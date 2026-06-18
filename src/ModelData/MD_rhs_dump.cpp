/* MD_rhs_dump.cpp — RHS snapshot dump impl.
 *
 * S0-7 (openmp issue #9) replaces the S0-6 stub with a real writer:
 * env-driven manifest probe, single-record / single-array snapshot
 * per call site match.
 *
 * When SHUD_DUMP_RHS is undefined the TU compiles to no emitted symbols
 * — DUMP=0 link binary stays functionally equivalent to a build that
 * omits this file entirely (B0 SHA256 invariant).
 *
 * SCHEMA DUPLICATION NOTE:
 *   The byte layout written here MUST match the writer at
 *     tools/rhs_snapshot/writer.cpp
 *   and the schema header at
 *     tools/rhs_snapshot/format.h
 *   We replicate the small writer here (rather than link against the
 *   outer-repo writer) because the SHUD submodule has no visibility of
 *   outer-repo objects at compile time. S0-9 CI MUST diff the two
 *   schema headers + replicate-writer logic on each PR.
 */
#include "MD_rhs_dump.h"

#ifdef SHUD_DUMP_RHS

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

struct DumpConfig {
    bool                initialized = false;
    bool                disabled    = false;  /* no t_values → silent no-op */
    std::string         output_dir;
    std::string         case_id;
    std::string         site;
    std::string         fname_suffix;  /* #43: empty = legacy snapshot_t<v>.bin;
                                        * non-empty = snapshot_t<v>_<suffix>.bin */
    double              tol = 0.5;
    std::vector<double> targets;
    std::vector<bool>   consumed;
};

DumpConfig& cfg() {
    static DumpConfig c;
    return c;
}

std::vector<double> parse_t_values(const char *raw) {
    std::vector<double> out;
    if (!raw || raw[0] == '\0') return out;
    const std::string s = raw;
    std::size_t i = 0;
    while (i < s.size()) {
        std::size_t j = s.find(',', i);
        if (j == std::string::npos) j = s.size();
        std::string tok = s.substr(i, j - i);
        /* trim ASCII whitespace */
        while (!tok.empty() && (tok.front() == ' ' || tok.front() == '\t')) tok.erase(tok.begin());
        while (!tok.empty() && (tok.back() == ' ' || tok.back() == '\t')) tok.pop_back();
        if (!tok.empty()) {
            char *end = nullptr;
            double v = std::strtod(tok.c_str(), &end);
            if (end != tok.c_str()) {
                out.push_back(v);
            }
        }
        i = j + 1;
    }
    return out;
}

void init_config() {
    DumpConfig& c = cfg();
    if (c.initialized) return;
    c.initialized = true;

    const char *dir = std::getenv("SHUD_DUMP_OUTPUT_DIR");
    c.output_dir = (dir && dir[0]) ? dir : ".";

    const char *cid = std::getenv("SHUD_DUMP_CASE_ID");
    c.case_id = (cid && cid[0]) ? cid : "unknown";

    const char *site = std::getenv("SHUD_DUMP_SITE");
    c.site = (site && site[0]) ? site : "f_update";

    /* #43: filename disambiguation suffix for coexistence of multiple
     * dump sites in same output dir (e.g. f_update vs
     * f_loop_before_passvalue). Empty = legacy filename, full bitwise
     * back-compat with PR #53 goldens. */
    const char *sfx = std::getenv("SHUD_DUMP_FNAME_SUFFIX");
    c.fname_suffix = (sfx && sfx[0]) ? sfx : "";
    /* Path-traversal guard: reject suffix containing '/' or '\\'. The
     * suffix lands in a snprintf("snapshot_t%.0f_%s.bin") format used
     * to compose a path under SHUD_DUMP_OUTPUT_DIR; a slash would
     * escape that dir. Failure mode: disable the dump and emit a
     * diagnostic so the calling test harness sees the rejection. */
    if (!c.fname_suffix.empty()) {
        if (c.fname_suffix.find('/')  != std::string::npos ||
            c.fname_suffix.find('\\') != std::string::npos) {
            std::fprintf(stderr,
                "shud_rhs_dump: SHUD_DUMP_FNAME_SUFFIX '%s' contains "
                "path separator; rejecting and disabling dump\n",
                c.fname_suffix.c_str());
            c.disabled = true;
            return;
        }
        /* F5 fix: length cap. Buffer math in shud_rhs_dump_point()
         * leaves ~95 chars headroom for the suffix; a 64-char limit
         * gives us a clear safety margin and keeps filenames
         * filesystem-friendly. Truncating silently would let a
         * misconfigured caller produce surprise paths, so we reject
         * + disable + diagnose. */
        if (c.fname_suffix.size() > 64) {
            std::fprintf(stderr,
                "shud_rhs_dump: SHUD_DUMP_FNAME_SUFFIX exceeds 64-char "
                "limit (got %zu); rejecting and disabling dump\n",
                c.fname_suffix.size());
            c.disabled = true;
            return;
        }
        /* F24 (PR #54 round-2): The two early-return paths above (path-
         * separator reject + F5 length-cap reject) intentionally leave
         * c.consumed unset.  This is safe by design: c.disabled = true
         * makes downstream call paths (shud_rhs_dump_point body, see
         * `if (c.disabled) return;` at the top of the dispatch) early-
         * out before any read of c.consumed[i].  c.consumed is only
         * assigned at line 172 below, after we've committed to a valid
         * targets vector — so any future reorder that touches consumed
         * before disabled-gate must re-audit these early returns. */
    }

    const char *tol = std::getenv("SHUD_DUMP_T_TOL");
    if (tol && tol[0]) {
        char *end = nullptr;
        double v = std::strtod(tol, &end);
        if (end != tol) c.tol = v;
    }

    const char *tv = std::getenv("SHUD_DUMP_T_VALUES");
    c.targets = parse_t_values(tv);

    /* F2-corr hygiene: drop non-finite targets (parse_t_values accepts
     * "nan"/"inf" tokens because strtod does; reject them here so they
     * cannot reach the filename buffer or the tolerance compare). */
    {
        auto fin = c.targets.begin();
        for (double v : c.targets) {
            if (std::isfinite(v)) *fin++ = v;
            else std::fprintf(stderr,
                "shud_rhs_dump: dropping non-finite target from SHUD_DUMP_T_VALUES\n");
        }
        c.targets.erase(fin, c.targets.end());
    }

    /* F1-corr hygiene: %%.0f filename precision means two targets within
     * 1.0 of each other collapse to the same snapshot_t<rounded>.bin file
     * and the second would silently overwrite the first. Reject at init. */
    if (c.targets.size() > 1) {
        std::vector<double> sorted(c.targets);
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t i = 1; i < sorted.size(); ++i) {
            if (std::fabs(sorted[i] - sorted[i - 1]) < 1.0) {
                std::fprintf(stderr,
                    "shud_rhs_dump: SHUD_DUMP_T_VALUES targets %g and %g "
                    "collide under %%.0f filename precision; disabling\n",
                    sorted[i - 1], sorted[i]);
                c.disabled = true;
                return;
            }
        }
    }

    if (c.targets.empty()) {
        c.disabled = true;
    }
    c.consumed.assign(c.targets.size(), false);
}

/* Replicated single-record writer. Layout: see SCHEMA DUPLICATION NOTE
 * above + format.h header. */
int write_one_snapshot(const std::string& path,
                       const std::string& case_id,
                       double             t_value,
                       const std::string& array_name,
                       const double      *data,
                       uint64_t           nelem) {
    FILE *fp = std::fopen(path.c_str(), "wb");
    if (!fp) return -1;

    ShudSnapshotFileHeader fh;
    std::memcpy(fh.magic, SHUD_RHS_SNAPSHOT_MAGIC, 4);
    fh.version = static_cast<uint32_t>(SHUD_RHS_SNAPSHOT_FORMAT_VERSION);
    std::memset(fh.case_id, 0, sizeof(fh.case_id));
    {
        const std::size_t n = case_id.size() < sizeof(fh.case_id)
                                ? case_id.size() : sizeof(fh.case_id);
        std::memcpy(fh.case_id, case_id.data(), n);
    }
    std::fwrite(&fh, sizeof(fh), 1, fp);

    ShudSnapshotRecordHeader rh;
    rh.t_value     = t_value;
    rh.array_count = 1;
    std::fwrite(&rh, sizeof(rh), 1, fp);

    const uint32_t name_len = static_cast<uint32_t>(array_name.size());
    std::fwrite(&name_len, sizeof(name_len), 1, fp);
    if (name_len > 0) {
        std::fwrite(array_name.data(), 1, name_len, fp);
    }
    const uint64_t nelem_le = nelem;
    std::fwrite(&nelem_le, sizeof(nelem_le), 1, fp);
    if (nelem > 0 && data != nullptr) {
        std::fwrite(data, sizeof(double), nelem, fp);
    }

    std::fclose(fp);
    return 0;
}

}  /* namespace */

void shud_rhs_dump_point(const char *site, double t,
                         const double *DY, int n) {
    init_config();
    DumpConfig& c = cfg();
    if (c.disabled) return;
    if (DY == nullptr || n <= 0) return;
    if (site == nullptr) return;
    if (c.site != site) return;

    /* Find first unconsumed target within tolerance. */
    int idx = -1;
    double best_dist = 0.0;
    for (std::size_t i = 0; i < c.targets.size(); ++i) {
        if (c.consumed[i]) continue;
        const double d = std::fabs(t - c.targets[i]);
        if (d <= c.tol) {
            if (idx < 0 || d < best_dist) {
                idx = static_cast<int>(i);
                best_dist = d;
            }
        }
    }
    if (idx < 0) return;
    c.consumed[idx] = true;

    /* #43: when fname_suffix is set, append `_<suffix>` between
     * t-stem and `.bin`. Filename buffer sized to fit the longest
     * suffix the path-traversal guard does not reject. Buffer math:
     * 128 - 10 ("snapshot_t") - 17 ("%.0f" max) - 1 ("_") - 4 (".bin")
     * - 1 (NUL) = 95 chars headroom for suffix. The F5 64-char
     * SHUD_DUMP_FNAME_SUFFIX guard in init_config() keeps any accepted
     * suffix well below this limit. */
    char fname[128];
    if (c.fname_suffix.empty()) {
        std::snprintf(fname, sizeof(fname), "snapshot_t%.0f.bin",
                      c.targets[idx]);
    } else {
        std::snprintf(fname, sizeof(fname), "snapshot_t%.0f_%s.bin",
                      c.targets[idx], c.fname_suffix.c_str());
    }
    std::string path = c.output_dir;
    if (!path.empty() && path.back() != '/') path.push_back('/');
    path += fname;

    write_one_snapshot(path, c.case_id, c.targets[idx],
                       std::string("DY"), DY,
                       static_cast<uint64_t>(n));
}

#endif /* SHUD_DUMP_RHS */
