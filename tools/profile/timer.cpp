/* timer.cpp — profile bucket accumulator + YAML dump.
 *
 * Only compiled into the SHUD binary when SHUD_ENABLE_PROFILE is defined
 * at translation-unit scope (Makefile passes -DSHUD_ENABLE_PROFILE=1).
 * The whole impl is wrapped in #ifdef so a stray build of this TU
 * without the flag is a harmless no-op TU.
 */
#include "timer.h"

#ifdef SHUD_ENABLE_PROFILE

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace shud_profile {

namespace {

/* Bucket = monotonically-accumulated nanoseconds, lock-free atomic add.
 * We store `unique_ptr<atomic>` (rather than the atomic by value) so
 * the unordered_map can rehash without UB on atomic copy/move. */
using BucketMap =
    std::unordered_map<std::string,
                       std::unique_ptr<std::atomic<int64_t>>>;

struct Registry {
    BucketMap  buckets;
    std::mutex mu;
};

Registry &registry() {
    static Registry r;
    return r;
}

}  /* namespace */

void add_elapsed(const char *bucket, std::chrono::nanoseconds dt) {
    if (bucket == nullptr) return;
    const int64_t ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(dt).count();
    Registry &r = registry();
    std::atomic<int64_t> *slot = nullptr;
    {
        std::lock_guard<std::mutex> lk(r.mu);
        auto it = r.buckets.find(bucket);
        if (it == r.buckets.end()) {
            auto p = std::make_unique<std::atomic<int64_t>>(0);
            slot   = p.get();
            r.buckets.emplace(std::string(bucket), std::move(p));
        } else {
            slot = it->second.get();
        }
    }
    /* Hot path: lock-free atomic add. */
    slot->fetch_add(ns, std::memory_order_relaxed);
}

namespace {

/* Buckets emitted by `dump`. Order preserved so the YAML output is
 * deterministic; a bucket missing from the registry is written as 0.0. */
const char *const kCanonicalBuckets[] = {
    "t_RHS_kernel",
    "t_RHS_total",
    "t_CVODE_internal",
    "t_forcing_io",
    "t_ET",
    "t_output",
    "t_other",
};
constexpr size_t kCanonicalBucketCount =
    sizeof(kCanonicalBuckets) / sizeof(kCanonicalBuckets[0]);

double load_seconds(Registry &r, const char *name) {
    /* Caller holds r.mu. */
    auto it = r.buckets.find(name);
    if (it == r.buckets.end()) return 0.0;
    const int64_t ns = it->second->load(std::memory_order_relaxed);
    return static_cast<double>(ns) * 1e-9;
}

}  /* namespace */

void dump(const char *yaml_path) {
    if (yaml_path == nullptr) return;
    FILE *fp = std::fopen(yaml_path, "w");
    if (fp == nullptr) {
        /* Non-fatal: profile dump is observability, not correctness. */
        std::fprintf(stderr,
                     "[shud_profile] fopen('%s') failed; profile not "
                     "written.\n",
                     yaml_path);
        return;
    }
    Registry &r = registry();
    std::lock_guard<std::mutex> lk(r.mu);

    /* Derived buckets.
     *
     * The instrumented code (src/Model/shud.cpp, f.cpp, MD_ET.cpp,
     * Model_Control.cpp) emits 7 raw buckets:
     *
     *   t_RHS_total      — wraps every f() callback
     *   t_RHS_kernel     — the RHS kernel inside f() only
     *   t_CVODE_raw      — wraps every CVode() solver call (INCLUDES
     *                       the RHS callbacks fired from inside CVODE)
     *   t_forcing_io     — wraps Model_Data::updateforcing
     *   t_ET             — wraps Model_Data::ET
     *   t_output         — wraps Control_Data::ExportResults
     *   t_wall_total     — wraps the main solver loop (NumSteps iters)
     *
     * The emitted buckets additionally need:
     *
     *   t_CVODE_internal = t_CVODE_raw - t_RHS_total   (net solver time)
     *   t_other          = t_wall_total
     *                      - t_CVODE_raw
     *                      - t_forcing_io
     *                      - t_ET
     *                      - t_output
     *                      (residual; t_RHS_* already inside CVODE_raw)
     *
     * Both clamped to 0 to defend against negative residual that could
     * arise from RAII timer overhead double-counting or thread-local
     * jitter (we don't have wall_total ≥ sum_components as a hard
     * invariant; clamp is the safe degradation).
     *
     * The raw t_CVODE_raw + t_wall_total values stay in the `extras`
     * section so verifiers can sanity-check the derivation. */
    const double sec_rhs_total   = load_seconds(r, "t_RHS_total");
    const double sec_rhs_kernel  = load_seconds(r, "t_RHS_kernel");
    const double sec_cvode_raw   = load_seconds(r, "t_CVODE_raw");
    const double sec_forcing_io  = load_seconds(r, "t_forcing_io");
    const double sec_et          = load_seconds(r, "t_ET");
    const double sec_output      = load_seconds(r, "t_output");
    const double sec_wall_total  = load_seconds(r, "t_wall_total");

    const double sec_cvode_internal =
        (sec_cvode_raw - sec_rhs_total) > 0.0
            ? (sec_cvode_raw - sec_rhs_total) : 0.0;
    const double sec_other_raw =
        sec_wall_total - sec_cvode_raw - sec_forcing_io
                       - sec_et - sec_output;
    const double sec_other =
        sec_other_raw > 0.0 ? sec_other_raw : 0.0;

    std::fprintf(fp, "# SHUD profile (S0-10, openMP #14)\n");
    std::fprintf(fp, "# Units: seconds. t_CVODE_internal and t_other are\n");
    std::fprintf(fp, "# derived; raw measurements in extras: section.\n");
    std::fprintf(fp, "buckets:\n");
    std::fprintf(fp, "  t_RHS_kernel: %.9f\n",     sec_rhs_kernel);
    std::fprintf(fp, "  t_RHS_total: %.9f\n",      sec_rhs_total);
    std::fprintf(fp, "  t_CVODE_internal: %.9f\n", sec_cvode_internal);
    std::fprintf(fp, "  t_forcing_io: %.9f\n",     sec_forcing_io);
    std::fprintf(fp, "  t_ET: %.9f\n",             sec_et);
    std::fprintf(fp, "  t_output: %.9f\n",         sec_output);
    std::fprintf(fp, "  t_other: %.9f\n",          sec_other);

    /* Surface any non-canonical buckets (raw measurements used in the
     * derivation above, plus any future bucket the source emits before
     * canonical-list churn catches up). Skipped if absent. */
    bool extras_header_emitted = false;
    auto emit_extra = [&](const char *name, double sec) {
        if (!extras_header_emitted) {
            std::fprintf(fp, "extras:\n");
            extras_header_emitted = true;
        }
        std::fprintf(fp, "  %s: %.9f\n", name, sec);
    };
    /* Surface raw timers used as derivation inputs (always present
     * when their corresponding code site fired). */
    if (sec_cvode_raw > 0.0) emit_extra("t_CVODE_raw", sec_cvode_raw);
    if (sec_wall_total > 0.0) emit_extra("t_wall_total", sec_wall_total);
    /* Any other bucket the source emits not in the canonical-or-raw
     * set above also gets surfaced. */
    static const char *const kKnownRawOrCanonical[] = {
        "t_RHS_kernel", "t_RHS_total", "t_CVODE_internal", "t_forcing_io",
        "t_ET", "t_output", "t_other", "t_CVODE_raw", "t_wall_total",
    };
    constexpr size_t kKnownCount =
        sizeof(kKnownRawOrCanonical) / sizeof(kKnownRawOrCanonical[0]);
    for (const auto &kv : r.buckets) {
        bool is_known = false;
        for (size_t i = 0; i < kKnownCount; ++i) {
            if (kv.first == kKnownRawOrCanonical[i]) {
                is_known = true;
                break;
            }
        }
        if (is_known) continue;
        const int64_t ns =
            kv.second->load(std::memory_order_relaxed);
        emit_extra(kv.first.c_str(),
                   static_cast<double>(ns) * 1e-9);
    }
    /* CVODE stats placeholder: the values are in cvode_stats.txt and
     * are not read back here. */
    std::fprintf(fp, "cvode_stats:\n");
    /* `ratio_nfeLS_over_nfe` is a fixed top-level key of the yaml. It is
     * written as 0.0; compute the real value from cvode_stats.txt. */
    std::fprintf(fp, "ratio_nfeLS_over_nfe: 0.0\n");
    std::fclose(fp);
}

}  /* namespace shud_profile */

#endif  /* SHUD_ENABLE_PROFILE */
