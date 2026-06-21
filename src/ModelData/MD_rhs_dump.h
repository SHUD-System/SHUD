/* MD_rhs_dump.h — RHS snapshot dump hook target.
 *
 * Hook insertion (S0-6, openmp issue #8): MD_f.cpp / MD_update.cpp.
 * Writer impl (S0-7, openmp issue #9): MD_rhs_dump.cpp.
 *
 * Header is unconditionally guarded so includes from MD_f.cpp /
 * MD_update.cpp are safe regardless of SHUD_DUMP_RHS state; the
 * call sites themselves are #ifdef SHUD_DUMP_RHS-guarded so
 * DUMP=0 builds emit zero added code.
 *
 * Format version + magic MUST stay equal to the outer-repo
 *   tools/rhs_snapshot/format.h
 * The outer repo header is authoritative; this copy is a sister
 * kept in sync (S0-9 build check enforces).
 *
 * Layout (must match tools/rhs_snapshot/format.h v1):
 *   FileHeader   (40 B packed)  : magic[4] + version + case_id[32]
 *   RecordHeader (12 B packed)  : t_value (double) + array_count
 *   Array entries: uint32 name_len + name + uint64 nelem + double[nelem]
 *
 * Endianness: file is little-endian; writer assumes host == LE.
 */
#ifndef MD_RHS_DUMP_H
#define MD_RHS_DUMP_H

#include <cstdint>

#define SHUD_RHS_SNAPSHOT_FORMAT_VERSION 1
#define SHUD_RHS_SNAPSHOT_MAGIC          "SHRH"   /* 4 chars exactly, no NUL */

#ifdef SHUD_DUMP_RHS

#pragma pack(push, 1)

struct ShudSnapshotFileHeader {
    char     magic[4];
    uint32_t version;
    char     case_id[32];
};

struct ShudSnapshotRecordHeader {
    double   t_value;
    uint32_t array_count;
};

#pragma pack(pop)

static_assert(sizeof(ShudSnapshotFileHeader)   == 40,
              "ShudSnapshotFileHeader must be exactly 40 bytes; struct layout drifted vs outer repo");
static_assert(sizeof(ShudSnapshotRecordHeader) == 12,
              "ShudSnapshotRecordHeader must be exactly 12 bytes; struct layout drifted vs outer repo");

/* Record one RHS snapshot point.
 *
 * Arguments:
 *   site : short literal naming the call site ("f_update", "f_loop",
 *          "f_loop_before_passvalue", "f_applyDY"); writer matches
 *          against SHUD_DUMP_SITE env.
 *   t    : model time (SHUD t-unit = minutes; absolute from epoch).
 *   DY   : derivative vector, may be NULL when call site does not own DY.
 *   n    : length of DY; 0 when DY is NULL.
 *
 * Runtime behaviour controlled by env vars (see MD_rhs_dump.cpp):
 *   SHUD_DUMP_OUTPUT_DIR    default "."
 *   SHUD_DUMP_CASE_ID       default "unknown"
 *   SHUD_DUMP_T_VALUES      comma-separated doubles; if unset/empty, no-op
 *   SHUD_DUMP_T_TOL         default "0.5" (half model time unit, i.e. min)
 *   SHUD_DUMP_SITE          default "f_update"
 *   SHUD_DUMP_FNAME_SUFFIX  default "" (empty); when non-empty filename
 *                           becomes snapshot_t<v>_<suffix>.bin instead of
 *                           snapshot_t<v>.bin. Used by #43 before-PassValue_legacy
 *                           probe to coexist with existing f_update goldens
 *                           in same output dir without collision. Suffix
 *                           MUST NOT contain '/' or '\\' (path traversal
 *                           guard, checked first), and max length is 64
 *                           chars (F5 length cap, checked second); rejected
 *                           suffixes disable the dump and emit a stderr
 *                           diagnostic.
 */
void shud_rhs_dump_point(const char *site, double t,
                         const double *DY, int n);
#endif /* SHUD_DUMP_RHS */

#endif /* MD_RHS_DUMP_H */
