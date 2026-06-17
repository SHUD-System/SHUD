/* MD_rhs_dump.h — RHS snapshot dump hook target (SHUD-OpenMP openmp issue #8).
 *
 * Stub declarations only. Real writer impl lands in S0-7 (issue #9).
 *
 * Header is unconditionally guarded so it's safe to include from
 * MD_f.cpp / MD_update.cpp regardless of SHUD_DUMP_RHS state; the
 * function-call sites themselves are #ifdef SHUD_DUMP_RHS-guarded
 * in the caller, so this header has no effect on a DUMP=0 build.
 *
 * Format version MUST stay equal to the outer-repo
 * `tools/rhs_snapshot/format.h` SHUD_RHS_SNAPSHOT_FORMAT_VERSION.
 * Outer repo header is authoritative; this copy is a sister kept in
 * sync (S0-7 build check enforces).
 */
#ifndef MD_RHS_DUMP_H
#define MD_RHS_DUMP_H

#define SHUD_RHS_SNAPSHOT_FORMAT_VERSION 1

#ifdef SHUD_DUMP_RHS
/* Record one RHS snapshot point. `site` is a short literal naming the
 * call site (e.g. "f_update", "f_loop", "f_applyDY"); `t` is the model
 * time at which the hook fires; `DY` is the optional derivative vector
 * of length `n` (n=0 / DY=NULL when the call site does not own DY).
 *
 * S0-6 ships a no-op stub (defined in MD_rhs_dump.cpp); S0-7 replaces
 * the body with the real writer (manifest-driven probe at the
 * configured t_values, binary payload per format.h).
 */
void shud_rhs_dump_point(const char *site, double t,
                         const double *DY, int n);
#endif /* SHUD_DUMP_RHS */

#endif /* MD_RHS_DUMP_H */
