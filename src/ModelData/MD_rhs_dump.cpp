/* MD_rhs_dump.cpp — RHS snapshot dump stub (SHUD-OpenMP openmp issue #8).
 *
 * S0-6 ships a no-op so SHUD_DUMP_RHS=1 builds link cleanly while
 * preserving keliya output SHA256 vs DUMP=0 (the stub does NOT touch
 * any model state or write any file). The real writer impl — manifest-
 * driven probe + binary payload per `tools/rhs_snapshot/format.h` —
 * lands in S0-7 (openmp issue #9).
 *
 * When SHUD_DUMP_RHS is undefined the TU compiles to no emitted symbols,
 * which keeps the DUMP=0 link binary-functionally equivalent to a build
 * that omits this file entirely.
 */
#include "MD_rhs_dump.h"

#ifdef SHUD_DUMP_RHS
void shud_rhs_dump_point(const char *site, double t,
                         const double *DY, int n) {
    /* Stub: deliberately empty. S0-7 fills in the manifest-driven write. */
    (void) site;
    (void) t;
    (void) DY;
    (void) n;
}
#endif
