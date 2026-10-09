/* MD_diagnostics.hpp — diagnostic timer accumulators.
 *
 * 7-bucket RHS wall-clock timer + forcing I/O wall-clock timer.
 * ALL declarations and helper utilities here are gated behind
 * `#ifdef SHUD_ENABLE_DIAGNOSTICS` so the default build (macro
 * undefined) emits ZERO code in the RHS hot path and ZERO new symbols
 * in the binary, and its results are unaffected.
 *
 * Buckets:
 *   0 update    — rhs_update() (zero-resets + BC reads)
 *   1 ET        — rhs_flux() pass-1 per-element ET/infil/recharge
 *   2 lateral   — rhs_flux() pass-2 per-element surf/sub flux
 *   3 segment   — rhs_flux() per-segment surf/sub
 *   4 river     — rhs_flux() per-river Flux_RiverDown + lake clamp
 *   5 gather    — rhs_flux() lake transitional gather + rhs_deterministic_gather()
 *   6 applyDY   — rhs_apply() write of DY
 *
 * Timer storage is a single global `long long` array (nanoseconds).
 * A plain global accumulator is race-free only while the RHS runs
 * serially. With the OpenMP RHS (SHUD_ENABLE_OPENMP_RHS=1) the bucket
 * timers inside rhs_flux() are updated by every team thread without
 * synchronisation, so the timings are approximate; use diagnostics
 * builds for profiling only.
 */
#ifndef MD_DIAGNOSTICS_HPP
#define MD_DIAGNOSTICS_HPP

#ifdef SHUD_ENABLE_DIAGNOSTICS

#include <chrono>

namespace shud_diag {

enum RhsBucket {
    RHS_BUCKET_UPDATE  = 0,
    RHS_BUCKET_ET      = 1,
    RHS_BUCKET_LATERAL = 2,
    RHS_BUCKET_SEGMENT = 3,
    RHS_BUCKET_RIVER   = 4,
    RHS_BUCKET_GATHER  = 5,
    RHS_BUCKET_APPLYDY = 6,
    RHS_BUCKET_COUNT   = 7
};

/* Per-bucket nanosecond accumulators (defined in MD_rhs_core.cpp). */
extern long long g_rhs_timer_ns[RHS_BUCKET_COUNT];

/* Forcing I/O accumulator (defined in TimeSeriesData.cpp). */
extern long long g_forcing_io_ns;

/* RAII timer: on construction records steady_clock::now(); on
 * destruction adds (now - start) nanoseconds to *target_ns.
 * Zero floating-point operations; chrono::duration_cast<ns> is integer
 * arithmetic in libstdc++ / libc++. */
class ScopeTimer {
public:
    explicit ScopeTimer(long long *target_ns)
      : target_ns_(target_ns),
        t_start_(std::chrono::steady_clock::now()) {}
    ~ScopeTimer() {
        auto t_end = std::chrono::steady_clock::now();
        *target_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            t_end - t_start_).count();
    }
private:
    long long *target_ns_;
    std::chrono::steady_clock::time_point t_start_;
};

}  // namespace shud_diag

#endif  /* SHUD_ENABLE_DIAGNOSTICS */

#endif  /* MD_DIAGNOSTICS_HPP */
