/* timer.h — wall-clock profile timers.
 *
 * Compile-switched (SHUD_ENABLE_PROFILE=0|1). When the macro is undefined
 * the public API expands to no-ops, so the call sites in SHUD compile
 * as if the timers were not there. Built with `make shud SHUD_ENABLE_PROFILE=1`.
 *
 * Usage:
 *   #include "tools/profile/timer.h"
 *   void f() {
 *       shud_profile::Timer _t("t_RHS_kernel");   // RAII bucket accum
 *       ...
 *   }                                              // dtor adds elapsed
 *   shud_profile::dump("output/ccw.out/profile_B0.yaml");
 *
 * Thread safety: the bucket map is guarded by a mutex on insert/read;
 * once a bucket exists, accumulation goes through std::atomic<int64_t>
 * (nanoseconds) so the hot path is lock-free. `dump` takes the mutex
 * to snapshot the map.
 *
 * Buckets emitted by `dump`:
 *   t_RHS_kernel, t_RHS_total, t_CVODE_internal,
 *   t_forcing_io, t_ET, t_output, t_other
 * A bucket that was never timed is written as 0.0.
 *
 * No SUNDIALS / SHUD dependency.
 */
#ifndef SHUD_PROFILE_TIMER_H
#define SHUD_PROFILE_TIMER_H

#include <chrono>
#include <cstdint>

namespace shud_profile {

#ifdef SHUD_ENABLE_PROFILE

/* Real implementations live in timer.cpp. */
void add_elapsed(const char *bucket, std::chrono::nanoseconds dt);
void dump(const char *yaml_path);

class Timer {
public:
    explicit Timer(const char *bucket)
        : bucket_(bucket),
          t0_(std::chrono::steady_clock::now()) {}
    ~Timer() {
        add_elapsed(bucket_,
                    std::chrono::steady_clock::now() - t0_);
    }
    /* Copying a live timer would double-count; deleting copy/move keeps
     * the RAII semantics honest. */
    Timer(const Timer &)            = delete;
    Timer &operator=(const Timer &) = delete;
    Timer(Timer &&)                 = delete;
    Timer &operator=(Timer &&)      = delete;

private:
    const char                                       *bucket_;
    std::chrono::steady_clock::time_point             t0_;
};

#else  /* !SHUD_ENABLE_PROFILE — header-only no-ops. */

inline void add_elapsed(const char * /*bucket*/,
                        std::chrono::nanoseconds /*dt*/) {}
inline void dump(const char * /*yaml_path*/) {}

class Timer {
public:
    /* Empty inline ctor/dtor compile to zero instructions at -O2; the
     * unused-bucket arg is `(void)bucket;` cast to silence -Wunused. */
    explicit Timer(const char *bucket) { (void)bucket; }
    ~Timer() = default;
    Timer(const Timer &)            = delete;
    Timer &operator=(const Timer &) = delete;
    Timer(Timer &&)                 = delete;
    Timer &operator=(Timer &&)      = delete;
};

#endif  /* SHUD_ENABLE_PROFILE */

}  /* namespace shud_profile */

#endif  /* SHUD_PROFILE_TIMER_H */
