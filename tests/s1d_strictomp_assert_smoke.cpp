/* tests/s1d_strictomp_assert_smoke.cpp -- S1d.1 (openMP #47).
 *
 * Regression guard for the decision to use `std::abort()` (NOT
 * `assert(false)`) inside the `ExecPolicy::StrictOMP` and
 * `ExecPolicy::ProductionOMP` cases of `Model_Data::rhs_core`.
 *
 * Why: under `-DNDEBUG` (release builds, and the
 * `EXTRA_CXXFLAGS=-DNDEBUG` smoke compile invoked by the
 * `smoke_strictomp` Makefile target), `assert(...)` expands to an
 * empty statement. If those switch cases used `assert(false)` they
 * would compile away into no-ops, and execution would silently fall
 * through to the next statement -- in this layout, to the closing
 * brace of `rhs_core` -- giving the caller a normal return rather
 * than the fail-fast abort the spec demands. That would let an
 * `ExecPolicy::StrictOMP` call impersonate Serial without anyone
 * noticing, defeating the entire S1d.1 dispatch contract.
 *
 * Test mechanism:
 *   - `fork()` a child process.
 *   - In the child: construct a Model_Data on the heap (members
 *     uninitialized but we never reach any of them; the StrictOMP
 *     case `std::abort()`s before reading any state), then call
 *     `rhs_core(Y, DY, t, ExecPolicy::StrictOMP)`. The abort
 *     terminates the child via SIGABRT; the heap Model_Data is
 *     never destructed (no FreeData() UB). If somehow the call
 *     returns, `_exit(0)` makes the assertion below fail.
 *   - In the parent: `waitpid` the child, check
 *     `WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT`.
 *     Exit 0 on PASS; exit 1 on FAIL (driving the Makefile target).
 *
 * Build / run is wired through the `smoke_strictomp` target in
 * SHUD/Makefile, which compiles with
 * `-DNDEBUG -DSHUD_ENABLE_OPENMP_RHS=1`. The `-DSHUD_ENABLE_OPENMP_RHS=1`
 * is needed because the StrictOMP / ProductionOMP cases are
 * `#ifdef SHUD_ENABLE_OPENMP_RHS`-gated out of the translation unit
 * by default to keep release binaries free of OMP-path symbols.
 *
 * Scope: this test deliberately does NOT exercise the Serial path
 * (that is covered by the standard 4-case bitwise validation under
 * tasks 4.6a/b). It tests one thing -- that the OMP-policy stubs
 * abort -- and nothing else.
 */

#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <cstdio>
#include <cstdlib>

#include "Model_Data.hpp"
#include "MD_rhs_core.hpp"

int main(void) {
    pid_t pid = fork();
    if (pid < 0) {
        std::fprintf(stderr, "FAIL: fork failed\n");
        return 1;
    }
    if (pid == 0) {
        /* Child: heap-allocate so the abort doesn't drag a stack
         * dtor + FreeData() into the signal-disposition window.
         * `new` may itself allocate but cannot fail in any
         * realistic test environment; if it throws we exit non-
         * zero and the parent reports it as not-SIGABRT.
         */
        Model_Data *MD = new Model_Data();
        double Y[1]  = {0.0};
        double DY[1] = {0.0};
        MD->rhs_core(Y, DY, 0.0, ExecPolicy::StrictOMP);
        /* If we reach here, the abort didn't fire. Exit non-zero
         * via _exit (not exit) to skip atexit handlers / static
         * dtors and keep the parent's waitpid signal interpretation
         * clean. _exit(0) -> WIFEXITED true, not WIFSIGNALED ->
         * parent reports FAIL. */
        _exit(0);
    }

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        std::fprintf(stderr, "FAIL: waitpid failed\n");
        return 1;
    }
    if (WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) {
        std::printf("PASS: StrictOMP stub aborted with SIGABRT under -DNDEBUG\n");
        return 0;
    }
    if (WIFEXITED(status)) {
        std::fprintf(stderr,
                     "FAIL: child exited normally with code %d "
                     "(expected SIGABRT) -- StrictOMP stub did not abort; "
                     "assert(false) regression?\n",
                     WEXITSTATUS(status));
    } else if (WIFSIGNALED(status)) {
        std::fprintf(stderr,
                     "FAIL: child died by signal %d (expected SIGABRT=%d)\n",
                     WTERMSIG(status), SIGABRT);
    } else {
        std::fprintf(stderr, "FAIL: unexpected wait status 0x%x\n", status);
    }
    return 1;
}
