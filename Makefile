# -----------------------------------------------------------------
# Makefile for SHUD — B0 baseline lock
# -----------------------------------------------------------------
# Programmer: Lele Shu (lele.shu@gmail.com)
# SHUD model is a heritage of Penn State Integrated Hydrologic Model (PIHM).
# B0 lockdown maintained per docs/build_manifest.md (top-level repo).
# Any change to the flag set requires an OpenSpec change against
#   openspec/changes/.../specs/build-environment-lockdown/spec.md
# -----------------------------------------------------------------
# Prerequisites:
#   - SUNDIALS 6.x installed at $(SUNDIALS_DIR) via ./configure
#   - For OpenMP on macOS: `brew install libomp`
#   - For OpenMP on Linux: GCC with libgomp
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# B0 baseline lock — flag variables (do not edit without OpenSpec change)
# -----------------------------------------------------------------
# Two-tier override: CXX_BASE_FLAGS is canonical; SHUD_BUILD_CFLAGS is an
# alias kept for backward-compat with readers that grep the recipe line.
# Both use GNU make's `override … :=` so a `make VAR=…` CLI override is
# silently ignored (per GNU make manual: override on `:=` immunizes against
# make-CLI assignment). Both must stay `override`-protected.
override CXX_BASE_FLAGS    := -O2 -g -ffp-contract=off -fno-fast-math -std=c++14
override SHUD_BUILD_CFLAGS := $(CXX_BASE_FLAGS)
# S1d.2 (openMP #48) — the legacy `CXX_OPENMP_DEFINE` variable +
# its single-concern `-D` define have been retired. The legacy
# triple-concern switch is replaced by two orthogonal macros
# (defined below), each gating exactly one concern:
#   SHUD_USE_OPENMP_NVECTOR — N_Vector backend (Serial vs OpenMP)
#   SHUD_ENABLE_OPENMP_RHS  — RHS execution policy stubs (StrictOMP /
#                             ProductionOMP); from #47
# S2 capstone (PR-8) — the third macro (legacy `_omp` RHS receiver
# compile inclusion) was retired together with the source file it
# gated; the receivers no longer exist in the codebase.
# The historical legacy define no longer exists in any source.

# CFLAGS is left as a non-override alias for backward-compat tooling that
# reads $(CFLAGS); recipes invoke $(SHUD_BUILD_CFLAGS) directly, so a user-
# supplied `make CFLAGS=…` cannot clobber the locked flag set even before
# the disallowed-flag scan below catches the injection.
CFLAGS            = $(CXX_BASE_FLAGS)

# -----------------------------------------------------------------
# B0 baseline lock — disallowed-flag guard
# -----------------------------------------------------------------
# Fast-fail if the user attempts to inject IEEE-754-violating flags via
# ANY user-controllable flag carrier. Two layers:
#
# Layer 1 (filter): scans the 6 standard carriers (CFLAGS / CXXFLAGS /
# CPPFLAGS / LDFLAGS / MAKEOVERRIDES / MAKEFLAGS) and the 2 project-local
# lock variables (SHUD_BUILD_CFLAGS / CXX_BASE_FLAGS). `filter` works
# word-level, so it catches `CXXFLAGS=-ffast-math …`. The two project-
# local variables are also `override`-protected above (so their values
# always equal the locked flag set); the scan extension is defense-in-
# depth in case a future refactor reverts the `override`.
#
# Layer 2 (anchored `=value` scan on MAKEOVERRIDES): iterates over the
# VAR=value tokens of $(MAKEOVERRIDES) and uses `filter %=$(f)` to detect
# exact-match `VAR=<disallowed-flag>` CLI assignments (e.g.
# `make SHUD_BUILD_CFLAGS=-Ofast`), which the `override :=` directive on
# the lock variables would otherwise silently ignore. The earlier
# `findstring -Ofast,$(MAKEOVERRIDES)` form was a literal-substring scan
# that false-positived on paths like `SUNDIALS_DIR=/opt/sundials-Ofast-tuned`;
# the anchored form pins to the `=` boundary and the end of the token, so
# only true `VAR=-Ofast` CLI assignments fire. Without this layer the
# user could think "the build accepted my flag" while in reality the
# locked flag set is still used — Layer 2 turns the silent rejection into
# a loud error.
#
# `DISALLOWED_FLAGS` itself is `override`-protected so a user-supplied
# `make DISALLOWED_FLAGS=` cannot disarm the scan list. Binary safety is
# also guaranteed by `override :=` on the lock variables above, but
# Layer 1/2 are the user-facing UX — keep them armed.
override DISALLOWED_FLAGS := -ffast-math -Ofast -funsafe-math-optimizations
ifneq (,$(filter $(DISALLOWED_FLAGS),$(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(MAKEOVERRIDES) $(MAKEFLAGS) $(SHUD_BUILD_CFLAGS) $(CXX_BASE_FLAGS)))
$(error disallowed flag detected (one of $(DISALLOWED_FLAGS)) in CFLAGS/CXXFLAGS/CPPFLAGS/LDFLAGS/MAKEOVERRIDES/MAKEFLAGS/SHUD_BUILD_CFLAGS/CXX_BASE_FLAGS; B0 baseline requires strict IEEE-754 — see docs/build_manifest.md §1)
endif
# Layer 2: anchored scan of $(MAKEOVERRIDES) for `VAR=<disallowed-flag>` CLI
# assignments. We split MAKEOVERRIDES into VAR=value tokens and, for each
# token, look for an exact-match `VAR=<flag>` (via `filter %=$(f)`). This
# avoids the previous literal-findstring false-positive on paths/values that
# legitimately contain a disallowed flag as a substring (e.g.
# `SUNDIALS_DIR=/opt/sundials-Ofast-tuned`), while still catching
# `make shud SHUD_BUILD_CFLAGS=-Ofast` / `CXX_BASE_FLAGS=-Ofast`, which the
# `override :=` directive would otherwise silently drop.
LAYER2_HITS := $(strip $(foreach tok,$(MAKEOVERRIDES),$(foreach f,$(DISALLOWED_FLAGS),$(if $(filter %=$(f),$(tok)),$(f)))))
ifneq (,$(LAYER2_HITS))
$(error disallowed flag detected ($(LAYER2_HITS)) in a make-CLI assignment (MAKEOVERRIDES=[$(MAKEOVERRIDES)]); attempts to inject via SHUD_BUILD_CFLAGS / CXX_BASE_FLAGS / etc are also rejected; B0 baseline requires strict IEEE-754 — see docs/build_manifest.md §1)
endif

# -----------------------------------------------------------------
# Optional RHS snapshot dump instrumentation (openmp issue #8 / S0-6)
# -----------------------------------------------------------------
# Off by default. Set `SHUD_DUMP_RHS=1` on the make CLI to compile in
# the dump hooks at f_update / f_loop / f_applyDY exits. Hook bodies
# are #ifdef SHUD_DUMP_RHS-guarded in the source, so DUMP=0 produces
# preprocessor output identical to the unmodified codebase — keliya
# `*.dat` SHA256 must match the pre-#8 baseline at the same flag set.
#
# Runtime env vars (consumed when SHUD_DUMP_RHS=1):
#   SHUD_DUMP_OUTPUT_DIR  — directory for snapshot files (default: cwd)
#   SHUD_DUMP_MANIFEST    — path to benchmark manifest.yaml driving
#                           probe `t_values` (see benchmarks/<case>/)
#
# Writer impl is S0-7 (#9). This stage ships a no-op stub
# (SHUD/src/ModelData/MD_rhs_dump.cpp) so SHUD_DUMP_RHS=1 builds link
# and the main keliya output stays bitwise-equal to DUMP=0.
SHUD_DUMP_RHS ?= 0
ifeq ($(SHUD_DUMP_RHS),1)
  SHUD_DUMP_DEFINE := -DSHUD_DUMP_RHS=1
else ifeq ($(SHUD_DUMP_RHS),0)
  SHUD_DUMP_DEFINE :=
else
$(error SHUD_DUMP_RHS must be 0 or 1, got '$(SHUD_DUMP_RHS)')
endif

# -----------------------------------------------------------------
# SHUD_ENABLE_OPENMP_RHS — S1d.1 (openMP issue #47); legacy fork retired in S2 capstone (PR-8 #152)
# -----------------------------------------------------------------
# The S2 capstone (PR-8) retired the legacy-vs-rhs_core fork in f.cpp:
# f() now unconditionally routes through `rhs_core(ExecPolicy::Serial)`,
# so the original `f_update/f_loop/f_applyDY` chain is dead source kept
# only as the PURE CARRY-OVER source of `rhs_update/rhs_flux/rhs_apply`.
# The legacy-routing macro that previously gated f.cpp has been removed.
#
# SHUD_ENABLE_OPENMP_RHS (default 0):
#   0 = `#ifdef SHUD_ENABLE_OPENMP_RHS` cases in `MD_rhs_core.cpp`
#       (StrictOMP / ProductionOMP) are excluded from the translation
#       unit; resulting binary has no OMP-path symbols.
#   1 = OMP cases compile in. They contain `std::abort()` stubs that
#       SIGABRT on any runtime call (verified by
#       `tests/s1d_strictomp_assert_smoke.cpp` under -DNDEBUG). Used
#       only for smoke compile / abort regression; NOT bitwise-validated.
#
# `assert(false)` is forbidden inside the OMP stubs because
# `-DNDEBUG` strips assert to a no-op and would let the switch case
# fall through silently to the next statement. `std::abort()` is
# unconditional; safe under `EXTRA_CXXFLAGS=-DNDEBUG` smoke compile.
# Release v1.0 default flip: `make shud_omp` defaults to
# SHUD_ENABLE_OPENMP_RHS=1 (Config C, Serial NVec + StrictOMP RHS,
# ADR-0002 Path 1 winner). `make shud` remains SHUD_ENABLE_OPENMP_RHS=0
# (Config A, canonical serial reference). User explicit override on the
# command line (e.g. `make shud_omp SHUD_ENABLE_OPENMP_RHS=0`) still
# wins via `?=` semantics — the research escape hatch is preserved.
# NB: `?=` late-evaluates SHUD_ENABLE_OPENMP_RHS; target-specific vars
# would NOT flow into the `ifeq` block below (parse-time vs recipe-time
# phases), so MAKECMDGOALS filtering is the correct hook.
ifneq (,$(filter shud_omp,$(MAKECMDGOALS)))
  SHUD_ENABLE_OPENMP_RHS ?= 1
else
  SHUD_ENABLE_OPENMP_RHS ?= 0
endif
ifeq ($(SHUD_ENABLE_OPENMP_RHS),0)
  SHUD_OMP_RHS_DEFINE :=
  SHUD_OMP_RHS_CK     :=
  SHUD_OMP_RHS_LK     :=
else ifeq ($(SHUD_ENABLE_OPENMP_RHS),1)
  SHUD_OMP_RHS_DEFINE := -DSHUD_ENABLE_OPENMP_RHS=1
  # P1e PR-G (#315 / design D3 + D6) — StrictOMP path needs the OpenMP
  # runtime so the outer `#pragma omp parallel` in
  # MD_rhs_core.cpp::ExecPolicy::StrictOMP actually spawns threads
  # (without -fopenmp the `_OPENMP` macro is undefined and the directive
  # collapses to serial code). Mirrors the SHUD_USE_OPENMP_NVECTOR=1
  # pattern: `=` (recursive expansion) is required because
  # CXX_OPENMP_CFLAGS / CXX_OPENMP_LFLAGS are defined later in the
  # Makefile (Platform-conditional OpenMP flags block); using `:=` here
  # would capture an empty value and silently drop -fopenmp / -lomp /
  # -lgomp from the link line, producing link errors at build time.
  # Wires Config C (`make shud SHUD_ENABLE_OPENMP_RHS=1`) into a real
  # multi-threaded binary; Config D (`make shud_omp SHUD_ENABLE_OPENMP_RHS=1`)
  # already inherits -fopenmp from the shud_omp recipe.
  SHUD_OMP_RHS_CK      = $(CXX_OPENMP_CFLAGS)
  SHUD_OMP_RHS_LK      = $(CXX_OPENMP_LFLAGS)
else
$(error SHUD_ENABLE_OPENMP_RHS must be 0 or 1, got '$(SHUD_ENABLE_OPENMP_RHS)')
endif

# -----------------------------------------------------------------
# S1d.2 (openMP #48) — SHUD_USE_OPENMP_NVECTOR
# -----------------------------------------------------------------
# This flag controls the N_Vector backend selection. S2 capstone
# (PR-8) retired the sibling `_omp` RHS receiver compile-inclusion
# switch together with the source file it gated; `SHUD_USE_OPENMP_NVECTOR`
# is now the only N_Vector concern remaining. Defaults OFF so the
# default build is Config A (bitwise vs B0).
#
# SHUD_USE_OPENMP_NVECTOR (default 0):
#   0 = nvector_serial backend. udata / du are allocated via
#       N_VNew_Serial; SET_VALUE / N_VGetArrayPointer dispatch to
#       the Serial ops table. No `-lsundials_nvecopenmp` link
#       dependency; the resulting binary does NOT pull in
#       libsundials_nvecopenmp.so / .dylib (verified via
#       `otool -L shud | grep nvecopenmp` on macOS or
#       `ldd shud | grep nvecopenmp` on Linux).
#   1 = nvector_openmp backend. udata / du via N_VNew_OpenMP
#       (NV_THREADS = MD->CS.num_threads). Pulls in
#       nvector_openmp.h (Macros.hpp) + links
#       libsundials_nvecopenmp. Independent of
#       SHUD_ENABLE_OPENMP_RHS (which gates the in-RHS
#       parallel-kernel stubs).
SHUD_USE_OPENMP_NVECTOR ?= 0
ifeq ($(SHUD_USE_OPENMP_NVECTOR),0)
  SHUD_NVEC_OMP_DEFINE :=
  SHUD_NVEC_OMP_CK     :=
  SHUD_NVEC_OMP_LK     :=
else ifeq ($(SHUD_USE_OPENMP_NVECTOR),1)
  SHUD_NVEC_OMP_DEFINE := -DSHUD_USE_OPENMP_NVECTOR=1
  # OpenMP NVector backend needs the OpenMP runtime (omp_set_num_threads
  # in shud.cpp + parallelized N_Vector ops inside libsundials_nvecopenmp).
  # When the user enables it via `make shud SHUD_USE_OPENMP_NVECTOR=1` we
  # implicitly pull in the platform OpenMP compile flag (sets `_OPENMP`)
  # and the OpenMP runtime link flag. This keeps `shud` + Config B/C/D
  # composable without forcing the user to also pass a separate -fopenmp.
  #
  # `=` (recursive expansion) is required here because CXX_OPENMP_CFLAGS
  # / CXX_OPENMP_LFLAGS are defined further down in the Makefile (in the
  # "Platform-conditional OpenMP flags" block). Using `:=` here would
  # capture an empty value (the variables aren't bound yet at this
  # point in parse order), silently dropping `-fopenmp` and `-lomp` /
  # `-lgomp` from the recipe and producing link errors at build time.
  SHUD_NVEC_OMP_CK      = $(CXX_OPENMP_CFLAGS)
  SHUD_NVEC_OMP_LK      = $(CXX_OPENMP_LFLAGS) -lsundials_nvecopenmp
else
$(error SHUD_USE_OPENMP_NVECTOR must be 0 or 1, got '$(SHUD_USE_OPENMP_NVECTOR)')
endif

# -----------------------------------------------------------------
# P12-nvec PR-N1 (#443) — SHUD_NVEC_HYBRID (Config E)
# -----------------------------------------------------------------
# Config E = Config C StrictOMP RHS + OpenMP NVector element-wise ops +
# SHUD-owned SERIAL reduction overrides (src/Model/MD_nvec_hybrid.cpp,
# whole TU #ifdef SHUD_NVEC_HYBRID). Element-wise ops keep the stock
# parallel implementation; every reduction slot is overwritten with a
# serial generic-API loop → bitwise-identical across thread counts AND vs
# Config C (Serial NVector). Design D2 / spec hybrid-nvec-tier1.
#
# SHUD_NVEC_HYBRID (default 0):
#   0 = MD_nvec_hybrid.cpp compiles to no-op fallbacks; the ops table is
#       untouched → behavior byte-identical to Config C (or Config D when
#       SHUD_USE_OPENMP_NVECTOR=1 without HYBRID). Default builds
#       (`make shud`, `make shud_omp`) leave this at 0 → CI build-and-
#       compare stays green.
#   1 = -DSHUD_NVEC_HYBRID=1: the serial reduction overrides install on
#       udata/du after N_VNew_OpenMP (before CVodeInit, before the
#       SHUD_NVEC_PROF wrap). Config E leg:
#         make shud_omp SHUD_USE_OPENMP_NVECTOR=1 SHUD_NVEC_HYBRID=1
#
# HYBRID=1 is meaningful ONLY with the OpenMP NVector backend — the
# overrides target N_VNew_OpenMP's ops table. If set without
# SHUD_USE_OPENMP_NVECTOR=1 the build MUST fail LOUD (not silently no-op
# on a Serial vector). The check is parse-time: like the release
# SHUD_ENABLE_OPENMP_RHS default-flip, target-specific variables do NOT
# reach a parse-time `ifeq`, so the guard reads the CLI/`?=` value of
# SHUD_USE_OPENMP_NVECTOR directly (which the Config E leg sets on the
# make CLI). `make shud_omp SHUD_NVEC_HYBRID=1` alone (no explicit
# NVector flag) therefore aborts, which is intended — Config E requires
# the NVector backend to be explicitly requested.
SHUD_NVEC_HYBRID ?= 0
ifeq ($(SHUD_NVEC_HYBRID),0)
  SHUD_NVEC_HYBRID_DEFINE :=
else ifeq ($(SHUD_NVEC_HYBRID),1)
  ifneq ($(SHUD_USE_OPENMP_NVECTOR),1)
$(error SHUD_NVEC_HYBRID=1 requires SHUD_USE_OPENMP_NVECTOR=1 (Config E = OpenMP NVector element-wise + serial reduction overrides); set both, e.g. `make shud_omp SHUD_USE_OPENMP_NVECTOR=1 SHUD_NVEC_HYBRID=1`)
  endif
  SHUD_NVEC_HYBRID_DEFINE := -DSHUD_NVEC_HYBRID=1
else
$(error SHUD_NVEC_HYBRID must be 0 or 1, got '$(SHUD_NVEC_HYBRID)')
endif

# EXTRA_CXXFLAGS: free-form append slot for caller-supplied defines
# the build doesn't otherwise know about. The S1d.1 smoke test uses
# `EXTRA_CXXFLAGS=-DNDEBUG` to verify that `std::abort()` stubs
# remain effective in release-build configuration (assert(false)
# would not). Layer-1/2 disallowed-flag guards above still apply:
# `-Ofast`/`-ffast-math`/`-funsafe-math-optimizations` in
# EXTRA_CXXFLAGS would be caught by the MAKEOVERRIDES scan.
EXTRA_CXXFLAGS ?=

# -----------------------------------------------------------------
# S5d.2-5a (#179) — ASan + UBSan build variant `shud_asan`
# -----------------------------------------------------------------
# Adds `-fsanitize=address,undefined -fno-omit-frame-pointer` on top
# of the locked CXX_BASE_FLAGS. Used to gate the S5d.2 jagged → flat
# refactor: any out-of-bounds write into QeleSurf_flat / QeleSub_flat
# or signed-overflow in the row-major index expression `3*i + j`
# surfaces under runtime sanitizer instead of silently corrupting
# memory. The sanitizer flags are NOT IEEE-754-relevant (they add
# instrumentation, not optimization) so they are kept OUT of the
# DISALLOWED_FLAGS scan; the `shud_asan` target uses them via a
# separate variable so the Layer-1/2 guards keep their original
# vigilance over `-Ofast`/`-ffast-math`.
#
# Wall-clock impact: ASan/UBSan-instrumented runs are 2-5x slower than
# default builds; keep test cases to keliya + qhh under 90-day
# truncation in CI. Use `make shud_asan` then run the binary normally
# (../../shud <case>); ASan output goes to stderr — redirect to
# sanitizer_report.txt for archival.
SHUD_ASAN_FLAGS := -fsanitize=address,undefined -fno-omit-frame-pointer
.PHONY: shud_asan
shud_asan: check_sundials $(MAIN_shud) $(SRC) $(SRC_H)
	@echo '...Compiling shud_asan (ASan + UBSan instrumented; S5d.2-5a #179) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(BUILDDIR)/shud_asan $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_NVEC_OMP_LK)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(BUILDDIR)/shud_asan $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_NVEC_OMP_LK)
	@echo
	@echo " $(BUILDDIR)/shud_asan is compiled successfully (ASan+UBSan)"
	@echo

# -----------------------------------------------------------------
# Optional profile timer instrumentation (openmp issue #10 / S0-8a)
# -----------------------------------------------------------------
# Off by default. Set `SHUD_ENABLE_PROFILE=1` on the make CLI to:
#   - compile in the call sites in shud.cpp that dump
#     `<outpath>/profile_B0.yaml` at end of run, and
#   - pull in the outer-repo profile timer library
#     (../tools/profile/timer.cpp + timer.h) on the compile line.
#
# When PROFILE=0 the shud.cpp call sites collapse to no-ops via the
# `#ifdef SHUD_ENABLE_PROFILE` block + header inline stubs, and the
# timer.cpp source is not added to the build. The PROFILE=0 binary is
# preprocessor-equivalent to a build that omits these files entirely.
# B0 bitwise invariant: keliya `*.dat` SHA256 MUST match between
# PROFILE=0 and PROFILE=1 builds at the same flag set.
#
# #10 ships infrastructure only; actual timer instrumentation inside
# SHUD source (timer Hooks inside MD_f.cpp / RHS core) is deferred
# to S0-10. Until then the YAML output is a skeleton with all bucket
# values = 0.0.
SHUD_ENABLE_PROFILE ?= 0
ifeq ($(SHUD_ENABLE_PROFILE),1)
  SHUD_PROFILE_DEFINE := -DSHUD_ENABLE_PROFILE=1
  # tools/profile lives in the outer repo, sibling to SHUD/. Pull in
  # both the include path (for shud.cpp's `#include "timer.h"`) and
  # the impl source so the build picks up the real timer body.
  SHUD_PROFILE_INC    := -I$(CURDIR)/../tools/profile
  SHUD_PROFILE_SRC    := $(CURDIR)/../tools/profile/timer.cpp
else ifeq ($(SHUD_ENABLE_PROFILE),0)
  SHUD_PROFILE_DEFINE :=
  # PROFILE=0: shud.cpp's #include "timer.h" is itself #ifdef-guarded
  # to SHUD_ENABLE_PROFILE, so the header path can stay empty. This
  # keeps the openmp-baseline branch SHUD checkout self-contained when
  # built standalone (without the outer Hydro-SHUD/openMP repo).
  SHUD_PROFILE_INC    :=
  SHUD_PROFILE_SRC    :=
else
$(error SHUD_ENABLE_PROFILE must be 0 or 1, got '$(SHUD_ENABLE_PROFILE)')
endif

# -----------------------------------------------------------------
# Platform-conditional OpenMP flags
# -----------------------------------------------------------------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LIBOMP_PREFIX     ?= $(shell brew --prefix libomp 2>/dev/null)
  # Guard libomp when building shud_omp (always needs OpenMP) OR when
  # building `shud` with SHUD_USE_OPENMP_NVECTOR=1 (S1d.2 #48 Config D
  # — needs omp_set_num_threads in shud.cpp + libsundials_nvecopenmp's
  # OpenMP runtime) OR when invoking the smoke_configd target (S1d.2
  # #49 — same OpenMP NVector runtime dependency, but bypasses
  # SHUD_USE_OPENMP_NVECTOR on the make CLI by hardcoding the flag
  # in the recipe). Serial `make shud` (Configs A/B/C) is fine without
  # libomp installed.
  # P1e PR-G (#315) — `SHUD_ENABLE_OPENMP_RHS=1` also brings in the
  # libomp runtime dependency (Config C/D need -fopenmp + -lomp at link
  # time; see SHUD_OMP_RHS_CK / SHUD_OMP_RHS_LK above), so extend the
  # libomp guard to fire for Config C builds as well.
  ifneq (,$(filter shud_omp smoke_configd,$(MAKECMDGOALS))$(filter 1,$(SHUD_USE_OPENMP_NVECTOR))$(filter 1,$(SHUD_ENABLE_OPENMP_RHS)))
    ifeq ($(LIBOMP_PREFIX),)
$(error libomp not found via 'brew --prefix libomp'; run 'brew install libomp' before make shud_omp, `make shud SHUD_USE_OPENMP_NVECTOR=1`, `make shud SHUD_ENABLE_OPENMP_RHS=1`, or make smoke_configd)
    endif
  endif
  INC_OMP           ?= $(LIBOMP_PREFIX)/include
  LIB_OMP           ?= $(LIBOMP_PREFIX)/lib
  CXX_OPENMP_CFLAGS = -Xpreprocessor -fopenmp
  CXX_OPENMP_LFLAGS = -L$(LIB_OMP) -lomp
else
  INC_OMP           ?= /usr/include
  LIB_OMP           ?=
  CXX_OPENMP_CFLAGS = -fopenmp
  CXX_OPENMP_LFLAGS = -lgomp
endif

# -----------------------------------------------------------------
# Paths
# -----------------------------------------------------------------
SUNDIALS_DIR ?= $(CURDIR)/InstallSundials

SHELL    = /bin/sh
BUILDDIR = .
SRC_DIR  = src

LIB_SYS  = /usr/local/lib/
LIB_SUN  = $(SUNDIALS_DIR)/lib

INC_MPI  = /usr/local/opt/open-mpi

TARGET_EXEC  = $(BUILDDIR)/shud
TARGET_OMP   = $(BUILDDIR)/shud_omp
TARGET_DEBUG = $(BUILDDIR)/shud_debug

MAIN_shud  = $(SRC_DIR)/main.cpp
MAIN_OMP   = $(SRC_DIR)/main.cpp
MAIN_DEBUG = $(SRC_DIR)/main.cpp

# -----------------------------------------------------------------
# Compiler — let user PATH resolve the toolchain.
# On macOS, `g++` is Apple Clang's wrapper. On Linux, it is GCC (use
# CXX=g++-12 to pin GCC 12 in CI per docs/build_manifest.md).
#
# NOTE: `CXX ?= g++` is a no-op — GNU make defines CXX = g++ as a built-in
# default with origin `default`, so `?=` (which only assigns when the
# variable is undefined) never fires, and reality defaults to `c++`.
# We check the origin explicitly so only the make built-in is overridden;
# any environment / CLI / Makefile-assigned value is preserved.
ifeq ($(origin CXX),default)
  CXX := g++
endif
MPICC ?= mpic++

# S2 capstone (PR-8) — MD_f_omp.cpp (legacy `_omp` RHS receivers) has
# been deleted from the tree. The prior compile-inclusion switch
# (sibling of SHUD_USE_OPENMP_NVECTOR) is retired together with the
# source file.
#
# S4 PR-10 (#154) — `$(SRC_DIR)/ModelData/*.cpp` glob already covers
# `SRC_DIR/ModelData/MD_adjacency.cpp` (new file in this PR). Spec
# `s4-adjacency-topology` Scenario L28 ("Makefile SHALL 把
# MD_adjacency.cpp 加入 SOURCE 列表（与 MD_f.cpp 等同级出现）") is
# satisfied by the wildcard — MD_adjacency.cpp lives in the same
# directory as MD_f.cpp and is compiled in the same translation pass
# (see `make shud` recipe below — SRC expands via wildcard at recipe
# time so the new .cpp ships in the link line automatically). Verified
# at build time: `make -n shud | grep MD_adjacency.cpp` shows the file
# on the compile command line.
SRC = $(SRC_DIR)/classes/*.cpp \
      $(SRC_DIR)/ModelData/*.cpp \
      $(SRC_DIR)/Model/*.cpp \
      $(SRC_DIR)/Equations/*.cpp

SRC_H = $(SRC_DIR)/classes/*.hpp \
        $(SRC_DIR)/ModelData/*.hpp \
        $(SRC_DIR)/Model/*.hpp \
        $(SRC_DIR)/Equations/*.hpp

# -----------------------------------------------------------------
# P8-tune.G0 PR-0 — Hypre / BoomerAMG link chain (additive)
# -----------------------------------------------------------------
# G0 wires the SUNLinSol_Hypre wrapper at
# $(SRC_DIR)/Equations/sunlinsol_hypre.{h,cpp}. The wrapper is
# linked into every `make shud` / `make shud_omp` binary
# (link-always, runtime-opt-in via SHUD_LINSOL=amg env var; see
# cvode_config.cpp factory dispatch). Default builds (SHUD_LINSOL
# unset or "spgmr") make zero Hypre calls and produce bit-identical
# SPGMR output vs the pre-G0 baseline.
#
# Install matrix (env-overridable defaults):
#   macOS brew  : HYPRE_INCDIR=/opt/homebrew/include          HYPRE_LIBDIR=/opt/homebrew/lib
#   Ubuntu apt  : HYPRE_INCDIR=/usr/include/hypre              HYPRE_LIBDIR=/usr/lib/x86_64-linux-gnu
#   Server      : HYPRE_INCDIR=/scratch/frd_muziyao/local/hypre-3.1.0/include
#                 HYPRE_LIBDIR=/scratch/frd_muziyao/local/hypre-3.1.0/lib
#
# Set per platform via `make HYPRE_INCDIR=... HYPRE_LIBDIR=... shud`
# or by exporting the vars in the shell before invoking make.
HYPRE_INCDIR ?= /opt/homebrew/include
HYPRE_LIBDIR ?= /opt/homebrew/lib

# MPI_INCDIR — mpi.h search path. Hypre's HYPRE_utilities.h unconditionally
# #include "mpi.h" so we need to expose the MPI headers at compile time even
# when SHUD itself does not invoke any MPI API.
#   macOS brew  : Hypre is built --without-MPI, so this is unused (but harmless)
#   Ubuntu apt  : libopenmpi-dev ships at /usr/lib/x86_64-linux-gnu/openmpi/include
#   Server      : /usr/lib/x86_64-linux-gnu/openmpi/include (cn-node OpenMPI 4.x)
MPI_INCDIR ?= /usr/lib/x86_64-linux-gnu/openmpi/include

# OPENBLAS_LIBDIR — openblas search path. On macOS brew openblas is
# keg-only (installed under /opt/homebrew/opt/openblas/lib, NOT in the
# default linker search path); on Ubuntu apt libopenblas-dev installs
# under the multiarch /usr/lib/x86_64-linux-gnu/ which IS in the default
# search path so the override is empty there. Default to the Mac brew
# path; CI overrides to empty (Ubuntu uses default linker path).
OPENBLAS_LIBDIR ?= /opt/homebrew/opt/openblas/lib

INCLUDES = -I $(SUNDIALS_DIR)/include \
           -I $(INC_OMP) \
           -I $(SRC_DIR)/Model \
           -I $(SRC_DIR)/ModelData \
           -I $(SRC_DIR)/classes \
           -I $(SRC_DIR)/Equations \
           -I $(HYPRE_INCDIR) \
           $(if $(wildcard $(MPI_INCDIR)/mpi.h),-I $(MPI_INCDIR))

# Use $(if …) so an empty $(LIB_OMP) / $(LIB_SYS) does NOT emit a bare `-L`
# token (which gobbles the next argument and breaks the link line on Linux,
# where LIB_OMP is empty by default).
LIBRARIES = $(if $(LIB_OMP),-L$(LIB_OMP)) \
            -L$(LIB_SUN) \
            $(if $(LIB_SYS),-L$(LIB_SYS))

RPATH = '-Wl,-rpath,$(LIB_SUN)'

# MPI_CXX_LIB — Ubuntu apt libhypre-dev is built with OpenMPI C++ bindings
# enabled, so libHYPRE.so has DT_NEEDED on libmpi_cxx.so.40. On Linux the
# linker requires it on the command line ("DSO missing from command line"
# error). Mac brew openmpi disabled C++ bindings (libmpi_cxx not shipped).
# Auto-detect via wildcard.
MPI_CXX_LIB ?= $(if $(wildcard /usr/lib/x86_64-linux-gnu/libmpi_cxx.so*),-lmpi_cxx)

LK_FLAGS = -lm -lsundials_cvode -lsundials_nvecserial \
           -L$(HYPRE_LIBDIR) -lHYPRE -lmpi $(MPI_CXX_LIB) \
           $(if $(OPENBLAS_LIBDIR),-L$(OPENBLAS_LIBDIR)) -lopenblas \
           -Wl,-rpath,$(HYPRE_LIBDIR)
# S1d.2 (openMP #48) — LK_OMP now only carries the platform OpenMP
# runtime flags (`-lgomp` on Linux, `-Xpreprocessor -fopenmp -lomp`
# on macOS via brew libomp). The historical hardcoded
# `-lsundials_nvecopenmp` link has been moved into the conditional
# `$(SHUD_NVEC_OMP_LK)` (driven by SHUD_USE_OPENMP_NVECTOR), so:
#   - serial `make shud` (Config A, default) — no nvecopenmp link
#   - `make shud SHUD_USE_OPENMP_NVECTOR=1` — adds -lsundials_nvecopenmp
#   - `make shud_omp` — implicitly sets SHUD_USE_OPENMP_NVECTOR=1
#     inside the recipe (back-compat with the pre-#48 shud_omp target,
#     which historically always pulled in nvecopenmp).
LK_OMP   = $(CXX_OPENMP_LFLAGS)
LK_DYLN  = "LD_LIBRARY_PATH=$(LIB_SUN)"

# -----------------------------------------------------------------
# SUNDIALS version + install-completeness guard
# -----------------------------------------------------------------
# - MAJOR pinned with `-Eq '^…6$'` (anchored regex) so substrings like
#   60 / 600 / 6X cannot pass as "6".
# - MINOR also pinned: 6.0.x is the supported series; 6.1+ is rejected.
# - PATCH unenforced: future 6.0.x patches are acceptable.
# - Library stat catches the "header present but libs missing" partial
#   install that the old grep-only guard silently allowed.
SUNDIALS_CFG_H = $(SUNDIALS_DIR)/include/sundials/sundials_config.h

.PHONY: check_sundials check_sundials_omp
check_sundials:
	@test -f $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: $(SUNDIALS_CFG_H) not found; run ./configure first"; exit 2)
	@grep -Eq '^#define SUNDIALS_VERSION_MAJOR 6$$' $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: SUNDIALS major != 6 in $(SUNDIALS_CFG_H); B0 requires exactly 6.x"; exit 2)
	@grep -Eq '^#define SUNDIALS_VERSION_MINOR 0$$' $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: SUNDIALS minor != 0; B0 requires 6.0.x; re-run ./configure to pin to 6.0.0"; exit 2)
	@ls $(SUNDIALS_DIR)/lib/libsundials_cvode.* >/dev/null 2>&1 || \
	  (echo "ERROR: libsundials_cvode.* not found under $(SUNDIALS_DIR)/lib; SUNDIALS install is incomplete; re-run ./configure"; exit 2)
	@ls $(SUNDIALS_DIR)/lib/libsundials_nvecserial.* >/dev/null 2>&1 || \
	  (echo "ERROR: libsundials_nvecserial.* not found under $(SUNDIALS_DIR)/lib; SUNDIALS install is incomplete; re-run ./configure"; exit 2)

# OpenMP build additionally needs the nvecopenmp variant.
check_sundials_omp: check_sundials
	@ls $(SUNDIALS_DIR)/lib/libsundials_nvecopenmp.* >/dev/null 2>&1 || \
	  (echo "ERROR: libsundials_nvecopenmp.* not found under $(SUNDIALS_DIR)/lib; SUNDIALS OpenMP variant missing; re-run ./configure"; exit 2)

# -----------------------------------------------------------------
# Targets
# -----------------------------------------------------------------
.PHONY: all check help cvode CVODE clean

all:
	$(MAKE) clean
	$(MAKE) shud
	@echo

check:
	ls $(SUNDIALS_DIR)
	ls $(SUNDIALS_DIR)/lib
	./shud
	@echo

help:
	@echo
	@echo "Usage:"
	@echo "       make all          - clean then build shud (serial)"
	@echo "       make cvode        - install SUNDIALS/CVODE 6.x to ./InstallSundials"
	@echo "       make shud         - build serial shud executable"
	@echo "       make shud_omp     - build OpenMP shud_omp executable"
	@echo "       make shud SHUD_DUMP_RHS=1     - serial build with RHS snapshot hooks compiled in"
	@echo "       make shud_omp SHUD_DUMP_RHS=1 - OpenMP build with RHS snapshot hooks compiled in"
	@echo "       make shud SHUD_ENABLE_PROFILE=1     - serial build with wall-clock profile timer compiled in"
	@echo "       make shud_omp SHUD_ENABLE_PROFILE=1 - OpenMP build with wall-clock profile timer compiled in"
	@echo "       make shud SHUD_ENABLE_OPENMP_RHS=1  - compile in StrictOMP/ProductionOMP std::abort stubs (smoke only; openMP #47)"
	@echo "       make shud SHUD_USE_OPENMP_NVECTOR=1 - serial build with OpenMP N_Vector backend (Config D dim; openMP #48)"
	@echo "       make shud EXTRA_CXXFLAGS=-DSHUD_ENABLE_DIAGNOSTICS - serial build with S5c CVODE diagnostic keys (hlast/qlast) added to cvode_stats.txt (S5c-A #173)"
	@echo "       make shud_asan    - serial build with -fsanitize=address,undefined (S5d.2-5a #179)"
	@echo "       make smoke_strictomp                - build + run StrictOMP SIGABRT regression smoke test (openMP #47)"
	@echo "       make smoke_configd                  - build + run Config D OpenMP NVector runtime probe (openMP #49)"
	@echo "       make check_sundials - verify SUNDIALS 6.x install"
	@echo "       make clean        - remove binary outputs (preserves InstallSundials)"
	@echo
	@echo "P8-tune.G0 Hypre/BoomerAMG link chain (additive, runtime-opt-in via SHUD_LINSOL=amg):"
	@echo "  macOS brew:  HYPRE_INCDIR=/opt/homebrew/include HYPRE_LIBDIR=/opt/homebrew/lib (defaults)"
	@echo "  Ubuntu apt:  HYPRE_INCDIR=/usr/include/hypre    HYPRE_LIBDIR=/usr/lib/x86_64-linux-gnu"
	@echo "  Server:      HYPRE_INCDIR=/scratch/frd_muziyao/local/hypre-3.1.0/include"
	@echo "               HYPRE_LIBDIR=/scratch/frd_muziyao/local/hypre-3.1.0/lib"
	@echo "  LK_FLAGS:    -lHYPRE -lmpi -lopenblas (openblas required by Hypre at link time)"
	@echo "  ColPack:     NOT required for G0 (per PRE0_SPIKE_NOTES.md §1.7)"
	@echo

cvode CVODE:
	@echo '...Install SUNDIALS/CVODE 6.0.0 ...'
	chmod +x configure
	./configure
	@echo

# S1d.2 (openMP #48) — `shud` recipe carries the remaining feature
# defines (SHUD_NVEC_OMP_DEFINE / SHUD_NVEC_OMP_LK), so Config C/D can
# build through the same target via `make shud SHUD_USE_OPENMP_NVECTOR=1`
# / `make shud SHUD_ENABLE_OPENMP_RHS=1` (or combinations). Config A =
# `make shud` (all flags default 0). The S2 capstone (PR-8) retired
# Config B (legacy `_omp` RHS receivers) together with its source file.
# S1d.2 (openMP #48) — pick check_sundials_omp when the user enables
# the OpenMP NVector backend (Config D needs libsundials_nvecopenmp);
# otherwise the basic check is enough (Configs A/B/C don't link nvecopenmp).
ifeq ($(SHUD_USE_OPENMP_NVECTOR),1)
  SHUD_CHECK_TARGET := check_sundials_omp
else
  SHUD_CHECK_TARGET := check_sundials
endif
shud SHUD: $(SHUD_CHECK_TARGET) $(MAIN_shud) $(SRC) $(SRC_H)
	@echo '...Compiling shud (B0 serial / Config A by default) ...'
	@echo  $(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_OMP_RHS_CK) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(SHUD_OMP_RHS_LK)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_OMP_RHS_CK) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(SHUD_OMP_RHS_LK)
	@echo
	@echo " $(TARGET_EXEC) is compiled successfully!"
	@echo

# Release v1.0 — `shud_omp` produces the production Config C binary by
# default: Serial NVector + StrictOMP RHS (ADR-0002 Path 1 winner,
# heihe_x4 sp@8 = 1.6–1.7×). The MAKECMDGOALS-conditional default flip
# above (around SHUD_ENABLE_OPENMP_RHS ?= 1 for shud_omp) means users
# get the parallel RHS without passing extra flags. Researcher escape
# hatches:
#   - `make shud_omp SHUD_ENABLE_OPENMP_RHS=0` → Config A/B (serial RHS)
#     for A/B/D reproducibility (P1c/d era build).
#   - `make shud_omp SHUD_USE_OPENMP_NVECTOR=1` → Config D (both OMP);
#     `SHUD_NVEC_OMP_DEFINE` + `SHUD_NVEC_OMP_LK` handle the link line.
# The recipe adds:
#   - $(CXX_OPENMP_CFLAGS)   : -fopenmp (sets _OPENMP compiler builtin)
#   - $(CXX_OPENMP_LFLAGS)   : -lgomp / -lomp (OpenMP runtime)
#   - $(SHUD_NVEC_OMP_DEFINE): -DSHUD_USE_OPENMP_NVECTOR=1 iff opted in
#   - $(SHUD_NVEC_OMP_LK)    : -lsundials_nvecopenmp iff opted in
# Historical note: pre-release, `shud_omp` hardcoded
# `-DSHUD_USE_OPENMP_NVECTOR=1` + `-lsundials_nvecopenmp` (P1c/d era
# Config B semantics). Release v1.0 flip aligns the default OpenMP
# target with the P1e-endorsed production build.
shud_omp: check_sundials_omp $(MAIN_OMP) $(SRC) $(SRC_H)
	@echo '...Compiling shud_OpenMP (Config C default: Serial NVec + StrictOMP RHS) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(LK_OMP)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(LK_OMP)
	@echo
	@echo " $(TARGET_OMP) is compiled successfully!"
	@echo

# -----------------------------------------------------------------
# S1d.1 smoke test — StrictOMP std::abort regression guard
# -----------------------------------------------------------------
# Builds + runs `tests/s1d_strictomp_assert_smoke.cpp`, which forks a
# child that invokes `Model_Data::rhs_core(..., ExecPolicy::StrictOMP)`
# under `-DNDEBUG` and asserts via `waitpid` that the child died by
# SIGABRT. This guards the decision to use `std::abort()` rather than
# `assert(false)` for the OMP-policy stubs (assert would be stripped
# under -DNDEBUG and the case would silently fall through to the
# next statement). Requires SHUD_ENABLE_OPENMP_RHS=1 to make the
# OMP cases compile-visible.
#
# The smoke binary supplies its own main(), so we must drop SHUD's
# main.cpp from the link line. Wildcard expansion happens at recipe
# time (SRC uses src/.../*.cpp globs), so the substring filter works
# only on the expanded list — wrap in $(filter-out ...) accordingly.
SHUD_SRC_NOMAIN := $(filter-out $(SRC_DIR)/main.cpp,$(wildcard $(SRC)))

.PHONY: smoke_strictomp
smoke_strictomp: check_sundials tests/s1d_strictomp_assert_smoke.cpp $(SRC) $(SRC_H)
	@echo '...Compiling s1d_strictomp_assert_smoke ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) -DNDEBUG -DSHUD_ENABLE_OPENMP_RHS=1 $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/s1d_strictomp_smoke tests/s1d_strictomp_assert_smoke.cpp $(SHUD_SRC_NOMAIN) $(LK_FLAGS)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) -DNDEBUG -DSHUD_ENABLE_OPENMP_RHS=1 $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/s1d_strictomp_smoke tests/s1d_strictomp_assert_smoke.cpp $(SHUD_SRC_NOMAIN) $(LK_FLAGS)
	@echo
	@echo '...Running s1d_strictomp_assert_smoke (expect SIGABRT in child) ...'
	@DYLD_LIBRARY_PATH=$(LIB_SUN) tests/s1d_strictomp_smoke && echo 'OK: child SIGABRT observed' || (echo 'FAIL: child did not SIGABRT'; exit 1)

# -----------------------------------------------------------------
# S1d.2 smoke test — Config D OpenMP N_Vector runtime probe (openMP #49)
# -----------------------------------------------------------------
# Builds + runs `tests/s1d_configd_nvec_smoke.cpp`, a standalone (no
# SHUD framework link) probe that:
#   1. constructs an OpenMP-backed N_Vector via N_VNew_OpenMP, and
#   2. asserts N_VGetVectorID == SUNDIALS_NVEC_OPENMP, and
#   3. exercises the generic N_VDestroy dispatch path (the §4.19
#      fix landed in #48 -- pre-fix SHUD used N_VDestroy_Serial on
#      an OpenMP vector, type-tag mismatch UB).
#
# Requires both:
#   - SHUD_USE_OPENMP_NVECTOR=1  (compile-line + nvecopenmp link;
#                                  also triggers the libomp guard
#                                  via the existing MAKECMDGOALS
#                                  filter below in this Makefile)
#   - SHUD_ENABLE_OPENMP_RHS=1   (Config D contract; symbolic --
#                                  this test does not touch rhs_core)
#
# We deliberately do NOT link $(SHUD_SRC_NOMAIN) and do NOT pull in
# $(LK_FLAGS) (which carries libsundials_cvode + libsundials_nvecserial):
# the test only exercises the OpenMP NVector + SUNContext API surface.
# `-lsundials_nvecopenmp` is bracketed by `-lsundials_generic` for
# SUNContext_Create / SUNContext_Free (nvecopenmp links generic at its
# own link time, but our standalone executable needs the generic
# symbols visible at link time too).
.PHONY: smoke_configd
smoke_configd: check_sundials_omp tests/s1d_configd_nvec_smoke.cpp
	@echo '...Compiling s1d_configd_nvec_smoke (Config D: OpenMP NVector runtime probe) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) -DSHUD_ENABLE_OPENMP_RHS=1 -DSHUD_USE_OPENMP_NVECTOR=1 $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/s1d_configd_nvec_smoke tests/s1d_configd_nvec_smoke.cpp -lsundials_nvecopenmp -lsundials_generic $(CXX_OPENMP_LFLAGS)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) -DSHUD_ENABLE_OPENMP_RHS=1 -DSHUD_USE_OPENMP_NVECTOR=1 $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/s1d_configd_nvec_smoke tests/s1d_configd_nvec_smoke.cpp -lsundials_nvecopenmp -lsundials_generic $(CXX_OPENMP_LFLAGS)
	@echo
	@echo '...Running s1d_configd_nvec_smoke (expect OK: N_VGetVectorID == SUNDIALS_NVEC_OPENMP) ...'
	@DYLD_LIBRARY_PATH=$(LIB_SUN) tests/s1d_configd_nvec_smoke || (echo 'FAIL: Config D NVector smoke did not print OK'; exit 1)

# -----------------------------------------------------------------
# S4 PR-10 (#154) — adjacency fallback unit test
# -----------------------------------------------------------------
# Builds tests/test_adjacency_fallback.cpp, which constructs a synthetic
# Model_Data mock whose entity `.index` fields violate the
# `index == array_index + 1` invariant. Verifies build_adjacency_lists()
# sets `adjacency_fallback_triggered = true` AND constructs lists in
# array-index order (NOT id-sort) — per spec s4-adjacency-topology
# Scenario "三条 assert 在所有 6 case 都 pass + fallback 单测 PASS".
#
# Test binary supplies its own main(); reuse the same `SHUD_SRC_NOMAIN`
# var from smoke_strictomp (Makefile L484) so we link in the rest of the
# SHUD framework objects.
.PHONY: test_adjacency_fallback
test_adjacency_fallback: check_sundials tests/test_adjacency_fallback.cpp $(SRC) $(SRC_H)
	@echo '...Compiling test_adjacency_fallback (S4 PR-10 fallback unit test) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/test_adjacency_fallback tests/test_adjacency_fallback.cpp $(SHUD_SRC_NOMAIN) $(LK_FLAGS)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/test_adjacency_fallback tests/test_adjacency_fallback.cpp $(SHUD_SRC_NOMAIN) $(LK_FLAGS)
	@echo
	@echo '...Running test_adjacency_fallback (expect PASS) ...'
	@DYLD_LIBRARY_PATH=$(LIB_SUN) tests/test_adjacency_fallback && echo 'OK: adjacency fallback unit test PASS' || (echo 'FAIL: adjacency fallback unit test'; exit 1)

# -----------------------------------------------------------------
# P8-tune.D KLU spike — additive carve-out (openspec change p8tune-klu-spike)
# -----------------------------------------------------------------
# `libshud.a` is the documented carve-out per spec
# `klu-pattern-spike-verdict` REQ-1 Scenario "Tool authoring with no
# SHUD source patch": the spike tool under `tools/p8tune.D/` links this
# archive instead of `main.cpp`, so the spike binary has the Model_Data
# API (loadinput / initialize / rhs_core / public AoS) without dragging
# in `shud.cpp::SHUD()` CVODE main loop. ZERO impact to existing `shud`
# / `shud_omp` / `shud_asan` / `smoke_*` / `test_adjacency_fallback`
# targets — uses the same `SHUD_SRC_NOMAIN` wildcard already defined at
# L548 for the s1d StrictOMP smoke test (proves the SHUD framework
# objects compile + link without `main.cpp`).
#
# Object files land under `_libshud_obj/` (gitignored via SHUD repo's
# existing `*.o` ignore + this directory's gitignored status under
# SHUD-OpenMP outer-repo `.gitignore` for `SHUD/_libshud_obj/`). The
# archive lands at SHUD-repo root next to the `shud` binary so the
# spike tool's `tools/p8tune.D/Makefile` can reference it via the
# relative path `../../SHUD/libshud.a`.
#
# Flag set: ZERO instrumentation defines (no SHUD_DUMP_RHS / no
# SHUD_ENABLE_PROFILE / no OMP defines / no DIAGNOSTICS) so the
# archive is a clean Config-A-equivalent object set. SUNDIALS includes
# are pulled in via $(INCLUDES) because `Model_Data.hpp` and friends
# include `nvector_serial.h` etc.
LIBSHUD_OBJ_DIR  := _libshud_obj
LIBSHUD_ARCHIVE  := libshud.a
LIBSHUD_SRC      := $(SHUD_SRC_NOMAIN)
LIBSHUD_OBJ      := $(patsubst $(SRC_DIR)/%.cpp,$(LIBSHUD_OBJ_DIR)/%.o,$(LIBSHUD_SRC))

$(LIBSHUD_OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: libshud.a
libshud.a: check_sundials $(LIBSHUD_OBJ)
	@echo '...Archiving libshud.a (P8-tune.D KLU spike — openspec change p8tune-klu-spike) ...'
	$(AR) rcs $(LIBSHUD_ARCHIVE) $(LIBSHUD_OBJ)
	@echo " $(LIBSHUD_ARCHIVE) is archived successfully (additive — Config A flag set, no instrumentation)"
	@echo

clean:
	@echo "Cleaning ... "
	@echo "  rm -f *.o"
	@rm -f *.o
	@echo "  rm -f $(TARGET_EXEC)"
	@rm -f $(TARGET_EXEC)
	@echo "  rm -f $(TARGET_OMP)"
	@rm -f $(TARGET_OMP)
	@echo "  rm -f $(TARGET_DEBUG)"
	@rm -f $(TARGET_DEBUG)
	@echo "  rm -f $(BUILDDIR)/shud_asan"
	@rm -f $(BUILDDIR)/shud_asan
	@echo "  rm -f $(LIBSHUD_ARCHIVE) + rm -rf $(LIBSHUD_OBJ_DIR)/"
	@rm -f $(LIBSHUD_ARCHIVE)
	@rm -rf $(LIBSHUD_OBJ_DIR)
	@echo "Done. (InstallSundials/ preserved)"
	@echo
