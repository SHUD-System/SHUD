# -----------------------------------------------------------------
# Makefile for SHUD
# -----------------------------------------------------------------
# Programmer: Lele Shu (lele.shu@gmail.com)
# SHUD model is a heritage of Penn State Integrated Hydrologic Model (PIHM).
# -----------------------------------------------------------------
# Prerequisites:
#   - SUNDIALS 6.x installed at $(SUNDIALS_DIR) (./configure installs 6.0.0)
#   - For OpenMP on macOS: `brew install libomp`
#   - For OpenMP on Linux: GCC with libgomp
# `make help` lists the targets and options. OpenMP_Guide.md explains
# the OpenMP builds.
# -----------------------------------------------------------------

# A plain `make` builds the serial model.
.DEFAULT_GOAL := all

# -----------------------------------------------------------------
# Compiler flags (fixed)
# -----------------------------------------------------------------
# The flag set is fixed because the results of the model must be
# reproducible bit for bit: -ffp-contract=off and -fno-fast-math keep
# strict IEEE-754 arithmetic. The bit-identity between the default
# `make shud_omp` build and the serial N_Vector build also depends on
# these flags (OpenMP_NVector_Determinism.md).
# `override ... :=` makes make ignore a command-line assignment such as
# `make CXX_BASE_FLAGS=...`. SHUD_BUILD_CFLAGS is the name used in the
# recipes; both variables must stay `override`-protected.
override CXX_BASE_FLAGS    := -O2 -g -ffp-contract=off -fno-fast-math -std=c++14
override SHUD_BUILD_CFLAGS := $(CXX_BASE_FLAGS)
# The OpenMP features are selected by two independent options, defined
# below:
#   SHUD_ENABLE_OPENMP_RHS   parallel evaluation of the right-hand side
#   SHUD_USE_OPENMP_NVECTOR  OpenMP N_Vector inside CVODE

# CFLAGS is kept only for tools that read $(CFLAGS). The recipes use
# $(SHUD_BUILD_CFLAGS), so `make CFLAGS=...` cannot change the flags.
CFLAGS            = $(CXX_BASE_FLAGS)

# -----------------------------------------------------------------
# Rejected flags
# -----------------------------------------------------------------
# Flags that break IEEE-754 arithmetic are rejected with an error instead
# of being silently ignored, so that nobody believes a build used them.
#
# Check 1 scans the usual flag variables word by word; it catches e.g.
# `CXXFLAGS=-ffast-math`.
#
# Check 2 (further down) looks for a command-line assignment of exactly
# such a flag, e.g. `make SHUD_BUILD_CFLAGS=-Ofast`, which the
# `override :=` above would drop without a message. It matches whole
# `VAR=<flag>` words of $(MAKEOVERRIDES); a substring search would wrongly
# reject a path such as `SUNDIALS_DIR=/opt/sundials-Ofast-tuned`.
#
# DISALLOWED_FLAGS is `override`-protected so that
# `make DISALLOWED_FLAGS=` cannot switch the checks off.
override DISALLOWED_FLAGS := -ffast-math -Ofast -funsafe-math-optimizations
ifneq (,$(filter $(DISALLOWED_FLAGS),$(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(MAKEOVERRIDES) $(MAKEFLAGS) $(SHUD_BUILD_CFLAGS) $(CXX_BASE_FLAGS)))
$(error disallowed flag detected (one of $(DISALLOWED_FLAGS)) in CFLAGS/CXXFLAGS/CPPFLAGS/LDFLAGS/MAKEOVERRIDES/MAKEFLAGS/SHUD_BUILD_CFLAGS/CXX_BASE_FLAGS; SHUD requires strict IEEE-754 arithmetic for reproducible results)
endif
# Check 2, described above.
LAYER2_HITS := $(strip $(foreach tok,$(MAKEOVERRIDES),$(foreach f,$(DISALLOWED_FLAGS),$(if $(filter %=$(f),$(tok)),$(f)))))
ifneq (,$(LAYER2_HITS))
$(error disallowed flag detected ($(LAYER2_HITS)) in a make-CLI assignment (MAKEOVERRIDES=[$(MAKEOVERRIDES)]); attempts to inject via SHUD_BUILD_CFLAGS / CXX_BASE_FLAGS / etc are also rejected; SHUD requires strict IEEE-754 arithmetic for reproducible results)
endif

# -----------------------------------------------------------------
# SHUD_DUMP_RHS — snapshots of the right-hand side (debugging aid)
# -----------------------------------------------------------------
# 0 (default): the dump code is not compiled.
# 1: snapshots of the RHS state are written at chosen times, controlled
#    at run time by environment variables (src/ModelData/MD_rhs_dump.cpp):
#      SHUD_DUMP_OUTPUT_DIR  directory of the snapshot files (default: cwd)
#      SHUD_DUMP_T_VALUES    times at which to write
#      SHUD_DUMP_T_TOL, SHUD_DUMP_CASE_ID, SHUD_DUMP_SITE,
#      SHUD_DUMP_FNAME_SUFFIX
# The model output is the same with and without this option.
SHUD_DUMP_RHS ?= 0
ifeq ($(SHUD_DUMP_RHS),1)
  SHUD_DUMP_DEFINE := -DSHUD_DUMP_RHS=1
else ifeq ($(SHUD_DUMP_RHS),0)
  SHUD_DUMP_DEFINE :=
else
$(error SHUD_DUMP_RHS must be 0 or 1, got '$(SHUD_DUMP_RHS)')
endif

# -----------------------------------------------------------------
# SHUD_ENABLE_OPENMP_RHS — parallel evaluation of the right-hand side
# -----------------------------------------------------------------
# 0: the OpenMP code of the RHS (src/Model/MD_rhs_core.cpp) is not
#    compiled. Default for `make shud`.
# 1: the RHS (fluxes of elements, river segments and lakes) runs in an
#    OpenMP parallel region. Default for `make shud_omp`. The number of
#    threads is taken from the environment variable SHUD_RHS_THREADS
#    (src/Model/shud.cpp).
# A value given on the command line always wins.
#
# The default depends on the goal, so it is chosen with MAKECMDGOALS. A
# target-specific variable would not work: the `ifeq` blocks below are
# evaluated while the Makefile is parsed, before any target-specific
# value exists.
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
  # The parallel RHS needs the OpenMP compiler and linker flags. `=`, not
  # `:=`, is required here: CXX_OPENMP_CFLAGS / CXX_OPENMP_LFLAGS are
  # defined further down, and `:=` would capture empty values and drop
  # -fopenmp and the OpenMP library from the build.
  SHUD_OMP_RHS_CK      = $(CXX_OPENMP_CFLAGS)
  SHUD_OMP_RHS_LK      = $(CXX_OPENMP_LFLAGS)
else
$(error SHUD_ENABLE_OPENMP_RHS must be 0 or 1, got '$(SHUD_ENABLE_OPENMP_RHS)')
endif

# -----------------------------------------------------------------
# SHUD_USE_OPENMP_NVECTOR — N_Vector implementation used by CVODE
# -----------------------------------------------------------------
# 0: serial N_Vector (N_VNew_Serial); libsundials_nvecopenmp is not
#    linked. Default for every goal except shud_omp.
# 1: OpenMP N_Vector (N_VNew_OpenMP); the number of threads is NUM_OPENMP
#    of the project's .cfg.para file. Links libsundials_nvecopenmp.
#    Default for `make shud_omp`.
# Independent of SHUD_ENABLE_OPENMP_RHS.
#
# Together with SHUD_NVEC_HYBRID=1 (also the default for shud_omp, see
# below) the output is bit-identical to the serial N_Vector at every
# thread count. `make shud_omp SHUD_USE_OPENMP_NVECTOR=0` builds the
# serial N_Vector with the parallel RHS.
#
# The goal-dependent default uses MAKECMDGOALS for the reason given
# above. Note that `make shud shud_omp` in one command applies the
# shud_omp defaults to both targets.
ifneq (,$(filter shud_omp,$(MAKECMDGOALS)))
  SHUD_USE_OPENMP_NVECTOR ?= 1
endif
SHUD_USE_OPENMP_NVECTOR ?= 0
ifeq ($(SHUD_USE_OPENMP_NVECTOR),0)
  SHUD_NVEC_OMP_DEFINE :=
  SHUD_NVEC_OMP_CK     :=
  SHUD_NVEC_OMP_LK     :=
else ifeq ($(SHUD_USE_OPENMP_NVECTOR),1)
  SHUD_NVEC_OMP_DEFINE := -DSHUD_USE_OPENMP_NVECTOR=1
  # The OpenMP N_Vector needs the OpenMP runtime, so the compiler and
  # linker flags are added here and `make shud SHUD_USE_OPENMP_NVECTOR=1`
  # works without a separate -fopenmp. `=`, not `:=`, for the same reason
  # as for SHUD_OMP_RHS_CK above.
  SHUD_NVEC_OMP_CK      = $(CXX_OPENMP_CFLAGS)
  SHUD_NVEC_OMP_LK      = $(CXX_OPENMP_LFLAGS) -lsundials_nvecopenmp
else
$(error SHUD_USE_OPENMP_NVECTOR must be 0 or 1, got '$(SHUD_USE_OPENMP_NVECTOR)')
endif

# -----------------------------------------------------------------
# SHUD_NVEC_HYBRID — order-preserving sums on the OpenMP N_Vector
# -----------------------------------------------------------------
# 0: the OpenMP N_Vector is used as SUNDIALS ships it. Its sums and norms
#    combine per-thread partial results, so the model output changes with
#    the number of threads. Not for production; see the guard below.
# 1: element-wise operations stay parallel; every accumulating operation
#    (dot product, norms, minimum, ...) is replaced by a serial loop
#    (src/Model/MD_nvec_hybrid.cpp). The output is bit-identical at every
#    thread count and to the serial N_Vector. The binary reports
#    "NVEC config: Config E" at startup. Default for `make shud_omp`.
# Requires SHUD_USE_OPENMP_NVECTOR=1; otherwise the build stops with an
# error instead of silently doing nothing on a serial vector.
# Background: OpenMP_NVector_Determinism.md.
#
# The default 1 applies only to the shud_omp goal and only when
# SHUD_USE_OPENMP_NVECTOR=1, so that
# `make shud_omp SHUD_USE_OPENMP_NVECTOR=0` leaves it at 0 and does not
# trip the requirement above.
ifneq (,$(filter shud_omp,$(MAKECMDGOALS)))
  ifeq ($(SHUD_USE_OPENMP_NVECTOR),1)
    SHUD_NVEC_HYBRID ?= 1
  endif
endif
SHUD_NVEC_HYBRID ?= 0
# Guard: `make shud_omp SHUD_NVEC_HYBRID=0` would build the OpenMP
# N_Vector without the order-preserving sums, whose output depends on the
# number of threads. It is refused unless SHUD_ALLOW_CONFIG_D=1 is given.
# The guard is limited to the shud_omp goal; the `smoke_configd` target
# sets its flags in its own recipe and does not pass through here.
ifneq (,$(filter shud_omp,$(MAKECMDGOALS)))
  ifeq ($(SHUD_USE_OPENMP_NVECTOR),1)
    ifeq ($(SHUD_NVEC_HYBRID),0)
      ifneq ($(SHUD_ALLOW_CONFIG_D),1)
$(error SHUD_NVEC_HYBRID=0 with the OpenMP N_Vector gives results that depend on the number of threads and is refused. Use SHUD_USE_OPENMP_NVECTOR=0 for the serial N_Vector, or add SHUD_ALLOW_CONFIG_D=1 to build it anyway (testing only))
      endif
    endif
  endif
endif
ifeq ($(SHUD_NVEC_HYBRID),0)
  SHUD_NVEC_HYBRID_DEFINE :=
else ifeq ($(SHUD_NVEC_HYBRID),1)
  ifneq ($(SHUD_USE_OPENMP_NVECTOR),1)
$(error SHUD_NVEC_HYBRID=1 requires SHUD_USE_OPENMP_NVECTOR=1 (the order-preserving sums are installed on the OpenMP N_Vector))
  endif
  SHUD_NVEC_HYBRID_DEFINE := -DSHUD_NVEC_HYBRID=1
else
$(error SHUD_NVEC_HYBRID must be 0 or 1, got '$(SHUD_NVEC_HYBRID)')
endif

# -----------------------------------------------------------------
# SHUD_NVEC_DETRED — parallel sums in a fixed tree
# -----------------------------------------------------------------
# 0 (default): the sums are the serial loops of SHUD_NVEC_HYBRID.
# 1: the summations run in parallel over fixed blocks of
#    SHUD_NVEC_DETRED_B entries (4096 by default, independent of the
#    number of threads). Each block is summed serially in index order and
#    the block results are combined in a fixed binary tree. The output is
#    bit-identical at every thread count, but the summation order differs
#    from DETRED=0, so the output is not bit-identical to the other
#    builds. The binary reports "NVEC config: Config E2" at startup.
#      make shud_omp SHUD_NVEC_DETRED=1
# Requires SHUD_NVEC_HYBRID=1 (and therefore SHUD_USE_OPENMP_NVECTOR=1);
# otherwise the build stops with an error.
#
# SHUD_NVEC_DETRED_B and SHUD_NVEC_DETRED_NEUMAIER (compensated
# summation, 0 by default) are compile-time values passed through
# EXTRA_CXXFLAGS, e.g. `EXTRA_CXXFLAGS=-DSHUD_NVEC_DETRED_B=256`.
SHUD_NVEC_DETRED ?= 0
ifeq ($(SHUD_NVEC_DETRED),0)
  SHUD_NVEC_DETRED_DEFINE :=
else ifeq ($(SHUD_NVEC_DETRED),1)
  ifneq ($(SHUD_NVEC_HYBRID),1)
$(error SHUD_NVEC_DETRED=1 requires SHUD_NVEC_HYBRID=1 (the fixed-tree sums replace the serial sums of SHUD_NVEC_HYBRID); use `make shud_omp SHUD_NVEC_DETRED=1`)
  endif
  SHUD_NVEC_DETRED_DEFINE := -DSHUD_NVEC_DETRED=1
else
$(error SHUD_NVEC_DETRED must be 0 or 1, got '$(SHUD_NVEC_DETRED)')
endif

# EXTRA_CXXFLAGS: additional compiler options, e.g. `-DNDEBUG` or
# `-DSHUD_ENABLE_DIAGNOSTICS`. The rejected-flags checks above apply to
# it as well.
EXTRA_CXXFLAGS ?=

# -----------------------------------------------------------------
# shud_asan — build with AddressSanitizer + UndefinedBehaviorSanitizer
# -----------------------------------------------------------------
# Adds `-fsanitize=address,undefined -fno-omit-frame-pointer` to the fixed
# flags. The sanitizers add checks and do not change floating-point
# optimization, so they are not in the rejected-flags list. Runs are 2-5
# times slower; use short simulations. Reports go to stderr.
SHUD_ASAN_FLAGS := -fsanitize=address,undefined -fno-omit-frame-pointer
.PHONY: shud_asan
shud_asan: check_sundials $(MAIN_shud) $(SRC) $(SRC_H)
	@echo '...Compiling shud_asan (ASan + UBSan instrumented) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_DETRED_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(BUILDDIR)/shud_asan $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_NVEC_OMP_LK)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_DETRED_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(BUILDDIR)/shud_asan $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_ASAN_FLAGS) $(SHUD_NVEC_OMP_LK)
	@echo
	@echo " $(BUILDDIR)/shud_asan is compiled successfully (ASan+UBSan)"
	@echo

# -----------------------------------------------------------------
# SHUD_ENABLE_PROFILE — wall-clock timers
# -----------------------------------------------------------------
# 0 (default): the timer calls are not compiled.
# 1: compiles tools/profile/timer.cpp and writes
#    `<output dir>/profile_B0.yaml` at the end of the run, with the time
#    spent in the RHS, in CVODE, in forcing input, in ET and in output.
# The model output is the same with and without this option.
SHUD_ENABLE_PROFILE ?= 0
ifeq ($(SHUD_ENABLE_PROFILE),1)
  SHUD_PROFILE_DEFINE := -DSHUD_ENABLE_PROFILE=1
  # Include path of "timer.h" and the implementation file.
  SHUD_PROFILE_INC    := -I$(CURDIR)/tools/profile
  SHUD_PROFILE_SRC    := $(CURDIR)/tools/profile/timer.cpp
else ifeq ($(SHUD_ENABLE_PROFILE),0)
  SHUD_PROFILE_DEFINE :=
  # No include path is needed: every `#include "timer.h"` in the sources
  # is inside `#ifdef SHUD_ENABLE_PROFILE`.
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
  # libomp is needed by shud_omp, by smoke_configd, and by `make shud`
  # with SHUD_USE_OPENMP_NVECTOR=1 or SHUD_ENABLE_OPENMP_RHS=1. A plain
  # `make shud` builds without it.
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
# Compiler. On macOS `g++` is Apple clang; on Linux it is GCC. Use
# `make CXX=...` to choose another one.
#
# `CXX ?= g++` would have no effect: GNU make predefines CXX (origin
# `default`), so `?=` never assigns. Testing the origin replaces only the
# built-in value and keeps a value from the environment or the command
# line.
ifeq ($(origin CXX),default)
  CXX := g++
endif
MPICC ?= mpic++

# The source lists are wildcards: a new .cpp file in one of these
# directories is compiled without editing this file.
SRC = $(SRC_DIR)/classes/*.cpp \
      $(SRC_DIR)/ModelData/*.cpp \
      $(SRC_DIR)/Model/*.cpp \
      $(SRC_DIR)/Equations/*.cpp

SRC_H = $(SRC_DIR)/classes/*.hpp \
        $(SRC_DIR)/ModelData/*.hpp \
        $(SRC_DIR)/Model/*.hpp \
        $(SRC_DIR)/Equations/*.hpp

# -----------------------------------------------------------------
# HYPRE — optional hypre BoomerAMG linear solver (experimental)
# -----------------------------------------------------------------
# 0 (default): hypre is not needed. src/Equations/sunlinsol_hypre.cpp
#    compiles to stubs, and running with SHUD_LINSOL=amg stops with an
#    error message.
# 1: compiles the wrapper that presents hypre BoomerAMG as a SUNDIALS
#    linear solver, and links hypre, MPI and OpenBLAS (hypre needs both
#    at link time). The solver is used only when the environment variable
#    SHUD_LINSOL=amg is set at run time; otherwise SPGMR is used and no
#    hypre function is called.
#
#   macOS (Homebrew): HYPRE_INCDIR=/opt/homebrew/include  HYPRE_LIBDIR=/opt/homebrew/lib
#   Ubuntu (apt)    : HYPRE_INCDIR=/usr/include/hypre     HYPRE_LIBDIR=/usr/lib/x86_64-linux-gnu
#
# Set the paths on the make command line or in the environment.
HYPRE ?= 0
HYPRE_INCDIR ?= /opt/homebrew/include
HYPRE_LIBDIR ?= /opt/homebrew/lib

# MPI_INCDIR — location of mpi.h. The hypre headers include it although
# SHUD calls no MPI function.
#   macOS (Homebrew): hypre is built without MPI; the variable is unused.
#   Ubuntu (apt)    : libopenmpi-dev, /usr/lib/<arch>-linux-gnu/openmpi/include
#                     (found automatically from `uname -m`)
MULTIARCH_DIR := /usr/lib/$(shell uname -m)-linux-gnu
MPI_INCDIR ?= $(MULTIARCH_DIR)/openmpi/include

# OPENBLAS_LIBDIR — location of the OpenBLAS library. Homebrew installs
# it outside the default linker path (/opt/homebrew/opt/openblas/lib),
# hence this default. On Ubuntu it is in the default path; set the
# variable to empty there.
OPENBLAS_LIBDIR ?= /opt/homebrew/opt/openblas/lib

INCLUDES = -I $(SUNDIALS_DIR)/include \
           -I $(INC_OMP) \
           -I $(SRC_DIR)/Model \
           -I $(SRC_DIR)/ModelData \
           -I $(SRC_DIR)/classes \
           -I $(SRC_DIR)/Equations \
           $(HYPRE_CK)

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
MPI_CXX_LIB ?= $(if $(wildcard $(MULTIARCH_DIR)/libmpi_cxx.so*),-lmpi_cxx)

# HYPRE_CK is part of INCLUDES and HYPRE_LK part of LK_FLAGS, so every
# recipe picks them up. Both are empty when HYPRE=0.
ifeq ($(HYPRE),1)
  HYPRE_CK = -DSHUD_USE_HYPRE=1 -I $(HYPRE_INCDIR) \
             $(if $(wildcard $(MPI_INCDIR)/mpi.h),-I $(MPI_INCDIR))
  HYPRE_LK = -L$(HYPRE_LIBDIR) -lHYPRE -lmpi $(MPI_CXX_LIB) \
             $(if $(OPENBLAS_LIBDIR),-L$(OPENBLAS_LIBDIR)) -lopenblas \
             -Wl,-rpath,$(HYPRE_LIBDIR)
else ifeq ($(HYPRE),0)
  HYPRE_CK =
  HYPRE_LK =
else
$(error HYPRE must be 0 or 1, got '$(HYPRE)')
endif

LK_FLAGS = -lm -lsundials_cvode -lsundials_nvecserial $(HYPRE_LK)
# LK_OMP holds only the OpenMP runtime library. libsundials_nvecopenmp is
# added through $(SHUD_NVEC_OMP_LK) when SHUD_USE_OPENMP_NVECTOR=1.
LK_OMP   = $(CXX_OPENMP_LFLAGS)
LK_DYLN  = "LD_LIBRARY_PATH=$(LIB_SUN)"

# -----------------------------------------------------------------
# SUNDIALS version and installation check
# -----------------------------------------------------------------
# - Major version 6 is required (tested with 6.0.0, 6.1.1, 6.4.1 and
#   6.7.0; ./configure installs 6.0.0). The pattern is anchored so that
#   60 or 600 does not pass as 6.
# - The libraries are checked too, to catch an installation where the
#   headers exist but the libraries were never built.
SUNDIALS_CFG_H = $(SUNDIALS_DIR)/include/sundials/sundials_config.h

.PHONY: check_sundials check_sundials_omp
check_sundials:
	@test -f $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: $(SUNDIALS_CFG_H) not found; run ./configure first"; exit 2)
	@grep -Eq '^#define SUNDIALS_VERSION_MAJOR 6$$' $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: SUNDIALS major != 6 in $(SUNDIALS_CFG_H); SHUD requires SUNDIALS 6.x"; exit 2)
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
	@echo "       make shud_omp     - build OpenMP shud_omp executable (default: parallel RHS + OpenMP N_Vector with serial sums; see OpenMP_Guide.md)"
	@echo "       make shud_omp SHUD_USE_OPENMP_NVECTOR=0 - serial N_Vector, parallel RHS only"
	@echo "       make shud_omp SHUD_NVEC_DETRED=1        - fastest: parallel sums in a fixed tree (output not bit-identical to the default build)"
	@echo "       make shud SHUD_DUMP_RHS=1     - serial build with RHS snapshot hooks compiled in"
	@echo "       make shud_omp SHUD_DUMP_RHS=1 - OpenMP build with RHS snapshot hooks compiled in"
	@echo "       make shud SHUD_ENABLE_PROFILE=1     - serial build with wall-clock profile timer compiled in"
	@echo "       make shud_omp SHUD_ENABLE_PROFILE=1 - OpenMP build with wall-clock profile timer compiled in"
	@echo "       make shud SHUD_ENABLE_OPENMP_RHS=1  - shud target with the parallel RHS compiled in"
	@echo "       make shud SHUD_USE_OPENMP_NVECTOR=1 - shud target with the OpenMP N_Vector backend"
	@echo "       make shud EXTRA_CXXFLAGS=-DSHUD_ENABLE_DIAGNOSTICS - serial build with extra CVODE keys (hlast/qlast) in cvode_stats.txt"
	@echo "       make shud_asan    - serial build with -fsanitize=address,undefined"
	@echo "       make smoke_configd                  - build + run the OpenMP N_Vector runtime probe"
	@echo "       make check_sundials - verify SUNDIALS 6.x install"
	@echo "       make clean        - remove binary outputs (preserves InstallSundials)"
	@echo
	@echo "       make shud HYPRE=1 - also build the experimental hypre BoomerAMG solver (run with SHUD_LINSOL=amg)."
	@echo "  It needs hypre, MPI and OpenBLAS:"
	@echo "  macOS brew:  HYPRE_INCDIR=/opt/homebrew/include HYPRE_LIBDIR=/opt/homebrew/lib (defaults)"
	@echo "  Ubuntu apt:  HYPRE_INCDIR=/usr/include/hypre    HYPRE_LIBDIR=/usr/lib/x86_64-linux-gnu"
	@echo "  Other:       HYPRE_INCDIR=<prefix>/include HYPRE_LIBDIR=<prefix>/lib"
	@echo

cvode CVODE:
	@echo '...Install SUNDIALS/CVODE 6.0.0 ...'
	chmod +x configure
	./configure
	@echo

# The option variables are part of the `shud` recipe too, so e.g.
# `make shud SHUD_ENABLE_OPENMP_RHS=1` works. With
# SHUD_USE_OPENMP_NVECTOR=1 the OpenMP N_Vector library must be
# installed, hence the choice of the check target.
ifeq ($(SHUD_USE_OPENMP_NVECTOR),1)
  SHUD_CHECK_TARGET := check_sundials_omp
else
  SHUD_CHECK_TARGET := check_sundials
endif
shud SHUD: $(SHUD_CHECK_TARGET) $(MAIN_shud) $(SRC) $(SRC_H)
	@echo '...Compiling shud (serial by default) ...'
	@echo  $(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_OMP_RHS_CK) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_DETRED_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(SHUD_OMP_RHS_LK)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_OMP_RHS_CK) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_DETRED_DEFINE) $(SHUD_NVEC_OMP_CK) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(SHUD_OMP_RHS_LK)
	@echo
	@echo " $(TARGET_EXEC) is compiled successfully!"
	@echo

# shud_omp — the OpenMP build. With the defaults chosen above for this
# goal (SHUD_ENABLE_OPENMP_RHS=1, SHUD_USE_OPENMP_NVECTOR=1,
# SHUD_NVEC_HYBRID=1) the RHS and the element-wise vector operations run
# in parallel, and the output is bit-identical at every thread count.
#   make shud_omp SHUD_NVEC_DETRED=1         parallel sums in a fixed tree (fastest)
#   make shud_omp SHUD_USE_OPENMP_NVECTOR=0  serial N_Vector, parallel RHS
#   make shud_omp SHUD_ENABLE_OPENMP_RHS=0   serial RHS
# Each option enters the recipe as a variable that is empty when the
# option is off.
shud_omp: check_sundials_omp $(MAIN_OMP) $(SRC) $(SRC_H)
	@echo '...Compiling shud_omp (OpenMP) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_DETRED_DEFINE) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(LK_OMP)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) $(SHUD_NVEC_OMP_DEFINE) $(SHUD_NVEC_HYBRID_DEFINE) $(SHUD_NVEC_DETRED_DEFINE) $(SHUD_DUMP_DEFINE) $(SHUD_OMP_RHS_DEFINE) $(SHUD_PROFILE_DEFINE) $(EXTRA_CXXFLAGS) $(INCLUDES) $(SHUD_PROFILE_INC) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(SHUD_PROFILE_SRC) $(LK_FLAGS) $(SHUD_NVEC_OMP_LK) $(LK_OMP)
	@echo
	@echo " $(TARGET_OMP) is compiled successfully!"
	@echo

# The SHUD sources without main.cpp, for programs that have their own
# main() (the unit test and libshud.a below). $(SRC) holds wildcards;
# they must be expanded with $(wildcard ...) before filter-out can match
# the file name.
SHUD_SRC_NOMAIN := $(filter-out $(SRC_DIR)/main.cpp,$(wildcard $(SRC)))

# -----------------------------------------------------------------
# smoke_configd — OpenMP N_Vector runtime probe
# -----------------------------------------------------------------
# Builds and runs tests/s1d_configd_nvec_smoke.cpp, a standalone program
# (not linked with the SHUD sources) that creates a vector with
# N_VNew_OpenMP, checks N_VGetVectorID == SUNDIALS_NVEC_OPENMP, and
# destroys it through the generic N_VDestroy.
#
# It links only libsundials_nvecopenmp and libsundials_generic (the
# latter for SUNContext_Create / SUNContext_Free), not $(LK_FLAGS).
.PHONY: smoke_configd
smoke_configd: check_sundials_omp tests/s1d_configd_nvec_smoke.cpp
	@echo '...Compiling s1d_configd_nvec_smoke (OpenMP N_Vector runtime probe) ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) -DSHUD_ENABLE_OPENMP_RHS=1 -DSHUD_USE_OPENMP_NVECTOR=1 $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/s1d_configd_nvec_smoke tests/s1d_configd_nvec_smoke.cpp -lsundials_nvecopenmp -lsundials_generic $(CXX_OPENMP_LFLAGS)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) -DSHUD_ENABLE_OPENMP_RHS=1 -DSHUD_USE_OPENMP_NVECTOR=1 $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/s1d_configd_nvec_smoke tests/s1d_configd_nvec_smoke.cpp -lsundials_nvecopenmp -lsundials_generic $(CXX_OPENMP_LFLAGS)
	@echo
	@echo '...Running s1d_configd_nvec_smoke (expect OK: N_VGetVectorID == SUNDIALS_NVEC_OPENMP) ...'
	@DYLD_LIBRARY_PATH=$(LIB_SUN) tests/s1d_configd_nvec_smoke || (echo 'FAIL: OpenMP N_Vector smoke test did not print OK'; exit 1)

# -----------------------------------------------------------------
# test_adjacency_fallback — unit test of the adjacency lists
# -----------------------------------------------------------------
# Builds and runs tests/test_adjacency_fallback.cpp. It sets up a
# Model_Data whose `.index` fields violate `index == array_index + 1`
# and checks that build_adjacency_lists() sets
# `adjacency_fallback_triggered` and builds the lists in array-index
# order.
#
# The test has its own main(), so it links $(SHUD_SRC_NOMAIN).
.PHONY: test_adjacency_fallback
test_adjacency_fallback: check_sundials tests/test_adjacency_fallback.cpp $(SRC) $(SRC_H)
	@echo '...Compiling test_adjacency_fallback ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/test_adjacency_fallback tests/test_adjacency_fallback.cpp $(SHUD_SRC_NOMAIN) $(LK_FLAGS)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) -I tests $(LIBRARIES) $(RPATH) -o tests/test_adjacency_fallback tests/test_adjacency_fallback.cpp $(SHUD_SRC_NOMAIN) $(LK_FLAGS)
	@echo
	@echo '...Running test_adjacency_fallback (expect PASS) ...'
	@DYLD_LIBRARY_PATH=$(LIB_SUN) tests/test_adjacency_fallback && echo 'OK: adjacency fallback unit test PASS' || (echo 'FAIL: adjacency fallback unit test'; exit 1)

# -----------------------------------------------------------------
# libshud.a — the SHUD sources without main.cpp, as a static library
# -----------------------------------------------------------------
# For programs that need the Model_Data interface (loadinput, initialize,
# rhs_core) without the CVODE main loop. Built with the fixed flags and
# none of the options above. Object files go to _libshud_obj/.
LIBSHUD_OBJ_DIR  := _libshud_obj
LIBSHUD_ARCHIVE  := libshud.a
LIBSHUD_SRC      := $(SHUD_SRC_NOMAIN)
LIBSHUD_OBJ      := $(patsubst $(SRC_DIR)/%.cpp,$(LIBSHUD_OBJ_DIR)/%.o,$(LIBSHUD_SRC))

$(LIBSHUD_OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: libshud.a
libshud.a: check_sundials $(LIBSHUD_OBJ)
	@echo '...Archiving libshud.a ...'
	$(AR) rcs $(LIBSHUD_ARCHIVE) $(LIBSHUD_OBJ)
	@echo " $(LIBSHUD_ARCHIVE) is archived successfully (fixed flags, no options)"
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

# -----------------------------------------------------------------
# Engineering gates (AGENTS.md, Enforcement Index)
# -----------------------------------------------------------------
# Thin entry points; the checks live in tools/ci/ and read their
# thresholds from constraints.yaml. `make ci` is what the CI runs.
#   make regress UPDATE=1   rewrite tests/reference/ (local only)
#   make test UPDATE=1      rewrite tests/io_contract/ (local only)
#   make pr-gates BASE=ref  checks against a PR base (default origin/master)
GATES := lint test regress coverage test-guardrails docs-check secrets pr-gates ci install-hooks
.PHONY: $(GATES)
lint test regress coverage test-guardrails docs-check secrets pr-gates ci install-hooks:
	@UPDATE=$(UPDATE) BASE=$(BASE) python3 tools/ci/gate.py $@
