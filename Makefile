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
CXX_OPENMP_DEFINE = -D_OPENMP_ON

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
# Layer 2 (findstring on MAKEOVERRIDES): catches CLI assignments to ANY
# project-local lock variable (e.g. `make SHUD_BUILD_CFLAGS=-Ofast`),
# which the `override :=` directive would otherwise silently ignore.
# Without this layer the user could think "the build accepted my flag"
# while in reality the locked flag set is still used — Layer 2 turns the
# silent rejection into a loud error.
DISALLOWED_FLAGS := -ffast-math -Ofast -funsafe-math-optimizations
ifneq (,$(filter $(DISALLOWED_FLAGS),$(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(MAKEOVERRIDES) $(MAKEFLAGS) $(SHUD_BUILD_CFLAGS) $(CXX_BASE_FLAGS)))
$(error disallowed flag detected (one of $(DISALLOWED_FLAGS)) in CFLAGS/CXXFLAGS/CPPFLAGS/LDFLAGS/MAKEOVERRIDES/MAKEFLAGS/SHUD_BUILD_CFLAGS/CXX_BASE_FLAGS; B0 baseline requires strict IEEE-754 — see docs/build_manifest.md §1)
endif
ifneq (,$(or $(findstring -ffast-math,$(MAKEOVERRIDES)),$(findstring -Ofast,$(MAKEOVERRIDES)),$(findstring -funsafe-math-optimizations,$(MAKEOVERRIDES))))
$(error disallowed flag detected (one of $(DISALLOWED_FLAGS)) in a make-CLI assignment (MAKEOVERRIDES=[$(MAKEOVERRIDES)]); attempts to inject via SHUD_BUILD_CFLAGS / CXX_BASE_FLAGS / etc are also rejected; B0 baseline requires strict IEEE-754 — see docs/build_manifest.md §1)
endif

# -----------------------------------------------------------------
# Platform-conditional OpenMP flags
# -----------------------------------------------------------------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LIBOMP_PREFIX     ?= $(shell brew --prefix libomp 2>/dev/null)
  # Guard libomp only when building shud_omp; serial `make shud` is fine
  # without libomp installed.
  ifeq ($(filter shud_omp,$(MAKECMDGOALS)),shud_omp)
    ifeq ($(LIBOMP_PREFIX),)
$(error libomp not found via 'brew --prefix libomp'; run 'brew install libomp' before make shud_omp)
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

SRC = $(SRC_DIR)/classes/*.cpp \
      $(SRC_DIR)/ModelData/*.cpp \
      $(SRC_DIR)/Model/*.cpp \
      $(SRC_DIR)/Equations/*.cpp

SRC_H = $(SRC_DIR)/classes/*.hpp \
        $(SRC_DIR)/ModelData/*.hpp \
        $(SRC_DIR)/Model/*.hpp \
        $(SRC_DIR)/Equations/*.hpp

INCLUDES = -I $(SUNDIALS_DIR)/include \
           -I $(INC_OMP) \
           -I $(SRC_DIR)/Model \
           -I $(SRC_DIR)/ModelData \
           -I $(SRC_DIR)/classes \
           -I $(SRC_DIR)/Equations

# Use $(if …) so an empty $(LIB_OMP) / $(LIB_SYS) does NOT emit a bare `-L`
# token (which gobbles the next argument and breaks the link line on Linux,
# where LIB_OMP is empty by default).
LIBRARIES = $(if $(LIB_OMP),-L$(LIB_OMP)) \
            -L$(LIB_SUN) \
            $(if $(LIB_SYS),-L$(LIB_SYS))

RPATH = '-Wl,-rpath,$(LIB_SUN)'

LK_FLAGS = -lm -lsundials_cvode -lsundials_nvecserial
# OpenMP build additionally needs sundials_nvecopenmp + platform LFLAGS.
LK_OMP   = $(CXX_OPENMP_LFLAGS) -lsundials_nvecopenmp
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
	@echo "       make check_sundials - verify SUNDIALS 6.x install"
	@echo "       make clean        - remove binary outputs (preserves InstallSundials)"
	@echo

cvode CVODE:
	@echo '...Install SUNDIALS/CVODE 6.0.0 ...'
	chmod +x configure
	./configure
	@echo

shud SHUD: check_sundials $(MAIN_shud) $(SRC) $(SRC_H)
	@echo '...Compiling shud (B0 serial) ...'
	@echo  $(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(LK_FLAGS)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(LK_FLAGS)
	@echo
	@echo " $(TARGET_EXEC) is compiled successfully!"
	@echo

shud_omp: check_sundials_omp $(MAIN_OMP) $(SRC) $(SRC_H)
	@echo '...Compiling shud_OpenMP ...'
	@echo $(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) $(CXX_OPENMP_DEFINE) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(LK_FLAGS) $(LK_OMP)
	@echo
	$(CXX) $(SHUD_BUILD_CFLAGS) $(CXX_OPENMP_CFLAGS) $(CXX_OPENMP_DEFINE) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo " $(TARGET_OMP) is compiled successfully!"
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
	@echo "Done. (InstallSundials/ preserved)"
	@echo
