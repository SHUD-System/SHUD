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
# B0 baseline lock — disallowed-flag guard
# -----------------------------------------------------------------
# Fast-fail if the user attempts to inject IEEE-754-violating flags.
DISALLOWED_FLAGS := -ffast-math -Ofast -funsafe-math-optimizations
ifneq (,$(filter $(DISALLOWED_FLAGS),$(CXXFLAGS) $(MAKEOVERRIDES) $(MAKEFLAGS)))
$(error disallowed flag detected (one of $(DISALLOWED_FLAGS)); B0 baseline requires strict IEEE-754 — see docs/build_manifest.md)
endif

# -----------------------------------------------------------------
# B0 baseline lock — flag variables (do not edit without OpenSpec change)
# -----------------------------------------------------------------
CXX_BASE_FLAGS    = -O2 -g -ffp-contract=off -fno-fast-math -std=c++14
CXX_OPENMP_DEFINE = -D_OPENMP_ON

# Existing rules read $(CFLAGS); keep them in sync with the locked flags.
CFLAGS = $(CXX_BASE_FLAGS)

# -----------------------------------------------------------------
# Platform-conditional OpenMP flags
# -----------------------------------------------------------------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LIBOMP_PREFIX     ?= $(shell brew --prefix libomp 2>/dev/null)
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
# -----------------------------------------------------------------
CXX   ?= g++
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

LIBRARIES = -L $(LIB_OMP) \
            -L $(LIB_SUN) \
            -L $(LIB_SYS)

RPATH = '-Wl,-rpath,$(LIB_SUN)'

LK_FLAGS = -lm -lsundials_cvode -lsundials_nvecserial
# OpenMP build additionally needs sundials_nvecopenmp + platform LFLAGS.
LK_OMP   = $(CXX_OPENMP_LFLAGS) -lsundials_nvecopenmp
LK_DYLN  = "LD_LIBRARY_PATH=$(LIB_SUN)"

# -----------------------------------------------------------------
# SUNDIALS major-version guard
# -----------------------------------------------------------------
SUNDIALS_CFG_H = $(SUNDIALS_DIR)/include/sundials/sundials_config.h

.PHONY: check_sundials
check_sundials:
	@test -f $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: $(SUNDIALS_CFG_H) not found; run ./configure first"; exit 2)
	@grep -q "define SUNDIALS_VERSION_MAJOR 6" $(SUNDIALS_CFG_H) || \
	  (echo "ERROR: wrong SUNDIALS major version (expected 6) in $(SUNDIALS_CFG_H)"; exit 2)

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
	@echo  $(CXX) $(CFLAGS) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(LK_FLAGS)
	@echo
	$(CXX) $(CFLAGS) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_EXEC) $(MAIN_shud) $(SRC) $(LK_FLAGS)
	@echo
	@echo " $(TARGET_EXEC) is compiled successfully!"
	@echo

shud_omp: check_sundials $(MAIN_OMP) $(SRC) $(SRC_H)
	@echo '...Compiling shud_OpenMP ...'
	@echo $(CXX) $(CFLAGS) $(CXX_OPENMP_CFLAGS) $(CXX_OPENMP_DEFINE) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(LK_FLAGS) $(LK_OMP)
	@echo
	$(CXX) $(CFLAGS) $(CXX_OPENMP_CFLAGS) $(CXX_OPENMP_DEFINE) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $(TARGET_OMP) $(MAIN_OMP) $(SRC) $(LK_FLAGS) $(LK_OMP)
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
