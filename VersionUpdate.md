# Simulator of Hydrologic Unstructured Domains (SHUD)

# Relation with PIHM family


## SHUD v1.0 (2019.12)

MODIFICATIONS/ADDITIONS in SHUD V1.0 from previous PIHM family.

  0. Change the language and structure of code from C to C++.
  1. Update the CVODE from v2.2 to v5.0.
  2. Support OpenMP Parrallel computing.
  3. Change the input/output format. Check the Manual of SHUD on github.
  4. Change the structure of River.
  5. The functions to handle the time-series data, including forcing, LAI,
     Roughness Length, Boundary Condition, Melting factor.
  6. Lake Module is added into the hydrological process.

## SHUD v2.0 (2022.04)

MODIFICATIONS/ADDITIONS from v1.0

1. Update to SUNDIALS 6.x
2. The units of forcing input.  
   1. Forcing data: Precipitation (mm/day), Temperature (C), Windspeed (m/s), Radiation (w/m2), Relative Humidity (0~1), Pressure (kPa).
   2. Landcover parameters: Rough (Manning's Roughness) from [day m^{1/3}] to [s m^{1/3}].
   3. River parameters: Rough (Manning's Roughness) from [day m^{1/3}] to [s m^{1/3}].
3. Add the Bucket Lake model. Water balance of a lake is: $ ds/dt = P + Q_surf + Q_sub + R_in - R_out - E $
4. Change of the names of inputfile .sp.rivseg(v2.0), instead of .sp.rivchn (v1.0)
5. The calculation of ET, particularly the Potential Evapotranspiration from Pennman-Monteith Equation.
6. More calibration parameters are open now. Total number is 38 or more.
7. Format of files:
    1. The number of columns of *.sp.att* file to 9 columns, that is "INDEX	SOIL	GEOL	LC	FORC	MF	BC	SS  iLake"
    2. Three table exist in the *.sp.riv* file: ggggRiver, parameter, points. Head of three tables are: **River**(Index	Down	Type	Slope	Length	BC), **Parameters**(Index	Depth	BankSlope	Width	Sinuosity	Manning	Cwr	KsatH	BedThick), **Points**(From.x	From.y	From.z	To.x	To.y	To.z)
    3. Change of the *.cfg.ic* file format, since the initial condition for lake stage is added. Three table (v2.0) (element, river reach and lake) exist within the file, instead of two tables (v1.0).
8. Temporary permafrost parameterization scheme is added; yet, the testing and validation is on the track.
9. Temperature decreases as elevation increases, dT/dz = 0.00065  Adiabatic Lapse Rate 6.5 [$K/km$]
99. Lots of bugs are fixed. 

 


## OpenMP CPU acceleration (2026.07, tag `cpu-accel-v1.1.1`)

MODIFICATIONS/ADDITIONS from v2.0. The physics and the input/output file formats are unchanged. How to build and run the parallel model is described in `OpenMP_Guide.md`.

1. Parallel computing.
   1. The right-hand side of the ODE system (fluxes of elements, river segments and lakes) is evaluated in parallel. Shared writes between threads were removed and all cross-element sums are accumulated in a fixed order, so the result does not depend on the number of threads.
   2. The vector operations inside CVODE are parallel. Element-wise operations use the SUNDIALS OpenMP vector; sums and norms are replaced by SHUD's own order-preserving versions (`src/Model/MD_nvec_hybrid.cpp`). See `OpenMP_NVector_Determinism.md`.
   3. Measured on a 40,046-element mesh with 16 threads: about 2.6 times faster than the serial model with the default build, about 3.6 times with `SHUD_NVEC_DETRED=1`.
   4. The former OpenMP code path (`src/ModelData/MD_f_omp.cpp`) is removed. `make shud_omp` now builds the new one.
   5. Internal changes that serve the above: the RHS is reorganized into `src/Model/MD_rhs_core.cpp`, the element data used by the RHS are stored as contiguous arrays (`src/ModelData/MD_layout.hpp`), and the element/river/lake neighbour lists are built once at start (`src/ModelData/MD_adjacency.cpp`).
2. Build.
   1. `./configure` installs SUNDIALS/CVODE 6.0.0 into `./InstallSundials` (before: `~/sundials`), and does nothing if it is already there. The Makefile checks for version 6.0.x and stops otherwise. Use `make SUNDIALS_DIR=...` for another location.
   2. The compiler flags are fixed to `-O2 -g -ffp-contract=off -fno-fast-math -std=c++14` (before: `-O3 -g -std=c++14`). `-ffast-math`, `-Ofast` and `-funsafe-math-optimizations` are rejected, because they break the reproducibility of the results.
   3. `HYPRE=1` builds an experimental hypre BoomerAMG linear solver and links hypre, MPI and OpenBLAS (paths: `HYPRE_INCDIR`, `HYPRE_LIBDIR`, `MPI_INCDIR`, `OPENBLAS_LIBDIR`). Off by default; the default build needs only SUNDIALS.
   4. Options of `make shud_omp`:
      - (none): parallel RHS, parallel element-wise vector operations, serial sums. Output is bit-identical to the serial-vector build.
      - `SHUD_NVEC_DETRED=1`: sums are parallel too, in a fixed tree. Fastest. Output is identical at every thread count, but not bit-identical to the other builds.
      - `SHUD_USE_OPENMP_NVECTOR=0`: serial vectors; only the RHS is parallel.
      - `SHUD_ENABLE_OPENMP_RHS=0`: serial RHS.
      - `SHUD_NVEC_HYBRID=0`: refused unless `SHUD_ALLOW_CONFIG_D=1`, because the output would depend on the number of threads.
   5. A plain `make` builds the serial model (`make all`).
   6. New targets: `make shud_asan` (AddressSanitizer + UndefinedBehaviorSanitizer build), `make libshud.a`, and the self-tests `make smoke_configd`, `make test_adjacency_fallback`.
3. Run.
   1. `NUM_OPENMP` in `.cfg.para` sets the number of threads of the CVODE vector layer.
   2. The environment variable `SHUD_RHS_THREADS` sets the number of threads of the RHS. If it is not set, `NUM_OPENMP` is used (`OMP_NUM_THREADS` in the serial-vector build).
   3. `OMP_PROC_BIND` should be set (e.g. `close`). If it is not, a warning is printed and the NUMA-aware memory initialization is skipped.
   4. `SHUD_SPGMR_MAXL` sets the maximum Krylov dimension of the linear solver. Accepted values: 5, 10, 15, 20, 30. If it is not set, the SUNDIALS default is used, as before. Setting it changes the solver's path, so the results differ slightly from the default.
4. Output. Two new text files in the output directory: `cvode_stats.txt` (CVODE solver statistics) and `nfcall.txt` (number of RHS evaluations).
5. Options for development and diagnosis. All are off by default, and the model output is the same as without them unless noted.
   1. Build options: `SHUD_DUMP_RHS=1` (snapshots of the RHS; controlled at run time by `SHUD_DUMP_OUTPUT_DIR`, `SHUD_DUMP_CASE_ID`, `SHUD_DUMP_SITE`, `SHUD_DUMP_FNAME_SUFFIX`, `SHUD_DUMP_T_VALUES`, `SHUD_DUMP_T_TOL`), `EXTRA_CXXFLAGS=-DSHUD_ENABLE_DIAGNOSTICS` (more keys in `cvode_stats.txt`), `SHUD_ENABLE_PROFILE=1` (wall-clock timers, written to `profile_B0.yaml` in the output directory; sources in `tools/profile/`).
   2. Environment variables: `SHUD_NVEC_PROF`, `SHUD_DUMP_CV_Y`, `SHUD_TELEMETRY_TSV`, `SHUD_ADJACENCY_LOG`.
   3. Experimental solver settings, which **change the results**: `SHUD_LINSOL=amg` (hypre BoomerAMG linear solver in place of SPGMR; needs a build with `HYPRE=1`; tested and not adopted), `SHUD_AMG_TOL`, `SHUD_CVODE_RELTOL`, `SHUD_CVODE_EPSLIN`.
6. Bugs fixed.
   1. Division by zero in the accumulated-temperature average when the queue is empty (`src/classes/AccTemperature.hpp`).
   2. Memory leaks at the end of a run in `Model_Data`, the lake and the time-series classes; uninitialized pointers in `FloodAlert`, `Model_Data` and the lake classes.

The full development record (benchmarks, validation tools, decisions) is in https://github.com/DankerMu/SHUD-OpenMP.
