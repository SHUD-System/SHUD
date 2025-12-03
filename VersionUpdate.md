# Simulator of Hydrologic Unstructured Domains (SHUD)


## SHUD v2.1 (2025.12)

MODIFICATIONS/ADDITIONS from v2.0

1. **Real-Time System (TimeManager) Implementation**
   - Introduced unified time management system with `TimeManager` class
   - Established `modelBaseDate` as the single authoritative time reference from forcing file (YYYYMMDD format)
   - Implemented comprehensive time variable tracking (t0, t1, dt, jd, timelocal, timeutc)
   - Added high-precision time handling using C++11 `std::chrono` library
   - Unified time formatting with 7 output formats (from year to microsecond precision)

2. **Time Stamp Validation System**
   - Automatic validation of timestamps across all input files (forcing, LAI, TSD files)
   - Detection and reporting of timestamp inconsistencies with detailed warnings
   - Interactive user confirmation mechanism for handling inconsistent timestamps
   - Clear warning messages showing file names, timestamps, and differences from base date

3. **Enhanced Output and User Experience**
   - All output files now include `modelBaseDate` timestamp in headers
   - Screen output enhanced with real-time information alongside relative time
   - Model summary displays base time, start time, and end time in human-readable format
   - Improved progress tracking with actual calendar dates during simulation

4. **Code Quality and Testing**
   - Added `TimeManager` class with comprehensive time management functionality
   - Extended `_TimeSeriesData` class with `getStartTime()` accessor method
   - Integrated time validation into model initialization workflow
   - Created comprehensive test suites for regression and timestamp validation
   - All existing test cases (ccw, heihe, qhh) pass with 100% success rate

5. **Technical Improvements**
   - Minimal code changes following "minimum modification, maximum stability" principle
   - Backward compatible with existing input/output formats
   - No performance degradation (time system overhead is negligible)
   - Clean separation of concerns with dedicated TimeManager class
   - Proper encapsulation and const-correctness in time-related APIs

6. **Documentation**
   - Complete specification documents (requirements, design, tasks)
   - Comprehensive verification reports for all implemented features
   - Detailed test scripts for integration and regression testing
   - User-friendly error messages and validation feedback

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



# Relation with PIHM family

SHUD is developed from PIHM family. The main differences are:

1. SHUD is developed in C++, while PIHM is developed in Fortran.
2. SHUD is developed in a modularized structure, while PIHM is developed in a monolithic structure.
3. SHUD is developed in a object-oriented structure, while PIHM is developed in a procedural structure.
4. SHUD is developed in a parallelized structure, while PIHM is developed in a sequential structure.
5. SHUD is developed in a platform-independent structure, while PIHM is developed in a platform-dependent structure.
6. SHUD is developed in a user-friendly structure, while PIHM is developed in a user-unfriendly structure.
7. SHUD is developed in a flexible structure, while PIHM is developed in a rigid structure.
8. SHUD is developed in a scalable structure, while PIHM is developed in a non-scalable structure.

