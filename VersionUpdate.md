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

 


