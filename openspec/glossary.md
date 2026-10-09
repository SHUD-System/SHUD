# Glossary

Terms of the SHUD hydrological model: how it steps in time, what it stores, how it is built and what its mesh is made of.

## Solving and time steps

**Macro step**:
The interval after which forcing, interception and snow are updated and the solver is called again; set by `MAX_SOLVER_STEP` (minutes).
_Avoid_: solver step, time step, coupling step

**Internal step**:
A step CVODE chooses by itself inside a macro step to meet the error tolerances.
_Avoid_: time step, sub-step

**Output interval**:
The period over which a variable is averaged before one record is written; one `DT_*` setting per output file.
_Avoid_: print step

**Right-hand side**:
The function that gives the time derivative of every state for a given time and state vector; the code calls it `f` and `rhs_core`.
_Avoid_: RHS function when a specific phase is meant, flux routine

**Coupled mode**:
The default run mode, in which all states are solved together in one CVODE system.
_Avoid_: implicit mode, full mode

**Uncoupled mode**:
The run mode selected with `-g`, in which surface, unsaturated zone, groundwater and river are solved by separate CVODE systems.
_Avoid_: explicit mode, decoupled mode

## States and fluxes

**State**:
A storage that CVODE integrates: surface water depth, unsaturated-zone storage and groundwater depth of each element, stage of each river reach, stage of each lake.
_Avoid_: variable, prognostic variable

**Interception storage**:
Water held on the canopy of an element; updated once per macro step and not part of the CVODE state vector.
_Avoid_: canopy state

**Snow storage**:
Snow water equivalent on an element; updated once per macro step and not part of the CVODE state vector.
_Avoid_: snow state

**Flux**:
A rate of water exchange between two storages or across the domain boundary.
_Avoid_: flow when a specific exchange is meant

**State output**:
An output file of a state; its records are the solution at the output time.
_Avoid_: y output

**Flux output**:
An output file of a flux; its records are averages of values taken at the end of each macro step, not time integrals.
_Avoid_: q output, cumulative flux

## Build variants

**Serial build**:
The model built by `make shud`: no OpenMP, serial vector operations.
_Avoid_: baseline build

**Default OpenMP build**:
The model built by `make shud_omp`: parallel right-hand side and element-wise vector operations, serial sums; the run log labels it "Config E".
_Avoid_: hybrid build, parallel build

**Fastest build**:
`make shud_omp SHUD_NVEC_DETRED=1`: sums are parallel too, in a fixed tree; the run log labels it "Config E2".
_Avoid_: DETRED build

**Serial-vector build**:
`make shud_omp SHUD_USE_OPENMP_NVECTOR=0`: parallel right-hand side, serial vector operations; the run log labels it "Config C".
_Avoid_: reference build

**Bit-identical**:
Two runs whose output files are equal byte for byte.
_Avoid_: identical, same results, equivalent

## Mesh and units

**Element**:
One triangle of the mesh; the unit that carries the surface, unsaturated-zone and groundwater states.
_Avoid_: cell, grid, triangle

**River reach**:
One entry of the river table; the unit that carries a river stage.
_Avoid_: river cell, channel

**River segment**:
The part of a river reach that lies in one element; the unit of river–element exchange.
_Avoid_: reach, sub-reach

**Lake element**:
An element flagged as part of a lake; it has no land-surface states of its own.
_Avoid_: lake cell

**Bank element**:
A land element that shares an edge with a lake element; the unit of lake–land exchange.
_Avoid_: shore cell

**Centroid**:
The centre of an element used for distances between elements; the mean of its three vertices.
_Avoid_: centre, circumcentre
