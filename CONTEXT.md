# CONTEXT.md

> Project identity, bounded contexts, and invariants for this repository.
> `AGENTS.md` is the operating contract; `openspec/glossary.md` defines the domain vocabulary; this file defines the boundaries and invariants that agents must respect when applying both.

## Project Identity

SHUD (Simulator for Hydrologic Unstructured Domains) is a multi-process, multi-scale hydrological model in which the major hydrological processes are fully coupled with the semi-discrete finite volume method on an unstructured triangular mesh.

- **Primary users / consumers**: hydrology researchers and watershed modellers; rSHUD prepares inputs and reads outputs.
- **Goal**: reproducible, physically consistent simulation of surface water, soil moisture, groundwater, rivers and lakes.
- **Lifecycle**: years.

## Domain Language

Canonical terms and their prohibited aliases live in `openspec/glossary.md`. Read it before naming a domain concept; never define a term here.

## Bounded Contexts

| Context | Owns | Key terms (defined in `openspec/glossary.md`) | Forbidden logic | Integration boundary |
|---------|------|-----------|-----------------|----------------------|
| Solver core (`src/Model/`) | time loop, CVODE calls, coupled right-hand side, vector sums | macro step, internal step, right-hand side, coupled mode | physical formulas; anything that depends on thread count | calls flux routines of the physics context; reads and writes only the state vector and flux arrays |
| Physical processes (`src/ModelData/MD_ET.cpp`, `MD_ElementFlux.cpp`, `MD_RiverFlux.cpp`, `src/classes/Element.cpp`, `River.cpp`, `Lake.cpp`, `src/Equations/`) | fluxes and storages of elements, rivers, lakes | state, flux, interception storage, element, river reach, lake element | solver settings; file reading and writing | called from the right-hand side or once per macro step |
| Input and output (`src/ModelData/MD_readin.cpp`, `MD_initialize.cpp`, `src/classes/IO.cpp`, `Model_Control.cpp`, `TabularData.cpp`, `TimeSeriesData.cpp`) | input file formats, output files, command line | state output, flux output, output interval | changing a computed value | formats are a public contract (`tests/io_contract/`) |
| Build and parallelism (`Makefile`, `configure`, `src/Model/MD_nvec_hybrid.cpp`) | compiler flags, build variants, SUNDIALS install | serial build, default OpenMP build, fastest build, bit-identical | options that trade reproducibility for speed without an explicit name | `OpenMP_Guide.md`, `OpenMP_NVector_Determinism.md` |

## Core Invariants

- The same input, commit, build variant and SUNDIALS build give the same output files, bit for bit, on the same platform.
- `shud` and the default `shud_omp` give bit-identical output at every thread count (with the SUNDIALS built by `./configure`). Checked by `make regress`.
- Inside the parallel region no two threads write the same memory, and sums across elements are accumulated in a fixed order.
- A change that alters model results says so: it updates `tests/reference/`, `VersionUpdate.md` and adds a record under `decisions/`.

Known violations of intended invariants, tracked as issues rather than stated here as facts: water is not conserved in the interception store (#15); the right-hand side is not a pure function of `(t, y)` (#16); flux outputs are not time integrals (#17). See epic #13.

## Public Interfaces and Contracts

| Interface | Contract source | Backward compatibility rule | Test seam |
|-----------|-----------------|-----------------------------|-----------|
| Command-line options | `src/classes/CommandIn.cpp` | never removed or re-purposed without approval | `tests/io_contract/cli_options.txt` |
| Output files (names, column counts, binary layout) | `src/classes/IO.cpp`, `src/classes/Model_Control.cpp` | unchanged without approval; rSHUD reads them | `tests/io_contract/ccw_outputs.txt` |
| Input file formats | `src/ModelData/MD_readin.cpp` | existing project folders must keep running | `make regress` runs `input/ccw` as shipped |
| Run-log lines read by scripts | `src/Model/shud.cpp` | wording kept | review-only |

## Forbidden Logic & Irreversible Operations

| Rule | Scope | Why |
|------|-------|-----|
| Do not change or bypass the reproducibility compiler flags, and do not introduce fast-math options | `Makefile`, `configure` | results would stop being reproducible |
| Do not update reference results or widen tolerances to make a gate pass | `tests/reference/`, `constraints.yaml` (regression) | the regression gate would confirm itself |
| Do not change input/output file formats or the run-log lines scripts depend on | input readers, output writers, `src/Model/shud.cpp` | existing projects and tools break |

## Open Terminology Questions

| Question | Why it matters | Candidate terms | Owner |
|----------|----------------|-----------------|-------|
| What is the user-facing name of the build variants? The run log prints "Config C / E / E2", the guide says serial-vector / default / fastest build | agents and users mix the two sets | keep the guide names, treat "Config X" as log labels only | maintainers |
| `LSM_STEP` / `ET_STEP` are read but unused (#24) | the term "land-surface step" has no behaviour behind it | remove the term, or implement it | maintainers |
