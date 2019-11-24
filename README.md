# Solver for Hydrological Unstructured Domain

## Brief

The Solver for Hydrologic Unstructured Domain  (SHUD - pronounced “SHOULD”) is a multi-process, multi-scale hydrological model where major hydrological processes are fully coupled using the semi-discrete **Finite Volume Method** (FVM). 

* **Maintainner**: Lele Shu (lele.shu@gmail.com)
* **Website (ongoing)**: https://SHUD-system.github.io
* **Document (ongoing)**: https://github.com/SHUD-System/SHUD_Doc
* **Programming**: C++
* **Platform**: Mac OS, Linux and Windows
* **Required library**:  SUNDIALS/CVODE V5.0+
* **Parallelization** : OpenMP
* **Support tools**: SHUD-tools in R (https://github.com/SHUD-System/SHUD-tool)

## Overview

The Solver for Hydrologic Unstructured Domain  (SHUD - pronounced “SHOULD”) is a multi-process, multi-scale hydrological model where major hydrological processes are fully coupled using the semi-discrete **Finite Volume Method** (FVM). 

SHUD encapsulates the strategy for the synthesis of multi-state distributed hydrological models using the integral representation of the underlying physical process equations and state variables. As a heritage of **Penn State Integrated Hydrologic Model (PIHM)**, the SHUD model is a continuation of 16 years of PIHM modeling in hydrology and related fields since the release of its first PIHM version (Qu, 2004). 

The SHUD’s design is based on a concise representation of a watershed and river basin’s hydrodynamics, which allows for interactions among major physical processes operating simultaneously, but with the flexibility to add or drop states/processes/constitutive relations depending on the objectives of the numerical experiment for research purpose. 

The SHUD is a distributed hydrological model in which the domain is discretized using an unstructured triangular irregular network (e.g., *Delaunay triangles*) generated with constraints (geometric and parametric). A local prismatic control volume is formed by the vertical projection of the Delaunay triangles forming each layer of the model. Given a set of constraints (river network, watershed boundary, elevation, and hydraulic properties), an *“optimized mesh”* is generated. The *“optimized mesh”* indicates the hydrological processes with the unstructured mesh can be calculated efficiently, stably and rationally (Farthing and Ogden, 2017; Vanderstraeten and Keunings, 1995; Kumar, Bhatt and Duffy, 2009). River volume elements are also prismatic, with trapezoidal or rectangular cross-section, and maintain the topological relation with the Delaunay triangles. The local control volumes encapsulate all equations to be solved and are herein referred to as the model kernel. 



### We now summarize the formulation and results from SHUD. 

- SHUD is a physically-based model, in which all equations used to emerge from the physics behind the hydrological processes within a catchment. The physical model can predict the water in an ungaged water system.  SHUD represents the spatial heterogeneity that influences the hydrology of the region. Consequently, it is practical to couple the SHUD model with models from biochemistry, reaction transport, geomorphology, limnology and other related research areas.
- SHUD is a fully-coupled hydrological model, where the conservative hydrological fluxes are calculated within the same time step. The state variables are the height of ponding water on the land surface, soil moisture, groundwater level, and river stage, while fluxes are infiltration, overland flow, groundwater recharge, lateral groundwater flow, river discharge, and exchange between river and elements. 
- The global ODE system solved in SHUD integrates all local ODE systems over the domain and solves with a state-of-the-art parallel ODE solver known as CVODE (Hindmarsh et al., 2005) developed at the Lawrence Livermore National Laboratory. 
- SHUD permits adaptable temporal and spatial resolution. The spatial resolution of the model varies from cen- timeters to kilometers based on modeling requirements computing resources. The internal time step of the iteration is adjustable; it can export the status of a catchment at time-intervals from a minutes to days. The flexible spatial and temporal resolution of the model is valuable for community model coupling. 
- SHUD can estimate either a long-term hydrological yield or a single-event flood. 
- SHUD is an open-source model --- anyone can access the source code and submit their modifications/improvements.

### Differences from PIHM

As a descendent of PIHM, SHUD inherits the fundamental idea of solving hydrological variables in CVODE. The code has been completely rewritten in a new programming language, with a new discretization and corresponding improvements to the underlying algorithms, adapting new mathematical schemes and a new user-friendly input/output data format. Although SHUD is forked from PIHM’s track, SHUD still inherits the use of CVODE for solving the ODE system, but modernizes and extends PIHM’s technical and scientific capabilities. The major differences are following: 

1. SHUD is written in C++, an object-oriented programming language with functionality to avoid risky memory leaks from C. Every functions in the code has been rewritten, so the functions, algorithm or data structure between SHUD and PIHM are incompatible. 
2. SHUD implements a re-design of the calculation of water exchange between hill slope and river. The PIHM defines the river channel as adjacent to bank elements – namely, the river channel shares the edges with bank elements. This design leads to sink problems in elements that share one node with a starting river channel. 
3. The mathematical equations used in infiltration, recharge, overland flow and river discharge are different among the two models. This change is so essential that the model results would be different with the same parameter set. 
4. SHUD adds mass-balance control within the calculation of each layer of elements and river channels, critical for long-term or micro-scale hydrologic modeling. 

We now briefly summarize the technical model improvements and technical capabilities of the model, compared to PIHM. This elaboration of the relevant technical features aims to assist future developers and advanced users with model coupling. Compared with PIHM, SHUD ... 

- supports the latest implicit Sundial/CVODE solver up to version 5.0.0 (the most recent version at time of writing), 
- supports OpenMP parallel computation, 
- redesigns the data structures and algorithm with object-oriented programming (C++), 
- supports human-readable input/output files and filenames, 
- exposes unified functions to handle the time-series data, including forcing, leaf area index, roughness length, boundary conditions and melt factor, 
- exports model initial condition at specific intervals that can be used for warm starts of continued simulation, 
- automatically checks the range of physical parameters and forcing data, 
- adds a debug mode that monitors potential errors in parameters and memory operations. 



