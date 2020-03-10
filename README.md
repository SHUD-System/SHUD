# Simulator of Hydrological Unstructured Domain

## Brief

The Simulator of Hydrological Unstructured Domain  (SHUD - pronounced “SHOULD”) is a multi-process, multi-scale hydrological model where major hydrological processes are fully coupled using the semi-discrete **Finite Volume Method** (FVM).

- **Maintainner**: Lele Shu ([lele.shu@gmail.com](mailto:lele.shu@gmail.com))
- **Website**: [shud.xyz](shud.xyz)
- **User Guide**: https://www.shud.xyz/_book/
- **Support tools**: SHUD-tools in R (<https://github.com/SHUD-System/SHUDtoolbox)
- **Programming**: C/C++
- **Platform**: Mac OS, Linux and Windows
- **Required library**:  SUNDIALS/CVODE V5.0+
- **Parallelization** : OpenMP

## Overview

The Simulator of Hydrological Unstructured Domain  (SHUD - pronounced “SHOULD”) is a multi-process, multi-scale hydrological model where major hydrological processes are fully coupled using the semi-discrete **Finite Volume Method** (FVM).

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

- 

## How to compile (Linux and Mac)

**Step 1: Install SUNDIALS/CVODE:**

```
./configure
```

This configure is to download teh SUNDIALS from GitHub and install it on your 

**Step 2: Compile SHUD with gcc**

```
make clean
make shud

```

**Step 3: Run the Cache Creek Watershed example**

```
./shud ccw
```

**Step4: Analysis the results of modeling.**

The output files from the SHUD model is save in `./output/ccw.out`.  The R package, SHUDtoolbox, helps to load the input/output files of SHUD. More details about prepare SHUD data, model input/output and visualization is available in SHUD website (https://www.shud.xyz) and help information of SHUDtoolbox.



