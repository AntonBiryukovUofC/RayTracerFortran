# Ray Tracing in Fortran with Tomography

This repository contains the codes for running ray tracing in 1D layered cake model, implemented in Fortran. There are multiple alternatives with better usage / more user-friendly interface to my codes; though they are also "heavier" to leverage MCMC seismic tomography.

The main purpose of Ray Tracing here is to be super quick for the particular case of 1D, so I could easily switch different velocity models, and still be able to quickly trace the rays with arrival times.

# Usage

The repository contains both the scripts that invoke the "workhorse" functions, as well as the functions themselves.
First, you have to compile them so R can recognise the functions from the compiled shared objects:
- compile the function codes with `R CMD SHLIB subroutineR-quiet.f90`



I will populate later this section, after I "clean" this repository.