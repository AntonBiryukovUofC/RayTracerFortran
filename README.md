![Image of Raytracing](https://github.com/AntonBiryukovUofC/RayTracerFortran/blob/master/pngs/ray_example.png)

# Ray Tracing in Fortran with Tomography

This repository contains the codes for running ray tracing in 1D layered cake model, implemented in Fortran. There are multiple alternatives with better usage / more user-friendly interface to my codes; though they are also "heavier" to leverage MCMC seismic tomography.

The main purpose of Ray Tracing here is to be super quick for the particular case of 1D, so I could easily switch different velocity models, and still be able to quickly trace the rays with arrival times.

# Usage

The repository contains both the scripts that invoke the "workhorse" functions, as well as the functions themselves.


1. Compile the function codes with `R CMD SHLIB subroutineR-quiet.f90`
2. In R, load the shared object: `dyn.load("./subroutineR-quiet.so")` 
3. Call the function in this fashion: `timeP <- .Fortran("dff" ,vels=as.numeric(v),
										depths=as.numeric(d),NLayers=as.integer(NLayers),
										src_offset=as.numeric(srco),
										src_depth=as.numeric(srcd),
										NSrc=as.integer(NSrc),
										timeP = as.numeric(timeP))`


I will populate later this section, after I "clean" this repository. For now, look at the details in **rayTracerR.R**.

# TODO:

## Parallel tempering Fortran codes:

1. READPARFILE() needs adapting for Ray Tracer (RT)
2. READDATA() needs adapting for Ray Tracer travel times 
3. Figure out what INTERPLAYER routine does 
4. Fix PROPOSAL_SDRT() for ray tracing