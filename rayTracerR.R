setwd('/home/geoanton/RayTracerFortran/RayTracer/')

dyn.load("./subroutineR-quiet.so")
src.df <- read.csv('/home/geoanton/RayTracerFortran/RayTracer/eqdf.csv',header=T)
srcd <- src.df$z
srco <- sqrt(src.df$x**2+src.df$y**2)
v = c(3100,3270,5000)
d= c(2000, 4000)
NSrc <- length(srcd)
NLayers <- length(d)
timeP <- rep(0,NSrc)
timeP <- .Fortran("dff" ,vels=as.numeric(v),depths=as.numeric(d),NLayers=as.integer(NLayers),
                  src_offset=as.numeric(srco),src_depth=as.numeric(srcd),NSrc=as.integer(NSrc),timeP = as.numeric(timeP))
