setwd('/home/geoanton/RayTracerFortran/RayTracer/')
library(dplyr)

read.rays <- function(file = "rays.dat",nrays=1)
{
  ind.delta <- 2*seq(nrays)-1
  df <- tibble()
  for (i in seq(nrays))
  {
    delta <- scan(file,skip = ind.delta[i]-1,nlines = 1)
    depth <- scan(file,skip = ind.delta[i],nlines = 1)  
    ray <- i
    df <- bind_rows(df,tibble(delta=delta,depth=depth,ray=ray))
  }
  df <- df %>% group_by(ray) %>% mutate(dp.cum = cumsum(depth), delta.cum = cumsum(delta) )
  return(df)
}



dyn.load("./subroutineR-quiet.so")
src.df <- read.csv('/home/geoanton/RayTracerFortran/RayTracer/eqdf.csv',header=T)
srcd <- src.df$z
srco <- sqrt(src.df$x**2+src.df$y**2)
v <- c(3100,3270,5000)
d <- c(2000, 4000)
NSrc <- length(srcd)
NLayers <- length(d)
timeP <- rep(0,NSrc)
kd <- 10 # keep rays in rays.dat ???
timeP <- .Fortran("dff" ,vels=as.numeric(v),depths=as.numeric(d),NLayers=as.integer(NLayers),
                  src_offset=as.numeric(srco),src_depth=as.numeric(srcd),NSrc=as.integer(NSrc),
                  timeP = as.numeric(timeP),keep_delta=as.integer(kd))


rays.df <- read.rays(nrays=NSrc)