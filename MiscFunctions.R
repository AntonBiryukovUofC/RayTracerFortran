# Miscellaneous functions here for later use:

#' Reads rays from the subroutine output file
#'
#' @param file 
#' @param nrays 
#' @param timeP 
#'
#' @return
#' @export
#'
#' @examples
read.rays <- function(file = "rays.dat",nrays=1,timeP= NA)
{
  ind.delta <- 2*seq(nrays)-1
  df <- tibble()
  for (i in seq(nrays))
  {
    delta <- c(0,scan(file,skip = ind.delta[i]-1,nlines = 1))
    depth <- c(0,scan(file,skip = ind.delta[i],nlines = 1))
    ray <- i
    df <- bind_rows(df,tibble(delta=delta,depth=depth,ray=ray,ray.label = round(timeP[i],digits=4)  ))
  }
  df <- df %>% group_by(ray) %>% mutate(depth.cum = cumsum(depth), delta.cum = cumsum(delta) )
  return(df)
}

#' Draws profiles from the chain 
#'
#' @param x 
#' @param NLayers 
#'
#' @return
#' @export
#'
#' @examples
drawProfiles <- function(x,NLayers=NA)
{
  d=as.numeric(x[(NLayers+2) : (2*(NLayers+1))])
  v=as.numeric(x[1:(NLayers+1)])
  f <- approxfun(d,v,method = "constant")
  dnew <- seq(1,8000,length.out = 200)
  vnew <- f(dnew)
  vnew[is.na(vnew)] = v[length(v)]
  dd <- d[-1]
  for (i in seq(length(dd)) )
  {
    ptsX=seq(from=min(v[i],v[i+1]),to = max(v[i],v[i+1]),by=40)
    ptsY=rep(dd[i],length(ptsX))
    dnew=c(dnew,ptsY)
    vnew=c(vnew,ptsX)
  }
  prof <-tbl_df(data.frame(d=dnew,v=vnew))
  return(prof)
}


#' Creates a valid parameter.dat file for MCMC given settings in param.list
#'
#' @param param.list 
#' @param example.file 
#' @param output.file 
#'
#' @return
#' @export
#'
#' @examples
makeParameterFile <- function(param.list=NA,example.file = '~/RayTracerFortran/RayTracer/example_parameter.dat',output.file ="~/RayTracerFortran/RayTracer/test_parameter.dat")
{
  if (is.na(param.list))
  {
    message('[INFO]: Param list is empty!...No files were changed.')
    return(NA)
  }
  
  fileContent <- readChar(example.file, file.info(example.file)$size)
  fileContent <- gsub("NPL_VAL",replacement = param.list$NPL,fileContent)
  fileContent <- gsub("ISMPPRIOR_VAL",replacement = param.list$ISMPPRIOR,fileContent)
  fileContent <- gsub("IEXCHANGE_VAL",replacement = param.list$IEXCHANGE,fileContent)
  fileContent <- gsub("NDAT_RT_VAL",replacement = param.list$NDAT_RT,fileContent)
  fileContent <- gsub("NMODE_VAL",replacement = param.list$NMODE,fileContent)
  fileContent <- gsub("NSRC_VAL",replacement = param.list$NSRC,fileContent)
  fileContent <- gsub("NLMN_VAL",replacement = param.list$NLMN,fileContent)
  fileContent <- gsub("NLMX_VAL",replacement = param.list$NLMX,fileContent)
  fileContent <- gsub("ICHAINTHIN_VAL",replacement = param.list$ICHAINTHIN,fileContent)
  fileContent <- gsub("NKEEP_VAL",replacement = param.list$NKEEP,fileContent)
  fileContent <- gsub("NPTCHAINS1_VAL",replacement = param.list$NPTCHAINS,fileContent)
  fileContent <- gsub("dTlog_VAL",replacement = param.list$dTlog,fileContent)
  fileContent <- gsub("lambda_VAL",replacement = param.list$lambda,fileContent)
  fileContent <- gsub("hmx_VAL",replacement = param.list$hmx,fileContent)
  fileContent <- gsub("hmn_VAL",replacement = param.list$hmn,fileContent)
  fileContent <- gsub("dVs_VAL",replacement = param.list$dVs,fileContent)
  fileContent <- gsub("dVpVs_VAL",replacement = param.list$dVpVs,fileContent)
  fileContent <- gsub("VpVsmin_VAL",replacement = param.list$VpVsmin,fileContent)
  fileContent <- gsub("VpVsmax_VAL",replacement = param.list$VpVsmax,fileContent)
  fileContent <- gsub("sdmn_VAL",replacement = paste0(param.list$sdmin,collapse = ' '),fileContent)
  fileContent <- gsub("sdmx_VAL",replacement = paste0(param.list$sdmax,collapse = ' '),fileContent)
  
  
  con <- file(output.file, "wb")
  writeChar(fileContent, con, nchar(fileContent), eos = NULL)
  close(con)
}


#' Creates a Maximum a posteriori file _MAP.dat given the param.list with MAP.* values of the most likely geological scenario
#'
#' @param param.list 
#' @param output.file 
#'
#' @return
#' @export
#'
#' @examples
makeMAPFile <- function(param.list = NA, output.file ="~/RayTracerFortran/RayTracer/test_MAP.dat")
{
  if (is.na(param.list))
  {
  message('Nothing was passed as a param list!') 
  return(NA)
  }
  NVoro <- param.list$MAP.k * param.list$NPL
  n.entries <- param.list$NLMX * param.list$NPL * param.list$NPL 
  MAP.row = rep(0.0,2*n.entries)
  map.vd<- rep(NA,param.list$MAP.k)
  # Make Voro Column names - this will be my velocity profiles
  map.vd[seq(1,NVoro-1,by=2)] <- param.list$MAP.d
  map.vd[seq(2,NVoro,by=2)] <- param.list$MAP.v

  MAP.row[1] <- param.list$MAP.k
  MAP.row[2:(2*param.list$MAP.k+2)] <- map.vd
  MAP.row[(2*param.list$MAP.k+2)] <- param.list$MAP.sd
  
  MAP.row[(param.list$NPL * param.list$NLMX+2)] <- param.list$MAP.sd
  
  fmt <- paste0(c('%d',rep(' %3.2f',2*n.entries-1)) ,collapse = '')
  MAP.row.char <- do.call(sprintf, c(list(fmt), MAP.row))
  
  
  fileContent <- MAP.row.char
  con <- file(output.file, "wb")
  writeChar(fileContent, con, nchar(fileContent), eos = NULL)
  close(con)
  
  return(MAP.row.char)
  
  
}


drawProfiles.TD <- function(x,NLayers=NA,NLMax = 10)
{
  d=as.numeric(x[(NLMax+1) : ( NLMax+NLayers )])
  v=as.numeric(x[1:(NLayers)])
  f <- approxfun(d,v,method = "constant")
  dnew <- seq(1,8000,length.out = 200)
  vnew <- f(dnew)
  vnew[is.na(vnew)] = v[length(v)]
  dd <- d[-1]
  for (i in seq(length(dd)) )
  {
    ptsX=seq(from=min(v[i],v[i+1]),to = max(v[i],v[i+1]),by=40)
    ptsY=rep(dd[i],length(ptsX))
    dnew=c(dnew,ptsY)
    vnew=c(vnew,ptsX)
  }
  prof <-tbl_df(data.frame(d=dnew,v=vnew))
  return(prof)
}


