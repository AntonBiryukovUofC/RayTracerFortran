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




