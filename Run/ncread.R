#Simple function to read .nc data:
library(ncdf4)
ncread  <- function(file, VAR, start = NA, count = NA){
  if(file.exists(file)){
    nc    <- nc_open(file)
  }else{
    stop(paste(file,'does not exist!'))
  }
  nvar    <- which(names(nc$var) == VAR)
  if (length(nvar) < 1) stop(paste0(VAR, ' not found in ',file))
    v4    <- nc$var[[nvar]] # The index of data to be extracted
  data    <- ncvar_get(nc, v4, start = start, count = count)
  nc_close(nc)
  return(data)
}