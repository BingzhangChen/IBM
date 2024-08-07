library(tidyverse)
library(ncdf4)
library(fields)
library(plot3D)
ncread  <- function(file, VAR, start = NA, count = NA){
  if(file.exists(file)){
    nc    <- nc_open(file)
  }else{
    stop(paste(file,'does not exist!'))
  }
  data    <- ncvar_get(nc, VAR, start = start, count = count)
  nc_close(nc)
  return(data)
}

File <- 'Euler.nc'

#Obtain phytoplankton carbon
 PC_MOD <- ncread(File, 'PC')
 PN_MOD <- ncread(File, 'PN')
Chl_MOD <- ncread(File, 'CHL')

#Obtain final year
Chl_MOD <- Chl_MOD[,(ncol(Chl_MOD)-365+1):ncol(Chl_MOD)]
 PC_MOD <-  PC_MOD[,(ncol( PC_MOD)-365+1):ncol(PC_MOD)]
 PN_MOD <-  PN_MOD[,(ncol( PN_MOD)-365+1):ncol(PN_MOD)]

#Calculate QN and theta
theta <- Chl_MOD/PC_MOD
QN <- PN_MOD/PC_MOD

nc   <- nc_open(File)
Z_r  <- nc$dim$Z_r$vals
nc_close(nc)

pdf('Fig7_QN_theta.pdf', width = 8, height = 5)
par(mfrow=c(1,2),
    mar = c(2,2.5,2,2),
    oma = c(2.5,2.5,2,1))
image2D(t(theta), 1:365, -Z_r,
        xlab = '',
        ylab = '')
mtext('(a) Chl:C', adj=0)

image2D(t(QN), 1:365, -Z_r,
        xlab = '',
        ylab = '')
mtext('(b) N:C', adj=0)
mtext('Day',       side = 1, outer=T)
mtext('Depth (m)', side = 2, outer=T)
dev.off()
