library(tidyverse)
library(ncdf4)
library(fields)
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

#Extract model output
File <- "Euler.nc"
NO3_MOD <- ncread(File, 'NO3')
Chl_MOD <- ncread(File, 'CHL')
NPP_MOD <- ncread(File, 'NPP')

#Obtain final year
NO3_MOD <- NO3_MOD[,(ncol(NO3_MOD)-365+1):ncol(NO3_MOD)]
Chl_MOD <- Chl_MOD[,(ncol(Chl_MOD)-365+1):ncol(Chl_MOD)]
NPP_MOD <- NPP_MOD[,(ncol(NPP_MOD)-365+1):ncol(NPP_MOD)]
  
newdata <- expand.grid(
  Depth = Z_r,
  DOY = 1:365
)

#Obtain BATS NO3 data
bats_NO3 <- read.csv('bats_NO3.csv') |>
  select(Depth, Date, NO3) |>
  filter(NO3 > 0) |>
  mutate(Date = as_date(ymd(Date))) |>
  mutate(DOY = as.numeric(strftime(Date, format='%j')))

#Model Chl as a function of Depth and DOY using knn
NO3_obs <- bats_NO3 |>
  mutate(Depth = -Depth)

#Estimate Chl using Knn
library(class)
NO3.knn <- knn( train = NO3_obs[,c('Depth', 'DOY')],
                test  = newdata[,c('Depth', 'DOY')],
                cl    = as.factor(as.character(NO3_obs[,'NO3'])),
                k     = 3)

newdata$NO3_knn <- as.numeric(as.character(NO3.knn))
NO3knn <- matrix(newdata$NO3_knn, nr = length(Z_r),
                nc = 365)

#Obtain new BATS Chl data
bats_pigment <- read.csv('bats_pigments.csv') |>
  select(Depth, Date, p16_Chl) |>
  mutate(Date = as_date(Date)) |>
  mutate(DOY = as.numeric(strftime(Date, format='%j'))) |>
  filter(p16_Chl > 0)

DOY_mth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
DOY_mth <- cumsum(DOY_mth)

#Model Chl as a function of Depth and DOY using knn
Chl_obs <- bats_pigment |>
  mutate(Depth = -Depth)

#Estimate Chl using Knn
Chl.knn <- knn( train = Chl_obs[,c('Depth', 'DOY')],
                test  = newdata[,c('Depth', 'DOY')],
                cl    = as.factor(as.character(Chl_obs[,'p16_Chl'])),
                k     = 3)

newdata$Chl_knn <- as.numeric(as.character(Chl.knn))
Chlknn <- matrix(newdata$Chl_knn, nr = length(Z_r),
                nc = 365)
Chlmax <- 1

#Obtain NPP data at BATS
NPP <- read.csv('BATS_Primary_Production.csv') |>
  select(depth, time_out, pp) |>
  mutate(Date = as_date(ymd_hm(time_out))) |>
  mutate(DOY = as.numeric(strftime(Date, format='%j'))) |>
  filter(pp > 0)

#Model Chl as a function of Depth and DOY using knn
NPP_obs <- NPP  |>
  mutate(Depth = -depth)

#Estimate NPP using Knn
NPP.knn <- knn( train = NPP_obs[,c('Depth', 'DOY')],
                test  = newdata[,c('Depth', 'DOY')],
                cl    = as.factor(as.character(NPP_obs[,'pp'])),
                k     = 3)

newdata$NPP_knn <- as.numeric(as.character(NPP.knn))
NPPknn <- matrix(newdata$NPP_knn, 
                 nr = length(Z_r),
                 nc = 365)

NPPmax <- 10
library(plot3D)
pdf('BATS_obs_mod_comparison.pdf', width = 6, height = 9)
par(mfrow=c(3,2),
    mar = c(3,2.5,2,2),
    oma = c(3,3,2,1))
NO3max <- 3
NO3knn[NO3knn > NO3max] <- NO3max
image2D(t(NO3knn), 1:365, -Z_r,
        zlim = c(0, NO3max),
        xlab = '',
        ylab = ''
        )
mtext('(a) DIN', adj=0)

NO3_MOD[NO3_MOD > NO3max] <- NO3max
image2D(t(NO3_MOD), 1:365, -Z_r,
        zlim = c(0, NO3max),
        xlab = '',
        ylab = ''
        )
mtext('(b) DIN', adj=0)

Chlknn[Chlknn > Chlmax] <- Chlmax
image2D(t(Chlknn), 1:365, -Z_r,
        zlim = c(0, Chlmax),
        xlab = '',
        ylab = ''
        )
mtext('(c) Chl', adj=0)

Chl_MOD[Chl_MOD > Chlmax] <- Chlmax
image2D(t(Chl_MOD), 1:365, -Z_r,
        zlim = c(0, Chlmax),
        xlab = '',
        ylab = ''
        )
mtext('(d) Chl', adj=0)

NPPknn[NPPknn > NPPmax] <- NPPmax
image2D(t(NPPknn), 1:365, -Z_r,
        zlim = c(0, NPPmax),
        xlab = '',
        ylab = ''
        )
mtext('(e) NPP', adj=0)

NPP_MOD[NPP_MOD > NPPmax] <- NPPmax
NPP_MOD[NPP_MOD < 0] <- 0
image2D(t(NPP_MOD), 1:365, -Z_r,
        zlim = c(0, NPPmax),
        xlab = '',
        ylab = ''
        )
mtext('(f) NPP', adj=0)
mtext('Observation', side = 3, outer=T, adj = .25)
mtext('Model',       side = 3, outer=T, adj = .75)
mtext('Day',       side = 1, outer=T)
mtext('Depth (m)', side = 2, outer=T)
dev.off()
