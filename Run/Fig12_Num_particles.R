source('get_FDiv.R')
Euler_20K <- 'Euler.nc'
Euler_5K  <- '../Run5000/Euler.nc'
Euler_10K <- '../Run10000/Euler.nc'
Euler_50K <- '../Run50K/Euler.nc'

par_20K <- 'ParY6.nc'
par_5K  <- '../Run5000/ParY6.nc'
par_10K <- '../Run10000/ParY6.nc'
par_50K <- '../Run50K/ParY6.nc'

#####Compute the 1D integrated functional richness and evenness############
FDiv_50K <- get_Fdiv(par_50K)
FDiv_20K <- get_Fdiv(par_20K)
FDiv_10K <- get_Fdiv(par_10K)
FDiv_5K  <- get_Fdiv(par_5K)

save(FDiv_5K, FDiv_10K, FDiv_20K, file='FDiv_num.Rdata')

#############Obtain 1D integrated total DIN, PN, PC, Chl, ZOO, and NPP
get_Euler_final <- function(File){
  nc   <- nc_open(File)
  Z_r  <- nc$dim$Z_r$vals
  Z_w  <- nc$dim$Z_w$vals
  Days <- nc$dim$Time$vals
  ZOO  <- nc$dim$ZOO$vals
  nc_close(nc)
  
  #Number of days
  N_Days <- length(Days)
  
  #Compute Hz
  Hz   <- diff(Z_w)
 
  NO3 <- ncread(File, 'NO3')
  PC  <- ncread(File, 'PC')
  PN  <- ncread(File, 'PN')
  CHL <- ncread(File, 'CHL')
  NPP <- ncread(File, 'NPP')
  
  #Retrieve final year
  wt   <- (N_Days-365+1):N_Days
  NO3  <- NO3[, wt]
  PC   <- PC[, wt]
  PN   <- PN[, wt]
  Chl  <- CHL[, wt]
  NPP  <- NPP[, wt]
  ZOO  <- ZOO[, , wt]
  TZOO <- apply(ZOO, c(2,3), sum) #Calculate total zooplankton biomass
  
  #Sum up all Eulerian fields
  Tnpp  <- rep(NA, length(wt))
  T_PC  <- Tnpp
  T_NO3 <- Tnpp
  T_Chl <- T_NO3
  T_PN  <- T_NO3
  T_ZOO <- T_NO3
  
  for (j in 1:length(wt)){
    #sum up across depth
    Tnpp[j] <- sum(NPP[,j] * Hz)
    T_PC[j] <- sum(PC[,j] * Hz)
    T_NO3[j]<- sum(NO3[,j] * Hz)
    T_PN[j] <- sum(PN[,j] * Hz)
    T_ZOO[j] <- sum(TZOO[,j] * Hz)
    T_Chl[j] <- sum(Chl[,j] * Hz)
  }
  
  return(list(Tnpp=Tnpp, 
              T_PC=T_PC, 
              T_NO3=T_NO3, 
              T_PN=T_PN,
              T_ZOO=T_ZOO,
              T_Chl=T_Chl
              )) 
}

#Euler_var_50K <- get_Euler_final(Euler_50K)
Euler_var_20K <- get_Euler_final(Euler_20K)
Euler_var_10K <- get_Euler_final(Euler_10K)
Euler_var_5K  <- get_Euler_final(Euler_5K)

###Plot three cases (different number of super-individuals together)
pdf('comp_num_sind.pdf', width = 6, height = 9)
par(mfrow=c(4,2),
    mar = c(2,4.8,2,1),
    oma = c(3,3,2,1),
    cex.lab = 1.5,
    cex.axis= 1.5,
    lwd=2)

Days <- 1:365
#Total DIN against time for three numbers of super-individuals
T_NO3_min <- min(c(Euler_var_5K$T_NO3, Euler_var_10K$T_NO3, Euler_var_20K$T_NO3))
T_NO3_max <- max(c(Euler_var_5K$T_NO3, Euler_var_10K$T_NO3, Euler_var_20K$T_NO3))
plot(Days, Euler_var_5K$T_NO3, 
     type = 'l',
     ylim = c(T_NO3_min, T_NO3_max),
     xlab = '',
     ylab = expression(paste('DIN (mmol '*m^-2*')')))
points(Days, Euler_var_10K$T_NO3, type='l', lty=2, col=2)
points(Days, Euler_var_20K$T_NO3, type='l', lty=3, col=3)
#points(Days, Euler_var_50K$T_NO4, type='l', lty=4, col=4)
mtext('(a) DIN', adj = 0)

T_PC_min <- min(c(Euler_var_5K$T_PC, Euler_var_10K$T_PC, Euler_var_20K$T_PC))
T_PC_max <- max(c(Euler_var_5K$T_PC, Euler_var_10K$T_PC, Euler_var_20K$T_PC))
plot(Days, Euler_var_5K$T_PC, 
     type = 'l',
     ylim = c(T_PC_min, T_PC_max),
     xlab = '',
     ylab = expression(paste(P[C]*' (mmol '*m^-2*')')))
points(Days, Euler_var_10K$T_PC, type='l',lty=2, col=2)
points(Days, Euler_var_20K$T_PC, type='l',lty=3, col=3)
mtext('(b) Phyto. carbon', adj = 0)
legend('topright',
       c('5000', '10000', '20000'),
       lty = 1:3,
       col = 1:3)

T_PN_min <- min(c(Euler_var_5K$T_PN, Euler_var_10K$T_PN, Euler_var_20K$T_PN))
T_PN_max <- max(c(Euler_var_5K$T_PN, Euler_var_10K$T_PN, Euler_var_20K$T_PN))
plot(Days, Euler_var_5K$T_PN, 
     type = 'l',
     ylim = c(T_PN_min, T_PN_max),
     xlab = '',
     ylab = expression(paste(P[N]*' (mmol '*m^-2*')')))
points(Days, Euler_var_10K$T_PN, type='l',lty=2, col=2)
points(Days, Euler_var_20K$T_PN, type='l',lty=3, col=3)
mtext('(c) Phyto. nitrogen', adj = 0)

T_Chl_min <- min(c(Euler_var_5K$T_Chl, Euler_var_10K$T_Chl, Euler_var_20K$T_Chl))
T_Chl_max <- max(c(Euler_var_5K$T_Chl, Euler_var_10K$T_Chl, Euler_var_20K$T_Chl))
plot(Days, Euler_var_5K$T_Chl, 
     type = 'l',
     ylim = c(T_Chl_min, T_Chl_max),
     xlab = '',
     ylab = expression(paste('Chl (mg '*m^-2*')')))
points(Days, Euler_var_10K$T_Chl, type='l',lty=2, col=2)
points(Days, Euler_var_20K$T_Chl, type='l',lty=3, col=3)
mtext('(d) Chl', adj = 0)

Tnpp_min <- min(c(Euler_var_5K$Tnpp, Euler_var_10K$Tnpp, Euler_var_20K$Tnpp))
Tnpp_max <- max(c(Euler_var_5K$Tnpp, Euler_var_10K$Tnpp, Euler_var_20K$Tnpp))
plot(Days, Euler_var_5K$Tnpp, 
     type = 'l',
     ylim = c(Tnpp_min, Tnpp_max),
     xlab = '',
     ylab = expression(paste('NPP (mg C '*m^-2*' '*d^-1*')')))
points(Days, Euler_var_10K$Tnpp, type='l',lty=2, col=2)
points(Days, Euler_var_20K$Tnpp, type='l',lty=3, col=3)
mtext('(e) NPP', adj = 0)

T_ZOO_min <- min(c(Euler_var_5K$T_ZOO, Euler_var_10K$T_ZOO, Euler_var_20K$T_ZOO))
T_ZOO_max <- max(c(Euler_var_5K$T_ZOO, Euler_var_10K$T_ZOO, Euler_var_20K$T_ZOO))
plot(Days, Euler_var_5K$T_ZOO, 
     type = 'l',
     ylim = c(T_ZOO_min, T_ZOO_max),
     xlab = '',
     ylab = expression(paste('ZOO (mmol N '*m^-2*')')))
points(Days, Euler_var_10K$T_ZOO, type='l', lty=2,col=2)
points(Days, Euler_var_20K$T_ZOO, type='l', lty=3,col=3)
mtext('(f) Zooplankton', adj = 0)

FRic_min <- min(min(FDiv_5K$FRic), 
                min(FDiv_10K$FRic), 
                min(FDiv_20K$FRic))
FRic_max <- max(c(FDiv_5K$FRic, FDiv_10K$FRic, FDiv_20K$FRic))
plot(Days, FDiv_5K$FRic, 
     type = 'l',
     ylim = c(FRic_min, FRic_max),
     xlab = '',
     ylab = 'Functional richness')
points(Days, FDiv_10K$FRic, type='l', lty=2, col=2)
points(Days, FDiv_20K$FRic, type='l', lty=3, col=3)
mtext('(g) Phyto. functional richness', adj = 0)

FEve_min <- min(min(FDiv_5K$FEve), 
                min(FDiv_10K$FEve), 
                min(FDiv_20K$FEve))
FEve_max <- max(c(FDiv_5K$FEve, FDiv_10K$FEve, FDiv_20K$FEve))
plot(Days, FDiv_5K$FEve, 
     type = 'l',
     ylim = c(FEve_min, FEve_max),
     xlab = '',
     ylab = 'Functional evenness')
points(Days, FDiv_10K$FEve, type='l', lty=2, col=2)
points(Days, FDiv_20K$FEve, type='l', lty=3, col=3)
mtext('(h) Phyto. functional evenness', adj = 0)

mtext('Days', side = 1, outer = T, line = .1, adj = .5)
dev.off()
