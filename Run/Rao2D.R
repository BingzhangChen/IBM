par_file <- 'ParY6.nc'
source('ncread.R')
library(tidyverse)
library(Rcpp)
sourceCpp('Rao.cpp')

#Obtain Z_r and Z_w######
nc   <- nc_open('Euler.nc')
Z_r  <- nc$dim$Z_r$vals
Z_w  <- nc$dim$Z_w$vals
Days <- nc$dim$Time$vals
nc_close(nc)
nlev <- length(Z_r) #Number of vertical layers
#########################

DOY  <- ncread(par_file, 'DOY')
Hour <- ncread(par_file, 'Hour')

#Calculate number of days
NT <- length(unique(DOY))

#Calculate Rao index for each date and each layer
Rao2D <- matrix(NA, nr = NT, nc = nlev)

#Extract particle traits
cell_Topt  <- ncread(par_file, 'Topt')
cell_CDiv  <- ncread(par_file, 'CDiv')
cell_alpha <- ncread(par_file, 'alpha')
cell_PC    <- ncread(par_file, 'PC')
cell_N     <- ncread(par_file, 'N_cell')
cell_IZ    <- ncread(par_file, 'IZ')

#Compute minimal and maximal traits
Topt_max <- max(cell_Topt)
Topt_min <- min(cell_Topt)
dTopt    <- Topt_max-Topt_min
CDiv_max <- max(cell_CDiv)
CDiv_min <- min(cell_CDiv)
dCDiv    <- log(CDiv_max)-log(CDiv_min)
alpha_max<- max(cell_alpha)
alpha_min<- min(cell_alpha)
dalpha   <- alpha_max - alpha_min
 
tindex <- 1:length(DOY)
for (i in 1:NT){
  wt  = DOY==unique(DOY)[i]
  wt  = tindex[wt]
  wt  = wt[1]
  IZ1 = cell_IZ[, wt]
  
  for (j in 1:nlev){
    
    #Extract the indeces corresponding to cells within each layer
    wz = IZ1 == j
    
    Topt1 <- cell_Topt[wz, wt]
    CDiv1 <- cell_CDiv[wz, wt]
    alpha1<- cell_alpha[wz, wt]
    PC1   <- cell_PC[wz, wt]
    N1    <- cell_N[wz, wt]
    cells <- data.frame(
      Topt  = Topt1,
      alpha = alpha1,
      CDiv  = log(CDiv1),
      PC    = PC1,
      N     = N1
    ) %>% 
      mutate(TC = PC * N/10e6)  %>% #unit: umol C ind-1
      mutate(wTC= TC/sum(.$TC)) %>% 
      mutate(Topt_norm = (Topt-Topt_min)/dTopt) %>% 
      mutate(CDiv_norm = (CDiv-log(CDiv_min))/dCDiv) %>% 
      mutate(alpha_norm= (alpha-alpha_min)/dalpha)
    
    #Retrieve only the trait data
    traits  <- cells  %>%
      select(CDiv_norm, Topt_norm, alpha_norm) %>% 
      as.matrix()
    
    #Calculate  Rao index of this community
    Rao2D[i,j] <- RaoCpp(traits, w = cells$wTC)
  }
}
save(Rao2D, file = 'Rao2D.Rdata')

write.csv(Rao2D, file='Rao2D.csv', row.names = F, col.names = F)

library(plot3D)
pdf('FIG9_Rao.pdf', width = 5, height = 5)
par(mfrow=c(1,1),
    mar = c(4,4,2,2),
    oma = c(3,3,2,1))
image2D(Rao2D, 1:NT, -Z_r,
        xlab = 'Days',
        ylab = 'Depth (m)'
        )
dev.off()
