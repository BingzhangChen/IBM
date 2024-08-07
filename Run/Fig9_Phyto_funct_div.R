par_file <- 'ParY6.nc'
source('ncread.R')
library(TPD)
library(tidyverse)

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

#Calculate Functional richness for each date and each layer
FRic <- matrix(NA, nr = NT, nc = nlev)
FEve <- FRic

#Extract particle traits
cell_Topt  <- ncread(par_file, 'Topt')
cell_CDiv  <- ncread(par_file, 'CDiv')
cell_alpha <- ncread(par_file, 'alpha')
cell_PC    <- ncread(par_file, 'PC')
cell_N     <- ncread(par_file, 'N_cell')
cell_IZ    <- ncread(par_file, 'IZ')

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
    ) 
    
    cells <- cells |>
              mutate(TC = PC * N/10e6) #unit: umol C ind-1
    
    #Compute TRic and Teve for this time point
    sps <- rep(paste0('Sp', 1), nrow(cells)) #Assume one single species
    traits_PHY <- cells |>
      select(Topt, alpha, CDiv)
    
    if (nrow(unique(traits_PHY)) > 3){
      TPDs_PHY <- TPDs(species = sps, 
                       traits_PHY, 
                       weight = cells$TC)
          
      #Calculate functional diversity of this community
      RED <- REND(TPDs = TPDs_PHY)
      
      FRic[i,j] <- RED$species$FRichness[[1]]
      FEve[i,j] <- RED$species$FEvenness[[1]]
    }
  }
}
save(FRic, FEve, file = 'Func_traits.Rdata')

load('Func_traits.Rdata')

#imputation
FD <- expand.grid(DOY = 1:NT, Z = -Z_r)
FD <- FD |>
  mutate(FRic = as.vector(FRic)) |>
  mutate(FEve = as.vector(FEve))
  
library(mice)
FD_imputed <- mice(FD, method = "cart")
FD$FRic_imp <- complete(FD_imputed)$FRic 
FD$FEve_imp <- complete(FD_imputed)$FEve 

FRic <- matrix(FD$FRic_imp,
               nr = NT,
               nc = nlev)
FEve <- matrix(FD$FEve_imp,
               nr = NT,
               nc = nlev)
FRic_max <- 50
FRic[FRic > FRic_max] <- FRic_max
FEve_max <- .7
FEve[FEve > FEve_max] <- FEve_max

#Plot out
library(plot3D)
pdf('Func_Div.pdf', width = 5, height = 8)
par(mfrow=c(2,1),
    mar = c(2,2.5,2,2),
    oma = c(3,3,2,1))
image2D(FRic, 1:NT, -Z_r,
        xlab = '',
        ylab = '',
        zlim = c(0, FRic_max)
        )
mtext('(a) Functional richness', adj=0)

image2D(FEve, 1:NT, -Z_r,
        xlab = '',
        ylab = ''
        )
mtext('(b) Functional evenness', adj=0)
mtext('Day',       side = 1, outer=T)
mtext('Depth (m)', side = 2, outer=T)
dev.off()
#head(TPDs_PHY$data$evaluation_grid)
#nrow(TPDs_PHY$data$evaluation_grid)
#names(TPDs_PHY$TPDs)


#How to calculate the annual mean FEve?

#Check the total probability
#sum(TPDs_PHY$TPDs$Sp1)
