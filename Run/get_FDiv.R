library(TPD)
library(tidyverse)
source('ncread.R')

#Extract particle traits and compute daily integrated functional trait diversity  
get_Fdiv <- function(par_file){
  DOY  <- ncread(par_file, 'DOY')
  Hour <- ncread(par_file, 'Hour')
  
  #Calculate number of days
  NT <- length(unique(DOY))

  #Calculate Functional richness of each date for the whole water column
  FRic <- numeric(NT)
  FEve <- FRic
 
  cell_Topt  <- ncread(par_file, 'Topt')
  cell_CDiv  <- ncread(par_file, 'CDiv')
  cell_alpha <- ncread(par_file, 'alpha')
  cell_PC    <- ncread(par_file, 'PC')
  cell_N     <- ncread(par_file, 'N_cell')
  cell_IZ    <- ncread(par_file, 'IZ')
  
  tindex <- 1:length(DOY)
  for (i in 1:NT){
    wt  <- DOY==unique(DOY)[i]
    wt  <- tindex[wt]
    wt  <- wt[1]
      
    Topt1 <- cell_Topt[, wt]
    CDiv1 <- cell_CDiv[, wt]
    alpha1<- cell_alpha[, wt]
    PC1   <- cell_PC[, wt]
    N1    <- cell_N[, wt]
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
      
      FRic[i] <- RED$species$FRichness[[1]]
      FEve[i] <- RED$species$FEvenness[[1]]
    }
  }
  return(list(FRic = FRic, FEve = FEve))
}
