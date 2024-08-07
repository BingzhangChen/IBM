source('ncread.R')
library(tidyverse)
library(ggpubr)
library(scales)

par_file   <- 'ParY6.nc'
Euler_file <- 'Euler.nc'

#Define a community
#Obtain Z_r and Z_w######
nc   <- nc_open(Euler_file)
Z_r  <- nc$dim$Z_r$vals
Z_w  <- nc$dim$Z_w$vals
Days <- nc$dim$Time$vals
ZOO0 <- nc$dim$ZOO$vals #Zooplantkon concentration
nc_close(nc)
nlev <- length(Z_r) #Number of vertical layers
Hz <- diff(Z_w) #Heights of each layer

#ESD of zooplankton
ZESD_min <- .8   #Minimal zooplankton ESD
ZESD_max <- 3600 #Maximal zooplankton ESD
NZOO     <- 20 #Number of zooplankton size classes
ZOOESD   <- seq(from=log(ZESD_min), 
              to=log(ZESD_max), 
              length.out=NZOO)

dESD <- diff(ZOOESD)[1]

#upper limit of ESD intervals
ESD_up <- ZOOESD+dESD/2

#Convert ESD to biovolume (um^3)
Vols <- pi/6*exp(ZOOESD)**3

#Volume upper limit
Vol_up <- pi/6*exp(ESD_up)**3

#Compute abundance of each zooplankton size class
#Following Harris (2000) to convert volume to carbon (pg)
ZOO_V2C <- function(V, pd=.2, pc=.45, rho=1.025) rho*V*pd*pc
redfield <- 106/16

#Select the surface grid in March
DOY_target <- c(90, 240)

#Extract carbon phytoplankton particles
DOY  <- ncread(par_file, 'DOY')
Hour <- ncread(par_file, 'Hour')

#Calculate number of days
NT <- length(unique(DOY))

#Extract phytoplankton super-individuals
cell_N  <- ncread(par_file, 'N_cell')
cell_IZ <- ncread(par_file, 'IZ')
cell_CDiv <- ncread(par_file, 'CDiv')

for (k in 1:length(DOY_target)){
  
  #Extract zooplankton
  wt  <- which(Days %% 365 == DOY_target[k])   
  ZOO <- ZOO0[,,wt[length(wt)]] 
  ZOO <- ZOO[,nlev]
  
  #Abundance of zooplankton
  Azoo <- ZOO*redfield *12/ZOO_V2C(Vols)*1e9 #Unit: ind/m3 

  wt <- which(DOY_target[k]==DOY)
  wt <- wt[1]
  
  IZ1 <- cell_IZ[, wt]
  
  wz <-  which(IZ1 == nlev)
  
  CDiv1 <- cell_CDiv[wz, wt] #Unit: pmol C cell-1
  
  #Convert CDiv to volume
  C2V <- function(CDiv, a=-.69, b=.88) (12*CDiv/10**a)**(1/b) 
  
  pVol <- C2V(CDiv1) #Cell volume of the phyto. particles
  N1 <- cell_N[wz, wt]
  
  #Count the total number of abundances within each size interval
  Aphy <- numeric(length(ZOOESD))
  
  for (i in 1:length(Aphy)){
    if (i == 1){
      wi = which(pVol <= Vol_up[i])
    }else{
      wi = which(pVol <= Vol_up[i] & pVol > Vol_up[i-1])
    }
    Aphy[i] <- sum(N1[wi])
  }
  
  #Convert to number of cells per m3
  Aphy <- Aphy/Hz[nlev]
  
  #Establish dataframe to hold phyto. and zoo. abundances
  #and Convert SS to the long-format
  SS <- data.frame(Vol = Vols,
                  APhy= Aphy,
                  AZoo= Azoo) |>
    mutate(Atot=APhy+AZoo) |>
    pivot_longer(
      cols = starts_with('A'),
      names_to = 'Group',
      names_prefix = 'A',
      values_to = 'Abun'
    )
  
  #Plot abundance of zooplankton against volume
  SS1 <- SS |>
    filter(Abun > 0) |>
    filter((Group == 'Zoo' & Vol > 3) | (Group == 'Phy')) 
  
  SS2 <- SS |>
    filter(Abun > 0) |>
    filter((Group == 'Zoo' & Vol > 3) | (Group == 'Phy') | (Group == 'tot')) 
  
  #Obtain the linear regression for phyto
  phy.lm <- lm(log10(Abun) ~ log10(Vol), SS1 |> filter(Group == 'Phy'))
  
  #Extract regression slope
  slope.phy <- round(coef(summary(phy.lm))[2,1],2)
  
  zoo.lm <- lm(log10(Abun) ~ log10(Vol), SS1 |> filter(Group == 'Zoo'))
  slope.zoo <- round(coef(summary(zoo.lm))[2,1],2)
  
  tot.lm <- lm(log10(Abun) ~ log10(Vol), SS2 |> filter(Group == 'tot'))
  slope.tot <- round(coef(summary(tot.lm))[2,1], 2)
  
  if (k == 1){
    p1 <- SS1 |>
      ggplot(aes(Vol, Abun, colour = Group, linetype=Group)) +
      geom_point() +
      geom_smooth(method = lm) + 
      xlim(1, 1e11) +
      ylim(0.1, 1e12) +
      scale_x_log10(breaks=breaks_log(base = 10),
                    labels=label_log(base = 10)) +
      scale_y_log10(breaks=breaks_log(base = 10),
                    labels=label_log(base = 10)) +
      xlab('') +
      ylab(bquote('Abundance (inds '*m^-3*')')) +
      annotate('text', label=paste('Slope=', slope.phy),
               x=10000, y =1e9, colour='red') +
      annotate('text', label=paste('Slope=', slope.zoo),
               x=1e6, y =1e2, colour='#00BF7D') +
      theme(legend.position = "top") +
      theme_light()
  }else{
     p2 <- SS1 |>
      ggplot(aes(Vol, Abun, colour = Group, linetype=Group)) +
      geom_point() +
      geom_smooth(method = lm) + 
      xlim(1, 1e11) +
      ylim(0.1, 1e12) +
      scale_x_log10(breaks=breaks_log(base = 10),
                    labels=label_log(base = 10)) +
      scale_y_log10(breaks=breaks_log(base = 10),
                    labels=label_log(base = 10)) +
      xlab(bquote('Volume ('*mu*m^3*')')) +
      ylab(bquote('Abundance (inds '*m^-3*')')) +
      annotate('text', label=paste('Slope=', slope.phy),
               x=10000, y =1e9, colour='red') +
      annotate('text', label=paste('Slope=', slope.zoo),
               x=1e6, y =1e2, colour='#00BF7D') +
      theme(legend.position="none") + 
      theme_light()  
  }
}

library(patchwork)
p1/p2 +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(',
                  tag_suffix = ')')

ggsave('Fig10Plankton_size_spectra.pdf', width = 4, height = 7)
