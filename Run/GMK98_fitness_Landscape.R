library(tidyverse)
library(patchwork)
source('GMK98_trait_functions.R')

Season <- 'summer'

if (Season == 'summer'){
  #summer surface
  Temp_ <- 28
  NO3_ <- 0.1
  PAR_ <- 250
}else if (Season == 'winter'){
  #winter surface
  Temp_ <- 18
  NO3_ <- 0.4
  PAR_ <- 150 
}

Topt     <- seq(0, 30,length.out=20)
alphaChl <- exp(seq(-3, -1, length.out=20))
ESD <- exp(seq(log(0.5), log(100), length.out = 20))
Vol <- pi/6*ESD^3
a <- 10^(-0.69)
b <- 0.88
PHY_CDiv <-  a*Vol^b /12 #pmol C cell-1
pGrid <- expand.grid(Topt = Topt, CDiv = PHY_CDiv, alphaChl = alphaChl)

C2ESD <- function(pC){
  Vol <- ((pC*12)/a)^(1/b)
  ESD <- (Vol*6/pi)^.33333
  return(ESD)
}

#Fix traits (Topt, alphaChl, CDiv) to compute the growth rate (dC/C, dN/N, and dChl/Chl)

#Function to estimate steady state of growth rate, N:C, and Chl:C for each environment
#and trait combination
GMK_ss <- function(Topt, CDiv, alphaChl, Temp = Temp_, NO3 = NO3_, PAR = PAR_){
  PHY_C <- CDiv/2 #pmol C cell-1
  PHY_N <- PHY_C/(106/16)
  PHY_Chl <- PHY_C * 12 /50 #pgChl cell-1
 
  dtdays <- 5/(60*24) #Timestep: 5 mins
  
  #Simulate for maximally 100 days
  Tot <- 100
  tol <- 1e-6 #numerical error that can be tolerated
  
  Nstep <- Tot/dtdays
  out <- matrix(NA, nr = Nstep, nc = 3)
  for (i in 1:Nstep){
   cff <- GMK98_TempSizeLight(Temp, PAR, NO3, Topt, PHY_C, 
                              PHY_N, PHY_Chl, CDiv, alphaChl)
   out[i, 1] <- cff[1]/PHY_C #Carbon-specific growth rate
   out[i, 2] <- PHY_N/PHY_C #N:C ratio
   out[i, 3] <- PHY_Chl/(PHY_C*12) #Chl:C ratio
   
   if (i > 1 &&
       abs(out[i, 1] - out[i-1,1]) < tol &&
       abs(out[i, 2] - out[i-1,2]) < tol &&
       abs(out[i, 3] - out[i-1,3]) < tol) break
   
   PHY_C   <- PHY_C   + dtdays*cff[1]
   PHY_N   <- PHY_N   + dtdays*cff[2]
   PHY_Chl <- PHY_Chl + dtdays*cff[3]
  }
  return(out[i,])
}

pGrid <- pGrid |>
  mutate(mu = NA) |>
  mutate(QN = NA) |>
  mutate(theta = NA)

#Compute the growth rate, QN and theta at steady state for each condition
for (k in 1:nrow(pGrid)){
  z <- GMK_ss(pGrid$Topt[k], 
              pGrid$CDiv[k], 
              pGrid$alphaChl[k])
  
  pGrid$mu[k]    <- z[1]
  pGrid$QN[k]    <- z[2]
  pGrid$theta[k] <- z[3]
}

out_name <- paste0('phyto_fitness_landscape_', 
                   Season,
                   '.Rdata')
save(pGrid, file = out_name)

#Plot growth rate
load(out_name)

#3 X 3 Graph
pGrid <- pGrid |>
  mutate(ESD = C2ESD(CDiv))

#Average responses across alphaChl
#pGrid1 <- 
#  pGrid |>
#  group_by(Topt, ESD) |>
#  summarise(mu = mean(mu),
#            QN = mean(QN),
#            C2Chl = mean(C2Chl))

p1 <- pGrid |>
  filter(alphaChl == alphaChl[10]) |>
  ggplot(aes(x = Topt, y = ESD, z = mu)) +
  geom_contour_filled(show.legend = F) +
  scale_y_log10() +
  xlab(bquote(T[opt] * ' (ºC)')) +
  ylab(bquote('Size (µm)')) + 
  annotate('text', x = 0, y = 120, 
    label = bquote('A) µ ('*d^-1*')'),
    hjust = 0) +
  theme_light()

p2 <- pGrid |>
  filter(Topt == Topt[length(Topt)]) |>
  ggplot(aes(x = ESD, y = alphaChl, z = mu)) +
  geom_contour_filled(show.legend = F) +
  scale_x_log10() +
  scale_y_log10() +
  xlab(bquote('Size (µm)')) + 
  ylab(bquote(alpha[Chl] * ' ('* m^2*' '*W^-1*' molC ' * gChl^-1*' '*d^-1*')')) +
  annotate('text', x = .6, y = .4, 
    label = bquote('B) µ ('*d^-1*')'),
    hjust = 0) +
  theme_light()

p3 <- pGrid |>
  filter(ESD == ESD[5]) |>
  ggplot(aes(x = Topt, y = alphaChl, z = mu)) +
  geom_contour_filled() +
  scale_y_log10() +
  xlab(bquote(T[opt] * ' (ºC)')) +
  ylab(bquote(alpha[Chl] * ' ('* m^2*' '*W^-1*' molC ' * gChl^-1*' '*d^-1*')')) +
  annotate('text', x = 0, y = .4, 
    label =  bquote('C) µ ('*d^-1*')'),
    hjust = 0) +
  theme_light()

#QN
p4 <- pGrid |>
  filter(alphaChl == alphaChl[10]) |>
  ggplot(aes(x = Topt, y = ESD, z = QN)) +
  geom_contour_filled(show.legend = F) +
  scale_y_log10() +
  xlab(bquote(T[opt] * ' (ºC)')) +
  ylab(bquote('Size (µm)')) + 
  annotate('text', x = 0, y = 120, 
    label = bquote('D) '*Q[N]*' (mol:mol)'),
    hjust = 0) +
  theme_light()

p5 <- pGrid |>
  filter(Topt == Topt[length(Topt)]) |>
  ggplot(aes(x = ESD, y = alphaChl, z = QN)) +
  geom_contour_filled(show.legend = F) +
  scale_x_log10() +
  scale_y_log10() +
  xlab(bquote('Size (µm)')) + 
  ylab(bquote(alpha[Chl] * ' ('* m^2*' '*W^-1*' molC ' * gChl^-1*' '*d^-1*')')) +
  annotate('text', x = .6, y = .4, 
    label = bquote('E) '*Q[N]*' (mol:mol)'),
    hjust = 0) +
  theme_light()

p6 <- pGrid |>
  filter(ESD == ESD[5]) |>
  ggplot(aes(x = Topt, y = alphaChl, z = QN)) +
  geom_contour_filled() +
  scale_y_log10() +
  xlab(bquote(T[opt] * ' (ºC)')) +
  ylab(bquote(alpha[Chl] * ' ('* m^2*' '*W^-1*' molC ' * gChl^-1*' '*d^-1*')')) +
  annotate('text', x = 0, y = .4, 
    label = bquote('F) '*Q[N]*' (mol:mol)'),
    hjust = 0) +
  theme_light()

#C:Chl
p7 <- pGrid |>
  filter(alphaChl == alphaChl[10]) |>
  ggplot(aes(x = Topt, y = ESD, z = theta)) +
  geom_contour_filled(breaks=c(0, 0.002, 0.005, 0.007, 0.01, 0.015, 0.02,0.03, 0.05), 
                      show.legend = F) +
  xlab(bquote(T[opt] * ' (ºC)')) +
  ylab(bquote('Size (µm)')) + 
  scale_y_log10() +
  annotate('text', x = 0, y = 120, 
    label = 'G) Chl:C (g:g)',
    hjust = 0) +
  theme_light()

p8 <- pGrid |>
  filter(Topt == Topt[length(Topt)]) |>
  ggplot(aes(x = ESD, y = alphaChl, z = theta)) +
  geom_contour_filled(breaks=c(0, 0.002, 0.005, 0.007, 0.01, 0.015, 0.02,0.03, 0.05), 
                      show.legend = F) +
  scale_x_log10() +
  scale_y_log10() +
  xlab(bquote('Size (µm)')) + 
  ylab(bquote(alpha[Chl] * ' ('* m^2*' '*W^-1*' molC ' * gChl^-1*' '*d^-1*')')) +
  annotate('text', x = .6, y = .4, 
    label = 'H) Chl:C (g:g)',
    hjust = 0) +
  theme_light()

p9 <- pGrid |>
  filter(ESD == ESD[5]) |>
  ggplot(aes(x = Topt, y = alphaChl, z = theta)) +
  geom_contour_filled(breaks=c(0, 0.002, 0.005, 0.007, 0.01, 0.015, 0.02,0.03, 0.05), 
                      show.legend = T) +
  xlab(bquote(T[opt] * ' (ºC)')) +
  ylab(bquote(alpha[Chl] * ' ('* m^2*' '*W^-1*' molC ' * gChl^-1*' '*d^-1*')')) +
  annotate('text', x = 0, y = .4, 
    label = 'I) Chl:C (g:g)',
    hjust = 0) +
  scale_y_log10() +
  theme_light()

p <- (p1+p2+p3)/(p4+p5+p6)/(p7+p8+p9)
ggsave(paste0('phyto_fitness_landscape_',
              Season,
              '.pdf'), p,
  width = 10, height = 10)

#Using plot3D
#library(plot3D)

#convert the mu
#mu3D <- array(pGrid$mu, dim = c(length(Topt),
#                                length(ESD),
#                                length(alphaChl)))
#jet.colors <- colorRampPalette( c("blue", "green") )
#
## Generate the desired number of colors from this palette
#nbcol <- 100
#color <- jet.colors(nbcol)
#
## Compute the z-value at the facet centres
#z <- mu3D
#zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
#
## Recode facet z-values into color indices
#facetcol <- cut(zfacet, nbcol)
#persp(x = Topt,
#        y = ESD,
#        z = alphaChl,
#        colvar = mu3D)
#Test whether steady state is achieved:
#Time <- seq(1, dtdays*i, length.out = i)
#out <- out[1:i,]
#plot(Time, out[,1],
#     xlab = 'Time',
#     ylab = expression(paste('Growth rate ' *'('* d^-1 * ')')),
#     type = 'l')
#plot(Time, out[,2],
#     xlab = 'Time',
#     ylab = 'N:C',
#     type = 'l')
#plot(Time, out[,3],
#     xlab = 'Time',
#     ylab = 'C:Chl',
#     type = 'l')
