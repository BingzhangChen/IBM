JOHNSON <- function(tC, mumax, Ea, Ed, Topt_){
#Temperature function following Dell et al. PNAS (2011) and Chen & Laws L&O (2017)
#Both tC and Topt_ are in ºC
  kb   = 8.62e-5
  T0   = 273.15
  Tref = 15
  
  if (Ed <=  0) stop("Ed must be greater than zero!")
  
  Eh = Ed+Ea
  x  = TK(tC)
  theta = TK(Topt_)
  b = x - theta
  return(mumax*(Ea/Ed + 1) * exp(Ea*b)/(1. + Ea/Ed*exp(Eh*b)) )
}

TEMPBOL <- function(Ea,tC){
  

  #The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
  # tC: in situ temperature
  # Tr: reference temperature
  
  # boltzman constant constant [ eV /K ]
  kb = 8.62e-5; Tr = 15e0
  
  return(exp(-(Ea/kb)*(1e0/(273.15 + tC)-1e0/(273.15 + Tr))))
}

TK <- function(tC){
#The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
# tC: in situ temperature
# Tr: reference temperature

kb = 8.62e-5; Tr = 15.0
return(-(1./kb)*(1./(273.15 + tC) - 1./(273.15 + Tr)))
}

alloscale <- function(Topt_, mu0p, alpha){
  return(mu0p * exp(TK(Topt_) * alpha) )
}

#------------------------------------------------------------------------------------------------
#Function to estimate photoinhibition following Nikolaou et al. (2016) (J. Theor. Biol.), and
#Han (2001) (J. Theor. Biol.)
#Assuming that acclimation to photoinhibition is at the time-scale of ms.
#------------------------------------------------------------------------------------------------
Ainf <- function(PAR_, alpha_, QN_, QNmin_, QNmax_, theta_){

  Tau   = 5.5e-3 #Turnover time of the electron transport chain [s]
  Beta  = 0.492  #Pre-exponential factor of effective cross-section eq [m2 uE-1 (g Chl)^(1/Kappa) (g C)^(-1/Kappa)]
  Kappa = 0.469  #Exponent of effective cross-section equation [nd]
  Kd    = 5e-6   #Damage constant of a photosynthetic unit [nd]
  WtouE = 4.57   #Constant to convert PAR units from [Wm-2] to [uE m-2 s-1]
  a_ = 2e-5      #The constant a in the equation relating Kr and alphaChl
  b_ = 5e-7      #The constant b in the equation relating Kr and alphaChl
  v_ = -6.64     #The constant v in the equation relating Kr and alphaChl

  #PAR, unit conversion:
  PARWm2 = PAR_ * WtouE    #[W m-2] to [uE m-2 s-1]

  #Convert the unit of alpha from molC/gChl (W m-2)-1 d-1 to  molC/gChl m2/uE
  alpha_new = alpha_ /WtouE/864e2

  #Carbon-specific chlorophyll quota, uit conversion:
  thetaA = theta_ / 12     #[mg Chl mmol C] to [g Chl g C-1]

  #Effective cross-section of the PSU [m2 uE-1] Nikolaou et al. (2016)
  Sigma = Beta * thetaA**Kappa

  #Repair constant of a photosynthetic unit [s-1], following Han et al. (2001):
  
  #Nutrient limitation index
  Lno3 = (QN_ - QNmin_) / (QNmax_ - QNmin_)

  #Kr0 depends on alpha_ using an empirical equation
  Kr0 = a_ * (alpha_new / b_)**v_
  Kr = Kr0 * Lno3

  if (Kr < 1e-10) {
    Y = 0
  }else{
    #Ratio of damage to repair constants [s]:
    K  = Kd / Kr
  
    #Calculate photoinhibition [nd]:
    Y = 1 / (1 + Tau * Sigma * PARWm2 + K * Tau * Sigma**2 * PARWm2**2)
  }
  return(Y) 
}

temp_Topt <- function(tC, mumax0, Topt_){
  #Function of a rate depending on Temperature and optimal temperature (Topt_) modified from Chen Ecol. Mod. (2022)
  #Maximal rate normalized to an optimal temperature of 15 ºC
  #Environmental temperature in ºC
  #Optimal temperature in ºC
  
  Ea0   = 0.98  
  Ed0   = 2.3
  Ei    = 0.22  
  beta  =-0.2  #Exponent for Ea0
  phi   = 0.27 #Exponent for Ed
  
  mumax = alloscale(Topt_, mumax0,  Ei) 
  Ea    = alloscale(Topt_, Ea0,  beta) 
  Ed    = alloscale(Topt_, Ed0,  phi) 
  y     = JOHNSON(tC, mumax, Ea, Ed, Topt_)
  return(y)
}

PHY_C2Vol <- function(p_C){
  
  #p_C : Phytoplankton carbon (pmol C/cell)
  a = -0.69
  b = 0.88
  
  #the parameters of a and b are derived from pg C cell-1 to volume; so we need to
  #convert the carbon unit from pmol C cell-1 to pg C cell-1
  
  y = (12 * p_C/10**a)**(1/b)
  return(y)
}

#Maximal growth rate of phytoplankton as a function of size following Wirtz (2011)
#Ignoring temperature

Pmax_size <- function(ESD, Pmax0 = 5){
  
  a <- .34
  rho_s <- .25
  rho0 <- .5
  a_p <- a*(rho_s/rho0)**(.33333)
  
  #Pmax0: Maximal photosynthesis rate (d-1)
  y <- Pmax0/(1+ a_p * ESD)
  
  return(y)
}

respiration <- function(ESD, r_s = .025){
  #Cell volume when rho_dia = rho_green
  V_s <- 8 #micron^3
  V <- pi/6*ESD^3 
  ESD_s <- (6*V_s/pi)^.333333
  
  b_rho <- 0 #Size scaling of C density
  return(r_s * ESD_s/ESD * (V/V_s)^b_rho)
}

GMK98_Size <- function(Temp, PAR, NO3, C, N, Chl, Cdiv, aI0){
  
  thetaNmax <- 3
  rhoChl_L  <- 0 #the value of rhochl before last sunset
  mu0       <- 5
  QNmin_a   <- 0.07 
  QNmin_b   <- -0.17
  
  #Scaling coefficients of half-saturation constant of nutrient uptake
  KN_a      <- .14
  KN_b      <- .33
  
  #Allometric scaling for Vm
  Vm_a <- .1 
  Vm_b <- .09 
  
  #Minimal NO3 concentration
  NO3_min <- .01
  Ep <- 0.32 #activation energy of phytoplankton growth
  # Temp             :Associated temperarure [degree C]
  # PAR              :Associated PAR [W m-2]
  # NO3              :Associated NO3 concentration [mmol N m-3]
  # C                :Current cellular carbon [pmol C cell-1]
  # N                :Current cellular nitrogen [pmol N cell-1]
  # Chl              :Current cellular Chl [pg C cell-1]
  # Cdiv             :Cellular carbon content threshold for division [pmol cell-1]
  # Topt_            :Optimal temperature [degree C]
  # alphaChl_        :Slope of the P-I curve [Unit the same as aI0]
  # dN               :Changes in the cellular nitrogen content [pmol N cell-1 d-1]
  # dC               :Changes in the cellular carbon content [pmol C cell-1 d-1]
  # dChl             :Changes in the cellular Chl content [pg Chl cell-1 d-1]
  
  Rc   <- 0.1   #Basic C respiration rate [d-1]
  RN   <- 0.1   #Basic N respiration rate [d-1]
  RChl <- 0.1   #Basic Chl respiration rate [d-1]
  zeta <- 3.0   #Cost of biosynthesis [mol C mol N-1]
  
  
  #Factor that governs the down regulation of nutrient uptake as Q approaches Qmax:
  nx <- 1L
  
  #Temperature coefficient
  tf_p <- TEMPBOL(Ep, Temp)  
  
  #Current N:C ratio [mmol N mmol C]:
  QN <- N/C
  
  #Current Chl:C ratio [mg Chl mmol C]
  theta <- Chl/C
  
  #Convert phytoplankton CDiv to Volume:
  Vol <- PHY_C2Vol(Cdiv)
  
  #Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
  KN <- KN_a * Vol**KN_b
  
  #Minimal N:C ratio [mmol N mmol C] following Ward et al. (2012):
  QNmin <- QNmin_a * Vol**QNmin_b
  
  #Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions:
  muT <- mu0 * tf_p
  
  #Apply the size-scaling relationship following Wirtz (2011)
  ESD_ <- (6*Vol/pi)^.3333  
  muT  <- Pmax_size(ESD_, muT)
 
  #Assuming the same temperature dependence of nutrient uptake rate as on photosynthesis rate.
  Vm <- Vm_a*tf_p*Vol**Vm_b
  
  #Maximal N:C ratio
  QNmax <- Vm / muT
  
  dQN <- QNmax - QNmin
  
  #Assume the same temperature dependence of respiration as on photosynthetic rate (long-term adaptation; Barton et al. 2020):
  RcT   <- Rc*tf_p
  RNT   <- RN*tf_p
  RChlT <- RChl*tf_p

  #Nutrient limitation [nd]:
  Lno3 <- (QN - QNmin) / dQN

  #Maximal photosynthesis rate (regulated by QN) [d-1]:
  PCmax <- muT * Lno3
  
  #Light saturation parameter [W m-2 d-1]:
  Ik <- PCmax / aI0 / theta
  
  #Light limitation index [nd]:
  SI <- - PAR / Ik
  
  if (abs(SI) < 1e-10) {
    SI <- 0.
  }else{
    SI <- 1. - exp(SI)
  }
  
  #Photosynthesis rate [d-1]:
  PC <- PCmax * SI

  #Define rhochl [g Chl mol C-1]: fraction of phytoplankton carbon production that is devoted to Chl synthesis.
  #If dark, assume that rhochl equaled the value calculated for the end of the preceding light period.
  if (PAR <= 0) {
    rhochl   <- rhoChl_L
  }else{
    rhochl   <- thetaNmax * PC / aI0 / theta / PAR
    rhoChl_L <- rhochl
  }

  #DIN uptake rate by phytoplankton [mol N mol C-1 d-1]:
  VCN <- Vm * NO3/ (NO3 + KN) * ((QNmax - QN) / dQN)**nx  #Vcref already temperature dependent

  #Changes of cellular carbon [d-1]:
  dC <- C * (PC - zeta * VCN - RcT)

  #Changes of cellular nitrogen [pmol N cell-1 d-1]:
  #RNT has to be zero to avoid continuous decline of N per cell
  dN <- N * (VCN / QN - RNT)

#Changes of cellular Chl [d-1]:
  dChl <- Chl * (rhochl * VCN / theta - RChlT)
  
  return(c(dC, dN, dChl))
}

GMK98_TempSizeLight <- function(Temp, PAR, NO3, Topt_, C, N, Chl, Cdiv, alphachl_){
  
  thetaNmax <- 3
  rhoChl_L <- 0 #the value of rhochl before last sunset
  QNmin_a <- 0.07 
  QNmin_b <- -0.17
  QNmax_a <- 0.167
  QNmax_b <- 0.0 
  mu0     <- 5
  KN_a    <- .14
  KN_b    <- .33
  NO3_min <- .01
  # Temp             :Associated temperature [ºC]
  # PAR              :Associated PAR [W m-2]
  # NO3              :Associated NO3 concentration [mmol N m-3]
  # C                :Current cellular carbon [pmol C cell-1]
  # N                :Current cellular nitrogen [pmol N cell-1]
  # Chl              :Current cellular Chl [pg C cell-1]
  # Cdiv             :Cellular carbon content threshold for division [pmol cell-1]
  # Topt_            :Optimal temperature [degree C]
  # alphaChl_        :Slope of the P-I curve [Unit the same as aI0]
  # dN               :Changes in the cellular nitrogen content [pmol N cell-1 d-1]
  # dC               :Changes in the cellular carbon content [pmol C cell-1 d-1]
  # dChl             :Changes in the cellular Chl content [pg Chl cell-1 d-1]
  
  Rc   = 0.025   #Basic C respiration rate [d-1]
  RN   = 0.025   #Basic N respiration rate [d-1]
  RChl = 0.025   #Basic Chl respiration rate [d-1]
  zeta = 3.0   #Cost of biosynthesis [mol C mol N-1]
  
  
  #Factor that governs the down regulation of nutrient uptake as Q approaches Qmax:
  nx = 1.
  a1 = 0. # Allometric exponent between mumax and alphaChl
  
  if (C <= 0) {
    stop("Non positive carbon!")
  }
  
  #Current N:C ratio [mmol N mmol C]:
  QN = N/C
  
  #Current Chl:C ratio [mg Chl mmol C]
  theta = Chl/C
  
  #Convert phytoplankton CDiv to Volume:
  Vol = PHY_C2Vol(Cdiv)
  ESD_ = (6*Vol/pi)^.3333  
    
  #Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
  KN = KN_a * Vol**KN_b
  
  #Minimal N:C ratio [mmol N mmol C] following Ward et al. (2012):
  QNmin = QNmin_a * Vol**QNmin_b
  
  #Maximal N:C ratio [mmol N mmol C] following Ward et al. (2012):
  QNmax = QNmax_a * Vol**QNmax_b
  
  #Constrain QN between QNmin and QNmax due to numerical issues
  QN = max(min(QN, QNmax), QNmin)
  
  #(Qmax - Qmin) [mmol N mmol C]:
  dQN = QNmax - QNmin
  
  #Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions:
    #mu0 should be a function of alphaChl
  muT = mu0 * exp(a1 * (alphachl_ - .1)) #0.1 is the average alphaChl value
  
  #Temperature dependent maximal growth rate at 1 um
  muT = temp_Topt(Temp, muT, Topt_)
  
  #Apply the size-scaling relationship following Wirtz (2011)
  muT = Pmax_size(ESD_, muT)
    
  #Assuming the same temperature dependence of nutrient uptake rate as on photosynthesis rate.
  Vcref = muT * QNmax
  
  #Assume the same temperature dependence of respiration as on photosynthetic rate (long-term adaptation; Barton et al. 2020):
  RcT <- temp_Topt(Temp, Rc,   Topt_)
  RcT <- respiration(ESD_, RcT)
  
  RNT = temp_Topt(Temp, RN,   Topt_)
  RNT = respiration(ESD_, RNT)
  
  RChlT = temp_Topt(Temp, RChl, Topt_)
  RChlT = respiration(ESD_, RChlT)

  #Nutrient limitation [nd]:
  Lno3 = (QN - QNmin) / dQN

  if (Lno3 <= 0) {
    PC = 0
  }else{
    #Maximal photosynthesis rate (regulated by QN) [d-1]:
    PCmax = muT * Lno3
  
    #Light saturation parameter [W m-2 d-1]:
    Ik = PCmax / alphachl_ / theta
  
    #Calculate the fraction of open PSU [nd]:
    if (PAR > 0.) { #Photoinhibition
      #Photoinhibition, following Nikolau et al. (2016)
      A = Ainf(PAR, alphachl_, QN, QNmin, QNmax, theta)
    }else{
      A = 1
    }
  
    #Light limitation index [nd]:
    SI <- - A * PAR / Ik
  
    if (abs(SI) < 1e-10) {
      SI <- 0. #zero light
    }else{
      SI <- 1. - exp(SI) 
    }
  
    #Photosynthesis rate [d-1]:
    PC <- PCmax * SI
  }

  #Define rhochl [g Chl mol C-1]: fraction of phytoplankton carbon production that is devoted to Chl synthesis.
  #If dark, assume that rhochl equaled the value calculated for the end of the preceding light period.
  if (PAR <= 0) {
    rhochl   <- rhoChl_L
  }else{
    rhochl   <- thetaNmax * PC / alphachl_ / theta / PAR
    rhoChl_L <- rhochl
  }

  #DIN uptake rate by phytoplankton [mol N mol C-1 d-1]:
  VCN = Vcref * (NO3 - NO3_min)/ (NO3 + KN) * ((QNmax - QN) / dQN)**nx  #Vcref already temperature dependent
  VCN = max(VCN, 0)

  #Changes of cellular carbon [d-1]:
  dC = C * (PC - zeta * VCN - RcT)

  #Changes of cellular nitrogen [pmol N cell-1 d-1]:
  #RNT has to be zero to avoid continuous decline of N per cell
  dN = N * (VCN / QN - RNT)

#Changes of cellular Chl [d-1]:
  dChl = Chl * (rhochl * VCN / theta - RChlT)
  
  return(c(dC, dN, dChl))
}