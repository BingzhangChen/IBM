Module Trait_functions
use params
!This module provides several functions calculating phytoplankton physiological rates as a function of environmental conditions (e.g., temperature) and traits
implicit none

private

public :: TEMPBOL, temp_Topt, palatability, PHY_C2Vol, PHY_ESD2C, Ainf


CONTAINS

!Function converting phytoplankton ESD (micron) to Carbon (unit: pmol C per cell) with the parameters obtained from Maranon et al. (2013)
real function PHY_ESD2C(ESD) result(y)
implicit none
real, intent(in)  :: ESD !Phytoplankton ESD (micron)

real, parameter :: a = -0.69
real, parameter :: b = 0.88

real :: Vol = 0d0

!Calculate volume from ESD
Vol = pi/6d0 * ESD**3

!the parameters of a and b are derived from pg C cell-1 to volume; so we need to
!convert the carbon unit from pmol C cell-1 to pg C cell-1

!Calculate carbon (pmol C per cell) from volume
y = 10.d0**a * Vol**b / 12.d0

return
end function PHY_ESD2C

!Function converting phytoplankton carbon to volume (unit: micron^3) with the parameters obtained from Maranon et al. (2013)
pure real function PHY_C2Vol(p_C) result(y)
implicit none
real, intent(in)  :: p_C !Phytoplankton carbon (pmol C/cell)

real, parameter :: a = -0.69
real, parameter :: b = 0.88

!the parameters of a and b are derived from pg C cell-1 to volume; so we need to
!convert the carbon unit from pmol C cell-1 to pg C cell-1

y = (12d0 * p_C/10.d0**a)**(1d0/b)
return
end function PHY_C2Vol

!Function calculating the prey palatability based on Ward et al. L&O 2012 (Eq. A21) and Banas Ecol. Mod. (2011) Table 2
real function palatability(Vpred, Vprey, SDpref_Z) result(y)
implicit none

!Predator Volume
real, intent(in)  :: Vpred

!Prey Volume
real, intent(in)  :: Vprey

! Standard deviation of log zooplankton feeding preference
real, intent(in)  :: SDpref_Z

!The actual predator:prey volume ratio
real  :: R_real = 0d0

!Optimal prey Volume of the predator
real :: Vprey_opt = 1d3 

real :: Xpred, Xprey !ESD of predator and prey
real :: Xprey_opt = 0. !Optimal prey ESD

!Optimal predator:prey volume ratio
real :: R_opt = 1d3 

real :: cff = 0d0

!First calculate prey and predator ESD (micron) from volume
Xpred = (6d0*Vpred/pi)**(.33333333333333)
Xprey = (6d0*Vprey/pi)**(.33333333333333)

!Then calculate optimal prey ESD (micron)
Xprey_opt = 0.65 * Xpred**0.56

!Convert prey ESD to volume
Vprey_opt = pi/6d0 * Xprey_opt**3

R_opt  = Vpred/Vprey_opt
R_real = Vpred/Vprey

cff = log(R_real/R_opt)
cff = cff**2/(2.d0 * SDpref_Z**2)

!To avoid underflow
if (cff > 5d2) then
   y = 0d0
else
   y = exp(-cff)
endif
return
end function palatability

pure real function TEMPBOL(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: Ea, tC

! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15D0

TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
return 
end function TEMPBOL

REAL function temp_Topt(tC, mumax0, Topt_) result(y)
!Function of a rate depending on Temperature and optimal temperature (Topt_) modified from Chen Ecol. Mod. (2022)
IMPLICIT NONE
real, intent(in) :: mumax0    !Maximal rate normalized to an optimal temperature of 15 ºC
real, intent(in) :: tC         !Environmental temperature in ºC
real, intent(in) :: Topt_   !Optimal temperature in ºC

real, parameter   :: Ea0   = 0.98  
real, parameter   :: Ed0   = 2.3
real, parameter   :: Ei      = 0.22  
real, parameter   :: beta  =-0.2  !Exponent for Ea0
real, parameter   :: phi    = 0.27  !Exponent for Ed
!real, parameter   :: mumax0 = 0.59  !Normalized growth rate for mumax (d-1)

real :: Ed, Ea, mumax

mumax = alloscale(Topt_, mumax0,  Ei) 
Ea    = alloscale(Topt_, Ea0,  beta) 
Ed    = alloscale(Topt_, Ed0,  phi) 
y     = JOHNSON(tC, mumax, Ea, Ed, Topt_)

END function temp_Topt

REAL FUNCTION JOHNSON(tC, mumax, Ea, Ed, Topt_) RESULT(y)
!Temperature function following Dell et al. PNAS (2011) and Chen & Laws L&O (2017)
IMPLICIT NONE
!Both tC and Topt_ are in ºC
real,   intent(in)     :: tC, mumax, Ea, Ed, Topt_
real,   parameter   :: kb   = 8.62D-5
real,   parameter   :: T0   = 273.15D0
real,   parameter   :: Tref = 15D0
real                         :: Eh, x, theta, b

if (Ed .le. 0d0) stop "Ed must be greater than zero!"
Eh = Ed+Ea
x    = TK(TC)
theta = TK(Topt_)
b = x - theta
y = mumax*(Ea/Ed + 1.d0) * exp(Ea*b)/(1.D0+Ea/ED*exp(Eh*b))   
return
END FUNCTION JOHNSON

PURE REAL FUNCTION TK(TC)
IMPLICIT NONE
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
REAL, INTENT (IN) :: TC
! boltzman constant constant [ eV /K ]
REAL, PARAMETER   :: kb = 8.62d-5, Tr = 15.0

TK = -(1./kb)*(1./(273.15 + tC) - 1./(273.15 + Tr))
return 
END FUNCTION TK

PURE REAL FUNCTION alloscale(Topt_, mu0p, alpha)
IMPLICIT NONE
real, intent(in) :: Topt_     !Topt in ºC
real, intent(in) :: mu0p  !Normalized growth rate
real, intent(in) :: alpha    !Exponent of thermal traits normalized to z
alloscale =  mu0p * exp(TK(Topt_) * alpha) 
END FUNCTION alloscale

!------------------------------------------------------------------------------------------------
!Function to estimate photoinhibition following Nikolaou et al. (2016) (J. Theor. Biol.), and
!Han (2001) (J. Theor. Biol.)
!Assuming that acclimation to photoinhibition is at the time-scale of ms.
!------------------------------------------------------------------------------------------------
PURE REAL FUNCTION Ainf(PAR_, alpha_, QN_, QNmin_, QNmax_, theta_)

implicit none

!Declaration of variables:
real, intent(in) :: PAR_           !Irradiance [W m-2]
real, intent(in) :: alpha_         !Slope of the P-I curve [Unit: molC/gChl m2/uE]
real, intent(in) :: QN_            !N:C ratio of the phyto. super-individual [mol N mol C-1]
real, intent(in) :: QNmin_         !Minimal N:C ratio
real, intent(in) :: QNmax_         !Maximal N:C ratio
real, intent(in) :: theta_         !Chl:C ratio [mg Chl mmol C]

real, parameter  :: Tau   = 5.5d-3 !Turnover time of the electron transport chain [s]
real, parameter  :: Beta  = 0.492  !Pre-exponential factor of effective cross-section eq [m2 uE-1 (g Chl)^(1/Kappa) (g C)^(-1/Kappa)]
real, parameter  :: Kappa = 0.469  !Exponent of effective cross-section equation [nd]
real, parameter  :: Kd    = 5d-6   !Damage constant of a photosynthetic unit [nd]

real, parameter  :: WtouE = 4.57   !Constant to convert PAR units from [Wm-2] to [uE m-2 s-1]
real, parameter  :: a_ = 2d-5      !The constant a in the equation relating Kr and alphaChl
real, parameter  :: b_ = 5d-7      !The constant b in the equation relating Kr and alphaChl
real, parameter  :: v_ = -6.64     !The constant v in the equation relating Kr and alphaChl
real             :: Kr0            !Repair constant of a photosynthetic unit [s-1] under nutrient saturated conditions which depends on alpha to impose a tradeoff
real             :: Kr             !Nutrient dependent Repair constant of a photosynthetic unit [s-1]
real             :: K              !Ratio of damage to repair constants [s]
real             :: Sigma          !Effective cross-section of the PSU [m2 uE-1]

real             :: thetaA         !Chl:C ratio (g Chl g C-1) to be consistent with Han (2001)
real             :: PARWm2
real             :: Lno3           !Nutrient limitation index
real             :: alpha_new      !alphaChl with the correct unit
!End of declaration

!PAR, unit conversion:
PARWm2 = PAR_ * WtouE    ![W m-2] to [uE m-2]

!Convert the unit of alpha from molC/gChl (W m-2)-1 d-1 to  molC/gChl m2/uE
alpha_new = alpha_ /WtouE/864d2

!Carbon-specific chlorophyll quota, uit conversion:
thetaA = theta_ / 12d0     ![mg Chl mmol C] to [g Chl g C-1]

!Effective cross-section of the PSU [m2 uE-1] Nikolaou et al. (2016)
Sigma = Beta * thetaA**Kappa

!Repair constant of a photosynthetic unit [s-1], following Han et al. (2001):

!Nutrient limitation index
Lno3 = (QN_ - QNmin_) / (QNmax_ - QNmin_)

!Kr0 depends on alpha_ using an empirical equation
Kr0 = a_ * (alpha_new / b_)**v_

Kr = Kr0 * Lno3

if (Kr < 1d-10) then
  Ainf = 0.d0
else
  !Ratio of damage to repair constants [s]:
  K  = Kd / Kr
  
  !Calculate photoinhibition [nd]:
  Ainf = 1d0 / (1d0 + Tau * Sigma * PARWm2 + K * Tau * Sigma**2 * PARWm2**2)
endif

return 
END FUNCTION Ainf
!------------------------------------------------------------------------------------------------

END MODULE
