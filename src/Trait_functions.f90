Module Trait_functions
!This module provides several functions calculating phytoplankton physiological rates as a function of environmental conditions (e.g., temperature) and traits
implicit none

private

public :: TEMPBOL, temp_Topt, palatability, PHY_C2Vol, PHY_ESD2C

real, parameter :: pi= 3.1415926535

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
end function 

!Function calculating the prey palatability based on Ward et al. L&O 2012 (Eq. A21) and Banas Ecol. Mod. (2011) Table 2
real function palatability(Vpred, Vprey) result(y)
implicit none

!Predator Volume
real, intent(in)  :: Vpred

!Prey Volume
real, intent(in)  :: Vprey

!The actual predator:prey volume ratio
real  :: R_real = 0d0

! Standard deviation of log zooplankton feeding preference
real, parameter  :: SDpref_Z = 0.1

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
y       = JOHNSON(tC, mumax, Ea, Ed, Topt_)

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

END MODULE
