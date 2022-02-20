Module Trait_functions
!This module provides several functions calculating phytoplankton physiological rates as a function of environmental conditions (e.g., temperature) and traits
implicit none

private

public :: TEMPBOL, temp_Topt

CONTAINS

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