SUBROUTINE GMK98_Ind_TempSize(Temp, PAR, NO3, Topt_, C, N, Chl, Cdiv, dC, dN, dChl)
USE Trait_functions, only : temp_Topt, PHY_C2Vol
USE params,          only : Ep, aI0, thetaNmax, mu0, rhoChl_L, nx

implicit none

! Declaration of variables:
real, intent(in)  :: Temp          ! Associated temperarure [degree C]
real, intent(in)  :: PAR           ! Associated PAR [W m-2]
real, intent(in)  :: NO3           ! Associated NO3 concentration [mmol N m-3]

real, intent(in)  :: C             ! Current cellular carbon [mmol C cell-1]
real, intent(in)  :: N             ! Current cellular nitrogen [mmol N cell-1]
real, intent(in)  :: Chl           ! Current cellular Chl [mg C cell-1]
real, intent(in)  :: Topt_         ! Optimal temperature [degree C]
real, intent(in)  :: Cdiv          ! Cellular carbon content threshold for division [pmol]

real              :: KN    = 0.50  ! Nitrate half-saturation constant of phyto growth [uM]
real              :: QNmin = 0.05  ! Minimal N:C ratio [mmol N mmol C]
real              :: QNmax = 0.18  ! Maximal N:C ratio [mmol N mmol C]
real              :: dQN   = 0.13  ! (Qmax - Qmin) [mmol N mmol C]
real              :: QN    = 0.    ! Current N:C ratio [mmol N mmol C]
real              :: theta = 0.    ! Current Chl:C ratio [mg Chl mmol C]

! Changes in the cellular nitrogen content [pmol N cell-1]:
real, intent(out) :: dN

! Changes in the cellular carbon content [pmol N cell-1]:
real, intent(out) :: dC

! Changes in the cellular Chl content [pg Chl cell-1]:
real, intent(out) :: dChl

! Maximal growth rate as a function of temperature under resource (nutrient and light)
! replete conditions [uM]:
real :: muT = 0.

! Maximal specific nitrogen uptake rate as a function of temperature under resource
! (nutrient and light) replete conditions [mol N mol C-1 d-1]:
! Vcref = Qmax * muT
real :: Vcref = 0.

! DIN uptake rate by phytoplankton [mol N mol C-1 d-1]:
real :: VCN = 0.

! Indices for nutrient, light and temperature limitation
real :: Lno3 = 0.   ! Nutrient limitation [pmol N m-3]
real :: SI   = 0.   ! Light limitation index [fpar]

real :: PCmax  = 0. ! Maximal photosynthesis rate (regulated by QN) [d-1]
real :: PC     = 0. ! Carbon specific rate of photosynthesis [d-1]
real :: rhoChl = 0. ! Phyto C production devoted to Chl synthesis [mg Chl mmol C-1]
real :: Ik     = 0. ! Saturation parameter for the PI curve [umol photons m-2 s-1]

real, parameter   :: RC   = 0.025d0  ! Basic respiration rate [d-1]
real, parameter   :: RN   = 0.025d0  ! Basic respiration rate [d-1]
real, parameter   :: RChl = 0.025d0  ! Basic respiration rate [d-1]
real, parameter   :: zeta = 3.0d0    ! Cost of biosynthesis [mol C mol N-1]

!Kn is an allometric function of Vol
real, parameter   :: KN_a = 10**(-0.84)  ! Normalization constant for KN
real, parameter   :: KN_b = 0.33d0  ! Allometric exponent for KN

!QNmin and QNmax are allometric functions of CDiv
real, parameter   :: QNmin_a = 0.12d0  ! Normalization constant for QNmin (pmol N per cell)
real, parameter   :: QNmin_b = 0.95d0  ! Allometric exponent for QNmin
real, parameter   :: QNmax_a = 0.2d0   ! Normalization constant for QNmax (pmol N per cell)
real, parameter   :: QNmax_b = 1d0  ! Allometric exponent for QNmax

real :: RcT   = 0. ! Temperature dependent respiration rate
real :: RNT   = 0. ! Temperature dependent N-based respiration rate
real :: RChlT = 0. ! Temperature dependent Chl-based respiration rate
real :: Vol = 0d0 !Cell volume of phytoplankton
! End of declaration
   
! Check input
if (NO3 .le. 0.d0) stop "Negative Nitrate concentration!"
if (PAR .lt. 0.d0) stop "Negative PAR!"

if (C .le. 0d0) then
   dN   = 0.d0
   dC   = 0.d0
   dChl = 0.d0
   return
endif

QN      = N / C
theta   = Chl / C
   
!Convert phytoplankton CDiv to Volume
Vol = PHY_C2Vol(CDiv)

! Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
KN = KN_a * Vol**KN_b

! Minimal N:C ratio [mmol N mmol C]:
QNmin = QNmin_a * Cdiv**(QNmin_b-1d0)
   
! Maximal N:C ratio [mmol N mmol C]:
QNmax = QNmax_a * Cdiv**(QNmax_b-1d0)
dQN   = QNmax - QNmin

! Temperature coefficient (what value of mu0 should be?):
muT   = temp_Topt(Temp, mu0, Topt_)
Vcref = muT * QNmax

!Assume the same temperature dependence of respiration as on photosynthetic 
!rate (long-term adaptation; Barton et al. 2020):
RcT   = temp_Topt(Temp, RC,   Topt_)
RNT   = temp_Topt(Temp, RN,   Topt_)
RChlT = temp_Topt(Temp, RChl, Topt_)

! Nutrient limitation [pmol N m-3]:
Lno3 = (QN - QNmin) / dQN

! Maximal photosynthesis rate (regulated by QN) [d-1]:
PCmax = muT * Lno3
Ik    = PCmax / aI0 / theta

! Light limitation index [fpar]:
SI = -PAR/Ik

if (abs(SI) < 1d-10) then
   SI = 0d0
else
   SI = 1.d0 - exp(SI)
endif

PC = PCmax * SI

! Define rhoChl [g Chl mol C-1]:
! If dark, assume that rhoChl equaled the value calculated for the end of the 
! preceding light period .
if (PAR <= 0d0) then
   rhoChl   = rhoChl_L
else
   rhoChl   = thetaNmax * PC / aI0 / theta / PAR
   rhoChl_L = rhoChl
endif

! DIN uptake rate by phytoplankton [mol N mol C-1 d-1]:
VCN = Vcref * NO3 / (NO3 + KN) * ((QNmax - QN) / dQN)**nx ! Vcref already temperature dependent

! Changes of cellular carbon [pmol C cell-1 d-1]:
dC  = C * (PC - zeta * VCN - RcT)

! Changes of cellular nitrogen [pmol N cell-1 d-1]:
dN  = N * (VCN / QN - RNT)

! Changes of cellular Chl [pg Chl cell-1 d-1]:
dChl= Chl * (rhoChl * VCN / theta - RChlT)

return
END subroutine GMK98_Ind_TempSize
