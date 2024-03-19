!Adding optimal temperature into the Geider model
SUBROUTINE GMK98_Ind_Temp(Temp, PAR, NO3, Topt_, C, N, Chl, dC, dN, dChl)
USE Trait_functions, only : temp_Topt
USE params,          only : Ep, aI0, thetaNmax, mu0, KN, rhoChl_L, nx
implicit none

real, intent(in)  :: Temp, PAR, NO3
real, intent(in)  :: C    !Current cellular carbon
real, intent(in)  :: N    !Current cellular nitrogen
real, intent(in)  :: Chl  !Current cellular Chl
real, intent(in)  :: Topt_  !Optimal temperature in C

! Minimal and maximal N:C ratio
real, parameter   :: QNmin = 0.05, QNmax = 0.18, dQN = 0.13  

! Current N:C ratio
real    :: QN    = 0.

! Current Chl:C ratio
real    :: theta = 0.

! Changes in the cellular nitrogen content (pmol/cell)
real, intent(out) :: dN

! Changes in the cellular carbon content
real, intent(out) :: dC  

! Changes in the cellular Chl content
real, intent(out) :: dChl  

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions
real :: muT = 0.

! Maximal specific nitrogen uptake rate as a function of temperature under resource (nutrient and light) replete conditions (molN molC-1 d-1; Vcref = Qmax * muT)
real :: Vcref    = 0.

! DIN uptake rate by phytoplankton (molN/molC/d)
real :: VCN = 0.

! Indices for nutrient, light and temperature limitation
real :: Lno3, SI

real :: PCmax, PC, rhochl, Ik

real, PARAMETER   :: Rc   = 0.025  !Basic respiration rate at the reference temperature 15 C (d-1)
real, PARAMETER   :: RN   = 0.025  !Basic N-based respiration rate at the reference temperature 15 C (d-1)
real, PARAMETER   :: RChl = 0.025  !Basic Chl-based respiration rate at the reference temperature 15 C (d-1)
real, PARAMETER   :: zeta = 3.0    !cost of biosynthesis (molC molN-1)

real :: RcT  = 0. !Temperature dependent respiration rate
real :: RNT = 0. !Temperature dependent N-based respiration rate
real :: RChlT = 0. !Temperature dependent Chl-based respiration rate
!End of declaration

!Check input
if (NO3 .le. 0.d0) stop "Negative Nitrate concentration!"
if (PAR .lt. 0.d0) stop "Negative PAR!"

if (C .le. 0d0) then
   dN=0.d0
   dC=0.d0
   dChl=0.d0
   return
endif

QN      = N/C
theta   = Chl/C

!Temperature coefficient (what value of mu0 should be?)
muT    = temp_Topt(Temp, mu0, Topt_)
Vcref  = muT * QNmax  !Assuming the same temperature dependence of nutrient uptake rate as on photosynthesis rate. QNmax/QNmin may be also a function of temperature which needs being further investigated. 

!Assume the same temperature dependence of respiration as on photosynthetic rate (long-term adaptation; Barton et al. 2020)
RcT     = temp_Topt(Temp, Rc, Topt_)
RNT    = temp_Topt(Temp, RN, Topt_)
RChlT= temp_Topt(Temp, RChl, Topt_)

! N limitation
Lno3 = (QN-QNmin)/dQN

!Maximal photosynthesis rate (regulated by QN)
PCmax = muT*Lno3

Ik = PCmax/aI0/theta

!The light limitation index (fpar)
SI = 1.d0 - exp(-par/Ik)

PC = PCmax * SI

! define rhochl (gChl/molC; the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
!If dark, assume that rhochl equaled the value calculated for the end of the preceding light period 
if (PAR <= 0d0) then
   rhochl   = rhoChl_L
else
   rhochl   = thetaNmax*PC/aI0/theta/PAR
   rhoChl_L = rhochl
endif

! DIN uptake rate by phytoplankton (molN/molC/d)
VCN = Vcref * NO3/(NO3 + KN)* ((QNmax-QN)/dQN)**nx  !Vcref already temperature dependent

! Changes of cellular carbon
dC  = C*(PC - zeta*VCN - RcT)

! Changes of cellular nitrogen (pmol N /d)
dN  = N*(VCN/QN - RNT)

! Changes of cellular Chl
dChl= Chl*(rhochl*VCN/theta - RChlT)

return
END SUBROUTINE GMK98_Ind_Temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GMK98_Ind_Size(Temp, PAR, NO3, C, N, Chl, Cdiv, dC, dN, dChl)

USE Trait_functions, only : PHY_C2Vol
USE Trait_functions, only : TEMPBOL
USE params,          only : Ep, aI0, thetaNmax, mu0, rhoChl_L
USE params,          only : QNmin_a, QNmin_b, nx
USE params,          only : QNmax_a, QNmax_b, KN_a, KN_b
USE state_variables, only : NO3_min

implicit none

! Declaration of variables:
real, intent(in)  :: Temp          ! Associated temperarure [degree C]
real, intent(in)  :: PAR           ! Associated PAR [W m-2]
real, intent(in)  :: NO3           ! Associated NO3 concentration [mmol N m-3]
real, intent(in)  :: C             ! Current cellular carbon [mmol C cell-1]
real, intent(in)  :: N             ! Current cellular nitrogen [mmol N cell-1]
real, intent(in)  :: Chl           ! Current cellular Chl [mg C cell-1]
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
real :: tf_p = 0.   ! Light limitation index [fpar]

real :: PCmax  = 0. ! Maximal photosynthesis rate (regulated by QN) [d-1]
real :: PC     = 0. ! Carbon specific rate of photosynthesis [d-1]
real :: rhoChl = 0. ! Phyto C production devoted to Chl synthesis [mg Chl mmol C-1]
real :: Ik     = 0. ! Saturation parameter for the PI curve [umol photons m-2 s-1]

real, parameter   :: RC   = 0.0d0  ! Basic respiration rate [d-1]
real, parameter   :: RN   = 0.0d0  ! Basic respiration rate [d-1]
real, parameter   :: RChl = 0.0d0  ! Basic respiration rate [d-1]
real, parameter   :: zeta = 3.0d0    ! Cost of biosynthesis [mol C mol N-1]

real, parameter   :: b0 = 4.3d0       ! Intercept between mumax and size 
real, parameter   :: b1 = -0.65d0     ! Allometric exponent between mumax and size
real, parameter   :: b2 = -0.08d0     ! Allometric exponent between mumax and size

!temperature correction
real, parameter   :: Tempcorr = exp(0.32/(8.62d-5 * 288.15))
real              :: CDiv1 = 0.       ! CDiv value with ug C per cell
!
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

tf_p = TEMPBOL(Ep, Temp)  

QN      = N / C
theta   = Chl / C
   
!Convert phytoplankton CDiv to Volume
Vol = PHY_C2Vol(CDiv)

! Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
KN = KN_a * Vol**KN_b

! Minimal N:C ratio [mmol N mmol C]:
QNmin = QNmin_a * Vol**QNmin_b
   
! Maximal N:C ratio [mmol N mmol C]:
QNmax = QNmax_a * Vol**QNmax_b

QN = max(min(QN, QNmax), QNmin)

dQN= QNmax - QNmin

! Maximal growth rate increases with temperature 
muT   = tf_P * mu0

!Apply unimodal size-scaling relationship
!First, convert the unit of CDiv to ugC cell-1
CDiv1  = CDiv * 12.d0/1d6

!Applying the size scaling following Chen and Liu (2011) Fig. 1B
muT = muT*10**(b0 + b1 * log10(Cdiv1) + b2 * log10(Cdiv1)**2)/Tempcorr

Vcref = muT * QNmax

!Assume the same temperature dependence of respiration as on photosynthetic 
!rate (long-term adaptation; Barton et al. 2020):
RcT   = RC*tf_p
RNT   = RN*tf_p
RChlT = RChl*tf_p

! Nutrient limitation [pmol N m-3]:
Lno3 = (QN - QNmin) / dQN

! Maximal photosynthesis rate (regulated by QN) [d-1]:
PCmax = muT * Lno3
Ik    = PCmax / aI0 / theta

! Light limitation index [fpar]:
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
END subroutine GMK98_Ind_Size
!