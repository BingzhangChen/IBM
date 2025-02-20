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