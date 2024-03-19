SUBROUTINE GMK98_Ind(Temp, PAR, NO3, C, N, Chl, dC, dN, dChl)
USE Trait_functions, only : TEMPBOL
USE params,          only : Ep, aI0, thetaNmax, mu0, KN, rhoChl_L
implicit none

real, intent(in)  :: Temp, PAR, NO3
real, intent(in)  :: C    !Current cellular carbon
real, intent(in)  :: N    !Current cellular nitrogen
real, intent(in)  :: Chl  !Current cellular Chl

! Minimal and maximal N:C ratio
real, parameter   :: QNmin = 0.05, QNmax = 0.18, dQN = 0.13  

! Current N:C ratio
real    :: QN    = 0.

! Current Chl:C ratio
real    :: theta = 0.

! Maximal specific nitrogen uptake rate (molN molC-1 d-1; Vcref = Qmax * Pcref)
real    :: Vcref    = 0.

! Changes in the cellular nitrogen content (pmol/cell)
real, intent(out) :: dN

! Changes in the cellular carbon content
real, intent(out) :: dC  

! Changes in the cellular Chl content
real, intent(out) :: dChl  

! DIN uptake rate by phytoplankton (molN/molC/d)
real :: VCN = 0.

! Indices for nutrient, light and temperature limitation
real :: Lno3, SI, tf_p

real :: PCmax, PC, rhochl, Ik

real, PARAMETER   :: Rc   = 0.025  !Basic respiration rate (d-1)
real, PARAMETER   :: RN   = 0.025  !Basic respiration rate (d-1)
real, PARAMETER   :: RChl = 0.025  !Basic respiration rate (d-1)
real, PARAMETER   :: zeta = 3.0    !cost of biosynthesis (molC molN-1)
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

QN       = N/C

!Force QN to be between QNmin and QNmax
QN = max(min(QNmax, QN), QNmin)

theta    = Chl/C
Vcref    = mu0 * QNmax

!Temperature coefficient
tf_p = TEMPBOL(Ep, Temp)  

! N limitation
Lno3 = (QN-QNmin)/dQN

!Maximal photosynthesis rate (regulated by QN)
PCmax = mu0 *tf_P*Lno3

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
VCN = Vcref * NO3/(NO3 + KN)* (QNmax-QN)/dQN*tf_p

! Changes of cellular carbon
dC  = C*(PC - zeta*VCN - Rc*tf_p)

! Changes of cellular nitrogen (pmol N /d)
dN  = N*(VCN/QN - RN*tf_p)

! Changes of cellular Chl
dChl= Chl*(rhochl*VCN/theta - RChl*tf_p)

return
end subroutine GMK98_Ind
