SUBROUTINE GMK98_Ind_TempSizeLight(Temp, PAR, NO3, Topt_, C, N, Chl, Cdiv, alphaChl_, dC, dN, dChl)
USE Trait_functions, only : temp_Topt, PHY_C2Vol, Ainf, Pmax_size, respiration
USE params,          only : thetaNmax, mu0, rhoChl_L, QNmin_a, QNmin_b
USE params,          only : QNmax_a, QNmax_b, KN_a, KN_b, nx, pi
USE state_variables, only : NO3_min

implicit none

!Declaration on variables:
real, intent(in)  :: Temp             !Associated temperarure [degree C]
real, intent(in)  :: PAR              !Associated PAR [W m-2]
real, intent(in)  :: NO3              !Associated NO3 concentration [mmol N m-3]
real, intent(in)  :: C                !Current cellular carbon [pmol C cell-1]
real, intent(in)  :: N                !Current cellular nitrogen [pmol N cell-1]
real, intent(in)  :: Chl              !Current cellular Chl [pg C cell-1]
real, intent(in)  :: Cdiv             !Cellular carbon content threshold for division [pmol cell-1]

real, intent(in)  :: Topt_            !Optimal temperature [degree C]
real, intent(in)  :: alphaChl_        !Slope of the P-I curve [Unit the same as aI0]
real, intent(out) :: dN               !Changes in the cellular nitrogen content [pmol N cell-1 d-1]
real, intent(out) :: dC               !Changes in the cellular carbon content [pmol C cell-1 d-1]
real, intent(out) :: dChl             !Changes in the cellular Chl content [pg Chl cell-1 d-1]
real              :: Vol    = 0d0     !Cell volume of phytoplankton [um3]
real              :: ESD_   = 0d0     !ESD of phytoplankton [um]
real              :: QNmin  = 0.05    !Minimal N:C ratio [mmol N mmol C]
real              :: QNmax  = 0.18    !Maximal N:C ratio [mmol N mmol C]
real              :: dQN    = 0.13    !(Qmax - Qmin) [mmol N mmol C]
real              :: QN     = 0.      !Current N:C ratio [mmol N mmol C]
real              :: theta  = 0.      !Current Chl:C ratio [mg Chl mmol C]
real              :: VCN    = 0.      !DIN uptake rate by phytoplankton [mol N mol C-1 d-1]
real              :: Lno3   = 0.      !Nutrient limitation [pmol N m-3]
real              :: SI     = 0.      !Light limitation index [fpar]
real              :: PCmax  = 0.      !Maximal photosynthesis rate (regulated by QN) [d-1]
real              :: PC     = 0.      !Carbon specific rate of photosynthesis [d-1]
real              :: rhoChl = 0.      !Phyto C production devoted to Chl synthesis [mg Chl mmol C-1]
real              :: Ik     = 0.      !Saturation parameter for the PI curve [W m-2 s-1]
real              :: A      = 0.      !Photoinhibition, following Nikolau et al. (2016)

real, parameter   :: RC     = 0.1d0   !Basic C respiration rate [d-1]
real, parameter   :: RN     = 0.1d0   !Basic N respiration rate [d-1]
real, parameter   :: RChl   = 0.1d0   !Basic Chl respiration rate [d-1]
real              :: RcT    = 0.d0    !Temperature dependent C-based respiration rate [d-1]
real              :: RNT    = 0.d0    !Temperature dependent N-based respiration rate [d-1]
real              :: RChlT  = 0.d0    !Temperature dependent Chl-based respiration rate [d-1]

real, parameter   :: zeta   = 3.0d0   !Cost of biosynthesis [mol C mol N-1]

! Maximal specific N uptake as a function of temp. under resource (nutrient and light) replete conditions [mol N mol C-1 d-1]:
real              :: Vcref  = 0.

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions [uM]:
real              :: muT    = 0.

!Kn is an allometric function of Vol (Cdiv) (Edwards et al. 2012) [uM]:
real              :: KN     = 0.      !Half-saturation constant [uM]
real, parameter   :: a1 = 0.d0        !Allometric exponent between mumax and alphaChl
real              :: CDiv1 = 0.       !CDiv value with ug C per cell
!End of declaration

if (C .le. 0d0) then
   dN   = 0.d0
   dC   = 0.d0
   dChl = 0.d0
   return
endif

!Current N:C ratio [mmol N mmol C]:
QN = N/C

!Current Chl:C ratio [mg Chl mmol C]
theta = Chl/C

!Convert phytoplankton CDiv to Volume:
Vol = PHY_C2Vol(CDiv)

!Convert Volume to ESD:
ESD_ = (6.d0*Vol/pi)**0.3333333  

!Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
KN = KN_a * Vol**KN_b

!Minimal N:C ratio [mmol N mmol C] following Ward et al. (2012):
QNmin = QNmin_a * Vol**QNmin_b

!Maximal N:C ratio [mmol N mmol C] following Maranon et al. (2013)
QNmax = QNmax_a * Vol**QNmax_b

!Constrain QN between QNmin and QNmax due to numerical issues
QN = max(min(QN, QNmax), QNmin)

!(Qmax - Qmin) [mmol N mmol C]:
dQN = QNmax - QNmin

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions:
!mu0 should be a function of alphaChl
muT = mu0 * exp(a1 * (alphaChl_ - .1)) !0.1 is the average alphaChl value

!Temperature dependent maximal growth rate at 1 ug C cell-1
muT = temp_Topt(Temp, muT, Topt_)

!Apply the size-scaling relationship following Wirtz (2011)
muT = Pmax_size(ESD_, muT)

!Assuming the same temperature dependence of nutrient uptake rate as on photosynthesis rate.
!QNmax/QNmin may be also a function of temperature which needs being further investigated.
Vcref  = muT * QNmax

!Assume the same temperature dependence of respiration as on photosynthetic rate (long-term adaptation; Barton et al. 2020):
RcT   = temp_Topt(Temp, Rc,   Topt_)
RcT   = respiration(ESD_, RcT)

RNT   = temp_Topt(Temp, RN,   Topt_)
RNT   = respiration(ESD_, RNT)

RChlT = temp_Topt(Temp, RChl, Topt_)
RChlT = respiration(ESD_, RChlT)

!Nutrient limitation [nd]:
Lno3 = (QN - QNmin) / dQN

if (Lno3 .le. 0d0) then
   PC = 0d0
else
   !Maximal photosynthesis rate (regulated by QN) [d-1]:
   PCmax = muT * Lno3

   !Light saturation parameter [W m-2 d-1]:
   Ik = PCmax / alphachl_ / theta

   !Calculate the fraction of open PSU [nd]:
   if (PAR > 0.) then !Photoinhibition
      A = Ainf(PAR, alphachl_, QN, QNmin, QNmax, theta)
   else
      A = 1d0
   endif

   !Light limitation index [nd]:
   SI = -A*PAR/Ik

   if (abs(SI) < 1d-10) then
      SI = 0.d0
   else
      SI = 1.d0 - exp(SI)
   endif

   !Photosynthesis rate [d-1]:
   PC = PCmax * SI
Endif

!Define rhochl [g Chl mol C-1]: fraction of phytoplankton carbon production that is devoted to Chl synthesis.
!If dark, assume that rhochl equaled the value calculated for the end of the preceding light period.
if (PAR <= 0d0) then
   rhochl   = rhoChl_L
else
   rhochl   = thetaNmax * PC / alphachl_ / theta / PAR
   rhoChl_L = rhochl
endif

!DIN uptake rate by phytoplankton [mol N mol C-1 d-1]:
VCN = Vcref * (NO3 - NO3_min)/ (NO3 + KN) * ((QNmax - QN) / dQN)**nx  !Vcref already temperature dependent
VCN = max(VCN, 0d0)

!Changes of cellular carbon [d-1]:
dC = C * (PC - zeta*VCN - RcT)

!Changes of cellular nitrogen [pmol N cell-1 d-1]:
!RNT has to be zero to avoid continuous decline of N per cell
dN = N * (VCN/QN - RNT)

!Changes of cellular Chl [d-1]:
dChl = Chl * (rhochl*VCN / theta - RChlT)

return
END subroutine GMK98_Ind_TempSizeLight