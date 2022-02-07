SUBROUTINE BIOLOGY
!1D lagrangian model using the model of Geider et al. L&O (1998)
USE params
USE state_variables
USE forcing,           only : Temp, TEMPBOL
USE grid,                only : Hz, nlev
USE Time_setting, only : dtdays, sec_of_day
implicit none
INTEGER :: k, i, j, m
INTEGER :: N_ = 0   !Number of particles in each grid
real    :: NO3 = 0.
real    :: PHY = 0.
real    :: PHYC=0. 
real    :: ZOO = 0. 
real    :: DET = 0. 
real    :: CHL = 0.
real    :: tf_p  = 0.
real    :: tf_z  = 0.
real    :: Graz  = 0.
real    :: dC_   = 0.
real    :: dN_   = 0.
real    :: dChl_ = 0.
real    :: uptake= 0.   !Total NO3 uptake
real    :: NPPc_(nlev)  = 0.  !C-based phytoplankton production (mg C m-3 d-1)
real    :: pp_ZP = 0.   
real    :: pp_DZ = 0.   
real    :: pp_ND = 0.   
real    :: pp_NZ = 0.   
real    :: Pmort  = 0.   
real    :: RES     = 0.   
real    :: EGES  = 0.   
real    :: gbar     = 0.   
real    :: INGES = 0.   
real    :: Zmort = 0.   
real,    parameter   :: GGE   = 0.30d0   !Zooplankton Gross Growth Efficiency
real,    parameter   :: unass = 0.24d0   !The fraction of unassimilated food ingested by zooplankton

! cellular carbon content threshold for division (pmol)
INTEGER, ALLOCATABLE :: index_(:)    !The indexes of particles in each grid
INTEGER, ALLOCATABLE :: scratch(:)    !The scratch indexes of particles in each grid
INTEGER                                 :: Allocatestatus = 0
!End of declaration

!Eulerian model for NO3, ZOO and DET and lagrangian model for PHY

DO k =  nlev, 1, -1
   tf_p  = TEMPBOL(Ep, Temp(k))  
   NO3 = t(iNO3, k)
   ZOO= t(iZOO, k)
   DET = t(iDET, k)

   if (sec_of_day == 0) then
       Varout(oNPP, k) = NPPc_(k)  !NPP of the past day; this is real NPP (mg C d-1 m-3)
       NPPc_(k)   = 0d0      !Reset NPPc_
   endif
   !The codes from Line 54-92 calculate the total amount of concentrations of 
   !phytoplankton carbon, nitrogen, and chl based on the cells present.
   PHYC= 0.d0
   PHY = 0.d0
   CHL = 0.d0
   N_  = 0

   !Get the indexes of particles in this grid (for grazing loss)
   allocate(index_(0), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating index_***"
   do i = 1, N_PAR
      if (p_PHY(i)%iz == k .and. p_PHY(i)%alive) then !Ignore dead super-individuals
         N_       = N_   + 1
         PHYC = PHYC + p_PHY(i)%num *1d-9*p_PHY(i)%C 
         PHY    = PHY  + p_PHY(i)%num *1d-9*p_PHY(i)%N
         CHL    = CHL  + p_PHY(i)%num *1d-9*p_PHY(i)%Chl

         if (N_ == 1) then
	         allocate(scratch(1), stat=Allocatestatus)
	         IF (AllocateStatus /= 0) STOP "*** Problem in allocating scratch***"
            scratch(1) = i

            !Update index_ and scratch has been deallocated
            call move_alloc(scratch, index_)   
         else
	         allocate(scratch(size(index_) + 1), stat=Allocatestatus)
	         IF (AllocateStatus /= 0) STOP "*** Problem in allocating scratch***"
	         scratch(1:size(index_))   = index_
            scratch(size(index_) + 1) = i

            !Update index_ and scratch has been deallocated
            call move_alloc(scratch, index_)   
         endif
      endif
   enddo

   PHYC = PHYC /Hz(k)   !Convert Unit to mmol/m^3
   PHY  = PHY  /Hz(k)   !Convert Unit to mmol/m^3
   CHL  = CHL  /Hz(k)   !Convert Unit to mmol/m^3
   !write(6,*) 'Total N before  = ', PHY + ZOO + NO3 + DET

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   ! In the NPZD model, phytoplankton cells utilize DIN and are eaten by zooplankton. 
   !The ingested food by zooplankton has three fates: 
   !1) being recycled to DIN; 2) being converted to detritus; and 3) supporting zooplankton growth. 
   !When zooplankton die, they are converted to detritus which is recycled to DIN and also sinks.

   tf_z      = TEMPBOL(Ez,Temp(k))
   gbar     = PHY**2/(PHY**2 + Kp**2)
   INGES = ZOO*gmax*tf_z*gbar
   Zmort  = ZOO*ZOO*mz*tf_z         !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES    = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES   = INGES*unass

   ! For production/destruction matrix:
   pp_ND  = RDN*DET * tf_z   
   pp_NZ  = ZOO*RES        
   pp_DZ  = ZOO*EGES + Zmort 
   pp_ZP  = ZOO*INGES      
  
   t(iZOO,k) = ZOO + dtdays*(pp_ZP - pp_DZ - pp_NZ)
   Varout(iZOO, k) = t(iZOO, k)

   !Now calculate new cell numbers associated with each particle
   if (PHY > 0d0) then
      Graz = pp_ZP/PHY*dtdays  !Specific mortality rate induced by zooplankton grazing per time step
   else
      Graz = 0d0 
   endif

   ! Impose the zooplankton grazing (the number of cells associated with each superindividual changes)
   IF (N_ > 0) THEN
      do j = 1, N_
         i = index_(j)
         p_PHY(i)%num = p_PHY(i)%num*(1d0 - Graz)   !Apply grazing
      enddo
   ENDIF

   ! Calculate total phytoplankton nitrogen uptake, mortality,  and PP (must after calculation of num(t+dt))
   uptake  = 0d0
   Pmort   = 0d0
   if (N_ > 0) then
      do j = 1, N_
         i = index_(j)
         call GMK98_Ind(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,         &
                        p_PHY(i)%C,    p_PHY(i)%N,   p_PHY(i)%Chl, dC_, dN_, dChl_)

         uptake      = uptake + dN_ * p_PHY(i)%num 
         NPPc_(k) = NPPc_(k) + dC_ * p_PHY(i)%num *1d-9/Hz(k)*12.d0   !Unit: mgC m-3 d-1

         ! Update cellular C, N, and Chl
         p_PHY(i)%C   =  p_PHY(i)%C   + dC_   * dtdays
         p_PHY(i)%N   =  p_PHY(i)%N   + dN_ * dtdays
         p_PHY(i)%Chl=  p_PHY(i)%Chl + dChl_ * dtdays

         ! If celular carbon is lower than the susbsistence threshold, it dies:
         p_PHY(i)%Cmin = 0.25d0 * p_PHY(i)%Cdiv
      
         if (p_PHY(i)%C < p_PHY(i)%Cmin) then  ! Dies
            Pmort = Pmort + p_PHY(i)%N * p_PHY(i)%num !Natural mortality of phytoplankton ==> DET
            p_PHY(i)%C     = 0d0
            p_PHY(i)%N     = 0d0
            p_PHY(i)%Chl  = 0d0
            p_PHY(i)%num = 0d0
            p_PHY(i)%alive = .false.

         endif
      enddo
  endif
  uptake= uptake*1d-9/Hz(k)!Convert uptake to mmol N m-3
  Pmort = Pmort*1d-9/Hz(k) !Convert Pmort to mmol N m-3

  !Now calculate NO3 and DET
  t(iNO3,k) = NO3 + dtdays*(pp_ND + pp_NZ - uptake)
  Varout(iNO3, k) = t(iNO3, k)
  t(iDET,k) = DET + Pmort + dtdays*( pp_DZ - pp_ND)
  Varout(iDET,k) = t(iDET,k)

   if (allocated(index_)) deallocate(index_)
ENDDO

!Needs to update t(iPN, :), t(iPC,:), and t(iChl,:)
call Par2PHY

END SUBROUTINE BIOLOGY

SUBROUTINE GMK98_Ind(Temp, PAR, NO3, C, N, Chl, dC, dN, dChl)
USE forcing, only : TEMPBOL
USE params,  only : Ep, aI0, thetaNmax, mu0, KN, rhoChl_L
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

subroutine Par2PHY
use state_variables, only : t, N_PAR, iPC, iPN, iChl, p_PHY, Varout
use grid,                   only : Hz, nlev
implicit none

!This subroutine calculate the total amount of concentrations of 
!phytoplankton carbon, nitrogen, and chl based on the cells present.
real      :: PHYC, PHY, CHL
integer :: k,i

do k = 1, nlev
   PHYC= 0.d0
   PHY = 0.d0
   CHL = 0.d0

   do i = 1, N_PAR
      if (p_PHY(i)%iz == k) then
         PHYC = PHYC + p_PHY(i)%num *1d-9*p_PHY(i)%C 
         PHY  = PHY  + p_PHY(i)%num *1d-9*p_PHY(i)%N
         CHL  = CHL  + p_PHY(i)%num *1d-9*p_PHY(i)%Chl
      endif

      ! If cellular carbon is above the division threshold, it divides
      if (p_PHY(i)%C >= p_PHY(i)%Cdiv) then  !Divide
         p_PHY(i)%C   = p_PHY(i)%C/2d0
         p_PHY(i)%N   = p_PHY(i)%N/2d0
         p_PHY(i)%Chl = p_PHY(i)%Chl/2d0
         p_PHY(i)%num = p_PHY(i)%num*2d0
      endif

   enddo

   t(iPC, k) = PHYC /Hz(k)   !Convert Unit to mmol/m^3
   t(iPN, k) = PHY/Hz(k)   !Convert Unit to mmol/m^3
   t(iChl,k) = CHL/Hz(k)   !Convert Unit to mmol/m^3

   Varout(iPC, k) = t(iPC, k)
   Varout(iPN, k) = t(iPN, k)
   Varout(iChl, k) = t(iChl, k)
enddo
return
end subroutine Par2PHY