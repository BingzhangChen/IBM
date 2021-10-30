SUBROUTINE BIOLOGY
!1D lagrangian model using the model of Geider et al. L&O (1998)
USE params
USE state_variables
USE forcing,      only : Temp, TEMPBOL
USE grid,         only : Hz, nlev
USE Time_setting, only : dtdays
implicit none
INTEGER :: k, i, j
INTEGER :: N_ = 0   !Number of particles in each grid
real    :: NO3, PHY, PHYC, ZOO, DET, CHL
real    :: tf_p  = 0.
real    :: tf_z  = 0.
real    :: Graz  = 0.
real    :: dC_   = 0.
real    :: dN_   = 0.
real    :: dChl_ = 0.
real    :: uptake= 0.   !Total NO3 uptake
real    :: NPP_  = 0.   !Net primary production
real    :: pp_ZP = 0.   
real    :: pp_DZ = 0.   
real    :: pp_ND = 0.   
real    :: pp_NZ = 0.   
real    :: RES   = 0.   
real    :: EGES  = 0.   
real    :: gbar  = 0.   
real    :: INGES = 0.   
real    :: Zmort = 0.   
real,    parameter   :: GGE   = 0.30d0   !Zooplankton Gross Growth Efficiency
real,    parameter   :: unass = 0.24d0   !The fraction of unassimilated food ingested by zooplankton

! cellular carbon content threshold for division (pmol)
INTEGER, ALLOCATABLE :: index_(:)    !The indexes of particles in each grid
INTEGER, ALLOCATABLE :: scratch(:)    !The scratch indexes of particles in each grid
INTEGER              :: Allocatestatus = 0
!End of declaration

!Eulerian model for NO3, ZOO and DET and lagrangian model for PHY
DO k = 1, nlev
   tf_p= TEMPBOL(Ep, Temp(k))  
   NO3 = t(iNO3, k)
   ZOO = t(iZOO, k)
   DET = t(iDET, k)

   !The codes from Line 45-88 calculate the total amount of concentrations of 
   !phytoplankton carbon, nitrogen, and chl based on the cells present.
   PHYC= 0.d0
   PHY = 0.d0
   CHL = 0.d0
   N_  = 0

   !Get the indexes of particles in this grid (for grazing loss)
   allocate(index_(0), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating index_***"
   do i = 1, N_PAR
      if (p_PHY(i)%iz == k) then
         N_   = N_   + 1
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

         PHYC = PHYC + p_PHY(i)%num *1d-9*p_PHY(i)%C 
         PHY  = PHY  + p_PHY(i)%num *1d-9*p_PHY(i)%N
         CHL  = CHL  + p_PHY(i)%num *1d-9*p_PHY(i)%Chl
      endif
   enddo

   PHYC = PHYC /Hz(k)   !Convert Unit to mmol/m^3
   PHY  = PHY  /Hz(k)   !Convert Unit to mmol/m^3
   CHL  = CHL  /Hz(k)   !Convert Unit to mmol/m^3

   ! Save total PHYtoplankton biomass and CHL
   Varout(iPN, k) = PHY
   Varout(iPC, k) = PHYC
   Varout(iCHL,k) = CHL

   ! Loss rate of phytoplankton to detritus equals to zero

   ! Calculate total phytoplankton uptake and PP
   uptake = 0d0
   NPP_   = 0d0

   if (N_ > 0) then
      do j = 1, N_
         i = index_(j)
         call GMK98_Ind(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,         &
                        p_PHY(i)%C,    p_PHY(i)%N,   p_PHY(i)%Chl,         &
                        dC_, dN_, dChl_)
         uptake = uptake + dN_ * p_PHY(i)%num *1d-9/Hz(k)
         NPP_   = NPP_   + dC_ * p_PHY(i)%num *1d-9/Hz(k)

         ! Update cellular C, N, and Chl
         p_PHY(i)%C  =  p_PHY(i)%C   + dC_   * dtdays
         p_PHY(i)%N  =  p_PHY(i)%N   + dN_   * dtdays
         p_PHY(i)%Chl=  p_PHY(i)%Chl + dChl_ * dtdays
      enddo
   else
      uptake = 0d0
      NPP_   = 0d0
   endif

   Varout(oNPP,k) = NPP_*12.     !Unit: mgC m-3 d-1

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   ! In the NPZD model, phytoplankton cells utilize DIN and are eaten by zooplankton. 
   !The ingested food by zooplankton has three fates: 
   !1) being recycled to DIN; 2) being converted to detritus; and 3) supporting zooplankton growth. 
   !When zooplankton die, they are converted to detritus which is recycled to DIN and also sinks.

   tf_z   = TEMPBOL(Ez,Temp(k))
   gbar   = PHY**2/(PHY**2 + Kp**2)
   INGES  = ZOO*gmax*tf_z*gbar
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
  
   Varout(iDET,k) = DET + dtdays*(pp_DZ - pp_ND)
   Varout(iNO3,k) = NO3 + dtdays*(pp_ND + pp_NZ - uptake)
   Varout(iZOO,k) = ZOO + dtdays*(pp_ZP - pp_DZ - pp_NZ)

   if (N_ > 0) then
      Graz = pp_ZP/PHY  !Specific grazing rate
   else
      Graz = 0d0 
   endif

   ! Impose the zooplankton grazing (the number of cells associated with each superindividual changes)
   ! Update Graz
   Graz = Graz * dtdays

   IF (N_ > 0) THEN
      do j = 1, N_
         p_PHY(j)%num = p_PHY(j)%num*(1d0 - Graz)
      enddo
   ENDIF

   if (allocated(index_)) deallocate(index_)

ENDDO

! If cellular carbon is above the division threshold, it divides
DO i = 1, N_PAR
   if (p_PHY(i)%C >= Cdiv) then  !Divide
      p_PHY(i)%C   = p_PHY(i)%C/2d0
      p_PHY(i)%N   = p_PHY(i)%N/2d0
      p_PHY(i)%Chl = p_PHY(i)%Chl/2d0
      p_PHY(i)%num = p_PHY(i)%num*2d0
   endif
ENDDO
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
real    :: alphaChl = 0.

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

QN       = N/C
theta    = Chl/C
Vcref    = mu0 * QNmax
alphaChl = aI0

!Temperature coefficient
tf_p = TEMPBOL(Ep, Temp)  

! N limitation
Lno3 = (QN-QNmin)/dQN

!Maximal photosynthesis rate (regulated by QN)
PCmax = mu0 *tf_P*Lno3

Ik = PCmax/alphaChl/theta

!The light limitation index (fpar)
SI = 1.d0 - exp(-par/Ik)

PC = PCmax * SI

! define rhochl (gChl/molC; the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
!If dark, assume that rhochl equaled the value calculated for the end of the preceding light period 
if (PAR <= 0d0) then
   rhochl   = rhoChl_L
else
   rhochl   = thetaNmax*PC/alphaChl/theta/PAR
   rhoChl_L = rhochl
endif

! DIN uptake rate by phytoplankton (molN/molC/d)
VCN = Vcref * NO3/(NO3 + KN)* (QNmax-QN)/dQN*tf_p

! Changes of cellular carbon
dC  = C*(PC - zeta*VCN - Rc*tf_p)

! Changes of cellular nitrogen
dN  = N*(VCN/QN - RN*tf_p)

! Changes of cellular Chl
dChl= Chl*(rhochl*VCN/theta - RChl*tf_p)

return
end subroutine GMK98_Ind