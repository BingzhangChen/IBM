SUBROUTINE BIOLOGY
!1D lagrangian model using the model of Geider et al. L&O (1998)
USE params
USE state_variables
USE forcing,          only : Temp, PAR
USE Trait_functions,  only : TEMPBOL, PHY_C2Vol, palatability
USE grid,             only : Hz, nlev
USE Time_setting,     only : dtdays, sec_of_day
implicit none
INTEGER :: k, i, j, m,kk
INTEGER :: N_ = 0   !Number of particles in each grid
real    :: Abun_ = 0   !Total abundance in each grid
real    :: NO3 = 0.
real    :: ZOO(NZOO) = 0. 
real    :: DET = 0. 
real    :: tf_z  = 0.
real    :: Graz  = 0.
real    :: dC_   = 0.
real    :: dN_   = 0.
real    :: dChl_ = 0.
real    :: uptake= 0.   !Total NO3 uptake
real    :: NPPc_(nlev)  = 0.  !C-based phytoplankton production (mg C m-3 d-1)
real    ::  IPAR(nlev)  = 0.  !Daily integrated PAR at each depth
real    :: pp_DZ = 0.   
real    :: pp_ND = 0.   
real    :: Pmort  = 0.   
real    :: Cmin  = 0. !Phytoplankton subsistence carbon quota below which the cell will die   
real    :: FZoo(NZOO) = 0.   !The total amount of palatable prey (in Nitrogen)
                             !for each zooplankton size class
real    :: phyV = 0.   !Phytoplankton cell volume
real    :: RES     = 0.   
real    :: EGES  = 0.   
real    :: gbar  = 0.   
real    :: INGES(NZOO) = 0.   
real    :: Zmort = 0.   
real    :: Gmatrix(NZOO,NZOO) = 0.d0     !Grazer biomass specific grazing rate matrix
real,    allocatable :: BN(:)                !The amount of nitrogen in each super-individual
real,    allocatable :: Pmatrix(:,:)     !Phytoplankton mortality rates by each zooplankton size class for each superindividual
real,    parameter   :: eta     = -1.d0*6.6  !Prey refuge parameter for nitrogen
real,    parameter   :: A_g    = 21.9   !Intercept of the allometric equation of maximal zooplankton grazing rate (Ward et al. 2012)
real,    parameter   :: B_g    = -0.16  !Slope of the allometric equation of maximal zooplankton grazing rate (Ward et al. 2012)
real,    parameter   :: mz_g   = 0.d0  !Power of zooplankton mortality following Ward et al. (2013)
real   :: Nt_min = 0.0 !Minimal amount carbon of each super-individual (number of cells * N content per cell)

! cellular carbon content threshold for division (pmol)
INTEGER, ALLOCATABLE :: index_(:)    !The indexes of particles in each grid
INTEGER, ALLOCATABLE :: scratch(:)    !The scratch indexes of particles in each grid
INTEGER              :: Allocatestatus = 0
real :: PHY_t = 0d0  !Total phytoplankton N

!End of declaration

!Update Nt_min
!Compute total phytoplankton nitrogen
PHY_t = 0d0
DO k = nlev, 1, -1
   PHY_t = PHY_t + t(iPN, k) * Hz(k)
ENDDO

!Minimal N of each superindividual should be 0.1% of the average
Nt_min = PHY_t*1d9/dble(N_PAR) * 0.001

!Eulerian model for NO3, ZOO and DET and lagrangian model for PHY
DO k = nlev, 1, -1
   NO3 = t(iNO3, k)
   DET = t(iDET, k)

   !Convert DET to NO3 at bottom
   if (k .eq. 1) then
      NO3 = NO3 + DET
      DET = 0.d0
   endif
   Varout(oTEMP,k) = Temp(k)
   IPAR(k) = IPAR(k) + PAR(k)*dtdays   !Unit: W m-2

   do kk = 1, NZOO
      ZOO(kk)= t(iZOO(kk), k)
   enddo

   if (sec_of_day == 0) then
       Varout(oNPP, k) = NPPc_(k)  !NPP of the past day; this is real NPP (mg C d-1 m-3)
       Varout(oPAR, k) = IPAR(k)   !Integrated PAR of the past day (W m-2)
       NPPc_(k)   = 0d0              !Reset NPPc_
       IPAR(k)    = 0d0              !Reset IPAR
   endif

   !Calculate the number of super-individuals (N_) in this vertical layer and obtain their indexes (index_)
   N_  = 0

   !Get the indexes of particles in this grid (for grazing loss)
   allocate(index_(0), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating index_***"

   do i = 1, N_PAR
      if (p_PHY(i)%iz == k .and. p_PHY(i)%alive) then !Ignore dead super-individuals
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
      endif
   enddo

   !Save number of super-individuals per m3
   Varout(oN_ind, k) = dble(N_)/Hz(k)

   !Allocate Pmatrix (matrix for super-individual grazing mortality)
   
   !Reset total abundance
   Abun_ = 0d0

   IF (N_ > 0) THEN
      allocate(Pmatrix(N_, NZOO), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating Pmatrix***"
      Pmatrix(:,:) = 0d0

      allocate(BN(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating BN***"
      BN(:) = 0d0

      ! calculate the amount of nitrogen and total abundances (cells/m3) in each super-individual
      DO m = 1, N_

         i = index_(m)

         !The amount of  N in the super-individual m
         BN(m) = p_PHY(i)%N * p_PHY(i)%num*1d-9/Hz(k)

         !Count the number of cells in each layer
         Abun_ = Abun_ + p_PHY(i)%num
      ENDDO

   ENDIF
   Varout(oN_cell, k) = Abun_/Hz(k) !Abundances (cells m-3)

   !The multiple zooplankton size class model follows Ward et al. L&O 2012
   ! In the NPZD model, phytoplankton cells utilize DIN and are eaten by zooplankton. 
   !The ingested food by zooplankton has three fates: 
   !1) being recycled to DIN; 2) being converted to detritus; and 3) supporting zooplankton growth. 
   !The natural mortality of zooplankton are converted to detritus which is recycled to DIN and also sinks.

   tf_z = TEMPBOL(Ez,Temp(k))  !Temperature coefficient of zooplankton

   !Calculate the total amount of prey N biomass available to each size class of zooplankton
   DO kk = 1, NZOO
      gmax = A_g * VolZOO(kk)**B_g 
      FZoo(kk) = 0d0

      !First calculate total phyto. prey from super-individuals
      IF (N_ > 0) THEN
         do m = 1, N_

           i = index_(m)

           !Volume of phytoplankton super-individual
           phyV = PHY_C2Vol(p_PHY(i)%C)

           !The amount of patalable prey in the super-individual m
           Pmatrix(m,kk) = palatability(VolZOO(kk), phyV, SDZoo) * BN(m)

           !Calculate the palatability of each prey superindividual and add to the total amount palatable prey
           FZoo(kk) = FZoo(kk) + Pmatrix(m,kk)
         enddo
      ENDIF

	  !Second, calculate the total zooplankton prey
	  IF (kk > 1) THEN
	    do m = 1, (kk - 1)

       !Save the palatability into Gmatrix
       Gmatrix(m,kk) = palatability(VolZOO(kk), VolZOO(m), SDZoo) 

		 !Calculate the palatability of each zoo. prey and add to the total amount palatable prey
		 FZoo(kk) = FZoo(kk) + Gmatrix(m,kk) * ZOO(m)

	    enddo
	  ENDIF

     !Save the total available prey for zooplankton
     !Varout(oFZ(kk),k) = FZOO(kk)

     gbar = FZoo(kk)/(FZoo(kk) + Kp)*(1.d0 - exp(eta *FZoo(kk)))

     !Total ingestion of zooplankton kk (mmol N m-3 d-1)
     INGES(kk) = ZOO(kk)*gmax*tf_z*gbar

	  IF (kk > 1) THEN
	    do m = 1, (kk - 1)

         !Calculate the total ingestion rate (mmol N m-3 d-1) of zooplankton kk on zooplankton m
         if (FZOO(kk) > 0d0) then
            Gmatrix(m,kk) = Gmatrix(m,kk)*ZOO(m)/FZoo(kk)*INGES(kk)
         else
            Gmatrix(m,kk) = 0d0
         endif

        enddo
	  ENDIF
   ENDDO !End of the zooplankton loop

   !Computing zooplankton mortality
   RES   = 0d0  !Total amount of nitrogen that is excreted by zooplankton and becomes DIN
   EGES  = 0d0  !Total egestion by zooplankton to detritus
   pp_DZ = 0d0 !Flux from ZOO to DET 

   DO kk = 1, NZOO

      !Zooplankton excretion rate (-> DIN)
      RES  = RES + INGES(kk)*(1d0-GGE-unass)

      !ZOOPLANKTON EGESTION (-> Detritus)
      EGES = EGES + INGES(kk)*unass

      pp_DZ = pp_DZ + INGES(kk)*unass

      !Calculate zooplankton mortality
      Zmort = ZOO(kk)* mz *tf_z * VolZOO(kk)**mz_g   !zooplankton Mortality term due to natural death

      pp_DZ = pp_DZ + Zmort

      !Loop through all predators
      if (kk .lt. NZOO) then
        do m = (kk + 1), NZOO
          Zmort = Zmort + Gmatrix(kk,m)    !Natural Mortality + grazing by other ZOO
        enddo
      endif

      !Update the biomass of ZOOplankton kk
      t(iZOO(kk),k) = max(ZOO(kk) + dtdays*(GGE*INGES(kk) - Zmort), 0d0)
      Varout(iZOO(kk), k) = t(iZOO(kk), k)

   ENDDO !End of the zooplankton loop

   ! For production/destruction matrix:
   pp_ND = RDN*DET*tf_z                   !Flux from DET to DIN

   !Now calculate new cell numbers associated with each particle
   ! Impose the zooplankton grazing (the number of cells associated with each superindividual changes)
   IF (N_ > 0) THEN
 
      do j = 1, N_
         i = index_(j)

         Graz = 0d0 

         !Calculate all the zooplankton ingestion for this superindividual (unit: mmol N m-3 d-1)
         do m = 1, NZOO
              if (FZOO(m) > 0d0) Graz = Graz + INGES(m) * Pmatrix(j, m)/FZoo(m)
         enddo

         p_PHY(i)%num = p_PHY(i)%num*(1d0 - Graz*dtdays/BN(j))   !Apply grazing
      enddo
   ENDIF

   ! Calculate total phytoplankton nitrogen uptake, mortality,  and PP (must after calculation of num(t+dt))
   uptake  = 0d0
   Pmort   = 0d0
   if (N_ > 0) then
      do j = 1, N_
         i = index_(j)

         !Check forcing fields:
         if (p_PHY(i)%NO3 .le. 0.d0) then
           write(6,*) "Nonpositive Nitrate found for this particle!"
           write(6,*) 'Particle at grid ',  p_PHY(i)%iz
           write(6,*) 'Particle at depth ', p_PHY(i)%rz
           write(6,*) 'NO3 of this particle ', p_PHY(i)%NO3
           stop
         endif

         if (p_PHY(i)%PAR .lt. 0.d0) then
           write(6,*) "Negative PAR found for this particle!"
           write(6,*) 'Particle at grid ',  p_PHY(i)%iz
           write(6,*) 'Particle at depth ', p_PHY(i)%rz
           write(6,*) 'PAR of this particle ', p_PHY(i)%PAR
           stop
         endif
 
         SELECTCASE (Model_ID)
         CASE(GMK98_simple)
            call GMK98_Ind(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,         &
                           p_PHY(i)%C,    p_PHY(i)%N,   p_PHY(i)%Chl, dC_, dN_, dChl_)
         CASE(GMK98_Topt)
            call GMK98_Ind_Temp(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,  p_PHY(i)%Topt, &
                           p_PHY(i)%C,    p_PHY(i)%N,   p_PHY(i)%Chl,  dC_, dN_, dChl_)

         CASE(GMK98_Light)
            stop "To be developed..."
         CASE(GMK98_Size)
            stop "To be developed..."
         CASE(GMK98_SizeLight)
            stop "To be developed..."
         CASE(GMK98_ToptLight)
            stop "To be developed..."
         CASE(GMK98_ToptSize)
            call GMK98_Ind_TempSize(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,  p_PHY(i)%Topt,&
                           p_PHY(i)%C, p_PHY(i)%N, p_PHY(i)%Chl, p_PHY(i)%CDiv, dC_, dN_, dChl_)
         CASE(GMK98_ToptSizeLight)

            call GMK98_Ind_TempSizeLight(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,  p_PHY(i)%Topt,&
                 p_PHY(i)%C, p_PHY(i)%N, p_PHY(i)%Chl, p_PHY(i)%CDiv, exp(p_PHY(i)%LnalphaChl),&
                 dC_, dN_, dChl_)

         CASE DEFAULT
            stop "Model choice is wrong!!"
         END SELECT

         uptake   =   uptake + dN_ * p_PHY(i)%num 
         NPPc_(k) = NPPc_(k) + dC_ * p_PHY(i)%num *1d-9/Hz(k)*12.d0*dtdays !Unit: mgC m-3 d-1

         ! Update cellular C, N, and Chl
         p_PHY(i)%C   =  p_PHY(i)%C   + dC_   * dtdays
         p_PHY(i)%N   =  p_PHY(i)%N   + dN_   * dtdays
         p_PHY(i)%Chl =  p_PHY(i)%Chl + dChl_ * dtdays

         ! If celular carbon is lower than the susbsistence threshold (Cmin), it dies:
         Cmin = 0.25d0 * p_PHY(i)%Cdiv

         if (p_PHY(i)%C < Cmin .or. (p_PHY(i)%num*p_PHY(i)%N) < Nt_min) then  ! The superindividual Dies
            N_death(k) = N_death(k) + 1
            Pmort = Pmort + p_PHY(i)%N * p_PHY(i)%num !Natural mortality of phytoplankton ==> DET
            p_PHY(i)%C   = 0d0
            p_PHY(i)%N   = 0d0
            p_PHY(i)%Chl = 0d0
            p_PHY(i)%num = 0d0
            p_PHY(i)%alive = .FALSE.
         endif
      enddo
  endif !End if of N_ > 0

  uptake= uptake*1d-9/Hz(k) !Convert uptake to mmol N m-3
  Pmort =  Pmort*1d-9/Hz(k) !Convert Pmort to mmol N m-3

  !Now calculate NO3 and DET
  t(iNO3,k) = NO3 + dtdays*(pp_ND + RES - uptake)

  if (t(iNO3, k) < 0.) then
    write(6,*) 'Current grid = ', k
    write(6,*) 'Number of PHY particles = ', N_
    write(6,*) 'Previous NO3 = ', NO3
    write(6,*) 'Phyto. uptake = ', dtdays*uptake
    write(6,*) 'DET reg. = ', dtdays*pp_ND
    write(6,*) 'ZOO excretion = ', dtdays*RES

    do j = 1, N_
      i = index_(j)
      call GMK98_Ind_TempSizeLight(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,  p_PHY(i)%Topt,&
                 p_PHY(i)%C, p_PHY(i)%N, p_PHY(i)%Chl, p_PHY(i)%CDiv, exp(p_PHY(i)%LnalphaChl),&
                 dC_, dN_, dChl_)

      write(6,*) 'Particle ', i, ' num= ', p_PHY(i)%num
      write(6,*) 'Nutrient change ',  dN_
    enddo
    stop
  endif

  Varout(iNO3, k) = t(iNO3, k)

  t(iDET,k) = DET + Pmort + dtdays*(pp_DZ - pp_ND)

  Varout(iDET,k) = t(iDET,k)

  if (allocated(index_))  deallocate(index_)
  if (allocated(Pmatrix)) deallocate(Pmatrix)
  if (allocated(BN)) deallocate(BN)
ENDDO

!update t(iPN, :), t(iPC,:), and t(iChl,:) and compute mean trait and trait variance
call Par2PHY

END SUBROUTINE BIOLOGY

SUBROUTINE GMK98_Ind(Temp, PAR, NO3, C, N, Chl, dC, dN, dChl)
USE Trait_functions, only : TEMPBOL, palatability
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

!Adding optimal temperature into the Geider model
SUBROUTINE GMK98_Ind_Temp(Temp, PAR, NO3, Topt_, C, N, Chl, dC, dN, dChl)
USE Trait_functions, only : temp_Topt
USE params,              only : Ep, aI0, thetaNmax, mu0, KN, rhoChl_L
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

real, PARAMETER   :: nx = 1.d0  !Exponent of nutrient uptake function in GMK98

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

subroutine GMK98_Ind_TempSize(Temp, PAR, NO3, Topt_, C, N, Chl, Cdiv, dC, dN, dChl)

USE Trait_functions, only : temp_Topt, PHY_C2Vol
USE params,          only : Ep, aI0, thetaNmax, mu0, rhoChl_L

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

real, parameter   :: nx   = 1.d0   ! Exponent of nutrient uptake function in GMK98

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
SI = 1.d0 - exp(-PAR/Ik)
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

!------------------------------------------------------------------------------------------------
SUBROUTINE GMK98_Ind_TempSizeLight(Temp, PAR, NO3, Topt_, C, N, Chl, Cdiv, alphaChl_, dC, dN, dChl)
USE Trait_functions, only : temp_Topt, PHY_C2Vol, Ainf
USE params,          only : thetaNmax, mu0, rhoChl_L, QNmin_a, QNmin_b
USE params,          only : QNmax_a, QNmax_b, KN_a, KN_b
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

real, parameter   :: RC     = 0.0d0   !Basic C respiration rate [d-1]
real, parameter   :: RN     = 0.0d0   !Basic N respiration rate [d-1]
real, parameter   :: RChl   = 0.0d0   !Basic Chl respiration rate [d-1]
real              :: RcT    = 0.d0    !Temperature dependent C-based respiration rate [d-1]
real              :: RNT    = 0.d0    !Temperature dependent N-based respiration rate [d-1]
real              :: RChlT  = 0.d0    !Temperature dependent Chl-based respiration rate [d-1]

real, parameter   :: zeta   = 3.0d0   !Cost of biosynthesis [mol C mol N-1]

! Maximal specific N uptake as a function of temp. under resource (nutrient and light) replete conditions [mol N mol C-1 d-1]:
real              :: Vcref  = 0.

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions [uM]:
real              :: muT    = 0.

!Exponent of nutrient uptake function in GMK98.
!Factor that governs the down regulation of nutrient uptake as Q approaches Qmax:
real, parameter   :: nx     = 1.d0

!Kn is an allometric function of Vol (Cdiv) (Edwards et al. 2012) [uM]:
real              :: KN     = 0.          !Half-saturation constant [uM]

real, parameter   :: a1 = 0.d0        ! Allometric exponent between mumax and alphaChl
real, parameter   :: b0 = 0.d0        ! Allometric exponent between mumax and size 
real, parameter   :: b1 = 0.d0        ! Allometric exponent between mumax and size
real, parameter   :: b2 = 0.d0        ! Allometric exponent between mumax and size
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

!Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
KN = KN_a * Vol**KN_b

!Minimal N:C ratio [mmol N mmol C] following Ward et al. (2012):
QNmin = QNmin_a * Vol**QNmin_b

!Maximal N:C ratio [mmol N mmol C] following Ward et al. (2012):
QNmax = QNmax_a * Vol**QNmax_b

!Constrain QN between QNmin and QNmax due to numerical issues
QN = max(min(QN, QNmax), QNmin)

!(Qmax - Qmin) [mmol N mmol C]:
dQN   = QNmax - QNmin

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions:
!(what value of mu0 should be?): mu0 should be a function of alphaChl
muT    = mu0 * exp(a1 * (alphaChl_ - .1)) !0.1 is the average alphaChl value
muT    = temp_Topt(Temp, muT, Topt_)
muT    = muT * exp(b0 + b1 * log(Cdiv) + b2 * log(Cdiv)**2 )

!Assuming the same temperature dependence of nutrient uptake rate as on photosynthesis rate.
!QNmax/QNmin may be also a function of temperature which needs being further investigated.
Vcref  = muT * QNmax

!Assume the same temperature dependence of respiration as on photosynthetic rate (long-term adaptation; Barton et al. 2020):
RcT   = temp_Topt(Temp, Rc,   Topt_)
RNT   = temp_Topt(Temp, RN,   Topt_)
RChlT = temp_Topt(Temp, RChl, Topt_)

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
   SI = - A * PAR / Ik

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
dC = C * (PC - zeta * VCN - RcT)

!Changes of cellular nitrogen [pmol N cell-1 d-1]:
!RNT has to be zero to avoid continuous decline of N per cell
dN = N * (VCN / QN - RNT)

!Changes of cellular Chl [d-1]:
dChl = Chl * (rhochl * VCN / theta - RChlT)

return
END subroutine GMK98_Ind_TempSizeLight
!------------------------------------------------------------------------------------------------

SUBROUTINE Par2PHY
use state_variables, only : t, N_PAR, iPC, iPN, iChl, p_PHY, Varout, nu, sigma, iTopt, iSize, ialphaChl
use state_variables, only : NTrait, N_birth, N_death, N_mutate
use state_variables, only : oN_cell, oCDiv_avg, oCDiv_var, oTopt_avg, oTopt_var 
use state_variables, only : oLnalpha_var, oLnalpha_avg
use state_variables, only : oTalp_cov, oALnV_cov, oTLnV_cov
use grid,            only : Hz, nlev
use mGf90,           only : srand_mtGaus
IMPLICIT NONE

!This subroutine calculate the total amount of concentrations of 
!phytoplankton carbon, nitrogen, and chl based on the cells present.
real      :: PHYC(nlev) = 0d0
real      ::    PHY(nlev) = 0d0
real      ::    CHL(nlev) = 0d0
real      ::   mTopt_(nlev) = 0d0
real      ::   vTopt_(nlev) = 0d0
real      ::   mCDiv_(nlev) = 0d0
real      ::   vCDiv_(nlev) = 0d0
real      ::   mlnalpha_(nlev) = 0d0
real      ::   vlnalpha_(nlev) = 0d0
real      ::   cov_TA(nlev) = 0d0
real      ::   cov_TL(nlev) = 0d0
real      ::   cov_AL(nlev) = 0d0
real      ::  nu_ = 0d0   !Basic Mutation rate
real      :: cff = 0.d0   !Random number [0,1]
real      :: oldtt(1) = 0.   !Scratch variable for storing the old trait
real      :: newtt(1) = 0.   !Scratch variable for storing the new trait
real      :: vartt(1,1) = 0.   !Variance of the mutating trait
integer :: k,i,m, ipar

!End of declaration

PHYC(:)= 0.d0
PHY(:) = 0.d0
CHL(:) = 0.d0
mCDiv_(:) = 0d0
mTopt_(:) = 0d0
mlnalpha_(:) = 0d0

!loop through all super-individuals to convert into Eulerian concentrations
DO i = 1, N_PAR
   ipar = p_PHY(i)%iz  !The current grid of super-individual i
   !Handle cell division and mutation
   ! If cellular carbon is above the division threshold, it divides
   IF (p_PHY(i)%C >= p_PHY(i)%Cdiv) THEN  !Divide
      N_birth(ipar)= N_birth(ipar) + 1
      p_PHY(i)%C   = p_PHY(i)%C/2d0
      p_PHY(i)%N   = p_PHY(i)%N/2d0
      p_PHY(i)%Chl = p_PHY(i)%Chl/2d0
      p_PHY(i)%num = p_PHY(i)%num*2d0

      !Mutation
      DO m = 1, NTrait
         nu_ = p_PHY(i)%num*nu(m)
         call random_number(cff)

         IF (cff < nu_) THEN !Mutation occurs
            N_mutate(ipar) = N_mutate(ipar) + 1

            select case(m)
            case(iTopt)
                oldtt(1) = p_PHY(i)%Topt
            case(iSize)
                oldtt(1) = log(p_PHY(i)%CDiv)
            case(ialphaChl)
                oldtt(1) = p_PHY(i)%LnalphaChl
            case DEFAULT
                stop "Trait index wrong!"
            end select

            vartt(1,1)= sigma(m)**2   !Construct the covariance matrix for the selected trait

            !A new Topt is randomly sampled from a Gaussian distribution with mean of previous Topt and SD of sigma
            newtt = srand_mtGaus(1, oldtt, vartt)
            select case(m)
            case(iTopt)
                p_PHY(i)%Topt = newtt(1)
            case(iSize)
                p_PHY(i)%CDiv = exp(newtt(1))
            case(ialphaChl)
                p_PHY(i)%LnalphaChl = newtt(1)
            case DEFAULT
                stop "Trait index wrong!"
            end select
         ENDIF
      ENDDO
   ENDIF

   !Calculate Eulerian concentrations of phyto C, N, and Chl, mean trait  for each layer
   PHYC(ipar) = PHYC(ipar) + p_PHY(i)%num *p_PHY(i)%C 
   PHY(ipar)  = PHY(ipar)  + p_PHY(i)%num *p_PHY(i)%N
   CHL(ipar)  = CHL(ipar)  + p_PHY(i)%num *p_PHY(i)%Chl

   mCDiv_(ipar) = mCDiv_(ipar) + p_PHY(i)%num * log(p_PHY(i)%Cdiv)
   mTopt_(ipar) = mTopt_(ipar) + p_PHY(i)%num * p_PHY(i)%Topt
   mlnalpha_(ipar) = mlnalpha_(ipar) + p_PHY(i)%num * p_PHY(i)%LnalphaChl
ENDDO

do k = 1, nlev
  t(iPC, k) = PHYC(k) *1d-9/Hz(k)!Convert Unit to mmol/m^3
  t(iPN,k) = PHY(k) *1d-9/Hz(k)   !Convert Unit to mmol/m^3
  t(iChl,k) = CHL(k)*1d-9/Hz(k) !Convert Unit to mmol/m^3

  if (Varout(oN_cell,k) > 0.) then
    mCDiv_(k)      = mCDiv_(k)/(Varout(oN_cell, k) * Hz(k))
    mTopt_(k)      = mTopt_(k)/(Varout(oN_cell, k) * Hz(k))
    mlnalpha_(k) = mlnalpha_(k)/(Varout(oN_cell, k) * Hz(k))
  endif
enddo

Varout(iPC, :) = t(iPC, :)
Varout(iPN, :) = t(iPN, :)
Varout(iChl, :) = t(iChl,:)
Varout(oCDiv_avg, :) = mCDiv_(:)
Varout(oTopt_avg, :) = mTopt_(:)
Varout(oLnalpha_avg, :) = mLnalpha_(:)

!Compute trait variances
vCDiv_(:) = 0d0
vTopt_(:) = 0d0
vlnalpha_(:) = 0d0
cov_AL(:) = 0d0
cov_TA(:) = 0d0
cov_TL(:) = 0d0

DO i = 1, N_PAR
   ipar = p_PHY(i)%iz  !The current grid of super-individual i
   vCDiv_(ipar) = vCDiv_(ipar) + p_PHY(i)%num * (log(p_PHY(i)%Cdiv) - mCDiv_(ipar))**2
   vTopt_(ipar) = vTopt_(ipar) + p_PHY(i)%num * (p_PHY(i)%Topt - mTopt_(ipar))**2
   vlnalpha_(ipar) = vlnalpha_(ipar) + p_PHY(i)%num*(p_PHY(i)%LnalphaChl - mlnalpha_(ipar))**2
   cov_TL(ipar) = cov_TL(ipar) +p_PHY(i)%num * (log(p_PHY(i)%Cdiv) - mCDiv_(ipar))*(p_PHY(i)%Topt - mTopt_(ipar))
   cov_AL(ipar) = cov_AL(ipar)+p_PHY(i)%num * (log(p_PHY(i)%Cdiv) - mCDiv_(ipar))*(p_PHY(i)%LnalphaChl - mlnalpha_(ipar))
   cov_TA(ipar) = cov_TA(ipar)+p_PHY(i)%num * (p_PHY(i)%LnalphaChl - mlnalpha_(ipar))*(p_PHY(i)%Topt - mTopt_(ipar))
ENDDO

do k = 1, nlev
  if (Varout(oN_cell,k) > 0.) then
    vCDiv_(k)      = vCDiv_(k)/(Varout(oN_cell, k) * Hz(k))
    vTopt_(k)      = vTopt_(k)/(Varout(oN_cell, k) * Hz(k))
    vlnalpha_(k) = vlnalpha_(k)/(Varout(oN_cell, k) * Hz(k))
    cov_TL(k)      = cov_TL(k)/(Varout(oN_cell, k) * Hz(k))
    cov_TA(k)      = cov_TA(k)/(Varout(oN_cell, k) * Hz(k))
    cov_AL(k)      = cov_AL(k)/(Varout(oN_cell, k) * Hz(k))
  endif
enddo

Varout(oCDiv_var, :) = vCDiv_(:)
Varout(oTopt_var, :) = vTopt_(:)
Varout(oLnalpha_var, :) = vLnalpha_(:)
Varout(oALnV_cov, :) = cov_AL(:)
Varout(oTalp_cov, :) = cov_TA(:)
Varout(oTLnV_cov, :) = cov_TL(:)
return
END subroutine Par2PHY
