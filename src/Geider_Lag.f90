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
real    :: Pmort = 0.   
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
real,    allocatable :: BN(:)            !The amount of nitrogen in each super-individual
real,    allocatable :: BC(:)            !The amount of carbon in each super-individual
real,    allocatable :: Pmatrix(:,:)     !Phytoplankton mortality rates by each zooplankton size class for each superindividual
real,    parameter   :: eta    = -1.d0*6.6  !Prey refuge parameter for nitrogen
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
       NPPc_(k) = 0d0              !Reset NPPc_
       IPAR(k)  = 0d0              !Reset IPAR
   endif

   !Calculate the number of super-individuals (N_) in this vertical layer and obtain their indexes (index_)
   N_ = 0

   !Get the indexes of particles in this grid (for grazing loss)
   allocate(index_(0), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating index_***"

   do i = 1, N_PAR
      if (p_PHY(i)%iz == k .and. p_PHY(i)%alive) then !Ignore dead super-individuals
         N_ = N_ + 1

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
   
   !Reset total abundance
   Abun_ = 0d0

   IF (N_ > 0) THEN

      !Allocate Pmatrix (matrix for super-individual grazing mortality) for size-structured model
      If (NZOO > 1) Then
         allocate(Pmatrix(N_, NZOO), stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP "*** Problem in allocating Pmatrix***"
         Pmatrix(:,:) = 0d0
      Endif

      allocate(BN(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating BN***"
      BN(:) = 0d0

      allocate(BC(N_), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating BC***"
      BC(:) = 0d0

      ! calculate the amount of nitrogen and total abundances (cells/m3) in each super-individual
      DO m = 1, N_

         i = index_(m)

         !The amount of N in the super-individual m
         BN(m) = p_PHY(i)%N * p_PHY(i)%num*1d-9/Hz(k)

         !The amount of C in the super-individual m
         BC(m) = p_PHY(i)%C * p_PHY(i)%num*1d-9/Hz(k)

         !Count the number of cells in each layer
         Abun_ = Abun_ + p_PHY(i)%num
      ENDDO
   ENDIF

   Varout(oN_cell, k) = Abun_/Hz(k) !Abundances (cells m-3)

   !The multiple zooplankton size class model follows Ward et al. L&O 2012
   !In the NPZD model, phytoplankton cells utilize DIN and are eaten by zooplankton. 
   !The ingested food by zooplankton has three fates: 
   !1) being recycled to DIN; 2) being converted to detritus; and 3) supporting zooplankton growth. 
   !The natural mortality of zooplankton are converted to detritus which is recycled to DIN and also sinks.

   tf_z = TEMPBOL(Ez,Temp(k))  !Temperature coefficient of zooplankton

   !Calculate the total amount of prey N biomass available to each size class of zooplankton
   IF (NZOO > 1) THEN
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
         Varout(oFZ(kk),k) = FZOO(kk)

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
   ELSE !Only one zooplankton

      !Calculate the sum of total phytoplankton biomass
      FZOO(1) = 0d0
      DO m = 1, N_
         i = index_(m)
         FZOO(1) = BN(m) + FZOO(1)
      ENDDO

      !Save the total available prey for zooplankton
      Varout(oFZ(1),k) = FZOO(1)

      gbar = FZoo(1)/(FZoo(1) + Kp)*(1.d0 - exp(eta *FZoo(1)))

      !Total ingestion of zooplankton (mmol N m-3 d-1)
      INGES(1) = ZOO(1)*gmax*tf_z*gbar

   ENDIF

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
      If (NZOO > 1) Then
         if (kk .lt. NZOO) then
           do m = (kk + 1), NZOO
             Zmort = Zmort + Gmatrix(kk,m)    !Natural Mortality + grazing by other ZOO
           enddo
         endif
      Endif

      Varout(oZmort(kk), k) = Zmort

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

         !Calculate all the zooplankton ingestion for this superindividual (Graz, unit: mmol N m-3 d-1)
         if (NZOO > 1) then
            do m = 1, NZOO
                 if (FZOO(m) > 0d0) Graz = Graz + INGES(m) * Pmatrix(j, m)/FZoo(m)
            enddo
            p_PHY(i)%num = p_PHY(i)%num*(1d0 - Graz*dtdays/BN(j))   !Apply grazing to super-individual j
         else
            Graz = Graz + INGES(1)/FZoo(1)
            p_PHY(i)%num = p_PHY(i)%num*(1d0 - Graz*dtdays)   !Apply grazing to super-individual j
         endif

      enddo
   ENDIF

   ! Calculate total phytoplankton nitrogen uptake, mortality, and PP (must after calculation of num(t+dt))
   uptake  = 0d0
   Pmort   = 0d0
   if (N_ > 0) then
      do j = 1, N_
         i = index_(j)

         SELECTCASE (Model_ID)
         CASE(GMK98_simple)
            call GMK98_Ind(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3, &
                           p_PHY(i)%C,    p_PHY(i)%N,   p_PHY(i)%Chl, &
                           dC_, dN_, dChl_)

         CASE(GMK98_Topt)
            call GMK98_Ind_Temp(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3,  p_PHY(i)%Topt, &
                           p_PHY(i)%C,    p_PHY(i)%N,   p_PHY(i)%Chl,  dC_, dN_, dChl_)

         CASE(GMK98_Light)
            stop "To be developed..."
         CASE(GMK98_Size)
            call GMK98_Ind_Size(p_PHY(i)%Temp, p_PHY(i)%PAR, p_PHY(i)%NO3, p_PHY(i)%C, &
                                p_PHY(i)%N, p_PHY(i)%Chl, p_PHY(i)%Cdiv, dC_, dN_, dChl_)
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
         ENDSELECT

         uptake   =   uptake + dN_ * p_PHY(i)%num 
         NPPc_(k) = NPPc_(k) + dC_ * p_PHY(i)%num *1d-9/Hz(k)*12.d0*dtdays !Unit: mgC m-3 d-1

         ! Save carbon-specific growth rate
         p_PHY(i)%mu_C = dC_/p_PHY(i)%C

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

         !Save fitness of each particle (per day)
         p_PHY(i)%fitness = (p_PHY(i)%C*p_PHY(i)%num*1d-9/Hz(k) - BC(j))/BC(j)/dtdays

      enddo !End of looping through particles
  endif !End if of N_ > 0

  uptake = uptake*1d-9/Hz(k) !Convert uptake to mmol N m-3
  Pmort  =  Pmort*1d-9/Hz(k) !Convert Pmort to mmol N m-3

  !Now calculate NO3 and DET
  t(iNO3,k) = NO3 + dtdays*(pp_ND + RES - uptake)

  Varout(iNO3, k) = t(iNO3, k)

  t(iDET,k) = DET + Pmort + dtdays*(pp_DZ - pp_ND)

  Varout(iDET,k) = t(iDET,k)

  if (allocated(index_))  deallocate(index_)
  if (allocated(Pmatrix)) deallocate(Pmatrix)
  if (allocated(BN))      deallocate(BN)
  if (allocated(BC))      deallocate(BC)
ENDDO !End of looping across different depths

!update t(iPN, :), t(iPC,:), and t(iChl,:) and compute mean trait and trait variance
call Par2PHY

END SUBROUTINE BIOLOGY