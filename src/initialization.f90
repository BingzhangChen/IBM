SUBROUTINE INITIALIZE
USE params
USE state_variables
USE Time_setting
USE grid,             only: Z_w, hmax
USE Trait_functions,  only: PHY_ESD2C, PHY_C2Vol
USE NETCDF_IO
USE forcing
USE MPI_Setting
IMPLICIT NONE

integer    :: rc    = 0
integer    :: k     = 0
integer    :: j_    = 0
real       :: cff   = 0.d0
real       :: Z_avg = 0.d0
real       :: Vol   = 0.
real       :: QN    = 0.
integer, parameter :: namlst         = 20   !Unit time for namelist files
integer            :: AllocateStatus = 0
logical            :: exists         = .true.
!==========================================================

!Namelist definition of time settings
namelist /timelist/  NDay_Run, dtsec, nsave, Nrand, read_previous_output

!Namelist definition of model choice and parameters
namelist /paramlist/ Model_ID, N_pass, N_Par, mu0, aI0, KN, gmax, Kp, mz,GGE, unass, RDN, &
                     wDET, SDZoo

! Check whether the namelist file exists.
inquire (file='time.nml', iostat=rc)

if (rc /= 0) then
    write (6, '(a)') 'Error: namelist file time.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='time.nml',status='old',action='read')
read(namlst,nml=timelist)
close(namlst)

!Total Number of time steps
Nstep = NDay_Run * INT(d_per_s)/INT(dtsec) 

!Calculate dtdays
dtdays = dtsec/d_per_s
if (taskid==0) then
   write(6,'(A13,1x,1pe12.2,A12)') 'Timestepping: ', dtdays, 'of one day.'
   write(6,'(A30,1x,I0)') 'Total number of simulation days: ', NDay_Run
endif
!==========================================================
!Read parameter namelist
!Check whether the namelist file exists.
inquire (file='param.nml', exist = exists)

if (.not. exists) then
    write (6, '(a)') 'Error: namelist file Model.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='param.nml',status='old',action='read')
read(namlst,nml=paramlist)
close(namlst)

IF (taskid .EQ. 0) WRITE(6,'(A15,1x,I1)') 'Select Model ID', Model_ID

! Check if need to read previous model output
IF (read_previous_output .EQ. 1) THEN
  inquire (file=restart_fname, exist = exists)
  
  if (.not. exists) then
    write (6, '(a)') 'Error: restart.nc does not exist.'
    stop
  else
    !Read restart.nc
    call read_restart
    write(6, '(a)') 'Restart file read successfully!'

    !Update time
    call update_time
  end if
ELSE
  restart_step = 1
ENDIF

!Allocate vectors for MPI data transfer
if (mod(N_Par + N_Pass, numtasks) .ne. 0) stop "Number of CPUs incorrect!"
N_chunk = (N_Par + N_Pass)/numtasks

ALLOCATE(Zr_mpi(N_PAR + N_Pass), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating Zr_mpi ***"
Zr_mpi(:) = 0.d0

ALLOCATE(Zi_mpi(N_PAR + N_Pass), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating Zi_mpi ***"
Zi_mpi(:) = 0

ALLOCATE(C_mpi(N_Par + N_Pass), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating C_mpi ***"
C_mpi(:) = 0d0

ALLOCATE(tag6(numtasks-1), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating tag6 ***"

ALLOCATE(tag7(numtasks-1), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating tag7 ***"

do k = 1, numtasks-1
   tag6(k) = 5 + k
   tag7(k) =15 + k
enddo

IF (TASKID==0) THEN

  ! Prepare forcing
  ! Prepare temperature for temporal interpolation
  call extract_WOAtemp
  
  ! Prepare Kv
  call extract_Kv

  !Initialize ZOOplankton size
  If (Model_ID .eq. GMK98_Size          .or. Model_ID .eq. GMK98_ToptSize .or.  &
      Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_SizeLight) then

    do k = 1, NZOO
    
       !Initialize zooplankton size (logESD)
       ESDZOO(k) = MinSzoo + dble(k-1)*dZOOESD
    
       !Make the model write out the zooplankton size
       write(6,1001) "ZOO", k, "ESD = ", exp(ESDZOO(k)), "micron"
    
       !Compute volume of zooplankton
       VolZOO(k) = pi/6d0*exp(ESDZOO(k))**3
    enddo
  Else

    !Check if NZOO is consistent with Model_ID (i.e., without size classes, NZOO should be one)
    if(NZOO > 1) then
      stop "Number of zooplankton size classes should be ONE if size is not modelled!"
    endif

    VolZOO(1) = pi/6d0*20.0**3  !Assume zooplankton volume 20 micron
  Endif
 
  IF (read_previous_output .EQ. 1) THEN

    !Check whether Euler.nc exists
    inquire (file=Euler_FNAME, exist = exists)
    
    if (exists) then
        write (6, '(a)') 'Warning: Euler.nc already exists and will be overwritten!'
        stop
    end if

  ELSE
 
    ! Initialize initial NO3:
    t(iNO3,:) = 0.5d0
  
    If (Model_ID .eq. GMK98_Size          .or. Model_ID .eq. GMK98_ToptSize .or.  &
        Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_SizeLight) then

      do k = 1, NZOO
         t(iZOO(k),:) = 0.1d0/dble(NZOO) !Assuming an initial condition of uniform biomass among different zoo. size classes
      enddo
    Else

      !Check if NZOO is consistent with Model_ID (i.e., without size classes, NZOO should be one)
      if(NZOO > 1) then
        stop "Number of zooplankton size classes should be ONE if size is not modelled!"
      endif

      t(iZOO(1),:) = 0.1d0
    Endif
  
    !Following Verity et al. AME (1996)
    t(iDET,:) = .1d0
  
    !For lagrangian model, initialize individual particles
  
    IF (.not. allocated(p_pass)) ALLOCATE(p_pass(N_Pass), stat=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_pass ***"
    
    IF (.not. allocated(p_PHY)) ALLOCATE(p_PHY(N_PAR), stat=AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_PHY ***"

    !Initialize passive particles
    Z_avg = 0d0
    DO k = 1, N_Pass
       p_pass(k)%ID = k
  
       !Compute rx (randomly distributed between Z_w(0) and Z_w(nlev)
       CALL RANDOM_NUMBER(cff)
       p_pass(k)%rz = cff *(Z_w(nlev)-Z_w(0)) + Z_w(0)
  
       !Update average position
       Z_avg = Z_avg + p_pass(k)%rz
  
       ! Find iz
       DO j_ = 1, nlev
          IF (Z_w(j_-1) .le. p_pass(k)%rz .AND. Z_w(j_) .gt. p_pass(k)%rz) THEN
             p_pass(k)%iz = j_
             EXIT            
          ENDIF
       ENDDO
    ENDDO
  
    !Print out the number of passive particles
    write(6, 1002)  N_Pass
    write(6, 1003)  Z_avg/dble(N_Pass)
    !!End of initializing passive particles
  
    !Print out the number of phyto particles
    write(6, 1004)  N_Par

    !Initialize phytoplankton super-individuals
    DO k = 1, N_PAR
       p_PHY(k)%ID = k        
       p_PHY(k)%alive = .true.
  
       if (Model_ID .eq. GMK98_Topt .or. Model_ID .eq. GMK98_ToptLight .or. &
           Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_ToptSize) then
          !Initialize phytoplankton optimal temperature (Topt) from a uniform distribution between 2 and 30 degree celcius
          call random_number(cff)
          cff = 2. + cff * (30. - 2.)
          p_PHY(k)%Topt = cff
       endif

       If (Model_ID .eq. GMK98_Size .or. Model_ID .eq. GMK98_ToptSize .or.  &
           Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_SizeLight) then

           !Initialize phytoplankton size from a uniform distribution between 0.8 and 60 um
           call random_number(cff)
           cff = log(0.8d0) + cff * (log(60.d0) - log(0.8d0))
           cff = exp(cff) !ESD
           p_PHY(k)%C = PHY_ESD2C(cff)      !Unit: pmol C cell-1
           p_PHY(k)%CDiv = p_PHY(k)%C*2d0   !Maximal size, Unit: pmol C cell-1

           !Assume nutrient replete for initial condition
           Vol= PHY_C2Vol(p_PHY(k)%CDiv)
           QN = QNmax_a * Vol**QNmax_b
           p_PHY(k)%N = p_PHY(k)%C*QN !Unit: pmol N cell-1

       else
           p_PHY(k)%C = PHY_ESD2C(1.0)      !Unit: pmol C cell-1
           p_PHY(k)%N = p_PHY(k)%C/106.*16. !Unit: pmol N cell-1
           p_PHY(k)%CDiv = p_PHY(k)%C* 2d0     !Unit: pmol C cell-1
       endif

       p_PHY(k)%Chl  = p_PHY(k)%C* 12./50. !Unit: pgChl cell-1

       if (Model_ID .eq. GMK98_Light .or. Model_ID .eq. GMK98_ToptLight .or. &
           Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_SizeLight) then
 
           !Initialize log(alphaChl) from a uniform distribution from log(0.01) to log(0.5) (W m-2)-1 (gChl molC)-1 d-1
           call random_number(cff)
           cff = log(0.01d0) + cff * (log(0.5) - log(0.01d0))
           p_PHY(k)%LnalphaChl = cff
        endif
  
       !Initialize the number of cells associated with each super-individual (assuming initial phytoplankton nitrogen is 0.1 mmol m-3)
       p_PHY(k)%num = 0.1 * hmax/dble(N_PAR)/p_PHY(k)%N * 1d9
  
       !Compute rx (randomly distributed between Z_w(0) and Z_w(nlev)
       CALL RANDOM_NUMBER(cff)
       p_PHY(k)%rz = cff *(Z_w(nlev)-Z_w(0)) + Z_w(0)
  
       ! Find iz
       DO j_ = 1, nlev
          IF (Z_w(j_-1) .le. p_PHY(k)%rz .AND. Z_w(j_) .gt. p_PHY(k)%rz) THEN
             p_PHY(k)%iz = j_
             EXIT            
          ENDIF
       ENDDO  !End of finding iz for each particle
    ENDDO !End of initializing PHY particles
  ENDIF !==> End of if statement of reading restart.nc
  
  !Pass positions of passive and phyto. particles to Zi_mpi, Zr_mpi and C_mpi
  do k = 1, N_par + N_Pass
    if (k .le. N_Pass) then !Passive particles
      Zi_mpi(k)=p_pass(k)%iz 
      Zr_mpi(k)=p_pass(k)%rz 
    else
      j_ = k - N_pass  !The index for phytoplankton particles
      Zi_mpi(k)=p_PHY(j_)%iz 
      Zr_mpi(k)=p_PHY(j_)%rz 
      C_mpi(j_)=p_PHY(j_)%C
    endif
  enddo

  !Estimating IDmax
  IDmax = N_par

  IF (read_previous_output .eq. 1) THEN
    do k = 1, N_PAR
      IDmax = MAX(IDmax, p_PHY(k)%ID)
    enddo
  ENDIF

  call UPDATE_PHYTO  !Initialize Eulerian concentrations of phytoplankton carbon, nitrogen and Chl
  call UPDATE_PARTICLE_FORCING
  
  !Initialize Varout
  do k = 1, NVAR
     Varout(k,:) = t(k,:)
  enddo
  
  !Initialize time
  it = restart_step
  call update_time

  !Save initial state to external file
  call create_Eulerian_file
  irec_Euler = 1
  call write_Eulerfile(irec_Euler, current_day, current_hour)

  
  !Name the initial phyto. particle file
  write(par_file, 1005) 'ParY', current_year, '.nc'
  irec_PHY = 1
  
  call Create_PHY_particlefile(par_file)
  call write_PHY_particlefile(par_file, irec_PHY, current_DOY, current_hour)
  
  !Name the initial passive particle file
  write(passive_file, 1006) 'PassY', current_year, '.nc'
  irec_Pass = 1
  
  call Create_Pass_particlefile(passive_file)
  call write_Pass_particlefile(passive_file, irec_Pass, current_DOY, current_hour)
 
ENDIF   !End of parent process (taskid == 0)

!Sinking rate
!Initialize sinking rate (UNIT: m/s !):
ww(:,:) = 0d0
do k = 0,nlev-1
   !Phytoplankton no sinking 
   !Detritus sinking rate (convert to UNIT: m/s)
   ww(k,NVsinkterms) = -wDET/dble(d_per_s) 
enddo
return

1001 format(A3, 1x, I0, 1x, A6, F8.2, 1x, A6 )
1002 format('Number of passive particles is:', 1x, I0)
1003 format('Their average position is:', 1x, F12.2, 1x, 'm')
1004 format('Number of phyto. particles is:', 1x, I0)
1005 format(A4,I0,A3)
1006 format(A5,I0,A3)
END SUBROUTINE INITIALIZE
