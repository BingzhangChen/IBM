SUBROUTINE INITIALIZE
USE params
USE state_variables
USE Time_setting, only: d_per_s, dtsec, dtdays, Nstep, nsave,NDay_Run,Nrand
USE grid,                   only: Z_w, hmax
USE Trait_functions,  only: PHY_ESD2C
USE IO
USE forcing
USE MPI_Setting
IMPLICIT NONE

integer    :: rc    = 0
integer    :: k     = 0
integer    :: j_     = 0
real          :: cff      = 0.d0
real          :: Z_avg = 0.d0
integer, parameter :: namlst   = 20   !Unit time for namelist files
character(len=10)  :: par_file = 'Filename'
character(len=12)  :: passive_file = 'Passive'
character(len=10), parameter :: format_string = "(A3,I0)"
integer           :: AllocateStatus = 0
!==========================================================

!Namelist definition of time settings
namelist /timelist/  NDay_Run, dtsec, nsave, Nrand

!Namelist definition of model choice and parameters
namelist /paramlist/ Model_ID, N_pass, N_Par, mu0, aI0, KN, gmax, Kp, mz,GGE, unass, RDN, &
                     wDET, SDZoo, nu,sigma

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
Nstep  = NDay_Run*INT(d_per_s)/INT(dtsec) 

!Calculate dtdays
dtdays = dtsec/d_per_s
if (taskid==0) write(6,'(A13,1x,1pe12.2,A12)') 'Timestepping: ', dtdays, 'of one day.'

!==========================================================
!Read parameter namelist
!Check whether the namelist file exists.
inquire (file='param.nml', iostat=rc)

if (rc /= 0) then
    write (6, '(a)') 'Error: namelist file Model.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='param.nml',status='old',action='read')
read(namlst,nml=paramlist)
close(namlst)

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
  write(6,'(A15,1x,I1)') 'Select Model ID', Model_ID
  
  !Initialize state variables
  ALLOCATE(p_pass(N_Pass), stat=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_pass ***"
  
  ALLOCATE(p_PHY(N_Par), stat=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_PHY ***"
 
  ! Initialize initial NO3:
  t(iNO3,:) = 0.5d0
  
  do k = 1, NZOO
     t(iZOO(k),:) = 0.1d0/dble(NZOO) !Assuming an initial condition of uniform biomass among different zoo. size classes
  
     !Initialize zooplankton size (logESD)
     ESDZOO(k) = MinSzoo + dble(k-1)*dZOOESD
  
     !Make the model write out the zooplankton size
     write(6,1001) "ZOO", k, "ESD = ", exp(ESDZOO(k)), "micron"
  
     !Compute volume of zooplankton
     VolZOO(k) = pi/6d0*exp(ESDZOO(k))**3
  enddo
  
  !Following Verity et al. AME (1996)
  t(iDET,:) = .1d0
  
  !For lagrangian model, initialize individual particles
  
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
  
  !Print out the number of phyto particles
  write(6, 1004)  N_Par

  !Initialize phytoplankton super-individuals
  DO k = 1, N_PAR
     p_PHY(k)%ID = k        
     p_PHY(k)%alive = .true.
  
     !Initialize phytoplankton optimal temperature (Topt) from a uniform distribution between 2 and 30 degree celcius
     call random_number(cff)
     cff = 2. + cff * (30. - 2.)
     p_PHY(k)%Topt = cff
  
     !Initialize phytoplankton size from a uniform distribution between 0.8 and 60 um
     call random_number(cff)
     cff = log(0.8d0) + cff * (log(60.d0) - log(0.8d0))
     cff = exp(cff) !ESD
  
     p_PHY(k)%C    = PHY_ESD2C(cff)      !Unit: pmol C cell-1
     p_PHY(k)%N    = p_PHY(k)%C/106.*16. !Unit: pmol N cell-1
     p_PHY(k)%Chl  = p_PHY(k)%C* 12./50. !Unit: pgChl cell-1
     p_PHY(k)%CDiv= p_PHY(k)%C* 2d0 !Unit: pmol C cell-1
  
     !Initialize log(alphaChl) from a uniform distribution from log(0.01) to log(0.5) (W m-2)-1 (gChl molC)-1 d-1
     call random_number(cff)
     cff = log(0.01d0) + cff * (log(0.5) - log(0.01d0))
     p_PHY(k)%LnalphaChl = cff
  
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
     ENDDO
  ENDDO
  
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

  IDmax=N_par
  
  call UPDATE_PHYTO  !Initialize Eulerian concentrations of phytoplankton carbon, nitrogen and Chl
  
  !Initialiaze labels for model output
  Labelout(iNO3) = 'NO3'
  Labelout(iDET) = 'DET'
  Labelout(iPC)  = 'PC '
  Labelout(iPN)  = 'PN '
  Labelout(iCHL) = 'CHL'
  Labelout(oNPP) = 'NPP'
  Labelout(oTEMP)= 'TEMP'
  Labelout(oPAR) = 'PAR'
  
  do k = 1, NZOO
     write(Labelout(iZOO(k)),  format_string) 'ZOO', k
     write(Labelout( oFZ(k)),  format_string) 'FZO',k
  enddo
  
  write(6, '(a)') 'Write out the labels for output for validation:'
  do k = 1, Nout
     write(6, 1000) 'Labelout(',k,') = ',trim(Labelout(k))
  enddo

  !Initialize Varout
  do k = 1, NVAR
     Varout(k,:) = t(k,:)
  enddo
  
  !Save initial state to external file
  call create_Eulerian_file
  call save_Eulerian
  
  !Save Kv
  call create_Kv_file
  call save_Kv
  
  !Save N_birth, N_death, and N_mutate
  call create_Particle_file(Death_file, Death)
  call       save_particles(Death_file, death)
  
  !Name the initial phyto. particle file
  par_file = 'ParY1_D0'
  
  call create_Particle_file(par_file, phyto)
  call         save_particles(par_file, phyto)
  
  !Name the initial passive particle file
  passive_file = 'PassY1_D0'
  
  call create_Particle_file(passive_file, passive)
  call         save_particles(passive_file, passive)
  
  ! Prepare forcing
  ! Prepare temperature for temporal interpolation
  call extract_WOAtemp
  
  ! Prepare Kv
  call extract_Kv
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

1000 format(A9, I0, A4, A5)
1001 format(A3, 1x, I0, 1x, A6, F8.2, 1x, A6 )
1002 format('Number of passive particles is:', 1x, I0)
1003 format('Their average position is:', 1x, F12.2, 1x, 'm')
1004 format('Number of phyto. particles is:', 1x, I0)
END SUBROUTINE INITIALIZE
