SUBROUTINE TIMESTEP
USE MPI_Setting
use forcing
use Time_setting
use state_variables
use grid
use IO
use Trait_functions,  only : PHY_C2Vol
implicit none

character(LEN=20)  :: par_file   = 'Phyto'
character(LEN=20)  :: passive_file = 'Passive'
real,       parameter  :: cnpar      = 0.6d0
real,       parameter  :: Taur(nlev) = 1D12  !Relaxation time
real,       parameter  :: zero       = 0.d0  !Vectors of zero
real,       parameter  :: Vec0(nlev) = zero  !Vectors of zero
real :: par_save_freq = 0d0           !scratch variable for saving frequency of particles
integer,parameter  :: mode0      = 0
integer,parameter  :: mode1      = 1
integer :: i,j

!Benchmarking
real(4) :: dt1, dt2, dt3, dt4, dt5, t1, t2, t3, t4, t5, t6

! 'START TIME STEPPING'
dt1 = 0.d0
dt2 = 0.d0
dt3 = 0.d0
dt4 = 0.d0
dt5 = 0.d0

DO it = 1, Nstep+1

  IF (TASKID .EQ. 0) THEN
    call cpu_time(t1) 
    call update_time

    !For each time step, read in external environmental data%%%%%%%%%%%%%%%
    !Interpolate vertical profile of temperature at each timestep
    call time_interp(int(current_sec), N_time_temp, nlev, obs_time_temp, VTemp, Temp)

    !Calculate PAR
    call VERTICAL_LIGHT(current_DOY, sec_of_day, t(iChl,:))

    !Directly use the ROMS model output
    call time_interp(int(current_sec), N_time_Kv, nlev+1, obs_time_Kv, VKv, Kv)

    !Start biology
    !Update environmental variables associated with each particle
    call UPDATE_PARTICLE_FORCING

    call cpu_time(t2) 
    dt1 = t2 - t1 + dt1 !The time for interpolating environmental data

    call BIOLOGY
    call cpu_time(t3) 
    dt2 = t3 - t2 + dt2 !The time for biology

    !Save the Eulerian output every day
    IF (mod(it, nsave) == 1) THEN

       ! Add calculations of total nitrogen and save to Eulerian output files
       call Cal_total_N !Including randomly split cells

       write(6, 101) "Day", current_day, ": Total Nitrogen =", Ntot

       call save_Eulerian
       call save_Kv

    ENDIF

    !Save the model output of particles to a separate file every day
    !And save the particles every hour
    If (current_year < NDay_Run/365) then
      par_save_freq = d_per_s
    Else
     par_save_freq = s_per_h
    Endif

    If (mod(current_sec, par_save_freq) == 0) then  !Can be modified to save the particles at daily frequency
        if (current_hour == 0) then

           !Create the phyto. particle file
           write(par_file, 100) 'ParY', current_year, '_D', current_DOY
           call create_Particle_file(par_file, phyto)

           !Create the passive particle file
           write(passive_file, 102) 'PassY', current_year, '_D', current_DOY
           call create_Particle_file(passive_file, passive)

        endif

        call save_particles(par_file, phyto)
        call save_particles(passive_file, passive)

        !Save N_birth, N_death, and N_mutate into a single file every hour
        call save_particles(Death_file, death)
    Endif

    call cpu_time(t4) 
    dt3 = t4 - t3 + dt3 !The time for saving data
  ENDIF  !! <-- End of main thread

  ! Vertical random walk for both phyto. particles (that are not dead) and passive particles
  !Use openmpi to allocate particles to ntasks threads
  call LAGRANGE

  IF (TASKID .EQ. 0) THEN
    call cpu_time(t5) 
    dt4 = t5 - t4 + dt4  !The time for random walk

    ! Update 
    ! Diffusion for Eulerian fields except phytoplankton
    Do j = 1,NVAR
       ! Zero flux at bottom
       if (j .ne. iPC .and. j .ne. iPN .and. j .ne. iCHL) then
          call diff_center(nlev,dtsec,cnpar,1,Hz, Neumann, Neumann, &
                        zero, zero, Kv, Vec0,Vec0,Taur, t(j,:), t(j,:), t(j,:))
       endif
    Enddo

    ! Sinking:
    do j = 1,NVsinkterms
       SELECT CASE (bot_bound)
       case(Neumann)  ! closed at bottom (Conserve total N mass)
          call adv_center(nlev,dtsec,Hz,Hz,ww(:,j),1,1,zero,zero,    6,mode1,t(Windex(j),:))
       case(Dirichlet)
          ! Open bottom boundary
          call adv_center(nlev,dtsec,Hz,Hz,ww(:,j),1,2,zero,t(Windex(j),1),6,mode1,t(Windex(j),:))
       case default
          stop "The boundary conditions incorrect! STOP!"
       ENDSELECT
    enddo

    call cpu_time(t6) 
    dt5 = t6 - t5 + dt5  !The time for diffusion and sinking
  ENDIF 

ENDDO

if (taskid == 0) then
  print '("Environmental interpolation costs ",f8.3," hours.")', dt1/3600.0 
  print '("Biology costs ",f8.3," hours.")', dt2/3600.0 
  print '("Saving data costs ",f8.3," hours.")', dt3/3600.0 
  print '("Random walk costs ",f8.3," hours.")', dt4/3600.0 
  print '("Diffusion and detritus sinking cost ",f8.3," hours.")', dt5/3600.0 
endif

100 format(A4,I0,A2,I0)
101 format(A3,1x, I0, 1x, A25, 1x, F10.4)
102 format(A5,I0,A2,I0)
END SUBROUTINE TIMESTEP

SUBROUTINE UPDATE_PARTICLE_FORCING
use forcing,         only : PARw, Temp
use grid,              only : Z_w, Z_r, nlev
use state_variables, only : p_PHY, t, iNO3
implicit none

! interpolation results for particles
real, ALLOCATABLE :: zp(:), zout(:,:), dat_in(:,:)
integer   :: AllocateStatus = 0
integer   :: k              = 0
!==================================================================

ALLOCATE(zp(1), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating zp ***"
zp(1) = 0.
ALLOCATE(zout(1,1), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating zout ***"
zout(1,1) = 0.

ALLOCATE(dat_in(0:nlev, 1), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating dat_in ***"
dat_in(:,1) = PARw(:)  !Observed profile of PAR

do k = 1, size(p_PHY) 
   zp(1) = p_PHY(k)%rz
   call gridinterpol(nlev+1, 1, Z_w, dat_in, 1, zp, zout)
   p_PHY(k)%PAR = zout(1,1)
enddo
DEALLOCATE(dat_in)

ALLOCATE(dat_in(nlev, 1), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating dat_in ***"
dat_in(:,1) = Temp(:)  !Observed profile of temperature
do k = 1, size(p_PHY) 
   zp(1) = p_PHY(k)%rz
   call gridinterpol(nlev, 1, Z_r, dat_in, 1, zp, zout)
   p_PHY(k)%temp = zout(1,1)
enddo

dat_in(:,1) = t(iNO3,:)  !Observed profile of NO3

do k = 1, size(p_PHY) 
   zp(1) = p_PHY(k)%rz
   call gridinterpol(nlev, 1, Z_r, dat_in, 1, zp, zout)

   !Needs to set an upper limit to the NO3 concentration of the particle which cannot be greater than the current concentration of the grid
   p_PHY(k)%NO3 = min(zout(1,1), t(iNO3, p_PHY(k)%iz))

enddo

DEALLOCATE(dat_in)
DEALLOCATE(zout)
DEALLOCATE(zp)
END SUBROUTINE UPDATE_PARTICLE_FORCING

subroutine Cal_total_N
use grid, only : nlev, Hz, Z_r
use state_variables, only : t, Ntot, iPC, iCHL,iPN, iZOO, iNO3, iDET,N_PAR, p_PHY, IDmax, NZOO
implicit none
integer :: k,i,m
real      :: rnd
real      :: max_num = 0.d0

Ntot = 0d0
do k = 1, nlev
   if (t(iPC,k) .ne. t(iPC,k)) then
      write(6,*) "Phyto Carbon is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iPC,k) < 0d0) then
      write(6,*) "Phyto Carbon is negative at depth", Z_r(k)
      stop 
   endif

   if (t(iCHL,k) .ne. t(iCHL,k)) then
      write(6,*) "Chl is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iCHL,k) < 0d0) then
      write(6,*) "Chl is negative at depth", Z_r(k)
      stop 
   endif

   if (t(iPN,k) .ne. t(iPN,k)) then
      write(6,*) "Phyto N is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iPN,k) < 0d0) then
      write(6,*) "Phyto N is negative at depth", Z_r(k)
      stop 
   endif

   do m = 1, NZOO
     if (t(iZOO(m),k) .ne. t(iZOO(m),k)) then
        write(6,*) "ZOO", m, " is NaN at depth", Z_r(k)
        stop 
     endif

     if (t(iZOO(m),k) < 0d0) then
        write(6,*) "ZOO", m, " is negative at depth", Z_r(k)
        stop 
     endif
   enddo

   if (t(iDET,k) .ne. t(iDET,k)) then
      write(6,*) "DET is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iDET,k) < 0d0) then
      write(6,*) "DET is negative at depth", Z_r(k)
      stop 
   endif

   if (t(iNO3,k) .ne. t(iNO3,k)) then
      write(6,*) "NO3 is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iNO3,k) < 0d0) then
      write(6,*) "NO3 is negative at depth", Z_r(k)
      stop 
   endif

   !The code below computes the maximal number of cells per super-individual
   max_num = 0.d0
   DO i = 1, N_PAR
      max_num = max(max_num, p_PHY(i)%num)
   END DO

   !The following code tries to find dead superindividuals and split
   !If there is a dead superindividual, choose a random live superindividual 
   !with a good number of cells (greater than 1/10 of maximal numbers of cells 
   !per super-individual) and split it
   DO i = 1, N_PAR
      IF (.not. p_PHY(i)%alive) THEN
         m = i
         do while ((.not. p_PHY(m)%alive) .or. (p_PHY(m)%num .le. max_num/1d1))
            call random_number(rnd)
            m = max(int(dble(N_PAR-1)*rnd + 1.d0), 1)
         enddo

         !Split it in to two identical superindividuals
         !The first one is identical with the parent, except that the number of cells is halved
         p_PHY(m)%num = p_PHY(m)%num/2d0

         !The second one is identical with its twin, but its ID needs to change to a new number (not overlapping with any ID of current superindividuals)
         p_PHY(i)        = p_PHY(m)
         p_PHY(i)%ID     = IDmax+1
         IDmax           = IDmax+1 !Update maximal ID
      ENDIF
   ENDDO

   !Update total N
   Ntot = Ntot + Hz(k)*(t(iPN,k)  + t(iNO3, k) + t(iDET, k))

   !Add total ZOOplankton N into total N
   do i = 1, NZOO
      Ntot = Ntot + t(iZOO(i), k) * Hz(k)
   enddo
enddo
end subroutine Cal_total_N
