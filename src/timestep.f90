SUBROUTINE TIMESTEP
use forcing
use Time_setting
use state_variables
use grid
use IO
implicit none

character(LEN=20)  :: par_file   = 'Particles'
real,   parameter  :: cnpar      = 0.6d0
real,   parameter  :: Taur(nlev) = 1D12  !Relaxation time
real,   parameter  :: zero       = 0.d0  !Vectors of zero
real,   parameter  :: Vec0(nlev) = zero  !Vectors of zero
integer,parameter  :: mode0      = 0
integer,parameter  :: mode1      = 1
integer            :: j, iit

! 'START TIME STEPPING'
DO it = 1, Nstep+1

   call update_time

   !For each time step, read in external environmental data%%%%%%%%%%%%%%%
   !Interpolate vertical profile of temperature at each timestep
   call time_interp(int(current_sec), N_time_temp, nlev, obs_time_temp, VTemp, Temp)

   !Calculate PAR
   call VERTICAL_LIGHT(current_DOY, sec_of_day, t(iChl,:))

   !Interpolate MLD, KV0, KVmax
   call time_interp(int(current_sec), N_time_Kv, 1, obs_time_Kv, obs_MLD,   MLD)
   call time_interp(int(current_sec), N_time_Kv, 1, obs_time_Kv, obs_Kv0,   Kv0)
   call time_interp(int(current_sec), N_time_Kv, 1, obs_time_Kv, obs_Kvmax, Kvmax)

   !Calculate vertical Kv (Kbg follows Christina L&O 2005)
   call analytic_Kv(nlev, Kv0(1), Kvmax(1), 3d-5, MLD(1), Kv)

   !Start biology
   !Update environmental variables associated with each particle
   call UPDATE_PARTICLE_FORCING

   call BIOLOGY

   !Save the Eulerian output every day
   IF (mod(it, nsave) == 1) THEN

      ! Add calculations of total nitrogen and save to Eulerian output files
      call Cal_total_N
      write(6, 101) "Day", current_day, ": Total Nitrogen =", Ntot
      call save_Eulerian
      call save_Kv

   ENDIF

   !Save the model output of particles to a separate file every day
   !And save the particles every hour
   If (mod(current_sec, s_per_h) == 0) then
       if (current_hour == 0) then

          !Name the par_file
          write(par_file, 100) 'ParY', current_year, '_D', current_DOY
          call create_Particle_file(par_file)
       endif

       call save_particles(par_file)
   Endif

   ! Vertical random walk for particles that are not dead
   DO iit = 1, Nrand
     Do j = 1, N_par
        !Assume closed boundary for particles
        if (p_PHY(j)%alive) then
           CALL LAGRANGE(nlev, Z_w, Kv, w, p_PHY(j)%iz, p_PHY(j)%rz)
        endif
     Enddo
   ENDDO

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

ENDDO

100 format(A4,I0,A2,I0)
101 format(A3,1x, I0, 1x, A25, 1x, F10.4)
END SUBROUTINE TIMESTEP

SUBROUTINE UPDATE_PARTICLE_FORCING
use forcing,         only : PARw, Temp
use grid,              only : Z_w, Z_r, nlev
use state_variables, only : p_PHY, t, iNO3
implicit none

! interpolation results for particles
real, ALLOCATABLE :: zp(:), zout(:,:), dat_in(:,:)
integer           :: AllocateStatus = 0
integer           :: k              = 0
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

   !The following code tries to find dead superindividuals and split
   !If there is a dead superindividual, choose a random live superindividual and split it
   DO i = 1, N_PAR
      IF (.not. p_PHY(i)%alive) THEN
         m = i
         do while (.not. p_PHY(m)%alive)
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

   Ntot = Ntot + Hz(k)*(t(iPN,k)  + t(iNO3, k) + t(iDET, k))

   do i = 1, NZOO
      Ntot = Ntot + t(iZOO(i), k) * Hz(k)
   enddo
enddo
end subroutine Cal_total_N
