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
integer            :: j

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

   !Calculate vertical Kv
   call analytic_Kv(nlev, Kv0(1), Kvmax(1), 1d-5, MLD(1), Kv)

   !Start biology
   !Update environmental variables associated with each particle
   call UPDATE_PARTICLE_FORCING

   !Update Eulerian state variables and phytoplankton particles
   call BIOLOGY

   !Save output
   IF (mod(it, nsave) == 1) THEN

      !Save the model output of particles to a separate file
      if (current_hour == 0) then

         !Name the par_file
         write(par_file, 100) 'ParY', current_year, '_D', current_DOY
         call create_Particle_file(par_file)
         call       save_particles(par_file)
      endif
   ELSE
      call save_Eulerian
      call save_particles(par_file)
   ENDIF

   ! Pass the new state variables to Vars
   do j = 1,NVAR
      t(j,:) = Varout(j,:)
   enddo
 
   ! Vertical random walk for particles
   Do j = 1, N_par
      !Assume closed boundary for particles
      CALL LAGRANGE(nlev, Z_w, Kv, w, p_PHY(j)%iz, p_PHY(j)%rz)
   Enddo

   ! Diffusion for Eulerian fields except phytoplankton
   do j = 1,NVAR
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
END SUBROUTINE TIMESTEP

SUBROUTINE UPDATE_PARTICLE_FORCING
use forcing,         only : PARw, Temp
use state_variables, only : p_PHY, t, iNO3
use grid,            only : Z_w, Z_r, nlev
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
   p_PHY(k)%NO3 = zout(1,1)
enddo

DEALLOCATE(dat_in)
DEALLOCATE(zout)
DEALLOCATE(zp)
END SUBROUTINE UPDATE_PARTICLE_FORCING