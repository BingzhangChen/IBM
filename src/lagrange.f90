! Phytoplankton and passive particles random walk
SUBROUTINE LAGRANGE
USE MPI_Setting
USE Time_setting, only: dtsec, Nrand
USE grid,         only: nlev, Z_w, Hz
USE forcing,    only: Kv, w, Kbg
USE STATE_VARIABLES, only: p_PHY, p_pass, N_PAR, N_pass, Zi_mpi, C_mpi, Zr_mpi
USE Trait_functions, only: PHY_C2Vol
IMPLICIT NONE

REAL     :: bio_w(0:nlev) = 0d0  !Scratch sinking rate in 1D

! !LOCAL VARIABLES:
INTEGER         :: i,j, zi, iit
integer           :: m = 0
REAL               :: rnd, rnd_var_inv, dt, step1, step2, zloc, zp
REAL               :: dzn(nlev),step
REAL               :: visc,rat,dt_inv, vs
REAL               :: wint  !interpolated w
REAL               :: dKvdz(nlev) = 0.d0  !Gradient of Kv
REAL,parameter :: rnd_var=0.333333333
REAL, external :: phyto_sinking
integer, parameter  :: tag1 = 1
integer, parameter  :: tag2 = 2
integer, parameter  :: tag3 = 3
integer, parameter  :: tag4 = 4
integer, parameter  :: tag5 = 5
integer :: stat(MPI_STATUS_SIZE)   ! required variable for receive routines
!--------------End of declaration-------------------------------------------------------

!The parent thread passes information to child processes

! Synchronize all processes:
call MPI_BARRIER (MPI_COMM_WORLD,ierr)

IF (TASKID .EQ. 0) THEN
   !Allocate p_Pass and p_PHY to Zi_Mpi and Zr_MPI in the parent process
   do m = 1, N_Pass + N_PAR
      if (m .le. N_Pass) then
          Zr_MPI(m) = p_Pass(m)%rZ
          Zi_MPI(m) = p_Pass(m)%iz
      else
          Zr_MPI(m) = p_PHY(m-N_Pass)%rZ
          Zi_MPI(m) = p_PHY(m-N_Pass)%iZ
          C_MPI(m) =  p_PHY(m-N_Pass)%C
      endif
   enddo

   Do j = 1, numtasks-1

      ! Send  Kv to child processes:
      call MPI_SEND(Kv, 1+nlev, MPI_REAL8, j, tag4, MPI_COMM_WORLD, ierr)

      ! Send  w to child processes:
      call MPI_SEND(w, 1+nlev, MPI_REAL8, j, tag5, MPI_COMM_WORLD, ierr)

      !Send Zi_Mpi, Zr_MPI and C_MPI from parent to child processes
      call MPI_SEND(Zi_MPI, N_Pass + N_PAR, MPI_INT, j, tag1, MPI_COMM_WORLD, ierr)
      call MPI_SEND(Zr_MPI, N_Pass + N_PAR, MPI_REAL8, j, tag2, MPI_COMM_WORLD, ierr)
      call MPI_SEND(C_MPI, N_PAR + N_Pass, MPI_REAL8, j, tag3, MPI_COMM_WORLD, ierr)
   Enddo
ELSE

   ! Receive  Kv from root (child processes):
   call MPI_RECV(Kv, 1+nlev, MPI_REAL8, 0, tag4, MPI_COMM_WORLD, stat, ierr)

   ! Receive  w from root (child processes):
   call MPI_RECV(w, 1+nlev, MPI_REAL8, 0, tag5, MPI_COMM_WORLD, stat, ierr)

   !Receive Zi_Mpi, Zr_MPI and C_MPI from root
   call MPI_RECV(Zi_MPI, N_Par + N_Pass, MPI_INT, 0, tag1, MPI_COMM_WORLD, stat, ierr)
   call MPI_RECV(Zr_MPI, N_Par + N_Pass, MPI_REAL8, 0, tag2, MPI_COMM_WORLD, stat, ierr)
   call MPI_RECV(C_MPI, N_Par + N_Pass, MPI_REAL8, 0, tag3, MPI_COMM_WORLD, stat, ierr)
ENDIF

!Calculate dKv/dz
dKvdz(:) = 0d0
DO i = 1,nlev
   dKvdz(i)= (Kv(i)-Kv(i-1))/Hz(i)
ENDDO

dt = dtsec/dble(Nrand)
dt_inv = 1.d0/dt
rnd_var_inv = 1.d0/rnd_var

! Vertical random walk for passive particles
! Real parallel computing
DO j = (taskid*N_chunk + 1), ((taskid+1)*N_chunk)

   !Compute sinking rate (negative) based on Durante et al. J. Phycol 2019
   vs = phyto_sinking(C_MPI(j))

   !Add biological sinking to w
   bio_w(:)    = vs
   bio_w(0)    = 0d0
   bio_w(nlev) = 0d0
   bio_w = w + bio_w

   zp= Zr_MPI(j)
   zi= Zi_MPI(j)

   do iit = 1, Nrand

     !Judge whether the input is correct
     IF (zp < Z_w(0) .or. zp > Z_w(nlev) .or. zi < 1 .or. zi > nlev) THEN
        write(6,*) "zp = ", zp
        write(6,*) "zi = ", zi
        write(6,*) "Particle index = ", j
        stop "The particle out of the domain!"
     ENDIF
     
     ! local viscosity calculation
     ! correction suggested by Visser [1997]
     rat  = (zp-Z_w(zi-1))/Hz(zi)
     wint= rat*bio_w(zi) +(1.-rat)*bio_w(zi-1)  ! interpolate w
     
     ! The new position after dt/2
     zloc = zp + 0.5*(dKvdz(zi)+wint)*dt 
     
     ! First judge whether the particle jumps out of the domain
     if (zloc > Z_w(nlev)) then  !Reflective boundary at surface
        zloc = 2.d0*Z_w(nlev) - zloc 
     elseif (zloc < Z_w(0)) then  !Reflective boundary at bottom
        zloc = Z_w(0) + (Z_w(0) - zloc)
     endif      
     
     ! The following is to find the index for the grid where the particle currently resides
     do i = 1,nlev
        if (Z_w(i) .ge. zloc .and. Z_w(i-1) .le. zloc) then
            zi = i
            EXIT
        endif
     end do   !End of diffusivity correction
     
     rat  = (zloc-Z_w(i-1))/Hz(i)
     visc = rat*Kv(i)+(1.-rat)*Kv(i-1)  ! interpolate Kv
     wint= rat*bio_w(i) +(1.-rat)*bio_w(i-1)     ! interpolate w
     
     IF (visc .lt. Kbg) visc=Kbg   !Kv background
     
     CALL RANDOM_NUMBER(rnd)
     rnd  = 2.d0*rnd - 1.d0
     step1= SQRT(2d0 * rnd_var_inv * dt * visc)*rnd
     step2= dt*(wint+dKvdz(i))
     step = step1 + step2
     
     zp = zp + step
     
     ! First judge whether the particle jumps out of the domain
     IF (zp > Z_w(nlev)) THEN
           zp = 2.d0*Z_w(nlev) - zp   !Reflective boundary at surface
     ELSEIF (zp < Z_w(0)) THEN      !Reflective boundary at bottom
           zp = Z_w(0) + (Z_w(0) - zp)
     ENDIF      
     
     !Compute new zi
     do i=1,nlev
        if (Z_w(i) .ge. zp .and. Z_w(i-1) .le. zp) then
            zi = i
            EXIT
        endif
     end do
   enddo  !End of iit

   !Save zi and zp back to the particle
   Zi_MPI(j) = zi
   Zr_MPI(j) = zp
ENDDO  !End of particle j

IF (TASKID .GT. 0) THEN
  !Return the information of Zi_MPI and Zr_MPI back to the parent process
  call MPI_SEND(Zi_MPI((taskid*N_chunk + 1):(taskid+1)*N_chunk), N_chunk, MPI_INT, 0, tag6(TASKID), MPI_COMM_WORLD, ierr)
  call MPI_SEND(Zr_MPI((taskid*N_chunk + 1):(taskid+1)*N_chunk), N_chunk, MPI_REAL8,0,tag7(TASKID), MPI_COMM_WORLD, ierr)

ELSE
  do j = 1, numtasks-1
    call MPI_RECV(Zi_MPI((j*N_chunk + 1):(j+1)*N_chunk), N_chunk, MPI_INT, j, tag6(j), MPI_COMM_WORLD, stat, ierr)
    call MPI_RECV(Zr_MPI((j*N_chunk + 1):(j+1)*N_chunk), N_chunk, MPI_REAL8,j,tag7(j), MPI_COMM_WORLD, stat, ierr)
  enddo

  !Transfer back the data of Zi and Zr to p_PHY and p_Pass
  do m = 1, N_Pass + N_PAR
     if (m .le. N_Pass) then
         p_Pass(m)%rZ=Zr_MPI(m)
         p_Pass(m)%iz =Zi_MPI(m)
     else
         p_PHY(m-N_Pass)%rZ=Zr_MPI(m) 
         p_PHY(m-N_Pass)%iZ=Zi_MPI(m)
     endif
  enddo

ENDIF
return
END SUBROUTINE LAGRANGE

!Phytoplankton sinking rate follow Durante et al. JPR 2019
pure real function phyto_sinking(Vol) result(y)
implicit none
real, intent(in) :: Vol  !Cell volume (um^3) 
  y =  -(0.0019 * Vol**0.43)/86400.d0 !Unit: m per second (downwards)
END function phyto_sinking
