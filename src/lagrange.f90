! Lagrange particle tracking on one particle
! !INTERFACE:
!TODO: Needs to implement Ross & Sharples (2004)
SUBROUTINE LAGRANGE(nlev,zlev, nuh,w,zi,zp)
USE Time_setting, only: dtdays, d_per_s, Nrand
! dtdays: timestep in days
! d_per_s: how many seconds in one day
!
IMPLICIT NONE
!
! !INPUT PARAMETERS:
INTEGER, intent(in)                :: nlev          !Number of vertical grids
REAL,    intent(in)                :: zlev(0:nlev)  !Depth of vertical grids at w points
REAL,    intent(in)                :: nuh(0:nlev)   !Vertical diffusivity at w points
REAL,    intent(in)                :: w(0:nlev)     !Vertical advection at w points

! INPUT/OUTPUT PARAMETERS:
! the ith grid where the particle currently resides
INTEGER, intent(inout)             :: zi

! Z position of particles
REAL,    intent(inout)             :: zp
!
! !LOCAL VARIABLES:
INTEGER            :: i
REAL               :: rnd, rnd_var_inv, dt, step1, step2
REAL,parameter     :: visc_back=0.e-6,rnd_var=0.333333333
REAL               :: dz(nlev),dzn(nlev),step
REAL               :: visc,rat,dt_inv,zloc
REAL               :: wint  !interpolated w
!EOP
!-----------------------------------------------------------------------
!BOC

dt = dtdays * d_per_s/dble(Nrand)

!Judge whether the input is correct
IF (zp < zlev(0) .or. zp > zlev(nlev) .or. zi < 1 .or. zi > nlev) THEN
   write(6,*) "zp = ", zp
   write(6,*) "zi = ", zi
   stop "The particle out of the domain!"
ENDIF

dt_inv      = 1.d0/dt
rnd_var_inv = 1.d0/rnd_var

DO i = 1,nlev
   
   ! Depth of each grid
   dz(i) = zlev(i)-zlev(i-1)

   ! gradient of Kv
   dzn(i)= (nuh(i)-nuh(i-1))/dz(i)

END DO

!  local viscosity calculation
! correction suggested by Visser [1997]
rat  = (zp-zlev(zi-1))/dz(zi)
wint = rat*w(zi) +(1.-rat)*w(zi-1)  ! interpolate w

! The new position after dt/2
zloc = zp + 0.5*(dzn(zi)+wint)*dt 

! First judge whether the particle jumps out of the domain
if (zloc > zlev(nlev)) then  !Reflective boundary at surface
   zloc = 2.d0*zlev(nlev) - zloc 
elseif (zloc < zlev(0)) then  !Reflective boundary at bottom
   zloc = zlev(0) + (zlev(0) - zloc)
endif      

! The following is to find the index for the grid where the particle currently resides
do i = 1,nlev
   if (zlev(i) .ge. zloc .and. zlev(i-1) .le. zloc) then
       zi = i
       EXIT
   endif
end do   !End of diffusivity correction

rat  = (zloc-zlev(i-1))/dz(i)
visc = rat*nuh(i)+(1.-rat)*nuh(i-1)  ! interpolate Kv
wint = rat*w(i) +(1.-rat)*w(i-1)     ! interpolate w

IF (visc.lt.visc_back) visc=visc_back   !Kv background

CALL RANDOM_NUMBER(rnd)
rnd  = 2.d0*rnd - 1.d0
step1= SQRT(2d0 * rnd_var_inv * dt * visc)*rnd
step2= dt*(wint+dzn(i))
step = step1 + step2

zp = zp + step

! First judge whether the particle jumps out of the domain
IF (zp > zlev(nlev)) THEN
      zp = 2.d0*zlev(nlev) - zp   !Reflective boundary at surface
ELSEIF (zp < zlev(0)) THEN      !Reflective boundary at bottom
      zp = zlev(0) + (zlev(0) - zp)
ENDIF      

!Compute new zi
do i=1,nlev
   if (zlev(i) .ge. zp .and. zlev(i-1) .le. zp) then
       zi = i
       EXIT
   endif
end do

return
END SUBROUTINE LAGRANGE
