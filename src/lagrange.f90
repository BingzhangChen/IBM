! Phytoplankton and passive particles random walk
!TODO: Needs to implement Ross & Sharples (2004)
SUBROUTINE LAGRANGE
USE Time_setting, only: dtsec, Nrand
USE grid,                only: nlev, Z_w, Hz
USE forcing,           only: Kv, dKvdz, w, Kbg
USE STATE_VARIABLES, only: p_PHY, p_pass, N_PAR, N_pass
USE Trait_functions,      only : PHY_C2Vol
IMPLICIT NONE

REAL     :: bio_w(0:nlev) = 0d0  !Scratch sinking rate in 1D

! !LOCAL VARIABLES:
INTEGER         :: i,j, zi, iit
integer           :: m = 0
REAL               :: rnd, rnd_var_inv, dt, step1, step2, zloc, zp
REAL               :: dzn(nlev),step
REAL               :: visc,rat,dt_inv, vs
REAL               :: wint  !interpolated w
REAL,parameter :: rnd_var=0.333333333
REAL, external  :: phyto_sinking
!-----------------------------------------------------------------------

dt = dtsec/dble(Nrand)
dt_inv = 1.d0/dt
rnd_var_inv = 1.d0/rnd_var

! Vertical random walk for passive particles
DO j = 1, (N_Pass + N_par)

   if (j .le. N_Pass) then  !Passive particles
      zp = p_pass(j)%rz
      zi = p_pass(j)%iz 
      bio_w = w  !No additional vertical movement

   else  !Phytoplankton particles

      m = j - N_pass  !The index for phytoplankton particles

      if (.not. p_PHY(m)%alive) cycle  !Only select live phyto. particles

      zp = p_PHY(m)%rz
      zi = p_PHY(m)%iz 

      !Compute sinking rate (negative) based on Durante et al. J. Phycol 2019
      vs = phyto_sinking(PHY_C2Vol(p_PHY(m)%C)) 
      
      !Add biological sinking to w
      bio_w(:)    = vs
      bio_w(0)    = 0d0
      bio_w(nlev) = 0d0
      bio_w = w + bio_w
   endif

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
   
   if (j .le. N_Pass) then  !Passive particles
      p_pass(j)%rz = zp
      p_pass(j)%iz = zi
   else   !Phyto. particles
      p_PHY(m)%rz = zp
      p_PHY(m)%iz = zi
   endif
ENDDO  !End of particle j

return
END SUBROUTINE LAGRANGE

!Phytoplankton sinking rate follow Durante et al. JPR 2019
pure real function phyto_sinking(Vol) result(y)
implicit none
real, intent(in) :: Vol  !Cell volume (um^3) 
  y =  -(0.0019 * Vol**0.43)/86400.d0 !Unit: m per second (downwards)
END function phyto_sinking
