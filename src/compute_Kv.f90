!Calculates analytic vertical eddy diffusivity profiles
subroutine analytic_Kv(nlev, Ksurf, Km, Kbg, H, Kv)
use grid, only : Z_w
implicit none
integer, intent(in) :: nlev
real, intent(in)  :: Ksurf  ! Surface diffusivity
real, intent(in)  :: Km     ! maximal diffusivity
real, intent(in)  :: Kbg    ! Background diffusivity
real, intent(in)  :: H      ! MLD
real, intent(out) :: Kv(0:nlev)  ! Output Kv file

real :: z !Water depth

real :: cff = 0. !Scratch variable

integer :: k  

Kv(nlev) = Ksurf
do k = 1, nlev
   z = Z_w(k) !Current depth
   if (z .gt. -H/2.d0) then
      cff  = -abs(z)**2.3
      if (cff .lt. -200.) then
         Kv(k)=Km
      else 
         Kv(k)=Ksurf + (Km - Ksurf)*( 1.d0 - 1.5**cff )
      endif

   elseif (z .le. -H/2d0 .and. z .ge. -H) then
      cff  = -(H+z)**2.3

      if (cff .lt. -200.) then
         Kv(k)=Km
      else 
         Kv(k)=Kbg + (Km - Kbg)*( 1.d0 - 1.5**cff )
      endif
   else
      Kv(k)=Kbg
   endif
enddo
END subroutine analytic_Kv
