Module GRID
integer, parameter :: nlev = 100  !Number of vertical layers

real,    parameter :: hmax = 250d0  !Maximal depth (positive)

!Grid position (m, negative) at r points 
real               :: Z_r(nlev)  = 0.d0 

!Grid position (m, negative) at w points 
real               :: Z_w(0:nlev)= 0.d0 

!Grid height (m, positive)
real               :: Hz(nlev)   = 0.d0

contains

subroutine setup_grid
implicit none

real,    parameter  :: thetaS = 2.d0!Surface stretching parameter

! Local scratch variable:
integer                :: i
real                   :: sc_r, C_sig

Z_w(0) = -hmax

!Following Song and Haidvogel (1994). sinh is the hyperbolic sin function
do i   = 1,nlev
    sc_r  = (float(i-nlev) - 0.5)/float(nlev)
   C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
   Z_r(i) = C_sig*hmax
    sc_r  = (float(i-nlev))/float(nlev)
   C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
   Z_w(i) = C_sig*hmax
    Hz(i) = Z_w(i) - Z_w(i-1)
enddo
end subroutine setup_grid
end module GRID