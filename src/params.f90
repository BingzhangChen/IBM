module params
implicit NONE

! Parameters
real :: mu0     =  0d0  !PCref: maximal growth rate normalized to 15 C
real :: aI0     =  0d0  !Slope of P-I curve
real :: KN      =  0d0  !Half saturation constant of phytoplankton growth on N
real :: gmax    =  0d0  !Maximal zooplankton grazing rate
real :: Kp      =  0d0  !Half saturation constant of zooplankton grazing
real :: mz      =  0d0  !Zooplankton mortality grazing coefficient
real :: RdN     =  0d0  !Rate of detritus converting to DIN
real :: wDET    =  0d0 ! Sinking rate of detritus

! Activation energy for phytoplankton growth
real, parameter :: Ep  = 0.32d0

! Activation energy for zooplankton grazing
real, parameter :: Ez  = 0.65d0

! Maixmal Chl:N ratio (g:mol) for phytoplankton
real, parameter :: thetaNmax = 3d0

!the value of rhochl before last sunset
real            :: rhochl_L  = 0d0
end module params
