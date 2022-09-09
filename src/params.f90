MODULE PARAMS
implicit NONE

! Parameters
real :: mu0     =  0d0  !Maximal growth rate normalized to 15 C (d-1)
real :: aI0     =  0d0  !Chl-specific Slope of P-I curve ((W m-2)-1 (gChl molC)-1 d-1)
real :: KN      =  0d0  !Half saturation constant of phytoplankton growth on N
real :: gmax    =  0d0  !Maximal zooplankton grazing rate
real :: Kp      =  0d0  !Half saturation constant of zooplankton grazing
real :: mz      =  0d0  !Zooplankton mortality grazing coefficient
real :: GGE     =  0d0  !Zooplankton Gross Growth Efficiency [dimensionless]
real :: unass   =  0d0  !Fraction of unassimilated food ingested by zooplankton [nd]
real :: RdN     =  0d0  !Rate of detritus converting to DIN
real :: wDET    =  0d0 ! Sinking rate of detritus

! Standard deviation of log zooplankton feeding preference
real :: SDZoo   =  0d0 

! Activation energy for phytoplankton growth
real, parameter :: Ep  = 0.32d0

! Activation energy for zooplankton grazing
real, parameter :: Ez  = 0.65d0

! Maixmal Chl:N ratio (g:mol) for phytoplankton
real, parameter :: thetaNmax = 3d0

!the value of rhochl before last sunset
real            :: rhochl_L  = 0d0
end module params
