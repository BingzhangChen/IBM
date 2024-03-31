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

! Acclimation parameter
real, parameter :: nx = 1d0

! Activation energy for phytoplankton growth
real, parameter :: Ep  = 0.32d0

! Activation energy for zooplankton grazing
real, parameter :: Ez  = 0.65d0

! Maixmal Chl:N ratio (g:mol) for phytoplankton
real, parameter :: thetaNmax = 3d0

!the value of rhochl before last sunset
real            :: rhochl_L  = 0d0

!QNmin and QNmax are allometric functions of Vol (Ward et al. 2012) [mol N: mol C]:
real, parameter   :: QNmin_a = 0.07d0     !Normalization constant for QNmin [molN:molC]
real, parameter   :: QNmin_b = -0.17d0    !Allometric exponent for QNmin
real, parameter   :: QNmax_a = 0.25d0     ! Normalization constant for QNmax [molN:molC]
real, parameter   :: QNmax_b = -0.07d0    ! Allometric exponent for QNmax

!Kn is an allometric function of Vol (Cdiv) (Edwards et al. 2012) [uM]:
real, parameter   :: KN_a   = 0.14d0      !Normalization constant for KN
real, parameter   :: KN_b   = 0.33d0      !Allometric exponent for KN

!Mutation parameters
integer, parameter :: NTrait = 3
real,    parameter :: nu(NTrait)    = [1d-12, 1d-12, 1d-12] !Probability per generation per cell
real,    parameter :: sigma(NTrait) = [0.1, 0.1, 0.1]       !Standard deviation of mutation of the three traits

real, parameter   :: pi= 3.1415926535

!The variance of random noise added to VKv
real, parameter :: VAR_noise_Kv = 0d0
end module params
