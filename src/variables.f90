MODULE STATE_VARIABLES
use grid, only : nlev, Hz
IMPLICIT NONE

!State variables
integer, private   :: i
integer, parameter :: NZOO = 20   !Number of zooplankton size classes

!Indexes for the state variables
integer, parameter :: iNO3 = 1 
integer, parameter :: iPC  = iNO3 + 1 
integer, parameter :: iPN  = iPC  + 1 
integer, parameter :: iCHL = iPN  + 1 
integer, parameter :: iZOO(NZOO) = (/ (iCHL + i , i = 1, NZOO) /)
integer, parameter :: iDET = iZOO(NZOO) + 1 
integer, parameter :: nvar = iDET !Total number of state variables
real               :: t(nvar, nlev) = 0.d0
real               :: Ntot = 0d0  !Total nitrogen in the domain
real, parameter :: MinSzoo = log(0.8d0)!Minimal zooplankton log ESD (micron)
real, parameter :: MaxSzoo = log(25d2) !Maximal zooplankton log ESD (micron)
real, parameter :: dZOOESD =  (MaxSzoo - MinSzoo)/dble(NZOO-1)! ESD difference between adjacent zooplankton size class (log)

!Log ESD of each zoo. size class
real                    :: ESDZOO(NZOO) = 0d0

!Volume of each zooplankton size class
real                    :: VolZOO(NZOO) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  1  ! DET

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [iDET]

!Sinking rate
real               :: ww(0:nlev, NVsinkterms) = 0.d0

!Labels for each state variable (for saving model output)
integer, parameter   :: oNPP = nvar + 1   !Daily net primary production integrated over a whole day

! Total available prey biomass for each zoo. size class
integer, parameter :: oFZ(NZOO) = (/ (oNPP + i , i = 1, NZOO) /)
integer, parameter :: oTEMP     = oFZ(NZOO) + 1
integer, parameter :: oPAR      = oTEMP     + 1
integer, parameter :: Nout      = oPAR
real                              :: Varout(Nout, nlev) = 0.d0
character(LEN=7)   :: Labelout(Nout)     = 'Unknown'

!Particles
!Declaration for phyto particles
TYPE Particle

    ! Particle ID
    integer :: ID = 1

    ! Grid indices for particles (range from nlev to 1)
    integer :: iz = 1
    
    ! Z coordinates for particles
    real    :: rz = 0.

    ! associated PAR 
    real    :: PAR = 0.01

    ! associated Temperature
    real    :: temp= 20.

    ! associated nutrient
    real    :: NO3 = 0.1

    ! cellular carbon content (pmol; assuming a 1 micron cell)
    real    :: C   = 0.02

    ! cellular nitrogen content (pmol per cell; assuming a 1 micron cell)
    real    :: N   = 0.02/106.*16.

    ! Cellular Chl content (pg Chl)
    real    :: Chl = 0.02 * 12/50

    ! Number of cells per superindividual
    real    :: num = 5d9 

    !cellular carbon content threshold for division (pmol/cell), can be used as a proxy for size and can be converted to ESD; Phytoplankton half-saturation constant, minimal N:C and maximal N:C ratios are allometric functions of this parameter
    !This trait will vary with mutation 
    real :: Cdiv= 0.04d0
    
    !Dead or alive
    logical :: alive = .true.

    !Optimal temperature
    real :: Topt = 20.d0

    !Ln alphaChl (slope of the P-I curve; unit: (W m-2)-1 (gChl molC)-1 d-1 instead of micro mol quanta m-2 s-1)
    real :: LnalphaChl = -2.3  !log(0.1)

END TYPE Particle

!Passive particles for diagosing random walk
TYPE Passive_Particle

    ! Particle ID
    integer :: ID = 1
    ! Grid indices for particles (range from nlev to 1)
    integer :: iz = 1
    
    ! Z coordinates for particles
    real    :: rz = 0.

End Type Passive_Particle

integer :: IDmax = 0 !Maximal ID number

!Fix the number of individuals in the system
INTEGER           :: N_PAR = 10000  !Number of phyto. particles
INTEGER           :: N_Pass= 1000    !Number of passive particles
Type (Particle), ALLOCATABLE       :: p_PHY(:)
Type (Passive_Particle), ALLOCATABLE :: p_pass(:)

!Model choices
integer, parameter :: GMK98_simple                 = 1 
integer, parameter :: GMK98_Topt                    = 2 
integer, parameter :: GMK98_Size                     = 3 
integer, parameter :: GMK98_Light                  = 4 
integer, parameter :: GMK98_ToptLight          = 5 
integer, parameter :: GMK98_ToptSize            = 6 
integer, parameter :: GMK98_SizeLight          = 7 
integer, parameter :: GMK98_ToptSizeLight  = 8 

!Current model selection
integer :: Model_ID = 8

!Number of traits
integer, parameter :: NTrait = 3
integer, parameter :: iTopt = 1     !Trait index for Topt
integer, parameter :: iSize = 2     !Trait index for Size (ESD)
integer, parameter :: ialphaChl = 3    !Trait index for optimal light

real :: nu(NTrait) = [1d-12, 1d-12, 1d-12] !Probability per generation per cell
real :: sigma(NTrait) = [0.1, 0.1, 0.1]    !Standard deviation of mutation of the three traits

integer :: N_birth(nlev) = 0  !Number of birth events during one hour
integer :: N_death(nlev) = 0  !Number of death events during one hour
integer :: N_mutate(nlev)= 0  !Number of mutation events during one hour at each grid

CONTAINS

SUBROUTINE UPDATE_PHYTO
IMPLICIT NONE
integer :: k, j_

!Update concentrations of phytoplankton carbon, nitrogen, and Chl
t(iPC, :) = 0.d0   !Unit: mmolC/m3
t(iPN, :) = 0.d0   !Unit: mmolN/m3
t(iChl,:) = 0.d0   !Unit: mgChl/m3

DO k = 1, N_PAR
   do j_ = 1, nlev
      if (p_PHY(k)%iz == j_) then
         t(iPC, j_) = t(iPC, j_) + p_PHY(k)%num * p_PHY(k)%C   * 1d-9  !Unit: mmolC 
         t(iPN, j_) = t(iPN, j_) + p_PHY(k)%num * p_PHY(k)%N   * 1d-9  !Unit: mmolN 
         t(iChl, j_)= t(iChl, j_)+ p_PHY(k)%num * p_PHY(k)%Chl * 1d-9  !Unit: mgChl 
         exit
      endif
   enddo
ENDDO

do j_ = 1, nlev
   t(iPC, j_) = t(iPC, j_) /Hz(j_) !Change unit to mmolC/m3
   t(iPN, j_) = t(iPN, j_) /Hz(j_) !Change unit to mmolN/m3
   t(iCHL,j_) = t(iCHL, j_)/Hz(j_) !Change unit to mgChl/m3
Enddo
END SUBROUTINE UPDATE_PHYTO

END MODULE STATE_VARIABLES
