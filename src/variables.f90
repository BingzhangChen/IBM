module state_variables
use grid, only : nlev, Hz
implicit NONE

!State variables

!Indexes for the state variables
integer, parameter :: iNO3 = 1 
integer, parameter :: iPC  = iNO3 + 1 
integer, parameter :: iPN  = iPC  + 1 
integer, parameter :: iCHL = iPN  + 1 
integer, parameter :: iZOO = iCHL + 1 
integer, parameter :: iDET = iZOO + 1 
integer, parameter :: nvar = iDET !Total number of state variables
real               :: t(nvar, nlev) = 0.d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  1  ! DET

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [iDET]

!Sinking rate
real               :: ww(0:nlev, NVsinkterms) = 0.d0

!Labels for each state variable (for saving model output)
integer, parameter :: oNPP               = nvar + 1
integer, parameter :: Nout               = oNPP
real               :: Varout(Nout, nlev) = 0.d0
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

    ! cellular nitrogen content (pmol; assuming a 1 micron cell)
    real    :: N   = 0.02/106.*16.

    ! Cellular Chl content (pg Chl)
    real    :: Chl = 0.02 * 12/50

    ! Number of cells per superindividual
    real    :: num = 5d9 

END TYPE Particle

!Fix the number of individuals in the system
INTEGER, PARAMETER           :: N_PAR = 2000
Type (Particle)              :: p_PHY(N_PAR)

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

END MODULE