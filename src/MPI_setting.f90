MODULE MPI_Setting
! This module calls the mpi module and defines general variables for OpenMPI
USE MPI
IMPLICIT NONE

! MPI variables:
integer  :: numtasks, ierr, MPIRUN, taskid
integer, allocatable :: tag6(:)
integer, allocatable :: tag7(:)

integer :: N_chunk = 0  !The number of particle units to be transferred between parent and child processes, initialized at initialize.f90
End MODULE

