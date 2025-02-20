PROGRAM IBM
USE GRID
USE MPI_Setting
IMPLICIT NONE
real(4) :: start,finish

! ***** Initialize MPI *****
call MPI_INIT(ierr)

! Returns the size of the group associated with a communicator.
call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

! Determines the rank of the calling process in the communicator
call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)

!Count time
call cpu_time(start) 

!Start the pseudo random number generator
call random_seed()

!Define grid
call setup_grid

!Initialization
call initialize

!Timestep
call Timestep

!End mpi
call MPI_finalize(ierr)

call cpu_time(finish)
if (taskid.eq.0) print '("Whole simulation time = ",f8.3," hours.")', (finish-start)/3600.0 

end program
