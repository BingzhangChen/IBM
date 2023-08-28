PROGRAM IBM
USE GRID
USE MPI_Setting
IMPLICIT NONE
real(4) :: start,finish, t1

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
call cpu_time(t1) 

if (taskid .eq. 0) then
  print '("Initialization costs ",f8.3," hours.")', (t1-start)/3600.0 
endif

stop

!Timestep
call Timestep

!End mpi
call MPI_finalize(ierr)

call cpu_time(finish)
if (taskid.eq.0) print '("Whole simulation time = ",f8.3," hours.")', (finish-start)/3600.0 

end program
