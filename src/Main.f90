PROGRAM IBM
USE GRID
IMPLICIT NONE
real(4) :: start,finish, t1

!Count time
call cpu_time(start) 

!Start the pseudo random number generator
call random_seed()

!Define grid
call setup_grid

!Initialization
call initialize

call cpu_time(t1) 
print '("Initialization costs ",f8.3," hours.")', (t1-start)/3600.0 

!Timestep
call Timestep

call cpu_time(finish)
print '("Whole simulation time = ",f8.3," hours.")', (finish-start)/3600.0 

end program