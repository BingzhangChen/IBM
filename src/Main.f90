program IBM
use GRID
implicit none
real(4) :: start,finish

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

call cpu_time(finish)
print '("Time = ",f8.3," hours.")', (finish-start)/3600.0 

end program