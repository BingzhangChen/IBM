Module Time_setting
implicit NONE

real,    parameter :: d_per_s = 864d2  ! how many seconds in one day
real,    parameter :: s_per_h = 3600d0 ! how many seconds in one hour
integer, parameter :: h_per_d = 24     ! how many hours in one day
integer, parameter :: y_per_d = 365    ! how many days of one year
integer            :: NDay_Run= 365    !Total number of days for simulation

!Time of the day
integer  :: sec_of_day   = 0
integer  :: current_day  = 0
integer  :: current_DOY  = 0
integer  :: current_hour = 0
integer  :: current_year = 0
real       :: current_sec  = 0.d0

!The current step
integer  :: it           = 0

!Total number of biological time steps
integer  :: Nstep        = 0

!Number of time steps of random walk per biological time step
integer  :: Nrand        = 0

!Timestep as a fraction of the day
real     :: dtdays       = 0.d0

!Timestep as seconds
real     :: dtsec        = 0.d0

!Number of iterations to save. Save everyday (86400/dtsec); save every hour (3600/dtsec)
integer  :: nsave        = 0 ! 

!Settings of reading previous outputs
integer :: read_previous_output = 0
integer :: NR_Euler = 0
integer :: Last_timestep = 0
integer :: Last_Day = 0
integer :: Last_hr = 0

CONTAINS

subroutine update_time
!This subroutine converts current time step to the current second, second of the day, current DOY etc.
implicit none

! Calculate current timing (zero is starting time):
current_sec = dble(it-1)*dtsec

! Calculate time of the day
sec_of_day = mod(nint(current_sec), nint(d_per_s))

! Calculate DOY
! Calculate model time in days:
current_day = int(current_sec/d_per_s) + 1

! Calculate current year
current_year= int(current_day/y_per_d) + 1

! Calculate DATE OF the YEAR (DOY)
Current_DOY = mod(current_day, y_per_d)

! Calculate current hour
current_hour= int(current_sec/s_per_h)
current_hour= mod(current_hour, h_per_d)

end subroutine
end module
