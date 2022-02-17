Module IO
implicit none

private

character(len=9), parameter :: Euler_file    = "Euler.out"
integer,          parameter ::    Euler_unit = 10
integer,          parameter :: Particle_unit = 11

public :: create_Eulerian_file, create_Particle_file
public :: save_Eulerian, save_particles

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_Eulerian_file
use GRID, only : Z_r
implicit none

! Create data out file
open (unit=Euler_unit, file = Euler_file, status = 'replace')
write(Euler_unit, 101) 'Type','Timestep','Day', 'Hour', Z_r
close(Euler_unit)

101 format(4(A10,1x), 200(2x,F12.5))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine save_Eulerian
use state_variables, only : Labelout, Varout, Nout
use Time_setting,    only : it, current_day, current_hour
implicit none
integer :: i

! save data into the Eulerian file
open (unit=Euler_unit, file = Euler_file, status = 'old', action='write', position='append')
do i = 1, Nout
   write(Euler_unit, 103) trim(Labelout(i)), it, current_day, current_hour, Varout(i,:)
enddo
close(Euler_unit)

103 format(A3, 1x, I7, 1x, 2(I0, 1x), 200(1pe12.3, 1x))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_Particle_file(fname)
implicit none
character(len=*), intent(in) :: fname  !Need to create multiple files for storing particles to avoid huge files

open (unit=Particle_unit, file = fname, status = 'replace')
write(Particle_unit, 102) 'Timestep','Day','Hour','ID','Grid', &
                          'Z', 'C', 'N', 'Chl', 'NO3', 'PAR',  &
                          'Temp', 'Num', 'Topt' 
close(Particle_unit)

return

102 format(14(A8, 2x))
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_particles(fname)
use state_variables, only : p_PHY
use Time_setting,    only : it, current_day, current_hour
implicit none
character(len=*), intent(in) :: fname  !The file for saving particles
integer :: i

! save data into the particle file
open (unit=Particle_unit, file=fname, status='old', action='write', position='append')
do i = 1, size(p_PHY)
   write(Particle_unit, 104) it,                             &
                  current_day, current_hour, p_PHY(i)%ID, &
                  p_PHY(i)%iz, p_PHY(i)%rz,  p_PHY(i)%C,  &
                  p_PHY(i)%N,  p_PHY(i)%Chl, p_PHY(i)%NO3,&
                  p_PHY(i)%PAR,p_PHY(i)%temp,p_PHY(i)%num, p_PHY(i)%Topt
enddo
close(Particle_unit)

104 format(I0, 1x, 2(I8, 1x), I0, 1x, I3, 1x, F12.3, 1x, 8(1pe12.3, 1x))
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE IO