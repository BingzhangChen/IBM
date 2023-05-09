Module IO
implicit none

private

character(len=9), parameter :: Euler_file    = "Euler.out"
character(len=6), parameter :: Kv_file       = "Kv.out"
character(len=9), parameter :: Death_file    = "Death.out"
integer,          parameter :: Euler_unit    = 10
integer,          parameter ::    Kv_unit    = 12
integer,          parameter :: Phyto_unit    = 13
integer,          parameter :: Passive_unit  = 14
integer,          parameter ::  Death_unit   = 15
integer,          parameter :: phyto         = 1
integer,          parameter :: passive       = 2
integer,          parameter :: death         = 3

public :: create_Eulerian_file, create_Particle_file, create_Kv_file
public :: save_Eulerian, save_particles, save_Kv, phyto, passive, death
public :: Death_file

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_Eulerian_file
use GRID, only : Z_r
implicit none

! Create data out file
open (unit=Euler_unit, file = Euler_file, status = 'replace')
write(Euler_unit, 101) 'Type','Timestep','Day', 'Hour', Z_r
close(Euler_unit)

101 format(4(A10,1x), 200(2x,F12.5))
end subroutine create_Eulerian_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_Kv_file
use GRID, only : Z_w
implicit none

! Create data out file
open (unit=Kv_unit, file = Kv_file, status = 'replace')
write(Kv_unit, 111) 'Timestep','Day', 'Hour', Z_w
close(Kv_unit)

111 format(3(A10,1x), 200(2x,F12.5))
end subroutine create_Kv_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine save_Eulerian
use state_variables, only : Labelout, Varout, Nout, oN_ind
use Time_setting,    only : it, current_day, current_hour
implicit none
integer :: i,j

! save data into the Eulerian file 
open (unit=Euler_unit, file = Euler_file, status = 'old', action='write', position='append')
do i = 1, Nout
   write(Euler_unit, 103) trim(Labelout(i)), it, current_day, current_hour, Varout(i,:)
enddo
close(Euler_unit)

103 format(A5, 1x, I7, 1x, 2(I0, 1x), 200(1pe12.3, 1x))
end subroutine save_Eulerian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine save_Kv
use forcing,         only : Kv
use Time_setting,    only : it, current_day, current_hour
implicit none
integer :: i

! save data into the Kv file
open (unit=Kv_unit, file = Kv_file, status = 'old', action='write', position='append')
write(Kv_unit, 113) it, current_day, current_hour, Kv(:)
close(Kv_unit)

113 format(I7, 1x, 2(I0, 1x), 200(1pe12.3, 1x))
end subroutine save_Kv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_Particle_file(fname, type)
implicit none
character(len=*), intent(in) :: fname  !Need to create multiple files for storing particles to avoid huge files
integer, intent(in) :: type  !type = 1, phyto; !type = 2, passive particles;
                             !type = 3, Birth and death file
integer :: fileunit = 0

select case(type)
case (phyto)
   fileunit = Phyto_unit
   open (unit=fileunit, file = fname, status = 'replace')
   write(fileunit, 102) 'Timestep','Day','Hour','ID','Grid', &
                               'Z', 'C', 'N', 'Chl', 'NO3', 'PAR',  &
                               'Temp', 'Num', 'Topt', 'CDiv', 'LnAlpha' 
case (passive)
   fileunit = passive_unit
   open (unit=fileunit, file = fname, status = 'replace')
   write(fileunit, 102) 'Timestep','Day','Hour','ID','Grid', 'Z'
case (death)
   fileunit = death_unit
   open (unit=fileunit, file = fname, status = 'replace')
   write(fileunit, 102) 'Timestep','Day','Hour','Grid','N_birth','N_death','N_mutate'

case default
   stop "Particle file type incorrect!"
end select
close(fileunit)
return

102 format(16(A8, 2x))
END subroutine create_Particle_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_particles(fname, type)
use state_variables, only : p_PHY, p_pass, N_birth, N_death, N_mutate
use Time_setting,    only : it, current_day, current_hour
use Grid,            only : nlev, Z_r

implicit none
character(len=*), intent(in) :: fname  !The file for saving particles
integer, intent(in) :: type
integer :: i, fileunit

selectcase (type)
case (phyto)
   ! save data into the phyto. file
   fileunit = Phyto_unit
   open (unit=fileunit, file=fname, status='old', action='write', position='append')
   do i = 1, size(p_PHY)
      write(fileunit, 104) it,                             &
                     current_day, current_hour, p_PHY(i)%ID, &
                     p_PHY(i)%iz, p_PHY(i)%rz,  p_PHY(i)%C,  &
                     p_PHY(i)%N,  p_PHY(i)%Chl, p_PHY(i)%NO3,&
                     p_PHY(i)%PAR,p_PHY(i)%temp,p_PHY(i)%num, &
                     p_PHY(i)%Topt, p_PHY(i)%CDiv, p_PHY(i)%LnalphaChl
   enddo
case (passive)
   fileunit = Phyto_unit
   open (unit=fileunit, file=fname, status='old', action='write', position='append')
   do i = 1, size(p_pass)
      write(fileunit, 104) it,                             &
                     current_day, current_hour, p_pass(i)%ID, &
                     p_pass(i)%iz, p_pass(i)%rz
   enddo

case (death)
   fileunit = death_unit
   open (unit=fileunit, file=fname, status='old', action='write', position='append')
   do i = 1, nlev
     write(fileunit, 105) it, current_day, current_hour, Z_r(i), &
                          N_birth(i), N_death(i), N_mutate(i)

     !Reset N_birth, N_death, and N_mutate
     N_birth(i) = 0
     N_death(i) = 0
     N_mutate(i)= 0
   enddo

case default
   stop "Particle file type incorrect!"

endselect
close(fileunit)

104 format(I0, 1x, 2(I8, 1x), I0, 1x, I3, 1x, F12.3, 1x, 10(1pe12.3, 1x))
105 format(I0, 1x, 2(I8, 1x), F12.3, 1x, 3(I0, 1x))
END SUBROUTINE save_particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE IO
