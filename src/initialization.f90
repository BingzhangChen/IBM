SUBROUTINE INITIALIZE
use params
use state_variables
use Time_setting, only: d_per_s, dtsec, dtdays, Nstep, nsave, NDay_Run
use grid,         only: Z_w
use IO
use forcing
implicit NONE

integer            :: rc       = 0
integer            :: k        = 0
integer            :: j_       = 0
real               :: cff      = 0.d0
integer, parameter :: namlst   = 20   !Unit time for namelist files
character(len=10)  :: par_file = 'Filename'
!==========================================================

!Namelist definition of time settings
namelist /timelist/  NDay_Run, dtsec, nsave

!Namelist definition of model choice and parameters
namelist /paramlist/ Model_ID, mu0, aI0, KN, gmax, Kp, mz,RDN, wDET 

write(6,*) 'Initialize model simulation time...'
! Check whether the namelist file exists.
inquire (file='time.nml', iostat=rc)

if (rc /= 0) then
    write (6, '(a)') 'Error: namelist file time.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='time.nml',status='old',action='read')
read(namlst,nml=timelist)
close(namlst)


!Update total number of time ! Number of time steps
Nstep  = NDay_Run*INT(d_per_s)/INT(dtsec) 

!Calculate dtdays
dtdays = dtsec/d_per_s
write(6,'(A13,1x,1pe12.2,A12)') 'Timestepping: ', dtdays, 'of one day.'

!==========================================================
write(6,*) 'Initialize model parameters...'
!Read parameter namelist

! Check whether the namelist file exists.
inquire (file='param.nml', iostat=rc)

if (rc /= 0) then
    write (6, '(a)') 'Error: namelist file Model.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='param.nml',status='old',action='read')
read(namlst,nml=paramlist)
close(namlst)
write(6,'(a)') 'Successfully initialize model parameters.'

write(6,'(A15,1x,I1)') 'Select Model ID', Model_ID

!Initialize state variables
write(6,'(a)') 'Initialize state variables...'

! Initialize initial NO3:
t(iNO3,:) = 0.5d0
t(iZOO,:) = 0.1d0

!Following Verity et al. AME (1996)
t(iDET,:) = t(iZOO,:) 

!For lagrangian model, initialize individuals
DO k = 1, N_PAR
   p_PHY(k)%ID = k        
   p_PHY(k)%alive = .true.

   !Initialize optimal temperature (Topt) to 20 C
   p_PHY(k)%Topt = 20.d0

   !Initialize ESD to 2 micron
   p_PHY(k)%LnESD = log(2.d0)

   !Initialize Iopt to 1000 umol photons m-2 s-1
   p_PHY(k)%LnIopt = log(1d3)

   !Compute rx (randomly distributed between Z_w(0) and Z_w(nlev)
   CALL RANDOM_NUMBER(cff)
   p_PHY(k)%rz = cff *(Z_w(nlev)-Z_w(0)) + Z_w(0)

   ! Find iz
   DO j_ = 1, nlev
      IF (Z_w(j_-1) .le. p_PHY(k)%rz .AND. Z_w(j_) .gt. p_PHY(k)%rz) THEN
         p_PHY(k)%iz = j_
         EXIT            
      ENDIF
   ENDDO
ENDDO

IDmax=N_par

call UPDATE_PHYTO  !Initialize Eulerian concentrations of phytoplankton carbon, nitrogen and Chl

!Initialiaze labels for model output
Labelout(iNO3) = 'NO3'
Labelout(iZOO) = 'ZOO'
Labelout(iDET) = 'DET'
Labelout(iPC)  = 'PC '
Labelout(iPN)  = 'PN '
Labelout(iCHL) = 'CHL'
Labelout(oNPP) = 'NPP'

write(6, '(a)') 'Write out the labels for output:'
do k = 1, Nout
   write(6, 1000) 'Labelout(',k,') = ',trim(Labelout(k))
enddo

1000 format(A9, I0, A4, A3)

!Initialize Varout
do k = 1, NVAR
   Varout(k,:) = t(k,:)
enddo

!Save initial state to external file
call create_Eulerian_file
call save_Eulerian

!Name the initial particle file
par_file = 'ParY1_D0'

call create_Particle_file(par_file)
call       save_particles(par_file)

! Prepare forcing
! Prepare temperature for temporal interpolation
call extract_WOAtemp

! Prepare Kv
call extract_Kv

!Sinking rate
!Initialize sinking rate (UNIT: m/s !):
ww(:,:) = 0d0
do k = 0,nlev-1
   !Phytoplankton no sinking 
   !Detritus sinking rate (convert to UNIT: m/s)
   ww(k,NVsinkterms) = -wDET/dble(d_per_s) 
enddo

return
END SUBROUTINE INITIALIZE