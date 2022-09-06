SUBROUTINE INITIALIZE
use params
use state_variables
use Time_setting, only: d_per_s, dtsec, dtdays, Nstep, nsave, NDay_Run
use grid,         only: Z_w, hmax
use Trait_functions,  only: PHY_ESD2C
use IO
use forcing
implicit NONE

integer            :: rc       = 0
integer            :: k        = 0
integer            :: j_       = 0
real               :: cff      = 0.d0
integer, parameter :: namlst   = 20   !Unit time for namelist files
character(len=10)  :: par_file = 'Filename'
character(len=10), parameter :: format_string = "(A3,I0)"
!==========================================================

!Namelist definition of time settings
namelist /timelist/  NDay_Run, dtsec, nsave

!Namelist definition of model choice and parameters
namelist /paramlist/ Model_ID, mu0, aI0, KN, gmax, Kp, mz,RDN, &
                     wDET, SDZoo, nu,sigma

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

do k = 1, NZOO
   t(iZOO(k),:) = 0.1d0/dble(NZOO) !Assuming an initial condition of uniform biomass among different zoo. size classes

   !Initialize zooplankton size (logESD)
   ESDZOO(k) = MinSzoo + dble(k-1)*dZOOESD

   !Make the model write out the zooplankton size
   write(6,1001) "ZOO", k, "ESD = ", exp(ESDZOO(k))

   !Compute volume of zooplankton
   VolZOO(k) = pi/6d0*exp(ESDZOO(k))**3
enddo

1001 format(A3, 1x, I0, 1x, A6, F12.3 )

!Following Verity et al. AME (1996)
t(iDET,:) = .1d0

!For lagrangian model, initialize individuals
DO k = 1, N_PAR
   p_PHY(k)%ID = k        
   p_PHY(k)%alive = .true.

   !Initialize optimal temperature (Topt) to 20 C
   p_PHY(k)%Topt = 20.d0

   !Initialize phytoplankton size from a uniform distribution between 0.8 and 60 um
   call random_number(cff)
   cff = log(0.8d0) + cff * (log(60.d0) - log(0.8d0))
   cff = exp(cff) !ESD

   p_PHY(k)%C    = PHY_ESD2C(cff)      !Unit: pmol C cell-1
   p_PHY(k)%N    = p_PHY(k)%C/106.*16. !Unit: pmol N cell-1
   p_PHY(k)%Chl  = p_PHY(k)%C* 12./50. !Unit: pgChl cell-1
   p_PHY(k)%CDiv = p_PHY(k)%C* 2d0 !Unit: pmol C cell-1

   !Initialize alphaChl to 0.1 (W m-2)-1 (gChl molC)-1 d-1
   p_PHY(k)%LnalphaChl = -2.3

   !Initialize the number of cells associated with each super-individual (assuming initial phytoplankton nitrogen is 0.1 mmol m-3)
   p_PHY(k)%num = 0.1 * hmax/dble(N_PAR)/p_PHY(k)%N * 1d9

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
Labelout(iDET) = 'DET'
Labelout(iPC)  = 'PC '
Labelout(iPN)  = 'PN '
Labelout(iCHL) = 'CHL'
Labelout(oNPP) = 'NPP'
Labelout(oTEMP)= 'TEMP'
Labelout(oPAR) = 'PAR'

do k = 1, NZOO
   write(Labelout(iZOO(k)),  format_string) 'ZOO', k
   write(Labelout( oFZ(k)),  format_string) 'FZO',k
enddo

write(6, '(a)') 'Write out the labels for output for validation:'
do k = 1, Nout
   write(6, 1000) 'Labelout(',k,') = ',trim(Labelout(k))
enddo

1000 format(A9, I0, A4, A5)

!Initialize Varout
do k = 1, NVAR
   Varout(k,:) = t(k,:)
enddo

!Save initial state to external file
call create_Eulerian_file
call save_Eulerian

!Save Kv
call create_Kv_file
call save_Kv

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
