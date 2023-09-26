Module forcing
!This module provides several analytic functions calculating temperature, surface PAR, MLD, maximal eddy diffusivity as a function of time
use Time_setting, only: d_per_s, y_per_d, s_per_h
use Grid,         only: nlev, Z_r, Hz, Z_w
use params,       only: pi
implicit none

private

!Temperature forcing
character(len=20), parameter :: temp_time_file = 'BATS_temp_time.dat'
character(len=20), parameter :: temp_file      = 'BATS_temp.dat'

!Number of time points of observed WOA temperature
integer,           parameter :: N_time_temp    = 12

!Number of vertical points of observed WOA temperature
integer,           parameter :: N_V_Temp       = 57

!Observed time of temperature (unit: second)
real                         :: obs_time_temp(N_time_temp)        = 0.d0

!Real temperature data (the first column is the depth)
real                         :: obs_temp(N_V_Temp, 1+N_time_temp) = 0.d0

!Full Temperature data interpolated to the model grid
real                         :: VTemp(nlev, N_time_temp)          = 0.d0

!Final temperature data at each model grid at the target time
real                         ::  Temp(nlev)                       = 0.d0

!Kv forcing
character(len=20), parameter :: Kv_time_file = 'BATS_Aks_time.dat'
character(len=20), parameter :: Kv0_file     = 'BATS_KV0.dat'
character(len=20), parameter :: Kv_file      = 'BATS_Aks.dat'
character(len=20), parameter :: Kvmax_file   = 'BATS_KVmax.dat'
character(len=20), parameter :: MLD_file     = 'BATS_MLD.dat'

!Number of time points of observed Kv
integer,           parameter :: N_time_Kv    = 12

!Observed time of Kv (unit: second)
real                         :: obs_time_Kv(N_time_Kv)            = 0.d0

!Surface Kv
real                         :: obs_Kv0(N_time_Kv)                = 0.d0
real                         :: Kv0(1)                            = 0.d0

!Water column Kv directly estimated from ROMS model output
integer,       parameter     :: N_obs_Kv                          = 41
real                         :: obs_Kv(N_obs_Kv, 1+N_time_Kv)     = 0.d0

!MLD
real                         :: obs_MLD(N_time_Kv)                = 0.d0
real                         :: MLD(1)                            = 0.d0

!Kvmax
real                         :: obs_Kvmax(N_time_Kv)              = 0.d0
real                         :: Kvmax(1)                          = 0.d0

!Full Kv data interpolated to the model grid
real                         :: VKv(0:nlev, N_time_Kv)          = 0.d0

!Final Kv at w points of each grid at the target time
real                         :: Kv(0:nlev)                      = 0.d0

!Background diffusivity
real, parameter :: Kbg = 3d-5 !Kbg follows Christina L&O 2005

!PAR
real                         :: PARw(0:nlev)                      = 0.d0
real                         :: PAR(   nlev)                      = 0.d0

!Final w at w points of each grid at the target time
real                         :: w(0:nlev)                         = 0.d0

! Bottom boundary condition:
integer,           parameter :: Dirichlet                         = 0
integer,           parameter :: Neumann                           = 1
integer,           parameter :: bot_bound                         = Neumann

!Calculating diel light
public :: diurnal_light, season, extract_WOAtemp, VERTICAL_LIGHT !Subroutines and functions
public :: Temp, VTemp, N_time_temp, obs_time_temp      !Scalars and vectors
public :: PARw, PAR
public :: w, Kv, VKv, N_time_Kv, obs_time_Kv, obs_MLD, obs_Kv0, obs_Kv, obs_Kvmax, MLD, Kvmax, Kv0, Kbg   !Scalars and vectors
public :: Dirichlet, Neumann, bot_bound
public :: extract_Kv

! Station latitude (Convert latitude to radians)
real, parameter :: Stn_lat = 31.67d0 *pi/180.d0  !BATS

CONTAINS

PURE REAL FUNCTION cosine(day) result(y) 
IMPLICIT NONE
real, intent(in) :: day
real, parameter  :: offset = 60.
  y   = cos(2. * pi * (day-offset) / dble(y_per_d) )
END FUNCTION cosine

pure real function season(mean, Amp, day)result(y)
IMPLICIT NONE
real, intent(in) :: mean, Amp, day

y = mean - Amp*cosine(day)/2.d0
end function season

!Function following Anderson et al. GMD 2015 to calculate PAR below the ocean surface
! Noon PAR, W m-2

pure real function noonparnow(jday, latradians)
implicit none
integer, intent(in) :: jday      !Day of year
real,    intent(in) :: latradians!Latitude in radianns

real, parameter :: solarconst = 1368.0!Solar constant (W m-2)
real, parameter :: e0    = 12.0  !partial pressure of water vapour in the atmosphere (mb)
real, parameter :: parrac= 0.43  !Fraction of solar radiation
real, parameter :: albedo= 0.04  !Ocean albedo
real, parameter :: clouds= 6d0   !cloud fraction in oktas
real            :: cfac     !Effect of clouds on atmospheric transmission
real            :: Iclear   !Short-wave irradiance at the ocean surface on a clear day
real            :: coszen    !Cosine of solar zenith angle
real            :: declin    !Solar declination angle (in radians)
real            :: zen       !solar zenith angle
real            :: Rvector, Inoon

declin = 23.45*sin(2.*pi*(284. + dble(jday))*0.00274)*pi/180.     ! solar declination angle
coszen = sin(latradians)*sin(declin)+cos(latradians)*cos(declin)  ! cosine of zenith angle
zen    = acos(coszen)*180./pi                                     !zenith angle, degrees
Rvector= 1d0/sqrt(1. + 0.033*cos(2.*pi*dble(jday)*0.00274))       ! Earth's radius vector
Iclear = solarconst*coszen**2/(Rvector**2)/(1.2*coszen+e0*(1.0+coszen)*0.001+0.0455) ! irradiance at ocean surface, clear sky
cfac   = (1d0-0.62*clouds*0.125 + 0.0019*(90d0-zen))                   ! cloud factor (atmospheric transmission)
Inoon  = Iclear*cfac*(1d0 - albedo)                                    ! noon irradiance: total solar
noonparnow = parrac*Inoon   
return
end function noonparnow

! ----------------------------------------------------------------------- #
! Calculation of day length as function of day of year and latitude       #
! ----------------------------------------------------------------------- #
pure real function FNdaylcalc(jday,latradians) result(daylnow)
implicit none
integer, intent(in) :: jday      !Day of year
real,    intent(in) :: latradians!Latitude in radianns

real                :: declin    !Solar declination angle (in radians)

declin = 23.45*sin(2.*pi*(284. + jday)*0.00274)*pi/180.      ! solar declination angle
daylnow= 2.*acos(-1. * tan(latradians)*tan(declin))*12./pi   ! hours
end function FNdaylcalc

! ----------------------------------------------------------------------- #
! Calculation of surface PAR as function of mid-day light, time of the day, and day length
! ----------------------------------------------------------------------- #
pure REAL FUNCTION surface_par(Inoon, timeofday, daylength)
IMPLICIT NONE
REAL,    INTENT(IN) :: Inoon         ! Surface par at noon Unit: W m-2
INTEGER, INTENT(IN) :: timeofday     ! Unit: second
REAL,    INTENT(IN) :: daylength     ! Unit: hours

real,    parameter  :: Tnoon = 43200d0 ! seconds at noon
real    :: Trise   !Unit: seconds
real    :: Tset    !Unit: seconds
real    :: Day_len_sec  ! Convert day length into seconds

Day_len_sec = daylength*s_per_h

Trise = Tnoon - Day_len_sec/2d0 
Tset  = Tnoon + Day_len_sec/2d0 

! Initialize surface PAR based on diel cycle
if (dble(timeofday) .le. Trise .or. dble(timeofday) .ge. Tset) then
    surface_par = 0d0
else
    surface_par = sin(pi*(dble(timeofday) - Trise)/Day_len_sec)*Inoon
endif

return
END FUNCTION surface_par

real function diurnal_light(jday,latradians, time_of_day) result(y)
IMPLICIT NONE
integer, intent(in) :: jday         !Julian date
real,    INTENT(in) :: latradians   !Latitude in radians
integer, INTENT(in) :: time_of_day  !Unit: second
real                :: Day_len = 0.d0
real                :: I_noon_ = 0.d0

if (time_of_day .eq. 0) then
   !Daylength
   Day_len = FNdaylcalc(jday, latradians)
   
   !I_noon
   I_noon_ = noonparnow(jday, latradians)
endif

!Diel light
y = surface_par(I_noon_, time_of_day, Day_len)

return
end function diurnal_light

subroutine VERTICAL_LIGHT(jday, time_of_day, Chl) 
IMPLICIT NONE
INTEGER, intent(in) :: jday         !Julian date
INTEGER, INTENT(in) :: time_of_day  !Unit: second
REAL,    INTENT(in) :: Chl(nlev)    !Vertical distributions of Chl

!Surface PAR 
real                ::  PAR0= 0.d0

PAR0 = diurnal_light(jday, Stn_lat, time_of_day)

call Calculate_PAR(PAR0, nlev, Hz, Chl, PAR, PARw)

return
END subroutine VERTICAL_LIGHT

subroutine extract_Kv
!This subroutine extracts Kv0, Kvmax and MLD from external files and should be called only once during initialization
implicit none

!Obtain time
call Readcsv(temp_time_file, 1, N_time_Kv, obs_time_Kv) 
obs_time_Kv(:) = obs_time_Kv(:)*3.d1*dble(d_per_s) !unit: second

!Obtain MLD
call Readcsv(MLD_file, 1, N_time_Kv, obs_MLD) 

!Obtain Kvmax
call Readcsv(Kvmax_file, 1, N_time_Kv, obs_Kvmax) 

!Obtain Kv0
call Readcsv(Kv0_file, 1, N_time_Kv, obs_Kv0) 

!Obtain Kv
call Readcsv(Kv_file, N_obs_Kv, N_time_Kv, obs_Kv) 

! Interpolate external Aks data:
call gridinterpol(N_obs_Kv, N_time_Kv, obs_Kv(:, 1), &
                  obs_Kv(:, 2:(N_time_Kv + 1)), nlev+1, Z_w, VKv)

end subroutine extract_Kv

subroutine extract_WOAtemp
!This subroutine extracts temperature data from WOA13 observations and should be called only once during initialization
implicit none

!Obtain time
call Readcsv(temp_time_file, 1, N_time_temp, obs_time_temp) 
obs_time_temp(:) = obs_time_temp(:)*3.d1*dble(d_per_s) !unit: second

!Obtain WOA temperature data
call Readcsv(temp_file, size(obs_Temp,1), size(obs_Temp,2), obs_Temp) 

! Interpolate external Temp data:
call gridinterpol(N_V_Temp, N_time_temp, obs_Temp(:,1),                   &
      obs_Temp(:,2:(N_time_temp+1)), nlev, Z_r, VTemp) 
 
return
end subroutine extract_WOAtemp

END MODULE
