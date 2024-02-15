MODULE NETCDF_IO
USE NETCDF
USE Grid,            only : nlev,  Z_r, Z_w
USE State_variables, only : t, NZOO, iNO3, iPC, iPN, iCHL, iZOO, iDET
USE State_variables, only : Nout, Varout, oNPP, oTEMP, oPAR, oCDiv_avg
USE State_variables, only : oCDiv_var, oTalp_cov, oALnV_cov, oTLnV_cov
USE State_variables, only : oTopt_avg, oTopt_var, oLnalpha_avg, oLnalpha_var
USE State_variables, only : p_pass, p_PHY, N_PAR, N_Pass, oN_cell, oN_ind
USE State_variables, only : N_birth, N_death, N_mutate
USE Time_setting,    only : it, Nstep, restart_step
IMPLICIT NONE

private

public :: irec_Euler, irec_Pass, irec_PHY
public :: create_Eulerian_file, write_Eulerfile
public :: Create_PHY_particlefile, write_PHY_particlefile
public :: Create_Pass_particlefile, write_Pass_particlefile
public :: read_restart, create_restart, write_restart

character (len=20), public :: par_file      = 'Par.nc'
character (len=20), public :: Passive_file  = 'Passive.nc'
character (len=8),  public ::   Euler_FNAME = 'Euler.nc'
character (len=20), public :: restart_FNAME = 'restart.nc'

character (len=2),  parameter  ::   ID_NAME = 'ID'
character (len=7),  parameter  ::   Pass_ID_NAME = 'Pass_ID'
character (len=7),  parameter  ::   PHY_ID_NAME = 'PHY_ID'
character (len=2),  parameter  ::   IZ_NAME = 'IZ'
character (len=7),  parameter  ::   Pass_IZ_NAME = 'Pass_IZ'
character (len=7),  parameter  ::   PHY_IZ_NAME = 'PHY_IZ'
character (len=1),  parameter  ::    Z_NAME = 'Z'
character (len=7),  parameter  ::   Pass_Z_NAME = 'Pass_Z'
character (len=7),  parameter  ::   PHY_Z_NAME = 'PHY_Z'
character (len=3),  parameter  ::   Zr_NAME = 'Z_r'
character (len=3),  parameter  ::   Zw_NAME = 'Z_w'
character (len=4),  parameter  :: hr_NAME = 'Hour'
character (len=3),  parameter  ::  DAY_NAME = 'Day'
character (len=3),  parameter  ::  DOY_NAME = 'DOY'
character (len=4),  parameter  ::  timestep_NAME = 'Step'   !For restart file
character (len=6),  parameter  ::  N_Pass_NAME   = 'N_Pass' !For restart file
character (len=5),  parameter  ::   N_PHY_NAME   = 'N_PHY'  !For restart file
character (len=3),  parameter  ::    NO3_NAME = 'NO3'
character (len=5),  parameter  ::  P_NO3_NAME = 'p_NO3'
character (len=2),  parameter  ::    PC_NAME = 'PC'
character (len=2),  parameter  ::    PN_NAME = 'PN'
character (len=3),  parameter  ::  CHL_NAME = 'CHL'
character (len=3),  parameter  ::  ZOO_NAME = 'ZOO'
character (len=3),  parameter  ::  DET_NAME = 'DET'
character (len=3),  parameter  ::  NPP_NAME = 'NPP'
character (len=2),  parameter  ::  Kv_NAME = 'Kv'
character (len=4),  parameter  ::  Temp_NAME = 'Temp'
character (len=5),  parameter  ::  N_ind_NAME = 'N_ind'
character (len=6),  parameter  ::  N_cell_NAME = 'N_cell'
character (len=3),  parameter  ::  PAR_NAME = 'PAR'
character (len=4),  parameter  ::  CDiv_NAME = 'CDiv'
character (len=4),  parameter  ::  Topt_NAME = 'Topt'
character (len=5),  parameter  ::  alp_NAME = 'alpha'
character (len=8),  parameter  ::  CDiv_avg_NAME = 'CDiv_avg'
character (len=8),  parameter  ::  CDiv_var_NAME = 'CDiv_var'
character (len=8),  parameter  ::  Topt_avg_NAME = 'Topt_avg'
character (len=8),  parameter  ::  Topt_var_NAME = 'Topt_var'
character (len=11), parameter  ::  Lnalpha_avg_NAME = 'Lnalpha_avg'
character (len=11), parameter  ::  Lnalpha_var_NAME = 'Lnalpha_var'
character (len=8),  parameter  ::  TLnV_cov_NAME = 'TLnV_cov'
character (len=8),  parameter  ::  Talp_cov_NAME = 'Talp_cov'
character (len=8),  parameter  ::  ALnV_cov_NAME = 'ALnV_cov'
character (len=8),  parameter  ::  N_birth_NAME = 'N_birth'
character (len=8),  parameter  ::  N_death_NAME = 'N_death'
character (len=8),  parameter  ::  N_mutate_NAME = 'N_mutate'

!Number of dimensions of particles
integer,  parameter :: NDIM_PARTICLE =  2

character (len=5),  parameter  ::  UNITS = 'units'

INTEGER :: Temp_varid, PAR_varid, Kv_varid, NPP_varid
INTEGER :: N_Cell_varid, N_ind_varid, N_Pass_varid, N_PHY_varid
INTEGER :: N_death_varid, N_birth_varid, N_mutate_varid
INTEGER :: CDiv_avg_varid, CDiv_var_varid, Topt_avg_varid, Topt_var_varid
INTEGER :: Lnalpha_avg_varid, Lnalpha_var_varid
INTEGER :: TLnV_cov_varid, Talp_cov_varid, ALnV_cov_varid
INTEGER :: Zr_varid, Zw_varid, step_varid, DAY_varid, DOY_varid, hr_varid
INTEGER :: Z_varid, IZ_varid, ID_varid, P_PAR_varid, P_temp_varid, P_NO3_varid
INTEGER :: C_varid, N_varid, P_CHL_varid, P_num_varid, CDiv_varid, Topt_varid
INTEGER :: alpha_varid, CHL_varid, DET_varid, NO3_varid, PC_varid
INTEGER :: PN_varid, ZOO_varid
INTEGER :: Pass_ID_varid, Pass_IZ_varid, Pass_Z_varid
INTEGER :: PHY_ID_varid, PHY_IZ_varid, PHY_Z_varid

!Records of Euler, Passive, and PHY files
integer :: irec_Euler = 0
integer :: irec_Pass = 0
integer :: irec_PHY = 0

CONTAINS

!Write current status onto restart.nc
SUBROUTINE create_restart
IMPLICIT NONE
!Number of Dimensions for normal tracers
integer,  parameter :: NDIM_res =  1
integer,  parameter :: NDIM_ZOO_res = 2
character (len=10), parameter :: Pass_NAME = 'Pass'
character (len=10), parameter :: PHY_NAME = 'PHY'
character (len=8)    ::  date

INTEGER :: ncid, rec_dimid, NZr_dimid, NZOO_dimid, Pass_dimid, PHY_dimid
INTEGER :: N_Pass_dimid, N_PHY_dimid
INTEGER :: dimid_r(NDIM_res),  dimid_zoo(NDIM_ZOO_res)
INTEGER :: dimid_Pass(NDIM_res), dimid_PHY(NDIM_res)

! Create restart file (nc file)
CALL check(nf90_create(restart_FNAME, nf90_clobber, ncid) )
 
! Define the dimensions. 
CALL check(nf90_def_dim(ncid, timestep_NAME, 1, rec_dimid))
CALL check(nf90_def_dim(ncid, N_Pass_NAME,   1, N_Pass_dimid)) !Dimension for number of particles
CALL check(nf90_def_dim(ncid, N_PHY_NAME,    1, N_PHY_dimid)) !Dimension for number of particles
CALL check(nf90_def_dim(ncid, Zr_NAME,  nlev,   NZr_dimid))
CALL check(nf90_def_dim(ncid, ZOO_NAME, NZOO,   NZOO_dimid))
CALL check(nf90_def_dim(ncid, Pass_NAME,N_Pass, Pass_dimid))
CALL check(nf90_def_dim(ncid, PHY_NAME, N_PAR,  PHY_dimid))

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
dimid_r   = (/ NZr_dimid/)
dimid_zoo = (/NZOO_dimid, NZr_dimid/)
dimid_Pass= (/ Pass_dimid/)
dimid_PHY = (/ PHY_dimid/)

! Define the variables in the restart file
CALL check(nf90_def_var(ncid, timestep_NAME, NF90_INT, rec_dimid, step_varid) )
CALL check(nf90_def_var(ncid, N_Pass_NAME, NF90_REAL,N_Pass_dimid, N_Pass_varid) )
CALL check(nf90_def_var(ncid, N_PHY_NAME, NF90_REAL,N_PHY_dimid, N_PHY_varid) )

! Define the netCDF variables for variables that need to be saved into the Euler.nc
CALL check( nf90_def_var(ncid, NO3_NAME, NF90_REAL, dimid_r, NO3_varid) )
CALL check( nf90_def_var(ncid, DET_NAME, NF90_REAL, dimid_r, DET_varid) )
CALL check( nf90_def_var(ncid, ZOO_NAME, NF90_REAL, dimid_zoo, ZOO_varid) )

! Define the netCDF variables for passive particle variables 
CALL check( nf90_def_var(ncid, Pass_ID_NAME, NF90_INT, dimid_Pass, Pass_ID_varid) )
CALL check( nf90_def_var(ncid, Pass_IZ_NAME, NF90_INT, dimid_Pass, Pass_IZ_varid) )
CALL check( nf90_def_var(ncid, Pass_Z_NAME,  NF90_REAL,dimid_Pass, Pass_Z_varid) )

! Define the netCDF variables for PHY particle variables 
CALL check( nf90_def_var(ncid, PHY_ID_NAME, NF90_INT, dimid_PHY, PHY_ID_varid) )
CALL check( nf90_def_var(ncid, PHY_IZ_NAME, NF90_INT, dimid_PHY, PHY_IZ_varid) )
CALL check( nf90_def_var(ncid, PHY_Z_NAME, NF90_REAL, dimid_PHY, PHY_Z_varid) )
CALL check( nf90_def_var(ncid, PAR_NAME, NF90_REAL, dimid_PHY, P_PAR_varid) )
CALL check( nf90_def_var(ncid, Temp_NAME, NF90_REAL, dimid_PHY, P_Temp_varid) )
CALL check( nf90_def_var(ncid, P_NO3_NAME, NF90_REAL, dimid_PHY, P_NO3_varid) )
CALL check( nf90_def_var(ncid, CHL_NAME, NF90_REAL, dimid_PHY, P_CHL_varid) )
CALL check( nf90_def_var(ncid, PC_NAME, NF90_REAL, dimid_PHY, C_varid) )
CALL check( nf90_def_var(ncid, PN_NAME, NF90_REAL, dimid_PHY, N_varid) )
CALL check( nf90_def_var(ncid, CDiv_NAME, NF90_REAL, dimid_PHY, CDiv_varid) )
CALL check( nf90_def_var(ncid, Topt_NAME, NF90_REAL, dimid_PHY, Topt_varid) )
CALL check( nf90_def_var(ncid, alp_NAME, NF90_REAL, dimid_PHY, alpha_varid) )
CALL check( nf90_def_var(ncid, N_cell_NAME, NF90_REAL, dimid_PHY, P_num_varid) )

! Assign units attributes to the netCDF variables.
! Get the date of the nc file
CALL date_and_time(DATE = date)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))
CALL check( nf90_put_att(ncid, step_varid,  UNITS, 'Number'))
CALL check( nf90_put_att(ncid, N_Pass_varid,UNITS, 'Number of passive particles'))
CALL check( nf90_put_att(ncid, N_PHY_varid, UNITS, 'Number of phytoplankton particles'))
CALL check( nf90_put_att(ncid, NO3_varid,   UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, DET_varid,   UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, ZOO_varid,   UNITS, 'mmol N m-3'))

CALL check( nf90_put_att(ncid, Pass_Z_varid,  UNITS, 'm'))
CALL check( nf90_put_att(ncid, PHY_Z_varid,   UNITS, 'm'))
CALL check( nf90_put_att(ncid, Pass_IZ_varid, UNITS, 'no. of grid'))
CALL check( nf90_put_att(ncid, PHY_IZ_varid,  UNITS, 'no. of grid'))
CALL check( nf90_put_att(ncid, Pass_ID_varid, UNITS, 'No'))
CALL check( nf90_put_att(ncid, PHY_ID_varid,  UNITS, 'No'))
CALL check( nf90_put_att(ncid, C_varid,       UNITS, 'pmol C cell-1'))
CALL check( nf90_put_att(ncid, N_varid,       UNITS, 'pmol N cell-1'))
CALL check( nf90_put_att(ncid, Topt_varid,    UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, P_PAR_varid,   UNITS, 'W m-2'))
CALL check( nf90_put_att(ncid, P_num_varid,   UNITS, 'cells ind-1'))
CALL check( nf90_put_att(ncid, P_CHL_varid,   UNITS, 'pg Chl cell-1'))
CALL check( nf90_put_att(ncid, P_NO3_varid,   UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, P_Temp_varid,  UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, alpha_varid,   UNITS, 'log((W m-2)-1 (gChl molC)-1 d-1)'))
CALL check( nf90_put_att(ncid, CDiv_varid,    UNITS, 'pmol C cell-1'))

! End define mode.
CALL check( nf90_enddef(ncid) )
  
CALL check( nf90_close(ncid) )
return
END subroutine create_restart

!!!Save data into restart file
SUBROUTINE write_restart(current_step)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: current_step  !The current time step to be written
INTEGER              :: i = 0
INTEGER              :: ncid = 0
REAL                 :: cff(NZOO, nlev) = 0.

!Open the nc file for writing
CALL check(NF90_OPEN(restart_FNAME, NF90_WRITE, ncid))

!Check Eulerian field ID
CALL check(NF90_INQ_VARID(ncid, timestep_NAME, step_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, N_Pass_NAME,   N_Pass_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, N_PHY_NAME,    N_PHY_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, NO3_NAME, NO3_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DET_NAME, DET_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ZOO_NAME, ZOO_varid))    ! get variable IDs

!Check Lagrangian ID
!!Check passive particles
CALL check(NF90_INQ_VARID(ncid, Pass_ID_NAME, Pass_ID_varid))
CALL check(NF90_INQ_VARID(ncid, Pass_IZ_NAME, Pass_IZ_varid))
CALL check(NF90_INQ_VARID(ncid, Pass_Z_NAME,  Pass_Z_varid ))

!!Check phyto. particles
CALL check(NF90_INQ_VARID(ncid, PHY_ID_NAME, PHY_ID_varid))
CALL check(NF90_INQ_VARID(ncid, PHY_IZ_NAME, PHY_IZ_varid))
CALL check(NF90_INQ_VARID(ncid, PHY_Z_NAME,  PHY_Z_varid ))
CALL check(NF90_INQ_VARID(ncid, PAR_NAME,    P_PAR_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Temp_NAME,   P_Temp_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, p_NO3_NAME,  P_NO3_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PC_NAME,     C_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PN_NAME,     N_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, CHL_NAME,    P_CHL_varid)) 
CALL check(NF90_INQ_VARID(ncid, CDiv_NAME,   CDiv_varid)) 
CALL check(NF90_INQ_VARID(ncid, Topt_NAME,   Topt_varid)) 
CALL check(NF90_INQ_VARID(ncid, alp_NAME,    alpha_varid)) 
CALL check(NF90_INQ_VARID(ncid, N_cell_NAME, P_num_varid)) 

!Add data into the restart.nc

!!Add Eulerian data
CALL check(NF90_PUT_VAR(ncid, step_varid, current_step, start = (/1/)))
CALL check(NF90_PUT_VAR(ncid, N_Pass_varid,N_Pass, start = (/1/)))
CALL check(NF90_PUT_VAR(ncid, N_PHY_varid, N_PAR,   start = (/1/)))
CALL check(NF90_PUT_VAR(ncid, NO3_varid, t(iNO3,:), start=[1], count=[nlev]))
CALL check(NF90_PUT_VAR(ncid, DET_varid, t(iDET,:), start=[1], count=[nlev]))

!save zooplankton into a temporary matrix                                                   
do i = 1, NZOO                                                
  cff(i,:) = t(iZOO(i),:)
enddo

CALL check(NF90_PUT_VAR(ncid, ZOO_varid, cff, start=[1,1], count=[NZOO,nlev]))

!!Add lagrangian data
!!!Add passive particles
DO i = 1, N_Pass
   CALL check(NF90_PUT_VAR(ncid, Pass_ID_varid, p_Pass(i)%ID, start=[i]))
   CALL check(NF90_PUT_VAR(ncid, Pass_IZ_varid, p_Pass(i)%iz, start=[i]))
   CALL check(NF90_PUT_VAR(ncid, Pass_Z_varid,  p_Pass(i)%rz, start=[i]))
ENDDO

!!!Add Phy particles
DO i = 1, N_PAR
   CALL check(NF90_PUT_VAR(ncid, PHY_ID_varid, p_PHY(i)%ID,         start=[i]))
   CALL check(NF90_PUT_VAR(ncid, PHY_IZ_varid, p_PHY(i)%iz,         start=[i]))
   CALL check(NF90_PUT_VAR(ncid, PHY_Z_varid,  p_PHY(i)%rz,         start=[i]))
   CALL check(NF90_PUT_VAR(ncid, p_PAR_varid,  p_PHY(i)%PAR,        start=[i]))
   CALL check(NF90_PUT_VAR(ncid, p_temp_varid, p_PHY(i)%temp,       start=[i]))
   CALL check(NF90_PUT_VAR(ncid, p_NO3_varid,  p_PHY(i)%NO3,        start=[i]))
   CALL check(NF90_PUT_VAR(ncid, C_varid,      p_PHY(i)%C,          start=[i]))
   CALL check(NF90_PUT_VAR(ncid, N_varid,      p_PHY(i)%N,          start=[i]))
   CALL check(NF90_PUT_VAR(ncid, P_CHL_varid,  p_PHY(i)%CHL,        start=[i]))
   CALL check(NF90_PUT_VAR(ncid, P_num_varid,  p_PHY(i)%num,        start=[i]))
   CALL check(NF90_PUT_VAR(ncid, Topt_varid,   p_PHY(i)%Topt,       start=[i]))
   CALL check(NF90_PUT_VAR(ncid, CDiv_varid,   p_PHY(i)%CDiv,       start=[i]))
   CALL check(NF90_PUT_VAR(ncid, alpha_varid,  p_PHY(i)%LnalphaChl, start=[i]))
ENDDO

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))

!Reset N_birth, N_death, and N_mutate
N_birth(:) = 0
N_death(:) = 0
N_mutate(:)= 0

RETURN
END SUBROUTINE write_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Create Eulerian file for saving Eulerian outputs
SUBROUTINE create_Eulerian_file
IMPLICIT NONE
CHARACTER (LEN=10), PARAMETER :: REC_NAME = 'Time'
CHARACTER (LEN=10), PARAMETER :: UNIT_dist = 'm'

!Number of Dimensions for normal tracers
integer,  parameter :: NDIMS    = 2
integer,  parameter :: NDIM_ZOO = 3
character (len=8)   :: date

INTEGER :: ncid, rec_dimid, NZr_dimid, NZw_dimid, NZOO_dimid
INTEGER :: dimid_r(NDIMS), dimid_w(NDIMS), dimid_zoo(NDIM_ZOO)

! Create Euler file (nc file)
CALL check(nf90_create(Euler_FNAME, nf90_clobber, ncid) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed.
CALL check(nf90_def_dim(ncid, Zr_NAME,  nlev,   NZr_dimid))
CALL check(nf90_def_dim(ncid, Zw_NAME,  nlev+1, NZw_dimid))
CALL check(nf90_def_dim(ncid, ZOO_NAME, NZOO,   NZOO_dimid))
CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
dimid_r = (/ NZr_dimid, rec_dimid /)
dimid_w = (/ NZw_dimid, rec_dimid /)
dimid_zoo = (/NZOO_dimid, NZr_dimid, rec_dimid /)

! Define the variables of Eulerian fields
CALL check(nf90_def_var(ncid, Zr_NAME,  NF90_REAL, NZr_dimid, Zr_varid) )
CALL check(nf90_def_var(ncid, Zw_NAME,  NF90_REAL, NZw_dimid, Zw_varid) )
CALL check(nf90_def_var(ncid, DAY_NAME, NF90_INT, rec_dimid, DAY_varid) )
CALL check(nf90_def_var(ncid, hr_NAME, NF90_INT, rec_dimid, hr_varid) )

! Define the netCDF variables for variables that need to be saved into the Euler.nc
CALL check( nf90_def_var(ncid, Temp_NAME,   NF90_REAL, dimid_r, Temp_varid) )
CALL check( nf90_def_var(ncid, PAR_NAME,    NF90_REAL, dimid_r, PAR_varid) )
CALL check( nf90_def_var(ncid, N_cell_NAME, NF90_REAL, dimid_r, N_cell_varid) )
CALL check( nf90_def_var(ncid, N_ind_NAME,  NF90_REAL, dimid_r, N_ind_varid) )
CALL check( nf90_def_var(ncid, Kv_NAME,     NF90_REAL, dimid_w, Kv_varid) )
CALL check( nf90_def_var(ncid, NPP_NAME,    NF90_REAL, dimid_r, NPP_varid) )
CALL check( nf90_def_var(ncid, NO3_NAME,    NF90_REAL, dimid_r, NO3_varid) )
CALL check( nf90_def_var(ncid, PC_NAME,     NF90_REAL, dimid_r, PC_varid) )
CALL check( nf90_def_var(ncid, PN_NAME,     NF90_REAL, dimid_r, PN_varid) )
CALL check( nf90_def_var(ncid, CHL_NAME,    NF90_REAL, dimid_r, CHL_varid) )
CALL check( nf90_def_var(ncid, DET_NAME,    NF90_REAL, dimid_r, DET_varid) )
CALL check( nf90_def_var(ncid, ZOO_NAME,   NF90_REAL, dimid_zoo, ZOO_varid) )
CALL check( nf90_def_var(ncid, CDiv_avg_NAME, NF90_REAL, dimid_r, CDiv_avg_varid) )
CALL check( nf90_def_var(ncid, CDiv_var_NAME, NF90_REAL, dimid_r, CDiv_var_varid) )
CALL check( nf90_def_var(ncid, Topt_avg_NAME, NF90_REAL, dimid_r, Topt_avg_varid) )
CALL check( nf90_def_var(ncid, Topt_var_NAME, NF90_REAL, dimid_r, Topt_var_varid) )
CALL check( nf90_def_var(ncid, Lnalpha_avg_NAME, NF90_REAL, dimid_r, Lnalpha_avg_varid) )
CALL check( nf90_def_var(ncid, Lnalpha_var_NAME, NF90_REAL, dimid_r, Lnalpha_var_varid) )
CALL check( nf90_def_var(ncid, TLnV_cov_NAME,  NF90_REAL, dimid_r, TLnV_cov_varid) )
CALL check( nf90_def_var(ncid, Talp_cov_NAME,  NF90_REAL, dimid_r, Talp_cov_varid) )
CALL check( nf90_def_var(ncid, ALnV_cov_NAME,  NF90_REAL, dimid_r, ALnV_cov_varid) )
CALL check( nf90_def_var(ncid, N_birth_NAME, NF90_INT, dimid_r, N_birth_varid) )
CALL check( nf90_def_var(ncid, N_death_NAME, NF90_INT, dimid_r, N_death_varid) )
CALL check( nf90_def_var(ncid, N_mutate_NAME, NF90_INT, dimid_r, N_mutate_varid) )

! Assign units attributes to the netCDF variables.
CALL date_and_time(DATE=date)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))
CALL check( nf90_put_att(ncid,   Zr_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid,   Zw_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, DAY_varid, UNITS, 'days'))
CALL check( nf90_put_att(ncid, hr_varid, UNITS, 'hour'))
CALL check( nf90_put_att(ncid, NPP_varid, UNITS, 'mg C m-3 d-1'))
CALL check( nf90_put_att(ncid, Kv_varid, UNITS, 'm^2 s-1'))
CALL check( nf90_put_att(ncid, PAR_varid, UNITS, 'W m-2'))
CALL check( nf90_put_att(ncid, Temp_varid, UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, NO3_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, PN_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, PC_varid, UNITS, 'mmol C m-3'))
CALL check( nf90_put_att(ncid, CHL_varid, UNITS, 'mg Chl m-3'))
CALL check( nf90_put_att(ncid, DET_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, ZOO_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, N_Cell_varid, UNITS, 'No. cells m-3'))
CALL check( nf90_put_att(ncid, N_Ind_varid, UNITS, 'No. inds m-3'))
CALL check( nf90_put_att(ncid, CDiv_avg_varid, UNITS, 'log(pmol C cell-1)'))
CALL check( nf90_put_att(ncid, CDiv_var_varid, UNITS, '(log(pmol C cell-1))^2'))
CALL check( nf90_put_att(ncid, Topt_avg_varid, UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, Topt_var_varid, UNITS, 'ºC^2'))
CALL check( nf90_put_att(ncid, Lnalpha_avg_varid, UNITS, 'ln (W m-2)-1 (g Chl mol C)-1 d-1'))
CALL check( nf90_put_att(ncid, Lnalpha_var_varid, UNITS, '(ln (W m-2)-1 (g Chl mol C)-1 d-1)^2'))
CALL check( nf90_put_att(ncid, TLnV_cov_varid, UNITS, 'ºC*(ln µm^3 cell-1)'))
CALL check( nf90_put_att(ncid, Talp_cov_varid, UNITS, 'ºC*(ln (W m-2)-1 (g Chl mol C)-1 d-1)'))
CALL check( nf90_put_att(ncid, ALnV_cov_varid, UNITS, '(ln µm^3 cell-1)*(ln (W m-2)-1 (g Chl mol C)-1 d-1)'))
CALL check( nf90_put_att(ncid, N_birth_varid, UNITS, 'No. births)'))
CALL check( nf90_put_att(ncid, N_death_varid, UNITS, 'No. deaths)'))
CALL check( nf90_put_att(ncid, N_mutate_varid, UNITS, 'No. mutations)'))

! End define mode.
CALL check( nf90_enddef(ncid) )
  
! Write the Z coordinate data
CALL check( nf90_put_var(ncid, Zr_varid, Z_r) )
CALL check( nf90_put_var(ncid, Zw_varid, Z_w) )
CALL check( nf90_close(ncid) )
return

END SUBROUTINE create_Eulerian_file

SUBROUTINE write_Eulerfile(rec, day, hour)
use Forcing, only : Kv
use state_variables, only : N_birth, N_mutate, N_death
IMPLICIT NONE
INTEGER, INTENT(IN)  :: rec  !The time index to be written
INTEGER, INTENT(IN)  :: day 
INTEGER, INTENT(IN)  :: hour 
INTEGER              :: i, j
INTEGER              :: ncid = 0
real :: cff(NZOO, nlev) = 0.

!Open the nc file for writing
CALL check(NF90_OPEN(Euler_FNAME, NF90_WRITE, ncid))

CALL check(NF90_INQ_VARID(ncid, N_cell_NAME, N_cell_varid))   ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, N_ind_NAME,  N_ind_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Kv_NAME,  Kv_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Temp_NAME,Temp_varid))   ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DAY_NAME, DAY_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Hr_NAME,  Hr_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PAR_NAME, PAR_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, NO3_NAME, NO3_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PC_NAME, PC_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PN_NAME, PN_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, CHL_NAME, CHL_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DET_NAME, DET_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ZOO_NAME, ZOO_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, CDiv_avg_NAME, CDiv_avg_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, CDiv_var_NAME, CDiv_var_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Topt_avg_NAME, Topt_avg_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Topt_var_NAME, Topt_var_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Lnalpha_avg_NAME, Lnalpha_avg_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Lnalpha_var_NAME, Lnalpha_var_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, TLnV_cov_NAME, TLnV_cov_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Talp_cov_NAME, Talp_cov_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ALnV_cov_NAME, ALnV_cov_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, N_birth_NAME, N_birth_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, N_death_NAME, N_death_varid))          ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, N_mutate_NAME, N_mutate_varid))          ! get variable IDs

!Add data into the Euler.nc
CALL check(NF90_PUT_VAR(ncid, DAY_varid, day, start = (/rec/)))

CALL check(NF90_PUT_VAR(ncid, Hr_varid, hour, start = (/rec/)))

CALL check(NF90_PUT_VAR(ncid, Kv_varid, Kv(:), start=[1,rec],           &
                                               count=[1+nlev,1]))

CALL check(NF90_PUT_VAR(ncid, PAR_varid, Varout(oPAR,:),start=[1,rec],  &
                                                        count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Temp_varid,Varout(oTemp,:),start=[1,rec],  &
                                                         count=[nlev,1]))
 
CALL check(NF90_PUT_VAR(ncid, NO3_varid, t(iNO3,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, PC_varid, t(iPC,:),  start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, PN_varid, t(iPN,:),  start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, CHL_varid, t(iCHL,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, DET_varid, t(iDET,:),start=[1,rec],  &
                                                   count=[nlev,1]))

!save zooplankton into a temporary matrix                                                   
do i = 1, NZOO                                                
  cff(i,:) = t(iZOO(i),:)
enddo

CALL check(NF90_PUT_VAR(ncid, ZOO_varid, cff, start=[1,1,rec],  &
                                              count=[NZOO,nlev,1]))

CALL check(NF90_PUT_VAR(ncid, NPP_varid, Varout(oNPP,:),start=[1,rec],  &
                                                        count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, N_cell_varid, Varout(oN_cell,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, N_ind_varid, Varout(oN_ind,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, TLnV_cov_varid, Varout(oTLnV_cov,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Talp_cov_varid, Varout(oTalp_cov,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, ALnV_cov_varid, Varout(oALnV_cov,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Topt_avg_varid, Varout(oTopt_avg,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Topt_var_varid, Varout(oTopt_var,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, CDiv_avg_varid, Varout(oCDiv_avg,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, CDiv_var_varid, Varout(oCDiv_var,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Lnalpha_avg_varid, Varout(oLnalpha_avg,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Lnalpha_var_varid, Varout(oLnalpha_var,:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, N_birth_varid, N_birth(:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, N_death_varid, N_death(:),start=[1,rec],  &
                                                   count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, N_mutate_varid, N_mutate(:),start=[1,rec],  &
                                                   count=[nlev,1]))

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))

!Reset N_birth, N_death, and N_mutate
N_birth(:) = 0
N_death(:) = 0
N_mutate(:)= 0

RETURN
END SUBROUTINE write_Eulerfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Create_Pass_particlefile(fname)
IMPLICIT NONE
character (len=*), intent(in) :: fname
character (len=10), parameter ::     PAT_NAME = 'ParticleID'
character (len=10), parameter ::     REC_NAME = 'Time'
character (len=10), parameter ::    UNIT_dist = 'm'
character (len=8) ::  date

integer :: ncid, rec_dimid, PAT_dimid
integer :: dimids(NDIM_PARTICLE)

! Create particle file (nc file)
CALL check(nf90_create(FNAME, nf90_clobber, ncid) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
CALL check(nf90_def_dim(ncid, PAT_NAME, N_Pass, PAT_dimid) )
CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

CALL check(nf90_def_var(ncid, DOY_NAME, NF90_INT, rec_dimid, DOY_varid) )
CALL check(nf90_def_var(ncid,   Hr_NAME, NF90_INT, rec_dimid, Hr_varid) )

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
dimids = (/ PAT_dimid, rec_dimid /)

! Define the netCDF variables for the position data.
CALL check( nf90_def_var(ncid, ID_NAME, NF90_INT, dimids, ID_varid) )
CALL check( nf90_def_var(ncid, IZ_NAME, NF90_INT, dimids, IZ_varid) )
CALL check( nf90_def_var(ncid, Z_NAME, NF90_REAL, dimids, Z_varid) )

! Assign units attributes to the netCDF variables.
! Get the date of the nc file
CALL date_and_time(DATE = date)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))

CALL check( nf90_put_att(ncid, Z_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, IZ_varid,  UNITS, 'no. of grid'))
CALL check( nf90_put_att(ncid, ID_varid,  UNITS, 'No'))
CALL check( nf90_put_att(ncid, DOY_varid, UNITS, 'Date of the Year'))
CALL check( nf90_put_att(ncid, hr_varid, UNITS, 'hour'))

! End define mode.
CALL check( nf90_enddef(ncid) )
CALL check( nf90_close(ncid) )
RETURN
END SUBROUTINE Create_Pass_particlefile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write_Pass_particlefile(fname, rec, DOY, hour)
IMPLICIT NONE
character (len=*), intent(in) :: fname
INTEGER, INTENT(IN)  :: rec  !The time index to be written
INTEGER, INTENT(IN)  :: DOY 
INTEGER, INTENT(IN)  :: hour
INTEGER :: start(NDIM_PARTICLE)
INTEGER :: ncid = 0
INTEGER :: j = 0

! These settings tell netcdf to write one timestep of data.

!Open the nc file for writing
CALL check(NF90_OPEN(FNAME, NF90_WRITE, ncid))

! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ID_NAME, ID_varid))    
CALL check(NF90_INQ_VARID(ncid, IZ_NAME, IZ_varid))    
CALL check(NF90_INQ_VARID(ncid, Z_NAME,     Z_varid))     
CALL check(NF90_INQ_VARID(ncid, DOY_NAME, DOY_varid))
CALL check(NF90_INQ_VARID(ncid, hr_NAME, hr_varid))

CALL check(NF90_PUT_VAR(ncid, DOY_varid, DOY, start=(/rec/))) !Add data into the nc file

DO j = 1, N_Pass
   start = (/j, rec/)
   CALL check(NF90_PUT_VAR(ncid, ID_varid, p_Pass(j)%ID, start=start))
   CALL check(NF90_PUT_VAR(ncid, IZ_varid, p_Pass(j)%iz, start=start))
   CALL check(NF90_PUT_VAR(ncid, Z_varid,  p_Pass(j)%rz, start=start))
ENDDO

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))
return
END SUBROUTINE write_Pass_particlefile

!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Create_PHY_particlefile(fname)
IMPLICIT NONE
character (len=*), intent(in) :: fname
character (len=10), parameter ::     PAT_NAME = 'ParticleID'
character (len=10), parameter ::     REC_NAME = 'Time'
character (len=10), parameter ::    UNIT_dist = 'm'
character (len=8) ::  date

integer :: ncid, rec_dimid, PAT_dimid
integer :: dimids(NDIM_PARTICLE)

! Create particle file (nc file)
CALL check(nf90_create(FNAME, nf90_clobber, ncid) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
CALL check(nf90_def_dim(ncid, PAT_NAME, N_PAR, PAT_dimid) )
CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

! Define the variables of particle group and ID. 
! Ordinarily we would need to provide an array of dimension IDs for each variable's dimensions, but
! since group and particle variables only have one dimension, we can
! simply provide the address of that dimension ID
CALL check(nf90_def_var(ncid, DOY_NAME, NF90_INT, rec_dimid, DOY_varid) )
CALL check(nf90_def_var(ncid,  Hr_NAME, NF90_INT, rec_dimid, Hr_varid) )

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
dimids = (/ PAT_dimid, rec_dimid /)

! Define the netCDF variables for the position data.
CALL check( nf90_def_var(ncid, ID_NAME, NF90_INT, dimids, ID_varid) )
CALL check( nf90_def_var(ncid, IZ_NAME, NF90_INT, dimids, IZ_varid) )
CALL check( nf90_def_var(ncid, Z_NAME, NF90_REAL, dimids, Z_varid) )
CALL check( nf90_def_var(ncid, PAR_NAME, NF90_REAL, dimids, P_PAR_varid) )
CALL check( nf90_def_var(ncid, Temp_NAME, NF90_REAL, dimids, P_Temp_varid) )
CALL check( nf90_def_var(ncid, NO3_NAME, NF90_REAL, dimids, P_NO3_varid) )
CALL check( nf90_def_var(ncid, CHL_NAME, NF90_REAL, dimids, P_CHL_varid) )
CALL check( nf90_def_var(ncid, PC_NAME, NF90_REAL, dimids, C_varid) )
CALL check( nf90_def_var(ncid, CDiv_NAME, NF90_REAL, dimids, CDiv_varid) )
CALL check( nf90_def_var(ncid, Topt_NAME, NF90_REAL, dimids, Topt_varid) )
CALL check( nf90_def_var(ncid, alp_NAME, NF90_REAL, dimids, alpha_varid) )
CALL check( nf90_def_var(ncid, PN_NAME, NF90_REAL, dimids, N_varid) )
CALL check( nf90_def_var(ncid, N_cell_NAME, NF90_REAL, dimids, P_num_varid) )

! Assign units attributes to the netCDF variables.
! Get the date of the nc file
CALL date_and_time(DATE = date)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))

CALL check( nf90_put_att(ncid, Z_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, IZ_varid,  UNITS, 'no. of grid'))
CALL check( nf90_put_att(ncid, ID_varid,  UNITS, 'No'))
CALL check( nf90_put_att(ncid, DOY_varid, UNITS, 'Date of the Year'))
CALL check( nf90_put_att(ncid, hr_varid, UNITS, 'hour'))
CALL check( nf90_put_att(ncid, C_varid, UNITS, 'pmol C cell-1'))
CALL check( nf90_put_att(ncid, N_varid, UNITS, 'pmol N cell-1'))
CALL check( nf90_put_att(ncid, Topt_varid, UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, P_PAR_varid, UNITS, 'W m-2'))
CALL check( nf90_put_att(ncid, P_num_varid, UNITS, 'cells ind-1'))
CALL check( nf90_put_att(ncid, P_CHL_varid, UNITS, 'pg Chl cell-1'))
CALL check( nf90_put_att(ncid, P_NO3_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, P_Temp_varid, UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, alpha_varid, UNITS, 'log((W m-2)-1 (gChl molC)-1 d-1)'))
CALL check( nf90_put_att(ncid, CDiv_varid, UNITS, 'pmol C cell-1'))

! End define mode.
CALL check( nf90_enddef(ncid) )
CALL check( nf90_close(ncid) )
RETURN
END SUBROUTINE Create_PHY_particlefile
  
SUBROUTINE write_PHY_particlefile(fname, rec, DOY, hour)
IMPLICIT NONE
character (len=*), intent(in) :: fname
INTEGER, INTENT(IN)  :: rec  !The time index to be written
INTEGER, INTENT(IN)  :: DOY 
INTEGER, INTENT(IN)  :: hour
INTEGER :: start(NDIM_PARTICLE)
INTEGER :: ncid
INTEGER :: j = 0

! These settings tell netcdf to write one timestep of data.

!Open the nc file for writing
CALL check(NF90_OPEN(FNAME, NF90_WRITE, ncid))

! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ID_NAME,   ID_varid))    
CALL check(NF90_INQ_VARID(ncid, IZ_NAME,   IZ_varid))    
CALL check(NF90_INQ_VARID(ncid, Z_NAME,     Z_varid))     
CALL check(NF90_INQ_VARID(ncid, DOY_NAME,  DOY_varid))
CALL check(NF90_INQ_VARID(ncid, hr_NAME,   hr_varid))
CALL check(NF90_INQ_VARID(ncid, PAR_NAME,  P_PAR_varid))
CALL check(NF90_INQ_VARID(ncid, Temp_NAME, p_Temp_varid))
CALL check(NF90_INQ_VARID(ncid, NO3_NAME,  p_NO3_varid))
CALL check(NF90_INQ_VARID(ncid, PC_NAME,   C_varid))
CALL check(NF90_INQ_VARID(ncid, PN_NAME,   N_varid))
CALL check(NF90_INQ_VARID(ncid, CHL_NAME,  p_CHL_varid))
CALL check(NF90_INQ_VARID(ncid, N_cell_NAME, p_num_varid))
CALL check(NF90_INQ_VARID(ncid, Topt_NAME, Topt_varid))
CALL check(NF90_INQ_VARID(ncid, CDiv_NAME, CDiv_varid))
CALL check(NF90_INQ_VARID(ncid, alp_NAME, alpha_varid))

CALL check(NF90_PUT_VAR(ncid, DOY_varid, DOY, start=(/rec/))) !Add DOY  into the nc file
CALL check(NF90_PUT_VAR(ncid,  hr_varid, hour,start=(/rec/))) !Add Hour into the nc file

DO j = 1, N_PAR
   start = (/j, rec/)
   CALL check(NF90_PUT_VAR(ncid, ID_varid,     p_PHY(j)%ID,   start=start))
   CALL check(NF90_PUT_VAR(ncid, IZ_varid,     p_PHY(j)%iz,   start=start))
   CALL check(NF90_PUT_VAR(ncid, Z_varid,      p_PHY(j)%rz,   start=start))
   CALL check(NF90_PUT_VAR(ncid, p_PAR_varid,  p_PHY(j)%PAR,  start=start))
   CALL check(NF90_PUT_VAR(ncid, p_temp_varid, p_PHY(j)%temp, start=start))
   CALL check(NF90_PUT_VAR(ncid, p_NO3_varid,  p_PHY(j)%NO3,  start=start))
   CALL check(NF90_PUT_VAR(ncid, C_varid,      p_PHY(j)%C,    start=start))
   CALL check(NF90_PUT_VAR(ncid, N_varid,      p_PHY(j)%N,    start=start))
   CALL check(NF90_PUT_VAR(ncid, p_CHL_varid,  p_PHY(j)%CHL,  start=start))
   CALL check(NF90_PUT_VAR(ncid, p_num_varid,  p_PHY(j)%num,  start=start))
   CALL check(NF90_PUT_VAR(ncid, Topt_varid,   p_PHY(j)%Topt, start=start))
   CALL check(NF90_PUT_VAR(ncid, CDiv_varid,   p_PHY(j)%CDiv, start=start))
   CALL check(NF90_PUT_VAR(ncid, alpha_varid,  p_PHY(j)%LnalphaChl, start=start))
ENDDO

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))
return
END SUBROUTINE write_PHY_particlefile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_restart
IMPLICIT NONE

! This will be the netCDF ID for the file and data variable.
integer :: ncid, i

!Scratch variable for temporally storing zooplankton data
real    :: cff(NZOO, nlev) = 0.
real,    allocatable :: pffr(:)
integer, allocatable :: pffi(:)

integer :: new_N_Pass = 0
integer :: new_N_PAR  = 0
integer :: allocateStatus = 0

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.
CALL check(nf90_open(restart_FNAME, NF90_NOWRITE, NCID) )

! Get the varid of the data variable, based on its name.
CALL check(NF90_INQ_VARID(ncid, timestep_NAME, step_varid))
CALL check(NF90_INQ_VARID(ncid, N_Pass_NAME, N_Pass_varid))
CALL check(NF90_INQ_VARID(ncid, N_PHY_NAME,  N_PHY_varid))
CALL check(NF90_INQ_VARID(ncid, NO3_NAME,      NO3_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DET_NAME,      DET_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ZOO_NAME,      ZOO_varid))    ! get variable IDs

!Check Lagrangian ID
!!Check passive particles
CALL check(NF90_INQ_VARID(ncid, Pass_ID_NAME, Pass_ID_varid))
CALL check(NF90_INQ_VARID(ncid, Pass_IZ_NAME, Pass_IZ_varid))
CALL check(NF90_INQ_VARID(ncid, Pass_Z_NAME,  Pass_Z_varid ))

!!Check phyto. particles
CALL check(NF90_INQ_VARID(ncid, PHY_ID_NAME, PHY_ID_varid))
CALL check(NF90_INQ_VARID(ncid, PHY_IZ_NAME, PHY_IZ_varid))
CALL check(NF90_INQ_VARID(ncid, PHY_Z_NAME,  PHY_Z_varid ))
CALL check(NF90_INQ_VARID(ncid, PAR_NAME,    P_PAR_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Temp_NAME,   P_Temp_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, p_NO3_NAME,  P_NO3_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PC_NAME,     C_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PN_NAME,     N_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, CHL_NAME,    P_CHL_varid)) 
CALL check(NF90_INQ_VARID(ncid, CDiv_NAME,   CDiv_varid)) 
CALL check(NF90_INQ_VARID(ncid, Topt_NAME,   Topt_varid)) 
CALL check(NF90_INQ_VARID(ncid, alp_NAME,    alpha_varid)) 
CALL check(NF90_INQ_VARID(ncid, N_cell_NAME, P_num_varid)) 

! Read the data.
call check(NF90_GET_VAR(ncid, step_varid, restart_step))

if (restart_step > Nstep) then
   stop "Initial step beyond the preset total number of steps!"
endif

call check(NF90_GET_VAR(ncid, N_Pass_varid, new_N_Pass) )

if (new_N_Pass .ne. N_Pass) stop "N_Pass read from restart.nc is unequal to that from Param.nml!"

call check(NF90_GET_VAR(ncid, N_PHY_varid,  new_N_PAR) )
if (new_N_PAR .ne. N_PAR) stop "N_PAR read from restart.nc is unequal to that from Param.nml!"

! Read the Eulerian data
call check(NF90_GET_VAR(ncid, NO3_varid, t(iNO3,:), start=[1], count=[nlev]) )
call check(NF90_GET_VAR(ncid, DET_varid, t(iDET,:), start=[1], count=[nlev]) )
call check(NF90_GET_VAR(ncid, ZOO_varid, cff, start=[1,1], count=[NZOO, nlev]))

do i = 1, NZOO
   t(iZOO(i), :) = cff(i, :)
enddo

!Initialize particles
IF (.not. allocated(p_pass)) ALLOCATE(p_pass(N_Pass), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_pass ***"

IF (.not. allocated(p_PHY)) ALLOCATE(p_PHY(N_Par), stat=AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_PHY ***"
 
!Read data of passive particles
allocate(pffi(N_Pass))
call check(NF90_GET_VAR(ncid, Pass_ID_varid, pffi))

DO i = 1, N_Pass
  p_Pass(i)%ID = pffi(i)
ENDDO

call check(NF90_GET_VAR(ncid, Pass_IZ_varid, pffi))

DO i = 1, N_Pass
  p_Pass(i)%iz = pffi(i)
ENDDO

allocate(pffr(N_Pass))
call check(NF90_GET_VAR(ncid, Pass_Z_varid,  pffr))
DO i = 1, N_Pass
  p_Pass(i)%rz = pffr(i)
ENDDO

deallocate(pffi)
deallocate(pffr)

!Read data of phytoplankton particles
allocate(pffi(N_PAR))
call check(NF90_GET_VAR(ncid, PHY_ID_varid, pffi))

DO i = 1, N_PAR
  p_PHY(i)%ID = pffi(i)
ENDDO

call check(NF90_GET_VAR(ncid, PHY_IZ_varid, pffi))

DO i = 1, N_PAR
  p_PHY(i)%iz = pffi(i)
ENDDO

allocate(pffr(N_PAR))
call check(NF90_GET_VAR(ncid, PHY_Z_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%rz = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, p_PAR_varid,  pffr))
DO i = 1, N_PAR
  p_PHY(i)%PAR = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, p_temp_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%temp = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, p_NO3_varid,  pffr))
DO i = 1, N_PAR
  p_PHY(i)%NO3 = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, C_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%C = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, N_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%N = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, P_CHL_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%CHL = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, P_num_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%num = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, Topt_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%Topt = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, CDiv_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%CDiv = pffr(i)
ENDDO

call check(NF90_GET_VAR(ncid, alpha_varid, pffr))
DO i = 1, N_PAR
  p_PHY(i)%LnalphaChl = pffr(i)
ENDDO

!Finish reading data from nc file and deallocate pffi and pffr
deallocate(pffi)
deallocate(pffr)

! Close the file, freeing all resources.
call check( nf90_close(ncid) )
return

END SUBROUTINE read_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check(status)
integer, intent (in) :: status

if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine check

END MODULE
