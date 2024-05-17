! !ROUTINE: Interpolate from time-series observation to model time
subroutine time_interp(time, N, nrow, obs_time, obs_data, mod_data)
! time: model time (in seconds)
!    N: number of time-series observations in observational data
! nrow: number of vertical points in observation data
! obs_time: time of observational data (in seconds), should be a vector of length N
! obs_data: observational data, should be a vector of length N
! mod_data: model data as output
implicit none
integer, intent(in)    :: N, time, nrow
real   , intent(in)    :: obs_time(N), obs_data(nrow,N)
real   , intent(out)   :: mod_data(nrow)

integer, parameter     :: y_per_s = 31536000
integer                :: i,j
real                   :: timeMOD, rat

! Get time of the year (mod):
timeMOD = float(mod(time,y_per_s))

! Deal with the case between two years:
IF ((timeMOD .lt. obs_time(1))) THEN 

   rat = (timeMOD    -obs_time(N)+float(y_per_s))   &
       / (obs_time(1)-obs_time(N)+float(y_per_s))

   do i=1,nrow
      mod_data(i)=obs_data(i,N)*(1d0-rat)+obs_data(i,1)*rat
   enddo

ELSEIF (timeMOD .ge. obs_time(N)) THEN

   rat = (timeMOD    -obs_time(N))   &
       / (obs_time(1)-obs_time(N)+float(y_per_s))

   do i=1,nrow
      mod_data(i)=obs_data(i,N)*(1d0-rat)+obs_data(i,1)*rat
   enddo

ELSE

   do j = 1, N-1
      if ((timeMOD.lt.obs_time(j+1)) .and. (timeMOD.GE.obs_time(j))) then
         
         rat=(timeMOD-obs_time(j))/(obs_time(j+1)-obs_time(j))
    
         do i=1,nrow
            mod_data(i)=obs_data(i,j)*(1d0-rat)+obs_data(i,j+1)*rat 
         enddo
         exit
      endif
   enddo 

ENDIF
return
end subroutine time_interp

