! !ROUTINE: Interpolate from observation space to model grid
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !INTERFACE:
subroutine gridinterpol(N,cols,obs_z,obs_prof,nlev_,model_z,model_prof)
!
! !DESCRIPTION:
!
!  This is a utility subroutine in which observational data, which might
!  be given on an arbitrary, but structured grid, are linearly interpolated and
!  extrapolated to the actual (moving) model grid.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   !N    : number of vertical layers in observation data
   !cols : number of profiles in obs. data
   integer,  intent(in)                :: N,cols
   REAL,     intent(in)                :: obs_z(N),obs_prof(N,cols)
   !nlev_: number of vertical layers in model
   integer,  intent(in)                :: nlev_
   REAL,     intent(in)                :: model_z(nlev_)
!
! !OUTPUT PARAMETERS:
   REAL,     intent(out)               :: model_prof(nlev_,cols)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii
   REAL                      :: rat
!-----------------------------------------------------------------------
!BOC
   model_prof(:,:) = 0d0
!  Set surface values to uppermost input value
   do i=nlev_,1,-1
      if(model_z(i) .ge. obs_z(N)) then
         do j=1,cols
            model_prof(i,j) = obs_prof(N,j)
         end do
      end if
   end do

!  Set bottom values to lowest input value
   do i=1,nlev_
      if(model_z(i) .le. obs_z(1)) then
         do j=1,cols
            model_prof(i,j) = obs_prof(1,j)
         end do
      end if
   end do

!  Interpolate inner values linearly
   do i=1, nlev_, 1
      if ((model_z(i) .lt. obs_z(N)) .and. (model_z(i) .gt. obs_z(1))) then
         ii=0
224      ii=ii+1
         !write(6,*) 'model_z(',i,') = ',model_z(i)
         !write(6,*) 'obs_z(',ii,') = ',obs_z(ii)
         if (obs_z(ii) .le. model_z(i)) goto 224

         rat=(model_z(i)-obs_z(ii-1))/(obs_z(ii)-obs_z(ii-1))

         do j=1,cols
            model_prof(i,j)=(1d0-rat)*obs_prof(ii-1,j)+rat*obs_prof(ii,j)
         enddo
      endif
   enddo
end subroutine
!--------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org

