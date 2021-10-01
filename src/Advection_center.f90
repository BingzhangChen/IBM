! !ROUTINE: Advection schemes --- grid centers\label{sec:advectionMean}
!  Original author(s): Karsten Bolding & Hans Burchard
subroutine adv_center(N,dt,h,ho,ww,Bcup,Bcdw,Yup,Ydw,method,mode,Y)
   IMPLICIT NONE
! !INPUT PARAMETERS:
!  number of vertical layers
   integer,  intent(in)                :: N
!  time step (s)
   real,     intent(in)                 :: dt
!  layer thickness (m)
   real,     intent(in)                 :: h(N)

!  old layer thickness (m)
   real,     intent(in)                 :: ho(N)

!  vertical advection speed (m/s)
   real,     intent(in)                 :: ww(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   real, intent(in)                 :: Yup

!  value of lower BC
   real,     intent(in)            :: Ydw

!  type of advection scheme
   integer,  intent(in)            :: method

!  advection mode (0: non-conservative, 1: conservative)
   integer,  intent(in)            :: mode
!
! !INPUT/OUTPUT PARAMETERS:
   real,  intent(inout)            :: Y(N)
!
! !DEFINED PARAMETERS:
   real,     parameter             :: one6th=1d0/6d0
   integer,  parameter             :: itmax=100

!  type of advection scheme
   integer,  parameter             :: UPSTREAM       = 1
   integer,  parameter             :: P1             = 2
   integer,  parameter             :: P2             = 3
   integer,  parameter             :: Superbee       = 4
   integer,  parameter             :: MUSCL          = 5
   integer,  parameter             :: P2_PDM         = 6

!  boundary condition type
!  for advection schemes
   integer,  parameter             :: flux           = 1
   integer,  parameter             :: value          = 2
   integer,  parameter             :: oneSided       = 3
   integer,  parameter             :: zeroDivergence = 4

! !LOCAL VARIABLES:
   integer                         :: i,k,it
   real                            :: x,r,Phi,limit
   real                            :: Yu,Yc,Yd
   real                            :: c,cmax
   real                            :: cu(0:N)
!-----------------------------------------------------------------------
!  initialize interface fluxes with zero
   cu   = 0.0

!  initialize maximum Courant number
   cmax = 0.0

!  compute maximum Courant number
   do k = 1,N-1
      c = abs(ww(k))*dt/(0.5*(h(k)+h(k+1)))
      if (c .gt. cmax) cmax=c
   enddo

   it   = min(itmax,int(cmax)+1)

!#ifdef DEBUG
   if (it .gt. 1) then
      write(0,*) 'In adv_center():'
      write(0,*) 'Maximum Courant number is ',cmax
      write(0,*) it,' iterations used for vertical advection'
   endif
!#endif

!  splitting loop
   do i=1,it

!     vertical loop
      do k=1,N-1

!        compute the slope ration
         if (ww(k) .gt. 0.0) then

!           compute Courant number
            c=ww(k)/float(it)*dt/(0.5*(h(k)+h(k+1)))

            if (k .gt. 1) then
               Yu=Y(k-1)                              ! upstream value
            else
               Yu=Y(k)
            end if
            Yc=Y(k  )                                 ! central value
            Yd=Y(k+1)                                 ! downstream value

!           compute slope ratio
            if (abs(Yd-Yc) .gt. 1e-10) then
               r=(Yc-Yu)/(Yd-Yc)
            else
               r=(Yc-Yu)*1.e10
            end if

!        negative speed
         else

!           compute Courant number
            c=-ww(k)/float(it)*dt/(0.5*(h(k)+h(k+1)))

            if (k .lt. N-1) then
               Yu=Y(k+2)                              ! upstream value
            else
               Yu=Y(k+1)
            end if
            Yc=Y(k+1)                                 ! central value
            Yd=Y(k  )                                 ! downstream value

!           compute slope ratio
            if (abs(Yc-Yd) .gt. 1e-10) then
               r=(Yu-Yc)/(Yc-Yd)
            else
               r=(Yu-Yc)*1.e10
            end if
         end if

!        compute the flux-factor phi
         x    =  one6th*(1d0-2d0*c)
         Phi  =  (5d-1+x)+(5d-1-x)*r

!   limit the flux according to different suggestions (Pietrzak 1998)
         select case (method)
            case (UPSTREAM)
               limit= 0.0
            case (P1)
               write(0,*) "Fatal error: P1 advection method not yet implemented, choose other method"
               stop  "adv_center.F90"
            case ((P2),(P2_PDM))
               if (method .eq. P2) then
                  limit=Phi
               else
                  limit=max(0.0,min(Phi,2d0/(1d0-c),2d0*r/(c+1.e-10)))
               end if
            case (Superbee)
               limit=max(0.0, min(1d0, 2.0*r), min(r,2.*1d0) )
            case (MUSCL)
               limit=max(0.0,min(2.*1d0,2.0*r,0.5*(1.0+r)))
            case default
               write(0,*) method
               write(0,*) 'Fatal error: unkown advection method in adv_center()'
               stop
          end select

!        compute the limited flux
         cu(k)=ww(k)*(Yc+0.5*limit*(1-c)*(Yd-Yc))

      end do

!     do the upper boundary conditions
      select case (Bcup)
      case (flux)
         cu(N) = - Yup              ! flux into the domain is positive
      case (value)
         cu(N) =  ww(N)*Yup
      case (oneSided)
         if (ww(N).ge.0.0) then
            cu(N) =  ww(N)*Y(N)
         else
            cu(N) = 0.0
         end if
      case (zeroDivergence)
         cu(N) = cu(N-1)
      case default
         write(0,*) 'Fatal error: unkown upper boundary condition type in adv_center()'
         stop
      end select

!     do the lower boundary conditions
      select case (Bcdw)
      case (flux)
         cu(0) =   Ydw               ! flux into the domain is positive
      case (value)
         cu(0) =  ww(0)*Ydw
      case (oneSided)
         if(ww(0).le.0.0) then
            cu(0) =  ww(0)*Y(1)
         else
            cu(0) = 0.0
         end if
      case (zeroDivergence)
         cu(0) = cu(1)
      case default
         write(0,*) 'Fatal error: unkown lower boundary condition type in adv_center()'
         stop
      end select

!     do the vertical advection step which will be used for prescribed
!     vertical flow velocity and for settling of suspended matter.

      if (mode .eq. 0) then ! non-conservative, water vertical velocity
         do k=1,N
            Y(k)=Y(k)-1d0/float(it)*dt*((cu(k)-cu(k-1))/        &
                 h(k)-Y(k)*(ww(k)-ww(k-1))/h(k))
         enddo
      else                ! conservative, for PHY and DET sinking
         do k=1,N
            Y(k)=Y(k)-1d0/float(it)*dt*((cu(k)-cu(k-1))/h(k))
         enddo
      end if

   end do ! end of the iteration loop
  return
end subroutine adv_center
!--------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org

