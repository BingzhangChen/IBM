!-----------------------------------------------------------------------
subroutine tridiagonal(N,au,bu,cu,du,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix structure is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}.
! The method used here is the simplified Gauss elimination, also called
! \emph{Thomas algorithm}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N,fi,lt
   real(8), dimension(N), intent(in)   :: au, bu, cu, du
!
! !OUTPUT PARAMETERS:
   real(8), intent(out)                :: value(N)
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: i
   real(8), dimension(N)               :: ru, qu
!
!-----------------------------------------------------------------------
!BOC
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)

!   write(6,*) 'fi = ',fi
!   write(6,*) 'value(fi) = ',value(fi)

   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do

!   write(6,*) 'value(lt) = ',value(lt)

   return
end subroutine tridiagonal

