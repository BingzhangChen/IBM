module mGf90
implicit none

private
public srand_mtGaus
contains

real function srand_mtGaus(N, mean, cvm) result(y)
implicit none

!The dimension of the multivariate vector
integer, intent(in)   :: N

real,    intent(in)   :: mean(N)  !Mean values of the multivariate distribution

real,    intent(in)   :: cvm(N,N) !Covariance matrix of the multivariate distribution
dimension             :: y(N)
integer    :: nullty  ! i/o for cholesky subroutine
integer    :: error  
integer    :: i, j
real       :: cvm_(N*(N+1)/2)
real       :: chol(N*(N+1)/2)

!Convert cvm to a vector with length N*(N+1)/2
do i = 1, N
  do j = 1,i
    cvm_(i*(i-1)/2+j) = cvm(i,j)
  enddo
enddo

call cholesky(cvm_,N,N*(N+1)/2,chol,nullty,error)

y = multiGauss(chol,mean,N)
end function srand_mtGaus

real function gasdev() result(gval)
!This function generates an array of independent Gaussian variates
implicit none
real :: FAC, R, V1, V2, X
real, save :: GSET
integer, save :: ISET = 0

IF (ISET.EQ.0) THEN
   R = 99 
   do while( R .ge. 1.0d0 )
      call random_number(V1)
      V1= 2d0*V1 - 1d0

      call random_number(V2)
      V2= 2d0*V2 - 1d0

      R = V1**2 + V2**2
   end do
   FAC  = SQRT( -2.0d0*LOG(R)/R )
   GSET = V1*FAC
   gval = V2*FAC
   ISET = 1
ELSE
   gval=GSET
   ISET=0
ENDIF

RETURN
END function gasdev

subroutine cholesky ( a, n, nn, u, nullty, ifault )
!*****************************************************************************80
!
!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
!
!  Discussion:
!
!    For a positive definite symmetric matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!  Modified:
!
!    01 February 2008
!
!  Author:
!
!    Michael Healy
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix 
!    stored by rows in lower triangular form as a one dimensional array, 
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, integer NN, the dimension of the array used to store A, 
!    which should be at least (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!    It is also the same as lower triangular matrix stored by rows (by B. Chen on 12 Aug 2019)
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY is zero,
!    the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, NN < (N*(N+1))/2.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: n
  integer ( kind = 4 ), intent(in) :: nn

  real    ( kind = 8 ), intent(in) :: a(nn)
  real    ( kind = 8 ), parameter :: eta = 1.0D-19
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ), intent(out) :: ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), intent(out) :: nullty
  real    ( kind = 8 ), intent(out) :: u(nn)
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
end subroutine cholesky

!-----------------------------------------------------------------------------
subroutine subchl ( a, b, n, u, nullty, ifault, ndim, det )

!*****************************************************************************80
! 
!! SUBCHL computes the Cholesky factorization of a (subset of a) PDS matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy, PR Freeman
!    FORTRAN90 version by  John Burkardt
!
!  Reference:
!
!    PR Freeman,
!    Remark AS R44:
!    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
!    and Inversion of a positive semi-definite symmetric matrix,
!    Applied Statistics,
!    Volume 31, Number 3, 1982, pages 336-339.
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((M*(M+1))/2), a positive definite matrix 
!    stored by rows in lower triangular form as a one dimensional array, 
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.  
!    In the simplest case, M, the order of A, is equal to N.
!
!    Input, integer ( kind = 4 ) B(N), indicates the order in which the
!    rows and columns of A are to be used.  In the simplest case, 
!    B = (1,2,3...,N).
!
!    Input, integer ( kind = 4 ) N, the order of the matrix, that is, 
!    the matrix formed by using B to select N rows and columns of A.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of A and U, which might 
!    be presumed to be (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndim

  real    ( kind = 8 ) a(ndim)
  integer ( kind = 4 ) b(n)
  real    ( kind = 8 ) det
  real    ( kind = 8 ), parameter :: eta = 1.0D-09
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nullty
  real    ( kind = 8 ) u(ndim)
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x

  ifault = 0
  nullty = 0
  det = 1.0D+00

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  ifault = 2
  j = 1
  k = 0

  do icol = 1, n

    ij = ( b(icol) * ( b(icol) - 1 ) ) / 2
    ii = ij + b(icol)
    x = eta * eta * a(ii)
    l = 0

    do irow = 1, icol

      kk = ( b(irow) * ( b(irow) + 1 ) ) / 2
      k = k + 1
      jj = ij + b(irow)
      w = a(jj)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        if ( abs ( x * a(kk) ) < w * w ) then
          ifault = 2
          return
        end if

        u(k) = 0.0D+00

      end if

    end do

    if ( abs ( eta * a(kk) ) <= abs ( w ) ) then

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    else

      u(k) = 0.0D+00
      nullty = nullty + 1

    end if

    j = j + icol
    det = det * u(k) * u(k)

  end do

  return
end subroutine subchl

real function multigauss(R,mu,n) result(G)
!!$ Generates a one dimensional array of length n containing multivariate Gaussian pseudo-random numbers
!!$ where R is the lower-triangular Cholesky factor of the covariance matrix of the desired distribution 
!!$ and mu is a one dimensional array of length n, containing the mean value for each random variable. 
!!$ R must be a one dimensional array, containing the lower-triangular matrix stored by row, starting from the 
!!$ uppermost, leftmost entry (first row, first column). 
   implicit none
   integer, intent(in) :: n
   real, intent(in) :: R(n*(n+1)/2)
   real, intent(in) :: mu(n)
   real :: Nu(n)
   dimension :: G(n)
   integer :: i, j
   
   if( n*(n+1)/2 .ne. size(R) ) then
      write(6,*) ' n*(n+1)/2 != size(R) in multiGauss '
      write(6,*) ' Cholesky factor matrix has size ', size(R) 
      write(6,*) ' but you are requesting a vector of size ',n
      stop
   else
      
!!$ generate an array of independent Gaussian variates
      do j = 1, n
         Nu(j) = gasdev()
      end do

!!$ start with the desired mean, and add the product of 
!!$ [the lower-triangular Cholesky factor of the covariance matrix] x [ the vector of iid Gaussian variates]
      do i = 1, n
         G(i) = mu(i)
         do j =  1,i
            G(i) = G(i) + R( i*(i-1)/2 + j ) * Nu(j)
         end do
      end do

   end if
   return
 end function multigauss
end module mGf90
