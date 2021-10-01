subroutine Calculate_PAR(I_0, nlev_, Hz, Chl, PAR, PARw)
  !top level is nlev_, bottom layer is 1, following ROMS convention
  implicit none
  real   , intent(in) :: I_0
  integer, intent(in) :: nlev_    ! Total number of vertical layers
  real   , intent(in) :: Hz(nlev_), Chl(nlev_)
  real   , intent(out):: PAR(nlev_),PARw(0:nlev_)
  integer             :: i
  real                :: par0, attn    ! Scratch variable
  real   , parameter  :: kw = 0.04, kc = 0.025

  par0         = I_0   !Light at the grid surface
  parw(nlev_)  = I_0
  do i = nlev_,1,-1

     attn   = exp(-0.5*Hz(i)*(kw + kc*Chl(i)))
     PAR(i) = par0*attn
     par0   = PAR(i)*attn
     PARw(i-1) = par0
  enddo
end subroutine Calculate_PAR

