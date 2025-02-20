SUBROUTINE Par2PHY
use state_variables, only : t, N_PAR, iPC, iPN, iChl, p_PHY, Varout, iTopt, iSize, ialphaChl
use state_variables, only : N_birth, N_mutate
use state_variables, only : oN_cell, oCDiv_avg, oCDiv_var, oTopt_avg, oTopt_var 
use state_variables, only : oLnalpha_var, oLnalpha_avg
use state_variables, only : oTalp_cov, oALnV_cov, oTLnV_cov
use grid,            only : Hz, nlev
use params,          only : NTrait, nu, sigma 
use mGf90,           only : srand_mtGaus
IMPLICIT NONE

!This subroutine calculate the total amount of concentrations of 
!phytoplankton carbon, nitrogen, and chl based on the cells present.
real      :: PHYC(nlev) = 0d0
real      ::    PHY(nlev) = 0d0
real      ::    CHL(nlev) = 0d0
real      ::   mTopt_(nlev) = 0d0
real      ::   vTopt_(nlev) = 0d0
real      ::   mCDiv_(nlev) = 0d0
real      ::   vCDiv_(nlev) = 0d0
real      ::   mlnalpha_(nlev) = 0d0
real      ::   vlnalpha_(nlev) = 0d0
real      ::   cov_TA(nlev) = 0d0
real      ::   cov_TL(nlev) = 0d0
real      ::   cov_AL(nlev) = 0d0
real      ::  nu_ = 0d0   !Basic Mutation rate
real      :: cff = 0.d0   !Random number [0,1]
real      :: oldtt(1) = 0.   !Scratch variable for storing the old trait
real      :: newtt(1) = 0.   !Scratch variable for storing the new trait
real      :: vartt(1,1) = 0.   !Variance of the mutating trait
integer :: k,i,m, ipar

!End of declaration

PHYC(:)= 0.d0
PHY(:) = 0.d0
CHL(:) = 0.d0
mCDiv_(:) = 0d0
mTopt_(:) = 0d0
mlnalpha_(:) = 0d0

!loop through all super-individuals to convert into Eulerian concentrations
DO i = 1, N_PAR
   ipar = p_PHY(i)%iz  !The current grid of super-individual i
   !Handle cell division and mutation
   ! If cellular carbon is above the division threshold, it divides
   IF (p_PHY(i)%C >= p_PHY(i)%Cdiv) THEN  !Divide
      N_birth(ipar)= N_birth(ipar) + 1
      p_PHY(i)%C   = p_PHY(i)%C/2d0
      p_PHY(i)%N   = p_PHY(i)%N/2d0
      p_PHY(i)%Chl = p_PHY(i)%Chl/2d0
      p_PHY(i)%num = p_PHY(i)%num*2d0

      !Mutation

      If (NTrait > 0) Then
         DO m = 1, NTrait
            nu_ = p_PHY(i)%num*nu(m)
            call random_number(cff)

            IF (cff < nu_) THEN !Mutation occurs
               N_mutate(ipar) = N_mutate(ipar) + 1

               select case(m)
               case(iTopt)
                   oldtt(1) = p_PHY(i)%Topt
               case(iSize)
                   oldtt(1) = log(p_PHY(i)%CDiv)
               case(ialphaChl)
                   oldtt(1) = p_PHY(i)%LnalphaChl
               case DEFAULT
                   stop "Trait index wrong!"
               end select

               vartt(1,1)= sigma(m)**2   !Construct the covariance matrix for the selected trait

               !A new Topt is randomly sampled from a Gaussian distribution with mean of previous Topt and SD of sigma
               newtt = srand_mtGaus(1, oldtt, vartt)
               select case(m)
               case(iTopt)
                   p_PHY(i)%Topt = newtt(1)
               case(iSize)
                   p_PHY(i)%CDiv = exp(newtt(1))
               case(ialphaChl)
                   p_PHY(i)%LnalphaChl = newtt(1)
               case DEFAULT
                   stop "Trait index wrong!"
               end select
            ENDIF
         ENDDO !End of looping through Traits
      ENDIF !End of if (NTrait > 0)
   ENDIF !End of division

   !Calculate Eulerian concentrations of phyto C, N, and Chl, mean trait  for each layer
   PHYC(ipar) = PHYC(ipar) + p_PHY(i)%num *p_PHY(i)%C 
   PHY(ipar)  = PHY(ipar)  + p_PHY(i)%num *p_PHY(i)%N
   CHL(ipar)  = CHL(ipar)  + p_PHY(i)%num *p_PHY(i)%Chl

   mCDiv_(ipar) = mCDiv_(ipar) + p_PHY(i)%num *p_PHY(i)%C * log(p_PHY(i)%Cdiv)
   mTopt_(ipar) = mTopt_(ipar) + p_PHY(i)%num *p_PHY(i)%C * p_PHY(i)%Topt
   mlnalpha_(ipar) = mlnalpha_(ipar) + p_PHY(i)%num * p_PHY(i)%C *p_PHY(i)%LnalphaChl
ENDDO !End of iterating over all super-individuals

DO k = 1, nlev
  t(iPC, k) = PHYC(k) *1d-9/Hz(k)!Convert Unit to mmol/m^3
  t(iPN,k) = PHY(k) *1d-9/Hz(k)   !Convert Unit to mmol/m^3
  t(iChl,k) = CHL(k)*1d-9/Hz(k) !Convert Unit to mmol/m^3

  if (Varout(oN_cell,k) > 0. .and. PHYC(k) > 0.) then
    mCDiv_(k)    = mCDiv_(k)/PHYC(k)
    mTopt_(k)    = mTopt_(k)/PHYC(k)
    mlnalpha_(k) = mlnalpha_(k)/PHYC(k)
  endif
ENDDO

Varout(iPC, :) = t(iPC, :)
Varout(iPN, :) = t(iPN, :)
Varout(iChl, :) = t(iChl,:)
Varout(oCDiv_avg, :) = mCDiv_(:)
Varout(oTopt_avg, :) = mTopt_(:)
Varout(oLnalpha_avg, :) = mLnalpha_(:)

!Compute trait covariances
vCDiv_(:) = 0d0
vTopt_(:) = 0d0
vlnalpha_(:) = 0d0
cov_AL(:) = 0d0
cov_TA(:) = 0d0
cov_TL(:) = 0d0

DO i = 1, N_PAR
   ipar = p_PHY(i)%iz  !The current grid of super-individual i
   vCDiv_(ipar) = vCDiv_(ipar) + p_PHY(i)%num*p_PHY(i)%C*(log(p_PHY(i)%Cdiv) - mCDiv_(ipar))**2
   vTopt_(ipar) = vTopt_(ipar) + p_PHY(i)%num*p_PHY(i)%C*(p_PHY(i)%Topt - mTopt_(ipar))**2
   vlnalpha_(ipar) = vlnalpha_(ipar) + p_PHY(i)%num*p_PHY(i)%C*(p_PHY(i)%LnalphaChl - mlnalpha_(ipar))**2
   cov_TL(ipar) = cov_TL(ipar)+p_PHY(i)%num*p_PHY(i)%C*(log(p_PHY(i)%Cdiv)-mCDiv_(ipar))*(p_PHY(i)%Topt - mTopt_(ipar))
   cov_AL(ipar) = cov_AL(ipar)+p_PHY(i)%num*p_PHY(i)%C*(log(p_PHY(i)%Cdiv)-mCDiv_(ipar))*(p_PHY(i)%LnalphaChl - mlnalpha_(ipar))
   cov_TA(ipar) = cov_TA(ipar)+p_PHY(i)%num*p_PHY(i)%C*(p_PHY(i)%LnalphaChl-mlnalpha_(ipar))*(p_PHY(i)%Topt - mTopt_(ipar))
ENDDO

do k = 1, nlev
  if (Varout(oN_cell,k) > 0. .and. PHYC(k) > 0d0) then
    vCDiv_(k)    = vCDiv_(k)   / PHYC(k)
    vTopt_(k)    = vTopt_(k)   / PHYC(k)
    vlnalpha_(k) = vlnalpha_(k)/ PHYC(k)
    cov_TL(k)    = cov_TL(k)   / PHYC(k)
    cov_TA(k)    = cov_TA(k)   / PHYC(k)
    cov_AL(k)    = cov_AL(k)   / PHYC(k)
  endif
enddo

Varout(oCDiv_var, :) = vCDiv_(:)
Varout(oTopt_var, :) = vTopt_(:)
Varout(oLnalpha_var, :) = vLnalpha_(:)
Varout(oALnV_cov, :) = cov_AL(:)
Varout(oTalp_cov, :) = cov_TA(:)
Varout(oTLnV_cov, :) = cov_TL(:)
return
END subroutine Par2PHY