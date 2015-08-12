!**************************************************************
!*
!*  STATE
!*
!**************************************************************
subroutine state
implicit none
include 'runhydro.h'
!**************************************************************
!*
!    Using an ideal gas equation of state, assign the
!    pressure as a function of the density and internal
!    energy per unit mass
!
!      p  =  (gamma - 1) tau ** gamma
!
!    where tau = (rho * eps) ** 1 / gamma
!*
!**************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(4) :: np
real, dimension(4) :: kappa
real, dimension(num_species) :: gammainit
real :: rho_c1, rho_c2
common /bipoly/ np, kappa, gammainit, rho_c1, rho_c2

real, dimension(numr_dd,numz_dd,numphi,num_species) :: species
real, dimension(numr_dd,numz_dd,numphi) :: gammaeff
common /multispecies/ species, gammaeff
!*
!**************************************************************
!*
!*   Local variables
integer :: J, K, L
!*
!**************************************************************

call gamma_eff

do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         p(J,K,L) = (gammaeff(J,K,L) - 1.0) * tau(J,K,L)**gammaeff(J,K,L)
      enddo
   enddo
enddo

return
end subroutine state
