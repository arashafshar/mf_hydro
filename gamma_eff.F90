!********************************************************************************
!*
!*  GAMMA_EFF
!*
!********************************************************************************
subroutine gamma_eff
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************************
!*
!   Calculates effective gamma in each cell 
!   
!                        SUM(i=1,num_species-1)  species(J,K,L,I) * gamma(I)
!   gammaeff(J,K,L)  =  ---------------------------------------------------
!                           SUM(i=1,num_species-1)   species(J,K,L,I)
!
!   gamma(num_species) is assumed to be vacuum and neglected in the calculations
!
    
!********************************************************************************
!*
!*  Global Variable

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(4) :: np
real, dimension(4) :: kappa
real, dimension(num_species) :: gammainit
real :: rho_c1, rho_c2
common /bipoly/ np, kappa, gammainit, rho_c1, rho_c2

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr_dd,numz_dd,numphi,num_species) :: species
real, dimension(numr_dd,numz_dd,numphi) :: gammaeff
common /multispecies/ species, gammaeff

!*
!***********************************************************************
!*
!*  Local Variable
real :: num, den
integer :: I, J, K, L

gammaeff = 0.0

do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         num = 0.0
         den = 0.0
! Sum over all species excluding vacuum
         do I = 1, num_species-1
            num = num + species(J,K,L,I) * gammainit(I)
            den = den + species(J,K,L,I)
         enddo
         if (den.le.0.0) then
            gammaeff(J,K,L) = gammainit(num_species)
         else
            gammaeff(J,K,L) = num / den
         endif
      enddo
   enddo
enddo


end subroutine gamma_eff
