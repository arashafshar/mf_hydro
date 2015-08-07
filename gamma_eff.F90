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
!   Calculates effective gamma in each cell, which is a mass weighed average
!   
!                            species(J,K,L,I) 
!   f(J,K,L,I)  =  --------------------------------------  
!                  SUM(i=1,num_species) species(J,K,L,I) 
!
!   m(J,K,L,I)  =  rho (I,J,K) *  vol(I,J,K) * f(I,J,K,L)
! 
!                        SUM(i=1,num_species)  m(J,K,L,I) * gamma(I)
!   gamma_eff(J,K,L)  =  -------------------------------------------
!                           SUM(i=1,num_species)   m(J,K,L,I)
!
    
!********************************************************************************
!*
!*  Global Variable

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real :: kappac1, kappac2, rho_c1, rho_c2, np1, np2, gamma1, gamma2  !bipoly 
real, dimension(num_species) :: gammainit
common /bipoly/ kappac1, kappac2, rho_c1, rho_c2, np1, np2, gamma1, &
       gamma2, gammainit

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
real, dimension(numr_dd,numz_dd,numphi,num_species) :: massys 
real, dimension(numr_dd,numz_dd,numphi) :: sumspecies

real :: fractional, num, den
integer :: I, J, K, L

!array for the mass of each passive scalar in species


massys = 0.0
sumspecies = 0.0
gammaeff = 0.0

do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         do I = 1, num_species
            sumspecies(J,K,L)=sumspecies(J,K,L)+species(J,K,L,I)            
         enddo
      enddo
   enddo
enddo


do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         do I = 1, num_species
            fractional = species(J,K,L,I)/sumspecies(J,K,L)
            massys(J,K,L,I) = rho(J,K,L) * rhf(J) * dr * dz * dphi * fractional
         enddo
      enddo
   enddo
enddo


do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         num = 0.0
         den = 0.0
         do I = 1, num_species
            num = num + massys(J,K,L,I) * gammainit(I)
            den = den + massys(J,K,L,I)
         enddo
         gammaeff(J,K,L) = num / den
      enddo
   enddo
enddo


end subroutine gamma_eff
