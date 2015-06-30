!************************************************************************
!* 
!*  DIAG_MI
!*
!************************************************************************
subroutine diag_mi(factor,star1com,star2com,mi1,mi2,mask)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!
!   diag_mi calculates the moment of inertia of each star wrt its
!   center of mass
!
!   the formula is:
!          
!   mi_i = SUM [ dm * ( (x - x_com-i)^2 - (y - y_com-i)^2 ) ]
!
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag, potr, vr, vz, vphi, x, y, xin, yin

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                      &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                           &
             in_neighbor, out_neighbor, root,                           &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                             &
           in_neighbor, out_neighbor, root,                             &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                      &
                        iam_on_bottom, iam_on_axis,                     &
                        iam_on_edge, down_neighbor,                     &
                        up_neighbor, in_neighbor,                       &
                        out_neighbor, root, column_num,                 &
                        row_num, pe_grid, iam_root,                     &
                        REAL_SIZE, INT_SIZE

!*
!************************************************************************
!*
!*   Local Variables

real, dimension(2) :: pass, passed

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: J, K, L

!*
!************************************************************************
!*
!*  Subroutine Arguments

real, dimension(3) :: star1com, star2com

real :: factor, mi1, mi2

integer, dimension(numr_dd,numz_dd,numphi) :: mask

!*
!************************************************************************
!  initialize the local variables
pass = 0.0
passed = 0.0
ierror = 0
 

!------------------------------------------------------------------------------
! sum up the z component of spin angular momentum

mi1 = 0.0
mi2 = 0.0
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( mask(J,K,L) < 0 ) then
            mi2 = mi2 + rho_diag(J,K,L) * rhf(J) *      &
                        ((x(J,L) - star2com(1))**2 +    &
                        (y(J,L) - star2com(2))**2)

         else if( mask(J,K,L) > 0 ) then
            mi1 = mi1 + rho_diag(J,K,L) * rhf(J) *      &
                        ((x(J,L) - star1com(1))**2 +    &
                        (y(J,L) - star1com(2))**2)     
          endif
      enddo
   enddo
enddo

!------------------------------------------------------------------------------
pass(1) = factor * mi1
pass(2) = factor * mi2

call mpi_reduce(pass,passed,2,REAL_SIZE,MPI_SUM,root,        &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then
   mi1 = passed(1)
   mi2 = passed(2)
endif
       
return
end subroutine diag_mi
