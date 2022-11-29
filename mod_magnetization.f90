MODULE MOD_MAGNETIZATION

USE MOD_INPUT
use MOD_INITIAL

IMPLICIT NONE

PRIVATE

PUBLIC :: MAGNETIZATION_FACTOR

CONTAINS



!***********************************************************************
! Subroutine MAGNETIZATION_FACTOR
! This subroutine computes magnetization correction factors
! 
! INPUT: 
! - r: the radius used to calculate n4
! - gam3: Lorentz factor gam3
! - sigma: global variable
!
! OUTPUT: magnetization correction factors fa, fb, fc, hgam3
!
! REMINDER: 
! - the accurate adiabatic index is derived by iterations 
! - external library FGSL is used for bessel function kn1 and kn2
!***********************************************************************

SUBROUTINE MAGNETIZATION_FACTOR(R, GAM3)

use fgsl

IMPLICIT NONE

DOUBLE PRECISION                         :: r, gam3
DOUBLE PRECISION                         :: hgam0, a1, b1, c1, d1, n3, e3, p3, th
DOUBLE PRECISION                         :: bet4, bet3, n4, gam34
DOUBLE PRECISION                         :: rmin, rmid, rmax, exch
DOUBLE PRECISION                         :: M, N, O
DOUBLE PRECISION,dimension(4)            :: a

integer(fgsl_int)                        :: l, j, status
complex(fgsl_double)                     :: z(3)
type(fgsl_poly_complex_workspace)        :: w

!***********************************************************************

bet4   = sqrt(1.D0-1.D0/gam4**2)
bet3   = sqrt(1.D0-1.D0/gam3**2)
gam34  = (1.D0-bet3*bet4)*gam4*gam3
hgam = 5./3
hgam0 = 5./3

do j = 1,100

  ! ZK05 Expression
  !a1  = hgam*(2.D0-hgam)*(gam34-1.D0)+2.D0
  !b1  = -(gam34+1.D0)*((2.D0-hgam)*(hgam*gam34**2+1.D0)+hgam*(hgam-1.D0)*gam34)*sigma-&
  !      (gam34-1.D0)*(hgam*(2.D0-hgam)*(gam34**2-2.D0)+(2.D0*gam34+3.D0))
  !c1  = (gam34+1.D0)*(hgam*(1.D0-hgam/4.D0)*(gam34**2-1.D0)+1.D0)*sigma**2+&
  !      (gam34**2-1.D0)*(2.D0*gam34-(2.D0-hgam)*(hgam*gam34-1.D0))*sigma + &
  !      (gam34+1.D0)*(gam34-1.D0)**2*(hgam-1.D0)**2
  !d1  = -(gam34-1.D0)*(gam34+1.D0)**2*(2.D0-hgam)**2.D0*sigma**2/4.D0
  ! QC Expression
  a1  = hgam*(hgam-2.0)*gam34**2.0-2.0*(hgam-1.0)**2.0*gam34+(hgam-1.0)**2.0+1.0
  b1  = hgam*(2.0-hgam)*(1.0+sigma)*gam34**4.0&
    +((2.0+sigma)*hgam-2.0)*(hgam-1.0)*gam34**3.0&
    +((1.0+sigma)*hgam**2.0-(3.0*sigma+2.0)*hgam+2.0*sigma-1.0)*gam34**2.0&
    +(hgam-1.0)*(4.0-hgam*(sigma+4.0))*gam34+2.0*hgam**2.0+(sigma-4.0)*hgam-2.0*sigma+3.0
  c1 = (0.25*(sigma**2.0-4.0*sigma-4.0)*hgam**2.0-(sigma**2.0-2.0*sigma-2.0)*hgam-2.0*sigma-1.0)*gam34**4.0&
    +(hgam-1.0)*(hgam*(sigma+2.0)-2.0)*gam34**3.0&
    +(-0.5*sigma*(sigma-2.0)*hgam**2.0+sigma*(2.0*sigma-3.0)*hgam-sigma*(sigma-4.0))*gam34**2.0&
    +(hgam-1)*(2.0-sigma*hgam-2.0*hgam)*gam34&
    +0.25*(sigma**2.0+4.0)*hgam**2.0-(sigma**2.0-sigma+2.0)*hgam+(1.0-sigma)**2.0
  d1  = (gam34**2-1.0)**2.0*(2.0-hgam)**2.0*sigma**2.0/4.0

  a = (/d1/c1,c1/c1,b1/c1,a1/c1/)
  w = fgsl_poly_complex_workspace_alloc (4_fgsl_size_t)
  status = fgsl_poly_complex_solve (a, 4_fgsl_size_t, w, z)
  call fgsl_poly_complex_workspace_free (w)
  rmin = real(z(1))
  rmid = real(z(2))
  rmax = real(z(3))

  if (rmin >= rmid) then
    exch = rmid
    rmid = rmin
    rmin = exch
  elseif (rmax >= rmid) then
  elseif (rmax >= rmin) then
      rmid = rmax
  endif
  u3s = sqrt(rmid)

  M  = sqrt(u3s*u3s+1.0)*sqrt(gam34*gam34-1.0)
  fa = 1.0-(gam34+1.0)/2.0/(u3s*u3s*gam34+u3s*M)*sigma

  N  = (hgam*gam34+1.0)/(hgam-1.0)
  fb = (gam34+M/u3s)/N

  O  = 1.0/(hgam-1.0)/2.0/(gam34-1.0)
  fc = 1.0+O*fb*N/fa*sigma

  !u4s    = u3s*gam34+dsqrt(xxx+1.D0)*dsqrt(gam34**2-1.D0)
  n4   = mej/(4.D0*pi*r*r*gam4*max(Delta0,r/gam4**2)*mp*(1.D0+sigma))

  n3 = n4*(hgam*gam34+1.0)/(hgam-1.0)*fb
  e3 = n3*mp*c*c*(gam34-1.0)*fa
  p3 = (hgam-1.0)*e3*fc
  th = p3/(n3*mp*c*c)

  ! use fgsl_sf_bessel_kcn(1,1.d0) or fgsl_sf_bessel_knu(1.d0,1.d0)
  ! to solve Modified Bessel function of second kind

  if (th>2e-3)then
    hgam = 1.0+th/(3.0*th+fgsl_sf_bessel_kcn(1,1./th)/fgsl_sf_bessel_kcn(2,1./th)-1.0)
    if ((hgam >5./3).or.(hgam<4./3)) then
      hgam = 4./3
    endif
  else
    hgam = 5./3
  endif

  if (abs(hgam-hgam0)<1e-12)then
    MAX_ITERATIONS = j
    exit
  endif
    hgam0 = hgam

end do

END SUBROUTINE MAGNETIZATION_FACTOR



END MODULE MOD_MAGNETIZATION

