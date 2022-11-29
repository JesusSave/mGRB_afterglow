MODULE MOD_DYNAMICS

USE FGSL
USE MOD_INPUT
USE MOD_INITIAL
USE MOD_MAGNETIZATION

IMPLICIT NONE

PRIVATE

PUBLIC :: FS_BEFORE_CROSS
PUBLIC :: FS_AFTER_CROSS
PUBLIC :: RS_AFTER_CROSS

PUBLIC :: RKDUMB ! DRIVER FOR RK4 
PUBLIC :: RK4 ! RUNGE KUTTA 4 ORDER ODE SOLVER

CONTAINS



!***********************************************************************
! Subroutine FS_BEFORE_CROSS
! This subroutine computes FS RS dynamics before crossings
! 
! INPUT: 
! - x: the logrithm of radius
! - y: the differential variables (obs_t2, gamma2, t3, gamma3, m2, m3)
!
! OUTPUT: differential of dy/dx
!
! REMINDER: 
!***********************************************************************

subroutine FS_BEFORE_CROSS(x,y,dydx)

IMPLICIT NONE

DOUBLE PRECISION                         :: x, y(*), dydx(*)
DOUBLE PRECISION                         :: gam2, gam3, bet4, bet2, bet3, betrs, gam34, bet34, bet3s
DOUBLE PRECISION                         :: tempQ, tempP, trans, temp
DOUBLE PRECISION                         :: W, W2, P2, P3

!***********************************************************************

gam2   = y(2)
gam3  = gam2
bet4   = sqrt(1.D0-1.D0/gam4**2)
bet2   = sqrt(1.D0-1.D0/gam2**2)
bet3   = bet2
gam34  = (1.D0-bet3*bet4)*gam4*gam3
bet34  = sqrt(1.0-1.0/gam34**2.0)

!bet3s  = (hgam-1.0)*(gam34-1.0)/bet34/gam34 ! QC: valid for hydrodynamics
!betrs  = (gam2*bet2*(4.D0*gam34+3.D0)*fb-gam4*bet4)/(gam2*(4.D0*gam34+3.D0)*fb-gam4)
!bet3s  = (bet3-betrs)/(1.0-betrs*bet3) ! Calculate from betrs is equivalent from u3s

bet3s = sqrt(u3s*u3s/(u3s*u3s+1.0))
betrs  = (bet3-bet3s)/(1.0-bet3*bet3s)

!temp   = mej/(4.D0*pi*10**x*10**x*gam4*max(Delta0,10**x/gam4**2)*mp)/Fparameter !n4
temp   = mej/(4.D0*pi*10**x*10**x*gam4*max(Delta0,10**x/gam4**2)*mp) !n4

P2     = y(2)**2-1.0
P3     = y(2)-1.0+y(2)*(gam34-1.D0)*fa*(1.D0+(fc-1.D0)*(hgam-1.0))-(gam4-1.0)*(1.D0 + sigma)

tempQ  = P2*en1(10**x)*mp + P3*gam4*temp*mp*(bet4-betrs)/betrs

W      = (1.D0-eps3+(fc-1.D0)*(hgam-1.0))*y(2)*(gam34-1.D0)*y(6)
W2     = y(2)*(hgam-1.0)*(gam34-1.D0)*y(6)*fa

tempP  = y(5)+y(6)&
         +(1.D0-eps2)*(2.D0*y(2)-1.D0)*y(5)&
         +(fc-1.D0)*(hgam-1.0)*(gam34-1.D0)*y(6)*fa&
         +(1-eps3)*y(6)*fa*(gam34-1.D0+y(2)*(gam4-gam4*bet4/bet2))&
         +y(2)*(fc-1.D0)/3.D0*y(6)*fa*(gam4-gam4*bet4/bet2) 

trans  = (log(10.D0)*10**x)                          !dR/dx

dydx(1) = trans*((1.D0-bet2)*(1.D0+redsh)/bet2/c)    !dt/dx
!dydx(2) = trans*min(0.0,-4.D0*pi*10**x*10**x*tempQ/tempP/(1.D0+W*faa/tempP+W2*fca/tempP))*(1.0-cos(thetaj))/2.
dydx(2) = trans*(-4.D0*pi*10**x*10**x*tempQ/tempP/(1.D0+W*faa/tempP+W2*fca/tempP))*(1.0-cos(thetaj))/2.
dydx(3) = dydx(1)
dydx(4) = dydx(2)
dydx(5) = trans*4.D0*pi*10**x*10**x*en1(10**x)*mp*(1.0-cos(thetaj))/2.
dydx(6) = trans*4.D0*pi*10**x*10**x*(bet4-betrs)/betrs*gam4*temp*mp*(1.0-cos(thetaj))/2.

return
end subroutine FS_BEFORE_CROSS



!***********************************************************************
! Subroutine FS_AFTER_CROSS
! This subroutine computes FS RS dynamics before crossings
! 
! INPUT: 
! - x: the logrithm of radius
! - y: the differential variables (obs_t2, gamma2, t3, gamma3, m2, m3)
!
! OUTPUT: differential of dy/dx
!
! REMINDER: 
!***********************************************************************

SUBROUTINE FS_AFTER_CROSS(x,y,dydx)

IMPLICIT NONE

DOUBLE PRECISION                         :: x, y(*), dydx(*)
DOUBLE PRECISION                         :: bet2, trans

!***********************************************************************

bet2    = sqrt(1.0D0-1.0D0/y(2)**2)
trans   = (log(10.D0)*10**x)
dydx(1) = trans*((1.D0-bet2)*(1.D0+redsh)/bet2/c)
dydx(2) = trans*(-4.D0*pi*10**x*10**x*en1(10**x)*mp*(y(2)**2-1.D0)/  &
     ((eps2+2.0D0*(1.0D0-eps2)*y(2))*y(5)))*(1.0-cos(thetaj))/2.
dydx(3) = 0.D0
dydx(4) = 0.D0
dydx(5) = trans*4.D0*pi*10**x*10**x*en1(10**x)*mp*(1.0-cos(thetaj))/2.
dydx(6) = 0.D0

return
end subroutine FS_AFTER_CROSS



!***********************************************************************
! Subroutine FS_BEFORE_CROSS
! This subroutine computes FS RS dynamics before crossings
! 
! INPUT: 
! - x: the logrithm of radius
! - y: the differential variables (obs_t2, gamma2, t3, gamma3, m2, m3)
!
! OUTPUT: differential of dy/dx
!
! REMINDER: 
! - after cross, region 4 vanish, gam34, u3s, fa, fb, fc should vanish
! - however we use dummy values by suspending the cross time values
! - indeed after cross FS is more important than RS
! - we use Blandford-McKee self-similar solutions, only valid for relativistic stage!
!***********************************************************************

SUBROUTINE RS_AFTER_CROSS(nstepD,x1,x2)

IMPLICIT NONE

DOUBLE PRECISION                         :: x1, x2
DOUBLE PRECISION                         :: f, g, ratio
INTEGER*8                                :: i, nstepD

DOUBLE PRECISION                         :: pp, nn, ee, th
DOUBLE PRECISION                         :: hgam1, hgam0
integer(fgsl_int)                        :: l, j
!***********************************************************************

f = n4(ncross)/en1(y(4,ncross))

do i=ncross+1, nstepD
  xx(i) = xx(ncross)+(i-ncross)*(x2-x1)/nstepD
  ratio = 10**(xx(i)-xx(ncross))                 !ratio = r/rx

  if( kk==0.D0 )then
    if(f >= gam4*gam4) then !Thin
      g = 2.D0
      y(3,i)= y(3,ncross)*ratio**(2.D0*g+1.D0)       !t3:(r/rx)**(2g+1)
      y(4,i)= y(4,ncross)*ratio**(-g)                !gam3:(r/rx)**(-g)
      e3(i) = e3(ncross)*ratio**(-8.D0*(3.D0+g)/7.D0)
      n3(i) = n3(ncross)*ratio**(-6.D0*(3.D0+g)/7.D0)
    else !Thick shell
      g = 3.5D0
      y(3,i)= y(3,ncross)*ratio**(2.D0*g+1.D0)       !t3:(r/rx)**(2g+1)
      y(4,i)= y(4,ncross)*ratio**(-g)                !gam3:(r/rx)**(-g)
      n3(i) = n3(ncross)*ratio**((2.0*g+1.0)*(-13./16))
      e3(i) = e3(ncross)*ratio**((2.0*g+1.0)*(-13./12))
    endif
  elseif( kk==2.D0 ) then
    if(f >= gam4*gam4)then
      g = 1.D0
      y(3,i)= y(3,ncross)*ratio**(2.D0*g+1.D0)       !t3:(r/rx)**(2g+1)
      y(4,i)= y(4,ncross)*ratio**(-g)                !gam3:(r/rx)**(-g)
      e3(i) = e3(ncross)*ratio**(-8.D0*(3.D0+g)/7.D0)
      n3(i) = n3(ncross)*ratio**(-6.D0*(3.D0+g)/7.D0)
    else
      g = 1.5D0
      y(3,i) = y(3,ncross)*ratio**(2.D0*g+1.D0)       !t3:(r/rx)**(2g+1)
      y(4,i) = y(4,ncross)*ratio**(-g)                !gam3:(r/rx)**(-g)
      n3(i) = n3(ncross)*ratio**((2.0*g+1.0)*(-9./8))
      e3(i) = e3(ncross)*ratio**((2.0*g+1.0)*(-3./2))
    endif
  else
    print*,'I can''t deal with this case. kk=',kk
    stop
  endif

  u3s_dump(i) = u3s_dump(ncross) 
  fa_dump(i) = fa_dump(ncross)
  fb_dump(i) = fb_dump(ncross)
  fc_dump(i) = fc_dump(ncross)

  !---iterations for hgam3 and p3, since should fc vanish, the results then not accurate
  nn = n3(i)
  ee = e3(i)
  hgam = 5./3
  hgam0 = 5./3
  do j = 1,100
    pp = (hgam-1.0)*ee*fc_dump(i)
    th = pp/(nn*mp*c*c)
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
 p3(i) = pp
 hgam3_dump(i) = hgam
end do

return
end subroutine RS_AFTER_CROSS



!***********************************************************************
! Subroutine RKDUMB
! This subroutine is the DRIVER for rk4 ODEs solver algorithm
!***********************************************************************

SUBROUTINE RKDUMB(vstartD,nvarD,x1,x2,nstepD,derivs,derivs0)

IMPLICIT NONE

INTEGER*8                                :: nvarD, nstepD, i, j
DOUBLE PRECISION                         :: vstartD(nvarD), v(nvarD), dv(nvarD)
DOUBLE PRECISION                         :: x1, x2, h, x
DOUBLE PRECISION                         :: gam34, bet4, bet2, betrs

DOUBLE PRECISION                         :: gam, p, n, e, th, r
DOUBLE PRECISION                         :: hgam1, hgam0

!***********************************************************************

do i=1,nvarD
  v(i)=vstartD(i)
  y(i,1)=v(i)
end do

xx(1)=x1
x=x1
h=(x2-x1)/nstepD

do j=1,nstepD-1

  call MAGNETIZATION_FACTOR(10**xx(j),y(2,j))
  u3s_dump(j) = u3s
  fa_dump(j) = fa
  fb_dump(j) = fb
  fc_dump(j) = fc
  hgam3_dump(j) = hgam
  if (j==1) then
    faa = fa/gam4 ! guess
    fca = fc/gam4 ! guess
  else
    faa = (fa-fa_dump(j-1))/(y(2,j)-y(2,j-1))
    fca = (fc-fc_dump(j-1))/(y(2,j)-y(2,j-1))
  endif

  call derivs(x,v,dv)
  call rk4(v,dv,nvarD,x,h,v,derivs)

  x=x+h
  xx(j+1)=x

  gam34   = gam4*y(2,j)-dsqrt((gam4*gam4-1.D0)*(y(2,j)*y(2,j)-1.D0))
  n4(j) = mej/(4.D0*pi*10**xx(j)*10**xx(j)*gam4*max(Delta0,10**xx(j)/gam4**2)*mp) !n4
  !n4(j) = mej/(4.D0*pi*10**xx(j)*10**xx(j)*gam4*max(Delta0,10**xx(j)/gam4**2)*mp)/fa/fb/fc 
  n3(j) = (hgam*gam34+1.D0)/(hgam-1.D0)*n4(j)*fb
  e3(j) = n3(j)*mp*c**2*(gam34-1.D0)*fa
  p3(j) = (hgam-1.0)*e3(j)*fc

  if( v(6).gt.mej )then   ! Deceleration time when swept mass .eq. initial mass
    ncross=j            ! This is CROSS TIME point
    exit
  end if

  do i=1,nvarD
   y(i,j+1)=v(i)
  end do

end do


do j=ncross,nstepD-1

  call derivs0(x,v,dv)
  call rk4(v,dv,nvarD,x,h,v,derivs0)

  n4(j) = mej/(4.D0*pi*10**xx(j)*10**xx(j)*gam4*max(Delta0,10**xx(j)/gam4**2)*mp)
  x=x+h
  xx(j+1)=x

  do i=1,nvarD
    y(i,j+1)=v(i)
  end do

end do


do i=1,nstep

  gam = y(2,i)
  r = 10**xx(i)
  hgam = 5./3
  hgam0 = 5./3
  do j = 1,100
    n = (hgam*gam+1.0)/(hgam-1.0)*en1(r)
    e = (gam-1.0)*n*mp*c*c
    p = (hgam-1.0)*e
    th = p/(n*mp*c*c)
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

  hgam2_dump(i) = hgam
  p2(i) = p
  e2(i) = e
  n2(i) = n

end do

return
END SUBROUTINE RKDUMB



!***********************************************************************
! Subroutine RK4
! This subroutine is rk4 ODEs solver algorithm
!***********************************************************************

SUBROUTINE RK4(y,dydx,n,x,h,yout,derivs)

IMPLICIT NONE

INTEGER*8                                :: i, n
DOUBLE PRECISION                         :: y(n),dydx(n),yout(n),yt(nvar),dyt(nvar),&
       dym(nvar)
DOUBLE PRECISION                         ::  hh, h6, xh, x, h

!***********************************************************************

hh=h*0.5
h6=h/6.
xh=x+hh

do i=1,n
  yt(i)=y(i)+hh*dydx(i)
end do
call derivs(xh,yt,dyt)
do i=1,n
  yt(i)=y(i)+hh*dyt(i)
end do
call derivs(xh,yt,dym)
do i=1,n
  yt(i)=y(i)+h*dym(i)
  dym(i)=dyt(i)+dym(i)
end do
call derivs(x+h,yt,dyt)
do i=1,n
  yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
end do

END SUBROUTINE RK4



END MODULE MOD_DYNAMICS
