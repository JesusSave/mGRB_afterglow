MODULE MOD_INITIAL

USE MOD_INPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: EN1       ! Circum number density
PUBLIC :: DC
PUBLIC :: eq_gamma_3
PUBLIC :: zbrent

CONTAINS



!***********************************************************************
! Function EN1: The number density ISM (k=0) or Wind (k=2)
! INPUT: r, radius
! OUTPUT: n
!***********************************************************************

FUNCTION EN1(R)
  real*8 :: en1, r
  if( kk==0 )then
     en1 = n0
  else
     en1 = (3.D35)*Aast/r**kk
  endif
  RETURN
END FUNCTION EN1



!***********************************************************************
! Function Dc: luminosity distance Dc
! INPUT: x, redshift
! OUTPUT: Dc
!***********************************************************************

function Dc( x, params ) bind(c)
use, intrinsic :: iso_c_binding
implicit none
  
real( kind = c_double )        :: Dc
real( kind = c_double ), value :: x
type( c_ptr ), value           :: params

Dc = 10.*c/H0/sqrt(OmegaM*(1+x)**3+OmegaL)
end function Dc 



!***********************************************************************
! Function EQ_GAMMA_3: initial gamma3
! INPUT: 
! - gam3, gam4
! - n4/n1
! OUTPUT: 
! - gamm3
! REMINDER: 
! - jump condition density relation func(gam4, gam4) = n3/n4
! - I want the MHD solution, so I multiply (1+sigma)
!***********************************************************************

FUNCTION eq_gamma_3(gamma_3)
implicit none

DOUBLE PRECISION                         :: eq_gamma_3
DOUBLE PRECISION                         :: gamma_3,gamma_34

gamma_34  = gamma_3*gam4 &
  -dsqrt((gamma_3*gamma_3-1.D0)*(gam4*gam4-1.D0))
eq_gamma_3 = (gamma_34-1.D0)*(4.D0*gamma_34+3.D0)*f*(1.D0+sigma) &
  -(gamma_3-1.D0)*(4.D0*gamma_3+3.D0)

END FUNCTION eq_gamma_3



!***********************************************************************
! Function ZBRENT: rooter algorithm
! INPUT: 
! - func : function
! - x1, x2 : possible range for the root
! - tol : toloration precision
! OUTPUT: 
! - root of func
!***********************************************************************

FUNCTION zbrent(func,x1,x2,tol)
EXTERNAL func
DOUBLE PRECISION                         :: zbrent,tol,x1,x2,func,EPS
DOUBLE PRECISION                         :: a,b,c,c1,c2,d,e,fa,fb,fc,p,q,r,s,tol1,xm
INTEGER*8                                :: ITMAX, iter
PARAMETER (ITMAX=1000,EPS=3.D-12)

  a=x1
  b=x2
  fa=func(a)
  fb=func(b)

  if((fa.gt.0.d0.and.fb.gt.0.d0).or.(fa.lt.0.d0.and.fb.lt.0.d0)) &
       pause 'root must be bracketed for zbrent'
  c=b
  fc=fb
  do iter=1,ITMAX
     if((fb.gt.0.d0.and.fc.gt.0.d0).or.(fb.lt.0.d0.and.fc.lt.0.d0))then
        c=a
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc).lt.abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.d0*EPS*abs(b)+0.5d0*tol
     xm=0.5d0*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.d0)then
        zbrent=b
        return
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa
        if(a.eq.c) then
           p=2.d0*xm*s
           q=1.d0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
           q=(q-1.d0)*(r-1.d0)*(s-1.d0)
        endif
        if(p.gt.0.d0) q=-q
        p=abs(p)
        if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        endif
     else
        d=xm
        e=d
     endif
     a=b
     fa=fb
     if(abs(d) .gt. tol1) then
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb=func(b)

  end do
  pause 'zbrent exceeding maximum iterations'
  zbrent=b
  return
END FUNCTION zbrent


END MODULE MOD_INITIAL
