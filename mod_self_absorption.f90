MODULE MOD_SELF_ABSORPTION
USE FGSL
USE MOD_INPUT
USE MOD_INITIAL

IMPLICIT NONE

PRIVATE

PUBLIC :: nu_absorb
PUBLIC :: nu_syn_absorb
PUBLIC :: eq_nu_a_3
PUBLIC :: gammln

CONTAINS



!***********************************************************************
! FUNCTION nu_absorb
! This function is used for
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function nu_absorb(Ra, n, B, gm, gc, Ntot)

IMPLICIT NONE

DOUBLE PRECISION                         :: nu_absorb
DOUBLE PRECISION                         :: Ra, B, gm, gc, Ntot, n, Ma
DOUBLE PRECISION                         :: length, nuc, num

!***********************************************************************

  nuc = 0.45D0*gc*gc*e*B/(2*pi*me*c)
  num = 0.45D0*gm*gm*e*B/(2*pi*me*c)
  !print*,'Ra, n, B, gm, gc, Ntot',Ra, n, B, gm, gc, Ntot
  !print*,'nuc',nuc
  !print*,'num',num

  length = Ntot/(4.D0*pi*Ra*Ra*n)
  if( num<nuc ) then
     nu_absorb = nu_syn_absorb(B,gm,gc,length,n,num,nuc,p)
  elseif( num>nuc ) then
     nu_absorb = nu_syn_absorb(B,gc,gm,length,n,nuc,num,2.D0)
  elseif( num==nuc ) then
     print*,'nu_m = nu_c!!!'
     print*,'It''s too strange, I can''t deal with it.:('
     pause
  else
  endif
  return
end function nu_absorb



!***********************************************************************
! FUNCTION nu_absorb
! This function is used for
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

FUNCTION nu_syn_absorb(B,gamma_1,gamma_2,length,n,nu_1,nu_2,p)

IMPLICIT NONE

DOUBLE PRECISION                         :: nu_syn_absorb
DOUBLE PRECISION                         :: A_1,A_2,A_2of2,nu_a_1,nu_a_2,nu_a_3,nu_a_2of2
DOUBLE PRECISION                         :: B,gamma_1,gamma_2,length,n,nu_1,nu_2,p
DOUBLE PRECISION                         :: N_gamma, tau_0

!***********************************************************************

  if( nu_1>=nu_2 )then
     print*,'Error in calculating the nu_a.',&
          '(nu_1 should be less than nu_2)'
     pause
  endif
  tau_0     = 0.35D0
  N_gamma   = n*(p-1.D0)*gamma_1**(p-1.D0)

  C_N_gamma = N_gamma
  C_B       = B
  C_gamma_2 = gamma_2
  C_length  = length
  C_p       = p  
  C_nu_2    = nu_2
  C_tau_0   = tau_0

  A_1       = (pi*pi*2.D0**(8.D0/3.D0)/9.D0/dexp(gammln(1.D0/3.D0))   &
       *(pi/2.D0))*(p+2.D0)/(p+2.D0/3.D0)
  A_2       = (8.D0*2.D0**(p/2.D0)*pi/9.D0*(pi/2.D0))                 &
       *dsqrt(3.D0)/4.D0 &
       *dexp(gammln((3.D0*p+2.D0)/12.D0))                             &
       *dexp(gammln((3.D0*p+22.D0)/12.D0))
  A_2of2    = (8.D0*2.D0**((p+1.D0)/2.D0)*pi/9.D0*(pi/2.D0))          &
       *dsqrt(3.D0)/4.D0 &
       *dexp(gammln((3.D0*(p+1.D0)+2.D0)/12.D0))&
       *dexp(gammln((3.D0*(p+1.D0)+22.D0)/12.D0))
  nu_a_1    = (A_1*e*N_gamma*nu_1**(5.D0/3.D0)*length            &
       /(B*gamma_1**(p+4.D0)*tau_0))**(0.6D0)
  nu_a_2    = (A_2*e*N_gamma*nu_1**(p/2.D0+2.D0)*length          &
       /(B*gamma_1**(p+4.D0)*tau_0))**(2.D0/(p+4.D0))
  !nu_a_3    = zbrent(eq_nu_a_3, 10.D0, nu_2*10.D0, 1.D-8) ! sometimes no root
  nu_a_3    = zbrent(eq_nu_a_3, 1.D-5, nu_2*1.D5, 1.D-8) ! QChen
  nu_a_2of2 = (A_2of2*e*N_gamma*nu_2**(p/2.D0+2.5D0)*length  &
       /(B*gamma_2**(p+4.D0)*tau_0))**(2.D0/(p+5.D0))
  if(.false.)then
     print*,' nu_1:',nu_1,'  nu_2:',nu_2
     print*,' nu_a_1:',nu_a_1,'  nu_a_2:',nu_a_2
     print*,' nu_a_3:',nu_a_3,'  nu_a_2of2:',nu_a_2of2
     print*
  endif
  if(nu_a_1 .gt. nu_1) nu_a_1 = 0.D0
  if((nu_a_2 .lt. nu_1) .or. (nu_a_2 .gt. nu_2)) nu_a_2 = 0.D0
  if(nu_a_3 .lt. nu_2) nu_a_3 = 0.D0
  if(nu_a_2of2 .lt. nu_2) nu_a_2of2 = 0.D0
  nu_syn_absorb = dmax1(nu_a_1,nu_a_2,nu_a_3,nu_a_2of2)
  if(nu_syn_absorb .eq. nu_a_3) print*,'nu_syn_absorb = nu_a_3.', &
       'In this case, I can''t deal with the emission.'
  if(nu_syn_absorb .eq. 0.D0) pause 'Error in nu_syn_absorb.'
END FUNCTION nu_syn_absorb



!***********************************************************************
! FUNCTION nu_absorb
! This function is used for
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

FUNCTION eq_nu_a_3(nu_a_3)

IMPLICIT NONE

DOUBLE PRECISION                         :: nu_a_3,eq_nu_a_3

!***********************************************************************

eq_nu_a_3 = (dsqrt(6.D0)/36.D0)*pi**2.5D0*(C_p+2.D0)*C_N_gamma*(e/C_B)&
       *C_gamma_2**(-C_p-4.D0)*(nu_a_3/C_nu_2)**(-2.5D0)                &
       *dexp(-nu_a_3/C_nu_2)*C_length-C_tau_0
end FUNCTION eq_nu_a_3



!***********************************************************************
! FUNCTION gammln
! This function is used for
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

FUNCTION gammln(xx)

IMPLICIT NONE

DOUBLE PRECISION                         :: gammln,xx
DOUBLE PRECISION                         :: ser,stp,tmp,x,y,cof(6)
INTEGER*8                                :: j

!***********************************************************************
SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,      & 
       24.01409824083091d0,  -1.231739572450155d0,            &
       .1208650973866179d-2, -.5395239384953d-5,              &
       2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammln=tmp+log(stp*ser/x)
  return
END function gammln

END MODULE MOD_SELF_ABSORPTION

