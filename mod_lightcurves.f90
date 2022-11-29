MODULE MOD_LIGHTCURVES
USE FGSL
USE MOD_INPUT
USE MOD_INITIAL
USE MOD_SELF_ABSORPTION

IMPLICIT NONE

PRIVATE

PUBLIC :: EqTfluxFs
PUBLIC :: EqTfluxRs
PUBLIC :: anglength
PUBLIC :: interpoly
PUBLIC :: interpolyne
PUBLIC :: obsfluxfs
PUBLIC :: obsfluxrs
PUBLIC :: frequencyrs
PUBLIC :: frequencyfs
PUBLIC :: eq_slow_gamma_c

CONTAINS



!***********************************************************************
! FUNCTION EqTfluxFs
! This function is used for FS EMISSION
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function EqTfluxFs(te)

IMPLICIT NONE

DOUBLE PRECISION                         :: EqTfluxFs, te, tedyn
DOUBLE PRECISION                         :: Ga, MA, THA, Ra
DOUBLE PRECISION                         :: dtheta,a,cosa,teth,dr,beta,dte,ang,ang1,ang2,Awidth,Alength
DOUBLE PRECISION                         :: partflux, totalflux, dw, totalw, Ntot, Npart
integer*8                                :: thstep,j,k
integer*8                                :: one, two, five

!***********************************************************************

  one = 1
  two = 2
  five= 5
  EqTfluxFs=0.D0
  if( deltaA < thetaj ) then
     ang1=0.0
  else 
     ang1 = deltaA - thetaj 
  end if

  ang2=thetaj+deltaA
  thstep=500
  dtheta=(ang2-ang1)/thstep
  do j=1, thstep
     a=ang1+dtheta*j
     cosa=cos(a)
     teth=y(1,1)  
     do k=1, nstep-1
        dr=(10**xx(k+1)-10**xx(k))
        beta=sqrt(1.0-1.0/y(2,k)/y(2,k))
        dte=(1.0-beta*cosa)*(1.D0+redsh)*dr/beta/c
        if( (teth+dte)>te ) then
           exit
        end if
        teth=teth+dte
     end do
     dte=te-teth
     dr=c*beta*dte/(1.0-beta*cosa)/(1.D0+redsh)
     Ra=10**xx(k)+dr             !/* this is the radius at angle a  */  

     Ga=interpoly(Ra, xx, y, two, k)
     Ma=interpoly(Ra, xx, y, five, k)
     THa=thetaj                 !/* this is the jet opening angle */
     tedyn=interpoly(Ra, xx, y, one, k)

     if( a < (deltaA-THa) ) then
        cycle
     end if
    
     if( (a-dtheta) < (deltaA-THa)) then
        Awidth=a+THa-deltaA
        ang=a
        Alength=anglength(deltaA,ang,THa)
        goto 330
     end if
    
     if( (a>(THa+deltaA)) .and. ((a-dtheta) >= (THa+deltaA)) )then
        return 
     endif
     
     if( (a>(THa+deltaA)) .and. ((a-dtheta)<(THa+deltaA)))then
        Awidth=THa+deltaA-a+dtheta
        ang=THa+deltaA
        Alength=anglength(deltaA,ang,THa)
        goto 330
     end if
     Awidth=dtheta
     ang=a
     Alength=anglength(deltaA,ang,THa)
     
330  dw=Alength*sin(ang)*Awidth
     totalw=2.0*PI*(1.0-cos(THa))
     Npart=Ma/Mp*dw/totalw
     Ntot=Ma/Mp
     partflux=obsfluxFs(tedyn, Ga, Ra, Npart, a, Ntot)
     EqTfluxFs=EqTfluxFs+partflux
  end do
  return
end function EqTfluxFs



!***********************************************************************
! FUNCTION EqTfluxRs
! This function is used for RS EMISSION
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function EqTfluxRs(te)

IMPLICIT NONE

DOUBLE PRECISION                         :: EqTfluxRs, te,tedyn
DOUBLE PRECISION                         :: Ga, MA, THA, Ra, E3a, N3a, Fca
DOUBLE PRECISION                         :: dtheta,a,cosa,teth,dr,beta,dte,ang,ang1,ang2,Awidth,Alength
DOUBLE PRECISION                         :: partflux, totalflux, dw, totalw, Ntot, Npart
integer*8                                :: thstep,j,k
integer*8                                :: six, three, four

!***********************************************************************

  six=6
  three=3
  four=4
  EqTfluxRs=0.D0
  if( deltaA < thetaj ) then
     ang1=0.0
  else 
     ang1 = deltaA - thetaj 
  end if

  ang2=thetaj+deltaA
  thstep=500
  dtheta=(ang2-ang1)/thstep
  do j=1, thstep
     a=ang1+dtheta*j
     cosa=cos(a)
     teth=y(3,1) !-gam4*wid*(1-cosa)*(1.D0+redsh)/c/2.D0/gam4/gam4  
     do k=1, nstep-1
        dr=10**xx(k+1)-10**xx(k)
        beta=sqrt(1.0-1.0/y(4,k)/y(4,k))
        dte=(1.0-beta*cosa)*(1.D0+redsh)*dr/beta/c
        if( (teth+dte)>te ) then
           exit
        end if
        teth=teth+dte
     end do

     dte=te-teth
     dr=c*beta*dte/(1.0-beta*cosa)/(1.D0+redsh)
     Ra=10**xx(k)+dr             !/* this is the radius at angle a  */  
     Ma = interpoly(Ra, xx, y, six, k) 
     Ga = interpoly(Ra, xx, y, four, k)
     N3a= interpolyne(Ra, xx, n3, k)
     E3a= interpolyne(Ra, xx, e3, k)
     Fca= interpolyne(Ra, xx, fc_dump, k) !QC
     THa=thetaj                      !/* this is the jet opening angle */
     tedyn=interpoly(Ra, xx, y, three, k)

     if( a < (deltaA-THa) ) then
        cycle
     end if
    
     if( (a-dtheta) < (deltaA-THa)) then
        Awidth=a+THa-deltaA
        ang=a
        Alength=anglength(deltaA,ang,THa)
        goto 330
     end if
    
     if( (a>(THa+deltaA)) .and. ((a-dtheta) >= (THa+deltaA)) )then
        return 
     endif
     
     if( (a>(THa+deltaA)) .and. ((a-dtheta)<(THa+deltaA)))then
        Awidth=THa+deltaA-a+dtheta
        ang=THa+deltaA
        Alength=anglength(deltaA,ang,THa)
        goto 330
     end if
     Awidth=dtheta
     ang=a
     Alength=anglength(deltaA,ang,THa)
     
330  dw=Alength*sin(ang)*Awidth
     totalw=2.0*PI*(1.0-cos(THa))
     Npart=Ma/Mp*dw/totalw
     Ntot=Ma/Mp
     partflux=obsfluxRs(tedyn, Ga, Ra, Npart, a, N3a, E3a, Ntot, Fca)
     EqTfluxRs=EqTfluxRs+partflux
  end do
  return
end function EqTfluxRs



!***********************************************************************
! FUNCTION interpolyne
! This function is used for RING ANGLE
!
! INPUT: 
! -
! -
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function anglength(deltA,THobs,alphA)

IMPLICIT NONE

DOUBLE PRECISION                         :: deltA, THobs, alphA
DOUBLE PRECISION                         :: temp, x1, x2, anglength

!***********************************************************************

  if( THobs >= (deltA+alphA)) then
     anglength=0.0
     return
  end if
  if( THobs <= (deltA-alphA)) then
     anglength=0.0
     return
  end if
  if( (deltA<alphA) .and. (THobs <= (alphA-deltA))) then
     anglength=2.0*PI
     return
  end if
  x1=cos(alphA)-cos(deltA)*cos(THobs)
  x2=sin(deltA)*sin(THobs)
  temp=x1/x2
  anglength=acos(temp)*2.0
  return
end function anglength




!***********************************************************************
! FUNCTION interpoly
! This function is used for INTERPOLATING OF GAM2,M2,M3
!
! INPUT: 
! - x  :
! - xx :
! - y :
! - k :
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function interpoly(x, xx, y, iX, k)

IMPLICIT NONE

DOUBLE PRECISION                         :: x, y(nvar, nstep), interpoly, xx(nstep)
INTEGER*8                                :: i, k, iX
DOUBLE PRECISION                         :: A0, A1, A2

!***********************************************************************

  i  = k
  if( k==1 )then
     k=k+1
  end if
  A0 = (x-10**xx(k))*(x-10**xx(k+1))/(10**xx(k-1)-10**xx(k))/          &
       (10**xx(k-1)-10**xx(k+1))
  A1 = (x-10**xx(k-1))*(x-10**xx(k+1))/(10**xx(k)-10**xx(k-1))/        &
       (10**xx(k)-10**xx(k+1))
  A2 = (x-10**xx(k-1))*(x-10**xx(k))/(10**xx(k+1)-10**xx(k-1))/        &
       (10**xx(k+1)-10**xx(k))
  interpoly = A0*y(iX,k-1)+A1*y(iX,k)+A2*y(iX,k+1)

  !print*,'t3(1)',y(3,1)
  !print*,'t3(2)',y(3,2)
  !print*,'t3(3)',y(3,3)

  !print*,'iX, k', iX,k
  !print*,'A0,A1,A2',A0,A1,A2
  !print*,'y(iX,k-1),y(iX,k),y(iX,k+1)',y(iX,k-1),y(iX,k),y(iX,k+1)
  !print*,'interpoly',interpoly
  return
  k = i
end function interpoly




!***********************************************************************
! FUNCTION interpolyne
! This function is used to INTERPOLATING N3,E3
!
! INPUT: 
! - x  :
! - xx :
! - y :
! - k :
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function interpolyne(x, xx, y, k)

IMPLICIT NONE

DOUBLE PRECISION                         :: x, y(nstep), interpolyne, xx(nstep)
INTEGER*8                                :: i, k, iX
DOUBLE PRECISION                         :: A0, A1, A2

!***********************************************************************

  i  = k
  if( k==1 )then
     k=k+1
  end if
  A0 = (x-10**xx(k))*(x-10**xx(k+1))/(10**xx(k-1)-10**xx(k))/          &
       (10**xx(k-1)-10**xx(k+1))
  A1 = (x-10**xx(k-1))*(x-10**xx(k+1))/(10**xx(k)-10**xx(k-1))/        &
       (10**xx(k)-10**xx(k+1))
  A2 = (x-10**xx(k-1))*(x-10**xx(k))/(10**xx(k+1)-10**xx(k-1))/        &
       (10**xx(k+1)-10**xx(k))
  interpolyne = A0*y(k-1)+A1*y(k)+A2*y(k+1)
  return
  k = i
end function interpolyne



!***********************************************************************
! FUNCTION obsfluxfs
! This 
!
! INPUT: 
! - te  :
! - gam :
! - Ra :
! - Npart :
! - ang : 
! - Ntot: 
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function obsfluxfs(te, gam2, Ra, Npart, ang, Ntot)

IMPLICIT NONE

DOUBLE PRECISION                         :: obsfluxfs, te, gam2, Ra, ang, Npart, Ntot
DOUBLE PRECISION                         :: trans, beta
DOUBLE PRECISION                         :: gm, gc, Pmax, Fmax, LamorFre, Bp, Fnu
DOUBLE PRECISION                         :: nuc, num, nua, nu
DOUBLE PRECISION                         :: gc_syn, YY, epsilon_rad !For Compton Effect
DOUBLE PRECISION                         :: dgec, Y, eta, gec2

DOUBLE PRECISION                         :: hgam, hgam1, hgam0, e2, p2,n2, th
integer(fgsl_int)                        :: i,l, j

!***********************************************************************

  ! Accurate adiabatic index hgam2 is through iteration:
  hgam = 5./3
  hgam0 = 5./3
  do j = 1,100
    n2 = (hgam*gam2+1.0)/(hgam-1.0)*en1(Ra)
    e2 = (gam2-1.0)*n2*mp*c*c
    p2 = (hgam-1.0)*e2
    th = p2/(n2*mp*c*c)

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
  !----End Adiabatic Index ------
  beta = dsqrt(1.D0-1.D0/gam2/gam2)
  trans= gam2*(1.D0-beta*cos(ang))
  gm   = epsefs*(gam2-1.D0)*mp/me*(p-2.0)/(p-1.0)+1.D0
  Bp   = sqrt(8.D0*pi*epsBfs*e2)
  gc   = 6.0D0*pi*me*c*(1.D0+redsh)/(sigmaT*gam2*Bp*Bp*te)
  LamorFre= e*Bp/(2.D0*pi*me*c)
  nu   = obsFre*trans*(1.D0+redsh)      !ÓîÖæµ±µØÏµµÄÆµÂÊ
  num  = 3.D0*gm*gm*LamorFre/2.D0
  nuc  = 3.D0*gc*gc*LamorFre/2.D0
  nua  = nu_absorb(Ra, n2, Bp, gm, gc, Ntot)
  Pmax = sqrt(3.0)*e**3.0*Bp/(me*c**2.0)*1.2 !*0.918
  Fmax = (1.D0+redsh)*Npart*Pmax/(4.D0*pi*D28**2.D0*1.D56)/trans**3.D0

  !-----------------------Compton Effect-------------------------------
  !----------------------Two method is 'ok'-----------------------------
  epse=epsefs
  epsB=epsBfs
  gc_syn = gc
  gamma_sync = gc_syn   !for module 'compton'
  gamma_m    = gm       !
  epsilon_rad = 1.D0
  YY = (dsqrt(1.D0+4.D0*epsilon_rad*epse/epsB)-1.0)/2.0D0
  gc = gc_syn/(1.D0+YY)
  if( gc_syn<gm .or. gc<gm )then
     gc = gc
  else
     gc = zbrent(eq_slow_gamma_c,1.D0,gamma_sync,1.D-8)
  end if

IF (ICS.EQV..TRUE.) THEN
  nuc = 1.5D0*gc*gc*LamorFre
ELSE
  nuc = 1.5D0*gc_syn*gc_syn*LamorFre
ENDIF

!--------------------µüŽú·œ·š,Ö»Ðè×¢ÊÍÈ¥µôŒŽ¿ÉÓÃ----------------------------
!  gc = gc_syn
!  dgec = 1.0
!  do while(dgec.ge.0.0001)
!     if(gm.le.gc)then
!	eta  = (gc/gm)**(2.-p)
!     else
!	eta  = 1.0
!     endif
!     Y    = (sqrt(1.+4.*eta*epsefs/epsBfs)-1.)/2.
!     gec2 = gc
!     gc   = gc_syn/(1+Y)
!     dgec = abs(gc-gec2)/gc
!  end do
!  nuc=1.5D0*gc*gc*LamorFre

  !--------------------------------------------------------------------
  if( (num>=nuc) .and. (nuc>=nua) )then
     !fast cooling, nu_a < nu_c < nu_m
     if(nu < nua)then
        Fnu = (nua/nuc)**(1.D0/3.D0)*(nu/nua)**2.D0
     elseif(nu < nuc)then
        Fnu = (nu/nuc)**(1.D0/3.D0)
     elseif(nu < num)then
        Fnu = (nu/nuc)**(-0.5D0)
     else
        Fnu = (num/nuc)**(-0.5D0)*(nu/num)**(-p/2.D0)
     endif
  elseif((num >= nua) .and. (nua >= nuc))then
     !     fast cooling, nu_c < nu_a < nu_m
     if(nu < nuc)then
        Fnu = (nuc/nua)**(2.5D0)*(nu/nuc)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     elseif(nu < num)then
        Fnu = (nu/nua)**(-0.5D0)
     else
        Fnu = (num/nua)**(-0.5D0)*(nu/num)**(-p/2.D0)
     endif
     Fnu = Fnu*(nua/nuc)**(-0.5D0)
  elseif((nua >= num) .and. (num >= nuc))then
     !     fast cooling, nu_c < nu_m < nu_a
     if(nu < nuc)then
        Fnu = (nuc/nua)**(2.5D0)*(nu/nuc)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     else
        Fnu = (nu/nua)**(-p/2.D0)
     endif
     Fnu = Fnu*(num/nuc)**(-0.5D0)*(nua/num)**(-p/2.D0)
  elseif((nua <= num) .and. (num <= nuc))then
     !     slow cooling, nu_a < nu_m < nu_c
     if(nu < nua)then
        Fnu = (nua/num)**(1.D0/3.D0)*(nu/nua)**2.D0
     elseif(nu < num)then
        Fnu = (nu/num)**(1.D0/3.D0)
     elseif(nu < nuc)then
        Fnu = (nu/num)**(-(p-1.D0)/2.D0)
     else
        Fnu = (nuc/num)**(-(p-1.D0)/2.D0) &
             *(nu/nuc)**(-p/2.D0)
     endif
  elseif((num <= nua) .and. (nua <= nuc))then
     !     slow cooling, nu_m < nu_a < nu_c
     if(nu < num)then
        Fnu = (num/nua)**(2.5D0)*(nu/num)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     elseif(nu < nuc)then
        Fnu = (nu/nua)**(-(p-1.D0)/2.D0)
     else
        Fnu = (nuc/nua)**(-(p-1.D0)/2.D0) &
             *(nu/nuc)**(-p/2.D0)
     endif
     Fnu = Fnu*(nua/num)**(-(p-1.D0)/2.D0)
  elseif((num <= nuc) .and. (nuc <= nua))then
     !     slow cooling, nu_m < nu_c < nu_a
     if(nu < nuc)then
        Fnu = (num/nua)**(2.5D0)*(nu/num)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     else
        Fnu = (nu/nua)**(-p/2.D0)
     endif
     Fnu = Fnu*(nuc/num)**(-(p-1.D0)/2.D0) &
          *(nua/nuc)**(-p/2.D0)
  else
     print*,'In subroutine emission, I find something impossible you should see it youself!'
     print*, 'nu_a:', nua, ' nu_c:', nuc, ' nu_m:', num
     pause
  endif
  obsfluxfs = Fnu*Fmax
  return 
end function obsfluxfs



!***********************************************************************
! FUNCTION obsfluxrs
! This 
!
! INPUT: 
! - te  :
! - gam :
! - Ra :
! - Npart :
! - ang : 
! - n3D : n3
! - e3D : e3
! - Ntot: 
! - fc
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

function obsfluxrs(te, gam, Ra, Npart, ang, n3D, e3D, Ntot, fc)

IMPLICIT NONE

DOUBLE PRECISION                         :: obsfluxrs, te, gam, Ra, ang, Ntot, n3D, e3D, Fnu, Npart
DOUBLE PRECISION                         :: trans, beta
DOUBLE PRECISION                         :: gm, gc, Pmax, LamorFre, ep, Bp, nu
DOUBLE PRECISION                         :: BpT, LamorFreT, nu_m3T, nu_c3T, nu_cut
DOUBLE PRECISION                         :: gc_syn, YY, epsilon_rad !For Compton Effect

DOUBLE PRECISION                         :: hgam, hgam0, p3, th, fc
integer(fgsl_int)                        :: i,l, j

!***********************************************************************

  ! Accurate adiabatic index hgam obtained by iteration
  hgam = 5./3
  hgam0 = 5./3
  do j = 1,100
    p3 = (hgam-1.0)*e3D*fc
    th = p3/(n3D*mp*c*c)
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
  !----End Adiabatic Index ------

  if (sigma > 0.) then
    epsBrs = (hgam-1.)*(fc-1.0)
  endif

  beta = dsqrt(1.D0-1.D0/gam/gam)
  trans= gam*(1.D0-beta*cos(ang))
  gm   = epsers*e3D/(n3D*me*c*c)*(p-2.D0)/(p-1.D0)+1.D0
  Bp   = sqrt(8.D0*pi*epsBrs*e3D)
  gc   = 6.0D0*pi*me*c/(sigmaT*gam*Bp*Bp*te)
  !print*,' pi, me, c, sigmaT, gam, Bp, te',pi, me, c, sigmaT, gam, Bp, te
  LamorFre= e*Bp/(2.D0*pi*me*c)
  nu=obsFre*trans*(1.D0+redsh)
  num = 1.5D0*gm*gm*LamorFre
  nuc = 1.5D0*gc*gc*LamorFre
  !print*,'3 Ra, n3D, Bp, gm, gc, Ntot',Ra, n3D, Bp, gm, gc, Ntot
  nua  = nu_absorb(Ra, n3D, Bp, gm, gc, Ntot)
  !print*,'nua',nua
  Pmax= (1.D0+redsh)*me*c*c*sigmaT*Bp/(3.D0*e)
  Fmax= Npart*Pmax/(4.D0*pi*D28**2*1.D56)/trans**3.D0

  !---------------------------Compton Effect----------------------------
  epse=epsers
  epsB=epsBrs
  gc_syn = gc
  gamma_sync = gc_syn   !for module 'compton'
  gamma_m    = gm       !
  epsilon_rad = 1.D0
  YY = (dsqrt(1.D0+4.D0*epsilon_rad*epse/epsB)-1.0)/2.0D0
  gc = gc_syn/(1.D0+YY)
  if( gc_syn<gm .or. gc<gm )then
     gc = gc
  else
     gc = zbrent(eq_slow_gamma_c,1.D0,gamma_sync,1.D-8)
  end if

IF (ICS.EQV..TRUE.) THEN
  nuc = 1.5D0*gc*gc*LamorFre
ELSE
  nuc = 1.5D0*gc_syn*gc_syn*LamorFre
ENDIF
  
  !---------------------------------------------------------------------
  !ÒÔÏÂµÄÓïŸäÊÇŽÓKobayashiÎÄÕÂÖÐÀŽµÄ£¬ÊµŒÊÉÏ²¢ÎÞ±ØÒª£¬·Ž¶ø»áÊ¹µÃÇúÏß±äµÃ²»¹â»¬£¡
!  BpT       = sqrt(8.D0*pi*epsBrs*e3(ncross))
!  LamorFreT = e*BpT/(2.D0*pi*me*c)
!  nu_m3T=(epsers*e3(ncross)/(n3(ncross)*me*c*c)*(p-2.D0)/       &
!       (p-1.D0)+1.D0)**2*LamorFreT
!  nu_c3T=(6.0D0*pi*me*c/(sigmaT*y(4,ncross)*BpT*BpT*y(3,ncross)))**2*&
!       LamorFreT

!  if( Ra>=10**xx(ncross) )then
!     nu_cut = nu_c3T*num/nu_m3T
!     nuc    = min(nuc, nu_cut)
!     if( nuc<nua )nua=nuc
!     if( nuc<num )num=nuc
!  end if

  !---------------------------------------------------------------------

  if( (num>=nuc) .and. (nuc>=nua) )then
     !fast cooling, nu_a < nu_c < nu_m
     if(nu < nua)then
        Fnu = (nua/nuc)**(1.D0/3.D0)*(nu/nua)**2.D0
     elseif(nu < nuc)then
        Fnu = (nu/nuc)**(1.D0/3.D0)
     elseif(nu < num)then
        Fnu = (nu/nuc)**(-0.5D0)
     else
        Fnu = (num/nuc)**(-0.5D0)*(nu/num)**(-p/2.D0)
     endif
  elseif((num >= nua) .and. (nua >= nuc))then
     !     fast cooling, nu_c < nu_a < nu_m
     if(nu < nuc)then
        Fnu = (nuc/nua)**(2.5D0)*(nu/nuc)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     elseif(nu < num)then
        Fnu = (nu/nua)**(-0.5D0)
     else
        Fnu = (num/nua)**(-0.5D0)*(nu/num)**(-p/2.D0)
     endif
     Fnu = Fnu*(nua/nuc)**(-0.5D0)
  elseif((nua >= num) .and. (num >= nuc))then
     !     fast cooling, nu_c < nu_m < nu_a
     if(nu < nuc)then
        Fnu = (nuc/nua)**(2.5D0)*(nu/nuc)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     else
        Fnu = (nu/nua)**(-p/2.D0)
     endif
     Fnu = Fnu*(num/nuc)**(-0.5D0)*(nua/num)**(-p/2.D0)
  elseif((nua <= num) .and. (num <= nuc))then
     !     slow cooling, nu_a < nu_m < nu_c
     if(nu < nua)then
        Fnu = (nua/num)**(1.D0/3.D0)*(nu/nua)**2.D0
     elseif(nu < num)then
        Fnu = (nu/num)**(1.D0/3.D0)
     elseif(nu < nuc)then
        Fnu = (nu/num)**(-(p-1.D0)/2.D0)
     else
        Fnu = (nuc/num)**(-(p-1.D0)/2.D0) &
             *(nu/nuc)**(-p/2.D0)
     endif
  elseif((num <= nua) .and. (nua <= nuc))then
     !     slow cooling, nu_m < nu_a < nu_c
     if(nu < num)then
        Fnu = (num/nua)**(2.5D0)*(nu/num)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     elseif(nu < nuc)then
        Fnu = (nu/nua)**(-(p-1.D0)/2.D0)
     else
        Fnu = (nuc/nua)**(-(p-1.D0)/2.D0) &
             *(nu/nuc)**(-p/2.D0)
     endif
     Fnu = Fnu*(nua/num)**(-(p-1.D0)/2.D0)
  elseif((num <= nuc) .and. (nuc <= nua))then
     !     slow cooling, nu_m < nu_c < nu_a
     if(nu < nuc)then
        Fnu = (num/nua)**(2.5D0)*(nu/num)**2.D0
     elseif(nu < nua)then
        Fnu = (nu/nua)**(2.5D0)
     else
        Fnu = (nu/nua)**(-p/2.D0)
     endif
     Fnu = Fnu*(nuc/num)**(-(p-1.D0)/2.D0) &
          *(nua/nuc)**(-p/2.D0)
  else
     print*,'In subroutine emission, I find something impossible, you should see it youself!'
     print*, 'nu_a:', nua, ' nu_c:', nuc, ' nu_m:', num
     pause
  endif
  obsfluxrs = Fnu*Fmax
  return 
end function obsfluxrs



!***********************************************************************
! Subroutine frequencyrs
! This subroutine computes RS frequency
!
! INPUT: 
! - te  :
! - gam :
! - R :
! - Ntot :
! - k
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

subroutine frequencyrs(te, gam, R, Ntot, k)

IMPLICIT NONE

DOUBLE PRECISION                         :: te, gam, R, Ntot, nn
DOUBLE PRECISION                         :: num, nuc, gm, gc,ee, B, trans, Pmax, hgam, nua, gmax, numax
DOUBLE PRECISION                         :: BpT, LamorFreT, nu_m3T, nu_c3T, nu_cut, LamorFre
DOUBLE PRECISION                         :: gc_syn, YY, epsilon_rad !For Compton Effect
INTEGER*8                                :: k

!***********************************************************************

  ee   = e3(k)
  nn   = n3(k)
  trans=gam*(1.D0-sqrt(1.D0-1.D0/gam/gam))*(1+redsh)
  B  = sqrt(8.D0*pi*ee*epsBrs)
  gm = epsers*ee/(nn*me*c*c)*(p-2.D0)/(p-1.D0)+1.D0
  gc = 6.0D0*pi*me*c*(1.D0+redsh)/(sigmaT*gam*B*B*te)
  LamorFre= e*B/(2.D0*pi*me*c)
  gmax = 1.D8/sqrt(B)
  num= gm**2.D0*3.D0*e*B/(4.D0*pi*me*c)
  nuc= gc**2.D0*3.D0*e*B/(4.D0*pi*me*c)
  numax=gmax**2.D0*3.D0*e*B/(4.D0*pi*me*c)
  nua = nu_absorb(R, nn, B, gm, gc, Ntot)
  Pmax= me*c*c*sigmaT*B/(3.D0*e)*0.7

  epse=epsers
  epsB=epsBrs

!  if( R>=10**xx(ncross) )then
!     BpT       = sqrt(8.D0*pi*epsBrs*e3(ncross))
!     LamorFreT = e*BpT/(2.D0*pi*me*c)
!     nu_m3T=(epsers*e3(ncross)/(n3(ncross)*me*c*c)*(p-2.D0)/(p-1.D0)+   &
!          1.D0)**2*LamorFreT
!     nu_c3T=(6.0D0*pi*me*c/(sigmaT*y(4,ncross)*BpT*BpT*y(3,ncross) &
!          ))**2*LamorFreT
!     nu_cut = nu_c3T*num/nu_m3T
!     nuc    = min(nuc, nu_cut)
!     if( nuc<nua )nua=nuc
!     if( nuc<num )num=nuc
!     gc = sqrt(nuc/1.5D0/gc/gc)
!     gm = sqrt(num/1.5D0/gm/gm)
!  end if

  gc_syn = gc
  gamma_sync = gc_syn   !for module 'compton'
  gamma_m    = gm       !
  epsilon_rad = 1.D0
  YY = (dsqrt(1.D0+4.D0*epsilon_rad*epse/epsB)-1.0)/2.0D0
  gc = gc_syn/(1.D0+YY)

  if( gc_syn<gm .or. gc<gm )then
     gc = gc
  else
     gc = zbrent(eq_slow_gamma_c,1.D0,gamma_sync,1.D-8)
  end if
  
IF (ICS.EQV..TRUE.) THEN
  nuc = 1.5D0*gc*gc*LamorFre
ELSE
  nuc = 1.5D0*gc_syn*gc_syn*LamorFre
ENDIF


!  if( R>=10**xx(ncross-1) )then
!     BpT       = sqrt(8.D0*pi*epsBrs*e3(ncross))
!     LamorFreT = e*BpT/(2.D0*pi*me*c)
!     nu_m3T=(epsers*e3(ncross)/(n3(ncross)*me*c*c)*(p-2.D0)/(p-1.D0)+   &
!          1.D0)**2*LamorFreT
!     nu_c3T=(6.0D0*pi*me*c/(sigmaT*y(4,ncross)*BpT*BpT*y(3,ncross) &
!          ))**2*LamorFreT
!     nu_cut = nu_c3T*num/nu_m3T
!     nuc    = min(nuc, nu_cut)
!     if( nuc<nua )nua=nuc
!     if( nuc<num )num=nuc
!  end if
!  write(200,*)log10(te),gam,R,log10(Ntot),gm,gc,gmax,num,nuc,nua,numax,Pmax
  write(outunit,*)te,gam,R,Ntot*mp,gm,gc,gmax,num,nuc,nua,numax,Pmax,B
end subroutine frequencyrs



!***********************************************************************
! Subroutine frequencyfs
! This subroutine computes FS frequency
! 
! INPUT: 
! - te  :
! - gam :
! - R :
! - Ntot :
!
! OUTPUT: 
!
! REMINDER: 
!***********************************************************************

subroutine frequencyfs(te, gam, R, Ntot, ee, nn)

IMPLICIT NONE

DOUBLE PRECISION                         :: te, gam, R, Ntot, nn
DOUBLE PRECISION                         :: num, nuc, gm, gc,ee, B, trans, Pmax
DOUBLE PRECISION                         :: nua, LamorFre, gmax, numax
DOUBLE PRECISION                         :: gc_syn, YY, epsilon_rad !For Compton Effect

!***********************************************************************

  trans=gam*(1.D0-sqrt(1.D0-1.D0/gam/gam))*(1+redsh)
  B  = sqrt(8.D0*pi*ee*epsBfs)
  gm = (gam-1.0D0)*epsefs*(p-2.0)*mp/((p-1.0)*me)+1.0       
  gc = 6.0D0*pi*me*c*(1.D0+redsh)/(sigmaT*gam*B*B*te)
  gmax= 1.D8/sqrt(B)
  LamorFre= e*B/(2.D0*pi*me*c)
  num= gm**2.D0*3.D0*e*B/(4.D0*pi*me*c)
  nuc= gc**2.D0*3.D0*e*B/(4.D0*pi*me*c)
  numax=gmax**2.D0*3.D0*e*B/(4.D0*pi*me*c)
  nua = nu_absorb(R, nn, B, gm, gc, Ntot)
  Pmax= me*c*c*sigmaT*B/(3.D0*e)*0.7

  epse=epsefs
  epsB=epsBfs
  gc_syn = gc
  gamma_sync = gc_syn   !for module 'compton'
  gamma_m    = gm       !
  epsilon_rad = 1.D0
  YY = (dsqrt(1.D0+4.D0*epsilon_rad*epse/epsB)-1.0)/2.0D0
  gc = gc_syn/(1.D0+YY)
  if( gc_syn<gm .or. gc<gm )then
     gc = gc
  else
     gc = zbrent(eq_slow_gamma_c,1.D0,gamma_sync,1.D-8)
  end if

IF (ICS.EQV..TRUE.) THEN
  nuc = 1.5D0*gc*gc*LamorFre
ELSE
  nuc = 1.5D0*gc_syn*gc_syn*LamorFre
ENDIF

write(outunit,*)te,gam,R,Ntot*mp,gm,gc,gmax,num,nuc,nua,numax,Pmax,B
end subroutine frequencyfs



!***********************************************************************
! FUNCTION EqTfluxFs
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

function eq_slow_gamma_c(gamma)

IMPLICIT NONE

DOUBLE PRECISION                         :: eq_slow_gamma_c,gamma

!***********************************************************************

eq_slow_gamma_c = epse/epsB*(gamma/gamma_m)**(4.D0-p)&
       +gamma*gamma_sync/gamma_m/gamma_m-(gamma_sync/gamma_m)**2.D0
end function eq_slow_gamma_c



END MODULE MOD_LIGHTCURVES

