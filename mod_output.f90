MODULE MOD_OUTPUT

USE FGSL
USE MOD_INPUT
USE MOD_INITIAL

IMPLICIT NONE

PRIVATE

PUBLIC :: DUMP_DYNAMICS ! Dynamics results
PUBLIC :: DUMP_MAGNETIZATION ! Lightcurves results
PUBLIC :: DUMP_DENSITY ! Density and pressure results
PUBLIC :: PARAMETERS ! Write parameter values to a file

CONTAINS


!***********************************************************************
! Subroutine PARAMETERS
! This subroutine writes parameter values to a file
!
!***********************************************************************

SUBROUTINE PARAMETERS()

IMPLICIT NONE

OPEN(9,FILE="./data/parameters.dat")

WRITE(9,*) "============================================"
WRITE(9,*) "NSTEP number of iteration steps             ", NSTEP
WRITE(9,*) "NVAR number of differential variables       ", NVAR
WRITE(9,*) "--------------------------------------------"
WRITE(9,*) "Eiso isotropic energy/1.D52         ", Eiso
WRITE(9,*) "eta initial Lorentz factor          ", eta
WRITE(9,*) "Delta0 ejecta shell width           ", Delta0
WRITE(9,*) "eps2 radiation loss in FS           ", eps2
WRITE(9,*) "eps2 radiation loss in RS           ", eps3
WRITE(9,*) "epsefs FS lepton dissipation gaining", epsefs
WRITE(9,*) "epsers RS lepton dissipation gaining", epsers
WRITE(9,*) "epsBfs FS mag field energy fraction ", epsBfs
WRITE(9,*) "epsBrs RS mag field energy fraction ", epsBrs
WRITE(9,*) "--------------------------------------------"
WRITE(9,*) "sigma ejecta magnetization          ", sigma
WRITE(9,*) "--------------------------------------------"
!WRITE(9,*) "sigma magnetization               ", 
!WRITE(9,*) "rhoc nominal Larmor radius        ", real(rhoc)
WRITE(9,*) "Rdec deceleration radius            ", Rdec
WRITE(9,*) "x1 radius lower limit in log        ", x1
WRITE(9,*) "x2 radius upper limit in log        ", x2
WRITE(9,*) "mej initial mass                    ", mej
WRITE(9,*) "gam4 initial Lorentz factor         ", gam4
WRITE(9,*) "t0 initial time                     ", t0
WRITE(9,*) "theta_j jet half opening angle      ", thetaj
WRITE(9,*) "--------------------------------------------"
WRITE(9,*) "kk circum medium index          ", kk
WRITE(9,*) "n0 electron number density      ", n0
WRITE(9,*) "Aast electron number density    ", Aast
WRITE(9,*) "--------------------------------------------"
WRITE(9,*) "redsh redshift                  ", redsh
IF (ICS.EQV..TRUE.) THEN 
  WRITE(9,*) "ICS include                      "
ELSE
  WRITE(9,*) "ICS not include                  " 
END IF
WRITE(9,*) "obsfre observed frequency       ", obsfre
WRITE(9,*) "luminosity distance             ", D28
WRITE(9,*) "============================================"
 CLOSE(9)

END SUBROUTINE PARAMETERS



!***********************************************************************
! Subroutine WRITE_DYN
!
! INPUT:
! 
!
! OUTPUT: 
! - dynfs: t2, gamma2, m2, r2, Ek2, U2, mom2
! - t2
! - gamma2
! - m2
! - r2
! - Ek2 : kinetic energy
! - U2 : thermal energy
! - mom2 : momentum
!
! - dynrs: t3, gamma3, m3, r3, Ek3, U3, EB3, mom3, mom4, E4
! - t3
! - gamma3
! - m3
! - r3
! - Ek3 : kinetic energy
! - U3 : thermal energy
! - EB3 : magnetic energy
! - mom3
! - mom4
! - E4
!
!***********************************************************************

SUBROUTINE DUMP_DYNAMICS
IMPLICIT NONE

DOUBLE PRECISION                         :: gam2, gam3, gam34, u2,u3,u4
DOUBLE PRECISION                         :: mom2, mom3, mom4
DOUBLE PRECISION                         :: E2, E3, E4
DOUBLE PRECISION                         :: Eint2com, Eint2, Eint3com, Eint3
DOUBLE PRECISION                         :: Ek2,Ek3
DOUBLE PRECISION                         :: EB3com, EB3
DOUBLE PRECISION                         :: m2,m3
INTEGER                                  :: i

!===================================================================
! PARALLEL DUMPING OF THE DATA
!===================================================================

IF (WRITE_DYNAMICS.EQV..TRUE.) THEN
  write(filename,'("./data/dynamics/dynfs.dat")')
  open(11,file=filename, form='formatted')
  write(filename,'("./data/dynamics/dynrs.dat")')
  open(21,file=filename, form='formatted')
  do i=1,nstep
    IF(xx(i).EQ.0)THEN
      EXIT
    ELSE
      gam2 = y(2,i)
      m2 = y(5,i)
      gam3 = y(4,i)
      m3 = y(6,i)
      u2 = sqrt(gam2*gam2-1.0)
      u3 = sqrt(gam3*gam3-1.0)
      u4 = sqrt(gam4*gam4-1.0)
      gam34 = gam4*gam3-u4*u3       
      Ek2 = (gam2-1.)*m2*c*c 
      Eint2com = (gam2-1.)*m2*c*c
      Eint2 = (1.-eps2)*gam2*Eint2com
      E2 = Ek2+Eint2
      Ek3 = (gam3-1.)*m3*c*c
      Eint3com = fa_dump(i)*(gam34-1.)*m3*c*c
      Eint3 = (1.-eps3)*gam3*Eint3com
      EB3com = (hgam3_dump(i)-1.)*(fc_dump(i)-1.)*Eint3com
      EB3 = gam3*EB3com
      E3 = Ek3+Eint3+EB3
      E4 = (gam4-1.)*(1.+sigma)*(mej-m3)*c*c
      mom2 = (m2+Eint2com/c/c)*u2 ! Piran 1999 Eq. 38
      mom3 = (m3+Eint3com/c/c+EB3com/c/c)*u3
      mom4 = (mej-m3)*u4*(1.+sigma)
      write(11,*)y(1,i),y(2,i),y(5,i),10**xx(i),Ek2, Eint2, mom2
      write(21,*)y(3,i),y(4,i),y(6,i),10**xx(i),Ek3, Eint3, EB3, mom3, mom4, E4
    ENDIF
  end do
  close(11)
  close(21)
ENDIF

END SUBROUTINE DUMP_DYNAMICS



!***********************************************************************
! Subroutine WRITE_MAGNETIZATION
!
! INPUT:
! 
! - t3, gamma3, gamma34, u3s, fa, fb, fc, hgam3, ncross
!
! OUTPUT: 
!***********************************************************************

SUBROUTINE DUMP_MAGNETIZATION
IMPLICIT NONE

DOUBLE PRECISION                         :: bet4, bet3, gam3, gam34
INTEGER                                  :: i

!===================================================================
! PARALLEL DUMPING OF THE DATA
!===================================================================

IF (WRITE_DYNAMICS.EQV..TRUE.) THEN
  write(filename,'("./data/magnetization/mag.dat")')
  open(11,file=filename, form='formatted')
  do i=1,nstep
     IF(xx(i).EQ.0)THEN
       EXIT
     ELSE
       gam3   = y(4,i)
       bet4   = sqrt(1.D0-1.D0/gam4**2)
       bet3   = sqrt(1.D0-1.D0/gam3**2)
       gam34  = (1.D0-bet3*bet4)*gam4*gam3
       write(11,*)y(3,i),gam3,gam34,u3s_dump(i),fa_dump(i),fb_dump(i),fc_dump(i), hgam3_dump(i), ncross
    ENDIF
  end do
  close(11)
ENDIF

END SUBROUTINE DUMP_MAGNETIZATION



!***********************************************************************
! Subroutine WRITE_DENSITY
!
! OUTPUT:
! 
! - n1, n2, p2, e2, hgam2
! - n4, n3, p3, e3, hgam3
!
!***********************************************************************

SUBROUTINE DUMP_DENSITY
IMPLICIT NONE

INTEGER                                  :: i

!***********************************************************************

!===================================================================
! PARALLEL DUMPING OF THE DATA
!===================================================================

IF (WRITE_DENSITY.EQV..TRUE.) THEN
  write(filename,'("./data/densities/denfs.dat")')
  open(11,file=filename, form='formatted')
  write(filename,'("./data/densities/denrs.dat")')
  open(21,file=filename, form='formatted')
  do i=1,nstep
    IF(xx(i).EQ.0)THEN
      EXIT
    ELSE
      write(11,*) en1(10**xx(i)),n2(i),p2(i),e2(i),hgam2_dump(i)
      write(21,*) n4(i), n3(i), p3(i), e3(i), hgam3_dump(i)
    ENDIF
  end do
  close(11)
  close(21)
ENDIF

END SUBROUTINE DUMP_DENSITY


END MODULE MOD_OUTPUT
