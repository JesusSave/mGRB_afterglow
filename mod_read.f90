MODULE MOD_READ

USE MOD_INPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: READ_DYNAMICS ! Dynamics results
PUBLIC :: READ_MAGNETIZATION ! Lightcurves results
PUBLIC :: READ_DENSITY ! Density and pressure results

CONTAINS



!***********************************************************************
! Subroutine READ_DYN
!
!
! INPUT: 
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
! REMEMBER:
! - radius is calculated again from x1, x2 and nstep
!***********************************************************************

SUBROUTINE READ_DYNAMICS
IMPLICIT NONE

DOUBLE PRECISION                         :: r2, Ek2, Eint2, mom2
DOUBLE PRECISION                         :: r3, Ek3, Eint3, EB3, mom3, mom4, E4
DOUBLE PRECISION                         :: h, x
INTEGER                                  :: i

!===================================================================

  OPEN (unit=11, file="./data/dynamics/dynfs.dat", status='old',    &
        access='sequential', form='formatted', action='read' )
  OPEN (unit=21, file="./data/dynamics/dynrs.dat", status='old',    &
        access='sequential', form='formatted', action='read' )

  xx(1)=x1
  x=x1
  h=(x2-x1)/nstep

  do i=1,nstep-1
    x=x+h
    xx(i+1)=x
  enddo

  do i=1,nstep
    READ(11,*)y(1,i),y(2,i),y(5,i),r2,Ek2, Eint2, mom2
    READ(21,*)y(3,i),y(4,i),y(6,i),r3,Ek3, Eint3, EB3, mom3, mom4, E4
  end do
  CLOSE(11)
  CLOSE(21)

END SUBROUTINE READ_DYNAMICS



!***********************************************************************
! Subroutine READ_MAGNETIZATION
!
! INPUT:
! 
! - t3, gamma3, gamma34, u3s, fa, fb, fc, hgam3
!
! OUTPUT: 
!***********************************************************************

SUBROUTINE READ_MAGNETIZATION
IMPLICIT NONE

DOUBLE PRECISION                         :: gam3, gam34
INTEGER                                  :: i

!===================================================================

  OPEN (unit=11, file="./data/magnetization/mag.dat", status='old',    &
        access='sequential', form='formatted', action='read' )
  do i=1,nstep
      READ(11,*)y(3,i),gam3,gam34,u3s_dump(i),fa_dump(i),fb_dump(i),fc_dump(i), hgam3_dump(i)
  end do
  CLOSE(11)

END SUBROUTINE READ_MAGNETIZATION



!***********************************************************************
! Subroutine READ_DENSITY
!
! INPUT:
! 
! - n1, n2, p2, e2, hgam2
! - n4, n3, p3, e3, hgam3
!
!***********************************************************************

SUBROUTINE READ_DENSITY
IMPLICIT NONE

DOUBLE PRECISION                         :: n1
INTEGER                                  :: i

!***********************************************************************

!===================================================================

  OPEN (unit=11, file="./data/densities/denfs.dat", status='old',    &
        access='sequential', form='formatted', action='read' )
  OPEN (unit=21, file="./data/densities/denrs.dat", status='old',    &
        access='sequential', form='formatted', action='read' )
  do i=1,nstep
      READ(11,*)n1,n2(i),p2(i),e2(i),hgam2_dump(i)
      READ(21,*)n4(i), n3(i), p3(i), e3(i), hgam3_dump(i)
  end do
  CLOSE(11)
  CLOSE(21)

END SUBROUTINE READ_DENSITY


END MODULE MOD_READ
