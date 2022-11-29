!=======================================================================
!
!                              TORCH 
!
! This is a GRB afterglow dynamics and radiation package.
! 
! FEATURE:
! - read/calculate dynamics
! - energy conserved GRB afterglow dynamics
! - magnetized ejecta
! - synchrotron radiation, self-absorption
! - inverse compton scattering
! 
! UPDATE DATE: 01 Aug 2020
!
! DEPEDENCY:
! GSL v2.7.1, FGSL v1.5.0 (1/5/2022)
! 
! UNDER DEVELOPMENT:
! - 
!
! AUTHORS:
!  - Prof. Xuewen Liu liuxuew@scu.edu.cn 
!  - Mr.   Qiang Chen chen@camk.edu.pl
!
! REMEMBER:
! - external library FGSL is needed. 
!
! USAGE: 
! - Settings: set initial parameters in mod_input.f90
! - Compile : execute ./run.sh to compile with Makefile
! - Manually: gfortran -I/usr/local/include/fgsl *.f90 &
!                             -lfgsl -lgsl -lgslcblas -lm
!=======================================================================

PROGRAM MAIN

USE MOD_INPUT
use MOD_INITIAL
USE MOD_MAGNETIZATION
USE MOD_READ
USE MOD_DYNAMICS
USE MOD_LIGHTCURVES
USE MOD_OUTPUT

use fgsl
use, intrinsic :: iso_c_binding

IMPLICIT NONE

!=======================================================================
!
! STEP 1: DECLARATION OF VARIABLES
!         LOCAL VARIABLES (VALID IN THE DOMAIN OF EACH PROCESS)
!
!=======================================================================

DOUBLE PRECISION :: CPU_T0,CPU_T1
DOUBLE PRECISION :: vstart(nvar)
DOUBLE PRECISION :: bet4
DOUBLE PRECISION :: fluxsyn ,flux_fs, flux_rs, sumsyn, hlogfre, hfre, logfre
integer*8        :: i, j

real( kind = fgsl_double )   :: result, error
integer( kind = fgsl_size_t) :: neval
integer( kind = fgsl_int)    :: i_int
type( fgsl_function ) :: func


CALL CPU_TIME(CPU_T0)
PRINT*,'INITIALIZING SIMULATION'
PRINT*,'  MAGNETIZATION = ',sigma
MAX_ITERATIONS = 0

if (kk ==0)then
  SedL  = ((3.D0-kk)*Eiso*1.D52/(4.D0*pi*n0*mp*c**2.0))**(1./(3-kk))
else
  SedL  = ((3.D0-kk)*Eiso*1.D52/(4.D0*pi*Aast*3.D35*mp*c**2.0))**(1./(3-kk))
endif

bet4  = sqrt(gam4*gam4-1.0)/gam4
Rdec  = SedL*eta**(-2./(3-kk))
t0    = R0/c ! source frame
tini  = (1+redsh)*t0*(1.-bet4) ! observer frame
f = mej/(4.D0*pi*R0*R0*gam4*max(Delta0,R0/gam4**2)*mp)/en1(R0)

vstart(1) = tini
vstart(2) = zbrent(eq_gamma_3, 1.D0, gam4,1.D-8)
vstart(3) = vstart(1)
vstart(4) = vstart(2)
vstart(5) = 4.D0/3.D0*pi*R0*R0*R0*en1(R0)*mp
vstart(6) = 0.D0

x1        = log10(R0)  
x2        = log10(1.D6*R0)

if (Delta0 > SedL/(gam4**(8./3)))then
    print*,'  THICK SHELL, Tx = ',Delta0/c
else
    print*,'  THIN SHELL, Tx = ',SedL/(2**(2./3)*c*gam4**(8./3))
end if

func = fgsl_function_init( Dc, c_null_ptr )
i_int = fgsl_integration_qng ( func,          &
                       0.0_fgsl_double,       &
                       redsh,                 &
                       1e-9_fgsl_double,      & 
                       1e-9_fgsl_double,      &   
                       result, error, neval ) 
! result = Dc, for Omega_k = 0, DM = Dc 
! DL = (1+z)DM ! unit pc
D28 = 3.085678*1.D18*(1+redsh)*result/1.D28

IF (ICS.EQV..TRUE.) THEN
  print*, '  INCLUDE Inverse Compton Scattering'
ENDIF

CALL PARAMETERS()

!=======================================================================
!
! STEP 2: READ DYNAMICS, OR SOLVE ODEs WITH RK4 ALGORITHM
!
!=======================================================================

IF (READ_DYNAMICS_FROM_FILE.EQV..TRUE.) THEN
  print*, 'READING DYNAMICS FROM FILE'
  CALL READ_DYNAMICS()
  CALL READ_MAGNETIZATION()
  CALL READ_DENSITY()
ELSE
  print*, 'CALCULATING DYNAMICS'
  CALL RKDUMB(VSTART,NVAR,X1,X2,NSTEP,FS_BEFORE_CROSS,FS_AFTER_CROSS)
  call RS_AFTER_CROSS(nstep, x1, x2)
ENDIF

!=======================================================================
!
! STEP 3: FORWARD SHOCK SPECTRUM AND LIGHTCURVES
!
!=======================================================================

print*, 'INITIALIZING RADIATION SECTION'
IF (WRITE_FORWARD_SHOCK_FREQUENCY.EQV..TRUE.) THEN
  print*, 'START FS Frequencies'
  !write(filename,'("./data/lightcurves/frequencyfs",i0,".dat")') filenum
  write(filename,'("./data/lightcurves/frequencyfs.dat")')
  open(unit=outunit,file=filename, form='formatted')
  j=0
  do i=1,nstep
     call frequencyfs(y(1,i),y(2,i),10**xx(i),y(5,i)/mp, e2(i), n2(i))
     j=j+1
  end do
  close(outunit)
ENDIF

IF (WRITE_FORWARD_SHOCK_LIGHTCURVES.EQV..TRUE.) THEN
  print*, 'START FS Light Curve! Please wait!'
  write(filename,'("./data/lightcurves/fluxfs.dat")')
  open(unit=outunit,file=filename, form='formatted')
  do i=4,nstep, 40
     flux_fs=EqTfluxFs(y(1,i))
     if( obsfre<1.D18 )then
        write(outunit,*)y(1,i),flux_fs*CorJet
     else
        sumsyn = 0.D0
        hlogfre = (log10(10.D0*2.417D17)-log10(0.3D0*2.417D17))/20.D0
        logfre = log10(0.3*2.417D17)
        do j = 1, 20
           obsfre = 10**logfre
           hfre   = 10**(logfre+hlogfre)-10**logfre
           fluxsyn = flux_fs
           sumsyn = sumsyn + fluxsyn*hfre
           logfre = logfre+hlogfre
        end do
        write(outunit,*)y(1,i),sumsyn*CorJet
     end if
     !if( y(1,i)>3.D5 )then
     !   exit
     !end if
  end do
ENDIF

!=======================================================================
!
! STEP 4: REVERSE SHOCK SPECTRUM AND LIGHTCURVES
!
!=======================================================================

IF (WRITE_REVERSE_SHOCK_FREQUENCY.EQV..TRUE.) THEN
  print*, 'START RS Frequencies'
  write(filename,'("./data/lightcurves/frequencyrs.dat")')
  open(unit=outunit,file=filename, form='formatted')
  j=0
  do i=2,nstep
     call frequencyrs(y(3,i),y(4,i),10**xx(i),y(6,i)/mp,i)
     j=j+1
     if( y(3,i)>1.D5 .or. y(4,i)<1. )then
        exit
     end if
  end do
  close(outunit)
ENDIF

IF (WRITE_REVERSE_SHOCK_LIGHTCURVES.EQV..TRUE.) THEN
  print*, 'START RS Light Curve! Please wait!'
  write(filename,'("./data/lightcurves/fluxrs.dat")')
  open(unit=outunit,file=filename, form='formatted')
    do i=4,nstep, 16
    ! nstep=10000, successful test: N = 25
    ! nstep=5000, N <= 16, why?
     flux_rs=EqTfluxRs(y(3,i))
      write(outunit,*)y(3,i),flux_rs*CorJet  !unit in Jule, = 10.**23*10**3 mJy
     if( y(3,i)>2.D6 .or. y(4,i)<1. )then
        exit
     end if
    end do
  close(outunit)
ENDIF

!=======================================================================
!
! STEP 5: DATA ANALYSIS AND OUTPUT
!
!=======================================================================

print*, 'DUMPING DATA'
IF (READ_DYNAMICS_FROM_FILE.EQV..FALSE.) THEN
  IF (WRITE_DYNAMICS.EQV..TRUE.) THEN
    CALL DUMP_DYNAMICS()
    CALL DUMP_DENSITY()
    CALL DUMP_MAGNETIZATION()
  ENDIF
ENDIF

PRINT*,'ADIABATIC INDEX ITERATIONS (<=100):', MAX_ITERATIONS
CALL CPU_TIME(CPU_T1)
print*,'HALLELUJAH! :-) ELAPSED TIME',CPU_T1-CPU_T0
END PROGRAM MAIN
