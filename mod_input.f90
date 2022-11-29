MODULE MOD_INPUT

IMPLICIT NONE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ CONSTANTS +++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Speed of light [cm/s]
DOUBLE PRECISION, PARAMETER, PUBLIC :: c=299792458d2
! Fundamental charge [esu]
DOUBLE PRECISION, PARAMETER, PUBLIC :: e=4.8032068d-10
! Boltzmann constant [erg/K]
DOUBLE PRECISION, PARAMETER, PUBLIC :: k=1.380658d-16
! Planck constant [erg.s]
DOUBLE PRECISION, PARAMETER, PUBLIC :: h=6.6260755d-27
! Mass of the electron [g]
DOUBLE PRECISION, PARAMETER, PUBLIC :: me=9.1093897d-28
! Mass of the proton [g]
DOUBLE PRECISION, PARAMETER, PUBLIC :: mp=1.6726231d-24
! 1 eV in erg [erg]
DOUBLE PRECISION, PARAMETER, PUBLIC :: evtoerg=1.602177d-12
! pi
DOUBLE PRECISION, PARAMETER, PUBLIC :: pi=acos(-1.0d0)
! Thomson cross section [cm^2]
DOUBLE PRECISION, PARAMETER, PUBLIC :: sigmaT=8.0*pi*e*e*e*e/&
                                               (3.0*me*me*c*c*c*c)
DOUBLE PRECISION, PARAMETER, PUBLIC :: hp=6.582D-22

! Mass ratio and ion mass
DOUBLE PRECISION, PARAMETER, PUBLIC :: mass_ratio=me/(mp+me)
DOUBLE PRECISION, PARAMETER, PUBLIC :: mec2=me*c*c

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++ COMPUTATIONAL DOMAIN PARAMETERS ++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of time steps
INTEGER*8, PARAMETER, PUBLIC :: nstep=5000
! Number of differential varibles
INTEGER*8, PARAMETER, PUBLIC :: nvar=6
! Adiabatic index hgam is from relativistic temperature T=p/(nmc^2)
! Maximum iteration of adiabatic index derivation
INTEGER*8, PUBLIC :: MAX_ITERATIONS
! Reverse shock crossing step 
INTEGER*8, PUBLIC :: ncross
! Radius in log10
DOUBLE PRECISION, PUBLIC :: xx(nstep)
! Results of differential varibles
DOUBLE PRECISION, PUBLIC :: y(nvar,nstep)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++ DYNAMICS MODELS AND SOLVERS ++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Choose one of the dynamic models
!CHARACTER(LEN=10), PARAMETER, PUBLIC :: MODEL='FS_BEFORE_CROSS'
!CHARACTER(LEN=10), PARAMETER, PUBLIC :: MODEL='FS_AFTER_CROSS'
!CHARACTER(LEN=10), PARAMETER, PUBLIC :: MODEL='RS_AFTER_CROSS'

!CHARACTER(LEN=10), PARAMETER, PUBLIC :: INTEGRATION='SIMPSON'
!CHARACTER(LEN=10), PARAMETER, PUBLIC :: INTEGRATION='TRAPEZOID'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++ INITIAL CIRCUMSTANCE PARAMETERS ++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! Circum medium index, ISM: kk=0
DOUBLE PRECISION, PARAMETER, PUBLIC :: kk=0.0
DOUBLE PRECISION, PARAMETER, PUBLIC :: Aast=1.D-2
DOUBLE PRECISION, PARAMETER, PUBLIC :: n0=1.0
! particle number density ratio n4/n1
DOUBLE PRECISION, PUBLIC :: f
! particle number density of region 2
DOUBLE PRECISION, PUBLIC :: n2(nstep)
! particle number density of region 3
DOUBLE PRECISION, PUBLIC :: n3(nstep)
! particle number density of region 4
DOUBLE PRECISION, PUBLIC :: n4(nstep)
! energy density of region 2
DOUBLE PRECISION, PUBLIC :: e2(nstep)
! energy density of region 3
DOUBLE PRECISION, PUBLIC :: e3(nstep)
! pressure of region 2
DOUBLE PRECISION, PUBLIC :: p2(nstep)
! pressure of region 3
DOUBLE PRECISION, PUBLIC :: p3(nstep)
! adiabatic index of region 3
DOUBLE PRECISION, PUBLIC :: hgam, hgam3_dump(nstep)
! adiabatic index of region 2
DOUBLE PRECISION, PUBLIC :: hgam2_dump(nstep)

! Cosmological parameters
DOUBLE PRECISION, PARAMETER, PUBLIC :: H0 = 68.0D0
DOUBLE PRECISION, PARAMETER, PUBLIC :: OmegaL = 0.69D0
DOUBLE PRECISION, PARAMETER, PUBLIC :: OmegaM = 0.31D0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++ INITIAL DYNAMICS PARAMETERS ++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! Isotropic energy/1.D52
DOUBLE PRECISION, PARAMETER, PUBLIC :: Eiso=1.0
! Ejecta initial Lorentz factor
DOUBLE PRECISION, PARAMETER, PUBLIC :: eta=300.0
! Ejecta Lorentz factor
DOUBLE PRECISION, PARAMETER, PUBLIC :: gam4=eta
! relative Lorentz factor between region 3 and RS
DOUBLE PRECISION, PUBLIC :: u3s, u3s_dump(nstep)
! Ejecta shell width
DOUBLE PRECISION, PARAMETER, PUBLIC :: Delta0=1.D13

DOUBLE PRECISION, PARAMETER, PUBLIC :: p=2.5

! Energy loss in Region 2
DOUBLE PRECISION, PARAMETER, PUBLIC :: eps2=0.
! Energy loss in Region 3
DOUBLE PRECISION, PARAMETER, PUBLIC :: eps3=0.

! fraction of energy gained by leptons in region 2
DOUBLE PRECISION, PARAMETER, PUBLIC :: epsefs=0.1
! fraction of energy gained by mag field in region 2
DOUBLE PRECISION, PARAMETER, PUBLIC :: epsBfs=0.01
! fraction of energy gained by leptons in region 3
DOUBLE PRECISION, PARAMETER, PUBLIC :: epsers=0.1
! fraction of energy gained by mag field in region 3
DOUBLE PRECISION, PUBLIC :: epsBrs=0.01

! Jet Half opening angle
DOUBLE PRECISION, PARAMETER, PUBLIC :: thetaj=pi
! Viewing angle
DOUBLE PRECISION, PARAMETER, PUBLIC :: deltaA=0.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ MAGNETIC FIELD PARAMETERS ++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! Ejecta magnetization
DOUBLE PRECISION, PARAMETER, PUBLIC :: sigma=0.0
! before cross magnetization correction factor with derivative
DOUBLE PRECISION, PUBLIC :: fa, fa_dump(nstep), faa
DOUBLE PRECISION, PUBLIC :: fb, fb_dump(nstep)
DOUBLE PRECISION, PUBLIC :: fc, fc_dump(nstep), fca

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ VARIABLE PARAMETERS ++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! Mass of ejecta
DOUBLE PRECISION, PARAMETER, PUBLIC :: mej=Eiso*1.D52/eta/c/c/(1.+sigma)

DOUBLE PRECISION, PARAMETER, PUBLIC :: CorJet = (1.D0-cos(thetaj))/2.D0

! Initial value of time in source frame
DOUBLE PRECISION, PUBLIC :: t0
! Initial value of time in observer frame
DOUBLE PRECISION, PUBLIC :: tini
! Sedov length
DOUBLE PRECISION, PUBLIC :: SedL
! Deceleration radius
DOUBLE PRECISION,PUBLIC :: Rdec
! Initial radius
DOUBLE PRECISION,PUBLIC :: R0 = Delta0*eta
! Radius lower limit
DOUBLE PRECISION,PUBLIC :: x1
! Radius upper limit
DOUBLE PRECISION, PUBLIC :: x2


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++ INITIAL SPECTRUM PARAMETERS ++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! observed frequency
DOUBLE PRECISION, PUBLIC :: obsfre=4.67D14 !optic red
!DOUBLE PRECISION, PUBLIC :: obsfre=5.D9 ! radio
DOUBLE PRECISION, PUBLIC :: D28

! character_frequency
DOUBLE PRECISION, PUBLIC :: num, nuc, nua, Fmax

! redshift: optic red 2.D14, radio 5.D9
DOUBLE PRECISION, PARAMETER, PUBLIC :: redsh=1.
!DOUBLE PRECISION, PUBLIC :: flux(nstep)
!DOUBLE PRECISION, PUBLIC :: nuc(nstep)
!DOUBLE PRECISION, PUBLIC :: num(nstep)

! Self-absorption
DOUBLE PRECISION, PUBLIC :: C_length,C_N_gamma, C_B, C_gamma_2, C_nu_2, C_p, C_tau_0

! Switch ON/OFF Include Inverse Compton Scattering
LOGICAL, PARAMETER, PUBLIC :: ICS=.TRUE.
! Inverse Compton 
DOUBLE PRECISION, PUBLIC :: gamma_sync, gamma_m, epse, epsB

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++ ANALYSIS AND OUTPUT +++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Output file name and index
character(len=200) :: filename
integer, parameter :: outunit = 44

! Switch ON/OFF reading dynamics from existed files
LOGICAL, PARAMETER, PUBLIC :: READ_DYNAMICS_FROM_FILE=.false.

! Switch ON/OFF writing data for dynamics
LOGICAL, PARAMETER, PUBLIC :: WRITE_DYNAMICS=.TRUE.
LOGICAL, PARAMETER, PUBLIC :: WRITE_MAGNETIZATION=.TRUE.
LOGICAL, PARAMETER, PUBLIC :: WRITE_DENSITY=.TRUE.
LOGICAL, PARAMETER, PUBLIC :: WRITE_FORWARD_SHOCK_FREQUENCY=.false.
LOGICAL, PARAMETER, PUBLIC :: WRITE_FORWARD_SHOCK_LIGHTCURVES=.false.
LOGICAL, PARAMETER, PUBLIC :: WRITE_REVERSE_SHOCK_FREQUENCY=.false.
LOGICAL, PARAMETER, PUBLIC :: WRITE_REVERSE_SHOCK_LIGHTCURVES=.false.

END MODULE MOD_INPUT
