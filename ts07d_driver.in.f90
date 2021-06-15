!======================================================================================================
! July 2017: G.K.Stephens (Grant.Stephens@jhuapl.edu)
! This program demonstrates how to setup and evaluate the Tsyganenko and Sitnov [2007] (TS07D)
! empirical geomagnetic field model and how to setup and evaluate the IGRF-12 model included
! in Tsyganenko's GEOPACK package.
!
! This July 2017 update incorporates an upgraded version of the Bessel function evaluator,
! provided by Jay Albert, which significantly speeds up the model. We thank Jay Albert for these
! contributions.
!
! Additionally, the model subroutine was updated to be a double precision subroutine.
! To indicate the update, the subroutine was renamed from EXTMODEL to TS07D_JULY_2017.
!
! To compile this code with GFortran, use:
! gfortran ts07d_geopack_example_july2017update.for -o ts07d_geopack_example_july2017update
! and to run: ./ts07d_geopack_example_july2017update
!======================================================================================================
Program   TS07D
implicit none

!     Inputs to the model
INTEGER   IOPT
REAL*8  PARMOD(10), XGSM,YGSM,ZGSM

!     Input parameters that will be passed through the common block
INTEGER    NTOT
PARAMETER (NTOT=101)
REAL*8 PDYN,PARAMS(NTOT),TSS,TSO,TSE
COMMON /PARAM/ PARAMS
COMMON /INPUT/  PDYN
COMMON /TSS/ TSS(80,5)
COMMON /TSO/ TSO(80,5,4)
COMMON /TSE/ TSE(80,5,4)

!     Set up the Geopack common block
REAL*8 AAA,BBB,PSI ! the dipole tilt
COMMON /GEOPACK1/ AAA(15),PSI,BBB(18)

!     Output parameters for TS07D and IGRF models
REAL*8 BXGSM,BYGSM,BZGSM, HXGSM,HYGSM,HZGSM

!     STATICDIR is the directory that contains the static shielding coefficients
CHARACTER(1024) :: STATICDIR
!    VARIABLEDIR is the directory where the variable coefficient directories will be
CHARACTER(1024) :: VARIABLEDIR

!     The meta data included in the .par file
REAL*8 COEFF_Q, COEFF_B_RMS,TILT
INTEGER M_INX,N_INX

CHARACTER(1024) :: FILENAME

character(16) :: argv
integer :: ierr

REAL*8 VXGSE,VYGSE,VZGSE,PI
INTEGER KK,IREAD,KREAD
INTEGER IYEAR,IDAY,IHOUR,MIN,ISEC

staticdir = "@ts07d_static_dir@"  !'./ts07d_tail_par/'
variabledir = "@ts07d_variable_dir@" ! './ts07_coeff/'

!     The following reads the static coefficients that are used for the Equatorial Field
!     shielding fields. These only need to be loaded once as they are common for all times.
      DO 1001 IREAD=1,5
        WRITE(filename,'(A,A,I0,A)') TRIM(STATICDIR),'/tailamebhr',IREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSS(KK,IREAD),KK=1,80)
 200  FORMAT(G17.10)
 1001 CLOSE(1)

      DO 1002 IREAD=1,5
        DO 1003 KREAD=1,4
        WRITE(filename,'(A,A,I0,I0,A)') TRIM(STATICDIR),'/tailamhr_o_',IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSO(KK,IREAD,KREAD),KK=1,80)
 1003   CONTINUE
 1002 CLOSE(1)

      DO 1004 IREAD=1,5
        DO 1005 KREAD=1,4
        WRITE(filename,'(A,A,I0,I0,A)') TRIM(STATICDIR),'/tailamhr_e_',IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSE(KK,IREAD,KREAD),KK=1,80)
 1005   CONTINUE
 1004 CLOSE(1)

!     Now, initialize the TS07D and IGRF models for the supplied date and time. If multiple times
!     are wanted, the following must be recalled for the new date and time.

!     The date and time is specified by the year, day of year (1-366), hour (0-23), minute (0-59),
!     and the second (0-59), note that the second field is not used for accessing the variable
!     parameter and coefficients file.
IYEAR = 2015
IDAY = 75
IHOUR = 12
MIN = 35
ISEC = 0

!     Use the date and time fields to construct the variable parameter and coefficients file.
!     Or, if preferred, you could hardcode the path as:
!     filename = '/Users/username/Downloads/2015_075/2015_075_12_35.par'
WRITE(filename,'(A,I0.4,A,I0.3,A,I0.4,A,I0.3,A,I0.2,A,I0.2,A)') &
      TRIM(VARIABLEDIR), IYEAR,'_',IDAY,'/',IYEAR,'_',IDAY,'_',IHOUR,'_',MIN,'.par'

!     Now read variable parameter and coefficients file, the variable parameters and coefficients
!     are stored in the array PARAMS. Additionaly, the dynamic pressure and the dipole tilt angle
!     are stored in PDYN and TILT respectively. All of these are variable or time-dependent values.
!     The other values stored in the file are not currently being used.
OPEN (UNIT=1,FILE=filename,action='read') ! open the filed
READ (1,100) (PARAMS(IREAD),IREAD=1,NTOT) ! read the variable coefficients and parameters
READ (1,101) COEFF_Q ! the Q factor, related to chi squared, a measure of the goodness of fit
READ (1,101) COEFF_B_RMS
READ (1,102) M_INX ! the number of azimuthal expansions in the equatorial field module (M=4)
READ (1,102) N_INX ! the number of radial expansions in the equatorial field module (N=5)
READ (1,101) PDYN ! the dynamic pressure for this time
READ (1,101) TILT ! the dipole tilt for this time
 100  FORMAT(G15.6)
 101  FORMAT(7x,G15.6)
 102  FORMAT(7x,I15)
CLOSE(1)

!     To get the IGRF model in GSM coordinates, set the solar wind velocity vector to the following
VXGSE=-400.d0
VYGSE=   0.d0
VZGSE=   0.d0
!     initializes the Geopack package for the supplied date and time
call RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VXGSE,VYGSE,VZGSE)


!     Now the TS07D and the IGRF model can be evaluated

!     or just code in the values, e.g. XGSM = 5.0d0, units are in Earth radii (1 RE=6371.2 KM).
if (command_argument_count() /= 3) error stop '  Enter S/C position: XGSM,YGSM,ZGSM '
call get_command_argument(1, argv)
read(argv, '(F9.3)', iostat=ierr) xgsm
call get_command_argument(2, argv)
read(argv, '(F9.3)', iostat=ierr) ygsm
call get_command_argument(3, argv)
read(argv, '(F9.3)', iostat=ierr) zgsm

if (ierr/=0) error stop '  Enter S/C position: XGSM,YGSM,ZGSM '


!     Evaluate the IGRF model in GSM coordinates
CALL IGRF_GSW_08 (XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)

IOPT = 0 ! the whole external magnetic field model (almost always what you want)
!     Now evaluate the TS07D model, note that we have two variables holding the dipole tilt,
!     PSI which comes from Geopack and TILT which comes from the variables file, which
!     ultimately also was computed using Geopack. So use either one.
CALL TS07D_JULY_2017 (IOPT,PARMOD,PSI,XGSM,YGSM,ZGSM, BXGSM,BYGSM,BZGSM)

PRINT *,'     Main field (nT):'
PRINT '(3A10)','HXGSM','HYGSM','HZGSM'
print '(3F10.3,/)',HXGSM,HYGSM,HZGSM
PRINT *,'   External field (nT):'
PRINT '(3A10)','BXGSM','BYGSM','BZGSM'
print '(3F10.3,/)', BXGSM,BYGSM,BZGSM

PI = 4*ATAN(1.)
PRINT '(A, F10.2)', 'Geodipole tilt (degrees):', PSI*180/PI

END program
