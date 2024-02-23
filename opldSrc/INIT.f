!***********************************************************************!

SUBROUTINE READ_INPUT

  USE SYSTEM

  !!PARAMETER INPUT FILE,MASTER CONTROL FILE
  OPEN (UNIT=10, FILE='par_init.dat')

  !!READING SYSTEM PARAMTER, SEE PAR_INIT.DAT FOR DETAILS
  READ (10, '( 11(I10/),5(F10.3/),2(I10/),8(F14.6/),(I10/),(ES15.5/),(L/),2(F10.6/),(F10.6) )') &
    REST, TJOB, TDAM, TPBC, TCRY, NTSTEP, NFREQ, NDEMAX, NAVMAX, NDIMER, TLAT, &
    RLAT, RBOX, PCUTR, DCUTR, DECUT, NDSMAX, NDRMAX, RDIMER, FNRMAX, FNRMIN, DSTEPF, &
    DSTEPM, TEMP, PKAE, DRATE, NFREKL, TIME, SIGRELON, SIGX, SIGY, SIGZ

  CLOSE (10)

END SUBROUTINE READ_INPUT
!***********************************************************************!
!!READING SYSTEM INPUT FROM VARIOUS FILES

SUBROUTINE SYS_INPUT

  USE SYSTEM
  USE LMP_DRIVER

  IMPLICIT NONE

  !!READING INTITIAL STRUCTURE

  NATOMC = int(lammps_get_natoms(lmp))
  if (rank .eq. 0) print *, NATOMC, " number of atoms"
  ALLOCATE (RA(3, NATOMC + NDEMAX))
  CALL RECV_LMP_POS

  !!CREATING REFERENCE STRUCTURE

  CALL CREATE_REF_LATTICE

  RETURN

END SUBROUTINE SYS_INPUT
!***********************************************************************!
!!SYSTEM INITITIALIZATION

SUBROUTINE SYS_INIT

  USE MPI
  USE SYSTEM
  USE DEFECT
  USE KMCVAR

  IMPLICIT NONE
  
  INTEGER, ALLOCATABLE :: SEED(:)
  INTEGER :: N

  !!EITHER RANDOM SEED OR FIXED RANDOM SEED

  CALL RANDOM_SEED()
  CALL RANDOM_SEED(SIZE=N)
  ALLOCATE (SEED(N))
  CALL RANDOM_SEED(GET=SEED)

  BSTEP = 0

  FDCONV = 1.0D-5

  !!HALF OF THE BOX SIZE, FOR PBC USE

  RBH1 = RBOX/2.0D0
  RBH2 = -RBOX/2.0D0

  !!POTENTIAL CUTOFF:PAIR TERM

  PCUTR = PCUTR/RLAT

  !!POTENTIAL CUTOFF:DENSITY TERM

  DCUTR = DCUTR/RLAT
  PCUTS = PCUTR*PCUTR
  DCUTS = DCUTR*DCUTR

  !!DEFECT CUTOFF

  DECUT = DECUT*DECUT
  NDVAC = 0
  NDINT = 0
  NNEIGH = 100
  UPNEIGH = .TRUE.

  !!IRRADIATION VARIABLES

  TIME_LAST = TIME
  DTIME = 0
  AVDPA = 0
  DPA = 0
  DPAC = 0
  CALC = .FALSE.

  !!INITIALIZING DEFECT VARIABLES

  ONDSTOT = 0
  NDSTOT = 0
  ONDTOT = 0
  NDTOT = 0

  !! KMC ENERGY BARRIER

  EMSKIP = -1*8.617E-5*TEMP*LOG(DRATE/1E12)

  !! SYSTEM PROPERTIES

  PRESS = 0.0D0
  EPSX = 1.0
  EPSY = 1.0
  EPSZ = 1.0

  !!VARIABLES ASSSOCIATED WITH MPI
  !! DIMER LIMIT VARIABLES BEING SET UP FOR J LOOP

  DIMERPR = NDIMER/NPROCS
  LDIMERI = DIMERPR*RANK + 1
  UDIMERI = DIMERPR*(RANK + 1)

  !!ALLOCATE DYNAMIC MATRIX
  !!FROM HERE,SEE DEFINE.F FOR THE MEANING OF VARIABLES

  ALLOCATE (FQ(NDIMER)); FQ = 0.0D0
  ALLOCATE (EM(NDIMER)); EM = 0.0D0
  ALLOCATE (OE(NDIMER)); OE = 0.0D0
  ALLOCATE (FF(3, NATOMC)); FF = 0.0D0
  ALLOCATE (EE(3, NATOMC)); EE = 0.0D0
  ALLOCATE (TV(NATOMR)); TV = 0
  ALLOCATE (TI(NATOMC)); TI = 0
  ALLOCATE (RNEIGH(NNEIGH, NATOMC)); RNEIGH = 0
  ALLOCATE (IV(3, NDEMAX)); IV = 0
  ALLOCATE (II(3, NDEMAX)); II = 0
  ALLOCATE (OC(10, NATOMR)); OC = 0
  ALLOCATE (OCC(NATOMR)); OCC = 0
  ALLOCATE (PV(3, NDEMAX)); PV = 0.0D0
  ALLOCATE (PI(3, NDEMAX)); PI = 0.0D0
  ALLOCATE (PX(7, NDEMAX)); PX = 0.0D0
  ALLOCATE (OX(7)); OX = 0.0D0
  ALLOCATE (ND(NDIMER)); ND = 0
  ALLOCATE (ON(NDIMER)); ON = 0
  ALLOCATE (SR(4, NAVMAX, NDIMER)); SR = 0.0D0
  ALLOCATE (OR(4, NAVMAX, NDIMER)); OR = 0.0D0

  IF (RANK .EQ. 0) THEN
    ALLOCATE (TVBUFF(NATOMR)); TVBUFF = 0
    ALLOCATE (TIBUFF(NATOMC)); TIBUFF = 0
    ALLOCATE (SRBUFF(4, NAVMAX, NDIMER)); SRBUFF = 0.0D0
    ALLOCATE (NDBUFF(NDIMER)); NDBUFF = 0
    ALLOCATE (EMBUFF(NDIMER)); EMBUFF = 0.0D0
  END IF

  !!VARIABLES ASSOCIATED WITH DIMER

  RDIMER = RDIMER/RLAT
  DSTEPF = DSTEPF/RLAT
  DSTEPM = DSTEPM/RLAT

  !!SETUP FOR PRE-SELECTION PROCESS

  SKIP = .FALSE.

  RETURN

END SUBROUTINE SYS_INIT

!***********************************************************************!
!BCC REFERENCE LATTICE CREATION ONLY
SUBROUTINE CREATE_REF_LATTICE

  USE SYSTEM, ONLY: RR, RBOX, NATOMR, TCRY, TLAT, RLAT, RANK

  IMPLICIT NONE

  LOGICAL :: EXISTS
  INTEGER :: I, THROWAWAY

  INQUIRE (FILE='pos_ref.dat', EXIST=EXISTS)
  IF (.NOT. EXISTS) THEN
    WRITE (*, *) "File 'pos_ref.dat' cannot be found."
    WRITE (*, *) "A reference file is needed for polycrystal."
    STOP
  END IF

  OPEN (UNIT=21, FILE='pos_ref.dat')
  READ (21, *) NATOMR
  READ (21, *)
  ALLOCATE (RR(3, NATOMR))
  READ (21, *) (THROWAWAY, RR(1, I), RR(2, I), RR(3, I), I=1, NATOMR)

  DO I = 1, NATOMR
    RR(1, I) = RR(1, I)/RLAT
    RR(2, I) = RR(2, I)/RLAT
    RR(3, I) = RR(3, I)/RLAT
  END DO

END SUBROUTINE CREATE_REF_LATTICE
