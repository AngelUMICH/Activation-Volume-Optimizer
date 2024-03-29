!***********************************************************************!
!***********************************************************************!
!***       SELF EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE          ***!
!***                        AUTHOR:XU,HAIXUAN                        ***!
!***             MATERIALS SCIENCE AND TECHNOLOGY DIVISION           ***!
!***                      OAK RIDGE NATIONAL LAB                     ***!
!***                          XUH1@ORNL.GOV                          ***!
!***                           VERSION 3.0                           ***!
!***                       ALL RIGHTS RESERVED                       ***!
!***********************************************************************!
!***********************************************************************!

!***********************************************************************!
!***********************************************************************!
!***********************************************************************!
MODULE CLASS_DIMER

  IMPLICIT NONE

  TYPE DIMER

    REAL*8 :: PI = 3.1415926535897931d0

!!NATOMT:TOTAL NUMBER OF ATOMS IN AN AV,LESS THAN NATOMM
!!NATOMS:TOTAL NUMBER OF ATOMS IN RC2
!!NATOMM:THE MAXIMUM NUMBER OF ATOMS IN AN AV, USER DEFINED

    INTEGER :: NATOMT, NATOMS, NATOMM

!!NATOM1:NUMBER OF ATOMS IN RC1
!!NATOM2:NUMBER OF ATOMS IN RC2
!!NATOM3:NUMBER OF ATOMS IN RC3

    INTEGER :: NATOM1, NATOM2, NATOM3

!!NROMAX:MAXIMUM NUMBER OF ROTATION DURING DIMER SEARCH
!!NROITR:THE ITERATION NUMBER OF ROTATION

    INTEGER :: NROMAX, NROITR

!!NTSMAX:MAXIMUM NUMBER OF TRANSLATION STEP IN DIMER SEARCH
!!NTSITR:THE ITERATION NUMBER OF TRANSLATION

    INTEGER :: NTSMAX, NTSITR

!!IX,IY,IZ: THE X,Y,Z POSITION OF THE AV, WHICH IS THE CENTER OF THE AV
!!RC1:CUTOFF RADIUS OF AV,THE SIZE OF THE AV
!!RC2:RC1+POTENTIAL CUTOFF, TO CALCULATE THE FORCE
!!RC3:RC1+POTENTIAL CUTOFF+ DENSITY CUTOFF,NEEDED FOR THE DENSITY

    REAL*8 :: IX, IY, IZ, RC1, RC2, RC3

!!DRATIO:A USER DEFINED RATIO,MAKE THE DIFFERENCE BETWEEN TWO DIMER IMAGE LARGER
!!DRATIO:ALSO AFFECT THE SADDLE POINT CONFIGURATION

    REAL*8 :: DRATIO

!!DISTDR:DISTANCE BETWEEN TWO DIMER IMAGE
!!ROFMIN:MINIMUM FORCE OF ROTATION CRITERIA
!!ROFMAX:MAXIMUM FORCE OF ROTATION CRITERIA

    REAL*8 :: DISTDR, ROFMIN, ROFMAX

!!FTOTAL:TOTAL FORCE OF ATOMS IN AV
!!FCONVG:FORCE CONVERGENCE CRITERIA
!!FNUNIT:SEE LINE 1300

    REAL*8 :: FTOTAL, FCONVG, FNUNIT

!!ENINIT:INITIAL ENERGY
!!ENCALC:CURRENT ENERGY
!!ENLAST:LAST STEP ENERGY

    REAL*8 :: ENINIT, ENCALC, ENLAST

!!EMDIFF:MIGRATION ENERGY BARRIER,DIFFERENCE BETWEEN INTITAL AND SADDLE
!!DISPS1:SUM OF THE DISTANCE FOR ONE DIMER IMAGE FROM THE INTITIAL
!!DISPS2:SUM OF THE DISTANCE FOR THE OTHER DIMER IMAGE FROM THE INITIAL

    REAL*8 :: EMDIFF, DISPS1, DISPS2

!!RTHETA:THETA, SEE REFERENCE FOR ALOGITHM

    REAL*8 :: RTHETA

!!ROCGA1:ROTATION_CG STEP_TEMP VARIABLE_A1
!!ROCGA2:ROTATION_CG STEP_TEMP VARAIBLE_A2
!!ROGAMN:ROTATION_VARIABLE_GAMN

    REAL*8 :: ROCGA1, ROCGA2, ROGAMN

!!ROTFN1:ROTATION_FORCE_FN1
!!ROTFN2:ROTATION_FORCE_FN2
!!ROCURV:ROTATION_CURVATURE

    REAL*8 :: ROTFN1, ROTFN2, ROCURV

!!ROFNRS:ROTATION FORCE_FN_STORAGE

    REAL*8 :: ROFNRS

!!TSTEPC:TRANSLATION STEP SIZE
!!TSTEPM:MAXIMUM TRANSLATION STEP SIZE
!!TSTEPF:TRIAL STEP SIZE IN THE TRANSLATION

    REAL*8 :: TSTEPC, TSTEPM, TSTEPF

!!TRTFP1:TRANSLATION_TEMP VARIABLE_FP1
!!TRTFP2:TRANSLATION_TEMP VARIABLE_FP2

    REAL*8 :: TRTFP1, TRTFP2, TRCURV

!!TRGAMN:TRANSLATION_GAMN
!!TRCGA1:TRANSLATION_CG_A1
!!TRCGA2:TRANSLATION_CG_A2

    REAL*8 :: TRGAMN, TRCGA1, TRCGA2

!!DI_MAX:DIMER REACH MAXIMUM STEP TAG

    LOGICAL :: DI_MAX

!!DI_CON:DIMER CONVERGENCE TAG
!!TR_FWD:TRANSLATION FORWARD STEP
!!RO_FWD:ROTATION FORWARD STEP

    LOGICAL :: DI_CON, TR_FWD, RO_FWD

!!RO_NEW:ROTATION NEW
!!RO_OPT:ROTATION OPT
!!RO_CGI:ROTATION CG

    LOGICAL :: RO_NEW, RO_OPT, RO_CGI

!!FO_CON:FORCE CONVERGENCE
!!ST_CON:STEP CONVERGENCE
!!EN_CON:ENERGY CONVERGENCE

    LOGICAL :: FO_CON, ST_CON, EN_CON

!!ID:ATOM IDENTIFICATION,3XN MATRIX, 1:ATOM INDEX,2:REGION INDEX,3:NOT USED

    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ID

!!SI:STORAGE MATRIX OF ID

    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: SI

!!XA:INITIAL POSITION OF ALL ATOMS RELAVANT

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: XA

!!XS:POSITION STORAGE MATRIX

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: XS

!!X0:MIDDLE POINT OF THE DIMER IMAGE

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: X0

!!X1:POSITION OF ONE DIMER IMAGE

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: X1

!!X2:POSITION OF THE OTHER DIMER IMAGE

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: X2

!!XC:POSITION AT THE CURRENT MOMENT

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: XC

!!FO:FORCE AT THE MIDDLE POINT

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: F0

!!F1:FORCE OF ONE DIMER IMAGE

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: F1

!!F2:FORCE OF THE OTHER DIMER IMAGE

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: F2

!!FC:FORCE AT THE CURRENT MOMENT

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: FC

!!FE:EFFECTIVE FORCE

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: FE

!!FN:FORCE ALONG N

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: FN

!!FL:FORCE OF THE LAST STEP

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: FL

!!MATRIX ASSOCIATED WITH CG

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: GN

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: GU

!!TEMPRARY VARIABLES

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: TF

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: TG

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: TU

!!VN:UNIT VECTOR

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: VN

!!SX:POSTION AT THE SADDLE POINTS,JUSTED OR UNJUSTED

    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: SX

  END TYPE DIMER

!***********************************************************************!

CONTAINS

!***********************************************************************!
!!INITIALIZATION OF THE DIMER

  SUBROUTINE DIMER_INIT(THIS, IX, IY, IZ, RC)

    USE LMP_DRIVER

    USE SYSTEM

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    REAL*8, INTENT(IN) :: IX, IY, IZ, RC

    REAL*8 RX, RY, RZ, RS, PBC

    INTEGER I, NTEMP

    CHARACTER(LEN=200) S_NATOM1

!!PASSING THE POSITION

    THIS%IX = IX

    THIS%IY = IY

    THIS%IZ = IZ

!!POSITION THE CUTOFF

    THIS%RC1 = RC**2

    THIS%RC2 = (RC + PCUTR)**2

    THIS%RC3 = (RC + PCUTR + DCUTR)**2

!!MAMIMUM TRANSLATION STEP:USER DEFINE

    THIS%NTSMAX = NDSMAX

!!RATIO:USED DEFINE (1.0DO MEANS NO EFFECTS)

    THIS%DRATIO = 1.0D0

!!CONVERENCE

    THIS%FCONVG = FDCONV

!!MAXIMUM ROTATION STEPS

    THIS%NROMAX = NDRMAX

!!MAXIMUM NUMBER OF ATOMS IN AV

    THIS%NATOMM = NAVMAX

!!DIMER DISTANCE

    THIS%DISTDR = RDIMER

    THIS%TSTEPM = DSTEPM

    THIS%TSTEPF = DSTEPF

    THIS%ROFMIN = FNRMIN

    THIS%ROFMAX = FNRMAX

    THIS%DISTDR = RDIMER

    THIS%NATOM1 = 0

    THIS%NATOM2 = 0

    THIS%NATOM3 = 0

!!IDENTIFY ATOMS IN THE AV BASED ON THEIR POSITIONS

    DO I = 1, NATOMC

      RX = PBC(IX - RA(1, I), 1)

      RY = PBC(IY - RA(2, I), 2)

      RZ = PBC(IZ - RA(3, I), 3)

      RS = RX*RX + RY*RY + RZ*RZ

      IF (RS .LE. THIS%RC1) THEN

        THIS%NATOM1 = THIS%NATOM1 + 1

      ELSEIF ((RS .GT. THIS%RC1) .AND. (RS .LE. THIS%RC2)) THEN

        THIS%NATOM2 = THIS%NATOM2 + 1

      ELSEIF ((RS .GT. THIS%RC2) .AND. (RS .LE. THIS%RC3)) THEN

        THIS%NATOM3 = THIS%NATOM3 + 1

      END IF

    END DO

!      WRITE(*,*) THIS%NATOM1,THIS%NATOM2,THIS%NATOM3

!!CALCULATE THE NUMBER OF ATOMS IN EACH REGION

    THIS%NATOMT = THIS%NATOM1 + THIS%NATOM2 + THIS%NATOM3

    THIS%NATOMS = THIS%NATOM1 + THIS%NATOM2

    IF (THIS%NATOMT .GT. THIS%NATOMM) WRITE (*, *) 'INCREASE NAVMAX'

!!$OMP CRITICAL

    ALLOCATE (THIS%ID(3, THIS%NATOMT)); THIS%ID = 0

    ALLOCATE (THIS%SI(3, THIS%NATOMT)); THIS%SI = 0

    ALLOCATE (THIS%XA(3, THIS%NATOMT)); THIS%XA = 0.0D0

    ALLOCATE (THIS%XS(3, THIS%NATOMT)); THIS%XS = 0.0D0

    ALLOCATE (THIS%X0(3, THIS%NATOM1)); THIS%X0 = 0.0D0

    ALLOCATE (THIS%X1(3, THIS%NATOM1)); THIS%X1 = 0.0D0

    ALLOCATE (THIS%X2(3, THIS%NATOM1)); THIS%X2 = 0.0D0

    ALLOCATE (THIS%XC(3, THIS%NATOM1)); THIS%XC = 0.0D0

    ALLOCATE (THIS%F0(3, THIS%NATOM1)); THIS%F0 = 0.0D0

    ALLOCATE (THIS%F1(3, THIS%NATOM1)); THIS%F1 = 0.0D0

    ALLOCATE (THIS%F2(3, THIS%NATOM1)); THIS%F2 = 0.0D0

    ALLOCATE (THIS%FC(3, THIS%NATOM1)); THIS%FC = 0.0D0

    ALLOCATE (THIS%FE(3, THIS%NATOM1)); THIS%FE = 0.0D0

    ALLOCATE (THIS%FN(3, THIS%NATOM1)); THIS%FN = 0.0D0

    ALLOCATE (THIS%FL(3, THIS%NATOM1)); THIS%FL = 0.0D0

    ALLOCATE (THIS%GN(3, THIS%NATOM1)); THIS%GN = 0.0D0

    ALLOCATE (THIS%GU(3, THIS%NATOM1)); THIS%GU = 0.0D0

    ALLOCATE (THIS%TF(3, THIS%NATOM1)); THIS%TF = 0.0D0

    ALLOCATE (THIS%TG(3, THIS%NATOM1)); THIS%TG = 0.0D0

    ALLOCATE (THIS%TU(3, THIS%NATOM1)); THIS%TU = 0.0D0

    ALLOCATE (THIS%VN(3, THIS%NATOM1)); THIS%VN = 0.0D0

    ALLOCATE (THIS%SX(3, THIS%NATOM1)); THIS%SX = 0.0D0

!!$OMP END CRITICAL

    NTEMP = 0

!!SORTING ATOMS BASED ON THEIR POSITIONS AND GIVE THEM INDEX

    DO I = 1, NATOMC

      RX = PBC(IX - RA(1, I), 1)

      RY = PBC(IY - RA(2, I), 2)

      RZ = PBC(IZ - RA(3, I), 3)

      RS = RX*RX + RY*RY + RZ*RZ

      IF (RS .LE. THIS%RC1) THEN

        NTEMP = NTEMP + 1

        THIS%ID(1, NTEMP) = I

        THIS%ID(2, NTEMP) = 1

        THIS%XS(1, NTEMP) = RA(1, I)

        THIS%XS(2, NTEMP) = RA(2, I)

        THIS%XS(3, NTEMP) = RA(3, I)

      END IF

    END DO

    DO I = 1, NATOMC

      RX = PBC(IX - RA(1, I), 1)

      RY = PBC(IY - RA(2, I), 2)

      RZ = PBC(IZ - RA(3, I), 3)

      RS = RX*RX + RY*RY + RZ*RZ

      IF ((RS .GT. THIS%RC1) .AND. (RS .LE. THIS%RC2)) THEN

        NTEMP = NTEMP + 1

        THIS%ID(1, NTEMP) = I

        THIS%ID(2, NTEMP) = 2

        THIS%XS(1, NTEMP) = RA(1, I)

        THIS%XS(2, NTEMP) = RA(2, I)

        THIS%XS(3, NTEMP) = RA(3, I)

      END IF

    END DO

    DO I = 1, NATOMC

      RX = PBC(IX - RA(1, I), 1)

      RY = PBC(IY - RA(2, I), 2)

      RZ = PBC(IZ - RA(3, I), 3)

      RS = RX*RX + RY*RY + RZ*RZ

      IF ((RS .GT. THIS%RC2) .AND. (RS .LE. THIS%RC3)) THEN

        NTEMP = NTEMP + 1

        THIS%ID(1, NTEMP) = I

        THIS%ID(2, NTEMP) = 3

        THIS%XS(1, NTEMP) = RA(1, I)

        THIS%XS(2, NTEMP) = RA(2, I)

        THIS%XS(3, NTEMP) = RA(3, I)

      END IF

    END DO

    THIS%XA = THIS%XS

    DO I = 1, THIS%NATOM1

      THIS%XC(1, I) = THIS%XA(1, I)

      THIS%XC(2, I) = THIS%XA(2, I)

      THIS%XC(3, I) = THIS%XA(3, I)

    END DO

!!INITIALIZE LOGICAL TAGS

    THIS%DI_CON = .FALSE.

    THIS%DI_MAX = .FALSE.

    THIS%FO_CON = .FALSE.

    THIS%ST_CON = .FALSE.

    THIS%EN_CON = .FALSE.

    THIS%ROCURV = 0.0D0

    THIS%NTSITR = 1

    THIS%RO_NEW = .TRUE.

    THIS%RO_FWD = .TRUE.

    THIS%TR_FWD = .TRUE.

!$OMP CRITICAL

    CALL RANDOM_NUMBER(THIS%VN)

!$OMP END CRITICAL

    THIS%VN = THIS%VN - 0.5D0

!!CREATE A RANDOM DISPLACEMENT

    CALL DIMER_UNIT(THIS%VN, THIS%NATOM1)

!SEND DEFECT POSITIONS TO LMPDIMER

    CALL CREATE_ATOMS(THIS)

!COMPUTING PER ATOM POTENTIAL ENERGY

    CALL LAMMPS_COMMAND(lmpDimer, 'compute peratom all pe/atom')

!SUM REDUCE ONLY NATOM1 POTENTIAL ENERGY

    WRITE (S_NATOM1, '(1I0)') THIS%NATOM1

    CALL LAMMPS_COMMAND(lmpDimer, 'group saddle id 1:'//TRIM(S_NATOM1))

    CALL LAMMPS_COMMAND(lmpDimer, 'compute pesum saddle reduce sum c_peratom')

    RETURN

  END SUBROUTINE DIMER_INIT
!***********************************************************************!
!!CLEAN UP THE DIMER

  SUBROUTINE DIMER_FINISH(THIS)

    USE LMP_DRIVER

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    !REMOVE THE SYSTEM FROM THE DIMER LAMMPS INSTANCE

    CALL LAMMPS_COMMAND(lmpDimer, 'uncompute pesum')

    CALL LAMMPS_COMMAND(lmpDimer, 'uncompute peratom')

    CALL LAMMPS_COMMAND(lmpDimer, 'group saddle delete')

    CALL LAMMPS_COMMAND(lmpDimer, 'delete_atoms group all')

!$OMP CRITICAL

    DEALLOCATE (THIS%ID)

    DEALLOCATE (THIS%XA)

    DEALLOCATE (THIS%XS)

    DEALLOCATE (THIS%X0)

    DEALLOCATE (THIS%X1)

    DEALLOCATE (THIS%X2)

    DEALLOCATE (THIS%XC)

    DEALLOCATE (THIS%F0)

    DEALLOCATE (THIS%F1)

    DEALLOCATE (THIS%F2)

    DEALLOCATE (THIS%FC)

    DEALLOCATE (THIS%FE)

    DEALLOCATE (THIS%FN)

    DEALLOCATE (THIS%FL)

    DEALLOCATE (THIS%GN)

    DEALLOCATE (THIS%GU)

    DEALLOCATE (THIS%TF)

    DEALLOCATE (THIS%TG)

    DEALLOCATE (THIS%TU)

    DEALLOCATE (THIS%VN)

    DEALLOCATE (THIS%SI)

    DEALLOCATE (THIS%SX)

!$OMP END CRITICAL

    RETURN

  END SUBROUTINE DIMER_FINISH
!***********************************************************************!
!!MAJOR SUBROUTINE OF THE DIMER
!!PLEASE SEE THE FOLLOWING REFERNCES TO UNDERSTAND THE ALOGITHM
!!J.CHEM.PHYS,VOL.111,NO.15,7010,(1999) BY G.HENKELMAN ET AL.
!!J.CHEM.PHYS,VOL.123.224101,(2005) BY A.HEYDEN ET AL.

  SUBROUTINE DIMER_SEARCH(THIS)

    IMPLICIT NONE

    INTEGER I

    REAL*8 RX, RY, RZ, RS

    TYPE(DIMER), INTENT(INOUT) :: THIS

!!CALCULATE INITIAL ENERGY

    CALL DIMER_FORCE_LMP(THIS)

    THIS%ENINIT = THIS%ENCALC

!!CHECK CONVERGENCE

    DO WHILE (THIS%DI_CON .EQV. .FALSE.)

      CALL DIMER_FORCE_LMP(THIS)

!!AVOID THE FIRST SEVERAL STEPS

      IF ((THIS%TR_FWD) .AND. (THIS%NROITR .EQ. 0)) THEN

        THIS%EN_CON = (ABS(THIS%ENCALC - THIS%ENLAST) .LT. 1.0D-5)

        THIS%ENLAST = THIS%ENCALC

      END IF

      CALL DIMER_ROTATE(THIS)

      THIS%FTOTAL = 0.0D0

      DO I = 1, THIS%NATOM1

        THIS%FTOTAL = THIS%FTOTAL + THIS%FC(1, I)*THIS%FC(1, I)

        THIS%FTOTAL = THIS%FTOTAL + THIS%FC(2, I)*THIS%FC(2, I)

        THIS%FTOTAL = THIS%FTOTAL + THIS%FC(3, I)*THIS%FC(3, I)

      END DO

!!CALCULATE TOTAL FORCE OF THE SYSTEM

      THIS%FTOTAL = DSQRT(THIS%FTOTAL)

!!CONVERGENCE CRITERIA

      THIS%FO_CON = (THIS%FTOTAL .LT. THIS%FCONVG)

      THIS%ST_CON = (THIS%NTSITR .GT. 5)

      THIS%DI_CON = ((THIS%FO_CON) .AND. THIS%ST_CON)

!!IF CONVERGED

      IF (THIS%DI_CON) THEN

        THIS%DI_CON = .TRUE.

!!OTHERWISE

      ELSE

        IF (THIS%RO_OPT) CALL DIMER_TRANSLATE(THIS)

      END IF

!!IF MAXIMUM STEPS ARE REACHED

      IF (THIS%NTSITR .GE. THIS%NTSMAX) THEN

        THIS%DI_MAX = .TRUE.

        THIS%DI_CON = .TRUE.

        EXIT

      END IF

    END DO

!!PASSING POSITIONS AND VARIABLES BACK

    THIS%X1 = THIS%XC

    THIS%X2 = 2.0D0*THIS%X0 - THIS%X1

    THIS%DISPS1 = 0.0D0

    THIS%DISPS2 = 0.0D0

!!CALCULATE WHICH IMAGE IS FAR FROM THE INITIAL CONFIGURATION

    DO I = 1, THIS%NATOM1

      RX = THIS%X1(1, I) - THIS%XS(1, I)

      RY = THIS%X1(2, I) - THIS%XS(2, I)

      RZ = THIS%X1(3, I) - THIS%XS(3, I)

      RS = RX*RX + RY*RY + RZ*RZ

      THIS%DISPS1 = THIS%DISPS1 + RS

      RX = THIS%X2(1, I) - THIS%XS(1, I)

      RY = THIS%X2(2, I) - THIS%XS(2, I)

      RZ = THIS%X2(3, I) - THIS%XS(3, I)

      RS = RX*RX + RY*RY + RZ*RZ

      THIS%DISPS2 = THIS%DISPS2 + RS

    END DO

!!STORE THE POSITION OF THE SADDLE POINT CONFIGURATION
!!IF ANY RATIO IS NEEDED TO MAKE THE DIFFERENCE LARGER

    IF (THIS%DISPS1 .GE. THIS%DISPS2) THEN

      THIS%SX = (THIS%X1 - THIS%X0)*THIS%DRATIO + THIS%X0

    ELSE IF (THIS%DISPS2 .GT. THIS%DISPS1) THEN

      THIS%SX = (THIS%X2 - THIS%X0)*THIS%DRATIO + THIS%X0

    ELSE

!!CAUTION IS NEEDED IF THE TWO IMAGES ARE THE SAME

      WRITE (*, *) 'DIMER IMAGES ARE SAME'

    END IF

    THIS%SI = THIS%ID

    THIS%XC = THIS%X0

!!CALCULATE THE FORCE AND ENERGY OF THE SADDLE POINT

    CALL DIMER_FORCE_LMP(THIS)

!!CALCULATE ENERGY BARRIER

    THIS%EMDIFF = THIS%ENCALC - THIS%ENINIT

!!IF MAXIMUM STEP IS REACHED, DO NOT USE THIS CONFIGURATION

    IF (THIS%DI_MAX) THIS%EMDIFF = 100.0D0

    RETURN

  END SUBROUTINE DIMER_SEARCH
!***********************************************************************!
!!PLEASE THE ABOVE REFERENCE TO UNDERSTAND THE ALOGRITHM

  SUBROUTINE DIMER_ROTATE(THIS)

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    INTEGER I

    IF (THIS%RO_OPT) THEN

      THIS%F0 = THIS%FC

      CALL DIMER_PROJECTION(THIS)

      THIS%FC = THIS%FE

      RETURN

    END IF

    IF (THIS%RO_NEW) THEN

      THIS%RO_NEW = .FALSE.

      THIS%RO_FWD = .TRUE.

      THIS%NROITR = 1

      THIS%X0 = THIS%XC

      THIS%F0 = THIS%FC

      CALL DIMER_PROJECTION(THIS)

      THIS%XC = THIS%X0 + THIS%VN*THIS%DISTDR

    ELSE

      CALL DIMER_UPDATE(THIS)

      IF (THIS%RO_FWD) THEN

        THIS%ROCURV = (SUM(THIS%F2*THIS%VN) - SUM(THIS%F1*THIS%VN))/ &
                      (2.0D0*THIS%DISTDR)

        CALL DIMER_PROJECTION(THIS)

      END IF

      IF (THIS%FNUNIT .LT. THIS%ROFMIN) THEN

        THIS%RO_OPT = .TRUE.

      ELSE

        CALL DIMER_ROTATESTEP(THIS)

        IF (THIS%RO_FWD) CALL DIMER_PROJECTION(THIS)

        IF (THIS%NROITR .GT. THIS%NROMAX) THEN

          THIS%RO_OPT = .TRUE.

        ELSE IF (THIS%RO_FWD .AND. (THIS%ROFNRS .LT. THIS%ROFMAX)) THEN

          THIS%RO_OPT = .TRUE.

        END IF

      END IF

    END IF

    IF (THIS%RO_OPT) THEN

      THIS%XC = THIS%X0

      THIS%RO_NEW = .TRUE.

      THIS%RO_FWD = .TRUE.

      THIS%NROITR = 0

      THIS%NTSITR = THIS%NTSITR + 1

    END IF

    DO I = 1, THIS%NATOM1

      THIS%XA(1, I) = THIS%XC(1, I)

      THIS%XA(2, I) = THIS%XC(2, I)

      THIS%XA(3, I) = THIS%XC(3, I)

    END DO

    THIS%FC = THIS%FE

    RETURN

  END SUBROUTINE DIMER_ROTATE
!***********************************************************************!
!!ROTATION STEP

  SUBROUTINE DIMER_ROTATESTEP(THIS)

    USE SYSTEM, ONLY: ISTEP, NFREQ

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    IF (THIS%RO_FWD) THEN

      THIS%RO_FWD = .FALSE.

      IF (THIS%RO_CGI) THEN

        THIS%RO_CGI = .FALSE.

        THIS%FL = THIS%FN

        THIS%GN = THIS%FN

      END IF

      THIS%ROCGA1 = ABS(SUM(THIS%FN*THIS%FL))

      THIS%ROCGA2 = SUM(THIS%FL*THIS%FL)

      IF ((THIS%ROCGA1 .LE. (0.5D0*THIS%ROCGA2)) .AND. &
          (THIS%ROCGA2 .NE. 0.0D0)) THEN

        THIS%ROGAMN = SUM(THIS%FN*(THIS%FN - THIS%FL))/THIS%ROCGA2

      ELSE

        THIS%ROGAMN = 0.0D0

      END IF

      IF (ABS(THIS%ROGAMN) .GT. 1.0D0) THEN

        IF (MOD(ISTEP, NFREQ) .EQ. 0) THEN
          WRITE (*, *) 'ROGAMN.GT.1.00'
        END IF

        THIS%ROGAMN = 0.0D0

      END IF

      THIS%GN = THIS%FN + THIS%GN*THIS%ROGAMN

      THIS%GU = THIS%GN/DSQRT(SUM(THIS%GN*THIS%GN))

      THIS%ROTFN1 = SUM(THIS%FN*THIS%GU)

      CALL DIMER_TRANSFORM(THIS%VN, THIS%GU, THIS%PI/4.0D0, THIS%NATOM1)

      THIS%XC = THIS%X0 + THIS%VN*THIS%DISTDR

      THIS%ROFNRS = THIS%FNUNIT

    ELSE

      THIS%RO_FWD = .TRUE.

      THIS%NROITR = THIS%NROITR + 1

      THIS%ROTFN2 = SUM(THIS%FN*THIS%GU)

      IF (THIS%ROTFN2 .NE. 0.0D0) THEN

        THIS%RTHETA = ATAN(THIS%ROTFN1/THIS%ROTFN2)/(-2.0D0)

      ELSE

        THIS%RTHETA = THIS%PI/(-2.0D0)

      END IF

      IF (THIS%ROTFN2 .GT. 0.0D0) THIS%RTHETA = THIS%RTHETA + THIS%PI/2.0D0

      THIS%RTHETA = THIS%RTHETA - THIS%PI/4.0D0

      CALL DIMER_TRANSFORM(THIS%VN, THIS%GU, THIS%RTHETA, THIS%NATOM1)

      CALL DIMER_UNIT(THIS%VN, THIS%NATOM1)

      THIS%XC = THIS%X0 + THIS%VN*THIS%DISTDR

    END IF

    RETURN

  END SUBROUTINE DIMER_ROTATESTEP
!***********************************************************************!
!!TRANSLATION STEP

  SUBROUTINE DIMER_TRANSLATE(THIS)

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    INTEGER I

!!FORWARD TRIP STEP

    IF (THIS%TR_FWD) THEN

      THIS%TR_FWD = .FALSE.

      THIS%RO_OPT = .TRUE.

      THIS%TRCGA1 = ABS(SUM(THIS%FC*THIS%TF))

      THIS%TRCGA2 = SUM(THIS%TF*THIS%TF)

      IF ((THIS%TRCGA1 .LE. (0.5D0*THIS%TRCGA2)) .AND. &
          (THIS%TRCGA2 .NE. 0.0D0)) THEN

        THIS%TRGAMN = SUM(THIS%FC*(THIS%FC - THIS%TF))/THIS%TRCGA2

      ELSE

        THIS%TRGAMN = 0.0D0

      END IF

      THIS%TG = THIS%FC + THIS%TG*THIS%TRGAMN

      THIS%TU = THIS%TG/DSQRT(SUM(THIS%TG*THIS%TG))

      THIS%TF = THIS%FC

      THIS%XC = THIS%XC + THIS%TU*THIS%TSTEPF

    ELSE

!!REAL STEP

      THIS%TR_FWD = .TRUE.

      THIS%RO_OPT = .FALSE.

      THIS%TRTFP1 = SUM(THIS%TF*THIS%TU)

      THIS%TRTFP2 = SUM(THIS%FC*THIS%TU)

      THIS%TRCURV = (THIS%TRTFP1 - THIS%TRTFP2)/THIS%TSTEPF

      IF (THIS%TRCURV .LT. 0.0D0) THEN

        THIS%TSTEPC = THIS%TSTEPM

      ELSE

        THIS%TSTEPC = 0.5D0*(THIS%TRTFP1 + THIS%TRTFP2)/THIS%TRCURV

        IF (ABS(THIS%TSTEPC) .GT. THIS%TSTEPM) THEN

          THIS%TSTEPC = SIGN(THIS%TSTEPM, THIS%TSTEPC) - &
                        SIGN(THIS%TSTEPF, THIS%TSTEPC)

        ELSE

          THIS%TSTEPC = THIS%TSTEPC - 0.5D0*THIS%TSTEPF

        END IF

      END IF

      THIS%XC = THIS%XC + THIS%TU*THIS%TSTEPC

    END IF

    DO I = 1, THIS%NATOM1

      THIS%XA(1, I) = THIS%XC(1, I)

      THIS%XA(2, I) = THIS%XC(2, I)

      THIS%XA(3, I) = THIS%XC(3, I)

    END DO

    RETURN

  END SUBROUTINE DIMER_TRANSLATE
!***********************************************************************!
!ENERGY AND FORCES ARE NEEDED FOR ATOMS UP TO NATOM1
  SUBROUTINE DIMER_FORCE_LMP(THIS)

    USE LMP_DRIVER

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS
    REAL(C_double), POINTER :: ENERGY => NULL()
    REAL*8, DIMENSION(:), ALLOCATABLE :: F
    INTEGER :: I

    THIS%FC = 0.0D0

    THIS%ENCALC = 0.0D0

    !UPDATE LAMMPS INSTANCE POSITIONS
    CALL SEND_POS(THIS)

    !UPDATE ALL VALUES
    CALL LAMMPS_COMMAND(lmpDimer, 'run 0')

    !CALCULATE ENERGY
    CALL LAMMPS_EXTRACT_COMPUTE(ENERGY, lmpDimer, 'pesum', &
                                LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)

    !CALCULATE FORCES
    CALL LAMMPS_GATHER_ATOMS(lmpDimer, 'f', 3, F)

    DO I = 1, THIS%NATOM1

      THIS%FC(1, I) = F(3*(I - 1) + 1)
      THIS%FC(2, I) = F(3*(I - 1) + 2)
      THIS%FC(3, I) = F(3*(I - 1) + 3)

    END DO

    THIS%ENCALC = ENERGY

  END SUBROUTINE DIMER_FORCE_LMP

!***********************************************************************!
!!UPDATE THE FORCE OF THE DIMER

  SUBROUTINE DIMER_UPDATE(THIS)

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    THIS%F1 = THIS%FC

    THIS%F2 = THIS%F0*2.0D0 - THIS%F1

    THIS%FN = (THIS%F1 - THIS%F2 - THIS%VN*(SUM(THIS%F1*THIS%VN) - &
                                            SUM(THIS%F2*THIS%VN)))/(2.0D0*THIS%DISTDR)

    THIS%FNUNIT = DSQRT(SUM(THIS%FN*THIS%FN))

    RETURN

  END SUBROUTINE DIMER_UPDATE
!***********************************************************************!
!!CALCULATE THE PROJECTION FORCE OF DIMER BASED ON THE CURVATURE

  SUBROUTINE DIMER_PROJECTION(THIS)

    IMPLICIT NONE

    TYPE(DIMER), INTENT(INOUT) :: THIS

    IF (THIS%ROCURV .LT. 0.0D0) THEN

      THIS%FE = THIS%F0 - THIS%VN*SUM(THIS%F0*THIS%VN)*2.0D0

    ELSE

      THIS%FE = -THIS%VN*SUM(THIS%F0*THIS%VN)

    END IF

    RETURN

  END SUBROUTINE DIMER_PROJECTION
!***********************************************************************!
!!TRANSFORM A VECTOR INTO A UNIT VECTOR

  SUBROUTINE DIMER_UNIT(V1, N)

    INTEGER, INTENT(IN) :: N

    REAL*8, INTENT(INOUT) :: V1(3, N)

    V1 = V1*(1.0D0/DSQRT(SUM(V1*V1)))

    RETURN

  END SUBROUTINE DIMER_UNIT
!***********************************************************************!
!!VECTOR TRANSFROMATION OPERATION,SIMPLE MATH

  SUBROUTINE DIMER_TRANSFORM(V1, V2, THETA, N)

    INTEGER, INTENT(IN) :: N

    REAL*8, INTENT(IN) :: THETA

    REAL*8, INTENT(INOUT) :: V1(3, N), V2(3, N)

    REAL*8 COST, SINT, V3(3, N)

    COST = COS(THETA)

    SINT = SIN(THETA)

    V3 = V1

    V1 = V1*COST + V2*SINT

    V2 = V2*COST - V3*SINT

    RETURN

  END SUBROUTINE DIMER_TRANSFORM
!***********************************************************************!
!USED TO INITIATE ATOMS IN THE DIMER LAMMPS INSTANCE
  SUBROUTINE CREATE_ATOMS(THIS)

    USE LMP_DRIVER
    USE SYSTEM, ONLY: RLAT

    IMPLICIT NONE

    TYPE(DIMER), INTENT(IN) :: THIS

    REAL*8, DIMENSION(THIS%NATOMT, 3) :: X
    INTEGER, DIMENSION(THIS%NATOMT) :: TYPE, ID
    INTEGER :: I

    !ALL ATOMS ARE OF TYPE 1
    TYPE = 1

    DO I = 1, THIS%NATOMT
      X(I, :) = THIS%XA(:, I)*RLAT
      ID(I) = I
    END DO

    CALL LAMMPS_CREATE_ATOMS(lmpDimer, ID=ID, TYPE=TYPE, X=X, shrinkexceed=.FALSE.)

  END SUBROUTINE CREATE_ATOMS
!***********************************************************************!
  SUBROUTINE SEND_POS(THIS)

    USE LMP_DRIVER
    USE SYSTEM, ONLY: RLAT

    IMPLICIT NONE

    TYPE(DIMER), INTENT(IN) :: THIS
    INTEGER :: I

    REAL*8, DIMENSION(:) :: POS(3*THIS%NATOMT)

    DO I = 1, THIS%NATOMT
      POS(3*(I - 1) + 1:3*(I - 1) + 3) = THIS%XA(:, I)*RLAT
    END DO

    CALL lammps_scatter_atoms(lmpDimer, 'x', POS)

  END SUBROUTINE SEND_POS
!***********************************************************************!
END MODULE CLASS_DIMER

!***********************************************************************!
!***********************************************************************!
!***********************************************************************!

!***********************************************************************!
!***********************************************************************!
!***       SELF EVOLVING ATOMISTIC KINETIC MONTE CARLO CODE          ***!
!***                        AUTHOR:XU,HAIXUAN                        ***!
!***             MATERIALS SCIENCE AND TECHNOLOGY DIVISION           ***!
!***                      OAK RIDGE NATIONAL LAB                     ***!
!***                          XUH1@ORNL.GOV                          ***!
!***                           VERSION 3.0                           ***!
!***                       ALL RIGHTS RESERVED                       ***!
!***********************************************************************!
!***********************************************************************!
