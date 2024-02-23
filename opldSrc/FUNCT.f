!***********************************************************************!
!!APPLY THE PERIODIC BOUNDARY CONDITION TO THE SYSTEM
REAL*8 FUNCTION PBC(X, DIM)

  USE SYSTEM, ONLY: RBOX, RBH1, RBH2, EPSX, EPSY, EPSZ

  IMPLICIT NONE

  INTEGER DIM
  REAL*8 X, EPS

  IF (DIM .EQ. 1) THEN
    EPS = EPSX
  ELSEIF (DIM .EQ. 2) THEN
    EPS = EPSY
  ELSEIF (DIM .EQ. 3) THEN
    EPS = EPSZ
  END IF

  IF ((X .LE. RBH1*EPS) .AND. (X .GE. RBH2*EPS)) THEN
    PBC = X
  ELSEIF (X .LT. RBH2*EPS) THEN
    PBC = X + RBOX*EPS
  ELSE
    PBC = X - RBOX*EPS
  END IF

  RETURN

END FUNCTION PBC
!***********************************************************************!
SUBROUTINE FIND_DEFECT

  USE SYSTEM
  USE DEFECT

  IMPLICIT NONE

  INTEGER I

  NDINT = 0
  NDVAC = 0
  TI = 0
  TV = 0

  CALL POLY_CRYSTAL_WS

  !!TALLY NUMBER OF VACANCIES AND INTERSTITIALS

  DO I = 1, NATOMC
    IF (TI(I) .NE. 1) NDINT = NDINT + 1
  END DO

  DO I = 1, NATOMR
    IF (TV(I) .NE. 1) NDVAC = NDVAC + 1
  END DO

  RETURN
END SUBROUTINE FIND_DEFECT
!***********************************************************************!
SUBROUTINE POLY_CRYSTAL_WS

  USE OMP_LIB
  USE SYSTEM
  USE DEFECT
  USE MPI

  IMPLICIT NONE

  REAL*8 RX, RY, RZ, RS, RMIN, PBC
  INTEGER I, J, MIN, N, ATOM
  INTEGER IDHI, IDLO, RID
  INTEGER OCCBUFF(SIZE(OCC))

  IDLO = NATOMC/NPROCS*(RANK) + 1
  IDHI = NATOMC/NPROCS*(RANK + 1)
  IF (RANK .eq. NPROCS - 1) IDHI = IDHI + MOD(NATOMC, NPROCS)

  CALL MPI_BCAST(RA, SIZE(RA), MPI_REAL8, 0, MPI_COMM_WORLD, IERR)
  CALL MPI_BCAST(UPNEIGH, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)

  IF (MOD(ISTEP - 1, 50) .EQ. 0 .and. UPNEIGH) CALL DEFECT_NEIGH_UPDATE(IDLO, IDHI)

  OC = 0
  OCC = 0
  OCCBUFF = 0

!!DETERMINE THE OCCUPANCY OF EACH SITE

!$OMP PARALLEL PRIVATE(RS,RX,RY,RZ,RMIN,MIN,J,N,ATOM) SHARED(OC,OCC)
!$OMP DO
  DO I = IDLO, IDHI

    RS = 0; RMIN = HUGE(RMIN); MIN = 0
    DO J = 1, NNEIGH
      RID = RNEIGH(J, I)
      RX = PBC(RA(1, I) - RR(1, RID), 1)
      RY = PBC(RA(2, I) - RR(2, RID), 2)
      RZ = PBC(RA(3, I) - RR(3, RID), 3)
      RS = RX*RX + RY*RY + RZ*RZ
      IF (RS .LT. RMIN) THEN
        MIN = RID
        RMIN = RS
      END IF
    END DO

!$OMP CRITICAL
    OCC(MIN) = OCC(MIN) + 1
    OC(OCC(MIN), MIN) = I
!$OMP END CRITICAL

  END DO
!$OMP END DO

!$OMP SINGLE
  CALL MPI_ALLREDUCE(OCC, OCCBUFF, SIZE(OCC), MPI_INT, MPI_SUM, MPI_COMM_WORLD, IERR)
!$OMP END SINGLE

!!USING OCCUPANCY TO DETERMINE DEFECTS
!$OMP DO
  DO I = 1, NATOMR

    N = OCC(I)
    IF (N .GT. 0) THEN
      TV(I) = N
      DO J = 1, N
        ATOM = OC(J, I)
        TI(ATOM) = OCCBUFF(I)
      END DO
    END IF

  END DO
!$OMP END DO
!$OMP END PARALLEL

  CALL MPI_REDUCE(TV, TVBUFF, SIZE(TV), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
  CALL MPI_REDUCE(TI, TIBUFF, SIZE(TI), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, IERR)

  IF (RANK .eq. 0) THEN
    TV = TVBUFF
    TI = TIBUFF
  END IF

END SUBROUTINE POLY_CRYSTAL_WS
!***********************************************************************!
SUBROUTINE DEFECT_NEIGH_UPDATE(LO, HI)

  USE SYSTEM
  USE DEFECT
  USE BINS, ONLY: maxNeighborSize

  IMPLICIT NONE

  !LO & HI ARE THE IDS FOR WHICH NEIGH LIST WILL BE UPDATED FOR
  INTEGER :: LO, HI
  INTEGER :: I, J, K
  REAL*8 :: RX, RY, RZ, RS, PBC
  REAL*8 START, FINISH
  REAL*8 :: RSLIST(NNEIGH)
  INTEGER :: RFNEIGH(maxNeighborSize)
  if (rank .eq. 0) WRITE (*, '(a)') 'Updating defect neighbor lists'

  UPNEIGH = .FALSE.

  CALL CPU_TIME(START)

  RX = 0; RY = 0; RZ = 0; RS = 0
  RSLIST = 0; RNEIGH = 0

  IF (.NOT. RFGRID%isInitialized) THEN
    CALL RFGRID%initBinGrid(20, 20, 20, RBOX*EPSX, RBOX*EPSY, RBOX*EPSZ)
    CALL RFGRID%assignAtomsToBins(RR, NATOMR)
  END IF

  !FOR EACH ATOM IN THE SYSTEM FIND N NEAREST REFERENCE ATOMS
  DO I = LO, HI
    RSLIST = HUGE(RS)
    RFNEIGH = RFGRID%getNearestNeighbors(RA(1, I), RA(2, I), RA(3, I))
    J = 1
    DO WHILE (RFNEIGH(J) .NE. 0 .AND. J .LE. SIZE(RFNEIGH))
      RX = PBC(RA(1, I) - RR(1, RFNEIGH(J)), 1)
      RY = PBC(RA(2, I) - RR(2, RFNEIGH(J)), 2)
      RZ = PBC(RA(3, I) - RR(3, RFNEIGH(J)), 3)
      RS = RX*RX + RY*RY + RZ*RZ
      IF (RS .lt. RSLIST(NNEIGH)) THEN
        RSLIST(NNEIGH) = RS
        RNEIGH(NNEIGH, I) = RFNEIGH(J)
        CALL quicksort(RSLIST, RNEIGH(:, I), 1, NNEIGH)
      END IF
      J = J + 1
    END DO
  END DO

  CALL CPU_TIME(FINISH)
  IF (RANK .eq. NPROCS - 1) WRITE (*, '(a,f14.5,a)') 'Neighbor list update time= ', FINISH - START, ' s'

END SUBROUTINE DEFECT_NEIGH_UPDATE
!***********************************************************************!
recursive subroutine quicksort(a, b, first, last)
  ! This quicksort routine sorts both a and b arrays based on the values in a
  implicit none
  real*8 a(*), x, t
  integer b(*)
  integer first, last
  integer i, j

  x = a((first + last)/2)
  i = first
  j = last
  do
    do while (a(i) .lt. x)
      i = i + 1
    end do
    do while (x .lt. a(j))
      j = j - 1
    end do
    if (i .ge. j) exit
    t = a(i); a(i) = a(j); a(j) = t
    t = b(i); b(i) = b(j); b(j) = t
    i = i + 1
    j = j - 1
  end do
  if (first .lt. i - 1) call quicksort(a, b, first, i - 1)
  if (j + 1 .lt. last) call quicksort(a, b, j + 1, last)
end subroutine quicksort
!***********************************************************************!
recursive subroutine quicksort_int(array1, array2, lo, hi)
    implicit none
    integer :: array1(*), array2(*)
    integer :: pivot, right, left, temp, n, lo, hi

    n = size(array1(lo:hi))

    if (n > 1) then
      pivot = array1(lo)
      left = lo
      right = hi
      do while (left .lt. right)
        do while (array1(left) < pivot)
          left = left + 1
        end do
        do while (array1(right) > pivot)
          right = right - 1
        end do
        if (left <= right) then
          temp = array1(left)
          array1(left) = array1(right)
          array1(right) = temp
          temp = array2(left)
          array2(left) = array2(right)
          array2(right) = temp
          left = left + 1
          right = right - 1
        end if
      end do

      call quicksort_int(array1, array2, lo, right)
      call quicksort_int(array1, array2, left, hi)

    end if

end subroutine quicksort_int
