PROGRAM MAIN

  USE MPI
  USE LMP_DRIVER
  USE SYSTEM, ONLY: TJOB, IERR, RANK, NPROCS, NDIMER

  IMPLICIT NONE

  call getarg(1, inputLAMMPSFilename)
  call getarg(2, configurationFilename)

  CALL MPI_INIT(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)

  CALL READ_INPUT
  CALL LMP_INIT
  CALL SYS_INPUT
  CALL SYS_INIT

  IF (MOD(NDIMER, NPROCS) .NE. 0) THEN

    WRITE (*, *) "# MOD(NDIMER, NPROCS) .NE. 0"
    WRITE (*, *) "PLEASE CHANGE NDIMER OR NUM OF PROCESSES"
    CALL MPI_FINALIZE(IERR)
    WRITE (*, *) "STOPPING PROGRAM"
    STOP

  ELSE

    CALL OPLD

  END IF

  CALL MPI_FINALIZE(IERR)

END PROGRAM MAIN
