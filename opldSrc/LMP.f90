!***********************************************************************!
module LMP_DRIVER

  use LAMMPS

  implicit none

!lmp: interactable LAMMPS object
!lmpDimer: LAMMPS object for dimer searches
  type(C_ptr) :: lmp
  type(C_ptr) :: lmpDimer
!lmpTemp: current measurement of system temperature
  real(C_double), pointer :: lmpTemp => NULL()
!lmpEnergy: current measurement of system energy
  real(C_double), pointer :: lmpEnergy => NULL()
!lmpPressure: current measurement of system pressure tensor
  real(C_double), dimension(:), pointer :: lmpPressure => NULL()
!lmpL*: length of the system in the respective dimensions
  real*8, dimension(3) :: lmpBoxLo = 0, lmpBoxHi = 0
  real*8 :: xy = 0, yz = 0, xz = 0
  logical, dimension(3) :: periodicity = 0
  logical :: box_change = 0

  character(len=200) :: inputLAMMPSFilename
  character(len=200) :: configurationFilename
!***********************************************************************!
contains
!***********************************************************************!
  subroutine LMP_INIT

    use MPI
    use SYSTEM

    implicit none
    character(len=200) :: str

    !attaches lammps object to pointer variable lmp
    call lammps_open('lmp -log none -screen none', MPI_COMM_WORLD, lmp)

    !Setting up system in LAMMPS space

    !Necessary commands for OPLD to work
    call lammps_command(lmp, 'units metal')
    call lammps_command(lmp, 'boundary p p p')
    call lammps_command(lmp, 'atom_style atomic')
    call lammps_command(lmp, 'atom_modify map array')

    !Read from input file using LAMMPS
    call lammps_file(lmp, trim(inputLAMMPSFilename))

    !Delete atoms in the system and replace with 'configurationFilename'
    call lammps_command(lmp, 'delete_atoms group all')
    call lammps_command(lmp, 'read_dump ' // trim(configurationFilename) // ' 0 x y z box no add yes format xyz')
  
    !Minimize the imported system
    call lammps_command(lmp, 'minimize 1.0e-6 1.0e-6 1000 10000')
    !Including necessary computes for OPLD
    call lammps_command(lmp, 'compute TotEn all pe')
    call lammps_command(lmp, 'compute TotPress all pressure NULL virial')

    !Initiating system values
    call lammps_command(lmp, 'run 0')

    !Extracting MD values for cascade
    call lammps_extract_compute(lmpTemp, lmp, 'thermo_temp', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
    call lammps_extract_compute(lmpEnergy, lmp, 'thermo_pe', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
    call lammps_extract_compute(lmpPressure, lmp, 'thermo_press', LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR)

    !!SETTING UP COMMS FOR DIMER LAMMPS OBJECTS
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, RANK, RANK, DIMER_COMM, IERR)

    CALL LMP_DIMER_INIT

  end subroutine LMP_INIT
!***********************************************************************!
  subroutine LMP_DIMER_INIT

    use MPI
    use SYSTEM

    implicit none

    character(len=200) :: str

    !attaches lammps object to pointer variable lmp
    call lammps_open('lmp -log none -screen none', DIMER_COMM, lmpDimer)

    !Necessary commands for OPLD to work
    call lammps_command(lmpDimer, 'units metal')
    call lammps_command(lmpDimer, 'boundary p p p')
    call lammps_command(lmpDimer, 'atom_style atomic')
    call lammps_command(lmpDimer, 'atom_modify map array')

    !Read from input file using LAMMPS
    call lammps_file(lmpDimer, trim(inputLAMMPSFilename))

    !Dimer lammps object initially has no atoms in it
    call lammps_command(lmpDimer, 'delete_atoms group all')

  end subroutine LMP_DIMER_INIT
!***********************************************************************!
  subroutine RECV_LMP_POS

    use SYSTEM, only: NATOMC, RA, RLAT

    implicit none

    integer :: i
    real*8, dimension(:), allocatable :: pos

    call lammps_gather_atoms(lmp, 'x', 3, pos)

    DO I = 1, NATOMC
      RA(:, i) = pos(3*(i - 1) + 1:3*(i - 1) + 3)/RLAT
    end do

  end subroutine RECV_LMP_POS
!***********************************************************************!
  subroutine SEND_OPLD_POS

    use SYSTEM, only: NATOMC, RA, RLAT, IERR
    use MPI

    implicit none

    integer :: i
    real*8, dimension(:) :: pos(3*NATOMC)

    CALL MPI_BCAST(RA, size(RA), MPI_REAL8, 0, MPI_COMM_WORLD, IERR)

    DO I = 1, NATOMC
      pos(3*(i - 1) + 1:3*(i - 1) + 3) = RA(:, i)*RLAT
    END DO

    CALL lammps_scatter_atoms(lmp, 'x', pos)

  end subroutine SEND_OPLD_POS
!***********************************************************************!
end module LMP_DRIVER
