!***********************************************************************
!                                                                      *
MODULE continuum_C
!                                                                      *
! Paweł Syty, 17/04/2019 (work started...)                             *
!             22/06/2024 (last modification)                           *
!                                                                      *
!***********************************************************************
    USE vast_kind_param, ONLY: DOUBLE
    USE parameter_def,    ONLY: NNNP

! Perform or not calculations of continuum wave function
    LOGICAL :: CO_CALCULATE = .FALSE.
! Perform or not final normalization of the continuum wave function
    LOGICAL :: CO_NORMALIZE = .FALSE.

! Serial number of the continuum orbital, determined automatically
    INTEGER :: CO_ORBITAL
! Kinetic energy (hartree) of the continuum electron (negative!)
    REAL(DOUBLE) :: CO_ENERGY

! Polarization potential parameters
    LOGICAL :: CO_INCLUDE_POLARIZATION = .TRUE.
    REAL(DOUBLE), DIMENSION(NNNP) :: CO_VPOL  ! Polarization potential
! Flags indicating how the parameters of the potential will be calculated
    LOGICAL :: CO_POLARIZATION_AUTO = .FALSE.
    LOGICAL :: CO_POLARIZATION_MANUAL = .FALSE.
    LOGICAL :: CO_POLARIZATION_FROM_FILE = .FALSE.
! Dipole polarization potential parameters (polarizability, <cutoff^3>)
    REAL(DOUBLE) :: CO_ALPHA_D = 0.0D0, CO_R0_3 = 0.0D0
! Quadrupole polarization potential parameters (polarizability, <cutoff^5>)
    REAL(DOUBLE) :: CO_ALPHA_Q = 0.0D0, CO_R0_5 = 0.0D0

! Calculated properties of the continuum wave function
    REAL(DOUBLE) :: CO_SL ! Scattering length, calculated if CO_ENERGY = 0
    REAL(DOUBLE) :: CO_PHASE_SHIFT ! Phase shift (TO BE DONE)

! Other variables
    CHARACTER(LEN=200) :: CO_DUMMY ! Dummy string
    INTEGER :: CO_ANSWER ! For collecting data

END MODULE continuum_C
