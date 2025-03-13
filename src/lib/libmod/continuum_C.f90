!***********************************************************************
!                                                                      *
MODULE continuum_C
!                                                                      *
!  This module contains definitions needed for calculations of         *
!  continuum orbitals. All names starts with 'CO_' prefix.             *
!                                                                      *
!   Paweł Syty                                                         *
!                                               Gdańsk, 2019 - 2025    *
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

! Array for the non-disturbed wave
    REAL(DOUBLE), DIMENSION(NNNP) :: CO_CSP_ZERO

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
    ! Octupole polarization potential parameters (polarizability, <cutoff^7>)
    REAL(DOUBLE) :: CO_ALPHA_O = 0.0D0, CO_R0_7 = 0.0D0

! Properties of the continuum wave function
    REAL(DOUBLE) :: CO_SL ! Scattering length
    REAL(DOUBLE) :: CO_PS, CO_PS_SHIFTED ! Phase shifts in two ranges
    LOGICAL :: CO_GRID_TOO_SHORT = .TRUE. ! True, if radial grid is too short for phase shift calculations

! Other variables
    CHARACTER(LEN=200) :: CO_DUMMY ! Dummy string
    INTEGER :: CO_ANSWER ! For collecting data

END MODULE continuum_C
