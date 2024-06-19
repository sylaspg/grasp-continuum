!***********************************************************************
!                                                                      *
    MODULE continuum_C
!                                                                      *
! Paweł Syty, 17/04/2019 (work started...)                             *
!             13/05/2024 (last modification)                           *
!                                                                      *
!***********************************************************************
    USE vast_kind_param, ONLY: DOUBLE
    USE parameter_def,    ONLY: NNNP

    LOGICAL :: CO_CALCULATE     ! Perform or not calculations for continuum orbital
    LOGICAL :: CO_NORMALIZE     ! Perform or not normalization of the continuum orbital
    LOGICAL :: CO_SAVE = .TRUE. ! Save or not continuum orbital to a file

    INTEGER :: CO_ORBITAL     ! Serial number of the continuum orbital
    REAL(DOUBLE) :: CO_ENERGY ! Kinetic energy (hartree) of the continuum electron (negative!)

! Polarization potental parameters
    LOGICAL :: CO_INCLUDE_POLARIZATION
    REAL(DOUBLE) :: CO_ALPHA_D, CO_R0_3       ! Dipole polarization potential parameters (polarizability, <cutoff^3>)
    REAL(DOUBLE) :: CO_ALPHA_Q, CO_R0_5       ! Quadrupole polarization potential parameters (polarizability, <cutoff^5>)
    REAL(DOUBLE), DIMENSION(NNNP) :: CO_VPOL  ! Polarization potential

! Other variables
    CHARACTER(LEN=200) :: CO_DUMMY ! Dummy string


    CONTAINS

!***********************************************************************
!                                                                      *
SUBROUTINE co_set_polarization_potential_parameters(atomic_number)
!                                                                      *
!   Sets polarization potential parameters for given atomic number     *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       atomic_number : (input)                                        *
!                                                                      *
!   Paweł Syty                                                         *
!                                           Last update: 10 May 2024   *
!                                                                      *
!***********************************************************************

    IMPLICIT none

    REAL(DOUBLE) atomic_number

    SELECT CASE (INT(atomic_number))

    CASE (4) ! Neon
        CO_ALPHA_D = 37.74D0
        CO_R0_3    = 28.49D0

    CASE (10) ! Neon
        CO_ALPHA_D = 2.6611D0
        CO_R0_3    = 2.02385D0

    CASE (18) ! Argon
        CO_ALPHA_D = 11.083D0
        CO_R0_3    = 7.65696D0

    CASE (30) ! Zinc
        CO_ALPHA_D = 38.67D0
        CO_R0_3    = 33.43D0

    CASE (36) ! Krypton
        CO_ALPHA_D = 16.78D0
        CO_R0_3    = 10.953D0

    CASE (38) ! Strontium
        CO_ALPHA_D = 197.2D0
        CO_R0_3    = 119.383D0   !139.641D0

    CASE (54) ! Xenon
        CO_ALPHA_D = 27.32D0
        CO_R0_3    = 19.3023D0

    CASE (56) ! Barium
        CO_ALPHA_D = 272.0D0        ! 272 += 10
        CO_R0_3    = 189.0D0

    CASE (80) ! Mercury
        CO_ALPHA_D = 33.91D0
        CO_R0_3    = 33.0D0        ! SL = 1.8 dla r0^3 = 600, 0.65

    CASE (86) ! Radon
        CO_ALPHA_D = 35.0D0
        CO_R0_3    = 25.132D0

    CASE (118) ! Oganesson
        CO_ALPHA_D = 58.0D0
        CO_R0_3    = 38.898D0

    CASE DEFAULT
        STOP "co_set_polarization_potential_parameters() - Unknown element for polarization parameters."

    END SELECT

    END SUBROUTINE co_set_polarization_potential_parameters


    END MODULE continuum_C
