!***********************************************************************
!                                                                      *
      SUBROUTINE co_initrwf(INDEX, NSUBS)
!                                                                      *
!   This subroutine estimates the wave function for calculations of    *
!   the continuum orbital.                                             *
!                                                                      *
!   If there are convergence problems, tune values of the parameters   *
!   PF_INIT, QF_INIT and PZ_INIT, or use the other estimation method.  *
!
!   Call(s) to: [LIB92]: START                                         *
!                                                                      *
!   Paweł Syty                                                         *
!                                                          July 2024   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE DEF_C
    USE parameter_def, ONLY: NNNW
    USE LEFT_C, ONLY: SET
    USE WAVE_C, ONLY: PF, QF, PZ
    USE WHFROM_C, ONLY: SOURCE

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: INDEX(NNNW)
    INTEGER, INTENT(IN) :: NSUBS

    REAL(DOUBLE), PARAMETER :: PF_INIT = 1.0D-7
    REAL(DOUBLE), PARAMETER :: QF_INIT = 1.0D-12
    REAL(DOUBLE), PARAMETER :: PZ_INIT = 10.0D0

    INTEGER :: J, LOC
    REAL(DOUBLE) :: Q0

    DO J = 1, NSUBS
        LOC = INDEX(J)
        IF (SET(LOC)) CYCLE
        PZ(LOC) = PZ_INIT

        ! Generate 6 starting points
        CALL START (J, 1, PZ(J), PF(:,J), Q0, QF(:,J))

        ! Constant values to the end of grid
        PF(7:,J) = PF_INIT
        QF(7:,J) = QF_INIT

        SOURCE(LOC) = 'con'
        SET(LOC) = .TRUE.
    END DO

END SUBROUTINE co_initrwf
