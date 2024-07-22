!***********************************************************************
!                                                                      *
SUBROUTINE co_normalization(J)
!                                                                      *
!   Normalizes the continuum spinor, acoording to                      *
!   Robert D. Cowan, The Theory of Atomic Structure and Spectra,       *
!   Univ. California Press, Berkeley 1981, page 519.                   *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!        J : (input) Continuum orbital number                          *
!                                                                      *
!   Call(s) to: [RMCDH90]: co_getmax                                   *
!                                                                      *
!   Paweł Syty                                                         *
!                                                         April 2019   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE continuum_C
    USE orb_C
    USE wave_C
    USE grid_C
    USE iounit_C

    IMPLICIT none

    INTEGER, INTENT(in)  :: J

    INTEGER :: i, i_max, mfj
    REAL(DOUBLE) :: e_rydbergs, r0, new_amplitude, old_amplitude, factor
    REAL(DOUBLE), PARAMETER :: PI = 4.0D0*ATAN(1.0D0)

    IF (CO_ENERGY == 0) THEN
        WRITE(ISTDE,*) "*** Skipping normalization of continuum wave function (energy = 0)"
        RETURN
    END IF

    mfj = MF(J)

    ! Determine the last maximum of the continuum wave function
    i_max = -1
    DO i = mfj-4, 1, -1
        IF (PF(i-1, J) < PF(i, J) .AND. PF(i+1, J) < PF(i, J))  THEN
            i_max = i
            EXIT
        END IF
    END DO

    IF (i_max > 0) THEN

        e_rydbergs = -CO_ENERGY * 2  ! Hartree to Rydbergs

        CALL co_getmax(R(i_max-4:i_max+4), PF(i_max-4:i_max+4, J), 9, r0, old_amplitude)

        IF (r0 == 0 .AND. old_amplitude == 0) THEN
            WRITE(ISTDE,*) "*** Skipping normalization of continuum wave function (no maximum)"
            RETURN
        ENDIF

        IF (old_amplitude < 0) old_amplitude = abs(old_amplitude)
        new_amplitude = PI**(-0.5) * e_rydbergs**(-0.25)
        factor = new_amplitude / old_amplitude
        WRITE(ISTDE,*) "Performing normalization of the continuum wave function."
        PF(:, J) = PF(:, J) * factor
        QF(:, J) = QF(:, J) * factor
    ELSE

        WRITE(ISTDE,*) "*** Skipping normalization of continuum wave function (no maximum)"

    ENDIF

END SUBROUTINE co_normalization
