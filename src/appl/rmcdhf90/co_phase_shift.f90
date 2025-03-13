!***********************************************************************
!                                                                      *
SUBROUTINE co_phase_shift(J)
!                                                                      *
!   Calculates phase shift of the continuum orbital, by comparing      *
!   it to the analytical behaviour for large r.
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       J : (input) Continuum orbital serial number                    *
!                                                                      *
!   Call(s) to: co_spherical_bessel                                    *
!                                                                      *
!   Paweł Syty                                                         *
!                                                           July 2024  *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE def_C, ONLY: C
    USE continuum_C
    USE wave_C, ONLY: PF, QF, MF
    USE grid_C, ONLY: R
    USE orb_C, ONLY: NAK
    USE iounit_C

    IMPLICIT none

    INTEGER, INTENT(in)  :: J

    INTEGER :: l, mfj, i, i_max, i_max1, i0_max
    INTEGER :: maxima_counter, maxima_counter0
    REAL(DOUBLE) :: r_max, y_max, r_max1, y_max1, r0_max, y0_max
    REAL(DOUBLE) :: bj, by
    REAL(DOUBLE) :: wavenumber, wavelength, wavelength_csp, wavelength_diff
    REAL(DOUBLE), PARAMETER :: PI = 4.0D0*ATAN(1.0D0)

    WRITE(ISTDE,*)

    IF (CO_ENERGY == 0.0D0) THEN
        WRITE(ISTDE,*) "Zero energy case, skipping phase shift calculation."
        RETURN
    END IF

    ! Determine the orbital quantum number l from kappa
    IF (NAK(J) > 0) THEN
        l = NAK(J)
    ELSE
        l = -(NAK(J) + 1)
    END IF

    mfj = MF(J)
    wavenumber = SQRT( 2*(-CO_ENERGY) + (-CO_ENERGY)/(C*C) )

    ! Determine the last maximum of the continuum wave function
    i_max = -1
    DO i = mfj-4, 1, -1
        IF (PF(i-1, J) < PF(i, J) .AND. PF(i+1, J) < PF(i, J)) THEN
            i_max = i
            EXIT
        END IF
    END DO
    IF (i_max <= 0) THEN
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (no maximum)"
        RETURN
    ENDIF
    CALL co_getmax(R(i_max-4:i_max+4), PF(i_max-4:i_max+4, J), 9, r_max, y_max)
    IF (r_max == 0 .AND. y_max == 0) THEN
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (no maximum)"
        RETURN
    ENDIF
    ! Enumerate the maximum
    maxima_counter = 0
    DO i = 1, i_max-1
        IF (R(i) <= CO_R0_3) CYCLE
        IF (PF(i-1, J) < PF(i, J) .AND. PF(i+1, J) < PF(i, J)) maxima_counter = maxima_counter + 1
    END DO

    ! Determine the one before last maximum of the continuum wave function
    ! (for estimation of the quality of calculated continuum spinor)
    i_max1 = -1
    DO i = i_max-4, 1, -1
        IF (PF(i-1, J) < PF(i, J) .AND.  PF(i+1, J) < PF(i, J)) THEN
            i_max1 = i
            EXIT
        END IF
    END DO
    IF (i_max1 <= 0) THEN
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (no maximum)"
        RETURN
    ENDIF
    CALL co_getmax(R(i_max1-4:i_max1+4), PF(i_max1-4:i_max1+4, J), 9, r_max1, y_max1)
    IF (r_max1 == 0 .AND. y_max1 == 0) THEN
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (no maximum)"
        RETURN
    ENDIF

    ! Check if the continuum spinor has been calculated far enough
    ! from the origin, by comparing the analytical and numerical wavelengths
    wavelength = 2 * PI / wavenumber
    wavelength_csp = r_max - r_max1
    wavelength_diff = abs(wavelength - wavelength_csp) / wavelength * 100
    WRITE(*,'(A,F13.8,A)') " Difference between analytical and numerical &
                   wavelength of the continuum spinor = ",wavelength_diff,"%"
    WRITE(ISTDE,'(A,F10.2,A,F7.2)') " Radial grid: Rmax = ", R(mfj), ", Rstep = ", R(mfj)-R(mfj-1)
    IF (wavelength_diff > 0.01D0) THEN
        CO_GRID_TOO_SHORT = .TRUE.
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (too short radial grid)"
        RETURN
    ELSE
        WRITE(ISTDE,*) "Radial grid is long enough to calculate the phase shift."
        CO_GRID_TOO_SHORT = .FALSE.
    END IF


    ! Construct the analytical "zero potential" solution on the grid
    DO i = 1, mfj
        CALL co_spherical_bessel(l, wavenumber*R(i), bj, by)
        CO_CSP_ZERO(i) =  R(i) * bj
    END DO

    ! Find the maximum of the "zero potential" solution
    i0_max = -1
    DO i = mfj-4, 1, -1
        IF (CO_CSP_ZERO(i-1) < CO_CSP_ZERO(i) .AND. CO_CSP_ZERO(i+1) < CO_CSP_ZERO(i)) THEN
            i0_max = i
            EXIT
        END IF
    END DO
    IF (i0_max <= 0) THEN
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (no maximum)"
        RETURN
    END IF
    CALL co_getmax(R(i0_max-4:i0_max+4), CO_CSP_ZERO(i0_max-4:i0_max+4), 9, r0_max, y0_max)
    IF (r0_max == 0 .AND. y0_max == 0) THEN
        WRITE(ISTDE,*) "*** Warning: skipping phase shift calculation (no maximum)"
        RETURN
    ENDIF
    ! Enumerate the maxima
    maxima_counter0 = 0
    DO i = 1, i0_max-1
        IF (R(i) <= CO_R0_3) CYCLE
        IF (CO_CSP_ZERO(i-1) < CO_CSP_ZERO(i) .AND. CO_CSP_ZERO(i+1) < CO_CSP_ZERO(i)) &
            maxima_counter0 = maxima_counter0 + 1
    END DO

    ! Phase shift calculation (with adjustments in case the maxima do not match)
    IF (maxima_counter0 == maxima_counter) CO_PS = wavenumber * (r0_max - r_max)
    IF (maxima_counter0 > maxima_counter)  CO_PS = wavenumber * ((r0_max-wavelength) - r_max)
    IF (maxima_counter0 < maxima_counter)  CO_PS = wavenumber * (r0_max - (r_max-wavelength))

    ! Shift the phase shift
    DO WHILE (CO_PS < -PI  .OR. CO_PS > PI)
        IF (CO_PS < -PI) CO_PS = CO_PS + 2*PI
        IF (CO_PS > PI)  CO_PS = CO_PS - 2*PI
    END DO

    ! Shift to the other range
    CO_PS_SHIFTED = CO_PS
    DO WHILE (CO_PS_SHIFTED < -PI/2  .OR. CO_PS_SHIFTED > PI/2)
        IF (CO_PS_SHIFTED < -PI/2) CO_PS_SHIFTED = CO_PS_SHIFTED + PI
        IF (CO_PS_SHIFTED >  PI/2) CO_PS_SHIFTED = CO_PS_SHIFTED - PI
    END DO

    IF (CO_PS /= CO_PS_SHIFTED) THEN
        WRITE(ISTDE,*) "Phase shift = ",CO_PS," = ", CO_PS_SHIFTED
    ELSE
        WRITE(ISTDE,*) "Phase shift = ",CO_PS
    END IF

    RETURN

END SUBROUTINE co_phase_shift
