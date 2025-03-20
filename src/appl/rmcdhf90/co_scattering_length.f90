!***********************************************************************
!                                                                      *
SUBROUTINE co_scattering_length(J)
!                                                                      *
!   Calculates scattering length.                                      *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       J : (input) Continuum orbital serial number                    *
!                                                                      *
!   Paweł Syty                                                         *
!                                                          July 2024   *
!                                            Last update: March 2025   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE continuum_C
    USE def_C, ONLY: C
    USE wave_C, ONLY: PF, QF, MF
    USE grid_C, ONLY: R
    USE iounit_C

    IMPLICIT none

    INTEGER, INTENT(in) :: J

    INTEGER :: mfj
    REAL(DOUBLE) :: a, b, sl2, sl_diff
    REAL(DOUBLE) :: wavenumber
    REAL(DOUBLE), PARAMETER :: PI = 4.0D0*ATAN(1.0D0)

    ! For linear regression
    INTEGER, PARAMETER :: n_points = NNNP/10
    REAL(DOUBLE) :: x(n_points), y(n_points)
    REAL(DOUBLE) :: aa, bb, sum_x, sum_y, sum_xx, sum_xy, RMSE
    INTEGER :: i

    WRITE(ISTDE,*)

    mfj = MF(J)

    IF (CO_ENERGY == 0) THEN

        ! Estimating the error of assuming that the tail is a straight line

        ! Collect last 'n_points' from the 'zero energy' wave function
        !IF (n_points > mfj) n_points = mfj
        DO i = 1, n_points
            x(i) = R(mfj-i+1)
            y(i) = PF(mfj-i+1,J)
        END DO

        ! Determine a and b coefficients (y = aa*x + bb)
        sum_x  = sum(x)
        sum_y  = sum(y)
        sum_xx = sum(x * x)
        sum_xy = sum(x * y)
        aa = (n_points * sum_xy - sum_x * sum_y) / (n_points * sum_xx - sum_x * sum_x)
        bb = (sum_xx * sum_y - sum_xy * sum_x) / (n_points * sum_xx - sum_x * sum_x)

        ! Calculate RMSE
        RMSE = 0.0D0
        DO i = 1, n_points
            RMSE = RMSE + (PF(mfj-i+1,J) - (aa*R(mfj-i+1) + bb))**2
        END DO
        RMSE = SQRT(RMSE / n_points)

        ! Calculation using formula dervided from the straight line for zero energy.
        ! It is based on last two points of continuum orbital: (MFJ-1, MFJ).
        ! Result is compared to that calculated using points (MFJ-2, MFJ-1).

        a = ( (PF(mfj,J) - PF(mfj-1,J)) / (R(mfj) - R(mfj-1)) )
        b = PF(mfj-1,J) - a*R(mfj-1)
        CO_SL = -b / a

        a = ( (PF(mfj-1,J) - PF(mfj-2,J)) / (R(mfj-1) - R(mfj-2)) )
        b = PF(mfj-2,J) - a*R(mfj-2)
        sl2 = -b / a
        sl_diff = ABS( (CO_SL - sl2) / CO_SL * 100 )

        WRITE(ISTDE,*) "Scattering length = ", CO_SL
        WRITE(ISTDE,*) "Information for error estimation:"
        WRITE(ISTDE,'(A,F10.2,A,F7.2)') "   Radial grid: Rmax = ", R(mfj), ", Rstep = ", R(mfj)-R(mfj-1)
        WRITE(ISTDE,'(A,F11.8,A1)') "   Relative error of straight line determination at the tail         = ", &
                                     ABS(RMSE/PF(mfj,J)*100), "%"
        WRITE(ISTDE,'(A,F11.8,A1)') "   Difference between two consecutive estimates of scattering length = ", &
                                     sl_diff, "%"
        WRITE(ISTDE,*) "  (calculated from the last two points on the grid, and the two penultimate points)"

    ELSE

        ! Ordinary formula based on the phase shift (used, if grid is long enough)
        IF (.NOT. CO_GRID_TOO_SHORT) THEN

            wavenumber = SQRT(2*(-CO_ENERGY) + (-CO_ENERGY)/(C*C))
            CO_SL = -1/wavenumber * TAN(CO_PS)
            WRITE(ISTDE,*) "Scattering length = ", CO_SL
            WRITE(ISTDE,*) "(this is just an estimate, calculated from the phase shift of non-zero energy continuum spinor;"
            WRITE(ISTDE,*) "to have a more accurate result, repeat the calculation with zero energy)"

        END IF

    END IF

    RETURN

END SUBROUTINE co_scattering_length