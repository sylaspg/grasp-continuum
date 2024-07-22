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

    IF (CO_ENERGY == 0) THEN

        ! Calculation using formula dervided from the straight line for zero energy.
        ! It is based on last two points of continuum orbital: (MFJ-1, MFJ).
        ! Result is compared to that calculated using points (MFJ-2, MFJ-1).

        mfj = MF(J)

        a = ( (PF(mfj,J) - PF(mfj-1,J)) / (R(mfj) - R(mfj-1)) )
        b = PF(mfj-1,J) - a*R(mfj-1)
        CO_SL = -b / a

        a = ( (PF(mfj-1,J) - PF(mfj-2,J)) / (R(mfj-1) - R(mfj-2)) )
        b = PF(mfj-2,J) - a*R(mfj-2)
        sl2 = -b / a
        sl_diff = ABS( (CO_SL - sl2) / CO_SL * 100 )

        WRITE(ISTDE,'(A,F10.2,A,F7.2)') " Radial grid: Rmax = ", R(mfj), ", Rstep = ", R(mfj)-R(mfj-1)
        WRITE(ISTDE,*) "Scattering length = ", CO_SL, "(diff = ", sl_diff, "% )"

    ELSE

        ! Ordinary formula based on the phase shift (used, if grid is long enough)
        IF (.NOT. CO_GRID_TOO_SHORT) THEN

            wavenumber = SQRT(2*(-CO_ENERGY) + (-CO_ENERGY)/(C*C))
            CO_SL = -1/wavenumber * TAN(CO_PS)
            WRITE(ISTDE,*) "Scattering length = ", CO_SL, "(estimation only - non-zero energy)"

        END IF

    END IF

    RETURN

END SUBROUTINE co_scattering_length