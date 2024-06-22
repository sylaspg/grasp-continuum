!***********************************************************************
!                                                                      *
SUBROUTINE co_scattering_length(J)
!                                                                      *
!   Calculates scattering length, based at last two points of          *
!   continuum orbital: (MFJ-1, MFJ).                                   *
!   Result is compared to that calculated using points (MFJ-2, MFJ-1). *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       J : (input) Continuum orbital serial number                    *
!                                                                      *
!   Paweł Syty                                                         *
!                                           Last update: 10 May 2024   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE continuum_C
    USE wave_C
    USE grid_C

    IMPLICIT none

    INTEGER, INTENT(in)  :: J
    INTEGER :: mfj

    REAL(DOUBLE) a, b, sl2, sl_diff

    mfj = MF(J)

    a = ( (PF(mfj,J) - PF(mfj-1,J)) / (R(mfj) - R(mfj-1)) )
    b = PF(mfj-1,J) - a*R(mfj-1)
    CO_SL = -b / a

    a = ( (PF(mfj-1,J) - PF(mfj-2,J)) / (R(mfj-1) - R(mfj-2)) )
    b = PF(mfj-2,J) - a*R(mfj-2)
    sl2 = -b / a
    sl_diff = ABS((CO_SL - sl2)/CO_SL *100)

    PRINT*,"Scattering length = ", CO_SL, "(diff = ", sl_diff, "% )"
    WRITE(*,'(A,F10.2,A,F7.2)') " Grid: Rmax = ", R(mfj), ", Rstep = ", &
                                R(mfj)-R(mfj-1)

    RETURN

END SUBROUTINE co_scattering_length