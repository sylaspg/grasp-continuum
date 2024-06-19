!***********************************************************************
!                                                                      *
    SUBROUTINE co_scattering_length(J)
!                                                                      *
!   Calculates scattering length, based at last two points of          *
!   continuum orbital: (MFJ-1, MFJ).                                    *
!   Result is compared to that calculated using points (MFJ-2, MFJ-1). *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       J : (input) Continuum orbital number                           *
!                                                                      *
!   Paweł Syty                                                         *
!                                           Last update: 10 May 2024   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE wave_C
    USE grid_C

    IMPLICIT none

    INTEGER, INTENT(in)  :: J
    INTEGER :: MFJ

    REAL(DOUBLE) a, b, sl1, sl2, sl_diff

    MFJ = MF(J)
    ! MFJ = 4050

    a = ( (PF(MFJ,J) - PF(MFJ-1,J)) / (R(MFJ) - R(MFJ-1)) )
    b = PF(MFJ-1,J) - a*R(MFJ-1)
    sl1 = -b / a

    a = ( (PF(MFJ-1,J) - PF(MFJ-2,J)) / (R(MFJ-1) - R(MFJ-2)) )
    b = PF(MFJ-2,J) - a*R(MFJ-2)
    sl2 = -b / a
    sl_diff = sl1 - sl2

    PRINT*," Rmax = ",R(MFJ), "Rstep = ", R(MFJ)-R(MFJ-1)
    PRINT*," Scattering length = ", sl1, "(diff = ", sl_diff, " )"

    RETURN

END SUBROUTINE co_scattering_length