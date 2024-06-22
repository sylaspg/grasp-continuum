!***********************************************************************
!                                                                      *
SUBROUTINE co_save_to_file(J)
!                                                                      *
!   Saves continuum orbital to file                                    *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       J : (input) Continuum orbital number                           *
!                                                                      *
!   Paweł Syty                                                         *
!                                           Last update: 13 May 2024   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE wave_C
    USE grid_C

    IMPLICIT none

    INTEGER, INTENT(in)  :: J
    INTEGER :: i
    CHARACTER(20), PARAMETER :: FNAME = "continuum.csp"

    OPEN(61,FILE=FNAME)
! For gnuplot
    WRITE(61,*) 'plot "-" index  0 using 1:2 with lines'
    WRITE(61,*) '#     r                   P(r)                   Q(r)'
    DO i=1, MF(J)
        WRITE(61,*) R(i), pf(i,J), qf(i,J)
    END DO
    CLOSE(61)

    PRINT*,"Continuum orbital saved to '", TRIM(FNAME), "' file."

    RETURN

END SUBROUTINE co_save_to_file