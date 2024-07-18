!***********************************************************************
!                                                                      *
SUBROUTINE co_save_to_file(J)
!                                                                      *
!   Saves continuum orbital to a file.                                 *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       J : (input) Continuum orbital number                           *
!                                                                      *
!   Paweł Syty                                                         *
!                                                           May 2024   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE continuum_C
    USE wave_C, ONLY: PF, QF, MF
    USE grid_C, ONLY: R
    USE orb_C, ONLY: NAK

    IMPLICIT none

    INTEGER, INTENT(in)  :: J

    INTEGER :: i
    CHARACTER(20), PARAMETER :: FNAME = "continuum.csp"


    OPEN(61,FILE=FNAME)
    WRITE(61,*) 'plot "-" index  0 using 1:2 with lines'    ! For gnuplot
    WRITE(61,*) '#'
    WRITE(61,*) '# Energy (hartree)  = ',-CO_ENERGY
    WRITE(61,*) '# kappa             = ',NAK(J)
    IF (.NOT. CO_GRID_TOO_SHORT) WRITE(61,*) '# Phase shift       = ',CO_PS
    IF (.NOT. CO_GRID_TOO_SHORT .OR. CO_ENERGY == 0.0D0) WRITE(61,*) '# Scattering length = ',CO_SL
    WRITE(61,*) '#'
    WRITE(61,*) '# r                         P(r)                      Q(r)'

    DO i=1, MF(J)
        WRITE(61,*) R(i), PF(i,J), QF(i,J)
    END DO
    CLOSE(61)

    PRINT*,"Continuum orbital saved to '", TRIM(FNAME), "' file."

    RETURN

END SUBROUTINE co_save_to_file