!***********************************************************************
!                                                                      *
SUBROUTINE co_spherical_bessel(L,Z,J,Y)
!                                                                      *
!   Calculates spherical bessel functions of first and second kind.    *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!       L, Z: (input) arguments                                        *
!       J, Y: (output) calculated J and Y values                       *
!                                                                      *
!   Paweł Syty                                                         *
!                                                           July 2024  *
!                                                                      *
!***********************************************************************

        USE vast_kind_param, ONLY: DOUBLE

        IMPLICIT none

        INTEGER, INTENT(in)       :: L
        REAL(DOUBLE), INTENT(in)  :: Z
        REAL(DOUBLE), INTENT(out) :: J, Y

        INTEGER      :: K
        REAL(DOUBLE) :: J0, Y0, J1, Y1

        J0 =  SIN(Z)/Z
        Y0 = -COS(Z)/Z

        IF (L <= 0) THEN
            J = J0
            Y = Y0
            RETURN
        END IF

        J1 =  Y0 + J0/Z
        Y1 = -J0 + Y0/Z

        IF(L <= 1) THEN
            J = J1
            Y = Y1
            RETURN
        END IF

        DO K=2,L
            J = (2*K - 1)*J1/Z - J0
            Y = (2*K - K)*Y1/Z - Y0
            J0 = J1
            Y0 = Y1
            J1 = J
            Y1 = Y
        END DO

        RETURN

    END SUBROUTINE co_spherical_bessel
