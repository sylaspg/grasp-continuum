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
!                                           Last update: 19 Apr 2019   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE continuum_C
    USE orb_C
    USE wave_C
    USE grid_C

    IMPLICIT none

    INTEGER, INTENT(in)  :: J
    INTEGER :: i, i_mtp_max, MFJ
    REAL(DOUBLE) :: E_RYDBERGS, r0, new_amplitude, old_amplitude, factor

    REAL (DOUBLE), PARAMETER :: PI = 4.0D00*ATAN(1.0D00)

    print*, "Normalization of the continuum orbital"
    MFJ = MF(J)
    ! Determine the last maximum
    i_mtp_max = -1
    do  i = MFJ-4,1,-1
        IF (PF(i-1,CO_ORBITAL) < PF(i,CO_ORBITAL)  .AND.  &
            PF(i+1,CO_ORBITAL) < PF(i,CO_ORBITAL))  THEN
            i_mtp_max = i
            EXIT
        END IF
    END DO

    IF (i_mtp_max > 0) THEN

        E_RYDBERGS = -CO_ENERGY * 2  ! Hartree to Rydbergs

        CALL co_getmax(R(i_mtp_max-2:i_mtp_max+2), &
                         PF(i_mtp_max-2:i_mtp_max+2,CO_ORBITAL), &
                         5,r0,old_amplitude)

        IF (r0 == 0 .AND. old_amplitude == 0) THEN
            PRINT*,"No maximum found, skipping continuum normalization"
            RETURN
        ENDIF

        IF (old_amplitude < 0) old_amplitude = abs(old_amplitude)
        new_amplitude = PI**(-0.5) * E_RYDBERGS**(-0.25)
        factor = new_amplitude / old_amplitude
        !print*, "  MTP = ", i_mtp_max
        !print*, "  Old amplitude = ", old_amplitude
        !print*, "  New amplitude = ", new_amplitude
        !print*, "  Normalization factor (applied to both P and Q) = ", factor

        PF(:,CO_ORBITAL) = PF(:,CO_ORBITAL) * factor
        QF(:,CO_ORBITAL) = QF(:,CO_ORBITAL) * factor
    ELSE
        PRINT*,"No maximum found, skipping continuum normalization"
    ENDIF

END SUBROUTINE co_normalization
