!***********************************************************************
!                                                                      *
SUBROUTINE co_add_polarization_potential()
!                                                                      *
!   Adds polarization potential to the direct potential.               *
!   The polarization potential can be either:                          *
!   - modeled as                                                       *
!       VPOL(R) = -0.5 * ALPHA_D * R**2 / (R**3 + <R^3>)**2 -          *
!                  0.5 * ALPHA_Q * R**4 / (R**5 + <R^5>)**2            *
!                                                                      *
!   - read from a text 'vpol' file in format                           *
!       R1 VPOL(R1)                                                    *
!       R2 VPOL(R2)                                                    *
!         ...                                                          *
!                                                                      *
!   Call(s) to: [LIB92] RINT, INTERP                                   *
!                                                                      *
!   Paweł Syty                                                         *
!                                                          June 2024   *
!                                                                      *
!***********************************************************************

    USE vast_kind_param, ONLY: DOUBLE
    USE continuum_C
    USE def_C, ONLY: Z
    USE rint_I
    USE pote_C
    USE orb_C
    USE grid_C

    IMPLICIT none

    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: rpol
    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: vpol
    CHARACTER(LEN=256) :: vpol_file = "vpol"
    INTEGER            :: n_vpol, i
    LOGICAL            :: vpol_file_exists = .FALSE.

    IF (CO_POLARIZATION_FROM_FILE) THEN

        INQUIRE (FILE=vpol_file, EXIST=vpol_file_exists)

        ! Try to read potential from file
        IF (vpol_file_exists) THEN
            PRINT*, "Include polarization potential: read from '", TRIM(vpol_file), "' file."

            ! Test for number of datapoints of numerical potential
            n_vpol = 0
            OPEN(76, FILE=vpol_file, FORM="FORMATTED")
            DO
                READ (76, *, END=10)
                n_vpol = n_vpol + 1
            END DO
            10 CLOSE(76)

            ! Read data
            OPEN(UNIT=76, FILE=vpol_file, FORM="FORMATTED")

            IF (.NOT. ALLOCATED(rpol)) ALLOCATE(rpol(1:n_vpol))
            IF (.NOT. ALLOCATED(vpol)) ALLOCATE(vpol(1:n_vpol))
            DO i = 1, n_vpol
                READ(76, *) rpol(i), vpol(i)
            END DO

            CLOSE(76)

            DO i=1, N
                IF (R(i) <= rpol(1) .OR. R(i) >= rpol(n_vpol)) THEN
                    CO_VPOL(i) = 0.0D0
                ELSE
                    CALL INTERP(rpol, vpol, n_vpol, R(i), CO_VPOL(i), 0.1D0)
                END IF
            END DO

        ELSE

            STOP "*** ERROR *** No 'vpol` file found."

        END IF

    END IF

    IF (CO_POLARIZATION_AUTO) THEN

        ! Use model potential
        ! Take dipole polarizabilities CO_ALPHA_D from
        ! P. Schwerdtfeger and J.K. Nagle, Molecular Physics_ 117 9-12 (2019)
        SELECT CASE (INT(Z))

        CASE (1) ! H
            CO_ALPHA_D = 4.50711D0
        CASE (2) ! He
            CO_ALPHA_D = 1.38375D0
        CASE (3) ! Li
            CO_ALPHA_D = 164.1125D0
        CASE (4) ! Be
            CO_ALPHA_D = 37.74D0
        CASE (5) ! B
            CO_ALPHA_D = 20.5D0
        CASE (6) ! C
            CO_ALPHA_D = 11.3D0
        CASE (7) ! N
            CO_ALPHA_D = 7.4D0
        CASE (8) ! O
            CO_ALPHA_D = 5.3D0
        CASE (9) ! F
            CO_ALPHA_D = 3.74D0
        CASE (10) ! Neon
            CO_ALPHA_D = 2.6611D0
        CASE (11) ! Na
            CO_ALPHA_D = 162.7D0
        CASE (12) ! Mg
            CO_ALPHA_D = 71.2D0
        CASE (13) ! Al
            CO_ALPHA_D = 57.8D0
        CASE (14) ! Si
            CO_ALPHA_D = 37.3D0
        CASE (15) ! P
            CO_ALPHA_D = 25.0D0
        CASE (16) ! S
            CO_ALPHA_D = 19.4D0
        CASE (17) !
            CO_ALPHA_D = 14.6D0
        CASE (18) ! Argon
            CO_ALPHA_D = 11.083D0
        CASE (19) ! K
            CO_ALPHA_D = 289.7D0
        CASE (20) ! Ca
            CO_ALPHA_D = 160.8D0
        CASE (21) ! Sc
            CO_ALPHA_D = 97.0D0
        CASE (22) ! Ti
            CO_ALPHA_D = 100.0D0
        CASE (23) ! V
            CO_ALPHA_D = 87.0D0
        CASE (24) ! Cr
            CO_ALPHA_D = 83.0D0
        CASE (25) ! Mn
            CO_ALPHA_D = 68.0D0
        CASE (26) ! Fe
            CO_ALPHA_D = 62.0D0
        CASE (27) ! Co
            CO_ALPHA_D = 55.0D0
        CASE (28) ! Ni
            CO_ALPHA_D = 49.0D0
        CASE (29) ! Cu
            CO_ALPHA_D = 46.5D0
        CASE (30) ! Zinc
            CO_ALPHA_D = 38.67D0
        CASE (31) ! Ga
            CO_ALPHA_D = 50.0D0
        CASE (32) ! Ge
            CO_ALPHA_D = 40.0D0
        CASE (33) ! As
            CO_ALPHA_D = 30.0D0
        CASE (34) ! Se
            CO_ALPHA_D = 28.9D0
        CASE (35) ! Br
            CO_ALPHA_D = 21.0D0
        CASE (36) ! Krypton
            CO_ALPHA_D = 16.78D0
        CASE (37) ! Rb
            CO_ALPHA_D = 319.8D0
        CASE (38) ! Strontium
            CO_ALPHA_D = 197.2D0
        CASE (39) ! Y
            CO_ALPHA_D = 162.0D0
        CASE (40) ! Zr
            CO_ALPHA_D = 112.0D0
        CASE (41) ! Nb
            CO_ALPHA_D = 98.0D0
        CASE (42) ! Mo
            CO_ALPHA_D = 87.0D0
        CASE (43) ! Tc
            CO_ALPHA_D = 79.0D0
        CASE (44) ! Ru
            CO_ALPHA_D = 72.0D0
        CASE (45) ! Th
            CO_ALPHA_D = 66.0D0
        CASE (46) ! Pd
            CO_ALPHA_D = 26.14D0
        CASE (47) ! Ag
            CO_ALPHA_D = 55.0D0
        CASE (48) ! Cd
            CO_ALPHA_D = 46.0D0
        CASE (49) ! In
            CO_ALPHA_D = 65.0D0
        CASE (50) ! Sn
            CO_ALPHA_D = 53.0D0
        CASE (51) ! Sb
            CO_ALPHA_D = 43.0D0
        CASE (52) ! Te
            CO_ALPHA_D = 38.0D0
        CASE (53) ! I
            CO_ALPHA_D = 32.9D0
        CASE (54) ! Xenon
            CO_ALPHA_D = 27.32D0
        CASE (55) ! Cs
            CO_ALPHA_D = 400.9D0
        CASE (56) ! Barium
            CO_ALPHA_D = 272.0D0
        CASE (57) ! La
            CO_ALPHA_D = 215.0D0
        CASE (58) ! Ce
            CO_ALPHA_D = 205.0D0
        CASE (59) ! Pr
            CO_ALPHA_D = 216.0D0
        CASE (60) ! Nd
            CO_ALPHA_D = 208.0D0
        CASE (61) ! Pm
            CO_ALPHA_D = 200.0D0
        CASE (62) ! Sm
            CO_ALPHA_D = 192.0D0
        CASE (63) ! Eu
            CO_ALPHA_D = 184.0D0
        CASE (64) ! Gd
            CO_ALPHA_D = 158.0D0
        CASE (65) ! Tb
            CO_ALPHA_D = 170.0D0
        CASE (66) ! Dy
            CO_ALPHA_D = 163.0D0
        CASE (67) ! Ho
            CO_ALPHA_D = 156.0D0
        CASE (68) ! Er
            CO_ALPHA_D = 150.0D0
        CASE (69) ! Tm
            CO_ALPHA_D = 144.0D0
        CASE (70) ! Yb
            CO_ALPHA_D = 139.0D0
        CASE (71) ! Lu
            CO_ALPHA_D = 137.0D0
        CASE (72) ! Hf
            CO_ALPHA_D = 103.0D0
        CASE (73) ! Ta
            CO_ALPHA_D = 74.0D0
        CASE (74) ! W
            CO_ALPHA_D = 68.0D0
        CASE (75) ! Re
            CO_ALPHA_D = 62.0D0
        CASE (76) ! Os
            CO_ALPHA_D = 57.0D0
        CASE (77) ! Ir
            CO_ALPHA_D = 54.0D0
        CASE (78) ! Pt
            CO_ALPHA_D = 58.0D0
        CASE (79) ! Au
            CO_ALPHA_D = 36.0D0
        CASE (80) ! Mercury
            CO_ALPHA_D = 33.91D0
        CASE (81) ! Ti
            CO_ALPHA_D = 50.0D0
        CASE (82) ! Pb
            CO_ALPHA_D = 47.0D0
        CASE (83) ! Bi
            CO_ALPHA_D = 48.0D0
        CASE (84) ! Po
            CO_ALPHA_D = 44.0D0
        CASE (85) ! At
            CO_ALPHA_D = 42.0D0
        CASE (86) ! Radon
            CO_ALPHA_D = 35.0D0
        CASE (87) ! Fr
            CO_ALPHA_D = 317.8D0
        CASE (88) ! Ra
            CO_ALPHA_D = 246.0D0
        CASE (89) ! Ac
            CO_ALPHA_D = 203.0D0
        CASE (90) ! Th
            CO_ALPHA_D = 217.0D0
        CASE (91) ! Pa
            CO_ALPHA_D = 154.0D0
        CASE (92) ! U
            CO_ALPHA_D = 129.0D0
        CASE (93) ! Np
            CO_ALPHA_D = 151.0D0
        CASE (94) ! Pu
            CO_ALPHA_D = 132.0D0
        CASE (95) ! Am
            CO_ALPHA_D = 131.0D0
        CASE (96) ! Cm
            CO_ALPHA_D = 144.0D0
        CASE (97) ! Bk
            CO_ALPHA_D = 125.0D0
        CASE (98) ! Cf
            CO_ALPHA_D = 122.0D0
        CASE (99) ! Es
            CO_ALPHA_D = 118.0D0
        CASE (100) ! Fm
            CO_ALPHA_D = 113.0D0
        CASE (101) ! Md
            CO_ALPHA_D = 109.0D0
        CASE (102) ! No
            CO_ALPHA_D = 110.0D0
        CASE (103) ! Lr
            CO_ALPHA_D = 320.0D0
        CASE (104) ! Rf
            CO_ALPHA_D = 112.0D0
        CASE (105) ! Db
            CO_ALPHA_D = 42.0D0
        CASE (106) ! Sg
            CO_ALPHA_D = 40.0D0
        CASE (107) ! Bh
            CO_ALPHA_D = 38.0D0
        CASE (108) ! Hs
            CO_ALPHA_D = 36.0D0
        CASE (109) ! Mt
            CO_ALPHA_D = 34.0D0
        CASE (110) ! Ds
            CO_ALPHA_D = 32.0D0
        CASE (111) ! Rg
            CO_ALPHA_D = 32.0D0
        CASE (112) ! Cn
            CO_ALPHA_D = 28.0D0
        CASE (113) ! Nh
            CO_ALPHA_D = 29.0D0
        CASE (114) ! Fl
            CO_ALPHA_D = 31.0D0
        CASE (115) ! Mc
            CO_ALPHA_D = 71.0D0
        CASE (116) ! Lv
            ! No data for Lv, 116
            CO_ALPHA_D = 0.0
            CO_INCLUDE_POLARIZATION = .FALSE.
            PRINT*,"No dipole polarizability is known for Lv, polarization potential will not be included."
        CASE (117) ! Ts
            CO_ALPHA_D = 76.0D0
        CASE (118) ! Oganesson
            CO_ALPHA_D = 58.0D0
        CASE DEFAULT
            CO_ALPHA_D = 0.0
            CO_INCLUDE_POLARIZATION = .FALSE.
            PRINT*,"No dipole polarization is found for the element, polarization potential will not be included."
        END SELECT

        ! Take cut-off from bound calculations (as size of the outermost bound orbital)
        CO_R0_3 = RINT(CO_ORBITAL-1,CO_ORBITAL-1, 3)
        PRINT*, "  <R0^3> is taken from size of the ", NP(CO_ORBITAL-1), NH(CO_ORBITAL-1), "orbital"

    ! In auto mode, there are no quadrupole terms available, so setting
    ! CO_ALPHA_Q and CO_R0_5 is omitted (they defaults to zero)

    END IF

    IF (CO_POLARIZATION_AUTO .OR. CO_POLARIZATION_MANUAL) THEN
        PRINT*,"Include model polarization potential, alpha_d = ", CO_ALPHA_D, ", <R0^3> = ", CO_R0_3
        CO_VPOL(:N) = -0.5 * CO_ALPHA_D * R(:N)**2 / (R(:N)**3 + CO_R0_3)**2

        IF (CO_ALPHA_Q /= 0 .AND. CO_R0_5 /= 0) THEN
            PRINT*,"Include quadrupole term, alpha_q = ", &
            CO_ALPHA_Q, ", <R0^5> = ", CO_R0_5
            CO_VPOL(:N) = &
            CO_VPOL(:N) -0.5 * CO_ALPHA_Q * R(:N)**4 / (R(:N)**5 + CO_R0_5)**2
        END IF

    END IF

    ! Add polarizaton potential to the direct one
    YP(:N) = YP(:N) - CO_VPOL(:N) * R(:N)

END SUBROUTINE co_add_polarization_potential