!***********************************************************************
!                                                                      *
subroutine co_getmax(x,y,n,xmax,ymax)
!                                                                      *
!   Returns the location xmax and the corresponding value ymax         *
!   of a local a local maximum by interpolation.                       *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!           x    : (input) Table with the x-values to interpolate      *
!           y    : (input) Table with the y-values to interpoate       *
!           n    : (input) Size of the x and y tables                  *
!           xmax : (output) The X value of the local maximum           *
!           ymax : (output) The Y value of the local maximum           *
!                                                                      *
!   Call(s) to: [LIB92]: INTERP                                        *
!                                                                      *
!   Paweł Syty                                                         *
!                                           Last update: 25 Apr 2024   *
!                                                                      *
!***********************************************************************
    USE vast_kind_param, ONLY: DOUBLE

    IMPLICIT none

    INTEGER, INTENT(IN)                      :: n
    REAL(DOUBLE), DIMENSION(1:n), INTENT(IN) :: x, y
    REAL(DOUBLE), INTENT(OUT)                :: xmax, ymax
    REAL(DOUBLE) :: x1, x2, x3, dx, y1, y2, y3
    REAL(DOUBLE), PARAMETER :: accy = 1.0D-3

    dx = (x(n) - x(1)) / 1000
    x1 = x(1) - dx/2

    DO
        x1 = x1 + dx
        x2 = x1 + dx
        x3 = x2 + dx
        CALL INTERP(x,y,n,x1,y1,accy)
        CALL INTERP(x,y,n,x2,y2,accy)
        CALL INTERP(x,y,n,x3,y3,accy)
        IF (y1 < y2 .and. y2 > y3) then
            xmax = x2
            ymax = y2
            return
        ELSE IF (x1 < x(1)  .or.  x3 < x(1)  .OR. &
                 x1 > x(n)  .or.  x3 > x(n)) THEN
            xmax = 0
            ymax = 0 ! Fail
            RETURN
        END IF
    END DO

end subroutine co_getmax
