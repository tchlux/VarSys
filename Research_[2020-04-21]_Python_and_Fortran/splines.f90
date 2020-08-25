! KNOTS -> T
! BREAKPOINTS -> XI
! STATUS -> INFO


MODULE REAL_PRECISION  ! module for 64-bit real arithmetic
  INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

MODULE SPLINES
  USE REAL_PRECISION
  IMPLICIT NONE

  ! PRIVATE :: EVAL_BSPLINE
  PUBLIC

CONTAINS


SUBROUTINE PMQSI(X, Y, KNOTS, COEFF, INFO)
  ! Compute a piecewise monotone quintic spline interpolant for data
  ! in terms of a B-spline basis represented with knots and
  ! coefficients.
  ! 
  ! INPUT:
  !   X(1:Z) -- 
  !   Y(1:Z) -- 
  ! 
  ! OUTPUT:
  !   KNOTS(1:Z+6) -- The knots for the PMQSI B-spline basis.
  !   COEFF(1:Z) -- The coefficients for the PMQSI B-spline basis.
  !   INFO -- Integer representing the subroutine execution status.
  !   
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X, Y
  ! 
  ! REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: KNOTS, COEFF
  REAL(KIND=R8), INTENT(OUT), DIMENSION(3*SIZE(X)+6) :: KNOTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(3*SIZE(X)) :: COEFF
  ! ^^^^ THESE ABOVE LINES ARE TEMPORARY, SHOULD NOT BE IN FINAL CODE.
  !      THEY ONLY EXIST TO MAKE AUTOMATIC WRAPPING EASIER.
  INTEGER, INTENT(OUT) :: INFO
  ! 
  ! Local variables.
  REAL(KIND=R8), DIMENSION(SIZE(X),3) :: VALUES ! Spline derivative values.
  REAL(KIND=R8), DIMENSION(SIZE(X),2) :: IDEAL_VALUES ! Estimated derivatives.
  REAL(KIND=R8), DIMENSION(3) :: AS, BS
  REAL(KIND=R8) :: ACCURACY, STEP_SIZE, A, B
  LOGICAL, DIMENSION(SIZE(X)) :: FIXING, MODIFIED
  LOGICAL :: SEARCHING
  INTEGER :: I, J, Z

  Z = SIZE(X)
  ! Check the shape of incoming arrays.
  IF      (Z .LT. 3)                 THEN ; INFO = 1 ; RETURN
  ELSE IF (SIZE(Y) .NE. Z)           THEN ; INFO = 2 ; RETURN
  ELSE IF (SIZE(KNOTS) .LT. 3*Z + 6) THEN ; INFO = 4 ; RETURN
  ELSE IF (SIZE(COEFF) .LT. 3*Z)     THEN ; INFO = 5 ; RETURN
  END IF

  ! Verify that X are increasing.
  DO I = 1, Z-1
     IF (X(I+1) - X(I) .LE. SQRT(EPSILON(1.0_R8))) THEN
        INFO = 6 ; RETURN
     END IF
  END DO

  ! Copy the "Y" values into the "VALUES" array, initialize others to zero.
  VALUES(:,1) = Y(:)
  FIXING(1) = .FALSE.
  FIXING(Z) = .FALSE.
  ! Store extreme points in the "FIXING" array.
  DO I = 2, Z-1
     FIXING(I) = (((Y(I) - Y(I-1)) * (Y(I+1) - Y(I))) .LT. 0.0_R8)
  END DO

  ! Use local quadratic interpolants to estimate slopes and second
  ! derivatives at all points. Use zero-slope quadratic interpolant
  ! at extrema to estimate curvature.
  !  - store outputs in "VALUES"
  DO I = 1, Z
     ! If this is an extreme point, construct quadratic interpolants
     ! that have zero slope here and hit left/right neighbors.
     IF (FIXING(I)) THEN
        AS(1) = (Y(I-1) - Y(I)) / (X(I-1) - X(I))**2
        BS(1) = - 2 * X(I) * AS(1)
        AS(2) = (Y(I+1) - Y(I)) / (X(I+1) - X(I))**2
        BS(2) = - 2 * X(I) * AS(2)
        AS(3) = HUGE(1.0_R8)
        BS(3) = 0.0_R8
     ELSE
        ! If there are no left quadratics, then skip.
        IF (I .LE. 2) THEN ; AS(1) = HUGE(1.0_R8) ; BS(1) = 0.0_R8
        ! If there's an extreme point to the left, use it's right interpolant.
        ELSE IF (FIXING(I-1)) THEN
           AS(1) = (Y(I) - Y(I-1)) / (X(I) - X(I-1))**2
           BS(1) = - 2 * X(I-1) * AS(1)
        ! Otherwise use the standard quadratic on the left.
        ELSE ; CALL QUADRATIC(I-1, AS(1), BS(1), INFO)
        END IF
        IF ((I .GT. 1) .AND. (I .LT. Z)) THEN
           ! Construct the quadratic interpolant through this point and neighbors.
           CALL QUADRATIC(I, AS(2), BS(2), INFO)
        ELSE ; AS(2) = HUGE(1.0_R8) ; BS(2) = 0.0_R8
        END IF
        ! If there are no right quadratics, then skip.
        IF (I .GE. Z-1) THEN ; AS(3) = HUGE(1.0_R8) ; BS(3) = 0.0_R8
        ! If there's an extreme point to the right, use it's left interpolant.
        ELSE IF (FIXING(I+1)) THEN
           AS(3) = (Y(I) - Y(I+1)) / (X(I) - X(I+1))**2
           BS(3) = - 2 * X(I+1) * AS(3)
        ! Otherwise use the standard quadratic on the right.
        ELSE ; CALL QUADRATIC(I+1, AS(3), BS(3), INFO)
        END IF
     END IF
     ! Get the best quadratic for this index (uses AS and BS)
     ! PRINT *, ''
     ! PRINT *, '-------------------------------------------------------'
     ! PRINT *, I
     ! PRINT *, '  As:', AS(:)
     ! PRINT *, '  Bs:', BS(:)
     J = BEST_QUADRATIC(I)
     ! PRINT *, '  Picked:', J
     IF (J .GT. 0) THEN
        ! Set the first and second derivatives.
        VALUES(I,2) = 2 * AS(J) * X(I) + BS(J)
        VALUES(I,3) = 2 * AS(J)
        ! Make sure there are no "infinity" values.
        IF (ABS(VALUES(I,2)) .GT. HUGE(1.0_R8)) THEN
           IF (VALUES(I,2) .GT. 0.0_R8) THEN ; VALUES(I,2) = HUGE(1.0_R8)
           ELSE ; VALUES(I,2) = -HUGE(1.0_R8)
           END IF
        END IF
        IF (ABS(VALUES(I,3)) .GT. HUGE(1.0_R8)) THEN
           IF (VALUES(I,3) .GT. 0.0_R8) THEN ; VALUES(I,3) = HUGE(1.0_R8)
           ELSE ; VALUES(I,3) = -HUGE(1.0_R8)
           END IF
        END IF
     ELSE
        ! PRINT *, 'WARNING: Found no appropriate quadratic. Using 0s at ',I
        VALUES(I,2) = 0.0_R8
        VALUES(I,3) = 0.0_R8
     END IF
  END DO

  ! Store the initially estimated values as "ideal".
  IDEAL_VALUES(:,1:2) = VALUES(:,2:3)

  ! PRINT *, ''
  ! print *, "IDEAL_VALUES(:,1): ",IDEAL_VALUES(:,1)
  ! print *, "IDEAL_VALUES(:,2): ",IDEAL_VALUES(:,2)
  ! PRINT *, ''

  ! Identify which intervals are not monotone and need to be fixed.
  FIXING(:) = .FALSE.
  MODIFIED(:) = .FALSE.
  DO I = 1, Z-1
     IF (.NOT. IS_MONOTONE(I,I+1)) THEN
        FIXING(I) = .TRUE.
        FIXING(I+1) = .TRUE.
     END IF
  END DO
  ! PRINT *, ''
  ! PRINT *, 'VALUES:'
  ! PRINT *, '  ', VALUES(:,1)
  ! PRINT *, '  ', VALUES(:,2)
  ! PRINT *, '  ', VALUES(:,3)
  ! PRINT *, ''
  ! PRINT *, 'FIXING:', FIXING
  ! PRINT *, ''
  ! Initialize step size to 1.0 (will be halved at beginning of loop.
  STEP_SIZE = 1.0_R8
  ! Define the accuracy of the solution (will be an optional parameter).
  ACCURACY = SQRT(SQRT(EPSILON(1.0_R8)))
  SEARCHING = .TRUE.
  ! Loop until the accuracy is achieved and *all* intervals are monotone.
  DO WHILE (SEARCHING .OR. ANY(FIXING(:)))
     ! Compute the step size for this iteration.
     IF (SEARCHING) THEN
        STEP_SIZE = STEP_SIZE / 2.0_R8
        IF (STEP_SIZE .LT. ACCURACY) THEN
           SEARCHING = .FALSE.
           STEP_SIZE = ACCURACY
        END IF
     END IF
     ! PRINT *, 'step size', STEP_SIZE
     ! Cycle through all intervals and make modifications where appropriate.
     DO I = 1, Z
        IF (FIXING(I)) THEN
           ! PRINT *, '     fixing', I, VALUES(I,2:3)
           FIXING(I) = .FALSE.
           MODIFIED(I) = .TRUE.
           ! Shrink those values that are causing nonmonotonicity.
           VALUES(I,2) = VALUES(I,2) - STEP_SIZE * IDEAL_VALUES(I,1)
           VALUES(I,3) = VALUES(I,3) - STEP_SIZE * IDEAL_VALUES(I,2)
           ! Make sure the first derivative does not pass zero.
           IF ((IDEAL_VALUES(I,1) .LT. 0.0_R8) .AND. &
                (VALUES(I,2) .GT. 0.0_R8)) THEN
              VALUES(I,2) = 0.0_R8
           ELSE IF ((IDEAL_VALUES(I,1) .GT. 0.0_R8) .AND. &
                (VALUES(I,2) .LT. 0.0_R8)) THEN
              VALUES(I,2) = 0.0_R8
           END IF
           ! Make sure the second derivative does not pass zero.
           IF ((IDEAL_VALUES(I,2) .LT. 0.0_R8) .AND. &
                (VALUES(I,3) .GT. 0.0_R8)) THEN
              VALUES(I,3) = 0.0_R8
           ELSE IF ((IDEAL_VALUES(I,2) .GT. 0.0_R8) .AND. &
                (VALUES(I,3) .LT. 0.0_R8)) THEN
              VALUES(I,3) = 0.0_R8
           END IF
        ELSE IF (MODIFIED(I) .AND. SEARCHING) THEN
           ! Grow those values that have been modified previously.
           VALUES(I,2) = VALUES(I,2) + STEP_SIZE * IDEAL_VALUES(I,1)
           VALUES(I,3) = VALUES(I,3) + STEP_SIZE * IDEAL_VALUES(I,2)
           ! Make sure the first derivative does not pass original values.
           IF ((IDEAL_VALUES(I,1) .LT. 0.0_R8) .AND. &
                (VALUES(I,2) .LT. IDEAL_VALUES(I,1))) THEN
              VALUES(I,2) = IDEAL_VALUES(I,1)
           ELSE IF ((IDEAL_VALUES(I,1) .GT. 0.0_R8) .AND. &
                (VALUES(I,2) .GT. IDEAL_VALUES(I,1))) THEN
              VALUES(I,2) = IDEAL_VALUES(I,1)
           END IF
           ! Make sure the second derivative does not pass original values.
           IF ((IDEAL_VALUES(I,2) .LT. 0.0_R8) .AND. &
                (VALUES(I,3) .LT. IDEAL_VALUES(I,2))) THEN
              VALUES(I,3) = IDEAL_VALUES(I,2)
           ELSE IF ((IDEAL_VALUES(I,2) .GT. 0.0_R8) .AND. &
                (VALUES(I,3) .GT. IDEAL_VALUES(I,2))) THEN
              VALUES(I,3) = IDEAL_VALUES(I,2)
           END IF
        END IF
     END DO
     ! Go through and identify which values are associated with
     ! nonmonotone intervals after all the updates.
     DO I = 1, Z-1
        IF (.NOT. IS_MONOTONE(I,I+1)) THEN
           FIXING(I) = .TRUE.
           FIXING(I+1) = .TRUE.
        END IF
     END DO
  END DO
  ! PRINT *, 'Fitting spline..'
  ! Use "FIT_SPLINE" to fit the final interpolating spline.
  CALL FIT_SPLINE(X, VALUES, KNOTS, COEFF, INFO)

CONTAINS
  FUNCTION BEST_QUADRATIC(I) RESULT(LOC)
    INTEGER, INTENT(IN) :: I
    INTEGER :: LOC, ORDER(3)
    REAL(KIND=R8) :: SLOPE, DIRECTION
    ! Construct the sorted-order of preference for chosen quadratics.
    IF (ABS(AS(1)) .LE. ABS(AS(2))) THEN
       IF (ABS(AS(2)) .LE. ABS(AS(3))) THEN
          ORDER(:) = (/ 1, 2, 3 /)
       ELSE IF (ABS(AS(1)) .LE. ABS(AS(3))) THEN
          ORDER(:) = (/ 1, 3, 2 /)
       ELSE
          ORDER(:) = (/ 3, 1, 2 /)
       END IF
    ELSE
       IF (ABS(AS(2)) .GT. ABS(AS(3))) THEN
          ORDER(:) = (/ 3, 2, 1 /)
       ELSE IF (ABS(AS(1)) .GT. ABS(AS(3))) THEN
          ORDER(:) = (/ 2, 3, 1 /)
       ELSE
          ORDER(:) = (/ 2, 1, 3 /)
       END IF
    END IF
    ! Determine the direction of change at the point "I".
    IF (I .GE. Z) THEN
       IF (I .LE. 1) THEN ; DIRECTION = 0.0_R8
       ELSE IF (Y(I-1) .LE. Y(I)) THEN ; DIRECTION = 1.0_R8
       ELSE ; DIRECTION = -1.0_R8
       END IF
    ELSE IF (Y(I) .LE. Y(I+1)) THEN
       IF (I .LE. 1) THEN ; DIRECTION = 1.0_R8
       ELSE IF (Y(I-1) .LE. Y(I)) THEN ; DIRECTION = 1.0_R8
       ELSE ; DIRECTION = 0.0_R8
       END IF
    ELSE
       IF (I .LE. 1) THEN ; DIRECTION = -1.0_R8
       ELSE IF (Y(I-1) .GT. Y(I)) THEN ; DIRECTION = -1.0_R8
       ELSE ; DIRECTION = 0.0_R8
       END IF
    END IF
    ! PRINT *, '  DIRECTION:', DIRECTION
    ! PRINT *, '  ORDER:    ', ORDER
    ! Pick the quadratic that has the least curvature and an agreeing slope.
    LOC = ORDER(1)
    SLOPE = 2 * AS(LOC) * X(I) + BS(LOC)
    ! PRINT *, '    SLOPE 1:', SLOPE
    IF (SLOPE * DIRECTION .LT. 0.0_R8) THEN
       LOC = ORDER(2)
       SLOPE = 2 * AS(LOC) * X(I) + BS(LOC)
       ! PRINT *, '    SLOPE 2:', SLOPE
       IF (SLOPE * DIRECTION .LT. 0.0_R8) THEN
          LOC = ORDER(3)
          SLOPE = 2 * AS(LOC) * X(I) + BS(LOC)
          ! PRINT *, '    SLOPE 3:', SLOPE
          IF (SLOPE * DIRECTION .LT. 0.0_R8) THEN
             LOC = -1
          END IF
       END IF
    END IF
  END FUNCTION BEST_QUADRATIC

  SUBROUTINE QUADRATIC(I2, A, B, INFO)
    ! Given an index "I", compute the quadratic interpolant of points
    ! X(I-1:I+1), Y(I-1:I+1). Return the slope and curvature.
    ! 
    ! INPUT:
    !   I2 -- Integer index of the middle point defining the quadratic.
    ! 
    ! OUTPUT:
    !   A -- Real valued coefficient on x^2 for interpolating quadratic.
    !   B -- Real valued coefficient on x for interpolating quadratic.
    !   INFO -- Integer value denoting execution info, on exit:
    !     0     Successful execution.
    !     1     Points X(I2-1), X(I2), X(I2+1) are too close together
    !           to compute good estimates of slope and curvature.
    INTEGER, INTENT(IN) :: I2
    REAL(KIND=R8), INTENT(OUT) :: A, & ! Curvature of the quadratic.
         B ! Slope of the quadratic.
    INTEGER, INTENT(OUT) :: INFO ! 
    ! Local variables.
    REAL(KIND=R8) :: D, & ! Denominator for computing "A" and "B".
         C1, C2, C3 ! Intermediate terms for computation.
    INTEGER :: I1, I3
    I1 = I2-1
    I3 = I2+1
    IF ((Y(I1) .EQ. Y(I2)) .OR. (Y(I2) .EQ. Y(I3))) THEN
       A = 0.0_R8
       B = 0.0_R8
    END IF
    ! Compute the shared denominator and check for numerical stability.
    D = (X(I1) - X(I2)) * (X(I1) - X(I3)) * (X(I2) - X(I3))
    IF (ABS(D) .LT. SQRT(EPSILON(1.0_R8))) THEN ; INFO = 1 ; RETURN ; END IF
    ! Compute coefficients A and B in "Ax^2 + Bx + C" quadratic interpolant.
    C1 = X(I1) * (Y(I3) - Y(I2))
    C2 = X(I2) * (Y(I1) - Y(I3))
    C3 = X(I3) * (Y(I2) - Y(I1))
    A = (C1 + C2 + C3) / D
    B = - (X(I1)*C1 + X(I2)*C2 + X(I3)*C3) / D
  END SUBROUTINE QUADRATIC


  FUNCTION IS_MONOTONE(I1, I2)
    ! Given an interval and function values (value, first, and seocond
    ! derivative), make the quintic over this interval monotone by
    ! modifying the first and second derivatives according to quintic
    ! monotonicity theory. Return the change applied to first and second
    ! derivative values.
    INTEGER, INTENT(IN) :: I1, I2
    LOGICAL IS_MONOTONE
    ! Local variables.
    REAL(KIND=R8) :: A, A0, A1, B, BOUND, B0, B1, D, DDX0, DDX1, &
         DX0, DX1, G, G0, G1, SIGN, TAU, U0, U1, X0, X1
    U0 = X(I1)
    U1 = X(I2)
    X0 = VALUES(I1,1)
    X1 = VALUES(I2,1)
    DX0 = VALUES(I1,2)
    DX1 = VALUES(I2,2)
    DDX0 = VALUES(I1,3)
    DDX1 = VALUES(I2,3)
    ! Always consider the monotone increasing case.
    IF (X1 .LT. X0) THEN
       X0 = -X0     ; X1 = -X1
       DX0 = -DX0   ; DX1 = -DX1
       DDX0 = -DDX0 ; DDX1 = -DDX1
    END IF
    ! Make sure the slopes point in the correct direction.
    IF (DX0 .LT. 0.0_R8) THEN ; IS_MONOTONE = .FALSE. ; RETURN ; END IF
    IF (DX1 .LT. 0.0_R8) THEN ; IS_MONOTONE = .FALSE. ; RETURN ; END IF
    ! Compute A and B.
    A = (U1 - U0) * DX0 / (X1 - X0)
    B = (U1 - U0) * DX1 / (X1 - X0)
    ! Simplified cubic monotone case.
    IF (A*B .LE. 0) THEN
       A = (4*DX1 + DDX1*(U0-U1)) * (U0-U1) / (X0-X1)
       B = 30 + ((-24*(DX0+DX1) + 3*(DDX0-DDX1)*(U0-U1))*(U0-U1)) / (2 * (X0-X1))
       G = ((U0-U1) * (4*DX0 + DDX0*(U1-U0))) / (X0-X1)
       D = DX0 * (U0-U1) / (X0-X1)
       IS_MONOTONE = (A .GE. 0.0_R8) .AND. (D .GE. 0.0_R8) .AND. &
            (B .GE. A - (4.0_R8 * A * D)**(0.5_R8)) .AND. &
            (G .GE. D - (4.0_R8 * A * D)**(0.5_R8))
       RETURN
    END IF
    ! Full quintic monotone case.
    TAU = 24 + 2*(A*B)**(0.5_R8) - 3*(A+B)
    IF (TAU .LT. 0.0_R8) THEN ; IS_MONOTONE = .FALSE. ; RETURN ; END IF
    ! Compute alpha, gamma ,beta from theorems, determine monotonicity.
    A0 = 4 * (B**(1/4) / A**(1/4))
    A1 = ((U0-U1) / DX1) * B**(1/4) / A**(1/4)
    G0 = 4 * (DX0 / DX1) * (B**(3/4) / A**(3/4))
    G1 = ((U1-U0) / DX1) * (B**(3/4) / A**(3/4))
    B0 = (12 * (DX0+DX1) * (U1-U0) + 30 * (X0-X1)) / ((X0-X1) * A**(1/2) * B**(1/2))
    B1 = (3 * (U0-U1)**2) / (2 * (X0-X1) * A**(1/2) * B**(1/2))
    ! Compute the monotonicity condition.
    A = A0 + A1 * DDX1
    G = G0 + G1 * DDX0
    B = B0  + B1  * (DDX0 - DDX1)
    IF (B .LE. 6.0_R8) THEN ; BOUND = - (B + 2) / 2
    ELSE ;               BOUND = -2 * (B - 2)**(1/2)
    END IF
    IS_MONOTONE = (A > BOUND) .AND. (G > BOUND)
  END FUNCTION IS_MONOTONE

END SUBROUTINE PMQSI


SUBROUTINE FIT_SPLINE(BREAKPOINTS, VALUES, KNOTS, COEFF, INFO)
  ! Subroutine for computing a linear combination of B-splines that
  ! interpolates the given function value (and function derivatives)
  ! at the given breakpoints.
  ! 
  ! INPUT:
  !   BREAKPOINTS(1:NB) -- The increasing real-valued locations of the
  !                        breakpoints for the interpolating spline.
  !   VALUES(1:NB,1:NCC) -- VALUES(I,J) contains the (J-1)st derivative
  !                         at BREAKPOINTS(I) to be interpolated.
  ! 
  ! OUTPUT:
  !   KNOTS(1:NB*NCC+2*NCC) -- The nondecreasing real-valued locations
  !                            of the knots for the B-spline basis.
  !   COEFF(1:NB*NCC) -- The coefficients for the B-splines
  !                      that define the interpolating spline.
  !   INFO -- Integer representing the subroutine execution status:
  !       0    Successful execution.
  !       1    SIZE(BREAKPOINTS) is less than 1.
  !       2    SIZE(VALUES) is less then 1.
  !       3    SIZE(VALUES,1) does not equal SIZE(BREAKPOINTS).
  !       4    SIZE(KNOTS) too small, should be at least NB*NCC + 2*NCC.
  !       5    SIZE(COEFF) too small, should be at least NB*NCC.
  !       6    Elements of BREAKPOINTS are not strictly increasing.
  !       7    The computed spline does not match the provided VALUES
  !            and this fit should be disregarded. This arises when
  !            the scaling of function values and derivative values
  !            causes the resulting linear system to have a
  !            prohibitively large condition number.
  !     >10    10 plus the info flag as returned by DGBSV from LAPACK.
  ! 
  !   
  ! DESCRIPTION:
  ! 
  !   This subroutine uses a B-spline basis to interpolate given
  !   function values (and derivative values) at unique breakpoints.
  !   The interpolating spline is returned in terms of KNOTS and
  !   COEFF that define the underlying B-splines and the
  !   corresponding linear combination that interpolates given data
  !   respectively. This function uses the subroutine EVAL_BSPLINE to
  !   evaluate the B-splines at all knots and the LAPACK routine DGBSV
  !   to compute the coefficients of all component B-splines. The
  !   difference between the provided function (and derivative) values
  !   and the actual values produced by this code can vary depending
  !   on the spacing of the knots and the magnitudes of the values
  !   provided. When the condition number of the intermediate linear
  !   system grows prohibitively large, this routine may fail to
  !   produce a correct set of coefficients and return INFO code 7.
  !   For very high levels of continuity, or when this routine fails,
  !   a Newton form of polynomial representation should be used
  !   instead.
  ! 
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:)   :: BREAKPOINTS
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: VALUES
  ! 
  ! REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: KNOTS, COEFF
  REAL(KIND=R8), INTENT(OUT), DIMENSION(SIZE(VALUES)+2*SIZE(VALUES,2)) :: KNOTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(SIZE(VALUES)) :: COEFF
  ! ^^^^ THE ABOVE 2 LINES ARE TEMPORARY, SHOULD NOT BE IN FINAL CODE.
  !      THEY ONLY EXIST TO MAKE AUTOMATIC WRAPPING EASIER.
  ! 
  INTEGER, INTENT(OUT) :: INFO
  ! Local variables.
  INTEGER, DIMENSION(SIZE(COEFF)) :: IPIV
  REAL(KIND=R8), DIMENSION(1 + 3*(2*SIZE(VALUES,2)-1), SIZE(VALUES)) :: AB
  REAL(KIND=R8) :: MAX_ERROR
  INTEGER :: DEGREE, & ! Degree of B-spline.
       DERIV, & ! Derivative loop index.       
       I, I1, I2, J, J1, J2, & ! Miscellaneous loop indices.
       K, &    ! Order of B-splines = 2*NCC.           
       NB, &   ! Number of breakpoints.
       NCC, &  ! Number of continuity conditions.
       NK, &   ! Number of knots = NSPL + 2*NCC.
       NSPL    ! Dimension of spline space = number of B-spline 
               ! coefficients = NB * NCC.
  ! LAPACK subroutine for solving banded linear systems.
  EXTERNAL :: DGBSV

  ! Define some local variables for notational convenience.
  NB = SIZE(BREAKPOINTS)
  NCC = SIZE(VALUES,2)
  NSPL = SIZE(VALUES)
  NK = NSPL + 2*NCC
  K = 2*NCC
  DEGREE = K - 1
  INFO = 0

  ! Check the shape of incoming arrays.
  IF      (NB .LT. 1)              THEN ; INFO = 1 ; RETURN
  ELSE IF (NSPL .LT. 1)            THEN ; INFO = 2 ; RETURN
  ELSE IF (SIZE(VALUES,1) .NE. NB) THEN ; INFO = 3 ; RETURN
  ELSE IF (SIZE(KNOTS) .LT. NK)    THEN ; INFO = 4 ; RETURN
  ELSE IF (SIZE(COEFF) .LT. NSPL)  THEN ; INFO = 5 ; RETURN
  END IF
  ! Verify that BREAKPOINTS are increasing.
  DO I = 1, NB - 1
     IF (BREAKPOINTS(I) .GE. BREAKPOINTS(I+1)) THEN
        INFO = 6 ; RETURN
     END IF
  END DO

  ! Copy over the knots that will define the B-spline representation.
  ! Each knot will be repeataed NCC times to maintain the necessary
  ! level of continuity for this spline.
  KNOTS(1:K) = BREAKPOINTS(1)
  DO I = 2, NB-1
     KNOTS(I*NCC+1 : (I+1)*NCC) = BREAKPOINTS(I)
  END DO
  ! Assign the last knot to exist a small step outside the supported
  ! interval to ensure the B-spline basis functions are nonzero at the
  ! rightmost breakpoint.
  KNOTS(NK-DEGREE:NK) = BREAKPOINTS(NB) * (1.0_R8 + &
       SQRT(EPSILON(BREAKPOINTS(NB))))

  ! The next block of code evaluates each B-spline and it's
  ! derivatives at all breakpoints. The first and last elements of
  ! BREAKPOINTS will be repeated K times and each internal breakpoint
  ! will be repeated NCC times. As a result, the support of each
  ! B-spline spans at most three breakpoints. The coefficients for the
  ! B-spline basis are determined by a linear solve of a matrix with
  ! NB columns (for each B-spline) and NB*NCC rows (for each value of
  ! the interpolating spline).
  ! 
  ! For example, a C^1 interpolating spline over three breakpoints
  ! will match function value and first derivative at each breakpoint
  ! requiring six fourth order (third degree) B-splines each composed
  ! from five knots. Below, the six B-splines are numbered (first
  ! number, columns) and may be nonzero at the three breakpoints
  ! (middle letter, rows) for each function value (odd rows, terms end
  ! with 0) and first derivative (even rows, terms end with 1). The
  ! linear system will look like:
  ! 
  !       B-SPLINE VALUES AT BREAKPOINTS     SPLINE           VALUES
  !        1st  2nd  3rd  4th  5th  6th    COEFICIENTS
  !      _                              _     _   _           _    _ 
  !     |                                |   |     |         |      |
  !   B |  1a0  2a0  3a0  4a0            |   |  1  |         |  a0  |
  !   R |  1a1  2a1  3a1  4a1            |   |  2  |         |  a1  |
  !   E |  1b0  2b0  3b0  4b0  5b0  6b0  |   |  3  |   ===   |  b0  |
  !   A |  1b1  2b1  3b1  4b1  5b1  6b0  |   |  4  |   ===   |  b1  |
  !   K |            3c0  4c0  5c0  6c0  |   |  5  |         |  c0  |
  !   S |            3c1  4c1  5c1  6c0  |   |  6  |         |  c1  |
  !     |_                              _|   |_   _|         |_    _|
  !   
  ! Notice this matrix is banded with lower / upper bandwidths equal
  ! to (one less than the maximum number of breakpoints for which a
  ! spline takes on a nonzero value) times (the number of continuity
  ! conditions) minus (one). In general KL = KU = DEGREE = K - 1.
  
  ! Initialize all values in AB to zero.
  AB(:,:) = 0_R8
  ! Evaluate all B-splines at all breakpoints (walking through columns).
  DO I = 1, NSPL
     ! Compute indices of the first and last knot for the current B-spline.
     J  = I + K ! Last knot.
     ! Compute the row indices in "A" that would be accessed.
     I1 = ((I-1)/NCC - 1) * NCC + 1 ! First row.
     I2  = I1 + 3*NCC - 1           ! Last row.
     ! Only two breakpoints will be covered for the first NCC
     ! B-splines and the last NCC B-splines.
     IF (I .LE. NCC)     I1 = I1 + NCC
     IF (I+NCC .GT. NSPL)  I2 = I2 - NCC
     ! Compute the indices of the involved breakpoints.
     J1 = I1 / NCC + 1 ! First breakpoint.
     J2 = I2  / NCC    ! Last breakpoint.
     ! Convert the "i,j" indices in "A" to the banded storage scheme.
     ! The mapping is looks like   AB[KL+KU+1+i-j,j] = A[i,j]
     I1 = 2*DEGREE+1 + I1 - I
     I2 = 2*DEGREE+1 + I2  - I
     ! Evaluate this B-spline, computing function value and derivatives.
     DO DERIV = 0, NCC-1
        ! Place the evaluations into a block of a column in AB, shift
        ! according to which derivative is being evaluated and use a
        ! stride determined by the number of continuity conditions.
        AB(I1+DERIV:I2:NCC,I) = BREAKPOINTS(J1:J2)
        CALL EVAL_BSPLINE(KNOTS(I:J), AB(I1+DERIV:I2:NCC,I), INFO, D=DERIV)
        ! ^ Correct usage is inherently enforced, only extrapolation
        !   warnings will be produced by this call (INFO=1). These
        !   extrapolation warnings are expected because underlying
        !   B-splines may not support the full interval.
     END DO
  END DO
  ! Copy the VALUES into the COEFF (output) variable.
  DO I = 1, NCC
     COEFF(I:NSPL:NCC) = VALUES(:,I)
  END DO
  ! Call the LAPACK subroutine to solve the banded linear system.
  CALL DGBSV(NSPL, DEGREE, DEGREE, 1, AB, SIZE(AB,1), IPIV, &
       COEFF, NSPL, INFO)
  ! Check for errors in the execution of DGBSV, (this should not happen).
  IF (INFO .NE. 0) THEN ; INFO = INFO + 10 ; RETURN ; END IF
  ! Check to see if the linear system was correctly solved by looking at
  ! the difference between prouduced B-spline values and provided values.
  MAX_ERROR = SQRT(SQRT(EPSILON(1.0_R8)))
  DO DERIV = 0, NCC-1
     ! Reuse the first row of AB as scratch space (the first column
     ! might not be large enough, but the first row certainly is).
     AB(1,1:NB) = BREAKPOINTS(:)
     ! Evaluate this spline at all breakpoints. Correct usage is
     ! enforced here, so it is expected that INFO=0 always.
     CALL EVAL_SPLINE(KNOTS, COEFF, AB(1,1:NB), INFO, D=DERIV)
     ! Check the precision of the reproduced values.
     ! Return an error if the precision is too low.
     IF (MAXVAL(ABS((AB(1,1:NB) - VALUES(:,DERIV+1)) &
          / (1.0_R8 + ABS(VALUES(:,DERIV+1))))) .GT. MAX_ERROR) THEN
        INFO = 7
        RETURN
     END IF
  END DO
END SUBROUTINE FIT_SPLINE


SUBROUTINE EVAL_SPLINE(KNOTS, COEFF, XY, INFO, D)
  ! Evaluate a spline construced with FIT_SPLINE. Similar interface
  ! to EVAL_BSPLINE. Evaluate D derivative at all XY, result in XY.
  ! 
  ! INPUT:
  !   KNOTS(1:N+2*NCC) -- The nondecreasing real-valued locations of
  !                       the knots for the underlying B-splines, where
  !                       the spline has NCC-1 continuous derivatives.
  !   COEFF(1:N) -- The coefficients assigned to each B-spline
  !                 that underpins this interpolating spline.
  ! 
  ! INPUT / OUTPUT:
  !   XY(1:Z) -- The locations at which the spline is evaluated on
  !              input, on output holds the value of the spline with
  !              KNOTS and COEFF evaluated at the given locations.
  ! 
  ! OUTPUT:
  !   INFO -- Integer representing subroutine exeuction status.
  !       0    Successful execution.
  !       1    Extrapolation warning, some X are outside of spline support.
  !       2    KNOTS contains at least one decreasing interval.
  !       3    KNOTS has size less than or equal to 1.
  !       4    KNOTS has an empty interior (KNOTS(1) = KNOTS(N+2*C)).
  !       5    Invalid COEFFICEINTS, size smaller than or equal to KNOTS.
  ! 
  ! OPTIONAL INPUT:
  !   D [= 0]  --  The derivative to take of the evaluated spline.
  !                When negative, this subroutine integrates the spline.
  !                The higher integrals of this spline are capped at
  !                the rightmost knot, using constant-valued extrapolation.
  ! 
  ! 
  ! DESCRIPTION:
  ! 
  !    This subroutine serves as a convenient wrapper to the
  !    underlying calls to EVAL_BSPLINE that need to be made to
  !    evaluate the full spline. Internally this uses a matrix-vector
  !    multiplication of the B-spline evaluations with the assigned
  !    coefficients. This requires O(Z*N) memory, meaning single XY
  !    points should be evaluated at a time when memory-constrained.
  ! 
  REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: KNOTS, COEFF
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
  INTEGER, INTENT(IN), OPTIONAL :: D
  INTEGER, INTENT(OUT) :: INFO
  ! Local variables.
  INTEGER :: DERIV, I, ORDER
  REAL(KIND=R8), DIMENSION(SIZE(XY), SIZE(COEFF)) :: VALUES

  ! Check for size-related errors and for an empty knot interior.
  IF (SIZE(KNOTS) .LE. 1)                    THEN ; INFO = 3 ; RETURN
  ELSE IF (KNOTS(1) .EQ. KNOTS(SIZE(KNOTS))) THEN ; INFO = 4 ; RETURN
  ELSE IF (SIZE(KNOTS) .LE. SIZE(COEFF))     THEN ; INFO = 5 ; RETURN
  END IF
  ! Check for valid (nondecreasing) knot sequence.
  DO I = 1, SIZE(KNOTS)-1
     IF (KNOTS(I) .GT. KNOTS(I+1)) THEN ; INFO = 2 ; RETURN ; END IF
  END DO

  ! Compute the ORDER (number of knots minus one) for each B-spline.
  ORDER = SIZE(KNOTS) - SIZE(COEFF)
  ! Assign the local value of the optional derivative "D" argument.
  set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
  ELSE ; DERIV = 0
  END IF set_derivative

  ! Evaluate all splines at all the X positions.
  DO I = 1, SIZE(COEFF)
     IF (KNOTS(I) .EQ. KNOTS(I+ORDER)) CYCLE
     ! ^ If this constituent B-spline has no support, skip it.
     VALUES(:,I) = XY(:)
     CALL EVAL_BSPLINE(KNOTS(I:I+ORDER), VALUES(:,I), &
          INFO, D=DERIV)
     ! ^ Correct usage is inherently enforced, only extrapolation
     !   warnings will be produced by this call. These
     !   extrapolation warnings are expected because underlying
     !   B-splines may not support the full interval.
  END DO
  ! Set the EXTRAPOLATION status flag.
  IF ((MINVAL(XY(:)) .LT. KNOTS(1)) .OR. &
       (MAXVAL(XY(:)) .GE. KNOTS(SIZE(KNOTS)))) THEN ; INFO = 1
  ELSE ; INFO = 0
  END IF
  ! Store the values into Y as the weighted sums of B-spline evaluations.
  XY(:) = MATMUL(VALUES(:,:), COEFF(:))
END SUBROUTINE EVAL_SPLINE


SUBROUTINE EVAL_BSPLINE(KNOTS, XY, INFO, D)
  ! Subroutine for evaluating a B-spline with provided knot sequence.
  ! 
  ! INPUT:
  !   KNOTS(1:N) -- The nondecreasing sequence of knots for the B-spline.
  ! 
  ! INPUT / OUTPUT:
  !   XY(1:Z) -- The locations at which the B-spline is evaluated on
  !              input, on output holds the value of the B-spline with
  !              prescribed knots evaluated at the given X locations.
  ! 
  ! OUTPUT:
  !   INFO -- Execution status of this subroutine on exit.
  !       0    Successful execution.
  !       1    Extrapolation warning, some points were outside all knots.
  !       2    Invalid knot sequence (not entirely nondecreasing).
  !       3    Invalid size for KNOTS (less than or equal to 1).
  ! 
  ! OPTIONAL INPUT:
  !   D [= 0]  --  The derivative to take of the evaluated B-spline.
  !                When negative, this subroutine integrates the B-spline.
  ! 
  ! DESCRIPTION:
  ! 
  !    This function uses the recurrence relation defining a B-spline:
  ! 
  !      B_{K,1}(X)   =   1     if KNOTS(K) <= X < KNOTS(K+1),
  !                       0     otherwise,
  ! 
  !    where K is the knot index, I = 2, ..., N-MAX(D,0)-1, and
  ! 
  !                               X - KNOTS(K)                      
  !      B_{K,I}(X) =      ------------------------- B_{K,I-1}(X)   
  !                         KNOTS(K+I-1) - KNOTS(K)                 
  !                                                                   
  !                             KNOTS(K+I) - X                    
  !                     +  ------------------------- B_{K+1,I-1}(X).
  !                         KNOTS(K+I) - KNOTS(K+1)                 
  ! 
  !    All of the intermediate steps (I) are stored in a single block
  !    of memory that is reused for each step.
  ! 
  !    The computation of the integral of the B-spline proceeds from
  !    the above formula one integration step at a time by adding a
  !    duplicate of the last knot, raising the order of all
  !    intermediate B-splines, summing their values, and rescaling
  !    each sum by the width of the supported interval divided by the
  !    degree plus the integration coefficient.
  ! 
  !    For the computation of the derivative of the B-spline, the
  !    continuation of the standard recurrence relation is used that
  !    builds from I = N-D, ..., N-1 as
  ! 
  !                           (I-1) B_{K,I-1}(X)     
  !      B_{K,I}(X) =      -------------------------
  !                         KNOTS(K+I-1) - KNOTS(K)  
  !                                                   
  !                           (I-1) B_{K+1,I-1}(X)    
  !                     -  -------------------------.
  !                         KNOTS(K+I) - KNOTS(K+1) 
  ! 
  !     The final B-spline is right continuous and has support over
  !     the interval [KNOTS(1), KNOTS(N)).
  ! 
  REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: KNOTS
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
  INTEGER, INTENT(OUT) :: INFO
  INTEGER, INTENT(IN), OPTIONAL :: D
  ! Local variables.
  REAL(KIND=R8), DIMENSION(SIZE(XY), SIZE(KNOTS)) :: VALUES
  INTEGER :: K, N, DERIV, ORDER, I
  REAL(KIND=R8) :: LEFT, RIGHT, LAST_KNOT
  ! Assign the local value of the optional derivative "D" argument.
  set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
  ELSE ; DERIV = 0
  END IF set_derivative
  ! Collect local variables for notational convenience.
  N = SIZE(KNOTS)
  ORDER = N - 1
  LAST_KNOT = KNOTS(N)
  INFO = 0
  ! Check for valid knot sequence.
  IF (N .LE. 1) THEN ; INFO = 3 ; RETURN ; END IF
  DO K = 1, N-1
     IF (KNOTS(K) .GT. KNOTS(K+1)) THEN ; INFO = 2 ; RETURN ; END IF
  END DO
  ! Check for extrapolation, set status if it is happening, but continue.
  IF ((MINVAL(XY(:)) .LT. KNOTS(1)) .OR. (MAXVAL(XY(:)) .GE. LAST_KNOT)) &
       INFO = 1
  ! If this is a large enough derivative, we know it is zero everywhere.
  IF (DERIV+1 .GE. N) THEN ; XY(:) = 0.0_R8 ; RETURN
  ! ---------------- Performing standard evaluation ------------------
  ! This is a standard B-spline with multiple unique knots, right continuous.
  ELSE
     ! Initialize all values to 0.
     VALUES(:,:) = 0.0_R8
     ! Assign the first value for each knot index.
     first_b_spline : DO K = 1, ORDER
        IF (KNOTS(K) .EQ. KNOTS(K+1)) CYCLE
        ! Compute all right continuous order 1 B-spline values.
        WHERE ( (KNOTS(K) .LE. XY(:)) .AND. (XY(:) .LT. KNOTS(K+1)) )
           VALUES(:,K) = 1.0_R8
        END WHERE
     END DO first_b_spline
  END IF

  ! Compute the remainder of B-spline by building up from the first.
  ! Omit the final steps of this computation for derivatives.
  compute_spline : DO I = 2, N-1-MAX(DERIV,0)
     ! Cycle over each knot accumulating the values for the recurrence.
     DO K = 1, N - I
        ! Check divisors, intervals with 0 width add 0 value to the B-spline.
        LEFT = (KNOTS(K+I-1) - KNOTS(K))
        RIGHT = (KNOTS(K+I) - KNOTS(K+1))
        ! Compute the B-spline recurrence relation (cases based on divisor).
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = &
                   ((XY(:) - KNOTS(K)) / LEFT) * VALUES(:,K) + &
                   ((KNOTS(K+I) - XY(:)) / RIGHT) * VALUES(:,K+1)
           ELSE
              VALUES(:,K) = ((XY(:) - KNOTS(K)) / LEFT) * VALUES(:,K)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = ((KNOTS(K+I) - XY(:)) / RIGHT) * VALUES(:,K+1)
           END IF
        END IF
     END DO
  END DO compute_spline

  ! -------------------- Performing integration ----------------------
  int_or_diff : IF (DERIV .LT. 0) THEN
  ! Integrals will be nonzero on [LAST_KNOT, \infty).
  WHERE (LAST_KNOT .LE. XY(:))
     VALUES(:,N) = 1.0_R8
  END WHERE
  ! Loop through starting at the back, raising the order of all
  ! constituents to match the order of the first.
  raise_order : DO I = 1, ORDER-1
     DO K = N-I, ORDER
        LEFT = (LAST_KNOT - KNOTS(K))
        RIGHT = (LAST_KNOT - KNOTS(K+1))
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = &
                   ((XY(:) - KNOTS(K)) / LEFT) * VALUES(:,K) + &
                   ((LAST_KNOT - XY(:)) / RIGHT) * VALUES(:,K+1)
           ELSE
              VALUES(:,K) = ((XY(:) - KNOTS(K)) / LEFT) * VALUES(:,K)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = ((LAST_KNOT - XY(:)) / RIGHT) * VALUES(:,K+1)
           END IF
        END IF
     END DO
  END DO raise_order
  
  ! Compute the integral(s) of the B-spline.
  compute_integral : DO I = 1, -DERIV
     ! Do a forward evaluation of all constituents.
     DO K = 1, ORDER
        LEFT = (LAST_KNOT - KNOTS(K))
        RIGHT = (LAST_KNOT - KNOTS(K+1))
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = &
                   ((XY(:) - KNOTS(K)) / LEFT) * VALUES(:,K) + &
                   ((LAST_KNOT - XY(:)) / RIGHT) * VALUES(:,K+1)
           ELSE
              VALUES(:,K) = ((XY(:) - KNOTS(K)) / LEFT) * VALUES(:,K)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = ((LAST_KNOT - XY(:)) / RIGHT) * VALUES(:,K+1)
           END IF
        END IF
     END DO
     ! Sum the constituent functions at each knot (from the back).
     DO K = ORDER, 1, -1
        VALUES(:,K) = (VALUES(:,K) + VALUES(:,K+1))
     END DO
     ! Divide by the degree plus the integration coefficient.
     VALUES(:,1) = VALUES(:,1) / (ORDER-1+I)
     ! Rescale then integral by its width.
     VALUES(:,1) = VALUES(:,1) * (LAST_KNOT - KNOTS(1))
     ! Extend the previous two computations if more integrals need
     ! to be computed after this one.
     IF (I+DERIV .LT. 0) THEN
        VALUES(:,2:N) = VALUES(:,2:N) / (ORDER-1+I)
        DO K = 2, N
           VALUES(:,K) = VALUES(:,K) * (LAST_KNOT - KNOTS(K))
        END DO
     END IF
  END DO compute_integral

  ! ------------------ Performing differentiation --------------------
  ELSE IF (DERIV .GT. 0) THEN
  ! Compute the derivative of the B-spline (if D > 0).
  compute_derivative : DO I = N-DERIV, ORDER
     ! Cycle over each knot, following the same structure with the
     ! derivative computing relation instead of the B-spline one.
     DO K = 1, N-I
        ! Assure that the divisor will not cause invalid computations.
        LEFT = (KNOTS(K+I-1) - KNOTS(K))
        RIGHT = (KNOTS(K+I) - KNOTS(K+1))
        ! Compute the derivative recurrence relation.
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) =  (I-1)*(VALUES(:,K)/LEFT - VALUES(:,K+1)/RIGHT)
           ELSE
              VALUES(:,K) = (I-1)*(VALUES(:,K)/LEFT)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              VALUES(:,K) = (I-1)*(-VALUES(:,K+1)/RIGHT)
           END IF
        END IF
     END DO
  END DO compute_derivative
  END IF int_or_diff

  ! Assign the values into the "Y" output.
  XY(:) = VALUES(:,1)
END SUBROUTINE EVAL_BSPLINE

END MODULE SPLINES