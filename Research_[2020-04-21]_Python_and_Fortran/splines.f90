
MODULE REAL_PRECISION  ! module for 64-bit real arithmetic
  INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

MODULE SPLINES
  USE REAL_PRECISION
  IMPLICIT NONE

  PRIVATE :: MAKE_MONOTONE
  PUBLIC

CONTAINS

FUNCTION L2(KNOTS, COEFF, D)
  REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: KNOTS, COEFF
  INTEGER, INTENT(IN), OPTIONAL :: D
  REAL(KIND=R8) :: L2
  ! Local variables.
  REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: Y, KS, CS, MIDS, AB
  INTEGER :: DERIV, K, I, NB, NK, STATUS, NCC, NSPL
  ! Set derivative.
  set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
  ELSE ; DERIV = 0
  END IF set_derivative
  ! Compute the order of the spline pieces for L2.
  K = 1 + 2*(SIZE(KNOTS) - SIZE(COEFF) - DERIV - 1)
  ! Determine the number of unique breakpoints.
  NB = 1
  DO I = 2, SIZE(KNOTS)
     IF (KNOTS(I) .NE. KNOTS(I-1)) NB = NB + 1
  END DO
  ! Allocate the knots and coefficients of the intermediate spline.
  !   NK = NB + (K-2) * (NB-1)
  !      = NB + K*NB - K - 2*NB + 2
  !      = NB*(1 - 2 + K) - K + 2
  !      = NB*(K - 1) - K + 2
  !      = NB*K - NB - K + 2
  !      = K*(NB - 1) - NB + 2
  NK = NB + (K-2)*(NB-1)
  ALLOCATE(KS(1:NK))
  ALLOCATE(CS(1:NK))
  ALLOCATE(Y(1:NK))
  ALLOCATE(MIDS(1:K-2))
  ! Compute the spacing of the internal points on each interval.
  DO I = 1, K-2
     MIDS(I) = REAL(I,KIND=R8) / REAL(K-1, KIND=R8)
  END DO
  ! Compute all the new knots.
  DO I = 1, NB-1
     KS(I*())
  END DO
  ! Compute the values of the spline.
  Y(:) = KS(:)
  CALL EVAL_SPLINE(KNOTS, COEFF, Y, STATUS, D=DERIV)
  Y(:) = Y(:)*Y(:)
  ! Fit the coefficients for the spline that achieve these valeus.
  

  
END FUNCTION L2


SUBROUTINE PMQSI(X, Y, KNOTS, COEFF, STATUS)
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
  !   STATUS -- Integer representing the subroutine execution status.
  !   
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X, Y
  ! 
  ! REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: KNOTS, COEFF
  REAL(KIND=R8), INTENT(OUT), DIMENSION(3*SIZE(X)+6) :: KNOTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(3*SIZE(X)) :: COEFF
  ! ^^^^ THESE ABOVE LINES ARE TEMPORARY, SHOULD NOT BE IN FINAL CODE.
  !      THEY ONLY EXIST TO MAKE AUTOMATIC WRAPPING EASIER.
  INTEGER, INTENT(OUT) :: STATUS
  INTEGER, PARAMETER :: MAX_STEPS=1000
  ! 
  ! Local variables.
  REAL(KIND=R8), DIMENSION(SIZE(X),3) :: VALUES ! Spline derivative values.
  REAL(KIND=R8), DIMENSION(3) :: SLOPES, RESIDUALS ! Linear fit information.
  REAL(KIND=R8), DIMENSION(2) :: V1, V2 ! Violation amounts.
  REAL(KIND=R8), DIMENSION(SIZE(X)) :: MIN_STEP
  INTEGER :: I, Z, STEP
  LOGICAL :: PREV_LINEAR, THIS_LINEAR
  LOGICAL, DIMENSION(SIZE(X)) :: DD_CHANGED
  LOGICAL, DIMENSION(SIZE(X)-1) :: TO_CHECK

  Z = SIZE(X)
  ! Check the shape of incoming arrays.
  IF      (Z .LT. 1)                 THEN ; STATUS = 1 ; RETURN
  ELSE IF (SIZE(Y) .NE. Z)           THEN ; STATUS = 2 ; RETURN
  ELSE IF (SIZE(KNOTS) .LT. 3*Z + 6) THEN ; STATUS = 4 ; RETURN
  ELSE IF (SIZE(COEFF) .LT. 3*Z)     THEN ; STATUS = 5 ; RETURN
  END IF

  ! Verify that X are increasing.
  DO I = 1, Z-1
     IF (X(I) .GE. X(I+1)) THEN ; STATUS = 6 ; RETURN ; END IF
  END DO

  ! Copy the "Y" values into the "VALUES" array, initialize others to zero.
  VALUES(:,1) = Y(:)
  VALUES(:,2:3) = 0_R8
  ! Initialize slopes, residuals, and "is linear" booleans.
  SLOPES(:) = 0.0_R8
  RESIDUALS(:) = HUGE(1.0_R8)
  PREV_LINEAR = .FALSE.
  THIS_LINEAR = .FALSE.

  ! Construct local derivative and second derivative fits for each point.
  DO I = 2, Z-1
     ! If this point has slope 0 between it and either neighbor..
     IF ((Y(I-1) .EQ. Y(I)) .OR. (Y(I) .EQ. Y(I+1))) THEN ; CYCLE
     ! Else if this is an extreme point (with different slopes on either side)..
     ELSE IF ((Y(I)-Y(I-1)) * (Y(I+1)-Y(I)) .LT. 0.0_R8) THEN
        ! Construct a local quadratic to determine the second
        ! derivative value.
        VALUES(I,3) = 0.0_R8
     ! Else if this point is inside a monotone segment..
     ELSE
        ! Get local linear fit slope and residual information.
        CALL LOCAL_LINEAR(I, SLOPES(3), RESIDUALS(3))
        THIS_LINEAR = .TRUE.
     END IF
     ! Use the slope of the best-fit linear model of the three for
     ! the previous interval, if applicable.
     IF (PREV_LINEAR) THEN
        VALUES(I-1, 2) = SLOPES(MINLOC(RESIDUALS(:), 1))
     END IF
     PREV_LINEAR = THIS_LINEAR
     THIS_LINEAR = .FALSE.
     ! Rotate existing slopes and residuals to free a space for the
     ! slope of the next interval (provides a sliding window).
     SLOPES(2) = SLOPES(3)
     RESIDUALS(1) = RESIDUALS(2)
     RESIDUALS(2) = RESIDUALS(3)
     SLOPES(3) = 0.0_R8
     RESIDUALS(3) = HUGE(1.0_R8)
  END DO
  ! Use the slope of the best-fit linear model of the three for
  ! the previous interval, if applicable.
  IF (PREV_LINEAR) THEN
     VALUES(I-1, 2) = SLOPES(MINLOC(RESIDUALS(:), 1))
  END IF
  ! Estimate the first and second derivative at edges using a
  ! quadratic interpolant over the neighboring slope and edge value.
  IF (Z .GT. 2) THEN
     VALUES(1,2) = (Y(2) - Y(1)) / (X(2) - X(1))
     VALUES(Z,2) = (Y(Z) - Y(Z-1)) / (X(Z) - X(Z-1))
  END IF

  ! Cycle through and enforce monotonicity on the derivative and
  ! second derivative values.
  TO_CHECK(:) = .TRUE.
  DD_CHANGED(:) = .FALSE.
  FORALL (I=1:Z) MIN_STEP(I) = ABS(VALUES(I,2)) / REAL(MAX_STEPS,KIND=R8)
  STEP = 1
  DO WHILE (ANY(TO_CHECK(:)))
     STEP = STEP + 1
     IF (STEP .GT. 100) EXIT
     ! Check all intervals.
     interval_loop : DO I = 1, Z-1
        monotone_check : IF (TO_CHECK(I)) THEN
           TO_CHECK(I) = .FALSE.
           CALL MAKE_MONOTONE(X(I), X(I+1), VALUES(I,:), VALUES(I+1,:), V1(:), V2(:))
           ! Convert the "violation" amounts to strictly positive numbers.
           V1(:) = ABS(V1(:))
           V2(:) = ABS(V2(:))
           print *, ''
           print *, "checking interval ", I
           print *, "   ", VALUES(I,2:)
           print *, "   ", V1(:)
           PRINT *, ''
           print *, "   ", VALUES(I+1,2:)
           print *, "   ", V2(:)
           ! Check the left side.
           check_left : IF (MAXVAL(V1(:)) .GT. 0.0_R8) THEN
              ! Mark the previous and next interval as "need to be checked".
              IF (I .GT. 1) TO_CHECK(I-1) = .TRUE.
              ! If this value has not been changed before, then mark
              ! it and move on.
              IF (.NOT. DD_CHANGED(I)) THEN ; DD_CHANGED(I) = .TRUE.
              ! The second derivative value was changed and this
              ! breakpoint is being revisited, shrink first derivative
              ! if it is greater than zero.
              ELSE IF (ABS(VALUES(I,2)) .GT. 0.0_R8) THEN
                 ! The first derivative was modified or will be,
                 ! reset "DD changed" variable.
                 DD_CHANGED(I) = .FALSE.
                 ! If the first derivative was not changed, it needs to be.
                 IF (V1(1) .EQ. 0.0_R8) THEN
                    ! Compute the "multiplier" that determines the shrinkage.
                    V1(1) = (1.0_R8 - 1.0_R8 / &
                         (1.0_R8 + (V1(2) + 1.0_R8/ABS(VALUES(I,2)))))
                    ! Make sure the shrinkage is at least MIN_STEP(I).
                    V1(1) = MIN(V1(1), MIN_STEP(I) / ABS(VALUES(I,2)))
                    ! Shrink the derivative value.
                    VALUES(I,2) = V1(1) * VALUES(I,2)
                 END IF
              ! The first derivative is zero and the second derivative
              ! has been changed twice, so there must be a first
              ! derivative conflict with the neighboring interval.
              ELSE IF (I .GT. 1) THEN
                 DD_CHANGED(I-1) = .TRUE.
              END IF
           END IF check_left

           ! Check the right side.
           check_right : IF (MAXVAL(V2(:)) .GT. 0.0_R8) THEN
              ! Mark the next interval as "need to be checked".
              IF (I .LT. Z-1) TO_CHECK(I+1) = .TRUE.
              ! If this value has not been changed before, then mark
              ! it and move on.
              IF (.NOT. DD_CHANGED(I+1)) THEN ; DD_CHANGED(I+1) = .TRUE.
              ! The second derivative value was changed and this
              ! breakpoint is being revisited, shrink first derivative
              ! if it is greater than zero.
              ELSE IF (ABS(VALUES(I+1,2)) .GT. 0.0_R8) THEN
                 ! The first derivative was modified or will be,
                 ! reset "DD changed" variable.
                 DD_CHANGED(I+1) = .FALSE.
                 ! If the first derivative was not changed, it needs to be.
                 IF (V2(1) .EQ. 0.0_R8) THEN
                    ! Compute the "multiplier" that determines the shrinkage.
                    V2(1) = (1.0_R8 - 1.0_R8 / &
                         (1.0_R8 + (V2(2) + 1.0_R8/ABS(VALUES(I+1,2)))))
                    ! Make sure the shrinkage is at least MIN_STEP(I+1).
                    V2(1) = MIN(V2(1), MIN_STEP(I) / ABS(VALUES(I+1,2)))
                    ! Shrink the derivative value.
                    VALUES(I+1,2) = V2(1) * VALUES(I+1,2)
                 END IF
              ! The first derivative is zero and the second derivative
              ! has been changed twice, so there must be a first
              ! derivative conflict with the neighboring interval.
              ELSE IF (I .LT. Z-1) THEN
                 DD_CHANGED(I+1) = .TRUE.
              END IF
           END IF check_right
        END IF monotone_check
     END DO interval_loop
  END DO

  ! Use "FIT_SPLINE" to fit the final interpolating spline.
  CALL FIT_SPLINE(X, VALUES, KNOTS, COEFF, STATUS)

CONTAINS
  SUBROUTINE LOCAL_LINEAR(I, A, R)
    ! Given an index "I", compute the linear regression about the
    ! points X(I-1:I+1), Y(I-1:I+1). Return the slope and the sum of
    ! squared errors of the regression line.
    INTEGER, INTENT(IN) :: I
    REAL(KIND=R8), INTENT(OUT) :: A, & ! Slope of the linear regression.
         R ! Sum of squared errors of the linear regression.
    ! Local variables.
    REAL(KIND=R8) :: B, & ! Y-intercept of the linear regression.
         AN, & ! Numerator for computing "A".
         BN, & ! Numerator for computing "B".
         D ! Denominator for computing "A" and "B".
    
    AN = 2*(X(I-1)*Y(I-1) + X(I)*Y(I) + X(I+1)*Y(I+1)) &
         - X(I-1)*(Y(I)+Y(I+1)) &
         - X(I)*(Y(I-1)+Y(I+1)) &
         - X(I+1)*(Y(I-1)+Y(I))
    BN = Y(I-1)*(X(I)**2 + X(I+1)**2 - X(I-1)*(X(I) + X(I+1))) &
         + Y(I)*(X(I-1)**2 + X(I+1)**2 - X(I)*(X(I-1) + X(I+1))) &
         + Y(I+1)*(X(I-1)**2 + X(I)**2 - X(I+1)*(X(I-1) + X(I)))
    D = 2*(X(I-1)**2 + X(I)**2 + X(I+1)**2) &
         - 2*(X(I)*X(I+1) + X(I-1)*(X(I)+X(I+1)))
    ! Compute the regression coefficient and intercept.
    A = AN / D
    B = BN / D
    ! Compute the sum of squared errors.
    R = (A*X(I-1) + B - Y(I-1))**2 &
         + (A*X(I) + B - Y(I))**2 &
         + (A*X(I+1) + B - Y(I+1))**2
    ! Compute the local slope, intercept
  END SUBROUTINE LOCAL_LINEAR
END SUBROUTINE PMQSI


SUBROUTINE MAKE_MONOTONE(X1, X2, Y1, Y2, D1, D2)
  ! Given an interval and function values (value, first, and seocond
  ! derivative), make the quintic over this interval monotone by
  ! modifying the first and second derivatives according to quintic
  ! monotonicity theory. Return the change applied to first and second
  ! derivative values.
  REAL(KIND=R8), INTENT(IN) :: X1, X2
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(3) :: Y1, Y2
  REAL(KIND=R8), INTENT(OUT), DIMENSION(2) :: D1, D2
  ! Local variables.
  REAL(KIND=R8) :: A, A0, A1, B, B0, B1, ETA1, ETA2, G0, G1, &
       SIGN, SLOPE, TEMP, V, W

  ! Compute useful variables, store original derivative values.
  V = Y2(1) - Y1(1)
  W = X2 - X1
  SLOPE = V / W
  D1(1:2) = Y1(2:3)
  D2(1:2) = Y2(2:3)
  ! Handle unchanging interval.
  IF (SLOPE .EQ. 0.0_R8) THEN
     ! Zero out the first and second derivative for flat intervals and return.
     Y1(2:3) = 0.0_R8
     Y2(2:3) = 0.0_R8
     D1(1:2) = Y1(2:3) - D1(1:2)
     D2(1:2) = Y2(2:3) - D2(1:2)
     RETURN
  ! This is not a flat interval, set the "SIGN" variable.
  ELSE IF (SLOPE .GT. 0.0_R8) THEN ; SIGN =  1.0_R8
  ELSE                             ; SIGN = -1.0_R8
  END IF

  ! Consider this interval to be monotone increasing.
  V = V * SIGN
  Y1(:) = Y1(:) * SIGN
  Y2(:) = Y2(:) * SIGN
  SLOPE = SLOPE * SIGN
  ! Set derivative to be the median of {0, Y(2), 14*SLOPE}.
  Y1(2) = MIN(14.0_R8*SLOPE, MAX(0.0_R8, Y1(2)))
  Y2(2) = MIN(14.0_R8*SLOPE, MAX(0.0_R8, Y2(2)))
  ! Compute "A" (left ratio) and "B" (right ratio).
  A = Y1(2) / SLOPE
  B = Y2(2) / SLOPE

  ! Use a (simplified) monotone cubic over this region if AB = 0.
  ! Only consider the coefficients less than the x^4, because that
  ! term is strictly positive.
  simplified_monotone : IF (A*B .LT. SQRT(EPSILON(0.0_R8))) THEN
     ! Ensure that DDf(X1) has a nonempty feasible region (given Df).
     TEMP = MAX(0.0_R8, (20.0_R8*V / W) / (5.0_R8*Y1(2) + 4.0_R8*Y2(2)))
     IF (TEMP .LT. 1) THEN
        Y1(2) = Y1(2) * TEMP
        Y2(2) = Y2(2) * TEMP
     END IF
     ! Cap DDf(X2) so that DDf(X1) feasible region is nonempty.
     Y1(3) = MIN(Y1(3), (4.0_R8*(2.0_R8*Y1(2) + Y2(2)) + 20.0_R8*V/W) / W)
     ! Enforce gamma >= delta.
     Y1(3) = MAX(Y1(3), 3.0_R8*Y1(2) / W)
     ! Enforce \alpha >= 0.
     Y2(3) = MIN(Y2(3), -4.0_R8*Y2(2) / W)
     ! Enforce \beta >= \alpha.
     Y2(3) = MAX(Y2(3), (3.0_R8*Y1(3)*W - &
          (24.0_R8*Y1(2) + 32.0_R8*Y2(2)) - 60.0_R8*V / W) / (5.0_R8*W))
     ! Undo the sign change, compute the derivative changes, and return.
     Y1(:) = Y1(:) * SIGN
     Y2(:) = Y2(:) * SIGN
     D1(1:2) = Y1(2:3) - D1(1:2)
     D2(1:2) = Y2(2:3) - D2(1:2)
     RETURN
  END IF simplified_monotone

  ! Clip derivative values that are too large (to ensure that shrinking
  ! the derivative vectors on either end will not break monotonicity).
  !   (clipping at 6 box is enough, using 8 box requires more steps)
  TEMP = 6.0_R8 / MAX(A,B)
  IF (TEMP .LT. 1.0_R8) THEN
     Y1(2) = Y1(2) * TEMP
     Y2(2) = Y2(2) * TEMP     
     A = A * TEMP
     B = B * TEMP
  END IF

  IF (24.0_R8 + 2.0_R8*SQRT(A*B) - 3.0_R8*(A+B) .LT. 0.0_R8) THEN
     PRINT *, "ERROR: Something went wrong, Tau_1 is negative."
  END IF

  ! Compute the terms needed to simplify the monotonicity check in
  ! terms of the second derivative values.
  TEMP = (B / A)**(0.75_R8) / Y2(2)
  G0 = 4.0_R8 * Y1(2) * TEMP
  G1 = W * TEMP
  TEMP = (B / A) ** (0.25)
  A0 = 4.0_R8 * TEMP
  A1 = - W / Y2(2) * TEMP
  TEMP = W / (V * SQRT(A) * SQRT(B))
  B0 = 30.0_R8 - 12.0_R8 * (Y1(2) + Y2(2)) * TEMP
  B1 = (-3.0_R8 * W / 2.0_R8) * TEMP

  ! Perform a binary search for allowable DDf values.
  A = SQRT(A) ; B = SQRT(B)
  ETA1 = - A * (7.0_R8*A + 3.0_R8*B) * SLOPE / W
  ETA2 =   B * (3.0_R8*A + 7.0_R8*B) * SLOPE / W
  nonmonotone_DDf : IF (.NOT. IS_MONOTONE(Y1(3), Y2(3))) THEN
     ! Reuse variables "A" and "B" as lower and upper bounds for
     ! binary search.
     A = 0.0_R8 ; B = 1.0_R8
     binary_search : DO WHILE ((B - A) .GT. SQRT(EPSILON(0.0_R8)))
        TEMP = (A + B) / 2.0_R8
        ! Compute current DDf estimates.
        Y1(3) = (1.0_R8 - TEMP) * D1(2)  +  TEMP * ETA1
        Y2(3) = (1.0_R8 - TEMP) * D2(2)  +  TEMP * ETA2
        ! Depending on the monotonicity, shrink search region.
        IF (IS_MONOTONE(Y1(3), Y2(3))) THEN ; B = TEMP
        ELSE                                ; A = TEMP
        END IF
     END DO binary_search
     ! Compute the final (monotone) estimate.
     Y1(3) = (1.0_R8 - B) * D1(2)  +  B * ETA1
     Y2(3) = (1.0_R8 - B) * D2(2)  +  B * ETA2
  END IF nonmonotone_DDf

  ! Undo the sign change, compute the derivative changes, and return.
  Y1(:) = Y1(:) * SIGN
  Y2(:) = Y2(:) * SIGN
  D1(1:2) = Y1(2:3) - D1(1:2)
  D2(1:2) = Y2(2:3) - D2(1:2)
  RETURN

CONTAINS
  FUNCTION IS_MONOTONE(DDY1, DDY2)
    ! Compute the condition that determines if monotonicity is
    ! achieved by current second derivative values.
    REAL(KIND=R8), INTENT(IN) :: DDY1, DDY2
    REAL(KIND=R8) :: G, A, B
    LOGICAL :: IS_MONOTONE
    G = G0 + G1*DDY1
    A = A0 + A1*DDY2
    B = B0 + B1*(DDY1 - DDY2)
    IF (B .LT. 6.0_R8) THEN ; IS_MONOTONE = (A .GT. -(B+2)/2.0_R8)
    ELSE                    ; IS_MONOTONE = (G .GT. -2.0_R8*SQRT(B-2.0_R8))
    END IF
  END FUNCTION IS_MONOTONE
END SUBROUTINE MAKE_MONOTONE


SUBROUTINE FIT_SPLINE(BREAKPOINTS, VALUES, KNOTS, COEFF, STATUS)
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
  !   STATUS -- Integer representing the subroutine execution status:
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
  !   produce a correct set of coefficients and return STATUS code 7.
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
  INTEGER, INTENT(OUT) :: STATUS
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
  STATUS = 0

  ! Check the shape of incoming arrays.
  IF      (NB .LT. 1)              THEN ; STATUS = 1 ; RETURN
  ELSE IF (NSPL .LT. 1)            THEN ; STATUS = 2 ; RETURN
  ELSE IF (SIZE(VALUES,1) .NE. NB) THEN ; STATUS = 3 ; RETURN
  ELSE IF (SIZE(KNOTS) .LT. NK)    THEN ; STATUS = 4 ; RETURN
  ELSE IF (SIZE(COEFF) .LT. NSPL)  THEN ; STATUS = 5 ; RETURN
  END IF
  ! Verify that BREAKPOINTS are increasing.
  DO I = 1, NB - 1
     IF (BREAKPOINTS(I) .GE. BREAKPOINTS(I+1)) THEN
        STATUS = 6 ; RETURN
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
        CALL EVAL_BSPLINE(KNOTS(I:J), AB(I1+DERIV:I2:NCC,I), STATUS, D=DERIV)
        ! ^ Correct usage is inherently enforced, only extrapolation
        !   warnings will be produced by this call (STATUS=1). These
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
       COEFF, NSPL, STATUS)
  ! Check for errors in the execution of DGBSV, (this should not happen).
  IF (STATUS .NE. 0) THEN ; STATUS = STATUS + 10 ; RETURN ; END IF
  ! Check to see if the linear system was correctly solved by looking at
  ! the difference between prouduced B-spline values and provided values.
  MAX_ERROR = SQRT(SQRT(EPSILON(1.0_R8)))
  DO DERIV = 0, NCC-1
     ! Reuse the first row of AB as scratch space (the first column
     ! might not be large enough, but the first row certainly is).
     AB(1,1:NB) = BREAKPOINTS(:)
     ! Evaluate this spline at all breakpoints. Correct usage is
     ! enforced here, so it is expected that STATUS=0 always.
     CALL EVAL_SPLINE(KNOTS, COEFF, AB(1,1:NB), STATUS, D=DERIV)
     ! Check the precision of the reproduced values.
     ! Return an error if the precision is too low.
     IF (MAXVAL(ABS((AB(1,1:NB) - VALUES(:,DERIV+1)) &
          / (1.0_R8 + ABS(VALUES:,DERIV+1)))) .GT. MAX_ERROR) THEN
        STATUS = 7
        RETURN
     END IF
  END DO
END SUBROUTINE FIT_SPLINE


SUBROUTINE EVAL_SPLINE(KNOTS, COEFF, XY, STATUS, D)
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
  !   STATUS -- Integer representing subroutine exeuction status.
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
  INTEGER, INTENT(OUT) :: STATUS
  ! Local variables.
  INTEGER :: DERIV, I, ORDER
  REAL(KIND=R8), DIMENSION(SIZE(XY), SIZE(COEFF)) :: VALUES

  ! Check for size-related errors and for an empty knot interior.
  IF (SIZE(KNOTS) .LE. 1)                    THEN ; STATUS = 3 ; RETURN
  ELSE IF (KNOTS(1) .EQ. KNOTS(SIZE(KNOTS))) THEN ; STATUS = 4 ; RETURN
  ELSE IF (SIZE(KNOTS) .LE. SIZE(COEFF))     THEN ; STATUS = 5 ; RETURN
  END IF
  ! Check for valid (nondecreasing) knot sequence.
  DO I = 1, SIZE(KNOTS)-1
     IF (KNOTS(I) .GT. KNOTS(I+1)) THEN ; STATUS = 2 ; RETURN ; END IF
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
          STATUS, D=DERIV)
     ! ^ Correct usage is inherently enforced, only extrapolation
     !   warnings will be produced by this call. These
     !   extrapolation warnings are expected because underlying
     !   B-splines may not support the full interval.
  END DO
  ! Set the EXTRAPOLATION status flag.
  IF ((MINVAL(XY(:)) .LT. KNOTS(1)) .OR. &
       (MAXVAL(XY(:)) .GE. KNOTS(SIZE(KNOTS)))) THEN ; STATUS = 1
  ELSE ; STATUS = 0
  END IF
  ! Store the values into Y as the weighted sums of B-spline evaluations.
  XY(:) = MATMUL(VALUES(:,:), COEFF(:))
END SUBROUTINE EVAL_SPLINE


SUBROUTINE EVAL_BSPLINE(KNOTS, XY, STATUS, D)
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
  !   STATUS -- Execution status of this subroutine on exit.
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
  INTEGER, INTENT(OUT) :: STATUS
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
  STATUS = 0
  ! Check for valid knot sequence.
  IF (N .LE. 1) THEN ; STATUS = 3 ; RETURN ; END IF
  DO K = 1, N-1
     IF (KNOTS(K) .GT. KNOTS(K+1)) THEN ; STATUS = 2 ; RETURN ; END IF
  END DO
  ! Check for extrapolation, set status if it is happening, but continue.
  IF ((MINVAL(XY(:)) .LT. KNOTS(1)) .OR. (MAXVAL(XY(:)) .GE. LAST_KNOT)) &
       STATUS = 1
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
