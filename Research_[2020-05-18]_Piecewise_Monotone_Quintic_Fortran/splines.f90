
MODULE REAL_PRECISION  ! module for 64-bit real arithmetic
  INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

MODULE SPLINES
  USE REAL_PRECISION
  IMPLICIT NONE

  PRIVATE :: EVAL_BSPLINE
  PUBLIC

CONTAINS

SUBROUTINE PMQSI(X, Y, T, BCOEF, INFO)
  ! Compute a piecewise monotone quintic spline interpolant (PMQSI)
  ! for data in terms of a B-spline basis represented with knots "T"
  ! and spline coefficients "BCOEF".
  ! 
  ! INPUT:
  !   X(1:Z) -- An array of increasing real values.
  !   Y(1:Z) -- An array of (real-valued) function values such that 
  !             Y(I) = F(X(I)) for some piecewise monotone function F.
  ! 
  ! OUTPUT:
  !   T(1:3*Z+6) -- The knots for the PMQSI B-spline basis.
  !   BCOEF(1:3*Z) -- The coefficients for the PMQSI B-spline basis.
  !   INFO -- Integer representing the subroutine execution status.
  !   
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X, Y
  ! REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: T, BCOEF
  REAL(KIND=R8), INTENT(OUT), DIMENSION(3*SIZE(X)+6) :: T
  REAL(KIND=R8), INTENT(OUT), DIMENSION(3*SIZE(X)) :: BCOEF
  ! ^^^^ THESE ABOVE LINES ARE TEMPORARY, SHOULD NOT BE IN FINAL CODE.
  !      THEY ONLY EXIST TO MAKE AUTOMATIC PYTHON WRAPPING EASIER.
  INTEGER, INTENT(OUT) :: INFO
  ! 
  ! Local variables.
  REAL(KIND=R8), DIMENSION(SIZE(X),3) :: FX ! Spline function values.
  REAL(KIND=R8), DIMENSION(SIZE(X),2) :: FHATX ! Estimated function values.
  LOGICAL, DIMENSION(SIZE(X)) :: CHECKING, GROWING, SHRINKING ! Execution queues
  INTEGER, DIMENSION(SIZE(X)) :: TO_CHECK, TO_GROW, TO_SHRINK !  ... cont'd
  REAL(KIND=R8) :: A, ACCURACY, B, DX, DIRECTION, STEP_SIZE, MAX_REAL
  INTEGER :: I, J, K, NC, NG, NS, Z
  LOGICAL :: SEARCHING
  ! Declare constants for computation.
  Z = SIZE(X)
  MAX_REAL = HUGE(1.0_R8)
  ACCURACY = SQRT(EPSILON(1.0_R8))
  ! Check the shape of incoming arrays.
  IF      (Z .LT. 3)             THEN ; INFO = 1 ; RETURN
  ELSE IF (SIZE(Y) .NE. Z)       THEN ; INFO = 2 ; RETURN
  ELSE IF (SIZE(T) .LT. 3*Z + 6) THEN ; INFO = 4 ; RETURN
  ELSE IF (SIZE(BCOEF) .LT. 3*Z) THEN ; INFO = 5 ; RETURN
  END IF
  ! Verify that X are increasing and appropriately spaced.
  DO I = 1, Z-1
     IF (X(I+1) - X(I) .LE. ACCURACY) THEN ; INFO = 6 ; RETURN ; END IF
  END DO
  ! ==================================================================
  !            Estimate initial derivatives by using a
  !            minimum curvature quadratic facet model.
  ! 
  ! Copy the "Y" values into the "FX" array.
  FX(:,1) = Y(:)
  ! Identify local extreme points and flat points. Denote location
  ! of flats in "GROWING" and extrema in "SHRINKING".
  GROWING(1) = (Z .GT. 1) .AND. (Y(1) .EQ. Y(2))
  GROWING(Z) = (Z .GT. 1) .AND. (Y(Z-1) .EQ. Y(Z))
  SHRINKING(1) = .FALSE.
  SHRINKING(Z) = .FALSE.
  DO I = 2, Z-1
     DIRECTION = (Y(I) - Y(I-1)) * (Y(I+1) - Y(I))
     GROWING(I) = DIRECTION .EQ. 0.0_R8
     SHRINKING(I) = DIRECTION .LT. 0.0_R8
  END DO
  ! Use local quadratic interpolants to estimate slopes and second
  ! derivatives at all points. Use zero-sloped quadratic interpolant
  ! of left/right neighbors at extrema to estimate curvature.
  estimate_derivatives : DO I = 1, Z
     ! Precompute these integers, they are used repeatedly.
     J = I+1
     K = I-1
     ! Determine the direction of change at the point "I".
     IF (I .EQ. 1) THEN
        IF (I .EQ. Z) THEN ;            DIRECTION =  0.0_R8
        ELSE IF (Y(I) .LT. Y(J)) THEN ; DIRECTION =  1.0_R8
        ELSE IF (Y(I) .GT. Y(J)) THEN ; DIRECTION = -1.0_R8
        ELSE ;                          DIRECTION =  0.0_R8
        END IF
     ELSE
        IF      (Y(K) .LT. Y(I)) THEN ; DIRECTION =  1.0_R8
        ELSE IF (Y(K) .GT. Y(I)) THEN ; DIRECTION = -1.0_R8
        ELSE ;                          DIRECTION =  0.0_R8
        END IF
     END IF
     ! Initialize the curvature to be maximally large.
     FX(I,3) = MAX_REAL
     ! If this is a local flat, first and second derivatives are zero valued.
     pick_quadratic : IF (GROWING(I)) THEN ; FX(I,2:3) = 0.0_R8
     ! If this is an extreme point, construct quadratic interpolants
     ! that have zero slope here and hit left/right neighbors.
     ELSE IF (SHRINKING(I)) THEN
        ! Set the first derivative to zero.
        FX(I,2) = 0.0_R8
        ! Compute the coefficient A in "Ax^2 + Bx + C" that interpolates X(I-1).
        FX(I,3) = (Y(K) - Y(I)) / (X(K) - X(I))**2
        ! Compute the coefficient A in "Ax^2 + Bx + C" that interpolates X(I+1).
        A = (Y(J) - Y(I)) / (X(J) - X(I))**2
        IF (ABS(A) .LT. ABS(FX(I,3))) THEN ; FX(I,3) = A ; END IF
        ! Compute the actual second derivative (instead of coefficient A).
        FX(I,3) = MAX(MIN(2.0_R8 * FX(I,3), MAX_REAL), -MAX_REAL)
     ELSE
        ! --------------------
        ! Quadratic left of I.
        IF (K .GT. 1) THEN
           ! If there's a zero-derivative to the left, use it's right interpolant.
           IF (SHRINKING(K) .OR. GROWING(K)) THEN
              A = (Y(I) - Y(K)) / (X(I) - X(K))**2
              B = -2.0_R8 * X(K) * A
           ! Otherwise use the standard quadratic on the left.
           ELSE ; CALL QUADRATIC(K)
           END IF
           DX = MAX(MIN(2.0_R8*A*X(I) + B, MAX_REAL), -MAX_REAL)
           IF (DX*DIRECTION .GE. 0.0_R8) THEN
              FX(I,2) = DX
              FX(I,3) = A
           END IF
        END IF
        ! ------------------------
        ! Quadratic centered at I (require that it has at least one
        ! neighbor that is not forced to zero slope).
        IF ((I .GT. 1) .AND. (I .LT. Z) .AND. .NOT. &
             ((SHRINKING(K) .OR. GROWING(K)) .AND. &
             (SHRINKING(J) .OR. GROWING(J)))) THEN
           ! Construct the quadratic interpolant through this point and neighbors.
           CALL QUADRATIC(I)
           DX = MAX(MIN(2.0_R8*A*X(I) + B, MAX_REAL), -MAX_REAL)
           ! Keep this new quadratic if it has less curvature.
           IF ((DX*DIRECTION .GE. 0.0_R8) .AND. (ABS(A) .LT. ABS(FX(I,3)))) THEN
              FX(I,2) = DX
              FX(I,3) = A
           END IF
        END IF
        ! ---------------------
        ! Quadratic right of I.
        IF (J .LT. Z) THEN
           ! If there's an extreme point to the right, use it's left interpolant.
           IF (SHRINKING(J) .OR. GROWING(J)) THEN
              A = (Y(I) - Y(J)) / (X(I) - X(J))**2
              B = -2.0_R8 * X(J) * A
           ! Otherwise use the standard quadratic on the right.
           ELSE ; CALL QUADRATIC(J)
           END IF
           DX = MAX(MIN(2.0_R8*A*X(I) + B, MAX_REAL), -MAX_REAL)
           ! Keep this new quadratic if it has less curvature.
           IF ((DX*DIRECTION .GE. 0.0_R8) .AND. (ABS(A) .LT. ABS(FX(I,3)))) THEN
              FX(I,2) = DX
              FX(I,3) = A
           END IF
        END IF
        ! Set the final quadratic.
        IF (FX(I,3) .EQ. MAX_REAL) THEN ; FX(I,2:3) = 0.0_R8
        ! Compute the actual curvature of the quadratic, instead of coefficient on x^2.        
        ELSE ; FX(I,3) = MAX(MIN(2.0_R8 * FX(I,3), MAX_REAL), -MAX_REAL)
        END IF
     END IF pick_quadratic
  END DO estimate_derivatives
  ! ------------------------------------------------------------------

  ! ==================================================================
  !           Identify viable piecewise monotone derivative
  !             values by doing a quasi-bisection search.
  ! 
  ! Store the initially estimated values.
  FHATX(:,1:2) = FX(:,2:3)
  ! Identify which intervals are not monotone and need to be fixed.
  CHECKING(:) = .FALSE.
  GROWING(:) = .FALSE.
  SHRINKING(:) = .FALSE.
  NC = 0 ; NG = 0 ; NS = 0
  DO I = 1, Z-1
     J = I+1
     IF (.NOT. IS_MONOTONE(X(I), X(J), FX(I,1), FX(J,1), &
          FX(I,2), FX(J,2), FX(I,3), FX(J,3))) THEN
        IF (.NOT. SHRINKING(I)) THEN
           SHRINKING(I) = .TRUE. ; NS = NS+1 ; TO_SHRINK(NS) = I
        END IF
        IF (.NOT. SHRINKING(J)) THEN
           SHRINKING(J) = .TRUE. ; NS = NS+1 ; TO_SHRINK(NS) = J
        END IF
     END IF
  END DO
  ! Initialize step size to 1.0 (will be halved at beginning of loop).
  STEP_SIZE = 1.0_R8
  SEARCHING = .TRUE.
  ! Loop until the accuracy is achieved and *all* intervals are monotone.
  DO WHILE (SEARCHING .OR. (NS .GT. 0))
     ! Compute the step size for this iteration.
     IF (SEARCHING) THEN
        STEP_SIZE = STEP_SIZE / 2.0_R8
        IF (STEP_SIZE .LT. ACCURACY) THEN
           SEARCHING = .FALSE. ; STEP_SIZE = ACCURACY ; NG = 0
        END IF
     ! Grow the step size (at a slower rate) if there are still intervals to fix.
     ELSE ; STEP_SIZE = STEP_SIZE + STEP_SIZE / 2.0_R8
     END IF
     ! Grow all those first and second derivatives that were previously 
     ! shrunk, but are not currently affecting monotonicity.
     grow_values : DO J = 1, NG
        I = TO_GROW(J)
        ! Do not grow values that are actively related to a nonmonotone segment.
        IF (SHRINKING(I)) CYCLE grow_values
        ! Otherwise, grow those values that have been modified previously.
        FX(I,2) = FX(I,2) + STEP_SIZE * FHATX(I,1)
        FX(I,3) = FX(I,3) + STEP_SIZE * FHATX(I,2)
        ! Make sure the first derivative does not exceed its original value.
        IF ((FHATX(I,1) .LT. 0.0_R8) .AND. (FX(I,2) .LT. FHATX(I,1))) THEN
           FX(I,2) = FHATX(I,1)
        ELSE IF ((FHATX(I,1) .GT. 0.0_R8) .AND. (FX(I,2) .GT. FHATX(I,1))) THEN
           FX(I,2) = FHATX(I,1)
        END IF
        ! Make sure the second derivative does not exceed its original value.
        IF ((FHATX(I,2) .LT. 0.0_R8) .AND. (FX(I,3) .LT. FHATX(I,2))) THEN
           FX(I,3) = FHATX(I,2)
        ELSE IF ((FHATX(I,2) .GT. 0.0_R8) .AND. (FX(I,3) .GT. FHATX(I,2))) THEN
           FX(I,3) = FHATX(I,2)
        END IF
        ! Set this value and its neighbors to be checked for monotonicity.
        IF ((I .GT. 1) .AND. (.NOT. CHECKING(I-1))) THEN
           CHECKING(I-1) = .TRUE. ; NC = NC+1 ; TO_CHECK(NC) = I-1
        END IF
        IF ((I .LT. Z) .AND. (.NOT. CHECKING(I))) THEN
           CHECKING(I) = .TRUE. ; NC = NC+1 ; TO_CHECK(NC) = I
        END IF
     END DO grow_values
     ! Shrink all those first and second derivatives that cause nonmonotonicity.
     shrink_values : DO J = 1, NS
        I = TO_SHRINK(J)
        SHRINKING(I) = .FALSE.
        IF (SEARCHING .AND. (.NOT. GROWING(I))) THEN
           GROWING(I) = .TRUE. ; NG = NG+1 ; TO_GROW(NG) = I
        END IF
        ! Shrink those values that are causing nonmonotonicity.
        FX(I,2) = FX(I,2) - STEP_SIZE * FHATX(I,1)
        FX(I,3) = FX(I,3) - STEP_SIZE * FHATX(I,2)
        ! Make sure the first derivative does not pass zero.
        IF ((FHATX(I,1) .LT. 0.0_R8) .AND. (FX(I,2) .GT. 0.0_R8)) THEN
           FX(I,2) = 0.0_R8
        ELSE IF ((FHATX(I,1) .GT. 0.0_R8) .AND. (FX(I,2) .LT. 0.0_R8)) THEN
           FX(I,2) = 0.0_R8
        END IF
        ! Make sure the second derivative does not pass zero.
        IF ((FHATX(I,2) .LT. 0.0_R8) .AND. (FX(I,3) .GT. 0.0_R8)) THEN
           FX(I,3) = 0.0_R8
        ELSE IF ((FHATX(I,2) .GT. 0.0_R8) .AND. (FX(I,3) .LT. 0.0_R8)) THEN
           FX(I,3) = 0.0_R8
        END IF
        ! Set this point and its neighbors to be checked.
        IF ((I .GT. 1) .AND. (.NOT. CHECKING(I-1))) THEN
           CHECKING(I-1) = .TRUE. ; NC = NC+1 ; TO_CHECK(NC) = I-1
        END IF
        IF ((I .LT. Z) .AND. (.NOT. CHECKING(I))) THEN
           CHECKING(I) = .TRUE. ; NC = NC+1 ; TO_CHECK(NC) = I
        END IF
     END DO shrink_values
     ! Go through and identify which values are associated with
     ! nonmonotone intervals after the updates.
     NS = 0
     check_monotonicity : DO K = 1, NC
        I = TO_CHECK(K)
        J = I+1
        CHECKING(I) = .FALSE.
        IF (.NOT. IS_MONOTONE(X(I), X(J), FX(I,1), FX(J,1), &
             FX(I,2), FX(J,2), FX(I,3), FX(J,3))) THEN
           IF (.NOT. SHRINKING(I)) THEN
              SHRINKING(I) = .TRUE. ; NS = NS+1 ; TO_SHRINK(NS) = I
           END IF
           IF (.NOT. SHRINKING(J)) THEN
              SHRINKING(J) = .TRUE. ; NS = NS+1 ; TO_SHRINK(NS) = J
           END IF
        END IF
     END DO check_monotonicity
     NC = 0
  END DO
  ! ------------------------------------------------------------------
  
  ! Use "FIT_SPLINE" to fit the final PMQSI.
  CALL FIT_SPLINE(X, FX, T, BCOEF, INFO)


CONTAINS
  FUNCTION IS_MONOTONE(U0, U1, X0, X1, DX0, DX1, DDX0, DDX1)
    ! Given an interval and function values (value, first, and seocond derivative),
    ! return TRUE if the quintic piece is monotone, FALSE if it is nonmonotone.
    REAL(KIND=R8), INTENT(IN) :: U0, U1, X0, X1, DX0, DX1, DDX0, DDX1
    LOGICAL IS_MONOTONE
    ! Local variables.
    REAL(KIND=R8) :: A, ALPHA, B, BETA, GAMMA, INV_SLOPE, SIGN, TAU, TEMP, W
    ! Identify the direction of change of the function (increasing or decreasing).
    IF      (X1 .GT. X0) THEN ; SIGN =  1.0_R8
    ELSE IF (X1 .LT. X0) THEN ; SIGN = -1.0_R8
    ! When the function values are flat, everything *must* be 0.
    ELSE
       IS_MONOTONE = (DX0  .EQ. 0.0_R8) .AND. (DX1  .EQ. 0.0_R8) .AND. &
                     (DDX0 .EQ. 0.0_R8) .AND. (DDX1 .EQ. 0.0_R8)
       RETURN
    END IF
    ! Make sure the slopes point in the correct direction.
    IF (SIGN*DX0 .LT. 0.0_R8) THEN ; IS_MONOTONE = .FALSE. ; RETURN ; END IF
    IF (SIGN*DX1 .LT. 0.0_R8) THEN ; IS_MONOTONE = .FALSE. ; RETURN ; END IF
    ! Compute some useful constants related to the theory.
    W = U1 - U0
    INV_SLOPE = W / (X1 - X0)
    A = INV_SLOPE * DX0
    B = INV_SLOPE * DX1
    ! Simplified cubic monotone case.
    IF (A*B .LE. 0.0_R8) THEN
       ! The check for delta >= 0 has already been performed above,
       ! next check for alpha >= 0.
       IF (SIGN*DDX1 .GT. SIGN*4.0_R8 * DX1 / W) THEN ; IS_MONOTONE = .FALSE.
       ELSE
          ! Now compute the value of 2 * \sqrt{ alpha delta }, store in "TEMP".
          TEMP = SIGN * 2.0_R8 * INV_SLOPE * SQRT(DX0 * (4*DX1 - DDX1*W))
          ! Check for gamma >= delta - 2 * \sqrt{ alpha delta }
          IF (TEMP + INV_SLOPE*(3.0_R8*DX0 + DDX0*W) .LT. 0.0_R8) THEN
             IS_MONOTONE = .FALSE.
          ! Check for beta >= alpha - 2 * \sqrt{ alpha delta }
          ELSE IF (60.0_R8 + 2.0_R8*TEMP + INV_SLOPE*((5.0_R8*DDX1 - 3.0_R8*DDX0) * &
               W - 8.0_R8*(3.0_R8*DX0 + 4.0_R8*DX1)) .LT. 0.0_R8) THEN
             IS_MONOTONE = .FALSE.
          ELSE ; IS_MONOTONE = .TRUE.
          END IF
       END IF
    ELSE
       ! Full quintic monotone case.
       TAU = 24.0_R8 + 2.0_R8*SQRT(A*B) - 3.0_R8*(A+B)
       IF (TAU .LT. 0.0_R8) THEN ; IS_MONOTONE = .FALSE.
       ELSE
          ! Compute alpha, gamma, beta from theorems to determine monotonicity.
          TEMP = (A / B)**(0.75_R8)
          ALPHA = TEMP * (4.0_R8*DX1 - DDX1*W) / DX0
          GAMMA = (4.0_R8*DX0 + DDX0*W) / (TEMP * DX1)
          BETA = (3.0_R8 * INV_SLOPE * ((DDX1-DDX0)*W - 8.0_R8*(DX0+DX1)) + 60.0_R8) &
               / (2.0 * SQRT(A*B))
          IF (BETA .LE. 6.0_R8) THEN ; TEMP = -(BETA + 2.0_R8) / 2.0_R8
          ELSE ;                       TEMP = -2.0_R8 * SQRT(BETA - 2.0_R8)
          END IF
          IS_MONOTONE = (ALPHA .GT. TEMP) .AND. (GAMMA .GT. TEMP)
       END IF
    END IF
  END FUNCTION IS_MONOTONE

  SUBROUTINE QUADRATIC(I2)
    ! Given an index "I", compute the coefficients "A" on x^2 and "B"
    ! on x for the quadratic interpolant of X(I-1:I+1), Y(I-1:I+1).
    INTEGER, INTENT(IN) :: I2
    ! Local variables.
    REAL(KIND=R8) :: D, & ! Denominator for computing "A" and "B".
         C1, C2, C3 ! Intermediate terms for computation.
    INTEGER :: I1, I3
    I1 = I2-1
    I3 = I2+1
    ! Compute the shared denominator.
    D = (X(I1) - X(I2)) * (X(I1) - X(I3)) * (X(I2) - X(I3))
    ! Prevent division by zero.
    IF (D .EQ. 0.0_R8) THEN
       A = MAX_REAL
       B = 0.0_R8
    END IF
    ! Compute coefficients A and B in "Ax^2 + Bx + C" quadratic interpolant.
    C1 = X(I1) * (Y(I3) - Y(I2))
    C2 = X(I2) * (Y(I1) - Y(I3))
    C3 = X(I3) * (Y(I2) - Y(I1))
    A = (C1 + C2 + C3) / D
    B = - (X(I1)*C1 + X(I2)*C2 + X(I3)*C3) / D
  END SUBROUTINE QUADRATIC

END SUBROUTINE PMQSI


SUBROUTINE FIT_SPLINE(XI, FX, T, BCOEF, INFO)
  ! Subroutine for computing a linear combination of B-splines that
  ! interpolates the given function value (and function derivatives)
  ! at the given breakpoints.
  ! 
  ! INPUT:
  !   XI(1:NB) -- The increasing real-valued locations of the
  !               breakpoints for the interpolating spline.
  !   FX(1:NB,1:NCC) -- FX(I,J) contains the (J-1)st derivative
  !                     at XI(I) to be interpolated.
  ! 
  ! OUTPUT:
  !   T(1:NB*NCC+2*NCC) -- The nondecreasing real-valued locations
  !                        of the knots for the B-spline basis.
  !   BCOEF(1:NB*NCC) -- The coefficients for the B-splines
  !                      that define the interpolating spline.
  !   INFO -- Integer representing the subroutine execution status:
  !       0    Successful execution.
  !       1    SIZE(XI) is less than 1.
  !       2    SIZE(FX) is less then 1.
  !       3    SIZE(FX,1) does not equal SIZE(XI).
  !       4    SIZE(T) too small, should be at least NB*NCC + 2*NCC.
  !       5    SIZE(BCOEF) too small, should be at least NB*NCC.
  !       6    Elements of XI are not strictly increasing.
  !       7    The computed spline does not match the provided FX
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
  !   The interpolating spline is returned in terms of knots "T" and
  !   BCOEF that define the underlying B-splines and the
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
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:)   :: XI
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: FX
  ! REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: T, BCOEF
  REAL(KIND=R8), INTENT(OUT), DIMENSION(SIZE(FX)+2*SIZE(FX,2)) :: T
  REAL(KIND=R8), INTENT(OUT), DIMENSION(SIZE(FX)) :: BCOEF
  ! ^^^^ THE ABOVE 2 LINES ARE TEMPORARY, SHOULD NOT BE IN FINAL CODE.
  !      THEY ONLY EXIST TO MAKE AUTOMATIC WRAPPING EASIER.
  INTEGER, INTENT(OUT) :: INFO
  ! Local variables.
  INTEGER, DIMENSION(SIZE(BCOEF)) :: IPIV
  REAL(KIND=R8), DIMENSION(1 + 3*(2*SIZE(FX,2)-1), SIZE(FX)) :: AB
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
  NB = SIZE(XI)
  NCC = SIZE(FX,2)
  NSPL = SIZE(FX)
  NK = NSPL + 2*NCC
  K = 2*NCC
  DEGREE = K - 1
  INFO = 0

  ! Check the shape of incoming arrays.
  IF      (NB .LT. 1)             THEN ; INFO = 1 ; RETURN
  ELSE IF (NSPL .LT. 1)           THEN ; INFO = 2 ; RETURN
  ELSE IF (SIZE(FX,1) .NE. NB)    THEN ; INFO = 3 ; RETURN
  ELSE IF (SIZE(T) .LT. NK)       THEN ; INFO = 4 ; RETURN
  ELSE IF (SIZE(BCOEF) .LT. NSPL) THEN ; INFO = 5 ; RETURN
  END IF
  ! Verify that breakpoints are increasing.
  DO I = 1, NB - 1
     IF (XI(I) .GE. XI(I+1)) THEN ; INFO = 6 ; RETURN ; END IF
  END DO

  ! Copy over the knots that will define the B-spline representation.
  ! Each internal knot will be repeataed NCC times to maintain the
  ! necessary level of continuity for this spline.
  T(1:K) = XI(1)
  DO I = 2, NB-1
     T(I*NCC+1 : (I+1)*NCC) = XI(I)
  END DO
  ! Assign the last knot to exist a small step outside the supported
  ! interval to ensure the B-spline basis functions are nonzero at the
  ! rightmost breakpoint.
  T(NK-DEGREE:NK) = XI(NB) * (1.0_R8 + SQRT(EPSILON(XI(NB))))

  ! The next block of code evaluates each B-spline and it's
  ! derivatives at all breakpoints. The first and last elements of
  ! XI will be repeated K times and each internal breakpoint
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
  AB(:,:) = 0.0_R8
  ! Evaluate all B-splines at all breakpoints (walking through columns).
  DO I = 1, NSPL
     ! Compute indices of the first and last knot for the current B-spline.
     J = I + K ! Last knot.
     ! Compute the row indices in "A" that would be accessed.
     I1 = ((I-1)/NCC - 1) * NCC + 1 ! First row.
     I2 = I1 + 3*NCC - 1            ! Last row.
     ! Only two breakpoints will be covered for the first NCC
     ! B-splines and the last NCC B-splines.
     IF (I .LE. NCC)      I1 = I1 + NCC
     IF (I+NCC .GT. NSPL) I2 = I2 - NCC
     ! Compute the indices of the involved breakpoints.
     J1 = I1 / NCC + 1 ! First breakpoint.
     J2 = I2 / NCC     ! Last breakpoint.
     ! Convert the "i,j" indices in "A" to the banded storage scheme.
     ! The mapping is looks like   AB[KL+KU+1+i-j,j] = A[i,j]
     I1 = 2*DEGREE+1 + I1 - I
     I2 = 2*DEGREE+1 + I2 - I
     ! Evaluate this B-spline, computing function value and derivatives.
     DO DERIV = 0, NCC-1
        ! Place the evaluations into a block of a column in AB, shift
        ! according to which derivative is being evaluated and use a
        ! stride determined by the number of continuity conditions.
        AB(I1+DERIV:I2:NCC,I) = XI(J1:J2)
        CALL EVAL_BSPLINE(T(I:J), AB(I1+DERIV:I2:NCC,I), D=DERIV)
     END DO
  END DO
  ! Copy the FX into the BCOEF (output) variable.
  DO I = 1, NCC
     BCOEF(I:NSPL:NCC) = FX(:,I)
  END DO
  ! Call the LAPACK subroutine to solve the banded linear system.
  CALL DGBSV(NSPL, DEGREE, DEGREE, 1, AB, SIZE(AB,1), IPIV, BCOEF, NSPL, INFO)
  ! Check for errors in the execution of DGBSV, (this should not happen).
  IF (INFO .NE. 0) THEN ; INFO = INFO + 10 ; RETURN ; END IF
  ! Check to see if the linear system was correctly solved by looking at
  ! the difference between prouduced B-spline values and provided values.
  MAX_ERROR = SQRT(SQRT(EPSILON(1.0_R8)))
  DO DERIV = 0, NCC-1
     ! Reuse the first row of AB as scratch space (the first column
     ! might not be large enough, but the first row certainly is).
     AB(1,1:NB) = XI(:)
     ! Evaluate this spline at all breakpoints. Correct usage is
     ! enforced here, so it is expected that INFO=0 always.
     CALL EVAL_SPLINE(T, BCOEF, AB(1,1:NB), INFO, D=DERIV)
     ! ^ Assumed shape "XY" argument of "EVAL_SPLINE" can be given a
     !   stride by its implicit array descriptor, no copy is necessary.
     ! Check the precision of the reproduced values.
     ! Return an error if the precision is too low.
     IF (MAXVAL(ABS((AB(1,1:NB) - FX(:,DERIV+1)) &
          / (1.0_R8 + ABS(FX(:,DERIV+1))))) .GT. MAX_ERROR) THEN
        INFO = 7
        RETURN
     END IF
  END DO
END SUBROUTINE FIT_SPLINE


SUBROUTINE EVAL_SPLINE(T, BCOEF, XY, INFO, D)
  ! Evaluate a spline construced with FIT_SPLINE. Points to estimate
  ! provided in XY, evaluates D derivative at all XY, result in XY.
  !   
  ! INPUT:
  !   T(1:N+2*NCC) -- The nondecreasing real-valued locations of
  !                       the knots for the underlying B-splines, where
  !                       the spline has NCC-1 continuous derivatives.
  !   BCOEF(1:N) -- The coefficients assigned to each B-spline
  !                 that underpins this interpolating spline.
  ! 
  ! INPUT / OUTPUT:
  !   XY(1:Z) -- The locations at which the spline is evaluated on
  !              input, on output holds the value of the spline with
  !              knots T and BCOEF evaluated at the given locations.
  ! 
  ! OUTPUT:
  !   INFO -- Integer representing subroutine exeuction status.
  !       0    Successful execution.
  !       1    Extrapolation error, some X are outside of spline support.
  !       2    T contains at least one decreasing interval.
  !       3    T has size less than or equal to 1.
  !       4    T has an empty interior (T(1) = T(N+2*C)).
  !       5    Invalid BCOEF, size smaller than or equal to T.
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
  REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: T, BCOEF
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
  INTEGER, INTENT(IN), OPTIONAL :: D
  INTEGER, INTENT(OUT) :: INFO
  ! Local variables.
  INTEGER :: DERIV, I, K, NT, NB
  REAL(KIND=R8), DIMENSION(SIZE(XY), SIZE(BCOEF)) :: BIATX

  NT = SIZE(T)
  NB = SIZE(BCOEF)
  ! Check for size-related errors and for an empty knot interior.
  IF      (NT .LE. 1)       THEN ; INFO = 3 ; RETURN
  ELSE IF (T(1) .EQ. T(NT)) THEN ; INFO = 4 ; RETURN
  ELSE IF (NT .LE. NB)      THEN ; INFO = 5 ; RETURN
  ! Check for extrapolation errors.
  ELSE IF ((MINVAL(XY(:)) .LT. T(1)) .OR. &
       (MAXVAL(XY(:)) .GE. T(NT))) THEN ; INFO = 1 ; RETURN
  END IF
  ! Check for valid (nondecreasing) knot sequence.
  DO I = 1, NT-1
     IF (T(I) .GT. T(I+1)) THEN ; INFO = 2 ; RETURN ; END IF
  END DO
  ! Assign the local value of the optional derivative "D" argument.
  set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
  ELSE ; DERIV = 0
  END IF set_derivative
  ! Compute the order (number of knots minus one) for each B-spline.
  K = NT - NB
  ! Evaluate all splines at all the X positions.
  DO I = 1, NB
     ! Copy potentially strided "XY" into contiguous block of memory.
     BIATX(:,I) = XY(:)
     CALL EVAL_BSPLINE(T(I:I+K), BIATX(:,I), D=DERIV)
     ! ^ Contiguous knots "T", as well as a contiguous "BIATX" along
     !   columns ensures no copy will be necessary for this call.
  END DO
  ! Store the values into Y as the weighted sums of B-spline evaluations.
  XY(:) = MATMUL(BIATX(:,:), BCOEF(:))
  ! Set "INFO" to a successful execution value.
  INFO = 0
END SUBROUTINE EVAL_SPLINE


SUBROUTINE EVAL_BSPLINE(T, XY, D)
  ! Subroutine for evaluating a B-spline with provided knot sequence.
  ! 
  ! INPUT:
  !   T(1:N) -- The nondecreasing sequence of knots for the B-spline.
  ! 
  ! INPUT / OUTPUT:
  !   XY(1:Z) -- The locations at which the B-spline is evaluated on
  !              input, on output holds the value of the B-spline with
  !              prescribed knots evaluated at the given X locations.
  ! 
  ! OPTIONAL INPUT:
  !   D [= 0]  --  The derivative to take of the evaluated B-spline.
  !                When negative, this subroutine integrates the B-spline.
  ! 
  ! DESCRIPTION:
  ! 
  !    This function uses the recurrence relation defining a B-spline:
  ! 
  !      B_{I,1}(X)   =   1     if T(I) <= X < T(I+1),
  !                       0     otherwise,
  ! 
  !    where I is the knot index, J = 2, ..., N-MAX(D,0)-1, and
  ! 
  !                        X - T(I)                         T(I+J) - X                    
  !      B_{I,J}(X) =  ----------------- B_{I,J-1}(X) +  ----------------- B_{I+1,J-1}(X).
  !                     T(I+J-1) - T(I)                   T(I+J) - T(I+1)                 
  !                                                                   
  !    All of the intermediate steps (J) are stored in a single block
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
  !    builds from J = N-D, ..., N-1 as
  ! 
  !                     (J-1) B_{I,J-1}(X)       (J-1) B_{I+1,J-1}(X)  
  !      B_{I,J}(X) =  --------------------  -  ----------------------.
  !                      T(I+J-1) - T(I)           T(I+J) - T(I+1)     
  !                                                   
  !     The final B-spline is right continuous and has support over
  !     the interval [T(1), T(N)).
  ! 
  REAL(KIND=R8), INTENT(IN),    DIMENSION(:) :: T
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
  INTEGER, INTENT(IN), OPTIONAL :: D
  ! Local variables.
  REAL(KIND=R8), DIMENSION(SIZE(XY), SIZE(T)) :: BIATX
  INTEGER :: DERIV, I, J, K, L, N
  REAL(KIND=R8) :: LEFT, RIGHT, TN
  ! Assign the local value of the optional derivative "D" argument.
  set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
  ELSE ; DERIV = 0
  END IF set_derivative
  ! Collect local variables for notational convenience.
  N = SIZE(T) ! Number of knots
  K = N - 1 ! Order of B-spline
  L = K - 1 ! One less than the order of the B-spline.
  TN = T(N) ! Value of the last knot, T(N).

  ! If this is a large enough derivative, we know it is zero everywhere.
  IF (DERIV+1 .GE. N) THEN ; XY(:) = 0.0_R8 ; RETURN

  ! ---------------- Performing standard evaluation ------------------
  ! This is a standard B-spline with multiple unique knots, right continuous.
  ELSE
     ! Initialize all values to 0.
     BIATX(:,:) = 0.0_R8
     ! Assign the first value for each knot index.
     first_b_spline : DO J = 1, K
        IF (T(J) .EQ. T(J+1)) CYCLE first_b_spline
        ! Compute all right continuous order 1 B-spline values.
        WHERE ( (T(J) .LE. XY(:)) .AND. (XY(:) .LT. T(J+1)) )
           BIATX(:,J) = 1.0_R8
        END WHERE
     END DO first_b_spline
  END IF
  ! Compute the remainder of B-spline by building up from the first.
  ! Omit the final steps of this computation for derivatives.
  compute_spline : DO I = 2, K-MAX(DERIV,0)
     ! Cycle over each knot accumulating the values for the recurrence.
     DO J = 1, N - I
        ! Check divisors, intervals with 0 width add 0 value to the B-spline.
        LEFT = (T(J+I-1) - T(J))
        RIGHT = (T(J+I) - T(J+1))
        ! Compute the B-spline recurrence relation (cases based on divisor).
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = &
                   ((XY(:) - T(J)) / LEFT) * BIATX(:,J) + &
                   ((T(J+I) - XY(:)) / RIGHT) * BIATX(:,J+1)
           ELSE
              BIATX(:,J) = ((XY(:) - T(J)) / LEFT) * BIATX(:,J)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = ((T(J+I) - XY(:)) / RIGHT) * BIATX(:,J+1)
           END IF
        END IF
     END DO
  END DO compute_spline

  ! -------------------- Performing integration ----------------------
  int_or_diff : IF (DERIV .LT. 0) THEN
  ! Integrals will be nonzero on [TN, \infty).
  WHERE (TN .LE. XY(:))
     BIATX(:,N) = 1.0_R8
  END WHERE
  ! Loop through starting at the back, raising the order of all
  ! constituents to match the order of the first.
  raise_order : DO I = 1, L
     DO J = N-I, K
        LEFT = (TN - T(J))
        RIGHT = (TN - T(J+1))
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = &
                   ((XY(:) - T(J)) / LEFT) * BIATX(:,J) + &
                   ((TN - XY(:)) / RIGHT) * BIATX(:,J+1)
           ELSE
              BIATX(:,J) = ((XY(:) - T(J)) / LEFT) * BIATX(:,J)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = ((TN - XY(:)) / RIGHT) * BIATX(:,J+1)
           END IF
        END IF
     END DO
  END DO raise_order
  
  ! Compute the integral(s) of the B-spline.
  compute_integral : DO I = 1, -DERIV
     ! Do a forward evaluation of all constituents.
     DO J = 1, K
        LEFT = (TN - T(J))
        RIGHT = (TN - T(J+1))
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = &
                   ((XY(:) - T(J)) / LEFT) * BIATX(:,J) + &
                   ((TN - XY(:)) / RIGHT) * BIATX(:,J+1)
           ELSE
              BIATX(:,J) = ((XY(:) - T(J)) / LEFT) * BIATX(:,J)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = ((TN - XY(:)) / RIGHT) * BIATX(:,J+1)
           END IF
        END IF
     END DO
     ! Sum the constituent functions at each knot (from the back).
     DO J = K, 1, -1
        BIATX(:,J) = (BIATX(:,J) + BIATX(:,J+1))
     END DO
     ! Divide by the degree plus the integration coefficient.
     BIATX(:,1) = BIATX(:,1) / (L+I)
     ! Rescale then integral by its width.
     BIATX(:,1) = BIATX(:,1) * (TN - T(1))
     ! Extend the previous two computations if more integrals need
     ! to be computed after this one.
     IF (I+DERIV .LT. 0) THEN
        BIATX(:,2:N) = BIATX(:,2:N) / (L+I)
        DO J = 2, N
           BIATX(:,J) = BIATX(:,J) * (TN - T(J))
        END DO
     END IF
  END DO compute_integral

  ! ------------------ Performing differentiation --------------------
  ELSE IF (DERIV .GT. 0) THEN
  ! Compute the derivative of the B-spline (if D > 0).
  compute_derivative : DO I = N-DERIV, K
     ! Cycle over each knot, following the same structure with the
     ! derivative computing relation instead of the B-spline one.
     DO J = 1, N-I
        ! Assure that the divisor will not cause invalid computations.
        LEFT = (T(J+I-1) - T(J))
        RIGHT = (T(J+I) - T(J+1))
        ! Compute the derivative recurrence relation.
        IF (LEFT .GT. 0) THEN
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) =  (I-1)*(BIATX(:,J)/LEFT - BIATX(:,J+1)/RIGHT)
           ELSE
              BIATX(:,J) = (I-1)*(BIATX(:,J)/LEFT)
           END IF
        ELSE
           IF (RIGHT .GT. 0) THEN
              BIATX(:,J) = (I-1)*(-BIATX(:,J+1)/RIGHT)
           END IF
        END IF
     END DO
  END DO compute_derivative
  END IF int_or_diff

  ! Assign the values into the "Y" output.
  XY(:) = BIATX(:,1)
END SUBROUTINE EVAL_BSPLINE

END MODULE SPLINES
