SUBROUTINE PMQSI(X, Y, T, BCOEF, INFO)
! Computes a piecewise monotone quintic spline interpolant (PMQSI) Q(X) to
! data in terms of spline coefficients BCOEF for a B-spline basis defined
! by knots T.  Q(X) is theoretically guaranteed to be monotone increasing
! (decreasing) over exactly the same intervals [X(I), X(J)] that the data
! Y(.) is monotone increasing (decreasing).
! 
! INPUT:
!   X(1:ND) -- A real array of increasing values.
!   Y(1:ND) -- A real array of function values Y(I) = F(X(I)).
! 
! OUTPUT:
!   T(1:3*ND+6) -- The knots for the PMQSI B-spline basis.
!   BCOEF(1:3*ND) -- The coefficients for the PMQSI B-spline basis.
!   INFO -- Integer representing the subroutine return status.
!    0  Normal return.
!    1  There are fewer than three data points, so there is nothing to do.
!    2  X(:) and Y(:) have different sizes.
!    3  The values in X(:) are not increasing or not separated by at least
!       the machine precision.
!    4  The size of T(:) must be at least the number of knots 3*ND+6.
!    5  The size of BCOEF(:) must be at least the spline space dimension 3*ND.
!    7  The computed spline does not match the provided data in Y(:) and
!       this result should not be used. This arises when the scaling of
!       function values and derivative values causes the linear system used
!       to compute the final spline interpolant to have a prohibitively
!       large condition number.
!  >10  10 plus the info flag as returned by DGBSV from LAPACK when computing
!       the final spline interpolant.
!   
USE REAL_PRECISION, ONLY: R8
IMPLICIT NONE
REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X, Y
REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: T, BCOEF
INTEGER, INTENT(OUT) :: INFO
! 
! Local variables.
REAL(KIND=R8), DIMENSION(SIZE(X),2) :: FHATX ! Estimated function values.
REAL(KIND=R8), DIMENSION(SIZE(X),3) :: FX ! Spline function values.
LOGICAL, DIMENSION(SIZE(X)) :: CHECKING, GROWING, SHRINKING ! Execution
INTEGER, DIMENSION(SIZE(X)) :: TO_CHECK, TO_GROW, TO_SHRINK !   queues.
REAL(KIND=R8) :: A, ACCURACY, B, DIRECTION, DX, REAL_MAX, STEP_SIZE
INTEGER :: I, J, K, NC, ND, NG, NS
LOGICAL :: SEARCHING
INTERFACE
   SUBROUTINE FIT_SPLINE(XI, FX, T, BCOEF, INFO)
     USE REAL_PRECISION, ONLY: R8
     REAL(KIND=R8), INTENT(IN),  DIMENSION(:)   :: XI
     REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: FX
     REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: T, BCOEF
     INTEGER, INTENT(OUT) :: INFO
   END SUBROUTINE FIT_SPLINE
END INTERFACE
! Declare constants for computation.
ND = SIZE(X)
REAL_MAX = HUGE(1.0_R8)
ACCURACY = SQRT(EPSILON(1.0_R8))
! Check the shape of incoming arrays.
IF      (ND .LT. 3)             THEN; INFO = 1; RETURN
ELSE IF (SIZE(Y) .NE. ND)       THEN; INFO = 2; RETURN
ELSE IF (SIZE(T) .LT. 3*ND + 6) THEN; INFO = 4; RETURN
ELSE IF (SIZE(BCOEF) .LT. 3*ND) THEN; INFO = 5; RETURN
END IF
! Verify that X values are increasing and appropriately spaced.
DO I = 1, ND-1
   IF (X(I+1) - X(I) .LE. ACCURACY) THEN; INFO = 3; RETURN; END IF
END DO
! ==================================================================
!            Estimate initial derivatives by using a
!            minimum curvature quadratic facet model.
! 
! Copy the Y values into the FX array.
FX(:,1) = Y(:)
! Identify local extreme points and flat points. Denote location
! of flats in GROWING and extrema in SHRINKING.
GROWING(1) = Y(1) .EQ. Y(2)
GROWING(ND) = Y(ND-1) .EQ. Y(ND)
SHRINKING(1) = .FALSE.
SHRINKING(ND) = .FALSE.
DO I = 2, ND-1
   DIRECTION = (Y(I) - Y(I-1)) * (Y(I+1) - Y(I))
   GROWING(I) = DIRECTION .EQ. 0.0_R8
   SHRINKING(I) = DIRECTION .LT. 0.0_R8
END DO
! Use local quadratic interpolants to estimate slopes and second
! derivatives at all points. Use zero-sloped quadratic interpolants
! at extrema with left/right neighbors to estimate curvature.
estimate_derivatives : DO I = 1, ND
   ! Precompute these frequently used integers.
   J = I+1
   K = I-1
   ! Determine the direction of change at the point I.
   IF (GROWING(I) .OR. SHRINKING(I)) THEN ; DIRECTION = 0.0_R8
   ELSE IF (I .EQ. 1) THEN
      IF      (Y(I) .LT. Y(J)) THEN; DIRECTION =  1.0_R8
      ELSE IF (Y(I) .GT. Y(J)) THEN; DIRECTION = -1.0_R8
      END IF
   ELSE
      IF      (Y(K) .LT. Y(I)) THEN; DIRECTION =  1.0_R8
      ELSE IF (Y(K) .GT. Y(I)) THEN; DIRECTION = -1.0_R8
      END IF
   END IF
   ! Initialize the curvature to be maximally large.
   FX(I,3) = REAL_MAX
   ! If this is a local flat, first and second derivatives are zero valued.
   pick_quadratic : IF (GROWING(I)) THEN; FX(I,2:3) = 0.0_R8
   ! If this is an extreme point, construct quadratic interpolants
   ! that have zero slope here and hit left/right neighbors.
   ELSE IF (SHRINKING(I)) THEN
      ! Set the first derivative to zero.
      FX(I,2) = 0.0_R8
      ! Compute the coefficient A in  Ax^2 + Bx + C  that interpolates X(I-1).
      FX(I,3) = (Y(K) - Y(I)) / (X(K) - X(I))**2
      ! Compute the coefficient A in  Ax^2 + Bx + C  that interpolates X(I+1).
      A = (Y(J) - Y(I)) / (X(J) - X(I))**2
      IF (ABS(A) .LT. ABS(FX(I,3))) THEN; FX(I,3) = A; END IF
      ! Compute the actual second derivative (instead of coefficient A).
      FX(I,3) = MAX(MIN(2.0_R8 * FX(I,3), REAL_MAX), -REAL_MAX)
   ELSE
      ! --------------------
      ! Quadratic left of I.
      IF (K .GT. 1) THEN
         ! If there's a zero-derivative to the left, use it's right interpolant.
         IF (SHRINKING(K) .OR. GROWING(K)) THEN
            A = (Y(I) - Y(K)) / (X(I) - X(K))**2
            B = -2.0_R8 * X(K) * A
         ! Otherwise use the standard quadratic on the left.
         ELSE; CALL QUADRATIC(K)
         END IF
         DX = MAX(MIN(2.0_R8*A*X(I) + B, REAL_MAX), -REAL_MAX)
         IF (DX*DIRECTION .GE. 0.0_R8) THEN
            FX(I,2) = DX
            FX(I,3) = A
         END IF
      END IF
      ! ------------------------
      ! Quadratic centered at I (require that it has at least one
      ! neighbor that is not forced to zero slope).
      IF ((I .GT. 1) .AND. (I .LT. ND) .AND. .NOT. &
           ((SHRINKING(K) .OR. GROWING(K)) .AND. &
           (SHRINKING(J) .OR. GROWING(J)))) THEN
         ! Construct the quadratic interpolant through this point and neighbors.
         CALL QUADRATIC(I)
         DX = MAX(MIN(2.0_R8*A*X(I) + B, REAL_MAX), -REAL_MAX)
         ! Keep this new quadratic if it has less curvature.
         IF ((DX*DIRECTION .GE. 0.0_R8) .AND. (ABS(A) .LT. ABS(FX(I,3)))) THEN
            FX(I,2) = DX
            FX(I,3) = A
         END IF
      END IF
      ! ---------------------
      ! Quadratic right of I.
      IF (J .LT. ND) THEN
         ! If there's an extreme point to the right, use it's left interpolant.
         IF (SHRINKING(J) .OR. GROWING(J)) THEN
            A = (Y(I) - Y(J)) / (X(I) - X(J))**2
            B = -2.0_R8 * X(J) * A
         ! Otherwise use the standard quadratic on the right.
         ELSE; CALL QUADRATIC(J)
         END IF
         DX = MAX(MIN(2.0_R8*A*X(I) + B, REAL_MAX), -REAL_MAX)
         ! Keep this new quadratic if it has less curvature.
         IF ((DX*DIRECTION .GE. 0.0_R8) .AND. (ABS(A) .LT. ABS(FX(I,3)))) THEN
            FX(I,2) = DX
            FX(I,3) = A
         END IF
      END IF
      ! Set the final quadratic.
      IF (FX(I,3) .EQ. REAL_MAX) THEN; FX(I,2:3) = 0.0_R8
      ! Compute the actual curvature of the quadratic, instead of coefficient on x^2.        
      ELSE; FX(I,3) = MAX(MIN(2.0_R8 * FX(I,3), REAL_MAX), -REAL_MAX)
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
NC = 0; NG = 0; NS = 0
DO I = 1, ND-1
   J = I+1
   IF (.NOT. IS_MONOTONE(X(I), X(J), FX(I,1), FX(J,1), &
        FX(I,2), FX(J,2), FX(I,3), FX(J,3))) THEN
      IF (.NOT. SHRINKING(I)) THEN
         SHRINKING(I) = .TRUE.; NS = NS+1; TO_SHRINK(NS) = I
      END IF
      IF (.NOT. SHRINKING(J)) THEN
         SHRINKING(J) = .TRUE.; NS = NS+1; TO_SHRINK(NS) = J
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
         SEARCHING = .FALSE.; STEP_SIZE = ACCURACY; NG = 0
      END IF
   ! Grow the step size (at a slower rate) if there are still intervals to fix.
   ELSE; STEP_SIZE = STEP_SIZE + STEP_SIZE / 2.0_R8
   END IF
   ! Grow all those first and second derivatives that were previously 
   ! shrunk, but are currently monotone.
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
         CHECKING(I-1) = .TRUE.; NC = NC+1; TO_CHECK(NC) = I-1
      END IF
      IF ((I .LT. ND) .AND. (.NOT. CHECKING(I))) THEN
         CHECKING(I) = .TRUE.; NC = NC+1; TO_CHECK(NC) = I
      END IF
   END DO grow_values
   ! Shrink all those first and second derivatives that cause nonmonotonicity.
   shrink_values : DO J = 1, NS
      I = TO_SHRINK(J)
      SHRINKING(I) = .FALSE.
      IF (SEARCHING .AND. (.NOT. GROWING(I))) THEN
         GROWING(I) = .TRUE.; NG = NG+1; TO_GROW(NG) = I
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
         CHECKING(I-1) = .TRUE.; NC = NC+1; TO_CHECK(NC) = I-1
      END IF
      IF ((I .LT. ND) .AND. (.NOT. CHECKING(I))) THEN
         CHECKING(I) = .TRUE.; NC = NC+1; TO_CHECK(NC) = I
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
            SHRINKING(I) = .TRUE.; NS = NS+1; TO_SHRINK(NS) = I
         END IF
         IF (.NOT. SHRINKING(J)) THEN
            SHRINKING(J) = .TRUE.; NS = NS+1; TO_SHRINK(NS) = J
         END IF
      END IF
   END DO check_monotonicity
   NC = 0
END DO
! ------------------------------------------------------------------
! Use FIT_SPLINE to fit the final PMQSI.
CALL FIT_SPLINE(X, FX, T, BCOEF, INFO)

CONTAINS
FUNCTION IS_MONOTONE(U0, U1, X0, X1, DX0, DX1, DDX0, DDX1)
  ! Given an interval and function values (value, first, and seocond derivative),
  ! return TRUE if the quintic piece is monotone, FALSE if it is nonmonotone.
  REAL(KIND=R8), INTENT(IN) :: U0, U1, X0, X1, DX0, DX1, DDX0, DDX1
  LOGICAL :: IS_MONOTONE
  ! Local variables.
  REAL(KIND=R8) :: A, ALPHA, B, BETA, GAMMA, INV_SLOPE, SIGN, TAU, TEMP, W
  ! Identify the direction of change of the function (increasing or decreasing).
  IF      (X1 .GT. X0) THEN; SIGN =  1.0_R8
  ELSE IF (X1 .LT. X0) THEN; SIGN = -1.0_R8
  ! When the function values are flat, everything *must* be 0.
  ELSE
     IS_MONOTONE = (DX0  .EQ. 0.0_R8) .AND. (DX1  .EQ. 0.0_R8) .AND. &
                   (DDX0 .EQ. 0.0_R8) .AND. (DDX1 .EQ. 0.0_R8)
     RETURN
  END IF
  ! Make sure the slopes point in the correct direction.
  IF (SIGN*DX0 .LT. 0.0_R8) THEN; IS_MONOTONE = .FALSE.; RETURN; END IF
  IF (SIGN*DX1 .LT. 0.0_R8) THEN; IS_MONOTONE = .FALSE.; RETURN; END IF
  ! Compute some useful constants related to the theory.
  W = U1 - U0
  INV_SLOPE = W / (X1 - X0)
  A = INV_SLOPE * DX0
  B = INV_SLOPE * DX1
  ! Simplified cubic monotone case.
  IF (A*B .LE. 0.0_R8) THEN
     ! The check for delta >= 0 has already been performed above,
     ! next check for alpha >= 0.
     IF (SIGN*DDX1 .GT. SIGN*4.0_R8 * DX1 / W) THEN; IS_MONOTONE = .FALSE.
     ELSE
        ! Now compute the value of 2 * \sqrt{ alpha delta }, store in TEMP.
        TEMP = SIGN * 2.0_R8 * INV_SLOPE * SQRT(DX0 * (4*DX1 - DDX1*W))
        ! Check for gamma >= delta - 2 * \sqrt{ alpha delta }
        IF (TEMP + INV_SLOPE*(3.0_R8*DX0 + DDX0*W) .LT. 0.0_R8) THEN
           IS_MONOTONE = .FALSE.
        ! Check for beta >= alpha - 2 * \sqrt{ alpha delta }
        ELSE IF (60.0_R8 + 2.0_R8*TEMP + INV_SLOPE*((5.0_R8*DDX1 - 3.0_R8*DDX0) * &
             W - 8.0_R8*(3.0_R8*DX0 + 4.0_R8*DX1)) .LT. 0.0_R8) THEN
           IS_MONOTONE = .FALSE.
        ELSE; IS_MONOTONE = .TRUE.
        END IF
     END IF
  ELSE
     ! Full quintic monotone case.
     TAU = 24.0_R8 + 2.0_R8*SQRT(A*B) - 3.0_R8*(A+B)
     IF (TAU .LT. 0.0_R8) THEN; IS_MONOTONE = .FALSE.
     ELSE
        ! Compute alpha, gamma, beta from theorems to determine monotonicity.
        TEMP = (A / B)**(0.75_R8)
        ALPHA = TEMP * (4.0_R8*DX1 - DDX1*W) / DX0
        GAMMA = (4.0_R8*DX0 + DDX0*W) / (TEMP * DX1)
        BETA = (3.0_R8 * INV_SLOPE * ((DDX1-DDX0)*W - 8.0_R8*(DX0+DX1)) + 60.0_R8) &
             / (2.0 * SQRT(A*B))
        IF (BETA .LE. 6.0_R8) THEN; TEMP = -(BETA + 2.0_R8) / 2.0_R8
        ELSE;                       TEMP = -2.0_R8 * SQRT(BETA - 2.0_R8)
        END IF
        IS_MONOTONE = (ALPHA .GT. TEMP) .AND. (GAMMA .GT. TEMP)
     END IF
  END IF
END FUNCTION IS_MONOTONE

SUBROUTINE QUADRATIC(I2)
  ! Given an index I, compute the coefficients A on x^2 and B
  ! on x for the quadratic interpolant of X(I-1:I+1), Y(I-1:I+1).
  INTEGER, INTENT(IN) :: I2
  ! Local variables.
  REAL(KIND=R8) :: D, & ! Denominator for computing A and B.
       C1, C2, C3 ! Intermediate terms for computation.
  INTEGER :: I1, I3
  I1 = I2-1
  I3 = I2+1
  ! Compute the shared denominator.
  D = (X(I1) - X(I2)) * (X(I1) - X(I3)) * (X(I2) - X(I3))
  ! Prevent division by zero.
  IF (D .EQ. 0.0_R8) THEN
     A = REAL_MAX
     B = 0.0_R8
  END IF
  ! Compute coefficients A and B in  Ax^2 + Bx + C  quadratic interpolant.
  C1 = X(I1) * (Y(I3) - Y(I2))
  C2 = X(I2) * (Y(I1) - Y(I3))
  C3 = X(I3) * (Y(I2) - Y(I1))
  A = (C1 + C2 + C3) / D
  B = - (X(I1)*C1 + X(I2)*C2 + X(I3)*C3) / D
END SUBROUTINE QUADRATIC

END SUBROUTINE PMQSI
