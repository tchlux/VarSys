
SUBROUTINE MQSI(X, Y, T, BCOEF, INFO)
! Computes a monotone quintic spline interpolant (MQSI) Q(X) to data
! in terms of spline coefficients BCOEF for a B-spline basis defined
! by knots T.  Q(X) is theoretically guaranteed to be monotone
! increasing (decreasing) over exactly the same intervals [X(I), X(J)]
! that the data Y(.) is monotone increasing (decreasing).
! 
! INPUT:
!   X(1:ND) -- A real array of increasing values.
!   Y(1:ND) -- A real array of function values Y(I) = F(X(I)).
! 
! OUTPUT:
!   T(1:3*ND+6) -- The knots for the MQSI B-spline basis.
!   BCOEF(1:3*ND) -- The coefficients for the MQSI B-spline basis.
!   INFO -- Integer representing the subroutine return status.
!    0  Normal return.
!    1  There are fewer than three data points, so there is nothing to do.
!    2  X(:) and Y(:) have different sizes.
!    3  The size of T(:) must be at least the number of knots 3*ND+6.
!    4  The size of BCOEF(:) must be at least the spline space dimension 3*ND.
!    5  The values in X(:) are not increasing or not separated by at least
!       their size times the square root of machine precision.
!    6  The values in X(:) are larger than the fourth root of the
!       largest representable real number.
!    7  The computed spline does not match the provided data in Y(:) and
!       this result should not be used. This arises when the scaling of
!       function values and derivative values causes the linear system used
!       to compute the final spline interpolant to have a prohibitively
!       large condition number.
!    8  The values in X(:) are separated by more than the fourth root
!       of the largest representable real number.
!    9  The values in Y(:) are separated by more than the fourth root
!       of the largest representable real number.
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
REAL(KIND=R8) :: A, ACCURACY, B, DIRECTION, DX, REAL_MAX, SEP_MAX, STEP_SIZE
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
ACCURACY = SQRT(EPSILON(1.0_R8))
REAL_MAX = HUGE(1.0_R8)
SEP_MAX = SQRT(SQRT(REAL_MAX))
! Check the shape of incoming arrays.
IF      (ND .LT. 3)             THEN; INFO = 1; RETURN
ELSE IF (SIZE(Y) .NE. ND)       THEN; INFO = 2; RETURN
ELSE IF (SIZE(T) .LT. 3*ND + 6) THEN; INFO = 3; RETURN
ELSE IF (SIZE(BCOEF) .LT. 3*ND) THEN; INFO = 4; RETURN
END IF
! Verify that X values are increasing and appropriately spaced.
DO I = 1, ND-1
   J = I + 1
   A = X(J) - X(I)
   IF (ABS(A) .LT. ACCURACY)               THEN; INFO = 5; RETURN
   ELSE IF (ABS(X(I)) .GE. SEP_MAX)        THEN; INFO = 6; RETURN
   ELSE IF (A .GE. SEP_MAX)                THEN; INFO = 8; RETURN
   ELSE IF (ABS(Y(J) - Y(I)) .GT. SEP_MAX) THEN; INFO = 9; RETURN
   END IF
END DO
! ==================================================================
!            Estimate initial derivatives by using a
!            minimum curvature quadratic facet model.
! 
! Copy the Y values into the FX array.
FX(:,1) = Y(:)
! Identify local extreme points and flat points. Denote location
! of flats in GROWING and extrema in SHRINKING.
GROWING(1) = ABS(Y(1) - Y(2)) .LT. ACCURACY
GROWING(ND) = ABS(Y(ND-1) - Y(ND)) .LT. ACCURACY
SHRINKING(1) = .FALSE.
SHRINKING(ND) = .FALSE.
DO I = 2, ND-1
   J = I+1
   K = I-1
   IF ((ABS(Y(K)-Y(I)) .LT. ACCURACY) .OR. (ABS(Y(I)-Y(J)) .LT. ACCURACY)) THEN
      GROWING(I) = .TRUE.
      SHRINKING(I) = .FALSE.
   ELSE
      SHRINKING(I) = (Y(I) - Y(K)) * (Y(J) - Y(I)) .LT. 0.0_R8
   END IF
END DO
! Use local quadratic interpolants to estimate slopes and second
! derivatives at all points. Use zero-sloped quadratic interpolants
! at extrema with left/right neighbors to estimate curvature.
estimate_derivatives : DO I = 1, ND
   ! Precompute these frequently used integers.
   J = I+1
   K = I-1
   ! Initialize the curvature to be maximally large.
   FX(I,3) = REAL_MAX
   ! If this is a local flat, first and second derivatives are zero valued.
   pick_quadratic : IF (GROWING(I)) THEN; FX(I,2:3) = 0.0_R8
   ! If this is an extreme point, construct quadratic interpolants
   ! that have zero slope here and hit left/right neighbors.
   ELSE IF (SHRINKING(I)) THEN
      ! Set the first derivative to zero.
      FX(I,2) = 0.0_R8
      ! Compute the coefficient A in  Ax^2+Bx+C  that interpolates at X(I-1).
      FX(I,3) = (Y(K) - Y(I)) / (X(K) - X(I))**2
      ! Compute the coefficient A in  Ax^2+Bx+C  that interpolates at X(I+1).
      A = (Y(J) - Y(I)) / (X(J) - X(I))**2
      IF (ABS(A) .LT. ABS(FX(I,3))) THEN; FX(I,3) = A; END IF
      ! Compute the actual second derivative (instead of coefficient A).
      FX(I,3) = 2.0_R8 * FX(I,3)
   ELSE
      ! Determine the direction of change at the point I.
      IF (I .EQ. 1) THEN
         IF      (Y(I) .LT. Y(J)) THEN; DIRECTION =  1.0_R8
         ELSE IF (Y(I) .GT. Y(J)) THEN; DIRECTION = -1.0_R8
         END IF
      ELSE
         IF      (Y(K) .LT. Y(I)) THEN; DIRECTION =  1.0_R8
         ELSE IF (Y(K) .GT. Y(I)) THEN; DIRECTION = -1.0_R8
         END IF
      END IF
      ! --------------------
      ! Quadratic left of I.
      IF (K .GT. 1) THEN
         ! If a zero derivative at left point, use its right interpolant.
         IF (SHRINKING(K) .OR. GROWING(K)) THEN
            A = (Y(I) - Y(K)) / (X(I) - X(K))**2
            B = -2.0_R8 * X(K) * A
         ! Otherwise use the standard quadratic on the left.
         ELSE; CALL QUADRATIC(K)
         END IF
         DX = 2.0_R8*A*X(I) + B
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
         ! Construct quadratic interpolant through this point and neighbors.
         CALL QUADRATIC(I)
         DX = 2.0_R8*A*X(I) + B
         ! Keep this new quadratic if it has less curvature.
         IF ((DX*DIRECTION .GE. 0.0_R8) .AND. (ABS(A) .LT. ABS(FX(I,3)))) THEN
            FX(I,2) = DX
            FX(I,3) = A
         END IF
      END IF
      ! ---------------------
      ! Quadratic right of I.
      IF (J .LT. ND) THEN
         ! If a zero derivative at right point, use its left interpolant.
         IF (SHRINKING(J) .OR. GROWING(J)) THEN
            A = (Y(I) - Y(J)) / (X(I) - X(J))**2
            B = -2.0_R8 * X(J) * A
         ! Otherwise use the standard quadratic on the right.
         ELSE; CALL QUADRATIC(J)
         END IF
         DX = 2.0_R8*A*X(I) + B
         ! Keep this new quadratic if it has less curvature.
         IF ((DX*DIRECTION .GE. 0.0_R8) .AND. (ABS(A) .LT. ABS(FX(I,3)))) THEN
            FX(I,2) = DX
            FX(I,3) = A
         END IF
      END IF
      ! Set the final quadratic.
      IF (FX(I,3) .EQ. REAL_MAX) THEN; FX(I,2:3) = 0.0_R8
      ! Compute curvature of quadratic from coefficient of x^2.
      ELSE; FX(I,3) = 2.0_R8 * FX(I,3)
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
! Identify which spline segments are not monotone and need to be fixed.
CHECKING(:) = .FALSE.
SHRINKING(:) = .FALSE.
NC = 0; NG = 0; NS = 0
DO I = 1, ND-1
   J = I+1
   ! Check for monotonicity on all segments that are not flat.
   IF ((.NOT. (GROWING(I) .AND. GROWING(J))) .AND. &
        (.NOT. IS_MONOTONE(X(I), X(J), FX(I,1), FX(J,1), &
        FX(I,2), FX(J,2), FX(I,3), FX(J,3)))) THEN
      ! Record points bounding nonomonotone segments in the TO_SHRINK queue.
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
GROWING(:) = .FALSE.
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
   ! shrunk, and correspond to currently monotone spline pieces.
   grow_values : DO J = 1, NG
      I = TO_GROW(J)
      ! Do not grow values that are actively related to a nonmonotone
      ! spline segment.
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
      ! Set this point and its neighboring intervals to be checked for
      ! monotonicity.
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
      ! Set this point and its neighboring intervals to be checked for
      ! monotonicity.
      IF ((I .GT. 1) .AND. (.NOT. CHECKING(I-1))) THEN
         CHECKING(I-1) = .TRUE.; NC = NC+1; TO_CHECK(NC) = I-1
      END IF
      IF ((I .LT. ND) .AND. (.NOT. CHECKING(I))) THEN
         CHECKING(I) = .TRUE.; NC = NC+1; TO_CHECK(NC) = I
      END IF
   END DO shrink_values
   ! Identify which spline segments are nonmonotone after the updates.
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
! Use FIT_SPLINE to fit the final MQSI.
CALL FIT_SPLINE(X, FX, T, BCOEF, INFO)

CONTAINS

FUNCTION IS_MONOTONE(U0, U1, F0, F1, DF0, DF1, DDF0, DDF1)
! Given an interval [U0, U1] and function values F0, F1, first derivatives
! DF0, DF1, and second derivatives DDF0, DDF1 at U0, U1, respectively,
! return TRUE if the quintic polynomial matching these values is monotone
! over [U0, U1], and FALSE otherwise.
!
REAL(KIND=R8), INTENT(IN) :: U0, U1, F0, F1, DF0, DF1, DDF0, DDF1
LOGICAL :: IS_MONOTONE
! Local variables.
REAL(KIND=R8) :: A, ALPHA, B, BETA, GAMMA, INV_SLOPE, SIGN, TAU, TEMP, W
! Identify the direction of change of the function (increasing or decreasing).
IF      (F1 .GT. F0+ACCURACY) THEN; SIGN =  1.0_R8
ELSE IF (F1 .LT. F0-ACCURACY) THEN; SIGN = -1.0_R8
! When the function values are flat, everything *must* be 0.
ELSE
   IS_MONOTONE = (ABS(DF0)  .LT. ACCURACY) .AND. (ABS(DF1)  .LT. ACCURACY) .AND. &
                 (ABS(DDF0) .LT. ACCURACY) .AND. (ABS(DDF1) .LT. ACCURACY)
   RETURN
END IF
W = U1 - U0
INV_SLOPE = W / (F1 - F0)
! Determine which set of monotonicity conditions to use based on the
! assigned first derivative values at either end of the interval.
IF ((ABS(DF0) .LT. ACCURACY) .OR. (ABS(DF1) .LT. ACCURACY)) THEN
   ! Simplified monotone case, whichredcues to a test of cubic
   ! positivity studied in
   ! 
   ! J. W. Schmidt and W. He{\ss}, ``Positivity of cubic polynomials on
   ! intervals and positive spline interpolation'', {\sl BIT Numerical
   ! Mathematics}, 28 (1988) 340--352.
   ! 
   ! Notably, monotonicity results when the following conditions hold:
   !    alpha >= 0,
   !    delta >= 0,
   !    gamma >= delta - 2 * \sqrt{ alpha delta },
   !    beta  >= alpha - 2 * \sqrt{ alpha delta },
   ! 
   ! where alpha, beta, delta, and gamma are defined in the paper. The
   ! equations that follow are the result of algebreic simplifications
   ! of the terms as they are defined by Schmidt and He{\ss}.
   ! 
   ! The check for delta >= 0 has already been implicitly performed
   ! above, next check for alpha >= 0.
   IF (SIGN*DDF1 .GT. SIGN*4.0_R8 * DF1 / W) THEN; IS_MONOTONE = .FALSE.
   ELSE
      ! Now compute the value of 2 * \sqrt{ alpha delta }, store in TEMP.
      TEMP = SIGN * 2.0_R8 * INV_SLOPE * SQRT(DF0 * (4*DF1 - DDF1*W))
      ! Check for gamma >= delta - 2 * \sqrt{ alpha delta }
      IF (TEMP + INV_SLOPE*(3.0_R8*DF0 + DDF0*W) .LT. 0.0_R8) THEN
         IS_MONOTONE = .FALSE.
      ! Check for beta >= alpha - 2 * \sqrt{ alpha delta }
      ELSE IF (60.0_R8 + 2.0_R8*TEMP + INV_SLOPE*((5.0_R8*DDF1 - 3.0_R8*DDF0) &
         * W - 8.0_R8*(3.0_R8*DF0 + 4.0_R8*DF1)) .LT. 0.0_R8) THEN
         IS_MONOTONE = .FALSE.
      ELSE; IS_MONOTONE = .TRUE.
      END IF
   END IF
ELSE
   ! Full quintic monotoncity case related to the theory in
   !
   ! G. Ulrich and L. T. Watson, ``Positivity conditions for quartic
   ! polynomials'', {\sl SIAM J. Sci. Comput.}, 15 (1994) 528--544.
   ! 
   ! The following terms A, B, TAU, ALPHA, GAMMA, and BETA all
   ! directly correspond to notation defined by Ulrich and Watson.
   A = INV_SLOPE * DF0
   B = INV_SLOPE * DF1
   ! Full monotone quintic case.
   TAU = 24.0_R8 + 2.0_R8*SQRT(A*B) - 3.0_R8*(A+B)
   IF (TAU .LT. 0.0_R8) THEN; IS_MONOTONE = .FALSE.
   ELSE
      ! Compute alpha, gamma, beta from theorems to determine monotonicity.
      TEMP = (DF0 / DF1)**(0.75_R8)
      ALPHA = TEMP * (4.0_R8*DF1 - DDF1*W) / DF0
      GAMMA = (4.0_R8*DF0 + DDF0*W) / (TEMP * DF1)
      BETA = (3.0_R8 * INV_SLOPE * ((DDF1-DDF0)*W - 8.0_R8*(DF0+DF1)) + &
         60.0_R8) / (2.0 * SQRT(A*B))
      IF (BETA .LE. 6.0_R8) THEN; TEMP = -(BETA + 2.0_R8) / 2.0_R8
      ELSE;                       TEMP = -2.0_R8 * SQRT(BETA - 2.0_R8)
      END IF
      IS_MONOTONE = (ALPHA .GT. TEMP) .AND. (GAMMA .GT. TEMP)
   END IF
END IF
END FUNCTION IS_MONOTONE

SUBROUTINE QUADRATIC(I2)
! Given an index I2, compute the coefficients A of x^2 and B
! of x for the quadratic interpolating Y(I2-1:I2+1) at X(I2-1:I2+1),.
INTEGER, INTENT(IN) :: I2
! Local variables.
REAL(KIND=R8) :: C1, C2, C3, & ! Intermediate terms for computation.
  D ! Denominator for computing A and B via Cramer's rule.
INTEGER :: I1, I3
I1 = I2-1
I3 = I2+1
! Compute the shared denominator.
D = (X(I1) - X(I2)) * (X(I1) - X(I3)) * (X(I2) - X(I3))
! Compute coefficients A and B in quadratic interpolant  Ax^2 + Bx + C.
C1 = X(I1) * (Y(I3) - Y(I2))
C2 = X(I2) * (Y(I1) - Y(I3))
C3 = X(I3) * (Y(I2) - Y(I1))
A = (C1 + C2 + C3) / D
B = - (X(I1)*C1 + X(I2)*C2 + X(I3)*C3) / D
END SUBROUTINE QUADRATIC

END SUBROUTINE MQSI
