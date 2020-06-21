! This file (sample_main.f90) contains a sample main program that
! calls MQSI to interpolate handcrafted nonmonotone data.

PROGRAM SAMPLE_MAIN
USE REAL_PRECISION, ONLY: R8  
IMPLICIT NONE

! Define the interfaces for relevant MQSI package subroutines.
INTERFACE
 SUBROUTINE MQSI(X, Y, T, BCOEF, INFO)
   USE REAL_PRECISION, ONLY: R8
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X, Y
   REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: T, BCOEF
   INTEGER, INTENT(OUT) :: INFO
 END SUBROUTINE MQSI
 SUBROUTINE EVAL_SPLINE(T, BCOEF, XY, INFO, D)
   USE REAL_PRECISION, ONLY: R8
   REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: T, BCOEF
   REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
   INTEGER, INTENT(OUT) :: INFO
   INTEGER, INTENT(IN), OPTIONAL :: D
 END SUBROUTINE EVAL_SPLINE
END INTERFACE

! Definitions of data.
INTEGER, PARAMETER :: ND = 18
REAL(KIND=R8), DIMENSION(ND) :: X
REAL(KIND=R8), DIMENSION(ND) :: &
     Y = (/ 0.0_R8, 1.0_R8, 1.0_R8, 1.0_R8, 0.0_R8, 20.0_R8, &
            19.0_R8, 18.0_R8, 17.0_R8, 0.0_R8, 0.0_R8, 3.0_R8, &
            0.0_R8, 1.0_R8, 6.0_R8, 16.0_R8, 16.1_R8, 1.0_R8 /)
! Definition of holder for spline knots, coefficients, and temporary data.
REAL(KIND=R8) :: SC(1:3*ND), SK(1:3*ND+6), U(ND)
! Loop index and info integer for checking exeuction status.
INTEGER :: I, INFO
! Initialize "X" values.
DO I = 1, ND ; X(I) = REAL(I, KIND=R8) ; END DO
! ------------------------------------------------------------------------------
!                   Rescaling operations (for testing)
! X(:) = X(:) * (2.0 * SQRT(EPSILON(1.0_R8))) ! Low separation X values.
! X(:) = X(:) * (10.0_R8**38 / (MAXVAL(X(:)) * 2.0_R8)) ! High separation X values.
! Y(:) = Y(:) * (10.0_R8**(-38)) ! Low separation Y values.
! Y(:) = Y(:) * (10.0_R8**20 / (MAXVAL(Y(:)) * 2.0_R8)) ! High separation Y values.
! ------------------------------------------------------------------------------
WRITE (*,*) ''
WRITE (*,"('     X: ',100ES10.2)") X(:)
WRITE (*,"('     Y: ',100ES10.2)") Y(:)
! Construct a MQSI.
CALL MQSI(X, Y, SK, SC, INFO)
IF (INFO .NE. 0) THEN; WRITE (*,"('MQSI error ',I3)") INFO; END IF
! Copy the "evaluation points" into the input/output buffer.
U(:) = X(:)
! Evaluate the spline at all points (result is updated in-place in U).
CALL EVAL_SPLINE(SK, SC, U, INFO)
IF (INFO .NE. 0) THEN; WRITE (*,"('EVAL_SPLINE error ',I3)") INFO; END IF
! Show the evaluations.
WRITE (*,"('  Q(X): ',100ES10.2)") U(:)
! Evaluate the first derivative of the spline at all points.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=1)
IF (INFO .NE. 0) THEN; WRITE (*,"('EVAL_SPLINE error ',I3)") INFO; END IF
WRITE (*,"(' DQ(X): ',100ES10.2)") U(:)
! Evaluate the second derivative of the spline at all points.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=2)
IF (INFO .NE. 0) THEN; WRITE (*,"('EVAL_SPLINE error ',I3)") INFO; END IF
WRITE (*,"('DDQ(X): ',100ES10.2)") U(:)
! Show the error of the values.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=0)
IF (INFO .NE. 0) THEN; WRITE (*,"('EVAL_SPLINE error ',I3)") INFO; END IF
WRITE (*,*) ''
WRITE (*,"('Machine precision:  ',ES18.10)") EPSILON(1.0_R8)
WRITE (*,"('   square root ^^:  ',ES18.10)") SQRT(EPSILON(1.0_R8))
WRITE (*,"('Max relative error: ',ES18.10)") MAXVAL(ABS(U(:) - Y(:)) / ABS(1.0_R8 + ABS(Y(:))))
WRITE (*,"('Max absolute error: ',ES18.10)") MAXVAL(ABS(U(:) - Y(:)))
WRITE (*,*) '' 

END PROGRAM SAMPLE_MAIN
  
