! This file (sample_main.f90) contains a sample main program that calls
! PMQSI to interpolate handcrafted nonmonotone data.

PROGRAM SAMPLE_MAIN
USE REAL_PRECISION, ONLY: R8  
IMPLICIT NONE

! Define the interfaces for relevant PMQSI package subroutines.
INTERFACE
 SUBROUTINE PMQSI(X, Y, T, BCOEF, INFO)
   USE REAL_PRECISION, ONLY: R8
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X, Y
   REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: T, BCOEF
   INTEGER, INTENT(OUT) :: INFO
 END SUBROUTINE PMQSI
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
REAL(KIND=R8), PARAMETER, DIMENSION(ND) :: &
     Y = (/ 0.0_R8, 1.0_R8, 1.0_R8, 1.0_R8, 0.0_R8, 20.0_R8, &
            19.0_R8, 18.0_R8, 17.0_R8, 0.0_R8, 0.0_R8, 3.0_R8, &
            0.0_R8, 1.0_R8, 6.0_R8, 16.0_R8, 16.1_R8, 1.0_R8 /)
! Definition of holder for spline knots, coefficients, and temporary data.
REAL(KIND=R8) :: SC(1:3*ND), SK(1:3*ND+6), U(ND)
! Loop index and info integer for checking exeuction status.
INTEGER :: I, INFO

! Initialize "X" values.
DO I = 1, ND ; X(I) = REAL(I, KIND=R8) ; END DO
! Construct a PMQSI.
CALL PMQSI(X, Y, SK, SC, INFO)
IF (INFO .NE. 0) PRINT "('PMQSI returned info 'I3)", INFO
! Copy the "evaluation points" into the input/output buffer.
U(:) = X(:)
! Evaluate the spline at all points (result is updated in-place in U).
CALL EVAL_SPLINE(SK, SC, U, INFO)
IF (INFO .NE. 0) PRINT "('EVAL_SPLINE returned info 'I3)", INFO
! Show the evaluations.
PRINT "('')"
PRINT "('    X:  '100F7.2)", X(:)
PRINT "('    Y:  '100F7.2)", Y(:)
PRINT "('  F(X): '100F7.2)", U(:)
PRINT "('')"
PRINT "('Machine precision:  'F18.16)", EPSILON(1.0_R8)
PRINT "('Max relative error: 'F18.16)", MAXVAL(ABS(U(:) - Y(:)) / ABS(1.0_R8 + ABS(Y(:))))
PRINT "('')"
! Evaluate the first derivative of the spline at all points.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=1)
IF (INFO .NE. 0) PRINT "('EVAL_SPLINE returned info 'I3)", INFO
PRINT "(' DF(X): '100F7.2)", U(:)
! Evaluate the second derivative of the spline at all points.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=2)
IF (INFO .NE. 0) PRINT "('EVAL_SPLINE returned info 'I3)", INFO
PRINT "('DDF(X): '100F7.2)", U(:)

END PROGRAM SAMPLE_MAIN
  
