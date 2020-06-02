! Test file for Piecewise Monotone Quintic Spline Interpolant (PMQSI) construction.

MODULE TESTING
  USE REAL_PRECISION, ONLY : R8
CONTAINS

SUBROUTINE SORT(VALUES)
  ! Insertion sort an array of values.
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: VALUES
  ! Local variables.
  REAL(KIND=R8) :: TEMP_VAL
  INTEGER :: I, J, K
  ! Return for the base case.
  IF (SIZE(VALUES) .LE. 1) RETURN
  ! Put the smallest value at the front of the list.
  I = MINLOC(VALUES,1)
  TEMP_VAL = VALUES(1)
  VALUES(1) = VALUES(I)
  VALUES(I) = TEMP_VAL
  ! Insertion sort the rest of the array.
  DO I = 3, SIZE(VALUES)
     TEMP_VAL = VALUES(I)
     ! Search backwards in the list until the insertion location is found. 
     J = I - 1
     K  = I
     DO WHILE (TEMP_VAL .LT. VALUES(J))
        VALUES(K) = VALUES(J)
        J = J - 1
        K  = K - 1
     END DO
     ! Put the value into its place (where it is greater than the
     ! element before it, but less than all values after it).
     VALUES(K) = TEMP_VAL
  END DO
END SUBROUTINE SORT

SUBROUTINE PIECEWISE_POLYNOMIAL(X, Y)
  !  "piecewise polynomial" that is formed from a C^1 spline.
  USE SPLINES, ONLY : FIT_SPLINE, EVAL_SPLINE
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  ! ------------------------------------
  ! Define the piecewise polynomial.
  INTEGER, PARAMETER :: Z = 18
  REAL(KIND=R8), DIMENSION(1:Z) :: XI
  REAL(KIND=R8), DIMENSION(1:Z,1:2) :: FX
  REAL(KIND=R8), DIMENSION(1:Z*2+2*2) :: SK
  REAL(KIND=R8), DIMENSION(1:Z*2) :: SC
  INTEGER :: I, INFO
  !   xi               f(x)      Df(x)  
  FX(1, 1:2)  = (/  0.0_R8,  1.0_R8 /)
  FX(2, 1:2)  = (/  1.0_R8,  0.0_R8 /)
  FX(3, 1:2)  = (/  1.0_R8,  0.0_R8 /)
  FX(4, 1:2)  = (/  1.0_R8,  0.0_R8 /)
  FX(5, 1:2)  = (/  0.0_R8,  0.0_R8 /)
  FX(6, 1:2)  = (/ 20.0_R8, -1.0_R8 /)
  FX(7, 1:2)  = (/ 19.0_R8, -1.0_R8 /)
  FX(8, 1:2)  = (/ 18.0_R8, -1.0_R8 /)
  FX(9, 1:2)  = (/ 17.0_R8, -1.0_R8 /)
  FX(10, 1:2) = (/  0.0_R8,  0.0_R8 /)  
  FX(11, 1:2) = (/  0.0_R8,  0.0_R8 /)  
  FX(12, 1:2) = (/  3.0_R8,  0.0_R8 /)
  FX(13, 1:2) = (/  0.0_R8,  0.0_R8 /)  
  FX(14, 1:2) = (/  1.0_R8,  3.0_R8 /)
  FX(15, 1:2) = (/  6.0_R8,  9.0_R8 /)
  FX(16, 1:2) = (/ 16.0_R8,  0.1_R8 /)
  FX(17, 1:2) = (/ 16.1_R8,  0.1_R8 /)
  FX(18, 1:2) = (/ 1.0_R8, -15.0_R8 /)
  DO I = 1, Z ; XI(I) = REAL(I-1, KIND=R8) ; END DO
  ! Rescale XI and FX to be on the unit interval.
  XI(:) = XI(:) / REAL(Z-1, KIND=R8)
  FX(:,2) = FX(:,2) * REAL(Z-1, KIND=R8)
  CALL FIT_SPLINE(XI, FX, SK, SC, INFO)
  ! ------------------------------------
  ! Make sure all provided values are inside the support of the spline.
  Y(:) = MIN(MAX(X(:), 0.0_R8), 1.0_R8)
  CALL EVAL_SPLINE(SK, SC, Y, INFO)
END SUBROUTINE PIECEWISE_POLYNOMIAL

SUBROUTINE SIGNAL(X, Y)
  !  "signal" function that is a decaying magnitude sine wave.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  REAL(KIND=R8) :: PI
  PI = ACOS(-1.0_R8)
  Y(:) = SIN(8.0_R8 * PI * X(:)) / (X(:)**2 + 0.1_R8)
END SUBROUTINE SIGNAL

SUBROUTINE LARGE_TANGENT(X, Y)
  !  "large tangent" function that rapidly grows from near 0 to 99 on right.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  Y(:) = -(1.0_R8 + ( 1.0_R8 / (X(:) - 1.01_R8) ))
END SUBROUTINE LARGE_TANGENT

SUBROUTINE RANDOM(X, Y)
  !  "random" function that is just that, purely random data.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  CALL RANDOM_NUMBER(Y)
END SUBROUTINE RANDOM

SUBROUTINE RANDOM_MONOTONE(X, Y)
  !  "random monotone" function that generates random monotone data.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  CALL RANDOM_NUMBER(Y)
  CALL SORT(Y)
END SUBROUTINE RANDOM_MONOTONE

END MODULE TESTING


! ====================================================================
PROGRAM TEST_SPLINES
  USE SPLINES, ONLY : PMQSI, EVAL_SPLINE
  USE REAL_PRECISION, ONLY : R8
  USE TESTING
  ! -------------------------------------------------
  !              Testing parameters
  INTEGER, PARAMETER :: Z = 18
  LOGICAL, PARAMETER :: EQUALLY_SPACED_X = .TRUE.
  PROCEDURE(TEST_FUNC), POINTER :: F => PIECEWISE_POLYNOMIAL
  ! -------------------------------------------------
  REAL(KIND=R8) :: U(Z), X(Z), Y(Z), SK(1:3*Z+6), SC(1:3*Z)
  INTEGER :: INFO, I
  ABSTRACT INTERFACE 
     SUBROUTINE TEST_FUNC(X, Y)
       USE REAL_PRECISION, ONLY : R8
       REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
       REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
     END SUBROUTINE TEST_FUNC
  END INTERFACE
  ! -------------------------------------------------

  ! Initialize X values.
  IF (EQUALLY_SPACED_X) THEN
     DO I = 1, Z ; X(I) = REAL(I-1, KIND=R8) / REAL(Z-1, KIND=R8) ; END DO
  ELSE
     CALL RANDOM_NUMBER(X)
     CALL SORT(X)
  END IF

  ! Initialize Y values.
  CALL F(X, Y)
 
  ! Show the X and Y values.
  PRINT "('')"
  PRINT "(' X: ',100ES12.3)", X
  PRINT "(' Y: ',100ES12.3)", Y

  ! Construct the piecewise monotone quintic spline interpolant.
  CALL PMQSI(X, Y, SK, SC, INFO)

  ! Print out some evaluations of the PMQSI.
  U(:) = X(:)
  CALL EVAL_SPLINE(SK, SC, U, INFO)
  PRINT "(' F: '100ES12.3)", U(:)
  U(:) = X(:)

  PRINT "('')"
  CALL EVAL_SPLINE(SK, SC, U, INFO, D=1)
  PRINT "(' DF:  '100ES12.3)", U(:)
  U(:) = X(:)
  CALL EVAL_SPLINE(SK, SC, U, INFO, D=2)
  PRINT "(' D2F: '100ES12.3)", U(:)
  U(:) = X(:)
  CALL EVAL_SPLINE(SK, SC, U, INFO, D=6)
  PRINT "(' D6F: '100ES12.3)", U(:)

  PRINT "('')"
  PRINT "(' MAX(ABS(F(X) - Y)) ')"
  U(:) = X(:)
  CALL EVAL_SPLINE(SK, SC, U, INFO)
  PRINT "(100ES12.3)", MAXVAL(ABS(U(:) - Y(:)))

END PROGRAM TEST_SPLINES

