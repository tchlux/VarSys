! This file contains code relevant to testing a local installation of
! the MQSI package. Example command to compile and run all tests:
! 
!  $F08 $OPTS REAL_PRECISION.f90 EVAL_BSPLINE.f90 SPLINE.f90 MQSI.f90 \
!    test.f90 -o test lapack.f blas.f && ./test
! 
! where '$F08' is the name of the Fortran 2008 compiler, '$OPTS' are
! compiler options such as '-O3'.
! 
! All tests will be run and an error message will be printed for any
! failed test cases. Otherwise a message saying tests have passed will
! be directed to standard output. A small timing test will also be
! run to display the expected time required to fit and evaluate a
! MQSI to data of a given size.
! 
! 
! ====================================================================
!                 Test program for MQSI package.
! 
PROGRAM TEST_ALL
USE REAL_PRECISION, ONLY: R8
IMPLICIT NONE
! ------------------------------------------------------------------
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
 SUBROUTINE FIT_SPLINE(XI, FX, T, BCOEF, INFO)
   USE REAL_PRECISION, ONLY: R8
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:)   :: XI
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: FX
   REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: T, BCOEF
   INTEGER, INTENT(OUT) :: INFO
 END SUBROUTINE FIT_SPLINE
END INTERFACE
! ------------------------------------------------------------------
INTEGER :: TRIALS
REAL(KIND=R8) :: ERROR_TOLERANCE
INTEGER :: NS(6)
LOGICAL :: EQ_SPACED(2)
CHARACTER(LEN=20) :: TEST_FUNC_NAMES(5)
INTEGER :: TIME_SIZE
! Iteration indices.
INTEGER :: I, J, K, N
LOGICAL :: ALL_PASSED
! Time recording variables.
REAL :: EVAL_TIMES(5,5,2), FIT_TIMES(5,5,2)
! ------------------------------------------------------------------
!                         Testing parameters
! 
TRIALS = 100
ERROR_TOLERANCE = SQRT(EPSILON(1.0_R8))
NS(:) = (/ 2**3, 2**4, 2**5, 2**6, 2**7, 2**8 /)
EQ_SPACED(:) = (/ .TRUE., .FALSE. /)
TEST_FUNC_NAMES(:) = (/ "Large tangent       ", "Piecewise polynomial", &
   "Random              ", "Random monotone     ", "Signal decay        " /)
TIME_SIZE = 100
! ------------------------------------------------------------------
ALL_PASSED = .TRUE.
! Run all tests.
WRITE (*,*) ''
WRITE (*,"('-------------------------------')")
WRITE (*,"('Running tests, ',I3,' trials each.')") TRIALS
main_loop : DO I = 1, SIZE(NS)
  N = NS(I)
  WRITE (*,*) ''
  WRITE (*,"('N: ',I5)") N
  ! Run a test with very nearby points.
  WRITE (*,"('  Tiny values')")
  IF (.NOT. EPSILON_TEST(N)) THEN
        ! Describe the test that failed and exit.
        WRITE (*,"('__________________________________')")
        WRITE (*,*) ''
        ALL_PASSED = .FALSE.
        EXIT main_loop
  END IF

  ! Run a test with maximally spaced points.
  WRITE (*,"('  Huge values')")
  IF (.NOT. HUGE_TEST(N)) THEN
        ! Describe the test that failed and exit.
        WRITE (*,"('__________________________________')")
        WRITE (*,*) ''
        ALL_PASSED = .FALSE.
        EXIT main_loop
  END IF

  ! Run a test on each test function (with this number of points).
  DO J = 1, SIZE(TEST_FUNC_NAMES)
    WRITE (*,"('  ',A)") TEST_FUNC_NAMES(J)
    DO K = 1, SIZE(EQ_SPACED)
      IF (.NOT. PASSES_TEST()) THEN
        ! Describe the test that failed and exit.
        WRITE (*,*) ''
        WRITE (*,"('Test configuration:')")
        WRITE (*,"('  function:   ',A)") TEST_FUNC_NAMES(J)
        IF (EQ_SPACED(K)) THEN ; WRITE (*,"('  spacing:    equally spaced')")
        ELSE ;                   WRITE (*,"('  spacing:    randomly spaced')")
        END IF
        WRITE (*,"('  data size: ',1I5)") N
        WRITE (*,"('  num trials:',1I5)") TRIALS
        WRITE (*,"('__________________________________')")
        WRITE (*,*) ''
        ALL_PASSED = .FALSE.
        EXIT main_loop
      END IF
    END DO
  END DO
END DO main_loop
run_timing_test : IF (ALL_PASSED) THEN
WRITE (*,*) '' ; WRITE (*,"('All tests PASSED.')") ; 
! End of testing code, beginning of timing code.
WRITE (*,*) ''
WRITE (*,"('--------------------------------------------------------------')")
N = TIME_SIZE
WRITE (*,"('Computing timing data for size ',I5,'.')") N
DO J = 1, SIZE(TEST_FUNC_NAMES)
   DO K = 1, SIZE(EQ_SPACED)
      CALL RUN_TIME_TEST(EVAL_TIMES(:,J,K), FIT_TIMES(:,J,K))
   END DO
END DO
! Average the percentiles over all executed tests.
EVAL_TIMES(:,1,1) = SUM(SUM(EVAL_TIMES(:,:,:), DIM=3), DIM=2) / &
     REAL(SIZE(EVAL_TIMES,2) * SIZE(EVAL_TIMES,3))
FIT_TIMES(:,1,1) = SUM(SUM(FIT_TIMES(:,:,:), DIM=3), DIM=2) / &
     REAL(SIZE(FIT_TIMES,2) * SIZE(FIT_TIMES,3))
! Print out the results.
WRITE (*,*) ''
WRITE (*,"(' Fit time percentiles for MQSI of ',I5,' points:')") N
J = SIZE(EVAL_TIMES,1)
DO I = 1, J
  IF (I .EQ. 1) THEN
    WRITE (*,"('      min  ',F8.6,' seconds')") FIT_TIMES(I,1,1)
  ELSE IF ((I-1 .EQ. (J-1)/2) .AND. (MOD(I,2) .EQ. 1)) THEN
    WRITE (*,"('   median  ',F8.6,' seconds')") FIT_TIMES(I,1,1)
  ELSE IF (I .EQ. J) THEN
    WRITE (*,"('      max  ',F8.6,' seconds')") FIT_TIMES(I,1,1)
  ELSE
    K = INT(1.0 + REAL((99)*(I-1))/REAL(J-1))
    WRITE (*,"('      ',I3,'  ',F8.6,' seconds')") K, FIT_TIMES(I,1,1)
  END IF
END DO
WRITE (*,*) ''
WRITE (*,"(' Evaluation time per point for MQSI built from ',I5,' points:')") N
J = SIZE(EVAL_TIMES,1)
DO I = 1, J
  IF (I .EQ. 1) THEN
    WRITE (*,"('      min  ',F8.6,' seconds')") EVAL_TIMES(I,1,1)
  ELSE IF ((I-1 .EQ. (J-1)/2) .AND. (MOD(I,2) .EQ. 1)) THEN
    WRITE (*,"('   median  ',F8.6,' seconds')") EVAL_TIMES(I,1,1)
  ELSE IF (I .EQ. J) THEN
    WRITE (*,"('      max  ',F8.6,' seconds')") EVAL_TIMES(I,1,1)
  ELSE
    K = INT(1.0 + REAL((99)*(I-1))/REAL(J-1))
    WRITE (*,"('      ',I3,'  ',F8.6,' seconds')") K, EVAL_TIMES(I,1,1)
  END IF
END DO
END IF run_timing_test

CONTAINS

FUNCTION PASSES_TEST()
  ! Wrapper function for converting "J" into a specific test function
  ! and executing the subsequent test.
  LOGICAL :: PASSES_TEST
  ! Large tangent
  IF (J .EQ. 1) THEN
     CALL RUN_TEST(N, EQ_SPACED(K), LARGE_TANGENT, &
          TEST_FUNC_NAMES(J), PASSES_TEST)
  ! Piecewise polynomial
  ELSE IF (J .EQ. 2) THEN
     CALL RUN_TEST(N, EQ_SPACED(K), PIECEWISE_POLYNOMIAL, &
          TEST_FUNC_NAMES(J), PASSES_TEST)
  ! Random data
  ELSE IF (J .EQ. 3) THEN
     CALL RUN_TEST(N, EQ_SPACED(K), RANDOM, &
          TEST_FUNC_NAMES(J), PASSES_TEST)
  ! Random monotone data
  ELSE IF (J .EQ. 4) THEN
     CALL RUN_TEST(N, EQ_SPACED(K), RANDOM_MONOTONE, &
          TEST_FUNC_NAMES(J), PASSES_TEST)
  ! Decaying signal
  ELSE IF (J .EQ. 5) THEN
     CALL RUN_TEST(N, EQ_SPACED(K), SIGNAL_DECAY, &
          TEST_FUNC_NAMES(J), PASSES_TEST)
  ! Unknown test number
  ELSE
     WRITE (*,"('Unknown test number ',I2)") J
     PASSES_TEST = .FALSE.
  END IF
END FUNCTION PASSES_TEST

SUBROUTINE RUN_TIME_TEST(EVAL_TIME, FIT_TIME)
  ! Wrapper function for converting "J" into a specific test function
  ! and executing the subsequent test.
  REAL, DIMENSION(:), INTENT(OUT) :: EVAL_TIME, FIT_TIME
  ! Large tangent
  IF (J .EQ. 1) THEN
     CALL TIME_TEST(N, EQ_SPACED(K), LARGE_TANGENT, &
          TEST_FUNC_NAMES(J), FIT_TIME, EVAL_TIME)
  ! Piecewise polynomial
  ELSE IF (J .EQ. 2) THEN
     CALL TIME_TEST(N, EQ_SPACED(K), PIECEWISE_POLYNOMIAL, &
          TEST_FUNC_NAMES(J), FIT_TIME, EVAL_TIME)
  ! Random data
  ELSE IF (J .EQ. 3) THEN
     CALL TIME_TEST(N, EQ_SPACED(K), RANDOM, &
          TEST_FUNC_NAMES(J), FIT_TIME, EVAL_TIME)
  ! Random monotone data
  ELSE IF (J .EQ. 4) THEN
     CALL TIME_TEST(N, EQ_SPACED(K), RANDOM_MONOTONE, &
          TEST_FUNC_NAMES(J), FIT_TIME, EVAL_TIME)
  ! Decaying signal
  ELSE IF (J .EQ. 5) THEN
     CALL TIME_TEST(N, EQ_SPACED(K), SIGNAL_DECAY, &
          TEST_FUNC_NAMES(J), FIT_TIME, EVAL_TIME)
  ! Unknown test number
  ELSE
     WRITE (*,"('Unknown test number ',I2)") J
  END IF
END SUBROUTINE RUN_TIME_TEST


! ====================================================================
!                       Test execution routines.
! 
FUNCTION EPSILON_TEST(ND) RESULT(PASSES)
! Run a test on data that is at the limit of being as close together
! as possible.
INTEGER, INTENT(IN) :: ND
LOGICAL :: PASSES
! -------------------------------------------------
!               Local variables

REAL(KIND=R8) :: MAX_ERROR, SK(1:3*ND+6), SC(1:3*ND), &
     T(TRIALS), U(ND), X(ND), Y(ND), Z(TRIALS)
INTEGER :: I, INFO, J
INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
! -------------------------------------------------
DO I = 1, ND ; X(I) = (I-1) ; END DO
X(:) = X(:) / REAL(ND-1, KIND=R8)
CALL PIECEWISE_POLYNOMIAL(X, Y)
! Construct the smallest spacing of "X" values that is allowed by the
! MQSI routine (then we will check for correct monotonicity).
X(:) = (X(:) * (ND-1)) * SQRT(EPSILON(1.0_R8))
INFO = 5
grow_x_gap : DO WHILE ((INFO .EQ. 5) .OR. (INFO .EQ. 7))
   CALL MQSI(X, Y, SK, SC, INFO)
   IF ((INFO .NE. 5) .AND. (INFO .NE. 7)) EXIT grow_x_gap
   X(:) = X(:) * 1.1_R8 ! <- Growth factor that roughly matches runtime of other tests.
END DO grow_x_gap
IF (INFO .NE. 0) THEN
   WRITE (*,*) ''   
   WRITE (*,"('Failed to construct MQSI, code ',I3,'.')") INFO
   PASSES = .FALSE.
   RETURN
END IF
! Check that the spline reproduces the function values correctly.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=0)
IF (INFO .NE. 0) THEN
   WRITE (*,*) ''   
   WRITE (*,"('Failed to evaluate spline, code ',I3,'.')") INFO
   PASSES = .FALSE.
   RETURN
END IF
MAX_ERROR = MAXVAL( ABS((U(:) - Y(:))) / (1.0_R8 + ABS(Y(:))) )
IF (MAX_ERROR .GT. ERROR_TOLERANCE) THEN
   WRITE (*,*) ''
   WRITE (*,"('Value test:        FAILED')")
   WRITE (*,"('  relative error:  ',ES10.3)") MAX_ERROR
   WRITE (*,"('  error tolerance: ',ES10.3)") ERROR_TOLERANCE
   PASSES = .FALSE.
   RETURN
END IF
! Check for monotonicity over all intervals.
check_monotonicity :DO I = 1, ND-1
   DO J = 1, TRIALS ; Z(J) = J-1 ; END DO
   Z(:) = X(I) + (Z(:) / REAL(TRIALS-1, KIND=R8)) * (X(I+1) - X(I))
   CALL EVAL_SPLINE(SK, SC, Z, INFO, D=1)
   IF (.NOT. VALUES_ARE_MONOTONE(X(I), X(I+1), Y(I), Y(I+1), Z)) THEN
      PASSES = .FALSE.; RETURN
   END IF
END DO check_monotonicity
PASSES = .TRUE.
! End of test subroutine.
END FUNCTION EPSILON_TEST

!=====================================================================
FUNCTION HUGE_TEST(ND) RESULT(PASSES)
! Run a test on data that is at the limit of being as close together
! as possible.
INTEGER, INTENT(IN) :: ND
LOGICAL :: PASSES
! -------------------------------------------------
!               Local variables
REAL(KIND=R8) :: MAX_ERROR, SK(1:3*ND+6), SC(1:3*ND), &
     T(TRIALS), U(ND), X(ND), Y(ND), Z(TRIALS)
INTEGER :: I, INFO, J
INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
! -------------------------------------------------
DO I = 1, ND ; X(I) = (I-1) ; END DO
X(:) = X(:) / REAL(ND-1, KIND=R8)
CALL PIECEWISE_POLYNOMIAL(X, Y)
! Construct the largest spacing of X and Y values that is allowed
! by the MQSI routine (then we will check for correct monotonicity).
X(:) = (X(:) * (ND-1)) * (SQRT(SQRT(HUGE(1.0_R8))))
Y(:) = Y(:) * (SQRT(SQRT(HUGE(1.0_R8))))
INFO = 6
grow_x_gap : DO WHILE ((INFO .EQ. 6) .OR. (INFO .EQ. 7))
   CALL MQSI(X, Y, SK, SC, INFO)
   IF ((INFO .NE. 6) .AND. (INFO .NE. 7)) EXIT grow_x_gap
   X(:) = X(:) / 10.0_R8
   Y(:) = Y(:) / 10.0_R8
END DO grow_x_gap
IF (INFO .NE. 0) THEN
   WRITE (*,*) ''   
   WRITE (*,"('Failed to construct MQSI, code ',I3,'.')") INFO
   PASSES = .FALSE.
   RETURN
END IF
! Check that the spline reproduces the function values correctly.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=0)
IF (INFO .NE. 0) THEN
   WRITE (*,*) ''   
   WRITE (*,"('Failed to evaluate spline, code ',I3,'.')") INFO
   PASSES = .FALSE.
   RETURN
END IF
MAX_ERROR = MAXVAL( ABS((U(:) - Y(:))) / (1.0_R8 + ABS(Y(:))) )
IF (MAX_ERROR .GT. ERROR_TOLERANCE) THEN
   WRITE (*,*) ''
   WRITE (*,"('Value test:        FAILED')")
   WRITE (*,"('  relative error:  ',ES10.3)") MAX_ERROR
   WRITE (*,"('  error tolerance: ',ES10.3)") ERROR_TOLERANCE
   PASSES = .FALSE.
   RETURN
END IF
! Check for monotonicity over all intervals.
check_monotonicity :DO I = 1, ND-1
   DO J = 1, TRIALS ; Z(J) = J-1 ; END DO
   Z(:) = X(I) + (Z(:) / REAL(TRIALS-1, KIND=R8)) * (X(I+1) - X(I))
   CALL EVAL_SPLINE(SK, SC, Z, INFO, D=1)
   IF (.NOT. VALUES_ARE_MONOTONE(X(I), X(I+1), Y(I), Y(I+1), Z)) THEN
      PASSES = .FALSE.; RETURN
   END IF
END DO check_monotonicity
PASSES = .TRUE.
! End of test subroutine.
END FUNCTION HUGE_TEST

!=====================================================================
SUBROUTINE RUN_TEST(ND, EQUALLY_SPACED_X, F, F_NAME, PASSES)
! Runs a batch of tests and prints information to standard output
! on the status of the test (expected output, execution time, etc.)
INTEGER, INTENT(IN) :: ND
LOGICAL, INTENT(IN) :: EQUALLY_SPACED_X
CHARACTER(LEN=20), INTENT(IN) :: F_NAME
LOGICAL, INTENT(INOUT) :: PASSES
! -------------------------------------------------
!               Local variables
INTERFACE
 SUBROUTINE F(X, Y)
   USE REAL_PRECISION, ONLY : R8
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
   REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
 END SUBROUTINE F
END INTERFACE
REAL(KIND=R8) :: MAX_ERROR, SK(1:3*ND+6), SC(1:3*ND), &
     T(TRIALS), U(ND), X(ND), Y(ND), Z(TRIALS)
INTEGER :: I, INFO, J
INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
! -------------------------------------------------
! Initialize X values.
IF (EQUALLY_SPACED_X) THEN
   DO I = 1, ND ; X(I) = I-1; END DO
   X(:) = X(:) / REAL(ND-1, KIND=R8)
ELSE
   ! Initialize random seed.
   CALL RANDOM_SEED(SIZE=J)
   ALLOCATE(SEED(J))
   SEED(:) = 0
   CALL RANDOM_SEED(PUT=SEED)
   ! Generate random data.
   CALL RANDOM_NUMBER(X)
   CALL SORT(X)
   ! Make sure the random points have ample spacing.
   DO I = 1, ND-1
      X(I+1) = MAX(X(I+1), X(I)+ERROR_TOLERANCE*REAL(2**8,KIND=R8))
   END DO
END IF
! Initialize Y values.
CALL F(X, Y)
! Show the X and Y values.
CALL MQSI(X, Y, SK, SC, INFO)
IF (INFO .NE. 0) THEN
   WRITE (*,*) ''   
   WRITE (*,"('Failed to construct MQSI, code ',I3,'.')") INFO
   PASSES = .FALSE.
   RETURN
END IF
! Check that the spline reproduces the function values correctly.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=0)
IF (INFO .NE. 0) THEN
   WRITE (*,*) ''   
   WRITE (*,"('Failed to evaluate spline, code ',I3,'.')") INFO
   PASSES = .FALSE.
   RETURN
END IF
MAX_ERROR = MAXVAL( ABS((U(:) - Y(:))) / (1.0_R8 + ABS(Y(:))) )
IF (MAX_ERROR .GT. ERROR_TOLERANCE) THEN
   WRITE (*,*) ''
   WRITE (*,"('Value test:        FAILED')")
   WRITE (*,"('  relative error:  ',ES10.3)") MAX_ERROR
   WRITE (*,"('  error tolerance: ',ES10.3)") ERROR_TOLERANCE
   PASSES = .FALSE.
   RETURN
END IF
! Check for monotonicity over all intervals.
check_monotonicity :DO I = 1, ND-1
   DO J = 1, TRIALS ; Z(J) = J-1 ; END DO
   Z(:) = X(I) + (Z(:) / REAL(TRIALS-1, KIND=R8)) * (X(I+1) - X(I))
   CALL EVAL_SPLINE(SK, SC, Z, INFO, D=1)
   IF (.NOT. VALUES_ARE_MONOTONE(X(I), X(I+1), Y(I), Y(I+1), Z)) THEN
      PASSES = .FALSE.; RETURN
   END IF
END DO check_monotonicity
PASSES = .TRUE.
! End of test subroutine.
END SUBROUTINE RUN_TEST

!=====================================================================
SUBROUTINE TIME_TEST(ND, EQUALLY_SPACED_X, F, F_NAME, TFIT, TEVAL)
! Runs a batch of tests and records equally spaced percentiles of
! times into "TFIT" (fit time) and "TEVAL" (evaluation time).
INTEGER, INTENT(IN) :: ND
LOGICAL, INTENT(IN) :: EQUALLY_SPACED_X
CHARACTER(LEN=20), INTENT(IN) :: F_NAME
REAL, INTENT(OUT), DIMENSION(:) :: TFIT, TEVAL
! -------------------------------------------------
!               Local variables
INTERFACE
 SUBROUTINE F(X, Y)
   USE REAL_PRECISION, ONLY : R8
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
   REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
 END SUBROUTINE F
END INTERFACE
REAL(KIND=R8) :: MAX_ERROR, SK(1:3*ND+6), SC(1:3*ND), &
     T(TRIALS), U(ND), X(ND), Y(ND), Z(TRIALS)
REAL :: FINISH_TIME_SEC, START_TIME_SEC
INTEGER :: I, INFO, J
INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
! -------------------------------------------------
! Initialize random seed.
CALL RANDOM_SEED(SIZE=J)
ALLOCATE(SEED(J))
SEED(:) = 0
CALL RANDOM_SEED(PUT=SEED)
! Initialize X values.
IF (EQUALLY_SPACED_X) THEN
   DO I = 1, ND ; X(I) = I-1; END DO
   X(:) = X(:) / REAL(ND-1, KIND=R8)
ELSE
   CALL RANDOM_NUMBER(X)
   CALL SORT(X)
   ! Make sure the random points have ample spacing.
   DO I = 1, ND-1
      X(I+1) = MAX(X(I+1), X(I)+ERROR_TOLERANCE*REAL(2**8,KIND=R8))
   END DO
END IF
! Initialize Y values.
CALL F(X, Y)
DO I = 1, TRIALS
   ! Construct the piecewise monotone quintic spline interpolant.
   CALL CPU_TIME(START_TIME_SEC)
   CALL MQSI(X, Y, SK, SC, INFO)
   CALL CPU_TIME(FINISH_TIME_SEC)
   T(I) = REAL(FINISH_TIME_SEC - START_TIME_SEC, KIND=R8)
   IF (INFO .NE. 0) THEN
      WRITE (*,*) ''   
      WRITE (*,"('Failed to construct MQSI, code ',I3,'.')") INFO
      RETURN
   END IF
END DO
CALL SORT(T)
J = SIZE(TFIT)
DO I = 0, J-1
   TFIT(I+1) = T( INT(1.0 + REAL((TRIALS-1)*I)/REAL(J-1)) )
END DO
! Print out some evaluations of the MQSI.
DO I = 1, TRIALS
   CALL RANDOM_NUMBER(Z)
   CALL CPU_TIME(START_TIME_SEC)
   CALL EVAL_SPLINE(SK, SC, Z, INFO, D=0)
   IF (INFO .NE. 0) THEN
      WRITE (*,*) ''   
      WRITE (*,"('Failed to evaluate spline, code ',I3,'.')") INFO
      RETURN
   END IF
   CALL CPU_TIME(FINISH_TIME_SEC)
   T(I) = REAL(FINISH_TIME_SEC - START_TIME_SEC, KIND=R8)
END DO
CALL SORT(T)
! Rescale the timings to be relative to a single approximation point.
T(:) = T(:) / REAL(ND, KIND=R8)
! Grab percentiles.
J = SIZE(TEVAL)
DO I = 0, J-1
   TEVAL(I+1) = T( INT(1.0 + REAL((TRIALS-1)*I)/REAL(J-1)) )
END DO
! End of test subroutine.
END SUBROUTINE TIME_TEST


! ====================================================================
!                       Testing functions.
! 
SUBROUTINE LARGE_TANGENT(X, Y)
  !  "large tangent" function that rapidly grows from near 0 to 99 on right.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  Y(:) = -(1.0_R8 + ( 1.0_R8 / (X(:) - 1.01_R8) ))
END SUBROUTINE LARGE_TANGENT

SUBROUTINE PIECEWISE_POLYNOMIAL(X, Y)
  !  "piecewise polynomial" that is formed from a C^1 spline.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  ! ------------------------------------
  ! Define the piecewise polynomial.
  INTEGER, PARAMETER :: ND = 18
  REAL(KIND=R8), DIMENSION(1:ND) :: XI
  REAL(KIND=R8), DIMENSION(1:ND,1:2) :: FX
  REAL(KIND=R8), DIMENSION(1:ND*2+2*2) :: SK
  REAL(KIND=R8), DIMENSION(1:ND*2) :: SC
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
  DO I = 1, ND ; XI(I) = REAL(I-1, KIND=R8) ; END DO
  ! Rescale XI and FX to be on the unit interval.
  XI(:) = XI(:) / REAL(ND-1, KIND=R8)
  FX(:,2) = FX(:,2) * REAL(ND-1, KIND=R8)
  CALL FIT_SPLINE(XI, FX, SK, SC, INFO)
  ! Make sure all provided values are inside the support of the spline.
  Y(:) = MIN(MAX(X(:), 0.0_R8), 1.0_R8)
  CALL EVAL_SPLINE(SK, SC, Y, INFO, D=0)
END SUBROUTINE PIECEWISE_POLYNOMIAL

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

SUBROUTINE SIGNAL_DECAY(X, Y)
  !  "signal" function that is a decaying magnitude sine wave.
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: Y
  REAL(KIND=R8) :: PI
  PI = ACOS(-1.0_R8)
  Y(:) = SIN(8.0_R8 * PI * X(:)) / (X(:)**2 + 0.1_R8)
END SUBROUTINE SIGNAL_DECAY


! ====================================================================
!                       Utility functions.
! 
FUNCTION VALUES_ARE_MONOTONE(U0, U1, F0, F1, DQ)
REAL(KIND=R8), INTENT(IN) :: U0, U1, F0, F1
REAL(KIND=R8), DIMENSION(:) :: DQ
LOGICAL :: VALUES_ARE_MONOTONE
REAL(KIND=R8) :: TM76
TM76 = 10.0_R8**(-76)
! Rescale the DQ (derivative) values by the function change on the interval.
DQ(:) = DQ(:) / (1.0_R8 + ABS(F1 - F0))
! The computation of derivatives necessitates a higher tolerance for error.
IF ( ((F1 .GT. F0+TM76) .AND. (ANY(DQ(:) .LT. -SQRT(ERROR_TOLERANCE)))) .OR. &
     ((F1 .LT. F0-TM76) .AND. (ANY(DQ(:) .GT.  SQRT(ERROR_TOLERANCE)))) .OR. &
     ((ABS(F1-F0) .LT. TM76) .AND. (MAXVAL(ABS(DQ(:))) .GT. SQRT(ERROR_TOLERANCE))) ) THEN
   WRITE (*,*) ''
   WRITE (*,"('Monotonicity test: FAILED')")
   WRITE (*,"('  interval [',ES10.3,', ',ES10.3,']')") U0, U1
   WRITE (*,"('  interval function change: ',ES10.3)") F1 - F0
   WRITE (*,"('  minimum derivative value: ',ES10.3)") MINVAL(DQ(:))
   WRITE (*,"('  maximum derivative value: ',ES10.3)") MAXVAL(DQ(:))
   WRITE (*,"('  error tolerance:          ',ES10.3)") ERROR_TOLERANCE
   VALUES_ARE_MONOTONE = .FALSE.
ELSE ; VALUES_ARE_MONOTONE = .TRUE.
END IF
END FUNCTION VALUES_ARE_MONOTONE

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

END PROGRAM TEST_ALL

