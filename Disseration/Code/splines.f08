MODULE SPLINES
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    IMPLICIT NONE

CONTAINS
  
  FUNCTION IS_MONOTONE(KNOTS, COEFFICIENTS)
    ! Use a computational method to check for the monotonicity of a
    ! spline over its sequence of knots.
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: COEFFICIENTS
    LOGICAL :: IS_MONOTONE
    ! Local variables.
    REAL(KIND=REAL64), DIMENSION(1000) :: EVAL_POINTS, DERIV_VALS
    ! Evaluate the first derivative of the spline at all points.
    CALL EVAL_SPLINE(KNOTS, COEFFICIENTS, EVAL_POINTS, DERIV_VALS, 1)
    ! See if the function is monotone.
    IS_MONOTONE = (MAXVAL(DERIV_VALS) .GT. 0_REAL64) .AND. (MINVAL(DERIV_VALS) .LT. 0_REAL64)
  END FUNCTION IS_MONOTONE

  SUBROUTINE FIT_MONOTONE_SPLINE(KNOTS, VALUES, SPLINE_KNOTS, &
       SPLINE_COEFFICIENTS, ACCURACY)
    ! Given a sequence of KNOTS and VALUES, find the monotone spline
    ! that fits some values as close as possible to VALUES along the
    ! line between provided VALUES(:,2:) and VALUES(:,2:) = 0.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: VALUES
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(VALUES) + & 
         2*SIZE(VALUES,2)) :: SPLINE_KNOTS
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(VALUES)) :: SPLINE_COEFFICIENTS
    REAL(KIND=REAL64), INTENT(IN),  OPTIONAl :: ACCURACY
    ! Local variables.
    REAL(KIND=REAL64), DIMENSION(SIZE(VALUES,1), SIZE(VALUES,2)) :: TEST_VALUES
    REAL(KIND=REAL64) :: LOW, UPP, MID, ACC
    ! Check for improper usage. Compute the difference between neighbors.
    TEST_VALUES(2:,1) = VALUES(2:,1) - VALUES(:SIZE(VALUES,1)-1,1)
    IF ((MAXVAL(TEST_VALUES(2:,1)) .GT. 0) .AND. &
         (MINVAL(TEST_VALUES(2:,1)) .LT. 0)) THEN
       PRINT '("ERROR: The given VALUES were not monotone.")'
       RETURN
    END IF
    ! First test the end point.
    CALL FIT_SPLINE(KNOTS, VALUES, SPLINE_KNOTS, SPLINE_COEFFICIENTS)
    ! Check to see if this is monotone.
    IF (IS_MONOTONE(SPLINE_KNOTS, SPLINE_COEFFICIENTS)) RETURN
    ! Initialize the "ACCURACY" local variable.
    IF (PRESENT(ACCURACY)) THEN ; ACC = ACCURACY
    ELSE ; ACC = SQRT(EPSILON(ACC))
    END IF
    ! Initialize the endpoints for a binary search.
    LOW = 0_REAL64
    UPP = 1_REAL64
    ! Perform the binary search.
    DO WHILE ((UPP - LOW) .GT. ACC)
       ! Assign all the of the derivative values to be the midpoint.
       MID = (LOW + UPP) / 2_REAL64
       TEST_VALUES(:,2:) = VALUES(:,2:) * MID
       ! Fit a spline to the data.
       CALL FIT_SPLINE(KNOTS, TEST_VALUES, SPLINE_KNOTS, SPLINE_COEFFICIENTS)
       ! Set the new bound based on the monotonicity of the spline.
       IF (IS_MONOTONE(SPLINE_KNOTS, SPLINE_COEFFICIENTS)) THEN ; LOW = MID
       ELSE ; UPP = MID
       END IF
    END DO
    ! If the final values tested were not monotone, then redo
    ! the fit for the "LOW" value.
    IF (.NOT. IS_MONOTONE(SPLINE_KNOTS, SPLINE_COEFFICIENTS)) THEN
       TEST_VALUES(:,2:) = VALUES(:,2:) * LOW
       CALL FIT_SPLINE(KNOTS, TEST_VALUES, SPLINE_KNOTS, SPLINE_COEFFICIENTS)
    END IF
  END SUBROUTINE FIT_MONOTONE_SPLINE


  SUBROUTINE FIT_SPLINE(KNOTS, VALUES, SPLINE_KNOTS, SPLINE_COEFFICIENTS)
    ! Subroutine for fitting a spline composed of B-splines that
    ! interpolates the given function value (and function derivatives)
    ! at provided knots.
    ! 
    ! INPUT:
    !   KNOTS(N)    -- The non-decreasing real-valued locations of the
    !                  break-points bfor the interpolating spline.
    !   VALUES(N,C) -- The real-valued function values (column 1), and
    !                  derivative values (columns 2 onwards) that the
    !                  interpolating spline should produce.
    ! 
    ! OUTPUT:
    !   SPLINE_KNOTS(K) -- The non-decreasing real-valued locations
    !                      of the breakpoints for the B-splines that
    !                      compose the resulting spline.
    !   SPLINE_COEFFICIENTS(NC) -- The coefficients for each of the 
    !                      B-splines that compose the interpolating
    !                      spline.
    ! 
    !   
    ! DESCRIPTION:
    ! 
    !   This function uses the EVAL_BSPLINE method to evaluate the
    !   B-splines at all locations and the LAPACK routine
    ! 
    !     DGBSV
    ! 
    !   to compute the coefficients of all component B-splines. The
    !   resulting fit over general knots in the range [0,1] with
    !   values in [0,1] has roughly EPSILON(1.0_REAL64) error when
    !   evaluating the function values, and SQRT(EPSILON(1.0_REAL64))
    !   error when evaluating the 3rd derivative of a spline.
    ! 
    !   This code becomes numerically unstable when evaluating a
    !   derivative past the 3rd for a C^{>=10} spline. Use of an exact
    !   arithmetic implementation is recommended for any functions
    !   beyond C5.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: VALUES
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(VALUES) + & 
         2*SIZE(VALUES,2)) :: SPLINE_KNOTS
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(VALUES)) :: SPLINE_COEFFICIENTS
    ! Use the LAPACK function for solving banded linear systems.
    EXTERNAL :: DGBSV
    INTEGER, DIMENSION(SIZE(SPLINE_COEFFICIENTS)) :: IPIV
    REAL(KIND=REAL64), DIMENSION(1 + 3*(2*SIZE(VALUES,2)-1), SIZE(SPLINE_COEFFICIENTS)) :: AB
    ! Local variables.
    INTEGER :: STEP, DERIV, DEGREE, CONT, INFO
    INTEGER :: START, END, FIRST, LAST

    ! Define some local variables for notional convenience.
    CONT = SIZE(VALUES,2)
    DEGREE = 2*CONT - 1

    ! Copy over the knots that will define the B-spline representation.
    SPLINE_KNOTS(:DEGREE+1) = KNOTS(1)
    DO STEP = 2, SIZE(KNOTS)-1
       SPLINE_KNOTS(DEGREE+(STEP-2)*CONT + 2 : &
                    DEGREE+(STEP-1)*CONT + 1) = KNOTS(STEP)
    END DO
    SPLINE_KNOTS(SIZE(SPLINE_KNOTS)-DEGREE:) = KNOTS(SIZE(KNOTS))

    ! Copy the VALUES into the SPLINE_COEFFICIENTS (output) variable.
    DO STEP = 1, SIZE(VALUES,1)
       SPLINE_COEFFICIENTS(1+(STEP-1)*CONT : STEP*CONT) = VALUES(STEP,:)
    END DO

    ! Evaluate each B-spline and it's appropriate derivatives at all
    ! knots. Each B-spline will have nonzero values at knots for at
    ! most (2*SIZE(VALUES,2) - 1) knots. That suggests a column of at
    ! most (2*SIZE(VALUES,2) - 1) * SIZE(VALUES,2) total values.
    ! 
    ! For a C1 interpolating spline with 4 knots (each providing
    ! value and first derivative), matrix A will look like:
    ! 
    !   a10 b10 c10 d10 
    !   a11 b11 c11 d11 
    !       b20 c20 d20 e20 f20 
    !       b21 c21 d21 e21 f21 
    !           c30 d30 e30 f30 g30
    !           c31 d31 e31 f31 g31
    !                   e40 f40 g40 h40
    !                   e41 f41 g41 h41
    !               
    ! This comes from:
    ! 
    ! (2 x K x N) (N x 1) = (2 x K x 1) ! Evaluation at all points.
    ! (2K x N) (N x 1)    = (2K x 1)    ! N = num splines, K = num knots
    ! 
    ! Notice this matrix is banded with equal lower / upper bandwidths
    ! being the maximum number of knots for which a spline takes on a
    ! nonzero value. In general KL = KU = 2*SIZE(VALUES,2) - 1.

    ! Initialize all untouched values in AB to zero.
    AB(:,:) = 0_REAL64

    FIRST = 1
    LAST  = 1
    ! Evaluate the B-splines at all of the knots.
    DO STEP = 1, SIZE(SPLINE_COEFFICIENTS)
       ! Increment the "FIRST" knot index until it is contained.
       DO WHILE ((FIRST .LT. SIZE(KNOTS)) .AND. &
            (SPLINE_KNOTS(MIN(SIZE(SPLINE_KNOTS),STEP+1)) .GT. KNOTS(FIRST)))
          FIRST = FIRST + 1
       END DO
       ! Increment the "LAST" knot index until it is contained.
       DO WHILE ((LAST .LT. SIZE(KNOTS)) .AND. (SPLINE_KNOTS( &
            MIN(SIZE(SPLINE_KNOTS),STEP+DEGREE)) .GT. KNOTS(LAST)))
          LAST = LAST + 1
       END DO
       ! Find the range of indices that will be written to in AB.
       ! The mapping is looks like   AB[LK+KU+1+i-j,j] = A[i,j]
       START = 2*DEGREE+1 + (1 + (FIRST-1)*CONT - STEP)
       END = MIN(SIZE(AB,1), START + (LAST-FIRST+1)*CONT - 1)
       ! Evaluate each derivative of this B-spline at relevant knots.
       DO DERIV = 0, CONT-1
          ! Place the evaluations into a block out of a column in AB,
          ! shift according to which derivative is being evaluated
          ! and use a stride appropriate for the continuity.
          CALL EVAL_BSPLINE(SPLINE_KNOTS(STEP:STEP+DEGREE+1), &
               KNOTS(FIRST:LAST), AB(START+DERIV:END:CONT,STEP), DERIV)
       END DO
    END DO

    ! Call the function to solve the system.
    CALL DGBSV(SIZE(SPLINE_COEFFICIENTS), DEGREE, DEGREE, 1, AB, &
         SIZE(AB,1), IPIV, SPLINE_COEFFICIENTS, SIZE(VALUES), INFO)
    ! Check for errors.
    IF (INFO .NE. 0) &
         PRINT '("WARNING: DGBSV INFO flag",I3," on output.")', INFO
  END SUBROUTINE FIT_SPLINE


  SUBROUTINE EVAL_SPLINE(KNOTS, COEFFICIENTS, X, Y, D)
    ! Evaluate a spline construced with FIT_SPLINE. Similar interface
    ! to EVAL_BSPLINE. Evaluate D derivative at all X, store in Y.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: COEFFICIENTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: X
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(X)) :: Y
    INTEGER, INTENT(IN), OPTIONAL :: D
    ! Local variables.
    INTEGER :: DERIV, STEP, NUM_KNOTS
    REAL(KIND=REAL64), DIMENSION(SIZE(Y), SIZE(COEFFICIENTS)) :: VALUES
    ! Compute the NUM_KNOTS (number of knots) for each B-spline.
    NUM_KNOTS = SIZE(KNOTS) - SIZE(COEFFICIENTS)
    ! Assign the local value of the optional derivative "D" argument.
    set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
    ELSE ; DERIV = 0
    END IF set_derivative
    ! Evaluate all splines at all the X positions.
    DO STEP = 1, SIZE(COEFFICIENTS)
       CALL EVAL_BSPLINE(KNOTS(STEP:STEP+NUM_KNOTS), X, VALUES(:,STEP), DERIV)
    END DO
    ! Store the values into Y as the weighted sums of B-spline evaluations.
    Y(:) = MATMUL(VALUES, COEFFICIENTS)
  END SUBROUTINE EVAL_SPLINE


  SUBROUTINE EVAL_BSPLINE(KNOTS, X, Y, D)
    ! Subroutine for evaluating a B-spline with provided knot sequence.
    ! 
    ! INPUT:
    !   KNOTS -- The nondecreasing sequence of break points for the B-spline.
    !   X     -- The locations at which the B-spline is evaluated.
    ! 
    ! OUTPUT:
    !   Y  --  The value of the B-spline with prescribed knots evaluated
    !          at the given X locations.
    ! 
    ! OPTIONAL INPUT:
    !   D [= 0]  --  The derivative to take of the evaluated B-spline. 
    ! 
    ! 
    ! DESCRIPTION:
    ! 
    !    This function uses the recurrence relation defining a B-spline:
    ! 
    !      B_{K,1}(X)   =   1     if KNOTS(K) <= X < KNOTS(K+1),
    !                       0     otherwise,
    ! 
    !    where K is the knot index, I = 2, ..., SIZE(KNOTS)-D-1, and
    ! 
    !                               X - KNOTS(K)                      
    !      B_{K,I}(X) =      ------------------------- B_{K,I-1}(X)   
    !                         KNOTS(K+I-1) - KNOTS(K)                 
    !                                                                   
    !                             KNOTS(K+I) - X                    
    !                     +  ------------------------- B_{K+1,I-1}(X).
    !                         KNOTS(K+I) - KNOTS(K+1)               
    ! 
    !    However, all of the intermediate steps (I) are stored in a
    !    single block of memory that is overwritten at each step.
    ! 
    !    For the computation of the derivative of the B-spline, the
    !    continuation of the above recurrence relation is used that
    !    builds from I = SIZE(KNOTS)-D, ..., SIZE(KNOTS)-1 as
    ! 
    !                           (I-1) B_{K,I-1}(X)     
    !      B_{K,I}(X) =      -------------------------
    !                         KNOTS(K+I-1) - KNOTS(K)  
    !                                                   
    !                           (I-1) B_{K+1,I-1}(X)    
    !                     -  -------------------------.
    !                         KNOTS(K+I) - KNOTS(K+1) 
    ! 
    !     The B-spline is left continuous everywhere except at the
    !     last knot, at which it is both left and right continuous.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: X
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(X)) :: Y
    INTEGER, INTENT(IN), OPTIONAL :: D
    ! Local variables.
    REAL(KIND=REAL64), DIMENSION(SIZE(X), SIZE(KNOTS)) :: VALUES
    INTEGER :: K, STEP, DERIV
    REAL(KIND=REAL64) :: DIV_LEFT, DIV_RIGHT

    ! Assign the local value of the optional derivative "D" argument.
    set_derivative : IF (PRESENT(D)) THEN ; DERIV = D
    ELSE ; DERIV = 0
    END IF set_derivative

    ! Initialize all values to 0.
    VALUES(:,:) = 0.0_REAL64

    ! Assign the first value for each knot index.
    ! Make the last knot inclusive on the right.
    DO K = 1, SIZE(KNOTS)-1
       IF (KNOTS(K+1) .NE. KNOTS(SIZE(KNOTS))) THEN
          WHERE ( (KNOTS(K) .LE. X(:)) .AND. (X(:) .LT. KNOTS(K+1)) )
             VALUES(:,K) = 1.0_REAL64
          END WHERE
       ELSE
          WHERE ( (KNOTS(K) .LE. X(:)) .AND. (X(:) .LE. KNOTS(K+1)) )
             VALUES(:,K) = 1.0_REAL64
          END WHERE
       END IF
    END DO

    ! Compute the value of the B-spline from the knot sequence.
    compute_spline : DO STEP=2, SIZE(KNOTS)-DERIV-1
       ! Cycle over each knot accumulating the values for the recurrence
       ! computation of the B-spline (skip last steps for derivatives).
       DO K = 1, SIZE(KNOTS)-STEP
          ! Assure that the divisor will not cause invalid computations.
          DIV_LEFT = (KNOTS(K+STEP-1) - KNOTS(K))
          DIV_RIGHT = (KNOTS(K+STEP) - KNOTS(K+1))
          ! Compute the B-spline recurrence relation.
          IF (DIV_LEFT .GT. 0) THEN
             IF (DIV_RIGHT .GT. 0) THEN
                VALUES(:,K) = &
                     ((X(:) - KNOTS(K)) / DIV_LEFT) * VALUES(:,K) + &
                     ((KNOTS(K+STEP) - X(:)) / DIV_RIGHT) * VALUES(:,K+1)
             ELSE
                VALUES(:,K) = ((X(:) - KNOTS(K)) / DIV_LEFT) * VALUES(:,K)
             END IF
          ELSE
             IF (DIV_RIGHT .GT. 0) THEN
                VALUES(:,K) = ((KNOTS(K+STEP) - X(:)) / DIV_RIGHT) * VALUES(:,K+1)
             END IF
          END IF
       END DO
    END DO compute_spline

    ! Compute the derivative of the B-spline (if D > 0).
    compute_derivative : DO STEP = SIZE(KNOTS)-DERIV, SIZE(KNOTS)-1
       ! Cycle over each knot, following the same structure with the
       ! derivative computing relation instead of the B-spline one.
       DO K = 1, SIZE(KNOTS) - STEP
          ! Assure that the divisor will not cause invalid computations.
          DIV_LEFT = (KNOTS(K+STEP-1) - KNOTS(K))
          DIV_RIGHT = (KNOTS(K+STEP) - KNOTS(K+1))
          ! Compute the derivative recurrence relation.
          IF (DIV_LEFT .GT. 0) THEN
             IF (DIV_RIGHT .GT. 0) THEN
                VALUES(:,K) =  (STEP-1) * (&
                     VALUES(:,K) / DIV_LEFT - VALUES(:,K+1) / DIV_RIGHT )
             ELSE
                VALUES(:,K) =  (STEP-1) * (VALUES(:,K) / DIV_LEFT)
             END IF
          ELSE
             IF (DIV_RIGHT .GT. 0) THEN
                VALUES(:,K) =  (STEP-1) * ( - VALUES(:,K+1) / DIV_RIGHT )
             END IF
          END IF
       END DO
    END DO compute_derivative

    ! Assign the values into the "Y" output.
    Y(:) = VALUES(:,1)

  END SUBROUTINE EVAL_BSPLINE

END MODULE SPLINES
