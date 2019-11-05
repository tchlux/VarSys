MODULE SPLINES
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    IMPLICIT NONE

CONTAINS
  
  SUBROUTINE FIT_SPLINE(KNOTS, VALUES, SPLINE_KNOTS, SPLINE_COEFFICIENTS)
    ! Subroutine for fitting a spline composed of B-splines that
    ! interpolate the given function value (and derivatives) at knots.
    ! 
    !   -> KNOTS( knot values )
    !   -> VALUES( knot index, func val )
    ! 
    !      this function
    ! 
    !   <-  SPLINE_KNOTS
    !   <-  SPLINE_COEFFICIENTS
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
    PRINT '()'
    PRINT '("K(1:",I2,") = ",f6.2)', DEGREE+1, KNOTS(1)
    SPLINE_KNOTS(:DEGREE+1) = KNOTS(1)
    DO STEP = 2, SIZE(KNOTS)-1
       PRINT '("K(",I2,":",I2,") = ",F6.2)', &
            DEGREE+(STEP-2)*CONT + 2, DEGREE+(STEP-1)*CONT + 1, KNOTS(STEP)
       SPLINE_KNOTS(DEGREE+(STEP-2)*CONT + 2 : &
                    DEGREE+(STEP-1)*CONT + 1) = KNOTS(STEP)
    END DO
    SPLINE_KNOTS(SIZE(SPLINE_KNOTS)-DEGREE:) = KNOTS(SIZE(KNOTS))
    PRINT '("K(",I2,":",I2,") = ",F6.2)', &
         SIZE(SPLINE_KNOTS)-DEGREE, SIZE(SPLINE_KNOTS), KNOTS(SIZE(KNOTS))
    PRINT '()'

    ! Copy the VALUES into the SPLINE_COEFFICIENTS (output) variable.
    DO STEP = 1, SIZE(VALUES,1)
       SPLINE_COEFFICIENTS(1+(STEP-1)*CONT : STEP*CONT) = VALUES(STEP,:)
    END DO

    ! Evaluate each B-spline and it's appropriate derivatives at all
    ! knots. Each B-spline will have nonzero values at knots for at
    ! most (2*SIZE(VALUES,2) - 1) knots. That suggests a column of at
    ! most (2*SIZE(VALUES,2) - 1) * SIZE(VALUES,2) total values.
    ! 
    ! (K x 3 x N) (N x 1) = (K x 3 x 1) ! Evaluation at all points.
    ! (3K x N) (N x 1)    = (3K x 1)    ! N = num splines, K = num knots
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
    ! Evaluate a spline construced with FIT_SPLINE.
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: COEFFICIENTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: X
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(X)) :: Y
    REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: D
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
    !     It should be noted that in the case of repeated knots, a minimum
    !     allowed denominator value is substituted to prevent a NaN result.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: X
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(X)) :: Y
    INTEGER, INTENT(IN), OPTIONAL :: D
    ! Set the minimum distance between knots before they're considered repeated.
    REAL(KIND=REAL64), PARAMETER :: MIN_DIV = SQRT(EPSILON(1.0_REAL64))
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
          IF (DIV_LEFT .GE. MIN_DIV) THEN
             IF (DIV_RIGHT .GE. MIN_DIV) THEN
                VALUES(:,K) =  (STEP-1) * (&
                     VALUES(:,K) / DIV_LEFT - VALUES(:,K+1) / DIV_RIGHT )
             ELSE
                VALUES(:,K) =  (STEP-1) * (VALUES(:,K) / DIV_LEFT)
             END IF
          ELSE
             IF (DIV_RIGHT .GE. MIN_DIV) THEN
                VALUES(:,K) =  (STEP-1) * ( - VALUES(:,K+1) / DIV_RIGHT )
             END IF
          END IF
       END DO
    END DO compute_derivative

    ! Assign the values into the "Y" output.
    Y(:) = VALUES(:,1)

  END SUBROUTINE EVAL_BSPLINE

END MODULE SPLINES
