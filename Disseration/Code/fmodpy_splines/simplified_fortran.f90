
MODULE SPLINES
USE ISO_FORTRAN_ENV , ONLY : REAL64


FUNCTION L2 ( KNOTS , VALUES1 , VALUES2 , DIVISIONS )
! Given a set of knots and two sets of coefficients, where each
! set defines a spline, compute the L2 difference between the two
! splines using Simpson's rule.
!
!   KNOTS(K) -- The non-decreasing sequence of unique real-valued
!               break points, with at least 2 unique values.
!   VALUES1(K,C1) -- The function (and derivative) values of a polynomial at knots.
!   VALUES2(K,C2) -- The function (and derivative) values of a polynomial at knots.
!   DIVISIONS -- The number of divisions used to estimate the
!                integral over EACH interval of the splines.
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : KNOTS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : VALUES1 , VALUES2
INTEGER , INTENT ( IN ) , OPTIONAL : : DIVISIONS
! Local variables.
! First define variables for the local spline fits of the values.
REAL ( KIND = REAL64 ) : : KNOTS1 ( SIZE ( VALUES1 ) + 2 * SIZE ( VALUES1 , 2 ) ) , KNOTS2 ( SIZE ( VALUES2 ) + 2 * SIZE ( VALUES2 , 2 ) ) , COEFS1 ( SIZE ( VALUES1 ) ) , COEFS2 ( SIZE ( VALUES2 ) )
! Now declare extra auxillary variables needed for computing L2.
INTEGER : : D , STEP
REAL ( KIND = REAL64 ) : : SQ_SUM , L2
REAL ( KIND = REAL64 ) , ALLOCATABLE , DIMENSION ( : ) : : REL_X , X , EVALS1 , EVALS2
! Declare the (even) number of divisions.
! Allocate space for storing the spline evaluations on each interval.
! Initialize a relative set of X to be evenly spaced on [0,1].
! Fit the splines to the two sets of knots and values.
! Cycle over all the intervals tracking the squared sum.
! Evaluate the (squared) difference between the splines at points.
! Compute the squared L2 over this interval with Simpson's rule.
! Compute the final L2 as the square root of the sum of the squares.
END FUNCTION L2


SUBROUTINE MULTIPLY_SPLINES ( KNOTS1 , COEFS1 , KNOTS2 , COEFS2 , SPLINE_KNOTS , SPLINE_COEFFICIENTS )
! Given two splines defined by their knot and coefficient
! sequences, compute the spline that is the pointwise
! multiplication of the two splines.
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : KNOTS1 , KNOTS2 , COEFS1 , COEFS2
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( KNOTS1 ) + SIZE ( KNOTS2 ) ) : : SPLINE_KNOTS
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( ( SIZE ( KNOTS1 ) + SIZE ( KNOTS2 ) ) * 2 * MAX ( SIZE ( KNOTS1 ) - SIZE ( COEFS1 ) , SIZE ( KNOTS2 ) - SIZE ( COEFS2 ) ) ) : : SPLINE_COEFFICIENTS
! Determine the knots that will be kept.
! Evaluate at appropriate points.
! Solve the linear system determining the spline coefficients.
END SUBROUTINE MULTIPLY_SPLINES

SUBROUTINE ADD_SPLINES ( KNOTS1 , COEFS1 , KNOTS2 , COEFS2 , SPLINE_KNOTS , SPLINE_COEFFICIENTS , NUM_KNOTS , NUM_COEFS )
! Given two splines defined by their knot and coefficient
! sequences, compute the spline that is the pointwise
! addition of the two splines.
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : KNOTS1 , KNOTS2 , COEFS1 , COEFS2
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( KNOTS1 ) + SIZE ( KNOTS2 ) ) : : SPLINE_KNOTS
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( COEFS1 ) + SIZE ( COEFS2 ) ) : : SPLINE_COEFFICIENTS
INTEGER , INTENT ( OUT ) : : NUM_KNOTS , NUM_COEFS
! Local variables.
INTEGER : : I , I1 , I2 , K , ORDER
! Determine the knots that will be kept.
!
! knots1:  5 [0, 1, 1.5, 2, 3]
! knots2:  6 [1.25, 1.5, 1.75, 2, 2.5, 5]
! answer:  9 [0, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 5]
!
! K1 size:  7  K2 size: 16
!
! knots1:  5 [0, 1, 1.5, 2, 3]
! knots2:  6 [1.25, 1.5, 1.75, 2, 2.5, 5]
! answer:  9 [0, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 5]
!
! K1 size:  7  K2 size: 16
! K1:   0.00  0.00  1.00  1.50  2.00  3.00  3.00
! K2:   1.25  1.25  1.25  1.25  1.50  1.50  1.75  1.75  2.00  2.00  2.50  2.50  5.00  5.00  5.00  5.00
!
! 1 I1:   1   I2:   1  0.00  1.25 -1.00
! 2 I1:   4   I2:   1  1.50  1.25 -1.00
! 1 I1:   4   I2:   7  1.50  1.75 -1.00
! 2 I1:   5   I2:   7  2.00  1.75 -1.00
! 1 I1:   5   I2:  11  2.00  2.50 -1.00
! 2 I1:   6   I2:  11  3.00  2.50 -1.00
! 1 I1:   6   I2:  13  3.00  5.00 -1.00
! 2 I1:   7   I2:  13  3.00  5.00 -1.00
! 1 I1:   7   I2:  13  3.00  5.00 -1.00
! 2 I1:   7   I2:  13  3.00  5.00 -1.00
!
! Cycle I1 until it is greater or done
! Count the number of equal valued knots.
! Skip over shared knots.
! Cycle I2 until it is greater or done
! Count the number of equal valued knots.
! If the copying is done, then exit the loop.
! Skip over shared knots.
! Get the size of the knots, set the last knot to be repeated filling the spline.
! Identify the maximum order of the two different splines.
! Allocate space for all of the function evaluations that will
! determine this spline fit over the data.
!
!   Number of evaluations = (ORDER - 2) * (NUM_KNOTS - 1)
!

! ! Solve the linear system determining the spline coefficients.
! !
! ! ----------------------------------------------------------------
! ! Initialize all untouched values in AB to zero.
! AB(:,:) = 0_REAL64

! FIRST = 1
! LAST  = 1
! ! Evaluate the B-splines at all of the knots.
! DO STEP = 1, SIZE(SPLINE_COEFFICIENTS)
!    ! Increment the "FIRST" knot index until it is contained.
!    DO WHILE ((FIRST .LT. SIZE(KNOTS)) .AND. &
!         (SPLINE_KNOTS(MIN(SIZE(SPLINE_KNOTS),STEP+1)) .GT. KNOTS(FIRST)))
!       FIRST = FIRST + 1
!    END DO
!    ! Increment the "LAST" knot index until it is contained.
!    DO WHILE ((LAST .LT. SIZE(KNOTS)) .AND. (SPLINE_KNOTS( &
!         MIN(SIZE(SPLINE_KNOTS),STEP+DEGREE)) .GT. KNOTS(LAST)))
!       LAST = LAST + 1
!    END DO
!    ! Find the range of indices that will be written to in AB.
!    ! The mapping is looks like   AB[LK+KU+1+i-j,j] = A[i,j]
!    START = 2*DEGREE+1 + (1 + (FIRST-1)*CONT - STEP)
!    END = MIN(SIZE(AB,1), START + (LAST-FIRST+1)*CONT - 1)
!    ! Evaluate each derivative of this B-spline at relevant knots.
!    DO DERIV = 0, CONT-1
!       ! Place the evaluations into a block out of a column in AB,
!       ! shift according to which derivative is being evaluated
!       ! and use a stride appropriate for the continuity.
!       AB(START+DERIV:END:CONT,STEP) = KNOTS(FIRST:LAST)
!       CALL EVAL_BSPLINE(SPLINE_KNOTS(STEP:STEP+DEGREE+1), &
!            AB(START+DERIV:END:CONT,STEP), DERIV)
!    END DO
! ! ----------------------------------------------------------------


END SUBROUTINE ADD_SPLINES



SUBROUTINE FIT_SPLINE ( KNOTS , VALUES , SPLINE_KNOTS , SPLINE_COEFFICIENTS )
! Subroutine for fitting a spline composed of B-splines that
! interpolates the given function value (and function derivatives)
! at provided knots.
!
! INPUT:
!   KNOTS(N)    -- The non-decreasing real-valued locations of the
!                  break-points bfor the interpolating spline.
!   VALUES(N,C) -- The real-valued function values (column 1), and
!                  derivative values (columns 2 onwards) that the
!                  interpolating spline should reproduce.
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
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : KNOTS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : VALUES
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( VALUES ) + 2 * SIZE ( VALUES , 2 ) ) : : SPLINE_KNOTS
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( VALUES ) ) : : SPLINE_COEFFICIENTS
! Use the LAPACK function for solving banded linear systems.
EXTERNAL : : DGBSV
INTEGER , DIMENSION ( SIZE ( SPLINE_COEFFICIENTS ) ) : : IPIV
REAL ( KIND = REAL64 ) , DIMENSION ( 1 + 3 * ( 2 * SIZE ( VALUES , 2 ) - 1 ) , SIZE ( SPLINE_COEFFICIENTS ) ) : : AB
! Local variables.
INTEGER : : STEP , DERIV , DEGREE , CONT , INFO
INTEGER : : START , END , FIRST , LAST

! Define some local variables for notional convenience.

! Copy over the knots that will define the B-spline representation.

! Copy the VALUES into the SPLINE_COEFFICIENTS (output) variable.

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

! Evaluate the B-splines at all of the knots.
! Increment the "FIRST" knot index until it is contained.
! Increment the "LAST" knot index until it is contained.
! Find the range of indices that will be written to in AB.
! The mapping is looks like   AB[LK+KU+1+i-j,j] = A[i,j]
END = MIN ( SIZE ( AB , 1 ) , START + ( LAST - FIRST + 1 ) * CONT - 1 )
! Evaluate each derivative of this B-spline at relevant knots.
! Place the evaluations into a block out of a column in AB,
! shift according to which derivative is being evaluated
! and use a stride appropriate for the continuity.

! Call the function to solve the system.
! Check for errors.
END SUBROUTINE FIT_SPLINE


SUBROUTINE EVAL_SPLINE ( KNOTS , COEFFICIENTS , XY , D , EXTRAP )
! Evaluate a spline construced with FIT_SPLINE. Similar interface
! to EVAL_BSPLINE. Evaluate D derivative at all XY, result in XY.
!
! TODO: Change extrapolation to keep fit values and derivatives
!       a constant, use Newton form extrapolation from value
!       and derivatives (0's for integration). Give extrapolation
!       parameter that makes extrapolation extend E derivatives.
!
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : KNOTS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : COEFFICIENTS
REAL ( KIND = REAL64 ) , INTENT ( INOUT ) , DIMENSION ( : ) : : XY
INTEGER , INTENT ( IN ) , OPTIONAL : : D , EXTRAP
! Local variables.
INTEGER : : DERIV , STEP , NUM_KNOTS , E , E_STEP
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( XY ) , SIZE ( COEFFICIENTS ) ) : : VALUES
! Variables for extrapolation-related computations.
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( XY ) ) : : CONTAINED_XY
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( KNOTS ) - SIZE ( COEFFICIENTS ) - 2 , 2 , SIZE ( COEFFICIENTS ) ) : : EDGE_VALUES
! Compute the NUM_KNOTS (number of knots) for each B-spline.
! Assign the local value of the optional derivative "D" argument.
! Assign the local value of the optional extrapolation argument.
! This value is the number of derivatives maintained constant.
! Cap the values to either end of the spline at the first and last knots.
! Evaluate all splines at all the X positions.
! ! Compute the value of the spline at the ends (for extrapolation).
! IF (E .GT. 0) THEN
!    EDGE_VALUES(:,1,:) = KNOTS(1)
!    EDGE_VALUES(:,2,:) = KNOTS(SIZE(KNOTS))
!    DO E_STEP = 1, MIN(E, SIZE(KNOTS) - SIZE(COEFFICIENTS) - 2)
!       DO STEP = 1, SIZE(COEFFICIENTS)
!          CALL EVAL_BSPLINE(KNOTS(STEP:STEP+NUM_KNOTS), &
!               EDGE_VALUES(:,:,STEP), E_STEP)
!       END DO
!    END DO
!    EDGE_VALUES(:,:,1) = MATMUL(EDGE_VALUES, COEFFICIENTS)
!    ! For each extrapolation point, convert it's position
! END IF
! Store the values into Y as the weighted sums of B-spline evaluations.
END SUBROUTINE EVAL_SPLINE


SUBROUTINE EVAL_BSPLINE ( KNOTS , XY , D )
! Subroutine for evaluating a B-spline with provided knot sequence.
!
! INPUT:
!   KNOTS -- The nondecreasing sequence of break points for the B-spline.
!
! INPUT / OUTPUT:
!   XY    -- The locations at which the B-spline is evaluated on
!            input, on output holds the value of the B-spline with
!            prescribed knots evaluated at the given X locations.
!
! OPTIONAL INPUT:
!   D [= 0]  --  The derivative to take of the evaluated B-spline.
!                When negative, this subroutine integrates the B-spline.
!
!
! DESCRIPTION:
!
!    This function uses the recurrence relation defining a B-spline:
!
!      B_{K,1}(X)   =   1     if KNOTS(K) <= X < KNOTS(K+1),
!                       0     otherwise,
!
!    where K is the knot index, I = 2, ..., SIZE(KNOTS)-MAX(D,0)-1, and
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
!    For the computation of the integral of the B-spline, the
!    continuation of the above formula proceeds one step at a
!    time by adding a duplicate of the last knot, raising the
!    order of all intermediate B-splines, summing their values,
!    and dividing the sums by the width of the supported interval.
!
!    For the computation of the derivative of the B-spline, the
!    continuation of the standard recurrence relation is used that
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
!     The B-spline is right continuous everywhere except at the
!     last knot, at which it is both left and right continuous.
!
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : KNOTS
REAL ( KIND = REAL64 ) , INTENT ( INOUT ) , DIMENSION ( : ) : : XY
INTEGER , INTENT ( IN ) , OPTIONAL : : D
! Local variables.
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( XY ) , SIZE ( KNOTS ) ) : : VALUES
INTEGER : : K , STEP , DERIV
REAL ( KIND = REAL64 ) : : DIV_LEFT , DIV_RIGHT
! Assign the local value of the optional derivative "D" argument.

! For integration, cap all X values to be the last knot.

! Initialize all values to 0.

! Assign the first value for each knot index.
! Make the last knot inclusive on the right.
! Make all knots before the last right continuous.
END WHERE
! Make the last interval left AND right continuous.
END WHERE

! Compute the value of the B-spline from the knot sequence.
! Cycle over each knot accumulating the values for the recurrence
! computation of the B-spline (skip last steps for derivatives).
! Assure that the divisor will not cause invalid computations.
! Compute the B-spline recurrence relation.

! Compute the derivative of the B-spline (if D > 0).
! Cycle over each knot, following the same structure with the
! derivative computing relation instead of the B-spline one.
! Assure that the divisor will not cause invalid computations.
! Compute the derivative recurrence relation.

! If integration is being performed, then we need to raise the
! order of all (sub) B-splines to be the same order as the first.
! Zero out repetitions of the last knot to prevent an integral
! that is too large at exactly the right endpoint.

! Loop through starting at the back, raising the order by
! normal evaluation of a B-spline.

! Compute the integral of the B-spline (if D < 0).
! Finish the integral by summing the constituent basis
! functions, then rescaling them according to their domain.

! Assign the values into the "Y" output.
END SUBROUTINE EVAL_BSPLINE

END MODULE SPLINES

