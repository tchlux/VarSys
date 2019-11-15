MODULE MONOTONE_SPLINES
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  USE SPLINES
  IMPLICIT NONE

CONTAINS

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

  FUNCTION IS_MONOTONE(KNOTS, COEFFICIENTS)
    ! Use a computational method to check for the monotonicity of a
    ! spline over its sequence of knots.
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: COEFFICIENTS
    LOGICAL :: IS_MONOTONE
    ! Local variables.
    REAL(KIND=REAL64), DIMENSION(1000) :: DERIV_VALS
    INTEGER :: I
    ! Assign the locations to evaluate the spline.
    DO I = 1, SIZE(DERIV_VALS)
       DERIV_VALS(I) = I
    END DO
    DERIV_VALS(:) = (DERIV_VALS(:) / SIZE(DERIV_VALS)) * (KNOTS(SIZE(KNOTS))-KNOTS(1))
    ! Evaluate the first derivative of the spline at all points.
    CALL EVAL_SPLINE(KNOTS, COEFFICIENTS, DERIV_VALS, 1)
    ! See if the function is monotone.
    IS_MONOTONE = (MAXVAL(DERIV_VALS) .GT. 0_REAL64) .AND. (MINVAL(DERIV_VALS) .LT. 0_REAL64)
  END FUNCTION IS_MONOTONE

END MODULE MONOTONE_SPLINES
