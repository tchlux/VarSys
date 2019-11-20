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
    ! Use the monotonicity of the coefficients to determine spline monotonicity.
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: KNOTS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:) :: COEFFICIENTS
    LOGICAL :: IS_MONOTONE
    ! Local variables.
    INTEGER :: I, SIGN
    ! Check for monotonicity.
    SIGN = 0
    IS_MONOTONE = .TRUE.
    monotone_check : DO I = 1, SIZE(COEFFICIENTS)-1
       ! Assign a "SIGN" if it has not yet been assigned.
       IF (SIGN .EQ. 0) THEN
          IF (COEFFICIENTS(I+1) .LT. COEFFICIENTS(I)) THEN
             SIGN = -1
          ELSE IF (COEFFICIENTS(I+1) .GT. COEFFICIENTS(I)) THEN
             SIGN = 1
          END IF
       ! Check that the function is monotone increasing.
       ELSE IF (SIGN .EQ. 1) THEN
          IF (COEFFICIENTS(I+1) .LT. COEFFICIENTS(I)) THEN
             IS_MONOTONE = .FALSE.
             EXIT monotone_check
          END IF
       ! Check that the function is monotone decreasing.
       ELSE IF (SIGN .EQ. -1) THEN
          IF (COEFFICIENTS(I+1) .GT. COEFFICIENTS(I)) THEN
             IS_MONOTONE = .FALSE.
             EXIT monotone_check
          END IF
       END IF
    END DO
  END FUNCTION IS_MONOTONE

END MODULE MONOTONE_SPLINES
