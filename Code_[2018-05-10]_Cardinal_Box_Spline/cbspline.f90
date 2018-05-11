! This file (cbspline.f90) contains the functions 
! 
!  CBSPLEV  -- Evaluate an arbitrary degree cardinal b-spline with
!              uniformly spaced knots on given range [lower, upper].
! 
!  BOXSPLEV -- Evaluate a box spline with a direction vector set
!              composed of <degree> repetitions of the identity
!              matrix and [lower,upper] corners of box.
!  
RECURSIVE FUNCTION CBSPLEV(X, DEGREE, LOWER, UPPER, IDX) RESULT(Y)
  ! Function for evaluating a cardinal b-spline of prescribed degree,
  ! which has evenly spaced knots located over [LOWER, UPPER] and
  ! nonzero support over the range (LOWER,UPPER).
  ! 
  ! Inputs:
  !   X      -- The (real-valued) location to evaluate a cardinal b-spline at.
  !   DEGREE -- The degree of the cardinal b-spline to evaluate.
  ! 
  ! Optional Input:
  !   LOWER  -- The location of the minimum-valued knot for this
  !             b-spline, default value of "0".
  !   UPPER  -- The location of the maximum-valued knot for this
  !             b-spline, default value of "DEGREE+1".
  !   IDX    -- A private input used for recurively evaluating the
  !             cardinal b-spline. *NOT FOR EXTERNAL USE*
  ! 
  ! Output:
  !   Y      -- The value of a cardinal b-spline of prescribed degree
  !             evaluated at the given position.
  ! 
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE

  ! Inputs
  REAL(KIND=REAL64), INTENT(IN)           :: X
  INTEGER, INTENT(IN)                     :: DEGREE
  REAL(KIND=REAL64), OPTIONAL, INTENT(IN) :: LOWER, UPPER  
  ! Optional Input
  INTEGER, INTENT(IN), OPTIONAL :: IDX
  ! Output
  REAL(KIND=REAL64) :: Y
  ! Local variables for optionals
  INTEGER :: INDEX
  REAL(KIND=REAL64) :: LOCAL_LOWER, LOCAL_UPPER, LOCAL_X

  LOCAL_X = X
  ! Initialize optional inputs
  DEFINE_LOWER : IF (PRESENT(LOWER)) THEN
     LOCAL_LOWER = LOWER
     ! Rescale X to be supported by the cardinal b-spline over [0, DEGREE+1]
     LOCAL_X = LOCAL_X - LOCAL_LOWER
  ELSE
     LOCAL_LOWER = 0
  END IF DEFINE_LOWER
  DEFINE_UPPER : IF (PRESENT(UPPER)) THEN
     LOCAL_UPPER = UPPER
     ! Rescale X to be supported by the cardinal b-spline over [0, DEGREE+1]
     LOCAL_X = (DEGREE+1) * LOCAL_X / (LOCAL_UPPER - LOCAL_LOWER)
  ELSE
     LOCAL_UPPER = DEGREE + 1
  END IF DEFINE_UPPER
  DEFINE_INDEX : IF (PRESENT(IDX)) THEN
     INDEX = IDX
  ELSE
     INDEX = 1
  END IF DEFINE_INDEX

  ! Body of cardinal b-spline evaluation function.
  IF (DEGREE .EQ. 0) THEN
     ! Base case, 1 for supported points.
     IF ((INDEX-1 .LT. LOCAL_X) .AND. (LOCAL_X .LE. INDEX)) THEN
        Y = 1
     ELSE
        Y = 0
     END IF
  ELSE
     ! Recursive evaluation of Cardinal B-Spline
     Y = (LOCAL_X - (INDEX-1)) / (DEGREE) * &
          CBSPLEV(LOCAL_X, DEGREE-1, IDX=INDEX) + &
          (INDEX + DEGREE - LOCAL_X) / (INDEX+DEGREE - INDEX) * &
          CBSPLEV(LOCAL_X, DEGREE-1, IDX=INDEX+1)
  END IF
END FUNCTION CBSPLEV

FUNCTION CBOXSPLEV(X, DEGREE, LOWER, UPPER) RESULT(Y)
  ! Function for evaluating a product of axis-aligned cardinal b-splines
  ! of arbitrary degree at a point X in [LOWER,UPPER]^SIZE(X).
  ! This will be referred to as a "cardinal box spline".
  ! 
  ! Inputs:
  !   X      -- Point in R^d (d-dimensional real vector space) to
  !             evaluate this box spline at.
  !   DEGREE -- The degree of the box spline being evaluated.
  ! 
  ! Optional Inputs:
  !   LOWER  -- Point in R^d defining the lower corner of this box 
  !             spline, default value of "0".
  !   UPPER  -- Point in R^d defining the upper corner of this box
  !             spline, default value of "DEGREE+1".
  ! 
  ! Output:
  !   Y      -- The value of this cardinal box spline at the given
  !             position for the given degree and bounds.
  ! 
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE

  ! Inputs
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: X
  INTEGER, INTENT(IN) :: DEGREE
  ! Optional Inputs
  REAL(KIND=REAL64), INTENT(IN), OPTIONAL, DIMENSION(:) :: LOWER, UPPER
  ! Output
  REAL(KIND=REAL64) :: Y
  ! Local variable (iterator)
  INTEGER :: I
  REAL(KIND=REAL64), DIMENSION(SIZE(X)) :: LOCAL_LOWER, LOCAL_UPPER
  ! Defined interface for the CBSPLEV function being called internally.
  INTERFACE 
     RECURSIVE FUNCTION CBSPLEV(X, DEGREE, LOWER, UPPER, IDX) RESULT(Y)
       USE ISO_FORTRAN_ENV, ONLY: REAL64
       REAL(KIND=REAL64), INTENT(IN)           :: X
       REAL(KIND=REAL64), OPTIONAL, INTENT(IN) :: LOWER, UPPER  
       INTEGER, INTENT(IN)                     :: DEGREE
       INTEGER, INTENT(IN), OPTIONAL :: IDX
       REAL(KIND=REAL64) :: Y
     END FUNCTION CBSPLEV
  END INTERFACE

  ! Defining optionals
  DEFINE_LOWER : IF (PRESENT(LOWER)) THEN
     LOCAL_LOWER = LOWER
  ELSE
     LOCAL_LOWER = 0
  END IF DEFINE_LOWER
  DEFINE_UPPER : IF (PRESENT(UPPER))  THEN
     LOCAL_UPPER = UPPER
  ELSE
     LOCAL_UPPER = DEGREE + 1
  END IF DEFINE_UPPER

  ! Body of function
  Y = 1
  EVALUATE_BOX_SPLINE : DO I = 1, SIZE(X)
     IF (Y .LT. EPSILON(Y)) EXIT EVALUATE_BOX_SPLINE
     Y = Y * CBSPLEV(X(I), DEGREE, LOCAL_LOWER(I), LOCAL_UPPER(I))
  END DO EVALUATE_BOX_SPLINE
END FUNCTION CBOXSPLEV
