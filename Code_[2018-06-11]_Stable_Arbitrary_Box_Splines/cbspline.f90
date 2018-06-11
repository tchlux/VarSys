! This file (cbspline.f90) contains the functions 
! 
!  CBSPLEV  -- Evaluate an arbitrary degree cardinal B-spline with
!              uniformly spaced knots on given range [lower, upper].
! 
!  BOXSPLEV -- Evaluate a box spline with a direction vector set
!              composed of <degree> repetitions of the identity
!              matrix supporting the box with corners [lower,upper].
!  
RECURSIVE FUNCTION CBSPLEV(X, DEGREE, LOWER, UPPER, IDX) RESULT(Y)
  ! Function for evaluating a cardinal B-spline of prescribed degree,
  ! which has evenly spaced knots located over [LOWER, UPPER] and
  ! nonzero support over the range (LOWER,UPPER).
  ! 
  ! Inputs:
  !   X      -- The (real-valued) location to evaluate a cardinal B-spline at.
  !   DEGREE -- The degree of the cardinal B-spline to evaluate.
  ! 
  ! Optional Input:
  !   LOWER  -- The location of the minimum-valued knot for this
  !             B-spline, default value of "0".
  !   UPPER  -- The location of the maximum-valued knot for this
  !             B-spline, default value of "DEGREE+1".
  !   IDX    -- A private input used for recursively evaluating the
  !             cardinal B-spline. *NOT FOR EXTERNAL USE*
  ! 
  ! Output:
  !   Y      -- The value of a cardinal B-spline of prescribed degree
  !             evaluated at the given position.
  ! 
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE

  ! Inputs
  REAL(KIND=REAL64), INTENT(IN)           :: X
  INTEGER,           INTENT(IN)           :: DEGREE
  REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: LOWER, UPPER  
  INTEGER,           INTENT(IN), OPTIONAL :: IDX
  ! Output
  REAL(KIND=REAL64) :: Y
  ! Local variables for optionals
  INTEGER :: INDEX
  REAL(KIND=REAL64) :: LOCAL_LOWER, LOCAL_UPPER, LOCAL_X

  LOCAL_X = X
  ! Initialize optional inputs
  DEFINE_LOWER : IF (PRESENT(LOWER)) THEN
     LOCAL_LOWER = LOWER
     ! Rescale X to be supported by the cardinal B-spline over [0, DEGREE+1]
     LOCAL_X = LOCAL_X - LOCAL_LOWER
  ELSE
     LOCAL_LOWER = 0
  END IF DEFINE_LOWER
  DEFINE_UPPER : IF (PRESENT(UPPER)) THEN
     LOCAL_UPPER = UPPER
     ! Rescale X to be supported by the cardinal B-spline over [0, DEGREE+1]
     LOCAL_X = (DEGREE+1) * LOCAL_X / (LOCAL_UPPER - LOCAL_LOWER)
  ELSE
     LOCAL_UPPER = DEGREE + 1
  END IF DEFINE_UPPER
  DEFINE_INDEX : IF (PRESENT(IDX)) THEN
     INDEX = IDX
  ELSE
     INDEX = 1
  END IF DEFINE_INDEX

  ! Body of cardinal B-spline evaluation function.
  IF (DEGREE .EQ. 0) THEN
     ! Base case, 1 for supported points.
     IF ((INDEX-1 .LT. LOCAL_X) .AND. (LOCAL_X .LE. INDEX)) THEN
        Y = 1
     ELSE
        Y = 0
     END IF
  ELSE
     ! Recursive evaluation of Cardinal B-spline
     Y = (LOCAL_X - (INDEX-1)) / (DEGREE) * &
          CBSPLEV(LOCAL_X, DEGREE-1, IDX=INDEX) + &
          (INDEX + DEGREE - LOCAL_X) / (INDEX+DEGREE - INDEX) * &
          CBSPLEV(LOCAL_X, DEGREE-1, IDX=INDEX+1)
  END IF
END FUNCTION CBSPLEV

FUNCTION CBOXSPLEV(X, DEGREE, LOWER, UPPER) RESULT(Y)
  ! Function for evaluating a product of axis-aligned cardinal B-splines
  ! of arbitrary degree at a point X in the box [LOWER,UPPER]^SIZE(X).
  ! This will be referred to as a "cardinal box spline".
  ! 
  ! Inputs:
  !   X      -- Point in R^d (d-dimensional real vector space) at
  !             which this box spline is evaluated.
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
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)           :: X
  INTEGER,           INTENT(IN)                         :: DEGREE
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:), OPTIONAL :: LOWER, UPPER
  ! Output
  REAL(KIND=REAL64) :: Y
  ! Local variables (iterator)
  INTEGER :: I
  REAL(KIND=REAL64), DIMENSION(SIZE(X)) :: LOCAL_LOWER, LOCAL_UPPER
  ! Defined interface for the CBSPLEV function being called internally.
  INTERFACE 
     RECURSIVE FUNCTION CBSPLEV(X, DEGREE, LOWER, UPPER, IDX) RESULT(Y)
       USE ISO_FORTRAN_ENV, ONLY: REAL64
       IMPLICIT NONE
       REAL(KIND=REAL64), INTENT(IN)           :: X
       INTEGER,           INTENT(IN)           :: DEGREE
       REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: LOWER, UPPER  
       INTEGER,           INTENT(IN), OPTIONAL :: IDX
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





! ====================================================================
!          Original C. De Boor code for evaluating a B-spline         
! ====================================================================
! 
! 
!       SUBROUTINE BSPLVB ( T, JHIGH, INDEX, X, LEFT, BIATX )
! CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT  X  OF ORDER
! C
! C               JOUT  =  MAX( JHIGH , (J+1)*(INDEX-1) )
! C
! C  WITH KNOT SEQUENCE  T .
! C
! C******  I N P U T  ******
! C  T.....KNOT SEQUENCE, OF LENGTH  LEFT + JOUT  , ASSUMED TO BE NONDE-
! C        CREASING.  A S S U M P T I O N . . . .
! C                       T(LEFT)  .LT.  T(LEFT + 1)   .
! C   D I V I S I O N  B Y  Z E R O  WILL RESULT IF  T(LEFT = T(LEFT+1)
! C  JHIGH,
! C  INDEX.....INTEGERS WHICH DETERMINE THE ORDER  JOUT = MAX(JHIGH,
! C        (J+1)*(INDEX-1))  OF THE B-SPLINES WHOSE VALUES AT  X  ARE TO
! C        BE RETURNED.  INDEX  IS USED TO AVOID RECALCULATIONS WHEN SEVE-
! C        RAL COLUMNS OF THE TRIANGULAR ARRAY OF B-SLPINE VALUES ARE NEE-
! C        DED (E.G., IN  BVALUE  OR IN  BSPLVD ). PRECISELY,
! C                     IF  INDEX = 1 ,
! C        THE CALCULATION STARTS FROM SCRATCH AND THE ENTIRE TRIANGULAR
! C        ARRAY OF B-SPLINE VALUES OF ORDERS 1,2,...,JHIGH  IS GENERATED
! C        ORDER BY ORDER , I.E., COLUMN BY COLUMN .
! C                     IF  INDEX = 2 ,
! C        ONLY THE B-SPLINE VALUES OF ORDER  J+1, J+2, ..., JOUT  ARE GE-
! C        NERATED, THE ASSUMPTION BEING THAT  BIATX , J , DELTAL , DELTAR
! C        ARE, ON ENTRY, AS THEY WERE ON EXIT AT THE PREVIOUS CALL.
! C           IN PARTICULAR, IF  JHIGH = 0, THEN  JOUT = J+1, I.E., JUST
! C        THE NEXT COLUMN OF B-SPLINE VALUES IS GENERATED.
! C
! C  W A R N I N G . . .  THE RESTRICTION   JOUT .LE. JMAX (= 20)  IS IM-
! C        POSED ARBITRARILY BY THE DIMENSION STATEMENT FOR  DELTAL  AND
! C        DELTAR  BELOW, BUT IS  N O W H E R E  C H E C K E D  FOR .
! C
! C  X.....THE POINT AT WHICH THE B-SPLINES ARE TO BE EVALUATED.
! C  LEFT.....AN INTEGER CHOSEN (USUALLY) SO THAT
! C                  T(LEFT) .LE. X .LE. T(LEFT+1)  .
! C
! C******  O U T P U T  ******
! C  BIATX.....ARRAY OF LENGTH  JOUT , WITH  BIATX(I)  CONTAINING THE VAL-
! C        UE AT  X  OF THE POLYNOMIAL OF ORDER  JOUT WHICH AGREES WITH
! C        THE B-SPLINE  B(LEFT-JOUT+I,JOUT,T)  ON THE INTERVAL (T(LEFT),
! C        T(LEFT+1)) .
! C
! C******  M E T H O D  ******
! C  THE RECURRENCE RELATION
! C
! C                       X - T(I)              T(I+J+1) - X
! C     B(I,J+1)(X)  =  -----------B(I,J)(X) + ---------------B(I+1,J)(X)
! C                     T(I+J)-T(I)            T(I+J+1)-T(I+1)
! C
! C  IS USED (REPEATEDLY) TO GENERATE THE (J+1)-VECTOR  B(LEFT-J,J+1)(X),
! C  ...,B(LEFT,J+1)(X)  FROM THE J-VECTOR  B(LEFT-J+1,J)(X),...,
! C  B(LEFT,J)(X), STORING THE NEW VALUES IN  BIATX  OVER THE OLD.  THE
! C  FACTS THAT
! C            B(I,1) = 1  IF  T(I) .LE. X .LT. T(I+1)
! C  AND THAT
! C            B(I,J)(X) = 0  UNLESS  T(I) .LE. X .LT. T(I+J)
! C  ARE USED. THE PARTICULAR ORGANIZATION OF THE CALCULATIONS FOLLOWS AL-
! C  GORITHM  (8)  IN CHAPTER X OF THE TEXT.
! C
!       INTEGER INDEX,JHIGH,LEFT,   I,J,JP1
!       REAL*8 BIATX(JHIGH),T(1),X,   DELTAL(20),DELTAR(20),SAVED,TERM
! C     DIMENSION BIATX(JOUT), T(LEFT+JOUT)
! CURRENT FORTRAN STANDARD MAKES IT IMPOSSIBLE TO SPECIFY THE LENGTH OF
! C  T  AND OF  BIATX  PRECISELY WITHOUT THE INTRODUCTION OF OTHERWISE
! C  SUPERFLUOUS ADDITIONAL ARGUMENTS.
!       DATA J/1/
! C     SAVE J,DELTAL,DELTAR (VALID IN FORTRAN 77)
! C
!                                         GO TO (10,20), INDEX
!    10 J = 1
!       BIATX(1) = 1.
!       IF (J .GE. JHIGH)                 GO TO 99
! C
!    20    JP1 = J + 1
!          DELTAR(J) = T(LEFT+J) - X
!          DELTAL(J) = X - T(LEFT+1-J)
!          SAVED = 0.
!          DO 26 I=1,J
!             TERM = BIATX(I)/(DELTAR(I) + DELTAL(JP1-I))
!             BIATX(I) = SAVED + DELTAR(I)*TERM
!    26       SAVED = DELTAL(JP1-I)*TERM
!          BIATX(JP1) = SAVED
!          J = JP1
!          IF (J .LT. JHIGH)              GO TO 20
! C
!    99                                   RETURN
!       END
