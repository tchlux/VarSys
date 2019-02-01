! This file (boxspline.f90) contains the subroutine BOXSPLEV that can
! be used to evaluate a box-spline, defined by its associated
! direction vector set, at given evaluation points, as well as the
! module REAL_PRECISION defining the real precision.
! 
! ====================================================================
MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

SUBROUTINE BOXSPLEV(DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  ! BOXSPLEV evaluates a box-spline, defined by a direction vector set in
  ! dimension S, at given evaluation points.
  ! 
  ! The implementation uses the numerically consistent algorithm for
  ! evaluating box-splines originally presented in [1]. Most notably,
  ! the evaluation of the box-spline near the boundaries of polynomial
  ! pieces does not exhibit the random behavior that is seen in the
  ! naive recursive implementation. Furthermore, the computational
  ! complexity for direction vector sets with repeated direction
  ! vectors is reduced from the naive recursive implementation.
  ! 
  ! [1] Kobbelt, Leif. "Stable Evaluation of Box-Splines." 
  !     Numerical Algorithms 14.4 (1997): 377-382.
  ! 
  ! The subroutine BOXSPLEV utilizes the LAPACK routines DGECON, DGELS,
  ! DGESVD, and DGETRF.
  ! 
  ! On input:
  ! 
  ! DVECS(:,:) is a real S x M array whose columns are the unique direction
  !    vectors used in defining the box-spline.  S <= M, and DVECS(1:S,1:S)
  !    must be invertible.
  !
  ! DVEC_MULTS(:) is an integer array of length M containing the multiplicity
  !    of each corresponding direction vector. 
  !
  ! EVAL_PTS(:,:) is a real S x L array whose columns are the points at which
  !    the box-spline is to be evaluated. 
  !
  ! On output:
  ! 
  ! BOX_EVALS(:) is a real array of length L containing the box-spline
  !    values at the evaluation points.
  !
  ! ERROR is an integer error flag of the form
  !             100*(LAPACK INFO flag) + 10*T + U.
  !    ERROR .EQ. 0 is a normal return.  For ERROR .NE. 0, the 
  !    meanings of the tens digit T and the units digit U are:
  !
  !    Tens digit T:
  !    0  Improper usage.
  !    1  Error computing normal vectors.
  !    2  Error computing box-spline.
  !    3  Error computing an orthogonal vector.
  !    4  Error computing a minimum norm representation.
  !    5  Error computing a matrix rank.
  !    6  Error computing a matrix determinant.
  !    7  Error finding nonzero entries in an array.
  !    8  Error regarding memory.
  ! 
  !    Units digit U, for T = 0:
  !    1  Mismatched dimension, SIZE(DVEC_MULTS) .NE. SIZE(DVECS,2).
  !    2  Mismatched dimension, SIZE(EVAL_PTS,1) .NE. SIZE(DVECS,1).
  !    3  Mismatched dimension, SIZE(BOX_EVALS) .NE. SIZE(EVAL_PTS,2).
  !    4  One of the multiplicity values provided was < 1.
  !    5  Columns of DVECS are not unique.
  !    6  M < S or DVECS(1:S,1:S) is near singular.
  ! 
  !    Units digit U, for T /= 0:
  !    0  Work array allocation failed.
  !    1  Work array size query failed.
  !    2  DGELS computation error, see 'INFO' for details.
  !    3  DGESVD computation error, see 'INFO' for details.
  !    4  DGETRF computation error, see 'INFO' for details.
  ! 
  ! The calling program should include the following interface to BOXSPLEV:
  ! INTERFACE 
  !   SUBROUTINE BOXSPLEV(DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  !     USE REAL_PRECISION, ONLY: R8
  !     REAL(KIND=R8), INTENT(IN), DIMENSION(:,:):: DVECS
  !     INTEGER, INTENT(IN), DIMENSION(:):: DVEC_MULTS
  !     REAL(KIND=R8), INTENT(IN), DIMENSION(:,:):: EVAL_PTS
  !     REAL(KIND=R8), INTENT(OUT), DIMENSION(:):: BOX_EVALS
  !     INTEGER, INTENT(OUT):: ERROR
  !   END SUBROUTINE BOXSPLEV
  ! END INTERFACE

  USE REAL_PRECISION, ONLY: R8
  IMPLICIT NONE
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: DVECS
  INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
  INTEGER,       INTENT(OUT)                 :: ERROR
  ! Local variables.
  INTEGER :: DIM, IDX_1, IDX_2, NUM_DVECS, INFO
  INTEGER,       DIMENSION(SIZE(DVECS,2))                  :: LOCATION, LOOKUP
  INTEGER,       DIMENSION(SIZE(DVECS,2),SIZE(DVECS,2))    :: IDENTITY
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1),2**SIZE(DVECS,2)) :: NORMAL_VECTORS
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE               :: MIN_NORM_DVECS
  REAL(KIND=R8), DIMENSION(:),   ALLOCATABLE               :: LAPACK_WORK
  REAL(KIND=R8), PARAMETER :: SQRTEPS = SQRT(EPSILON(1.0_R8))
  REAL(KIND=R8) :: RCOND

  ERROR = 0 ! Initialize error flag.
  ! Allocate a work array large enough for all LAPACK calls.
  LAPACK_WORK = ALLOCATE_MAX_LAPACK_WORK(SIZE(DVECS,1), SIZE(DVECS,2))
  IF (ERROR .NE. 0) RETURN
  ! Store 'global' constants for box-spline evaluation.
  DIM = SIZE(DVECS, 1)
  NUM_DVECS = SIZE(DVECS, 2)
  ! Create an identity matrix (for easy matrix-vector multiplication)
  ! and compute the index lookup indices (for binary style hashing).
  IDENTITY = 0
  FORALL (IDX_1 = 1:NUM_DVECS)
     IDENTITY(IDX_1,IDX_1) = 1
     LOOKUP(IDX_1) = 2**(IDX_1-1)
  END FORALL

  ! Check for usage errors.
  ! Check for dimension mismatches.
  IF (SIZE(DVEC_MULTS) .NE. SIZE(DVECS,2))    THEN; ERROR = 1; RETURN; END IF
  IF (SIZE(EVAL_PTS,1) .NE. SIZE(DVECS,1))    THEN; ERROR = 2; RETURN; END IF
  IF (SIZE(BOX_EVALS)  .NE. SIZE(EVAL_PTS,2)) THEN; ERROR = 3; RETURN; END IF
  ! Check for invalid multiplicity.
  IF (MINVAL(DVEC_MULTS) .LT. 1) THEN; ERROR = 4; RETURN; END IF
  ! Check uniqueness of DVECS columns.
  DO IDX_1 = 1, NUM_DVECS-1
     DO IDX_2 = IDX_1+1, NUM_DVECS
        IF (SUM(ABS(DVECS(:,IDX_1) - DVECS(:,IDX_2))) .LT. SQRTEPS) THEN
           ERROR = 5; RETURN
        END IF
     END DO
  END DO
  ! Check if NUM_DVECS < DIM or rank DVECS(1:DIM,1:DIM) < DIM.
  IF (SIZE(DVECS,2) .LT. SIZE(DVECS,1)) THEN 
     ERROR = 6; RETURN
  ELSE   
     ! Compute condition number COND(DVECS(1:DIM,1:DIM)) and test for near
     ! rank deficiency: 1/COND(DVECS(1:DIM,1:DIM)) < SQRT(EPSILON(1.0_R8)).
     RCOND = MATRIX_CONDITION(DVECS(1:SIZE(DVECS,1), 1:SIZE(DVECS,1)))
     IF (RCOND .LT. SQRTEPS) THEN
        ERROR = 6; RETURN
     ENDIF
  ENDIF

  ! Get the minimum norm representation of the direction vectors.
  MIN_NORM_DVECS = MATRIX_MINIMUM_NORM(DVECS)
  IF (ERROR .NE. 0) RETURN
  
  ! Compute the table of normal vectors defining boundaries of
  ! polynomial pieces of box-spline.
  LOCATION = 0
  NORMAL_VECTORS = 0
  CALL COMPUTE_NORMALS(DIM-1, NUM_DVECS, LOCATION, NORMAL_VECTORS)
  IF (ERROR .NE. 0) RETURN

  ! Recursively evaluate the box-spline.
  LOCATION = 0.
  BOX_EVALS = 0.
  CALL EVALUATE_BOX_SPLINE(DVEC_MULTS, LOCATION, MIN_NORM_DVECS, &
       MATMUL(EVAL_PTS, MIN_NORM_DVECS), BOX_EVALS)
  IF (ERROR .NE. 0) RETURN

CONTAINS

  ! ==================================================================
  RECURSIVE SUBROUTINE COMPUTE_NORMALS(DIM, NUM, LOC, NORMAL_VECTORS_SUBSET)
    ! 1) COMPUTE_NORMALS:
    ! 
    !   Evaluate and store the normal vectors that define polynomial
    !   piece boundaries for the box-spline with given direction
    !   vector set. (O(2^"NUM") calls in recursive subtree).
    ! 
    ! Inputs:
    !   DIM, NUM -- Integers used for proceeding with recursion.
    !   LOC      -- Array of integers storing binary-style 
    !               representation of active unique direction vectors.
    ! 
    ! Output:
    !   NORMAL_VECTORS_SUBSET -- The subset of real-valued normal
    !                            column vectors being computed by this
    !                            branch of the recursive tree. At base
    !                            case there is only one column vector.
    ! 
    INTEGER,           INTENT(IN)                    :: DIM, NUM
    INTEGER,           INTENT(IN),    DIMENSION(:)   :: LOC
    REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: NORMAL_VECTORS_SUBSET
    ! Handle recursion
    IF (DIM .GT. 0) THEN
       ! Get the left and right sub-trees (left is 0 for current bit,
       ! right is 1 for current bit) and store evaluated subtrees.
       ! 
       ! Only compute the left subtree if it will not go out of bounds.
       left_subtree : IF (NUM .GT. DIM) THEN
          CALL COMPUTE_NORMALS(DIM, NUM-1, LOC, &
               NORMAL_VECTORS_SUBSET(:,:1+2**(NUM-1)))
          IF (ERROR .NE. 0) RETURN
       END IF left_subtree
       ! Always compue the right subtree, because bounds are size-constrained.
       CALL COMPUTE_NORMALS(DIM-1, NUM-1, LOC+IDENTITY(:,NUM), &
            NORMAL_VECTORS_SUBSET(:,1+2**(NUM-1):))
    ELSE
       ! Set the column vector to be a vector orthogonal to all
       ! selected vectors of the box-spline direction vector set.
       CALL MATRIX_ORTHOGONAL(TRANSPOSE(DVECS(:,NONZERO(LOC))), &
            NORMAL_VECTORS_SUBSET(:,1))
    END IF
  END SUBROUTINE COMPUTE_NORMALS
  
  ! ==================================================================
  RECURSIVE SUBROUTINE EVALUATE_BOX_SPLINE(MULTS, LOC, SUB_DVECS,&
       SHIFTED_EVAL_PTS, EVALS_AT_PTS)
    ! 2) EVALUATE_BOX_SPLINE:
    !   
    !   Evaluate the box-spline defined by "DVECS" recursively, where
    !   this iteration has the remaining vectors "SUB_DVECS", at all
    !   points in "SHIFTED_EVAL_PTS" and store box-spline evaluations
    !   in "EVALS_AT_PTS".
    ! 
    ! Inputs:
    !   MULTS            -- Integer array of direction vector multiplicities.
    !   LOC              -- Integer array tracking current recursion position.
    !   SUB_DVECS        -- Real matrix of (subset of) (transformed)
    !                       direction vectors used to evaluated at
    !                       provided evaluation points.
    !   SHIFTED_EVAL_PTS -- Real matrix of (shifted) evaluation points
    !                       at which the box-spline value will be computed.
    ! 
    ! Output:
    !   EVALS_AT_PTS -- Real array of box-spline evaluations at each
    !                   corresponding (shifted) evaluation point.
    ! 
    INTEGER,           INTENT(IN), DIMENSION(NUM_DVECS) :: MULTS, LOC
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:)       :: SUB_DVECS
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:)       :: SHIFTED_EVAL_PTS
    REAL(KIND=R8), INTENT(OUT), &
         DIMENSION(SIZE(EVAL_PTS,1)) :: EVALS_AT_PTS
    ! Temporary variables for recursive computations of eval points,
    ! shifted evaluation points, and the actual evaluations of box splines.
    REAL(KIND=R8), DIMENSION(SIZE(SHIFTED_EVAL_PTS,1),&
         SIZE(SHIFTED_EVAL_PTS,2))                   :: TEMP_SHIFTED_EVAL_PTS
    REAL(KIND=R8), DIMENSION(SIZE(EVAL_PTS,1),&
         SIZE(EVAL_PTS,2))                           :: TEMP_EVAL_PTS
    REAL(KIND=R8), DIMENSION(SIZE(EVALS_AT_PTS)) :: TEMP_EVALS_AT_PTS
    REAL(KIND=R8), DIMENSION(SIZE(EVAL_PTS,1))   :: LOCATIONS
    REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE   :: NEXT_DVECS
    INTEGER,           DIMENSION(NUM_DVECS)          :: NEXT_MULTS, NEXT_LOC
    INTEGER,           DIMENSION(:),   ALLOCATABLE   :: REMAINING_DVECS
    REAL(KIND=R8), DIMENSION(DIM) :: PT_SHIFT
    REAL(KIND=R8) :: POSITION, 
    INTEGER :: IDX_1, IDX_2, IDX_3

    ! Recursion case ...
    IF (SUM(MULTS) > DIM) THEN
       EVALS_AT_PTS = 0.
       IDX_2 = 1
       ! Sum over all direction vectors.
       DO IDX_1 = 1, NUM_DVECS
          ! Update multiplicity of directions and position in recursion tree.
          NEXT_MULTS = MULTS - IDENTITY(:, IDX_1) ! Reduce multiplicity.
          NEXT_LOC   = LOC   + IDENTITY(:, IDX_1) ! Track recursion position.
          ! Make recursive calls.
          IF (MULTS(IDX_1) .GT. 1) THEN
             ! Perform recursion with only reduced multiplicity.
             CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, SUB_DVECS, &
                  SHIFTED_EVAL_PTS, TEMP_EVALS_AT_PTS)
             IF (ERROR .NE. 0) RETURN
             EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                  SHIFTED_EVAL_PTS(:,IDX_2)
             ! Perform recursion with transformed set of direction
             ! vectors and evaluation points.
             compute_shift_1 : FORALL (IDX_3 = 1:SIZE(SUB_DVECS,2))
                PT_SHIFT(1) = SUM(DVECS(:,IDX_1) * SUB_DVECS(:,IDX_3))
                TEMP_EVALS_AT_PTS(:,IDX_3) = TEMP_EVALS_AT_PTS(:,IDX_3) - PT_SHIFT(1)
             END FORALL compute_shift_1
             CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, NEXT_LOC, SUB_DVECS, &
                  TEMP_SHIFTED_EVAL_PTS, TEMP_EVALS_AT_PTS)
             IF (ERROR .NE. 0) RETURN
             EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                  (MULTS(IDX_1) - SHIFTED_EVAL_PTS(:,IDX_2))
             IDX_2 = IDX_2 + 1
          ELSE IF (MULTS(IDX_1) .GT. 0) THEN
             ! Find the next direction vectors (ones with nonzero multiplicities).
             NEXT_DVECS = DVECS(:,NONZERO(NEXT_MULTS))
             IF (ERROR .NE. 0) RETURN
             IF (MATRIX_RANK(TRANSPOSE(NEXT_DVECS)) .EQ. DIM) THEN
                IF (ERROR .NE. 0) RETURN
                ! Update Least norm representation
                NEXT_DVECS = MATRIX_MINIMUM_NORM(NEXT_DVECS)
                IF (ERROR .NE. 0) RETURN
                ! Perform recursion with only reduced multiplicity.
                compute_shift_2 : FORALL (IDX_3 = 1:DIM)
                   PT_SHIFT(1) = SUM(LOC * DVECS(IDX_3,:))
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - PT_SHIFT(1)
                END FORALL compute_shift_2
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     SHIFTED_EVAL_PTS(:,IDX_2)
                ! Perform recursion with transformed set of direction vectors.
                compute_shift_3 : FORALL (IDX_3 = 1:DIM)
                   PT_SHIFT(1) = SUM(NEXT_LOC * DVECS(IDX_3,:))
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - PT_SHIFT(1)
                END FORALL compute_shift_3
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, NEXT_LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     (MULTS(IDX_1) - SHIFTED_EVAL_PTS(:,IDX_2))
             END IF
             IDX_2 = IDX_2 + 1
          END IF
       END DO
       ! Normalize by number of direction vectors in computation.
       EVALS_AT_PTS = EVALS_AT_PTS / (SUM(MULTS) - DIM)
    ELSE
       ! Base case ... compute characteristic function.
       EVALS_AT_PTS = 1.
       ! Delayed translations (this is what makes the algorithm more stable).
       REMAINING_DVECS = NONZERO(MULTS)
       IF (ERROR .NE. 0) RETURN
       compute_shift_4 : FORALL (IDX_1 = 1:DIM)
          PT_SHIFT(1) = SUM(LOC * DVECS(IDX_1,:))
          TEMP_EVAL_PTS(:,IDX_1) = EVAL_PTS(:,IDX_1) - PT_SHIFT(1)
       END FORALL compute_shift_4
       ! Check against all *precomputed* hyperplanes (also contributes to stability).
       DO IDX_1 = 1, DIM
          ! Lookup normal vector to current hyperplane.
          IDX_3 = 1 + SUM(LOOKUP * ( &
               MULTS - IDENTITY(:,REMAINING_DVECS(IDX_1)) ))
          ! Compute shifted position (relative to normal ector).
          POSITION = SUM(DVECS(:,REMAINING_DVECS(IDX_1)) * NORMAL_VECTORS(:,IDX_3))
          ! Compute shifted evaluation locations. (NUM_PTS, DIM) x (DIM, 1)
          compute_shifted_point_1 : FORALL (IDX_2 = 1:SIZE(EVAL_PTS,1))
             LOCATIONS(IDX_2) = SUM(TEMP_EVAL_PTS(IDX_2,:) * NORMAL_VECTORS(:,IDX_3))
          END FORALL compute_shifted_point_1
          ! Identify those points that are outside of this box (0-side).
          IF (POSITION .GT. 0) THEN
             WHERE (LOCATIONS .LT. 0) EVALS_AT_PTS = 0.
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (LOCATIONS .GE. 0) EVALS_AT_PTS = 0.
          END IF
          ! Recompute shifted locations based on remaining direction vectors.
          NEXT_LOC = LOC + IDENTITY(:,REMAINING_DVECS(IDX_1))
          compute_shift_5 : FORALL (IDX_2 = 1:DIM)
             PT_SHIFT(IDX_2) = SUM(NEXT_LOC * DVECS(IDX_2,:))
          END FORALL compute_shift_5
          compute_shifted_point_2 : FORALL (IDX_2 = 1:SIZE(EVAL_PTS,1))
             LOCATIONS(IDX_2) = SUM((EVAL_PTS(IDX_2,:) - PT_SHIFT) * &
                  NORMAL_VECTORS(:,IDX_3))
          END FORALL compute_shifted_point_2
          ! Identify those shifted points that are outside of this box (DVEC-side).
          IF (POSITION .GT. 0) THEN
             WHERE (LOCATIONS .GE. 0) EVALS_AT_PTS = 0.
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (LOCATIONS .LT. 0) EVALS_AT_PTS = 0.
          END IF
       END DO
       ! Normalize evaluations by determinant of box.
       EVALS_AT_PTS = EVALS_AT_PTS / ABS( &
            MATRIX_DET(TRANSPOSE(DVECS(:,REMAINING_DVECS(1:DIM)))) )
       IF (ERROR .NE. 0) RETURN
    END IF
  END SUBROUTINE EVALUATE_BOX_SPLINE

  !===============================================================
  !             Mathematical Convenience Operations     
  !===============================================================

  ! ==================================================================
  SUBROUTINE MATRIX_ORTHOGONAL(A, ORTHOGONAL)
    ! 3) MATRIX_ORTHOGONAL
    ! 
    !   Given a matrix A of row vectors, compute a vector orthogonal to
    !   the row vectors in A and store it in ORTHOGONAL using DGESVD.
    ! 
    ! Input:
    !   A -- Real dense matrix.
    ! 
    ! Output:
    !   ORTHOGONAL -- Real vector orthogonal to given matrix.
    ! 
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:)       :: A
    REAL(KIND=R8), INTENT(OUT), DIMENSION(SIZE(A,2)) :: ORTHOGONAL
    ! Local variables for computing the orthogonal vector
    REAL(KIND=R8), DIMENSION(MIN(SIZE(A,1),SIZE(A,2))) :: S
    REAL(KIND=R8), DIMENSION(SIZE(A,2),SIZE(A,2)) :: VT
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: IDX
    LOGICAL :: FOUND_ZERO
    ! Unused parameters
    REAL(KIND=R8), DIMENSION(1) :: U

    ! Query the size of the work array to construct.
    CALL DGESVD('N', 'A', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), S, &
         U, SIZE(U,1), VT, SIZE(VT,1), U, -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 81 + 100*INFO; RETURN; END IF
    ! Allocate the work array.
    ALLOCATE(WORK(INT(U(1))))

    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N','A',SIZE(A,1),SIZE(A,2),A,SIZE(A,1),S,U,SIZE(U), &
         VT, SIZE(VT,1), WORK, SIZE(WORK), INFO)
    !   'N'        -- No columns of U are computed
    !   'A'        -- All N rows of V**T are returned in the array VT
    !   SIZE(A,1)  -- M, number of rows in A
    !   SIZE(A,2)  -- N, number of cols in A
    !   A          -- Inputs matrix A (destroyed during computation)
    !   SIZE(A,1)  -- Leading dimension of A
    !   S          -- Container for singular values of A
    !   U          -- Double precision array for storing U
    !   SIZE(U,1)  -- Leading dimension of U
    !   VT         -- Double precision array for storing V^T
    !   SIZE(VT,1) -- Leading dimension of V^T
    !   WORK       -- Work array for computing SVD
    !   SIZE(WORK) -- Size of the work array
    !   INFO       -- Info message parameter

    IF (INFO .NE. 0) THEN; ERROR = 33 + 100*INFO; RETURN; END IF

    ORTHOGONAL = 0.
    FOUND_ZERO = .FALSE.
    ! Find the first vector in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value).
    find_null : DO IDX = 1, SIZE(S)
       IF (S(IDX) .LE. (EPSILON(U(1)) * SIZE(S))) THEN
          ! If we found a singular value, copy out the orthogonal vector.
          ORTHOGONAL = VT(IDX,:)
          FOUND_ZERO = .TRUE.
          EXIT find_null
       END IF
    END DO find_null
    ! If no orthogonal was found and the matrix does not contain
    ! enough vectors to span the space, use the first vector in VT at
    ! (RANK(A) + 1).
    IF ((SIZE(VT,1) > SIZE(S)) .AND. (.NOT. FOUND_ZERO))THEN
       ORTHOGONAL = VT(SIZE(S)+1,:)
    END IF
  END SUBROUTINE MATRIX_ORTHOGONAL

  ! ==================================================================
  FUNCTION MATRIX_MINIMUM_NORM(MATRIX) RESULT(MIN_NORM)
    ! 4) MATRIX_MINIMUM_NORM
    ! 
    !   Compute the minimum norm representation of 'MARTIX' and store
    !   it in 'MIN_NORM', use DGELS to find the least squares
    !   solution to the problem (AX = I).
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   MIN_NORM -- Real dense matrix that is the minimum norm
    !               representation of MATRIX.
    ! 
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: MIN_NORM
    ! Local variables
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: DGELS_WORK_ARRAY
    REAL(KIND=R8), DIMENSION(1) :: DGELS_SIZE_HOLDER
    INTEGER :: DGELS_DIM
    ! Allocate the output matrix.
    ALLOCATE(MIN_NORM(1:SIZE(MATRIX,1),1:SIZE(MATRIX,2)))
    ! Make "MIN_NORM" the identity matrix
    MIN_NORM = 0.
    FORALL (DGELS_DIM = 1 : SIZE(MIN_NORM
    ! Get the size of the work array necessary.
    CALL DGELS('N', SIZE(MATRIX,1), SIZE(MATRIX,2), SIZE(MIN_NORM,2),&
         MATRIX, SIZE(MATRIX,1), MIN_NORM, SIZE(MIN_NORM,1), &
         DGELS_SIZE_HOLDER, -1, INFO)

    IF (INFO .NE. 0) THEN; ERROR = 42; RETURN; END IF

    DGELS_DIM = DGELS_SIZE_HOLDER(1)
    ALLOCATE( DGELS_WORK_ARRAY(1:DGELS_DIM) )

    ! Call DGELS for actual solve.
    CALL DGELS('N', SIZE(MATRIX,1), SIZE(SQUARE,2), SIZE(MIN_NORM,2),&
         SQUARE, SIZE(SQUARE,1), MIN_NORM, SIZE(MIN_NORM,1), &
         DGELS_WORK_ARRAY, DGELS_DIM, INFO)
    !   'N'              -- No transpose of A is necessary
    !   SIZE(SQUARE,1)   -- number of rows in A
    !   SIZE(SQUARE,2)   -- number of columns in A
    !   SIZE(MIN_NORM,2) -- Number of right hand sides
    !   SQUARE           -- A
    !   SIZE(SQUARE,1)   -- Leading dimension of A (number of rows)
    !   MIN_NORM         -- B (will be overwritten to hold X after call)
    !   SIZE(MIN_NORM,1) -- Leading dimension of B (number of rows)
    !   DGELS_WORK_ARRAY -- Workspace array for DGELS
    !   DGELS_DIM        -- Size of dgels_work_array
    !   INFO             -- For verifying successful execution

    IF (INFO .NE. 0) THEN; ERROR = 43; RETURN; END IF
  END FUNCTION MATRIX_MINIMUM_NORM

  ! ==================================================================
  FUNCTION MATRIX_RANK(MATRIX)
    ! 5) MATRIX_RANK
    ! 
    !   Get the rank of the provided matrix using the SVD.
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   MATRIX_RANK -- The integer rank of MATRIX.
    ! 
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    ! Local variables for computing the orthogonal vector.
    REAL(KIND=R8), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: A
    REAL(KIND=R8), DIMENSION(MIN(SIZE(MATRIX,1),SIZE(MATRIX,2))) :: S
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: IDX, MATRIX_RANK
    ! Unused DGESVD parameters.
    REAL(KIND=R8), DIMENSION(1) :: U, VT

    A = MATRIX
    ! Query the size of the work array to construct.
    CALL DGESVD('N', 'N', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), S, &
         U, SIZE(U), VT, SIZE(VT), U, -1, INFO)

    IF (INFO .NE. 0) THEN; ERROR = 54; RETURN; END IF

    ! Allocate the work array.
    ALLOCATE(WORK(INT(U(1))))

    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N','N',SIZE(A,1),SIZE(A,2),A,SIZE(A,1),S,U,SIZE(U), &
         VT, SIZE(VT), WORK, SIZE(WORK), INFO)
    !   'N'        -- No columns of U are computed
    !   'N'        -- No vectors in V**T are computed
    !   SIZE(A,1)  -- M, number of rows in A
    !   SIZE(A,2)  -- N, number of cols in A
    !   A          -- Inputs matrix A
    !   SIZE(A,1)  -- Leading dimension of A
    !   S          -- Container for singular values of A
    !   U          -- Double precision array for storing U
    !   SIZE(U,1)  -- Leading dimension of U
    !   VT         -- Double precision array for storing V^T
    !   SIZE(VT)   -- Leading dimension of V^T
    !   WORK       -- Work array for computing SVD
    !   SIZE(WORK) -- Size of the work array
    !   INFO       -- Info message parameter

    IF (INFO .NE. 0) THEN; ERROR = 55; RETURN; END IF

    MATRIX_RANK = SIZE(S)
    ! Find the first singular value in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value)
    find_null : DO IDX = 1, SIZE(S)
       IF (S(IDX) .LE. (EPSILON(U(1)) * MIN(SIZE(A,1),SIZE(A,2)))) THEN
          MATRIX_RANK = IDX - 1
          EXIT find_null
       END IF
    END DO find_null
  END FUNCTION MATRIX_RANK

  ! ==================================================================
  FUNCTION MATRIX_DET(MATRIX)
    ! 6) MATRIX_DET
    ! 
    !   Compute the determinant of a matrix without modifying it
    !   using the LU decomposition (and a copy).
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   MATRIX_DET -- The real-valued determinant of MATRIX.
    ! 
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=R8) :: MATRIX_DET
    ! Local variable
    REAL(KIND=R8), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: M
    INTEGER, DIMENSION(MIN(SIZE(MATRIX,1), SIZE(MATRIX,2))) :: IPIV
    INTEGER :: IDX
    ! Copy into a local matrix.
    M = MATRIX
    ! Do the LU decomposition.
    CALL DGETRF(SIZE(M,1), SIZE(M,2), M, SIZE(M,1), IPIV, INFO)

    IF (INFO .NE. 0) THEN
       MATRIX_DET = 1.
       ERROR = 66
       RETURN
    END IF

    ! Compute the determinant (product of diagonal of U).
    MATRIX_DET = 1.
    DO IDX = 1, MIN(SIZE(M,1), SIZE(M,2))
       MATRIX_DET = MATRIX_DET * M(IDX,IDX)
    END DO
  END FUNCTION MATRIX_DET

  ! ==================================================================
  FUNCTION MATRIX_CONDITION(MATRIX) RESULT(COND)
    ! 7) MATRIX_CONDITION
    ! 
    ! Compute the condition number (for testing near rank deficiency)
    ! using DGETRF and DGECON.
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   COND -- Real value corresponding to the condition number of
    !           the provided matrix.
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=R8) :: COND
    ! Use LAPACK

  END FUNCTION MATRIX_CONDITION

  ! ==================================================================
  FUNCTION NONZERO(ARRAY) RESULT(NE_ZERO)
    ! 7) NONZERO
    ! 
    ! Return a new array of the indices of 'ARRAY' that contain
    ! nonzero elements. Set error and return array with 1 if ALL(ARRAY==0).
    ! 
    ! Input:
    !   ARRAY -- Real array of numbers.
    ! 
    ! Output:
    !   NE_ZERO -- Integer array corresponding to those indices of
    !              ARRAY that contain nonzero elements.
    ! 
    INTEGER, INTENT(IN), DIMENSION(:) :: ARRAY
    INTEGER, DIMENSION(:), ALLOCATABLE :: NE_ZERO
    ! Local variables
    INTEGER, DIMENSION(SIZE(ARRAY)) :: INDICES
    INTEGER :: COUNT_NONZERO, IDX
    ! Identify nonzero elements of ARRAY
    COUNT_NONZERO = 0
    DO IDX = 1, SIZE(ARRAY)
       IF (ARRAY(IDX) .NE. 0) THEN
          COUNT_NONZERO = COUNT_NONZERO + 1
          INDICES(COUNT_NONZERO) = IDX
       END IF
    END DO
    IF (COUNT_NONZERO .LE. 0) THEN
       ERROR = 70
       ALLOCATE(NE_ZERO(1))
       NE_ZERO(1) = 1
    ELSE
       ! Allocate the smaller output array and copy in the values
       ALLOCATE(NE_ZERO(1:COUNT_NONZERO))
       NE_ZERO = INDICES(:COUNT_NONZERO)
    END IF
  END FUNCTION NONZERO

END SUBROUTINE BOXSPLEV
