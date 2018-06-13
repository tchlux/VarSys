! This file (boxspline.f90) contains the subroutine BOXSPLEV that can
! be used to evaluate a box-spline defined by its associated
! direction vector set at provided evaluation points.
! 
! The BOXSPLEV subroutine utilizes the following LAPACK routines:
! 
!     DGELS   -- Computing a minimum norm representation with LS.
!     DGESVD  -- Computing SVD to find matrix rank / an orthogonal vector.
!     DGETRF  -- Computing the determinant of matrix.
! 
! The implementation uses the numerically consistent algorithm for
! evaluating box-splines originally presented in [1]. Most notably,
! the evaluation of the box-spline near the boundaries of polynomial
! pieces does *not* exhibit the random behavior that is seen in the
! naive recursive implementation. Furthermore, the computational
! complexity for direction vector sets with repeated direction
! vectors is reduced from the naive recursive implementation.
! 
! [1] Kobbelt, Leif. "Stable evaluation of box‐splines." 
!     Numerical Algorithms 14.4 (1997): 377-382.
! 
! ====================================================================
SUBROUTINE BOXSPLEV(DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR, INFO)
  ! Subroutine for evaluating a box-spline defined by a direction
  ! vector set at provided evaluation points.
  ! 
  ! Inputs:
  !   DVECS      -- Double precision real column vectors providing the
  !                 *unique* direction vectors that define this box-spline.
  !   DVEC_MULTS -- Integer array containing the multiplicity
  !                 of each corresponding direction vector, of length
  !                 SIZE(DVECS,2).
  !   EVAL_PTS   -- Double precision real row vectors providing the
  !                 points at which the box-spline is evaluated. 
  !                 SIZE(EVAL_PTS,2) = SIZE(DVECS,1)
  ! 
  ! Outputs:
  !   BOX_EVALS -- Double precision real array of box-spline
  !                evaluations at each corresponding evaluation point.
  !                SIZE(BOX_EVALS) = SIZE(EVAL_PTS,1)
  !   ERROR     -- Integer error flag with corresponding meanings
  !                listed under "Error flags" section. Any error at or
  !                above 10 is likely caused by hardware limitations.
  !   INFO      -- For all ** error codes relating to LAPACK, this
  !                integer contains the associated "INFO" provided by
  !                LAPACK subroutine evaluation.
  ! 
  ! Error flags (asterisks represent nonzero digits):
  ! 
  !       00 -> Successful execution.
  !       0* -> Improper usage.
  !       1* -> Error computing normal vectors.
  !       2* -> Error computing box-spline.
  !       3* -> Error computing a orthogonal vector.
  !       4* -> Error computing a minimum norm representation.
  !       5* -> Error computing a matrix rank.
  !       6* -> Error computing a matrix determinant.
  !       7* -> Error finding nonzero entries in an array.
  ! 
  !       Least significant digit meanings for errors "0*":
  !         01 - Mismatched dimension, SIZE(EVAL_PTS,2) /= SIZE(DVECS,1).
  !         02 - One of the multiplicity values provided was < 1.
  !         03 - At least 1 non-unique direction vector was provided.
  ! 
  !       Least significant digit meanings for errors "**":
  !         *1 - Memory allocation error (failed allocation).
  !         *2 - DGELS allocation error, see 'INFO' for details.
  !         *3 - DGELS computatoin error, see 'INFO' for details.
  !         *4 - DGESVD allocation error, see 'INFO' for details.
  !         *5 - DGESVD computation error, see 'INFO' for details.
  !         *6 - DGETRF computation error, see 'INFO' for details.
  ! 
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  ! Inputs and outputs
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)              :: DVECS
  INTEGER,           INTENT(IN),  DIMENSION(SIZE(DVECS,2))    :: DVEC_MULTS
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)              :: EVAL_PTS
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(EVAL_PTS,1)) :: BOX_EVALS
  INTEGER,           INTENT(OUT)                              :: ERROR, INFO
  ! Local variables
  INTEGER :: DIM, NUM_DVECS, IDX_1, IDX_2
  INTEGER, DIMENSION(SIZE(DVECS,2))                            :: LOCATION
  REAL(KIND=REAL64), DIMENSION(SIZE(DVECS,2))                  :: LOOKUP
  INTEGER,           DIMENSION(SIZE(DVECS,2),SIZE(DVECS,2))    :: IDENTITY
  REAL(KIND=REAL64), DIMENSION(SIZE(DVECS,1),2**SIZE(DVECS,2)) :: NORMAL_VECTORS
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE               :: MIN_NORM_DVECS

  ! Initialize error flag and info flag
  INFO = 0
  ERROR = 0
  ! Store 'global' constants for box-spline evaluation.
  DIM = SIZE(DVECS, 1)
  NUM_DVECS = SIZE(DVECS, 2)
  ! Create an identity matrix (for easy matrix-vector multiplication)
  ! and compute the index lookup indices (for binary style hashing).
  IDENTITY = 0
  init_ident_lookup : DO IDX_1 = 1, NUM_DVECS
     IDENTITY(IDX_1,IDX_1) = 1
     LOOKUP(IDX_1) = 2**(IDX_1-1)
  END DO init_ident_lookup

  ! Error checking
  mismatched_dim_check : IF (SIZE(DVECS,1) .NE. SIZE(EVAL_PTS,2)) THEN
     ERROR = 1
     RETURN
  END IF mismatched_dim_check
  bad_multiplicity_check : IF (MINVAL(DVEC_MULTS) .LT. 1) THEN
     ERROR = 2
     RETURN
  END IF bad_multiplicity_check
  nonunique_dvec_check : DO IDX_1 = 1, NUM_DVECS
     DO IDX_2 = IDX_1+1, NUM_DVECS
        IF (SUM(ABS(DVECS(:,IDX_1) - DVECS(:,IDX_2))) .LT. &
             SQRT(EPSILON(DVECS(1,1)))) THEN
           ERROR = 3
           RETURN
        END IF
     END DO
  END DO nonunique_dvec_check
  
  ! Get the minimum norm representation of the direction vectors
  MIN_NORM_DVECS = MINIMUM_NORM_REPR(DVECS)
  failed_min_check : IF (ERROR .NE. 0) THEN
     RETURN
  END IF failed_min_check
  
  ! Compute the table of normal vectors defining boundaries of
  ! polynomial pieces of box-spline.
  LOCATION = 0
  NORMAL_VECTORS = 0
  CALL COMPUTE_NORMALS(DIM-1, NUM_DVECS, LOCATION, NORMAL_VECTORS)
  failed_norm_check : IF (ERROR .NE. 0) THEN
     RETURN
  END IF failed_norm_check

  ! Recursive evaluation of box-spline
  LOCATION = 0.
  BOX_EVALS = 0.
  CALL EVALUATE_BOX_SPLINE(DVEC_MULTS, LOCATION, MIN_NORM_DVECS, &
       MATMUL(EVAL_PTS, MIN_NORM_DVECS), BOX_EVALS)
  failed_rec_check : IF (ERROR .NE. 0) THEN
     RETURN
  END IF failed_rec_check

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
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:) :: NORMAL_VECTORS_SUBSET
    ! Handle recursion
    IF (DIM .GT. 0) THEN
       ! Get the left and right sub-trees (left is 0 for current bit,
       ! right is 1 for current bit) and store evaluated subtrees.
       ! 
       ! Only compute the left subtree if it will not go out of bounds
       left_subtree : IF (NUM .GT. DIM) THEN
          CALL COMPUTE_NORMALS(DIM, NUM-1, LOC, &
               NORMAL_VECTORS_SUBSET(:,:1+2**(NUM-1)))
          IF (ERROR .NE. 0) RETURN
       END IF left_subtree
       ! Always compue the right subtree (because bounds are size-constrained)
       CALL COMPUTE_NORMALS(DIM-1, NUM-1, LOC+IDENTITY(NUM,:), &
            NORMAL_VECTORS_SUBSET(:,1+2**(NUM-1):))
    ELSE
       ! Set the column vector to be a vector orthogonal to all
       ! selected vectors of the box-spline direction vector set.
       CALL COMPUTE_ORTHOGONAL(TRANSPOSE(DVECS(:,FIND(REAL(LOC,REAL64)))), &
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
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:)       :: SUB_DVECS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:)       :: SHIFTED_EVAL_PTS
    REAL(KIND=REAL64), INTENT(OUT), &
         DIMENSION(SIZE(EVAL_PTS,1)) :: EVALS_AT_PTS
    ! Temporary variables for recursive computations of eval points,
    ! shifted evaluation points, and the actual evaluations.
    REAL(KIND=REAL64), DIMENSION(SIZE(SHIFTED_EVAL_PTS,1),&
         SIZE(SHIFTED_EVAL_PTS,2)) :: TEMP_SHIFTED_EVAL_PTS
    REAL(KIND=REAL64), DIMENSION(SIZE(EVAL_PTS,1),&
         SIZE(EVAL_PTS,2)) :: TEMP_EVAL_PTS
    REAL(KIND=REAL64), DIMENSION(SIZE(EVALS_AT_PTS)) :: TEMP_EVALS_AT_PTS
    ! For computing evaluations in the base case
    REAL(KIND=REAL64), DIMENSION(SIZE(EVAL_PTS,1)) :: LOCATIONS
    ! For computing recursive call values
    INTEGER,           DIMENSION(NUM_DVECS) :: NEXT_MULTS, NEXT_LOC
    REAL(KIND=REAL64), DIMENSION(NUM_DVECS) :: PT_SHIFT
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: NEXT_DVECS
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: REMAINING_DVECS
    ! Constants (integers are for indexing)
    REAL(KIND=REAL64) :: POSITION
    INTEGER :: IDX_1, IDX_2, IDX_3

    ! Recursion case ...
    IF (SUM(MULTS) > DIM) THEN
       EVALS_AT_PTS = 0.
       IDX_2 = 1
       ! Sum over all direction vectors.
       DO IDX_1 = 1, NUM_DVECS
          ! Update multiplicity of directions and position in recursion tree
          NEXT_MULTS = MULTS - IDENTITY(:, IDX_1) ! Reduce multiplicity
          NEXT_LOC   = LOC   + IDENTITY(:, IDX_1) ! Track recursion position
          ! Recursive calls
          IF (MULTS(IDX_1) .GT. 1) THEN
             ! Recursion with only reduced multiplicity
             CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, SUB_DVECS, &
                  SHIFTED_EVAL_PTS, TEMP_EVALS_AT_PTS)
             IF (ERROR .NE. 0) RETURN
             EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                  SHIFTED_EVAL_PTS(:,IDX_2)
             ! Recursion with different set of direction vectors
             compute_shift_1 : DO IDX_3 = 1, SIZE(SUB_DVECS,2)
                PT_SHIFT(IDX_3) = SUM(DVECS(:,IDX_1) * SUB_DVECS(:,IDX_3))
             END DO compute_shift_1
             shift_pts_1 : DO IDX_3 = 1, SIZE(EVAL_PTS,1)
                TEMP_SHIFTED_EVAL_PTS(IDX_3,:) = SHIFTED_EVAL_PTS(IDX_3,:) - &
                     PT_SHIFT(:SIZE(SUB_DVECS,2))
             END DO shift_pts_1
             CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, NEXT_LOC, SUB_DVECS, &
                  TEMP_SHIFTED_EVAL_PTS, TEMP_EVALS_AT_PTS)
             IF (ERROR .NE. 0) RETURN
             EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                  (MULTS(IDX_1) - SHIFTED_EVAL_PTS(:,IDX_2))
             IDX_2 = IDX_2 + 1
          ELSE IF (MULTS(IDX_1) .GT. 0) THEN
             ! Find the next direction vectors (ones with nonzero multiplicities)
             NEXT_DVECS = DVECS(:,FIND(REAL(NEXT_MULTS,REAL64)))
             IF (ERROR .NE. 0) RETURN
             IF (MATRIX_RANK(TRANSPOSE(NEXT_DVECS)) .EQ. DIM) THEN
                IF (ERROR .NE. 0) RETURN
                ! Update Least norm representation
                NEXT_DVECS = MINIMUM_NORM_REPR(NEXT_DVECS)
                IF (ERROR .NE. 0) RETURN
                ! Recursion with only reduced multiplicity
                compute_shift_2 : DO IDX_3 = 1, DIM
                   PT_SHIFT(IDX_3) = SUM(LOC * DVECS(IDX_3,:))
                END DO compute_shift_2
                shift_pts_2 : DO IDX_3 = 1, SIZE(EVAL_PTS,1)
                   TEMP_EVAL_PTS(IDX_3,:) = EVAL_PTS(IDX_3,:) - PT_SHIFT(:DIM)
                END DO shift_pts_2
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     SHIFTED_EVAL_PTS(:,IDX_2)
                ! Recursion with different set of direction vectors
                compute_shift_3 : DO IDX_3 = 1, DIM
                   PT_SHIFT(IDX_3) = SUM(NEXT_LOC * DVECS(IDX_3,:))
                END DO compute_shift_3
                shift_pts_3 : DO IDX_3 = 1, SIZE(EVAL_PTS,1)
                   TEMP_EVAL_PTS(IDX_3,:) = EVAL_PTS(IDX_3,:) - PT_SHIFT(:DIM)
                END DO shift_pts_3
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, NEXT_LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     (MULTS(IDX_1) - SHIFTED_EVAL_PTS(:,IDX_2))
             END IF
             IDX_2 = IDX_2 + 1
          END IF
       END DO
       ! Normalization by number of direction vectors.
       EVALS_AT_PTS = EVALS_AT_PTS / (SUM(MULTS) - DIM)
    ELSE
       ! Base case ... compute characteristic function
       EVALS_AT_PTS = 1.
       ! Delayed translations (this is what makes the algorithm more stable)
       REMAINING_DVECS = FIND(REAL(MULTS,REAL64))
       IF (ERROR .NE. 0) RETURN
       compute_shift_4 : DO IDX_1 = 1, DIM
          PT_SHIFT(IDX_1) = SUM(LOC * DVECS(IDX_1,:))
       END DO compute_shift_4
       shift_pts_4 : DO IDX_1 = 1, SIZE(EVAL_PTS,1)
          TEMP_EVAL_PTS(IDX_1,:) = EVAL_PTS(IDX_1,:) - PT_SHIFT(:DIM)
       END DO shift_pts_4
       ! Check against all hyperplanes
       DO IDX_1 = 1, DIM
          ! Lookup normal vector to current hyperplane
          IDX_3 = 1 + SUM(LOOKUP * ( &
               MULTS - IDENTITY(:,REMAINING_DVECS(IDX_1)) ))
          ! Compute shifted position (relative to normal ector)
          POSITION = SUM(DVECS(:,REMAINING_DVECS(IDX_1)) * NORMAL_VECTORS(:,IDX_3))
          ! Compute shifted evaluation locations = (NUM_PTS, DIM) x (DIM, 1)
          compute_shifted_point_1 : DO IDX_2 = 1, SIZE(EVAL_PTS,1)
             LOCATIONS(IDX_2) = SUM(TEMP_EVAL_PTS(IDX_2,:) * NORMAL_VECTORS(:,IDX_3))
          END DO compute_shifted_point_1
          ! Identify those points that are outside of this box
          IF (POSITION .GT. 0) THEN
             WHERE (LOCATIONS .LT. 0) EVALS_AT_PTS = 0.
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (LOCATIONS .GE. 0) EVALS_AT_PTS = 0.
          END IF
          ! Recompute shifted locations
          NEXT_LOC = LOC + IDENTITY(:,REMAINING_DVECS(IDX_1))
          compute_shift_5 : DO IDX_2 = 1, DIM
             PT_SHIFT(IDX_2) = SUM(NEXT_LOC * DVECS(IDX_2,:))
          END DO compute_shift_5
          compute_shifted_point_2 : DO IDX_2 = 1, SIZE(EVAL_PTS,1)
             LOCATIONS(IDX_2) = SUM((EVAL_PTS(IDX_2,:) - PT_SHIFT(:DIM)) * &
                  NORMAL_VECTORS(:,IDX_3))
          END DO compute_shifted_point_2
          ! Identify those points that are outside of this box
          IF (POSITION .GT. 0) THEN
             WHERE (LOCATIONS .GE. 0) EVALS_AT_PTS = 0.
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (LOCATIONS .LT. 0) EVALS_AT_PTS = 0.
          END IF
       END DO
       ! Normalization of evaluations by determinant of box.
       EVALS_AT_PTS = EVALS_AT_PTS / ABS( &
            DET(TRANSPOSE(DVECS(:,REMAINING_DVECS(1:DIM)))) )
       IF (ERROR .NE. 0) RETURN
    END IF
  END SUBROUTINE EVALUATE_BOX_SPLINE

  !===============================================================
  !     Mathematical Operations (Following Matlab Intrinsics)     
  !===============================================================

  ! ==================================================================
  SUBROUTINE COMPUTE_ORTHOGONAL(A, ORTHOGONAL)
    ! 3) COMPUTE_ORTHOGONAL
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
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)       :: A
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(A,2)) :: ORTHOGONAL
    ! Local variables for computing the orthogonal vector
    REAL(KIND=REAL64), DIMENSION(MIN(SIZE(A,1),SIZE(A,2))) :: S
    REAL(KIND=REAL64), DIMENSION(SIZE(A,2),SIZE(A,2)) :: VT
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: IDX
    LOGICAL :: FOUND_ZERO
    ! Unused parameters
    REAL(KIND=REAL64), DIMENSION(1) :: U

    ! Query the size of the work array to construct
    CALL DGESVD('N', 'A', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), S, &
         U, SIZE(U,1), VT, SIZE(VT,1), U, -1, INFO)

    IF (INFO .NE. 0) THEN
       ERROR = 34
       RETURN
    END IF

    ! Allocate the work array
    ALLOCATE(WORK(INT(U(1))))

    ! Use the SVD to get the orthogonal vectors
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

    IF (INFO .NE. 0) THEN
       ERROR = 35
       RETURN
    END IF

    ORTHOGONAL = 0.
    FOUND_ZERO = .FALSE.
    ! Find the first vector in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value)
    find_null : DO IDX = 1, SIZE(S)
       IF (S(IDX) .LE. (EPSILON(U(1)) * SIZE(S))) THEN
          ! If we found a singular value, copy out the orthogonal vector.
          ORTHOGONAL = VT(IDX,:)
          FOUND_ZERO = .TRUE.
          EXIT find_null
       END IF
    END DO find_null
    ! If no orthogonal was found and the matrix is not a span of the
    ! space, use the first vector in VT at (RANK(A) + 1).
    IF ((SIZE(VT,1) > SIZE(S)) .AND. (.NOT. FOUND_ZERO))THEN
       ORTHOGONAL = VT(SIZE(S)+1,:)
    END IF
  END SUBROUTINE COMPUTE_ORTHOGONAL

  ! ==================================================================
  FUNCTION MINIMUM_NORM_REPR(MATRIX) RESULT(MIN_NORM)
    ! 4) MINIMUM_NORM_REPR
    ! 
    !   Compute the minimum norm representation of 'MARTIX' and store
    !   it in 'MIN_NORM', use DGELS for least squares computation.
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   MIN_NORM -- Real dense matrix that is the minimum norm
    !               representation of MATRIX computed by finding the
    !               least squares solution to the problem (AA^T)X = A.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: MIN_NORM
    ! Local variables
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: SQUARE
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: DGELS_WORK_ARRAY
    REAL(KIND=REAL64), DIMENSION(1) :: DGELS_SIZE_HOLDER
    INTEGER :: DGELS_DIM
    ! Allocate the output matrix
    ALLOCATE(MIN_NORM(1:SIZE(MATRIX,1),1:SIZE(MATRIX,2)))
    ! Store the transpose and square of the matrix
    SQUARE = MATMUL(MATRIX, TRANSPOSE(MATRIX))
    MIN_NORM = MATRIX
    ! Get the size of the work array necessary
    CALL DGELS('N', SIZE(SQUARE,1), SIZE(SQUARE,2), SIZE(MIN_NORM,2),&
         SQUARE, SIZE(SQUARE,1), MIN_NORM, SIZE(MIN_NORM,1), &
         DGELS_SIZE_HOLDER, -1, INFO)

    IF (INFO .NE. 0) THEN
       ERROR = 42
       RETURN
    END IF

    DGELS_DIM = DGELS_SIZE_HOLDER(1)
    ALLOCATE( DGELS_WORK_ARRAY(1:DGELS_DIM) )

    ! Call DGELS for actual solve
    CALL DGELS('N', SIZE(SQUARE,1), SIZE(SQUARE,2), SIZE(MIN_NORM,2),&
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

    IF (INFO .NE. 0) THEN
       ERROR = 43
       RETURN
    END IF
  END FUNCTION MINIMUM_NORM_REPR

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
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    ! Local variables for computing the orthogonal vector
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: A
    REAL(KIND=REAL64), DIMENSION(MIN(SIZE(MATRIX,1),SIZE(MATRIX,2))) :: S
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: IDX, MATRIX_RANK
    ! Unused parameters
    REAL(KIND=REAL64), DIMENSION(1) :: U, VT

    A = MATRIX
    ! Query the size of the work array to construct
    CALL DGESVD('N', 'N', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), S, &
         U, SIZE(U), VT, SIZE(VT), U, -1, INFO)

    IF (INFO .NE. 0) THEN
       ERROR = 54
       RETURN
    END IF

    ! Allocate the work array
    ALLOCATE(WORK(INT(U(1))))

    ! Use the SVD to get the orthogonal vectors
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

    IF (INFO .NE. 0) THEN
       ERROR = 55
       RETURN
    END IF

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
  FUNCTION DET(MATRIX)
    ! 6) DET
    ! 
    !   Compute the determinant of a matrix without modifying it
    !   using the QR decomposition (and a copy).
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   DET -- The real-valued determinant of matrix M.
    ! 
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=REAL64) :: DET
    ! Local variable
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: M
    INTEGER, DIMENSION(MIN(SIZE(MATRIX,1), SIZE(MATRIX,2))) :: IPIV
    INTEGER :: IDX
    ! Copy into a local matrix
    M = MATRIX
    ! Do the LU factorization
    CALL DGETRF(SIZE(M,1), SIZE(M,2), M, SIZE(M,1), IPIV, INFO)

    IF (INFO .NE. 0) THEN
       ERROR = 66
    END IF

    ! Compute the determinant (product of diagonal of U)
    DET = 1.
    DO IDX = 1, MIN(SIZE(M,1), SIZE(M,2))
       DET = DET * M(IDX,IDX)
    END DO
  END FUNCTION DET

  ! ==================================================================
  FUNCTION FIND(ARRAY) RESULT(NE_ZERO)
    ! 7) FIND
    ! 
    ! Return a new array of the indices of 'ARRAY' that contain
    ! nonzero elements. Print error and return array with 0 of ARRAY=0.
    ! 
    ! Input:
    !   ARRAY -- Real array of numbers.
    ! 
    ! Output:
    !   NE_ZERO -- Integer array corresponding to those indices of
    !              ARRAY that contain nonzero elements.
    ! 
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: ARRAY
    INTEGER, DIMENSION(:), ALLOCATABLE :: NE_ZERO
    ! Local variables
    INTEGER, DIMENSION(SIZE(ARRAY)) :: INDICES
    INTEGER :: COUNT_NONZERO, IDX
    ! Identify nonzero elements of ARRAY
    COUNT_NONZERO = 0
    DO IDX = 1, SIZE(ARRAY)
       ! IF (ABS(ARRAY(IDX)) .GT. EPSILON(ARRAY(IDX))) THEN
       IF (ARRAY(IDX) .NE. 0) THEN
          COUNT_NONZERO = COUNT_NONZERO + 1
          INDICES(COUNT_NONZERO) = IDX
       END IF
    END DO
    usage_check : IF (COUNT_NONZERO .LE. 0) THEN
       ERROR = 70
       ALLOCATE(NE_ZERO(1))
       NE_ZERO(1) = 0
    ELSE
       ! Allocate the smaller output array and copy in the values
       ALLOCATE(NE_ZERO(1:COUNT_NONZERO))
       NE_ZERO = INDICES(:COUNT_NONZERO)
    END IF usage_check
  END FUNCTION FIND

END SUBROUTINE BOXSPLEV