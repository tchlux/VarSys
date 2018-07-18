! ====================================================================
SUBROUTINE BOX_EVAL(DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR, INFO)
  ! Evaluate a box spline given direction vector set and
  ! multiplicity of direction vectors (output stored in B)
  ! 
  ! This subroutine utilizes the following LAPACK routines:
  !     DGELS   -- Computing a minimum norm representation with LS.
  !     DGESVD  -- Computing SVD to find orthogonal vector.
  !     DGETRF  -- Computing the determinant of matrix with QR.
  ! 
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  ! Direction column vector matrix (k x s)
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)              :: DVECS
  INTEGER,           INTENT(IN),  DIMENSION(SIZE(DVECS,2))    :: DVEC_MULTS
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)              :: EVAL_PTS
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(EVAL_PTS,1)) :: BOX_EVALS
  INTEGER,           INTENT(OUT)                              :: ERROR, INFO

  ! Local variables
  INTEGER :: DIM, NUM_DVECS, IDX
  CHARACTER(LEN=5) :: CHARS
  INTEGER, DIMENSION(SIZE(DVECS,2))                            :: LOCATION
  REAL(KIND=REAL64), DIMENSION(SIZE(DVECS,2))                  :: LOOKUP
  INTEGER,           DIMENSION(SIZE(DVECS,2),SIZE(DVECS,2))    :: IDENTITY
  REAL(KIND=REAL64), DIMENSION(SIZE(DVECS,1),2**SIZE(DVECS,2)) :: NORMAL_VECTORS
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE               :: MIN_NORM_DVECS

  ! Initialize error flag and info flag
  INFO = 0
  ERROR = 0
  DIM = SIZE(DVECS, 1)
  NUM_DVECS = SIZE(DVECS, 2)
  ! Create an identity matrix (for easy matrix-vector multiplication)
  ! and compute the index lookup indices (for binary style hashing)
  IDENTITY = 0
  init_ident_lookup : DO IDX = 1, NUM_DVECS
     IDENTITY(IDX,IDX) = 1
     LOOKUP(IDX) = 2**(IDX-1)
  END DO init_ident_lookup
  ! Error checking
  bad_multiplicity_check : IF (MINVAL(DVEC_MULTS) .LT. 1) THEN
     ERROR = 1
     RETURN
  END IF bad_multiplicity_check

  ! Get the minimum norm representation of the direction vectors
  MIN_NORM_DVECS = MINIMUM_NORM_REPR(DVECS)
  ! Error checking
  failed_dgels_check : IF (ERROR .NE. 0) THEN
     INFO = ERROR
     ERROR = 10
     RETURN
  END IF failed_dgels_check
  
  ! Compute the table of normal vectors defining boundaries of
  ! polynomial pieces of box spline.
  LOCATION = 0
  NORMAL_VECTORS = 0
  CALL COMPUTE_NORMALS(DIM-1, NUM_DVECS, LOCATION, NORMAL_VECTORS)
  failed_norm_check : IF (ERROR .NE. 0) THEN
     RETURN
  END IF failed_norm_check

  ! C-x C-k e -- Edit currently defined keyboard macro.
  ! C-x C-k n -- Give name to the most recently defined keyboard macro (session).
  ! C-x C-k b -- Bind most recent keyboard macro to a key (session).

  ! Recursive evaluation of box spline
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
    !   piece boundaries for the box spline with given direction
    !   vector set. (O(2^"NUM") calls in recursive subtree).
    ! 
    ! INPUTS


    ! Inputs and output
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
       END IF left_subtree
       ! Always compue the right subtree (because bounds are size-constrained)
       CALL COMPUTE_NORMALS(DIM-1, NUM-1, LOC+IDENTITY(NUM,:), &
            NORMAL_VECTORS_SUBSET(:,1+2**(NUM-1):))
    ELSE
       ! Set the column vector to be a vector orthogonal to all
       ! selected vectors of the box spline direction vector set
       CALL COMPUTE_ORTHOGONAL(TRANSPOSE(DVECS(:,FIND(REAL(LOC,REAL64)))), &
            NORMAL_VECTORS_SUBSET(:,1))
    END IF
  END SUBROUTINE COMPUTE_NORMALS
  
  ! ==================================================================
  RECURSIVE SUBROUTINE EVALUATE_BOX_SPLINE(MULTS, LOC, SUB_DVECS,&
       SHIFTED_EVAL_PTS, EVALS_AT_PTS)
    ! 2) EVALUATE_BOX_SPLINE:
    !   
    !   Evaluate the box spline defined by "DVECS" recursively, where
    !   this iteration has the remaining vectors "SUB_DVECS", at all
    !   points in "SHIFTED_EVAL_PTS" and store box spline evaluations
    !   in "EVALS_AT_PTS".
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
    REAL(KIND=REAL64), DIMENSION(SIZE(EVAL_PTS,1)) :: LOCATIONS, SUBTRACTOR
    ! For computing recursive call values
    INTEGER,           DIMENSION(NUM_DVECS) :: NEXT_MULTS, NEXT_LOC
    REAL(KIND=REAL64), DIMENSION(NUM_DVECS) :: PT_SHIFT
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: NEXT_DVECS
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: REMAINING_DVECS
    ! Constants (integers are for indexing)
    REAL(KIND=REAL64) :: POSITION
    INTEGER :: IDX, IDX_1, IDX_2, NV_IDX

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
             ! Left recursion
             CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, SUB_DVECS, &
                  SHIFTED_EVAL_PTS, TEMP_EVALS_AT_PTS)
             EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                  SHIFTED_EVAL_PTS(:,IDX_2)
             ! Right recursion
             compute_shift_1 : DO IDX = 1, SIZE(SUB_DVECS,2)
                PT_SHIFT(IDX) = SUM(DVECS(:,IDX_1) * SUB_DVECS(:,IDX))
             END DO compute_shift_1
             shift_pts_1 : DO IDX = 1, SIZE(EVAL_PTS,1)
                TEMP_SHIFTED_EVAL_PTS(IDX,:) = SHIFTED_EVAL_PTS(IDX,:) - &
                     PT_SHIFT(:SIZE(SUB_DVECS,2))
             END DO shift_pts_1
             CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, NEXT_LOC, SUB_DVECS, &
                  TEMP_SHIFTED_EVAL_PTS, TEMP_EVALS_AT_PTS)
             EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                  (MULTS(IDX_1) - SHIFTED_EVAL_PTS(:,IDX_2))
             IDX_2 = IDX_2 + 1
          ELSE IF (MULTS(IDX_1) .GT. 0) THEN
             ! Find the next direction vectors (ones with nonzero multiplicities)
             NEXT_DVECS = DVECS(:,FIND(REAL(NEXT_MULTS,REAL64)))
             IF (MATRIX_RANK(TRANSPOSE(NEXT_DVECS)) .EQ. DIM) THEN
                ! Update Least norm representation
                NEXT_DVECS = MINIMUM_NORM_REPR(NEXT_DVECS)
                ! Left recursion
                compute_shift_2 : DO IDX = 1, DIM
                   PT_SHIFT(IDX) = SUM(LOC * DVECS(IDX,:))
                END DO compute_shift_2
                shift_pts_2 : DO IDX = 1, SIZE(EVAL_PTS,1)
                   TEMP_EVAL_PTS(IDX,:) = EVAL_PTS(IDX,:) - PT_SHIFT(:DIM)
                END DO shift_pts_2
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     SHIFTED_EVAL_PTS(:,IDX_2)
                ! Right recursion
                compute_shift_3 : DO IDX = 1, DIM
                   PT_SHIFT(IDX) = SUM(NEXT_LOC * DVECS(IDX,:))
                END DO compute_shift_3
                shift_pts_3 : DO IDX = 1, SIZE(EVAL_PTS,1)
                   TEMP_EVAL_PTS(IDX,:) = EVAL_PTS(IDX,:) - PT_SHIFT(:DIM)
                END DO shift_pts_3
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, NEXT_LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     (MULTS(IDX_1) - SHIFTED_EVAL_PTS(:,IDX_2))
             END IF
             IDX_2 = IDX_2 + 1
          END IF
       END DO
       ! Normalization
       EVALS_AT_PTS = EVALS_AT_PTS / (SUM(MULTS) - DIM)
    ELSE
       ! Base case ... compute characteristic function
       EVALS_AT_PTS = 1.
       ! Delayed translations (this is what makes the algorithm more stable)
       REMAINING_DVECS = FIND(REAL(MULTS,REAL64))
       compute_shift_4 : DO IDX = 1, DIM
          PT_SHIFT(IDX) = SUM(LOC * DVECS(IDX,:))
       END DO compute_shift_4
       shift_pts_4 : DO IDX = 1, SIZE(EVAL_PTS,1)
          TEMP_EVAL_PTS(IDX,:) = EVAL_PTS(IDX,:) - PT_SHIFT(:DIM)
       END DO shift_pts_4
       ! Check against all hyperplanes
       DO IDX_1 = 1, DIM
          ! Lookup normal vector to current hyperplane
          NV_IDX = 1 + SUM(LOOKUP * ( &
               MULTS - IDENTITY(:,REMAINING_DVECS(IDX_1)) ))
          ! Compute position
          POSITION = SUM(DVECS(:,REMAINING_DVECS(IDX_1)) * NORMAL_VECTORS(:,NV_IDX))
          ! Compute shifted locations = (NUM_PTS, DIM) x (DIM, 1)
          compute_shifted_point_1 : DO IDX = 1, SIZE(EVAL_PTS,1)
             LOCATIONS(IDX) = SUM(TEMP_EVAL_PTS(IDX,:) * NORMAL_VECTORS(:,NV_IDX))
          END DO compute_shifted_point_1
          ! Identify those points that are outside of this box
          IF (POSITION .GT. 0) THEN
             WHERE (LOCATIONS .LT. 0) EVALS_AT_PTS = 0.
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (LOCATIONS .GE. 0) EVALS_AT_PTS = 0.
          END IF
          ! Recompute shifted locations
          NEXT_LOC = LOC + IDENTITY(:,REMAINING_DVECS(IDX_1))
          compute_shift_5 : DO IDX = 1, DIM
             PT_SHIFT(IDX) = SUM(NEXT_LOC * DVECS(IDX,:))
          END DO compute_shift_5
          compute_shifted_point_2 : DO IDX = 1, SIZE(EVAL_PTS,1)
             LOCATIONS(IDX) = SUM((EVAL_PTS(IDX,:) - PT_SHIFT(:DIM)) * &
                  NORMAL_VECTORS(:,NV_IDX))
          END DO compute_shifted_point_2
          ! Identify those points that are outside of this box
          IF (POSITION .GT. 0) THEN
             WHERE (LOCATIONS .GE. 0) EVALS_AT_PTS = 0.
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (LOCATIONS .LT. 0) EVALS_AT_PTS = 0.
          END IF
       END DO
       ! Normalization of evaluations
       EVALS_AT_PTS = EVALS_AT_PTS / ABS( &
            DET(TRANSPOSE(DVECS(:,REMAINING_DVECS(1:DIM)))) )
    END IF
  END SUBROUTINE EVALUATE_BOX_SPLINE

  !===============================================================
  !     Mathematical Operations (Following Matlab Intrinsics)     
  !===============================================================

  ! ==================================================================
  SUBROUTINE COMPUTE_ORTHOGONAL(A, ORTHOGONAL)
    ! 3) COMPUTE_ORTHOGONAL
    ! 
    !   Given a matrix "A" of row vectors, compute a vector orthogonal
    !   to the row vectors in "A" and store it in "ORTHOGONAL".
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
         U, SIZE(U,1), VT, SIZE(VT,1), U, -1, ERROR)

    ! Allocate the work array
    ALLOCATE(WORK(INT(U(1))))

    ! Use the SVD to get the orthogonal vectors
    CALL DGESVD('N','A',SIZE(A,1),SIZE(A,2),A,SIZE(A,1),S,U,SIZE(U), &
         VT, SIZE(VT,1), WORK, SIZE(WORK), ERROR)
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
    !   it in 'MIN_NORM'.
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
         DGELS_SIZE_HOLDER, -1, ERROR)
    allocate_dgels_array : IF (ERROR .EQ. 0) THEN
       DGELS_DIM = DGELS_SIZE_HOLDER(1)
       ALLOCATE( DGELS_WORK_ARRAY(1:DGELS_DIM) )
    ELSE
       INFO = ERROR
       ERROR = 20
       RETURN
    END IF allocate_dgels_array
    ! Call DGELS for actual solve
    CALL DGELS('N', SIZE(SQUARE,1), SIZE(SQUARE,2), SIZE(MIN_NORM,2),&
         SQUARE, SIZE(SQUARE,1), MIN_NORM, SIZE(MIN_NORM,1), &
         DGELS_WORK_ARRAY, DGELS_DIM, ERROR)
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
    !   ERROR            -- For verifying successful execution
    ls_error_check : IF (ERROR .NE. 0) THEN
       INFO = ERROR
       ERROR = 30
       RETURN
    END IF ls_error_check
  END FUNCTION MINIMUM_NORM_REPR

  ! ==================================================================
  FUNCTION MATRIX_RANK(MATRIX)
    ! Get the rank of the provided matrix.
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
         U, SIZE(U), VT, SIZE(VT), U, -1, ERROR)

    ! Allocate the work array
    ALLOCATE(WORK(INT(U(1))))

    ! Use the SVD to get the orthogonal vectors
    CALL DGESVD('N','N',SIZE(A,1),SIZE(A,2),A,SIZE(A,1),S,U,SIZE(U), &
         VT, SIZE(VT), WORK, SIZE(WORK), ERROR)
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
  FUNCTION DET(M) RESULT(D)
    ! Compute the determinant of a matrix without modifying it.
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: M
    REAL(KIND=REAL64) :: D
    ! Local variable
    REAL(KIND=REAL64), DIMENSION(SIZE(M,1),SIZE(M,2)) :: MATRIX
    INTEGER, DIMENSION(MIN(SIZE(MATRIX,1), SIZE(MATRIX,2))) :: IPIV
    INTEGER :: INFO, IDX
    ! Copy into a local matrix
    MATRIX = M
    ! Do the LU factorization
    CALL DGETRF(SIZE(MATRIX,1), SIZE(MATRIX,2), MATRIX, &
         SIZE(MATRIX,1), IPIV, INFO)
    ! Compute the determinant (product of diagonal of U)
    D = 1
    DO IDX = 1, MIN(SIZE(MATRIX,1), SIZE(MATRIX,2))
       D = D * MATRIX(IDX,IDX)
    END DO
  END FUNCTION DET

  ! ==================================================================
  FUNCTION FIND(ARRAY) RESULT(NE_ZERO)
    ! Return a new array of the indices of 'ARRAY' that contain
    ! nonzero elements. Print error and return array with 0 of ARRAY=0.
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
       PRINT *, 'ERROR: "FIND" was given 0 array.'
       ERROR = 60
       ALLOCATE(NE_ZERO(1))
       NE_ZERO(1) = 0
    ELSE
       ! Allocate the smaller output array and copy in the values
       ALLOCATE(NE_ZERO(1:COUNT_NONZERO))
       NE_ZERO = INDICES(:COUNT_NONZERO)
    END IF usage_check
  END FUNCTION FIND

END SUBROUTINE BOX_EVAL


! ====================================================================
!       Matlab output
! -------------------------
! 
! BoxEv_N:
!    0.00000   0.00000  -1.00000   0.00000   0.70711   0.00000   0.00000   0.00000  -0.70711   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
!    0.00000   1.00000   0.00000   0.00000   0.70711   0.00000   0.00000   0.00000   0.70711   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
! 
! b:      
!  0.25005
! ====================================================================
