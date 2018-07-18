! ====================================================================
SUBROUTINE BOX_EVAL(X, NU, P, B, ERROR)
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
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)       :: X
  ! Multiplicities of vectors
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(SIZE(X,2)) :: NU
  ! Locations of row-vector points at which the box spline is evaluated
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)       :: P
  ! Output box evaluations
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(P,1)) :: B
  ! Output error flag
  !  0 - Successful execution
  !  1 - Improper usage, one of the elements of NU was <= 0, only NU > 0 allowed.
  !  10 - Call to DGELS querying work array size failed.
  !  11 - Call to DGESVD querying work array size failed.
  !  20 - Computation of DGELS failed.
  !  21 - Computation of DGESVD failed.
  INTEGER,           INTENT(OUT)                       :: ERROR

  ! Allocate 'global' variables for the recursive computations

  ! Local copies of provided P and X
  REAL(KIND=REAL64), DIMENSION(SIZE(X,2),SIZE(X,1)) :: BOXEV_X
  REAL(KIND=REAL64), DIMENSION(SIZE(P,1),SIZE(P,2)) :: BOXEV_P
  ! Identity matrix
  REAL(KIND=REAL64), DIMENSION(SIZE(X,2),SIZE(X,2)) :: BOXEV_I
  ! matrix - BOXEV_J * row_vector
  REAL(KIND=REAL64), DIMENSION(SIZE(P,1),1)         :: BOXEV_J
  ! Table of normal vectors.
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE    :: BOXEV_N
  ! Hashing function for table of normal vectors.
  REAL(KIND=REAL64), DIMENSION(1,SIZE(X,2))         :: BOXEV_U
  ! Local copy of NU
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE    :: Y
  ! Number of direction vectors, dimension of domain space,
  !  number of box spline evaluation points, and an iterator respectively.
  INTEGER :: BOXEV_K, BOXEV_S, N, IDX
  REAL(KIND=REAL64), DIMENSION(SIZE(X,2)) :: BOX_WORK_ARRAY

  ! Initialize error flag
  ERROR = 0
  ! Store locals
  BOXEV_X = TRANSPOSE(X)
  BOXEV_P = P
  ! Store constants
  BOXEV_K = SIZE(X,2)
  BOXEV_S = SIZE(X,1)
  N = SIZE(P,1)
  ! Set BOXEV_J to be the 1-vector
  BOXEV_J = 1
  ! Make BOXEV_I into the identity matrix
  BOXEV_I = 0
  set_identity : DO IDX = 1, SIZE(BOXEV_I,1)
     BOXEV_I(IDX,IDX) = 1
  END DO set_identity
  ! Table lookup function (binary digits to select direction vectors)
  calculate_u : DO IDX = 0, BOXEV_K-1
     BOXEV_U(1,IDX+1) = 2 ** IDX
  END DO calculate_u

  ! Error checking
  bad_nu_check : IF (MINVAL(NU) .LT. 1) THEN
     ERROR = 1
     RETURN
  END IF bad_nu_check

  ! Get the minimum norm representation of X and keep it in Y
  Y = MINIMUM_NORM_REPR(X, ERROR)

  ! gfortran box_eval.f90 -llapack -lblas
  failed_dgels_check : IF (ERROR .NE. 0) THEN
     RETURN
  END IF failed_dgels_check

  ! Compute the table of normal vectors defining boundaries of
  ! polynomial pieces of box spline.
  BOX_WORK_ARRAY = 0
  BOXEV_N = BOX_NORM(BOXEV_S-1, BOXEV_K, BOX_WORK_ARRAY, ERROR)

  ERROR = 0
  ! Recursive evaluation of box spline
  BOX_WORK_ARRAY = 0
  B = BOX_REC(NU, BOX_WORK_ARRAY, Y, MATMUL(P,Y), ERROR)

CONTAINS

  ! ==================================================================
  RECURSIVE FUNCTION BOX_NORM(T, K, M, ERROR) RESULT(N)
    ! Evaluate and store the normal vectors that define polynomial
    ! piece boundaries for the box spline with given direction vector
    ! set.

    ! Inputs and output
    INTEGER,           INTENT(IN)                  :: T, K
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:)   :: M
    INTEGER,           INTENT(OUT)                 :: ERROR
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: N
    ! Allocate space for N
    ALLOCATE(N(1:BOXEV_S, 1:2**K))
    ! Initialize N to be zeros
    N = 0
    IF (K .GE. T) THEN
       IF (T .GT. 0) THEN
          ! Get the left and right sub-trees (left is 0 for current bit,
          ! right is 1 for current bit) and store evaluated subtrees.
          N(:,:2**(K-1)+1) = BOX_NORM(T, K-1, M, ERROR)
          N(:,1+2**(K-1):) = BOX_NORM(T-1, K-1, M+BOXEV_I(K,:), ERROR)
       ELSE
          ! Set the column vector to be the vector orthogonal to all
          ! selected rows of the box spline direction vector set
          N(:,1) = ORTHOGONAL(BOXEV_X(INT(FIND(M)),:), ERROR)
       END IF
    END IF
  END FUNCTION BOX_NORM

  ! ==================================================================
  RECURSIVE FUNCTION BOX_REC(N, M, Y, T, ERROR) RESULT(B)
    IMPLICIT NONE
    ! Multiplicities of direction vectors
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(BOXEV_K,1) :: N
    ! Box work array
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(1,BOXEV_K) :: M
    ! Direction vectors transposed
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:)       :: Y
    ! Evaluation points
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:)       :: T
    ! Error code
    INTEGER,           INTENT(OUT)                      :: ERROR
    ! Evaluations of box spline
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:)        :: B
    ! Holder for checking outside/inside conditions
    REAL(KIND=REAL64), DIMENSION(SIZE(T,1))             :: SUBTRACTOR
    ! Local variables
    INTEGER :: I, J
    REAL(KIND=REAL64), DIMENSION(BOXEV_K,1) :: NN
    REAL(KIND=REAL64), DIMENSION(1,BOXEV_K) :: MM
    REAL(KIND=REAL64), DIMENSION(:),   ALLOCATABLE :: V
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: Z, TEMP_NN, P, Q

    ! Allocate B (for return)
    ALLOCATE(B(1:SIZE(T,1)))

    ! Recursion case ...
    IF (SUM(N(:,1)) > BOXEV_S) THEN
       B = 0.
       J = 1
       ! Sum over the remaining directions in BoxEv_X ...
       DO I = 1, BOXEV_K
          ! Update multiplicity of directions and position in recursion tree
          NN = N - BOXEV_I(:,I:I) ! Reduce multiplicity
          MM = M + BOXEV_I(I:I,:) ! Track recursion position
          ! Recursive calls
          IF (N(I,1) .GT. 1) THEN
             B = B + T(:,J)*BOX_REC(NN, M, Y, T, ERROR)
             B = B + (N(I,1)-T(:,J)) * BOX_REC(NN, MM, Y, &
                  T - MATMUL(BOXEV_J,MATMUL(BOXEV_X(I:I,:),Y)), ERROR)
             J = J + 1
          ELSE IF (N(I,1) .GT. 0) THEN
             ! Update Least norm representation
             Z = BOXEV_X(INT(FIND(NN(:,1))),:)
             IF (MATRIX_RANK(Z,ERROR) .EQ. BOXEV_S) THEN
                Z = MINIMUM_NORM_REPR(TRANSPOSE(Z), ERROR)
                B = B +          T(:,J) * BOX_REC(NN,M ,Z,MATMUL((BOXEV_P-&
                     MATMUL(BOXEV_J,MATMUL(M, BOXEV_X))),Z), ERROR)
                B = B + (N(I,1)-T(:,J)) * BOX_REC(NN,MM,Z,MATMUL((BOXEV_P-&
                     MATMUL(BOXEV_J,MATMUL(MM,BOXEV_X))),Z), ERROR)
             END IF
             J = J+1
          END IF
       END DO
       ! Normalization
       B = B / REAL(SUM(N(:,1))-BOXEV_S, REAL64)
    ELSE
       ! Base case ... compute characteristic function
       B = 1.
       ! Delayed translations
       V = FIND(N(:,1))
       Z = BOXEV_P - MATMUL(BOXEV_J, MATMUL(M, BOXEV_X))
       ! Check against all hyperplanes
       DO I = 1, BOXEV_S
          ! Lookup normal vector to current hyperplane
          TEMP_NN = 1 + MATMUL(BOXEV_U, (N-BOXEV_I(:,INT(V(I)):INT(V(I)))))
          TEMP_NN = BOXEV_N(:,INT(TEMP_NN(1,1)):INT(TEMP_NN(1,1)))
          ! Box is half-open!!!
          P = MATMUL(BOXEV_X(INT(V(I)):INT(V(I)),:), TEMP_NN)
          Q  = MATMUL(Z, TEMP_NN)
          ! Change some elements of subtractor to 1 based on location 
          ! relative to normal vector
          SUBTRACTOR = 0
          IF (P(1,1) .GT. 0) THEN
             WHERE (Q(:,1) .LT. 0) SUBTRACTOR = MAX(SUBTRACTOR,1.)
          END IF
          IF (P(1,1) .LT. 0) THEN
             WHERE (Q(:,1) .GE. 0) SUBTRACTOR = MAX(SUBTRACTOR,1.)
          END IF
          B = MIN(B,REAL(1. - SUBTRACTOR,REAL64))
          ! Change a few more elements based on location.
          Q = MATMUL(BOXEV_P - MATMUL(BOXEV_J, MATMUL(M + &
               BOXEV_I(INT(V(I)):INT(V(I)),:), BOXEV_X)), TEMP_NN)
          SUBTRACTOR = 0
          IF (P(1,1) .GT. 0) THEN
             WHERE (Q(:,1) .GE. 0) SUBTRACTOR = MAX(SUBTRACTOR,1.)
          END IF
          IF (P(1,1) .LT. 0) THEN
             WHERE (Q(:,1) .LT. 0) SUBTRACTOR = MAX(SUBTRACTOR,1.)
          END IF
          B = MIN(B,REAL(1. - SUBTRACTOR,REAL64))
       END DO
       ! Normalization
       B = B / ABS( DET(BOXEV_X(INT(V(1:BOXEV_S)),:)) )
    END IF
  END FUNCTION BOX_REC

  !===============================================================
  !     Mathematical Operations (Following Matlab Intrinsics)     
  !===============================================================

  ! ==================================================================
  FUNCTION MINIMUM_NORM_REPR(MATRIX, ERROR) RESULT(MIN_NORM)
    ! Compute the new minimum norm representation of 'MARTIX' and
    ! store it in 'MIN_NORM', report any errors in 'ERROR'.
    ! 
    !  Error Codes:
    !    0   -- Successful execution
    !    2   -- Allocation error
    !    >10 -- DGELS error 'ERROR - 10'
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: MIN_NORM
    INTEGER,           INTENT(OUT)                 :: ERROR
    ! Local variables
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: SQUARE
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: DGELS_WORK_ARRAY
    REAL(KIND=REAL64), DIMENSION(1) :: DGELS_SIZE_HOLDER
    INTEGER :: DGELS_DIM
    ! Allocate the output matrix
    ALLOCATE(MIN_NORM(1:SIZE(MATRIX,1),1:SIZE(MATRIX,2)))
    ! Initialize error flag
    ERROR = 0
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
       ERROR = 2
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
  END FUNCTION MINIMUM_NORM_REPR

  ! ==================================================================
  FUNCTION ORTHOGONAL(MATRIX, ERROR) RESULT(ORTHO)
    ! Get a vector orthogonal to the row vectors stored in matrix
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    INTEGER,           INTENT(OUT)                 :: ERROR
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE   :: ORTHO 
    ! Local variables for computing the orthogonal vector
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: A
    REAL(KIND=REAL64), DIMENSION(MIN(SIZE(MATRIX,1),SIZE(MATRIX,1))) :: S
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,2),SIZE(MATRIX,2)) :: VT
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: IDX
    ! Unused parameters
    REAL(KIND=REAL64), DIMENSION(1) :: U

    A = MATRIX
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
    !   A          -- Inputs matrix A
    !   SIZE(A,1)  -- Leading dimension of A
    !   S          -- Container for singular values of A
    !   U          -- Double precision array for storing U
    !   SIZE(U,1)  -- Leading dimension of U
    !   VT         -- Double precision array for storing V^T
    !   SIZE(VT,1) -- Leading dimension of V^T
    !   WORK       -- Work array for computing SVD
    !   SIZE(WORK) -- Size of the work array
    !   INFO       -- Info message parameter

    ! Find the first vector in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value)
    find_null : DO IDX = 1, SIZE(S)
       IF (S(IDX) .LE. (EPSILON(U(1)) * MIN(SIZE(A,1),SIZE(A,2)))) THEN
          EXIT find_null
       END IF
    END DO find_null

    ! Allocate space for the output vector
    ALLOCATE(ORTHO(1:SIZE(MATRIX,2)))
    ! Copy out the vector
    ORTHO = VT(IDX,:)
  END FUNCTION ORTHOGONAL

  ! ==================================================================
  FUNCTION MATRIX_RANK(MATRIX, ERROR)
    ! Get a vector orthogonal to the row vectors stored in matrix
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    INTEGER,           INTENT(OUT)                 :: ERROR
    ! Local variables for computing the orthogonal vector
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: A
    REAL(KIND=REAL64), DIMENSION(MIN(SIZE(MATRIX,1),SIZE(MATRIX,1))) :: S
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
    ! Return a new array with only the nonzero elements of ARRAY.
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: ARRAY
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: NE_ZERO
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
       ALLOCATE(NE_ZERO(1))
       NE_ZERO(1) = 0
    ELSE
       ! Allocate the smaller output array and copy in the values
       ALLOCATE(NE_ZERO(1:COUNT_NONZERO))
       NE_ZERO = INDICES(:COUNT_NONZERO)
    END IF usage_check
  END FUNCTION FIND


END SUBROUTINE BOX_EVAL
