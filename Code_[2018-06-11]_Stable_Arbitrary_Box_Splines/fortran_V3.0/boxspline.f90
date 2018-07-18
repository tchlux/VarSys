! - Allocate all necessary storage for vectors and add an integer to
!   the recursion for "DEPTH" use that to determine storage slices.
! - Remove all local variable allocation from all functions and
!   subroutines


! This file (boxspline.f90) contains the subroutine BOXSPLEV that can
! be used to evaluate a box-spline, defined by its associated
! direction vector set, at given evaluation points, as well as the
! module REAL_PRECISION defining the real precision.
! 
! ====================================================================
! MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
!   INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
! END MODULE REAL_PRECISION

SUBROUTINE BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  ! BOXSPLEV evaluates a box-spline, defined by a direction vector set
  ! in dimension S, at given evaluation points.
  ! 
  ! This implementation uses the numerically consistent algorithm for
  ! evaluating box-splines originally presented in [1]. Most notably,
  ! the evaluation of the box-spline near the boundaries of polynomial
  ! pieces does not exhibit the random behavior that is seen in the
  ! naive recursive implementation. As described in [1], the
  ! computational complexity for direction vector sets with repeated
  ! direction vectors is reduced from the naive recursive
  ! implementation complexity
  !     O( 2^{n - s} {n! / s!} ),
  ! where "n" is the number of direction vectors and "s" is the
  ! dimension, to
  !     O( (2k)^{n - s} ),
  ! where "k" is the number of unique direction vectors, "n" is the
  ! total number of direction vectors, and "s" is the dimension.
  ! 
  ! The memory complexity of the implementation provided here is
  !     O( s k + n m s + n^2 ),
  ! where "s" is the dimension, "k" is the number of unique direction
  ! vectors, "n" is the number of direction vectors, and "m" is the
  ! number of evaluation points. This memory complexity is achieved by
  ! not precomputing the 2^k normal vectors, which require O( 2^k s )
  ! space in [1].
  ! 
  ! 
  ! [1] Kobbelt, Leif. "Stable Evaluation of Box-Splines." 
  !     Numerical Algorithms 14.4 (1997): 377-382.
  ! 
  ! The subroutine BOXSPLEV utilizes the LAPACK routines DGECON, DGELS,
  ! DGESVD, and DGETRF.
  ! 
  ! On input:
  ! 
  ! DVECS(:,:) 
  !    is a real S x M array whose columns are the unique direction
  !    vectors used in defining the box-spline.  S <= M, and DVECS(1:S,1:S)
  !    must be invertible.
  !
  ! DVEC_MULTS(:) 
  !    is an integer array of length M containing the multiplicity
  !    of each corresponding direction vector. 
  !
  ! EVAL_PTS(:,:)
  !    is a real S x L array whose columns are the points at which
  !    the box-spline is to be evaluated. 
  !
  ! On output:
  ! 
  ! BOX_EVALS(:)
  !    is a real array of length L containing the box-spline
  !    values at the evaluation points.
  !
  ! ERROR 
  !    is an integer error flag of the form
  !             100*(LAPACK INFO flag) + 10*T + U.
  !    ERROR .EQ. 0 is a normal return.  For ERROR .NE. 0, the 
  !    meanings of the tens digit T and the units digit U are:
  !
  !    Tens digit T:
  !      0  Improper usage.
  !      1  Error computing box-spline.                               
  !      2  Error computing an orthogonal vector.                     
  !      3  Error computing a minimum norm representation.            
  !      4  Error computing a matrix rank.                            
  !      5  Error computing a matrix determinant.                     
  !      6  Error computing the reciprocal condition number of matrix.
  !      7  Error finding nonzero entries in an array.                
  !      8  Error preparing memory.                                   
  ! 
  !    Units digit U, for T = 0:
  !      1  Mismatched dimension, SIZE(DVEC_MULTS) .NE. SIZE(DVECS,2).
  !      2  Mismatched dimension, SIZE(EVAL_PTS,1) .NE. SIZE(DVECS,1).
  !      3  Mismatched dimension, SIZE(BOX_EVALS)  .NE. SIZE(EVAL_PTS,1).
  !      4  One of the multiplicity values provided was < 1.
  !      5  Column vectors of DVECS are not unique.
  !      6  M < S or DVECS(1:S,1:S) is near singular.
  ! 
  !    Units digit U, for T /= 0:
  !      0  Work array allocation failed.
  !      1  Work array size query failed.
  !      2  DGELS computation error, see 'INFO' for details.
  !      3  DGESVD computation error, see 'INFO' for details.
  !      4  DGETRF computation error, see 'INFO' for details.
  !      5  DGECON computation error, see 'INFO' for details.
  ! 
  ! ------------------------------------------------------------------
  ! The calling program should include the following interface to BOXSPLEV:
  ! 
  ! INTERFACE 
  !   SUBROUTINE BOXSPLEV(DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  !     USE REAL_PRECISION, ONLY: R8
  !     REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: DVECS
  !     INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
  !     REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
  !     REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
  !     INTEGER,       INTENT(OUT)                 :: ERROR
  !   END SUBROUTINE BOXSPLEV
  ! END INTERFACE
  ! ------------------------------------------------------------------
  ! 
  USE REAL_PRECISION, ONLY: R8
  IMPLICIT NONE
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: UNIQUE_DVECS
  INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
  INTEGER,       INTENT(OUT)                 :: ERROR
  ! Single-usage recursion variables and global variables.
  INTEGER :: DIM, NUM_DVECS, NUM_PTS, DEPTH, INFO, IDX, IDX_1, IDX_2
  INTEGER,       DIMENSION(SIZE(UNIQUE_DVECS,2),SIZE(UNIQUE_DVECS,2))  :: IDENTITY
  REAL(KIND=R8), DIMENSION(SIZE(UNIQUE_DVECS,2),SIZE(UNIQUE_DVECS,2))  :: ORTHO_MIN_WORK
  REAL(KIND=R8), DIMENSION(SIZE(UNIQUE_DVECS,1) - 1)                   :: ORTHO_S
  REAL(KIND=R8), PARAMETER :: SQRTEPS = SQRT(EPSILON(1.0_R8))
  REAL(KIND=R8), PARAMETER :: ONE  = 1.0_R8
  REAL(KIND=R8), PARAMETER :: ZERO = 0.0_R8
  ! Recursion variables for evaluating box-spline. Last dimension is depth.
  REAL(KIND=R8), DIMENSION(:,:,:), ALLOCATABLE :: DVECS
  INTEGER,       DIMENSION(:,:),   ALLOCATABLE :: LOC
  INTEGER,       DIMENSION(:,:),   ALLOCATABLE :: MULTS
  real(KIND=R8), DIMENSION(:,:),   ALLOCATABLE :: EVALS_AT_PTS
  REAL(KIND=R8), DIMENSION(:,:,:), ALLOCATABLE :: SHIFTED_EVAL_PTS
  ! Single-usage variables for evaluating the box-spline.
  REAL(KIND=R8) :: POSITION
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: TEMP_EVAL_PTS
  REAL(KIND=R8), DIMENSION(:),   ALLOCATABLE :: PT_SHIFT, NORMAL_VECTOR
  REAL(KIND=R8), DIMENSION(:),   ALLOCATABLE :: LAPACK_WORK
  INTEGER,       DIMENSION(:),   ALLOCATABLE :: REMAINING_DVECS, NONZERO_DVECS

  ! Initialize error and info flags.
  ERROR = 0; INFO = 0; 
  ! Initialize parameters 
  ! Store 'global' constants for box-spline evaluation.
  DIM       = SIZE(UNIQUE_DVECS, 1)
  NUM_DVECS = SIZE(UNIQUE_DVECS, 2)
  NUM_PTS   = SIZE(EVAL_PTS, 2)
  ! Check for usage errors (dimension mismatches, invalid multiplicity)
  IF (SIZE(DVEC_MULTS)   .NE. NUM_DVECS) THEN; ERROR = 1; RETURN; END IF
  IF (SIZE(EVAL_PTS,1)   .NE. DIM)       THEN; ERROR = 2; RETURN; END IF
  IF (SIZE(BOX_EVALS)    .NE. NUM_PTS)   THEN; ERROR = 3; RETURN; END IF
  IF (MINVAL(DVEC_MULTS) .LT. 1)         THEN; ERROR = 4; RETURN; END IF
  ! Check uniqueness of DVECS columns.
  DO IDX_1 = 1, NUM_DVECS-1
     DO IDX_2 = IDX_1+1, NUM_DVECS
        IF (SUM(ABS(UNIQUE_DVECS(:,IDX_1) - UNIQUE_DVECS(:,IDX_2))) .LT. SQRTEPS) THEN
           ERROR = 5; RETURN
        END IF
     END DO
  END DO
  ! Check if NUM_DVECS < DIM or rank DVECS(1:DIM,1:DIM) < DIM.
  IF (NUM_DVECS .LT. DIM) THEN 
     ERROR = 6; RETURN
  ELSE
     ! Compute condition number COND(DVECS(1:DIM,1:DIM)) and test for near
     ! rank deficiency: 1/COND(DVECS(1:DIM,1:DIM)) < SQRT(EPSILON(1.0_R8)).
     IF (MATRIX_CONDITION_INV(UNIQUE_DVECS(1:DIM, 1:DIM)) .LT. SQRTEPS) THEN
        ERROR = 6; RETURN
     ENDIF
  ENDIF
  ! Allocate a work array large enough for all LAPACK calls.
  CALL RESERVE_MEMORY()
  IF (ERROR .NE. 0) THEN; CALL FREE_MEMORY(); RETURN; END IF
  ! Create an identity matrix (for easy matrix-vector multiplication).
  IDENTITY = 0
  FORALL (IDX = 1:NUM_DVECS) IDENTITY(IDX,IDX) = 1
  ! Calculate the DVECS for the beginning of box-spline evaluation.
  DEPTH = 0; DVECS(:,:,1) = UNIQUE_DVECS(:,:);
  CALL MAKE_DVECS_MIN_NORM()
  IF (ERROR .NE. 0) THEN; CALL FREE_MEMORY(); RETURN; END IF
  ! Initialize the multiplicities and location.
  MULTS(:,1) = DVEC_MULTS(:); LOC(:,1) = 0;
  ! Compute the shifted evaluation points.
  CALL DGEMM('T', 'N', SIZE(DVECS,2), SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
       ONE, DVECS(:,:,1), SIZE(DVECS,1), EVAL_PTS, SIZE(EVAL_PTS,1), &
       ZERO, SHIFTED_EVAL_PTS(:,:,1), SIZE(SHIFTED_EVAL_PTS,1))
  ! Recursively evaluate the box-spline.
  CALL EVALUATE_BOX_SPLINE()
  IF (ERROR .NE. 0) THEN; CALL FREE_MEMORY(); RETURN; END IF
  BOX_EVALS(:) = EVALS_AT_PTS(1:SIZE(BOX_EVALS),1)
CONTAINS

  ! ==================================================================
  RECURSIVE SUBROUTINE EVALUATE_BOX_SPLINE()
    ! 1) EVALUATE_BOX_SPLINE:
    !   
    !   Evaluate the box-spline defined by "DVECS" recursively, where
    !   this iteration handles the remaining direction vectors, at all
    !   points in "SHIFTED_EVAL_PTS" and store box-spline evaluations
    !   in "EVALS_AT_PTS".
    ! 
    INTEGER :: IDX_1, IDX_2, IDX_3, LOCAL_COUNT
    ! Adjust the global variable "DEPTH" to use appropriate memory slices.
    DEPTH = DEPTH + 1
    ! Recursion case ...
    IF (SUM(MULTS(:,DEPTH)) > DIM) THEN
       EVALS_AT_PTS(:,DEPTH) = 0_R8
       IDX_2 = 1
       ! Sum over all direction vectors.
       DO IDX_1 = 1, SIZE(MULTS,1)
          ! Update multiplicity of directions and position in recursion tree.
          MULTS(:,DEPTH+1) = MULTS(:,DEPTH) - IDENTITY(:,IDX_1) ! Reduce multiplicity.
          ! Make recursive calls.
          IF (MULTS(IDX_1,DEPTH) .GT. 1) THEN
             ! Copy LOC, DVECS, and SHIFTED_EVAL_PTS down one level for next recursion
             LOC(:,DEPTH+1)     = LOC(:,DEPTH)
             DVECS(:,:,DEPTH+1) = DVECS(:,:,DEPTH) ! Pass the DVECS down
             SHIFTED_EVAL_PTS(:,:,DEPTH+1) = SHIFTED_EVAL_PTS(:,:,DEPTH)
             ! Perform recursion with only reduced multiplicity.
             CALL EVALUATE_BOX_SPLINE()
             IF (ERROR .NE. 0) RETURN
             EVALS_AT_PTS(:,DEPTH) = EVALS_AT_PTS(:,DEPTH) + &
                  EVALS_AT_PTS(:,DEPTH+1) * SHIFTED_EVAL_PTS(IDX_2,:,DEPTH)
             ! Perform recursion with transformed set of direction
             ! vectors and evaluation points. (store at next 'DEPTH')
             PT_SHIFT(:) = MATMUL(DVECS(:,:,DEPTH), UNIQUE_DVECS(:,IDX_1))
             compute_shift_1 : DO IDX_3 = 1, NUM_PTS
                SHIFTED_EVAL_PTS(:,IDX_3,DEPTH+1) = &
                     SHIFTED_EVAL_PTS(:,IDX_3,DEPTH) - PT_SHIFT(:)
             END DO compute_shift_1
             ! Update location for next level of recursion.
             LOC(:,DEPTH+1) = LOC(:,DEPTH) + IDENTITY(:,IDX_1) 
             CALL EVALUATE_BOX_SPLINE()
             IF (ERROR .NE. 0) RETURN
             EVALS_AT_PTS(:,DEPTH) = EVALS_AT_PTS(:,DEPTH) + &
                  EVALS_AT_PTS(:,DEPTH+1) * &
                  (MULTS(IDX_1,DEPTH) - SHIFTED_EVAL_PTS(IDX_2,:,DEPTH))
             IDX_2 = IDX_2 + 1
          ELSE IF (MULTS(IDX_1,DEPTH) .GT. 0) THEN
             ! Find the next direction vectors (ones with nonzero multiplicities).
             CALL NONZERO(MULTS(:,DEPTH+1), NONZERO_DVECS)
             IF (ERROR .NE. 0) RETURN
             ! Initialize next DVECS to only those active direction vectors.
             DVECS(:,:NUM_DVECS,DEPTH+1) = UNIQUE_DVECS(:,NONZERO_DVECS(:NUM_DVECS))
             LOCAL_COUNT = NUM_DVECS
             IF (MATRIX_RANK(TRANSPOSE(DVECS(:,:NUM_DVECS,DEPTH+1))) .EQ. DIM) THEN
                IF (ERROR .NE. 0) RETURN
                ! Assign location for next recursion level to be unchanged.
                LOC(:,DEPTH+1) = LOC(:,DEPTH)
                ! Update Least norm representation of the direction vectors.
                CALL MAKE_DVECS_MIN_NORM()
                IF (ERROR .NE. 0) RETURN
                ! Perform recursion with only reduced multiplicity.
                PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), LOC(:,DEPTH))
                compute_shift_2 : DO IDX_3 = 1, NUM_PTS
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - PT_SHIFT(:)
                END DO compute_shift_2
                ! Next dvecs x temp eval points
                CALL DGEMM('T', 'N', NUM_DVECS, SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
                     ONE, DVECS(:,:NUM_DVECS,DEPTH+1), SIZE(DVECS,1), TEMP_EVAL_PTS, SIZE(EVAL_PTS,1), &
                     ZERO, SHIFTED_EVAL_PTS(:NUM_DVECS,:,DEPTH+1), NUM_DVECS)
                ! Perform recursion with only reduced multiplicity.
                CALL EVALUATE_BOX_SPLINE()
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS(:,DEPTH) = EVALS_AT_PTS(:,DEPTH) + &
                     EVALS_AT_PTS(:,DEPTH+1) * SHIFTED_EVAL_PTS(IDX_2,:,DEPTH)
                ! Update location for next level of recursion.
                LOC(:,DEPTH+1) = LOC(:,DEPTH) + IDENTITY(:,IDX_1) 
                PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), LOC(:,DEPTH+1))
                compute_shift_3 : DO IDX_3 = 1, NUM_PTS
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - PT_SHIFT(:)
                END DO compute_shift_3
                ! Recalculate "NUM_DVECS" since it was destroyed in recursion.
                NUM_DVECS = LOCAL_COUNT
                ! Next dvecs x temp eval points
                CALL DGEMM('T', 'N', NUM_DVECS, SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
                     ONE, DVECS(:,:NUM_DVECS,DEPTH+1), SIZE(DVECS,1), TEMP_EVAL_PTS, SIZE(EVAL_PTS,1), &
                     ZERO, SHIFTED_EVAL_PTS(:NUM_DVECS,:,DEPTH+1), NUM_DVECS)
                ! Perform recursion with transformed set of direction vectors.
                CALL EVALUATE_BOX_SPLINE()
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS(:,DEPTH) = EVALS_AT_PTS(:,DEPTH) + &
                     EVALS_AT_PTS(:,DEPTH+1) * &
                     (MULTS(IDX_1,DEPTH) - SHIFTED_EVAL_PTS(IDX_2,:,DEPTH))
             END IF
             IDX_2 = IDX_2 + 1
          END IF
       END DO
       ! Normalize by number of direction vectors in computation.
       EVALS_AT_PTS(:,DEPTH) = EVALS_AT_PTS(:,DEPTH) / &
            REAL(SUM(MULTS(:,DEPTH)) - DIM, R8)
    ELSE
       ! Base case ... compute characteristic function.
       EVALS_AT_PTS(:,DEPTH) = 1_R8

       ! -------------------------------------------------------------
       ! Get the set of remaining direction vectors
       CALL NONZERO(MULTS(:,DEPTH), REMAINING_DVECS)
       ! 
       ! This should be more efficient. Pack the remaining direction
       ! vectors into an array (that has them transposed). Then when
       ! the iterations are performed below, simply cycle out single
       ! direction vectors at a time.
       ! -------------------------------------------------------------

       ! Delayed translations (this is what makes the algorithm more stable).
       IF (ERROR .NE. 0) RETURN
       PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), LOC(:,DEPTH))
       compute_shift_4 : DO IDX_1 = 1, NUM_PTS
          SHIFTED_EVAL_PTS(:DIM,IDX_1,DEPTH) = EVAL_PTS(:,IDX_1) - PT_SHIFT(:)
       END DO compute_shift_4
       ! Check evaluation point locations against all remaining direction vectors.
       DO IDX_1 = 1, DIM

          ! ----------------------------------------------------------
          ! Identify the only nonzero direction vectors (except 'current').
          CALL NONZERO(MULTS(:,DEPTH) - IDENTITY(:,REMAINING_DVECS(IDX_1)), &
               NONZERO_DVECS)
          IF (ERROR .NE. 0) RETURN
          ! Calculate normal vector to current selection of direction vectors.
          CALL MATRIX_ORTHOGONAL( &
               TRANSPOSE(UNIQUE_DVECS(:,NONZERO_DVECS(:NUM_DVECS))), NORMAL_VECTOR)
          IF (ERROR .NE. 0) RETURN
          ! 
          !  This needs to be rewritten to conserve more memory. The
          !  operation is, find a normal vector to a selection of the
          !  original direction vectors.
          ! ----------------------------------------------------------

          ! Compute shifted position (relative to normal vector).
          POSITION = DOT_PRODUCT(UNIQUE_DVECS(:,REMAINING_DVECS(IDX_1)), NORMAL_VECTOR(:))
          ! Compute shifted evaluation locations. (NUM_PTS, DIM) x (DIM, 1)
          EVALS_AT_PTS(:,DEPTH+1) = MATMUL(NORMAL_VECTOR, SHIFTED_EVAL_PTS(:DIM,:,DEPTH))
          ! Identify those points that are outside of this box (0-side).
          IF (POSITION .GT. 0) THEN
             WHERE (EVALS_AT_PTS(:,DEPTH+1) .LT. 0) EVALS_AT_PTS(:,DEPTH) = 0_R8
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (EVALS_AT_PTS(:,DEPTH+1) .GE. 0) EVALS_AT_PTS(:,DEPTH) = 0_R8
          END IF
          ! Recompute shifted location (other side of box) based on selected direction vector.
          LOC(:,DEPTH+1) = LOC(:,DEPTH) + IDENTITY(:,REMAINING_DVECS(IDX_1))
          PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), LOC(:,DEPTH+1))
          compute_shifted_point_2 : DO IDX_2 = 1, NUM_PTS
             EVALS_AT_PTS(IDX_2,DEPTH+1) = SUM((EVAL_PTS(:,IDX_2) - PT_SHIFT(:)) * &
                  NORMAL_VECTOR(:))
          END DO compute_shifted_point_2
          ! Identify those shifted points that are outside of this box on REMAINING_DVEC(IDX_1)-side.
          IF (POSITION .GT. 0) THEN
             WHERE (EVALS_AT_PTS(:,DEPTH+1) .GE. 0) EVALS_AT_PTS(:,DEPTH) = 0_R8
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (EVALS_AT_PTS(:,DEPTH+1) .LT. 0) EVALS_AT_PTS(:,DEPTH) = 0_R8
          END IF
       END DO
       ! -------------------------------------------------------------
       ! Normalize evaluations by determinant of box.
       EVALS_AT_PTS(:,DEPTH) = EVALS_AT_PTS(:,DEPTH) / ABS( &
            MATRIX_DET(TRANSPOSE(UNIQUE_DVECS(:,REMAINING_DVECS(1:DIM)))) )
       IF (ERROR .NE. 0) RETURN
       ! 
       ! This needs to be cleaned up. We want the determinant of the
       ! remaining direction vectors. If packed correctly (as
       ! suggested before), then this operation will be simple.
       ! -------------------------------------------------------------
    END IF
    DEPTH = DEPTH - 1
  END SUBROUTINE EVALUATE_BOX_SPLINE

  !===============================================================
  !             Mathematical Convenience Operations               
  !===============================================================

  ! ==================================================================
  SUBROUTINE MATRIX_ORTHOGONAL(A, ORTHOGONAL)
    ! 2) MATRIX_ORTHOGONAL
    ! 
    !   Given a matrix A of row vectors, compute a vector orthogonal to
    !   the row vectors in A and store it in ORTHOGONAL using DGESVD.
    !   If there are any singular values [value < SQRT(epsilon)],
    !   the vector associated with the largest such singular values is
    !   returned.
    ! 
    ! Input:
    !   A(:,:) -- Real dense matrix.
    ! 
    ! Output:
    !   ORTHOGONAL(:) -- Real vector orthogonal to given matrix.
    ! 
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:)       :: A
    REAL(KIND=R8), INTENT(OUT), DIMENSION(SIZE(A,2)) :: ORTHOGONAL
    ! Local variables for computing the orthogonal vector
    INTEGER :: IDX
    LOGICAL :: FOUND_ZERO
    REAL(KIND=R8) :: TOL
    ! Unused DGESVD parameter
    REAL(KIND=R8), DIMENSION(1) :: U

    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N', 'A', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), &
         ORTHO_S, U, SIZE(U), ORTHO_MIN_WORK(:SIZE(A,2),:SIZE(A,2)), &
         SIZE(A,2), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
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

    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 23; RETURN; END IF

    ORTHOGONAL(:) = 0._R8
    TOL = EPSILON(0._R8) * SIZE(ORTHO_S) ! <- Compute tolerance for considering a singular value '0'
    ! Find a vector in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value).
    IF (ORTHO_S(SIZE(ORTHO_S)) .LE. TOL) THEN
       ORTHOGONAL = ORTHO_MIN_WORK(SIZE(ORTHO_S),:SIZE(A,2))
       FOUND_ZERO = .TRUE.
    ! If no orthogonal vector was found and the matrix does not contain
    ! enough vectors to span the space, use the first vector in VT at
    ! (RANK(A) + 1).
    ELSE IF ((SIZE(A,2) > SIZE(ORTHO_S)) .AND. (.NOT. FOUND_ZERO))THEN
       ORTHOGONAL = ORTHO_MIN_WORK(SIZE(ORTHO_S)+1,:SIZE(A,2))
    END IF
  END SUBROUTINE MATRIX_ORTHOGONAL

  ! ==================================================================
  SUBROUTINE MAKE_DVECS_MIN_NORM()
    ! 3) MAKE_DVECS_MIN_NORM
    ! 
    !   Compute the minimum norm representation of 'MATRIX' and store
    !   it in 'MIN_NORM', use DGELS to find the least squares
    !   solution to the problem (AX = I). This is a more numerically
    !   stable solution to the linear system (A^T A) X = A^T.

    ! Make "ORTHO_MIN_WORK" the identity matrix
    ORTHO_MIN_WORK(:,:) = 0_R8
    FORALL (IDX = 1:NUM_DVECS) ORTHO_MIN_WORK(IDX,IDX) = 1_R8
    ! Call DGELS for actual solve.
    CALL DGELS('T', SIZE(DVECS,1), NUM_DVECS, SIZE(DVECS,2),&
         DVECS(:,:NUM_DVECS,DEPTH+1), SIZE(DVECS,1), ORTHO_MIN_WORK, &
         SIZE(ORTHO_MIN_WORK,1), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    ! Check for error.
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 32; RETURN; END IF
    ! Extract the minimum norm representation from the output of DGELS.
    DVECS(:,:NUM_DVECS,DEPTH+1) = ORTHO_MIN_WORK(:SIZE(DVECS,1),:NUM_DVECS)
  END SUBROUTINE MAKE_DVECS_MIN_NORM

  ! ==================================================================
  FUNCTION MATRIX_RANK(MATRIX)
    ! 4) MATRIX_RANK
    ! 
    !   Get the rank of the provided matrix using the SVD.
    ! 
    ! Input:
    !   MATRIX(:,:) -- Real dense matrix.
    ! 
    ! Output:
    !   MATRIX_RANK -- The integer rank of MATRIX.
    ! 
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    ! Local variables for computing the orthogonal vector.
    REAL(KIND=R8), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2))      :: A
    REAL(KIND=R8), DIMENSION(MIN(SIZE(MATRIX,1),SIZE(MATRIX,2))) :: S
    INTEGER :: IDX, MATRIX_RANK
    ! Unused DGESVD parameters.
    REAL(KIND=R8), DIMENSION(1) :: U, VT

    A = MATRIX
    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N','N',SIZE(A,1),SIZE(A,2),A,SIZE(A,1),S,U,SIZE(U), &
         VT, SIZE(VT), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
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
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 43; RETURN; END IF

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
    ! 5) MATRIX_DET
    ! 
    !   Compute the determinant of a matrix without modifying it
    !   using the LU decomposition (and a copy).
    ! 
    ! Input:
    !   MATRIX(:,:) -- Real dense matrix.
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
       ERROR = 100*INFO + 54
       RETURN
    END IF

    ! Compute the determinant (product of diagonal of U).
    MATRIX_DET = 1.
    DO IDX = 1, MIN(SIZE(M,1), SIZE(M,2))
       MATRIX_DET = MATRIX_DET * M(IDX,IDX)
    END DO
  END FUNCTION MATRIX_DET

  ! ==================================================================
  FUNCTION MATRIX_CONDITION_INV(MATRIX) RESULT(RCOND)
    ! 6) MATRIX_CONDITION_INV
    ! 
    ! Compute the condition number (for testing near rank deficiency)
    ! using DGECON (which computes 1 / CONDITION).
    ! 
    ! Input:
    !   MATRIX(:,:) -- Real dense matrix.
    ! 
    ! Output:
    !   RCOND -- Real value corresponding to the reciprocal of the
    !            condition number of the provided matrix.
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=R8) :: RCOND
    ! Local arrays for work.
    REAL(KIND=R8), DIMENSION(4*SIZE(MATRIX,2)) :: WORK
    INTEGER,       DIMENSION(  SIZE(MATRIX,2)) :: IWORK
    ! Use LAPACK
    CALL DGECON('I', SIZE(MATRIX,2), MATRIX, SIZE(MATRIX,1),&
         SUM(ABS(MATRIX)), RCOND, WORK, IWORK, INFO)
    ! Store output of 'info' flag.
    IF (INFO .NE. 0) ERROR = 100 * INFO + 65
  END FUNCTION MATRIX_CONDITION_INV

  ! ==================================================================
  SUBROUTINE NONZERO(ARRAY, NE_ZERO)
    ! 7) NONZERO
    ! 
    ! Return a new array of the indices of 'ARRAY' that contain
    ! nonzero elements. Set error and return array with 1 if ALL(ARRAY==0).
    ! 
    ! Input:
    !   ARRAY(:)   -- Integer array of numbers.
    ! 
    ! Output:
    !   NE_ZERO(:)    -- Integer array corresponding to those indices of
    !                    ARRAY that contain nonzero elements.
    ! 
    INTEGER, INTENT(IN), DIMENSION(:) :: ARRAY
    INTEGER, INTENT(OUT), DIMENSION(:) :: NE_ZERO
    ! Local variables
    INTEGER :: IDX
    ! Identify nonzero elements of ARRAY.
    NE_ZERO = 0
    NUM_DVECS = 0
    DO IDX = 1, SIZE(ARRAY)
       IF (ARRAY(IDX) .NE. 0) THEN
          NUM_DVECS = NUM_DVECS + 1
          NE_ZERO(NUM_DVECS) = IDX
       END IF
    END DO
    IF (NUM_DVECS .LE. 0) THEN; ERROR = 70; END IF
  END SUBROUTINE NONZERO

  ! ==================================================================
  SUBROUTINE RESERVE_MEMORY()
    ! 8) RESERVE_MEMORY
    ! 
    ! Return an allocated real array that has a size equal to the
    ! maximum requested LAPACK work array size across all routines
    ! that will be executed in the evaluation of a box-spline.
    ! 
    ! Output:
    !   WORK -- Real array with size large enough to accommodate all
    !           LAPACK subroutines that will be used to evaluated a
    !           box-spline given the dimension and number of vectors.
    ! 
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: WORK
    REAL(KIND=R8), DIMENSION(3) :: SIZES
    ! Initialize all sizes to 0.
    SIZES = 0_R8
    ! Retrieve expected size from each of the LAPACK routines.

    ! Query the size of the work array for DGESVD. (MATRIX_ORTHOGONAL)
    CALL DGESVD('N', 'A', NUM_DVECS, DIM, LAPACK_WORK, NUM_DVECS, &
         LAPACK_WORK, LAPACK_WORK, 1, LAPACK_WORK, NUM_DVECS, &
         SIZES(1:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 81; RETURN; END IF

    ! Query the size of the work array to (MATRIX_RANK).
    CALL DGESVD('N', 'N', NUM_DVECS, DIM, LAPACK_WORK, NUM_DVECS, &
         LAPACK_WORK, LAPACK_WORK, 1, LAPACK_WORK, 1, SIZES(2:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 81; RETURN; END IF

    ! Get the size of the work array for DGELS (MAKE_DVECS_MIN_NORM).
    CALL DGELS('T', DIM, NUM_DVECS, NUM_DVECS, LAPACK_WORK, DIM, &
         LAPACK_WORK, NUM_DVECS, SIZES(3:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 81; RETURN; END IF

    ! Allocate the work array by rounding the max value in 'SIZES'.
    ALLOCATE(LAPACK_WORK(1:INT(.5_R8 + MAXVAL(SIZES))))

    ! Allocate various global variables for working.
    DEPTH = SUM(DVEC_MULTS) - DIM + 1 ! <- Max recursion depth.
    ALLOCATE(&
         SHIFTED_EVAL_PTS (NUM_DVECS, NUM_PTS,   DEPTH), &
         DVECS            (DIM,       NUM_DVECS, DEPTH),&
         LOC              (NUM_DVECS, DEPTH+1), &
         MULTS            (NUM_DVECS, DEPTH), &
         EVALS_AT_PTS     (NUM_PTS,   DEPTH+1), &
         TEMP_EVAL_PTS    (DIM, NUM_PTS), &
         REMAINING_DVECS  (NUM_DVECS), &
         NONZERO_DVECS    (NUM_DVECS), &
         NORMAL_VECTOR    (DIM), &
         PT_SHIFT         (DIM), &
         )

  END SUBROUTINE RESERVE_MEMORY

  ! ==================================================================
  SUBROUTINE FREE_MEMORY()
    ! Deallocate storage reserved for box-spline evaluation.
    IF (ALLOCATED(LAPACK_WORK))      DEALLOCATE(LAPACK_WORK)
    IF (ALLOCATED(SHIFTED_EVAL_PTS)) DEALLOCATE(SHIFTED_EVAL_PTS)
    IF (ALLOCATED(DVECS))            DEALLOCATE(DVECS)
    IF (ALLOCATED(LOC))              DEALLOCATE(LOC)
    IF (ALLOCATED(MULTS))            DEALLOCATE(MULTS)
    IF (ALLOCATED(EVALS_AT_PTS))     DEALLOCATE(EVALS_AT_PTS)
    IF (ALLOCATED(TEMP_EVAL_PTS))    DEALLOCATE(TEMP_EVAL_PTS)
    IF (ALLOCATED(REMAINING_DVECS))  DEALLOCATE(REMAINING_DVECS)
    IF (ALLOCATED(NONZERO_DVECS))    DEALLOCATE(NONZERO_DVECS)
    IF (ALLOCATED(NORMAL_VECTOR))    DEALLOCATE(NORMAL_VECTOR)
    IF (ALLOCATED(PT_SHIFT))         DEALLOCATE(PT_SHIFT)
  END SUBROUTINE FREE_MEMORY


END SUBROUTINE BOXSPLEV

