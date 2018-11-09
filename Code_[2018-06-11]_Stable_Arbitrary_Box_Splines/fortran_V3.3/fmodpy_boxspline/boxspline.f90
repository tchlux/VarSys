! This file (boxspline.f90) contains the subroutine BOXSPLEV that can
! be used to evaluate a box-spline, defined by its associated
! direction vector set, at given evaluation points, as well as the
! module REAL_PRECISION defining the real precision.
! 
! ====================================================================
! MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
!   INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
! END MODULE REAL_PRECISION

SUBROUTINE BOXSPLEV(DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
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
  !     O( (k)^{n - s} ),
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
  ! [1] Kobbelt, Leif. "Stable Evaluation of Box-Splines." 
  !     Numerical Algorithms 14.4 (1997): 377-382.
  ! 
  ! 
  ! The subroutine BOXSPLEV utilizes LAPACK routines DGECON, DGELS,
  ! DGESVD, DGETRF, and BLAS routine DGEMM.
  ! 
  ! On input:
  ! 
  ! DVECS(:,:) 
  !    is a real S x M array whose columns are the unique direction
  !    vectors used in defining the box-spline.  S <= M, and 
  !    DVECS(1:S,1:S) must be invertible.
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
  !      2  Error computing a minimum norm representation of direction vectors.
  !      3  Error copying direction vectors with nonzero multiplicity.
  !      4  Error computing an orthogonal vector.
  !      5  Error computing the reciprocal condition number of matrix.
  !      6  Error computing a matrix rank.                            
  !      7  Error computing a matrix determinant.                     
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
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: DVECS
  INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
  INTEGER,       INTENT(OUT)                 :: ERROR
  ! Reusable recursion variables and global variables.
  INTEGER :: DIM, NUM_DVECS, NUM_PTS, INFO
  REAL(KIND=R8) :: POSITION
  REAL(KIND=R8), PARAMETER :: SQRTEPS = SQRT(EPSILON(1.0_R8))
  REAL(KIND=R8), PARAMETER :: ONE  = 1.0_R8
  REAL(KIND=R8), PARAMETER :: ZERO = 0.0_R8
  ! ------------------------------------------------------------------
  INTEGER,       DIMENSION(SIZE(DVECS,2),  SIZE(DVECS,2)) :: IDENTITY
  INTEGER,       DIMENSION(SIZE(DVECS,1),  SIZE(DVECS,1)) :: DET_WORK
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,2),  SIZE(DVECS,2)) :: ORTHO_MIN_WORK
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1)-1,SIZE(DVECS,1)) :: TRANS_DVECS
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1),  SIZE(EVAL_PTS,2)) :: TEMP_PTS
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1))  :: PT_SHIFT
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1))  :: SING_VALS
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1))  :: NORMAL_VECTOR
  INTEGER,       DIMENSION(SIZE(DVECS,2))  :: REMAINING_DVECS
  INTEGER,       DIMENSION(SIZE(DVECS,2))  :: NONZERO_DVECS
  REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: LAPACK_WORK
  ! Unique recursion variables for evaluating box-spline. Last dimension is depth.
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,1),SIZE(DVECS,2)   ) :: MIN_NORM_DVECS
  INTEGER,       DIMENSION(SIZE(DVECS,2)                 ) :: LOC
  REAL(KIND=R8), DIMENSION(SIZE(DVECS,2),SIZE(EVAL_PTS,2)) :: SHIFTED_PTS
  ! ------------------------------------------------------------------
  ! Local variables for preparation.
  INTEGER :: IDX_1, IDX_2
  ! Initialize error and info flags.
  ERROR = 0; INFO = 0; 
  ! Initialize parameters 
  ! Store 'global' constants for box-spline evaluation.
  DIM       = SIZE(DVECS, 1)
  NUM_DVECS = SIZE(DVECS, 2)
  NUM_PTS   = SIZE(EVAL_PTS, 2)
  ! Check for usage errors (dimension mismatches, invalid multiplicity)
  IF (SIZE(DVEC_MULTS)   .NE. NUM_DVECS) THEN; ERROR = 1; RETURN; END IF
  IF (SIZE(EVAL_PTS,1)   .NE. DIM)       THEN; ERROR = 2; RETURN; END IF
  IF (SIZE(BOX_EVALS)    .NE. NUM_PTS)   THEN; ERROR = 3; RETURN; END IF
  IF (MINVAL(DVEC_MULTS) .LT. 1)         THEN; ERROR = 4; RETURN; END IF
  ! Check uniqueness of DVECS columns.
  DO IDX_1 = 1, NUM_DVECS-1
     DO IDX_2 = IDX_1+1, NUM_DVECS
        IF (SUM(ABS(DVECS(:,IDX_1) - DVECS(:,IDX_2))) .LT. SQRTEPS) THEN
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
     IF (MATRIX_CONDITION_INV(DVECS(1:DIM, 1:DIM)) .LT. SQRTEPS) THEN
        ERROR = 6; RETURN
     ENDIF
  ENDIF
  ! Allocate a work array large enough for all LAPACK calls.
  IF (ERROR .NE. 0) THEN; RETURN; END IF
  ! Create an identity matrix (for easy matrix-vector multiplication).
  IDENTITY = 0
  FORALL (IDX_1 = 1:NUM_DVECS) IDENTITY(IDX_1,IDX_1) = 1
  ! Calculate the DVECS for the beginning of box-spline evaluation.
  MIN_NORM_DVECS = DVECS;
  CALL MAKE_DVECS_MIN_NORM(MIN_NORM_DVECS)
  IF (ERROR .NE. 0) THEN; RETURN; END IF
  ! Compute the shifted evaluation points.
  CALL DGEMM('T', 'N', SIZE(DVECS,2), SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
       ONE, MIN_NORM_DVECS, SIZE(DVECS,1), EVAL_PTS, SIZE(EVAL_PTS,1), &
       ZERO, SHIFTED_PTS, SIZE(SHIFTED_PTS,1))
  ! Recursively evaluate the box-spline.
  CALL EVALUATE_BOX_SPLINE(MIN_NORM_DVECS, LOC, DVEC_MULTS, SHIFTED_PTS, BOX_EVALS)
  IF (ERROR .NE. 0) THEN; RETURN; END IF
CONTAINS

  ! ==================================================================
  RECURSIVE SUBROUTINE EVALUATE_BOX_SPLINE(C_DVECS, C_LOC, C_MULTS, &
       C_SHIFTED_PTS, C_EVALS)
    ! 1) EVALUATE_BOX_SPLINE:
    !   
    !   Evaluate the box-spline defined by "C_DVECS" recursively, where
    !   this iteration handles the remaining direction vectors, at all
    !   points in "SHIFTED_PTS" and store box-spline evaluations
    !   in "C_EVALS".
    ! 
    !   "C_" stands for "current"
    !   "N_" stands for "next"
    ! 
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: C_DVECS
    INTEGER,       INTENT(IN),  DIMENSION(:)   :: C_LOC, C_MULTS
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: C_SHIFTED_PTS
    REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: C_EVALS
    INTEGER,       DIMENSION(SIZE(C_LOC))      :: N_LOC, N_MULTS
    REAL(KIND=R8), DIMENSION(SIZE(C_EVALS))    :: N_EVALS
    REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: N_DVECS
    ! Local variables
    INTEGER :: IDX_1, IDX_2, IDX_3
    ! Recursion case ...
    IF (SUM(C_MULTS) > DIM) THEN
       C_EVALS = 0._R8
       IDX_2 = 1
       ! Sum over all direction vectors.
       DO IDX_1 = 1, SIZE(C_MULTS,1)
          ! Update multiplicity of directions and position in recursion tree.
          N_MULTS = C_MULTS - IDENTITY(:,IDX_1) ! Reduce multiplicity.
          ! Make recursive calls.
          IF (C_MULTS(IDX_1) .GT. 1) THEN
             ! Copy LOC, DVECS, and SHIFTED_PTS down one level for next recursion.
             ! Perform recursion with only reduced multiplicity.
             CALL EVALUATE_BOX_SPLINE( &
                  C_DVECS, C_LOC, N_MULTS, &
                  C_SHIFTED_PTS, N_EVALS)
             IF (ERROR .NE. 0) RETURN
             C_EVALS = C_EVALS + &
                  N_EVALS * C_SHIFTED_PTS(IDX_2,:)
             ! Perform recursion with transformed set of direction
             ! vectors and evaluation points.
             PT_SHIFT_R = MATMUL(DVECS(:,IDX_1), C_DVECS)
             compute_shift_1 : DO IDX_3 = 1, NUM_PTS
                TEMP_PTS(:,IDX_3) = &
                     C_SHIFTED_PTS(:,IDX_3) - PT_SHIFT_R
             END DO compute_shift_1
             ! Update location for next level of recursion.
             N_LOC = C_LOC + IDENTITY(:,IDX_1) 
             CALL EVALUATE_BOX_SPLINE( &
                  C_DVECS, N_LOC, N_MULTS, &
                  TEMP_PTS, N_EVALS)
             IF (ERROR .NE. 0) RETURN
             C_EVALS = C_EVALS + N_EVALS * &
                  (C_MULTS(IDX_1) - C_SHIFTED_PTS(IDX_2,:))
             IDX_2 = IDX_2 + 1
          ELSE IF (C_MULTS(IDX_1) .GT. 0) THEN
             ! Get the remaining direction vectors.
             N_DVECS = PACK_DVECS(N_MULTS, DVECS)
             IF (ERROR .NE. 0) RETURN
             IF (FULL_RANK(TRANSPOSE(N_DVECS))) THEN
                IF (ERROR .NE. 0) RETURN
                ! Update Least norm representation of the direction vectors.
                CALL MAKE_DVECS_MIN_NORM(N_DVECS)
                IF (ERROR .NE. 0) RETURN
                ! Perform recursion with only reduced multiplicity.
                ! Compute the appropriate point shift according the new DVECS.
                PT_SHIFT = MATMUL(DVECS, C_LOC)
                DO IDX_3 = 1, NUM_PTS
                   TEMP_PTS(:,IDX_3) = PTS(:,IDX_3) - PT_SHIFT
                END DO
                ! Perform recursion with only reduced multiplicity.
                CALL EVALUATE_BOX_SPLINE( &
                     N_DVECS, C_LOC, N_MULTS, &
                     MATMUL(TRANSPOSE(N_DVECS), TEMP_PTS), N_EVALS)
                IF (ERROR .NE. 0) RETURN
                C_EVALS = C_EVALS + &
                     N_EVALS * C_SHIFTED_PTS(IDX_2,:)
                ! Update location for next level of recursion.
                N_LOC = C_LOC + IDENTITY(:,IDX_1) 
                PT_SHIFT = MATMUL(DVECS, N_LOC)
                compute_shift_3 : DO IDX_3 = 1, NUM_PTS
                   TEMP_PTS(:,IDX_3) = PTS(:,IDX_3) - PT_SHIFT
                END DO compute_shift_3
                ! Perform recursion with transformed set of direction vectors.
                CALL EVALUATE_BOX_SPLINE( &
                     N_DVECS, N_LOC, N_MULTS, &
                     MATMUL(TRANSPOSE(N_DVECS), TEMP_PTS), N_EVALS)
                IF (ERROR .NE. 0) RETURN
                C_EVALS = C_EVALS + &
                     N_EVALS * &
                     (C_MULTS(IDX_1) - C_SHIFTED_PTS(IDX_2,:))
             END IF
             IDX_2 = IDX_2 + 1
          END IF
       END DO
       ! Normalize by number of direction vectors in computation.
       C_EVALS = C_EVALS / &
            REAL(SUM(C_MULTS) - DIM, R8)
    ELSE
       ! Base case ... compute characteristic function.
       C_EVALS = 1.0_R8
       ! Pack the unique direction vectors into the memory location
       ! for the current set of direction vectors (since the 'current'
       ! are not needed for base case evaluation).
       N_DVECS = PACK_DVECS(C_MULTS, DVECS)
       IF (ERROR .NE. 0) RETURN
       ! Delayed translations (this is what makes the algorithm more stable).
       PT_SHIFT = MATMUL(DVECS, C_LOC)
       compute_shift_4 : DO IDX_1 = 1, NUM_PTS
          TEMP_PTS(:,IDX_1) = PTS(:,IDX_1) - PT_SHIFT
       END DO compute_shift_4
       ! Check evaluation point locations against all remaining direction vectors.
       DO IDX_1 = 1, DIM
          ! Get the active set of direction vectors (exluding current vector).
          TRANS_DVECS = TRANSPOSE(N_DVECS(:,:DIM-1))
          IF (IDX_1 .LT. DIM) TRANS_DVECS(IDX_1,:) = N_DVECS(:,DIM)
          ! Calculate the orthogonal vector to the remaining direction vectors.
          CALL COMPUTE_ORTHOGONAL(TRANS_DVECS, NORMAL_VECTOR)
          IF (ERROR .NE. 0) RETURN
          ! Compute shifted position (relative to normal vector).
          POSITION = DOT_PRODUCT(N_DVECS(:,IDX_1), NORMAL_VECTOR)
          ! Compute shifted evaluation locations. (1, DIM) x (DIM, NUM_PTS)
          N_EVALS = MATMUL(NORMAL_VECTOR, TEMP_PTS)
          ! Identify those points that are outside of this box (0-side).
          IF (POSITION .GT. 0) THEN
             WHERE (N_EVALS .LT. 0) C_EVALS = 0._R8
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (N_EVALS .GE. 0) C_EVALS = 0._R8
          END IF
          ! Recompute shifted location (other side of box) based on selected direction vector.
          compute_shifted_point_2 : DO IDX_2 = 1, NUM_PTS
             N_EVALS(IDX_2) = DOT_PRODUCT( &
                  PTS(:,IDX_2) - PT_SHIFT - N_DVECS(:,IDX_1), &
                  NORMAL_VECTOR)
          END DO compute_shifted_point_2
          ! Identify those shifted points that are outside of this box on REMAINING_DVEC(IDX_1)-side.
          IF (POSITION .GT. 0) THEN
             WHERE (N_EVALS .GE. 0) C_EVALS = 0._R8
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (N_EVALS .LT. 0) C_EVALS = 0._R8
          END IF
       END DO
       ! Normalize evaluations by determinant of box.
       C_EVALS = C_EVALS / &
            ABS(MATRIX_DET(TRANSPOSE(N_DVECS(:,:DIM))))
       IF (ERROR .NE. 0) RETURN
    END IF
  END SUBROUTINE EVALUATE_BOX_SPLINE

  !===================================================================
  !             Supporting code for computing box-spline              
  !===================================================================

  ! ==================================================================
  SUBROUTINE MAKE_DVECS_MIN_NORM(DVECS)
    ! 2) MAKE_DVECS_MIN_NORM
    ! 
    !   Compute the minimum norm representation of 'MATRIX' and store
    !   it in 'MIN_NORM', use DGELS to find the least squares
    !   solution to the problem (AX = I). This is a more numerically
    !   stable solution to the linear system (A^T A) X = A^T.
    REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: DVECS
    INTEGER :: IDX
    ! Make "ORTHO_MIN_WORK" the identity matrix
    ORTHO_MIN_WORK = 0._R8
    FORALL (IDX = 1:SIZE(DVECS,2)) ORTHO_MIN_WORK(IDX,IDX) = 1.0_R8
    ! Call DGELS for actual solve.
    CALL DGELS('T', SIZE(DVECS,1), SIZE(DVECS,2), SIZE(DVECS,2),&
         DVECS, SIZE(DVECS,1), ORTHO_MIN_WORK, &
         SIZE(ORTHO_MIN_WORK,1), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    ! Check for error.
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 22; RETURN; END IF
    ! Extract the minimum norm representation from the output of DGELS.
    DVECS = ORTHO_MIN_WORK(:SIZE(DVECS,1),:SIZE(DVECS,2))
  END SUBROUTINE MAKE_DVECS_MIN_NORM

  ! ==================================================================
  FUNCTION PACK_DVECS(MULTS, SRCE_DVECS) RESULT(DEST_DVECS)
    ! 3) PACK_DVECS
    ! 
    !   Given multiplicities and source direction vectors, pack all
    !   of the nonzero-multiplicity source direction vectors into the
    !   storage location provided for remaining direction vectors.
    ! 
    ! Inputs:
    !   MULTS(:)        -- Integer array of multiplicities.
    !   SRCE_DVECS(:,:) -- Real dense matrix of (source) direction
    !                      vectors.
    !   DEST_DVECS(:,:) -- Real dense matrix for storing direction
    !                      vectors.
    ! 
    INTEGER,       INTENT(IN),  DIMENSION(:)   :: MULTS
    REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: SRCE_DVECS
    REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: DEST_DVECS
    INTEGER :: IDX, NUM_DVECS
    ! Loop through once to find out how many DVECS there will be.
    NUM_DVECS = 0
    DO IDX = 1, SIZE(MULTS)
       IF (MULTS(IDX) .GT. 0) (NUM_DVECS = NUM_DVECS + 1)
    END DO
    ! Allocate the new direction vectors and copy them over.
    ALLOCATE(DEST_DVECS(1:SIZE(SRCE_DVECS,1), 1:NUM_DVECS))
    NUM_DVECS = 0
    DO IDX = 1, SIZE(MULTS)
       IF (MULTS(IDX) .GT. 0) THEN
          NUM_DVECS = NUM_DVECS + 1
          DEST_DVECS(:,NUM_DVECS) = SRCE_DVECS(:,IDX)
       END IF
    END DO
  END SUBROUTINE PACK_DVECS

  ! ==================================================================
  SUBROUTINE COMPUTE_ORTHOGONAL(A, ORTHOGONAL)
    ! 4) COMPUTE_ORTHOGONAL
    ! 
    !   Given a matrix A of row vectors, compute a vector orthogonal to
    !   the row vectors in A and store it in ORTHOGONAL using DGESVD.
    !   If there any near-zero singular values are identified, the
    !   vector associated with the smallest such singular values is
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
    ! Local variables.
    REAL(KIND=R8) :: TOL
    REAL(KIND=R8), DIMENSION(1) :: U
    ! Use the SVD to get the orthogonal vector(s).
    CALL DGESVD('N', 'A', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), &
         SING_VALS, U, SIZE(U), ORTHO_MIN_WORK(:SIZE(A,2),:SIZE(A,2)), &
         SIZE(A,2), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 43; RETURN; END IF
    ! Compute the appropriate tolerance based on calculations.
    TOL = EPSILON(1.0_R8) * MAXVAL(SHAPE(A)) * MAXVAL(SING_VALS)
    ORTHOGONAL = 0._R8
    ! Check for a vector in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value).
    IF (SING_VALS(SIZE(SING_VALS)-1) .LE. TOL) THEN
       ORTHOGONAL = ORTHO_MIN_WORK(SIZE(SING_VALS)-1,:SIZE(A,2))
    ! If no orthogonal vector was found and the matrix does not contain
    ! enough vectors to span the space, use the first vector in VT at
    ! (RANK(A) + 1).
    ELSE IF (SIZE(A,2) > SIZE(SING_VALS)-1) THEN
       ORTHOGONAL = ORTHO_MIN_WORK(SIZE(SING_VALS),:SIZE(A,2))
    END IF
  END SUBROUTINE COMPUTE_ORTHOGONAL

  ! ==================================================================
  FUNCTION MATRIX_CONDITION_INV(MATRIX) RESULT(RCOND)
    ! 5) MATRIX_CONDITION_INV
    ! 
    !   Compute the condition number (for testing near rank deficiency)
    !   using DGECON (which computes 1 / CONDITION).
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
    ! Use LAPACK to compue recirpocal of matrix codition number.
    CALL DGECON('I', SIZE(MATRIX,2), MATRIX, SIZE(MATRIX,1),&
         SUM(ABS(MATRIX)), RCOND, WORK, IWORK, INFO)
    ! Check for errors during exeuction.
    IF (INFO .NE. 0) ERROR = 100 * INFO + 55
  END FUNCTION MATRIX_CONDITION_INV

  ! ==================================================================
  FUNCTION FULL_RANK(MATRIX)
    ! 6) FULL_RANK
    ! 
    !   Get the rank of the provided matrix using the SVD.
    ! 
    ! Input:
    !   MATRIX(:,:) -- Real dense matrix.
    ! 
    ! Output:
    !   MATRIX_RANK -- The integer rank of MATRIX.
    ! 
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: MATRIX
    LOGICAL :: FULL_RANK
    REAL(KIND=R8) :: TOL
    ! Early return if possible, when we know the rank must be low.
    IF (SIZE(MATRIX,1) .LT. SIZE(SING_VALS)) THEN; FULL_RANK = .FALSE.; RETURN; END IF
    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N', 'N', SIZE(MATRIX,1), SIZE(MATRIX,2), MATRIX, &
         SIZE(MATRIX,1), SING_VALS, NULL(), 0, NULL(), 0, &
         LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 63; RETURN; END IF
    ! Compute a reasonable singular value tolerance based on expected
    ! numerical error ["Numerical Recipes" by W. H. Press et al., Matlab].
    TOL = EPSILON(1.0_R8) * MAXVAL(SHAPE(MATRIX)) * MAXVAL(SING_VALS)
    FULL_RANK = (SING_VALS(SIZE(SING_VALS)) > TOL)
  END FUNCTION FULL_RANK

  ! ==================================================================
  FUNCTION MATRIX_DET(MATRIX)
    ! 7) MATRIX_DET
    ! 
    !   Compute the determinant of a matrix using the LU decomposition.
    ! 
    ! Input:
    !   MATRIX(:,:) -- Real dense matrix.
    ! 
    ! Output:
    !   MATRIX_DET -- The real-valued determinant of MATRIX.
    ! 
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=R8) :: MATRIX_DET
    INTEGER :: IDX
    ! Do the LU decomposition.
    CALL DGETRF(SIZE(MATRIX,1), SIZE(MATRIX,2), MATRIX, SIZE(MATRIX,1), DET_WORK, INFO)
    ! Check for errors.
    IF (INFO .NE. 0) THEN
       MATRIX_DET = 1.0_R8 ! <- Set a value that will not break caller.
       ERROR = 100*INFO + 74
       RETURN
    END IF
    ! Compute the determinant (product of diagonal of U).
    MATRIX_DET = 1.
    DO IDX = 1, MIN(SIZE(MATRIX,1), SIZE(MATRIX,2))
       MATRIX_DET = MATRIX_DET * MATRIX(IDX,IDX)
    END DO
  END FUNCTION MATRIX_DET

  ! ==================================================================
  SUBROUTINE PRINT_OPTIMAL_BLOCK_SIZE()
    ! 8) PRINT_OPTIMAL_BLOCK_SIZE
    ! 
    !   Use the LAPACK routine ILAENV to print out the optimal block
    !   size for the current architecture (requires optimal LAPACK and
    !   BLAS installation).
    ! 
    INTEGER :: BLOCK_SIZE
    INTERFACE
       FUNCTION ILAENV(ISPEC, NAME, OPTS, N1, N2, N3, N4)
         INTEGER :: ISPEC
         CHARACTER :: NAME, OPTS
         INTEGER, OPTIONAL :: N1,N2,N3,N4
         INTEGER :: ILAENV
       END FUNCTION ILAENV
    END INTERFACE
    ! 
    ! Compute storage needed by each of the LAPACK routines.
    ! DGELS optimal work array size:
    !   DIM * NUM_DVECS * (1 +  BLOCK_SIZE)
    ! DGESVD optimal work array size:
    !   MAX(3*DIM + NUM_DVECS, 5*DIM)
    ! 
    ! BLOCK_SIZE = 65536 ! <- 64KB
    ! BLOCK_SIZE = 8192  ! <- 8KB
    ! 
    BLOCK_SIZE = ILAENV(1, 'DGELS', 'T', DIM, NUM_DVECS, NUM_DVECS)
    PRINT *, 'BLOCK SIZE:', BLOCK_SIZE
  END SUBROUTINE PRINT_OPTIMAL_BLOCK_SIZE

END SUBROUTINE BOXSPLEV
