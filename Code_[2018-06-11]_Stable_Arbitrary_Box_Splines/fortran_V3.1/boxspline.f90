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
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: UNIQUE_DVECS
  INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
  REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
  REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
  INTEGER,       INTENT(OUT)                 :: ERROR
  ! Reusable recursion variables and global variables.
  INTEGER :: DIM, NUM_DVECS, NUM_PTS, DEPTH, INFO
  REAL(KIND=R8) :: POSITION
  REAL(KIND=R8), PARAMETER :: SQRTEPS = SQRT(EPSILON(1.0_R8))
  REAL(KIND=R8), PARAMETER :: ONE  = 1.0_R8
  REAL(KIND=R8), PARAMETER :: ZERO = 0.0_R8
  INTEGER,       DIMENSION(SIZE(UNIQUE_DVECS,2),  SIZE(UNIQUE_DVECS,2)) :: IDENTITY
  INTEGER,       DIMENSION(SIZE(UNIQUE_DVECS,1),  SIZE(UNIQUE_DVECS,1)) :: DET_WORK
  REAL(KIND=R8), DIMENSION(SIZE(UNIQUE_DVECS,2),  SIZE(UNIQUE_DVECS,2)) :: ORTHO_MIN_WORK
  REAL(KIND=R8), DIMENSION(SIZE(UNIQUE_DVECS,1)-1,SIZE(UNIQUE_DVECS,1)) :: TRANS_DVECS
  REAL(KIND=R8), DIMENSION(SIZE(UNIQUE_DVECS,1))                        :: SING_VALS
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: TEMP_EVAL_PTS
  REAL(KIND=R8), DIMENSION(:),   ALLOCATABLE :: PT_SHIFT, PT_SHIFT_R, NORMAL_VECTOR, LAPACK_WORK
  INTEGER,       DIMENSION(:),   ALLOCATABLE :: REMAINING_DVECS, NONZERO_DVECS
  ! Unique recursion variables for evaluating box-spline. Last dimension is depth.
  REAL(KIND=R8), DIMENSION(:,:,:), ALLOCATABLE :: DVECS
  INTEGER,       DIMENSION(:,:),   ALLOCATABLE :: LOC, MULTS
  REAL(KIND=R8), DIMENSION(:,:),   ALLOCATABLE :: EVALS_AT_PTS
  REAL(KIND=R8), DIMENSION(:,:,:), ALLOCATABLE :: SHIFTED_EVAL_PTS
  ! Local variables for preparation.
  INTEGER :: IDX_1, IDX_2

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
  FORALL (IDX_1 = 1:NUM_DVECS) IDENTITY(IDX_1,IDX_1) = 1
  ! Calculate the DVECS for the beginning of box-spline evaluation.
  DEPTH = 0; DVECS(:,:,1) = UNIQUE_DVECS(:,:);
  CALL MAKE_DVECS_MIN_NORM(DVECS(:,:,1))
  IF (ERROR .NE. 0) THEN; CALL FREE_MEMORY(); RETURN; END IF
  ! Initialize the multiplicities and location.
  MULTS(:,1) = DVEC_MULTS(:); LOC(:,1) = 0;
  ! Compute the shifted evaluation points.
  CALL DGEMM('T', 'N', SIZE(DVECS,2), SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
       ONE, DVECS(:,:,1), SIZE(DVECS,1), EVAL_PTS, SIZE(EVAL_PTS,1), &
       ZERO, SHIFTED_EVAL_PTS(:,:,1), SIZE(SHIFTED_EVAL_PTS,1))
  ! Recursively evaluate the box-spline.
  CALL EVALUATE_BOX_SPLINE(&
       DVECS(:,:,1), LOC(:,1), MULTS(:,1), SHIFTED_EVAL_PTS(:,:,1), EVALS_AT_PTS(:,1), &
       DVECS(:,:,2), LOC(:,2), MULTS(:,2), SHIFTED_EVAL_PTS(:,:,2), EVALS_AT_PTS(:,2))
  IF (ERROR .NE. 0) THEN; CALL FREE_MEMORY(); RETURN; END IF
  BOX_EVALS(:) = EVALS_AT_PTS(:,1)
CONTAINS

  ! ==================================================================
  RECURSIVE SUBROUTINE EVALUATE_BOX_SPLINE(&
       R_DVECS, R_LOC, R_MULTS, R_SHIFTED_EVAL_PTS, R_EVALS_AT_PTS, &
       R_NEXT_DVECS, R_NEXT_LOC, R_NEXT_MULTS, R_NEXT_SHIFTED_EVAL_PTS, &
       R_NEXT_EVALS_AT_PTS)
    ! 1) EVALUATE_BOX_SPLINE:
    !   
    !   Evaluate the box-spline defined by "DVECS" recursively, where
    !   this iteration handles the remaining direction vectors, at all
    !   points in "SHIFTED_EVAL_PTS" and store box-spline evaluations
    !   in "EVALS_AT_PTS".
    ! 
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: R_DVECS
    INTEGER,       INTENT(IN), DIMENSION(:)   :: R_LOC, R_MULTS
    REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: R_SHIFTED_EVAL_PTS
    REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: R_NEXT_DVECS
    INTEGER,       INTENT(OUT), DIMENSION(:)   :: R_NEXT_LOC, R_NEXT_MULTS
    REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: R_NEXT_SHIFTED_EVAL_PTS
    REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: R_EVALS_AT_PTS, R_NEXT_EVALS_AT_PTS
    ! Local variables
    INTEGER :: IDX_1, IDX_2, IDX_3, TEMP_NUM_DVECS_1, TEMP_NUM_DVECS_2
    ! Adjust the global variable "DEPTH" to use appropriate memory slices.
    DEPTH = DEPTH + 1
    ! Recursion case ...
    IF (SUM(R_MULTS) > DIM) THEN
       R_EVALS_AT_PTS(:) = 0._R8
       IDX_2 = 1
       ! Sum over all direction vectors.
       DO IDX_1 = 1, SIZE(MULTS,1)
          ! Update multiplicity of directions and position in recursion tree.
          R_NEXT_MULTS(:) = R_MULTS(:) - IDENTITY(:,IDX_1) ! Reduce multiplicity.
          ! Make recursive calls.
          IF (R_MULTS(IDX_1) .GT. 1) THEN
             ! Copy LOC, DVECS, and SHIFTED_EVAL_PTS down one level for next recursion.
             R_NEXT_LOC(:) = R_LOC(:)
             R_NEXT_DVECS(:,:NUM_DVECS) = R_DVECS(:,:NUM_DVECS) ! Pass the DVECS down
             R_NEXT_SHIFTED_EVAL_PTS(:,:) = R_SHIFTED_EVAL_PTS(:,:)
             ! Perform recursion with only reduced multiplicity.
             CALL EVALUATE_BOX_SPLINE( &
                  R_NEXT_DVECS, R_NEXT_LOC, R_NEXT_MULTS, &
                  R_NEXT_SHIFTED_EVAL_PTS, R_NEXT_EVALS_AT_PTS, &
                  DVECS(:,:,DEPTH+2), LOC(:,DEPTH+2), MULTS(:,DEPTH+2), &
                  SHIFTED_EVAL_PTS(:,:,DEPTH+2), EVALS_AT_PTS(:,DEPTH+2))
             IF (ERROR .NE. 0) RETURN
             R_EVALS_AT_PTS(:) = R_EVALS_AT_PTS(:) + &
                  R_NEXT_EVALS_AT_PTS(:) * R_SHIFTED_EVAL_PTS(IDX_2,:)
             ! Perform recursion with transformed set of direction
             ! vectors and evaluation points. (store at next 'DEPTH')
             PT_SHIFT_R(:) = MATMUL(UNIQUE_DVECS(:,IDX_1), R_DVECS(:,:))
             compute_shift_1 : DO IDX_3 = 1, NUM_PTS
                R_NEXT_SHIFTED_EVAL_PTS(:,IDX_3) = &
                     R_SHIFTED_EVAL_PTS(:,IDX_3) - PT_SHIFT_R(:)
             END DO compute_shift_1
             ! Update location for next level of recursion.
             R_NEXT_LOC(:) = R_LOC(:) + IDENTITY(:,IDX_1) 
             CALL EVALUATE_BOX_SPLINE( &
                  R_NEXT_DVECS, R_NEXT_LOC, R_NEXT_MULTS, &
                  R_NEXT_SHIFTED_EVAL_PTS, R_NEXT_EVALS_AT_PTS, &
                  DVECS(:,:,DEPTH+2), LOC(:,DEPTH+2), MULTS(:,DEPTH+2), &
                  SHIFTED_EVAL_PTS(:,:,DEPTH+2), EVALS_AT_PTS(:,DEPTH+2))
             IF (ERROR .NE. 0) RETURN
             R_EVALS_AT_PTS(:) = R_EVALS_AT_PTS(:) + &
                  R_NEXT_EVALS_AT_PTS(:) * &
                  (R_MULTS(IDX_1) - R_SHIFTED_EVAL_PTS(IDX_2,:))
             IDX_2 = IDX_2 + 1
          ELSE IF (R_MULTS(IDX_1) .GT. 0) THEN
             TEMP_NUM_DVECS_1 = NUM_DVECS
             ! Get the remaining direction vectors.
             CALL PACK_DVECS(R_NEXT_MULTS, UNIQUE_DVECS, R_NEXT_DVECS)
             IF (ERROR .NE. 0) RETURN
             ORTHO_MIN_WORK(:NUM_DVECS,:DIM) = TRANSPOSE(R_NEXT_DVECS(:,:NUM_DVECS))
             IF (MATRIX_RANK(ORTHO_MIN_WORK(:NUM_DVECS,:DIM)) .EQ. DIM) THEN
                IF (ERROR .NE. 0) RETURN
                ! Store the number of direction vectors at this level before recursion.
                TEMP_NUM_DVECS_2 = NUM_DVECS
                ! Assign location for next recursion level to be unchanged.
                R_NEXT_LOC(:) = R_LOC(:)
                ! Update Least norm representation of the direction vectors.
                CALL MAKE_DVECS_MIN_NORM(R_NEXT_DVECS(:,:NUM_DVECS))
                IF (ERROR .NE. 0) RETURN
                ! Perform recursion with only reduced multiplicity.
                PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), R_LOC(:))
                compute_shift_2 : DO IDX_3 = 1, NUM_PTS
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - PT_SHIFT(:)
                END DO compute_shift_2
                ! Next dvecs x temp eval points
                CALL DGEMM('T', 'N', NUM_DVECS, SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
                     ONE, R_NEXT_DVECS(:,:NUM_DVECS), SIZE(DVECS,1), TEMP_EVAL_PTS, SIZE(EVAL_PTS,1), &
                     ZERO, R_NEXT_SHIFTED_EVAL_PTS(:NUM_DVECS,:), NUM_DVECS)
                ! Perform recursion with only reduced multiplicity.
                CALL EVALUATE_BOX_SPLINE( &
                     R_NEXT_DVECS, R_NEXT_LOC, R_NEXT_MULTS, &
                     R_NEXT_SHIFTED_EVAL_PTS, R_NEXT_EVALS_AT_PTS, &
                     DVECS(:,:,DEPTH+2), LOC(:,DEPTH+2), MULTS(:,DEPTH+2), &
                     SHIFTED_EVAL_PTS(:,:,DEPTH+2), EVALS_AT_PTS(:,DEPTH+2))
                IF (ERROR .NE. 0) RETURN
                R_EVALS_AT_PTS(:) = R_EVALS_AT_PTS(:) + &
                     R_NEXT_EVALS_AT_PTS(:) * R_SHIFTED_EVAL_PTS(IDX_2,:)
                ! Update location for next level of recursion.
                R_NEXT_LOC(:) = R_LOC(:) + IDENTITY(:,IDX_1) 
                PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), R_NEXT_LOC(:))
                compute_shift_3 : DO IDX_3 = 1, NUM_PTS
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - PT_SHIFT(:)
                END DO compute_shift_3
                ! Recalculate "NUM_DVECS" since it was destroyed in recursion.
                NUM_DVECS = TEMP_NUM_DVECS_2
                ! Next dvecs x temp eval points.
                CALL DGEMM('T', 'N', NUM_DVECS, SIZE(EVAL_PTS,2), SIZE(DVECS,1), &
                     ONE, R_NEXT_DVECS(:,:NUM_DVECS), SIZE(DVECS,1), TEMP_EVAL_PTS, SIZE(EVAL_PTS,1), &
                     ZERO, R_NEXT_SHIFTED_EVAL_PTS(:NUM_DVECS,:), NUM_DVECS)
                ! Perform recursion with transformed set of direction vectors.
                CALL EVALUATE_BOX_SPLINE( &
                     R_NEXT_DVECS, R_NEXT_LOC, R_NEXT_MULTS, &
                     R_NEXT_SHIFTED_EVAL_PTS, R_NEXT_EVALS_AT_PTS, &
                     DVECS(:,:,DEPTH+2), LOC(:,DEPTH+2), MULTS(:,DEPTH+2), &
                     SHIFTED_EVAL_PTS(:,:,DEPTH+2), EVALS_AT_PTS(:,DEPTH+2))
                IF (ERROR .NE. 0) RETURN
                R_EVALS_AT_PTS(:) = R_EVALS_AT_PTS(:) + &
                     R_NEXT_EVALS_AT_PTS(:) * &
                     (R_MULTS(IDX_1) - R_SHIFTED_EVAL_PTS(IDX_2,:))
             END IF
             ! Reset "NUM_DVECS" after executing the above steps.
             NUM_DVECS = TEMP_NUM_DVECS_1
             IDX_2 = IDX_2 + 1
          END IF
       END DO
       ! Normalize by number of direction vectors in computation.
       R_EVALS_AT_PTS(:) = R_EVALS_AT_PTS(:) / &
            REAL(SUM(R_MULTS(:)) - DIM, R8)
    ELSE
       ! Base case ... compute characteristic function.
       R_EVALS_AT_PTS(:) = 1_R8
       ! Pack the unique direction vectors into the memory location
       ! for the current set of direction vectors (since the 'current'
       ! are not needed for base case evaluation).
       CALL PACK_DVECS(R_MULTS(:), UNIQUE_DVECS, R_NEXT_DVECS(:,:))
       IF (ERROR .NE. 0) RETURN
       ! Delayed translations (this is what makes the algorithm more stable).
       PT_SHIFT(:) = MATMUL(UNIQUE_DVECS(:,:), R_LOC(:))
       compute_shift_4 : DO IDX_1 = 1, NUM_PTS
          R_NEXT_SHIFTED_EVAL_PTS(:DIM,IDX_1) = EVAL_PTS(:,IDX_1) - PT_SHIFT(:)
       END DO compute_shift_4
       ! Check evaluation point locations against all remaining direction vectors.
       DO IDX_1 = 1, DIM
          ! Get the active set of direction vectors (exluding current vector).
          TRANS_DVECS(:,:) = TRANSPOSE(R_NEXT_DVECS(:,:DIM-1))
          IF (IDX_1 .LT. DIM) TRANS_DVECS(IDX_1,:) = R_NEXT_DVECS(:,DIM)
          ! Calculate the orthogonal vector to the remaining direction vectors.
          CALL COMPUTE_ORTHOGONAL(TRANS_DVECS, NORMAL_VECTOR)
          IF (ERROR .NE. 0) RETURN
          ! Compute shifted position (relative to normal vector).
          POSITION = DOT_PRODUCT(R_NEXT_DVECS(:,IDX_1), NORMAL_VECTOR(:))
          ! Compute shifted evaluation locations. (1, DIM) x (DIM, NUM_PTS)
          R_NEXT_EVALS_AT_PTS(:) = MATMUL(NORMAL_VECTOR, R_NEXT_SHIFTED_EVAL_PTS(:DIM,:))
          ! Identify those points that are outside of this box (0-side).
          IF (POSITION .GT. 0) THEN
             WHERE (R_NEXT_EVALS_AT_PTS(:) .LT. 0) R_EVALS_AT_PTS(:) = 0._R8
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (R_NEXT_EVALS_AT_PTS(:) .GE. 0) R_EVALS_AT_PTS(:) = 0._R8
          END IF
          ! Recompute shifted location (other side of box) based on selected direction vector.
          compute_shifted_point_2 : DO IDX_2 = 1, NUM_PTS
             R_NEXT_EVALS_AT_PTS(IDX_2) = DOT_PRODUCT( &
                  EVAL_PTS(:,IDX_2) - PT_SHIFT(:) - R_NEXT_DVECS(:,IDX_1), &
                  NORMAL_VECTOR(:))
          END DO compute_shifted_point_2
          ! Identify those shifted points that are outside of this box on REMAINING_DVEC(IDX_1)-side.
          IF (POSITION .GT. 0) THEN
             WHERE (R_NEXT_EVALS_AT_PTS(:) .GE. 0) R_EVALS_AT_PTS(:) = 0._R8
          ELSE IF (POSITION .LT. 0) THEN
             WHERE (R_NEXT_EVALS_AT_PTS(:) .LT. 0) R_EVALS_AT_PTS(:) = 0._R8
          END IF
       END DO
       ! Normalize evaluations by determinant of box.
       R_EVALS_AT_PTS(:) = R_EVALS_AT_PTS(:) / &
            ABS(MATRIX_DET(TRANSPOSE(R_NEXT_DVECS(:,:DIM))))
       IF (ERROR .NE. 0) RETURN
    END IF
    DEPTH = DEPTH - 1
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
    ORTHO_MIN_WORK(:,:) = 0._R8
    FORALL (IDX = 1:NUM_DVECS) ORTHO_MIN_WORK(IDX,IDX) = 1_R8
    ! Call DGELS for actual solve.
    CALL DGELS('T', SIZE(DVECS,1), SIZE(DVECS,2), SIZE(UNIQUE_DVECS,2),&
         DVECS, SIZE(DVECS,1), ORTHO_MIN_WORK, &
         SIZE(ORTHO_MIN_WORK,1), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    ! Check for error.
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 22; RETURN; END IF
    ! Extract the minimum norm representation from the output of DGELS.
    DVECS(:,:) = ORTHO_MIN_WORK(:SIZE(DVECS,1),:SIZE(DVECS,2))
  END SUBROUTINE MAKE_DVECS_MIN_NORM

  ! ==================================================================
  SUBROUTINE PACK_DVECS(MULTS, SRCE_DVECS, DEST_DVECS)
    ! 3) PACK_DVECS
    ! 
    !   Given multiplicities and source direction vectors, pack all
    !   of the nonzero-multiplicity source direction vectors into the
    !   storage location provided for remaining direction vectors.
    !   Update global variable "NUM_DVECS" in the process.
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
    REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: DEST_DVECS
    INTEGER :: IDX
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
    TOL = EPSILON(1._R8) * MAXVAL(SHAPE(A)) * MAXVAL(SING_VALS)
    ORTHOGONAL(:) = 0._R8
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
  FUNCTION MATRIX_RANK(MATRIX)
    ! 6) MATRIX_RANK
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
    ! Local variables for computing the orthogonal vector.
    INTEGER :: IDX, MATRIX_RANK
    REAL(KIND=R8) :: TOL
    ! Unused DGESVD parameters.
    REAL(KIND=R8), DIMENSION(1) :: U, VT
    ! Early return if possible, when we know the rank must be low.
    IF (SIZE(MATRIX,1) .LT. SIZE(SING_VALS)) THEN; MATRIX_RANK = 0; RETURN; END IF
    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N', 'N', SIZE(MATRIX,1), SIZE(MATRIX,2), MATRIX, SIZE(MATRIX,1), SING_VALS,&
         U, SIZE(U), VT, SIZE(VT), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 63; RETURN; END IF
    ! Compute a reasonable singular value tolerance based on expected
    ! numerical error ["Numerical Recipes" by W. H. Press et al., Matlab].
    TOL = EPSILON(1._R8) * MAXVAL(SHAPE(MATRIX)) * MAXVAL(SING_VALS)
    ! Find the first near-0 singular value starting from smallest
    ! value, assumes high rank is more likely than low rank.
    find_null : DO IDX = SIZE(SING_VALS), 1, -1
       IF (SING_VALS(IDX) .GT. TOL) THEN
          MATRIX_RANK = IDX
          EXIT find_null
       END IF
    END DO find_null
  END FUNCTION MATRIX_RANK

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
       MATRIX_DET = 1._R8 ! <- Set a value that will not break caller.
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
  SUBROUTINE RESERVE_MEMORY()
    ! 8) RESERVE_MEMORY
    ! 
    !   Allocate a real array that has a size equal to the
    !   maximum requested LAPACK work array size across all routines
    !   that will be executed in the evaluation of a box-spline.
    !   Allocate all global arrays for recursive evaluation of box-spline.
    ! 
    REAL(KIND=R8), DIMENSION(:), ALLOCATABLE :: WORK
    REAL(KIND=R8), DIMENSION(3) :: SIZES
    ! Initialize all sizes to 0.
    SIZES = 0._R8
    ! Compute storage needed by each of the LAPACK routines.

    ! Get the size of the work array for DGELS (MAKE_DVECS_MIN_NORM).
    CALL DGELS('T', DIM, NUM_DVECS, NUM_DVECS, LAPACK_WORK, DIM, &
         LAPACK_WORK, NUM_DVECS, SIZES(3:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 81; RETURN; END IF
    ! Query the size of the work array to (MATRIX_RANK).
    CALL DGESVD('N', 'N', NUM_DVECS, DIM, LAPACK_WORK, NUM_DVECS, &
         LAPACK_WORK, LAPACK_WORK, 1, LAPACK_WORK, 1, SIZES(2:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 81; RETURN; END IF
    ! Query the size of the work array for DGESVD. (COMPUTE_ORTHOGONAL)
    CALL DGESVD('N', 'A', NUM_DVECS, DIM, LAPACK_WORK, NUM_DVECS, &
         LAPACK_WORK, LAPACK_WORK, 1, LAPACK_WORK, NUM_DVECS, &
         SIZES(1:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO + 81; RETURN; END IF
    ! Allocate the work array by rounding the max value in 'SIZES'.
    ALLOCATE(LAPACK_WORK(1:INT(.5_R8 + MAXVAL(SIZES))))

    ! Allocate various global variables.
    DEPTH = SUM(DVEC_MULTS) - DIM + 1 ! <- Max recursion depth.
    ALLOCATE(&
         SHIFTED_EVAL_PTS (NUM_DVECS, NUM_PTS,   DEPTH+1), &
         DVECS            (DIM,       NUM_DVECS, DEPTH+1), &
         EVALS_AT_PTS     (NUM_PTS,   DEPTH+1),            &
         LOC              (NUM_DVECS, DEPTH+1),            &
         MULTS            (NUM_DVECS, DEPTH),              &
         TEMP_EVAL_PTS    (DIM, NUM_PTS),                  &
         REMAINING_DVECS  (NUM_DVECS),                     &
         NONZERO_DVECS    (NUM_DVECS),                     &
         PT_SHIFT_R       (NUM_DVECS),                     &
         NORMAL_VECTOR    (DIM),                           &
         PT_SHIFT         (DIM),                           &
         )
  END SUBROUTINE RESERVE_MEMORY

  ! ==================================================================
  SUBROUTINE FREE_MEMORY()
    ! 9) FREE_MEMORY
    ! 
    !   Deallocate storage reserved for box-spline evaluation.
    ! 
    IF (ALLOCATED(LAPACK_WORK))      DEALLOCATE(LAPACK_WORK)
    IF (ALLOCATED(SHIFTED_EVAL_PTS)) DEALLOCATE(SHIFTED_EVAL_PTS)
    IF (ALLOCATED(DVECS))            DEALLOCATE(DVECS)
    IF (ALLOCATED(EVALS_AT_PTS))     DEALLOCATE(EVALS_AT_PTS)
    IF (ALLOCATED(LOC))              DEALLOCATE(LOC)
    IF (ALLOCATED(MULTS))            DEALLOCATE(MULTS)
    IF (ALLOCATED(TEMP_EVAL_PTS))    DEALLOCATE(TEMP_EVAL_PTS)
    IF (ALLOCATED(REMAINING_DVECS))  DEALLOCATE(REMAINING_DVECS)
    IF (ALLOCATED(NONZERO_DVECS))    DEALLOCATE(NONZERO_DVECS)
    IF (ALLOCATED(PT_SHIFT_R))       DEALLOCATE(PT_SHIFT_R)
    IF (ALLOCATED(NORMAL_VECTOR))    DEALLOCATE(NORMAL_VECTOR)
    IF (ALLOCATED(PT_SHIFT))         DEALLOCATE(PT_SHIFT)
  END SUBROUTINE FREE_MEMORY


END SUBROUTINE BOXSPLEV
