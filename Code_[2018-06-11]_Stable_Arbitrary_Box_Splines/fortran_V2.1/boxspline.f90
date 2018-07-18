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
  !      1  Error computing normal vectors.
  !      2  Error computing box-spline.
  !      3  Error computing an orthogonal vector.
  !      4  Error computing a minimum norm representation.
  !      5  Error computing a matrix rank.
  !      6  Error computing a matrix determinant.
  !      7  Error computing the reciprocal condition number of matrix.
  !      8  Error finding nonzero entries in an array.
  !      9  Error preparing memory.
  ! 
  !    Units digit U, for T = 0:
  !      1  Mismatched dimension, SIZE(DVEC_MULTS) .NE. SIZE(DVECS,2).
  !      2  Mismatched dimension, SIZE(EVAL_PTS,1) .NE. SIZE(DVECS,1).
  !      3  Mismatched dimension, SIZE(BOX_EVALS)  .NE. SIZE(EVAL_PTS,2).
  !      4  One of the multiplicity values provided was < 1.
  !      5  Columns of DVECS are not unique.
  !      6  M < S or DVECS(1:S,1:S) is near singular.
  ! 
  !    Units digit U, for T /= 0:
  !      0  Work array allocation failed.
  !      1  Work array size query failed.
  !      2  DGELS computation error, see 'INFO' for details.
  !      3  DGESVD computation error, see 'INFO' for details.
  !      4  DGETRF computation error, see 'INFO' for details.
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
  ! USE REAL_PRECISION, ONLY: R8
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: DVECS
  INTEGER,           INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
  REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(EVAL_PTS,1)) :: BOX_EVALS
  INTEGER,           INTENT(OUT)                 :: ERROR
  ! Local variables.
  INTEGER :: DIM, NUM_DVECS, NUM_PTS, INFO, IDX_1, IDX_2
  INTEGER,           DIMENSION(SIZE(DVECS,2))                  :: LOCATION, LOOKUP
  INTEGER,           DIMENSION(SIZE(DVECS,2),SIZE(DVECS,2))    :: IDENTITY
  REAL(KIND=REAL64), DIMENSION(SIZE(DVECS,1),SIZE(DVECS,2))    :: DVECS_COPY
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE               :: MIN_NORM_DVECS
  REAL(KIND=REAL64), DIMENSION(:),   ALLOCATABLE               :: LAPACK_WORK
  REAL(KIND=REAL64), PARAMETER :: SQRTEPS = SQRT(EPSILON(1.0_REAL64))

  ERROR = 0; INFO = 0; ! Initialize error and info flags.
  ! Store 'global' constants for box-spline evaluation.
  DIM = SIZE(DVECS, 1)
  NUM_DVECS = SIZE(DVECS, 2)
  NUM_PTS = SIZE(EVAL_PTS, 1)
  ! Allocate a work array large enough for all LAPACK calls.
  LAPACK_WORK = ALLOCATE_MAX_LAPACK_WORK()
  IF (ERROR .NE. 0) RETURN
  ! Create an identity matrix (for easy matrix-vector multiplication)
  ! and compute the index lookup indices (for binary style hashing).
  IDENTITY = 0
  FORALL (IDX_1 = 1:NUM_DVECS)
     IDENTITY(IDX_1,IDX_1) = 1
     LOOKUP(IDX_1) = 2**(IDX_1-1)
  END FORALL
  ! Check for usage errors.
  ! Check for dimension mismatches.
  IF (SIZE(DVEC_MULTS) .NE. NUM_DVECS) THEN; ERROR = 1; RETURN; END IF
  IF (SIZE(EVAL_PTS,2) .NE. DIM)       THEN; ERROR = 2; RETURN; END IF
  IF (SIZE(BOX_EVALS)  .NE. NUM_PTS)   THEN; ERROR = 3; RETURN; END IF
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
  IF (NUM_DVECS .LT. DIM) THEN 
     ERROR = 6; RETURN
  ELSE
     ! Compute condition number COND(DVECS(1:DIM,1:DIM)) and test for near
     ! rank deficiency: 1/COND(DVECS(1:DIM,1:DIM)) < SQRT(EPSILON(1.0_REAL64)).
     IF (MATRIX_CONDITION_INV(DVECS(1:DIM, 1:DIM)) .LT. SQRTEPS) THEN
        ERROR = 6; RETURN
     ENDIF
  ENDIF

  ! Get the minimum norm representation of the direction vectors, use
  ! a copy to ensure the original DVECS are not modified.
  DVECS_COPY(:,:) = DVECS(:,:)
  MIN_NORM_DVECS = MATRIX_MINIMUM_NORM(DVECS_COPY)
  IF (ERROR .NE. 0) RETURN
  
  ! Recursively evaluate the box-spline.
  LOCATION(:) = 0
  BOX_EVALS(:) = 0_REAL64
  CALL EVALUATE_BOX_SPLINE(DVEC_MULTS, LOCATION, MIN_NORM_DVECS, &
       MATMUL(EVAL_PTS, MIN_NORM_DVECS), BOX_EVALS)
  IF (ERROR .NE. 0) RETURN

CONTAINS
  
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
    INTEGER,           INTENT(IN),  DIMENSION(SIZE(DVECS,2))    :: MULTS, LOC
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)              :: SUB_DVECS
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:)              :: SHIFTED_EVAL_PTS
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(EVAL_PTS,1)) :: EVALS_AT_PTS
    ! Temporary variables for recursive computations of eval points,
    ! shifted evaluation points, and the actual evaluations of box splines.
    REAL(KIND=REAL64), DIMENSION(SIZE(SHIFTED_EVAL_PTS,1),&
         SIZE(SHIFTED_EVAL_PTS,2))                 :: TEMP_SHIFTED_EVAL_PTS
    REAL(KIND=REAL64), DIMENSION(SIZE(EVAL_PTS,1), &
         SIZE(EVAL_PTS,2))                         :: TEMP_EVAL_PTS
    REAL(KIND=REAL64), DIMENSION(SIZE(EVAL_PTS,1)) :: TEMP_EVALS_AT_PTS
    REAL(KIND=REAL64), DIMENSION(SIZE(EVAL_PTS,1)) :: LOCATIONS
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: NEXT_DVECS
    INTEGER,           DIMENSION(:),   ALLOCATABLE :: REMAINING_DVECS
    INTEGER,           DIMENSION(SIZE(DVECS,2))    :: NEXT_MULTS, NEXT_LOC
    REAL(KIND=REAL64), DIMENSION(SIZE(DVECS,1))    :: PT_SHIFT, NORMAL_VECTOR
    REAL(KIND=REAL64) :: POSITION
    INTEGER :: IDX_1, IDX_2, IDX_3
    ! Recursion case ...
    IF (SUM(MULTS) > DIM) THEN
       EVALS_AT_PTS = 0.
       IDX_2 = 1
       ! Sum over all direction vectors.
       DO IDX_1 = 1, NUM_DVECS
          ! Update multiplicity of directions and position in recursion tree.
          NEXT_MULTS(:) = MULTS(:) - IDENTITY(:, IDX_1) ! Reduce multiplicity.
          NEXT_LOC(:)   = LOC(:)   + IDENTITY(:, IDX_1) ! Track recursion position.
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
                TEMP_SHIFTED_EVAL_PTS(:,IDX_3) = SHIFTED_EVAL_PTS(:,IDX_3) - &
                     SUM(DVECS(:,IDX_1) * SUB_DVECS(:,IDX_3))
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
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - SUM(LOC * DVECS(IDX_3,:))
                END FORALL compute_shift_2
                CALL EVALUATE_BOX_SPLINE(NEXT_MULTS, LOC, NEXT_DVECS,&
                     MATMUL(TEMP_EVAL_PTS, NEXT_DVECS), TEMP_EVALS_AT_PTS)
                IF (ERROR .NE. 0) RETURN
                EVALS_AT_PTS = EVALS_AT_PTS + TEMP_EVALS_AT_PTS * &
                     SHIFTED_EVAL_PTS(:,IDX_2)
                ! Perform recursion with transformed set of direction vectors.
                compute_shift_3 : FORALL (IDX_3 = 1:DIM)
                   TEMP_EVAL_PTS(:,IDX_3) = EVAL_PTS(:,IDX_3) - SUM(NEXT_LOC * DVECS(IDX_3,:))
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
       EVALS_AT_PTS = EVALS_AT_PTS / REAL(SUM(MULTS) - DIM, REAL64)
    ELSE
       ! Base case ... compute characteristic function.
       EVALS_AT_PTS = 1.
       ! Delayed translations (this is what makes the algorithm more stable).
       REMAINING_DVECS = NONZERO(MULTS)
       IF (ERROR .NE. 0) RETURN
       compute_shift_4 : FORALL (IDX_1 = 1:DIM)
          TEMP_EVAL_PTS(:,IDX_1) = EVAL_PTS(:,IDX_1) - SUM(LOC * DVECS(IDX_1,:))
       END FORALL compute_shift_4
       ! Check against all *precomputed* hyperplanes (also contributes to stability).
       DO IDX_1 = 1, DIM
          ! Calculate normal vector to current hyperplane.
          CALL MATRIX_ORTHOGONAL(TRANSPOSE(DVECS(:, &
               NONZERO(MULTS - IDENTITY(:,REMAINING_DVECS(IDX_1))))), &
               NORMAL_VECTOR)
          ! Compute shifted position (relative to normal ector).
          POSITION = SUM(DVECS(:,REMAINING_DVECS(IDX_1)) * NORMAL_VECTOR(:))
          ! Compute shifted evaluation locations. (NUM_PTS, DIM) x (DIM, 1)
          compute_shifted_point_1 : FORALL (IDX_2 = 1:NUM_PTS)
             LOCATIONS(IDX_2) = SUM(TEMP_EVAL_PTS(IDX_2,:) * NORMAL_VECTOR(:))
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
          compute_shifted_point_2 : FORALL (IDX_2 = 1:NUM_PTS)
             LOCATIONS(IDX_2) = SUM((EVAL_PTS(IDX_2,:) - PT_SHIFT) * &
                  NORMAL_VECTOR(:))
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
    !   If there are any singular values [value < SQRT(epsilon)],
    !   the vector associated with the largest such singular values is
    !   returned.
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
    INTEGER :: IDX
    LOGICAL :: FOUND_ZERO
    REAL(KIND=REAL64) :: TOL
    ! Unused parameters
    REAL(KIND=REAL64), DIMENSION(1) :: U

    ! Use the SVD to get the orthogonal vectors.
    CALL DGESVD('N','A',SIZE(A,1),SIZE(A,2),A,SIZE(A,1),S,U,SIZE(U), &
         VT, SIZE(VT,1), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
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

    ORTHOGONAL(:) = 0.
    FOUND_ZERO = .FALSE.
    TOL = EPSILON(0._REAL64) * SIZE(S) ! <- Compute tolerance for considering a singular value '0'
    ! Find the first vector in the orthonormal basis for the null
    ! space of A (a vector with an extremely small singular value).
    find_null : DO IDX = 1, SIZE(S)
       IF (S(IDX) .LE. TOL) THEN
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
    !   solution to the problem (AX = I). This is a more numerically
    !   stable solution to the linear system (A^T A) X = A^T.
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   MIN_NORM -- Real dense matrix that is the minimum norm
    !               representation of MATRIX.
    ! 
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: MIN_NORM
    ! Local variables
    INTEGER :: IDX

    ! Allocate the output matrix.
    ALLOCATE(MIN_NORM(1:SIZE(MATRIX,2),1:SIZE(MATRIX,2)))
    ! Make "MIN_NORM" the identity matrix
    MIN_NORM(:,:) = 0.
    FORALL (IDX = 1:SIZE(MIN_NORM,1)) MIN_NORM(IDX,IDX) = 1.

    ! Call DGELS for actual solve.
    CALL DGELS('T', SIZE(MATRIX,1), SIZE(MATRIX,2), SIZE(MIN_NORM,1),&
         MATRIX, SIZE(MATRIX,1), MIN_NORM, SIZE(MIN_NORM,1), &
         LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
    !   'T'               -- Transpose of A is necessary
    !   SIZE(MATRIX,1)    -- number of rows in A
    !   SIZE(MATRIX,2)    -- number of columns in A
    !   SIZE(MIN_NORM,2)  -- Number of right hand sides
    !   MATRIX            -- A
    !   SIZE(MATRIX,1)    -- Leading dimension of A (number of rows)
    !   MIN_NORM          -- B (will be overwritten to hold X after call)
    !   SIZE(MIN_NORM,1)  -- Leading dimension of B (number of rows)
    !   LAPACK_WORK       -- Workspace array for DGELS
    !   SIZE(LAPACK_WORK) -- Size of dgels_work_array
    !   INFO              -- For verifying successful execution

    IF (INFO .NE. 0) THEN; ERROR = 43 + 100*INFO; RETURN; END IF
    ! Extract the minimum norm representation from the output of DGELS.
    MIN_NORM = MIN_NORM(1:SIZE(MATRIX,1),:)
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
    REAL(KIND=REAL64), INTENT(IN),  DIMENSION(:,:) :: MATRIX
    ! Local variables for computing the orthogonal vector.
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: A
    REAL(KIND=REAL64), DIMENSION(MIN(SIZE(MATRIX,1),SIZE(MATRIX,2))) :: S
    INTEGER :: IDX, MATRIX_RANK
    ! Unused DGESVD parameters.
    REAL(KIND=REAL64), DIMENSION(1) :: U, VT

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
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=REAL64) :: MATRIX_DET
    ! Local variable
    REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: M
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
  FUNCTION MATRIX_CONDITION_INV(MATRIX) RESULT(RCOND)
    ! 7) MATRIX_CONDITION_INV
    ! 
    ! Compute the condition number (for testing near rank deficiency)
    ! using DGECON (which computes 1 / CONDITION).
    ! 
    ! Input:
    !   MATRIX -- Real dense matrix.
    ! 
    ! Output:
    !   RCOND -- Real value corresponding to the reciprocal of the
    !            condition number of the provided matrix.
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: MATRIX
    REAL(KIND=REAL64) :: RCOND
    ! Local arrays for work.
    REAL(KIND=REAL64), DIMENSION(4*SIZE(MATRIX,2)) :: WORK
    INTEGER,       DIMENSION(  SIZE(MATRIX,2)) :: IWORK
    ! Use LAPACK
    CALL DGECON('I', SIZE(MATRIX,2), MATRIX, SIZE(MATRIX,1),&
         SUM(ABS(MATRIX)), RCOND, WORK, IWORK, INFO)
    ! Store output of 'info' flag.
    ERROR = 100 * INFO
  END FUNCTION MATRIX_CONDITION_INV

  ! ==================================================================
  FUNCTION NONZERO(ARRAY) RESULT(NE_ZERO)
    ! 8) NONZERO
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

  ! ==================================================================
  FUNCTION ALLOCATE_MAX_LAPACK_WORK() RESULT(WORK)
    ! 9) ALLOCATE_MAX_LAPACK_WORK
    ! 
    ! Return an allocated real array that has a size equal to the
    ! maximum requested LAPACK work array size across all routines
    ! that will be executed in the evaluation of a Box Spline.
    ! 
    ! Output:
    !   WORK -- Real array with size large enough to accomadate all
    !           LAPACK subroutines that will be used to evaluated a
    !           Box Spline given the dimension and number of vectors.
    ! 
    REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: WORK
    REAL(KIND=REAL64), DIMENSION(3) :: SIZES
    ! Initialize all sizes to 1
    SIZES = 0.
    ! Retrieve expected size from each of the LAPACK routines.

    ! Query the size of the work array for DGESVD. (MATRIX_ORTHOGONAL)
    CALL DGESVD('N', 'A', NUM_DVECS, DIM, WORK, NUM_DVECS, &
         WORK, WORK, 1, WORK, NUM_DVECS, SIZES(1:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO; RETURN; END IF

    ! Query the size of the work array to (MATRIX_RANK)
    CALL DGESVD('N', 'N', NUM_DVECS, DIM, WORK, NUM_DVECS, &
         WORK, WORK, 1, WORK, 1, SIZES(2:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO; RETURN; END IF

    ! Get the size of the work array for DGELS. (MATRIX_MINIMUM_NORM)
    CALL DGELS('T', DIM, NUM_DVECS, NUM_DVECS, WORK, DIM, WORK, &
         NUM_DVECS, SIZES(3:), -1, INFO)
    IF (INFO .NE. 0) THEN; ERROR = 100*INFO; RETURN; END IF

    ! Allocate the work array by rounding the max value in 'SIZES'.
    ALLOCATE(WORK(1:INT(.5 + MAXVAL(SIZES))))
  END FUNCTION ALLOCATE_MAX_LAPACK_WORK

END SUBROUTINE BOXSPLEV

    
! python3 -c "import fmodpy; fmodpy.wrap('boxspline.f90', module_link_args=['-lblas','-llapack','-lgfortran'], verbose=True)"
