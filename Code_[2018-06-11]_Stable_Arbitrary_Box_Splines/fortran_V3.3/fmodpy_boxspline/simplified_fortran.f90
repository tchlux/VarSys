
! This file (boxspline.f90) contains the subroutine BOXSPLEV that can
! be used to evaluate a box-spline, defined by its associated
! direction vector set, at given evaluation points, as well as the
! module REAL_PRECISION defining the real precision.
!
! ====================================================================
! MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
!   INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
! END MODULE REAL_PRECISION

SUBROUTINE BOXSPLEV ( DVECS , DVEC_MULTS , EVAL_PTS , BOX_EVALS , ERROR )
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
USE REAL_PRECISION , ONLY : R8
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : DVECS
INTEGER , INTENT ( IN ) , DIMENSION ( : ) : : DVEC_MULTS
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : EVAL_PTS
REAL ( KIND = R8 ) , INTENT ( OUT ) , DIMENSION ( : ) : : BOX_EVALS
INTEGER , INTENT ( OUT ) : : ERROR
! Reusable recursion variables and global variables.
INTEGER : : DIM , NUM_DVECS , NUM_PTS , INFO
REAL ( KIND = R8 ) : : POSITION
REAL ( KIND = R8 ) , PARAMETER : : SQRTEPS = SQRT ( EPSILON ( 1.0_R8 ) )
REAL ( KIND = R8 ) , PARAMETER : : ONE = 1.0_R8
REAL ( KIND = R8 ) , PARAMETER : : ZERO = 0.0_R8
! ------------------------------------------------------------------
INTEGER , DIMENSION ( SIZE ( DVECS , 2 ) , SIZE ( DVECS , 2 ) ) : : IDENTITY
INTEGER , DIMENSION ( SIZE ( DVECS , 1 ) , SIZE ( DVECS , 1 ) ) : : DET_WORK
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 2 ) , SIZE ( DVECS , 2 ) ) : : ORTHO_MIN_WORK
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 1 ) - 1 , SIZE ( DVECS , 1 ) ) : : TRANS_DVECS
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 1 ) , SIZE ( EVAL_PTS , 2 ) ) : : TEMP_PTS
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 1 ) ) : : PT_SHIFT
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 1 ) ) : : SING_VALS
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 1 ) ) : : NORMAL_VECTOR
INTEGER , DIMENSION ( SIZE ( DVECS , 2 ) ) : : REMAINING_DVECS
INTEGER , DIMENSION ( SIZE ( DVECS , 2 ) ) : : NONZERO_DVECS
REAL ( KIND = R8 ) , DIMENSION ( : ) , ALLOCATABLE : : LAPACK_WORK
! Unique recursion variables for evaluating box-spline. Last dimension is depth.
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 1 ) , SIZE ( DVECS , 2 ) ) : : MIN_NORM_DVECS
INTEGER , DIMENSION ( SIZE ( DVECS , 2 ) ) : : LOC
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( DVECS , 2 ) , SIZE ( EVAL_PTS , 2 ) ) : : SHIFTED_PTS
! ------------------------------------------------------------------
! Local variables for preparation.
INTEGER : : IDX_1 , IDX_2
! Initialize error and info flags.
! Initialize parameters
! Store 'global' constants for box-spline evaluation.
! Check for usage errors (dimension mismatches, invalid multiplicity)
! Check uniqueness of DVECS columns.
! Check if NUM_DVECS < DIM or rank DVECS(1:DIM,1:DIM) < DIM.
! Compute condition number COND(DVECS(1:DIM,1:DIM)) and test for near
! rank deficiency: 1/COND(DVECS(1:DIM,1:DIM)) < SQRT(EPSILON(1.0_R8)).
! Allocate a work array large enough for all LAPACK calls.
! Create an identity matrix (for easy matrix-vector multiplication).
! Calculate the DVECS for the beginning of box-spline evaluation.
! Compute the shifted evaluation points.
! Recursively evaluate the box-spline.

! ==================================================================
RECURSIVE SUBROUTINE EVALUATE_BOX_SPLINE ( C_DVECS , C_LOC , C_MULTS , C_SHIFTED_PTS , C_EVALS )
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
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : C_DVECS
INTEGER , INTENT ( IN ) , DIMENSION ( : ) : : C_LOC , C_MULTS
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : C_SHIFTED_PTS
REAL ( KIND = R8 ) , INTENT ( OUT ) , DIMENSION ( : ) : : C_EVALS
INTEGER , DIMENSION ( SIZE ( C_LOC ) ) : : N_LOC , N_MULTS
REAL ( KIND = R8 ) , DIMENSION ( SIZE ( C_EVALS ) ) : : N_EVALS
REAL ( KIND = R8 ) , DIMENSION ( : , : ) , ALLOCATABLE : : N_DVECS
! Local variables
INTEGER : : IDX_1 , IDX_2 , IDX_3
! Recursion case ...
! Sum over all direction vectors.
! Update multiplicity of directions and position in recursion tree.
! Make recursive calls.
! Copy LOC, DVECS, and SHIFTED_PTS down one level for next recursion.
! Perform recursion with only reduced multiplicity.
! Perform recursion with transformed set of direction
! vectors and evaluation points.
! Update location for next level of recursion.
! Get the remaining direction vectors.
! Update Least norm representation of the direction vectors.
! Perform recursion with only reduced multiplicity.
! Compute the appropriate point shift according the new DVECS.
! Perform recursion with only reduced multiplicity.
! Update location for next level of recursion.
! Perform recursion with transformed set of direction vectors.
! Normalize by number of direction vectors in computation.
! Base case ... compute characteristic function.
! Pack the unique direction vectors into the memory location
! for the current set of direction vectors (since the 'current'
! are not needed for base case evaluation).
! Delayed translations (this is what makes the algorithm more stable).
! Check evaluation point locations against all remaining direction vectors.
! Get the active set of direction vectors (exluding current vector).
! Calculate the orthogonal vector to the remaining direction vectors.
! Compute shifted position (relative to normal vector).
! Compute shifted evaluation locations. (1, DIM) x (DIM, NUM_PTS)
! Identify those points that are outside of this box (0-side).
! Recompute shifted location (other side of box) based on selected direction vector.
! Identify those shifted points that are outside of this box on REMAINING_DVEC(IDX_1)-side.
! Normalize evaluations by determinant of box.
END SUBROUTINE EVALUATE_BOX_SPLINE

!===================================================================
!             Supporting code for computing box-spline
!===================================================================

! ==================================================================
SUBROUTINE MAKE_DVECS_MIN_NORM ( DVECS )
! 2) MAKE_DVECS_MIN_NORM
!
!   Compute the minimum norm representation of 'MATRIX' and store
!   it in 'MIN_NORM', use DGELS to find the least squares
!   solution to the problem (AX = I). This is a more numerically
!   stable solution to the linear system (A^T A) X = A^T.
REAL ( KIND = R8 ) , INTENT ( INOUT ) , DIMENSION ( : , : ) : : DVECS
INTEGER : : IDX
! Make "ORTHO_MIN_WORK" the identity matrix
! Call DGELS for actual solve.
! Check for error.
! Extract the minimum norm representation from the output of DGELS.
END SUBROUTINE MAKE_DVECS_MIN_NORM

! ==================================================================
FUNCTION PACK_DVECS ( MULTS , SRCE_DVECS ) RESULT ( DEST_DVECS )
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
INTEGER , INTENT ( IN ) , DIMENSION ( : ) : : MULTS
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : SRCE_DVECS
REAL ( KIND = R8 ) , DIMENSION ( : , : ) , ALLOCATABLE : : DEST_DVECS
INTEGER : : IDX , NUM_DVECS
! Loop through once to find out how many DVECS there will be.
! Allocate the new direction vectors and copy them over.
END SUBROUTINE PACK_DVECS

! ==================================================================
SUBROUTINE COMPUTE_ORTHOGONAL ( A , ORTHOGONAL )
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
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : A
REAL ( KIND = R8 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( A , 2 ) ) : : ORTHOGONAL
! Local variables.
REAL ( KIND = R8 ) : : TOL
REAL ( KIND = R8 ) , DIMENSION ( 1 ) : : U
! Use the SVD to get the orthogonal vector(s).
! Compute the appropriate tolerance based on calculations.
! Check for a vector in the orthonormal basis for the null
! space of A (a vector with an extremely small singular value).
! If no orthogonal vector was found and the matrix does not contain
! enough vectors to span the space, use the first vector in VT at
! (RANK(A) + 1).
END SUBROUTINE COMPUTE_ORTHOGONAL

! ==================================================================
FUNCTION MATRIX_CONDITION_INV ( MATRIX ) RESULT ( RCOND )
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
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : MATRIX
REAL ( KIND = R8 ) : : RCOND
! Local arrays for work.
REAL ( KIND = R8 ) , DIMENSION ( 4 * SIZE ( MATRIX , 2 ) ) : : WORK
INTEGER , DIMENSION ( SIZE ( MATRIX , 2 ) ) : : IWORK
! Use LAPACK to compue recirpocal of matrix codition number.
! Check for errors during exeuction.
END FUNCTION MATRIX_CONDITION_INV

! ==================================================================
FUNCTION FULL_RANK ( MATRIX )
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
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : MATRIX
LOGICAL : : FULL_RANK
REAL ( KIND = R8 ) : : TOL
! Early return if possible, when we know the rank must be low.
! Use the SVD to get the orthogonal vectors.
! Compute a reasonable singular value tolerance based on expected
! numerical error ["Numerical Recipes" by W. H. Press et al., Matlab].
END FUNCTION FULL_RANK

! ==================================================================
FUNCTION MATRIX_DET ( MATRIX )
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
REAL ( KIND = R8 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : MATRIX
REAL ( KIND = R8 ) : : MATRIX_DET
INTEGER : : IDX
! Do the LU decomposition.
! Check for errors.
! Compute the determinant (product of diagonal of U).
END FUNCTION MATRIX_DET

! ==================================================================
SUBROUTINE PRINT_OPTIMAL_BLOCK_SIZE ( )
! 8) PRINT_OPTIMAL_BLOCK_SIZE
!
!   Use the LAPACK routine ILAENV to print out the optimal block
!   size for the current architecture (requires optimal LAPACK and
!   BLAS installation).
!
INTEGER : : BLOCK_SIZE
INTERFACE
FUNCTION ILAENV ( ISPEC , NAME , OPTS , N1 , N2 , N3 , N4 )
INTEGER : : ISPEC
CHARACTER : : NAME , OPTS
INTEGER , OPTIONAL : : N1 , N2 , N3 , N4
INTEGER : : ILAENV
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
END SUBROUTINE PRINT_OPTIMAL_BLOCK_SIZE

END SUBROUTINE BOXSPLEV

