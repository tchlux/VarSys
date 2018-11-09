'''! This file (boxspline.f90) contains the subroutine BOXSPLEV that can
! be used to evaluate a box-spline, defined by its associated
! direction vector set, at given evaluation points, as well as the
! module REAL_PRECISION defining the real precision.
!
! ====================================================================
! MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
!   INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
! END MODULE REAL_PRECISION'''


import cython
import numpy

class NotFortranCompatible(Exception): pass

#      Wrapper for fortran function boxsplev     
# =================================================

cdef extern:
    void c_boxsplev( int* dvecs_0, int* dvecs_1, double* dvecs, int* dvec_mults_0, int* dvec_mults, int* eval_pts_0, int* eval_pts_1, double* eval_pts, int* box_evals_0, double* box_evals, int* error )

@cython.boundscheck(False)
@cython.wraparound(False)
def boxsplev( double[:,:] dvecs, int[:] dvec_mults, double[:,:] eval_pts, double[:] box_evals, error=None ):
    '''! BOXSPLEV evaluates a box-spline, defined by a direction vector set
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
    !'''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(dvecs).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int dvecs_0 = dvecs.shape[0]
    cdef int dvecs_1 = dvecs.shape[1]
    
    cdef int dvec_mults_0 = dvec_mults.shape[0]
    
    if (not numpy.asarray(eval_pts).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int eval_pts_0 = eval_pts.shape[0]
    cdef int eval_pts_1 = eval_pts.shape[1]
    
    cdef int box_evals_0 = box_evals.shape[0]
    
    if (type(error) == type(None)):
        error = 1
    cdef int local_error = error
    
    
    # Make fortran function call
    c_boxsplev(&dvecs_0, &dvecs_1, &dvecs[0][0], &dvec_mults_0, &dvec_mults[0], &eval_pts_0, &eval_pts_1, &eval_pts[0][0], &box_evals_0, &box_evals[0], &local_error)
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(box_evals, order='F'), local_error



