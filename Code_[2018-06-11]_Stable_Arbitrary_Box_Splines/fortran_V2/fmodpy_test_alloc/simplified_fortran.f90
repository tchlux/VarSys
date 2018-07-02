
FUNCTION ALLOCATE_MAX_LAPACK_WORK ( DIM , NUM_DVECS ) RESULT ( SIZE )
USE ISO_FORTRAN_ENV , ONLY : REAL64
! 8) ALLOCATE_MAX_LAPACK_WORK
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
INTEGER , INTENT ( IN ) : : DIM , NUM_DVECS
INTEGER : : SIZE , INFO
!
REAL ( KIND = REAL64 ) , DIMENSION ( : ) , ALLOCATABLE : : WORK
REAL , DIMENSION ( 3 ) : : SIZES
! Initialize all sizes to 1
! Retrieve expected size from each of the LAPACK routines.

! Query the size of the work array for DGESVD. (MATRIX_ORTHOGONAL)

! Query the size of the work array to (MATRIX_RANK)

! Get the size of the work array for DGELS. (MATRIX_MINIMUM_NORM)

END FUNCTION ALLOCATE_MAX_LAPACK_WORK

