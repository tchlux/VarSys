! ==================================================================
FUNCTION MATRIX_RANK(MATRIX)
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
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
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: MATRIX
  ! Local variables for computing the orthogonal vector.
  REAL(KIND=REAL64), DIMENSION(SIZE(MATRIX,1),SIZE(MATRIX,2)) :: A
  REAL(KIND=REAL64), DIMENSION(MIN(SIZE(A,1),SIZE(A,2))) :: SING_VALS
  REAL(KIND=REAL64), DIMENSION(100) :: LAPACK_WORK

  INTEGER :: IDX, MATRIX_RANK, INFO
  ! Unused DGESVD parameters.
  REAL(KIND=REAL64), DIMENSION(1) :: U, VT
  REAL(KIND=REAL64) :: TOL

  A = MATRIX
  ! Use the SVD to get the orthogonal vectors.
  CALL DGESVD('N', 'N', SIZE(A,1), SIZE(A,2), A, SIZE(A,1), SING_VALS,&
       U, SIZE(U), VT, SIZE(VT), LAPACK_WORK, SIZE(LAPACK_WORK), INFO)
  IF (INFO .NE. 0) PRINT *, "ERROR:", INFO

  TOL = EPSILON(1._REAL64) * MAXVAL(SHAPE(A)) * MAXVAL(SING_VALS)
  PRINT *, "TOL:       ", TOL
  PRINT *, "SING_VALS: ", SING_VALS

  ! Find the first singular value in the orthonormal basis for the null
  ! space of A (a vector with an extremely small singular value)
  find_null : DO IDX = SIZE(SING_VALS), 1
     IF (SING_VALS(IDX) .GT. TOL) THEN
        MATRIX_RANK = IDX
        EXIT find_null
     ELSE
  END DO find_null
END FUNCTION MATRIX_RANK

