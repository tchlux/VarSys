! ==================================================================
SUBROUTINE MATRIX_MINIMUM_NORM(MATRIX)
  USE ISO_FORTRAN_ENV, ONLY: REAL64
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
  REAL(KIND=REAL64), DIMENSION(300) :: LAPACK_WORK
  ! Local variables
  INTEGER :: IDX, INFO
  ! Allocate the output matrix.
  ALLOCATE(MIN_NORM(1:SIZE(MATRIX,2),1:SIZE(MATRIX,2)))
  ! Make "MIN_NORM" the identity matrix
  MIN_NORM = 0.
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

  ! Extract the minimum norm representation from the output of DGELS.
  PRINT *, "MIN_NORM: "
  DO INFO = 1, SIZE(MIN_NORM,1)
     PRINT *, MIN_NORM(INFO,:)
  END DO
  INFO = 0
  MIN_NORM = MIN_NORM(1:SIZE(MATRIX,1),:)
  PRINT *, "MIN_NORM: "
  DO INFO = 1, SIZE(MIN_NORM,1)
     PRINT *, MIN_NORM(INFO,:)
  END DO
END SUBROUTINE MATRIX_MINIMUM_NORM
