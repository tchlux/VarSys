! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!                          sample_main.f90
! 
! DESCRIPTION:
!   This file (sample_main.f90) contains a sample main program that calls
!   MQSI to interpolate handcrafted nonmonotone data. Compile with:
! 
!    $F03 $OPTS REAL_PRECISION.f90 EVAL_BSPLINE.f90 SPLINE.f90 MQSI.f90 \
!       sample_main.f90 -o main $LIB
!
!   where '$F03' is the name of the Fortran 2003 compiler, '$OPTS' are
!   compiler options such as '-O3', and '$LIB' provides a flag to link
!   BLAS and LAPACK. If the BLAS and LAPACK libraries are not
!   available on your system, then replace $LIB with the filenames
!   'blas.f lapack.f'; these files contain the routines from the BLAS
!   and LAPACK libraries that are necessary for this package.
!
! CONTAINS:
!   PROGRAM SAMPLE_MAIN
! 
! DEPENDENCIES:
!   MODULE REAL_PRECISION
!     INTEGER, PARAMETER :: R8
!   END MODULE REAL_PRECISION
! 
!   SUBROUTINE MQSI(X, Y, T, BCOEF, INFO)
!     USE REAL_PRECISION, ONLY: R8
!     REAL(KIND=R8), INTENT(IN),    DIMENSION(:) :: X
!     REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: Y
!     REAL(KIND=R8), INTENT(OUT),   DIMENSION(:) :: T, BCOEF
!     INTEGER, INTENT(OUT) :: INFO
!   END SUBROUTINE MQSI
! 
!   SUBROUTINE EVAL_SPLINE(T, BCOEF, XY, INFO, D)
!     USE REAL_PRECISION, ONLY: R8
!     REAL(KIND=R8), INTENT(IN),    DIMENSION(:) :: T, BCOEF
!     REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
!     INTEGER, INTENT(OUT) :: INFO
!     INTEGER, INTENT(IN), OPTIONAL :: D
!   END SUBROUTINE EVAL_SPLINE
! 
! CONTRIBUTORS:
!   Thomas C.H. Lux (tchlux@vt.edu)
!   Layne T. Watson (ltwatson@computer.org)
!   William I. Thacker (thackerw@winthrop.edu)
! 
! VERSION HISTORY:
!   June 2020 -- (tchl) Created file, (ltw / wit) reviewed and revised.
! 
PROGRAM SAMPLE_MAIN
USE REAL_PRECISION, ONLY: R8  
IMPLICIT NONE

! Define the interfaces for relevant MQSI package subroutines.
INTERFACE
 SUBROUTINE MQSI(X, Y, T, BCOEF, INFO)
   USE REAL_PRECISION, ONLY: R8
   REAL(KIND=R8), INTENT(IN),  DIMENSION(:) :: X
   REAL(KIND=R8), INTENT(INOUT),  DIMENSION(:) :: Y
   REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: T, BCOEF
   INTEGER, INTENT(OUT) :: INFO
 END SUBROUTINE MQSI
 SUBROUTINE EVAL_SPLINE(T, BCOEF, XY, INFO, D)
   USE REAL_PRECISION, ONLY: R8
   REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: T, BCOEF
   REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: XY
   INTEGER, INTENT(OUT) :: INFO
   INTEGER, INTENT(IN), OPTIONAL :: D
 END SUBROUTINE EVAL_SPLINE
END INTERFACE

! Definitions of data.
INTEGER, PARAMETER :: ND = 18  ! Define the size of the problem.

!Define the X and Y for interpolation.
REAL(KIND=R8), DIMENSION(ND) :: X
REAL(KIND=R8), DIMENSION(ND) :: &
     Y = (/ 0.0_R8, 1.0_R8, 1.0_R8, 1.0_R8, 0.0_R8, 20.0_R8, &
            19.0_R8, 18.0_R8, 17.0_R8, 0.0_R8, 0.0_R8, 3.0_R8, &
            0.0_R8, 1.0_R8, 6.0_R8, 16.0_R8, 16.1_R8, 1.0_R8 /)
! Definition of variables to hold the spline coefficients SC, knots 
! SK, and temporary data storage U.
REAL(KIND=R8) :: SC(1:3*ND), SK(1:3*ND+6), U(ND)
! Loop iterator I and integer INFO for checking execution status.
INTEGER :: I, INFO

! Initialize "X" values.
DO I = 1, ND 
   X(I) = REAL(I, KIND=R8) 
END DO
! Construct a MQSI.
CALL MQSI(X, Y, SK, SC, INFO)
IF (INFO .NE. 0) WRITE(*,101) INFO
101 FORMAT('MQSI returned info =',I4)
! Copy the "evaluation points" into the temporary variable U to
! save the X values since that argument is overwritten by EVAL_SPLINE.
U(:) = X(:)
! Evaluate the spline at all points (result is updated in-place in U).
CALL EVAL_SPLINE(SK, SC, U, INFO)
IF (INFO .NE. 0) WRITE(*,111) INFO
111 FORMAT(/,'EVAL_SPLINE returned info= ',I4)
! Show the evaluations.
WRITE (*,121)
121 FORMAT('     I', 5X, 'X(I)', 5X, 'Y(I)', 3X, 'Spline Value')
DO I=1,ND
   WRITE(*,131) I, X(I), Y(I), U(I)
END DO
131 FORMAT( I6, 3F9.2)

WRITE (*,141)  EPSILON(1.0_R8)
141 FORMAT (/,/,'Machine precision:  ',F18.16)
WRITE(*,151) MAXVAL(ABS(U(:) - Y(:)) / ABS(1.0_R8 + ABS(Y(:))))
151 FORMAT(/,'Max relative error: ',F18.16)
! Evaluate the first derivative of the spline at all points.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=1)
IF (INFO .NE. 0) WRITE(*,161) INFO
161 FORMAT(/,'EVAL_SPLINE returned info ',I4)

WRITE(*,171)
171 FORMAT(/,'   I   First Derivative')
DO I=1,ND
   WRITE(*,181) I,U(I)
END DO
181 FORMAT(I4, 2X, F9.2)

! Evaluate the second derivative of the spline at all points.
U(:) = X(:)
CALL EVAL_SPLINE(SK, SC, U, INFO, D=2)
IF (INFO .NE. 0) WRITE(*,191) INFO
191 FORMAT (/,'EVAL_SPLINE returned info =',I4)
WRITE(*,201)
201 FORMAT(/,'   I   Second Derivative')
DO I=1,ND
   WRITE(*,181) I,U(I)
END DO

END PROGRAM SAMPLE_MAIN
  
