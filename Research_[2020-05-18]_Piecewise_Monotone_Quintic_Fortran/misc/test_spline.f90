MODULE FS_TIME
  INTEGER, PARAMETER :: MAX_LAPS = 100
  REAL :: TIMES(1:MAX_LAPS)
  REAL :: START_TIME_SEC, STOP_TIME_SEC
  INTEGER :: FS_I

 CONTAINS
   SUBROUTINE TIME_RESET()
     ! Reset all recorded times.
     FS_I = 1
     TIMES(:) = 0.0
   END SUBROUTINE TIME_RESET

   SUBROUTINE START()
     ! Initialize a start time.
     CALL CPU_TIME(START_TIME_SEC)
   END SUBROUTINE START

   SUBROUTINE LAP()
     ! Measure the time since START_TIME_SEC, store, and print.
     CALL CPU_TIME(STOP_TIME_SEC)
     ! Give error for unexpected usage, handle elegantly.
     IF (FS_I .EQ. MAX_LAPS+1) THEN
        WRITE (*,'(A,/,A)') &
             'WARNING (FS_TIME): Called LAP too many times, record stops here.', &
             ' Increase MAX_LAPS parameter to store more lap times.'
     ELSE
        TIMES(FS_I) = STOP_TIME_SEC - START_TIME_SEC
     END IF
     WRITE (*,'(I3,A,F10.6)') FS_I, '  at  ', TIMES(FS_I)
     FS_I = FS_I + 1
     START_TIME_SEC = STOP_TIME_SEC
   END SUBROUTINE LAP
END MODULE FS_TIME


PROGRAM TIME_SPLINE
  USE REAL_PRECISION, ONLY: R8
  USE FS_TIME, ONLY: TIME_RESET

  ! Program constants.
  INTEGER, PARAMETER :: NB_MIN = 100
  INTEGER, PARAMETER :: NB_MAX = 1000
  INTEGER, PARAMETER :: NB_STEP = 100
  INTEGER, PARAMETER :: NCC = 3
  
  ! FIT_SPLINE routine interface.
  INTERFACE
     SUBROUTINE FIT_SPLINE(XI, FX, T, BCOEF, INFO)
       USE REAL_PRECISION, ONLY: R8
       REAL(KIND=R8), INTENT(IN),  DIMENSION(:)   :: XI
       REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: FX
       REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: T, BCOEF
       INTEGER, INTENT(OUT) :: INFO
     END SUBROUTINE FIT_SPLINE
  END INTERFACE

  ! Local arguments for FIT_SPLINE.
  REAL(KIND=R8), DIMENSION(1:NB_MAX)   :: XI
  REAL(KIND=R8), DIMENSION(1:NB_MAX,1:NCC) :: FX
  REAL(KIND=R8), DIMENSION(1:NB_MAX*NCC + 2*NCC) :: T
  REAL(KIND=R8), DIMENSION(1:NB_MAX*NCC) :: BCOEF
  INTEGER :: INFO

  ! Loop index.
  INTEGER :: I, J, NB

  ! Run all timing tests.
  DO NB = NB_MIN, NB_MAX, NB_STEP
     CALL TIME_RESET()
     ! Initialize arguments to FIT_SPLINE.
     DO I = 1, NB
        XI(I) = REAL(I, KIND=R8)
        DO J = 1, NCC
           FX(I,J) = REAL(I**2,R8)**(1.0_R8 / REAL(J,R8))
        END DO
     END DO
     ! Describe this test.
     WRITE (*,'(/,35A1,/,A,I4,/)') ('_',I=1,35), 'NB = ', NB
     ! Fit the spline (will record and print timings).
     CALL FIT_SPLINE(XI(1:NB), FX(1:NB,1:NCC), T, BCOEF, INFO)
     WRITE (*,'(A,I3)') ' INFO = ', INFO
  END DO

END PROGRAM TIME_SPLINE
