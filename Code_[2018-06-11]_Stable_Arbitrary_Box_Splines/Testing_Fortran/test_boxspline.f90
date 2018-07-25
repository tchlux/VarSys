
PROGRAM TEST_BOXSPLINE
  ! Define local variables for testing the box-spline code.
  INTEGER :: DIM, NUM_DVECS, NUM_PTS
  REAL :: START, FINISH
  ! Define local variables for reading direction vectors,
  ! multiplicities, and evaluation points.
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: UNIQUE_DVECS
  INTEGER,       DIMENSION(:),   ALLOCATABLE :: DVEC_MULTS
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: EVAL_PTS
  REAL(KIND=R8), DIMENSION(:),   ALLOCATABLE :: BOX_EVALS
  INTEGER :: ERROR
  ! Define an interface for the BOXSPLEV subroutine.
  INTERFACE 
    SUBROUTINE BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
      REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: UNIQUE_DVECS
      INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
      REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
      REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
      INTEGER,       INTENT(OUT)                 :: ERROR
    END SUBROUTINE BOXSPLEV
  END INTERFACE
  ! Body code for running all tests.

  
  ! __________________________________________________________________
  ! 
  ! TenP element, 1 points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_TenP_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_TenP_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_TenP_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,TenP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 2K points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_TenP_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_TenP_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_TenP_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,TenP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 4K points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_TenP_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_TenP_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_TenP_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,TenP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 1 points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_Diag_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_Diag_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_Diag_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,Diag,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 2K points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_Diag_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_Diag_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_Diag_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,Diag,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 4K points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_Diag_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_Diag_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_Diag_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,Diag,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! ZP element, 1 points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_ZP_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_ZP_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_ZP_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,ZP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! ZP element, 2K points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_ZP_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_ZP_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_ZP_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,ZP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! ZP element, 4K points, multiplicity 1
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_ZP_1.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_ZP_1.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_ZP_1.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '1,ZP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 1 points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_TenP_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_TenP_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_TenP_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,TenP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 2K points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_TenP_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_TenP_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_TenP_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,TenP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 4K points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_TenP_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_TenP_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_TenP_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,TenP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 1 points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_Diag_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_Diag_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_Diag_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,Diag,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 2K points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_Diag_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_Diag_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_Diag_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,Diag,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 4K points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_Diag_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_Diag_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_Diag_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,Diag,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! ZP element, 1 points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_ZP_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_ZP_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_ZP_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,ZP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! ZP element, 2K points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_ZP_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_ZP_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_ZP_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,ZP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! ZP element, 4K points, multiplicity 2
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_ZP_2.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_ZP_2.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_ZP_2.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '2,ZP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 1 points, multiplicity 3
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_TenP_3.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_TenP_3.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_TenP_3.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '3,TenP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 2K points, multiplicity 3
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_TenP_3.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_TenP_3.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_TenP_3.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '3,TenP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 4K points, multiplicity 3
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_TenP_3.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_TenP_3.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_TenP_3.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '3,TenP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 1 points, multiplicity 3
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_Diag_3.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_Diag_3.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_Diag_3.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '3,Diag,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 2K points, multiplicity 3
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_Diag_3.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_Diag_3.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_Diag_3.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '3,Diag,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! Diag element, 4K points, multiplicity 3
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_Diag_3.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_Diag_3.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_Diag_3.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '3,Diag,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 1 points, multiplicity 4
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_TenP_4.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_TenP_4.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_TenP_4.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '4,TenP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 2K points, multiplicity 4
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_TenP_4.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_TenP_4.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_TenP_4.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '4,TenP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 4K points, multiplicity 4
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_TenP_4.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_TenP_4.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_TenP_4.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '4,TenP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 1 points, multiplicity 5
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_1_TenP_5.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_1_TenP_5.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_1_TenP_5.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '5,TenP,1,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 2K points, multiplicity 5
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_2K_TenP_5.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_2K_TenP_5.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_2K_TenP_5.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '5,TenP,2K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------

  ! __________________________________________________________________
  ! 
  ! TenP element, 4K points, multiplicity 5
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs_4K_TenP_5.csv')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults_4K_TenP_5.csv')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts_4K_TenP_5.csv')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '5,TenP,4K,', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------


END PROGRAM TEST_BOXSPLINE
