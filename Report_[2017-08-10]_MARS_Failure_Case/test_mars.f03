PROGRAM TEST_MARS
  IMPLICIT NONE
  INTEGER :: N, P, NK, MI
  REAL,             DIMENSION(10,3) :: X
  REAL,             DIMENSION(10)    :: Y
  REAL,             DIMENSION(10)    :: W
  INTEGER,          DIMENSION(3)    :: LX
  REAL,             DIMENSION(219)    :: FM
  INTEGER,          DIMENSION(191)    :: IM
  REAL,             DIMENSION(269)    :: SP
  DOUBLE PRECISION, DIMENSION(161)    :: DP
  INTEGER,          DIMENSION(36)    :: MM
  
  ! See descriptions of the above variables in 'marspack.f'
  ! documentation for subroutine 'mars' starting at line 126

  N = 10    ! Number of points
  P = 3    ! Number of dimensions
  NK = 10   ! Maximum number of basis functions
  MI = 3   ! Maximum interaction between dimensions
  W = 1.0   ! Weights for each point
  LX = 1    ! Flags for each dimension (1 = ordinal)
  FM = 1.0  ! Holder for MARS model
  IM = 1    ! Holder for MARS model
  SP = 1.0  ! Real workspace array for MARS
  DP = 1.0  ! Double precision workspace array for MARS
  MM = 1    ! Integer workspace array for MARS
  
  ! Size of fm = 3+nk*(5*mi+nmcv+6)+2*p+ntcv
  !            = 3+10*(5*3+0+6)+2*3+0)
  !            = 219
  ! Size of im = 21+nk*(3*mi+8)
  !            = 21+10*(3*3+8)
  !            = 191
  ! Size of sp = n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk
  !            = n*(max(10+1,2)+3)+max(3*10+5*10+3,2*3,4*10)+2*3+4*10
  !            = 269
  ! Size of dp = max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)
  !            = max(10*10,(10+1)*(10+1))+max((10+2)*(0+3),4*10)
  !            = 161
  ! Size of mm = n*p+2*max(mi,nmcv)
  !            = 10*3+2*max(3,0)
  !            = 36
  
  ! Initialization of X that causes segmentation fault
  X(1,:) = (/ 0.0, 0.0, 0.0 /)
  X(2,:) = (/ 0.0, 0.0, 1.0 /)
  X(3,:) = (/ 1.0, 0.0, 0.0 /)
  X(4,:) = (/ 1.0, 0.0, 1.0 /)
  X(5,:) = (/ 1.0, 1.0, 0.0 /)
  X(6,:) = (/ 1.0, 1.0, 1.0 /)
  ! X(7,:) = (/ 0.0, 0.0, 0.1110494169905608 /)
  X(8,:) = (/ 0.0, 0.0, 0.8883953359244864 /)
  X(9,:) = (/ 0.2, 0.2, 0.4441976679622432 /)
  X(10,:) = (/ 1.0, 0.2, 0.9994447529150472 /)
  ! Y values do not matter for this bug
  Y = 1.0

  PRINT *, "SHAPE(X): ", SHAPE(X)
  PRINT *,
  PRINT *, 'Calling MARS...'

  ! This should build a model out of the data and store the output in
  ! "FM" and "IM" to be used later by FMOD.
  CALL MARS(N,P,X,Y,W,NK,MI,LX,FM,IM,SP,DP,MM)

  PRINT *, 'Done calling MARS.'

END PROGRAM TEST_MARS
