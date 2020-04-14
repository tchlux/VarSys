      REAL FUNCTION FMIN(AX,BX,F,TOL)
C***BEGIN PROLOGUE  FMIN
C***DATE WRITTEN   730101  (YYMMDD)
C***REVISION DATE  730101  (YYMMDD) 
C***CATEGORY NO.  G1A2
C***KEYWORDS  ONE-DIMENSIONAL MINIMIZATION, UNIMODAL FUNCTION
C***AUTHOR  BRENT, R.
C***PURPOSE  An approximation to the point where F attains a minimum on
C            the interval (AX,BX) is determined as the value of the 
C            function FMIN.
C***DESCRIPTION
C
C     From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C     The method used is a combination of golden section search and
C     successive parabolic interpolation.  Convergence is never much 
C     slower than that for a Fibonacci search.  If F has a continuous 
C     second derivative which is positive at the minimum (which is not
C     at AX or BX), then convergence is superlinear, and usually of the 
C     order of about 1.324....
C
C     The function F is never evaluated at two points closer together
C     than EPS*ABS(FMIN) + (TOL/3), where EPS is approximately the 
C     square root of the relative machine precision.  If F is a unimodal
C     function and the computed values of F are always unimodal when
C     separated by at least EPS*ABS(XSTAR) + (TOL/3), then FMIN 
C     approximates the abcissa of the global minimum of F on the 
C     interval AX,BX with an error less than 3*EPS*ABS(FMIN) + TOL.  
C     If F is not unimodal, then FMIN may approximate a local, but 
C     perhaps non-global, minimum to the same accuracy.
C
C     This function subprogram is a slightly modified version of the
C     ALGOL 60 procedure LOCALMIN given in Richard Brent, Algorithms for
C     Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
C
C INPUT PARAMETERS
C
C  AX    (real)  left endpoint of initial interval
C  BX    (real) right endpoint of initial interval
C  F     Real function of the form REAL FUNCTION F(X) which evaluates 
C          F(X)  for any  X in the interval  (AX,BX)
C        Must be declared EXTERNAL in calling routine.
C  TOL   (real) desired length of the interval of uncertainty of the 
C        final result ( .ge. 0.0)
C
C
C OUTPUT PARAMETERS
C
C FMIN   abcissa approximating the minimizer of F
C AX     lower bound for minimizer
C BX     upper bound for minimizer
C
C***REFERENCES  RICHARD BRENT, ALGORITHMS FOR MINIMIZATION WITHOUT 
C                 DERIVATIVES, PRENTICE-HALL, INC. (1973).
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FMIN
      REAL  AX,BX,F,TOL
      REAL  A,B,C,D,E,EPS,XM,P,Q,R,TOL1,TOL2,U,V,W
      REAL  FU,FV,FW,FX,X
      REAL  ABS,SQRT,SIGN
C***FIRST EXECUTABLE STATEMENT  FMIN
      C = 0.5*(3. - SQRT(5.0))
C
C  C is the squared inverse of the golden ratio
C
C  EPS is approximately the square root of the relative machine
C  precision.
C
      EPS = 1.0
   10 EPS = EPS/2.0
      TOL1 = 1.0 + EPS
      IF (TOL1 .GT. 1.0) GO TO 10
      EPS = SQRT(EPS)
C
C  initialization
C
      A = AX
      B = BX
      V = A + C*(B - A)
      W = V
      X = V
      E = 0.0
      FX = F(X)
      FV = FX
      FW = FX
C
C  main loop starts here
C
   20 XM = 0.5*(A + B)
      TOL1 = EPS*ABS(X) + TOL/3.0
      TOL2 = 2.0*TOL1
C
C  check stopping criterion
C
      IF (ABS(X - XM) .LE. (TOL2 - 0.5*(B - A))) GO TO 90
C
C is golden-section necessary
C
      IF (ABS(E) .LE. TOL1) GO TO 40
C
C  fit parabola
C
      R = (X - W)*(FX - FV)
      Q = (X - V)*(FX - FW)
      P = (X - V)*Q - (X - W)*R
      Q = 2.0*(Q - R)
      IF (Q .GT. 0.0) P = -P
      Q = ABS(Q)
      R = E
      E = D
C
C  is parabola acceptable
C
   30 IF (ABS(P) .GE. ABS(0.5*Q*R)) GO TO 40
      IF (P .LE. Q*(A - X)) GO TO 40
      IF (P .GE. Q*(B - X)) GO TO 40
C
C  a parabolic interpolation step
C
      D = P/Q
      U = X + D
C
C  F must not be evaluated too close to AX or BX
C
      IF ((U - A) .LT. TOL2) D = SIGN(TOL1, XM - X)
      IF ((B - U) .LT. TOL2) D = SIGN(TOL1, XM - X)
      GO TO 50
C
C  a golden-section step
C
   40 IF (X .GE. XM) E = A - X
      IF (X .LT. XM) E = B - X
      D = C*E
C
C  F must not be evaluated too close to X
C
   50 IF (ABS(D) .GE. TOL1) U = X + D
      IF (ABS(D) .LT. TOL1) U = X + SIGN(TOL1, D)
      FU = F(U)
C
C  update  A, B, V, W, and X
C
      IF (FU .GT. FX) GO TO 60
      IF (U .GE. X) A = X
      IF (U .LT. X) B = X
      V = W
      FV = FW
      W = X
      FW = FX
      X = U
      FX = FU
      GO TO 20
   60 IF (U .LT. X) A = U
      IF (U .GE. X) B = U
      IF (FU .LE. FW) GO TO 70
      IF (W .EQ. X) GO TO 70
      IF (FU .LE. FV) GO TO 80
      IF (V .EQ. X) GO TO 80
      IF (V .EQ. W) GO TO 80
      GO TO 20
   70 V = W
      FV = FW
      W = U
      FW = FU
      GO TO 20
   80 V = U
      FV = FU
      GO TO 20
C
C  end of main loop
C
   90 FMIN = X
      RETURN
      END
