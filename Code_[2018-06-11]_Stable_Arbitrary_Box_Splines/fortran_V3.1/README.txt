This version updates 3.1 by removing all the "DEPTH" and "DEPTH+1" access operations from EVALUATE_BOX_SPLINE. Rather, the slices of memory at "DEPTH" and "DEPTH+1" are provided as input arguments to the subroutine.

Speedup over Matlab ~ 8.3x
