This version implements almost all of Dr. Watson's initial suggestsion
from V1.1 (in the V1.1 folder). Particularly, the MATRIX_MINIMUM_NORM
function was rewritten correctly with LAPACK (not computing SQUARE at
all).

Speedup over Matlab ~ 15.5x