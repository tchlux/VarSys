This version responds to Watson's suggestsions from 2.1 and
reimplements column-vector storage (even though row vector storage
appears faster for the average use-case based on 2.2). This version
also makes the first attempt at performing all necessary array
allocation at the top-level subroutine and removing all allocations
from internal subroutine headers.