''''''


import cython
import numpy

class NotFortranCompatible(Exception): pass

#      Wrapper for fortran function allocate_max_lapack_work     
# =================================================

cdef extern:
    void c_allocate_max_lapack_work( int* dim, int* num_dvecs, int* size_output )

@cython.boundscheck(False)
@cython.wraparound(False)
def allocate_max_lapack_work( int dim, int num_dvecs, size_output=None ):
    '''! 8) ALLOCATE_MAX_LAPACK_WORK
    !
    ! Return an allocated real array that has a size equal to the
    ! maximum requested LAPACK work array size across all routines
    ! that will be executed in the evaluation of a Box Spline.
    !
    ! Output:
    !   WORK -- Real array with size large enough to accomadate all
    !           LAPACK subroutines that will be used to evaluated a
    !           Box Spline given the dimension and number of vectors.
    !'''
    # Prepare for fortran function call (initialize optionals)
    if (type(size_output) == type(None)):
        size_output = 1
    cdef int local_size_output = size_output
    
    
    # Make fortran function call
    c_allocate_max_lapack_work(&dim, &num_dvecs, &local_size_output)
    # Return appropriate values based on "INTENT" from fortran code
    
    return local_size_output



