SUBROUTINE nearest(point, neighbor, data, sort_by_dim)
  ! Given "data" and "sort_by_dim" of equal size such that the values
  ! in "sorted_by_dim" are the indices of the sorted "data" along
  ! each dimension, lookup the nearest datum to "point"
  ! 
  USE ISO_FORTRAN_ENV
  REAL(KIND=REAL64), DIMENSION(:), INTENT(IN) :: point
  INTEGER, INTENT(OUT) :: neighbor
  REAL(KIND=REAL64), DIMENSION(:,SIZE(point)), INTENT(IN) :: data
  INTEGER, DIMENSION(SIZE(data)), INTENT(OUT) :: sorted_by_dim
  ! 
  

END SUBROUTINE nearest

SUBROUTINE lookup_index(val, array, index)
  USE ISO_FORTRAN_ENV
  REAL(KIND=REAL64), INTENT(IN) :: val
  REAL(KIND=REAL64), DIMENSION(:), INTENT(IN) :: array
  INTEGER, INTENT(OUT) :: index
  ! Local variables
  INTEGER :: below, above
  below = 1
  above = SIZE(array)
  find_index : DO WHILE ((above - below) .LE. 1)
     index = (below+above) / 2
     ! If the middle index was truncated to the lower value,
     IF (index .EQ. below) THEN
        below = above
     ELSE IF (array(index) .LT. val) THEN
        below = index
     ELSE IF (array(index) .GT. val) THEN
        above = index
     END IF
  END DO find_index
  index = below
END SUBROUTINE lookup_index

SUBROUTINE sort_by_dimension(data, sorted_by_dim)
  ! Given "data" (dimensions along columns), make the values of an
  ! equally sized integer array "sorted_by_dim" the indices in data of
  ! the points in sorted order. Then, data(sorted_by_dim) = sorted(data).
  ! 
  USE ISO_FORTRAN_ENV
  REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: data
  INTEGER, DIMENSION(SIZE(data)), INTENT(OUT) :: sorted_by_dim
  INTEGER :: dim
  ! Sort along each dimension
  sorting_dims : DO CONCURRENT (dim=1:SIZE(data,2))
     QSORTC( data(:,dim), sorted_by_dim(:,dim) )
  END DO sorting_dims
END SUBROUTINE sort_by_dimension


SUBROUTINE QSORTC(TO_SORT, SORTED_INDICES)
  ! This is a QuickSort routine adapted from Orderpack 2.0.
  !
  ! Also, this implementation incorporates ideas from "A Practical Introduction
  ! to Data Structures and Algorithm Analysis", by Clifford Shaffer.
  !
  ! It sorts real numbers into ascending numerical order
  ! and keeps an index of the value's original array position.
  !
  ! Author: Will Thacker, Winthrop University, July 2013.
  !
  ! QSORTC sorts the real array A and keeps the index of the value's
  ! original array position along with the value (integer array IDX).
  !
  ! On input:
  !
  ! TO_SORT(:) is the array to be sorted.
  !
  ! On output:
  !
  ! SORTED_INDICES(:) is a list of indices such that 
  !   TO_SORT(SORTED_INDICES(i)) = the i-th point in ascending sorted order.
  !
  ! This has been modified in 2016 by Thomas Lux (tchlux@vt.edu) to
  ! not have an effect on the input array and rather to only return
  ! the indices of input elements sorted.
  ! 
  REAL(KIND=REAL64), DIMENSION(:), INTENT(IN):: TO_SORT
  INTEGER, DIMENSION(SIZE(A)), INTENT(OUT):: SORTED_INDICES
  ! Local variables
  REAL(KIND=REAL64), DIMENSION(:) :: LOCAL_COPY
  INTEGER:: I   ! Loop iteration variable.
  ! Initialize the copy array of values
  LOCAL_COPY = TO_SORT
  ! Initialize the array of original positions.
  DO CONCURRENT (I=1:SIZE(A))
     SORTED_INDICES(I)=I
  END DO
  CALL QSORTC_HELPER(LOCAL_COPY, SORTED_INDICES, 1, SIZE(TO_SORT))
  RETURN
CONTAINS
  RECURSIVE SUBROUTINE QSORTC_HELPER(A, IDX, ISTART, ISTOP)
    ! This internal recursive subroutine performs the recursive quicksort
    ! algorithm.  It is needed because the boundaries of the part of the
    ! array being sorted change and the initial call to the sort routine
    ! does not need to specify boundaries since, generally, the user will
    ! want to sort the entire array passed.
    !
    ! On input:
    !
    ! A(:) contains a subpart to be sorted.
    !
    ! IDX(i) contains the initial position of the value A(i) before sorting.
    !
    ! ISTART is the starting position of the subarray to be sorted.
    !
    ! ISTOP is the ending position of the subarray to be sorted.
    !
    ! On output:
    !
    ! A(ISTART:ISTOP) will be sorted.
    !
    ! IDX(i) contains the original position for the value at A(i).
    !
    !
    REAL(KIND=REAL64), DIMENSION(:), INTENT (IN OUT):: A
    INTEGER, DIMENSION(SIZE(A)), INTENT(IN OUT):: IDX
    INTEGER, INTENT(IN):: ISTART, ISTOP
    !  Local variables
    INTEGER:: ILEFT ! A position on the left to be swapped with value at IRIGHT.
    INTEGER:: IMID ! The middle position used to select the pivot.
    INTEGER:: IRIGHT ! A position on the right to be swapped with value at ILEFT.
    INTEGER:: ITEMP  ! Used for swapping within IDX.
    REAL(KIND=REAL64):: ATEMP ! Used for swapping.
    REAL(KIND=REAL64):: PIVOT ! Holds the temporary pivot.
    ! INSMAX is used to stop recursively dividing the array and to instead
    ! use a sort that is more efficient for small arrays than quicksort.
    !
    ! The best cutoff point is system dependent.
    INTEGER, PARAMETER:: INSMAX=24
    ! Check to see if we have enough values to make quicksort useful.
    ! Otherwise let the insertion sort handle it.
    IF ((ISTOP - ISTART) .LT. INSMAX) THEN
       CALL INSERTION(A, IDX, ISTART, ISTOP)
    ELSE
       ! Use the median of the first, middle and last items for the pivot
       ! and place the median (pivot) at the end of the list.
       ! Putting it at the end of the list allows for a guard value to keep
       ! the loop from falling off the right end of the array (no need to
       ! check for at the end of the subarray EACH time through the loop).
       IMID = (ISTART + ISTOP)/2
       IF (A(ISTOP) .LT. A(ISTART)) THEN
          ATEMP = A(ISTART)
          A(ISTART) = A(ISTOP)
          A(ISTOP) = ATEMP
          ITEMP = IDX(ISTART)
          IDX(ISTART) = IDX(ISTOP)
          IDX(ISTOP) = ITEMP
       END IF
       IF (A(IMID) .LT. A(ISTOP)) THEN
          ATEMP = A(ISTOP)
          A(ISTOP) = A(IMID)
          A(IMID) = ATEMP
          ITEMP = IDX(ISTOP)
          IDX(ISTOP) = IDX(IMID)
          IDX(IMID) = ITEMP
          IF (A(ISTOP) .LT. A(ISTART)) THEN
             ATEMP = A(ISTOP)
             A(ISTOP) = A(ISTART)
             A(ISTART) = ATEMP
             ITEMP = IDX(ISTOP)
             IDX(ISTOP) = IDX(ISTART)
             IDX(ISTART) = ITEMP
          END IF
       END IF
       ! Now, the first position has a value that is less or equal to the
       ! partition. So, we know it belongs in the left side of the partition
       ! and we can skip it. Also, the pivot is at the end.  So, there is
       ! no need to compare the pivot with itself.
       PIVOT = A(ISTOP)
       ILEFT = ISTART + 1
       IRIGHT = ISTOP - 1
       DO WHILE (ILEFT < IRIGHT)
          ! Find a value in the left side that is bigger than the pivot value.
          ! Pivot is at the right end so ILEFT will not fall off the end
          ! of the subarray.
          DO WHILE (A(ILEFT) .LT. PIVOT)
             ILEFT = ILEFT + 1
          END DO
          DO WHILE (IRIGHT .NE. ILEFT)
             IF (A(IRIGHT) .LT. PIVOT) EXIT
             IRIGHT = IRIGHT - 1
          END DO
          ! Now we have a value bigger than pivot value on the left side that can be
          ! swapped with a value smaller than the pivot on the right side.
          !
          ! This gives us all values less than pivot on the left side of the
          ! array and all values greater than the pivot on the right side.
          ATEMP = A(IRIGHT)
          A(IRIGHT) = A(ILEFT)
          A(ILEFT) = ATEMP
          ITEMP = IDX(IRIGHT)
          IDX(IRIGHT) = IDX(ILEFT)
          IDX(ILEFT) = ITEMP
       END DO
       !
       ! The last swap was in error (since the while condition is not checked
       ! until after the swap is done) so we swap again to fix it.
       ! This is done (once) rather than having an if (done many times) in the
       ! loop to prevent the swapping.
       ATEMP = A(IRIGHT)
       A(IRIGHT) = A(ILEFT)
       A(ILEFT) = ATEMP
       ITEMP = IDX(IRIGHT)
       IDX(IRIGHT) = IDX(ILEFT)
       IDX(ILEFT) = ITEMP
       ! Put the pivot value in its correct spot (between the 2 partitions)
       ! When the WHILE condition finishes, ILEFT is greater than IRIGHT.
       ! So, ILEFT has the position of the first value in the right side.
       ! This is where we can put the pivot (and where it will finally rest,
       ! so no need to look at it again).  Also, place the first value of
       ! the right side (being displaced by the pivot) at the end of the
       ! subarray (since it is bigger than the pivot).
       ATEMP = A(ISTOP)
       A(ISTOP) = A(ILEFT)
       A(ILEFT) = ATEMP
       ITEMP = IDX(ISTOP)
       IDX(ISTOP) = IDX(ILEFT)
       IDX(ILEFT) = ITEMP
       CALL QSORTC_HELPER(A, IDX, ISTART, ILEFT-1)
       CALL QSORTC_HELPER(A, IDX, ILEFT+1, ISTOP)
    END IF
  END SUBROUTINE QSORTC_HELPER

  SUBROUTINE INSERTION(A, IDX, ISTART, ISTOP)
    ! This subroutine performs an insertion sort used for sorting
    ! small subarrays efficiently.
    !
    ! This subroutine sorts a subarray of A (between positions ISTART
    ! and ISTOP) keeping a record of the original position (array IDX).

    ! On input:
    !
    ! A(:) contains a subpart to be sorted.
    !
    ! IDX(i) contains the initial position of the value A(i) before sorting.
    !
    ! ISTART is the starting position of the subarray to be sorted.
    !
    ! ISTOP is the ending position of the subarray to be sorted.
    !
    ! On output:
    !
    ! A(ISTART:ISTOP) will be sorted.
    !
    ! IDX(i) contains the original position for the value at A(i).
    !
    REAL(KIND=REAL64), DIMENSION(:), INTENT (IN OUT):: A
    INTEGER, DIMENSION(SIZE(A)), INTENT(IN OUT)::IDX
    INTEGER, INTENT(IN):: ISTART, ISTOP
    ! Local variables.
    REAL(KIND=REAL64):: AMIN  ! Temporary minimum.
    REAL(KIND=REAL64):: ATEMP  ! The value to be inserted.
    INTEGER:: I    ! Index variable.
    INTEGER:: IABOVE ! Index to find insertion point.
    INTEGER:: IMIN ! Temporary minimum position.
    INTEGER:: ITEMP ! Temporary for swapping.
    IF (ISTOP .EQ. ISTART) THEN
       RETURN
    END IF
    ! Find the smallest and put it at the top as a "guard" so there is
    ! no need for the DO WHILE to check if it is going past the top.
    AMIN = A(ISTART)
    IMIN = ISTART
    DO I=ISTOP,ISTART+1,-1
       IF (A(I) .LT. AMIN) THEN
          AMIN = A(I)
          IMIN = I
       END IF
    END DO
    A(IMIN) = A(ISTART)
    A(ISTART) = AMIN
    ITEMP = IDX(ISTART)
    IDX(ISTART) = IDX(IMIN)
    IDX(IMIN) = ITEMP
    ! Insertion sort the rest of the array.
    DO I=ISTART+2,ISTOP
       ATEMP = A(I)
       ITEMP = IDX(I)
       IABOVE = I - 1
       IF (ATEMP .LT. A(IABOVE)) THEN
          A(I) = A(IABOVE)
          IDX(I) = IDX(IABOVE)
          IABOVE = IABOVE - 1
          ! Stop moving items down when the position for "insertion" is found.
          !
          ! Do not have to check for "falling off" the beginning of the
          ! array since the smallest value is a guard value in the first position.
          DO WHILE (ATEMP .LT. A(IABOVE))
             A(IABOVE+1) = A(IABOVE)
             IDX(IABOVE+1) = IDX(IABOVE)
             IABOVE = IABOVE - 1
          END DO
       END IF
       A(IABOVE+1) = ATEMP
       IDX(IABOVE+1) = ITEMP
    END DO
  END SUBROUTINE INSERTION
END SUBROUTINE QSORTC
