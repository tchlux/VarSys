! Given minimum values for the separation of X.
! 
! D = ACCURACY^3                          
! C. = ACCURACY * SEP_MAX                 
! A = (3 * ACCURACY * SEP_MAX) / ACCURACY^3
! B = (3 * ACCURACY^2 * SEP_MAX) / ACCURACY^3
! 
! DX = 2 * ((3 * ACCURACY * SEP_MAX) / ACCURACY^3) * ACCURACY + (3 * ACCURACY^2 * SEP_MAX / ACCURACY^3)
!    = 6 * SEP_MAX / ACCURACY + 3 * SEP_MAX / ACCURACY
!    = 9 * SEP_MAX / ACCURACY
! 
! W = ACCURACY
! INV_SLOPE = ACCURACY / ACCURACY
! A = (ACCURACY / ACCURACY) * (9 * SEP_MAX / ACCURACY)
! B = (ACCURACY / ACCURACY) * (9 * SEP_MAX / ACCURACY)
! 
! 
! 
! Given maximum values for the separation of X.
! 
! D = MAX_VAL^3
! C. = MAX_VAL * SEP_MAX
! A = (3 * MAX_VAL * SEP_MAX) / MAX_VAL^3
! B = 3 * MAX_VAL^2 * SEP_MAX / MAX_VAL^3
! 
! DX = 2 * ((3 * MAX_VAL * SEP_MAX) / MAX_VAL^3) * MAX_VAL + (3 * MAX_VAL^2 * SEP_MAX / MAX_VAL^3)
!    = 6 * SEP_MAX / MAX_VAL + 3 * SEP_MAX / MAX_VAL
!    = 9 * SEP_MAX / MAX_VAL
! 
! W = SEP_MAX
! INV_SLOPE = SEP_MAX / ACCURACY
! A = (SEP_MAX / ACCURACY) * (9 * SEP_MAX / MAX_VAL)
!   = 9 * SEP_MAX^2 / (ACCURACY * MAX_VAL)
! B = (SEP_MAX / ACCURACY) * (9 * SEP_MAX / MAX_VAL)
!   = 9 * SEP_MAX^2 / (ACCURACY * MAX_VAL)
! 
! With "MAX_VAL" = "SEP_MAX" we have:
! 
!  A*B = 81 * (SEP_MAX / ACCURACY) * (SEP_MAX / ACCURACY)
!      ~ 81 * SEP_MAX^2 / ACCURACY^2
! TEMP = 2 * SEP_MAX / ACCURACY * SQRT(81*(SEP_MAX / SEP_MAX)^2 - <same>)
!      = 2 * SEP_MAX / ACCURACY
! (2 * SEP_MAX / ACCURACY) + (SEP_MAX / ACCURACY) * 162
! 60 + 2*(2 * SEP_MAX / ACCURACY) + (SEP_MAX / ACCURACY) * SEP_MAX
! 
! 
! Given a minimum value for one derivative and maximum for the other..
! 
! DX = SEP_MAX / ACCURACY
! DX = ACCURACY
! A = SEP_MAX / ACCURACY * ACCURACY
!   = SEP_MAX
! B = SEP_MAX / ACCURACY * ACCURACY
!   = SEP_MAX
! TAU = SQRT( SEP_MAX^2 ) - SEP_MAX
!     = SEP_MAX
! TEMP = ((SEP_MAX / ACCURACY) / ACCURACY) ^ (3/4)
!      = (SEP_MAX / ACCURACY^2) ^ (3/4)
! ALPHA = (SEP_MAX / ACCURACY^2) * (SEP_MAX / ACCURACY ...) / ACCURACY
!       = SEP_MAX^2 / ACCURACY^4
! GAMMA = (SEP_MAX / ACCURACY) / ((SEP_MAX / ACCURACY^2) * ACCURACY)
!       = ACCURACY / SEP_MAX
! BETA = (SEP_MAX / ACCURACY) * (SEP_MAX^2 / ACCURACY) / (ACCURACY)
!      = (SEP_MAX^3 / ACCURACY^3)
! 
