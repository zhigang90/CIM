module kinds

  implicit none

          ! named constants for integer kinds

  integer, parameter ::         &
    i_1 = selected_int_kind(2), & ! integer*1
    i_2 = selected_int_kind(4), & ! integer*2
    i_4 = selected_int_kind(9), & ! integer*4
    i_8 = selected_int_kind(18)   ! integer*8

          ! sizes of the various integer kinds, that is, 
          ! how many integers fit in a double precision word

  integer, parameter ::             &
    intsize = 64 / bit_size(0),     & ! size of default integer
    i1size  = 64 / bit_size(0_i_1), & ! size of integer*1
    i2size  = 64 / bit_size(0_i_2), & ! size of integer*2
    i4size  = 64 / bit_size(0_i_4), & ! size of integer*4
    i8size  = 64 / bit_size(0_i_8)    ! size of integer*8

  integer, parameter :: bytesize = 8  ! size of characters

         ! named constants for real kinds

  integer, parameter ::                            &
    r_4 = kind(1.0),                               & ! real*4 
    r_8 = selected_real_kind(2*precision(1.0_r_4))   ! real*8

end module kinds
