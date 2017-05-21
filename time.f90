!
! Fortran 95 version of the time routines.
! This version uses only calls to the Fortran 95
! intrinsic functions:
!
! cpu_time
! system_clock
! date_and_time
!
! thus it should be quite portable
!
! The only nonstandard call is to the subroutine
! system, for executing the command 'hostname'
!
! the function perpetual_calendar returns the day of the week.
!
subroutine secund( t )

           ! returns the cpu time in seconds

  implicit none
  real( kind = 8 ), intent( out ) :: t
  real :: cputime

  call cpu_time( cputime )
  t = dble( cputime )
  
end subroutine secund

subroutine elapsec( t )

 ! returns the elapsed time in seconds

 ! this routines calls the standard Fortran 95 intrinsic procedure
 ! system_clock, which return a clock counter, a rate (number of counts per
 ! second), and a maxrate (maximum value that the counter can assume)
 !
 ! these values depend on the default integer size:
 !
 !                  rate          maxrate          counter restart time
 !
 !    integer*4     10000       2147483647              59.6 hours
 !    integer*8    1000000  9223372036854775807     292471.2 years
 !
 ! as one can see, for large 32-bit jobs (integer*4) there is the chance
 ! that the clock counter can reach the maximum and be restarted. The routine
 ! tries to detect that and take care of it, but it can still fail if two 
 ! consecutive call to system clock are done more than 59.6 hours apart, 
 ! which hopefully is unlikely. For 64-bit this will be a problem only
 ! for very large (geologic time) jobs.
 !
 ! There is also a different behavior associated with different compilers:
 ! with g95 the first call to system_clock always returns a count of 0,
 ! while for the intel compiler the clock count is relative to 00.00 CUT
 ! of January 1 1970. This means that the counter in Intel-compiled 32-bit
 ! jobs can be restarted in the middle of a job also for small jobs. 
 ! The routine should be able to handle that if it happens.
 
  implicit none
  real( kind = 8 ), intent( out ) :: t
  integer :: count, rate, maxr
  integer, save :: ifirst = 0, lastcount, rate100
  real( kind = 8 ), save :: t0, tprec, tmax

  if( ifirst .eq. 0 ) then
    call system_clock( count, rate, maxr )
    t0 = dble( count ) / dble( rate )
    lastcount = count
    rate100 = rate / 100
    tprec = 0.0d0
    tmax = dble( maxr ) / dble( rate )
    ifirst = 1
  endif

  call system_clock( count, rate, maxr )

  if( count .lt. ( lastcount  - rate100) ) tprec = tprec + tmax  ! the counter has been restarted

  lastcount = count
  t = dble( count ) / dble( rate ) + tprec - t0

end subroutine elapsec

subroutine date1( x )

          ! returns the date as in 'Wed Oct  5 18:36:02 2005'

  implicit none
  character( 24 ), intent( out ) :: x
  integer :: dt(8)

  character(3) :: weekday, month
  character(3), parameter :: months(12) = (/ 'Jan', 'Feb', 'Mar', &
  'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)
  character(3), external :: perpetual_calendar

  call date_and_time( values = dt )

  weekday = perpetual_calendar( dt(2), dt(3), dt(1) )

  month = '???'
  if( dt(2) .gt. 0 .and. dt(2) .le. 12 ) month = months(dt(2))

  write( x, &
        '(a3,1x,a3,1x,i2,1x,i2,'':'',i2,'':'',i2,1x,i4)' ) &
        weekday, month, dt(3), dt(5), dt(6), dt(7), dt(1)

end subroutine date1

subroutine chartime(x)

  character(24), intent( out ) :: x

  call date1( x )

end subroutine chartime

subroutine get_host(hostname)
   
  implicit none
  character(80), intent( out ) :: hostname

  call system('hostname>hname.tmp')
  open(1,file='hname.tmp')
  read(1,'(a80)',end=200) hostname

200 close(1,status='delete')

end subroutine get_host

function perpetual_calendar( month, day, year ) result ( weekday )

  implicit none
  character( 3 ) :: weekday
  integer, intent( in ) :: month, day, year

           ! the reference day is Saturday Jan 1 2000.

  character( 3 ), dimension( 7), parameter :: days = &
               (/'Sat', 'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri' /)

  integer, parameter :: mday( 12 ) = (/ 31, 28, 31, 30, 31, 30, 31, &
                        31, 30, 31, 30, 31 /)
  integer :: ndays, i

  weekday = '???'
  if( month .le. 0 .or. month .gt. 12 ) return
  if( day .le. 0 ) return
  if( month .ne. 2 ) then
    if( day .gt. mday( month ) ) return
  else
    if( leap( year ) )then
      if( day .gt. 28 ) return
    else
      if( day .gt. 29 ) return
    endif
  endif

  if ( year .ge. 2000 ) then

           ! counting forward
    ndays = -1

    do i = 2000, year - 1, 1
      ndays = ndays + 365
      if( leap( i ) ) ndays = ndays + 1
    enddo

    do i = 1, month - 1, 1
      ndays = ndays + mday( i )
      if( i .eq. 2 .and. leap( year ) ) ndays = ndays + 1
    enddo

    ndays = ndays + day
    weekday = days( mod( ndays, 7 )  + 1 )
  else
          ! counting backwards (to cover the case of time travel)
    ndays = 0

    do i = 1999, year + 1, -1
      ndays = ndays + 365
      if( leap( i ) ) ndays = ndays + 1
    enddo

    do i = 12, month + 1, -1
      ndays = ndays + mday(i)
      if( i .eq.2 .and. leap( year ) ) ndays = ndays + 1
    enddo

    ndays = ndays + mday( month ) - day
    if( month .eq.2 .and. leap( year ) ) ndays = ndays + 1

    weekday = days( 7 - mod( ndays, 7 ) )

  endif

  contains

  function leap( year ) result ( isleap )

    implicit none
    integer :: year
    logical :: isleap

    if ( mod( year, 400 ) .eq. 0 ) then
      isleap = .true.
    else if( mod( year, 100 ) .eq. 0 ) then
      isleap = .false.
    else if( mod( year, 4 ) .eq. 0 ) then
      isleap = .true.
    else
      isleap = .false.
    endif

  end function leap

end function perpetual_calendar
