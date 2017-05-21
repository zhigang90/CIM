module memory

  use kinds
  implicit none
  private

  real (kind=r_8), pointer :: bl(:) => null()  ! pointer for the memory system
! real (kind=8) :: bl
! common /big/bl(1)

  integer :: lcore = 0! memory that has been allocated


  integer, parameter :: maxmark = 100, maxall = 20000
  integer :: nreq, mark, maxmem, marks(maxmark), nadr(maxall)
! common /mmanag/lcore,nreq,maxmem,mark,marks(maxmark),nadr(maxall)

  public :: retall, getmem, getint, getint_1, getint_2, getint_4,&
            getint_8, getbyte, retmem, mmark, retmark, memstat, &
            mem_alloc, mem_dealloc

  public :: bl

  contains

!================================================================

   subroutine mem_alloc( nwords, ierr )
 
     implicit none
     integer, intent( inout ) :: nwords
     integer, intent( out ) :: ierr
     integer :: istat = 0
 
 
     allocate ( bl( nwords ), stat = istat )
 
     if ( istat .eq. 0 ) then
       lcore = nwords
       ierr = 0
     else
       ierr = -1
       lcore = 0
     endif
 
     return
 
   end subroutine mem_alloc

!================================================================

   subroutine mem_dealloc( ierr )
 
     implicit none
     integer, intent( out ) :: ierr
     integer :: istat
 
     ierr = 1
 
     deallocate( bl, stat = istat )
 
     if ( istat .eq. 0 ) ierr = 0
 
     return
 
   end subroutine mem_dealloc
  
!================================================================

  subroutine retall(ioffset)
          
  implicit none
  integer :: ioffset
           ! (re)initializes the memory allocation system

    nreq = 0
    nadr(1) = ioffset
    maxmem = 0
    mark = 0

  end subroutine retall

!================================================================

  subroutine getmem( amount, addr )

           ! reserves amount words of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    implicit none
    integer, intent( in ) :: amount
    integer, intent( out ) :: addr
    integer :: ixx

    nreq = nreq + 1
    if ( nreq > maxall - 1 ) call nerror( 1, 'getmem', & 
                'number of requests exceeds maximum', nreq, maxall )
!   if ( nreq > maxall - 1 ) call abort

    addr = nadr( nreq ) + 1
    ixx  = addr + amount - 1
    if ( ixx > lcore ) call nerror( 1, 'getmem', &
      ' memory overflow request no., amount needed', nreq, ixx-nadr(1) )
!   if ( ixx > lcore ) call abort

    nadr( nreq + 1 ) = ixx
    if ( ixx > maxmem ) maxmem = ixx
!     write(*,*) 'memory allocation request no=',nreq,'addr=',addr,&
!                 'ceiling=',ixx

  end subroutine getmem

!================================================================

  subroutine getint( amint, addr )

           ! reserves amint default integer words of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    integer, intent( in ) :: amint
    integer, intent( out ) :: addr
    integer :: amount

    amount = amint / intsize 
    if( mod( amint, intsize ) .ne. 0 ) amount = amount + 1

    call getmem( amount, addr )
!    write(*,*)'memory  from getint amount= ', amount, ' addr= ', addr

  end subroutine getint

!================================================================

  subroutine getint_1( amint, addr )

           ! reserves amint integer*1 words of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    integer, intent( in ) :: amint
    integer, intent( out ) :: addr
    integer :: amount

    amount = amint / i1size 
    if( mod( amint, i1size ) .ne. 0 ) amount = amount + 1

    call getmem( amount, addr )
!    write(*,*)'memory  from getint amount= ', amount, ' addr= ', addr

  end subroutine getint_1

!================================================================

  subroutine getint_2( amint, addr )

           ! reserves amint integer*2 words of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    integer, intent( in ) :: amint
    integer, intent( out ) :: addr
    integer :: amount

    amount = amint / i2size 
    if( mod( amint, i2size ) .ne. 0 ) amount = amount + 1

    call getmem( amount, addr )
!    write(*,*)'memory  from getint amount= ', amount, ' addr= ', addr

  end subroutine getint_2

!================================================================

  subroutine getint_4( amint, addr )

           ! reserves amint integer*4 words of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    integer, intent( in ) :: amint
    integer, intent( out ) :: addr
    integer :: amount

    amount = amint / i4size 
    if( mod( amint, i4size ) .ne. 0 ) amount = amount + 1

    call getmem( amount, addr )
!    write(*,*)'memory  from getint amount= ', amount, ' addr= ', addr

  end subroutine getint_4

!================================================================

  subroutine getint_8( amint, addr )

           ! reserves amint integer*8 words of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    integer, intent( in ) :: amint
    integer, intent( out ) :: addr
    integer :: amount

    amount = amint / i8size 
    if( mod( amint, i8size ) .ne. 0 ) amount = amount + 1

    call getmem( amount, addr )
!    write(*,*)'memory  from getint amount= ', amount, ' addr= ', addr

  end subroutine getint_8

!================================================================

  subroutine getbyte( ambyte, addr )

           ! reserves ambyte bytes of space in the bl array
           ! and returns the starting address of the reserved
           ! block

    integer, intent( in ) :: ambyte
    integer, intent( out ) :: addr
    integer :: amount

    amount = ambyte / bytesize
    if( mod( ambyte, bytesize ) .ne. 0 ) amount = amount + 1

    call getmem( amount, addr )

  end subroutine getbyte

!================================================================

  subroutine retmem ( n )

           ! removes the reservation for the last n blocks

    integer, intent( in ) :: n
    integer i

    nreq = max( 1, nreq - n ) ! the first assignement must stay

!    write(*,*) 'memory memory deallocation, last valid request and ceiling=',&
!               nreq,nadr(nreq+1)

  end subroutine retmem

!================================================================

  subroutine mmark

!  this subroutine puts down a mark in the memory management if called.
!  its purpose is to facilitate the return of unneded memory
!  retmem serves well if a few block have to be returned but sometimes
!  one does not know the exact number of blocks called. retall is
!  less practical in such  cases because it return everything.
!  mmark is used in conjunction with retmark. If retmark is called,
!  it returns memory allocated after the last call to mmark.
!  mmark and retmark may be used multiply like parentheses, i.e.
!     call mmark
!       call getmem ( ... array1..)
!       call getmem (.. array 2..)
!       call getmem (... array 3)
!       call mmark
!         call getmem(... array 4)
!         call getmem (...array 5)
!       call retmark !returns array locations 4 and 5 to the memory pool
!     call retmark ! returns arrays 1 ,2 and 3
!
!  The maximum depth of these mmark-retmark parentheses is given
!  by the parameter maxmark

    implicit none

    mark = mark + 1
    if( mark > 100 ) call nerror( 1, 'mmark', & 
                    'too many marks, max=100', mark, nreq )
    marks( mark ) = nreq

  end subroutine mmark

!================================================================

  subroutine retmark

!  see the comments at mmark. Returns memory to the last mark.
!  do not use this in inner loops
!  if the memory has been released by using retmem, this routine
!  will not do anything

    implicit none
    integer :: newallo

    if(mark.gt.0) then
      newallo = nreq - marks(mark)
      if( newallo .gt. 0) then
        call retmem( newallo )
      endif
      mark = mark-1
    endif

  end subroutine retmark

!================================================================

  subroutine memstat( nreque, nmark, lastadr, memtot, iceiling, ioffset)

!   this routine returns the memory status
!   it returns the current memory request number, i.e.
!   the number of allocated blocks, the number of current marks,
!   the last used address, and the total available memory

    implicit none

    integer, intent( out ) :: nreque, nmark, lastadr, memtot, & 
                              iceiling, ioffset
      nreque  = nreq
      nmark   = mark
      if( nreq .gt. maxall) return
      lastadr  = nadr(nreq+1)
      memtot   = lcore
      iceiling = maxmem
      ioffset  = nadr(1) - 1

  end subroutine memstat

end module memory

subroutine int_to_bl( iarray, ipos, ival )

  implicit none
  integer :: iarray(*)
  integer :: ipos, ival

  iarray(ipos) = ival

end subroutine int_to_bl

subroutine int_from_bl( iarray, ipos, ival )

  implicit none
  integer :: iarray(*)
  integer :: ipos, ival

  ival = iarray(ipos)

end subroutine int_from_bl
