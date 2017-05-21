module messages

! Module for composing and displaying messages
!
! this is a system dependent code. On linux systems
! the message is written to a fortran output unit,
! so the whole thing is mostly redundant. On the Windows
! intel version, however, the message is displayed on
! a message box.
!
! The message to be displayed is stored in a message buffer, wbuf. 
! Routines are provided to add a line of text to the message buffer
! and to display the contents of the message buffer. There is also
! a routine for displaying the contents of an I/O channel that works
! only on Windows
!
! addbuf(string)     adds string as a line to the message buffer
!
! print_buf(n,title)                        print buffer on I/O unit n.
!                                           using the given title line.

! display_buf(n,title) on Linux:            calls print_buf
!                      on Windows (intel) : display the buffer on a 
!                                           message box with the given title.
!                                           N is ignored
!
! cat_channel(ni,title) on LInux:           do nothing
!                       on Windows (intel): display the contents of I/O
!                                           unit ni on a message box
!                                           with the given title.
!                                           The content of the unit is
!                                           deleted.

  implicit none
  private

  integer, parameter :: lenbuf=10000
  character(len=lenbuf), public :: wbuf

  public :: addbuf, print_buf, display_buf, cat_channel

  contains

  subroutine addbuf(string)
  
    implicit none
    character(len=*) :: string
    integer leni,lens
    
    leni=len_trim(wbuf)
    lens=len_trim(string)
    if(leni+lens+2.le.lenbuf)then
      if(leni.gt.0)then
#ifdef __INTEL_COMPILER

            ! using C strings extension of the Intel compiler

        wbuf(leni+1:leni+1)="\n"C
        wbuf(leni+2:leni+lens+2)=string(1:len_trim(string))
#else
        wbuf(leni+1:leni+lens+2)="\n"//string(1:len_trim(string))
#endif
      else
        wbuf=trim(string)
      endif
    endif
  
  end subroutine addbuf

  subroutine print_buf(n,title)

    use kinds

    implicit none
    integer n
    character(len=*) :: title
    integer leni

    leni=len_trim(wbuf)
    if(leni.gt.0)then
      write(n,'(/a/)')title
      write(n,'(a)')trim(wbuf)
    endif

  end subroutine print_buf

  subroutine display_buf(n,title)

#ifdef WINDOWS
#ifdef  __INTEL_COMPILER
      use IFQWIN ! Quick windows libraries
#endif
#endif
    use kinds

    implicit none
    integer n
    character(len=*) :: title
    integer(kind=i_4) :: iresp

#ifdef WINDOWS

#ifdef  __INTEL_COMPILER
    iresp=MESSAGEBOXQQ(trim(wbuf)//''C,trim(title)//''C,MB$OK)
#else
    call print_buf(n,title)
#endif

#else
    call print_buf(n,title)
#endif

  end subroutine display_buf

  subroutine cat_channel( ni, title )

    implicit none
    integer ni
    character(len=*) :: title
    character(len=256) :: string

#ifdef WINDOWS
#ifdef  __INTEL_COMPILER
    rewind(ni)
    wbuf=''
    do
      string=''
      read(ni,'(a)',end=99)string
      call addbuf(string)
    enddo

    99 rewind(ni)
    endfile(ni)

    call display_buf(6,title)
#endif
#endif

    return

  end subroutine cat_channel

end Module messages
