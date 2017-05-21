!
! collection of (potentially) system-dependent  utility function
! and subroutines.
! 
! the idea is to have local version of these functions in case
! a particular hardware or operating system needs it
!
!================================================================

subroutine printexetype( exename, iout )

         ! This subroutine prints out information on the
         ! executable currently running:
         !
         ! 1) the type of executable we are currently running
         !    ( 32- or 64-bit, ecc)
         !
         ! 2) the size of the default integer type:
         !   intsize = 1 for 64-bit integers
         !   intsize = 2 for 32-bit integers

         ! uses the nonstandard routine system, also depends
         ! on the unix command "file" and windows command "objdump"

#ifdef __INTEL_COMPILER
  use ifport
#endif

  use kinds

  implicit none
  integer, intent(in) :: iout
  character, intent(in) :: exename*256
  integer, parameter :: lenout = 58
  character(len=lenout) :: line
  character :: freeout*256, a*128, pre*14
  integer :: i, ic, il, len, info, nbits, lenf, ierr

  nbits = 0

#ifdef WINDOWS
         ! no executable type on Windows
#else
         !  generate a file name in freeout

  call getchval('jobname',freeout) 
  lenf=len_trim(freeout)
  freeout(lenf+1:lenf+5)='.finf'

         ! execute the "file" command, writing the output on the
         ! freeout file

#ifdef __INTEL_COMPILER
  ierr = system( 'file -b '//trim(exename)//' > "'//trim(freeout)//'"' )
#else
  ierr = 0
  call system( 'file -b '//trim(exename)//' > "'//trim(freeout)//'"' )
#endif

         ! read in freeout

  if( ierr .ge. 0 ) then
    open( 1, file = trim(freeout), status = 'unknown', err = 1000 )
    a = ''
    read( 1,'(a)', err = 1000 ) a
    close( 1, err = 1000, status = 'delete' )

           ! Print the type information with some primitive
           ! formatting (indentation and line wrapping)
           
    len = len_trim( a )
    pre = ' Type       : '

    il = 0
    ic = 0
    i = index( a, ':', back=.true. )
    do
      if ( i .ge. len ) exit
      i = i + 1
      il = il + 1
      if ( il .gt. lenout ) then
        if( ic .ne. 0 ) then
          write( iout, '(a,a)') pre, adjustl(line(1:ic))
          pre = ''
          if( ic .lt. lenout )line(1:lenout-ic) = line(ic+1:lenout)
          il = lenout - ic + 1
        else
          write( iout, '(a,a)') pre, adjustl(line(1:lenout))
          pre = ''
          il = 1
        endif
        ic = 0
      endif
      line(il:il) = a(i:i)
      if( line(il:il) .eq. ',' ) ic = il
    enddo
    if(il .ne. 0 ) &
       write( iout, '(a,a)' ) pre, adjustl(line(1:min(il,lenout)))

           ! find out if we are compiled 32-bit or 64-bit

    i = index( a, '-bit' )
    if( i .le. 2 ) i = index( a, '-BIT' )
    if( i .gt. 2 ) then
      read( a(i-2:i-1), '(i2)', err=1000 ) nbits
    endif
  endif

1000 continue
#endif

         ! print intsize

  write( iout, '('' Intsize    : '',i2)' ) intsize

         ! store nbits into depository

  call setival( 'nbits', nbits )

end subroutine printexetype

!================================================================

subroutine init_random( iseed )

           ! initializes the random number generator.
           ! iseed is the seed value. If it is zero, a new
           ! seed will be computed based on the process id
           ! or on the system date, depending on the compiler
           ! used

#ifdef __INTEL_COMPILER
  use ifport
#endif

  use kinds

  implicit none
  integer, intent(inout) :: iseed
  integer, allocatable :: seeda(:)
  integer :: nseed = 1, i, ipid
#ifdef __INTEL_COMPILER
 !integer(i_4), external :: getpid
#else
  character*24 date
  integer ibuf
#endif

#ifdef __INTEL_COMPILER
  if( iseed .eq. 0 ) iseed = getpid()   ! get process id
#else
  if( iseed .eq. 0 ) then

          ! do not know how to get pid on Mac, so use date and time instead

    call date1( date ) ! returns the date as in 'Wed Oct  5 18:36:02 2005' 

    read(date,'(8x,i2)') ibuf  ! days
    ipid = ibuf
    read(date,'(11x,i2)') ibuf  ! hours
    ipid = ipid + ibuf
    read(date,'(14x,i2)') ibuf  ! minutes
    ipid = ipid + ibuf
    read(date,'(17x,i2)') ibuf  ! seconds
    ipid = ipid + ibuf
    read(date,'(20x,i4)') ibuf  ! years
    ipid = ipid + ibuf
    iseed = ipid

  endif
#endif

  call random_seed( size = nseed ) ! get size of seed
  allocate( seeda(nseed) )
  seeda = iseed 

  call random_seed( put = seeda ) ! initialize random generator

  deallocate( seeda )
       
end subroutine init_random

!================================================================

subroutine checkfile(check)

!  This routine copies the files of a previous run to the current
!  filenames, enabling thus the checkpoint option
!  uses the nonstandard routine system and the unix commands
!  cp and tar

#ifdef __INTEL_COMPILER
  use ifport
#endif

  character check*256,jobname*256
  logical ex,exx
  parameter(nexten=10)
  dimension  Len(nexten)
  character*8 ext(nexten)
  integer ierr
  data ext/'.coord  ','.control','.basis  ','.hess   ','.mos    ', &
           '.mob    ','.sym    ','.grad   ','.deriv  ','.zmat   '/
  data Len/6,8,6,5,4,4,4,5,6,5/

  iout=igetival('iout')
  call rmblan2(check,256,LenC)
!   write message
  write(iout,*) 'Check file used=',check(1:LenC)
  call getchval('jobname',jobname)
  call rmblan2(jobname,256,LenJ)
  inquire(file=check(1:LenC)//ext(1)(1:Len(1)),exist=exx)
  if(exx) then
    write(iout,*) check(1:LenC)//'.coord',' file found; ', &
             'copying it and other check files to ',jobname(1:LenJ)
    do i=1,nexten
      inquire(file=check(1:LenC)//ext(i)(1:Len(i)),exist=ex)
      if(ex) then

#ifdef WINDOWS

#ifdef __INTEL_COMPILER
        ierr = system('copy /y "'//check(1:LenC)//ext(i)(1:Len(i)) &
             //'"  "'//jobname(1:LenJ)//ext(i)(1:Len(i))//'" > nul')
#else
        call system('copy /y "'//check(1:LenC)//ext(i)(1:Len(i)) &
             //'"  "'//jobname(1:LenJ)//ext(i)(1:Len(i))//'" > nul')
#endif

#else 

#ifdef __INTEL_COMPILER
        ierr = system('cp -f "'//check(1:LenC)//ext(i)(1:Len(i)) &
             //'" "'//jobname(1:LenJ)//ext(i)(1:Len(i))//'"')
#else
        call system('cp -f "'//check(1:LenC)//ext(i)(1:Len(i)) &
             //'" "'//jobname(1:LenJ)//ext(i)(1:Len(i))//'"')
#endif

#endif

      end if
    end do
  else

#ifdef WINDOWS
    write(iout,*) 'Could not locate checkpoint files ',check(1:LenC)
#else
    write(iout,*) 'File ',check(1:LenC)//'.coord not found'
    write(iout,*) 'Try to reconstruct check files from the ', &
                  'compressed archive'
    inquire(file=check(1:LenC)//'.tgz',exist=exx)
    if(exx) then
      write(iout,*) 'Archive file ',check(1:LenC)//'.tgz',' found', &
                   '; decompressing and renaming check files'
#ifdef __INTEL_COMPILER
      ierr = system('tar -xvzf "'//check(1:LenC)//'.tgz" > /dev/null')
#else
      call system('tar -xvzf "'//check(1:LenC)//'.tgz" > /dev/null')
#endif

      do i=1,nexten
        inquire(file=check(1:LenC)//ext(i)(1:Len(i)),exist=ex)
        if(ex) then

#ifdef __INTEL_COMPILER
          ierr =  system('cp -f "'//check(1:LenC)//ext(i)(1:Len(i)) &
               //'" "'//jobname(1:LenJ)//ext(i)(1:Len(i))//'"')
#else
          call system('cp -f "'//check(1:LenC)//ext(i)(1:Len(i)) &
               //'" "'//jobname(1:LenJ)//ext(i)(1:Len(i))//'"')
#endif

        end if
      end do
    else
      write(iout,*) 'Could not locate checkpoint or archive files ',check(1:LenC)
    end if
#endif

  end if

end subroutine checkfile

!=======================================================================

subroutine runchelp(jobname, lenj, natoms, atsymb, iout, icon)

  use sysdef

#ifdef __INTEL_COMPILER
  use ifport

#ifdef WINDOWS
  use kernel32, only: CreateProcess, WaitForsingleObject, CloseHandle
  use ifwinty
#endif

#endif

  implicit none

!  CHELP charges. This routine runs the external program "pqspoint" in
!  order to compute the chelp charges, then reads the results from file
!  jobname.chelp and prints out the chelp charges and dipole moment.
!  usees the nonstandard routine system. Paths and directory separators
!  are also involved

  character*256 jobname,field
  character*8 atsymb(natoms)
  character*8 readat
  character*9 nulla
  character*256 path,args
  character*256 pqsroot
  character*500 cmdline
  logical found,isthere
  integer lenj,natoms,iout,icon
  integer lenp,ierr,i
  real*8 charge
#ifdef WINDOWS
#ifdef __INTEL_COMPILER

        !  declarations needed for call to Windows API function CreateProcess

  type(T_STARTUPINFO) :: sinfo
  integer :: sinfosize=9*DWORD+3*LPSTR+2*WORD+LPBYTE+3*HANDLE
  type(T_PROCESS_INFORMATION) :: pinfo
  integer :: pinfosize=2*HANDLE+2*DWORD

#endif
#endif

!  print banner

  write(iout,1000)
  write(icon,1000)

!  to locate the pqspoint executable we try:
!
!  1) $PQS_ROOT/UTILS
!  2) $PQS_ROOT
!
! if pqspoint is not found in any of these locations, we give up

  call getchval('pqs_root',pqsroot)
  call rmblan2(pqsroot,256,lenp)
  found=.false.

!  1) $PQS_ROOT/UTILS

#ifdef WINDOWS
  path='\UTILS\pqspoint.exe'
#else
  path='/UTILS/pqspoint'
#endif
  inquire( file=pqsroot(1:lenp)//path(1:len_trim(path)),exist=isthere)
  if(isthere) found=.true.
  if(.not.found)then

!  2) $PQS_ROOT

#ifdef WINDOWS
    path='\pqspoint.exe'
#else
    path='/pqspoint'
#endif
    inquire( file=path(1:lenp)//path(1:len_trim(path)),exist=isthere)
    if(isthere) found=.true.
  endif
  if(.not.found)then

!  we give up

    write(iout,*)'***WARNING: pqspoint executable not found in'
    write(iout,*)'       '//pqsroot(1:lenp)//DIR_SEP//'UTILS'
    write(iout,*)'       '//pqsroot(1:lenp)
    write(iout,*)
    write(iout,*)'        CHELP calculation aborted.'
    write(iout,*)
    return
  endif

!  pqspoint has been located. Now we run the command
!  pqspoint jobname -xp 0 -yp 0 -zp 0 -thresh 1e-8 -chelp -end

#ifdef WINDOWS
  nulla = 'nul'
#else
  nulla = '/dev/null'
#endif

  cmdline = '"'//pqsroot(1:lenp)//path(1:len_trim(path)) &
            //'" "'//jobname(1:lenj)//'" -chelp -end  > '//nulla

#ifdef __INTEL_COMPILER

#ifdef WINDOWS
    
      ! this code uses the WINDOWS API

   sinfo=T_STARTUPINFO(sinfosize,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   pinfo=T_PROCESS_INFORMATION(0,0,0,0)
   ierr = CreateProcess(        &
            NULL,               & ! no module line (use command line)
            trim(cmdline)//''C, & ! command line
            NULL,               & ! process security attributes
            NULL,               & ! thread security attributes
            FALSE,              & ! do not inherit handles
            0,                  & ! no creation flags
            NULL,               & ! use parent's environment
            NULL,               & ! use parent's current directory
            sinfo,              & ! startup info
            pinfo               & ! process info
          )

          ! Wait until child process exits

  if( ierr .eq. int(FALSE) ) goto 666
  ierr = WaitForSingleObject( pinfo%hProcess, INFINITE )

          ! Close process and thread handles

  ierr = CloseHandle( pinfo%hProcess )
  ierr = CloseHandle( pinfo%hThread )

#else
  ierr =  system( trim(cmdline) )
#endif

#else
  call system( trim(cmdline) )
#endif

!  now read the chelp file

  open( unit=1, file=jobname(1:lenj)//'.chelp', status='old', form='formatted', err=666 )

!  read file until the string 'Atom' is located

  10 read(1,'(a)',end=666)field
    field=adjustl(field)
    if(field(1:1).ne.'a'.and.field(1:1).ne.'A')goto 10
    if(field(2:2).ne.'t'.and.field(2:2).ne.'T')goto 10
    if(field(3:3).ne.'o'.and.field(3:3).ne.'O')goto 10
    if(field(4:4).ne.'m'.and.field(4:4).ne.'M')goto 10

!  read and print chelp charges

  read(1,*,end=666)
  do i = 1, natoms
    read(1,'(a8,f20.6)',end=666)readat,charge
    write(iout,1100),i,atsymb(i),charge
    write(icon,1100),i,atsymb(i),charge
  enddo

!  read and print chelp dipole moment

  read(1,*,end=666)
  read(1,*,end=666)
  read(1,*,end=666)
  read(1,'(a)',end=666)field
  call rmblan(field,256,lenp)
  write(iout,'(/a)') field(1:lenp)
  write(icon,'(/a)') field(1:lenp)

  close(unit=1)
  return

  666 write(iout,*)'***WARNING: '//jobname(1:lenj)//'.chelp'
      write(iout,*)'           the file does not exist or is incomplete'
      write(iout,*)
      write(iout,*)'            CHELP calculation aborted.'
      write(iout,*)
      close(unit=1)

  1000 format(/,' ---- CHELP charges ----',/)
  1100 format(1X,I5,1X,A8,2X,F12.6)

end subroutine runchelp

!=======================================================================

subroutine set_pqs_root

! set the pqs_root variable, either from the environment or from the
! defaults. The default values are system dependent

  implicit none

  character(len=256) :: pqs_root
  integer :: len
  logical :: isthere

     !   1) try the environmental variable

  pqs_root=''
  call get_environment_variable('PQS_ROOT',pqs_root)
  call rmblan(pqs_root,256,len)

  if(len.eq.0) then

           !   2) try the default values

#ifdef WINDOWS
   pqs_root='c:\Program Files\PQS\PQS 3.3'
   inquire(file=pqs_root(1:28)//'\PQSv33.exe', EXIST=isthere)
#else
   pqs_root='/usr/local/share/PQS'
   inquire(file=pqs_root(1:20)//'/pqs.x', EXIST=isthere)
#endif

    if(.not.isthere) then

#ifdef WINDOWS
     pqs_root='c:\PQS\PQS 3.3'
     inquire(file=pqs_root(1:14)//'\PQSv33.exe', EXIST=isthere)
#else
     pqs_root='/usr/local/PQS'
     inquire(file=pqs_root(1:14)//'/pqs.x', EXIST=isthere)
#endif

     if(.not.isthere) then
       call nerror(1,'set_pqs_root','PQS_ROOT could not be located',0,0)
     end if

    endif

  endif
  call setchval('pqs_root',pqs_root)

end subroutine set_pqs_root

!=======================================================================

subroutine set_scrdir

! set the scrdir variable, either from the environment or from the
! defaults. The default values are system dependent

  use sysdef
  implicit none

  character(len=256) :: scrdir,user
  integer :: len

     !   1) try the environmental variable

  scrdir=''
  call get_environment_variable('PQS_SCRDIR',scrdir)
  call rmblan(scrdir,256,len)

  if(len.eq.0) then

           !   2) use the default value

#ifdef WINDOWS
    call get_environment_variable('USERPROFILE',user)
    call rmblan(user,256,len)
    scrdir=user(1:len)//'\Local Settings\Temp'
#else
    call get_environment_variable('USER',user)
    call rmblan(user,256,len)
    scrdir='/scr/'//user(1:len)
#endif

  endif

  ! if needed, append a directory separator at the end of scrdir

  call rmblan(scrdir,256,len)
  if(scrdir(len:len).ne. DIR_SEP) scrdir(len+1:len+1)=DIR_SEP

  call setchval('scrdir',scrdir)

end subroutine set_scrdir

!=======================================================================

subroutine get_hostname( hname )

! get the hostname.
!
! uses environmental variables instead of calling non standard and not
! much portable functions

  implicit none
  character(len=256) :: hname

  hname = ''
  call get_environment_variable('COMPUTERNAME',hname) ! this should work on Windows
  if ( len_trim( hname ) .eq. 0 ) call get_environment_variable('HOSTNAME',hname)
  if ( len_trim( hname ) .eq. 0 ) call get_environment_variable('HOST',hname)
  if ( len_trim( hname ) .eq. 0 ) hname = 'unknown'  ! we give up

end subroutine get_hostname

!=======================================================================

subroutine f_lush( ichan )

! try to flush the I/O channel ichan
! calls the utility routine flush which is non standard
! the funny name is for fooling grep

#ifdef __INTEL_COMPILER
  use ifport
#endif

  use kinds

  implicit none
  integer :: ichan
  integer(kind=i_4) :: i4chan

  i4chan=ichan
  call flush( i4chan )

end subroutine f_lush

!=======================================================================

subroutine open_tim

! handles the opening of Krzysztof's timing file (unit 91).
! Now by default this file will be connected to the
! "big bit bucket in the sky" (/dev/null). It is possible
! to instruct the program to connect the unit to the <jobname>.tim
! file by setting the depository integer variable 'tim' to a value
! greater than zero. 
! This can be activated by the command line option -tim
!
! Now the tim file is redirected  to /dev/null also in the single processor
! version in preparation for the parallel thread-based version, to avoid
! having different threads writing to the same unit at the same time.
!
! A system dependent routine is needed because /dev/null works
! only for Linux. For windows the syntax is different.
!
  implicit none
  character*256 jobname
  integer itim, isthere

                         !  test for 'tim' depository variable

  call tstival('tim',isthere)
  if(isthere.ne.0)then
    call getival('tim',itim)
  else
    itim=0
  endif
  if(itim.eq.0)then         ! send output to nowhere

#ifdef WINDOWS
      open(91,file='nul')
#else
      open(91,file='/dev/null')
#endif

  else                      ! send output to <jobname>.tim

    call getchval('jobname',jobname)
    open(unit=91,file=jobname(1:len_trim(jobname))//'.tim',status='unknown')

  endif

end subroutine open_tim

!=======================================================================

subroutine set_pqs_env

! sets the environmental variable PQS_ROOT. this is used by the slaves,
! to avoid potential problem with license checking if the environment
! is not properly set.

#ifdef __INTEL_COMPILER
  use ifport
#endif
  implicit none
  character*254 pqsroot,name
  logical success
  integer iret, pqs_setenv

  name = ''
  pqsroot=''
  call getchval('pqs_root',pqsroot)

#ifdef __INTEL_COMPILER
  success = setenvqq('PQS_ROOT='//pqsroot(1:len_trim(pqsroot)))
  if (.not. success ) call message( 'set_pqs_env', 'WARNING: cannot set PQS_ROOT', 0, 0)
#else
  name = 'PQS_ROOT'
  iret = pqs_setenv( name, len_trim(name), pqsroot, len_trim(pqsroot) )
  if (iret .ne. 0 ) call message( 'set_pqs_env', 'WARNING: cannot set PQS_ROOT', 0, 0)
#endif

end subroutine set_pqs_env
