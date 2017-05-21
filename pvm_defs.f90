
               ! pvm definitions

module pvm_defs

  use kinds
  implicit none
  private

!
! $Id: fpvm3.h,v 1.15 2000/02/10 23:53:08 pvmsrc Exp $
!
!  -------------------------------------------------------------------
!          PVM version 3.4:  Parallel Virtual Machine System
!                University of Tennessee, Knoxville TN.
!            Oak Ridge National Laboratory, Oak Ridge TN.
!                    Emory University, Atlanta GA.
!       Authors:  J. J. Dongarra, G. E. Fagg, M. Fischer
!           G. A. Geist, J. A. Kohl, R. J. Manchek, P. Mucci,
!          P. M. Papadopoulos, S. L. Scott, and V. S. Sunderam
!                    (C) 1997 All Rights Reserved
!
!                               NOTICE
!
!  Permission to use, copy, modify, and distribute this software and
!  its documentation for any purpose and without fee is hereby granted
!  provided that the above copyright notice appear in all copies and
!  that both the copyright notice and this permission notice appear in
!  supporting documentation.
!
!  Neither the Institutions (Emory University, Oak Ridge National
!  Laboratory, and University of Tennessee) nor the Authors make any
!  representations about the suitability of this software for any
!  purpose.  This software is provided ``as is'' without express or
!  implied warranty.
!
!  PVM version 3 was funded in part by the U.S. Department of Energy,
!  the National Science Foundation and the State of Tennessee.
!  -------------------------------------------------------------------
!     ----------------------------------
!     fpvm3.h
!
!     Definitions to be included with
!     User Fortran application
!     ----------------------------------

! --------------------
! spawn 'flag' options
! --------------------
  integer(kind=i_4), public, parameter :: PVMTASKDEFAULT    =  0
  integer(kind=i_4), public, parameter :: PVMTASKHOST       =  1
  integer(kind=i_4), public, parameter :: PVMTASKARCH       =  2
  integer(kind=i_4), public, parameter :: PVMTASKDEBUG      =  4
  integer(kind=i_4), public, parameter :: PVMTASKTRACE      =  8
  integer(kind=i_4), public, parameter :: PVMMPPFRONT       = 16
  integer(kind=i_4), public, parameter :: PVMHOSTCOMPL      = 32
  integer(kind=i_4), public, parameter :: PVMNOSPAWNPARENT  = 64

! --------------------------------
! old option names still supported
! --------------------------------
  integer(kind=i_4), public, parameter :: PVMHOST  =  1
  integer(kind=i_4), public, parameter :: PVMARCH  =  2
  integer(kind=i_4), public, parameter :: PVMDEBUG =  4
  integer(kind=i_4), public, parameter :: PVMTRACE =  8

! -------------------------
! buffer 'encoding' options
! -------------------------
  integer(kind=i_4), public, parameter :: PVMDATADEFAULT = 0
  integer(kind=i_4), public, parameter :: PVMDATARAW     = 1
  integer(kind=i_4), public, parameter :: PVMDATAINPLACE = 2
  integer(kind=i_4), public, parameter :: PVMDATATRACE   = 4

! --------------------------------
! old option names still supported
! --------------------------------
  integer(kind=i_4), public, parameter :: PVMDEFAULT = 0
  integer(kind=i_4), public, parameter :: PVMRAW     = 1
  integer(kind=i_4), public, parameter :: PVMINPLACE = 2

! ----------------------
! notify 'about' options
! ----------------------
  integer(kind=i_4), public, parameter :: PVMTASKEXIT     = 1 
  integer(kind=i_4), public, parameter :: PVMHOSTDELETE   = 2 
  integer(kind=i_4), public, parameter :: PVMHOSTADD      = 3 
  integer(kind=i_4), public, parameter :: PVMROUTEADD     = 4 
  integer(kind=i_4), public, parameter :: PVMROUTEDELETE  = 5 
  integer(kind=i_4), public, parameter :: PVMNOTIFYCANCEL = 256

! --------------------------------
! packing/unpacking 'what' options
! --------------------------------
  integer(kind=i_4), public, parameter :: STRING   = 0
  integer(kind=i_4), public, parameter :: BYTE1    = 1
  integer(kind=i_4), public, parameter :: INTEGER2 = 2
  integer(kind=i_4), public, parameter :: INTEGER4 = 3
  integer(kind=i_4), public, parameter :: REAL4    = 4
  integer(kind=i_4), public, parameter :: COMPLEX8 = 5
  integer(kind=i_4), public, parameter :: REAL8    = 6
  integer(kind=i_4), public, parameter :: COMPLEX16= 7
  integer(kind=i_4), public, parameter :: INTEGER8 = 8

! --------------------------------
! setopt/getopt options for 'what'
! --------------------------------
  integer(kind=i_4), public, parameter :: PVMROUTE         = 1
  integer(kind=i_4), public, parameter :: PVMDEBUGMASK     = 2
  integer(kind=i_4), public, parameter :: PVMAUTOERR       = 3
  integer(kind=i_4), public, parameter :: PVMOUTPUTTID     = 4
  integer(kind=i_4), public, parameter :: PVMOUTPUTCODE    = 5
  integer(kind=i_4), public, parameter :: PVMTRACETID      = 6
  integer(kind=i_4), public, parameter :: PVMTRACECODE     = 7
  integer(kind=i_4), public, parameter :: PVMTRACEBUFFER   = 8
  integer(kind=i_4), public, parameter :: PVMTRACEOPTIONS  = 9
  integer(kind=i_4), public, parameter :: PVMFRAGSIZE      = 10
  integer(kind=i_4), public, parameter :: PVMRESVTIDS      = 11
  integer(kind=i_4), public, parameter :: PVMSOUTPUTTID    = 12
  integer(kind=i_4), public, parameter :: PVMSOUTPUTCODE   = 13
  integer(kind=i_4), public, parameter :: PVMSTRACETID     = 14
  integer(kind=i_4), public, parameter :: PVMSTRACECODE    = 15
  integer(kind=i_4), public, parameter :: PVMSTRACEBUFFER  = 16
  integer(kind=i_4), public, parameter :: PVMSTRACEOPTIONS = 17
  integer(kind=i_4), public, parameter :: PVMSHOWTIDS      = 18
  integer(kind=i_4), public, parameter :: PVMPOLLTYPE      = 19
  integer(kind=i_4), public, parameter :: PVMPOLLTIME      = 20
  integer(kind=i_4), public, parameter :: PVMOUTPUTCTX     = 21
  integer(kind=i_4), public, parameter :: PVMTRACECTX      = 22
  integer(kind=i_4), public, parameter :: PVMSOUTPUTCTX    = 23
  integer(kind=i_4), public, parameter :: PVMSTRACECTX     = 24
  integer(kind=i_4), public, parameter :: PVMNORESET       = 25

! --------------------------------------------
! tracing option values for setopt function
! --------------------------------------------
  integer(kind=i_4), public, parameter :: PVMTRACEFULL     = 1
  integer(kind=i_4), public, parameter :: PVMTRACETIME     = 2
  integer(kind=i_4), public, parameter :: PVMTRACECOUNT    = 3

! --------------------------------------------
! poll type options for 'how' in setopt function
! --------------------------------------------
  integer(kind=i_4), public, parameter :: PVMPOLLCONSTANT = 1
  integer(kind=i_4), public, parameter :: PVMPOLLSLEEP    = 2

! --------------------------------------------
! for message mailbox operations
! --------------------------------------------
  integer(kind=i_4), public, parameter :: PVMMBOXDEFAULT       =  0
  integer(kind=i_4), public, parameter :: PVMMBOXPERSISTENT    =  1
  integer(kind=i_4), public, parameter :: PVMMBOXMULTIINSTANCE =  2
  integer(kind=i_4), public, parameter :: PVMMBOXOVERWRITABLE  =  4
  integer(kind=i_4), public, parameter :: PVMMBOXFIRSTAVAIL    =  8
  integer(kind=i_4), public, parameter :: PVMMBOXREADANDDELETE = 16
  integer(kind=i_4), public, parameter :: PVMMBOXWAITFORINFO   = 32

! --------------------------------------------
! routing options for 'how' in setopt function
! --------------------------------------------
  integer(kind=i_4), public, parameter :: PVMDONTROUTE  = 1
  integer(kind=i_4), public, parameter :: PVMALLOWDIRECT= 2
  integer(kind=i_4), public, parameter :: PVMROUTEDIRECT= 3

! --------------------------
! error 'info' return values
! --------------------------
  integer(kind=i_4), public, parameter :: PvmOk           =   0
  integer(kind=i_4), public, parameter :: PvmBadParam     =  -2
  integer(kind=i_4), public, parameter :: PvmMismatch     =  -3
  integer(kind=i_4), public, parameter :: PvmOverflow     =  -4
  integer(kind=i_4), public, parameter :: PvmNoData       =  -5
  integer(kind=i_4), public, parameter :: PvmNoHost       =  -6
  integer(kind=i_4), public, parameter :: PvmNoFile       =  -7
  integer(kind=i_4), public, parameter :: PvmDenied       =  -8
  integer(kind=i_4), public, parameter :: PvmNoMem        = -10
  integer(kind=i_4), public, parameter :: PvmBadMsg       = -12
  integer(kind=i_4), public, parameter :: PvmSysErr       = -14
  integer(kind=i_4), public, parameter :: PvmNoBuf        = -15
  integer(kind=i_4), public, parameter :: PvmNoSuchBuf    = -16
  integer(kind=i_4), public, parameter :: PvmNullGroup    = -17
  integer(kind=i_4), public, parameter :: PvmDupGroup     = -18
  integer(kind=i_4), public, parameter :: PvmNoGroup      = -19
  integer(kind=i_4), public, parameter :: PvmNotInGroup   = -20
  integer(kind=i_4), public, parameter :: PvmNoInst       = -21
  integer(kind=i_4), public, parameter :: PvmHostFail     = -22
  integer(kind=i_4), public, parameter :: PvmNoParent     = -23
  integer(kind=i_4), public, parameter :: PvmNotImpl      = -24
  integer(kind=i_4), public, parameter :: PvmDSysErr      = -25
  integer(kind=i_4), public, parameter :: PvmBadVersion   = -26
  integer(kind=i_4), public, parameter :: PvmOutOfRes     = -27
  integer(kind=i_4), public, parameter :: PvmDupHost      = -28
  integer(kind=i_4), public, parameter :: PvmCantStart    = -29
  integer(kind=i_4), public, parameter :: PvmAlready      = -30
  integer(kind=i_4), public, parameter :: PvmNoTask       = -31
  integer(kind=i_4), public, parameter :: PvmNotFound     = -32
  integer(kind=i_4), public, parameter :: PvmExists       = -33
  integer(kind=i_4), public, parameter :: PvmHostrNMstr   = -34
  integer(kind=i_4), public, parameter :: PvmParentNotSet = -35
  integer(kind=i_4), public, parameter :: PvmIPLoopback   = -36

! --------------------------
! these are going away in the next version.
! use the replacements
! --------------------------
  integer(kind=i_4), public, parameter :: PvmNoEntry    = -32
  integer(kind=i_4), public, parameter :: PvmDupEntry   = -33


! Uncomment this include for use with the WIN32 WATCOM fortran compiler
!
!     include '../include/fpvm3_watcom.h'
!


  integer, public :: PVM_NH  ! No. of hosts in the virtual machine
  character(len=79), public :: GROUP_NAME
  integer(kind=i_4), public ,allocatable :: DMNID(:)

  public :: gethstnm, para_gettids

 contains

  subroutine para_gettids

! setup list of slave tids (consistent on all processes)

    use newpara
    implicit none

    integer islv
    do islv=0,NSLV
       call pvmfgettid(GROUP_NAME,islv,SLVID(islv+1))
    end do
    
  end subroutine para_gettids
      
  subroutine gethstnm

!
!     get PVM host name for this machine
!

    use newpara
    implicit none

    integer*4 ihostid,narch,itd,ispeed,nh,info
    integer i
    character*20 arch
    call blankit(MY_HOSTNM,256)
    call pvmftidtohost(MY_ID,ihostid)
!
!     loop until (ihostid.eq.itd) which means we have found the hostname for
!     this machine .or. until we have run out of hosts (i.gt.PVM_NH).
!
    call pvmfconfig(PVM_NH,narch,itd,MY_HOSTNM,arch,ispeed,info)
    i=2
    do while((itd.ne.ihostid).and.(i.le.PVM_NH))
      call pvmfconfig(nh,narch,itd,MY_HOSTNM,arch,ispeed,info)
      i=i+1
    end do
    
  end subroutine gethstnm

end module pvm_defs
