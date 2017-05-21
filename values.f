c Aug 19 PP added the following functionality:
c  integer function   igetival(name)
c  real*8 function    rgetrval(name)
c  character*256 function chgetchval(name)
c  logical functions ltstival(name), ltstrval(name), ltstchval(name)
cJUN10,93 PP the value names have been changed to 8-character words
cJUN10,93 the names have been changed from 'options' to 'values'
      subroutine valinit
c  this routine has to be called first before the value subsystem
c  is called
c
      implicit real*8 (a-h,o-z)
      character*8 names,rnames,cnames
      character*256 cvals
      parameter (nvali=1000,nvalr=500,nvalc=100)
      common /ivalns/icur,ivals(nvali),names(nvali)
      common /rvalns/ircur,jcur,rvals(nvalr),rnames(nvalr)
      common /chvalns/iccur,cvals(nvalc),cnames(nvalc)
      save ifirst
      data ifirst/0/
      if(ifirst.ne.0) return
      icur=0
      ircur=0
      iccur=0
      ifirst=1
      end
c========================================================================
      subroutine getival(namei,ival)
c  this routine recovers an integer number which corresponds to the
c  8-character name "namei" in ival. The routine will fail and call
c  nerror if the value has not been yet defined. If you are uncertain,
c  call first testival
c  example for calling
c      call getival('na',na)
c
      character*(*) namei
      character*8 name,names
      parameter (nval=1000)
      common /ivalns/icur,ivals(nval),names(nval)
      name=namei
      do i=1,icur
         if (name.eq.names(i)) then
            ival=ivals(i)
            return
         end if
      end do
      call nerror(1,'getival','no such value exists '//name,i,i)
      end
c======================================================================
      subroutine tstival(namei,iyes)
c  this routine tests whether the integer value with the 4-character
c  name "namei" hasd been defined earlier. If yes, it returns 1 in iyes,
c  otherwise 0.
c  example for using it:
c      call tstival('na',ipresent)
c      if(ipresent.eq.0) then
c        write(*,*) 'number of atoms undefined
c        go to ....
c      else
c        call getival('na',na)
c      end if
      character*(*) namei
      character*8 name,names
      parameter (nval=1000)
      common /ivalns/icur,ivals(nval),names(nval)
      name=namei
      iyes=0
      do 100 i=1,icur
      if(name.eq.names(i)) then
        iyes=1
       return
      end if
 100  continue
      end
c=======================================================================
      subroutine setival(namei,ival)
c  this routine sets the value of the value with the name namei
c  to the value ival (an integer). If the value was undefined before,
c  it will be defined after this call. If it was defined, the value
c  changes to the new value. It is very important not to redefine
c  accidentally an value. Therefore, consult the file opvalues before
c  introducing a new name
c  example for calling
c       call setival('ncf',ncf)
c
      character*(*) namei
      character*8 name,names
      parameter (nval=1000)
      common /ivalns/icur,ivals(nval),names(nval)
      name=namei
      do 100 i=1,icur
      if(name.eq.names(i)) then
        ivals(i)=ival
       return
      end if
 100  continue
      icur=icur+1
      if(icur.gt.nval) then
      call nerror(1,'setival','too many values, max=',nval,icur)
      else
      names(icur)=name
      ivals(icur)=ival
      end if
      end
c=======================================================================
      subroutine getrval(namei,rval)
c  this is the real equivalent of getival
c
      character*(*) namei
      real*8 rval,rvals
      character*8 name,rnames
      parameter (nval=500)
      common /rvalns/ircur,jcur,rvals(nval),rnames(nval)
      name=namei
      do 100 i=1,ircur
      if(name.eq.rnames(i)) then
        rval=rvals(i)
       return
      end if
 100  continue
      call nerror(1,'getrval','no such value exists '//name,i,i)
      end
c========================================================================
      subroutine tstrval(namei,iyes)
c  this is the real equivalent of tstival
c
      character*(*) namei
      real*8 rvals
      character*8 name,rnames
      parameter (nval=500)
      common /rvalns/ircur,jcur,rvals(nval),rnames(nval)
      name=namei
      iyes=0
      do 100 i=1,ircur
      if(name.eq.rnames(i)) then
        iyes=1
       return
      end if
 100  continue
      end
c=======================================================================
      subroutine setrval(namei,rval)
c  this is the real equivalent of setival
c
      character*(*) namei
      real*8 rval,rvals
      character*8 name,rnames
      parameter (nval=500)
      common /rvalns/ircur,jcur,rvals(nval),rnames(nval)
      name=namei
      do 100 i=1,ircur
      if(name.eq.rnames(i)) then
        rvals(i)=rval
       return
      end if
 100  continue
      ircur=ircur+1
      if(ircur.gt.nval) then
      call nerror(1,'setrval','too many values, max=',nval,ircur)
      else
      rnames(ircur)=name
      rvals(ircur)=rval
      end if
      end
c=======================================================================
      subroutine getchval(namei,cval)
c  this is the character equivalent of getival
c  the string length stored is limited to 256 characters
c
      character*(*) namei,cval
      character*256 cvals
      character*8 name,cnames
      parameter (nval=100)
      common /chvalns/iccur,cvals(nval),cnames(nval)
      name=namei
      do 100 i=1,iccur
      if(name.eq.cnames(i)) then
        cval=cvals(i)
       return
      end if
 100  continue
      call nerror(1,'getchval','no such value exists '//name,i,i)
      end
c========================================================================
      subroutine tstchval(namei,iyes)
c  this is the character equivalent of tstival
c
      character*(*) namei
      character*8 name,cnames
      character*256 cvals
      parameter (nval=100)
      common /chvalns/iccur,cvals(nval),cnames(nval)
      name=namei
      iyes=0
      do 100 i=1,iccur
      if(name.eq.cnames(i)) then
        iyes=1
       return
      end if
 100  continue
      end
c=======================================================================
      subroutine setchval(namei,cvali)
c  this is the character equivalent of setival
c  cval is limited to 256 characters
c
      character*(*) namei,cvali
      character*256 cval,cvals
      character*8 name,cnames
      parameter (nval=100)
      common /chvalns/iccur,cvals(nval),cnames(nval)
      name=namei
      cval=cvali
      do 100 i=1,iccur
      if(name.eq.cnames(i)) then
        cvals(i)=cval
       return
      end if
 100  continue
      iccur=iccur+1
      if(iccur.gt.nval) then
      call nerror(1,'setchval','too many values, max=',nval,iccur)
      else
      cnames(iccur)=name
      cvals(iccur)=cval
      end if
      end
c========================================================================
      integer function igetival(namei)
c  this routine recovers an integer number which corresponds to the
c  8-character name "namei" in ival. The routine will fail and call
c  nerror if the value has not been yet defined. If you are uncertain,
c  call first testival
c  example for calling
c      na=igetival('na')
c
      character*(*) namei
      character*8 name,names
      parameter (nval=1000)
      common /ivalns/icur,ivals(nval),names(nval)
      name=namei
      do 100 i=1,icur
      if(name.eq.names(i)) then
        igetival=ivals(i)
       return
      end if
 100  continue
      call nerror(1,'getival','no such value exists '//name,i,i)
      end
c========================================================================
      logical function ltstival(namei)
c  this routine tests whether the integer value with the 4-character
c  name "namei" hasd been defined earlier. If yes, it returns 1 in iyes,
c  otherwise 0.
c  example for using it:
c      logical present
c      present=ltstival('na')
c      if(present) then
c        write(*,*) 'number of atoms undefined
c        go to ....
c      else
c        call getival('na',na)
c      end if
      character*(*) namei
      character*8 name,names
      parameter (nval=1000)
      common /ivalns/icur,ivals(nval),names(nval)
      name=namei
      ltstival=.false.
      do 100 i=1,icur
      if(name.eq.names(i)) then
        ltstival=.true.
       return
      end if
 100  continue
      end
c========================================================================
      real*8 function rgetrval(namei)
c  this is the real equivalent of igetival
c
      character*(*) namei
      real*8 rvals
      character*8 name,rnames
      parameter (nval=500)
      common /rvalns/ircur,jcur,rvals(nval),rnames(nval)
      name=namei
      do 100 i=1,ircur
      if(name.eq.rnames(i)) then
        rgetrval=rvals(i)
       return
      end if
 100  continue
      call nerror(1,'getrval','no such value exists '//name,i,i)
      end
c====================================================================
      logical function ltstrval(namei,iyes)
c  this is the real equivalent of tstival
c
      character*(*) namei
      real*8 rvals
      character*8 name,rnames
      parameter (nval=500)
      common /rvalns/ircur,jcur,rvals(nval),rnames(nval)
      name=namei
      ltstrval=.false.
      do 100 i=1,ircur
      if(name.eq.rnames(i)) then
        ltstrval=.true.
       return
      end if
 100  continue
      end
c=======================================================================
      character*256 function chgetchval(namei,cval)
c  this is the character equivalent of getival
c  the string length stored is limited to 256 characters
c
      character*(*) namei
      character*256 cvals
      character*8 name,cnames
      parameter (nval=100)
      common /chvalns/iccur,cvals(nval),cnames(nval)
      name=namei
      do 100 i=1,iccur
      if(name.eq.cnames(i)) then
        chgetchval=cvals(i)
        return
      end if
 100  continue
      call nerror(1,'getchval','no such value exists '//name,i,i)
      end
c========================================================================
      logical function ltstchval(namei)
c  this is the character equivalent of tstival
c
      character*(*) namei
      character*8 name,cnames
      character*256 cvals
      parameter (nval=100)
      common /chvalns/iccur,cvals(nval),cnames(nval)
      name=namei
      ltstchval=.false.
      do 100 i=1,iccur
      if(name.eq.cnames(i)) then
        ltstchval=.true.
        return
      end if
 100  continue
      end
c=======================================================================

