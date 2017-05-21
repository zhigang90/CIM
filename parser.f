      subroutine inpcopy(iraw,inew)
c
c     this subroutine massages the raw input in a form acceptable
c     for the later program parts. It is important to copy the input to
c     a new files because of the possibility of backspaces
c    
c
      implicit none

      integer iraw, inew
                                      ! the width of the input line
      integer, parameter :: linew=300 ! is set to 300, because on
                                      ! windows paths can get quite long
      character(len=linew) :: input
      integer ios

      do
        input=''
        read(iraw,'(a)',end=200,iostat=ios)input
        call writeline(input,inew,linew)
      enddo

 200  call writeline(input,inew,linew)

      end

c=====================================================================

      subroutine writeline(input,ichan,linew)
c
c     this routine writes the contents of input
c     into channel ichan.
c
c     input lines beginning with ! or ? are skipped, and
c     the text following an unquoted exclamation mark is
c     removed. An exclamation mark is considered quoted if it
c     is found between a pair of quote characters (',`,").
c     This is to allow the use of exclamation marks as 
c     part of file names.
c
c     linew is the maximum width of an imput line
c
      implicit none

      integer linew,ichan
      character(len=linew) :: input,input1
      character(len=1), parameter:: excl='!', ques='?'
      character(len=3), parameter :: quote='"''`'
      integer i,iq,iqq,len
      logical quoted

      len=len_trim(input)

      if(input(1:1).eq.excl.or.input(1:1).eq.ques) return ! skip lines beginning
                                                          ! with ! or ?
      input1=''
      quoted=.false.
      do i=1,len
        iqq=scan(quote,input(i:i))
        if(iqq.ne.0)then
          if(.not.quoted)then
             quoted=.true.
             iq=iqq
          else
             if(iq.eq.iqq)quoted=.false.
          endif
        endif
        if(input(i:i).eq.excl.and..not.quoted)exit  ! skip comments
        input1(i:i)=input(i:i)
      enddo

      write(ichan,'(a)')input1(1:len_trim(input1))

      end subroutine writeline

c=====================================================================
      subroutine rmblan(string,len1,len)
      character*1 string,blank
      dimension string(len1)
      data blank/' '/
c     counts the trailing blanks from a string
c     len1 is the original length of the string, len is the
c     length after the trailing blanks have been removed
      len=0
      do 100 k=len1,1,-1
        if(string(k).eq.blank) go to 100
        len=k
        go to 200
 100  continue
 200  return
      end
c=====================================================================
      subroutine readopt(inpf,nopt,opnames,ioptyp,iopval,
     $                   ropval,chopval,ifound)
c   this subroutine scans inpf
c   if it finds one of the opnames, it consults ioptyp corresponding
c  to this particular option name. If ioptyp=0, than this is essentially
c   a logical option, and iopval corresponding to this option is set
c   to 1
c   if ioptyp=k, (k<4) then this is an integer option, and the program
c   reads  the next k integers and returns them in iopval(k,ii)
c   note that if fewer than k integers are present than fewer are read
c   if ioptyp=10+k, then k reaL VALUES ARE READ
c   if ioptyp=21 than a string is read (terminated by the first blank)
c   the results of the k-th option are returned in iopval(3,k),
c   ropval(3,k), or chopval(k), depending on whether the option was
c   integer, real or character type. Integer and real options can have
c   up to 3 values
c  What was missing is the kind of option which may or may not have a
c  value. We can assume that a value is not necessary for an integer
c  option, the value being zero be default
c   ifound is 0 if the option in question was not found, and 1 otherwise
c
c   MM 04/18/2007 gave the routine a clenup. Deleted unused variables,
c   set line width to 300 to allow for long filenames (Windows).
c
      implicit real*8 (a-h,o-z)
      character*4 opnames(nopt),b
      character*256 chopval
      logical first
      parameter (linew=300)
      character(len=linew) :: a
      character*256 c
      dimension ioptyp(nopt),iopval(3,nopt),ropval(3,nopt),chopval(nopt)
      dimension ifound(nopt)
c
      call izeroit(ifound,nopt)
c   if a wrong keyword is encountered on the first card at the first
c   position, do not backspace the input
      first=.false.
 30   continue
      a=''
      read(inpf,'(a)',end=1200)a
c  first card - skip the first word until the first  space or comma or (
c  initial skip. do not get back here again
      i=1
      call rdword(a,i,linew,c,256,j,ival)
 200  call rdword(a,i,linew,b,4,j,ival)
      call lowercas(b,4)
      if (j.ne.0) then
c  a word was found
c  analyze  the option name
        do 700 k=1,nopt
          if(opnames(k).eq.b) go to 800
 700    continue
c  option name does not ring a bell
        go to 1200
 800    ifound(k)=1
        if (ioptyp(k).lt.10) then
         iopval(1,k)=0
         iopval(2,k)=0
         iopval(3,k)=0
        else if (ioptyp(k).lt.20) then
          ropval(1,k)=0.0d0
          ropval(2,k)=0.0d0
          ropval(3,k)=0.0d0
        else if (ioptyp(k).eq.21) then
          call blankit(chopval(k),256)
        end if
        if (ioptyp(k).eq.0) then
          iopval(1,k)=1
        end if
        if(ival.eq.0) then
          l=0
        else
          call rdword(a,i,linew,c,256,l,ival)
        end if
        if (l.eq.0.and.ioptyp(k).gt.1) then
          call message('readopt','option value is missing,zero assumed',
     $                 k,l)
        end if
        if(l.gt.0) then
          if (ioptyp(k).eq.1) then
            read(c,*,end=900) iopval(1,k)
          else if (ioptyp(k).eq.2) then
            read(c,*,end=900) iopval(1,k),iopval(2,k)
          else if (ioptyp(k).eq.3) then
            read(c,*,end=900) iopval(1,k),iopval(2,k),iopval(3,k)
          else if (ioptyp(k).eq.11) then
            read(c,*,end=900) ropval(1,k)
          else if (ioptyp(k).eq.12) then
            read(c,*,end=900) ropval(1,k),ropval(2,k)
          else if (ioptyp(k).eq.13) then
            read(c,*,end=900) ropval(1,k),ropval(2,k),ropval(3,k)
          else if (ioptyp(k).eq.21) then
            chopval(k)=c
          else if (ioptyp(k).ne.0) then
            call nerror(2,'readopt','Wrong ioptyp', ioptyp(k),0)
          end if
        end if
        first=.false.
c  look at the next option on this card
 900    continue
        go to 200
      else
        first=.true.
        a=''
        read(inpf,'(a)',end=1200) a
        i=1
        go to 200
      end if
 1200 if(first) then
         backspace inpf
       else
         call nerror(2,'readopt','wrong keyword= '//b,0,0)
       end if
       end
c=====================================================================
      subroutine readop1(inpf,nopt,opnames,ioptyp,iopval,
     $                   ropval,chopval,ifound)
c This is similar to readopt but it does not need an extra name parameter
c  E.g. instead of BASIS STAN=6-31G* it reads BASIS=6-31G*
c  or MEMO=4.5
c  It also allows parentheses in the file name arguments
c   this subroutine scans inpf
c   if it finds one of the opnames, it consults ioptyp corresponding
c  to this particular option name. If ioptyp=0, than this is essentially
c   a logical option, and iopval corresponding to this option is set
c   to 1
c   if ioptyp=k, (k<4) then this is an integer option, and the program
c   reads  the next k integers and returns them in iopval(k,ii)
c   note that if fewer than k integers are present than fewer are read
c   if ioptyp=10+k, then k reaL VALUES ARE READ
c   if ioptyp=21 than a string is read (terminated by the first blank)
c   the results of the k-th option are returned in iopval(3,k),
c   ropval(3,k), or chopval(k), depending on whether the option was
c   integer, real or character type. Integer and real options can have
c   up to 3 values
c  What was missing is the kind of option which may or may not have a
c  value. We can assume that a value is not necessary for an integer
c  option, the value being zero be default
c   ifound is 0 if the option in question was not found, and 1 otherwise
c
c   MM 04/18/2007 gave the routine a clenup. Deleted unused variables,
c   set line width to 300 to allow for long filenames (Windows).
c
      implicit real*8 (a-h,o-z)
      character*4 opnames(nopt),b
      character*256 chopval
      logical first
      parameter (linew=300)
      character(len=linew) :: a
      character*256 c
      dimension ioptyp(nopt),iopval(3,nopt),ropval(3,nopt),chopval(nopt)
      dimension ifound(nopt)
c
      call izeroit(ifound,nopt)
c   if a wrong keyword is encountered on the first card at the first
c   position, do not backspace the input
      first=.false.
c I think "first" shows whether the parser has found a card already
c which makes sense. If not, it just reads in.
c I need to remove the error message if the very first option
c does not have an equal sign
      a=''
      read(inpf,'(a)',end=1200) a
c  in readop1, do not skip the first word until the first space or comma or (
      i=1
 200  call rdword(a,i,linew,b,4,j,ival)
      call lowercas(b,4)
      IF (j.ne.0) THEN
c  a word was found
c  analyze  the option name
      do 700 k=1,nopt
        if(opnames(k).eq.b) go to 800
 700  continue
c  option name does not ring a bell
      go to 1200
 800  ifound(k)=1
      if (ioptyp(k).lt.10) then
       iopval(1,k)=0
       iopval(2,k)=0
       iopval(3,k)=0
      else if (ioptyp(k).lt.20) then
        ropval(1,k)=0.0d0
        ropval(2,k)=0.0d0
        ropval(3,k)=0.0d0
      else if (ioptyp(k).eq.21) then
        call blankit(chopval(k),256)
      end if
      if (ioptyp(k).eq.0) then
        iopval(1,k)=1
      end if
      if(ival.eq.0) then
        l=0
      else
        call rdword1(a,i,linew,c,256,l,ival)
      end if
cc      if (l.eq.0.and.ioptyp(k).gt.1.and.izeropt.gt.1) then
cc      call message('readop1','option value is missing,zero assumed',k,l)
cc      end if
      if(l.gt.0) then
        if (ioptyp(k).eq.1) then
          read(c,*,end=900) iopval(1,k)
        else if (ioptyp(k).eq.2) then
          read(c,*,end=900) iopval(1,k),iopval(2,k)
        else if (ioptyp(k).eq.3) then
          read(c,*,end=900) iopval(1,k),iopval(2,k),iopval(3,k)
        else if (ioptyp(k).eq.11) then
          read(c,*,end=900) ropval(1,k)
        else if (ioptyp(k).eq.12) then
          read(c,*,end=900) ropval(1,k),ropval(2,k)
        else if (ioptyp(k).eq.13) then
          read(c,*,end=900) ropval(1,k),ropval(2,k),ropval(3,k)
        else if (ioptyp(k).eq.21) then
            chopval(k)=c(1:len_trim(c))
        else if (ioptyp(k).ne.0) then
          call nerror(2,'readop1','Wrong ioptyp', ioptyp(k),0)
        end if
      end if
      first=.false.
c  look at the next option on this card
 900   continue
       go to 200
      ELSE
        first=.true.
        a=''
        read(inpf,'(a)',end=1200) a
        i=1
        go to 200
      ENDIF
 1200  if(first) then
         backspace inpf
       else
        call nerror(3,'readop1','wrong keyword= '//b,0,0)
       end if
       end
c=====================================================================

       subroutine rdword(a,i,lg,b,khar,j1,ival)
c   this routine parses the string a from position i. It skips
c  the spaces, equal signs, commas, parentheses until the first
c  non-space, non =, non comma or parenthesis character. From this
c  point, it copies the content of a into b, until
c  (a) it encounters a space, = or paranthesis, or
c  (b) it fills up khar characters in b
c   j1 is the number of returned characters in b
c   if the option name is terminated by an equal sign, ival is set to 1
c
c   MM 04/18/2007  now the routine allows for quoting: if a quotation
c      character (",',`) is encountered, the routine will continue
c      copying until a matching quote is encountered (the quotation
c      characters themselves are skipped) then it will resume normal
c      operation, e.g., it will look for a space, comma, equal sign or
c      parenthesis. This is to allow spaces (and other characters as
c      well) into file names (mainly for Windows)
c
      character(len=lg) :: a
      character(len=khar) :: b
      character(len=3), parameter :: quote='"''`' ! quotation characters
      character(len=5), parameter :: sepa=' =(),' ! separators
      character(len=4), parameter :: sepb=' =()'  ! separators for b string
      logical quoted

      quoted=.false.
      ival=0
      init=i
      i=i-1
c   j is the subscript to b
      b=''
      j=1
c  initial skipping
 100  i=i+1
      if(i.le.lg)then
        if(scan(sepa,a(i:i)).ne.0) goto 100
c  copying
 200    continue
                        ! check for quotes
        iqq=scan(quote,a(i:i))
        if(iqq.ne.0)then
          if(.not.quoted)then
             quoted=.true.
             iq=iqq
             i=i+1
             goto 200
          else if(iq.eq.iqq) then
             quoted=.false.
             i=i+1
             goto 200
          endif
        endif

c  commas are permitted in b!

        if(i.le.lg.and.j.le.khar)then
  
          if(quoted)then
            b(j:j)=a(i:i)
            j=j+1
            i=i+1
            if(i.le.lg)go to 200
          else if(scan(sepb,a(i:i)).eq.0)then
            b(j:j)=a(i:i)
            j=j+1
            i=i+1
            if(i.le.lg)go to 200
          endif

        else if(i.gt.lg) then

          if(quoted) 
     $       call nerror(1,'Rdword','Unmatched quote in'//a(init:i),0,0)
          j1=j-1
          return

        else if(i.le.lg.and.j.gt.khar)then

          if(quoted) 
     $       call nerror(1,'Rdword','Unmatched quote in'//a(init:i),0,0)
          if(scan(sepb,a(i:i)).eq.0)then
            i=i+1
            go to 200
          endif

        end if

      endif

      if(quoted) 
     $   call nerror(1,'Rdword','Unmatched quote in'//a(init:i),0,0)
      j1=j-1
      if(i.le.lg)then
        if(a(i:i).eq.'=') ival=1
      endif
      end

c=====================================================================
      subroutine rdword1(a,i,lg,b,khar,j1,ival)
c   this routine parses the string a from position i. Itv skips
c  the spaces and equal signs but retains commas and parentheses
c  until the first non-space, non-equal sign, non-comma, non-
c  parenthesis character From this
c  point, it copies the content of a into b, until
c  (a) it encounters a space or equal sign, or
c  (b) it fills up khar characters in b
c   j1 is the number of returned characters in b
c   if the option name is terminated by an equal sign, ival is set to 1
c  this differs from rdword in that the string can contain commas and/or
c  parenthese
c
c   MM 04/18/2007  now the routine allows for quoting: if a quotation
c      character (",',`) is encountered, the routine will continue
c      copying until a matching quote is encountered (the quotation
c      characters themselves are skipped) then it will resume normal
c      operation, e.g., it will look for a space, comma, equal sign or
c      parenthesis. This is to allow spaces (and other characters as
c      well) into file names (mainly for Windows)
c
      character(len=lg) :: a
      character(len=khar) :: b
      character(len=3), parameter :: quote='"''`' ! quotation characters
      character(len=5), parameter :: sepa=' =(),' ! separators
      character(len=2), parameter :: sepb=' ='    ! separators for b string
      logical quoted

      quoted=.false.
      ival=0
      init=i
      i=i-1
c   j is the subscript to b
      b=''
      j=1
c  initial skipping
 100  i=i+1
      if(i.le.lg)then
        if(scan(sepa,a(i:i)).ne.0) goto 100
c  copying
 200    continue
                        ! check for quotes
        iqq=scan(quote,a(i:i))
        if(iqq.ne.0)then
          if(.not.quoted)then
             quoted=.true.
             iq=iqq
             i=i+1
             goto 200
          else if(iq.eq.iqq) then
             quoted=.false.
             i=i+1
             goto 200
          endif
        endif

c  commas and parentheses are permitted in b!

        if(i.le.lg.and.j.le.khar)then

          if(quoted)then
            b(j:j)=a(i:i)
            j=j+1
            i=i+1
            if(i.le.lg)go to 200
          else if(scan(sepb,a(i:i)).eq.0)then
            b(j:j)=a(i:i)
            j=j+1
            i=i+1
            if(i.le.lg)go to 200
          endif

        else if(i.gt.lg) then

          if(quoted) 
     $      call nerror(1,'Rdword1','Unmatched quote in'//a(init:i),0,0)
          j1=j-1
        if(a(i:i).eq.'=') ival=1
          return

        else if(i.le.lg.and.j.gt.khar)then
          if(quoted) 
     $      call nerror(1,'Rdword1','Unmatched quote in'//a(init:i),0,0)
          if(scan(sepb,a(i:i)).eq.0)then
            i=i+1
            go to 200
          endif

        endif

      endif

      if(quoted) 
     $  call nerror(1,'Rdword1','Unmatched quote in'//a(init:i),0,0)
      j1=j-1
      if(a(i:i).eq.'=') ival=1
      end

c=====================================================================

      subroutine blankit(b,khar)
      character*1 b(khar),space
      data space/' '/
      do 100 i=1,khar
        b(i)=space
 100  continue
      end
c=====================================================================
      subroutine lowercas(string,n)
c  converts string to lower case
      character*1 string(n)
      character*1 lower(26),capit(26)
      data capit/'A','B','C','D','E','F','G','H','I','J','K','L','M',
     1 'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1 'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      do 100 i=1,n
       do 100 j=1,26
         if(string(i).eq.capit(j)) then
           string(i)=lower(j)
         end if
 100  continue
      end
c=====================================================================
      subroutine lowerca2(string,n)
c  converts the fist n characters of a string to lower case
c  n must be smaller or equal than the length of the string
      character*(*) string
      character*1 lower(26),capit(26)
      data capit/'A','B','C','D','E','F','G','H','I','J','K','L','M',
     1 'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1 'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      do 100 i=1,n
       do 100 j=1,26
         if(string(i:i).eq.capit(j)) then
           string(i:i)=lower(j)
         end if
 100  continue
      end
c====================================================================
      subroutine uppercase(string,n)
c  converts the first n characters of a string to upper case
c  n must be smaller or equal than the length of the string
      character*(*) string
      character*1 lower(26),capit(26)
      data capit/'A','B','C','D','E','F','G','H','I','J','K','L','M',
     1 'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1 'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      do 100 i=1,n
       do 100 j=1,26
         if(string(i:i).eq.lower(j)) then
           string(i:i)=capit(j)
         end if
 100  continue
      end
c====================================================================
      subroutine findrootd(string,length)

      use sysdef

      character*(*) string
      character*1 blank
      data blank/' '/
      do i=length,1,-1
        if(string(i:i).eq.DIR_SEP) exit
        string(i:i)=blank
      end do
      end
c====================================================================
      subroutine leadblan2(string,len1,len)
c     removes the leading blanks from a string
c     len1 is the original length of the string, len is the
c     length after the leading blanks have been removed
      character*(*) string
      character*1 blank
      data blank/' '/
c  count the number of leading blanks
      k1=0
      do k=1,len1
        if(string(k:k).ne.blank.or.k1.gt.0) k1=k1+1
        if(k.gt.k1) string(k1:k1)=string(k:k)
      end do
      len=k1
      end
c
c====================================================================
c
      subroutine rmblan2(string,len1,len)
      character*(*) string
      character*1 blank
      data blank/' '/
c     counts the trailing blanks from a string
c     len1 is the original length of the string, len is the
c     length after the trailing blanks have been removed
      len=0
      do 100 k=len1,1,-1
        k1=k
        if(string(k1:k1).eq.blank) go to 100
        len=k
        go to 200
 100  continue
 200  return
      end
c====================================================================
      subroutine basename(filnm,LenIn,Basnm,LenOut)
c  This routine finds the basename (i.e. the name without the extension)
c  of a  filename. E.g. the basename of butane.inp is butane. The basename of
c  butane is butane
c  Arguments:
c  INTENT(INOUT):
c  string*(*) = a string at most 256 characters long. NOTE: the program
c  will remove any leading blanks in filnm, and also remove the last
c  periosd and characters following it
c  INTENT(IN):
c  LenIn (no connection to Vladimir Ilyich): the length of the input string
c  INTENT(OUT):
c  length     = the length of the resulting basename string. The basename id
c  filnm(1:length)
c
c  MM 04/18/2007  Now the routine allows for directory separators in the
c   input name, in other words, the search for a file extension is
c   limited to the last part of the file name. This is to allow input
c   names containing a path, and in particular input names with a
c   relative path (which has dot characters on it)
c
      use sysdef

      character*(*) filnm,basnm

      basnm=''
      len1=len_trim(filnm)
      if(len1.gt.256)len1=256
      basnm=adjustl(filnm(1:len1))
      len1=len_trim(basnm)

               ! scan for last directory separator

      isep=scan(basnm(1:len1),DIR_SEP,back=.true.)
      if(isep.ge.len1) call nerror(1,'Basename',
     $      'Empty basename: '//basnm(1:len1),0,0)

               ! scan for last dot in the filename part of basnm

      idot=scan(basnm(isep+1:len1),'.',back=.true.)
      if(idot.eq.0)then
        lenout=len1
      else
        lenout=isep+idot-1
      endif

      if(lenout.lt.len1)basnm(lenout+1:len1)=''

      end
