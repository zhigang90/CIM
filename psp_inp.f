      subroutine pspupdateab(nalpha,nbeta)
      implicit real*8 (a-h,o-z)
c...
c...  MM August 2004
c...
c...  this routine updates the number of alpha and beta
c...  electrons in the <control> file
c...
      character*256 jobname
      logical busy
      real*8 rdum
      character*20 cdum
c...
c...
c...  look for a free I/O channel
c...
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
        endif
      enddo
      call nerror(1,'Pspupdateab','Cannot find a free channel',0,0)
 10   continue
c...
c...  open <control> file
c...
      call getchval('jobname',jobname)
      open(unit=ich,file=jobname(1:len_trim(jobname))//'.control',
     $     status='old',form='formatted')
c...
c...  replace alpha an beta electrons
c...
      call wrcntrl(ich,7,'$nalpha',1,nalpha,rdum,cdum)
      call wrcntrl(ich,7,'$nbeta ',1,nbeta,rdum,cdum)
c...
c...  done
c...
      close(unit=ich,status='keep')
      end
c======================================================================
      subroutine pspupdatech(na,xnuc)

      use memory

      implicit real*8 (a-h,o-z)
c...
c...  MM August 2004
c...
c...  this routine updates the <coord> file with the
c...  psp corrected atomic charges
c...
      dimension xnuc(5,na)
c...
c...  we need an area for the atomic symbols. hopefully
c...  10 000 atoms will be enough for a while
c...
      character*256 jobname
      parameter (maxat=10000)
      character*8 atsymb(maxat)
c     common /big/bl(100)
      logical busy
c...
c...
      if(na.gt.maxat)then
        call nerror(1,'Pspupdatech','Too many atoms!',na,maxat)
      endif
c...
c...  we need some memory for reading stuff
c...
      call getival('NMol',nmol)
      call mmark
      call getmem(3*na,ixyz)
      call getmem(na,icharge)
      call getmem(na,ixm)
      call getmem(na+1,imol)
c...
c...  read <coord> file
c...
      call getchval('jobname',jobname)
c
c
c   look for a free I/O channel
c
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(2,'Pspupdatech','Cannot find a free channel',0,0)
 10   continue
      open(unit=ich,file=jobname(1:len_trim(jobname))//'.coord',
     $     status='old',form='formatted')
      call RdCoordF(ich,na,atsymb,bl(ixyz),nmol,bl(imol),bl(icharge),
     $              bl(ixm))
      close(unit=ich,status='keep')
c...
c...  update charges
c...
      do i=1,na
        bl(icharge+i-1)=xnuc(1,i)
      enddo
c...
c...  write back <coord> file
c...
      call WrCoord(na,atsymb,bl(ixyz),nmol,bl(imol),bl(icharge),bl(ixm))
c...
c...  done
c...
      call retmark
      end
c======================================================================
      subroutine psptmpread(icpsp,icradp,npsp,nradp,nelpsp,iecpdat,
     $                   nrrad,aradp,cradp,maxl)
      implicit real*8 (a-h,o-z)
c...
c...  MM August 2004
c...
c...  this routine reads ithe pseudopotential data from temporary
c...  storage
c...
c...  icpsp      input channel for general psp data
c...  icradp     input channel for radial psp data
c...  npsp       number of pseudopotential
c...  nradp      number of radial terms
c...  nelpsp     number of electron represented by psp
c...  iecpdat    storage for pseudopotential data
c...  nrrad      storage for powers of r in radial terms
c...  aradp      storage for exponent of radial terms
c...  cradp      storage for coefficient of radial terms
c...  maxl       maximum angular momentum we can treat
c...
      dimension iecpdat(maxl+6,npsp),nrrad(nradp),aradp(nradp),
     $          cradp(nradp)
      rewind icpsp
      rewind icradp
c...
c...  skip file headers
c...
      read(icpsp,*)
      read(icradp,*)
      read(icpsp,*)
      read(icradp,*)
      read(icpsp,*)
      read(icradp,*)
c...
c...  read data
c...
      do i=1,npsp
        read(icpsp,*)(iecpdat(j,i),j=1,maxl+6)
      enddo
      do i=1,nradp
        read(icradp,*)nrrad(i),aradp(i),cradp(i)
      enddo
c...
c...  consistency check
c...
      nelr=0
      do i=1,npsp
        nelr=nelr+iecpdat(2,i)
      enddo
      if(nelr.ne.nelpsp)
     $  call nerror(1,'Psptmpread','Inconsistency detected',nelr,nelpsp)
c...
      end
c======================================================================
      subroutine rpspinp(icin,icpsp,icradp,atsy,na,npsp,nradp,nelpsp,
     $                   maxl,lprint,iout)
      implicit real*8 (a-h,o-z)
c...
c...  MM August 2004
c...
c...  this routine reads pseudopotentials (psp) from the user input
c...  and store the information on the two files icsps and icradp.
c...
c...  it is an incremental read, i.e. it can be called several times
c...
c...  icin       input channell
c...  icpsp      output channel for general psp data
c...  icradp     output channel for radial psp data
c...  atsy       atomic symbol of psp being read
c...  na         number of atoms
c...  npsp       number of pseudopotential
c...  nradp      number of radial terms
c...  nelpsp     number of electron represented by psp
c...  maxl       maximum angular momentum we can treat
c...  lprint     print flag
c...  iout       output channell
c...
      character*8 atsy,atsyex
      character*2 atsym
      parameter (lwrd=120)
      character*128 line
      dimension inica(lwrd/2+1),ifica(lwrd/2+1)
c...
c...  maximum number of angular term we can read for a single psp
c...
      parameter (maxradp=100)
      parameter (maxldim=10)
      dimension iecpdat(maxldim+6),iradp(maxradp),radp(2,maxradp)
c...
      character*1 pot(0:maxldim)
      data pot/'s','p','d','f','g','h','i',4*'x'/
c...
c...  consistency check
c...
      if(maxl.gt.maxldim)
     $  call nerror(1,'Rpspinp','Inconsistency in dimensioning',
     $              maxl,maxldim)
c
      iout = igetival('iout')
c...
c...  the first line must be of type
c...
c...  'psp ' izcore lmax+1
c...
c...  or
c...
c...  'ecp ' izcore lmax+1
c...
c...  that is, at least three fields
c...
      call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
      call lowerca2(line,lwrd)
      if(noca.ge.3.and.iend.eq.0)then  ! read  the rest of this line
        read(line(inica(2):ifica(2)),*)izcore
        if(izcore.lt.0.or.mod(izcore,2).ne.0)then
          write(iout,*)'Rpspinp: incorrect core charge'
          write(iout,*)line
          call nerror(1,'rpspinp','incorrect core charge',0,0)
        endif
        read(line(inica(3):ifica(3)),*)lmax
        if(lmax.le.0.or.lmax.gt.maxl+1)then
          write(iout,*)'Rpspinp: incorrect angular momentum'
          write(iout,*)line
          call nerror(2,'rpspinp','incorrect angular momentum',0,0)
        endif
c...
c...  optional printing
c...
        if(lprint.gt.0)then
          write(iout,100)izcore,lmax
100     format(1x,'Pseudopotential: izcore=',i3,' lmax+1=',i3)
        endif
        call izeroit(iecpdat,maxl+6)
        iecpdat(1)=lmax-1
        iecpdat(2)=izcore
        call zeroit(radp,2*maxradp)
        call izeroit(iradp,maxradp)
        irad=0
c...
c...   now we expect a line of type
c...
c...   ngpot
c...
c...   that is, at least one field
c...
        call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
        if(noca.ge.1.and.iend.eq.0)then
          read(line(inica(1):ifica(1)),*)ngpot
          if(ngpot.le.0)then
            write(iout,*)'Rpspinp: incorrect number of radial terms'
            write(iout,*)line
          call nerror(1,'rpspinp','incorrect no, of radial terms',0,0)
          endif
        else
          write(iout,*)'Rpspinp: error while reading',
     $              ' number of radial terms'
          write(iout,*)line
          call nerror(1,'rpspinp','incorrect no, of radial terms',0,0)
        endif
c...
c...  optional printing
c...
        if(lprint.gt.0)then
          write(iout,101)ngpot,pot(lmax)
        endif
101     format(1x,i2,6x,'-----',1x,a1,' potential',5x,'-----')
c...
c...   now we expect ngpot lines of the type
c...
c...   clp nlp zlp
c...
c...   that is, at least three fields
c...
        ngpotl=ngpot
        iecpdat(4)=1
        do nlin=1,ngpotl
          call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
          if(noca.ge.3.and.iend.eq.0)then
            irad=irad+1
            read(line(inica(1):ifica(3)),*)radp(1,irad),iradp(irad),
     $                                    radp(2,irad)
          else
            write(iout,*)'Rpspinp: error while reading',
     $                ' radial term'
            write(iout,*)line
            call nerror(1,'rpspinp','incorrect no, of radial terms',0,0)
          endif
        enddo
c...
c...  optional printing
c...
        if(lprint.gt.0)then
          irpr=0
          do nlin=1,ngpotl
            irpr=irpr+1
            write(iout,102)radp(1,irpr),iradp(irpr),radp(2,irpr)
          enddo
        endif
102     format(1x,'coeff=',f14.8,' n=',i2,' exp=',f14.8)
c...
c... loop over angular momentum
c...
        do l=0,lmax-1
          iecpdat(5+l)=irad+1
c...
c...   now we expect a line of type
c...
c...   ngpot
c...
          call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
          if(noca.ge.1.and.iend.eq.0)then
            read(line(inica(1):ifica(1)),*)ngpot
            if(ngpot.le.0)then
              write(iout,*)'Rpspinp: incorrect number of radial terms'
              write(iout,*)line
            call nerror(1,'rpspinp','incorrect no, of radial terms',0,0)
            endif
          else
            write(iout,*)'Rpspinp: error while reading',
     $                ' number of radial terms'
            write(iout,*)line
            call nerror(1,'rpspinp','incorrect no, of radial terms',0,0)
          endif
c...
c...  optional printing
c...
          if(lprint.gt.0)then
            write(iout,103)ngpot,pot(l),pot(lmax)
          endif
103       format(1x,i2,6x,'-----',1x,a1,'-',a1,' potential',3x,'-----')
c...
c...   now we expect ngpot lines of the type
c...
c...   clp nlp zlp
c...
          do nlin=1,ngpot
            call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
            if(noca.ge.3.and.iend.eq.0)then
              irad=irad+1
              read(line(inica(1):ifica(3)),*)radp(1,irad),
     $               iradp(irad),radp(2,irad)
            else
              write(iout,*)'Rpspinp: error while reading',
     $                  ' radial term'
              write(iout,*)line
            call nerror(1,'rpspinp','incorrect no, of radial terms',0,0)
            endif
          enddo
c...
c...  optional printing
c...
          if(lprint.gt.0)then
            do nlin=1,ngpot
              irpr=irpr+1
              write(iout,102)radp(1,irpr),iradp(irpr),radp(2,irpr)
            enddo
          endif
        enddo
        iecpdat(maxl+6)=irad
c...
c...   the reading has been successfull, now we are ready
c...   to store the data on the files
c...
        do n=1,na
          call getatsymex(n,na,atsyex)
          if(atsy.eq.atsyex)then
            write(icpsp,'(i7,50(1x,i7))')(iecpdat(i),i=1,2),n,
     $      (nradp+iecpdat(l),l=4,lmax+4),
     $      (iecpdat(l),l=lmax+5,maxl+5),nradp+iecpdat(maxl+6)
            npsp=npsp+1
            nelpsp=nelpsp+izcore
            do i=1,irad
              write(icradp,'(i7,1x,e25.17,1x,e25.17)')
     $                   iradp(i),radp(2,i),radp(1,i)
            enddo
            nradp=nradp+irad
          endif
        enddo
      else
        write(iout,*)'Rpspinp: incorrect pseudopotential specification'
        write(iout,*)line
        call nerror(1,'rpspinp','incorrect psp specification',0,0)
      endif
      end
c======================================================================
      subroutine rpsplib(icin,icpsp,icradp,atsy,na,npsp,nradp,nelpsp,
     $                   maxl,lprint,iout)
      implicit real*8 (a-h,o-z)
c...
c...  MM August 2004
c...
c...  this routine reads pseudopotentials (psp) from library files
c...  and store the information on the two files icsps and icradp.
c...
c...  if atsy is not blank, only the pseudopotentials for atom
c...  atsy will be read, otherwise all atoms that can be matched
c...  will be read
c...
c...  it is an incremental read, i.e. it can be called several times
c...
c...  icin       input channell
c...  icpsp      output channel for general psp data
c...  icradp     output channel for radial psp data
c...  atsy       atomic symbol if just one atom type is to be
c...             considered, otherwise blank
c...  na         number of atoms
c...  npsp       number of pseudopotential
c...  nradp      number of radial terms
c...  nelpsp     number of electron represented by psp
c...  maxl       maximum angular momentum we can treat
c...  lprint     print flag
c...  iout       output channell
c...
      character*8 atsy,atsyex,reasy
      character*2 atsym
      parameter (lwrd=120)
      character*128 line
      character*24 pspline
      dimension inica(lwrd/2+1),ifica(lwrd/2+1)
      logical special,match
c...
c...  maximum number of angular term we can read for a single psp
c...
      parameter (maxradp=100)
      parameter (maxldim=10)
      dimension iecpdat(maxldim+6),iradp(maxradp),radp(2,maxradp)
c...
      character*1 pot(0:maxldim)
      data pot/'s','p','d','f','g','h','i',4*'x'/
c...
c...  consistency check
c...
      if(maxl.gt.maxldim)
     $  call nerror(1,'Rpsplib','Inconsistency in dimensioning',
     $              maxl,maxldim)
c
      iout = igetival('iout')
c...
      if(atsy.ne.'        ')then
        special=.true.
        atsym=atsy(1:2)
        jo=ichar(atsym(2:2))
        if(jo.lt.97.or.jo.gt.122)atsym(2:2)=' '
      else
        special=.false.
      endif
c...
c...  the effective core potential section starts with the line
c...
c...  'Effective Core Potential'
c...
  50  continue
      read(icin,'(a)',end=999)pspline
      call lowerca2(pspline,24)
      if(pspline.ne.'effective core potential')goto 50
c...
c...  look for a line of type:
c...
c...  symbol ptype izcore lmax+1
c...
c...  that is, at lest four fields
c...
 100  continue
      call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
      if(iend.ne.0)goto 999
      call lowerca2(line,lwrd)
      if(noca.ge.4)then
        call blankit(reasy,2)
        reasy=line(inica(1):min(ifica(1),inica(1)+1)) ! get symbol
        jo=ichar(reasy(2:2))
        if(jo.lt.97.or.jo.gt.122)reasy(2:2)=' '
        match=.false.
        if(special)then
          match=reasy(1:2).eq.atsym
        else
          do i=1,na
            call getatsymex(i,na,atsyex)
            if(reasy.eq.atsyex)match=.true.
          enddo
        endif
        if(match)then  ! read  the rest of this line
          read(line(inica(3):ifica(3)),*,err=100)izcore
          if(izcore.lt.0.or.mod(izcore,2).ne.0)then
            write(iout,*)'Rpsplib: incorrect core charge'
            write(iout,*)line
            call nerror(1,'rpsplib','incorrect core charge',0,0)
          endif
          read(line(inica(4):ifica(4)),*,err=100)lmax
          if(lmax.le.0.or.lmax.gt.maxl+1)then
            write(iout,*)'Rpsplib: incorrect angular momentum'
            write(iout,*)line
            call nerror(1,'rpsplib','incorrect angular momentum',0,0)
          endif
          call izeroit(iecpdat,maxl+6)
          iecpdat(1)=lmax-1
          iecpdat(2)=izcore
          call zeroit(radp,2*maxradp)
          call izeroit(iradp,maxradp)
          irad=0
c...
c...      optional printing
c...
          if(lprint.gt.0)then
            if(special)then
              write(iout,199)atsy
            else
              write(iout,199)reasy
            endif
            write(iout,200)izcore,lmax
          endif
199       format(' For atom: ',a4)
200       format(1x,'Pseudopotential: izcore=',i3,' lmax+1=',i3)
c...
c...   now we expect a line of type
c...
c...   ngpot
c...
c...   that is, at least one field
c...
          call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
          if(noca.ge.1.and.iend.eq.0)then
            read(line(inica(1):ifica(1)),*)ngpot
            if(ngpot.le.0)then
              write(iout,*)'Rpsplib: incorrect number of radial terms'
              write(iout,*)line
              call nerror(1,'rpsplib','wrong no. of radial terms',0,0)
            endif
          else
            write(iout,*)'Rpsplib: error while reading',
     $                ' number of radial terms'
            write(iout,*)line
            call nerror(1,'rpsplib','wrong no. of radial terms',0,0)
          endif
c...
c...  optional printing
c...
          if(lprint.gt.0)then
            write(iout,201)ngpot,pot(lmax)
          endif
201       format(1x,i2,6x,'-----',1x,a1,' potential',5x,'-----')
c...
c...   now we expect ngpot lines of the type
c...
c...   clp nlp zlp
c...
c...   that is, at least three fields
c...
          ngpotl=ngpot
          iecpdat(4)=1
          do nlin=1,ngpotl
            call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
            if(noca.ge.3.and.iend.eq.0)then
              irad=irad+1
              read(line(inica(1):ifica(3)),*)radp(1,irad),iradp(irad),
     $                                      radp(2,irad)
            else
              write(iout,*)'Rpsplib: error while reading',
     $                  ' radial term'
              write(iout,*)line
              call nerror(1,'rpsplib','wrong no. of radial terms',0,0)
            endif
          enddo
c...
c...  optional printing
c...
          if(lprint.gt.0)then
            irpr=0
            do nlin=1,ngpotl
              irpr=irpr+1
              write(iout,202)radp(1,irpr),iradp(irpr),radp(2,irpr)
            enddo
          endif
202       format(1x,'coeff=',f14.8,' n=',i2,' exp=',f14.8)
c...
c... loop over angular momentum
c...
          do l=0,lmax-1
            iecpdat(5+l)=irad+1
c...
c...   now we expect a line of type
c...
c...   ngpot
c...
            call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
            if(noca.ge.1.and.iend.eq.0)then
              read(line(inica(1):ifica(1)),*)ngpot
              if(ngpot.le.0)then
                write(iout,*)'Rpsplib: incorrect number of radial terms'
                write(iout,*)line
                call nerror(1,'rpsplib','wrong no. of radial terms',0,0)
              endif
            else
              write(iout,*)'Rpsplib: error while reading',
     $                  ' number of radial terms'
              write(iout,*)line
              call nerror(1,'rpsplib','wrong no. of radial terms',0,0)
            endif
c...
c...  optional printing
c...
             if(lprint.gt.0)then
               write(iout,203)ngpot,pot(l),pot(lmax)
             endif
203        format(1x,i2,6x,'-----',1x,a1,'-',a1,' potential',3x,'-----')
c...
c...   now we expect ngpot lines of the type
c...
c...   clp nlp zlp
c...
            do nlin=1,ngpot
              call realin(icin,line,noca,inica,ifica,lwrd,iend,0)
              if(noca.ge.3.and.iend.eq.0)then
                irad=irad+1
                read(line(inica(1):ifica(3)),*)radp(1,irad),
     $                 iradp(irad),radp(2,irad)
              else
                write(iout,*)'Rpsplib: error while reading',
     $                    ' radial term'
                write(iout,*)line
                call nerror(1,'rpsplib','wrong no. of radial terms',0,0)
              endif
            enddo
c...
c...  optional printing
c...
            if(lprint.gt.0)then
              do nlin=1,ngpot
                irpr=irpr+1
                write(iout,202)radp(1,irpr),iradp(irpr),radp(2,irpr)
              enddo
            endif
          enddo
          iecpdat(maxl+6)=irad
c...
c...   the reading has been successfull, now we are ready
c...   to store the data on the files
c...
          if(special)then
            do n=1,na
              call getatsymex(n,na,atsyex)
              if(atsy.eq.atsyex)then
                write(icpsp,'(i7,50(1x,i7))')(iecpdat(i),i=1,2),n,
     $          (nradp+iecpdat(l),l=4,lmax+4),
     $          (iecpdat(l),l=lmax+5,maxl+5),nradp+iecpdat(maxl+6)
                npsp=npsp+1
                nelpsp=nelpsp+izcore
                do i=1,irad
                  write(icradp,'(i7,1x,e25.17,1x,e25.17)')
     $                       iradp(i),radp(2,i),radp(1,i)
                enddo
                nradp=nradp+irad
              endif
            enddo
          else
            do n=1,na
              call getatsymex(n,na,atsyex)
              if(reasy.eq.atsyex)then
                write(icpsp,'(i7,50(1x,i7))')(iecpdat(i),i=1,2),n,
     $          (nradp+iecpdat(l),l=4,lmax+4),
     $          (iecpdat(l),l=lmax+5,maxl+5),nradp+iecpdat(maxl+6)
                npsp=npsp+1
                nelpsp=nelpsp+izcore
                do i=1,irad
                  write(icradp,'(i7,1x,e25.17,1x,e25.17)')
     $                       iradp(i),radp(2,i),radp(1,i)
                enddo
                nradp=nradp+irad
              endif
            enddo
          endif
        endif
      endif
c...
c...  read next line
c...
      goto 100
c...
 999  continue
      end
c======================================================================
      subroutine psptmpsetup(icpsp,icradp,maxl)
      implicit real*8 (a-h,o-z)
c...
c...  MM August 2004
c...
c...  file setup for temporary storage of psp data
c...
c...  icpsp  file channel for psp data
c...  icradp file channel for radial data
c...  maxl   maximum angular momentum of psp we can treat
c...
      character*256 jobname,pspfile,radpfile
      character*100 char
c...
      call getchval('scrf',jobname)
      pspfile=jobname(1:len_trim(jobname))//'.psp'
      open(unit=icpsp,file=pspfile(1:len_trim(pspfile)),
     $     status='unknown',form='formatted')
      radpfile=jobname(1:len_trim(jobname))//'.radp'
      open(unit=icradp,file=radpfile(1:len_trim(radpfile)),
     $     status='unknown',form='formatted')
c...
c...  write file headers
c...
      write(icpsp,*)
      char(1:32) = '   Lmax  Izcore  Center   Local'
      do l=0,maxl
        write(char(32+8*l:39+8*l),1000) l
 1000   format(7X,I1)
      enddo
      char(33+8*(maxl+1):40+8*(maxl+1)) = '    End'
      write(icpsp,'(a)') char(1:39+8*(maxl+1))
      write(icpsp,*)
      write(icradp,*)
      write(icradp,*)'   Nr   ',
     $ '         Exponent         ','       Coefficient        '
      write(icradp,*)
c...
      end

c-----------------------------------------------------------------------

      subroutine get_psp_centre( iecpdat, maxlpsp, npsp, n, ic)

           ! returns the center of pseudopotential n

      implicit none

      integer maxlpsp, npsp, n, ic
      integer iecpdat(maxlpsp+6, npsp)

      if( n .le. 0 .or. n .gt. npsp ) then
        call nerror( 1, 'get_psp_centre', 'n is out of range',
     $              n, npsp)
      endif

      ic = iecpdat(3,n)

      end

c-----------------------------------------------------------------------

      subroutine realin(nchan,linea,noca,inica,ifica,lunmax,iend,icap)
c...
c...  this subroutine reads one line of characters from channel nchan.
c...  if nchan.eq.0 it reads from standard input channel.
c...  every line will contain at most lunmax characters.
c...  the line will be dived into noca fields, each field is defined
c...  as a sequence of characters other than ' '(blank).
c...  the index of the initial character of every field is stored in
c...  the inica array, and the index of the final character in the
c...  ifica array. so inica an ifica must be dimensioned at least
c...  lunmax/2. if the end of file is found, the variable iend is set
c...  to a value different from zero.
c...  if icap.ne.0 every small letter is rewritten as capital.
c...
c...  parameters:
c...
c...  name    type       contents
c...
c...  nchan  integer    chnannel to reas from (nchan.eq.0 correspond to
c...                    standard input)
c...  linea  char*127   in output will contain the line read from nchan
c...  noca   integer    in output will contain the number of fields
c...                    contained in linea
c...  inica  integer    array containing the index of the first chracter
c...                    of every field  (must be dimensioned at least
c...                    lunmax/2)
c...  ifica  integer    array containing the index of the last chracter
c...                    of every field  (must be dimensioned at least
c...                    lunmax/2)
c...  lumnax integer    maximum number of characters to be read
c...  iend   integer    in output it will be =0 if no errors occurs,
c...                    different from 0 if the end of file is found
c...  icap   integer    if icap.ne.0, every small letter is rewritten
c...                    as capital
c...
c...  questa subroutine legge una linea di caratteri dal canale nchan.
c...  se nchan=0 legge dal canale di input standard.
c...  ogni linea puo' contenere al massimo lunmax caratteri.
c...  la linea viene suddivisa in noca campi costituiti da sequenze
c...  di caratteri <> ' ', l'indice del carattere iniziale di ogni campo
c...  e' memorizzato nel vettore inica e il finale in ifica.
c...  inica ed ifica devono essere dimensionati almeno lunmax/2.
c...  se viene letto l'eof ritorna iend <> 0.
c...  se icap <> 0  ogni lettera minuscola viene riscritta in maiuscolo.
c...
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      character*127 linea
      dimension inica(1),ifica(1)
      character*1 bl,car
      data bl/' '/
c...
      iend=0
c...
c...  azzera inica e ifica.
c...
      do 5 i=1,lunmax/2
      inica(i)=0
5     ifica(i)=0
c...
c...  legge la linea di input in formato alfanumerico.
c...
      if(nchan.ne.0)then
      read(nchan,100,end=200)linea
      elsE
      read(*,100,end=200)linea
      endif
c...
c...  determina i caratteri iniziali e finali di ogni campo.
c...
      car=bl
      noca=0
      do 1 i=1,lunmax
      if(car.eq.bl.and.linea(i:i).ne.bl)then
      noca=noca+1
      inica(noca)=i
      endiF
      car=linea(i:i)
1     continue
c...
      if(noca.ne.0)then
      do 2 i=1,noca
      j=inica(i)
      do 3 k=j,lunmax
      if(linea(k:k).eq.bl)goto 4
3     continue
      k=lunmax+1
4     ifica(i)=k-1
2     continue
c...
c...  cambia da minuscolo a maiuscolo.
c...
      if(icap.ne.0)then
      do 10 i=1,noca
      do 10 j=inica(i),ifica(i)
      jo=ichar(linea(j:j))
      if(jo.ge.97.and.jo.le.122)linea(j:j)=char(jo-32)
10    continue
      endiF
      endif
c...
      return
200   iend=1
      return
c...
100   format(a)
      end
