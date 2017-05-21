      subroutine gaussinp(inpf)
      implicit real*8 (a-h,o-z)
c  this routine scans the input and determines if it is a Gaussian-style
c  input or not
c  Gaussian (Pople) style input is identified by having a pound (hash)
c  symbol as the first non-blank character of a line
c  if the input file (unit inpf) is non-Gaussian, nothing is done
c  if it is Gaussian, some preliminary input processing is performed
c  ARGUMENTS:
c  INPUT:
c  inpf : unit number of the input file (usually 30)
c  OUTPUT:
c  none but the file inpf may be changed
      parameter (linew=80)
      logical ginp,memc,chkc,rwfc,printc,uhf,dft,prelim,tight,zmat
c  stopit will be set .true. if any of the Pople-style input cards
c  starts with the string CONV (or conv). In this case, the program
c  will not continue
      logical stopit
      character*1 a(80),blank,hash,slash,card(5)
      character aa*80,title*80,text*80,job*256,bb*4,geomfile*256
      equivalence (a(1),aa)
      character chk*5,mem*5,rwf*5
      character chkfile*256, rwfile*256
      character*15 method,basis
      character*20 dftstr,thrcard,itercard,shiftcard,ishift
      character printcard*8
      character*80 cchrg
      character*4 unit,iter
      character*10 coorc,typec,optc
      parameter(nopt=30)
      character*4 opnames(nopt)
      character*256 chopval
      dimension ioptyp(nopt),iopval(3,nopt),ropval(3,nopt)
      dimension chopval(nopt),ifound(nopt)
      data blank,hash,slash/' ','#','/'/
      data mem,chk,rwf/'%mem=','%chk=','%rwf='/
      data opnames /'sp  ','forc','opt ','freq','nmr',
     1              'mp2 ','dyna','pop ','symm','xxxx',
     2              'geom','nosy','unit','gues','scf',
     3              'scfc','vshi','coor','optc','ts  ',
     4              'maxd','gdii','time','maxc','temp',
     5              'step','test','file','xxxx','card'/
c  ioptyp= 0=logical, 1= one integer, 2=two integers, 3= three int.
c          11=one real, 12=  two reals, 13=three reals, 21=character
      data ioptyp / 0,     0,     21,    0,     0,
     1              0,     0,     21,   11,     0,
     2              21,    0,     21,    21,    21,
     3              1,     11,    21,    1,     0,
     4              1,     1,     1,     1,     11,
     5              1,     0,     21,    0,     0/
      data card/'c','a','r','d',' '/
CPP
      ginp=.false.
      memc=.false.
      chkc=.false.
      rwfc=.false.
      printc=.false.
      printcard=' '
      dft=.false.
      dftstr=' '
      prelim=.false.
      zmat=.false.
      inp1=1
      inp2=2
c  determine if this is a conversion-only job
      stopit=.false.
 50   continue
c     jump-back address
        call blankit(a,linew)
        read(inpf,'(80a1)',end=70)  (a(k),k=1,linew)
        do i=1,linew
          if(a(i).ne.blank) then
            call lowercas(a,80)
            if(aa(i:i+3).eq.'conv') then
                stopit=.true.
              go to 70
            end if
          else
            go to 60
          end if
        end do
  60    continue
      go to 50
  70  continue
      rewind inpf
c  now decide if this is a gaussian-style input
 100  continue
c     jump-back address
      call blankit(a,linew)
      read(inpf,'(80a1)',end=1000) (a(k),k=1,linew)
      do i=1,linew
        if(a(i).ne.blank) then
          if(a(i).eq.hash) then
            If(a(2).eq.'c'.and.a(3).eq.'o'.and.a(4).eq.'m') go to 1000
            If(a(2).eq.'e'.and.a(3).eq.'n'.and.a(4).eq.'d') go to 1000
            ginp=.true.
            go to 1000
          else if(a(i).eq.'$') then
c -- cannot be gaussian-style if any line contains '$'
            go to 1000
          else
            exit
          end if
        else
          go to 200
        end if
  200   continue
      end do
      go to 100
 1000 continue
      if(.not.ginp) then
        rewind inpf
        return
      end if
c  do some preliminary input processing - eliminate leading blanks
      rewind inpf
 1100 continue
c  jump back
      call blankit(a,linew)
      read(inpf,'(300a1)',end=1000) (a(k),k=1,linew)
      call leadblan2(a,linew,len1)
      call rmblan(a,len1,len)
      call lowercas(a,len)
      if(aa(1:5).eq.mem) then
c memory allocation card
        read(aa(6:len),*) memsize
          if(memsize.gt.2000) memsize=memsize/1 000 000+1
        memc=.true.
        go to 1100
      end if
      if(aa(1:5).eq.chk) then
c CHK card - set the CHK option. This can differ slightly from Gaussian,
c  as PQS CHK affects only the file WE READ FROM
        chkfile=aa(6:len)
        lenchk=len-5
        chkc=.true.
        go to 1100
      end if
      if(aa(1:5).eq.rwf) then
c RWF card - set the SCR option.
        rwfile=aa(6:len)
        lenrw=len-5
        rwfc=.true.
        go to 1100
      end if
      if(a(1).eq.hash) then
c  read Gaussian Job card (route card)
CPP
c        print *,aa
c  the route card is the most important card
        istart=3
        if(a(2).eq.'p') then
          printc=.true.
          printcard=' PRINT=3'
          istart=4
        end if
c  search for the second  and third words
       call blankit(basis,15)
c  first locate the first non-blank character past the second
        do i=istart,linew
          if(a(i).ne.blank) then
c  determine the method
            do j=i,len
              if(aa(j:j).eq.slash) then
                method=aa(i:j-1)
                go to 2000
              end if
            end do
 2000       continue
            uhf=.false.
            if(method(1:1).eq.'u') uhf=.true.
            if(method(1:1).eq.'r'.or.method(1:1).eq.'u') then
              do k=2,15
                method(k-1:k-1)=method(k:k)
              end do
            end if
c extract the DFT potential from the method
            dft=.true.
            if(method(1:2).eq.'hf'.or.method(1:2).eq.'mp') dft=.false.
            if(method(1:3).eq.'hfs') dft=.true.
            if(dft) then
              dftstr=' DFTP='//method
            else
              dftstr=' '
            end if
            do k=j+1,linew
              if(a(k).eq.blank) then
c   determine the basis
                basis=aa(j+1:k-1)
                go to 3000
              end if
            end do
          end if
        end do
 3000   continue
        i1=0
        do i=k+1,len
          i1=i1+1
          a(i1)=a(i)
        end do
        len=len-k
      end if
c  For testing:
c      print *, 'Mem? ',memc
c      if(memc) print *, 'memory size in double words=',memsize
c      print *,'Chk? ',chkc
c      if(chkc) print *,'Check file=',chkfile(1:lenchk)
c      print *,'Rwf? ',rwfc
c      if(rwfc) print *,'RWF file=',rwfile(1:lenrw)
c      print *,'Printc? ',printc
c
c      print *, 'Method=',method,'uhf? ',uhf
c      print *, 'Basis=',basis
c      print *,'Remaining options=',aa(1:len)
 4000 continue
      open(inp1,status='scratch',form='formatted')
      write(inp1,'(80a1)') (card(k),k=1,5),(a(k),k=1,len)
      rewind inp1
      call readop2(inp1,nopt,opnames,ioptyp,iopval,ropval,
     1  chopval,ifound)
c      data opnames /'sp  ','forc','opt ','freq','nmr',
c     1              'mp2 ','dyna','pop ','symm','xxxx',
c     2              'geom','nosy','unit','gues','scf',
c     3              'scfc','vshi','coor','optc','ts  ',
c     4              'maxd','gdii','time','maxc','temp',
c     5              'step','test','file','xxxx','card'/
c check for continuation lines
      close(inp1,status='delete')
      if(method(1:3).eq.'mp2') ifound(6)=1
 4100 continue
      call blankit(a,linew)
      read(inpf,'(300a1)',end=1000) (a(k),k=1,linew)
      call leadblan2(a,linew,len1)
c      if(a(1).eq.'\r') a(1)=' '
      call rmblan(a,len1,len)
      call lowercas(a,len)
c      print *,'line=',aa
c  continuation found
      if(len.gt.0) then
        open(1,status='scratch',form='formatted')
        write(1,'(80a1)') (a(k),k=1,len)
        rewind 1
      call readop2(1,nopt,opnames,ioptyp,iopval,ropval,chopval,ifound)
          close(1,status='delete')
        goto 4100
      end if
c determine the level shift card
      if(ifound(17).gt.0) then
        write(ishift,'(f6.2)') ropval(1,17)
        call leadblan2(ishift,6,lshift)
        shiftcard=' LVSH='//ishift
      else
        shiftcard=' '
      end if
      call rmblan(shiftcard,20,lshift)
c blank card found - read title
      i1=0
 5000 continue
      call blankit(a,linew)
      read(inpf,'(300a1)',end=9999) (a(k),k=1,linew)
      call leadblan2(a,linew,len1)
      call rmblan(a,len1,len)
      if(len.gt.0) then
        i1=i1+1
        if(i1.eq.1) then
          title=aa
          lentit=len
        else
          text=aa
          lentxt=len
        end if
        goto 5000
      end if
c read charge+multiplicity
      read(inpf,*,end=9999) ichrg,mult
c copy the rest of the output to PQS input later
c  Generate PQS input
c      data opnames /'sp  ','forc','opt ','freq','nmr',
c     1              'mp2 ','dyna','pop ','symm','xxxx',
c     2              'geom','nosy','unit','gues','scf',
c     3              'scfc','vshi','coor','optc','ts  ',
c     4              'maxd','gdii','time','maxc','temp',
c     5              'step','test','file','xxxx','card'/
c  Do first the %MEM and FILE cards
      open(unit=inp2,status='scratch',form='formatted')
      rewind inp2
      if(memc) then
        if(memsize.lt.10) then
          write(inp2,6100) memsize
        else if(memsize.lt.100) then
          write(inp2,6200) memsize
        else
          write(inp2,6300) memsize
        end if
 6100 format('%MEM=',i1)
 6200 format('%MEM=',i2)
 6300 format('%MEM=',i3)
      end if
c  check card
      if(chkc.and.rwfc) then
        write(inp2,*) 'FILE SAVE='//chkfile(1:lenchk),
     1   ' CHK='//chkfile(1:lenchk),' SCR='//rwfile(1:lenrw)
      else if(chkc) then
        write(inp2,*) 'FILE SAVE='//chkfile(1:lenchk),
     1   ' CHK='//chkfile(1:lenchk)
      else if(rwfc) then
        write(inp2,*) 'FILE  SCR='//rwfile(1:lenrw)
      end if
      if(lentit.lt.1) lentit=1
      write(inp2,*) 'TITLE='//title(1:lentit)
c  if the title was 2 lines, write the 2nd line as TEXT
      if(i1.gt.1.and.lentxt.gt.0) then
        write(inp2,*)'TEXT '//text(1:lentxt)
      end if
c  if no GEOM card is found then assume Z matrix
c  if there is a NOSYM card, set symm=0
      symm=1.0d-5
      if(ifound(9).eq.1) then
        symm=ropval(1,9)
      end if
      unit='    '
      if(ifound(13).eq.1.and.
     1 (chopval(13)(1:2).eq.'au'.or.chopval(13)(1:4).eq.'bohr')) then
        unit='bohr'
      end if
      if(ifound(12).eq.1) SYMM=0.0d0
      if(ifound(28).eq.0) then
        geomfile=' '
      else
        geomfile=' FILE='//chopval(28)
      end if
      call rmblan2(geomfile,256,lgeom)
c  transform ichrg, the charge as an integer, into a char. string
      write(cchrg,'(i3)') ichrg
      call leadblan2(cchrg,4,lchrg)
      call rmblan2(cchrg,lchrg,l1chrg)
      if(ifound(11).eq.0) then
c     write(*,'(a4)') unit
        write(inp2,'(a20,f8.6,a8,a4,a6,i1,1x,a4,1x,a80)')
     1 'GEOM=ZMAT GEOP SYMM=',symm,' CHARge=',cchrg(1:l1chrg),
     2 ' MULT=',mult,unit,geomfile
        zmat=.true.
      else
        call lowercas(chopval(11),4)
        if(chopval(11)(1:4).eq.'cart') chopval(11)(1:4)='pqs '
        write(inp2,'(a20,f8.6,a8,a4,a6,i1,1x,a4,1x,a80)')
     1  'GEOM='//chopval(11)(1:4)//' GEOP SYMM=',symm,
     2   ' CHARge=',cchrg(1:l1chrg),' MULT=',mult,unit,geomfile
         if(chopval(11)(1:4).eq.'zmat') zmat=.true.
      end if
c  Copy the geometry
      if(ifound(28).ne.0) go to 7500
      nblank=0
 7000 continue
c  jump-back address
      read(inpf,'(80a1)',end=7500) (a(k),k=1,linew)
      call leadblan2(aa,linew,len1)
      call rmblan(aa,len1,lena)
      call lowercas(a,lena)
      if(aa(1:8).eq.'variable') then
        do i=1,lena
          aa(i:i)=' '
        end do
        lena=0
      end if
      if(lena.eq.0) nblank=nblank+1
      if (lena.eq.0) then
        if(nblank.eq.1.and.zmat) then
          write(inp2,*) 'VARIABLES'
        end if
          go to 7000
      end if
      if(lena.ne.0) then
        call leadblan2(aa,lena,lena1)
        lena=lena1
        bb(1:4)=aa(1:4)
        call lowercas(bb,4)
        if(bb.ne.'conv') write(inp2,*) aa(1:lena)
        go to 7000
      end if
c   geometry input is terminated by empty line or a Variables: card
c  if the basis is not very small or semiempirical, first put in the
c  3-21G basis, unless GUESS=READ was specified
 7500 continue
      if(ifound(14).ne.0) then
        if(chopval(14)(1:4).eq.'read') go to 8000
      else
        chopval(14)(1:4)='    '
      end if
      if(basis(1:3).eq.'sto'.or.basis(1:2).eq.'3-'.or.
     1  basis(1:3).eq.'min'.or.basis(1:4).eq.'semi') go to 8000
c Do first a 3-21 calculation
      if(ifound(14).gt.0.and.chopval(14)(1:4).eq.'read') go to 8000
        prelim=.true.
        write(inp2,*) 'BASIS=3-21G'
c        write(inp2,*) 'GUESS'
        chopval(14)(1:4)='read'
c  Determine if this is a DFT calculation
      write(inp2,*) 'SCF ITER=6 ',dftstr,shiftcard
 8000 continue
c  Now put in the real basis
      write(inp2,*) 'BASIS='//basis//printcard
c Put in the OPTIM and DYNAM cards, if any
c First gather the OPT options
c  In general, one must jump back by 4 cards
      jumpback=4
      if(basis(1:4).ne.'semi'.and.ifound(14).gt.0) then
        write(inp2,*) 'GUESS='//chopval(14)(1:4)
      end if
c  OPTIMIZATION
      if(ifound(3).gt.0) then
        typec=chopval(3)(1:4)
        if(chopval(3)(1:4).ne.'    ') then
          ifound(20)=1
          chopval(20)(1:4)=chopval(3)(1:4)
        end if
        coorc=' '
        if(ifound(18).gt.0) then
          coorc=' COOR='//chopval(18)(1:4)
        end if
        typec=' '
        if(ifound(20).gt.0) typec=' TYPE='//chopval(20)(1:4)
        optc=' '
        if(ifound(19).gt.0) then
          write(optc,*) iopval(1,19)
          call leadblan2(optc,10,lenl)
          call rmblan2(optc,lenl,lenll)
          optc=' OPTC='//optc(1:lenll)
        end if
        call rmblan2(coorc,10,lenc)
        call rmblan2(typec,10,lent)
        call rmblan2(optc,10,leno)
        if(lenc.eq.0) lenc=1
        if(lent.eq.0) lent=1
        if(leno.eq.0) leno=1
        write(inp2,*) 'OPTI'//coorc(1:lenc)//typec(1:lent)//
     1   optc(1:leno)//printcard
        if(ifound(18).gt.0.and.chopval(18)(2:4).eq.'zmat') then
          write(inp2,*) 'GEOM=ZMAT'
          jumpback=jumpback+1
        end if
c        write(inp2,*) 'GUESS=READ'
c  END OPT
      end if
c  DYNAM
      if(ifound(7).gt.0) then
        coorc=' '
        if(ifound(23).gt.0) write(coorc,*) iopval(1,23)
        call leadblan2(coorc,10,lenc)
        typec=' '
        if(ifound(25).gt.0) write(typec,'(f5.0)') ropval(1,25)
        call leadblan2(typec,10,lent)
        optc=' '
        if(ifound(24).gt.0) write(optc,'(i6)') iopval(1,24)
        call leadblan2(optc,10,lenl)
        call rmblan2(optc,lenl,lenll)
        write(inp2,*) 'DYNA STEP='//coorc//' MAXC='//optc//' TEMP='
     1    //typec
c        write(inp2,*) 'GUESS=READ'
c  END DYNAM
      end if
c
c  OPT or DYNAM
      if(ifound(3).gt.0.or.ifound(7).gt.0) then
        if(dft) then
          write(inp2,*) 'SCF LOCA=PIPEK'//printcard//dftstr
        else
          write(inp2,*) 'SCF LOCA=PIPEK'//printcard
        end if
        write(inp2,*) 'FORCE'
c  if NMR is specified together with dynam, do an NMR calculation in
c  each cycle. If NMR is specified w/OPTI, do only 1 NMR later
        if(ifound(7).gt.0.and.ifound(5).gt.0) then
          write(inp2,*) 'NMR'
          jumpback=jumpback+1
        end if
c        write(inp2,*) 'JUMP ',jumpback
        write(inp2,*) 'JUMP'
c if it is an OPT run, calculate geometry parameters
        if(ifound(3).gt.0) then
          write(inp2,*) 'GEOM GEOP'
        end if
        go to 9000
c  END OF OPT OR DYNAM
      end if
c  Determine first the highest property to be calculated
c  For this purpose, they are ranked as
c    sp<force<opt<pop<nmr<mp2<freq
c  DETERMINE SCF OPTIONS
      tight=.false.
      if(ifound(2).gt.0.or.ifound(3).gt.0.or.ifound(4).gt.0.
     1  .or.ifound(5).gt.0.or.ifound(6).gt.0) then
c  set the "tight" option
          tight=.true.
      end if
      if(ifound(15).gt.0) then
         if(chopval(15)(1:4).eq.'tigh') then
           tight=.true.
         end if
      end if
      if(.not.tight) then
        write(inp2,*) 'INTE THRE=7,9'
      end if
c  SCF
c  Put in the SCF card
      if(.not.tight) then
        thrcard=' THRE=4'
      else
        thrcard=' '
      end if
      if(ifound(16).gt.0) then
        write(iter,*) iopval(1,16)
        call leadblan2(iter,4,liter)
        itercard=' ITER='//iter
      else
        itercard=' '
      end if
      call rmblan(thrcard,20,lthre)
      call rmblan(itercard,20,liter)
      call rmblan(dftstr,20,ldft)
      write(inp2,*)
     1         'SCF LOCA=PIPEK'//itercard(1:liter)//
     1   ' '//thrcard(1:lthre)//' '//dftstr(1:ldft)//shiftcard(1:lshift)
c FORCE
      if(ifound(2).ge.1) then
        write(inp2,*) 'FORCE'
      end if
 9000 continue
c  Population analysis and NBO
c  do population analysis except in a dynamics run
      if(ifound(8).gt.0.and.ifound(7).eq.0) then
        if(chopval(8)(1:3).eq.'nbo') then
          write(inp2,*) 'NBO'
        else
          write(inp2,*) 'POP='//chopval(8)(1:4)
        end if
      end if
c  NMR
c  do NMR except in a dynam run
      if(ifound(5).gt.0.and.ifound(7).eq.0) then
        write(inp2,*) 'NMR'//printcard
      end if
c  do MP2 energies
      optc=' '
c  MP2 should not be a property - it is a method. It was put
c  accidentaly there. It can be called both as
c  # MP2/basis  or as
c  # HF/basis  MP2
      if(ifound(6).gt.0) then
        if(ifound(21).gt.0) then
          write(optc,*) iopval(1,21)
          call leadblan2(optc,10,lenl)
          call rmblan2(optc,lenl,lenll)
          write(inp2,*) 'MP2 MAXD='//optc(1:lenll)//printcard
        else
          write(inp2,*) 'MP2'//printcard
        end if
      end if
c  FREQ
c  make frequencies last
      if(ifound(4).gt.0) then
        optc=' '
        if(ifound(26).gt.0) then
          fdstep=dble(iopval(1,26))*1.8897269D-4
          write(optc,'(f10.6)') fdstep
          call leadblan2(optc,10,lenl)
          call rmblan2(optc,lenl,lenll)
          write(inp2,*) 'NUMH FDSTep='//optc(1:lenll)
        else
          write(inp2,*) 'NUMH'
        end if
        write(inp2,*) 'GEOM NOORient Print=1'
c        write(inp2,*) 'GUESS=READ'
        write(inp2,*) 'SCF '//dftstr
        write(inp2,*) 'FORCe'
c        write(inp2,*) 'JUMP 5'
        write(inp2,*) 'JUMP'
        write(inp2,*) 'GEOM print=1'
        write(inp2,*) 'FREQ'
      end if
 9999 continue
      rewind inp2
      rewind inpf
 9500 continue
      call getchval('jobname',job)
      call rmblan2(job,256,len)
      open(1,file=job(1:len)//'.pqs',form='formatted')
c
c     status='unknown' seems to be redundant here,
c     and gfortran does not like it
c
c     open(1,file=job(1:len)//'.pqs',form='formatted',
c    1        status='unknown')
      read(inp2,'(80a1)',end=9600) (a(k),k=1,linew)
      call leadblan2(aa,80,lena)
      call rmblan2(aa,lena,lenaa)
      write(1,'(80a1)') (a(k),k=1,lenaa)
      write(inpf,'(80a1)') (a(k),k=1,lenaa)
c        print *, (a(k),k=1,lenaa)
      go to 9500
 9600 continue
      close(1)
      close(2)
      rewind inpf
      if(stopit) then
        close(1)
        close(2)
        iout=igetival('iout')
        write(iout,*) 'Pople style input converter found CONV keyword'
        write(iout,*) 'Only input conversion was done, PQS input is in '
     2   ,job(1:len)//'.pqs'
        STOP
      end if
      end
c===================================================================
      subroutine readop2(inpf,nopt,opnames,ioptyp,iopval,
     $                   ropval,chopval,ifound)
c  Almost identical with readop1 but it does not zero ifound
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
      implicit real*8 (a-h,o-z)
      character*4 opnames(nopt),b
      character*256 chopval,c,lbl
      logical first,paren
      parameter (linew=80)
      character*1 a(linew),space,comma,lparen,rparen,eqsgn
      dimension ioptyp(nopt),iopval(3,nopt),ropval(3,nopt),chopval(nopt)
      dimension ifound(nopt)
c  this is the maximum width of a line
      data space,comma,lparen,rparen,eqsgn /' ',',','(',')','='/
c
c   if a wrong keyword is encountered on the first card at the first
c   position, do not backspace the input
      call izeroit(ifound,nopt)
      first=.false.
      izeropt=0
c I think "first" shows whether the parser has found a card already
c which makes sense. If not, it just reads in.
c I need to remove the error message if the very first option
c does not have an equal sign
      call blankit(a,linew)
      read(inpf,100,end=1200) (a(i),i=1,linew)
c  in readop1, do not skip the first word until the first space or comma or (
 100  format(80a1)
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
      izeropt=izeropt+1
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
            chopval(k)=c
        else if (ioptyp(k).ne.0) then
          call nerror(2,'readop2','Wrong ioptyp', ioptyp(k),0)
        end if
      end if
      first=.false.
c  look at the next option on this card
 900   continue
       go to 200
      ELSE
        first=.true.
        call blankit(a,linew)
        read(inpf,100,end=1200) a
        i=1
        go to 200
      ENDIF
 1200  if(first) then
         backspace inpf
       else
        call message('readop2','wrong keyword= '//b,0,0)
       end if
       end
c=============================================================
