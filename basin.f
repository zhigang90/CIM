c  MM August 2004  added pseudopotentials
c ......................................................................
c  JB July 25 1997    basis set data now written to <basis> file
c ......................................................................
C PP JUL 7, 97 I have removed the basis set printing for basis
c   sets which are read in.
C PP JUL 7  97 write error message if basis set file is missing, ln. 226,232-34
C PP May 29 97 do not go in basin with a stored basis if none is given
c PP line 118
C PP May 29 97 Basin should return if there are no functions
c added in basin, at label 190:   190 if (ncs.eq.0) return
C PP 05/23/97 Reset the integral thresholds to the default values
c    in any case, not only if the value is not defined. lines
c    commented out: 86-88, 90-92,94
c ......................................................................
      subroutine basrea(inp)
c  this routine reads the basis set

      use memory

      implicit real*8 (a-h,o-z)
      parameter (nopt=6)
      parameter (nbases=52,ipribas=30)
      character*4 options(nopt)
      character*8 basname,bases,blank,AtSymb
      character*256 chopval,filname,basdir
      character*256 baspath
      logical isthere
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     1 ifound(nopt)
c  bases is the basis set name, length is the length
c (i.e. number of card images), icol is the number of columns in the
c  card images. For most bases, this is smaller than the maximum 100
c  the primitive bases are: 1=sto-3g, 2=3-21g, 3=4-21g, 4=4-31g, 5=6-31g
c  6=6-311g, 7=star, 8=star5(5-comp.d),9=dstar (the second star),
c  10=d(for 6-311g), 11=2d(for 6311g),12=f (6-311g),13=3d(6-311g),
c  14=p (for 6-311g),15=2p(6311g),16=+, 17=++,
c  18=mini, 19=midi, 20=dz, 21=dzp (only polar.), 22=tz2p(polar only),
c  23=tzext(3d2p1f), 24=dgcdz, 25=dgctz, 26=dgcqz, 27=dgc5z,
c  24=dgc+, 25=dgcp, 26=dgc2p, 27=dgc3d2f,28=dgc4d3f2g
c  dgc stands for Dunning General Contraction
      dimension bases(nbases)
c  dummy is usually false. If true, the program lets you do a
c  calculation with no basis function on an atom
      logical dummy
      integer*8 ntwoel ! integer*8 (always) to hold no. of 2el integrals
      common /symm/nsym,nsy(7)
      character*40 sto3gb,g321b,g421b,g431b,g631b,g6311b,starb,star5b,
     1 dstarb,d6311b,d26311b,f6311b,d36311b,p6311b,p26311b,plusb,
     2 starpab
      common /sto3gc/sto3gb(176)
      common /g321c/ g321b(144)
      common /g421c/g421b(44)
      common /g431c/g431b(98)
      common /g631c/g631b(234)
      common /g6311c/g6311b(66)
      common /starc/starb(32)
      common /star5c/star5b(32)
      common /dstarc/dstarb(4)
      common /d6311c/d6311b(30)
      common /d26311c/d26311b(45)
      common /f6311c/f6311b(34)
      common /d36311c/d36311b(68)
      common /p6311c/p6311b(4)
      common /p26311c/p26311b(4)
      common /plusc/plusb(30)
      common /plussc/plussb(4)
      common /starpac/starpab(22)
      data options/'basi','file','next','dumm','prin','gene'/
      data bases
     1 /'sto-3g  ','3-21g   ','4-21g   ','4-31g   ','6-31g   ',
     2 '3-21g(*)','4-21g*  ','4-21g(*)','4-31g*  ','4-31g(*)',
     3 '6-31g*  ','6-31g** ','6-311g  ','6-311gd ','6-311gdp',
     4 '6-311g2d','4-21+g  ','4-31+g  ','6-31+g  ','6-31+g* ',
     5 '6-311+g ','4-21+g* ','4-31+g* ','6-31+g* ','6-311+gd',
     6 '6-311+gp','6-311gex','mini    ','midi    ','dz      ',
     7 'dzp     ','dz2p    ','tz      ','tzp     ','tz2p    ',
     8 'tz3p    ','qz3p    ','dgcdz   ','dgctz   ','dgcqz   ',
     9 'dgcdz+  ','dgctz+  ','dgcqz+  ','dgcdzp  ','dgctzp  ',
     x 'dgcdz2p ','dgctz2p ','dgcqz2p ','dgctz3p ','dgcqz3p ',
     1 'dgc5zext','xxxxxxx '/
      data blank/'        '/
c
c  standard is one of the internally stored basis sets (e.g. 6-31G*)
c
c  file is an external filename on which the basis set info resides
c
c  next means that the basis is given following this card
c  these options are not mutually exclusive. However, the user
c  should not specify the same basis function twice. It is legal
c  to have all three present: the basis set will be the sum of
c  all three
      data ioptyp/21,21,0,0,1,0/
c
c  pseudopotential stuff:
c
c  maxl  is the maximum angular momentum we can treat
c  icpsp and icpradp are unit numbers for psp temporary storage
c
      parameter (maxl=6)
      parameter (icpsp=55,icradp=56,icbasinfo=57)
c
      call getival('iout',iout)
c .......................................................
c  set default for instability
      call setival('stab',+1)       ! assume everything is STABLE
c .......................................................
c  set default integral thresholds
      thres=1.0d-10
      thres1=1.0d-7
c   integral thresholds: main threshold first, secondary follows
      call setrval('ithr',thres)
      call setrval('ith1',thres1)
c ........................................................
c  integral printing option
      call setival('iprn',0)
      call getival('inuc',inuc)
c  reserve memory for the basis data and for the block inx
c  contraction info (integers)
c  (reserve tentatively a lot of memory)
      call getmem(78 000,ibas)
c  (good for 6000 primitive shells, i.e. about 24000 functions)   ??
      call getint(5000,ictr)
c
c  zero out the reading buffers, to avoid uninitialized crap
c  
      call zeroit(bl(ibas),78 000)
      call izeroit(bl(ictr),5000)
c
      call setival('ibas',ibas)
      na=igetival('na')
      dummy=.false.
c   ictr is the beginning point in inx of the contraction info (old inx)
      call setival('ictr',ictr)
c  (good for 2000 contracted shells)
c   zero out the basis set data
      nbf=0
      nsh=0
      ncf=0
      ncs=0
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,nopt,options,ioptyp,iopval,ropval,chopval,ifound)
c
c   initialize pseudopotential counters:
c
c    npsp    total no. of pseudopotentials
c    nradp   total no. of radial terms
c    nelpsp  total no. of electrons treated by psp
c
      npsp=0
      nradp=0
      nelpsp=0
c
c   setup temporary storage files for pseudopotential data
c
      call psptmpsetup(icpsp,icradp,maxl)
c
c  print level. Note the PRINt is equivalent to PRINT=1
c  print levels higher  than 1 are not used currently
      if(ifound(5).gt.0) then
        lprint=1
        if(iopval(1,5).gt.1) lprint=iopval(1,5)
      else
        lprint=0
      end if
c
c -- first check we have a standard basis input
      IF(chopval(1)(1:4).NE.'    ') THEN
c  test if this basis set is in the basis library
      call getchval('BASDIR',basdir)
      call rmblan2(basdir,256,len)
      filname=chopval(1)
      call rmblan2(filname,256,len1)
      call lowerca2(filname,len1)
c  massagefn massages the file name
      call massagefn(filname,len1)
      len2=len+len1
      if(len2.gt.256) then
        call nerror(1,'basrea',
     1  'basis set full pathname exceeds 256 characters',len,len1)
      end if
      baspath=basdir(1:len)//filname(1:len1)//'.bas'
      len2=len2+4
c       print *, baspath(1:len2)
      iunit=45
      inquire(file=baspath(1:len2),exist=isthere)
      if(isthere) then
        open(iunit,file=baspath(1:len2), status='old')
        rewind iunit
c  ibastyp is the basis number for internally stored bases.
c  It is 0 for basis sets read in. In this case, the basis is printed.
        call basin(iunit,  blank,  lprint, bl(inuc),nbf,
     $             nsh,    ncf,    ncs,    bl(ibas),bl(ictr),
     $             npsp,   nradp,  nelpsp, maxl,   icpsp,
     $             icradp, na)
        write(iout,*)'Basis set ',filname(1:len1),' from library'
c
c  look for pseudopotential input in the library file
c
        rewind iunit
        call rpsplib(iunit,icpsp,icradp,blank,na,npsp,nradp,nelpsp,maxl,
     $               lprint,iout)
        close (iunit,status='keep')
      else
c  standard basis - stored internally
        basname=chopval(1)(1:8)
        filname = basname
        if(ifound(1).eq.1.and.basname.ne.blank) then
          call lowerca2(basname,8)
          ibastyp=0
          do 100 i=1,nbases
            if (basname.eq.bases(i)) then
              ibastyp=i
              go to 105
            end if
 100      continue
 105      iscrf=99
          open(iscrf,err=110,status='unknown')
          rewind iscrf
 110      continue
         if(ibastyp.ne.0)
     1       write(iout,*) 'Internally stored basis ',basname
         if (ibastyp.eq.0) then
            call nerror(2,'basrea','no such basis '//basname,ibastyp,
     1      nbases)
          else if (ibastyp.eq.1) then
            call sto3g
            call copyinf(sto3gb,iscrf,40,176)
          else if (ibastyp.eq.2) then
            call g321
            call copyinf(g321b,iscrf,40,144)
          else if (ibastyp.eq.3) then
            call g421
            call copyinf(g421b,iscrf,40,44)
          else if (ibastyp.eq.4) then
            call g431
            call copyinf(g431b,iscrf,40,98)
          else if (ibastyp.eq.5) then
            call g631
            call copyinf(g631b,iscrf,40,234)
          else if (ibastyp.eq.6) then
            call g321
            call starpa
            call copyinf(g321b,iscrf,40,144)
            call copyinf(starpab,iscrf,40,22)
          else if (ibastyp.eq.7) then
            call g421
            call star5
            call copyinf(g421b,iscrf,40,44)
            call copyinf(star5b,iscrf,40,32)
          else if (ibastyp.eq.8) then
            call g421
            call starpa
            call copyinf(g421b,iscrf,40,44)
            call copyinf(starpab,iscrf,40,22)
          else if (ibastyp.eq.9) then
            call g431
            call star
            call copyinf(g431b,iscrf,40,44)
            call copyinf(starb,iscrf,40,32)
          else if (ibastyp.eq.10) then
            call g431
            call starpa
            call copyinf(g431b,iscrf,40,44)
            call copyinf(starpab,iscrf,40,22)
              else if (ibastyp.eq.11) then
                call g631
                call star
                call copyinf(g631b,iscrf,40,234)
                call copyinf(starb,iscrf,40,32)
              else if (ibastyp.eq.12) then
                call g631
                call star
                call dstar
                call copyinf(g631b,iscrf,40,234)
                call copyinf(starb,iscrf,40,32)
                call copyinf(dstarb,iscrf,40,4)
              else if (ibastyp.eq.13) then
                call g6311
                call copyinf(g6311b,iscrf,40,66)
              else if (ibastyp.eq.14) then
                call g6311
                call d6311
                call copyinf(g6311b,iscrf,40,66)
                call copyinf(d6311b,iscrf,40,30)
                else if (ibastyp.eq.15) then
                call g6311
                call d6311
                call p6311
                call copyinf(g6311b,iscrf,40,66)
                call copyinf(d6311b,iscrf,40,30)
                call copyinf(p6311b,iscrf,40,4)
              else if (ibastyp.eq.16) then
                call g6311
                call d26311
                call copyinf(g6311b,iscrf,40,66)
                call copyinf(d26311b,iscrf,40,45)
              else if (ibastyp.eq.17) then
              else if (ibastyp.eq.18) then
              else if (ibastyp.eq.19) then
              else if (ibastyp.eq.20) then
              else if (ibastyp.eq.21) then
              else if (ibastyp.eq.22) then
              else if (ibastyp.eq.23) then
              else if (ibastyp.eq.24) then
              else if (ibastyp.eq.25) then
            end if
            rewind iscrf
            call basin(iscrf,  blank,  lprint, bl(inuc),nbf,
     $                 nsh,    ncf,    ncs,    bl(ibas),bl(ictr),
     $                 npsp,   nradp,  nelpsp, maxl,    icpsp,
     $                 icradp, na)
           close(iscrf,status='delete')
        end if
      end if
      ENDIF
cc
      IF(ifound(2).eq.1) THEN
        call rmblan(chopval(2),256,len)
        iunit=45
        open(iunit,file=chopval(2)(1:len),err=200,status='old')
        rewind iunit
        call basin(iunit,  blank,  lprint, bl(inuc),nbf,
     $             nsh,    ncf,    ncs,    bl(ibas),bl(ictr),
     $             npsp,   nradp,  nelpsp, maxl,    icpsp,
     $             icradp, na)
c
c  look for pseudopotential input in the library file
c
        rewind iunit
        call rpsplib(iunit,icpsp,icradp,blank,na,npsp,nradp,nelpsp,maxl,
     $               lprint,iout)
      write(iout,*) 'Basis set augmented from file ',chopval(2)(1:len)
      ENDIF
      go to 300
  200 call nerror(3,'basrea',
     1 'basis set file '//trim(chopval(2))//' does not exist',0,0)
      stop
  300 continue
cc
      IF(ifound(3).eq.1) THEN
c -- two options
c    (1) additional basis functions read direct from input file
c    (2) basis for special atoms read from file or standard basis library
c
  310   Continue
        Call bastyp2(inp,iout,AtSymb,filname,len)
        If(AtSymb.Eq.blank) Then
         iunit=inp
         call basin(iunit,  blank,  lprint, bl(inuc),nbf,
     $              nsh,    ncf,    ncs,    bl(ibas),bl(ictr),
     $              npsp,   nradp,  nelpsp, maxl,    icpsp,
     $              icradp, na)
         write(iout,*) 'Basis set augmented from input'
        Else If(AtSymb.Ne.'Finish  ') Then
         iunit=45
         open (iunit,file=filname(1:len),err=400,status='old')
         rewind iunit
         call basin(iunit,  AtSymb, lprint, bl(inuc),nbf,
     $              nsh,    ncf,    ncs,    bl(ibas),bl(ictr),
     $              npsp,   nradp,  nelpsp, maxl,    icpsp,
     $              icradp, na)
c
c  look for pseudopotential input in the library file
c
        rewind iunit
        call rpsplib(iunit,icpsp,icradp,AtSymb,na,npsp,nradp,nelpsp,
     $               maxl,lprint,iout)
         close (iunit,status='keep')
        Else
         Go To 500
        EndIf
        Go To 310
  400   call nerror(4,'basrea',
     1   'basis set file '//filname(1:len)//' does not exist',0,0)
        stop
      ENDIF
  500 Continue
cc
      if(ifound(4).eq.1) dummy=.true.
cc
c--------------------------------------------------------------------

                ! temporarily dump the basis info into a file

      open( unit=icbasinfo, status='scratch', form='unformatted')
      call dump_basinfo(bl(ictr),bl(ibas),ncs,nsh,icbasinfo)

                ! deaallocate ictr and ibas

      call retmem(2)
c     call retint(1)
      
                ! now reallocate the right amount of memory

      call getmem(13*nsh,ibas)
      call setival('ibas',ibas)
      call getint(12*ncs,ictr)
      call setival('ictr',ictr)

                ! read back basis set info

      call read_basinfo(bl(ictr),bl(ibas),ncs,nsh,icbasinfo)
      close( unit=icbasinfo)

      call setival('ncs',ncs)
      call setival('nbf',nbf)
      call setival('ncf',ncf)
      call setival('nsh',nsh)
c  check for atoms without basis functions
c  ncenter is a temporary array used to check if all centers have basis
c  functions on them
      call getint(na,ncenter)
      if(.not.dummy) call chkbasis(na,bl(ncenter),ncs,bl(ictr),
     1                             bl(inuc))
      call retmem(1)
c ................................................................
c
c  check if there are any generally contracted shells in basis
c
      call check4lgc(ncs,bl(ictr),lshells,lgcshell,ltype_gc,ldeep_gc)
c
c  default is now to segment all generally contracted basis sets
c  i.e., split into separate basis functions
c  (Do NOT do this if GENEral keyword is found)      ! JB  March 2011
c
      If(ifound(6).GT.0) THEN
c
        if(lshells.gt.0 .and. lgcshell.gt.0) then
          write(iout,*)'                                              '
          write(iout,*)'               *** CAUTION ***               '
          write(iout,*)
          write(iout,*)
     *    '   There are ',lshells,' l-shells and ',lgcshell,' gc-shells'
          write(iout,*)
     *                 '   program cannot handle a basis set containing'
          write(iout,*)'   both: L-shells and General Contracted shells'
          write(iout,*)'   L-shells have been segmented into S,P-shells'
          call trans_l_2_sp(.true.,bl,bl(ictr),ictr,lshell,'basiss')
          call getival('ncs ',ncs)
          call getival('ibas',ibas)
          call getival('nsh',nsh)
          call getival('ictr',ictr)
        endif
cc
      ELSE
cc
        if(lgcshell.gt.0) then
          if(lshells.gt.0)
     $     call trans_l_2_sp(.true.,bl,bl(ictr),ictr,lshell,'basiss')
          call trans_gc_2_sg(.true.,bl,bl(ictr),ictr,nogcsi,'basiss')
          call getival('ncs ',ncs)
          call getival('ibas',ibas)
          call getival('nsh',nsh)
          call getival('ictr',ictr)
          write(iout,*)
         write(iout,*)'   General Contracted shells have been segmented'
        endif
      ENDIF
c ................................................................
      ntwoel=ncf*(ncf+1)/2
      ntwoel=ntwoel*(ntwoel+1)/2
      write (iout,1000) nbf,nsh,ncf,ncs,ntwoel
 1000 format(/1x,i0,' gaussians ',i0,' shells ',i0,' contr. gaussians ',
     1  i0,' contr. shells',/1x,i0,' integrals (less symmetry or',
     2  ' neglect)',/)
c ................................................................
c  reserve memory for the array giving the relation between the
c  old and new orders of the contracted functions. This array
c  gives, at the i-th position, the serial number of the basis function
c  i in the new order
      call getint(ncf,ireor)
      call setival('reor',ireor)
      call getint(ncs,iswap)
      call getint(12*ncs,ncso)
c ................................................................
c  write basis set data to <control> file
c  (do before reorder; use "iswap" and "ireor" as scratch)
      if(ifound(2)+ifound(3).gt.0) filname = 'nonstandard'
      CALL wrbasis(ncf,    ncs,    nsh,    npsp, bl(ireor),
     $             bl(ictr),bl(ibas),filname, icpsp)
c ................................................................
c  now read basis back - this is a desperate attempt to ensure
c  that <reorder> produces the same Wolinski ordering throughout
c
      CALL RdBackBasis(na,ncs,nsh,bl(ictr),bl(ibas))
c ................................................................
      call normaliz(ncs,bl(ictr),bl(ibas))
      call getmem(3*na,ix)
      call getmem(na,ian)
      call getnucdat(na,bl(ix),bl(ian))
      call reorder(ncs,bl(ictr),bl(ireor),bl(iswap),bl(ncso),
     $             bl(ian))
      call retmem(2)
c
c ................................................................
c2002 KW check for possible numerical instability in integral calc.
c
      call check4stab(na,ncs,bl(ictr),bl(ibas),bl(inuc),istab)
      If(istab.eq.0) Then
        write(iout,*) 'Basis might be numerically unstable',
     *                ' - check results carefully'
        istab = +1
      Else If(istab.eq.-1) Then
        write(iout,*) 'Basis may be numerically unstable',
     *                ' - integral stability switched on'
      EndIf
      call setival('stab',istab)
c ................................................................
c     call retint(2)
      call symfunc(.false.)
c
c   if pseudopotential were read in, setup the storage and
c   variables for pseudopotential calculations
c
      if(npsp.ne.0)then
        call getint(npsp*(maxl+6),iiecpdat)
        call getint(nradp,inrrad)
        call getmem(nradp,iarad)
        call getmem(nradp,icrad)
c
c   read in pseudopotential data from temporary files
c
        call psptmpread(icpsp,icradp,npsp,nradp,nelpsp,bl(iiecpdat),
     $              bl(inrrad),bl(iarad),bl(icrad),maxl)

           ! update atomic charges. need to set some
           ! variables for routine setnucchg

        call setival('npsp',npsp)
        call setival('iiecpdat',iiecpdat)
        call setival('maxlpsp',maxl)
        call setnucchg( 62, bl(inuc), na )

c
c  now we need to update the <coord> file with the corrected
c  atomic charges
c
        call pspupdatech(na,bl(inuc))
c
c  update total charge and nuclear repulsion energy
c
        call getival('ndum',ndum)
        call getival('mass',iadr2)
        call nucparam(na-ndum,bl(inuc),bl(iadr2))
        call calcnuc(na-ndum,bl(inuc),enuc)

        write(iout,2000)nelpsp,npsp,nradp
2000    format(/,' Pseudopotentials will be used:',/,
     $         i5,' electrons simulated by ',i5,
     $         ' pseudopotentials (',i5,' radial terms)'/)
        write(iout,2001)enuc
2001    format(' Nuclear repulsion energy after psp correction: ',
     $         f16.9,' au'/)
c
c  update no. of alpha and beta electrons
c
        call getival('Multip',imult)
        call getrval('charge',charg)
        call GetNAB(na-ndum,bl(inuc),Charg,IMult,NAlpha,NBeta)
        call pspupdateab(nalpha,nbeta)
        call setival('NAlpha',NAlpha)
        call setival('NBeta',NBeta)
c
c  store remaining psp parameters into depository
c
        call setival('nelpsp',nelpsp)
        call setival('nradp',nradp)
        call setival('inrrad',inrrad)
        call setival('iarad',iarad)
        call setival('icrad',icrad)
      endif
c
c  delete the (temporay) pseudopotential files
C
      close(unit=icpsp,status='delete')
      close(unit=icradp,status='delete')
c
c   this has to be always set (npsp=0 in a full electron run)
c
      call setival('npsp',npsp)
c
      end
c============================================================
      subroutine copyinf(array,ifile,linew,numli)
      character*(*) array(numli)
      character*1 line(100)
      do 200 i=1,numli
      call blankit(line,100)
      read (array(i),100,end=200) (line(k),k=1,linew)
 100  format(100a1)
      write(ifile,100) (line(k),k=1,linew)
 200  continue
      end
c=============================================================
      subroutine massagefn(filname,len)
      character*(*) filname
c   this routine replaces ** with (d,p); * with (d)
c   and then replaces opening parentheses with a hyphen and
c   deletes closing parentheses and commas
c  Arguments:
c   - filname: file name (input and output)
c   - len    : its length (input and output)
      character*256 name,blank
      character*1 star,opar,cpar,comma,hyphen
      data star,opar,cpar,comma,hyphen/'*','(',')',',','-'/
      name=filname
      len1=len
      ii=1
      do i=2,len1
        if(name(i-1:i).ne.star//star) then
          ii=ii+1
          if(name(i:i+1).eq.star//star) then
            filname(ii:ii+4)='(d,p)'
            ii=ii+4
          else if (name(i:i).eq.star) then
            filname(ii:ii+2)='(d)'
            ii=ii+2
          else
            filname(ii:ii)=name(i:i)
          end if
        end if
      end do
      name=filname
      call blankit(filname,256)
      len1=ii
      ii=0
      do i=1,len1
        if(name(i:i).ne.cpar.and.name(i:i).ne.comma) then
          ii=ii+1
          if(name(i:i).eq.opar) then
            filname(ii:ii)=hyphen
          else
            filname(ii:ii)=name(i:i)
          end if
        end if
      end do
      len=ii
      end

c==========================================================

      subroutine dump_basinfo(inx,basdat,ncs,nsh,icbasinfo)

               ! dumps the basis set info into channel icbasinfo

      implicit none

      integer ncs, nsh, icbasinfo
      integer inx(12,ncs)
      real*8 basdat(13,nsh)

      rewind(icbasinfo)
      write(icbasinfo)inx
      write(icbasinfo)basdat

      end

c==========================================================

      subroutine read_basinfo(inx,basdat,ncs,nsh,icbasinfo)

               ! reads the basis set info from channel icbasinfo

      implicit none

      integer ncs, nsh, icbasinfo
      integer inx(12,ncs)
      real*8 basdat(13,nsh)

      rewind(icbasinfo)
      read(icbasinfo)inx
      read(icbasinfo)basdat

      end

c==========================================================
      subroutine basin(inp,    AtSymb, lprint, xa,     nbf,
     $                 nsh,    ncf,    ncs,    basdat, inx,
     $                 npsp,   nradp,  nelpsp, maxl,   icpsp,
     $                 icradp, na)
c    this routine reads in the basis set info
c    this is an incremental read, i.e. it can be called several times
c     parameters
c     inp (in): input file
c     AtSymb (in): symbol for atom if basis for just ONE atom type
C                  is to be considered; otherwise blank
c     lprint (in): print level. If > 0, the basis is printed
c     xa (in): nuclear info
c     nbf (in and out): total no. of primitive functions
c     nsh (in and out): total no. of primitive shells
c     ncf (in and out): total no. of contracted functions
c     ncs (in and out): total no. of contracted shells
c     basdat (out): an array which holds the primitive basis set info
c     inx (out): integer general storage for basis set info
c     npsp (in and out) total no. of pseudopotentials
c     nrad (in and out) total no. of radial terms in pseudopotentials
c     nelpsp (in and out) total no. of electrons treated by pseudopot.
c     maxl (in)  maximum psp angular momentum that can be treated
c     icpsp (in) channel number for temporary psp storage
c     icradp (in) ditto
c     na (in)  number of atoms
c
      implicit real*8 (a-h,o-z)
c
c     .... maximum 30 functions per contraction
c
      dimension xa(5,na), basdat(13,*), iat(na), eta(30), cs(30),
     1          cgen(30,9),cp(30),inx(12,*)
      dimension nfun(13)
      character*4 for,blank,word,atyp,conti,ftyp(13)
      character*8 AtSymb,xnam,xat,Symb,blank8
      character*1 atchar,chnum(0:9)
      character*1 line(100)
      logical ltype,found
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
c *** on the ibm this is needed
c
      equivalence(xxnam,xnam)
      equivalence (xxat,xat)
      data for,blank,conti,ftyp /'for ','    ','c   ','s   ','p   ',
     1 'l   ','d   ','d6  ','f   ','f10 ','g15 ','h21 ','i28 ',
     2 'g   ','h   ','i   '/
      data nfun/1,3,4,5,6,7,10,15,21,28,9,11,13/
      data chnum/'0','1','2','3','4','5','6','7','8','9'/
      Data blank8/'        '/
c
      ltype = .false.
      found = .False.
      call getival('iout',iout)
      if(lprint.gt.0) write(iout,210)
  210 format (1x,15hbasis functions,/)
c     search for the first 'for' card, disregard cards before it
    5 read (inp,220,end=190) word
  220 format (a4)
      call lowerca2(word,4)
      if (word.eq.for) then
        backspace inp
        go to 10
      else
        go to 5
      end if
   10 continue
      If(found.AND.AtSymb.NE.blank8) GO TO 190
c
c -- At this point, if AtSymb is not blank, we are attempting
c -- to add a basis to an atom with a "special symbol"
      If(AtSymb.NE.blank8) Call GetSymbM(AtSymb,Symb)   ! get special symbol
c  xxnam is the same as xnam but is of type real
      read (inp,230) word,xnam
  230 format (a4,6x,a4)
      call lowerca2(xnam,4)
c -- if appropriate, add special symbol to atom
      If(AtSymb.NE.blank8) Call AddSymbM(Symb,xnam)
   30 nata=0
      do 40 j=1,na
        xxat=xa(5,j)
c  xat is equivalenced to xxat
        jchar=0
        do ichar=1,8
          atchar=xat(ichar:ichar)
          if(atchar.eq.' ') go to 32
c    delete numerals from the atom name
          do num=0,9
          if(atchar.eq.chnum(num)) go to 32
          end do
          jchar=jchar+1
          xat(jchar:jchar)=atchar
  32      continue
        end do
        do kchar=jchar+1,8
        xat(kchar:kchar)=' '
        end do
        if (xxat.ne.xxnam) go to 40
        If(AtSymb.NE.blank8.AND.AtSymb.NE.xnam) go to 40
c
c *** on the ibm-style machines the preceding statement should be
c *** if (xa(5,j).ne.xxnam(i)) go to 40 because of the difficulty
c *** of comparing characters with reals
c
        nata=nata+1
        iat(nata)=j
        found = .True.
   40 continue
      if (nata.eq.0) go to 5
      if(lprint.gt.0) write(iout,45) xnam,(iat(kk),kk=1,nata)
 45   format(' For atom: ',a4,2x,10i4,/,(17x,10i4))
   60 read (inp,260,end=190,err=190) atyp
      backspace inp
  260 format (a4,6x,11f10.6)
      call lowerca2(atyp,4)
      if (atyp.eq.for) go to 10
      if(atyp.eq.'psp ')then  ! user pseudopotential input
        call rpspinp(inp,icpsp,icradp,xnam,na,npsp,nradp,nelpsp,maxl,
     $               lprint,iout)
        goto 60
      endif
      if(atyp.eq.'ecp ')then  ! user pseudopotential input
        call rpspinp(inp,icpsp,icradp,xnam,na,npsp,nradp,nelpsp,maxl,
     $               lprint,iout)
        goto 60
      endif
      do 70 j=1,13
         if (atyp.eq.ftyp(j)) go to 80
   70 continue
      go to 190
   80 continue
      mxtyp=j
c .........................................................................
c  WARNING -  spherical i-functions currently DO NOT work
      If(mxtyp.EQ.13) call nerror(1,'basin',
     $  'Sorry - cannot use spherical i-functions at this time',0,0)
c .........................................................................
c
c     the type specification may contain only the characters s,p,l,d,d6,
c     f,f10,g15,h21,i28,g,h,i, blank, or the word for
c     end-of-basis set information is assumed if another character is en
c     countered
c
      if(mxtyp.eq.3) then
        ltype=.true.
        ngcmax=0
        ngc=0
c  l type function -no general contraction
        read (inp,260) atyp,eta(1),cs(1),cp(1)
        if(cs(1).eq.zero) cs(1)=one
        if(cp(1).eq.zero) cp(1)=one
ccc        if (abs(cp(1)).le.acc) cp(1)=one
      else
c      other types -gen. contraction possible
        do 85 igc=1,9
          cgen(1,igc)=zero
 85     continue
        read (inp,260,end=86) atyp,eta(1),cs(1),(cgen(1,igc),igc=1,9)
          if(cs(1).eq.zero) cs(1)=one
 86     continue
      end if
c     if (abs(cp(1)).lt.acc.and.abs(cs(1)).ge.acc) cp(1)=cs(1)
      if (eta(1).le.zero) then
        call nerror(2,'basin','non-positive exponent',ncs,1)
      end if
      ngc=0
      if(.not.ltype) then
c     count last general contraction
        ngcmax=0
        do igc=1,9
          if(abs(cgen(1,igc)).gt.zero) ngcmax=igc
        end do
      end if
      isum=nfun(mxtyp)
c        write (iout,270) atyp
      if(lprint.gt.0) write (iout,270) atyp
      if(ltype) then
c        write (iout,320) eta(1),cs(1),cp(1)
        if(lprint.gt.0) write (iout,320) eta(1),cs(1),cp(1)
  320   format (1x,'eta=',f14.7,3x,'c(s)=',f14.7,3x,'c(p)=',f14.7)
      else
c        write (iout,321) eta(1),cs(1),(cgen(1,igc),igc=1,ngc)
        if(lprint.gt.0) write (iout,321) eta(1),cs(1),(cgen(1,igc),
     1     igc=1,ngc)
  321   format (1x,'eta=',f14.7,3x,'c=',6f14.7)
      end if
      i=1
  100 read (inp,280,end=140) word
  280 format (a4)
      backspace inp
      if (word.ne.blank.and.word.ne.conti) go to 140
c  if word = blank, check if the whole line is empty
      if(word.eq.blank) then
        read(inp,'(100a1)',end=62) line
  62    backspace inp
        call rmblan(line,100,lgt)
        if(lgt.eq.0) go to 140
      end if
c
      i=i+1
      if (i.gt.30) then
        call nerror(3,'basin','maximum contr. length 30 exceeded',i,i)
      end if
  110 continue
      if (mxtyp.eq.3) then
        read (inp,300) eta(i),cs(i),cp(i)
  300   format (10x,11f10.6)
        if(cs(1).eq.zero) cs(1)=one
        if(cp(1).eq.zero) cp(1)=one
ccc        if (abs(cp(i)).eq.zero) cp(i)=one
      else
        do 115 igc=1,9
          cgen(i,igc)=zero
 115    continue
        read (inp,300,end=116) eta(i),cs(i),(cgen(i,igc),igc=1,9)
 116    continue
      end if
c count number of gen. contraction coeff. for this primitive
      nngc=0
      do igc=1,9
        if(abs(cgen(i,igc)).ne.zero) nngc=igc
      end do
      if(nngc.gt.ngcmax.and.mxtyp.ne.3) ngcmax=nngc
c      if (abs(cs(i)).eq.zero) cs(i)=one
      if (eta(i).gt.0.0d0) go to 120
      call nerror(4,'basin','non-positive exponent',ncs,i)
  120 continue
      if(mxtyp.eq.3) then
c        write (iout,320) eta(i),cs(i),cp(i)
        if(lprint.gt.0) write (iout,320) eta(i),cs(i),cp(i)
        ngcmax=0
      else
c        write (iout,321) eta(i),cs(i),(cgen(i,igc),igc=1,nngc)
        if(lprint.gt.0) write (iout,321) eta(i),cs(i),(cgen(i,igc),
     1     igc=1,nngc)
      end if
      go to 100
  140 continue
c  massage the coefficients etc
c if there are no general contractions then a missing (zero) coefficient
c is interpreted as 1
      nprim=i
      ngc=ngcmax
      if(ngc.eq.0) then
        do k=1,nprim
          if(cs(k).eq.zero) cs(k)=one
        end do
      end if
      inum=i
      do 180 i=1,nata
         ncf=ncf+isum*(ngc+1)
         ncs=ncs+1
c        note that ncs is the number of sets of general contractions
c        say, 3 general contractions over the same primitives
c        increase ncs only by 1
c        if (ncs.ge.maxsh) then
c          call nerror(5,'basin','too many contracted shells',ncs,maxsh)
c        end if
c **   contraction starts
         inx(1,ncs)=nsh
         noat=iat(i)
c **    atom
         inx(2,ncs)=noat
c **    shell size (1 2or s, 3 for p, 4 for l 5 for d,..)
         inx(3,ncs)=isum
c **    number of general contractions
c       note that this is zero if there is only one contraction present
         inx(4,ncs)=ngc
c **
         inx(11,ncs)=ncf-isum*(ngc+1)
         inx(12,ncs)=mxtyp
         do 170 j=1,inum
            nsh=nsh+1
c           basdat(last+1)=eta(j)
c           basdat(last+2)=cs(j)
            basdat(1,nsh)=eta(j)
            basdat(2,nsh)=cs(j)
            if(mxtyp.eq.3) then
c             basdat(last+3)=cp(j)
              basdat(3,nsh)=cp(j)
            else
              do 160 igc=1,ngc
c               basdat(last+igc+2)=cgen(j,igc)
                basdat(igc+2,nsh)=cgen(j,igc)
 160          continue
            end if
c           last=last+10
            do 165 ii=1,3
              basdat(10+ii,nsh)=xa(ii+1,noat)
  165       continue
            nbf=nbf+isum
  170    continue
  180 continue
      go to 60
  190 if (ncs.eq.0) return
      do 195 i=1,ncs-1
        inx(5,i)=inx(1,i+1)
        inx(10,i)=inx(11,i+1)
 195  continue
      inx(5,ncs)=nsh
      inx(10,ncs)=ncf
c
c *** data on the basis functions are stored in the array basdat
c *** beginning with the address ibas. 13 consecutive numbers are
c *** stored for each (uncontracted) shell.  these are the numbers..
c *** (1) orbital exponent
c *** (2) contraction coefficient (for the s function in an l shell)
c *** (3) contraction coefficient for p functions in an l shell or
c ***    contraction coeff. for the next general contraction
c *** (4)-(10): further contr. coeff. for general contractions
c *** (11)-(13): x, y and z coordinates of the center (not needed
c *** absolutely but convenient
c
c *** in the array inx 12 integers are stored for each contracted
c *** shell. the contraction goes from inx(1,k)+1 to inx(5,k)
c *** if inx is dimensioned as inx(12,ncs). inx(2,k) contains the
c *** atom corresponding to contraction k. inx(3,k) contains the
c *** shell size.. 1 for s,3 for p, 4 for l, 5 for d, 6 for d6,
c *** 7 for f, 10 for f10, 15 for g15. inx(4) contains the number of
c *** general contraction: note that in a segmented scheme it is 0
c *** inx(5) is the last primitive of the (general) contraction
c *** inx(6) to inx (nsym+6-1) (i.e. inx(6) or inx(6),inx(7) and (8)
c *** contain the symmetry pair of the contraction, i.e. another
c *** (general) contraction. This is not very useful, as it is the
c *** symmetry pairs of the contracted functions which is  needed.
c *** inx(9) is unused
c *** inx(11,k)+1 to inx(10,k) gives the indices of the
c *** contracted functions making up the contraction. note that
c *** the contracted functions for this gen. contraction go
c *** from inx(11,k)+1 to inx(10,k).
c *** inx(12,k) gives the function type for contracted shell k.
c
c
c
  240 format (1x,a4,6x,7(a4,6x))
  270 format (1x,a4)
      end
c============================================================================
      subroutine normaliz(ncs,inx,basdat)
c   this routine normalizes the basis functions
c   ncs=number of contracted shells (input)
c   inx= contraction information (input)
c   basdat= shell information (input and output)
      implicit real*8 (a-h,o-z)
      character*4 betu(8)
      dimension basdat(13,*),inx(12,*)
      data betu/'s','p','l','d','f','g','h','i'/
c                1   2   3   4   5   6   7   8
      call getrval('acc',acc)
      data zero,one,two,three,four,five /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,
     1  5.0d0/
      eight=two*four
      eigsq3=sqrt(eight/three)
      xfnorm=sqrt(1.6d0)
      gnorm=four*four/sqrt(three*five)
      hnorm=32.0d0/sqrt(105.0d0)
      xinorm=64.0d0/sqrt(945.0d0)
      sgnorm=1.0d0
      shnorm=1.0d0
      sinorm=1.0d0
      do 290 i=1,ncs
         ityp=inx(12,i)
         ifu=inx(11,i)+1
         ifu2=inx(10,i)
c        begin and end of the contraction
         ig=inx(1,i)+1
         ie=inx(5,i)
         ngc=inx(4,i)
         do 285 igg=0,ngc
           igc=igg+2
           s1=zero
           s2=zero
           do 260 ii=ig,ie
              a=basdat(1,ii)
              do 260 jj=ig,ie
                b=basdat(1,jj)
c                        s  p   l   d   d6  f  f10 g15 h21 i28
c                        g  h   i
                go to (230,240,245,250,250,251,251,252,253,254,
     1                 252,253,254), ityp
c 230    s1=s1+bl(ja+igc)*bl(ib+igc)*two*sqrt(two*(sqrt(a*b)/(a+b))**3)
  230           s1=s1+basdat(igc,ii)*basdat(igc,jj)*two*
     1           sqrt(two*(sqrt(a*b)/(a+b))**3)
                go to 260
c 240    s1=s1+bl(ja+igc)*bl(ib+igc)*four*sqrt(two*(sqrt(a*b)/(a+b))**5)
  240           s1=s1+basdat(igc,ii)*basdat(igc,jj)*four*
     1           sqrt(two*(sqrt(a*b)/(a+b))**5)
            go to 260
c 245    s1=s1+bl(ja+1)*bl(ib+1)*two*sqrt(two*(sqrt(a*b)/(a+b))**3)
c    l shell igg=0, igc=2
  245           s1=s1+basdat(igc,ii)*basdat(igc,jj)*two*
     1           sqrt(two*(sqrt(a*b)/(a+b))**3)
c        s2=s2+bl(ja+2)*bl(ib+2)*four*sqrt(two*(sqrt(a*b)/(a+b))**5)
                s2=s2+basdat(3,ii)*basdat(3,jj)*four*
     1           sqrt(two*(sqrt(a*b)/(a+b))**5)
            go to 260
c 250   s1=s1+bl(ja+igc)*bl(ib+igc)*8.0d0*sqrt(two*(sqrt(a*b)/(a+b))**7)
  250           s1=s1+basdat(igc,ii)*basdat(igc,jj)*eight*
     1           sqrt(two*(sqrt(a*b)/(a+b))**7)
            go to 260
c 251  s1=s1+bl(ja+igc)*bl(ib+igc)*16.0d0*sqrt(two*(sqrt(a*b)/(a+b))**9)
  251           s1=s1+basdat(igc,ii)*basdat(igc,jj)*16.0d0*
     1           sqrt(two*(sqrt(a*b)/(a+b))**9)
            go to 260
c 252 s1=s1+bl(ja+igc)*bl(ib+igc)*32.0d0*sqrt(two*(sqrt(a*b)/(a+b))**11)
  252           s1=s1+basdat(igc,ii)*basdat(igc,jj)*32.0d0*
     1           sqrt(two*(sqrt(a*b)/(a+b))**11)
            go to 260
  253           s1=s1+basdat(igc,ii)*basdat(igc,jj)*64.0d0*
     1           sqrt(two*(sqrt(a*b)/(a+b))**13)
            go to 260
c 254 s1=s1+bl(ja+igc)*bl(ib+igc)*128.d0*sqrt(two*(sqrt(a*b)/(a+b))**15)
  254           s1=s1+basdat(igc,ii)*basdat(igc,jj)*128.d0*
     1           sqrt(two*(sqrt(a*b)/(a+b))**15)
  260    continue
         if (s1.lt.acc) s1=one
         if (s2.lt.acc) s2=one
         s1=one/sqrt(s1)
         s2=one/sqrt(s2)
         do 280 ii=ig,ie
            basdat(igc,ii)=basdat(igc,ii)*s1
            if(ityp.eq.3) then
              basdat(3,ii)=basdat(3,ii)*s2
            end if
            i1=inx(2,i)
            if(ityp.eq.3) then
              kk=2
            else
              kk=1
            end if
c           for some reason, the p functions are extra scaled by
c          sqrt(a)*two, and the d's by 4a. The f's are scaled by
c           a**(3/2)*sqrt(1.6) and the g's by a**2*16
 270        a=basdat(1,ii)
            if(ityp.eq.2) then
              basdat(igc,ii)=basdat(igc,ii)*sqrt(a)*two
            else if (ityp.eq.3) then
              basdat(igc+1,ii)=basdat(igc+1,ii)*sqrt(a)*two
            else if (ityp.eq.4.or.ityp.eq.5) then
              basdat(igc,ii)=basdat(igc,ii)*four*a
            else if(ityp.eq.6) then
              basdat(igc,ii)=basdat(igc,ii)*xfnorm*a*sqrt(a)
            else if(ityp.eq.7)  then
              basdat(igc,ii)=basdat(igc,ii)*eigsq3*a*sqrt(a)
            else if(ityp.eq.8)  then
              basdat(igc,ii)=basdat(igc,ii)*gnorm*a**2
            else if(ityp.eq.9)  then
              basdat(igc,ii)=basdat(igc,ii)*hnorm*a**2*sqrt(a)
            else if(ityp.eq.10)  then
              basdat(igc,ii)=basdat(igc,ii)*xinorm*a**3
            else if(ityp.eq.11)  then
              basdat(igc,ii)=basdat(igc,ii)*sgnorm*a**2
            else if(ityp.eq.12)  then
              basdat(igc,ii)=basdat(igc,ii)*shnorm*a**2*sqrt(a)
            else if(ityp.eq.13)  then
              basdat(igc,ii)=basdat(igc,ii)*sinorm*a**3
            end if
c  take out this print
cc          call getival('iout',iout)
cc           write (iout,470) betu(ityp),i1,basdat(1,ii),
cc     1      (basdat(k+igc-1,ii),k=1,kk)
cc 470  format (3x,a1,3x,i5,2x,f13.7,12f12.7)
  280    continue
c       loop over general contraction
  285    continue
  290 continue
      end
c ===================================================================
c
      subroutine preorder(na,ncs,inx1,inx)
      implicit integer(a-z)
c
c  Preorders the INX array so that basis functions are ordered
c  per ATOM, in the same order as the original input geometry.
c  This should hopefully ensure a consistent Wolinski order
c  when input to the routine <reorder> (which does not give a
c  unique ordering independent of the original input order).
c    <CURRENTLY NOT USED>
c
c  ARGUMENTS
c
c  na      -  number of atoms
c  ncs     -  number of contracted shells
c  inx1    -  temporary array for reordering
c  inx     -  INX array (12,ncs)
c             on exit is reordered per atom
c
c
      Dimension inx1(12,ncs),inx(12,ncs)
c
c
c -- loop over all atoms
c
      it = 0
      Do iat=1,na
c
c -- loop over shells
c
      Do j=1,ncs
      If(inx(2,j).eq.iat) Then
c
c -- shell is on current atom
c
       it = it+1
       inx1(1,it) = inx(1,j)
       inx1(2,it) = iat
       Do i=3,12
       inx1(i,it) = inx(i,j)
       EndDo
      EndIf
      EndDo
c
      EndDo
c
c -- copy inx1 back into inx
c
      Call ICpyVEC(12*ncs,inx1,inx)
c
      end
c ===================================================================

      subroutine reorder(ncs,inx,iorder,iswap,inxo,IAN)
      implicit real*8 (a-h,o-z)
      dimension inx(12,ncs),inxo(12,ncs),iorder(*),iswap(ncs)
      dimension IAN(*)
c----------------------------------------------------------------
c   inx is the contraction info. iorder(i) gives the position of
c   the contracted basis function (not shell) with index i in the
c   original ordering in the reordered basis set new order.
c   iswap is a temporary array used to keep track of
c   the order of shells. iswap(i) is the position of the i-th
c   rearranged contracted shell in the original order
c
c   reorder the  basis set from the lowest to the highest ang.
c   momentum. Within one ang. momentum, reorder them so that
c   the contraction depth decreases. Within one ang. mom. and
c   contr. depth, reorder them so that the basis functions
c   for the same atom TYPE follow each other.
c
c -- IMPORTANT:  Changed argument           ! JB   Nov 2004
c    IAN has actual (integer) atomic numbers. This used to be datnuc
c    the stupid attempt at a C structure in which the atomic charges
c    were embedded as reals. Using atomic charges instead of atomic
c    numbers messes up if the charges change, e.g., with ECPs
c----------------------------------------------------------------
c find the lowest ang. momentum, and within that the highest contr.
c length and within that the lowest atom number in the new
c atomic ordering
cc      write(6,*) 'original basis set'
cc      do 100 ics=1,ncs
cc      write(6,*) (inx(i,ics),i=1,12)
cc 100  continue
c----------------------------------------------------------------
      if(ncs.le.1) return
      do 200 i=1,ncs
        iswap(i)=i
 200  continue
c save the original inx
      do 300 ics=1,ncs
      do 300 i=1,12
        inxo(i,ics)=inx(i,ics)
 300  continue
c----------------------------------------------------------------
c return address for the next sweep of reordering
c
 400  continue
      iexch=0
      do 500 i=1,ncs-1
        ics0=i
        ics1=i+1
        ict0=inx(12,ics0)
        ict1=inx(12,ics1)
        ina0=inx(2,ics0)
        ina1=inx(2,ics1)
        nza0=0
        if(ina0.gt.0) nza0=IAN(ina0)
        nza1=0
        if(ina1.gt.0) nza1=IAN(ina1)
c
        iclen0=inx(5,ics0)-inx(1,ics0)-1
        iclen1=inx(5,ics1)-inx(1,ics1)-1
        if (ict1.lt.ict0) then
          iexch=1
          do 410 k=1,12
             itemp=inx(k,ics0)
             inx(k,ics0)=inx(k,ics1)
             inx(k,ics1)=itemp
 410      continue
          itemp=iswap(ics0)
          iswap(ics0)=iswap(ics1)
          iswap(ics1)=itemp
        else if (ict1.eq.ict0.and.iclen0.lt.iclen1) then
          iexch=1
          do 420 k=1,12
             itemp=inx(k,ics0)
             inx(k,ics0)=inx(k,ics1)
             inx(k,ics1)=itemp
 420      continue
          itemp=iswap(ics0)
          iswap(ics0)=iswap(ics1)
          iswap(ics1)=itemp
c       else if(ict1.eq.ict0.and.iclen0.eq.iclen1.and.ina0.gt.ina1) then
       else if(ict1.eq.ict0.and.iclen0.eq.iclen1.and.nza1.gt.nza0) then
          iexch=1
          do 430 k=1,12
             itemp=inx(k,ics0)
             inx(k,ics0)=inx(k,ics1)
             inx(k,ics1)=itemp
 430      continue
          itemp=iswap(ics0)
          iswap(ics0)=iswap(ics1)
          iswap(ics1)=itemp
        end if
 500  continue
c     now we have to re-generate the beginning-ending contraction
c     arrays (inx(11,i) and inx(10,i))
      icf=0
      do 700 i=1,ncs
c  the functions ibgo to iendo will be permuted in the new order to
c  ibg to iend
        inx(11,i)=icf
c       end of contr= beginning+(shell-size)*(1+number of gen. contr.)
        inx(10,i)=icf+inx(3,i)*(1+inx(4,i))
        icf=inx(10,i)
 700  continue
      if (iexch.eq.1) go to 400
c  now create the array iorder
      do 800 ics=1,ncs
        ibeg=inx(11,ics)+1
        iend=inx(10,ics)
c  these are the indices of the functions belonging to the present
c  contraction. this contraction used to be iswap(ics) in the
c  original order. we may need to have a copy of the original inx
c  array
        icso=iswap(ics)
        ibego=inxo(11,icso)+1
        iendo=inxo(10,icso)
        idiff=(iend-ibeg)-(iendo-ibego)
        if (idiff.ne.0) then
c          write(*,*) ics,icso,ibeg,ibego,iend,iendo
          call nerror(1,'reorder','shell size wrong',ics,icso)
        end if
        ii=0
        do 850 i=ibego,iendo
          iorder(i)=ibeg+ii
          ii=ii+1
 850    continue
 800  continue
c     write(6,*) 'reordered basis set'
c     do 600 ics=1,ncs
c     write(6,*) (inx(i,ics),i=1,12)
c 600 continue
c     write(6,*) 'new addresses of old functions'
c     ncf=inx(10,ncs)
c     write(6,*) (iorder(i),i=1,ncf)
c
      end
c===============================================================
      subroutine chkbasis(na,ncen,ncs,inx,datnuc)
      implicit real*8 (a-h,o-z)
      dimension ncen(na),inx(12,ncs),datnuc(5,na)
      do i=1,na
        ncen(i)=0
      end do
      do ics=1,ncs
        iat=inx(2,ics)
        ncen(iat)=ncen(iat)+1
      end do
      do i=1,na
        if (ncen(i).le.0) go to 100
      end do
      return
  100 continue
      if(Abs(datnuc(1,i)).gt.1.0d-5) then
        call nerror(1,'chkbasis',
     1 'check basis or use DUMMY - no basis defined for atom',i,0)
      end if
      end
c================================================================
      SUBROUTINE BasTyp2(IUnit,IOut,AtSymb,Filnam,len)
      IMPLICIT INTEGER(A-Z)
C
C  reads additional basis input
C  checks for possible read from basis set library or file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to read from
C  IOut    -  output unit (should be unit 6)
C  AtSymb  -  atomic symbol basis is being input for
C  Filnam  -  name of file from which to read basis set
C  len     -  number of characters in file name
C
C
      CHARACTER*256 Filnam,basdir,Char
      CHARACTER*8 AtSymb,blank,Finish
c
      DATA blank/'        '/, Finish/'Finish  '/
C
C  Read current line in input file
C  We are looking for cards of the form
C   for     AtSymb
C   for     AtSymb     basis=<basis name>
C   for     AtSymb     file=<filename>
C  (the format for the first two fields is A4,6X,A4)
C
      Read(IUnit,'(A256)') Char
C
C  if the first three characters on the current card are
C  NOT 'for', then we are done
C
      call lowerca2(Char(1:3),3)
      IF(Char(1:3).NE.'for') THEN
       AtSymb = Finish
       BackSpace IUnit
       RETURN
      ELSE
       AtSymb = Char(11:14)       ! WARNING - assumes fixed format
       call lowerca2(AtSymb,4)
C
C  check for "basis" or "file" directive
C
       DO I=14,75
       call lowerca2(Char(I:I+4),5)
       If(Char(I:I+3).EQ.'file') Then
        Filnam = Char(I+5:256)
        Call rmblan(Filnam,251-I+1,len)
        RETURN
       Else If(Char(I:I+4).EQ.'basis') Then
C
C  standard basis
C  get basis library
C
        call getchval('BASDIR',basdir)
        call rmblan(basdir,256,len0)
        Filnam = Char(I+6:256)
        call rmblan(Filnam,250-I+1,len)
        call lowerca2(Filnam,len)
        call massagefn(Filnam,len)
        write(IOut,*)'Basis set ',Filnam(1:len),' from library'
        Filnam = basdir(1:len0)//Filnam(1:len)//'.bas'
        len = len0+len+4
        RETURN
       EndIf
       EndDO
C
C  no basis/file marker found
C
       AtSymb = blank
       len = 0
       BackSpace IUnit
      ENDIF
C
      END
c================================================================
      subroutine rawbasis(ncs,inx,nsh,basdat)
      implicit real*8 (a-h,o-z)
c  This routine write a file with the primitive basis set data
c  Arguments:
c  ncs      = number of contracted shells
c  inc      = contraction array
c  nsh      = number of primitive shells
c  basdata  = basis set data
c
c  See the definition of these arrays in SUBROUTINE BASIN
      dimension inx(12,ncs),basdat(13,nsh),istart(10)
      data istart /0,1,4,8,13,19,26,36,51,73/
      character*256 jobname
      character*6 comp(51),typex,functype
      data comp /'S     ','Px    ','Py    ','Pz    ',
     1           'S     ','Px    ','Py    ','Pz    ',
     2  'Dz2-r2','Dx2-y2','Dxy   ','Dxz   ','Dyz   ',
     3  'Dxx   ','Dyy   ','Dzz   ','Dxy   ','Dxz   ','Dyz   ',
     4  'Fx2y  ','Fx2z  ','Fxy2  ','Fy2z  ','Fxz2  ','Fyz2  ','Fxyz  ',
     6  'Fxxx  ','Fxxy  ','Fxxz  ','Fxyy  ','Fxyz  ','Fxzz  ','Fyyy  ',
     7  'Fyyz  ','Fyzz  ','Fzzz  ',
     8  'Gxxxx ','Gxxxy ','Gxxxz ','Gxxyy ','Gxxyz ','Gxxzz ','Gxyyy ',
     9  'Gxyyz ','Gxyzz ','Gxzzz ','Gyyyy ','Gyyyz ','Gyyzz ','Gyzzz ',
     &  'Gzzzz '/
      rfnorm=1.0d0/sqrt(1.6d0)
      rgnorm=sqrt(15.0d0)/16.0d0
      sq38=sqrt(0.375d0)
c  First open the raw basis file
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,len)
      open(1,file=jobname(1:len)//'.raw',form='formatted',
     1       status='unknown')
      write(1,100)
  100 format
     1 ('Type   Exponent     Coefficient      Atom         ',
     2  'x             y            z')
c  Loop over contr. shells
      do ics=1,ncs
      ityp=inx(12,ics)
c  nfun= number of Cartesian components (3 for P, 5 for D)
      ncomp=inx(3,ics)
c  Atom associated with the basis function\
      iatom=inx(2,ics)
c  The position of the first component in the component array istart
      icomp1=istart(ityp)
c  first and last primitive
      ifirst=inx(1,ics)+1
      ilast=inx(5,ics)
c  First and last genrally contracted functions (test only)
      ifirstc=inx(11,ics)
      ilastc =inx(10,ics)
c  number of general contractions=ngc
      ngc=inx(4,ics)
      do igc=0,ngc
      do icomp=1,ncomp
      functype=comp(icomp1+icomp)
      do ish=ifirst,ilast
      if(ish.gt.ifirst) functype='      '
c  basline contains the type for the first member of a contraction ONLY
      coef=basdat(2+igc,ish)
c  For L type functions, the P contr. coefficient differs from the S
      if(ityp.eq.3.and.icomp.gt.1) coef=basdat(3,ish)
      xc=basdat(11,ish)
      yc=basdat(12,ish)
      zc=basdat(13,ish)
      expo=basdat(1,ish)
c coefficients are scaled by (4*eta)**(L/2)
      if(ityp.eq.2.or.ityp.eq.3.and.icomp.gt.1) then
        coef=coef*0.5d0/sqrt(expo)
      else if (ityp.eq.4.or.ityp.eq.5.and.icomp.gt.4) then
        coef=coef*0.25d0/expo
      else if(ityp.eq.5.and.icomp.lt.4) then
        coef=coef*0.433012702d0/expo
      else if (ityp.eq.6) then
        coef=coef*rfnorm/expo/sqrt(expo)
      else if (ityp.eq.7) then
        coef=coef*sq38/expo**2
      end if
      write(1,200)  functype,expo,coef,iatom,xc,yc,zc
 200  format(a6,1x,2f14.8,i5,1x,3f14.8)
      end do  ! ish
      end do  ! Cartesian components icomp
      end do  ! general contractions igc
      end do  !  contracted shells ics
      close(1)
      end
c================================================================
      subroutine check4stab(na,ncs,inx,datbas,xnuc,istab)
      implicit real*8 (a-h,o-z)
      dimension inx(12,ncs),datbas(13,*),xnuc(5,na)
c
c     correspondance between ityp and ityp1:
c     ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28)  11   12   13
c     ityp1 1    2    3    4    4     5    5      6      7       8        6    7    8
c
c
c  check for integral stability
c  returns istab  (now modified  June 2002   JB)
c      1  -  stable
c      0  -  might be unstable, print warning message only
c     -1  -  switch on stability checking
c
c assuming , everything is stable
c
      istab=+1
c
      call get_maxtype(inx, ncs, itype_max) ! returns ityp1
      if(itype_max.le.3) return       ! only s,p,l shells
c     in Devin's toluene-xe job (dp|dp) was unstable !!!
c
      expon_max=0.d0
      expon_min=1.d+10
      do ics=1,ncs
         itype=inx(12,ics)
         ia=inx(1,ics)+1
         ie=inx(5,ics)
         do ips=ia,ie
            aa=datbas(1,ips)
            expon_max=max(expon_max,aa)
            expon_min=min(expon_min,aa)
         enddo
      enddo
c
      boamax=expon_max/max(1.d0,expon_min)
      boamax=2.d0*boamax
c
c     10**4.5=31622.8
      if(boamax.gt.31622.8d0 .and. itype_max.ge.4) istab=-1  ! d and higher
c
c     6-311++g-dp up to Ne should be stable
c     6-311++g-dp above Ne should be unstable
c     6-311++g-3df3pd up to Ne should be stable
c     6-311++g-3df3pd above Ne should be unstable
c     cc-pqvTz up to Ne should be stable
c     cc-pqvQz up to B  should be stable
c     cc-pqvQz above B  should be unstable
c
c     in Devin's toluene-xe job (dp|dp) was unstable !!!
c
c -- additional testing
c -- only heavy elements appear to be really unstable
c -- atomic numbers : 1 -  2    I row (H-He)
c                     3 - 10   II row (Li-Ne)
c                    11 - 18  III row (Na-Ar)
c                    19 - 36   IV row (K -Kr)
c                    37 - 54    V row (Rb-Xe)
c                    55 - 86   VI row (Cs-Rn)
c     correspondance between ityp and ityp1:
c     ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28)  11   12  13
c     ityp1 1    2    3    4    4     5    5      6      7       8        6    7   8
c
c find heavy elements : maximum atomic number
c
      If(istab.EQ.-1) Then
        iatno_max=0
        do i=1,na
          iatno=nint(xnuc(1,i))
          iatno_max=max(iatno_max,iatno)
        enddo
c -- iatno_max lowered from 18 to 10 - JB Now 2004
        if(iatno_max.le.10.AND.itype_max.le.5) istab=0  ! III row element & f-orbitals
        if(iatno_max.le.36.AND.itype_max.le.4) istab=0  !  IV row element & d-orbitals
      EndIf
c
      end
c================================================================
      subroutine check4lgc(ncs,inx,lshells,lgcshell,ltype_gc,ldeep_gc)
      dimension inx(12,ncs)
c
      ltype_gc=0     ! max type of gc shell
      ldeep_gc=0     ! max deep of gc shell
      lshells=0
      lgcshell=0
      do ics=1,ncs
         itype=inx(12,ics)
         if(itype.eq.3) lshells=lshells+1
         igc=inx(4,ics)   ! 0 if segmented
         if(igc.gt.0) then
            lgcshell=lgcshell+1
            ltype_gc=max(ltype_gc,itype)   ! s,p,d .....
            ldeep_gc=max(ldeep_gc,igc)     ! 0,1,2
         endif
      enddo
c
      end
c================================================================
      subroutine ChkAngMom(ncs,inx,iwvf)
      implicit integer(a-z)
      integer inx(12,ncs)
c
c  checks basis set
c  currently DFT is limited to F functions
c  correlated wavefunctions are limited to G functions
c  check this limit is not exceeded
c
c  input argument iwvf should be either
c     1 for DFT
c     2 for correlated wavefunctions
c
      IF(iwvf.EQ.1) THEN
        Do ics=1,ncs
        If(inx(12,ics).gt.7) Then
          call nerror(99,'DFT Calculation',
     $      'No higher than F functions can be used with DFT',0,0)
        EndIf
        EndDo
      ELSE
        Do ics=1,ncs
        If(inx(12,ics).gt.8.and.inx(12,ics).ne.11) Then
          call nerror(99,'Correlated Calculation',
     $      'No higher than G functions with Correlated Wavefunctions',
     $       0,0)
        EndIf
        EndDo
      ENDIF
c
      return
      end
