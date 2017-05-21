      Subroutine ump2(NAlpha,NBeta)
c
      use memory
c
c-----------------------------------------------------------------------
c
c The MP2 program for the Unrestricted Hartree-Fock wavefunction
c
c JB 2009/10 and KW 2010
c-----------------------------------------------------------------------
C
C    UMP2 MODULE
C
C    FORMULA
C
C        alpha/alpha term
C    EAA  =  Sum(i<=j) Sum(a<b) [(ia|jb)-(ib|ja)]**2/(ei+ej-ea-eb)
C
C         beta/beta term
C    EBB  =  Sum(i<=j) Sum(a<b) [(ia|jb)-(ib|ja)]**2/(ei+ej-ea-eb)
C
C        alpha/beta term       (I,A alpha; j,b beta)
C    EAB  =  Sum(I,j) Sum(A,b) (IA|jb)**2/(eI+ej-eA-eb)
C
C    E2  =  EAA  +  EBB  +  EAB
C
C
C    The algorithm is as follows (courtesy of PP):
C
C    Loop over Lam, Mu<=Lam
C      Calculate the X matrix     (i.e., the integrals)
C      Calculate  Ya = Ca(T)*X  and  Yb = Cb(T)*X
C      Calculate
C        Ya*Ca  to give the alpha-alpha half-transformed integrals
C        Ya*Cb  to give half of the alpha-beta half-transformed integrals
C        Yb*Ca  to give the other half of the alpha-beta    "
C        Yb*Cb  to give the beta-beta half-transformed integrals
C        Store all of these half-transformed integrals
C      End
C    End
C
C    Sort all of the half-transformed integrals to reorder indices
C
C    Then
C      EAA  =  Va(T)*(Ya*Ca)*Va
C      EAB  =  Va(T)*(Ya*Cb)*Vb  +  Vb(T)*(Ya*Cb)*Va  +
C              Va(T)*(Yb*Ca)*Vb  +  Vb(T)*(Yb*Ca)*Va
C      EBB  =  Vb(T)*(Yb*Cb)*Vb
C
c-----------------------------------------------------------------------
C
      implicit real*8 (a-h,o-z)
      integer*4 restrt
      character*256 jobname,scrfile,filename,
     $              filname1,filname2,filname3,filname4
      logical LastSymPair,emp2only,scs,incore0,incore
c
      dimension xintxx(9)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      common /job/jobname,lenJ
      Parameter (mdisk1=51,mdisk2=52,mdisk3=53,mdisk4=54)    ! UMP2 unit nos.
c
      parameter(nopt=16)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      data smal/.TRUE./
      character*256 chopv(nopt)
      character*3 ch3
      character opnames*4,scftype*11,wvfnc*20,cdum*20
      logical nofr,exst,dualbasis,smal,Test,alphabeta
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
ckw
      dimension Ipass(2,28), Kpass(2,28) !28 is 28 comp.of cartisian I-function
c-----------------------------------------------------------------------

      data opnames  /'nofr','orba','orbb','thre','core',
     1               'rest','dual','maxd','keep','prin',
     2               'grad','pmij','scs ','sos ','noio','test'/
c
c  0=logical option, 1=single integer, 2=two integers, 3=three integers
c  11=single real, 12=two reals, 13=three reals, 21=character option
c
      data ioptyp / 0, 2, 2,11,11, 1, 1,11, 0, 1, 0, 2,12, 0, 0, 1/
c
c  NOFR = no frozen core correlate all orbitals
c  ORBA = correlate alpha orbitals from i to j
c  ORBB = correlate beta orbitals from i to j
c  THRE = integral thresold, n pH form, default is  9, meaning 1.0d-9
c  CORE = the orbital energy separating the cores from the valence, in
c         atomic units, default is -3 Eh; if set in RHF, that value is
c         used
c  REST = restart option
c          1 - restart after alpha-alpha half-transformed integrals
c          2 - restart after beta-beta half-transformed integrals
c          3 - restart after alpha-beta half-transformed integrals
c  DUAL = dual basis MP2 (with optional read-in of dual type)
c  MAXD = maximum disk space in gigabytes (default 50 GB)
c  KEEP = do not delete the half-transformed integral files and the
c         sort file after the job
c  PRIN = print level, default=2
c  GRAD = calculate the extra quantities needed for an MP2 gradient
c         and write them to disk
c  PMIJ = print par correlation matrices (beware!!) from i to j
c  SCS  = use scaled MP2 (with optional read-in of 2 scale factors)
c  SOS  = use Head-Gordon scaled MP2, with alpha-beta term only
c  TYPE = dual-basis type (1, 2 or 3)
c  NOIO = use in-core algorithm for 2nd half-transformation
c  TEST = undocumented feature. Will stop the job after the first
c         half-transformation. Used to test the restart.
c
c
      call secund(tt0)
      call elapsec(elaps0)
c
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c  get symmetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c  basis function symmetry pairs are in inx(ifp)
      iprnt=1
c .............................................................
c -- CHECK! Only do UMP2 if current wavefunction is UHF
c --   (at the same time get the lowest eigenvalue of the overlap
c --    and the SCF energy)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,cdum)
      Call rdcntrl(IUnit,5,'$escf',2,idum,euhf,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(wvfnc(1:3).EQ.'UHF' .or. wvfnc(1:3).EQ.'UMP') then
      else
        call nerror(1,'UMP2 Module',
     $  'UMP2 Requested But Current Wavefun. is NOT UHF or UMP2 !',0,0)
      endif
c .............................................................
c
      call readopt(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
c
c  'NOFR' no frozen core
      nofr=.false.
      if(ifound(1).gt.0) nofr=.true.
c  ORBA/ORBB - see below
c  THRESHOLD (integral threshold)
      if(ifound(4).gt.0) then
        thresh=10.0d0**(-ropv(1,4))
      else
c -- check lowest eigenvalue of overlap (from SCF)
        thresh = MIN(1.0d-10,xlow**2)
        If(thresh.LT.1.0d-12) then
          write(iout,2000) xlow
 2000 Format('**WARNING** Smallest eigenvalue of overlap is ',d12.4,/,
     $       '  Final MP2 energy may have greater than mhartree error')
          thresh = 1.0d-12
        EndIf
      end if
      call setrval('thresh',thresh)
c  core is the limit separating core from valence;
c  it is normally set in scf but it can be set here too
      call tstrval('core',iyes)
      if(iyes.eq.1) then
        core=rgetrval('core')
      else
        core=-3.0d0
      end if
      if(ifound(5).gt.0)  core=ropv(1,5)
c  restart
      restrt=0
      if(ifound(6).gt.0) restrt=iopv(1,6)
      if(restrt.LT.0.OR.restrt.GT.3) call nerror(2,'UMP2 Module',
     $    'Impossible Restart Option - Check your Input',0,0)
c .............................................................................
c -- File handling --
c  filname1-filname3 = half-transformed integrals
c  filname4 = sorted bins
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filname1=scrfile(1:len)//'.htr.AA'
      filname2=scrfile(1:len)//'.htr.AB'
      filname3=scrfile(1:len)//'.htr.BB'
      len1 = len+7
      filname4=scrfile(1:len)//'.bins'
      len4 = len+5
c .........................................................................
      dualbasis=ifound(7).gt.0
c -- dual type (1, 2 or 3)
      itype=1
      If(iopv(1,7).gt.0) Then
        itype=iopv(1,7)
        if(itype.gt.3) itype=1
      EndIf
c .........................................................................
      natoms=igetival('na')
c  *************************************************
c   WARNING!  I don't think this routine properly
c    takes into account dummy atoms   ! JB Nov 07
c  *************************************************
      if(ifound(8).gt.0) then
c  Disk allotment in gigabytes
        xmaxdisk=ropv(1,8)
        xmaxdisk=125 000 000.0d0*xmaxdisk     ! convert to doublewords
      else
        xmaxdisk=6 250 000 000.0d0      ! default is 50 GB
      end if
      if(ifound(10).gt.0) iprnt=iopv(1,10)
c  'PMIJ' print matrix; first number is 1 (AO Xchange) or 2
c   (both AO & MO) second number: pair to be printed
      ipmij=0
      ipairpr=0
      if(ifound(12).gt.0) then
        ipmij=iopv(1,12)
        ipairpr=iopv(2,12)
      end if
c -- SCS scaling --
      scs=.false.
      iscs=0
      p1=1.2d0
      p2=1.0d0/3.0d0
      alphabeta = ifound(14).gt.0      ! Head-Gordon scaled MP2
      If(ifound(13).gt.0) Then
        scs=.true.
        iscs=1
        if(ropv(1,13).gt.0) then
          p11=ropv(1,13)
          p22=ropv(2,13)
          if(abs(p1-p11).gt.1.0d-05) then
            iscs=2
            p1=p11
            p2=p22
            write(iout,1234) p1,p2
 1234       format(' Scaled MP2 with p1: ',F8.5,' p2: ',F8.5)
            if(p2.eq.zero) alphabeta = .True.
          endif
        else
          write(iout,*) ' Conventional SCS-MP2'
        endif
      EndIf
c
      incore0 = ifound(15).gt.0
      if(incore0) write(iout,*) ' ** Forcing in-core Algorithm **'
c
      Test = ifound(16).gt.0
c
      call setival('iscs',iscs)
      call setrval('p1',p1)
      call setrval('p2',p2)
c-----------------------------------------------------------
c -- determine logical variable emp2only
c    normally FORCE is the next line for an MP2 gradient and
c    we are going to check for this automatically here
      READ(inp,'(A4)',END=95) chopv(12)(1:4)
      call lowercas(chopv(12),4)
      If(chopv(12)(1:4).EQ.'forc') ifound(12)=1
      BACKSPACE inp
c
 95   emp2only=.true.
      mp2only = 1
      If(ifound(12).gt.0) Then
        emp2only=.false.
        mp2only=-1
        filname3=scrfile(1:len)//'.Tij'
        len3 = len+4
        filname4=scrfile(1:len)//'.Kov'
        len4 = len+4
      EndIf
      call setival('mp2only',mp2only)
c-----------------------------------------------------------
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
c-----------------------------------------------------------
c allocate memory for an array showing if a big bs shell
c was present (0 or 1) in a small bs
c
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
c
c this is needed for dual basis set but is also used
c in blocking for ordinary MP2 (must have all 1)
c------------------------------------------------------------
  45  format(72('-'))
  46  format(/72('='))
      write(iout,46)
      write(iout,*) '                           The UMP2 Module  '
      write(iout,*)' '
c-----------------------------------------------------------
      if(dualbasis) then
c        the current "extended" basis set was already read-in
c        after SCF. Get the smaller bais set info (one used for SCF):
c
         call get_small_basis(bl,ncs,ncf,bl(ictr), ncs_sm,ncf_sm)
c
         write(iout,45)
         write(iout,*)
     *   '              Dual basis set will be used in MP2 calculations'
         write(iout,*) ' '
         write(iout,47) ncs_sm,ncf_sm
         write(iout,48) ncs   ,ncf
   47    format(7x,i5,' shells & ',i5,
     *    ' contracted functions for occupied MOs')
   48    format(7x,i5,' shells & ',i5,
     *    ' contracted functions for virtual  MOs')
         write(iout,45)
         write(iout,49) itype
         write(iout,45)
   49    format(22x,'dual basis set approach type = ',i2)
c-----------------------------------------------------------
      else
         ncs_sm=ncs
         ncf_sm=ncf
         call set_mpres(bl(mpres_in_bg),ncs)
      endif
c-----------------------------------------------------------
      write(iout,*) ' MP2 integral thresh    = ',thresh
      write(iout,*) ' core orbitals with eps < ', core
      write(iout,*) ' max disk storage (GB)  = ',8*xmaxdisk*1.d-9
      write(iout,*) ' '
      call f_lush(iout)
c-----------------------------------------------------------
c  put down a memory marker
      call matreset
c
      call mmark
      call matmark
c
c-----------------------------------------------------------
c zero out irrelevant dft stuff
      nfock=1
      idft=0
      ax=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
c-------------------------------------------------------
c allocate memory for an array mapping contr.func. to contr.shells :
c
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
c
c get memory for MOs, orbital energies and screening density
c
      np4=4
      call matdef('cmoA','q',ncf,ncf)
      call matdef('cmoB','q',ncf,ncf)
      call matdef('eorbA','d',ncf,ncf)
      call matdef('eorbB','d',ncf,ncf)
      call matdef('dsmxA','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      call matdef('dsmxB','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      icanoA=mataddr('cmoA')
      ieorbA=mataddr('eorbA')
      idsmxA=mataddr('dsmxA')
      icanoB=mataddr('cmoB')
      ieorbB=mataddr('eorbB')
      idsmxB=mataddr('dsmxB')
c-------------------------------------------------------
      If(dualbasis) Then
c
c get eigenvectors & eigenvalues
c
c the Fock operator is defined in a small basis set
c but its matrix representation goes over the big one
c It will be diagonalized WITHOUT allowing occupied
c and virtuals orbitals to mix .
c
        call get_mix_eigen2(bl,nblocks,bl(ictr),labels,thresh,
     *                      ncf,ncf_sm,ncs,ncs_sm,nmo,itype)
c
c  integral threshold can be changed on return if "big"
c  basis set overlap shows linear dependency
c
      Else
c
c -- set up the MOs filename
        call tstchval('mosfname',iyes)
        if(iyes.eq.1) then
          call getchval('mosfname',filename)
        else
          filename = jobname(1:lenJ)
        endif
        call rmblan2(filename,256,lenM)
        if(lenM.le.252) filename(lenM+1:lenM+4)='.mos'
        lenM=lenM+4
c
c -- read in the MOs
        itype = 1
        CALL ReadMOS(ncf,bl(icanoA),bl(ieorbA),.True.,lenM,
     $               filename(1:lenM),itype,IErr)
        If(IErr.NE.0) Call nerror(4,'UMP2 Module',
     $     'MOS File Does Not Exist',0,0)
        If(NBeta.GT.0) Then
          filename(lenM:lenM)='b'
          CALL ReadMOS(ncf,bl(icanoB),bl(ieorbB),.True.,lenM,
     $                 filename(1:lenM),itype,IErr)
          If(IErr.NE.0) Call nerror(4,'UMP2 Module',
     $       'MOB File Does Not Exist',0,0)
        endif
      EndIf
c-------------------------------------------------------
c initialize the two-el. integral program
c
      thint=thresh
      iforwhat=5
      call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
     *              ax,     nrad,   nang,   lrad,   lang,
     *              Iradq,  NBatch, .false.,NAlpha, NBeta,
     *              scftype,xintxx, nblocks,maxbuffer,maxlabels)
c
c-------------------------------------------------------
c
c  determine the orbitals to be correlated
c
      if(ifound(2).gt.0) then
        nfirstA=iopv(1,2)
        nlastA=iopv(2,2)
        nvalA=nlastA-nfirstA+1
        ncoreA=NAlpha-nvalA
      else
        if(nofr) then
          nfirstA=1
          nlastA=NAlpha
          nvalA=NAlpha
          ncoreA=0
        else
          ncoreA=ncores(ncf,bl(ieorbA),core)
          nfirstA=ncoreA+1
          nlastA=NAlpha
          nvalA=NAlpha-ncoreA
        endif
      endif
c
      if(ifound(3).gt.0) then
        nfirstB=iopv(1,3)
        nlastB=iopv(2,3)
        nvalB=nlastB-nfirstB+1
        ncoreB=NBeta-nvalB
      else
        if(nofr) then
          nfirstB=1
          nlastB=NBeta
          nvalB=NBeta
          ncoreB=0
        else
          ncoreB=ncores(ncf,bl(ieorbB),core)
          nfirstB=ncoreB+1
          nlastB=NBeta
          nvalB=NBeta-ncoreB
        endif
      endif
c---------------------------------------------------------------------
c Construct "screening density" matrix to be used to estimate the
c importance of integrals. It can be based either on canonical or
c localized orbitals. The latter has some advantage : integral
c prescreening based on localizability arguments is supposed to be
c more efficient .
c
c For ordinary one-basis set MP2 the "screening density" is made
c out of canonical or localized vectors taken from the disk as
c "evec_uhf" or "loca_uhf" .
c
c For dual-basis set the "screening density is formed using current
c MOS content. In the input deck, the GUESS=READ statement MUST precede
c the MP2 or LMP2 keyword in order to ensure the "corresponding orbital"
c transformation (projection from small to big basis set) has been
c preformed.
c
c For dual basis set make screening density using CANONICAL orbitals
c
      if(.not.dualbasis) then
c
c check if Localization was performed
c -- (1) alpha spin
        call sudat(np4,'loca_rhf',ni)
        if(ni.gt.0) then
          loca_done=1
      write(iout,*)' Localized orbitals are used for integral screening'
      write(iout,*)' '
          if(restrt.EQ.0) Then
            call DmxMakeL(dualbasis, nvalA, NAlpha, ncf, ncs,
     $                    np4,.false.,inx, iprnt,bl(idsmxA))
          endif
        else
          loca_done=0
          write(iout,*)
     1 ' Canonical orbitals are used for integral screeing; results',
     2 ' may be inaccurate for large molecules. Please use the LOCA',
     3 ' option in SCF (say LOCA=PIPEK) to improve accuracy and',
     4 ' performance for large systems'
          call nerror(5,'UMP2 Module',
     1  'Please use localized orbitals in SCF for screening',0,0)
          if(restrt.EQ.0)
     1    call DmxMakeC(dualbasis,nvalA,NAlpha,ncf,ncs,np4,inx,iprnt)
        endif
c
c -- (2) beta spin
        call sudat(np4,'loca_uhf',ni)
        if(ni.gt.0) then
          if(restrt.EQ.0) then
            call DmxMakeL(dualbasis, nvalB, NBeta,  ncf, ncs,
     $                    np4,.true., inx, iprnt,bl(idsmxB))
          endif
        else
          call nerror(13,'UMP2 Module',
     $      'Problems with localized beta spin MOs',0,0)
        endif
c
      else      !  dual-basis
        call DmxMakeC(dualbasis,nvalA,NAlpha,ncf,ncs,np4,inx,iprnt)
      endif
c---------------------------------------------------------------------
c make ONE screening density matrix; put it in alpha location
c
       call makeOneScreen(ncs*ncs,bl(idsmxA),bl(idsmxB))
c---------------------------------------------------------------------
c
      iorbprA=1
      iorbprB=1
c  establish orbital symmetry characteristics
      if(nsym.gt.0) then
        call getint(nvalA*nsym,iorbprA)
        call getint(nvalB*nsym,iorbprB)
        call getmem(ncf,itmp)
        call OrbSymPair(nsym,  ncf,  bl(icanoA), bl(ifp) ,nfirstA,
     1                  nlastA, iprnt, bl(itmp), bl(iorbprA),iout)
        call OrbSymPair(nsym,  ncf,  bl(icanoB), bl(ifp) ,nfirstB,
     1                  nlastB, iprnt, bl(itmp), bl(iorbprB),iout)
        call retmem(1)
      endif
c
      lastvirt=igetival('nonredun')
c
c  print out some informsation
      write(iout,'("  Number of contracted functions       =",i8,/,
     1             "  Number of alpha correlated orbitals  =",i8,/,
     2             "  Number of alpha virtual orbitals     =",i8,/,
     3             "  Number of beta correlated orbitals   =",i8,/,
     4             "  Number of beta virtual orbitals      =",i8)')
     5                ncf,nvalA,lastvirt-NAlpha,nvalB,lastvirt-NBeta
      write(iout,*) ' '
      call f_lush(iout)
c.................................................
c check if there will be a split in integral calculations
c
      call check_sizes(bl(ictr),ncs,ncf,nvalA)
c.................................................
cc      if(ncoreA.gt.0) call matsub('occaA','cmoA',1,NAlpha)
      call  matsub('occuA','cmoA',nfirstA,nlastA)
      call  matsub('virtA','cmoA',NAlpha+1,lastvirt)
      ioccA=mataddr('occuA')
      ivrtA=mataddr('virtA')
c
cc      if(ncoreB.gt.0) call matsub('occaB','cmoB',1,NBeta)
      call  matsub('occuB','cmoB',nfirstB,nlastB)
      call  matsub('virtB','cmoB',NBeta+1,lastvirt)
      ioccB=mataddr('occuB')
      ivrtB=mataddr('virtB')
c.................................................
      if(iprnt.gt.3) then
        call matprint ('occuA',6)
        call matprint ('virtA',6)
        call matprint ('occuB',6)
        call matprint ('virtB',6)
      end if
c
c***********************************************************************
c
c  FIRST HALF-TRANSFORMED INTEGRALS
c
c***********************************************************************
c
c  open direct access scratch files to store half-transformed integrals
c
      lrec1 = 5*(nvalA**2) + 16
      lrec2 = 5*(2*nvalA*nvalB) + 16
      lrec3 = 5*(nvalB**2) + 16
c
      OPEN (UNIT=mdisk2,FILE=filname2(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec2)
      If(.NOT.alphabeta) Then
        OPEN (UNIT=mdisk1,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec1)
        OPEN (UNIT=mdisk3,FILE=filname3(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec3)
      EndIf
c
      write(iout,1235) lrec1,lrec2,lrec3
 1235 format(' Record length for integral files (bytes): ',3(I10,1x))
c-----------------------------------------------------------------------
c
c  reserve space for one AO integral matrix
      call matdef('xmat','q',ncf,ncf)
      ixadr=mataddr('xmat')
c
c  put down a memory marker
      call mmark
c
c  space for the half-transformed exchange operators
      nval = MAX(nvalA,nvalB)
      call getmem(nval**2,ihalf)
c
c  reserve space for the various half-transformed integrals
c  ints     holds the half-transformed integrals in integer form
c  int1     holds a 1-byte integer
c          (mimic 5-byte integers to allow lower integral threshold)
c
      call getint_4(nvalA**2,intsAA)
      call getint_1(nvalA**2,int1AA)
      call getint_4(nvalA*nvalB,intsAB)
      call getint_1(nvalA*nvalB,int1AB)
      call getint_4(nvalB*nvalA,intsBA)
      call getint_1(nvalB*nvalA,int1BA)
      call getint_4(nvalB**2,intsBB)
      call getint_1(nvalB**2,int1BB)
c-----------------------------------------------------------------------
c
      thint=thresh
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
c  nrec is the total number of records written on all files
c  irec is the current counter in the file
      nrec=0
      irec=0
      if(restrt.EQ.1) then
        mulam=ncf*(ncf-1)/2
        mulamd=ncf
        call symmoff
        go to 50
      end if
c---------------------------------------------------
c     if(dualbasis.and.itype.eq.1) then
c        update content of "basis sets" according to current needs
c        i.e. put small basis set into "basis 2" & "basis 4"
c
c        call update_basis(bl,ncs_sm)
c     endif
c---------------------------------------------------
c
      call secund(tt)
      if(iprnt.gt.3) write(iout,*) 'Startup time=',tt-tt0
      call secund(tt3)
      call elapsec(telap3)
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c  iktot is the number of ics,kcs pairs really calculated & transformed
c   if an (ics,kcs) pair is skipped because there are no integrals in it
c   then it is NOT incremented
      iktot=0
c  Count the actual number of (mu, lam) pairs, both the diagonal
c   and the non-diagonal ones
      mulam=0
      mulamd=0
      ENonZero=zero
ckw..........
      call secund(txxx1)
      call init_info('fock')
ckw..........
c
c turn off symmetry for mp2 integrals :
      call symmoff
c
c  bl(icol) and bl(jcol) store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
      call getint(ncf,irow)
      call getint(ncf,icol)
      call getint(ncf,irow1)
      call getint(ncf,icol1)
      call getint(ncf,lzero)
c--------------------------------------------------------------------
      ncf2=ncf*ncf
c--------------------------------------------------------------------
      nskipped=0
      do ics=ncs,1,-1
      call get_shell_size(bl(ictr),ics,ics_size)
      lmp2_siz1=ncf*ncf*ics_size
      do kcs=ics,1,-1
      if(LastSymPair(ics,kcs,nsym,bl(ifp1),inegl,iret)) cycle
c  of each (ics,kcs) pair, calculate only the last one
c  (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
c
      call get_shell_size(bl(ictr),kcs,kcs_size)
c
      lmp2_size=lmp2_siz1*kcs_size
c
c check if size of the mp2 integral buffer is not too big
c if it is then split over kcs (and possibly ics)
c
      call check_size1(lmp2_size,ncf,nvalA,nvalA,ntimes)
c
      call set_passes(bl(ictr),ics,kcs,ics_size,kcs_size,
     *                ntimes,Ipass,Kpass,Itimes,Ktimes)
c
      do itime=1,itimes
      icf1=ipass(1,itime)
      icf2=ipass(2,itime)
      iatonce=icf2-icf1+1
      do ktime=1,ktimes
c
        call secund(tt2)
        call elapsec(telap2)
c
        kcf1=kpass(1,ktime)
        kcf2=kpass(2,ktime)
        katonce=kcf2-kcf1+1
c
        lmp2_size=iatonce*katonce*ncf2
c
        call getmem(lmp2_size,lmp2int)
        call mmark
        call int_lmp2(bl, bl(ictr), thresh, ics, icf1,
     1                icf2, kcs, kcf1, kcf2, bl(mapf2s),
     2                bl(idsmxA),iprnt,bl(lmp2int),nintotal,nrow,
     3                ncol, bl(irow),bl(icol),bl(lzero))
        call retmark
        if(nintotal.eq.0) then
          nskipped=nskipped+1
          call retmem(1)
          cycle
        else
          iktot=iktot+1
        endif
c
        call secund(tt3)
        call elapsec(telap3)
        tinteg=tinteg+tt3-tt2
        elapint=elapint+telap3-telap2
c......................................................................
        call TransOneShlAB(ncf,    ncs,    nvalA,  nvalB,  ics,
     1                     kcs,    icf1,   icf2,   kcf1,   kcf2,
     2             bl(ictr),bl(lmp2int),bl(ioccA),bl(ioccB),alphabeta,
     3                   iprnt,thresh,bl(ixadr),bl(ihalf), bl(ihalf),
     4             irec, bl(intsAA),bl(int1AA),bl(intsAB),bl(int1AB),
     5                   bl(intsBA),bl(int1BA),bl(intsBB),bl(int1BB),
     6                     mulam,  mulamd, nrow,   ncol, bl(irow),
     7                bl(icol),bl(irow1),bl(icol1),ENonZero, smal,
     8                     ihalf,  mdisk1, mdisk2, mdisk3)
c
        call secund(tt4)
        call elapsec(telap4)
        ttrans=ttrans+tt4-tt3
        elaptrans=elaptrans+telap4-telap3
c......................................................................
        call retmem(1)          ! lmp2int
c......................................................................
        enddo    ! over ktime (selected kcf belonging to kcs shell )
        enddo    ! over itime (selected icf belonging to ics shell )
c......................................................................
      end do     !  over kcs shell
      end do     !  over ics shell
c-------------------------------------------------------------------
      if(nskipped.gt.0) then
        write(iout,*) nskipped,
     1   ' pairs were skipped because no integrals were calculated'
      end if
c-----------------------------------------------------------------------
c timing here is irrelevant :
      call secund(txxx2)
      txxx0=txxx1
      call term_info(thresh,txxx2-txxx1,txxx1-txxx0,'fock')
c
      if(iprnt.gt.2) then
         if(nsym.gt.0) then
            write(iout,*) 'Number of shell pairs omitted by symmetry=',
     1       inegl, ' retained=',iret,'transformed',iktot
         endif
      endif
c-----------------------------------------------------------------------
      nrec=nrec+irec
      if(iprnt.ge.1) then
         write(iout,'(a,i10)')' Total number of records written=',nrec
      endif
c-----------------------------------------------------------------------
c temporary close files with half-transformed integrals
c
      close (unit=mdisk2,status='keep')
      If(.NOT.alphabeta) Then
        close (unit=mdisk1,status='keep')
        close (unit=mdisk3,status='keep')
      EndIf
c-----------------------------------------------------------------------
c
      if(iprnt.ge.1) then
c        Calculate statistics
         PerCent=100.0d0*ENonZero/(dble(mulam+mulamd)*ncf**2)
         write(iout,40) PerCent
      endif
c
  40  format(
     *' Percentage of the AO matrix used in the transformation=',f7.3/)
c......................................................................
c  release memory
  50  continue
      call retmark
      call matrem('xmat')
c......................................................................
      write(iout,61)
   61 format('  CPU & Elapsed timings in the MP2 module ')
      write(iout,*) '  '
c
      write(iout,*)' ----------- First  half transformation -----------'
c
      write(iout,62) tinteg/sixty, elapint/sixty
   62 format('  Integrals for MP2  = ',f8.2,'  and ',f8.2,' mins')
      write(iout,64) ttrans/sixty, elaptrans/sixty
   64 format('  First half transf. = ',f8.2,'  and ',f8.2,' mins')
      write(iout,*) '  '
      call f_lush(iout)
ckw---------------------------------------------------------
      if(restrt.EQ.0 .and. iprnt.ge.2) then
         call getrval('timeblks',timeblks)
         write(91,123) timeblks/sixty
  123    format('CPU time for new blocking =',f8.2)
      endif
ckw---------------------------------------------------------
      call secund(tt1)
      call elapsec(elaps1)
c
c .....................................................
      If(Test) Then
        write(iout,65)
   65   format('  == Calculation Halted at User Request ==')
        STOP
       EndIf
c .....................................................
c
c
c *******************************************************************
c
c  SECOND HALF-TRANSFORMATION
c
c *******************************************************************
c
      call elapsec(elaps0)
      call secund(tt0)
c
      lcore=igetival('lcore')
      ioffset=igetival('ioffset')
c --  lastvirt is less than ncf if redundant basis functions suppressed
      lastvirt=igetival('nonredun')
      nvirtA = lastvirt-NAlpha
      nvirtB = lastvirt-NBeta
c
c-----------------------------------------------------------------------
c  put down a memory marker
      call mmark
c
c-----------------------------------------------------------------------
c  reserve space for the various half-transformed integrals
c  ints     holds the half-transformed integrals in integer form
c  int1     holds a 1-byte integer
c          (mimic 5-byte integers to allow lower integral threshold)
c
      call getint_4(nvalA**2,intsAA)
      call getint_1(nvalA**2,int1AA)
      call getint_4(nvalA*nvalB,intsAB)
      call getint_1(nvalA*nvalB,int1AB)
c
c-----------------------------------------------------------------------
c  reserve space for batch of fully transformed integrals plus intermediate
c  storage
c
      nvirt = MAX(nvirtA,nvirtB)
      call getmem(nvirt**2,ixtrn)
      call getmem(ncf*nvirt,iz)
c
c-----------------------------------------------------------------------
c  Integrals are stored as ALL IJ values for each Mu,Lam (Mu >= Lam)
c  They need to be resorted as ALL Mu,Lam for each IJ
c  This will be done by a modified bin sort
c
c-----------------------------------------------------------------------
c  determine memory available for the bins. They should be as big
c  as possible but leave space for other allocations
c-----------------------------------------------------------------------
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      memoccu=lastaddr-ioffset
      memavail=memtot-memoccu
      membin=0.8*memavail               ! units are doublewords (8-bytes)
c-----------------------------------------------------------------------
c
      If(membin.LT.1 000 000) call nerror(14,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
c-----------------------------------------------------------------------
C  For the alpha-alpha and beta-beta contributions we need an IJ and its
C  corresponding JI pair at the same time. Any split of the total number
C  of IJ pairs MUST ensure that both IJ and JI pairs are handled together.
C  This applies in particular to the parallel algorithm where IJ pairs
C  are split between the slaves. It is accomplised by reordering all the
C  indices so that each IJ pair is followed by its corresponding JI pair
C  with the diagonal pairs (I=J) coming at the end. Thus
C    1 1   1 2   1 3   1 4   .....
C    2 1   2 2   2 3   2 4   .....
C    3 1   3 2   3 3   3 4   .....
C  becomes
C    1 2   2 1   1 3   3 1   1 4   4 1   .....
C    2 3   3 2   2 4   4 2   2 5   5 2   .....
C    3 4   4 3   3 5   5 3   .....
C  with    1 1   2 2   3 3   4 4   .....      at the end
C  The various IJ batches (if they cannot all be dealt with in the same pass)
C  should as far as possible contain an even number of IJ pairs. The last
C  batch can be odd as this contains the diagonal pairs. Similar comments
C  apply to any division among the slaves in the parallel algorithm
c-----------------------------------------------------------------------
c
c  initialize
c
      eAA = zero
      eAB = zero
      eBB = zero
      tnorm= one 
c
      tbin = zero
      tsort = zero
      ttran = zero
c
      itA = ieorbA+nfirstA-1     ! start of alpha orbital energies
      itB = ieorbB+nfirstB-1     ! start of beta  orbital energies
c
c-----------------------------------------------------------------------
c
c -- we are going to handle each type separately
c
c -- determine the disk storage used by current integral files
      halfspace1 = float(nrec)*float(lrec1)/8.0d0      ! file size in doublewords
      halfspace2 = float(nrec)*float(lrec2)/8.0d0
      halfspace3 = float(nrec)*float(lrec3)/8.0d0
c
      IF(alphabeta) GO TO 999
c
c
c  (1)  ALPHA-ALPHA
c
c----------------------------------------------------------------------
c  reopen the Alpha-Alpha half-transformed integral file
c
      OPEN (UNIT=mdisk1,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec1,STATUS='OLD')
c-----------------------------------------------------------------------
c
      call mmark
c
c -- allocate the Wolinski reordering array
      call getint(nvalA*nvalA,listAA)
      call getint(nvalA*nvalA,listAA_back)
      call makeListAA(nvalA,bl(listAA),bl(listAA_back))
c
c -- determine the disk storage currently being used
      halfspace = halfspace1 + halfspace2 + halfspace3
c
c determine how many passes will be needed to do the second half 
c transformation and the number of IJ indices per pass.
c this depends on the available disk storage
c
      npairs = nvalA**2             !  total number of IJ pairs
c
      If(incore0) Then
        incore = .True.
      Else
      CALL MemSortUMP2(membin, xmaxdisk, halfspace, halfspace1, ncf2,
     $                 npairs,   1,      1,         npass,      nIJ,
     $              nIntPerBin,nBinPerRec,  nrem,   numrec,     incore)
c
      If(incore) Then
        If(iprnt.GT.0) write(iout,*) ' ** Insufficient Disk Storage ',
     $           '-- Switching to In-Core Algorithm **'
      EndIf
      EndIf
c
      IF(.NOT.incore) THEN
cc
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half alpha-alpha transformation ---'
        write(iout,1236) npass,nIJ,nIntPerBin,nBinPerRec,numrec
 1236 format(/,' Number of passes (bin sort files written): ',I4,/,
     $         ' Maximum number of IJ indices per pass:   ',I6,/,
     $         ' Maximum number of integrals per bin:   ',I8,/,
     $         ' Maximum number of IJ indices in memory:  ',I6,/,
     $         ' Number of bins per record on file:       ',I6,/)
      EndIf
c
      If(nBinPerRec.LT.5.AND.nIJ.GT.5) call nerror(15,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
c----------------------------------------------------------------------
c  open the Alpha-Alpha bin sort file
c
      lrec4 = nIntPerBin*(5*nBinPerRec + 4) + 4      ! record length in bytes
      OPEN (UNIT=mdisk4,FILE=filname4(1:len4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec4,STATUS='NEW')
c-----------------------------------------------------------------------
c
c -- allocate the memory for JMu and JLam
      call getmem(nIntPerBin/4+1,imu)          ! integer*2
      call getmem(nIntPerBin/4+1,ilam)         ! integer*2
c
c  Loop over the number of passes
c
      lpass=0
c
c  return address for multipasses
c
 2100 CONTINUE
c
c  get the range of pairs this pass
c
      lpass = lpass+1
      istrtAA = 1 + (lpass-1)*nIJ
      iendAA = MIN(istrtAA+nIJ-1,npairs)
c
c -- allocate memory for the bins
      call getint_4(nIntPerBin*(nIJ+nrem),intsIJ)
      call getint_1(nIntPerBin*(nIJ+nrem),int1IJ)
c
      call elapsec(t1)
      Call BinSendAA(nvalA,   nrec,    mdisk1, nIntPerBin,nBinPerRec,
     $               numrec,  nrem,    istrtAA,  iendAA,   bl(intsAA),
     $            bl(ints1AA),bl(imu), bl(ilam),bl(intsIJ),bl(int1IJ),
     $            bl(listAA), mdisk4,  mrec)
      call elapsec(t2)
      tbin = tbin + t2-t1
c
c -- on exit from <BinSendAA> the bins file has been written for this pass
c -- and contains mrec records in total
c
      If(iprnt.gt.0) write(iout,1237) mrec
 1237 format(1X,' Number of records on bin sort file: ',I6)
c
c -- allocate memory for sorted integrals
      call retmem(2)      !  deallocate memory for integral bins
      call getint_4(nIntPerBin*nBinPerRec,intsIJ)
      call getint_1(nIntPerBin*nBinPerRec,int1IJ)
      call getmem(ncf2*nBinPerRec,ibAA)
c
      DO iloop=1,numrec
c
      jstrt = istrtAA + (iloop-1)*nBinPerRec
      jend = MIN(jstrt+nBinPerRec-1,iendAA)
c
      call elapsec(t1)
      Call BinSrt1AA(ncf,      nvalA,    nsym,     mdisk4, nIntPerBin,
     1             nBinPerRec, numrec,   mrec,     iloop, bl(iorbprA),
     2               bl(ifp),  thresh,   jstrt,    jend,  bl(intsIJ),
     3             bl(int1IJ), bl(imu), bl(ilam), bl(ibAA), bl(listAA),
     4                bl(listAA_back))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAA(ncf,      nvirtA,   nvalA,    jstrt,    jend,
     1                bl(ibAA), bl(ivrtA),bl(itA),  bl(iz),  bl(ixtrn),
     2                eAA,      tnorm,  bl(listAA), iprnt,    ibAA,
     3                ivrtA,    ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(3)
c
c     go back to make next pass
c
      If(iendAA.LT.npairs) GO TO 2100
c
      call retmark
c
c  -- close integral files
      close (unit=mdisk1,status='delete')
      close (unit=mdisk4,status='delete')
cc
      ELSE
c
c=====================================================================
c  IN-CORE ALGORITHM
c=====================================================================
c
      CALL MemSrt1UMP2(membin, ncf2, npairs, nIJ, npass)
c
      If(nIJ.LT.2) call nerror(15,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half alpha-alpha transformation ---'
        write(iout,1238) npairs,npass,nIJ
 1238 format(/,10X'  ** IN-CORE ALGORITHM **',/,
     $         ' Total number of IJ indices:            ',I8,/,
     $         ' Number of passes:                          ',I4,/,
     $         ' Maximum number of IJ indices per pass:   ',I6,/)
      EndIf
c
c -- allocate memory for sorted integrals
      call getmem(ncf2*nIJ,ibAA)
c
      Do lpass=1,npass
c
c  get the range of pairs this pass
c
      istrtAA = 1 + (lpass-1)*nIJ
      iendAA = MIN(istrtAA+nIJ-1,npairs)
c
      call elapsec(t1)
      Call BinSrtAA(ncf,      nvalA,    nsym,     mdisk1,   nrec,
     1           bl(iorbprA), bl(ifp),  thresh,   istrtAA,  iendAA,
     2           bl(intsAA), bl(int1AA),bl(ibAA), bl(listAA),
     3           bl(listAA_back))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAA(ncf,      nvirtA,   nvalA,    istrtAA,  iendAA,
     1                bl(ibAA), bl(ivrtA),bl(itA),  bl(iz),  bl(ixtrn),
     2                eAA,      tnorm,  bl(listAA), iprnt,    ibAA,
     3                ivrtA,    ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(1)
cc
      ENDIF
c
c
c  (2) BETA-BETA
c
c----------------------------------------------------------------------
c  reopen the Beta-Beta half-transformed integral file
c
      OPEN (UNIT=mdisk3,FILE=filname3(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec3,STATUS='OLD')
c-----------------------------------------------------------------------
c
      call mmark
c
c -- allocate the Wolinski reordering array
      call getint(nvalB*nvalB,listBB)
      call getint(nvalB*nvalB,listBB_back)
      call makeListAA(nvalB,bl(listBB),bl(listBB_back))
c
c -- determine the disk storage currently being used
      halfspace =  halfspace2 + halfspace3
c
c determine how many passes will be needed to do the second half
c transformation and the number of IJ indices per pass.
c this depends on the available disk storage
c
      npairs = nvalB**2             !  total number of IJ pairs
c
      If(incore0) Then
        incore = .True.
      Else
      CALL MemSortUMP2(membin, xmaxdisk, halfspace, halfspace1, ncf2,
     $                 npairs,   1,      1,         npass,      nIJ,
     $              nIntPerBin,nBinPerRec,  nrem,   numrec,     incore)
c
      If(incore) Then
        If(iprnt.GT.0) write(iout,*) ' ** Insufficient Disk Storage ',
     $           '-- Switching to In-Core Algorithm **'
      EndIf
      EndIf
c
      IF(.NOT.incore) THEN
cc
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half beta-beta transformation ---'
        write(iout,1236) npass,nIJ,nIntPerBin,nBinPerRec,numrec
      EndIf
c
      If(nBinPerRec.LT.5.AND.nIJ.GT.5) call nerror(15,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
c----------------------------------------------------------------------
c  open the Beta-Beta bin sort file
c
      lrec4 = nIntPerBin*(5*nBinPerRec + 4) + 4      ! record length in bytes
      OPEN (UNIT=mdisk4,FILE=filname4(1:len4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec4,STATUS='NEW')
c-----------------------------------------------------------------------
c
c -- allocate the memory for JMu and JLam
      call getmem(nIntPerBin/4+1,imu)          ! integer*2
      call getmem(nIntPerBin/4+1,ilam)         ! integer*2
c
c  Loop over the number of passes
c
      lpass=0
c
c     return address for multipasses
c
 2200 CONTINUE
c
c  get the range of pairs this pass
c
      lpass = lpass+1
      istrtBB = 1 + (lpass-1)*nIJ
      iendBB = MIN(istrtBB+nIJ-1,npairs)
c
c -- allocate memory for the bins
      call getint_4(nIntPerBin*(nIJ+nrem),intsIJ)
      call getint_1(nIntPerBin*(nIJ+nrem),int1IJ)
c
      call elapsec(t1)
      Call BinSendAA(nvalB,   nrec,    mdisk3, nIntPerBin,nBinPerRec,
     $               numrec,  nrem,    istrtBB,  iendBB,   bl(intsAA),
     $            bl(ints1AA),bl(imu), bl(ilam),bl(intsIJ),bl(int1IJ),
     $            bl(listBB), mdisk4,  mrec)
      call elapsec(t2)
      tbin = tbin + t2-t1
c
c -- on exit from <BinSendAA> the bins file has been written for this pass
c -- and contains mrec records in total
c
      If(iprnt.gt.0) write(iout,1237) mrec
c
c -- allocate memory for sorted integrals
      call retmem(2)      !  deallocate memory for integral bins
      call getint_4(nIntPerBin*nBinPerRec,intsIJ)
      call getint_1(nIntPerBin*nBinPerRec,int1IJ)
      call getmem(ncf2*nBinPerRec,ibBB)
c
      DO iloop=1,numrec
c
      jstrt = istrtBB + (iloop-1)*nBinPerRec
      jend = MIN(jstrt+nBinPerRec-1,iendBB)
c
      call elapsec(t1)
      Call BinSrt1AA(ncf,      nvalB,    nsym,     mdisk4, nIntPerBin,
     1             nBinPerRec, numrec,   mrec,     iloop, bl(iorbprB),
     2               bl(ifp),  thresh,   jstrt,    jend,  bl(intsIJ),
     3             bl(int1IJ), bl(imu), bl(ilam), bl(ibBB), bl(listBB),
     4                bl(listBB_back))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAA(ncf,      nvirtB,   nvalB,    jstrt,    jend,
     1                bl(ibBB), bl(ivrtB),bl(itB),  bl(iz),  bl(ixtrn),
     2                eBB,      tnorm,  bl(listBB), iprnt,    ibBB,
     3                ivrtB,    ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(3)
c
c     go back to make next pass
c
      If(iendBB.LT.npairs) GO TO 2200
c
      call retmark
c
c  -- close integral files
      close (unit=mdisk3,status='delete')
      close (unit=mdisk4,status='delete')
c
      ELSE
c
c=====================================================================
c  IN-CORE ALGORITHM
c=====================================================================
c
      CALL MemSrt1UMP2(membin, ncf2, npairs, nIJ, npass)
c
      If(nIJ.LT.2) call nerror(15,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half beta-beta transformation ---'
        write(iout,1238) npairs,npass,nIJ
      EndIf
c
c -- allocate memory for sorted integrals
      call getmem(ncf2*nIJ,ibBB)
c
      Do lpass=1,npass
c
c  get the range of pairs this pass
c
      istrtBB = 1 + (lpass-1)*nIJ
      iendBB = MIN(istrtBB+nIJ-1,npairs)
c
      call elapsec(t1)
      Call BinSrtAA(ncf,      nvalB,    nsym,     mdisk3,   nrec,
     1           bl(iorbprB), bl(ifp),  thresh,   istrtBB,  iendBB,
     2           bl(intsAA), bl(int1AA),bl(ibBB), bl(listBB),
     3           bl(listBB_back))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAA(ncf,      nvirtB,   nvalB,    istrtBB,  iendBB,
     1                bl(ibBB), bl(ivrtB),bl(itB),  bl(iz),  bl(ixtrn),
     2                eBB,      tnorm,  bl(listBB), iprnt,    ibBB,
     3                ivrtB,    ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(1)
cc
      ENDIF
c
  999 CONTINUE
c
c
c  (3) ALPHA-BETA
c
c----------------------------------------------------------------------
c  reopen the Alpha-Beta half-transformed integral file
c
      OPEN (UNIT=mdisk2,FILE=filname2(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec2,STATUS='OLD')
c-----------------------------------------------------------------------
c
c -- determine the disk storage used by current integral files
      halfspace =  halfspace2
c
c determine how many passes will be needed to do the second-half
c transformation and the number of IJ indices per pass.
c This depends on the available disk storage
c
      npairs = nvalA*nvalB          !  total number of IJ pairs
c
      If(incore0) Then
        incore = .True.
      Else
      CALL MemSortUMP2(membin, xmaxdisk, halfspace, halfspace1, ncf2,
     $                 npairs,   2,      1,         npass,      nIJ,
     $              nIntPerBin,nBinPerRec,  nrem,   numrec,     incore)
c
      If(incore) Then
        If(iprnt.GT.0) write(iout,*) ' ** Insufficient Disk Storage ',
     $           '-- Switching to In-Core Algorithm **'
      EndIf
      EndIf
c
      IF(.NOT.incore) THEN
cc
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half alpha-beta transformation ---'
        write(iout,1236) npass,nIJ,nIntPerBin,nBinPerRec,numrec
      EndIf
c
      If(nBinPerRec.LT.5.AND.nIJ.GT.5) call nerror(15,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
c----------------------------------------------------------------------
c  open the Alpha-Beta bin sort file
c
      lrec4 = nIntPerBin*(10*nBinPerRec + 4) + 4      ! record length in bytes
      OPEN (UNIT=mdisk4,FILE=filname4(1:len4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec4,STATUS='NEW')
c----------------------------------------------------------------------
c
c -- allocate the memory for JMu and JLam
      call getmem(nIntPerBin/4+1,imu)          ! integer*2
      call getmem(nIntPerBin/4+1,ilam)         ! integer*2
c
c  Loop over the number of passes
c
      lpass = 0
c
c     return address for multipasses
c
 2300 CONTINUE
c
c  get the range of pairs this pass
c
      lpass = lpass+1
      istrtAB = 1 + (lpass-1)*nIJ
      iendAB = MIN(istrtAB+nIJ-1,npairs)
c
c -- allocate memory for the bins
      call getint_4(nIntPerBin*(nIJ+nrem),intsIJ)
      call getint_1(nIntPerBin*(nIJ+nrem),int1IJ)
      call getint_4(nIntPerBin*(nIJ+nrem),intsJI)
      call getint_1(nIntPerBin*(nIJ+nrem),int1JI)
c
      call elapsec(t1)
      Call BinSendAB(nvalA,   nvalB,   nrec,     mdisk2, nIntPerBin,
     $            nBinPerRec, numrec,  nrem,     istrtAB,  iendAB,
     $            bl(intsAA),bl(ints1AA),bl(intsAB),bl(int1AB),bl(imu),
     $            bl(ilam),bl(intsIJ),bl(int1IJ),bl(intsJI),bl(int1JI),
     $                mdisk4,  mrec)
      call elapsec(t2)
      tbin = tbin + t2-t1
c
c -- on exit from <BinSendAB> the bins file has been written for this pass
c -- and contains mrec records in total
c
      If(iprnt.gt.0) write(iout,1237) mrec
c
c -- allocate memory for sorted integrals
      call retmem(4)      !  deallocate memory for integral bins
      call getint_4(nIntPerBin*nBinPerRec,intsIJ)
      call getint_1(nIntPerBin*nBinPerRec,int1IJ)
      call getint_4(nIntPerBin*nBinPerRec,intsJI)
      call getint_1(nIntPerBin*nBinPerRec,int1JI)
      call getmem(ncf2*nBinPerRec,ibAB)
c
      Do iloop=1,numrec
c
      jstrt = istrtAB + (iloop-1)*nBinPerRec
      jend = MIN(jstrt+nBinPerRec-1,iendAB)
c
      call elapsec(t1)
      Call BinSrt1AB(ncf,      nvalA,    nvalB,    nsym,    mdisk4,
     1          nIntPerBin, nBinPerRec,  numrec,   mrec,    iloop,
     2          bl(iorbprA),bl(iorbprB),bl(ifp), thresh,    jstrt,
     3               jend, bl(intsIJ),bl(int1IJ),bl(intsJI),bl(int1JI),
     4               bl(imu),  bl(ilam), bl(ibAB))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAB(ncf,      nvirtA,   nvirtB,   nvalA,    nvalB,
     1                jstrt,    jend,     bl(ibAB), bl(ivrtA),bl(ivrtB),
     2                bl(itA),  bl(itB),  bl(iz),   bl(ixtrn),eAB,
     3                tnorm,    iprnt,    ibAB,     ivrtA,    ivrtB,
     4                ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(5)
c
c     go back to make next pass
c
      If(iendAB.LT.npairs) GO TO 2300
c
c  -- close integral files
      close (unit=mdisk2,status='delete')
      close (unit=mdisk4,status='delete')
cc
      ELSE
c
c=====================================================================
c  IN-CORE ALGORITHM
c=====================================================================
c
      CALL MemSrt1UMP2(membin, ncf2, npairs, nIJ, npass)
c
      If(nIJ.LT.2) call nerror(15,'UMP2 Module',
     $   'Insufficient Memory for Second Half-Transformation',0,membin)
c
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half alpha-beta transformation ---'
        write(iout,1238) npairs,npass,nIJ
      EndIf
c
c -- allocate memory for sorted integrals
      call getmem(ncf2*nIJ,ibAB)
c
      Do lpass=1,npass
c
c  get the range of pairs this pass
c
      istrtAB = 1 + (lpass-1)*nIJ
      iendAB = MIN(istrtAB+nIJ-1,npairs)
c
      call elapsec(t1)
      Call BinSrtAB(ncf,      nvalA,    nvalB,    nsym,    mdisk2,
     1              nrec, bl(iorbprA),bl(iorbprB),bl(ifp), thresh,
     2              istrtAB,  iendAB,  bl(intsAA), bl(int1AA),
     3              bl(intsAB), bl(int1AB), bl(ibAB))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAB(ncf,      nvirtA,   nvirtB,   nvalA,    nvalB,
     1                istrtAB,  iendAB,   bl(ibAB), bl(ivrtA),bl(ivrtB),
     2                bl(itA),  bl(itB),  bl(iz),   bl(ixtrn),eAB,
     3                tnorm,    iprnt,    ibAB,     ivrtA,    ivrtB,
     4                ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(1)
cc
      ENDIF
c
c-----------------------------------------------------------------------
c
c -- that's it!
c
      If(iprnt.gt.0) WRITE(iout,2500) tbin/sixty,tsort/sixty,ttran/sixty
 2500 Format(/,' Time spent in bin sort:         ',F10.2,' mins',/,
     $         ' Time spent unpacking integrals: ',F10.2,' mins',/,
     $         ' Time in integral transformation:',F10.2,' mins',/)
c
      call secund(tt1)
      call elapsec(elaps1)
c
      WRITE(iout,705) (tt1-tt0)/sixty,(elaps1-elaps0)/sixty
 705  Format(1X,' Second half-transformation = ',f8.2,' and ',f8.2,
     $          ' mins')
c
c *******************************************************
c
c Print out the final results :
c
      emp2 = eAA + eBB + eAB
      tnorm=sqrt(tnorm)
c
      write(iout,*) ' '
      write(iout,*) '-------------------------------------------------'
      write(iout,*) '                Final SCF/MP2 results '
      write(iout,*) ' '
      write(iout,3001) euhf
      write(iout,*) '-------------------------------------------------'
      write(iout,3002) eAA,eAB,eBB
      write(iout,*) '-------------------------------------------------'
      write(iout,3003) emp2
      write(iout,*) '-------------------------------------------------'
      write(iout,3004) euhf+emp2
      write(iout,3005) tnorm
      If(alphabeta) then
      write(iout,*) '-------------------------------------------------'
      write(iout,*) '              Scaled Opposite Spin MP2    '
      write(iout,*) '   [See Head-Gordon, J.Chem.Phys. 121 (2004) 9793]'
      write(iout,*) ' '
        emp2 = 1.3d0*eAB
        write(iout,3039) emp2
        write(iout,3040) euhf+emp2
      else
      write(iout,*) '-------------------------------------------------'
      write(iout,*) '                    Scaled MP2    '
      write(iout,*) '      [See S.Grimme, J.Chem.Phys. 118 (2003) 9095]'
      write(iout,*) ' '
        emp2scs = 1.2d0*eAB + (eAA+eBB)/3.0d0
        write(iout,3035) emp2scs
        write(iout,3036) euhf+emp2scs
        if(iscs.eq.2) then
          emp2scs = p1*eAB + p2*(eAA+eBB)
          write(iout,3037)emp2scs,p1,p2
          write(iout,3038) euhf+emp2scs
        endif
      endif
      write(iout,*) '-------------------------------------------------'
      write(iout,*) ' '
c
 3001 format('        Total SCF energy         =',f15.8)
 3002 format('        alpha-alpha contrib.     =',f15.8,/,
     $       '        alpha-beta  contrib.     =',f15.8,/,
     $       '        beta -beta  contrib.     =',f15.8)
 3003 format('        MP2 correlation energy   =',f15.8)
 3004 format('        Total MP2 energy         =',f15.8)
 3005 format('        Norm of the wavefunction =',f15.8)
c
 3035 format('        SCS-MP2 energy           =',f15.8)
 3036 format('        Total SCS-MP2 energy     =',f15.8)
 3037 format('        Scaled MP2 energy        =',f15.8,' p1=',f5.3,
     $        'p2=',f5.3)
 3038 format('        Total Scaled MP2 energy  =',f15.8)
 3039 format('        SOS-MP2 energy           =',f15.8)
 3040 format('        Total SOS-MP2 energy     =',f15.8)
c
      call f_lush(6)
c -- set wavefunction type
      wvfnc = 'UMP2'
      If(dualbasis) wvfnc = 'UMP2-dual'
      If(scs) wvfnc = 'UMP2-scs'
      If(scs.AND.dualbasis) wvfnc = 'UMP2-dual-scs'
      If(alphabeta) wvfnc = 'UMP2-sos'
c
c -- write final MP2 energy to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      If(iscs.eq.0) Then
        Call wrcntrl(IUnit,7,'$energy',2,idum,euhf+emp2,wvfnc)
      Else
        Call wrcntrl(IUnit,7,'$energy',2,idum,euhf+emp2scs,wvfnc)
      Endif
      If(dualbasis) Call wrcntrl(IUnit,6,'$edual',2,idum,edual,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
      call symmon        ! turn symmetry back on for integrals
c---------------------------------------------------------------------
      call retmark
c---------------------------------------------------------------------
      call setival('ncoreA',ncoreA)
      call setival('ncoreB',ncoreB)
      call setrval('thrsx',thresh)
c---------------------------------------------------------------------
c
c -- release allocated memory
      call matremark
      call retmark
      call memory_status('end of mp2')
c---------------------------------------------------------------------
      end
c=======================================================================
c
      subroutine TransOneShlAB(ncf,    ncs,    nvalA,  nvalB,  ics,
     1                         kcs,    icf1,   icf2,   kcf1,   kcf2,
     2                         inx,    xint,   CMOA,   CMOB, alphabeta,
     3                         iprnt,  thrsh,  xmat,   halfA,  halfB,
     4                         nrec,   intsAA, int1AA, instAB, int1AB,
     5                         intsBA, int1BA, instBB, int1BB, mulam,
     6                         mulamd, nrow,   ncol,   irow,   icol,
     7                         irow1,  icol1, ENonZero,smal,   ihalf,
     8                         mdisk1, mdisk2, mdisk3)

      use memory

      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c
c  transform a batch of integrals
c
c  ARGUMENTS
c  ncf     number of contracted basis functions
c  ncs     number of contracted shells
c  nvalA   number of alpha MOs to be transformed (usually the valence ones)
c  nvalB   number of beta  MOs to be transformed
c  ics,kcs: contracted shells for which the transformation is carried out
c  inx:    array containing contraction info
c  xint:   AO integrals ordered as
c          (l=1:ncf,j=1:ncf,k=kcf1:kcf2,i=icf1:icf2)   where icf1 & icf2
c          are the first and the last contr. funct. in the ICS shell ,
c          similarly for kcs .
c          NOTE THAT THEY ARE SCALED by 1/thresh
c  CMOA    alpha MO coefficients
c  CMOB    beta  MO coefficients
c  iprnt:  print level
c  thresh  integral threshold
c  nrow    number of non-zero rows of the AO integral matrices
c  ncol      ditto for columns
c  irow(k),k=1,nrow = the indices of the non-zero rows
c  icol(k),k=1,ncol = ditto for columns
c  xmat  (ncf,ncf) storage for 1/4 transformed integrals
c  halfA   storage for half-transformed alpha integrals
c  halfB   storage for half-transformed beta integrals
c          ** NOTE: halfA and halfB can share storage **
c  nrec    number of records written  (input/output)
c  intsAA  4-byte storage for alpha-alpha half-transformed integrals
c  int1AA    ditto  1-byte overflow storage
c  intsAB  4-byte storage for alpha-beta  half-transformed integrals
c  int1AB    ditto  1-byte overflow storage
c  intsBA  4-byte storage for beta-alpha  half-transformed integrals
c  int1BA    ditto  1-byte overflow storage
c  intsBB  4-byte storage for  beta-beta  half-transformed integrals
c  int1BB    ditto  1-byte overflow storage
c  mulam: number of non-diagonal (mu, lambda) pairs transformed
c  integrals are {mu,nu|lam,isig);  mu.ne.lam
c  mulamd: number of diagonal pairs (mu,mu) transformed
c  ENonZero = number of the elements in non-zero columns and rows
c  ihalf- address of halfA and halfB matrices
c  mdisk1   unit number for alpha-alpha integrals
c  mdisk2   unit number for alpha-beta  integrals
c  mdisk3   unit number for  beta-beta  integrals
c-----------------------------------------------------------------------
c
      dimension irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
      dimension xint(*),xmat(ncf,ncf),halfA(nvalA,*),halfB(nvalB,*)
      dimension CMOA(ncf,*),CMOB(ncf,*)
      INTEGER*4 intsAA(nvalA,nvalA),intsAB(nvalA,nvalB),
     $          intsBA(nvalB,nvalA),intsBB(nvalB,nvalB)
      INTEGER*1 int1AA(nvalA,nvalA),int1AB(nvalA,nvalB),
     $          int1BA(nvalB,nvalA),int1BB(nvalB,nvalB)
      dimension inx(12,*)
      logical smal,alphabeta
c
c  return if there are no integrals
      if(nrow.eq.0.or.ncol.eq.0) RETURN
c
c  icol and jcol store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
c  itrunc serves as temporary storage for a compacted coefficient matrix
c  only the rows (or columns in the second quarter transformation) of
c  X which are non-zero are present in bl(itrunc)
c
      nval = MAX(nvalA,nvalB)
      call getmem(ncf*nval,itrunc)
c
c     kcf1=inx(11,kcs)+1
c     kcf2=inx(10,kcs)
c     icf1=inx(11,ics)+1
c     icf2=inx(10,ics)
c
      do icf=icf1,icf2
      do kcf=kcf1,kcf2
      if(kcf.gt.icf) exit
      if(smal) then
        call ExtractOne(ncf,   icf,   kcf,   icf1,  icf2,
     2                  kcf1,  kcf2,  xint,  xmat,  nrow,
     3                  ncol,  irow,  icol,  nrow1, ncol1,
     4                  irow1, icol1)
      else
        call mmark
        call Extract1KI(ncf,   icf,   kcf,  icf1,  icf2,
     2                  kcf1,  kcf2,  xint, xmat,  nrow,
     3                  ncol,  irow,  icol, nrow1, ncol1,
     4                  irow1, icol1)
        call retmark
      endif
c  Determine some statistics:
      ENonZero=ENonZero+dble(nrow*ncol)
c
      if(nrow1.eq.0.or.ncol1.eq.0) cycle
c
      call getmem(nrow1*ncol1,iymat)
      call getmem(nvalA*ncol1,iymA)
      call getmem(nvalB*ncol1,iymB)
c
      call IntTranAB(ncf,    nvalA,  nvalB,  icf,    kcf,
     1               iprnt,          nrow1,  ncol1,  irow1,
     2               icol1,  thrsh,  xmat,   CMOA,   CMOB, 
     3               intsAA, int1AA, intsAB, int1AB, bl(iymA),
     4               intsBA, int1BA, intsBB, int1BB, bl(iymB),
     5            bl(iymat),bl(itrunc),bl(itrunc),halfA,halfB,
     6                nrec, alphabeta, mdisk1, mdisk2, mdisk3,
     7               iymA,iymB,iymat,itrunc,ihalf)
c
      call retmem(3)
      if(icf.ne.kcf) then
        mulam=mulam+1
      else
        mulamd=mulamd+1
      end if
c
      end do
      end do
      call retmem(1)
      return
      end
c=====================================================================
c
      subroutine IntTranAB(ncf,    nvalA,  nvalB,  mu,      lam,
     1                     iprnt,          nrow,   ncol,    irow,
     2                     icol,   thrsh,  xmat,   coefA,   coefB,
     3                     intsAA, int1AA, intsAB, int1AB,  ymatA,
     4                     intsBA, int1BA, intsBB, int1BB,  ymatB,
     5                     ymat,   Trunc1, Trunc2, halfA,   halfB,
     6                     nrec, alphabeta, mdisk1, mdisk2, mdisk3,
     7                     iymA,iymB,iymat,itrunc,ihalf)
c---------------------------------------------------------------------
c
c  This routine transforms the integrals for a given mu, lam to MO basis,
c  yielding (mu,i|lam,j) where mu,lam are AO indices and i,j are MO
c  ones. The resulting integrals are written to disk as integers
c
c  ARGUMENTS:
c
c  ncf     number of contracted basis functions
c  nvalA   number of alpha MOs to be transformed (valence orbitals)
c  nvalB   number of beta  MOs to be transformed
c  mu,lam: fixed basis function indices. The transformed integrals are
c          (mu,i|lam,j) where i,j are MO indices (i-alpha; j-beta)
c  iprnt   print level
c  nrow, ncol: the number of nonzero rows and columns of xmat
c  irow, icol = the indices of nonzero rows + columns
c  thrsh   integral threshold
c  xmat    xmat(nu.isig)=(mu.nu | lam,isig) originally. However,
c          it is compacted here: xmat(p,q)=(mu,irow(p)|lam,icol(q))
c  coefA   alpha MO coefficients
c  coefB   beta  MO coefficients
c  intsAA  4-byte storage for alpha-alpha half-transformed integrals
c  int1AA    ditto  1-byte overflow storage
c  intsAB  4-byte storage for alpha-beta  half-transformed integrals
c  int1AB    ditto  1-byte overflow storage
c  ymatA   storage for coefA(T)*xmat
c  intsBA  4-byte storage for beta-alpha  half-transformed integrals
c  int1BA    ditto  1-byte overflow storage
c  intsBB  4-byte storage for  beta-beta  half-transformed integrals
c  int1BB    ditto  1-byte overflow storage
c  ymatB    storage for coefB(T)*xmat
c  ymat     storage for non-zero (compact) block of xmat
c  Trunc1   storage for truncated MO coefficients (left side)
c  Trunc2   storage for truncated MO coefficients (right side)
c           ** NOTE:  Trunc1 & Trunc2 can share storage
c  halfA    storage for half-transformed alpha integrals
c  halfB    storage for half-transformed beta integrals
c           ** NOTE: halfA and halfB can share storage **
c  nrec     number of records written to ndisk1,ndisk2,ndisk3
c  alphabeta  logical flag for alpha-beta contribution only
c  mdisk1   unit number for alpha-alpha integrals
c  mdisk2   unit number for alpha-beta  integrals
c  mdisk3   unit number for  beta-beta  integrals
c
c in order to use the matrix system (matconn):
c
c  iymA  -  address of the ymatA matrix
c  iymB  -  address of the ymatB matrix
c  iymat -  address of the ymat  matrix
c  itrunc-  address of the trunc1 and trunc2 matrices 
c  ihalf -  address of the halfA and halfB matrices
c---------------------------------------------------------------------
c
      use memory
c
      implicit real*8 (a-h,o-z)
      INTEGER*4 intsAA(nvalA,nvalA),intsAB(nvalA,nvalB),
     $          intsBA(nvalB,nvalA),intsBB(nvalB,nvalB)
      INTEGER*1 int1AA(nvalA,nvalA),int1AB(nvalA,nvalB),
     $          int1BA(nvalB,nvalA),int1BB(nvalB,nvalB)
      dimension halfA(nvalA,*),halfB(nvalB,*),irow(nrow),icol(ncol)
      dimension xmat(ncf,ncf),coefA(ncf,nvalA),coefB(ncf,nvalB)
      dimension ymat(nrow,ncol),ymatA(ncol,nvalA),ymatB(ncol,nvalB),
     $          Trunc1(nrow,*),Trunc2(ncol,*)
      Logical alphabeta
c
      common/marks/m(7)
c
c -- initialize the marks for the 5-byte integral storage
cc      mm = 0
cc      m(1)=ibset(mm,5)
cc      m(2)=ibset(mm,6)
cc      m(3)=ibset(m(2),5)
cc      m(4)=ibset(mm,7)
cc      m(5)=ibset(m(4),5)
cc      m(6)=ibset(m(4),6)
cc      m(7)=ibset(m(6),5)
c
c---------------------------------------------------------------------
c     write(6,*) ' for mu,lam=',mu,lam,' input xmat matrix is:'
c     do i=1,nrow
c     write(6,1111) (xmat(i,j)*thrsh,j=1,ncol)
c1111 format(1x,6f12.6)
c     enddo
c---------------------------------------------------------------------
      iout = 6     ! do NOT change this   JB
c---------------------------------------------------------------------
c
c  xmat has the COMPACT AO exchange matrix
c  perform first quarter transformation
c    ymatA = Ca(T)*xmat;   ymatB = Cb(T)*xmat
c
c  copy compacted xmat into ymat
c
      call PutTrunc(ncf,nrow,ncol,xmat,ymat)
c
c  prepare the compact alpha MO coefficients for the left side
c
      call CompactCoef(ncf,nvalA,coefA,nrow,irow,Trunc1)
c
c---------------------------------------------------------
c     DO 10 I=1,nvalA
c     DO 10 J=1,ncol
c     Temp = 0.0d0
c     DO 5 K=1,nrow
c     Temp = Temp + Trunc1(K,I)*ymat(K,J)
c 5   CONTINUE
c     ymatA(J,I) = Temp
c10   CONTINUE
c---------------------------------------------------------
c Above stuff replaced by below
c    ymatA=ymat(T)*Trunc1
c
      call matconn('ymat' ,'r',nrow,ncol,iymat)
      call matconn('trunc','r',nrow,nvalA,itrunc)
      call matconn('ymatA','r',ncol,nvalA,iymA)
c
      call matmmul2('ymat','trunc','ymatA','t','n','n')
c
      call matdisc('ymatA')
      call matdisc('trunc')
      call matdisc('ymat')
c---------------------------------------------------------
c
c  prepare the compact beta MO coefficients for the left side
c
      call CompactCoef(ncf,nvalB,coefB,nrow,irow,Trunc1)
c
c---------------------------------------------------------
c     DO 20 I=1,nvalB
c     DO 20 J=1,ncol
c     Temp = 0.0d0
c     DO 15 K=1,nrow
c     Temp = Temp + Trunc1(K,I)*ymat(K,J)
c15   CONTINUE
c     ymatB(J,I) = Temp
c20   CONTINUE
c---------------------------------------------------------
c Above stuff replaced by below
c    ymatB=ymat(T)*Trunc1
c
      call matconn('ymat' ,'r',nrow,ncol,iymat)
      call matconn('trunc','r',nrow,nvalB,itrunc)
      call matconn('ymatB','r',ncol,nvalB,iymB)
c
      call matmmul2('ymat','trunc','ymatB','t','n','n')
c
      call matdisc('ymatB')
      call matdisc('trunc')
      call matdisc('ymat')
c--------------------------------------------------
c  now complete first half transformation by forming
c    ymatA*Ca    ymatA*Cb    ymatB*Ca    ymatB*Cb
c
c  prepare the compact alpha MO coefficients for the right side
c
      call CompactCoef(ncf,nvalA,coefA,ncol,icol,Trunc2)
c
c---------------------------------------------------------
c
c     DO 30 I=1,nvalA
c     DO 30 J=1,nvalA
c     Temp = 0.0d0
c     DO 25 K=1,ncol
c     Temp = Temp + ymatA(K,I)*Trunc2(K,J)
c25   CONTINUE
c     halfA(I,J) = Temp
c30   CONTINUE
c---------------------------------------------------------
c Above stuff replaced by below
c    halfA=ymatA(T)*Trunc2
c
      IF(.NOT.alphabeta) THEN
        call matconn('ymatA','r',ncol,nvalA,iymA)
        call matconn('trunc','r',ncol,nvalA,itrunc)
        call matconn('halfA','r',nvalA,nvalA,ihalf)
c
        call matmmul2('ymatA','trunc','halfA','t','n','n')
c
        call matdisc('halfA')
        call matdisc('trunc')
        call matdisc('ymatA')
c---------------------------------------------------------
c
        CALL IntTRAN(nvalA,nvalA,halfA,intsAA,int1AA)
c
c--------------------------------------------------
c print out halfA
        if(iprnt.gt.3) then
          write(6,*) ' for mu,lam=',mu,lam,' AA-half-trans matrix is:'
          call dscal(nvalA*nvalA,thrsh,halfA,1)
          call prntmat(nvalA,nvalA,nvalA,halfA)
          call dscal(nvalA*nvalA,1.0d0/thrsh,halfA,1)
        endif
      ENDIF
c--------------------------------------------------
c
c     DO 40 I=1,nvalA
c     DO 40 J=1,nvalB
c     Temp = 0.0d0
c     DO 35 K=1,ncol
c     Temp = Temp + ymatB(K,J)*Trunc2(K,I)
c35   CONTINUE
c     halfB(J,I) = Temp
c40   CONTINUE
c---------------------------------------------------------
c Above stuff replaced by below
c    halfB=ymatB(T)*Trunc2
c
      call matconn('ymatB','r',ncol,nvalB,iymB)
      call matconn('trunc','r',ncol,nvalA,itrunc)
      call matconn('halfB','r',nvalB,nvalA,ihalf)
c
      call matmmul2('ymatB','trunc','halfB','t','n','n')
c
      call matdisc('halfB')
      call matdisc('trunc')
      call matdisc('ymatB')
c---------------------------------------------------------
c
      CALL IntTRAN(nvalB,nvalA,halfB,intsBA,int1BA)
c---------------------------------------------------------
c print out halfA
      if(iprnt.gt.3) then
        write(6,*) ' for mu,lam=',mu,lam,' BA-half-trans matrix is:'
        call dscal(nvalA*nvalB,thrsh,halfB,1)
        call prntmat(nvalB,nvalA,nvalB,halfB)
        call dscal(nvalA*nvalB,1.0d0/thrsh,halfB,1)
      endif
c---------------------------------------------------------
c
c  prepare the compact beta MO coefficients for the right side
c
      call CompactCoef(ncf,nvalB,coefB,ncol,icol,Trunc2)
c
c---------------------------------------------------------
c     DO 50 I=1,nvalB
c     DO 50 J=1,nvalA
c     Temp = 0.0d0
c     DO 45 K=1,ncol
c     Temp = Temp + ymatA(K,J)*Trunc2(K,I)
c45   CONTINUE
c     halfA(J,I) = Temp
c50   CONTINUE
c---------------------------------------------------------
c Above stuff replaced by below
c    halfA=ymatA(T)*Trunc2
c
      call matconn('ymatA','r',ncol,nvalA,iymA)
      call matconn('trunc','r',ncol,nvalB,itrunc)
      call matconn('halfA','r',nvalA,nvalB,ihalf)
c
      call matmmul2('ymatA','trunc','halfA','t','n','n')
c
      call matdisc('halfA')
      call matdisc('trunc')
      call matdisc('ymatA')
c---------------------------------------------------------
      CALL IntTRAN(nvalA,nvalB,halfA,intsAB,int1AB)
c---------------------------------------------------------
c print out halfA
      if(iprnt.gt.3) then
        write(6,*) ' for mu,lam=',mu,lam,' AB-half-trans matrix is:'
        call dscal(nvalA*nvalB,thrsh,halfA,1)
        call prntmat(nvalA,nvalB,nvalA,halfA)
        call dscal(nvalA*nvalB,1.0d0/thrsh,halfA,1)
      endif
c---------------------------------------------------------
c
c     DO 60 I=1,nvalB
c     DO 60 J=1,nvalB
c     Temp = 0.0d0
c     DO 55 K=1,ncol
c     Temp = Temp + ymatB(K,I)*Trunc2(K,J)
c55   CONTINUE
c     halfB(I,J) = Temp
c60   CONTINUE
c---------------------------------------------------------
c Above stuff replaced by below
c    halfB=ymatB(T)*Trunc2
c
      IF(.NOT.alphabeta) THEN
        call matconn('ymatB','r',ncol,nvalB,iymB)
        call matconn('trunc','r',ncol,nvalB,itrunc)
        call matconn('halfB','r',nvalB,nvalB,ihalf)
c
        call matmmul2('ymatB','trunc','halfB','t','n','n')
c
        call matdisc('halfB')
        call matdisc('trunc')
        call matdisc('ymatB')
c---------------------------------------------------------
        CALL IntTRAN(nvalB,nvalB,halfB,intsBB,int1BB)
c---------------------------------------------------------
c print out halfB
        if(iprnt.gt.3) then
          write(6,*) ' for mu,lam=',mu,lam,' BB-half-trans matrix is:'
          call dscal(nvalB*nvalB,thrsh,halfB,1)
          call prntmat(nvalB,nvalB,nvalB,halfB)
          call dscal(nvalB*nvalB,1.0d0/thrsh,halfB,1)
        endif
      ENDIF
c---------------------------------------------------------
c
      nrec = nrec+1
c
c------------------------------------------------------------------
c
      write(unit=mdisk2,rec=nrec) mu,lam,int1AB,intsAB,int1BA,intsBA
      If(.NOT.alphabeta) Then
        write(unit=mdisk1,rec=nrec) mu,lam,int1AA,intsAA
        write(unit=mdisk3,rec=nrec) mu,lam,int1BB,intsBB
      EndIf
c------------------------------------------------------------------
c
      return
      end
c=====================================================================
c
      SUBROUTINE IntTRAN(nval1,nval2,halftra,intint,int1)
      implicit real*8(a-h,o-z)
C
C  Converts half-transformed integrals from real*8 to integer*5
C  prior to writing to disk
C
      REAL*8 halftra(nval1,nval2)
      INTEGER*4 intint(nval1,nval2)
      INTEGER*1 int1(nval1,nval2),i1
      common /marks/ m(7)
c
c  the number below is the largest number to fit in a 4 byte signed integer
      parameter (zero=0.0d0,one=1.0d0,dblmax=2 147 483 648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
c
      do j=1,nval2
      do i=1,nval1
        x=halftra(i,j)
cc        call integerst(m,x,intint(i,j),int1(i,j))
        IF(abs(x).ge.dblmax) THEN
          b = x*dblinv
          if(abs(b).ge.d1max) then
            dfac = abs(x)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
            dfac = LOG10(dfac)
            i1 = -NINT(dfac+0.5d0)
            x = x*10.0d0**i1
            int1(i,j) = i1
            intint(i,j) = x
cc            write(iout,100) mu,i,lam,j,thrsh*10**-i1
cc  100   format(' half-transformed integral (mu i|lam j) = ',
cc     1      '(',i4,2x,i4,'|',i4,2x,i4,')',/,
cc     2      ' too large.  Threshold reduced to ',e20.6)
cc             call nerror(1,'IntTransp',
cc     1      'transformed integral is too large', 0,0)
          else
            i1 = abs(b)
            b = x - SIGN(i1*dblmax,x)
            int1(i,j)  = i1
            intint(i,j) = b
          endif
        ELSE
          int1(i,j) = 0
          intint(i,j) = x
        ENDIF
      end do
      end do
c
      return
      end
c=====================================================================
c
      subroutine BinSrtAA(ncf,    nvalA,  nsym,   ndisk1,   nrec,
     1                    iorbprA,ifp,    thresh, istrtAA,  iendAA,
     2                    intsAA, int1AA, xmatAA, listAA,  listAA_back) 
      implicit real*8(a-h,o-z)
C
C
C  This routine sorts the Alpha-Alpha or Beta-Beta half-transformed integrals
C
C  ARGUMENTS
C
C  ncf       -  number of AO basis functions
C  nvalA     -  number of MOs to be transformed (usually the valence ones)
C  nsym      -  number of Abelian symmetry operations (less identity)
C  ndisk1    -  unit number of the half-transformed integral file
C  nrec      -  number of records written on the half-transformed integral file
C  iorbprA   -  symmetry characteristics of each correlated MO
C               under each Abelian symmetry operation
C  ifp       -    ditto but for the basis functions
C  thresh    -  threshold for integral neglect
C  istrtAA   -  starting IJ index for current batch integrals
C  iendAA    -    ditto ending index
C  intsAA    -  4-byte storage for half-transformed integrals
C  int1AA    -    ditto  1-byte overflow storage
C  xmatAA    -  storage for reordered integrals
C  listAA    -  Wolinski reordering arrays
C  listAA_back
c-----------------------------------------------------------------------
c
      INTEGER*4 intsAA(nvalA,nvalA)
      INTEGER*1 int1AA(nvalA,nvalA)
      REAL*8 xmatAA(ncf,ncf,istrtAA:iendAA)
      DIMENSION iorbprA(nsym,nvalA),ifp(7,ncf),
     $          ipair(7),jpair(7)
c
c-----------------------------------------------------------------------
      dimension listAA(*)        ! nvalA*nvalA
      dimension listAA_back(*)   ! nvalA*nvalA
c-----------------------------------------------------------------------
      Parameter (dblmax=2 147 483 648.0d0)
c-----------------------------------------------------------------------
      dblcmp = dblmax*thresh
c-----------------------------------------------------------------------
c  initialize
c
      call ZeroIT(xmatAA,(iendAA-istrtAA+1)*ncf**2)
c
c  read all records on the half-transformed integral file and select the range
c  of IJ indices being processed this pass
c
      DO 20 irec=1,nrec
      read(unit=ndisk1,rec=irec) mu,lam,int1AA,intsAA
c
      DO 10 lp0=istrtAA,iendAA
      IJ=listAA(lp0)
      call get_ij_full(IJ,nvalA,I,J)
c--------------------------------------
        ji=(j-1)*nvalA+i
        lp1=listAA_back(ji)
c--------------------------------------
        do isym=1,nsym
        ipair(isym)=iorbprA(isym,I)
        jpair(isym)=iorbprA(isym,J)
        end do
c
c -- decompress the integral ----------
cc        call reconstruct(intsAA(I,J),int1AA(I,J),xx)
cc        xx = xx*thresh
        If(int1AA(I,J).eq.0) Then
          xx = intsAA(I,J)*thresh
        Else If(int1AA(I,J).gt.0) Then
          xx = intsAA(I,J)*thresh
          xx = xx + SIGN(int1AA(I,J)*dblcmp,xx)
        Else
          xx = intsAA(I,J)*thresh*10.0d0**(-int1AA(I,J))
        EndIf
c ------------------------------------------------------
        xmatAA(mu,lam,lp0)=xx
        xmatAA(lam,mu,lp1)=xx
c
        do isym=1,nsym
          xxx=xx
          mu1=ifp(isym,mu)
          mu2=abs(mu1)
          if(mu1.lt.0) xxx = -xxx
          lam1=ifp(isym,lam)
          lam2=abs(lam1)
          if(mu2.eq.mu.and.lam2.eq.lam) cycle
          if(lam1.lt.0) xxx = -xxx
          if(ipair(isym).lt.0) xxx = -xxx
          if(jpair(isym).lt.0) xxx = -xxx
          xmatAA(mu2,lam2,lp0)=xxx
          xmatAA(lam2,mu2,lp1)=xxx
        end do
c
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c=====================================================================
c
      Subroutine BinSrtAB(ncf,      nvalA,    nvalB,    nsym,    ndisk2,
     1                    nrec,     iorbprA,  iorbprB,  ifp,     thresh,
     2                    istrtAB,  iendAB,   intsAB,   int1AB,  intsBA,
     3                    int1BA,   xmatAB)
      implicit real*8(a-h,o-z)
C
C
C  This routine sorts the Alpha-Beta half-transformed integrals
C
C  ARGUMENTS
C  ncf       -  number of AO basis functions
C  nvalA     -  number of alpha MOs to be transformed (usually the valence ones)
C  nvalB     -  number of beta  MOs to be transformed
C  nsym      -  number of Abelian symmetry operations (less identity)
C  ndisk2    -  unit number of half transformed integral file
C  nrec      -  number of records written on the half-transformed integral file
C  iorbprA   -  symmetry characteristics of each correlated alpha MO
C               under each Abelian symmetry operation
C  iorbprB   -    ditto for the beta MOs
C  ifp       -    ditto but for the basis functions
C  thresh    -  threshold for integral neglect
C  istrtAB   -  starting IJ index for current batch alpha-beta integrals
C  iendAB    -    ditto ending index
C  intsAB    -  4-byte storage for alpha-beta  half-transformed integrals
C  int1AB    -    ditto  1-byte overflow storage
C  intsBA    -  4-byte storage for beta-alpha  half-transformed integrals
C  int1BA    -    ditto  1-byte overflow storage
C  xmatAB    -  storage for reordered alpha-beta  integrals
c-----------------------------------------------------------------------
c
      INTEGER*4 intsAB(nvalA,nvalB),intsBA(nvalB,nvalA)
      INTEGER*1 int1AB(nvalA,nvalB),int1BA(nvalB,nvalA)
      REAL*8 xmatAB(ncf,ncf,istrtAB:iendAB)
      DIMENSION iorbprA(nsym,nvalA),iorbprB(nsym,nvalB),ifp(7,ncf),
     $          ipair(7),jpair(7)
c
c-----------------------------------------------------------------------
      Parameter (dblmax=2 147 483 648.0d0)
c-----------------------------------------------------------------------
      dblcmp = dblmax*thresh
c-----------------------------------------------------------------------
c  initialize
c
      call ZeroIT(xmatAB,(iendAB-istrtAB+1)*ncf**2)
c
c  read all records on the half-transformed integral file and select the range
c  of IJ indices being processed this pass
c
      DO 20 irec=1,nrec
      read(unit=ndisk2,rec=irec) mu,lam,int1AB,intsAB,int1BA,intsBA
c
c -- alpha-beta
      DO 10 IJ=istrtAB,iendAB
      call get_ij_full(IJ,nvalB,I,J)
c
      do isym=1,nsym
      ipair(isym)=iorbprA(isym,I)
      jpair(isym)=iorbprB(isym,J)
      end do
c
c -- decompress the integral ----------
cc        call reconstruct(intsAB(I,J),int1AB(I,J),xx)
cc        xx = xx*thresh
        If(int1AB(I,J).eq.0) Then
          xx = intsAB(I,J)*thresh
        Else If(int1AB(I,J).gt.0) Then
          xx = intsAB(I,J)*thresh
          xx = xx + SIGN(int1AB(I,J)*dblcmp,xx)
        Else
          xx = intsAB(I,J)*thresh*10.0d0**(-int1AB(I,J))
        EndIf
c ------------------------------------------------------
        xmatAB(mu,lam,IJ)=xx
        do isym=1,nsym
          xxx=xx
          mu1=ifp(isym,mu)
          mu2=abs(mu1)
          if(mu1.lt.0) xxx = -xxx
          lam1=ifp(isym,lam)
          lam2=abs(lam1)
          if(mu2.eq.mu.and.lam2.eq.lam) cycle
          if(lam1.lt.0) xxx = -xxx
          if(ipair(isym).lt.0) xxx = -xxx
          if(jpair(isym).lt.0) xxx = -xxx
          xmatAB(mu2,lam2,IJ)=xxx
        end do
c
 10   CONTINUE
c
c -- beta-alpha
c
      DO 11 IJ=istrtAB,iendAB
      call get_ij_full(IJ,nvalB,I,J)
c
        do isym=1,nsym
        ipair(isym)=iorbprB(isym,J)
        jpair(isym)=iorbprA(isym,I)
        end do
c
c -- decompress the integral ----------
cc        call reconstruct(intsBA(J,I),int1BA(J,I),xx)
cc        xx = xx*thresh
        If(int1BA(J,I).eq.0) Then
          xx = intsBA(J,I)*thresh
        Else If(int1BA(J,I).gt.0) Then
          xx = intsBA(J,I)*thresh
          xx = xx + SIGN(int1BA(J,I)*dblcmp,xx)
        Else
          xx = intsBA(J,I)*thresh*10.0d0**(-int1BA(J,I))
        EndIf
c ------------------------------------------------------
        xmatAB(lam,mu,ij)=xx
        do isym=1,nsym
          xxx=xx
          mu1=ifp(isym,mu)
          mu2=abs(mu1)
          if(mu1.lt.0) xxx = -xxx
          lam1=ifp(isym,lam)
          lam2=abs(lam1)
          if(mu2.eq.mu.and.lam2.eq.lam) cycle
          if(lam1.lt.0) xxx = -xxx
          if(ipair(isym).lt.0) xxx = -xxx
          if(jpair(isym).lt.0) xxx = -xxx
          xmatAB(lam2,mu2,IJ)=xxx
        end do
c
 11   CONTINUE
 20   CONTINUE
cccccccccc
cc      write(6,*) ' About to leave <BinSrtAB>  xmatAB matrix is:'
cc      do ij=istrtAB,iendAB
cc      write(6,*) ' IJ is:',ij
cc      call prntmat(ncf,ncf,ncf,xmatAB(1,1,ij))
cc      enddo
cccccccccc
C
      RETURN
      END
c=====================================================================
c
      Subroutine TrnsVirtAA(ncf,    nvirtA,  nvalA,   istrtAA, iendAA,
     1                      xmatAA, virtA,   eorbA,   Z,       XMOS,
     2                      eAA,    tnorm,   listAA,  iprnt,   ibAA,
     3                      ivrtA,  ixtrn,   iz)
c-----------------------------------------------------------------------
C
C  Transforms the half-transformed integrals for the current range of
C  MO indices IJ to virtual orbitals
C
C  ARGUMENTS
C
C  ncf     -  number of basis functions (AOs)
C  nvirtA  -  number of alpha virtual MOs
C  nvalA   -  number of alpha occupied MOs
C  istrtAA -  starting IJ index for current batch alpha-alpha integrals
C  iendAA  -    ditto ending index
C  xmatAA  -  batch of reordered alpha-alpha integrals
C  virtA   -  alpha virtual MOs
C  eorbA   -  alpha orbital energies
C  Z       -  working storage:  ncf*nvirtA
C  XMOS    -  intermediate storage for fully transformed integrals
C
C  on exit
C
C  eAA     -  cumulative alpha-alpha pair energies
C  tnorm   -  cumulative wave function norm (squared)
C
c
c matrix addresses to use the matrix system
c
c  ibAA  -  xmatAA
c  ivrtA -  virtA
c  ixtrn -  XMOS
c  iz    -  scratch Z
c-----------------------------------------------------------------------
c
      use memory
c
      implicit real*8 (a-h,o-z)
c
      REAL*8 xmatAA(ncf,ncf,istrtAA:iendAA)
c
      Dimension virtA(ncf,nvirtA),eorbA(ncf)
      Dimension XMOS(nvirtA,*),Z(ncf,*)
c-----------------------------------------------------------------------
      dimension listAA(*)        ! nvalA*nvalA
c-----------------------------------------------------------------------
c
c -- OK, let's go!
c
      DO lp=istrtAA,iendAA
c
      IJ=listAA(lp)
      call get_ij_full(IJ,nvalA,I,J)
c
c------------------------------------------------------------
      IF(I.LT.J) THEN
c------------------------------------------------------------
        if(iprnt.gt.4) then
          write(6,*) ' IJ:',ij,' I:',i,' J:',j
          write(6,*) ' Input xmatAA matrix is:'
          call prntmat(ncf,ncf,ncf,xmatAA(1,1,lp))
          call f_lush(6)
        endif
c------------------------------------------------------------
c -- matrix multiply:  XMOS = virtA(T)*xmatAA*virtA
c -- BEWARE!  using matrix system
c
        call matconn('virtualA','r',ncf,nvirtA, ivrtA)
        lp1=lp-istrtAA
        call matconn('xmtAA','q',ncf,ncf, ibAA+lp1*ncf*ncf)
        call matconn('xmoS', 'q',nvirtA,nvirtA, ixtrn)
c
        call matsimtr('xmtAA','virtualA','xmoS')
c
        call matdisc('xmoS')
        call matdisc('xmtAA')
        call matdisc('virtualA')
c------------------------------------------------------------
        if(iprnt.gt.4) then
          write(6,*) ' Fully Transformed alpha-alpha matrix is:'
          call prntmat(nvirtA,nvirtA,nvirtA,XMOS)
          call f_lush(6)
        endif
c------------------------------------------------------------
c
        call PairEnAA(nvalA, nvirtA, I, J, XMOS, eorbA, xnorm, eIJ)
c
c------------------------------------------------------------
        if(iprnt.gt.2) then
          write(6,1234) ij,I,J,eIJ
          call f_lush(6)
        endif
 1234 format(' alpha-alpha: IJ=',i5,' i,j=',2(i5,1x),'eij= ',F12.6)
c------------------------------------------------------------
c
        eAA = eAA+eIJ
        tnorm= tnorm + xnorm
c
c------------------------------------------------------------
      ENDIF     !   IF(I.LE.J) THEN
c------------------------------------------------------------
      EndDO
C
      RETURN
      END
c=====================================================================
c
      Subroutine TrnsVirtAB(ncf,     nvirtA,  nvirtB,  nvalA,   nvalB,
     1                      istrtAB, iendAB,  xmatAB,  virtA,   virtB,
     2                      eorbA,   eorbB,   Z,       XMOS,    eAB,
     3                      tnorm,   iprnt,   ibAB,    ivrtA,   ivrtB,
     4                      ixtrn,   iz)
c-----------------------------------------------------------------------
C
C  Transforms the half-transformed integrals for the current range of
C  MO indices IJ to virtual orbitals
C
C  ARGUMENTS
C
C  ncf     -  number of basis functions (AOs)
C  nvirtA  -  number of alpha virtual MOs
C  nvirtB  -  number of beta  virtual MOs
C  nvalA   -  number of alpha occupied MOs
C  nvalB   -  number of beta  occupied MOs
C  istrtAB -  starting IJ index for current batch alpha-beta integrals
C  iendAB  -    ditto ending index
C  xmatAB  -  batch of reordered alpha-beta  integrals
C  virtA   -  alpha virtual MOs
C  virtB   -  beta  virtual MOs
C  eorbA   -  alpha orbital energies
C  eorbB   -  beta  orbital energies
C  Z       -  working storage:  ncf*nvirtA
C  XMOS    -  intermediate storage for fully transformed integrals
C
C  on exit
C
C  eAB     -  cumulative alpha-beta  pair energies
C  tnorm   -  cumaulative wave function norm (squerd)
C
c
c matrix addresses to use the matrix system
c
c  ibAB         -  xmatAB
c  ivrtA,ivrtB  -  virtA and virtB
c  ixtrn        -  XMOS
c  iz           -  scratch Z
c-----------------------------------------------------------------------
c
      use memory
c
      implicit real*8 (a-h,o-z)
c
      REAL*8 xmatAB(ncf,ncf,istrtAB:iendAB)
c
      Dimension virtA(ncf,nvirtA),virtB(ncf,nvirtB)
      dimension eorbA(ncf), eorbB(ncf)
      dimension XMOS(nvirtA,*),Z(ncf,*)
c-----------------------------------------------------------------------
c
c -- OK, let's go!
c
      DO IJ=istrtAB,iendAB
      call get_ij_full(IJ,nvalB,I,J)
c
c------------------------------------------------------------
      if(iprnt.gt.4) then
        write(6,*) ' IJ:',ij,' I:',i,' J:',j
        write(6,*) ' Input xmatAB matrix is:'
        call prntmat(ncf,ncf,ncf,xmatAB(1,1,IJ))
        call f_lush(6)
      endif
c------------------------------------------------------------
c -- matrix multiply:  XMOS = virtA(T)*xmatAB*virtB
c -- BEWARE!  using matrix system
c
      call matconn('zzzzz',   'r',nvirtA,ncf, iz)
      call matconn('virtualA','r',ncf,nvirtA, ivrtA)
      ij1=ij-istrtAB
      call matconn('xmtAB',   'q',ncf,ncf, ibAB+ij1*ncf*ncf)
c
      call matmmul2('virtualA','xmtAB','zzzzz','t','n','n')
c
      call matdisc('xmtAB')
      call matdisc('virtualA')
c
      call matconn('virtualB','r',ncf,nvirtB, ivrtB)
      call matconn('xmoS','r',nvirtA,nvirtB, ixtrn)
c
      call matmmult('zzzzz','virtualB','xmoS')
c
      call matdisc('xmoS')
      call matdisc('virtualB')
c
      call matdisc('zzzzz')
c------------------------------------------------------------
      if(iprnt.gt.4) then
        write(6,*) ' Fully Transformed alpha-Beta  matrix is:'
        call prntmat(nvirtA,nvirtB,nvirtA,XMOS)
        call f_lush(6)
      endif
c------------------------------------------------------------
c
      call PairEnAB(nvalA,  nvalB,
     *             nvirtA, nvirtB,
     *               I,      J,      XMOS,
     *             eorbA,  eorbB,  xnorm, eIJ)
c
c------------------------------------------------------------
      if(iprnt.gt.2) then
        write(6,1235) ij,I,J,0.5d0*eIJ
        call f_lush(6)
      endif
 1235 format(' alpha- beta: IJ=',i5,' i,j=',2(i5,1x),'eij= ',F12.6)
c------------------------------------------------------------
c
      eAB = eAB+eIJ
      tnorm= tnorm+xnorm
c
      EndDO
C
      RETURN
      END
c=====================================================================
c
      Subroutine PairEnAA(NAlpha,nvirtA,I,J,X,eorbA, tnorm, eIJ)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
C  Compute the pair energy for the alpha-beta term for the UMP2 energy
C
C  ARGUMENTS
C
C  NAlpha  -  number of occupied alpha MOs
C  nvirtA  -  number of correlated alpha virtual MOs
C  i,j     -  occupied MOs defining the ij pair
C  X       -  fully transformed integrals in matrix form
C             X(A,B) = (AI}BJ)
C  eorbA   -  alpha orbital energies
c  eij     -  on exit pair energy
c  tnorm   -  on exit partial wavefunction norm
c---------------------------------------------------------------------
c
      Dimension X(nvirtA,nvirtA),eorbA(*)
      parameter(zero=0.0d0,one=1.0d0)
c
      ei=eorbA(i)
      ej=eorbA(j)
      eocc=ei+ej
      eij=zero
      tnorm=zero
c
      do ia=1,nvirtA
      ea=eocc-eorbA(NAlpha+ia)
      do ib=ia+1,nvirtA
        eb=ea-eorbA(NAlpha+ib)
        rdenom=one/eb
        xmo=x(ia,ib)-x(ib,ia)
        tij=rdenom*xmo
        tnorm=tnorm+tij*tij
        eij=eij+tij*xmo
      end do
      end do
c
      end
c=====================================================================
c
      Subroutine PairEnAB(NAlpha, NBeta,  nvirtA, nvirtB, 
     *                    I,        J,      X,  
     *                    eorbA,  eorbB,  tnorm , eIJ)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
C  Compute the pair energy for the alpha-beta term for the UMP2 energy
C
C  ARGUMENTS
C
C  NAlpha  -  number of occupied alpha MOs
C  NBeta   -  number of occupied beta  MOs
C  nvirtA  -  number of correlated alpha virtual MOs
C  nvirtB  -  number of correlated beta  virtual MOs
C  i,j     -  occupied MOs defining the ij pair
C  X       -  fully transformed integrals in matrix form
C             X(A,B) = (AI}BJ)
C  eorbA   -  alpha orbital energies
C  eorbB   -  beta  orbital energies
c  eij     -  on exit pair energy
c  tnorm   -  on exit partial wavefunction norm
c---------------------------------------------------------------------
c
      Dimension X(nvirtA,nvirtB),eorbA(*),eorbB(*)
      parameter(zero=0.0d0,one=1.0d0)
c---------------------------------------------------------------------
c
      ei=eorbA(i)
      ej=eorbB(j)
      eocc=ei+ej
      eij=zero
      tnorm=zero
c
      do ib=1,nvirtB
      eb=eocc-eorbB(NBeta+ib)
      do ia=1,nvirtA
        rdenom=one/(eb-eorbA(NAlpha+ia))
        xmo=x(ia,ib)
        tij=rdenom*xmo
        tnorm=tnorm+tij*tij
        eij=eij+tij*xmo
      end do
      end do
c
      end
c=====================================================================
c
      subroutine makeOneScreen(ncs2,dsmxA,dsmxB)
      implicit real*8 (a-h,o-z)
      dimension dsmxA(ncs2), dsmxB(ncs2)
c
      do i=1,ncs2
         dsmxA(i)=max( dsmxA(i),dsmxB(i))
      enddo
c
      end
c=====================================================================
c
      subroutine makeListAA(nvalA,listAA,listAA_back)
      implicit real*8 (a-h,o-z)
      dimension listAA(nvalA*nvalA)
      dimension listAA_back(nvalA*nvalA)
C
C  This is a modification of Krys Wolinski's original routine
C
C  When doing the alpha-alpha or beta-beta IJ pairs, IJ and JI must be done
C  together, i.e., handled in the same batch.  This is because a given integral
C  with MO indices I and J contributes to both the IJ and the JI metrices.
C  This affects the way the IJ indices are divided up not just in the serial
C  algorithm, but also during the parallel algorithm when IJ indices are sent
C  to different slaves.
C
C  This is achieved by making sure that the IJ pairs are ordered such that JI
C  immediately follows IJ and further ensuring that however the IJ pairs are
C  divided up, each "batch" contains an even number of indices. (The indices
C  where I=J are put at the end, so it's OK for the last batch to have an odd
C  number of IJ indices.)
C
C
c -- off-diagonal IJ pairs
      call IZeroIT(listAA_back,nvalA**2)
      IT = 0
      Do 11 I=1,nvalA
      Do 10 J=1,nvalA
      If(J.EQ.I) GO TO 10
      IJ = (I-1)*nvalA + J
      If(listAA_back(IJ).EQ.0) Then
        IT = IT+1
        listAA(IT) = IJ
        listAA_back(IJ) = 1
      EndIf
      JI = (J-1)*nvalA + I
      If(listAA_back(JI).EQ.0) Then
        IT = IT+1
        listAA(IT) = JI
        listAA_back(JI) = 1
      EndIf
 10   CONTINUE
 11   CONTINUE
c
c -- diagonal IJ pairs (I=J)
      DO 12 I=1,nvalA
      IT = IT+1
      IJ = (I-1)*nvalA + I
      listAA(IT) = IJ
 12   CONTINUE
c
      do l=1,nvalA*nvalA
         ij=listAA(l)
         listAA_back(ij)=l
      enddo
cccccccc
cc      write(6,*) ' In <makeListAA>'
cc      write(6,*) ' listAA arrays are:'
cc      l=0
cc      do i=1,nvalA
cc      do j=1,nvalA
cc      l=l+1
cc      write(6,1234) l,i,j,listAA(l),listAA_back(l)
cc 1234 format(I6,' i:',I4,' j:',I4,' listAA:',i6,' listAA-back:',i6)
cc      enddo
cc      enddo
cccccccc
c
      return
      end
c=====================================================================
c
      SUBROUTINE BinSendAA(nval,   nrec,   ndisk,nIntPerBin,nBinPerRec,
     $                     numrec, nrem,   istrt,   iend,    intsAA,
     $                     int1AA, JMu,    JLam,    intsIJ,  int1IJ,
     $                     listAA, mdisk4, lrec)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine sorts the half-transformed integrals from mu,lam all IJ to
C  IJ all mu,lam and writes them to the bin sort file
C
C
C  ARGUMENTS
C
C  nval        -  number of occupied MOs to be correlated
C  nrec        -  number of records on the half-transformed integral file
C  ndisk       -  unit number of half-transformed integral file
C  nIntPerBin  -  number of 5-byte integrals that can be stored in each IJ bin
C                 (i.e., the bin size)
C  nBinPerRec  -  number of IJ bins written to the bin sort file per record
C  numrec      -  number of records written to bin sort file at one time
C  nrem        -  number of "padded" IJ bins in last record
C  istrt       -  starting IJ index on this pass
C  iend        -  ending IJ index on this pass
C  intsAA      -  4-byte storage for half-transformed integrals
C  int1AA      -    ditto  1-byte overflow storage
C  JMu         -  2-byte storage for mu AO-indices
C  JLam        -  2-byte storage for lamda AO-indices
C  intsIJ      -  4-byte storage for reordered integrals
C  int1IJ      -    ditto  1-byte overflow storage
C  listAA      -  Wolinski reordering array
C  mdisk4      -  unit number of bin sort file
C  lrec        -  on exit number of records on bin sort file
C  
C
      INTEGER*4 intsAA(nval,nval),intsIJ(nIntPerBin,istrt:iend+nrem)
      INTEGER*1 int1AA(nval,nval),int1IJ(nIntPerBin,istrt:iend+nrem)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
c-----------------------------------------------------------------------
      dimension listAA(*)        ! nvalA*nvalA
c-----------------------------------------------------------------------
C
C  Read all records on the half-transformed integral file and select the range of
C  IJ indices being processed this pass
C
      IT = 0
      lrec = 0      !  number of records on bin sort file
c
      DO 20 irec=1,nrec
      READ (UNIT=ndisk,rec=irec) mu,lam,int1AA,intsAA
c
      IT = IT+1
      JMu(IT) = mu
      JLam(IT) = lam
c
      DO 10 lp=istrt,iend
      IJ = listAA(lp)
      call get_ij_full(IJ,nval,I,J)
      intsIJ(IT,lp) = intsAA(I,J)
      int1IJ(IT,lp) = int1AA(I,J)
 10   CONTINUE
c
      If(IT.EQ.nIntPerBin) Then
c
c -- write out the full bins
        NonZero = nIntPerBin
        IJ = 1
        Do K=1,numrec
        CALL BinWrtAA(mdisk4, nIntPerBin, nBinPerRec, JMu,    JLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),  NonZero,  lrec)
        IJ = IJ+NBinPerRec
        EndDO
c
        IT = 0
      EndIf
 20   CONTINUE
C
C  All integrals processed
C  write out the last bins
C
      If(IT.GT.0) Then
c
        NonZero = IT
        IJ = 1
        Do K=1,numrec
        CALL BinWrtAA(mdisk4, nIntPerBin, nBinPerRec, JMu,    JLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),  NonZero,  lrec)
        IJ = IJ+NBinPerRec
        EndDO
      EndIf
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE BinWrtAA(ndisk,  nIntPerBin, nBinPerRec, JMu,    JLam,
     $                    intsIJ, int1IJ,     NonZero,   lrec)
      IMPLICIT INTEGER(A-Z)
C
C  Writes a record to the bin sort file
C
      INTEGER*4 intsIJ(nIntPerBin,nBinPerRec)
      INTEGER*1 int1IJ(nIntPerBin,nBinPerRec)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
c
      lrec = lrec+1
      WRITE(unit=ndisk,rec=lrec) JMu,JLam,int1IJ,intsIJ,NonZero
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE BinRdAA(ndisk,  nIntPerBin, nBinPerRec, JMu,    JLam,
     $                   intsIJ, int1IJ,     NonZero,   lrec)
      IMPLICIT INTEGER(A-Z)
C
C  Reads a record from the bin sort file
C
      INTEGER*4 intsIJ(nIntPerBin,nBinPerRec)
      INTEGER*1 int1IJ(nIntPerBin,nBinPerRec)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
c
      READ(unit=ndisk,rec=lrec) JMu,JLam,int1IJ,intsIJ,NonZero
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE BinSrt1AA(ncf,    nval,    nsym,   ndisk, nIntPerBin,
     $                 nBinPerRec, numrec,  mrec,   iloop,  iorbpr,
     $                     ifp,    thresh,  istrt,  iend,   intsIJ,
     $                     int1IJ, JMu,     JLam,   xmatAA, listAA,
     $                   listAA_back)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads back the alpha-alpha (or beta-beta) sorted half-transformed integrals
C  and puts then into the xmatAA matrix as reals
C
C  ARGUMENTS
C
C  ncf         -  number of basis functions
C  nval        -  number of occupied MOs to be correlated
C  nsym        -  number of Abelian symmetry operations (less identity)
C  ndisk       -  unit number of bin sort file
C  nIntPerBin  -  number of 5-byte integrals that can be stored in each IJ bin
C                 (i.e., the bin size)
C  nBinPerRec  -  number of IJ bins written to the bin sort file per record
C  numrec      -  number of records written to bin sort file at one time
C  mrec        -  total number of records on the bin sort file
C  iloop       -  pointer to the record we are currently dealing with
C  iorbpr      -  symmetry characteristics of each correlated alpha/beta MO
C                 under each Abelian symmetry operation
C  ifp         -    ditto but for the basis functions
C  thresh      -  threshold for integral neglect
C  istrt       -  starting IJ index on this pass
C  iend        -  ending IJ index on this pass
C  intsIJ      -  4-byte storage for reordered integrals
C  int1IJ      -    ditto  1-byte overflow storage
C  JMu         -  2-byte storage for mu AO-indices
C  JLam        -  2-byte storage for lamda AO-indices
C  xmatAA      -  reordered alpha-alpha integrals ready for transformation
C  listAA      -  Wolinski reordering arrays
C  listAA_back
C
C
      INTEGER*4 intsIJ(nIntPerBin,istrt:iend) 
      INTEGER*1 int1IJ(nIntPerBin,istrt:iend)
      INTEGER*2 JMu(NIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
      REAL*8 xmatAA(ncf,ncf,istrt:iend)
      DIMENSION iorbpr(nsym,nval),ifp(7,ncf),ipair(7),jpair(7)
c
c-----------------------------------------------------------------------
      dimension listAA(*)        ! nvalA*nvalA
      dimension listAA_back(*)   ! nvalA*nvalA
c-----------------------------------------------------------------------
      Parameter (dblmax=2 147 483 648.0d0)
c-----------------------------------------------------------------------
      dblcmp = dblmax*thresh
c-----------------------------------------------------------------------
C
C  initialize
C
      numIJ = iend-istrt+1
      Call ZeroIT(xmatAA,numIJ*ncf**2)
c
      lrec = iloop      !  starting record number for this batch of bins
 100  CONTINUE
c
      Call BinRdAA(ndisk,  nIntPerBin, nBinPerRec, JMu,    JLam,
     $             intsIJ(1,istrt), int1IJ(1,istrt), NonZero, lrec)
c
      DO 11 lp0=istrt,iend
c
      IJ = listAA(lp0)
      call get_ij_full(IJ,nval,I,J)
c--------------------------------------
      ji=(j-1)*nval+i
      lp1=listAA_back(ji)
c--------------------------------------
      do isym=1,nsym
      ipair(isym)=iorbpr(isym,I)
      jpair(isym)=iorbpr(isym,J)
      end do
c
      DO 10 kk=1,NonZero
      mu = JMu(kk)
      lam = JLam(kk)
c
c -- decompress the integral ----------
cc      call reconstruct(intsIJ(kk,lp0),int1IJ(kk,lp0),xx)
cc      xx = xx*thresh
        If(int1IJ(kk,lp0).eq.0) Then
          xx = intsIJ(kk,lp0)*thresh
        Else If(int1IJ(kk,lp0).gt.0) Then
          xx = intsIJ(kk,lp0)*thresh
          xx = xx + SIGN(int1IJ(kk,lp0)*dblcmp,xx)
        Else
          xx = intsIJ(kk,lp0)*thresh*10.0d0**(-int1IJ(kk,lp0))
        EndIf
c ------------------------------------------------------
        xmatAA(mu,lam,lp0)=xx
        xmatAA(lam,mu,lp1)=xx
c
        do isym=1,nsym
          xxx=xx
          mu1=ifp(isym,mu)
          mu2=abs(mu1)
          if(mu1.lt.0) xxx = -xxx
          lam1=ifp(isym,lam)
          lam2=abs(lam1)
          if(mu2.eq.mu.and.lam2.eq.lam) cycle
          if(lam1.lt.0) xxx = -xxx
          if(ipair(isym).lt.0) xxx = -xxx
          if(jpair(isym).lt.0) xxx = -xxx
          xmatAA(mu2,lam2,lp0)=xxx
          xmatAA(lam2,mu2,lp1)=xxx
        end do
c
 10   CONTINUE
 11   CONTINUE
C
C  finished this batch of integrals
C  are there any more?
C
      lrec = lrec + numrec
      If(lrec.LE.mrec) GO TO 100
C
cccccccccc
cc      write(6,*) ' About to leave <BinSrt1AA>  xmatAA matrix is:'
cc      do ij=istrt,iend
cc      write(6,*) ' IJ is:',ij
cc      call prntmat(ncf,ncf,ncf,xmatAA(1,1,ij))
cc      enddo
cccccccccc
      RETURN
      END
c=====================================================================
c
      SUBROUTINE BinSendAB(nvalA,  nvalB,  nrec,    ndisk, nIntPerBin,
     $                  nBinPerRec,numrec, nrem,    istrt,   iend,
     $                     intsAB, int1AB, intsBA,  int1BA,  JMu,
     $                     JLam,   intsIJ, int1IJ,  intsJI,  int1JI,
     $                     mdisk4, lrec)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine sorts the half-transformed integrals from mu,lam all IJ to
C  IJ all mu,lam and writes them to the bin sort file
C
C
C  ARGUMENTS
C
C  nvalA       -  number of occupied alpha MOs to be correlated
C  nvalB       -  number of occupied beta  MOs to be correlated
C  nrec        -  number of records on the half-transformed integral file
C  ndisk       -  unit number of half-transformed integral file
C  nIntPerBin  -  number of 5-byte integrals that can be stored in each IJ bin
C                 (i.e., the bin size)
C  nBinPerRec  -  number of IJ bins written to the bin sort file per record
C  numrec      -  number of records written to bin sort file at one time
C  nrem        -  number of "padded" IJ bins in last record
C  istrt       -  starting IJ index on this pass
C  iend        -  ending IJ index on this pass
C  intsAB      -  4-byte storage for half-transformed alpha-beta integrals
C  int1AB      -    ditto  1-byte overflow storage
C  intsBA      -  4-byte storage for half-transformed beta-alpha integrals
C  int1BA      -    ditto  1-byte overflow storage
C  JMu         -  2-byte storage for mu AO-indices
C  JLam        -  2-byte storage for lamda AO-indices
C  intsIJ      -  4-byte storage for reordered alpha-beta integrals
C  int1IJ      -    ditto  1-byte overflow storage
C  intsJI      -  4-byte storage for reordered beta-alpha integrals
C  int1JI      -    ditto  1-byte overflow storage
C  mdisk4      -  unit number of bin sort file
C  lrec        -  on exit number of records on bin sort file
C  
C
      INTEGER*4 intsAB(nvalA,nvalB),intsIJ(nIntPerBin,istrt:iend+nrem),
     $          intsBA(nvalB,nvalA),intsJI(nIntPerBin,istrt:iend+nrem)
      INTEGER*1 int1AB(nvalA,nvalB),int1IJ(nIntPerBin,istrt:iend+nrem),
     $          int1BA(nvalB,nvalA),int1JI(nIntPerBin,istrt:iend+nrem)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
C
C  Read all records on the half-transformed integral file and select the range of
C  IJ indices being processed this pass
C
      IT = 0
      lrec = 0      !  number of records on bin sort file
c
      DO 20 irec=1,nrec
      READ (UNIT=ndisk,rec=irec) mu,lam,int1AB,intsAB,int1BA,intsBA
c
      IT = IT+1
      JMu(IT) = mu
      JLam(IT) = lam
c
      DO 10 IJ=istrt,iend
      call get_ij_full(IJ,nvalB,I,J)
      intsIJ(IT,IJ) = intsAB(I,J)
      int1IJ(IT,IJ) = int1AB(I,J)
      intsJI(IT,IJ) = intsBA(J,I)
      int1JI(IT,IJ) = int1BA(J,I)
 10   CONTINUE
c
      If(IT.EQ.nIntPerBin) Then
c
c -- write out the full bins
        NonZero = nIntPerBin
        IJ = 1
        Do K=1,numrec
        CALL BinWrtAB(mdisk4, nIntPerBin, nBinPerRec, JMu,    JLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),
     $                intsJI(1,IJ),  int1JI(1,IJ),  NonZero,  lrec)
        IJ = IJ+NBinPerRec
        EndDO
c
        IT = 0
      EndIf
 20   CONTINUE
C
C  All integrals processed
C  write out the last bins
C
      If(IT.GT.0) Then
c
        NonZero = IT
        IJ = 1
        Do K=1,numrec
        CALL BinWrtAB(mdisk4, nIntPerBin, nBinPerRec, JMu,    JLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),
     $                intsJI(1,IJ),  int1JI(1,IJ),  NonZero,  lrec)
         IJ = IJ+NBinPerRec
        EndDO
      EndIf
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE BinWrtAB(ndisk,  nIntPerBin, nBinPerRec, JMu,    JLam,
     $                    intsIJ, int1IJ, intsJI, int1JI,NonZero, lrec)
      IMPLICIT INTEGER(A-Z)
C
C  Writes a record to the bin sort file
C
      INTEGER*4 intsIJ(nIntPerBin,nBinPerRec)
      INTEGER*1 int1IJ(nIntPerBin,nBinPerRec)
      INTEGER*4 intsJI(nIntPerBin,nBinPerRec)
      INTEGER*1 int1JI(nIntPerBin,nBinPerRec)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
c
      lrec = lrec+1
      WRITE(unit=ndisk,rec=lrec) JMu,JLam,int1IJ,intsIJ,int1JI,intsJI,
     $                           NonZero
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE BinRdAB(ndisk,  nIntPerBin, nBinPerRec, JMu,    JLam,
     $                   intsIJ, int1IJ, intsJI, int1JI,NonZero, lrec)
      IMPLICIT INTEGER(A-Z)
C
C  Reads a record from the bin sort file
C
      INTEGER*4 intsIJ(nIntPerBin,nBinPerRec)
      INTEGER*1 int1IJ(nIntPerBin,nBinPerRec)
      INTEGER*4 intsJI(nIntPerBin,nBinPerRec)
      INTEGER*1 int1JI(nIntPerBin,nBinPerRec)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
c
      READ(unit=ndisk,rec=lrec) JMu,JLam,int1IJ,intsIJ,int1JI,intsJI,
     $                          NonZero
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE BinSrt1AB(ncf,    nvalA,   nvalB,  nsym,   ndisk,
     $                 nIntPerBin,nBinPerRec,numrec,mrec,   iloop,
     $                   iorbprA, iorbprB,  ifp,    thresh, istrt,
     $                     iend,   intsIJ,  int1IJ, intsJI, int1JI,
     $                     JMu,    JLam,    xmatAB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads back the alpha-alpha (or beta-beta) sorted half-transformed integrals
C  and puts then into the xmatAA matrix as reals
C
C  ARGUMENTS
C
C  ncf         -  number of basis functions
C  nvalA       -  number of occupied alpha MOs to be correlated
C  nvalB       -  number of occupied beta  MOs to be correlated
C  nsym        -  number of Abelian symmetry operations (less identity)
C  ndisk       -  unit number of bin sort file
C  nIntPerBin  -  number of 5-byte integrals that can be stored in each IJ bin
C                 (i.e., the bin size)
C  nBinPerRec  -  number of IJ bins written to the bin sort file per record
C  numrec      -  number of records written to bin sort file at one time
C  mrec        -  total number of records on the bin sort file
C  iloop       -  pointer to the record we are currently dealing with
C  iorbprA     -  symmetry characteristics of each correlated alpha MO
C                 under each Abelian symmetry operation
C  iorbprA     -    ditto but for the beta MOs
C  ifp         -    ditto but for the basis functions
C  thresh      -  threshold for integral neglect
C  istrt       -  starting IJ index on this pass
C  iend        -  ending IJ index on this pass
C  intsIJ      -  4-byte storage for reordered alpha-beta integrals
C  int1IJ      -    ditto  1-byte overflow storage
C  intsJI      -  4-byte storage for reordered beta-alpha integrals
C  int1JI      -    ditto  1-byte overflow storage
C  JMu         -  2-byte storage for mu AO-indices
C  JLam        -  2-byte storage for lamda AO-indices
C  xmatAB      -  reordered alpha-beta integrals ready for transformation
C
C
      INTEGER*4 intsIJ(nIntPerBin,istrt:iend) 
      INTEGER*1 int1IJ(nIntPerBin,istrt:iend)
      INTEGER*4 intsJI(nIntPerBin,istrt:iend) 
      INTEGER*1 int1JI(nIntPerBin,istrt:iend)
      INTEGER*2 JMu(NIntPerBin),JLam(nIntPerBin)
      INTEGER*4 NonZero
      REAL*8 xmatAB(ncf,ncf,istrt:iend)
      DIMENSION iorbprA(nsym,nvalA),iorbprB(nsym,nvalB),
     $          ifp(7,ncf),ipair(7),jpair(7)
c
c-----------------------------------------------------------------------
      Parameter (dblmax=2 147 483 648.0d0)
c-----------------------------------------------------------------------
      dblcmp = dblmax*thresh
c-----------------------------------------------------------------------
C
C  initialize
C
      numIJ = iend-istrt+1
      Call ZeroIT(xmatAB,numIJ*ncf**2)
c
      lrec = iloop      !  starting record number for this batch of bins
 100  CONTINUE
c
      Call BinRdAB(ndisk,  nIntPerBin, nBinPerRec, JMu,    JLam,
     $             intsIJ(1,istrt), int1IJ(1,istrt),
     $             intsJI(1,istrt), int1JI(1,istrt), NonZero, lrec)
c
c -- alpha-beta
c
      DO 11 IJ=istrt,iend
      call get_ij_full(IJ,nvalB,I,J)
c
      do isym=1,nsym
      ipair(isym)=iorbprA(isym,I)
      jpair(isym)=iorbprB(isym,J)
      end do
c
      DO 10 kk=1,NonZero
      mu = JMu(kk)
      lam = JLam(kk)
c
c -- decompress the integral ----------
cc        call reconstruct(intsIJ(kk,IJ),int1IJ(kk,IJ),xx)
cc        xx = xx*thresh
        If(int1IJ(kk,IJ).eq.0) Then
          xx = intsIJ(kk,IJ)*thresh
        Else If(int1IJ(kk,IJ).gt.0) Then
          xx = intsIJ(kk,IJ)*thresh
          xx = xx + SIGN(int1IJ(kk,IJ)*dblcmp,xx)
        Else
          xx = intsIJ(kk,IJ)*thresh*10.0d0**(-int1IJ(kk,IJ))
        EndIf
c ------------------------------------------------------
        xmatAB(mu,lam,IJ)=xx
c
        do isym=1,nsym
          xxx=xx
          mu1=ifp(isym,mu)
          mu2=abs(mu1)
          if(mu1.lt.0) xxx = -xxx
          lam1=ifp(isym,lam)
          lam2=abs(lam1)
          if(mu2.eq.mu.and.lam2.eq.lam) cycle
          if(lam1.lt.0) xxx = -xxx
          if(ipair(isym).lt.0) xxx = -xxx
          if(jpair(isym).lt.0) xxx = -xxx
          xmatAB(mu2,lam2,IJ)=xxx
        end do
c
 10   CONTINUE
 11   CONTINUE
c
c -- beta-alpha
c
      DO 21 IJ=istrt,iend
      call get_ij_full(IJ,nvalB,I,J)
c
      do isym=1,nsym
      ipair(isym)=iorbprB(isym,J)
      jpair(isym)=iorbprA(isym,I)
      end do
c
      DO 20 kk=1,NonZero
      mu = JMu(kk)
      lam = JLam(kk)
c
c -- decompress the integral ----------
cc        call reconstruct(intsJI(kk,IJ),int1JI(kk,IJ),xx)
cc        xx = xx*thresh
        If(int1JI(kk,IJ).eq.0) Then
          xx = intsJI(kk,IJ)*thresh
        Else If(int1JI(kk,IJ).gt.0) Then
          xx = intsJI(kk,IJ)*thresh
          xx = xx + SIGN(int1JI(kk,IJ)*dblcmp,xx)
        Else
          xx = intsJI(kk,IJ)*thresh*10.0d0**(-int1JI(kk,IJ))
        EndIf
c ------------------------------------------------------
        xmatAB(lam,mu,IJ)=xx
c
        do isym=1,nsym
          xxx=xx
          mu1=ifp(isym,mu)
          mu2=abs(mu1)
          if(mu1.lt.0) xxx = -xxx
          lam1=ifp(isym,lam)
          lam2=abs(lam1)
          if(mu2.eq.mu.and.lam2.eq.lam) cycle
          if(lam1.lt.0) xxx = -xxx
          if(ipair(isym).lt.0) xxx = -xxx
          if(jpair(isym).lt.0) xxx = -xxx
          xmatAB(lam2,mu2,IJ)=xxx
        end do
c
 20   CONTINUE
 21   CONTINUE
C
C  finished this batch of integrals
C  are there any more?
C
      lrec = lrec + numrec
      If(lrec.LE.mrec) GO TO 100
C
cccccccccc
cc      write(6,*) ' About to leave <BinSrt1AB>  xmatAB matrix is:'
cc      do ij=istrt,iend
cc      write(6,*) ' IJ is:',ij
cc      call prntmat(ncf,ncf,ncf,xmatAB(1,1,ij))
cc      enddo
cccccccccc
      RETURN
      END
c=====================================================================
c
      SUBROUTINE MemSortUMP2(membin,xmaxdisk, spaceU, spaceN,  ncf2,
     $                       npairs,  nfact,  nslv,    npass,  nIJ,
     $                   nIntPerBin,nBinPerRec, nrem,  numrec, incore)
      IMPLICIT INTEGER(A-Z)
      REAL*8 xmaxdisk,spaceU,spaceN,fraction
      Logical incore
C
C  This routine calculates various paramemers for the binsort step during
C  the second half-transformation based on the amount of disk storage and
C  memory available
C
C  ARGUMENTS
C
C  membin      -  available memory in double words (8-bytes)
C  xmaxdisk    -  available disk storage in double words
C  spaceU      -  size of existing integral files in double words
C  spaceN      -  estimated space needed for bins file in double words
C  ncf2        -  number of basis functions squared
C  npairs      -  total number of IJ pairs
C  nfact       -  1 for alpha-alpha & beta-beta term
C                 2 if dealing with alpha-beta term
C                 (latter has twice as many integrals per record)
C  nslv        -  number of slaves (must be 1 for serial)
C
C  on exit
C
C  npass       -  number of passes needed
C  nIJ         -  number of IJ pairs per pass
C  nIntPerBin  -  number of integrals per bin
C  nBinPerRec  -  number of bins per record on bin sort file
C  nrem        -  number of empty "padded" bins on final record
C  numrec      -  number of records written for each set of IJ bins
C  incore      -  logical flag for in-core execution
C
C
      incore = .False.
c
      If(nslv.GT.1) Then      !  parallel implementation
c -- for the time being we will assume that there is enough disk storage
c -- to write a full bin sort file across all the nodes
        npass = 1
        nIJ = npairs
        GO TO 95
      EndIf
c
      xmaxd = xmaxdisk-spaceU       !  estimated amount of disk storage left
cc      write(6,*) ' xmaxd is: ',xmaxd
c
      If(xmaxd.LT.zero) Then
        incore = .True.
        RETURN
      EndIf
c
      fraction = MIN(1.0d0,xmaxd/spaceN)
      If(fraction.LT.1.0d0) Then
        npass = 1.0d0/fraction + 1
        nIJ = npairs/npass          !  number of IJ indices per pass
        nrem = npairs - npass*nIJ
        If(nrem.GT.0) nIJ = nIJ+1
c -- nIJ must be even
        nn = nIJ/2
        nn = nIJ - 2*nn
        If(nn.GT.0) nIJ = nIJ-1
      Else
        npass = 1
        nIJ = npairs
      EndIF
c
   95 CONTINUE
c
c determine how many integrals can be written into each of the nIJ bins at any
c one time, i.e., the size of each bin
c
      nIntPerBin = (8*membin)/(5*nfact*nIJ+4)
      If(nIntPerBin.GT.ncf2) nIntPerBin = ncf2
c
c determine the number of IJ indices that can be read back in the available
c memory at any one time - this is the number of bins to write per record
c
      nBinPerRec = (8*membin-4*nIntPerBin)/(8*ncf2+5*nfact*nIntPerBin)
      If(nBinPerRec.GT.nIJ) nBinPerRec = nIJ
c -- nBinPerRec should be even
      If(nBinPerRec.LT.nIJ) Then
        nn = nBinPerRec/2
        nn = nBinPerRec - 2*nn
        If(nn.GT.0) nBinPerRec = nBinPerRec-1
        If(nBinPerRec.LT.2) call nerror(15,'UMP2 Module',
     $  'Too few IJ indices can be dealt with at once',0,nBinPerRec)
      EndIf
c
c if we are writing nBinPerRec bins at at time and there are nIJ bins determine
c the number of records that will be written and the remainder
c (i.e., the last record may not have its full complement of bins)
c
      numrec = nIJ/nBinPerRec
      nrem = nIJ - numrec*nBinPerRec
      If(nrem.GT.0) Then
        numrec = numrec+1
        nrem = nBinPerRec-nrem
      EndIf
c
c  The maximum amount of memory we are actually going to allocate is
c    5*nfact*nIntPerBin*(nIJ+nrem) + 4*nIntPerBin + 8*ncf2*nBinPerRec
c  Make sure this amount is not exceeded as we have not allowed for nrem
c  If it is, reduce nIntPerBin and recalculate
c
   96 CONTINUE
      memreal = 5*nfact*nIntPerBin*(nIJ+nrem) + 4*nIntPerBin
     $                 + 8*ncf2*nBinPerRec
cc      write(6,*) ' Memory available is: ',8*membin
cc      write(6,*) ' Memory allocated is: ',memreal
cc      write(6,*) ' nIntPerBin is: ',nIntPerBin
cc      write(6,*) ' nBinPerRec is: ',nBinPerRec
cc      write(6,*) ' numrec is: ',numrec
cc      write(6,*) ' nrem is: ',nrem
      If(memreal.GT.8*membin) Then
        nIntPerBin = nIntPerBin-1
        GO TO 96
      EndIf
c
      RETURN
      END
c=====================================================================
c
      SUBROUTINE MemSrt1UMP2(membin, ncf2, npairs, nIJ, npass)
      IMPLICIT INTEGER(A-Z)
C
C  This routine calculates the number of IJ indices that can be handled in core
C  at once based on the amount of memory available
C
C  ARGUMENTS
C
C  membin  -  available memory in double words (8-bytes)
C  ncf2    -  number of basis functions squared
C  npairs  -  total number of IJ pairs
C
C  on exit
C
C  nIJ     -  number of IJ pairs per pass
C  npass   -  number of passes needed
C
      nIJ = membin/ncf2
      If(nIJ.GT.npairs) nIJ = npairs
c -- nIJ must be even (unless last pass)
      If(nIJ.LT.npairs) Then
        nn = nIJ/2
        nn = nIJ - 2*nn
        If(nn.GT.0) nIJ = nIJ-1
      EndIf
c
      npass = npairs/nIJ
      nn = npairs - nIJ*npass
      If(nn.GT.0) npass = npass+1
c
      RETURN
      END
