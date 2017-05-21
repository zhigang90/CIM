      SUBROUTINE para_ump2(NAlpha,NBeta)
c
      use memory
      use newpara
c
c------------------------------------------------------------------------------
c
c The Parallel UMP2 program for the Unrestricted Hartree-Fock wavefunction
c
c Based on the serial code by JB 2009 and KW 2010
c
c JB 2010
c------------------------------------------------------------------------------
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
      integer*1 i01
      integer*4 i0,restrt
      INTEGER*4 SLAVE_ID(nslv)
      character*256 jobname,scrfile,filename,filname1,filname2
      logical LastSymPair,emp2only,scs,incore0,incore
c
      dimension xintxx(9)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      common /job/jobname,lenJ
c
      parameter(nopt=16)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      data smal/.TRUE./
      character*256 chopv(nopt)
      character*3 ch3
      character opnames*4,scftype*11,wvfnc*20,cdum*20
      logical nofr,exst,dualbasis,smal,keep,Test,alphabeta
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
ckw
      dimension Ipass(2,28), Kpass(2,28) !28 is 28 comp.of cartisian I-function
c-----------------------------------------------------------------------

      data opnames  /'nofr','orba','orbb','thre','core',
     1               'rest','dual','maxd','keep','prin',
     2               'grad','pmij','scs ','sos ','noio','test'/
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c  11=a single real, 12=two reals, 13=three reals, 21=character option
c
      data ioptyp / 0, 2, 2,11,11, 1, 1,11, 0, 1, 0, 2,12, 0, 0, 1/
c
c  Explanation of the options:
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
c  MAXD = maximum disk space in gigabytes (default 20 GB)
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
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c  get symetry info: number of symm. op. starting address of contr. info
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
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filname1=scrfile(1:len)//'.htr'
      len1 = len+4
      filname2=scrfile(1:len)//'.bins'
      len2 = len+5
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
        xmaxdisk=125 000 000.0d0*xmaxdisk     ! convert to megawords
      else
c  default disk storage is 20 GB
        xmaxdisk=2 500 000 000.0d0
      end if
      keep = ifound(9).gt.0
      if(ifound(10).gt.0) iprnt=iopv(1,10)
c  'PMIJ' print matrix; first number is 1 (AO Xchange) or 2
c   (both AO & MO) second number: pair to be printed
      ipmij=0
      ipairpr=0
      if(ifound(12).gt.0) then
        ipmij=iopv(1,12)
        ipairpr=iopv(2,12)
      end if
c -- SCS scaling
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
cc      If(ifound(12).gt.0) Then
cc        emp2only=.false.
cc        mp2only=-1
cc        filname3=scrfile(1:len)//'.Tij'
cc        len3 = len+4
cc        filname4=scrfile(1:len)//'.Kov'
cc        len4 = len+4
cc      EndIf
      call setival('mp2only',mp2only)
c-----------------------------------------------------------
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
c
      ncf2 = ncf**2
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
      write(iout,*) ' max disk storage (MW)  = ',xmaxdisk*1.d-6
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
      iforwhat=5
      thint=thresh
      call para_jobinit(iforwhat)
      call ptwoint(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .false.,NAlpha, NBeta,
     *             0,      scftype,xintxx, nblocks,maxbuffer,
     *             maxlabels)
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
c  print out some information
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
      call  matsub('occuA','cmoA',nfirstA,nlastA)
      call  matsub('virtA','cmoA',NAlpha+1,lastvirt)
      ioccA=mataddr('occuA')
      ivrtA=mataddr('virtA')
c
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
c ......................................................................
c -- determine length of direct access file records
c
      lrec1 = 5*(nvalA**2) + 16
      lrec2 = 5*(2*nvalA*nvalB) + 16
      lrec3 = 5*(nvalB**2) + 16
c
      write(iout,1235) lrec1,lrec2,lrec3
 1235 format(' Record length for integral files (bytes): ',3(I10,1x))
c ......................................................................
c
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
c
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c
      ENonZero=zero
c
      call init_info('fock')
      call elapsec(tstart)
c ....................................................................
c
c -- send initial UMP2 data to slaves
      call para_initsend
      call para_pack_int(nfirstA,1)
      call para_pack_int(nlastA,1)
      call para_pack_int(nvalA,1)
      call para_pack_int(lrec1,1)
      call para_pack_int(nfirstB,1)
      call para_pack_int(nlastB,1)
      call para_pack_int(nvalB,1)
      call para_pack_int(lrec2,1)
      call para_pack_int(lrec3,1)
      call para_pack_byte(alphabeta,1)
      call para_pack_int(smal,1)
      call para_pack_int(keep,1)
      call para_pack_int(restrt,1)
      call para_pack_int(Test,1)
      call para_pack_int(mp2only,1)
      call para_pack_int(bl(mpres_in_bg),ncs)
      call para_pack_real(thresh,1)
      call para_pack_real(bl(ieorbA),ncf)
      call para_pack_real(bl(icanoA),ncf*ncf)
      call para_pack_real(bl(ieorbB),ncf)
      call para_pack_real(bl(icanoB),ncf*ncf)
      call para_bcast_pack(TxMP2Dat)
c
c -- send screening density, symm.info and filename separately
      call para_initsend
      call para_pack_real(bl(idsmxA),ncs*ncs)
      If(nsym.gt.0) Then
        call para_pack_int(bl(iorbprA),nvalA*nsym)
        call para_pack_int(bl(iorbprB),nvalB*nsym)
      EndIf
c
c -- take care with file handling
c -- send each slave the base filenames
c -- they will make it unique based on their gid (1...nslv)
      call para_pack_string(filname1,256)
      call para_pack_int(len1,1)
      call para_pack_string(filname2,256)
      call para_pack_int(len2,1)
cc      If(.NOT.emp2only) Then
cc         call para_pack_string(filname3,256)
cc         call para_pack_int(len3,1)
cc         call para_pack_string(filname4,256)
cc         call para_pack_int(len4,1)
cc      End If
      call para_bcast_pack(TxMP2File)
c .....................................................................
c
      if(restrt) go to 50
c .....................................................................
c
c    F I R S T   H A L F - T R A N S F O R M A T I O N
c    -------------------------------------------------
c
      do ics=ncs,1,-1
      do kcs=ics,1,-1
c  of each (ics,kcs) pair, calculate only the last one
c  (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
      if(LastSymPair(ics,kcs,nsym,bl(ifp1),inegl,iret)) cycle
c
c -- send ics,kcs pair to slave
        call para_recv(imsg,islave,TxDftReq)
        call para_initsend
        call para_pack_int(ics,1)
        call para_pack_int(kcs,1)
        call para_send_pack(islave,TxMP2Assign)
      end do
      end do
c
c -- all done  stop all slaves
      call elapsec(timo)
      call getrval('imbal',tlost)
      Do ics=1,nslv
        call para_recv(imsg,islave,TxDftReq)
        call para_initsend
        call para_pack_int(0,1)
        call para_pack_int(0,1)
        call para_send_pack(islave,TxMP2Assign)
      call elapsec(timb)
      tlost=tlost+(ics-1)*(timb-timo)
      timo=timb
      EndDo
      call setrval('imbal',tlost)
      call elapsec(tend)
c
      If(iprnt.gt.0) Then
        telap = (tend-tstart)/sixty
        write(iout,1000) telap
 1000   Format(/,' End of First Half-transformation    elapsed time = ',
     $           F10.2,' min',/)
        call f_lush(iout)
      EndIf
c
c .....................................................
      If(Test) Then
        write(iout,65)
   65   format('  == Calculation Halted at User Request ==')
        STOP
       EndIf
c .....................................................
c
 50   CONTINUE
c
c
c *******************************************************************
c
c  SECOND HALF-TRANSFORMATION
c
c *******************************************************************
c
c  Integrals from the first half-transformation are stored as ALL IJ values for
c  each Mu,Lam (Mu >= Lam). They need to be resorted as ALL Mu,Lam for each IJ.
c  This will be done by a modification to the existing bin sort algorithm.
c  
      call elapsec(tstart)
c
c -- check if the slaves are still OK
      call para_check
c
c  tell the slaves to start
c
      IStart = 1
      call para_initsend
      call para_pack_int(IStart,1)
      call para_bcast_pack(TxBinDat)
c
c -- get from each slave the size of the half-transformed integral file
      DO islv=1,nslv
      call para_recv(nrec,islave,TxMP2S4)
      If(iprnt.gt.0) write(iout,1236) islave,nrec
      EndDo
 1236 Format(' Slave: ',I4,' Half-transformed integral file has ',I8,
     $         ' records')
c
c -- get back from each slave the amount of memory available for the bin sort
C    (this should be the same on each slave)
c -- at the same time store the slave IDs
c
      call para_recv(membin,islave,TxDftReq)
      SLAVE_ID(1) = islave
c
      DO islv=2,nslv
      call para_recv(memb,islave,TxDftReq)
      SLAVE_ID(islv) = islave
      If(memb.LT.membin) membin = memb
      EndDO
c
c -- broadcast the slave ID array
      call para_initsend
      call para_pack_int4(SLAVE_ID,nslv)
      call para_bcast_pack(TxMP2S5)
c
c -- the master now does nothing except wait for the slaves to finish
c -- each term after which the partial UMP2 energies are summed up
c
      If(alphabeta) GO TO 999
C
C
C  (1)  ALPHA-ALPHA
C
C  based on the memory available on the slaves, compute the various
C  parameters for the bin sort
C
      call elapsec(t0)
c
      npairs = nvalA**2
c
      CALL MemSortUMP2(membin,   rjnk,   rjnk,      rjnk,     ncf2,
     $                 npairs,   1,      nslv,      npass,    nIJ,
     $              nIntPerBin,nBinPerRec,  nrem,   numrec,   incore)
c
c -- now broadcast the bin sort parameters to the slaves
      call para_initsend
      call para_pack_int(nIJ,1)
      call para_pack_int(nIntPerBin,1)
      call para_pack_int(nBinPerRec,1)
      call para_bcast_pack(TxMP2S4)
c
      If(iprnt.gt.0) Then
        write(iout,*)
        write(iout,*) ' --- Second half alpha-alpha transformation ---'
        write(iout,1237) npass,nIJ,nIntPerBin,nBinPerRec
 1237 format(/,' Number of passes (bin sort files written): ',I4,/,
     $         ' Maximum number of IJ indices per pass:   ',I6,/,
     $         ' Maximum number of integrals per bin:   ',I8,/,
     $         ' Maximum number of IJ indices in memory:  ',I6,/)
      EndIf
c
c -- accumulate alpha-alpha term
      eAA = zero
      call para_reduce(eAA,   1,TxMP2S1)
c
      call elapsec(t1)
      ttime = (t1-t0)/sixty
c
      If(iprnt.gt.0) WRITE(iout,2500) ttime
 2500 Format(1X,' Elapsed time to compute alpha-alpha energy is: ',
     $            F10.2,' mins',/)
C
C
C  (2)  BETA-BETA
C
C  based on the memory available on the slaves, compute the various
C  parameters for the bin sort
C
      call elapsec(t0)
c
      npairs = nvalB**2
c
      CALL MemSortUMP2(membin,   rjnk,   rjnk,      rjnk,     ncf2,
     $                 npairs,   1,      nslv,      npass,    nIJ,
     $              nIntPerBin,nBinPerRec,  nrem,   numrec,   incore)
c
c -- now broadcast the bin sort parameters to the slaves
      call para_initsend
      call para_pack_int(nIJ,1)
      call para_pack_int(nIntPerBin,1)
      call para_pack_int(nBinPerRec,1)
      call para_bcast_pack(TxMP2S5)
c
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half beta-beta transformation ---'
        write(iout,1237) npass,nIJ,nIntPerBin,nBinPerRec
      EndIf
c
c -- accumulate beta-beta term
      eBB = zero
      call para_reduce(eBB,   1,TxMP2S2)
c
      call elapsec(t1)
      ttime = (t1-t0)/sixty
c
      If(iprnt.gt.0) WRITE(iout,2600) ttime
 2600 Format(1X,' Elapsed time to compute beta-beta energy is:   ',
     $            F10.2,' mins',/)
c
  999 CONTINUE
C
C
C  (3)  ALPHA-BETA
C
C  based on the memory available on the slaves, compute the various
C  parameters for the bin sort
C
      call elapsec(t0)
c
      npairs = nvalA*nvalB
c
      CALL MemSortUMP2(membin,   rjnk,   rjnk,      rjnk,     ncf2,
     $                 npairs,   2,      nslv,      npass,    nIJ,
     $              nIntPerBin,nBinPerRec,  nrem,   numrec,   incore)
c
c -- now broadcast the bin sort parameters to the slaves
      call para_initsend
      call para_pack_int(nIJ,1)
      call para_pack_int(nIntPerBin,1)
      call para_pack_int(nBinPerRec,1)
      call para_bcast_pack(TxMP2S6)
c
      If(iprnt.gt.0) Then
        write(iout,*) ' --- Second half alpha-beta transformation ---'
        write(iout,1237) npass,nIJ,nIntPerBin,nBinPerRec
      EndIf
c
c -- accumulate alpha-beta term
      eAB = zero
      call para_reduce(eAB,   1,TxMP2S3)
c
      call elapsec(t1)
      ttime = (t1-t0)/sixty
c
      If(iprnt.gt.0) WRITE(iout,2700) ttime
 2700 Format(1X,' Elapsed time to compute alpha-beta energy is:  ',
     $            F10.2,' mins',/)
c
c-----------------------------------------------------------------------
c
c -- that's it!
c
c -- get detailed timings from all the slaves
      tbin = zero
      tsort = zero
      ttran = zero
      Do islv=1,nslv
      call para_recv_pack(ifrom,TxMP2S5)
      call para_unpack_real(tb,1)
      call para_unpack_real(ts,1)
      call para_unpack_real(tt,1)
      If(tb.GT.tbin) tbin = tb
      If(ts.GT.tsort) tsort = ts
      If(tt.GT.ttran) ttran = tt
      EndDo
c
      If(iprnt.gt.0) WRITE(iout,2800) tbin/sixty,tsort/sixty,ttran/sixty
 2800 Format(
     $    ' Maximum time spent in bin sort:           ',F10.2,' mins',/,
     $    ' Maximum time spent unpacking integrals:   ',F10.2,' mins',/,
     $    ' Maximum time in integral transformation:  ',F10.2,' mins',/)
C
C  now sum up the wavefunction norm
C
      tnorm = one
      call para_reduce(tnorm, 1,TxMP2S4)
c
      call elapsec(tend)
c
      If(iprnt.gt.0) Then
        telap = (tend-tstart)/sixty
        write(iout,1200) telap
 1200   Format(' End of Second Half-transformation    elapsed time = ',
     $           F10.2,' min')
        call f_lush(iout)
      EndIf
c
c -- timing stuff
      call para_next(-1)
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
        Call wrcntrl(IUnit,7,'$energy',2,idum,euhf+emp2scal,wvfnc)
      Endif
      If(dualbasis) Call wrcntrl(IUnit,6,'$edual',2,idum,edual,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
      call retmark
c---------------------------------------------------------------------
      call setival('ncoreA',ncoreA)
      call setival('ncoreB',ncoreB)
      call setrval('thrsx',thresh)
c---------------------------------------------------------------------
      call memory_status('end of mp2')
c---------------------------------------------------------------------
      end
c=====================================================================
c
      SUBROUTINE Para_BinSendAA(
     $                     nval,   nrec,   ndisk,nIntPerBin,nBinPerRec,
     $                     numrec, nrem,   istrt,   iend,    intsAA,
     $                     int1AA, JMu,    JLam,    intsIJ,  int1IJ,
     $                     listAA, mdisk4, SLAVE_ID,IBIN,    IJSTRT,
     $                     MyID,   lrec,    jstrt,  jend)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine sorts the half-transformed integrals from mu,lam all IJ to
C  IJ all mu,lam and writes them to the bin sort file
C
C
C  ARGUMENTS
C
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
C  SLAVE_ID    -  slave ID array
C  IBIN        -  number of IJ indices per slave
C  IJSTRT      -  starting IJ index for each slave
C  MyID        -  position in SLAVE_ID array of THIS slave
C
C  on exit
C
C  lrec        -  number of records on bin sort file
C  jstrt       -  starting IJ index to be transformed on this slave
C  jend        -    ditto  ending IJ index
C
C
      INTEGER*4 intsAA(nval,nval),intsIJ(nIntPerBin,istrt:iend+nrem)
      INTEGER*1 int1AA(nval,nval),int1IJ(nIntPerBin,istrt:iend+nrem)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*2 KMu(nIntPerBin),KLam(nIntPerBin)
      INTEGER*4 SLAVE_ID(nslv),IBIN(nslv),IJSTRT(nslv)
      INTEGER*4 NonZero,NonZero1
      LOGICAL   DONE(nslv),DONE1(nslv)
c-----------------------------------------------------------------------
      dimension listAA(*)        ! nvalA*nvalA
c-----------------------------------------------------------------------
C
C
C  initialize the DONE array
C
      Do islv=1,nslv
      DONE(islv) = .False.
      DONE1(islv) = .False.
      EndDo
C
C  Read all records on the half-transformed integral file and select the range of
C  IJ indices being processed this pass
C
      IT = 0
      irec = 0          !  number of records on integral file
      lrec = 0          !  number of records on bin sort file
      NLeft = nslv-1    !  number of slaves still with data to send
C
C  Rewrite due to failure of original algorithm under MPI
C  Worked fine under PVM
C
  5   CONTINUE
C
C  Fill One complete set of bins and then process that
C
      DO 20 IT=1,nIntPerBin
c
      irec = irec+1
      READ (UNIT=ndisk,rec=irec) mu,lam,int1AA,intsAA
c
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
      If(irec.EQ.nrec) GO TO 25      ! all integrals read
c
 20   CONTINUE
      IT = nIntPerBin
 25   CONTINUE
C
C  At this point we have a complete set of bins ready to send
C
      NonZero = IT                   ! no. of integrals per bin
c
c -- first write out the bins that will be processed on this slave
c
      IJ = IJSTRT(MyID)
      Do K=1,numrec
      CALL BinWrtAA(mdisk4, nIntPerBin, nBinPerRec, JMu,    JLam,
     $              intsIJ(1,IJ),  int1IJ(1,IJ),  NonZero,  lrec)
      IJ = IJ+NBinPerRec
      EndDo
c
c -- now send the other bins to the slave that will process them
c    and get the bins for this slave
c
      DO 40 islv=1,nslv
      If(MY_GID.EQ.SLAVE_ID(islv)) Then
c
c -- receive
c
        DO 30 jslv=1,nslv
        If(DONE(jslv)) GO TO 30
        kslave = SLAVE_ID(jslv)
        If(kslave.EQ.MY_GID) GO TO 30
        IJ = IJSTRT(MyID)
        nIJ = IBIN(MyID)
c -- receive data from slave
        Call para_recv(imsg,islave,TxBlockReq)
        Call para_recv_pack(islave,TxBinDat)
        Call para_unpack_int4(NonZero1,1)
        Call para_unpack_int2(KMu,nIntPerBin)
        Call para_unpack_int2(KLam,nIntPerBin)
        Call para_unpack_byte(int1IJ(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_int4(intsIJ(1,IJ),nIntPerBin*nIJ)
        Do K=1,numrec
        CALL BinWrtAA(mdisk4, nIntPerBin, nBinPerRec, KMu,    KLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),  NonZero1, lrec)
        IJ = IJ+NBinPerRec
        EndDo
 30     CONTINUE
      Else
c
c -- send
c
        kslave = SLAVE_ID(islv)
        IJ = IJSTRT(islv)
        nIJ = IBIN(islv)
c
        Call para_send(MY_GID,kslave,TxBlockReq)
        Call para_initsend
        Call para_pack_int4(NonZero,1)
        Call para_pack_int2(JMu,nIntPerBin)
        Call para_pack_int2(JLam,nIntPerBin)
        Call para_pack_byte(int1IJ(1,IJ),nIntPerBin*nIJ)
        Call para_pack_int4(intsIJ(1,IJ),nIntPerBin*nIJ)
        Call para_send_pack(kslave,TxBinDat)
      EndIf
 40   CONTINUE
C
C  At this point we have sent/received one bin
C  Communicate to all slaves to continue with next bin
C  but let all slaves know if THIS slave is done
C
      If(irec.EQ.nrec) Then
        DONE(MyID) = .True.
        DONE1(MyID) = .True.
        ISend = 0
      Else
        ISend = 1
      EndIf
c
      DO 60 islv=1,nslv
      If(MY_GID.EQ.SLAVE_ID(islv)) Then
c
c -- receive
c
        DO 50 jslv=1,nslv
        If(DONE(jslv)) GO TO 50
        kslave = SLAVE_ID(jslv)
        If(kslave.EQ.MY_GID) GO TO 50
c -- receive continuation data from slave
        Call para_recv(ISend1,islave,TxMP2S4)
c -- check if slave has finished sending
        If(ISend1.EQ.0) Then
          Call Set_DONE(islave,nslv,SLAVE_ID,DONE1)
          NLeft = NLeft-1
        EndIf
 50     CONTINUE
c -- update DONE array
        Do jslv=1,nslv
        DONE(jslv) = DONE1(jslv)
        EndDo
      Else
c
c -- send
c
        kslave = SLAVE_ID(islv)
        Call para_send(ISend,kslave,TxMP2S4)
      EndIf
 60   CONTINUE
C
C  Continuation data sent
C  Either continue or drop into receive-only mode
C
      If(ISend.NE.0) GO TO 5
C
C  THIS slave is now done
C  Are any bins left to get from other slaves?
C  This will be the case if NLeft > 0
C
 75   CONTINUE
      If(NLeft.GT.0) Then
        DO 80 jslv=1,nslv
        If(DONE(jslv)) GO TO 80
        IJ = IJSTRT(MyID)
        nIJ = IBIN(MyID)
c -- receive data from slave
        Call para_recv(imsg,islave,TxBlockReq)
        Call para_recv_pack(islave,TxBinDat)
        Call para_unpack_int4(NonZero1,1)
        Call para_unpack_int2(KMu,nIntPerBin)
        Call para_unpack_int2(KLam,nIntPerBin)
        Call para_unpack_byte(int1IJ(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_int4(intsIJ(1,IJ),nIntPerBin*nIJ)
        Do K=1,numrec
        CALL BinWrtAA(mdisk4, nIntPerBin, nBinPerRec, KMu,    KLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),  NonZero1, lrec)
        IJ = IJ+NBinPerRec
        EndDo
 80     CONTINUE
c
c
c -- receive
c
        DO 85 jslv=1,nslv
        If(DONE(jslv)) GO TO 85
        kslave = SLAVE_ID(jslv)
        If(kslave.EQ.MY_GID) GO TO 85
c -- receive continuation data from slave
        Call para_recv(ISend1,islave,TxMP2S4)
c -- check if slave has finished sending
        If(ISend1.EQ.0) Then
          Call Set_DONE(islave,nslv,SLAVE_ID,DONE1)
          NLeft = NLeft-1
        EndIf
 85     CONTINUE
c -- update DONE array
        Do jslv=1,nslv
        DONE(jslv) = DONE1(jslv)
        EndDo
      EndIf
c
      If(NLeft.GT.0) GO TO 75
C
C  return starting and ending IJ indices for THIS slave
C
      jstrt = IJSTRT(MyID)
      If(MyID.EQ.nslv) Then
        jend = iend
      Else
        jend = IJSTRT(MyID+1)-1
      EndIf
C
      RETURN
      END
c=====================================================================
c
      SUBROUTINE Para_BinSendAB(
     $                     nvalA,  nvalB,  nrec,    ndisk,nIntPerBin,
     $                  nBinPerRec,numrec, nrem,    istrt,   iend,
     $                     intsAB, int1AB, intsBA,  int1BA,  JMu,
     $                     JLam,   intsIJ, int1IJ,  intsJI,  int1JI,
     $                     mdisk4, SLAVE_ID,IBIN,   IJSTRT,  MyID,
     $                     lrec,   jstrt,  jend)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine sorts the half-transformed integrals from mu,lam all IJ to
C  IJ all mu,lam and sends them to the slave that will ultimately transform them
C
C
C  ARGUMENTS
C
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
C  SLAVE_ID    -  slave ID array
C  IBIN        -  number of IJ indices per slave
C  IJSTRT      -  starting IJ index for each slave
C  MyID        -  position in SLAVE_ID array of THIS slave
C
C  on exit
C
C  lrec        -  number of records on bin sort file
C  jstrt       -  starting IJ index to be transformed on this slave
C  jend        -    ditto  ending IJ index
C
C
      INTEGER*4 intsAB(nvalA,nvalB),intsIJ(nIntPerBin,istrt:iend+nrem)
      INTEGER*1 int1AB(nvalA,nvalB),int1IJ(nIntPerBin,istrt:iend+nrem)
      INTEGER*4 intsBA(nvalB,nvalA),intsJI(nIntPerBin,istrt:iend+nrem)
      INTEGER*1 int1BA(nvalB,nvalA),int1JI(nIntPerBin,istrt:iend+nrem)
      INTEGER*2 JMu(nIntPerBin),JLam(nIntPerBin)
      INTEGER*2 KMu(nIntPerBin),KLam(nIntPerBin)
      INTEGER*4 SLAVE_ID(nslv),IBIN(nslv),IJSTRT(nslv)
      INTEGER*4 NonZero,NonZero1
      LOGICAL   DONE(nslv),DONE1(nslv)
C
C
C  initialize the DONE array
C
      Do islv=1,nslv
      DONE(islv) = .False.
      DONE1(islv) = .False.
      EndDo
C
C  Read all records on the half-transformed integral file and select the range of
C  IJ indices being processed this pass
C
      IT = 0
      irec = 0          !  number of records on integral file
      lrec = 0          !  number of records on bin sort file
      NLeft = nslv-1    !  number of slaves still with data to send
C
C  Rewrite due to failure of original algorithm under MPI
C  Worked fine under PVM
C
  5   CONTINUE
C
C  Fill One complete set of bins and then process that
C
      DO 20 IT=1,nIntPerBin
c
      irec = irec+1
      READ (UNIT=ndisk,rec=irec) mu,lam,int1AB,intsAB,int1BA,intsBA
c
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
      If(irec.EQ.nrec) GO TO 25      ! all integrals read
c
 20   CONTINUE
      IT = nIntPerBin
 25   CONTINUE
C
C  At this point we have a complete set of bins ready to send
C
      NonZero = IT                   ! no. of integrals per bin
c
c -- first write out the bins that will be processed on this slave
c
      IJ = IJSTRT(MyID)
      Do K=1,numrec
      CALL BinWrtAB(mdisk4, nIntPerBin, nBinPerRec, JMu,    JLam,
     $              intsIJ(1,IJ),  int1IJ(1,IJ),
     $              intsJI(1,IJ),  int1JI(1,IJ),  NonZero,  lrec)
      IJ = IJ+NBinPerRec
      EndDo
c
c -- now send the other bins to the slave that will process them
c    and get the bins for this slave
c
      DO 40 islv=1,nslv
      If(MY_GID.EQ.SLAVE_ID(islv)) Then
c
c -- receive
c
        DO 30 jslv=1,nslv
        If(DONE(jslv)) GO TO 30
        kslave = SLAVE_ID(jslv)
        If(kslave.EQ.MY_GID) GO TO 30
        IJ = IJSTRT(MyID)
        nIJ = IBIN(MyID)
c -- receive data from slave
        Call para_recv(imsg,islave,TxBlockReq)
        Call para_recv_pack(islave,TxBinDat)
        Call para_unpack_int4(NonZero1,1)
        Call para_unpack_int2(KMu,nIntPerBin)
        Call para_unpack_int2(KLam,nIntPerBin)
        Call para_unpack_byte(int1IJ(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_int4(intsIJ(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_byte(int1JI(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_int4(intsJI(1,IJ),nIntPerBin*nIJ)
        Do K=1,numrec
        CALL BinWrtAB(mdisk4, nIntPerBin, nBinPerRec, KMu,    KLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),
     $                intsJI(1,IJ),  int1JI(1,IJ),  NonZero1, lrec)
        IJ = IJ+NBinPerRec
        EndDo
 30     CONTINUE
      Else
c
c -- send
c
        kslave = SLAVE_ID(islv)
        IJ = IJSTRT(islv)
        nIJ = IBIN(islv)
c
        Call para_send(MY_GID,kslave,TxBlockReq)
        Call para_initsend
        Call para_pack_int4(NonZero,1)
        Call para_pack_int2(JMu,nIntPerBin)
        Call para_pack_int2(JLam,nIntPerBin)
        Call para_pack_byte(int1IJ(1,IJ),nIntPerBin*nIJ)
        Call para_pack_int4(intsIJ(1,IJ),nIntPerBin*nIJ)
        Call para_pack_byte(int1JI(1,IJ),nIntPerBin*nIJ)
        Call para_pack_int4(intsJI(1,IJ),nIntPerBin*nIJ)
        Call para_send_pack(kslave,TxBinDat)
      EndIf
 40   CONTINUE
C
C  At this point we have sent/received one bin
C  Communicate to all slaves to continue with next bin
C  but let all slaves know if THIS slave is done
C
      If(irec.EQ.nrec) Then
        DONE(MyID) = .True.
        DONE1(MyID) = .True.
        ISend = 0
      Else
        ISend = 1
      EndIf
c
      DO 60 islv=1,nslv
      If(MY_GID.EQ.SLAVE_ID(islv)) Then
c
c -- receive
c
        DO 50 jslv=1,nslv
        If(DONE(jslv)) GO TO 50
        kslave = SLAVE_ID(jslv)
        If(kslave.EQ.MY_GID) GO TO 50
c -- receive continuation data from slave
        Call para_recv(ISend1,islave,TxMP2S1)
c -- check if slave has finished sending
        If(ISend1.EQ.0) Then
          Call Set_DONE(islave,nslv,SLAVE_ID,DONE1)
          NLeft = NLeft-1
        EndIf
 50     CONTINUE
c -- update DONE array
        Do jslv=1,nslv
        DONE(jslv) = DONE1(jslv)
        EndDo
      Else
c
c -- send
c
        kslave = SLAVE_ID(islv)
        Call para_send(ISend,kslave,TxMP2S1)
      EndIf
 60   CONTINUE
C
C  Continuation data sent
C  Either continue or drop into receive-only mode
C
      If(ISend.NE.0) GO TO 5
C
C  THIS slave is now done
C  Are any bins left to get from other slaves?
C  This will be the case if NLeft > 0
C
 75   CONTINUE
      If(NLeft.GT.0) Then
        DO 80 jslv=1,nslv
        If(DONE(jslv)) GO TO 80
        IJ = IJSTRT(MyID)
        nIJ = IBIN(MyID)
c -- receive data from slave
        Call para_recv(imsg,islave,TxBlockReq)
        Call para_recv_pack(islave,TxBinDat)
        Call para_unpack_int4(NonZero1,1)
        Call para_unpack_int2(KMu,nIntPerBin)
        Call para_unpack_int2(KLam,nIntPerBin)
        Call para_unpack_byte(int1IJ(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_int4(intsIJ(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_byte(int1JI(1,IJ),nIntPerBin*nIJ)
        Call para_unpack_int4(intsJI(1,IJ),nIntPerBin*nIJ)
        Do K=1,numrec
        CALL BinWrtAB(mdisk4, nIntPerBin, nBinPerRec, KMu,    KLam,
     $                intsIJ(1,IJ),  int1IJ(1,IJ),
     $                intsJI(1,IJ),  int1JI(1,IJ),  NonZero1, lrec)
        IJ = IJ+NBinPerRec
        EndDo
 80     CONTINUE
c
c
c -- receive
c
        DO 85 jslv=1,nslv
        If(DONE(jslv)) GO TO 85
        kslave = SLAVE_ID(jslv)
        If(kslave.EQ.MY_GID) GO TO 85
c -- receive continuation data from slave
        Call para_recv(ISend1,islave,TxMP2S1)
c -- check if slave has finished sending
        If(ISend1.EQ.0) Then
          Call Set_DONE(islave,nslv,SLAVE_ID,DONE1)
          NLeft = NLeft-1
        EndIf
 85     CONTINUE
c -- update DONE array
        Do jslv=1,nslv
        DONE(jslv) = DONE1(jslv)
        EndDo
      EndIf
c
      If(NLeft.GT.0) GO TO 75
C
C  return starting and ending IJ indices for THIS slave
C
      jstrt = IJSTRT(MyID)
      If(MyID.EQ.nslv) Then
        jend = iend
      Else
        jend = IJSTRT(MyID+1)-1
      EndIf
C
      RETURN
      END
c=====================================================================
c
      SUBROUTINE Set_DONE(islave,nslv,SLAVE_ID,DONE)
      IMPLICIT INTEGER(A-Z)
C
C  Sets the logical entry in DONE for the slave islave
C
C  ARGUMENTS
C
C  islave    -  slave under consideration
C  nslv      -  total number of slaves
C  SLAVE_ID  -  array of slave IDs
C  DONE      -  logical array (entry set to .TRUE. if that slave is finished)
C
      INTEGER*4 SLAVE_ID(nslv)
      LOGICAL DONE(nslv)
c
      DO islv=1,nslv
      If(islave.EQ.SLAVE_ID(islv)) Then
        DONE(islv) = .True.
        RETURN
      EndIf
      EndDO
c
c -- should NEVER get here
      call nerror(17,'Parallel UMP2 Module',
     $      'Something is badly wrong -- lost track of Slave IDs',0,0)
c
      END
c=====================================================================
c
      SUBROUTINE Split_IJ(istrt,  iend,   nslv,   MY_GID, SLAVE_ID,
     $                nBinPerRec, IBIN,   IJSTRT, numrec, nrem,
     $                    MyID)
      IMPLICIT INTEGER(A-Z)
C
C  This routine divides the IJ indices being done this pass among the slaves
C  and determines numrec & nrem for THIS slave
C
C  ARGUMENTS
C
C  istrt       -  starting IJ index
C  iend        -  ending IJ index
C  nslv        -  number of slaves
C  MY_GID      -  ID of THIS slave
C  SLAVE_ID    -  ID array for all slaves
C  nBinPerRec  -  number of IJ bins written to the bin sort file per record
C
C  on exit
C
C  IBIN        -  number of IJ indices per slave
C  IJSTRT      -  starting IJ index on each slave
C  numrec      -  number of records written to bin sort file at one time on THIS slave
C  nrem        -  number of "padded" IJ bins in last record
C  MyID        -  position in SLAVE_ID array of THIS slave
C
C
      INTEGER*4 SLAVE_ID(nslv),IBIN(nslv),IJSTRT(nslv)
C
C  First determine which slaves do which IJ indices
C
      npairs = iend-istrt+1
      nIJ = npairs/nslv
c -- nIJ must be even
      nn = nIJ/2
      nn = nIJ - 2*nn
      If(nn.GT.0) nIJ = nIJ-1
c -- how many IJ pairs are left?
      nn = npairs - nslv*nIJ
c -- remaining IJ pairs are added to each bin 2 at a time until there are none
c -- or one left; if the latter this is added to the last bin
      Do islv=1,nslv
      If(MY_GID.EQ.SLAVE_ID(islv)) MYID = islv      ! ID of THIS slave
      IBIN(islv) = nIJ
      EndDo
 5    CONTINUE
      Do islv=1,nslv
      If(nn.GE.2) Then
        IBIN(islv) = IBIN(islv)+2
        nn = nn-2
      Else
        If(nn.EQ.1) IBIN(nslv) = IBIN(nslv)+1
        nn = nn-1
        EXIT
      EndIf
      EndDO
      If(nn.GT.0) GO TO 5
c
c -- fill the IJSTRT array
      IJSTRT(1) = 1
      Do islv=1,nslv-1
      IJSTRT(islv+1) = IJSTRT(islv) + IBIN(islv)
      EndDo
C
C  determine nrem (for last bin)
C
      nIJ = IBIN(nslv)
      numrec = nIJ/nBinPerRec
      nrem = nIJ - numrec*nBinPerRec
      If(nrem.GT.0) Then
        numrec = numrec+1
        nrem = nBinPerRec-nrem
      EndIf
C
C  determine numrec for THIS slave
C
      nIJ = IBIN(MyID)
      numrec = nIJ/nBinPerRec
      nn = nIJ - numrec*nBinPerRec
      If(nn.GT.0) numrec = numrec+1
C
      RETURN
      END
