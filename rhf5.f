      subroutine rhfmain(inp,GUES)

      use memory

      implicit real*8 (a-h,o-z)
c ....................................................................
c  Restricted Self-Consistent Field Module
c    revised     August 1998    (MultiFock REMOVED)
c ....................................................................
      Logical conver,conv,itlimit,interpol,switch,locali,localv
      Logical lpseudo,lpwave,escfonly,noabc,tz
      character lmeth*5,lmethv*5,scftype*11,dmethod*1,Sflies*4
      parameter (nopt=39)
      character*4 opnames(nopt)
      parameter(ndft=27)
      character*6 dftype,dftxc(ndft)
      character*20 wvfnc,GUES
      character*256 chopv(nopt),char,jobname,MOS,LOS
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension xintxx(9),dip(3),field(3)
      dimension xnmo(2),values(5)
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
c     common /mmanag/ lcore
      common /symm/nsym,nsy(7)
      common /errflag/ierrflag(10)
c   These are sums: nintsum1 sums up nintez, nintsum2 nblsiz*lsh and
c   nquartsum nblsiz   They must be zeroed in the SCF program
c
      common /counters/nintsum1,nintsum2,nquartsum
      common /job/jobname,lenJ
cccccc
      common/ftctime/tinit,einit,tsmooth,esmooth,tftc1,eftc1,
     $               tmix,emix,tftc2,eftc2,tfock1,efock1,tmpole,empole
cccccc
c
c  character variables for COSMO
c
      character*100 cerm
c
      PARAMETER (IUnit=1)          ! unit for checkpoint read
c  ...................................................................
      data nrad/400/, nang/434/, IradQ/0/, lrad/0/, lang/0/
      parameter (MaxEL=94,MaxC=5)  ! DFT dispersion
      common /dftpnt/ ird,iaij,idst,iixx,iiwt,ipre,iexp,
     $                ixgd,iwt,iscr,icofm,nscr
c  ...................................................................
      data opnames / 'gapm','nodd','prin','sthr','thre','diis',
     1               'ipol','punc','iter','lvsh','virt','pseu',
     2               'nfoc','loca','star','gran','dftp','grid',
     3               'core','dft ','nbat','semi','aufb','pwav',
     4               'pwgd','rang','rest','lrad','lang','nomp',
     5               'expc','dmul','ccdd','fact','elpo','locv',
     6               'disp','nabc','tz  ' /
c
      data ioptyp /11, 1, 1,11,11,11,  11, 1, 1,11, 1,11,
     $              1,21, 1, 1,21,11,  11,21, 1, 0, 1, 0,
     $             12, 2, 0, 1, 1, 0,   1,11, 0,11, 0,21,
     $              1, 0, 0/
c
      data dftxc/'hfs   ','svwn  ','svwn5 ','hfb   ','bvwn  ','bvwn5 ',
     1           'bp86  ','bpw91 ','blyp  ','bvp86 ','optx  ','ovwn  ',
     2           'ovwn5 ','op86  ','opw91 ','olyp  ','pw91  ','pbe   ',
     3           'o3lyp ','b3lyp ','b3pw91','wah   ','user  ','b97   ',
     4           'b97-1 ','b97-2 ','hcth  '/
c
c
      call elapsec(et0)
c
c ------------------------------------------------
c  set memory marker
      call mmark
c ------------------------------------------------
c  first set the defaults
      np1=1             ! file used to store integral quantities
      np4=4             ! file used to store SCF quantities
      nfock = 1         ! no multifock
      nscr = 0          ! DFT scratch memory
      nodelta = 30      ! switch off delta density after 30 cycles
      gapmin = 0.3d0    ! minimum homo-lumo gap  (for level shift)
      laufbau = 0       ! do not force aufbau occupancy
      lgrid = 1         ! force grid recalculation for DFT
      idisp = -1        ! NO dispersion in DFT unless requested
      wvfnc = 'RHF '
      escfonly = .True.
c
c ::::::::::::::::::::::::::::::::::::::::::
c -- parallel?
      call getival('nslv',nslv)
c ::::::::::::::::::::::::::::::::::::::::::
      call getrval('zero',zero)
      call getrval('half',half)
      call getrval('one',one)
      call getrval('two',two)
      call getival('iout',iout)
      call getival('isumscf',isumscf)
      call getival('ncf',ncf)
      call getival('ncs',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('ictr',ictr)
      call getrval('qtot',qtot)
      call getival('na  ',na)
      call getival('ndum',ndum)           ! number of dummy atoms
      call getival('ngener',ngen)         ! number of generators
      call getival('SymFunPr',ifp)        ! basis function symmetry pairs
      call getival('npsp',npsp)           ! number of pseudopotentials
c ------------------------------------------------------------------
c    get Cartesian coordinates from the <coord> file
c    put into bl(inuc) [XNuc] and bl(ibas) [BASDAT] arrays
      CALL GetCoord(na,nsh,ncs,bl(ictr),bl(inuc),bl(ibas))
c
c    if part of a geometry optimization, the initial integral
c    threshold is dependent on the rms gradient
      CALL GetIntThrsh
c ------------------------------------------------------------------
c  total charge
      call getrval('charge',charge)
      iprscf=0                            ! print level
      IVirt=5                             ! # virtual MOs to print
      sthre=6.0d0                         ! linear-dependency threshold
      scfthre=5.0d0                       ! SCF convergence threshold
      maxiter=50                          ! maximum iterations
      diisthr=two                         ! diis parameter
c   default point for switching to pseudo-diagonalization
c   on smaller molecules, set pseudo to zero
      pseudo=zero
      if(ncf.gt.250) pseudo=0.002d0
c   this value governs the Jacobi rotations. Only the rotations which
c   are larger than cutoff*(max. rot. angle) are performed.
c   maybe only the rotations in the same row and column should be
c   considered
c   cutoff for orbital rotation in pseudo-diagonalization
      cutoff=0.001d0
      nstart=0
      igran=20                            ! integral granularity
      xlvsh = one                         ! level shift
c dft type: 0=HF, 1=HFS, 2=SVWN, 3=SVWN5, 4=HFB, 5=BVWN, 6=BVWN5,
c           7=BP86, 8=BPW91, 9=BLYP, 10=BVP86, 11=OPTX, 12=OVWN,
c           13=OVWN5, 14=OP86, 15=OPW91, 16=OLYP, 17=PW91, 18=PBE,
c           19=O3LYP, 20=B3LYP, 21=B3PW91,22=WAH, 23=USER (user-defined)
c           24=B97, 25=B97-1, 26=B97-2, 27=HCTH
      idft=0
c  scale factor for dft grid quality (principally radial grid)
      fac0=1.25d0
c  default core (limiting energy separating core from valence)
c   (only used in localization procedure)
      core = -3.0d0
      NBatch = 50                         ! # points per batch in DFT
      lsemi = 0                           ! fully direct is default
      lpwgrid = 0                         ! initialization of plane-wave
      edisp = zero                        ! distersion energy
c ..................................................................
c  input options
c
      call readopt(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
      if(ifound(1).gt.0) gapmin=ropv(1,1)
      if(ifound(2).gt.0) nodelta=iopv(1,2)
      if(ifound(3).gt.0) iprscf=iopv(1,3)
      if(ifound(4).gt.0) sthre=ropv(1,4)
      if(ifound(5).gt.0) scfthre=ropv(1,5)
      if(ifound(6).gt.0) diisthr=ropv(1,6)
      if(ifound(9).gt.0) maxiter=iopv(1,9)
      if(ifound(10).gt.0) xlvsh=ropv(1,10)
      if(ifound(11).gt.0) ivirt=iopv(1,11)
      if(ifound(12).gt.0) pseudo=ropv(1,12)
c -- orbital localization
      locali=.false.
      if(ifound(14).gt.0) then
        locali=.true.
        lmeth=chopv(14)(1:5)
        call lowerca2(lmeth,5)
      end if
c -- virtual orbital localization
      localv=.false.
      if(ifound(36).gt.0) then
        localv=.true.
        lmethv=chopv(36)(1:5)
        call lowerca2(lmethv,5)
      end if
c -- if MP2 is next in the input file, force localization
      If(locali.eq.false) Then
        READ(inp,'(A4)',END=94) chopv(14)(1:4)
        call lowercas(chopv(14),4)
        If(chopv(14)(1:3).EQ.'mp2') Then
          locali=.true.
          lmeth='pipek'
          write(iout,*) ' == Pipek Localization forced for MP2 Job =='
        EndIf
        BACKSPACE inp
      EndIf
 94   CONTINUE
c -------------------------
      if(ifound(16).gt.0) igran=iopv(1,16)
      call setival('gran',igran)
      if(ifound(22).gt.0) lsemi=-1
c ...................................................................
c -- DFT
      If(ifound(17).gt.0.or.ifound(20).gt.0) then
c -- check basis set (F is highest angular momentum basis function coded)
        call ChkAngMom(ncs,bl(ictr),1)
        if(ifound(17).gt.0) dftype=chopv(17)(1:6)
        if(ifound(20).gt.0) dftype=chopv(20)(1:6)
        if(ifound(37).gt.0) idisp=iopv(1,37)      ! dispersion term
        noabc = ifound(38).gt.0                   ! omit 3-body dispersion term
        tz = ifound(39).gt.0                      ! special parametrization
        call lowerca2(dftype,6)
        wvfnc = 'RDFT'
        do i=1,ndft
          if(dftype.eq.dftxc(i)) then
            idft=i
            exit
          end if
          idft=ndft+1
        end do
      End if
      If(idft.gt.0) then
        if(idft.gt.ndft) then
         call nerror(9,'scfmain','undefined dft XC potential',idft,idft)
        else
         write(iout,*)  'DFT exchange-correlation potential=',dftype
         if(idisp.ge.0) write(iout,*)
     $                  'Dispersion correction will be included'
         if(idisp.ge.50) call RdDispUser(inp,values)
        end if
        call setival('dft',idft)
c -- set DFT exchange-correlation coefficient common block
        call setDFTblock(inp,iout,idft,ax)
      End if
c ...................................................................
      if(ifound(18).gt.0) fac0=ropv(1,18)
      if(ifound(34).gt.0) fac0=ropv(1,34)   ! backwards compatability
      if(ifound(19).gt.0) core=ropv(1,19)
      call setrval('core',core)
      if(ifound(21).gt.0) NBatch=iopv(1,21)
      if(ifound(23).gt.0) then      ! enforce occupancy of lowest MOs
        laufbau=iopv(1,23)
        if(laufbau.eq.0) laufbau=1
      endif
c ...................................................................
c -- plane wave options
c
      lpwave = ifound(24).gt.0      ! set plane wave flag
      call setival('ftc0',0)        ! switch off FTC (may be set later)
c
      If(lpwave) Then
c -- can we do plane wave?
       If(idft.eq.0.OR.aX.ne.zero) Then
        call nerror(10,'scfmain',
     $         'Cannot use FTC algorithm with this functional',0,0)
       else
        write(6,*) 'Using FTC plane wave algorithm'
       endif
      Else
        write(6,*) 'Using traditional integral algorithm'
      EndIf
c
      If(lpwave) Then
c -- check for segmented basis
        call check4lgc(ncs,bl(ictr),lshells,lgcshell,ltype_gc,ldeep_gc)
        if(lgcshell.gt.0) then
          call nerror(13,'scfmain',
     $          'This Basis MUST be segmented for FTC energy',0,0)
        endif
c
        griddensf=3.75d0
        griddensi=2.75d0
        icorecut = 1
        call setival('iroute',2)      ! force route=2 for integrals
        iexpl = 6                     ! exponential accuracy limit
        if(ifound(25).gt.0) then      ! grid density (final+initial)
          griddensf=ropv(1,25)
          if(ropv(2,25).gt.0.1d0) griddensi=ropv(2,25)
        endif
        irange1 = 7
        irange2 = 7
        irange1act=6
        irange2act=6
        if(ifound(26).gt.0) then      ! ranges (1 -> final and 2 -> start)
          irange1=iopv(1,26)
          irange2=irange1
          irange1act=iopv(2,26)
          irange2act=irange1act
        endif
cc        if(ifound(30) .gt. 0) then
cc          nomultipole=1
cc          write(6,*) 'multipoles will not be used in FTC algorithm'
cc        else
cc          nomultipole=0
cc          dist2mp = 6.0d0
cc          if(ifound(32).gt.0) dist2mp=ropv(1,32)
cc          If(dist2mp.LT.6.0d0) Write(iout,1150)
cc 1150  Format('**WARNING** Multipole Distance Cutoff ',f9.5,
cc     $        ' too low',/'   Accuracy will be compromised')
cc        end if
c - temporary new option 'ccdd' - include/exclude (cc|dd) integrals
c   in FTC space
cc        iccdd = 1
cc        if(ifound(33).gt.0) then
cc           iccdd = 0
cc           write(6,*) ' (cc|dd) integrals evaluated in plane wave space'
cc        endif
        nomultipole = 1   ! WARNING - CURRENTLY DISABLED (JB)
        call setival('expl',iexpl)
        call setival('nomultipole',nomultipole)
        call setival('ccdd',1)   ! force (cc!dd) integrals in classical space
      EndIf
c
c -- determine logical variable escfonly
c    normally FORCE is the next line for an SCF gradient and
c    we are going to check for this automatically here
        READ(inp,'(A4)',END=95) char(1:4)
        call lowercas(char(1:4),4)
        If(char(1:4).EQ.'forc') Then
          escfonly=.FALSE.
          call setival('expl',iexpl)
          call setival('irange',irange1)
          call setrval('griddf',griddensf)
        EndIf
        BACKSPACE inp
 95     CONTINUE
c ...................................................................
c
      if(ifound(28).gt.0) lrad=iopv(1,27)
      if(ifound(29).gt.0) lang=iopv(1,28)
      if(ifound(31).gt.0) iexpl=iopv(1,31)
c
c -- see if field defined earlier (e.g. in polarizability module)
      call getival('field',ifield)
      If(ifield.ne.0) then
        call getrval('fieldX',field(1))
        call getrval('fieldY',field(2))
        call getrval('fieldZ',field(3))
        call zeroit(dip,3)
        call nucdipole(na,dip,bl(inuc))
        enucF = -ddot(3,dip,1,field,1)   ! nuclear contribution in field
      else
        field(1)=zero
        field(2)=zero
        field(3)=zero
        enucF = zero
      end if
c
c    check if the COSMO flag has been defined, otherwise
c    define it now.
c
      call tstival('cosmo',icosmo)
      If(icosmo.EQ.0) Then
        call setival('cosmo',0)
      Else
        call getival('cosmo',icosmo)
      EndIf
c
      natom = na-ndum                        ! # real atoms
c .................................................................
c
c -- initialize timing
      call secund(t0)
      oneltime=zero
      twoeltime=zero
      timeloca=zero
      totcosmo=zero
      totdiag=zero
      totdiis=zero
      totftc=zero
c  nuclear repulsion energy
      call calcnuc(na,bl(inuc),enuc)
      call setrval('enuc',enuc)
      ntri=ncf*(ncf+1)/2
c  total number of electrons
      nel=qtot-charge+0.0001d0
      nmo=nel/2
c  nmo1 is the number of virtuals orbitals calculated
c  (except for the last cycle, this should be nmo+1)
c  **WARNING**  Do NOT set this TO ANY OTHER VALUE
      nmo1 = nmo+1
      if(nmo1.gt.ncf) nmo1=ncf
      call setival('nmo',nmo)
c   convert ph-like quantities to final numbers
      scfthre=10.0d0**(-scfthre)
      sthre=10.0d0**(-sthre)
c
c .................................................................
c  allocate memory
c
      call matdef('hmat','s',ncf,ncf)
      call matdef('diag','d',ncf,ncf)
      call matdef('coef','q',ncf,ncf)
      call matdef('dens','s',ncf,ncf)
      call matdef('fock','s',ncf,ncf)
c
c   differential Fock matrix builder
      call matdef('oldfock','s',ncf,ncf)
      call matdef('olddens','s',ncf,ncf)
c
      call matdef('u','q',ncf,ncf)
      call matdef('uinv','q',ncf,ncf)
      call matdef('commut','a',ncf,ncf)
c
      ihmat = mataddr('hmat')                 ! pointer to H0 (1-el)
      icof = mataddr('coef')                  ! pointer to MO coeffs
      iadm = mataddr('dens')                  ! pointer to density
      iafm = mataddr('fock')                  ! pointer to Fock matrix
      iod = mataddr('olddens')                ! pointer to old density
      iof = mataddr('oldfock')                ! pointer to old alpha Fock
      icom = mataddr('commut')                ! pointer to commutator
      idiag = mataddr('diag')                 ! pointer to orbital energies
c  ....................................................................
c   initialize DFT:1

      totdft = zero                           ! total dft cpu time
      IF(idft.NE.0) THEN
        call elapsec(tdft1)
c
c -- use differential exchange-correlation matrix
        If(lsemi.EQ.-1) Then
          call matdef('oldxc','s',ncf,ncf)
          ioxc = mataddr('oldxc')             ! pointer to old XC matrix
        EndIf
c
c -- restore coordinates, atomic charges and get symmetry data
        call mmark
        call getmem(3*natom,ixnc)             ! nuclear coordinates
        call getmem(natom,ian)                ! atomic charges/numbers
        call getmem(natom,iuq)                ! symmetry-unique atoms
        call getmem(ncf,idmx)                 ! max. density/column
c -- **WARNING** use old density matrix storage for XC matrix
        ixc = iod                             ! XC matrix
c
        call getnucdat(natom,bl(ixnc),bl(ian))
c
c -- prepare symmetry data (done even if no symmetry)
c -- be careful as DFT only needs real atoms
        nupr = 1                              ! not set if symmetry not used
        If(nsym.gt.0) call getival('SymNuPr1',nupr)
c
        call GetUNQ(natom,  nsym,   bl(ixnc),bl(nupr),NQ,
     $              bl(iuq),bl(icof))         ! icof used as scratch
c  ....................................................................
c -- if requested add Grimme's dispersion correction term
        If(idisp.GE.0) Then
          call mmark
c -- need to remove any ghost atoms
          call getmem(3*natom,ixnc1)
          call getmem(natom,ian1)
          call getnucdat(natom,bl(ixnc1),bl(ian1))
          call getmem(natom,ian2)
          call delete_ghosts(natom,bl(inuc),latom,bl(ixnc1),bl(ian2))
          call retmem(1)
CPP
c         call print_atomic_numbers(natom,iout,bl(ian1))
CPP
          call getmem(MaxEL**2,ir0ab)
          call getmem(MaxEL*MaxEL*MaxC*MaxC*3,ic6ab)
          call getmem(latom,icn)
          call getmem(MaxEL,imxc)
          call getmem(MaxEL,itmp)
          call dftd3(idft,     idisp,   iprscf,   latom,    bl(ian1),
     $               bl(ixnc1),values,  MaxEL,    MaxC,     bl(ir0ab),
     $               bl(ic6ab),bl(icn), bl(imxc), bl(itmp), noabc,
     $               tz,       edisp)
          call retmark
        EndIf
c  ....................................................................
c :::::::::::::::::::::::::::::::::::::::::::::
c -- parallel?
        If(nslv.eq.0) Then
c  -  set up DFT pointers and arrays
          call setup_dft(idft,   0,    nfock,  nrad,   nang,
     $                   NBatch, bl(ian),bl(ixnc))
          islv = 1
          lslv = 1
          ifni = 1
        Else
c -- set up arrays for controlled parallelism
          call getmem(na,islv)                ! atom/slave array
          call getmem(na,lslv)                !  ditto logical
          call getmem(na,ifni)                !  ditto premature
        EndIf
c :::::::::::::::::::::::::::::::::::::::::::::
        call elapsec(tdft2)
        tdft = tdft2-tdft1
        totdft = totdft + tdft
      ENDIf
c  .....................................................................
c
c   initialize COSMO
c
      if(icosmo.ne.0)then
        write(iout,*) ' COSMO Solvation Model Calculation'
        cosmotime=zero
        cosmotone=zero
        call secund(ct0)
        call elapsec(cet0)
c
c   compute maximum number of surface elements
c   and initialize COSMO radii.
c
        call mmark
        if(idft.EQ.0) then                      ! already done if DFT
          call getmem(3*natom,ixnc)             ! nuclear coordinates
          call getmem(natom,ian)                ! atomic charges/numbers
          call getnucdat(natom,bl(ixnc),bl(ian))
        endif 
        call getmem(natom,israd)
        call setival('c_israd',israd)
        call getint(natom,icharge)
        call setival('c_inuc',icharge)
        call cosmo_radii(natom,bl(ian),bl(israd),bl(icharge),
     $                   maxnps)
c
c   these arrays are needed for the whole COSMO calculation
c
        call getmem(3*2*maxnps,icosurf)
        call setival('c_icosur',icosurf)
        call getmem(2*maxnps,iar)
        call setival('c_iar',iar)
        call getmem(2*maxnps,iiatsp)
        call setival('c_iiatsp',iiatsp)
c
c   scratch memory for surface construction and A matrix setup
c
        call mmark
        call getival('c_nspa',nspa)
        call getival('c_nsph',nsph)
        call getival('c_nppa',nppa)
        call getmem(3*nspa,idirsm)
        call getmem(3*nsph,idirsmh)
        call getmem(3*nppa,idirvec)
        call getmem(3*nppa,idirtm)
        call getmem(3*natom,ixyzpert)
        call getmem(natom,iradtmp)
        call getmem(9*natom,itm)
        call getmem(2*maxnps,inar)
        call getmem(2*maxnps,insetf)
        call getint(natom*nppa,inset)
c
c   surface construction
c
        call cosmo_surf(natom,maxnps,nspa,nsph,nppa,bl(ixnc),
     $                  bl(icharge),bl(israd),bl(icosurf),
     $                  bl(iar),bl(iiatsp),bl(idirsm),bl(idirsmh),
     $                  bl(idirvec),bl(idirtm),bl(ixyzpert),bl(iradtmp),
     $                  bl(itm),bl(inar),bl(insetf),bl(inset))
c
c   temporary location for A matrix
c
        call getival('c_nps',nps)
        call getmem(nps*(nps+1)/2,ia1tmp)
c
c   A Matrix setup
c
        call getrval('c_disex2',disex2)
        call getrval('c_rsolv',rsolv)
        call getrval('c_routf',routf)
        call getival('c_npsphe',npspher)
        call getival('c_npsd',npsd)
        call getmem(npsd,ia23mat)
        ierr=0
        call setamat(bl(ixyzpert),natom,bl(iiatsp),
     $               bl(idirvec),bl(inar),bl(iar),bl(ia1tmp),
     $               bl(insetf),bl(inset),bl(icosurf),disex2,
     $               rsolv,routf,nppa,nps,npspher,bl(iradtmp),bl(itm),
     $               ierr,cerm,npsd,bl(ia23mat),jobname,lenj)
        call retmem(1)
cc
        If(ierr.NE.0) Then
          write(iout,'(a)') cerm
          call nerror(11,'SCF Module',
     $                   'Error in COSMO A-matrix initialization',0,0)
        EndIf
c
c   Temporarily save the A matrix to disk
c
        call cosmo_dump(jobname(1:lenj)//'.a1tmp',lenj+6,
     $                  bl(ia1tmp),nps*(nps+1)/2)
c
c   release scratch memory
c
        call retmark
c
c   memory storage for A matrix, potential and charges
c
        call getmem(nps*(nps+1)/2,ia1mat)
        call setival('c_ia1mat',ia1mat)
        call getmem(nps,iphi)
        call setival('c_iphi',iphi)
        call getmem(nps,iphin)
        call setival('c_iphin',iphin)
        call getmem(nps,iqcos)
        call setival('c_iqcos',iqcos)
c
c    memory for electron-surface repulsion for one surface
c    element
c
        call getmem(ntri,ivimat)
c
c    memory for cosmo contribution to one-electron Hamiltonian
c
        call getmem(ntri,ih0cos)
c
c   restore the A matrix
c
        call cosmo_rest(jobname(1:lenj)//'.a1tmp',lenj+6,
     $                  bl(ia1mat),nps*(nps+1)/2)
c
c   compute nuclear part of surface potential
c
        call cosmo_potn(bl(ixnc),bl(icharge),bl(iphin),
     $                  bl(icosurf),natom,nps)
        call secund(ct2)
        call elapsec(cet1)
c
        cosmotime=cosmotime+ct2-ct0
        totcosmo=totcosmo+cet1-cet0
      endif
c  .....................................................................
c
c -- now get general scratch storage
c -- this will be used in DIIS (for F, FDS-SDF and D) and for
c -- temporary and scratch storage in both SCF and DFT
      mscr = MAX(nscr+nmo*ncf,3*ntri)
      call getmem(mscr,istor)
c
c  parameters: itype: 0=overlap, 1=h0 mtx, see the rest in the routine
c  na=number of atoms,integrals,contraction info,2 directions for
c  x,y,z-dependent quantities,basis set info,nucl..info, numer of
c  contracted shells
c  determine the overlap and h0 matrices
c
      call elapsec(t1)
      call matredef('olddens','smat','s',ncf,ncf)
      call inton(0,na,bl(iod),bl(ictr),0,0,bl(ibas),bl(inuc),ncs)
c
c  write out the S matrix
      call matwrite('smat',np1,1,'s matrix')
c
c  If Cholesky decomposition is used, calculate the LU factorization of S
c  and find the inverse of the upper triangular matrix U
c  Replace cholesky with <seigvec> if the lowest eigenvalue of the overlap
c  matrix is less than sthre
      call lowesteigv('smat',ncf,xlow)
      write(iout,*) 'Lowest eigenvalue of the overlap matrix =', xlow
      xlow1=xlow
      if(xlow.gt.sthre) then
        call cholesky1('smat',ncf,'u')
        newn=ncf
        call u_inverse('u','uinv',ncf)
      else
        call seigvec('smat',ncf,sthre,'u','uinv',xlow1,newn)
        write(iout,*) 'Lowest eigenvalue of S that exceeds ',sthre,
     1                ' is ',xlow1
      end if
c
      call setival('nonredun',newn)
c
c  number of MOs to save in MOS file
      nwmo = MIN(nmo+20,newn)
c
      if(ncf.gt.newn) then
        write(iout,*)
     1     'There are ',ncf-newn,' suppressed basis functions'
c -- redefine the dimension of the MO coefficients
        write(iout,*) 'Number of non-redundant basis functions=',newn
        call matredef('commut','commut','a',newn,newn)
        call matredef('diag','diag','d',newn,newn)
        call matredef('coef','coef','r',ncf,newn)
      end if
c
c -- save xlow for nmr and MP2
c   (If <seigvec> is used, the lowest eigenvalue of S is, by definition, sthre)
      call setrval('xlows',max(xlow,sthre))
      call matredef('smat','olddens','s',ncf,ncf)
      call matwrite('olddens',np4,0,'ovlap')
c
c -- check integral thresholds
      elow = xlow
      if(elow.lt.sthre) elow = sthre
      call chkthrsh(elow,thint,thresh,iout)
c -- on return from <chkthrsh>
c      thint  - initial (loose) integral threshold
c      thresh - final   (tight) integral threshold
c
      If(thint.EQ.thresh) Then
c -- check whether to switch on semi-direct immediately for DFT
       If(lsemi.EQ.-1) lsemi=1
c -- check whether to use tight threshold immediately for FTC
       If(lpwave) Then
         griddensi = griddensf
         irange1act=irange1
         irange2act=irange2
       EndIf
      EndIf
c
      call elapsec(tt)
      oneltime=oneltime+tt-t1
c
c   put down a memory marker
      call mmark
c
      iforwhat=1
      if(lpwave) iforwhat=11
c
c ::::::::::::::::::::::::::::::::::::::::::
c -- parallel?
      If(nslv.GT.0) call para_JobInit(iforwhat)
c ::::::::::::::::::::::::::::::::::::::::::
c
c
c---------------------------------------------------
c  Initialize plane wave
      If(lpwave) Then
       griddens = griddensi
       call elapsec(tt)
       call Para_preinit_FTC(
     &        iexpl,     dist2mp,   isharpgrd, isharpness, ncspl,
     &        ncfpl,    iicsplsize, icssharps, iicsdiff,   ilistsd,
     &       iicspltype, iplbasdat, imultipolecut,nomultipole)
       call elapsec(t2)
       totftc=totftc+t2-tt
c
       tinit = zero
       einit = zero
       tsmooth = zero
       esmooth = zero
       tftc1 = zero
       eftc1 = zero
       tmix = zero
       emix = zero
       tftc2 = zero
       eftc2 = zero
       tfock1 = zero
       efock1 = zero
       tmpole=zero
       empole=zero
c
      EndIf
c-----------------------------------------
c ** TWO-ELECTRON INITIALIZATION **
c -- moved to here as 1-electron integrals are now parallel
c -- and the parallel initialization is done in <ptwoint>
c
c  iforwhat shows the type of integrals to be calculated;
c  1) for ordinary two-el.integrals (iforwhat=1)
c  2) for GIAO two-el. derivatives  (iforwhat=2)
c  3) for gradient derivatives      (iforwhat=3)
c  4) for second derivatives        (iforwhat=4)
c  6) for FTC "core" integrals      (iforwhat=11)
c  (used only in blocking routines)
c
      call ptwoint(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             lsemi,  scftype,xintxx, nblocks,maxbuffer,
     *             maxlabels)
c
      if(scftype.ne.'full-direct')then
        call getival('incore',incore)
        if(incore.gt.0)then
          write(6,*)'Integral calculation is ',scftype,
     $              ' with core storage'
        else
          write(6,*)'Integral calculation is ',scftype,
     $              ' with disk storage'
        endif
      endif
c
      call elapsec(t1)
      twoeltime=twoeltime+t1-tt
c---------------------------------------------------
c
c Construct the modified H0 matrix - needed for the very first SCF iter.
c the type of it is 10 (so you call Inton with 10-1=9 )
c
      call inton(9,na,bl(iafm),bl(ictr),0,0,bl(ibas),
     1           bl(inuc),ncs)
c
c Construct the full H0 matrix
      call elapsec(tnuc0)
      call para_oneint(1,na,bl(ihmat),bl(ictr),0,0,bl(ibas),
     1                 bl(inuc),ncs)
      call elapsec(tnuc1)
cc      write(6,*) ' one-electron time is:',tnuc1-tnuc0
c
c  pseudopotential contribution to one electron hamiltonian
c  use 'oldfock' as temporary matrix
c
      if(npsp.ne.0)then
cc        call secund(pspt0)
cc        call elapsec(pspet0)
        call getrval('ithr',psptol)
        psptol=psptol*0.01d0       ! psp thres.=0.01*main int. thres
        call matredef('oldfock','hpsp','s',ncf,ncf)
        call psph0(na,npsp,ncf,ncs,bl(ictr),bl(inuc),bl(ibas),
     $             psptol,bl(iof))
        if(iprscf.gt.3) call matprint('hpsp',iout)
        call matadd('hpsp','hmat')
        call matredef('hpsp','oldfock','s',ncf,ncf)
cc        call secund(pspt1)
cc        call elapsec(pspet1)
cc        psptime=pspt1-pspt0
cc        pspelap=pspet1-pspet0
      endif
c
c Add electric field to H0 if needed. Use 'fock' as temporary matrix
      if(ifield.ne.0) then
        do icomp=1,3
          write(iout,53) icomp,field(icomp)
 53   Format(' The ',I2,' component of electric field is added:',F10.6)
          if(abs(field(icomp)).gt.1.0d-6) then
            call inton(3,na,bl(iafm),bl(ictr),icomp,0,bl(ibas),
     1                 bl(inuc),ncs)
            call matadd1('fock',field(icomp),'hmat')
          end if
        end do
      end if
c
c  Old Penalty function code now redundant     ! PP  Sep 2008
      alpha=zero
c
c ................. penalty function ........................
c  Calculate S**-2 for supressing near-redundant linear combinations
c  S**(-2) is calculated as Uinv*Univ**(-1)*Uinv*Uinv**(-1)
c  determine the coefficient of S**(-2) from the lowest eigenvalue of S
cc      if(xlow.ge.sthre) then
cc        alpha=zero
cc      else
cc        alpha=sthre**2*1.0d-8
cc      write(iout,*) 'Penalty function is used to remove redundant ',
cc     1 'basis function linear combinations'
cc      write(iout,*) 'alpha=',alpha
c  Use 'u' and 'fock' as temporary matrices
cc        call matwrite('u',np1,0,'u mtrx')
cc        iumx=mataddr('u')
cc        iuinv=mataddr('uinv')
cc        call matcopy('uinv','u')
cc        call dtrmm('r','U','t','n',ncf,ncf,alpha,bl(iuinv),ncf,
cc     1              bl(iumx),ncf)
cc        call dtrmm('r','U','n','n',ncf,ncf,one,bl(iuinv),ncf,
cc     1              bl(iumx),ncf)
cc        call dtrmm('r','U','t','n',ncf,ncf,one,bl(iuinv),ncf,
cc     1              bl(iumx),ncf)
cc        call matcopy('u','fock')
cc        call matread('u',np1,'u mtrx')
c  write alpha*S^-2 to the disk and add it to H0
cc        call matwrite('fock',np1,0,'s**-2')
cc        call matadd('fock','hmat')
cc      end if
c ..................end penalty function ..........................
c
c -- write to binary file for possible NBO analysis
      call matwrite('hmat',np1,0,'h0-mtx ')
c one-electron timing
      call elapsec(tt)
      oneltime=oneltime+tt-t1
c
      call matzero('oldfock')
      call matzero('olddens')
c =====================================================================
c -- get the initial orbitals
c -- should be available from the GUESS module on the binary file
c -- <jobname>.mos
      call tstchval('mos-file',iyes)
      If(iyes.eq.1) Then
        call getchval('mos-file',MOS)
      Else
        MOS = jobname(1:lenJ)//'.mos'
      EndIf
      call rmblan2(MOS,256,lenM)
      itype = 1
      If(GUES(1:4).EQ.'atom') itype=3
      call ReadMOS(ncf,bl(icof),jnk,.False.,lenM,MOS,itype,IErr)
      If(IErr.EQ.0) then
        nstart = 1
c -- switch off level shift if HF and guess MOs available
        If(idft.eq.0.and.ifound(10).eq.0) xlvsh = zero
        If(itype.EQ.3) Then
c -- atomic density guess
          call tri(bl(icof),bl(iadm),ncf)
        Else
          call densma('coef','dens',nmo,.true.)
        EndIf
        call DensSymm(ngen,ncf,bl(iadm),bl(ifp),devmx)
      Else
        call elapsec(tdiag1)
c -- check the symmetry of the Fock matrix before diagonalization
        call DensSymm(ngen,ncf,bl(iafm),bl(ifp),devmx)
        If(devmx.gt.scfthre) then
          char = '  Fock matrix breaks symmetry - maximum deviation: '
          write(char(52:62),'(d10.4)') devmx
          call message('**WARNING** from SCF module',char,0,0)
        endif
c ---------------------------------------------------------------------
        If(newn.eq.ncf) Then
          call geneig('fock','u','uinv','diag','coef','olddens',zero,
     *                 nmo,'U')
        Else
          call geneig1('fock','u','uinv','diag','coef','olddens',zero,
     *                 nmo,'U')
        EndIf
c ---------------------------------------------------------------------
        call elapsec(tdiag2)
        totdiag=totdiag+tdiag2-tdiag1
        call densma('coef','dens',nmo,.true.)
        call DensSymm(ngen,ncf,bl(iadm),bl(ifp),devmx)
      End if
c -- print warning if symmetry deviation in initial density matrix
      If(devmx.gt.scfthre) then
        char =
     1  '  Initial density matrix breaks symmetry - maximum deviation: '
        write(char(62:72),'(d10.4)') devmx
        call message('**WARNING** from SCF module',char,0,0)
      endif
c
c -- set minimum level shift
      xlvsh0 = 0.3d0*xlvsh
c
      call elapsec(ett)
c
      write(iout,50) charge,nel,iprscf,sthre,scfthre,diisthr,
     $               xlvsh,xlvsh0,ett/60.0d0
 50   format(' SCF parameters:',/,
     1 ' wave function type              =',2x,'rhf',/,
     2 ' charge                          =',f8.3,/,
     3 ' number of electrons             =',i5,/,
     4 ' print level                     =',i5,/,
     5 ' threshold for linear dependency =',e10.3,/,
     6 ' SCF threshold                   =',e10.3,/,
     7 ' diis switch-on                  =',f10.5,/,
     8 ' initial level shift             =',f8.3,/,
     9 ' minimum level shift             =',f8.3,/,
     * ' Elapsed time before SCF (min)   =',f10.2)
c
      if(laufbau.gt.0) write(iout,51) laufbau
 51   format(' aufbau occupancy enforced for first ',i2,' iterations')
      if(ifield.gt.0) write(iout,52) field
 52   format(' electric field                  =',3f10.6,/)
c ======================================================================
c
c  reserve memory for labels in two-electron evaluation
c
      call getmem(maxlabels,labels)
c
c---------------------------------------------------
c Output from the call to <ptwoint>:
c nblocks   - total number of blocks of c.s.q.
c maxbuffer - maximum size of the integral buffer
c maxlabels - maximum size of 3 integer arrays for inegral labels :
c             lsize(4), lindex(4,ngcd,nbls),lgenc(nbls)
c             lsize(1)=ics_len, lsize(2)=jcs_len etc.
c             lindex(1,iqu,ijkl)=icf_start,
c             lindex(2,iqu,ijkl)=jcf_start etc.
c scftype   - scf mode to be run
c xintxx(9)  contains number of integrals more expensive than :
c 0.0%, 0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75% of maximum price
c---------------------------------------------------
c GENERAL REMARKS about run-mode:
c  (non-, semi- or full-direct)
c
c 1. the user decides which mode will be performed
c 2. it may be changed by the program only if requested
c    mode can not be performed , for instance , requested
c    mode is non -direct but there is not enough room (in
c    core or on disk) to store all integrals. In such
c    a case program is expected to switch the mode to
c    semi- or even full-direct.
c 3. if there are integrals which suppose to be stored
c    the program tries to put them in-core storage first.
c    The in-core storage will be used only if at least
c    10% of all integrals can be kept there. If not then
c    disk storage will be used.
c 4. All information above mode and storage are set up
c    in twoint.
c
c
      iter=0
      if(nstart.gt.0) iter=1
      iteration=0
c  iter2 counts the iterations following the switch to sharp integral
c  threshold
      iter2=0
      if(thint.eq.thresh) iter2=1
      conver=.false.
      conv=.false.
      itlimit=.false.
      lpseudo=.false.
      scferr=1.0d4
c  header for iteration printing
      write(iout,101)
 101  format(67x,'timing/min',/,
     1 'SCFiter    etot           e2     Brillouin',
     2 '  Delta-dens  Errsq   master  elapsed')
c
c  initialize DIIS  -- arrays do not matter
      call elapsec(tdiis1)
      ii = 1
      call diis(1,      bl(ii), bl(ii), bl(ii), bl(ii),
     1          bl(ii), bl(ii), bl(ii), ntri, .true.,
     2          .true., iprscf, xlam,   ndiis)
      call elapsec(tdiis2)
      totdiis=totdiis+tdiis2-tdiis1
c
c  Initialize DMAX
c  copy density matrix into scratch storage
      call tfer(bl(iadm),bl(istor),ntri)
      call absmax(ntri,bl(iadm),iiii,dmax)
c
c   compute and store COSMO electron-surface repulsion matrix elements
c
      if(icosmo.ne.0)then
        call elapsec(cet1)
        call secund(ct1)
        if(nslv.eq.0)then
          call cosmo_surfrep(bl(ivimat),bl(icosurf),bl(ictr),bl(ibas),
     $                       nps,ncs,ncf,ntri)
        else
          call para_cosmo_data(bl(icosurf),bl(ixnc),
     $                         natom,nps,npspher)
          call para_cosmo_surfrep(nps)
        endif
        call secund(ct2)
        call elapsec(cet2)
c
        cosmotime=cosmotime+ct2-ct1
        cosmotone=cosmotone+ct2-ct1
        totcosmo=totcosmo+cet2-cet1
      endif
c
c---------------------------------------------------
c if scftype = 'semi-direct' or 'non -direct'
c More expensive integrals will be kept in core
c OR written on a disk   !!! WARNING - needs checking with parallel  JB
c
      If(lpwave) then
        call update_run_mode(bl,bl(ictr),ncf,lsemi,scftype)
        if(scftype.ne.'full-direct')then
          call getival('incore',incore)
          if(incore.gt.0)then
            write(6,*)'Non-FTC integral calculation is ',scftype,
     $                ' with core storage'
          else
            write(6,*)'Non-FTC integral calculation is ',scftype,
     $                ' with disk storage'
          endif
        endif
      endif
c
      If(scftype.ne.'full-direct') Then
        if(lpwave) call symmoff
        call pstore(nblocks,bl,bl(ictr),thresh,bl(labels) )
        if(lpwave) call symmon
      EndIf
c---------------------------------------------------
c
      switch=.false.
      mgo=1
c
c  RETURN POINT FOR ITERATION
c
 200  continue
      call f_lush(iout)
      nintsum2=0
      nquartsum=0
      call mmark
c -- zero out the Fock matrix
      call matzero('fock')
c ..................................................................
      call elapsec(tt)
c -- switch off symmetry if FTC
      If(lpwave) call symmoff
      call pfock(idft,     ax,       nblocks,  nfock,    .true.,
     $           ncf,      bl,       bl(ictr), thint,    mgo,
     $           bl(istor),bl(iafm), rdum,     bl(iadm), rdum,
     $           bl(labels),iforwhat)
      If(lpwave) call symmon
      call elapsec(t2)
c ..................................................................
      twoeltime=twoeltime+t2-tt
      call retmark
c ..................................................................
      IF(lpwave) THEN
        If(lpwgrid.EQ.0) Then
        lpwgrid=1
        call elapsec(tt)
        call mmark
        call Para_FTC_INITIAL(
     &         griddens, irange1act, irange2act, iexpl,     ncspl,
     &         ncfpl,     isharpgrd, isharpness, iicsplsize,icssharps,
     &         iicspltype,iplbasdat, rLxmin,     rLxmax,    rLymin,
     &         rLymax,    rLzmin,    rLzmax,     rLxo,      rLyo,
     &         rLzo,      rLxe,      rLye,       rLze,      PLDmax,
     &         npwx,      npwy,      npwz,       npwxe,     npwye,
     &         npwze,     iinit,   igridranges,igridranges2,maxovs,
     &         isharpovs, insharpovs, ii4sharps, icorecut,
     &      iisharpgridrange1,iisharpgridrange2, griddensf)
        call elapsec(t2)
        totftc=totftc+t2-tt
        EndIf
        call elapsec(tt)
        call mmark
        call Para_FTC_ENERGY(
     &     griddens,   rLxmin,     rLymin,     rLzmin,     rLxo,
     &     rLyo,       rLzo,       rLxe,       rLye,       rLze,
     &     PLDmax,     npwx,       npwy,       npwz,       npwxe,
     &     npwye,      npwze,      isharpness, iinit,      ncspl,
     &     igridranges,igridranges2,icssharps, iicsdiff,   iicspltype,
     &     maxovs,     isharpovs,  insharpovs, ii4sharps,  icorecut,
     &iisharpgridrange1,iisharpgridrange2,bl(iadm),bl(ii), bl(iafm),
     &     bl(ii),     thint,      dist2mp,  imultipolecut,nomultipole,
     &     conver)
        call retmark
        call elapsec(t2)
        totftc=totftc+t2-tt
      ENDIF
c-------------------------------------kwol----------------------
      iter=iter+1
      if(iter2.gt.0) iter2=iter2+1
c  add the old Fock matrix
      call matadd('oldfock','fock')
c  Restore the full density
      if(lsemi.gt.0.and.lsemi.lt.3) call tfer(bl(iadm),bl(iof),ntri)
      call matadd('olddens','dens')
      if(lsemi.le.0.or.lsemi.eq.3) call tfer(bl(iadm),bl(iof),ntri)
c
      if(icosmo.ne.0)then
        call secund(ct0)
        call elapsec(cet1)
c
c  compute the potential on the COSMO surface
c
        if(nslv.eq.0)then
          call cosmo_pot(bl(iadm),bl(ictr),bl(ibas),bl(iphin),
     $                   bl(iphi),bl(icosurf),bl(ivimat),
     $                   nps,ncf,ncs,ntri)
        else
          call para_cosmo_pot(bl(iadm),bl(iphin),bl(iphi),nps,ntri)
        endif
        call secund(ct1)
c
c  compute the COSMO screening charges
c
        call coschol2(bl(ia1mat),bl(iqcos),bl(iphi),nps,-1.0d0)
c
c  compute COSMO contribution to one electron hamiltonian
c
        call getrval('c_fepsi',fepsi)
        if(nslv.eq.0)then
          call cosmo_h0(bl(ih0cos),bl(ivimat),bl(ictr),bl(ibas),
     $                bl(icosurf),bl(iqcos),fepsi,nps,ncf,ncs,ntri)
        else
          call para_cosmo_h0(bl(ih0cos),bl(iqcos),fepsi,nps,ntri)
        endif
c
c  the one electron COSMO Hamiltonian need to be symmetrized.
c  I call denssymm because in fact I just need to average the
c  symmetry-equivalent elements.
c  Print a warning message if things look bad.
c
        call DensSymm(ngen,ncf,bl(ih0cos),bl(ifp),devmx)
        if(devmx.gt.1.0d-4) then
          char='COSMO H0 matrix  breaks symmetry - maximum deviation: '
          write(char(60:70),'(d10.4)') devmx
          call message('**WARNING** from SCF module',char,0,0)
        end if
        call secund(ct2)
        call elapsec(cet2)
        cosmotime=cosmotime+ct2-ct0
        cosmotr=ct1-ct0
        totcosmo=totcosmo+cet2-cet1
      endif
c
c  calculate the total energy
c
      call matprodtr('fock','dens',e2)
      call matprodtr('hmat','dens',e1)
c
c  COSMO correction to one-electron energy
c
      if(icosmo.ne.0)then
        call secund(ct0)
        call elapsec(cet1)
        call cosmo_ediel(bl(iqcos),bl(iphi),fepsi,nps,ediel)
        e1=e1+ediel
c
c     store ediel in depository
c
        call setrval('c_ediel',ediel)
        call secund(ct1)
        call elapsec(cet2)
        cosmotime=cosmotime+ct1-ct0
        totcosmo=totcosmo+cet2-cet1
      endif
c ........................................................
c  subtract any penalty function from the energy
      if(alpha.gt.zero) then
        call matconn('s**-2','s',ncf,ncf,icom)
        call matread('s**-2',np1,'s**-2')
        call matprodtr('s**-2','dens',penalty)
        if(iprscf.ge.1) write(iout,*) 'penalty=',penalty
        e1=e1-penalty
        call matdisc('s**-2')
      end if
      e2 = e2*half                             ! JB
      etot=e1+e2+enuc+enucF
      call matadd('hmat','fock')
c
c  add COSMO contribution to total Fock matrix
c
      if(icosmo.ne.0) then
        call secund(ct0)
        call elapsec(cet1)
        call addvec(ntri,bl(iafm),bl(ih0cos),bl(iafm))
        call secund(ct1)
        call elapsec(cet2)
        cosmotime=cosmotime+ct1-ct0
        totcosmo=totcosmo+cet2-cet1
      endif
c  ...................................................................
      If(idft.NE.0) Then
        call elapsec(tdft1)
        thrsh = thint
c -- determine grid quality factor based on scf convergence
c    (intermediate value factor=1.0 eliminated due to diis)
        If(thint.GT.thresh) factor=0.6d0*fac0
        If(thint.LE.thresh) factor=fac0
c -- form maximum density matrix element per column
        call FormMaxDen(ncf,bl(iof),bl(idmx))
c -- zero exchange-correlation matrix
        call zeroit(bl(ixc),ntri)
c ========================================================================
        lsemit = lsemi
        if(lsemi.eq.3) lsemit=0
        CALL Para_Dft(idft,   natom, bl(ixnc), bl(ian), nsym,
     $                ngen,   NSY,  bl(nupr), NQ,    bl(iuq),
     $                iprscf, lrad,   lang,    IradQ, factor,
     $                NBatch, lgrid,  lsemit, bl(idst), bl(iaij),
     $                bl(ird),bl(iixx),bl(iiwt),bl(ixgd),bl(iwt),
     $                thrsh,  ncf,    ncs,   bl(ibas), bl(ictr),
     $                bl(ipre),bl(iexp),.true., bl(iof),bl(ii),
     $                bl(idmx),bl(islv),bl(lslv),bl(ifni),nscr,
     $                bl(istor),bl(ixc), bl(ii), edft,    el)
c ========================================================================
c -- if symmetry, need to symmetrize exchange-correlation matrix
        If(nsym.GT.0) Then
          CALL FockSymm_SCF(ngen,ncf,1,bl(ixc),bl(ifp))
c -- correct edft and el
          el = el*DFloat(nsym+1)
          edft = edft*DFloat(nsym+1)
        EndIf
        If(lsemi.EQ.2) Then
c -- form full exchange-correlation matrix from difference matrix
          call matadd('oldxc','olddens')
          call cpyvec(ntri,bl(ixc),bl(ioxc))
        Else If(lsemi.EQ.1) Then
          lsemi = 2
          call cpyvec(ntri,bl(ixc),bl(ioxc))
        EndIf
c -- add exchange-correlation matrix to Fock matrix
        call matadd('olddens','fock')
        If(iprscf.gt.1) write(iout,1100) el,edft,edisp
 1100 FORMAT(/,' .................................................',
     $       /,'    Call From DFT Routines',
     $       /,'  Number of electrons over grid:   ',F12.6,
     $       /,'  DFT Exchange-Correlation energy: ',F15.9,
     $       /,'  Dispersion energy correction:    ',F15.9,
     $       /,'  ................................................',/)
c ========================================================================
        lgrid = 0
        call elapsec(tdft2)
        tdft = tdft2-tdft1
        totdft = totdft + tdft
cc        write(iout,1234) tdft,totdft
 1234   format(1x,' dft time this cycle: ',f9.2,'  total: ',f9.2)
c -- add exchange-correlation energy to 2-el & total energy
        e2 = e2 + edft
        etot = etot + edft + edisp
      EndIf
c  ....................................................................
c
c -- check for convergence
c -- commutator of F and D stored in 'commut'
      call converged1('u','uinv','dens','fock','commut',scferr,errsq)
c -- have we converged?
      conv = scferr.lt.scfthre
c -- set diagonalization method (D=diag; P=pseudodiag)
      dmethod='P'
      if(scferr.gt.pseudo.or.conver.or.itlimit.or.iter.le.2) dmethod='D'
      if(newn.ne.ncf) dmethod='D'      ! pseudodiagonalization likely NOT working
c                                        with suppressed basis functions   JB Nov 2009
c
      call secund(cumt)
      cumt=(cumt-t0)/60.0d0
      call elapsec(ecumt)
      ecumt=(ecumt-et0)/60.0d0
      iteration = iteration+1
      write(iout,301) iteration,etot,e2,scferr,dmax,errsq,
     $                cumt,ecumt,dmethod
 301  format(i3,f16.9,f13.6,e10.3,f9.5,e10.3,2f8.2,1x,a1)
c
      IF(.not.conv.AND..not.itlimit) THEN
c
c -- save current density in <olddens> and current Fock matrix
c -- (minus 1-electron/DFT part) in <oldfock> ready for next SCF cycle
        If(idft.NE.0) Then
          do i=0,ntri-1
          bl(iof+i) = bl(iafm+i) - bl(ihmat+i) - bl(ixc+i)
          bl(iod+i) = bl(iadm+i)
          enddo
        Else
          do i=0,ntri-1
          bl(iof+i) = bl(iafm+i)  - bl(ihmat+i)
          bl(iod+i) = bl(iadm+i)
          enddo
        EndIf
        if(icosmo.ne.0)then
          do i=0,ntri-1
          bl(iof+i) = bl(iof+i)  - bl(ih0cos+i)
          enddo
        endif
c
c  perform the DIIS extrapolation if sufficiently converged
        If(scferr.lt.diisthr.and.iter.gt.1.and.
     $     iteration.gt.laufbau) Then
          call elapsec(tdiis1)
          call diis(0,      bl(icom), bl(ii),   bl(iafm),bl(ii),
     1               bl(iadm),bl(ii), bl(istor),ntri,    .true.,
     2              .true., iprscf,   xlam,     ndiis)
          call elapsec(tdiis2)
          totdiis=totdiis+tdiis2-tdiis1
        EndIf
      ENDIF
c
c .........................................
      if(conver.or.itlimit) then
         mgo=5                   ! tell slaves SCF is over, but do not
         call para_next(mgo)     ! collect timing data now
         go to 400
      else
         call para_next(mgo)
      end if
c .........................................
      if(conv) conver=.true.
c
c  if sufficiently converged, sharpen the integral threshold
c  do this anyway if there is no rough convergence by the 10th iteration
c  discontinue delta density if the iteration is almost finished
      IF( (thint.gt.thresh.and.(scferr.lt.1.0d-3.or.iter.gt.10)) .or.
     1    (iter2.gt.0.and.mod(iter2,10).eq.0.and.scferr.gt.two*scfthre)
     2    .or.(iter.gt.nodelta).or.(conv.or.itlimit) ) THEN
        iter2=1
        mgo=2
        write(iout,150) thresh
 150  Format(' Switching to Full Fock Evaluation',
     $       '  Integral Threshold: ',D12.4)
c
c  do NOT keep any old DIIS vectors if integral accuracy changed
        If(thint.gt.thresh) Then
          call elapsec(tdiis1)
          call diis(-1,     bl(ii), bl(ii), bl(ii), bl(ii),
     1              bl(ii), bl(ii), bl(ii), ntri,   .true.,
     2             .true.,  iprscf, xlam,   ndiis)
          call elapsec(tdiis2)
          totdiis=totdiis+tdiis2-tdiis1
          thint=thresh
          lgrid = 2        ! force grid recalculation + delete old grid files
          If(lpwave) Then
            griddens=griddensf
            lpwgrid=0
            irange1act=irange1
            irange2act=irange2
c  return the memory allocated in the plane wave initialization routine
            call retmark
          EndIf
        EndIf
c
c  do not continue with delta density if the threshold was changed
        switch=.true.
        If(lsemi.NE.0) lsemi = 1
      ENDIF
c
c  diagonalize the Fock matrix and add the level shift.
c  It can be shown that UDU(t) has to be added to Uinv(t)FUinv
c  use regular diagonalization at the beginning and in the
c  last step. Otherwise use pseudo-diagonalization
c  in the first step,add only 90% of the charge if started from scratch
      if(iter.eq.1) then
        call matscal('fock',0.9d0)
        call matadd1('hmat',0.1d0,'fock')
c
c  take care of a possible COSMO contribution here?
c
        if(icosmo.ne.0)then
          call secund(ct0)
          call elapsec(cet1)
          call daxpy(ntri,0.1d0,bl(ih0cos),1,bl(iafm),1)
          call elapsec(cet2)
          call secund(ct1)
          cosmotime=cosmotime+ct1-ct0
          totcosmo=totcosmo+cet2-cet1
        endif
      end if
c
c  factor of half because of double occupancy in dens
      xlv=half*xlvsh
      if(laufbau.ge.iteration) xlv=zero    ! no level shift while aufbau on
      IF(scferr.gt.pseudo.or.iter.le.2.or..not.lpseudo) THEN
        nmop = nmo1
        If(scferr.le.pseudo.or.lpseudo) Then
         nmop = 0      ! all MOs needed for pseudodiag next cycle
         lpseudo = .true.
        EndIf
        call elapsec(tdiag1)
c ---------------------------------------------------------------------
        If(newn.eq.ncf) Then
          call geneig('fock','u','uinv','diag','coef',
     $                'dens',xlv, nmop, 'U')
        Else
          call geneig1('fock','u','uinv','diag','coef',
     $                 'dens',xlv, nmop, 'U')
        EndIf
c ---------------------------------------------------------------------
        call elapsec(tdiag2)
        totdiag=totdiag+tdiag2-tdiag1
c  restore the level shifted eigenvalues to proper values
        if(xlv.ne.zero) call fixdiag('diag',nmo,xlvsh)
c -- save the MOs
        call WriteMOS(ncf,    nwmo,   bl(icof),jnk, .False.,
     $                lenM,   MOS,    1)
c  calculate the HOMO-LUMO gap
        call matelem('diag',nmo,nmo,homo)
        call matelem('diag',nmo1,nmo1,xlumo)
        gap=xlumo-homo
      ELSE
c -- pseudodiagonalization
        call pdiag('fock','coef','diag',nmo,ncf,cutoff,xlv)
c -- save the MOs
        call WriteMOS(ncf,    nwmo,   bl(icof),jnk, .False.,
     $                lenM,   MOS,    1)
      ENDIF
c
      if(iteration.ge.maxiter-1) itlimit=.true.
c
c  calculate the density matrix
      call densma('coef','dens',nmo,.true.)
c -- check symmetry of density matrices
      call DensSymm(ngen,ncf,bl(iadm),bl(ifp),devmx)
      if(devmx.gt.scfthre) then
        char = '  Density matrix breaks symmetry - maximum deviation: '
        write(char(60:70),'(d10.4)') devmx
        call message('**WARNING** from SCF module',char,0,0)
      end if
c  subtract the old density except if a switch to full density was made
      if(.not.switch.and..not.(conv.and.iter2.gt.1)) then
        call matadd1('olddens',-one,'dens')
      end if
c -- copy delta density into scratch storage
      call tfer(bl(iadm),bl(istor),ntri)
      call absmax(ntri,bl(istor),iiii,dmax)
c
      if(switch.or.(conver.and.iter2.gt.1)) then
        call matzero('oldfock')
        call matzero('olddens')
        switch=.false.
        if(conver.and.lsemi.gt.0) lsemi=3
      end if
c
      IF(xlvsh0.NE.zero) THEN
c -- set the level shift
c -- mimimum level shift is xlvsh0
       if(gap.gt.gapmin) then
          xlvsh=xlvsh0
          if(ifound(10).eq.0) xlvsh=zero
        else
          xlvsh=MAX(gapmin-gap,xlvsh0)
        end if
        if(iprscf.gt.1) write(iout,*) ' level shift is ',xlvsh
      ENDIF
c
c -- debug printout
      if(iprscf.gt.3) then
        call matprint('fock',iout)
        call matprint('diag',iout)
        call matprint('coef',iout)
        call matprint('dens',iout)
      end if
c
c --  LOOP BACK TO THE BEGINNING TO DO SCF
      go to 200
c
 400  continue
c
      If(lpwave) call retmark      ! return memory in plane wave setup
c
c -- Make the final diagonalization complete to get all virtuals
c -- if converged use big level shift to ensure retention of MO order
      If(conver) xlv=5.0d0
      call elapsec(tdiag1)
c ---------------------------------------------------------------------
      If(newn.eq.ncf) Then
        call geneig('fock','u','uinv','diag','coef','dens',xlv,0,'U')
      Else
        call geneig1('fock','u','uinv','diag','coef','dens',xlv,0,'U')
      EndIf
c ---------------------------------------------------------------------
      call elapsec(tdiag2)
      totdiag=totdiag+tdiag2-tdiag1
      if(xlv.ne.zero) call fixdiag('diag',nmo,two*xlv)
c
c -- calculate the maximum absolute value of the occupied coefficients
      call absmax(nmo*ncf,bl(icof),imax,xmax)
      if(xmax.gt.5.0d0) then
       write(iout,*) '** Warning - largest occupied MO coefficient is ',
     1                   xmax
      endif
      xmax0 = xmax
      call absmax(newn*ncf,bl(icof),imax,xmax)
      if(xmax.gt.10.0d0.AND.xmax.gt.xmax0) then
       write(iout,*) '** Warning - largest virtual MO coefficient 1s ',
     1                   xmax
      endif
c
c -- save all MOs (AND the orbital energies)
c
c first zero out the non-used SCF coefficients and put in huge numbers
c for the corresponding orbital energies
      if(ncf.gt.newn) then
        istcoef=icof+ncf*newn
        istorbe=idiag+newn
        do ii=newn+1,ncf
          call zeroit(bl(istcoef),ncf)
          istcoef=istcoef+ncf
          bl(istorbe)=1.0d15
          istorbe=istorbe+1
        end do
      end if
c
c  Symmetry-adapt the orbitals to the abelian subgroup of the full point group
c  This is only necessary if there are degenerate symmetry representations
c  (** NOTE: Does not work for very high symmetry **)
      call getchval('Sflies',Sflies)
      if(nsym.gt.0.and.Sflies.NE.'ih  ') then
c -- read the overlap matrix into 'olddens'
        call matread('olddens',np1,'s matrix')
c -- expand it into square matrix form using existing scratch storage
        CALL EXPAND(ncf,bl(iod),bl(istor))
c -- reuse existing storage for scratch arrays
        iiu1 = mataddr('u')
        iiu2 = mataddr('uinv')
        call getmem(newn**2,icsc)
        call SymORB(ngen,   ncf,    newn,   bl(ifp),  bl(icof),
     1           bl(istor),bl(iiu1),bl(iiu2),bl(icsc),bl(idiag))
        call retmem(1)
      end if
c
      call WriteMOS(ncf,    newn,   bl(icof),bl(idiag),.True.,
     $              lenM,   MOS,    1)
c
      if(itlimit.and..not.conver) then
        write(iout,*) 'Attention: NO CONVERGENCE'
        write(iout,550) Iteration,etot,scferr
  550 format(' No convergence in ',i4,' steps, Energy= ',F18.9,' Eh',/,
     $       '    Brillouin=',E12.5)
      end if
c
c  if we made it to this point then we are done. print
c  results and exit.
c
c  calculate the kinetic energy matrix for NBO analysis
c  inton with ityp=2 is the kinetic energy
c
      call inton(2,na,bl(ihmat),bl(ictr),0,0,bl(ibas),
     1           bl(inuc),ncs)
      call matwrite('hmat',np1,0,'kinetic ')
c
c  write out a few SCF quantities
      xnmo(1)=nmo
      xnmo(2)=etot
      call wri(xnmo,2,np4,1,'nocc_rhf')
      call matwrite('dens',np4,0,'den0_rhf')
      call matwrite('fock',np4,0,'fock_rhf')
      call matwrite('coef',np4,0,'evec_rhf')
      call matwrite('diag',np4,0,'eval_rhf')
c
c  COSMO outlying charge correction and output
c
      if(icosmo.ne.0)then
        call elapsec(cet1)
        call secund(ct0)
        call mmark
        call getival('c_npsphe',npspher)
        call getmem(npspher,iphio)
        call setival('c_iphio',iphio)
        call getmem(nps,iphic)
        call setival('c_iphic',iphic)
        call getmem(npspher,iqcoso)
        call getmem(nps,iqcosc)
        call setival('c_iqcosc',iqcosc)
        call getmem(npspher*nps,ia2mat)
        call getmem(npspher*(npspher+1)/2,ia3mat)
c
c   COSMO matrix elements for outer surface
c
        if(nslv.eq.0)then
          call cosmo_surfrep(bl(ivimat),bl(icosurf+3*nps),bl(ictr),
     $                     bl(ibas),npspher,ncs,ncf,ntri)
        else
          call para_cosmo_surfrep(npspher)
        endif
c
c   COSMO potential for outer surface
c
        call cosmo_potn(bl(ixnc),bl(icharge),bl(iphin),
     $                  bl(icosurf+3*nps),natom,npspher)
        if(nslv.eq.0)then
          call cosmo_pot(bl(iadm),bl(ictr),bl(ibas),bl(iphin),
     $                   bl(iphio),bl(icosurf+3*nps),bl(ivimat),
     $                   npspher,ncf,ncs,ntri)
        else
          call para_cosmo_pot(bl(iadm),bl(iphin),bl(iphio),
     $                        npspher,ntri)
        endif
c
c  compute COSMO outlying charge correction
c
        ierr=0
        call cosmo_oc(ediel,bl(iphi),bl(iqcos),bl(ia1mat),natom,
     $                bl(iphio),fepsi,nps,npspher,qsum,dq,de,bl(iqcosc),
     $                bl(iphic),ierr,cerm,bl(iqcoso),bl(ia2mat),
     $                bl(ia3mat),jobname,lenj)
cc
        if(ierr.ne.0)then
          write(iout,'(a)')'Error in COSMO outlying charge correction'
          write(iout,'(a)')cerm
          call nerror(1,'rhfmain','error in COSMO OC correction',0,0)
        endif
        call setrval('c_de',de)
        call setrval('c_dq',dq)
        call setrval('c_qsum',qsum)
c
c  write COSMO output file
c
        ierr=0
        call cosmo_write(bl,bl(icharge),bl(iiatsp),etot,bl(ixnc),
     $                   ierr,cerm,natom)
cc
        if(ierr.ne.0)then
          write(iout,'(a)')'Error in COSMO output'
          write(iout,'(a)')cerm
          call nerror(1,'rhfmain','error in COSMO output',0,0)
        endif
c
c   store on file COSMO surface data, to be used in a later part
c   of the calculation (e.g. gradient)
c
        call cosmo_store(bl(ixnc),bl(icosurf),bl(iar),
     $                   bl(iqcos),bl(iiatsp),bl(icharge),
     $                   natom,nps)
c
        call retmark
        call secund(ct1)
        call elapsec(cet2)
        cosmotime=cosmotime+ct1-ct0
        totcosmo=totcosmo+cet2-cet1
c
c  print COSMO results
c
        call cosmo_print(iout,etot,cosmotime,totcosmo)
      endif
c
c  now it is really over, we can collect timing data from slaves
c
      call para_next(0)
c
c -- start freeing up memory pointers
      call retmark
c
c  delete the DIIS files
      call elapsec(tdiis1)
      ii=mataddr('fock')
      call diis(3,      bl(ii), bl(ii), bl(ii), bl(ii),
     1          bl(ii), bl(ii), bl(ii), ntri,    .true.,
     2          .true., iprscf, xlam,   ndiis)
      call elapsec(tdiis2)
      totdiis=totdiis+tdiis2-tdiis1
c
c -- delete scratch storage
      call retmem(1)
c
      If(icosmo.ne.0) Then
c -- delete cosmo matrix elements file
        call cosmo_del('c_onel  ')
c -- deallocate memory for COSMO
        call retmark
      EndIf
c  ...................................................................
c -- now deallocate memory allocated in dft
      if(idft.ne.0) call retmark
c  ...................................................................
c  remove matrices
      if(idft.ne.0.AND.lsemi.ne.0) call matrem('oldxc')
      call matrem('commut')
      call matrem('uinv')
      call matrem('u')
      call matrem('olddens')
      call matrem('oldfock')
c
      call elapsec(tt)
      call scfresu(nmo,    0,      ncf,    ncs,    natom,
     $             etot,   e1,     e2,     enuc,   iout,
     1             ivirt,  bl(ictr),bl(ifp),ibas,iprscf, dip)
      call elapsec(t1)
      call setrval('erhf',etot)
      oneltime=oneltime+t1-tt
c.......................................................
c
      loca_done=0
      if(locali.or.localv) then
        call matdef('dipx','s',ncf,ncf)
        call matdef('dipy','s',ncf,ncf)
        call matdef('dipz','s',ncf,ncf)
        call matread('dipx',np1,'aoX     ')
        call matread('dipy',np1,'aoY     ')
        call matread('dipz',np1,'aoZ     ')
        write(6,*) 'nmo before loca',nmo
        if(locali) call loca(nmo,lmeth ,iprscf,core,'fock','coef')
c -- virtual localization (use zero for the Fermi level)
        if(localv) call loca(ncf,lmethv,iprscf,zero,'fock','coef')
        loca_done=1
c ...................................................................
c -- write the localized orbitals to the LOS file
c -- this are needed for dual-basis Local MP2    ?? - JB
        call matdef('eloc','d',nmo,nmo)
        call matread('eloc',np4,'eloc_rhf')
        iloc = mataddr('eloc')
c
        LOS = jobname(1:lenJ)//'.los'
        call WriteMOS(ncf,    newn,   bl(icof),bl(iloc),.True.,
     $                lenM,   LOS,    1)
        call matrem('eloc')
c ...................................................................
        call matrem('dipz')
        call matrem('dipy')
        call matrem('dipx')
      end if
c...................................................
      call setival('loca_done',loca_done)
      call secund(t3)
      timeloca=timeloca+t3-t1
c
c -- kw   for dual basis set mp2 & lmp2
      call savebas4mp2(bl,bl(ibas),bl(ictr),ncs,ncf,nsh,nbf)
c...................................................
c
c -- If DFT, and single node, delete all DFT grid files
      If(idft.GT.0.AND.nslv.EQ.0)
     $   Call TidyGRID(-1,NQ,bl(iuq),bl(ii),lsemi)
c
c --- clear ALL memory ---------------------
      call matreset
      call retmark
c -------------------------------------------
c
c -- reset the values of some FTC flags in order not to screw up
c -- subsequent job steps
      If(escfonly) call setival('ftc0',0)     ! switch off FTC
      call setival('nomultipole',1)           ! switch off multipoles
c
c -- write SCF energy/number of iterations to LOG file
      If(isumscf.EQ.0) CALL SumSCF(etot,iteration)
c ==========================================================
c -- write dft flag, final energy, dispersion parameters,xlow and dipole moment
c -- to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,4,'$dft',1,idft,rdum,chopv)
      If(idft.GT.0) Call wrcntrl(IUnit,7,'$factor',2,idum,factor,chopv)
      Call wrcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call wrcntrl(IUnit,7,'$energy',2,idum,etot,chopv)
      Call wrcntrl(IUnit,5,'$escf',2,idum,etot,chopv)
c -- save xlow1 in case redundant basis functions have been eliminated
      Call wrcntrl(IUnit,5,'$xlow',2,idum,xlow1,chopv)
      Call WrDIP(IUnit,dip)
      If(idisp.gt.0) Then
        Call wrcntrl(IUnit,11,'$dispersion',1,idisp,rdum,chopv)
        Call wrcntrl(IUnit,6,'$edisp',2,idum,edisp,chopv)
        Call WrDISP(IUnit,noabc,tz,values)
      EndIf
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c ==========================================================
c
c -- if semidirect, remove integral files
      If(scftype.ne.'full-direct') then
        call getival('indisk',indisk)
        if(indisk.gt.0) call clos4int()
      EndIf
c
c timings
      call secund(t1)
      t0=t1-t0
ccccccccccccccccccccccccccccccccccccccccccccccccc
c ** TEMPORARY PRINT-OUT OF FTC TIMINGS
c
cc      if(lpwave) then
cc      write(6,*) ' '
cc      write(6,3333) nslv,tinit,einit,tsmooth,esmooth,tftc1,eftc1,
cc     $              tmix,emix,tftc2,eftc2,tfock1,efock1,tmpole,empole
 3333 format (' ***** Timings for FTC steps (seconds) on ',I2,
     $        ' nodes *****'/,
     $        ' initialization    -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' smooth density    -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' FFT 1st part      -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' mixed density     -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' FFT 2nd part      -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' Smooth Fock       -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' Multipole <cd|cd> -  cpu:',f12.2,'  elapsed:',f12.2,/,
     $        ' *************************************************')
cc      endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
c----------------------------------------------------------------------
      write(iout,5001)
 5001 format(
     * /'============================================================='/
     * '                    JOB INFORMATION                          '/)
      t0=t0/60.0d0
      t1=t1/60.0d0
      oneltime=oneltime/60.0d0
      twoeltime=twoeltime/60.0d0
      totdft=totdft/60.0d0
      totftc=totftc/60.0d0
      totdiag=totdiag/60.0d0
      totdiis=totdiis/60.0d0
      timeloca=timeloca/60.0d0
      write(iout,1000) t0,t1
 1000 format(' Time for SCF and total time=',2f10.2,' min'/)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
      call elapsec(ett)
      telap=ett-et0
      telap=telap/60.0d0
      tmisc=telap-oneltime-twoeltime-totdft-totftc-totcosmo
     $           -totdiag-totdiis
      write(iout,1050) oneltime,twoeltime,totdft,totftc,totcosmo,
     $                 totdiag,totdiis,tmisc,telap
 1050 format(' Total SCF elapsed timings in minutes:',/,
     1 ' 1-el=   ',f10.2,'  2-el=   ',f10.2,'  DFT=    ',f10.2,/,
     2 ' ftc=    ',f10.2,'  cosmo=  ',f10.2,'  diag=   ',f10.2,/,
     3 ' diis=   ',f10.2,'  misc=   ',f10.2,'  total=  ',f10.2)
      write(iout,1200) nreq,nmark,lastadr-ioffset,mxmem,memtot
 1200 format(/' Memory status:'/
     *' request number=',i4,' memory marks=',i3,
     *' last used address=',i9,/,
     *' high water=',i9,' total available memory=',i9)
c----------------------------------------------------------------------
      If(.not.conver) then
        call message('SCFMAIN',
     1   '**WARNING** Incomplete SCF convergence ',0,0)
        ierrflag(1)=1
      Else
        ierrflag(1)=0
      End if
c
      RETURN
      END
c=======================================================================
        subroutine print_atomic_numbers(na,iout,ian)
       implicit real*8 (a-h,o-z)
       dimension ian(na)
       write(iout,*) 'Atomic numbers before calling DFT-D3'
       write(iout,'((i3,i5))') (i,ian(i),i=1,na)
       end
c=======================================================================
