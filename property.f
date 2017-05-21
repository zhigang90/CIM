c ==================================================================
c  PROPERTY MODULE            JB   April 2000
c ==================================================================
c
      subroutine preprop(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval
c
c  reads the PROP line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=11)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character cdum*20
      character*256 jobname
c
      parameter (IUnit=1)
c
      data options/'prop','spin','efg ','thre','grid','radf',
     $             'lmax','nrad','nang','prin','fact'/
      data ioptyp/21,0,0,11,11,11,1,1,1,1,11/
c
c -- deal with jobname header
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c -- what properties are being requested?
      if(ifound(2).eq.1) cdum = 'spin'
      if(ifound(3).eq.1) then
        if(ifound(2).eq.1) then
          cdum = 'spin   efg'
        else
          cdum = 'efg'
        endif
      endif
c
c -- grid factor
      If(ifound(5).eq.1) then
        factor = ropval(1,5)
      else
        factor = 1.25d0
      endif
c ......................................................
c -- grid factor (obsolescent)
      If(ifound(11).eq.1) then
        factor = ropval(1,11)
      else
        factor = 1.25d0
      endif
c ......................................................
c
c -- Chipman radial factor
      If(ifound(6).eq.1) then
        r0f = ropval(1,6)
      else
        r0f = 0.35d0
      endif
c
c - Maximum angular momentum in spherical harmonics
      If(ifound(7).eq.1) then
        LMax = iopval(1,7)
      Else
        LMax = 4
      EndIf
c
c -- number of radial/angular grid points
      NRad = 0
      NAng = 0
      If(ifound(8).eq.1) NRad = iopval(1,8)
      If(ifound(9).eq.1) NAng = iopval(1,9)
c
      If(NAng.gt.7) call nerror(1,'PROPERTY Module',
     $  'Input error - NAng must be +ve integer between 1 and 7',0,0)
c
      if(ifound(10).eq.1) then
        IPRNT = iopval(1,10)
      else
        IPRNT = 1          ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,5,'$prop',3,idum,rdum,cdum)
      call WrPropG(IUnit,factor,r0f,LMax,NRad,NAng)
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c  =======================================================================
c
      SUBROUTINE PROPERTY(NMem,Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** PROPERTY MODULE **
C
C  This module computes various atomic and molecular properties
C  Currently these are limited to
C    charge density at the nucleus
C    spin density at the nucleus (open shell)
C    electric field gradient
C  all of which are calculated numerically over atomic-centred
C  grids, similar to those used in DFT
C
C  ..........................................................
C
C  PROPERTY itself is simply a "wrapper" which reads the
C  <control> file to determine the size of the current
C  system and the program option requested and, based on
C  this information, allocates the necessary memory.
C  It then calls PROPMAIN which is itself another "wrapper"
C  which completes job input and is responsible for all
C  other I/O. PROPMAIN calls various routines to compute
C  the desired properties
C  .............................................................
C
C
      CHARACTER jobname*256,cdum*20
      Logical rhf
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    number of atoms
C    multiplicity
C    number of alpha/closed-shell orbitals
C    number of beta orbitals
C    number of radial grid points
C    maximum angular momentum in spherical harmonics
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,rdum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,rdum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,rdum,cdum)
      call RdPropG(IUnit,factor,r0f,LMax,NRad,NAng)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read from the <sym> file
C    number of symmetry operations
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c ............................................................
c -- for the properties, dummy atoms are ignored
c
      NAtoms = NAtoms-Ndum1-Ndum2
c ............................................................
c
      rhf = (NBeta.EQ.0.AND.IMult.EQ.1)
C
C  Read from the <basis> file
C    number of basis functions, primitives and shells
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$nbasis',1,NBas,rdum,cdum)
      call rdcntrl(IUnit,7,'$nshell',1,NShell,rdum,cdum)
      call rdcntrl(IUnit,6,'$nprim',1,NPrm,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  get number of primitive shells from old Texas depository
C
      call getival('nsh ',nsh)
C
C  set maximum number of radial grid points
C
      numR = MAX(NRad,400)
C
C  get number of slaves (if parallel)
C
      call getival('nslv',nslv)
C
C
C  Now get the memory
C
      NScr = NBas*NBas + NBas + 13*NShell + 5*NAtoms
      IMem = 19*NAtoms + 13*NPrm + 12*NShell + 9*NTrans +
     $          NAtoms*NTrans + NAlpha*NBas + 7*NAlpha
      If(.NOT.rhf) Then
       NScr = NScr + NBas + 3*NBeta + 3
       IMem = IMem + NBeta*NBas + 7*NBeta
      EndIf
c
      If(nslv.eq.0) Then
       IMem = IMem + NAtoms + NShell + nsh + 2*NAtoms**2 +
     $          2*numR + 4*1130 + (NAlpha+1)*(LMax+1)*(2*LMax+1)
       If(.NOT.rhf) IMem = IMem +  NBeta*(LMax+1)*(2*LMax+1)
      EndIf
      IMem = IMem + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,8,'PROPERTY')
C
C
C  Allocate permanent memory pointers
C
      ICHG  = iptr                         !  charge density per atom
      ISPN  = ICHG + 2*NAtoms              !  spin density per atom
      IEFG  = ISPN + 2*NAtoms              !  EFG per atom
      IAN   = IEFG + 9*NAtoms              !  atomic numbers
      IQA   = IAN  + NAtoms                !  actual atomic charges
      IXC   = IQA  + NAtoms                !  geometry
      IUQ   = IXC  + 3*NAtoms              !  list of symmetry-unique atoms
      ITN   = IUQ  + NAtoms                !  symmetry operations as 3x3 matrices
      INQ   = ITN  + 9*NTrans              !  list of atomic equivalences
      IBAS  = INQ  + NAtoms*NTrans         !  basis function data
      INX   = IBAS + 13*NPrm               !  basis indexing array
      ICMO  = INX  + 12*NShell             !  alpha/closed-shell MOs
      ISYA  = ICMO + NAlpha*NBas           !  MO symmetry transformation
      IEnd  = ISYA + NAlpha*7
c
      If(.NOT.rhf) Then
       ICMB = IEnd                         !  beta MOs
       ISYB = ICMB + NBeta*NBas            !  beta MO symmetry transformation
       IEnd = ISYB + NBeta*7
      Else
       ICMB = 1
       ISYB = 1
      EndIf
C
C  Memory that is needed only on the slaves
C
      IF(nslv.EQ.0) THEN
       IEXP  = IEnd                        !  minimum exponents per shell
       IBL   = IEXP + NShell               !  precomputed normalization factors
       IAIJ  = IBL  + nsh                  !  precomputed array of Becke weights
       IDST  = IAIJ + NAtoms*NAtoms        !  nearest neighbour distances
       IRST  = IDST + NAtoms               !  inverse interatomic distances
       IRAD  = IRST + NAtoms*NAtoms        !  radial points
       IRWT  = IRAD + numR                 !  radial weights
       IXXA  = IRWT + numR                 !  angular quadrature grid
       IWTA  = IXXA + 3*1130               !  angular weights
       IYLM  = IWTA + 1130                 !  spherical harmonics
       ICA   = IYLM + (LMax+1)*(2*LMax+1)  !  spherically-averaged MOs
       IEnd  = ICA  + NAlpha*(LMax+1)*(2*LMax+1)
c
       If(.NOT.rhf) Then
        ICB  = IEnd                        !  spherically-averaged beta MOs
        IEnd = ICB  + NBeta*(LMax+1)*(2*LMax+1)
       Else
        ICB = 1
       EndIf
      ENDIF
C
C  scratch memory
C
      IScr = IEnd
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,8,'PROPMAIN')
c
c -- set memory mark for TEXAS memory manager
      call getmem(IMem,lastx)
C
C
C  ----------------------------------------------------------------------
C
      CALL Para_PROPMAIN(NAtoms,  Z(IAN),  Z(IQA),  Z(IXC),  NBas,
     $                   NShell,  nsh,     Z(IBAS), Z(INX),  factor,
     $                   r0f,     LMax,    NRad,    NAng,    Z(IEXP),
     $                   Z(IBL),  Z(IAIJ), Z(IDST), Z(IRST), Z(IRAD),
     $                   Z(IRWT), Z(IXXA), Z(IWTA), NAlpha,  NBeta,
     $                   rhf,     Z(IYLM), Z(ICMO), Z(ICA),  Z(ICMB),
     $                   Z(ICB),  NTrans,  Z(IUQ),  Z(ITN),  Z(INQ),
     $                   Z(ISYA), Z(ISYB), Z(ICHG), Z(ISPN), Z(IEFG),
     $                   NScr,    Z(IScr), IErr)
C
C  ----------------------------------------------------------------------
C
C  free memory used in this routine
C
        call retmem(1)
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE PROPMAIN(NAtoms, IAN,    QA,     XC,     NBas,
     $                    NShell, nsh,    BASDAT, INX,    factor,
     $                    r0f,    LMax,   NRad,   NAng,   ExpMIN,
     $                    BL,     AIJ,    DISTN,  RDIST,  radii,
     $                    radwght,XXA,    WTA,    NAlpha, NBeta,
     $                    rhf,    Ylm,    CMO,    CA,     CMOB,
     $                    CB,     NTrans, IUNQ,   TRANS,  NEqATM,
     $                    SIGNA,  SIGNB,  XCharg, XSpin,  XEFG,
     $                    NScr,   Z,      IErr)

      use memory, block => bl

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for property module
C  PROPMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION IAN(NAtoms),QA(NAtoms),XC(3,NAtoms),BASDAT(13,*),
     $          INX(12,*),ExpMIN(NShell),BL(*),DISTN(NAtoms),
     $          AIJ(NAtoms,NAtoms),RDIST(NAtoms,NAtoms),
     $          radii(*),radwght(*),XXA(3,1130),WTA(1130),
     $          Ylm(0:LMax,-LMax:LMax),CMO(NAlpha,NBas),
     $          CA(NAlpha,0:LMax,-LMax:LMax),CMOB(NBeta,NBas),
     $          CB(NAlpha,0:LMax,-LMax:LMax)
      DIMENSION IUNQ(NAtoms),TRANS(9,NAtoms),NEqATM(NAtoms,NTrans)
      INTEGER   SIGNA(NAlpha,7),SIGNB(NAlpha,7)
      DIMENSION XCharg(NAtoms,2),XSpin(NAtoms,2),XEFG(9,NAtoms)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
      DIMENSION RM(3,3)
      CHARACTER GROUP*4,cdum*20
      Character*256 jobname,MOS,MOB
      Logical Symflag,abelian,rhf,spin,efg,full
c  ...................................................................
c     common /intbl/maxsh,ixx(100)     ! old texas integer depository
      common /symm/nsym,nsy(7)         ! old texas symmetry data
c ....................................................................
      DIMENSION Z(NScr)
c
      PARAMETER (Zero=0.0d0,One=1.0d0,Three=3.0d0)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      Data IradQ/0/
c
      Common /job/jobname,lenJ
C
C
C  initialize
C
      IOut=ioutfil('iout')
      ICon=ioutfil('icon')
      thrsh = 1.0d-12
C
C set up the MOs filenames MOS and MOB
C
      call tstchval('mosfname',iyes)
      if(iyes.eq.1) then
        call getchval('mosfname',MOS)
      else
        MOS=jobname(1:lenJ)
      end if
      call rmblan2(MOS,256,lenM)
      if(lenM.le.76) MOS(lenM+1:lenM+4)='.mos'
      lenM=lenM+4
      MOB=MOS
      MOB(lenM:lenM)='b'
C
C ...............................................................
C  Read from the <sym> file
C    rotation matrix
C    point group symbol
C    number of degrees of freedom
C    number of symmetry-unique atoms
C    symmetry operations
C    symmetry-equivalent atoms array
C  This information read even if system is C1
C
      Symflag = NTrans.GT.1
      CALL RdSYM(Symflag,NAtoms, RM,     GROUP,  NTrans,
     $           NDEG,   NQ,     IUNQ,   TRANS,  NEqATM)
C
C  determine if system has abelian symmetry
C
      abelian = GROUP.EQ.'cs  '.OR.GROUP.EQ.'ci  '.OR.
     $          GROUP.EQ.'c2  '.OR.GROUP.EQ.'c2v '.OR.
     $          GROUP.EQ.'c2h '.OR.GROUP.EQ.'d2  '.OR.
     $          GROUP.EQ.'d2h '
cc      abelian = .true.
C
C  Read from the <control> file
C     print flag
C     which properties to be computed
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call rdcntrl(IUnit,5,'$prop',3,idum,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  which properties?
C
      spin = cdum(1:4).eq.'spin'
      efg = cdum(1:3).eq.'efg'.OR.cdum(5:7).eq.'efg'
C
C  print header
C
      WRITE(IOut,1000)
      WRITE(ICon,1000)
c
      If(IPrnt.gt.1) Then
       write(iout,*) ' Rassolov/Chipman radial factor: ',r0f
       write(iout,*) ' maximum angular order (LMax):    ',LMax
       write(iout,*) ' radial grid factor:             ',factor
       write(iout,*) ' properties to be computed:  ',cdum(1:4),
     $               '  ',cdum(5:7)
      EndIf
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XC,-1,jnk)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read the basis set data
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdbasis(IUnit,  NAtoms, AtSymb, XC,     Z(1),
     $             Z(1+NAtoms),BASDAT)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C ....................................................................
C  Be Careful with basis order!
C  We will use Wolinski order for property module
C  Wolinski ordering (for integral evaluation and hence for the MOs)
C  is PER SHELL (e.g. all S functions together regardless of atom)
C
C  The following routines are relevant here
C   SortBAS1  -  sorts basis per shell (but NOT full Wolinski)
C   reorder   -  orders basis from SortBAS1 into Wolinski
C ....................................................................

      CALL SortBAS1(NAtoms,NShell,Z(1+NAtoms),INX)
      CALL normaliz(NShell,INX,BASDAT)
C
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
C  get actual atomic charges (there may, e.g., be ghost atoms)
C
      CALL GetAtChrg(NAtoms,QA)
C
C  need to reorder basis functions to Wolinski order
C
      i1 = 1
      i2 = i1 + NBas
      i3 = i2 + NShell
      i4 = i3 + 12*NShell
      IEnd = i4 + NBas*NBas - 1
      CALL MemCHK(NScr,IEnd,8,'PROPMAIN')
c
      Call reorder(NShell,INX,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinski order
C
C ================================================================
C  get location of abelian symmetry data in old texas depository
C
      call getival('ngener',ngen)         ! number of generators
      nupr = 1                            ! not set if symmetry not used
      ifp = 1                             ! ditto
      If(nsym.gt.0) Then
       call getival('SymNuPr1',nupr)      ! equivalent atoms
       call getival('SymFunPr',ifp)       ! equivalent basis functions
      EndIf
C ================================================================
C
C  now read the old MOs (from the existing binary file)
C  and determine how MOs transform over the abelian symmetry
C  operations (i.e., whether they change sign or not)
C
      itype = 1
      CALL ReadMOS(NBas,Z(i4),jnk,.False.,lenM,MOS(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(1,'Properties module',
     $   'MOS File Does Not Exist',0,0)
c
      If(abelian)
     $ CALL FindMOSym(NBas,   NAlpha, NSym,   nsy,    block(ifp),
     $                Z(i4),  Z(i1),  Z(i2),  SIGNA,  full)
C
C  now transpose MOs
C
      call trspmo(Z(i4),NBas,CMO,NAlpha)
C
C  repeat for beta MOs (if unrestricted wavefunction)
C
      If(NBeta.GT.0) Then
        CALL ReadMOS(NBas,Z(i4),jnk,.False.,lenM,MOB(1:lenM),itype,IErr)
        If(IErr.NE.0) Call nerror(1,'Properties module',
     $     'MOB File Does Not Exist',0,0)
c
        If(abelian.AND.full)
     $   CALL FindMOSym(NBas,   NBeta,  NSym,   nsy,  block(ifp),
     $                  Z(i4),  Z(i1),  Z(i2),  SIGNB,  full)
        call trspmo(Z(i4),NBas,CMOB,NBeta)
      EndIf
C
C
C  Precalculate data for the numerical grid
c  ----------------------------------------
c
c -- calculate inverse atomic distances, Becke aij parameters
c -- and angular grid and weights prior to full grid construction
       call PreGRID(NAtoms, XC,     IAN,    rdum,   rdum,
     $              DISTN,  RDIST,  AIJ,    XXA,    WTA)
c
c -- get array of smallest exponent per shell
       call GetEXP(nsh,NShell,BASDAT,INX,ExpMIN)
c
c -- precompute shell normalization factors
       call AOInit(NShell,BASDAT,INX,BL)
C
C  switch off symmetry for charge/spin density calculation for
C  non-abelian point group and if MO symmetry analysis incomplete
C
      MSym = 0
      If(abelian.AND.full) MSym=NSym
C
C
C  Now calculate the charge and spin densities
C
C  ------------------------------------------------------------
      IF(spin) THEN
C
C  Charge/Spin density
C  -------------------
C
        DO 100 IAtm=1,NQ
        ICntr = IUNQ(IAtm)
        If(rhf) Then
          CALL NucSPINC(ICntr,  NAtoms, XC,     IAN,    QA,
     $                  MSym,   NGen,   nsy, block(nupr), IPrnt,
     $                  NRad,   NAng,   factor, r0f,    DISTN,
     $                  AIJ,    RDIST,  XXA,    WTA,    thrsh,
     $                  radii,  radwght,NBas,   NShell, BASDAT,
     $                  INX,    BL,     ExpMIN, LMax,   Ylm,
     $                  NAlpha, CMO,    CA,     SIGNA,  NScr,
     $                  Z,      XCharg)
c
        Else
          CALL NucSPINU(ICntr,  NAtoms, XC,     IAN,    QA,
     $                  MSym,   NGen,   nsy, block(nupr), IPrnt,
     $                  NRad,   NAng,   factor, r0f,    DISTN,
     $                  AIJ,    RDIST,  XXA,    WTA,    thrsh,
     $                  radii,  radwght,NBas,   NShell, BASDAT,
     $                  INX,    BL,     ExpMIN, LMax,   Ylm,
     $                  NAlpha, NBeta,  CMO,    CMOB,   CA,
     $                  CB,     SIGNA,  SIGNB,  NScr,   Z,
     $                  XCharg, XSpin)
c
        EndIf
 100    CONTINUE
C
C  print out results
C  -----------------
C
        IF(NTrans.GT.1) THEN
          DO 200 IAtm=1,NQ
          ICntr = IUNQ(IAtm)
C
C  generate all symmetry-related centres from current
C  symmetry-unique centre
C
          DO 190 IOP=1,NTrans
          JAtm = NEqATM(ICntr,IOP)
          If(JAtm.NE.ICntr) Then
            XCharg(JAtm,1) = XCharg(ICntr,1)
            XCharg(JAtm,2) = XCharg(ICntr,2)
            If(.NOT.rhf) Then
              XSpin(JAtm,1) = XSpin(ICntr,1)
              XSpin(JAtm,2) = XSpin(ICntr,2)
            EndIf
          EndIf
 190      CONTINUE
 200      CONTINUE
        ENDIF
C
C  now print out the results
C
        WRITE(IOut,1100)
        WRITE(ICon,1100)
c
        If(rhf) Then
          WRITE(IOut,1200)
          WRITE(ICon,1200)
        Else
          WRITE(IOut,1300)
          WRITE(ICon,1300)
        EndIf
c
        DO 300 IAtm=1,NAtoms
        If(rhf) Then
          WRITE(IOut,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2)
          WRITE(ICon,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2)
        Else
          WRITE(IOut,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2),
     $                   XSpin(IAtm,1),XSpin(IAtm,2)
          WRITE(ICon,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2),
     $                   XSpin(IAtm,1),XSpin(IAtm,2)
        EndIf
 300    CONTINUE
cc
      ENDIF
C
      IF(efg) THEN
C
C  Electric Field Gradient
C  -----------------------
C
        CALL ZeroIT(XEFG,9*NAtoms)
c
c *******************************************
c -- Symmetry not working properly
c -- switch off until fixed
        MSym = NSym
        call TempSYM(NAtoms,NSym,NGen,NQ,IUNQ)
c *******************************************
C
        DO 400 IAtm=1,NQ    ! WARNING  needs to be changed to abelian
        ICntr = IUNQ(IAtm)
        CALL NucEFG(ICntr,  NAtoms, XC,     IAN,    QA,
     $              NSym,   NGen,   nsy, block(nupr), NQ,
     $              IUNQ,   IPrnt,  NRad,   NAng,   IradQ,
     $              factor, DISTN,  AIJ,    RDIST,  XXA,
     $              WTA,    thrsh,  radii,  radwght,NBas,
     $              NShell, BASDAT, INX,    BL,     ExpMIN,
     $              rhf,    NAlpha, NBeta,  CMO,    CMOB,
     $              NScr,   Z,      XEFG)
 400    CONTINUE
C
C  print out results
C  -----------------
C
        IF(NSym.GT.0) THEN
          DO 500 IAtm=1,NQ
          ICntr = IUNQ(IAtm)  ! WARNING - This is WRONG  needs to be abelian
C
C  scale EFG tensor by symmetry factor
C
          CALL RedEFG(NAtoms,ICntr,NGen,block(nupr),XEFG(1,ICntr))
C
C  generate all symmetry-related centres from current
C  symmetry-unique centre
C
          CALL SymEFG(NAtoms,ICntr,NGen,nsy,block(nupr),XEFG)
c
 500      CONTINUE
        ENDIF
c
c *******************************************
c -- restore real NSym
        NSym = MSym
c *******************************************
C
C  now print out the results
C
        WRITE(IOut,1500)
        WRITE(ICon,1500)
c
cc        DO 600 IAtm=1,NAtoms
cc        XEFG(2,IAtm) = Three*XEFG(2,IAtm)
cc        XEFG(4,IAtm) = Three*XEFG(4,IAtm)
cc        XEFG(5,IAtm) = Three*XEFG(5,IAtm)
c
cc        WRITE(IOut,1600) IAtm,AtSymb(IAtm),
cc     $                   XEFG(7,IAtm),XEFG(8,IAtm),XEFG(9,IAtm),
cc        WRITE(ICon,1600) IAtm,AtSymb(IAtm),
cc     $                   XEFG(7,IAtm),XEFG(8,IAtm),XEFG(9,IAtm),
cc     $                   XEFG(2,IAtm),XEFG(4,IAtm),XEFG(5,IAtm)
cc 600    CONTINUE
C
C  Now do Tensor representation
C
        DO 650 IAtm=1,NAtoms
        XEFG(2,IAtm) = Three*XEFG(2,IAtm)
        XEFG(4,IAtm) = Three*XEFG(4,IAtm)
        XEFG(5,IAtm) = Three*XEFG(5,IAtm)
c
        WRITE(IOut,1600) IAtm,AtSymb(IAtm),
     $                   XEFG(1,IAtm),XEFG(3,IAtm),XEFG(6,IAtm),
     $                   XEFG(2,IAtm),XEFG(4,IAtm),XEFG(5,IAtm)
        WRITE(ICon,1600) IAtm,AtSymb(IAtm),
     $                   XEFG(1,IAtm),XEFG(3,IAtm),XEFG(6,IAtm),
     $                   XEFG(2,IAtm),XEFG(4,IAtm),XEFG(5,IAtm)
 650    CONTINUE
C
C  diagonalize tensor for each atom and print out eigenvalues
C
        WRITE(IOut,1800)
        WRITE(ICon,1800)
c
        I1 = 1
        I2 = I1 + 9
        I3 = I2 + 9
c
        DO 700 IAtm=1,NAtoms
        CALL EXPAND(3,XEFG(1,IAtm),Z(I1))
        CALL DiagMAT(Z(I1),3,Z(I2),Z(I3),DISTN,IErr)
c
        If(IErr.NE.0) Then
          WRITE(IOut,2000) IAtm
          GO TO 700
        EndIf
c
        WRITE(IOut,1600) IAtm,AtSymb(IATM),DISTN(1),DISTN(2),DISTN(3)
        WRITE(ICon,1600) IAtm,AtSymb(IATM),DISTN(1),DISTN(2),DISTN(3)
 700    CONTINUE
cc
      ENDIF
C  ------------------------------------------------------------
C
      RETURN
c
 1000 FORMAT(/,'                          ** NUCLEAR PROPERTIES **')
 1100 FORMAT(/,' <== CHARGE AND SPIN DENSITY AT THE NUCLEUS',
     $          '  (atomic units)  ==>')
 1200 FORMAT(/,'                    CHARGE  DENSITY',/,
     $         '    Atom     Delta Function  Rassolov/Chipman')
 1300 FORMAT(/,20X,'CHARGE  DENSITY',20X,'SPIN  DENSITY',/,
     $         '  Atom     Delta Function  Rassolov/Chipman',
     $         '  Delta Function  Rassolov/Chipman')
 1400 FORMAT(I4,1X,A4,2X,4(F12.5,5X))
 1500 FORMAT(/,' ELECTRIC FIELD GRADIENT COMPONENTS in atomic units',
     $         ' (Tensor Representation)',
     $        /,'  Atom       3XX-RR     3YY-RR     3ZZ-RR',
     $          '      3XY        3XZ        3YZ')
 1600 FORMAT(I4,1X,A4,1X,6(F10.5,1X))
 1800 FORMAT(/,' INVARIANT EIGENVALUES OF THE EFG TENSOR',
     $       /,'  Atom',15X,'Eigenvalues')
 2000 FORMAT(2X,'***ERROR*** Unable to diagonalize EFG tensor',
     $          ' for atom: ',I4)
c
      END
c  =======================================================================
c
      SUBROUTINE NucSPINC(ICntr,  NAtoms, XNuc,   IAN,    QA,
     $                    NSym,   NGen,   ISYM,   NEqATM, IPrnt,
     $                    NRad,   NAng,   factor, r0f,    DISTN,
     $                    AIJ,    RDIST,  XXA,    WTA,    thrsh,
     $                    radii,  radwght,NBas,   NShell, BASDAT,
     $                    INX,    BL,     ExpMIN, LMax,   Ylm,
     $                    NOcc,   CMO,    CA,     SIGN,   NScr,
     $                    Z,      XCharg)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates spin/charge density at nucleus using density operator method
C  of Rassolov and Chipman (J.Chem.Phys. 105 (1996) 1470)
C
C  ** CLOSED SHELL **
C
C    originally written by Bing Wang based on modified JB DFT code
C    revised by JB April 2000
C
C  ARGUMENTS
C
C  ICntr   -  current atom for which spin density is being calculated
C  NAtoms  -  number of atoms
C  XNuc    -  nuclear coordinates
C  IAN     -  atomic number
C  QA      -  actual nuclear charge (there may, e.g., be ghost atoms)
C  NSym    -  number of abelian symmetry operations
C  NGen    -  number of generators
C  ISYM    -  list of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  IPrnt   -  print flag
C  NRad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  NAng    -  angular quadrature order for numerical grid
C             (if zero, will be determined in grid generation routines)
C  factor  -  determines  grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  r0f     -  radial factor for Rassolov/Chipman function
C  DISTN   -  distance to nearest neighbour for each atom
C  AIJ     -  Becke grid parameters (see A2 in Becke article)
C  RDIST   -  inverse interatomic distance array
C  XXA     -  storage for angular quadrature grid
C  WTA     -  storage for angular quadrature weights
C  thrsh   -  threshold for neglecting contribution
C  radii   -  radial grid points
C  radwght -  radial weights
C  NBas    -  total number of basis functions
C  NShell  -  total number of shells
C  BASDAT  -  basis set data for TEXAS
C  INX     -  more basis set data for TEXAS
C  BL      -  precomputed normalization factors per primitive shell
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  LMax    -  Maximum L value for spherical harmonics
C  Ylm     -  spherical harmonics
C  NOcc    -  number of occupied MOs
C  CMO     -  occupied MO coefficients (transposed)
C  CA      -  spherically averaged MOs
C  SIGN    -  how MOs transform over abelian symmetry operations
C  NScr    -  size of general scratch array
C             (5*NBas + 4*NOcc + 6)
C  Z       -  general scratch array
C
C  on exit
C
C  XCharg  -  charge density per atom
C
C  ----------------------------------------------------------------
C  references
C
C     "A multicenter numerical integration scheme for polyatomic
C      molecules"
C      A.D.Becke      J.Chem.Phys.  88 (1988) 2547
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),QA(NAtoms),ExpMIN(NShell),
     $          BL(*),BASDAT(13,*),INX(12,*),CMO(NOcc,NBas)
      DIMENSION AIJ(NAtoms,NAtoms),RDIST(NAtoms,NAtoms)
      DIMENSION DISTN(NAtoms),radii(*),radwght(*)
      DIMENSION XXA(3,1130),WTA(1130),INDEX(8),XXAN(3)
      DIMENSION ISYM(NSym),NEqATM(NAtoms,NSym)
      DIMENSION Ylm(0:LMax,-LMax:LMax)
      DIMENSION CA(NOcc,0:LMax,-LMax:LMax)
      DIMENSION XCharg(NAtoms,2)
      DIMENSION Z(NScr),LF(NSym),nbf(2)
      INTEGER SIGN(NOcc,7)
      Integer SFac
C
      PARAMETER (Zero=0.0d0,One=1.0d0,Two=2.0d0,Four=4.0d0)
      PARAMETER (MaxRad=400)
      Data INDEX/1,15,41,91,201,395,697,1131/
C
C
C  This method involves numerical integration over an atom-centred
C  grid, similar to current DFT implementations. However, the numerical
C  integration scheme is different. Angular integration is done using
C  spherical harmonics. Radial integration uses standard Gauss-Legendre
C  quadrature.
C
C
      PI = Four*ATAN(One)
      XX = XNuc(1,ICntr)
      YY = XNuc(2,ICntr)
      ZZ = XNuc(3,ICntr)
c
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
C
C  allocate scratch pointers
C
      i1  = 1
      i2  = i1  + NAtoms
      i3  = i2  + NAtoms
c
      iao = 1                         ! AO values over grid points
      imo = iao + NBas                ! MO values over grid points
      inb = imo + NOcc                ! index array for "non-zero" AOs
      iaox = inb + NBas               ! AO derivatives
      imox = iaox + 3*NBas            ! MO derivatives
      idx  = imox + 3*NOcc            ! density derivatives
      IEnd = idx  + 3 - 1
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'NucSPINC')
C
C ===============================================================
C  First, calculate charge density at nucleus using
C  delta function method
C
      CALL DeltaFUNC(ICntr,  NBas,   NShell, NOcc,   XNuc,
     $               BASDAT, INX,    BL,     ExpCUT, ExpMIN,
     $               CMO,    Z(iao), Z(imo), DCharg)
C
C ===============================================================
C
C
C  Begin calculation for this centre
C  ---------------------------------
C
      ino = IAN(ICntr)      ! atomic number
      qai = QA(ICntr)       ! atomic charge
c
c -- set radial factor
      IF(ino.LE.10) THEN
        r0 = r0f/dble(ino)
      ELSE IF(11.LE.ino.AND.ino.LE.18) THEN
        r0 = 1.5d0*r0f/dble(ino)
      ELSE
        r0 = 4.5d0*r0f/dble(ino)
      ENDIF
      r02 = r0**2
C
C  -----------------------------------
C    RADIAL GRID
C  -----------------------------------
C  determine number of radial grid points
C
      thenum = factor*50.0d0
c
      IF(NRad.LE.0) THEN
        IF(ino.LE.2) THEN
         numR = thenum*0.6d0                     ! 30 points
        ELSE IF(3.LE.ino.AND.ino.LE.10) THEN
         numR = thenum                           ! 50 points
        ELSE IF(11.LE.ino.AND.ino.LE.18) THEN
         numR = thenum*1.4d0                     ! 70 points
        ELSE IF(19.LE.ino.AND.ino.LE.36) THEN
         numR = thenum*1.8d0                     ! 90 points
        ELSE
         numR = thenum*2.4d0                     ! 120 points
        ENDIF
      ELSE
        numR = NRad
      ENDIF
      call f_lush(6)
c
      If(numR.GT.MaxRad) numR=MaxRad
C
C  now get the radial grid
C
      rmin = 0.0d0
      rmax = 15.0d0*r0
      If(ino.GT.10) rmax = 45.0d0*r0
      CALL RadGRID(numR,rmin,rmax,radii,radwght)
c
      If(IPrnt.gt.2) Then
       write(6,*) ' Radial Grid for atom: ',icntr
       write(6,*) ' RMin: ',rmin,' RMax: ',rmax
       write(6,*) ' Number of radial grid points: ',numr
       If(IPrnt.gt.3) Then
        do i=1,numR
        write(6,*) i,radii(i),radwght(i)
        enddo
       EndIf
      EndIf
C
C  ------------------------------------------
C
      XD1A = Zero
      XD2A = Zero
      XD3A = Zero
      XD4A = Zero
      XUC = Zero
      XLC = Zero
c
      NP  = 1
      IVal = 1
      SFac = 1
      LNum = 0
c
      iang = NAng
      If(iang.EQ.0) iang=4        ! default
      IA1 = INDEX(iang)
      IA2 = INDEX(iang+1)-1
C
C  loop over radial grid points on atom IAtom
C
      DO 50 irad=1,numR
      rr = radii(irad)
      rr2 = rr**2
      rwght = radwght(irad)
c
      DenA = Zero
      DenRAX = Zero
      Call ZeroIT(CA,NOcc*(LMax+1)*(2*LMax+1))
C
C  loop over angular grid points on atom IAtom
C
      DO 30 iang=IA1,IA2
      xap = XXA(1,iang)
      yap = XXA(2,iang)
      zap = XXA(3,iang)
      XXAN(1) = xap*rr + XX
      XXAN(2) = yap*rr + YY
      XXAN(3) = zap*rr + ZZ
C
C  get symmetry factor
C
      If(NSym.GT.0) Then
       CALL SymFAC(NAtoms, ICntr,  XX,     YY,     ZZ,
     $             XXAN(1),XXAN(2),XXAN(3),NSym,   NGen,
     $             ISYM,   NEqATM, SFac,   LNum,   LF)
       If(SFac.EQ.0) GO TO 30
      EndIf
C
C  calculate the Becke weight
C
      call WBeckeG(ICntr, XXAN(1), XXAN(2), XXAN(3),
     $             NAtoms, XNuc, RDIST, AIJ,
     $             Z(i1), Z(i2), Z(i3), Wbecke)
c
      Wtot = WTA(iang)*Wbecke*SFac
C
C  form AO and derivative values over grid points
C  and sort "non-zero" AOs
C
      CALL zeroit(Z(iao),NBas)
      CALL zeroit(Z(iaox),3*NBas)
      CALL AOGrad(NP,     NShell, NBas,   XXAN,
     $            XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $            ExpMIN, Z(iao), Z(iaox))
      CALL AOXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $             nbf,    Z(inb), Z(iao), Z(iaox),VMX)
C
C  form density at current grid point
C
      CALL GetDENS(NP,     NBas,   nbf,    NOcc,   CMO,
     $             Z(iao), Z(inb), Z(imo), DVal)
      DVal = DVal+DVal               ! closed-shell density
C
C  ..........................................................
C  on return from <GetDENS>  DVal holds current density
C  ..........................................................
C
C  prepare for numerical integration
C  neglect point if density below threshold
C
      If(DVal.LT.thrsh) GO TO 30
C
C  form density derivatives
C
      CALL GetDENDeriv(NP,    NBas,   thrsh,  nbf,    NOcc,
     $                 CMO,  Z(iaox), Z(inb), DVal,   Z(imo),
     $                Z(imox),Z(idx))
C
c -- spherically averaged electron density
      DenA = DenA + DVal*Wtot
c
c -- spherically averaged MOs
      call spherical_harmonics(xap,yap,zap,LMax,Ylm)
      call SphAvMO(NOcc,   LMax,   Wtot,   LNum,   LF,
     $             SIGN,   Ylm,    Z(imo), CA)
c
c -- spherically averaged MO derivatives
      call xyz_ang(xap,yap,zap,th,phi)
      CALL SphAvMOX(NOcc,th,phi,Z(imox),Z(imo),DenRX)
      DenRAX = DenRAX + Two*DenRX*Wtot
cc
 30   CONTINUE
C
C  numerical integration for charge density
C
c -- U and D part
      XHSFAD = DenA*rwght*exp(-rr2/r02)
      XD1A = XD1A + XHSFAD*(-Two*rr2**1.5d0/r02**3)
      XD2A = XD2A + XHSFAD*(7.0d0*sqrt(rr2)/r02**2)
      XD3A = XD3A + XHSFAD*(-Two/(sqrt(rr2)*r02))
      XD4A = XD4A + DenRAX*exp(-rr2/r02)*rwght
     $                    *(Two*sqrt(rr2)/r02)
      XUC = XUC + XHSFAD/rr2
      XDC = XD1A + XD2A + XD3A + XD4A
c
c -- L part
      CCA = Zero
      DO 40 J=1,NOcc
      DO 40 L=1,LMax
      DO 40 M=-L,L
      CCA = CCA - L*(L+1)*CA(J,L,M)**2
 40   CONTINUE
c
      XLC = XLC + Two*CCA*rwght*EXP(-rr2/r02)/(rr2**1.5d0)
cc
 50   CONTINUE
C
      XUC = XUC*qai
      XLC = XLC*(Four*PI)
      RCharg = (XUC+XDC+XLC)/(Two*PI)
c
      If(IPrnt.GT.2) Then
       write(6,*)' U part CHARGE integral is', XUC
       write(6,*)' D1 part CHARGE integral is', XD1A
       write(6,*)' D2 part CHARGE integral is', XD2A
       write(6,*)' D3 part CHARGE integral is', XD3A
       write(6,*)' D4 part CHARGE integral is', XD4A
       write(6,*)' D part CHARGE integral is', XDC
       write(6,*)' L part CHARGE integral is', XLC
       write(6,*)' Total charge density is',RCharg
      EndIf
c
      XCharg(ICntr,1) = DCharg
      XCharg(ICntr,2) = RCharg
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE DeltaFUNC(ICntr,  NBas,   NShell, NOcc,   XNuc,
     $                     BASDAT, INX,    BL,     ExpCUT, ExpMIN,
     $                     CMO,    VAO,    VMOA,   DCharg)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine calculates the charge density at the nucleus
C  by the delta function method
C
      Dimension XNuc(3,*), BL(*), ExpMIN(NShell), BASDAT(13,*),
     $          INX(12,*), CMO(NOcc,NBas)
      Dimension VAO(NBas),VMOA(NOcc)
C
C  get basis function values at the nucleus
C
      CALL ZeroIT(VAO,NBas)
      CALL AOVal(1,      NShell, NBas,   XNuc(1,ICntr),
     $           XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $           ExpMIN, VAO)
C
C  now evaluate the charge density
C
      CALL zeroit(VMOA,NOcc)
      DO 20 I=1,NBas
      Val = VAO(I)
      DO 10 J=1,NOcc
      VMOA(J) = VMOA(J) + CMO(J,I)*Val
 10   CONTINUE
 20   CONTINUE
c
      DCharg = DDOT(NOcc,VMOA,1,VMOA,1)
      DCharg = DCharg + DCharg
C
      RETURN
      END
c ===================================================================
c
      SUBROUTINE SphAvMOX(NOcc,   th,     phi,    VMOX,   VMORX,
     $                    DenRX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  convert Cartesian MO derivatives to spherically-averaged
C  MO derivatives
C
      DIMENSION VMOX(NOcc,3),VMORX(NOcc)
      Parameter (One=1.0d0)
c
      s = sqrt(One-th*th)
      cosp = s*cos(phi)
      sinp = s*sin(phi)
c
      DO 10 I=1,NOcc
      VMORX(I) = cosp*VMOX(I,1)
     $         + sinp*VMOX(I,2)
     $         + th*VMOX(I,3)
 10   CONTINUE
c
      DENRX = DDOT(NOcc,VMORX,1,VMORX,1)
c
      RETURN
      END
c ===================================================================
c
      SUBROUTINE SphAvMO(NOcc,   LMax,   ww,     LNum,   LF,
     $                   SIGN,   Ylm,    VMO,    CA)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  form spherically-averaged MOs over current grid point
C
C  ARGUMENTS
C
C  NOcc    -  number of occupied orbitals
C  LMax    -  maximum L value for spherical harmonics
C  ww      -  weight for current angular grid point
C  LNum    -  number of symmetry operations
C  LF      -  list of operations
C  SIGN    -  array indicating whether or not MOs change
C             sign under symmetry operations
C  Ylm     -  spherical harmonics
C  VMO     -  current MO values
C
C  on exit
C
C  CA      -  spherically averaged MOs
C
C
      DIMENSION VMO(NOcc),Ylm(0:LMax,-LMax:LMax)
      DIMENSION CA(NOcc,0:LMax,-LMax:LMax)
      INTEGER LF(LNum),SIGN(NOcc,*)
      Logical leven,meven
C
C
      IF(LNum.EQ.0) THEN
C
C  NO SYMMETRY
C  -----------
C
       DO 20 L=0,LMax
       DO 20 M=-L,L
       www = Ylm(L,M)*ww
       DO 10 J=1,NOcc
       CA(J,L,M) = CA(J,L,M) + VMO(J)*www
 10    CONTINUE
 20    CONTINUE
cc
      ELSE
C
C  SYMMETRY
C  --------
C
C  For "standard" grid-based integration procedures, symmetry
C  readily be incorporated by computing the relevant quantities
C  at symmetry-unique grid points only, multiplying the normal
C  weight factor by the number of symmetry-equivalent points.
C
C  This will not work for CA, above, as symmetry causes some
C  contributions to CA at symmetry-related grid points to
C  cancel instead of always summing. Which terms cancel and
C  which sum depends on the values of L and M, the symmetry
C  operation and the sign of the MO under the symmetry operation.
C
C  The following terms cancel for the various symmetry operations:
C  (shown in parentheses are Cartesian coordinates which change
C   sign under that symmetry operation)
C
C ..................................
C  (1) MOs which do NOT change sign
C ..................................
C
C       operation 1 (X)      operation 2 (Y)      operation 3 (Z)
C       ---------------      ---------------      ---------------
C   L         M                    M                    M
C   0         -                    -                    -
C   1         1                   -1                    0
C   2       -2,1                 -2,-1                -1,1
C   3      -2,1,3               -3,-2,-1             -2,0,2
C   4     -4,-2,1,3            -4,-3,-2,-1          -3,-1,1,3
C
C
C       operation 4 (XY)     operation 5 (XZ)     operation 6 (YZ)
C       ----------------     ----------------     ----------------
C   L         M                    M                    M
C   0         -                    -                    -
C   1       -1,1                  0,1                 -1,0
C   2       -1,1                 -2,-1                -2,1
C   3    -3,-1,1,3              0,1,2,3            -3,-1,0,2
C   4    -3,-1,1,3            -4,-3,-2,-1          -4,-2,1,3
C
C
C       operation 7 (XYZ)
C       -----------------
C   L         M                  The effects of each operation can
C   0         -                  be derived by combinations of the
C   1      -1,0,1                first 3 operations, eliminating
C   2         -                  any M value that appears twice
C   3  -3,-2,-1,0,1,2,3
C   4         -
C
C ............................
C  (2) MOs which CHANGE sign
C ............................
C
C       operation 1 (X)      operation 2 (Y)      operation 3 (Z)
C       ---------------      ---------------      ---------------
C   L         M                    M                    M
C   0         0                    0                    0
C   1       -1,0                  0,1                 -1,1
C   2      -1,0,2                0,1,2               -2,0,2
C   3     -3,-1,0,2             0,1,2,3             -3,-1,1,3
C   4    -3,-1,0,2,4           0,1,2,3,4           -4,-2,0,2,4
C
C
C       operation 4 (XY)     operation 5 (XZ)     operation 6 (YZ)
C       ----------------     ----------------     ----------------
C   L         M                    M                    M
C   0         0                    0                    0
C   1         0                   -1                    1
C   2      -2,0,2                0,1,2               -1,0,2
C   3      -2,0,2               -3,-2,-1             -2,1,3
C   4   -4,-2,0,2,4            0,1,2,3,4          -3,-1,0,2,4
C
C
C       operation 7 (XYZ)
C       -----------------
C   L         M
C   0         0
C   1         -
C   2    -2,-1,0,1,2
C   3         -
C   4  -4,-3,-2,-1,0,1,2,3,4
C
C
       DO 60 J=1,NOcc
       DO 80 L=0,LMax
       leven = (L/2).EQ.((L+1)/2)
       DO 90 M=-L,L
       meven = (Abs(M)/2).EQ.((Abs(M)+1)/2)
       DO 70 IOP=1,LNum
       LFac = LF(IOP)
c
c -- determine sign of MO under operation IOP
       ISign = SIGN(J,LFac)
cc       write(6,*) ' LFac:',lfac,' J:',j,' SIGN:',sign(j,lfac)
c
       IF(ISign.GT.0) THEN
c
c -- no sign change
        If(LFac.EQ.7.AND..not.leven) GO TO 80
        IF(LFac.EQ.1) THEN
         If( (M.GT.0.AND..not.meven) .OR.
     $       (M.LT.0.AND.meven) ) GO TO 90
        ELSE IF(LFac.EQ.2) THEN
         If(M.LT.0) GO TO 90
        ELSE IF(LFac.EQ.3) THEN
         If( (leven.AND..not.meven) .OR.
     $       (.not.leven.AND.meven) ) GO TO 90
        ELSE IF(LFac.EQ.4) THEN
         If(.not.meven) GO TO 90
        ELSE IF(LFac.EQ.5) THEN
         If( (leven.AND.M.LT.0) .OR.
     $       (.not.leven.AND.M.GE.0) ) GO TO 90
        ELSE
         If( (leven.AND.((meven.AND.M.LT.0).OR.
     $                   (.not.meven.AND.M.GT.0))) .OR.
     $       (.not.leven.AND.((meven.AND.M.GE.0).OR.
     $                        (.not.meven.AND.M.LT.0))) ) GO TO 90
        ENDIF
cc
       ELSE
c
c -- sign change
        If(LFac.EQ.7.AND.leven) GO TO 80
        IF(LFac.EQ.1) THEN
         If( (M.GE.0.AND.meven) .OR.
     $       (M.LT.0.AND..not.meven) ) GO TO 90
        ELSE IF(LFac.EQ.2) THEN
         If(M.GE.0) GO TO 90
        ELSE IF(LFac.EQ.3) THEN
         If( (leven.AND.meven) .OR.
     $       (.not.leven.AND..not.meven) ) GO TO 90
        ELSE IF(LFac.EQ.4) THEN
         If(meven) GO TO 90
        ELSE IF(LFac.EQ.5) THEN
         If( (leven.AND.M.GE.0) .OR.
     $       (.not.leven.AND.M.LT.0) ) GO TO 90
        ELSE
         If( (leven.AND.((meven.AND.M.GE.0).OR.
     $                   (.not.meven.AND.M.LT.0))) .OR.
     $       (.not.leven.AND.((meven.AND.M.LT.0).OR.
     $                        (.not.meven.AND.M.GT.0))) ) GO TO 90
        ENDIF
cc
       ENDIF
 70    CONTINUE
       www = Ylm(L,M)*ww
       CA(J,L,M) = CA(J,L,M) + VMO(J)*www
 90    CONTINUE
 80    CONTINUE
 60    CONTINUE
cc
      ENDIF
C
      RETURN
      END
c ===================================================================
c
      SUBROUTINE NucEFG(ICntr,  NAtoms, XNuc,   IAN,    QA,
     $                  NSym,   NGen,   ISYM,   NEqATM, NQ,
     $                  IUNQ,   IPrnt,  NRad,   NAng,   IradQ,
     $                  factor, DISTN,  AIJ,    RDIST,  XXA,
     $                  WTA,    thrsh,  radii,  radwght,NBas,
     $                  NShell, BASDAT, INX,    BL,     ExpMIN,
     $                  rhf,    NAlpha, NBeta,  CMO,    CMOB,
     $                  NScr,   Z,      XEFG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates Electric Field Gradient at the nucleus
C  Written by JB - May 2000
C
C  ARGUMENTS
C
C  ICntr   -  current atom over which integration is being carried out
C  NAtoms  -  number of atoms
C  XNuc    -  nuclear coordinates
C  IAN     -  atomic numbers
C  QA      -  actual nuclear charge (there may, e.g., be dummy atoms)
C  NSym    -  number of abelian symmetry operations
C  NGen    -  number of generators
C  ISYM    -  list of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of symmetry-unique atoms
C  IPrnt   -  print flag
C  NRad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  NAng    -  angular quadrature order for numerical grid
C             (if zero, will be determined in grid generation routines)
C  IradQ   -  radial quadrature type
C              0 - Euler-Maclaurin (Handy)          default
C              1 - Chebyshev second type (Becke mapping)
C  factor  -  determines grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  DISTN   -  distance to nearest neighbour for each atom
C  AIJ     -  Becke grid parameters (see A2 in Becke article)
C  RDIST   -  inverse interatomic distance array
C  XXA     -  storage for angular quadrature grid
C  WTA     -  storage for angular quadrature weights
C  thrsh   -  threshold for neglecting contribution
C  radii   -  radial grid points
C  radwght -  radial weights
C  NBas    -  total number of basis functions
C  NShell  -  total number of shells
C  BASDAT  -  basis set data for TEXAS
C  INX     -  more basis set data for TEXAS
C  BL      -  precomputed normalization factors per primitive shell
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  rhf     -  logical flag for closed-shell
C  NAlpha  -  number of alpha/closed shell occupied MOs
C  NBeta   -  number of beta MOs
C  CMO     -  occupied alpha/closed shell MO coefficients (transposed)
C  CMOB    -  occupied beta MO coefficients (transposed)
C  NScr    -  size of general scratch array
C             (2*NBas + NAlpha)
C  Z       -  general scratch array
C
C  on exit
C
C  XEFG    -  electric field gradient
C
C  ----------------------------------------------------------------
C  references
C
C     "A multicenter numerical integration scheme for polyatomic
C      molecules"
C      A.D.Becke      J.Chem.Phys.  88 (1988) 2547
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),QA(NAtoms),ExpMIN(NShell),
     $          BL(*),BASDAT(13,*),INX(12,*),IUNQ(NQ),
     $          CMO(NAlpha,NBas),CMOB(NBeta,NBas)
      DIMENSION AIJ(NAtoms,NAtoms),RDIST(NAtoms,NAtoms)
      DIMENSION DISTN(NAtoms),radii(*),radwght(*)
      DIMENSION XXA(3,1130),WTA(1130),INDEX(8),XXAN(9)
      DIMENSION XEFG(9,NAtoms)
      DIMENSION Z(NScr),LF(NSym),nbf(2)
      Integer SFac
      Logical rhf
C
      PARAMETER (Zero=0.0d0,One=1.0d0,Three=3.0d0)
      PARAMETER (MaxRad=400)
      Data INDEX/1,15,41,91,201,395,697,1131/
C
C
C  This method involves numerical integration over an atom-centred
C  grid, similar to current DFT implementations.
C
      XX = XNuc(1,ICntr)
      YY = XNuc(2,ICntr)
      ZZ = XNuc(3,ICntr)
c
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
C
C  allocate scratch pointers
C
      i1  = 1
      i2  = i1  + NAtoms
      i3  = i2  + NAtoms
c
      iao = 1                         ! AO values over grid points
      imo = iao + NBas                ! MO values over grid points
      inb = imo + NAlpha              ! index array for "non-zero" AOs
      IEnd = inb + NBas
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'NucEFGC')
C
C
C  Begin calculation for this centre
C  ---------------------------------
C
      ino = IAN(ICntr)      ! atomic number
C
C  -----------------------------------
C    RADIAL GRID
C  -----------------------------------
C  determine number of radial grid points
C
      thenum = factor*50.0d0
c
      IF(NRad.LE.0) THEN
        IF(ino.LE.2) THEN
         numR = thenum*0.6d0                     ! 30 points
        ELSE IF(3.LE.ino.AND.ino.LE.10) THEN
         numR = thenum                           ! 50 points
        ELSE IF(11.LE.ino.AND.ino.LE.18) THEN
         numR = thenum*1.4d0                     ! 70 points
        ELSE IF(19.LE.ino.AND.ino.LE.36) THEN
         numR = thenum*1.8d0                     ! 90 points
        ELSE
         numR = thenum*2.4d0                     ! 120 points
        ENDIF
      ELSE
        numR = NRad
      ENDIF
c
      If(numR.GT.MaxRad) numR=MaxRad
C
C  get Bragg-Slater atomic radius
C
      CALL braggslaterradius(ino,rbs)            ! rbs in au
      RMax = 20.0d0*rbs
      rm = rbs
C
C  now get the radial grid
C
      CALL RadialGRID(numR,IradQ,ino,rm,radii,radwght)
c
      If(IPrnt.gt.2) Then
       write(6,*) ' Radial Grid for atom: ',icntr
       write(6,*) ' Number of radial grid points: ',numr
       If(IPrnt.gt.3) Then
        do i=1,numR
        write(6,*) i,radii(i),radwght(i)
        enddo
       EndIf
      EndIf
C
C  "Prune" the grid by dividing into regions with different
C   numbers of angular grid points per radial shell
C
      CALL PruneGRID(ino,    rbs,    factor, R1,     R2,
     $               R3,     lang)
C
C  ===============================================================
C
      efgx = Zero
      efgy = Zero
      efgz = Zero
c
      NP  = 1
      IVal = 1
      SFac = 1
C
C  loop over radial grid points on atom IAtom
C
      DO 50 irad=1,numR
      rr = radii(irad)
c ------------------------------------------------------------------
      If(rr.GT.RMax) exit               ! ignore if radius too large
c ------------------------------------------------------------------
      rr2 = rr**2
      rwght = radwght(irad)
C
C  determine angular grid order
C  (i.e., number of grid points in this shell)
C
      IF(NAng.EQ.0) THEN
        If(rr.LT.R1) Then
          iang = lang+2
        Else If(rr.LT.R2.OR.rr.GT.R3) Then
          iang = lang+4
        Else
          iang = lang+6
        EndIf
      ELSE
        iang = NAng
      ENDIF
c
      IA1 = INDEX(iang)
      IA2 = INDEX(iang+1)-1
C
C  loop over angular grid points on atom IAtom
C
      DO 30 iang=IA1,IA2
      xap = XXA(1,iang)
      yap = XXA(2,iang)
      zap = XXA(3,iang)
      XXAN(1) = xap*rr + XX
      XXAN(2) = yap*rr + YY
      XXAN(3) = zap*rr + ZZ
C
C  get symmetry factor
C
      If(NSym.GT.0) Then
       CALL SymFAC(NAtoms, ICntr,  XX,     YY,     ZZ,
     $             XXAN(1),XXAN(2),XXAN(3),NSym,   NGen,
     $             ISYM,   NEqATM, SFac,   LNum,   LF)
       If(SFac.EQ.0) GO TO 30
      EndIf
C
C  calculate the Becke weight
C
      call WBeckeG(ICntr, XXAN(1), XXAN(2), XXAN(3),
     $             NAtoms, XNuc, RDIST, AIJ,
     $             Z(i1), Z(i2), Z(i3), Wbecke)
c
      Wtot = WTA(iang)*Wbecke*SFac
C
C  form AO values over grid points
C  and sort "non-zero" AOs
C
      CALL zeroit(Z(iao),NBas)
      CALL AOVAL(NP,     NShell, NBas,   XXAN,
     $           XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $           ExpMIN, Z(iao))
      CALL AOSort(NP,     NBas,   thint,  Z(iao), nbf,
     $            Z(inb), Z(iao), VMX)
cc      CALL AOSort(NP,NBas,thint,Z(iao),nbf,Z(inb))
C
C  form density at current grid point
C
      CALL GetDENS(NP,     NBas,   nbf,    NAlpha, CMO,
     $             Z(iao), Z(inb), Z(imo), DA)
      If(rhf) Then
       DVal = DA + DA                 ! closed-shell density
      Else
       CALL GetDENS(NP,     NBas,   nbf,    NBeta,  CMOB,
     $              Z(iao), Z(inb), Z(imo), DB)
       DVal = DA + DB                 ! total open-shell density
      EndIf
C
C  ..........................................................
C  on return from <GetDENS>  DVal holds current density
C  ..........................................................
C
C  neglect point if density below threshold
C
      If(DVal.LT.thrsh) GO TO 30
c
c -- electric field gradient
      GFac = Wtot*rwght*DVal
      DO 20 IQ=1,NQ
      IAtm = IUNQ(IQ)
      XXX = XNuc(1,IAtm)
      YYY = XNuc(2,IAtm)
      ZZZ = XNuc(3,IAtm)
      xp = XXAN(1) - XXX
      yp = XXAN(2) - YYY
      zp = XXAN(3) - ZZZ
      x2 = xp**2
      y2 = yp**2
      Z2 = zp**2
      rrr = x2 + y2 + z2
      rr1 = sqrt(rrr)
      Fac = GFac/(rr1**5)
      XEFG(1,IAtm) = XEFG(1,IAtm) + (Three*x2-rrr)*Fac
      XEFG(3,IAtm) = XEFG(3,IAtm) + (Three*y2-rrr)*Fac
      XEFG(6,IAtm) = XEFG(6,IAtm) + (Three*z2-rrr)*Fac
      XEFG(2,IAtm) = XEFG(2,IAtm) + xp*yp*Fac
      XEFG(4,IAtm) = XEFG(4,IAtm) + xp*zp*Fac
      XEFG(5,IAtm) = XEFG(5,IAtm) + yp*zp*Fac
cc      XEFG(7,IAtm) = XEFG(7,IAtm) + x2*Fac
cc      XEFG(8,IAtm) = XEFG(8,IAtm) + y2*Fac
cc      XEFG(9,IAtm) = XEFG(9,IAtm) + z2*Fac
      If(NGen.GT.0) CALL SymGenEFG(XXX,    YYY,  ZZZ, XXAN(1),XXAN(2),
     $                           XXAN(3),NGen, ISYM, GFac, XEFG(1,IAtm))
 20   CONTINUE
cc
 30   CONTINUE
C
C  end of angular integration
C
 50   CONTINUE
C
C  Calculate the nuclear contribution to the EFG
C  for current integration centre
C
      CALL EFG_NUC(ICntr,NAtoms,XNuc,QA,XXAN)
c
      DO 60 I=1,6
      XEFG(I,ICntr) = XEFG(I,ICntr) - XXAN(I)
 60   CONTINUE
C
      RETURN
      END
c ===================================================================
c
      SUBROUTINE EFG_NUC(ICntr,  NAtoms, XNuc,   QA,     EFGN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculates the nuclear contribution to the electric
C  field gradient at the nucleus ICntr
C
      DIMENSION XNuc(3,NAtoms),QA(NAtoms),EFGN(9)
      Parameter (Zero=0.0d0,One=1.0d0,Three=3.0d0)
C
      CALL ZeroIT(EFGN,9)
      XX = XNuc(1,ICntr)
      YY = XNuc(2,ICntr)
      ZZ = XNuc(3,ICntr)
c
      DO 10 IAtm=1,NAtoms
      If(IAtm.EQ.ICntr) GO TO 10
      X1 = (XX - XNuc(1,IAtm))
      Y1 = (YY - XNuc(2,IAtm))
      Z1 = (ZZ - XNuc(3,IAtm))
      X2 = X1**2
      Y2 = Y1**2
      Z2 = Z1**2
      r2 = X2 + Y2 + Z2
      rp = One/(sqrt(r2))**5
      qq = QA(IAtm)
      EFGN(1) = EFGN(1) + qq*(Three*X2 - r2)*rp
      EFGN(3) = EFGN(3) + qq*(Three*Y2 - r2)*rp
      EFGN(6) = EFGN(6) + qq*(Three*Z2 - r2)*rp
      EFGN(2) = EFGN(2) + qq*X1*Y1*rp
      EFGN(4) = EFGN(4) + qq*X1*Z1*rp
      EFGN(5) = EFGN(5) + qq*Y1*Z1*rp
cc      EFGN(7) = EFGN(7) + qq*X2*rp
cc      EFGN(8) = EFGN(8) + qq*Y2*rp
cc      EFGN(9) = EFGN(9) + qq*Z2*rp
 10   CONTINUE
C
      RETURN
      END
c ===================================================================
c
      SUBROUTINE NucSPINU(ICntr,  NAtoms, XNuc,   IAN,    QA,
     $                    NSym,   NGen,   ISYM,   NEqATM, IPrnt,
     $                    NRad,   NAng,   factor, r0f,    DISTN,
     $                    AIJ,    RDIST,  XXA,    WTA,    thrsh,
     $                    radii,  radwght,NBas,   NShell, BASDAT,
     $                    INX,    BL,     ExpMIN, LMax,   Ylm,
     $                    NAlpha, NBeta,  CMOA,   CMOB,   CA,
     $                    CB,     SIGNA,  SIGNB,  NScr,   Z,
     $                    XCharg, XSpin)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates spin/charge density at nucleus using density operator method
C  of Rassolov and Chipman (J.Chem.Phys. 105 (1996) 1470)
C
C  ** UNRESTRICTED OPEN SHELL **
C
C    originally written by Bing Wang based on modified JB DFT code
C    revised by JB April 2000
C
C  ARGUMENTS
C
C  ICntr   -  current atom for which spin density is being calculated
C  NAtoms  -  number of atoms
C  XNuc    -  nuclear coordinates
C  IAN     -  atomic numbers
C  QA      -  actual atomic charge (there may, e.g., be ghost atoms)
C  NSym    -  number of abelian symmetry operations
C  NGen    -  number of generators
C  ISYM    -  list of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  IPrnt   -  print flag
C  NRad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  NAng    -  angular quadrature order for numerical grid
C             (if zero, will be determined in grid generation routines)
C  factor  -  determines  grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  r0f     -  radial factor for Rassolov/Chipman function
C  DISTN   -  distance to nearest neighbour for each atom
C  AIJ     -  Becke grid parameters (see A2 in Becke article)
C  RDIST   -  inverse interatomic distance array
C  XXA     -  storage for angular quadrature grid
C  WTA     -  storage for angular quadrature weights
C  thrsh   -  threshold for neglecting contribution
C  radii   -  radial grid points
C  radwght -  radial weights
C  NBas    -  total number of basis functions
C  NShell  -  total number of shells
C  BASDAT  -  basis set data for TEXAS
C  INX     -  more basis set data for TEXAS
C  BL      -  precomputed normalization factors per primitive shell
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  LMax    -  Maximum L value for spherical harmonics
C  Ylm     -  spherical harmonics
C  NAlpha  -  number of occupied alpha MOs
C  NBeta   -  number of occupied beta MOs
C  CMOA    -  alpha MO coefficients (transposed)
C  CMOB    -  beta MO coefficients (transposed)
C  CA      -  spherically averaged MOs (alpha)
C  CB      -  spherically averaged MOs (beta)
C  SIGNA   -  how alpha MOs transform over abelian symmetry operations
C  SIGNB   -  how beta MOs transform over abelian symmetry operations
C  NScr    -  size of general scratch array
C             (5*NBas + 4*NOcc + 6)
C  Z       -  general scratch array
C
C  on exit
C
C  XCharg  -  charge density per atom
C  XSpin   -  spin density per atom
C
C  ----------------------------------------------------------------
C  references
C
C     "A multicenter numerical integration scheme for polyatomic
C      molecules"
C      A.D.Becke      J.Chem.Phys.  88 (1988) 2547
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),QA(NAtoms),ExpMIN(NShell),
     $          BL(*),BASDAT(13,*),INX(12,*),CMOA(NAlpha,NBas),
     $          CMOB(NBeta,NBas)
      DIMENSION AIJ(NAtoms,NAtoms),RDIST(NAtoms,NAtoms)
      DIMENSION DISTN(NAtoms),radii(*),radwght(*)
      DIMENSION XXA(3,1130),WTA(1130),INDEX(8),XXAN(3)
      DIMENSION Ylm(0:LMax,-LMax:LMax)
      DIMENSION CA(NAlpha,0:LMax,-LMax:LMax)
      DIMENSION CB(NBeta,0:LMax,-LMax:LMax)
      DIMENSION XCharg(NAtoms,2),XSpin(NAToms,2)
      DIMENSION Z(NScr),LF(NSym),nbf(2)
      INTEGER SIGNA(NAlpha,7),SIGNB(NBeta,7)
      Integer SFac
C
      PARAMETER (Zero=0.0d0,One=1.0d0,Two=2.0d0,Four=4.0d0)
      PARAMETER (MaxRad=400)
      Data INDEX/1,15,41,91,201,395,697,1131/
C
C
C  This method involves numerical integration over an atom-centred
C  grid, similar to current DFT implementations. However, the numerical
C  integration scheme is different. Angular integration is done using
C  spherical harmonics. Radial integration uses standard Gauss-Legendre
C  quadrature.
C
C
      PI = Four*ATAN(One)
      XX = XNuc(1,ICntr)
      YY = XNuc(2,ICntr)
      ZZ = XNuc(3,ICntr)
c
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
C
C  allocate scratch pointers
C
      i1  = 1
      i2  = i1  + NAtoms
      i3  = i2  + NAtoms
c
      iao = 1                         ! AO values over grid points
      imoA = iao + NBas               ! alpha MO values over grid points
      imoB = imoA + NAlpha            ! beta MO values over grid points
      inb = imoB + NBeta              ! index array for "non-zero" AOs
      iaox = inb + NBas               ! AO derivatives
      imoxA = iaox + 3*NBas           ! alpha MO derivatives
      imoxB = imoxA + 3*NAlpha        ! beta MO derivatives
      idxA = imoxB + 3*NBeta          ! alpha density derivatives
      idxB = idxA + 3                 ! beta density derivatives
      IEnd = idxB + 3 - 1
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'NucSPINU')
C
C ===============================================================
C  First, calculate charge/spin density at nucleus using
C  delta function method
C
      CALL DeltaFUNU(ICntr,  NBas,   NShell, NAlpha, NBeta,
     $               XNuc,   BASDAT, INX,    BL,     ExpCUT,
     $               ExpMIN, CMOA,   CMOB,   Z(iao), Z(imoA),
     $               Z(imoB),DCharg, DSpin)
c
      If(IPrnt.gt.2) Then
       write(6,*) ' Charge density from Delta Function method:',dcharg
       write(6,*) '   Spin density from Delta Function method:',dspin
      EndIf
C
C ===============================================================
C
C
C  Begin calculation for this centre
C  ---------------------------------
C
      ino = IAN(ICntr)      ! atomic number
      qai = QA(ICntr)       ! atomic charge
c
c -- set radial factor
      IF(ino.LE.10) THEN
        r0 = r0f/dble(ino)
      ELSE IF(11.LE.ino.AND.ino.LE.18) THEN
        r0 = 1.5d0*r0f/dble(ino)
      ELSE
        r0 = 4.5d0*r0f/dble(ino)
      ENDIF
      r02 = r0**2
C
C  -----------------------------------
C    RADIAL GRID
C  -----------------------------------
C  determine number of radial grid points
C
      thenum = factor*50.0d0
c
      IF(NRad.LE.0) THEN
        IF(ino.LE.2) THEN
         numR = thenum*0.6d0                     ! 30 points
        ELSE IF(3.LE.ino.AND.ino.LE.10) THEN
         numR = thenum                           ! 50 points
        ELSE IF(11.LE.ino.AND.ino.LE.18) THEN
         numR = thenum*1.4d0                     ! 70 points
        ELSE IF(19.LE.ino.AND.ino.LE.36) THEN
         numR = thenum*1.8d0                     ! 90 points
        ELSE
         numR = thenum*2.4d0                     ! 120 points
        ENDIF
      ELSE
        numR = NRad
      ENDIF
c
      If(numR.GT.MaxRad) numR=MaxRad
C
C  now get the radial grid
C
      rmin = 0.0d0
      rmax = 15.0d0*r0
      If(ino.GT.10) rmax = 45.0d0*r0
      CALL RadGRID(numR,rmin,rmax,radii,radwght)
c
      If(IPrnt.gt.2) Then
       write(6,*) ' Radial Grid for atom: ',icntr
       write(6,*) ' RMin: ',rmin,' RMax: ',rmax
       write(6,*) ' Number of radial grid points: ',numr
       If(IPrnt.gt.3) Then
        do i=1,numR
        write(6,*) i,radii(i),radwght(i)
        enddo
       EndIf
      EndIf
C
C  ===============================================================
C
      XD1A = Zero
      XD2A = Zero
      XD3A = Zero
      XD4A = Zero
      XUC = Zero
      XLC = Zero
c
      XD1B = Zero
      XD2B = Zero
      XD3B = Zero
      XD4B = Zero
      XUS = Zero
      XLS = Zero
c
      NP  = 1
      IVal = 1
      SFac = 1
      LNum = 0
c
      iang = NAng
      If(iang.EQ.0) iang=4        ! default
      IA1 = INDEX(iang)
      IA2 = INDEX(iang+1)-1
C
C  loop over radial grid points on atom IAtom
C
      DO 50 irad=1,numR
      rr = radii(irad)
      rr2 = rr**2
      rwght = radwght(irad)
c
      DenA = Zero
      DenB = Zero
      DenRAX = Zero
      DenRBX = Zero
      Call ZeroIT(CA,NAlpha*(LMax+1)*(2*LMax+1))
      Call ZeroIT(CB,NBeta*(LMax+1)*(2*LMax+1))
C
C  loop over angular grid points on atom IAtom
C
      DO 30 iang=IA1,IA2
      xap = XXA(1,iang)
      yap = XXA(2,iang)
      zap = XXA(3,iang)
      XXAN(1) = xap*radii(irad) + XX
      XXAN(2) = yap*radii(irad) + YY
      XXAN(3) = zap*radii(irad) + ZZ
C
C  get symmetry factor
C
      If(NSym.GT.0) Then
       CALL SymFAC(NAtoms, ICntr,  XX,     YY,     ZZ,
     $             XXAN(1),XXAN(2),XXAN(3),NSym,   NGen,
     $             ISYM,   NEqATM, SFac,   LNum,   LF)
       If(SFac.EQ.0) GO TO 30
      EndIf
C
C  calculate the Becke weight
C
      call WBeckeG(ICntr, XXAN(1), XXAN(2), XXAN(3),
     $             NAtoms, XNuc, RDIST, AIJ,
     $             Z(i1), Z(i2), Z(i3), Wbecke)
c
      Wtot = WTA(iang)*Wbecke*SFac
C
C  form AO and derivative values over grid points
C  and sort "non-zero" AOs
C
      CALL zeroit(Z(iao),NBas)
      CALL zeroit(Z(iaox),3*NBas)
      CALL AOGrad(NP,     NShell, NBas,   XXAN,
     $            XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $            ExpMIN, Z(iao), Z(iaox))
      CALL AOXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $             nbf,    Z(inb), Z(iao), Z(iaox),VMX)
C
C  form alpha/beta density at current grid point
C
      CALL GetDENS(NP,     NBas,   nbf,    NAlpha, CMOA,
     $             Z(iao), Z(inb), Z(imoA),DVA)
      CALL GetDENS(NP,     NBas,   nbf,    NBeta,  CMOB,
     $             Z(iao), Z(inb), Z(imoB),DVB)
C
C  prepare for numerical integration
C  neglect point if density below threshold
C
      If(DVA+DVB.LT.thrsh) GO TO 30
C
C  form density derivatives
C
      CALL GetDENDeriv(NP,    NBas,   thrsh,  nbf,    NAlpha,
     $                 CMOA,  Z(iaox),Z(inb), DVA,    Z(imoA),
     $               Z(imoxA),Z(idxA))
      CALL GetDENDeriv(NP,    NBas,   thrsh,  nbf,    NBeta,
     $                 CMOB,  Z(iaox),Z(inb), DVB,    Z(imoB),
     $               Z(imoxB),Z(idxB))
C
c -- spherically averaged electron density
      DenA = DenA + DVA*Wtot
      DenB = DenB + DVB*Wtot
c
c -- spherically averaged MOs
      call spherical_harmonics(xap,yap,zap,LMax,Ylm)
      call SphAvMO(NAlpha, LMax,   Wtot,   LNum,   LF,
     $             SIGNA,  Ylm,    Z(imoA),CA)
      call SphAvMO(NBeta,  LMax,   Wtot,   LNum,   LF,
     $             SIGNB,  Ylm,    Z(imoB),CB)
c
c -- spherically averaged MO derivatives
      call xyz_ang(xap,yap,zap,th,phi)
      CALL SphAvMOX(NAlpha,th,phi,Z(imoxA),Z(imoA),DenARX)
      CALL SphAvMOX(NBeta,th,phi,Z(imoxB),Z(imoB),DenBRX)
      DenRAX = DenRAX + DenARX*Wtot
      DenRBX = DenRBX + DenBRX*Wtot
cc
 30   CONTINUE
C
C  numerical integration
C
c -- U and D part
      XHSFAD = DenA*rwght*exp(-rr2/r02)
      XD1A = XD1A + XHSFAD*(-Two*rr2**1.5d0/r02**3)
      XD2A = XD2A + XHSFAD*(7.0d0*sqrt(rr2)/r02**2)
      XD3A = XD3A + XHSFAD*(-Two/(sqrt(rr2)*r02))
      XD4A = XD4A + DenRAX*exp(-rr2/r02)*radwght(irad)
     $                    *(Two*sqrt(rr2)/r02)
c
      XHSFBD = DenB*rwght*exp(-rr2/r02)
      XD1B = XD1B + XHSFBD*(-Two*rr2**1.5d0/r02**3)
      XD2B = XD2B + XHSFBD*(7.0d0*sqrt(rr2)/r02**2)
      XD3B = XD3B + XHSFBD*(-Two/(sqrt(rr2)*r02))
      XD4B = XD4B + DenRBX*exp(-rr2/r02)*radwght(irad)
     $                    *(Two*sqrt(rr2)/r02)
c
      XUC = XUC + (XHSFAD+XHSFBD)/rr2
      XDC = XD1A + XD2A + XD3A + XD4A + XD1B + XD2B + XD3B + XD4B
      XUS = XUS + (XHSFAD-XHSFBD)/rr2
      XDS = XD1A + XD2A + XD3A + XD4A - XD1B - XD2B - XD3B - XD4B
c
c -- L part
      CCA = Zero
      DO 40 J=1,NAlpha
      DO 40 L=1,LMax
      DO 40 M=-L,L
      CCA = CCA - L*(L+1)*CA(J,L,M)**2
 40   CONTINUE
      CCB = Zero
      DO 41 J=1,NBeta
      DO 41 L=1,LMax
      DO 41 M=-L,L
      CCB = CCB - L*(L+1)*CB(J,L,M)**2
 41   CONTINUE
c
      XLC = XLC + (CCA+CCB)*rwght*EXP(-rr2/r02)/(rr2**1.5d0)
      XLS = XLS + (CCA-CCB)*rwght*EXP(-rr2/r02)/(rr2**1.5d0)
 50   CONTINUE
C
      XUC = XUC*qai
      XLC = XLC*(Four*PI)
      RCharg = (XUC+XDC+XLC)/(Two*PI)
      XUS = XUS*qai
      XLS = XLS*(Four*PI)
      RSpin = (XUS+XDS+XLS)/(Two*PI)
cc      RSpin = (XUS+XDS+XLS)/(Two*PI*(NAlpha-NBeta))
c
      If(IPrnt.GT.2) Then
       write(6,*)' U part CHARGE integral is', XUC
       write(6,*)' D1 part CHARGE integral is', XD1A
       write(6,*)' D2 part CHARGE integral is', XD2A
       write(6,*)' D3 part CHARGE integral is', XD3A
       write(6,*)' D4 part CHARGE integral is', XD4A
       write(6,*)' D part CHARGE integral is', XDC
       write(6,*)' L part CHARGE integral is', XLC
       write(6,*)' Total charge density is',RCharg
       write(6,*)' U part SPIN integral is', XUS
       write(6,*)' D1 part SPIN integral is', XD1B
       write(6,*)' D2 part SPIN integral is', XD2B
       write(6,*)' D3 part SPIN integral is', XD3B
       write(6,*)' D4 part SPIN integral is', XD4B
       write(6,*)' D part SPIN integral is', XDS
       write(6,*)' L part SPIN integral is', XLS
       write(6,*)' Total spin density is',RSpin
      EndIf
c
      XCharg(ICntr,1) = DCharg
      XCharg(ICntr,2) = RCharg
      XSpin(ICntr,1)  = DSpin
      XSpin(ICntr,2)  = RSpin
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE DeltaFUNU(ICntr,  NBas,   NShell, NAlpha, NBeta,
     $                     XNuc,   BASDAT, INX,    BL,     ExpCUT,
     $                     ExpMIN, CMOA,   CMOB,   VAO,    VMOA,
     $                     VMOB,   DCharg, DSpin)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine calculates the charge/spin density at the nucleus
C  by the delta function method
C
      Dimension XNuc(3,*),BL(*),ExpMIN(NShell),BASDAT(13,*),
     $          INX(12,*),CMOA(NAlpha,NBas),CMOB(NBeta,NBas)
      Dimension VAO(NBas),VMOA(NAlpha),VMOB(NBeta)
C
C  get basis function values at the nucleus
C
      CALL ZeroIT(VAO,NBas)
      CALL AOVal(1,      NShell, NBas,   XNuc(1,ICntr),
     $           XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $           ExpMIN, VAO)
C
C  now evaluate the charge density
C
      Call ZeroIT(VMOA,NAlpha)
      Call ZeroIT(VMOB,NBeta)
      DO 20 I=1,NBas
      Val = VAO(I)
      DO 10 J=1,NAlpha
      VMOA(J) = VMOA(J) + CMOA(J,I)*Val
 10   CONTINUE
      DO 11 J=1,NBeta
      VMOB(J) = VMOB(J) + CMOB(J,I)*Val
 11   CONTINUE
 20   CONTINUE
c
      DenA = DDOT(NAlpha,VMOA,1,VMOA,1)
      DenB = DDOT(NBeta,VMOB,1,VMOB,1)
      DCharg = DenA + DenB
      DSpin  = DenA - DenB
cc      DSpin = (DenA-DenB)/(NAlpha-NBeta)
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE FindMOSym(NBas,   NOcc,   NSym,   ISYM,   IFP,
     $                     CMO,    VNorm,  VSym,   SIGN,   full)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines whether or not an MO changes sign under the
C  (Abelian) symmetry operations of the point group
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  NOcc    -  number of MOs
C  NSym    -  number of abelian symmetry operations
C  ISYM    -  list of operations
C              1 - reflection in YZ plane (X changes sign)
C              2 - reflection in XZ plane (Y changes sign)
C              3 - reflection in XY plane (Z changes sign)
C              4 - rotation about Z axis  (X,Y change sign)
C              5 - rotation about Y axis  (X,Z change sign)
C              6 - rotation about X axis  (Y,Z change sign)
C              7 - inversion             (X,Y,Z change sign)
C  IFP     -  basis function symmetry pairs under
C             symmetry operations
C  CMO     -  molecular orbital coefficients
C  VNorm   -  scratch vector for "normalized" MO
C  VSym    -  scratch vector for symmetry image
C
C  on exit
C
C  SIGN    -  sign of MO under symmetry operation
C              SIGN(I,J)=1   sign of MO I unchanged under operation J
C              SIGN(I,J)=-1     ditto  sign changed
C  full    -  logical flag for full symmetry analysis
C
C
      DIMENSION CMO(NBas,NOcc),VNorm(NBas),VSym(NBas)
      INTEGER ISYM(NSym),IFP(7,NBas),SIGN(NOcc,7)
      Logical full
C
      Parameter (Zero=0.0d0,One=1.0d0,TolZero=1.0d-5)
C
C
C  Assume the worst
C
      full = .FALSE.
C
      DO 30 I=1,NOcc
C
C  "normalize" the current MO
C
      Ovlp = SProd(NBas,CMO(1,I),CMO(1,I))
      Ovlp = One/SQRT(Ovlp)
      DO 25 J=1,NBas
      VNorm(J) = CMO(J,I)*Ovlp
 25   CONTINUE
c
      DO 20 J=1,NSym
      JJ = ISYM(J)
      DO 10 K=1,NBas
C
C  generate the symmetry image of basis function K under
C  symmetry operation J
C
      IS = IFP(J,K)
      IT = Abs(IS)
      If(IS.GT.0) Then
       VSym(IT) = VNorm(K)
      Else If(IS.LT.0) Then
       VSym(IT) = -VNorm(K)
      Else
       VSym(IT) = Zero
      EndIf
 10   CONTINUE
C
C  determine the overlap between MO I and its symmetry image
C
      Ovlp = SProd(NBas,VNorm,VSym)
cc      write(6,*) ' MO:',i,' Sym:',j,' Ovlp is ',ovlp
C
C  now get sign
C
      If(Abs(Ovlp-One).LT.TolZero) Then
       SIGN(I,JJ) = 1
      Else If(Abs(Ovlp+One).LT.TolZero) Then
       SIGN(I,JJ) = -1
      Else
c
c -- cannot determine sign
       RETURN
      EndIf
 20   CONTINUE
 30   CONTINUE
C
C  completed full analysis
C
      full = .TRUE.
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE TempSYM(NAtoms,NSym,NGen,NQ,IUNQ)
      IMPLICIT INTEGER(A-Z)
      DIMENSION IUNQ(NAtoms)
C
C  Symmetry not working properly
C  Switch off until fixed
C
      NSym = 0
      NGen = 0
      NQ = NAtoms
      DO I=1,NQ
      IUNQ(I) = I
      EndDO
C
      RETURN
      END
