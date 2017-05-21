c ==============================================================
c  REACTION PATH MODULE          JB   Sep. 2000
c ==============================================================
c
      subroutine prepath(inp,scyc)
      implicit real*8(a-h,o-z)
      character*256 chopval
      integer scyc
c
c  reads the PATH line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=7)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character*20 coord,cdum
      character*256 jobname
c
      parameter (IUnit=1)
      Common /job/jobname,lenJ
c
      data options/'coor','dmax','sign','iter','dtol','prin',
     $             'path'/
      data ioptyp/21,11,1,1,11,1,0/
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c ...........................................................
      If(scyc.gt.0) RETURN        ! only read first time
c ...........................................................
c
      call setival('isumscf',1)      ! switch off SCF summary print
c
c -- coordinate type
      if(ifound(1).eq.1) then
        coord = chopval(1)
      else
        coord = 'mwgt'             ! default is mass-weighted Cartesians
      endif
c
c -- maximum step size
      if(ifound(2).eq.1) then
        dmax = ropval(1,2)
      else
        dmax = 0.15d0              ! default step length
      endif
c
c -- direction of first step (from TS eigenvector)
      if(ifound(3).eq.1) then
        isign = iopval(1,3)
      else
        isign = +1
      endif
c
c -- number of steps
      if(ifound(4).eq.1) then
        iter = iopval(1,4)
      else
        iter = 20                  ! default number of steps
      endif
c
c -- step size for convergence
      if(ifound(5).eq.1) then
        dtol = ropval(1,5)
      else
        dtol = 0.005d0             ! default step size for convergence
      endif
c
c -- print flag
      if(ifound(6).eq.1) then
        IPRNT = iopval(1,6)
      else
        IPRNT = 2                  ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call wrpath(IUnit,coord,dmax,isign,iter,dtol)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c .............................................................................
c
      SUBROUTINE PATH(NMem,Z,back,Cnvgd)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** REACTION PATH PROGRAM **
C
C  This Program follows a "reaction path" downhill starting from a
C  transition state. The first step is in the direction of the TS
C  eigenmode; thereafter the initial step is in the direction of the
C  gradient. After the first step (which in general does not follow
C  the reaction path), a line search is carried out in the direction
C  of the gradient bisector to get back on the "true" reaction path.
C
C  The reaction path can be defined in Cartesian coordinates, Z-matrix
C  coordinates or Mass-Weighted Cartesians (the default). In the latter
C  case, the path is equivalent to the Intrinsic Reaction Coordinate
C  (IRC) defined by Fukui.
C
C
C  References
C  ----------
C
C  "A Formulation of the Reaction Coordinate"
C   K.Fukui  J.Phys.Chem.  74 (1970) 4161
C
C  "The Intrinsic Reaction Coordinate: An Ab Initio Calculation for
C   HCN --> HNC and H- + CH4 --> CH4 + H-"
C   K.Ishida, K.Morokuma and A.Komornicki,  J.Chem.Phys. 66 (1977) 2153
C
C  "The Intrinsic Reaction Coordinate and the Rotational Barrier
C   in Silaethylene"
C   M.W.Schmidt, M.S.Gordon and M.Dupuis  J.Am.Chem.Soc. 107 (1985) 2585
C
C  ----------------------------------------------------------------------
      CHARACTER cdum*20,jobname*256,coord*4
      LOGICAL back,Cnvgd,mwght
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/              !  unit number for checkpoint I/O
      Common /job/jobname,lenJ
C
      IOut = ioutfil('iout')
      mwght = .FALSE.
C
C --------------------------------------------------------------
C  ** CHECK LICENSE **
C
cc      CALL ChkLicense
C --------------------------------------------------------------
C
C  Read from the <control> file
C    total number of atomic centres, including dummies
C    total number of molecules
C    data for the reaction path
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      call rdcntrl(IUnit,5,'$nmol',1,NMol,rdum,cdum)
      call rdpath(IUnit,coord,DMax,ISign,MaxCyc,DTol)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Read from the <sym> file
C    number of real atoms
C    number of symmetry operations
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,rdum,cdum)
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  For a Z-matrix reaction path, determine the number of centers
C  in the Z matrix
C
      IF(coord(1:4).EQ.'zmat') THEN
       CALL ScanZMAT(NZ,IErr)
       If(IErr.EQ.-1) Call nerror(1,'REACTION PATH module',
     $    'Z-Matrix Path Requested but there is NO Z Matrix!!',0,0)
      ELSE
       NZ = 0
       If(coord(1:4).EQ.'mwgt') mwght=.TRUE.
      ENDIF
C
      NAT3 = 3*NAtoms
C
C  Allocate the maximum scratch storage
C
      NScr = MAX(3*NAT3**2,36*NZ**2 + 12*NZ)
C
C  Now get the memory
C
      IMem = 5*NAtom + NMol+1 + 9*NTrans + NAtoms*NTrans + 19*NAtoms
     $        + 9*NAtoms**2 + 16*NZ + 9*NZ**2 + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,4,'PATH')
C
C  Allocate memory pointers
C
      IXS = iptr                     !  TS coordinates
      IXC =  IXS + 3*NAtoms          !  current coordinates
      IMOL = IXC + 3*NAtom           !  pointer to start/end of molecules
      IXCG = IMOL + NMol+1           !  atomic charges
      IXCM = IXCG + NAtom            !  atomic masses
      IXO = IXCM + NAtom             !  old coordinates
      IGC = IXO + 3*NAtoms           !  current gradient
      IGO = IGC + 3*NAtoms           !  old gradient
      IGS = IGO + 3*NAtoms           !  gradient bisector
      IHS = IGS + 3*NAtoms           !  Hessian matrix
      ID  = IHS + NAT3*NAT3          !  current displacement
      ITN = ID  + 3*NAtoms           !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*NTrans           !  list of atomic equivalences
      IUQ = INQ + NAtoms*NTrans      !  list of symmetry-unique atoms
      IZO = IUQ + NAtoms             !  Z-matrix parameter values
      IZS = IZO + 3*NZ               !  Z-matrix connectivity
      IZG = IZS + 4*NZ               !  parameter optimization array
      IXT = IZG + 3*NZ               !  values of internal coordinates
      IGT = IXT + 3*NZ               !  internal gradient
      IHT = IGT + 3*NZ               !  internal Hessian
      IEnd = IHT + 9*NZ**2
c
      IScr = IEnd + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(NMem,IEnd,4,'SCAN')
C
C
c memory status
c -- assign memory for high water mark (old TEXAS)
      call getmem(IEnd,lastx)
C
C  ----------------------------------------------------------------------
C
      CALL PATHMAIN(NAtom,   NAtoms,  NMol,    Z(IXC),  Z(IMOL),
     $              Z(IXCG), Z(IXCM), mwght,   DMax,    ISign,
     $              MaxCyc,  DTol,    Z(IXS),  Z(IXO),  Z(IGC),
     $              Z(IGO),  Z(IGS),  Z(IHS),  Z(ID),   NTrans,
     $              Z(ITN),  Z(INQ),  Z(IUQ),  NZ,      Z(IZO),
     $              Z(IZS),  Z(IZG),  Z(IXT),  Z(IGT),  Z(IHT),
     $              NScr,    Z(IScr), back,    Cnvgd)
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
      write(IOut,1100) IEnd,mxmem,memtot
 1100 format(/,' Reaction Path memory status:',/,
     $         ' memory needed=',i15,' high water=',i15,/,
     $         ' total available memory=',i15)
c-------------------------------------------
C
C  Exit procedure
C
      RETURN
      END
c ..............................................................................
c
      SUBROUTINE PATHMAIN(NAtom,  NAtoms, NMol,   XC,     IMOL,
     $                    XCharg, XMass,  mwght,  DMax,   ISign,
     $                    MaxCyc, DTol,   XOrig,  XOld,   GC,
     $                    GOld,   GS,     HESS,   D,      NTrans,
     $                    TRANS,  NEqATM, IUNQ,   NZ,     IGEO,
     $                    GEO,    IG,     XINT,   GINT,   HINT,
     $                    NMem,   Z,      back,   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for Reaction Path module
C  PATHMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION XC(3*NAtom),XCharg(NAtom),XMass(NAtom),
     $          XOld(3,NAtoms),XOrig(3,NAtoms)
      DIMENSION GC(3*NAtoms),GOld(3*NAtoms),GS(3*NAtoms),
     $          HESS(9*NAtoms**2),D(3*NAtoms)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans),
     $          IUNQ(NAtoms),IMOL(NAtoms),RM(3,3)
      DIMENSION GEO(NZ,3),IGEO(NZ,4),IG(3*NZ)
      DIMENSION XINT(3*NZ),GINT(3*NZ),HINT(9*NZ**2)
c ............................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtom),ZSymb(NZ),VARNAM(3*NZ)
      CHARACTER*8 CZ(14*NZ+3*NAtoms)
c ............................................................
      CHARACTER GROUP*4,cdum*20,jobname*256
      LOGICAL Symflag,mwght,parab,found,back,Cnvgd
C
      DIMENSION Z(NMem)
C
      Common /job/jobname,lenJ
      Data IUnit/1/
C
C
      IOut = ioutfil('iout')
      ICond = ioutfil('icond')
c
      NAT3 = 3*NAtoms
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
C  Read from the <control> file
C    current energy (if available)
C    print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call fdcntrl(IUnit,7,'$energy',idum)
      If(idum.EQ.0) call rdcntrl(IUnit,7,'$energy',2,idum,EC,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  read initial coordinates from <coord> file
C  **WARNING** There may not be a <coord> file available for
C    a Z-Matrix reaction path
C
      INQUIRE(FILE=jobname(1:lenJ)//'.coord',EXIST=found)
      If(found) Then
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $        FORM='FORMATTED',STATUS='OLD')
        CALL RdCoordF(IUnit,  NAtoms, AtSymb, XC,     NMol,
     $                IMOL,   XCharg, XMass)
        CLOSE (UNIT=IUnit,STATUS='KEEP')
      EndIf
C
C  Attempt to read <path> file
C  If no <path> file is available, this is the first
C  step of a new reaction path
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.path',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(40) NStep,NVar,DMax,ISign,MaxCyc,DTol,IPRNT,IStep,LStep,
     $         E0,d0,E1,d1,E2,d2,dm,parab,XOrig,XOld,GOld,GS,HESS
      If(NZ.GT.0) READ(40) IGEO,GEO,IG,XINT,ZSymb,VARNAM
      CLOSE (UNIT=40,STATUS='KEEP')
      GO TO 20
C
C  Nothing on file
C  Initialize
C
 10   CONTINUE
      NStep = 0
      IStep = 0
      LStep = 0
      EC = 0.0d0
      Cnvgd = .False.
      NVar = NAT3
C
      IF(NZ.GT.0) THEN
C
C  Read Z-Matrix data
C  allocate scratch memory pointers
C
       ismb = 1
       is1 = ismb + NZ
       is2 = is1 + 7*NZ
       IEnd = is2 + 6*NZ - 1
c
       CALL RdZMAT(NZ,  CZ(ismb), CZ(is1),CZ(is2),Z(1),
     $             GEO,   IGEO,   IG,     VARNAM, NVar)
C
C  get atomic symbols from Z-matrix symbols
C  and assign dummy atoms
C
       CALL GetAtSym(NZ,CZ(ismb),ZSymb)
      ENDIF
C
C  At the start of a new reaction path search there should
C  be a Hessian we can use
C
      CALL RdHESS(jobname(1:lenJ)//'.hess',lenJ+5,NAT3,IPRNT,
     $             HESS,IHess)
      If(IHess.EQ.-1) Call nerror(2,'REACTION PATH module',
     $    'No Hessian Matrix Available at Start of Path',0,0)
C
 20   CONTINUE
C
C  If appropriate, read gradient from <grad> file
C
      If( (IStep.EQ.2.AND.NStep.GT.1).OR.IStep.EQ.4 ) Then
       CALL RdGRAD(NAtoms,GC,'save')
       CALL VScal(NAT3,-1.0d0,GC)            ! file contains forces
      EndIf
C
C  Do the Reaction Path
C
C ---------------------------------------------------------------------
      CALL RPATH(NStep,  NAtoms, AtSymb, XC,     XMass,
     $           mwght,  DMax,   ISign,  MaxCyc, DTol,
     $           XOrig,  XOld,   EC,     GC,     GOld,
     $           GS,     HESS,   D,      NTrans, TRANS,
     $           NEqATM, IPRNT,  NZ,     NVar,   IGEO,
     $           GEO,    IG,     ZSymb,  RM,     XINT,
     $           GINT,   HINT,   IStep,  LStep,  E0,
     $           d0,     E1,     d1,     E2,     d2,
     $           dm,     parab,  NMem,   Z,      Cnvgd)
C ---------------------------------------------------------------------
C
      If(Cnvgd) Then
c -- delete <path> file
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.path',FORM='FORMATTED')
       CLOSE (UNIT=40,STATUS='DELETE')
       RETURN
      EndIf
C
C  Another step needed
C
C  Write <path> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.path',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(40) NStep,NVar,DMax,ISign,MaxCyc,DTol,IPRNT,IStep,LStep,
     $          E0,d0,E1,d1,E2,d2,dm,parab,XOrig,XOld,GOld,GS,HESS
      If(NZ.GT.0) WRITE(40) IGEO,GEO,IG,XINT,ZSymb,VARNAM
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  Write necessary data to <coord>/<zmat> file
C
      IF(NZ.GT.0) THEN
       If(IStep.NE.2.AND.IStep.NE.4)
     $    CALL WrZMAT(NZ,GEO,IG,VARNAM,.true.)       ! <coord> file deleted
      ELSE
c -- check if there are dummy atoms whose coordinates need to be written
       If(NAtom.GT.NAtoms) Then
         OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $         FORM='FORMATTED',STATUS='OLD')
         CALL RdCoordF(IUnit,  NAtom,  AtSymb, Z,      NMol,
     $                 IMOL,   XCharg, XMass)
         CLOSE (UNIT=IUnit,STATUS='KEEP')
         CALL CpyVEC(3*(NAtom-NAtoms),Z(3*NAtoms+1),XC(3*NAtoms+1))
       EndIF
c
       CALL WrCoord(NAtom,  AtSymb, XC,     NMol,   IMOL,
     $              XCharg, XMass)
      ENDIF
C
C  ...............................................
C     Maximum Reaction Path Steps Reached
C
      IF(NStep.GT.MaxCYC) THEN
       WRITE(IOut,2000)
       WRITE(ICond,2000)
       Cnvgd = .TRUE.
      ENDIF
C  ...............................................
C
C  Set logical flag back
C  If linesearch, energy only
C
      back = (2*((LStep+1)/2)).NE.LStep      ! LStep odd
      back = (back).OR.(NStep.EQ.1.AND.IStep.EQ.2)
C
      RETURN
c
 2000 FORMAT(//,' *******************************************',/,
     $          ' **  MAXIMUM REACTION PATH STEPS REACHED  **',/,
     $          ' *******************************************')
c
      END
c ..............................................................................
c
      SUBROUTINE RPATH(NStep,  NAtoms, AtSymb, XC,     XMass,
     $                 mwght,  DMax,   ISign,  MaxCyc, DTol,
     $                 XOrig,  XOld,   EC,     GC,     GOld,
     $                 GS,     HESS,   D,      NTrans, TRANS,
     $                 NEqATM, IPRNT,  NZ,     NVar,   IGEO,
     $                 GEO,    IG,     ZSymb,  RM,     XINT,
     $                 GINT,   HINT,   IStep,  LStep,  E0,
     $                 d0,     E1,     d1,     E2,     d2,
     $                 dm,     parab,  NMem,   Z,      Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out a Reaction Path Search
C  This is done EITHER in:
C   (1) Cartesian coordinates
C   (2) Mass-Weighted Cartesians     (default)
C   (3) Z-Matrix coordinates
C
C  ARGUMENTS
C
C  NStep    -  step number
C              i.e. number of times this routine has been called
C              initial value should be -1 (for coordinate generation)
C              incremented by 1 on exit
C  NAtoms   -  number of real atoms
C  AtSymb   -  atomic symbols (used for nice printout)
C  XC       -  current Cartesian coordinates
C              contains new coordinates on exit
C  XMass    -  atomic masses
C  mwght    -  logical flag for use of mass-weighted Cartesians
C  DMax     -  maximum stepsize during Reaction Path search
C  ISign    -  direction of initial step from TS
C  MaxCyc   -  maximum number of points to find on Reaction Path
C  DTol     -  convergence criterion on step size
C  XOld     -  old coordinates
C  EC       -  current energy
C  GC       -  current Cartesian gradient
C  GOld     -  old gradient
C  GS       -  gradient bisector
C  HESS     -  Cartesian Hessian matrix at TS
C  D        -  current step
C  NTrans   -  number of symmetry operations
C  TRANS    -  symmetry operations as 3x3 transformation matrices
C  NEqATM   -  list of atomic equivalences under symmetry operations
C  IPRNT    -  print flag
C  NZ       -  number of Z-matrix centers for a Z-matrix scan
C              (set to zero for scan in Cartesians)
C  NVar     -  number of Z-matrix variables
C  IGEO     -  Z-matrix connectivity
C  GEO      -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IG       -  array determining what to do with Z-matrix parameter
C                  0 - optimize it
C                  J - assign same value as previous (Jth) variable
C                 -J - assign same value, opposite sign
C               1000 - fixed
C  ZSymb    -  Z-Matrix atomic symbols
C  RM       -  Rotation matrix
C  XINT     -  array to hold Z-Matrix values
C  GINT     -  internal gradient
C  HINT     -  internal Hessian
C  IStep    -  current/previous step time
C               0 - first step of a new search
C               1 - calculate gradient for new step
C               2 - just got gradient
C               3 - calculate gradient for bisector step
C               4 - just got gradient
C               5 - energy-only line search
C  LStep    -  line search step
C  E0       -  energy at start of new Reaction Path step or
C              start of line search
C  d0       -  line search parameter corresponding to E0
C  E1       -  energy after first gradient or line search step
C  d1       -  line search parameter corresponding to E1
C  E2       -  energy during line search
C  d2       -  line search parameter corresponding to E2
C  dm       -  current estimate of line search minimum
C  parab    -  logical flag for parabolic fitting in line search
C  NMem     -  amount of available scratch storage
C  Z        -  scratch array
C  Cnvgd    -  Logical flag (set to TRUE for last step of scan)
C
C
      DIMENSION XC(3*NAtoms),GC(3*NAtoms),GOld(3*NAtoms),
     $          GS(3*NAtoms),HESS(3*NAtoms,3*NAtoms),
     $          D(3*NAtoms),XOld(3*NAtoms),XMass(NAtoms),
     $          XOrig(3*NAtoms)
      DIMENSION TRANS(9,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION IGEO(NZ,4),GEO(NZ,3),IG(3*NZ)
      DIMENSION XINT(3*NZ),GINT(3*NZ),HINT(3*NZ,3*NZ)
      CHARACTER*8 AtSymb(NAtoms),ZSymb(NZ)
      LOGICAL mwght,parab,Step,Cnvgd
c
      Dimension Z(NMem)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
      PARAMETER (small=1.0d-8,delta=0.025d0)
C
C
      IOut = ioutfil('iout')
      ICond = ioutfil('icond')
c
      NAT3 = 3*NAtoms
      NVib = MAX(1,3*NZ-6)
      If(NZ.EQ.0) NVar = NAT3
      IDB = 0
      Step = .True.
C
C  ....................................................................
C    SYMMETRY SECTION
C    make sure the geometry and gradient conform fully
C    to the molecular symmetry
C
C  **WARNING** If Z-matrix path and IStep=2 or 4 there are NO
C    coordinates available, so skip
C
      If(NZ.GT.0.AND.(IStep.EQ.2.OR.IStep.EQ.4)) GO TO 5
        CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z,
     $              small,  XC)
 5    CONTINUE
c
      If(IStep.EQ.2.OR.IStep.EQ.4)
     $  CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z,
     $              small,  GC)
c
      IF(IStep.EQ.0.AND.NStep.EQ.0) THEN
       IV1 = 1
       IV2 = IV1 + NAT3**2
       IV3 = IV2 + NAT3**2
c
       CALL CpyVEC(NAT3*NAT3,HESS,Z(IV1))
       CALL SymHES(NAtoms, NTrans, NEqATM, TRANS,  Z(IV1),
     $             Z(IV2), Z(IV3), small,  HESS)
      ENDIF
C
C .....................................................................
C    Z-MATRIX COORDINATES
C    transform gradient, if available, to Z-Matrix variables
C
      IF(NZ.GT.0) THEN
        If(NStep.EQ.0.AND.IStep.EQ.0) Then
         IHess = 0
         CALL ZeroIT(GC,3*NAtoms)     ! no gradient available
        Else If(IStep.EQ.2.OR.IStep.EQ.4) Then
         IHess = 1
        Else
         CALL CpyVEC(NVar,XINT,XC)
         GO TO 100                    ! no transformation needed
        EndIf
c
        IXZ = 1
        INT = IXZ + 3*NZ
        IEnd = INT +  NVib
        IMem = NMem - IEnd
c
        CALL GetZINT(NZ,     NAtoms, NVib,   NVar,   IPRNT,
     $               IHess,  ZSymb,  GEO,    IGEO,   IG,
     $               IDB,    GC,     HESS,   RM,     Z(IXZ),
     $               Z(INT), XINT,   GINT,   HINT,   IMem,
     $               Z(IEnd) )
C
C  we will use Cartesian arrays for all manipulations
C
        CALL CpyVEC(NVar,XINT,XC)
        If(IStep.EQ.0) Then
          CALL CpyVEC(NVar**2,HINT,HESS)
        Else
          CALL CpyVEC(NVar,GINT,GC)
        EndIf
      ENDIF
C .....................................................................
C
C
C  REACTION PATH SEARCH
C  --------------------
C
      If(IPRNT.GT.0) Then
        If(NZ.GT.0) Then
          WRITE(IOut,1000)
        Else If(mwght) Then
          WRITE(IOut,1010)
        Else
          WRITE(IOut,1020)
        EndIf
      EndIf
C
 100  CONTINUE
C
C  Determine what action to take
C  This is controlled mainly by IStep
C
      IF(IStep.EQ.0) THEN
C
C  First step of a new Reaction Path Search
C  Make copy of starting geometry, increment counters and exit
C
        NStep = NStep + 1
        IStep = IStep + 1
        If(NZ.GT.0) Then
          CALL BackToZ(NZ,     NVib,   NVar,   IG,     XC,
     $                 XOrig,  GEO)
        Else
          CALL CpyVEC(NVar,XC,XOrig)
        EndIf
        If(IPRNT.GT.0) WRITE(IOut,1050)
        Step = .False.
        GO TO 200
C
      ELSE IF(IStep.EQ.1) THEN
C
C  We need to compute a gradient for the first step
C  simply increment IStep and exit UNLESS we are on the
C  very first step in which case use the Hessian mode
C
        IF(NStep.GT.1) THEN
          If(IPRNT.GT.0) WRITE(IOut,1100)
          IStep = IStep + 1
          Step = .False.
          GO TO 200
        ELSE
C
C  Very first step
C  ---------------
C  It is assumed that we have a geometry corresponding to the
C  transition state with an exact Hessian which will be diagonalized
C  to get the initial search direction
C
        If(IPRNT.GT.0) WRITE(IOut,1200)
        If(NZ.GT.0) Then
C
C  convert new geometry (internal coordinates on XC)
C  back to Z-Matrix internals
C
          CALL BackToZ(NZ,     NVib,   NVar,   IG,     XC,
     $                 XINT,   GEO)
c
          CALL PrntZPATH(IOut,   NStep,  NZ,     ZSymb,  GEO,
     $                   IGEO,   NVib,   XOrig,  XINT,   EC)
          CALL PrntZPATH(ICond,  NStep,  NZ,     ZSymb,  GEO,
     $                   IGEO,   NVib,   XOrig,  XINT,   EC)
        Else
          CALL PrntPATH(IOut,   NStep,  NAtoms, mwght,  AtSymb,
     $                  XOrig,  XC,     EC)
          CALL PrntPATH(ICond,  NStep,  NAtoms, mwght,  AtSymb,
     $                  XOrig,  XC,     EC)
        EndIf
c
        If(mwght) Then
C
C  mass-weight the Hessian
C
          CALL HessWT(NAtoms,XMass,HESS)
c
          If(IPRNT.GT.3) Then
           WRITE(IOut,1210)
           CALL PrntMAT(NVar,NVar,NVar,HESS)
          EndIf
C
C  project out translations/rotations
C
          I1 = 1
          I2 = I1 + NVar**2
c
          CALL ProjTRM(NAtoms, IPRNT,  XC,     XMass,  Z(I1),
     $                 Z(I2),  HESS)
c
          If(IPRNT.GT.3) Then
           WRITE(IOut,1220)
           CALL PrntMAT(NVar,NVar,NVar,HESS)
          EndIf
c
        Else If(NZ.GT.0) Then
C
C  Z-Matrix Hessian already available
C
          If(IPRNT.GT.3) Then
           WRITE(IOut,1230)
           CALL PrntMAT(NVar,NVar,NVar,HESS)
          EndIf
c
        Else
C
C  Cartesian coordinates
C  project out translations/rotations
C
          I1 = 1
          I2 = I1 + NVar**2
c
          CALL ProjTR(NAtoms, XC,     2,      IPRNT,  Z(I1),
     $                Z(I2),  HESS)
c
          If(IPRNT.GT.3) Then
           WRITE(IOut,1220)
           CALL PrntMAT(NVar,NVar,NVar,HESS)
          EndIf
        EndIf
C
C  Diagonalize Hessian to get initial step direction
C  (GC used as scratch, GS for eigenvalues)
C
        CALL DIAGMAT(HESS,NVar,Z,GC,GS,IErr)
c
        If(IPRNT.GT.2) Then
         WRITE(6,1240)
         WRITE(6,1250) (GS(I),I=1,NVar)
        EndIf
C
C  partial check
C  Lowest eigenvalue should be negative, next positive
C
        If(GS(1).GT.Zero.AND.GS(2).LT.Zero) Then
           Call nerror(3,'REACTION PATH module',
     $      'Starting Geometry Does NOT Correspond to TS',0,0)
        EndIf
C
C  Prepare for search
C  copy lowest Hessian mode into D
C
        CALL CpyVEC(NVar,HESS,D)
        If(mwght) CALL MassWT(-1,NAtoms,XMass,D)  ! not sure about this  JB
        CALL CpyVEC(NVar,XC,XOld)
cc        write(6,*) ' Hessian mode is:'
cc        do i=1,nvar
cc        write(6,*) i,d(i)
cc        enddo
C
C  calculate step and copy -D into GOld for use with gradient bisector
C
        DO I=1,NVar
        GOld(I) = -ISign*D(I)
        D(I) = ISign*DMax*D(I)
        EndDO
C
C  increment IStep
C
        IStep = IStep + 1
        ENDIF
C
C ----------------
C
      ELSE IF(IStep.EQ.2) THEN
C
C  We have a new gradient so take the first step
C
        If(NStep.EQ.1) Then
C
C  already taken Hessian step last time
C  increment IStep and exit
C
          If(IPRNT.GT.0) WRITE(IOut,1300)
          IStep = IStep + 1
          Step = .False.
          GO TO 200
        EndIf
c
        If(IPRNT.GT.0) WRITE(IOut,1310)
c
        If(mwght) CALL MassWT(-1,NAtoms,XMass,GC)
        G0 = SQRT(SProd(NVar,GC,GC))
cc        write(6,*) ' gradient norm is: ',g0
C
c ....................................................
C  if gradient step too small stop
C
        IF(G0.LT.DTol) THEN
          If(IPRNT.GT.0) WRITE(IOut,1330)
          Cnvgd = .TRUE.
          Step = .False.
          GO TO 200
        ENDIF
c ....................................................
C
C  set scale factor for steepest descent
C
        Skal = DMax/G0
        If(IPRNT.GT.0) WRITE(IOut,1320) DMax
C
C  put step into D
C
        DO I=1,NVar
        D(I) = -Skal*GC(I)
        EndDO
c
        CALL CpyVEC(NVar,GC,GOld)
        CALL CpyVEC(NVar,XC,XOld)
        IStep = IStep + 1
C
      ELSE IF(IStep.EQ.3) THEN
C
C  We need to compute a gradient for the bisector step
C  inilialize E), increment IStep and exit
C
        E0 = EC
c
        If(IPRNT.GT.0) WRITE(IOut,1100)
        IStep = IStep + 1
        Step = .False.
        GO TO 200
C
      ELSE IF(IStep.EQ.4) THEN
C
C  Calculate the gradient bisector step and initialize
C  the line search
C
        If(IPRNT.GT.0) WRITE(IOut,1400)
c
        If(mwght) CALL MassWT(-1,NAtoms,XMass,GC)
        G1 = SQRT(SProd(NVar,GC,GC))
        G0 = SQRT(SProd(NVar,GOld,GOld))
cc        write(6,*) ' G1:',g1,' G0:',g0
c
        DO I=1,NVar
        GS(I) = -GC(I)/G1 + GOld(I)/G0
        EndDO
C
C  normalize
C
        Skal = SProd(NVar,GS,GS)
        Skal = One/SQRT(Skal)
        CALL VScal(NVar,Skal,GS)
C
C  project gradient at origin along search direction to
C  get initial step length
C
        d1 = -SProd(NVar,GC,GS)
cc        write(6,*) ' gradient projected step is: ',d1
        if(d1.lt.zero) call nerror(1,'rpath',
     $              'd1 .lt. zero',0,0)
c
        If(IPRNT.GT.0) Then
         WRITE(IOut,1500)
         WRITE(IOut,1510) d1
        EndIf
c
        If(d1.GT.2*DMax) Then
          If(IPRNT.GT.0) WRITE(IOut,1520) DMax
          d1 = DMax
        EndIf
C
C  put step into d
C
        DO I=1,NVar
        D(I) = d1*GS(I)
        EndDO
c
        CALL CpyVEC(NVar,XC,XOld)
        CALL CpyVEC(NVar,GC,GOld)
c
        d0 = Zero
        LStep = 0
        IStep = IStep + 1
C
      ELSE IF(IStep.EQ.5) THEN
C
C  Line search step
C  ----------------
C  Line search involves an energy loop only
C  If the step counter, LStep, is odd, simply exit
C  to compute another energy
C
        LL = 2*((LStep+1)/2)
cc        write(6,*) ' In Line Search  LStep:',lstep,' LL:',ll
c
        IF(LL.NE.LStep) THEN       ! L is odd
          If(IPRNT.GT.0) WRITE(IOut,1600)
          LStep = LStep + 1
          RETURN
        ENDIF
C
C  LStep is even
C
        If(IPrnt.GT.0) WRITE(IOut,1700) LL+1
c
        IF(LStep.EQ.0) THEN
C
          E1 = EC
c
          If(d1.NE.DMax) Then
C
C  just taken gradient projector step
C  fit a parabola using two energies and the projected gradient
C
            CALL FitParab(E0,d0,E1,d1,E2,d2,.true.,dm)
            parab = .True.
            If(IPRNT.GT.0) WRITE(IOut,1710) dm
            if(dm.lt.zero) call nerror(1,'rpath',
     $              'dm .lt. zero',0,0)
C
C  dm is the estimated minimum
C  if dm > DMax, reject
C
            If(dm.GT.DMax) Then
             d2 = d1 + d1
             parab = .False.
             If(IPRNT.GT.0) WRITE(IOut,1720) d2
            Else
             d2 = dm
            EndIf
          Else
C
C  gradient projector step rejected
C  take another line search step
C
            If(E1.GT.E0) Then
             d2 = Half*d1
            Else
             d2 = d1 + d1
            EndIf
            parab = .False.
            If(IPRNT.GT.0) WRITE(IOut,1730) d2
          EndIf
C
C  put step into D
C
          DO I=1,NVar
          D(I) = d2*GS(I)
          EndDO
c
          LStep = LStep + 1
cc
        Else If(LStep.EQ.2) Then
C
C  just taken line search step
C
          E2 = EC
c
cc          write(6,*) ' LStep:',lstep,' dm:',dm
cc          write(6,*) ' d0:',d0,' E0:',e0
cc          write(6,*) ' d1:',d1,' E1:',e1
cc          write(6,*) ' d2:',d2,' E2:',e2
          If(parab.AND.(E2.LT.E1.AND.E2.LT.E0)) Then
C
C  terminate line search with d2 as minimum
C
            If(IPRNT.GT.0) WRITE(IOut,1800)
            E0 = E2
            NStep = NStep + 1
            IStep = 1
            LStep = 0
            If(NZ.GT.0) Then
C
C  convert new geometry (internal coordinates on XC)
C  back to Z-Matrix internals
C
              CALL BackToZ(NZ,     NVib,   NVar,   IG,     XC,
     $                     XINT,   GEO)
c
              CALL PrntZPATH(IOut,   NStep,  NZ,     ZSymb,  GEO,
     $                       IGEO,   NVib,   XOrig,  XINT,   EC)
              CALL PrntZPATH(ICond,  NStep,  NZ,     ZSymb,  GEO,
     $                       IGEO,   NVib,   XOrig,  XINT,   EC)
            Else
              CALL PrntPATH(IOut,   NStep,  NAtoms, mwght,  AtSymb,
     $                      XOrig,  XC,     EC)
              CALL PrntPATH(ICond,  NStep,  NAtoms, mwght,  AtSymb,
     $                      XOrig,  XC,     EC)
            EndIf
            CALL CpyVEC(NVar,XC,XOld)
            GO TO 100
C
          Else If(parab.AND.(d1.LT.DTol.AND.d2.LT.DTol)) Then
C
C  we are very close to the minimum
C  simply take lowest energy of E1 or E2 as mimimum
C
            If(E1.LT.E2) Then
              E0 = E1
              d2 = d1
            Else
              E0 = E2
            EndIf
c
            If(IPRNT.GT.0) WRITE(IOut,1810) d2
            NStep = NStep + 1
            IStep = 1
            LStep = 0
            If(NZ.GT.0) Then
C
C  convert new geometry (internal coordinates on XC)
C  back to Z-Matrix internals
C
              CALL BackToZ(NZ,     NVib,   NVar,   IG,     XC,
     $                     XINT,   GEO)
c
              CALL PrntZPATH(IOut,   NStep,  NZ,     ZSymb,  GEO,
     $                       IGEO,   NVib,   XOrig,  XINT,   EC)
              CALL PrntZPATH(ICond,  NStep,  NZ,     ZSymb,  GEO,
     $                       IGEO,   NVib,   XOrig,  XINT,   EC)
            Else
              CALL PrntPATH(IOut,   NStep,  NAtoms, mwght,  AtSymb,
     $                      XOrig,  XC,     EC)
              CALL PrntPATH(ICond,  NStep,  NAtoms, mwght,  AtSymb,
     $                      XOrig,  XC,     EC)
            EndIf
            CALL CpyVEC(NVar,XC,XOld)
            GO TO 100
C
          Else
C
C  order the steps
C
            If(d1.GT.d2) Then
              temp = d1
              d1 = d2
              d2 = temp
              temp = E1
              E1 = E2
              E2 = temp
            EndIf
C
C  If E1 is the lowest energy then
C  fit parabola through (d0,E0), (d1,E1), (d2,E2)
C
            If(E1.LT.E0.AND.E1.LT.E2) Then
              CALL FitParab(E0,d0,E1,d1,E2,d2,.false.,dm)
              parab = .True.
              If(IPRNT.GT.0) WRITE(IOut,1820) dm
            Else
C
C  take another line search step
C
              If(E0.LT.E1) Then
                dm = Half*d1
              Else
                dm = d2 + (d2-d1)
              EndIf
              parab = .False.
            EndIf
          EndIf
C
C  put step into D
C
          If(IPRNT.GT.0) WRITE(IOut,1730) dm
          DO I=1,NVar
          D(I) = dm*GS(I)
          EndDO
c
          LStep = LStep + 1
cc
        Else If(LStep.GT.2) Then
C
C  just taken line search step
C
cc          write(6,*) ' LStep:',lstep,' dm:',dm,' EC:',ec
cc          write(6,*) ' d0:',d0,' E0:',e0
cc          write(6,*) ' d1:',d1,' E1:',e1
cc          write(6,*) ' d2:',d2,' E2:',e2
          If(parab.AND.(EC.LT.E0.AND.EC.LT.E1.AND.EC.LT.E2)) Then
C
C  terminate line search with dm as minimum
C
            If(IPRNT.GT.0) WRITE(IOut,1800)
            E0 = EC
            NStep = NStep + 1
            IStep = 1
            LStep = 0
            If(NZ.GT.0) Then
C
C  convert new geometry (internal coordinates on XC)
C  back to Z-Matrix internals
C
              CALL BackToZ(NZ,     NVib,   NVar,   IG,     XC,
     $                     XINT,   GEO)
c
              CALL PrntZPATH(IOut,   NStep,  NZ,     ZSymb,  GEO,
     $                       IGEO,   NVib,   XOrig,  XINT,   EC)
              CALL PrntZPATH(ICond,  NStep,  NZ,     ZSymb,  GEO,
     $                       IGEO,   NVib,   XOrig,  XINT,   EC)
            Else
              CALL PrntPATH(IOut,   NStep,  NAtoms, mwght,  AtSymb,
     $                      XOrig,  XC,     EC)
              CALL PrntPATH(ICond,  NStep,  NAtoms, mwght,  AtSymb,
     $                      XOrig,  XC,     EC)
            EndIf
            CALL CpyVEC(NVar,XC,XOld)
            GO TO 100
          Else
C
C  order the steps
C
            CALL StepOrder(EC,dm,E0,d0,E1,d1,E2,d2)
cc          write(6,*) ' back from <StepOrder>'
cc          write(6,*) ' d0:',d0,' E0:',e0
cc          write(6,*) ' d1:',d1,' E1:',e1
cc          write(6,*) ' d2:',d2,' E2:',e2
C
C  If E1 is the lowest energy then
C  fit parabola through (d0,E0), (d1,E1), (d2,E2)
C
            If(E1.LT.E0.AND.E1.LT.E2) Then
              CALL FitParab(E0,d0,E1,d1,E2,d2,.false.,dm)
              parab = .True.
              If(IPRNT.GT.0) WRITE(IOut,1820) dm
            Else
C
C  take another line search step
C
              If(E0.LT.E1) Then
                dm = Half*d1
              Else
                dm = d2 + (d2-d1)
              EndIf
              parab = .False.
            EndIf
          EndIf
C
C  put step into D
C
          If(IPRNT.GT.0) WRITE(IOut,1730) dm
          DO I=1,NVar
          D(I) = dm*GS(I)
          EndDO
c
          LStep = LStep + 1
        EndIf
      ENDIF
C
C
C  Take the new step
C
 200  CONTINUE
c
      IF(.NOT.Step) THEN
        If(NZ.GT.0) CALL BackToZ(NZ,     NVib,   NVar,   IG,     XINT,
     $                           Z,      GEO)
        RETURN
      ELSE
        IF(NZ.EQ.0) THEN
          If(mwght) CALL MassWT(-1,NAtoms,XMass,D)
          CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z,
     $                small,  D)
          CALL AddVEC(NVar,XOld,D,XC)
          CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z,
     $                small,  XC)
        ELSE
          CALL AddVEC(NVar,XOld,D,XINT)
C
C  convert new geometry (internal coordinates on XC)
C  back to Z-Matrix internals
C
          CALL BackToZ(NZ,     NVib,   NVar,   IG,     XINT,
     $                 Z,      GEO)
        ENDIF
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,'** REACTION PATH SEARCH IN Z-MATRIX COORDINATES **',/)
 1010 FORMAT(/,'** REACTION PATH SEARCH IN MASS-WEIGHTED',
     $           ' CARTESIAN COORDINATES **',/)
 1020 FORMAT(/,'** REACTION PATH SEARCH IN CARTESIAN COORDINATES **',/)
 1050 FORMAT(' Start of New Reaction Path Search')
 1100 FORMAT(' Computing Gradient')
 1200 FORMAT(' Taking Hessian mode as Initial Reaction Path Step')
 1210 FORMAT(/,' Mass-Weighted Hessian Matrix:')
 1220 FORMAT(/,' Hessian Matrix after Projection:')
 1230 FORMAT(/,' Hessian Matrix in Z-Matrix Variables:')
 1240 FORMAT(/,' Hessian Eigenvalues (atomic units)')
 1250 FORMAT(1X,6F12.6)
 1300 FORMAT(' Taken Hessian mode as Initial Reaction Path Step')
 1310 FORMAT(' Taking Steepest Descent as Initial Reaction Path Step')
 1320 FORMAT(' Step size is ',F9.6)
 1330 FORMAT(' Steepest Descent step below threshold - SEARCH HALTED')
 1400 FORMAT(' Computing Gradient Bisector')
 1500 FORMAT(' Initializing Line Search')
 1510 FORMAT(' gradient projector step has length ',F9.6)
 1520 FORMAT(' gradient projector step too large - reducing to ',F9.6)
 1600 FORMAT(' Computing Energy for Line Search')
 1700 FORMAT(' Line Search  Step:',I4)
 1710 FORMAT(' Fit Parabola using two energies and projected gradient',/,
     $       '  estimated minimum is ',F9.6)
 1720 FORMAT(' Estimated minimum too far - taking step of length ',F9.6)
 1730 FORMAT(' Continuing Line Search - taking step of length ',F9.6)
 1800 FORMAT(' Successful Termination of Line Search')
 1810 FORMAT(' Energies very close - minimum taken to be ',F9.6)
 1820 FORMAT(' Fit Parabola using three energies',/,
     $       '  estimated minimum is ',F9.6)
c
      END
c ...........................................................................
c
      SUBROUTINE PrntPATH(IOut,   NStep,  NAtoms, mwght,  AtSymb,
     $                    XOrig,  XC,     EC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out current step on reaction path
C
C  ARGUMENTS
C
C  IOut    -  unit to print to
C  NStep   -  current step on reaction path
C  NAtoms  -  number of atoms
C  mwght   -  logical flag for mass-weighted coordinates
C  AtSymb  -  atomic symbols
C  XOrig   -  Original Cartesian coordinates (at TS)
C  XC      -  Cartesian coordinates
C  EC      -  current energy
C
C
      DIMENSION XC(3*NAtoms),XOrig(3*NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
      Logical mwght
C
      WRITE(IOut,1000)
      If(mwght) Then
        WRITE(IOut,1001)
      Else
        WRITE(IOut,1002)
      EndIf
c
      WRITE(IOut,1100) NStep,EC
      CALL PrntCAR(IOut,0,NAtoms,AtSymb,XC)
C
C  determine distance from origin
C
      Dist = 0.0d0
      DO 10 I=1,3*NAtoms
      Dist = Dist + (XOrig(I)-XC(I))**2
 10   CONTINUE
c
      Dist = SQRT(Dist)
      WRITE(IOut,1200) Dist
c
      WRITE(IOut,1000)
C
      RETURN
c
 1000 FORMAT(1X,70('-'))
 1001 FORMAT(' REACTION PATH IN MASS-WEIGHTED CARTESIAN COORDINATES')
 1002 FORMAT(' REACTION PATH IN CARTESIAN COORDINATES')
 1100 FORMAT('  Step: ',I4,'  Energy is ',F18.9)
 1200 FORMAT('  Distance from Origin of Search: ',F9.6,'  atomic units')
c
      END
c ...........................................................................
c
      SUBROUTINE PrntZPATH(IOut,   NStep,  NZ,     ZSymb,  GEO,
     $                     IGEO,   NVib,   XOrig,  ZINT,   EC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out current step on reaction path
C
C  ARGUMENTS
C
C  IOut    -  unit to print to
C  NStep   -  current step on reaction path
C  NZ      -  number of Z-Matrix centres
C  ZSymb   -  Z-Matrix symbols
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO    -  Z-matrix connectivity
C  NVar    -  number of independent Z-matrix variables
C  XOrig   -  Original Z-Matrix coordinates (at TS)
C  ZINT    -  current Z-Matrix coordinates
C  EC      -  current energy
C
C
      DIMENSION GEO(NZ,3),IGEO(NZ,4),ZINT(NVib),XOrig(NVib)
      CHARACTER*8 ZSymb(NZ)
C
      WRITE(IOut,1000)
      WRITE(IOut,1010)
c
      WRITE(IOut,1100) NStep,EC
      CALL PrntZMAT(IOut,NZ,ZSymb,GEO,IGEO)
C
C  determine distance from origin
C
      Dist = 0.0d0
      DO 10 I=1,NVib
      Dist = Dist + (XOrig(I)-ZINT(I))**2
 10   CONTINUE
c
      Dist = SQRT(Dist)
      WRITE(IOut,1200) Dist
c
      WRITE(IOut,1000)
C
      RETURN
c
 1000 FORMAT(1X,70('-'))
 1010 FORMAT(' REACTION PATH IN Z-MATRIX COORDINATES')
 1100 FORMAT('  Step: ',I4,'  Energy is ',F18.9)
 1200 FORMAT('  Distance from Origin of Search: ',F9.6,
     $       '  Z-matrix units')
c
      END
c ...........................................................................
c
      SUBROUTINE StepOrder(EC,dm,E0,d0,E1,d1,E2,d2)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Orders line search steps according to distance
C  and selects three values for next line search step
C  choosing values d0, d1, d2 such that E0 > E1 < E2
C  wherever possible
C  NOTE: steps should already be ordered such that d0 < d1 < d2
C        so it is just a question of inserting dm
C
      IF(dm.LT.d1) THEN
        temp = d2
        d2 = d1
        d1 = dm
        dm = temp
        temp = E2
        E2 = E1
        E1 = EC
        EC = temp
      ELSE IF(dm.GT.d1.AND.dm.LT.d2) THEN
        temp = d2
        d2 = dm
        dm = temp
        temp = E2
        E2 = EC
        EC = temp
      ENDIF
C
C  distance ordering is now d0 < d1 < d2 < dm
C  select three values depending on energies
C
      IF(E2.LT.E1) THEN
        d0 = d1
        d1 = d2
        d2 = dm
        E0 = E1
        E1 = E2
        E2 = EC
      ENDIF
C
      RETURN
      END
c ..........................................................................
c
      SUBROUTINE FitParab(E0,d0,E1,d1,E2,d2,grd,dm)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Fits a parabola through three points or two points and the
C  projected gradient
C  For three-point fit points should satisfy
C   d0 < d1 < d2    and   E0 > E1 < E2
C
C  ARGUMENTS
C
C  (E0,d0)  -  first point
C  (E1,d1)  -  second point (d1 may be projected gradient)
C  (E2,d2)  -  third point
C  grd      -  logical flag for projected gradient
C
C  on exit
C
C  dm       -  estimate of minimum
C
      Logical grd
C
cc      write(6,*) ' In <FitParab>'
cc      write(6,*) ' d0:',d0,' E0:',e0
cc      write(6,*) ' d1:',d1,' E1:',e1
cc      write(6,*) ' d2:',d2,' E2:',e2
      IF(grd) Then
C
C  two points plus projected gradient
C  only done, if at all, from origin of line search
C  (so d0 = 0)
C
        B = -d1
        A = (E1-E0-B*d1)/(d1**2)
c
        dm = - B/(A+A)
C
      ELSE
C
C  three-point fit
C
        ed1 = (E2-E1)/(d2-d1)
        ed0 = (E2-E0)/(d2-d0)
c
        A =  (ed1-ed0)/(d1-d0)
        B =  ed1 - A*(d2+d1)
c
        dm = - B/(A+A)
cc        write(6,*) ' In <FitParab>  Three-Point fit'
cc        write(6,*) ' d0: ',d0,' d1:',d1,' d2:',d2
cc        write(6,*) ' estimated minimum is ',dm
c
      ENDIF
C
      RETURN
      END
