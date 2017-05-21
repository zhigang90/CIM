c ==============================================================
c  GENERAL QM/MM ONIOM MODULE          JB   Jan. 2001
c ==============================================================
c
      subroutine preqmmm(inp,scyc)
      implicit real*8(a-h,o-z)
      character*256 chopval
      integer scyc
c
c  reads the QMMM line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=2)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character*20 cdum
      character*256 jobname
      logical found
c
      parameter (IUnit=1)
      Common /job/jobname,lenJ
c
      data options/'prin','qmmm'/
      data ioptyp/1,0/
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c ...........................................................
cc      If(scyc.gt.0) RETURN        ! only read first time
c ...........................................................
c
c -- Currently the only option is the print flag
c
c -- print flag
      if(ifound(1).eq.1) then
        IPRNT = iopval(1,1)
      else
        IPRNT = 2                  ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c .............................................................................
c
      SUBROUTINE QMMM(NMem,Z,IQM)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** GENERAL QM/MM PROGRAM **
C
C  This module implements a general ONIOM QM/MM procedure.
C  Energy/gradient computation is done elsewhere using standard
C  PQS modules; all that is done here is to form the model QM
C  system (by adding "link" hydrogen atoms to terminate the QM
C  region) and energy/gradient manipulations to form the final
C  ONIOM energy and gradient, e.g., prior to the next optimization
C  step.
C
C
C  References
C  ----------
C
C  "A New ONIOM Implementation in Gaussian98. Part I. The Calculation
C   of Energies, Gradients, Vibrational Frequencies and Electric
C   Field Derivatives"
C   S.Dapprich, I.Komaromi, K.S.Byun, K.Morokuma and M.Frisch
C    J.Mol.Struct.(Theochem)  461-462 (1999) 1
C
C  ----------------------------------------------------------------------
      CHARACTER cdum*20,jobname*256
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/              !  unit number for checkpoint I/O
      Common /job/jobname,lenJ
C
      IOut = ioutfil('iout')
C
C  Read from the <control> file
C    current number of atomic centres, NO DUMMIES ALLOWED
C    current number of molecules
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      call rdcntrl(IUnit,5,'$nmol',1,NMol,rdum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Read from the <sym> file
C    total number of real atoms
C    number of symmetry operations
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
      
C
C  Get the total number of atoms in the FULL system from the depository
C
      call tstival('nqmmm',iexist)
      if(iexist.eq.1) then
       call getival('nqmmm',NAtoms)
      else
       NAtoms = NAtom   ! on first entry dealing with FULL system
       call setival('nqmmm',NAtoms)
      endif
C
C ...................................................................
C  ** IMPORTANT **
C  The number of atoms will change depending on whether we are
C  dealing with the FULL or the MODEL system. The TOTAL number
C  of atoms in the FULL system is written to the depository and
C  should be accessed from there. The number of atoms on the
C  <control> and <sym> files will vary, depending on whether we
C  are dealing with the FULL or the MODEL system.
C  Thus    NAtom  - # atoms in current system
C          NAtoms - # atoms in FULL system
C  At the start of a QM/MM "loop", NAtom = NAtoms
C ...................................................................
C
C  Allocate the maximum scratch storage
C
      NScr = 10                      ! no scratch needed so far
C
C  Now get the memory
C
      IMem = 22*NAtoms + 9*NTrans + NAtoms*NTrans
     $        + 18*NAtoms**2 + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,5,'QM/MM')
C
C  Allocate memory pointers
C
      IXC = iptr                     !  current coordinates
      IMOL = IXC + 3*NAtoms          !  pointer to start/end of molecules
      IAN = IMOL + NAtoms            !  atomic numbers
      IXCG = IAN + NAtoms            !  current atomic charges
      IXCM = IXCG + NAtoms           !  current atomic masses
      ITN = IXCM + NAtoms            !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*NTrans           !  list of atomic equivalences
      IUQ  = INQ + NAtoms*NTrans     !  list of symmetry-unique atoms
      IXF =  IUQ + NAtoms            !  full coordinates
      ICGF = IXF + 3*NAtoms          !  full atomic charges
      ICMF = ICGF + NAtoms           !  full atomic masses
      IGC = ICMF + NAtoms            !  current gradient
      IGF = IGC + 3*NAtoms           !  full gradient (QM/MM)
      IHS = IGF + 3*NAtoms           !  current Hessian matrix
      IHF = IHS + 9*NAtoms**2        !  full Hessian  (QM/MM)
      ICT = IHF + 9*NAtoms**2        !  storage for Link atom connectivity
      IG  = ICT + 2*NAtoms           !  bond factors for link atoms
      IEnd = IG + NAtoms
c
      IScr = IEnd + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(NMem,IEnd,5,'QM/MM')
C
C
c memory status
c -- assign memory for high water mark (old TEXAS)
      call getmem(IEnd,lastx)
C
C  ----------------------------------------------------------------------
C
      CALL QMMMMAIN(NAtom,   NAtoms,  NMol,    Z(IMOL), Z(IAN),
     $              Z(IXC),  Z(IXCG), Z(IXCM), NTrans,  Z(ITN),
     $              Z(INQ),  Z(IUQ),  Z(IXF),  Z(ICGF), Z(ICMF),
     $              Z(IGC),  Z(IGF),  Z(IHS),  Z(IHF),  Z(ICT),
     $              Z(IG),   NScr,    Z(IScr), IQM)
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
      write(IOut,1100) IEnd,mxmem,memtot
 1100 format(/,' QM/MM memory status:',/,
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
      SUBROUTINE QMMMMAIN(NAtom,  NAtoms, NMol,   IMOL,   IAN,
     $                    XC,     XCharg, XMass,  NTrans, TRANS,
     $                    NEqATM, IUNQ,   XCF,    XChrgF, XMassF,
     $                    GC,     GF,     HESS,   HF,     ICnnct,
     $                    g,      NMem,   Z,      IQM)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for general QM/MM module
C  QMMMMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION XC(3*NAtoms),XCharg(NAtoms),XMass(NAtoms),
     $          XCF(3,NAtoms),XChrgF(NAtoms),XMassF(NAtoms),
     $          GC(3*NAtoms),GF(3*NAtoms),
     $          HESS(9*NAtoms**2),HF(9*NAtoms**2)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans),
     $          IUNQ(NAtoms),IMOL(NAtoms),IAN(NAtoms),
     $          ICnnct(2*NAtoms),g(NAtoms),RM(3,3)
c ............................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms),AtSymF(NAtoms)
c ............................................................
      CHARACTER GROUP*4,cdum*20,jobname*256
      LOGICAL Symflag
C
      DIMENSION Z(NMem)
C
      PARAMETER (thrbnd=0.85d0)  !  distance ratio for bonding
C
      Common /job/jobname,lenJ
      Data IUnit/1/
C
C
c     common/big/bl(100)         !  old texas depository
      IOut = ioutfil('iout')
      ICond = ioutfil('icond')
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
      CALL RdSYM(Symflag,NAtom,  RM,     GROUP,  NTrans,
     $           NDEG,   NQ,     IUNQ,   TRANS,  NEqATM)
C
C  Read from the <control> file
C    current energy
C    print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$energy',2,idum,EC,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  read initial coordinates from <coord> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoordF(IUnit,  NAtom,  AtSymb, XC,     NMol,
     $              IMOL,   XCharg, XMass)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Attempt to read <qmmm> file
C  If no <qmmm> file is available, this is the first
C  step of a new QM/MM job
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.qmmm',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(40) IQM
      READ(40) NAtoms,MAtom,MAtoms,ICnnct,g,IAN,XCF,
     $         AtSymF,XChrgF,XMassF,EQmm,GF,HF
      CLOSE (UNIT=40,STATUS='KEEP')
      GO TO 20
C
C  Nothing on file
C  Initialize
C
 10   CONTINUE
      IQM = 0
C
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
C  save original atomic symbols,charges & masses
C
      DO 15 IAtm=1,NAtoms
      AtSymF(IAtm) = AtSymb(IAtm)
      XChrgF(IAtm) = XCharg(IAtm)
      XMassF(IAtm) = XMass(IAtm)
 15   CONTINUE
C
C  determine connectivity between model and full system
C
      CALL QMConnect(NAtoms, NMol,   IMOL,   IAN,    XC,
     $               thrbnd, ICnnct, g,      MAtom,  MAtoms)
C
 20   CONTINUE
C
C  NOTE: Throughout this module
C    NAtoms is the total number of atoms in the FULL system
C    MAtom  is the number of atoms in the QM region
C    MAtoms is the number of atoms in the MODEL system
C           (MAtom + link atoms)
C    NAtom  is the number of atoms in the current system
C           (could be NAtoms OR MAtoms)
C
C  read gradient from <grad> file
C
      CALL RdGRAD(NAtom,GC,'save')
      CALL VScal(3*NAtom,-1.0d0,GC)            ! file contains forces
C
C  Is there a Hessian file to read?
C    <Hessian stuff here>
C
C  Determine what QM/MM action to take
C
C ---------------------------------------------------------------------
      CALL QMMMF(NAtoms, MAtom,  MAtoms, IAN,    AtSymb,
     $           XCharg, XMass,  XC,     EC,     EQmm,
     $           IPRNT,  XCF,    GC,     GF,     HESS,
     $           HF,     ICnnct, IQM,    g)
C ---------------------------------------------------------------------
C
C  Write <qmmm> file for next step
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.qmmm',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(40) IQM
      WRITE(40) NAtoms,MAtom,MAtoms,ICnnct,g,IAN,XCF,
     $          AtSymF,XChrgF,XMassF,EQmm,GF,HF
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  Write new <coord> file and write current size of system
C  to control file
C    If IQM=0   dealing with FULL system
C       IQM=1   dealing with MODEL system, MM part
C       IQM=2   dealing with MODEL system, QM part
C
      If(IQM.EQ.0) Then
       NAtom = NAtoms
       NMol = 2
       IMOL(2) = MAtom
       IMOL(3) = NAtoms
C
C  restore original atomic symbols, charges & masses
C
       DO 25 IAtm=1,NAtoms
       AtSymb(IAtm) = AtSymF(IAtm)
       XCharg(IAtm) = XChrgF(IAtm)
       XMass(IAtm) = XMassF(IAtm)
 25    CONTINUE
C
C  write modified gradient to <grad> file
C
       CALL VScal(3*NAtoms,-1.0d0,GF)            ! file contains forces
       CALL WrGRAD(NAtoms,GF)
      Else
       NAtom = MAtoms
       NMol = 1
      EndIf
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      If(IQM.EQ.0) Call wrcntrl(IUnit,7,'$energy',2,idum,EQmm,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
cc      write(6,*) ' About to call <WrCoord>  NAtom:',natom
cc      do i=1,natom
cc      write(6,*) i,'  ',atsymb(i),XCharg(i),xmass(i)
cc      enddo
      CALL WrCoord(NAtom,  AtSymb, XC,     NMol,   IMOL,
     $             XCharg, XMass)
C
C  save current number of atoms and current atomic masses
C  to depository
C
      call setival('na',NAtom)
      call getival('mass',iadr)
      do i=1,natom
      bl(iadr+i-1) = XMass(i)
      enddo
C
      RETURN
      END
c ..............................................................................
c
      SUBROUTINE QMMMF(NAtoms, MAtom,  MAtoms, IAN,    AtSymb,
     $                 XCharg, XMass,  XC,     EC,     EQmm,
     $                 IPRNT,  XCF,    GC,     GF,     HESS,
     $                 HF,     ICnnct, IQM,    g)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main QM/MM driver
C
C  ARGUMENTS
C
C  NAtoms   -  number of atoms in FULL system
C  MAtom    -  number of atoms in QM region
C  MAtoms   -  number of atoms in MODEL system
C  IAN      -  atomic numbers
C  AtSymb   -  atomic symbols for current system (used for nice printout)
C  XCharg   -  atomic charges
C  XMass    -  atomic masses
C  XC       -  Cartesian coordinates of current system
C  EC       -  current energy
C  EQmm     -  QM/MM energy
C  IPRNT    -  print flag
C  XCF      -  Cartesian coordinates of FULL system
C  GC       -  gradient of current system
C  GF       -  QM/MM gradient
C  HESS     -  Hessian matrix for current system
C  HF       -  QM/MM Hessian
C  ICnnct   -  connectivity data for link atoms
C  IQM      -  QM/MM step (on entry)
C               0 - dealing with FULL system
C               1 - dealing with MODEL system, MM part
C               2 - dealing with MODEL system, QM part
C              will be reset on exit ready for next entry
C  g        -  bond factors for link atoms
C
C
      DIMENSION XCharg(NAtoms),XMass(NAtoms),XC(3,NAtoms),
     $          XCF(3*NAtoms),GC(3,NAtoms),GF(3,NAtoms),
     $          HESS(9*NAtoms**2),HF(3*NAtoms,3*NAtoms)
      DIMENSION ICnnct(2*NAtoms),g(NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
C
C
C .................................................................
C  What needs to be done in the routine is controlled by IQM
C
C  IQM=0:  Start of a new QM/MM calculation
C          Should have MM energy and gradient for FULL system
C          Save details of FULL system and prepare MODEL system
C          On exit, dealing with MODEL system - increment IQM
C  IQM=1:  Starting second part of QM/MM calculation
C          Should have MM energy and gradient for MODEL system
C          Modify QM/MM energy and gradient
C          On exit, dealing with MODEL system - increment IQM
C  IQM=2:  Starting final part of QM/MM calculation
C          Should have QM energy and gradient for MODEL system
C          Modify QM/MM energy and gradient
C          Restore FULL system
C          On exit, dealing with FULL system - reset IQM=0
C .................................................................
C
      IOut = ioutfil('iout')
      ICond = ioutfil('icond')
c
      If(IPRNT.GT.0) WRITE(IOut,1000)
c
      NLink = MAtoms - MAtom       ! number of link atoms
c
      IF(IQM.EQ.0) THEN
cc
        If(IPRNT.GT.0) WRITE(IOut,1100)
C
C  save current geometry, energy and gradient
C
        CALL CpyVEC(3*NAtoms,XC,XCF)
        CALL CpyVEC(3*NAtoms,GC,GF)
        EQmm = EC
c
        If(IPRNT.GT.3)
     $     write(IOut,*) ' MM energy:',ec,'  Total energy:',EQmm
        If(IPRNT.GT.4) Then
          write(IOut,*) ' MM Gradient (au) for full system is:'
          do i=1,natoms
          write(IOut,1234) i,gc(1,i),gc(2,i),gc(3,i)
          enddo
        EndIf
C
C  prepare MODEL system
C
        CALL GetMODEL(NAtoms, MAtom,  MAtoms, IAN,    XCF,
     $                ICnnct, g,      IPRNT,  XC,     AtSymb,
     $                XCharg, XMass)
c
        IQM = IQM+1
cc
      ELSE IF(IQM.EQ.1) THEN
cc
        If(IPRNT.GT.0) WRITE(IOut,1200)
C
C  should have MM energy and gradient of model system
C  incorporate into FULL system
C
        EQmm = EQmm - EC
c
        If(IPRNT.GT.3)
     $     write(IOut,*) ' MM energy:',ec,'  Total energy:',EQmm
        If(IPRNT.GT.4) Then
          write(IOut,*) ' MM Gradient (au) for model system is:'
          do i=1,matoms
          write(IOut,1234) i,gc(1,i),gc(2,i),gc(3,i)
          enddo
        EndIf
C
C  modify gradient
C
        DO 10 I=1,MAtom
        GF(1,I) = GF(1,I) - GC(1,I)
        GF(2,I) = GF(2,I) - GC(2,I)
        GF(3,I) = GF(3,I) - GC(3,I)
 10     CONTINUE
C
C  now project link atoms
C
        DO 11 I=1,NLink
        I1 = ICnnct(2*I-1)     ! atom in QM region
        I3 = ICnnct(2*I)       ! atom it is connect to in MM region
        II = MAtom+I           ! link H atom between these 2 atoms
c
        GF(1,I1) = GF(1,I1) - (1-g(I))*GC(1,II)
        GF(2,I1) = GF(2,I1) - (1-g(I))*GC(2,II)
        GF(3,I1) = GF(3,I1) - (1-g(I))*GC(3,II)
c
        GF(1,I3) = GF(1,I3) - g(I)*GC(1,II)
        GF(2,I3) = GF(2,I3) - g(I)*GC(2,II)
        GF(3,I3) = GF(3,I3) - g(I)*GC(3,II)
 11     CONTINUE
c
        IQM = IQM+1
cc
      ELSE IF(IQM.EQ.2) THEN
cc
        If(IPRNT.GT.0) WRITE(IOut,1300)
C
C  should have QM energy and gradient of model system
C  incorporate into FULL system
C
        EQmm = EQmm + EC
c
        If(IPRNT.GT.3)
     $     write(IOut,*) ' QM energy:',ec,'  Total energy:',EQmm
        If(IPRNT.GT.4) Then
          write(IOut,*) ' QM Gradient (au) for model system is:'
          do i=1,matoms
          write(IOut,1234) i,gc(1,i),gc(2,i),gc(3,i)
          enddo
        EndIf
C
C  modify gradient
C
        DO 20 I=1,MAtom
        GF(1,I) = GF(1,I) + GC(1,I)
        GF(2,I) = GF(2,I) + GC(2,I)
        GF(3,I) = GF(3,I) + GC(3,I)
 20     CONTINUE
C
C  now project link atoms
C
        DO 21 I=1,NLink
        I1 = ICnnct(2*I-1)     ! atom in QM region
        I3 = ICnnct(2*I)       ! atom it is connect to in MM region
        II = MAtom+I           ! link H atom between these 2 atoms
c
        GF(1,I1) = GF(1,I1) + (1-g(I))*GC(1,II)
        GF(2,I1) = GF(2,I1) + (1-g(I))*GC(2,II)
        GF(3,I1) = GF(3,I1) + (1-g(I))*GC(3,II)
c
        GF(1,I3) = GF(1,I3) + g(I)*GC(1,II)
        GF(2,I3) = GF(2,I3) + g(I)*GC(2,II)
        GF(3,I3) = GF(3,I3) + g(I)*GC(3,II)
 21     CONTINUE
cc
        If(IPRNT.GT.4) Then
          write(IOut,*) ' Final QM/MM gradient is:'
          do i=1,natoms
          write(IOut,1234) i,gf(1,i),gf(2,i),gf(3,i)
          enddo
        EndIf
C
C  restore geometry of FULL system
C
        If(IPRNT.GT.0) WRITE(IOut,1350)
        CALL CpyVEC(3*NAtoms,XCF,XC)
c
        IQM = 0
cc
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,'** QM/MM MODULE **',/)
 1100 FORMAT(' MM Calculation on FULL system',/,
     $       ' Generating Link Atoms and Preparing Model System')
 1200 FORMAT(' MM Calculation on MODEL system - modifying gradient')
 1300 FORMAT(' QM Calculation on MODEL system - modifying gradient')
 1350 FORMAT(' Restoring FULL system')
 1234 format(1X,i4,2X,3(F15.7))
c
      END
c ..............................................................................
c
      SUBROUTINE QMConnect(NAtoms, NMol,   IMOL,   IAN,    XC,
     $                     thrbnd, ICnnct, g,      MAtom,  MAtoms)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines which atoms in the inner QM layer are connected
C  to the outer MM layer for a QM/MM calculation
C  Based on simple interatomic distance criterion
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms in FULL system
C  NMol    -  number of "molecules" in system
C             (should be 2 in this implementation)
C  IMOL    -  pointers to start/end of molecules in XC array
C  IAN     -  atomic numbers
C  XC      -  current Cartesian coordinates
C  thrbnd  -  distance ratio above which two atoms are
C             considered to be connected
C
C  on exit
C
C  ICnnct  -  list of atoms in QM region that are connected
C             to the outer MM region stored as pairs of atoms
C               atom in inner region   atom in outer region
C             if there is more than one connection to a given
C             atom, it will have multiple entries
C  g       -  bond factors for link atoms
C  MAtom   -  number of atoms in QM region
C  MAtoms  -  total number of atoms in MODEL system
C             (number of atoms in QM region + number of connections)
C
C
      DIMENSION IMOL(NMol+1),IAN(NAtoms),XC(3,NAtoms)
      DIMENSION ICnnct(2*NAtoms),g(NAtoms)
C
      REAL*8 BndLen(92)      ! homonuclear bond distances in au
C
c                    H-H  He-He Li-Li Be-Be  B-B   C-C   N-N   O-O
      DATA BndLen / 1.42, 4.00, 5.10, 4.16, 3.21, 2.93, 2.74, 2.74,
c
c                    F-F  Ne-Ne Na-Na Mg-Mg Al-Al Si-Si  P-P   S-S
     $              2.65, 6.00, 5.86, 5.67, 5.01, 4.35, 4.16, 3.87,
c
c                   Cl-Cl Ar-Ar  K-K  Ca-Ca Sc-Sc Ti-Ti  V-V  Cr-Cr
     $              3.78, 8.00, 7.37, 6.58, 5.44, 5.00, 4.62, 4.42,
c
c                   Mn-Mn Fe-Fe Co-Co Ni-Ni Cu-Cu Zn-Zn Ga-Ga Ge-Ge
     $              4.40, 4.38, 4.38, 4.36, 4.35, 4.72, 4.72, 4.62,
c
c                   As-As Se-Se Br-Br Kr-Kr Rb-Rb Sr-Sr  Y-Y  Zr-Zr
     $              4.58, 4.44, 4.35, 9.00, 8.16, 7.22, 6.12, 5.50,
c
c                   Nb-Nb Mo-Mo Tc-Tc Ru-Ru Rh-Rh Pd-Pd Ag-Ag Cd-Cd
     $              5.06, 4.88, 4.80, 4.70, 4.72, 4.84, 5.06, 5.33,
c
c                   In-In Sn-Sn Sb-Sb Te-Te  I-I  Xe-Xe Cs-Cs Ba-Ba
     $              5.67, 5.30, 5.33, 5.18, 5.03, 9.00, 8.88, 7.48,
c
c                   La-La Ce-Ce Pr-Pr Nd-Nd Pm-Pm Sm-Sm Eu-Eu Gd-Gd
     $              6.40, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00,
c
c                   Tb-Tb Dy-Dy Ho-Ho Er-Er Tm-Tm Yb-Yb Lu-Lu Hf-Hf
     $              6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 5.44,
c
c                   Ta-Ta  W-W  Re-Re Os-Os Ir-Ir Pt-Pt Au-Au Hg-Hg
     $              5.06, 4.92, 4.84, 4.76, 4.76, 4.88, 5.06, 5.44,
c
c                   Tl-Tl Pb-Pb Bi-Bi Po-Po At-At Rn-Rn Fr-Fr Ra-Ra
     $              5.86, 5.82, 5.75, 6.00, 6.00, 6.00, 6.00, 6.00,
c
c                   Ac-Ac Th-Th Pa-Pa  U-U
     $              6.00, 6.00, 6.00, 6.00 /
C
C
C  loop over atoms in QM layer
C
      IT = 0
      MAtom = IMOL(2)
c
      DO 10 I=1,MAtom
      DO 10 J=MAtom+1,NAtoms
C
C  get interatomic distance
C
      Dist = SQRT( (XC(1,I) - XC(1,J))**2 +
     $             (XC(2,I) - XC(2,J))**2 +
     $             (XC(3,I) - XC(3,J))**2 )
C
C  get "standard" I-J bond length from BndLen table
C
      bondL = SQRT(BndLen(IAN(I))*BndLen(IAN(J)))
C
C  If the "standard" bond length divided by the actual interatomic
C  distance is greater than thrbnd, the atoms are considered to
C  be bonded
C
      IF(bondL/Dist.GT.thrbnd) THEN
       If(IT+2.GT.2*NAtoms) GO TO 95    ! error check
       IT = IT+1
       ICnnct(IT) = I
       IT = IT+1
       ICnnct(IT) = J
C
C  get bond factor
C
       bondH = SQRT(BndLen(IAN(I))*BndLen(1))
       g(IT/2) = bondH/bondL
      ENDIF
c
 10   CONTINUE
C
C  include the link atoms
C
      NLink = IT/2
      MAtoms = MAtom + NLink
c
c -- save data in depository
      call setival('nlink',NLink)
C
      RETURN
C
 95   CONTINUE
c -- should never get here!
      Call nerror(1,'QM/MM module',
     $              'Too many bonds between QM and MM layers!',0,0)
C
      END
c ..............................................................................
c
      SUBROUTINE GetMODEL(NAtoms, MAtom,  MAtoms, IAN,    XCF,
     $                    ICnnct, g,      IPRNT,  XC,     AtSymb,
     $                    XCharg, XMass)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Form the MODEL system from the QM-MM connectivity and bond factors
C  This is done by appending link atoms to the list of atoms in the QM region
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms in FULL system
C  MAtom   -  number of atoms in QM region
C  MAtoms  -  number of atoms in MODEL system
C  IAN     -  atomic numbers
C  XCF     -  Cartesian coordinates for FULL system
C  ICnnct  -  list of atoms in QM region that are connected
C             to the outer MM region stored as pairs of atoms
C               atom in inner region   atom in outer region
C             if there is more than one connection to a given
C             atom, it will have multiple entries
C  g       -  bond factors for link atoms
C  IPRNT   -  print flag
C  XC      -  on input Cartesian coordinates for FULL system
C             on exit  Cartesian coordinates for MODEL system
C  AtSymb  -  on input atomic symbols for FULL system
C             on exit  atomic symbols for MODEL system
C  XCharg  -    ditto atomic charges
C  XMass   -    sitto atomic masses
C
C
      DIMENSION IAN(NAtoms),ICnnct(2*NAtoms),g(NAtoms)
      DIMENSION XCF(3,NAtoms),XC(3,NAtoms),XCharg(NAtoms),XMass(NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
C
      PARAMETER (One=1.0d0,HMass=1.00794d0)
C
C
      NLink = MAtoms - MAtom       ! number of link atoms
cc      write(6,*) ' MAtoms:',matoms,' MAtom:',matom,' NLink:',nlink
c
      If(IPRNT.GT.2) Then
        WRITE(6,1000)
        DO I=1,NLink
        WRITE(6,1100) ICnnct(2*I-1),ICnnct(2*I),g(I)
        EndDO
      EndIf
C
C  loop over each connected pair
C
      DO 10 I=1,NLink
      I1 = ICnnct(2*I-1)
      I3 = ICnnct(2*I)
      II = MAtom+I
C
C  put coordinates of link atom into XC array
C  at the same time add default charge and mass for H atom
C
      XC(1,II) = XCF(1,I1) + g(I)*(XCF(1,I3) - XCF(1,I1))
      XC(2,II) = XCF(2,I1) + g(I)*(XCF(2,I3) - XCF(2,I1))
      XC(3,II) = XCF(3,I1) + g(I)*(XCF(3,I3) - XCF(3,I1))
c
      XCharg(II) = One
      XMass(II)  = HMass
      AtSymb(II) = 'h       '
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(' Bond Connectivity Data Between QM and MM Regions:')
 1100 FORMAT(' Atom ',I4,' QM Connected to Atom ',I4,' MM',
     $       '   bond factor ',F10.6)
c
      END
