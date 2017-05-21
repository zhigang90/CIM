c ==============================================================
c  POTENTIAL SCAN MODULE                   JB   July 2000
c      updated for external force scan     JB    Feb 2012
c ==============================================================
c
      subroutine prescan(inp,nline)
      implicit real*8(a-h,o-z)
C
C  ARGUMENTS
C
C  inp     -  unit number for input file
C  nline   -  number of lines following SCAN line
C             (e.g., for definition of composite coordinate)
C             such lines are skipped in all entries to this routine
C             EXCEPT the first
C             on first entry, nline should be 0
C             on subsequent entries, number of lines to skip
C
      character*80 Char
      character*256 jobname
      character*4 Word,Coord
      Logical zmat,comp,force
c
c  reads the SCAN line in the input file and writes options
c  (must be some) to the <control> file
c
c  Input reading is handled differently in SCAN
c  input card must be EITHER of the form
c         SCAN  coord  <atom list>  FROM <range> <step>
c  e.g.   SCAN  tors   I  J  K  L   FROM -30.0 30.0 5.0
c  OR of the form (for a composite scan)
c         SCAN  comp  FROM <range> <step>
c         #composite
c         <definition of composite coordinate>
c         #endcomposite
c  e.g.   SCAN COMP  FROM  5.0 4.0 -0.1
c         #composite
c         stre  5  11
c         stre  6  13
c         #endcomposite
c  OR of the form (for a Z-matrix scan)
c         SCAN ZMAT <Z-matrix variable>  FROM <range> <step>
c  e.g.   SCAN ZMAT  L1  FROM  1.0 1.6 0.05
c  OR of the form (for an external force scan)
c         SCAN  EXTF  <2 atoms>  FROM <range> <step>
c  e.g.   SCAN  extf   9  16  FROM  0.3  0.6  0.02
c         If applied force is +ve, assumed PULL
c         If applied force is -ve, assumed PUSH
c
c  Primitive internals that can be scanned include:
c     stre  -  distance between 2 atoms
c     bend  -  bond angle between 3 atoms
c     outp  -  out-of-plane bend involving 4 atoms  (L is central atom)
c     tors  -  proper torsion involving 4 atoms
c     comp  -  composite coordinate (any linear combination of the above)
c
      parameter (IUnit=1,MaxP=20,One=1.0d0)
      DIMENSION IPCOMP(5,MaxP),PCWght(MaxP)
      Common /job/jobname,lenJ
      Data zmat/.False./, comp/.False./, NPComp/0/, force/.False./
c
c
c -- read SCAN line
      Read(inp,900) Char
c ...........................................................
      If(nline.EQ.-1) RETURN        ! only read first time
      IF(nline.gt.0) THEN
        DO I=1,nline
        Read(inp,900) Char
        EndDO
        RETURN
      ENDIF
c ...........................................................
c
c -- convert to lower case
      call lowerca2(Char,80)
c
c -- remove any leading blanks
      call leadblan2(Char,80,len)
c
      call setival('isumscf',1)     ! switch off SCF summary print
c
c -- initialize
      nline = -1
      n = 1
c
c -- start to interpret string
c -- read first keyword - this must be SCAN
c
      call rdword(Char,n,80,Word,4,nk,ival)
c
c -- read second keyword - this is the internal coordinate type to scan
c -- or ZMAT, indicating a Z-Matrix scan, or EXTF
c
      call rdword(Char,n,80,Coord,4,nk,ival)
c
      IF(coord.EQ.'zmat') THEN
c
c -- Z-Matrix scan
c -- read Z-Matrix variable to scan
c
       call rdword(Char,n,80,Coord,4,nk,ival)
       zmat = .True.
c
      ELSE IF(coord.EQ.'extf') THEN
c -- we have an external force scan
        force = .True.
        GO TO 90
c
      ELSE IF(coord.EQ.'comp') THEN
c
c -- scan of composite coordinate
       comp = .True.
c
      ELSE
c
c -- scan using delocalized internals
c -- coord should be one of:
c      stre   stretch
c      bend   planar bend
c      outp   out-of-plane bend
c      tors   torsion
c
       K = 0
       L = 0
c
       If(Coord.eq.'stre') Then
        Read(Char(n:80),*,END=95,ERR=95) I,J
c -- read two "words" (corresponding to integer atom positions)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
       Else If(Coord.eq.'bend') Then
        Read(Char(n:80),*,END=95,ERR=95) I,J,K
c -- read three "words" (corresponding to integer atom positions)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
       Else
        Read(Char(n:80),*,END=95,ERR=95) I,J,K,L
c -- read four "words" (corresponding to integer atom positions)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
       EndIf
       zmat = .False.
      ENDIF
c
c -- read third keyword - should be FROM
      call rdword(Char,n,80,Word,4,nk,ival)
c
c -- now read range and step length
      Read(Char(n:80),*,END=95,ERR=95) R1,R2,Step
c
c -- if composite scan, read coordinate definition
      IF(comp) THEN
        READ(inp,900) Char
        If(Char(1:4).NE.'#com') Then
         Call nerror(9,'SCAN Module',
     $    'Composite Coordinate expected but None found.',
     $    ' Check Your Input',0,0)
        EndIf
c
        nt = 0
 25     CONTINUE
        READ(inp,900) Char
c
c -- how many "parameters" are on this card?
        CALL NumFIELD(Char,80,NParam)
        Word = CHAR(1:4)
c
        If(Word.EQ.'#end') Then
          NPComp = nt
          GO TO 30
cc
        Else If(Word.EQ.'stre') Then
cc
          nt = nt+1
          IPComp(1,nt) = 1
          If(NParam.EQ.4) Then
           READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),PCWght(nt)
          Else If(NParam.EQ.3) Then
           READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt)
           PCWght(nt) = One
          Else
           Call nerror(11,'SCAN Module',
     $      'Problems with input for Distance Composite Coordinate',0,0)
          EndIf
cc
        Else If(Word.EQ.'bend') Then
cc
          nt = nt+1
          IPComp(1,nt) = 2
          If(NParam.EQ.5) Then
            READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                         PCWght(nt)
          Else If(NParam.EQ.4) Then
            READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt)
            PCWght(nt) = One
          Else
            Call nerror(12,'SCAN Module',
     #       'Problems with input for Angle Composite Coordinate',0,0)
          EndIf
cc
        Else If(Word.EQ.'outp') Then
cc
          nt = nt+1
          IPComp(1,nt) = 3
          If(NParam.EQ.6) Then
            READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                         IPComp(5,nt),PCWght(nt)
          Else If(NParam.EQ.5) Then
            READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                         IPComp(5,nt)
            PCWght(nt) = One
          Else
            Call nerror(13,'SCAN Module',
     $       'Problems with input for Out-of-Plane Bend Composite',
     $       ' Coordinate',0,0)
          EndIf
cc
        Else If(Word.EQ.'tors') Then
cc
          nt = nt+1
          IPComp(1,nt) = 4
          If(NParam.EQ.6) Then
            READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                         IPComp(5,nt),PCWght(nt)
          Else If(NParam.EQ.5) Then
            READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                         IPComp(5,nt)
            PCWght(nt) = One
          Else
            Call nerror(14,'SCAN Module',
     $       'Problems with input for Torsion Composite Coordinate',0,0)
          EndIf
cc
        Else
cc
          CHAR = 'Unknown or Unsupported Composite Coordinate Type: '
     $                //Word
          Call nerror(15,'SCAN Module',CHAR,0,0)
        EndIf
c
        GO TO 25
 30     CONTINUE
c
c  At this point we should have
c   NPComp - number of primitives in composite coordinate
c   IPCOMP - type of each primitive
c   PCWght - weight of each primitive
c
        nline = NPComp+2
      ENDIF
c
 90   CONTINUE
      IF(force) THEN
        Read(Char(n:80),*,END=95,ERR=95) I,J
c -- read two "words" (corresponding to two atoms)
        call rdword(Char,n,80,Word,4,nk,ival)
        call rdword(Char,n,80,Word,4,nk,ival)
c -- read third keyword - should be FROM
        call rdword(Char,n,80,Word,4,nk,ival)
c -- now read range and step length
        Read(Char(n:80),*,END=95,ERR=95) R1,R2,Step
      ENDIF
c
c -- write SCAN data to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call WrSCAN(IUnit,zmat,Coord,I,J,K,L,NPComp,
     $            IPCOMP,PCWght,R1,R2,Step,force)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
c
 95   CONTINUE
      Call nerror(1,'SCAN Module',
     $  'Problems Reading Data for Potential Scan - check input file',
     $   0,0)
c
  900 format(A80)
c
      end
c .............................................................................
c
      SUBROUTINE PSCAN(inp,NMem,Z,Cnvgd,natoms)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** POTENTIAL SCAN PROGRAM **
C
C  This Program carries out a potential scan by varying a given internal
C  coordinate between two values, R1 and R2, with a steplength, Step.
C
C  The scan is handled in delocalized internal coordinates by isolating
C  the internal coordinate to be scanned as if it were a constraint, and
C  then constraining all the other coordinates, moving only the isolated
C  primitive.
C
C  The initial geometry can be given in Cartesian coordinates and the
C  scanned coordinate does not need to have value R1 in the starting
C  geometry; it will be set to R1 on first entry to this module.
C
C  The program can also carry out a scan over an applied external force.
C
C
      Parameter (MaxP=20)
      DIMENSION IPCOMP(5,MaxP),PCWght(MaxP)
      CHARACTER cdum*20,jobname*256
      LOGICAL zmat,Tors,force,Cnvgd
C
      DIMENSION Z(NMem)
c ............................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(natoms)
      DIMENSION IAN(natoms),IC(natoms,natoms)
c ............................................................
      Data NIC/0/, thrbnd/0.85d0/
c
      Data IUnit/1/              !  unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
      IOut = ioutfil('iout')
      NZ = 0
      Tors=.TRUE.
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
C    data for the scan
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      call rdcntrl(IUnit,5,'$nmol',1,NMol,rdum,cdum)
      call rdscan0(IUnit,zmat,NPComp,force,Cnvgd)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(.NOT.Cnvgd) call nerror(28,'SCAN Module',
     $                  'Problems reading Scan data',0,0)
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
C  For a Z-matrix scan, determine the number of centers
C  in the Z matrix
C
      IF(zmat) THEN
       CALL ScanZMAT(NZ,IErr)
       If(IErr.EQ.-1) Call nerror(2,'SCAN module',
     $    'Z-Matrix Scan Requested but there is NO Z Matrix!!',0,0)
       NIC = 0
      ELSE IF(force) THEN
       NIC = 0
      ELSE
C
C  ==============================================================
C  WARNING - Cannot do scan in delocalized internals for diatomic
C            Print an error message and exit
C
        If(NAtoms.EQ.2) Then
          Call nerror(3,'SCAN module',
     $    'Sorry! Scans on Diatomic Molecules MUST Use a Z Matrix',0,0)
        EndIf
C  ==============================================================
C
C  Determine number of primitives for a scan using
C  delocalized internals
C
       Call tstival('npscan',iexist)
       If(iexist.EQ.1) Then
         call getival('npscan',NIC)
       Else
c -- (over) estimate the maximum number of primitives
c --  first entry only
         OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.coord',
     $        FORM='FORMATTED',STATUS='OLD')
         CALL RdCoord(40,NAtoms,AtSymb,Z,-1,idum)
         CLOSE (UNIT=40,STATUS='KEEP')
         CALL GetAtNo(NAtoms,AtSymb,IAN)
         CALL IZeroIT(IC,NAtoms*NAtoms)
         CALL CONNECTM(NAtoms,IAN,Z,thrbnd,IC,IErr)
         If(IErr.NE.0) STOP '**ERROR FORMING CONNECTIVITY MATRIX**'
         CALL GetNPrim(NAtoms,Tors,IC,NIC)
       EndIf
      ENDIF
C
      NAT3 = 3*NAtoms
      NScan = 1
C
C  Allocate the maximum scratch storage
C
      NScr = 2*NAT3*NIC + 2*NIC*NIC + 3*NIC + 4*NAT3
C
C  Now get the memory
C
      IMem = 6*NAtom + NMol+1 + 6*NScan + 9*NTrans + NAtoms*NTrans +
     $       NAtoms**2 + 10*NZ + 8*NIC + MAX(NIC,3*NZ) + 3*NAtoms*NIC
     $          + 7*NPComp + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,4,'SCAN')
C
C  Allocate memory pointers
C
      IXC = iptr                     !  coordinates
      IMOL = IXC + 3*NAtom           !  pointer to start/end of molecules
      IXCG = IMOL + NMol+1           !  atomic charges
      IXCM = IXCG + NAtom            !  atomic masses
      ICT = IXCM + NAtom             !  primitive type to scan
      ICN = ICT + NScan              !  list of atoms involved in scan primitive
      ICNM = ICN + 4*NScan           !  positive of scan primitive
      ITN = ICNM + NScan             !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*NTrans           !  list of atomic equivalences
      IUQ  = INQ + NAtoms*NTrans     !  list of symmetry-unique atoms
      iccn = IUQ + NAtoms            !  atomic connectivity matrix
      ityp = iccn + NAtoms**2        !  type of primitive internal
      ilst = ityp + NIC              !  list of atoms in primitives
      icof = ilst + 4*NIC            !  primitive weights
      isav = icof + NIC              !  values of primitive torsions
      iut = isav + NIC               !  transformation vectors
      ixpm = iut + NIC*NAT3          !  values of primitives
      IZO = ixpm + NIC               !  Z-matrix parameter values
      IZS = IZO + 3*NZ               !  Z-matrix connectivity
      IZG = IZS + 4*NZ               !  parameter optimization array
      IXT = IZG + 3*NZ               !  values of internal coordinates
      icmp = IXT + MAX(NIC,3*NZ)     !  composite coordinate definition
      iwt = icmp + 5*NPComp          !  primitive weights
      ipcn = iwt + NPComp            !  primitive numbers
      IEnd = ipcn + NPComp
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
      CALL SCANMAIN(NAtom,   NAtoms,  NMol,    Z(IXC),  Z(IMOL),
     $              Z(IXCG), Z(IXCM), NScan,   Z(ICT),  Z(ICN),
     $              Z(ICNM), NTrans,  Z(ITN),  Z(INQ),  Z(IUQ),
     $              NIC,     Z(iccn), Z(ityp), Z(ilst), Z(icof),
     $              Z(isav), Z(iut),  Z(ixpm), NZ,      Z(IZO),
     $              Z(IZS),  Z(IZG),  Z(IXT),  NPComp,  Z(icmp),
     $              Z(iwt),  Z(ipcn), NScr,    Z(IScr), force,
     $              inp,     Cnvgd)
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
      write(IOut,1100) IEnd,mxmem,memtot
 1100 format(/,' Potential Scan memory status:',/,
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
      SUBROUTINE SCANMAIN(NAtom,  NAtoms, NMol,   XC,     IMOL,
     $                    XCharg, XMass,  NScan,  ICTYP,  ICON,
     $                    ICNUM,  NTrans, TRANS,  NEqATM, IUNQ,
     $                    NIC,    ICNNCT, ktyp,   klist,  Coeff,
     $                    SavTOR, UT,     XPRIM,  NZ,     IGEO,
     $                    GEO,    IG,     XINT,   NPComp, IPCOMP,
     $                    PCWght, IPCNUM, NMem,   Z,      force,
     $                    inp,    Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for Potential Scan module
C  SCANMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(3*NAtom),XCharg(NAtom),XMass(NAtom),XINT(*)
      DIMENSION ktyp(NIC),klist(4,NIC),Coeff(NIC),SavTOR(NIC),
     $          UT(3*NAtoms*NIC),XPRIM(NIC)
      DIMENSION GEO(NZ,3),IGEO(NZ,4),IG(3*NZ)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans),
     $          IUNQ(NAtoms),IMOL(NAtoms),RM(3,3)
      DIMENSION ICTYP(NScan),ICON(4,NScan),ICNUM(NScan),
     $          ICNNCT(NAtoms,NAtoms)
      DIMENSION IPCOMP(5,NPComp),PCWght(NPComp),IPCNUM(NPComp)
c ............................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtom),ZSymb(NZ),VARNAM(3*NZ)
      CHARACTER*8 CZ(14*NZ+3*NAtoms)
c ............................................................
      CHARACTER GROUP*4,Coord*4,cdum*20,jobname*256
      LOGICAL Symflag,zmat,force,Cnvgd,check
C
      DIMENSION Z(NMem)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      Common /job/jobname,lenJ
C
      Data IUnit/1/
C
C
      TORAD = PI/180.0d0
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
C    scanning data
C    current energy (if available)
C    print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdscan(IUnit,zmat,Coord,I,J,K,L,NPComp,
     $            IPCOMP,PCWght,R1,R2,Step,force)
      call fdcntrl(IUnit,7,'$energy',idum)
      If(idum.EQ.0) call rdcntrl(IUnit,7,'$energy',2,idum,EC,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C .........................................................................
C
      IF(zmat) THEN
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
cc
      ELSE IF(.NOT.force) THEN
cc
       If(Coord.EQ.'comp') Then
        ICTYP(1) = 9
        GO TO 5
       EndIf
c
c -- store definition of scanned primitive
       ICON(1,1) = I
       ICON(2,1) = J
       ICON(3,1) = K
       ICON(4,1) = L
c -- what type of primitive
       If(Coord.EQ.'stre') Then
        ICTYP(1) = 1
       Else If(Coord.EQ.'bend') Then
        ICTYP(1) = 2
       Else If(Coord.EQ.'outp') Then
        ICTYP(1) = 3
       Else If(Coord.EQ.'tors') Then
        ICTYP(1) = 4
       Else
        Call nerror(4,'SCAN Module',
     $       'Unrecognizable Type of Primitive Internal to Scan',0,0)
       EndIf
c -- convert scan data to atomic units
c    (composite scan should already be in au)
       If(Coord.EQ.'stre') Then
        R1 = R1*ANTOAU
        R2 = R2*ANTOAU
        Step = Step*ANTOAU
       Else
        R1 = R1*TORAD
        R2 = R2*TORAD
        Step = Step*TORAD
       EndIf
      ENDIF
c
  5   CONTINUE
C
C  read initial coordinates from <coord> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoordF(IUnit,  NAtoms, AtSymb, XC,     NMol,
     $              IMOL,   XCharg, XMass)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Attempt to read <scan> file
C  If no <scan> file is available, this is the first
C  step of a new scan
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.scan',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      READ(40) NStep
      If(.NOT.zmat.AND..NOT.force)
     $    READ(40) NPrim,NDEG,ktyp,klist,SavTOR,UT
      CLOSE (UNIT=40,STATUS='KEEP')
      GO TO 20
C
C  Nothing on file
C  Initialize
C
 10   CONTINUE
      NStep = -1
      EC = 0.0d0
C
 20   CONTINUE
C
C  Do the Scan
C
C ---------------------------------------------------------------------
      IF(force) THEN
        CALL ForceScan(inp,    NStep,  I,      J,      R1,
     $                 R2,     Step,   NAtoms, AtSymb, XC,
     $                 EC,     Cnvgd)
      ELSE
        CALL PotScan(NStep,  NAtoms, AtSymb, XC,     NMol,
     $               IMOL,   GROUP,  NDEG,   NScan,  ICTYP,
     $               R1,     R2,     Step,   ICON,   ICNUM,
     $               NPComp, IPCOMP, PCWght, IPCNUM,
     $               NQ,     NTrans, TRANS,  NEqATM, IPRNT,
     $               NIC,    NPrim,  ICNNCT, ktyp,   klist,
     $               Coeff,  SavTOR, UT,     XPRIM,  NZ,
     $               IGEO,   GEO,    IG,     VARNAM, Coord,
     $               XINT,   EC,     NMem,   Z,      Cnvgd)
      ENDIF
C ---------------------------------------------------------------------
C
      If(Cnvgd) Then
c -- delete <scan> file
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.scan',FORM='FORMATTED')
       CLOSE (UNIT=40,STATUS='DELETE')
       RETURN
      EndIf
C
C  Another step needed
C
C  Write <scan> file
C
      check = (zmat.OR.force)
      CALL WrtSCAN(check,  NAtoms, NPrim,  NDEG,   NStep,
     $             ktyp,   klist,  SavTOR, UT)
C
C  Write necessary data to <coord>/<zmat> file
C
      IF(NZ.GT.0) THEN
       CALL WrZMAT(NZ,GEO,IG,VARNAM,.true.)         ! <coord> file DELETED
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
C
C  on first entry write number of primitives to TEXAS depository
C
       If(.NOT.force.AND.NStep.eq.0) call setival('npscan',nprim)
      ENDIF
C
      RETURN
      END
c ..............................................................................
c
      SUBROUTINE PotScan(NStep,  NAtoms, AtSymb, XC,     NMol,
     $                   IMOL,   GROUP,  NDEG,   NScan,  ICTYP,
     $                   R1,     R2,     Step,   ICON,   ICNUM,
     $                   NPComp, IPCOMP, PCWght, IPCNUM,
     $                   NQ,     NTrans, TRANS,  NEqATM, IPRNT,
     $                   NIC,    NPrim,  ICNNCT, ktyp,   klist,
     $                   Coeff,  SavTOR, UT,     XPRIM,  NZ,
     $                   IGEO,   GEO,    IG,     VARNAM, Coord,
     $                   XINT,   EC,     NMem,   Z,      Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out a potential scan
C  This is done EITHER in delocalized internal coordinates
C  OR Z-Matrix coordinates
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
C  NMol     -  number of molecules (for, e.g., cluster optimizations)
C  IMOL     -  pointers to start/end of molecules in XC array
C  GROUP    -  point group symbol
C  NDEG     -  number of non-zero delocalized internal coordinates
C              (this is GREATER THAN OR EQUAL to the real number of
C               degrees of freedom)
C  NScan    -  number of primitives to scan (should be 1 in this implementation)
C  ICTYP    -  scanned primitive type
C               1 - scan distance               stre
C               2 - scan bond angle             bend
C               3 - scan out-of-plane bend      outp
C               4 - scan dihedral angle         tors
C               9 - scan composite coordinate   comp
C  R1       -  initial value of scanned primitive (atomic units)
C  R2       -  final value of scanned primitive (atomic units)
C  Step     -  step length for scan (atomic units)
C  ICON     -  atoms involved in scanned primitive
C                IC1-IC2           distance
C                IC1-IC2-IC3       bond angle
C                IC1-IC2-IC3-IC4   all other primitives
C  ICNUM    -  indicates which primitive in list corresponds
C              to the scanned primitive
C  NPComp   -  number of primitives in composite coordinate
C  IPComp   -  primitive type and atoms involved in composite coordinate
C               IPComp(1,I) - primitive type
C               IPComp(2,I) to IPComp(5,I) - atoms in primitive
C  PCWght   -  weight of each primitive in composite constraint
C             (if no value given, assumed to be unity)
C  IPCNUM   -  indicates which primitive in list corresponds
C              to the desired composite coordinate
C  NQ       -  number of symmetry unique atoms
C  NTrans   -  number of symmetry operations
C  TRANS    -  symmetry operations as 3x3 transformation matrices
C  NEqATM   -  list of atomic equivalences under symmetry operations
C  IPRNT    -  print flag
C  NIC      -  estimate of the maximum number of primitive internal
C              coordinates (stretches, bends and torsions) in system
C              (set to zero for Z-Matrix scan)
C  NPrim    -  number of primitive internals generated
C  ICNNCT   -  atomic connectivity matrix
C  ktyp     -  integer array containing internal coordinate type
C  klist    -  list of atoms involved in each primitive
C  Coeff    -  primitive weighting factors
C  SavTOR   -  array for storing primitive torsions
C              (possible sign changes near limiting values)
C  UT       -  full set of delocalized internal coordinates
C  XPRIM    -  values of primitive internals
C  NZ       -  number of Z-matrix centers for a Z-matrix scan
C              (set to zero for scan in delocalized internals)
C  IGEO     -  Z-matrix connectivity
C  GEO      -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IG       -  array determining what to do with Z-matrix parameter
C                  0 - optimize it
C                  J - assign same value as previous (Jth) variable
C                 -J - assign same value, opposite sign
C               1000 - fixed
C  VARNAM   -  Z-Matrix variable names
C  Coord    -  Z-matrix variable to be scanned
C  XINT     -  array to hold delocalized internal coordinate or
C              Z-Matrix values
C  EC       -  current energy
C  NMem     -  amount of available scratch storage
C  Z        -  scratch array
C  Cnvgd    -  Logical flag (set to TRUE for last step of scan)
C
C
      DIMENSION XC(3,NAtoms),IMOL(NMol+1),ICNNCT(NAtoms,NAtoms)
      DIMENSION ICTYP(NScan),ICON(4,NScan),ICNUM(NScan)
      DIMENSION IPCOMP(5,NPComp),PCWght(NPComp),IPCNUM(NPComp)
      DIMENSION TRANS(9,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION ktyp(NIC),klist(4,NIC),Coeff(NIC),SavTOR(NIC),
     $          XPRIM(NIC),XINT(*),UT(*)
      DIMENSION IGEO(NZ,4),GEO(NZ,3),IG(3*NZ)
      CHARACTER*8 AtSymb(NAtoms),VARNAM(3*NZ)
      CHARACTER*4 GROUP,Coord
      LOGICAL Cnvgd
c
      Dimension Z(NMem)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,TolZero=1.0d-8)
      DATA thrsh/1.0d-7/
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      TORAD = PI/180.0d0
C
C
      IOut = ioutfil('iout')
      ICond = ioutfil('icond')
C
C .....................................................................
C
      IF(NZ.GT.0) THEN
C
C  Z-Matrix Scan
C  -------------
C ** WARNING **  All values in angstroms/degrees
C
        NVib = MAX(1,3*NZ-6)
C
C  determine scanned variable type (stretch, bend or torsion)
C
        DO 15 I=1,NVib
        IF(VARNAM(I)(1:4).EQ.Coord) THEN
          If(I.LT.NZ) Then
           ICTYP(1) = 1
          Else If(I.LT.2*NZ-2) Then
           ICTYP(1) = 2
          Else
           ICTYP(1) = 4
          EndIf
          Exit
        ENDIF
 15     CONTINUE
C
C ------------------------------------------------------------------------
C  Print Out for current step
C
        If(IPrnt.GT.0)
     $     CALL PrntZSCAN(IOut,   NStep,  NScan,  R1,     R2,
     $                    Step,   ICTYP,  Coord,  EC,     NAtoms,
     $                    AtSymb, XC)
        If(NStep.GE.0)
     $     CALL PrntZSCAN(ICond,  NStep,  NScan,  R1,     R2,
     $                    Step,   ICTYP,  Coord,  EC,     NAtoms,
     $                    AtSymb, XC)
C
C ------------------------------------------------------------------------
C
C  NEW STEP
C  --------
C
C  set value of scanned primitive
C
        NStep = NStep+1
        RVal = R1 + NStep*Step
c
        If( (Step.GT.Zero.AND.(RVal-R2).GT.TolZero) .OR.
     $      (Step.LT.Zero.AND.(R2-RVal).GT.TolZero) ) Then
         Cnvgd = .True.
         RETURN
        Else
         Cnvgd = .False.
        EndIf
C
C  get the array of internal coordinates
C
        CALL GetZVAR(NZ,NVib,GEO,XINT)
C
C  set the scanned variable in the internal coordinate array
C
        DO 18 I=1,NVib
        IF(VARNAM(I)(1:4).EQ.Coord) THEN
          If(IG(I).GE.0) Then
           XINT(I) = RVal
          Else
           XINT(I) = -RVal
          EndIf
        ENDIF
 18     CONTINUE
C
C  update Z-Matrix GEO array
C
        DO 20 I=1,NZ-1
        GEO(I+1,1) = XINT(I)
 20     CONTINUE
c
        IT = NZ-1
c
        DO 30 I=1,NZ-2
        IT = IT+1
        GEO(I+2,2) = XINT(IT)
 30     CONTINUE
c
        DO 40 I=1,NZ-3
        IT = IT+1
        GEO(I+3,3) = XINT(IT)
 40     CONTINUE
C
      ELSE
C
C  Delocalized Internal Coordinate Scan
C  ------------------------------------
C
        NVib = 3*NAtoms-6
        If(GROUP.EQ.'c*v '.OR.GROUP.EQ.'d*h ') NVib=NAtoms-1
c
        IF(NStep.EQ.-1) THEN
C
C  First Entry
C  check if scan of single variable does not break symmetry
C
          CALL SymSCAN(NAtoms, ICTYP,  ICON,   NTrans, NEqATM)
C
C  generate full set of delocalized internal coordinates
C
          NFix = 0
          NDrive = 0
          LGen = 2
          IType = 0
          ITors = 0
          IOrd = -1
          PThrsh = 170.0d0*TORAD    ! threshold for near-linear bend
          Call IZeroIT(ICNNCT,NAtoms*NAtoms)
C
C  make sure scanned variable is connected
C
          CALL CnnctSCAN(NAtoms, NScan,  ICTYP,  ICON,   NPComp,
     $                   IPCOMP, ICNNCT)
c
          IB = 1
          INB = IB + 12*NIC
          IEnd = INB + 12*NIC - 1
          JMem = NMem - IEnd
c
          CALL MakeINTC(NAtoms, AtSymb, XC,     NMol,   IMOL,
     $                  GROUP,  NVib,   NScan,  ICTYP,  R1,
     $                  ICON,   ICNUM,  NPComp, IPComp, IPCNUM,
     $                  NCon,   NFix,   IFIX,   NDrive, IDRTYP,
     $                  IDRIVE, FDRIVE, LDRIVE, NQ,     NTrans,
     $                  TRANS,  NEqATM, IPRNT,  NIC,    NPrim,
     $                  NPrim0, ICNNCT, LGen,   IType,  ITors,
     $                  CutOff, BSkal,  PThrsh, ktyp,   klist,
     $                  Z(IB),  Z(INB), Coeff,  SavTOR, NCmp,
     $                  NP1,    INT1,   UT,     XPRIM,  IGEO,
     $                  GEO,    IG,     IOrd,   JMem,   Z(IEnd),
     $                  IErr)
c
          If(IErr.NE.0) Call nerror(5,'SCAN module',
     $      'Unable to Generate Coordinates for Potential Scan',0,0)
C
C  Partial symmetry sort
C  Eliminate all delocalized internals with zero value
C  (This removes most - but not all - symmetry-breaking coordinates)
C
          CALL SymDIC(NPrim,  NVib,   XPRIM,  UT,     thrsh,
     $                XINT,   NDGG)
C
C  ** WARNING ** The number of non-zero delocalized internals is
C  returned in NDGG
C  Note that this is GREATER THAN OR EQUAL to the real number
C  of degrees of freedom, NDEG.
C
          NDEG = NDGG
C
C  Now isolate the primitive to be scanned
C  (don't isolate if only one variable, e.g., a diatomic)
C
          IF(NDEG.GT.1) THEN
            IB1 = 1
            IB2 = IB1 + NPrim
            IEnd = IB2 + (NDEG+1)*NPrim - 1
            CALL MemCHK(NMem,IEnd,7,'PotSCAN')
c
            CALL ConVEC(NDEG,   NPrim,  NScan,  ICNUM,  NPComp,
     $                  PCWght, IPCNUM, UT,     thrsh,  IPRNT,
     $                  LCON,   NCon,   Z(IB1), Z(IB2))
c
            CALL SCHMIDT(NDEG,   NPrim,  NCon,   Z(IB1), thrsh,
     $                   IPRNT,  Z(IB2), NS,     LCON,   UT,
     $                   IErr)
c
            If(IErr.NE.0) Then
              If(NS+NCon.GT.NDEG) Then
                WRITE(IOut,*) ' ** IGNORING ERROR - CONTINUING **'
                NDEG = NS+NCon
              Else
                Call nerror(6,'SCAN module',
     $            'Problems Isolating Scanned Variable in SCHMIDT',0,0)
              EndIf
            EndIf
          ENDIF

C
C  scanned primitive is isolated in the last vector in UT
C
        ENDIF
C
C ------------------------------------------------------------------------
C  Print Out for current step
C
        If(IPrnt.GT.0)
     $     CALL PrntSCAN(IOut,   NStep,  NScan,  R1,     R2,
     $                   Step,   ICTYP,  ICON,   NPComp, IPCOMP,
     $                   PCWght, EC,     NAtoms, AtSymb, XC)
        If(NStep.GE.0)
     $     CALL PrntSCAN(ICond,  NStep,  NScan,  R1,     R2,
     $                   Step,   ICTYP,  ICON,   NPComp, IPCOMP,
     $                   PCWght, EC,     NAtoms, AtSymb, XC)
C
C ------------------------------------------------------------------------
C
C  NEW STEP
C  --------
C
        NStep = NStep+1
C
C  set value of scanned primitive
C
        RVal = R1 + NStep*Step
c
        If( (Step.GT.Zero.AND.(RVal-R2).GT.TolZero) .OR.
     $      (Step.LT.Zero.AND.(R2-RVal).GT.TolZero) ) Then
         Cnvgd = .True.
         RETURN
        Else
         Cnvgd = .False.
        EndIf
C
C  get values of primitives
C
        IF(NStep.GT.0) THEN
C
C  initialize primitive weights
C
          DO 50 I=1,NPrim
          Coeff(I) = One
 50       CONTINUE
c
          CALL BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $               Coeff,  SavTOR, 1,      0,      1,
     $               BSkal,  IPRNT,  NDEG,   UT,     BJnk,
     $               XPRIM,  IErr)
        ENDIF
C
C  transform to internal coordinates
C
        CALL TranINT(NPrim,NDEG,XPRIM,UT,XINT)
C
C  set value of corresponding internal coordinate
C  (If it is a composite coordinate, it will need to be normalized)
C
        If(ICTYP(1).EQ.9) Then
          Fac = 0.0d0
          DO I=1,NPComp
          Fac = Fac + PCWght(I)**2
          EndDO
          RVal = RVal/SQRT(Fac)
        EndIf
c
        ID = 1
        IEnd = ID + NDEG
        JMem = NMem - IEnd
c
        CALL ZeroIT(Z(ID),NDEG)
        Z(ID+NDEG-1) = RVal - XINT(NDEG)
        XINT(NDEG) = RVal
C
C  if the constraint is a torsion make sure signs are the same
C
        If(ICTYP(1).EQ.4) Then
          CALL LocateTors(NPrim,ICON,ktyp,klist,I)
          SavTOR(I) = SIGN(SavTOR(I),XINT(NDEG))
          Z(ID+NDEG-1) = RVal - SavTOR(I)
        EndIf
C
C  transform back to Cartesians
C
        IGen = 0    ! delocalized internal coordinates
c
        CALL GetCART(NAtoms, AtSymb, NDEG,   NPrim,  IGen,
     $               XINT,   Z(ID),  ktyp,   klist,  Coeff,
     $               SavTOR, UT,     NTrans, TRANS,  NEqATM,
     $               BSkal,  IPRNT,  Scal,   XC,     JMem,
     $               Z(IEnd),IErr)
c
        If(IErr.NE.0.OR.Scal.NE.One) Call nerror(7,'SCAN module',
     $    'Problems with the Cartesian Back-Transformation',0,0)
      ENDIF
C
C  we are done
C
      RETURN
      END
c ..............................................................................
c
      SUBROUTINE CnnctSCAN(NAtoms, NScan,  ICTYP,  ICON,   NPComp,
     $                     IPCOMP, ICNNCT)
      IMPLICIT INTEGER(A-Z)
C
C  Ensures that scanned variable has suitable atomic connectivity
C  for successful generation of delocalized internal coordinates
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NScan   -  number of variables to scan
C             (should be 1 in this implementation)
C  ICTYP   -  scanned variable type
C  ICON    -  atoms involved in scanned variable
C  NPComp  -  number of primitives (if any) in composite coordinate
C  IPCOMP  -  type and atoms involved in composite coordinate
C  ICNNCT  -  connectivity matrix
C
      DIMENSION ICTYP(NScan),ICON(4,NScan),IPCOMP(5,NPComp),
     $          ICNNCT(NAtoms,NAtoms)
C
C
C  Loop over scanned variables
C
      IT = 1
      I = ICON(1,IT)
      J = ICON(2,IT)
      K = ICON(3,IT)
      L = ICON(4,IT)
c
      IF(ICTYP(IT).EQ.1) THEN
c
c -- stretch
        ICNNCT(I,J) = 1
        ICNNCT(J,I) = 1
c
      ELSE IF(ICTYP(IT).EQ.2) THEN
c
c -- bend
        ICNNCT(I,J) = 1
        ICNNCT(J,I) = 1
        ICNNCT(J,K) = 1
        ICNNCT(K,J) = 1
c
      ELSE IF(ICTYP(IT).EQ.3) THEN
c
c -- out-of-plane bend
        ICNNCT(I,L) = 1
        ICNNCT(L,I) = 1
        ICNNCT(J,L) = 1
        ICNNCT(L,J) = 1
        ICNNCT(K,L) = 1
        ICNNCT(L,K) = 1
c
      ELSE IF(ICTYP(IT).EQ.4) THEN
c
c -- torsion
        ICNNCT(I,J) = 1
        ICNNCT(J,I) = 1
        ICNNCT(J,K) = 1
        ICNNCT(K,J) = 1
        ICNNCT(K,L) = 1
        ICNNCT(L,K) = 1
c
      ELSE IF(ICTYP(IT).EQ.9) THEN
c
c -- composite coordinate
        DO IP=1,NPComp
        ITYP = IPCOMP(1,IP)
        I = IPCOMP(2,IP)
        J = IPCOMP(3,IP)
        K = IPCOMP(4,IP)
        L = IPCOMP(5,IP)
c
        If(ITYP.EQ.1) Then
c
c -- stretch
          ICNNCT(I,J) = 1
          ICNNCT(J,I) = 1
c
        Else If(ITYP.EQ.2) Then
c
c -- bend
          ICNNCT(I,J) = 1
          ICNNCT(J,I) = 1
          ICNNCT(J,K) = 1
          ICNNCT(K,J) = 1
c
        Else If(ITYP.EQ.3) Then
c
c -- out-of-plane bend
          ICNNCT(I,L) = 1
          ICNNCT(L,I) = 1
          ICNNCT(J,L) = 1
          ICNNCT(L,J) = 1
          ICNNCT(K,L) = 1
          ICNNCT(L,K) = 1
c
        Else IF(ITYP.EQ.4) Then
c
c -- torsion
          ICNNCT(I,J) = 1
          ICNNCT(J,I) = 1
          ICNNCT(J,K) = 1
          ICNNCT(K,J) = 1
          ICNNCT(K,L) = 1
          ICNNCT(L,K) = 1
c
        EndIf
        EndDO
c
      ENDIF
C
      RETURN
      END
c ..............................................................................
c
      SUBROUTINE LocateTors(NPrim,ITors,ktyp,klist,IT)
      IMPLICIT INTEGER(A-Z)
C
C  locates a given torsion in the primitive list
C
C  ARGUMENTS
C
C  NPrim   -  number of primitives
C  ITors   -  atoms in torsion
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  IT      -  on exit, position of torsion in primitive list
C
      DIMENSION ITors(4),ktyp(NPrim),klist(4,NPrim)
C
      DO 10 I=1,NPrim
      IF(ktyp(I).EQ.4) THEN
        If((klist(1,I).EQ.ITors(1).AND.klist(2,I).EQ.ITors(2).AND.
     $      klist(3,I).EQ.ITors(3).AND.klist(4,I).EQ.ITors(4)) .OR.
     $     (klist(1,I).EQ.ITors(4).AND.klist(2,I).EQ.ITors(3).AND.
     $      klist(3,I).EQ.ITors(2).AND.klist(4,I).EQ.ITors(1))) Then
          IT = I
          GO TO 20
        EndIf
      ENDIF
 10   CONTINUE
C
C  should never get here
C
      Call nerror(8,'SCAN module',
     $    'Cannot Locate Scanned Torsion in Primitive List!',0,0)
c
 20   RETURN
      END
c ..............................................................................
c
      SUBROUTINE PrntSCAN(IOut,   NStep,  NScan,  R1,     R2,
     $                    Step,   ICTYP,   IC,    NPComp, IPCOMP,
     $                    PCWght, EC,     NAtoms, AtSymb, XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out the scanned parameter values and current geometry
C
C  ARGUMENTS
C
C  IOut    -  output unit
C  NStep   -  step number
C  NScan   -  number of primitives to scan
C             (should be 1 in this implementation)
C  R1      -  starting value for scanned primitive (atomic units)
C  R2      -  final value for scanned primitive (atomic units)
C  Step    -  scanning step length (atomic units)
C  ICTYP   -  primitive type
C  IC      -  atoms defining primitive
C  NPComp  -  number of primitives in composite coordinate
C  IPCOMP  -  type and atoms involved in composite coordinate
C  PCWght  -  weight of each primitive in composite coordinate
C  EC      -  current energy
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C
      DIMENSION IC(4), XC(3,NAtoms)
      DIMENSION IPComp(5,NPComp), PCWght(NPComp)
      CHARACTER*8 AtSymb(NAtoms)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      AUTOAN = 1.0d0/ANTOAU
      TOANG = 180.0d0/PI
C
C
      WRITE(IOut,1000)
      If(NStep.EQ.-1) Then
       WRITE(IOut,1001)
       RVal = R1
      Else
       WRITE(IOut,1002) NStep+1
       RVal = R1 + NStep*Step
      EndIf
c
      IF(ICTYP.EQ.1) THEN
C
C  scanned distance
C
        Val = RVal*AUTOAN
        WRITE(IOut,1100) IC(1),IC(2),R1*AUTOAN,R2*AUTOAN,Step*AUTOAN
        If(NStep.GT.-1) WRITE(IOut,1010) Val,EC
cc
      ELSE IF(ICTYP.EQ.2) THEN
C
C  scanned bend
C
        Val = RVal*TOANG
        WRITE(IOut,1200) IC(1),IC(2),IC(3),R1*TOANG,R2*TOANG,Step*TOANG
        If(NStep.GT.-1) WRITE(IOut,1010) Val,EC
cc
      ELSE IF(ICTYP.EQ.3) THEN
C
C  scanned out-of-plane bend
C
        Val = RVal*TOANG
        WRITE(IOut,1300) IC(1),IC(2),IC(3),IC(4),R1*TOANG,R2*TOANG,
     $                   Step*TOANG
        If(NStep.GT.-1) WRITE(IOut,1010) Val,EC
cc
      ELSE IF(ICTYP.EQ.4) THEN
C
C  scanned dihedral angle
C
        Val = RVal*TOANG
        WRITE(IOut,1400) IC(1),IC(2),IC(3),IC(4),R1*TOANG,R2*TOANG,
     $                   Step*TOANG
        If(NStep.GT.-1) WRITE(IOut,1010) Val,EC
cc
      ELSE IF(ICTYP.EQ.9) THEN
C
C  scanned composite coordinate
C
        WRITE(IOut,1700) R1,R2,Step
        If(NStep.GT.-1) WRITE(IOut,1010) RVal,EC
c
        WRITE(IOut,2000)
        val = 0.0d0
c
c -- get values of individual primitives
        DO IP=1,NPComp
        ITYP = IPCOMP(1,IP)
        Wght = PCWght(IP)
        I = IPCOMP(2,IP)
        J = IPCOMP(3,IP)
        K = IPCOMP(4,IP)
        L = IPCOMP(5,IP)
c
        IF(ITYP.EQ.1) THEN
cc
C  distance constraint
C
          XIJ = XC(1,I) - XC(1,J)
          YIJ = XC(2,I) - XC(2,J)
          ZIJ = XC(3,I) - XC(3,J)
          R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
          val = val+R*Wght
          R = R*AUTOAN
c
          WRITE(IOut,2100) I,J,R,Wght
cc
        ELSE IF(ITYP.EQ.2) THEN
cc
C  bond angle constraint
C
          CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
          val = val+Th*Wght
          Th = Th*TOANG
c
          WRITE(IOut,2200) I,J,K,Th,Wght
cc
        ELSE IF(ITYP.EQ.3) THEN
cc
C  out-of-plane bend constraint
C
          CALL OutpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
          val = val+Th*Wght
          Th = Th*TOANG
c
          WRITE(IOut,2300) I,J,K,L,Th,Wght
cc
        ELSE IF(ITYP.EQ.4) THEN
cc
C  dihedral angle constraint
C
          CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
          val = val+Th*Wght
          Th = Th*TOANG
c
          WRITE(IOut,2400) I,J,K,L,Th,Wght
cc
        ENDIF
        EndDO
      ENDIF
c
      CALL PrntCAR(IOut,0,NAtoms,AtSymb,XC)
      WRITE(IOut,1000)
C
      RETURN
c
 1000 FORMAT(1X,70('-'))
 1001 FORMAT(' STARTING POTENTIAL SCAN')
 1002 FORMAT(' POTENTIAL SCAN   Step: ',I4)
 1010 FORMAT(/,'  Current value: ',F10.4,'  Energy is ',F18.9)
 1100 FORMAT('  Scanning Distance:     ',2I4,12X,/,
     $     '  From ',F8.4,' to ',F8.4,'  Angstrom  step length: ',F8.4)
 1200 FORMAT('  Scanning Angle:        ',3I4,8X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1300 FORMAT('  Scanning Out-of-Plane: ',4I4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1400 FORMAT('  Scanning Dihedral:     ',4I4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1500 FORMAT('  Scanning Linear Plane: ',4I4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1600 FORMAT('  Scanning Linear Perp.: ',4I4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1700 FORMAT('  Scanning composite coordinate',/,
     $     '  From ',F8.3,' to ',F8.3,'  au        step length: ',F8.3)
 2000 FORMAT('  Composite Coordinate',18X,'value',8X,'weight')
 2100 FORMAT('   Distance:     ',2I4,12X,F9.6,4X,F9.6)
 2200 FORMAT('   Angle:        ',3I4,8X,F9.3,4X,F9.6)
 2300 FORMAT('   Out-of-Plane: ',4I4,4X,F9.3,4X,F9.6)
 2400 FORMAT('   Dihedral:     ',4I4,4X,F9.3,4X,F9.6)
c
      END
c ..............................................................................
c
      SUBROUTINE PrntZSCAN(IOut,   NStep,  NScan,  R1,     R2,
     $                     Step,   ICTYP,  Coord,  EC,     NAtoms,
     $                     AtSymb, XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out the scanned Z-matrix parameter value
C
C  ARGUMENTS
C
C  IOut    -  output unit
C  NStep   -  step number
C  NScan   -  number of primitives to scan
C             (should be 1 in this implementation)
C  R1      -  starting value for scanned primitive (angstrom/degrees)
C  R2      -  final value for scanned primitive (angstrom/degrees)
C  Step    -  scanning step length (angstrom/degrees)
C  ICTYP   -  Z-Matrix parameter type
C  Coord   -  Z-Matrix parameter name
C  EC      -  current energy
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C
      CHARACTER*4 Coord
C
C
      WRITE(IOut,1000)
      If(NStep.EQ.-1) Then
       WRITE(IOut,1001)
       RVal = R1
      Else
       WRITE(IOut,1002) NStep+1
       RVal = R1 + NStep*Step
      EndIf
c
      IF(ICTYP.EQ.1) THEN
C
C  scanned distance
C
        WRITE(IOut,1100) Coord,R1,R2,Step
        If(NStep.GT.-1) WRITE(IOut,1010) RVal,EC
cc
      ELSE IF(ICTYP.EQ.2) THEN
C
C  scanned bend
C
        WRITE(IOut,1200) Coord,R1,R2,Step
        If(NStep.GT.-1) WRITE(IOut,1010) RVal,EC
cc
      ELSE IF(ICTYP.EQ.3) THEN
C
C  scanned out-of-plane bend
C
        WRITE(IOut,1300) Coord,R1,R2,Step
        If(NStep.GT.-1) WRITE(IOut,1010) RVal,EC
cc
      ELSE IF(ICTYP.EQ.4) THEN
C
C  scanned dihedral angle
C
        WRITE(IOut,1400) Coord,R1,R2,Step
        If(NStep.GT.-1) WRITE(IOut,1010) RVal,EC
cc
      ENDIF
c
      CALL PrntCAR(IOut,0,NAtoms,AtSymb,XC)
      WRITE(IOut,1000)
C
      RETURN
c
 1000 FORMAT(1X,70('-'))
 1001 FORMAT(' STARTING POTENTIAL SCAN')
 1002 FORMAT(' POTENTIAL SCAN   Step: ',I4)
 1010 FORMAT(/,'  Current value: ',F10.4,'  Energy is ',F18.9)
 1100 FORMAT('  Scanning Distance:     ',A4,4X,/,
     $     '  From ',F8.4,' to ',F8.4,'  Angstrom  step length: ',F8.4)
 1200 FORMAT('  Scanning Angle:        ',A4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1300 FORMAT('  Scanning Out-of-Plane: ',A4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
 1400 FORMAT('  Scanning Dihedral:     ',A4,4X,/,
     $     '  From ',F8.3,' to ',F8.3,'  degrees   step length: ',F8.3)
c
      END
c .........................................................................
c
      SUBROUTINE SymDIC(NPrim,  NVib,   XPRIM,  UT,     thrsh,
     $                  XINT,   NDEG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Removes delocalized internal coordinates with zero value
C  Coordinates are typically only zero by symmetry
C  This will eliminate most, but not necessarily all,
C  symmetry-breaking internal coordinates
C
C  ARGUMENTS
C
C  NPrim   -  number of primitives
C  NVib    -  number of delocalized internals
C  XPRIM   -  values of primitive internals
C  UT      -  full set of delocalized internal coordinates
C  thrsh   -  zero threshold
C  XINT    -  values of delocalized internals
C
C  on exit
C
C  NDEG    -  number of non-zero internals remaining
C
C
      DIMENSION XPRIM(NPrim),UT(NPrim,NVib),XINT(NVib)
C
C
C  transform to internal coordinates
C
      CALL TranINT(NPrim,NVIB,XPRIM,UT,XINT)
C
C  now eliminate zero coordinates
C
      NDEG = 0
      DO 10 I=1,NVib
      Val = Abs(XINT(I))
      IF(Val.GT.thrsh) THEN
       NDEG = NDEG+1
       CALL CpyVEC(NPrim,UT(1,I),UT(1,NDEG))
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c .........................................................................
c
      SUBROUTINE SymSCAN(NAtoms, ICTYP,  ICON,   NTrans, NEqATM)
      IMPLICIT INTEGER(A-H,O-Z)
C
C
C  Checks whether the scanned parameter will preserve molecular symmetry.
C  If not, prints out the symmetry-related variables and Exits
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  ICTYP   -  coordinate type
C              1 - scanned distance
C              2 - scanned bond angle
C              3 - scanned out-of-plane-bend
C              4 - scanned dihedral angle
C              9 - composite coordinate (not checked)
C  ICON    -  atoms involved in scanned parameter
C               ID1-ID2           scanned distance
C               ID1-ID2-ID3       scanned bond angle
C               ID1-ID2-ID3-ID4   scanned out-of-plane bend or dihedral
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C
      DIMENSION ICON(4)
      DIMENSION NEqATM(NAtoms,NTrans)
C
C
      If(NTrans.EQ.1) RETURN     ! there is no symmetry to conserve
      If(ICTYP.EQ.9) RETURN      ! don't check composite coordinate
c
      IOut = ioutfil('iout')
      IErr = 0
C
      IF(ICTYP.EQ.1) THEN
cc
C  scanning a distance (I-J)
C
       I = ICON(1)
       J = ICON(2)
C
C  check for equivalent atoms
C  these must also be scanned for symmetry
C  to be maintained
C
       DO 10 IOP=2,NTrans
       IT = NEqATM(I,IOP)
       JT = NEqATM(J,IOP)
       If( (IT.EQ.I.AND.JT.EQ.J).OR.
     $     (JT.EQ.I.AND.IT.EQ.J) ) GO TO 10
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1000) IT,JT,I,J
C
 10    CONTINUE
cc
      ELSE IF(ICTYP.EQ.2) THEN
cc
C  scanning an angle (I-J-K)
C
       I = ICON(1)
       J = ICON(2)
       K = ICON(3)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C
       DO 20 IOP=2,NTrans
       IT = NEqATM(I,IOP)
       JT = NEqATM(J,IOP)
       KT = NEqATM(K,IOP)
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K).OR.
     $     (KT.EQ.I.AND.JT.EQ.J.AND.IT.EQ.K) ) GO TO 20
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1100) IT,JT,KT,I,J,K
C
 20    CONTINUE
cc
      ELSE IF(ICTYP.EQ.3) THEN
cc
C  scanning an out-of-plane bend (I-J-K-L)
C
       I = ICON(1)
       J = ICON(2)
       K = ICON(3)
       L = ICON(4)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C
       DO 30 IOP=2,NTrans
       IT = NEqATM(I,IOP)
       JT = NEqATM(J,IOP)
       KT = NEqATM(K,IOP)
       LT = NEqATM(L,IOP)
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K.AND.LT.EQ.L).OR.
     $     (IT.EQ.I.AND.KT.EQ.J.AND.JT.EQ.K.AND.LT.EQ.L) ) GO TO 30
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1200) IT,JT,KT,LT,I,J,K,L
C
 30    CONTINUE
cc
      ELSE IF(ICTYP.EQ.4) THEN
cc
C  scanninmg a dihedral angle (I-J-K-L)
C
       I = ICON(1)
       J = ICON(2)
       K = ICON(3)
       L = ICON(4)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C
       DO 40 IOP=2,NTrans
       IT = NEqATM(I,IOP)
       JT = NEqATM(J,IOP)
       KT = NEqATM(K,IOP)
       LT = NEqATM(L,IOP)
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K.AND.LT.EQ.L).OR.
     $     (LT.EQ.I.AND.KT.EQ.J.AND.JT.EQ.K.AND.IT.EQ.L) ) GO TO 40
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1300) IT,JT,KT,LT,I,J,K,L
C
 40    CONTINUE
cc
      ENDIF
C
      IF(IErr.NE.0) THEN
       WRITE(IOut,1400)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(' Distances ',2I4,8X,' and ',2I4,8X,' Symmetry-Related')
 1100 FORMAT(' Angles    ',3I4,4X,' and ',3I4,4X,' Symmetry-Related')
 1200 FORMAT(' OOP bends ',4I4,   ' and ',4I4,   ' Symmetry-Related')
 1300 FORMAT(' Dihedrals ',4I4,   ' and ',4I4,   ' Symmetry-Related')
 1400 FORMAT(/,2X,'***ERROR*** Potential Scan Will Break Symmetry',/,
     $     '      Either scan a composite coordinate or distort your',/,
     $     '      starting geometry and switch off symmetry',/)
c
      END
c .........................................................................
c
      SUBROUTINE WrtSCAN(check,  NAtoms, NPrim,  NDEG,   NStep,
     $                   ktyp,   klist,  SavTOR, UT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  writes the <scanchk> file for the potential scan
C
      DIMENSION ktyp(NPrim),klist(4,NPrim),SavTOR(NPrim)
      DIMENSION UT(3*NAtoms*NPrim)
      Character*256 jobname
      Logical check
c
      Common /job/jobname,lenJ
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.scan',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(40) NStep
      If(.NOT.check) WRITE(40) NPrim,NDEG,ktyp,klist,SavTOR,UT
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
      END
c .........................................................................
c
      SUBROUTINE ForceScan(inp,    NStep,  I,      J,      R1,
     $                     R2,     Step,   NAtoms, AtSymb, XC,
     $                     EC,     Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine does a scan with respect to an applied external force.
C  It is done by modifying the input file at each step of the scan (NOT
C  during any optimization step) by embedding the applied force directly
C  into the input file.
C
C  ARGUMENTS
C
C  inp     -  unit number of input file
C  NStep   -  step number of force scan
C  I       -  first atom involved in external force
C  J       -  second atom involved in external force
C  R1      -  starting value for applied external force
C  R2      -  ending value for applied external force
C  Step    -  applied force increment
C  Cnvgd   -  Logical flag (set to TRUE for last step of scan)
C
C  ** IMPORTANT **
C     If applied force is +ve, assumed PULL
C     If applied force is -ve, assumed PUSH
C
C
      REAL*8 XC(3,NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
c
      Character jobname*256,Char*80,Direction*4
      Logical Cnvgd,found
      Common /job/jobname,lenJ
      Parameter (Zero=0.0d0)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      IOut = ioutfil('iout')
      ICond = ioutfil('icond')
c
      Direction = 'PULL'
      NStep = NStep+1
      force = R1 + NStep*Step
      If(force.LT.Zero) Direction = 'PUSH'
C
C  get the distance between the two atoms
C
      dist = SQRT( (XC(1,I)-XC(1,J))**2 +
     $             (XC(2,I)-XC(2,J))**2 +
     $             (XC(3,I)-XC(3,J))**2 )
      dist = dist/ANTOAU      ! convert to angstroms
c
      If( (Step.GT.Zero.AND.force.GT.R2) .OR.
     $    (Step.LT.Zero.AND.force.LT.R2) ) Then
       Cnvgd = .True.
      Else
       Cnvgd = .False.
      EndIf
c
      WRITE(IOut,1000)
      WRITE(ICond,1000)
      If(NStep.EQ.0) Then
       WRITE(IOut,1001)
       WRITE(ICond,1001)
      Else
       WRITE(IOut,1002) NStep
       WRITE(ICond,1002) NStep
      EndIf
      WRITE(IOut,1100) I,J,R1,R2,Step
      WRITE(ICond,1100) I,J,R1,R2,Step
      If(NStep.GT.0) Then
       WRITE(IOut,1010) force-step,EC
       WRITE(ICond,1010) force-step,EC
      EndIf
      WRITE(IOut,1020) I,J,dist
      WRITE(ICond,1020) I,J,dist
c
      CALL PrntCAR(IOut,0,NAtoms,AtSymb,XC)
      CALL PrntCAR(ICond,0,NAtoms,AtSymb,XC)
      WRITE(IOut,1000)
      WRITE(ICond,1000)
c
      If(Cnvgd) RETURN
C
C  modify input file to include applied force
C
  5   CONTINUE
      READ(inp,900) Char
      call lowercas(Char,4)
      If(Char(1:4).NE.'forc') GO TO 5
C
C  At this point found FORCE line
C  copy input file beyond FORCE line to unit 41
C
      OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
c -- don't copy any external force added on a previous step
      READ(inp,900) Char
      call lowercas(Char,6)
      If(Char(1:6).EQ.'$force') Then
       READ(inp,900) Char
       READ(inp,900) Char
       found = .True.
      Else
       backspace inp
       found = .False.
      EndIf
c
      line = 0
 10   CONTINUE
      READ(inp,900,End=20) Char
      line = line+1
      WRITE(41,900) Char
      GO TO 10
c
 20   CONTINUE
C
C  add external force data to input file
C
      If(found) line = line+3
      DO L=1,line+1
      backspace inp
      EndDO
c
      WRITE(inp,'(a)') '$force'
      WRITE(inp,910) I,J,force,Direction
      WRITE(inp,'(a)') '$endforce'
C
C  copy back rest of file
C
      REWIND 41
 30   CONTINUE
      READ(41,900,END=40) Char
      call rmblan(Char,80,Len)
      If(Len.GT.0) WRITE(inp,'(a)') Char(1:Len)
      GO TO 30
c
 40   CONTINUE
      CLOSE (UNIT=41,STATUS='DELETE')
C
C  now reposition read head in input file
C
      REWIND inp
 50   CONTINUE
      READ(inp,900) Char
      call lowercas(Char,4)
      If(Char(1:4).NE.'scan') GO TO 50
c
      RETURN
c
 900  Format(A80)
 910  Format(2I4,2X,F10.6,2X,A4)
 1000 FORMAT(1X,70('-'))
 1001 FORMAT(' STARTING EXTERNAL FORCE SCAN')
 1002 FORMAT(' EXTERNAL FORCE SCAN   Step: ',I4)
 1010 FORMAT(/,'  Current value: ',F10.4,'  Energy is ',F18.9)
 1020 FORMAT('  Distance between atoms ',I4,' and ',I4,' is ',F10.6,
     $       '  Angstroms')
 1100 FORMAT('  Applied External force between atoms: ',2I4,/,
     $       '  From ',F8.4,' to ',F8.4,'  au  step length: ',F8.4)
c
      END
