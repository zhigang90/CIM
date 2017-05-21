c ==================================================================
c  GEOMETRY MODULE            JB   August 1997
c ==================================================================
c
      SUBROUTINE GEOMETRY(NMem,Z,NAtoms,NDum)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** MOLECULAR GEOMETRY PROGRAM **
C
C  This Program reads the input geometry in Cartesian
C  coordinates from the <coord> file and optionally
C  determines the molecular point group symmetry.
C
C  The Cartesian coordinates are repositioned into a
C  centre-of-mass frame (with all atoms assumed to have
C  unit mass). The molecular point group symmetry is then
C  determined and symmetry data written to the <sym> file.
C  This may result in a reorientation of the Cartesian axes.
C  ............................................................
C
      DIMENSION Z(NMem),Field(3)
      Character cdum*20
      Character*256 jobname
      Logical FFlag
C
      PARAMETER (MSymOP=120)     !  maximum number of symmetry operations
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
      PARAMETER (Zero=0.0d0)
c
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    total number of atoms (including dummies)
C    total number of molecules
C    external electric field (if there is one)
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call rdcntrl(IUnit,5,'$nmol',1,NMol,dum,cdum)
      call RdField(IUnit,Field,FFlag)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C .............................................................
C  If there is an external field we are going to MIMIC its
C  effects (for symmetry purposes) by using point charges
C
      NFDum = 0
      If(FFlag) Then
        If(Field(1).NE.Zero) NFDum = NFDum + 2
        If(Field(2).NE.Zero) NFDum = NFDum + 2
        If(Field(3).NE.Zero) NFDum = NFDum + 2
      EndIf
      NAtoms = NAtoms + NFDum
C .............................................................
C
C
C  Now get the memory
C
      NScr = MAX(9*NAtoms + 7*MSymOP,9*NAtoms**2 + 9)
      IMem = 9*NAtoms + 9*MSymOP + NAtoms*MSymOP + NMol+1 + NScr
c
cc      CALL falloc(Z(0),8*IMem,iptr,IErr)
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,8,'GEOMETRY')
      CALL ZeroIT(Z,IMem)      ! clear the memory
C
C  Allocate memory pointers
C
      IXC = iptr                    !  reoriented geometry
      IXO = IXC + 3*NAtoms          !  input geometry
      IML = IXO + 3*NAtoms          !  pointer to molecules
      ICH = IML + NMol+1            !  atomic charges
      ICM = ICH + NAtoms            !  atomic masses
      IUQ = ICM + NAtoms            !  list of symmetry unique atoms
      ITN = IUQ + NAtoms            !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*MSymOP          !  list of atomic equivalences
      IScr = INQ + NAtoms*MSymOP    !  general scratch storage
      IEnd = IScr + NScr - 1
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,8,'GEOMETRY')
C
      NTrans = MSymOP
      NAtoms = NAtoms - NFdum       !  restore original value
C
C  ----------------------------------------------------------------------
C
      CALL GEOMAIN(NAtoms,  NTrans,  Z(IXO),  Z(IXC),  NMol,
     $             Z(IML),  Z(ICH),  Z(ICM),  Z(IUQ),  Z(ITN),
     $             Z(INQ),  NFDum,   Field,   NDum,    NScr,
     $             Z(IScr))
C
C  ----------------------------------------------------------------------
C
C  free memory used in this routine
C
cc      CALL ffree(Z(iptr))
C
C  Exit procedure
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE GEOMAIN(NAtoms, NTrans, XOld,   XC,     NMol,
     $                   IMOL,   XCharg, XMass,  IUNQ,   TRANS,
     $                   NEqATM, NFdum,  Field,  Ndum,   IMem,   Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for symmetry determination
C  GEOMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION XC(3,NAtoms),XOld(3,NAtoms),XCharg(NAtoms),XMass(NAtoms)
      DIMENSION TRANS(9,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION IMOL(NMol+1),IUNQ(NAtoms),Field(3),RM(3,3)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms+NFdum),Symb(NAtoms)
c ..................................................
      CHARACTER*4 GROUP,group1
      Character cdum*20
      Character*256 jobname
      LOGICAL Symflag,adjst,com,Error,rotate
C
      DIMENSION Z(IMem)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=1.0d-5)
      PARAMETER (IUnit=1)
C
      Data IPRNT/2/                ! default print flag
c
      Common /job/jobname,lenJ
cc
      iout=igetival('iout')
      icond=igetival('icond')
C
c -- set up constants for conversion to atomic units etc..
      call setconst
C
C
      adjst = .TRUE.          ! default is to reorient to symmetry axes
      com   = .TRUE.          ! default is to use "centre-of-mass" frame
      RThrsh = Zero           ! no distance printout
      Error = .FALSE.         ! flag for potential symmetry-breaking dummy atoms
C
C ...............................................................
C  Read from the <control> file
C    symmetry threshold
C    print flag
C    number of dummy atoms (if any)
C    threshold for distance printing (if any)
C    axes reorientation flag
C    cms flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,14,'$sym_threshold',2,idum,thrsym,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      call fdcntrl(IUnit,18,'$dist_print_thresh',idum)
      If(idum.eq.0) then
       call rdcntrl(IUnit,18,'$dist_print_thresh',2,idum,RThrsh,cdum)
      endIf
      call fdcntrl(IUnit,9,'$noorient',idum)
      If(idum.eq.0) adjst = .false.
      call fdcntrl(IUnit,6,'$nocms',idum)
      If(idum.eq.0) com = .false.
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(.NOT.adjst.AND.IPRNT.GT.0) WRITE(IOut,1000)
      If(.NOT.com.AND.IPRNT.GT.0) WRITE(IOut,1100)
C
C  Read the atomic symbols, coordinates, charges & masses
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call RdCoordF(IUnit,  NAtoms, Symb,   XOld,   NMol,
     $              IMOL,   XCharg, XMass)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  initialize the rotation matrix
C
      CALL SetDiagMat(3,One,RM)
C
C  Do the symmetry analysis
C
      Symflag = thrsym.gt.Zero
      if(Natoms.eq.1)Symflag=.false. ! no symmetry in atomic calculation
C
C  Check for incompatibility (if symmetry MUST use cms frame)
C
      If(.NOT.com.AND.Symflag) Call nerror(7,'GEOMETRY module',
     $ 'MUST Use "Centre-of-Mass" Axis System With Symmetry',0,0)
C
C
C  **NOTE**
C  NAtoms  -   total number of atomic centres (real + ALL dummies)
C  NAtom   -   number of real atoms + charged dummies
C  MAtom   -   number of real atoms only
C
      If(.NOT.Symflag) NFdum=0      ! no need if no symmetry
      NAtom = NAtoms-Ndum2
      MAtom = NAtom-Ndum1
      Ndum  = Ndum1+Ndum2
C
C  Are the masses defined?
C  If not, use the default values
C
      DO 5 I=1,MAtom
      Call nugrep(Symb(I)(1:2),INum)
      If(XMass(I).EQ.Zero) CALL DefMASS(1,INum,XMass(I))
 5    CONTINUE
C
C ----------------------------------------------------------------
C
      IF(Symflag) THEN
C
C  Get the actual atomic symbols from the user-defined atom symbols
C  ** atoms may be treated differently depending on symbol **
C     (see comments inside subroutine <SymbolM>)
C
       CALL SymbolM(NAtoms, Symb, AtSymb)
C
C  .............................................................
       IF(NFdum.GT.0) THEN
C
C  mimic field using charged dummy atoms
C  first move type 2 dummy atoms to end
C
        DO 10 I=1,Ndum2
        AtSymb(NAtom+NFdum+I) = AtSymb(NAtom+I)
        XOld(1,NAtom+NFdum+I) = XOld(1,NAtom+I)
        XOld(2,NAtom+NFdum+I) = XOld(2,NAtom+I)
        XOld(3,NAtom+NFdum+I) = XOld(3,NAtom+I)
 10     CONTINUE
C
C  now insert type 1 dummy atoms to represent field
C
        KK = 0
        If(Field(1).NE.Zero) Then
          KK = KK+1
          AtSymb(NAtom+KK) = 'Q'
          XOld(1,NAtom+KK) = 57.0d0
          XOld(2,NAtom+KK) = Zero
          XOld(3,NAtom+KK) = Zero
          KK = KK+1
          AtSymb(NAtom+KK) = 'X'
          XOld(1,NAtom+KK) = -57.0d0
          XOld(2,NAtom+KK) = Zero
          XOld(3,NAtom+KK) = Zero
        EndIf
        If(Field(2).NE.Zero) Then
          KK = KK+1
          AtSymb(NAtom+KK) = 'Q1'
          If(Field(2).EQ.Field(1)) AtSymb(NAtom+KK) = 'Q'
          If(Field(2).EQ.-Field(1)) AtSymb(NAtom+KK) = 'X'
          XOld(1,NAtom+KK) = Zero
          XOld(2,NAtom+KK) = 57.0d0
          XOld(3,NAtom+KK) = Zero
          KK = KK+1
          AtSymb(NAtom+KK) = 'X1'
          If(Field(2).EQ.Field(1)) AtSymb(NAtom+KK) = 'X'
          If(Field(2).EQ.-Field(1)) AtSymb(NAtom+KK) = 'Q'
          XOld(1,NAtom+KK) = Zero
          XOld(2,NAtom+KK) = -57.0d0
          XOld(3,NAtom+KK) = Zero
        EndIf
        If(Field(3).NE.Zero) Then
          KK = KK+1
          AtSymb(NAtom+KK) = 'Q2'
          If(Field(3).EQ.Field(1)) AtSymb(NAtom+KK) = 'Q'
          If(Field(3).EQ.-Field(1)) AtSymb(NAtom+KK) = 'X'
          If(Field(3).EQ.Field(2))
     $       AtSymb(NAtom+KK) = AtSymb(NAtom+KK-2)
          If(Field(3).EQ.-Field(2))
     $       AtSymb(NAtom+KK) = AtSymb(NAtom+KK-1)
          XOld(1,NAtom+KK) = Zero
          XOld(2,NAtom+KK) = Zero
          XOld(3,NAtom+KK) = 57.0d0
          KK = KK+1
          AtSymb(NAtom+KK) = 'X2'
          If(Field(3).EQ.Field(1)) AtSymb(NAtom+KK) = 'X'
          If(Field(3).EQ.-Field(1)) AtSymb(NAtom+KK) = 'Q'
          If(Field(3).EQ.Field(2))
     $       AtSymb(NAtom+KK) = AtSymb(NAtom+KK-2)
          If(Field(3).EQ.-Field(2))
     $       AtSymb(NAtom+KK) = AtSymb(NAtom+KK-3)
          XOld(1,NAtom+KK) = Zero
          XOld(2,NAtom+KK) = Zero
          XOld(3,NAtom+KK) = -57.0d0
        EndIf
       ENDIF
C  .............................................................
C
C  determine molecular symmetry
C  for symmetry purposes, type 2 dummy atoms (uncharged)
C  should be ignored
C
C  .............................................................
C
       KAtom = NAtom + NFdum
       CALL SYMTRY(KAtom,  AtSymb, XOld,   thrsym, adjst,
     $             IPRNT,  IMem,   Z,      XC,     RM,
     $             GROUP,  NTrans, TRANS,  NEqATM, NDEG)
C  .............................................................
C
C  ** Symmetry Check for atomic charges **
C
       CALL SymCHK(IOut,KAtom,NAtom,NTrans,NEqATM,XCharg)
C
C  reorient any type 2 dummy atoms
C
       If(Ndum2.GT.0) Then
        DO 20 I=KAtom+1,NAtoms+NFdum
        DO 20 J=1,3
        XC(J,I) = RM(J,1)*XOld(1,I) + RM(J,2)*XOld(2,I) +
     $            RM(J,3)*XOld(3,I)
 20     CONTINUE
       EndIf
C  .............................................................
      ELSE
C
C  simply put molecule in CMS frame
C
       CALL CpyVEC(3*NAtoms,XOld,XC)
       If(com) CALL CMS(NAtoms,X,Y,Z,XC)
C
C  set point group to C1
C
       GROUP = 'c1  '
       NTrans = 1
       NDEG = 3*MAtom-6
       If(MAtom.EQ.2) NDEG = 1        ! simple linear
       If(MAtom.EQ.1) NDEG = 0        ! single atom
       KAtom = NAtom
C
C  set unit symmetry operation
C
       CALL SetDiagMat(3,One,TRANS)
c
       DO 30 IAtm=1,MAtom
       NEqATM(IAtm,1) = IAtm
 30    CONTINUE
      ENDIF
C
C  ----------------------------------------------------------------
C  ** ONLY PUT INFORMATION ABOUT REAL ATOMS ON <sym> FILE **
C
      If(Ndum1+NFdum.GT.0) Then
        CALL ICpyVEC(KAtom*NTrans,NEqATM,Z)
        CALL SymFilter(KAtom,MAtom,NTrans,Z,NEqATM)
      EndIf
C
C  .................................................................
C  remove field dummy atoms
C
      IF(NFdum.GT.0) THEN
C
C   restore type 2 dummy atoms
C
        DO 40 I=1,Ndum2
        AtSymb(NAtom+I) = AtSymb(NAtom+NFdum+I)
        XOld(1,NAtom+I) = XOld(1,NAtom+NFdum+I)
        XOld(2,NAtom+I) = XOld(2,NAtom+NFdum+I)
        XOld(3,NAtom+I) = XOld(3,NAtom+NFdum+I)
 40     CONTINUE
      ENDIF
C  .................................................................
C
C  determine the number of degrees of freedom
C
      CALL NumDEG(MAtom,NTrans,TRANS,GROUP,NEqATM,NDEG)
C
C  find symmetry-unique atoms
C
      CALL GetUNQ(MAtom,  NTrans, XC,     NEqATM, NQ,
     $            IUNQ,   Z)
C
C  display
C
      If(IPRNT.GT.0) Then
       WRITE(IOut,1200)
       CALL PrntCAR(IOut,0,NAtoms,Symb,XC)
       group1=group
       call uppercase(group1,1)
       WRITE(IOut,1300) group1,NDEG,NQ
c -- check for, and alert user to, the presence of ghost atoms
       CALL ChkGhost(IOut,MAtom,Symb,XC,XCharg,IPRNT)
       If(IPRNT.GT.2)
     $    CALL PrntGEOM(IOut,NAtom,Symb,XC,Z(1),Z(1+NAtoms))
       If(RThrsh.GT.Zero) CALL PrntDIST(IOut,NAtoms,Symb,XC,RThrsh)
c -- write to summary output
       If(IPRNT.GT.1) Then
        WRITE(ICond,1200)
        CALL PrntCAR(ICond,0,NATOMS,Symb,XC)
        group1=group
        call uppercase(group1,1)
        WRITE(ICond,1300) group1,NDEG,NQ
        call f_lush(icond)
       EndIf
      EndIf
C
C  now put symmetry information onto <sym> file
C
      CALL WrSYM(MAtom,  RM,     GROUP,  NTrans, NDEG,
     $           NQ,     IUNQ,   TRANS,  NEqATM)
C
C  write new coordinates to <coord> file
C
      CALL WrCoord(NAtoms, Symb,   XC,     NMol,   IMOL,
     $             XCharg, XMass)
C
C  -- NEW --
C  rotate existing Hessian to new coordinate frame
C  ** WARNING ** Don't rotate if Z-matrix is being used
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.zmat',
     $      FORM='FORMATTED',STATUS='OLD',ERR=94)
      CLOSE (UNIT=IUNit,STATUS='KEEP')
      RETURN
 94   CONTINUE
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.hess',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      CALL ChkROT(RM,IPRNT,thrsh,rotate)
      If(rotate) Then
       IH1 = 1
       IH2 = IH1 + 9*MAtom**2
       CALL RdHESS(jobname(1:lenJ)//'.hess',lenJ+5,3*MAtom,IPRNT,
     $             Z(IH1),IHess)
       CALL RotHES(MAtom,RM,Z(IH2),Z(IH1))
       CALL WrHESS(jobname(1:lenJ)//'.hess',lenJ+5,0,3*MAtom,Z(IH1))
       WRITE(IOut,1400)
      EndIf
C
 95   CONTINUE
      RETURN
c
  900 Format(9X,I4,2X,I4)
 1000 FORMAT(' WARNING - coordinate axes will NOT be reoriented')
 1100 FORMAT(' WARNING - system will NOT use "centre-of-mass" frame')
 1200 FORMAT(/,4X,' Cartesian Coordinates in Standard Orientation')
 1300 FORMAT('   Point Group: ',A4,'  Number of degrees of',
     $         ' freedom: ',I3,/,
     $       '   Number of symmetry-unique atoms: ',I3,/)
 1400 FORMAT(' WARNING - Existing Hessian rotated to new orientation')
c
      END
c ========================================================================
c
      SUBROUTINE SYMTRY(NATOMS, AtSymb, XOld,   thrsym0,adjst,
     $                  IPRNT,  IMem,   Z,      XNew,   RM,
     $                  GROUP,  NTrans, TRANS,  NEqATM, NDEG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines Molecular Point Group Symmetry and number
C  of degrees of freedom
C
C  ARGUMENTS
C
C  NATOMS  -  total number of atoms
C  AtSymb  -  atomic symbols
C  XOld    -  current cartesian coordinates
C             (put into cms frame on exit)
C  thrsym0 -  accuracy threshold for symmetry determination
C  adjst   -  logical flag for reorienting axes
C  IPRNT   -  flag for controlling printout (mainly debug)
C  IMem    -  amount of scratch memory
C  Z       -  scratch storage array
C
C  on exit ...
C
C  XNew    -  reoriented coordinates
C  RM      -  rotation matrix effecting reorientation
C              CNew = RM * COld
C  GROUP   -  molecular point group
C  NTrans  -  number of symmetry operations found
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NDEG    -  number of internal degrees of freedom
C
C
      REAL*8 XOld(3,NATOMS),XNew(3,NATOMS),RM(3,3)
      DIMENSION TRANS(3,3,*),NEqATM(NATOMS,*)
      CHARACTER*8 AtSymb(NATOMS)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 Symb(NAtoms)
c ..................................................
      CHARACTER*4 GROUP
      LOGICAL adjst
C
      DIMENSION Z(IMem)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0,MSymOP=120,small=1.0d-8)
      PARAMETER (HALF=0.5d0,RV1=0.86602540378444d0)
      PARAMETER (GR1=0.80901699437495d0,GR2=GR1-HALF)  ! for Ih
C
C
C  Put into cms frame
C
      CALL CMS(NATOMS,CX,CY,CZ,XOld)
C
C  assign scratch storage
C
      ISH = 1
      IR = ISH + 3*NATOMS
      IRS = IR + NATOMS
      ILB = IRS + NATOMS
      INP = ILB + NATOMS
      IUN = INP + NATOMS
      IST = IUN + 3*MSymOP
      IIN = IST + MSymOP
      IIP = IIN + MSymOP
      IJP = IIP + MSymOP
      ILG = IJP + MSymOP
      IAT = ILG + NATOMS        ! this now allocated to avoid PC error
      IEnd = IAT + NATOMS - 1
      CALL MemCHK(IMem,IEnd,6,'SYMTRY')
c
      thrsym = thrsym0
 50   CONTINUE
      CALL CpyVEC(3*NATOMS,XOld,XNew)
cc
      CALL GetSYM(XNew,   Z(ISH), Z(IR),  Z(IRS), TRANS,
     $            Z(IUN), Z(ILB), Z(INP), Z(IST), Z(IIN),
     $            Z(IIP), Z(IJP), NATOMS, NTrans, Z(ILG),
     $            AtSymb, Symb,   GROUP,  MSymOP, thrsym,
     $            adjst,  IErr)
cc
      IF(GROUP.NE.'c1  '.AND.adjst) THEN
       adjst = .FALSE.
       CALL GetSYM(XNew,   Z(ISH), Z(IR),  Z(IRS), TRANS,
     $             Z(IUN), Z(ILB), Z(INP), Z(IST), Z(IIN),
     $             Z(IIP), Z(IJP), NATOMS, NTrans, Z(ILG),
     $             AtSymb, Symb,   GROUP,  MSymOP, thrsym,
     $             adjst,  IErr)
      ELSE
C
C  Put molecule in "Standard Orientation"
C  ** WARNING **  This needs to be done for non-C1 molecules too **
C
cc       CALL ORIENT(NATOMS,XNew,Z(ISH),Z(IR),Z(IST),RM)
      ENDIF
c
      IF(IErr.NE.0) THEN
       WRITE(6,1000)
       GROUP = 'c1  '
       NTrans = 1
      ENDIF
C
C  Make sure symmetry operations are all OK
C  (assuming "standard" orientation)
C
      DO 10 I=1,3
      DO 10 J=1,3
      DO 10 K=1,NTrans
      IF(Abs(TRANS(I,J,K)).LT.thrsym) THEN
        TRANS(I,J,K) = ZERO
      ELSE IF(Abs(Abs(TRANS(I,J,K))-ONE).LT.thrsym) THEN
        TRANS(I,J,K) = SIGN(ONE,TRANS(I,J,K))
      ELSE IF(Abs(Abs(TRANS(I,J,K))-HALF).LT.thrsym) THEN
        TRANS(I,J,K) = SIGN(HALF,TRANS(I,J,K))
      ELSE IF(Abs(Abs(TRANS(I,J,K))-RV1).LT.thrsym) THEN
        TRANS(I,J,K) = SIGN(RV1,TRANS(I,J,K))
      ELSE IF(Abs(Abs(TRANS(I,J,K))-GR1).LT.thrsym) THEN
        TRANS(I,J,K) = SIGN(GR1,TRANS(I,J,K))
      ELSE IF(Abs(Abs(TRANS(I,J,K))-GR2).LT.thrsym) THEN
        TRANS(I,J,K) = SIGN(GR2,TRANS(I,J,K))
      ENDIF
 10   CONTINUE
C
C  determine the equivalent atoms array
C
      CALL EQVATM(NATOMS, XNew,   NTrans, TRANS,  thrsym,
     $            NEqATM, IErr)
c
      IF(IErr.NE.0) THEN
       If(thrsym.GT.1.0d-5) Then
c -- if equivalent atoms cannot be found, try reducing threshold
        thrsym = thrsym*0.1d0
        WRITE(6,1050) thrsym
cc        adjst = .True.
        GO TO 50
       Else
        WRITE(6,1100)
        GROUP = 'c1  '
        NTrans = 1
        CALL CpyVEC(3*NATOMS,XOld,XNew)
       EndIf
      ENDIF
c
      If(IPRNT.GT.3) CALL PrntSYM(NATOMS,GROUP,NTrans,TRANS,NEqATM)
C
C  symmetrize the coordinates
C
      CALL SymVEC(NATOMS, NTrans, NEqATM, TRANS,  Z(IRS),
     $            small,  XNew)
C
C  find the rotation matrix connecting the original and
C  reoriented coordinates
C
      If(thrsym.LE.1.0d-5) Then
        CALL GetROT(NATOMS, IPRNT,  XOld,  XNew,  thrsym,
     $              RM,     IErr)
        If(IErr.NE.0) Call nerror(3,'GEOMETRY module',
     $   'Problems with Rotation Matrix in <GetROT>',0,0)
      EndIf
C
      RETURN
c
 1000 FORMAT('**WARNING** Error determining Point Group - Switching',
     $       ' off Symmetry')
 1050 FORMAT(' Equivalent atoms could not be determined.',
     $       ' Lowering Symmetry threshold to: ',F10.6)
 1100 FORMAT('**WARNING** Error determining Equivalent Atoms -',
     $       ' Switching off Symmetry')
c
      END
c ========================================================================
c
      SUBROUTINE EQVATM(NATOMS, XC,     NTrans, TRANS,  thrsym,
     $                  NEqATM, IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine atomic equivalences under the various
C  symmetry operations
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  XC      -  atomic coordinates
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  thrsym  -  threshold for determining equivalences
C             a given symmetry operation must move atom I
C             to within thrsym of atom J for atoms I and J
C             to be equivalent under this operation
C  NEqATM  -  on exit contains atom equivalence array
C  IErr    -  error flag   0 - success
C                         -1 - problems with equivalences
C
C
      REAL*8 XC(3,NATOMS),V(3)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NATOMS,NTrans)
      LOGICAL CmpVEC,same
C
C
      IErr = -1
C
C  first operation is always the identity
C
      DO 10 IAtm=1,NATOMS
      NEqATM(IAtm,1) = IAtm
 10   CONTINUE
C
      If(NTrans.EQ.1) GO TO 95
C
C  now loop over the other operations
C
      DO 40 IAtm=1,NATOMS
      DO 30 IOP=2,NTrans
      NEqATM(IAtm,IOP) = 0
      CALL MatVEC(3,TRANS(1,1,IOP),XC(1,IAtm),V)
      DO 20 JAtm=1,NATOMS
      same = CmpVEC(3,V,XC(1,JAtm),thrsym)
      If(same) Then
       NEqATM(IAtm,IOP) = JAtm
       GO TO 30
      ENDIF
 20   CONTINUE
      RETURN                     ! error
 30   CONTINUE
 40   CONTINUE
C
 95   CONTINUE
      IErr = 0
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE GetROT(NATOMS, IPRNT,  XOld,   XNew,   thrsh,
     $                  RM,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine the rotation matrix linking two sets of
C  Cartesian coordinates
C              XNew = RM * XOld
C  IMPORTANT:  It is assumed that both coordinate sets
C  are in the centre of mass frame (with all atoms
C  having unit mass).
C  On exit, RM contains the rotation matrix
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  IPRNT   -  controls level of print out
C  XOld    -  old coordinates (before reorientation)
C  XNew    -  new coordinates (after reorientation)
C  thrsh   -  accuracy threshold
C  RM      -  on exit contains rotation matrix
C  IErr    -  error flag
C               0 - rotation matrix found
C              -1 - unable to determine rotation matrix
C
C
      REAL*8 XOld(3,NATOMS),XNew(3,NATOMS),RM(3,3)
      REAL*8 C1(3,3),C2(3,3)
      LOGICAL same,CmpVEC
C
      PARAMETER (ZERO=0.0d0)
C
C  initialize error flag
C
      IErr = -1
C
C  get old orientation matrix (by schmidt orthogonalization)
C
      CALL GRAM(IL,NATOMS,XOld,C1)
      IF(IL.EQ.4) THEN
       WRITE(6,1000)
       RETURN
      ENDIF
c
      IF(IPRNT.GT.5) THEN
       WRITE(6,1100)
       Call PrntMAT(3,3,3,C1)
      ENDIF
C
C  get new orientation matrix
C
      CALL GRAM(IL,NATOMS,XNew,C2)
      IF(IL.EQ.4) THEN
       WRITE(6,1200)
       RETURN
      ENDIF
c
      IF(IPRNT.GT.5) THEN
       WRITE(6,1300)
       Call PrntMAT(3,3,3,C2)
      ENDIF
C
C  Get rotation matrix RM
C    CNew = RM*COld  ==>  RM = C2*C1(t)
C
      DO 20 I=1,3
      DO 20 J=1,3
      TEMP = ZERO
      DO 10 K=1,3
      TEMP = TEMP + C2(I,K)*C1(J,K)
 10   CONTINUE
      RM(I,J) = TEMP
 20   CONTINUE
C
C  Check accuracy of transformation
C
      DO 30 I=1,NATOMS
      C1(1,1) = RM(1,1)*XOld(1,I) + RM(1,2)*XOld(2,I)
     $             + RM(1,3)*XOld(3,I)
      C1(2,1) = RM(2,1)*XOld(1,I) + RM(2,2)*XOld(2,I)
     $             + RM(2,3)*XOld(3,I)
      C1(3,1) = RM(3,1)*XOld(1,I) + RM(3,2)*XOld(2,I)
     $             + RM(3,3)*XOld(3,I)
c
      same = CmpVEC(3,C1,XNew(1,I),thrsh)
      IF(.NOT.same) THEN
       WRITE(6,1400) thrsh
       write(6,*) ' fails for atom: ',i
       write(6,*) ' old coordinates'
       call prntmat(natoms,3,natoms,xold)
       write(6,*) ' new coordinates'
       call prntmat(natoms,3,natoms,xnew)
       write(6,*) ' rotation matrix'
       call prntmat(3,3,3,rm)
       write(6,*) ' C1 vector'
       call prntmat(1,3,1,c1)
       RETURN
      ENDIF
c
 30   CONTINUE
C
C  Zero out small elements
C
      DO 40 I=1,3
      DO 40 J=1,3
      If(Abs(RM(I,J)).LT.thrsh) RM(I,J) = ZERO
 40   CONTINUE
c
      IF(IPRNT.GT.5) THEN
       WRITE(6,1500)
       Call PrntMAT(3,3,3,RM)
      ENDIF
C
      IErr = 0
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unable to construct Orientation',
     $            ' Matrix for Old Coordinates in <GetROT>',/)
 1100 FORMAT(/,' Orientation Matrix for old Coordinates is:')
 1200 FORMAT(/,2X,'***ERROR*** Unable to construct Orientation',
     $            ' Matrix for New Coordinates in <GetROT>',/)
 1300 FORMAT(/,' Orientation Matrix for new Coordinates is:')
 1400 FORMAT(/,2X,'***ERROR*** Coordinates do not transform',
     $            ' within specified threshold of',D10.5,
     $ '     Either increase symmetry threshold or turn off symmetry')
 1500 FORMAT(/,' Rotation Matrix is:')
c
      END
c ========================================================================
c
      SUBROUTINE GetUNQ(NATOMS, NTrans, XC,     NEqATM, NQ,
     $                  IUNQ,   ISYM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines which atoms in the current system are
C  symmetry unique
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  NTrans  -  number of symmetry operations
C  XC      -  Cartesian coordinates of all atoms
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C  on exit
C
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of symmetry unique atoms in ascending order
C  ISYM    -  list of all atoms
C              ISYM(I) = 1 if atom I symmetry-unique
C              ISYM(I) = 0 otherwise
C
C
      REAL*8 XC(3,NATOMS)
      INTEGER NEqATM(NATOMS,NTrans),IUNQ(NATOMS),ISYM(NATOMS)
C
      PARAMETER (ZERO=0.0d0)
C
      CALL IZeroIT(ISYM,NATOMS)
      NQ = 0
c
      DO 20 I=1,NATOMS
c
      If(ISYM(I).EQ.1) GO TO 20
      ISYM(I) = 1
C
C  atom I is potentially symmetry unique
C  check number of zero coordinates for this atom
C
      NQ = NQ+1
      IUNQ(NQ) = I
c
      NZero = 0
      If(XC(1,I).EQ.ZERO) NZero = NZero+1
      If(XC(2,I).EQ.ZERO) NZero = NZero+1
      If(XC(3,I).EQ.ZERO) NZero = NZero+1
      MZero = NZero
C
C  check all symmetry operations
C  locate all operations which convert atom I to a new atom
C  If this atom has more zero coordinates than I then
C  take this atom as symmetry unique
C
      DO 10 IOP=1,NTrans
      J = NEqATM(I,IOP)
      IF(ISYM(J).NE.1) THEN
       ISYM(J) = 1
       NZero = 0
       If(XC(1,J).EQ.ZERO) NZero = NZero+1
       If(XC(2,J).EQ.ZERO) NZero = NZero+1
       If(XC(3,J).EQ.ZERO) NZero = NZero+1
       If(NZero.GT.MZero) IUNQ(NQ) = J
      ENDIF
 10   CONTINUE
c
 20   CONTINUE
C
C  Now sort symmetry-unique atoms into ascending order
C  and fill ISYM array
C
      CALL IZeroIT(ISYM,NATOMS)
      DO 30 I=1,NQ
      IQ = IUNQ(I)
      ISYM(IQ) = 1
 30   CONTINUE
c
      IQ = 0
      DO 40 I=1,NATOMS
      IF(ISYM(I).EQ.1) THEN
       IQ = IQ+1
       IUNQ(IQ) = I
      ENDIF
 40   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE NumDEG(NATOMS,NTrans,TRANS,GROUP,NEqATM,NDEG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine determines the number of symmetry-independent
C  degrees of freedom
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 permutation matrices
C  GROUP   -  molecular point group
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NDEG    -  on exit number of degrees of freedom
C
C
      DIMENSION TRANS(9,NTrans),NEqATM(NATOMS,NTrans)
      CHARACTER*4 GROUP
C
      PARAMETER (ZERO=0.0d0)
C
C
      NDEG = 0
      If(NATOMS.EQ.1) RETURN
C
C  special cases for linear molecules
C
      IF(GROUP.EQ.'c*v ') THEN
       NDEG = NATOMS-1
       RETURN
      ELSE IF(GROUP.EQ.'d*h ') THEN
       NDEG = NATOMS/2
       RETURN
      ENDIF
C
      XYZSUM = ZERO
      ROTSUM = ZERO
      CRTSUM = ZERO
c
      DO 20 IOP=1,NTrans
      XX = TRANS(1,IOP)
      XY = TRANS(2,IOP)
      XZ = TRANS(3,IOP)
      YX = TRANS(4,IOP)
      YY = TRANS(5,IOP)
      YZ = TRANS(6,IOP)
      ZX = TRANS(7,IOP)
      ZY = TRANS(8,IOP)
      ZZ = TRANS(9,IOP)
      XYZ = XX+YY+ZZ
      XYZSUM = XYZSUM + XYZ
c
      DO 10 IATM=1,NATOMS
      If(NEqATM(IATM,IOP).EQ.IATM) CRTSUM = CRTSUM + XYZ
 10   CONTINUE
c
      XXYYZZ = XX*(YY*ZZ - YZ*ZY) - XY*(YX*ZZ - YZ*ZX)
     $            + XZ*(YX*ZY - YY*ZX)
      If(XXYYZZ.GT.ZERO) XYZ = -XYZ
      ROTSUM = ROTSUM - XYZ
 20   CONTINUE
C
      NT = NINT(XYZSUM/NTrans)
      NR = NINT(ROTSUM/NTrans)
      NC = NINT(CRTSUM/NTrans)
C
      NDEG = NC - NT - NR
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE PrntSYM(NATOMS,GROUP,NTrans,TRANS,NEqATM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out symmetry information
C
      INTEGER NEqATM(NATOMS,NTrans)
      REAL*8 TRANS(3,3,NTrans)
      CHARACTER*4 GROUP,groupc,group1
C
      groupc=group
      if(groupc(1:1).eq.'c') groupc(1:1)='C'
      if(groupc(1:1).eq.'d') groupc(1:1)='D'
      if(groupc(1:1).eq.'s') groupc(1:1)='S'
      if(groupc(1:1).eq.'t') groupc(1:1)='T'
      if(groupc(1:1).eq.'o') groupc(1:1)='O'
      if(groupc(1:1).eq.'i') groupc(1:1)='I'
C
      iout=igetival('iout')
      group1=groupc
      call uppercase(group1,1)
      WRITE(iout,1000) group1,NTrans
 1000 FORMAT(/,' Point group: ',A4,'  Number of symmetry',
     $         ' operations: ',I3)
c
      DO 10 I=1,NTrans
        WRITE(iout,1100) I
 1100   FORMAT(/,' Symmetry operation: ',I3)
        CALL PrntMat(3,3,3,TRANS(1,1,I))
 10   CONTINUE
c
      WRITE(iout,1200)
 1200 FORMAT(/,' Equivalent atoms array')
      DO 20 I=1,NATOMS
        WRITE(iout,1300) (NEqATM(I,J),J=1,NTrans)
 1300   FORMAT(20I4)
 20   CONTINUE
C
      END
c ======================================================================
c
      SUBROUTINE SymCHK(IOut,NAtoms,NAtom,NSym,NEqATM,XCharg)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Checks if atomic charges are consistent with symmetry
C  This will pick up ghost atoms and charged dummy atoms
C  that may break symmetry
C  (quick and dirty and not the best way to do this - should include
C   charge as a part of the symmetry analysis)
C
C  ARGUMENTS
C
C  IOut    -  output unit for error printout
C  NAtoms  -  number of real + charged dummy atoms
C  NAtom   -  number of atoms to check
C  NSym    -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  XCharg  -  charges on atomic centres
C
C
      DIMENSION NEqATM(NAtoms,NSym),XCharg(NAtom)
      Logical Error
C
      Data Error/.False./
C
C
      DO 20 IAtm=1,NAtom
      Charge = XCharg(IAtm)
C
C  check all symmetry-equivalent atoms
C
      DO 10 IOP=2,NSym
      JAtm = NEqATM(IAtm,IOP)
      If(XCharg(JAtm).NE.Charge) Then
       WRITE(IOut,1000) IAtm,JAtm
       Error = .True.
      EndIf
 10   CONTINUE
 20   CONTINUE
c
        If(Error) Call nerror(8,'GEOMETRY MODULE',
     $     'Atomic Charges Break Symmetry',0,0)
C
      RETURN
c

 1000 FORMAT(' Atoms ',2I4,' Should have the same charge for',
     $       ' this symmetry'/,' Either make charges equal, change',
     $       ' atomic symbol or switch off symmetry')
c
      END
c ========================================================================
c
      SUBROUTINE SymEqAtm(NAtoms,NSym,NEqATM,IVec,ISYM)
      IMPLICIT INTEGER(A-Z)
C
C  determines how many symmetry partners each atom has
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NSym    -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  IVec    -  work array
C  ISYM    -  on exit number of symmetry partners per atom
C
C
      DIMENSION NEqATM(NAtoms,NSym),IVec(NAtoms),ISYM(NAtoms)
C
C  loop over atoms
C
      DO 20 IAtm=1,NAtoms
      Call IZeroIT(IVec,NAtoms)
      IVec(IAtm) = 1            ! for identity (which may not be included)
C
C  loop over symmetry operations
C
      DO 10 IOP=1,NSym
      II = NEqATM(IAtm,IOP)
      IVec(II) = 1
 10   CONTINUE
C
C  now add up all IVec entries to determine number of symmetry partners
C  (1 means atom is symmetry unique)
C
      II = 0
      Do I=1,NAtoms
      II = II + IVec(I)
      EndDo
c
      ISYM(IAtm) = II
 20   CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE SymFILTER(NATOMS,NAtom,NTrans,NEqOLD,NEqNEW)
      IMPLICIT INTEGER(A-Z)
C
C
C  Filters Equivalent atoms array to exclude dummy atoms
C
C  ARGUMENTS
C
C  NATOMS  -  number of real atoms + charged dummy atoms
C  NAtom   -  number of real atoms only
C  NTrans  -  number of symmetry operations
C  NEqOLD  -  equivalent atoms array including dummies
C  NEqNEW  -  on exit equivalent atoms array for real atoms only
C
C
      DIMENSION NEqOLD(NATOMS,NTrans),NEqNEW(NAtom,NTrans)
C
      DO 10 IOP=1,NTrans
      DO 10 IAtm=1,NAtom
      NEqNEW(IAtm,IOP) = NEqOLD(IAtm,IOP)
 10   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE PrntDIST(IOut,NAtoms,AtSymb,XC,RThrsh)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  print out all interatomic distances below threshold
C
C  ARGUMENTS
C
C  IOut    -  output unit
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  RThrsh  -  threshold for distance print out
C
C
      DIMENSION XC(3,NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
      Character*80 Char
      Logical special
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      AUTOAN = 1.0d0/ANTOAU
C
C
C  First determine if there are any user-defined atomic labels
C  If so, then DO NOT number atoms; otherwise atoms are numbered
C  for easier identification during print out
C
C  Initialize position of zero in collating series
C
      JZero = ICHAR('0')
      special = .True.
c
      DO 5 I=1,NATOMS
      DO 4 j=2,8
      If(AtSymb(I)(j:j).EQ.' ') exit
      ICh = ICHAR(AtSymb(I)(j:j)) - JZero
      If(ICh.GE.0.AND.ICh.LE.9) GO TO 6     ! found a number
 4    CONTINUE
 5    CONTINUE
C
C  if reached here, there are NO numbers in ANY atomic symbols
C
      special = .False.
c
 6    CONTINUE
      WRITE(IOut,1000) RThrsh
c
      DO 10 I=2,NATOMS
      DO 10 J=1,I-1
      XIJ = XC(1,I) - XC(1,J)
      YIJ = XC(2,I) - XC(2,J)
      ZIJ = XC(3,I) - XC(3,J)
      R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
      R = R*AUTOAN
      IF(R.LT.RThrsh) THEN
       If(special) Then
        WRITE(Char,1140) AtSymb(I),'-',AtSymb(J)
       Else
        WRITE(Char,1150) AtSymb(I),I,'-',AtSymb(J),J
       EndIf
       Call RmBlank(80,Char,Len)
       WRITE(IOut,1200) Char(1:10),R
      ENDIF
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,' InterAtomic Distances Below ',F9.6,' Angstroms')
 1140 FORMAT(4(A8,A1))
 1150 FORMAT(4(A8,I4,A1))
 1200 FORMAT(2X,A10,F10.6)
c
      END
c ========================================================================
c
      SUBROUTINE PrntGEOM(IOut,NAtoms,AtSymb,XC,IAN,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  print out molecular geometry in terms of internal coordinates
C  (bond lengths, bond angles and torsions)
C
C  ARGUMENTS
C
C  IOut    -  output unit
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  IAN     -  atomic numbers
C  IC      -  connectivity matrix
C
C
      DIMENSION XC(3,NAtoms),IAN(NAtoms),IC(NAtoms,NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
C
      PARAMETER (thrbnd=0.85d0)       ! distance ratio for bonding
C
C
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NATOMS,AtSymb,IAN)
C
C  get connectivity matrix
C
      CALL IZeroIT(IC,NAtoms*NAtoms)
      CALL CONNECTM(NAtoms,IAN,XC,thrbnd,IC,IErr)
C
C  print out geometries
C
      CALL PrntINTL(IOut,NAtoms,AtSymb,XC,IC)
C
      RETURN
      END
c =========================================================================
c
      SUBROUTINE PrntINTL(IOut,NAtoms,AtSymb,XC,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates bond stretches and primitive bends and torsions
C  from the atomic connectivity matrix
C
C  ARGUMENTS
C
C  IOut    -  output unit
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  IC      -  atomic connectivity matrix
C
      DIMENSION XC(3,NAtoms),IC(NAtoms,NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
      LOGICAL Bend,Outp,Tors,special
      Character*80 Char
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      AUTOAN = 1.0d0/ANTOAU
      TOANG = 180.0d0/PI
C
      DATA Bend/.TRUE./, Outp/.FALSE./, Tors/.TRUE./
C
C
C  First determine if there are any user-defined atomic labels
C  If so, then DO NOT number atoms; otherwise atoms are numbered
C  for easier identification during print out
C
C  Initialize position of zero in collating series
C
      JZero = ICHAR('0')
      special = .True.
c
      DO 5 I=1,NATOMS
      DO 4 j=2,8
      If(AtSymb(I)(j:j).EQ.' ') exit
      ICh = ICHAR(AtSymb(I)(j:j)) - JZero
      If(ICh.GE.0.AND.ICh.LE.9) GO TO 6     ! found a number
 4    CONTINUE
 5    CONTINUE
C
C  if reached here, there are NO numbers in ANY atomic symbols
C
      special = .False.
c
 6    CONTINUE
      WRITE(IOut,1000)
C
C  (a) Stretches
C
      WRITE(IOut,1100)
c
      DO 10 I=2,NATOMS
      DO 10 J=1,I-1
      IF(IC(I,J).NE.0) THEN
       XIJ = XC(1,I) - XC(1,J)
       YIJ = XC(2,I) - XC(2,J)
       ZIJ = XC(3,I) - XC(3,J)
       R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
       R = R*AUTOAN
       If(special) Then
        WRITE(Char,1140) AtSymb(I),'-',AtSymb(J)
       Else
        WRITE(Char,1150) AtSymb(I),I,'-',AtSymb(J),J
       EndIf
       Call RmBlank(80,Char,Len)
       WRITE(IOut,1200) Char(1:10),R
      ENDIF
 10   CONTINUE
c
      If(.NOT.Bend.OR.NATOMS.LT.3) GO TO 25
C
C  (b) Bends
C
      WRITE(IOut,1300)
c
      DO 20 I=3,NATOMS
      DO 20 J=2,I-1
      DO 20 K=1,J-1
c
      IF(IC(I,J).NE.0.AND.IC(J,K).NE.0) THEN
       CALL AngGRAD(NAtoms,I,J,K,XC,Ang,.false.,jnk)
       Ang = Ang*TOANG
       If(special) Then
        WRITE(Char,1140) AtSymb(I),'-',AtSymb(J),'-',AtSymb(K)
       Else
        WRITE(Char,1150) AtSymb(I),I,'-',AtSymb(J),J,'-',AtSymb(K),K
       EndIf
       Call RmBlank(80,Char,Len)
       WRITE(IOut,1210) Char(1:15),Ang
      ENDIF
c
      IF(IC(J,K).NE.0.AND.IC(K,I).NE.0) THEN
       CALL AngGRAD(NAtoms,J,K,I,XC,Ang,.false.,jnk)
       Ang = Ang*TOANG
       If(special) Then
        WRITE(Char,1140) AtSymb(J),'-',AtSymb(K),'-',AtSymb(I)
       Else
        WRITE(Char,1150) AtSymb(J),J,'-',AtSymb(K),K,'-',AtSymb(I),I
       EndIf
       Call RmBlank(80,Char,Len)
       WRITE(IOut,1210) Char(1:15),Ang
      ENDIF
c
      IF(IC(K,I).NE.0.AND.IC(I,J).NE.0) THEN
       CALL AngGRAD(NAtoms,K,I,J,XC,Ang,.false.,jnk)
       Ang = Ang*TOANG
       If(special) Then
        WRITE(Char,1140) AtSymb(K),'-',AtSymb(I),'-',AtSymb(J)
       Else
        WRITE(Char,1150) AtSymb(K),K,'-',AtSymb(I),I,'-',AtSymb(J),J
       EndIf
       Call RmBlank(80,Char,Len)
       WRITE(IOut,1210) Char(1:15),Ang
      ENDIF
c
 20   CONTINUE
c
 25   CONTINUE
      If(.NOT.Outp.OR.NATOMS.LT.4) GO TO 35
C
C  (c) Out-of-Plane Bends
C
      WRITE(IOut,1400)
c
cc      DO 30 I=4,NATOMS
cc      DO 30 J=3,I-1
cc      DO 30 K=2,J-1
cc      DO 30 L=1,K-1
c
cc      IF(IC(I,J).NE.0.AND.IC(I,K).NE.0.AND.IC(I,L).NE.0) THEN
cc       CALL OutpGRAD(NAtoms,J,K,L,I,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(J),J,'-',AtSymb(K),K,'-',AtSymb(L),L,'-',
cc     $                  AtSymb(I),I
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,K,L,J,I,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(K),K,'-',AtSymb(L),L,'-',AtSymb(J),J,'-',
cc     $                  AtSymb(I),I
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,L,J,K,I,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(L),L,'-',AtSymb(J),J,'-',AtSymb(K),K,'-',
cc     $                  AtSymb(I),I
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
cc      ENDIF
c
cc      IF(IC(J,I).NE.0.AND.IC(J,K).NE.0.AND.IC(J,L).NE.0) THEN
cc       CALL OutpGRAD(NAtoms,I,K,L,J,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(I),I,'-',AtSymb(K),K,'-',AtSymb(L),L,'-',
cc     $                  AtSymb(J),J
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,K,L,I,J,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(K),K,'-',AtSymb(L),L,'-',AtSymb(I),I,'-',
cc     $                  AtSymb(J),J
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,L,I,K,J,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(L),L,'-',AtSymb(I),I,'-',AtSymb(K),K,'-',
cc     $                  AtSymb(J),J
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
cc      ENDIF
c
cc      IF(IC(K,I).NE.0.AND.IC(K,J).NE.0.AND.IC(K,L).NE.0) THEN
cc       CALL OutpGRAD(NAtoms,I,J,L,K,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(I),I,'-',AtSymb(J),J,'-',AtSymb(L),L,'-',
cc     $                  AtSymb(K),K
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,J,L,I,K,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(J),J,'-',AtSymb(L),L,'-',AtSymb(I),I,'-',
cc     $                  AtSymb(K),K
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,L,I,J,K,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(L),L,'-',AtSymb(I),I,'-',AtSymb(J),J,'-',
cc     $                  AtSymb(K),K
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
cc      ENDIF
c
cc      IF(IC(L,I).NE.0.AND.IC(L,J).NE.0.AND.IC(L,K).NE.0) THEN
cc       CALL OutpGRAD(NAtoms,I,J,K,L,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(I),I,'-',AtSymb(J),J,'-',AtSymb(K),K,'-',
cc     $                  AtSymb(L),L
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,J,K,I,L,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(J),J,'-',AtSymb(K),K,'-',AtSymb(I),I,'-',
cc     $                  AtSymb(L),L
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
c
cc       CALL OutpGRAD(NAtoms,K,I,J,L,XC,Ang,.false.,jnk)
cc       Ang = Ang*TOANG
cc       WRITE(Char,1150) AtSymb(K),K,'-',AtSymb(I),I,'-',AtSymb(J),J,'-',
cc     $                  AtSymb(L),L
cc       Call RmBlank(80,Char,Len)
cc       WRITE(IOut,1220) Char(1:20),Ang
cc      ENDIF
c
cc 30   CONTINUE
c
c .........................................................
c  ** WARNING   PP Version    Has this been tested?
      DO I=1,NATOMS
        do j=1,natoms
          if(ic(j,i).ne.0) then
            do k=j+1,natoms
              if(ic(k,i).ne.0) then
                do l=k+1,natoms
                  if(ic(l,i).ne.0) then
                    CALL OutpGRAD(NAtoms,J,K,L,I,XC,Ang,.false.,jnk)
                    Ang = Ang*TOANG
       WRITE(Char,1150) AtSymb(J),J,'-',AtSymb(K),K,'-',AtSymb(L),L,'-',
     $                  AtSymb(I),I
                    Call RmBlank(80,Char,Len)
                    WRITE(IOut,1220) Char(1:20),Ang
c
                    CALL OutpGRAD(NAtoms,K,L,J,I,XC,Ang,.false.,jnk)
                    Ang = Ang*TOANG
       WRITE(Char,1150) AtSymb(K),K,'-',AtSymb(L),L,'-',AtSymb(J),J,'-',
     $                  AtSymb(I),I
                    Call RmBlank(80,Char,Len)
                    WRITE(IOut,1220) Char(1:20),Ang
c
                    CALL OutpGRAD(NAtoms,L,J,K,I,XC,Ang,.false.,jnk)
                    Ang = Ang*TOANG
       WRITE(Char,1150) AtSymb(L),L,'-',AtSymb(J),J,'-',AtSymb(K),K,'-',
     $                  AtSymb(I),I
                    Call RmBlank(80,Char,Len)
                    WRITE(IOut,1220) Char(1:20),Ang
      If(.NOT.Tors.OR.NATOMS.LT.4) GO TO 45
                  end if
                end do
              end if
            end do
          end if
        end do
      end do
c ............................................................
 35   CONTINUE
      If(.NOT.Tors.OR.NATOMS.LT.4) GO TO 45
C
C  (d) Torsions
C
      WRITE(IOut,1500)
c
      DO 42 I=2,NAtoms
      DO 42 J=1,I-1
      IF(IC(I,J).NE.0) THEN
c -- consider I-J as middle 2 atoms in proper torsion
        DO 41 K=1,NAtoms
        IF(IC(I,K).NE.0.AND.K.NE.J) THEN
c -- have K-I-J
          DO 40 L=1,NAtoms
          IF(IC(J,L).NE.0.AND.I.NE.L.AND.K.NE.L) THEN
            CALL DihGRAD(NAtoms,K,I,J,L,XC,Ang,.false.,jnk)
            Ang = Ang*TOANG
            If(special) Then
             WRITE(Char,1140) AtSymb(K),'-',AtSymb(I),'-',
     $                        AtSymb(J),'-',AtSymb(L)
            Else
             WRITE(Char,1150) AtSymb(K),K,'-',AtSymb(I),I,'-',
     $                        AtSymb(J),J,'-',AtSymb(L),L
            EndIf
            Call RmBlank(80,Char,Len)
            WRITE(IOut,1220) Char(1:20),Ang
          ENDIF
 40       CONTINUE
        ENDIF
 41     CONTINUE
      ENDIF
 42   CONTINUE
c
 45   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,' ---- MOLECULAR GEOMETRY (Angstroms/Degrees) ----')
 1100 FORMAT(/,' ** Bond Lengths **')
 1140 FORMAT(4(A8,A1))
 1150 FORMAT(4(A8,I4,A1))
 1200 FORMAT(2X,A10,F10.6)
 1210 FORMAT(2X,A15,F10.4)
 1220 FORMAT(2X,A20,F10.4)
 1300 FORMAT(/,' ** Bond Angles **')
 1400 FORMAT(/,' ** Out of Plane Bends **')
 1500 FORMAT(/,' ** Dihedral Angles **')
c
      END
c ========================================================================
c
      SUBROUTINE SortGEOM(NAtoms, XNuc,   AtMass, XScr,   XVec,
     $                    INDX,   Ndum1,  Ndum2)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Moves ALL dummy atoms to end of Coordinate list
C  -----------------------------------------------------------------
C  There are two types of dummy atom: Those with point charges
C  (to mimic the effects of an applied field) and those with no
C  charge (used to calculate a molecular property, e.g., chemical
C  shift, at a given point). Charged dummies are included when
C  determining the molecular symmetry; uncharged dummies are not.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XNuc    -  atomic nuclei (old TEXAS, includes charge and symbol)
C  AtMass  -  atomic masses
C  XScr    -  scratch array for nuclei
C  XVec    -  scratch vector for atomic masses
C  INDX    -  index array of dummy atoms
C  Ndum1   -  on exit number of charged dummy atoms
C  Ndum2   -  on exit number of uncharged dummy atoms
C
C  ** NOTE:  This routine assumes that dummy atom symbols
C            start with x (but not xe), q or du
C
C
      DIMENSION XNuc(5,NAtoms),XScr(5,NAtoms),INDX(NAtoms),
     $          AtMass(NAtoms),XVec(NAtoms)
      Character*8 name
      Equivalence (xname,name)
C
      PARAMETER (Zero=0.0d0)
C
C
      Ndum = 0
      Ndum1 = 0
      Ndum2 = 0
c
      DO 10 I=1,NAtoms
      xname = XNuc(5,I)
      If( (name(1:1).EQ.'x'.AND.name(2:2).NE.'e') .OR.
     $    (name(1:1).EQ.'d'.AND.name(2:2).EQ.'u') .OR.
     $     name(1:1).EQ.'q') Then
C
C  we have a dummy atom.
C  charged or uncharged?
C
       Ndum = Ndum+1
       If(XNuc(1,I).EQ.Zero) Then
        Ndum2 = NDum2+1
        INDX(Ndum) = I
       Else
        Ndum1 = Ndum1 + 1
        INDX(Ndum) = -I
       EndIf
      EndIf
 10   CONTINUE
c
      IF(Ndum.GT.0) THEN
C
C  There are dummy atoms, so need to sort coordinates
C
        NS = NAtoms-Ndum+1
        NS1 = NAtoms-Ndum2+1
c
        Ndum = 1
        Ndum1 = 0
        Ndum2 = 0
c
        II = 0
        Call CpyVEC(5*NAtoms,XNuc,XScr)
        Call CpyVEC(NAtoms,AtMass,XVec)
c
        DO 20 I=1,NAtoms
        IF(I.EQ.INDX(Ndum)) THEN
          CALL CpyVEC(5,XScr(1,I),XNuc(1,NS1+Ndum2))
          AtMass(NS1+Ndum2) = XVec(I)
          Ndum = Ndum+1
          Ndum2 = Ndum2+1
        ELSE IF(I.EQ.-INDX(Ndum)) THEN
          CALL CpyVEC(5,XScr(1,I),XNuc(1,NS+Ndum1))
          AtMASS(NS+Ndum1) = XVec(I)
          Ndum = Ndum+1
          Ndum1 = Ndum1+1
        ELSE
          II = II+1
          CALL CpyVEC(5,XScr(1,I),XNuc(1,II))
          AtMass(II) = XVec(I)
        ENDIF
 20     CONTINUE
      ENDIF
C
      RETURN
      END
c ========================================================================
c
      subroutine postgeom(natoms,natom,trans,meqatm,neqatm)

      use memory

      implicit real*8(a-h,o-z)
c
c  detects abelian symmetry operations and stores this information in
c  common /symm/ and in depository for compatability with old texas code
c
c  ARGUMENTS
c
c  natoms  -  total number of atomic centres (real and dummy)
c  natom   -  number of real atoms only
c  tran    -  3x3 array for symmetry operation(s)
c  meqatm  -  list of atomic equivalences under abelian symmetry operations
c  neqatm  -  list of atomic equivalences under all symmetry operations
c
c  data set here is
c  ----------------
c  nsym    -  number of abelian symmetry operations
c  nsy     -  list of abelian operations
c               1 - reflection in YZ plane
c               2 - reflection in XZ plane
c               3 - reflection in XY plane
c               4 - C2 rotation about Z-axis
c               5 - C2 rotation about Y-axis
c               6 - C2 rotation about X-axis
c               7 - inversion through origin
c  nupair  -  symmetry-equivalent atoms array for abelian operations
c             (extracted from array for all operations)
c
c     common /big/bl(3000)
c     common /intbl/maxsh,inx(100)
      common /symm/nsym,nsy(7)
      dimension trans(3,3),meqatm(natoms,7),neqatm(natom,*)
      character jobname*256,cdum*20,group*4
      character*3 symop(7)
c
      parameter (one=1.0d0)
      parameter (IUnit=1)
      data symop/' x ',' y ',' z ','xy ','xz ','yz ','xyz'/
c
c
c -- deal with filename header
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
c -- make sure final coordinates are in xnuc (meqatm used as scratch)
      call getival('inuc',inuc)
      call rdgeom(natoms,ndum1,ndum2,bl(inuc),meqatm,
     $            -1,    idum, idum, idum)
c
c  open <sym> file and read in operations one at a time
      open (unit=IUnit,file=jobname(1:lenJ)//'.sym',
     $      form='formatted',status='old')
      call rdcntrl(IUnit,6,'$group',3,ndum,dum,cdum)
      group = cdum(1:4)                        ! point group symbol
      call rdcntrl(IUnit,7,'$ntrans',1,ntrans,dum,cdum)
c
c -- locate start of symmetry operations
      call fdcntrl(IUnit,18,'symmetry_operation',ierr)
      If(ierr.eq.1) Call nerror(4,'GEOMETRY module',
     $    'No Symmetry Operations found on <sym> file!',0,0)
c
      backspace IUnit
c
      call izeroit(nsy,7)
c
c -- loop over symmetry operations
      Do iop=1,ntrans
      read(IUnit,'(a20)') cdum
c
      do i=1,3
      read(IUnit,900) (trans(i,j),j=1,3)
      enddo
c
c -- now check if operation is abelian
      If(trans(1,1).eq.-one.and.trans(2,2).eq.one.and.
     $        trans(3,3).eq.one) Then
        nsy(1) = iop
      Else If(trans(1,1).eq.one.and.trans(2,2).eq.-one.and.
     $        trans(3,3).eq.one) Then
        nsy(2) = iop
      Else If(trans(1,1).eq.one.and.trans(2,2).eq.one.and.
     $        trans(3,3).eq.-one) Then
        nsy(3) = iop
      Else If(trans(1,1).eq.-one.and.trans(2,2).eq.-one.and.
     $        trans(3,3).eq.one) Then
        nsy(4) = iop
      Else If(trans(1,1).eq.-one.and.trans(2,2).eq.one.and.
     $        trans(3,3).eq.-one) Then
        nsy(5) = iop
      Else If(trans(1,1).eq.one.and.trans(2,2).eq.-one.and.
     $        trans(3,3).eq.-one) Then
        nsy(6) = iop
      Else If(trans(1,1).eq.-one.and.trans(2,2).eq.-one.and.
     $        trans(3,3).eq.-one) Then
        nsy(7) = iop
      EndIf
c
      EndDo
c
c -- read equivalent atoms array
      read(IUnit,'(a20)') cdum
      If(cdum.ne.'$equivalent_atoms_ar')
     $   Call nerror(5,'GEOMETRY module',
     $    'Equivalent atoms array missing from <sym> file!',0,0)
c
      do i=1,ntrans
      read(IUnit,910) (neqatm(j,i),j=1,natom)
      enddo
c
      close (unit=IUnit,status='keep')
c
      call izeroit(meqatm,7*natoms)
c
c -- put data in common /symm/
      nsym = 0
      do i=1,7
      iop = nsy(i)
      if(iop.gt.0) then
       nsy(i) = 0
       nsym = nsym+1
       nsy(nsym) = i
       do j=1,natom
       meqatm(j,nsym) = neqatm(j,iop)
       enddo
      endif
      enddo
c
c -- determine the Abelian subgroup of the point group
      call abelian(nsym,nsy,group)
      call abeliansymm(natoms,natom,nsym,nsy,meqatm,bl(inuc))
c
c -- print the Abelian symmetry operations
      iout = igetival('iout')
      if(nsym.gt.0) then
        write(iout,1000) (symop(nsy(k)),k=1,nsym)
      end if
c
c -- determine the number of group generators
      ngener=nsym
      if(nsym.gt.0) then
        if(nsym.eq.3) ngener=2
        if(nsym.eq.7) ngener=3
      endif
c
c -- also put data in depository
c
c ****************************************************************
c   Normally the GEOMETRY module, and this subroutine, will be
c   invoked only ONCE per job. For a geometry optimization in
c   Z-matrix coordinates the GEOMETRY module will be called EACH
c   optimization cycle. After the first call the memory pointers
c   for the symmetry data are already set and are simply reused.
c *****************************************************************
c
      call tstival('nsyo',iexist)
      if(iexist.eq.0) then
c -- get maximum possible storage in case symmetry changes later
       call getint(natoms*7,nupair)
       nupair1 = nupair
       if(natoms.gt.natom) call getint(natom*7,nupair1)
       call getint(7,iadr)
      else
       call getival('SymNuPr',nupair)
       call getival('SymNuPr1',nupair1)
       call getival('nsyo',iadr)
      endif
c
c  -- set equivalent atom array for ALL atoms
c     ij = nupair-1
      ij = 0
      do j=1,nsym
      do i=1,natoms
      ij=ij+1
c     bl(ij) = meqatm(i,j)
      call int_to_bl( bl(nupair), ij, meqatm(i,j) )
      enddo
      enddo
c
c -- if necessary, set equivalent atoms array for REAL atoms
      If(natoms.gt.natom) Then
       Call SymFilter(natoms,natom,nsym,bl(nupair),bl(nupair1))
      EndIf
c
      call izeroit(bl(iadr), 7 )
      do i=1,nsym            ! why do this? data already in common /symm/
      call int_to_bl( bl(iadr), i, nsy(i) )
      enddo
c
c -- set values in depository
      call setival('nsym',nsym)
      call setival('ngener',ngener)
      If(iexist.eq.0) Then
c -- don't set if already there
       call setival('SymNuPr',nupair)
       call setival('SymNuPr1',nupair1)
       call setival('nsyo',iadr)
      EndIf
c
c --  reset basis symmetry data
       call tstival('ncs',ifind)
       if(ifind.eq.1) call symfunc(.true.)
c
      return
c
 900  Format(3F20.14)
 910  Format(20I4)
 920  Format(9X,I4,2X,I4)
 1000 Format(' Abelian symmetry operations: ',7(A3,2x))
c
      end
c============================================================================
      subroutine abelian(nsym,nsy,sflies)
      implicit integer(a-z)
c  This routine determines the Abelian subgroup of the molecular
c  point group and writes point group symbols to depository
c  Parameters:
c  nsym:   number of idempotent symmetry operations
c  nsy:    symmetry operations
c  sflies: character symbol of full molecular point group
      character*3 groups(8),group,group1
      character*4 sflies
      dimension nsy(7)
      data groups/'d2h','d2 ','c2v','c2h','c2 ','cs ','ci ','c1 '/
c
      iout=igetival('iout')
      if(nsym.eq.7) then
c  D2h
        group=groups(1)
      else if(nsym.eq.3) then
        if(nsy(1).le.3.and.nsy(2).le.3) then
c generators are two symmetry planes - C2v
          group=groups(3)
        else if(nsy(1).le.3.or.nsy(2).le.3) then
c  generators are one plane and one C2 - C2h
          group=groups(4)
        else if(nsy(1).gt.3.and.nsy(2).gt.3.and.nsy(2).lt.7) then
c  generators are two C2 axes - D2
          group=groups(2)
        end if
      else if(nsym.eq.1) then
        if(nsy(1).le.3) then
c one plane only - Cs
          group=groups(6)
        else if(nsy(1).lt.7) then
c one C2 axis
          group=groups(5)
        else if(nsy(1).eq.7) then
          group=groups(7)
        end if
      else
        group=groups(8)
      end if
      group1=group
      call uppercase(group1,1)
      write(iout,1400) group1
c      write(icond,1400) group
 1400 format('   Largest Abelian subgroup of the molecular',
     $         ' point group is ',A3)
c  store the point group information
      call setchval('Sflies',sflies)
      call setchval('Abelian',group)
c
      end
c============================================================================
      subroutine abeliansymm(natoms,natom,nsym,isym,meqatom,xnuc)
c  This routine enforces any Abelian symmetry that was detected
c  Arguments:
c  INTENT(IN)
c  natoms = total numer of atomic centers
c  natom  = number of real atoms
c  nsym   = number of Abelian symmetry operations (excluding identity); 1,3 or 7
c  isym   = Abelian symmetry operations
c               1 - reflection in YZ plane
c               2 - reflection in XZ plane
c               3 - reflection in XY plane
c               4 - C2 rotation about Z-axis
c               5 - C2 rotation about Y-axis
c               6 - C2 rotation about X-axis
c               7 - inversion through origin
c meqatom = list of symmetry equivalent atoms under the Abelian subgroup
c              meqatom(natoms,7)
c   INTENT(IN-OUT)
c  xnuc   = nuclear data array [xnuc(5,*): charge,x,y,z,name]
c
      implicit real*8 (a-h,o-z)
      logical last
      parameter (zero=0.0d0)
      dimension isym(7),meqatom(natoms,7),xnuc(5,natom),xx(3)
c  mirror(k,isy) gives the Abelian symmetry behavior for coordinate k
c  under symmetry operation isy: it is +/- 1.
      dimension mirror(3,7)
      data mirror(1,1),mirror(2,1),mirror(3,1),mirror(1,2),mirror(2,2),
     1 mirror(3,2),mirror(1,3),mirror(2,3),mirror(3,3),mirror(1,4),
     2 mirror(2,4),mirror(3,4),mirror(1,5),mirror(2,5),mirror(3,5),
     3 mirror(1,6),mirror(2,6),mirror(3,6),mirror(1,7),mirror(2,7),
     4 mirror(3,7)
     5 /-1,1,1, 1,-1,1, 1,1,-1, -1,-1,1, -1,1,-1, 1,-1,-1, -1,-1,-1/
      if(nsym.lt.1) RETURN
CPP
cc      iout=igetival('iout')
cc      write(iout,*) 'initial nuclear coordinates and symmetry pairs'
cc      do i=1,natom
cc        write(iout,100) i, xnuc(2,i),xnuc(3,i),xnuc(4,i),
cc     1    (meqatom(i,k),k=1,nsym)
cc  100 format(i3,2x,3f15.8,2x,7i4)
cc      end do
cc      write(iout,*) 'nsym,isym',nsym,'   ',(isym(k),k=1,nsym)
CPP
      do i=1,natom
        last=.true.
        do j=1,nsym
          if(meqatom(i,j).gt.i) last=.false.
c  perform this only for the last of the set
        end do
        if(last) then
          do k=1,3
            xx(k)=xnuc(k+1,i)
            do iop=1,nsym
              isy=isym(iop)
              ipair=meqatom(i,iop)
              xx(k)=xx(k)+mirror(k,isy)*xnuc(k+1,ipair)
            end do
            xx(k)=xx(k)/dble(nsym+1)
c  xx(k) is the averaged position over symmetry images. Now put them back
            xnuc(k+1,i)=xx(k)
            do isy=1,nsym
              ipair=meqatom(i,isy)
              xnuc(k+1,ipair)=xx(k)*mirror(k,isy)
            end do
          end do
        end if
      end do
CPP
cc      do i=1,natom
cc        write(iout,100) i, xnuc(2,i),xnuc(3,i),xnuc(4,i)
cc      end do
c
      end
c ====================================================================
c
      SUBROUTINE ChkGhost(IOut,NAtoms,AtSymb,XC,XCharg,IPRNT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine checks for, and alerts the user to, the presence of ghost atoms
C
C  ARGUMENTS
C
C  IOut    -  unit number to write to
C  NAtoms  -  total number of real atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates (assumed in atomic units)
C  XCharg  -  atomic charges
C  IPrnt   -  print flag
C
      DIMENSION XC(3,NAtoms),XCharg(NAtoms)
      Character*8 AtSymb(NAtoms)
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
      Parameter (TolZero=1.0d-8)
c
      nghost = 0
c
      DO 10 I=1,NAtoms
      If(Abs(XCharg(I)).LT.TolZero) Then
        nghost = nghost+1
        If(IPrnt.GE.2) Then
          If(nghost.EQ.1) WRITE(IOut,1000)
          X = XC(1,I)/ANTOAU
          Y = XC(2,I)/ANTOAU
          Z = XC(3,I)/ANTOAU
          WRITE(IOut,1100) I,AtSymb(I),X,Y,Z
        EndIf
      EndIf
 10   CONTINUE
C
      If(nghost.GT.0) WRITE(IOut,1200) nghost
C
      RETURN
c
 1000 FORMAT(1X,' The following are ghost atoms:')
 1100 FORMAT(1X,I4,2X,A8,3F12.6)
 1200 FORMAT(1X,' There are ',I4,' ghost atoms',/)
c
      END
