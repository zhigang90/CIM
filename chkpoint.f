c ==============================================================
c  routines used to transfer data to and from the checkpoint
c  <control> file      ** WARNING - NOT IN FINAL FORM **
c ==============================================================
c
      SUBROUTINE RdBASIS(IUnit,  NAtoms, AtSymb, XC,     ISCR,
     $                   ILST,   BASDAT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads basis set from <basis> file and loads data into
C  BASDAT and precursor INX arrays
C
C  ARGUMENTS
C
C  IUnit   -  unit number of basis file
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  ISCR    -  integer scratch array
C  ILST    -  presort of INX array
C  BASDAT  -  basis data array for primitive exponents and
C             contraction coefficients
C
C  ** NOTE: This routine reads each shell and then duplicates
C           the shell data for every atom of that type
C
C
      DIMENSION XC(3,NAtoms),ISCR(NAtoms),ILST(4,*)
      DIMENSION BASDAT(13,*)
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 Symb(NAtoms)
c ..................................................
      Character*3 TYPE(13)
      Character Char*136,Atm*4
C
      Parameter (Zero=0.0d0)
C
      Data  TYPE
     $   / 's  ','p  ','l  ','d  ','d6 ','f  ','f10','g15','h21','i28',
     $     'g  ','h  ','i  '/
C
C
C  first massage atomic symbols to remove numbers and
C  catch special symbols
C
      CALL SymbolM(NAtoms,AtSymb,Symb)
C
C  locate <$basis> keyword
C
      CALL rdcntrl(IUnit,6,'$basis',0,idum,dum,char)
C
C  initialize
C
      nprm1 = 0             ! starting counter for primitive shells/atom
      nprm = 0              ! number of primitive shells
      ncs = 0               ! number of contracted shells
      nprms = 0             ! starting primitive for current shell
      NewAtom = 0           ! counter for new atom
C
C
C  loop over atoms in basis set definition
C
 10   CONTINUE
      READ(IUnit,900) Char
c
      IF(Char(1:3).EQ.'for') THEN
cc
        If(nprm.NE.0) Then
C
C  finished a shell
C  so store shell data in ILST array
C
          ncs = ncs+1
          ILST(1,ncs) = IType                   ! shell type
          ILST(2,ncs) = nprm - nprms            ! no. of primitives
          ILST(3,ncs) = ISCR(1)                 ! atom
          ILST(4,ncs) = ngr                     ! no. of general contractions
C
C  are there more of these atoms in molecule?
C  if so, need to duplicate shell
C
          nprm2 = nprm
          ncs0 = ncs                            ! shell to be duplicated
c
          num = MAX(3,ngr+2)
          Do I=2,NA
          IAtm = ISCR(I)
c -- duplicate BASDAT (primitive shells)
          Do J=nprm1,nprm2
          nprm = nprm+1
          Do K=1,num
          BASDAT(K,nprm)  = BASDAT(K,J)
          EndDo
          BASDAT(11,nprm) = XC(1,IAtm)
          BASDAT(12,nprm) = XC(2,IAtm)
          BASDAT(13,nprm) = XC(3,IAtm)
          EndDo
c -- duplicate ILST (contracted shells)
          ncs = ncs+1
          ILST(1,ncs) = ILST(1,ncs0)
          ILST(2,ncs) = ILST(2,ncs0)
          ILST(3,ncs) = IAtm
          ILST(4,ncs) = ILST(4,ncs0)
          EndDo
c
          nprms = nprm
        EndIf
C
C  started a new atom
C  how many of the new atom type are there in the system?
C
        NewAtom = 1
        Atm = Char(11:14)                   ! this implies fixed format
        NA = 0
        Do IAtm=1,NAtoms
        If(Symb(IAtm)(1:4).EQ.Atm) Then
          NA = NA+1
          ISCR(NA)=IAtm
        EndIf
        EndDo
        GO TO 10
cc
      ELSE IF(Char(1:1).NE.' ') THEN
C
C  start of a new shell on current atom type
C  (or end of basis)
C
        If(nprm.NE.0.AND.NewAtom.EQ.0) Then
C
C  finished a shell
C  so store shell data in ILST array
C
          ncs = ncs+1
          ILST(1,ncs) = IType                   ! shell type
          ILST(2,ncs) = nprm - nprms            ! no. of primitives
          ILST(3,ncs) = ISCR(1)                 ! atom
          ILST(4,ncs) = ngr                     ! no. of general contractions
C
C  are there more of these atoms in molecule?
C  if so, need to duplicate shell
C
          nprm2 = nprm
          ncs0 = ncs                            ! shell to be duplicated
c
          num = MAX(3,ngr+2)
          Do I=2,NA
          IAtm = ISCR(I)
c -- duplicate BASDAT (primitive shells)
          Do J=nprm1,nprm2
          nprm = nprm+1
          Do K=1,num
          BASDAT(K,nprm)  = BASDAT(K,J)
          EndDo
          BASDAT(11,nprm) = XC(1,IAtm)
          BASDAT(12,nprm) = XC(2,IAtm)
          BASDAT(13,nprm) = XC(3,IAtm)
          EndDo
c -- duplicate ILST (contracted shells)
          ncs = ncs+1
          ILST(1,ncs) = ILST(1,ncs0)
          ILST(2,ncs) = ILST(2,ncs0)
          ILST(3,ncs) = IAtm
          ILST(4,ncs) = ILST(4,ncs0)
          EndDo
c
          nprms = nprm
        EndIf
c
        If(Char(1:1).EQ.'$') RETURN         ! we're done
c...............................
c
        nprm1 = nprm+1
        nprms = nprm
C
C  started new shell
C  what type is it?
C
        Do I=1,13
        If(Char(1:3).EQ.TYPE(I)) Then
         IType = I
         GO TO 20
        EndIf
        EndDo
C
C  should not get here
C
        Char = 'Unknown basis function type: '//Char(1:3)
        Call nerror(1,'File IO routine <RdBASIS>',Char(1:32),0,0)
      ENDIF
c
 20   CONTINUE
C
C  read exponent and coefficient
C
      ngr = 0
      nprm = nprm+1
      If(IType.EQ.3) Then
       READ(Char(5:80),*) BASDAT(1,nprm),BASDAT(2,nprm),BASDAT(3,nprm)
      Else
       READ(Char(5:38),*) BASDAT(1,nprm),BASDAT(2,nprm)
C
C  now check for general contractions
C  i.e. third BASDAT entry non-zero for other than L shell
C
       Do igr=3,10
       ii = 39+12*(igr-3)
       READ(Char(ii:ii+11),*,END=25,ERR=25) BASDAT(igr,nprm)
       If(BASDAT(igr,nprm).EQ.Zero) GO TO 25
       EndDo
 25    ngr = igr-3
      EndIf
C
C  put coordinates of atomic centre in BASDAT
C  (apparently obsolescent, but still used)
C
      IAtm = ISCR(1)
      BASDAT(11,nprm) = XC(1,IAtm)
      BASDAT(12,nprm) = XC(2,IAtm)
      BASDAT(13,nprm) = XC(3,IAtm)
C
C  read another card
C
      NewAtom = 0
      GO TO 10
c
  900 Format(A136)
c
      END
c ====================================================================
c
      SUBROUTINE RdBASIS1(IUnit,  NAtoms, AtSymb, XC,     ILST,
     $                    BASDAT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads basis set from <basis> file and loads data into
C  BASDAT and precursor INX arrays
C
C  ARGUMENTS
C
C  IUnit   -  unit number of basis file
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  ILST    -  presort of INX array
C  BASDAT  -  basis data array for primitive exponents and
C             contraction coefficients
C
C  ** NOTE: This routine reads in the basis for each atom
C           in the original geometry input order
C     <CURRENTLY NOT USED>
C
C
      DIMENSION XC(3,NAtoms),ILST(4,*)
      DIMENSION BASDAT(13,*)
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 Symb(NAtoms)
c ..................................................
      Character*3 TYPE(10)
      Character Char*120,Atm*8
C
      Parameter (Zero=0.0d0)
C
      Data  TYPE
     $   / 's  ','p  ','l  ','d  ','d6 ','f  ','f10','g15','h21','i28'/
C
C
C  first massage atomic symbols to remove numbers and
C  catch special symbols
C
      CALL SymbolM(NAtoms,AtSymb,Symb)
C
C  initialize
C
      nprm = 0              ! number of primitive shells
      ncs = 0               ! number of contracted shells
C
C
C  loop over atoms in system
C
      DO 50 IAtm=1,NAtoms
      Atm = Symb(IAtm)     ! this is the atom we are dealing with
      ncs0 = ncs
      nprms = nprm
C
C  locate <$basis> keyword
C
      CALL rdcntrl(IUnit,6,'$basis',0,idum,dum,char)
C
 10   CONTINUE
      READ(IUnit,900,END=40) Char
      If(Char(1:3).NE.'for'.OR.Char(11:18).NE.Atm) GO TO 10
C
C  this is the atom we want
C  read basis
C
 20    CONTINUE
       READ(IUnit,900) Char
c
       If(Char(1:3).EQ.'for'.OR.Char(1:4).EQ.'$end') GO TO 40
c
       IF(Char(1:1).NE.' ') THEN
C
C  new shell
C  check if this atom already has shells
C  if so, store shell data in ILST array
C
        If(nprms.NE.nprm) Then
         ncs = ncs+1
       ILST(1,ncs) = IType                   ! shell type
       ILST(2,ncs) = nprm - nprms            ! no. of primitives
       ILST(3,ncs) = IAtm                    ! atom
       ILST(4,ncs) = ngr
      EndIf
C
C  started new shell
C  what type is it?
C
        nprms = nprm
        Do I=1,10
        If(Char(1:3).EQ.TYPE(I)) Then
         IType = I
         GO TO 30
        EndIf
        EndDo
C
C  should not get here
C
        Char = 'Unknown basis function type: '//Char(1:3)
        Call nerror(1,'File IO routine <RdBASIS>',Char(1:32),0,0)
       ENDIF
c
 30   CONTINUE
      ngr = 0
      nprm = nprm+1
      If(IType.EQ.3) Then
       READ(Char(5:80),*) BASDAT(1,nprm),BASDAT(2,nprm),BASDAT(3,nprm)
      Else
       READ(Char(5:32),*) BASDAT(1,nprm),BASDAT(2,nprm)
C
C  now check for general contractions
C  i.e. third BASDAT entry non-zero for other than L shell
C
       Do igr=3,10
       ii = 33+12*(igr-3)
       READ(Char(ii:ii+11),*,END=25) BASDAT(igr,nprm)
       If(BASDAT(igr,nprm).EQ.Zero) GO TO 25
       EndDo
 25    ngr = igr-3
      EndIf
C
C  put coordinates of atomic centre in BASDAT
C  (apparently obsolescent, but still used)
C
      BASDAT(11,nprm) = XC(1,IAtm)
      BASDAT(12,nprm) = XC(2,IAtm)
      BASDAT(13,nprm) = XC(3,IAtm)
C
C  read another card
C
      GO TO 20
c
 40   CONTINUE
C
C  finished this atom
C  write final ILST array data
C
      If(nprms.NE.0) Then
       ncs = ncs+1
       ILST(1,ncs) = IType                   ! shell type
       ILST(2,ncs) = nprm - nprms            ! no. of primitives
       ILST(3,ncs) = IAtm                    ! atom
       ILST(4,ncs) = ngr
      EndIf
c
      REWIND IUnit
c
      If(ncs.EQ.ncs0) Then
C
C  this atom has no basis!
C
       Char = 'Atom Shown has no basis!: '//Atm
       Call nerror(1,'File IO routine <RdBASIS>',Char(1:34),0,0)
      EndIf
c
 50   CONTINUE
C
      RETURN
c
  900 Format(A120)
c
      END
c ====================================================================
c
      SUBROUTINE RdBackBasis(NAtoms, NShell, NPrim, INX, BASDAT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  HOPEFULLY A TEMPORARY ROUTINE UNTIL <reorder> IS REWRITTEN
C  Reads back the basis set just written using <WrBASIS>
C  This is to ensure that routine <reorder> gets the SAME
C  input each time it is called and hence should generate
C  the SAME Wolinski ordering
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NShell  -  number of shells
C  NPrim   -  number of primitive gaussians
C  INX     -  Texas INX array (integer basis data)
C  BASDAT  -  exponents and contraction coefficients
C
C
      DIMENSION INX(12,NShell),BASDAT(13,NPrim)
c ..................................................
c -- automatic allocation of arrays in F90
      REAL*8 XC(3,NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
      INTEGER ISCR(NAtoms),ILST(4,NShell)
c ..................................................
C
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      CHARACTER*256 jobname
      Common /job/jobname,lenJ
C
C
C  read the Cartesian coordinates and atomic symbols
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoord(IUnit,NAtoms,AtSymb,XC,-1,jnk)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  open the <basis> file and read back the basis
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdBASIS(IUnit,  NAtoms, AtSymb, XC,     ISCR,
     $             ILST,   BASDAT)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  sort the basis (by shell)
C
      CALL SortBAS1(NAtoms,NShell,ILST,INX)
C
      RETURN
      END
c ====================================================================
c
      SUBROUTINE rdcntrl(IUnit,Length,Key,IType,INum,RNum,CNum)
      IMPLICIT INTEGER(A-Z)
C
C  Finds location of keyword "Key" in <control> file on unit IUnit,
C  optionally reads a value or string and returns with read
C  pointer at beginning of next record
C
C  ARGUMENTS
C
C  IUnit   -  unit number of <control> file (assumed open)
C  Length  -  length of string
C  Key     -  string to be found
C  IType   -  control of read
C              0 - locate string and return
C              1 - read integer value
C              2 - read real value
C              3 - read character string (unspecified length)
C
C  on exit
C
C  INum    -  contains integer value     (if IType=1)
C             contains length of string  (if IType=3)
C  RNum    -  contains real value        (if IType=2)
C  CNum    -  contains character string  (if IType=3)
C             (maximum 20 characters)
C
C
      REAL*8 RNum
      CHARACTER*(*) Key
      CHARACTER CNum*20,Char*80
C
C
      REWIND IUnit
 10   CONTINUE
      READ(IUnit,900,End=95) Char
      If(Char(1:Length).EQ.Key) GO TO 20
      GO TO 10
c
 20   CONTINUE
      IF(IType.EQ.0) THEN
      ELSE IF(IType.EQ.1) THEN
        READ(Char(Length+1:80),*) INum
      ELSE IF(IType.EQ.2) THEN
        READ(Char(Length+1:80),*) RNum
      ELSE
        CNum = Char(Length+1:Length+20)
        CALL RmBlank(20,CNum,INum)
      ENDIF
c
      RETURN
c
 95   CONTINUE
      Char = 'Unable to find Keyword '//Key//' in <control> file'
      Call nerror(2,'File IO routine <rdcntrl>',Char(1:41+Length),0,0)
c
  900 Format(A80)
c
      END
c ===================================================================
c
      SUBROUTINE fdcntrl(IUnit,Length,Key,IEnd)
      IMPLICIT INTEGER(A-Z)
C
C  Finds location of keyword "Key" in <control> file on unit IUnit,
C  if cannot find then stops at end-of-file ($end)
C
C  ARGUMENTS
C
C  IUnit   -  unit number of <control> file (assumed open)
C  Length  -  length of string
C  Key     -  string to be found
C  IEnd    -  integer flag
C              0 - found keyword
C              1 - found end-of-file
C
C
      CHARACTER*(*) Key
      CHARACTER Char*80
C
C
      IEnd = 1
c
      REWIND IUnit
 10   CONTINUE
      READ(IUnit,900,End=95) Char
      If(Char(1:Length).EQ.Key) GO TO 20
      If(Char(1:4).EQ.'$end') GO TO 30
      GO TO 10
c
 20   CONTINUE
      IEnd = 0
 30   RETURN
c
 95   CONTINUE
      Call nerror(3,'File IO routine <fdcntrl>',
     $       'Unable to find $end marker in <control> file',0,0)
c
  900 Format(A80)
c
      END
c ===================================================================
c
      SUBROUTINE RmBlank(Len,CNum,RealLen)
      IMPLICIT INTEGER(A-Z)
C
C  Removes ALL blanks from a character string
C
C  ARGUMENTS
C
C  Len     -  total length of initial string
C  CNum    -  character string
C  RealLen -  actual length after ALL blanks removed
C
      Character*(*) CNum
      Character*1 Blank
c
      Data Blank/' '/
C
      RealLen = 0
c
      DO 10 I=1,Len
      If(CNum(I:I).NE.Blank) Then
      RealLen = RealLen+1
      CNum(RealLen:RealLen) = CNum(I:I)
      EndIf
 10   CONTINUE
c
      DO 20 I=RealLen+1,Len
      CNum(I:I) = Blank
 20   CONTINUE
C
      RETURN
      END
c ====================================================================
c
      SUBROUTINE RdCoord(IUnit,NAtoms,AtSymb,XC,NMol,IMOL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Routine to read Cartesian coordinates from <control> file
C  ** Short Read (No charges/masses) **
C
C  ARGUMENTS
C
C  IUnit   -  unit number of control file (assumed open)
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  current geometry (Cartesian coordinates)
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C             If -1 on entry, IMOL array will NOT be filled
C  IMOL    -  pointers to start/end of molecules in XC array
C
C
      DIMENSION XC(3,NAtoms),IMOL(*)
      CHARACTER*8 AtSymb(NAtoms)
      Character*80 Char
C
C
C  locate the $coordinate section
C
      CALL rdcntrl(IUnit,12,'$coordinates',0,idum,dum,char)
C
C  now read coordinates
C
      IMOL(1) = 0         ! always true
c
      IF(NMol.EQ.1) THEN
       DO 10 IAtm=1,NAtoms
       READ(IUnit,910) AtSymb(IAtm),XC(1,IAtm),XC(2,IAtm),XC(3,IAtm)
 10    CONTINUE
       IMOL(2) = NAtoms
      ELSE
       IAtm = 0
       LMol = 0
 20    READ(IUnit,920,END=95) Char
       If(Char(1:1).EQ.'$') Then
        LMol = LMol+1
        If(NMol.GE.0) IMOL(LMol+1) = IAtm
        If(Char(2:4).EQ.'end') GO TO 90
       Else
        IAtm = IAtm+1
        If(IAtm.GT.NAtoms) GO TO 90
        READ(Char,910) AtSymb(IAtm),XC(1,IAtm),XC(2,IAtm),XC(3,IAtm)
       EndIf
       GO TO 20
      ENDIF
C
 90   CONTINUE
cc      write(6,*) ' Debug check in <RdCoord>  NAtoms:',natoms,' IAtm:',
cc     $             iatm,' NMol:',nmol,' LMol:',lmol
cc      if(nmol.ne.-1) then
cc       write(6,*) ' IMOL array is:'
cc       write(6,*) (imol(i),i=1,nmol+1)
cc      endif
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no $end section marker found
      Call nerror(4,'File IO routine <RdCoord>',
     $  'No End-Of-File marker ($end) found on <coord> file',0,0)
c
  910 Format(A8,2X,3F20.14)
  920 Format(A80)
c
      END
c ======================================================================
c
      SUBROUTINE RdCoordF(IUnit,  NAtoms, AtSymb, XC,     NMol,
     $                    IMOL,   XCharg, XMass)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Routine to read Cartesian coordinates from <control> file
C  ** Full Read (including charges/masses) **
C
C  ARGUMENTS
C
C  IUnit   -  unit number of control file (assumed open)
C  NATOMS  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  current geometry (Cartesian coordinates)
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C             If -1 on entry, IMOL array will NOT be filled
C  IMOL    -  pointers to start/end of molecules in XC array
C  XCharg  -  atomic charges (can be user-defined)
C  XMass   -  atomic masses
C
C
      REAL*8 XC(3,NAtoms),XCharg(NAtoms),XMass(NAtoms)
      DIMENSION IMOL(*)
      CHARACTER*8 AtSymb(NAtoms)
      Character*110 Char
C
C
C  locate the $coordinate section
C
      CALL rdcntrl(IUnit,12,'$coordinates',0,idum,dum,char)
C
C  now read coordinates
C
      IMOL(1) = 0         ! always true
c
      IF(NMol.EQ.1) THEN
       DO 10 IAtm=1,NAtoms
       READ(IUnit,910) AtSymb(IAtm),XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),
     $                 XCharg(IAtm),XMass(IAtm)
 10    CONTINUE
       IMOL(2) = NAtoms
      ELSE
       IAtm = 0
       LMol = 0
       IMOL(1) = 0
 20    READ(IUnit,920,END=95) Char
       If(Char(1:1).EQ.'$') Then
        LMol = LMol+1
        If(NMol.GE.0) IMOL(LMol+1) = IAtm
        If(Char(2:4).EQ.'end') GO TO 90
       Else
        IAtm = IAtm+1
        If(IAtm.GT.NAtoms) GO TO 90
        READ(Char,910) AtSymb(IAtm),XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),
     $                 XCharg(IAtm),XMass(IAtm)
       EndIf
       GO TO 20
      ENDIF
C
 90   CONTINUE
cc      write(6,*) ' Debug check in <RdCoordF>  NAtoms:',natoms,' IAtm:',
cc     $             iatm,' NMol:',nmol,' LMol:',lmol
cc      if(nmol.ne.-1) then
cc       write(6,*) ' IMOL array is:'
cc       write(6,*) (imol(i),i=1,nmol+1)
cc      endif
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no $end section marker found
      Call nerror(4,'File IO routine <RdCoordF>',
     $  'No End-Of-File marker ($end) found on <coord> file',0,0)
c
  910 Format(A8,2X,5F20.14)
  920 Format(A110)
c
      END
c ======================================================================
c
      SUBROUTINE RdDeriv(NAt3,DipD,PolD,Dipole,Polar)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Attempt to read Cartesian derivatives from <deriv> file
C
C  ARGUMENTS
C
C  NAt3    -  3 x number of atoms
C  DipD    -  dipole derivatives in Cartesian coordinates
C  PolD    -  polarizability derivatives in Cartesian coordinates
C  Dipole  -  logical flag
C             set to .TRUE. on exit if dipole derivatives found
C  Polar   -  logical flag
C             set to .TRUE. on exit if polarizability derivatives found
C
C
      DIMENSION DipD(NAt3,3),PolD(NAt3,6)
      CHARACTER CHAR*80
      LOGICAL Dipole,Polar
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
      Dipole = .FALSE.
      Polar = .FALSE.
C
C  open file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.deriv',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  look for dipole derivatives
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      If(CHAR(1:7).NE.'$dipole') GO TO 10
C
C  dipole derivatives found
C
      DO 20 I=1,NAt3
      READ(40,910) DipD(I,1),DipD(I,2),DipD(I,3)
 20   CONTINUE
      Dipole = .TRUE.
C
C  look for polarizability derivatives
C
 30   CONTINUE
      READ(40,900,END=95) CHAR
      If(CHAR(1:15).NE.'$polarizability') GO TO 30
C
C  polarizability derivatives found
C
      DO 40 I=1,NAt3
      READ(40,910) PolD(I,1),PolD(I,2),PolD(I,3)
      READ(40,910) PolD(I,4),PolD(I,5),PolD(I,6)
 40   CONTINUE
      Polar = .TRUE.

C  close file and exit
C
      IErr = 0
      CLOSE (UNIT=40,STATUS='KEEP')
C
 95   CONTINUE
      RETURN
c
  900 Format(A80)
  910 Format(10X,3F20.14)
c
      END
c ====================================================================
c
      SUBROUTINE RdDFT(IUnit,dft)
      IMPLICIT INTEGER(A-Z)
      Character*20 cdum
C
C  reads dft flag from control file
C
      dft = 0
      call fdcntrl(IUnit,4,'$dft',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,'(a20)') cdum
        READ(cdum(5:20),*) dft
C
C  ** WARNING **
C     need to read ACM coefficients here if set by user
C     (maybe put into common block /dft/ ?
C
      EndIf
C
      RETURN
      END
c ====================================================================
c
      SUBROUTINE RdDISP(IUnit,noabc,tz,VALUES,disp)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VALUES(5)
      LOGICAL noabc,tz,disp
C
C  reads dispersion parameters from <control> file
C
      disp = .false.
      call fdcntrl(IUnit,10,'$paramdisp',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,910) (VALUES(J),J=1,5)
       READ(IUnit,920) iabc,itz
       noabc = iabc.eq.1
       tz = itz.eq.1
       disp = .true.
      EndIf
C
      RETURN
c
  910 Format(10X,5(2X,F10.5))
  920 Format(7X,I2,8X,I2)
c
      END
c ====================================================================
c
      SUBROUTINE RdDIP(IUnit,Dip,Dipole)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Dip(3)
      LOGICAL Dipole
C
C  reads dipole moment (X, Y, Z components) from <control> file
C
      dipole = .false.
      call fdcntrl(IUnit,7,'$dipole',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,910) Dip(1),Dip(2),Dip(3)
       Dipole = .true.
      EndIf
C
      RETURN
c
  910 Format(7X,3F20.14)
c
      END
c ====================================================================
c
      SUBROUTINE RdField(IUnit,Field,efld)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Field(3)
      LOGICAL Efld
C
C  reads field components (X, Y, Z components) from <control> file
C
      efld = .false.
      call fdcntrl(IUnit,6,'$field',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,910) Field(1),Field(2),Field(3)
       efld = .true.
      EndIf
C
      RETURN
c
  910 Format(8X,3F12.6)
c
      END
c ====================================================================
c
      SUBROUTINE RdGEOM(NAtoms, Ndum1,  Ndum2,  XNuc,   AtMASS,
     $                  NMol,   IMOL,   Charg,  IMult)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads geometry data from an old <control> file into XNuc array
C  (involving old texas format)
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms (including dummies)
C  Ndum1   -  number of dummy atoms with charge
C  Ndum2   -  number of dummy atoms without charge
C   (the above three arguments may be modified on exit)
C  XNuc    -  Cartesian coordinates      (old Texas format)
C             (also includes atomic symbols and numbers)
C  AtMASS  -  atomic masses
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C             If -1 on entry, IMOL array will NOT be filled
C  IMOL    -  pointers to start/end of molecules in XC array
C  Charg   -  charge
C  IMult   -  multiplicity
C
C
      DIMENSION XNuc(5,NAtoms),AtMASS(NAtoms),IMOL(*)
      Character Char*110,cdum*20,name*8
      Equivalence (xname,name)
C
      PARAMETER (IUnit=1)
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  first read from the <control> file
C    total number of atoms
C    number of dummy atoms (if any)
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      call rdcntrl(IUnit,7,'$charge',2,Idum,Charg,cdum)
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  now open and read the <coord> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD',ERR=96)
C
C  first locate the $coordinate section
C
      CALL rdcntrl(IUnit,12,'$coordinates',0,idum,dum,cdum)
C
C  read coordinates
C
      IMOL(1) = 0          ! always true
c
      IF(NMol.EQ.1) THEN
       DO 10 I=1,NAtoms
       READ(IUnit,910) name,XNuc(2,I),XNuc(3,I),XNuc(4,I),XNuc(1,I),
     $                 AtMASS(I)
       XNuc(5,I) = xname
 10    CONTINUE
       IMOL(2) = NAtoms
      ELSE
       I = 0
       LMol = 0
 20    READ(IUnit,920) Char
       If(Char(1:1).EQ.'$') Then
        LMol = LMol+1
        If(NMol.GE.0) IMOL(LMol+1) = I
        If(Char(2:4).EQ.'end') GO TO 90
       Else
        I = I+1
        If(I.GT.NAtoms) GO TO 90
        READ(Char,910) name,XNuc(2,I),XNuc(3,I),XNuc(4,I),XNuc(1,I),
     $                 AtMASS(I)
        XNuc(5,I) = xname
       EndIf
       GO TO 20
      ENDIF
c
 90   CONTINUE
      CLOSE (UNIT=IUnit,STATUS='KEEP')
      If(NMol.NE.-1) NMol=LMol
C
C  Make sure that the nuclear charges are set correctly,
C  in case pseudopotentials are or were involved
C
      call setnucchg(IUnit,XNuc,NAtoms)
c
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE
      Call nerror(9,'File IO routine <RdGEOM>',
     $     'GEOM=READ Specified but old <control> File Does Not Exist',
     $      0,0)
c
 96   CONTINUE
      Call nerror(10,'File IO routine <RdGEOM>',
     $     'GEOM=READ Specified but old <coord> File Does Not Exist',
     $      0,0)
c
  900 Format(9X,I4,2X,I4)
  910 Format(A8,2X,5F20.14)
  920 Format(A110)
c
      END

c ====================================================================

      subroutine setnucchg( iunit, xnuc, natoms )
      
           ! This routine makes sure that the atomic charges
           ! are consistent with the basis set currently in use.
           ! This is needed in case pseudopotentials were used in
           ! the previous run, or are currently in use

      use memory
      implicit none

      integer iunit, natoms
      real*8 xnuc(5,natoms)

      integer lenj
      character*256 jobname
      Common /job/jobname,lenJ

      logical exx
      integer idum, npspold, npsp, iiecpdat, idef, i, icentr
      integer maxlpsp, ndum, iadr2
      real*8 rdum
      character*1 cdum

           ! get number of ECPs in the current basis

      call tstival( 'npsp', idef )
      if( idef .eq. 0 ) call setival( 'npsp', 0 )
      call getival( 'npsp', npsp )
           
           ! set charges for current ECP centers

      if( npsp .ne. 0 )  then
        call getival( 'iiecpdat', iiecpdat )
        call getival( 'maxlpsp', maxlpsp )
        do i=1, npsp
          call get_psp_centre( bl(iiecpdat), maxlpsp, npsp, i, icentr)
          call setcharge( icentr, bl(iiecpdat), maxlpsp, npsp,
     &                    xnuc, natoms )
        enddo
      endif

           ! get number of ECPs in the old basis

      npspold = 0
      inquire( file=jobname(1:lenj)//'.basis', exist=exx )
      if ( exx ) then
        open( unit=iunit, file=jobname(1:lenj)//'.basis',
     $        form='formatted', status='old', err=96 )
        call fdcntrl( iunit, 7, '$npseud', idum )
        if ( idum .eq. 0 ) then
          call rdcntrl( iunit, 7, '$npseud', 1, npspold, rdum, cdum )
        endif
      endif

           ! set charges for old ECP centers

      if( npspold .ne. 0 ) then
        do i=1, npspold
          read(iunit,*) icentr  ! read ECP center from .basis file

          call setcharge( icentr, bl(iiecpdat), maxlpsp, npsp,
     &                    xnuc, natoms )
        enddo
      endif

      close( iunit, err=96 )
96    continue

           ! for good measure, repeat the procedure for the .basis2 file
           ! get number of ECPs in the old basis

      npspold = 0
      inquire( file=jobname(1:lenj)//'.basis2', exist=exx )
      if ( exx ) then
        open( unit=iunit, file=jobname(1:lenj)//'.basis2',
     $        form='formatted', status='old', err=97 )
        call fdcntrl( iunit, 7, '$npseud', idum )
        if ( idum .eq. 0 ) then
          call rdcntrl( iunit, 7, '$npseud', 1, npspold, rdum, cdum )
        endif
      endif

           ! set charges for old ECP centers

      if( npspold .ne. 0 ) then
        do i=1, npspold
          read(iunit,*) icentr  ! read ECP center from .basis file

          call setcharge( icentr, bl(iiecpdat), maxlpsp, npsp,
     &                    xnuc, natoms )
        enddo
      endif

      close( iunit, err=97 )
97    continue

           ! update total charge (if possible)

      call tstival( 'ndum', idef )
      if ( idef .ne. 0 ) then
        call getival( 'ndum', ndum )
        call tstival( 'mass', idef )
        if ( idef .ne. 0 ) then
          call getival( 'mass', iadr2 )
          call nucparam( natoms-ndum, xnuc, bl(iadr2) )
        endif
      endif

      end

c ====================================================================

      subroutine setcharge( icentr, iecpdat, maxlpsp, npsp,
     &                      xnuc, natoms )

           ! set the nuclear charge for centre icentr
 
           ! icentr       centre to be examined
           ! iecpdat      ecp information, used only if npsp .ne. 0
           ! maxlpsp      first dimension of iecpdat array
           ! npsp         number of pseudopotentials
           ! xnuc         nuclear data, in output will contain
           !              the appropriate nuclear charge for
           !              centre icentr
           ! natoms       number of atoms

      implicit none

      integer icentr, npsp, maxlpsp, natoms
      integer iecpdat(maxlpsp+6,npsp)
      real*8 xnuc(5,natoms)

      integer  i, izcore, iatn
      character*8 atsymb
      real*8 xname
      Equivalence (xname,atsymb)

           ! get number of core electrons

      izcore = 0
      if ( npsp .ne. 0 ) then
        do i = 1, npsp
          if( iecpdat(3,i) .eq. icentr ) izcore = iecpdat(2,i)
        enddo
      endif

          ! get the atomic number from the atomic symbol

      xname = xnuc(5,icentr)
      call nugrep( atsymb(1:2), iatn )

          ! set charge

      xnuc(1,icentr) = real( iatn - izcore )

      end

c ====================================================================
c
      SUBROUTINE RdGRAD(NATOMS,GC,FStat)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads gradient from <grad> file
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  GC      -  on exit contains current gradient
C  FStat   -  status of file after read (save/delete)
C
C
      REAL*8 GC(3,NATOMS)
      CHARACTER*80 CHAR
      Character*4 FStat
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open GRAD file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.grad',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $gradient section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      IF(CHAR(1:9).NE.'$gradient') GO TO 10
C
C  read gradient
C
      DO 20 IAtm=1,NATOMS
      READ(40,910) GC(1,IAtm),GC(2,IAtm),GC(3,IAtm)
 20   CONTINUE
C
C  check for end of data
C  (just look for "$")
C
      READ(40,900) CHAR
      If(CHAR(1:1).NE.'$') GO TO 97
c
      IF(FStat.EQ.'save'.OR.FStat.EQ.'keep') THEN
        CLOSE (UNIT=40,STATUS='KEEP')
      ELSE
        CLOSE (UNIT=40,STATUS='DELETE')
      ENDIF
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE
      Call nerror(6,'File IO routine <RdGRAD>',
     $     'Unable to find <grad> file!',0,0)
c
 96   CONTINUE
      Call nerror(7,'File IO routine <RdGRAD>',
     $      'No $gradient section found on <grad> file!',0,0)
c
 97   CONTINUE
      Call nerror(8,'File IO routine <RdGRAD>',
     $  'No End-Of-File marker ($end) found on <grad> file',0,0)
c
 900  Format(A80)
 910  Format(10X,3F20.14)
c
      END
c =====================================================================
c
      SUBROUTINE RdPATH(IUnit,coord,DMax,ISign,MaxCyc,DTol)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads Reaction Path information from <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to read from
C  coord   -  coordinate type
C              cart - Cartesian coordinates
C              zmat - z-matrix coordinates
C              mwgt - mass-weighted Cartesians
C  DMax    -  maximum step size
C  ISign   -  sign of Hessian mode on first step
C  MaxCyc  -  maximum number of steps
C  DTol    -  displacement criterion for stopping search
C
      CHARACTER*4 coord
C
      CALL fdcntrl(IUnit,5,'$path',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,1000) coord,DMax,ISign,MaxCyc,DTol
      EndIf
C
      RETURN
c
 1000 FORMAT(7X,A4,2X,F9.6,2X,I2,2X,I4,2X,F9.6)
c
      END
c =====================================================================
c
      SUBROUTINE RdPOL(IUnit,field,IType)
      IMPLICIT INTEGER(A-Z)
C
C  reads numerical polarizability information from <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to read from
C  field   -  field strength (atomic units)
C  IType   -  job type
C               1 - calculate polarizabilities
C               2 - calculate polarizability derivatives
C
      REAL*8 field
C
      CALL fdcntrl(IUnit,4,'$pol',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,1000) field,IType
      EndIf
C
      RETURN
c
 1000 FORMAT(13X,F12.6,7X,I2)
c
      END
c =====================================================================
c
      SUBROUTINE RdPropG(IUnit,factor,r0f,LMax,NRad,NAng)
      IMPLICIT INTEGER(A-Z)
C
C  Reads numerical grid information from <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to read from
C  factor  -  radial grid factor
C  r0f     -  Chipman radial factor
C  LMax    -  maximum angular momentum for spherical harmonics
C  NRad    -  number of radial points
C  NAng    -  number of angular points
C
      REAL*8 factor,r0f
C
      CALL fdcntrl(IUnit,10,'$prop_grid',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,1000) factor,r0f,LMax,NRad,NAng
      EndIf
C
      RETURN
c
 1000 FORMAT(16x,f7.4,6x,f7.4,7x,i5,2x,2i5)
c
      END
c ====================================================================
c
      SUBROUTINE RdSCAN0(IUnit,  zmat,   NPComp, force,  found)
C
C  Reads scan type and number of primitives in composite constraint
C  (if any) from <control> file
C
      Character Char*80
      Logical zmat,force,found
c
      zmat = .False.
      force = .False.
      NPComp = 0
c
      call fdcntrl(IUnit,5,'$scan',IEnd)
      IF(IEnd.EQ.0) THEN
        found = .True.
        backspace IUnit
        READ(IUnit,900) Char
        If(Char(8:11).EQ.'zmat') Then
          zmat = .True.
          RETURN
        Else If(Char(8:11).EQ.'extf') Then
          force = .True.
          RETURN
        Else If(Char(8:11).EQ.'comp') Then
          READ(IUnit,900) Char
 10       CONTINUE
          READ(IUnit,900) Char
          If(Char(1:4).EQ.'#end') RETURN
          NPComp = NPComp+1
          GO TO 10
        Else
          found = .True.
          RETURN
        EndIf
      ENDIF
c
c -- if reached here then no $scan keyword found
      found = .False.
      RETURN
c
 900  Format(A80)
c
      END
c ====================================================================
c
      SUBROUTINE RdSCAN(IUnit,  zmat,   Coord,  I,      J,
     $                  K,      L,      NPComp, IPCOMP, PCWght,
     $                  R1,     R2,     Step,   force)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads Data for potential scan from <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number for (open) control file
C  zmat    -  logical flag indicating Z-Matrix scan
C  Coord   -  internal coordinate type or Z-Matrix variable to scan
C  I,J,K,L -  atoms defining internal coordinate
C  NPComp  -  number of primitives in composite coordinate
C  IPCOMP  -  constraint type and atoms involved in constraint
C               IPComp(1,I) - constraint type
C               IPComp(2,I) to IPComp(5,I) - atoms in constraint
C  PCWght  -  weight of each primitive in composite constraint
C  R1      -  starting value for coordinate scan
C  R2      -  final value for coordinate scan
C  Step    -  step length for scan
C  force   -  logical variable indicating external force scan
C
C
      DIMENSION IPCOMP(5,*),PCWght(*)
      Character Char*80,Coord*4
      Logical zmat,force
C
      zmat = .False.
      force = .False.
      nt = 0
c
      call fdcntrl(IUnit,5,'$scan',IEnd)
      IF(IEnd.EQ.0) THEN
       backspace IUnit
       READ(IUnit,900) Char
       If(Char(8:11).EQ.'zmat') Then
        READ(Char,910) Coord,R1,R2,Step
        zmat = .True.
       Else If(Char(8:11).EQ.'extf') Then
        READ(Char,930) I,J,R1,R2,Step
        force = .True.
       Else If(Char(8:11).EQ.'comp') Then
        Coord = Char(8:11)
        READ(Char,915) R1,R2,Step
        READ(IUnit,900) Char
 10     CONTINUE
        READ(IUnit,900) Char
        If(Char(1:4).EQ.'#end') GO TO 95
        nt = nt+1
        READ(Char,925) IPCOMP(1,nt),(IPCOMP(J,nt),J=2,5),PCWght(nt)
        GO TO 10
       Else
        READ(Char,920) Coord,I,J,K,L,R1,R2,Step
        zmat = .False.
       EndIf
      ENDIF
C
 95   CONTINUE
      NPComp = nt
      RETURN
c
  900 Format(A80)
  910 Format(13X,A4,2X,3F12.6)
  915 Format(15X,3F12.6)
  920 Format(7X,A4,2X,4I5,2X,3F12.6)
  925 Format(1X,I4,2X,4I4,F12.6)
  930 Format(13X,2I4,2X,3F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE RdSYM(Symflag,NATOMS, RM,     GROUP,  NTrans,
     $                 NDEG,   NQ,     IUNQ,   TRANS,  NEqATM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads symmetry information from the <jobname>.sym file
C
C  Symflag -  logical flag indicating use of symmetry
C  NATOMS  -  number of real atoms
C  RM      -  rotation matrix used to reorient to "standard orientation"
C  NTrans  -  number of symmetry operations
C  GROUP   -  molecular point group
C  NDEG    -  number of internal degrees of freedom
C  NQ      -  number of symmetry unique atoms
C  IUNQ    -  list of symmetry unique atoms
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C
      REAL*8 RM(3,3),TRANS(3,3,NTrans)
      DIMENSION IUNQ(NQ),NEqATM(NATOMS,NTrans)
      CHARACTER GROUP*4,CHAR*80
      LOGICAL Symflag
c
      character*256 jobname
      Common /job/jobname,LenJ
C
C
C  Try to Open the <sym> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD',ERR=96)
C
C  locate the $symmetry section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:4).NE.'$sym') GO TO 10
c
      READ(40,910) GROUP
      call lowercas(GROUP,4)
      READ(40,920) NATOMS
      READ(40,920) NTrans
      READ(40,920) NDEG
C
C  read rotation matrix
C
      READ(40,900) CHAR
      DO 20 I=1,3
      READ(40,930) (RM(I,J),J=1,3)
 20   CONTINUE
C
      IF(.NOT.Symflag) THEN
C
C  symmetry is NOT being used
C  setup for C1 and exit
C
       GROUP = 'c1  '
       NDEG = MAX(3*NATOMS-6,1)            ! be careful if linear
       NQ = NATOMS
       DO 30 I=1,NATOMS
       IUNQ(I) = I
 30    CONTINUE
       GO TO 95
      ENDIF
C
C  read symmetry-unique atoms
C
      READ(40,940) NQ
      READ(40,950) (IUNQ(I),I=1,NQ)
C
C  read symmetry operations
C
      DO 50 IOP=1,NTrans
      READ(40,900) CHAR
      DO 40 I=1,3
      READ(40,930) (TRANS(I,J,IOP),J=1,3)
 40   CONTINUE
 50   CONTINUE
C
C  read equivalent atoms array
C
      READ(40,900) CHAR
      DO 60 I=1,NTrans
      READ(40,950) (NEqATM(J,I),J=1,NATOMS)
 60   CONTINUE
C
C  close file and exit
C
 95   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C  ..............................................
C    ERROR SECTION
C
 96   CONTINUE
      Call nerror(5,'File IO routine <RdSYM>',
     $      'Unable to find or Unexpected End of <sym> file!',0,0)
c
 900  Format(A80)
 910  Format(8X,A4)
 920  Format(8X,I4)
 930  Format(3F20.14)
 940  Format(19X,I4)
 950  Format(20I4)
c
      END
c =====================================================================
c
      SUBROUTINE RdS2(IUnit,S2,XMult,IsThere)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL IsThere
C
C  reads dipole moment (X, Y, Z components) from <control> file
C
      IsThere = .false.
      call fdcntrl(IUnit,5,'$<s2>',IEnd)
      If(IEnd.EQ.0) Then
       backspace IUnit
       READ(IUnit,1000) S2,XMult
       IsThere = .true.
      EndIf
C
      RETURN
c
 1000 FORMAT(6X,2F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE WrBASIS(NBasis, NShell, NPrim,  NPsp,   NAT,
     $                   INX,    BASDAT, basname, icpsp)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Write basis set data to the <jobname>.basis file
C  If there are pseudopotentials, the atomic center and the number of
C  core electrons will be read from a temporary file
C
C  ARGUMENTS
C
C  NBasis  -  total number of contracted basis functions  (ncf)
C  NShell  -  total number of contracted shells           (ncs)
C  NPrim   -  total number of primitive shells            (nsh)
C  NPsp    -  number of pseudopotentials
C  NAT     -  storage to indicate different atoms
C  INX     -  integer basis set information
C  BASDAT  -  basis exponents and contraction coefficients
C  basname -  basis set name
C  icpsp   -  unit number for (temporary) pseudopotential file
C
C
      DIMENSION INX(12,NShell),BASDAT(13,NPrim),NAT(*)
      CHARACTER*20 basname,cdum
c ..................................................
c -- automatic allocation of arrays in F90
c -- **WARNING**  assume NBasis > NAtoms
      CHARACTER*8 AtSym(NBasis),Symb(NBasis)
c ..................................................
      CHARACTER*8 Atom
      Character*3 TYPE(13),Shell
      character*80 Char
      Logical found
C
      PARAMETER (IUnit=1)
C
      Data TYPE
     $   /'s  ','p  ','l  ','d  ','d6 ','f  ','f10','g15','h21','i28',
     $    'g  ','h  ','i  '/
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open <control> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $    FORM='FORMATTED',STATUS='OLD')
c -- read in number of atoms and number of molecules
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,rdum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,901) Ndum1,Ndum2
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      call rdcntrl(IUnit,5,'$nmol',1,NMol,rdum,cdum)
c -- write basis name
      call wrcntrl(IUnit,6,'$basis',3,idum,dum,basname)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c -- dummy atoms are ignored
      NAtoms = NAtoms-Ndum1-Ndum2
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $   FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,12,'$coordinates',0,idum,rdum,cdum)
c
      IF(NMol.EQ.1) THEN
       DO 10 IAtm=1,NAtoms
       READ(IUnit,900) AtSym(IAtm)
 10    CONTINUE
      ELSE
       IAtm = 0
 20    READ(IUnit,910) Char
       If(Char(1:1).EQ.'$') Then
        If(Char(2:4).EQ.'end') GO TO 25
c -- skip dummy atoms
       Else If( (Char(1:1).eq.'x'.and.Char(2:2).ne.'e') .OR.
     $          (Char(1:1).eq.'d'.and.Char(2:2).eq.'u') .OR.
     $           Char(1:1).eq.'q' ) Then
        GO TO 20
       Else
        IAtm = IAtm+1
        READ(Char,900) AtSym(IAtm)
       EndIf
       GO TO 20
      ENDIF
c
 25   CONTINUE
C
C  massage atomic symbols to remove numbers and
C  catch special symbols
C
      CALL SymbolM(NAtoms,AtSym,Symb)
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  determine the number of different atoms according to atomic
C  symbols, and list their locations in NAT array
C
      CALL GetDiffAtoms(NAtoms,Symb,NDiff,NAT)
C
C  inquire regarding existence of basis file
C  if already exists, make copy as basis2
C
      INQUIRE(FILE=jobname(1:lenJ)//'.basis',EXIST=found)
      If(found) CALL CopyFile(jobname(1:lenJ)//'.basis',lenJ+6,
     $      jobname(1:lenJ)//'.basis2',lenJ+7)
C
C  open <basis> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $    FORM='FORMATTED',STATUS='UNKNOWN')
C
C  start writing basis set data
C
      WRITE(IUnit,1000) NBasis
      WRITE(IUnit,1100) NShell
      WRITE(IUnit,1200) NPrim
c
c ...................................................................
c   PSEUDOPOTENTIAL STUFF HERE
c ...................................................................
      IF(NPsp.GT.0) THEN
C
C  read which centers involve pseudopotentials and the number
C  of core electrons from the (temporary) <psp> file
C
        REWIND icpsp
        READ (icpsp,*)     ! skip 3 lines
        READ (icpsp,*)
        READ (icpsp,*)
c
        WRITE(IUnit,1250) NPsp
c
        DO 30 I=1,NPsp
        READ(icpsp,*) idum,NCore,IAtm
        WRITE(IUnit,1260) IAtm,NCore
 30     CONTINUE
      ENDIF
c ..................................................................
c
      WRITE(IUnit,'(a)') '$basis'
C
C  Loop over the different atoms, extract the basis and write
C  in user-interpretable form
C
      DO 70 I=1,NDiff
      IAtm = NAT(I)
      Atom = Symb(IAtm)
      found = .FALSE.
      DO 60 ics=1,NShell
      LAtm = INX(2,ics)
      IF(LAtm.EQ.IAtm) THEN
        If(.NOT.found) WRITE(IUnit,1300) Atom
        found = .TRUE.
        IType = INX(12,ics)          ! shell type
        Shell = TYPE(IType)          ! shell
        istart = INX(1,ics)+1        ! start of contraction
        iend = INX(5,ics)            ! end of contraction
        ngr = INX(4,ics)             ! no. of general contractions
        IF(IType.EQ.3) THEN          ! l shell
          DO 40 L=istart,iend
          If(L.EQ.istart) Then
           WRITE(IUnit,1400) Shell,BASDAT(1,L),BASDAT(2,L),BASDAT(3,L)
          Else
           WRITE(IUnit,1410) BASDAT(1,L),BASDAT(2,L),BASDAT(3,L)
          EndIf
 40       CONTINUE
        ELSE
          DO 50 L=istart,iend
          If(L.EQ.istart) Then
            WRITE(IUnit,1400) Shell,BASDAT(1,L),
     $                        (BASDAT(igr+2,L),igr=0,ngr)
          Else
            WRITE(IUnit,1410) BASDAT(1,L),
     $                        (BASDAT(igr+2,L),igr=0,ngr)
          EndIf
 50       CONTINUE
        ENDIF
      ENDIF
 60   CONTINUE
 70   CONTINUE
c
      WRITE(IUnit,'(a)') '$end'
      CLOSE (UNIT=IUnit,STATUS='KEEP')
      RETURN
c
  900 Format(A8)
  901 Format(9X,I4,2X,I4)
  910 Format(A80)
 1000 FORMAT('$nbasis',1X,I6)
 1100 FORMAT('$nshell',1X,I6)
 1200 FORMAT('$nprim',2X,I6)
 1250 FORMAT('$npseud',1X,I6)
 1260 FORMAT(2(1X,I4))
 1300 FORMAT('for       ',A8)
 1400 FORMAT(A3,3X,F20.7,9(2X,F10.7))
 1410 FORMAT(6X,F20.7,9(2X,F10.7))
c
      END
c =====================================================================
c
      SUBROUTINE GetDiffAtoms(NAtoms,AtSymb,NDiff,NAT)
      IMPLICIT INTEGER(A-Z)
C
C  how many different atoms are there?
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  NDiff   -  on exit number of different atoms
C  NAT     -  list of first different atoms
C
C
      DIMENSION NAT(NAtoms)
      CHARACTER*8 AtSymb(NAtoms),Symb
C
C
      CALL IZeroIT(NAT,NAtoms)
      NDiff = 0
c
      DO 20 I=1,NAtoms
      Symb = AtSymb(I)
      IF(NAT(I).EQ.0) THEN
        NDiff = NDiff+1
        NAT(I) = NDiff
        DO 10 J=I+1,NAtoms
        If(AtSymb(J).EQ.Symb) NAT(J) = NDiff
 10     CONTINUE
      ENDIF
 20   CONTINUE
c
      DO 40 I=1,NDiff
      DO 30 J=I,NAtoms
      If(NAT(J).EQ.I) Then
       NAT(I) = J
       GO TO 40
      EndIf
 30   CONTINUE
 40   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE wrcntrl(IUnit,Length,Key,IType,INum,RNum,CNum)
      IMPLICIT INTEGER(A-Z)
C
C  Writes keyword "Key" (optionally with value) to <control> file.
C  If keyword already exists, replaces its value
C
C  ARGUMENTS
C
C  IUnit   -  unit number of <control> file (assumed open)
C  Length  -  length of string
C  Key     -  string to be written/replaced
C  IType   -  control of write
C              0 - write keyword
C              1 - write keyword + integer value
C              2 - write keyword + real value
C              3 - write keyword + character string
C             -1 - remove keyword entirely
C  INum    -  integer value to write     (if IType=1)
C  RNum    -  real value to write        (if IType=2)
C  CNum    -  character string to write  (if IType=3)
C
C
      REAL*8 RNum
      CHARACTER*(*) Key
      CHARACTER CNum*20,Char*80
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
      CALL fdcntrl(IUnit,Length,Key,IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,Length,Key,IEnd)
        backspace IUnit
        If(IType.EQ.-1) Then
c  -- do nothing (keyword is to be eliminated)
        Else If(IType.EQ.0) Then
         WRITE(IUnit,'(a)') Key
        Else If(IType.EQ.1) Then
         WRITE(Char,'(I6)') INum
         Char = Key//' '//Char(1:6)
         WRITE(IUnit,'(a)') Char(1:Length+7)
        Else If(IType.EQ.2) Then
         WRITE(Char,'(F20.10)') RNum
         Char = Key//' '//Char(1:20)
         WRITE(IUnit,'(a)') Char(1:Length+21)
        Else
         Char = Key//'  '//CNum
         call rmblan(Char,80,Len)
         WRITE(IUnit,'(a)') Char(1:Len)
        EndIf
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write keyword direct to <control> file
C
        backspace IUnit
        If(IType.EQ.0) Then
         WRITE(IUnit,'(a)') Key
        Else If(IType.EQ.1) Then
         WRITE(Char,'(I6)') INum
         Char = Key//' '//Char(1:6)
         WRITE(IUnit,'(a)') Char(1:Length+7)
        Else If(IType.EQ.2) Then
         WRITE(Char,'(F20.10)') RNum
         Char = Key//' '//Char(1:20)
         WRITE(IUnit,'(a)') Char(1:Length+21)
        Else
         Char = Key//'  '//CNum
         call rmblan(Char,80,Len)
         WRITE(IUnit,'(a)') Char(1:Len)
        EndIf
c
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
c
      END
c =====================================================================
c
      SUBROUTINE WrCoord(NAtoms, AtSymb, XC,     NMol,   IMOL,
     $                   XCharg, XMass)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Routine to write new coordinates to <coord> file
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  current geometry (Cartesian coordinates)
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C  IMOL    -  pointers to start/end of molecules in XC array
C  XCharg  -  atomic charges (can be user-defined)
C  XMass   -  atomic masses
C
C
      REAL*8 XC(3,NAtoms),XCharg(NAtoms),XMass(NAtoms)
      DIMENSION IMOL(NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open the <coord> file
C  (note: a new file will be created if one does not exist)
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='UNKNOWN')
C
C  write the new coordinates
C
      WRITE(40,'(a)') '$coordinates'
c
      IF(NMol.EQ.1) THEN
cc
       DO 10 IAtm=1,NAtoms
       WRITE(40,910) AtSymb(IAtm),XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),
     $               XCharg(IAtm),XMass(IAtm)
 10    CONTINUE
cc
      ELSE
cc
       DO 30 I=1,NMol
       IStrt = IMOL(I)+1
       IEnd = IMOL(I+1)
       DO 20 IAtm=IStrt,IEnd
       WRITE(40,910) AtSymb(IAtm),XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),
     $               XCharg(IAtm),XMass(IAtm)
 20    CONTINUE
       If(IEnd.LT.NAtoms) WRITE(40,'(a)') '$molecule'
 30    CONTINUE
cc
      ENDIF
c
      WRITE(40,'(a)') '$end'
c
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
  910 Format(A8,2X,5F20.12)
c
      END
c ======================================================================
c
      SUBROUTINE WrDeriv(NAt3,DipD,PolD,Dipole,Polar)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  write Cartesian dipole and polarizability derivatives
C  to <deriv> file
C
C  ARGUMENTS
C
C  NAt3    -  3 x number of atoms
C  DipD    -  dipole derivatives in Cartesian coordinates
C  PolD    -  polarizability derivatives in Cartesian coordinates
C  Dipole  -  logical flag
C             set to .TRUE. if dipole derivatives to be written
C  Polar   -  logical flag
C             set to .TRUE. if polarizability derivatives to be written
C
C
      DIMENSION DipD(NAt3,3),PolD(NAt3,6)
      character*80 CHAR
      LOGICAL Dipole,Polar
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.deriv',
     $    FORM='FORMATTED',STATUS='UNKNOWN')
c
      If(Dipole) Then
C
C  write dipole derivatives
C
       WRITE(40,'(a)') '$dipole derivatives'
       DO 20 I=1,NAt3
       WRITE(40,910) DipD(I,1),DipD(I,2),DipD(I,3)
 20    CONTINUE
      EndIf
C
C  make sure any existing dipole derivatives are NOT
C  overwritten
C
      CLOSE (UNIT=40,STATUS='KEEP')
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.deriv',
     $    FORM='FORMATTED',STATUS='OLD',ACCESS='APPEND')
c
      If(Polar) Then
C
C  write polarizability derivatives
C
       WRITE(40,'(a)') '$polarizability derivatives'
       DO 40 I=1,NAt3
       WRITE(40,910) PolD(I,1),PolD(I,2),PolD(I,3)
       WRITE(40,910) PolD(I,4),PolD(I,5),PolD(I,6)
 40    CONTINUE
      EndIf

C  close file and exit
C
      WRITE(40,'(a)') '$end'
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
c
  910 Format(10X,3F20.14)
c
      END
c ====================================================================
c
      SUBROUTINE WrDISP(IUnit,noabc,tz,VALUES)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VALUES(5)
      Logical noabc,tz
      Character*80 Char
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  writes dispersion parameters to control file
C
      iabc = 0
      itz = 0
      If(noabc) iabc=1
      If(tz) itz=1
c
      CALL fdcntrl(IUnit,10,'$paramdisp',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        If(Char(1:6).EQ.'$noabc') GO TO 10    ! skip this line as rewritten here
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,10,'$paramdisp',IEnd)
        backspace IUnit
        WRITE(IUnit,910) (VALUES(J),J=1,5)
        WRITE(IUnit,920) iabc,itz
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write dipole direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) (VALUES(J),J=1,5)
        WRITE(IUnit,920) iabc,itz
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$paramdisp',5(2X,F10.5))
  920 Format('$noabc ',I2,4X,'$tz ',I2)
c
      END
c ====================================================================
c
      SUBROUTINE WrDIP(IUnit,Dip)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Dip(3)
      Character*80 Char
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  writes dipole moment (X, Y, Z components) to <control> file
C
      CALL fdcntrl(IUnit,7,'$dipole',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,7,'$dipole',IEnd)
        backspace IUnit
        WRITE(IUnit,910) Dip(1),Dip(2),Dip(3)
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write dipole direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) Dip(1),Dip(2),Dip(3)
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$dipole',3F20.14)
c
      END
c =====================================================================
c
      SUBROUTINE WrField(IUnit,Field)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Field(3)
      Character*80 Char
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  writes field components (X, Y, Z components) to <control> file
C
      CALL fdcntrl(IUnit,6,'$field',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,6,'$field',IEnd)
        backspace IUnit
        WRITE(IUnit,910) Field(1),Field(2),Field(3)
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write field direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) Field(1),Field(2),Field(3)
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$field',2X,3F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE WrGEOM(NAtoms, NDum1,  Ndum2,  XNuc,   NMol,
     $                  IMOL,   XMass,  Charg,  IMult,  NAlpha,
     $                  NBeta,  IPRNT,  rthrsh, sthrsh, adjst,
     $                  com,    efld,   field)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Write geometry data (including charge, multiplicity
C  and symmetry threshold) to the <control> file
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  NDum1   -  number of charged dummy atoms
C  NDum2   -  number of uncharged dummy atoms
C  XNuc    -  Cartesian coordinates      (old Texas format)
C             (also includes atomic symbols and numbers)
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C  IMOL    -  pointers to start/end of molecules in XNuc array
C  XMass   -  atomic masses
C  Charg   -  total molecular charge
C  IMult   -  multiplicity
C  NAlpha  -  number of occupied closed-shell/alpha MOs
C  NBeta   -  number of occupied beta MOs
C             (may change later if UHF singlet)
C  IPRNT   -  print flag
C  rthrsh  -  threshold for printing interatomic distances
C  sthrsh  -  threshold used to determine symmetry
C  adjst   -  logical flag for reorienting axes
C  com     -  logical flag for putting in "centre of mass" frame
C  efld    -  logical flag for external electric field
C  field   -  field components
C
C
      DIMENSION XNuc(5,NAtoms),XMass(NAtoms),IMOL(NAtoms),field(3)
      LOGICAL adjst,com,efld
      REAL*8 rdum
      Character*20 cdum
C
      PARAMETER (IUnit=1, Zero=0.0d0)
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  Two Operational Modes
C   (1) <control> file already exists
C   (2) new <control> file
C  In the second case we can simply write the data we've got;
C  in the first we have to overwrite any new data, preserving
C  the existing data
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD',ERR=50)
C
C -- (1) <control> file extant
C  replace existing data
C
      CALL wrcntrl(IUnit,7,'$natoms',1,NAtoms,rdum,cdum)
      If(NDum1+NDum2.GT.0) Then
        WRITE(IUnit,1050) NDum1,NDum2
      Else
        CALL fdcntrl(IUnit,7,'$ndummy',IEnd)
        If(IEnd.EQ.0)
     $    CALL wrcntrl(IUnit,7,'$ndummy',-1,idum,rdum,cdum)
      EndIf
      CALL wrcntrl(IUnit,5,'$nmol',1,NMol,rdum,cdum)
      CALL wrcntrl(IUnit,7,'$charge',2,Idum,Charg,cdum)
      CALL wrcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      CALL wrcntrl(IUnit,7,'$nalpha',1,NAlpha,rdum,cdum)
      CALL wrcntrl(IUnit,7,'$nbeta ',1,NBeta,rdum,cdum)
      CALL wrcntrl(IUnit,14,'$sym_threshold',2,idum,sthrsh,cdum)
      If(.not.adjst) Then
        CALL wrcntrl(IUnit,9,'$noorient',0,idum,rdum,cdum)
      Else
        CALL fdcntrl(IUnit,9,'$noorient',IEnd)
        If(IEnd.EQ.0)
     $    CALL wrcntrl(IUnit,9,'$noorient',-1,idum,rdum,cdum)
      EndIf
      If(.not.com) Then
        CALL wrcntrl(IUnit,6,'$nocms',0,idum,rdum,cdum)
      Else
        CALL fdcntrl(IUnit,6,'$nocms',IEnd)
        If(IEnd.EQ.0)
     $    CALL wrcntrl(IUnit,6,'$nocms',-1,idum,rdum,cdum)
      EndIf
      If(rthrsh.GT.Zero) Then
        CALL wrcntrl(IUnit,18,'$dist_print_thresh',2,idum,rthrsh,cdum)
      Else
        CALL fdcntrl(IUnit,18,'$dist_print_thresh',IEnd)
        If(IEnd.EQ.0)
     $    CALL wrcntrl(IUnit,18,'$dist_print_thresh',-1,idum,rdum,cdum)
      EndIf
      CALL wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      If(efld) Then
        CALL WrField(IUnit,field)
      Else
        CALL fdcntrl(IUnit,4,'$pol',IEnd)
c
c -- do NOT remove field if polarizability option is set
        If(IEnd.NE.0) Then
          CALL fdcntrl(IUnit,6,'$field',IEnd)
          If(IEnd.EQ.0)
     $      CALL wrcntrl(IUnit,6,'$field',-1,idum,rdum,cdum)
        EndIf
      EndIf
c
c  write type of executable ( 32 or 64-bit ) and intsize
c
      nbits = 0
      call tstival('nbits',ifound)
      if(ifound.ne.0) call getival('nbits',nbits)
      CALL wrcntrl(IUnit,6,'$nbits',1,nbits,rdum,cdum)
      intsize = 0
      call tstival('ints',ifound)
      if(ifound.ne.0) call getival('ints',intsize)
      CALL wrcntrl(IUnit,8,'$intsize',1,intsize,rdum,cdum)

      GO TO 95
 50   CONTINUE
C
C -- (2) new <control> file
C  this is the first block of data written to the <control> file
C  start writing from the top
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='NEW')
c
      WRITE(IUnit,1000) NAtoms
      If(NDum1+NDum2.GT.0) WRITE(IUnit,1050) NDum1,NDum2
      WRITE(IUnit,1150) NMol
      WRITE(IUnit,1200) Charg
      WRITE(IUnit,1300) IMult
      WRITE(IUnit,1400) NAlpha
      WRITE(IUnit,1500) NBeta
      WRITE(IUnit,1600) sthrsh
      If(.not.adjst) WRITE(IUnit,'(a)') '$noorient'
      If(.not.com) WRITE(IUNIT,'(a)') '$nocms'
      WRITE(IUnit,1700) rthrsh
      WRITE(IUnit,1800) IPRNT
      If(efld) WRITE(IUnit,1900) field
c
c  write type of executable ( 32 or 64-bit ) and intsize
c
      nbits = 0
      call tstival('nbits',ifound)
      if(ifound.ne.0) call getival('nbits',nbits)
      WRITE(IUnit,'(''$nbits'',2x,i4)') nbits
      intsize = 0
      call tstival('ints',ifound)
      if(ifound.ne.0) call getival('ints',intsize)
      WRITE(IUnit,'(''$intsize'',2x,i4)') intsize

      WRITE(IUnit,'(a)') '$end'
c
 95   CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  write the coordinates to the <coord> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(IUnit,'(a)') '$coordinates'
c
      IF(NMol.EQ.1) THEN
       DO 10 I=1,NAtoms
       WRITE(IUnit,1100) XNuc(5,I),XNuc(2,I),XNuc(3,I),XNuc(4,I),
     $                   XNuc(1,I),XMass(I)
 10    CONTINUE
      ELSE
       DO 30 J=1,NMol
       IStrt = IMOL(J)+1
       IEnd = IMOL(J+1)
       DO 20 I=IStrt,IEnd
       WRITE(IUnit,1100) XNuc(5,I),XNuc(2,I),XNuc(3,I),XNuc(4,I),
     $                   XNuc(1,I),XMass(I)
 20    CONTINUE
       If(IEnd.LT.NAtoms) WRITE(IUnit,'(a)') '$molecule'
 30    CONTINUE
      ENDIF
c
      WRITE(IUnit,'(a)') '$end'
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
      RETURN
c
 1000 FORMAT('$natoms',2X,I4)
 1050 FORMAT('$ndummy',2X,I4,2X,I4)
 1100 FORMAT(A8,2X,5F20.12)
 1150 FORMAT('$nmol',1X,I4)
 1200 FORMAT('$charge',1X,F12.6)
 1300 FORMAT('$multiplicity',1X,I4)
 1400 FORMAT('$nalpha',1X,I4)
 1500 FORMAT('$nbeta ',1X,I4)
 1600 FORMAT('$sym_threshold',1X,F9.6)
 1700 FORMAT('$dist_print_thresh',1X,F9.6)
 1800 FORMAT('$print',2X,I2)
 1900 FORMAT('$field',2X,3F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE WrGRAD(NATOMS,GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Writes gradient to <grad> file
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  GC      -  Cartesian gradient
C
C
      REAL*8 GC(3,NATOMS)
      Character*80 CHAR
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open GRAD file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.grad',
     $     FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE(40,'(a)') '$gradient'
c
      DO 10 IAtm=1,NATOMS
      WRITE(40,910) GC(1,IAtm),GC(2,IAtm),GC(3,IAtm)
 10   CONTINUE
      WRITE(40,'(a)') '$end'
c
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
c
 910  Format(10X,3F20.14)
c
      END
c =====================================================================
c
      SUBROUTINE WrGUESS(IUnit,  GUESS,  uhfs,   iop,    angle)
      IMPLICIT INTEGER(A-Z)
C
C  write SCF guess options to <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number of <control> file  (should already be open)
C  GUESS   -  main guess optin (character string)
C  uhfs    -  integer flag indicating UHF singlet
C              1 = yes   0 = no
C  iop     -  array of integer options
C             relevant entries are as follows:
C              iop(1,2)  - alpha orbital to swap, if any (occupied)
C              iop(2,2)  -  ditto  corresponding orbital (virtual)
C              iop(1,3)  -  beta orbital to swap, if any (occupied)
C              iop(2,3)  -  ditto  corresponding orbital (virtual)
C              iop(1,4)  - alpha orbital involved in orbital mixing
C              iop(2,4)  -  ditto  beta orbital
C              iop(1,9)  - first alpha orbital in multiple swap
C              iop(2,9)  -  ditto  corresponding orbital
C              iop(1,10) - second alpha orbital in multiple swap
C              iop(2,10) -  ditto  corresponding orbital
C              iop(1,11) - third alpha orbital in multiple swap
C              iop(2,11) -  ditto  corresponding orbital
C              iop(1,12) - first beta orbital in multiple swap
C              iop(2,12) -  ditto  corresponding orbital
C              iop(1,13) - second beta orbital in multiple swap
C              iop(2,13) -  ditto  corresponding orbital
C              iop(1,14) - third beta orbital in multiple swap
C              iop(2,14) -  ditto  corresponding orbital
C  angle   -  rotation angle for orbital mixing
C
C
      REAL*8 angle,rdum
      INTEGER iop(3,*)
      CHARACTER GUESS*80,cdum*20
C
C
C  assign variables from iop array before we forget
C
      swap1 = iop(1,2)
      swap2 = iop(2,2)
      swab1 = iop(1,3)
      swab2 = iop(2,3)
      mix1  = iop(1,4)
      mix2  = iop(2,4)
c -- multiple swaps
      swa1  = iop(1,9)
      swa2  = iop(2,9)
      swa3  = iop(1,10)
      swa4  = iop(2,10)
      swa5  = iop(1,11)
      swa6  = iop(2,11)
      swb1  = iop(1,12)
      swb2  = iop(2,12)
      swb3  = iop(1,13)
      swb4  = iop(2,13)
      swb5  = iop(1,14)
      swb6  = iop(2,14)
C
C  locate $guess line OR end-of-file ($end)
C
      call fdcntrl(IUnit,6,'$guess',IEnd)
C
C  if found $guess then remove it
C  (always write guess options at end of <control> file)
C
      If(IEnd.EQ.0) Then
       call wrcntrl(IUnit,6,'$guess',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,12,'$uhf_singlet',IEnd)
       If(IEnd.EQ.0)
     $    call wrcntrl(IUnit,12,'$uhf_singlet',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,5,'$swap',IEnd)
       If(IEnd.EQ.0) call wrcntrl(IUnit,5,'$swap',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,10,'$beta_swap',IEnd)
       If(IEnd.EQ.0)
     $    call wrcntrl(IUnit,10,'$beta_swap',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,4,'$mix',IEnd)
       If(IEnd.EQ.0) call wrcntrl(IUnit,4,'$mix',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,6,'$guess',IEnd)
c -- multiple swaps
       call fdcntrl(IUnit,5,'$swp1',IEnd)
       If(IEnd.EQ.0) call wrcntrl(IUnit,5,'$swp1',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,5,'$swp2',IEnd)
       If(IEnd.EQ.0) call wrcntrl(IUnit,5,'$swp2',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,5,'$swp3',IEnd)
       If(IEnd.EQ.0) call wrcntrl(IUnit,5,'$swp3',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,10,'$beta_swp1',IEnd)
       If(IEnd.EQ.0)
     $    call wrcntrl(IUnit,10,'$beta_swp1',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,10,'$beta_swp2',IEnd)
       If(IEnd.EQ.0)
     $    call wrcntrl(IUnit,10,'$beta_swp2',-1,idum,rdum,cdum)
       call fdcntrl(IUnit,10,'$beta_swp3',IEnd)
       If(IEnd.EQ.0)
     $    call wrcntrl(IUnit,10,'$beta_swp3',-1,idum,rdum,cdum)
      EndIf
c
      backspace IUnit
C
C  now write guess option
C
      CALL RmBlank(80,GUESS,Len)
      call lowerca2(GUESS,Len)
      GUESS = '$guess  '//GUESS(1:Len)
      WRITE(IUnit,'(a)') GUESS(1:Len+8)
C
C uhf singlet?
C
      If(uhfs.EQ.1) WRITE(IUnit,'(a)') '$uhf_singlet'
C
C  write swap option if MOs are to be swapped
C
      If(swap1.NE.0) WRITE(IUnit,1000) swap1,swap2
      If(swab1.NE.0) WRITE(IUnit,1100) swab1,swab2
      If(mix1.NE.0)  WRITE(IUnit,1200) mix1,mix2,angle
c -- multiple swaps (actually up to 3) now a possibity
      If(swa1.NE.0)  WRITE(IUnit,1301) swa1,swa2
      If(swa3.NE.0)  WRITE(IUnit,1302) swa3,swa4
      If(swa5.NE.0)  WRITE(IUnit,1303) swa5,swa6
      If(swb1.NE.0)  WRITE(IUnit,1304) swb1,swb2
      If(swb3.NE.0)  WRITE(IUnit,1305) swb3,swb4
      If(swb5.NE.0)  WRITE(IUnit,1306) swb5,swb6
c
      WRITE(IUnit,'(a)') '$end'
c .............................................................
c -- Special: Swap/Mix option currently unavailable for GUESS=CORE
c
      If((swap1.NE.0.OR.swab1.NE.0).AND.GUESS(9:12).EQ.'core') Then
        Call nerror(9,'File IO routine <WrGUESS>',
     $       'Swap option currently unavailable for GUESS=CORE',0,0)
      EndIf
      If(mix1.NE.0.AND.GUESS(9:12).EQ.'core') Then
        Call nerror(10,'File IO routine <WrGUESS>',
     $       'Mix option currently unavailable for GUESS=CORE',0,0)
      EndIf
C  ...................................................
C
      RETURN
c
 1000 FORMAT('$swap  ',I4,1X,I4)
 1100 FORMAT('$beta_swap  ',I4,1X,I4)
 1200 FORMAT('$mix   ',I4,1X,I4,2X,F8.4)
 1301 FORMAT('$swp1  ',I4,1X,I4)
 1302 FORMAT('$swp2  ',I4,1X,I4)
 1303 FORMAT('$swp3  ',I4,1X,I4)
 1304 FORMAT('$beta_swp1  ',I4,1X,I4)
 1305 FORMAT('$beta_swp2  ',I4,1X,I4)
 1306 FORMAT('$beta_swp3  ',I4,1X,I4)
c
      END
c =====================================================================
c
      SUBROUTINE WrHESS(HFile,Len,IHflag,NDim,HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Routine to write Hessian to file <HFile>
C
C  ARGUMENTS
C
C  HFile   -  File name
C  Len     -  Length of character string HFile
C  IHflag  -  Hessian type
C              0 - Exact Cartesian Hessian
C              1 - approximate Cartesian Hessian
C              2 - Hessian in internal coordinates
C              3 - primitive Hessian
C              4 - Force Field Hessian
C              5 - SQM scaled Hessian
C  NDim    -  Dimension of Hessian
C  HESS    -  Hessian matrix
C
C
      REAL*8 HESS(NDim,NDim)
      CHARACTER*(*) HFile
C
C
C  Open Hessian file
C
      OPEN (UNIT=40,FILE=HFile,FORM='FORMATTED',STATUS='UNKNOWN')
C
C  write title
C
      IF(IHflag.EQ.0) THEN
       WRITE(40,'(a)') '$hessian'
      ELSE IF(IHflag.EQ.1) THEN
       WRITE(40,'(a)') '$hessian approximate'
      ELSE IF(IHflag.EQ.2) THEN
       WRITE(40,'(a)') '$hessint'
      ELSE IF(IHflag.EQ.3) THEN
       WRITE(40,'(a)') '$hessprim'
      ELSE IF(IHflag.EQ.4) THEN
       WRITE(40,'(a)') '$hessian forcefield'
      ELSE IF(IHflag.EQ.5) THEN
       WRITE(40,'(a)') '$hessian SQM-scaled'
      ELSE
       Call nerror(42,'OPTIMIZE module',
     $    'Unknown Hessian Type in routine <WrHESS>',0,0)
      ENDIF
C
C  write Hessian dimension
C
      WRITE(40,900) NDim
C
C  write the Hessian
C
      DO 10 I=1,NDim
      WRITE(40,*) (HESS(I,J),J=1,I)
 10   CONTINUE
c
      WRITE(40,'(a)') '$end'
c
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
 900  Format('Dimension  ',I5)
c
      END
c =====================================================================
c
      SUBROUTINE WrPATH(IUnit,coord,DMax,ISign,MaxCyc,DTol)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads Reaction Path information from <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to read from
C  coord   -  coordinate type
C              cart - Cartesian coordinates
C              zmat - z-matrix coordinates
C              mwgt - mass-weighted Cartesians
C  DMax    -  maximum step size
C  ISign   -  sign of Hessian mode on first step
C  MaxCyc  -  maximum number of steps
C  DTol    -  displacement criterion for stopping search
C
C
      Character*80 Char
      Character*4 coord
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  locate path keyword
C
      CALL fdcntrl(IUnit,5,'$path',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,5,'$path',IEnd)
        backspace IUnit
        WRITE(IUnit,910) coord,DMax,ISign,MaxCyc,DTol
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write scan data direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) coord,DMax,ISign,MaxCyc,DTol
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$path  ',A4,2X,F9.6,2X,I2,2X,I4,2X,F9.6)
c
      END
c =====================================================================
c
      SUBROUTINE WrPropG(IUnit,factor,r0f,LMax,NRad,NAng)
      IMPLICIT INTEGER(A-Z)
C
C  Writes numerical grid information to <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to write to
C  factor  -  radial grid factor
C  r0f     -  Chipman radial factor
C  LMax    -  maximum angular momentum for spherical harmonics
C  NRad    -  number of radial points
C  NAng    -  number of angular points
C
      REAL*8 factor,r0f
      Character*80 Char
c
      character*256 jobname
      Common /job/jobname,lenJ
C
      CALL fdcntrl(IUnit,10,'$prop_grid',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,10,'$prop_grid',IEnd)
        backspace IUnit
        WRITE(IUnit,910) factor,r0f,LMax,NRad,NAng
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write dipole direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) factor,r0f,LMax,NRad,NAng
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$prop_grid  fac=',f7.4,'  r0f=',f7.4,'  LMax=',I5,2X,2I5)
c
      END
c =====================================================================
c
      SUBROUTINE WrPOL(IUnit,field,IType)
      IMPLICIT INTEGER(A-Z)
C
C  Writes numerical polarizability information to <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number to write to
C  field   -  field strength (atomic units)
C  IType   -  job type
C               0 - calculate polarizabilities
C               1 -  ditto + dipole derivatives
C               2 -  ditto + polarizability derivatives
C
      REAL*8 field
      Character*80 Char
c
      character*256 jobname
      Common /job/jobname,lenJ
C
      CALL fdcntrl(IUnit,4,'$pol',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,4,'$pol',IEnd)
        backspace IUnit
        WRITE(IUnit,910) field,IType
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write dipole direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) field,IType
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$pol   field=',f12.6,'  type=',I2)
c
      END
c =====================================================================
c
      SUBROUTINE WrSCAN(IUnit,  zmat,   Coord,  I,      J,
     $                  K,      L,      NPComp, IPCOMP, PCWght,
     $                  R1,     R2,     Step,   force)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Writes Data for potential scan onto <control> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number for (open) control file
C  zmat    -  logical flag indicating Z-Matrix scan
C  Coord   -  internal coordinate type or Z-Matrix variable to scan
C  I,J,K,L -  atoms defining internal coordinate
C  NPComp  -  number of primitives in composite coordinate
C  IPCOMP  -  constraint type and atoms involved in constraint
C               IPComp(1,I) - constraint type
C               IPComp(2,I) to IPComp(5,I) - atoms in constraint
C  PCWght  -  weight of each primitive in composite constraint
C  R1      -  starting value for coordinate scan
C  R2      -  final value for coordinate scan
C  Step    -  step length for scan
C  force   -  logical flag for external force scan
C
C
      DIMENSION IPCOMP(5,NPComp),PCWght(NPComp)
      Character*80 Char
      Character*4 Coord
      Logical zmat,force
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  locate scan keyword
C
      CALL fdcntrl(IUnit,5,'$scan',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,5,'$scan',IEnd)
        backspace IUnit
        If(zmat) Then
         WRITE(IUnit,910) Coord,R1,R2,Step
        Else If(force) Then
         WRITE(IUnit,930) I,J,R1,R2,Step
        Else If(NPComp.GT.0) Then
         WRITE(IUnit,915) R1,R2,Step
         WRITE(IUnit,'(a)') '#composite'
         DO I=1,NPComp
         WRITE(IUnit,925) (IPCOMP(I,J),J=1,5),PCWght(I)
         EndDO
         WRITE(IUnit,'(a)') '#endcomposite'
        Else
         WRITE(IUnit,920) Coord,I,J,K,L,R1,R2,Step
        EndIf
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write scan data direct to <control> file
C
        backspace IUnit
        If(zmat) Then
         WRITE(IUnit,910) Coord,R1,R2,Step
        Else If(force) Then
         WRITE(IUnit,930) I,J,R1,R2,Step
        Else If(NPComp.GT.0) Then
         WRITE(IUnit,915) R1,R2,Step
         WRITE(IUnit,'(a)') '#composite'
         DO I=1,NPComp
         WRITE(IUnit,925) IPCOMP(1,I),(IPCOMP(J,I),J=2,5),PCWght(I)
         EndDO
         WRITE(IUnit,'(a)') '#endcomposite'
        Else
         WRITE(IUnit,920) Coord,I,J,K,L,R1,R2,Step
        EndIf
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$scan  zmat  ',A4,2X,3F12.6)
  915 Format('$scan  comp  ',2X,3F12.6)
  920 Format('$scan  ',A4,2X,4I5,2X,3F12.6)
  925 Format(1X,I4,2X,4I4,F12.6)
  930 Format('$scan  extf  ',2I4,2X,3F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE WrSYM(NAtoms, RM,     GROUP,  NTrans, NDEG,
     $                 NQ,     IUNQ,   TRANS,  NEqATM)
      IMPLICIT INTEGER(A-Z)
C
C
C  Writes symmetry information to <sym> file
C
C  ARGUMENTS
C
C  NATOMS  -  number of real atoms
C  RM      -  rotation matrix used to reorient to "standard orientation"
C  GROUP   -  point group
C  NTrans  -  number of symmetry operations
C  NDEG    -  number of internal degrees of freedom
C  NQ      -  number of symmetry unique atoms
C  IUNQ    -  list of symmetry unique atoms
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C
      REAL*8 RM(3,3),TRANS(3,3,NTrans)
      DIMENSION IUNQ(NQ),NEqATM(NATOMS,NTrans)
      CHARACTER*4 GROUP
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  open <sym> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.sym',
     $    FORM='FORMATTED',STATUS='UNKNOWN')
C
      WRITE(40,'(a)') '$symmetry'
      WRITE(40,900) '$group  ',GROUP
      WRITE(40,910) '$natoms ',NAtoms
      WRITE(40,910) '$ntrans ',NTrans
      WRITE(40,910) '$ndeg   ',NDEG
C
      WRITE(40,'(a)') '$rotation_matrix:'
      DO 10 I=1,3
      WRITE(40,920) (RM(I,J),J=1,3)
 10   CONTINUE
C
      WRITE(40,930) NQ
      WRITE(40,940) (IUNQ(I),I=1,NQ)
C
      DO 20 IOP=1,NTrans
      WRITE(40,950) IOP
      DO 19 I=1,3
      WRITE(40,920) (TRANS(I,J,IOP),J=1,3)
 19   CONTINUE
 20   CONTINUE
C
      WRITE(40,'(a)') '$equivalent_atoms_array:'
      DO 30 I=1,NTrans
      WRITE(40,940) (NEqATM(J,I),J=1,NAtoms)
 30   CONTINUE
C
      WRITE(40,'(a)') '$end'
C
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
c
  900 FORMAT(A8,A4)
  910 FORMAT(A8,I4)
  920 FORMAT(3F20.14)
  930 FORMAT('$sym_unique_atoms  ',I4)
  940 FORMAT(20I4)
  950 FORMAT('symmetry_operation:  ',I4)
c
      END
c =====================================================================
c
      SUBROUTINE CopyFile(File1,Len1,File2,Len2)
      IMPLICIT INTEGER(A-Z)
C
C  makes a duplicate of file1 as file2
C
C  ARGUMENTS
C
C  File1   -  name of file to be duplicated
C  Len1    -  length of character string "File1"
C  File2   -  name of duplicate file
C  Len2    -  length of character string "File2"
C
C
      CHARACTER File1*(*),File2*(*)
      Character*80 Char
      Character*1 Char1(80)
      Equivalence(Char,Char1(1))
C
C
      OPEN (UNIT=40,File=File1,FORM='FORMATTED',STATUS='OLD')
      OPEN (UNIT=41,File=File2,FORM='FORMATTED',STATUS='UNKNOWN')
c
 10   CONTINUE
      call blankit(Char,80)
      READ(40,'(80a1)',END=95) Char1
      call rmblan(Char,80,Len)
      WRITE(41,910) (Char(I:I),I=1,Len)
      GO TO 10
c
 95   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
      CLOSE (UNIT=41,STATUS='KEEP')
c
      RETURN
c
  910  Format(80A1)
c
      END
c =====================================================================
c
      SUBROUTINE WrS2(IUnit,S2,XMult)
      IMPLICIT REAL*8(A-H,O-Z)
      Character*80 Char
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C  writes <S**2> and multiplicity to <control> file
C
      CALL fdcntrl(IUnit,5,'$<s2>',IEnd)
c
      IF(IEnd.EQ.0) THEN
C
C  keyword already exists on file
C  copy <control> file beyond keyword to unit 41
C
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.tempppxx',
     $        FORM='FORMATTED',STATUS='NEW')
c
 10     CONTINUE
        READ(IUnit,900,End=20) Char
        WRITE(41,900) Char
        GO TO 10
c
 20     CONTINUE
C
C  now replace keyword
C
        CALL fdcntrl(IUnit,5,'$<s2>',IEnd)
        backspace IUnit
        WRITE(IUnit,910) S2,XMult
C
C  copy back rest of file
C
        REWIND 41
 30     CONTINUE
        READ(41,900,END=40) Char
        call rmblan(Char,80,Len)
        WRITE(IUnit,'(a)') Char(1:Len)
        GO TO 30
c
 40     CONTINUE
        CLOSE (UNIT=41,STATUS='DELETE')
      ELSE
C
C  found end-of-file
C  simply write dipole direct to <control> file
C
        backspace IUnit
        WRITE(IUnit,910) S2,XMult
        WRITE(IUnit,'(a)') '$end'
      ENDIF
C
      RETURN
c
  900 Format(A80)
  910 Format('$<s2> ',2F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE AppendFile(IUnit1,IUnit2,NRec)
      IMPLICIT INTEGER(A-Z)
C
C  Appends NRec records from IUnit1 to IUnit2
C  (both files most be opened and positioned for read and write
C   respectively)
C
C  ARGUMENTS
C
C  IUnit1  -  unit to be read
C  IUnit2  -  unit to be written
C  NRec    -  number of records to append
C             (if -1 append all records in IUnit1)
C
C
      CHARACTER*80 Char
C
C
      ICount = 0
 10   CONTINUE
      READ(IUnit1,900,END=95) Char
      call rmblan(Char,80,Len)
      WRITE(IUnit2,910) (Char(I:I),I=1,Len)
      ICount = ICount+1
      If(ICount.NE.NRec) GO TO 10
c
 95   CONTINUE
      If(NRec.GT.0.AND.ICount.LT.NRec)
     $  Call nerror(11,'File IO routine <AppendFile>',
     $   'EOF reached before requested number of records appended',
     $    NRec,0)
c
      RETURN
c
 900  Format(A80)
 910  Format(80A1)
c
      END
c =====================================================================
c
      SUBROUTINE summary(iflag,Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This subroutine writes a full or a partial job summary on the
C  short output (formerly unit icond in old texas code)
C
C  ARGUMENTS
C
C  iflag   -  controls what is output to summary output on this call
C              -2 - just write basis set name or data
C              -1 - same as 0 below, but SCF did not converge
C               0 - IMPLIES JOB IS A SINGLE POINT ENERGY + PROPERTIES
C                   write basis set data, energy, dipole and <S**2>
C              >0 - IMPLIES JOB IS AN OPTIMIZATION
C                   write final geometry and any molecular properties
C  na      -  number of atoms
C             (for automatic allocation of character array AtSymb)
C  Z       -  general working storage
C
C
      DIMENSION Z(*),Dip(3)
      CHARACTER Char*80,cdum*20,dftxc(27)*6,wvfnc*20,
     $          AtSymb*8
      Logical Lflag
      PARAMETER (IUnit=1,icond=8)
      PARAMETER (Debye=0.39342658d0)
      data dftxc/'hfs   ','svwn  ','svwn5 ','hfb   ','bvwn  ','bvwn5 ',
     1           'bp86  ','bpw91 ','blyp  ','bvp86 ','optx  ','ovwn  ',
     2           'ovwn5 ','op86  ','opw91 ','olyp  ','pw91  ','pbe   ',
     3           'o3lyp ','b3lyp ','b3pw91','wah   ','user  ','b97   ',
     4           'b97-1 ','b97-2 ','hcth  '/
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open the <control file>
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  always read the wavefunction type
C
      call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
C  ...................................................................
C
      IF(iflag.EQ.-2) THEN
C
C  -----------------------------------
C    basis set summary
C  -----------------------------------
C
        WRITE(icond,*)
c
        call rdcntrl(IUnit,7,'$charge',2,idum,Chrge,cdum)
        call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
        call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,rdum,cdum)
        call rdcntrl(IUnit,6,'$nbeta',1,NBeta,rdum,cdum)
        CALL RdDFT(IUnit,idft)
c
        WRITE(icond,1000) Chrge,IMult
        If(idft.EQ.0) WRITE(icond,1100) wvfnc
        If(idft.GT.0) WRITE(icond,1200) wvfnc,dftxc(idft)
c
        call rdcntrl(IUnit,6,'$basis',3,idum,rdum,cdum)
        call RmBlank(20,cdum,len)
c
c -- write basis set name and size
        Char = '  Basis set: '//cdum(1:len)
        WRITE(icond,'(a)') Char(1:len+13)
c
        call tstival('ncf',iexist)
        If(iexist.eq.1) Then
          ncf=igetival('ncf')
          WRITE(icond,1250) ncf
        EndIf
c
        If(cdum(1:11).EQ.'nonstandard') Then
c -- print out basis as appears on <basis> file
          OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.basis',
     $          FORM='FORMATTED',STATUS='OLD')
          CALL rdcntrl(40,6,'$basis',0,idum,dum,cdum)
 10       READ(40,'(a80)') Char
          If(Char(1:4).NE.'$end') Then
           call rmblan(Char,80,Len)
           WRITE(icond,'(a)') '  '//Char(1:Len)
           GO TO 10
          Else
           CLOSE(UNIT=40,STATUS='KEEP')
          EndIf
        EndIf
c
        WRITE(icond,*)
c
      ELSE IF(iflag.LE.0) THEN
C
C  --------------------------------------
C    single-point energy
C  --------------------------------------
C
        call fdcntrl(IUnit,7,'$energy',IEnd)
c
        If(IEnd.EQ.0) Then
          If(wvfnc(2:3).NE.'HF'.AND.wvfnc(2:3).NE.'DF'.AND.
     $       wvfnc(1:5).NE.'Force') Then
            call rdcntrl(IUnit,7,'$energy',2,idum,E,cdum)
            call rdcntrl(IUnit,5,'$escf',2,idum,ESCF,cdum)
            WRITE(icond,1300) wvfnc,E,ESCF
          Else
            call rdcntrl(IUnit,7,'$energy',2,idum,E,cdum)
            WRITE(icond,1310) E
          EndIf
        EndIf
c
        If(iflag.EQ.-1)
     $     WRITE(icond,*) ' *** WARNING  SCF DID NOT CONVERGE ***'
c
      ELSE IF(iflag.GT.0) THEN
C
C  --------------------------------------
C    final geometry printout
C  --------------------------------------
C
        call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.coord',
     $        FORM='FORMATTED',STATUS='OLD')
C
        WRITE(icond,1400)
C
C  locate the $coordinate section
C
        CALL rdcntrl(40,12,'$coordinates',0,idum,dum,char)
C
C  now read coordinates
C
 20     READ(40,'(a80)') Char
        If(Char(1:4).NE.'$end') Then
         call rmblan(Char,70,Len)
         If(Char(1:1).NE.'$') Then
          READ(Char,910) AtSymb,XX,YY,ZZ
          XX = XX/ANTOAU
          YY = YY/ANTOAU
          ZZ = ZZ/ANTOAU
          WRITE(icond,910) AtSymb,XX,YY,ZZ
         Else
          WRITE(icond,'(a)') Char(1:Len)
         EndIf
         GO TO 20
        EndIf
        CLOSE(UNIT=40,STATUS='KEEP')
        WRITE(icond,*)
      ENDIF
C
C  --------------------------------------
C    dipole moment and <S**2>
C  --------------------------------------
C
      IF(iflag.GE.0.AND.
     $    (wvfnc(2:3).EQ.'HF'.OR.wvfnc(2:3).EQ.'DF')) THEN
c
        call RdDIP(IUnit,Dip,Lflag)
        If(Lflag) Then
          Dip(1) = Dip(1)/debye
          Dip(2) = Dip(2)/debye
          Dip(3) = Dip(3)/debye
          diptot=sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
          WRITE(icond,1500) Dip(1),Dip(2),Dip(3),diptot
        EndIf
c
        call RdS2(IUnit,s2,xmult,Lflag)
        If(Lflag) Write(icond,1600) s2,xmult
        WRITE(icond,*)
      ENDIF
C
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
 95   CONTINUE
      RETURN
c
  910 Format(A8,2X,3F20.12)
 1000 FORMAT('  Charge: ',F9.6,'  Multiplicity: ',I3)
 1100 FORMAT('  Wavefunction:  ',A20)
 1200 FORMAT('  Wavefunction:  ',A20,'    XC potential: ',A6)
 1250 FORMAT('  Number of contracted basis functions: ',I5)
 1300 FORMAT(2X,A8,' Energy is: ',F20.9,/,
     $        '  SCF      Energy is: ',F20.9)
 1310 FORMAT('  Energy is: ',F20.9,' au',/)
 1400 FORMAT(/,25X,' CONVERGED GEOMETRY',/,
     $         22X,' Coordinates (Angstroms)',/,
     $         22X,'X',19X,'Y',19X,'Z')
 1500 FORMAT('  dipole/D = ',2X,3F10.6,'  total=',F10.6)
 1600 FORMAT('  Expectation value of S**2:',F11.7,'  Multiplicity:',
     $          F11.7)
c
      END
c======================================================================
      SUBROUTINE GetCoord(NAtoms,nsh,ncs,INX,XNuc,BASDAT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Gets Cartesian coordinates from <coord> file and puts them
C  in XNuc and BASDAT arrays for compatability with SCF module
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  nsh     -  number of primitive shells
C  ncs     -  number of contracted shells
C  INX     -  integer basis function data
C  XNuc    -  on exit contains current nuclear data
C  BASDAT  -  basis function exponents and coefficients
C             (updated here with current nuclear coordinates - obsolescent)
C
C
      Dimension INX(12,ncs),XNuc(5,NAtoms),BASDAT(13,nsh)
c .........................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
      DIMENSION XC(3,NAtoms),IMOL(NAtoms)
c .........................................................
      Character Char*20
      PARAMETER (IUnit=1)
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open <coord> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoord(IUnit,NAtoms,AtSymb,XC,-1,jnk)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  update XNuc array
C
      DO 10 IAtm=1,NAtoms
      XNuc(2,IAtm) = XC(1,IAtm)
      XNuc(3,IAtm) = XC(2,IAtm)
      XNuc(4,IAtm) = XC(3,IAtm)
 10   CONTINUE
C
C  now try and sort out basis crap for that Wolinski guy  (sigh!)
C
      DO 30 ics=1,ncs
      IAtm = INX(2,ics)
      X = XNuc(2,IAtm)
      Y = XNuc(3,IAtm)
      Z = XNuc(4,IAtm)
      istart = INX(1,ics)+1
      iend = INX(5,ics)
      DO 20 I=istart,iend
      BASDAT(11,I) = X
      BASDAT(12,I) = Y
      BASDAT(13,I) = Z
 20   CONTINUE
 30   CONTINUE
C
C  actually, that wasn't as bad as I first thought ....
C
      RETURN
c
 910  Format(10X,4F20.14)
c
      END
c ======================================================================
c
      SUBROUTINE WrAAT(NAt,AAT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  write atomic axial tensors
C  to <deriv> file
C
C  ARGUMENTS
C
C  NAt    -   number of atoms
C  AAT    -   atomic axial tensors in Cartesian coordinates
C
C
      DIMENSION AAT(3,3,NAt)
      character*80 CHAR
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.aat',
     $    FORM='FORMATTED',STATUS='UNKNOWN')
c
       WRITE(40,'(a)') '$aa tensors'
       DO 20 I=1,NAt
       WRITE(40,910) (AAt(k,1,I),k=1,3)
       WRITE(40,910) (AAt(k,2,I),k=1,3)
       WRITE(40,910) (AAt(k,3,I),k=1,3)
 20    CONTINUE

C  close file and exit
C
      WRITE(40,'(a)') '$end'
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
c
  910 Format(10X,3F20.14)
c
      END
c ======================================================================
c
      SUBROUTINE RdAAT(NAt3,AAT,AtomAx)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Attempt to read Cartesian derivatives from <deriv> file
C
C  ARGUMENTS
C
C  NAt     -  number of atoms
C  AAT     -  atomic axial tensors in Cartesian coordinates
C  AtomAx  -  logical flag
C             set to .TRUE. on exit if atomic axial tensors found
C
C
      DIMENSION AAt(NAt3,3)
      CHARACTER CHAR*80
      LOGICAL AtomAx
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
      AtomAx = .FALSE.
C
C  open file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.aat',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  look for the tensors
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      If(CHAR(1:11).NE.'$aa tensors') GO TO 10
C
C  found
C
      DO 20 I=1,NAt3
      READ(40,910) (AAT(I,k),k=1,3)
 20   CONTINUE
      AtomAx = .TRUE.
C
C  close file and exit
C
      IErr = 0
      CLOSE (UNIT=40,STATUS='KEEP')
C
 95   CONTINUE
      RETURN
c
  900 Format(A80)
  910 Format(10X,3F20.14)
c
      END
