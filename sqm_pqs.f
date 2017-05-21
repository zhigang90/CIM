c ==================================================================
c  SQM SCALING MODULE FOR PQS                    JB   Feb 2010
c ==================================================================
c
      SUBROUTINE SQM_PQS(inp,NMem,Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..................................................................
C  ** SQM Scaling Module **
C
C  This module scales a primitive Hessian matrix using a set of
C  standard SQM scale factors. It can read in additional scale
C  factors (either directly embedded in the input or from an
C  external file).  It reads in an existing Cartesian Hessian
C  matrix, converts the Hessian to primitive internals, scales
C  it and converts back to Cartesians. The result should be an
C  SQM Hessian that can be used in the frequency module to give
C  a full SQM vibrational analysis.
C
C  NOTE: This module does NOT have the capability to optimize the
C  scale factors. It simply scales a Hessian using existing (known)
C  scale factors. Scale factor optimization is retained in SQM proper.
C  ....................................................................
C
      DIMENSION Z(NMem)
      Character*256 chopval,jobname
      Character cdum*20,Char*80,Filename*24
      Logical File,dskal
      Common /job/jobname,lenJ
c
      parameter (nopt=4)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
c
      data options/'sqm ','file','prin','node'/
      data ioptyp/0,21,1,0/
c
      PARAMETER (IUnit=1)      !  unit number for checkpoint I/O
c
c -- parse input line
      File = .False.
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c -- filename for scaling data
      if(ifound(2).eq.1) then
        call rmblan(chopval(2),256,LenF)
        If(LenF.GT.24) call nerror(1,'SQM module',
     $     'Filename too long - 24 characters maximum',0,0)
        Filename(1:LenF) = chopval(2)(1:LenF)
        File = .True.
      endif
c
c -- print flag
      if(ifound(3).eq.1) then
        IPrnt = iopval(1,3)
      else
        IPrnt = 1
      endif
c
c -- no default scaling
      if(ifound(4).eq.1) then
        dskal = .False.
      else
        dskal = .True.
      endif
C
C  Read the number of atoms from the <control> file
C  including the number of dummy atoms (if any)
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,MAtoms,rdum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,910) Ndum1,Ndum2
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  only real atoms in this module
C
      NAtoms = MAtoms-Ndum1-Ndum2
C
C  determine the number of lines of scaling data, if any
C
      NLines = 0
      imp = inp
c
      If(File) Then
        OPEN(Unit=IUnit,FILE=Filename(1:LenF),FORM='FORMATTED',
     $       STATUS='OLD')
        imp = IUnit
      EndIF
c
      READ(imp,900) Char
c
      If(Char(1:5).NE.'$scal') Then
        BACKSPACE imp
        If(dskal) GO TO 20
        call nerror(1,'SQM module',
     $    'No Scaling Factors found - Check Your Input!',0,0)
      EndIF
c
 10   CONTINUE
      READ(imp,900) Char
      If(Char(1:5).NE.'$ends') Then
        NLines = NLines+1
        GO TO 10
      EndIf
c
      If(File) Then
        CLOSE (UNIT=IUnit,STATUS='KEEP')
      Else
c -- reposition the input file
        Do I=1,NLines+2
        BACKSPACE imp
        EndDO
      EndIf
c
      If(NLines.EQ.0) Then
        call nerror(1,'SQM module',
     $    'No Scaling Factors found - Check Your Input!',0,0)
      EndIf
c
 20   CONTINUE
      NAt3 = 3*NAtoms
      NInt = 6*NAt3        ! maximum number of primitive internals
      MSkal = NLines+18    ! 18 is the number of lines of default
C                            defined lines of scaling data
C
C  Now allocate the memory pointers
C
      NScr = 2*NInt**2
c
      IXC  = 1
      IHES = IXC  + NAt3
      ISKL = IHES + NAt3**2
      ICTP = ISKL + NInt
      ISTP = ICTP + MSkal
      IVKL = ISTP + MSkal
      ICC  = IVKL + MSkal
      IB   = ICC  + MIN(NAtoms**2,MSkal)
      IGNV = IB   + NAt3*NInt
      IG   = IGNV + NAt3*NInt
      IHPM = IG   + NInt**2
      ktyp = IHPM + NInt**2
      klst = ktyp + NInt
      IScr = klst + 4*NInt
      IEnd = IScr + NScr - 1
C
C  Check memory storage not exceeded
C
      CALL MemChk(NMem,IEnd,7,'SQM_PQS')
C
C
c memory status
c -- assign memory for high water mark (old TEXAS)
      call getmem(IEnd,lastx)
C
C  ------------------------------------------------------------------
C
      CALL SQMM(inp,     NAtoms,  NInt,    File,    Filename,
     $          LenF,    dskal,   NLines,  IPrnt,   Z(IXC),
     $          Z(IHES), Z(ISKL), Z(ICTP), Z(ISTP), Z(IVKL),
     $          Z(ICC),  Z(IB),   Z(IGNV), Z(IG),   Z(IHPM),
     $          Z(ktyp), Z(klst), NScr,    Z(IScr), IErr)
C
C  -------------------------------------------------------------------
C
      call retmem(1)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
      write(6,1100) IEnd,mxmem,memtot
 1100 format(/,' SQM memory status:',/,
     $         ' memory needed=',i15,' high water=',i15,/,
     $         ' total available memory=',i15)
c
c  ----------------------------------------------------------------------
C
C  Exit procedure
C
      RETURN
c
  900 Format(A80)
  910 Format(9X,I4,2X,I4)
c
      END
c ======================================================================
c
      SUBROUTINE SQMM(inp,    NAtoms, NInt,   File,   Filename,
     $                LenF,   dskal,  NLines, IPrnt,  XC,
     $                HESS,   ISKAL,  ctyp,   styp,   VSKAL,
     $                IC,     B,      GINV,   G,      HPRIM,
     $                ktyp,   klist,  NScr,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for SQM treatment
C  SQMM is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(3,NAtoms),HESS(3*NATOMS,3*NATOMS),B(3*NAtoms,NInt),
     $       GINV(3*NAtoms*NInt),G(NInt*NInt),HPRIM(NInt*NInt),
     $       VSKAL(NLines+11),Z(NScr)
      INTEGER IC(NAtoms*NAtoms),ktyp(NInt),klist(4,NInt)
      INTEGER ISKAL(NInt),ctyp(NLines+18),styp(NLines+18)
      Character*4 Group
      Logical File,dskal
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms),clist(4,NLines+18)
c ..................................................
c
      Character cdum*20,Char*80,Filename*24,jobname*256
      Common /job/jobname,lenJ
c
      PARAMETER (IUnit=1)      !  unit number for checkpoint I/O
C
C 
      NAt3 = 3*NAtoms
      NSkal = 0
C
C  Read in the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoord(IUnit,NAtoms,AtSymb,XC,-1,Jnk)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read in the molecular point group
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
c
c -- locate the $symmetry section
  10  CONTINUE
      READ(IUnit,900) Char
      If(Char(1:4).NE.'$sym') GO TO 10
      READ(IUnit,910) Group
      call lowercas(Group,4)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read in Cartesian Hessian matrix
C
      CALL RdHESS(jobname(1:lenJ)//'.hess',lenJ+5,NAt3,IPrnt,
     $            HESS,IHess)
C
C  Read in any additional scale factors
C
      IF(NLines.GT.0) THEN
        If(File) Then
          OPEN(Unit=IUnit,FILE=Filename(1:LenF),FORM='FORMATTED',
     $         STATUS='OLD')
          imp = IUnit
        Else
          imp = inp
        EndIF
c
        CALL RdSkal(imp,    NLines, ctyp,   clist,  styp,
     $              NSkal,  VSKAL)
c
c -- should be unnecessary check
        If(dskal.AND.NSkal.EQ.0) call nerror(2,'SQM module',
     $     'No Default Scaling But NO Scale factors! Check Input.',0,0)
c
        If(File) CLOSE(Unit=IUnit,STATUS='KEEP')
      ENDIF
C
C  Now call the main routine
C
C  ------------------------------------------------------------------
C
      CALL SQM_Scale(NAtoms, Group,  NInt,   NLines, NSkal,
     $               dskal,  IPrnt,  AtSymb, XC,     HESS,
     $               ISKAL,  ctyp,   clist,  styp,   VSKAL,
     $               IC,     B,      GINV,   G,      HPRIM,
     $               ktyp,   klist,  NScr,   Z,      IErr)
C
C  -------------------------------------------------------------------
C
C  Following a successful return from <SQM_Scale> HESS contains
C  the SQM scaled Hessian in Cartesian coordinates
C  Write to file
C
      If(IErr.EQ.0) Then
        CALL WrHESS(jobname(1:lenJ)//'.hess',lenJ+5,5,NAt3,HESS)
      EndIf
c
      RETURN
c
 900  Format(A80)
 910  Format(8X,A4)
c
      END
c ======================================================================
c
      SUBROUTINE SQM_Scale(NAtoms, Group,  NInt,   MSkal,  NSkal,
     $                     dskal,  IPrnt,  AtSymb, XC,     HESS,
     $                     ISKAL,  ctyp,   clist,  styp,   VSKAL,
     $                     IC,     B,      GINV,   G,      HPRIM,
     $                     ktyp,   klist,  NScr,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main routine for SQM module
C  Currently limited functionality as a part of PQS
C  Scales Hessian matrix using default scale factors plus any additional
C  scale factors read in via input file.  Does NOT optimize scale factors
C
C  ARGUMENTS
C
C  NAtoms  -  total number of real atoms (NO dummy atoms)
C  Group   -  molecular point group
C  NInt    -  upper bound to nimber of primitive internals
C  MSkal   -  number of lines of scaling data, if any
C  NSkal   -  number of additional scale factors, if any
C  dskal   -  logical flag for default scaling
C              .True.  - Apply default scaling factors
C              .False. - Apply only user-defined scale factors
C  IPrnt   -  print flag
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  HESS    -  Cartesian Hessian
C  ISKAL   -  vector associating primitives with scale factors
C             empty on entry
C  ctyp    -  array indicating each user-defined primitive type to scale
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  clist    -  list of atoms involved in user-defined scaling coordinate
C               (as atomic symbols, i.e. this is a character array)
C  styp     -  array indicating which scaling parameter to apply
C  VSKAL    -  array of scaling parameter values
C  IC       -  integer storage for connectivity matrix
C  B        -  storage for B matrix
C  GINV     -  storage for transformation matrix (B*B(t))**-1 B
C  U        -  general NInt*NInt scratch storage
C  HPRIM    -  primitive Hessian matrix
C  ktyp     -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C  NScr    -  amount of available scratch storage
C  Z       -  general scratch array
C  IErr    -  error flag      0 - success
C                            -1 - something went wrong
C
C
      REAL*8 XC(3,NAtoms),HESS(3*NATOMS,3*NATOMS),B(3*NAtoms*NInt),
     $       GINV(3*NAtoms*NInt),G(NInt*NInt),HPRIM(NInt*NInt),
     $       VSKAL(NSkal+11),Z(NScr)
      INTEGER IC(NAtoms*NAtoms),ktyp(NInt),klist(4,NInt)
      INTEGER ISKAL(NInt),ctyp(MSkal+18),styp(MSkal+18)
      Character*8 AtSymb(NAtoms),clist(4,MSkal+18)
      Character*4 Group,type(6)
      Logical dskal,First
c
      DATA type /'stre','bend','outp','tors','linc','linp'/
C
C
      IOut = igetival('iout')
      Icon = igetival('icond')
      WRITE(IOut,1200)
      WRITE(ICon,1200)
      If(NSkal.GT.0) WRITE(IOut,1250) NSkal
c
      NAt3 = 3*NAtoms
C
C  generate the primitive internals
C
      CALL GetPRIM(NAtoms, AtSymb, XC,     Group,  NInt,
     $             IC,     Z,      IPrnt,  NPrim,  ktyp,
     $             klist,  IErr)
c
      If(IErr.NE.0) call nerror(3,'SQM module',
     $    'Problems generating primitive internals',0,0)
C
C  assign scaling parameters to primitives
C
      Do I=1,11+NSkal
      IC(I) = 1
      EndDO
c
      If(dskal) CALL DefaultSKAL(MSkal,  NSkal,  ctyp,   clist,  styp,
     $                           VSKAL)
      CALL AssignSKAL(MSkal,  NSkal,  NPrim,  AtSymb, ktyp,
     $                klist,  ctyp,   clist,  styp,   IC,
     $                ISKAL,  NOpt)
C
C  print out scale factors used in this run
C
      WRITE(IOut,1300)
      WRITE(ICon,1300)
      WRITE(IOut,1400)
      WRITE(ICon,1400)
c
      is = 0
      DO 30 I=1,NSkal
      If(IC(I).EQ.-1) GO TO 30      ! scale factor unused
      is = is+1
      first = .True.
      DO 20 J=1,MSkal
      II = ctyp(J)
      If(styp(J).EQ.I) Then
        If(first) Then
          WRITE(IOut,1500) is,type(II),(clist(L,J),L=1,4),VSKAL(I)
          WRITE(ICon,1500) is,type(II),(clist(L,J),L=1,4),VSKAL(I)
          first = .False.
        Else
          WRITE(IOut,1600) type(II),(clist(L,J),L=1,4)
          WRITE(ICon,1600) type(II),(clist(L,J),L=1,4)
        EndIf
      EndIf
  20  CONTINUE
  30  CONTINUE
C
C  form the B-matrix and the generalized inverse
C
      ieig = 1
      iu   = ieig + NPrim
      iv   = iu + NPrim**2
c
      CALL GenINVT(NAtoms, XC,     NPrim,  ktyp,   klist,
     $             IPrnt,  B,      Z(iu),  Z(iv),  Z(ieig),
     $             G,      GINV,   IErr)
c
      If(IErr.NE.0) call nerror(4,'SQM module',
     $    'Problems generating Generalized Inverse matrix',0,0)

C
C  form the primitive Hessian
C
      CALL GetHPRIM(NAt3,   NPrim,  HESS,   GINV,   IPrnt,
     $              Z(iu),  HPRIM)
C
C  scale it
C
      CALL SkalHPRIM(NPrim,  NSkal,  IPrnt,  ISKAL,  VSKAL,
     $               Z(iv),  HPRIM)
C
C  now transform back to Cartesians
C
      CALL GetHCART(NAt3,   NPrim,  IPrnt,  HPRIM,  B,
     $              Z,      HESS)
C
      RETURN
c
 1200 FORMAT(/,17X,'The SQM Scaling Module')
 1250 FORMAT(/,I4,' Additional Scaling Parameters Read in')
 1300 FORMAT(/,'  The Following Scaling Parameters Will Be Used')
 1400 FORMAT('  Parameter',16X,'type',18X,'value')
 1500 FORMAT(4X,I3,4X,A4,3X,4(1X,A6),2X,F7.4,4X,A8)
 1600 FORMAT(11X,A4,3X,4(1X,A6))
c
      END
c ======================================================================
c
      SUBROUTINE AssignSKAL(NEntry, NSkal,  NPrim,  AtSymb, ktyp,
     $                      klist,  ctyp,   clist,  styp,   IFIX,
     $                      ISKAL,  NOpt)
      IMPLICIT INTEGER(A-Z)
C
C
C  Assigns scaling factors to primitives
C
C  ARGUMENTS
C
C  NEntry  -  number of user-defined scalings to assign
C  NSkal   -  number of independent scaling parameters
C  NPrim   -  number of primitives
C  AtSymb  -  atomic symbols
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C  ctyp    -  array indicating each user-defined primitive type
C              (same definition as ktyp)
C  clist   -  list of atoms involved in user-defined scaling coordinate
C              (as atomic symbols, i.e. this is a character array)
C  styp    -  array indicating which scaling factor to apply
C  IFIX    -  array indicating which scaling parameters are fixed
C               0 - scaling parameter to be optimized
C               1 - scaling parameter fixed
C              -1 - scaling factor ignored (primitive not present)
C             ** note: -1 entries set in this routine **
C
C  on exit
C
C  ISKAL   -  array assigning primitives to scaling parameters
C              ISKAL(I) = J  -  primitive I scaled by parameter J
C              ISKAL(I) = 0  -  primitive I unscaled
C  NOpt    -  number of scaling parameters to be optimized
C
C
      DIMENSION ktyp(NPrim),klist(4,NPrim),ctyp(*),styp(*),
     $          IFix(NSkal),ISKAL(NPrim)
      CHARACTER*8 AtSymb(*),clist(4,*)
      CHARACTER*8 a1,a2,a3,a4,c1,c2,c3,c4,X
C
      DATA X/'x       '/       ! X is any atom EXCEPT hydrogen
C
C
C  initialize
C
      CALL IZeroIT(ISKAL,NPrim)
C
C  loop over number of entries in clist
C
      DO 70 I=1,NEntry
c
      Ityp = ctyp(I)
      c1 = clist(1,I)
      c2 = clist(2,I)
      c3 = clist(3,I)
      c4 = clist(4,I)
c
      IF(Ityp.EQ.1) THEN
C
C  stretch
C
       DO 10 J=1,NPrim
       IF(ktyp(J).EQ.1) THEN
        I1 = klist(1,J)
        I2 = klist(2,J)
        a1 = AtSymb(I1)
        a2 = AtSymb(I2)
c
        If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $      (c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) ) Then
             ISKAL(J) = styp(I)
        Else If( (c1.EQ.a2.OR.(c1.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $           (c2.EQ.a1.OR.(c2.EQ.X.AND.a1(1:1).NE.'h')) ) Then
             ISKAL(J) = styp(I)
        EndIf
c
       ENDIF
 10    CONTINUE
C
      ELSE IF(Ityp.EQ.2) THEN
C
C  bend
C
       DO 20 J=1,NPrim
       IF(ktyp(J).EQ.2) THEN
        I1 = klist(1,J)
        I3 = klist(2,J)
        I2 = klist(3,J)
        a1 = AtSymb(I1)
        a2 = AtSymb(I2)
        a3 = AtSymb(I3)
c
        IF(c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) THEN
C
C  central atom matches, check ends
C
         If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $       (c3.EQ.a3.OR.(c3.EQ.X.AND.a3(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         Else If( (c1.EQ.a3.OR.(c1.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $            (c3.EQ.a1.OR.(c3.EQ.X.AND.a1(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         EndIf
        ENDIF
c
       ENDIF
 20    CONTINUE
C
      ELSE IF(Ityp.EQ.3) THEN
C
C  out-of-plane-bend
C
       DO 30 J=1,NPrim
       IF(ktyp(J).EQ.3) THEN
        I1 = klist(1,J)
        I2 = klist(2,J)
        I3 = klist(3,J)
        I4 = klist(4,J)
        a1 = AtSymb(I1)
        a2 = AtSymb(I2)
        a3 = AtSymb(I3)
        a4 = AtSymb(I4)
c
        IF(c4.EQ.a4.OR.(c4.EQ.X.AND.a4(1:1).NE.'h')) THEN
C
C  central atom of out-of-plane bend matches, check others
C
         If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $       (c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $       (c3.EQ.a3.OR.(c3.EQ.X.AND.a3(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         Else If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $            (c2.EQ.a3.OR.(c2.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $            (c3.EQ.a2.OR.(c3.EQ.X.AND.a2(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         Else If( (c1.EQ.a2.OR.(c1.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $            (c2.EQ.a1.OR.(c2.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $            (c3.EQ.a3.OR.(c3.EQ.X.AND.a3(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         Else If( (c1.EQ.a2.OR.(c1.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $            (c2.EQ.a3.OR.(c2.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $            (c3.EQ.a1.OR.(c3.EQ.X.AND.a1(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         Else If( (c1.EQ.a3.OR.(c1.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $            (c2.EQ.a1.OR.(c2.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $            (c3.EQ.a2.OR.(c3.EQ.X.AND.a2(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         Else If( (c1.EQ.a3.OR.(c1.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $            (c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $            (c3.EQ.a1.OR.(c3.EQ.X.AND.a1(1:1).NE.'h')) ) Then
              ISKAL(J) = styp(I)
         EndIf
        ENDIF
c
       ENDIF
 30    CONTINUE
C
      ELSE IF(Ityp.EQ.4) THEN
C
C  torsion
C
       DO 40 J=1,NPrim
       IF(ktyp(J).EQ.4) THEN
        I1 = klist(1,J)
        I2 = klist(2,J)
        I3 = klist(3,J)
        I4 = klist(4,J)
        a1 = AtSymb(I1)
        a2 = AtSymb(I2)
        a3 = AtSymb(I3)
        a4 = AtSymb(I4)
c
        IF( (c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $      (c3.EQ.a3.OR.(c3.EQ.X.AND.a3(1:1).NE.'h')) ) THEN
C
C  central atoms match, check ends
C
         If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $       (c4.EQ.a4.OR.(c4.EQ.X.AND.a4(1:1).NE.'h')) )
     $        ISKAL(J) = styp(I)
c
        ENDIF
cc
        IF( (c2.EQ.a3.OR.(c2.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $      (c3.EQ.a2.OR.(c3.EQ.X.AND.a2(1:1).NE.'h')) ) THEN
C       
C  central atoms match, check ends
C
         If( (c1.EQ.a4.OR.(c1.EQ.X.AND.a4(1:1).NE.'h')) .AND.
     $       (c4.EQ.a1.OR.(c4.EQ.X.AND.a1(1:1).NE.'h')) )
     $        ISKAL(J) = styp(I)
c
        ENDIF
c
       ENDIF
 40    CONTINUE
C
      ELSE IF(Ityp.EQ.5) THEN
C
C  linear coplanar bend
C
       DO 50 J=1,NPrim
       IF(ktyp(J).EQ.5) THEN
        I1 = klist(1,J)
        I2 = klist(2,J)
        I3 = klist(3,J)
        I4 = klist(4,J)
        a1 = AtSymb(I1)
        a2 = AtSymb(I2)
        a3 = AtSymb(I3)
        a4 = AtSymb(I4)
c
        If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $      (c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $      (c3.EQ.a3.OR.(c3.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $      (c4.EQ.a4.OR.(c4.EQ.X.AND.a4(1:1).NE.'h')) )
     $       ISKAL(J) = styp(I)
c
       ENDIF
 50    CONTINUE
C
      ELSE IF(Ityp.EQ.6) THEN
C
C  linear perpendicular bend
C
       DO 60 J=1,NPrim
       IF(ktyp(J).EQ.6) THEN
        I1 = klist(1,J)
        I2 = klist(2,J)
        I3 = klist(3,J)
        I4 = klist(4,J)
        a1 = AtSymb(I1)
        a2 = AtSymb(I2)
        a3 = AtSymb(I3)
        a4 = AtSymb(I4)
c
        If( (c1.EQ.a1.OR.(c1.EQ.X.AND.a1(1:1).NE.'h')) .AND.
     $      (c2.EQ.a2.OR.(c2.EQ.X.AND.a2(1:1).NE.'h')) .AND.
     $      (c3.EQ.a3.OR.(c3.EQ.X.AND.a3(1:1).NE.'h')) .AND.
     $      (c4.EQ.a4.OR.(c4.EQ.X.AND.a4(1:1).NE.'h')) )
     $       ISKAL(J) = styp(I)
c
       ENDIF
 60    CONTINUE
C
      ENDIF
 70   CONTINUE
C
C
C  now check that all scaling parameters have been utilized
C  if not, set corresponding entry in IFix array to -1
C
      DO 80 I=1,NSkal
      DO 75 J=1,NPrim
      If(ISKAL(J).EQ.I) GO TO 80
 75   CONTINUE
C
C  if get here, scaling parameter not found
C
      IFix(I) = -1
 80   CONTINUE
C
C  how many scaling parameters need optimizing?
C
      NOpt = 0
      DO 90 I=1,NSkal
      If(IFix(I).EQ.0) NOpt = NOpt+1
 90   CONTINUE
C
C  warn user if some primitives are unscaled
C
      ITot = 0
      DO 95 I=1,NPrim
      If(ISKAL(I).EQ.0) ITot = ITot+1
 95   CONTINUE
      If(ITot.NE.0) WRITE(6,1000) ITot
cccccccccc
cc      write(6,*) ' IFIX array is:'
cc      do i=1,nskal
cc      write(6,*)i,ifix(i)
cc      enddo
cccccccccc
C
      RETURN
c
 1000 FORMAT('**WARNING** There are ',I4,' Unscaled Primitives in',
     $       ' This Molecule')
c
      END
c ======================================================================
c
      SUBROUTINE DefaultSKAL(MSkal,  NSkal,  ctyp,   clist,  styp,
     $                       VSKAL)
      IMPLICIT INTEGER(A-Z)
C
C  Load up the standard default SQM scale factors
C  These are taken from the original primitive SQM paper
C  J. Baker, A. Jarzecki and P. Pulay, J. Phys. Chem. A 102 (1998) 1412
C
C  ARGUMENTS
C
C  MSkal   -  number of additional lines of scaling data (if any)
C  NSkal   -  number of additional scale factors, if any
C  ctyp    -  array indicating each user-defined primitive type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  clist   -  list of atoms involved in user-defined scaling coordinate
C              (as atomic symbols, i.e. this is a character array)
C  styp    -  array indicating which scaling factor to apply
C  VSKAL   -  array of scale factors
C
C
      DIMENSION ctyp(MSkal+18),styp(MSkal+18)
      REAL*8 VSKAL(NSkal+11)
      CHARACTER*8 clist(4,MSkal+18)
      Character*8 X,blank8
C
C  What we are going to do is pad the existing ctyp, clist and styp arrays
C  containing user-defined scaling data with all the default scaling data.
C  The standard default scaling data goes first.
C
      DATA X/'x       '/       ! X is any atom EXCEPT hydrogen
      DATA blank8/'        '/
C
C  If there is any other stuff, move it to the end
C
      DO I=1,MSkal
      II = I+18
      ctyp(II) = ctyp(I)
      styp(II) = styp(I)+11
      clist(1,II) = clist(1,I)
      clist(2,II) = clist(2,I)
      clist(3,II) = clist(3,I)
      clist(4,II) = clist(4,I)
      EndDO
c
      DO I=1,NSkal
      VSKAL(I+11) = VSKAL(I)
      EndDO
C
C  now load up the defaults
C
c -- X-X Stretch
      ctyp(1) = 1
      styp(1) = 1
      clist(1,1) = X
      clist(2,1) = X
      clist(3,1) = blank8
      clist(4,1) = blank8
      VSKAL(1) = 0.9207d0
c
c -- C-Cl Stretch
      ctyp(2) = 1
      styp(2) = 2
      clist(1,2) = 'c       '
      clist(2,2) = 'cl      '
      clist(3,2) = blank8
      clist(4,2) = blank8
      VSKAL(2) = 1.0438d0
c
c -- C-H Stretch
      ctyp(3) = 1
      styp(3) = 3
      clist(1,3) = 'c       '
      clist(2,3) = 'h       '
      clist(3,3) = blank8
      clist(4,3) = blank8
      VSKAL(3) = 0.9164d0
c
c -- N-H Stretch
      ctyp(4) = 1
      styp(4) = 4
      clist(1,4) = 'n       '
      clist(2,4) = 'h       '
      clist(3,4) = blank8
      clist(4,4) = blank8
      VSKAL(4) = 0.9242d0
c
c -- O-H Stretch
      ctyp(5) = 1
      styp(5) = 5
      clist(1,5) = 'o       '
      clist(2,5) = 'h       '
      clist(3,5) = blank8
      clist(4,5) = blank8
      VSKAL(5) = 0.9527d0
c
c -- X-X-X Bend
      ctyp(6) = 2
      styp(6) = 6
      clist(1,6) = X
      clist(2,6) = X
      clist(3,6) = X
      clist(4,6) = blank8
      VSKAL(6) = 1.0144d0
c
c -- X-X-H Bend
      ctyp(7) = 2
      styp(7) = 7
      clist(1,7) = X
      clist(2,7) = X
      clist(3,7) = 'h       '
      clist(4,7) = blank8
      VSKAL(7) = 0.9431d0
c
c -- H-C-H Bend
      ctyp(8) = 2
      styp(8) = 8
      clist(1,8) = 'h       '
      clist(2,8) = 'c       '
      clist(3,8) = 'h       '
      clist(4,8) = blank8
      VSKAL(8) = 0.9016d0
c
c -- H-N-H Bend
      ctyp(9) = 2
      styp(9) = 9
      clist(1,9) = 'h       '
      clist(2,9) = 'n       '
      clist(3,9) = 'h       '
      clist(4,9) = blank8
      VSKAL(9) = 0.8753d0
c
c -- X-X-X-X Torsion
      ctyp(10) = 4
      styp(10) = 10
      clist(1,10) = X
      clist(2,10) = X
      clist(3,10) = X
      clist(4,10) = X
      VSKAL(10) = 0.9523d0
c
c -- H-X-X-X Torsion
      ctyp(11) = 4
      styp(11) = 10
      clist(1,11) = 'h       '
      clist(2,11) = X
      clist(3,11) = X
      clist(4,11) = X
c
c -- H-X-X-H Torsion
      ctyp(12) = 4
      styp(12) = 10
      clist(1,12) = 'h       '
      clist(2,12) = X
      clist(3,12) = X
      clist(4,12) = 'h       '
c
c -- X-X-X-X  linear coplanar bend
      ctyp(13) = 5
      styp(13) = 11
      clist(1,13) = X
      clist(2,13) = X
      clist(3,13) = X
      clist(4,13) = X
      VSKAL(11) = 0.8847d0
c
c -- H-X-X-X  linear coplanar bend
      ctyp(14) = 5
      styp(14) = 11
      clist(1,14) = 'h       '
      clist(2,14) = X
      clist(3,14) = X
      clist(4,14) = X
c
c -- H-X-X-H  linear coplanar bend
      ctyp(15) = 5
      styp(15) = 11
      clist(1,15) = 'h       '
      clist(2,15) = X
      clist(3,15) = X
      clist(4,15) = 'h       '
c
c -- X-X-X-X  linear perpendicular bend
      ctyp(16) = 6
      styp(16) = 11
      clist(1,16) = X
      clist(2,16) = X
      clist(3,16) = X
      clist(4,16) = X
c
c -- H-X-X-X  linear perpendicular bend
      ctyp(17) = 6
      styp(17) = 11
      clist(1,17) = 'h       '
      clist(2,17) = X
      clist(3,17) = X
      clist(4,17) = X
c
c -- H-X-X-H  linear perpendicular bend
      ctyp(18) = 6
      styp(18) = 11
      clist(1,18) = 'h       '
      clist(2,18) = X
      clist(3,18) = X
      clist(4,18) = 'h       '
c
c -- increment MSkal and NSkal
      MSkal = MSKal+18
      NSkal = NSkal+11
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE GetHCART(NAt3,   NPrim,  IPrnt,  HPRIM,  B,
     $                    VM,     HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine transforms the force constant matrix from
C  primitive internals back to Cartesian coordinates using the B-matrix
C
C     HESS  =  B HPRIM Bt
C
C  NOTE: It is assumed that the Hessian has been computed at a stationary
C        point, so the gradient is zero. If this is not the case, then
C        a second term is needed in the transformation involving the
C        derivative of the B-matrix.
C
C  ARGUMENTS
C
C  NAt3    -  3*number of atoms (dimension of FC)
C  NPrim   -  number of primitive internals
C  IPRNT   -  print flag
C  HPRIM   -  Hessian over primitive internals
C  B       -  B matrix 
C  VM      -  intermediate storage
C  HESS    -  on exit Hessian matrix in Cartesian coordinates
C
C
      REAL*8 HPRIM(NPrim,NPrim),B(NAt3,NPrim),VM(NAt3,NPrim),
     $       HESS(NAt3,NAt3)
C
C
      If(IPRNT.GT.1) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
      EndIf
C
C  now transform the primitive Hessian
C    HESS  =  B HPRIM Bt
C
      DO 40 I=1,NAt3
      DO 40 J=1,NPrim
      Val = 0.0d0
      DO 35 K=1,NPrim
      Val = Val + B(I,K)*HPRIM(K,J)
  35  CONTINUE
      VM(I,J) = Val
  40  CONTINUE
c
      DO 50 I=1,NAt3
      DO 50 J=1,I
      Val = 0.0d0
      DO 45 K=1,NPrim
      Val = Val + VM(I,K)*B(J,K)
  45  CONTINUE
      HESS(I,J) = Val
      HESS(J,I) = Val
  50  CONTINUE
C
C  see what we've got?
C
      If(IPRNT.GT.5) Then
       WRITE(IOut,1100)
       CALL PrntMAT(NAt3,NAt3,NAt3,HESS)
      EndIf
C
      RETURN
c
 1000 FORMAT(/,' Transforming Primitive Hessian Matrix back to',
     $         ' Cartesian Coordinates')
 1100 FORMAT(/,'  Cartesian Force Constant Matrix')
c
      END
c ======================================================================
c
      SUBROUTINE GetPARAM(Char,LMax,NMax,Param,NParam)
      IMPLICIT INTEGER(A-Z)
C
C
C  Decompose a character string into its separate components
C  NOTE:  Each component assumed to be a maximum of 8 characters
C
C  ARGUMENTS
C
C  Char    -  input character string
C  LMax    -  maximum length of each component
C  NMax    -  maximum number of components expected
C
C  on exit
C
C  NParam  -  actual number of separate components
C  Param   -  array to hold each separate component
C
C
      CHARACTER*8 Param(NMax),blank8
      CHARACTER Char*80
C
      DATA blank8/'        '/
C
C
C  initialize
C
      DO 10 I=1,NMax
      Param(I) = blank8
 10   CONTINUE
C
C  split up input string
C  individual components separated by blanks
C
      NParam = 0
      I = 1
 20   CONTINUE
C
C  look for non-blank entry to start component
C
      IF(Char(I:I).NE.' ') THEN
       IStart = I
C
C  now look for blank entry to terminate
C
 30    CONTINUE
       I = I+1
       If(Char(I:I).NE.' ') GO TO 30
C
C  we have a new component
C
       IEnd = I-1
       NParam = NParam+1
       If(NParam.GT.NMax) GO TO 96
       Length = IEnd - IStart + 1
       If(Length.GT.LMax) GO TO 97
       Param(NParam)(1:Length) = Char(IStart:IEnd)
      ENDIF
c
      I = I+1
      If(I.LE.80) GO TO 20
C
C  at this point string length exceeded
C
      RETURN
C
C  ..............................................................
C    ERROR SECTION
C
 96   CONTINUE
      WRITE(0,1000) NMax
      CALL Exit
c
 97   CONTINUE
      WRITE(0,1100) NParam,LMax
      CALL Exit
c
 1000 FORMAT(/,2X,'***ERROR*** Input String Contains More than the',
     $       ' Maximum',/,14X,'Number of ',I2,' Components Expected')
 1100 FORMAT(/,2X,'***ERROR*** Component ',I2,' is Longer than the',
     $       ' Maximum',/,14X,'Component Length ',I2,' Allowed')
c
      END
c ======================================================================
c
      SUBROUTINE GetPRIM(NAtoms, AtSymb, XC,     Group,  NInt,
     $                   ICNNCT, Z,      IPRNT,  NPrim,  ktyp,
     $                   klist,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  -----------------------------------------------------------------
C    P R I M I T I V E   C O O R D I N A T E   G E N E R A T I O N
C  -----------------------------------------------------------------
C
C  This routine is responsible for generating the set of primitive
C  internals directly from input Cartesians based on the atomic
C  connectivity.
C
C  Primitive internals can be:
C      1.   Stretches
C      2.   Planar Bends
C      3.   Out-of-Plane Bends
C      4.   Proper Torsions
C      5.   Linear Coplanar Bend
C      6.   Linear Perpendicular Bend
C
C  For normal molecules typically only stretches, bends and torsions
C  are used. Coplanar and perpendicular bends are used for linear
C  or near-linear arrangements of three or more atoms.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates (in cms and symmetry alligned)
C  Group   -  molecular point group
C  NInt    -  maximum allowed number of primitive internals
C  ICNNCT  -  atomic connectivity matrix
C  Z       -  scratch storage
C             ** WARNING **  At least 8*NAtoms integer storage needed
C  IPRNT   -  flag controlling print out
C
C  on exit
C
C  NPrim   -  actual number of primitives generated
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C
C
      DIMENSION XC(3,NAtoms),ICNNCT(NAtoms,NAtoms)
      DIMENSION ktyp(NInt),klist(4,NInt),Z(*)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 Group
      Dimension IMOL(2)
C
      PARAMETER (thrbnd=0.85)       ! distance ratio for bonding
      PARAMETER (PThrsh=3.05)       ! threshold for near-linear bond angle
C                                     (approximately 175 degrees)
C
C
      If(IPRNT.GT.1) WRITE(6,1000)
c
      IErr = -1
C
C  can't handle linear molecules
C  so trap this and exit
C
      IF(Group.EQ.'c*v '.OR.Group.EQ.'d*h ') THEN
       WRITE(6,1100)
       RETURN
      ENDIF
C
C  generate connectivity data from interatomic distances
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,Z)
C
C  determine connectivity matrix
C
      CALL IZeroIT(ICNNCT,NAtoms*NAtoms)
      CALL CONNECTM(NAtoms,Z,XC,thrbnd,ICNNCT,IErr)
c
      IF(IErr.NE.0) THEN
       If(IPRNT.GT.2) WRITE(6,1200)
       RETURN
      ENDIF
C
C  -----------------------------------
C  Generation of primitive internals
C  -----------------------------------
C
C  attempt automatic assignment of primitive internal
C  coordinates by topological analysis
C
      ITors = 0
      IMOL(1) = 0
      IMOL(2) = NAtoms
      NQ = NAtoms      ! number of unique atoms (not known)
c                        only needed if <TOPOLOGY> includes oop bends
c
      CALL TOPOLOGY(NAtoms, 1,      IMOL,   ITors,  NInt,
     $              Group,  NQ,     ICNNCT, IPrnt,  NPrim,
     $              ktyp,   klist)
C
C  designate scratch pointers
C
      I1 = 1
      I2 = I1 + 4*NAtoms
c
      CALL ChkANG(NAtoms, 1,      IMOL,   XC,     PThrsh,
     $            0,      ICNNCT, IPrnt,  Z(I1),  Z(I2),
     $            NPrim,  ktyp,   klist,  IErr)
C
      RETURN
c
 1000 FORMAT(/,' Attempting to Generate Primitive Internal',
     $         ' Coordinates')
 1100 FORMAT('  SORRY! Cannot Handle Linear Molecules')
 1200 FORMAT(/,2X,'***ERROR*** Unable to Fully Determine Atomic',
     $            ' Connectivity',/,5X,'Missing Homonuclear Bond',
     $            ' Distance in Bond Length Table')
c
      END
c ======================================================================
c
      SUBROUTINE RdSkal(inp,    MSkal,  ctyp,   clist,  styp,
     $                  NSkal,  VSKAL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Read any scaling data from unit inp
C  (This is either from the input or an external file)
C
C  ARGUMENTS
C
C  inp     -  unit number to read
C  MSkal   -  number of lines of scaling information expected
C
C  on exit
C
C  ctyp    -  array indicating each user-defined primitive type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  clist   -  list of atoms involved in user-defined scaling coordinate
C              (as atomic symbols, i.e. this is a character array)
C  styp    -  array indicating which scaling parameter to apply
C  NSkal   -  number of independent scaling parameters
C  VSkal   -  array of scaling parameter values
C
C
      REAL*8 VSkal(MSkal)
      INTEGER ctyp(MSkal),styp(MSkal)
      CHARACTER*8 clist(4,MSkal)
c
      CHARACTER*8 Param(8),blank8
      CHARACTER Char*80,type*4
C
      PARAMETER (ONE=1.0d0)
      DATA blank8/'        '/
C
C
C  Expected Input Format
C  ---------------------
C
C  $scale
C  stre  C  C  1  <value>
C  bend  H  C  H   2  <value>
C  tors  X  X  X  X  3   <value>
C  $endscale
C
C  There should be MSkal lines of input between the delimiters
C
C  initialize clist array
C
      DO 5 I=1,MSkal
      clist(1,I) = blank8
      clist(2,I) = blank8
      clist(3,I) = blank8
      clist(4,I) = blank8
 5    CONTINUE
C
C  now read in scaling parameters
C
 50   CONTINUE
      READ(inp,900,END=97) Char
      If(Char(1:5).NE.'$scal') GO TO 50
c
      DO 60 I=1,MSkal
      READ(inp,900) Char
C
C  decompose this line
C
      CALL GetPARAM(Char,8,8,Param,NParam)
C
C  convert to lower case
C
      Do k=1,NParam
      call lowercas(Param(k),8)
      EndDo
c
      type = Param(1)(1:4)
      clist(1,I) = Param(2)
      clist(2,I) = Param(3)
c
      If(type.EQ.'stre') Then
       ctyp(I) = 1
       II = 4
      Else If(type.EQ.'bend') Then
       ctyp(I) = 2
       clist(3,I) = Param(4)
       II = 5
      Else If(type.EQ.'outp') Then
       ctyp(I) = 3
       clist(3,I) = Param(4)
       clist(4,I) = Param(5)
       II = 6
      Else If(type.EQ.'tors') Then
       ctyp(I) = 4
       clist(3,I) = Param(4)
       clist(4,I) = Param(5)
       II = 6
      Else If(type.EQ.'linc') Then
       ctyp(I) = 5
       clist(3,I) = Param(4)
       clist(4,I) = Param(5)
       II = 6
      Else If(type.EQ.'linp') Then
       ctyp(I) = 6
       clist(3,I) = Param(4)
       clist(4,I) = Param(5)
       II = 6
      Else
       WRITE(6,1000) type
       Call Exit
      EndIf
c
      If(NParam.GE.II) Then
       READ(Param(II),*) JJ
      Else
       JJ = I
      EndIf
      styp(I) = JJ
      If(NParam.GE.II+1) Then
       READ(Param(II+1),*) VSkal(JJ)
      Else
       VSkal(JJ) = ONE
      EndIf
c
 60   CONTINUE
C
C  check for end of scaling input
C
      READ(inp,900) Char
      If(Char(1:1).NE.'$') GO TO 98
C
C  determine actual number of scaling parameters
C  (this is maximum entry in styp array)
C
      NSkal = styp(1)
      DO 65 I=1,MSkal
      If(styp(I).GT.NSkal) NSkal=styp(I)
 65   CONTINUE
c
      RETURN
C
C  ................................................
C    ERROR SECTION
C
 97   CONTINUE
      WRITE(6,1400)
      Call Exit
c
 98   CONTINUE
      WRITE(6,1500)
      Call Exit
c
  900 Format(A80)
 1000 FORMAT(/,2X,'***ERROR*** Unknown Coordinate type found in',
     $              ' input: ',A6)
 1400 FORMAT(/,2X,'***ERROR*** No Scaling Data Found in input')
 1500 FORMAT(/,2X,'***ERROR*** No Terminator Found for Scaling Data')
c
      END
c ======================================================================
c
      SUBROUTINE SkalHPRIM(NPrim,  NSkal,  IPrnt,  ISKAL,  VSkal,
     $                     BSkal,  HPRIM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Scale the primitive Hessian
C
C  ARGUMENTS
C
C  NPrim   -  number of primitives
C  NSkal   -  number of scaling parameters
C  IPrnt   -  print flag
C  ISKAL   -  array indicating which scaling parameters are
C             associated with which primitives
C  VSkal   -  scaling parameters
C  BSkal   -  storage for expanded vector of scaling parameters
C             (unit entries if no scaling for that primitive)
C  HPRIM   -  primitive Hessian (scaled on exit)
C
C
      DIMENSION ISKAL(NPrim),VSkal(NSkal),BSkal(NPrim)
      REAL*8 HPRIM(NPrim,NPrim)
C
      PARAMETER (ONE=1.0d0)
C
C
C  get the expanded scaling vector
C
      DO 10 I=1,NPrim
      IS = ISKAL(I)
      IF(IS.EQ.0) THEN
       BSkal(I) = ONE
      ELSE
       BSkal(I) = VSkal(IS)
      ENDIF
 10   CONTINUE
C
C  scale the Hessian
C
      DO 20 J=1,NPrim
      DO 20 I=1,NPrim
      HPRIM(I,J) = HPRIM(I,J)*SQRT(BSkal(I)*BSkal(J))
 20   CONTINUE
c
      If(IPRNT.GT.6) Then
        WRITE(6,1000)
        CALL PrntMAT(NPrim,NPrim,NPrim,HPRIM)
      EndIf
c
 1000 FORMAT(/,' Scaled Primitive Force Constant Matrix')
c
      RETURN
      END
