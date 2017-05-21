c ==================================================================
c  GEOMETRY OPTIMIZATION ROUTINES A-C          JB   October 1999
c ==================================================================
c
      SUBROUTINE AddCON(IType,  ISign,  II,     JJ,     KK,
     $                  LL,     intcor, ktyp,   klist,  CV)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Finds the primitive associated with a given constraint
C  that is symmetry related to the current constraint and
C  includes it in the constraint vector
C
C  ARGUMENTS
C
C  IType   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane bend
C              4 - fixed dihedral angle
C              5 - fixed linear coplanar bend
C              6 - fixed linear perpendicular bend
C  ISign   -  reflects sign change for oop bend and torsion
C  II, JJ  -  atoms involved in constraint
C  KK, LL
C  intcor  -  number of primitives internals
C  ktyp    -  integer array containing internal coordinate type
C             ** NOTE:  These types are:
C                       1 - stretch
C                       2 - bend
C                       3 - out-of-plane bend
C                       4 - torsion
C                       5 - linear coplanar bend
C                       6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive
C  CV      -  constraint vector
C
C
      DIMENSION ktyp(intcor),klist(4,intcor),CV(intcor)
C
      PARAMETER (One=1.0d0)
C
C
C  Loop over all primitives and look for constraint
C
      DO 10 I=1,intcor
      JType = ktyp(I)
      If(JType.NE.IType) GO TO 10
c
      IF(IType.EQ.1) THEN
cc
       IT = klist(1,I)
       JT = klist(2,I)
       If( (IT.EQ.II.AND.JT.EQ.JJ) .OR.
     $     (IT.EQ.JJ.AND.JT.EQ.II) ) Then
        CV(I) = ISign*One
        RETURN
       EndIf
cc
      ELSE IF(IType.EQ.2) THEN
cc
       IT = klist(1,I)
       JT = klist(2,I)
       KT = klist(3,I)
       If( (IT.EQ.II.AND.JT.EQ.JJ.AND.KT.EQ.KK) .OR.
     $     (KT.EQ.II.AND.JT.EQ.JJ.AND.IT.EQ.KK) ) Then
        CV(I) = ISign*One
        RETURN
       EndIf
cc
      ELSE IF(IType.EQ.3) THEN
cc
       IT = klist(1,I)
       JT = klist(2,I)
       KT = klist(3,I)
       LT = klist(4,I)
       IF(IT.EQ.II.AND.JT.EQ.JJ.AND.KT.EQ.KK.AND.LT.EQ.LL) THEN
        CV(I) = ISign*One
        RETURN
       ELSE IF(IT.EQ.II.AND.KT.EQ.JJ.AND.JT.EQ.KK.AND.LT.EQ.LL) THEN
        CV(I) = -ISign*One
        RETURN
       ENDIF
cc
      ELSE IF(IType.EQ.4) THEN
cc
       IT = klist(1,I)
       JT = klist(2,I)
       KT = klist(3,I)
       LT = klist(4,I)
       If( (IT.EQ.II.AND.JT.EQ.JJ.AND.KT.EQ.KK.AND.LT.EQ.LL).OR.
     $     (LT.EQ.II.AND.KT.EQ.JJ.AND.JT.EQ.KK.AND.IT.EQ.LL) ) Then
        CV(I) = ISign*One
        RETURN
       EndIf
cc
      ELSE IF(IType.EQ.5.OR.IType.EQ.6) THEN
cc
       IT = klist(1,I)
       JT = klist(2,I)
       KT = klist(3,I)
       LT = klist(4,I)
       If(IT.EQ.II.AND.JT.EQ.JJ.AND.KT.EQ.KK.AND.LT.EQ.LL) Then
        CV(I) = ISign*One
        RETURN
       EndIf
cc
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE AddSCAN(MaxCON, NCons,  NAtoms, XC,     Coord,
     $                   I,      J,      K,      L,      NComp,
     $                   NPComp, ICOMP,  IPCOMP, PCWght, npcs,
     $                   IPCS,   PCWS,   ICTYP,  RCON,   ICON)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Adds the primitive being scanned from a potential scan to the
C  list of primitives to be constrained (if any)
C
C  ARGUMENTS
C
C  MaxCON  -  maximum allowed number of constraints
C  NCons   -  current number of constraints
C             (incremented by one on exit)
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  Coord   -  character string indicating scanned primitive type
C  I,J,K,L -  atoms defining scanned primitive
C  NComp   -  number of composite constraints
C             (may be incremented on exit)
C  NPComp  -  number of primitives in composite constraints
C             (may be incremented on exit)
C  ICOMP   -  number of primitives in each composite constraint
C  IPComp  -  constraint type and atoms involved in constraint
C  PCWght  -  weight of each primitive in composite constraint
C  npsc    -  number of primitives in SCAN coordinate
C  IPCS    -  primitive type and atoms involved in coordinate
C  PCWS    -  weight of each primitive in composite coordinate
C  ICTYP   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane-bend
C              4 - fixed dihedral angle
C              5 - fixed linear coplanar bend
C              6 - fixed linear perpendicular bend
C              7 - fixed inverse distance
C              9 - fixed composite coordinate
C  RCON    -  constraint values
C  IC      -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C
      DIMENSION XC(3,NAtoms),ICTYP(NCons),RCON(NCons),ICON(4,NCons)
      DIMENSION ICOMP(*),IPCOMP(5,*),PCWght(*),IPCS(5,npcs),PCWS(npcs)
      Character*4 Coord
C
C
C  increment NCons
C
      NCons = NCons+1
c
      If(NCons.GT.MaxCON) Call nerror(28,'OPTIMIZE module',
     $   'Trying to Impose Too Many Constraints.  Maximum Allowed is',
     $    MaxCON,0)
c
      ICON(1,NCons) = I
      ICON(2,NCons) = J
      ICON(3,NCons) = K
      ICON(4,NCons) = L
C
C  what type of primitive?
C
      If(Coord.EQ.'stre') Then
        ICTYP(NCons) = 1
        CALL StreGRAD(NAtoms,I,J,XC,Th,.false.,jnk)
      Else If(Coord.EQ.'bend') Then
        ICTYP(NCons) = 2
        CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
      Else If(Coord.EQ.'outp') Then
        ICTYP(NCons) = 3
        CALL OutpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
      Else If(Coord.EQ.'tors') Then
        ICTYP(NCons) = 4
        CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
      Else If(Coord.EQ.'linc') Then
        ICTYP(NCons) = 5
        CALL LincGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
      Else If(Coord.EQ.'linp') Then
        ICTYP(NCons) = 6
        CALL LinpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
      Else If(Coord.EQ.'comp') Then
        NComp = NComp+1
        ICON(1,NCons) = 0
        ICON(2,NCons) = 0
        ICON(3,NCons) = 0
        ICON(4,NCons) = 0
        ICOMP(NComp) = npcs
        ICTYP(NCons) = 9
        Th = 0.0d0           ! total value of scanned coordinate
        DO IP=1,npcs
        NPComp = NPComp+1
        ityp = IPCS(1,IP)
        Wght = PCWS(IP)
        I = IPCS(2,IP)
        J = IPCS(3,IP)
        K = IPCS(4,IP)
        L = IPCS(5,IP)
c
        IPCOMP(1,NPComp) = ityp
        IPCOMP(2,NPComp) = I
        IPCOMP(3,NPComp) = J
        IPCOMP(4,NPComp) = K
        IPCOMP(5,NPComp) = L
        PCWght(NPComp) = Wght
c
        If(ityp.EQ.1) Then
c -- distance 
          CALL StreGRAD(NAtoms,I,J,XC,va,.false.,jnk)
        Else If(ityp.EQ.2) Then
c -- bond angle
          CALL AngGRAD(NAtoms,I,J,K,XC,va,.false.,jnk)
        Else If(ityp.EQ.3) Then
c -- out-of-plane-bend
          CALL OutpGRAD(NAtoms,I,J,K,L,XC,va,.false.,jnk)
        Else If(ityp.EQ.4) Then
c -- dihedral angle
          CALL DihGRAD(NAtoms,I,J,K,L,XC,va,.false.,jnk)
        EndIf
        Th = Th+va*Wght
        EndDO
c
      EndIf
C
C  assign value of constraint (scanned coordinate)
C
      RCON(NCons) = Th
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE AddSCANZ(MaxCON, NCons,  NZ,     XC,     Coord,
     $                    VARNAM, IG,     IGEO,   ICTYP,  RCON,
     $                    ICON)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Adds the Z-matrix variable being scanned from a potential scan to the
C  list of primitives to be constrained (if any)
C
C  ARGUMENTS
C
C  MaxCON  -  maximum allowed number of constraints
C  NCons   -  current number of constraints
C             (incremented by one on exit)
C  NZ      -  number of Z-matrix rows (same as number of atoms)
C  XC      -  Cartesian coordinates
C  Coord   -  character string indicating scanned Z-matrix variable
C  VARNAM  -  names of all "variables" in Z matrix
C  IG      -  array determining what to do with Z-matrix parameter
C             0 - optimize it
C             J - assign same value as previous (Jth) variable
C            -J - assign same value  opposite sign
C          1000 - fixed
C  IGEO    -  Z-matrix connectivity
C              Atoms are usually defined with respect to previously
C              defined atoms by a stretch, a bend and a torsion;
C              an alternative is a stretch and two bends  or a stretch
C              a bend and an out-of-plane bend.  The fourth column
C              of IGEO is used to distinguish these cases
C                0 - bend + torsion
C                1 - 2 bends
C                2 - bend + out-of-plane bend
C      ** WARNING - CURRENTLY SET UP FOR TORSIONS ONLY **
C
C  on exit
C
C  ICTYP   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane-bend
C              4 - fixed dihedral angle
C              5 - fixed linear coplanar bend
C              6 - fixed linear perpendicular bend
C  RCON    -  constraint values
C  ICON    -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C
      DIMENSION XC(3,NZ),ICTYP(NCons),RCON(NCons),ICON(4,NCons)
      DIMENSION IGEO(NZ,4),IG(3*NZ)
      CHARACTER*8 VARNAM(3*NZ)
      Character*4 Coord
C
C
C  increment NCons
C
      NCons = NCons+1
c
      If(NCons.GT.MaxCON) Call nerror(28,'OPTIMIZE module',
     $   'Trying to Impose Too Many Constraints.  Maximum Allowed is',
     $    MaxCON,0)
C
C  locate the scanned Z-matrix variable
C
      DO 10 IV=1,3*NZ-6
      If(VARNAM(IV)(1:4).EQ.Coord) Exit
 10   CONTINUE
C
C  which atoms are involved in the scanned variable?
C
      If(IV.LE.NZ-1) Then
c
c -- stretch
        I = IV+1
        J = IGEO(I,1)
        K = 0
        L = 0
        ICTYP(NCons) = 1
        CALL StreGRAD(NZ,I,J,XC,Th,.false.,jnk)
      Else If(IV.GE.NZ.AND.IV.LE.2*NZ-3) Then
c
c -- bend
        IV = IV-(NZ-1)
        I = IV+2
        J = IGEO(I,1)
        K = IGEO(I,2)
        L = 0
        ICTYP(NCons) = 2
        CALL AngGRAD(NZ,I,J,K,XC,Th,.false.,jnk)
      Else If(IV.GT.2*NZ-3) Then
c
c -- torsion or out-of-plane-bend (cannot handle 2 bends)
        IV = IV-(2*NZ-3)
        I = IV+3
        J = IGEO(I,1)
        K = IGEO(I,2)
        L = IGEO(I,3)
c
        If(IGEO(I,4).EQ.0) Then
          ICTYP(NCons) = 4
          CALL DihGRAD(NZ,I,J,K,L,XC,Th,.false.,jnk)
        Else
          ICTYP(NCons) = 3
          CALL OutpGRAD(NZ,I,J,K,L,XC,Th,.false.,jnk)
        EndIF
      EndIf
c
      ICON(1,NCons) = I
      ICON(2,NCons) = J
      ICON(3,NCons) = K
      ICON(4,NCons) = L
C
C  assign value of constraint (scanned primitive)
C
      RCON(NCons) = Th
cc      write(6,*) ' scanned z-matrix variable is: ',I,J,K,L
cc      write(6,*) ' value: ',Th
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE BackToZ(NZ,     NIC,    NVar,    IG,     XINT,
     $                   ZINT,   GEO)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Restore internal coordinates (changed by the optimization
C  step) to the Z-Matrix GEO array.
C  (bohr/radians are converted back to angstrom/degrees)
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  NIC     -  total number of internal coordinates
C  NVar    -  number of variables in Z-Matrix
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  XINT    -  current internal parameter values
C  ZINT    -  full set of previous internal coordinate values
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C
C
      DIMENSION IG(NIC),XINT(NVar),ZINT(NIC),GEO(NZ,3)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
C
C  expand variables to full set of internal coordinates
C
      ToRAD = PI/180.0d0
      NV = 0
c
      DO 10 I=1,NIC
      IF(IG(I).EQ.0) THEN
       NV = NV+1
       ZINT(I) = XINT(NV)
      ELSE IF(IG(I).GT.0.AND.IG(I).NE.1000) THEN
       ZINT(I) = ZINT(IG(I))
      ELSE IF(IG(I).LT.0) THEN
       ZINT(I) = -ZINT(Abs(IG(I)))
      ENDIF
 10   CONTINUE
C
C  Check on number of variables
C
      IF(NV.NE.NVar) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
C  update Z-Matrix GEO array
C
      DO 20 I=1,NZ-1
      GEO(I+1,1) = ZINT(I)/ANTOAU
 20   CONTINUE
c
      IT = NZ-1
c
      DO 30 I=1,NZ-2
      IT = IT+1
      GEO(I+2,2) = ZINT(IT)/ToRAD
 30   CONTINUE
c
      DO 40 I=1,NZ-3
      IT = IT+1
      GEO(I+3,3) = ZINT(IT)/ToRAD
 40   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Wrong Number of Variables in <BackToZ>')
c
      END
c =====================================================================
c
      SUBROUTINE BCHECK(Nat3,   NPrim,  IPRNT,  NInt,   thrsh,
     $                  UT,     B,      BOrth,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Schmidt orthogonalizes the (symmetrized) B-Matrix columns and
C  eliminates linear-dependent coordinates (those coordinates
C  whose B-matrix columns have norms less than thrsh before
C  or after Schmidt orthogonalization)
C
C  ARGUMENTS
C
C  Nat3    -  3*number of atoms
C  NPrim   -  number of primitives
C  thrsh   -  threshold for eliminating coordinate
C  IPRNT   -  flag controlling print out
C  NInt    -  on entry number of natural internal coordinates
C             on exit  number of linearly independent coordinates
C  UT      -  vectors of coefficients defining each internal coordinate
C  B       -  B matrix
C  BOrth   -  Schmidt-Orthogonalized B-matrix columns
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      REAL*8 UT(NPrim,NInt),B(Nat3,NInt),BOrth(Nat3,NInt)
C
C
      IErr = 0
      IC = 0               ! counter for linearly-independent internals
C
      DO 30 I=1,NInt
      IC = IC+1
C
C  check norm of B-matrix column
C
      snorm = SProd(Nat3,B(1,IC),B(1,IC))
      If(IPRNT.GT.5) WRITE(6,1000) I,ic,ic,snorm
      IF(snorm.LT.thrsh) THEN
       If(IPRNT.GT.5) WRITE(6,1010) I
       IC = IC-1
       GO TO 30
      ENDIF
C
C  check for linear dependency
C
      CALL CpyVEC(Nat3,B(1,I),BOrth(1,IC))
C
      DO 20 JC=IC-1,1,-1
      snorm = SProd(Nat3,BOrth(1,JC),BOrth(1,JC))
c
      IF(snorm.LT.thrsh) THEN
       If(IPRNT.GT.5) WRITE(6,2000)
       IErr = -1
       RETURN
      ENDIF
c
      cf = SProd(Nat3,BOrth(1,JC),BOrth(1,IC))/snorm
      DO 10 L=1,Nat3
      BOrth(L,IC) = BOrth(L,IC) - cf*BOrth(L,JC)
 10   CONTINUE
c
      snorm = SProd(Nat3,BOrth(1,IC),BOrth(1,IC))
      If(IPRNT.GT.5) WRITE(6,1000) I,IC,JC,snorm
      IF(snorm.LT.thrsh) THEN
       If(IPRNT.GT.5) WRITE(6,1020) I
       IC = IC-1
       GO TO 30
      ENDIF
 20   CONTINUE
C
C  at this point internal coordinate IC is linearly
C  independent of all previous coordinates
C  store it appropriately in UT
C
      IF(IC.LT.I) CALL CpyVEC(NPrim,UT(1,I),UT(1,IC))

 30   CONTINUE
C
C  at this point we have a total of IC acceptable coordinates
C
      If(IC.LT.NInt.AND.IPRNT.GT.1) WRITE(6,1030) NInt-IC
      NInt = IC
      RETURN
c
 1000 FORMAT(' B-Matrix Column ',I4,' (',I3,') has overlap with',
     $       ' Column (',I3,') of ',G15.8)
 1010 FORMAT(' Internal Coordinate ',I4,' was found to have a Zero',
     $       ' B-Matrix Column'/,
     $       ' probably due to Symmetry - Coordinate Rejected')
 1020 FORMAT(' Internal Coordinate ',I4,' was found to be Linearly',
     $       ' Dependent',/,
     $       ' on previous coordinates - Coordinate Rejected')
 1030 FORMAT(' Eliminated ',I4,' Coordinates due to Linear Dependency')
 2000 FORMAT(/,2X,'***ERROR*** Zero Norm encountered in <BCHECK>')
c
      END
c =====================================================================
c
      SUBROUTINE BINVT(NAtoms, intcor, IPRNT,  B,      BmBt,
     $                 BTinv,  BINV,   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  "Invert" the B-matrix
C  i.e. set up  B**(-1) = Bt*(B*Bt)**(-1)
C  (the m-matrix usual in this formalism is taken as a unit matrix)
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  intcor  -  number of symmetry-independent internal coordinates
C  B       -  symmetrized B-matrix
C  IPRNT   -  print flag (debug)
C
C  on exit ....
C
C  BmBt    -  contains (B*Bt)
C  BTinv   -  contains (B*Bt)**(-1)
C  BINV    -  inverse of B-matrix
C  IErr    -  error flag   0 - success
C                         -1 - could not invert B-matrix
C
      REAL*8 B(3*NAtoms,intcor),BINV(3*NAtoms,intcor)
      REAL*8 BmBt(intcor,intcor),BTinv(intcor,intcor)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=1.0d-5)
C
C
      IErr = -1
C
C  set up (B*Bt)
C
      DO 10 I=1,intcor
      DO 10 J=1,I
      BmBt(I,J) = SProd(3*NAtoms,B(1,I),B(1,J))
      BmBt(J,I) = BmBt(I,J)
 10   CONTINUE
C
C  now invert it
C  (BINV used as scratch)
C
      IF(intcor.EQ.1) THEN
       If(Abs(BmBt(1,1)).LT.thrsh) RETURN
       BTinv(1,1) = One/BmBt(1,1)
       IErr = 0
      ELSE
       CALL INVMAT(BmBt,   intcor, BINV(1,1), BINV(1,2),
     $             BTinv,  IErr)
       If(IErr.NE.0) RETURN
      ENDIF
C
C  set up B**(-1) = Bt*(B*Bt)**(-1)
C
      DO 30 I=1,3*NAtoms
      DO 30 J=1,intcor
      SUM = Zero
      DO 20 K=1,intcor
      SUM = SUM + B(I,K)*BTinv(K,J)
 20   CONTINUE
      BINV(I,J) = SUM
 30   CONTINUE
c
      If(IPRNT.GT.6) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL PrntMAT(intcor,3*NAtoms,intcor,BINV)
      EndIf
C
      RETURN
c
 1000 FORMAT(/,'  Inverse B-Matrix',/)
c
      END
c =====================================================================
c
      SUBROUTINE BMAT2(I,      NAtoms, XC,     NIC,    ktyp,
     $                 klist,  SavTOR, makeq,  makeb,  itor,
     $                 BSkal,  IPRNT,  INDX,   BVal,   QQ,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates value and non-zero entries to B-Matrix column
C  for primitive internal coordinate IC
C
C  Coded from the description in
C  S.Califano, "Vibrational States" (Wiley, London, 1976)
C
C  This routine recognizes several different primitive
C  internal coordinates:
C
C  stre    stretch: a bond distance between 2 atoms
C  bend    a bond angle involving 3 atoms (a planar bend)
C  outp    out-of-plane bend involving 4 atoms
C  tors    standard torsion involving 4 atoms
C  linc    special coordinate to describe bending in a
C          near-linear system (involves 4 atoms)
C          colinear bending of a-b-c in the plane bcd
C  linp    similar to linc but describes bending of a-b-c
C          perpendicular to the plane bcd
C  inv6    inverse stretch to the sixth power multiplied by 2000
C          (experimental coordinate for VderW clusters)
C
C  ARGUMENTS
C
C  I       -  current primitive internal coordinate
C  NATOMS  -  number of atoms
C  XC      -  Cartesian coordinates
C  NIC     -  total number of primitive internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C             (1-7 for the seven primitive types given above)
C  klist   -  list of atoms involved in each primitive
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  makeq   -  integer flag for determining current value of
C             internal coordinate
C               0 - NO;  1 - YES
C  makeb   -  integer flag for calculating B-Matrix column
C               0 - NO;  1 - YES
C  itor    -  integer flag for saving/checking values of primitive
C             torsions (this is done due to possible convergence
C             problems in iterative generation of new Cartesians)
C              0 - no action
C              1 - save initial primitive torsions in SavTOR
C              2 - check current torsions against initial
C                  (may need to change sign of angles > PI)
C  BSkal   -  scaling factor for inverse-distance coordinates
C  IPRNT   -  flag controlling print out
C
C  on exit
C
C  INDX    -  index of non-zero entries in B-Matrix column
C  BVal    -  non-zero B Matrix values
C  QQ      -  value of primitive internal coordinate IC
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      DIMENSION XC(3,NAtoms),ktyp(NIC),klist(4,NIC),SavTOR(NIC)
      DIMENSION INDX(12),BVal(12)
c
      dimension T(3),U(3),V(3),W(3),X(3),Y(3),Z(3),SS(3),TT(3),YY(3),
     1          UU(3),VV(3),WW(3),ZZ(3),UV(12)
      equivalence (UV(1),UU(1)),(UV(4),VV(1)),(UV(7),WW(1)),
     1            (UV(10),ZZ(1))
c
      PARAMETER (ZERO=0.0d0,HALF=0.5d0,ONE=1.0d0,TWO=2.0d0)
      PARAMETER (small=1.0d-6)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      dimension limtyp(7)
      data limtyp/ 2,3,4,4,4,4,2  /
C
C
      IErr = -1
      IOut = ioutfil('iout')
c
      If(IPRNT.GT.6) write(IOut,*) ' Calculating B-Matrix Column:',i
C
C  initialize for this coordinate
C
      NTyp = limtyp(ktyp(I))
c
      IT1 = klist(1,I)
      IT2 = klist(2,I)
      IT3 = klist(3,I)
      IT4 = klist(4,I)
c
      IF(ktyp(I).EQ.1) THEN
cc
C ------STRETCH------
C
       CALL VecDIF(UU,QQ,XC(1,IT1),XC(1,IT2))
       VV(1) = -UU(1)
       VV(2) = -UU(2)
       VV(3) = -UU(3)
cc
      ELSE IF(ktyp(I).EQ.7) THEN
cc
C ------INVERSE-DISTANCE------
C
       UU(1)=XC(1,IT1)-XC(1,IT2)
       UU(2)=XC(2,IT1)-XC(2,IT2)
       UU(3)=XC(3,IT1)-XC(3,IT2)
       RM2 = ONE/(UU(1)*UU(1) + UU(2)*UU(2) + UU(3)*UU(3))
cc       QQ = 2000.0d0*(RM2**3)
cc       RM8 = -QQ*RM2*6.0d0
       QQ = BSkal*SQRT(RM2)
       RM8 = -QQ*RM2
       UU(1) = UU(1)*RM8
       UU(2) = UU(2)*RM8
       UU(3) = UU(3)*RM8
       VV(1) = -UU(1)
       VV(2) = -UU(2)
       VV(3) = -UU(3)
cc
      ELSE IF(ktyp(I).EQ.2) THEN
cc
C ------BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT3))
       CALL VecDIF(V,R2,XC(1,IT2),XC(1,IT3))
       CO = SProd(3,U,V)
       SI = S2(CO)
       SIR1 = SI*R1
       SIR2 = SI*R2
c
       IF(Abs(SIR1).LT.small.OR.Abs(SIR2).LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2200)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       DO 10 L=1,3
       UU(L) = (CO*U(L) - V(L))/SIR1
       VV(L) = (CO*V(L) - U(L))/SIR2
       WW(L) = - (UU(L) + VV(L))
 10    CONTINUE
       QQ = ACOS(CO)
cc
      ELSE IF(ktyp(I).EQ.3) THEN
cc
C ------OUT OF PLANE BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT4))
       CALL VecDIF(V,R2,XC(1,IT2),XC(1,IT4))
       CALL VecDIF(W,R3,XC(1,IT3),XC(1,IT4))
       CO = SProd(3,V,W)
       SI = S2(CO)
c
       IF(SI.LT.small) THEN
c -- three atoms defining plane are linear
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2300)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       CALL Normal(V,W,Z)
       CP = SProd(3,U,Z)
       SJ = S2(CP)
       CQ = SProd(3,W,U)
       CR = SProd(3,V,U)
       D = SJ*SI*SI
c
       IF(SJ.LT.small) THEN
c -- out-of-plane bend is 90 degrees
c -- cannot define B-matrix column
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2400)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ELSE
        ST2 = (CO*CQ-CR)/(R2*D)
        ST3 = (CO*CR-CQ)/(R3*D)
       ENDIF
c
       DO 20 L=1,3
       VV(L) = Z(L)*ST2
       WW(L) = Z(L)*ST3
 20    CONTINUE
       CALL Normal(Z,U,X)
       CALL Normal(U,X,Z)
       DO 30 L=1,3
       UU(L) = Z(L)/R1
       ZZ(L) = - (UU(L) + VV(L) + WW(L))
 30    CONTINUE
c
       CX = -ONE
       If(CP.LT.ZERO) CX = ONE
       QQ = -CX*ACOS(SJ)
C
C  As defined here the out-of-plane bend cannot have a magnitude
C  greater than 90 degrees; we want it to go up to 180 degrees.
C  We use the torsion ILJK to decide this; if magnitude of torsion
C  is greater than 90 degrees then so is the magnitude of the
C  out-of-plane bend
C
cc      Call DihGRAD(NAtoms,IT1,IT4,IT2,IT3,XC,Dih,.False.,Jnk)
cc      If(ABS(Dih).GT.0.5d0*PI) Then
cc        If(QQ.GT.ZERO) QQ = PI-QQ
cc        If(QQ.LT.ZERO) QQ = -(PI+QQ)
cc      EndIf
cc
      ELSE IF(ktyp(I).EQ.4) THEN
cc
C ------TORSION------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT2))
       CALL VecDIF(V,R2,XC(1,IT3),XC(1,IT2))
       CALL VecDIF(W,R3,XC(1,IT3),XC(1,IT4))
       CO = SProd(3,U,V)
       CP = SProd(3,V,W)
       SI = S2(CO)
       SJ = S2(CP)
       SIR1 = SI*R1
       SJR3 = SJ*R3
c
       IF(Abs(SIR1).LT.small.OR.Abs(SJR3).LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2500)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       CALL Normal(U,V,Z)
       CALL Normal(W,V,X)
       DO 40 L=1,3
       UU(L) = Z(L)/SIR1
       ZZ(L) = X(L)/SJR3
       VV(L) = (R1*CO/R2 - ONE)*UU(L) - (R3*CP/R2)*ZZ(L)
       WW(L) = - (UU(L) + VV(L) + ZZ(L))
 40    CONTINUE
       CO = SProd(3,Z,X)
       CALL Cross(Z,X,U)
       SI = SQRT(SProd(3,U,U))
       CP = SProd(3,U,V)
c
       S = ARC1(-CO,SI)
       If(CP.LT.ZERO) S = -S
       QQ = -S
c  .............................................................
c    ** saving/checking of primitive torsions **
c
       IF(itor.EQ.1) THEN
        SavTOR(I) = QQ
       ELSE IF(itor.EQ.2) THEN
        If(Abs(QQ).GT.ONE.AND.SIGN(ONE,QQ).NE.
     $                        SIGN(ONE,SavTOR(I))) Then
         If(IPRNT.GT.3) Then
         write(IOut,*) ' Torsional sign change  primitive torsion',i
         write(IOut,*) '   previous value: ',savtor(i)
         write(IOut,*) ' calculated value: ',qq
         EndIf
         d = PI - Abs(QQ)
         QQ = SIGN(PI+d,SavTOR(I))
         If(IPRNT.GT.3) write(IOut,*) ' reassigned value: ',qq
        EndIf
       ENDIF
c  .............................................................
cc
      ELSE IF(ktyp(I).EQ.5) THEN
cc
C ------LINEAR COPLANAR BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT3))
       CALL VecDIF(V,R2,XC(1,IT4),XC(1,IT3))
       CALL VecDIF(X,R3,XC(1,IT2),XC(1,IT3))
       CO = SProd(3,V,U)
       CP = SProd(3,X,V)
c
       QQ = PI - ACOS(CO) - ACOS(CP)
c
       CALL Normal(V,U,W)
       CALL Normal(U,W,Z)
       CALL Normal(W,V,Y)
       CALL Normal(X,V,W)
       CALL Normal(W,X,U)
       CALL Normal(V,W,T)
C
C  internal coordinate +Ve if atom b moves towards atom d
C
       DO 50 L=1,3
       UU(L) = Z(L)/R1
       VV(L) = U(L)/R3
       ZZ(L) = (Y(L)+T(L))/R2
       WW(L) = - (UU(L) + VV(L) + ZZ(L))
 50    CONTINUE
cc
      ELSE IF(ktyp(I).EQ.6) THEN
cc
C ------LINEAR PERPENDICULAR BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT3))
       CALL VecDIF(V,R2,XC(1,IT4),XC(1,IT3))
       CALL VecDIF(Z,R3,XC(1,IT2),XC(1,IT3))
       CALL Cross(V,U,W)
       WNorm = SQRT(SProd(3,W,W))
c
       IF(WNorm.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2600)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       WNorm = ONE/WNorm
       CALL Vscal(3,WNorm,W)
       CALL Cross(Z,V,X)
       XNorm = SQRT(SProd(3,X,X))
c
       IF(XNorm.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2600)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       XNorm = ONE/XNorm
       CALL Vscal(3,XNorm,X)
       CO = SProd(3,U,X)
       CP = SProd(3,Z,W)
c
       QQ = HALF*(PI - ACOS(CO) - ACOS(CP))
c
       SI = TWO*S2(CO)
       SJ = TWO*S2(CP)
c
       IF(SI.LT.small.OR.SJ.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2700)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       CALL Cross(Z,W,Y)
       CALL Cross(W,Y,T)
       CALL Cross(U,T,Y)
       CALL Cross(V,T,YY)
       CALL Cross(U,X,T)
       CALL Cross(X,T,SS)
       CALL Cross(Z,SS,T)
       CALL Cross(V,SS,TT)
       DO 60 L=1,3
       UU(L) = X(L)/(R1*SI) - YY(L)*WNorm/(R1*SJ)
       VV(L) = W(L)/(R3*SJ) + TT(L)*XNorm/(R3*SI)
       ZZ(L) = Y(L)*WNorm/(R2*SJ) - T(L)*XNorm/(R2*SI)
 60    CONTINUE
       CALL Cross(UU,U,T)
       CALL Cross(U,T,UU)
       CALL Cross(VV,Z,T)
       CALL Cross(Z,T,VV)
       CALL Cross(ZZ,V,T)
       CALL Cross(V,T,ZZ)
       DO 70 L=1,3
       WW(L) = - (UU(L) + VV(L) + ZZ(L))
 70    CONTINUE
cc
      ELSE
cc
C  Unknown Internal coordinate type
C
       WRITE(IOut,2000) I
       WRITE(IOut,2800)
       RETURN
cc
      ENDIF
C
C
C  Set up column of B-Matrix for this coordinate
C
      IF(makeb.EQ.1) THEN
       DO 80 K=1,NTyp
       KK = 3*(K-1)
       L = 3*(klist(K,I)-1)
       INDX(KK+1) = L+1
       INDX(KK+2) = L+2
       INDX(KK+3) = L+3
       BVal(KK+1) = UV(KK+1)
       BVal(KK+2) = UV(KK+2)
       BVal(KK+3) = UV(KK+3)
 80    CONTINUE
      ENDIF
C
C
      IErr = 0
C
      RETURN
c
 2000 FORMAT(/,2X,'***ERROR*** B-Matrix Construction  Internal ',
     $            'Coordinate ',I3)
 2100 FORMAT(5X,'Internal Coordinate ',I3,' became Ill-Conditioned')
 2200 FORMAT(5X,'Bending coordinate is virtually Linear')
 2300 FORMAT(5X,'Atoms in base plane for out-of-plane bend are ',
     $          'virtually Linear')
 2400 FORMAT(5X,'Out-of-plane bend is virtually 90 degrees')
 2500 FORMAT(5X,'Three or more atoms for torsion are virtually ',
     $          'Linear')
 2600 FORMAT(5X,'Linear perpendicular bend is degenerate')
 2700 FORMAT(5X,'Linear perpendicular bend ill-defined')
 2800 FORMAT(5X,'Unknown internal coordinate type')
c
      END
c =====================================================================
c
      SUBROUTINE BMATRIX(NAtoms, XC,     NIC,    ktyp,   klist,
     $                   coeff,  SavTOR, makeq,  makeb,  itor,
     $                   BSkal,  IPRNT,  B,      XINT,   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates internal coordinate values and the B-Matrix
C  over the primitive internal coordinates
C
C  Coded from the description in
C  S.Califano, "Vibrational States" (Wiley, London, 1976)
C
C  This routine recognizes several different primitive
C  internal coordinates:
C
C  stre    stretch: a bond distance between 2 atoms
C  bend    a bond angle involving 3 atoms (a planar bend)
C  outp    out-of-plane bend involving 4 atoms
C  tors    standard torsion involving 4 atoms
C  linc    special coordinate to describe bending in a
C          near-linear system (involves 4 atoms)
C          colinear bending of a-b-c in the plane bcd
C  linp    similar to linc but describes bending of a-b-c
C          perpendicular to the plane bcd
C  inv6    inverse stretch to the sixth power multiplied by 2000
C          (experimental coordinate for VderW clusters)
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NIC     -  total number of primitive internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C             (1-7 for the seven primitive types given above)
C  klist   -  list of atoms involved in each primitive
C  coeff   -  weighting of primitive in non-redundant
C             natural internal coordinate space
C             note: if the weighting of any given primitive is
C                   zero, it can be ignored in B-Matrix construction
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  makeq   -  integer flag for determining current value of
C             internal coordinate
C               0 - NO;  1 - YES
C  makeb   -  integer flag for calculating B-Matrix column
C               0 - NO;  1 - YES
C  itor    -  integer flag for saving/checking values of primitive
C             torsions (this is done due to possible convergence
C             problems in iterative generation of new Cartesians)
C              0 - no action
C              1 - save initial primitive torsions in SavTOR
C              2 - check current torsions against initial
C                  (may need to change sign of angles > PI)
C  BSkal   -  scaling factor for inverse-distance coordinates
C  IPRNT   -  flag controlling print out
C
C  on exit
C
C  B       -  B-Matrix
C  XINT    -  internal coordinate values
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      REAL*8 XC(3,NAtoms),B(3*NAtoms,NIC),XINT(NIC)
      DIMENSION ktyp(NIC),klist(4,NIC),coeff(NIC),SavTOR(NIC)
c
      dimension T(3),U(3),V(3),W(3),X(3),Y(3),Z(3),SS(3),TT(3),YY(3),
     1          UU(3),VV(3),WW(3),ZZ(3),UV(12)
      equivalence (UV(1),UU(1)),(UV(4),VV(1)),(UV(7),WW(1)),
     1            (UV(10),ZZ(1))
c
      PARAMETER (ZERO=0.0d0,HALF=0.5d0,ONE=1.0d0,TWO=2.0d0)
      PARAMETER (small=1.0d-6)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      dimension limtyp(7)
      data limtyp/ 2,3,4,4,4,4,2  /
C
C
      IOut = ioutfil('iout')
      IErr = -1
      If(makeb.EQ.1) CALL ZeroIT(B,3*NAtoms*NIC)
      If(makeq.EQ.1) CALL ZeroIT(XINT,NIC)
C
C  Loop over all primitive internal coordinates
C
      DO 200 I=1,NIC
C
C  ignore this primitive if its weighting in the
C  symmetry non-redundant space is zero
C
      If(coeff(I).EQ.ZERO) GO TO 200
      If(IPRNT.GT.6) write(IOut,*) ' Calculating B-Matrix Column:',i
C
C  initialize for this coordinate
C
      NTyp = limtyp(ktyp(I))
c
      IT1 = klist(1,I)
      IT2 = klist(2,I)
      IT3 = klist(3,I)
      IT4 = klist(4,I)
c
      IF(ktyp(I).EQ.1) THEN
cc
C ------STRETCH------
C
       CALL VecDIF(UU,QQ,XC(1,IT1),XC(1,IT2))
       VV(1) = -UU(1)
       VV(2) = -UU(2)
       VV(3) = -UU(3)
cc
      ELSE IF(ktyp(I).EQ.7) THEN
cc
C ------INVERSE-DISTANCE------
C
       UU(1)=XC(1,IT1)-XC(1,IT2)
       UU(2)=XC(2,IT1)-XC(2,IT2)
       UU(3)=XC(3,IT1)-XC(3,IT2)
       RM2 = ONE/(UU(1)*UU(1) + UU(2)*UU(2) + UU(3)*UU(3))
cc       QQ = 2000.0d0*(RM2**3)
cc       RM8 = -QQ*RM2*6.0d0
       QQ = BSkal*SQRT(RM2)
       RM8 = -QQ*RM2
       UU(1) = UU(1)*RM8
       UU(2) = UU(2)*RM8
       UU(3) = UU(3)*RM8
       VV(1) = -UU(1)
       VV(2) = -UU(2)
       VV(3) = -UU(3)
cc
      ELSE IF(ktyp(I).EQ.2) THEN
cc
C ------BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT3))
       CALL VecDIF(V,R2,XC(1,IT2),XC(1,IT3))
       CO = SProd(3,U,V)
       SI = S2(CO)
       SIR1 = SI*R1
       SIR2 = SI*R2
c
       IF(Abs(SIR1).LT.small.OR.Abs(SIR2).LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2200)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       DO 10 L=1,3
       UU(L) = (CO*U(L) - V(L))/SIR1
       VV(L) = (CO*V(L) - U(L))/SIR2
       WW(L) = - (UU(L) + VV(L))
 10    CONTINUE
       QQ = ACOS(CO)
cc
      ELSE IF(ktyp(I).EQ.3) THEN
cc
C ------OUT OF PLANE BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT4))
       CALL VecDIF(V,R2,XC(1,IT2),XC(1,IT4))
       CALL VecDIF(W,R3,XC(1,IT3),XC(1,IT4))
       CO = SProd(3,V,W)
       SI = S2(CO)
c
c -- test for near-linearity of three atoms defining plane
       IF(SI.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2300)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       CALL Normal(V,W,Z)
       CP = SProd(3,U,Z)
       SJ = S2(CP)
       CQ = SProd(3,W,U)
       CR = SProd(3,V,U)
       D = SJ*SI*SI
c
       IF(SJ.LT.small) THEN
c -- out-of-plane bend is 90 degrees
c -- cannot define B-matrix column
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2400)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ELSE
        ST2 = (CO*CQ-CR)/(R2*D)
        ST3 = (CO*CR-CQ)/(R3*D)
       ENDIF
c
       DO 20 L=1,3
       VV(L) = Z(L)*ST2
       WW(L) = Z(L)*ST3
 20    CONTINUE
       CALL Normal(Z,U,X)
       CALL Normal(U,X,Z)
       DO 30 L=1,3
       UU(L) = Z(L)/R1
       ZZ(L) = - (UU(L) + VV(L) + WW(L))
 30    CONTINUE
c
       CX = -ONE
       If(CP.LT.ZERO) CX = ONE
       QQ = -CX*ACOS(SJ)
C
C  As defined here the out-of-plane bend cannot have a magnitude
C  greater than 90 degrees; we want it to go up to 180 degrees.
C  We use the torsion ILJK to decide this; if magnitude of torsion
C  is greater than 90 degrees then so is the magnitude of the
C  out-of-plane bend
C
cc      Call DihGRAD(NAtoms,IT1,IT4,IT2,IT3,XC,Dih,.False.,Jnk)
cc      If(ABS(Dih).GT.0.5d0*PI) Then
cc        If(QQ.GT.ZERO) QQ = PI-QQ
cc        If(QQ.LT.ZERO) QQ = -(PI+QQ)
cc      EndIf
cc
      ELSE IF(ktyp(I).EQ.4) THEN
cc
C ------TORSION------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT2))
       CALL VecDIF(V,R2,XC(1,IT3),XC(1,IT2))
       CALL VecDIF(W,R3,XC(1,IT3),XC(1,IT4))
       CO = SProd(3,U,V)
       CP = SProd(3,V,W)
       SI = S2(CO)
       SJ = S2(CP)
       SIR1 = SI*R1
       SJR3 = SJ*R3
c
       IF(Abs(SIR1).LT.small.OR.Abs(SJR3).LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2500)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       CALL Normal(U,V,Z)
       CALL Normal(W,V,X)
       DO 40 L=1,3
       UU(L) = Z(L)/SIR1
       ZZ(L) = X(L)/SJR3
       VV(L) = (R1*CO/R2 - ONE)*UU(L) - (R3*CP/R2)*ZZ(L)
       WW(L) = - (UU(L) + VV(L) + ZZ(L))
 40    CONTINUE
       CO = SProd(3,Z,X)
       CALL Cross(Z,X,U)
       SI = SQRT(SProd(3,U,U))
       CP = SProd(3,U,V)
c
       S = ARC1(-CO,SI)
       If(CP.LT.ZERO) S = -S
       QQ = -S
c  .............................................................
c    ** saving/checking of primitive torsions **
c
       IF(itor.EQ.1) THEN
        SavTOR(I) = QQ
       ELSE IF(itor.EQ.2) THEN
        If(Abs(QQ).GT.ONE.AND.SIGN(ONE,QQ).NE.
     $                        SIGN(ONE,SavTOR(I))) Then
         If(IPRNT.GT.3) Then
         write(IOut,*) ' Torsional sign change  primitive torsion',i
         write(IOut,*) '   previous value: ',savtor(i)
         write(IOut,*) ' calculated value: ',qq
         EndIf
         d = PI - Abs(QQ)
         QQ = SIGN(PI+d,SavTOR(I))
         If(IPRNT.GT.3) write(IOut,*) ' reassigned value: ',qq
        EndIf
       ENDIF
c  .............................................................
cc
      ELSE IF(ktyp(I).EQ.5) THEN
cc
C ------LINEAR COPLANAR BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT3))
       CALL VecDIF(V,R2,XC(1,IT4),XC(1,IT3))
       CALL VecDIF(X,R3,XC(1,IT2),XC(1,IT3))
       CO = SProd(3,V,U)
       CP = SProd(3,X,V)
c
       QQ = PI - ACOS(CO) - ACOS(CP)
c
       CALL Normal(V,U,W)
       CALL Normal(U,W,Z)
       CALL Normal(W,V,Y)
       CALL Normal(X,V,W)
       CALL Normal(W,X,U)
       CALL Normal(V,W,T)
C
C  internal coordinate +Ve if atom b moves towards atom d
C
       DO 50 L=1,3
       UU(L) = Z(L)/R1
       VV(L) = U(L)/R3
       ZZ(L) = (Y(L)+T(L))/R2
       WW(L) = - (UU(L) + VV(L) + ZZ(L))
 50    CONTINUE
cc
      ELSE IF(ktyp(I).EQ.6) THEN
cc
C ------LINEAR PERPENDICULAR BEND------
C
       CALL VecDIF(U,R1,XC(1,IT1),XC(1,IT3))
       CALL VecDIF(V,R2,XC(1,IT4),XC(1,IT3))
       CALL VecDIF(Z,R3,XC(1,IT2),XC(1,IT3))
       CALL Cross(V,U,W)
       WNorm = SQRT(SProd(3,W,W))
c
       IF(WNorm.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2600)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       WNorm = ONE/WNorm
       CALL Vscal(3,WNorm,W)
       CALL Cross(Z,V,X)
       XNorm = SQRT(SProd(3,X,X))
c
       IF(XNorm.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2600)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       XNorm = ONE/XNorm
       CALL Vscal(3,XNorm,X)
       CO = SProd(3,U,X)
       CP = SProd(3,Z,W)
c
       QQ = HALF*(PI - ACOS(CO) - ACOS(CP))
c
       SI = TWO*S2(CO)
       SJ = TWO*S2(CP)
c
       IF(SI.LT.small.OR.SJ.LT.small) THEN
        IF(IPRNT.GT.3) THEN
         WRITE(IOut,2000) I
         WRITE(IOut,2700)
         WRITE(IOut,2100) I
        ENDIF
        RETURN
       ENDIF
c
       CALL Cross(Z,W,Y)
       CALL Cross(W,Y,T)
       CALL Cross(U,T,Y)
       CALL Cross(V,T,YY)
       CALL Cross(U,X,T)
       CALL Cross(X,T,SS)
       CALL Cross(Z,SS,T)
       CALL Cross(V,SS,TT)
       DO 60 L=1,3
       UU(L) = X(L)/(R1*SI) - YY(L)*WNorm/(R1*SJ)
       VV(L) = W(L)/(R3*SJ) + TT(L)*XNorm/(R3*SI)
       ZZ(L) = Y(L)*WNorm/(R2*SJ) - T(L)*XNorm/(R2*SI)
 60    CONTINUE
       CALL Cross(UU,U,T)
       CALL Cross(U,T,UU)
       CALL Cross(VV,Z,T)
       CALL Cross(Z,T,VV)
       CALL Cross(ZZ,V,T)
       CALL Cross(V,T,ZZ)
       DO 70 L=1,3
       WW(L) = - (UU(L) + VV(L) + ZZ(L))
 70    CONTINUE
cc
      ELSE
cc
C  Unknown Internal coordinate type
C
       WRITE(IOut,2000) I
       WRITE(IOut,2800)
       RETURN
cc
      ENDIF
C
C
C  Set up column of B-Matrix for this coordinate
C
      If(makeq.EQ.1) XINT(I) = QQ
      IF(makeb.EQ.1) THEN
       DO 80 K=1,NTyp
       KK = 3*(K-1)
       L = 3*(klist(K,I)-1)
       B(L+1,I) = B(L+1,I) + UV(KK+1)
       B(L+2,I) = B(L+2,I) + UV(KK+2)
       B(L+3,I) = B(L+3,I) + UV(KK+3)
 80    CONTINUE
      ENDIF
C
 200  CONTINUE
C
C  end of loop over internal coordinates
C
C
      IF(makeb.EQ.1.AND.IPRNT.GT.6) THEN
       WRITE(IOut,1000)
       CALL PrntMAT(NIC,3*NATOMS,NIC,B)
      ENDIF
C
      IErr = 0
C
      RETURN
c
 1000 FORMAT(/,6X,'B-Matrix',/)
 2000 FORMAT(/,2X,'***ERROR*** B-Matrix Construction  Internal ',
     $            'Coordinate ',I3)
 2100 FORMAT(5X,'Internal Coordinate ',I3,' became Ill-Conditioned')
 2200 FORMAT(5X,'Bending coordinate is virtually Linear')
 2300 FORMAT(5X,'Atoms in base plane for out-of-plane bend are ',
     $          'virtually Linear')
 2400 FORMAT(5X,'Out-of-plane bend is virtually 90 degrees')
 2500 FORMAT(5X,'Three or more atoms for torsion are virtually ',
     $          'Linear')
 2600 FORMAT(5X,'Linear perpendicular bend is degenerate')
 2700 FORMAT(5X,'Linear perpendicular bend ill-defined')
 2800 FORMAT(5X,'Unknown internal coordinate type')
c
      END
c =====================================================================
c
      SUBROUTINE BMatZ(NZ,XZ,IGEO,IPRNT,B)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates B-Matrix for optimization in Z-matrix
C  internal coordinates
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres (including dummy atoms)
C             in Z-matrix
C  XZ      -  Cartesian coordinates of all centres
C  IGEO    -  Z-matrix connectivity
C              IGEO(1,I) = J:  I bonded to J
C              IGEO(2,I) = K:  bond angle I-J-K
C              IGEO(3,I) = L:  torsion I-J-K-L
C  IPRNT   -  print flag
C  B       -  on exit contains B-matrix
C
C
      DIMENSION XZ(3,NZ),IGEO(NZ,3),B(3*NZ,3*NZ)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0)
C
C
C  Zero the B-matrix
C
      CALL ZeroIT(B,9*NZ*NZ)
C
C ------STRETCH------
C
      DO 10 N=1,NZ-1
      I = N+1
      J = IGEO(I,1)
      II = 3*(I-1)
      JJ = 3*(J-1)
      RX = XZ(1,I) - XZ(1,J)
      RY = XZ(2,I) - XZ(2,J)
      RZ = XZ(3,I) - XZ(3,J)
      RIJ = ONE/SQRT(RX*RX+RY*RY+RZ*RZ)
c
      B(II+1,N) =  RX*RIJ
      B(II+2,N) =  RY*RIJ
      B(II+3,N) =  RZ*RIJ
      B(JJ+1,N) = -RX*RIJ
      B(JJ+2,N) = -RY*RIJ
      B(JJ+3,N) = -RZ*RIJ
c
 10   CONTINUE
C
C ------BEND------
C
      NBond = NZ-1
c
      DO 20 M=1,NZ-2
      N = NBond+M
      I = M+2
      J = IGEO(I,1)
      K = IGEO(I,2)
      II = 3*(I-1)
      JJ = 3*(J-1)
      KK = 3*(K-1)
      RIJX = XZ(1,I) - XZ(1,J)
      RIJY = XZ(2,I) - XZ(2,J)
      RIJZ = XZ(3,I) - XZ(3,J)
      RKJX = XZ(1,K) - XZ(1,J)
      RKJY = XZ(2,K) - XZ(2,J)
      RKJZ = XZ(3,K) - XZ(3,J)
      RIJ = ONE/SQRT(RIJX*RIJX+RIJY*RIJY+RIJZ*RIJZ)
      RKJ = ONE/SQRT(RKJX*RKJX+RKJY*RKJY+RKJZ*RKJZ)
      DIJX = RIJX*RIJ
      DIJY = RIJY*RIJ
      DIJZ = RIJZ*RIJ
      DKJX = RKJX*RKJ
      DKJY = RKJY*RKJ
      DKJZ = RKJZ*RKJ
      DOTJ = DIJX*DKJX + DIJY*DKJY + DIJZ*DKJZ
      RSINJ = ONE/SQRT(ONE-DOTJ*DOTJ)
c
      B(II+1,N) = RIJ*RSINJ*(DOTJ*DIJX - DKJX)
      B(II+2,N) = RIJ*RSINJ*(DOTJ*DIJY - DKJY)
      B(II+3,N) = RIJ*RSINJ*(DOTJ*DIJZ - DKJZ)
      B(KK+1,N) = RKJ*RSINJ*(DOTJ*DKJX - DIJX)
      B(KK+2,N) = RKJ*RSINJ*(DOTJ*DKJY - DIJY)
      B(KK+3,N) = RKJ*RSINJ*(DOTJ*DKJZ - DIJZ)
      B(JJ+1,N) = -B(II+1,N) - B(KK+1,N)
      B(JJ+2,N) = -B(II+2,N) - B(KK+2,N)
      B(JJ+3,N) = -B(II+3,N) - B(KK+3,N)
c
 20   CONTINUE
C
C ------TORSION------
C
      NBend = NZ-2
      IOff = NBond + NBend
c
      DO 30 M=1,NZ-3
      N = IOff+M
      I = M+3
      J = IGEO(I,1)
      K = IGEO(I,2)
      L = IGEO(I,3)
      II = 3*(I-1)
      JJ = 3*(J-1)
      KK = 3*(K-1)
      LL = 3*(L-1)
      RJIX = XZ(1,J) - XZ(1,I)
      RJIY = XZ(2,J) - XZ(2,I)
      RJIZ = XZ(3,J) - XZ(3,I)
      RKJX = XZ(1,K) - XZ(1,J)
      RKJY = XZ(2,K) - XZ(2,J)
      RKJZ = XZ(3,K) - XZ(3,J)
      RLKX = XZ(1,L) - XZ(1,K)
      RLKY = XZ(2,L) - XZ(2,K)
      RLKZ = XZ(3,L) - XZ(3,K)
      RJI = ONE/SQRT(RJIX*RJIX+RJIY*RJIY+RJIZ*RJIZ)
      RKJ = ONE/SQRT(RKJX*RKJX+RKJY*RKJY+RKJZ*RKJZ)
      RLK = ONE/SQRT(RLKX*RLKX+RLKY*RLKY+RLKZ*RLKZ)
      DJIX = RJIX*RJI
      DJIY = RJIY*RJI
      DJIZ = RJIZ*RJI
      DKJX = RKJX*RKJ
      DKJY = RKJY*RKJ
      DKJZ = RKJZ*RKJ
      DLKX = RLKX*RLK
      DLKY = RLKY*RLK
      DLKZ = RLKZ*RLK
      DOTJ = -(DJIX*DKJX + DJIY*DKJY + DJIZ*DKJZ)
      DOTK = -(DKJX*DLKX + DKJY*DLKY + DKJZ*DLKZ)
      RSSQJ = ONE/(ONE-DOTJ*DOTJ)
      RSSQK = ONE/(ONE-DOTK*DOTK)
c
      C1X = DJIY*DKJZ - DJIZ*DKJY
      C1Y = DJIZ*DKJX - DJIX*DKJZ
      C1Z = DJIX*DKJY - DJIY*DKJX
      C2X = DKJY*DLKZ - DKJZ*DLKY
      C2Y = DKJZ*DLKX - DKJX*DLKZ
      C2Z = DKJX*DLKY - DKJY*DLKX
c
      COF = -RJI*RSSQJ
      B(II+1,N) = COF*C1X
      B(II+2,N) = COF*C1Y
      B(II+3,N) = COF*C1Z
c
      COF1 = RSSQJ*(RJI - DOTJ*RKJ)
      COF2 = RSSQK*DOTK*RKJ
      B(JJ+1,N) = COF1*C1X - COF2*C2X
      B(JJ+2,N) = COF1*C1Y - COF2*C2Y
      B(JJ+3,N) = COF1*C1Z - COF2*C2Z
c
      COF = RLK*RSSQK
      B(LL+1,N) = COF*C2X
      B(LL+2,N) = COF*C2Y
      B(LL+3,N) = COF*C2Z
c
      B(KK+1,N) = -(B(II+1,N) + B(JJ+1,N) + B(LL+1,N))
      B(KK+2,N) = -(B(II+2,N) + B(JJ+2,N) + B(LL+2,N))
      B(KK+3,N) = -(B(II+3,N) + B(JJ+3,N) + B(LL+3,N))
c
 30   CONTINUE
C
C
C  Pad out the B-matrix with translations and rotations
C
      NIC = 3*NZ-6
C
C ---Translations---
C
      DO 40 I=1,NZ
      II = 3*(I-1)
      B(II+1,NIC+1) = ONE
      B(II+2,NIC+2) = ONE
      B(II+3,NIC+3) = ONE
 40   CONTINUE
C
C ---Rotations---
C
      DO 50 I=1,NZ
      II = 3*(I-1)
C
C  rotate about X
C
      B(II+1,NIC+4) =  ZERO
      B(II+2,NIC+4) =  XZ(3,I)
      B(II+3,NIC+4) = -XZ(2,I)
C
C  rotate about Y
C
      B(II+1,NIC+5) = -XZ(3,I)
      B(II+2,NIC+5) =  ZERO
      B(II+3,NIC+5) =  XZ(1,I)
C
C  rotate about Z
C
      IF(NZ.GT.2) THEN
       B(II+1,NIC+6) =  XZ(2,I)
       B(II+2,NIC+6) = -XZ(1,I)
       B(II+3,NIC+6) =  ZERO
      ENDIF
 50   CONTINUE
C
      If(IPRNT.GT.6) Then
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL PrntMAT(3*NZ,3*NZ,3*NZ,B)
      EndIf
C
      RETURN
c
 1000 FORMAT(/,6X,'B-Matrix (Z-matrix coordinates)')
c
      END
c =====================================================================
c
      SUBROUTINE BSCHMIDT(NAtoms, XC,     NPrim,  ktyp,   klist,
     $                    Coeff,  SavTOR, IPRNT,  NInt,   UT,
     $                    NTrans, TRANS,  NEqATM, thrsh,  XPrim,
     $                    NMem,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for linear dependency of natural internal coordinates
C  by forming the B matrix for the whole set of coordinates and
C  Schmidt-orthogonalizing the B-matrix columns against each
C  other. Internal coordinates with linearly independent B-matrix
C  columns are removed
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NPrim   -  total number of primitive internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  Coeff   -  weighting of primitive in non-redundant
C             natural internal coordinate space
C             note: if the weighting of any given primitive is
C                   zero, it can be ignored in B-Matrix construction
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  IPRNT   -  flag controlling print out
C  NInt    -  number of natural internal coordinates
C             on exit, contains number of linearly independent coordinates
C  UT      -  vectors of coefficients defining each internal coordinate
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  thrsh   -  threshold for eliminating linearly dependent coordinates
C  XPrim   -  contains primitive internal coordinate values on exit
C  NMem    -  amount of working storage available
C             (should be at least 6*NAtoms*NInt)
C  Z       -  work array
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      DIMENSION XC(3,NAtoms),ktyp(NPrim),klist(4,NPrim),Coeff(NPrim),
     $          SavTOR(NPrim),UT(NPrim,NInt),XPrim(NPrim)
      DIMENSION TRANS(9,NTrans),NEqATM(NAtoms,NTrans)
      Dimension Z(NMem)
C
      PARAMETER (One=1.0d0)
C
C
      Nat3 = 3*NAtoms
C
C  get scratch pointers
C
      IB = 1
      IBO = IB + Nat3*NInt
      IEnd = IBO + Nat3*NInt - 1
      CALL MemCHK(NMem,IEnd,8,'BSCHMIDT')
C
C  initialize primitive weights
C
      DO 10 I=1,NPrim
      Coeff(I) = One
 10   CONTINUE

C  form full B matrix
C
      CALL BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $           Coeff,  SavTOR, 1,      1,      1,
     $           BSkal,  IPRNT,  NInt,   UT,     Z(IB),
     $           XPrim,  IErr)
C
C  symmetrize all B matrix columns
C
      If(NTrans.GT.1) CALL BSYM(NAtoms, NInt,   NTrans, TRANS,  NEqATM,
     $                          Z(IBO), Z(IB))
C
C  now check for linear dependency
C
      CALL BCHECK(Nat3,   NPrim,  IPRNT,  NInt,   thrsh,
     $            UT,     Z(IB),  Z(IBO), IErr)
C
C  on exit from <BCHECK> NInt is the number of surviving internal
C  coordinates stored in the first NInt columns of UT
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE BSPARSE(NAtoms, XC,     NPrim,  ktyp,   klist,
     $                   Coeff,  SavTOR, makeq,  makeb,  itor,
     $                   BSkal,  IPRNT,  B,      INDX,   XPRIM,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sparse-matrix construction of B-Matrix in primitive internals
C  storing ONLY the non-zero elements of B (maximum 12 per column)
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NPrim   -  total number of primitive internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  Coeff   -  weighting of primitive in non-redundant
C             natural internal coordinate space
C             note: if the weighting of any given primitive is
C                   zero, it can be ignored in B-Matrix construction
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  makeq   -  integer flag for determining current value of
C             primitive internal coordinate
C               0 - NO;  1 - YES
C  makeb   -  integer flag for calculating B-Matrix column
C               0 - NO;  1 - YES
C  itor    -  integer flag for saving/checking values of primitive
C             torsions (this is done due to possible convergence
C             problems in iterative generation of new Cartesians)
C              0 - no action
C              1 - save initial primitive torsions in SavTOR
C              2 - check current torsions against initial
C                  (may need to change sign of angles > PI)
C  BSkal   -  scale factor for inverse-distance coordinates
C  IPRNT   -  flag controlling print out
C
C  on exit
C
C  B       -  non-zero elements of B-Matrix
C  INDX    -  index to non-zero elements of B-Matrix
C  XPRIM   -  primitive internal coordinate values
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      DIMENSION XC(3,NAtoms),ktyp(NPrim),klist(4,NPrim),Coeff(NPrim),
     $          SavTOR(NPrim),B(12,NPrim),INDX(12,NPrim),XPRIM(NPrim)
C
      PARAMETER (Zero=0.0d0)
C
C
C  initialize
C
      If(makeb.eq.1) Then
       CALL ZeroIT(B,12*NPrim)
       CALL IZeroIT(INDX,12*NPrim)
      EndIf
      If(makeq.eq.1) CALL ZeroIT(XPRIM,NPrim)
C
C  Loop over all NPrim primitive internals
C
      DO 50 IC=1,NPrim
c
      If(Coeff(IC).EQ.Zero) GO TO 50       ! skip if no weight
C
C  construct primitive B-Matrix column for coordinate IC
C
      CALL BMAT2(IC,     NAtoms, XC,     NPrim,  ktyp,
     $           klist,  SavTOR, makeq,  makeb,  itor,
     $           BSkal,  IPRNT, INDX(1,IC),B(1,IC),QQ,IErr)
c
      If(IErr.NE.0) CALL OptExit(9)
      If(makeq.eq.1) XPRIM(IC) = QQ
C
C  on exit from <BMAT2>
C     INDX - index of non-zero entries to B-Matrix column
C     B    - non-zero B Matrix values  (maximum 12)
C
 50   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE BSYM(NAtoms, NInt,   NTrans, TRANS,  NEqATM,
     $                VV,     B)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  symmetrizes each column of the B-Matrix
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NInt    -  number of internal coordinates
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  VV      -  scratch vector
C  B       -  B matrix
C
C
      DIMENSION TRANS(9,NTrans),NEqATM(NAtoms,NTrans)
      REAL*8 VV(3*NAtoms),V(3),B(3*NAtoms,NInt)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      Nat3 = 3*NAtoms
C
C  loop over all B-matrix columns
C
      DO 30 I=1,NInt
      CALL CpyVEC(Nat3,B(1,I),VV)
C
      DO 20 IOP=2,NTrans
C
C  symmetry mapping of B-vector
C
      DO 10 IAT=1,NAtoms
      II = 3*(IAT-1) + 1
      JJ = 3*(NEqATM(IAT,IOP)-1) + 1
      CALL MatVEC(3,TRANS(1,IOP),VV(II),V)
cc      CALL MATAB(3,1,3,TRANS(1,IOP),VV(II),V,0)
      B(JJ,I)   = B(JJ,I)   + V(1)
      B(JJ+1,I) = B(JJ+1,I) + V(2)
      B(JJ+2,I) = B(JJ+2,I) + V(3)
 10   CONTINUE
 20   CONTINUE
 30   CONTINUE
C
C  normalize
C
      BNorm = One/Float(NTrans)
      CALL Vscal(Nat3*NInt,BNorm,B)
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $                 Coeff,  SavTOR, makeq,  makeb,  itor,
     $                 BSkal,  IPRNT,  NDEG,   UT,     B,
     $                 XPRIM,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Direct construction of B-Matrix in delocalized internal coordinates
C  [primitive B matrix is NOT explicitly formed - the non-zero elements
C   of each column (maximum 12) multiply UT directly to form B(UT)]
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NPrim   -  total number of primitive internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  Coeff   -  weighting of primitive in non-redundant
C             natural internal coordinate space
C             note: if the weighting of any given primitive is
C                   zero, it can be ignored in B-Matrix construction
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  makeq   -  integer flag for determining current value of
C             primitive internal coordinate
C               0 - NO;  1 - YES
C  makeb   -  integer flag for calculating B-Matrix column
C               0 - NO;  1 - YES
C  itor    -  integer flag for saving/checking values of primitive
C             torsions (this is done due to possible convergence
C             problems in iterative generation of new Cartesians)
C              0 - no action
C              1 - save initial primitive torsions in SavTOR
C              2 - check current torsions against initial
C                  (may need to change sign of angles > PI)
C  BSkal   -  scale factor for inverse-distance coordinates
C  IPRNT   -  flag controlling print out
C  NDEG    -  number of active delocalized internal coordinates
C  UT      -  vectors of coefficients defining each internal coordinate
C
C  on exit
C
C  B       -  B-Matrix
C  XPRIM   -  primitive internal coordinate values
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      DIMENSION XC(3,NAtoms),ktyp(NPrim),klist(4,NPrim),Coeff(NPrim),
     $          SavTOR(NPrim),UT(NPrim,NDEG),B(3*NAtoms,NDEG),
     $          XPRIM(NPrim)
      Dimension INDX(12),BVal(12),kval(7)
C
      PARAMETER (Zero=0.0d0)
C
      Data kval/ 6, 9, 12, 12, 12, 12, 6 /
C
C
C  initialize
C
      If(makeb.eq.1) CALL ZeroIT(B,3*NAtoms*NDEG)
      If(makeq.eq.1) CALL ZeroIT(XPRIM,NPrim)
C
C  Loop over all NPrim primitive internals
C
      DO 50 IC=1,NPrim
c
      If(Coeff(IC).EQ.Zero) GO TO 50       ! skip if no weight
C
C  construct primitive B-Matrix column for coordinate IC
C
      CALL BMAT2(IC,     NAtoms, XC,     NPrim,  ktyp,
     $           klist,  SavTOR, makeq,  makeb,  itor,
     $           BSkal,  IPRNT,  INDX,   BVal,   QQ,  IErr)
c
      If(IErr.NE.0) CALL OptExit(9)
      If(makeq.eq.1) XPRIM(IC) = QQ
      if(makeb.eq.0) GO TO 50
C
C  on exit from <BMAT2>
C     INDX - index of non-zero entries to B-Matrix column
C     BVal - non-zero B Matrix values  (maximum 12)
C
      nval = kval(ktyp(IC))
      DO 20 I=1,nval
      II = INDX(I)
      BB = BVal(I)
      DO 10 J=1,NDEG
      B(II,J) = B(II,J) + UT(IC,J)*BB
 10   CONTINUE
 20   CONTINUE
c
 50   CONTINUE
c
      If(makeb.EQ.1.AND.IPRNT.GT.6) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL PrntMAT(NDEG,3*NAtoms,NDEG,B)
      EndIf
C
      RETURN
c
 1000 FORMAT(/,6X,'B-Matrix',/)
c
      END
c =====================================================================
c
      SUBROUTINE ChkANG(NAtoms, NMol,   IMOL,   XC,     PThrsh,
     $                  IType,  ICnnct, IPRNT,  klin,   ktor,
     $                  NPrim,  ktyp,   klist,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for and eliminates/replaces near-linear bond angles
C  by colinear bends and modifies torsions.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NMol    -  number of molecules (for, e.g., cluster optimizations)
C  IMOL    -  pointers to start/end of molecules in XC array
C  XC      -  Cartesian coordinates
C  PThrsh  -  threshold for near-linearity
C             (angles > PThrsh are replaced)
C  IType   -  flag for cluster/surface optimizations
C               0 - standard optimization
C               1 - adsorbate-surface optimization
C               2 - VdW cluster  (use 1/R distance coordinates)
C  (If IType=1; then simply eliminate near-linear bond angles)
C  ICnnct  -  connectivity matrix
C  IPRNT   -  flag controlling printout
C  klin    -  scratch array for linear bend coordinates
C  ktor    -  scratch array for additional torsions
C  NPrim   -  total number of primitives so far
C             WARNING: could be modified by this routine
C  ktyp    -  integer array containing internal coordinate type
C             ** NOTE:  These types are:
C                       1 - stretch
C                       2 - bend
C                       3 - out-of-plane bend
C                       4 - torsion
C                       5 - linear coplanar bend
C                       6 - linear perpendicular bend
C                assumed to be in this order
C             ** IMPORTANT: NO COLINEAR BENDS DEFINED ON ENTRY **
C  klist   -  list of atoms involved in each primitive
C  IErr    -  error flag
C               0 - success
C              -1 - linear bend found, but could not be replaced
C                   by colinear bend pair
C
      REAL*8 XC(3,NAtoms)
      DIMENSION IMOL(NMol+1),ICnnct(NAtoms,NAtoms),klin(4,NAtoms)
      DIMENSION ktor(4,NAtoms),ktyp(NPrim),klist(4,NPrim)
      Logical next
C
      Data next/.FALSE./
C
C
      IOut = ioutfil('iout')
      IErr = -1
      NStre = 0
C
C
C  Loop over each molecule separately
C
      DO 100 Mol=1,NMol
      IStrt = IMOL(Mol)+1
      IEnd  = IMOL(Mol+1)
C
C  determine how many primitives of each type there are for
C  this molecule (assume primitives are in order, stretches
C  first, for each molecule)
C
      NBend = 0
      NOutp = 0
      NTor = 0
c
      DO 10 I=NStre+1,NPrim
      If(ktyp(I).EQ.4) Then
       NTor = NTor+1
      Else If(ktyp(I).EQ.2) Then
       NBend = NBend+1
       next = .TRUE.
      Else If(ktyp(I).EQ.1) Then
       If(next) Exit
       NStre = NStre+1
      Else
       NOutp = NOutp+1
      EndIf
 10   CONTINUE
c
      next = .FALSE.
      MPrim = I-1
      NPrim0 = NPrim
C
C  Loop over all bends
C
      mlin = 0
      mtor = 0
      melim = 0
c
      DO 50 IQ=NStre+1,NStre+NBend
      I = klist(1,IQ)
      J = klist(3,IQ)
      K = klist(2,IQ)
      CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
c
      IF(Th.GT.PThrsh) THEN
C
C  mark this bend for elimination
C
       melim = melim+1
       ktyp(IQ) = 0
       If(IPRNT.GT.2) WRITE(IOut,1000) I,J,K
C
C  check all torsions
C  mark any torsions involving IJK for elimination
C
       DO 15 JQ=NStre+NBend+NOutp+1,NPrim
       IF(J.EQ.klist(2,JQ)) THEN
        If( (I.EQ.klist(1,JQ).AND.K.EQ.klist(3,JQ)) .OR.
     $      (I.EQ.klist(3,JQ).AND.K.EQ.klist(1,JQ)) ) ktyp(JQ) = 0
       ELSE IF(J.EQ.klist(3,JQ)) THEN
        If( (I.EQ.klist(2,JQ).AND.K.EQ.klist(4,JQ)) .OR.
     $      (I.EQ.klist(4,JQ).AND.K.EQ.klist(2,JQ)) ) ktyp(JQ) = 0
       ENDIF
 15    CONTINUE
C
C  replace planar bend by appropriate colinear bends
C  (for surface, simply eliminate without replacing)
C
       If(IType.EQ.1.AND.Mol.EQ.1) GO TO 50
C ................................................................
C
C  linear bend
C  look for forth (non-linear) atom connected to atoms I or K
C  to define planar/perpendicular bends
C
       mlin0 = mlin
       DO 20 IAtm=IStrt,IEnd
       If(IAtm.EQ.J) GO TO 20
       IF(ICnnct(I,IAtm).EQ.1) THEN
        CALL AngGRAD(NAtoms,IAtm,I,J,XC,Th,.false.,jnk)
        IF(Th.LE.PThrsh) THEN
         mlin = mlin + 1
         klin(1,mlin) = K
         klin(3,mlin) = J
         klin(2,mlin) = I
         klin(4,mlin) = IAtm
         If(IPRNT.GT.2) WRITE(IOut,1010) K,J,I,IAtm
        ENDIF
       ELSE IF(ICnnct(K,IAtm).EQ.1) THEN
        CALL AngGRAD(NAtoms,J,K,IAtm,XC,Th,.false.,jnk)
        IF(Th.LE.PThrsh) THEN
         mlin = mlin + 1
         klin(1,mlin) = I
         klin(3,mlin) = J
         klin(2,mlin) = K
         klin(4,mlin) = IAtm
         If(IPRNT.GT.2) WRITE(IOut,1010) I,J,K,IAtm
        ENDIF
       ENDIF
 20    CONTINUE
C
C  If we haven't found a colinear bend then we either have 5+
C  atoms near-linear OR we have a terminal linear atom which
C  situation may be salvagable
C
       IF(mlin.EQ.mlin0) THEN
        II = 0
        KK = 0
        DO 25 IAtm=IStrt,IEnd
        II = II + ICnnct(I,IAtm)
        KK = KK + ICnnct(K,IAtm)
 25     CONTINUE
        If(II.EQ.1.AND.KK.EQ.1) GO TO 29     ! both atoms are terminal
C
C  we have one terminal atom
C  relabel so that I is terminal and try to define colinear bend
C
        IF(II.GT.1) THEN
         II = I
         I = K
         K = II
        ENDIF
c
        DO 27 IAtm=IStrt,IEnd
        If(IAtm.EQ.J) GO TO 27
        IF(ICnnct(K,IAtm).EQ.1) THEN
C
C  look for other atoms connected to IAtm
C
         DO 26 JAtm=IStrt,IEnd
         If(JAtm.EQ.K) GO TO 26
         IF(ICnnct(IAtm,JAtm).EQ.1) THEN
          CALL AngGRAD(NAtoms,JAtm,IAtm,K,XC,Th,.false.,jnk)
          IF(Th.LE.PThrsh) THEN
           mlin = mlin + 1
           klin(1,mlin) = -I
           klin(3,mlin) = K
           klin(2,mlin) = IAtm
           klin(4,mlin) = JAtm
           If(IPRNT.GT.2) WRITE(IOut,1000) I,J,K,I,K,IAtm,JAtm
          ENDIF
         ENDIF
 26      CONTINUE
        ENDIF
 27     CONTINUE
       ENDIF
C
 29    CONTINUE
       IF(mlin.EQ.mlin0) THEN
        WRITE(IOut,1100) I,J,K
cc        RETURN
       ENDIF
C
C  mark this bend for elimination
C
       ktyp(IQ) = 0
C
C  check all torsions
C  mark any torsions involving IJK for elimination
C
       DO 30 JQ=NStre+NBend+NOutp+1,NPrim
       IF(J.EQ.klist(2,JQ)) THEN
        If( (I.EQ.klist(1,JQ).AND.K.EQ.klist(3,JQ)) .OR.
     $      (I.EQ.klist(3,JQ).AND.K.EQ.klist(1,JQ)) ) ktyp(JQ) = 0
       ELSE IF(J.EQ.klist(3,JQ)) THEN
        If( (I.EQ.klist(2,JQ).AND.K.EQ.klist(4,JQ)) .OR.
     $      (I.EQ.klist(4,JQ).AND.K.EQ.klist(2,JQ)) ) ktyp(JQ) = 0
       ENDIF
 30    CONTINUE
C
C  If there is anything attached to the terminal linear atom
C  of the colinear bend, need to add the torsion
C
       DO 40 JQ=mlin0+1,mlin
       I = klin(1,JQ)
ccc
       If(I.LT.0) THEN
C
C  I is a terminal atom - no torsions possible
C
        klin(1,JQ) = -I
        GO TO 40
       EndIf
ccc
       J = klin(3,JQ)
       K = klin(2,JQ)
       L = klin(4,JQ)
       DO 39 IAtm=IStrt,IEnd
       IF(ICnnct(I,IAtm).EQ.1.AND.IAtm.NE.J) THEN
        CALL AngGRAD(NAtoms,J,I,IAtm,XC,Th,.false.,jnk)
        IF(Th.LE.PThrsh) THEN
C
C  check existing torsions for duplication
C
         DO 34 itor=1,mtor
         If(ktor(1,itor).EQ.IAtm.AND.ktor(2,itor).EQ.I.AND.
     $      ktor(3,itor).EQ.K.AND.ktor(4,itor).EQ.L) GO TO 39
 34      CONTINUE
c
         mtor = mtor+1
         ktor(1,mtor) = L
         ktor(2,mtor) = K
         ktor(3,mtor) = I
         ktor(4,mtor) = IAtm
         If(IPRNT.GT.4) WRITE(IOut,1200) L,K,I,IAtm
        ELSE
C
C  -- SPECIAL CASE --
C  see if IAtm has other connected atoms that can be used
C  to define the torsion
C
         DO 36 JAtm=IStrt,IEnd
         IF(ICnnct(IAtm,JAtm).EQ.1.AND.JAtm.NE.I) THEN
          CALL AngGRAD(NAtoms,I,IAtm,JAtm,XC,Th,.false.,jnk)
          IF(Th.LE.PThrsh) THEN
C
C  check existing torsions for duplication
C
           DO 35 itor=1,mtor
           If(ktor(1,itor).EQ.JAtm.AND.ktor(2,itor).EQ.IAtm.AND.
     $        ktor(3,itor).EQ.K.AND.ktor(4,itor).EQ.L) GO TO 36
 35        CONTINUE
c
           mtor = mtor+1
           ktor(1,mtor) = L
           ktor(2,mtor) = K
           ktor(3,mtor) = IAtm
           ktor(4,mtor) = JAtm
           If(IPRNT.GT.4) WRITE(IOut,1200) L,K,IAtm,JAtm
          ENDIF
         ENDIF
 36      CONTINUE
        ENDIF
       ENDIF
 39    CONTINUE
 40    CONTINUE
      ENDIF
c
 50   CONTINUE
C
      If(melim.EQ.0) GO TO 95          ! no change in primitives
C
C  update ktyp/klist arrays
C  remove all bends & torsions marked for elimination
C
C  no change to stretches
C  bends
C
      intcor = NStre
      DO 60 I=NStre+1,NStre+NBend
      IF(ktyp(I).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = ktyp(I)
       klist(1,intcor) = klist(1,I)
       klist(2,intcor) = klist(2,I)
       klist(3,intcor) = klist(3,I)
      ENDIF
 60   CONTINUE
C
C  torsions and ALL primitives for unprocessed molecules
C
      DO 70 I=NStre+NBend+1,NPrim
      IF(ktyp(I).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = ktyp(I)
       klist(1,intcor) = klist(1,I)
       klist(2,intcor) = klist(2,I)
       klist(3,intcor) = klist(3,I)
       klist(4,intcor) = klist(4,I)
      ELSE IF(IPRNT.GT.2) THEN
       WRITE(IOut,1300) klist(1,I),klist(2,I),klist(3,I),klist(4,I)
      ENDIF
 70   CONTINUE
C
C  set the number of original primitives for this molecule
C  that have been eliminated
C
      melim = NPrim - intcor
C
C  now append the new primitives for the current molecule to the
C  END of the existing list of ALL primitives (for ALL molecules)
C
C  include the extra torsions
C
      DO 80 I=1,mtor
      intcor = intcor+1
      ktyp(intcor) = 4
      klist(1,intcor) = ktor(1,I)
      klist(2,intcor) = ktor(2,I)
      klist(3,intcor) = ktor(3,I)
      klist(4,intcor) = ktor(4,I)
 80   CONTINUE
C
C  include the new planar/perpendicular bends
C
      DO 90 I=1,mlin
      intcor = intcor+1
      ktyp(intcor) = 5
      klist(1,intcor) = klin(1,I)
      klist(2,intcor) = klin(2,I)
      klist(3,intcor) = klin(3,I)
      klist(4,intcor) = klin(4,I)
      intcor = intcor+1
      ktyp(intcor) = 6
      klist(1,intcor) = klin(1,I)
      klist(2,intcor) = klin(2,I)
      klist(3,intcor) = klin(3,I)
      klist(4,intcor) = klin(4,I)
 90   CONTINUE
c
      NPrim = intcor
      If(IPRNT.GT.2) WRITE(IOut,1400) NPrim
c
 95   CONTINUE
C
C  reset NStre to number of ORIGINAL primitives remaining so far
C  (ignore any NEW primitives that have been appended to the list)
C
      NStre = MPrim - melim
cc      write(6,*) ' End of molecule:',mol,' MPrim:',mprim,
cc     $           ' melim:',melim
cc      write(6,*) ' New value for NStre:',nstre
cc      write(6,*) ' New value for NPrim:',nprim
 100  CONTINUE
c
      IErr = 0
C
      RETURN
c
 1000 FORMAT(' Eliminating near-linear Angle ',3I5)
 1010 FORMAT(' Adding colinear bends ',4I5)
cc 1100 FORMAT(/,2X,'***ERROR*** Angle ',3I5,'   is near-linear',/,
cc     $         14X,'But No atom available to define colinear bend')
 1100 FORMAT('**WARNING** Angle ',3I5,'   is near-linear',/,
     $       '  No atom available to define colinear bend - Eliminated')
 1200 FORMAT(' Adding torsion ',4I5)
 1300 FORMAT(' Eliminating torsion ',4I5)
 1400 FORMAT(' There are now ',I6,' Primitive Internals')
c
      END
c =====================================================================
c
      SUBROUTINE ChkCART(NAT3,IPRNT,DCart,XOld,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks Cartesian displacement during back-transformation
C  If too large, scale back
C
C  ARGUMENTS
C
C  NAT3    -  dimension of Cartesian space (3*NAtoms)
C  IPRNT   -  print flag
C  DCart   -  storage for Cartesian displacement
C  XOld    -  Cartesian coordinates previous cycle
C  XC      -  current Cartesian coordinates
C
C
      DIMENSION DCart(NAT3),XOld(NAT3),XC(NAT3)
      PARAMETER (DMax=0.05d0)
C
C
      IOut = ioutfil('iout')
c
      Diff = 0.0d0
      DO 10 I=1,NAT3
      DCart(I) = XC(I) - XOld(I)
      If(Abs(DCart(I)).GT.Diff) Diff = Abs(DCart(I))
 10   CONTINUE
      If(IPRNT.GT.2) WRITE(IOut,1000) Diff
c
      IF(Diff.GT.DMax) THEN
       Skal = DMax/Diff
       DO 20 I=1,NAT3
       XC(I) = XOld(I) + Skal*DCart(I)
 20    CONTINUE
       If(IPRNT.GT.2) WRITE(IOut,1100) Skal
      ENDIF
c
      RETURN
c
 1000 FORMAT(' Maximum Cartesian Displacement is ',F10.6)
 1100 FORMAT('**WARNING** Scaling Cartesian Displacement by ',F10.6)
c
      END
c =====================================================================
c
      SUBROUTINE ChkCMPLT(NAtoms,IC,XC,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks that the connectivity matrix forms a "closed loop"
C  with all atoms connected
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IC      -  connectivity matrix
C  XC      -  Cartesian coordinates
C  IErr    -  error flag on exit
C              0 - OK; -1 - problem
C
C
      DIMENSION XC(3,NAtoms),IC(NAtoms,NAtoms)
c .........................................................
c -- F90 dynamically allocated memory
      INTEGER I1(NAtoms),I2(NAtoms),I3(NAtoms)
c .........................................................
C
      Parameter (Big=1000.0d0)
C
C
      IErr = 0      ! should not be any problems
c
      num = 0
      CALL IZeroIT(I1,NAtoms)
C
C  First find a connectivity
C
      DO 10 I=2,NAtoms
      DO 10 J=1,I
      If(IC(I,J).EQ.1) GO TO 11
 10   CONTINUE
C
C  should not get here unless NOTHING is bonded
C
      If(NAtoms.EQ.2) Then
        IC(1,2) = 1
        IC(2,1) = 1
      Else
        IErr = -1
      EndIf
      RETURN
C
 11   CONTINUE
      num = num+1
      I1(I) = 1
      I2(num) = I
      num = num+1
      I1(J) = 1
      I2(num) = J
C
C  now search among known connected atoms for further atoms
C  connected to any of the existing atoms
C
      num0 = 0
 20   CONTINUE
      numS = num        ! store number of atoms in group
      DO 30 inum=num0+1,numS
      I = I2(inum)
      DO 25 J=1,NAtoms
      If(IC(I,J).EQ.1.AND.I1(J).EQ.0) Then
c
c -- found a new atom not previously connected to existing group
        num = num+1
        I1(J) = 1
        I2(num) = J
      EndIf
 25   CONTINUE
 30   CONTINUE
C
C  check if num has changed
C
      num0 = numS
      If(num.GT.numS.AND.num.LT.NAtoms) GO TO 20
C
C  At this point we are finished
C  Either everything is connected, in which case we are done OR
C  not all atoms are connected
C
      If(num.EQ.NAtoms) RETURN
C
C  find which atoms are not connected and store them in I3
C
      num3 = 0
      DO 40 I=1,NAtoms
      If(I1(I).EQ.0) Then
        num3 = num3+1
        I3(num3) = I
      EndIf
 40   CONTINUE
C
C  now find the shortest distance between the atoms in the
C  two groups and connect them
C
      DistM = Big
      DO 50 inum=1,num
      I = I2(inum)
      DO 45 jnum=1,num3
      J = I3(jnum)
C
C  get interatomic distance I-J
C
      Dist = SQRT( (XC(1,I) - XC(1,J))**2 +
     $             (XC(2,I) - XC(2,J))**2 +
     $             (XC(3,I) - XC(3,J))**2 )
c
      If(Dist.LT.DistM) Then
        DistM = Dist
        IMin = I
        JMin = j
      EndIf
c
 45   CONTINUE
 50   CONTINUE
C
C  connect the new pair
      IC(IMin,JMin) = 1
      IC(JMin,IMin) = 1
C
C  check for other atoms (most likely symmetry-related) that
C  are the same distance as (IMin,JMin) pair and connect them
C
      DO 60 inum=1,num
      I = I2(inum)
      DO 55 jnum=1,num3
      J = I3(jnum)
C
C  get interatomic distance I-J
C
      Dist = SQRT( (XC(1,I) - XC(1,J))**2 +
     $             (XC(2,I) - XC(2,J))**2 +
     $             (XC(3,I) - XC(3,J))**2 )
c
      If(Dist.EQ.DistM) Then
        IC(I,J) = 1
        IC(J,I) = 1
      EndIf
c
 55   CONTINUE
 60   CONTINUE
C
C  now go back and continue the search
C
      num = num+1
      I2(num) = JMin
      I1(JMin) = 1
      GO TO 20
C
      END
c =====================================================================
c
      SUBROUTINE ChkCnnct(File,NAtoms,MaxC,IC)
      IMPLICIT INTEGER(A-Z)
C
C  Checks if there is any additional connectivity data
C
C  ARGUMENTS
C
C  File    -  name of options file
C  NAtoms  -  number of atoms
C  MaxC    -  maximum atomic connectivity
C  IC      -  connectivity matrix
C
      INTEGER IC(NAtoms,NAtoms)
c ..................................................
c -- automatic allocation of arrays in F90
      INTEGER IAC(MaxC)
c ..................................................
      CHARACTER*256 File,CHAR
C
C  initialize
C
      NCnnct = 0
      call rmblan(File,256,Len)
C
C  open options file
C
      OPEN (UNIT=40,FILE=File(1:Len),FORM='FORMATTED',STATUS='OLD',
     $      ERR=95)
C
C  is there a connectivity section?

 10   CONTINUE
      READ(40,900,END=30) CHAR
      If(CHAR(1:8).NE.'$connect') GO TO 10
C
C  there is atom connectivity data
C  determine how many lines
C
 20    CONTINUE
       READ(40,900,END=30) CHAR
       If(CHAR(1:8).EQ.'$endconn') GO TO 30
c
       NCnnct = NCnnct+1
       GO TO 20
C
 30   CONTINUE
      If(NCnnct.EQ.0) Then
       CLOSE (UNIT=40,STATUS='KEEP')
       RETURN
      EndIf
C
C  there is connectivity data - read it
C  EXPECTED FORMAT
C   atom number     list of connected atoms
C
      REWIND 40
      READ(40,900) CHAR
c
      DO 50 I=1,NCnnct
c
      CALL IZeroIT(IAC,MaxC)
      READ(40,*,ERR=35,END=35) IAtom,(IAC(L),L=1,MaxC)
c
  35  continue
      DO 40 J=1,MaxC
      JJ = IAC(J)
      If(JJ.EQ.0) GO TO 50
c
      If(IAtom.GT.NAtoms.OR.JJ.GT.NAtoms) Then
       Call nerror(37,'OPTIMIZE module',
     $ 'Connectivity Data (Line # given) Contains too Large Entry',
     $  I,NAtoms)
      EndIF
c
      IC(IAtom,JJ) = 1
      IC(JJ,IAtom) = 1
 40   CONTINUE
c
 50   CONTINUE
c
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  All data read
C  check that no atom is connected to itself
C
      DO 60 I=1,NAtoms
      IF(IC(I,I).NE.0) THEN
       Call nerror(38,'OPTIMIZE module',
     $  'Atom (# given) is Connected to itself!  Check input data',
     $   I,0)
      EndIF
 60   CONTINUE
C
C  If we get here, we're OK
C
      RETURN
C
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no options file found
      Call nerror(39,'OPTIMIZE module',
     $  'There Should be an options file but none found!',
     $   0,0)
 97   CONTINUE          ! error reading connectivity data
      Call nerror(41,'OPTIMIZE module',
     $      'Problem Reading Atomic Connectivity',IAtom,0)
c
  900 Format(A80)
  910 Format(I4,2X,8I4)
c
      END
c =====================================================================
c
      SUBROUTINE ChkD(N,DMAX,IPRNT,D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Checks the length of vector D
C  If |D| > DMAX, scales D down accordingly
C
      REAL*8 D(N)
C
C
      IOut = ioutfil('iout')
      DD = SQRT(SProd(N,D,D))
C
      IF(DD.GT.DMAX) THEN
       skal = DMAX/DD
       If(IPRNT.GT.1) WRITE(IOut,1000) skal
       CALL Vscal(N,skal,D)
       DD = DMAX
      ENDIF
C
      If(IPRNT.GT.1) WRITE(IOut,1100) DD
      RETURN
c
 1000 FORMAT(' Calculated Step too Large.  Step scaled by ',F9.6)
 1100 FORMAT(' Step Taken.  Stepsize is ',F9.6)
c
      END
c =====================================================================
c
      SUBROUTINE ChkDIH(RCon,Dih)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Checks for possible dihedral angle changes near 180 degrees
C
C  ARGUMENTS
C
C  RCon    -  desired dihedral angle constraint
C  Dih     -  current value
C
C
      PARAMETER (ZERO=0.0d0)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
C
      IF(RCon.GT.ZERO.AND.Dih.LT.ZERO.OR.
     $   RCon.LT.ZERO.AND.Dih.GT.ZERO) THEN
c
       d = Abs(Dih) - PI
       d = -SIGN((PI-d),Dih)
c
       If(Abs(RCon-d).LT.Abs(RCon-Dih)) Dih = d
c
      ENDIF
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE ChkDUM(NAtoms, NCons,  IDum,   XC,     ICTYP,
     $                  RCON,   IC,     Switch)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks that dummy atoms are positioned appropriately for
C  correct definition of near-limiting dihedral constraints
C
C  ARGUMENTS
C
C  NAtoms  -  number of real atoms
C  NCons   -  number of constraints
C  IDum    -  designated centre number of dummy atom
C  XC      -  Cartesian coordinates for all real atoms
C  ICTYP    -  constraint type
C               1 - fixed distance
C               2 - fixed bond angle
C               3 - fixed out-of-plane bend
C               4 - fixed dihedral angle
C  RCON     -  constraint values
C  IC       -  atoms involved in constraint
C                IC1-IC2           distance constraint
C                IC1-IC2-IC3       bond angle constraint
C                IC1-IC2-IC3-IC4   all other constraints
C  Switch  -  Logical flag assigning dummy atom position
C             on exit      .true.  - switch dummy atom
C                          .false. - leave as is
C
C
      REAL*8 XC(3,NAtoms),RCON(NCons),Ninety
      DIMENSION ICTYP(NCons),IC(4,NCons)
      LOGICAL Switch
C
      PARAMETER (ToANG=180.0d0/3.14159265358979323844d0)
      PARAMETER (TolZERO=1.0d-5,Ninety=90.0d0)
C
C
      Switch = .FALSE.
c
      DO 10 IQ=1,NCons
C
C  Look for dihedral constraints involving dummy atom IDum
C  (there should be two paired constraints)
C
      If(ICTYP(IQ).NE.4) GO TO 10
C
C  does this constraint involve atom IDum?
C  (IDum should the last atom in the first constraint
C   and the first atom in the second and one of the pair
C   should be either +90.0 or -90.0 degrees)
C
      DihVal = Abs(RCON(IQ))*ToANG
c
      IF( (IC(4,IQ).EQ.IDum.OR.IC(1,IQ).EQ.IDum).AND.
     $     Abs(DihVal-Ninety).LT.TolZERO) THEN
cc
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       L = IC(4,IQ)
       CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       val = RCON(IQ)*ToANG
       Th = Th*ToANG
C
C  Val and Th should be the same
C  i.e. either both +90.0 or both -90.0
C  If there is a sign incompatibility, then set Switch
C
       If(Abs(val-Th).GT.Ninety) Switch = .TRUE.
cc
      ENDIF
c
 10   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE ChkHES(N,      NCons,  IPRNT,  U,      EigVal,
     $                  Symflag,G,      thrsh,  NCon,   NC )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Removes eigenvectors with "zero" eigenvalues from a set of
C  N Hessian modes. Additionally, if symmetry is switched on,
C  removes all symmetry-breaking modes including, if applicable,
C  constraint modes
C
C  ARGUMENTS
C
C  N       -  number of modes to check
C             (i.e. total dimension of problem)
C  NCons   -  initial number of constraints
C  IPRNT   -  controls level of printout
C  U       -  eigenvectors stored as columns
C  EigVal  -  eigenvalues (assumed to be ordered)
C  Symflag -  Logical flag determining whether symmetry
C             checking is on or not
C  G       -  current gradient vector
C             (if Symflag is .true. all modes whose overlap
C              with the gradient is less than thrsh will be
C              rejected)
C  thrsh   -  threshold below which eigenvalues will be
C             rejected as "zero"
C  NCon    -  number of constraints remaining
C  NC      -  number of modes remaining
C
C
      REAL*8 U(N,N),EigVal(N),G(N)
      LOGICAL Symflag
C
      PARAMETER (Zero=0.0d0)
      PARAMETER (CTol=1.0d-8)
C
C
      IOut = ioutfil('iout')
      NAT3 = N-NCons
      If(NCons.GT.0) NCon = NCons
C
C  first remove zeroes
C
      IMIN = 0
      DO 10 I=1,N
      IF(EigVal(I).LT.Zero.AND.Abs(EigVal(I)).GT.thrsh) THEN
       IMIN = I
       GO TO 10
      ENDIF
      If(Abs(EigVal(I)).GT.thrsh) GO TO 20
 10   CONTINUE
C
C  should never get here accept for di- and maybe tri-atomics
C
      IF(N.GT.9) THEN
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
 20   CONTINUE
      IMAX = I-1
C
      IF(IMIN.GE.IMAX) THEN
cc
C  there are no zero eigenvalues
C
       NC = N
cc
      ELSE
cc
C  All eigenvalues between IMIN+1 and IMAX to be deleted
C
       DO 30 I=1,N-IMAX
       EigVal(IMIN+I) = EigVal(IMAX+I)
       CALL CpyVEC(N,U(1,IMAX+I),U(1,IMIN+I))
 30    CONTINUE
C
       NC = N - (IMAX-IMIN)
cc
      ENDIF
C
C  "zero" eigenvalues removed
C  now check for symmetry
C
      IF(Symflag) THEN
C
C  systematically remove all modes that break symmetry
C  (whose overlap with the gradient is less than thrsh)
C
       IT = 0
 40    CONTINUE
       IT = IT+1
       If(IT.GT.NC) GO TO 70
C
C  form scalar product with gradient
C  ** WARNING:  Look into normalizing gradient to unity to prevent  **
C  ** modes being removed near convergence when gradient very small **
C
       GDV = SProd(NAT3,G,U(1,IT))
       IF(Abs(GDV).LT.thrsh) THEN
C
C  mode has no overlap with gradient - remove
C
        NC = NC-1
C
C  if we are doing a constrained optimization and the eigenvalue
C  is negative, check if a constraint mode is being removed
C  **  WARNING  **  THIS CHECK IS NOT BULLET-PROOF
C
        IF(NCons.GT.0.AND.EigVal(IT).LT.Zero) THEN
C
C  check if the mode has any constraint component
C
         DO 50 I=NAT3+1,N
         IF(Abs(U(I,IT)).GT.CTol) THEN
          NCon = NCon - 1
          GO TO 51
         ENDIF
 50      CONTINUE
 51      CONTINUE
        ENDIF
C
C  shift all eigenvalues and vectors down
C
        DO 60 J=IT,NC
        EigVal(J) = EigVal(J+1)
        CALL CpyVEC(N,U(1,J+1),U(1,J))
 60     CONTINUE
c
        IT = IT-1
       ENDIF
c
       GO TO 40
 70    CONTINUE
c
      ENDIF
C
C  NC is the total number of modes remaining
C
      IF(IPRNT.GT.1) THEN
       WRITE(IOut,1100) NC
       WRITE(IOut,1200)
       WRITE(IOut,1300) (EigVal(I),I=1,NC)
       If(NCons.GT.NCon) WRITE(IOut,1400) NCons-NCon
      ENDIF
      IF(IPRNT.GT.3) THEN
       WRITE(IOut,1500)
       CALL PrntMAT(NC,N,NC,U)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Hessian Appears to have all zero or',
     $            ' negative eigenvalues')
 1100 FORMAT(/,I3,' Hessian modes will be used to form the next step')
 1200 FORMAT('  Hessian Eigenvalues:')
 1300 FORMAT(1X,6F12.6)
 1400 FORMAT(/,I3,' constraint modes removed due to symmetry')
 1500 FORMAT(/,'  Hessian Eigenvectors:')
c
      END
c =====================================================================
c
      SUBROUTINE ChkHssINT(NAtoms, IAN,    XC,     NPrim,  NDEG,
     $                     HESS,   IPRNT,  ktyp,   klist,  UT,
     $                     ILST,   IMAP,   HP,     HPRIM,  HINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks if Cartesian Hessian used to generate internal
C  coordinate Hessian is only partial and, if so, uses
C  default diagonal primitive force constants and reforms
C  internal coordinate Hessian
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C  NPrim   -  number of primitive internals
C  NDEG    -  number of internal degrees of freedom
C  HESS    -  Hessian in Cartesian coordinates
C  IPRNT   -  print flag
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C                7 - 1/R distance coordinates
C  klist   -  list of atoms involved in each primitive internal
C  UT      -  full set of delocalized internal coordinates
C  ILST    -  list of active atoms in partial Cartesian Hessian
C  IMAP    -  list of active primitives
C  HP      -  storage for partial transformed Hessian
C  HPRIM   -  primitive Hessian
C  HINT    -  final Hessian over delocalized internal coordinates
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),HESS(3*NAtoms,3*NAtoms),
     $          ktyp(NPrim),klist(4,NPrim),UT(NPrim,NDEG),
     $          ILST(NAtoms),IMAP(NPrim),HP(NPrim,NDEG),
     $          HPRIM(NPrim,NPrim),HINT(NDEG,NDEG)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,TolZero=1.0d-8)
C
C
C  Initialize
C
      NP = 0
      CALL IZeroIT(ILST,NAtoms)
      CALL IZeroIT(IMAP,NPrim)
C
C  First determine if we have a partial Cartesian Hessian
C  This is detected by the presence of 3x3 unit diagonal blocks
C
      DO 10 I=1,NAtoms
      II = 3*(I-1) + 1
      H1 = ABS(One-HESS(II,II))
      H2 = ABS(One-HESS(II+1,II+1))
      H3 = ABS(One-HESS(II+2,II+2))
      If(H1.LT.TolZero.AND.H2.LT.TolZero.AND.H3.LT.TolZero) Then
       NP = NP+1
       ILST(NP) = I
      EndIf
 10   CONTINUE
cc      write(6,*) ' ILST array is:'
cc      do i=1,np
cc      write(6,*) i,'  ',ilst(i)
cc      enddo
C
C  At this point we have NP inactive atoms stored in ILST
C
      If(NP.EQ.0) RETURN
c
      If(IPRNT.GT.1) WRITE(6,1000)
C
C  Determine which primitives are inactive
C  If a primitive has any involved atoms inactive, then that
C  primitive is inactive and a default force constant will be used
C
      DO 30 I=1,NPrim
C
       DO 21 J=1,NP
       If(klist(1,I).EQ.ILST(J)) GO TO 30
 21    CONTINUE
       DO 22 J=1,NP
       If(klist(2,I).EQ.ILST(J)) GO TO 30
 22    CONTINUE
       If(ktyp(I).EQ.1) Then     ! active stretch
        IMAP(I) = 1
        GO TO 30
       EndIf
cc
       DO 23 J=1,NP
       If(klist(3,I).EQ.ILST(J)) GO TO 30
 23    CONTINUE
       If(ktyp(I).EQ.2) Then     ! active bend
        IMAP(I) = 1
        GO TO 30
       EndIf
cc
       DO 24 J=1,NP
       If(klist(4,I).EQ.ILST(J)) GO TO 30
 24    CONTINUE
       IMAP(I) = 1               ! active 4-atom internal
C
 30   CONTINUE
cc      write(6,*) ' Mapping Array is:'
cc      do i=1,nprim
cc      write(6,*) i,'  ',imap(i)
cc      enddo
C
C  At this point the IMAP array has 0 entries for all primitives
C  that will be given a default force constant
C
C  Transform the Hessian to primitive internals
C
C  transform columns
C
      DO 40 I=1,NPrim
      DO 39 J=1,NDEG
      Val = Zero
      DO 38 K=1,NDEG
      Val = Val + UT(I,K)*HINT(K,J)
 38   CONTINUE
      HP(I,J) = Val
 39   CONTINUE
 40   CONTINUE
cc      write(6,*) ' Partially transformed Hessian is:'
cc      call prntmat(ndeg,nprim,ndeg,hp)
C
C  transform rows
C  at the same time include the default force constants
C
      CALL ZeroIT(HPRIM,NPrim*NPrim)
      DO 50 I=1,NPrim
      IF(IMAP(I).EQ.0) THEN
       HPRIM(I,I) =  ForceCnst(ktyp(I),klist(1,I),klist(2,I),klist(3,I),
     $                         klist(4,I),0,IAN,XC)
      ELSE
       DO 49 J=1,NPrim
       If(IMAP(J).EQ.0) GO TO 49
       Val = Zero
       DO 48 K=1,NDEG
       Val = Val + HP(I,K)*UT(J,K)
 48    CONTINUE
       HPRIM(I,J) = Val
 49    CONTINUE
      ENDIF
 50   CONTINUE
c
      If(IPRNT.GT.5) Then
       write(6,*) ' Primitive Hessian is:'
       call prntmat(nprim,nprim,nprim,hprim)
      EndIf
C
C  now transform back to non-redundant coordinates
C
      CALL TrnHPRIM(NDEG,NPrim,HPRIM,UT,HP,HINT)
cc      write(6,*) ' Non-Redundant Hessian is'
cc      call prntmat(ndeg,ndeg,ndeg,hint)
C
      RETURN
c
 1000 FORMAT(' Partial Cartesian Hessian only - using default',
     $       ' primitive force constants')
c
      END
c =====================================================================
c
      SUBROUTINE ChkMAG(N,EigMIN,EigMAX,EigVAL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks absolute magnitudes of a set of N eigenvalues.
C  All eigenvalues below EigMIN in magnitude are increased
C  to this value while all those above EigMAX are reduced
C  to this value.
C  IMPORTANT:  It is assumed that the eigenvalues are ordered.
C
C
      REAL*8 EigVAL(N)
C
      PARAMETER (Zero=0.0d0)
C
C
      IOut = ioutfil('iout')
C
C  First check if magnitudes are below EigMIN
C
      DO 10 I=1,N
      If(EigVAL(I).GT.EigMIN) GO TO 20
      IF(Abs(EigVAL(I)).LT.EigMIN) THEN
       If(EigVAL(I).LT.Zero) EigVAL(I) = -EigMIN
       If(EigVAL(I).GT.Zero) EigVAL(I) =  EigMIN
       WRITE(IOut,1000) I,EigVAL(I)
      ENDIF
 10   CONTINUE
      RETURN
c
 20   CONTINUE
C
C  Now check if magnitudes exceed EigMAX
C
      DO 30 I=1,N
      If(Abs(EigVAL(I)).LT.EigMAX) GO TO 40
      If(EigVAL(I).LT.Zero) EigVAL(I) = -EigMAX
      If(EigVAL(I).GT.Zero) EigVAL(I) =  EigMAX
      WRITE(IOut,1100) I,EigVAL(I)
 30   CONTINUE
      RETURN
c
 40   CONTINUE
      DO 50 I=N,1,-1
      If(EigVAL(I).LT.EigMAX) RETURN
      EigVAL(I) = EigMAX
      WRITE(IOut,1100) I,EigVAL(I)
 50   CONTINUE
C
      RETURN
c
 1000 FORMAT('**WARNING** Magnitude of eigenvalue ',I3,' too small.',
     $       '  Replaced by ',F12.6)
 1100 FORMAT('**WARNING** Magnitude of eigenvalue ',I3,' too large.',
     $       '  Replaced by ',F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE ChkOPT(NCNTR,  NAtoms, NDum,   IOptC,  IType,
     $                  ICons,  NCons,  ICTYP,  IC,     RCON,
     $                  LThrsh, NFix,   IFix,   MaxL,   IFunc,
     $                  MaxDiis,TSflag, Symflag)
      IMPLICIT INTEGER(A-Z)
C
C
C  Checks for errors and incompatible optimization options
C
C  ARGUMENTS
C
C  NCNTR    -  total number of atomic centres (real and dummy)
C  NAtoms   -  number of real atoms
C  NDum     -  number of dummy atoms
C  IOptC    -  coordinate type flag
C               0 - optimize in Cartesian coordinates
C               1 - generate and optimize in internal coordinates
C                   if this fails abort
C              -1 - generate and optimize in internal coordinates
C                   if this fails during any stage of the optimization
C                   switch to Cartesians and continue
C               2 - optimize in Z-Matrix internal coordinates
C                   if this fails abort
C              -2 - optimize in Z-Matrix internal coordinates
C                   if this fails during any stage of the optimization
C                   switch to Cartesians and continue
C  IType    -  internal coordinate optimization type flag
C               0 - standard ab initio/semiempirical optimization
C               1 - surface adsorption/reaction using delocalized internals
C               2 - cluster coordinates (delocalized/inverse-distance)
C               3 - distance-only coordinates
C  ICons    -  constraint flag for Cartesian optimization
C               0 - full optimization (no constraints)
C               1 - constrained optimization using Lagrange multipliers
C               2 - constrained optimization using Penalty functions
C              -1 - attempt a constrained optimization using
C                   Lagrange multipliers; if this fails switch
C                   to Penalty functions. At convergence, tidy up
C                   using Lagrange multipliers if this failed originally
C              -2 - constrained optimization using Penalty functions
C                   tidy up at convergence using Lagrange multipliers
C  NCons    -  number of constraints
C  ICTYP    -  constraint type
C               1 - fixed distance
C               2 - fixed bond angle
C               3 - fixed out-of-plane bend
C               4 - fixed dihedral angle
C               5 - fixed linear coplanar bend
C               6 - fixed linear perpendicular bend
C               7 - fixed 1/R inverse distance
C               9 - composite constraint
C  IC       -  atoms involved in constraint
C                IC1-IC2           distance constraint
C                IC1-IC2-IC3       bond angle constraint
C                IC1-IC2-IC3-IC4   all other constraints
C  RCON     -  constraint values
C  LThrsh   -  threshold for near-linear bond angle constraints
C  NFix     -  number of fixed Cartesian coordinates
C  IFix     -  list of fixed/active coordinates
C               0 - coordinate active
C               1 - coordinate will be fixed during optimization
C  MaxL     -  maximum number of real atoms involved in definition
C              of dummy atom position
C  IFunc    -  list for each dummy atom of the real atoms (if any)
C              defining its position
C  MaxDiis  -  maximum size of GDIIS subspace
C  TSflag   -  Logical flag   .true.  - TS search
C           -                 .false. - minimization
C  Symflag  -  Logical flag   .true.  - use symmetry
C                             .false. - do not use symmetry
C
      DIMENSION ICTYP(NCons),IC(4,NCons),IFix(3,NAtoms),IFunc(NDum,MaxL)
      REAL*8 RCON(NCons),TORAD,PThrsh
      LOGICAL TSflag,Symflag,Error
C
      PARAMETER (TORAD=3.14159265358979323844d0/180.0d0)
C
C
      IOut = ioutfil('iout')
      Error = .FALSE.
      PThrsh = DFloat(LThrsh)*TORAD
C
C  Is GDIIS switched on for TS search?
C
      IF(TSflag.AND.MaxDiis.NE.0) THEN
       WRITE(IOut,1100)
       MaxDiis = 0
      ENDIF
C
C  Cannot use internal coordinates with fixed Cartesian coordinates
C  unless ALL coordinates of atom are fixed
C
      IF(NFix.NE.0.AND.IOptC.NE.0) THEN
       DO I=1,NAtoms
       MFix = IFix(1,I) + IFix(2,I) + IFix(3,I)
       If(MFix.EQ.1.OR.MFix.EQ.2) Then
        If(IOptC.EQ.-1) Then
         IOptC = 0
        Else
         WRITE(IOut,1150)
         Error = .TRUE.
        EndIf
        GO TO 5
       EndIf
       EndDO
      ENDIF
c
  5   CONTINUE
C
C  We're not allowing GDIIS with fixed atoms
C
      IF(NFix.NE.0.AND.IOptC.EQ.0.AND.MaxDiis.NE.0) THEN
       WRITE(IOut,1200)
       MaxDiis = 0
      ENDIF
C
C  Must use symmetry with internals
C
      IF(.NOT.Symflag.AND.Abs(IOptC).EQ.1) THEN
       WRITE(IOut,1300)
       Error = .TRUE.
      ENDIF
C
C  Constrained TS search using Penalty functions?
C
      IF(TSflag.AND.Abs(ICons).EQ.2) THEN
       WRITE(IOut,1400)
       Error = .TRUE.
      ENDIF
C
C  Dummy atoms but no constraints?
C
      IF(NDum.GT.0.AND.NCons.EQ.0) THEN
       WRITE(IOut,1500)
       Error = .TRUE.
      ENDIF
C
C  Check dummy atom dependencies
C
      DO 20 I=1,NDum
      DO 10 J=1,MaxL
      If(IFunc(I,J).EQ.0) GO TO 20
      IF(IFunc(I,J).LT.1.OR.IFunc(I,J).GT.NATOMS) THEN
       WRITE(IOut,1600) I,IFunc(I,J)
       Error = .TRUE.
      ENDIF
 10   CONTINUE
 20   CONTINUE
C
C  Check constraints
C
      If(NCons.EQ.0) ICons=0
c
      DO 70 I=1,NCons
C
C  if optimization in Cartesian coordinates
C  cannot handle out-of-plane and colinear bends
C
      IF(IOptC.EQ.0.AND.(ICTYP(I).EQ.3.OR.ICTYP(I).EQ.5.OR.
     $                   ICTYP(I).EQ.6) ) THEN
       WRITE(IOut,1700) I
       Error = .TRUE.
      ENDIF
C
C  if optimization in delocalized internal coordinates
C  best if near-linear bend constraint replaced by equivalent
C  colinear bend constraint
C
      IF(Abs(IOptC).EQ.1.AND.ICTYP(I).EQ.2.AND.
     $                       RCON(I).GT.PThrsh) THEN
       WRITE(IOut,1800) I,DFloat(LThrsh)
       Error = .TRUE.
      ENDIF
C
C  do constraints involve non-existent centres?
C
      IF(IC(1,I).GT.NCNTR.OR.IC(2,I).GT.NCNTR.OR.
     $   IC(3,I).GT.NCNTR.OR.IC(4,I).GT.NCNTR) THEN
       WRITE(IOut,1900) I
       Error = .TRUE.
      ENDIF
C
C  check constraints that involve dummy atoms
C  constraints should
C    (i)  involve ONE dummy atom only
C    (ii) involve at least one of the atoms defining the
C         position of the dummy atom
C
      num = 0
      DO 30 J=1,4
      IF(IC(J,I).GT.NAtoms) Then
       num = num+1
       idum = IC(J,I)
      ENDIF
 30   CONTINUE
c
      IF(num.GT.1) THEN
       WRITE(IOut,2000) I
       Error = .TRUE.
      ENDIF
c
      IF(num.EQ.1) THEN
C
C  constraint involves a dummy atom
C  see if it also involves an atom in the dummy defining list
C
       idum = idum-NAtoms
       DO 50 J=2,MaxL
       If(IFunc(idum,J).EQ.0) GO TO 50
       DO 40 K=1,4
       If(IC(K,I).EQ.IFunc(idum,J)) GO TO 60
 40    CONTINUE
 50    CONTINUE
C
C  If we get here, it means constraint does NOT involve
C  a dummy defining atom
C
       WRITE(IOut,2100) I
       Error = .TRUE.
 60    CONTINUE
c
       If(IFunc(idum,1).EQ.1.AND..NOT.Error) WRITE(IOut,2200) I
c
      ENDIF
C
 70   CONTINUE
C
C
      If(Error) CALL OptExit(9)
C
      RETURN
c
 1100 FORMAT('**WARNING** DIIS not applicable for TS search',/,
     $       '  GDIIS Switched Off')
 1150 FORMAT(/,2X,'***ERROR*** Cannot Fix Cartesians with',
     $              ' Internal Coordinates')
 1200 FORMAT('**WARNING** DIIS not available with Fixed Atoms',/,
     $       '  GDIIS Switched Off')
 1300 FORMAT(/,2X,'***ERROR*** Cannot Switch Off Symmetry with',
     $              ' Internal Coordinates')
 1400 FORMAT(/,2X,'***ERROR*** Cannot Locate Constrained TS Using',
     $              ' Penalty Functions')
 1500 FORMAT(/,2X,'***ERROR*** There are dummy atoms but no',
     $              ' Constraints!',/,
     $          '     This is Pointless.  Check your input')
 1600 FORMAT(/,2X,'***ERROR*** Dummy atom ',I2,' is defined with',
     $              ' respect to centre ',I2,/,
     $          '     This is not a Real Atom!')
 1700 FORMAT(/,2X,'***ERROR*** Constraint ',I2,' CANNOT be handled',
     $              ' in Cartesian Coordinates',/,
     $          '     (out-of-plane or colinear bend)')
 1800 FORMAT(/,2X,'**WARNING** Angle Constraint ',I2,' is Greater',
     $             ' than Near-Linear Threshold ',F5.1,/,
     $          '    Please Replace by Appropriate Colinear Bend',
     $             ' Constraint')
 1900 FORMAT(/,2X,'***ERROR*** Constraint ',I2,' Involves a',
     $              ' Non-existent Centre!')
 2000 FORMAT(/,2X,'***ERROR*** Constraint ',I2,' Involves more than',
     $              ' One Dummy Atom',/,
     $          '     Such Constraints are NOT Supported')
 2100 FORMAT(/,2X,'***ERROR*** Constraint ',I2,' Which Involves a',
     $              ' Dummy Atom Must also',/,
     $          '     Involve at Least One of the Atoms in the',
     $              ' Dummy Atom Definition')
 2200 FORMAT('**WARNING** Constraint ',I2,' Involves a Non-Fixed',
     $       ' Dummy Atom',/,
     $       '  Such Constraints are NOT Guaranteed to Converge!',/,
     $       '  The Best Way to Define a Constraint with Respect',
     $       ' to a Dummy Atom',/,'  is Via a Z-Matrix')
c
      END
c =====================================================================
c
      SUBROUTINE ChkOSCIL(NCount, NDim,   NCons,  ICTYP,  IC,
     $                    IPRNT,  TolG,   EC,     XC,     XOld,
     $                    D,      GC,     EOS,    GOS,    IOS,
     $                    oscil )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine checks for oscillatory behaviour and modifies
C  the current Lagrange multipliers accordingly
C
C  ARGUMENTS
C
C  NCount  -  number of times this routine entered
C             incremented by 1 on entry
C  NDim    -  total dimension of problem
C  NCons   -  number of constraints
C  ICTYP   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane bend    CURRENTLY NOT IMPLEMENTED
C              4 - fixed dihedral angle
C  IC      -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C  IPRNT   -  controls level of print out
C  TolG    -  convergence criterion on maximum gradient component
C  EC      -  current energy
C  XC      -  current geometry  (including Lagrange multipliers)
C  XOld    -  scratch space for copy of current geometry
C  D       -  previous displacement (including Lagrange multipliers)
C  GC      -  current gradient  (including Lagrange multipliers)
C  EOS     -  storage for current and two previous energies
C  GOS     -  storage for current and two previous gradients
C  IOS     -  oscillation counter
C  oscil   -  Logical flag indicating oscillation status
C             on exit   .true.  - we are oscillating
C                       .false. - we are not oscillating
C
C
      REAL*8 XC(NDim),XOld(NDim),D(NDim),GC(NDim),EOS(3),GOS(NCons,3)
      INTEGER ICTYP(NCons),IC(4,NCons)
      LOGICAL oscil
C
      PARAMETER (Zero=0.0d0)
C
C
      NAT3 = NDim-NCons
C
C  increment counter and store energy and Lagrange gradients
C
      NCount = NCount+1
      EOS(1) = EOS(2)
      EOS(2) = EOS(3)
      EOS(3) = EC
      DO 10 IQ=1,NCons
      GOS(IQ,1) = GOS(IQ,2)
      GOS(IQ,2) = GOS(IQ,3)
      GOS(IQ,3) = GC(NAT3+IQ)
 10   CONTINUE
C
C  If oscillation occured on the previous cycle, exit
C
      IF(oscil) THEN
       oscil = .FALSE.
       RETURN
      ENDIF
C
C  If we haven't done at least 3 cycles, exit
C
      If(NCount.LT.3) RETURN
C
C  Check for oscillatory behaviour
C
      IF( (EOS(1).GT.EOS(2).AND.EOS(3).GT.EOS(2)) .OR.
     $    (EOS(1).LT.EOS(2).AND.EOS(3).LT.EOS(2)) ) THEN
cccc
C  first save current geometry
C
       CALL CpyVEC(NDim,XC,XOld)
C
C  now construct previous geometry in D
C
       DO 20 I=1,NDim
       D(I) = XC(I)-D(I)
 20    CONTINUE
C
C  now check each constraint in turn
C
       DO 60 IQ=1,NCons
C
C  skip if Lagrange gradient below convergence criterion
C
       If(Abs(GOS(IQ,3)).LT.TolG) GO TO 60
c
       IF( (GOS(IQ,1).LT.Zero.AND.GOS(IQ,2).GT.Zero.AND.
     $      GOS(IQ,3).LT.Zero) .OR.
     $     (GOS(IQ,1).GT.Zero.AND.GOS(IQ,2).LT.Zero.AND.
     $      GOS(IQ,3).GT.Zero) ) THEN
ccc
C  estimate zero Lagrange gradient by linear interpolation
C
        oscil = .TRUE.
        A0 = -GOS(IQ,3)/(GOS(IQ,2)-GOS(IQ,3))
        A1 =  GOS(IQ,2)/(GOS(IQ,2)-GOS(IQ,3))
c
        XC(NAT3+IQ) = A0*D(NAT3+IQ) + A1*XC(NAT3+IQ)
C
C  interpolate the positions of all atoms involved in
C  the constraints
C
        IF(ICTYP(IQ).EQ.1) THEN
cc
C  distance constraint
C
         IAT = 3*(IC(1,IQ)-1)
         JAT = 3*(IC(2,IQ)-1)
         DO 30 I=1,3
         II = IAT+I
         JJ = JAT+I
         XC(II) = A0*D(II) + A1*XC(II)
         XC(JJ) = A0*D(JJ) + A1*XC(JJ)
 30      CONTINUE
cc
        ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  bond angle constraint
C
         IAT = 3*(IC(1,IQ)-1)
         JAT = 3*(IC(2,IQ)-1)
         KAT = 3*(IC(3,IQ)-1)
         DO 40 I=1,3
         II = IAT+I
         JJ = JAT+I
         KK = KAT+I
         XC(II) = A0*D(II) + A1*XC(II)
         XC(JJ) = A0*D(JJ) + A1*XC(JJ)
         XC(KK) = A0*D(KK) + A1*XC(KK)
 40      CONTINUE
cc
        ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral angle constraint
C
         IAT = 3*(IC(1,IQ)-1)
         JAT = 3*(IC(2,IQ)-1)
         KAT = 3*(IC(3,IQ)-1)
         LAT = 3*(IC(4,IQ)-1)
         DO 50 I=1,3
         II = IAT+I
         JJ = JAT+I
         KK = KAT+I
         LL = LAT+I
         XC(II) = A0*D(II) + A1*XC(II)
         XC(JJ) = A0*D(JJ) + A1*XC(JJ)
         XC(KK) = A0*D(KK) + A1*XC(KK)
         XC(LL) = A0*D(LL) + A1*XC(LL)
 50      CONTINUE
cc
        ENDIF
ccc
       ENDIF
 60    CONTINUE
C
C  now modify D to reflect geometry changes
C
       DO 70 I=1,NDim
       D(I) = XC(I)-XOld(I)
 70    CONTINUE
cccc
      ENDIF
C
C  increment oscillation counter
C
      IF(oscil) THEN
       IOS = IOS+1
       IOut = ioutfil('iout')
       If(IPRNT.GT.1) WRITE(IOut,1000)
      ENDIF
C
      RETURN
c
 1000 FORMAT('**WARNING** Energy is Oscillating  Using Linear',
     $       ' Interpolation')
c
      END
c =====================================================================
c
      SUBROUTINE ChkROT(RM,IPRNT,thrsh,rotate)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines whether any rotation of the coordinate
C  system occured following e.g. symmetry checking.
C  On exit the logical flag rotate is set .true. if
C  RM is anything other than a unit matrix
C
C  ARGUMENTS
C
C  RM      -  rotation matrix
C  IPRNT   -  controls level of print out
C  thrsh   -  accuracy threshold
C  rotate  -  logical flag
C             set .false. if RM is a unit matrix
C             i.e. coordinates are identical
C             set .true. otherwise
C
C
      REAL*8 RM(3,3),U(3,3)
      LOGICAL rotate,CmpVEC
C
      DATA U / 1.0d0, 0.0d0, 0.0d0,
     $         0.0d0, 1.0d0, 0.0d0,
     $         0.0d0, 0.0d0, 1.0d0 /
C
      If(IPRNT.GT.6) Then
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       Call PrntMAT(3,3,3,RM)
      EndIf
C
C  Check for Unit Matrix
C
      rotate = CmpVEC(9,RM,U,thrsh)
      if(rotate) then
       rotate = .FALSE.
      else
       rotate = .TRUE.
      endif
C
      RETURN
c
 1000 FORMAT(/,' Rotation Matrix is:')
c
      END
c =====================================================================
c
      SUBROUTINE ChkSCAL(NAT3,IPRNT,EigVal,HESS,scale)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks whether it is likely to be profitable to scale
C  the Hessian in an attempt to calculate a suitable
C  shift parameter for EF during a constrained optimization
C  (The original value for the shift parameter was unacceptable)
C
C  The assumption here is that the lowest (positive) Hessian
C  eigenvalue scales linearly whilst the (initially too high)
C  shift parameter scales at a slower rate; thus it may be
C  possible to find a suitable shift parameter (i.e. one lower
C  than the Hessian eigenvalue) by scaling the Hessian.
C
C  ARGUMENTS
C
C  NAT3       -  dimension of Hessian (3*number of atoms)
C  IPRNT      -  controls level of printout
C  EigVal(1)  -  lowest "positive" Hessian Eigenvalue
C  Eigval(2)  -  previously calculated shift parameter
C                * NOTE *  These values set in CONFormD
C  HESS       -  Hessian matrix
C  scale      -  contains Hessian scale factor on exit
C
C
      REAL*8 HESS(NAT3,NAT3),EigVal(2)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
      PARAMETER (EigMIN=0.001d0,ScalMAX=5.0d0)
C
C
      IOut = ioutfil('iout')
c
      IF(scale.NE.One) THEN
cc
C  already made an attempt at scaling which failed
C  restore original Hessian and exit
C
       CALL Vscal(NAT3*NAT3,One/scale,HESS)
       scale = One
cc
      ELSE
cc
C  calculate likely scale factor
C
       If(EigVal(1).LT.EigMIN.OR.EigVal(2).LT.ZERO) GO TO 95
c
       scale = EigVal(2)/EigVal(1)
       scale = scale+scale
       IF(scale.LE.ScalMAX) THEN
        If(IPRNT.GT.1) WRITE(IOut,1000) scale
        CALL Vscal(NAT3*NAT3,scale,HESS)
       ELSE
        scale = One
        GO TO 95
       ENDIF
cc
      ENDIF
C
      RETURN
C
 95   CONTINUE
      If(IPRNT.GT.2) WRITE(IOut,1100)
      RETURN
c
 1000 FORMAT(' ** Scaling Hessian by ',F9.6,' **',/)
 1100 FORMAT(' ** No Hessian Scaling Attempted **',/)
c
      END
c =====================================================================
c
      SUBROUTINE ChkXSC(NAtoms,XC,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for and removes "excess connectivity"
C  Principally to prevent in the near-linear arrangement A--B--C
C  that if A is connected to B and B is connected to C then
C  A is not also connected to C
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  IC      -  connectivity matrix (may be modified on exit)
C
C
      DIMENSION XC(3,NAtoms),IC(NAtoms,NAtoms)
c
      PARAMETER (TollZero=0.002d0)
C
C
C  First find a connectivity
C
      DO 20 I=2,NAtoms
      DO 20 J=1,I
      IF(IC(I,J).EQ.1) THEN
C
C  I and J are connected
C  What else is J connected to?
C
        DO 10 K=J+1,I
        IF(IC(K,J).EQ.1) THEN
C
C  J and K are connected
C  Is K also connected to I?
C
          If(IC(I,K).EQ.1) Then
C
C  check the angle I-J-K
C  if near-linear, remove longest bond
C
            CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
            If(Th.LT.TollZero) Then
              DIJ = SQRT( (XC(1,I) - XC(1,J))**2 +
     $                    (XC(2,I) - XC(2,J))**2 +
     $                    (XC(3,I) - XC(3,J))**2 )
              DIK = SQRT( (XC(1,I) - XC(1,K))**2 +
     $                    (XC(2,I) - XC(2,K))**2 +
     $                    (XC(3,I) - XC(3,K))**2 )
              If(DIJ.GT.DIK) Then
                IC(I,J) = 0
                IC(J,I) = 0
              Else
                IC(I,K) = 0
                IC(K,I) = 0
              EndIf
            EndIf
          EndIf
        ENDIF
 10     CONTINUE
cc
      ENDIF
 20   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE CnvBACK(NCycle, NIC,    XINT,   QQ,     thrsh,
     $                   IPRNT,  QSQR,   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for convergence during iterative back-transformation
C  of natural internal coordinates to Cartesians.
C  Convergence is attained when each internal coordinate
C  generated from the current estimate of the corresponding
C  Cartesian coordinates differs by less than thrsh from its
C  actual (known) value.
C
C
C  ARGUMENTS
C
C  NCycle  -  current cycle number
C  NIC     -  number of internal coordinates
C  XINT    -  actual internal coordinate values
C  QQ      -  internals as generated from current estimate
C             of Cartesians
C  thrsh   -  convergence threshold
C  IPRNT   -  flag controlling printout
C  QSQR    -  sum of squares of total error
C  Cnvgd   -  logical flag indicating convergence
C             on exit    .true.  - convergence achieved
C                        .false. - another step needed
C
C
      REAL*8 XINT(NIC),QQ(NIC)
      LOGICAL Cnvgd
C
      PARAMETER (Zero=0.0d0)
C
C
      IOut = ioutfil('iout')
c
      If(IPRNT.GT.4) WRITE(IOut,1000)
C
      DMax = Abs(XINT(1)-QQ(1))
      QSQR = Zero
c
      DO 10 I=1,NIC
      DX = XINT(I)-QQ(I)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
      If(IPRNT.GT.4) WRITE(IOut,1100) I,XINT(I),QQ(I),DX
c
      QSQR = QSQR + DX**2
 10   CONTINUE
C
      QSQR = SQRT(QSQR/NIC)
C
      If(IPRNT.GT.2) WRITE(IOut,1200) NCycle,DMax,QSQR
C
C  Converged?
C
      Cnvgd = DMax.LT.thrsh
C
      RETURN
c
 1000 FORMAT(/,5X,'Iterative generation of Cartesian Coordinates',/,
     $       '     Internal   True Value     Estimate     Difference')
 1100 FORMAT(5X,I5,2X,3(2X,F12.8))
 1200 FORMAT(5X,'Cycle: ',I3,'  Maximum deviation: ',F12.8,
     $          '  RMS deviation: ',F12.8)
c
      END
c =====================================================================
c
      SUBROUTINE CnvCART(NAtoms, AtSymb, XC,     D,      GC,
     $                   EC,     EOld,   IPRNT,  TolG,   TolD,
     $                   TolE,   NEG,    NegReq, Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for convergence during geometry optimization in Cartesian coordinates
C  The flag Cnvgd is returned .true. if the maximum component
C  of the gradient vector, GC, is below TolG AND either the
C  maximum component of the displacement vector, D, is below
C  TolD OR the energy change is below TolE.
C  (TolG, TolD and TolE are currently in atomic units)
C  Additionally the Hessian must have the correct
C  eigenvalue structure.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  current Cartesian coordinates
C  D       -  displacement vector
C  GC      -  current gradient
C  EC      -  current energy
C  EOld    -  previous energy
C  IPRNT   -  print flag
C  TolG    -  convergence criterion on maximum gradient component
C  TolD    -  convergence criterion on maximum component of
C             displacement vector
C  TolE    -  convergence criterion on maximum energy change
C             from previous cycle
C  NEG     -  number of negative Hessian eigenvalues
C  NegReq  -  number of negative Hseeian eigenvalues required
C  Cnvgd   -  logical flag     .true. - converged
C                             .false. - not converged
C
      REAL*8 XC(3,NAtoms),D(3,NAtoms),GC(3,NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
      LOGICAL Cnvgd
C
      CHARACTER*3 DCnvgd,GCnvgd,ECnvgd
      CHARACTER*1 C(3)
      DATA C(1)/'X'/, C(2)/'Y'/, C(3)/'Z'/
C
C
      IOut = ioutfil('iout')
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
      GMax = Abs(GC(1,1))
      DMax = Abs(D(1,1))
c
      DO 10 IAtm=1,NATOMS
      DO 10 J=1,3
      GX = GC(J,IAtm)
      DX = D(J,IAtm)
      If(Abs(GX).GT.GMax) GMax = Abs(GX)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
      IF(IPRNT.GT.2) THEN
       CX = XC(J,IAtm)
       CN = CX + DX
       If(J.EQ.1) WRITE(IOut,1100) IAtm,AtSymb(IAtm),C(J),CX,GX,DX,CN
       If(J.NE.1) WRITE(IOut,1200)                   C(J),CX,GX,DX,CN
      ENDIF
c
 10   CONTINUE
c
      Edif = EC - EOld
c
      IF(IPRNT.GT.1) THEN
       GCnvgd = ' NO'
       DCnvgd = ' NO'
       ECnvgd = ' NO'
       If(GMax.LT.TolG) GCnvgd = 'YES'
       If(DMax.LT.TolD) DCnvgd = 'YES'
       If(Abs(Edif).LT.TolE) ECnvgd = 'YES'
       WRITE(IOut,1300)
       WRITE(IOut,1400) GMax,TolG,GCnvgd
       WRITE(IOut,1500) DMax,TolD,DCnvgd
       WRITE(IOut,1600) Edif,TolE,ECnvgd
      ENDIF
C
C  Converged?
C
      Cnvgd = GMax.LT.TolG.AND.(DMax.LT.TolD.OR.Abs(Edif).LT.TolE)
     $                    .AND.NEG.EQ.NegReq
C
      RETURN
c
 1000 FORMAT(//,21X,' Coordinates and Displacements in Atomic Units',/,
     $       '   ATOM            Current Value    Gradient  ',
     $       'Displacement   New Value')
 1100 FORMAT(I3,2X,A8,2X,A1,5X,F10.6,3X,F10.6,3X,F10.6,3X,F10.6)
 1200 FORMAT(15X,A1,5X,F10.6,3X,F10.6,3X,F10.6,3X,F10.6)
 1300 FORMAT(/,29X,'Maximum     Tolerance    Cnvgd?')
 1400 FORMAT(9X,'Gradient           ',F8.6,6X,F8.6,5X,A3)
 1500 FORMAT(9X,'Displacement       ',F8.6,6X,F8.6,5X,A3)
 1600 FORMAT(9X,'Energy change     ',F9.6,6X,F8.6,5X,A3,/)
c
      END
c =====================================================================
c
      SUBROUTINE CnvCON(NAtoms, NCNTR,  NCon,   AtSymb, XC,
     $                  D,      GC,     EC,     EOld,   IPRNT,
     $                  TolG,   TolD,   TolE,   NEG,    NegReq,
     $                  Cnvgd )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for convergence during constrained geometry optimization
C  in Cartesian coordinates.
C  The flag Cnvgd is returned .true. if the maximum component
C  of the gradient vector, GC, is below TolG AND either the
C  maximum component of the displacement vector, D, is below
C  TolD OR the energy change is below TolE.
C  (TolG, TolD and TolE are currently in atomic units)
C  Additionally the Hessian must have the correct
C  eigenvalue structure.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NCNTR   -  number of centres
C  NCon    -  number of constraints
C  AtSymb  -  atomic symbols
C  XC      -  current Cartesian coordinates
C  D       -  displacement vector
C  GC      -  current gradient
C  EC      -  current energy
C  EOld    -  previous energy
C  IPRNT   -  print flag
C  TolG    -  convergence criterion on maximum gradient component
C  TolD    -  convergence criterion on maximum component of
C             displacement vector
C  TolE    -  convergence criterion on maximum energy change
C             from previous cycle
C  NEG     -  number of negative Hessian eigenvalues
C  NegReq  -  number of negative Hseeian eigenvalues required
C  Cnvgd   -  logical flag     .true. - converged
C                             .false. - not converged
C

      REAL*8 XC(*),D(*),GC(*)
      CHARACTER*8 AtSymb(NAtoms)
      LOGICAL Cnvgd
C
      CHARACTER*3 DCnvgd,GCnvgd,ECnvgd
      CHARACTER*1 C(3)
      DATA C(1)/'X'/, C(2)/'Y'/, C(3)/'Z'/
C
C
      IOut = ioutfil('iout')
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
      GMax = Abs(GC(1))
      DMax = Abs(D(1))
c
      DO 20 IATM=1,NAtoms
      II = 3*(IATM-1)
      DO 10 J=1,3
      IT = II+J
      GX = GC(IT)
      DX = D(IT)
      If(Abs(GX).GT.GMax) GMax = Abs(GX)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
      IF(IPRNT.GT.2) THEN
       CX = XC(IT)
       CN = CX + DX
       If(J.EQ.1) WRITE(IOut,1100) IATM,AtSymb(IATM),C(J),CX,GX,DX,CN
       If(J.NE.1) WRITE(IOut,1200)                   C(J),CX,GX,DX,CN
      ENDIF
c
 10   CONTINUE
 20   CONTINUE
C
      If(IPRNT.GT.2) WRITE(IOut,1700)
C
      DO 30 J=1,NCon
      IT = 3*NCNTR+J
      GX = GC(IT)
      DX = D(IT)
      If(Abs(GX).GT.GMax) GMax = Abs(GX)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
      IF(IPRNT.GT.2) THEN
       CX = XC(IT)
       CN = CX + DX
       WRITE(IOut,1800) J,CX,GX,DX,CN
      ENDIF
c
 30   CONTINUE
c
      Edif = EC - EOld
c
      IF(IPRNT.GT.1) THEN
       GCnvgd = ' NO'
       DCnvgd = ' NO'
       ECnvgd = ' NO'
       If(GMax.LT.TolG) GCnvgd = 'YES'
       If(DMax.LT.TolD) DCnvgd = 'YES'
       If(Abs(Edif).LT.TolE) ECnvgd = 'YES'
       WRITE(IOut,1300)
       WRITE(IOut,1400) GMax,TolG,GCnvgd
       WRITE(IOut,1500) DMax,TolD,DCnvgd
       WRITE(IOut,1600) Edif,TolE,ECnvgd
      ENDIF
C
C  Converged?
C
      Cnvgd = GMax.LT.TolG.AND.(DMax.LT.TolD.OR.Abs(Edif).LT.TolE)
     $                    .AND.NEG.EQ.NegReq
C
      RETURN
c
 1000 FORMAT(//,21X,' Coordinates and Displacements in Atomic Units',/,
     $       '   ATOM            Current Value    Gradient  ',
     $       'Displacement   New Value')
 1100 FORMAT(I3,2X,A8,2X,A1,5X,F10.6,3X,F10.6,3X,F10.6,3X,F10.6)
 1200 FORMAT(15X,A1,5X,F10.6,3X,F10.6,3X,F10.6,3X,F10.6)
 1300 FORMAT(/,29X,'Maximum     Tolerance    Cnvgd?')
 1400 FORMAT(9X,'Gradient           ',F8.6,6X,F8.6,5X,A3)
 1500 FORMAT(9X,'Displacement       ',F8.6,6X,F8.6,5X,A3)
 1600 FORMAT(9X,'Energy change     ',F9.6,6X,F8.6,5X,A3,/)
 1700 FORMAT(/,21X,' Lagrange Multipliers for Constraints',/,
     $       '       Constraint  Current Value    Gradient  ',
     $       'Displacement   New Value')
 1800 FORMAT(10X,I3,9X,F9.6,4X,F9.6,4X,F9.6,4X,F9.6)
c
      END
c =====================================================================
c
      SUBROUTINE CnvINT(NVar,   XINT,   D,      GC,     EC,
     $                  EOld,   IPRNT,  TolG,   TolD,   TolE,
     $                  NEG,    NegReq, Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for convergence during geometry optimization in Internal coordinates
C  The flag Cnvgd is returned .true. if the maximum component
C  of the gradient vector, GC, is below TolG AND either the
C  maximum component of the displacement vector, D, is below
C  TolD OR the energy change is below TolE.
C  (TolG, TolD and TolE are currently in atomic units)
C  Additionally the Hessian must have the correct
C  eigenvalue structure.
C
C  ARGUMENTS
C
C  NVar    -  number of variables
C  XINT    -  current parameter values
C  D       -  displacement vector
C  GC      -  current gradient
C  EC      -  current energy
C  EOld    -  previous energy
C  IPRNT   -  print flag
C  TolG    -  convergence criterion on maximum gradient component
C  TolD    -  convergence criterion on maximum component of
C             displacement vector
C  TolE    -  convergence criterion on maximum energy change
C             from previous cycle
C  NEG     -  number of negative Hessian eigenvalues
C  NegReq  -  number of negative Hseeian eigenvalues required
C  Cnvgd   -  logical flag     .true. - converged
C                             .false. - not converged
C
      REAL*8 XINT(NVar),D(NVar),GC(NVar)
      LOGICAL Cnvgd
C
      CHARACTER*3 DCnvgd,GCnvgd,ECnvgd
C
C
      IOut = ioutfil('iout')
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
      GMax = Abs(GC(1))
      DMax = Abs(D(1))
c
      DO 10 I=1,NVar
      GX = GC(I)
      DX = D(I)
      If(Abs(GX).GT.GMax) GMax = Abs(GX)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
      IF(IPRNT.GT.2) THEN
       CX = XINT(I)
       CN = CX + DX
       WRITE(IOut,1100) I,CX,GX,DX,CN
      ENDIF
c
 10   CONTINUE
c
      Edif = EC - EOld
c
      IF(IPRNT.GT.1) THEN
       GCnvgd = ' NO'
       DCnvgd = ' NO'
       ECnvgd = ' NO'
       If(GMax.LT.TolG) GCnvgd = 'YES'
       If(DMax.LT.TolD) DCnvgd = 'YES'
       If(Abs(Edif).LT.TolE) ECnvgd = 'YES'
       WRITE(IOut,1300)
       WRITE(IOut,1400) GMax,TolG,GCnvgd
       WRITE(IOut,1500) DMax,TolD,DCnvgd
       WRITE(IOut,1600) Edif,TolE,ECnvgd
      ENDIF
C
C  Converged?
C
      Cnvgd = GMax.LT.TolG.AND.(DMax.LT.TolD.OR.Abs(Edif).LT.TolE)
     $                    .AND.NEG.EQ.NegReq
C
      RETURN
c
 1000 FORMAT(//,5X,'Parameter Values and Displacements in Internal',
     $              ' Coordinates',/,
     $       '  Coordinate   Current Value    Gradient  ',
     $       'Displacement   New Value')
 1100 FORMAT(5X,I3,9X,F10.6,3X,F10.6,3X,F10.6,3X,F10.6)
 1300 FORMAT(/,29X,'Maximum     Tolerance    Cnvgd?')
 1400 FORMAT(9X,'Gradient           ',F8.6,6X,F8.6,5X,A3)
 1500 FORMAT(9X,'Displacement       ',F8.6,6X,F8.6,5X,A3)
 1600 FORMAT(9X,'Energy change     ',F9.6,6X,F8.6,5X,A3,/)
c
      END
c =====================================================================
c
      SUBROUTINE CONFormD(NCycle, NC,     N,      NCons,  NegReq,
     $                    mode,   IPRNT,  U,      EigVal, FX,
     $                    VMODE,  D,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculate the new RFO step for constrained optimization
C  MIN Search:  Forms a step by P-RFO that attempts to maximize
C               along the NCons constraint modes and minimize
C               along all other Hessian modes
C  TS Search:   Forms a step by P-RFO that attempts to maximize
C               along the NCons constraint modes plus one other
C               chosen mode (stored in VMODE) and minimize
C               along all other Hessian modes
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C  NC      -  number of modes used to form new step
C  N       -  actual number of coordinates
C             (i.e. full dimension of problem)
C  NCons   -  number of constraints
C             (for constrained Cartesian coordinate optimization)
C  NegReq  -  number of negative eigenvalues desired
C  mode    -  mode being followed during Transition State search
C             (zero if mode following switched off)
C  IPRNT   -  controls level of print out
C  U       -  Hessian eigenvectors
C  EigVal  -  Hessian eigenvalues
C  FX      -  gradient along each mode
C  VMODE   -  eigenvector followed on previous cycle
C             (if mode is non-zero)
C
C  On exit
C
C  VMODE   -  eigenvector followed this cycle
C             (if mode is non-zero)
C  D       -  new coordinate displacement (next step)
C  IErr    -  error flag    0 - new step calculated successfully
C                          -1 - something went wrong
C
C
      REAL*8 U(N,*),EigVal(NC),FX(NC),VMODE(N),D(N)
C
      REAL*8 Lamda,Lamda1,Lamda2
C
      PARAMETER (ZERO=0.0d0,HALF=0.5d0,FOUR=4.0d0)
      PARAMETER (TOLL=1.0d-8,STEP=0.05d0,BIG=1.0d+3,MAXIT=999)
C
C
      IOut = ioutfil('iout')
      IErr = -1
c
      NUMIT = 0
      IT = NCons
      CALL ZeroIT(D,N)
C
      If(NCons.EQ.0) GO TO 55
C
C  (a)  Maximize along NegReq Hessian modes
C
      IF(mode.NE.0) THEN
cc
       CALL ModeOverlap(NCycle, N,      mode,   Nmode,  IPRNT,
     $                  U,      VMODE)
C
C  On return from ModeOverlap, Nmode is the additional mode
C  along which the energy is to be maximized
C
       If(IPRNT.GT.1.AND.Nmode.NE.mode) WRITE(IOut,1000) mode,Nmode
       mode = Nmode
       If(IPRNT.GT.1) WRITE(IOut,1100) mode
       IT = mode
C
C  If the mode now being followed is now amongst the NCons+1
C  lowest modes, switch off mode following
C
       IF(mode.LE.NegReq) THEN
        mode = 0
        If(IPRNT.GT.1) WRITE(IOut,1200)
       ENDIF
cc
      ELSE
cc
       IF(IPRNT.GT.1) THEN
        If(NegReq.EQ.NCons) WRITE(IOut,1300)
        If(NegReq.GT.NCons) WRITE(IOut,1350)
       ENDIF
       If(NegReq.GT.NCons) IT = NegReq
cc
      ENDIF
C
      IF(NegReq.EQ.1) THEN
cc
C  There is only one constraint and we are minimizing
C
      Lamda = EigVal(IT) + SQRT( EigVal(IT)**2 + FOUR*FX(IT)*FX(IT) )
      Lamda = HALF*Lamda
c
      If(IPRNT.GT.1) WRITE(IOut,1400) Lamda
      IF(EigVal(IT+1).GT.ZERO.AND.Lamda.GT.EigVal(IT+1)) THEN
C
C  Set up for potential Hessian scaling in subroutine ChkSCAL
C
       EigVal(1) = EigVal(IT+1)
       EigVal(2) = Lamda
       GO TO 97
      ENDIF
cc
      ELSE
cc
C  We are maximizing along the lowest NCons modes
C  and possibly mode IT
C
C  Solve iteratively for Lamda
C  Initial guess for Lamda is zero EXCEPT Lamda > EigVal(IT)
C
       Lamda = ZERO
       If(EigVal(IT).GT.ZERO) Lamda = EigVal(IT) + STEP
c
 10    CONTINUE
       NUMIT = NUMIT+1
       TEMP = ZERO
       DO 20 I=1,NCons
       TEMP = TEMP + ( FX(I)*FX(I) )/( Lamda-EigVal(I) )
 20    CONTINUE
       If(IT.GT.NCons)
     $    TEMP = TEMP + ( FX(IT)*FX(IT) )/( Lamda-EigVal(IT) )
       if(iprnt.gt.4) write(IOut,1111) lamda,temp
C
C  Check for Convergence of Lamda
C
       If(Abs(Lamda-TEMP).LT.TOLL) GO TO 30
C
C  Check for Maximum Iterations Exceeded
C
       If(NUMIT.GT.MAXIT) GO TO 96
c
       Lamda = TEMP
       GO TO 10
c
 30    CONTINUE
C
C  At this point we should have an acceptable value for Lamda
C  Make final check
C
       If(IPRNT.GT.1) WRITE(IOut,1400) Lamda
       If(Lamda.LT.ZERO) GO TO 97
       IF(EigVal(IT+1).GT.ZERO.AND.Lamda.GT.EigVal(IT+1)) THEN
C
C  Set up for potential Hessian scaling in subroutine ChkSKAL
C
        EigVal(1) = EigVal(IT+1)
        EigVal(2) = Lamda
        GO TO 97
       ENDIF
cc
      ENDIF
C
C  Calculate the Maximization Step
C
      DO 40 I=1,NCons
      TEMP = FX(I)/( Lamda-EigVal(I) )
      DO 39 J=1,N
      D(J) = D(J) + TEMP*U(J,I)
 39   CONTINUE
 40   CONTINUE
      IF(IT.GT.NCons) THEN
       TEMP = FX(IT)/( Lamda-EigVal(IT) )
       DO 50 J=1,N
       D(J) = D(J) + TEMP*U(J,IT)
 50    CONTINUE
      ENDIF
C
C
      NUMIT = 0
 55   CONTINUE
      JT = 1+IT
      If(JT.GT.NCons+2) JT = NCons+1
C
C  (b)  Minimize along all other Hessian modes
C
      IF(IPRNT.GT.1) THEN
       If(NCons.EQ.0) WRITE(IOut,1500)
       If(NCons.GT.0) WRITE(IOut,1600)
      ENDIF
C
C  Solve Iteratively for Lamda
C  Initial guess for Lamda is zero EXCEPT Lamda < EigVal(JT)
C
      Lamda = ZERO
      IF(EigVal(JT).LT.ZERO) THEN
       Lamda = EIGVAL(JT)-STEP
       Lamda1 = EIGVAL(JT)
       Lamda2 = -BIG
      ENDIF
C
 60   CONTINUE
      NUMIT = NUMIT+1
      TEMP = ZERO
      DO 70 I=NCons+1,NC
      If(I.EQ.IT) GO TO 70
      TEMP = TEMP + ( FX(I)*FX(I) )/( Lamda-EigVal(I) )
 70   CONTINUE
      if(iprnt.gt.4) write(IOut,1111) lamda,temp
C
C  Check for Convergence of Lamda
C
      If(Abs(Lamda-TEMP).LT.TOLL) GO TO 80
C
C  Check for Maximum Iterations Exceeded
C
      If(NUMIT.GT.MAXIT) GO TO 96
C
      IF(EigVal(JT).GT.ZERO) THEN
cc
C  (i)  Simple Iterative Scheme
C
       Lamda = TEMP
       GO TO 60
cc
      ELSE
cc
C  (ii) Cautious Bracketing Scheme
C
       If(TEMP.LT.Lamda) Lamda1 = Lamda
       If(TEMP.GT.Lamda) Lamda2 = Lamda
       If(Lamda2.GT.-BIG) Lamda = HALF*(Lamda1+Lamda2)
       If(Lamda2.EQ.-BIG) Lamda = Lamda-STEP
       GO TO 60
cc
      ENDIF
C
C  At this point we should have an acceptable value for Lamda
C  Make final check
C
 80   CONTINUE
      If(IPRNT.GT.1) WRITE(IOut,1400) Lamda
      If(Lamda.GT.EigVal(JT)) GO TO 97
      If(EigVal(JT).GT.ZERO.AND.Lamda.GT.ZERO) GO TO 97
C
C  Calculate the Minimization Step
C
      DO 90 I=NCons+1,NC
      If(I.EQ.IT) GO TO 90
      TEMP = FX(I)/( Lamda-EigVal(I) )
      DO 89 J=1,N
      D(J) = D(J) + TEMP*U(J,I)
 89   CONTINUE
 90   CONTINUE
C
      IErr = 0
      RETURN
C
C
C  .......................................................
C    ** ERROR SECTION **
C
 96   CONTINUE
      WRITE(IOut,1700)
      RETURN
 97   CONTINUE
      WRITE(IOut,1800)
      RETURN
c
 1000 FORMAT('**WARNING** Mode Switching: Was Following mode ',I3,
     $       ' Now Following mode ',I3)
 1100 FORMAT(' Searching for Lambda that Maximizes along mode ',I3)
 1200 FORMAT(' Mode Following Switched Off')
 1300 FORMAT(' Searching for Lambda that Maximizes Along the',
     $       ' Constraint modes Only')
 1350 FORMAT(' Searching for Lambda that Maximizes Along the',
     $       ' Lowest Non-Constraint mode')
 1400 FORMAT(' Value Taken    Lamda = ',F12.8)
 1500 FORMAT(' Searching for Lambda that Minimizes Along All modes')
 1600 FORMAT(' Searching for Lambda that Minimizes Along All',
     $       ' other modes')
 1700 FORMAT(//,' ********************************************',/,
     $          ' ** UNABLE TO DETERMINE Lambda IN CONFormD **',/,
     $          ' ********************************************',//)
 1800 FORMAT(//,' *********************************************',/,
     $          ' ** ERROR IN DETERMINING Lambda IN CONFormD **',/,
     $          ' *********************************************',//)
c
 1111 format(' in iterative cycle:  Lambda = ',f12.8,' TEMP = ',f12.8)
c
      END
c =====================================================================
c
      SUBROUTINE ConGRAD(NCNTR,  NCons,  method, ICTYP,  RCON,
     $                   IC,     XC,     GC,     GCC )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine modifies the standard Cartesian gradient to
C  incorporate geometric constraints (involving interatomic
C  distances, angles and dihedral angles) using either
C  Lagrange Multipliers or Penalty Functions
C
C  ARGUMENTS
C
C  NCNTR   -  total number of atomic centers (including dummy atoms)
C  NCons   -  number of constraints
C  method  -  method flag     1 - Lagrange multipliers
C                             2 - Penalty functions
C  ICTYP   -  integer array indicating constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane bend    CURRENTLY NOT IMPLEMENTED
C              4 - fixed dihedral angle
C  RCON    -  value of constraint
C  IC      -  list of atoms involved in the constraints
C              IC1-IC2           distance constraint
C              IC1-IC2-IC3       bond angle constraint
C              IC1-IC2-IC3-IC4   all other constraints
C  XC      -  coordinate vector (including Lagrange multipliers)
C  GC      -  gradient vector (including dL/dLambda)
C  GCC     -  matrix of constraint normals
C
C  ** WARNING **  This routine forms the negative of GCC  **
C
C
      DIMENSION ICTYP(NCons),RCON(NCons),XC(*),GC(*),
     $          IC(4,NCons),GCC(3*NCNTR,NCons)
C
      PARAMETER (One=1.0d0,sigma=10.0d0)
C
C
      NCTR3 = 3*NCNTR
c
      CALL ZeroIT(GCC,NCTR3*NCons)
c
      DO 10 IQ=1,NCons
c
      IF(ICTYP(IQ).EQ.1) THEN
cc
C  distance constraint (I-J)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       II = 3*(I-1)
       JJ = 3*(J-1)
c
       XIJ = XC(II+1) - XC(JJ+1)
       YIJ = XC(II+2) - XC(JJ+2)
       ZIJ = XC(II+3) - XC(JJ+3)
       R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
       DCS = One/R
c
       GCC(II+1,IQ) = -DCS*XIJ
       GCC(JJ+1,IQ) =  DCS*XIJ
       GCC(II+2,IQ) = -DCS*YIJ
       GCC(JJ+2,IQ) =  DCS*YIJ
       GCC(II+3,IQ) = -DCS*ZIJ
       GCC(JJ+3,IQ) =  DCS*ZIJ
C
C  now append the dL/dLambda term (Lagrange multipliers)
C  or determine the current constraint multiplier
C
       IT = NCTR3 + IQ
       ROff = RCON(IQ) - R
       If(method.EQ.1) GC(IT) = ROff
       If(method.EQ.2) XC(IT) = sigma*ROff
cc
      ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  angle constraint (I-J-K)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       CALL AngGRAD(NCNTR,I,J,K,XC,Th,.true.,GCC(1,IQ))
C
C  now append the dL/dLambda term (Lagrange multipliers)
C  or determine the current constraint multiplier
C
       IT = NCTR3 + IQ
       ROff = RCON(IQ) - Th
       If(method.EQ.1) GC(IT) = ROff
       If(method.EQ.2) XC(IT) = sigma*ROff
cc
      ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral constraint (I-J-K-L)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       L = IC(4,IQ)
       CALL DihGRAD(NCNTR,I,J,K,L,XC,Dih,.true.,GCC(1,IQ))
       CALL ChkDIH(RCON(IQ),Dih)
C
C  now append the dL/dLambda term (Lagrange multipliers)
C  or determine the current constraint multiplier
C
       IT = NCTR3 + IQ
       ROff = RCON(IQ) - Dih
       If(method.EQ.1) GC(IT) = ROff
       If(method.EQ.2) XC(IT) = sigma*ROff
cc
      ENDIF
c
 10   CONTINUE
C
C  now modify the Cartesian gradient to include
C  the constraint normals
C
      DO 30 IQ=1,NCons
      RLambda = XC(NCTR3+IQ)
      DO 20 I=1,NCTR3
      GC(I) = GC(I) + RLambda*GCC(I,IQ)
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE ConHESS(NCNTR,  NCons,  method, ICTYP,  IC,
     $                   XC,     GCC,    GU,     GD,     HESS )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine modifies the standard Cartesian Hessian to
C  incorporate geometric constraints (involving interatomic
C  distances, angles and dihedral angles) using either
C  Lagrange multipliers or Penalty functions
C
C  ARGUMENTS
C
C  NCNTR   -  total number of atomic centres (including dummy atoms)
C  NCons   -  number of constraints
C  method  -  method flag     1 - Lagrange multipliers
C                             2 - Penalty functions
C  ICTYP   -  integer array indicating constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane bend    CURRENTLY NOT IMPLEMENTED
C              4 - fixed dihedral angle
C  IC      -  list of atoms involved in the constraints
C              IC1-IC2           distance constraint
C              IC1-IC2-IC3       bond angle constraint
C              IC1-IC2-IC3-IC4   all other constraints
C  XC      -  coordinate vector (including Lagrange multipliers)
C  GU      -  scratch space for constraint normal (step up)
C  GD      -  scratch space for constraint normal (step down)
C  HESS    -  Cartesian Hessian (expanded for dummy atoms)
C
C
      DIMENSION ICTYP(NCons),XC(*),GU(3*NCNTR),GD(3*NCNTR),
     $          IC(4,NCons),GCC(3*NCNTR,NCons),HESS(3*NCNTR,3*NCNTR)
      DIMENSION IV(4)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,delta=0.005d0,sigma=10.0d0)
C
C
      NCTR3 = 3*NCNTR
C
C  There are two possible terms modifying the Cartesian Hessian matrix.
C  The first involves the constraint second derivatives and is the same
C  for both the Lagrange and Penalty function methods. The second involves
C  the constraint normals (first derivatives) "squared" (the dot product)
C  and is only included when using Penalty functions.
C
C  1. FIRST TERM
C  -------------
C
      DO 50 IQ=1,NCons
      RLambda = XC(NCTR3+IQ)
c
      IF(ICTYP(IQ).EQ.1) THEN
cc
C  distance constraint (I-J)
C  ** ANALYTICAL **
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       II = 3*(I-1)
       JJ = 3*(J-1)
c
       XIJ = XC(II+1) - XC(JJ+1)
       YIJ = XC(II+2) - XC(JJ+2)
       ZIJ = XC(II+3) - XC(JJ+3)
       R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
c
       RLambda = RLambda/R
c
       XIJD = RLambda*(One-(XIJ/R)**2)
       YIJD = RLambda*(One-(YIJ/R)**2)
       ZIJD = RLambda*(One-(ZIJ/R)**2)
       XYIJ = RLambda*XIJ*YIJ/(R**2)
       XZIJ = RLambda*XIJ*ZIJ/(R**2)
       YZIJ = RLambda*YIJ*ZIJ/(R**2)
C
C  incorporate corrections into Hessian
C
       HESS(II+1,II+1) = HESS(II+1,II+1) - XIJD
       HESS(II+1,JJ+1) = HESS(II+1,JJ+1) + XIJD
       HESS(II+1,II+2) = HESS(II+1,II+2) + XYIJ
       HESS(II+1,JJ+2) = HESS(II+1,JJ+2) - XYIJ
       HESS(II+1,II+3) = HESS(II+1,II+3) + XZIJ
       HESS(II+1,JJ+3) = HESS(II+1,JJ+3) - XZIJ
c
       HESS(JJ+1,II+1) = HESS(II+1,JJ+1)
       HESS(JJ+1,JJ+1) = HESS(JJ+1,JJ+1) - XIJD
       HESS(JJ+1,II+2) = HESS(JJ+1,II+2) - XYIJ
       HESS(JJ+1,JJ+2) = HESS(JJ+1,JJ+2) + XYIJ
       HESS(JJ+1,II+3) = HESS(JJ+1,II+3) - XZIJ
       HESS(JJ+1,JJ+3) = HESS(JJ+1,JJ+3) + XZIJ
c
       HESS(II+2,II+1) = HESS(II+1,II+2)
       HESS(II+2,JJ+1) = HESS(JJ+1,II+2)
       HESS(II+2,II+2) = HESS(II+2,II+2) - YIJD
       HESS(II+2,JJ+2) = HESS(II+2,JJ+2) + YIJD
       HESS(II+2,II+3) = HESS(II+2,II+3) + YZIJ
       HESS(II+2,JJ+3) = HESS(II+2,JJ+3) - YZIJ
c
       HESS(JJ+2,II+1) = HESS(II+1,JJ+2)
       HESS(JJ+2,JJ+1) = HESS(JJ+1,JJ+2)
       HESS(JJ+2,II+2) = HESS(II+2,JJ+2)
       HESS(JJ+2,JJ+2) = HESS(JJ+2,JJ+2) - YIJD
       HESS(JJ+2,II+3) = HESS(JJ+2,II+3) - YZIJ
       HESS(JJ+2,JJ+3) = HESS(JJ+2,JJ+3) + YZIJ
c
       HESS(II+3,II+1) = HESS(II+1,II+3)
       HESS(II+3,JJ+1) = HESS(JJ+1,II+3)
       HESS(II+3,II+2) = HESS(II+2,II+3)
       HESS(II+3,JJ+2) = HESS(JJ+2,II+3)
       HESS(II+3,II+3) = HESS(II+3,II+3) - ZIJD
       HESS(II+3,JJ+3) = HESS(II+3,JJ+3) + ZIJD
c
       HESS(JJ+3,II+1) = HESS(II+1,JJ+3)
       HESS(JJ+3,JJ+1) = HESS(JJ+1,JJ+3)
       HESS(JJ+3,II+2) = HESS(II+2,JJ+3)
       HESS(JJ+3,JJ+2) = HESS(JJ+2,JJ+3)
       HESS(JJ+3,II+3) = HESS(II+3,JJ+3)
       HESS(JJ+3,JJ+3) = HESS(JJ+3,JJ+3) - ZIJD
cc
      ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  angle constraint (I-J-K)
C  ** NUMERICAL **
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       IV(1) = 3*(I-1)
       IV(2) = 3*(J-1)
       IV(3) = 3*(K-1)
       RLambda = RLambda/(delta+delta)
c
       DO 20 IATM=1,3
       II = IV(IATM)
       DO 20 M=1,3
       JJ = II + M
C
C   step up
C
       XC(JJ) = XC(JJ) + delta
       CALL AngGRAD(NCNTR,I,J,K,XC,Ang,.true.,GU)
C
C  step down
C
       XC(JJ) = XC(JJ) - delta - delta
       CALL AngGRAD(NCNTR,I,J,K,XC,Ang,.true.,GD)
C
C  incorporate corrections into Hessian
C
       DO 10 IT=1,NCTR3
       TEMP = RLambda*(GU(IT) - GD(IT))
       HESS(IT,JJ) = HESS(IT,JJ) - TEMP
       HESS(JJ,IT) = HESS(IT,JJ)
 10    CONTINUE
C
C  restore original coordinate
C
       XC(JJ) = XC(JJ) + delta
c
 20    CONTINUE
cc
      ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral constraint
C  ** NUMERICAL **
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       L = IC(4,IQ)
       IV(1) = 3*(I-1)
       IV(2) = 3*(J-1)
       IV(3) = 3*(K-1)
       IV(4) = 3*(L-1)
       RLambda = RLambda/(delta+delta)
c
       DO 40 IAtm=1,4
       II = IV(IAtm)
       DO 40 M=1,3
       JJ = II + M
C
C  step up
C
       XC(JJ) = XC(JJ) + delta
       CALL DihGRAD(NCNTR,I,J,K,L,XC,Ang,.true.,GU)
C
C  step down
C
       XC(JJ) = XC(JJ) - delta - delta
       CALL DihGRAD(NCNTR,I,J,K,L,XC,Ang,.true.,GD)
C
C  incorporate corrections into Hessian
C
       DO 30 IT=1,NCTR3
       TEMP = RLambda*(GU(IT) - GD(IT))
       HESS(IT,JJ) = HESS(IT,JJ) - TEMP
       HESS(JJ,IT) = HESS(IT,JJ)
 30    CONTINUE
C
C  restore original coordinates
C
       XC(JJ) = XC(JJ) + delta
C
 40    CONTINUE
cc
      ENDIF
c
 50   CONTINUE
C
C
      If(method.EQ.1) RETURN
C
C
C  2.  SECOND TERM
C  ---------------
C
      DO 70 I=1,NCTR3
      DO 70 J=1,I
c
      VAL = Zero
      DO 60 IQ=1,NCons
      VAL = VAL + GCC(I,IQ)*GCC(J,IQ)
 60   CONTINUE
c
      HESS(I,J) = HESS(I,J) + sigma*VAL
      HESS(J,I) = HESS(I,J)
 70   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE ConINT(NAtoms, NCons,  XC,     ICTYP,  RCON,
     $                  IC,     ICOMP,  IPComp, PCWght, LCON,
     $                  CTol,   IPRNT,  ICHK,   RCHK,   IAlg)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks all the constraints to see if they are satisfied in
C  the current geometry
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NCons   -  number of constraints
C  XC      -  current cartesian coordinates
C             contains new coordinates on exit
C  ICTYP   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane bend
C              4 - fixed dihedral angle
C              5 - fixed coplanar bend
C              6 - fixed perpendicular bend
C              7 - fixed inverse-power distance
C              9 - composite constraint
C  RCON    -  constraint value (in atomic units)
C  IC      -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C  ICOMP   -  number of primitives in each composite constraint
C  IPComp  -  constraint type and atoms involved in constraint
C              IPComp(1,I) - constraint type (same definition as ICTYP array)
C              IPComp(2,I) to IPComp(5,I) - atoms in constraint
C  PCWght  -  weight of each primitive in composite constraint
C  LCON    -  which constrains are "active" (independent)
C              LCON(I) = 0  constraint is active
C              LCON(I) = 1  constraint should be eliminated
C  CTol    -  accuracy to which constraint must be satisfied
C  IPRNT   -  flag controlling degree of printout
C  ICHK    -  array indicating which constrains are satisfied
C               ICHK(I) = 0 - ith constraint is satisfied
C               ICHK(I) = 1 - ith constraint is not satisfied
C  RCHK    -  values for constraint differences
C             (current value minus desired value)
C  IAlg    -  exit flag    0 - all constraints satisfied
C                         -1 - one or more constraints unsatisfied
C
      REAL*8 XC(3,NAtoms),RCON(NCons),RCHK(NCons)
      INTEGER ICTYP(NCons),IC(4,NCons),LCON(NCons),ICHK(NCons)
      DIMENSION ICOMP(*),IPCOMP(5,*),PCWght(*)
      LOGICAL Satisfy
C
      PARAMETER (One=1.0d0,GMax=0.3d0)
C
C
      IAlg = 0
      If(NCons.EQ.0) RETURN
c
      nc = 0
      np1 = 1
c
      Satisfy = .TRUE.
C
C  Loop over all constraints
C
      DO 10 IQ=1,NCons
c
      ICHK(IQ) = 0
c
      I = IC(1,IQ)
      J = IC(2,IQ)
      K = IC(3,IQ)
      L = IC(4,IQ)
C
C  check that the constraint value is satisfied
C  to within CTol
C
      IF(ICTYP(IQ).EQ.1) THEN
cc
C  distance constraint
C
       XIJ = XC(1,I) - XC(1,J)
       YIJ = XC(2,I) - XC(2,J)
       ZIJ = XC(3,I) - XC(3,J)
       R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
c
       diff = RCON(IQ)-R
cc
      ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  bond angle constraint
C
       CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
c
       diff = RCON(IQ)-Th
cc
      ELSE IF(ICTYP(IQ).EQ.3) THEN
cc
C  out-of-plane bend constraint
C
       CALL OutpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       diff = RCON(IQ)-Th
cc
      ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral angle constraint
C
       CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
       CALL ChkDIH(RCON(IQ),Th)
c
       diff = RCON(IQ)-Th
cc
      ELSE IF(ICTYP(IQ).EQ.5) THEN
cc
C  linear coplanar bend constraint
C
       CALL LincGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       diff = RCON(IQ)-Th
cc
      ELSE IF(ICTYP(IQ).EQ.6) THEN
cc
C  linear perpendicular bend constraint
C
       CALL LinpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       diff = RCON(IQ)-Th
cc
      ELSE IF(ICTYP(IQ).EQ.7) THEN
cc
C  inverse-power distance constraint
C
       XIJ = XC(1,I) - XC(1,J)
       YIJ = XC(2,I) - XC(2,J)
       ZIJ = XC(3,I) - XC(3,J)
       R = One/SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
c
       diff = RCON(IQ)-R
cc
      ELSE
cc
C  composite constraint
C
       nc = nc+1
       np2 = np1 + ICOMP(nc) - 1
       val = 0.0d0
c
       DO 9 IP=np1,np2
       ITYP = IPCOMP(1,IP)
       Wght = PCWGHT(IP)
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
cc
       ELSE IF(ITYP.EQ.2) THEN
cc
C  bond angle constraint
C
        CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
        Th = Th
        val = val+Th*Wght
cc
       ELSE IF(ITYP.EQ.3) THEN
cc
C  out-of-plane bend constraint
C
        CALL OutpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
        Th = Th
        val = val+Th*Wght
cc
       ELSE IF(ITYP.EQ.4) THEN
cc
C  dihedral angle constraint
C
        CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
        Th = Th
        val = val+Th*Wght
cc
       ENDIF
  9    CONTINUE
c
       diff = RCON(IQ)-val
       np1 = np2+1
cc
      ENDIF
C
      RCHK(IQ) = diff
c
      IF(Abs(diff).GT.CTol) THEN      ! constraint not satisfied
       ICHK(IQ) = 1
       Satisfy = .FALSE.
      ENDIF
C
 10   CONTINUE
C
C  now remove dependent constraints
C
      II = 0
      DO 20 I=1,NCons
      IF(LCON(I).EQ.0) THEN
       II = II+1
       ICHK(II) = ICHK(I)
       RCHK(II) = RCHK(I)
       If(Abs(RCHK(II)).GT.GMax) Then
        IOut = ioutfil('iout')
        If(IPRNT.GT.0) WRITE(IOut,1000) I
        RCHK(II) = SIGN(GMax,RCHK(II))
       EndIf
      ENDIF
 20   CONTINUE
c
      If(.NOT.Satisfy) IAlg = -1
C
      RETURN
c
 1000 FORMAT('**WARNING** Constraint Gradient Reduced for',
     $        ' constraint ',I3)
c
      END
c =====================================================================
c
      SUBROUTINE CONNECTM(NAtoms,IAN,XC,thrbnd,IC,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms the connectivity matrix from the interatomic distances
C  and a list of "standard" single bond distances as obtained
C  from the array BndLen.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  cartesian coordinates
C  thrbnd  -  distance ratio above which two atoms are
C             considered to be bonding
C  IC      -  on exit contains atomic connectivity
C  IErr    -  error flag   0 - success
C                         -1 - no default bond length available
C
C
      REAL*8 XC(3,NAtoms)
      INTEGER IAN(NAtoms),IC(NAtoms,NAtoms)
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
C  set error flag
C
      IErr = -1
c
      DO 10 I=2,NAtoms
      DO 10 J=1,I-1
C
C  get interatomic distance I-J
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
       IC(I,J) = 1
       IC(J,I) = 1
      ENDIF
c
 10   CONTINUE
C
C  Check for and remove "excess connectivity"
C
      CALL ChkXSC(NAtoms,XC,IC)
C
C  Now check that the connectivity matrix forms a "closed loop"
C  i.e., that all atoms are connected to the whole
C
      CALL ChkCMPLT(NAtoms,IC,XC,IErr)
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE CONJUGATE(NAtoms, N,      B,      BT,     C,
     $                     D,      X,      thrsh,  MaxIT,  IPRNT,
     $                     PV,     RV,     YV,     ZV,     ST,
     $                     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine solves the linear equation set  (B*B(t))X = C by the
C  method of conjugate gradients
C  Adapted from PP routine <CONGRAD> used in sparse
C  matrix transformation
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  N       -  number of active delocalized internal coordinates
C  B       -  B matrix
C  BT      -  B(t) matrix
C  C       -  right-hand side vector  (B*gcart)
C  D       -  inverse diagonals of A matrix
C  X       -  on entry initial guess for solution
C             on exit  solution (if converged)
C  thrsh   -  convergence criterion
C  MaxIT   -  maximum number of cycles allowed
C  IPRNT   -  print flag
C  PV      -  scratch vector
C  RV      -   ditto
C  YV      -   ditto
C  ZV      -   ditto
C  ST      -   ditto
C  IErr    -  error flag    0 - converged
C                          -1 - failed
C
C
      DIMENSION B(3*NAtoms,N),BT(N,3*NAtoms),C(N),D(N),X(N),
     $          PV(N),RV(N),YV(N),ZV(N),ST(3*NAtoms)
      Dimension cc(2)
C
      PARAMETER (Zero=0.0d0,thrs=1.0d-30)
C
C
      IOut = ioutfil('iout')
C
C  initialize circular counters
C
      NAT3 = 3*NAtoms
      IErr = 0
c
      i0 = 1
      i1 = 2
C
C  set YV = B*(BT*X(0))
C
      DO 10 I=1,NAT3
      ST(I) = SProd(N,BT(1,I),X)
 10   CONTINUE
      DO 11 I=1,N
      YV(I) = SProd(NAT3,B(1,I),ST)
 11   CONTINUE
C
C  set RV = C - B*BT*X(0)
C      ZV = D*RV
C      PV = ZV
C      cc(1) = RV(t)*ZV
C
      Sum = Zero
      DO 20 I=1,N
      RV(I) = C(I) - YV(I)
      ZV(I) = D(I)*RV(I)
      PV(I) = ZV(I)
      Sum = Sum + RV(I)*ZV(I)
 20   CONTINUE
C
C  if sum is very small already means initial gradient
C  is essentially exact
C
      If(Abs(Sum).LT.thrs) RETURN
C
C  start of iterative loop
C
      cc(i0) = Sum
      k = 0
c
 100  CONTINUE
      k = k+1
C
C  set YV = B*(BT*PV(k-1))
C
      DO 30 I=1,NAT3
      ST(I) = SProd(N,BT(1,I),PV)
 30   CONTINUE
      DO 31 I=1,N
      YV(I) = SProd(NAT3,B(1,I),ST)
 31   CONTINUE
c
      Sum = SProd(N,PV,YV)
      bb = cc(i0)/Sum
C
C  X(k)  = X(k-1)  + bb*PV(k-1)
C  RV(k) = RV(k-1) - bb*YV
C
      DO 40 I=1,N
      X(I) = X(I) + bb*PV(I)
      RV(I) = RV(I) - bb*YV(I)
 40   CONTINUE
C
C  check for convergence
C  Absolute value of RV < thrsh
C
      DMax = Zero
      QSQR = Zero
C
      DO 50 I=1,N
      DX = RV(I)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
      QSQR = QSQR + DX**2
 50   CONTINUE
C
      QSQR = SQRT(QSQR/N)
c
      If(IPRNT.GT.2) WRITE(IOut,1200) k,DMax,QSQR
C
C  converged?
C
cc      If(DMax.LT.thrsh) RETURN
      If(DMax.LT.thrsh) Then
       write(IOut,'('' Gradient converged in'',i4,'' cycles'')')k
       return
      endif
c
      If(k.GT.MaxIT) Then
        WRITE(IOut,1300)
        IErr = -1
        RETURN
      EndIf
C
C  prepare for next cycle
C
      DO 60 I=1,N
      ZV(I) = RV(I)*D(I)
 60   CONTINUE
c
      Sum = SProd(N,ZV,RV)
      cc(i1) = Sum
      bb = cc(i1)/cc(i0)
C
C  PV(k) = ZV + bb*PV(k-1)
C
      DO 70 I=1,N
      PV(I) = ZV(I) + bb*PV(I)
 70   CONTINUE
c
      ii = i1
      i1 = i0
      i0 = ii
      GO TO 100
c
 1200 FORMAT(5X,'Cycle: ',I3,'  Maximum deviation: ',F12.8,
     $          '  RMS deviation: ',F12.8)
 1300 FORMAT(/,2X,'***ERROR*** Exceeded allowed number of iterative',
     $            ' cycles in <CONJUGATE>')
c
      END
c =====================================================================
c
      SUBROUTINE ConORD(NAT3,   NCon,   NDim,   NEG,    EigVAL,
     $                  U,      IPRNT,  V,      VM,     INDX,
     $                  ISYM )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ** WARNING **  THIS CODE NOT PARTICULARLY WELL CHECKED  **
C  When there are more Lagrangian modes with negative eigenvalues
C  than there are constraint modes, this routine partitions the
C  negative modes into their respective subspaces
C
C  ARGUMENTS
C
C  NAT3    -  3*number of atoms
C  NCon    -  number of constraints
C  NDim    -  total dimension of each mode
C  NEG     -  number of negative eigenvalues
C  EigVal  -  Lagrangian eigenvalues
C  U       -  corresponding eigenvectors in columns
C  IPRNT   -  controls level of print out
C  V       -  scratch space for eigenvalues
C  VM      -  scratch space for eigenvectors
C  INDX    -  scratch space for indexing array
C  ISYM    -  scratch space for tracking symmetry equivalent
C             constraint components (if any)
C
      REAL*8 EigVal(*),U(NDim,*),V(*),VM(NDim,*)
      INTEGER INDX(*),ISYM(*)
      LOGICAL remove
C
      PARAMETER (Zero=0.0d0,TolW=1.0d-6)
C
C
      IOut = ioutfil('iout')
      NCons = NDim-NAT3           ! number of original constraints
      remove = NCons.NE.NCon      ! any constraints removed by symmetry?
      CALL IZeroIT(ISYM,NCons)
      IT = 0
C
C  Check each of the NEG negative modes in turn and
C  select for the constraint space the mode with the
C  largest weighting for each constraint
C
      DO 40 IC=1,NCons
      If(ISYM(IC).NE.0) GO TO 40
c
      WMax = Zero
      IT = IT+1
c
      DO 10 I=1,NEG-IT+1
      Weight = U(NAT3+IC,I)**2
      IF(Weight.GT.WMax) THEN
       IMax = I
       WMax = Weight
      ENDIF
 10   CONTINUE
C
C  IMax is the selected mode
C  put it in the constraint space and shift
C  all other negative modes down
C
      V(IT) = EigVal(IMax)
      CALL CpyVEC(NDim,U(1,IMax),VM(1,IT))
c
      DO 20 J=IMax+1,NEG-IT+1
      EigVal(J-1) = EigVal(J)
      CALL CpyVEC(NDim,U(1,J),U(1,J-1))
 20   CONTINUE
C
      IF(remove) THEN
C
C  check the constraint weightings of the selected mode and
C  set ISYM so as to skip over components with identical
C  weightings since these are symmetry related
C
       DO 30 J=1,NCons
       Weight = VM(NAT3+J,IT)**2
       If(Abs(WMax-Weight).LT.TolW) ISYM(J) = 1
 30    CONTINUE
      ENDIF
c
 40   CONTINUE
C
C  check
C  we should have selected NCon modes
C
      IF(IT.NE.NCon) THEN
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
C  At this point we should have NCon constraints in VM
C  with their corresponding eigenvalues in V.
C  The lowest NEG-NCon modes in U will be the negative
C  non-constraint modes - move these up to join the
C  positive modes
C
      DO 50 I=1,NEG-NCon
      EigVal(NCon+I) = EigVal(I)
      CALL CpyVEC(NDim,U(1,I),U(1,NCon+I))
 50   CONTINUE
C
C  now put back the constraint modes after sorting
C  them into ascending energy ordering
C
      CALL INDEXA(NCon,V,INDX)
C
      DO 60 I=1,NCon
      EigVal(I) = V(INDX(I))
      CALL CpyVEC(NDim,VM(1,INDX(I)),U(1,I))
 60   CONTINUE
C
      IF(IPRNT.GT.3) THEN
       WRITE(IOut,1100)
       WRITE(IOut,1200)
       WRITE(IOut,1300) (EigVal(I),I=1,NEG)
       WRITE(IOut,1400)
       CALL PrntMAT(NEG,NDim,NDim,U)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Something Wrong in Assigning Constraint',
     $            ' modes')
 1100 FORMAT(/,'  ReOrdering Negative modes')
 1200 FORMAT('  Hessian Eigenvalues:')
 1300 FORMAT(1X,6F12.6)
 1400 FORMAT(/,'  Hessian Eigenvectors:')
c
      END
c =====================================================================
c
      SUBROUTINE ConVEC(NDEG,   NPrim,  NCons,  ICNUM,  ICOMP,
     $                  PCWght, IPCNUM, UT,     thrsh,  IPRNT,
     $                  MCON,   NCon,   VC,     VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Set up and symmetrize initial constraint vectors and form
C  projection onto the current space of active delocalized internals
C
C  ARGUMENTS
C
C  NDEG    -  number of active internals
C             (number of degrees of freedom)
C  NPrim   -  number of primitive internals
C  NCons   -  number of constrained (fixed) primitives
C  ICNUM   -  array indicating which primitives correspond
C             to the desired constraints
C  ICOMP   -  array indicating number of primitives in each
C             composite constraint
C  PCWght  -  weight of each primitive in composite constraint
C  IPCNUM  -  array indicating which primitives correspond
C             to the desired composite constraints
C  UT      -  natural internal coordinate coefficients
C  thrsh   -  zero threshold
C  IPRNT   -  print flag
C  MCON    -  integer sort array
C
C  on exit
C
C  NCon    -  number of symmetrized constraints
C  VC      -  symmetrized constraint vectors
C  VM      -  projected constraint vectors
C
C
      DIMENSION UT(NPrim,NDEG),VC(NPrim,NCons),VM(NPrim,NCons)
      DIMENSION ICNUM(NCons),MCON(NCons)
      DIMENSION ICOMP(*),PCWght(*),IPCNUM(*)
c
      PARAMETER (One=1.0d0)
C
C
      IOut = ioutfil('iout')
c
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
C  start by setting up unit vectors
C
      CALL ZeroIT(VC,NPrim*NCons)
c
      nc = 0
      JT = 0
      KT = 0
c
      DO 10 I=1,NCons
      If(ICNUM(I).NE.0) Then
       II = ICNUM(I)
       VC(II,I) = One
      Else
c -- composite constraint
       nc = nc+1
       np = ICOMP(nc)
c -- calculate normalization factor
       vnorm = 0.0d0
       DO 8 J=1,np
       KT = KT+1
       vnorm = vnorm + PCWght(KT)**2
  8    CONTINUE
       vnorm = One/DSQRT(vnorm)
c --  set up composite constraint vector
       DO 9 J=1,np
       JT = JT+1
       II = IPCNUM(JT)
       VC(II,I) = vnorm*PCWght(JT)
  9    CONTINUE
      EndIf
 10   CONTINUE
c
      CALL ZeroIT(VM,NPrim*NCons)
C
C  project the constraints onto the active space
C
      DO 30 I=1,NCons
      DO 30 J=1,NDEG
      proj = SProd(NPrim,VC(1,I),UT(1,J))
      DO 20 K=1,NPrim
      VM(K,I) = VM(K,I) + proj*UT(K,J)
 20   CONTINUE
 30   CONTINUE
c
      If(IPRNT.GT.5) Then
       WRITE(IOut,1100)
       CALL PrntMAT(NCons,NPrim,NCons,VM)
      EndIf
C
C  preliminary symmetrization and elimination of constraints
C
C   first vector always kept
C
      CALL IZeroIT(MCON,NCons)
      MCON(1) = 1
      IC = 1
c
 100  CONTINUE
C
C  find vectors whose components are identical in magnitude
C  to current constraint vector
C
      DO 50 I=IC+1,NCons
      IF(MCON(I).EQ.0) THEN
       DO 40 J=1,NPrim
       If(Abs(VM(J,IC))-Abs(VM(J,I)).GT.thrsh) GO TO 50
 40    CONTINUE
C
C  projected constraint vectors are the same within a sign
C  determine sign
C
       DO 45 J=1,NPrim
       If(Abs(VM(J,IC)).GT.thrsh) Then
        If(Abs(VM(J,IC)-VM(J,I)).GT.thrsh) Then
         MCON(I) = -IC
        Else
         MCON(I) = IC
        EndIf
        GO TO 50
       EndIf
 45    CONTINUE
      ENDIF
c
 50   CONTINUE
cc      write(6,*) ' IC:',ic,' MCON array is:'
cc      do i=1,ncons
cc      write(6,*) i,'  ',mcon(i)
cc      enddo
C
C  collect all vectors that are the same
C
      II = 0
      DO 70 I=1,NCons
      If(Abs(MCON(I)).EQ.IC) Then
       II = II+1
       If(II.GT.1) Then
        ISign = SIGN(1,MCON(I))
        DO 60 J=1,NPrim
        VC(J,IC) = VC(J,IC) + ISign*VC(J,I)
        VM(J,IC) = VM(J,IC) + ISign*VM(J,I)
 60     CONTINUE
       EndIf
      EndIf
 70   CONTINUE
C
C  prepare next vector
C
      IC = IC+1
      DO 80 I=IC,NCons
      If(MCON(I).EQ.0) Then
       MCON(I) = IC
       CALL CpyVEC(NPrim,VM(1,I),VM(1,IC))
       CALL CpyVEC(NPrim,VC(1,I),VC(1,IC))
       GO TO 100
      EndIf
 80   CONTINUE
C
C  if we get here there we are all done
C
      NCon = IC-1
C
C  Normalize the original constraint vectors
C
      DO 90 I=1,NCon
      snorm = SProd(NPrim,VC(1,I),VC(1,I))
      snorm = One/SQRT(snorm)
      CALL VScal(NPrim,snorm,VC(1,I))
 90   CONTINUE
C
      RETURN
c
 1000 FORMAT(/' Constructing Initial Set of Constraint Vectors')
 1100 FORMAT(' Projected Constraint Vectors')
c
      END
c =====================================================================
c
      SUBROUTINE CpyMAT(N,N1,N2,A,B)
      REAL*8 A(N1,N1),B(N2,N2)
C
C  Copies the upper left NxN block of matrix A into matrix B
C
C  ARGUMENTS
C
C  N   -  size of block to copy
C  N1  -  dimension of matrix A
C  N2  -  dimension of matrix B
C  A   -  matrix being copied from
C  B   -  matrix being copied to
C
C
      IF(N.GT.N1.OR.N.GT.N2) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
      DO 10 J=1,N
      DO 10 I=1,N
      B(I,J) = A(I,J)
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Matrix block too big to copy in CpyMAT')
c
      END
