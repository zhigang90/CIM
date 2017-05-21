c ==================================================================
c  GEOMETRY OPTIMIZATION ROUTINES D-G          JB   October 1999
c ==================================================================
c
      SUBROUTINE DefHES(NDEG,   intcor, NCon,   ktyp,   klist,
     $                  XC,     IAN,    FC,     UT,     HINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate Hessian by forming UT(t)*H(prim)*UT
C  where H(prim) is a "Hessian" diagonal matrix in the
C  primitive space with appropriate force constants
C  given for each primitive type
C
C  ARGUMENTS
C
C  NDEG    -  number of degrees of freedom
C             (and number of natural internal coordinates)
C  intcor  -  number of primitive internals
C  NCon    -  number of constraints
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C  XC      -  Cartesian coordinates
C  IAN     -  list of atomic numbers
C  FC      -  space for primitive force constants
C  UT      -  set of active delocalized internal coordinates
C  HINT    -  on exit contains scaled-guess Hessian/Inverse Hessian
C
C
      REAL*8 XC(3,*),FC(intcor),UT(intcor,NDEG),HINT(NDEG,NDEG)
      DIMENSION ktyp(intcor),klist(4,intcor),IAN(*)
C
      PARAMETER (Zero=0.0d0)
C
C
C  form primitive force constant array
C
      DO 10 I=1,intcor
      FC(I) = ForceCnst(ktyp(I),klist(1,I),klist(2,I),klist(3,I),
     $                  klist(4,I),NCon,IAN,XC)
 10   CONTINUE
C
C  now form the Hessian
C  this is simply UT(t)*UT scaled by force constants
C
      DO 30 I=1,NDEG
      DO 30 J=1,I
      Val = Zero
      DO 20 K=1,intcor
      Val = Val + UT(K,I)*UT(K,J)*FC(K)
 20   CONTINUE
      HINT(I,J) = Val
      HINT(J,I) = Val
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE DefHESZ(NZ,NVar,IG,HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Set up a default diagonal Hessian suitable for a minimization
C  in Z-Matrix internal coordinates
C
C  Diagonal elements are assigned depending on the coordinate
C  type according to:
C
C  stretch                 0.5
C  all other               0.2
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  NVar    -  number of variables in Z-Matrix
C             i.e. parameters to be optimized
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  HESS    -  Hessian matrix
C
C
      DIMENSION IG(*),HESS(NVar,NVar)
C
C
      CALL ZeroIT(HESS,NVar*NVar)
      IT = 0
c
      DO 10 I=1,NZ-1
      IF(IG(I).EQ.0) THEN
       IT = IT+1
       HESS(IT,IT) = 0.5d0
      ENDIF
 10   CONTINUE
c
      DO 20 I=NZ,3*NZ-6
      IF(IG(I).EQ.0) THEN
       IT = IT+1
       HESS(IT,IT) = 0.2d0
      ENDIF
 20   CONTINUE
C
C  Check that correct number of elements set
C
      IF(IT.NE.NVar) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Wrong Number of Variables in <DefHESZ>')
c
      END
c =====================================================================
c
      SUBROUTINE DefHPRIM(NPrim,  NCon,   ktyp,   klist,  XC,
     $                    IAN,    HPRIM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Set up default Hessian in full primitive space
C  with appropriate force constants for each primitive type
C
C  ARGUMENTS
C
C  NPrim   -  number of primitive internals
C  NCon    -  number of constraints
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C                7 - 1/R**6
C  klist   -  list of atoms involved in each primitive internal
C  XC      -  Cartesian coordinates
C  IAN     -  list of atomic numbers
C  HPRIM   -  on exit contains primitive Hessian
C
C
      REAL*8 XC(3,*),HPRIM(NPrim,NPrim)
      DIMENSION ktyp(NPrim),klist(4,NPrim),IAN(*)
C
C
C  now form the Hessian
C
      CALL ZeroIT(HPRIM,NPrim*NPrim)
      DO 10 I=1,NPrim
      HPRIM(I,I) = ForceCnst(ktyp(I),klist(1,I),klist(2,I),klist(3,I),
     $                       klist(4,I),NCon,IAN,XC)
 10   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE DTOPOLOGY(NAtoms, NMol,   IMOL,   XC,     NIC,
     $                     NVib,   Cutoff, IPRNT,  IC,     intcor,
     $                     ktyp,   klist)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates "bond" stretches between separate molecules for
C  use as a set of inverse power distance-only coordinates
C  Used to connect different molecules in cluster optimizations
C  ** TEMPORARY BRUTE FORCE ALGORITHM **
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  NMol    -  number of molecules (for, e.g., cluster optimizations)
C  IMOL    -  pointers to start/end of molecules in XC array
C  XC      -  Cartesian coordinates
C  NIC     -  maximum number of primitive internal coordinates
C             that can be generated
C             ** NOTE:  This is typically far more than 3*NATOMS-6
C                       as many redundancies are found that are
C                       dealt with later
C  NVib    -  total number of degrees of freedom
C             (at least this many coordinates must be generated in total)
C  CutOff  -  distance cutoff for bonding (in Angstroms)
C             atoms less than a distance cutoff from one another
C             are considered as bonded
C  IPRNT   -  print flag
C  IC      -  scratch storage for connectivity matrix
C  INDX    -  scratch storage for index sorting
C  Dist    -  scratch storage for interatomic distances
C
C  on exit
C
C  intcor  -  number of primitive internals found
C  ktyp    -  integer array containing internal coordinate type
C             ** NOTE:  These types are:
C                       1 - stretch
C                       2 - bend
C                       3 - out-of-plane bend
C                       4 - torsion
C                       5 - linear coplanar bend
C                       6 - linear perpendicular bend
C                       7 - inverse power stretch
C             This routine only assigns inverse power stretches
C  klist   -  list of atoms involved in each primitive
C
C
      DIMENSION XC(3,NAtoms),IMOL(NMol+1)
      DIMENSION IC(NAtoms,NAtoms),ktyp(NIC),klist(4,NIC)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      PARAMETER (One=1.0d0)
C
C
      IOut = ioutfil('iout')
      CutOff0 = CutOff                 ! save initial cutoff
      intcor0 = intcor                 !  ditto
c
 500  CONTINUE
      CALL IZeroIT(IC,NAtoms*NAtoms)
c
      If(IPRNT.GT.1) Then
       WRITE(IOut,1000)
       WRITE(IOut,1100) CutOff
      EndIf
C
C  convert cutoff to distance squared in au
C
      CutOf2 = (CutOff*ANTOAU)**2
C
C  Loop over all molecules
C  generate primitives for intermolecular bonding
C
      DO 200 Mol1=2,NMol
      IStrt = IMOL(Mol1)+1
      IEnd  = IMOL(Mol1+1)
      DO 100 Mol2=1,Mol1-1
      JStrt = IMOL(Mol2)+1
      JEnd  = IMOL(Mol2+1)
C
C  Loop over atoms in the two molecules
C
      DO 20 I=IStrt,IEnd
      DO 10 J=JStrt,JEnd
C
C  get interatomic distance I-J
C
      X = XC(1,I) - XC(1,J)
      Y = XC(2,I) - XC(2,J)
      Z = XC(3,I) - XC(3,J)
      R2 = X*X + Y*Y + Z*Z
C
C  atoms "bonded" if distance less than cutoff
C
      IF(R2.LT.CutOf2) THEN
       IC(I,J) = 1
       IC(J,I) = 1
      ENDIF
 10   CONTINUE
 20   CONTINUE
C
C  End of loop over molecules
C
 100  CONTINUE
 200  CONTINUE
C
C  assign inverse-power stretches
C
      DO 40 I=2,NAtoms
      DO 30 J=1,I-1
      IF(IC(I,J).EQ.1) THEN
       intcor = intcor+1
       ktyp(intcor) = 7
       klist(1,intcor) = I
       klist(2,intcor) = J
      ENDIF
 30   CONTINUE
 40   CONTINUE
C
C  check
C  do we have enough coordinates?
C
      IF(intcor.LT.NVib) THEN
       CutOff = CutOff + One
       intcor = intcor0
       If(IPRNT.GT.1) WRITE(IOut,2000)
       if(CutOff.gt.40.0d0) call nerror(1,'dtopology',
     $              'cutoff is too large',0,0)
       GO TO 500
      ENDIF
c
      If(IPRNT.GT.4) Then
       WRITE(IOut,1200)
       Do I=intcor0+1,intcor
       WRITE(IOut,1300) klist(1,I),klist(2,I)
       EndDo
      EndIf
c
      If(IPRNT.GT.2) WRITE(IOut,1400) intcor-intcor0
c
      CutOff = CutOff0             ! restore original cutoff
      RETURN
c
 1000 FORMAT(' Generating Primitive Inverse-Stretch Coordinates')
 1100 FORMAT(' Cutoff for bonding is ',F10.6,' Angstroms')
 1200 FORMAT(' Primitive Inverse Stretches:')
 1300 FORMAT(5X,4I5)
 1400 FORMAT(' There are ',I6,' Inverse Stretches')
 2000 FORMAT('**WARNING**  Insufficient Primitives - Increasing Cutoff')
c
      END
c =====================================================================
c
      SUBROUTINE EstNPrim(NAtoms, XC,     CutOff, NPrim)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Estimates maximum number of "bonds" between all atoms in the
C  system based on number of interatomic distances less than CutOff.
C  This should be an upper bound to the total number of primitives
C  in a cluster optimization using inverse power distance coordinates.
C  ** TEMPORARY BRUTE FORCE ALGORITHM **
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  XC      -  Cartesian coordinates
C  CutOff  -  distance cutoff for bonding (in Angstroms)
C             atoms less than a distance cutoff from one another
C             are considered as bonded
C  NPrim   -  on exit number of "n\bonds" found
C
C
      DIMENSION XC(3,NAtoms)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
C
      NPrim = 0
C
C  convert cutoff to distance squared in au
C
      CutOf2 = (CutOff*ANTOAU)**2
C
C  Loop over all atoms
C  generate primitives for intermolecular bonding
C
      DO 20 I=2,NAtoms
      DO 10 J=1,I-1
C
C  get interatomic distance I-J
C
      X = XC(1,I) - XC(1,J)
      Y = XC(2,I) - XC(2,J)
      Z = XC(3,I) - XC(3,J)
      R2 = X*X + Y*Y + Z*Z
C
C  atoms "bonded" if distance less than cutoff
C
      If(R2.LT.CutOf2) NPrim=NPrim+1
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE ExpandV(NS,NCons,NCon,ICHK,V,VX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine expands a vector from the current to the
C  "full" constrained optimization space, including padding
C  for satisfied constraints
C
C  ARGUMENTS
C
C  NS      -  number of unconstrained degrees of freedom
C  NCons   -  total number of constraints
C  NCon    -  number of unsatisfied constraints
C  ICHK    -  array indicating which constrains are satisfied
C               ICHK(I) = 0 - ith constraint is satisfied
C               ICHK(I) = 1 - ith constraint is not satisfied
C  V       -  vector to be expanded
C  VX      -  expanded vector
C
C
      DIMENSION ICHK(NCons),V(NS+2*NCon),VX(NS+2*NCons)
C
      PARAMETER (Zero=0.0d0)
C
C
C  The first NS (the active) coordinates are always kept
C  Pad out the NCon unsatisfied constraints to the full NCons
C
      IT = 0
      DO 10 I=1,NCons
      IF(ICHK(I).EQ.1) THEN
       IT = IT+1
       VX(NS+I) = V(NS+IT)
       VX(NS+NCons+I) = V(NS+NCon+IT)
      ELSE
       VX(NS+I) = Zero
       VX(NS+NCons+I) = Zero
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE FillZ(NZ,XPrim,INDX,GEO)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Fills Z-Matrix with selected primitives according to
C  indexing array, INDX
C
C  ARGUMENTS
C
C  NZ      -  number of Z-matrix centres
C  XPrim   -  values of primitives
C  INDX    -  indexing array into primitives
C  GEO     -  on exit contains Z-matrix values
C
C
      DIMENSION XPrim(*),INDX(3*NZ),GEO(NZ,3)
C
C
      IT = 0
C
C  fill bond lengths
C
      DO 10 I=2,NZ
      IT = IT+1
      GEO(I,1) = XPrim(INDX(IT))
 10   CONTINUE
C
C  fill bond angles
C
      DO 20 I=3,NZ
      IT = IT+1
      GEO(I,2) = XPrim(INDX(IT))
 20   CONTINUE
C
C  fill torsions
C
      DO 30 I=4,NZ
      IT = IT+1
      GEO(I,3) = XPrim(INDX(IT))
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE FixCON(NAtoms, XC,     NFix,   IFIX,   NPrim,
     $                  ktyp,   klist,  IPRNT,  NCons,  ICTYP,
     $                  ICON,   RCON)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sets up constraints for fixed atoms
C  Implemented by fixing ALL primitives that comprise
C  ONLY fixed atoms
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NFix    -  on entry simply indicates there are fixed surface atoms
C             on exit number of fixed atoms
C  IFIX    -  list of fixed surface atoms (considered as X,Y,Z
C              components for compatibility with fixed Cartesians)
C              0 - coordinate active
C              1 - coordinate inactive (will be "fixed")
C  NPrim   -  number of primitive internals
C  ktyp    -  integer array containing internal coordinate type
C             ** NOTE:  These types are:
C                       1 - stretch
C                       2 - bend
C                       3 - out-of-plane bend
C                       4 - torsion
C                       5 - linear coplanar bend
C                       6 - linear perpendicular bend
C                       7 - inverse-power distance
C  klist   -  list of atoms involved in each primitive
C  IPRNT   -  print flag
C  NCons   -  number of existing constraints
C             (on exit includes additional surface constraints)
C  ICTYP   -  constraint type (same as in ktyp array)
C  ICON    -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C  RCON    -  constraint values
C
C
      DIMENSION XC(3,NAtoms),IFIX(3,NAtoms),ktyp(NPrim),klist(4,NPrim),
     $          ICTYP(NCons),ICON(4,NCons),RCON(NCons)
C
C
      IOut = ioutfil('iout')
C
C  How many fixed surface atoms are there?
C
      NFix = 0
      DO 5 I=1,NAtoms
      NFix = NFix + IFIX(1,I)
 5    CONTINUE
      If(NFix.EQ.0) RETURN
c
      If(IPRNT.GT.2) WRITE(IOut,1000) NFix
C
C  Constrain ALL primitives that comprise totally fixed atoms
C
      NCon0 = NCons
c
      DO 10 I=1,NPrim
      II = klist(1,I)
      JJ = klist(2,I)
      KK = klist(3,I)
      LL = klist(4,I)
c
      IF(ktyp(I).EQ.1.OR.ktyp(I).EQ.7) THEN
cc
       If(IFIX(1,II).EQ.1.AND.IFIX(1,JJ).EQ.1) Then
        NCons = NCons+1
        ICTYP(NCons) = ktyp(I)
        ICON(1,NCons) = II
        ICON(2,NCons) = JJ
        RCON(NCons) = SQRT( (XC(1,II)-XC(1,JJ))**2 +
     $                      (XC(2,II)-XC(2,JJ))**2 +
     $                      (XC(3,II)-XC(3,JJ))**2 )
       EndIf
cc
      ELSE IF(ktyp(I).EQ.2) THEN
cc
       If(IFIX(1,II).EQ.1.AND.IFIX(1,JJ).EQ.1.AND.IFIX(1,KK).EQ.1) Then
        NCons = NCons+1
        ICTYP(NCons) = ktyp(I)
        ICON(1,NCons) = II
        ICON(2,NCons) = KK
        ICON(3,NCons) = JJ
        Call AngGRAD(MAtoms,II,KK,JJ,XC,Th,.false.,jnk)
        RCON(NCons) = Th
       EndIF
cc
      ELSE IF(ktyp(I).EQ.3.OR.ktyp(I).EQ.4) THEN
cc
       If(IFIX(1,II).EQ.1.AND.IFIX(1,JJ).EQ.1.AND.
     $    IFIX(1,KK).EQ.1.AND.IFIX(1,LL).EQ.1) Then
        NCons = NCons+1
        ICTYP(NCons) = ktyp(I)
        ICON(1,NCons) = II
        ICON(2,NCons) = JJ
        ICON(3,NCons) = KK
        ICON(4,NCons) = LL
        If(ktyp(I).EQ.3) Then
         Call OutpGRAD(MAtoms,II,JJ,KK,LL,XC,Th,.false.,jnk)
        Else
         Call DihGRAD(MAtoms,II,JJ,KK,LL,XC,Th,.false.,jnk)
        EndIf
        RCON(NCons) = Th
       EndIf
cc
      ELSE IF(ktyp(I).EQ.5.OR.ktyp(I).EQ.6) THEN
cc
       If(IFIX(1,II).EQ.1.AND.IFIX(1,JJ).EQ.1.AND.
     $    IFIX(1,KK).EQ.1.AND.IFIX(1,LL).EQ.1) Then
        NCons = NCons+1
        ICTYP(NCons) = ktyp(I)
        ICON(1,NCons) = II
        ICON(2,NCons) = KK
        ICON(3,NCons) = JJ
        ICON(4,NCons) = LL
        If(ktyp(I).EQ.5) Then
         Call LincGRAD(MAtoms,II,KK,JJ,LL,XC,Th,.false.,jnk)
        Else
         Call LinpGRAD(MAtoms,II,KK,JJ,LL,XC,Th,.false.,jnk)
        EndIf
        RCON(NCons) = Th
       EndIf
cc
      ENDIF
 10   CONTINUE
c
      If(IPRNT.GT.2) WRITE(IOut,1100) NCons-NCon0
C
      RETURN
c
 1000 FORMAT(' There are ',I4,' Fixed Atoms')
 1100 FORMAT(' Added ',I4,' Fixed Atom Constraints')
c
      END
c =====================================================================
c
      SUBROUTINE FndNEG(N,EigVal,NEG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Find the number of negative Hessian eigenvalues
C  (it is assumed that the eigenvalues are ordered)
C
C  ARGUMENTS
C
C  N       -  number of eigenvalues to search
C  EigVal  -  Hessian eigenvalues
C  NEG     -  on return contains the number of negative eigenvalues
C
C
      REAL*8 EigVal(N)
      PARAMETER (Zero=0.0d0)
C
      DO 10 I=1,N
      If(EigVal(I).GT.Zero) GO TO 20
 10   CONTINUE
 20   NEG = I-1
C
      RETURN
      END
c =====================================================================
c
      DOUBLE PRECISION FUNCTION ForceCnst(ITyp,I,J,K,L,NCon,IAN,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Estimates the diagonal force constant for a particular
C  primitive internal coordinate
C
C  ARGUMENTS
C
C  ITyp    -  primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C                7 - inverse-power distance
C  I,J     -  atoms involved in each primitive internal
C  K,L         (maximum of 4)
C  NCon    -  number of constraints
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C              (for determining distance)
C
C  Based on:
C  "Estimating the Hessian for gradient-type geometry optimizations"
C   H.B.Schlegel,  Theoret.Chim.Acta  66 (1984) 333
C
C
      DIMENSION IAN(*),XC(3,*)
c
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
      IF(ITyp.EQ.1) THEN
C
C  stretch
C
       R = SQRT( (XC(1,I) - XC(1,J))**2 +
     $           (XC(2,I) - XC(2,J))**2 +
     $           (XC(3,I) - XC(3,J))**2 )
c
       IAT = IAN(I)
       JAT = IAN(J)
c
       If(JAT.GT.IAT) Then
        ITmp = IAT
        IAT = JAT
        JAT = ITmp
       EndIf
c
       A = 1.734d0
       If(IAT.EQ.1) Then
        B = -0.244d0
       Else If(IAT.LT.11) Then
        B = 1.085d0
        If(JAT.EQ.1) B = 0.352d0
       Else
        B = 2.068d0
        If(JAT.LT.11) B = 1.522d0
        If(JAT.EQ.1) B = 0.660d0
       EndIf
c
       ForceCnst = A/((R-B)**3)
cc
      ELSE IF(ITyp.EQ.2) THEN
C
C  bend
C
       IAT = IAN(I)
       JAT = IAN(J)
       KAT = IAN(K)
c
       If(IAT.GT.1.AND.JAT.GT.1.AND.KAT.GT.1) THEN
        ForceCnst = 0.250d0
       Else
        ForceCnst = 0.160d0
       EndIf
cc
      ELSE IF(ITyp.EQ.3) THEN
C
C  out-of-plane-bend
C
       ForceCnst = 0.045d0
cc
      ELSE IF(ITyp.EQ.4) THEN
C
C  torsion
C
       IF(NCon.EQ.0) THEN
        R = SQRT( (XC(1,J) - XC(1,K))**2 +
     $            (XC(2,J) - XC(2,K))**2 +
     $            (XC(3,J) - XC(3,K))**2 )
c
        bondL =  SQRT(BndLen(IAN(J))*BndLen(IAN(K)))
c
        FTor = 0.0025d0 - 0.07d0*(R-bondL)
        ForceCnst = MAX(FTor,0.0020d0)
       ELSE
        ForceCnst = 0.045d0    ! no small force constants with constraints
       ENDIF
cc
      ELSE IF(ITyp.EQ.7) THEN
C
C  inverse-power distance
C
       ForceCnst = 0.2d0
cc
      ELSE
C
C  linear bends
C
       ForceCnst = 0.2d0
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
c
      DOUBLE PRECISION FUNCTION ForceCnst1(ITyp,I,J,K,L,NCon,IAN,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Estimates the diagonal force constant for a particular
C  primitive internal coordinate
C
C  ARGUMENTS
C
C  ITyp    -  primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C                7 - inverse-power distance
C  I,J     -  atoms involved in each primitive internal
C  K,L         (maximum of 4)
C  NCon    -  number of constraints
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C              (for determining distance)
C
C  Based on:
C  "Estimating the Hessian for gradient-type geometry optimizations"
C   H.B.Schlegel,  Theoret.Chim.Acta  66 (1984) 333
C ** APPEARS TO BE MUCH WORSE THAN THE ORIGINAL ALGORITHM **
C
C
      DIMENSION IAN(*),XC(3,*)
c
      REAL*8 VdW(103)      !atomic VdW radii in au
C
c                H       He      Li      Be      B       C       N
      DATA VdW/2.886d0,2.362d0,2.451d0,2.745d0,4.083d0,3.851d0,3.660d0,
c                O       F       Ne      Na      Mg      Al      Si
     $         3.500d0,3.364d0,3.243d0,2.983d0,3.021d0,4.499d0,4.295d0,
c                P       S       Cl      Ar      K       Ca      Sc
     $         4.147d0,4.035d0,3.947d0,3.868d0,3.812d0,3.399d0,3.295d0,
c                Ti      V       Cr      Mn      Fe      Co      Ni
     $         3.175d0,3.144d0,3.023d0,2.961d0,2.912d0,2.872d0,2.834d0,
c                Cu      Zn      Ga      Ge      As      Se      Br
     $         3.495d0,2.763d0,4.383d0,4.280d0,4.230d0,4.205d0,4.189d0,
c                Kr      Rb      Sr      Y       Zr      Nb      Mo
     $         4.141d0,4.114d0,3.641d0,3.345d0,3.124d0,3.165d0,3.052d0,
c                Tc      Ru      Rh      Pd      Ag      Cd      In   
     $         2.998d0,2.963d0,2.929d0,2.899d0,3.148d0,2.848d0,4.463d0,
c                Sn      Sb      Te      I_      Xe      Cs      Ba   
     $         4.392d0,4.420d0,4.470d0,4.50d0, 4.404d0,4.517d0,3.703d0,
c                La      Ce      Pr      Nd      Pm      Sm      Eu   
     $         3.522d0,3.556d0,3.606d0,3.575d0,3.547d0,3.520d0,3.493d0,
c                Gd      Tb      Dy      Ho      Er      Tm      Yb   
     $         3.368d0,3.451d0,3.428d0,3.409d0,3.391d0,3.374d0,3.355d0,
c                Lu      Hf      Ta      W       Re      Os      Ir
     $         3.640d0,3.141d0,3.170d0,3.069d0,2.954d0,3.120d0,2.840d0,
c                Pt      Au      Hg      Tl      Pb      Bi      Po
     $         2.754d0,3.293d0,2.705d0,4.347d0,4.297d0,4.370d0,4.709d0,
c                At      Rn      Fr      Ra      Ac      Th      Pa
     $         4.750d0,4.765d0,4.90d0, 3.677d0,3.478d0,3.396d0,3.424d0,
c                U       Np      Pu      Am      Cm      Bk      Cf
     $         3.395d0,3.424d0,3.424d0,3.381d0,3.326d0,3.339d0,3.313d0,
c                Es      Fm      Md      No      Lw
     $         3.299d0,3.286d0,3.274d0,3.248d0,3.236d0/
C
C
      IF(ITyp.EQ.1) THEN
C
C  stretch
C
       RIJ = SQRT( (XC(1,I) - XC(1,J))**2 +
     $             (XC(2,I) - XC(2,J))**2 +
     $             (XC(3,I) - XC(3,J))**2 )
c
       IAT = IAN(I)
       JAT = IAN(J)
c
       CIJ = VdW(IAT) + VdW(JAT)
c
       rho = EXP(1.0d0-(RIJ/CIJ))
c
       ForceCnst1 = 0.40d0*rho
cc
      ELSE IF(ITyp.EQ.2) THEN
C
C  bend   (warning - K is the central atom)
C
       RIK = SQRT( (XC(1,I) - XC(1,K))**2 +
     $             (XC(2,I) - XC(2,K))**2 +
     $             (XC(3,I) - XC(3,K))**2 )
       RKJ = SQRT( (XC(1,K) - XC(1,J))**2 +
     $             (XC(2,K) - XC(2,J))**2 +
     $             (XC(3,K) - XC(3,J))**2 )
c
       IAT = IAN(I)
       KAT = IAN(K)
       JAT = IAN(J)
c
       CIK = VdW(IAT) + VdW(KAT)
       CKJ = VdW(KAT) + VdW(JAT)
c
       rho1 = EXP(1.0d0-(RIK/CIK))
       rho2 = EXP(1.0d0-(RKJ/CKJ))
c
       ForceCnst1 = 0.20d0*rho1*rho2
cc

      ELSE IF(ITyp.EQ.3) THEN
C
C  out-of-plane-bend
C
       ForceCnst1 = 0.045d0
cc
      ELSE IF(ITyp.EQ.4) THEN
C
C  torsion
C
       IF(NCon.EQ.0) THEN
        RIJ = SQRT( (XC(1,I) - XC(1,J))**2 +
     $              (XC(2,I) - XC(2,J))**2 +
     $              (XC(3,I) - XC(3,J))**2 )
        RJK = SQRT( (XC(1,J) - XC(1,K))**2 +
     $              (XC(2,J) - XC(2,K))**2 +
     $              (XC(3,J) - XC(3,K))**2 )
        RKL = SQRT( (XC(1,K) - XC(1,L))**2 +
     $              (XC(2,K) - XC(2,L))**2 +
     $              (XC(3,K) - XC(3,L))**2 )
c
        IAT = IAN(I)
        JAT = IAN(J)
        KAT = IAN(K)
        LAT = IAN(L)
c
        CIJ = VdW(IAT) + VdW(JAT)
        CJK = VdW(JAT) + VdW(KAT)
        CKL = VdW(KAT) + VdW(LAT)
c
        rho1 = EXP(1.0d0-(RIJ/CIJ))
        rho2 = EXP(1.0d0-(RJK/CJK))
        rho3 = EXP(1.0d0-(RKL/CKL))
c
        ForceCnst1 = 0.01d0*rho1*rho2*rho3
       ELSE
        ForceCnst = 0.045d0    ! no small force constants with constraints
       ENDIF
cc
      ELSE IF(ITyp.EQ.7) THEN
C
C  inverse-power distance
C
       ForceCnst1 = 0.2d0
cc
      ELSE
C
C  linear bends
C
       ForceCnst1 = 0.2d0
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE FormD(NCycle, NC,     N,      NegReq, mode,
     $                 IPRNT,  U,      EigVal, FX,     VMODE,
     $                 D,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculate the new RFO step for standard optimization
C  MIN Search:  Forms a step by simple RFO that attempts to
C               minimize along all Hessian modes
C  TS Search:   Forms a step by P-RFO that attempts to maximize
C               along the direction of a chosen Hessian mode
C               (stored in VMODE) and minimize along all other modes
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C  NC      -  number of modes used to form new step
C  N       -  actual number of coordinates
C             (i.e. full dimension of problem)
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
C  D       -  new coordinate displacement (next step)
C  IErr    -  error flag    0 - new step calculated successfully
C                          -1 - something went wrong
C
C
      REAL*8 U(N,*),EigVal(NC),FX(NC),VMODE(N),D(N)
C
      REAL*8 Lamda,Lamda0,Lamda1,Lamda2
C
      PARAMETER (ZERO=0.0d0,HALF=0.5d0,FOUR=4.0d0)
      PARAMETER (TOLL=1.0d-8,STEP=0.05d0,BIG=1.0d+3,MAXIT=999)
C
C
      IOut = ioutfil('iout')
      IErr = -1
c
      NUMIT = 0
      IT = 0
      CALL ZeroIT(D,N)
C
      If(NegReq.EQ.0) GO TO 10
C
C  (a)  Maximize along selected Hessian mode
C
      IF(mode.NE.0) THEN
cc
       CALL ModeOverlap(NCycle, N,      mode,   Nmode,  IPRNT,
     $                  U,      VMODE)
C
C  On return from ModeOverlap, Nmode is the new mode along which
C  the energy is to be maximized
C
       If(IPRNT.GT.1.AND.Nmode.NE.mode) WRITE(IOut,1000) mode,Nmode
       mode = Nmode
       If(IPRNT.GT.1) WRITE(IOut,1100) mode
       IT = mode
C
C  If the mode now being followed is the lowest mode
C  switch off mode following
C
       IF(mode.EQ.1) THEN
        mode = 0
        If(IPRNT.GT.1) WRITE(IOut,1200)
       ENDIF
cc
      ELSE
cc
       If(IPRNT.GT.1) WRITE(IOut,1300)
       IT = 1
cc
      ENDIF
C
C  Calculate the Shift
C
      Lamda0 = EigVal(IT) + SQRT( EigVal(IT)**2 + FOUR*FX(IT)*FX(IT) )
      Lamda0 = HALF*Lamda0
C
      If(IPRNT.GT.1) WRITE(IOut,1400) Lamda0
      If(NC.EQ.1) GO TO 50
C
C
C  (b)  Minimize along all other Hessian modes
C
 10   CONTINUE
      JT = 1+IT
      If(JT.GT.2) JT = 1
C
      IF(IPRNT.GT.1) THEN
       If(NegReq.EQ.0) WRITE(IOut,1500)
       If(NegReq.EQ.1) WRITE(IOut,1600)
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
 20   CONTINUE
      NUMIT = NUMIT+1
      TEMP = ZERO
      DO 30 I=1,NC
      If(I.EQ.IT) GO TO 30
      TEMP = TEMP + ( FX(I)*FX(I) )/( Lamda-EigVal(I) )
 30   CONTINUE
      if(iprnt.gt.4) write(IOut,1111) lamda,temp
C
C  Check for Convergence of Lamda
C
      If(Abs(Lamda-TEMP).LT.TOLL) GO TO 40
C
C  Check for Maximum Iterations Exceeded
C
      If(NUMIT.GT.MAXIT) GO TO 96
C
C
      IF(EigVal(JT).GT.ZERO) THEN
cc
C  (i)  Simple Iterative Scheme
C
       Lamda = TEMP
       GO TO 20
cc
      ELSE
cc
C  (ii) Cautious Bracketing Scheme
C
       If(TEMP.LT.Lamda) Lamda1 = Lamda
       If(TEMP.GT.Lamda) Lamda2 = Lamda
       If(Lamda2.GT.-BIG) Lamda = HALF*(Lamda1+Lamda2)
       If(Lamda2.EQ.-BIG) Lamda = Lamda-STEP
       GO TO 20
cc
      ENDIF
C
C  At this point we should have an acceptable value for Lamda
C  Make final check
C
 40   CONTINUE
      If(IPRNT.GT.1) WRITE(IOut,1400) Lamda
      If(Lamda.GT.EigVal(JT)) GO TO 97
      If(EigVal(JT).GT.ZERO.AND.Lamda.GT.ZERO) GO TO 97
C
C  Everything OK
C  Calculate the new Step
C
 50   CONTINUE
      DO 70 I=1,NC
      TEMP = FX(I)/( Lamda-EigVal(I) )
      If(I.EQ.IT) TEMP = FX(I)/( Lamda0-EigVal(I) )
      DO 60 J=1,N
      D(J) = D(J) + TEMP*U(J,I)
 60   CONTINUE
 70   CONTINUE
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
     $       ' Lowest mode')
 1400 FORMAT(' Value Taken    Lambda = ',F12.8)
 1500 FORMAT(' Searching for Lambda that Minimizes Along All modes')
 1600 FORMAT(' Searching for Lambda that Minimizes Along All',
     $       ' other modes')
 1700 FORMAT(//,' ****************************************',/,
     $          ' ** UNABLE TO DETERMINE Lambda IN FormD **',/,
     $          ' ****************************************',//)
 1800 FORMAT(//,' *****************************************',/,
     $          ' ** ERROR IN DETERMINING Lambda IN FormD **',/,
     $          ' *****************************************',//)
c
 1111 format(' in iterative cycle:  Lambda = ',f12.8,' TEMP = ',f12.8)
c
      END
c =====================================================================
c
      SUBROUTINE FormNR(NC,N,U,EigVal,FX,D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Take the simple Newton-Raphson step
C
C  ARGUMENTS
C
C  NC      -  number of modes used to form new step
C  N       -  actual number of coordinates
C             (i.e. full dimension of problem)
C  U       -  Hessian eigenvectors
C  EigVal  -  Hessian eigenvalues
C  FX      -  gradient in Hessian eigenvector basis
C             (if mode is non-zero)
C  D       -  new coordinate displacement (next step)
C
C
      REAL*8 U(N,*),EigVal(NC),FX(NC),D(N)
C
      CALL ZeroIT(D,N)
      DO 20 I=1,NC
      TEMP = FX(I)/EigVal(I)
      DO 10 J=1,N
      D(J) = D(J) - TEMP*U(J,I)
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GDIIS(NDEG,   N,      MaxDiis,NDiis,  IPRNT,
     $                 U,      EigVal, HINV,   X,      G,
     $                 XSTR,   GSTR,   B,      ErrVec, C,
     $                 V1,     V2,     D )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  GDIIS - Variable Metric
C  See  P.Csaszar and P.Pulay  J.Mol.Struct.  114 (1984) 31
C  "Switched On" when RMS gradient is below a specified
C  tolerance - see calling routine
C
C  ARGUMENTS
C
C  NDEG    -  number of degrees of freedom
C  N       -  actual number of coordinates used
C  MaxDiis -  maximum allowed size of GDIIS subspace
C             (MaxDiis > 1; best if MaxDiis is around 4)
C  NDiis   -  current size of GDIIS subspace
C  IPRNT   -  flag controlling level of printout
C  U       -  Hessian eigenvectors
C  EigVal  -  Hessian eigenvalues
C  HINV    -  scratch space for Hessian inverse
C  X       -  current geometry
C  G       -  current gradient
C  XSTR    -  storage of all previous geometries in GDIIS subspace
C  GSTR    -  storage of all previous gradients  in GDIIS subspace
C  B       -  scratch space for GDIIS error matrix
C  ErrVec  -  scratch space for GDIIS error vectors
C  C       -  scratch space for GDIIS coefficient vector
C             (MAY USE SAME STORAGE AS ErrVec(1,1) - see calling routine)
C  V1      -  scratch vector (dimension MAX(N,MaxDiis+1))
C  V2      -  scratch vector (dimension MAX(N,MaxDiis+1))
C  D       -  on exit contains new step (unscaled)
C
C
      REAL*8 U(N,NDEG),EigVal(NDEG),HINV(N,N),X(N),G(N),
     $       XSTR(N,MaxDiis),GSTR(N,MaxDiis),B((MaxDiis+1)**2),
     $       ErrVec(N,MaxDiis),C(MaxDiis+1),V1(*),V2(*),D(N)
C
      PARAMETER (RMIN=0.001d0,RMAX=0.3d0,dtoll=1.0d-6)
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      IOut = ioutfil('iout')
C
C  On first entry simply increment counter and save current
C  geometry and gradient
C
      IF(NDiis.EQ.0) THEN
       CALL CpyVEC(N,X,XSTR(1,1))
       CALL CpyVEC(N,G,GSTR(1,1))
       NDiis = 1
       If(IPRNT.GT.1) WRITE(IOut,1000)
       RETURN
      ENDIF
c
      IDiis = NDiis
      IF(NDiis.GE.2) THEN
C
C  Check previous points relative to current point
C  if too far away or too near, reject
C
       ICount = 0
       I = 1
 10    CONTINUE
c
       DO 20 J=1,N
       V1(J) = XSTR(J,I) - X(J)
 20    CONTINUE
       RLEN = SQRT(SProd(N,V1,V1))
       If(IPRNT.GT.2) WRITE(IOut,1100) ICount+1,RLEN
c
       IF(RLEN.GT.RMAX.OR.RLEN.LT.RMIN) THEN
       If(IPRNT.GT.2) WRITE(IOut,1200)
C
C  reject old point
C  move remaining points down stack
C
        DO 30 J=I+1,NDiis
        CALL CpyVEC(N,XSTR(1,J),XSTR(1,J-1))
        CALL CpyVEC(N,GSTR(1,J),GSTR(1,J-1))
 30     CONTINUE
        I = I-1
        IDiis = IDiis-1
       ENDIF
c
       I = I+1
       ICount = ICount+1
       If(ICount.LT.(NDiis-1)) GO TO 10
      ENDIF
C
      NDiis = NDiis+1
      If(IDiis.LT.(NDiis-1)) NDiis = IDiis+1
C
C  Now we have checked all old points
C  make sure we have not exceeded maximum allowed size of
C  DIIS subspace. If so, remove the oldest point
C
      IF(NDiis.GT.MaxDiis) THEN
C
C  delete first stored geometry and gradient
C  and shift all others down
C
       NDiis = MaxDiis
       DO 40 I=2,NDiis
       CALL CpyVEC(N,XSTR(1,I),XSTR(1,I-1))
       CALL CpyVEC(N,GSTR(1,I),GSTR(1,I-1))
 40    CONTINUE
      ENDIF
C
C  store the current point at the top of the stack
C
      CALL CpyVEC(N,X,XSTR(1,NDiis))
      CALL CpyVEC(N,G,GSTR(1,NDiis))
C
C  Now start GDIIS proper
C  construct the inverse Hessian using the NDEG selected
C  eigenvectors and eigenvalues
C
      If(IPRNT.GT.1) WRITE(IOut,1300) NDiis
C
      CALL ZeroIT(HINV,N*N)
      DO 60 K=1,NDEG
      EVal = One/EigVal(K)
      DO 60 J=1,N
      SVal = U(J,K)*EVal
      DO 50 I=1,N
      HINV(I,J) = HINV(I,J) + U(I,K)*SVal
 50   CONTINUE
 60   CONTINUE
C
C  calculate the error vectors
C
      DO 70 I=1,NDiis
      CALL MatVEC(N,HINV,GSTR(1,I),ErrVec(1,I))
 70   CONTINUE
C
C  form the error matrix
C
      IDiis = NDiis+1
      DO 80 I=1,NDiis
      II = (I-1)*IDiis
      DO 80 J=1,I
      JJ = (J-1)*IDiis
      IJ = II+J
      JI = JJ+I
      B(IJ) = SProd(N,ErrVec(1,I),ErrVec(1,J))
      B(JI) = B(IJ)
 80   CONTINUE
      II = (IDiis-1)*IDiis
      DO 81 J=1,NDiis
      IJ = II+J
      JI = J*IDiis
      B(IJ) = One
      B(JI) = One
 81   CONTINUE
      B(IDiis*IDiis) = Zero
C
C  initialize the C vector
C
      CALL ZeroIT(C,NDiis)
      C(IDiis) = One
C
C  now solve the linear equations  BX = C
C  solution will be returned in C
C
      CALL LUDCMP(B,IDiis,V1,V2,ID,IErr)
      IF(IErr.NE.0) THEN
       NDiis = 0
       If(IPRNT.GT.1) WRITE(IOut,1400)
       RETURN
      ENDIF
      CALL LUBKSB(B,IDiis,V2,C)
C
      Sum = Zero
      DO 90 I=1,NDiis
      Sum = Sum + C(I)
 90   CONTINUE
      If(Abs(SUM-One).GT.dtoll) WRITE(IOut,1500) Sum
C
C  generate an intermediate geometry and gradient as a linear
C  combination of previous geometries and gradients
C
      CALL ZeroIT(V1,N)
      CALL ZeroIT(V2,N)
      DO 100 I=1,NDiis
      DO 100 J=1,N
      V1(J) = V1(J) + C(I)*XSTR(J,I)
      V2(J) = V2(J) + C(I)*GSTR(J,I)
 100  CONTINUE
C
C  now generate the new geometry
C
      CALL MatVEC(N,HINV,V2,D)
      DO 110 I=1,N
      V2(I) = V1(I) - D(I)
 110  CONTINUE
C
C  The new geometry is in V2
C  However we need a displacement vector in D to check
C  the stepsize and for possible Hessian update
C
      DO 120 I=1,N
      D(I) = V2(I) - XSTR(I,NDiis)
 120  CONTINUE
C
      RETURN
c
 1000 FORMAT('** Switching on GDIIS **')
 1100 FORMAT(' Previous point: ',I2,'  Distance from Current point: ',
     $         F12.8)
 1200 FORMAT(' point will be rejected')
 1300 FORMAT(' Minimizing Using GDIIS   Size of DIIS Subspace: ',I3)
 1400 FORMAT('**WARNING** Error Matrix Singular - Switching Off GDIIS')
 1500 FORMAT('**WARNING** Sum of GDIIS Coefficients Differs from',
     $       ' Unity',/,'  Sum is ',F12.8)
c
      END
c =====================================================================
c
      SUBROUTINE GDrive(NDrive, NPrim,  NS,     IPrnt,  LDRIVE,
     $                  FDRIVE, UT,     GINT,   GPRIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Coordinate Driving
C  This routine adds an "external force" to specified primitive internal
C  coordinates so as to drive that variable to a potential new structure.
C  (This is the internal coordinate extension of an applied external force
C   between pairs of atoms in Cartesian coordinates.)
C
C  ARGUMENTS
C
C  NDrive  -  number of separate primitives to drive
C  NPrim   -  total number of primitives
C  NS      -  number of delocalized internal coordinates (degrees of freedom)
C  IPrnt   -  print flag
C  LDRIVE  -  location of each primitive in list
C  FDRIVE  -  force to apply to each primitive
C  UT      -  current delocalized internal coordinate vectors
C  GINT    -  internal gradient
C  GPRIM   -  intermediate storage for gradient over primitives
C
      DIMENSION FDRIVE(NDrive),LDRIVE(NDrive)
      DIMENSION UT(NPrim,NS),GINT(NS),GPRIM(NPrim)
C
C  What we are going to do is to transform the gradient over delocalized
C  internals to primitives via GPRIM = UT*GINT, modify the gradient of the
C  NDrive selected primitives in GPRIM and then convert it back via
C  GINT = UT(t)*GPRIM
C
C  Get the primitive gradient
C
      write(6,*) ' GINT on entry is:'
      call prntmat(1,NS,1,gint)
cccccc
      DO 20 I=1,NPrim
      Val = 0.0d0
      DO 10 J=1,NS
      Val = Val + UT(I,J)*GINT(J)
 10   CONTINUE
      GPRIM(I) = Val
 20   CONTINUE
c
      If(IPrnt.GT.5) Then
        WRITE(6,*) ' gradient over primitives is:'
        DO I=1,NPrim
        WRITE(6,1000) I,GPRIM(I)
        EndDO
      EndIf
C
C  modify the primitive gradient
C
      DO I=1,NDrive
      IT = LDRIVE(I)
      GPRIM(IT) = GPRIM(IT) + FDRIVE(I)
      EndDO
cccccc
      WRITE(6,*) ' modified gradient over primitives is:'
      DO I=1,NPrim
      WRITE(6,1000) I,GPRIM(I)
      EndDO
cccccc
C
C  now transform back
C
      DO 40 I=1,NS
      Val = 0.0d0
      DO 30 J=1,NPrim
      Val = Val + UT(J,I)*GPRIM(J)
 30   CONTINUE
      GINT(I) = Val
 40   CONTINUE
cccccc
      write(6,*) ' GINT on exit is:'
      call prntmat(1,NS,1,gint)
cccccc
c
 1000 Format(1X,I4,2X,F10.6)
c
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GetCART(NAtoms, AtSymb, NDEG,   NPrim,  IGen,
     $                   XINT,   DINT,   ktyp,   klist,  Coeff,
     $                   SavTOR, UT,     NTrans, TRANS,  NEqATM,
     $                   BSkal,  IPRNT,  Scal,   XC,     NMem,
     $                   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Conversion of Internal Coordinates back into Cartesians
C  This is done iteratively according to:
C     XC(K) = XC(K-1) + B**(-1)*(XINT - XINT(K-1))
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols (for nice print out)
C  NDEG    -  number of degrees of freedom
C  NPrim   -  number of primitive internals
C  IGen    -  delocalized/natural internal flag
C              -1 denotes natural internals; otherwise delocalized
C  XINT    -  current internal coordinate values
C  DINT    -  step in internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  Coeff   -  primitive weights
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  UT      -  transformation matrix
C              i.e. which linear combination of primitive internals
C                   make up each delocalized internal coordinate
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  BSkal   -  scaling factor for inverse-distance coordinates
C  IPRNT   -  flag controlling print out
C  Scal    -  on exit    factor if step has been scaled
C  XC      -  on entry   previous Cartesian coordinates
C             on exit    new (current) Cartesian coordinates
C  ..................................................................
C  NMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  IErr    -  error flag    0 - new cartesians found
C                          -1 - unable to find new coordinates
C
C
      REAL*8 XC(3*NAtoms),XINT(NDEG),DINT(NDEG),UT(NPrim,NDEG)
      DIMENSION ktyp(NPrim),klist(4,NPrim),Coeff(NPrim),
     $          SavTOR(NPrim)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      CHARACTER*8 AtSymb(NAtoms)
      LOGICAL Cnvgd,invrF
C
      DIMENSION Z(NMem)
C
      PARAMETER (thrsh=5.0d-9,Maxit=30,One=1.0d0,Half=0.5d0)
C
C
      IOut = ioutfil('iout')
      NAT3 = 3*NAtoms
      invrF = .False.
C
C  Allocate memory for B-Matrix construction
C
      IMem = 2*NAT3*NPrim + 2*NPrim*NPrim + 2*NPrim
     $           + 3*NDEG + NAT3
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,7,'GetCART')
C
C  allocate scratch space
C
      IXC = iptr                         ! copy of old Cartesians
      IX0 = IXC + NAT3                   ! copy of old internals
      IXN = IX0 + NDEG                   ! copy of new internals
      Id0 = IXN + NDEG                   ! copy of step
      IB = Id0 + NDEG
      IXt = IB + NAT3*NPrim
      IG = IXt + NPrim
      IEig = IG + NAT3*NPrim
      ISc1 = IEig + NPrim
      ISc2 = ISc1 + NPrim*NPrim
      IEnd = ISc2 + NPrim*NPrim - iptr
c
      CALL MemCHK(IMem,IEnd,7,'GetCART')
c
      IT = 0
      JT = 0
      itor = 2
c
      If(IPRNT.GT.1) WRITE(IOut,1000)
C
C  make copy of original and new internals
C
      CALL MinusVEC(NDEG,XINT,DINT,Z(IX0))
      CALL CpyVEC(NDEG,XINT,Z(IXN))
c
 200  CONTINUE
C
C  make copy of original Cartesian coordinates and step
C
      CALL CpyVEC(NAT3,XC,Z(IXC))
      CALL CpyVEC(NDEG,DINT,Z(Id0))
      Scal = One
      QSQR = 1.0d+6
c
 100  CONTINUE
      IT = IT+1
      JT = JT+1
C
C  directly form the B-matrix over delocalized internals
C
      CALL BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $           Coeff,  SavTOR, 1,      1,      itor,
     $           BSkal,  IPRNT,  NDEG,   UT,     Z(IG),
     $           Z(IXt), IErr)
C
C  if we are using natural internals we need to symmetrize the B-matrix
C
      If(IGen.EQ.-1.AND.NTrans.GT.1)
     $   CALL BSYM(NAtoms, NPrim,  NTrans, TRANS,  NEqATM,
     $             Z(ISc1),Z(IG))
C
C  now invert the new B-Matrix
C
      CALL BINVT(NAtoms, NDEG,   IPRNT,  Z(IG),  Z(ISc1),
     $           Z(ISc2),Z(IB),  IErr)
c
      IF(IErr.NE.0) THEN
       If(Scal.EQ.One) Then
        invrF = .True.
        GO TO 50
       EndIf
       WRITE(IOut,1100)
       RETURN
      ENDIF
C
C  get the internal coordinates from the primitives
C
      CALL TranINT(NPrim,NDEG,Z(IXt),UT,Z(IEig))
C
C  check for convergence
C
      QSOld = QSQR
      CALL CnvBACK(IT,     NDEG,   XINT,   Z(IEig),thrsh,
     $             IPRNT,  QSQR,   Cnvgd)
c
      If(Cnvgd) GO TO 95
c
      IF(IT.GT.MaxIT) THEN
       WRITE(IOut,1200)
       IErr = -1
       RETURN
      ENDIF
C
C  if sum of squares displacement increases attempt to divide total
C  step into two halves and converge to each half separately
C
 50   If((QSQR.GT.QSOld.AND.JT.GT.3).OR.invrF) Then
       Scal = Half*Scal
       CALL CpyVEC(NAT3,Z(IXC),XC)
       CALL CpyVEC(NDEG,Z(Id0),DINT)
       CALL VScal(NDEG,Scal,DINT)
       CALL AddVEC(NDEG,Z(IX0),DINT,XINT)
       QSQR = 1.0d+6
       If(IPRNT.GT.2) WRITE(IOut,1300)
       invrF = .False.
       JT = 0
       GO TO 100
      EndIf
C
C  generate new Cartesian coordinates
C
      CALL NewCART(NAtoms,NDEG,Z(IB),XINT,Z(IEig),XC)
C
C  symmetrize
C
      CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IB),
     $            thrsh,  XC)
C
C  go back and get new internals
C
      GO TO 100
C
 95   CONTINUE
C
      If(Scal.NE.One) Then
C
C  only converged "part way"
C  reduce stepsize for this cycle
C
       WRITE(IOut,1400) Scal
       CALL CpyVEC(NDEG,Z(Id0),DINT)
      EndIF
c
      write(6,1500)it
C
      RETURN
c
 1000 FORMAT(/,' New Cartesian Coordinates Obtained by Inverse',
     $         ' Iteration')
 1100 FORMAT(/,2X,'***ERROR*** Unable to Invert B-Matrix')
 1200 FORMAT(/,2X,'***ERROR*** Exceeded allowed number of iterative',
     $            ' cycles in GetCART')
 1300 FORMAT(' Sum of Squares increases - halving stepsize')
 1400 FORMAT('**WARNING** Problems converging back-transformation',/,
     $       '  Stepsize reduced by factor of ',F12.6)
 1500 FORMAT(' Full backtransformation converged in',i4,' cycles')
c
      END
c =====================================================================
c
      SUBROUTINE GetCARTR(NAtoms, AtSymb, NDEG,   intcor, XINT,
     $                    DINT,   ktyp,   klist,  Coeff,  SavTOR,
     $                    UT,     NTrans, NEqATM, BSkal,  IPRNT,
     $                    XPrim,  Scal,   XC,     NMem,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Conversion of Inverse-power Distance Coordinates back into Cartesians
C  This is done by a distance matrix algorithm which currently uses
C  ALL interatomic distances.
C  (MUST SURELY be a way to do this with only unique distances)
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols (for nice print out)
C  NDEG    -  number of degrees of freedom
C  intcor  -  number of primitive internals
C  XINT    -  current values of delocalized internal coordinates
C  DINT    -  displacement (step) in delocalized internals
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  Coeff   -  primitive weights (set to unity in this routine)
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  UT      -  transformation matrix
C              i.e. which linear combination of primitive internals
C                   make up each delocalized internal coordinate
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  BSkal   -  scaling factor for inverse-distance coordinates
C  IPRNT   -  flag controlling print out
C  XPrim   -  values of primitives at start of optimization cycle
C             on exit contains new primitive values
C  Scal    -  on exit    factor if step has been scaled
C  XC      -  on exit    new (current) Cartesian coordinates
C  ..................................................................
C  NMem    -  Amount of available scratch space
C  Z       -  Scratch array
C  ..................................................................
C  IErr    -  error flag    0 - new cartesians found
C                          -1 - unable to find new coordinates
C
C
      REAL*8 XC(3,NAtoms),XINT(NDEG),DINT(NDEG),UT(intcor,NDEG)
      DIMENSION ktyp(intcor),klist(4,intcor),Coeff(intcor),
     $          SavTOR(intcor),XPrim(intcor),NEqATM(NAtoms,NTrans)
      CHARACTER*8 AtSymb(NAtoms)
      LOGICAL Cnvgd
C
      DIMENSION Z(NMem)
C
      PARAMETER (thrsh=5.0d-9,Maxit=999,One=1.0d0,Half=0.5d0)
C
C
      IOut = ioutfil('iout')
C
C  Allocate memory for main Distance Matrix algorithm
C
      IXp = 1                       ! copy of old primitives
      IX0 = IXp + intcor            ! copy of old internals
      IXN = IX0 + NDEG              ! copy of new internals
      Id0 = IXN + NDEG              ! copy of step
      IXt = Id0 + NDEG
      idm = IXt + intcor
      ide = idm + NAtoms*NAtoms
      idv = ide + NAtoms
      idw = idv + NAtoms*NAtoms
      iw  = idw + 32*NAtoms
      iif = iw  + 5*NAtoms
      IEnd = iif + NAtoms - 1
c
      CALL MemCHK(NMem,IEnd,8,'GetCARTR')
C
C  set all weights to unity
C
      DO 5 I=1,intcor
      Coeff(I) = One
 5    CONTINUE
c
      IT = 0
      JT = 0
      itor = 2
      itree = 1
c
      If(IPRNT.GT.1) WRITE(IOut,1000)
C
C  make copy of original and new internals
C
      CALL MinusVEC(NDEG,XINT,DINT,Z(IX0))
      CALL CpyVEC(NDEG,XINT,Z(IXN))
c
 200  CONTINUE
C
C  make copies of original primitives and step
C
      CALL CpyVEC(intcor,XPrim,Z(IXp))
      CALL CpyVEC(NDEG,DINT,Z(Id0))
      Scal = One
      QSQR = 1.0d+6
c
 100  CONTINUE
      IT = IT+1
      JT = JT+1
C
C  get the primitives from the delocalized internals
C
      DO 20 J=1,NDEG
      Val = DINT(J)
      DO 10 I=1,intcor
      XPrim(I) = XPrim(I) + UT(I,J)*Val
 10   CONTINUE
 20   CONTINUE
C
C  convert inverse-power distances to regular distances
C
      DO 30 I=1,intcor
      Z(IXt+I-1) = One/XPrim(I)
 30   CONTINUE
C
C  now get Cartesians from distance matrix
C
      CALL RIJ2CART(NAtoms, intcor, klist,  Z(IXt), Z(idm),
     $              Z(ide), Z(idv), Z(idw), Z(iw),  Z(iif),
     $              XC)
C
C  from Cartesians, generate new primitives
C
      CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $             Coeff,  SavTOR, 1,      0,      itor,
     $             BSkal,  IPRNT,  BJnk,   XPrim,  IErr)
C
C  get the internal coordinates from the primitives
C
      CALL TranINT(intcor,NDEG,XPrim,UT,Z(IXt))
C
C  check for convergence
C
      QSOld = QSQR
      CALL CnvBACK(IT,     NDEG,   XINT,   Z(IXt), thrsh,
     $             IPRNT,  QSQR,   Cnvgd)
C
C  if sum of squares displacement increases attempt to divide total
C  step into two halves and converge to each half separately
C
      If(QSQR.GT.QSOld.AND.JT.GT.3) Then
       Scal = Half*Scal
       CALL CpyVEC(intcor,Z(IXp),XPrim)
       CALL CpyVEC(NDEG,Z(Id0),DINT)
       CALL VScal(NDEG,Scal,DINT)
       CALL AddVEC(NDEG,Z(IX0),DINT,XINT)
       CALL CpyVEC(NDEG,DINT,Z(Id0))          ! save new stepsize
       QSQR = 1.0d+6
       If(IPRNT.GT.2) WRITE(IOut,1100)
       JT = 0
       GO TO 100
      EndIf
c
      If(Cnvgd) GO TO 95
c
      IF(IT.GT.MaxIT) THEN
       WRITE(IOut,1200)
       IErr = -1
       RETURN
      ENDIF
C
C  calculate new displacement
C
      CALL MinusVEC(NDEG,XINT,Z(IXt),DINT)
C
C  go back and get new internals
C
      GO TO 100
C
 95   CONTINUE
C
      If(Scal.NE.One) Then
C
C  only converged "part way"
C  reduce stepsize for this cycle
C
       WRITE(IOut,1300) Scal
       CALL CpyVEC(NDEG,Z(Id0),DINT)
      EndIF
c
      write(6,*) ' distance matrix converged in ',it,' cycles'
      IErr = 0
C
      RETURN
c
 1000 FORMAT(/,' New Cartesian Coordinates Obtained by Distance',
     $         ' Matrix algorithm')
 1100 FORMAT(' Sum of Squares increases - halving stepsize')
 1200 FORMAT(/,2X,'***ERROR*** Exceeded allowed number of iterative',
     $            ' cycles in GetCARTR')
 1300 FORMAT('**WARNING** Problems converging back-transformation',/,
     $       '  Stepsize reduced by factor of ',F12.6)
c
      END
c =====================================================================
c
      SUBROUTINE GetCARTZ(NAtoms, AtSymb, NDEG,   intcor, IGen,
     $                    XINT,   DINT,   ktyp,   klist,  Coeff,
     $                    SavTOR, UT,     NCmp,   NP1,    INT1,
     $                    UT1,    NTrans, NEqATM, IPRNT,  XPrim,
     $                    IGEO,   IMAP,   INDX,   IOrd,   XC,
     $                    NMem,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Conversion of Internal Coordinates back into Cartesians
C  This is done by selecting appropriate primitives and
C  constructing a Z-matrix, then using standard Z-matrix
C  code to do the Z-matrix ---> Cartesian transformation
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols (for nice print out)
C  NDEG    -  number of degrees of freedom
C  intcor  -  number of primitive internals
C  IGen    -  coordinate type
C              -1 - sparse natural internal coordinates
C              anything else - delocalized internal coordinates
C  XINT    -  current values of delocalized internal coordinates
C  DINT    -  displacement (step) in delocalized internals
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  Coeff   -  primitive weights (set to unity in this routine)
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  UT      -  transformation matrix
C              i.e. which linear combination of primitive internals
C                   make up each delocalized internal coordinate
C  NCmp    -  total number of non-zero natural internal components
C  NP1     -  number of primitives in each NIC
C  INT1    -  indices for all non-zero components per NIC
C  UT1     -  non-zero NIC components
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  IPRNT   -  flag controlling print out
C  XPrim   -  values of primitives at start of optimization cycle
C             on exit contains new primitive values
C  IGEO    -  Z-matrix connectivity for back-transformation
C  IMAP    -  mapping order from original to Z-matrix order
C  INDX    -  indexing array of which primitives are in Z-matrix
C  IOrd    -  flag indicating change of order between original
C             atom order and order in Z-matrix
C               0 - no change;  1 - change
C  XC      -  on exit    new (current) Cartesian coordinates
C  ..................................................................
C  NMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  IErr    -  error flag    0 - new cartesians found
C                          -1 - unable to find new coordinates
C
C
      REAL*8 XC(3*NAtoms),XINT(NDEG),DINT(NDEG),UT(intcor,NDEG),
     $       UT1(NCmp)
      DIMENSION ktyp(intcor),klist(4,intcor),Coeff(intcor),
     $          SavTOR(intcor),XPrim(intcor),NEqATM(NAtoms,NTrans),
     $          IGEO(NAtoms,4),IMAP(NAtoms,2),INDX(NDEG)
      INTEGER NP1(NDEG),INT1(NCmp)
      CHARACTER*8 AtSymb(NAtoms)
      LOGICAL Cnvgd
      dimension xppp(intcor)
C
      DIMENSION Z(NMem)
C
      PARAMETER (thrsh=5.0d-9,Maxit=30,One=1.0d0)
C
C
      IOut = ioutfil('iout')
C
C  Allocate memory for Z-Matrix construction
C
      NZ = NAtoms
      IMem =  7*NZ + MAX(intcor,3*NAtoms)
      If(NTrans.GT.1) IMem = IMem + intcor
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,8,'GetCARTZ')
C
C  allocate scratch space
C
      IXt = iptr
      IGZ = IXt + MAX(intcor,3*NAtoms)
      IEnd = IGZ + 3*NZ
c
      If(NTrans.GT.1) Then
       IG = IEnd
       IEnd = IG + intcor
      EndIf
c
      IEnd = IEnd - iptr
c
      CALL MemCHK(IMem,IEnd,8,'GetCARTZ')
C
C  set all weights to unity
C
      DO 5 I=1,intcor
      Coeff(I) = One
 5    CONTINUE
C
C  get the primitive symmetry equivalence vector if there is symmetry
C
      If(NTrans.GT.1)
     $   CALL GetSymPrim(NAtoms, intcor, NTrans, NEqATM, ktyp,
     $                   klist,  XPrim,  Z(IXt), Z(IG))
c
      IT = 0
      itor = 2
      itree = 1
c
      If(IPRNT.GT.1) WRITE(IOut,1000)
c
 100  CONTINUE
      IT = IT+1
C
C  get the primitives from the delocalized internals
C
      IF(IGen.EQ.-1) THEN
c -- sparse algorithm using natural internals
       np = 0
       DO 20 J=1,NDEG
       Val = DINT(J)
       DO 10 I=1,NP1(J)
       np = np+1
       II = INT1(np)
       XPrim(II) = XPrim(II) + UT1(np)*Val
 10    CONTINUE
 20    CONTINUE
      ELSE
c -- normal algorithm using delocalized internals
       DO 21 J=1,NDEG
       Val = DINT(J)
       DO 11 I=1,intcor
       XPrim(I) = XPrim(I) + UT(I,J)*Val
 11    CONTINUE
 21    CONTINUE
      ENDIF
C
C  fill Z-matrix with new primitives
C
      CALL FillZ(NZ,XPrim,INDX,Z(IGZ))
C
C  now get Cartesians from Z-Matrix
C
      CALL GMETRY(NZ,Z(IGZ),IGEO,XC)
C
C  Cartesians may need reordering
C
      If(IOrd.NE.0) Then
       CALL CpyVEC(3*NZ,XC,Z(IXt))
       CALL ReOrderZ(NZ,IMAP,Z(IXt),XC)
      EndIf
C
C  from Cartesians, generate new primitives
C
      CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $             Coeff,  SavTOR, 1,      0,      itor,
     $             BJnk,   IPRNT,  BJnk,   XPrim,  IErr)
C
C  symmetrize the primitives
C
      If(NTrans.GT.1) CALL SymPRIM(intcor,Z(IG),XPrim)
C
C  get the internal coordinates from the primitives
C
      IF(IGen.EQ.-1) THEN
       CALL S_TranINT(intcor, NDEG,   XPrim,  NCmp,   NP1,
     $                INT1,   UT1,    Z(IXt))
      ELSE
       CALL TranINT(intcor,NDEG,XPrim,UT,Z(IXt))
      ENDIF
C
C  check for convergence
C
      CALL CnvBACK(IT,     NDEG,   XINT,   Z(IXt), thrsh,
     $             IPRNT,  QSQR,   Cnvgd)
c
      If(Cnvgd) GO TO 95
c
      IF(IT.GT.MaxIT) THEN
       WRITE(IOut,1200)
       IErr = -1
       RETURN
      ENDIF
C
C  calculate new displacement
C
      CALL MinusVEC(NDEG,XINT,Z(IXt),DINT)
C
C  go back and get new internals
C
      GO TO 100
C
 95   CONTINUE
      write(6,*) ' z-matrix converged in ',it,' cycles'
      IErr = 0
C
      RETURN
c
 1000 FORMAT(/,' New Cartesian Coordinates Obtained by Z-Matrix',
     $         ' Construction')
 1200 FORMAT(/,2X,'***ERROR*** Exceeded allowed number of iterative',
     $            ' cycles in GetCARTZ')
c
      END
c =====================================================================
c
      SUBROUTINE GetD(NDim,IPRNT,XC,D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms the Cartesian displacement from new and old
C  Cartesian coordinates
C
C  ARGUMENTS
C
C  NDim    -  total dimension of Cartesian space
C             (typically 3*NATOMS)
C  IPRNT   -  print flag
C  XC      -  new Cartesian coordinates
C  D       -  on entry contains old Cartesians
C             on exit contains displacement
C
      REAL*8 XC(NDim),D(NDim)
C
      DO 10 I=1,NDim
      D(I) = XC(I) - D(I)
 10   CONTINUE
c
      IF(IPRNT.GT.1) THEN
C
C  print Cartesian displacement
C
       DD = SProd(NDim,D,D)
       DD = SQRT(DD)
       IOut = ioutfil('iout')
       WRITE(IOut,1000) DD
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,' Displacement from previous Coordinates is: ',F9.6)
c
      END
c =====================================================================
c
      SUBROUTINE GetHPRIM(NAt3,   NPrim,  HESS,   GINV,   IPRNT,
     $                    VM,     HPRIM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine transforms the force constant matrix from
C  Cartesian to primitive internal coordinates
C
C     HPRIM  =  GINV HESS GINV(t)
C
C  ARGUMENTS
C
C  NAt3    -  3*number of atoms (dimension of FC)
C  NPrim   -  number of primitive internals
C  HESS    -  Hessian matrix in Cartesian coordinates
C  GINV    -  transformation matrix
C  IPRNT   -  print flag
C  VM      -  general matrix storage
C  HPRIM   -  on exit contains Hessian over primitive internals
C
C
      REAL*8 HESS(NAt3,NAt3),GINV(NPrim,NAt3),VM(NPrim,NAt3),
     $       HPRIM(NPrim,NPrim)
C
C
      If(IPRNT.GT.1) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
      EndIf
C
C  first half of transformation
C
      CALL ZeroIT(VM,NPrim*NAt3)
      DO 20 J=1,NAt3
      DO 20 K=1,NAt3
      Val = HESS(K,J)
      DO 10 I=1,NPrim
      VM(I,J) = VM(I,J) + GINV(I,K)*Val
 10   CONTINUE
 20   CONTINUE
C
C  second half of transformation
C
      CALL ZeroIT(HPRIM,NPrim*NPrim)
      DO 40 K=1,NAt3
      DO 40 J=1,NPrim
      Val = GINV(J,K)
      DO 30 I=1,NPrim
      HPRIM(I,J) = HPRIM(I,J) + VM(I,K)*Val
 30   CONTINUE
 40   CONTINUE
C
C  see what we've got?
C
      If(IPRNT.GT.5) Then
       WRITE(IOut,1100)
       CALL PrntMAT(NPrim,NPrim,NPrim,HPRIM)
      EndIf
C
      RETURN
c
 1000 FORMAT(/,' Transforming Cartesian Hessian to Primitive',
     $         ' Internal Coordinates')
 1100 FORMAT(/,'  Primitive Force Constant Matrix')
c
      END
c =====================================================================
c
      SUBROUTINE GetINTL(NCycle, NAtoms, AtSymb, XC,     NMol,
     $                   IMOL,   GROUP,  NDEG,   NCons,  ICTYP,
     $                   RCON,   ICON,   LCON,   ICNUM,  NCon,
     $                   NComp,  NPComp, ICOMP,  PCWght, IPCNUM,
     $                   NDrive, LDRIVE, FDRIVE, NTrans, TRANS,
     $                   NEqATM, Steep,  IPRNT,  IHess,  IUpDat,
     $                   IDB,    GC,     HESS,   NPrim,  NPrim1,
     $                   MPrim,  IGen,   BSkal,  ktyp,   klist,
     $                   Coeff,  SavTOR, NS,     NCmp,   NP1,
     $                   INT1,   UT,     XPrim,  HPRIM,  UT0,
     $                   INDX,   Change, IOrd,   IHPrim, GOld,
     $                   DINT,   XINT,   GINT,   HINT,   NMem,
     $                   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  -----------------------------------------------------------------------
C    I N T E R N A L   C O O R D I N A T E   T R A N S F O R M A T I O N
C  -----------------------------------------------------------------------
C
C  This routine is responsible for transforming Cartesian quantities
C  into their equivalents over internal coordinates
C
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates (in cms and symmetry alligned)
C  NMol    -  number of molecules (for, e.g., cluster optimizations)
C  IMOL    -  pointers to start/end of molecules in XC array
C  GROUP   -  molecular point group
C  NDEG    -  number of internal degrees of freedom
C  NCons   -  number of constraints (fixed primitives)
C  ICTYP   -  list of constraint types
C  ICON    -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       angle constraint
C               IC1-IC2-IC3-IC4   dihedral constraint
C  RCON    -  constraint values (in atomic units)
C  LCON    -  on exit indicates which constraints are "active" (independent)
C              LCON(I) = 0  constraint is active
C              LCON(I) = 1  constraint should be eliminated
C  ICNUM   -  on exit list of which primitives are to be fixed
C  NCon    -  on exit contains actual number of constraints
C  NComp   -  number of composite constraints
C  NPComp  -  number of primitives in composite constraints
C  ICOMP   -  number of primitives in each composite constraint
C  PCWght  -  weight of each primitive in composite constraint
C  IPCNUM  -  on exit list of primitives in composite constraints
C  NDrive  -  number of primitives to drive
C  LDRIVE  -  which primitives to drive
C  FDRIVE  -  force to apply to each primitive
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  Steep   -  logical flag for steepest descent step
C  IPRNT   -  flag controlling print out
C  IHess   -  flag controlling Hessian transformation
C               0 - convert cartesian Hessian to internals
C              -1 - set up default diagonal Hessian matrix
C               1 - update internal coordinate Hessian
C  IUpDat   -  Hessian update flag
C               0 - do not update the Hessian (!?)
C              -1 - use default update
C               1 - Murtagh-Sargent update
C               2 - Powell update
C               3 - Powell/Murtagh-Sargent update (default for TS)
C               4 - BFGS update (default for minimization)
C               5 - BFGS with safeguards to ensure retention of
C                   positive definiteness
C  IDB     -  integer flag controlling dB/dcart term
C               0 - do not include this term in Hessian transformation
C                   (for approximate Hessian & retention of Cartesian
C                    Hessian eigenvalue structure)
C               1 - include it
C  GC      -  current Cartesian gradient
C  HESS    -  Hessian matrix in cartesians
C  NPrim   -  size of primitive space on current cycle
C  NPrim1  -  size of primitive space on previous cycle
C  MPrim   -  maximum size of (old+new) primitive space
C  IGen    -  delocalized/natural internal flag
C              -1 denotes natural internals; otherwise delocalized
C  BSkal   -  scale factor for inverse-distance coordinates
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
C  Coeff   -  weighting of primitives in non-redundant
C             natural internal coordinate space
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  NS      -  dimension of final optimization space
C  NCmp    -  total number of non-zero natural internal components
C  NP1     -  number of primitives in each NIC
C  INT1    -  indices for all non-zero components per NIC
C  UT      -  non-zero NIC components OR full set of
C             delocalized internal coordinates
C  XPrim   -  values of primitive internals
C  HPRIM   -  Hessian matrix over primitive internals
C             (only required if IGen = 1 or 2)
C  UT0     -  delocalized internals on previous cycle
C  INDX    -  array relating previous and current primitive sets
C             (specifically what to do to first set to get second)
C                0 - primitive common to both sets
C               -1 - primitive needs to be deleted from first set
C               +1 - primitive needs to be added to first set
C  Change  -  logical flag indicating change in primitive space
C  IOrd    -  flag for Z-matrix back-transformation
C              >=0 - Z-matrix back-transformation will be done
C               -1 - no Z-matrix back-transformation
C                    (can tidy up primitive space)
C  IHPrim  -  flag for primitive Hessian update
C                0 - use default when estimating new diagonal element
C                1 - use unity for new diagonal element
C  GOld    -  old gradient in internals
C  DINT    -  previous step in internals
C  XINT    -  on exit contains non-redundant internal coordinates
C  GINT    -  on exit contains gradient in internals
C  HINT    -  on exit contains Hessian in internals
C  ..................................................................
C  NMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  IErr    -  error flag    0 - full set of internals found
C                          -1 - unable to find full set
C
C
      DIMENSION XC(3,NAtoms),GC(3,NAtoms)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION ktyp(NPrim),klist(4,NPrim),Coeff(NPrim),SavTOR(NPrim),
     $          UT(*),XPrim(NPrim),HPRIM(NPrim**2)
      INTEGER NP1(NDEG),INT1(12*NDEG)
      Dimension UT0(3*NAtoms*NPrim1)
      DIMENSION GOld(NDEG),XINT(NDEG),GINT(NDEG),DINT(NDEG),
     $          HINT(NDEG,NDEG)
      DIMENSION IG(NDEG)
      DIMENSION ICTYP(NCons),ICON(4,NCons),RCON(NCons),LCON(NCons),
     $          ICOMP(NComp),PCWght(NPComp),IPCNUM(NPComp),
     $          ICNUM(NCons),IMOL(NMol+1)
      DIMENSION LDRIVE(NDrive),FDRIVE(NDrive)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 GROUP
      LOGICAL Steep,Change
C
      DIMENSION Z(NMem)
C
      PARAMETER (One=1.0d0)
C
      DATA thrsh/1.0d-7/
C
C
      IOut = ioutfil('iout')
C
C =================================================================
C  D E L O C A L I Z E D   I N T E R N A L   C O O R D I N A T E S
C =================================================================
C
      IF(IGen.EQ.-1) GO TO 100
C
C
C  initialization
C
      If(IHess.EQ.0) thrsh = 1.0d-9
      NAT3 = 3*NAtoms
      NVib = NDEG
c
      LGen = IGen
      If(NCycle.EQ.1.OR.IGen.GT.0) Then
       LGen = 1
       If(IGen.NE.-1) Then
        NVib = MAX(NDEG,NAT3-6)
        If(GROUP.EQ.'c*v '.OR.GROUP.EQ.'d*h ') NVib=NAtoms-1
       EndIf
       NCon = NCons           ! ** NEED TO CHECK THIS **
      EndIf
C
C  ...........................................................
C    B-MATRIX CONSTRUCTION FOR PRIMITIVE INTERNALS
C
C  initialize primitive weights
C
      DO 20 I=1,NPrim
      Coeff(I) = One
 20   CONTINUE
C  ------------------------------
C
C  Allocate memory for B-Matrix construction
C
      mem1 = NAT3*NVib + 2*NPrim
      mem2 = 0
      If(IHess.EQ.0) Then
       mem1 = mem1 + NAT3*NVib
       mem2 = NAT3*NPrim + 3*NAT3*NVib + 2*NAT3*NDEG
      EndIf
      IMem = MAX(mem1,mem2)
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,6,'GetINT')
C
C  allocate scratch pointers
C
      IB   = iptr
      IBv  = IB
      If(IHess.EQ.0) IBv = IB + NAT3*NVib
      ISc1 = IBv  + NAT3*NVib
      ISc2 = ISc1 + NPrim
      IEnd = ISc2 + NPrim - iptr
      CALL MemCHK(NMem,IEnd,8,'GetINT-1')
c
      itor = 1
C
 200  CONTINUE
C
      call secund(t1)
      CALL BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $           Coeff,  SavTOR, 1,      1,      itor,
     $           BSkal,  IPRNT,  NVib,   UT,     Z(IB),
     $           XPrim,  IErr)
C
      call secund(t2)
      t1 = t2-t1
cc      write(6,*) ' CPU time for Sparse B-Matrix construction',t1
C
C  .................................................................
C    At this point the correct B matrix is in Z(IB)
C    If the Cartesian Hessian needs to be transformed, we need
C    the inverse B-matrix.
C
      If(IHess.EQ.0) Then
c
        IB1 = IBv + NAT3*NVib
        IB2 = IB1 + NVib*NVib
        IEnd = IB2 + NVib*NVib - iptr
        CALL MemCHK(NMem,IEnd,8,'GetINT-2')
c
        CALL BINVT(NAtoms, NVib,   IPRNT,  Z(IB),  Z(IB1),
     $             Z(IB2), Z(IBv), IErr)
c
      EndIF
C
C  .................................................................
C    TRANSFORMATION SECTION - CARTESIANS TO INTERNALS
C
C -- coordinates
C
       CALL TranINT(NPrim,NVib,XPrim,UT,XINT)
C
C -- O(N) iterative transformation of gradient
C
       call secund(t1)
       i1 = IBv + NAT3*NVib
       i2 = i1  + NVib*NAT3
       i3 = i2  + NVib
       i4 = i3  + NVib
c
       CALL GrdINT2(NAtoms, NVib,   Z(IB),  GC,     IPRNT,
     $              Z(i1),  Z(i2),  Z(i3),  Z(i4),  GINT,
     $              IErr)
c
       If(IErr.NE.0) RETURN
       call secund(t2)
       t1 = t2-t1
cc       write(6,*) ' CPU time for iterative gradient transformation: ',t1
C
C  .................................................................
C    REMOVE SYMMETRY-REDUNDANT DELOCALIZED INTERNALS
C    Can any composite internal coordinates be
C    eliminated by symmetry?
C
      IF(LGen.GT.0) THEN
cc
       If(NVib.NE.NDEG) Then
        CALL SymNIC(NPrim,  NVib,   NDEG,   NAT3,   thrsh,
     $              IPRNT,  UT,     GINT,   XINT,   IHess,
     $              Z(IBv), IErr)
c
        If(IErr.NE.0) RETURN
       EndIf
C
C  ********************************
C  ** ARE THERE ANY CONSTRAINTS? **
C  ********************************
C
       If(NCons.GT.0) Then
c
        IB1  = IBv + NAT3*NVib
        IB2  = IB1 + NCons*NPrim
        IEnd = IB2 + (NCons+NDEG)*NPrim - iptr
        CALL MemCHK(NMem,IEnd,8,'GetINT-3')
c
        CALL ConVEC(NDEG,   NPrim,  NCons,  ICNUM,  ICOMP,
     $              PCWght, IPCNUM, UT,     thrsh,  IPRNT,
     $              LCON,   NCon,   Z(IB1), Z(IB2))
c
        CALL SCHMIDT(NDEG,   NPrim,  NCon,   Z(IB1), thrsh,
     $               IPRNT,  Z(IB2), NS,     LCON,   UT,
     $               IErr)
c
        If(IErr.NE.0) RETURN
C
C  get new coordinate values and gradients with respect to
C  the set of Schmidt-Orthogonalized natural internals
C
        LGen = 0
        NVib = NDEG
        GO TO 200
       Else
        NS = NDEG              ! size of optimization space
       EndIf
cc
      ENDIF
C  .................................................................
C  check for coordinate driving
C
       If(NDrive.GT.0)
     $    CALL GDrive(NDrive, NPrim,  NS,     IPrnt,  LDRIVE,
     $                FDRIVE, UT,     GINT,   Z(i1))
C  ..................................................................
C
C  find the weight of each primitive in the final
C  non-redundant optimization space
C
      If(IGen.GT.0.OR.NCycle.EQ.1)
     $ CALL GetWGHT(NPrim,NS,NCons,IPRNT,UT,Coeff)
C
C
C  avoid Hessian manipulations if steepest descent
C
      If(Steep.AND.IGen.GE.0) GO TO 95       !! not sure about this
C
C  .................................................................
C    HESSIAN TRANSFORMATION/GENERATION
C
      IF(IHess.EQ.0) THEN
C
C  transformation of Cartesian Hessian to internals
C
       IBt = IBv + NAT3*NVib
       IBp = IBt + NAT3*NPrim
       IBm = IBp + NAT3*NS
       IEnd = IBm + NAT3*NS - iptr
       CALL MemCHK(NMem,IEnd,8,'GetINT-4')
c
       CALL HssINT(NAtoms, NPrim,  NDEG,   IDB,    IPRNT,
     $             XC,     GINT,   HESS,   ktyp,   klist,
     $             Coeff,  Jnk,    BSkal,  UT,     Z(IBv),
     $             Z(IBt), Z(IBp), Z(IBm), HINT)
C
C  check if we started with only a partial Cartesian Hessian
C  If so, use partial default diagonal primitive force constants
C
       ILs1 = iptr
       ILs2 = ILs1 + NAtoms
       ian  = ILs2 + NPrim
       IHp  = ian  + NAtoms
       IH2  = IHp  + NPrim**2
       IEnd = IH2  + NPrim*NDEG - iptr
       CALL MemCHK(NMem,IEnd,8,'GetINT-5')
C
C  get atomic numbers from atomic symbols
C
       CALL GetAtNo(NAtoms,AtSymb,Z(ian))
       CALL ChkHssINT(NAtoms, Z(ian), XC,     NPrim,  NDEG,
     $                HESS,   IPRNT,  ktyp,   klist,  UT,
     $                Z(ILs1),Z(ILs2),Z(IH2), Z(IHp), HINT)
cc
      ELSE IF(IHess.EQ.-1) THEN
c
       If(IPRNT.GT.2) WRITE(IOut,1300)
C
C  set up default Hessian
C  get atomic numbers from atomic symbols
C
       CALL GetAtNo(NAtoms,AtSymb,Z(ISc1))
       CALL DefHES(NDEG,   NPrim,  NCon,   ktyp,   klist,
     $             XC,     Z(ISc1),Z(ISc2),UT,     HINT)
cc
      ELSE IF(IHess.EQ.-2) THEN
c
       If(IPRNT.GT.2) WRITE(IOut,1350)
C
C  set up unit Hessian
C
       CALL SetDiagMat(NDEG,ONE,HINT)
cc
      ELSE IF(IGen.GT.0) THEN
c
       If(IPRNT.GT.2) WRITE(IOut,1400)
C
C  update primitive Hessian and transform to new non-redundant space
C
       NPM = MAX(NPrim,NPrim1)
       IP1 = iptr
       IP2 = IP1  + NPM
       IP3 = IP2  + NPM
       IP4 = IP3  + NPM
       IP5 = IP4  + NPM
       IEnd = IP5 + NPrim*NPrim1 - iptr
       CALL MemCHK(NMem,IEnd,8,'GetINT-6')
c
       CALL UTPrim(NDEG,  NPrim,  NPrim1, GOld,   DINT,
     $             GINT,  UT,     UT0,    Z(IP1), Z(IP2),
     $             Z(IP3))
c
       If(Change) CALL MODPrim(NPrim,  NPrim1, MPrim,  IHPrim, INDX,
     $                         Z(IP4), Z(IP5), Z(IP5), Z(IP1), Z(IP2),
     $                         HPRIM,  Z(IP3), HPRIM)
       If(Steep) GO TO 95
c
       CALL UpdHES(NPrim,  IUpDat, IPRNT,  Z(IP2), Z(IP3),
     $             Z(IP1), Z(IP4), Z(IP5), HPRIM)
c
       CALL TrnHPRIM(NDEG,NPrim,HPRIM,UT,Z(IP1),HINT)
cc
      ELSE IF(IHess.EQ.1) THEN
C
C  update Hessian directly in non-redundant space
C
       CALL UpdHES(NDEG,   IUpDat, IPRNT,  DINT,   GINT,
     $             GOld,   Z(ISc1),Z(ISc2),HINT)
cc
      ENDIF
c
      If(IGen.GT.0.AND.NCycle.EQ.1) Then
C
C  need to generate primitive Hessian
C
 50    CONTINUE
       If(IHess.EQ.0) Then
C
C  Cartesian Hessian available, so transform that
C
        ignv = iptr
        ib = ignv + NAt3*NPrim
        ivm = ib + NAt3*NPrim
        iv = ivm + NPrim**2
        ieig = iv + NPrim
        igg = ieig + NPrim
        IEnd = igg + NPrim**2
        CALL MemCHK(NMem,IEnd,8,'GetINT-7')
c
        CALL GenINVT(NAtoms, XC,     NPrim,  ktyp,   klist,
     $               IPRNT,  Z(ib),  Z(ivm), Z(iv),  Z(ieig),
     $               Z(igg), Z(ignv),IErr)
c
        If(IErr.NE.0) Then
         WRITE(IOut,2000)
         IHess = -1
         GO TO 50
        EndIf
c
        CALL GetHPRIM(NAt3,   NPrim,  HESS,   Z(ignv),IPRNT,
     $                Z(ib),  HPRIM)
cc
       Else If(IHess.EQ.-2) Then
c
        If(IPRNT.GT.2) WRITE(IOut,1450)
C
C  set up unit primitive Hessian
C
        CALL SetDiagMat(NPrim,ONE,HPRIM)
        IHPrim = 1    ! for primitive Hessian update
cc
       Else
        If(IPRNT.GT.1) WRITE(IOut,1500)
        CALL GetAtNo(NAtoms,AtSymb,Z(ISc1))
        CALL DefHPRIM(NPrim,  NCon,   ktyp,   klist,  XC,
     $                Z(ISc1),HPRIM)
       EndIf
      EndIf
c
      IHess = 0       ! prevent Hessian update elsewhere
C
 95   CONTINUE
C
C  ................................................................
C    TIDY UP PRIMITIVE SPACE
C
      IF(IGen.LE.0.AND.NCycle.EQ.1.AND.IOrd.EQ.-1) THEN
       CALL TidyUP(NDEG,   NPrim,  thrsh,  IPRNT,  Z(ISc1),
     $             Coeff,  ktyp,   klist,  SavTOR, UT)
      ENDIF
c
      IErr = 0
      RETURN
C
C
 100  CONTINUE
C
C =========================================================
C  N A T U R A L   I N T E R N A L   C O O R D I N A T E S
C =========================================================
C
C  Large molecula, sparse matrix algorithm for
C  Unconstrained optimization
C
      NAT3 = 3*NAtoms
C
C  ...........................................................
C    B-MATRIX CONSTRUCTION FOR PRIMITIVE INTERNALS
C
C  Allocate memory for B-Matrix construction
C
      iibm = 1
      ibm = iibm + 52*NAT3
      nb1 = ibm + 52*NAT3
      ib1 = nb1 + NAT3
      ib  = ib1 + 12*NDEG
      nb2 = ib  + 12*NDEG
      ib2 = nb2 + NDEG
      ibt = ib2 + 12*NDEG
      IEnd = ibt + 12*NDEG
      CALL MemCHK(NMem,IEnd,8,'GetINT-8')
c
      itor = 1
c
      CALL S_BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $             SavTOR, 1,      1,      itor,   IPRNT,
     $             NDEG,   NCmp,   NP1,    INT1,   UT,
     $             Z(iibm),Z(ibm), Nb,     Z(nb1), Z(ib1),
     $             Z(ib),  Z(nb2), Z(ib2), Z(ibt), XPrim,  IErr)
C
C  .................................................................
C    TRANSFORMATION SECTION - CARTESIANS TO INTERNALS
C
C -- coordinates
C
      CALL S_TranINT(NPrim,  NDEG,   XPrim,  NCmp,   NP1,
     $               INT1,   UT,     XINT)
C
C -- O(N) iterative transformation of gradient
C
      call secund(t1)
      i1 = IEnd
      i2 = i1  + NDEG
      i3 = i2  + NDEG
      IEnd = i3 + 4*NDEG + NAT3
      CALL MemCHK(NMem,IEnd,8,'GetINT-9')
c
      CALL S_GrdINT2(NAT3,   NDEG,   GC,     IPRNT,  Nb,
     $               Z(nb1), Z(ib1), Z(ib),  Z(nb2), Z(ib2),
     $               Z(ibt), Z(i1),  Z(i2),  Z(i3),  GINT,
     $               IErr)
      call secund(t2)
      t1 = t2-t1
      write(6,*) ' gradient from first transformation is:'
      do i=1,ndeg
      write(6,'(1x,i4,f13.7)') i,gint(i)
      enddo
cc      write(6,*) ' CPU time for sparse iterative gradient',
cc     $           ' transformation: ',t1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c -  second transformation
c
      call secund(t1)
      i1 = IEnd
      i2 = i1  + NAT3
      i3 = i2  + NAT3
      IEnd = i3 + 4*NAT3 + NDEG
      CALL MemCHK(NMem,IEnd,8,'GetINT-0')
c
      CALL S_GrdINT3(NAT3,   NDEG,   GC,     IPRNT,  Nb,
     $               Z(nb1), Z(ib1), Z(ib),  Z(nb2), Z(ib2),
     $               Z(ibt), Z(i1),  Z(i2),  Z(i3),  GINT,
     $               IErr)
      call secund(t2)
      t1 = t2-t1
      write(6,*) ' gradient from second transformation is:'
      do i=1,ndeg
      write(6,'(1x,i4,f13.7)') i,gint(i)
      enddo
cc      write(6,*) ' CPU time for sparse iterative gradient',
cc     $           ' transformation: ',t1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C  .................................................................
C    HESSIAN TRANSFORMATION/GENERATION
C
      IF(IHess.EQ.-1) THEN
c
       If(IPRNT.GT.2) WRITE(IOut,1300)
C
C  set up default Hessian
C  get atomic numbers from atomic symbols
C
       CALL GetAtNo(NATOMS,AtSymb,Z(i1))
       CALL S_DefHES(NDEG,   NPrim,  NCon,   ktyp,   klist,
     $               XC,     Z(i1),  Z(i2),  NCmp,   NP1,
     $               INT1,   UT,     HINT)
cc
      ELSE IF(IHess.EQ.1) THEN
C
C  update Hessian directly in non-redundant space
C
       CALL UpdHES(NDEG,   IUpDat, IPRNT,  DINT,   GINT,
     $             GOld,   Z(i1),  Z(i2),  HINT)
cc
      ELSE
c
       If(IPRNT.GT.2) WRITE(IOut,1350)
C
C  set up unit Hessian
C
       CALL SetDiagMat(NDEG,One,HINT)
cc
      ENDIF
c
      IHess = 0       ! prevent Hessian update elsewhere
      NS = NDEG
      IErr = 0
      RETURN
c
 1300 FORMAT(' Setting Up Default Hessian')
 1350 FORMAT(' Initial Hessian is a Unit Matrix')
 1400 FORMAT(' Updating Hessian Over Primitive Internals')
 1450 FORMAT(' Initial Primitive Hessian is a Unit Matrix')
 1500 FORMAT(' Setting up Default Primitive Hessian')
 2000 FORMAT('**WARNING** Problems Transforming Hessian to Primitive',
     $       ' Internals')
c
      END
c =====================================================================
c
      SUBROUTINE GenINVT(NAtoms, XC,     NPrim,  ktyp,   klist,
     $                   IPRNT,  B,      VM,     V,      EigVal,
     $                   G,      GINV,   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  -------------------------------------------------------------
C    C O N S T R U C T   G E N E R A L I Z E D   I N V E R S E
C  -------------------------------------------------------------
C
C  This routine constructs
C    G        kinetic energy matrix
C    G**-1 B  used to transform the Cartesian force constant
C             matrix to primitive internals
C
C  ARGUMENTS
C
C  IUnit   -  unit to write mass-weighted G**-1 for later use
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NPrim   -  number of primitives
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C  IPRNT   -  flag controlling print out
C  B       -  scratch space for B-matrix
C  VM      -  matrix scratch storage
C  V       -  vector scratch storage
C  EigVal  -  storage for eigenvalues of G-matrix
C  G       -  storage for B*B(t)
C
C  on exit
C
C  GINV    -  transformation matrix
C  IErr    -  error return flag    0 - success
C                                 -1 - something went wrong
C
C
      DIMENSION XC(3,NAtoms)
      DIMENSION ktyp(NPrim),klist(4,NPrim)
      REAL*8 B(3*NAtoms,NPrim),VM(NPrim,NPrim),V(*),
     $       EigVal(NPrim),G(NPrim,NPrim),GINV(NPrim,3*NAtoms)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0)
      PARAMETER(thrsh=1.0d-7)         ! zero-eigenvalue threshold
C
C
      If(IPRNT.GT.4) WRITE(6,1000)
c
      IErr = -1
      NAt3 = 3*NAtoms
C
C  form the B-Matrix
C  make sure all primitives are evaluated
C
      DO 5 I=1,NPrim
      V(I) = ONE
 5    CONTINUE
c
      CALL BMATRIX(NAtoms, XC,     NPrim,  ktyp,   klist,
     $             V,      Jnk,    0,      1,      0,
     $             ONE,    IPRNT,  B,      Jnk,    IErr)
c
      If(IErr.NE.0) RETURN
C
C  form the G matrix [B*B(t)]
C
      DO 20 I=1,NPrim
      DO 20 J=1,I
      Val = ZERO
      DO 15 K=1,NAt3
      Val = Val + B(K,I)*B(K,J)
 15   CONTINUE
      G(I,J) = Val
      G(J,I) = Val
 20   CONTINUE
c
      If(IPRNT.GT.6) Then
       WRITE(6,1100)
       CALL PrntMAT(NPrim,NPrim,NPrim,G)
      EndIf
C
C  diagonalize G
C
      CALL DIAGMAT(G,NPrim,VM,V,EigVal,JErr)
c
      IF(JErr.NE.0) THEN
       WRITE(6,1200)
       RETURN
      ENDIF
C
C  see what we've got?
C
      If(IPRNT.GT.5) Then
       WRITE(6,1300)
       WRITE(6,1400) (EigVal(I),I=1,NPrim)
       WRITE(6,1500)
       CALL PrntMAT(NPrim,NPrim,NPrim,G)
      EndIf
C
C  invert (and count) the non-zero eigenvalues
C
      NZero = 0
      DO 30 I=1,NPrim
      If(Abs(EigVal(I)).GT.thrsh) Then
       EigVal(I) = ONE/EigVal(I)
      Else
       NZero = NZero+1
      EndIf
 30   CONTINUE
c
      NVib = NPrim - NZero
C
C  check on NVib
C
      If(NVib.NE.NAt3-6) Then
       WRITE(6,1600)
       RETURN
      EndIf
C
C  construct G**-1
C
      CALL ZeroIT(VM,NPrim*NPrim)
      DO 40 J=1,NPrim
      DO 40 K=NZero+1,NPrim
      Val = G(J,K)*EigVal(K)
      DO 35 I=1,NPrim
      VM(I,J) = VM(I,J) + G(I,K)*Val
 35   CONTINUE
 40   CONTINUE
C
C  now construct the generalized inverse [G**-1 B]
C
      CALL ZeroIT(GINV,NPrim*NAt3)
      DO 50 K=1,NPrim
      DO 50 J=1,NAt3
      Val = B(J,K)
      DO 45 I=1,NPrim
      GINV(I,J) = GINV(I,J) + VM(I,K)*Val
 45   CONTINUE
 50   CONTINUE
C
C  see what we've got?
C
      If(IPRNT.GT.6) Then
       WRITE(6,1700)
       CALL PrntMAT(NAt3,NPrim,NAt3,GINV)
      EndIf
C
 1000 FORMAT(/,' Constructing Generalized Inverse')
 1100 FORMAT(/,' B*B(t) Matrix')
 1200 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize B*B(t) in',
     $            ' <GenINV>')
 1300 FORMAT(/,'   Eigenvalues of B*B(t):')
 1400 FORMAT(1X,6F12.6)
 1500 FORMAT(/,'   Eigenvectors:')
 1600 FORMAT(/,2X,'***ERROR*** G-matrix has Wrong Number of Zero',
     $              ' Eigenvalues')
 1700 FORMAT(/,'  Generalized Inverse Matrix')
c
      END
c =====================================================================
c
      SUBROUTINE GetNPrim(NAtoms,Tors,IC,NPrim)
      IMPLICIT INTEGER(A-Z)
C
C
C  Calculates in advance the number of primitive internal
C  for the current system based on the connectivity matrix
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  Tors    -  logical flag for use of primitive torsions
C  IC      -  connectivity matrix
C  NPrim   -  on exit number of primitives
C
C
      DIMENSION IC(NAtoms,NAtoms)
      LOGICAL Bend,Outp,Tors
C
      DATA Bend/.TRUE./, Outp/.FALSE./
C
      NPrim = 0
C
C
C  (a) Stretches
C
      DO 10 I=2,NAtoms
      DO 10 J=1,I-1
      If(IC(I,J).NE.0) NPrim = NPrim+1
 10   CONTINUE
c
      If(.NOT.Bend) GO TO 25
C
C  (b) Bends
C
      DO 20 I=3,NAtoms
      DO 20 J=2,I-1
      DO 20 K=1,J-1
c
      If(IC(I,J).NE.0.AND.IC(J,K).NE.0) NPrim = NPrim+1
      If(IC(J,K).NE.0.AND.IC(K,I).NE.0) NPrim = NPrim+1
      If(IC(K,I).NE.0.AND.IC(I,J).NE.0) NPrim = NPrim+1
c
 20   CONTINUE
c
 25   CONTINUE
C
      IF(.NOT.Outp) THEN            ! usually skip out-of-plane bends
       If(NAtoms.NE.4) GO TO 35
C
C -- Special Case   Planar molecule --
C    With 3 atoms attached to a central atom there will
C    be no proper torsions
C    Use out-of-plane-bend instead
C
       I = 0
       J = 0
       K = 0
       L = 0
c
       DO 27 II=1,4
       I = I + IC(1,II)
       J = J + IC(2,II)
       K = K + IC(3,II)
       L = L + IC(4,II)
 27    CONTINUE
c
       If( .NOT.(I.EQ.3.AND.J.EQ.1.AND.K.EQ.1.AND.L.EQ.1) .AND.
     $     .NOT.(I.EQ.1.AND.J.EQ.3.AND.K.EQ.1.AND.L.EQ.1) .AND.
     $     .NOT.(I.EQ.1.AND.J.EQ.1.AND.K.EQ.3.AND.L.EQ.1) .AND.
     $     .NOT.(I.EQ.1.AND.J.EQ.1.AND.K.EQ.1.AND.L.EQ.3) ) GO TO 35
      ENDIF
C
C  (c) Out-of-Plane Bends
C
      DO 30 I=4,NAtoms
      DO 30 J=3,I-1
      DO 30 K=2,J-1
      DO 30 L=1,K-1
c
      If(IC(I,J).NE.0.AND.IC(I,K).NE.0.AND.IC(I,L).NE.0) NPrim = NPrim+3
      If(IC(J,I).NE.0.AND.IC(J,K).NE.0.AND.IC(J,L).NE.0) NPrim = NPrim+3
      If(IC(K,I).NE.0.AND.IC(K,J).NE.0.AND.IC(K,L).NE.0) NPrim = NPrim+3
      If(IC(L,I).NE.0.AND.IC(L,J).NE.0.AND.IC(L,K).NE.0) NPrim = NPrim+3
c
 30   CONTINUE
c
 35   CONTINUE
c
      If(.NOT.Tors) GO TO 45
C
C  (d) Torsions
C
      DO 42 I=2,NAtoms
      DO 42 J=1,I-1
      IF(IC(I,J).NE.0) THEN
C
C  consider I-J as middle 2 atoms in proper torsion
C
       DO 41 K=1,NAtoms
       IF(IC(I,K).NE.0.AND.K.NE.J) THEN
C
C  have K-I-J
C
        DO 40 L=1,NAtoms
        If(IC(J,L).NE.0.AND.I.NE.L.AND.K.NE.L) NPrim = NPrim+1
 40     CONTINUE
       ENDIF
 41    CONTINUE
      ENDIF
c
 42   CONTINUE
c
 45   CONTINUE
C
C  NPrim is the total number of primitives that will be generated
C  HOWEVER
C   (a) Some planar bends may need to be replaced by near-linear bends
C   (b) Additional internals may need to be generated if there are constraints
C  So add some extra just in case
C
      NPrim = NPrim + 2*NAtoms
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GetSymPrim(NAtoms, NPrim,  NTrans, NEqATM, ktyp,
     $                      klist,  XPrim,  IITM,   IG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prepares the symmetry equivalence vector for the primitives
C  indicating which primitives are equivalent by symmetry
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NPrim   -  number of primitives
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  XPrim   -  primitive values
C  IITM    -  work storage (NAtoms)
C  IG      -  on exit contains primitive equivalences
C                 0 - parent primitive
C                 J - has same value as previous (Jth) primitive
C                -J - has same value, opposite sign
C
C
      DIMENSION NEqATM(NAtoms,NTrans),ktyp(NPrim),klist(4,NPrim),
     $          XPrim(NPrim),IITM(NAtoms),IG(NPrim)
C
      PARAMETER (TollZero=1.0d-6)
C
C
C  form the IITM index vector
C  (which LOWER atom each atom is equivalent to by symmetry)
C
      CALL IZeroIT(IITM,NAtoms)
      CALL IZeroIT(IG,NPrim)
c
      DO 10 IAtm=1,NAtoms
      DO 10 IOP=1,NTrans
      DO 10 JAtm=IAtm,NAtoms
      If(IITM(JAtm).EQ.0.AND.NEqATM(JAtm,IOP).EQ.IAtm) IITM(JAtm)=IAtm
 10   CONTINUE
C
C  now loop over the primitives and determine which are equivalent
C
      DO 50 I=1,NPrim
      IType = ktyp(I)
c
      IF(IType.EQ.1) THEN
C
C  bond length
C
       IAtm = klist(1,I)
       JAtm = klist(2,I)
       IEQ = IITM(IAtm)
       JEQ = IITM(JAtm)
c
       DO 20 J=I+1,NPrim
       IF(IG(J).EQ.0.AND.ktyp(J).EQ.IType.AND.XPrim(I).EQ.XPrim(J)) THEN
        II = klist(1,J)
        JJ = klist(2,J)
        IIEQ = IITM(II)
        JJEQ = IITM(JJ)
c
        If( (IEQ.EQ.IIEQ.AND.JEQ.EQ.JJEQ) .OR.
     $      (IEQ.EQ.JJEQ.AND.JEQ.EQ.IIEQ) ) IG(J)=I
c
       ENDIF
 20    CONTINUE
C
      ELSE IF(IType.EQ.2) THEN
C
C  bond angle
C
       IAtm = klist(1,I)
       JAtm = klist(3,I)
       KAtm = klist(2,I)
       IEQ = IITM(IAtm)
       JEQ = IITM(JAtm)
       KEQ = IITM(KAtm)
c
       DO 30 J=I+1,NPrim
       IF(IG(J).EQ.0.AND.ktyp(J).EQ.IType.AND.XPrim(I).EQ.XPrim(J)) THEN
        II = klist(1,J)
        JJ = klist(3,J)
        KK = klist(2,J)
        IIEQ = IITM(II)
        JJEQ = IITM(JJ)
        KKEQ = IITM(KK)
c
        If( (IEQ.EQ.IIEQ.AND.JEQ.EQ.JJEQ.AND.KEQ.EQ.KKEQ) .OR.
     $      (IEQ.EQ.KKEQ.AND.JEQ.EQ.JJEQ.AND.KEQ.EQ.IIEQ) ) IG(J)=I
c
       ENDIF
 30    CONTINUE
C
      ELSE IF(IType.EQ.4) THEN
C
C  dihedral angle
C
       IAtm = klist(1,I)
       JAtm = klist(2,I)
       KAtm = klist(3,I)
       LAtm = klist(4,I)
       IEQ = IITM(IAtm)
       JEQ = IITM(JAtm)
       KEQ = IITM(KAtm)
       LEQ = IITM(LAtm)
c
       DO 40 J=I+1,NPrim
       IF(IG(J).EQ.0.AND.ktyp(J).EQ.IType.AND.Abs(XPrim(I)).EQ.
     $          Abs(XPrim(J)) ) THEN
        II = klist(1,J)
        JJ = klist(2,J)
        KK = klist(3,J)
        LL = klist(4,J)
        IIEQ = IITM(II)
        JJEQ = IITM(JJ)
        KKEQ = IITM(KK)
        LLEQ = IITM(LL)
c
        If( (IEQ.EQ.IIEQ.AND.JEQ.EQ.JJEQ.AND.KEQ.EQ.KKEQ.AND.
     $       LEQ.EQ.LLEQ) .OR.
     $      (IEQ.EQ.LLEQ.AND.JEQ.EQ.KKEQ.AND.KEQ.EQ.JJEQ.AND.
     $       LEQ.EQ.IIEQ) ) Then
C
C  for torsions, need to check sign
C  (do this by comparing original primitives)
C
         If(Abs(XPrim(I)+XPrim(J)).LT.TollZero) Then
            IG(J) = -I
         Else
            IG(J) = I
         EndIf
c
        EndIf
       ENDIF
 40    CONTINUE
c
      ELSE
C
C  should not get here
C
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
 50   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unknown primitive type in <GetSymPrim>')
c
      END
c =====================================================================
c
      SUBROUTINE GetWGHT(intcor, NDEG,   NCon,   IPRNT,  UT,
     $                   Coeff)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine the weight of each primitive in the final
C  symmetry non-redundant optimization space
C
C  ARGUMENTS
C
C  intcor  -  number of primitive internal coordinates
C  NDEG    -  number of degrees of freedom
C             (size of optimization space)
C  NCon    -  number of constraints (fixed primitives)
C  IPRNT   -  print flag
C  UT      -  transformation matrix
C              i.e. which linear combination of primitive internals
C                   make up each compound (natural) internal coordinate
C  Coeff   -  on exit contains primitive weights
C
C
      REAL*8 UT(intcor,NDEG),Coeff(intcor)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=1.0d-7)
C
C
      DO 20 I=1,intcor
      VAL = Zero
      DO 10 J=1,NDEG
      VAL = VAL + UT(I,J)*UT(I,J)
 10   CONTINUE
      If(Abs(VAL).LT.thrsh) VAL = Zero
      Coeff(I) = VAL
 20   CONTINUE
c
      If(IPRNT.GT.4) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       DO 30 I=1,intcor
       WRITE(IOut,1100) I,Coeff(I)
 30    CONTINUE
      EndIf
C
C  make sure any fixed primitives get a weighting or else
C  the Cartesian back-transformation will fail
C
      DO 40 I=NDEG+1,NDEG+NCon
      DO 40 J=1,intcor
      If(Coeff(J).EQ.Zero.AND.Abs(UT(J,I)).GT.thrsh)
     $   Coeff(J) = One
 40   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,'  Primitive Weights in Final Non-Redundant',
     $          ' Optimization Space',/)
 1100 FORMAT('  Weight for primitive ',I5,' is  ',F10.6)
c
      END
c =====================================================================
c
      SUBROUTINE GetZINT(NZ,     NAtoms, NIC,    NVar,   IPRNT,
     $                   IHess,  ZSymb,  GEO,    IGEO,   IG,
     $                   IDB,    GC,     HESS,   RM,     XZ,
     $                   ZINT,   XINT,   GINT,   HINT,   NMem,
     $                   Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Preparatory step for Optimization in Z-Matrix internal coordinates
C  Gets Cartesian coordinates from Z-Matrix, rotates to symmetry
C  axes and converts Cartesian gradient and (optionally) Hessian
C  to internal coordinate frame
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  NAtoms  -  number of real atoms
C  NIC     -  total number of internal coordinates (3*NZ-6)
C  NVar    -  number of variables in Z-Matrix
C             i.e. parameters to be optimized
C  IPRNT   -  print flag
C  IHess   -  flag controlling Hessian transformation
C               1 - Hessian will be updated later - do nothing
C               0 - convert cartesian Hessian to internals
C              -1 - set up default diagonal Hessian matrix
C              -2 - set up unit Hessian matrix
C              -3 - Hessian in Z-Matrix internals already available
C  ZSymb   -  Z-Matrix symbols (used to sort dummy atoms)
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO    -  Z-matrix connectivity
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  IDB     -  integer flag controlling dB/dcart term
C               0 - do not include this term in Hessian transformation
C                   (for approximate Hessian & retention of Cartesian
C                    Hessian eigenvalue structure)
C               1 - include it
C  GC      -  Cartesian gradient
C  HESS    -  Cartesian Hessian
C  RM      -  rotation matrix
C  XZ      -  Cartesian coordinates of all centres in
C             Z-matrix orientation
C  ZINT    -  on exit contains full set of internal coordinate values
C  XINT    -  on exit contains Z-Matrix variables
C  GINT    -  on exit contains Z-Matrix gradient
C  HINT    -  on exit contains Z-Matrix Hessian
C  ..................................................................
C  NMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C
C
      DIMENSION GEO(NZ,3),IGEO(NZ,4),IG(NIC),ZINT(NIC)
      DIMENSION GC(3,NAtoms),HESS(3*NAtoms,3*NAtoms),RM(3,3)
      DIMENSION XZ(3,NZ),XINT(NVar),GINT(NVar),HINT(NVar,NVar)
      CHARACTER*8 ZSymb(NZ)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      DIMENSION Z(NMem)
C
      PARAMETER (One=1.0d0)
C
C
      IOut=ioutfil('iout')
      ToRAD = PI/180.0d0
C
C  allocate scratch memory
C
      IMem = 2*(9*NZ*NZ) + 6*NZ
      If(IHess.EQ.0) IMem = IMem + 2*(9*NZ*NZ)
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,7,'GetZINT')
C
C  allocate scratch space
C
      isc1 = iptr
      isc2 = isc1 + 3*NZ
      ibv  = isc2 + 3*NZ
      ibvt = ibv  + 9*NZ*NZ
c
      If(IHess.EQ.0) THEN
       ihz = ibvt + 9*NZ*NZ
       ibm = ihz  + 9*NZ*NZ
      EndIf
C
C  convert angles and dihedrals to radians
C  and distances to Bohrs
C
      DO 10 I=1,NZ
      GEO(I,1) = ANTOAU*GEO(I,1)
      GEO(I,2) = ToRAD*GEO(I,2)
      GEO(I,3) = ToRAD*GEO(I,3)
 10   CONTINUE
C
C  get cartesian coordinates
C
      CALL GMETRY(NZ,GEO,IGEO,XZ)
c
      If(IPRNT.GT.2) Then
       WRITE(IOut,1000)
       CALL PrntCAR(IOut,0,NZ,ZSymb,XZ)
      EndIf
C
C  rotate coordinates to CMS frame
C
      CALL CMS(NZ,X,Y,Z,XZ)
      CALL RotVEC(NZ,RM,XZ)
c
      If(IPRNT.GT.2) Then
       WRITE(IOut,1100)
       CALL PrntCAR(IOut,0,NZ,ZSymb,XZ)
      EndIf
C
C  get the array of internal coordinates
C
      CALL GetZVAR(NZ,NIC,GEO,ZINT)
C
C  reduce the set of internal coordinates to the set
C  of parameters to be optimized
C
      CALL RedTOPX(NIC,NVar,IG,ZINT,XINT)
C
C  construct the B-Matrix
C
      CALL BMatZ(NZ,XZ,IGEO,IPRNT,Z(ibv))
C
C  invert it
C
      CALL INVMAT(Z(ibv),3*NZ,Z(isc1),Z(isc2),Z(ibvt),IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1200)
       RETURN
      ENDIF
c
      If(IPRNT.GT.6) Then
       WRITE(IOut,1300)
       CALL PrntMAT(3*NZ,3*NZ,3*NZ,Z(ibvt))
      EndIf
C
C  .................................................................
C    TRANSFORMATION SECTION - CARTESIANS TO INTERNALS
C
C  transform the gradient
C
      CALL GrdZMAT(NZ,     NIC,    ZSymb,  GC,     Z(ibvt),
     $             Z(isc1),Z(isc2))
C
C  reduce the set of internal gradients to the set
C  of parameters to be optimized
C
      CALL RedTOPG(NIC,NVar,IG,Z(isc2),Z(isc1),GINT)
C
C  transform the Hessian if requested
C
      IF(IHess.EQ.0) THEN

       CALL HssINTZ(NZ,     NATOMS, NIC,    NVar,   IPRNT,
     $              ZSymb,  GEO,    IGEO,   IG,     XZ,
     $              Z(isc2),IDB,    HESS,   Z(ihz), Z(isc1),
     $              Z(ibvt),Z(ibv), Z(ibm), Z(ibm), HINT)

      ELSE IF(IHess.EQ.-1) THEN
cc
C  set up default diagonal Hessian
C
       CALL DefHESZ(NZ,NVar,IG,HINT)
       IHess = 0
cc
      ELSE IF(IHess.EQ.-2) THEN
cc
C  set up unit matrix
C
       CALL SetDiagMat(NVar,One,HINT)
       IHess = 0
cc
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,4X,' Cartesian Coordinates in Z-Matrix Orientation')
 1100 FORMAT(/,4X,' Cartesian Coordinates in Standard Orientation')
 1200 FORMAT(/,2X,'***ERROR*** Problems Inverting B-Matrix in',
     $              ' <GetZINT>')
 1300 FORMAT(/,6X,'Inverse B-Matrix (Z-matrix coordinates)')
c
      END
c =====================================================================
c
      SUBROUTINE GetZVAR(NZ,NIC,GEO,ZINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Get the full array of internal coordinates from
C  the array GEO
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  NIC     -  total number of internal coordinates
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C             (**note:  already converted to bohr/radians)
C  ZINT    -  on exit contains vector of internal coordinates
C
C
      REAL*8 GEO(NZ,3),ZINT(NIC)
C
C
      DO 10 I=1,NZ-1
      ZINT(I) = GEO(I+1,1)
 10   CONTINUE
c
      IT = NZ-1
c
      DO 20 I=1,NZ-2
      IT = IT+1
      ZINT(IT) = GEO(I+2,2)
 20   CONTINUE
c
      DO 30 I=1,NZ-3
      IT = IT+1
      ZINT(IT) = GEO(I+3,3)
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GrdINT(NAtoms,NDEG,GC,BINV,GINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transforms cartesian gradient to gradient over the
C  set of internal coordinates
C
      REAL*8 GC(3*NAtoms),BINV(3*NAtoms,NDEG),GINT(NDEG)
C
      DO 10 I=1,NDEG
      GINT(I) = SProd(3*NAtoms,BINV(1,I),GC)
 10   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GrdINT2(NAtoms, N,      B,      gcart,  IPRNT,
     $                   BT,     Bg,     D,      Z,      gint,
     $                   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Iterative transformation of Cartesian gradient to gradient
C  over delocalized internal coordinates
C  Solve
C     Gint(k+1)  =  D**-1 { B*Gcart + (D - B*Bt)*Gint(k) }
C  where D is the diagonal of B*Bt
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  N       -  size of optimization space
C  B       -  B-Matrix over delocalized internals
C  gcart   -  gradient in Cartesian coordinates
C  IPRNT   -  print flag
C  BT      -  storage for B(t)
C  Bg      -  storage for B*gcart
C  D       -  storage for diagonal elements of BBT
C  Z       -  general scratch storage (at least 4*N + 3*NAtoms)
C  gint    -  on exit contains gradient in delocalized internals
C  IErr    -  error flag    0 - new gradient found
C                          -1 - unable to converge to new gradient
C
C
      DIMENSION B(3*NAtoms,N),gcart(3*NAtoms),Bg(N),
     $          BT(N,3*NAtoms),D(N),Z(*),gint(N)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=5.0d-9,MaxIT=999)
C
C
C  form B*gcart
C
      DO 10 I=1,N
      Bg(I) = SProd(3*NAtoms,B(1,I),gcart)
 10   CONTINUE
C
C  form and invert diagonal elements of B*B(t)
C
      DO 30 I=1,N
      Val = Zero
      DO 20 K=1,3*NAtoms
      Val = Val + B(K,I)**2
 20   CONTINUE
      D(I) = One/Val
      gint(I) = Bg(I)*D(I)       ! first estimate of gradient
 30   CONTINUE
C
C  form B-transpose
C
      DO 40 I=1,N
      DO 40 J=1,3*NAtoms
      BT(I,J) = B(J,I)
 40   CONTINUE
c
      If(IPRNT.GT.2) Then
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
      EndIf
C
C  now call <CONJUGATE> to solve linear equations
C
      i1 = 1
      i2 = i1 + N
      i3 = i2 + N
      i4 = i3 + N
      i5 = i4 + N
c
      CALL CONJUGATE(NAtoms, N,      B,      BT,     Bg,
     $               D,      gint,   thrsh,  MaxIT,  IPRNT,
     $               Z(i1),  Z(i2),  Z(i3),  Z(i4),  Z(i5),
     $               IErr)
C
      RETURN
c
 1000 FORMAT(/,5X,'Iterative generation of Internal Gradient')
c
      END
c =====================================================================
c
      SUBROUTINE GrdZMAT(NZ,     NIC,    ZSymb,   GC,     BINV,
     $                   GZ,     ZINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transforms Cartesian gradient to full Z-Matrix internal
C  coordinate gradient
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  NIC     -  total number of internal coordinates
C  ZSymb   -  Z-Matrix symbols (used to sort dummy atoms)
C  GC      -  Cartesian gradient
C  BINV    -  inverse B-Matrix
C  GZ      -  storage for expanded Cartesian gradient
C             (with zeros added for dummy atoms)
C  ZINT    -  on exit contains Z-Matrix gradient
C
C
      REAL*8 GC(3,*),BINV(3*NZ,3*NZ),GZ(3,NZ),ZINT(NIC)
      CHARACTER*8 ZSymb(NZ)
C
      PARAMETER (Zero=0.0d0)
C
C
C  first expand cartesian gradient to include dummy atoms
C
      IAtm = 0
      DO 10 I=1,NZ
      IF(ZSymb(I)(1:2).EQ.'Du') THEN
       GZ(1,I) = Zero
       GZ(2,I) = Zero
       GZ(3,I) = Zero
      ELSE
       IAtm = IAtm+1
       GZ(1,I) = GC(1,IAtm)
       GZ(2,I) = GC(2,IAtm)
       GZ(3,I) = GC(3,IAtm)
      ENDIF
 10   CONTINUE
C
C  now generate the internal coordinate gradient
C
      CALL ZeroIT(ZINT,NIC)
      DO 30 I=1,NZ
      II = 3*(I-1)
      DO 20 J=1,NIC
      ZINT(J) = ZINT(J) + BINV(J,II+1)*GZ(1,I)
     $                  + BINV(J,II+2)*GZ(2,I)
     $                  + BINV(J,II+3)*GZ(3,I)
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
