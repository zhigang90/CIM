      SUBROUTINE SymFAC(NAtoms, IAtom,  X,      Y,      Z,
     $                  XP,     YP,     ZP,     NSym,   NGen,
     $                  IGEN,   NEqATM, SFac,   LNum,   LF)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine checks the grid point on the current symmetry-unique
C  atom, eliminating grid points which are equivalent by symmetry
C  and adjusting the quadrature weights appropriately for those
C  points that are retained
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAtom   -  current atom
C  X       -  X-coordinate of current atom
C  Y       -  Y-coordinate of current atom
C  Z       -  Z-coordinate of current atom
C  XP      -  X-coordinate of current grid point
C  YP      -  Y-coordinate of current grid point
C  ZP      -  Z-coordinate of current grid point
C  NSym    -  total number of symmetry operations
C  NGen    -  number of generators of the group
C  IGEN    -  list of group generators
C             these are highest non-abelian subgroup and are:
C               1 - reflection in YZ plane
C               2 - reflection in XZ plane
C               3 - reflection in XY plane
C               4 - C2 rotation about Z-axis
C               5 - C2 rotation about Y-axis
C               6 - C2 rotation about X-axis
C               7 - inversion through origin
C             (contains ALL symmetry operations; first NGen are generators)
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C  on exit
C
C  SFac    -  symmetry weight factor
C  LNum    -  number of symmetry operations for current grid point
C  LF      -  list of symmetry operations for current grid point
C
C
      DIMENSION IGEN(NSym),NEqATM(NAtoms,*),LF(NSym)
      INTEGER SFac
c
      Dimension LTable(0:7,0:7)
      Data LTable / 0, 1, 2, 3, 4, 5, 6, 7,
     $              1, 0, 4, 5, 2, 3, 7, 6,
     $              2, 4, 0, 6, 1, 7, 3, 5,
     $              3, 5, 6, 0, 7, 1, 2, 4,
     $              4, 2, 1, 7, 0, 6, 5, 3,
     $              5, 3, 7, 1, 6, 0, 4, 2,
     $              6, 7, 3, 2, 5, 4, 0, 1,
     $              7, 6, 5, 4, 3, 2, 1, 0 /
C
      PARAMETER (Zero=0.0d0)
C
C
C  Loop over number of generators
C  ** WARNING - ASSUMES MOLECULE IS CORRECTLY ORIENTED **
C
      SFac = 1
      LNum = 0
cc
      DO 90 IOp=1,NGen
C
C  If the symmetry operation takes atom IAtom into a new atom
C  then do nothing
C
      If(NEqATM(IAtom,IOp).NE.IAtom) GO TO 90
      ISym = IGEN(IOp)
C
C  If molecule has been oriented correctly then we can prune
C  the grid on every atom which has a zero coordinate
C
      IF(ISYM.EQ.1.AND.X.EQ.Zero) THEN
C
C  atom is in YZ plane
C
        LNum = LNum+1
        LF(LNum) = 1
C
C  eliminate grid point with negative X
C  double weight if positive X
C
        If(XP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(XP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ELSE IF(ISYM.EQ.2.AND.Y.EQ.Zero) THEN
C
C  atom is in XZ plane
C
        LNum = LNum+1
        LF(LNum) = 2
C
C  eliminate grid point with negative Y
C  double weight if positive Y
C
        If(YP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(YP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ELSE IF(ISym.EQ.3.AND.Z.EQ.Zero) THEN
C
C  atom is in XY plane
C
        LNum = LNum+1
        LF(LNum) = 3
C
C  eliminate grid point with negative Z
C  double weight if positive Z
C
        If(ZP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(ZP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ELSE IF(ISym.EQ.4.AND.X.EQ.Zero.AND.Y.EQ.Zero) THEN
C
C  atom is on Z-axis
C
        LNum = LNum+1
        LF(LNum) = 4
C
C  eliminate grid point with negative X
C  double weight if positive X
C
        If(XP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(XP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ELSE IF(ISym.EQ.5.AND.X.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on Y-axis
C
        LNum = LNum+1
        LF(LNum) = 5
C
C  eliminate grid point with negative Z
C  double weight if positive Z
C
        If(ZP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(ZP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ELSE IF(ISym.EQ.6.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on X-axis
C
        LNum = LNum+1
        LF(LNum) = 6
C
C  eliminate grid point with negative Y
C  double weight if positive Y
C
        If(YP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(YP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ELSE IF(ISym.EQ.7.AND.X.EQ.Zero.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is at origin of axis system
C
        LNum = LNum+1
        LF(LNum) = 7
C
C  eliminate grid point with negative XYZ product
C  double weight if positive product
C
        XYZ = XP*YP*ZP
        If(XYZ.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(XYZ.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc        write(6,*) ' symmetry operation: ',isym,' LNum:',lnum
cc
      ENDIF
 90   CONTINUE
C
      If(SFac.GT.1) RETURN
C
C ------------------------------------------------------------
C  WARNING - EXPERIMENTAL CODE
C
C  No pruning has been done
C  HOWEVER, there MAY be a symmetry operation of the (Abelian) group
C  that takes an atom into itself, even though the generators do not.
C  In this case we need to prune
C
      DO 91 IOp=NGen+1,NSym
      If(NEqATM(IAtom,IOp).NE.IAtom) GO TO 91
      ISym = IGEN(IOp)
C
C  If molecule has been oriented correctly then we can prune
C  the grid on every atom which has a zero coordinate
C
      IF(ISYM.EQ.1.AND.X.EQ.Zero) THEN
C
C  atom is in YZ plane
C
        LNum = LNum+1
        LF(LNum) = 1
C
C  eliminate grid point with negative X
C  double weight if positive X
C
        If(XP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(XP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ELSE IF(ISYM.EQ.2.AND.Y.EQ.Zero) THEN
C
C  atom is in XZ plane
C
        LNum = LNum+1
        LF(LNum) = 2
C
C  eliminate grid point with negative Y
C  double weight if positive Y
C
        If(YP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(YP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ELSE IF(ISym.EQ.3.AND.Z.EQ.Zero) THEN
C
C  atom is in XY plane
C
        LNum = LNum+1
        LF(LNum) = 3
C
C  eliminate grid point with negative Z
C  double weight if positive Z
C
        If(ZP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(ZP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ELSE IF(ISym.EQ.4.AND.X.EQ.Zero.AND.Y.EQ.Zero) THEN
C
C  atom is on Z-axis
C
        LNum = LNum+1
        LF(LNum) = 4
C
C  eliminate grid point with negative X
C  double weight if positive X
C
        If(XP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(XP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ELSE IF(ISym.EQ.5.AND.X.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on Y-axis
C
        LNum = LNum+1
        LF(LNum) = 5
C
C  eliminate grid point with negative Z
C  double weight if positive Z
C
        If(ZP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(ZP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ELSE IF(ISym.EQ.6.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on X-axis
C
        LNum = LNum+1
        LF(LNum) = 6
C
C  eliminate grid point with negative Y
C  double weight if positive Y
C
        If(YP.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(YP.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ELSE IF(ISym.EQ.7.AND.X.EQ.Zero.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is at origin of axis system
C
        LNum = LNum+1
        LF(LNum) = 7
C
C  eliminate grid point with negative XYZ product
C  double weight if positive product
C
        XYZ = XP*YP*ZP
        If(XYZ.LT.Zero) Then
         SFac = Zero
         RETURN
        Else If(XYZ.GT.Zero) Then
         SFac = 2*SFac
        EndIf
cc
      ENDIF
 91   CONTINUE
C
C  ----------------------------------------------------------------
C
C  The number of generators (LNum) should be either
C  1, 2 or 3 (for D2h only). For the spherical harmonic
C  stuff, the effect of multiple reflections is the
C  accumulated effect of each reflection; however for
C  rotations the effect is NOT cumulative and any two
C  rotations can be replaced by the third
C
      IF(LNum.EQ.2) THEN
        L1 = LF(1)
        L2 = LF(2)
        If(L1.GT.L2) Then
         IOp = L1
         L1 = L2
         L2 = IOp
        EndIf                    ! L1 < L2
        If(L1.EQ.4.AND.L2.EQ.5) Then
         LNum = 1
         LF(1) = 6
        Else If(L1.EQ.4.AND.L2.EQ.6) Then
         LNum = 1
         LF(1) = 5
        Else If(L1.EQ.5.AND.L2.EQ.6) Then
         LNum = 1
         LF(1) = 4
        EndIf
      ENDIF
C
C  -----------------------------------------------------------------
C
      RETURN
      END
c ......................................................................
c
      SUBROUTINE SymGenEFG(X,      Y,      Z,      XP,     YP,
     $                     ZP,     NGen,   IGEN,   GFac,   XEFG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  From the symmetry unique grid points, this routine generates all
C  symmetry equivalent grid points and folds the contribution from
C  each grid point into the electric field gradient of the current atom
C
C  ARGUMENTS
C
C  X       -  X-coordinate of current atom
C  Y       -  Y-coordinate of current atom
C  Z       -  Z-coordinate of current atom
C  XP      -  X-coordinate of current grid point
C  YP      -  Y-coordinate of current grid point
C  ZP      -  Z-coordinate of current grid point
C  NGen    -  number of generators of the group
C  IGEN    -  list of group generators
C             these are highest non-abelian subgroup and are:
C               1 - reflection in YZ plane
C               2 - reflection in XZ plane
C               3 - reflection in XY plane
C               4 - C2 rotation about Z-axis
C               5 - C2 rotation about Y-axis
C               6 - C2 rotation about X-axis
C               7 - inversion through origin
C             (contains ALL symmetry operations; first NGen are generators)
C  GFac    -  grid weight density factor
C
C  on exit
C
C  XEFG    -  accumulated electric field gradient for current atom
C
C
      DIMENSION IGEN(NGen),XEFG(9)
      PARAMETER (Three=3.0d0)
c
      Dimension LFAC(3,7)
      Data LFAC / -1, 1, 1,   1,-1, 1,   1, 1,-1,  -1,-1, 1,
     $            -1, 1,-1,   1,-1,-1,  -1,-1,-1 /
c
C  LFAC indicates how the coordinates change sign with the different
C  abelian generators
C
C  We have already accumulated the contribution from the original
C  grid point. Generate the symmetry equivalent grid points and
C  their accumulations. Note that if there are
C   1 generator  there will be 1 more grid point (2 in total)
C   2 generators    "     "    3    "       "    (4 in total)
C   3 generators    "     "    7    "       "    (8 in total)
C
C
C  First generator
C
cc      write(6,*) ' In <SymEFG>  NGen:',ngen,' IGEN: ',(igen(i),i=1,ngen)
cc      write(6,*) ' X=',x,' Y=',y,' Z=',z
cc      write(6,*) ' XP=',xp,' YP=',yp,' ZP=',zp
      LL = IGEN(1)
cc      X1 = LFAC(1,LL)*X
cc      Y1 = LFAC(2,LL)*Y
cc      Z1 = LFAC(3,LL)*Z
      XP1 = LFAC(1,LL)*XP
      YP1 = LFAC(2,LL)*YP
      ZP1 = LFAC(3,LL)*ZP
cc      write(6,*) ' New atomic coords from 1st generator:'
cc      write(6,*) ' X1=',x1,' Y1=',y1,' Z1=',z1
cc      write(6,*) ' New grid point from 1st generator:'
cc      write(6,*) ' XP1=',xp1,' YP1=',yp1,' ZP1=',zp1
c -- form the contribution
      xx = XP1 - X
      yy = YP1 - Y
      zz = ZP1 - Z
      x2 = xx**2
      y2 = yy**2
      z2 = zz**2
      rr = x2 + y2 + z2
      r1 = sqrt(rr)
cc      write(6,*) ' r1 is ',r1
      Fac = GFac/(r1**5)
      XEFG(1) = XEFG(1) + (Three*x2-rr)*Fac
      XEFG(3) = XEFG(3) + (Three*y2-rr)*Fac
      XEFG(6) = XEFG(6) + (Three*z2-rr)*Fac
      XEFG(2) = XEFG(2) + xx*yy*Fac
      XEFG(4) = XEFG(4) + xx*zz*Fac
      XEFG(5) = XEFG(5) + yy*zz*Fac
cc      XEFG(7) = XEFG(7) + x2*Fac
cc      XEFG(8) = XEFG(8) + y2*Fac
cc      XEFG(9) = XEFG(9) + z2*Fac
cc
      If(NGen.EQ.1) RETURN
C
C  Second generator
C
      LL = IGEN(2)
cc      X2 = LFAC(1,LL)*X
cc      Y2 = LFAC(2,LL)*Y
cc      Z2 = LFAC(3,LL)*Z
      XP2 = LFAC(1,LL)*XP
      YP2 = LFAC(2,LL)*YP
      ZP2 = LFAC(3,LL)*ZP
cc      write(6,*) ' New atomic coords from 2nd generator:'
cc      write(6,*) ' X2=',x2,' Y2=',y2,' Z2=',z2
cc      write(6,*) ' New grid point from 2nd generator is:'
cc      write(6,*) ' XP2=',xp2,' YP2=',yp2,' ZP2=',zp2
c -- form the contribution
      xx = XP2 - X
      yy = YP2 - Y
      zz = ZP2 - Z
      x2 = xx**2
      y2 = yy**2
      z2 = zz**2
      rr = x2 + y2 + z2
      r1 = sqrt(rr)
cc      write(6,*) ' r1 is ',r1
      Fac = GFac/(r1**5)
      XEFG(1) = XEFG(1) + (Three*x2-rr)*Fac
      XEFG(3) = XEFG(3) + (Three*y2-rr)*Fac
      XEFG(6) = XEFG(6) + (Three*z2-rr)*Fac
      XEFG(2) = XEFG(2) + xx*yy*Fac
      XEFG(4) = XEFG(4) + xx*zz*Fac
      XEFG(5) = XEFG(5) + yy*zz*Fac
cc      XEFG(7) = XEFG(7) + x2*Fac
cc      XEFG(8) = XEFG(8) + y2*Fac
cc      XEFG(9) = XEFG(9) + z2*Fac
c
cc      X3 = LFAC(1,LL)*X1
cc      Y3 = LFAC(2,LL)*Y1
cc      Z3 = LFAC(3,LL)*Z1
      XP3 = LFAC(1,LL)*XP1
      YP3 = LFAC(2,LL)*YP1
      ZP3 = LFAC(3,LL)*ZP1
cc      write(6,*) ' New atomic coords from 3rd generator:'
cc      write(6,*) ' X3=',x3,' Y3=',y3,' Z3=',z3
cc      write(6,*) ' New grid point from 3rd generator is:'
cc      write(6,*) ' XP3=',xp3,' YP3=',yp3,' ZP3=',zp3
c -- form the contribution
      xx = XP3 - X
      yy = YP3 - Y
      zz = ZP3 - Z
      x2 = xx**2
      y2 = yy**2
      z2 = zz**2
      rr = x2 + y2 + z2
      r1 = sqrt(rr)
cc      write(6,*) ' r1 is ',r1
      Fac = GFac/(r1**5)
      XEFG(1) = XEFG(1) + (Three*x2-rr)*Fac
      XEFG(3) = XEFG(3) + (Three*y2-rr)*Fac
      XEFG(6) = XEFG(6) + (Three*z2-rr)*Fac
      XEFG(2) = XEFG(2) + xx*yy*Fac
      XEFG(4) = XEFG(4) + xx*zz*Fac
      XEFG(5) = XEFG(5) + yy*zz*Fac
cc      XEFG(7) = XEFG(7) + x2*Fac
cc      XEFG(8) = XEFG(8) + y2*Fac
cc      XEFG(9) = XEFG(9) + z2*Fac
cc
      If(NGen.EQ.2) RETURN
C
C  Third generator  (D2h only)
C
      RETURN
      END
c ......................................................................
c
      SUBROUTINE RedEFG(NAtoms, ICntr,  NSym,   NEqATM, XEFG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Scale back the symmetry-generated EFG components
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  ICntr   -  current symmetry-unique atom
C  NSym    -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  XEFG    -  electric field gradient tensor
C             on exit may be scaled down
C
C
      DIMENSION ISYM(NSym),NEqATM(NAtoms,NSym),XEFG(6)
C
      Parameter (One=1.0d0)
C
C
C  EFG-components for symmetry-unique atoms have been generated
C  in <SymEFG>. What we need to do here is to divide by the
C  correct integer symmetry factor. This factor is equal to the
C  number of symmetry operations which take a symmetry-unique
C  atom into itself.
C
      NFac = 1
      DO 10 IOP=1,NSym
      If(NEqATM(ICntr,IOP).EQ.ICntr) NFac = NFac+1
 10   CONTINUE
c
      If(NFac.EQ.1) RETURN       ! nothing to do
C
C  divide by symmetry factor
C
      skal = One/DBLE(NFac)
      CALL VScal(6,skal,XEFG)
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE SymEFG(NAtoms,ICntr,NGen,ISYM,NEqATM,XEFG)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Generates electric field gradient components for symmetry-related
C  atoms from symmetry unique atoms
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  ICntr   -  current symmetry-unique centre
C  NGen    -  number of abelian generators
C  ISYM    -  list of abelian symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  XEFG    -  electric field gradient tensor
C
C
      DIMENSION ISYM(NGen),NEqATM(NAtoms,NGen),XEFG(9,NAtoms)
C
C  loop over symmetry generators
C
      DO 10 IOP=1,NGen
      JAtm = NEqATM(ICntr,IOP)
      If(JAtm.NE.ICntr) Then
        CALL GetSymFac(ISYM(IOP),facXY,facXZ,facYZ)
        XEFG(1,JAtm) = XEFG(1,ICntr)
        XEFG(3,JAtm) = XEFG(3,ICntr)
        XEFG(6,JAtm) = XEFG(6,ICntr)
        XEFG(2,JAtm) = facXY*XEFG(2,ICntr)
        XEFG(4,JAtm) = facXZ*XEFG(4,ICntr)
        XEFG(5,JAtm) = facYZ*XEFG(5,ICntr)
cc        XEFG(7,JAtm) = XEFG(7,ICntr)
cc        XEFG(8,JAtm) = XEFG(8,ICntr)
cc        XEFG(9,JAtm) = XEFG(9,ICntr)
      EndIf
 10   CONTINUE
C
      RETURN
      END
c ......................................................................
c
      SUBROUTINE GetSymFac(IOP,facXY,facXZ,facYZ)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Gets (abelian) symmetry factors for possible sign changes when
C  generating XY, XZ and YZ electric field gradient components on
C  symmetry related atoms using values on symmetry-unique atom
C
C  ARGUMENTS
C
C  IOP     -  generating symmetry operation
C
C  on exit
C
C  facXY
C  facXZ   -  symmetry factors (either +1 or -1)
C  facYZ
C
C
      Dimension LFAC(3,7)
      Data LFAC / -1, 1, 1,   1,-1, 1,   1, 1,-1,  -1,-1, 1,
     $            -1, 1,-1,   1,-1,-1,  -1,-1,-1 /
C
C  LFAC indicates how the coordinates change sign with the different
C  abelian generators
C
C
      FacXY = LFAC(1,IOP)*LFAC(2,IOP)
      FacXZ = LFAC(1,IOP)*LFAC(3,IOP)
      FacYZ = LFAC(2,IOP)*LFAC(3,IOP)
C
      RETURN
      END
