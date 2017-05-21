c ==================================================================
c  GEOMETRY OPTIMIZATION ROUTINES H-P          JB   October 1999
c ==================================================================
c
      SUBROUTINE HSSCART(NAtoms, intcor, NDEG,   IDB,    IPRNT,
     $                   XC,     GINT,   HINT,   ktyp,   klist,
     $                   Coeff,  SavTOR, BSkal,  UT,     B,
     $                   BP,     BM,     BMT,    BT,     HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Back-Transformation of Hessian in internal coordinates to
C  Cartesian Hessian according to:
C     HESS = Bt * HINT * B  + GINT*dB/dcart
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  intcor  -  number of primitive internals
C  NDEG    -  number of degrees of freedom
C             (and number of active internal coordinates)
C  IDB     -  integer flag controlling dB/dcart term
C               0 - do not include this term in Hessian transformation
C                   (for approximate Hessian & retention of Hessian
C                    eigenvalue structure)
C               1 - include it
C  IPRNT   -  flag controlling printout
C  XC      -  cartesian coordinates
C  GINT    -  gradient in internal coordinates
C  HINT    -  Hessian in internal coordinates
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C  klist   -  list of atoms involved in each primitive internal
C  Coeff   -  weighting of primitive in non-redundant
C             delocalized internal coordinate space
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  BSkal   -  scale factor for inverse-distance coordinates
C  UT      -  transformation matrix
C               i.e. which linear combination of primitive internals
C                    make up each compound (natural) internal coordinate
C  B       -  will hold B-Matrix in primitive internals
C  BP      -  will hold transformed B-Matrix forward step
C  BM      -  will hold transformed B-Matrix backward step
C             (MAY USE SAME STORAGE AS B - see calling routine)
C  BMT     -  will hold intermediate half-transformed matrix
C             (MAY USE SAME STORAGE AS BP - see calling routine)
C  BT      -  will hold transposed B-Matrix
C             (MAY USE SAME STORAGE AS B - see calling routine)
C  HESS    -  on exit Hessian in Cartesian coordinates
C
C
      REAL*8 XC(3*NAtoms),GINT(NDEG),HINT(NDEG,NDEG),
     $       UT(intcor,NDEG),B(3*NAtoms,intcor),
     $       BP(3*NAtoms,NDEG),BM(3*NAtoms,NDEG),
     $       BT(NDEG,3*NAtoms),BMT(NDEG,3*NAtoms),
     $       HESS(3*NAtoms,3*NAtoms)
C
      PARAMETER (Zero=0.0d0,TolZero=1.0d-8,One=1.0d0)
      PARAMETER (delta=0.0005d0)
C
C
      IOut = ioutfil('iout')
c
      If(IPRNT.GT.1) Then
       WRITE(IOut,1000)
       If(IDB.EQ.0) WRITE(IOut,1100)
      EndIf
c
      NAT3 = 3*NAtoms
C
C  Transform the internal coordinate Hessian
C  -----------------------------------------
C  First form the B-matrix at the current geometry
C
      CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $             coeff,  Jnk,    0,      1,      0,
     $             BSkal,  0,      B,      Jnk,    IErr)
c
      If(IErr.NE.0) GO TO 95
C
C  transform the original B-Matrix according to UT
C
      CALL TranBM(NAtoms,intcor,NDEG,UT,B,BP)
C
C  transpose the transformed B-Matrix
C
      DO 10 I=1,NAT3
      DO 10 J=1,NDEG
      BT(J,I) = BP(I,J)
 10   CONTINUE
C
C  Now transform the Hessian
C  Transform rows
C
      DO 20 I=1,NAT3
      DO 20 J=1,NDEG
      BMT(J,I) = SProd(NDEG,BT(1,I),HINT(1,J))
 20   CONTINUE
C
C  Transform columns
C
      DO 30 I=1,NAT3
      DO 30 J=1,I
      HESS(I,J) = SProd(NDEG,BMT(1,I),BT(1,J))
      HESS(J,I) = HESS(I,J)
 30   CONTINUE
c
      IF(IDB.EQ.1) THEN
c
       If(IPRNT.GT.1) WRITE(IOut,1200)
cc
       DO 60 I=1,NAT3
C
C  evaluate the GINT*dB/cart term and incorporate
C  into HESS - this will be done numerically
C
C  forward step
C
       XC(I) = XC(I) + delta
c
       CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $              coeff,  Jnk,    0,      1,      0,
     $              BSkal,  0,      BT,     Jnk,    IErr)
c
       If(IErr.NE.0) GO TO 95
C
C  transform the original B-Matrix to the new coordinates
C
       CALL TranBM(NAtoms,intcor,NDEG,UT,BT,BP)
C
C  backward step
C
       XC(I) = XC(I) - delta - delta
c
       CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $              coeff,  Jnk,    0,      1,      0,
     $              BSkal,  0,      BT,     Jnk,    IErr)
c
       If(IErr.NE.0) GO TO 95
C
C  transform the original B-Matrix to the new coordinates
C
       CALL TranBM(NAtoms,intcor,NDEG,UT,BT,BM)
C
C  restore coordinate
C
       XC(I) = XC(I) + delta
C
C  form the derivative contribution
C
       DO 50 J=1,NAT3
       gdBdx = Zero
       DO 40 IC=1,NDEG
       gdBdx = gdBdx + GINT(IC)*(BP(J,IC) - BM(J,IC))
 40    CONTINUE
       gdBdx = gdBdx/(delta + delta)
C
C  add contribution to cartesian Hessian prior
C  to transformation
C
       HESS(I,J) = HESS(I,J) - gdBdx
 50    CONTINUE
c
 60    CONTINUE
cc
      ENDIF
C
C  replace small diagonal elements by unity
C  (if molecule has symmetry, may have zero diagonals)
C
      DO 70 I=1,NAT3
      If(Abs(HESS(I,I)).LT.TolZero) HESS(I,I) = One
 70   CONTINUE
c
      RETURN
C
C -------------------------------------------------------
C ** ERROR SECTION **
C
 95   CONTINUE
      WRITE(IOut,1300)
      CALL OptExit(9)
c
 1000 FORMAT(/,' Transforming Internal Coordinate Hessian to',
     $         ' Cartesian Coordinates')
 1100 FORMAT(' Hessian Transformation does not Include',
     $       ' Derivative of B-matrix')
 1200 FORMAT(' Hessian Transformation Includes Derivative of',
     $       ' B-matrix')
 1300 FORMAT(/,2X,'***ERROR*** Something Wrong in Transforming',
     $            ' Hessian to Cartesian Coordinates')
c
      END
c ======================================================================
c
      SUBROUTINE HSSINT(NAtoms, intcor, NDEG,   IDB,    IPRNT,
     $                  XC,     GINT,   HESS,   ktyp,   klist,
     $                  Coeff,  SavTOR, BSkal,  UT,     BINV,
     $                  BT,     BP,     BM,     HINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transformation of cartesian Hessian to Hessian over the set
C  of internal coordinates according to:
C     HINT = B**(-1)t * (HESS - GINT*dB/dcart) * B**(-1)
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  intcor  -  number of primitive internals
C  NDEG    -  number of degrees of freedom
C             (and number of active internal coordinates)
C  IDB     -  integer flag controlling dB/dcart term
C               0 - do not include this term in Hessian transformation
C                   (for approximate Hessian & retention of Cartesian
C                    Hessian eigenvalue structure)
C               1 - include it
C  IPRNT   -  flag controlling printout
C  XC      -  cartesian coordinates
C  GINT    -  gradient in internal coordinates
C  HESS    -  cartesian Hessian
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C  klist   -  list of atoms involved in each primitive internal
C  Coeff   -  weighting of primitive in non-redundant
C             delocalized internal coordinate space
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  BSkal   -  scale factor for inverse-distance coordinates
C  UT      -  transformation matrix
C               i.e. which linear combination of primitive internals
C                    make up each compound (natural) internal coordinate
C  BINV    -  current inverse B-Matrix
C  BT      -  will hold B-Matrix in primitive internals
C  BP      -  will hold transformed B-Matrix forward step
C  BM      -  will hold transformed B-Matrix backward step
C  HINT    -  on exit Hessian in internal coordinates
C             (MAY USE SAME STORAGE AS BM - see calling routine)
C
C
      REAL*8 XC(3*NATOMS),GINT(NDEG),HESS(3*NATOMS,3*NATOMS),
     $       UT(intcor,NDEG),BT(3*NATOMS,intcor),
     $       BINV(3*NATOMS,NDEG),BP(3*NATOMS,NDEG),
     $       BM(3*NATOMS,NDEG),HINT(NDEG,NDEG)
C
      PARAMETER (Zero=0.0d0,delta=0.0005d0)
C
C
      IOut = ioutfil('iout')
c
      If(IPRNT.GT.1) Then
       WRITE(IOut,1000)
       If(IDB.EQ.0) WRITE(IOut,1100)
      EndIf
c
      NAT3 = 3*NAtoms
c
      IF(IDB.EQ.1) THEN
c
       If(IPRNT.GT.1) WRITE(IOut,1200)
cc
       DO 40 I=1,NAT3
C
C  evaluate the GINT*dB/cart term and incorporate
C  into HESS - this will be done numerically
C
C  forward step
C
       XC(I) = XC(I) + delta
c
       CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $              coeff,  Jnk,    0,      1,      0,
     $              BSkal,  0,      BT,     Jnk,    IErr)
c
       If(IErr.NE.0) GO TO 95
C
C  transform the original B-Matrix to the new coordinates
C
       CALL TranBM(NAtoms,intcor,NDEG,UT,BT,BP)
C
C  backward step
C
       XC(I) = XC(I) - delta - delta
c
       CALL BMATRIX(NAtoms, XC,     intcor, ktyp,   klist,
     $              coeff,  Jnk,    0,      1,      0,
     $              BSkal,  0,      BT,     Jnk,    IErr)
c
       If(IErr.NE.0) GO TO 95
C
C  transform the original B-Matrix to the new coordinates
C
       CALL TranBM(NAtoms,intcor,NDEG,UT,BT,BM)
C
C  restore coordinate
C
       XC(I) = XC(I) + delta
C
C  form the derivative contribution
C
       DO 30 J=1,NAT3
       gdBdx = Zero
       DO 20 IC=1,NDEG
       gdBdx = gdBdx + GINT(IC)*(BP(J,IC) - BM(J,IC))
 20    CONTINUE
       gdBdx = gdBdx/(delta + delta)
C
C  add contribution to cartesian Hessian prior
C  to transformation
C
       HESS(I,J) = HESS(I,J) - gdBdx
 30    CONTINUE
c
 40    CONTINUE
cc
      ENDIF
C
C  Now do the transformation
C  Transform columns
C
      DO 50 I=1,NAT3
      DO 50 J=1,NDEG
      BP(I,J) = SProd(NAT3,HESS(1,I),BINV(1,J))
 50   CONTINUE
C
C  Transform rows
C
      DO 60 I=1,NDEG
      DO 60 J=1,I
      HINT(I,J) = SProd(NAT3,BINV(1,I),BP(1,J))
      HINT(J,I) = HINT(I,J)
 60   CONTINUE
c
      RETURN
c
 95   CONTINUE
      WRITE(IOut,1300)
      CALL OptExit(9)
c
 1000 FORMAT(/,' Transforming Cartesian Hessian to Internal',
     $         ' Coordinates')
 1100 FORMAT(' Hessian Transformation does not Include',
     $       ' Derivative of B-matrix')
 1200 FORMAT(' Hessian Transformation Includes Derivative of',
     $       ' B-matrix')
 1300 FORMAT(/,2X,'***ERROR*** Something Wrong in Transforming',
     $            ' Hessian to Internal Coordinates')
c
      END
c ======================================================================
c
      SUBROUTINE HssINTZ(NZ,     NAtoms, NIC,    NVar,   IPRNT,
     $                   ZSymb,  GEO,    IGEO,   IG,     XZ,
     $                   GINT,   IDB,    HESS,   HZ,     ITMP,
     $                   BINV,   BP,     BM,     HTMP,   HINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transformation of cartesian Hessian to Hessian over the set
C  of Z-Matrix internal coordinates according to:
C     HINT = B**(-1) * (HESS - GINT*dB/dcart) * B**(-1)t
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  NAtoms  -  number of real atoms
C  NIC     -  total number of internal coordinates (3*NZ-6)
C  NVar    -  number of variables in Z-Matrix
C             i.e. parameters to be optimized
C  IPRNT   -  print flag
C  ZSymb   -  Z-Matrix symbols (used to sort dummy atoms)
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO    -  Z-matrix connectivity
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  XZ      -  Cartesian coordinates of all centres in
C             Z-matrix orientation
C  GINT    -  current (fully expanded) internal coordinate gradient
C  IDB     -  integer flag controlling dB/dcart term
C               0 - do not include this term in Hessian transformation
C                   (for approximate Hessian & retention of Cartesian
C                    Hessian eigenvalue structure)
C               1 - include it
C  HESS    -  Cartesian Hessian
C  HZ      -  scratch storage for expanded Cartesian Hessian
C  ITMP    -  integer scratch vector (length NIC)
C  BINV    -  current inverse B-Matrix
C  BP      -  will hold B-Matrix forward step
C  BM      -  will hold B-Matrix backward step
C  HTMP    -  fully expanded internal coordinate Hessian
C             (MAY USE SAME STORAGE AS BM - see calling routine)
C  HINT    -  on exit contains Z-Matrix Hessian
C
C
      DIMENSION GEO(NZ,3),IGEO(NZ,3),IG(NIC),XZ(3*NZ),ITMP(NIC)
      REAL*8 GINT(NIC),HESS(3*NAtoms,3*NAtoms),HZ(3*NZ,3*NZ),
     $       BINV(3*NZ,3*NZ),BP(3*NZ,3*NZ),BM(3*NZ,3*NZ),
     $       HTMP(NIC,NIC),HINT(NVar,NVar)
      CHARACTER*8 ZSymb(NZ)
C
      PARAMETER (Zero=0.0d0,delta=0.0005d0)
C
C
      IOut = ioutfil('iout')
c
      If(IPRNT.GT.1) Then
       WRITE(IOut,1000)
       If(IDB.EQ.0) WRITE(IOut,1100)
      EndIf
c
      NZ3 = 3*NZ
C
C  expand HESS into HZ (taking account of dummy atoms)
C
      IF(NAtoms.EQ.NZ) THEN
cc
       CALL CpyVEC(NZ3*NZ3,HESS,HZ)
cc
      ELSE
cc
C  expand out each row/column
C
       IAtm = 0
       DO 20 I=1,NZ
       II = 3*(I-1)
       If(ZSymb(I)(1:2).NE.'Du') Then
        IAtm = IAtm+1
        IK = 3*(IAtm-1)
       EndIf
c
       JAtm = 0
       DO 10 J=1,NZ
       JJ = 3*(J-1)
       If(ZSymb(J)(1:2).NE.'Du') Then
        JAtm = JAtm+1
        JL = 3*(JAtm-1)
       EndIf
c
       IF(ZSymb(I)(1:2).EQ.'Du'.OR.ZSymb(J)(1:2).EQ.'Du') THEN
        HZ(II+1,JJ+1) = Zero
        HZ(II+2,JJ+1) = Zero
        HZ(II+3,JJ+1) = Zero
        HZ(II+1,JJ+2) = Zero
        HZ(II+2,JJ+2) = Zero
        HZ(II+3,JJ+2) = Zero
        HZ(II+1,JJ+3) = Zero
        HZ(II+2,JJ+3) = Zero
        HZ(II+3,JJ+3) = Zero
       ELSE
        HZ(II+1,JJ+1) = HESS(IK+1,JL+1)
        HZ(II+2,JJ+1) = HESS(IK+2,JL+1)
        HZ(II+3,JJ+1) = HESS(IK+3,JL+1)
        HZ(II+1,JJ+2) = HESS(IK+1,JL+2)
        HZ(II+2,JJ+2) = HESS(IK+2,JL+2)
        HZ(II+3,JJ+2) = HESS(IK+3,JL+2)
        HZ(II+1,JJ+3) = HESS(IK+1,JL+3)
        HZ(II+2,JJ+3) = HESS(IK+2,JL+3)
        HZ(II+3,JJ+3) = HESS(IK+3,JL+3)
       ENDIF
c
 10    CONTINUE
 20    CONTINUE
cc
      ENDIF
c
      If(IPRNT.GT.6) Then
       WRITE(IOut,1300)
       CALL PrntMAT(NZ3,NZ3,NZ3,HZ)
      EndIf
C
C  now evaluate the GINT*dB/cart term and incorporate
C  into HZ - this will be done numerically
C
      IF(IDB.EQ.1) THEN
c
       If(IPRNT.GT.1) WRITE(IOut,1200)
cc
       DO 50 I=1,NZ3
C
C  forward step
C
       XZ(I) = XZ(I) + delta
c
       CALL BMatZ(NZ,XZ,IGEO,IPRNT-1,BP)
C
C  backward step
C
       XZ(I) = XZ(I) - delta - delta
c
       CALL BMatZ(NZ,XZ,IGEO,IPRNT-1,BM)
C
C  restore coordinate
C
       XZ(I) = XZ(I) + delta
C
C  form the derivative contribution
C
       DO 40 J=1,NZ3
       gdBdx = Zero
       DO 30 IC=1,NIC
       gdBdx = gdBdx + GINT(IC)*(BP(J,IC) - BM(J,IC))
 30    CONTINUE
       gdBdx = gdBdx/(delta + delta)
C
C  add contribution to expanded cartesian Hessian prior
C  to transformation
C
       HZ(I,J) = HZ(I,J) - gdBdx
 40    CONTINUE
c
 50    CONTINUE
c
       If(IPRNT.GT.6) Then
        WRITE(IOut,1400)
        CALL PrntMAT(NZ3,NZ3,NZ3,HZ)
       EndIf
cc
      ENDIF
C
C  Now do the transformation
C  Transform columns
C
      DO 60 I=1,NIC
      DO 60 J=1,NZ3
      Val = Zero
      DO 59 K=1,NZ3
      Val = Val + BINV(I,K)*HZ(K,J)
 59   CONTINUE
      BP(I,J) = Val
 60   CONTINUE
C
C  Transform rows
C
      DO 70 I=1,NIC
      DO 70 J=1,I
      Val = Zero
      DO 69 K=1,NZ3
      Val = Val + BP(I,K)*BINV(J,K)
 69   CONTINUE
      HTMP(I,J) = Val
      HTMP(J,I) = Val
 70   CONTINUE
c
      If(IPRNT.GT.5) Then
       WRITE(IOut,1500)
       CALL PrntMAT(NIC,NIC,NIC,HTMP)
      EndIf
C  ....................................................................
C
C  Now transform from internal to Z-Matrix variables
C
C  form ITMP vector from IG
C
      IPar = 0
      DO 80 I=1,NIC
      IF(IG(I).EQ.1000) THEN
       ITMP(I) = 0
      ELSE IF(IG(I).EQ.0) THEN
       IPar = IPar+1
       ITMP(I) = IPar
      ELSE
       ITMP(I) = ITMP(IAbs(IG(I)))
      ENDIF
 80   CONTINUE
C
C  first transform the columns
C
      CALL ZeroIT(BP,NZ3*NZ3)
      DO 90 I=1,NIC
      DO 89 J=1,NIC
      If(ITMP(J).EQ.0) GO TO 89
      IBL = ITMP(J)
      FX = HTMP(J,I)
      If(IG(J).LT.0) FX=-FX
      BP(IBL,I) = BP(IBL,I) + FX
 89   CONTINUE
 90   CONTINUE
C
C  now transform the rows
C
      CALL ZeroIT(HINT,NVar*NVar)
      DO 100 I=1,NVar
      DO 99  J=1,NIC
      If(ITMP(J).EQ.0) GO TO 99
      IBL = ITMP(J)
      FX = BP(I,J)
      If(IG(J).LT.0) FX=-FX
      HINT(I,IBL) = HINT(I,IBL) + FX
 99   CONTINUE
 100  CONTINUE
C
      RETURN
c
 1000 FORMAT(/,' Transforming Cartesian Hessian to Z-Matrix',
     $         ' Coordinates')
 1100 FORMAT(' Hessian Transformation does not Include',
     $       ' Derivative of B-matrix')
 1200 FORMAT(' Hessian Transformation Includes Derivative of',
     $       ' B-matrix')
 1300 FORMAT(/,'   Expanded Cartesian Hessian Matrix')
 1400 FORMAT(/,'   Cartesian Hessian including derivative',
     $           ' B-Matrix Term')
 1500 FORMAT(/,'   Full Hessian in internal coordinates')
c
      END
c ======================================================================
c
      SUBROUTINE IntCnvCON(NS,     NCons,  NCon,   ICHK,   XINT,
     $                     D,      GC,     EC,     EOld,   IPRNT,
     $                     TolG,   TolD,   TolE,   NEG,    NegReq,
     $                     Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for convergence during constrained optimization in
C  delocalized internal coordinates
C
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
C  NS      -  number of unconstrained degrees of freedom
C  NCons   -  total number of constraints
C  NCon    -  number of unsatisfied constraints
C  ICHK    -  array indicating which constrains are satisfied
C               ICHK(I) = 0 - ith constraint is satisfied
C               ICHK(I) = 1 - ith constraint is not satisfied
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
      REAL*8 XINT(NS+2*NCons),D(NS+2*NCon),GC(NS+2*NCon)
      DIMENSION ICHK(NCons)
      LOGICAL Cnvgd
C
      CHARACTER*3 DCnvgd,GCnvgd,ECnvgd
C
C
      IOut = ioutfil('iout')
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
      NDEG = NS+NCons
      GMax = Abs(GC(1))
      DMax = Abs(D(1))
c
      DO 10 I=1,NS
      GX = GC(I)
      DX = D(I)
      If(Abs(GX).GT.GMax) GMax = Abs(GX)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
      If(IPRNT.GT.2) Then
       CX = XINT(I)
       CN = CX + DX
       WRITE(IOut,1100) I,CX,GX,DX,CN
      EndIf
c
 10   CONTINUE
c
      IT = NS
      DO 20 I=1,NCons
      IF(ICHK(I).EQ.1) THEN
       IT = IT+1
       GX = GC(IT)
       DX = D(IT)
       If(Abs(GX).GT.GMax) GMax = Abs(GX)
       If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
       If(IPRNT.GT.2) Then
        CX = XINT(NS+I)
        CN = CX + DX
        WRITE(IOut,1100) IT,CX,GX,DX,CN
       EndIf
c
      ENDIF
 20   CONTINUE
c
      If(IPRNT.GT.2) WRITE(IOut,1700)
c
      IT = NS+NCon
      DO 30 I=1,NCons
      IF(ICHK(I).EQ.1) THEN
       IT = IT+1
       GX = GC(IT)
       DX = D(IT)
       If(Abs(GX).GT.GMax) GMax = Abs(GX)
       If(Abs(DX).GT.DMax) DMax = Abs(DX)
c
       If(IPRNT.GT.2) Then
        CX = XINT(NDEG+I)
        CN = CX + DX
        WRITE(IOut,1100) I,CX,GX,DX,CN
       EndIf
c
      ENDIF
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
 1000 FORMAT(//,5X,'Parameter Values and Displacements in Internal',
     $              ' Coordinates',/,
     $       '  Coordinate  Current Value    Gradient  ',
     $       'Displacement   New Value')
 1100 FORMAT(5X,I3,9X,F10.6,3X,F10.6,3X,F10.6,3X,F10.6)
 1300 FORMAT(/,29X,'Maximum     Tolerance    Cnvgd?')
 1400 FORMAT(9X,'Gradient           ',F8.6,6X,F8.6,5X,A3)
 1500 FORMAT(9X,'Displacement       ',F8.6,6X,F8.6,5X,A3)
 1600 FORMAT(9X,'Energy change     ',F9.6,6X,F8.6,5X,A3,/)
 1700 FORMAT(/,15X,' Lagrange Multipliers for Constraints',/,
     $       '  Constraint  Current Value    Gradient  ',
     $       'Displacement   New Value')
c
      END
c ======================================================================
c
      SUBROUTINE IntConGRAD(NS,     NCons,  NCon,   ICHK,   RCHK,
     $                      XINT,   GC,     GINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine modifies the internal coordinate gradient
C  to incorporate geometrical constraints using the method
C  of Lagrange multipliers
C
C  ARGUMENTS
C
C  NS      -  number of unconstrained degrees of freedom
C  NCons   -  total number of constraints
C  NCon    -  number of unsatisfied constraints
C  ICHK    -  array indicating which constrains are satisfied
C               ICHK(I) = 0 - ith constraint is satisfied
C               ICHK(I) = 1 - ith constraint is not satisfied
C  RCHK    -  values for constraint differences
C             (current value minus desired value)
C  XINT    -  coordinate vector (including Lagrange multipliers)
C  GC      -  gradient over all degrees of freedom
C  GINT    -  gradient vector (including dL/dLambda on exit)
C
C
      DIMENSION ICHK(NCons),RCHK(NCons),XINT(*),GC(*),GINT(*)
C
C
C  The first NS (the active) coordinates are always kept
C  From the NCons total constraints, only the NCon unsatisfied
C  constraints are retained
C
      CALL CpyVEC(NS,GC,GINT)
c
      IT = 0
      DO 10 I=1,NCons
      IF(ICHK(I).EQ.1) THEN
       IT = IT+1
C
C  modify the gradient to incorporate Lambda
C
       GINT(NS+IT) = GC(NS+I) + XINT(NS+NCons+I)
C
C  now append the dL/dLambda term
C
       GINT(NS+NCon+IT) = -RCHK(I)
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE IntLgngeHESS(NS,NCons,NCon,ICHK,HESS,HESSC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine forms the Lagrangian Hessian from the standard
C  internal coordinate Hessian and the (unit) constraint normals
C
C  ARGUMENTS
C
C  NS      -  number of unconstrained degrees of freedom
C  NCons   -  total number of constraints
C  NCon    -  number of unsatisfied constraints
C  ICHK    -  array indicating which constrains are satisfied
C               ICHK(I) = 0 - ith constraint is satisfied
C               ICHK(I) = 1 - ith constraint is not satisfied
C  HESS    -  Hessian over all degrees of freedom
C  HESSC   -  on exit contains Lagrangian Hessian
C
C
      DIMENSION ICHK(NCons),HESS(NS+NCons,NS+NCons),
     $          HESSC(NS+2*NCon,NS+2*NCon)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
C  The first NS (the active) coordinates are always kept
C  From the NCons total constraints, only the NCon unsatisfied
C  constraints are retained
C
C  first copy the active part of HESS into HESSC
C
      DO 10 J=1,NS
      DO 10 I=1,NS
      HESSC(I,J) = HESS(I,J)
 10   CONTINUE
C
C  zero out the additional part
C
      DO 20 J=NS+NCon+1,NS+2*NCon
      DO 20 I=1,NS+2*NCon
      HESSC(I,J) = Zero
      HESSC(J,I) = Zero
 20   CONTINUE
C
C  now copy in the unsatisfied constraints
C  also the constraint (unit) normals
C
      IT = 0
      DO 50 IC=1,NCons
      IF(ICHK(IC).EQ.1) THEN
       IT = IT+1
       DO 30 I=1,NS
       HESSC(I,NS+IT) = HESS(I,NS+IC)
       HESSC(NS+IT,I) = HESS(I,NS+IC)
 30    CONTINUE
       JT = 0
       DO 40 JC=1,NCons
       IF(ICHK(JC).EQ.1) THEN
        JT = JT+1
        HESSC(NS+IT,NS+JT) = HESS(NS+IC,NS+JC)
       ENDIF
 40    CONTINUE
C
C  add the constraint normals
C
       HESSC(NS+NCon+IT,NS+IT) = One
       HESSC(NS+IT,NS+NCon+IT) = One
      ENDIF
 50   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE LgngeHESS(NCTR3,NCons,HESS,GCC,HESSC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine forms the Lagrangian Hessian from the modified
C  Cartesian Hessian and the constraint normals
C
C  ARGUMENTS
C
C  NCTR3   -  3*number of centres (for dimensioning)
C  NCons   -  number of constraints
C  HESS    -  Cartesian Hessian
C  GCC     -  matrix of constraint normals
C  HESSC   -  on exit contains Lagrangian Hessian
C
C
      REAL*8 HESS(NCTR3,NCTR3),GCC(NCTR3,NCons),
     $       HESSC(NCTR3+NCons,NCTR3+NCons)
C
C
C  first copy HESS into HESSC
C
      DO 10 J=1,NCTR3
      DO 10 I=1,NCTR3
      HESSC(I,J) = HESS(I,J)
 10   CONTINUE
C
C  now copy in the constraint normals
C  and zero the lower right block
C
      DO 30 IC=1,NCons
      IT = NCTR3 + IC
      DO 20 I=1,NCTR3
      HESSC(I,IT) = GCC(I,IC)
      HESSC(IT,I) = GCC(I,IC)
 20   CONTINUE
      CALL ZeroIT(HESSC(NCTR3+1,IT),NCons)
 30   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE LincGRAD(NAtoms,I,J,K,L,XC,linc,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the linear coplanar bend I-J-K-L
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  first atom in linear bend
C  J       -  central atom in linear bend
C  K       -  third atom in linear bend
C  L       -  non-linear atom (attached to K)
C  XC      -  Cartesian coordinates
C  linc    -  on exit contains linear coplanar bend
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate linc gradient   NOT IMPLEMENTED
C              .false. -  skip gradient calculation
C  G       -  on exit contains linear coplanar bend gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms),linc
      LOGICAL grd
C
      DIMENSION R1(3),R2(3),R3(3)
C
      PARAMETER (PI=3.14159 26535 89793d0)
C
C
      CALL VecDIF(R1,rn,XC(1,I),XC(1,J))
      CALL VecDIF(R2,rn,XC(1,L),XC(1,J))
      CALL VecDIF(R3,rn,XC(1,K),XC(1,J))
      CO = SProd(3,R2,R1)
      CP = SProd(3,R3,R2)
c
      linc = PI - ACOS(CO) - ACOS(CP)
C
C
      If(.NOT.grd) RETURN
C
C
C  set terms for linear coplanar bend derivatives
C  ** NOT IMPLEMENTED **
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE LinpGRAD(NAtoms,I,J,K,L,XC,linp,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the linear coplanar bend I-J-K-L
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  first atom in linear bend
C  J       -  central atom in linear bend
C  K       -  third atom in linear bend
C  L       -  non-linear atom (attached to K)
C  XC      -  Cartesian coordinates
C  linp    -  on exit contains linear perpendicular bend
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate linp gradient   NOT IMPLEMENTED
C              .false. -  skip gradient calculation
C  G       -  on exit contains linear perp. bend gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms),linp
      LOGICAL grd
C
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0,small=1.0d-6)
      PARAMETER (PI=3.14159 26535 89793d0)
C
C
      linp = Zero
c
      CALL VecDIF(R1,rn,XC(1,I),XC(1,J))
      CALL VecDIF(R2,rn,XC(1,L),XC(1,J))
      CALL VecDIF(R3,rn,XC(1,K),XC(1,J))
      CALL Cross(R2,R1,R4)
      WNorm = SQRT(SProd(3,R4,R4))
c
      If(WNorm.LT.small) RETURN
c
      WNorm = One/WNorm
      CALL Vscal(3,WNorm,R4)
      CALL Cross(R3,R2,R5)
      XNorm = SQRT(SProd(3,R5,R5))
c
      If(XNorm.LT.small) RETURN
c
      XNorm = One/XNorm
      CALL Vscal(3,XNorm,R5)
      CO = SProd(3,R1,R5)
      CP = SProd(3,R3,R4)
c
      linp = Half*(PI - ACOS(CO) - ACOS(CP))
C
C
      If(.NOT.grd) RETURN
C
C
C  set terms for linear coplanar bend derivatives
C  ** NOT IMPLEMENTED **
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE LinTOPO(NAtoms, IC,     IVdW,   intcor, ktyp,
     $                   klist,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms internal coordinates for a Linear molecule
C  These are taken as simple bond stretches
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IC      -  atomic connectivity
C  IVdW    -  flag for VdW cluster
C              0 - standard optimization
C              1 - surface adsorption/reaction
C              2 - VdW cluster  (use 1/R**6 distance coordinates)
C
C  on exit  ....
C
C  intcor  -  total number of internal coordinates
C  ktyp    -  integer array of internal coordinate types
C             (all entries set to 1 for bond stretch)
C  klist   -  list of atoms involved in each primitive
C  IErr    -  error flag   0 - success
C                         -1 - something wrong with IC matrix
C
      DIMENSION IC(NAtoms,NAtoms)
      DIMENSION ktyp(*),klist(4,*)
C
C
C  initialize
C
      IErr = -1
      intcor = 0
c
      DO 10 I=2,NAtoms
      DO 10 J=1,I-1
      IF(IC(I,J).GT.0) THEN
C
C  atoms I and J are bonded
C
       intcor = intcor + 1
       ktyp(intcor) = 1
       If(IVdW.EQ.2) ktyp(intcor) = 7
       klist(1,intcor) = I
       klist(2,intcor) = J
      ENDIF
 10   CONTINUE
C
C  check that we have the correct number of internals
C
      IF(intcor.NE.NAtoms-1) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       RETURN
      ENDIF
C
      IErr = 0
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Something Wrong with Atom Connectivity',
     $            ' in <LinTOPO>')
c
      END
c ======================================================================
c
      SUBROUTINE MakeDIC(NAtoms, intcor, NVib,   IPRNT,  ktyp,
     $                   B,      INDX,   G,      EigVal, V,
     $                   UT,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Constructs the lower triangle of the G = B*B(t) matrix and extracts
C  the highest NVib eigenvalues from it. The corresponding eigenvectors
C  are taken as the definition of our delocalized internal coordinates
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  intcor  -  total number of primitive internal coordinates
C  NVib    -  number of vibrational degrees of freedom
C             (should be same as number of non-redundant coordinates)
C  IPRNT   -  print flag
C  ktyp    -  integer array containing internal coordinate type
C  B       -  non-zero elements of primitive B Matrix
C  INDX    -  index to non-zero elements of B Matrix
C  G       -  storage for G Matrix (lower triangle)
C  EigVal  -  storage for eigenvalues of G
C  V       -  scratch storage
C  UT      -  transformation matrix for non-redundant internals
C  IErr    -  error flag
C              0 - non-redundant set of natural internals generated
C             -1 - something went wrong
C
C
      DIMENSION B(12,intcor),INDX(12,intcor),EigVal(intcor),
     $          G(intcor*(intcor+1)/2),V(intcor),
     $          UT(intcor,NVib),ktyp(intcor)
      Dimension kval(7)
      integer*4 i4err
C
      Data kval/ 6, 9, 12, 12, 12, 12, 6 /
C
      PARAMETER (Zero=0.0d0,thrsh=1.0d-7)
C
C
      IOut = ioutfil('iout')
      NAT3 = 3*NAtoms
c
      IErr = -1
c
      If(IPRNT.GT.3) WRITE(IOut,1000)
C
C  construct G (using sparse matrix multiply)
C  (Don't like this, has IF in inner loop, but it is usually
C   a "once only" and can't think of anything obviously better)
C
      IJ = 0
      DO 30 I=1,intcor
      IVal = kval(ktyp(I))
      DO 20 J=1,I
      IJ = IJ+1
      JVal = kval(ktyp(J))
      Val = Zero
      DO 11 K=1,IVal
      KK = INDX(K,I)
      DO 10 L=1,JVal
      LL = INDX(L,J)
      If(LL.EQ.KK) Val = Val + B(K,I)*B(L,J)
 10   CONTINUE
 11   CONTINUE
      G(IJ) = Val
 20   CONTINUE
 30   CONTINUE
C
C  "diagonalize" it
C  we only need the highest NVib vectors
C
      IRed = intcor-NVib
      istrt = IRed+1
      if(istrt.gt.intcor)istrt=intcor
c
      CALL DSPEVX('V',    'I',    'U',    intcor, G,
     $            jnk,    jnk,    istrt,  intcor, 1.0d-10,
     $            Nfound, EigVal, UT,     intcor, B,
     $            INDX,   V,      i4err)
      ierr=int(i4err)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1200)
       RETURN
      ENDIF
C
C
C  see what we've got?
C
      If(IPRNT.GT.5) Then
       WRITE(IOut,1400)
       WRITE(IOut,1500) (EigVal(I),I=1,NVib)
       WRITE(IOut,1600)
       CALL PrntMAT(NVib,intcor,NVib,UT)
      EndIf
C
C  check the lowest eigenvalue
C
      IF(Eigval(1).LT.thrsh) THEN
       WRITE(IOut,1300)
       IErr = -1
       RETURN
      ENDIF
c
      If(IPRNT.GT.2) WRITE(IOut,1700) IRed
c
      RETURN
c
 1000 FORMAT(/,' Generating Delocalized Internal Coordinates')
 1200 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize B*B(t) in',
     $            ' <MakeDIC>')
 1300 FORMAT(/,2X,'***ERROR*** B*B(t) Has Small or Zero Eigenvalues',
     $            ' in <MakeDIC>',/,
     $        '     Check the Atomic Connectivity')
 1400 FORMAT(/,'   Non-Zero Eigenvalues of B*B(t):')
 1500 FORMAT(1X,6F12.6)
 1600 FORMAT(/,'  Preliminary Set of Delocalized Internal Coordinates:')
 1700 FORMAT(/,' There are ',I6,' Redundant primitive internals')
c
      END
c ======================================================================
c
      SUBROUTINE MakeINTC(NAtoms, AtSymb, XC,     NMol,   IMOL,
     $                    GROUP,  NDEG,   NCons,  ICTYP,  RCON,
     $                    ICON,   ICNUM,  NPComp, IPComp, IPCNUM,
     $                    NCon,   NFix,   IFIX,   NDrive, IDRTYP,
     $                    IDRIVE, FDRIVE, LDRIVE, NQ,     NTrans,
     $                    TRANS,  NEqATM, IPRNT,  NIC,    NPrim,
     $                    NPrim0, ICNNCT, IGen,   IType,  ITors,
     $                    CutOff, BSkal,  PThrsh, ktyp,   klist,
     $                    BPRIM,  INDB,   Coeff,  SavTOR, NCmp,
     $                    NP1,    INT1,   UT,     XPrim,  IGEO,
     $                    IMAP,   IG,     IOrd,   NMem,   Z,
     $                    IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ---------------------------------------------------------------
C    I N T E R N A L   C O O R D I N A T E   G E N E R A T I O N
C  ---------------------------------------------------------------
C
C  This routine is responsible for generating internal coordinates
C  directly from input Cartesians. First of all, it generates a set
C  of redundant "primitive" internals based on atomic connectivity
C  and from these primitives generates the internal coordinates.
C
C  Primitive internals can be:
C      1.   Stretches
C      2.   Planar Bends
C      3.   Out-of-Plane Bends
C      4.   Proper Torsions
C      5.   Linear Coplanar Bend
C      6.   Linear Perpendicular Bend
C      7.   Inverse Stretches
C
C  These primitives are used to generate a set of either (i) natural
C  or (ii) delocalized internal coordinates. Natural internals are
C  obtained from an exhaustive topological analysis using essentially
C  the original code of Pulay and Fogarasi; delocalized internals are
C  obtained by constructing and diagonalizing the G [B*B(t)] matrix.
C
C  For normal molecules typically only stretches, bends and torsions
C  are used (out-of-plane bends are included for natural internals).
C  Coplanar and perpendicular bends are used for linear or near-linear
C  arrangements of three or more atoms. Inverse stretch coordinates
C  can be used for clusters (typically VdW clusters).
C
C
C  References
C  ----------
C
C  "Systematic ab initio Gradient Calculation of Molecular Geometries,
C   Force Constants and Dipole Moment Derivatives"
C   P.Pulay, G.Fogarasi, F.Pang and J.E.Boggs,  J.Am.Chem.Soc. 101 (1979) 2550
C
C  "The Calculation of ab initio Molecular Geometries: Efficient Optimization
C   by Natural Internal Coordinates and Empirical Corrections by Offset Forces"
C   G.Fogarasi, X.Zhou, P.W.Taylor and P.Pulay,  J.Am.Chem.Soc. 114 (1992) 8191
C
C  "Geometry Optimization in Redundant Internal Coordinates"
C   P.Pulay and G.Fogarasi  J.Chem.Phys. 96 (1992) 2856
C
C  "The Generation and Use of Delocalized Internal Coordinates
C   in Geometry Optimization"
C   J.Baker, A.Kessi and B.Delley  J.Chem.Phys. 105 (1996) 192
C
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates (in cms and symmetry alligned)
C  NMol    -  number of molecules (for, e.g., cluster optimizations)
C  IMOL    -  pointers to start/end of molecules in XC array
C  GROUP   -  molecular point group
C  NDEG    -  number of internal degrees of freedom
C  NCons   -  number of constraints (fixed primitives)
C  ICTYP   -  list of constraint types
C  RCON    -  constraint values (in atomic units)
C  ICON    -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       angle constraint
C               IC1-IC2-IC3-IC4   dihedral constraint
C  ICNUM   -  on exit list of which primitives are to be fixed
C  NPComp  -  number of primitives in composite constraints
C  IPComp  -  constraint type and atoms involved in constraint
C  IPCNUM  -  on exit indicates which primitive in list corresponds
C             to the desired composite constraints
C  NCon    -  on exit contains actual number of constraints
C  NFix    -  on entry simply indicates there are fixed surface atoms
C             on exit number of fixed atoms
C  IFIX    -  list of fixed surface atoms (for surface optimization)
C  NDrive  -  number of primitives to drive
C  IDRTYP  -  type of each primitive to drive
C               1 - distance               stre
C               2 - bond angle             bend
C               3 - out-of-plane bend      outp
C               4 - dihedral angle         tors
C  IDRIVE  -  definition of each primitive (i.e., list of atoms)
C  FDRIVE  -  force (in internal coordinates) to be applied
C  LDRIVE  -  on exit primitive number, i.e., which primitive we are
C             referring to in list of all primitives
C  NQ      -  number of symmetry unique atoms
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  IPRNT   -  flag controlling print out
C  NIC     -  maximum allowed number of primitive internals
C  NPrim   -  on entry size of primitive space on previous cycle
C             on exit size of primitive space on current cycle
C  NPrim0  -  number of intramolecular primitives (when IType=1)
C  ICNNCT  -  atomic connectivity matrix
C  IGen    -  flag controlling generation of internal coordinates
C               0 - use existing set of internal coordinates
C                   (this routine should not be called)
C               1 - generate a new set of non-redundant internals
C                   from an existing primitive space
C               2 - generate a new set of underlying primitives
C                   and a new set of delocalized internals
C              -1 - generate a set of natural internal coordinates
C  IType   -  flag for cluster/surface optimizations
C               0 - standard optimization
C               1 - adsorbate-surface optimization
C               2 - VdW cluster  (use 1/R distance coordinates)
C  ITors   -  use of torsions in delocalized internal coordinates
C               0 - use torsions throughout
C               1 - do not use torsions for FIRST molecule only
C             (in surface optimizations, the first "molecule" is the surface)
C  CutOff  -  bond distance threshold for inverse-power distance
C             coordinates
C  BSkal   -  scale factor for inverse-distance coordinates
C  PThrsh  -  smallest value (in radians) allowed for near-linear
C             angle constraint
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
C  BPRIM   -  non-zero elements of B-Matrix over primitive internals
C  INDB    -  index to non-zero elements of B-Matrix
C  Coeff   -  primitive weighting factors
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  NCmp    -  total number of non-zero natural internal components
C  NP1     -  number of primitives in each NIC
C  INT1    -  indices for all non-zero components per NIC
C  UT      -  non-zero NIC components OR full set of
C             delocalized internal coordinates
C  XPrim   -  values of primitive internals
C ................................................................
C  IGEO    -  Z-matrix connectivity for back-transformation
C  IMAP    -  mapping order from original to Z-matrix order
C  IG      -  indexing array of which primitives are in Z-matrix
C  IOrd    -  flag indicating change of order between original
C             atom order and order in Z-matrix
C              0 - no change
C              1 - change
C             -1 - Z-matrix mapping failed; need full back-transformation
C  ..................................................................
C  NMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  IErr    -  error flag    0 - full set of internals found
C                          -1 - unable to find full set
C
C
      DIMENSION XC(3,NAtoms),IMOL(NMol)
      DIMENSION ktyp(NIC),klist(4,NIC),Coeff(NIC),SavTOR(NIC),
     $          UT(*),XPrim(NIC),BPRIM(12*NIC),
     $          INDB(12*NIC),ICNNCT(NAtoms,NAtoms)
      INTEGER NP1(NIC),INT1(12*NIC)
      DIMENSION IFIX(3,NAtoms),IGEO(NAtoms,4),IMAP(NAtoms,2),IG(NDEG)
      DIMENSION ICTYP(NCons),ICON(4,NCons),RCON(NCons),ICNUM(NCons),
     $          IPComp(5,NPComp),IPCNUM(NPComp)
      DIMENSION IDRTYP(NDrive),IDRIVE(4,NDrive),FDRIVE(NDrive),
     $          LDRIVE(NDrive)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 GROUP
      LOGICAL Checked
C
      DIMENSION Z(NMem)
C
      PARAMETER (thrbnd=0.85d0)     !  distance ratio for bonding
      PARAMETER (One=1.0d0)
      PARAMETER (GWght=0.2d0)
C
C
      IOut = ioutfil('iout')
C
C  initialize
C
      NVib = 3*NAtoms-6
      IErr = -1
      If(NTrans.GT.1) IOrd = -1     ! no Z-matrix with symmetry for now
cc      Checked = .False.
      Checked = .True.
c
      IF(IGEN.GE.1) THEN
C
C
C =================================================================
C  D E L O C A L I Z E D   I N T E R N A L   C O O R D I N A T E S
C =================================================================
C
       IF(IGen.EQ.2) THEN
C
C  -----------------------------------
C  Generation of primitive internals
C  -----------------------------------
C
C  generate connectivity data from interatomic distances
C  get atomic numbers from atomic symbols
C
        CALL GetAtNo(NAtoms,AtSymb,Z)
C
C  determine connectivity matrix
C
        CALL CONNECTM(NAtoms,Z,XC,thrbnd,ICNNCT,IErr)
c
        If(IErr.NE.0) Then
         If(IPRNT.GT.2) WRITE(IOut,1100)
         RETURN
        EndIf
C
        IF(GROUP.EQ.'c*v '.OR.GROUP.EQ.'d*h ') THEN
C
C  special case - Linear molecule
C
         CALL LinTOPO(NAtoms, ICNNCT, IType,  NPrim,  ktyp,
     $                klist,  IErr)
c
         If(IErr.NE.0) RETURN
         NVib = NAtoms-1
c
        ELSE
C
C  attempt automatic assignment of primitive internal
C  coordinates by topological analysis
C
         CALL TOPOLOGY(NAtoms, NMol,   IMOL,   ITors,  NIC,
     $                 GROUP,  NQ,     ICNNCT, IPRNT,  NPrim,
     $                 ktyp,   klist)
c
         CALL ChkANG(NAtoms, NMol,   IMOL,   XC,     PThrsh,
     $               IType,  ICNNCT, IPRNT,  Z(1),  Z(1+NIC),
     $               NPrim,  ktyp,   klist,  IErr)
         If(IErr.NE.0) RETURN
c
         NPrim0  = NPrim      ! save # of intramolecular primitives
c
         If(IType.EQ.2) Then
C
C  inverse-distance intermolecular cluster coordinates
C
          If(IPRNT.GT.1) WRITE(IOut,1000) BSkal
          CALL DTOPOLOGY(NAtoms, NMol,   IMOL,   XC,     NIC,
     $                   NVib,   CutOff, IPRNT,  ICNNCT, NPrim,
     $                   ktyp,   klist)
c
         Else If(IType.EQ.1) Then
C
C  surface adsorption
C
          CALL SurfINT(NAtoms, NMol,   IMOL,   XC,     PThrsh,
     $                 NIC,    CutOff, IPRNT,  ICNNCT, NPrim,
     $                 ktyp,   klist)
c
         EndIf
C
C  check for fixed atoms
C
          If(NFix.GT.0)
     $       CALL FixCON(NAtoms, XC,     NFix,   IFIX,   NPrim,
     $                   ktyp,   klist,  IPRNT,  NCons,  ICTYP,
     $                   ICON,   RCON)
C
        ENDIF
C
C  relate any user-defined constraints to primitives generated
C
        If(NCons.GT.0)
     $     CALL MapCON(NCons,  ICTYP,  ICON,   RCON,   NPComp,
     $                 IPComp, IPRNT,  NPrim,  ktyp,   klist,
     $                 ICNUM,  IPCNUM)
C
C  relate any user-defined coordinate driving to primitives generated
C
        If(NDrive.GT.0)
     $     CALL MapDRIVE(NDrive, IDRTYP, IDRIVE, LDRIVE, IPrnt,
     $                   NPrim,  ktyp,   klist)
C
       ENDIF
C
C  On exit from topology section the actual number of primitive
C  internal coordinates found is NPrim
C
C  ...........................................................
C    B-MATRIX CONSTRUCTION FOR PRIMITIVE INTERNALS
C
 50   CONTINUE
C
C  initialize primitive weights
C
        DO 10 I=1,NPrim
        Coeff(I) = One
 10     CONTINUE
C  ------------------------------
C
C  Allocate memory for B-Matrix construction
C
       IMem = 2*NPrim + (NPrim*(NPrim+1))/2
c
       iptr = 1
       IErr = NMem - IMem
       If(IErr.LT.0) CALL MemERR(8*IMem,8,'MakeINTC')
C
C  allocate scratch pointers
C
       ISc1 = iptr
       ISc2 = ISc1 + NPrim
       IGG  = ISc2 + NPrim
       IEnd = IGG  + (NPrim*(NPrim+1))/2 - iptr
       CALL MemCHK(NMem,IEnd,8,'MakeINTC')
c
       itor = 1
C
C  .................................................................
C    MAKE THE DELOCALIZED INTERNAL COORDINATES
C    form B*B(t) and diagonalize
C
       call secund(t1)
       CALL BSPARSE(NAtoms, XC,     NPrim,  ktyp,   klist,
     $              Coeff,  SavTOR, 1,      1,      itor,
     $              BSkal,  IPRNT,  BPRIM,  INDB,   XPrim,  IErr)
       If(IErr.NE.0) RETURN
c
       CALL MakeDIC(NAtoms, NPrim,  NVib,   IPRNT,  ktyp,
     $              BPRIM,  INDB,   Z(IGG), Z(ISc1),Z(ISc2),
     $              UT,     IErr)
       If(IErr.NE.0) RETURN
C
C  Check weight of each primitive in the non-redundant space
C  If less than GWght, eliminate the primitive
C
       IF(.NOT.Checked) THEN
        Checked = .True.
        CALL GetWGHT(NPrim, NVib, NCons, IPRNT, UT, XPrim)
c
c -- XPrim used as scratch
c -- on exit, contains weight of each primitive
c
        MPrim = 0
        DO 20 I=1,NPrim
        IF(XPrim(I).GT.GWght) Then
          MPrim = MPrim+1
          ktyp(MPrim) = ktyp(I)
          klist(1,MPrim) = klist(1,I)
          klist(2,MPrim) = klist(2,I)
          klist(3,MPrim) = klist(3,I)
          klist(4,MPrim) = klist(4,I)
        ENDIF
 20     CONTINUE
c
        If(MPrim.LT.NPrim) Then
         WRITE(IOut,1200) NPrim-MPrim,GWght
         NPrim = MPrim
         GO TO 50
        Else
         WRITE(IOut,1300)
        EndIf
c
       ENDIF
C 
C  we're done
C
C
      ELSE
C
C =========================================================
C  N A T U R A L   I N T E R N A L   C O O R D I N A T E S
C =========================================================
C
C
C  ** NEED TO REDO THIS ENTIRE SECTION **
C     Will probably use redundant primitives
C
       WRITE(IOut,2000)
cc
      ENDIF
C
C  ................................................................
C
C  Coordinates generated
C  Can we do Z-Matrix back-transformation?
C
      If(IOrd.EQ.-1) RETURN
c
      IGZ = iptr
      IGI = IGZ + 3*NAtoms
      ICN = IGI + 5*NAtoms
      IEnd = ICN + NAtoms*NAtoms - iptr
      CALL MemCHK(NMem,IEnd,8,'MakeINTC')
C
      IF(IType.EQ.0) THEN
       CALL Tor2ZN(NAtoms, NPrim,  IGen,   ktyp,   klist,
     $             XPrim,  XC,     IPRNT,  Z(ICN), Z(IGI),
     $             Z(IGZ), IMAP,   IG,     IOrd)
C
C  copy Z(IGI) into IGEO array
C
       If(IOrd.NE.-1) CALL ICpyVEC(4*NAtoms,Z(IGI),IGEO)
      ENDIF
C
      RETURN
c
 1000 FORMAT(' Inverse Distance Scaling Factor is ',F10.6)
 1100 FORMAT(/,2X,'***ERROR*** Unable to Fully Determine Atomic',
     $            ' Connectivity',/,5X,'Missing Homonuclear Bond',
     $            ' Distance in Bond Length Table')
 1200 FORMAT(' Eliminated ',I6,' primitives with weight less than ',
     $         F8.6)
 1300 FORMAT(' No primitives eliminated')
 2000 FORMAT(/,2X,'***ERROR*** Redundant Internal Coordinate ',
     $            ' Optimization Not Yet Implemented')
c
      END
c ======================================================================
c
      SUBROUTINE MapCON(NCon,   ICTYP,  ICON,   RCON,   NPComp,
     $                  IPCOMP, IPRNT,  intcor, ktyp,   klist,
     $                  ICNUM,  IPCNUM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine which primitive internal corresponds to
C  the desired constraints
C
C  ARGUMENTS
C
C  NCon    -  number of constraints (fixed primitives)
C  ICTYP   -  constraint type
C               1 - fixed distance               stre
C               2 - fixed bond angle             bend
C               3 - fixed out-of-plane bend      outp
C               4 - fixed dihedral angle         tors
C               5 - fixed colinear bend          linc
C               6 - fixed perpendicular bend     linp
C               9 - composite constraint
C  ICON    -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C  RCON    -  constraint values
C  NPComp  -  number of primitives in composite constraints
C  IPCOMP  -  constraint type and atoms involved in constraint
C  IPRNT   -  print flag
C  intcor  -  total number of primitives so far
C             WARNING: could be increased by this routine
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
C  ICNUM   -  on exit indicates which primitive in list
C             correspond to the desired constraints
C  IPCNUM  -  on exit indicates which primitive in list
C             correspond to the desired composite constraints
C
C
      DIMENSION ICTYP(NCon),ICON(4,NCon),RCON(NCon),ICNUM(NCon),
     $          IPCOMP(5,NPComp),IPCNUM(NPComp),
     $          ktyp(intcor),klist(4,intcor)
C
C
      IOut = ioutfil('iout')
C
C  see if we can map the constraints
C  ** WARNING **  Primitives can now come in any order
C
      intcor0 = intcor
c
      DO 90 I=1,NCon
c
      IF(ICTYP(I).EQ.1) THEN
cc
       DO 20 J=1,intcor0
       IF(ktyp(J).EQ.1) THEN
       If( (ICON(1,I).EQ.klist(1,J).AND.ICON(2,I).EQ.klist(2,J)) .OR.
     $     (ICON(1,I).EQ.klist(2,J).AND.ICON(2,I).EQ.klist(1,J))) Then
        ICNUM(I) = J
        GO TO 90
       EndIf
       ENDIF
 20    CONTINUE
C
C  Couldn't find stretch among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1000) ICON(1,I),ICON(2,I)
       intcor = intcor+1
       ktyp(intcor) = 1
       klist(1,intcor) = ICON(1,I)
       klist(2,intcor) = ICON(2,I)
       ICNUM(I) = intcor
cc
      ELSE IF(ICTYP(I).EQ.2) THEN
cc
       DO 30 J=1,intcor0
       IF(ktyp(J).EQ.2) THEN
       If( (ICON(1,I).EQ.klist(1,J).AND.ICON(2,I).EQ.klist(3,J).AND.
     $      ICON(3,I).EQ.klist(2,J)) .OR.
     $     (ICON(1,I).EQ.klist(2,J).AND.ICON(2,I).EQ.klist(3,J).AND.
     $      ICON(3,I).EQ.klist(1,J)) ) Then
        ICNUM(I) = J
        GO TO 90
       EndIf
       ENDIF
 30    CONTINUE
C
C  Couldn't find bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1100) ICON(1,I),ICON(2,I),ICON(3,I)
       intcor = intcor+1
       ktyp(intcor) = 2
       klist(1,intcor) = ICON(1,I)
       klist(2,intcor) = ICON(3,I)
       klist(3,intcor) = ICON(2,I)
       ICNUM(I) = intcor
cc
      ELSE IF(ICTYP(I).EQ.3) THEN
cc
       DO 40 J=1,intcor0
       IF(ktyp(J).EQ.3) THEN
       If( (ICON(1,I).EQ.klist(1,J).AND.ICON(2,I).EQ.klist(2,J).AND.
     $      ICON(3,I).EQ.klist(3,J).AND.ICON(4,I).EQ.klist(4,J))) Then
        ICNUM(I) = J
        GO TO 90
       Else If((ICON(1,I).EQ.klist(1,J).AND.ICON(3,I).EQ.klist(2,J).AND.
     $       ICON(2,I).EQ.klist(3,J).AND.ICON(4,I).EQ.klist(4,J))) Then
        ICNUM(I) = J
        RCON(I) = -RCON(I)      ! reverse sign of out-of-plane bend
       EndIf
       ENDIF
 40    CONTINUE
C
C  Couldn't find out-of-plane bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1200) ICON(1,I),ICON(2,I),
     $                                 ICON(3,I),ICON(4,I)
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = ICON(1,I)
       klist(2,intcor) = ICON(2,I)
       klist(3,intcor) = ICON(3,I)
       klist(4,intcor) = ICON(4,I)
       ICNUM(I) = intcor
cc
      ELSE IF(ICTYP(I).EQ.4) THEN
cc
       DO 50 J=1,intcor0
       IF(ktyp(J).EQ.4) THEN
       If( (ICON(1,I).EQ.klist(1,J).AND.ICON(2,I).EQ.klist(2,J).AND.
     $      ICON(3,I).EQ.klist(3,J).AND.ICON(4,I).EQ.klist(4,J)) .OR.
     $     (ICON(1,I).EQ.klist(4,J).AND.ICON(2,I).EQ.klist(3,J).AND.
     $      ICON(3,I).EQ.klist(2,J).AND.ICON(4,I).EQ.klist(1,J))) Then
        ICNUM(I) = J
        GO TO 90
       EndIf
       ENDIF
 50    CONTINUE
C
C  Couldn't find torsion among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1300) ICON(1,I),ICON(2,I),
     $                                 ICON(3,I),ICON(4,I)
       intcor = intcor+1
       ktyp(intcor) = 4
       klist(1,intcor) = ICON(1,I)
       klist(2,intcor) = ICON(2,I)
       klist(3,intcor) = ICON(3,I)
       klist(4,intcor) = ICON(4,I)
       ICNUM(I) = intcor
cc
      ELSE IF(ICTYP(I).EQ.5.OR.ICTYP(I).EQ.6) THEN
cc
       DO 60 J=1,intcor0
       If(ktyp(J).EQ.ICTYP(I).AND.
     $    ICON(1,I).EQ.klist(1,J).AND.ICON(2,I).EQ.klist(3,J).AND.
     $    ICON(3,I).EQ.klist(2,J).AND.ICON(4,I).EQ.klist(4,J)) Then
        ICNUM(I) = J
        GO TO 90
       EndIf
 60    CONTINUE
C
C  Couldn't find linear coplanar bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2) Then
        If(ICTYP(I).EQ.5) WRITE(IOut,1400) ICON(1,I),ICON(2,I),
     $                                     ICON(3,I),ICON(4,I)
        If(ICTYP(I).EQ.6) WRITE(IOut,1500) ICON(1,I),ICON(2,I),
     $                                     ICON(3,I),ICON(4,I)
       EndIf
c
       intcor = intcor+1
       ktyp(intcor) = ICTYP(I)
       klist(1,intcor) = ICON(1,I)
       klist(2,intcor) = ICON(3,I)
       klist(3,intcor) = ICON(2,I)
       klist(4,intcor) = ICON(4,I)
       ICNUM(I) = intcor
cc
      ELSE IF(ICTYP(I).EQ.7) THEN       ! temporary
cc
       DO 70 J=1,intcor0
       IF(ktyp(J).EQ.7) THEN
       If( (ICON(1,I).EQ.klist(1,J).AND.ICON(2,I).EQ.klist(2,J)) .OR.
     $     (ICON(1,I).EQ.klist(2,J).AND.ICON(2,I).EQ.klist(1,J))) Then
        ICNUM(I) = J
        GO TO 90
       EndIf
       ENDIF
 70    CONTINUE
C
C  Couldn't find inverse stretch among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1700) ICON(1,I),ICON(2,I)
       intcor = intcor+1
       ktyp(intcor) = 7
       klist(1,intcor) = ICON(1,I)
       klist(2,intcor) = ICON(2,I)
       ICNUM(I) = intcor
cc
      ELSE IF(ICTYP(I).EQ.9) THEN       ! composite constraint
cc
       ICNUM(I) = 0
c
       DO 900 K=1,NPComp
c
       IF(IPCOMP(1,K).EQ.1) THEN
cc
       DO 200 J=1,intcor0
       IF(ktyp(J).EQ.1) THEN
       If( (IPCOMP(2,K).EQ.klist(1,J).AND.IPCOMP(3,K).EQ.klist(2,J)).OR.
     $   (IPCOMP(2,K).EQ.klist(2,J).AND.IPCOMP(3,K).EQ.klist(1,J))) Then
        IPCNUM(K) = J
        GO TO 900
       EndIf
       ENDIF
 200   CONTINUE
C
C  Couldn't find stretch among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1000) IPCOMP(2,K),IPCOMP(3,K)
       intcor = intcor+1
       ktyp(intcor) = 1
       klist(1,intcor) = IPCOMP(2,K)
       klist(2,intcor) = IPCOMP(3,K)
       IPCNUM(K) = intcor
cc
       ELSE IF(IPCOMP(1,K).EQ.2) THEN
cc
       DO 300 J=1,intcor0
       IF(ktyp(J).EQ.2) THEN
       If( (IPCOMP(2,K).EQ.klist(1,J).AND.IPCOMP(3,K).EQ.klist(3,J).AND.
     $      IPCOMP(4,K).EQ.klist(2,J)) .OR.
     $     (IPCOMP(2,K).EQ.klist(2,J).AND.IPCOMP(3,K).EQ.klist(3,J).AND.
     $      IPCOMP(4,K).EQ.klist(1,J)) ) Then
        IPCNUM(K) = J
        GO TO 900
       EndIf
       ENDIF
 300   CONTINUE
C
C  Couldn't find bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2)
     $    WRITE(IOut,1100) IPCOMP(2,K),IPCOMP(3,K),IPCOMP(4,K)
       intcor = intcor+1
       ktyp(intcor) = 2
       klist(1,intcor) = IPCOMP(2,K)
       klist(2,intcor) = IPCOMP(4,K)
       klist(3,intcor) = IPCOMP(3,K)
       IPCNUM(K) = intcor
cc
       ELSE IF(IPCOMP(1,K).EQ.3) THEN
cc
       DO 400 J=1,intcor0
       IF(ktyp(J).EQ.3) THEN
       If(IPCOMP(2,K).EQ.klist(1,J).AND.IPCOMP(3,K).EQ.klist(2,J).AND.
     $    IPCOMP(4,K).EQ.klist(3,J).AND.IPCOMP(5,K).EQ.klist(4,J)) Then
        IPCNUM(K) = J
        GO TO 900
       Else If(IPCOMP(2,K).EQ.klist(1,J).AND.IPCOMP(4,K).EQ.klist(2,J)
     $.AND.IPCOMP(3,K).EQ.klist(3,J).AND.IPCOMP(5,K).EQ.klist(4,J)) Then
        IPCNUM(K) = J
       EndIf
       ENDIF
 400   CONTINUE
C
C  Couldn't find out-of-plane bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1200) IPCOMP(2,K),IPCOMP(3,K),
     $                                 IPCOMP(4,K),IPCOMP(5,K)
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = IPCOMP(2,K)
       klist(2,intcor) = IPCOMP(3,K)
       klist(3,intcor) = IPCOMP(4,K)
       klist(4,intcor) = IPCOMP(5,K)
       IPCNUM(K) = intcor
cc
       ELSE IF(IPCOMP(1,K).EQ.4) THEN
cc
       DO 500 J=1,intcor0
       IF(ktyp(J).EQ.4) THEN
       If( (IPCOMP(2,K).EQ.klist(1,J).AND.IPCOMP(3,K).EQ.klist(2,J).AND.
     $     IPCOMP(4,K).EQ.klist(3,J).AND.IPCOMP(5,K).EQ.klist(4,J)) .OR.
     $     (IPCOMP(2,K).EQ.klist(4,J).AND.IPCOMP(3,K).EQ.klist(3,J).AND.
     $    IPCOMP(4,K).EQ.klist(2,J).AND.IPCOMP(5,K).EQ.klist(1,J))) Then
        IPCNUM(K) = J
        GO TO 900
       EndIf
       ENDIF
 500   CONTINUE
C
C  Couldn't find torsion among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1300) IPCOMP(2,K),IPCOMP(3,K),
     $                                 IPCOMP(4,K),IPCOMP(5,K)
       intcor = intcor+1
       ktyp(intcor) = 4
       klist(1,intcor) = IPCOMP(2,K)
       klist(2,intcor) = IPCOMP(3,K)
       klist(3,intcor) = IPCOMP(4,K)
       klist(4,intcor) = IPCOMP(5,K)
       IPCNUM(K) = intcor
cc
       ENDIF
 900   CONTINUE
cc
      ENDIF
 90   CONTINUE
c
      If(IPRNT.GT.2.AND.intcor.NE.intcor0) WRITE(IOut,1600) intcor
C
      RETURN
c
 1000 FORMAT(' Adding Distance ',2I4,9X,' to Primitive Set',
     $       ' in order to Constrain it')
 1100 FORMAT(' Adding Bend     ',3I4,5X,' to Primitive Set',
     $       ' in order to Constrain it')
 1200 FORMAT(' Adding OOP Bend ',4I4,1X,' to Primitive Set',
     $       ' in order to Constrain it')
 1300 FORMAT(' Adding Torsion  ',4I4,1X,' to Primitive Set',
     $       ' in order to Constrain it')
 1400 FORMAT(' Adding LINC Bend',4I4,1X,' to Primitive Set',
     $       ' in order to Constrain it')
 1500 FORMAT(' Adding LINP Bend',4I4,1X,' to Primitive Set',
     $       ' in order to Constrain it')
 1600 FORMAT(' There are now ',I6,' Primitive Internals')
 1700 FORMAT(' Adding Inv.Dist. ',2I4,9X,' to Primitive Set',
     $       ' in order to Constrain it')
c
      END
c ======================================================================
c
      SUBROUTINE MapDRIVE(NDrive, IDRTYP, IDRV,   LDRIVE, IPrnt,
     $                    intcor, ktyp,   klist)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine which primitive internal corresponds to
C  the desired driving coordinate
C
C  ARGUMENTS
C
C  NDrive  -  number of primitives being driven
C  IDRTYP  -  type of each primitive being driven
C               1 - distance               stre
C               2 - bond angle             bend
C               3 - out-of-plane bend      outp
C               4 - dihedral angle         tors
C  IDRV    -  definition of each primitive (i.e., list of atoms)
C  LDRIVE  -  on exit primitive number, i.e., which primitive we are
C              referring to in list of all primitives
C  IPrnt   -  print flag
C  intcor  -  total number of primitives so far
C             WARNING: could be increased by this routine
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
C
C
      DIMENSION IDRTYP(NDrive),IDRV(4,NDrive),LDRIVE(NDrive),
     $          ktyp(intcor),klist(4,intcor)
C
C
      IOut = ioutfil('iout')
C
C  see if we can map the coordinates being driven
C  ** WARNING **  Primitives can now come in any order
C
      intcor0 = intcor
c
      DO 90 I=1,NDrive
c
      IF(IDRTYP(I).EQ.1) THEN
cc
       DO 20 J=1,intcor0
       IF(ktyp(J).EQ.1) THEN
       If( (IDRV(1,I).EQ.klist(1,J).AND.IDRV(2,I).EQ.klist(2,J)) .OR.
     $     (IDRV(1,I).EQ.klist(2,J).AND.IDRV(2,I).EQ.klist(1,J))) Then
        LDRIVE(I) = J
        GO TO 90
       EndIf
       ENDIF
 20    CONTINUE
C
C  Couldn't find stretch among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1000) IDRV(1,I),IDRV(2,I)
       intcor = intcor+1
       ktyp(intcor) = 1
       klist(1,intcor) = IDRV(1,I)
       klist(2,intcor) = IDRV(2,I)
       LDRIVE(I) = intcor
cc
      ELSE IF(IDRTYP(I).EQ.2) THEN
cc
       DO 30 J=1,intcor0
       IF(ktyp(J).EQ.2) THEN
       If( (IDRV(1,I).EQ.klist(1,J).AND.IDRV(2,I).EQ.klist(3,J).AND.
     $      IDRV(3,I).EQ.klist(2,J)) .OR.
     $     (IDRV(1,I).EQ.klist(2,J).AND.IDRV(2,I).EQ.klist(3,J).AND.
     $      IDRV(3,I).EQ.klist(1,J)) ) Then
        LDRIVE(I) = J
        GO TO 90
       EndIf
       ENDIF
 30    CONTINUE
C
C  Couldn't find bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1100) IDRV(1,I),IDRV(2,I),IDRV(3,I)
       intcor = intcor+1
       ktyp(intcor) = 2
       klist(1,intcor) = IDRV(1,I)
       klist(2,intcor) = IDRV(3,I)
       klist(3,intcor) = IDRV(2,I)
       LDRIVE(I) = intcor
cc
      ELSE IF(IDRTYP(I).EQ.3) THEN
cc
       DO 40 J=1,intcor0
       IF(ktyp(J).EQ.3) THEN
       If( (IDRV(1,I).EQ.klist(1,J).AND.IDRV(2,I).EQ.klist(2,J).AND.
     $      IDRV(3,I).EQ.klist(3,J).AND.IDRV(4,I).EQ.klist(4,J))) Then
        LDRIVE(I) = J
        GO TO 90
       Else If((IDRV(1,I).EQ.klist(1,J).AND.IDRV(3,I).EQ.klist(2,J).AND.
     $       IDRV(2,I).EQ.klist(3,J).AND.IDRV(4,I).EQ.klist(4,J))) Then
        LDRIVE(I) = J
cc        RCON(I) = -RCON(I)      ! reverse sign of out-of-plane bend
       EndIf
       ENDIF
 40    CONTINUE
C
C  Couldn't find out-of-plane bend among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1200) IDRV(1,I),IDRV(2,I),
     $                                 IDRV(3,I),IDRV(4,I)
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = IDRV(1,I)
       klist(2,intcor) = IDRV(2,I)
       klist(3,intcor) = IDRV(3,I)
       klist(4,intcor) = IDRV(4,I)
       LDRIVE(I) = intcor
cc
      ELSE IF(IDRTYP(I).EQ.4) THEN
cc
       DO 50 J=1,intcor0
       IF(ktyp(J).EQ.4) THEN
       If( (IDRV(1,I).EQ.klist(1,J).AND.IDRV(2,I).EQ.klist(2,J).AND.
     $      IDRV(3,I).EQ.klist(3,J).AND.IDRV(4,I).EQ.klist(4,J)) .OR.
     $     (IDRV(1,I).EQ.klist(4,J).AND.IDRV(2,I).EQ.klist(3,J).AND.
     $      IDRV(3,I).EQ.klist(2,J).AND.IDRV(4,I).EQ.klist(1,J))) Then
        LDRIVE(I) = J
        GO TO 90
       EndIf
       ENDIF
 50    CONTINUE
C
C  Couldn't find torsion among existing primitives
C  so add it
C
       If(IPRNT.GT.2) WRITE(IOut,1300) IDRV(1,I),IDRV(2,I),
     $                                 IDRV(3,I),IDRV(4,I)
       intcor = intcor+1
       ktyp(intcor) = 4
       klist(1,intcor) = IDRV(1,I)
       klist(2,intcor) = IDRV(2,I)
       klist(3,intcor) = IDRV(3,I)
       klist(4,intcor) = IDRV(4,I)
       LDRIVE(I) = intcor
cc
      ENDIF
 90   CONTINUE
c
      If(IPRNT.GT.2.AND.intcor.NE.intcor0) WRITE(IOut,1600) intcor
C
      RETURN
c
 1000 FORMAT(' Adding Distance ',2I4,9X,' to Primitive Set',
     $       ' in order to Drive it')
 1100 FORMAT(' Adding Bend     ',3I4,5X,' to Primitive Set',
     $       ' in order to Drive it')
 1200 FORMAT(' Adding OOP Bend ',4I4,1X,' to Primitive Set',
     $       ' in order to Drive it')
 1300 FORMAT(' Adding Torsion  ',4I4,1X,' to Primitive Set',
     $       ' in order to Drive it')
 1600 FORMAT(' There are now ',I6,' Primitive Internals')
c
      END
c ======================================================================
c
      SUBROUTINE ModeOverlap(NCycle, N,      mode,   Nmode,  IPRNT,
     $                       U,      VMODE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine which Hessian mode to follow during a transition
C  state search.
C
C  On the first cycle this is either set by the user
C  or (default) is the lowest mode. On subsequent cycles
C  it is the mode with the maximum overlap with the mode
C  followed on the previous cycle.
C
C  ** NOTE: When the mode being followed becomes the lowest
C           mode, mode following is switched off
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C  N       -  full dimension of space
C  mode    -  mode followed on previous cycle
C             (or initial mode if first cycle)
C  Nmode   -  mode followed on this cycle
C  IPRNT   -  controls level of print out
C  U       -  Hessian eigenvectors
C  VMODE   -  eigenvector followed on previous cycle
C             (if mode is non-zero)
C             on exit contains mode followed on this cycle
C
      REAL*8 U(N,*),VMODE(N)
C
C
      IOut = ioutfil('iout')
c
      IF(NCycle.EQ.1) THEN
cc
       IT = mode
       If(IPRNT.GT.1) WRITE(IOut,1000) mode
cc
      ELSE
cc
       IT = 1
       TOVLP = SProd(N,U(1,1),VMODE)
       TOVLP = Abs(TOVLP)
C
       DO 10 I=2,N
       OVLP = SProd(N,U(1,I),VMODE)
       OVLP = Abs(OVLP)
       IF(OVLP.GT.TOVLP) THEN
        TOVLP = OVLP
        IT = I
       ENDIF
 10    CONTINUE
C
       If(IPRNT.GT.1) WRITE(IOut,1100) TOVLP
cc
      ENDIF
C
C  Store the mode to be followed in VMODE
C
      DO 20 I=1,N
      VMODE(I) = U(I,IT)
 20   CONTINUE
C
      Nmode = IT
      RETURN
c
 1000 FORMAT(' Hessian Mode Following Switched On',/,
     $       ' Following mode: ',I3)
 1100 FORMAT(' Overlap of Current mode with previous mode is ',F12.6)
c
      END
c ======================================================================
c
      SUBROUTINE MODPrim(NPrim,  NPrim1, MPrim,  IHPrim, INDX,
     $                   V1,     V2,     HPTmp,  GOld,   DOld,
     $                   HPOld,  GNew,   HPRIM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Modifies previous gradient, displacement and Hessian over
C  primitive internals to account for changes in the underlying
C  primitive space
C
C  ARGUMENTS
C
C  NPrim   -  size of current primitive internal space
C  NPrim1  -  size of previous primitive internal space
C  MPrim   -  total number of entries in INDX array
C  IHPrim  -  primitive Hessian flag
C             in particular if IHess=1, set new Hessian entries to unity
C  INDX    -  contains array of differences between two primitive sets
C             (specifically what to do to first set to get second)
C                0 - primitive common to both sets
C               -1 - primitive needs to be deleted from first set
C               +1 - primitive needs to be added to first set
C  V1      -  scratch vector for sorting gradient
C  V2      -  scratch vector for sorting displacement
C  HPTmp   -  scratch matrix for sorting Hessian
C  GOld    -  previous gradient over old/new primitives
C  DOld    -  displacement over old/new primitives
C  HPOld   -  Hessian over old primitives
C  GNew    -  current gradient over new primitives
C  HPRIM   -  Hessian over new primitives
C
C  **WARNING**
C    HPRIM and HPOld can share storage (see calling routine)
C
C
      REAL*8 V1(NPrim1),V2(NPrim1),HPTmp(NPrim1,NPrim),
     $       GOld(NPrim),DOld(NPrim),HPOld(NPrim1,NPrim1),
     $       GNew(NPrim),HPRIM(NPrim,NPrim)
      DIMENSION INDX(MPrim)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
C  initialize
C
      If(IHPrim.EQ.1) Then
       HDiag = One
      Else
       HDiag = 0.2d0
      EndIf
c
      CALL CpyVEC(NPrim1,GOld,V1)
      CALL CpyVEC(NPrim1,DOld,V2)
C
C  first sort gradient and displacement
C
      I1 = 0
      I2 = 0
      DO 10 I=1,MPrim
      IF(INDX(I).EQ.0) THEN
       I1 = I1+1
       I2 = I2+1
       GOld(I2) = V1(I1)
       DOld(I2) = V2(I1)
      ELSE IF(INDX(I).EQ.1) THEN
       I2 = I2+1
       GOld(I2) = Zero
       DOld(I2) = Zero
      ELSE
       I1 = I1+1
      ENDIF
 10   CONTINUE
C
C  now sort Hessian
C
C  (i)  Rows
C
      I1 = 0
      I2 = 0
      DO 20 J=1,MPrim
      IF(INDX(J).EQ.0) THEN
       I1 = I1+1
       I2 = I2+1
       CALL CpyVEC(NPrim1,HPOld(1,I1),HPTmp(1,I2))
      ELSE IF(INDX(J).EQ.1) THEN
       I2 = I2+1
       CALL ZeroIT(HPTmp(1,I2),NPrim1)
      ELSE
       I1 = I1+1
      ENDIF
 20   CONTINUE
C
C  (ii)  Columns
C
      DO 40 I=1,NPrim
      I1 = 0
      I2 = 0
      DO 30 J=1,MPrim
      IF(INDX(J).EQ.0) THEN
       I1 = I1+1
       I2 = I2+1
       HPRIM(I2,I) = HPTmp(I1,I)
      ELSE IF(INDX(J).EQ.1) THEN
       I2 = I2+1
       HPRIM(I2,I) = Zero
      ELSE
       I1 = I1+1
      ENDIF
 30   CONTINUE
 40   CONTINUE
C
C  now add diagonal elements for new primitives
C  and zero corresponding gradient entries
C
      I2 = 0
      DO 50 I=1,MPrim
      I2 = I2+1
      IF(INDX(I).EQ.1) THEN
       HPRIM(I2,I2) = HDiag
       GNew(I2) = Zero
      ELSE IF(INDX(I).EQ.-1) THEN
       I2 = I2-1
      ENDIF
 50   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE NewCART(NAtoms,NIC,BINV,XINT,QQ,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Gets new Cartesian coordinates from current estimate
C  of Cartesians and internals
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NIC     -  number of internal coordinates
C  BINV    -  inverse of B-Matrix
C  XINT    -  set of exact internal coordinates
C  QQ      -  estimate of internal coordinates
C  XC      -  estimate of corresponding Cartesian Coordinates
C             on exit contains new estimate
C
      REAL*8 BINV(3*NAtoms,NIC),XINT(NIC),QQ(NIC),XC(3*NAtoms)
C
      DO 20 I=1,NIC
      QDIF = XINT(I) - QQ(I)
      DO 10 J=1,3*NAtoms
      XC(J) = XC(J) + BINV(J,I)*QDIF
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE OPTEF(NCycle, NC,     N,      NCon,   NEG,
     $                 NegReq, mode,   IPRNT,  U,      Eigval,
     $                 G,      FX,     VMODE,  D,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Driving Routine for Eigenvector (mode) Following algorithm
C
C  Based on:
C
C   "On Finding Transition States"
C    Cerjan & Miller, J.Chem.Phys. 75 (1981) 2800
C
C   "Walking on Potential Energy Surfaces"
C    Simons, Jorgenson, Taylor & Ozment, J.Phys.Chem. 87 (1983) 2745
C
C   "Searching for Stationary Points on Surfaces"
C    Banerjee, Adams, Simons & Shepard, J.Phys.Chem. 89 (1985) 52
C
C  For full details see:
C
C   "An Algorithm for the Location of Transition States"
C    J.Baker, J.Comp.Chem. 7 (1986) 385
C
C
C  ARGUMENTS
C
C  NCycle  -  current optimization cycle number
C  NC      -  number of modes used to form new step
C  N       -  actual number of coordinates
C             (i.e. full dimension of problem)
C  NCon    -  number of constraints
C             (for constrained Cartesian coordinate optimization)
C  NEG     -  number of negative eigenvalues of Hessian
C  NegReq  -  actual number of negative eigenvalues desired
C  mode    -  mode being followed during Transition State search
C             (zero if mode following switched off)
C  IPRNT   -  controls level of print out
C  U       -  Hessian eigenvectors
C  EigVal  -  Hessian eigenvalues
C  G       -  gradient vector
C  VMODE   -  eigenvector followed on previous cycle
C             (if mode is non-zero)
C
C  On exit
C
C  FX      -  gradient in Hessian eigenvector basis
C  VMODE   -  eigenvector followed this cycle
C  D       -  new coordinate displacement (next step)
C  IErr    -  error flag    0 - new step calculated successfully
C                          -1 - something went wrong
C
C
      REAL*8 U(N,*),EigVal(NC),G(N),FX(NC),VMODE(N),D(N)
C
C
C  Use the selected NC eigenvectors of the Hessian (stored in U)
C  to form the next step
C
      IOut = ioutfil('iout')
C
C  First form the FX vector
C  (the components of G along the local Hessian modes)
C
      DO 10 I=1,NC
      FX(I) = SProd(N,U(1,I),G)
 10   CONTINUE
C
      IF(NCon.EQ.0) THEN
cc
C  Normal Unconstrained Optimization
C
C  Take the P-RFO step for a TS search
C  Take the simple RFO step for a minimum search
C
       IF(IPRNT.GT.1) THEN
        If(NEG.NE.NegReq) WRITE(IOut,1000)
        If(NegReq.EQ.0) WRITE(IOut,1100)
        If(NegReq.EQ.1) WRITE(IOut,1200)
       ENDIF
c
       CALL FormD(NCycle, NC,     N,      NegReq, mode,
     $            IPRNT,  U,      EigVal, FX,     VMODE,
     $            D,      IErr )
C
C  If the RFO step failed then fall back on the
C  simple Newton-Raphson step if possible
C
       IF(IErr.EQ.-1.AND.NEG.EQ.NegReq) THEN
        If(IPRNT.GT.1) WRITE(IOut,1300)
        CALL FormNR(NC,N,U,EigVal,FX,D)
        IErr = 0
       ENDIF
cc
      ELSE
cc
C  Constrained optimization in Cartesian coordinates
C
C  Take the P-RFO step
C
       IF(IPRNT.GT.1) THEN
        If(NEG.NE.NegReq) WRITE(IOut,1000)
        If(NegReq.EQ.NCon) WRITE(IOut,1400)
        If(NegReq.EQ.NCon+1) WRITE(IOut,1200)
       ENDIF
c
       CALL ConFormD(NCycle, NC,     N,      NCon,   NegReq,
     $               mode,   IPRNT,  U,      EigVal, FX,
     $               VMODE,  D,      IErr )
cc
      ENDIF
C
      RETURN
c
 1000 FORMAT('**WARNING** Hessian does not have the Desired',
     $       ' Local Structure')
 1100 FORMAT(/,' Minimum Search - Taking Simple RFO Step')
 1200 FORMAT(/,' Transition State Search - Taking P-RFO Step')
 1300 FORMAT(/,' Taking Simple Newton-Raphson Step')
 1400 FORMAT(/,' Minimum Search - Taking P-RFO Step')
c
      END
c ======================================================================
c
      SUBROUTINE PenGRAD(NAtoms, NCons,  ICTYP,  RCON,   IC,
     $                   XC,     GC,     GCC )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine modifies the standard Cartesian gradient to
C  incorporate geometric constraints (involving interatomic
C  distances, angles and dihedral angles) using quadratic
C  penalty functions
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NCons   -  number of constraints
C  ICTYP   -  integer array indicating constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane bend    CURRENTLY NOT IMPLEMENTED
C              4 - fixed dihedral angle
C  RCON    -  value of constraint
C  IC      -  list of atoms involved in the constraints
C              IC1-IC2           distance constraint
C              IC1-IC2-IC3       bond angle constraint
C              IC1-IC2-IC3-IC4   dihedral constraint
C  XC      -  coordinate vector
C  GC      -  gradient vector
C  GCC     -  matrix of constraint normals
C
C  ** WARNING **  This routine forms the negative of GCC  **
C
C
      DIMENSION ICTYP(NCons),RCON(NCons),XC(*),GC(*),
     $          IC(4,NCons),GCC(3*NAtoms,NCons)
C
      PARAMETER (One=1.0d0,sigma=10.0d0)
C
C
      NAT3 = 3*NAtoms
      CALL ZeroIT(GCC,NAT3*NCons)
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
C  now determine the current constraint multiplier
C
       IT = NAT3 + IQ
       XC(IT) = sigma*(RCON(IQ) - R)
cc
      ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  angle constraint (I-J-K)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       CALL AngGRAD(NAtoms,I,J,K,XC,Th,.true.,GCC(1,IQ))
C
C  now determine the current constraint value
C
       IT = NAT3 + IQ
       XC(IT) = sigma*(RCON(IQ) - Th)
cc
      ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral constraint (I-J-K-L)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       L = IC(4,IQ)
       CALL DihGRAD(NAtoms,I,J,K,L,XC,Dih,.true.,GCC(1,IQ))
       CALL ChkDIH(RCON(IQ),Dih)
C
C  now determine the current constraint value
C
       IT = NAT3 + IQ
       XC(IT) = sigma*(RCON(IQ) - Dih)
cc
      ENDIF
c
 10   CONTINUE
C
C  now modify the Cartesian gradient to include
C  the constraint normals
C
      DO 30 IQ=1,NCons
      RLambda = XC(NAT3+IQ)
      DO 20 I=1,NAT3
      GC(I) = GC(I) + RLambda*GCC(I,IQ)
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE PrntCON(IOut,   NAtoms, NCons,  XC,     ICTYP,
     $                   IC,     RCON,   ICOMP,  IPCOMP, PCWght)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out the constraints and their current values
C
      REAL*8 XC(3,NAtoms),RCON(NCons)
      DIMENSION ICTYP(NCons),IC(4,NCons),
     $          ICOMP(*),IPCOMP(5,*),PCWght(*)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      AUTOAN = 1.0d0/ANTOAU
      TOANG = 180.0d0/PI
C
C
      WRITE(IOut,1000)
c
      nc = 0
      np1 = 1
c
      DO 10 IQ=1,NCons
      I = IC(1,IQ)
      J = IC(2,IQ)
      K = IC(3,IQ)
      L = IC(4,IQ)
c
      IF(ICTYP(IQ).EQ.1) THEN
cc
C  distance constraint
C
       XIJ = XC(1,I) - XC(1,J)
       YIJ = XC(2,I) - XC(2,J)
       ZIJ = XC(3,I) - XC(3,J)
       R = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
c
       val = RCON(IQ)*AUTOAN
       R = R*AUTOAN
c
       WRITE(IOut,1100) I,J,R,val
cc
      ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  bond angle constraint
C
       CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
c
       val = RCON(IQ)*TOANG
       Th = Th*TOANG
c
       WRITE(IOut,1200) I,J,K,Th,val
cc
      ELSE IF(ICTYP(IQ).EQ.3) THEN
cc
C  out-of-plane bend constraint
C
       CALL OutpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       val = RCON(IQ)*TOANG
       Th = Th*TOANG
c
       WRITE(IOut,1300) I,J,K,L,Th,val
cc
      ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral angle constraint
C
       CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       val = RCON(IQ)*TOANG
       Th = Th*TOANG
c
       WRITE(IOut,1400) I,J,K,L,Th,val
cc
      ELSE IF(ICTYP(IQ).EQ.5) THEN
cc
C  linear coplanar bend constraint
C
       CALL LincGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       val = RCON(IQ)*TOANG
       Th = Th*TOANG
c
       WRITE(IOut,1500) I,J,K,L,Th,val
cc
      ELSE IF(ICTYP(IQ).EQ.6) THEN
cc
C  linear perpendicular bend constraint
C
       CALL LinpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
c
       val = RCON(IQ)*TOANG
       Th = Th*TOANG
c
       WRITE(IOut,1600) I,J,K,L,Th,val
cc
      ELSE IF(ICTYP(IQ).EQ.7) THEN
cc
C  inverse-power distance constraint
C
       XIJ = XC(1,I) - XC(1,J)
       YIJ = XC(2,I) - XC(2,J)
       ZIJ = XC(3,I) - XC(3,J)
       R = 1.0d0/SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
c
       val = RCON(IQ)*AUTOAN
       R = R*AUTOAN
c
       WRITE(IOut,1700) I,J,R,val
cc
      ELSE
cc
C  composite constraint
C
       WRITE(IOut,2000)
c
       nc = nc+1
       np2 = np1 + ICOMP(nc) - 1
       val = 0.0d0
cccccccccc
cc       write(6,*) ' In <PrntCON>'
cc       write(6,*) ' nc: ',nc
cc       write(6,*) ' np1:',np1,'  np2:',np2
cc       write(6,*) ' PCWght:',(pcwght(ii),ii=np1,np2)
cc       write(6,*) ' IPCOMP array:'
cc       do ii=np1,np2
cc       write(6,*) ' ii is ',ii
cc       write(6,*) ' IPCOMP(1,ii):',ipcomp(1,ii)
cc       write(6,*) ' IPCOMP(2,ii):',ipcomp(2,ii)
cc       write(6,*) ' IPCOMP(3,ii):',ipcomp(3,ii)
cc       write(6,*) ' IPCOMP(4,ii):',ipcomp(4,ii)
cc       write(6,*) ' IPCOMP(5,ii):',ipcomp(5,ii)
cc       enddo
cccccccccc
c
       DO 9 IP=np1,np2
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
        WRITE(IOut,1100) I,J,R,Wght
cc
       ELSE IF(ITYP.EQ.2) THEN
cc
C  bond angle constraint
C
        CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
        val = val+Th*Wght
        Th = Th*TOANG
c
        WRITE(IOut,1200) I,J,K,Th,Wght
cc
       ELSE IF(ITYP.EQ.3) THEN
cc
C  out-of-plane bend constraint
C
        CALL OutpGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
        val = val+Th*Wght
        Th = Th*TOANG
c
        WRITE(IOut,1300) I,J,K,L,Th,Wght
cc
       ELSE IF(ITYP.EQ.4) THEN
cc
C  dihedral angle constraint
C
        CALL DihGRAD(NAtoms,I,J,K,L,XC,Th,.false.,jnk)
        val = val+Th*Wght
        Th = Th*TOANG
c
        WRITE(IOut,1400) I,J,K,L,Th,Wght
cc
       ENDIF
  9    CONTINUE
c
       WRITE(IOut,2100) val,RCON(IQ)
       np1 = np2+1
cc
      ENDIF
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(14X,'Constraints and their Current Values',/,
     $       37X,'   Value     Constraint')
 1100 FORMAT('   Distance:     ',2I4,12X,F9.6,4X,F9.6)
 1200 FORMAT('   Angle:        ',3I4,8X,F9.3,4X,F9.3)
 1300 FORMAT('   Out-of-Plane: ',4I4,4X,F9.3,4X,F9.3)
 1400 FORMAT('   Dihedral:     ',4I4,4X,F9.3,4X,F9.3)
 1500 FORMAT('   Linear Plane: ',4I4,4X,F9.3,4X,F9.3)
 1600 FORMAT('   Linear Perp.: ',4I4,4X,F9.3,4X,F9.3)
 1700 FORMAT('   Inverse Dist: ',2I4,12X,F9.6,4X,F9.6)
 2000 FORMAT('  Composite Coordinate',18X,'value',8X,'weight')
 2100 FORMAT('   Current Constraint Value (au): ',F12.6,/,
     $       '   Desired Constraint Value (au): ',F12.6,/,
     $       '  End of Composite Constraint')
c
      END
c ======================================================================
c
      SUBROUTINE PrntDRIVE(IOut,NDrive,IDRTYP,IDRIVE,FDRIVE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out primitives involved in coordinate driving
C
      DIMENSION IDRTYP(NDrive),IDRIVE(4,NDrive),FDRIVE(NDrive)
C
      WRITE(IOut,1000)
c
      DO 10 I=1,NDrive
      IF(IDRTYP(I).EQ.1) THEN
        WRITE(IOut,1010) IDRIVE(1,I),IDRIVE(2,I),FDRIVE(I)
      ELSE IF(IDRTYP(I).EQ.2) THEN
        WRITE(IOut,1020) IDRIVE(1,I),IDRIVE(2,I),IDRIVE(3,I),FDRIVE(I)
      ELSE IF(IDRTYP(I).EQ.3) THEN
        WRITE(IOut,1030) IDRIVE(1,I),IDRIVE(2,I),IDRIVE(3,I),
     $                   IDRIVE(4,I),FDRIVE(I)
      ELSE IF(IDRTYP(I).EQ.4) THEN
        WRITE(IOut,1040) IDRIVE(1,I),IDRIVE(2,I),IDRIVE(3,I),
     $                   IDRIVE(4,I),FDRIVE(I)
      ENDIF
 10   CONTINUE
c
      RETURN
c
 1000 FORMAT(/,1X,' The Following Primitives are Involved in',
     $            ' Coordinate Driving'/,
     $       16X ' Primitive              Driving Force (au)')
 1010 FORMAT(4X,'stre',4X,2I6,16X,F12.6)
 1020 FORMAT(4X,'bend',4X,3I6,10X,F12.6)
 1030 FORMAT(4X,'outp',4X,4I6,6X,F12.6)
 1040 FORMAT(4X,'tors',4X,4I6,6X,F12.6)
c
      END
c ======================================================================
c
      SUBROUTINE PrntFIX(IOut,NAtoms,IFix,CFix,ISTR)
      IMPLICIT INTEGER(A-Z)
C
C  Prints out the "fixed" atoms
C
      DIMENSION IFix(3,NAtoms),ISTR(NAtoms)
      CHARACTER*1 CFix(3,NAtoms)
C
C
      WRITE(IOut,1000)
c
      DO 10 I=1,NAtoms
       ISTR(I) = 0
       CFix(1,I) = ' '
       CFix(2,I) = ' '
       CFix(3,I) = ' '
 10   CONTINUE
c
      DO 20 I=1,NAtoms
      IF(IFix(1,I).EQ.1) THEN
       ISTR(I) = 1
       CFix(1,I) = 'X'
      ENDIF
      IF(IFix(2,I).EQ.1) THEN
       ISTR(I) = 1
       CFix(2,I) = 'Y'
      ENDIF
      IF(IFix(3,I).EQ.1) THEN
       ISTR(I) = 1
       CFix(3,I) = 'Z'
      ENDIF
 20   CONTINUE
C
      DO 30 I=1,NAtoms
      If(ISTR(I).EQ.1) WRITE(IOut,1100) I,CFix(1,I),CFix(2,I),CFix(3,I)
 30   CONTINUE
C
      RETURN
c
 1000 FORMAT(12X,'Atoms with "Fixed" Coordinates',/,
     $       12X,'    Atom    Coordinates')
 1100 FORMAT(16X,I3,6X,3(2X,A1))
c
      END
c ======================================================================
c
      SUBROUTINE ProjTR(NAtoms, XC,     IProj,  IPRNT, P,
     $                  TRVec,  HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Projects out from the Hessian matrix in Cartesian
C  coordinates vectors corresponding to translations
C  and infinitesimal rotations
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  IProj   -  Flag controlling projection
C              0 - don't project
C              1 - project out translations
C              2 - project out translations and rotations
C  IPRNT   -  print flag
C  P       -  scratch space for projection matrix
C  TRVec   -  scratch space for TR vectors
C  HESS    -  on entry contains Cartesian Hessian matrix
C             on exit contains projected Hessian
C
      REAL*8 XC(3,NAtoms),P(3*NAtoms,*),TRVec(3*NAtoms,*),
     $       HESS(3*NAtoms,3*NAtoms)
      REAL*8 T(9),PMom(3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      IOut = ioutfil('iout')
C
C  Coordinates should already be in centre of mass
C  frame, but transform just in case
C
      CALL CMS(NAtoms,CX,CY,CZ,XC)
C
C  Find the principal moments & rotation generators
C
      CALL ZeroIT(T,9)
      DO 20 I=1,NAtoms
      X = XC(1,I)
      Y = XC(2,I)
      Z = XC(3,I)
      T(1) = T(1) + (Y*Y + Z*Z)
      T(5) = T(5) + (X*X + Z*Z)
      T(9) = T(9) + (X*X + Y*Y)
      T(2) = T(2) - X*Y
      T(3) = T(3) - X*Z
      T(6) = T(6) - Y*Z
 20   CONTINUE
      T(4) = T(2)
      T(7) = T(3)
      T(8) = T(6)
C
C  Diagonalize T
C
      CALL DIAGMAT(T,3,TRVec,P,PMom,IErr)    ! TRVec & P used as scratch
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
C  Set up Orthogonal coordinate vectors for translation and
C  rotation about principal axes of inertia
C
      NAT3 = 3*NAtoms
      CALL ZeroIT(TRVec,NAT3*6)
      CALL FormTR(NAtoms,XC,T,TRVec)
c
      IF(IPRNT.GT.5) THEN
       WRITE(IOut,1100)
       CALL PrntMAT(6,NAT3,6,TRVec)
      ENDIF
C
C  Now form the Projection Matrix
C
      CALL ZeroIT(P,NAT3*NAT3)
C
C  decide what to project
C
      If(IProj.EQ.1) NumP = 3         ! just project translations
      If(IProj.EQ.2) NumP = 6         ! translations and rotations
c
      DO 30 K=1,NumP
      DO 30 J=1,NAT3
      DO 30 I=1,NAT3
      P(I,J) = P(I,J) - TRVec(I,K)*TRVec(J,K)
 30   CONTINUE
      DO 40 I=1,NAT3
      P(I,I) = One + P(I,I)
 40   CONTINUE
c
      IF(IPRNT.GT.5) THEN
       WRITE(IOut,1200)
       CALL PrntMat(NAT3,NAT3,NAT3,P)
      ENDIF
C
C  Project out the translations/rotations from Hessian
C     HESS = P * HESS * P(t)
C
      CALL DGemm('N',    'N',    NAT3,   NAT3,   NAT3,
     $            One,    HESS,  NAT3,   P,      NAT3,
     $            Zero,   TRVec, NAT3)
      CALL DGemm('N',    'N',    NAT3,   NAT3,   NAT3,
     $            One,    P,     NAT3,   TRVec,  NAT3,
     $            Zero,   HESS,  NAT3)
cc      CALL MatAB(NAT3,NAT3,NAT3,HESS,P,TRVec,0)
cc      CALL MatAB(NAT3,NAT3,NAT3,P,TRVec,HESS,0)
c
      If(IPRNT.GT.1) Then
       If(IProj.EQ.1) WRITE(IOut,1300)
       If(IProj.EQ.2) WRITE(IOut,1400)
      EndIf
C
C  Restore original coordinates
C
      DO 50 I=1,NAtoms
      XC(1,I) = XC(1,I) + CX
      XC(2,I) = XC(2,I) + CY
      XC(3,I) = XC(3,I) + CZ
 50   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize inertia tensor')
 1100 FORMAT(/,' Vectors for Translations and Rotations')
 1200 FORMAT(/,' The Projection Matrix')
 1300 FORMAT(/,' Translations Projected Out of Hessian')
 1400 FORMAT(/,' Translations and Rotations Projected Out of Hessian')
c
      END

*Deck formtr
      SUBROUTINE FormTR(NAtoms,XC,T,V)
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8 XC(3,NAtoms),T(9),V(3,NAtoms,6)
C
      PARAMETER (One=1.0d0,TollZERO=1.0d-8)
C
C
C  This routine generates vectors corresponding to translations
C  and infinitesimal rotations given the coordinates (in centre
C  of mass frame) and the eigenvectors of the inertia tensor
C
C
      NAT3 = 3*NAtoms
c
      DO 10 I=1,NAtoms
      X = XC(1,I)
      Y = XC(2,I)
      Z = XC(3,I)
      CX = X*T(1) + Y*T(2) + Z*T(3)
      CY = X*T(4) + Y*T(5) + Z*T(6)
      CZ = X*T(7) + Y*T(8) + Z*T(9)
      V(1,I,1) = One
      V(2,I,2) = One
      V(3,I,3) = One
      V(1,I,4) = CY*T(7) - CZ*T(4)
      V(2,I,4) = CY*T(8) - CZ*T(5)
      V(3,I,4) = CY*T(9) - CZ*T(6)
      V(1,I,5) = CZ*T(1) - CX*T(7)
      V(2,I,5) = CZ*T(2) - CX*T(8)
      V(3,I,5) = CZ*T(3) - CX*T(9)
      V(1,I,6) = CX*T(4) - CY*T(1)
      V(2,I,6) = CX*T(5) - CY*T(2)
      V(3,I,6) = CX*T(6) - CY*T(3)
 10   CONTINUE
c
      DO 20 I=1,6
      skal = SProd(NAT3,V(1,1,I),V(1,1,I))
      IF(skal.GT.TollZERO) THEN
       skal = One/SQRT(skal)
       CALL VScal(NAT3,skal,V(1,1,I))
      ENDIF
 20   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE PutDUM(NCycle, NAtoms, NDum,   NCons,  XC,
     $                  MaxL,   IFunc,  ICTYP,  RCON,   IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates dummy atom coordinates from the coordinates of
C  the real atoms in its defining list
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C  NAtoms  -  number of real atoms
C  NDum    -  number of dummy atoms
C  NCons   -  number of constraints
C  XC      -  Cartesian coordinates for all real atoms
C             on exit coordinates of dummy atoms are appended
C  MaxL    -  maximum number of real atoms involved in
C             definition of dummy atom position + 1
C  IFunc   -  list for each dummy atom of the atom type
C             and the real atoms defining its position
C  ICTYP    -  constraint type
C               1 - fixed distance
C               2 - fixed bond angle
C               3 - fixed dihedral angle
C  RCON     -  constraint values
C  IC       -  atoms involved in constraint
C                IC1-IC2           distance constraint
C                IC1-IC2-IC3       angle constraint
C                IC1-IC2-IC3-IC4   dihedral constraint
C
C  All dummy atoms are defined with reference to a list of
C  real atoms and dummy atom coordinates will be generated
C  from the coordinates of the real atoms in its defining
C  list. There are three types of dummy atom:
C
C    1.  Those positioned at the arithmetic mean of the
C        up to MaxL-1 real atoms in the defining list
C    2.  Those positioned a unit distance along the normal to
C        a plane defined by 3 atoms, centred on the middle
C        atom of the 3 (for near-limiting dihedral constraints)
C    3.  Those positioned a unit distance along the bisector
C        of a given angle (for near-limiting angle constraints)
C
C  On the first cycle ALL dummy atom coordinates are generated.
C  On subsequent cycles only type 1 coordinates are regenerated;
C  type 2 & 3 dummies retain their previous coordinates except for
C  any repositioning due to movement of the system centre of mass
C
C
      REAL*8 XC(3,NAtoms+NDum),RCON(NCons),R1(3),R2(3),R3(3)
      DIMENSION IFunc(NDum,MaxL),ICTYP(NCons),IC(4,NCons)
      LOGICAL Switch
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      DO 30 I=1,NDum
c
      ITYP = IFunc(I,1)
c
      IF(ITYP.EQ.1) THEN
cc
       X = Zero
       Y = Zero
       Z = Zero
       DO 10 J=2,MaxL
       IAtm = IFunc(I,J)
       If(IAtm.EQ.0) GO TO 20
       X = X + XC(1,IAtm)
       Y = Y + XC(2,IAtm)
       Z = Z + XC(3,IAtm)
 10    CONTINUE
 20    CONTINUE
       J = J-2
       IDum = NAtoms+I
       XC(1,IDum) = X/DFloat(J)
       XC(2,IDum) = Y/DFloat(J)
       XC(3,IDum) = Z/DFloat(J)
cc
      ELSE IF(ITYP.EQ.2.AND.NCycle.EQ.1) THEN
cc
       IAtm = IFunc(I,2)
       JAtm = IFunc(I,3)
       KAtm = IFunc(I,4)
c
       XIJ = XC(1,IAtm) - XC(1,JAtm)
       YIJ = XC(2,IAtm) - XC(2,JAtm)
       ZIJ = XC(3,IAtm) - XC(3,JAtm)
       XJK = XC(1,JAtm) - XC(1,KAtm)
       YJK = XC(2,JAtm) - XC(2,KAtm)
       ZJK = XC(3,JAtm) - XC(3,KAtm)
       R1(1) = -XIJ
       R1(2) = -YIJ
       R1(3) = -ZIJ
       R2(1) = -XJK
       R2(2) = -YJK
       R2(3) = -ZJK
       CALL Normal(R1,R2,R3)
C
C  put the dummy atom along R3
C
       IDum = NAtoms+I
       XC(1,IDum) = XC(1,JAtm) + R3(1)
       XC(2,IDum) = XC(2,JAtm) + R3(2)
       XC(3,IDum) = XC(3,JAtm) + R3(3)
C
C  now check if dummy atom positioned correctly
C  (may need to be in opposite direction)
C
       CALL ChkDUM(NAtoms, NCons,  IDum,   XC,     ICTYP,
     $             RCON,   IC,     Switch)
c
       IF(Switch) THEN
        XC(1,IDum) = XC(1,JAtm) - R3(1)
        XC(2,IDum) = XC(2,JAtm) - R3(2)
        XC(3,IDum) = XC(3,JAtm) - R3(3)
       ENDIF
cc
      ELSE IF(ITYP.EQ.3.AND.NCycle.EQ.1) THEN
cc
       IAtm = IFunc(I,2)
       JAtm = IFunc(I,3)
       KAtm = IFunc(I,4)
c
       XIJ = XC(1,IAtm) - XC(1,JAtm)
       YIJ = XC(2,IAtm) - XC(2,JAtm)
       ZIJ = XC(3,IAtm) - XC(3,JAtm)
       XJK = XC(1,JAtm) - XC(1,KAtm)
       YJK = XC(2,JAtm) - XC(2,KAtm)
       ZJK = XC(3,JAtm) - XC(3,KAtm)
       R1(1) =  XIJ
       R1(2) =  YIJ
       R1(3) =  ZIJ
       R2(1) = -XJK
       R2(2) = -YJK
       R2(3) = -ZJK
       WNorm = SQRT(SProd(3,R1,R1))
       WNorm = One/WNorm
       CALL VScal(3,WNorm,R1)
       WNorm = SQRT(SProd(3,R2,R2))
       WNorm = One/WNorm
       CALL VScal(3,WNorm,R2)
C
C  calculate bisector of R1 and R2
C
       R3(1) = R1(1) + R2(1)
       R3(2) = R1(2) + R2(2)
       R3(3) = R1(3) + R2(3)
       WNorm = SQRT(SProd(3,R3,R3))
       WNorm = One/WNorm
       CALL VScal(3,WNorm,R3)
C
C  put the dummy atom along R3
C
       IDum = NAtoms+I
       XC(1,IDum) = XC(1,JAtm) + R3(1)
       XC(2,IDum) = XC(2,JAtm) + R3(2)
       XC(3,IDum) = XC(3,JAtm) + R3(3)
cc
      ENDIF
c
 30   CONTINUE
C
      RETURN
      END
