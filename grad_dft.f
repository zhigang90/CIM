      SUBROUTINE BasATOM(NAtoms,NShell,ILST,NBAtm)
      IMPLICIT INTEGER(A-Z)
C
C  determines the number of basis functions per atom
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NShell  -  total number of contracted shells
C  ILST    -  array containing data pertaining to each shell
C              (1,J) - shell type of Jth shell
C              (2,J) - # primitive shells in Jth shell
C              (3,J) - atom associated with Jth shell
C              (4,J) - # general contractions
C  NBAtm   -  on exit number of basis functions per atom
C
C
      DIMENSION ILST(4,NShell),NBAtm(NAtoms)
      Dimension ISIZE(10)
C
c                   S    P    L    D   D6    F  F10  G15  H21  I28
      Data ISIZE /  1,   3,   4,   5,   6,   7,  10,  15,  21,  28 /
C
C
      CALL IZeroIT(NBAtm,NAtoms)
C
C  loop over entries in ILST array
C
      DO 10 J=1,NShell
      IAtm = ILST(3,J)
      IType = ILST(1,J)
      ngr = ILST(4,J)
      nbas = ISIZE(IType)*(ngr+1)
      NBAtm(IAtm) = NBAtm(IAtm) + nbas
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE FormDEN(NBas,NOcc,CMO,DEN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  form the closed-shell density matrix (lower triangle)
C  from the transposed MOs
C
      DIMENSION CMO(NOcc,NBas),DEN(NBas*(NBas+1)/2)
      Parameter (Zero=0.0d0)
C
      IJ = 0
      DO 20 I=1,NBas
      DO 20 J=1,I
      IJ = IJ+1
      Val = Zero
      DO 10 K=1,NOcc
      Val = Val + CMO(K,I)*CMO(K,J)
 10   CONTINUE
      DEN(IJ) = Val + Val
 20   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE FormDENU(NBas,NAlpha,NBeta,CMOA,CMOB,DA,DB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  form the open-shell density matrices (lower triangle)
C  from the transposed MOs
C
      DIMENSION CMOA(NAlpha,NBas),CMOB(NBeta,NBas),DA(*),DB(*)
      Parameter (Zero=0.0d0)
C
      IJ = 0
      DO 40 I=1,NBas
      DO 40 J=1,I
      IJ = IJ+1
c -- alpha density
      Val = Zero
      DO 30 K=1,NAlpha
      Val = Val + CMOA(K,I)*CMOA(K,J)
 30   CONTINUE
      DA(IJ) = Val
c -- beta density
      Val = Zero
      DO 31 K=1,NBeta
      Val = Val + CMOB(K,I)*CMOB(K,J)
 31   CONTINUE
      DB(IJ) = Val
 40   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE GradMu(X,Y,Z,XB,YB,ZB,XA,YA,ZA,rAB,xaAB,xmuAB,
     $                  GradX,GradY,GradZ)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculates (Eq.B10) from Johnson, Gill & Pople
C  modified due to size-adjustment factor
C     Grad = Grad{xmuAB} * (1 - 2*xaAB*xmuAB)
c
c   Grad{A}
C
C  ARGUMENTS
C
C  X       -  x coordinate of grid point
C  Y       -  y coordinate of grid point
C  Z       -  z coordinate of grid point
C  XA      -  x coordinate of atom A
C  YA      -  y coordinate of atom A
C  ZA      -  z coordinate of atom A
C  XB      -  x coordinate of atom B
C  YB      -  y coordinate of atom B
C  ZB      -  z coordinate of atom B
C  rAB     -  inverse distance between atoms A and B
C  xaAB    -  atomic size adjustment factor
C  xmuAB   -  adjusted hyperbolic coordinate
C  GradX   -  on exit contains x gradient component
C  GradY   -  on exit contains y gradient component
C  GradZ   -  on exit contains z gradient component
C
C
      PARAMETER (One=1.0d0,Two=2.0d0)
C
C  unit vector of grid point to atom A
C
      XX = (X-XA)
      YY = (Y-YA)
      ZZ = (Z-ZA)
      rAI = SQRT(XX*XX + YY*YY + ZZ*ZZ)
      rAII = -rAB/rAI
      uAX = XX*rAII
      uAY = YY*rAII
      uAZ = ZZ*rAII
C
C  unit vector in direction of B to A
C
      rAII = -rAB*rAB*xmuAB
      uBX = (XB-XA)*rAII
      uBY = (YB-YA)*rAII
      uBZ = (ZB-ZA)*rAII
C
C  now put it together
C
      coef = (One - Two*xaAB*xmuAB)
      GradX = (uAX-uBX)*coef
      GradY = (uAY-uBY)*coef
      GradZ = (uAZ-uBZ)*coef
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GradWT(NPP,    XGRID,  WGHT,   INDX,   IATOM,
     $                  EVec,   NAtoms, XC,     AIJ,    rrij,
     $                  RDist,  P,      T,      xmuAB,  GWT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates derivative of weight at contributing grid points for
C  each atom and stores contributions in GWT array.
C
C  See:  Appendix B of Johnson, Gill and Pople
C        J.Chem.Phys.  98 (1993) 5612
C
C  ARGUMENTS
C
C  NPP     -  number of contributing (non-zero) grid points this batch
C  XGRID   -  X,Y,Z grid points
C  WGHT    -  grid quadrature weights
C  INDX    -  index into contributing grid points
C  IATOM   -  which atom grid point formerly associated with
C  EVec    -  Exchange-Correlation energy at grid points
C  NAtoms  -  total number of atoms in system
C  XC      -  Cartesian coordinates of all atoms
C  AIJ     -  Becke atomic size-adjustment factors
C  rrij    -  inverse of interatomic distances
C  RDist   -  scratch array for interatomic distances
C  P       -  scratch array for cell functions
C  T       -  scratch array for auxillary functions
C  xmuAB   -  scratch array for gradients of cell functions
C  GWT     -  on exit contains cumulative weight derivatives
C
C
      DIMENSION XGRID(3,*),WGHT(*),INDX(NPP),EVec(NPP)
      REAL*8 XC(3,NAtoms),GWT(3,NAtoms)
      REAL*8 AIJ(NAtoms,NAtoms),rrij(natoms,natoms)
      real*8 RDist(natoms),P(natoms),T(natoms,natoms),
     $       xmuAB(natoms,natoms)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0,ThreeHalf=3.0d0/2.0d0)
      PARAMETER (CON1=-27.0d0/16.0d0)
C
      data thrsh/1.0d-14/
C
C
C  loop over the grid points
C
      DO 90 IP=1,NPP
      IPP = INDX(IP)
c...
c...  we compute the prefactor and test it.
c...  If it is smaller than the threshold, we can skip
c...  this grid point
c...
      Wt = evec(IPP)*wght(IPP)
      if(abs(wt).lt.thrsh) goto 90
c...
      X = XGRID(1,IPP)
      Y = XGRID(2,IPP)
      Z = XGRID(3,IPP)
C
C  form the distance vector
C  (distances of current grid point to all atoms)
C
      DO 10 IAtm=1,NAtoms
      XX = X - XC(1,IAtm)
      YY = Y - XC(2,IAtm)
      ZZ = Z - XC(3,IAtm)
      RDist(IAtm) = SQRT(XX*XX + YY*YY + ZZ*ZZ)
 10   CONTINUE
C
C  form the auxiliary function array T  (Eq.B9)
C
      ZF = Zero
      DO 30 IAtm=1,NAtoms
      PSum = One
      DO 20 JAtm=1,NAtoms
      If(JAtm.EQ.IAtm) GO TO 20
      xmu = (RDist(IAtm)-RDist(JAtm))*rrij(IAtm,JAtm)
      xmum = xmu + AIJ(IAtm,JAtm)*(One-xmu*xmu)
      xmuAB(IAtm,JAtm) = xmu
c
      p1 = ThreeHalf*xmum - Half*(xmum**3)
      p2 = ThreeHalf*p1   - Half*(p1**3)
      p3 = ThreeHalf*p2   - Half*(p2**3)
      s  = Half*(One-p3)
c
      IF(s.LT.thrsh) THEN
       T(IAtm,JAtm) = Zero
      ELSE
       T(IAtm,JAtm) = CON1*(One-p2*p2)*(One-p1*p1)*(One-xmum*xmum)/s
      ENDIF
c
      PSum = PSum*s
c
 20   CONTINUE
      P(IAtm) = PSum
      ZF = ZF + PSum
 30   CONTINUE
C
      IF(ZF.LT.thrsh) THEN
       ZFF = Zero
      ELSE
       ZFF = One/ZF
      ENDIF
c...
c...  loop over the atom respect to which the gradient has
c...  to be computed (that is, all atoms except IATOM)
c...
      DO 60 IAtm=1,NATOMS
      If(IAtm.EQ.IATOM) GO TO 60
c...
c...  Eq. (B7) can be developed and divided into two contributions:
c...
c... D_I(W(A)) = -P(A)/Z [(1 + P(A)/Z ) T(A,I) + P(I)/Z T(I,A)] D_I(mu(I,A))
c...           - P(A)/Z**2 SUM_J [P(I) T(I,J) - P(J) T(J,I)] D_I(mu(I,J))
c...
c...   where the summation is over all atoms J different from A and I.
c...   D_I(mu(I,J)) is the gradient of mu(I,J) with respect atom I (this
c...   is computed by subroutine GradMu).
c...
c...   First part:
c...
c...  we compute the gradient of mu(I,A)
c...
      CALL GradMu(X,Y,Z,XC(1,IATOM),XC(2,IATOM),XC(3,IATOM),
     $            XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),rrij(IATOM,IAtm),
     $            AIJ(IATOM,IAtm),xmuAB(IATOM,IAtm),
     $            GradX,GradY,GradZ)
c...
c...  now we compute the corresponding coefficient and sum into GWT.
c...  here one term P(A)/Z  is omitted, as P(A)/Z = omega(A) has already
c...  been incorporated into the prefactor WT
c...
      coef= -wt*
     +   ((1-p(iatom)*zff)*t(iatom,iatm)+p(iatm)*t(iatm,iatom)*zff)
      GWT(1,IAtm) = GWT(1,IAtm) + coef*GradX
      GWT(2,IAtm) = GWT(2,IAtm) + coef*GradY
      GWT(3,IAtm) = GWT(3,IAtm) + coef*GradZ
c...
c...  second part:
c...
c...  loop over atom J, skipping J=A and J=I
c...
      DO 50 JAtm=1,NATOMS
c
      IF(JAtm.EQ.IATOM) goto 50
      IF(JAtm.EQ.IAtm) goto 50
c...
c...  we compute the gradient of mu(I,J)
c...

       CALL GradMu(X,Y,Z,XC(1,JAtm),XC(2,JAtm),XC(3,JAtm),
     $             XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),rrij(JAtm,IAtm),
     $             AIJ(JAtm,IAtm),xmuAB(JAtm,IAtm),
     $             GradX,GradY,GradZ)
c...
c...  now we compute the corresponding coefficient and sum into GWT.
c...  here one term P(A)/Z  is omitted, as P(A)/Z = omega(A) has already
c...  been incorporated into the prefactor WT
c...
       coef =-wt*
     +       (p(iatm)*t(iatm,jatm)-p(jatm)*t(jatm,iatm))*zff
       GWT(1,IAtm) = GWT(1,IAtm) + coef*GradX
       GWT(2,IAtm) = GWT(2,IAtm) + coef*GradY
       GWT(3,IAtm) = GWT(3,IAtm) + coef*GradZ
cc
 50   CONTINUE
 60   CONTINUE
c -- end loop over grid points
 90   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE MakeGXC(dft,    NPP,    NBas,   nbf,    thrsh,
     $                   WGHT,   Pot,    PotX,   VAO,    VAOX,
     $                   VAOXX,  INB,    INDX,   GX,     GY,     GZ)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into Exchange-Correlation gradient matrices
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All methods include Slater exchange)
C              1 - 3 - local correlation
C             >3 - nonlocal - need AO 2nd derivatives
C  NPP     -  number of contributing (non-zero) grid points this batch
C  NBas    -  total number of basis functions
C  nbf     -  indexing array to location of "non-zero" basis functions
C  thrsh   -  threshold for neglect of contribution
C  WGHT    -  grid quadrature weights
C  Pot     -  potential
C  PotX    -  potential derivatives
C  VAO     -  "non-zero" basis function values at grid points
C  VAOX    -  basis function 1st derivatives at grid points
C  VAOXX   -  basis function 2nd derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAO
C  INDX    -  index into contributing columns of VAO
C
C  on exit
C
C  GX      -  x DFT derivative matrix (closed shell)
C  GY      -  y DFT derivative matrix (closed shell)
C  GZ      -  z DFT derivative matrix (closed shell)
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VAOXX(6,*),
     $          Pot(NPP),PotX(3,NPP)
      DIMENSION GX(NBas,NBas),GY(NBas,NBas),GZ(NBas,NBas)
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
C
      PARAMETER (Half=0.5d0)
C
C
      IF(dft.LE.3) THEN
C
C  local only
C
        DO 30 IP=1,NPP
        IPP = INDX(IP)
        VXC = Pot(IP)*WGHT(IPP)
        DO 20 J=nbf(IPP)+1,nbf(IPP+1)
        JJ = INB(J)
        Val = VAO(J)*VXC
        IF(Abs(Val).GT.thrsh) THEN
          DO 10 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          GX(II,JJ) = GX(II,JJ) + VAOX(1,I)*Val
          GY(II,JJ) = GY(II,JJ) + VAOX(2,I)*Val
          GZ(II,JJ) = GZ(II,JJ) + VAOX(3,I)*Val
 10       CONTINUE
        ENDIF
 20     CONTINUE
 30     CONTINUE
cc
      ELSE
C
C  non-local
C
        DO 70 IP=1,NPP
        IPP = INDX(IP)
        Wt = WGHT(IPP)
        VXC = Pot(IP)*Wt
        VAx = PotX(1,IP)*Wt
        VAy = PotX(2,IP)*Wt
        VAz = PotX(3,IP)*Wt
        DO 60 J=nbf(IPP)+1,nbf(IPP+1)
        JJ = INB(J)
        ValJ = VAO(J)
        ValJA = ValJ*VXC +
     $    (VAOX(1,J)*VAx + VAOX(2,J)*VAy + VAOX(3,J)*VAz)
        IF(Abs(ValJA).GT.thrsh) THEN
          DO 50 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          GX(II,JJ) = GX(II,JJ) + ValJA*VAOX(1,I) + ValJ*
     $      (VAOXX(1,I)*VAx + VAOXX(2,I)*VAy + VAOXX(4,I)*VAz)
          GY(II,JJ) = GY(II,JJ) + ValJA*VAOX(2,I) + ValJ*
     $      (VAOXX(2,I)*VAx + VAOXX(3,I)*VAy + VAOXX(5,I)*VAz)
          GZ(II,JJ) = GZ(II,JJ) + ValJA*VAOX(3,I) + ValJ*
     $      (VAOXX(4,I)*VAx + VAOXX(5,I)*VAy + VAOXX(6,I)*VAz)
 50       CONTINUE
        ENDIF
 60     CONTINUE
 70     CONTINUE
      ENDIF
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE MakeGXU(dft,    NPP,    NBas,   nbf,    thrsh,
     $                   WGHT,   PotA,   PotB,   PotAX,  PotBX,
     $                   VAO,    VAOX,   VAOXX,  INB,    INDX,
     $                   GXA,    GYA,    GZA,    GXB,    GYB,   GZB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into Exchange-Correlation gradient matrices
C  ** UNRESTRICTED OPEN SHELL **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All methods include Slater exchange)
C              1 - 3 - local correlation
C             >3 - nonlocal - need AO 2nd derivatives
C  NPP     -  number of contributing (non-zero) grid points this batch
C  NBas    -  total number of basis functions
C  nbf     -  number of "non-zero" basis functions
C  thrsh   -  threshold for neglect of contribution
C  WGHT    -  grid quadrature weights
C  PotA    -  alpha potential
C  PotB    -  beta potential
C  PotAX   -  alpha potential derivatives
C  PotBX   -  beta potential derivatives
C  VAO     -  "non-zero" basis function values at grid points
C  VAOX    -  basis function 1st derivatives at grid points
C  VAOXX   -  basis function 2nd derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAO
C  INDX    -  index into contributing columns of VAO
C
C  on exit
C
C  GXA     -  x DFT derivative matrix (alpha spin)
C  GYA     -  y DFT derivative matrix (alpha spin)
C  GZA     -  z DFT derivative matrix (alpha spin)
C  GXB     -  x DFT derivative matrix (beta spin)
C  GYB     -  y DFT derivative matrix (beta spin)
C  GZB     -  z DFT derivative matrix (beta spin)
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VAOXX(6,*),
     $          PotA(NPP),PotB(NPP),PotAX(3,NPP),PotBX(3,NPP)
      DIMENSION GXA(NBas,NBas),GYA(NBas,NBas),GZA(NBas,NBas),
     $          GXB(NBas,NBas),GYB(NBas,NBas),GZB(NBas,NBas)
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
C
      PARAMETER (Half=0.5d0)
C
C
      IF(dft.LE.3) THEN
C
C  local only
C
        DO 30 IP=1,NPP
        IPP = INDX(IP)
        VXCA = PotA(IP)*WGHT(IPP)
        VXCB = PotB(IP)*WGHT(IPP)
        DO 20 J=nbf(IPP)+1,nbf(IPP+1)
        JJ = INB(J)
        VA = VAO(J)*VXCA
        VB = VAO(J)*VXCB
        IF(Abs(VA).GT.thrsh.OR.Abs(VB).GT.thrsh) THEN
          DO 10 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          GXA(II,JJ) = GXA(II,JJ) + VAOX(1,I)*VA
          GYA(II,JJ) = GYA(II,JJ) + VAOX(2,I)*VA
          GZA(II,JJ) = GZA(II,JJ) + VAOX(3,I)*VA
          GXB(II,JJ) = GXB(II,JJ) + VAOX(1,I)*VB
          GYB(II,JJ) = GYB(II,JJ) + VAOX(2,I)*VB
          GZB(II,JJ) = GZB(II,JJ) + VAOX(3,I)*VB
 10       CONTINUE
        ENDIF
 20     CONTINUE
 30     CONTINUE
cc
      ELSE
C
C  non-local
C
        DO 70 IP=1,NPP
        IPP = INDX(IP)
        Wt = WGHT(IPP)
        VXCA = PotA(IP)*Wt
        VXCB = PotB(IP)*Wt
        VAx = PotAX(1,IP)*Wt
        VAy = PotAX(2,IP)*Wt
        VAz = PotAX(3,IP)*Wt
        VBx = PotBX(1,IP)*Wt
        VBy = PotBX(2,IP)*Wt
        VBz = PotBX(3,IP)*Wt
        DO 60 J=nbf(IPP)+1,nbf(IPP+1)
        JJ = INB(J)
        ValJ = VAO(J)
        ValJA = ValJ*VXCA +
     $    (VAOX(1,J)*VAx + VAOX(2,J)*VAy + VAOX(3,J)*VAz)
        ValJB = ValJ*VXCB +
     $    (VAOX(1,J)*VBx + VAOX(2,J)*VBy + VAOX(3,J)*VBz)
        IF(Abs(ValJA).GT.thrsh.OR.Abs(ValJB).GT.thrsh) THEN
          DO 50 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          GXA(II,JJ) = GXA(II,JJ) + ValJA*VAOX(1,I) + ValJ*
     $      (VAOXX(1,I)*VAx + VAOXX(2,I)*VAy + VAOXX(4,I)*VAz)
          GYA(II,JJ) = GYA(II,JJ) + ValJA*VAOX(2,I) + ValJ*
     $      (VAOXX(2,I)*VAx + VAOXX(3,I)*VAy + VAOXX(5,I)*VAz)
          GZA(II,JJ) = GZA(II,JJ) + ValJA*VAOX(3,I) + ValJ*
     $      (VAOXX(4,I)*VAx + VAOXX(5,I)*VAy + VAOXX(6,I)*VAz)
          GXB(II,JJ) = GXB(II,JJ) + ValJB*VAOX(1,I) + ValJ*
     $      (VAOXX(1,I)*VBx + VAOXX(2,I)*VBy + VAOXX(4,I)*VBz)
          GYB(II,JJ) = GYB(II,JJ) + ValJB*VAOX(2,I) + ValJ*
     $      (VAOXX(2,I)*VBx + VAOXX(3,I)*VBy + VAOXX(5,I)*VBz)
          GZB(II,JJ) = GZB(II,JJ) + ValJB*VAOX(3,I) + ValJ*
     $      (VAOXX(4,I)*VBx + VAOXX(5,I)*VBy + VAOXX(6,I)*VBz)
 50       CONTINUE
        ENDIF
 60     CONTINUE
 70     CONTINUE
      ENDIF
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE PartGrad(NAtoms, NBas,   NBAtm,  IATOM,  DA,
     $                    GXA,    GYA,    GZA,    GWT,    GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Form the partial Cartesian gradient from all contributions due
C  to motion of atom IATOM (excluding the gradient for IATOM itself).
C  Done this way so that IATOM's contribution can be obtained by
C  invoking translational invarience
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NBas    -  total number of basis functions
C  NBAtm   -  number of basis functions per atom
C  IATOM   -  atom currently being "moved"
C  DA      -  closed shell density matrix
C  GXA     -  x DFT derivative matrix (closed shell)
C  GYA     -  y DFT derivative matrix (closed shell
C  GZA     -  z DFT derivative matrix (closed shell
C  GWT     -  grid weight derivatives
C  GC      -  accumulates total DFT contribution to Cartesian gradient
C
C
      DIMENSION GXA(NBas,NBas),GYA(NBas,NBas),GZA(NBas,NBas),
     $          DA(NBas*(NBas+1)/2),NBAtm(NAtoms),
     $          GWT(3,NAtoms),GC(3,NAtoms)
C
      PARAMETER (Zero=0.0d0)
C
C
C  Form partial Cartesian gradient from DFT gradient matrices
C  and accumulate with weight derivatives
C
      write(6,*) ' Weight derivatives are:'
      do i=1,natoms
      write(6,*) gwt(1,i),gwt(2,i),gwt(3,i)
      enddo
cc
      I1 = 0
      DO 50 IAtm=1,NAtoms
      I1 = I1+1
      I2 = I1+NBAtm(IAtm)-1
      IF(IAtm.NE.IATOM) THEN
        GradX = Zero
        GradY = Zero
        GradZ = Zero
        DO 49 J=1,NBas
        JJ = (J*(J-1))/2
        DO 48 I=I1,I2
        If(I.GT.J) Then
         IJ = (I*(I-1))/2 + J
        Else
         IJ = JJ + I
        EndIf
        GradX = GradX + DA(IJ)*GXA(I,J)
        GradY = GradY + DA(IJ)*GYA(I,J)
        GradZ = GradZ + DA(IJ)*GZA(I,J)
 48     CONTINUE
 49     CONTINUE
        GWT(1,IAtm) = GWT(1,IAtm) - GradX - GradX
        GWT(2,IAtm) = GWT(2,IAtm) - GradY - GradY
        GWT(3,IAtm) = GWT(3,IAtm) - GradZ - GradZ
      ENDIF
      I1 = I2
 50   CONTINUE
      write(6,*) ' DFT derivatives for this atom are:'
      do i=1,natoms
      write(6,*) gwt(1,i),gwt(2,i),gwt(3,i)
      enddo
C
C  ...................................................................
C    Get Derivative for atom IATOM via Translational Invarience
C  ...................................................................
C
      GX = Zero
      GY = Zero
      GZ = Zero
      DO 90 IAtm=1,NAtoms
      GX = GX + GWT(1,IAtm)
      GY = GY + GWT(2,IAtm)
      GZ = GZ + GWT(3,IAtm)
 90   CONTINUE
c
      GWT(1,IATOM) = - GX
      GWT(2,IATOM) = - GY
      GWT(3,IATOM) = - GZ
C
C  Accumulate into full DFT Cartesian gradient
C
      DO 95 IAtm=1,NAtoms
      GC(1,IAtm) = GC(1,IAtm) + GWT(1,IAtm)
      GC(2,IAtm) = GC(2,IAtm) + GWT(2,IAtm)
      GC(3,IAtm) = GC(3,IAtm) + GWT(3,IAtm)
 95   CONTINUE
      write(6,*) ' Cumulative Cartesian Gradient is:'
      do i=1,natoms
      write(6,*) gc(1,i),gc(2,i),gc(3,i)
      enddo
C
C  Zero out DFT gradient arrays
C
      CALL ZeroIT(GWT,3*NAtoms)
      CALL ZeroIT(GXA,NBas*NBas)
      CALL ZeroIT(GYA,NBas*NBas)
      CALL ZeroIT(GZA,NBas*NBas)
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE PartGradU(NAtoms, NBas,   NBAtm,  IATOM,  DA,
     $                     DB,     GXA,    GYA,    GZA,    GXB,
     $                     GYB,    GZB,    GWT,    GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Form the partial Cartesian gradient from all contributions due
C  to motion of atom IATOM (excluding the gradient for IATOM itself).
C  Done this way so that IATOM's contribution can be obtained by
C  invoking translational invarience
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NBas    -  total number of basis functions
C  NBAtm   -  number of basis functions per atom
C  IATOM   -  atom currently being "moved"
C  DA      -  alpha density matrix
C  DB      -  beta density matrix
C  GXA     -  alpha x DFT derivative matrix
C  GYA     -  alpha y DFT derivative matrix
C  GZA     -  alpha z DFT derivative matrix
C  GXB     -  beta x DFT derivative matrix
C  GYB     -  beta y DFT derivative matrix
C  GZB     -  beta z DFT derivative matrix
C  GWT     -  grid weight derivatives
C  GC      -  accumulates total DFT contribution to Cartesian gradient
C
C
      DIMENSION GXA(NBas,NBas),GYA(NBas,NBas),GZA(NBas,NBas),
     $          GXB(NBas,NBas),GYB(NBas,NBas),GZB(NBas,NBas),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),NBAtm(NAtoms),
     $          GWT(3,NAtoms),GC(3,NAtoms)
C
      PARAMETER (Zero=0.0d0)
C
C
C  Form partial Cartesian gradient from DFT gradient matrices
C  and accumulate with weight derivatives
C
      write(6,*) ' Weight derivatives are:'
      do i=1,natoms
      write(6,*) gwt(1,i),gwt(2,i),gwt(3,i)
      enddo
cc
      I1 = 0
      DO 50 IAtm=1,NAtoms
      I1 = I1+1
      I2 = I1+NBAtm(IAtm)-1
      IF(IAtm.NE.IATOM) THEN
        GradX = Zero
        GradY = Zero
        GradZ = Zero
        DO 49 J=1,NBas
        JJ = (J*(J-1))/2
        DO 48 I=I1,I2
        If(I.GT.J) Then
         IJ = (I*(I-1))/2 + J
        Else
         IJ = JJ + I
        EndIf
        GradX = GradX + DA(IJ)*GXA(I,J) + DB(IJ)*GXB(I,J)
        GradY = GradY + DA(IJ)*GYA(I,J) + DB(IJ)*GYB(I,J)
        GradZ = GradZ + DA(IJ)*GZA(I,J) + DB(IJ)*GZB(I,J)
 48     CONTINUE
 49     CONTINUE
        GWT(1,IAtm) = GWT(1,IAtm) - GradX - GradX
        GWT(2,IAtm) = GWT(2,IAtm) - GradY - GradY
        GWT(3,IAtm) = GWT(3,IAtm) - GradZ - GradZ
      ENDIF
      I1 = I2
 50   CONTINUE
      write(6,*) ' DFT derivatives for this atom are:'
      do i=1,natoms
      write(6,*) gwt(1,i),gwt(2,i),gwt(3,i)
      enddo
C
C  ...................................................................
C    Get Derivative for atom IATOM via Translational Invarience
C  ...................................................................
C
      GX = Zero
      GY = Zero
      GZ = Zero
      DO 90 IAtm=1,NAtoms
      GX = GX + GWT(1,IAtm)
      GY = GY + GWT(2,IAtm)
      GZ = GZ + GWT(3,IAtm)
 90   CONTINUE
c
      GWT(1,IATOM) = - GX
      GWT(2,IATOM) = - GY
      GWT(3,IATOM) = - GZ
C
C  Accumulate into full DFT Cartesian gradient
C
      DO 95 IAtm=1,NAtoms
      GC(1,IAtm) = GC(1,IAtm) + GWT(1,IAtm)
      GC(2,IAtm) = GC(2,IAtm) + GWT(2,IAtm)
      GC(3,IAtm) = GC(3,IAtm) + GWT(3,IAtm)
 95   CONTINUE
      write(6,*) ' Cumulative Cartesian Gradient is:'
      do i=1,natoms
      write(6,*) gc(1,i),gc(2,i),gc(3,i)
      enddo
C
C  Zero out DFT gradient arrays
C
      CALL ZeroIT(GWT,3*NAtoms)
      CALL ZeroIT(GXA,NBas*NBas)
      CALL ZeroIT(GYA,NBas*NBas)
      CALL ZeroIT(GZA,NBas*NBas)
      CALL ZeroIT(GXB,NBas*NBas)
      CALL ZeroIT(GYB,NBas*NBas)
      CALL ZeroIT(GZB,NBas*NBas)
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE ReorderMO(NBas,NOcc,INDX,COld,CNew)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reorder MOs so basis functions are ordered per atom
C  as opposed to Wolinski special ordering
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  NOcc    -  number of occupied MOs
C  INDX    -  index array relating atom ordering to Wolinski ordering
C             (output from subroutine <SortBAS2>)
C  COld    -  Full set of MO coefficients in Wolinski order
C  CNew    -  on output contains transposed MOs in atom order
C
C
      DIMENSION INDX(NBas),COld(NBas,NBas),CNew(NOcc,NBas)
C
C
      DO 20 I=1,NBas
      II = INDX(I)
      DO 10 J=1,NOcc
      CNew(J,II) = COld(I,J)
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE ReorderIFP(NBas,NSym,INDX,IFP,IFPP)
      IMPLICIT INTEGER(A-Z)
C
C  Reorder IFP so basis functions are ordered per atom
C  as opposed to Wolinski special ordering
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  NSym    -  number of (Abelian) symmetry operations
C  INDX    -  index array relating Wolinski ordering to atom ordering
C             (output from subroutine <SortBAS2>)
C  IFP     -  original shell pair symmetry ordering
C  IFPP    -  on exit reordered shell pair symmetry ordering
C
C
      DIMENSION INDX(NBas),IFP(7,NBas),IFPP(7,NBas)
C
C
      DO 20 I=1,NBas
      II = INDX(I)
      DO 10 J=1,NSym
      KK = IAbs(IFP(J,I))
      KK = INDX(KK)
      IFPP(J,II) = ISIGN(KK,IFP(J,I))
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
