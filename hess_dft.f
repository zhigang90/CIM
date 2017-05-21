      SUBROUTINE SortBATM(NBas,NShell,INX,NBAtm)
      IMPLICIT INTEGER(A-Z)
C
C  forms the basis function/atom ownership list
C  i.e., which atom each basis function is associated with
C  this data is derived from the INX array
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  INX     -  standard INX array with basis function data
C  NBAtm   -  on exit, basis function/atom ownership list
C
C
      DIMENSION INX(12,NShell),NBAtm(NBas)
C
      IT = 0
      DO 10 I=1,NShell
      nb = INX(10,I) - INX(11,I)
      IAtom = INX(2,I)
      DO 9 J=1,nb
      IT = IT+1
      NBAtm(IT) = IAtom
  9   CONTINUE
 10   CONTINUE
C
C  IT should equal NBas (number of contracted basis functions)
C
      If(IT.NE.NBas) write(6,*) ' Oh, My God!  We are DEAD!!'
C
      RETURN
      END
c =================================================================
      SUBROUTINE DENHess(NP,     NBas,   thrsh,  DA,     nbf,
     $                   VAO,    VAOX,   VAOXX,  INB,    DenX,
     $                   DenXX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the gradient and Hessian density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  DA      -  closed-shell density matrix (lower triangle)
C  nbf     -  index array to indices of "non-zero" basis functions
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  "non-zero" basis function 1st derivatives at grid point
C  VAOXX   -  "non-zero" basis function 2nd derivatives at grid point
C  INB     -  indices corresponding to entries in VAO
C
C  on exit
C
C  DenX    -  density gradient per grid point
C  DenXX   -  density Hessian per grid point
C
C
      DIMENSION DA(*),VAO(*),VAOX(3,*),VAOXX(6,*),nbf(*),INB(*)
      DIMENSION DenX(3,*),DenXX(6,*)
C
      PARAMETER (Zero=0.0d0)
C
C
      DO 30 IP=1,NP
c
      DX = Zero
      DY = Zero
      DZ = Zero
      DXX = Zero
      DXY = Zero
      DYY = Zero
      DXZ = Zero
      DYZ = Zero
      DZZ = Zero
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      DO 19 J=nbf(IP)+1,nbf(IP+1)
      JJ = INB(J)
      If(II.GE.JJ) Then
        IJ = IT + JJ
      Else
        IJ = (JJ*(JJ-1))/2 + II
      EndIf
c -- gradient density
      DAIJ = DA(IJ)*VAO(J)
      DX = DX + DAIJ*VAOX(1,I)
      DY = DY + DAIJ*VAOX(2,I)
      DZ = DZ + DAIJ*VAOX(3,I)
c -- Hessian density
      DXX = DXX + DA(IJ)*(VAOXX(1,I)*VAO(J) + VAOX(1,I)*VAOX(1,J))
      DXY = DXY + DA(IJ)*(VAOXX(2,I)*VAO(J) + VAOX(2,I)*VAOX(1,J))
      DYY = DYY + DA(IJ)*(VAOXX(3,I)*VAO(J) + VAOX(2,I)*VAOX(2,J))
      DXZ = DXZ + DA(IJ)*(VAOXX(4,I)*VAO(J) + VAOX(3,I)*VAOX(1,J))
      DYZ = DYZ + DA(IJ)*(VAOXX(5,I)*VAO(J) + VAOX(3,I)*VAOX(2,J))
      DZZ = DZZ + DA(IJ)*(VAOXX(6,I)*VAO(J) + VAOX(3,I)*VAOX(3,J))
 19   CONTINUE
 20   CONTINUE
      DenX(1,IP) = DX + DX
      DenX(2,IP) = DY + DY
      DenX(3,IP) = DZ + DZ
      DenXX(1,IP) = DXX + DXX
      DenXX(2,IP) = DXY + DXY
      DenXX(3,IP) = DYY + DYY
      DenXX(4,IP) = DXZ + DXZ
      DenXX(5,IP) = DYZ + DYZ
      DenXX(6,IP) = DZZ + DZZ
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE DENCphf1(NP,     NAtoms, NBas,   thrsh,  nbf,
     $                    DEN1,   DM1,    VAO,    INB,    VM,
     $                    PotD,   Den1P)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the gradient and Hessian density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NAtoms  -  number of atoms
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  nbf     -  index array to indices of "non-zero" basis functions
C  DEN1    -  current solution to CPHF, lower triangle (closed-shell)
C  DM1     -  array containing maximum 1st-order density element per column
C  VAO     -  "non-zero" basis function values at grid point
C  INB     -  indices corresponding to entries in VAO
C  VM      -  array containing maximum magnitude AO per grid point
C  PotD    -  density potential derivative
C  Den1P   -  on exit contains 1st-order density per grid point
C
C
      DIMENSION DEN1(3,NAtoms,NBas*(NBas+1)/2),VAO(NBas*NP),VM(NP),
     $          nbf(*),INB(*),DM1(NBas),PotD(NP),Den1P(3,NAtoms,NP)
C
      PARAMETER (Two=2.0d0)
C
C
      DO 30 IP=1,NP
      VMx = VM(IP)
      PD = PotD(IP)
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = PD*VAO(I)
      DMx = DM1(II)
      If(Two*DMx*VMx*Abs(Val).LT.thrsh) GO TO 20
      DO 19 J=nbf(IP)+1,I-1
      JJ = INB(J)
      IJ = IT + JJ
c -- first-order density
      VAIJ = Two*Val*VAO(J)
      If(Abs(VAIJ)*DMx.GT.thrsh) Then
        DO IAtm=1,NAtoms
        Den1P(1,IAtm,IP) = Den1P(1,IAtm,IP) + DEN1(1,IAtm,IJ)*VAIJ
        Den1P(2,IAtm,IP) = Den1P(2,IAtm,IP) + DEN1(2,IAtm,IJ)*VAIJ
        Den1P(3,IAtm,IP) = Den1P(3,IAtm,IP) + DEN1(3,IAtm,IJ)*VAIJ
        EndDO
      EndIf
 19   CONTINUE
      Val2 = Val**2/PD
      If(DMx*Val2.GT.thrsh) Then
        IJ = IT+II
        DO IAtm=1,NAtoms
        Den1P(1,IAtm,IP) = Den1P(1,IAtm,IP) + DEN1(1,IAtm,IJ)*Val2
        Den1P(2,IAtm,IP) = Den1P(2,IAtm,IP) + DEN1(2,IAtm,IJ)*Val2
        Den1P(3,IAtm,IP) = Den1P(3,IAtm,IP) + DEN1(3,IAtm,IJ)*Val2
        EndDO
      EndIf
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE DENHess1(NP,     NAtoms, NBas,   NBAtm,  thrsh,
     $                    DA,     nbf,    VAO,    VAOX,   VAOXX,
     $                    INB,    VM,     DenX,   DenXX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the gradient and Hessian density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NAtoms  -  number of atoms
C  NBas    -  total number of basis functions
C  NBAtm   -  basis functions --> atom index
C  thrsh   -  threshold for neglect of contribution
C  DA      -  closed-shell density matrix (lower triangle)
C  nbf     -  index array to indices of "non-zero" basis functions
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  "non-zero" basis function 1st derivatives at grid point
C  VAOXX   -  "non-zero" basis function 2nd derivatives at grid point
C  INB     -  indices corresponding to entries in VAO
C  VM      -  array containing maximum magnitude AO per grid point
C
C  on exit
C
C  DenX    -  density gradient per atom perturbation
C  DenXX   -  density Hessian per atom perturbation
C
C
      DIMENSION DA(*),VAO(*),VAOX(3,*),VAOXX(6,*),nbf(*),INB(*)
      DIMENSION NBAtm(NBas),VM(NP),DenX(3,NAtoms,NP),
     $          DenXX(3,NAtoms,3,NAtoms,NP)
C
      PARAMETER (Two=2.0d0)
C
C
      DO 30 IP=1,NP
      VMx = VM(IP)
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      IAtm = NBAtm(II)
      DO 18 J=nbf(IP)+1,I
      JJ = INB(J)
      IJ = IT + JJ
      JAtm = NBAtm(JJ)
      DAIJ = Two*DA(IJ)*VAO(J)
      If(Abs(DAIJ)*VMx.LT.thrsh) GO TO 185
c -- gradient density
      DenX(1,IAtm,IP) = DenX(1,IAtm,IP) - DAIJ*VAOX(1,I)
      DenX(2,IAtm,IP) = DenX(2,IAtm,IP) - DAIJ*VAOX(2,I)
      DenX(3,IAtm,IP) = DenX(3,IAtm,IP) - DAIJ*VAOX(3,I)
c -- Hessian density
c -- (a) IAtm with IAtm
      DenXX(1,IAtm,1,IAtm,IP)=DenXX(1,IAtm,1,IAtm,IP)+DAIJ*VAOXX(1,I)
      xyc=DAIJ*VAOXX(2,I)
      DenXX(1,IAtm,2,IAtm,IP)=DenXX(1,IAtm,2,IAtm,IP)+xyc
      DenXX(2,IAtm,1,IAtm,IP)=DenXX(2,IAtm,1,IAtm,IP)+xyc
      DenXX(2,IAtm,2,IAtm,IP)=DenXX(2,IAtm,2,IAtm,IP)+DAIJ*VAOXX(3,I)
      xzc=DAIJ*VAOXX(4,I)
      DenXX(1,IAtm,3,IAtm,IP)=DenXX(1,IAtm,3,IAtm,IP)+xzc
      DenXX(3,IAtm,1,IAtm,IP)=DenXX(3,IAtm,1,IAtm,IP)+xzc
      yzc=DAIJ*VAOXX(5,I)
      DenXX(2,IAtm,3,IAtm,IP)=DenXX(2,IAtm,3,IAtm,IP)+yzc
      DenXX(3,IAtm,2,IAtm,IP)=DenXX(3,IAtm,2,IAtm,IP)+yzc
      DenXX(3,IAtm,3,IAtm,IP)=DenXX(3,IAtm,3,IAtm,IP)+DAIJ*VAOXX(6,I)
 185  CONTINUE
c -- (b) IAtm with JAtm
      If(JAtm.LT.IAtm) GO TO 18
      DAIJ = Two*DA(IJ)*VAOX(1,J)
      DenXX(1,IAtm,1,JAtm,IP)=DenXX(1,IAtm,1,JAtm,IP)+DAIJ*VAOX(1,I)
      DenXX(2,IAtm,1,JAtm,IP)=DenXX(2,IAtm,1,JAtm,IP)+DAIJ*VAOX(2,I)
      DenXX(3,IAtm,1,JAtm,IP)=DenXX(3,IAtm,1,JAtm,IP)+DAIJ*VAOX(3,I)
      DAIJ = Two*DA(IJ)*VAOX(2,J)
      DenXX(1,IAtm,2,JAtm,IP)=DenXX(1,IAtm,2,JAtm,IP)+DAIJ*VAOX(1,I)
      DenXX(2,IAtm,2,JAtm,IP)=DenXX(2,IAtm,2,JAtm,IP)+DAIJ*VAOX(2,I)
      DenXX(3,IAtm,2,JAtm,IP)=DenXX(3,IAtm,2,JAtm,IP)+DAIJ*VAOX(3,I)
      DAIJ = Two*DA(IJ)*VAOX(3,J)
      DenXX(1,IAtm,3,JAtm,IP)=DenXX(1,IAtm,3,JAtm,IP)+DAIJ*VAOX(1,I)
      DenXX(2,IAtm,3,JAtm,IP)=DenXX(2,IAtm,3,JAtm,IP)+DAIJ*VAOX(2,I)
      DenXX(3,IAtm,3,JAtm,IP)=DenXX(3,IAtm,3,JAtm,IP)+DAIJ*VAOX(3,I)
 18   CONTINUE
      DO 19 J=I+1,nbf(IP+1)
      JJ = INB(J)
      IJ = (JJ*(JJ-1))/2 + II
      JAtm = NBAtm(JJ)
      DAIJ = Two*DA(IJ)*VAO(J)
      If(Abs(DAIJ)*VMx.LT.thrsh) GO TO 195
c -- gradient density
      DenX(1,IAtm,IP) = DenX(1,IAtm,IP) - DAIJ*VAOX(1,I)
      DenX(2,IAtm,IP) = DenX(2,IAtm,IP) - DAIJ*VAOX(2,I)
      DenX(3,IAtm,IP) = DenX(3,IAtm,IP) - DAIJ*VAOX(3,I)
c -- Hessian density
c -- (a) IAtm with IAtm
      DenXX(1,IAtm,1,IAtm,IP)=DenXX(1,IAtm,1,IAtm,IP)+DAIJ*VAOXX(1,I)
      xyc=DAIJ*VAOXX(2,I)
      DenXX(1,IAtm,2,IAtm,IP)=DenXX(1,IAtm,2,IAtm,IP)+xyc
      DenXX(2,IAtm,1,IAtm,IP)=DenXX(2,IAtm,1,IAtm,IP)+xyc
      DenXX(2,IAtm,2,IAtm,IP)=DenXX(2,IAtm,2,IAtm,IP)+DAIJ*VAOXX(3,I)
      xzc=DAIJ*VAOXX(4,I)
      DenXX(1,IAtm,3,IAtm,IP)=DenXX(1,IAtm,3,IAtm,IP)+xzc
      DenXX(3,IAtm,1,IAtm,IP)=DenXX(3,IAtm,1,IAtm,IP)+xzc
      yzc=DAIJ*VAOXX(5,I)
      DenXX(2,IAtm,3,IAtm,IP)=DenXX(2,IAtm,3,IAtm,IP)+yzc
      DenXX(3,IAtm,2,IAtm,IP)=DenXX(3,IAtm,2,IAtm,IP)+yzc
      DenXX(3,IAtm,3,IAtm,IP)=DenXX(3,IAtm,3,IAtm,IP)+DAIJ*VAOXX(6,I)
  195 CONTINUE
c -- (b) IAtm with JAtm
      If(JAtm.LT.IAtm) GO TO 19
      DAIJ = Two*DA(IJ)*VAOX(1,J)
      DenXX(1,IAtm,1,JAtm,IP)=DenXX(1,IAtm,1,JAtm,IP)+DAIJ*VAOX(1,I)
      DenXX(2,IAtm,1,JAtm,IP)=DenXX(2,IAtm,1,JAtm,IP)+DAIJ*VAOX(2,I)
      DenXX(3,IAtm,1,JAtm,IP)=DenXX(3,IAtm,1,JAtm,IP)+DAIJ*VAOX(3,I)
      DAIJ = Two*DA(IJ)*VAOX(2,J)
      DenXX(1,IAtm,2,JAtm,IP)=DenXX(1,IAtm,2,JAtm,IP)+DAIJ*VAOX(1,I)
      DenXX(2,IAtm,2,JAtm,IP)=DenXX(2,IAtm,2,JAtm,IP)+DAIJ*VAOX(2,I)
      DenXX(3,IAtm,2,JAtm,IP)=DenXX(3,IAtm,2,JAtm,IP)+DAIJ*VAOX(3,I)
      DAIJ = Two*DA(IJ)*VAOX(3,J)
      DenXX(1,IAtm,3,JAtm,IP)=DenXX(1,IAtm,3,JAtm,IP)+DAIJ*VAOX(1,I)
      DenXX(2,IAtm,3,JAtm,IP)=DenXX(2,IAtm,3,JAtm,IP)+DAIJ*VAOX(2,I)
      DenXX(3,IAtm,3,JAtm,IP)=DenXX(3,IAtm,3,JAtm,IP)+DAIJ*VAOX(3,I)
 19   CONTINUE
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE MakeCphfC(dft,    NPP,    NBas,   nbf,    NAtoms,
     $                     natdo,  thrsh,  WGHT,   prara,  prarb,
     $                     pga,    pgc,    praga,  pragb,  pragc,
     $                     pgaga,  pgagb,  pgagc,  pgcgc,  DEN1,
     $                     DM1,    DENX,   VAO,    VAOX,   INB,
     $                     VM,     INDX,   Den1P,  Den1X1, Den1X2,
     $                     Den1X3, G1P,    SV,     SW1,    SW2,
     $                     SW3,    FDA,    td1,    tg1,    tquad,
     $                     tsw)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms the gradient density at the current grid point and carries
C  out numerical integration and accumulates contribution from
C  current grid point into partial derivative Fock matrices
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft   -  method flag (NOTE: All methods include Slater exchange)
C           1 - 3 - local correlation
C           >3 - nonlocal
C  NPP   - number of contributing (non-zero) grid points this batch
C  NBas  - total number of basis functions
C  nbf   - indexing array to location of "non-zero" basis functions
C  NAtoms- number of (symmetry unique) atoms this pass
C  NatDo - number of atoms not yet converged, for wich quadrature
C          has to be performed
C  thrsh - threshold for neglect of contribution
C  WGHT  - grid quadrature weights
C  pra   - Functional derivative w.r.t. alpha density at grid points
C  prara - Functional 2nd deriv. w.r.t. alpha density at grid points
C  prarb - Funct. 2nd deriv. w.r.t. alpha and beta density (DFT > 3)
C  pga   - Funct. deriv. w.r.t. alpha gradient (DFT > 3)
C  pgc   - Funct. deriv. w.r.t. alpha beta gradient (DFT > 3)
C  praga - Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (DFT > 3)
C  pragb - Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C          (DFT > 3)
C  pragc - Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C          (DFT > 3)
C  pgaga - Funct. 2nd. deriv. w.r.t. alpha grad. (DFT > 3)
C  pgagb - Funct. 2nd. deriv. w.r.t. alpha and beta grad. (DFT > 3)
C  pgagc - Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C          (DFT > 3)
C  pgcgc - Funct. 2nd. deriv. w.r.t. alpha beta grad. (DFT > 3)
C  DEN1  - current solution to CPHF, lower triangle (closed-shell)
C  DM1   - array containing maximum 1st-order density element per column
C  DENX  - density gradient at each grid point (dft > 3)
C  VAO   - "non-zero" basis function values at grid points
C  VAOX  - basis function derivatives at grid points (dft > 3)
C  INB   - indexing array for non-zero entries to VAO
C  VM    - array containing maximum magnitude AO per grid point
C  INDX  - index into contributing columns of VAO
C  Den1P - storage for 1st-order density at grid point
C  Den1XP- storage for 1st-order density gradient
C          at grid point (dft > 3)
C  G1P   - storage for in situ 1st-order density gradient invariant
C          at grid point (dft > 3)
C  SV    - storage for coefficient of vao(i)*vao(j) (dft > 3)
C  SW    - storage for coefficient of
C            (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C
C  on exit
C
C  FDA     -  contribution to partial derivative Fock matrices
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VM(*)
      DIMENSION DEN1(3*NAtoms,*),DM1(*)
      DIMENSION DENX(3,*)
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
      DIMENSION Den1P(3*NAtoms), G1P(3*NAtoms)
      DIMENSION Den1X1(3*Natoms),Den1X2(3*Natoms),Den1X3(3*Natoms)
      DIMENSION SV(3*NAtoms),SW1(3*NAtoms),SW2(3*NAtoms),SW3(3*NAtoms)
      DIMENSION FDA(3*NAtoms,NBas*(NBas+1)/2)
      dimension pga(npp),pgc(npp)
      dimension prara(npp),prarb(npp),praga(npp),pragb(npp),pragc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0)
      PARAMETER  (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d+17)
C
C
      NAt3 = 3*NAtoms
      natd3=3*NatDo
c
      IF(dft.LE.3) THEN
C
C  local only
C
        DO 50 IP=1,NPP
        IPP = INDX(IP)
        VMx = VM(IPP)
        rara = prara(IP)*WGHT(IPP)
        CALL ZeroIT(Den1P,NAt3)
c
c -- form 1st-order density (incorporating effect of derivative potential)
c
c    NOTE: for the closed shell case, the array DEN1 contains the
c          total density matrix, while here DEN1P will contain
c          only half the first order density at the current
c          grid point.
c
        vmx2=vmx*vmx
        thrsh1=thrsh/vmx2
        call secund(t1)
        DO 20 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        Val = VAO(I)*rara
        DMx = DM1(II)   ! maximum element of 1st-order density
        If(DMx*VMx*Abs(Val).LT.thrsh1) GO TO 20
        DO 19 J=nbf(IPP)+1,I-1
        JJ = INB(J)
        IJ = IT + JJ
c -- first-order density
        VAIJ = Val*VAO(J)
        If(Abs(VAIJ)*DMx.GT.thrsh1) Then
          DO IA=1,natd3
            Den1P(IA) = Den1P(IA) + DEN1(IA,IJ)*VAIJ
          EndDO
        EndIf
 19     CONTINUE
        Val2 = half*val*VAO(I)
        If(DMx*Abs(Val2).GT.thrsh1) Then
          IJ = IT+II
          DO IA=1,natd3
            Den1P(IA) = Den1P(IA) + DEN1(IA,IJ)*Val2
          EndDO
        EndIf
 20     CONTINUE
        call secund(t2)
        td1=td1+t2-t1
c
c -- get maximum value of 1st-order density at this grid point
        call secund(t1)
        Call absmax(NAt3,Den1P,iixx,DMaxyz)
        DMax1 = DMaxyz*VMx
        If(DMax1*VMx.LT.thrsh) GO TO 49
c
c -- now do numerical integration
        thrsh1=thrsh/dmax1
        thrsh2=thrsh/dmaxyz
        DO 40 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        Val = VAO(I)
        IF(Abs(Val).GT.thrsh1) THEN
          DO 30 J=nbf(IPP)+1,I
          JJ = INB(J)
          IJ = IT + JJ
          VAIJ = Val*VAO(J)
          If(Abs(VAIJ).GT.thrsh2) Then
            DO 29 IA=1,natd3
              FDA(IA,IJ) = FDA(IA,IJ) + VAIJ*Den1P(IA)
 29         CONTINUE
          EndIf
 30       CONTINUE
        ENDIF
 40     CONTINUE
 49     continue
        call secund(t3)
        tquad=tquad+t3-t2
 50     CONTINUE
cc
      ELSE
cc
C
C   non-local dft
C
        DO 100 IP=1,NPP
          IPP = INDX(IP)
          WG = WGHT(IPP)
          VMx = VM(IPP)
          ga = pga(ip)*wg
          gc = pgc(ip)*wg
          rara = (prara(ip)+prarb(ip))*wg
          raga = praga(ip)*wg
          ragb = pragb(ip)*wg
          ragc = pragc(ip)*wg
          gaga = pgaga(ip)*wg
          gagb = pgagb(ip)*wg
          gagc = pgagc(ip)*wg
          gcgc = pgcgc(ip)*wg
c
c  some sums of the potentials that will be used later
c
          prr=rara
          prg=raga+ragb+ragc
          pgg=gaga+gagb+Two*gagc+Half*gcgc
          pg=ga+Half*gc
c
c  density gradient at current point
c
          DX=DENX(1,IPP)
          DY=DENX(2,IPP)
          DZ=DENX(3,IPP)
c
c  zero out first order densities and gradients
c
          CALL ZeroIT(Den1P,NAt3)
          CALL ZeroIT(Den1X1,NAt3)
          CALL ZeroIT(Den1X2,NAt3)
          CALL ZeroIT(Den1X3,NAt3)
c
c  cutoff for the first order densities
c
          call secund(t1)
          vmx2=vmx*vmx
          thrsh1=thrsh/(four*vmx2)
          fdxyz=(dx+dy+dz)*four
          fpg=four*pg
          pfprod=prg+fdxyz*pgg
          denomij=max(abs(prr+fdxyz*prg),abs(dx*pfprod+fpg),
     $                abs(dy*pfprod+fpg),abs(dz*pfprod+fpg))
          if(denomij.gt.epsi)then
             thrshm=thrsh1/denomij
          else
             goto 100
          endif
c
c -- form 1st-order density and density gradient.
c    Unlike the local density case above, we do not
c    multiply by the potential derivatives at this stage
c
c    NOTE: for the closed shell case, the array DEN1 contains the
c          total density matrix, while here DEN1P and DEN1XP will
c          contain only half the first order density and density
c          gradient at the current grid point.
c
          DO 70 I=nbf(IPP)+1,nbf(IPP+1)
            II = INB(I)
            DMx = DM1(II)       ! max. element of first order density
            if(dmx.gt.epsi)then
              thtest=thrshm/dmx
            else
              goto 70
            endif
            IT = (II*(II-1))/2
            ValI = VAO(I)
            ValXI = VAOX(1,I)
            ValYI = VAOX(2,I)
            ValZI = VAOX(3,I)
            abvali=abs(vali)
            valt=max(abvali,abs(valxi),abs(valyi),abs(valzi))
            if((valt+abvali)*vmx.lt.thtest)goto 70
            DO 69 J=nbf(IPP)+1,I-1
              JJ = INB(J)
              ValJ = VAO(J)
              VAIJ = ValI*ValJ
              V1X=ValXI*ValJ+ValI*VAOX(1,J)
              V1Y=ValYI*ValJ+ValI*VAOX(2,J)
              V1Z=ValZI*ValJ+ValI*VAOX(3,J)
              valt=max(Abs(vaij),Abs(V1X),Abs(V1Y),Abs(V1Z))
              if(valt.lt.thtest) goto 69
              IJ = IT + JJ
              DO IA=1,natd3
                Den1P(IA) = Den1P(IA) + DEN1(IA,IJ)*VAIJ
              EndDO
              DO IA=1,natd3
                Den1X1(IA) = Den1X1(IA)+DEN1(IA,IJ)*V1X
              EndDO
              DO IA=1,natd3
                Den1X2(IA) = Den1X2(IA)+DEN1(IA,IJ)*V1Y
              EndDO
              DO IA=1,natd3
                Den1X3(IA) = Den1X3(IA)+DEN1(IA,IJ)*V1Z
              EndDO
 69         CONTINUE
            ValI2=ValI*ValI*Half
            V1X=ValXI*ValI
            V1Y=ValYI*ValI
            V1Z=ValZI*ValI
            valt=max(Abs(vali2),Abs(V1X),Abs(V1Y),Abs(V1Z))
            If(Valt.LT.thtest) GO TO 70
            IJ = IT+II
            DO IA=1,natd3
              Den1P(IA) = Den1P(IA) + DEN1(IA,IJ)*ValI2
            EndDO
            DO IA=1,natd3
              Den1X1(IA) = Den1X1(IA)+DEN1(IA,IJ)*V1X
            EndDO
            DO IA=1,natd3
              Den1X2(IA) = Den1X2(IA)+DEN1(IA,IJ)*V1Y
            EndDO
            DO IA=1,natd3
              Den1X3(IA) = Den1X3(IA)+DEN1(IA,IJ)*V1Z
            EndDO
 70       CONTINUE
        call secund(t2)
        td1=td1+t2-t1
c
c -- form 1st-order gradient invariant.
c
          DO IA=1,natd3
            G1P(IA)=Two*(DX*Den1X1(IA)+DY*Den1X2(IA)+DZ*Den1X3(IA))
          EndDO
          call secund(t3)
          tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V and W at page 7434 of Johnson and Frisch).
c
          prg2=two*prg
          pgg2=two*pgg
          pg2=two*pg
          DO IA=1,natd3
            sv(ia)=rara*den1p(ia)+prg*g1p(ia)
          EndDO
          DO IA=1,natd3
            vd=prg2*den1p(ia)+pgg2*g1p(ia)
            SW1(IA)=pg2*den1x1(IA)+DX*VD
            SW2(IA)=pg2*den1x2(IA)+DY*VD
            SW3(IA)=pg2*den1x3(IA)+DZ*VD
          EndDO
          call secund(t4)
          tsw=tsw+t4-t3
c
c   test the maximum absolute value of coefficients V and W
c
          Call absmax(NAtd3,SV,isv,svmax)
          Call absmax(NAtd3,SW1,isw,swmax1)
          Call absmax(NAtd3,SW2,isw,swmax2)
          Call absmax(NAtd3,SW3,isw,swmax3)
          swmax=max(swmax1,swmax2,swmax3)
          tswmax=Three*swmax
          if((svmax+tswmax+tswmax)*VMx2.lt.thrsh) GO TO 99
          thtest=thrsh/vmx
c
c   we are now ready for the numerical integration
c
          DO 90 I=nbf(IPP)+1,nbf(IPP+1)
            ValI = VAO(I)
            ValIX = VAOX(1,I)
            ValIY = VAOX(2,I)
            ValIZ = VAOX(3,I)
            valm=abs(vali)
            valmx=max(abs(valix),abs(valiy),abs(valiz))
            if(svmax*valm+tswmax*(valmx+valm).LT.thtest) GO TO 90
            II = INB(I)
            IT = (II*(II-1))/2
            DO 80 J=nbf(IPP)+1,I
              VALJ = VAO(J)
              VIJ=VALI*VALJ
              VIJX=VALIX*VALJ+VAOX(1,J)*VALI
              VIJY=VALIY*VALJ+VAOX(2,J)*VALI
              VIJZ=VALIZ*VALJ+VAOX(3,J)*VALI
              if(svmax*abs(vij)+swmax*(abs(vijx)+abs(vijy)+
     $            abs(vijz)).lt.thrsh)  GO TO 80
              JJ = INB(J)
              IJ = IT + JJ
              DO IA=1,natd3
                FDA(IA,IJ) = FDA(IA,IJ) + (SV(IA)*VIJ+
     $           SW1(IA)*VIJX+SW2(IA)*VIJY+SW3(IA)*VIJZ)
              EndDO
 80         CONTINUE
 90       CONTINUE
c
 99     continue
        call secund(t5)
        tquad=tquad+t5-t4
c
 100    CONTINUE
C
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE MakeCphfU(dft,    NPP,    NBas,   nbf,    NAtoms,
     $                     NatDo,  thrsh,  WGHT,   prara,  prbrb,
     $                     prarb,  pga,    pgb,    pgc,    praga,
     $                     pragb,  pragc,  prbga,  prbgb,  prbgc,
     $                     pgaga,  pgagb,  pgagc,  pgbgb,  pgbgc,
     $                     pgcgc,  DEN1A,  DEN1B,  DM1,    DENXA,
     $                     DENXB,  VAO,    VAOX,   INB,    VM,
     $                     INDX,   Den1PA, den1PB, Den1XA1,Den1XA2,
     $                     Den1XA3,Den1XB1,Den1XB2,Den1XB3,G1PAA,
     $                     G1PBB,  G1PAB,  SVA,    SVB,    SWA1,
     $                     SWA2,   SWA3,   SWB1,   SWB2,   SWB3,
     $                     FDA,    FDB,    td1,    tg1,    tquad,
     $                     tsw)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms the gradient density at the current grid point and carries
C  out numerical integration and accumulates contribution from
C  current grid point into partial derivative Fock matrices
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft   -  method flag (NOTE: All methods include Slater exchange)
C           1 - 3 - local correlation
C           >3 - nonlocal
C  NPP   - number of contributing (non-zero) grid points this batch
C  NBas  - total number of basis functions
C  nbf   - indexing array to location of "non-zero" basis functions
C  NAtoms- number of (symmetry unique) atoms this pass
C  NatDo -  number of atoms not yet converged, for wich quadrature
C           has to be performed
C  thrsh - threshold for neglect of contribution
C  WGHT  - grid quadrature weights
C  prara - Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb - Functional 2nd deriv. w.r.t. abeta density at grid points
C  prarb - Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   - Funct. deriv. w.r.t. alpha gradient (DFT > 3)
C  pgb   - Funct. deriv. w.r.t. beta gradient (DFT > 3)
C  pgc   - Funct. deriv. w.r.t. alpha beta gradient (DFT > 3)
C  praga - Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (DFT > 3)
C  pragb - Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C          (DFT > 3)
C  pragc - Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C          (DFT > 3)
C  prbga - Funct. 2nd. deriv. w.r.t. beta dens. and  alpa grad.
C  prbgb - Funct. 2nd. deriv. w.r.t. beta dens. and  grad. (DFT > 3)
C  prbgc - Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C          (DFT > 3)
C  pgaga - Funct. 2nd. deriv. w.r.t. alpha grad. (DFT > 3)
C  pgagb - Funct. 2nd. deriv. w.r.t. alpha and beta grad. (DFT > 3)
C  pgagc - Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C          (DFT > 3)
C  pgbgb - Funct. 2nd. deriv. w.r.t. beta grad. (DFT > 3)
C  pgbgc - Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C          (DFT > 3)
C  pgcgc - Funct. 2nd. deriv. w.r.t. alpha beta grad. (DFT > 3)
C  DEN1A - current alpha solution to CPHF, lower triangle
C  DEN1B - current beta solution to CPHF, lower triangle
C  DM1   - array containing maximum 1st-order density element per column
C  DENXA - alpha density gradient at each grid point (dft > 3)
C  DENXB - beta density gradient at each grid point (dft > 3)
C  VAO   - "non-zero" basis function values at grid points
C  VAOX  - basis function derivatives at grid points (dft > 3)
C  INB   - indexing array for non-zero entries to VAO
C  VM    - array containing maximum magnitude AO per grid point
C  INDX  - index into contributing columns of VAO
C  Den1PA- storage for 1st-order alpha density at grid point
C  Den1PB- storage for 1st-order beta density at grid point
C  Den1XPA- storage for 1st-order alpha density gradient
C           at grid point (dft > 3)
C  Den1XPB- storage for 1st-order beta density gradient
C           at grid point (dft > 3)
C  G1PAA - storage for in situ 1st-order alpha density gradient
C          invariant at grid point (dft > 3)
C  G1PBB - storage for in situ 1st-order beta density gradient
C          invariant at grid point (dft > 3)
C  G1PAB - storage for in situ 1st-order alpha beta density gradient
C          invariant at grid point (dft > 3)
C  SVA   - storage for alpha coefficient of vao(i)*vao(j)
C  SVB   - storage for beta coefficient of vao(i)*vao(j)
C  SWA   - storage for alpha coefficient of
C            (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  SWB   - storage for beta coefficient of
C            (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C
C  on exit
C
C  FDA     -  contribution to partial alpha derivative Fock matrices
C  FDB     -  contribution to partial beta derivative Fock matrices
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VM(*)
      DIMENSION DEN1A(3*NAtoms,*),DM1(*)
      DIMENSION DEN1B(3*NAtoms,*)
      DIMENSION DENXA(3,*)
      DIMENSION DENXB(3,*)
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
      DIMENSION Den1PA(3*NAtoms), G1PAA(3*NAtoms)
      DIMENSION Den1PB(3*NAtoms), G1PBB(3*NAtoms)
      DIMENSION G1PAB(3*NAtoms)
      DIMENSION Den1XA1(3*NAtoms),Den1XA2(3*NAtoms),Den1XA3(3*NAtoms)
      DIMENSION Den1XB1(3*NAtoms),Den1XB2(3*NAtoms),Den1XB3(3*NAtoms)
      DIMENSION SVA(3*NAtoms),SWA1(3*NAtoms),SWA2(3*NAtoms),
     $          SWA3(3*NAtoms)
      DIMENSION SVB(3*NAtoms),SWB1(3*NAtoms),SWB2(3*NAtoms),
     $          SWB3(3*NAtoms)
      DIMENSION FDA(3*NAtoms,NBas*(NBas+1)/2)
      DIMENSION FDB(3*NAtoms,NBas*(NBas+1)/2)
      dimension pga(npp),pgb(npp),pgc(npp)
      dimension prara(npp),prarb(npp),praga(npp),pragb(npp),pragc(npp)
      dimension prbrb(npp),prbga(npp),prbgb(npp),prbgc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp)
      dimension pgbgb(npp),pgbgc(npp),pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0)
      PARAMETER  (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d+17)
C
C
      NAt3 = 3*NAtoms
      natd3=3*NatDo
c
      IF(dft.LE.3) THEN
C
C  local only
C
        DO 50 IP=1,NPP
        IPP = INDX(IP)
        VMx = VM(IPP)
        rara = prara(IP)*WGHT(IPP)
        rbrb = prbrb(IP)*WGHT(IPP)
        rarb = prarb(IP)*WGHT(IPP)
        CALL ZeroIT(Den1PA,NAt3)
        CALL ZeroIT(Den1PB,NAt3)
c
c -- form 1st-order density
c
        vmx2=vmx*vmx
        thrsh1=thrsh/vmx2
        call secund(t1)
        DO 20 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        Val = VAO(I)
        Val2 = Val+Val
        DMx2 = two*DM1(II)   ! maximum element of 1st-order density
        If(DMx2*VMx*Abs(Val).LT.thrsh1) GO TO 20
        DO 19 J=nbf(IPP)+1,I-1
        JJ = INB(J)
        IJ = IT + JJ
        VAIJ2 = Val2*VAO(J)
        If(Abs(VAIJ2)*DMx2.GT.thrsh1) Then
          DO IA=1,NatD3
            Den1PA(IA) = Den1PA(IA) + DEN1A(IA,IJ)*VAIJ2
            Den1PB(IA) = Den1PB(IA) + DEN1B(IA,IJ)*VAIJ2
          EndDO
        EndIf
 19     CONTINUE
        ValS = val*val
        If(DMx2*Abs(ValS).GT.thrsh1) Then
          IJ = IT+II
          DO IA=1,NatD3
            Den1PA(IA) = Den1PA(IA) + DEN1A(IA,IJ)*ValS
            Den1PB(IA) = Den1PB(IA) + DEN1B(IA,IJ)*ValS
          EndDO
        EndIf
 20     CONTINUE
        call secund(t2)
        td1=td1+t2-t1
c
c -- now form the coefficients that multiply the basis functions
c    in the expression for the Fock Matrix
c    (i.e., quantities V at page 7434 of Johnson and Frisch).
c
        DO IA=1,NatD3
          sva(ia)=rara*den1pa(ia)+rarb*den1pb(ia)
          svb(ia)=rbrb*den1pb(ia)+rarb*den1pa(ia)
        EndDO
        call secund(t3)
        td1=tsw+t3-t2
c
c -- get maximum value of sva and svb coefficienta
        Call absmax(NAtd3,sva,iixx,svamax)
        Call absmax(NAtd3,svb,iixx,svbmax)
        svmax=max(svamax,svbmax)
        DMax1 = svmax*VMx
        If(DMax1*VMx.LT.thrsh) GO TO 49
c
c -- now do numerical integration
        thrsh1=thrsh/dmax1
        thrsh2=thrsh/svmax
        DO 40 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        Val = VAO(I)
        IF(Abs(Val).GT.thrsh1) THEN
          DO 30 J=nbf(IPP)+1,I
          JJ = INB(J)
          IJ = IT + JJ
          VAIJ = Val*VAO(J)
          If(Abs(VAIJ).GT.thrsh2) Then
            DO 29 IA=1,NatD3
              FDA(IA,IJ) = FDA(IA,IJ) + VAIJ*sva(IA)
              FDB(IA,IJ) = FDB(IA,IJ) + VAIJ*svb(IA)
 29         CONTINUE
          EndIf
 30       CONTINUE
        ENDIF
 40     CONTINUE
 49     continue
        call secund(t4)
        tquad=tquad+t4-t3
 50     CONTINUE
cc
      ELSE
cc
C
C   non-local dft
C
        DO 100 IP=1,NPP
          IPP = INDX(IP)
          WG = WGHT(IPP)
          VMx = VM(IPP)
          ga = pga(ip)*wg
          gb = pgb(ip)*wg
          gc = pgc(ip)*wg
          rara = prara(ip)*wg
          rbrb = prbrb(ip)*wg
          rarb = prarb(ip)*wg
          raga = praga(ip)*wg
          ragb = pragb(ip)*wg
          ragc = pragc(ip)*wg
          rbga = prbga(ip)*wg
          rbgb = prbgb(ip)*wg
          rbgc = prbgc(ip)*wg
          gaga = pgaga(ip)*wg
          gagb = pgagb(ip)*wg
          gagc = pgagc(ip)*wg
          gbgb = pgbgb(ip)*wg
          gbgc = pgbgc(ip)*wg
          gcgc = pgcgc(ip)*wg
c
c  density gradient at current point
c
          DAX=DENXA(1,IPP)
          DAY=DENXA(2,IPP)
          DAZ=DENXA(3,IPP)
          DBX=DENXB(1,IPP)
          DBY=DENXB(2,IPP)
          DBZ=DENXB(3,IPP)
          DAX2=DAX+DAX
          DAY2=DAY+DAY
          DAZ2=DAZ+DAZ
          DBX2=DBX+DBX
          DBY2=DBY+DBY
          DBZ2=DBZ+DBZ
c
c  zero out first order densities and gradients
c
          CALL ZeroIT(Den1PA,NAt3)
          CALL ZeroIT(Den1PB,NAt3)
          CALL ZeroIT(Den1XA1,NAt3)
          CALL ZeroIT(Den1XA2,NAt3)
          CALL ZeroIT(Den1XA3,NAt3)
          CALL ZeroIT(Den1XB1,NAt3)
          CALL ZeroIT(Den1XB2,NAt3)
          CALL ZeroIT(Den1XB3,NAt3)
          call secund(t1)
c
c  cutoff for the first order densities
c
          call secund(t1)
          vmx2=vmx*vmx
          thrsh1=thrsh/vmx2
c
          daxyz2=dax2+day2+daz2
          dbxyz2=dbx2+dby2+dbz2
          daxyz4=daxyz2+daxyz2
          dbxyz4=dbxyz2+dbxyz2
          ga2=ga+ga
          gb2=gb+gb
          ga4=ga2+ga2
          gb4=gb2+gb2
          gc2=gc+gc
c
          d1aa=rara+dbxyz2*ragc+daxyz4*raga
          d1ab=rarb+daxyz2*ragc+dbxyz4*ragb
          d1ba=rarb+dbxyz2*rbgc+daxyz4*rbga
          d1bb=rbrb+daxyz2*rbgc+dbxyz4*rbgb
          den1=max(abs(d1aa),abs(d1ab),abs(d1ba),abs(d1bb))
c
          sumaa1=ragc+dbxyz2*gcgc+daxyz4*gagc
          sumaa2=raga+dbxyz2*gagc+daxyz4*gaga
          d1xaa=dbx*sumaa1+dax*sumaa2+ga4
          d1yaa=dby*sumaa1+day*sumaa2+ga4
          d1zaa=dbz*sumaa1+daz*sumaa2+ga4
c
          sumab1=rbgc+daxyz2*gcgc+dbxyz4*gbgc
          sumab2=rbga+daxyz2*gagc+dbxyz4*gagb
          d1xab=dbx*sumab1+dax*sumab2+gc2
          d1yab=dby*sumab1+day*sumab2+gc2
          d1zab=dbz*sumab1+daz*sumab2+gc2
          den2=max(abs(d1xaa),abs(d1yaa),abs(d1zaa),abs(d1xab),
     $             abs(d1yab),abs(d1zab))
c
          sumbb2=rbgb+daxyz2*gbgc+dbxyz4*gbgb
          d1xbb=dax*sumab1+dbx*sumbb2+gb4
          d1ybb=day*sumab1+dby*sumbb2+gb4
          d1zbb=daz*sumab1+dbz*sumbb2+gb4
c
          sumba2=ragb+dbxyz2*gbgc+daxyz4*gagb
          d1xba=dax*sumaa1+dbx*sumba2+gc2
          d1yba=day*sumaa1+dby*sumba2+gc2
          d1zba=daz*sumaa1+dbz*sumba2+gc2
          den3=max(abs(d1xbb),abs(d1ybb),abs(d1zbb),abs(d1xba),
     $             abs(d1yba),abs(d1zba))
c
          denij=max(den1,den2,den3)
          if(denij.gt.epsi)then
             thrshm=thrsh1/denij
          else
             goto 100
          endif
c
c -- form 1st-order density and density gradient.
c
          DO 70 I=nbf(IPP)+1,nbf(IPP+1)
            II = INB(I)
            DMx = DM1(II)       ! max. element of first order density
            if(dmx.gt.epsi)then
              thtest=thrshm/dmx
            else
              goto 70
            endif
            IT = (II*(II-1))/2
            ValI = VAO(I)
            ValXI = VAOX(1,I)
            ValYI = VAOX(2,I)
            ValZI = VAOX(3,I)
            abvali=abs(vali)
            valt=max(abvali,abs(valxi),abs(valyi),abs(valzi))
            if((valt+abvali)*vmx.lt.thtest)goto 70
            DO 69 J=nbf(IPP)+1,I-1
              JJ = INB(J)
              ValJ = VAO(J)
              VAIJ2 = two*ValI*ValJ
              V1X2=two*(ValXI*ValJ+ValI*VAOX(1,J))
              V1Y2=two*(ValYI*ValJ+ValI*VAOX(2,J))
              V1Z2=two*(ValZI*ValJ+ValI*VAOX(3,J))
              valt=max(Abs(vaij2),Abs(V1X2),Abs(V1Y2),Abs(V1Z2))
              if(valt.lt.thtest) goto 69
              IJ = IT + JJ
              DO IA=1,NatD3
                Den1PA(IA) = Den1PA(IA) + DEN1A(IA,IJ)*VAIJ2
                Den1PB(IA) = Den1PB(IA) + DEN1B(IA,IJ)*VAIJ2
              EndDO
              DO IA=1,NatD3
                Den1XA1(IA) = Den1XA1(IA) + DEN1A(IA,IJ)*V1X2
                Den1XB1(IA) = Den1XB1(IA) + DEN1B(IA,IJ)*V1X2
              EndDO
              DO IA=1,NatD3
                Den1XA2(IA) = Den1XA2(IA) + DEN1A(IA,IJ)*V1Y2
                Den1XB2(IA) = Den1XB2(IA) + DEN1B(IA,IJ)*V1Y2
              EndDO
              DO IA=1,NatD3
                Den1XA3(IA) = Den1XA3(IA) + DEN1A(IA,IJ)*V1Z2
                Den1XB3(IA) = Den1XB3(IA) + DEN1B(IA,IJ)*V1Z2
              EndDO
 69         CONTINUE
            ValIS=ValI*ValI
            V1X2=two*ValXI*ValI
            V1Y2=two*ValYI*ValI
            V1Z2=two*ValZI*ValI
            valt=max(Abs(valiS),Abs(V1X2),Abs(V1Y2),Abs(V1Z2))
            If(Valt.LT.thtest) GO TO 70
            IJ = IT+II
            DO IA=1,NatD3
              Den1PA(IA) = Den1PA(IA) + DEN1A(IA,IJ)*ValIS
              Den1PB(IA) = Den1PB(IA) + DEN1B(IA,IJ)*ValIS
            EndDO
            DO IA=1,NatD3
              Den1XA1(IA) = Den1XA1(IA) + DEN1A(IA,IJ)*V1X2
              Den1XB1(IA) = Den1XB1(IA) + DEN1B(IA,IJ)*V1X2
            EndDO
            DO IA=1,NatD3
              Den1XA2(IA) = Den1XA2(IA) + DEN1A(IA,IJ)*V1Y2
              Den1XB2(IA) = Den1XB2(IA) + DEN1B(IA,IJ)*V1Y2
            EndDO
            DO IA=1,NatD3
              Den1XA3(IA) = Den1XA3(IA) + DEN1A(IA,IJ)*V1Z2
              Den1XB3(IA) = Den1XB3(IA) + DEN1B(IA,IJ)*V1Z2
            EndDO
 70       CONTINUE
        call secund(t2)
        td1=td1+t2-t1
c
c -- form 1st-order gradient invariant.
c
          DO IA=1,NatD3
            G1PAA(IA)= Two*(DAX*Den1XA1(IA) +
     $                     DAY*Den1XA2(IA) + DAZ*Den1XA3(IA))
          EndDO
          DO IA=1,NatD3
            G1PBB(IA)= Two*(DBX*Den1XB1(IA) +
     $                     DBY*Den1XB2(IA) + DBZ*Den1XB3(IA))
          EndDO
          DO IA=1,NatD3
            G1PAB(IA)= DAX*Den1XB1(IA) + DBX*Den1XA1(IA) +
     $                 DAY*Den1XB2(IA) + DBY*Den1XA2(IA) +
     $                 DAZ*Den1XB3(IA) + DBZ*Den1XA3(IA)
          EndDO
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V and W at page 7434 of Johnson and Frisch).
c
c
c     coeff. V
c
        DO IA=1,NatD3
          sva(ia)=rara*Den1PA(ia)+rarb*Den1PB(ia)+
     $        raga*G1PAA(ia)+ragb*G1PBB(ia)+ragc*G1PAB(ia)
          svb(ia)=rarb*Den1PA(ia)+rbrb*Den1PB(ia)+
     $        rbga*G1PAA(ia)+rbgb*G1PBB(ia)+rbgc*G1PAB(ia)
        EndDO
c
c     coeff. W
c
        DO IA=1,NatD3
          sga=raga*Den1PA(ia)+rbga*Den1PB(ia)+
     $        gaga*G1PAA(ia)+gagb*G1PBB(ia)+gagc*G1PAB(ia)
          sgb=ragb*Den1PA(ia)+rbgb*Den1PB(ia)+
     $        gagb*G1PAA(ia)+gbgb*G1PBB(ia)+gbgc*G1PAB(ia)
          sgc=ragc*Den1PA(ia)+rbgc*Den1PB(ia)+
     $        gagc*G1PAA(ia)+gbgc*G1PBB(ia)+gcgc*G1PAB(ia)
          swa1(ia)=ga2*Den1XA1(ia)+gc*Den1XB1(ia)+DAX2*sga+DBX*sgc
          swa2(ia)=ga2*Den1XA2(ia)+gc*Den1XB2(ia)+DAY2*sga+DBY*sgc
          swa3(ia)=ga2*Den1XA3(ia)+gc*Den1XB3(ia)+DAZ2*sga+DBZ*sgc
          swb1(ia)=gb2*Den1XB1(ia)+gc*Den1XA1(ia)+DBX2*sgb+DAX*sgc
          swb2(ia)=gb2*Den1XB2(ia)+gc*Den1XA2(ia)+DBY2*sgb+DAY*sgc
          swb3(ia)=gb2*Den1XB3(ia)+gc*Den1XA3(ia)+DBZ2*sgb+DAZ*sgc
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c   test the maximum absolute value of coefficients V and W
c
          Call absmax(NAtd3,sva,iixx,svamax)
          Call absmax(NAtd3,svb,iixx,svbmax)
          svmax=max(svamax,svbmax)
          Call absmax(NAtd3,swa1,iixx,swamax1)
          Call absmax(NAtd3,swa2,iixx,swamax2)
          Call absmax(NAtd3,swa3,iixx,swamax3)
          Call absmax(NAtd3,swb1,iixx,swbmax1)
          Call absmax(NAtd3,swb2,iixx,swbmax2)
          Call absmax(NAtd3,swb3,iixx,swbmax3)
          swmax=max(swamax1,swamax2,swamax3,swbmax1,swbmax2,swbmax3)
          tswmax=Three*swmax
          if((svmax+tswmax+tswmax)*vmx*vmx.lt.thrsh) GO TO 99
          thtest=thrsh/vmx
c
c   we are now ready for the numerical integration
c
          DO 90 I=nbf(IPP)+1,nbf(IPP+1)
            ValI = VAO(I)
            ValIX = VAOX(1,I)
            ValIY = VAOX(2,I)
            ValIZ = VAOX(3,I)
            valm=abs(vali)
            valmx=max(abs(valix),abs(valiy),abs(valiz))
            if(svmax*valm+tswmax*(valmx+valm).LT.thtest) GO TO 90
            II = INB(I)
            IT = (II*(II-1))/2
            DO 80 J=nbf(IPP)+1,I
              VALJ = VAO(J)
              VIJ=VALI*VALJ
              VIJX=VALIX*VALJ+VAOX(1,J)*VALI
              VIJY=VALIY*VALJ+VAOX(2,J)*VALI
              VIJZ=VALIZ*VALJ+VAOX(3,J)*VALI
              if(svmax*abs(vij)+swmax*(abs(vijx)+abs(vijy)+
     $            abs(vijz)).lt.thrsh)  GO TO 80
              JJ = INB(J)
              IJ = IT + JJ
              DO IA=1,NatD3
               FDA(IA,IJ) = FDA(IA,IJ) + SVA(IA)*VIJ +
     $           SWA1(IA)*VIJX + SWA2(IA)*VIJY + SWA3(IA)*VIJZ
               FDB(IA,IJ) = FDB(IA,IJ) + SVB(IA)*VIJ +
     $           SWB1(IA)*VIJX + SWB2(IA)*VIJY + SWB3(IA)*VIJZ
              EndDO
 80         CONTINUE
 90       CONTINUE
c
 99     continue
        call secund(t5)
        tquad=tquad+t5-t4
c
 100    CONTINUE
C
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE POTSortF(dft,    lsemi,   NP,     thrsh,  nrec,
     $                    POld,   POTD,    INDX,   NPP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts this batch of potentials, determining which grid points
C  will be retained for the numerical integration
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  DFT flag
C             of importance here is if we have a non-local functional
C  lsemi   -  integer flag, if >0 use delta potential
C  NP      -  number of grid points this batch
C  thrsh   -  threshold for neglecting contribution
C  nrec    -  record counter for direct access file
C  POld    -  previous density potential derivatives at each grid point
C  POTD    -    ditto current density potential derivatives
C
C  on exit
C
C  INDX    -  index into contributing columns of POT
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION POld(NP),POTD(NP),INDX(NP)
      INTEGER dft
C
C
      NPP = 0
c
      IF(lsemi.GT.1) THEN
C
C  read old potential for current batch of grid points
C  replace on file by current potential derivatives
C
        READ(41,rec=nrec) POld
        WRITE(41,rec=nrec) POTD
C
C  form difference potential
C
        DO 20 IP=1,NP
        POTD(IP) = POTD(IP) - POld(IP)
        If(Abs(POTD(IP)).GT.thrsh) Then
          NPP = NPP+1
          INDX(NPP) = IP
          POTD(NPP) = POTD(IP)
        EndIf
 20     CONTINUE
cc
      ELSE
C
C  do not use difference potential
C  save current potential
C
        If(lsemi.GT.0) WRITE(41,rec=nrec) POTD
C
C  sort absolute potential
C
        DO 25 IP=1,NP
        If(Abs(POTD(IP)).GT.thrsh) Then
          NPP = NPP+1
          INDX(NPP) = IP
          POTD(NPP) = POTD(IP)
        EndIf
 25     CONTINUE
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE POTSortH(dft,    NP,     thrsh,  pra,    prara,
     $                    prarb,  pga,    pgc,    praga,  pragb,
     $                    pragc,  pgaga,  pgagb,  pgagc,  pgcgc,
     $                    INDX,   NPP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts this batch of potentials, determining which grid points
C  will be retained for the numerical integration
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  DFT flag
C             of importance here is if we have a non-local functional
C  NP      -  number of grid points this batch
C  thrsh   -  threshold for neglecting contribution
C  pra     -  Functional derivative w.r.t. alpha density at grid points
C  prara   -  Functional 2nd deriv. w.r.t. alpha density at grid points
C
C  for non local DFT (dft > 3) only
C
C  prarb -   Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -   Funct. deriv. w.r.t. alpha gradient
C  pgc   -   Funct. deriv. w.r.t. alpha beta gradient
C  praga -   Funct. 2nd. deriv. w.r.t. alpha dens. and  grad.
C  pragb -   Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C  pragc -   Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C  pgaga -   Funct. 2nd. deriv. w.r.t. alpha grad.
C  pgagb -   Funct. 2nd. deriv. w.r.t. alpha and beta grad.
C  pgagc -   Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C  pgcgc -   Funct. 2nd. deriv. w.r.t. alpha beta grad.
C
C  on exit
C
C  INDX    -  index into contributing columns of POT
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION INDX(NP)
      dimension pra(np),pga(np),pgc(np)
      dimension prara(np),prarb(np),praga(np),pragb(np),pragc(np)
      dimension pgaga(np),pgagb(np),pgagc(np),pgcgc(np)
      INTEGER dft
C
C
      NPP = 0
c
      IF(dft.GT.3) THEN
C
C  sort absolute potential & derivatives
C
        DO 20 IP=1,NP
        If(Abs(pra(IP)).GT.thrsh) Then
          NPP = NPP+1
          INDX(NPP) = IP
          pra(npp)=pra(ip)
          pga(npp)=pga(ip)
          pgc(npp)=pgc(ip)
          prara(npp)=prara(ip)
          prarb(npp)=prarb(ip)
          praga(npp)=praga(ip)
          pragb(npp)=pragb(ip)
          pragc(npp)=pragc(ip)
          pgaga(npp)=pgaga(ip)
          pgagb(npp)=pgagb(ip)
          pgagc(npp)=pgagc(ip)
          pgcgc(npp)=pgcgc(ip)
        EndIf
 20     CONTINUE
cc
      ELSE
C
C  sort absolute potential
C
        DO 25 IP=1,NP
        If(Abs(pra(IP)).GT.thrsh) Then
          NPP = NPP+1
          INDX(NPP) = IP
          pra(NPP) = pra(IP)
          prara(NPP) = prara(IP)
        EndIf
 25     CONTINUE
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE POTSortHU(dft,    NP,     thrsh,  pra,    prb,
     $                     prara,  prbrb,  prarb,  pga,    pgb,
     $                     pgc,    praga,  pragb,  pragc,  prbga,
     $                     prbgb,  prbgc,  pgaga,  pgagb,  pgagc,
     $                     pgbgb,  pgbgc,  pgcgc,  INDX,   NPP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts this batch of potentials, determining which grid points
C  will be retained for the numerical integration
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  DFT flag
C             of importance here is if we have a non-local functional
C  NP      -  number of grid points this batch
C  thrsh   -  threshold for neglecting contribution
C  pra     -  Functional derivative w.r.t. alpha density at grid points
C  prb     -  Functional derivative w.r.t. beta density at grid points
C  prara   -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb   -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb   -   Funct. 2nd deriv. w.r.t. alpha and beta density
C
C  for non local DFT (dft > 3) only
C
C  pga   -   Funct. deriv. w.r.t. alpha gradient
C  pgb   -   Funct. deriv. w.r.t. beta gradient
C  pgc   -   Funct. deriv. w.r.t. alpha beta gradient
C  praga -   Funct. 2nd. deriv. w.r.t. alpha dens. and  grad.
C  pragb -   Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C  pragc -   Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C  prbga -   Funct. 2nd. deriv. w.r.t. beta dens. and  alpha grad.
C  prbgb -   Funct. 2nd. deriv. w.r.t. beta dens. and  beta grad.
C  prbgc -   Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C  pgaga -   Funct. 2nd. deriv. w.r.t. alpha grad.
C  pgagb -   Funct. 2nd. deriv. w.r.t. alpha and beta grad.
C  pgagc -   Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C  pgbgb -   Funct. 2nd. deriv. w.r.t. beta grad.
C  pgbgc -   Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C  pgcgc -   Funct. 2nd. deriv. w.r.t. alpha beta grad.
C
C  on exit
C
C  INDX    -  index into contributing columns of POT
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION INDX(NP)
      dimension pra(np),prb(np),pga(np),pgb(np),pgc(np)
      dimension prara(np),prarb(np),praga(np),pragb(np),pragc(np)
      dimension prbrb(np),prbga(np),prbgb(np),prbgc(np)
      dimension pgaga(np),pgagb(np),pgagc(np),pgbgb(np)
      dimension pgbgc(np),pgcgc(np)
      INTEGER dft
C
C
      NPP = 0
c
      IF(dft.GT.3) THEN
C
C  sort absolute potential & derivatives
C
        DO 20 IP=1,NP
        If(Abs(pra(IP))+Abs(prb(IP)).GT.thrsh) Then
          NPP = NPP+1
          INDX(NPP) = IP
          pra(npp)=pra(ip)
          prb(npp)=prb(ip)
c
          pga(npp)=pga(ip)
          pgb(npp)=pgb(ip)
          pgc(npp)=pgc(ip)
c
          prara(npp)=prara(ip)
          prbrb(npp)=prbrb(ip)
          prarb(npp)=prarb(ip)
c
          praga(npp)=praga(ip)
          pragb(npp)=pragb(ip)
          pragc(npp)=pragc(ip)
          prbga(npp)=prbga(ip)
          prbgb(npp)=prbgb(ip)
          prbgc(npp)=prbgc(ip)
c
          pgaga(npp)=pgaga(ip)
          pgagb(npp)=pgagb(ip)
          pgagc(npp)=pgagc(ip)
          pgbgb(npp)=pgbgb(ip)
          pgbgc(npp)=pgbgc(ip)
          pgcgc(npp)=pgcgc(ip)
        EndIf
 20     CONTINUE
cc
      ELSE
C
C  sort absolute potential
C
        DO 25 IP=1,NP
        If(Abs(pra(IP))+Abs(prb(IP)).GT.thrsh) Then
          NPP = NPP+1
          INDX(NPP) = IP
          pra(NPP) = pra(IP)
          prb(NPP) = prb(IP)
          prara(NPP) = prara(IP)
          prbrb(NPP) = prbrb(IP)
          prarb(NPP) = prarb(IP)
        EndIf
 25     CONTINUE
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE hdersymm_hes(NTrans, NAtoms, ISYM,   NEqATM, TRANS,
     $                        HESS,   R,      U,      HSYM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine symmetrizes the DFT contribution to the analytical Hessian
C
C  ARGUMENTS
C
C  NTrans  -  number of (Abelian) symmetry operations
C  NAtoms  -  number of atoms
C  ISYM    -  list of Abelian symmetry operations
C               1 - reflection in YZ plane
C               2 - reflection in XZ plane
C               3 - reflection in XY plane
C               4 - C2 rotation about Z-axis
C               5 - C2 rotation about Y-axis
C               6 - C2 rotation about X-axis
C               7 - inversion through origin
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C             (will be generated internally)
C  HESS    -  current partial (unsymmetrized) Hessian
C  R       -  scratch space (for transformation matrix)
C  U       -  scratch space
C  HSYM    -  on exit contains symmetrized matrix
C
C
      DIMENSION ISYM(7),NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),
     $          HESS(3*NAtoms,3*NAtoms),R(3*NAtoms,3*NAtoms),
     $          U(3*NAtoms,3*NAtoms),HSYM(3*NAtoms,3*NAtoms)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=1.0d-7)
C
C
      NAt3 = 3*NAtoms
C
C  First form the symmetry operations
C
      CALL ZeroIT(TRANS,3*3*NTrans)
C
C   generate the symmetry operators
C
      DO 10 IOP=1,NTrans
      If(ISYM(IOP).EQ.1) Then
        TRANS(1,1,IOP) = -One
        TRANS(2,2,IOP) =  One
        TRANS(3,3,IOP) =  One
      Else If(ISYM(IOP).EQ.2) Then
        TRANS(1,1,IOP) =  One
        TRANS(2,2,IOP) = -One
        TRANS(3,3,IOP) =  One
      Else If(ISYM(IOP).EQ.3) Then
        TRANS(1,1,IOP) =  One
        TRANS(2,2,IOP) =  One
        TRANS(3,3,IOP) = -One
      Else If(ISYM(IOP).EQ.4) Then
        TRANS(1,1,IOP) = -One
        TRANS(2,2,IOP) = -One
        TRANS(3,3,IOP) =  One
      Else If(ISYM(IOP).EQ.5) Then
        TRANS(1,1,IOP) = -One
        TRANS(2,2,IOP) =  One
        TRANS(3,3,IOP) = -One
      Else If(ISYM(IOP).EQ.6) Then
        TRANS(1,1,IOP) =  One
        TRANS(2,2,IOP) = -One
        TRANS(3,3,IOP) = -One
      Else If(ISYM(IOP).EQ.7) Then
        TRANS(1,1,IOP) = -One
        TRANS(2,2,IOP) = -One
        TRANS(3,3,IOP) = -One
      EndIf
 10   CONTINUE
C
C  copy Hessian upper triangle into lower
C
      DO 20 I=1,NAt3
      DO 20 J=1,I
      HESS(I,J) = HESS(J,I)
 20   CONTINUE
C
C  now symmetrize
C
      CALL SymHES1(NAtoms, NTrans, NEqATM, TRANS,  HESS,
     $             R,      U,      thrsh,  HSYM)
C
C  symmetrized Hessian is in HSYM
C
      RETURN
      END
c =====================================================================
      SUBROUTINE SymHES1(NAtoms, NTrans, NEqATM, TRANS,  HOld,
     $                   R,      U,      thrsh,  HNew)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Symmetrizes a Cartesian Hessian matrix according to all
C  operations of the molecular point group
C  Additionally zeros all elements below thrsh
C  ** CURRENTLY USED FOR DFT CONTRIBUTION TO ANALYTICAL HESSIAN **
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  HOld    -  matrix (Hessian) to be symmetrized
C  R       -  scratch space (for transformation matrix)
C  U       -  scratch space
C  thrsh   -  threshold below which elements will be set to zero
C  HNew    -  on exit contains symmetrized matrix
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),
     $          HOld(3*NAtoms,3*NAtoms),R(3*NAtoms,3*NAtoms),
     $          U(3*NAtoms,3*NAtoms),HNew(3*NAtoms,3*NAtoms)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      NAT3 = 3*NAtoms
      CALL CpyVEC(NAT3*NAT3,HOld,HNew)
c
      DO 70 IOP=1,NTrans
C
C  Construct the R matrix for each symmetry operation
C
      CALL ZeroIT(R,NAT3*NAT3)
c
      DO 20 IAtm=1,NAtoms
      II = (IAtm-1)*3
      JAtm = NEqATM(IAtm,IOP)
      JJ = (JAtm-1)*3
c
      DO 10 K=1,3
      KK = II+K
      DO 10 L=1,3
      LL = JJ+L
      R(KK,LL) = TRANS(L,K,IOP)
 10   CONTINUE
c
 20   CONTINUE
C
C  Form R * HESS * R(t)
C
      CALL ZeroIT(U,NAT3*NAT3)
c
      DO 40 J=1,NAT3
      DO 40 K=1,NAT3
      VAL = HOld(K,J)
      DO 30 I=1,NAT3
      U(I,J) = U(I,J) + R(I,K)*VAL
 30   CONTINUE
 40   CONTINUE
c
      DO 60 J=1,NAT3
      DO 60 K=1,NAT3
      VAL = R(J,K)
      DO 50 I=1,NAT3
      HNew(I,J) = HNew(I,J) + U(I,K)*VAL
 50   CONTINUE
 60   CONTINUE
c
 70   CONTINUE
C
C  zero out elements below thrsh
C
      DO 80 I=1,NAT3
      DO 80 J=1,I
      If(Abs(HNew(I,J)).LT.thrsh) HNew(I,J) = Zero
      HNew(J,I) = HNew(I,J)
 80   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE GHMu(X,Y,Z,XB,YB,ZB,XA,YA,ZA,rAB,xaAB,xmuAB,
     $                  Gradx,Hess)
      IMPLICIT REAL*8(A-H,O-Z)
C
C    MM (01/09/2003), most of this subroutine has been generated
C    with maxima
C
C     computes the gradient and hessian of the hyperbolic coordinates
C     mu(A,B) used in the calculations of Becke's quadrature weights,
C     taking into account the modifications due to size adjustment
C
C     for the hyperbolic coordinates definition and their gradient
C     see: Johnson, Gill and Pople, J.Chem.Phys.  98 (1993) 5612,
C     eqs. (B6) and (B10).
C
C     For the modifications due to size adjustment, see:
C     A.D. Becke, J. Chem. Phys. 88, 2547 (1988);
C     O. Treutler and R. Ahlrichs, J. Chem. Phys. 102, 346 (1995).
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
C  xmuAB   -  hyperbolic coordinate
C  Gradx   -  on exit contains the 3 components of gradient
C             adjusted for size (only the gradient with respect
C             atom a is computed)
C  Hess    -  on exit contains the 15 components of the hessian
C             adjusted for size (the second derivatives involving
C             two times atom B are not computed)
C
C     Ordering of array Hess:
C
C                         1     XA  XA
C                         2     YA  XA
C                         3     YA  YA
C                         4     ZA  XA
C                         5     ZA  YA
C                         6     ZA  ZA
C                         7     XB  XA
C                         8     XB  YA
C                         9     XB  ZA
C                        10     YB  XA
C                        11     YB  YA
C                        12     YB  ZA
C                        13     ZB  XA
C                        14     ZB  YA
C                        15     ZB  ZA
C
      DIMENSION Gradx(3),Hess(15)
      Real*8  Grad(6)
      PARAMETER (One=1.0d0,Two=2.0d0,Three=3.0d0)
C
C  vectors from grid point to atoms A and B
C
      ax = (XA-X)
      ay = (YA-Y)
      az = (ZA-Z)
      bx = (XB-X)
      by = (YB-Y)
      bz = (ZB-Z)
c...
c...  vector lengths
c...
      A = SQRT(az*az+ay*ay+ax*ax)
      B = SQRT(bz*bz+by*by+bx*bx)
c...
c...  precomputes common quantities
c...
      abx = bx-ax
      aby = by-ay
      abz = bz-az
      arab= rab/A
      brab= rab/B
      ab3 = rab*rab*rab
      ambab3 = ab3*(A-B)
c...
c...  gradient components
c...
      GRAD(1) = ax*arab-ambab3*abx
      GRAD(2) = ay*arab-ambab3*aby
      GRAD(3) = az*arab-ambab3*abz
      GRAD(4) = -bx*brab+ambab3*abx
      GRAD(5) = -by*brab+ambab3*aby
      GRAD(6) = -bz*brab+ambab3*abz
c...
c...  modified for size adjustment
c...
      coef = (One - Two*xaAB*xmuAB)
      GRADX(1) = grad(1)*coef
      GRADX(2) = grad(2)*coef
      GRADX(3) = grad(3)*coef
c...
c...  additional common quantities for hessian
c...
      aab3=ab3/A
      bab3=ab3/B
      a3ab=rab/A**3
      tambab5=Three*ambab3*rab*rab
c...
c...  hessian components with size adjustment factors
c...
      TxaAB = Two*xaAB
      hess(1)= (abx*abx*tambab5-two*ax*abx*aab3-ambab3-ax*ax*a3ab+arab)
     1   *coef-TxaAB*grad(1)*grad(1)
      hess(2)= (abx*aby*tambab5-(ax*aby+ay*abx)*aab3-ax*ay*a3ab)
     1   *coef-TxaAB*grad(2)*grad(1)
      hess(3)=(aby*aby*tambab5-two*ay*aby*aab3-ambab3-ay*ay*a3ab+arab)
     1   *coef-TxaAB*grad(2)*grad(2)
      hess(4)= (abx*abz*tambab5-(ax*abz+az*abx)*aab3-ax*az*a3ab)
     1   *coef-TxaAB*grad(3)*grad(1)
      hess(5)= (aby*abz*tambab5-(ay*abz+az*aby)*aab3-ay*az*a3ab)
     1   *coef-TxaAB*grad(3)*grad(2)
      hess(6)= (abz*abz*tambab5-two*az*abz*aab3-ambab3-az*az*a3ab+arab)
     1   *coef-TxaAB*grad(3)*grad(3)
      hess(7) = (abx*bx*bab3-tambab5*abx*abx+ax*abx*aab3+ambab3)
     1   *coef-TxaAB*grad(4)*grad(1)
      hess(8) = (aby*bx*bab3-tambab5*abx*aby+ay*abx*aab3)
     1   *coef-TxaAB*grad(4)*grad(2)
      hess(9) = (abz*bx*bab3-tambab5*abx*abz+az*abx*aab3)
     1   *coef-TxaAB*grad(4)*grad(3)
      hess(10) = (abx*by*bab3-tambab5*abx*aby+ax*aby*aab3)
     1   *coef-TxaAB*grad(5)*grad(1)
      hess(11) = (aby*by*bab3-tambab5*aby*aby+ay*aby*aab3+ambab3)
     1   *coef-TxaAB*grad(5)*grad(2)
      hess(12) = (abz*by*bab3-tambab5*abz*aby+az*aby*aab3)
     1   *coef-TxaAB*grad(5)*grad(3)
      hess(13) = (abx*bz*bab3-tambab5*abx*abz+ax*abz*aab3)
     1   *coef-TxaAB*grad(6)*grad(1)
      hess(14) = (aby*bz*bab3-tambab5*aby*abz+ay*abz*aab3)
     1   *coef-TxaAB*grad(6)*grad(2)
      hess(15) = (abz*bz*bab3-tambab5*abz*abz+az*abz*aab3+ambab3)
     1   *coef-TxaAB*grad(6)*grad(3)
      RETURN
      END
c =====================================================================
c
      SUBROUTINE GWTONLY(NPP,    XGRID,  WGHT,   INDX,   IATOM,
     $                   NAtoms, XC,     AIJ,    rrij,   RDist,
     $                   P,      T,      xmuAB,  GWT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates derivative of weight at contributing grid points for
C  each atom and stores contributions in GWT array.
C
C  See:  Appendix B of Johnson, Gill and Pople
C        J.Chem.Phys.  98 (1993) 5612
C
C  MM December 2003 This is a modified version of <GradWT>, that
C  does not multiply by the exchange-correlation energy. It is
C  meant to be used in the Hessian code, during multiple passes
C  calculation of the derivative Fock matrices.
C
C  ARGUMENTS
C
C  NPP     -  number of contributing (non-zero) grid points this batch
C  XGRID   -  X,Y,Z grid points
C  WGHT    -  grid quadrature weights
C  INDX    -  index into contributing grid points
C  IATOM   -  which atom grid point formerly associated with
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
      DIMENSION XGRID(3,*),WGHT(*),INDX(NPP)
      REAL*8 XC(3,NAtoms),GWT(3,NAtoms,*)
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
      Wt = wght(IPP)
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
      GWT(1,IAtm,IPP) = GWT(1,IAtm,IPP) + coef*GradX
      GWT(2,IAtm,IPP) = GWT(2,IAtm,IPP) + coef*GradY
      GWT(3,IAtm,IPP) = GWT(3,IAtm,IPP) + coef*GradZ
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
       GWT(1,IAtm,IPP) = GWT(1,IAtm,IPP) + coef*GradX
       GWT(2,IAtm,IPP) = GWT(2,IAtm,IPP) + coef*GradY
       GWT(3,IAtm,IPP) = GWT(3,IAtm,IPP) + coef*GradZ
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
      SUBROUTINE GHWT(NPP,   XGRID,  WGHT,   INDX, IATOM,
     $                NAtoms,XC,     AIJ,    rrij, P,
     $                T,     Q,      gmuab,  hmuAB,Evec,
     $                GWT,    HWT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C    MM (01/09/2003).
C   The structure of this subroutine is derived from subroutine GradWT.
C
C  Computes first and second derivative of weight at contributing
C  grid points for each atom and stores contributions in GWT and HWT
C  arrays.
C
C  See:  Appendix B of Johnson, Gill and Pople
C        J.Chem.Phys.  98 (1993) 5612
C        and
C        Johnson and Frisch, J. Chem. Phys. 100, 7492 (1994).
C        (the actual working equations for the weight second
C        derivatives have been derived with maxima)
C
C  ARGUMENTS
C
C  NPP     -  number of contributing (non-zero) grid points this batch
C  XGRID   -  X,Y,Z grid points
C  WGHT    -  grid quadrature weights
C  INDX    -  index into contributing grid points
C  IATOM   -  which atom grid point formerly associated with
C  NAtoms  -  total number of atoms in system
C  XC      -  Cartesian coordinates of all atoms
C  AIJ     -  Becke atomic size-adjustment factors
C  rrij    -  inverse of interatomic distances
C  P       -  scratch array for cell functions
C  T       -  scratch array for auxiliary functions
C  Q       -  scratch array for auxiliary functions
C  Gmuab   -  scratch array for gradient of hyperbolic coordinates
C  HmuAB   -  scratch array for hessian of hyperbolic coordinates
C  Evec    -  exchange-correlation energy contribution per grid point
C  GWT     -  on exit contains weight first derivatives
C             over the grid points
C  HWT     -  on exit contains weight second derivatives over
C             the grid points
C
C
      DIMENSION XGRID(3,*),WGHT(*),INDX(NPP),
     $          GWT(3,Natoms,*),HWT(3,natoms,3,natoms,*),
     $          Evec(*)
      REAL*8 XC(3,NAtoms)
      REAL*8 AIJ(NAtoms,NAtoms),rrij(natoms,natoms)
      real*8 P(natoms),T(natoms,natoms),Q(natoms,natoms),
     $       gmuab(3,natoms,natoms),hmuab(15,natoms,natoms)
      real*8 gzb(3),gzc(3),hz(9)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0,ThreeHalf=3.0d0/2.0d0)
      PARAMETER (two=2.0d0)
      PARAMETER (CON1=-27.0d0/16.0d0,CON2=243.0d0/32.0d0)
      PARAMETER (CON3=81.0d0/16.0d0,CON4=27.0d0/8.0d0)
C
      data thrsh/1.0d-14/
      integer i12(3,3)
      data i12/1,2,4,2,3,5,4,5,6/
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
      Wt = wght(IPP)
      if(abs(wt).lt.thrsh) goto 90
c...
      X = XGRID(1,IPP)
      Y = XGRID(2,IPP)
      Z = XGRID(3,IPP)
C
C    loop over atom I
C
      ZF = Zero
      DO 30 IAtm=1,NAtoms
      PSum = One
      XI = X - XC(1,IAtm)
      YI = Y - XC(2,IAtm)
      ZI = Z - XC(3,IAtm)
      DistI = SQRT(XI*XI + YI*YI + ZI*ZI)
C
C   loop over atom J, skipping J=i
C
      DO 20 JAtm=1,NAtoms
      If(JAtm.EQ.IAtm) GO TO 20
      XJ = X - XC(1,JAtm)
      YJ = Y - XC(2,JAtm)
      ZJ = Z - XC(3,JAtm)
      DistJ = SQRT(XJ*XJ + YJ*YJ + ZJ*ZJ)
C
C   compute the hyperbolyc coordinate xmu and the
C   adjusted hyperbolic cordinate xmum for the current
C   atom pair
C
      xmu = (DistI-DistJ)*rrij(IAtm,JAtm)
      xmum = xmu + AIJ(IAtm,JAtm)*(One-xmu*xmu)
C
C  compute first and second derivatives of the adjusted
C  hyperbolyc coordinate and store them in the scratch
C  arrays gmuab and hmuab
C
      CALL GHMu(X,Y,Z,XC(1,JAtm),XC(2,JAtm),XC(3,JAtm),
     $       XC(1,IAtm),XC(2,IAtm),XC(3,IAtm),rrij(IAtm,JAtm),
     $       AIJ(IAtm,JAtm),xmu,gmuab(1,IAtm,JAtm),hmuab(1,IAtm,JAtm))

C
C  form the auxiliary function arrays T  (Eq.B9)
C  and Q. Q is defined as in Eq. B9, except that
C  s' is replaced by s''.
C

      p1 = ThreeHalf*xmum - Half*(xmum**3)
      p2 = ThreeHalf*p1   - Half*(p1**3)
      p3 = ThreeHalf*p2   - Half*(p2**3)
      s  = Half*(One-p3)
c
      IF(s.LT.thrsh) THEN
       T(IAtm,JAtm) = Zero
       Q(IAtm,JAtm) = Zero
      ELSE
       Om=(One-xmum*xmum)
       Op1=(One-p1*p1)
       Op2=(One-p2*p2)
       T(IAtm,JAtm) = CON1*Op2*Op1*Om/s
       Q(iAtm,jAtm) = (CON2*p2*Op1*Op1*Om*Om+CON3*Op2*p1*Om*Om+
     $                CON4*Op2*Op1*xmum)/s
      ENDIF
c
      PSum = PSum*s
c
 20   CONTINUE
C
C    compute the cell function
C
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
c...  gradient and hessian calculation.
c...
c...  Gradient (eq. B7):
c...
c...  d omega_a             1  d P(A)   P(A) d Z
c...  --------- (i.ne.a) = --- ------ - ---- ---
c...     d_i                Z   d_i     Z**2 d_i
c...
c...  Hessian:
c...
c...   d**2 omega_a                        1  d**2 P(A)
c...  ------------- (i.ne.a.and.j.ne.a) = --- --------- -
c...     d_i d_j                           Z   d_i d_j
c...
c...   1    d P(A) d Z    1    d P(A) d Z   2 P(A) d Z d Z
c...  ----  ------ --- - ----  ------ --- + ------ --- --- -
c...  Z**2   d_i   d_j   Z**2   d_j   d_i    Z**3  d_i d_j
c...
c...  P(A) d**2 Z
c...  ---- -------
c...  Z**2 d_i d_j
c...
c...  the derivatives of Z are computed in subroutines
c...  gradz, hbbz and hbcz. In the comment lines of these
c...  subroutines are also given the expressions for the
c...  derivatives of the cell functions P
c...
      oma=WGHT(IPP)
      omae=oma*Evec(IPP)
      omap=oma
      DO 60 IAtm=1,NATOMS
        If(IAtm.EQ.IATOM) GO TO 60
        qai=q(IATOM,iatm)
        tai=t(IATOM,iatm)
        call gradz(p,t,gmuab,gzb,iatm,natoms)
        call hbbz(p,t,q,gmuab,hmuab,hz,iatm,natoms)
c...
c...  gradient
c...
       gwt(1,iatm,IPP)=-omap*(tai*gmuab(1,iatm,iatom)+zff*gzb(1))
       gwt(2,iatm,IPP)=-omap*(tai*gmuab(2,iatm,iatom)+zff*gzb(2))
       gwt(3,iatm,IPP)=-omap*(tai*gmuab(3,iatm,iatom)+zff*gzb(3))
c...
c...  hessian with respect center iatm two times
c...
        do ic1=1,3
          do ic2=1,3
             hwt(ic1,iatm,ic2,iatm,IPP)=
     $      +omae*(qai*gmuab(ic1,iatm,iatom)*gmuab(ic2,iatm,iatom)
     $              -tai*hmuab(i12(ic1,ic2),iatm,iatom)
     $        +zff*tai*(gmuab(ic1,iatm,iatom)*gzb(ic2)
     $                        +gmuab(ic2,iatm,iatom)*gzb(ic1))
     $        +two*zff*zff*gzb(ic1)*gzb(ic2)
     $      -zff*hz(i12(ic1,ic2)))
          enddo
        enddo
        DO 50 JAtm=IAtm+1,NATOMS
          IF(JAtm.EQ.IATOM) goto 50
          call gradz(p,t,gmuab,gzc,jatm,natoms)
          call hbcz(p,t,q,gmuab,hmuab,hz,iatm,jatm,natoms)
          taj=t(IATOM,jatm)
c...
c...  hessian with respect centers iatm and jatm
c...
          index=0
          do ic1=1,3
            do ic2=1,3
            index=index+1
              hwt(ic1,iatm,ic2,jatm,IPP)=
     $       +omae*(tai*taj*gmuab(ic1,iatm,iatom)*gmuab(ic2,jatm,iatom)
     $        +zff*taj*gmuab(ic2,jatm,iatom)*gzb(ic1)-
     $         zff*hz(index)+
     $         zff*tai*gmuab(ic1,iatm,iatom)*gzc(ic2)+
     $         two*zff*zff*gzb(ic1)*gzc(ic2))
            enddo
          enddo
 50     CONTINUE
 60   CONTINUE
c -- end loop over grid points
 90   CONTINUE
C
      RETURN
      END
c======================================================================
      subroutine hbbz(p,t,q,gmuab,hmuab,hz,ib,natoms)
      IMPLICIT REAL*8(A-H,O-Z)
c...
c...  MM (01/09/2003)
c...
c...  computes the second derivatives (with respect atom ib
c...  two times), of Z, defined as the sum of the
c...  cell functions P, see eq B2 of Johnson, Gill and Pople,
c...  J. Chem Phys, 98, 5612 (1993).
c...
c...    d**2 Z     d**2 P(ib)                 d**2 P(i)
c...  ---------  = ---------- + sum_(i.ne.ib) ---------
c...  (d_ib)**2    (d_ib)**2                  (d_ib)**2
c...
c...  where the sum is over all centers except ib.
c...
c...  P(i) is the cell function of center i (eq. B3).
c...
c...  Developing further, one obtains:
c...
c...  d**2 P(ib)
c...  ---------- = P(ib) sum_(i.ne.ib) [ q(ib,i) gmuab(ib,i)**2 +
c...  (d_ib)**2
c...                  t(ib,i) hmuab(ib,i)-(t(ib,i) gmuab(ib,i)**2 ] +
c...
c...             + P(ib) [sum_(i.ne.ib) t(ib,i) gmuab(ib,i)] **2
c...
c...  and
c...
c...  d**2 P(i)
c...  --------- (i.ne.ib) = P(i) [ q(i,ib) gmuab(ib,i)**2
c...  (d_ib)**2
c...                              - t(i,ib) hmuab(ib,i) ]
c...
c...   t(i,ib)    auxiliary function         (eq  B9)
c...   q(i,ib)    auxiliary function (the same as t, but involving
c...              the second derivative of s)
c...   gmuab      gradient of adjusted hyperbolic coordinates
c...              (see subroutine ghmu)
c...   hmuab      hessian of adjusted hyperbolic coordinates
c...              (see subroutine ghmu)
c...
      real*8 hz(6)
      real*8 P(natoms),T(natoms,natoms),Q(natoms,natoms),
     $       gmuab(3,natoms,natoms),hmuab(15,natoms,natoms)
      real*8 psum(3)
      parameter (zero=0.0d0)
c...
      psum(1)=zero
      psum(2)=zero
      psum(3)=zero
      do ic=1,6
      hz(ic)=zero
      enddo
c...
c...              d**2 P(ib)
c... first part:  ----------
c...              (d_ib)**2
c...
      do 10 i=1,natoms
        if(i.eq.ib) goto 10
        qbi=q(ib,i)
        tbi=t(ib,i)
        index=0
        do ic1=1,3
          do ic2=1,ic1
            index=index+1
            hz(index)=hz(index)+(qbi-tbi*tbi)*gmuab(ic1,ib,i)
     $              *gmuab(ic2,ib,i)+tbi*hmuab(index,ib,i)
          enddo
          psum(ic1)=psum(ic1)+tbi*gmuab(ic1,ib,i)
        enddo
10    continue
      index=0
      do ic1=1,3
        do ic2=1,ic1
           index=index+1
           hz(index)=p(ib)*(hz(index)+psum(ic1)*psum(ic2))
        enddo
      enddo
c...
c...                            d**2 P(i)
c... second part: sum_(i.ne.ib) ---------
c...                            (d_ib)**2
c...
      do 20 i=1,natoms
        if(i.eq.ib) goto 20
        qib=q(i,ib)
        tib=t(i,ib)
        index=0
        do ic1=1,3
          do ic2=1,ic1
            index=index+1
           hz(index)=hz(index)+p(i)*(qib*gmuab(ic1,ib,i)*gmuab(ic2,ib,i)
     $                             -tib*hmuab(index,ib,i))
          enddo
        enddo
20    continue
      return
      END
c======================================================================
      subroutine gradz(p,t,gmuab,grzet,ib,natoms)
      IMPLICIT REAL*8(A-H,O-Z)
c...
c...  MM (01/09/2003)
c...
c...  computes the gradient of Z, defined as the sum of the
c...  cell functions P, see eq B2 of Johnson, Gill and Pople,
c...  J. Chem Phys, 98, 5612 (1993),
c...  with respect center ib:
c...
c...  gz(ic)=sum_(i.ne.ib) (p(ib)*t(ib,i)-p(i)*t(i,ib))*gmuab(ic,ib,i)
c...
c...  where the summation is over all centers except ib.
c...
c...   p(i)       cell function of center i  (eq. B3)
c...   t(i,ib)    auxiliary function         (eq  B9)
c...   gmuab      gradient of adjusted hyperbolic coordinates
c...              (see subroutine ghmu)
c...
c...
      real*8 grzet(3)
      real*8 P(natoms),T(natoms,natoms),gmuab(3,natoms,natoms)
      parameter (zero=0.0d0)
      grzet(1)=zero
      grzet(2)=zero
      grzet(3)=zero
      do 10 i=1,natoms
        if(i.eq.ib) goto 10
        piti=p(ib)*t(ib,i)-p(i)*t(i,ib)
        grzet(1)=grzet(1)+piti*gmuab(1,ib,i)
        grzet(2)=grzet(2)+piti*gmuab(2,ib,i)
        grzet(3)=grzet(3)+piti*gmuab(3,ib,i)
10    continue
      return
      END
c======================================================================
      subroutine cpartd(p,t,q,gmuab,hmuab,d,ib,ic,natoms)
      IMPLICIT REAL*8(A-H,O-Z)
c...
c...  MM (01/09/2003)
c...
c...  computes a partial sum needed in the calculation
c...  of the weight hessian. See the comments in subroutine
c...  hbcz for details
c...
      real*8 d(9)
      real*8 P(natoms),T(natoms,natoms),Q(natoms,natoms),
     $       gmuab(3,natoms,natoms),hmuab(15,natoms,natoms)
      real*8 psum(3)
      parameter (zero=0.0d0)
      tbc=t(ib,ic)
      qbc=q(ib,ic)
      psum(1)=zero
      psum(2)=zero
      psum(3)=zero
      do 10 i=1,natoms
        if(i.eq.ib) goto 10
        if(i.eq.ic) goto 10
        tbi=t(ib,i)
        psum(1)=psum(1)+tbi*gmuab(1,ib,i)
        psum(2)=psum(2)+tbi*gmuab(2,ib,i)
        psum(3)=psum(3)+tbi*gmuab(3,ib,i)
10    continue
c...
c...  this is done with ic2 in the external loop,
c...  so the hmuab array is properly addressed
c...
      index=0
      do ic2=1,3
        do ic1=1,3
          index=index+1
          d(index)=p(ib)*(tbc*hmuab(index+6,ib,ic)-qbc*gmuab(ic1,ib,ic)*
     $             gmuab(ic2,ic,ib)-tbc*gmuab(ic2,ic,ib)*psum(ic1))
        enddo
      enddo
      return
      END
c======================================================================
      subroutine hbcz(p,t,q,gmuab,hmuab,hz,ib,ic,natoms)
      IMPLICIT REAL*8(A-H,O-Z)
c...
c...  MM (01/09/2003)
c...
c...  computes the second derivatives (with respect atoms ib
c...  and ic), of Z, defined as the sum of the
c...  cell functions P, see eq B2 of Johnson, Gill and Pople,
c...  J. Chem Phys, 98, 5612 (1993).
c...
c...   d**2 Z     d**2 P(ib)   d**2 P(ic)
c...  --------- = ---------- + ---------- +
c...  d_ib d_ic   d_ib d_ic    d_ib d_ic
c...
c...                                             d**2 P(i)
c...                 + sum_(i.ne.ib.and.i.ne.ic) ---------
c...                                             d_ib d_ic
c...
c...  where the sum is over all centers except ib and ic.
c...
c...  P(i) is the cell function of center i (eq. B3).
c...
c...  Developing further, one obtains:
c...
c...  d**2 P(ib)
c...  ---------- = P(ib) [ t(ib,ic) hmuab(ib,ic) -
c...  d_ib d_ic
c...                    q(ib,ic) gmuab(ib,ic) gmuab(ic,ib) -
c...
c...       t(ib) gmuab(ic,ib) sum_(i.ne.ib.and.i.ne.ic) t(ib,i) gmuab(ib,i) ]
c...
c...  and
c...
c...  d**2 P(i)
c...  --------- (i.ne.ib.and.i.ne.ic) =
c...  d_ib d_ic
c...
c...                       = P(i) t(i,ib) t(i,ic) gmuab(ib,i) gmuab(ic,i)
c...
c...   t(i,ib)    auxiliary function         (eq  B9)
c...   q(i,ib)    auxiliary function (the same as t, but involving
c...              the second derivative of s)
c...   gmuab      gradient of adjusted hyperbolic coordinates
c...              (see subroutine ghmu)
c...   hmuab      hessian of adjusted hyperbolic coordinates
c...              (see subroutine ghmu)
c...
      real*8 hz(9),dcbb(9),dbcc(9)
      real*8 P(natoms),T(natoms,natoms),Q(natoms,natoms),
     $       gmuab(3,natoms,natoms),hmuab(15,natoms,natoms)
      parameter (zero=0.0d0)
c...
c...  d**2 P(ib)
c...  ----------
c...  d_ib d_ic
c...
c...  actually, this call to cpartd returns the derivative
c...  with respect ic as first index and ib as second. Thus
c...  the actual contribution is given by the transposed of
c...  the dcbb array
c...
      call cpartd(p,t,q,gmuab,hmuab,dcbb,ib,ic,natoms)
c...
c...  d**2 P(ic)
c...  ----------
c...  d_ib d_ic
c...
      call cpartd(p,t,q,gmuab,hmuab,dbcc,ic,ib,natoms)
c...
      hz(1)=dcbb(1)+dbcc(1)
      hz(2)=dcbb(4)+dbcc(2)
      hz(3)=dcbb(7)+dbcc(3)
      hz(4)=dcbb(2)+dbcc(4)
      hz(5)=dcbb(5)+dbcc(5)
      hz(6)=dcbb(8)+dbcc(6)
      hz(7)=dcbb(3)+dbcc(7)
      hz(8)=dcbb(6)+dbcc(8)
      hz(9)=dcbb(9)+dbcc(9)
c...
c...                           d**2 P(i)
c... sum_(i.ne.ib.and.i.ne.ic) ---------
c...                           d_ib d_ic
c...
      do 10 i=1,natoms
        if(i.eq.ib) goto 10
        if(i.eq.ic) goto 10
        ptt=p(i)*t(i,ib)*t(i,ic)
        index=0
        do ic1=1,3
          do ic2=1,3
            index=index+1
            hz(index)=hz(index)+ptt*gmuab(ic1,ib,i)*gmuab(ic2,ic,i)
          enddo
        enddo
10    continue
      return
      end
