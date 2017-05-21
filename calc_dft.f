      subroutine setup_dft(dft,    ityp,   nfock,  nrad,   nang,
     $                     NBatch, IAN,    XNuc)

      use memory

      implicit real*8(a-h,o-z)
c
c  Reserve memory and recompute the (constant) data needed
c  for the DFT routines and put into the depository
c  Additionally puts all dft pointers in common /dftpnt/
c
c  ARGUMENTS
c
c  dft    -  DFT flag
c  ityp    -  what type of DFT is it?
c              0 - closed-shell energy
c              1 - open-shell energy
c              2 - closed-shell gradient
c              3 - open-shell gradient
c              4 - closed-shell NMR (**WARNING** need to change NMR code)
c              5 - closed-shell Hessian
c              6 - closed-shell Hessian CPHF
c              7 - open-shell Hessian
c              8 - open-shell Hessian CPHF
c  nfock   -  number of simultaneous Fock matrices constructed
c               (for multifock procedure - normally 1)
c  nrad    -  maximum number of radial grid points on largest atom
c  nang    -   ditto for angular grid points
c  NBatch  -  maximum number of grid points per batch
c  IAN     -  atomic charges (numbers)
c  XNuc    -  nuclear coordinates
c
c     common /big/ bl(300)
c     common /intbl/ maxsh,inx(100)
      common /dftpnt/ ird,iaij,idst,iixx,iiwt,ipre,iexp,
     $                ixgd,iwt,iscr,icofm,nscr
      Dimension IAN(*),XNuc(3,*)
      Integer dft
c
c -- get some quantities from the depository
      call getival('na  ',na)              ! number of atoms
      call getival('ndum',ndum)            ! number of dummy atoms
      call getival('nsh ',nsh)             ! number of primitive shells
      call getival('ncs ',ncs)             ! number of contracted shells
      call getival('ncf ',ncf)             ! number of basis functions
      call getival('ibas',ibas)            ! pointer to basis set data
      call getival('ictr',ictr)            ! pointer to integer basis set data
      call getival('nsym',nsym)            ! number of abelian operations
      call getival('nmo ',nmo)             ! number of occupied MOs
c
      na = na-ndum                         ! real atoms only
c
c -- reserve memory in bl array
      call getmem(na*na,ird)               ! for grid
      call getmem(na*na,iaij)              ! for grid
      call getmem(na,idst)
      call getmem(3*1130,iixx)             ! precomputed angular grids
      call getmem(1130,iiwt)               ! precomputed angular weights
      call getmem(nsh,ipre)                ! precomputed normalization
      call getmem(ncs,iexp)                ! mimimum exponents per shell
c
c -- grid (coordinates and weights)
      mxpts = nrad*nang
      call getmem(mxpts*3,ixgd)            ! grid itself
      call getmem(mxpts,iwt)               ! grid weights
c
c -- general scratch storage
      If(ityp.le.1) Then
        nscr = 2*NBatch*ncf + 5*NBatch + 1
        If(ityp.eq.1) nscr = nscr + 2*NBatch
        If(dft.GT.3) Then
          nscr = nscr + 3*NBatch*ncf + 6*NBatch
          If(ityp.eq.1) nscr = nscr + 6*NBatch
        EndIf
      Else If(ityp.eq.2) Then
        nscr = 4*NBatch*ncf +
     $            MAX(NBatch*ncf + 5*NBatch, na*(na+3)) + 1
        If(dft.GT.3) nscr = nscr + 6*NBatch*ncf + 6*NBatch
      Else If(ityp.eq.3) Then
        nscr = 4*NBatch*ncf +
     $            MAX(NBatch*ncf + 7*NBatch, na*(na+3)) + 1
        If(dft.GT.3) nscr = nscr + 6*NBatch*ncf + 12*NBatch
      Else If(ityp.eq.4) Then
c -- for NMR (**WARNING** need to change NMR code)
        nscr = 2*NBatch*ncf + NBatch*nmo + 4*NBatch + 1
        If(dft.GT.3) nscr = nscr + 3*NBatch*ncf + 6*NBatch + 3*nmo
c -- add extra in case Malkin
        nscr = nscr + (ncf-nmo)*NBatch + NBatch
      Else If(ityp.eq.5) Then
c -- for closed-shell Hessian
        nscr = 11*ncf*NBatch + 6*NBatch + 1 + 3*na*(3*na+1)
        If(dft.GT.3) nscr = nscr + 10*ncf*NBatch + 13*NBatch +
     $                       24*na + 36*na**2
        call getival('wghtd',IdWt)
        If(IdWt.EQ.1) nscr = nscr + (3*na*(3*na+1)+2)*NBatch +
     $                        20*na**2 + na
c -- make sure there is enough scratch memory for derivative Fock buffer
        nscr = MAX(nscr,3*(ncf*(ncf+1))/2)
      Else If(ityp.eq.6) Then
c -- for closed-shell Hessian CPHF
cc        nscr = 6*ncf*NBatch + 6*NBatch + 1 + 3*na*NBatch
        nscr = 2*ncf*NBatch + 6*NBatch + 1 + 3*na
cc        If(dft.GT.3) nscr = nscr + 3*ncf*NBatch + 6*NBatch + 24*na
        If(dft.GT.3) nscr = nscr + 3*ncf*NBatch + 13*NBatch + 24*na
      Else If(ityp.eq.7) Then
c -- for open-shell Hessian
        nscr = 11*ncf*NBatch + 10*NBatch + 1 + 18*na**2 + 12*na
        If(dft.GT.3) nscr = nscr + 10*ncf*NBatch + 21*NBatch + 54*na +
     $                             3*27*na**2
        call getival('wghtd',IdWt)
        If(IdWt.EQ.1) nscr = nscr + (3*na*(3*na+1)+2)*NBatch +
     $                        20*na**2 + na
c -- make sure there is enough scratch memory for derivative Fock buffer
        nscr = MAX(nscr,3*(ncf*(ncf+1))/2)
      Else If(ityp.eq.8) Then
c -- for open-shell Hessian CPHF
        nscr = 2*ncf*NBatch + 10*NBatch + 1 + 12*na
        If(dft.GT.3) nscr = nscr + 3*ncf*NBatch + 21*NBatch + 45*na
      EndIf
cc      If(ityp.ne.1) call getmem(nscr,iscr)  ! don't allocate for UDFT
      If(ityp.gt.1) call getmem(nscr,iscr)  ! don't allocate for energy
c
c -- need to save occupied MOs if multifock
      if(nfock.gt.1) call getmem(nmo*ncf*(nfock-1),icofm)
c
c -- calculate inverse atomic distances, Becke aij parameters
c -- and angular grid and weights prior to full grid construction
      call PreGRID(na,    XNuc,   IAN,     bl(ixgd),bl(iwt),
     $           bl(idst),bl(ird),bl(iaij),bl(iixx),bl(iiwt))
c
c -- get array of smallest exponent per shell
      call GetEXP(nsh,ncs,bl(ibas),bl(ictr),bl(iexp))
c
c -- precompute shell normalization factors
      call AOInit(ncs,bl(ibas),bl(ictr),bl(ipre))
c
      return
      end
c ===============================================================
      SUBROUTINE GetEXP(nsh,ncs,BASDAT,INX,ExpMIN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Extracts and stores smallest exponent for each shell
C
C  ARGUMENTS
C
C  nsh     -  number of primitive shells
C  ncs     -  number of contracted shells
C  BASDAT  -  basis set data (TEXAS format)
C  INX     -    ditto
C  ExpMIN  -  on exit, contains smallest exponent per shell
C
C
      DIMENSION BASDAT(13,nsh),ExpMIN(ncs)
      INTEGER INX(12,ncs)
      parameter (huge=1.0d+25)
C
C
cc      write(6,*) ' In <GetEXP>'
cc      write(6,*) ' BASDAT array is:'
cc      do i=1,13
cc      write(6,*) (basdat(i,j),j=1,ncs)
cc      enddo
cc      write(6,*) ' INX array is'
cc      do i=1,12
cc      write(6,1111) (inx(i,j),j=1,ncs)
cc 1111 format(12i4)
cc      enddo
      DO 20 ics=1,ncs
        ifirstsh = INX(1,ics)+1
        ilastsh  = INX(5,ics)
c
        expm = huge
c
        DO 10 ish=ifirstsh,ilastsh      ! loop over primitive shells
          expon = BASDAT(1,ish)
          If(expon.LT.expm) expm = expon
 10     CONTINUE
c
        ExpMIN(ics) = expm
c
 20   CONTINUE
C
      RETURN
      END
c ====================================================================
      SUBROUTINE MakeXC(dft,    NPP,    NBas,   nbf,    thrsh,
     $                  WGHT,   Pot,    PotX,   VAO,    VAOX,
     $                  INB,    INDX,   XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into Exchange-Correlation Fock matrix
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All methods include Slater exchange)
C              1 - 3 - local correlation
C             >3 - nonlocal - need AO derivatives
C  NPP     -  number of contributing (non-zero) grid points this batch
C  NBas    -  total number of basis functions
C  nbf     -  indexing array to location of "non-zero" basis functions
C  thrsh   -  threshold for neglect of contribution
C  WGHT    -  grid quadrature weights
C  Pot     -  potential
C  PotX    -  potential derivatives
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  basis function derivatives at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  INDX    -  index into contributing columns of VAO
C  XC      -  exchange-correlation Fock matrix
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),Pot(NPP),
     $          PotX(3,NPP),XC(*)
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
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
          JT = (JJ*(JJ-1))/2
          DO 10 I=nbf(IPP)+1,J
          II = INB(I)
          IJ = JT + II
          XC(IJ) = XC(IJ) + Val*VAO(I)
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
          JT = (JJ*(JJ-1))/2
          DO 50 I=nbf(IPP)+1,J
          II = INB(I)
          IJ = JT + II
          XC(IJ) = XC(IJ) + VAO(I)*ValJA +
     $     (VAOX(1,I)*VAx + VAOX(2,I)*VAy + VAOX(3,I)*VAz)*ValJ
 50       CONTINUE
        ENDIF
 60     CONTINUE
 70     CONTINUE
cc
      ENDIF
C
      RETURN
      END
c ==================================================================
      SUBROUTINE MakeXU(dft,    NPP,    NBas,   nbf,    thrsh,
     $                  WGHT,   PotA,   PotB,   PotAX,  PotBX,
     $                  VAO,    VAOX,   INB,    INDX,   XCA,  XCB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into Exchange-Correlation Fock matrix
C  ** UNRESTRICTED OPEN SHELL **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All methods include Slater exchange)
C              1 - 3 - local correlation
C             >3 - nonlocal - need AO derivatives
C  NPP     -  number of contributing (non-zero) grid points this batch
C  NBas    -  total number of basis functions
C  nbf     -  indexing array to location of "non-zero" basis functions
C  thrsh   -  threshold for neglect of contribution
C  WGHT    -  grid quadrature weights
C  PotA    -  alpha potential
C  PotB    -  beta potential
C  PotAX   -  alpha potential derivatives
C  PotBX   -  beta potential derivatives
C  VAO     -  basis function values at grid point
C  VAOX    -  basis function derivatives at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  INDX    -  index into contributing columns of VAO
C  XCA     -  alpha exchange-correlation Fock matrix
C  XCB     -  beta exchange-correlation Fock matrix
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),PotA(NPP),PotB(NPP),
     $          PotAX(3,NPP),PotBX(3,NPP),XCA(*),XCB(*)
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
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
        ValA = VAO(J)*VXCA
        ValB = VAO(J)*VXCB
        IF(Abs(ValA).GT.thrsh.OR.Abs(ValB).GT.thrsh) THEN
          JT = (JJ*(JJ-1))/2
          DO 10 I=nbf(IPP)+1,J
          II = INB(I)
          IJ = JT + II
          XCA(IJ) = XCA(IJ) + ValA*VAO(I)
          XCB(IJ) = XCB(IJ) + ValB*VAO(I)
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
          JT = (JJ*(JJ-1))/2
          DO 50 I=nbf(IPP)+1,J
          II = INB(I)
          IJ = JT + II
          XCA(IJ) = XCA(IJ) + VAO(I)*ValJA +
     $     (VAOX(1,I)*VAx + VAOX(2,I)*VAy + VAOX(3,I)*VAz)*ValJ
          XCB(IJ) = XCB(IJ) + VAO(I)*ValJB +
     $     (VAOX(1,I)*VBx + VAOX(2,I)*VBy + VAOX(3,I)*VBz)*ValJ
 50       CONTINUE
        ENDIF
 60     CONTINUE
 70     CONTINUE
      ENDIF
C
      RETURN
      END
c ==================================================================
      SUBROUTINE DENSort(NP,thrsh,DEN,INDX,NPP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts this batch of densities, determining which grid points
C  will be retained for the numerical integration
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  thrsh   -  threshold for neglecting contribution
C  DEN     -  density at each grid point
C
C  on exit
C
C  INDX    -  index into contributing columns of VAO
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION DEN(NP),INDX(NP)
C
C
      NPP = 0
c
      DO 10 IP=1,NP
      IF(DEN(IP).GT.thrsh) THEN
        NPP = NPP+1
        INDX(NPP) = IP
        DEN(NPP) = DEN(IP)
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c ==================================================================
      SUBROUTINE DENSortU(NP,thrsh,DENA,DENB,INDX,NPP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts this batch of densities, determining which grid points
C  will be retained for the numerical integration
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  thrsh   -  threshold for neglecting contribution
C  DENA    -  alpha density at each grid point
C  DENB    -  beta density at each grid point
C
C  on exit
C
C  INDX    -  index into contributing columns of VAO
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION DENA(NP),DENB(NP),INDX(NP)
C
C
      NPP = 0
c
      DO 10 IP=1,NP
      IF(DENA(IP)+DENB(IP).GT.thrsh) THEN
        NPP = NPP+1
        INDX(NPP) = IP
        DENA(NPP) = DENA(IP)
        DENB(NPP) = DENB(IP)
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c ====================================================================
      SUBROUTINE GetDENDeriv(NP,     NBas,   thrsh,  nbf,    NOcc,
     $                       CMO,    VAOX,   INB,    DEN,    DVT,
     $                       DVTX,   DenX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density derivatives at the current batch of grid points
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  nbf     -  number of "non-zero" basis functions
C  NOcc    -  number of occupied MOs
C  CMO     -  occupied MO coefficients (transposed)
C  VAOX    -  basis function derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAOX
C  DEN     -  density values at grid points
C  DVT     -  array for storing MO values and partial derivatives
C             (MO values already calculated in <GetDENS>)
C  DVTX    -  storage for partial density derivatives
C  DenX    -  on exit contains density derivatives
C
C
      DIMENSION CMO(NOcc,NBas),VAOX(3,*),DVT(NOcc,*),nbf(NP+1),
     $          DVTX(NOcc,3),DenX(3,NP),INB(NBas*NP),DEN(NP)
C
      PARAMETER (Zero=0.0d0)
C
C
      DO 40 IP=1,NP
c
      If(DEN(IP).LT.thrsh) GO TO 40
c
      CALL ZeroIT(DVTX,3*NOcc)
      DO 30 IC=1,3
      DO 20 J=nbf(IP)+1,nbf(IP+1)
      JJ = INB(J)
      Val = VAOX(IC,J)
      DO 19 I=1,NOcc
      DVTX(I,IC) = DVTX(I,IC) + Val*CMO(I,JJ)
 19   CONTINUE
cc      call daxpy(NOcc,Val,CMO(1,JJ),1,DVTX(1,IC),1)
 20   CONTINUE
 30   CONTINUE
c
      DVx = DDOT(NOcc,DVT(1,IP),1,DVTX(1,1),1)
      DVy = DDOT(NOcc,DVT(1,IP),1,DVTX(1,2),1)
      DVz = DDOT(NOcc,DVT(1,IP),1,DVTX(1,3),1)
c
c  density derivatives
c  (closed-shell should be doubled, but taken into account elsewhere)
      DenX(1,IP) = DVx + DVx
      DenX(2,IP) = DVy + DVy
      DenX(3,IP) = DVz + DVz
c
 40   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENS(NP,     NBas,   nbf,    NOcc,   CMO,
     $                   VAO,    INB,    DVT,    DEN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NBas    -  total number of basis functions
C  nbf     -  index array to indices of "non-zero" basis functions
C  NOcc    -  number of occupied MOs
C  CMO     -  occupied MO coefficients (transposed)
C  VAO     -  "non-zero" basis function values at grid point
C  INB     -  indices corresponding to entries in VAO
C  DVT     -  array for storing MO values
C  DEN     -  on exit contains density values per grid point
C
C
      REAL*8 CMO(Nocc,NBas),VAO(NBas*NP),DVT(NOcc,NP),DEN(NP)
      INTEGER INB(NBas*NP),nbf(NP+1)
C
C
      CALL ZeroIT(DVT,NP*NOcc)
      DO 30 IP=1,NP
      DO 20 J=nbf(IP)+1,nbf(IP+1)
      JJ = INB(J)
      Val = VAO(J)
      DO 19 I=1,NOcc
      DVT(I,IP) = DVT(I,IP) + Val*CMO(I,JJ)
 19   CONTINUE
cc      call daxpy(NOcc,Val,CMO(1,JJ),1,DVT(1,IP),1)
 20   CONTINUE
      DEN(IP) = DDOT(NOcc,DVT(1,IP),1,DVT(1,IP),1)
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENS1(NP,     NBas,   thrsh,  nbf,    DA,
     $                    DM,     VAO,    INB,    VM,     DEN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect
C  nbf     -  index array to indices of "non-zero" basis functions
C  DA      -  density matrix (lower triangle)
C  DM      -  array containing maximum density matrix element per column
C  VAO     -  "non-zero" basis function values at grid point
C  INB     -  indices corresponding to entries in VAO
C  VM      -  array containing maximum magnitude AO per grid point
C  DEN     -  on exit contains density values per grid point
C
C
      REAL*8 DA(NBas*(NBas+1)/2),DM(NBas),VAO(NBas*NP),VM(NP),DEN(NP)
      INTEGER INB(NBas*NP),nbf(NP+1)
C
      Parameter (Zero=0.0d0,Two=2.0d0)
C
C
      DO 30 IP=1,NP
      DVal = Zero
      VMx = VM(IP)
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      If(DM(II)*VMx*Abs(Val).LT.thrsh) GO TO 20
      DO 19 J=nbf(IP)+1,I-1
      JJ = INB(J)
      IJ = IT + JJ
      DVal = DVal + Two*DA(IJ)*Val*VAO(J)
 19   CONTINUE
      DVal = DVal + DA(IT+II)*Val**2
 20   CONTINUE
      DEN(IP) = DVal
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENSF(NP,     NBas,   nbf,    DA,     VAO,
     $                    INB,    DEN)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density at the current grid point
C  ** Same <GetDENS1> but without threshold testing **
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NBas    -  total number of basis functions
C  nbf     -  index array to indices of "non-zero" basis functions
C  DA      -  density matrix (lower triangle)
C  VAO     -  "non-zero" basis function values at grid point
C  INB     -  indices corresponding to entries in VAO
C  DEN     -  on exit contains density values per grid point
C
C
      REAL*8 DA(NBas*(NBas+1)/2),VAO(NBas*NP),DEN(NP)
      INTEGER INB(NBas*NP),nbf(NP+1)
C
      Parameter (Zero=0.0d0,Two=2.0d0)
C
C
      DO 30 IP=1,NP
      DVal = Zero
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      DO 19 J=nbf(IP)+1,I-1
      JJ = INB(J)
      IJ = IT + JJ
      DVal = DVal + Two*DA(IJ)*Val*VAO(J)
 19   CONTINUE
      DVal = DVal + DA(IT+II)*Val**2
 20   CONTINUE
      DEN(IP) = DVal
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENSU(NP,     NBas,   nbf,    DA,     DB,
     $                    VAO,    INB,    DENA,   DENB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NBas    -  total number of basis functions
C  nbf     -  index array to indices of "non-zero" basis functions
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta  density matrix (lower triangle)
C  VAO     -  "non-zero" basis function values at grid point
C  INB     -  indices corresponding to entries in VAO
C  DENA    -  on exit contains alpha density values per grid point
C  DENB    -   ditto beta density values
C
C
      REAL*8 DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),VAO(NBas*NP),
     $       DENA(NP),DENB(NP)
      INTEGER INB(NBas*NP),nbf(NP+1)
C
      Parameter (Zero=0.0d0,Two=2.0d0)
C
C
      DO 30 IP=1,NP
      DENA(IP) = Zero
      DENB(IP) = Zero
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      DO 19 J=nbf(IP)+1,I-1
      JJ = INB(J)
      IJ = IT + JJ
      DENA(IP) = DENA(IP) + Two*DA(IJ)*Val*VAO(J)
      DENB(IP) = DENB(IP) + Two*DB(IJ)*Val*VAO(J)
 19   CONTINUE
      DENA(IP) = DENA(IP) + DA(IT+II)*Val**2
      DENB(IP) = DENB(IP) + DB(IT+II)*Val**2
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENS1U(NP,     NBas,   thrsh,  nbf,    DA,
     $                     DB,     DM,     VAO,    INB,    VM,
     $                     DENA,   DENB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points in this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect
C  nbf     -  index array to indices of "non-zero" basis functions
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta  density matrix (lower triangle)
C  DM      -  array containing maximum density matrix element per column
C  VAO     -  "non-zero" basis function values at grid point
C  INB     -  indices corresponding to entries in VAO
C  VM      -  array containing maximum magnitude AO per grid point
C  DENA    -  on exit contains alpha density values per grid point
C  DENB    -   ditto beta density values
C
C
      REAL*8 DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NP),
     $       VAO(NBas*NP),VM(NP),DENA(NP),DENB(NP)
      INTEGER INB(NBas*NP),nbf(NP+1)
C
      Parameter (Zero=0.0d0,Two=2.0d0)
C
C
      DO 30 IP=1,NP
      DValA = Zero
      DValB = Zero
      VMx = VM(IP)
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      If(DM(II)*VMx*Abs(Val).LT.thrsh) GO TO 20
      DO 19 J=nbf(IP)+1,I-1
      JJ = INB(J)
      IJ = IT + JJ
      DValA = DValA + Two*DA(IJ)*Val*VAO(J)
      DValB = DValB + Two*DB(IJ)*Val*VAO(J)
 19   CONTINUE
      DValA = DValA + DA(IT+II)*Val**2
      DValB = DValB + DB(IT+II)*Val**2
 20   CONTINUE
      DENA(IP) = DValA
      DENB(IP) = DValB
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENDeriv1(NP,     NBas,   thrsh,  nbf,    DA,
     $                        DM,     VAO,    VAOX,   INB,    VMX,
     $                        DEN,    DenX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density derivatives at the current batch of grid points
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  nbf     -  number of "non-zero" basis functions
C  DA      -  density matrix (lower triangle)
C  DM      -  array containing maximum density matrix element per column
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  basis function derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAOX
C  VMX     -  array containing maximum magnitude AO derivative per grid point
C  DEN     -  density values at grid points
C  DenX    -  on exit contains density derivatives
C
C
      DIMENSION DA(NBas*(NBas+1)/2),DM(NBas),VAO(*),VAOX(3,*),nbf(NP+1),
     $          INB(NBas*NP),VMX(NP),DEN(NP),DenX(3,NP)
C
      PARAMETER (Zero=0.0d0)
C
C
      DO 30 IP=1,NP
c
      Den1 = Zero
      Den2 = Zero
      Den3 = Zero
      If(Abs(DEN(IP)).LT.thrsh) GO TO 29
c
      VMxx = VMX(IP)
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      If(DM(II)*VMxx*Abs(Val).LT.thrsh) GO TO 20
      DO 19 J=nbf(IP)+1,I
      JJ = INB(J)
      IJ = IT+JJ
      DAIJ = DA(IJ)*Val
c  (closed-shell should be doubled, but taken into account elsewhere)
      Den1 = Den1 + DAIJ*VAOX(1,J)
      Den2 = Den2 + DAIJ*VAOX(2,J)
      Den3 = Den3 + DAIJ*VAOX(3,J)
 19   CONTINUE
      DO 18 J=I+1,nbf(IP+1)
      JJ = INB(J)
      IJ = (JJ*(JJ-1))/2 + II
      DAIJ = DA(IJ)*Val
c  (closed-shell should be doubled, but taken into account elsewhere)
      Den1 = Den1 + DAIJ*VAOX(1,J)
      Den2 = Den2 + DAIJ*VAOX(2,J)
      Den3 = Den3 + DAIJ*VAOX(3,J)
 18   CONTINUE
 20   CONTINUE
c
 29   CONTINUE
      DenX(1,IP) = Den1
      DenX(2,IP) = Den2
      DenX(3,IP) = Den3
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENDerivF(NP,     NBas,   thrsh,  nbf,    DA,
     $                        VAO,    VAOX,   INB,    Den,    DenX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density derivatives at the current batch of grid points
C  ** Same as <GetDENDeriv1> but without threshold testing **
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  nbf     -  number of "non-zero" basis functions
C  DA      -  density matrix (lower triangle)
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  basis function derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAOX
C  DEN     -  density values at grid points
C  DenX    -  on exit contains density derivatives
C
C
      DIMENSION DA(NBas*(NBas+1)/2),VAO(*),VAOX(3,*),nbf(NP+1),
     $          INB(NBas*NP),DEN(NP),DenX(3,NP)
C
      PARAMETER (Zero=0.0d0)
C
C
      DO 30 IP=1,NP
c
      Den1 = Zero
      Den2 = Zero
      Den3 = Zero
      If(Abs(DEN(IP)).LT.thrsh) GO TO 29
c
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      DO 19 J=nbf(IP)+1,I
      JJ = INB(J)
      IJ = IT+JJ
      DAIJ = DA(IJ)*Val
c  (closed-shell should be doubled, but taken into account elsewhere)
      Den1 = Den1 + DAIJ*VAOX(1,J)
      Den2 = Den2 + DAIJ*VAOX(2,J)
      Den3 = Den3 + DAIJ*VAOX(3,J)
 19   CONTINUE
      DO 18 J=I+1,nbf(IP+1)
      JJ = INB(J)
      IJ = (JJ*(JJ-1))/2 + II
      DAIJ = DA(IJ)*Val
c  (closed-shell should be doubled, but taken into account elsewhere)
      Den1 = Den1 + DAIJ*VAOX(1,J)
      Den2 = Den2 + DAIJ*VAOX(2,J)
      Den3 = Den3 + DAIJ*VAOX(3,J)
 18   CONTINUE
 20   CONTINUE
c
 29   CONTINUE
      DenX(1,IP) = Den1
      DenX(2,IP) = Den2
      DenX(3,IP) = Den3
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENDerivU(NP,     NBas,   thrsh,  nbf,    DA,
     $                        DB,     VAO,    VAOX,   INB,    DEN,
     $                        DenAX,  DenBX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density derivatives at the current batch of grid points
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  nbf     -  number of "non-zero" basis functions
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta  density matrix (lower triangle)
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  basis function derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAOX
C  DEN     -  density (alpha+beta) values at grid points
C  DenAX   -  on exit contains alpha density derivatives
C  DenBX   -   ditto beta density derivatives
C
C
      DIMENSION DA(*),DB(*),VAO(*),VAOX(3,*),nbf(NP+1),
     $          INB(NBas*NP),DEN(NP),DenAX(3,NP),DenBX(3,NP)
C
      PARAMETER (Zero=0.0d0)
C
C
      DO 30 IP=1,NP
c
      DenA1 = Zero
      DenA2 = Zero
      DenA3 = Zero
      DenB1 = Zero
      DenB2 = Zero
      DenB3 = Zero
      If(Abs(DEN(IP)).LT.thrsh) GO TO 29
c
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      DO 19 J=nbf(IP)+1,I
      JJ = INB(J)
      IJ = IT + JJ
      DAIJ = DA(IJ)*Val
      DBIJ = DB(IJ)*Val
      DenA1 = DenA1 + DAIJ*VAOX(1,J)
      DenA2 = DenA2 + DAIJ*VAOX(2,J)
      DenA3 = DenA3 + DAIJ*VAOX(3,J)
      DenB1 = DenB1 + DBIJ*VAOX(1,J)
      DenB2 = DenB2 + DBIJ*VAOX(2,J)
      DenB3 = DenB3 + DBIJ*VAOX(3,J)
 19   CONTINUE
      DO 18 J=I+1,nbf(IP+1)
      JJ = INB(J)
      IJ = (JJ*(JJ-1))/2 + II
      DAIJ = DA(IJ)*Val
      DBIJ = DB(IJ)*Val
      DenA1 = DenA1 + DAIJ*VAOX(1,J)
      DenA2 = DenA2 + DAIJ*VAOX(2,J)
      DenA3 = DenA3 + DAIJ*VAOX(3,J)
      DenB1 = DenB1 + DBIJ*VAOX(1,J)
      DenB2 = DenB2 + DBIJ*VAOX(2,J)
      DenB3 = DenB3 + DBIJ*VAOX(3,J)
 18   CONTINUE
 20   CONTINUE
c
 29   CONTINUE
      DenAX(1,IP) = DenA1 + DenA1
      DenAX(2,IP) = DenA2 + DenA2
      DenAX(3,IP) = DenA3 + DenA3
      DenBX(1,IP) = DenB1 + DenB1
      DenBX(2,IP) = DenB2 + DenB2
      DenBX(3,IP) = DenB3 + DenB3
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE GetDENDeriv1U(NP,     NBas,   thrsh,  nbf,    DA,
     $                         DB,     DM,     VAO,    VAOX,   INB,
     $                         VMX,    DEN,    DenAX,  DenBX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  forms the density derivatives at the current batch of grid points
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglect of contribution
C  nbf     -  number of "non-zero" basis functions
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta  density matrix (lower triangle)
C  DM      -  array containing maximum density matrix element per column
C  VAO     -  "non-zero" basis function values at grid point
C  VAOX    -  basis function derivatives at grid points
C  INB     -  indexing array for non-zero entries to VAOX
C  VMX     -  array containing maximum magnitude AO derivative per grid point
C  DEN     -  density (alpha+beta) values at grid points
C  DenAX   -  on exit contains alpha density derivatives
C  DenBX   -   ditto beta density derivatives
C
C
      DIMENSION DA(*),DB(*),DM(NBas),VAO(*),VAOX(3,*),nbf(NP+1),
     $          INB(NBas*NP),VMX(NP),DEN(NP),DenAX(3,NP),DenBX(3,NP)
C
      PARAMETER (Zero=0.0d0)
C
C
      DO 30 IP=1,NP
c
      DenA1 = Zero
      DenA2 = Zero
      DenA3 = Zero
      DenB1 = Zero
      DenB2 = Zero
      DenB3 = Zero
      If(Abs(DEN(IP)).LT.thrsh) GO TO 29
c
      VMxx = VMX(IP)
      DO 20 I=nbf(IP)+1,nbf(IP+1)
      II = INB(I)
      IT = (II*(II-1))/2
      Val = VAO(I)
      If(DM(II)*VMxx*Abs(Val).LT.thrsh) GO TO 20
      DO 19 J=nbf(IP)+1,I
      JJ = INB(J)
      IJ = IT+JJ
      DAIJ = DA(IJ)*Val
      DBIJ = DB(IJ)*Val
      DenA1 = DenA1 + DAIJ*VAOX(1,J)
      DenA2 = DenA2 + DAIJ*VAOX(2,J)
      DenA3 = DenA3 + DAIJ*VAOX(3,J)
      DenB1 = DenB1 + DBIJ*VAOX(1,J)
      DenB2 = DenB2 + DBIJ*VAOX(2,J)
      DenB3 = DenB3 + DBIJ*VAOX(3,J)
 19   CONTINUE
      DO 18 J=I+1,nbf(IP+1)
      JJ = INB(J)
      IJ = (JJ*(JJ-1))/2 + II
      DAIJ = DA(IJ)*Val
      DBIJ = DB(IJ)*Val
      DenA1 = DenA1 + DAIJ*VAOX(1,J)
      DenA2 = DenA2 + DAIJ*VAOX(2,J)
      DenA3 = DenA3 + DAIJ*VAOX(3,J)
      DenB1 = DenB1 + DBIJ*VAOX(1,J)
      DenB2 = DenB2 + DBIJ*VAOX(2,J)
      DenB3 = DenB3 + DBIJ*VAOX(3,J)
 18   CONTINUE
 20   CONTINUE
c
 29   CONTINUE
      DenAX(1,IP) = DenA1 + DenA1
      DenAX(2,IP) = DenA2 + DenA2
      DenAX(3,IP) = DenA3 + DenA3
      DenBX(1,IP) = DenB1 + DenB1
      DenBX(2,IP) = DenB2 + DenB2
      DenBX(3,IP) = DenB3 + DenB3
 30   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE POTSort(dft,    lsemi,  NP,     thrsh,  nrec,
     $                   POld,   POldX,  POT,    POTX,   INDX,
     $                   NPP)
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
C  POld    -  potential at previous SCF iteration
C  POldX   -    ditto  potential derivatives
C  POT     -  potential at each grid point
C  POTX    -    ditto  potential derivatives
C
C  on exit
C
C  INDX    -  index into contributing columns of POT
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION POld(NP),POldX(3,NP),POT(NP),POTX(3,NP),INDX(NP)
      INTEGER dft
C
C
      NPP = 0
c
      IF(lsemi.GT.1) THEN
C
C  read old potential for current batch of grid points
C  replace on file by current potential
C
        IF(dft.GT.3) THEN
          READ(41,rec=nrec) POld,POldX
          WRITE(41,rec=nrec) POT,POTX
C
C  form difference potential & derivatives
C
          DO 10 IP=1,NP
          POT(IP) = POT(IP) - POld(IP)
          POTX(1,IP) = POTX(1,IP) - POldX(1,IP)
          POTX(2,IP) = POTX(2,IP) - POldX(2,IP)
          POTX(3,IP) = POTX(3,IP) - POldX(3,IP)
          If(Abs(POT(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POT(NPP) = POT(IP)
            POTX(1,NPP) = POTX(1,IP)
            POTX(2,NPP) = POTX(2,IP)
            POTX(3,NPP) = POTX(3,IP)
          EndIf
 10       CONTINUE
cc
        ELSE
cc
          READ(41,rec=nrec) POld
          WRITE(41,rec=nrec) POT
C
C  form difference potential
C
          DO 15 IP=1,NP
          POT(IP) = POT(IP) - POld(IP)
          If(Abs(POT(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POT(NPP) = POT(IP)
          EndIf
 15       CONTINUE
        ENDIF
cc
      ELSE
C
C  do not use difference potential
C  save current potential
C
        IF(dft.GT.3) THEN
          If(lsemi.GT.0) WRITE(41,rec=nrec) POT,POTX
C
C  sort absolute potential & derivatives
C
          DO 20 IP=1,NP
          If(Abs(POT(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POT(NPP) = POT(IP)
            POTX(1,NPP) = POTX(1,IP)
            POTX(2,NPP) = POTX(2,IP)
            POTX(3,NPP) = POTX(3,IP)
          EndIf
 20       CONTINUE
cc
        ELSE
          If(lsemi.GT.0) WRITE(41,rec=nrec) POT
C
C  sort absolute potential
C
          DO 25 IP=1,NP
          If(Abs(POT(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POT(NPP) = POT(IP)
          EndIf
 25       CONTINUE
        ENDIF
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE POTSortU(dft,    lsemi,  NP,     thrsh,  nrec,
     $                    POldA,  POldB,  POldAX, POldBX, POTA,
     $                    POTB,   POTAX,  POTBX,  INDX,   NPP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts this batch of potentials, determining which grid points
C  will be retained for the numerical integration
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  dft     -  DFT flag
C             of importance here is if we have a non-local functional
C  lsemi   -  integer flag, if >0 use delta potential
C  NP      -  number of grid points this batch
C  thrsh   -  threshold for neglecting contribution
C  nrec    -  record counter for direct access file
C  POldA   -  alpha potential at previous SCF iteration
C  POldB   -  beta  potential at previous SCF iteration
C  POldAX  -  alpha potential derivatives at previous SCF iteration
C  POldBX  -  beta  potential derivatives at previous SCF iteration
C  POTA    -  alpha potential at each grid point
C  POTB    -  beta  potential at each grid point
C  POTAX   -  alpha potential derivatives at each grid point
C  POTBX   -  beta  potential derivatives at each grid point
C
C  on exit
C
C  INDX    -  index into contributing columns of POT/POTX
C  NPP     -  number of "non-zero" grid points
C
C
      DIMENSION POldA(NP),POldB(NP),POldAX(3,NP),POldBX(3,NP),
     $          POTA(NP),POTB(NP),POTAX(3,NP),POTBX(3,NP),INDX(NP)
      INTEGER dft
C
C
      NPP = 0
c
      IF(lsemi.GT.1) THEN
C
C  read old potential for current batch of grid points
C  replace on file by current potential
C
        IF(dft.GT.3) THEN
          READ(41,rec=nrec) POldA,POldB,POldAX,POldBX
          WRITE(41,rec=nrec) POTA,POTB,POTAX,POTBX
C
C  form difference potential & derivatives
C
          DO 10 IP=1,NP
          POTA(IP) = POTA(IP) - POldA(IP)
          POTB(IP) = POTB(IP) - POldB(IP)
          POTAX(1,IP) = POTAX(1,IP) - POldAX(1,IP)
          POTAX(2,IP) = POTAX(2,IP) - POldAX(2,IP)
          POTAX(3,IP) = POTAX(3,IP) - POldAX(3,IP)
          POTBX(1,IP) = POTBX(1,IP) - POldBX(1,IP)
          POTBX(2,IP) = POTBX(2,IP) - POldBX(2,IP)
          POTBX(3,IP) = POTBX(3,IP) - POldBX(3,IP)
          If(Abs(POTA(IP))+Abs(POTB(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POTA(NPP) = POTA(IP)
            POTB(NPP) = POTB(IP)
            POTAX(1,NPP) = POTAX(1,IP)
            POTAX(2,NPP) = POTAX(2,IP)
            POTAX(3,NPP) = POTAX(3,IP)
            POTBX(1,NPP) = POTBX(1,IP)
            POTBX(2,NPP) = POTBX(2,IP)
            POTBX(3,NPP) = POTBX(3,IP)
          EndIf
 10       CONTINUE
cc
        ELSE
cc
          READ(41,rec=nrec) POldA,POldB
          WRITE(41,rec=nrec) POTA,POTB
C
C  form difference potential
C
          DO 15 IP=1,NP
          POTA(IP) = POTA(IP) - POldA(IP)
          POTB(IP) = POTB(IP) - POldB(IP)
          If(Abs(POTA(IP))+Abs(POTB(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POTA(NPP) = POTA(IP)
            POTB(NPP) = POTB(IP)
          EndIf
 15       CONTINUE
        ENDIF
cc
      ELSE
C
C  do not use difference potential
C  save current potential
C
        IF(dft.GT.3) THEN
          If(lsemi.GT.0) WRITE(41,rec=nrec) POTA,POTB,POTAX,POTBX
C
C  sort absolute potential & derivatives
C
          DO 20 IP=1,NP
          If(Abs(POTA(IP))+Abs(POTB(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POTA(NPP) = POTA(IP)
            POTB(NPP) = POTB(IP)
            POTAX(1,NPP) = POTAX(1,IP)
            POTAX(2,NPP) = POTAX(2,IP)
            POTAX(3,NPP) = POTAX(3,IP)
            POTBX(1,NPP) = POTBX(1,IP)
            POTBX(2,NPP) = POTBX(2,IP)
            POTBX(3,NPP) = POTBX(3,IP)
          EndIf
 20       CONTINUE
cc
        ELSE
          If(lsemi.GT.0) WRITE(41,rec=nrec) POTA,POTB
C
C  sort absolute potential
C
          DO 25 IP=1,NP
          If(Abs(POTA(IP))+Abs(POTB(IP)).GT.thrsh) Then
            NPP = NPP+1
            INDX(NPP) = IP
            POTA(NPP) = POTA(IP)
            POTB(NPP) = POTB(IP)
          EndIf
 25       CONTINUE
        ENDIF
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE FullDEN(dft,    lsemi,  NP,     nrec,   DOld,
     $                   DOldX,  DEN,    DENX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generate the full density and density derivatives from difference
C  quantities and values from the previous cycle
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  DFT flag
C             of importance here is if we have a non-local functional
C  lsemi   -  integer flag, if >0 delta densities are being used
C  NP      -  number of grid points this batch
C  nrec    -  record counter for direct access file
C  DOld    -  difference density at each grid point
C  DOldX   -    ditto  difference density derivatives
C  DEN     -  old density  at each grid point
C             on exit contains full density
C  DENX    -    ditto  old density derivatives
C             on exit contains full density derivatives
C
C
      DIMENSION DOld(NP),DOldX(3,NP),DEN(NP),DENX(3,NP)
      INTEGER dft
C
C
      IF(lsemi.GT.1) THEN
C
C  read old density for current batch of grid points
C  replace on file by current density
C
        IF(dft.GT.3) THEN
          READ(43,rec=nrec) DEN,DENX
C
C  form full density and density derivatives
C
          DO 10 IP=1,NP
          DEN(IP) = Abs(DOld(IP) + DEN(IP))    ! must be +ve
          DENX(1,IP) = DOldX(1,IP) + DENX(1,IP)
          DENX(2,IP) = DOldX(2,IP) + DENX(2,IP)
          DENX(3,IP) = DOldX(3,IP) + DENX(3,IP)
 10       CONTINUE
          WRITE(43,rec=nrec) DEN,DENX
cc
        ELSE
cc
          READ(43,rec=nrec) DEN
C
C  form full density
C
          DO 15 IP=1,NP
          DEN(IP) = Abs(DOld(IP) + DEN(IP))    ! must be +ve
 15       CONTINUE
          WRITE(43,rec=nrec) DEN
        ENDIF
cc
      ELSE
C
C  do not use difference density
C  save current density
C
        DO 20 IP=1,NP
        DEN(IP) = Abs(DOld(IP))   ! must be +ve
 20     CONTINUE
cc
        IF(dft.GT.3) THEN
          CALL CpyVEC(3*NP,DOldX,DENX)
          If(lsemi.GT.0) WRITE(43,rec=nrec) DEN,DENX
        ELSE
          If(lsemi.GT.0) WRITE(43,rec=nrec) DEN
        ENDIF
cc
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE FullDENU(dft,    lsemi,  NP,     nrec,   DOldA,
     $                    DOldB,  DOldAX, DOldBX, DENA,   DENB,
     $                    DENAX,  DENBX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generate the full density and density derivatives from difference
C  quantities and values from the previous cycle
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  dft     -  DFT flag
C             of importance here is if we have a non-local functional
C  lsemi   -  integer flag, if >0 delta densities are being used
C  NP      -  number of grid points this batch
C  nrec    -  record counter for direct access file
C  DOldA   -  alpha difference density at each grid point
C  DOldB   -  beta  difference density at each grid point
C  DOldAX  -  alpha difference density derivatives at each grid point
C  DOldBX  -  beta  difference density derivatives at each grid point
C  DENA    -  old alpha density  at each grid point
C             on exit contains full alpha density
C  DENB    -  old beta  density  at each grid point
C             on exit contains full beta  density
C  DENAX   -  old alpha density derivatives at each grid point
C             on exit contains full alpha density derivatives
C  DENBX   -  old beta  density derivatives at each grid point
C             on exit contains full beta  density derivatives
C
C
      DIMENSION DOldA(NP),DOldB(NP),DOldAX(3,NP),DOldBX(3,NP),
     $          DENA(NP),DENB(NP),DENAX(3,NP),DENBX(3,NP)
      INTEGER dft
C
C
      IF(lsemi.GT.1) THEN
C
C  read old density for current batch of grid points
C  replace on file by current density
C
        IF(dft.GT.3) THEN
          READ(43,rec=nrec) DENA,DENB,DENAX,DENBX
C
C  form full density and density derivatives
C
          DO 10 IP=1,NP
          DENA(IP) = Abs(DOldA(IP) + DENA(IP))   ! must be +ve
          DENB(IP) = Abs(DOldB(IP) + DENB(IP))   ! must be +ve
          DENAX(1,IP) = DOldAX(1,IP) + DENAX(1,IP)
          DENAX(2,IP) = DOldAX(2,IP) + DENAX(2,IP)
          DENAX(3,IP) = DOldAX(3,IP) + DENAX(3,IP)
          DENBX(1,IP) = DOldBX(1,IP) + DENBX(1,IP)
          DENBX(2,IP) = DOldBX(2,IP) + DENBX(2,IP)
          DENBX(3,IP) = DOldBX(3,IP) + DENBX(3,IP)
 10       CONTINUE
          WRITE(43,rec=nrec) DENA,DENB,DENAX,DENBX
cc
        ELSE
cc
          READ(43,rec=nrec) DENA,DENB
C
C  form full density
C
          DO 15 IP=1,NP
          DENA(IP) = Abs(DOldA(IP) + DENA(IP))   ! must be +ve
          DENB(IP) = Abs(DOldB(IP) + DENB(IP))   ! must be +ve
 15       CONTINUE
          WRITE(43,rec=nrec) DENA,DENB
        ENDIF
cc
      ELSE
C
C  do not use difference density
C  save current density
C
        DO 20 IP=1,NP
        DENA(IP) = Abs(DOldA(IP))   ! must be +ve
        DENB(IP) = Abs(DOldB(IP))   ! must be +ve
 20     CONTINUE
cc
        IF(dft.GT.3) THEN
          CALL CpyVEC(3*NP,DOldAX,DENAX)
          CALL CpyVEC(3*NP,DOldBX,DENBX)
          If(lsemi.GT.0) WRITE(43,rec=nrec) DENA,DENB,DENAX,DENBX
        ELSE
          If(lsemi.GT.0) WRITE(43,rec=nrec) DENA,DENB
        ENDIF
cc
      ENDIF
C
      RETURN
      END
