      SUBROUTINE DFTGradU(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                    NSym,   NGen,   ISYM,   NEqATM, IPRNT,
     $                    NRad,   NAng,   IradQ,  factor, NBatch,
     $                    DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                    XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                    BASDAT, INX,    BL,     ExpMIN, IdWt,
     $                    NBAtm,  DA,     DB,     NScr,   Z,
     $                    GXA,    GYA,    GZA,    GXB,    GYB,
     $                    GZB,    GWT,    GC,     EL,     IEntry)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates DFT contribution to Exchange-Correlation gradient
C  ** UNRESTRICTED OPEN SHELL **
C
C    Jon Baker   April 1998
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All methods include Slater exchange)
C              0 - no dft contribution (this routine should NOT be called)
C              1 - Slater exchange only
C              2 - VWN (Slater exchange + local RPA correlation)
C              3 - VWN5 (ditto but RECOMMENDED Ceperley-Alder fit)
C              4 - Becke 88 exchange
C              5 - Becke 88 + VWN (local correlation)
C              6 - Becke 88 + VWN5 (local correlation)
C              7 - Becke 88 + Perdew 86
C              8 - Becke 88 + Perdew-Wang 91
C              9 - Becke 88 + LYP (correlation)
C             10 - Becke 88 + VWN5 (local) + Perdew 86 (nonlocal)
C             11 - Handy/Cohen exchange
C             12 - Handy/Cohen + VWN
C             13 - Handy/Cohen + VWN5
C             14 - Handy/Cohen + Perdew 86
C             15 - Handy/Cohen + Perdew-Wang 91
C             16 - Handy/Cohen + LYP
C             17 - PW91 exchange + correlation
C             18 - PBE exchange + correlation
C             19 - O3LYP   (VWN5 + OPTX     + LYP  + HF exchange)
C             20 - B3LYP   (VWN  + Becke 88 + LYP  + HF exchange)
C             21 - B3PW91  (VWN5 + Becke 88 + PW91 + HF exchange)
C             22 - WAH     (VWN5 + Becke 88 + LYP  + HF exchange;
C                           modified B3LYP specifically for NMR)
C             23 - user-defined functional
C -- The following functionals are self-contained and were downloaded from the web
C             24 - B97    (Becke 97 exchange-correlation functional)
C             25 - B97-1  ( ditto, but reparametrized)
C             26 - B97-2  ( ditto, but reparametrized again)
C             27 - HCTH   (Hamprecht-Cohen-Tozer-Handy)
C                  (this is similar in form to B97, but has no HF exchange)
C  ICntr   -  current atom for DFT contribution being calculated
C  NAtoms  -  number of atoms
C  XNuc    -  nuclear coordinates
C  IAN     -  nuclear charge (atomic number)
C  NSym    -  number of symmetry operations
C             (currently highest abelian subgroup - excludes identity)
C  NGen    -  number of group generators
C  ISYM    -  list of symmetry operations (generators)
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of NQ symmetry-unique atoms
C  IPrnt   -  print flag
C  NRad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  NAng    -  angular quadrature order for numerical grid
C             (if zero, will be determined in grid generation routines)
C  IradQ   -  radial quadrature type
C              0 - Euler-Maclaurin (Handy)          default
C              1 - Chebyshev second type (Becke mapping)
C  factor  -  determines  grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  NBatch  -  maximum number of grid points treated in one batch
C  DISTN   -  distance to nearest neighbour for each atom
C  AIJ     -  Becke grid parameters (see A2 in Becke article)
C  RDIST   -  inverse interatomic distance array
C  XXA     -  storage for angular quadrature grid
C  WTA     -  storage for angular quadrature weights
C  XGRID   -  X,Y,Z grid points
C             (GRID array contains enough space for maximum number
C              of grid points on "largest" atom)
C  WGHT    -  grid quadrature weights
C  thrsh   -  threshold for neglecting contribution
C  NBas    -  total number of basis functions
C  NShell  -  total number of shells
C  BASDAT  -  basis set data for TEXAS
C  INX     -  more basis set data for TEXAS
C  BL      -  precomputed normalization factors per primitive shell
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C  NBAtm   -  number of basis functions per atom
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta density matrix (lower triangle)
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  GXA     -  alpha x DFT derivative matrix
C  GYA     -  alpha y DFT derivative matrix
C  GZA     -  alpha z DFT derivative matrix
C  GXB     -  beta x DFT derivative matrix
C  GYB     -  beta y DFT derivative matrix
C  GZB     -  beta z DFT derivative matrix
C  GWT     -  grid weight derivatives (if calculated)
C  GC      -  DFT contribution to Cartesian gradient
C  EL      -  number of electrons over the grid (accumulated)
C
C  IEntry  -  indicates initial/final subroutine entry
C
C
C  ----------------------------------------------------------------
C  references
C
C     "A multicenter numerical integration scheme for polyatomic
C      molecules"
C      A.D.Becke      J.Chem.Phys.  88 (1988) 2547
C
C     "The performance of a family of density functionals"
C      B.G.Johnson, P.M.W.Gill and J.A.Pople
C      J.Chem.Phys.  98 (1993) 5612
C
C     "The effect of grid quality and weight derivatives in density
C      functional calculations"
C      J.Baker, J.Andzelm, A.Scheiner and B.Delley
C      J.Chem.Phys.  101 (1994) 8894
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),ISYM(NGen),NBAtm(NAtoms),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),
     $          GXA(NBas,NBas),GYA(NBas,NBas),GZA(NBas,NBas),
     $          GXB(NBas,NBas),GYB(NBas,NBas),GZB(NBas,NBas),
     $          GWT(3,NAtoms),GC(3,NAtoms)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
C
      PARAMETER (Zero=0.0d0)
C
C
      IF(IEntry.EQ.1) THEN
        DO 10 J=1,NBas
        DO 10 I=1,NBas
        GXA(I,J) = Zero
        GYA(I,J) = Zero
        GZA(I,J) = Zero
        GXB(I,J) = Zero
        GYB(I,J) = Zero
        GZB(I,J) = Zero
 10     CONTINUE
        CALL ZeroIT(GWT,3*NAtoms)
        CALL ZeroIT(GC,3*NAtoms)
cc
        tgrid = Zero
        tao = Zero
        tden = Zero
        tdenx = Zero
        tpot = Zero
        tnum = Zero
        twght = Zero
cc
        EL  = Zero
        IEntry = 0
      ENDIF
c
      NBasS = NBas*NBatch
      NBasT = 3*NBasS
      NBasV = 6*NBasS
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
      if(IPrnt.gt.2) then
        write(6,*) '** DFT thresholds  thrsh:',thrsh,' thint:',thint
      end if
C
C  allocate scratch pointers
C
      iao = 1                         ! AO values over grid points
      iaox = iao + NBasS              ! AO derivatives over grid points
      inb = iaox + NBasT              ! index array for "non-zero" AOs
      ibf = inb + NBasS               ! number of "non-zero" AOs per point
      idA = ibf + NBatch+1            ! alpha density per grid point
      idB = idA + NBatch              ! beta density per grid point
      ivA = idB + NBatch              ! alpha potential per grid point
      ivB = ivA + NBatch              ! beta potential per grid point
      ind = ivB + NBatch              ! index of "non-zero" densities
      iscr = ind + NBatch             ! scratch storage (NBatch)
      IEnd = iscr + NBatch - 1
cc
      i1 = inb                        ! scratch pointers for grid
      i2 = i1 + NAtoms
      i3 = i2 + NAtoms
      i4 = i3 + NAtoms*NAtoms
cc
      If(dft.GT.3) Then
C
C  additional scratch space for non-local functionals
C
        iaoxx = 1 + IEnd              ! AO 2nd derivatives over grid points
        idxA = iaoxx+ NBasV           ! alpha density derivatives
        idxB = idxA + 3*NBatch        ! beta density derivatives
        ivxA = idxB + 3*NBatch        ! alpha potential derivatives
        ivxB = ivxA + 3*NBatch        ! beta potential derivatives
        IEnd = ivxB + 3*NBatch - 1
      EndIf
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'DFTFockU')
C
C
C  Begin calculation for this centre
C  ---------------------------------
C
        call secund(t1)
        CALL GRID(NAtoms, XNuc,   IAN,    IPRNT,  NRad,
     $            NAng,   factor, IradQ,  Z(i1),  Z(i2),
     $            Z(i3),  DISTN,  RDIST,  AIJ,    ICntr,
     $            XXA,    WTA,    NPoint, XGRID,  WGHT)
C
C  ==============================================================
C  ** SYMMETRY SORT OF GRID **
        If(NSym.GT.0) Then
         XX = XNuc(1,ICntr)
         YY = XNuc(2,ICntr)
         ZZ = XNuc(3,ICntr)
         CALL SymGRID(NAtoms, ICntr,  XX,     YY,     ZZ,
     $                NSym,   NGen,   ISYM,   NEqATM, NPoint,
     $                XGRID,  WGHT)
        EndIf
C  ===============================================================
        call secund(t2)
        t1 = t2-t1
cc        write(6,*) ' Back from <GRID> for Atom: ',ICntr
cc        write(6,1111) NPoint,t1
        tgrid = tgrid + t1
 1111   format(' Number of grid points is ',I6,' time: ',f7.2)
C
C  ..........................................................
C  on return from <GRID>  NPoint is number of grid points
C  ..........................................................
C
C  loop over grid points on atom IAtom
C
        IP = 1
        NLeft = NPoint
 50     CONTINUE
        NP = NBatch
        If(NP.GT.NLeft) NP=NLeft
C
C  ...........................................................
C    NP is the number of points in one batch on this pass
C  ...........................................................
C
C  form AO and derivative values over grid points
C  and sort "non-zero" AOs
C
        call secund(t1)
        CALL zeroit(Z(iao),NBasS)
        CALL zeroit(Z(iaox),NBasT)
        IF(dft.LE.3) THEN
          CALL AOGrad(NP,     NShell, NBas,   XGRID(1,IP),
     $                XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $                ExpMIN, Z(iao), Z(iaox))
          CALL AOXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $                 Z(ibf), Z(inb), Z(iao), Z(iaox),Z(iscr))
        ELSE
          CALL zeroit(Z(iaoxx),NBasV)
          CALL AODer2(NP,     NShell, NBas,   XGRID(1,IP),
     $                XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $                ExpMIN, Z(iao), Z(iaox),Z(iaoxx))
          CALL AOXXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $                 Z(iaoxx),Z(ibf), Z(inb), Z(iao), Z(iaox),
     $                 Z(iaoxx),z(iscr))
        ENDIF
        call secund(t2)
        tao = tao + t2-t1
C
C  form alpha/beta density at current grid point
C
        call secund(t1)
        CALL GetDENSU(NP,     NBas,   Z(ibf), DA,   DB,
     $                Z(iao), Z(inb), Z(idA), Z(idB))
        call secund(t2)
        tden = tden + t2-t1
C
C  ..........................................................
C  on return from <GetDENC>  Z(id) holds current densities
C  ..........................................................
C
C  prepare for numerical integration
C  sort "non-zero" densities
C
cc        CALL DENSortU(NP,thrsh,Z(idA),Z(idB),Z(ind),NPP)
C
C  form density derivatives
C
        If(dft.GT.3) Then
          call secund(t1)
c
c -- add alpha & beta densities (zero threshold is on density sum)
          CALL AddVEC(NP,Z(idA),Z(idB),Z(ivA))
c
          CALL GetDENDerivU(NP,     NBas,   thrsh,  Z(ibf), DA,
     $                      DB,     Z(iao), Z(iaox),Z(inb), Z(ivA),
     $                      Z(idxA),Z(idxB))
          call secund(t2)
          tdenx = tdenx + t2-t1
        EndIf
C
C  calculate potential (and potential derivatives)
C
        call secund(t1)
        CALL VXCPOTU(dft,    NP,     thrsh,  Z(idA), Z(idB),
     $               Z(idxA),Z(idxB),WGHT(IP),Z(ivA), Z(ivB),
     $               Z(ivxA),Z(ivxB),EXC,    EL,     Z(iscr))
c ...................................................................
c  on exit from VXCPOTU, Z(imo) contains the exchange-correlation
c  energy contribution AT EACH OF NPP GRID POINTS
c ...................................................................
C
C  determine number of contributing points
C
        CALL POTSortU(dft,    0,      NP,     thint,  nrec,
     $                Z(idA), Z(idB), Z(idxA),Z(idxB),Z(ivA),
     $                Z(ivB), Z(ivxA),Z(ivxB),Z(ind), NPP)
        call secund(t2)
        tpot = tpot + t2-t1
C
C  Accumulate DFT contribution into gradient matrices
C
        call secund(t1)
        CALL MakeGXU(dft,    NPP,    NBas,   Z(ibf), thrsh,
     $              WGHT(IP),Z(ivA), Z(ivB), Z(ivxA),Z(ivxB),
     $              Z(iao), Z(iaox),Z(iaoxx),Z(inb), Z(ind),
     $               GXA,    GYA,    GZA,    GXB,    GYB,   GZB)
        call secund(t2)
        tnum = tnum + t2-t1
C
C  ..................................................................
C    Derivatives of Weight Section
C  ..................................................................
C
       If(IdWt.GT.0) Then
        call secund(t1)
        CALL GradWT(NPP, XGRID(1,IP), WGHT(IP), Z(ind), ICntr,
     $              Z(iscr), NAtoms,  XNuc,    AIJ,     RDIST,
     $              Z(i1),   Z(i2),   Z(i3),   Z(i4),   GWT)
        write(6,*) ' back from <GradWT>  GWT is:'
        call prntmat(natoms,3,natoms,gwt)
        call secund(t2)
        twght = twght + t2-t1
       EndIf
C  ..................................................................
C
        IP = IP+NP
        NLeft = NLeft-NP
        If(NLeft.GT.0) GO TO 50
C
C  end loop over grid points
C
C  ................................................................
C  Form full Cartesian gradient from DFT gradient matrices
C
      If(IdWt.GT.0) Then
       CALL PartGRADU(NAtoms, NBas,   NBAtm,  ICntr,  DA,
     $                DB,     GXA,    GYA,    GZA,    GXB,
     $                GYB,    GZB,    GWT,    GC)
      EndIf
C  ................................................................
C
      If(IEntry.EQ.-1.AND.IPrnt.GT.1)
     $   write(6,2222) tgrid,tao,tden,tdenx,tpot,tnum,twght
C
      RETURN
c
 2222   Format(/,' Time taken to calculate grid: ',f7.2,/,
     $           ' Time to calculate & sort AOs: ',f7.2,/,
     $           ' Time to calculate density:    ',f7.2,/,
     $           ' Time to calculate den.derivs: ',f7.2,/,
     $           ' Time to evaluate potential:   ',f7.2,/,
     $           ' Time numerical integration:   ',f7.2,/,
     $           ' Time for weight derivatives:  ',f7.2,/)
c
      END
