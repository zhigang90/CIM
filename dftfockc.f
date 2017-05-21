      SUBROUTINE DFTFockC(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                    NSym,   NGen,   ISYM,   NEqATM, IPRNT,
     $                    NRad,   NAng,   IradQ,  factor, NBatch,
     $                    DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                    XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                    BASDAT, INX,    BL,     ExpMIN, DA,
     $                    DM,     NScr,   Z,      XC,     EXC,
     $                    EL,     lgrid,  lsemi,  IEntry)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates DFT contribution to Exchange-Correlation Fock Matrix
C  and total Exchange-Correlation energy
C  ** CLOSED SHELL **
C
C    Jon Baker   March/April 1997    (revised March 2002)
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All standard methods include Slater exchange)
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
C  DA      -  closed-shell density matrix (lower triangle)
C  DM      -  maximum density matrix element per column
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  XC      -  DFT Exchange-Correlation matrix (accumulated)
C  EXC     -  Exchange-Correlation energy (accumulated)
C  EL      -  number of electrons over the grid (accumulated)
C
C  lgrid   -  integer flag for generation/regeneration of DFT grid
C              0 - do not generate grid, already exists
C              1 - first entry so generate grid
C              2 - regenerate grid and REMOVE old grid files
C  lsemi   -  integer flag for switching on difference potential
C  IEntry  -  indicates initial/final subroutine entry this SCF cycle
C
C  ----------------------------------------------------------------
C  references
C
C     "A multicenter numerical integration scheme for polyatomic
C      molecules"
C      A.D.Becke      J.Chem.Phys.  88 (1988) 2547
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),ISYM(NGen),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     $          DA(NBas*(NBas+1)/2),DM(NBas),XC(NBas*(NBas+1)/2)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
      CHARACTER*256 Filename,Filenam1,scrfile
      CHARACTER*3 ch3
      save tgrid,tao,tden,tdenx,tpot,tnum,mpoint,npptot
C
      PARAMETER (IdWt=0,Zero=0.0d0,Two=2.0d0)
C
C
      IF(IEntry.EQ.1) THEN
        EXC = Zero
        EL  = Zero
cc
        tgrid = Zero
        tao = Zero
        tden = Zero
        tdenx = Zero
        tpot = Zero
        tnum = Zero
        mpoint = 0
        npptot = 0
cc
        IEntry = 0
      ENDIF
c
      NBasS = NBas*NBatch
      NBasT = 3*NBasS
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
      if(IPrnt.gt.2) then
        write(6,*) '** DFT thresholds  thrsh:',thrsh,' thint:',thint
      end if
C
C  allocate scratch pointers
C
      iao = 1                         ! AO values over grid points
      inb = iao + NBasS               ! index array for "non-zero" AOs
      ibf = inb + NBasS               ! number of "non-zero" AOs per point
      id  = ibf + NBatch+1            ! density per grid point
      iv  = id  + NBatch              ! potential per grid point
      ind = iv  + NBatch              ! index of "non-zero" densities
      iscr = ind + NBatch             ! scratch storage (NBatch)
      IEnd = iscr + NBatch - 1
cc
      i1 = 1                          ! scratch pointers for grid
      i2 = i1 + NAtoms
      i3 = i2 + NAtoms
cc
      If(dft.GT.3) Then
C
C  additional scratch space for non-local functionals
C
        iaox = 1 + IEnd               ! AO derivatives
        idx  = iaox + NBasT           ! density derivatives
        ivx  = idx  + 3*NBatch        ! potential derivatives
        IEnd = ivx  + 3*NBatch - 1
      EndIf
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'DFTFockC')
C
C
C  Begin calculation for this centre
C  ---------------------------------
C
C  assign grid filename for this centre
C
      call secund(t1)
      Call getchval('scrf',scrfile)
      Call rmblan(scrfile,256,len)
      Write(ch3,'(I3)') ICntr
      If(ICntr.LT.100) ch3(1:1) = '0'
      If(ICntr.LT.10)  ch3(2:2) = '0'
      Filename = scrfile(1:len)//'.grid.'//ch3
      len1 = len+9
C
C  open grid potential direct access file and grid AO file
C
      If(lsemi.GT.0) Then
        Filenam1 = scrfile(1:len)//'.pots.'//ch3
        lrec = 32*NBatch
        If(dft.LE.3) lrec = 8*NBatch
        OPEN (UNIT=41,FILE=Filenam1(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec)
cc        Filenam1 = scrfile(1:len)//'.aovs.'//ch3
cc        OPEN (UNIT=42,FILE=Filenam1(1:len1),FORM='UNFORMATTED',
cc     $        STATUS='UNKNOWN')
        Filenam1 = scrfile(1:len)//'.dens.'//ch3
        OPEN (UNIT=43,FILE=Filenam1(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec)
      EndIf
c
      IF(lgrid.GT.0) THEN
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
C
C -- save grid to disk
        CALL WrGRID(Filename,len1,NPoint,XGRID,WGHT)
c
      ELSE
C
C -- read existing file back from disk
        CALL RdGRID(Filename,len1,NPoint,XGRID,WGHT)
c
      ENDIF
      call secund(t2)
      t1 = t2-t1
cc      write(6,*) ' Back from <GRID> for Atom: ',ICntr
cc      write(6,1111) NPoint,t1
      tgrid = tgrid + t1
 1111 format(' Number of grid points is ',I6,' time: ',f7.2)
C
C  ..........................................................
C  on return from <GRID>  NPoint is number of grid points
C  ..........................................................
C
C  loop over grid points on atom IAtom
C
      nrec = 0               ! record counter for DA file
      IP = 1
      NLeft = NPoint
      mpoint = mpoint+NPoint
 50   CONTINUE
      nrec = nrec+1
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
      call zeroit(Z(iao),NBasS)
cc      IF(lsemi.LE.1) THEN
C
C  calculate AO and derivatives values
C
        If(dft.LE.3) Then
          CALL AOVal(NP,     NShell, NBas,   XGRID(1,IP),
     $               XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $               ExpMIN, Z(iao))
          CALL AOSort(NP,     NBas,   thint,  Z(iao), Z(ibf),
     $                Z(inb), Z(iao), Z(iscr))
        Else
          CALL zeroit(Z(iaox),NBasT)
          CALL AOGrad(NP,     NShell, NBas,   XGRID(1,IP),
     $                XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $                ExpMIN, Z(iao), Z(iaox))
          CALL AOXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $                 Z(ibf), Z(inb), Z(iao), Z(iaox),Z(iscr))
        EndIf
c
cc        If(lsemi.EQ.1) Then
C
C  save AO and derivative values to disk
C
cc          If(dft.LE.3) Then
cc            CALL WrAOVal(NP+1, Z(ibf), Z(iao), Z(inb), Z(iscr))
cc          Else
cc            CALL WrAOGrad(NP+1, Z(ibf), Z(iao), Z(iaox), Z(inb),
cc     $                    Z(iscr))
cc          EndIf
cc        EndIf
cc      ELSE
C
C  Values already exist on disk - simply read back
C
cc        If(dft.LE.3) Then
cc          CALL RdAOVal(NP+1, Z(ibf), Z(iao), Z(inb), Z(iscr))
cc        Else
cc          CALL RdAOGrad(NP+1, Z(ibf), Z(iao), Z(iaox), Z(inb),
cc     $                  Z(iscr))
cc        EndIf
cc      ENDIF
      call secund(t2)
      tao = tao + t2-t1
C
C  form density at current grid point
C
      call secund(t1)
      CALL GetDENS1(NP,     NBas,   thrsh,  Z(ibf), DA,
     $              DM,     Z(iao), Z(inb), Z(iscr),Z(iv))
      call secund(t2)
      tden = tden + t2-t1
C
C  ..........................................................
C  on return from <GetDENC>  Z(iv) holds current densities
C  ..........................................................
C
C  prepare for numerical integration
C
C  form density derivatives if non-local
C
      If(dft.GT.3) Then
        call secund(t1)
        CALL GetDENDeriv1(NP,     NBas,   thrsh,  Z(ibf), DA,
     $                    DM,     Z(iao), Z(iaox),Z(inb), Z(iscr),
     $                    Z(iv),  Z(ivx))
        call secund(t2)
        tdenx = tdenx + t2-t1
      EndIf
C
C  form full density and density derivatives from difference quantities
C
      CALL FullDEN(dft,    lsemi,  NP,     nrec,   Z(iv),
     $             Z(ivx), Z(id),  Z(idx))
C
C  calculate potential (and potential derivatives)
C
      call secund(t1)
      CALL VXCPOT(dft,    NP,     thrsh,  IdWt,   Z(id),
     $           Z(idx), WGHT(IP),Z(iv),  Z(ivx), EXC,
     $            EL,    Z(iscr))
C
C  determine difference potential and number of contributing points
C
      CALL POTSort(dft,    lsemi,  NP,     thint,  nrec,
     $             Z(id),  Z(idx), Z(iv),  Z(ivx), Z(ind),
     $             NPP)
      npptot = npptot+npp
      call secund(t2)
      tpot = tpot + t2-t1
C
C  Accumulate DFT contribution into XC matrix
C
      call secund(t1)
      CALL MakeXC(dft,    NPP,    NBas,   Z(ibf), thrsh,
     $           WGHT(IP),Z(iv),  Z(ivx), Z(iao), Z(iaox),
     $           Z(inb),  Z(ind), XC)
      call secund(t2)
      tnum = tnum + t2-t1
cc
      IP = IP+NP
      NLeft = NLeft-NP
      If(NLeft.GT.0) GO TO 50
C
C  end loop over grid points
C
      If(IEntry.EQ.-1.AND.IPrnt.GT.1)
     $   write(6,2222) tgrid,tao,tden,tdenx,tpot,tnum,mpoint,npptot
C
C  close grid potential direct access file and grid AO file
C
      If(lsemi.GT.0) Then
        close( unit=41, status='keep' )
        close( unit=43, status='keep' )
      endif
C
      RETURN
c
 2222   Format(/,' Time taken to calculate grid:      ',f8.2,/,
     $           ' Time to calculate & sort AOs:      ',f8.2,/,
     $           ' Time to calculate density:         ',f8.2,/,
     $           ' Time to calculate den.derivs:      ',f8.2,/,
     $           ' Time to evaluate & sort potential: ',f8.2,/,
     $           ' Time numerical integration:        ',f8.2,/,
     $           ' Total number of grid points:       ',I8,/,
     $           ' Number of contributing points:     ',I8)
c
      END
