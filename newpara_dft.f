c ================================================================
c  MAIN DRIVING ROUTINES FOR PARALLEL DFT
c ================================================================
c
      SUBROUTINE Para_DFT(dft,    NAtoms, XNuc,   IAN,    NSym,
     $                    NGen,   ISYM,   NEqATM, NQ,     IUNQ,
     $                    IPRNT,  NRad,   NAng,   IradQ,  factor,
     $                    NBatch, lgrid,  lsemi,  DISTN,  AIJ,
     $                    RDIST,  XXA,    WTA,    XGRID,  WGHT,
     $                    thrsh,  NBas,   NShell, BASDAT, INX,
     $                    BL,     ExpMIN, rhf,    DA,     DB,
     $                    DM,     NSLAVE, LSlave, LFINI,  NScr,
     $                    Z,      XCA,    XCB,    EXC,    EL)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Drives the slaves for parallel computation of DFT contribution to
C  Exchange-Correlation Fock Matrix and total Exchange-Correlation energy
C
C    Jon Baker   Jan 1998
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
C  lgrid   -  integer flag for generation/regeneration of DFT grid
C              0 - do not generate grid, already exists
C              1 - first entry so generate grid
C              2 - regenerate grid and REMOVE old grid files
C  lsemi   -  integer flag for switching on difference potential
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
C  rhf     -  logical flag indicating if system is closed-shell
C  DA      -  alpha/closed-shell density matrix (lower triangle)
C  DB      -  beta density matrix (lower triangle)
C  DM      -  maximum density matrix element per column
C  NSLAVE  -  list of which atoms go to which slaves
C  LSlave  -  logical array (if that atom has been done)
C  LFINI   -  small array to handle slaves that finish prematurely
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  XCA     -  alpha/closed-shell DFT Exchange-Correlation matrix
C  XCB     -  beta DFT Exchange-Correlation matrix
C  EXC     -  Exchange-Correlation energy
C  EL      -  number of electrons over the grid
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),IUNQ(NQ),ISYM(NGen),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NBas),
     $          XCA(NBas*(NBas+1)/2),XCB(NBas*(NBas+1)/2)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      Dimension NSLAVE(NQ),LFINI(NQ)
      DIMENSION Z(NScr)
      INTEGER dft
      LOGICAL rhf,LSlave(NQ),FirstEntry
      SAVE FirstEntry
c
      DATA FirstEntry/.True./
c
      If(nslv.GT.NQ.AND.FirstEntry) Then
c -- inefficiency warning
        call message('**WARNING** in DFT',
     $  'Not possible to use slaves efficiently: Nslaves > Natoms ',
     *  nslv,NQ)
        FirstEntry = .False.
      EndIf
c -----------------------------------------------------------------
c
c -- check if the slaves are OK
      call para_check
      nsize = (NBas*(NBas+1))/2
c
c -- send variable data to slaves
      call para_initsend
      call para_pack_int(lgrid,1)
      call para_pack_int(lsemi,1)
      call para_pack_int(NQ,1)
      call para_pack_int(IUNQ,NQ)
      call para_pack_int(NSLAVE,NQ)
      call para_pack_real(factor,1)
      call para_pack_real(thrsh,1)
      call para_pack_real(DA,nsize)
      If(.NOT.rhf) call para_pack_real(DB,nsize)
      call para_pack_real(DM,NBas)
      call para_bcast_pack(TxDftDat)
c
      EXC = 0.0d0
      EL = 0.0d0
c
      nfini = 0
      Do IAtom=1,NQ
      LSlave(IAtom) = .False.
      EndDo
c
c -- loop over atomic centres, sending each atom to a different slave
      call para_check
c
      IF(lgrid.GT.0) THEN
c
c -- allow free-for-all and save atom-slave list
        DO 200 IAtom=1,NQ
        ICntr = IUNQ(IAtom)
        call para_recv(imsg,islave,TxDftReq)
        call para_send(ICntr,islave,TxDftAssign)
        NSlave(IAtom) = islave
 200    CONTINUE
cc        write(6,*) ' MASTER:  Finished free-for-all  atom-slave list:'
cc        do i=1,nq
cc        write(6,*) i,'  ',nslave(i)
cc        enddo
cc        call f_lush(6)
      ELSE
c
c -- use existing atom-slave list to determine which slaves get which atoms
c
c -- the time imbalance was not summed properly              
        IAtom = 0
        call elapsec(timo)
        call getrval('imbal',tlost)
 201    CONTINUE
        IAtom = IAtom+1
        call para_recv(imsg,islave,TxDftReq)
        CALL WhichSlave(islave,NQ,NSLAVE,LSlave,JAtom)
        If(JAtom.EQ.-1) Then
          nfini = nfini+1
          call para_send(0,islave,TxDftAssign)
          call elapsec(timb)
          tlost=tlost+(nfini-1)*(timb-timo)
          timo=timb          ! this slave already finished, idling
          IAtom = IAtom-1    ! decrement loop index
        Else
          ICntr = IUNQ(JAtom)
          call para_send(ICntr,islave,TxDftAssign)
        EndIf
        If(IAtom.LT.NQ) GO TO 201
        call setrval('imbal',tlost)
      ENDIF
c
c -- all done  stop all slaves
      call elapsec(timo)
      call getrval('imbal',tlost)
cc      write(6,*) ' MASTER:  About to stop all slaves.  Nfini:',nfini
cc      call f_lush(6)
      Do IAtom=nfini+1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_send(0,islave,TxDftAssign)
      call elapsec(timb)
      tlost=tlost+(IAtom-1)*(timb-timo)
      timo=timb
      EndDo
cc      write(6,*) ' MASTER:  Stopped all slaves'
cc      call f_lush(6)
      call setrval('imbal',tlost)
c
c -- now sum up exchange-correlation matrix and energy from each slave
      call para_reduce(XCA,nsize,TxDftF1)
      call para_reduce(EXC,1,TxDftF2)
      call para_reduce(EL,1,TxDftF3)
      If(.NOT.rhf) call para_reduce(XCB,nsize,TxDftF4)
c
c -- That's it!
c
      RETURN
      END
c =======================================================================
c
      SUBROUTINE Para_DFTG(dft,    NAtoms, XNuc,   IAN,    nsym,
     $                     ngen,   ISYM,   NEqATM, NEqBAS, NQ,
     $                     IUNQ,   IPRNT,  lrad,   lang,   IradQ,
     $                     factor, NBatch, DISTN,  AIJ,    RDIST,
     $                     XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                     NBas,   NShell, BASDAT, INX,    BL,
     $                     ExpMIN, NBAtm,  IdWt,   rhf,    DA,
     $                     DB,     NScr,   Z,      GX,     GY,
     $                     GZ,     GXB,    GYB,    GZB,    GWT,
     $                     GC,     EL)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Drives the slaves for parallel computation of DFT Gradient
C
C    Jon Baker         January 1998
C    (modified for UHF   April 1998)
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
C  NAtoms  -  number of atoms
C  XNuc    -  nuclear coordinates
C  IAN     -  nuclear charge (atomic number)
C  NSym    -  number of symmetry operations
C             (currently highest abelian subgroup - excludes identity)
C  NGen    -  number of group generators
C  ISYM    -  list of symmetry operations (generators)
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NEqBas  -  list of basis function equivalences under symmetry operations
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
C  NBAtm   -  number of basis functions per atom
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C  rhf     -  logical flag for RHF
C  DA      -  density matrix (closed-shell/alpha)
C  DB      -  density matrix (beta)
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  GX      -  x DFT derivative matrix (closed-shell/alpha)
C  GY      -  y DFT derivative matrix (closed-shell/alpha)
C  GZ      -  z DFT derivative matrix (closed-shell/alpha)
C  GXB     -  x DFT derivative matrix (beta)
C  GYB     -  y DFT derivative matrix (beta)
C  GZB     -  z DFT derivative matrix (beta)
C  GWT     -  grid weight derivatives (if calculated)
C  GC      -  DFT contribution to Cartesian gradient
C  EL      -  number of electrons over the grid (accumulated)
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
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),ISYM(NGen),IUNQ(NQ),
     $          XGRID(3,*),WGHT(*),BL(*),ExpMIN(NShell),NBAtm(NAtoms),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,NShell),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),
     $          GX(NBas,NBas),GY(NBas,NBas),GZ(NBas,NBas),
     $          GXB(NBas,NBas),GYB(NBas,NBas),GZB(NBas,NBas),
     $          NEqBAS(7,NBas),GWT(NAtoms,3),GC(3,NAtoms)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
      Logical rhf
c  ...................................................................
      Data nrad/400/, nang/434/     ! maximum values
c  ...................................................................
C
C
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  Main Driving Routine for Parallel Implementation of DFT Gradient
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      call getival('nsh',nsh)      ! needed because of Wolinski
C
C -- check if the Slaves are OK
C
      call para_check
C
C -- Send initial data to Slaves
C
      call para_initsend
      call para_pack_int(NRad,1)
      call para_pack_int(NAng,1)
      call para_pack_int(lrad,1)
      call para_pack_int(lang,1)
      call para_pack_int(IradQ,1)
      call para_pack_int(NBatch,1)
      call para_pack_int(IdWt,1)
c .........................................................
c -- these next lines needed because BASDAT array is
c -- corrupted in the slaves  (I blame Wolinski!)
      call para_pack_int(NShell,1)
      call para_pack_int(nsh,1)
      call para_pack_int(INX,12*NShell)
      call para_pack_real(BASDAT,13*nsh)
c ......................................................
      call para_pack_real(factor,1)
      call para_pack_real(thrsh,1)
      call para_bcast_pack(TxDftInit)
C
C -- wait awhile for Slaves to allocate memory
C
      nsq = NAtoms**2
      nsize = (NBas*(NBas+1))/2
      nbsq = NBas**2
c
      CALL ZeroIT(GX,nbsq)
      CALL ZeroIT(GY,nbsq)
      CALL ZeroIT(GZ,nbsq)
      If(.NOT.rhf) Then
        CALL ZeroIT(GXB,nbsq)
        CALL ZeroIT(GYB,nbsq)
        CALL ZeroIT(GZB,nbsq)
      EndIf
      If(IdWt.GT.0) Then
        CALL ZeroIT(GWT,3*NAtoms)
        CALL ZeroIT(GC,3*NAtoms)
      EndIf
      EL = 0.0d0
C
C -- Send second batch of data to Slaves
C
      call para_initsend
      call para_pack_int(NBAtm,NAtoms)
      call para_pack_int(NEqBAS,7*NBas)
      call para_pack_real(DA,nsize)
      If(.NOT.rhf) Then
        call para_pack_real(DB,nsize)
      EndIf
      call para_bcast_pack(TxDftDat)
C
C -- loop over atomic centres, sending each atom to a different slave
C
      call para_distr_atoms(iunq,nq)
c
c -- now sum up partial exchange-correlation gradient matrices
c -- and electron density over grid from each slave
      call para_reduce(GX,nbsq,TxDftF1)
      call para_reduce(GY,nbsq,TxDftF2)
      call para_reduce(GZ,nbsq,TxDftF3)
      call para_reduce(EL,1,TxDftF4)
      If(.NOT.rhf) Then
        call para_reduce(GXB,nbsq,TxDftF1b)
        call para_reduce(GYB,nbsq,TxDftF2b)
        call para_reduce(GZB,nbsq,TxDftF3b)
      EndIf
c
c -- That's it!
      call para_next(0)
c
      RETURN
      END
c ================================================================
c
      SUBROUTINE para_dftgiao(dft,    NAtoms, XNuc,   IAN,    NSym,
     *                        NGen,   ISYM,   NEqATM, NQ,     IUNQ,
     *                        IPRNT,  lrad,   lang,   IradQ,  factor,
     *                        NBatch, DISTN,  AIJ,    RDIST,  XXA,
     *                        WTA,    XGRID,  WGHT,   thrsh,  NBas,
     *                        NShell, BASDAT, INX,    BL,     ExpMIN,
     *                        NOcc,   CMO,    VCMO,   Malk,   NScr,
     *                        Z,      XC,     XMalkin)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Drives the slaves for parallel computation of DFT contribution to
C  NMR derivative Fock matrices (CLOSED SHELL)
C
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
C              7 - Becke 88 + Perdew-Wang 91
C              8 - Becke 88 + LYP (correlation)
C              9 - B3LYP   (VWN  + Becke 88 + LYP  + HF exchange)
C             10 - B3PW91  (VWN5 + Becke 88 + PW91 + HF exchange)
C             11 - Adiabatic connection method (user-defined)
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
C  lrad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  lang    -  angular quadrature order for numerical grid
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
C  NOcc    -  number of occupied MOs
C  CMO     -  occupied MO coefficients (transposed)
C  VCMO    -  virtual MO coefficients (transposed) used only if Malkin
C  Malk    -  method for doing malkin corrections
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  XC      -  3 deriv. DFT Exchange-Correlation matrix contibutions
C  XMalkin -  Malkin correction matrix
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),IUNQ(NQ),ISYM(NGen),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     $          CMO(NOcc,NBas),
     $          XC(3*NBas*(NBas+1)/2),
     $          VCMO(NBas-NOcc,NBas),XMalkin(NBas-NOcc,NOcc)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
      data nrad/400/, nang/434/
c
      call getival('ncphf',ncphf)
c -----------------------------------------------------------------
c
c -- check if the slaves are OK
      call para_check
c
c -- send data to slaves
      call para_initsend
      call para_pack_int(nrad,1)
      call para_pack_int(nang,1)
      call para_pack_int(lrad,1)
      call para_pack_int(lang,1)
      call para_pack_int(IradQ,1)
      call para_pack_int(NBatch,1)
      call para_pack_int(ncphf,1)
      call para_pack_int(Malk,1)
      call para_pack_real(factor,1)
      call para_pack_real(thrsh,1)
      call para_pack_real(CMO,NOcc*NBas)
      if(Malk.ne.0) call para_pack_real(VCMO,(NBas-NOcc)*NBas)
      call para_bcast_pack(TxDftNDat)
c
      nsize = 3*(NBas*(NBas+1))/2
c
c -- loop over atomic centres, sending each atom to a different slave
      call para_check
c
      DO 200 IAtom=1,NQ
      ICntr = IUNQ(IAtom)
      call para_recv(imsg,islave,TxDftNReq)
      call para_send(ICntr,islave,TxDftNAssign)
 200  CONTINUE
c
c -- all done  stop all slaves
      call elapsec(timo)
      call getrval('imbal',tlost)
      Do IAtom=1,nslv
      call para_recv(imsg,islave,TxDftNReq)
      call para_send(0,islave,TxDftNAssign)
      call elapsec(timb)
      tlost=tlost+(IAtom-1)*(timb-timo)
      timo=timb
      EndDo
      call setrval('imbal',tlost)
c
c -- now sum up exchange-correlation matrix and energy from each slave
      call para_reduce(XC,nsize,TxDftN)
      if(Malk.ne.0) call para_reduce(XMalkin,(NBas-NOcc)*NOcc,TxDftN2)
c
c -- That's it!
c
      RETURN
      END
c ================================================================
c
      SUBROUTINE Para_DFTH(dft,    NAtoms, ipass,  npass,  Nb,
     $                     Ne,     XNuc,   IAN,    nsym,   ngen,
     $                     ISYM,   NEqATM, NQ,     IUNQ,   IPRNT,
     $                     lrad,   lang,   IradQ,  factor, NBatch,
     $                     thrsh,  NBas,   NShell, BASDAT, INX,
     $                     NBAtm,  IdWt,   rhf,    DA,     DB,
     $                     DM,     FDA,    FDB,    HESS,   EL)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Drives the slaves for parallel computation of direct contribution to
C  Hessian matrix and contribution to fixed part of Fock derivative
C  matrices for CPHF
C
C    Jon Baker         June 2002
C                      November 2003    (ditto for UHF)
C    MM                January 2004     (multiple passes)
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
C  NAtoms  -  number of atoms
C  ipass   -  number of current pass (1 to npass)
C  npass   -  total number of passes needed to compute all Fock deriv.
C  Nb      -  first componet to be computed by current pass
C  Ne      -  last component to be computed by cureent pass
C  XNuc    -  nuclear coordinates
C  IAN     -  nuclear charge (atomic number)
C  NSym    -  number of symmetry operations
C             (currently highest abelian subgroup - excludes identity)
C  NGen    -  number of group generators
C  ISYM    -  list of symmetry operations (generators)
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of NQ symmetry-unique atoms
C  IPRNT   -  print flag
C  lrad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  lang    -  angular quadrature order for numerical grid
C             (if zero, will be determined in grid generation routines)
C  IradQ   -  radial quadrature type
C              0 - Euler-Maclaurin (Handy)          default
C              1 - Chebyshev second type (Becke mapping)
C  factor  -  determines  grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  NBatch  -  maximum number of grid points treated in one batch
C  thrsh   -  threshold for neglecting contribution
C  NBas    -  total number of basis functions
C  NShell  -  total number of shells
C  BASDAT  -  basis set data for TEXAS
C  INX     -  more basis set data for TEXAS
C  NBAtm   -  which atom each basis function belongs to (NBas)
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C  rhf     -  logical flag for RHF
C  DA      -  density matrix (closed-shell/alpha)
C  DB      -  density matrix (beta)
C  DM      -  maximum absolute density matrix element per row/column
C
C  on exit
C
C  FDA     -  derivative Fock matrices (alpha/closed-shell)
C  FDB     -  derivative Fock matrices (beta)
C  HESS    -  Hessian matrix
C  EL      -  number of electrons over the grid (accumulated)
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
C
C     "Analytic Second Derivatives of the Potential Energy Surface"
C      N.C.Handy, D.J.Tozer, G.J.Laming, C.W.Murray and R.D.Amos
C      Israel J.Chem.  33 (1993) 331
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),ISYM(NGen),IUNQ(NQ),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,NShell),
     $          NBAtm(NBas),DA(NBas*(NBas+1)/2),DB(*),DM(NBas),
     $          FDA(*),FDB(*),HESS(3*NAtoms,3*NAtoms)
      INTEGER dft
      Logical rhf
c  ...................................................................
      Data nrad/400/, nang/434/     ! maximum values
c  ...................................................................
C
C
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  Main Driving Routine for Parallel Implementation of DFT Hessian
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      call getival('nsh',nsh)      ! needed because of Wolinski
C
C -- check if the Slaves are OK
C
      call para_check
C
C -- Send initial data to Slaves
C    (only if this is the first pass)
C
      if(ipass.eq.1)then
        call para_initsend
        call para_pack_int(NRad,1)
        call para_pack_int(NAng,1)
        call para_pack_int(lrad,1)
        call para_pack_int(lang,1)
        call para_pack_int(IradQ,1)
        call para_pack_int(NBatch,1)
        call para_pack_int(IdWt,1)
c .........................................................
c -- these next lines needed because BASDAT array is
c -- corrupted in the slaves  (I blame Wolinski!)
        call para_pack_int(NShell,1)
        call para_pack_int(nsh,1)
        call para_pack_int(INX,12*NShell)
        call para_pack_real(BASDAT,13*nsh)
c ......................................................
        call para_pack_real(factor,1)
        call para_pack_real(thrsh,1)
        call para_bcast_pack(TxDftInit)
      endif
C
C -- wait awhile for Slaves to allocate memory
C
      ntri = (NBas*(NBas+1))/2
      nat3 = 3*NAtoms
c
C
C -- Send second batch of data to Slaves
C    (again, only in the first pass)
C
      if(ipass.eq.1)then
        EL = 0.0d0
        call para_initsend
        call para_pack_int(NBAtm,NBas)
        call para_pack_real(DA,ntri)
        If(.NOT.rhf) call para_pack_real(DB,ntri)
        call para_pack_real(DM,NBas)
        call para_bcast_pack(TxDftDat)
      endif
c
c  send total number of passes, current pass, first and last component
c
      natcurr=Ne-Nb+1
      call para_initsend
      call para_pack_int(npass,1)
      call para_pack_int(ipass,1)
      call para_pack_int(nb,1)
      call para_pack_int(ne,1)
      call para_bcast_pack(TxDftHmp)
C
C -- loop over atomic centres, sending each atom to a different slave
C
      call para_distr_atoms(iunq,nq)
c
c -- now sum up partial derivative Fock matrices, Hessian matrices
c -- and electron density over grid from each slave
      call para_reduce(HESS,nat3**2,TxHess4)
      call para_reduce(EL,1,TxDftF4)
c --------------------------------------------------------------------
      call para_reduce(FDA,3*ntri*natcurr,TxHess2)
c -- now accumulate beta Fock derivative matrices
      If(.NOT.rhf) call para_reduce(FDB,3*ntri*natcurr,TxHess3)
c ---------------------------------------------------------------------
c -- That's it!
cc      call para_next(0)
c
      RETURN
      END
c ================================================================
c
      SUBROUTINE Para_DFTCphf(dft,    NAtoms, NatOnce,NatDo,  XNuc,
     $                        XSym,   IAN,    NQ,     IUNQ,   lrad,
     $                        lang,   IradQ,  factor, NBatch, thrsh,
     $                        NBas,   NShell, BASDAT, INX,    IdWt,
     $                        rhf,    DA,     DB,     DEN1,   DEN1B,
     $                        DM1,    FDA,    FDB,    lgrid)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Drives the slaves for parallel computation of direct contribution to
C  Hessian matrix and contribution to fixed part of Fock derivative
C  matrices for CPHF
C
C    Jon Baker         June 2002
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
C  NAtoms  -  total number of symmetry-unique atoms
C  NatOnce -  number of atoms on this pass
C  NatDo   -  number of atoms for which the CPHF has not converged
C  XNuc    -  atomic coordinates
C  XSym    -   ditto symmetry-unique atoms
C  IAN     -  nuclear charge (atomic number) for symmetry-unique atoms
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of NQ symmetry-unique atoms
C  lrad    -  number of radial grid points for numerical grid
C             (if zero, will be determined in grid generation routines)
C  lang    -  angular quadrature order for numerical grid
C             (if zero, will be determined in grid generation routines)
C  IradQ   -  radial quadrature type
C              0 - Euler-Maclaurin (Handy)          default
C              1 - Chebyshev second type (Becke mapping)
C  factor  -  determines  grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  NBatch  -  maximum number of grid points treated in one batch
C  thrsh   -  threshold for neglecting contribution
C  NBas    -  total number of basis functions
C  NShell  -  total number of shells
C  BASDAT  -  basis set data for TEXAS
C  INX     -  more basis set data for TEXAS
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C  rhf     -  logical flag for RHF
C  DA      -  density matrix (closed-shell/alpha)
C  DB      -  density matrix (beta)
C  DEN1    -  current solution to CPHF, lower triangle (alpha/closed-shell)
C  DEN1B   -  current solution to CPHF, lower triangle (beta)
C  DM1     -  maximum absolute 1st-order density matrix element per row/column
C
C  on exit
C
C  FDA     -  derivative Fock matrices (alpha/closed-shell)
C  FDB     -  derivative Fock matrices (beta)
C  HESS    -  Hessian matrix
C  EL      -  number of electrons over the grid (accumulated)
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
C
C     "Analytic Second Derivatives of the Potential Energy Surface"
C      N.C.Handy, D.J.Tozer, G.J.Laming, C.W.Murray and R.D.Amos
C      Israel J.Chem.  33 (1993) 331
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),IUNQ(NQ),XSym(3,NQ),
     $          BASDAT(13,*),INX(12,NShell),
     $          DA(NBas*(NBas+1)/2),DB(*),DM1(NBas),
     $          DEN1(3,NatOnce,*),DEN1B(3,NatOnce,*),FDA(*),FDB(*)
      INTEGER dft
      LOGICAL rhf
c  ...................................................................
      Data nrad/400/, nang/434/     ! maximum values
c  ...................................................................
C
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  Main Driving Routine for Parallel Implementation of DFT Hessian CPHF
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      call getival('nsh',nsh)      ! needed because of Wolinski or Pulay
c                                    each of them blames the other
C
C -- check if the Slaves are OK
C
      call para_check
C
C --  NO NEED TO SEND ANY INITIAL DFT DATA TO SLAVES AS SHOULD
C --  ALREADY BE PRESENT FOLLOWING PREVIOUS DFT STEP
C
C -- wait awhile for Slaves to allocate memory
C
      ntri = (NBas*(NBas+1))/2
      nat3 = 3*NatOnce
C
C -- Send second batch of data to Slaves
C
      lgrid = 1     ! temporary (free-for-all)    JB
      call para_initsend
      call para_pack_int(lgrid,1)
      call para_pack_int(NatDo,1)
      call para_pack_real(DM1,NBas)
      call para_pack_real(BASDAT,13*nsh)
      call para_pack_int(NQ,1)
      call para_pack_real(XSym,3*NQ)
      call para_pack_int(IAN,NQ)
      call para_pack_real(DA,ntri)
      If(.NOT.rhf) Then
        call para_pack_real(DB,ntri)
      EndIf
      call para_bcast_pack(TxDftDat)
C
C -- loop over atomic centres, sending each atom to a different slave
C
      call para_distr_atoms(iunq,nq)
c
c
c summing up is done together with the normal integrals      
c ---------------------------------------------------------------------
c -- That's it!
cc      call para_next(0)
c
      RETURN
      END
c ================================================================
c
      SUBROUTINE WhichSlave(islave,NQ,NSLAVE,LSlave,JAtom)
      IMPLICIT INTEGER(A-Z)
C
C  Determines which atom should be sent to current available slave
C  for controlled parallelism
C
C  ARGUMENTS
C
C  islave  -  slave ID
C  NQ      -  number of symmetry-unique atoms
C  NSLAVE  -  array matching atoms to slaves
C  LSlave  -  logical array  .true. - is atom has been done
C                           .false. -  otherwise
C  JAtom   -  on exit the atomic centre to send to slave islave
C              -1 if all atoms done for that slave
C
C
      DIMENSION NSLAVE(NQ)
      LOGICAL   LSLAVE(NQ)
C
C
C  loop over all atoms
C
      DO 10 IAtom=1,NQ
      If(NSLAVE(IAtom).EQ.islave.AND..NOT.LSlave(IAtom)) Then
        JAtom = IAtom
        LSlave(IAtom) = .True.
        RETURN
      EndIf
 10   CONTINUE
C
C  no match found
C
      JAtom = -1
C
      RETURN
      END
