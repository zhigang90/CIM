      SUBROUTINE DFTHessU(dft,    ICntr,  NAtoms, Nb,     Ne,
     $                    XNuc,   IAN,    NSym,   NGen,   ISYM,
     $                    NEqATM, IPRNT,  NRad,   NAng,   IradQ,
     $                    factor, NBatch, DISTN,  AIJ,    RDIST,
     $                    XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                    NBas,   NShell, BASDAT, INX,    BL,
     $                    ExpMIN, IdWt,   NBAtm,  DA,     DB,
     $                    DM,     NScr,   Z,      FDA,    FDB,
     $                    HESS,   EL,     IEntry)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates DFT contribution to Hessian and derivative Fock matrices
C  ** OPEN SHELL **
C
C    MM September 2003 (based on DFTHessC)
C
C    Multiple passes         MM December 2003
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
C  NBAtm   -  which atom each basis function belongs to (NBas)
C  DA      -  alpha density matrix, lower triangle
C  DB      -  beta density matrix, lower triangle
C  DM      -  maximum density matrix element per column
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  FDA     -  alpha derivative Fock matrices
C  FDB     -  beta derivative Fock matrices
C  HESS    -  Hessian matrix
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
C     "Analytical Second Derivatives of the Potential Energy Surface"
C      N.C.Handy, D.J.Tozer, G.J.Laming, C.W.Murray and R.D.Amos
C      Israel J. Chem.  33 (1993) 331
C
C     "An implementation of analytic second derivatives of the
C      gradient-corrected density functional energy",
C      B.G. Johnson and M.J. Frisch,
C      J. Chem. Phys. 100 (1994) 7429 (*)
C
C     (*) NOTE: this paper has two errors/misprints, one in Fig.2 and
C               one in Fig.3
C
C  ----------------------------------------------------------------
C
C   Multiple passes inplementation:
C
C   The input parameters Nb and Ne hold the first and
C   last component of the derivative Fock matrix that
C   are to be computed in this pass.
C
C   There can be 3 cases:
C
C   Case 1:  one pass only: NAtonce=Ne-Nb+1=NAtoms
C            All contribution to Fock derivatives and Hessian
C            are computed at the same time. With weight derivatives,
C            The contribution of to the ICntr component of Fock
C            derivatives and Hessian is obtained by translational
C            invariance at the end of the quadrature.
C
C   Case 2:  Natonce=Ne-Nb+1 < Natoms, Ne <= Icntr <= Nb
C            Multiple passes, the Icntr component is computed
C            in the current pass. In this case we perform quadrature
C            for Natonce fock derivatives and all the Hessian at
C            the same time. With weight derivatives, the contribution
C            to Fock derivative Icntr is obtained applying
C            translational invariance at each grid point. The
C            contributions to the Icntr components of the Hessian
C            are computed at the end of the quadrature, as in Case 1
C
C   Case 3:  Natonce=Ne-Nb+1 < NAtoms, Icntr <= Nb Or ICntr >= Ne
C            Multiple passes, the Icntr component is not computed
C            in the current pass. In this case we perform quadrature
C            for NAtonce Fock derivatives only. The Hessian is not
C            computed at all, as it will be computed by the pass
C            belonging to Case 2. With weight derivatives, there is
C            no need for translational invariance in the current pass.
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),ISYM(NGen),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     $          DA(*),DB(*),DM(NBas),FDA(3*NAtoms*NBas*(NBas+1)/2),
     $          FDB(3*NAtoms*NBas*(NBas+1)/2),
     $          NBAtm(NBas),HESS(9*NAtoms**2)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
      CHARACTER*256 Filename,scrfile
      CHARACTER*3 ch3
C
      PARAMETER (Zero=0.0d0,Two=2.0d0)
      save mpoint,npptot
      save td1,tqf,tg1,tsw,tqh
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
        tdwt = Zero
c
        mpoint=0
        npptot=0
        td1=zero
        tg1=zero
        tsw=zero
        tqf=zero
        tqh=zero
cc
        IEntry = 0
      ENDIF
c
c compute wich case the current pass belongs to (see initial comments)
c
      NAtonce=Ne-Nb+1
      if(NAtonce.eq.NAtoms)then
        Icase=1  ! One pass only
      else if(Icntr.ge.Nb.and.ICntr.le.Ne)then
        Icase=2  ! Multiple passes, compute ICntr
      else
        Icase=3  ! Multiple passes, do not compute ICntr
      endif
c
c   sanity check
c
      if(NAtonce.le.0.or.NAtonce.gt.Natoms)then
        call nerror(1,'DFTHessU','something is terrybly wrong',
     $              natonce,natoms)
      endif
c
      ntri=nbas*(nbas+1)/2
      nb3=3*natoms*ntri
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
      iao = 1                      ! AO values over grid points
      iaox = iao + NBasS           ! AO deriv. over grid points
      iaoxx = iaox + NBasT         ! AO 2nd deriv. over grid points
      inb = iaoxx + NBasV          ! index array for "non-zero" AOs
      ibf = inb + NBasS            ! number of "non-zero" AOs per point
      ida  = ibf + NBatch+1        ! alpha density per grid point
      idb = ida  + NBatch          ! beta density per grid point
      idxa = idb  + NBatch         ! atomic grad. of alpha dens.
      idxb = idxa + 3*NAtoms       ! atomic grad. of beta dens.
      idxxa = idxb + 3*NAtoms      ! atomic Hess. of alpha dens.
      idxxb = idxxa + 9*NAtoms**2  ! atomic Hess. of beta dens.
      ipra  = idxxb + 9*NAtoms**2  ! Func. deriv. w.r.t. alpha dens.
      iprb = ipra + NBatch         ! Func. deriv. w.r.t. beta dens.
      iprara = iprb + NBatch       ! Funct. 2nd deriv. w.r.t. alpha dens.
      iprbrb = iprara + NBatch     ! Funct. 2nd deriv. w.r.t. beta dens.
      iprarb = iprbrb + NBatch     ! 2nd deriv. w.r.t. alpha beta dens.
      isva = iprarb + NBatch       ! alpha V coeff. in quadrature
      isvb = isva + 3*NAtoms       ! beta V coeff. in quadrature
      ind = isvb + 3*NAtoms        ! index of "non-zero" densities
      iscr = ind + NBatch          ! max AO value per grid point
      Iend = iscr + NBatch -1
cc
      i1 = inb                     ! scratch pointers for grid
      i2 = i1 + NAtoms
      i3 = i2 + NAtoms
cc
      If(dft.GT.3) Then
C
C  additional scratch space for non-local functionals
C
        iaoxxx=1+IEnd              !AO 3rd deriv. at grid points
        idga=iaoxxx+10*NBasS       !alpha density gradient
        idgb=idga+3*NBatch         !beta density gradient
        igdxa=idgb+3*NBatch        !in situ alpha atomic grad. of dens. grad.
        igdxb=igdxa+9*NAtoms       !in situ beta atomic grad. of dens. grad.
        igdxxa=igdxb+9*NAtoms      !in situ alpha atomic Hess. of dens. grad.
        igdxxb=igdxxa+27*NAtoms**2 !in situ beta atomic Hess. of dens. grad.
        igxaa=igdxxb+27*NAtoms**2  !in situ alpha atomic grad. of grad. invariant
        igxbb=igxaa+3*NAtoms       !in situ beta atomic grad. of grad. invariant
        igxab=igxbb+3*NAtoms       !in situ alpha beta atomic grad. of grad. invariant
        igxxaa=igxab+3*NAtoms      !in situ alpha atomic Hess. of grad. invariant
        igxxbb=igxxaa+9*NAtoms**2  !in situ beta atomic Hess. of grad. invariant
        igxxab=igxxbb+9*NAtoms**2  !in situ alpha beta atomic Hess. of grad. invariant
        iswa=igxxab+9*NAtoms**2    !alpha W coeff. in numerical quadrature
        iswb=iswa+9*NAtoms         !beta W coeff. in numerical quadrature
        isga=iswb+9*NAtoms         !partial sum for quadrature
        isgb=isga+3*NAtoms         !partial sum for quadrature
        isgc=isgb+3*NAtoms         !partial sum for quadrature
        ipga=isgc+3*NAtoms         !Func. deriv. w.r.t. alpha grad.
        ipgb=ipga+NBAtch           !Func. deriv. w.r.t. beta grad.
        ipgc=ipgb+NBAtch           !Func. deriv. w.r.t. alpha beta grad.
        ipraga=ipgc+NBAtch         !Func. 2nd deriv. w.r.t. alpha dens. and grad.
        ipragb=ipraga+NBAtch       !Func. 2nd deriv. w.r.t alpha dens. and beta grad.
        ipragc=ipragb+NBAtch       !Func. 2nd deriv. w.r.t. alpha dens and alpha beta grad.
        iprbga=ipragc+NBAtch       !func. 2nd deriv. w.r.t beta dens. and alpha grad.
        iprbgb=iprbga+NBAtch       !Func. 2nd deriv. w.r.t. beta dens. and grad.
        iprbgc=iprbgb+NBAtch       !Func. 2nd deriv. w.r.t. beta dens and alpha beta grad.
        ipgaga=iprbgc+NBAtch       !Func. 2nd. deriv. w.r.t. alpha grad.
        ipgagb=ipgaga+NBAtch       !Func. 2nd. deriv. w.r.t. alpha and beta grad.
        ipgagc=ipgagb+NBAtch       !Func. 2nd. deriv. w.r.t. alpha and alpha beta grad.
        ipgbgb=ipgagc+NBAtch       !Func. 2nd. deriv. w.r.t. beta grad.
        ipgbgc=ipgbgb+NBAtch       !Func. 2nd. deriv. w.r.t. beta and alpha beta grad.
        ipgcgc=ipgbgc+NBAtch       !Func. 2nd. deriv. w.r.t. alpha beta grad.
        IEnd=ipgcgc+NBAtch -1
      EndIf
cc
      if(IdWt.eq.1)then
c
c  scratch space for weight derivatives
c
        iexc  = 1 + IEnd        ! Exchange-correlation energy per grid point
        igwt  = iexc + NBatch   ! Gradient of quadrature weight per grid point
        ihwt = igwt + NBatch*3*NAtoms ! Hessian of quadrature weight per grid point
        icell = ihwt + NBatch*9*NAtoms*NAtoms ! Cell functions
        iauxt = icell+ NAtoms          ! Auxiliary function t
        iauxq = iauxt+ NAtoms*NAtoms   ! Auxiliary function q
        igmu  = iauxq+ NAtoms*NAtoms   ! Gradient of adjusted hyperbolic coordinates
        ihmu  = igmu+ 3*NAtoms*NAtoms  ! Hessian of adjusted hyperbolic coordinates
        IEnd = ihmu + 15*NAtoms*NAtoms -1
      endif
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'DFTHessU')
C
C
C  Begin calculation for this centre
C  ---------------------------------
C
c
c    if using weight derivatives, save current hessian and
c    Fock matrices to file, according to the value of icase
c
      if(idwt.eq.1)then
        Call getchval('scrf',scrfile)
        Call rmblan(scrfile,256,len)
        Write(ch3,'(I3)') ICntr
        If(ICntr.LT.100) ch3(1:1) = '0'
        If(ICntr.LT.10)  ch3(2:2) = '0'
c
        if(Icase.eq.1)then
          Filename = scrfile(1:len)//'.fda.'//ch3
          len1 = len+8
          OPEN (UNIT=41,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='UNKNOWN')
          do i=0,NAtoms-1
            ii = 1 + 3*ntri*i
            Call WriteBinary1(41,3*ntri,FDA(ii))
          enddo
          CLOSE(UNIT=41,STATUS='KEEP')
          call zeroit(fda,nb3)
c
          Filename = scrfile(1:len)//'.fdb.'//ch3
          len1 = len+8
          OPEN (UNIT=41,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='UNKNOWN')
          do i=0,NAtoms-1
            ii = 1 + 3*ntri*i
            Call WriteBinary1(41,3*ntri,FDB(ii))
          enddo
          CLOSE(UNIT=41,STATUS='KEEP')
          call zeroit(fdb,nb3)
        endif
c
        if(Icase.le.2)then
          Filename = scrfile(1:len)//'.hess.'//ch3
          len1 = len+9
          OPEN (UNIT=41,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='UNKNOWN')
          WRITE(41) HESS
          CLOSE(UNIT=41,STATUS='KEEP')
          call zeroit(HESS,9*NAtoms**2)
        endif
      endif
      call secund(t1)
      CALL GRID(NAtoms, XNuc,   IAN,    IPRNT,  NRad,
     $          NAng,   factor, IradQ,  Z(i1),  Z(i2),
     $          Z(i3),  DISTN,  RDIST,  AIJ,    ICntr,
     $          XXA,    WTA,    NPoint, XGRID,  WGHT)
C
C  ==============================================================
C  ** SYMMETRY SORT OF GRID **
      If(NSym.GT.0) Then
        XX = XNuc(1,ICntr)
        YY = XNuc(2,ICntr)
        ZZ = XNuc(3,ICntr)
        CALL SymGRID(NAtoms, ICntr,  XX,     YY,     ZZ,
     $               NSym,   NGen,   ISYM,   NEqATM, NPoint,
     $               XGRID,  WGHT)
      EndIf
C  ===============================================================
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
      IP = 1
      NLeft = NPoint
      mpoint=mpoint+npoint
 50   CONTINUE
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
      CALL zeroit(Z(iaoxx),NBasV)
      If(dft.LE.3) Then
        CALL AODer2(NP,     NShell, NBas,   XGRID(1,IP),
     $              XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $              ExpMIN, Z(iao), Z(iaox),Z(iaoxx))
        CALL AOXXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $               Z(iaoxx),Z(ibf), Z(inb), Z(iao), Z(iaox),
     $               Z(iaoxx),Z(iscr))
      Else
        CALL zeroit(Z(iaoxxx),10*NBasS)
        CALL AODer3(NP,     NShell, NBas,   XGRID(1,IP),
     $              XNuc,   BL,     BASDAT, INX,        ExpCUT,
     $              ExpMIN, Z(iao), Z(iaox),Z(iaoxx),   Z(iaoxxx))
        CALL AOXXXSort(NP,      NBas,     thint,    Z(iao), Z(iaox),
     $                 Z(iaoxx),Z(iaoxxx),Z(ibf),   Z(inb), Z(iao),
     $                 Z(iaox), Z(iaoxx), Z(iaoxxx),Z(iscr))
      EndIf
      call secund(t2)
      tao = tao + t2-t1
C
C  form alpha/beta density at current grid point
C
      call secund(t1)
      CALL GetDENSU(NP,     NBas,   Z(ibf), DA,   DB,
     $              Z(iao), Z(inb), Z(idA), Z(idB))
      call secund(t2)
      tden = tden + t2-t1
C
C  form density derivatives if non-local
C
      If(dft.GT.3) Then
        call secund(t1)
c
c  use the ipra storage to compute the total (alpha + beta)
c  density (the threshold is on density sum).
c
        CALL AddVEC(NP,Z(idA),Z(idB),Z(ipra))
c
        CALL GetDENDerivU(NP,     NBas,   thrsh,  z(ibf), DA,
     $                    DB,     Z(iao), z(iaox),z(inb), z(ipra),
     $                    z(idga),z(idgb))
        call secund(t2)
        tdenx = tdenx + t2-t1
      EndIf
C
C
C  calculate potentials (and potential derivatives)
C
      call secund(t1)
      CALL VXCPOTDU(dft,     NP,       thrsh,    Z(idA),   Z(idB),
     $             Z(idga),  Z(idgb),  WGHT(IP), Z(ipra),  Z(iprb),
     $             Z(iprara),Z(iprbrb),Z(iprarb),Z(ipga),  Z(ipgb),
     $             Z(ipgc),  Z(ipraga),Z(ipragb),Z(ipragc),Z(iprbga),
     $             Z(iprbgb),Z(iprbgc),Z(ipgaga),Z(ipgagb),Z(ipgagc),
     $             Z(ipgbgb),Z(ipgbgc),Z(ipgcgc),EXC,      EL0,
     $             Z(iexc),  IdWt)
c
c  update EL only in cases 1 and 2
c
      if(icase.le.2)EL=EL+EL0
c ...................................................................
c  If the weight derivatives section is active (IdWt=1),
c  on exit from VXCPOTD, Z(iexc) contains the exchange-correlation
c  energy contribution AT EACH OF NP GRID POINTS
c ...................................................................
C
C  determine number of contributing points
C
      call potsorthu(dft,      np,       thint,    z(ipra),  z(iprb),
     $               z(iprara),z(iprbrb),z(iprarb),z(ipga),  z(ipgb),
     $               z(ipgc),  z(ipraga),z(ipragb),z(ipragc),z(iprbga),
     $               z(iprbgb),z(iprbgc),z(ipgaga),z(ipgagb),z(ipgagc),
     $               z(ipgbgb),z(ipgbgc),z(ipgcgc),z(ind),   npp)
      npptot=npptot+npp
      call secund(t2)
      tpot = tpot + t2-t1
c
c  now we have to chose the appropriate routine for numerical
c  quadrature, according to Icase
c
      if(icase.eq.1)then       ! one pass only
        if(IdWt.eq.1)then
          call secund(t1)      ! weight derivatives
          CALL zeroit(Z(igwt),3*NAtoms*NBatch)
          CALL zeroit(Z(ihwt),9*NAtoms*NAtoms*NBatch)
          CALL GHWT(NPP,     XGRID(1,IP),WGHT(IP),Z(ind),  ICntr,
     $              NAtoms,  XNuc,       AIJ,     RDIST,   Z(icell),
     $              Z(iauxt),Z(iauxq),   Z(igmu), Z(ihmu), Z(iexc),
     $              Z(igwt),    Z(ihwt))
          call secund(t2)
          tdwt = tdwt + t2-t1
          call secund(t1) ! one pass quadrature with weight derivatives
          call cwopfhu(dft,    npp,      nbas,     z(ibf),   natoms,
     $               nbatm,    thrsh,    da,       db,       dm,
     $               z(idga),  z(idgb),  wght(ip), z(ipra),  z(iprb),
     $               z(iprara),z(iprbrb),z(iprarb),z(ipga),  z(ipgb),
     $               z(ipgc),  z(ipraga),z(ipragb),z(ipragc),z(iprbga),
     $               z(iprbgb),z(iprbgc),z(ipgaga),z(ipgagb),z(ipgagc),
     $               z(ipgbgb),z(ipgbgc),z(ipgcgc),z(iao),   z(iaox),
     $               z(iaoxx), z(iaoxxx),z(inb),   z(iscr),  z(ind),
     $               z(idxa),  z(idxb),  z(idxxa), z(idxxb), z(igdxa),
     $               z(igdxb), z(igdxxa),z(igdxxb),z(igxaa), z(igxbb),
     $               z(igxab), z(igxxaa),z(igxxbb),z(igxxab),z(isva),
     $               z(isvb),  z(iswa),  z(iswb),  z(isga),  z(isgb),
     $               z(isgc),  icntr,    z(igwt),  z(ihwt),  fda,
     $               fdb,      hess,     td1,      tg1,      tsw,
     $               tqf,      tqh)
        else      ! one pass quadrature without weight derivatives
          call secund(t1)
          call mkmpfhu(dft,  npp,      nbas,     z(ibf),    natoms,
     $             nb,       ne,       nbatm,    thrsh,     da,
     $             db,       dm,       z(idga),  z(idgb),   wght(ip),
     $             z(ipra),  z(iprb),  z(iprara),z(iprbrb), z(iprarb),
     $             z(ipga),  z(ipgb),  z(ipgc),  z(ipraga), z(ipragb),
     $             z(ipragc),z(iprbga),z(iprbgb),z(iprbgc), z(ipgaga),
     $             z(ipgagb),z(ipgagc),z(ipgbgb),z(ipgbgc), z(ipgcgc),
     $             z(iao),   z(iaox),  z(iaoxx), z(iaoxxx), z(inb),
     $             z(iscr),  z(ind),   z(idxa),  z(idxb),   z(idxxa),
     $             z(idxxb), z(igdxa), z(igdxb), z(igdxxa), z(igdxxb),
     $             z(igxaa), z(igxbb), z(igxab), z(igxxaa), z(igxxbb),
     $             z(igxxab),z(isva),  z(isvb),  z(iswa),   z(iswb),
     $             z(isga),  z(isgb),  z(isgc),  icntr,     fda,
     $             fdb,      hess,     td1,      tg1,       tsw,
     $             tqf,      tqh)
        endif
      else if(icase.eq.2)then ! multiple passes, compute Icntr
        if(IdWt.eq.1)then
          call secund(t1)     ! weight derivatives
          CALL zeroit(Z(igwt),3*NAtoms*NBatch)
          CALL zeroit(Z(ihwt),9*NAtoms*NAtoms*NBatch)
          CALL GHWT(NPP,     XGRID(1,IP),WGHT(IP),Z(ind),  ICntr,
     $              NAtoms,  XNuc,       AIJ,     RDIST,   Z(icell),
     $              Z(iauxt),Z(iauxq),   Z(igmu), Z(ihmu), Z(iexc),
     $              Z(igwt),    Z(ihwt))
          call secund(t2)
          tdwt = tdwt + t2-t1
          call secund(t1) ! n passess, Icntr, weight derivatives
          call cwmpfhu(dft,    npp,      nbas,     z(ibf),   natoms,
     $               natonce,  nb,       ne,       nbatm,    thrsh,
     $               da,       db,       dm,       z(idga),  z(idgb),
     $               wght(ip), z(ipra),  z(iprb),  z(iprara),z(iprbrb),
     $               z(iprarb),z(ipga),  z(ipgb),  z(ipgc),  z(ipraga),
     $               z(ipragb),z(ipragc),z(iprbga),z(iprbgb),z(iprbgc),
     $               z(ipgaga),z(ipgagb),z(ipgagc),z(ipgbgb),z(ipgbgc),
     $               z(ipgcgc),z(iao),   z(iaox),  z(iaoxx), z(iaoxxx),
     $               z(inb),   z(iscr),  z(ind),   z(idxa),  z(idxb),
     $               z(idxxa), z(idxxb), z(igdxa), z(igdxb), z(igdxxa),
     $               z(igdxxb),z(igxaa), z(igxbb), z(igxab), z(igxxaa),
     $               z(igxxbb),z(igxxab),z(isva),  z(isvb),  z(iswa),
     $               z(iswb),  z(isga),  z(isgb),  z(isgc),  icntr,
     $               z(igwt),  z(ihwt),  z(igmu),  z(ihmu),  fda,
     $               fdb,      hess,     td1,      tg1,      tsw,
     $               tqf,      tqh)
        else
          call secund(t1) ! n passess, Icntr, no weight derivatives
          call mkmpfhu(dft,  npp,      nbas,     z(ibf),    natoms,
     $             nb,       ne,       nbatm,    thrsh,     da,
     $             db,       dm,       z(idga),  z(idgb),   wght(ip),
     $             z(ipra),  z(iprb),  z(iprara),z(iprbrb), z(iprarb),
     $             z(ipga),  z(ipgb),  z(ipgc),  z(ipraga), z(ipragb),
     $             z(ipragc),z(iprbga),z(iprbgb),z(iprbgc), z(ipgaga),
     $             z(ipgagb),z(ipgagc),z(ipgbgb),z(ipgbgc), z(ipgcgc),
     $             z(iao),   z(iaox),  z(iaoxx), z(iaoxxx), z(inb),
     $             z(iscr),  z(ind),   z(idxa),  z(idxb),   z(idxxa),
     $             z(idxxb), z(igdxa), z(igdxb), z(igdxxa), z(igdxxb),
     $             z(igxaa), z(igxbb), z(igxab), z(igxxaa), z(igxxbb),
     $             z(igxxab),z(isva),  z(isvb),  z(iswa),   z(iswb),
     $             z(isga),  z(isgb),  z(isgc),  icntr,     fda,
     $             fdb,      hess,     td1,      tg1,       tsw,
     $             tqf,      tqh)
        endif
      else                    ! multiple passes, no Icntr
        if(IdWt.eq.1)then
          call secund(t1)     ! weight derivatives (gradient only)
          CALL zeroit(Z(igwt),3*NAtoms*NBatch)
          CALL GWTONLY(NPP,     XGRID(1,IP),WGHT(IP),Z(ind),  ICntr,
     $                 NAtoms,  XNuc,       AIJ,     RDIST,   Z(iauxq),
     $                 Z(icell),Z(iauxt),   Z(igmu), Z(igwt))
          call secund(t2)
          tdwt = tdwt + t2-t1
          call secund(t1) ! n passess, no Icntr, weight derivatives
          call cwmpfu(dft,    npp,      nbas,     z(ibf),   natoms,
     $              natonce,  nb,       ne,       nbatm,    thrsh,
     $              da,       db,       dm,       z(idga),  z(idgb),
     $              wght(ip), z(ipra),  z(iprb),  z(iprara),z(iprbrb),
     $              z(iprarb),z(ipga),  z(ipgb),  z(ipgc),  z(ipraga),
     $              z(ipragb),z(ipragc),z(iprbga),z(iprbgb),z(iprbgc),
     $              z(ipgaga),z(ipgagb),z(ipgagc),z(ipgbgb),z(ipgbgc),
     $              z(ipgcgc),z(iao),   z(iaox),  z(iaoxx), z(inb),
     $              z(iscr),  z(ind),   z(idxa),  z(idxb),  z(igdxa),
     $              z(igdxb), z(igxaa), z(igxbb), z(igxab), z(isva),
     $              z(isvb),  z(iswa),  z(iswb),  z(isga),  z(isgb),
     $              z(isgc),  icntr,    z(igwt),  fda,      fdb,
     $              td1,      tg1,      tsw,      tqf)
        else
          call secund(t1) ! n passess, no Icntr, no weight derivatives
          call mkmpfu(dft,  npp,      nbas,     z(ibf),    natoms,
     $            nb,       ne,       nbatm,    thrsh,     da,
     $            db,       dm,       z(idga),  z(idgb),   wght(ip),
     $            z(ipra),  z(iprb),  z(iprara),z(iprbrb), z(iprarb),
     $            z(ipga),  z(ipgb),  z(ipgc),  z(ipraga), z(ipragb),
     $            z(ipragc),z(iprbga),z(iprbgb),z(iprbgc), z(ipgaga),
     $            z(ipgagb),z(ipgagc),z(ipgbgb),z(ipgbgc), z(ipgcgc),
     $            z(iao),   z(iaox),  z(iaoxx), z(inb),    z(iscr),
     $            z(ind),   z(idxa),  z(idxb),  z(igdxa), z(igdxb),
     $            z(igxaa), z(igxbb), z(igxab), z(isva),  z(isvb),
     $            z(iswa),  z(iswb),  z(isga),  z(isgb),  z(isgc),
     $            icntr,    fda,      fdb,      td1,      tg1,
     $            tsw,      tqf)
        endif
      endif
c
      call secund(t2)
      tnum = tnum + t2-t1
c
      IP = IP+NP
      NLeft = NLeft-NP
      If(NLeft.GT.0) GO TO 50
C
C  end loop over grid points
C
C  ................................................................
C
c
c  if using weight derivatives, apply translational invariance
c  according to icase, and add current contribution to the existing one
c
      if(idwt.eq.1)then
        if(Icase.eq.1)then   ! Fock derivative
          call trinvaf(fda,natoms,nbas,icntr)
          Filename = scrfile(1:len)//'.fda.'//ch3
          len1 = len+8
          OPEN (UNIT=41,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='OLD')
          do i=0,NAtoms-1
            ii = 1 + 3*ntri*i
            CALL ReadBinary1(41,3*ntri,Z)
            CALL AddVec(3*ntri,FDA(ii),Z,FDA(ii))
          enddo
          CLOSE(UNIT=41,STATUS='DELETE')
c
          call trinvaf(fdb,natoms,nbas,icntr)
          Filename = scrfile(1:len)//'.fdb.'//ch3
          len1 = len+8
          OPEN (UNIT=41,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='OLD')
          do i=0,NAtoms-1
            ii = 1 + 3*ntri*i
            CALL ReadBinary1(41,3*ntri,Z)
            CALL AddVec(3*ntri,FDB(ii),Z,FDB(ii))
          enddo
          CLOSE(UNIT=41,STATUS='DELETE')
        endif
c
        if(Icase.le.2)then   ! Hessian
          call trinvah(hess,natoms,icntr)
          Filename = scrfile(1:len)//'.hess.'//ch3
          len1 = len+9
          OPEN (UNIT=41,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='OLD')
          CALL ReadBinary1(41,9*NAtoms**2,Z)
          CLOSE(UNIT=41,STATUS='DELETE')
          call addhess(HESS,Z,natoms)
        endif
      endif
c
      If(IEntry.EQ.-1.AND.IPrnt.GT.0)then
        write(6,2222) tgrid,tao,tden,tdenx,tpot,tdwt,tnum,
     $                td1,tg1,tsw,tqf,tqh,
     $                mpoint,npptot
      endif
C
      RETURN
c
 2222   Format(/,' Time taken to calculate grid:   ',f7.2,/,
     $           ' Time to calculate & sort AOs:   ',f7.2,/,
     $           ' Time to calculate density:      ',f7.2,/,
     $           ' Time to calculate den.derivs:   ',f7.2,/,
     $           ' Time to evaluate potential:     ',f7.2,/,
     $           ' Time weight derivatives:        ',f7.2,/,
     $           ' Time numerical integration:     ',f7.2,/,
     $           ' of which:',/,
     $           '      for density derivatives:   ',f7.2,/,
     $           '      for gradient invariant:    ',f7.2,/,
     $           '      for coefficients s and w   ',f7.2,/,
     $           '      for quadrature Fock:       ',f7.2,/,
     $           '      for quadrature Hessian:    ',f7.2,//
     $           ' Total number of grid points:    ',I8,/,
     $           ' Number of contributing points:  ',I8,/)
c
      END
