      SUBROUTINE DFTCphfC(dft,    ICntr,  NAtoms, NatOnce,NatDo,
     $                    XNuc,   XSym,   IAN,    IPRNT,  NRad,
     $                    NAng,   IradQ,  factor, NBatch, DISTN,
     $                    AIJ,    RDIST,  XXA,    WTA,    XGRID,
     $                    WGHT,   thrsh,  NBas,   NShell, BASDAT,
     $                    INX,    BL,     ExpMIN, IdWt,   DA,
     $                    DEN1,   DM1,    NScr,   Z,      FDA,
     $                    lgrid,  IEntry)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates DFT contribution to partial derivative Fock matrices
C  during solution of CPHF equations
C  ** CLOSED SHELL **
C
C    Jon Baker   April 2002
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
C  NAtoms  -  total number of symmetry-unique atoms
C  NAtOnce -  number of atoms on this pass
C  NatDo   -  number of atoms not yet converged, for wich quadrature
C             has to be performed (NatDo <= NatOnce)
C  XNuc    -  nuclear coordinates
C  XSym    -    ditto symmetry-unique atoms
C  IAN     -  nuclear charge (atomic number)
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
C              -1 - do not take weight derivatives
C               0 - take quadrature weight derivatives
C  DA      -  closed-shell density matrix (lower triangle)
C  DEN1    -  current solution to CPHF, lower triangle (closed-shell)
C  DM1     -  maximum element of DEN1 per column
C  NScr    -  size of general scratch array
C             (for dft < 3 only need MAX(2*NAtoms,NOcc), otherwise
C              need MAX(NBas*NBas,4*NOcc)
C  Z       -  general scratch array
C
C  on exit
C
C  FDA     -  derivative Fock matrices
C
C  lgrid   -  integer flag for generation/regeneration of DFT grid
C              0 - do not generate grid, already exists
C              1 - first entry so generate grid
C              2 - regenerate grid and REMOVE old grid files
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
C      J. Chem. Phys. 100 (1994) 7429
C
C  ----------------------------------------------------------------
C
C
      DIMENSION XNuc(3,NAtoms),XSym(3,NAtoms),IAN(NAtoms),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          BASDAT(13,*),INX(12,*),
     $          DA(NBas*(NBas+1)/2),DEN1(3,NatOnce,NBas*(NBas+1)/2),
     $          DM1(NBas),FDA(3,NatOnce,NBas*(NBas+1)/2)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
      CHARACTER*256 Filename,Filenam1,scrfile
      CHARACTER*3 ch3
C
      PARAMETER (Zero=0.0d0,Two=2.0d0)
      save mpoint,npptot
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
        tden1 = Zero
        tpot = Zero
        tnum = Zero
cc
        mpoint=0
        npptot=0
cc
        td1=zero
        tg1=zero
        tsw=zero
        tquad=zero
cc
        IEntry = 0
      ENDIF
c
      NBasS = NBas*NBatch
      NBasT = 3*NBasS
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
cc      if(IPrnt.gt.1) then
cc        write(6,*) '** DFT thresholds  thrsh:',thrsh,' thint:',thint
cc      end if
C
C  allocate scratch pointers
C
      iao = 1                 ! AO values over grid points
      inb = iao + NBasS       ! index array for "non-zero" AOs
      ibf = inb + NBasS       ! number of "non-zero" AOs per point
      id  = ibf + NBatch+1    ! density per grid point
      id1 = id  + NBatch      ! 1st-order density per grid point
      ipra  = id1 + 3*NatOnce ! potential per grid point
      iprara = ipra + NBatch  ! gradient of potential per grid point
      ind = iprara + NBatch   ! index of "non-zero" densities
      iscr = ind + NBatch     ! scratch storage (max AO value per point)
      IEnd = iscr + NBatch - 1
cc
      i1 = inb                  ! scratch pointers for grid
      i2 = i1 + NAtoms
      i3 = i2 + NAtoms
      i4 = i3 + NAtoms*NAtoms
cc
      If(dft.GT.3) Then
C
C  additional scratch space for non-local functionals
C
        iaox = 1 + IEnd         ! AO derivatives at grid point
        idx  = iaox + NBasT     ! density gradient at grid points
        id1x1= idx  + 3*NBatch  ! Functional gradient derivative
        id1x2= id1x1+ 3*NatOnce ! in situ first order gradient invariant
        id1x3= id1x2+ 3*NatOnce !
        idg1 = id1x3+ 3*NatOnce !
        isv  = idg1 + 3*NatOnce ! V coefficients in numerical quadrature
        isw1 = isv  + 3*NatOnce ! W coefficients in numerical quadrature
        isw2 = isw1 + 3*NatOnce ! Func. deriv. w.r.t. alpha grad.
        isw3 = isw2 + 3*NatOnce !
        ipga = isw3 + 3*NatOnce !
        ipgc = ipga + NBAtch    ! Func. deriv. w.r.t. alpha beta grad.
        iprarb = ipgc + NBAtch  ! Func. 2nd deriv. w.r.t. alpha and beta dens.
        ipraga = iprarb + NBAtch ! Func. 2nd deriv. w.r.t. alpha dens. and grad.
        ipragb = ipraga + NBAtch ! func. 2nd deriv. w.r.t alpha dens. and beta grad.
        ipragc = ipragb + NBAtch ! Func. 2nd deriv. w.r.t. alpha dens and alpha beta grad.
        ipgaga = ipragc + NBAtch ! Func. 2nd. deriv. w.r.t. alpha grad.
        ipgagb = ipgaga + NBAtch ! Func. 2nd. deriv. w.r.t. alpha and beta grad.
        ipgagc = ipgagb + NBAtch ! Func. 2nd. deriv. w.r.t. alpha and alpha beta grad.
        ipgcgc = ipgagc + NBAtch ! Func. 2nd. deriv. w.r.t. alpha beta grad.
        Iend = ipgcgc + NBAtch -1
      EndIf
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,8,'DFTCphfC')
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
c
      If(lgrid.GT.0) THEN
        CALL GRID(NAtoms, XSym,   IAN,    IPRNT,  NRad,
     $            NAng,   factor, IradQ,  Z(i1),  Z(i2),
     $            Z(i3),  DISTN,  RDIST,  AIJ,    ICntr,
     $            XXA,    WTA,    NPoint, XGRID,  WGHT)
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
      nrec = 0
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
C  form AO values over grid points
C  and sort "non-zero" AOs
C
      call secund(t1)
      CALL zeroit(Z(iao),NBasS)
      If(dft.LE.3) Then
        CALL AOVal(NP,     NShell, NBas,   XGRID(1,IP),
     $             XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $             ExpMIN, Z(iao))
        CALL AOSort(NP,     NBas,   thint,  Z(iao), Z(ibf),
     $              Z(inb), Z(iao), Z(iscr))
      Else
        CALL zeroit(Z(iaox),NBasT)
        CALL AOGrad(NP,     NShell, NBas,   XGRID(1,IP),
     $              XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $              ExpMIN, Z(iao), Z(iaox))
        CALL AOXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $               Z(ibf), Z(inb), Z(iao), Z(iaox),Z(iscr))
      EndIf
      call secund(t2)
      tao = tao + t2-t1
C
C  form density at current grid point
C
      call secund(t1)
      CALL GetDENSF(NP,     NBas,   Z(ibf), DA,     Z(iao),
     $              Z(inb), Z(id))
      call secund(t2)
      tden = tden + t2-t1
C
C  ..........................................................
C  on return from <GetDENSF>  Z(id) holds current densities
C  ..........................................................
C
C  form density derivatives if non-local
C
      If(dft.GT.3) Then
        call secund(t1)
        CALL GetDENDerivF(NP,     NBas,   thrsh,  z(ibf), DA,
     $                    Z(iao), z(iaox),z(inb), z(id),  z(idx))
        call secund(t2)
        tdenx = tdenx + t2-t1
      EndIf
C
C  ..........................................................
C  on return from <GetDENDerivF>  Z(idx) holds current density
C  derivatives
C  ..........................................................
C
C
C  calculate potential (and derivatives of potential)
C
      call secund(t1)
      CALL VXCPOTD(dft,      NP,       thrsh,    Z(id),    Z(idx),
     $             WGHT(IP), Z(ipra),  Z(iprara),Z(iprarb),Z(ipga),
     $             Z(ipgc),  Z(ipraga),Z(ipragb),Z(ipragc),Z(ipgaga),
     $             Z(ipgagb),Z(ipgagc),Z(ipgcgc),EXC,      EL,
     $             Jnk,      -1)
C
C  determine number of contributing points
C
      call potsorth(dft,      np,       thint,    z(ipra),  z(iprara),
     $              z(iprarb),z(ipga),  z(ipgc),  z(ipraga),z(ipragb),
     $              z(ipragc),z(ipgaga),z(ipgagb),z(ipgagc),z(ipgcgc),
     $              z(ind),   npp)
cc      write(6,*) ' Back from <POTSortH>  NPP is ',npp
cc      write(6,*) ' Sorted Pot and PotD are:'
cc      do i=0,npp-1
cc      write(6,*) i+1,'  ',z(iv+i),'  ',z(ivv+i)
cc      enddo
      npptot = npptot+npp
      call secund(t2)
      tpot = tpot + t2-t1
C
C  form "first-order density" at current grid point
C
cc      call secund(t1)
cc      CALL zeroit(Z(id1),3*NAtoms*NBatch)
cc      CALL DENCphf1(NP,     NAtoms, NBas,   thrsh,  Z(ibf),
cc     $              DEN1,   DM1,    Z(iao), Z(inb), Z(iscr),
cc     $              Z(ivv), Z(id1))
cc      call secund(t2)
cc      tden1 = tden1 + t2-t1
C
C  .....................................................................
C  on return from <DENCPHF1>  Z(id1) holds current 1st-order densities
C  .....................................................................
C
C  Accumulate DFT contribution into partial derivative Fock matrices
C
      call secund(t1)
      CALL makecphfc(dft,      npp,      nbas,     z(ibf),   natonce,
     $               natdo,    thrsh,    wght(ip), z(iprara),z(iprarb),
     $               Z(ipga),  z(ipgc),  z(ipraga),z(ipragb),z(ipragc),
     $               z(ipgaga),z(ipgagb),z(ipgagc),z(ipgcgc),DEN1,
     $               DM1,      z(idx),   z(iao),   z(iaox),  z(inb),
     $               z(iscr),  z(ind),   z(id1),   z(id1x1), z(id1x2),
     $               z(id1x3), z(idg1),  z(isv),   z(isw1),  z(isw2),
     $               z(isw3),  FDA,      td1,      tg1,      tquad,
     $               tsw)
      call secund(t2)
      tnum = tnum + t2-t1
C
C  ..................................................................
C    Derivatives of Weight Section     ** NOT ACTIVATED **
C  ..................................................................
C
      IP = IP+NP
      NLeft = NLeft-NP
      If(NLeft.GT.0) GO TO 50
C
C  end loop over grid points
C
C  ................................................................
C
c     If(IEntry.EQ.-1.AND.IPrnt.GT.0)then
c       write(6,2222) tgrid,tao,tden,tdenx,tpot,tnum,
c    $                td1,tg1,tsw,tquad,
c    $                mpoint,npptot
c     endif
C
      RETURN
 2222   Format(/,' Time taken to calculate grid:   ',f7.2,/,
     $           ' Time to calculate & sort AOs:   ',f7.2,/,
     $           ' Time to calculate density:      ',f7.2,/,
     $           ' Time to calculate den.derivs:   ',f7.2,/,
     $           ' Time to evaluate potential:     ',f7.2,/,
     $           ' Time numerical integration:     ',f7.2,/,
     $           ' of which:',/,
     $           '      for density derivatives:   ',f7.2,/,
     $           '      for gradient invariant:    ',f7.2,/,
     $           '      for coefficients s and w   ',f7.2,/,
     $           '      for quadrature Fock:       ',f7.2,/,
     $           ' Total number of grid points:    ',I8,/,
     $           ' Number of contributing points:  ',I8)
c
      END
