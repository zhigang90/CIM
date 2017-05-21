      subroutine dft_giao(idft,malk,lfxyz,na,ncf,ibas,iprint)
c
c  Adds exchange contribution to derivative fock matrices
c
c idft:   	dft type IN
c lfxyz:	pointer to derivative Fock matrices
c na:		number of atoms
c ncf:		number of contr. functions
c ibas:		pointer to basdat
c iprint:	print level
c
c  ...................................................................

      use memory

      implicit real*8 (a-h,o-z)
      data nrad/400/, nang/434/, IradQ/0/
      data NBatch/50/,lrad/0/,lang/0/, factor/1.25d0/
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
      common /symm/nsym,nsy(7)
      common /dftpnt/ ird,iaij,idst,iixx,iiwt,ipre,iexp,
     *                ixgd,iwt,iscr,icofm,nscr
c  ....................................................................
c initialize DFT
c restore coordinates, atomic charges and get all kinds of data
      ntri=ncf*(ncf+1)/2
      call mmark
      call getrval('ithr',thrsh)
      call getival('nslv',nslv)
      call getival('ictr',ictr)
      call getival('ncs',ncs)
      call getival('nmo',nmo)
      call getival('ngener',ngen)
      call getmem(3*na,ixnc)               ! nuclear coordinates
      call getmem(na,iqa)                  ! atomic charges/numbers
      call getmem(na,iuq)                  ! symmetry-unique atoms
      if (malk.eq.0) then
         call getmem(ncf*nmo,iadr1)        ! MO coeffs
      else
         call getmem(ncf*ncf,iadr1)        ! virtuals also
         call getival('imalk',imalk)
      endif
      iadr2=iadr1+ncf*nmo                  ! only used if malkin set
c
      call getnucdat(na,bl(ixnc),bl(iqa))
c
c prepare symmetry data (done even if no symmetry)
      nupr = 1                             ! not set if symmetry not used
      If(nsym.gt.0) call getival('SymNuPr1',nupr)
      call GetUNQ(na, nsym,   bl(ixnc),bl(nupr),NQ,
     *              bl(iuq),bl(iadr1))        ! iadr1 used as scratch
c if not parallel set up DFT pointers and arrays
      If(nslv.eq.0) Then
        call setup_dft(idft,   4,      1,  nrad,   nang,
     *                 NBatch, bl(iqa),bl(ixnc))
c I do not appreciate this way of memory allocation
        if(malk.ne.0) then
           nscr=nscr+NBatch
           nscr=nscr+NBatch*(ncf-nmo)
        endif
        call getmem(nscr,iscr)
      EndIf
c2006
c read in mo coefficients and transpose them
c       if(malk.eq.0)then
c          call getmem(ncf*nmo,ltmp)
c          call rea(bl(ltmp),ncf*nmo,4,'evec_rhf')
c      else           ! also the virtuals with malkin
c          call getmem(ncf*ncf,ltmp)
c          call rea(bl(ltmp),ncf*ncf,4,'evec_rhf')
c          call trspmo(bl(ltmp+ncf*nmo),ncf,bl(iadr2),ncf-nmo)
c       endif
c2006
        call getival('lvec',ltmp) ! evec location from CHSHIFT
c
        if(malk.eq.0)then
c          call getmem(ncf*nmo,ltmp)
c          call rea(bl(ltmp),ncf*nmo,4,'evec_rhf')
       else           ! also the virtuals with malkin
c          call getmem(ncf*ncf,ltmp)
c          call rea(bl(ltmp),ncf*ncf,4,'evec_rhf')
           call trspmo(bl(ltmp+ncf*nmo),ncf,bl(iadr2),ncf-nmo)
        endif
c2006
        call trspmo(bl(ltmp),ncf,bl(iadr1),nmo)
c2006   call retmem(1)
        call getmem(3*ntri,lfd1)  !exchange der. fock matrix
        call zeroit(bl(lfd1),3*ntri)
c
      CALL Para_Dftgiao(idft,   na,   bl(ixnc), bl(iqa), nsym,
     *                  ngen,   NSY,  bl(nupr), NQ,    bl(iuq),
     *                  iprint, lrad,   lang,    IradQ, factor,
     *                  NBatch,bl(idst),bl(iaij),bl(ird),bl(iixx),
     *                  bl(iiwt),bl(ixgd),bl(iwt),thrsh, ncf,
     *                  ncs, bl(ibas),bl(ictr),bl(ipre),bl(iexp),
     *                  nmo, bl(iadr1),bl(iadr2),malk,
     *                  nscr,  bl(iscr),bl(lfd1),bl(imalk))
c
c form vector product with nuclear coords and add to F1 matrices
c
      call vecnuc(na,ncs,ncf,ntri,bl(ictr),bl(ixnc),bl(lfd1),bl(lfxyz))
      call retmark
      end
c ====================================================================
      SUBROUTINE DFTNMRC(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                   NSym,   NGen,   ISYM,   NEqATM, IPRNT,
     $                   NRad,   NAng,   IradQ,  factor, NBatch,
     $                   DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                   XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                   BASDAT, INX,    BL,     ExpMIN, NOcc,
     $                   CMO,    VCMO,   Malk,   NScr,   Z,
     $                   XC,     XMalkin,IEntry)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates DFT contribution to the 3 derivative fock matrices on atom
C  ** CLOSED SHELL **
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
C  XC      -  derivative fock matrices, with exchange contribution added
C  XMalkin -  Malkin correction matrix
C
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
     $          CMO(NOcc,NBas),XC(3*NBas*(NBas+1)/2),
     $          VCMO(NBas-NOcc,NBas),XMalkin(NBas-NOcc,NOcc)

      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
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
cc
        IEntry = 0
      ENDIF
c
      lsemi = 0     ! needs to be checked later  JB
c
      NBasS = NBas*NBatch
      NBasT = 3*NBasS
      ntri = NBas*(NBas+1)/2
      ntri2 = 2*ntri
      ExpCUT = -LOG(thrsh)        ! exponent cutoff equivalent to thrsh
      thint = 0.1d0*SQRT(thrsh)
      if(IPrnt.gt.2) then
        write(6,*) '** DFT thresholds  thrsh:',thrsh,' thint:',thint
      end if
C
C  allocate scratch pointers
C
      iao = 1                         ! AO values over grid points
      imo = iao + NBasS               ! MO values over grid points
      inb = imo + NBatch*NOcc         ! index array for "non-zero" AOs
      ibf = inb + NBasS               ! number of "non-zero" AOs per point
      id  = ibf + NBatch+1            ! density per grid point
      iv  = id  + NBatch              ! potential per grid point
      ind = iv  + NBatch              ! index of "non-zero" densities
      IEnd = ind + NBatch - 1
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
        imox = iaox + NBasT           ! MO derivatives
        idx  = imox + 3*NOcc          ! density derivatives
        ivx  = idx  + 3*NBatch        ! potential derivatives
        IEnd = ivx  + 3*NBatch - 1
      EndIf
      if (Malk.ne.0) then
         ievec= 1 + IEnd                 ! energy contibution of each grid point
         ivmo = ievec + NBatch           ! virtual MO values over grid points
         IEnd = ivmo + (NBas-NOcc)*NBatch - 1
      else
         ievec = imo
      endif
C
C  make sure there is enough scratch memory
C
      CALL MemCHK(NScr,IEnd,7,'DFTNMRC')
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
        call zeroit(Z(iao),NBasS)
        IF(dft.LE.3) THEN
          CALL AOVal(NP,     NShell, NBas,   XGRID(1,IP),
     $               XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $               ExpMIN, Z(iao))
          CALL AOSort(NP,     NBas,   thint,  Z(iao), Z(ibf),
     $                Z(inb), Z(iao), Z(iv))
        ELSE
          CALL zeroit(Z(iaox),NBasT)
          CALL AOGrad(NP,     NShell, NBas,   XGRID(1,IP),
     $                XNuc,   BL,     BASDAT, INX,    ExpCUT,
     $                ExpMIN, Z(iao), Z(iaox))
          CALL AOXSort(NP,     NBas,   thint,  Z(iao), Z(iaox),
     $                 Z(ibf), Z(inb), Z(iao), Z(iaox),Z(iv))
        ENDIF
        call secund(t2)
        tao = tao + t2-t1
C
C  form density at current grid point
C
        if (Malk.ne.0) then          ! get virtual MO values too
           CALL GetDENS(NP,     NBas,   Z(ibf), NBas-NOcc,  VCMO,
     $                  Z(iao), Z(inb), Z(ivmo), Z(id))
        endif
        call secund(t1)
        CALL GetDENS(NP,     NBas,   Z(ibf), NOcc,   CMO,
     $               Z(iao), Z(inb), Z(imo), Z(id))
        CALL VScal(NP,Two,Z(id))             ! closed-shell density
        call secund(t2)
        tden = tden + t2-t1
C
C  ..........................................................
C  on return from <GetDENS>  Z(id) holds current densities
C  ..........................................................
C
C  prepare for numerical integration
C  sort "non-zero" densities
C
cc        CALL DENSort(NP,thrsh,Z(id),Z(ind),NPP)
C
C  form density derivatives if non-local
C
        If(dft.GT.3) Then
          call secund(t1)
          CALL GetDENDeriv(NP,     NBas,   thrsh,  Z(ibf), NOcc,
     $                     CMO,    Z(iaox),Z(inb), Z(id),  Z(imo),
     $                     Z(imox),Z(idx))
          call secund(t2)
          tdenx = tdenx + t2-t1
        EndIf
C
C  calculate potential (and potential derivatives)
C
        call secund(t1)
        CALL VXCPOT(dft,    NP,     thrsh,  IdWt,   Z(id),
     $             Z(idx), WGHT(IP),Z(iv),  Z(ivx), EXC,
     $              EL,  Z(ievec))
C
C  determine difference potential and number of contributing points
C
      CALL POTSort(dft,    lsemi,  NP,     thint,  nrec,
     $             Z(id),  Z(idx), Z(iv),  Z(ivx), Z(ind),
     $             NPP)
        call secund(t2)
        tpot = tpot + t2-t1
C
C  Accumulate DFT contribution into fock der. matrices
C
        call secund(t1)
        CALL MakeNMRXC(dft,    NPP,    NBas,   Z(ibf), thrsh,
     $                WGHT(IP),Z(iv),  Z(ivx), Z(iao), Z(iaox),
     $                 Z(inb), Z(ind), XGRID(1,IP), ntri, ntri2, XC)
C
C  Calculate the Malkin correction
C
        if(Malk.ne.0) then
         if(Malk.ne.dft) then
           CALL VXCPOT(Malk,    NP,    thrsh,  IdWt,   Z(id),
     $                Z(idx), WGHT(IP),Z(iv),  Z(ivx), EXC,
     $                 EL,  Z(ievec))
c no need to potsort again as the functional is most likely hfs
         endif
         call makemalkin(NP,    NPP,    NBas,   NOcc,   thrsh,
     *                   WGHT(IP), Z(iv),Z(imo),Z(ivmo),Z(ievec),
     *                   Z(id),Z(ind),XMalkin)
        endif
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
     $   write(6,2222) tgrid,tao,tden,tdenx,tpot,tnum
C
      RETURN
c
 2222   Format(/,' Time taken to calculate grid: ',f7.2,/,
     $           ' Time to calculate & sort AOs: ',f7.2,/,
     $           ' Time to calculate density:    ',f7.2,/,
     $           ' Time to calculate den.derivs: ',f7.2,/,
     $           ' Time to evaluate potential:   ',f7.2,/,
     $           ' Time numerical integration:   ',f7.2,/)
c
      END
c ====================================================================
      SUBROUTINE MakeNMRXC(dft,    NPP,    NBas,   nbf,    thrsh,
     $                     WGHT,   Pot,    PotX,   VAO,    VAOX,
     $                     INB,    INDX,   XGRID,  ntri, ntri2, XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid points into Exchange-Correlation Fock matrix derivatives
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All methods include Slater exchange)
C              1 - 3 - local correlation
C             >3 - nonlocal - need AO derivatives
C  NPP     -  number of contributing (non-zero) grid points this batch
C  NBas    -  total number of basis functions
C  nbf     -  number of "non-zero" basis functions
C  thrsh   -  threshold for neglect of contribution
C  WGHT    -  grid quadrature weights
C  Pot     -  potential
C  PotX    -  potential derivatives
C  VAO     -  basis function values at grid point
C  VAOX    -  basis function derivatives at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  INDX    -  index into contributing columns of VAO
C  XGRID   -  gridpoint coordinates
C  ntri,ntri2 - shift of y and z derivative matrices in XC
C  XC      -  exchange-correlation Fock matrix derivatives
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),Pot(NPP),
     $          PotX(3,NPP),XGRID(3,*),XC(*)
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
          S= Val*VAO(I)
          XC(IJ) = XC(IJ) + S*XGRID(1,IP)
          XC(IJ+ntri) = XC(IJ+ntri) + S*XGRID(2,IP)
          XC(IJ+ntri2) = XC(IJ+ntri2) + S*XGRID(3,IP)
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
        TX = ValJ*VAx
        TY = ValJ*VAy
        TZ = ValJ*VAz
        ValJA = ValJ*VXC +
     $    (VAOX(1,J)*VAx + VAOX(2,J)*VAy + VAOX(3,J)*VAz)
        IF(Abs(ValJA).GT.thrsh) THEN
          JT = (JJ*(JJ-1))/2
          DO 50 I=nbf(IPP)+1,J
          II = INB(I)
          IJ = JT + II
          S =  VAO(I)*ValJA +
     $     (VAOX(1,I)*VAx + VAOX(2,I)*VAy + VAOX(3,I)*VAz)*ValJ
          XC(IJ) = XC(IJ) + S*XGRID(1,IP) + TX*VAO(I)
          XC(IJ+ntri) = XC(IJ+ntri) + S*XGRID(2,IP) + TY*VAO(I)
          XC(IJ+ntri2) = XC(IJ+ntri2) + S*XGRID(3,IP) + TZ*VAO(I)
 50       CONTINUE
        ENDIF
 60     CONTINUE
 70     CONTINUE
      ENDIF
C
      RETURN
      END
c ==================================================================
      subroutine vecnuc(na,ncs,ncf,ntri,inx,xnuc,fder,fock)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER p,q
c
c forms dft exchange by taking vector product of the precalculated matrix
c (fder) and the nuclear coordinates, result added to (fock)
c
C  na   -  number of atoms
C  ncs  -  total number of shells
C  ncf  -  total number of basis functions
C  ntri - fock matrix dimension
C  inx     -  more basis set data
C  xnuc    -  nuclear coordinates
C  fder    -  exchange correlation deriv. fock matrix data
C  fock    -  GIAO derivative fock matrix (input - Coulomb, output ready)
      DIMENSION xnuc(3,na), inx(12,*), fder(3*ntri),
     *          fock(3*ntri),iatn(ncf)
      ntri2=2*ntri
      do 550 ics=1,ncs
        ia=inx(11,ics)+1
        if(ics.lt.ncs) then
          ie=inx(11,ics+1)
        else
          ie=ncf
        endif
        do 540 q=ia,ie
          iatn(q)=inx(2,ics)
 540    continue
 550  continue
      k=0
      do 700 q=1,ncf
        iaq=iatn(q)
        do 600 p=1,q
          iap=iatn(p)
          k=k+1
c  these are the center coordinates: Rp-Rq
          xx=xnuc(1,iap)-xnuc(1,iaq)
          yy=xnuc(2,iap)-xnuc(2,iaq)
          zz=xnuc(3,iap)-xnuc(3,iaq)
          vx=fder(k)
          vy=fder(k+ntri)
          vz=fder(k+ntri2)
          fock(k)=fock(k)-(yy*vz-zz*vy)
          fock(k+ntri)=fock(k+ntri)-(zz*vx-xx*vz)
          fock(k+ntri2)=fock(k+ntri2)-(xx*vy-yy*vx)
 600    continue
 700  continue
      end
c ======================================================================
      subroutine makemalkin(NP,     NPP,    NBas,   NOcc,   thrsh,
     $                      WGHT,   Pot,    VMO,    VMOV,   EVec,
     $                      DEN,    INDX,   XMalkin)
      IMPLICIT REAL*8(A-H,O-Z)
C
C Numerical integration to calculate the malkin correction matrices
C using a batch of gridpoints
C
C Arguments:
C
C  NP      -  number of points
C  NPP     -  number of contributing (non-zero) grid points this batch
C  NBas    -  total number of basis functions
C  NOcc    -  number of occupied orbitals
C  thrsh   -  threshold for neglect of contribution
C  WGHT    -  grid quadrature weights
C  Pot     -  potential
C  VMO     -  MO values at grid point
C  VMOV    -  virtual MO values at grid point
C  EVEC    -  vector of exchange-correlation energy contributions
C             per grid point (premultiplied by the weights)
C  DEN     -  density values per grid point
C  INB     -  indexing array for non-zero entries to VAO
C  INDX    -  index into contributing columns of VAO
C  XMalkin -  Malkin corr. matrix
C
C
      DIMENSION WGHT(*),Pot(NP),VMO(NOcc,NP),VMOV(NBas-NOcc,NP),
     *          EVec(NP),DEN(NP),XMalkin(NBas-NOcc,NOcc)
      integer dft,INDX(NP)
      PARAMETER (felha=0.492372511d0) !0.5 (3/pi)**(1/3)
c constant is different because the density is closed shell already
C
      DO 30 IP=1,NPP
c       if (DEN(ipp).gt.thrsh)then      ! taken care of with sorting
        IPP=INDX(IP)
        corr=(EVEC(IPP)-Pot(IPP))*WGHT(IPP)*2.0d0/DEN(IPP)
c       corr=WGHT(IPP)*felha*den(iPp)**(-0.6666666666)
        DO 20 J=1,NOcc
          denfrac=corr*VMO(J,IPP)**2
          DO 10 I=1,NBas-NOcc
            XMalkin(I,J)=XMalkin(I,J)+denfrac*VMOV(I,IPP)**2
 10       CONTINUE
 20     CONTINUE
c       endif
 30   CONTINUE
      RETURN
      END
