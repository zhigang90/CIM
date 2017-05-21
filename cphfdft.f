c ==================================================================
c  DFT HESSIAN CPHF MODULE            JB   April 2002
c ==================================================================
c
      SUBROUTINE CPHFDFT(NAtoms, NatOnce,natend, npos,   NQ,
     $                   IUNQ,   DEN1,   Den1B,  FDA,    FDB,
     $                   thrsh,  lsemi,  NMem,   Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** DFT HESSIAN PROGRAM **
C
C  This Program determines the DFT exchange-correlation
C  contribution to the partial derivative Fock matrices
C  combining it with the Coulomb contribution previously
C  calculated (prior to solving the CPHF equations)
C  ............................................................
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  NatOnce -  number of atoms on this pass
C  natend  -  array indicating for which atoms the CPHF is already
C             converged
C  npos    -  starting position of current first atom
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of symmetry-unique atoms
C  DEN1    -  current solution to CPHF equations (1st-order density)
C             (alpha density for uhf case)
C  DEN1B   -  beta density for uhf
C  FDA     -  partial derivative Fock matrices (closed-shell or alpha)
C  FDB     -  beta partial derivative Fock matrices
C  thrsh   -  integral neglect threshold
C  lsemi   -  flag for reuse of DFT grid
C             (should be reset to 0 for each new Hessian calculation)
C  NMem    -  available scratch storage
C  Z       -  scratch storage
C
C
      DIMENSION FDA(3,natonce,*),DEN1(3,natonce,*)
      DIMENSION FDB(3,natonce,*),DEN1B(3,natonce,*)
      DIMENSION IUNQ(NQ),Z(NMem)
c  ...................................................................
c     common /big/ bl(1000)            ! old texas depository
c     common /intbl/maxsh,ixx(100)     ! old texas integer depository
c ....................................................................
      Integer dft
      Character cdum*20
      Logical rhf
C
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
c  ...................................................................
      Data nrad/400/, nang/434/, NBatch/50/, IdWt/0/
c  ...................................................................
c
      Character*256 jobname
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    dft functional type
C    multiplicity
C    number of alpha/closed-shell orbitals
C    number of beta orbitals
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call RdDFT(IUnit,dft)
      If(dft.gt.0) Then
        call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
        call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
        call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      EndIf
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  do we have a dft contribution?
C
      If(dft.EQ.0) RETURN
c
      rhf = (NBeta.EQ.0.AND.IMult.EQ.1)
C
C --------------------------------------------------------------------
C  Get all relevant pointers from Texas depository
C
      call getival('ncf',NBas)            ! number of basis functions
      call getival('ncs',NShell)          ! number of shells
      call getival('ibas',ibas)           ! location of BASDAT array
      call getival('nsh',NPrim)           ! number of primitive shells
      call getival('ictr',ictr)           ! location of INX array
      call getival('inuc',inuc)           ! nuclear coordinates etc...
      call getival('ldensi',ida)          ! closed-shell/alpha density matrix
      if(.not.rhf)call getival('ldensB',idb) ! beta density matrix
      call getival('nslv',nslv)           ! number of slaves (parallel)
C --------------------------------------------------------------------
C
C  number of atoms to do on this pass is natonce
C  Now get the memory
C
      MaxGrid = NRad*NAng
c
      IMem = 7*NAtoms + NQ + NBas + NatOnce
      NScr = NAtoms
      If(nslv.eq.0) Then
       IMem = IMem + NAtoms + 2*NAtoms**2 + 4*1130 + 4*MaxGrid +
     $               NPrim + NBas
c -- scratch memory for closed-shell local DFT
       NScr = 2*NBas*NBatch + 6*NBatch + 1 + 3*NatOnce
       if(.not.rhf) NScr = NScr + 4*NBatch + 9*NatOnce !local uhf
c -- additional scratch memory for non-local DFT
       If(dft.GT.3)then
         NScr = NScr + 3*NBas*NBatch + 13*NBatch + 24*NatOnce
         if(.not.rhf) NScr = NScr + 8*NBatch + 24*NatOnce!non-local uhf
       endif
      EndIf
c
      IMem = IMem + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,7,'CPHFDFT')
      CALL ZeroIT(Z,IMem)      ! clear the memory
C
C  Allocate memory pointers
C
      IXC  = iptr                       !  nuclear coordinates
      IXS  = IXC + 3*NAtoms             !   ditto symmetry-unique atoms
      IAN  = IXS + 3*NAtoms             !  nuclear charge (atomic numbers)
      IUQ  = IAN + NAtoms               !  duplicate list of symmetry unique atoms
      IDM  = IUQ + NQ                   !  largest 1st-order density element per column
      incex = IDM + NBas                !  index array for contracting density/fock matrices
      IEnd = incex + NatOnce
c
c
c -- the following arrays are NOT needed in parallel mode
      IF(nslv.EQ.0) THEN
       IDST = IEnd                      !  distance to nearest neighbour
       IAIJ = IDST + NAtoms             !  Becke grid parameters
       IRST = IAIJ + NAtoms*NAtoms      !  inverse interatomic distances
       IXXA = IRST + NAtoms*NAtoms      !  angular quadrature grid
       IWTA = IXXA + 3*1130             !  angular quadrature weights
       IGRD = IWTA + 1130               !  X,Y,Z grid points
       IWGT = IGRD + MaxGrid*3          !  grid quadrature weights
       IBL  = IWGT + MaxGrid            !  precomputed normalization factors
       IExp = IBL  + NPrim              !  smallest exponent per basis function
       IEnd = IExp + NBas
      ENDIF
C
C  scratch memory
C
      IScr = IEnd
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,7,'CPHFDFT')
c -- assign memory for high water mark
      call getmem(IEnd,lastx)
c
c -- get coordinates and atomic numbers
c    from Texas format to decent format
      call getnucdat(Natoms,Z(IXC),Z(IAN))
c
c -- duplicate list of symmetry-unique atoms (original is modified)
      CALL ICpyVEC(NQ,IUNQ,Z(IUQ))
c
c -- form maximum 1st-order density element per column
c     if(rhf)then
c       CALL FormMaxDen1(3*NatOnce,NBas,DEN1,Z(IDM))
c     else
c       CALL FormMaxDen1U(3*NatOnce,NBas,DEN1,DEN1B,Z(IDM))
c     endif
C  ----------------------------------------------------------------------
C
      CALL CPHFMAIN(dft,     NAtoms,  NatOnce, npos,    Z(IXC),
     $              Z(IXS),  Z(IAN),  NQ,      Z(IUQ),  NRad,
     $              NAng,    NBatch,  Z(IDST), Z(IAIJ), Z(IRST),
     $              Z(IXXA), Z(IWTA), Z(IGRD), Z(IWGT), thrsh,
     $              NBas,    NShell,  NPrim,  bl(ibas),bl(ictr),
     $              Z(IBL),  Z(IExp), NAlpha,  NBeta,   rhf,
     $              bl(ida), bl(idb), DEN1,    DEN1B,   Z(IDM),
     $              natend, z(incex),FDA,     FDB,     NScr,
     $              Z(IScr), lsemi)
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      RETURN
      END
c ========================================================================
c
      SUBROUTINE CPHFMAIN(dft,    NAtoms, NatOnce,npos,   XNuc,
     $                    XSym,   IAN,    NQ,     IUNQ,   NRad,
     $                    NAng,   NBatch, DISTN,  AIJ,    RDIST,
     $                    XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                    NBas,   NShell, NPrim,  BASDAT, INX,
     $                    BL,     ExpMIN, NAlpha, NBeta,  rhf,
     $                    DA,     DB,     DEN1,   DEN1B,  DM1,
     $                    natend, ncex,   FDA,    FDB,    NScr,
     $                    Z,      lsemi)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for DFT contributions to Hessian
C  and derivative Fock matrices
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),XGRID(3,*),WGHT(*),
     $          BASDAT(13,*),INX(12,*),BL(*),ExpMIN(NShell),
     $          DA(*),FDA(3,NatOnce,(NBas*(NBas+1))/2),
     $          DB(*),FDB(3,NatOnce,(NBas*(NBas+1))/2),
     $          DEN1(3,NatOnce,(NBas*(NBas+1))/2),DM1(NBas),
     $          DEN1B(3,NatOnce,(NBas*(NBas+1))/2)
      DIMENSION IUNQ(NQ),XSym(3,NAtoms),natend(NatOnce),ncex(NatOnce)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      Character cdum*20
      INTEGER dft
      Logical rhf,compact
      DIMENSION Z(NScr)
c  ......................................................
      data lrad/0/, lang/0/, IradQ/0/, lpos/0/
c  ......................................................
C
      PARAMETER (IUnit=1)
      SAVE lpos
c
      Character*256 jobname
      Common /job/jobname,lenJ
C
C
C  initialization
C
      call secund(t1)
      call elapsec(elacphf1)
C     WRITE(6,1000)
ccccc
      if(lpos.ne.npos) then
        lsemi = 0
        lpos = npos
      endif
ccccc
c
c    if for some of the atoms the cphf has already converged,
c    we can save some time by doing the numerical quadrature
c    only for the not yet converged atoms.
c    this is done by compacting the densities and fock matrices,
c    so that the non converged components are contiguously stored.
c
c    this process uses the following variables:
c
c    natdo      number of atoms that have to be computed (<= natonce)
c    compact    logical, whether we need to compact the matrices
c    ncex       array for pointers needed in compacting/expanding
c
c   we compute natdo
c
      natdo=0
      do ia=1,natonce
        if(natend(ia).eq.0)natdo=natdo+1
      enddo
c
c   we check whether we need to compact and
c   we fill the ncex array
c
      compact=.false.
      if(natdo.ne.natonce)then
        ist=natdo+1
        do 10 ia=1,natdo
          if(natend(ia).ne.0)then
            do ib=ist,natonce
              if(natend(ib).eq.0)then
                compact=.true.
                ncex(ia)=ib
                ist=ib+1
                goto 10
              endif
            enddo
          else
            ncex(ia)=0
          endif
 10     continue
      endif
c
c     if needed, we compact the matrices
c
      if(compact)then
        call swap_df1(DEN1,ncex,NatOnce,NBas,natdo)
        call swap_df1(FDA,ncex,NatOnce,NBas,natdo)
        if(.not.rhf)then
          call swap_df1(DEN1B,ncex,NatOnce,NBas,natdo)
          call swap_df1(FDB,ncex,NatOnce,NBas,natdo)
        endif
      endif
c
c -- form maximum 1st-order density element per column
c    moved this here from <cphfdft> (MM, October 2003)
c
      if(rhf)then
        CALL FormMaxDen1(3*NatOnce,natdo,NBas,DEN1,DM1)
      else
        CALL FormMaxDen1U(3*NatOnce,natdo,NBas,DEN1,DEN1B,DM1)
      endif
C
C ...............................................................
C  read from the <control> file
C    grid factor
C    print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$factor',2,idum,factor,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C ................................................................
C  get data from old texas depository
C
      call getival('nslv',nslv)           ! number of slaves (parallel)
C
C  If there is symmetry, we are only going to use the symmetry-unique
C  atoms when generating the grid and solving the CPHF equations.
C  Put these atoms in XSym
C
      If(NQ.LT.NAtoms) Then
        DO I=1,NQ
        IAtm = IUNQ(I)
        XSym(1,I) = XNuc(1,IAtm)
        XSym(2,I) = XNuc(2,IAtm)
        XSym(3,I) = XNuc(3,IAtm)
        IAN(I) = IAN(IAtm)
        IUNQ(I) = I
        EndDO
      Else
        CALL CpyVEC(3*NAtoms,XNuc,XSym)
      EndIf
C
C  Now calculate the various DFT contributions
C
      IF(nslv.EQ.0) THEN
cc
c  SINGLE PROCESSOR MODE
c  ---------------------
C
C  calculate inverse atomic distances, Becke aij parameters
C  and angular grid and weights prior to full grid construction
C
      CALL PreGRID(NQ,     XSym,   IAN,    XGRID,  WGHT,
     $             DISTN,  RDIST,  AIJ,    XXA,    WTA)
C
C  get array of smallest exponent per shell
C
      CALL GetEXP(NPrim,NShell,BASDAT,INX,ExpMIN)
C
C  precompute shell normalization factors
C
      CALL AOInit(NShell,BASDAT,INX,BL)
c
        lsemi = lsemi+1
        lgrid = 0
        If(lsemi.EQ.1) lgrid = 1
c
        DO 100 IAtom=1,NQ
        If(IAtom.EQ.1) Then
          IEntry = 1
        Else If(IAtom.EQ.NQ) Then
          IEntry = -1
        EndIf
        ICntr = IUNQ(IAtom)
        If(rhf) Then
          CALL DFTCphfC(dft,    ICntr,  NQ,     NatOnce,NatDo,
     $                  XNuc,   XSym,   IAN,    IPRNT,  lrad,
     $                  lang,   IradQ,  factor, NBatch, DISTN,
     $                  AIJ,    RDIST,  XXA,    WTA,    XGRID,
     $                  WGHT,   thrsh,  NBas,   NShell, BASDAT,
     $                  INX,    BL,     ExpMIN, IdWt,   DA,
     $                  DEN1,   DM1,    NScr,   Z,      FDA,
     $                  lgrid,  IEntry)
        Else
          CALL DFTCphfU(dft,    ICntr,  NQ,     NatOnce,NatDo,
     $                  XNuc,   XSym,   IAN,    IPRNT,  lrad,
     $                  lang,   IradQ,  factor, NBatch, DISTN,
     $                  AIJ,    RDIST,  XXA,    WTA,    XGRID,
     $                  WGHT,   thrsh,  NBas,   NShell, BASDAT,
     $                  INX,    BL,     ExpMIN, IdWt,   DA,
     $                  DB,     DEN1,   DEN1B,  DM1,    NScr,
     $                  Z,      FDA,    FDB,    lgrid,  IEntry)
        EndIf
 100    CONTINUE
cc
      ELSE
cc
c  PARALLEL MODE
c  -------------
c
        CALL Para_DFTCphf(dft,    NQ,     NatOnce,NatDo,  XNuc,
     $                    XSym,   IAN,    NQ,     IUNQ,   lrad,
     $                    lang,   IradQ,  factor, NBatch, thrsh,
     $                    NBas,   NShell, BASDAT, INX,    IdWt,
     $                    rhf,    DA,     DB,     DEN1,   DEN1B,
     $                    DM1,    FDA,    FDB,    lgrid)
cc
      ENDIF
c
      lgrid = 0
c
c     we expand back densities and fock matrices, if needed
c
      if(compact)then
        call swap_df1(DEN1,ncex,NatOnce,NBas,natdo)
        call swap_df1(FDA,ncex,NatOnce,NBas,natdo)
        if(.not.rhf)then
          call swap_df1(DEN1B,ncex,NatOnce,NBas,natdo)
          call swap_df1(FDB,ncex,NatOnce,NBas,natdo)
        endif
      endif
C
C ------------------------------------------------------------------
C
C
      If(IPRNT.GT.4) Then
        call getival('iout',iout)
        write(iout,*) ' 1st-order density matrices are:'
        do i=1,natonce
        write(6,*) ' Atom:',i,' X-derivative density matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (den1(1,i,kk+l),l=1,k)
        enddo
        if(.not.rhf)then
        write(6,*) ' Atom:',i,' beta X-derivative density matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (den1B(1,i,kk+l),l=1,k)
        enddo
        endif
c
        write(6,*) ' Atom:',i,' Y-derivative density matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (den1(2,i,kk+l),l=1,k)
        enddo
        if(.not.rhf)then
        write(6,*) ' Atom:',i,' beta Y-derivative density matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (den1B(2,i,kk+l),l=1,k)
        enddo
        endif
c
        write(6,*) ' Atom:',i,' Z-derivative density matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (den1(3,i,kk+l),l=1,k)
        enddo
        if(.not.rhf)then
        write(6,*) ' Atom:',i,' beta Z-derivative density matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (den1B(3,i,kk+l),l=1,k)
        enddo
        endif
        enddo
c
 1112   format(1x,6(2x,F10.6))
        write(iout,*) ' Partial Fock matrices are:'
        do i=1,natonce
        write(6,*) ' Atom:',i,' X-derivative Fock matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (fda(1,i,kk+l),l=1,k)
        enddo
        if(.not.rhf)then
        write(6,*) ' Atom:',i,' beta X-derivative Fock matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (fdb(1,i,kk+l),l=1,k)
        enddo
        endif
c
        write(6,*) ' Atom:',i,' Y-derivative Fock matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (fda(2,i,kk+l),l=1,k)
        enddo
        if(.not.rhf)then
        write(6,*) ' Atom:',i,' beta Y-derivative Fock matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (fdb(2,i,kk+l),l=1,k)
        enddo
        endif
c
        write(6,*) ' Atom:',i,' Z-derivative Fock matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (fda(3,i,kk+l),l=1,k)
        enddo
        if(.not.rhf)then
        write(6,*) ' Atom:',i,' beta Z-derivative Fock matrix is:'
        do k=1,nbas
        kk = (k*(k-1))/2
        write(6,1112) (fdb(3,i,kk+l),l=1,k)
        enddo
        endif
        enddo
      EndIf
C  ----------------------------------------------------------------
C
      call secund(t2)
      call elapsec(elacphf2)
      t1 = (t2-t1)/60.0d0
      elaps = (elacphf2-elacphf1)/60.0d0
      If(IPRNT.GT.2) WRITE(6,1100) t1,elaps
c
      RETURN
c
 1000 FORMAT(/,' Calculating DFT Contribution to CPHF Fock matrices')
 1100 FORMAT('Master CPU time for XC part of CPHF Fock matrices  ='
     *,f8.2,' Elapsed = ',f8.2,' min')
c
      END
c ========================================================================
      SUBROUTINE FormMaxDen1(NAt3,natdo,NBas,DEN1,DM1)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum 1st-order density matrix element per column.
C  Used for density threshold in DFT code
C  (NOTE: Used for density AND difference-density matrices)
C
C  ARGUMENTS
C
C  NAt3    -  three*number of atoms
c  natdo   -  atoms that are to be used
C  NBas    -  number of basis functions
C  DEN1    -  1st-order density matrix (lower triangle)
C  DM1     -  on exit contains maximum 1st-order density element per column
C
C
      DIMENSION DEN1(NAt3,NBas*(NBas+1)/2),DM1(NBas)
C
      natd3=3*natdo
      DO 20 I=1,NBas
      DMx = 0.0d0
      II = (I*(I-1))/2
      DO 10 J=1,I
      IJ = II+J
      DO K=1,NAtd3
      If(Abs(DEN1(K,IJ)).GT.DMx) DMx = Abs(DEN1(K,IJ))
      EndDO
 10   CONTINUE
      DO 11 J=I+1,NBas
      IJ = (J*(J-1))/2 + I
      DO K=1,NAtd3
      If(Abs(DEN1(K,IJ)).GT.DMx) DMx = Abs(DEN1(K,IJ))
      EndDO
 11   CONTINUE
      DM1(I) = DMx
 20   CONTINUE
C
      RETURN
      END
c ========================================================================
      SUBROUTINE FormMaxDen1U(NAt3,natdo,NBas,DEN1A,DEN1B,DM1)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum 1st-order density matrix element per column.
C  Used for density threshold in DFT uhf code
C  (NOTE: Used for density AND difference-density matrices)
C
C  ARGUMENTS
C
C  NAt3    -  three*number of atoms
c  natdo   -  atoms that are to be used
C  NBas    -  number of basis functions
C  DEN1A   -  alpha 1st-order density matrix (lower triangle)
C  DEN1B   -  beta 1st-order density matrix (lower triangle)
C  DM1     -  on exit contains maximum 1st-order density element per column
C
C
      DIMENSION DEN1A(NAt3,NBas*(NBas+1)/2),DM1(NBas)
      DIMENSION DEN1B(NAt3,NBas*(NBas+1)/2)
C
      natd3=3*natdo
      DO 20 I=1,NBas
      DMx = 0.0d0
      II = (I*(I-1))/2
      DO 10 J=1,I
      IJ = II+J
      DO K=1,NAtd3
      DMx=max(DMx,Abs(DEN1A(K,IJ)),Abs(DEN1B(K,IJ)))
      EndDO
 10   CONTINUE
      DO 11 J=I+1,NBas
      IJ = (J*(J-1))/2 + I
      DO K=1,NAtd3
      DMx=max(Dmx,Abs(DEN1A(K,IJ)),Abs(DEN1B(K,IJ)))
      EndDO
 11   CONTINUE
      DM1(I) = DMx
 20   CONTINUE
C
      RETURN
      END
c ========================================================================
      subroutine swap_df1(df1,ncex,NatOnce,NBas,natdo)
      implicit real*8(a-h,o-z)
c
      dimension df1(3,NatOnce,Nbas*(Nbas+1)/2)
      integer ncex(natdo)
      do ia=1,natdo
        if(ncex(ia).ne.0)then
          ib=ncex(ia)
          do ij=1,Nbas*(Nbas+1)/2
            ax=df1(1,ib,ij)
            ay=df1(2,ib,ij)
            az=df1(3,ib,ij)
            df1(1,ib,ij)=df1(1,ia,ij)
            df1(2,ib,ij)=df1(2,ia,ij)
            df1(3,ib,ij)=df1(3,ia,ij)
            df1(1,ia,ij)=ax
            df1(2,ia,ij)=ay
            df1(3,ia,ij)=az
          enddo
        endif
      enddo
c
      return
      end
