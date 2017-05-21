C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C THE REAL PARA_CPHF ROUTINE WAS MOVED!!!!!!
c=====================================================================
c
      subroutine pcphf_xyz(idft,ax,nblocks,bl,inx,ntri,thres,mgo,
     *                     dens,fock,labels,lsemi)
c
c=====================================================================
c This routine calls two-el. int. block by block and constructes
c losed-shell Fock matrix using 1st-order density.
c---------------------------------------------------------------------
      use memory, block => bl
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),dens(*),fock(*)
      dimension labels(*)
      parameter (Zero=0.0d0)
      data ncall /1/
      save ncall
cccc  data lsemi /0/
cccc  save lsemi
c----------------------------------------------------------------
c
c if this is the first cycle then store integrals
c
      if(scftype.ne.'full-direct' .and. ncall.eq.1) then
         call getrval('thr2',thres1) ! final integral threshold
         call pstore(nblocks,bl,inx,thres1,labels)
         ncall=ncall+1
      endif
c----------------------------------------------------------------
      call getival('ncf',ncf)
cxxx  call getival('na',natoms)
      call getival('realat',natoms)
ccc   call getival('gran',igran)
c----------------------------------------------------------------
c
      call zeroit(fock,ntri*3)
      call para_cphf_xyz(idft,ax,nblocks,bl,inx,ntri,thres,dens,fock,
     *                   labels)
c
c----------------------------------------------------------------
      thresax=thres
c
c re-scale fock matrices :
c
      call rescale_fock(fock,ntri,bl(lind),thresax)
c
c----------------------------------------------------------------
c symmetry data :
c
      call getival('nsym',nsym)
      call getival('ngener',ngener)
      call getival('ngener',ngener)
      call getival('nsyo',nsyo)
      call getival('SymFunPr',ifp)
      nupr = 1      ! not set if symmetry not used
      If(nsym.gt.0) call getival('SymNuPr1',nupr)
c----------------------------------------------------------------
c  DFT CONTRIBUTIONS
c
c  Calculate DFT XC contributions to partial derivative Fock
c  matrices for CPHF equations
c
      IF(idft.GT.0) THEN
c
c get number of symmetry-unique atoms and pointer to the list
c (set up in polariz)
c
        call getival('nofuniq',nq)
        call getival('iunqato',iunq)
c
        call getival('lcore',lcore)
        call getmem(0,lastx)
        call retmem(1)
c
c set up irrelevent stuff :
c
         npos=0
c
         call cphfdft_ef(natoms,npos,
     *               bl(nupr),
     *               nsym,ngener,NQ, bl(nsyo),bl(iunq), 
     *               dens,   dens,  fock,   fock,
     $               thres,  lsemi,lcore-lastx,bl(lastx))
c----------------------------------------------------------------
      ENDIF
c----------------------------------------------------------------
c Symetrize the Ist-order Fock matrices (G(D1,g0)) for exac symmetry :
c
      if(nsym.gt.0) then
         call fdersymm_gen(ngener,ncf,fock,bl(ifp),bl(nsyo))
      endif
c---------------------------------------------------------------------
c
      end
c=====================================================================
c     symmetrization of the first-order Fock matrices
c            in a general CPHF case except NMR
c===================================================================
      subroutine fdersymm_gen(ngener,ncf,fock,ifp,nsymop)
c---------------------------------------------------------------------
c currently uesed in polarizability CPHF 
c---------------------------------------------------------------------
c _gen - stands for general (not for NMR)
c
c    symmetrize F(D0,g1) (with the integral derivatives)
c                 or
c    symmetrize F(D1,g0) (with the ordinary zeroth-order integrals)
c               if D1 does not break symmetry
c
c input : ngener= number of group generators
c         ncf   = number of basis functions
c         fock  = 3 Fock matrices (X,Y,Z perturbed)
c         ifp   = symmetry images of a given function for each sym.op.
c         nsymop= symmetry operations
c
c output: fock  = symmetrized Fock matrices
c
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension fock(*)
      dimension nsymop(7),ifp(7,ncf)
c local stuff:
      dimension mirror(3,7),negx(7),negy(7),negz(7)
      data mirror/-1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
      data half,one /0.5d0,1.d0/
c-----------------------------------------
      if (ncf.eq.1) return
c-----------------------------------------
      ntri=ncf*(ncf+1)/2
      lrix=0
      lriy=lrix+ntri
      lriz=lriy+ntri
c-----------------------------------------
      do 10 ns=1,ngener
         nop=nsymop(ns)
         negx(ns)=mirror(1,nop)
         negy(ns)=mirror(2,nop)
         negz(ns)=mirror(3,nop)
   10 continue
c-----------------------------------------
      do 350 ns=1,ngener
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if(i1.lt.0) then
               i1=-i1
               fct=-fct
            endif
            j1=ifp(ns,j)
            if (j1.lt.0) then
              j1=-j1
              fct=-fct
            endif
            ij1=i1*(i1-1)/2 +j1
            if (j1.gt.i1) then
              ij1=j1*(j1-1)/2 +i1
            endif
            ffx=fock(lrix+ij)+fct*fock(lrix+ij1)*negx(ns)
            ffy=fock(lriy+ij)+fct*fock(lriy+ij1)*negy(ns)
            ffz=fock(lriz+ij)+fct*fock(lriz+ij1)*negz(ns)
            if (ij.gt.ij1) then
               ffx=ffx*half
               ffy=ffy*half
               ffz=ffz*half
            endif
            fock(lrix+ij)=ffx
            fock(lrix+ij1)=fct*ffx*negx(ns)
            fock(lriy+ij)=ffy
            fock(lriy+ij1)=fct*ffy*negy(ns)
            fock(lriz+ij)=ffz
            fock(lriz+ij1)=fct*ffz*negz(ns)
  340    continue
  350 continue
c
      end
c===================================================================
c  DFT POALRIZABILITY CPHF MODULE            KW January 2008
c  modification of Jon's CPHFDFT 
c
c      _ef stands for electric field
c ==================================================================
      subroutine cphfdft_ef(natoms,npos,
     *                      neqatm,
     *                      nsym,ngen,NQ,nsymop,iunq,
     *                      DEN1,  Den1B,  FDA,    FDB,
     $                      thrsh,  lsemi,  NMem,   Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C  ..........................................................
C  ** DFT Poalriz PROGRAM **
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
C  npos    -  starting position of current first atom
c  nsym    -  number of symmetry operations (highest abelian)
c  ngen    -  number of symmetry generators
C  NQ      -  number of symmetry-unique atoms
c  nsymop  -  list of symmetry operations
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
      DIMENSION FDA(    *      ),DEN1(    *      )
      DIMENSION FDB(    *      ),DEN1B(    *      )
      DIMENSION IUNQ(NQ),Z(NMem)
c  ...................................................................
      dimension nsymop(*)     ! ngen
      dimension neqatm(natoms,nsym)
c  ...................................................................
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
C  Now get the memory
C
      MaxGrid = NRad*NAng
c
      IMem = 7*NAtoms + NQ + NBas + 1
      NScr = NAtoms
      If(nslv.eq.0) Then
       IMem = IMem + NAtoms + 2*NAtoms**2 + 4*1130 + 4*MaxGrid +
     $               NPrim + NBas
c -- scratch memory for closed-shell local DFT
       NScr = 2*NBas*NBatch + 6*NBatch + 1 + 3*1
       if(.not.rhf) NScr = NScr + 4*NBatch + 9*1 !local uhf
c -- additional scratch memory for non-local DFT
       If(dft.GT.3)then
         NScr = NScr + 3*NBas*NBatch + 13*NBatch + 24*1
         if(.not.rhf) NScr = NScr + 8*NBatch + 24*1!non-local uhf
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
      IEnd = incex + 1
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
C  ----------------------------------------------------------------------
      CALL CPHFMAIN_ef(dft, npos,    NAtoms, nq, Z(iuq),
     *                       nsym, ngen, nsymop,neqatm,
     *              Z(IXC),          Z(IAN), 
     *              NRad,
     $              NAng,    NBatch,  Z(IDST), Z(IAIJ), Z(IRST),
     $              Z(IXXA), Z(IWTA), Z(IGRD), Z(IWGT), thrsh,
     $              NBas,    NShell,  NPrim,  bl(ibas),bl(ictr),
     $              Z(IBL),  Z(IExp), NAlpha,  NBeta,   rhf,
     $              bl(ida), bl(idb), DEN1,    DEN1B,   Z(IDM),
     $                                FDA,     FDB,     NScr,
     $              Z(IScr), lsemi)
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      RETURN
      END
c ========================================================================
c
      SUBROUTINE CPHFMAIN_ef(dft, npos,   NAtoms,  nq, iunq,
     *                       nsym, ngen, nsymop,neqatm,
     *                    XNuc,           IAN, 
     &                    NRad,
     $                    NAng,   NBatch, DISTN,  AIJ,    RDIST,
     $                    XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                    NBas,   NShell, NPrim,  BASDAT, INX,
     $                    BL,     ExpMIN, NAlpha, NBeta,  rhf,
     $                    DA,     DB,     DEN1,   DEN1B,  DM1,
     $                                    FDA,    FDB,    NScr,
     $                    Z,      lsemi)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for DFT contributions to Hessian
C  and derivative Fock matrices
C
C
      dimension neqatm(natoms,nsym)
      dimension nsymop(*)
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),XGRID(3,*),WGHT(*),
     $          BASDAT(13,*),INX(12,*),BL(*),ExpMIN(NShell),
     $          da(*),fda(*),db(*),fdb(*),den1(*),den1b(*),dm1(*)
      DIMENSION IUNQ(NQ)
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
      natdo=1
      compact=.false.
c
c -- form maximum 1st-order density element per column
c    moved this here from <cphfdft> (MM, October 2003)
c
      if(rhf)then
        CALL FormMaxD1(NBas,DEN1,DM1)
      else
        CALL FormMaxD1U(NBas,DEN1,DEN1B,DM1)
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
c
C ................................................................
C  get data from old texas depository
C
      call getival('nslv',nslv)           ! number of slaves (parallel)
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
      CALL PreGRID(natoms, xnuc,   IAN,    XGRID,  WGHT,
     $             DISTN,  RDIST,  AIJ,    XXA,    WTA)
C
C  get array of smallest exponent per shell
C
      CALL GetEXP(NPrim,NShell,BASDAT,INX,ExpMIN)
C
C  precompute shell normalization factors
C
      CALL AOInit(NShell,BASDAT,INX,BL)
ckw
c         write(6,*)' nq=',nq,' natoms=',natoms
ckw
c
        lsemi = lsemi+1
        lgrid = 0
        If(lsemi.EQ.1) lgrid=1
c
        DO 100 IAtom=1,NQ
        If(IAtom.EQ.1) Then
          IEntry = 1
        Else If(IAtom.EQ.NQ) Then
          IEntry = -1
        EndIf
        ICntr = IUNQ(IAtom)
c
        If(rhf) Then
          CALL DFTCphfC_ef(dft,    ICntr,  natoms,
     *                  nsym, ngen, nsymop, neqatm,
     $                  XNuc,           IAN,    IPRNT,  lrad,
     $                  lang,   IradQ,  factor, NBatch, DISTN,
     $                  AIJ,    RDIST,  XXA,    WTA,    XGRID,
     $                  WGHT,   thrsh,  NBas,   NShell, BASDAT,
     $                  INX,    BL,     ExpMIN, IdWt,   DA,
     $                  DEN1,   DM1,    NScr,   Z,      FDA,
     $                  lgrid,  IEntry)
        Else
ckw       CALL DFTCphfU_ef(dft,    ICntr,  NQ,
ckw  $                  XNuc,   XSym,   IAN,    IPRNT,  lrad,
ckw  $                  lang,   IradQ,  factor, NBatch, DISTN,
ckw  $                  AIJ,    RDIST,  XXA,    WTA,    XGRID,
ckw  $                  WGHT,   thrsh,  NBas,   NShell, BASDAT,
ckw  $                  INX,    BL,     ExpMIN, IdWt,   DA,
ckw  $                  DB,     DEN1,   DEN1B,  DM1,    NScr,
ckw  $                  Z,      FDA,    FDB,    lgrid,  IEntry)
        EndIf
 100    CONTINUE
cc
      ELSE
cc
c  PARALLEL MODE
c  -------------
c
ckw     CALL Para_DFTCphf_ef(dft,    NQ,                  XNuc,
ckw  $                    XSym,   IAN,    NQ,     IUNQ,   lrad,
ckw  $                    lang,   IradQ,  factor, NBatch, thrsh,
ckw  $                    NBas,   NShell, BASDAT, INX,    IdWt,
ckw  $                    rhf,    DA,     DB,     DEN1,   DEN1B,
ckw  $                    DM1,    FDA,    FDB,    lgrid)
cc
      ENDIF
c
      lgrid = 0
c
C ------------------------------------------------------------------
c     If(IPRNT.GT.4) Then
c       call getival('iout',iout)
c       write(iout,*) ' 1st-order density matrices are:'
c       do i=1,natonce
c       write(6,*) ' Atom:',i,' X-derivative density matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (den1(1,i,kk+l),l=1,k)
c       enddo
c       if(.not.rhf)then
c       write(6,*) ' Atom:',i,' beta X-derivative density matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (den1B(1,i,kk+l),l=1,k)
c       enddo
c       endif
c
c       write(6,*) ' Atom:',i,' Y-derivative density matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (den1(2,i,kk+l),l=1,k)
c       enddo
c       if(.not.rhf)then
c       write(6,*) ' Atom:',i,' beta Y-derivative density matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (den1B(2,i,kk+l),l=1,k)
c       enddo
c       endif
c
c       write(6,*) ' Atom:',i,' Z-derivative density matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (den1(3,i,kk+l),l=1,k)
c       enddo
c       if(.not.rhf)then
c       write(6,*) ' Atom:',i,' beta Z-derivative density matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (den1B(3,i,kk+l),l=1,k)
c       enddo
c       endif
c       enddo
c
 1112   format(1x,6(2x,F10.6))
c       write(iout,*) ' Partial Fock matrices are:'
c       do i=1,natonce
c       write(6,*) ' Atom:',i,' X-derivative Fock matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (fda(1,i,kk+l),l=1,k)
c       enddo
c       if(.not.rhf)then
c       write(6,*) ' Atom:',i,' beta X-derivative Fock matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (fdb(1,i,kk+l),l=1,k)
c       enddo
c       endif
c
c       write(6,*) ' Atom:',i,' Y-derivative Fock matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (fda(2,i,kk+l),l=1,k)
c       enddo
c       if(.not.rhf)then
c       write(6,*) ' Atom:',i,' beta Y-derivative Fock matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (fdb(2,i,kk+l),l=1,k)
c       enddo
c       endif
c
c       write(6,*) ' Atom:',i,' Z-derivative Fock matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (fda(3,i,kk+l),l=1,k)
c       enddo
c       if(.not.rhf)then
c       write(6,*) ' Atom:',i,' beta Z-derivative Fock matrix is:'
c       do k=1,nbas
c       kk = (k*(k-1))/2
c       write(6,1112) (fdb(3,i,kk+l),l=1,k)
c       enddo
c       endif
c       enddo
c     EndIf
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
      SUBROUTINE FormMaxD1(NBas,DEN1,DM1)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum 1st-order density matrix element per column.
C  Used for density threshold in DFT code
C  (NOTE: Used for density AND difference-density matrices)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  DEN1    -  1st-order density matrix (lower triangle)
C  DM1     -  on exit contains maximum 1st-order density element per column
C
C
      dimension den1(nbas*(nbas+1)/2,3), dm1(nbas)
c
      do i=1,nbas
         dm1(i)=0.0d0
      enddo
c
      do k=1,3
       ij=0
       do i=1,nbas
          do j=1,i
             ij=ij+1
             dx=abs(den1(ij,k))
             if(dx.gt.dm1(i)) dm1(i)=dx
             if(dx.gt.dm1(j)) dm1(j)=dx
          enddo
       enddo
      enddo
C
      RETURN
      END
c ========================================================================
      SUBROUTINE FormMaxD1U(NBas,DEN1A,DEN1B,DM1)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum 1st-order density matrix element per column.
C  Used for density threshold in DFT uhf code
C  (NOTE: Used for density AND difference-density matrices)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  DEN1A   -  alpha 1st-order density matrix (lower triangle)
C  DEN1B   -  beta 1st-order density matrix (lower triangle)
C  DM1     -  on exit contains maximum 1st-order density element per column
C
C
c     DIMENSION DEN1A(NAt3,NBas*(NBas+1)/2),DM1(NBas)
c     DIMENSION DEN1B(NAt3,NBas*(NBas+1)/2)
c
      dimension den1a(nbas*(nbas+1)/2,3), dm1(nbas)
      dimension den1b(nbas*(nbas+1)/2,3)
c
c
      do i=1,nbas
         dm1(i)=0.0d0
      enddo
c
      DO K=1,3
      ij=0
      do i=1,nbas
         do j=1,i
            ij=ij+1
            dx=max( abs(den1a(ij,k)),abs(den1b(ij,k))) 
            if(dx.gt.dm1(i)) dm1(i)=dx
            if(dx.gt.dm1(j)) dm1(j)=dx
         enddo
      enddo
      ENDDO
C
      END
c======================================================================
      SUBROUTINE DFTCphfC_ef(dft,    ICntr,  NAtoms,               
     *                  nsym, ngen, nsymop, neqatm,
     $                    XNuc,           IAN,    IPRNT,  NRad,
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
C             19 - O3LYP   (VWN  + OPTX     + LYP  + HF exchange)
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
c
C  NAtoms  -  total number of symmetry-unique atoms
c
C  NAtoms  -  total number of  ALL atoms
c
C  XNuc    -  nuclear coordinates
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
      dimension neqatm(natoms,nsym) 
      dimension nsymop(nsym)
      DIMENSION XNuc(3,NAtoms),               IAN(NAtoms),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          BASDAT(13,*),INX(12,*),
     $          da(*), den1(*),dm1(*),fda(*)
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
      if(IPrnt.gt.1) then
        write(6,*) '** DFT thresholds  thrsh:',thrsh,' thint:',thint
      end if
C
C  allocate scratch pointers
C
      iao = 1                 ! AO values over grid points
      inb = iao + NBasS       ! index array for "non-zero" AOs
      ibf = inb + NBasS       ! number of "non-zero" AOs per point
      id  = ibf + NBatch+1    ! density per grid point
      id1 = id  + NBatch      ! 1st-order density per grid point
      ipra  = id1 + 3         ! potential per grid point
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
        id1x2= id1x1+ 3         ! in situ first order gradient invariant
        id1x3= id1x2+ 3         !
        idg1 = id1x3+ 3         !
        isv  = idg1 + 3         ! V coefficients in numerical quadrature
        isw1 = isv  + 3         ! W coefficients in numerical quadrature
        isw2 = isw1 + 3         ! Func. deriv. w.r.t. alpha grad.
        isw3 = isw2 + 3         !
        ipga = isw3 + 3         !
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
      IF(lgrid.GT.0) THEN
        CALL GRID(NAtoms, xnuc,   IAN,    IPRNT,  NRad,
     $            NAng,   factor, IradQ,  Z(i1),  Z(i2),
     $            Z(i3),  DISTN,  RDIST,  AIJ,    ICntr,
     $            XXA,    WTA,    NPoint, XGRID,  WGHT)
c
c
c
C       symmetry sort of grid
        If(NSym.GT.0) Then
           XX = XNuc(1,ICntr)
           YY = XNuc(2,ICntr)
           ZZ = XNuc(3,ICntr)
           CALL SymGRID(natoms, icntr,  xx,     yy,     zz,
     $               nsym,   ngen, nsymop, neqatm , NPoint,
     $               XGRID,  WGHT)
        EndIf
C
C -- save grid to disk
c
        CALL WrGRID(Filename,len1,NPoint,XGRID,WGHT)
      ELSE
C       read existing file back from disk
        CALL RdGRID(Filename,len1,NPoint,XGRID,WGHT)
      ENDIF
c
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
C  Accumulate DFT contribution into partial derivative Fock matrices
C
      call secund(t1)
      CALL makecphfc_ef(dft,      npp,      nbas,     z(ibf), 
     $                         thrsh,    wght(ip), z(iprara),z(iprarb),
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
c======================================================================
c January 2008, KW : this is modification of the MakeCphfC routine
c                    used for hessian
c This is used for CPHF with an electric field perturbation
c
c                    _ef : stends for electric field
c======================================================================
      SUBROUTINE MakeCphfC_ef(dft,    NPP,    NBas,   nbf,   
     $                             thrsh,  WGHT,   prara,  prarb,
     $                     pga,    pgc,    praga,  pragb,  pragc,
     $                     pgaga,  pgagb,  pgagc,  pgcgc,  DEN1,
     $                     DM1,    DENX,   VAO,    VAOX,   INB,
     $                     VM,     INDX,   Den1P,  Den1X1, Den1X2,
     $                     Den1X3, G1P,    SV,     SW1,    SW2,
     $                     SW3,    FDA,    td1,    tg1,    tquad,
     $                     tsw)
      IMPLICIT REAL*8(A-H,O-Z)
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
      INTEGER dft,nbf(*),INB(*),INDX(NPP)
      DIMENSION DENX(3,*)
ckw
c     DIMENSION DEN1(3*NAtoms,*),DM1(*)
c     DIMENSION Den1P(3*NAtoms), G1P(3*NAtoms)
c     DIMENSION Den1X1(3*Natoms),Den1X2(3*Natoms),Den1X3(3*Natoms)
c     DIMENSION SV(3*NAtoms),SW1(3*NAtoms),SW2(3*NAtoms),SW3(3*NAtoms)
c     DIMENSION FDA(3*NAtoms,NBas*(NBas+1)/2)
ckw
      dimension den1(nbas*(nbas+1)/2,3),DM1(*)
      dimension den1p(3), g1p(3)
      dimension den1x1(3),den1x2(3),den1x3(3)
      dimension sv(3),sw1(3),sw2(3),sw3(3)
      dimension fda(NBas*(NBas+1)/2,3)
ckw
c
      dimension pga(npp),pgc(npp)
      dimension prara(npp),prarb(npp),praga(npp),pragb(npp),pragc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0)
      PARAMETER  (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d+17)
C
C
      NAt3 = 3
      natd3=3
c
      IF(dft.LE.3) THEN
C
C  local only
C
        DO 50 IP=1,NPP
        IPP = INDX(IP)
        VMx = VM(IPP)
        rara = prara(IP)*WGHT(IPP)
        CALL ZeroIT(Den1P,3)
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
ckw       DO IA=1,natd3
          DO IA=1,3
            Den1P(IA) = Den1P(IA) + DEN1(IJ,ia)*VAIJ
          EndDO
        EndIf
 19     CONTINUE
        Val2 = half*val*VAO(I)
        If(DMx*Abs(Val2).GT.thrsh1) Then
          IJ = IT+II
          DO IA=1,3
            Den1P(IA) = Den1P(IA) + DEN1(IJ,ia)*Val2
          EndDO
        EndIf
 20     CONTINUE
        call secund(t2)
        td1=td1+t2-t1
c
c -- get maximum value of 1st-order density at this grid point
        call secund(t1)
        Call absmax(3,Den1P,iixx,DMaxyz)
        DMax1 = DMaxyz*VMx
        If(DMax1*VMx.LT.thrsh) GO TO 49
c
c -- now do numerical integration
c
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
            DO 29 IA=1,3
              FDA(IJ,ia) = FDA(IJ,ia) + VAIJ*Den1P(IA)
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
          CALL ZeroIT(Den1P,3)
          CALL ZeroIT(Den1X1,3)
          CALL ZeroIT(Den1X2,3)
          CALL ZeroIT(Den1X3,3)
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
              DO IA=1,3
                Den1P(IA) = Den1P(IA) + DEN1(IJ,ia)*VAIJ
              EndDO
              DO IA=1,3
                Den1X1(IA) = Den1X1(IA)+DEN1(IJ,ia)*V1X
              EndDO
              DO IA=1,3
                Den1X2(IA) = Den1X2(IA)+DEN1(IJ,ia)*V1Y
              EndDO
              DO IA=1,3
                Den1X3(IA) = Den1X3(IA)+DEN1(IJ,ia)*V1Z
              EndDO
 69         CONTINUE
            ValI2=ValI*ValI*Half
            V1X=ValXI*ValI
            V1Y=ValYI*ValI
            V1Z=ValZI*ValI
            valt=max(Abs(vali2),Abs(V1X),Abs(V1Y),Abs(V1Z))
            If(Valt.LT.thtest) GO TO 70
            IJ = IT+II
            DO IA=1,3
              Den1P(IA) = Den1P(IA) + DEN1(IJ,ia)*ValI2
            EndDO
            DO IA=1,3
              Den1X1(IA) = Den1X1(IA)+DEN1(IJ,ia)*V1X
            EndDO
            DO IA=1,3
              Den1X2(IA) = Den1X2(IA)+DEN1(IJ,ia)*V1Y
            EndDO
            DO IA=1,3
              Den1X3(IA) = Den1X3(IA)+DEN1(IJ,ia)*V1Z
            EndDO
 70       CONTINUE
        call secund(t2)
        td1=td1+t2-t1
c
c -- form 1st-order gradient invariant.
c
          DO IA=1,3
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
          DO IA=1,3
            sv(ia)=rara*den1p(ia)+prg*g1p(ia)
          EndDO
          DO IA=1,3
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
          Call absmax(3,SV,isv,svmax)
          Call absmax(3,SW1,isw,swmax1)
          Call absmax(3,SW2,isw,swmax2)
          Call absmax(3,SW3,isw,swmax3)
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
              DO IA=1,3
                FDA(IJ,ia) = FDA(IJ,ia) + (SV(IA)*VIJ+
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
c======================================================================
