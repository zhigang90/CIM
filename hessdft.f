      SUBROUTINE HESSDFT(HESS,NMem,Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** DFT HESSIAN PROGRAM **
C
C  This Program determines the DFT exchange-correlation
C  contribution to the total Hessian and the derivative
C  Fock matrices, combining it with the Coulomb contribution
C  previously calculated.
C  ............................................................
C
      DIMENSION HESS(*),Z(NMem)
c  ...................................................................
c     common /big/ bl(1000)            ! old texas depository
c     common /intbl/maxsh,ixx(100)     ! old texas integer depository
c ....................................................................
      Integer dft
      Character jobname*256,cdum*20
      Logical rhf
C
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
c  ...................................................................
      Data nrad/400/, nang/434/, NBatch/50/
c  ...................................................................
c
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    dft functional type
C    total number of atoms (including dummies)
C    number of dummy atoms (if any)
C    multiplicity
C    number of alpha/closed-shell orbitals
C    number of beta orbitals
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call RdDFT(IUnit,dft)
      If(dft.gt.0) Then
        call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
        call fdcntrl(IUnit,7,'$ndummy',idum)
        If(idum.EQ.0) Then
          backspace IUnit
          READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
          Ndum = Ndum1+Ndum2
        Else
          Ndum = 0
        EndIf
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
C  DFT uses real atoms only
C
      NAtoms = NAtoms-Ndum
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
      if(.not.rhf)call getival('ldensB',idb)! beta density matrix
      call getival('nsym',nsym)           ! number of symmetry operations - 1
      call getival('ngener',ngen)         ! number of generators
      If(nsym.gt.0) Then
       call getival('SymNuPr1',nupr)      ! equivalent atoms
       call getival('SymFunPr',ifp)       ! equivalent basis functions
      Else
       nupr = 1
       ifp = 1
      EndIf
      call getival('wghtd',IdWt)          ! use of weight derivatives
      call getival('nslv',nslv)           ! number of slaves (parallel)
C --------------------------------------------------------------------
C
C
C  Allocate the DFT scratch memory
C
      ntri = (NBas*(NBas+1))/2
      NAt3 = 3*NAtoms
      MaxGrid = NRad*NAng
c
c  we compute the scratch memory requirement for the quadrature.
c  this memory will be actually allocated only in serial runs.
c  In parallel runs, the memory is not allocated by the master,
c  but we need these numbers in order to estimate the meory
c  requirements of the slave for the calculation of the number
c  of passses needed to compute all Fock derivatives
c
      imem = 0
      IMem = IMem + NAtoms + 2*NAtoms**2 + 4*1130 + 4*MaxGrid +
     $               NPrim + NBas
      NScr = 11*NBas*NBatch + 6*NBatch + 1 + NAt3*(NAt3+1)
c -- above memory requirement for closed-shell local DFT only
      if(.not.rhf)NScr=NScr+4*NBatch+ NAt3*(Nat3+3) ! open-shell local
c
c    memory for Non Local DFT
c
      if(dft.gt.3)then
        NScr=NScr+10*NBas*NBatch+13*NBatch+24*Natoms+36*NAtoms**2
        if(.not.rhf)NScr=NScr+9*NBatch+33*NAtoms+45*NAtoms**2 ! uhf
      endif
c
c    memory for weight derivatives
c
      if(IdWt.eq.1)then
         NScr=NScr+(NAt3*(NAt3+1)+2)*NBatch+20*NAtoms*NAtoms+NAtoms
      endif
      imemp=imem
      nscrp=nscr
c
c   now we compute the memory we actually have to allocate.
c   imemp and nscrp will contain the memory needed by the slaves
c   in case of a parallel run
c
      if(nslv.eq.0)then
        IMem = Imem + 5*NAtoms + 2*NBas
      else
        IMem =  5*NAtoms + 2*NBas
        IMemp = imemp + 5*NAtoms + 2*NBas
        Nscr = Natoms
      endif
C
C  make sure there is enough scratch memory for symmetrizing Hessian
C  if necessary, and for a buffer for accumulating the fock matrices
C
      NScr = MAX(NScr,3*NAt3**2 + 3*3*8,3*ntri)
      NScrp = MAX(NScrp,3*NAt3**2 + 3*3*8,3*ntri)
c
      IMem = IMem + NScr
      IMemp = IMemp + NScrp
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,7,'HESSDFT')
      CALL ZeroIT(Z,IMem)      ! clear the memory
C
C  Allocate the memory pointers
C
      IXC  = iptr                       !  nuclear coordinates
      IAN  = IXC + NAt3                 !  nuclear charge (atomic numbers)
      IUQ  = IAN + NAtoms               !  list of symmetry unique atoms
      IAtm = IUQ + NAtoms               !  atom/basis function ownership
      IDM  = IAtm + NBas                !  largest density element per column
      IEnd = IDM + NBas
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
      CALL MemCHK(IMem,IEnd,7,'HESSDFT')
c
c   compute how many components of the Fock derivative matrix
c   can be kept in memory in one pass, and compute the number
c   of passes needed.
c
      call npass_dfth(rhf,NMem,IEnd,imemp,Nslv,Natoms,Ntri,
     $                npass,natonce)
c
c  allocate memory for fock derivative matrices
c
      IFDA=IEnd+1
      IEnd=IFDA+3*natonce*ntri
      if(.not.rhf)then
        IFDB=IEnd
        IEnd=IFDB+3*natonce*ntri
      endif
      IEnd=IEnd-1
c -- assign memory for high water mark (see <forces>)
      call getmem(IEnd,lastx)
c
c -- get coordinates and atomic numbers
c    from Texas format to decent format
      call getnucdat(NAtoms,Z(IXC),Z(IAN))
c
c -- form maximum density element per column
      if(rhf)then
        CALL FormMaxDen(NBas,bl(ida),Z(IDM))
      else
        CALL FormMaxDenUH(NBas,bl(ida),bl(idb),Z(IDM))
      endif
C  ----------------------------------------------------------------------
C
      CALL HESSMAIN(dft,     NAtoms,  NPass,   NAtonce, Z(IXC),
     $              Z(IAN),  Z(IUQ),  bl(ifp), bl(nupr),NRad,
     $              NAng,    NBatch,  Z(IDST), Z(IAIJ), Z(IRST),
     $              Z(IXXA), Z(IWTA), Z(IGRD), Z(IWGT), NBas,
     $              NShell,  NPrim,   bl(ibas),bl(ictr),Z(IBL),
     $              Z(IExp), Z(IAtm), NAlpha,  NBeta,   rhf,
     $              bl(ida), bl(idb), Z(IDM),  Z(IFDA), Z(IFDB),
     $              HESS,    NScr,    IdWt,    Z(IScr))
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      RETURN
      END
c ========================================================================
c
      SUBROUTINE HESSMAIN(dft,    NAtoms, NPass,  Natonce,XNuc,
     $                    IAN,    IUNQ,   NEqBAS, NEqATM, NRad,
     $                    NAng,   NBatch, DISTN,  AIJ,    RDIST,
     $                    XXA,    WTA,    XGRID,  WGHT,   NBas,
     $                    NShell, NPrim,  BASDAT, INX,    BL,
     $                    ExpMIN, NBAtm,  NAlpha, NBeta,  rhf,
     $                    DA,     DB,     DM,     FDA,    FDB,
     $                    HESS,   NScr,   IdWT,   Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for DFT contributions to Hessian
C  and derivative Fock matrices
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),XGRID(3,*),WGHT(*),
     $          BASDAT(13,*),INX(12,*),BL(NPrim),ExpMIN(NShell),
     $          NBAtm(NBas),DA(*),DB(*),FDA(*),FDB(*),DM(NBas),
     $          HESS(3*NAtoms,3*NAtoms)
      DIMENSION IUNQ(NAtoms),NEqBAS(7,NBas),NEqATM(NAtoms,*)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      Character jobname*256,cdum*20
      INTEGER dft
      Logical rhf
c  ...................................................................
      common /symm/nsym,nsy(7)         ! old texas symmetry data
c ..............................................................
      DIMENSION Z(NScr)
c  ...................................................................
      data lrad/0/, lang/0/, IradQ/0/, thrsh/1.0d-12/
c  ...................................................................
C
      PARAMETER (IUnit=1)
c
      Common /job/jobname,lenJ
C
C
C  initialization
C
      call secund(t1)
      call elapsec(elahess1)
      call getival('iout',iout)
ckw   WRITE(IOut,1000)
      NAt3 = 3*NAtoms
      ntri = (NBas*(NBas+1))/2
      EL=0.0d0
C
C  If there is symmetry, we need to construct and symmetrize the
C  direct DFT contribution to the Hessian separately
C
      If(nsym.GT.0) Then
C
C  save existing Hessian to file
C
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.temp',FORM='UNFORMATTED',
     $       STATUS='UNKNOWN')
       WRITE(40) HESS
       CLOSE (UNIT=40,STATUS='KEEP')
C
C  zero out Hessian
C
       CALL ZeroIT(HESS,NAt3**2)
      EndIf

C
C ...............................................................
C  read from the <control> file
C    grid factor
C    print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$factor',2,idum,factor,cdum)
c     call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      call getival('printh',IPRNT)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C ...............................................................
C
C ** WARNING **   POSSIBLE MODIFICATION OF GRID FACTOR
C    For accurate frequencies, need quadrature weight derivatives.
C    These are currently only available for systems with NO SYMMETRY.
C    If there is symmetry, we automatically increase the grid factor
C    (which controls the number of radial points) to compensate
C
c     If(nsym.GT.0) Then
c       IdWt = 0
c       factor0 = factor
c       If(factor.LT.2.0d0) Then
c         factor = 3.0d0*factor
c         WRITE(IOut,1100) factor0,factor
c       EndIf
c     EndIf
C .................................................................
C
C  get data from old texas depository
C
      call getival('nslv',nslv)           ! number of slaves (parallel)
      call getival('ngener',ngen)         ! number of generators
C
C  get number and location of symmetry-unique atoms under
C  highest abelian point group
C
      CALL GetUNQ(NAtoms, nsym,   XNuc,   NEqATM, NQ,
     $            IUNQ,   Z)
C
C  need basis function/associated atom list
C
      CALL SortBATM(NBas,NShell,INX,NBAtm)
cc      write(6,*) ' basis functions/atom ownership is:'
cc      do i=1,nbas
cc      write(6,*) i,nbatm(i)
cc      enddo
C ...............................................................
C
      IF(nslv.EQ.0) THEN
C
C  initialization for single processor mode
C
C  calculate inverse atomic distances, Becke aij parameters
C  and angular grid and weights prior to full grid construction
C
      CALL PreGRID(NAtoms, XNuc,   IAN,    XGRID,  WGHT,
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
      ENDIF
c
c  multiple passes implementation:
c
c  the variable npass contains the number of passes necessary
c  to compute all components of the derivativative Fock matrices.
c  the variable natonce contains the maximum number of Fock
c  derivative components that can be kept in memory at the same time
c
c  if more than one pass, try to have balanced passes
c
      np1=0
      natone=natonce
      if(npass.gt.1)then
        natone=natoms/npass
        np1=natoms-natone*npass
      endif
ckw
      if(npass.eq.1) write(iout,200)
  200 format(/'DFT contributions to Fder will be calculated in 1 pass')
      if(npass.gt.1) write(iout,201) npass
  201 format(/'DFT contributions to Fder will be calculated in ',i2,
     *      ' passes')
      write(iout,*) '  '
      call f_lush(iout)
ckw
C
C  Now calculate the various DFT contributions
C
C  Nb=first component to compute in the current pass
C  Ne=last component to compute in current pass
C  Natcurr=Ne-Nb+1=total number of components computed by current pass
C
      nb=1
      call getival('iout',iout)
      do ipass=1,npass
        if(ipass.eq.npass)then
          ne=natoms
        else if(ipass.le.np1)then
          ne=nb+natone
        else
          ne=nb+natone-1
        endif
        natcurr=ne-nb+1
        if(natcurr.gt.natonce)then
          write(6,*)'HESSDFT: too many atoms in current pass!'
          write(6,*)'ipass=',ipass,' natcurr=',natcurr,
     $              'natonce=',natonce
          call nerror(1,'HESSDFT',' too many atom in current pass',0,0)
        endif
ckw
c       write(iout,
c    $  '(''Pass'',i4,'', computing'',i4,'' Fock derivatives'//
c    $  ', atom'',i4,'' through'',i4)')ipass,natcurr,nb,ne
c       call f_lush(iout)
ckw
        call zeroit(FDA,3*NAtcurr*ntri)
        if(.not.rhf)call zeroit(FDB,3*NAtcurr*ntri)
        IF(nslv.EQ.0) THEN
cc
c  SINGLE PROCESSOR MODE
c  ---------------------
          DO 100 IAtom=1,NQ
          If(IAtom.EQ.1) Then
            IEntry = 1
          Else If(IAtom.EQ.NQ) Then
            IEntry = -1
          EndIf
          ICntr = IUNQ(IAtom)
c
          If(rhf) Then
            CALL DFTHessC(dft,    ICntr,  NAtoms, Nb,     Ne,
     $                    XNuc,   IAN,    nsym,   ngen,   nsy,
     $                    NEqATM, IPRNT,  lrad,   lang,   IradQ,
     $                    factor, NBatch, DISTN,  AIJ,    RDIST,
     $                    XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                    NBas,   NShell, BASDAT, INX,    BL,
     $                    ExpMIN, IdWt,   NBAtm,  DA,     DM,
     $                    NScr,   Z,      FDA,    HESS,   EL,
     $                    IEntry)
          Else
            CALL DFTHessU(dft,    ICntr,  NAtoms, Nb,     Ne,
     $                    XNuc,   IAN,    nsym,   ngen,   nsy,
     $                    NEqATM, IPRNT,  lrad,   lang,   IradQ,
     $                    factor, NBatch, DISTN,  AIJ,    RDIST,
     $                    XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                    NBas,   NShell, BASDAT, INX,    BL,
     $                    ExpMIN, IdWt,   NBAtm,  DA,     DB,
     $                    DM,     NScr,   Z,      FDA,    FDB,
     $                    HESS,   EL,     IEntry)
          EndIf
c
 100      CONTINUE
cc
        ELSE
cc
c  PARALLEL MODE
c  -------------
c
          CALL Para_DFTH(dft,    NAtoms, ipass,  npass,  Nb,
     $                   Ne,     XNuc,   IAN,    nsym,   ngen,
     $                   nsy,    NEqATM, NQ,     IUNQ,   IPRNT,
     $                   lrad,   lang,   IradQ,  factor, NBatch,
     $                   thrsh,  NBas,   NShell, BASDAT, INX,
     $                   NBAtm,  IdWt,   rhf,    DA,     DB,
     $                   DM,     FDA,    FDB,    HESS,   EL)
cc
        ENDIF
c
c  add DFT contribution to the fock derivatives saved on file
c
        if(rhf)then
          call add_trsp_dsk(61,ntri,Nb,Ne,fda,z)
        else
          call add_trsp_dsk(61,ntri,Nb,Ne,fda,z)
          call add_trsp_dsk(71,ntri,Nb,Ne,fdb,z)
        endif
        nb=ne+1                 ! ready for next pass
      ENDDO
C
C ------------------------------------------------------------------
C
      If(NSym.GT.0) Then
C
C  if symmetry was used, correct EL
C
       EL = EL*DFloat(NSym+1)
C
C  symmetrize the partial Fock derivative matrices that
c  are now stored on file
c  fdersym_hes1 needs a memory buffer of 6*ntri locations.
c  this memory is available, because one area of 3*ntri has
c  been reserved in the scratch memory Z, and at least another
c  area of 3*ntri is allocated in fda. As the memory of z and
c  fda is in fact contiguous, we use z as buffer.
C
       if(rhf)then
       call fdersymm_hes1(61,ngen,NBas,ntri,NAtoms,Z,nsy,NEqBas,NEqATM)
       else
       call fdersymm_hes1(61,ngen,NBas,ntri,NAtoms,Z,nsy,NEqBas,NEqATM)
       call fdersymm_hes1(71,ngen,NBas,ntri,NAtoms,Z,nsy,NEqBas,NEqATM)
       endif
C
C  symmetrize Hessian contribution
C
       itrn = 1
       ih   = itrn + 9*nsym
       ir   = ih + NAt3**2
       iu   = ir + NAt3**2
c
       CALL hdersymm_hes(nsym,   NAtoms, nsy,    NEqATM, Z(itrn),
     $                   HESS,   Z(ir),  Z(iu),  Z(ih))
C
C  symmetrized Hessian is in Z(ih)
C  now add to existing contribution
C
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.temp',FORM='UNFORMATTED',
     $       STATUS='OLD')
       READ(40) HESS
       CLOSE (UNIT=40,STATUS='DELETE')
       CALL AddVEC(NAt3**2,HESS,Z(ih),HESS)
      EndIf
c
      If(IPRNT.GT.3) Then
        call getival('iout',iout)
        write(iout,*) ' Direct Hessian contribution is:'
        do i=1,natoms
        ii = 3*(i-1)+1
        do j=1,natoms
        jj = 3*(j-1)+1
        write(iout,*) ' Atom: ',i,' with Atom: ',j
        write(iout,1112) hess(jj,ii),hess(jj,ii+1),hess(jj,ii+2)
        write(iout,1112) hess(jj+1,ii),hess(jj+1,ii+1),hess(jj+1,ii+2)
        write(iout,1112) hess(jj+2,ii),hess(jj+2,ii+1),hess(jj+2,ii+2)
 1112   format(1x,6(2x,F10.6))
        enddo
        enddo
      EndIf
      If(IPRNT.GT.4) Then
        if(rhf)then
          write(iout,*) ' Fock-Derivative matrices are:'
          do i=1,natoms
            call prt_fder(61,iout,i,nbas,ntri,fda)
          enddo
        else
          write(iout,*) ' Alpha Fock-Derivative matrices are:'
          do i=1,natoms
            call prt_fder(61,iout,i,nbas,ntri,fda)
          enddo
          write(iout,*) ' Beta Fock-Derivative matrices are:'
          do i=1,natoms
            call prt_fder(71,iout,i,nbas,ntri,fda)
          enddo
        endif
      EndIf
C  ----------------------------------------------------------------
C
      call secund(t2)
      call elapsec(elahess2)
      t1 = (t2-t1)/60.0d0
      elaps = (elahess2-elahess1)/60.0d0
      WRITE(6,1200) t1,elaps
      If(IPRNT.GT.1) write(6,*) ' Number of electrons over grid:',el
c
      RETURN
c
 1000 FORMAT(/,' Calculating DFT Contribution to Hessian')
 1100 FORMAT('**WARNING** DFT Weight Derivatives Currently Unavailable',
     $       ' with Symmetry',/,'  Increasing radial grid factor from ',
     $       F6.4,' to ',F6.4,' to Compensate')
 1200 FORMAT('Master CPU time for XC part of Hessian  ='
     *,f8.2,' Elapsed = ',f8.2,' min'/)
c
      END
c ========================================================================
      SUBROUTINE FormMaxDenUH(NBas,DA,DB,DM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum density matrix element per column.
C  This version stores the maximum element of alpha or beta
C  density, not their sum.
C  Used for density threshold in UHF/DFT Hessian code
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta density matrix (lower triangle)
C  DM      -  on exit contains maximum density element per column
C
C
      DIMENSION DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NBas)
C
      DO 20 I=1,NBas
      DMx = 0.0d0
      II = (I*(I-1))/2
      DO 10 J=1,I
      IJ = II+J
      DMx = max(DMx,Abs(DA(IJ)),Abs(DB(IJ)))
 10   CONTINUE
      DO 11 J=I+1,NBas
      IJ = (J*(J-1))/2 + I
      DMx = max(Dmx,Abs(DA(IJ)),Abs(DB(IJ)))
 11   CONTINUE
      DM(I) = DMx
 20   CONTINUE
C
      RETURN
      END
c ======================================================================
      subroutine npass_dfth(rhf,NMem,IEnd,imemp,Nslv,Natoms,Ntri,
     $                      npass,natonce)
      implicit real*8 (a-h,o-z)
c
c   compute how many components of the fock derivative matrix
c   can be kept in memory in one pass
c
c   rhf        logical flag for rhf/uhf
c   nmem       total memory available for DFT run
c   iend       end of assigned memory (free memory = nmem-iend-1)
c   imemp      memory needed by slaves in parallel run
c   Nslv       Number of slaves (different from zero in parallel runs)
c   natoms     number of atoms
c   ntri       dimension of a triangular matrix
c   npass      in exit will contain the number of passes necessary
c              to compute all the components
c   natonce    in exit will contain  the maximum number of atoms
c              to compute in one pass
c
      logical rhf
c
c   compute free memory:
c
      nfree=nmem-iend-1
c
c   in case of parallel runs, the free memory has to be computed
c   based on the memory requirements of the slaves also, as the master
c   does not allocate the scratch memory needed for the dft quadrature.
c
      if(nslv.ne.0)then
       nfreep=nmem-imemp
       if(nfreep.le.0)then
         write(6,*)'HessDFT: not enough memory for slaves:'
         write(6,*)'Available',nmem*8,', needed',imemp*8,' bytes'
         call nerror(2,'HESSDFT','not enough memory for slaves',0,0)
       endif
       nfree=min(nfree,nfreep)
      endif
c
c   memory needed for one atom
c
      if(rhf)then
        noneat=3*ntri
      else
        noneat=2*3*ntri
      endif
c
c   if we cannot do at least one atom, we are in trouble
c
      if(nfree.lt.noneat)then
        write(6,*)'Not enough memory in HessDFT!'
        write(6,*)'Memory needed for one atom (bytes):',noneat*8
        write(6,*)'Free memory (bytes)               :',nfree*8
        write(6,*)'Please, increase memory allocation by at least',
     $             (noneat-nfree)*8,' bytes'
        call nerror(3,'HESSDFT','not enough memory',0,0)
      endif
c
c   compute how many atoms we can do in one pass
c
      nonepass=int(nfree/noneat)
      if(nonepass.ge.natoms)then   ! one pass only
        npass=1
        natonce=natoms
      else                         ! more passes needed.
        npass=int(natoms/nonepass)
        nrest=natoms-npass*nonepass
        if(nrest.gt.0)npass=npass+1
        natonce=nonepass
      endif
c
      return
      end
c======================================================================
      subroutine add_trsp_dsk(iunit,ntri,nb,ne,fd,z)
      implicit real*8 (a-h,o-z)
      dimension fd(3,nb:ne,ntri),z(ntri,3)
c
      do iat=nb,ne
        read(unit=iunit,rec=iat)z
        do i=1,ntri
          z(i,1)=z(i,1)+fd(1,iat,i)
          z(i,2)=z(i,2)+fd(2,iat,i)
          z(i,3)=z(i,3)+fd(3,iat,i)
        enddo
        write(unit=iunit,rec=iat)z
      enddo
c
      end
c======================================================================
      subroutine prt_fder(iunit,iout,iat,nbas,ntri,fd)
      implicit real*8 (a-h,o-z)
      dimension fd(ntri,3)
      read(unit=iunit,rec=iat)fd
      write(iout,*) ' Atom:',iat,' X-derivative Fock matrix is:'
      do k=1,nbas
        kk = (k*(k-1))/2
        write(iout,1112) (fd(kk+l,1),l=1,k)
      enddo
      write(iout,*) ' Atom:',iat,' Y-derivative Fock matrix is:'
      do k=1,nbas
        kk = (k*(k-1))/2
        write(iout,1112) (fd(kk+l,2),l=1,k)
      enddo
      write(iout,*) ' Atom:',iat,' Y-derivative Fock matrix is:'
      do k=1,nbas
        kk = (k*(k-1))/2
        write(iout,1112) (fd(kk+l,3),l=1,k)
      enddo
 1112 format(1x,6(2x,F10.6))
      end
