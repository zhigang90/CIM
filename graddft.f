c ==================================================================
c  DFT GRADIENT MODULE            JB   September 1997
c ==================================================================
c
      SUBROUTINE GRADDFT(NMem,Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** DFT GRADIENT PROGRAM **
C
C  This Program determines the DFT exchange-correlation
C  contribution to the total gradient, combining it with
C  the coulomb contribution previously calculated.
C  ............................................................
C
      DIMENSION Z(NMem)
      Integer dft
      Character jobname*256,cdum*20
      Logical rhf
C
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
c  ...................................................................
      Data nrad/400/, nang/434/, NBatch/50/, IdWt/0/
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
        call rdcntrl(IUnit,7,'$weight',1,IdWt,dum,cdum)
      EndIf
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  do we have a dft contribution?
C
      If(dft.EQ.0) RETURN
c
      rhf = (NBeta.EQ.0.AND.IMult.EQ.1)
C
C  Read from the <basis> file
C    number of basis functions, primitives and shells
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$nbasis',1,NBas,dum,cdum)
      call rdcntrl(IUnit,7,'$nshell',1,NShell,dum,cdum)
      call rdcntrl(IUnit,6,'$nprim',1,NPrim,dum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  DFT uses real atoms only
C
      NAtoms = NAtoms-Ndum
C
C  Now get the memory
C
      MaxGrid = NRad*NAng
      NScr = 4*NBas*NBatch + NAlpha*NBatch + 1
      If(rhf) Then
        NScr = NScr + MAX(NBas*NBatch + 4*NBatch, NAtoms*(NAtoms+2))
      Else
        NScr = NScr + MAX(NBas*NBatch + 6*NBatch, NAtoms*(NAtoms+2))
      EndIf
      If(dft.gt.3) NScr = NScr + 6*NBas*NBatch + 6*NBatch + 3*NAlpha
      If(.NOT.rhf.AND.dft.gt.3) NScr = NScr + 6*NBatch + NBatch*NBeta
c
      IMem = 14*NAtoms + 2*NAtoms*NAtoms + 4*1130 + 4*MaxGrid +
     $       14*NPrim + 12*NShell + 3*NBas*NBas + (NBas*(NBas+1))/2 +
     $       8*NBas
      If(.NOT.rhf) IMem = IMem + 3*NBas*NBas + (NBas*(NBas+1))/2
c
      IMem = IMem + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,8,'GRADDFT')
      CALL ZeroIT(Z,IMem)      ! clear the memory
C
C  Allocate memory pointers
C  **WARNING**  Need to do something about memory allocation
C    density matrics NOT needed till later if IdWt=1
C
      IXC  = iptr                      !  nuclear coordinates
      IAN  = IXC + 3*NAtoms            !  atomic numbers
      IQA  = IAN + NAtoms              !  atomic charge
      IUQ  = IQA + NAtoms              !  list of symmetry unique atoms
      IFP  = IUQ + NAtoms              !  list of basis function pairs
      IBas = IFP + 7*NBas              !  basis set data
      INX  = IBas + 13*NPrim           !  integer basis set data
      IAtm = INX + 12*NShell           !  no. basis functions per atom
      IDA  = IAtm + NAtoms             !  density matrix (lower triangle)
      IGX  = IDA  + (NBas*(NBas+1))/2  !  x DFT derivative matrix
      IGY  = IGX + NBas*NBas           !  y DFT derivative matrix
      IGZ  = IGY + NBas*NBas           !  z DFT derivative matrix
      IGWT = IGZ + NBas*NBas           !  DFT weight-derivatives
      IGC  = IGWT + 3*NAtoms           !  final DFT gradient
      IEnd = IGC + 3*NAtoms
c
      If(.NOT.rhf) Then
        IDB  = IEnd                      !  beta density matrix (lower triangle)
        IGXB = IDB  + (NBas*(NBas+1))/2  !  beta x DFT derivative matrix
        IGYB = IGXB + NBas*NBas          !  beta y DFT derivative matrix
        IGZB = IGYB + NBas*NBas          !  beta z DFT derivative matrix
        IEnd = IGZB + NBas*NBas
      EndIf
c
c
c -- the following arrays are NOT needed in parallel mode
      IDST = IEnd                      !  distance to nearest neighbour
      IAIJ = IDST + NAtoms             !  Becke grid parameters
      IRST = IAIJ + NAtoms*NAtoms      !  inverse interatomic distances
      IXXA = IRST + NAtoms*NAtoms      !  angular quadrature grid
      IWTA = IXXA + 3*1130             !  angular quadrature weights
      IGRD = IWTA + 1130               !  X,Y,Z grid points
      IWGT = IGRD + MaxGrid*3          !  grid quadrature weights
      IBL  = IWGT + MaxGrid            !  precomputed normalization factors
      IExp = IBL + NPrim               !  smallest exponent per basis function
      IScr = IExp + NBas               !  general scratch storage
      IEnd = IScr + NScr - 1
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,7,'GRADDFT')
c -- assign memory for high water mark (see <forces>)
      call getmem(IEnd,lastx)
C
C
C  ----------------------------------------------------------------------
C
      CALL GRADMAIN(dft,     NAtoms,  Z(IXC),  Z(IAN),  Z(IQA),
     $              Z(IUQ),  Z(IFP),  NRad,    NAng,    NBatch,
     $              Z(IDST), Z(IAIJ), Z(IRST), Z(IXXA), Z(IWTA),
     $              Z(IGRD), Z(IWGT), NBas,    NShell,  Z(IBas),
     $              Z(INX),  Z(IBL),  Z(IExp), Z(IAtm), NAlpha,
     $              NBeta,   rhf,     Z(IDA),  Z(IDB),  Z(IGX),
     $              Z(IGY),  Z(IGZ),  Z(IGXB), Z(IGYB), Z(IGZB),
     $              IdWt,    Z(IGWT), Z(IGC),  NScr,    Z(IScr))
      call retmem(1)
C
C  ----------------------------------------------------------------------
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE GRADMAIN(dft,    NAtoms, XNuc,   IAN,    QA,
     $                    IUNQ,   NEqBAS, NRad,   NAng,   NBatch,
     $                    DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                    XGRID,  WGHT,   NBas,   NShell, BASDAT,
     $                    INX,    BL,     ExpMIN, NBAtm,  NAlpha,
     $                    NBeta,  rhf,    DA,     DB,     GX,
     $                    GY,     GZ,     GXB,    GYB,    GZB,
     $                    IdWt,   GWT,    GC,     NScr,   Z)
      use memory, block => bl
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for DFT gradients
C
C
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),QA(NAtoms),XGRID(3,*),
     $          WGHT(*),BASDAT(13,*),INX(12,*),BL(*),ExpMIN(NShell),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),
     $          GX(NBas,NBas),GY(NBas,NBas),GZ(NBas,NBas),
     $          GXB(NBas,NBas),GYB(NBas,NBas),GZB(NBas,NBas),
     $          GWT(NAtoms,3),GC(3,NAtoms),NBAtm(NAtoms)
      DIMENSION IUNQ(NAtoms),NEqBAS(7,NBas)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
      Character*20 cdum
      Character*256 jobname
      INTEGER dft
      Logical rhf
c  ...................................................................
c     common /intbl/maxsh,ixx(100)     ! old texas integer depository
      common /symm/nsym,nsy(7)         ! old texas symmetry data
c ..............................................................
      DIMENSION Z(NScr)
c  ...................................................................
      data lrad/0/, lang/0/, IradQ/0/, thrsh/1.0d-12/
c  ...................................................................
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
      PARAMETER (IUnit=1)
c
      Common /job/jobname,lenJ
C
C
C  initialization
C
C     WRITE(6,1000)
      call secund(t1)
      call elapsec(elagrad1)
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
C
C  read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XNuc,-1,jnk)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  read the basis set data
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdbasis(IUnit,  NAtoms, AtSymb, XNuc,   Z(1),
     $             Z(1+NAtoms),BASDAT)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C ....................................................................
C  Be Careful with basis order!
C  We need basis ordered PER ATOM for DFT gradient
C  Wolinski ordering (for integral evaluation and hence for the MOs)
C  is PER SHELL (e.g. all S functions together regardless of atom)
C
C  The following routines are relevant here
C   SortBAS1  -  sorts basis per shell (but NOT full Wolinski)
C   reorder   -  orders basis from SortBAS1 into Wolinski
C   SortBAS2  -  orders basis per atom and supplies integer
C                ordering array relating atom and Wolinski order
C ....................................................................
C
      CALL SortBAS1(NAtoms,NShell,Z(1+NAtoms),INX)
      CALL normaliz(NShell,INX,BASDAT)
C
C  get number of basis functions per atom
C
      CALL BasATOM(NAtoms,NShell,Z(1+NAtoms),NBAtm)
cc      write(6,*) ' Number of basis functions per atom is:'
cc      do i=1,natoms
cc      write(6,*) i,nbatm(i)
cc      enddo
C
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
C  get actual atomic charges (there may, e.g., be ghost atoms)
C
      CALL GetAtChrg(NAtoms,QA)
C
C  need to reorder MOs as need basis functions ordered per atom
C  and MOs have Wolinski special order
C
      i1 = 1
      i2 = i1 + NBas
      i3 = i2 + NShell
      IEnd = i3 + 12*NShell - 1
      CALL MemCHK(NScr,IEnd,8,'GRADMAIN')
c
cc      Call preorder(NAtoms,NShell,Z(i3),INX)
      Call reorder(NShell,INX,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinski order
      Call SortBAS2(NShell,Z(i3),Z(i2),Z(i1),INX)        ! per atom
C
C  now read the MOs
C
      IF(rhf) THEN
        call rea(GX,NBas*NBas,4,'evec_rhf')
C
C  now reorder the occupied MOs, at the same time transpose them
C
        CALL ReorderMO(NBas,NAlpha,Z(i1),GX,GY)
c -- form density matrix
        Call FormDEN(NBas,NAlpha,GY,DA)
      ELSE
        call rea(GX,NBas*NBas,4,'evea_uhf')
        call rea(GXB,NBas*NBas,4,'eveb_uhf')
C
C  now reorder the occupied MOs, at the same time transpose them
C
        CALL ReorderMO(NBas,NAlpha,Z(i1),GX,GY)
        CALL ReorderMO(NBas,NBeta,Z(i1),GXB,GYB)
c -- form density matrix now if using weight derivatives
        Call FormDENU(NBas,NAlpha,NBeta,GY,GYB,DA,DB)
      ENDIF
C
C ================================================================
C  get location of abelian symmetry data in old texas depository
C
      call getival('ngener',ngen)         ! number of generators
      nupr = 1                            ! not set if symmetry not used
      ifp = 1                             ! ditto
      If(nsym.gt.0) Then
       call getival('SymNuPr1',nupr)      ! equivalent atoms
       call getival('SymFunPr',ifp)       ! equivalent basis functions
C
C  reorder the equivalent basis function array from Wolinski
C  special order to per atom order
C
       CALL ICpyVEC(7*NBas,block(ifp),Z(i2))
       CALL ReorderIFP(NBas,nsym,Z(i1),Z(i2),NEqBAS)
      EndIf
      call getival('nslv',nslv)           ! number of slaves (parallel)
C ================================================================
C
C  get number and location of symmetry-unique atoms under
C  highest abelian point group
C
      CALL GetUNQ(NAtoms, nsym,   XNuc,   block(nupr), NQ,
     $            IUNQ,   Z)
C
C ----------------------------------------------------------------
C
C  Now calculate the gradient
C
      IF(nslv.EQ.0) THEN
cc
c  SINGLE PROCESSOR MODE
c  ---------------------
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
        DO 100 IAtom=1,NQ
        If(IAtom.EQ.1) Then
          IEntry = 1
        Else If(IAtom.EQ.NQ) Then
          IEntry = -1
        EndIf
        ICntr = IUNQ(IAtom)
        If(rhf) Then
          CALL DFTGradC(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                  nsym,   ngen,   nsy, block(nupr), IPRNT,
     $                  lrad,   lang,   IradQ,  factor, NBatch,
     $                  DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                  XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                  BASDAT, INX,    BL,     ExpMIN, IdWt,
     $                  NBAtm,  DA,     NScr,   Z,      GX,
     $                  GY,     GZ,     GWT,    GC,     EL,
     $                  IEntry)
        Else
          CALL DFTGradU(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                  nsym,   ngen,   nsy, block(nupr), IPRNT,
     $                  lrad,   lang,   IradQ,  factor, NBatch,
     $                  DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                  XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                  BASDAT, INX,    BL,     ExpMIN, IdWt,
     $                  NBAtm,  DA,     DB,     NScr,   Z,
     $                  GX,     GY,     GZ,     GXB,    GYB,
     $                  GZB,    GWT,    GC,     EL,     IEntry)
        EndIf
 100    CONTINUE
cc
      ELSE
cc
c  PARALLEL MODE
c  -------------
c
        CALL Para_DFTG(dft,    NAtoms, XNuc,   IAN,    nsym,
     $                 ngen,   nsy, block(nupr), NEqBAS, NQ,
     $                 IUNQ,   IPRNT,  lrad,   lang,   IradQ,
     $                 factor, NBatch, DISTN,  AIJ,    RDIST,
     $                 XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                 NBas,   NShell, BASDAT, INX,    BL,
     $                 ExpMIN, NBAtm,  IdWt,   rhf,    DA,
     $                 DB,     NScr,   Z,      GX,     GY,
     $                 GZ,     GXB,    GYB,    GZB,    GWT,
     $                 GC,     EL)
cc
      ENDIF
C
C ------------------------------------------------------------------
C
C  if symmetry was used, correct EL
C
       EL = EL*DFloat(NSym+1)
C
C  ................................................................
C  If Weight Derivatives were NOT used, form full gradient
C  from DFT gradient matrices
C
      If(NSym.GT.0) Then
C
C  symmetrize partial DFT gradient matrices
C
        CALL GradSYMM(NGen,nsy,NBas,NEqBas,1,GX)
        CALL GradSYMM(NGen,nsy,NBas,NEqBas,2,GY)
        CALL GradSYMM(NGen,nsy,NBas,NEqBas,3,GZ)
        If(.NOT.rhf) Then
          CALL GradSYMM(NGen,nsy,NBas,NEqBas,1,GXB)
          CALL GradSYMM(NGen,nsy,NBas,NEqBas,2,GYB)
          CALL GradSYMM(NGen,nsy,NBas,NEqBas,3,GZB)
        EndIf
      EndIf
c
      IF(IdWt.EQ.0) THEN
C
C  contract DFT gradient matrices with density
C
        If(rhf) Then
          I1 = 0
          DO 70 IAtom = 1,NAtoms
          I1 = I1+1
          I2 = I1+NBAtm(IAtom)-1
          GradX = Zero
          GradY = Zero
          GradZ = Zero
          DO 69 J=1,NBas
          JJ = (J*(J-1))/2
          DO 68 I=I1,I2
          If(I.GT.J) Then
           IJ = (I*(I-1))/2 + J
          Else
           IJ = JJ + I
          EndIf
          GradX = GradX + DA(IJ)*GX(I,J)
          GradY = GradY + DA(IJ)*GY(I,J)
          GradZ = GradZ + DA(IJ)*GZ(I,J)
 68       CONTINUE
 69       CONTINUE
          GC(1,IAtom) = - GradX - GradX
          GC(2,IAtom) = - GradY - GradY
          GC(3,IAtom) = - GradZ - GradZ
          I1 = I2
 70       CONTINUE
        Else
          I1 = 0
          DO 90 IAtom = 1,NAtoms
          I1 = I1+1
          I2 = I1+NBAtm(IAtom)-1
          GradX = Zero
          GradY = Zero
          GradZ = Zero
          DO 89 J=1,NBas
          JJ = (J*(J-1))/2
          DO 88 I=I1,I2
          If(I.GT.J) Then
           IJ = (I*(I-1))/2 + J
          Else
           IJ = JJ + I
          EndIf
          GradX = GradX + DA(IJ)*GX(I,J) + DB(IJ)*GXB(I,J)
          GradY = GradY + DA(IJ)*GY(I,J) + DB(IJ)*GYB(I,J)
          GradZ = GradZ + DA(IJ)*GZ(I,J) + DB(IJ)*GZB(I,J)
 88       CONTINUE
 89       CONTINUE
          GC(1,IAtom) = - GradX - GradX
          GC(2,IAtom) = - GradY - GradY
          GC(3,IAtom) = - GradZ - GradZ
          I1 = I2
 90       CONTINUE
        EndIf
      ENDIF
c
      If(IPRNT.GT.1) Then
        call getival('iout',iout)
        write(iout,*) ' DFT gradient contribution is:'
        do i=1,natoms
        write(iout,*) gc(1,i),gc(2,i),gc(3,i)
        enddo
      EndIf
C  ----------------------------------------------------------------
C
C  write total gradient to <grad> file
C
      CALL RdGRAD(NAtoms,XNuc,'save')    ! read coulomb gradient
      CALL VScal(3*NAtoms,-One,GC)       ! force module writes FORCES!!
      CALL AddVEC(3*NAtoms,XNuc,GC,GC)   ! add coulomb + exchange (DFT)
      CALL WrGRAD(NAtoms,GC)             ! write total gradient to <grad>
c
      call secund(t2)
      call elapsec(elagrad2)
      t1 = (t2-t1)/60.0d0
      elaps = (elagrad2-elagrad1)/60.0d0
      WRITE(6,1100) t1,elaps
      write(6,*) 'Number of electrons over grid:   ',EL
c
      RETURN
c
 1000 FORMAT(/,' Calculating DFT Contribution to Gradient')
 1100 FORMAT('Master CPU time for XC part of gradient ='
     *,f8.2,' Elapsed = ',f8.2,' min'/)
c
      END
