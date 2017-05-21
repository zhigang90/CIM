      SUBROUTINE ATOM_Density(NAtoms, NBas,   NShell, NPrim,  NPseud,
     $                        NBeta,  AtSymb, NMem,   Z,      IErr)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Interface to ATOMSCF atomic program
C
C  Solves SCF for each unique atom producing spherically symmetric
C  atomic density. Total density is then a superposition of all the
C  atomic densities.  Result is a block-diagonal initial guess
C  density matrix. (Needs to be reordered into Wolinski order.)
C
C
C  ARGUMENTS
C
C  NAtoms  -  Total number of (real) atoms
C  NBas    -  Total number of contracted basis functions
C  NShell  -  Total number of shells
C  NPrim   -  Total number of primitive shells
C  NPseud  -  Total number of pseudopotentials (ECPs)
C  NBeta   -  number of beta spin electrons
C             (zero if closed-shell)
C  AtSymb  -  character storage for atomic symbols
C  NMem    -  scratch memory in double words
C  Z       -  scratch array
C  IErr    -  error flag
C              0 - OK    -1 - something went wrong
C
C
      DIMENSION Z(NMem)
      CHARACTER*8 AtSymb(NAtoms)
C
C
C  Get the memory
C  BE CAREFUL - due to the segmentation of L (SP) shells, for
C  small systems involving a heavy atom with hydrogen, there may
C  be more primitives for the heavy atom than the total number of
C  primitives for the system before segmentation. So allocate
C  extra storage for ATOMSCF arrays
C
      NScr = MAX(13*NShell,70+NBas**2)
      IMem = 8*NAtoms + 17*NPrim + 14*NShell + 9*NBas
     $         + NBas**2 + 12 + NScr
C
C  Check memory
C
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,8,'ATOM_Den')
      CALL ZeroIT(Z,IMem)               ! clear the memory
C
C
C  Allocate memory pointers
C
      IAN   = iptr                     !  atomic numbers
      IXC   = IAN  + NAtoms            !  geometry
      ICHG  = IXC  + 3*NAtoms          !  atomic charges
      IXM   = ICHG + NAtoms            !  atomic masses
      IBAS  = IXM  + NAtoms            !  basis function data
      INX   = IBAS + 13*NPrim          !  basis indexing array
      IPSP  = INX  + 12*NShell         !  pseudopotential array
      IBTM  = IPSP + NAtoms            !  no. basis functions per atom
      IODR  = IBTM + NAtoms            !  reordering array for Wolinski
      ISDR  = IODR + NBas              !  ordering array for S shells
      IPDR  = ISDR + NBas              !  ordering array for P shells
      IDCT  = IPDR + NBas              !  Spherical/Cartesian D shells
      IFCT  = IDCT + NShell            !  Spherical/Cartesian F shells
      JC    = IFCT + NShell            !  basis function data for ATOMSCF
      JZTA  = JC   + 6*NBas            !  basis exponents for ATOMSCF
      JCTR  = JZTA + NPrim*2           !  contraction coeffs for ATOMSCF
      JBAS  = JCTR + NPrim*2           !  basis data for ATOMSCF
      JBC   = JBAS + 6                 !    ditto
      IDEN  = JBC  + 6                 !  final density matrix
      IZ    = IDEN + NBas**2           !  general scratch storage
      IEnd  = IZ   + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,8,'ATOM_Den')
c
c -- set memory mark for TEXAS memory manager
      call getmem(IMem,lastx)
C
C
C  ----------------------------------------------------------------------
C
      CALL ATOMMAIN(NAtoms, NBas,   NShell, NPrim,  NPseud,
     $              NBeta,  AtSymb, Z(IAN), Z(IXC), Z(ICHG),
     $              Z(IXM), Z(IBAS),Z(INX), Z(IPSP),Z(IBTM),
     $              Z(IODR),Z(ISDR),Z(IPDR),Z(IDCT),Z(IFCT),
     $              Z(JC),  Z(JBAS),Z(JBC), Z(JZTA),Z(JCTR),
     $              Z(IDEN),NScr,   Z(IZ),  IErr)

C
C  ----------------------------------------------------------------------
C  free memory used in this routine
C
        call retmem(1)
c
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE ATOMMAIN(NAtoms, NBas,   NShell, NPrim,  NPseud,
     $                    NBeta,  AtSymb, IAN,    XC,     XCharg,
     $                    XMass,  BASDAT, INX,    IPSP,   NBAtm,
     $                    IORDER, SOrder, POrder, DCart,  FCart,
     $                    JC,     JBAS,   JBC,    ZETA,   CONTR,
     $                    D,      NScr,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ARGUMENTS
C
C  NAtoms  -  Total number of (real) atoms
C  NBas    -  Total number of contracted basis functions
C  NShell  -  Total number of shells
C  NPrim   -  Total number of primitive shells
C  NPseud  -  Total number of pseudopotentials (ECPs)
C  NBeta   -  number of beta spin electrons
C             (zero if closed-shell)
C  AtSymb  -  atomic symbols
C  IAN     -  array of atomic numbers
C  XC      -  Cartesian coordinates
C  XCharg  -  atomic charges
C  XMass   -  atomic masses    (not used)
C  BASDAT  -  basis function data (Texas format)
C  INX     -    ditto
C  IPSP    -  number of ECPs on each atomic center
C  NBAtm   -  number of basis functions per atom
C  IORDER  -  reordering array (atom order to Wolinski)
c  SOrder  -  array indicating order of s shells in final density matrix
c  POrder  -  array indicating order of p shells in final density matrix
c    These arrays are needed because L shells have to be split into
c    separate s and p for ATOMSCF and then recombined
c  DCart   -  integer array for spherical/Cartesian d shells
C  FCart   -  integer array for spherical/Cartesian f shells
C
C  The following are basis set data reformatted for ATOMSCF program
C
C  JC(6,*) -  contraction pattern. The first subscipt is the angular
C             momentum: 1=s,2=p,3=d,4=f etc. The second subscript
C             goes over the contractions, and shows the length of
C             each contraction
C  JBAS(6) -  number of primitive shells for each angular momentum
C  JBC(6)  -  number of contracted shells for each angular momentum value
C  ZETA    -  Gaussian orbital exponents for all primitive shells.
C             First the JBAS(1) s exponents, then the JBAS(2) p exponents...
C             Total number: JBAS(1)+JBAS(2)+..JBAS(6)=NPrim  (per atom)
C  CONTR   -  contraction coefficients in the same order
C
C  D       -  on exit contains block density matrix
C  NScr    -  dimension of scratch storage (double words)
C  Z       -  scratch storage 
C  IErr    -  error flag
C              0 - success;    -1 - error
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),XCharg(NAtoms),XMass(NAtoms),
     $          BASDAT(13,NPrim),INX(12,NShell),IPSP(NAtoms),
     $          NBAtm(NAtoms),IORDER(NBas)
      INTEGER SOrder(NBas),POrder(NBas),DCart(NShell),FCart(NShell)
      LOGICAL CartD,CartF
      DIMENSION JC(6,NBas),JBAS(6),JBC(6),ZETA(NPrim*2),CONTR(NPrim*2)
      DIMENSION D(NBas,NBas),Z(NScr)
      CHARACTER*8 AtSymb(NAtoms)
      Character*256 jobname,MOS,MOB
      Character*20 cdum
      Logical pflag
c
      Data IUnit/1/
c
      Common /job/jobname,lenJ
C
C
C  Read the Cartesian coordinates and atomic charges
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoordF(IUnit,  NAtoms, AtSymb, XC,     -1,
     $              jnk,    XCharg, XMass)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read the basis set data
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdbasis(IUnit,  NAtoms, AtSymb, XC,     Z(1),
     $             Z(1+NAtoms),BASDAT)
C  ...................................................................
C  Are there pseudopotentials?
C
      If(NPseud.GT.0) Then
        CALL IZeroIT(IPSP,NAtoms)
        call rdcntrl(IUnit,7,'$npseud',1,NPseud,rdum,cdum)
        Do I=1,NPseud
        READ(IUnit,*) NAtm,Izcore
        IPSP(NAtm) = Izcore
        EndDO
      EndIf
C  ...................................................................
C
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C ....................................................................
C  Be Careful with basis order!
C  We need basis ordered PER ATOM for ATOMSCF
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
C  need to reorder final Density matrix as will produce Density
C  ordered per atom and will need Wolinski special order
C
      i1 = 1
      i2 = i1 + NShell
      IEnd = i2 + 12*NShell - 1
      CALL MemCHK(NScr,IEnd,8,'ATOMMAIN')
c
      Call reorder(NShell,INX,IORDER,Z(i1),Z(i2),IAN)     ! Wolinski order
      Call SortBAS2(NShell,Z(i2),Z(i1),IORDER,INX)        ! per atom
cc      write(6,*) ' IORDER array is:'
cc      do i=1,nbas
cc      write(6,*) i,iorder(i)
cc      enddo
C
C 
C ....................................................................
C  Now call ATOMSCF
C
      IErr = -1
      pflag = .False.
      IOut = igetival('iout')
      IStart = 0
c
      WRITE(IOut,1000)
c
      DO 20 IAtm=1,NAtoms
cc      write(6,*) ' Solving for atom: ',IAtm
      zn = XCharg(IAtm)
      nb = NBAtm(IAtm)
      If(NPseud.GT.0) Then
        zeff = zn - IPSP(IAtm)
      Else
        zeff = zn
      EndIf
c
c -- first reorder atomic basis for ATOMSCF
      Call primbas(IAtm,   NShell, INX,    BASDAT, JC,
     $             JBAS,   JBC,    ZETA,   CONTR,  SOrder,
     $             POrder, DCart,  FCart,  ISymax)
ccccccccccccccccc
cc      write(6,*) ' Back from <primbas>'
cc      write(6,*) ' ISymax is:',isymax
cc      write(6,*)' JC array is:'
cc      do i=1,nbas
cc      write(6,1234) (jc(j,i),j=1,6)
cc 1234 format(2X,6I4)
cc      enddo
cc      write(6,*) ' JBAS  and  JBC arrays:'
cc      do i=1,6
cc      write(6,*) jbas(i),jbc(i)
cc      enddo
cc      write(6,*) ' ZETA  and  CONTR arrays:'
cc      do i=1,nprim
cc      write(6,*) zeta(i),contr(i)
cc      enddo
cc      write(6,*) ' SOrder array is:'
cc      do i=1,jbc(1)
cc      write(6,*) i,'  ',sorder(i)
cc      enddo
cc      write(6,*) ' POrder array is:'
cc      do i=1,jbc(2)
cc      write(6,*) i,'  ',porder(i)
cc      enddo
cc      write(6,*) ' DCart:',(DCart(i),i=1,JBC(3))
cc      write(6,*) ' FCart:',(FCart(i),i=1,JBC(4))
cccccccccccccccccc
c
c -- now solve atomic SCF
      Call ZeroIT(Z,nb**2)    ! why does this make a difference?  JB
      Call atomscf(pflag, IOut,   zeff,   zn,     ISymax,
     $             JC,    JBAS,   JBC,    ZETA,   CONTR,
     $             Z)
C
C  On return from <atomscf> atomic density is in Z and is given
C  for only ONE angular momentum component. Thus s-orbitals are
C  fine, p-orbitals need to be divided by 3 and duplicated 3 times,
C  d-orbitals need to be divided by 5 and duplicated 5 times etc...
C
      IJ = 0
      nsi = 0
      npi = 0
      DO IAng=1,ISymax
      IShell = 2*IAng-1     ! shell size (1,3,5,7)
      IDim = JBC(IAng)      ! no. of shells per angular momenta
c
      If(IAng.EQ.1) Then     ! S-orbitals
        IS = IStart+1
        IE = IS+IDim-1
        DO 10 I=IS,IE
        nsi = nsi+1
        II = SOrder(nsi) + IStart
        nsj = 0
        DO 10 J=IS,I
        nsj = nsj+1
        JJ = SOrder(nsj) + IStart
        IJ = IJ+1
        D(II,JJ) = Z(IJ)
cc        write(6,*) ' S-loop  I:',ii,' J:',jj,z(ij)
 10     CONTINUE
      Else If(IAng.EQ.2) Then     ! P-orbitals
        IS = IStart+JBC(1)+1
        IE = IS+IDim-1
        DO 11 I=IS,IS+3*(IE-IS),3
        npi = npi+1
        II = POrder(npi) + IStart
        npj = 0
        DO 11 J=IS,I,3
        npj = npj+1
        JJ = POrder(npj) + IStart
        IJ = IJ+1
        D(II,JJ) = Z(IJ)/3.0d0
        D(II+1,JJ+1) = D(II,JJ)
        D(II+2,JJ+2) = D(II,JJ)
cc        write(6,*) ' P-loop  I:',ii,' J:',jj,z(ij)/3.0d0
cc        write(6,*) ' P-loop  I:',ii+1,' J:',jj+1
cc        write(6,*) ' P-loop  I:',ii+2,' J:',jj+2
 11     CONTINUE
      Else If(IAng.EQ.3) Then     ! D-orbitals
        IS = IStart+JBC(1)+3*JBC(2)+1
        IE = IS+IDim-1
        DO 12 I=IS,IS+5*(IE-IS),5
        DO 12 J=IS,I,5
        IJ = IJ+1
        D(I,J) = Z(IJ)/5.0d0
        D(I+1,J+1) = D(I,J)
        D(I+2,J+2) = D(I,J)
        D(I+3,J+3) = D(I,J)
        D(I+4,J+4) = D(I,J)
cc        write(6,*) ' D-loop  I:',i,' J:',j,z(ij)/5.0d0
cc        write(6,*) ' D-loop  I:',i+1,' J:',j+1
cc        write(6,*) ' D-loop  I:',i+2,' J:',j+2
cc        write(6,*) ' D-loop  I:',i+3,' J:',j+3
cc        write(6,*) ' D-loop  I:',i+4,' J:',j+4
 12     CONTINUE
      Else If(IAng.EQ.4) Then     ! F-orbitals
        IS = IStart+JBC(1)+3*JBC(2)+5*JBC(3)+1
        IE = IS+IDim-1
        DO 13 I=IS,IS+7*(IE-IS),7
        DO 13 J=IS,I,7
        IJ = IJ+1
        D(I,J) = Z(IJ)/7.0d0
        D(I+1,J+1) = D(I,J)
        D(I+2,J+2) = D(I,J)
        D(I+3,J+3) = D(I,J)
        D(I+4,J+4) = D(I,J)
        D(I+5,J+5) = D(I,J)
        D(I+6,J+6) = D(I,J)
cc        write(6,*) ' F-loop  I:',i,' J:',j,z(ij)/7.0d0
cc        write(6,*) ' F-loop  I:',i+1,' J:',j+1
cc        write(6,*) ' F-loop  I:',i+2,' J:',j+2
cc        write(6,*) ' F-loop  I:',i+3,' J:',j+3
cc        write(6,*) ' F-loop  I:',i+4,' J:',j+4
cc        write(6,*) ' F-loop  I:',i+5,' J:',j+5
cc        write(6,*) ' F-loop  I:',i+6,' J:',j+6
 13     CONTINUE
      EndIF
      EndDO
C
C  Density is given in terms of Spherical (pure) D and F
C  Basis may have Cartesian D/F
C  If so, need to transform
C
      nd = JBC(3)
      nf = JBC(4)
      CartD = .False.
      CartF = .False.
      If(nd.GT.0) Then
        DO I=1,nd
        If(DCart(I).EQ.6) Then
          CartD = .True.
          Exit
        EndIf
        EndDO
cc        write(6,*) (DCart(I),I=1,nd)
      EndIf
      If(nf.GT.0) Then
        DO I=1,nf
        If(FCart(I).EQ.10) Then
          CartF = .True.
          Exit
        EndIf
        EndDO
cc        write(6,*) (FCart(I),I=1,nf)
      EndIf
c
      If(CartD) Then
c -- locate and transform D5 blocks to D6
        ids = IStart + JBC(1) + JBC(2)*3     ! start of D block
        i1 = 1
        i2 = i1 + NBas**2
        CALL D5TOD6(NBas,nd,DCart,ids,nf,Z(i1),Z(i2),D)
      EndIf
c
      If(CartF) Then
c -- locate and transform F7 blocks to F10
        ifs = IStart + JBC(1) + JBC(2)*3 + JBC(3)*5
        DO I=1,nd
        If(DCart(I).EQ.6) ifs = ifs+1
        EndDO
        CALL F7TOF10(NBas,nf,FCart,ifs,Z(i1),Z(i2),D)
      EndIf
cc
      IStart = IStart+nb
 20   CONTINUE
cc      write(6,*) ' Density Matrix is:'
cc      call prntmat(NBas,NBas,NBas,D)
C
C  now reorder density to Wolinski
C
      IJ = 0
      DO 40 I=1,NBas
      II = IORDER(I)
      DO 30 J=1,I
      JJ = IORDER(J)
      IJ = IJ+1
      If(II.GE.JJ) Then
        Z(IJ) = D(II,JJ)
      Else
        Z(IJ) = D(JJ,II)
      EndIf
 30   CONTINUE
 40   CONTINUE
C
      CALL EXPAND(NBas,Z,D)
cc      write(6,*) ' Wolinski Density Matrix is:'
cc      call prntmat(NBas,NBas,NBas,D)
C
C  We are going to save the density matrix on the MOs file
C  set up the MOs filename
C
      MOS=jobname(1:lenJ)//'.mos'
      lenM=lenJ+4
      MOB=MOS
      MOB(lenM:lenM)='b'
      call setchval('mos-file',MOS)
      call setchval('mob-file',MOB)
c
      If(NBeta.GT.0) Then
        CALL VScal(NBas**2, 0.5d0, D)
        Call WriteMOS(NBas,   NBas,   D,      jnk,    .False.,
     $                lenM,   MOB,    3)
      EndIf
c
      Call WriteMOS(NBas,   NBas,   D,      jnk,    .False.,
     $              lenM,   MOS,    3)
c
      IErr = 0
c
      RETURN
c
 1000 FORMAT(' Superposition of Atomic Densities Guess')
c
      END
c  =======================================================================
c
      SUBROUTINE D5TOD6(NBas,   nd,     DCart,  ids,    nf,
     $                  DT,     Z,      D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transforms blocks of atomic density matrix from D5 to D6
C  (Spherical to Cartesian)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  nd      -  number of D shells
C  DCart   -  whether D shell is D5 or D6
C  ids     -  starting address of D block
C             note: may be changed on exit
C  nf      -  number of F shells (if any)
C  DT      -  intermediate Density storage
C  Z       -  scratch storage
C  D       -  on entry - density matrix with D5 blocks
C             on exit  - density matrix with D6 blocks
C
C
      DIMENSION DT(NBas,NBas),D(NBas,NBas),Z(6,5)
      INTEGER DCart(nd)
c
c -- transformation matrix
      Dimension DTRANS(6,5)
C
C  set up DTRANS
C
      CALL ZeroIT(DTRANS,30)
      DTRANS(1,1) = -1.0d0/SQRT(12.0d0)
      DTRANS(2,1) = -1.0d0/SQRT(12.0d0)
      DTRANS(3,1) =  2.0d0/SQRT(12.0d0)
      DTRANS(1,2) =  0.5d0
      DTRANS(2,2) = -0.5d0
      DTRANS(4,3) =  1.0d0
      DTRANS(5,4) =  1.0d0
      DTRANS(6,5) =  1.0d0
C
C  copy all S and P blocks into DT
      CALL ZeroIT(DT,NBas**2)
      DO 10 J=1,ids
      DO 10 I=1,ids
      DT(I,J) = D(I,J)
 10   CONTINUE
c
      ids0 = ids      ! save original pointer
      id5 = ids
      id6 = ids
      DO L=1,nd
      jd5 = ids0
      jd6 = ids0
      DO M=1,L
c
      If(DCart(L).EQ.6.AND.DCart(M).EQ.6) Then
c -- transform columns
        DO I=1,6
        DO J=1,5
        JJ = J+jd5
        Temp = 0.0d0
        DO K=1,5
        Temp = Temp + DTRANS(I,K)*D(K+id5,JJ)
        EndDO
        Z(I,J) = Temp
        EndDO
        EndDO
c -- transform rows
        DO I=1,6
        II = I+id6
        DO J=1,6
        JJ = J+jd6
        Temp = 0.0d0
        DO K=1,5
        Temp = Temp + Z(I,K)*DTRANS(J,K)
        EndDO
        DT(II,JJ) = Temp
        EndDO
        EndDO
cc
      Else If(DCart(L).EQ.6) Then
c -- transform columns
        DO I=1,6
        II = I+id6
        DO J=1,5
        JJ = J+jd5
        Temp = 0.0d0
        DO K=1,5
        Temp = Temp + DTRANS(I,K)*D(K+id5,JJ)
        EndDO
        DT(II,JJ) = Temp
        EndDO
        EndDO
cc
      Else If(DCart(M).EQ.6) Then
c -- transform rows
        DO I=1,5
        II = I+id5
        DO J=1,6
        JJ = J+jd6
        Temp = 0.0d0
        DO K=1,5
        Temp = Temp + D(II,K+jd5)*DTRANS(J,K)
        EndDO
        DT(II,JJ) = Temp
        EndDO
        EndDO
      EndIf
c
      jd5 = jd5 + 5
      jd6 = jd6 + 6
      EndDO
      id5 = id5 + 5
      id6 = id6 + 6
      EndDO
c
      If(nf.GT.0) Then
c -- copy all (untransformed) F blocks into DT
        DO I=1,nd
        ids = ids0 + DCart(I)
        EndDO
        ids0 = ids0 + 5*nd
        DO 20 J=1,7*nf
        JJ = J+ids
        JJ0 = J+ids0
        DO 20 I=1,7*nf
        DT(I+ids,JJ) = D(I+ids0,JJ0)
 20     CONTINUE
      EndIF
C
C  now copy DT back into D
C
      CALL CpyVEC(NBas**2,DT,D)
c
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE F7TOF10(NBas,   nf,     FCart,  ifs,    DT,
     $                   Z,      D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transforms blocks of atomic density matrix from F7 to F10
C  (Spherical to Cartesian)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  nf      -  number of F shells
C  FCart   -  whether F shell is F7 or F10
C  ifs     -  starting address of F block
C             note: may be changed on exit
C  DT      -  intermediate Density storage
C  Z       -  scratch storage
C  D       -  on entry - density matrix with F7 blocks
C             on exit  - density matrix with F10 blocks
C
C
      DIMENSION DT(NBas,NBas),D(NBas,NBas),Z(10,7)
      INTEGER FCart(nf)
c
c -- transformation matrix
      Dimension FTRANS(10,7)
C
C  set up FTRANS
C
      CALL ZeroIT(FTRANS,70)
      FTRANS(2,1)  =  4.0d0
      FTRANS(7,1)  = -1.0d0
      FTRANS(9,1)  = -1.0d0
      FTRANS(3,2)  =  4.0d0
      FTRANS(8,2)  = -1.0d0
      FTRANS(10,2) = -1.0d0
      FTRANS(1,3)  = -1.0d0
      FTRANS(4,3)  =  4.0d0
      FTRANS(6,3)  = -1.0d0
      FTRANS(3,4)  = -1.0d0
      FTRANS(8,4)  =  4.0d0
      FTRANS(10,4) = -1.0d0
      FTRANS(1,5)  = -1.0d0
      FTRANS(4,5)  = -1.0d0
      FTRANS(6,5)  =  4.0d0
      FTRANS(2,6)  = -1.0d0
      FTRANS(7,6)  = -1.0d0
      FTRANS(9,6)  =  4.0d0
      FTRANS(5,7)  =  1.0d0
C
C  copy all S and P and D blocks into DT
      CALL ZeroIT(DT,NBas**2)
      DO 10 J=1,ifs
      DO 10 I=1,ifs
      DT(I,J) = D(I,J)
 10   CONTINUE
c
      ifs0 = ifs      ! save original pointer
      if7 = ifs
      ift = ifs
      DO L=1,nf
      jf7 = ifs0
      jft = ifs0
      DO M=1,L
c
      If(FCart(L).EQ.10.AND.FCart(M).EQ.10) Then
c -- transform columns
        DO I=1,10
        DO J=1,7
        JJ = J+jf7
        Temp = 0.0d0
        DO K=1,7
        Temp = Temp + FTRANS(I,K)*D(K+if7,JJ)
        EndDO
        Z(I,J) = Temp
        EndDO
        EndDO
c -- transform rows
        DO I=1,10
        II = I+ift
        DO J=1,10
        JJ = J+jft
        Temp = 0.0d0
        DO K=1,7
        Temp = Temp + Z(I,K)*FTRANS(J,K)
        EndDO
        DT(II,JJ) = Temp
        EndDO
        EndDO
cc
      Else If(FCart(L).EQ.10) Then
c -- transform columns
        DO I=1,10
        II = I+ift
        DO J=1,7
        JJ = J+jf7
        Temp = 0.0d0
        DO K=1,7
        Temp = Temp + FTRANS(I,K)*D(K+if7,JJ)
        EndDO
        DT(II,JJ) = Temp
        EndDO
        EndDO
cc
      Else If(FCart(M).EQ.10) Then
c -- transform rows
        DO I=1,7
        II = I+if7
        DO J=1,10
        JJ = J+jft
        Temp = 0.0d0
        DO K=1,7
        Temp = Temp + D(II,K+jf7)*FTRANS(J,K)
        EndDO
        DT(II,JJ) = Temp
        EndDO
        EndDO
      EndIf
c
      jf7 = jf7 + 7
      jft = jft + 10
      EndDO
      if7 = if7 + 7
      ift = ift + 10
      EndDO
C
C  now copy DT back into D
C
      CALL CpyVEC(NBas**2,DT,D)
c
      RETURN
      END
c  =======================================================================
c
      subroutine primbas(iatom,  ncs,    inx,    basdat, ic,
     $                   nbas,   nbc,    zeta,   contr,  SOrder,
     $                   POrder, DCart,  FCart,  nsymax)
      implicit real*8 (a-h,o-z)
c
c  ARGUMENTS
c  INPUT
c  iatom: the atom, for which the basis is retrieved
c  ncs: number of contracted shells in the molecule
c  inx: contraction information, not identical with common inx
c  basdat(13,*): basis set data
c  OUTPUT
c  ic(6,*): contraction info for ATOMSCF, lengths of contractions
c    per angular momentum (s,p,d,f,g,h)
c  nbas(6): number of primitive shells for each angular momentum
c  nbc(6): number of contracted shells for each angular momentum
c  zeta(*): orbital exponents, for ALL types together, first the s
c    exponents, then the p ones etc.
c  contr(*): contraction coefficients, in the same order as the
c    exponents
c  SOrder: array indicating order of s shells in final density matrix
c  POrder: array indicating order of p shells in final density matrix
c    These arrays are needed because l shells have to be split into
c    separate s and p for ATOMSCF and then recombined
c    IMPORTANT - assumed that ALL s shells and ALL p shells on a
c                given atom are dealt with BEFORE any l shells
c  DCart: integer array for spherical/Cartesian d shells
C  FCart: integer array for spherical/Cartesian f shells
c  nsymax: maximum type of basis function which occurs
c
      dimension ic(6,*),nbas(6),nbc(6),zeta(*),contr(*)
      dimension inx(12,*),basdat(13,*)
      dimension iptr(6)
      integer SOrder(*),POrder(*),DCart(*),FCart(*)
c
      nsymax=1     ! changed - was originally 0!!   JB
      do i=1,6
        nbc(i)=0
        nbas(i)=0
        iptr(i)=0
      end do
      nss=0
      npp=0
      nd=0
      nf=0
c
      DO 10 ics=1,ncs
c
      If(inx(2,ics).eq.iatom) Then
        itype=inx(12,ics)
c
c -- sort out itype1, s and p ordering arrays, and Cartesian D/F
        If(itype.EQ.1) Then
          itype1 = itype
          nss = nss+1
          SOrder(nss) = nss
        Else If(itype.EQ.2) Then
          itype1 = itype
          npp = npp+1
          if(npp.eq.1) POrder(npp) = nss+1
          if(npp.gt.1) POrder(npp) = POrder(npp-1)+3
        Else If(itype.EQ.3) Then
c -- L shell here
          itype1 = 2
          nss = nss+1
          npp = npp+1
          If(npp.eq.1) SOrder(nss) = SOrder(nss-1)+1
          If(npp.gt.1) SOrder(nss) = POrder(npp-1)+3
          POrder(npp) = SOrder(nss)+1
        Else If(itype.eq.4) Then
          nd = nd+1
          itype1 = 3
          DCart(nd) = 5
        Else If(itype.eq.5) Then
          nd = nd+1
          itype1 = 3
          DCart(nd) = 6
        Else If(itype.eq.6) Then
          nf = nf+1
          itype1 = 4
          FCart(nf) = 7
        Else If(itype.eq.7) Then
          nf = nf+1
          itype1 = 4
          FCart(nf) = 10
        Else If(itype.gt.7) Then
c -- greater than f shell, simply ignore
          GO TO 10
        EndIf
        if(itype.lt.1) then
          call nerror(1,'primbas','impossible type',itype,0)
        end if
c
        nbcit=nbc(itype1)+1
        nbc(itype1)=nbcit
        nprim=inx(5,ics)-inx(1,ics)
cc        write(6,*) ' itype1:',itype1,' nprim:',nprim
        nbas(itype1)=nbas(itype1)+nprim
        ic(itype1,nbcit)=nprim
        if(itype.eq.3) then
c  L type shell
          nbc(1)=nbc(1)+1
          nbas(1)=nbas(1)+nprim
          ic(1,nbc(1))=nprim
        end if
      EndIf
 10   CONTINUE
c
      iptr(1)=0
      do ii=2,6
        iptr(ii)=iptr(ii-1)+nbas(ii-1)
c  maximum angular momentum  1=s, 2=p etc./        
        if(nbas(ii).gt.0) nsymax=ii
      end do
c
      DO 20 ics=1,ncs
      If(inx(2,ics).eq.iatom) Then
        itype=inx(12,ics)
        itype1=itype
        if(itype.eq.3) itype1=2
        if(itype.eq.4) itype1=3
        if(itype.eq.5) itype1=3
        if(itype.eq.6) itype1=4
        if(itype.eq.7) itype1=4
c -- greater than f shell, simply ignore
        if(itype.gt.7) GO TO 20
        nprim=inx(5,ics)-inx(1,ics)
cc        write(6,*) ' shell=',ics,' itype1:',itype1,' # prim=',nprim
c  deal with the S component of the L shell first
        if(itype.eq.3) itype1=1
        do k=inx(1,ics)+1,inx(5,ics)
          ipt=iptr(itype1)+1
          iptr(itype1)=ipt
          zeta(ipt)=basdat(1,k)
          factor=1.0d0
          if(itype1.eq.2) factor=0.5d0/sqrt(zeta(ipt))
          if(itype1.eq.3) factor=0.25d0/zeta(ipt)
          if(itype1.eq.4) factor=0.125d0/zeta(ipt)*sqrt(zeta(ipt))
cc          if(itype1.eq.5) factor=0.0625d0/zeta(ipt)**2
cc          write(6,*) ' primitive=',k,' exponent=',zeta(ipt)
          contr(ipt)=basdat(2,k)*factor
          if(itype.eq.3) then
c  L type shell
            ipt=iptr(2)+1
            iptr(2)=ipt
            zeta(ipt)=basdat(1,k)
            contr(ipt)=basdat(3,k)*0.5d0/sqrt(zeta(ipt))
          end if
        end do
      EndIf
 20   CONTINUE
c
      return
      end
c  =======================================================================
c
      subroutine atomscf(oprin,  iwr,    znps,   zn_in,  isymax,
     1                   ic_in,  nbas_in,nbc_in, zeta_in,cont_in,
     2                   dt)
      implicit real*8 (a-h,o-z), integer (i-n)
c
c...  calling routine for atomscf,
c...  jvl 1998
c  ARGUMENTS
c  INPUT
c  oprin(logical): print flag
c  iwr: unit number for printing
c  znps: effective nuclear charge
c  zn_in: actual nuclear charge
c  isymax: highest occupied angular momentum value (s=1,p=2,...)
c  ic_in: contraction pattern, see the description in atombasis
c  nbas_in: number of primitive shells/angular momentum
c  nbc_in:  number of contracted shells/angular momentum
c  zeta_in: Gaussian exponents
c  cont_in: contraction coefficients.
c     See atombasis for the latter two
c  OUTPUT
c  dt: density matrix  
c
      logical oprin
      real*8 zeta_in(*),cont_in(*)
      integer ic_in(6,*),nbc_in(*),nbas_in(*)
      real*8 dt(*)
c
      parameter (nbig=100,no=50)
      parameter (ntr=nbig*(nbig+1)/2)
      parameter (nsq=nbig*nbig)
      parameter (nno=nbig*no)
c
      integer ic(6,nbig)
      real*8 pcap(ntr),qcap(ntr),fc(ntr),fo(ntr),s(ntr),
     1 u(ntr),t(ntr),h(ntr),dc(ntr),dos(ntr),dold(ntr),ss(ntr),
     1 cvec(nsq), copn(nsq), smin(nno), qmin(nno), transf(nsq), cc(nsq)
      real*8 hatom(ntr)
c
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
c
      zn = zn_in
      nn = 0
      do i=1,isymax
       nbc(i) = nbc_in(i)
       nbas(i) = nbas_in(i)
       nn = nn + nbas(i)
       do j=1,nbc(i)
        ic(i,j) = ic_in(i,j)
       end do
      end do
      do i=1,nn
       zeta(i) = zeta_in(i)
       cont(i) = cont_in(i)
      end do
c
      call  atomd(oprin,iwr,znps,ic,isymax,hatom,
     + pcap,qcap,fc,fo,s,u,t,h,dc,dos,dt,dold,ss,
     + cvec, copn, smin, qmin, transf, cc , nbig)
c    
      return
      end
c  =======================================================================
c
      subroutine atomd(oprin,iwr,znps,ic,isymax,hatom,
     + pcap,qcap,fc,fo,s,u,t,h,dc,dos,dt,dold,ss,
     + cvec, copn, smin, qmin, transf, cc , nbb)
c
      implicit real*8  (a-h,o-z),integer   (i-n)
      parameter (nbig=100, no=50)
      dimension ic(6,*),hatom(*)
      dimension pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      dimension h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      dimension cvec(*),copn(*),smin(nbb,*),qmin(nbb,*),transf(*),cc(*)
c
      logical oprin
c.......................................................................
c     atomic r h f code for gto's. uses roothaan double diagonalization.
c.......................................................................
c
c      routine datoms .. tramad are a totally separate unit
c      they use commom junk to communicate
c      oprin,iwr are new parameters to control printing
c      hatom are the pseudo-corrections in the contracted basis from gam
c            they replace atom's own integrals
c.......................................................................
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
c
c
c     nsht      = total number of shells
c     n1(i)     = nbas(i) * (nbas(i) + 1 ) / 2
c     nbc(i)    = number of cont. orbitals in symmetry i
c     cont(i)   = contraction coeff. assosiated with primitive no. i
c     nstrt(i)  = number for first primitive in cont. no. i
c     nbct      = total number of cont. basis functions.
c.......................................................................
c     zn = znx
      if (dabs(zn).ge.1.d-8) then
c.......................................................................
c
c     distribute electrons according to aufbau
c.......................................................................
         call atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps)
crz
c     fix for use of pseudopotentials
         znsave = zn
         zn = znps
c.......................................................................
c
c     move basis set information from transfer variables to working
c     variables. ** not necessary in present version
c.......................................................................
         nbct = 0
         ndim = 0
         nsqt = 0
         nbcdim = 0
         jcount = 0
         nsht = 0
         k = 0
         do 40 i = 1 , nsym
            n1(i) = nbas(i)*(nbas(i)+1)/2
            nbcdim = nbcdim + nbc(i)*(nbc(i)+1)/2
            ndim = ndim + n1(i)
            nsht = nsht + ncsh(i) + nosh(i)
            nsqt = nsqt + nbas(i)**2
            nbct = nbct + nbc(i)
            do 20 j = 1 , nbas(i)
               jcount = jcount + 1
               nqn(jcount) = i
 20         continue
            do 30 j = 1 , n1(i)
               k = k + 1
               dold(k) = 0.0d0
 30         continue
 40      continue
         ns = 1
         nstrt(ns) = 1
         do 60 l = 1 , nsym
            iant = nbc(l)
            ksum = 0
            do 50 i = 1 , iant
               nstrt(ns+i) = nstrt(ns+i-1) + ic(l,i)
               ksum = ksum + ic(l,i)
 50         continue
            if (ksum.ne.nbas(l)) write (iwr,6020) l , ksum , nbas(l)
            ns = ns + iant
 60      continue
         do 70 i = 1 , nsqt
            cvec(i) = 0.d0
            copn(i) = 0.0d0
            cc(i) = 0.d0
 70      continue
         do 80 i = 1 , ndim
            dold(i) = 0.0d0
 80      continue
cjvl  few extra checks
         do 81 i=1,nsym
            if (nbc(i).lt.ncsh(i)+nosh(i)) then
               oprin = .true.
               ndim = 0
            end if
81       continue
cjvl
         if (oprin) write (iwr,6030) zn
         if (oprin) write (iwr,6040) (nbas(i),i=1,nsym)
         if (oprin) write (iwr,6050) (nbc(i),i=1,nsym)
         if (oprin) write (iwr,6060) (ncsh(i),i=1,nsym)
         if (oprin) write (iwr,6070) (nosh(i),i=1,nsym)
         if (oprin) write (iwr,6080) (nccup(i),i=1,nsym)
c        if (oprin) write(iwr,1700)(2*ajmn(i),i=1,20),(ajmn(i),i=21,24)
cjvl
         if (ndim.eq.0) call nerror(1,'atomd','ndim is 0',ndim,0)
cjvl
         maxitr = 100
c..
c..     calculate 1-electron ints
c..
         call oeigd(fc,s,u,t,h)
c.......................................................................
c
c     copy overlap matrix to ss
c.......................................................................
         do 90 i = 1 , ndim
            ss(i) = s(i)
            fc(i) = s(i)
 90      continue
c.......................................................................
c
c     now transform ss to contracted basis, then set up transformation
c     matrix to o.n. contracted basis.
c.......................................................................
         call trafsd(nsym,nbas,ndim,ss,nbc,cont,nstrt,dc)
         call trafsd(nsym,nbas,ndim,fc,nbc,cont,nstrt,dc)
c...
         nstep1 = 1
         nstep2 = 1
         do 100 i = 1 , nsym
crz      to surpress problems when there is one type of shell missing...
            if (nbc(i).ne.0) then
               call shalfd(fc(nstep1),transf(nstep2),nbc(i),nbb)
               call starcd(cc(nstep2),ss(nstep1),nbc(i),ncsh(i),
     +                     nosh(i))
               nstep1 = nstep1 + nbc(i)*(nbc(i)+1)/2
               nstep2 = nstep2 + nbc(i)**2
            end if
 100     continue
         nitscf = 0
         nconv = 0
         damp = .30d0
         if (znsave.eq.30.0d0) then
            damp = .9d0
            maxitr = 200
c         print *,' zn is special ',damp,maxitr
         end if
 110     nitscf = nitscf + 1
c.......................................................................
c
c     transform vectors and set up matrices in primitive basis,
c     then transform fock matrices to contracted basis.
c.......................................................................
         call tracd(cvec,cc,nsqt)
         call densid(dt,dold,dos,nsym,nosh,ncsh,nccup,cvec,damp,
     +               nconv,nbas,nitscf,tlarge)
c
c... check for convergence on tlarge (max change of d-matrix)
c
         if (tlarge.le.1.0d-5) nconv = 1
         if (nitscf+20.ge.maxitr) then
            write (iwr,6010) nitscf , energ , tlarge
         end if
c
         call hamild(pcap,qcap,fc,fo,s,u,t,h,dos,dt,cvec,smin,qmin,nbb)
         call trafsd(nsym,nbas,ndim,fc,nbc,cont,nstrt,dc)
         call trafsd(nsym,nbas,ndim,fo,nbc,cont,nstrt,dc)
c...
c...    now add the h-pseudo-contributions (from xpsnld)
c...
         do 120 i = 1 , nbcdim
            fc(i) = fc(i) + hatom(i)
            fo(i) = fo(i) + hatom(i)
 120     continue
c.......................................................................
c
c     do double diagonalization by symmetries:
c         1. transform block to o.n.basis (contracted).
c         2. store o.n. transformation matrix in vector matrix.
c         3. diagonalize.
c         4. order eigenvectors by eigenvalue.
c         5. if necessary, merge open and closed vectors.
c.......................................................................
         nstep1 = 1
         nstep2 = 1
         knteps = 0
         do 160 i = 1 , nsym
            nbc1 = nbc(i)
            nbc2 = nbc1**2
            nbc3 = (nbc2+nbc1)/2
            if (ncsh(i).ne.0) then
               call tramad(fc(nstep1),transf(nstep2),dc,nbc3,nbc1,dt)
               call dcopy(nbc2,transf(nstep2),1,cc(nstep2),1)
               call jacod(fc(nstep1),cc(nstep2),nbc1,n1(i),nbc2,1,nbc1,
     +                    dc,dt,nbc1)
            end if
            if (nosh(i).ne.0) then
               call tramad(fo(nstep1),transf(nstep2),dc,nbc3,nbc1,dt)
               call dcopy(nbc2,transf(nstep2),1,copn(nstep2),1)
               call jacod(fo(nstep1),copn(nstep2),nbc1,n1(i),nbc2,1,
     +                    nbc1,dc,dt,nbc1)
            end if
            icount = nstep1
            do 130 j = 1 , nbc1
               dc(j) = fc(icount)
               dos(j) = fo(icount)
               icount = icount + 1 + j
 130        continue
            call orderd(cc(nstep2),nbc1,nbc1,idum,idum,dt,dc,100)
            if (nosh(i).gt.0) then
               call orderd(copn(nstep2),nbc1,nbc1,idum,idum,dt,dos,100)
               call cmergd(cc(nstep2),copn(nstep2),ncsh(i),nbc1,nosh(i))
            end if
            nstep1 = nstep1 + nbc1*(nbc1+1)/2
            nstep2 = nstep2 + nbc2
            if (nconv.gt.0) then
               if (ncsh(i).gt.0) then
                  do 140 j = 1 , ncsh(i)
                     knteps = knteps + 1
                     eps(knteps) = dc(j)
 140              continue
               end if
               if (nosh(i).gt.0) then
                  do 150 j = 1 , nosh(i)
                     knteps = knteps + 1
                     eps(knteps) = dos(ncsh(i)+j)
 150              continue
               end if
            end if
 160     continue
         if (nitscf.ge.maxitr) nconv = 1
         if (nconv.le.0) go to 110
         if (oprin) then
            call outpud(copn,cc,1,iwr)
         else
            call outpud(copn,cc,0,iwr)
         end if
         call densid(dt,dold,dos,nsym,nosh,ncsh,nccup,cc,damp,nconv,nbc,
     +               nitscf,tlarge)
c
c...     correct energy for pseudopential or zora
c
         energ = energ + ddot(nbcdim,hatom,1,dt,1)
c        
         lm = 0
         do 190 i = 1 , nsym
            noff = lm
            nbci = nbc(i)
            do 180 l = 1 , nbci
               ll = noff + l*(l+1)/2
               do 170 m = 1 , l
                  mm = noff + m*(m+1)/2
                  lm = lm + 1
                  dt(lm) = dt(lm)*dsqrt(ss(ll)*ss(mm))
                  if (m.ne.l) dt(lm) = dt(lm)/2.0d0
 170           continue
 180        continue
 190     continue
      else
c.......................................................................
c
c     special section for handling the case of floating functions
c     on centers with no charge.
c.......................................................................
         nbcdim = 0
         do 200 i = 1 , isymax
            nbcdim = nbcdim + nbc(i)*(nbc(i)+1)/2
 200     continue
         do 210 i = 1 , nbcdim
            dt(i) = 0.0d0
 210     continue
         nsym = isymax
         energ = 0.0d0
      end if
c
      return
 6010 format (' it.',i4,'  energy',d19.10,'  div.',d13.5)
 6020 format ('-',' wrong contraction in symmetry   ',3i5)
 6030 format (/'      charge =',f10.6,//'      symmetry species',12x,
     +        's',5x,'p',5x,'d',5x,'f')
 6040 format (6x,'number of basis functions =',4(i2,4x))
 6050 format (6x,'number of cont. functions =',4(i2,4x))
 6060 format (6x,'number of closed shells   =',4(i2,4x))
 6070 format (6x,'number of open shells     =',4(i2,4x))
 6080 format (6x,'open shell occupation     =',4(i2,4x))
c1700 format('0     vector coupling coefficients',/(6f16.8))
      end
c  =======================================================================
c
      subroutine atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps)
      implicit real*8  (a-h,o-z),integer (i-n)
c.......................................................................
c
c     for atom of nuclear charge zn find electron configuration
c     and k(l,l,0) coupling coefficient. simple aufbau is
c     assumed throughout the periodic system. the algorithm used
c     will work for zn less than 119.
c
c     to allow for pseudo-potentials the effective charge znps
c     is used to determine the number of closed shells filled
c     jvl  (daresbury 1988): corrected July 93
c.......................................................................
      dimension ncsh(*),nccup(*),nosh(*),ajmn(*)
      dimension nelhw(103),nelcep(103)
      dimension isymhw(13)
      dimension isymce(13)
      data isymhw / 1,1,2,1,2, 3,1,2, 3,1,2, 4, 3/
      data isymce / 1,1,2,1,2, 3,1,2, 3, 4,1,2, 3/
      data nelhw/
     $         0, 0,
     $         2, 2, 2, 2, 2, 2, 2, 2,
     $         10, 10, 10, 10, 10, 10, 10, 10,
     $         18, 18,
     $                     18, 18, 18, 18, 18,
     $                     18, 18, 18, 18, 18,
     $                     28, 28, 28, 28, 28, 28,
     $ 36,36,36,36,36,36,36,36,36,36,36,36,
     $ 46,46,46,46,46,46,54,54,54,0,0,0,
     $ 0,0,0,0,0,0,0,0,0,0,0,68,
     $ 68,68,68,68,68,68,68,68,68,78,78,0,
     $ 0,0,0,0,0,0,0,0,0,0,0,0,
     $ 0,0,0,0,0,0,0   /
      data nelcep /
     $         0, 0,
     $         2, 2, 2, 2, 2, 2, 2, 2,
     $         10, 10, 10, 10, 10, 10, 10, 10,
     $         18, 18,
     $                     10, 10, 10, 10, 10,
     $                     10, 10, 10, 10, 10,
     $                     10, 28, 28, 28, 28, 28,
     $ 36,36,28,28,28,28,28,28,28,28,28,28,
     $ 28,46,46,46,46,46,54,54,46,46,46,46,
     $ 46,46,46,46,46,46,46,46,46,46,46,60,
     $ 60,60,60,60,60,60,60,60,60,78,78,78,
     $ 78,78,0,0,0,0,0,0,0,0,0,0,
     $ 0,0,0,0,0,0,0   /
c.......................................................................
c
c     initialize .
c.......................................................................
      nsym = 0
      maxbas = 4
      do 20 i = 1 , maxbas
         ncsh(i) = 0
         nosh(i) = 0
         nccup(i) = 0
 20   continue
c.......................................................................
      nz = zn + 0.1d0
      nzps = znps + 0.1d0
      idiff = nz - nzps
      ihw = 0
      if (idiff.ne.0) then
c     determine pseudopotential type (hw or cep)
       if (nelhw(nz).eq.idiff) then
         ihw = 1
       else if (nelcep(nz).eq.idiff) then
         ihw = 2
       else
         call nerror(1,'atcond','unrecognised pseudopotential',0,0)
       endif
      endif
      nzps = idiff
c.......................................................................
c
c     fill up shells - all electron first
c.......................................................................
      nlast = 0
       do 50 i = 1 , maxbas
          do 40 j = 1 , 2
             ksym = i
             do 30 k = 1 , i
                nshell = nlast + 4*ksym - 2
                if (nz.lt.nshell) go to 60
                nsym = max0(nsym,ksym)
                ncsh(ksym) = ncsh(ksym) + 1
                ksym = ksym - 1
                nlast = nshell
 30          continue
 40       continue
 50    continue
c.......................................................................
c
c     check if open shell atom. test for la and ac.
c.......................................................................
 60   if (nz.eq.57 .or. nz.eq.89) then
         ncsh(4) = ncsh(4) - ncsh(4)/2
         ksym = 3
      end if
c
c ... now consider pseudopotentials ...
c ... decrease ncsh according to nelcor (nzps)
c
      if (nzps.gt.0) then
       nelec = 0
       if(ihw.eq.1) then
        do 80 i = 1,13
        nelec = nelec + 4*isymhw(i) - 2
        if (nelec.le.nzps) then
         isymm = isymhw(i)
         ncsh(isymm) = ncsh(isymm) - 1
        else
        go to 90
        endif
 80    continue
      else
        do 800 i = 1,13
        nelec = nelec + 4*isymce(i) - 2
        if (nelec.le.nzps) then
         isymm = isymce(i)
         ncsh(isymm) = ncsh(isymm) - 1
        else
        go to 90
        endif
 800   continue
       endif
      endif
c
c     check for cases where pseudo-potential changes  state
c     (cu  pseudo incl d10 =>s1)
c
90    if (nz.eq.29.and.nzps.eq.28) then
         nlast = 28
         ksym = 1
      end if
c
      do 70 i = 1 , 24
         ajmn(i) = 0.0d0
 70   continue
      if (nz.ne.nlast) then
         nosh(ksym) = 1
         nsym = max0(nsym,ksym)
         nccup(ksym) = nz - nlast
c.......................................................................
c
c     set k(l,l,0).
c.......................................................................
         ind = ksym*(ksym+1)*(ksym+2)/6 - ksym + 1
         t = nccup(ksym)
         f = 4*ksym - 4
         ajmn(ind) = -(f+1.0d0)/(f+2.0d0) + (t-1.0d0)/t
      end if
      return
      end
c  =======================================================================
c
      subroutine cmergd(c1,c2,nc,nb,no)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     merge open and closed shell coefficient matrices.
c.......................................................................
      dimension c1(nb,nb), c2(nb,nb)
      nlast = nc + no
      nrow = nc + 1
      do 30 i = nrow , nlast
         do 20 j = 1 , nb
            c1(j,i) = c2(j,i)
 20      continue
 30   continue
      return
      end
c  =======================================================================
c
      subroutine denmad(d,c,ns,nb,occ,nrow)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     make actual density matrices.
c.......................................................................
      dimension d(*), c(nb,nb)
      icount = 0
      do 40 i = 1 , nb
         do 30 j = 1 , i
            icount = icount + 1
            klast = nrow + ns - 1
            sum = 0.0d0
            do 20 k = nrow , klast
               sum = sum + c(i,k)*c(j,k)
 20         continue
            sum = sum*occ
            if (i.ne.j) sum = 2*sum
            d(icount) = sum
 30      continue
 40   continue
      return
      end
c  =======================================================================
c
      subroutine densid(dt,dold,dos,nsym,nosh,ncsh,nccup,c,damp,nconv,
     x                 nbas,nitscf,tlarge)
c.......................................................................
c
c     driver routine for density matrix processing
c.......................................................................
      implicit real*8  (a-h,o-z),integer   (i-n)
      dimension dt(*),dold(*),dos(*),ncsh(*),nosh(*),nccup(*),c(*),
     x          nbas(*)
      nstep1 = 1
      nstep2 = 1
      k = 0
      vamp1 = 1.0d0
      vamp2 = 0.0d0
      if (nitscf.gt.1 .and. nconv.eq.0) vamp1 = 1.0d0 - damp
      if (nitscf.gt.1 .and. nconv.eq.0) vamp2 = damp
      do 40 i = 1 , nsym
         occucl = 4*i - 2
         occuop = nccup(i)
         nbas1 = nbas(i)
         do 30 m = 1 , nbas1
            do 20 n = 1 , m
               k = k + 1
               dt(k) = 0.0d0
               dos(k) = 0.0d0
 20         continue
 30      continue
         if (ncsh(i).ne.0) call denmad(dt(nstep1),c(nstep2),ncsh(i),
     +                                 nbas1,occucl,1)
         if (nosh(i).ne.0) call denmad(dos(nstep1),c(nstep2),nosh(i),
     +                                 nbas1,occuop,ncsh(i)+1)
         nstep1 = nstep1 + nbas1*(nbas1+1)/2
         nstep2 = nstep2 + nbas1**2
 40   continue
      tlarge = 0.0d0
      icount = 0
      do 70 i = 1 , nsym
         do 60 j = 1 , nbas(i)
            do 50 k = 1 , j
               icount = icount + 1
               dt(icount) = (dt(icount)+dos(icount))
     +                      *vamp1 + dold(icount)*vamp2
               ddiff = dabs(dt(icount)-dold(icount))
               dold(icount) = dt(icount)
               if (ddiff.gt.tlarge) then
c                 jmax = j
c                 kmax = k
                  tlarge = ddiff
               end if
 50         continue
 60      continue
 70   continue
      return
      end
c  =======================================================================
c
      subroutine hamild(pcap, qcap, fc, fo, s, u, t, h, dos, dt,
     +                  c, smin, qmin, nbb)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     construct fock matrices. this is a direct s c f procedure,
c     and the two-electron integrals are recalculated for every
c     iteration.
c.......................................................................
      parameter (nbig=100, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
      dimension c(*), smin(nbb,*), qmin(nbb,*)
      dimension pcap(*), qcap(*), fc(*), fo(*), s(*)
      dimension u(*), t(*), h(*), dos(*), dt(*)
c
      call teigd (pcap, qcap, u, t, dt, dos)
c.......................................................................
c
c     compute smin and qmin
c.......................................................................
      nstep1 = 0
      nstep2 = 0
      nstep = 0
      do 50 i = 1 , nsym
         nsh = ncsh(i) + nosh(i)
         nbas1 = nbas(i)
         do 40 j = 1 , nsh
            naddr = nstep2 + (j-1)*nbas1
            j1 = nstep + j
            do 30 m = 1 , nbas1
               smin(m,j1) = 0.d0
               qmin(m,j1) = 0.d0
               do 20 n = 1 , nbas1
                  k = max(m,n)*(max(m,n)-1)/2 + min(m,n) + nstep1
                  smin(m,j1) = smin(m,j1) + s(k)*c(n+naddr)
                  qmin(m,j1) = qmin(m,j1) + qcap(k)*c(n+naddr)
 20            continue
 30         continue
 40      continue
         nstep = nstep + nsh
         nstep2 = nstep2 + nbas1**2
         nstep1 = nstep1 + n1(i)
 50   continue
c.......................................................................
c
c     compute fc and fo
c.......................................................................
      k = 1
      nstep = 0
      do 90 i = 1 , nsym
         occucl = 4*i - 2
         nosh1 = nosh(i)
         nsh = ncsh(i) + nosh1
         fact1 = occucl/(occucl-nccup(i))
         fact2 = nccup(i)/(occucl-nccup(i))
         nbas1 = nbas(i)
         do 80 m = 1 , nbas1
            do 70 n = 1 , m
               term1 = 0.d0
               term2 = 0.d0
               do 60 j = 1 , nsh
                  j1 = j + nstep
                  if (j.ne.nsh .or. nosh1.eq.0) then
                     term1 = term1 + smin(m,j1)*qmin(n,j1) + qmin(m,j1)
     +                       *smin(n,j1)
                  else
                     term2 = smin(m,j1)*qmin(n,j1) + qmin(m,j1)
     +                       *smin(n,j1)
                  end if
 60            continue
               fo(k) = fact1*term1 + pcap(k) + h(k) - qcap(k)
               fc(k) = fact2*term2 + pcap(k) + h(k)
               k = k + 1
 70         continue
 80      continue
         nstep = nstep + nsh
 90   continue
      return
      end
c  =======================================================================
c
      subroutine jacod(f,v,nb,nb1,nb2,nmin,nmax,big,jbig,maxao)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c      f is the matrix to be diagonalized. f is stored triangular
c      v is the array of eigenvectors. quadratic array, dimension nb*nb
c      big and jbig are temporary scratch areas of dimension nb
c      the rotations among the first nmin basis functions are not
c      accounted for.
c      the rotations among the last nb-nmax basis functions are not
c      accounted for.
c.......................................................................
      dimension big(nb),jbig(nb),f(nb1),v(nb2)
      data root2 /0.707106781186548d0/
c     data c1/1.d-12/
      data c2,c3,c4,c5,c6/1.d-12,4.d-16,2.d-16,8.d-9,3.d-6/
c...    ligen ignored
c     if(iligen.ne.1) goto 17171
c     call ligen(f,v,jbig,nb,maxao,big,1.e-13)
c     return
c       17171 continue
      if (nb.eq.1) then
         v(1) = 1.0d0
         return
      end if
      ii = 0
c.......................................................................
c      loop over rows (i) of triangular matrix
c.......................................................................
      do 30 i = 1 , nb
         big(i) = 0.00d0
         jbig(i) = 0
         if (i.ge.nmin .and. i.ne.1) then
            j = min0(i-1,nmax)
c.......................................................................
c      loop over columns (k) of triangular matrix to determine
c      largest off-diagonal elements in row(i).
c.......................................................................
            do 20 k = 1 , j
               if (dabs(big(i)).lt.dabs(f(ii+k))) then
                  big(i) = f(ii+k)
                  jbig(i) = k
               end if
 20         continue
         end if
         ii = ii + i
 30   continue
 40   sd = 1.050d0
c.......................................................................
c      find smallest diagonal element and corresponding largest
c      off-diagonal element.
c.......................................................................
      jj = 0
      do 50 j = 1 , nb
         jj = jj + j
         sd = dmin1(sd,dabs(f(jj)))
 50   continue
      sd = dmax1(sd,c6)*c2
      t = 0.0d0
      i1 = max0(2,nmin)
      do 60 i = i1 , nb
         if (t.lt.dabs(big(i))) then
            t = dabs(big(i))
            ib = i
         end if
 60   continue
c.......................................................................
c      test for convergence, then determine rotation.
c.......................................................................
      ia = jbig(ib)
      if (t.lt.sd) then
         return
      else
         iaa = ia*(ia-1)/2
         ibb = ib*(ib-1)/2
         jaa = (ia-1)*maxao
         jbb = (ib-1)*maxao
         dif = f(iaa+ia) - f(ibb+ib)
         if (dabs(dif).gt.c3*t) then
            t2x2 = big(ib)/dif
            t2x25 = t2x2*t2x2
            if (t2x25.le.c4) then
               cx = 1.0d0
               sx = t2x2
            else if (t2x25.le.c5) then
               sx = t2x2*(1.0d0-1.5d0*t2x25)
               cx = 1.0d0 - 0.5d0*t2x25
            else if (t2x25.gt.c6) then
               t = 0.25d0/dsqrt(0.25d0+t2x25)
               cx = dsqrt(0.5d0+t)
               sx = dsign(dsqrt(0.5d0-t),t2x2)
            else
               cx = 1.0d0 + t2x25*(t2x25*1.375d0-0.5d0)
               sx = t2x2*(1.00d0+t2x25*(t2x25*3.875d0-1.5d0))
            end if
         else
            sx = root2
            cx = root2
         end if
         iear = iaa + 1
         iebr = ibb + 1
         do 90 ir = 1 , nb
            t = f(iear)*sx
            f(iear) = f(iear)*cx + f(iebr)*sx
            f(iebr) = t - f(iebr)*cx
            if (ir.lt.ia) then
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else if (ir.eq.ia) then
               tt = f(iebr)
               ieaa = iear
               ieab = iebr
               f(iebr) = big(ib)
               iear = iear + ir - 1
               if (jbig(ir).ne.0) go to 70
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else
               t = f(iear)
               it = ia
               iear = iear + ir - 1
               if (ir.lt.ib) then
               else if (ir.eq.ib) then
                  f(ieaa) = f(ieaa)*cx + f(ieab)*sx
                  f(ieab) = tt*cx + f(iebr)*sx
                  f(iebr) = tt*sx - f(iebr)*cx
                  iebr = iebr + ir - 1
                  go to 70
               else
                  if (dabs(t).lt.dabs(f(iebr))) then
                     if (ib.le.nmax) then
                        t = f(iebr)
                        it = ib
                     end if
                  end if
                  iebr = iebr + ir - 1
               end if
            end if
            if (dabs(t).ge.dabs(big(ir))) then
               big(ir) = t
               jbig(ir) = it
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else if (ia.ne.jbig(ir) .and. ib.ne.jbig(ir)) then
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            end if
 70         kq = iear - ir - ia + 1
            big(ir) = 0.00d0
            ir1 = min0(ir-1,nmax)
            do 80 i = 1 , ir1
               k = kq + i
               if (dabs(big(ir)).lt.dabs(f(k))) then
                  big(ir) = f(k)
                  jbig(ir) = i
               end if
 80         continue
            iear = iear + 1
            iebr = iebr + 1
 90      continue
         do 100 i = 1 , maxao
            t = v(jbb+i)*sx
            v(jbb+i) = v(jaa+i)*sx - v(jbb+i)*cx
            v(jaa+i) = v(jaa+i)*cx + t
 100     continue
         go to 40
      end if
      end
c  =======================================================================
c
      subroutine oeigd(fc, s, u, t, h)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     one-electron integrals. general.
c.......................................................................
      parameter (nbig=100, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
      dimension fc(*),s(*),u(*),t(*),h(*)
      ufacl = dsqrt(8.0d0/3.1415926536d0)
      nstep = 0
      k = 0
      do 40 l = 1 , nsym
         nbas1 = nbas(l)
         expfac = l + 0.5d0
         do 30 i = 1 , nbas1
            do 20 j = 1 , i
               k = k + 1
               zp = zeta(nstep+i)
               zq = zeta(nstep+j)
               zpq = 0.5d0*(zp+zq)
               term1 = dsqrt(zpq)
               ppq = zp*zq/zpq
               rpq = dsqrt(ppq/zpq)
               s(k) = rpq**expfac
               fc(k) = s(k)
               u(k) = ufacl*s(k)*term1
               t(k) = expfac*s(k)*ppq
               h(k) = t(k) - zn*u(k)
 20         continue
 30      continue
         nstep = nstep + nbas1
         ufacl = ufacl*2.0d0*l/(2.0d0*l+1.0d0)
 40   continue
      return
      end
c  =======================================================================
c
      subroutine orderd(amat,nrow,ncol,imemb,nblock,ind,vec,iflag)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c     this routine sorts a set of column vectors in amat(*,*) according
c     to increasing values in vec(*). care is taken that the output orde
c     of the vectors is as close to the input order as possible.
c.......................................................................
      dimension amat(nrow,ncol), ind(*), vec(*), imemb(*)
c.......................................................................
c
c     set tolerance for test, then determine for each element how
c     many smaller elements vec(*) contains.
c.......................................................................
      tol = 1.0d-10
      do 30 i = 1 , ncol
         test = vec(i) - tol
         indi = 1
         do 20 j = 1 , ncol
            if (vec(j).lt.test) indi = indi + 1
 20      continue
         ind(i) = indi
 30   continue
c.......................................................................
c
c     if desired,scan ind(*) to determine the number of different values
c     in vec(*) and the number of elements of each value.
c.......................................................................
      if (iflag.le.99) then
         do 40 i = 1 , ncol
            imemb(i) = 0
 40      continue
         do 50 i = 1 , ncol
            indi = ind(i)
            imemb(indi) = imemb(indi) + 1
 50      continue
         icount = 0
         do 60 i = 1 , ncol
            if (imemb(i).ne.0) then
               icount = icount + 1
               imemb(icount) = imemb(i)
            end if
 60      continue
         nblock = icount
      end if
c.......................................................................
c
c     establish order in degeneracies of the ordering vector.
c.......................................................................
      do 80 i = 2 , ncol
         itest = ind(i)
         im1 = i - 1
         do 70 j = 1 , im1
            if (ind(j).eq.itest) itest = itest + 1
 70      continue
         ind(i) = itest
 80   continue
c.......................................................................
c
c     ind(*) contains ordering indices for amat(*,*). sort following
c     input order as far as this is correct.
c.......................................................................
      ilow = 1
 90   do 100 i = 1 , ncol
         if (i.ne.ind(i)) go to 110
 100  continue
      go to 130
c.......................................................................
c
c     input order wrong. swap present vector into correct col.
c.......................................................................
 110  index = ind(i)
      do 120 j = 1 , nrow
         scra = amat(j,i)
         amat(j,i) = amat(j,index)
         amat(j,index) = scra
 120  continue
      scra = vec(i)
      vec(i) = vec(index)
      vec(index) = scra
      ind(i) = ind(index)
      ind(index) = index
      if (ind(i).ne.i) go to 110
c.......................................................................
c
c     amat(*,*) is ordered through i. go back to checking order.
c.......................................................................
      ilow = i + 1
      if (ilow.lt.ncol) go to 90
 130  return
      end
c  =======================================================================
c
      subroutine outpud(copn,cc,ntest,iwr)
      implicit real*8  (a-h,o-z),integer   (i-n)
      parameter (nbig=100, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
      dimension cc(*),copn(*)
      if (ntest.gt.0) write (iwr,6010) nitscf , energ , cin , vir
      if (ntest.gt.0) write (iwr,6030)
      nbc1 = nbc(1)
      do 20 i = 2 , nsym
         nbc1 = max0(nbc1,nbc(i))
 20   continue
      naddr = 0
      noddr = 0
      do 60 i = 1 , nsym
         do 50 j = 1 , ncsh(i) + nosh(i)
            do 30 k = 1 , nbc(i)
               naddr = naddr + 1
               noddr = noddr + 1
               copn(noddr) = cc(naddr)
 30         continue
            if (nbc(i).lt.nbc1) then
               do 40 k = nbc(i) + 1 , nbc1
                  noddr = noddr + 1
                  copn(noddr) = 0.0d0
 40            continue
            end if
 50      continue
         naddr = naddr + nbc(i)*(nbc(i)-ncsh(i)-nosh(i))
 60   continue
      jfirst = 1
 70   jlast = min(jfirst+7,nsht)
      if (ntest.gt.0) then
         write (iwr,6020) (eps(j),j=jfirst,jlast)
         write (iwr,6030)
         do 80 i = 1 , nbc1
            write (iwr,6040) (copn(i+(j-1)*nbc1),j=jfirst,jlast)
 80      continue
      end if
      if (jlast.eq.nsht) return
      jfirst = jfirst + 8
      if (ntest.gt.0) write (iwr,6030)
      go to 70
 6010 format (/,8x,'final scf results at iteration',i4,/,8x,
     +        'total hf energy',4x,'kinetic energy',4x,'virial theorem',
     +        /4x,3(e19.10),//,8x,'orbital energies and eigenvectors')
 6020 format (/,2x,8(2x,f12.5),1x)
 6030 format (' ')
 6040 format (2x,8f14.6)
      end
c  =======================================================================
c
      subroutine shalfd(s,v,n,maxao)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     make s**(-1/2) for transformation to orthonormal basis.
c.......................................................................
      common/blkin/scr1(100),scr2(100)
      dimension s(*), v(n,n)
      icount = 0
      nnd = n*(n+1)/2
      nn = n*n
      if (n.gt.maxao) call nerror(1,'shalfd','dimensioning error',
     $                            n,maxao)
      do 30 i = 1 , n
         do 20 j = 1 , i
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
 20      continue
         v(i,i) = 1.0d0
 30   continue
      call jacod(s,v,n,nnd,nn,1,n,scr1,scr2,n)
      icount = 1
      do 40 i = 1 , n
         scr1(i) = 1.0d0/dsqrt(s(icount))
         icount = icount + i + 1
 40   continue
      icount = 0
      do 70 i = 1 , n
         do 60 j = 1 , i
            icount = icount + 1
            hlp = 0.0d0
            do 50 k = 1 , n
               hlp = hlp + v(i,k)*scr1(k)*v(j,k)
 50         continue
            s(icount) = hlp
 60      continue
 70   continue
      icount = 0
      do 90 i = 1 , n
         do 80 j = 1 , i
            icount = icount + 1
            v(j,i) = s(icount)
            v(i,j) = s(icount)
 80      continue
 90   continue
      return
      end
c  =======================================================================
c
      subroutine starcd(c,ss,nbci,ncshi,noshi)
      implicit real*8  (a-h,o-z),integer   (i-n)
      dimension ss(*), c(nbci,nbci)
      common/blkin/a(100)
c.......................................................................
c
c     this routine provides a set of schmidt
c     orthogonalized start vectors.
c
c     inline indexing function for triangular matrices.
c
c.......................................................................
      index(i,j) = max0(i,j)*(max0(i,j)-1)/2 + min0(i,j)
c.......................................................................
c
c     uneducated guess of trial vectors.
c
c.......................................................................
      nshti = ncshi + noshi
      if (nshti.eq.0) return
      if (nshti.gt.nbci) call nerror(1,'starcd',
     +  'more atomic shells than basis functions in atomic startup',
     $   nsht1,nbci)
      ixx = nbci/nshti
      khi = 0
      do 30 j = 1 , nshti
         klo = khi + 1
         khi = khi + ixx
         do 20 k = klo , khi
            c(k,j) = 1.0d0
 20      continue
 30   continue
c.......................................................................
c
c     orthogonalize.
c
c.......................................................................
      do 130 j = 1 , nshti
c.......................................................................
c
c     take scalar products with preceding vectors of same symm.
c     first do a(*) = ss(*)*c(*,j).
c
c.......................................................................
         do 50 l = 1 , nbci
            sum = 0.0d0
            do 40 m = 1 , nbci
               indxlm = index(l,m)
               sum = sum + ss(indxlm)*c(m,j)
 40         continue
            a(l) = sum
 50      continue
c.......................................................................
c
c     if first vector of symmetry, normalize directly.
c
c.......................................................................
         if (j.ne.1) then
c.......................................................................
c
c     c(*,k)*a(*)
c
c.......................................................................
            jm1 = j - 1
            do 80 k = 1 , jm1
               sum = 0.0d0
               do 60 l = 1 , nbci
                  sum = sum + c(l,k)*a(l)
 60            continue
c.......................................................................
c
c     multiply c(*,k) by the scalar product and subtract from c(*,j)
c
c.......................................................................
               do 70 l = 1 , nbci
                  c(l,j) = c(l,j) - c(l,k)*sum
 70            continue
 80         continue
c.......................................................................
c
c     c(*,j) is orthogonal to all preceding vectors. normalize.
c
c.......................................................................
            do 100 l = 1 , nbci
               sum = 0.0d0
               do 90 m = 1 , nbci
                  indxlm = index(l,m)
                  sum = sum + ss(indxlm)*c(m,j)
 90            continue
               a(l) = sum
 100        continue
         end if
         sum = 0.0d0
         do 110 l = 1 , nbci
            sum = sum + c(l,j)*a(l)
 110     continue
         sum = dsqrt(sum)
         do 120 l = 1 , nbci
            c(l,j) = c(l,j)/sum
 120     continue
 130  continue
      return
      end
c  =======================================================================
c
      subroutine teigd(pcap,qcap,u,t,dt,dos)
      implicit real*8  (a-h,o-z),integer   (i-n)
      logical open,klnemn
c.......................................................................
c
c     two-electron integral routine for s,p,d, and f functions.
c.......................................................................
      parameter (nbig=100, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
      dimension twopow(14), pifac(14)
      dimension pcap(*),qcap(*),u(*),t(*),dos(*),dt(*)
c.......................................................................
c
c     angular factors for exchange integrals. obtained as sum of squares
c     of slater coefficients c(kappa;l1,m1;l2,m2) divided by
c     2*(2*l1+1)*(2*l2+1) for any kappa, l1, and l2.
c.......................................................................
      data ss0,sp1,pp0,pp2,sd2,pd1,pd3,dd0,dd2,
     x     sf3,pf2,pf4,df1,df3,df5,ff0,ff2,ff4,ff6
     x/ .50000000000d+00, .16666666667d+00, .16666666667d+00,
     x  .66666666667d-01, .10000000000d+00, .66666666667d-01,
     x  .42857142857d-01, .10000000000d+00, .28571428571d-01,
     x                    .71428571429d-01, .42857142857d-01,
     x  .31746031746d-01, .42857142857d-01, .19047619048d-01,
     x  .21645021645d-01, .71428571429d-01, .19047619048d-01,
     x  .12987012987d-01, .16650016650d-01/
      f0pol(a,b) = 3*(16*a**6+104*a**5*b+286*a**4*b**2+429*(a*b)
     +             **3+286*a**2*b**4+104*a*b**5+16*b**6)
      f2pol(a,b) = 8*a**4 + 52*a**3*b + 143*(a*b)**2 + 52*a*b**3 +
     +             8*b**4
      df1pol(a,b) = 8*a**4 + 44*a**3*b + 99*(a*b)**2 + 44*a*b**3 +
     +              8*b**4
      df3pol(a,b) = 2*a**2 + 11*a*b + 2*b**2
c.......................................................................
c     two-electron integral routine for lcgo atom scf.
c     restricted to principal quantum numbers 1,2 and 3 for respectively
c     s,p, and d orbitals.
c.......................................................................
c     pi = 3.14159265d0
c     on6 = 1.0d0/6.0d0
c     on15 = 1.0d0/15.d0
c     on35 = 2*on70
c     on70 = 1.0d0/70.d0
      pifac(1) = dsqrt(3.14159265d0)
      twopow(1) = 1
      do 20 i = 2 , 14
         twopow(i) = twopow(i-1)*0.5d0
         pifac(i) = pifac(i-1)*(2*i-1)
 20   continue
c.......................................................................
c
c     this part sets up the coefficients lambda,p,q and mu,r,s.
c.......................................................................
      j = 0
c     nstep1 = 0
      kmx = 0
      factkl = 1
      kl = 0
      pot = 0.0d0
      potn = 0.0d0
      cin = 0.0d0
      do 190 i = 1 , nsym
         prfac1 = twopow(i+1)*factkl
         kin = kmx + 1
         kmx = kin + nbas(i) - 1
         do 180 k = kin , kmx
            zp = zeta(k)
            do 170 l = kin , k
               kl = kl + 1
               pcap(kl) = 0.0d0
               qcap(kl) = 0.0d0
               zq = zeta(l)
               zpq = zp + zq
               prfac2 = prfac1*u(kl)
               xfac1 = prfac2*zpq**i*pifac(i)
c              nstep2 = 0
               mmx = 0
               factmn = 1.0d0
               mn = 0
               do 160 im = 1 , i
                  open = (nosh(i).ne.0 .and. nosh(im).ne.0)
                  prfac3 = prfac2*pifac(im)*twopow(im)*factmn
                  xfac2 = xfac1*factmn*twopow(im+1)
                  min = mmx + 1
                  mmx = min - 1 + nbas(im)
                  mmxp = mmx
                  if (im.eq.i) mmxp = k
                  do 150 m = min , mmxp
                     zr = zeta(m)
                     zpr = zp + zr
                     zqr = zq + zr
                     nmx = l
                     if (m.lt.k) nmx = m
                     do 140 n = min , nmx
                        mn = mn + 1
                        klnemn = (kl.ne.mn)
                        j = j + 1
c.......................................................................
c
c     j is the number label of the matrix elements to be calculated
c     i=lambda+1,k=p,l=q,im=mu+1,m=r,n=s
c.......................................................................
                        zs = zeta(n)
                        zrs = zr + zs
                        zqs = zq + zs
                        zps = zp + zs
                        zpqrs = zpq + zrs
                        zpqrs2 = 2*zpqrs**2
                        zprzqs = zpr*zqs
                        zpszqr = zps*zqr
                        xterm = (1.0d0/dsqrt(zpqrs))**(2*(i+im)-3)
                        prfac4 = prfac3*u(mn)*xterm
                        xfac3 = xfac2*u(mn)*xterm
                        xfac11 = xfac3*(zrs/zprzqs)**im
                        xfac21 = xfac3*(zrs/zpszqr)**im
                        xfsum = xfac11 + xfac21
                        ntest = i*(i-1)/2 + im
                        go to (30,40,50,60,70,80,90,100,110,120,140) ,
     +                         ntest
c.......................................................................
c
c     i=1,im=1,(ss)-loop. x0=j0(ss),y0=k0(ss)
c.......................................................................
 30                     x0 = prfac4
                        y0 = xfsum
                        pj = x0 - y0*ss0
                        qj = -ajmn(1)*y0
                        go to 130
c.......................................................................
c
c     i=2,im=1,(sp)-loop. x0=j0(sp),y1=k1(sp)
c.......................................................................
 40                     x0 = prfac4*(3*zpq+2*zrs)
                        y1 = xfsum
                        pj = x0 - y1*sp1
                        qj = -ajmn(2)*y1
                        go to 130
c.......................................................................
c
c     i=2,im=2,(pp)-loop. x0=j0(pp),y0=k0(pp),y2=k2(pp)
c.......................................................................
 50                     x0 = prfac4*(zpqrs2+zpq*zrs)
                        y0 = xfsum*zpqrs2
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
                        y0 = y0 + xfsum
                        y2 = xfsum*5
                        pj = x0 - y0*pp0 - y2*pp2
                        qj = -(ajmn(3)*y0+ajmn(4)*y2)
                        go to 130
c.......................................................................
c
c     i=3,im=1,(sd)-loop. x0=j0(sd),y2=k2(sd)
c.......................................................................
 60                     x0 = prfac4*(15*zpq**2+20*zpq*zrs+8*zrs**2)
                        y2 = xfsum
                        pj = x0 - y2*sd2
                        qj = -ajmn(5)*y2
                        go to 130
c.......................................................................
c
c     i=3,im=2,(pd)-loop. x0=j0(pd),y1=k1(pd),y3=k3(pd),x2=j2(pd)
c.......................................................................
 70                     x0 = prfac4*(10*zpq**3+35*zpq**2*zrs+
     +                       28*zpq*zrs**2+8*zrs**3)
                        y1 = xfsum*zpqrs2
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
                        y1 = y1 + 3*xfsum
                        y3 = 7*xfsum
                        pj = x0 - y1*pd1 - y3*pd3
                        qj = 0.d0
                        if (open) then
                           x2 = prfac4*5*zpq*zrs*(7*zpq+2*zrs)
                           qj = ajmn(21)*x2 - (ajmn(6)*y1+ajmn(7)*y3)
                        end if
                        go to 130
c.......................................................................
c
c     i=3,im=3,(dd)-loop. x0=j0(dd),y0=k0(dd),y2=k2(dd),y4=k4(dd)
c.......................................................................
 80                     zprzqs = zpr*zqs
                        zpszqr = zps*zqr
                        x0 = prfac4*((zpqrs2+zpq*zrs)
     +                       *zpqrs2*2+7*(zpq*zrs)**2)
                        y01 = xfac11*((zpqrs2+zprzqs)
     +                        *zpqrs2*2+7*zprzqs**2)
                        y02 = xfac21*((zpqrs2+zpszqr)
     +                        *zpqrs2*2+7*zpszqr**2)
                        y0 = y01 + y02
                        xfac11 = xfac11*7*zprzqs
                        xfac21 = xfac21*7*zpszqr
                        y21 = xfac11*(zpqrs2+5*zprzqs)
                        y22 = xfac21*(zpqrs2+5*zpszqr)
                        y2 = y21 + y22
                        y4 = (xfac11*zprzqs+xfac21*zpszqr)*9
                        pj = x0 - y0*dd0 - (y2+y4)*dd2
                        qj = -ajmn(8)*y0 - ajmn(9)*y2 - ajmn(10)*y4
                        go to 130
c.......................................................................
c
c     i=4,im=1,(sf)-loop. x0=j0(sf),y3=k3(sf)
c.......................................................................
 90                     x0 = prfac4*(105*zpq**3+210*zpq**2*zrs+
     +                       168*zpq*zrs**2+48*zrs**3)
                        y3 = xfsum
                        pj = x0 - y3*sf3
                        qj = -ajmn(11)*y3
                        go to 130
c.......................................................................
c
c     i=4,im=2,(pf)-loop. x0=j0(pf),y2=k2(pf),y4=k4(pf),x2=j2(pf)
c.......................................................................
 100                    x0 = prfac4*(70*zpq**4+315*zpq**3*zrs+
     +                       378*(zpq*zrs)**2+216*zpq*zrs**3+48*zrs**4)
c.......................................................................
c                       y2 = xfsum*zpqrs2
c.......................................................................
                        y2 = xfac11*(2*zpr**2+9*zpr*zqs+2*zqs**2)
     +                       + xfac21*(2*zps**2+9*zps*zqr+2*zqr**2)
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
c.......................................................................
c                       y2 = y2 + 5*xfsum
c.......................................................................
                        y4 = 9*xfsum
                        pj = x0 - y2*pf2 - y4*pf4
                        qj = 0.0d0
                        if (open) then
                           x2 = prfac4*5*zpq*zrs*(63*zpq**2+26*zpq*zrs+
     +                          8*zrs**2)
                           qj = ajmn(22)*x2 - ajmn(12)*y2 - ajmn(13)*y4
                        end if
                        go to 130
c.......................................................................
c
c     i=4,im=3,(df)-loop. x0=j0(df),y1=k1(df),y3=k3(df),y5=k5(df)
c                         x2=j2(df),x4=j4(df).
c.......................................................................
 110                    x0 = prfac4*(56*zpq**5+308*zpq**4*zrs+
     +                       693*zpq**3*zrs**2+594*zpq**2*zrs**3+
     +                       264*zpq*zrs**4+48*zrs**5)
                        y1 = xfac11*df1pol(zpr,zqs)
     +                       + xfac21*df1pol(zps,zqr)
                        xfac11 = 9*zpr*zqs*xfac11
                        xfac21 = 9*zps*zqr*xfac21
                        y3 = xfac11*df3pol(zpr,zqs)
     +                       + xfac21*df3pol(zps,zqr)
                        y5 = 11*(xfac11*zpr*zqs+xfac21*zps*zqr)
                        pj = x0 - y1*df1 - y3*df3 - y5*df5
                        qj = 0.0d0
                        if (open) then
                           prfac4 = prfac4*7*zpq*zrs
                           x2 = prfac4*(18*zpq**3+99*zpq**2*zrs+
     +                          44*zpq*zrs**2+8*zrs**3)
                           x4 = 9*prfac4*(11*zpq+2*zrs)*zpq*zrs
                           qj = x2*ajmn(23) + x4*ajmn(24) - y1*ajmn(14)
     +                          - y3*ajmn(15) - y5*ajmn(16)
                        end if
                        go to 130
c.......................................................................
c
c     i=4,im=4,(ff)-loop. x0=j0(ff),y0=k0(ff),y2=k2(ff),y4=k4(ff)
c              y6=k6(ff)
c.......................................................................
 120                    x0 = prfac4*f0pol(zpq,zrs)
                        y0 = xfac11*f0pol(zpr,zqs)
     +                       + xfac21*f0pol(zps,zqr)
                        xfac11 = xfac11*zpr*zqs*9
                        xfac21 = xfac21*zps*zqr*9
                        y2 = xfac11*f2pol(zpr,zqs)
     +                       + xfac21*f2pol(zps,zqr)
                        xfac11 = xfac11*zpr*zqs*11
                        xfac21 = xfac21*zps*zqr*11
                        y4 = xfac11*(2*zpr**2+13*zpr*zqs+2*zqs**2)
     +                       + xfac21*(2*zps**2+13*zps*zqr+2*zqr**2)
                        y6 = 13*(xfac11*zpr*zqs+xfac21*zps*zqr)
                        pj = x0 - y0*ff0 - y2*ff2 - y4*ff4 - y6*ff6
                        qj = -y0*ajmn(17) - y2*ajmn(18) - y4*ajmn(19)
     +                       - y6*ajmn(20)
 130                    pcap(kl) = pcap(kl) + pj*dt(mn)
                        if (klnemn) pcap(mn) = pcap(mn) + pj*dt(kl)
                        term = dt(kl)*pj*dt(mn)
                        if (open) then
                           qcap(kl) = qcap(kl) + qj*dos(mn)
                           if (klnemn) qcap(mn) = qcap(mn) + qj*dos(kl)
                           term = term - dos(kl)*qj*dos(mn)
                        end if
                        if (.not.klnemn) term = term*0.5d0
                        pot = pot + term
 140                 continue
 150              continue
                  factmn = factmn/im
 160           continue
               potn = potn + u(kl)*dt(kl)
               cin = cin + t(kl)*dt(kl)
 170        continue
 180     continue
         factkl = factkl/i
 190  continue
c     istart = 1
      pot = pot - zn*potn
      energ = cin + pot
      vir = pot/cin
      return
      end
c  =======================================================================
c
      subroutine tracd(c,cc,nsqt)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     transform vectors from contr. to primitive basis functions.
c.......................................................................
      parameter (nbig=100, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont
      dimension c(*),cc(*)
c
      do 20 i = 1 , nsqt
         c(i) = 0.0d0
 20   continue
c.......................................................................
c
c     cc(i,j)  contains vectors over cont. functions.
c              i runs over functions, j runs over orbitals
c.......................................................................
      nstep = 0
      nstep1 = 0
      nstep2 = 0
      do 60 l = 1 , nsym
         nbc1 = nbc(l)
         nbasl = nbas(l)
         nsh1 = ncsh(l) + nosh(l)
         ii = 1
         do 50 n = 1 , nbc1
            k1 = nstrt(nstep+n)
            k2 = nstrt(nstep+n+1) - 1
            do 40 k = k1 , k2
               do 30 j = 1 , nsh1
                  c(nstep2+ii+(j-1)*nbasl) = cc(nstep1+n+(j-1)*nbc1)
     +               *cont(k)
 30            continue
               ii = ii + 1
 40         continue
 50      continue
         nstep = nstep + nbc1
         nstep1 = nstep1 + nbc1**2
         nstep2 = nstep2 + nbasl**2
 60   continue
      return
      end
c  =======================================================================
c
      subroutine trafsd(nsym,nbas,ndim,a,nbc,cont,nstrt,ffc)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c     transform matrix a
c     from prim. to cont. basis.
c.......................................................................
      dimension ffc(*),nbas(*),a(*),nbc(*),cont(*),nstrt(*)
      do 20 i = 1 , ndim
         ffc(i) = a(i)
 20   continue
      mnstep = 0
      ijstep = 0
      k = 1
      ijbas = 0
      mnbas = 0
      do 70 l = 1 , nsym
         nbc1 = nbc(l)
         do 60 m = 1 , nbc1
            do 50 n = 1 , m
               m1 = nstrt(mnbas+m)
               n0 = nstrt(mnbas+n)
               m2 = nstrt(mnbas+m+1) - 1
               n2 = nstrt(mnbas+n+1) - 1
               sumc = 0.0d0
               do 40 i = m1 , m2
                  ij = ijstep + (i-ijbas-1)*(i-ijbas)/2
                  if (m.eq.n) then
                     n2 = i
                  end if
                  do 30 j = n0 , n2
                     jj = j - ijbas
                     if (m.ne.n) then
                        sumc = sumc + cont(i)*cont(j)*ffc(ij+jj)
                     else if (i.ne.j) then
                        sumc = sumc + 2.0d0*cont(i)*cont(j)*ffc(ij+jj)
                     else
                        sumc = sumc + cont(i)*cont(j)*ffc(ij+jj)
                     end if
 30               continue
 40            continue
               a(k) = sumc
               k = k + 1
 50         continue
 60      continue
         mnstep = mnstep + nbc(l)*(nbc(l)+1)/2
         ijstep = ijstep + nbas(l)*(nbas(l)+1)/2
         mnbas = mnbas + nbc(l)
         ijbas = ijbas + nbas(l)
 70   continue
      return
      end
c  =======================================================================
c
      subroutine tramad(a,b,c,mdima,mdim,scr)
      implicit real*8  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     this routine tranforms a to (b+)ab, where a is a lower
c     triangular symmetric matrix, and b is a square matrix.
c     the transformed matrix is returned in a.
c.......................................................................
      dimension a(mdima), b(mdim,mdim), c(mdima), scr(mdim)
      do 60 i = 1 , mdim
c.......................................................................
c
c     generate i'th column of ab.
c.......................................................................
         do 30 j = 1 , mdim
            sum = 0.0d0
            do 20 l = 1 , mdim
               maxjl = max0(j,l)
               jl = maxjl*(maxjl-3)/2 + j + l
               sum = sum + a(jl)*b(l,i)
 20         continue
            scr(j) = sum
 30      continue
c.......................................................................
c
c     multiply this by rows of b+
c.......................................................................
         do 50 j = 1 , i
            sum = 0.0d0
            do 40 k = 1 , mdim
               sum = sum + scr(k)*b(k,j)
 40         continue
            c(i*(i-1)/2+j) = sum
 50      continue
 60   continue
c.......................................................................
c
c     transfer to a for return.
c.......................................................................
      do 70 i = 1 , mdima
         a(i) = c(i)
 70   continue
      return
      end
