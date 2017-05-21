c ==================================================================
c  POPULATION ANALYSIS MODULE            JB   Dec 1999
c  modified for semiempirical            JB   Nov 2008
c ==================================================================
c
      subroutine prepop(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval
c
c  reads the POP line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=3)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character cdum*20
      character*256 jobname
c
      parameter (IUnit=1)
c
      data options/'pop ','pthr','prin'/
      data ioptyp/21,11,1/
c
c -- deal with jobname header
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c -- what type of population analysis?
c
      cdum = chopval(1)
      call lowerca2(cdum,20)
      if(cdum(1:4).eq.'    ') cdum='standard'
c
c -- bond order print threshold
      If(ifound(2).eq.1) then
        PThrsh = ropval(1,2)
      else
        PThrsh = 0.01d0
      endif
c
      if(ifound(3).eq.1) then
        IPRNT = iopval(1,3)
      else
        IPRNT = 2          ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,4,'$pop',3,idum,rdum,cdum)
      call wrcntrl(IUnit,7,'$pthrsh',2,idum,PThrsh,cdum)
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c  =======================================================================
c
      SUBROUTINE POPML(NMem,Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** POPULATION ANALYSIS MODULE **
C
C  This module carries out standard Mulliken and Lowdin population
C  analysis, forming and printing
C    gross orbital occupancies
C    atomic charges
C    bond orders
C    atomic valencies
C    free valencies (for open-shell wavefunctions)
C
C  ..........................................................
C
C  POPML itself is simply a "wrapper" which reads the
C  <control> file to determine the size of the current
C  system and the program option requested and, based on
C  this information, allocates the necessary memory.
C  It then calls POPMAIN which is itself another "wrapper"
C  which completes job input and is responsible for all
C  other I/O. POPMAIN calls <popanal> which is the main
C  analysis routine.
C  .............................................................
C
C
      CHARACTER jobname*256,wvfnc*20,cdum*20
      Logical Semi
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      Data Semi/.False./
c
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    number of atoms
C    number of alpha/closed-shell orbitals
C    number of beta orbitals
C    multiplicity
C    wavefunction type
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,dum,cdum)
      call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Check the wavefunction type
C
      If(wvfnc(1:7).EQ.'Semiemp') Semi=.True.
      If(wvfnc(2:3).EQ.'MP') call nerror(3,'Population Analysis module',
     $   'Population Analysis Unavailable for MP2 Wavefunctions',0,0)
c
c ............................................................
c -- for the population analysis, dummy atoms are ignored
c
      NAtoms = NAtoms-Ndum1-Ndum2
c ............................................................
c
      If(.NOT.Semi) Then
C
C  Read from the <basis> file
C    number of basis functions, primitives and shells
C
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $        FORM='FORMATTED',STATUS='OLD')
        call rdcntrl(IUnit,7,'$nbasis',1,NBas,dum,cdum)
        call rdcntrl(IUnit,7,'$nshell',1,NShl,dum,cdum)
        call rdcntrl(IUnit,6,'$nprim',1,NPrm,dum,cdum)
        CLOSE(UNIT=IUnit,STATUS='KEEP')
      Else
C
C  for semiempirical (over)estimate the number of basis functiona
C
        NBas = 4*NAtoms
        NShl = 0
        NPrm = 0
      EndIf
C
C
C  Now get the memory
C
      NScr = MAX(127008,NBas*(NBas+1))
      IMem = 6*NAtoms + 12*NShl + 13*NPrm + 4*NBas**2 + NBas +NAtoms**2
      If(NBeta.GT.0) IMem = IMem + 2*NBas**2
      IMem = IMem + NScr
c
cc      CALL falloc(Z(0),8*IMem,iptr,IErr)
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,5,'POPML')
C
C
C  Allocate memory pointers
C
      IAN   = iptr                      !  atomic numbers
      IQA   = IAN  + NAtoms             !  actual atomic charge
      IXC   = IQA  + NAtoms             !  geometry
      IBAS  = IXC  + 3*NAtoms           !  basis function data
      INX   = IBAS + 13*NPrm            !  basis indexing array
      IBA   = INX  + 12*NShl            !  number of basis functions per atom
      IS    = IBA  + NAtoms             !  overlap matrix
      ISH   = IS   + NBas*NBas          !  S**1/2
      IPA   = ISH  + NBas*NBas          !  density matrix
      IPS   = IPA  + NBas*NBas          !  PS matrix
      IG    = IPS  + NBas*NBas          !  gross orbital occupancies
      IBO   = IG   + NBas               !  bond order matrix
      IEnd  = IBO  + NAtoms*NAtoms
c
      IF(NBeta.GT.0) THEN
       IPB  = IEnd                      !  beta density matrix
       IEnd = IPB  + NBas*NBas
      ELSE
       IPB  = 1
      ENDIF
C
C  scratch memory
C
      IScr = IEnd                       !  sufficient for 1-el integral call
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,5,'POPML')
c
c -- set memory mark for TEXAS memory manager
      call getmem(IMem,lastx)
C
C
C  ----------------------------------------------------------------------
C
      CALL POPMAIN(NAtoms,  Z(IAN),  Z(IQA),  Z(IXC),  NBas,
     $             NShl,    Z(IBAS), Z(INX),  Z(IBA),  Z(IS),
     $             Z(ISH),  NAlpha,  NBeta,   Z(IPA),  Z(IPB),
     $             Z(IPS),  Z(IG),   Z(IBO),  NScr,    Z(IScr),
     $             IMult,   Semi,    IErr)
C
C  ----------------------------------------------------------------------
C
C  reset print flag to 1
C
      IPRNT = 1
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  free memory used in this routine
C
        call retmem(1)
cc      CALL ffree(Z(iptr))
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE POPMAIN(NAtoms, IAN,    QA,     XC,     NBas,
     $                   NShell, BASDAT, INX,    NBAtm,  S,
     $                   SHalf,  NAlpha, NBeta,  PA,     PB,
     $                   PS,     g,      BO,     NMem,   Z,
     $                   IMult,  Semi,   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for population analysis module
C  POPMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION IAN(NAtoms),QA(NAtoms),XC(3,NAtoms),BASDAT(13,*),
     $          INX(12,*),NBAtm(NAtoms),S(NBas,NBas),SHalf(NBas,NBas),
     $          PA(NBas,NBas),PB(NBas,NBas),PS(NBas,NBas),
     $          g(NBas),BO(NAtoms,NAtoms)
      DIMENSION Z(NMem)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
      CHARACTER Pop*8,cdum*20
      Character*256 jobname,MOS,MOB
      Logical Semi
c
      PARAMETER (Zero=0.0d0)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
C  initialize
C
      IOut=ioutfil('iout')
C
C set up the MOs filenames MOS and MOB
C
      call tstchval('mosfname',iyes)
      if(iyes.eq.1) then
        call getchval('mosfname',MOS)
      else
        MOS=jobname(1:lenJ)
      end if
      call rmblan2(MOS,256,lenM)
      if(lenM.le.252) MOS(lenM+1:lenM+4)='.mos'
      lenM=lenM+4
      MOB=MOS
      MOB(lenM:lenM)='b'
C
C  Read from the <control> file
C     type of population analysis
C     threshold for printing bond orders
C     print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$pthrsh',2,idum,PThrsh,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      call rdcntrl(IUnit,4,'$pop',3,idum,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
      Pop = cdum(1:8)
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XC,-1,jnk)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
C  get actual atomic charges (there may, e.g., be ghost atoms)
C
      CALL GetAtChrg(NAtoms,QA)
C ..........................................................
C
      IF(Semi) THEN
C
C  get the actual number of basis functions
C
      CALL GetNHeavy(NAtoms,IAN,NHeavy)
      CALL GetNGhosts(Natoms,IAN,QA,NHGhosts,NLGhosts)
c
      NLight = NAtoms-NHeavy-NLGhosts
      NHeavy = NHeavy-NHGhosts
      NBas = 4*NHeavy + NLight
C
C  modify the number of occupied MOs and the charge for semiempirical
C
      NBetaS = NBeta
      DO I=1,NAtoms
      II = IAN(I)
      If(II.GT.2.AND.II.LE.10) Then
       NAlpha = NAlpha-1
       NBeta = NBeta-1
       QA(I) = QA(I)-2.0d0
      Else If(II.GT.10.AND.II.LE.18) Then
       NAlpha = NAlpha-5
       NBeta = NBeta-5
       QA(I) = QA(I)-10.0d0
      Else If(II.EQ.18.OR.II.EQ.19) Then
       NAlpha = NAlpha-9
       NBeta = NBeta-9
       QA(I) = QA(I)-18.0d0
      Else If(II.GE.30.AND.II.LE.36) Then
       NAlpha = NAlpha-14
       NBeta = NBeta-14
       QA(I) = QA(I)-28.0d0
      Else If(II.EQ.37.OR.II.EQ.38) Then
       NAlpha = NAlpha-18
       NBeta = NBeta-18
       QA(I) = QA(I)-36.0d0
      Else If(II.GE.48.AND.II.LE.54) Then
       NAlpha = NAlpha-23
       NBeta = NBeta-23
       QA(I) = QA(I)-46.0d0
      EndIf
      EndDO
      If(NBetaS.EQ.0) NBeta=0
C
C  set number of basis functions per atom
C
      DO I=1,NAtoms
      If(IAN(I).GT.1) Then
        NBAtm(I) = 4
      Else
        NBAtm(I) = 1
      EndIf
      EndDO
C
C  now read the MOs (from the existing binary file)
C
      itype = 2
      CALL ReadMOS(NBas,PA,jnk,.False.,lenM,MOS(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(1,'Population module',
     $   'MOS File Does Not Exist',0,0)
      If(NBeta.GT.0) Then
        CALL ReadMOS(NBas,PB,jnk,.False.,lenM,MOB(1:lenM),itype,IErr)
        If(IErr.NE.0) Call nerror(1,'Population module',
     $     'MOB File Does Not Exist',0,0)
      EndIf
C
C  form the density matrix
C
      CALL MakeDENS(NBas,NAlpha,PA,Z)
      If(NBeta.EQ.0.AND.IMult.EQ.1) CALL VScal(NBas*(NBas+1)/2,2.0d0,Z)
      Call Expand(NBas,Z,PA)
      If(NBeta.GT.0) Then
        CALL MakeDENS(NBas,NBeta,PB,Z)
        Call Expand(NBas,Z,PB)
      EndIf
cc      write(6,*) ' Density matrix is:'
cc      call prntmat(nbas,nbas,nbas,pa)
C
C  Overlap matrix is a unit matrix for semiempirical
C
      CALL SetDiagMat(NBas,1.0d0,S)
      CALL SetDiagMat(NBas,1.0d0,SHalf)
cc
      ELSE
cc
C
C  Hartree-Fock or DFT wavefunction
C
C  Read the basis set data
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdbasis(IUnit,  NAtoms, AtSymb, XC,     Z(1),
     $             Z(1+NAtoms),BASDAT)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C ....................................................................
C  Be Careful with basis order!
C  We need basis ordered PER ATOM for population analysis
C  Wolinski ordering (for integral evaluation and hence for the MOs)
C  is PER SHELL (e.g. all S functions together regardless of atom)

C  The following routines are relevant here
C   SortBAS1  -  sorts basis per shell (but NOT full Wolinski)
C   reorder   -  orders basis from SortBAS1 into Wolinski
C   SortBAS2  -  orders basis per atom and supplies integer
C                ordering array relating atom and Wolinski order
C ....................................................................
      CALL SortBAS1(NAtoms,NShell,Z(1+NAtoms),INX)
      CALL normaliz(NShell,INX,BASDAT)

C  get number of basis functions per atom

      CALL BasATOM(NAtoms,NShell,Z(1+NAtoms),NBAtm)
cc      write(6,*) ' Number of basis functions per atom is:'
cc      do i=1,natoms
cc      write(6,*) i,nbatm(i)
cc      enddo
C
C  need to reorder MOs as need basis functions ordered per atom
C  and MOs have Wolinski special order
C
      i1 = 1
      i2 = i1 + NBas
      i3 = i2 + NShell
      IEnd = i3 + 12*NShell - 1
      CALL MemCHK(NMem,IEnd,5,'POPML')
c
      Call reorder(NShell,INX,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinski order
      Call SortBAS2(NShell,Z(i3),Z(i2),Z(i1),INX)          ! per atom
C
C  now read the old MOs (from the existing binary file)
C
      itype = 1
      CALL ReadMOS(NBas,PS,jnk,.False.,lenM,MOS(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(1,'Population module',
     $   'MOS File Does Not Exist',0,0)
      If(NBeta.GT.0) Then
        CALL ReadMOS(NBas,S,jnk,.False.,lenM,MOB(1:lenM),itype,IErr)
        If(IErr.NE.0) Call nerror(1,'Population module',
     $     'MOB File Does Not Exist',0,0)
      EndIf
C
C  need to sort MOs to per atom order
C  and form overlap and density matrices
C
      If(NBeta.EQ.0.AND.IMult.EQ.1) Then
        CALL ReorderMO(NBas,NAlpha,Z(i1),PS,PA)
        CALL FormDEN(NBas,NAlpha,PA,PS)
        Call Expand(NBas,PS,PA)
      Else
        CALL ReorderMO(NBas,NAlpha,Z(i1),PS,PA)
        CALL ReorderMO(NBas,NBeta,Z(i1),S,PB)
        CALL FormDENU(NBas,NAlpha,NBeta,PA,PB,PS,S)
        Call Expand(NBas,PS,PA)
        Call Expand(NBas,S,PB)
      EndIf
C  ...................................................................
C
C  get the overlap matrix
C
      CALL inton2(0,      NAtoms, S,      INX,    INX,
     $            0,      0,      BASDAT, BASDAT, XC,
     $            IAN,    NShell, NShell, NBas,   NBas,
     $            Z(i1))
C
C  diagonalize S
C
      CALL CpyVEC(NBas*NBas,S,PS)
      CALL DIAGMAT(PS,NBas,Z(i2),Z(i1),g,IErr)
c
      If(IErr.NE.0) Then
        Call nerror(2,'Population Analysis module',
     $     'Unable to Diagonalize Overlap Matrix!',0,0)
      EndIf
C
C  finally form S**1/2
C
      DO 10 I=1,NBas
      g(I) = SQRT(g(I))
 10   CONTINUE
c
      DO 20 I=1,NBas
      DO 20 J=1,I
      Val = Zero
      DO 19 K=1,NBas
      Val = Val + PS(I,K)*PS(J,K)*g(K)
 19   CONTINUE
      SHalf(I,J) = Val
      SHalf(J,I) = Val
 20   CONTINUE
cc      write(6,*) ' S**1/2 matrix is:'
cc      call prntmat(nbas,nbas,nbas,shalf)
cc
      ENDIF
C
C  Now do the population analysis
C
C  ------------------------------------------------------------
      CALL POPANAL(Pop,    NAtoms, AtSymb, QA,     NBas,
     $             NAlpha, NBeta,  PThrsh, IPRNT,  NBAtm,
     $             PA,     PB,     S,      SHalf,  PS,
     $             Z,      g,      BO)
C  ------------------------------------------------------------
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE POPANAL(Pop,    NAtoms, AtSymb, QA,     NBas,
     $                   NAlpha, NBeta,  PThrsh, IPRNT,  NBAtm,
     $                   PA,     PB,     S,      SHalf,  PS,
     $                   INDX,   g,      BO)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Mulliken and Lowdin Population analysis
C  Forms and prints out:
C      gross orbital occupancies (optional)
C      atomic charges
C      bond orders
C      atomic valencies
C      free atomic valencies (for unrestricted wavefunctions)
C
C  ARGUMENTS
C
C  Pop     -  type of population analysis
C  NAtoms  -  number of atoms
C  QA      -  atomic charge
C  NBas    -  number of basis functions
C  NAlpha  -  number of alpha spin/closed-shell occupied MOs
C  NBeta   -  number of beta spin occupied MOs
C  PThrsh  -  threshold for printing bond orders
C  IPRNT   -  print flag
C  NBAtm   -  number of basis functions per atom
C  PA      -  alpha/closed-shell density matrix
C  PB      -  beta density matrix
C  S       -  overlap matrix
C  SHalf   -  S**1/2
C  PS      -  array for storing (PA+PB)*S
C  INDX    -  general integer storage
C  g       -  general vector storage (for atomic valencies)
C  BO      -  bond order matrix
C
C
      DIMENSION QA(NAtoms),NBAtm(NAtoms),PA(NBas,NBas),PB(NBas,NBas),
     $          S(NBas,NBas),SHalf(NBas,NBas),PS(NBas,NBas)
      DIMENSION INDX(NAtoms),g(NBas),BO(NAtoms,NAtoms)
      CHARACTER*8 Pop,AtSymb(NAtoms)
      Logical mulliken,lowdin,full,chelp
      character*256 jobname
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
C
C
      IOut = ioutfil('iout')
      ICon = ioutfil('icon')
c
      mulliken =  Pop.EQ.'mulliken'.OR.Pop.EQ.'standard'.OR.
     $            Pop(1:4).EQ.'full'
      lowdin = Pop(1:6).EQ.'lowdin'.OR.Pop.EQ.'standard'.OR.
     $         Pop(1:4).EQ.'full'
      full = Pop(1:4).EQ.'full'
      chelp = Pop(1:5).EQ.'chelp'.OR.Pop(1:4).eq.'full'
c
c  -- form spinless density matrix if open shell
      If(NBeta.GT.0) Call AddVEC(NBas*NBas,PA,PB,PA)
C
C
      IF(mulliken) THEN
C
C  Mulliken Analysis
C  -----------------
C
      WRITE(IOut,1000)
      WRITE(ICon,1000)
C
C  Form PS matrix
C
      Call DGEMM('N',    'N',    NBas,   NBas,   NBas,
     $            One,    PA,    NBas,   S,      NBas,
     $            Zero,   PS,    NBas)
C
C  Get gross orbital occupancies
C
      If(full.OR.IPRNT.GT.2) Then
       WRITE(IOut,1100)
       Do I=1,NBas
       WRITE(IOut,1110) I,PS(I,I)
       EndDo
      EndIf
C
C  get atomic charges
C
      WRITE(IOut,1200)
      WRITE(ICon,1200)
      I1 = 0
c
      DO 15 I=1,NAtoms
      I2 = I1 + NBAtm(I)
      I1 = I1 + 1
      Val = Zero
      DO 14 II=I1,I2
      Val = Val + PS(II,II)
 14   CONTINUE
      Val = QA(I) - Val
      WRITE(IOut,1210) I,AtSymb(I),Val
      WRITE(ICon,1210) I,AtSymb(I),Val
      I1 = I2
 15   CONTINUE
C
C  Form Bond Order Matrix
C
      CALL ZeroIT(BO,NAtoms*NAtoms)
      I1 = 0
c
      DO 30 I=1,NAtoms
      I2 = I1 + NBAtm(I)
      I1 = I1 + 1
      J1 = 0
      DO 25 J=1,I-1
      J2 = J1 + NBAtm(J)
      J1 = J1 + 1
      DO 20 II=I1,I2
      DO 20 JJ=J1,J2
      BO(I,J) = BO(I,J) + PS(II,JJ)*PS(JJ,II)
 20   CONTINUE
      BO(J,I) = BO(I,J)
      J1 = J2
 25   CONTINUE
      I1 = I2
 30   CONTINUE
      BO(1,1) = Zero
C
C  print selected bond orders (greater than a threshold)
C
      WRITE(IOut,1300) PThrsh
      WRITE(ICon,1300) PThrsh
      DO 40 I=2,NAtoms
      nt = 0
      DO 39 J=1,I-1
      If(Abs(BO(J,I)).GT.PThrsh) Then
       nt = nt+1
       g(nt) = BO(J,I)
       INDX(nt) = J
      EndIF
 39   CONTINUE
      WRITE(IOut,1310) (I,INDX(J),g(J),J=1,nt)
      WRITE(ICon,1310) (I,INDX(J),g(J),J=1,nt)
 40   CONTINUE
C
C  get atomic valencies
C
      If(NBeta.EQ.0) WRITE(IOut,1400)
      If(NBeta.EQ.0) WRITE(ICon,1400)
      DO 50 I=1,NAtoms
      Val = Zero
      DO 45 J=1,NAtoms
      Val = Val + BO(J,I)
 45   CONTINUE
      g(I) = Val
      If(NBeta.EQ.0) WRITE(IOut,1210) I,AtSymb(I),Val
      If(NBeta.EQ.0) WRITE(ICon,1210) I,AtSymb(I),Val
 50   CONTINUE
C
C  for open shell get free valencies
C
      If(NBeta.GT.0) THEN
       WRITE(IOut,1410)
       WRITE(ICon,1410)
       I1 = 0
       DO 51 I=1,NAtoms
       I2 = I1 + NBAtm(I)
       I1 = I1 + 1
       Val = Zero
       DO 46 II=I1,I2
       Val = Val + PS(II,II) + PS(II,II)
       DO 46 JJ=I1,I2
       Val = Val - PS(II,JJ)*PS(JJ,II)
 46    CONTINUE
       Val = Val - g(I)
       WRITE(IOut,1420) I,AtSymb(I),g(I),Val
       WRITE(ICon,1420) I,AtSymb(I),g(I),Val
       I1 = I2
 51    CONTINUE
      ENDIF
cc
      ENDIF
C
C
      IF(lowdin) THEN
C
C  Lowdin Analysis
C  ---------------
C
      WRITE(IOut,1500)
      WRITE(ICon,1500)
C
C  form Lowdin Density Matrix  (SHalf*PA*SHalf)
C
      Call DGEMM('N',    'N',    NBas,   NBas,   NBas,
     $            One,    PA,    NBas,   SHalf,  NBas,
     $            Zero,   PS,    NBas)
c
      Call DGEMM('N',    'N',    NBas,   NBas,   NBas,
     $            One,    SHalf, NBas,   PS,     NBas,
     $            Zero,   PA,    NBas)
C
C  get gross orbital occupancies
C
      If(full.OR.IPRNT.GT.2) Then
       WRITE(IOut,1100)
       Do I=1,NBas
       WRITE(IOut,1110) I,PA(I,I)
       EndDo
      EndIf
C
C  get atomic charges
C
      WRITE(IOut,1200)
      WRITE(ICon,1200)
      I1 = 0
c
      DO 55 I=1,NAtoms
      I2 = I1 + NBAtm(I)
      I1 = I1 + 1
      Val = Zero
      DO 54 II=I1,I2
      Val = Val + PA(II,II)
 54   CONTINUE
      Val = QA(I) - Val
      WRITE(IOut,1210) I,AtSymb(I),Val
      WRITE(ICon,1210) I,AtSymb(I),Val
      I1 = I2
 55   CONTINUE
C
C  Form Bond Order Matrix
C
      CALL ZeroIT(BO,NAtoms*NAtoms)
      I1 = 0
c
      DO 70 I=1,NAtoms
      I2 = I1 + NBAtm(I)
      I1 = I1 + 1
      J1 = 0
      DO 65 J=1,I-1
      J2 = J1 + NBAtm(J)
      J1 = J1 + 1
      DO 60 II=I1,I2
      DO 60 JJ=J1,J2
      BO(I,J) = BO(I,J) + PA(II,JJ)*PA(JJ,II)
 60   CONTINUE
      BO(J,I) = BO(I,J)
      J1 = J2
 65   CONTINUE
      I1 = I2
 70   CONTINUE
      BO(1,1) = Zero
c
      WRITE(IOut,1300) PThrsh
      WRITE(ICon,1300) PThrsh
      DO 80 I=2,NAtoms
      nt = 0
      DO 79 J=1,I-1
      If(Abs(BO(J,I)).GT.PThrsh) Then
       nt = nt+1
       g(nt) = BO(J,I)
       INDX(nt) = J
      EndIF
 79   CONTINUE
      WRITE(IOut,1310) (I,INDX(J),g(J),J=1,nt)
      WRITE(ICon,1310) (I,INDX(J),g(J),J=1,nt)
 80   CONTINUE
C
C  get atomic valencies
C
      If(NBeta.EQ.0) WRITE(IOut,1400)
      If(NBeta.EQ.0) WRITE(ICon,1400)
      DO 90 I=1,NAtoms
      Val = Zero
      DO 85 J=1,NAtoms
      Val = Val + BO(J,I)
 85   CONTINUE
      g(I) = Val
      If(NBeta.EQ.0) WRITE(IOut,1210) I,AtSymb(I),Val
      If(NBeta.EQ.0) WRITE(ICon,1210) I,AtSymb(I),Val
 90   CONTINUE
C
C  for open shell get free valencies
C
      IF(NBeta.GT.0) THEN
       WRITE(IOut,1410)
       WRITE(ICon,1410)
       I1 = 0
       DO 91 I=1,NAtoms
       I2 = I1 + NBAtm(I)
       I1 = I1 + 1
       Val = Zero
       DO 86 II=I1,I2
       Val = Val + PA(II,II) + PA(II,II)
       DO 86 JJ=I1,I2
       Val = Val - PA(II,JJ)*PA(JJ,II)
 86    CONTINUE
       Val = Val - g(I)
       WRITE(IOut,1420) I,AtSymb(I),g(I),Val
       WRITE(ICon,1420) I,AtSymb(I),g(I),Val
       I1 = I2
 91    CONTINUE
      ENDIF
cc
      ENDIF
C
C  Charge from Electrostatic Potential (CHELP)
C  ---------------
      If(chelp) Then
        call getchval('jobname',jobname)
        call rmblan2(jobname,256,lenJ)
        call runchelp(jobname,lenj,NAtoms, AtSymb, IOut, ICon)
      EndIf
c
c -- atomic charges from dipole derivatives   ! GM   Sep. 2006
      call DipPop(NAtoms,AtSymb,IOut,ICon,BO,PS)
C
      RETURN
c
 1000 FORMAT(/,' ---- MULLIKEN POPULATION ANALYSIS ----')
 1100 FORMAT(/,'  gross orbital occupancies')
 1110 FORMAT(1X,I5,2X,F12.6)
 1200 FORMAT(/,18X,'atomic charges')
 1210 FORMAT(1X,I5,1X,A8,2X,F12.6)
 1300 FORMAT(/,'  Bond Orders of magnitude greater than ',F7.4)
 1310 FORMAT(1X,3(I4,'  -',I4,1X,F8.4,2X))
 1400 FORMAT(/,18X,'atomic valencies')
 1410 FORMAT(/,18X,'atomic and free valencies')
 1420 FORMAT(1X,I5,1X,A8,2(2X,F12.6))
 1500 FORMAT(/,' ---- LOWDIN POPULATION ANALYSIS ----')
c
      END
c  =======================================================================
c
      Subroutine DipPop(NAtoms,AtSymb,IOut,ICon,DipD,PolD)
      implicit real*8(A-H,O-Z)
C
C  Cioslowski population analysis (J.Am.Chem.Soc. 111 (1989) 8333)
C  Atomic charges are the trace of the dipole derivative tensor
C  These are calculated and printed if dipole derivatives are available
C
      DIMENSION DipD(3*NAtoms,3),PolD(3*NAtoms,6)
      character*8 AtSymb(NAtoms)
      logical dipole,polar
c
      call rdderiv(3*NAtoms,DipD,PolD,dipole,polar)
      If(.NOT.dipole) RETURN
c
      write(iout,1000)
      write(icon,1000)
      do iat=1,NAtoms
      chg=0.0d0
      ii = 3*(iat-1)
      do i=1,3
      chg = chg + DipD(ii+i,i)
      end do      
      chg=chg/3.0d0
      write(iout,1100) iat,AtSymb(iat),chg
      write(icon,1100) iat,AtSymb(iat),chg
      end do
c
      return
c
 1000 format(/,' ---- Atomic charges from DIPOLE DERIVATIVES ----',/)
 1100 format(1X,I5,1X,A8,2X,F12.6)
c
      end
c  =======================================================================
c
      SUBROUTINE MakeDENS(NBas,nmo,CMO,P)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  form density matrix (lower triangle) from MOs
C
      REAL*8 CMO(NBas,nmo),P(NBas*(NBas+1)/2)
c
      CALL ZeroIT(P,NBas*(NBas+1)/2)
c
      DO 30 K=1,nmo
      IJ = 0
      DO 20 I=1,NBas
      Val = CMO(I,K)
      DO 10 J=1,I
      IJ = IJ+1
      P(IJ) = P(IJ) + Val*CMO(J,K)
 10   CONTINUE
 20   CONTINUE
 30   CONTINUE
c
      RETURN
      END
