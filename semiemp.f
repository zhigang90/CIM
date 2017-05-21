c ==================================================================
c  SEMIEMPIRICAL  MODULE            JB   April 1999
c                              extended  July  2000
c ==================================================================
c
      subroutine presemi(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval,semi,jobname
c
c  reads the SEMI line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=10)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character*20 cdum
c
      parameter (IUnit=1)
c
      Data options / 'semi','nogu','ethr','dthr','lvsh','diis',
     $               'iter','pseu','uhfs','prin'/
      data ioptyp/21,0,11,11,11,11,1,11,0,1/
c
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c -- semiempirical option
      if(ifound(1).eq.1) then
        if(chopval(1)(1:1).eq.' ') then
          semi = 'pm3'
        else
          semi = chopval(1)
        endif
      endif
c put this in the depository
      call setchval('semi',semi)
c
c -- print flag
      if(ifound(10).eq.1) then
        IPRNT = iopval(1,10)
      else
        IPRNT = 2          ! default print flag
      endif
c
c *******************************************************
c *  All other options considered later in SCF routine  *
c *******************************************************
c
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call RmBlank(256,semi,Len)
      call lowerca2(semi,Len)
      call wrcntrl(IUnit,6,'$basis',3,idum,rdum,semi(1:Len))
      If(ifound(9).eq.1) Then
        call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,rdum,cdum)
        call wrcntrl(IUnit,6,'$nbeta',1,NAlpha,rdum,cdum)
        call wrcntrl(IUnit,12,'$uhf_singlet',0,idum,rdum,cdum)
      EndIf
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
      return
      end
c  =======================================================================
c
      SUBROUTINE SEMIEMP(inp,NMem,Z,natom)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** SEMIEMPIRICAL PROGRAM **
C
C  This Program calculates both the energy and gradient
C  for semiempirical wavefunctions
C
C  Hamiltonians supported:    MINDO/3    MNDO    AM1    PM3
C
C  References
C  ----------
C  MINDO/3:  "Ground States of Molecules. 25. MINDO/3. An Improved
C             Version of the MINDO Semiempirical SCF-MO Method"
C             R.C.Bingham, M.J.S.Dewar and D.H.Lo,
C             J.Am.Chem.Soc. 97 (1975) 1285
C
C     MNDO:  "Ground States of Molecules. 38. The MNDO Method
C             Approximations and Parameters"
C             M.J.S.Dewar and W.Thiel, J.Am.Chem.Soc. 99 (1977) 4899
C
C      AM1:  "AM1: A New General Purpose Quantum Mechanical
C             Molecular Model"
C             M.J.S.Dewar, E.Zoebisch, E.F.Healy and J.J.P.Stewart
C             J.Am.Chem.Soc. 107 (1985) 3902
C
C      PM3:  "Optimization of Parameters for Semiempirical Methods"
C             J.J.P.Stewart, J.Comput.Chem. 10 (1989) 209,221
C  .............................................................
C
C
      CHARACTER jobname*256,cdum*20,SEMI*20
      LOGICAL Symflag
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(natom)
c ..................................................
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
      Symflag = .FALSE.                 ! temporary
C
C  Read from the <control> file
C    number of atoms
C    number of alpha/closed-shell orbitals
C    number of beta orbitals
C    semiempirical method
C    if we have a UHF-singlet
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,rdum,cdum)
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
      Else
       Ndum1 = 0
       Ndum2 = 0
      EndIf
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,rdum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,rdum,cdum)
      call rdcntrl(IUnit,6,'$basis',3,idum,rdum,SEMI)
      call fdcntrl(IUnit,12,'$uhf_singlet',IEnd)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c ............................................................
c -- for semiempirical, dummy atoms are ignored
c
      NAtoms = NAtoms-Ndum1-Ndum2
c ............................................................
c
c -- check if uhf singlet
      If(IEnd.EQ.0) Then
        If(NBeta.EQ.0.OR.NBeta.EQ.NAlpha) Then
          NBeta = NAlpha
          call setival('NBeta',NBeta)     ! replace in depository
        Else
          Call nerror(1,'SEMIEMPIRICAL module',
     $         'UHF Singlet requested and System is Not a Singlet!',
     $          0,0)
        EndIf
      EndIf
C
C  Determine
C    number of heavy atoms
C    number of light (i.e., hydrogen) atoms
C    number of basis functions
C    total number of two-electron integrals
C
C  Read the atomic symbols and charges
C
      IXC = 1
      ICHG = IXC + 3*NAtoms
      IXM = ICHG + NAtoms
      IAN = IXM + NAtoms
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoordF(IUnit,  NAtoms, AtSymb, Z(IXC), -1,
     $              jnk,    Z(ICHG),Z(IXM))
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,Z(IAN))
C
C  Get number of heavy atoms
C
      CALL GetNHeavy(NAtoms,Z(IAN),NHeavy)
C
C  Get the number of ghost atoms
C
      CALL GetNGhosts(Natoms,Z(IAN),Z(ICHG),NHGhosts,NLGhosts)
c
      NLight = NAtoms-NHeavy-NLGhosts
      NHeavy = NHeavy-NHGhosts
      NBas = 4*NHeavy + NLight
      N2Elec = 50*NHeavy*(NHeavy-1) + 10*NHeavy*NLight
     $           + (NLight*(NLight-1))/2
      NLin = (NBas*(NBas+1))/2
c ..............................................................
C
C  determine amount of "scratch" memory
C
      MDiis = 6                  ! maximum size of Diis subspace
c
      NScr = 2*NLin
C
C  Now get the memory
C
      IMem = 13*NAtoms + 5*NBas + N2Elec + 5*NLin + NBas**2
     $         + (MDiis+1)**2
      If(NBeta.GT.0) IMem = IMem + NBas + 3*NLin + NBas**2
     $                        + (MDiis+1)**2
      If(Symflag) IMem = IMem + 9*NTrans + NAtoms*NTrans
      IMem = IMem + NScr
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,4,'SEMI')
cc      CALL ZeroIT(Z,IMem)             ! clear the memory
C
C
C  Allocate memory pointers
C
      IAN   = iptr                    !  atomic numbers
      IUQ   = IAN  + NAtoms           !  symmetry-unique atoms
      IXC   = IUQ  + NAtoms           !  geometry
      ICHG  = IXC  + 3*NAtoms         !  atomic charges
      IXM   = ICHG + NAtoms           !  atomic masses
      INF   = IXM  + NAtoms           !  NFIRST array
      INM   = INF  + NAtoms           !  NMIDLE array
      INL   = INM  + NAtoms           !  NLAST array
      IFAC  = INL  + NAtoms           !  triangular numbers  I*(I-1)/2
      IFC1  = IFAC + NBas             !  triangular numbers  I*(I+1)/2
      IUSP  = IFC1 + NBas             !  USPD array
      IPSP  = IUSP + NBas             !  PSPD array
      IW    = IPSP + NBas             !  two-electron integrals
      IH    = IW   + N2Elec           !  one-electron matrix
      IP    = IH   + NLin             !  total density matrix
      IF    = IP   + NLin             !  alpha/closed-shell Fock matrix
      IC    = IF   + NLin             !    ditto  MOs
      IEA   = IC   + NBas*NBas        !    ditto  orbital energies
      IPA   = IEA  + NBas             !    ditto  density matrix
      IPAO  = IPA  + NLin             !    ditto  old density matrix
      IBA   = IPAO + NLin             !    ditto  DIIS error matrix
      IEnd  = IBA  + (MDiis+1)**2
      If(NBeta.GT.0) Then
       IFB  = IEnd                    !  beta Fock matrix
       ICB  = IFB  + NLin             !    ditto  MOs
       IEB  = ICB  + NBas*NBas        !    ditto  orbital energies
       IPB  = IEB  + NBas             !    ditto  density matrix
       IPBO = IPB  + NLin             !    ditto  old density matrix
       IBB  = IPBO + NLin             !    ditto  DIIS error matrix
       IEnd = IBB  + (MDiis+1)**2
      Else
       IFB  = 1
       ICB  = 1
       IEB  = 1
       IPB  = 1
       IPBO = 1
       IBB  = 1
      EndIf
C
C  storage for symmetry
C
      IF(Symflag) THEN
c
       ITN = IEnd                     !  symmetry operations as 3x3 matrices
       INQ = ITN + 9*NTrans           !  list of atomic equivalences
       IEnd = INQ + NAtoms*NTrans
      ELSE
       ITN = 1
       INQ = 1
c
      ENDIF
c
      IGC  = IEnd                     !  gradient
      IEnd = IGC  + 3*NAtoms
C
C  general scratch storage
C
      IScr = IEnd
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,7,'SEMIEMP')
c
c -- set memory mark for TEXAS memory manager
      call getmem(IMem,lastx)
C
C
C  ----------------------------------------------------------------------
C
      CALL SEMIMAIN(inp,     NAtoms,  SEMI,    AtSymb,  Z(IAN),
     $              NQ,      Z(IUQ),  Z(IXC),  Z(ICHG), Z(IXM),
     $              NBas,    NAlpha,  NBeta,   N2Elec,  MDiis,
     $              Z(INF),  Z(INM),  Z(INL),  Z(IFAC), Z(IFC1),
     $              Z(IUSP), Z(IPSP), Z(IW),   Z(IH),   Z(IP),
     $              Z(IF),   Z(IC),   Z(IEA),  Z(IPA),  Z(IPAO),
     $              Z(IBA),  Z(IFB),  Z(ICB),  Z(IEB),  Z(IPB),
     $              Z(IPBO), Z(IBB),  NTrans,  Z(ITN),  Z(INQ),
     $              Z(IGC),  NScr,    Z(IScr), IErr)
C
C  ----------------------------------------------------------------------
C  free memory used in this routine
C
      call retmem(1)
C
C  Exit procedure
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE SEMIMAIN(inp,    NAtoms, SEMI,   AtSymb, IAN,
     $                    NQ,     IUNQ,   XC,     XCharg, XMass,
     $                    NBas,   NAlpha, NBeta,  N2Elec, MDiis,
     $                    NFIRST, NMIDLE, NLAST,  IFACT,  I1FACT,
     $                    USPD,   PSPD,   W,      H,      P,
     $                    F,      C,      EIGA,   PA,     POLD,
     $                    BA,     FB,     CB,     EIGB,   PB,
     $                    PBOLD,  BB,     NTrans, TRANS,  NEqATM,
     $                    GC,     NMem,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main driver for semiempirical energy + gradient
C
C  ARGUMENTS
C
C  inp     -  unit number for input file
C  NAtoms  -  number of atoms
C  SEMI    -  which semiempirical method
C  AtSymb  -  atomic symbols
C  IAN     -  atomic numbers
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of symmetry unique atoms in ascending order
C  XC      -  Cartesian coordinates
C  XCharg  -  atomic charges
C  XMass   -  atomic masses
C  NBas    -  number of basis functions (on input only upper bound)
C  NAlpha  -  number of occupied alpha/closed shell orbitals
C  NBeta   -  number of occupied beta orbitals
C  N2Elec  -  number of 2-electron repulsion integrals
C  MDiis   -  maximum size of DIIS subspace
C  NFIRST
C  NMIDLE
C  NLAST
C  IFACT   -  triangular numbers   I*(I-1)/2
C  I1FACT  -   ditto               I*(I+1)/2
C  USPD    -  one-electron diagonal terms
C  PSPD    -  two-electron diagonal terms
C  W       -  two-electron matrix
C  H       -  one-electron matrix
C  P       -  total density matrix
C  F       -  alpha or closed-shell Fock matrix
C  C       -    ditto molecular orbitals
C  EIGA    -    ditto orbital energies
C  PA      -    ditto density matrix
C  POLD    -    ditto old density matrix
C  BA      -  alpha or closed-shell DIIS error matrix
C  FB      -  beta Fock matrix
C  CB      -    ditto molecular orbitals
C  EIGB    -    ditto orbital energies
C  PB      -    ditto density matrix
C  PBOLD   -    ditto old density matrix
C  BB      -  beta DIIS error matrix (?)
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C             (currently symmetry is not used)
C  GC      -  Cartesian gradient
C  NMem    -  available scratch memory
C  Z       -  scratch storage
C  IErr    -  error flag
C
C
      DIMENSION IAN(NAtoms),IUNQ(NQ),XC(3,NAtoms),XCharg(NAtoms),
     $          XMass(NAtoms),GC(3,NAtoms)
      DIMENSION NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms),
     $          IFACT(NBas),I1FACT(NBas)
      DIMENSION USPD(NBas),PSPD(NBas),W(N2Elec),C(NBas,NBas),
     $          H(*),P(*),F(*),PA(*),FB(*),PB(*),POLD(*),PBOLD(*),
     $          EIGA(NBas),EIGB(NBas),CB(NBas,NBas),BA(*),BB(*)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      CHARACTER*8 AtSymb(NAtoms)
      Character*20 SEMI
      DIMENSION Z(NMem)
      Dimension DIP(3)
      Logical precise
c
      Character*256 scrf,jobname,MOS,MOB
      Character char*35,cdum*20,wvfnc*20
      Common /job/jobname,lenJ
c
      COMMON /CORES / CORE(54)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,clight,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,RGas
      PARAMETER (EVKCAL=23.0605423d0,CVRT=0.0015936d0)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      Data JUnit/15/, KUnit/16/         ! unit number for direct access files
C
C
      AUTOAN = 1.0d0/ANTOAU
c      CVRT=rgetrval('kcal/mol')
c      EVKCAL=RGETRVAL('evol')/CVRT
C
C  Read from the <control> file
C     total molecular charge
C     print flag
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$charge',2,idum,Charg,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoordF(IUnit, NAtoms, AtSymb, XC,     -1,
     $              jnk,   XCharg, XMass)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c -- convert to angstroms
      CALL VScal(3*NAtoms,AUTOAN,XC)
C
C  ...................................................................
C  Open DIIS direct access files
C
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,lnscr)
c
      LenRec = 8*NBas*(NBas+1)
c
      OPEN (UNIT=JUnit,FILE=scrf(1:lnscr)//'.da',
     $      ACCESS='DIRECT',RECL=LenRec,FORM='UNFORMATTED',
     $      STATUS='UNKNOWN')
c
      If(NBeta.GT.0)
     $ OPEN (UNIT=KUnit,FILE=scrf(1:lnscr)//'.db',
     $       ACCESS='DIRECT',RECL=LenRec,FORM='UNFORMATTED',
     $       STATUS='UNKNOWN')
C  ....................................................................
c
      IOut=ioutfil('iout')
      wvfnc = 'Semiemp-'//SEMI(1:5)
c
      If(IPRNT.GT.1) Then
        WRITE(IOut,1000)
        If(NAlpha.EQ.NBeta) WRITE(IOut,1100)
      EndIf
C
C  Get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
c -- remove ghost atoms from the semiempirical module
      call RemoveGhosts(NAtoms,AtSymb,IAN,XCharg,XC)
C
C  Can we do the calculation?
C
      CALL ChkSEMI(SEMI,NAtoms,IAN,NAlpha,NBeta,IErr)
c
      If(IErr.NE.0) Then
        char = SEMI(1:5)//' Not Possible for this System'
        Call nerror(3,'SEMIEMPIRICAL module',Char,0,0)
      EndIf
C
C  set up semiempirical common blocks
C
      ICharg = NINT(Charg)
      CALL MOLDAT(NAtoms, NBas,   IAN,    XC,     ICharg,
     $            SEMI,   NFIRST, NMIDLE, NLAST,  USPD,
     $            PSPD,   ATHEAT)
C
C  set up triangular number arrays
C
      DO 10 I=1,NBas
      IFACT(I) = (I*(I-1))/2
      I1FACT(I) = IFACT(I) + I
 10   CONTINUE
C
C  form the one-electron matrix and two-electron integrals
C
      CALL HCORE(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $           NFIRST, NMIDLE, NLAST,  USPD,   H,
     $           W,      ENuclr)
C
C  do the SCF
C
      CALL ITER(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $          NFIRST, NMIDLE, NLAST,  IFACT,  I1FACT,
     $          NBas,   PSPD,   NAlpha, NBeta,  MDiis,
     $          H,      W,      P,      F,      C,
     $          EIGA,   PA,     POLD,   BA,     FB,
     $          CB,     EIGB,   PB,     PBOLD,  BB,
     $          ATHEAT, ENuclr, inp,    Z,      EE,     IErr)
c
      If(IErr.NE.0) Then
        Call nerror(4,'SEMIEMPIRICAL module',
     $        'SCF Failed to Converge',0,0)
      EndIF
C
C  form total energy
C
      E = (EE + ENuclr)*EVKCAL + ATHEAT
c
      If(IPRNT.GT.1) WRITE(IOut,1200) SEMI,E
C
C  calculate the dipole moment
C
      CALL CHRGE(NAtoms, NFIRST, NLAST,  P,      Z)
c
      DO 20 I=1,NAtoms
      II = IAN(I)
      Z(I) = CORE(II) - Z(I)
 20   CONTINUE
c
      DipM = DIPOLE(NAtoms, IAN,    XC,     SEMI,   NFIRST,
     $              P,      Z,   Z(1+NAtoms), DIP)
c
      If(IPRNT.GT.1) WRITE(IOut,1300) DIP(1),DIP(2),DIP(3),DipM
C
C  calculate the gradient
C
      precise = .TRUE.
c
      CALL DCART(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $           NFIRST, NMIDLE, NLAST,  IFACT,  I1FACT,
     $           NBas,   NAlpha, NBeta,  P,      PA,
     $           PB,    precise, GC)
C
C  convert energy/gradient to atomic units
C
      E = E*CVRT
      CALL VScal(3*NAtoms,CVRT*AUTOAN,GC)
C
c ==========================================================
c -- write final energy and dipole moment to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,0,rdum,wvfnc)
      Call wrcntrl(IUnit,7,'$energy',2,0,E,cdum)
      Call wrcntrl(IUnit,5,'$escf',2,0,E,cdum)      ! tentative  JB May 2011
      Call WrDIP(IUnit,DIP)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c -- write forces to <grad> file
      CALL WrGRAD(NAtoms,GC)
c ==========================================================
C
C  Write converged MOs to MOS file
C
      nwmo = MIN(NAlpha+20,NBas)     ! number of MOS to save in MOS file
      call tstchval('mos-file',iyes)
      If(iyes.eq.1) Then
        call getchval('mos-file',MOS)
        If(NBeta.GT.0) call getchval('mob-file',MOB)
      Else
        MOS = jobname(1:lenJ)//'.mos'
        MOB = jobname(1:lenJ)//'.mob'
      EndIf
      call rmblan2(MOS,256,lenM)
      call WriteMOS(NBas,nwmo,C,EIGA,.True.,lenM,MOS,2)
      If(NBeta.GT.0) Then
        call rmblan2(MOB,256,lenM)
        call WriteMOS(NBas,nwmo,CB,EIGB,.True.,lenM,MOB,2)
      EndIf
C
C  Close DIIS direct access files
C
      CLOSE (UNIT=JUnit,STATUS='DELETE')
      CLOSE (UNIT=KUnit,STATUS='DELETE')
C
      RETURN
c
 1000 FORMAT(/,' ** SEMIEMPIRICAL MODULE **')
 1100 FORMAT(/,' System is an unrestricted open-shell singlet')
 1200 FORMAT(/,1X,A6,'  Heat of Formation: ',F18.9,' Kcal/mol')
 1300 FORMAT(/,' Dipole/au = ',3F10.6,' Total: ',F10.6)
c
      END
c =======================================================================
c
      SUBROUTINE ChkSEMI(SEMI,NAtoms,IAN,NAlpha,NBeta,IErr)
      IMPLICIT INTEGER(A-Z)
C
C  Checks to see if all atoms in the molecule have parameters for
C  the semiempirical model requested. Also modifies the number of
C  occupied orbitals to eliminate the core.
C
C  ARGUMENTS
C
C  SEMI    -  which semiempirical method
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  NAlpha  -  on entry number of alpha/closed-shell MOs
C             on exit ditto for semiempirical
C  NBeta   -  on entry number of beta MOs
C             on exit ditto for semiempirical
C  IErr    -  Error flag on exit
C               0 - all atoms have semiempirical parameters
C              -1 - 1 or more atoms DO NOT have parameters
C
C
C .......................................................
C   ELEMENTS THAT HAVE BEEN PARAMETRIZED ARE:
C             PM3   AM1  MNDO  MINDO
C   1   H      Y     Y     Y     Y
C   3   Li     Y     -     Y     -
C   4   Be     Y     Y     Y     -
C   5   B      Y     Y     Y     Y
C   6   C      Y     Y     Y     Y
C   7   N      Y     Y     Y     Y
C   8   O      Y     Y     Y     Y
C   9   F      Y     Y     Y     Y
C  11   Na     Y     -     -     -
C  12   Mg     Y     -     -     -
C  13   Al     Y     Y     Y     -
C  14   Si     Y     Y     Y     Y
C  15   P      Y     Y     Y     Y
C  16   S      Y     Y     Y     Y
C  17   Cl     Y     Y     Y     Y
C  19   K      Y     -     -     -
C  20   Ca     Y     -     -     -
C  30   Zn     Y     Y     Y     -
C  31   Ga     Y     -     -     -
C  32   Ge     Y     Y     Y     -
C  33   As     Y     -     -     -
C  34   Se     Y     -     -     -
C  35   Br     Y     Y     Y     -
C  37   Rb     Y     -     -     -
C  38   Sr     Y     -     -     -
C  48   Cd     Y     -     -     -
C  49   In     Y     -     -     -
C  50   Sn     Y     Y     Y     -
C  51   Sb     Y     -     -     -
C  52   Te     Y     -     -     -
C  53   I      Y     Y     Y     -
C ........................................................
C
      DIMENSION IAN(NAtoms)
      Character*20 SEMI
C
      IErr = -1
      NBetaS = NBeta
c
      IF(SEMI(1:3).EQ.'pm3') THEN
        DO 10 I=1,NAtoms
        II = IAN(I)
        If(II.EQ.2.OR.II.EQ.10.OR.II.EQ.18.OR.II.EQ.36.OR.
     $     II.GE.54) RETURN
        If( (II.GT.20.AND.II.LT.30) .OR.
     $      (II.GT.38.AND.II.LT.48) ) RETURN
c
        If(II.GT.2.AND.II.LE.10) Then
         NAlpha = NAlpha-1
         NBeta = NBeta-1
        Else If(II.GT.10.AND.II.LE.18) Then
         NAlpha = NAlpha-5
         NBeta = NBeta-5
        Else If(II.EQ.18.OR.II.EQ.19) Then
         NAlpha = NAlpha-9
         NBeta = NBeta-9
        Else If(II.GE.30.AND.II.LE.36) Then
         NAlpha = NAlpha-14
         NBeta = NBeta-14
        Else If(II.EQ.37.OR.II.EQ.38) Then
         NAlpha = NAlpha-18
         NBeta = NBeta-18
        Else If(II.GE.48.AND.II.LE.54) Then
         NAlpha = NAlpha-23
         NBeta = NBeta-23
        EndIf
 10     CONTINUE
      ELSE IF(SEMI(1:5).EQ.'mindo') THEN
        DO 20 I=1,NAtoms
        II = IAN(I)
        If(II.EQ.2.OR.II.EQ.3.OR.II.EQ.4.OR.II.EQ.10.OR.
     $     II.EQ.11.OR.II.EQ.12.OR.II.EQ.13.OR.II.GE.18) RETURN
c
        If(II.GT.2) Then
         NAlpha = NAlpha-1
         NBeta = NBeta-1
        EndIf
        If(II.GT.10) Then
         NAlpha = NAlpha-4
         NBeta = NBeta-4
        EndIf
 20     CONTINUE
      ELSE
        DO 30 I=1,NAtoms
        II = IAN(I)
c -- elements for MNDO/AM1
        If(II.EQ.2.OR.II.EQ.10.OR.II.EQ.11.OR.II.EQ.12) RETURN
        If(II.GT.17.AND.(II.NE.30.AND.II.NE.32.AND.II.NE.35.AND.
     $                   II.NE.50.AND.II.NE.53) ) RETURN
c -- for AM1, Li is not parametrized
        If(SEMI(1:3).EQ.'am1'.AND.II.EQ.3) RETURN
c
        If(II.GT.2.AND.II.LE.10) Then
         NAlpha = NAlpha-1
         NBeta = NBeta-1
        Else If(II.GT.10.AND.II.LE.18) Then
         NAlpha = NAlpha-5
         NBeta = NBeta-5
        Else If(II.EQ.18.OR.II.EQ.19) Then
         NAlpha = NAlpha-9
         NBeta = NBeta-9
        Else If(II.GE.30.AND.II.LE.36) Then
         NAlpha = NAlpha-14
         NBeta = NBeta-14
        Else If(II.EQ.37.OR.II.EQ.38) Then
         NAlpha = NAlpha-18
         NBeta = NBeta-18
        Else If(II.GE.48.AND.II.LE.54) Then
         NAlpha = NAlpha-23
         NBeta = NBeta-23
        EndIf
 30     CONTINUE
      ENDIF
c
      If(NBetaS.EQ.0) NBeta=0
      IErr = 0
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE GetNHeavy(NAtoms,IAN,NHeavy)
      IMPLICIT INTEGER(A-Z)
C
C  gets number of heavy (i.e., non-hydrogen) atoms
C  in the current system
C
      DIMENSION IAN(NAtoms)
C
      NHeavy = 0
      DO 10 I=1,NAtoms
c  PP changes to 2 - He is also light!!!
      If(IAN(I).GT.2) NHeavy = NHeavy+1
 10   CONTINUE
C
      END
c =======================================================================
c
      SUBROUTINE GetNGhosts(NAtoms,IAN,XCharg,NHGhosts,NLGhosts)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  ARGUMENTS
C
C  NAtoms  -  Total number of all atoms (ghosts + non-ghosts)
C  IAN     -  atomic numbers
C  XCharg  -  atomic charges
C  NHGhosts - on exit number of heavy (non-hydrogen) ghost atoms
C  NLGhosts - on exit number of light (hydrogen, helium) ghost atoms
C
      Dimension IAN(NAtoms),XCharg(NAtoms)
C
      NHGhosts = 0
      NLGhosts = 0
      do I=1,NAtoms
      If(XCharg(I).LT.0.99d0) Then
        If(IAN(I).gt.2) then
          NHGhosts=NHGhosts+1
        Else
          NLGhosts=NLGhosts+1
        EndIF
      EndIf
      end do
C
      END
c =======================================================================
c
      SUBROUTINE RemoveGhosts(NAtoms,AtSymb,IAN,XCharg,XC)
      IMPLICIT  real*8 (A-H,O-Z)
C
C  Removes the ghost (zero charge) atoms from the semiempirical calculation
C
C  ARGUMENTS
C
C  NAtoms  -  on input, total number of atoms
C             on exit,  number of real (non-ghost) atoms
C  AtSymb  -  on input, all atomic symbols
C             on exit,  replaced by list of non-ghost atomic symbols
C  IAN     -  atomic numbers, replaced by condensed list on exit
C  XCharg  -  atomic charges
C             (in particular, zero for ghost atoms)
C  XC      -  Cartesian coordinates, replaced by the condensed list on exit
C
      Character*8 AtSymb(NAtoms)
      DIMENSION IAN(NAtoms),XCharg(NAtoms),XC(3,NAtoms)
C
      nat=0
      do i=1,NAtoms
      If(XCharg(i).gt.0.99d0) Then
        nat=nat+1
        AtSymb(nat)=AtSymb(i)
        IAN(nat)=IAN(i)
        XC(1,nat)=XC(1,i)
        XC(2,nat)=XC(2,i)
        XC(3,nat)=XC(3,i)
      EndIf
      end do
c
      NAtoms = nat
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE CHRGE(NAtoms, NFIRST, NLAST,  P,      Q)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Stores in Q the total electron density on the atoms
C
      DIMENSION NFIRST(NAtoms),NLAST(NAtoms)
      DIMENSION P(*),Q(NAtoms)
C
      K=0
      DO 10 I=1,NAtoms
      IA=NFIRST(I)
      IB=NLAST(I)
      Q(I)=0.0d0
      DO 10 J=IA,IB
      K=K+J
      Q(I)=Q(I)+P(K)
 10   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE CNVG(NBas, PNEW, P, P1, NITER, PL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Two-point interpolation scheme for speeding convergence
C  of the density matrix
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  PNEW    -  current density matrix
C  P       -  on entry previous density matrix
C             on exit new density matrix
C  P1      -  storage for diagonal of old density matrix
C  NITER   -  SCF itertion counter
C  PL      -  on exit largest difference between old
C             and new density matrix
C
C
      DIMENSION P1(*), P(*), PNEW(*)
      LOGICAL EXTRAP
C
      PARAMETER (Zero=0.0d0)
C
      PL = Zero
      FACA = Zero
      DAMP = 1.D10
      If(NITER.GT.3) DAMP=0.05D0
      FACB = Zero
      FAC = Zero
      II=MOD(NITER,3)
      EXTRAP=II.NE.0
      K=0
      DO 30 I=1,NBas
      K=K+I
      A=PNEW(K)
      SA=ABS(A-P(K))
      If(SA.GT.PL) PL=SA
      If(EXTRAP) GO TO 20
      FACA=FACA+SA**2
      FACB=FACB+(A-2.D00*P(K)+P1(I))**2
   20 P1(I)=P(K)
   30 P(K)=A
      If (FACB.LE.Zero) GO TO 40
      If (FACA.LT.(100.D00*FACB)) FAC=SQRT(FACA/FACB)
   40 IE=0
      DO 80 I=1,NBas
      II=I-1
      DO 60 J=1,II
      IE=IE+1
      A=PNEW(IE)
      P(IE)=A+FAC*(A-P(IE))
      PNEW(IE)=P(IE)
   60 CONTINUE
      IE=IE+1
      IF(ABS(P(IE)-P1(I)) .GT. DAMP) THEN
        P(IE)=P1(I)+SIGN(DAMP,P(IE)-P1(I))
      ELSE
        P(IE)=P(IE)+FAC*(P(IE)-P1(I))
      ENDIF
      PNEW(IE)=P(IE)
  80  CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE DCART(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $                 NFIRST, NMIDLE, NLAST,  IFACT,  I1FACT,
     $                 NBas,   NAlpha, NBeta,  P,      PA,
     $                 PB,    precise, GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main driving routine for Cartesian first derivatives
C  This is done by finite differences
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C  IPRNT   -  print flag
C  NFIRST
C  NMIDLE
C  NLAST
C  IFACT   -  triangular numbers   I*(I-1)/2
C  I1FACT  -   ditto               I*(I+1)/2
C  NBas    -  number of basis functions
C  NAlpha  -  number of occupied alpha or closed-shell orbitals
C  NBeta   -  number of occupied beta orbitals
C  P       -  total density matrix
C  PA      -  alpha density matrix
C  PB      -  beta density matrix
C  precise -  logical flag for accurate derivative
C
C  on exit
C
C  GC      -  Cartesian gradient
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),GC(3,NAtoms)
      DIMENSION NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms),
     $          IFACT(NBas),I1FACT(NBas)
      DIMENSION P(*),PA(*),PB(*)
      Dimension PDI(171),PADI(171),PBDI(171),CDI(3,2),NDI(2)
      LOGICAL UHF,precise
      Character*20 SEMI
C
      PARAMETER (EVKCAL=23.0605423d0)
      DATA CHNGE,CHNGE2 /1.D-6,5.D-7/
C
c       EVKCAL=RGETRVAL('evol')/RGETRVAL('kcal/mol')
C
      UHF = NBeta.GT.0
      CALL ZeroIT(GC,3*NAtoms)
c
      DO 1 II=2,NAtoms
      IF=NFIRST(II)
      IM=NMIDLE(II)
      IL=NLAST(II)
      NDI(2)=IAN(II)
      DO 6 I=1,3
   6  CDI(I,2)=XC(I,II)
c
      DO 1 JJ=1,II-1
C
C  form diatomic matrices
C
      JF=NFIRST(JJ)
      JM=NMIDLE(JJ)
      JL=NLAST(JJ)
c
c -- first atom
      NDI(1)=IAN(JJ)
      DO 7 I=1,3
   7  CDI(I,1)=XC(I,JJ)
      IJ=0
      DO 2 I=JF,JL
      K=I*(I-1)/2+JF-1
      DO 2 J=JF,I
      IJ=IJ+1
      K=K+1
      PADI(IJ)=PA(K)
      PBDI(IJ)=PB(K)
   2  PDI(IJ)=P(K)
C
C  get second atom-first atom intersection
C
      DO 3 I=IF,IL
      L=I*(I-1)/2
      K=L+JF-1
      DO 4 J=JF,JL
      IJ=IJ+1
      K=K+1
      PADI(IJ)=PA(K)
      PBDI(IJ)=PB(K)
   4  PDI(IJ)=P(K)
      K=L+IF-1
      DO 5 L=IF,I
      K=K+1
      IJ=IJ+1
      PADI(IJ)=PA(K)
      PBDI(IJ)=PB(K)
   5  PDI(IJ)=P(K)
   3  CONTINUE
      IIJJ=IIJJ+1
      IF(.NOT.precise) THEN
        CDI(1,1)=CDI(1,1)+CHNGE2
        CDI(2,1)=CDI(2,1)+CHNGE2
        CDI(3,1)=CDI(3,1)+CHNGE2
        CALL DHC(NAtoms,SEMI,PDI,PADI,PBDI,CDI,NDI,JF,JM,JL,IF,IM,IL,
     +           UHF,NBas,IFACT,I1FACT,AA)
      ENDIF
      DO 8 K=1,3
      IF(precise)THEN
        CDI(K,2)=CDI(K,2)-CHNGE2
        CALL DHC(NAtoms,SEMI,PDI,PADI,PBDI,CDI,NDI,JF,JM,JL,IF,IM,IL,
     +           UHF,NBas,IFACT,I1FACT,AA)
      ENDIF
      CDI(K,2)=CDI(K,2)+CHNGE
      CALL DHC(NAtoms,SEMI,PDI,PADI,PBDI,CDI,NDI,JF,JM,JL,IF,IM,IL,
     +         UHF,NBas,IFACT,I1FACT,EE)
      CDI(K,2)=CDI(K,2)-CHNGE2
      If(.NOT.precise) CDI(K,2)=CDI(K,2)-CHNGE2
C +++++++++++ MODIFICACIO DE LES DERIVADES CARTESIANES  ++++++++++
C              DERIV=(AA-EE)*46.122D0/CHNGE
C +++++++++++ VENEN MULTIPLICADES PER UN FACTOR DE 2.0  ++++++++++
C +++++++++++ QUE ES INCORRECTE.                        ++++++++++
      DERIV=(AA-EE)*EVKCAL/CHNGE
      GC(K,II)=GC(K,II)+DERIV
      GC(K,JJ)=GC(K,JJ)-DERIV
   8  CONTINUE
   1  CONTINUE
c
      If(IPRNT.GT.2) Then
       WRITE(6,'(//12X,''Cartesian Gradient (Kcal/mol/A)'',//3X,
     +   ''ATOM  AT. NO.'',5X,''X'',12X,''Y'',12X,''Z'',/)')
       WRITE(6,'(2I6,F13.6,2F13.6)')
     +          (I,IAN(I),(GC(J,I),J=1,3),I=1,NAtoms)
      EndIf
C
      RETURN
      END
c ----------------------------------------------------------------------
c
      SUBROUTINE DHC(NAtoms,SEMI,P,PA,PB,XI,IAN,IF,IM,IL,JF,JM,JL,
     +               UHF,NBas,IFACT,I1FACT,DENER)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates the energy contributions from those pairs
C  of atoms that have been moved
C
C
      DIMENSION P(*), PA(*), PB(*)
      DIMENSION XI(3,*),NFIRST(2),NMIDLE(2),NLAST(2),IAN(*),
     $          IFACT(NBas),I1FACT(NBas)
      COMMON /ONELEC/ USS(54),UPP(54),UDD(54)
      LOGICAL UHF
      DIMENSION DCL(3), DIF(3), H(171), SHMAT(9,9), F(171),
     +          W(100), E1B(10), E2A(10)
      Character*20 SEMI
C
C  initialize
C
      NFIRST(1)=1
      NMIDLE(1)=IM-IF+1
      NLAST(1)=IL-IF+1
      NFIRST(2)=NLAST(1)+1
      NMIDLE(2)=NFIRST(2)+JM-JF
      NLAST(2)=NFIRST(2)+JL-JF
cc      WRITE(6,'(6I4)')(NFIRST(I),NMIDLE(I),NLAST(I),I=1,2)
      NLin=(NLAST(2)*(NLAST(2)+1))/2
      DO 10 I=1,NLin
      F(I)=0.0D0
 10   H(I)=0.0D0
c
      DO 56 I=1,2
      NI=IAN(I)
      J=NFIRST(I)
      H((J*(J+1))/2)=USS(NI)
      H(((J+1)*(J+2))/2)=UPP(NI)
      H(((J+2)*(J+3))/2)=UPP(NI)
      H(((J+3)*(J+4))/2)=UPP(NI)
      H(((J+4)*(J+5))/2)=UDD(NI)
      H(((J+5)*(J+6))/2)=UDD(NI)
      H(((J+6)*(J+7))/2)=UDD(NI)
      H(((J+7)*(J+8))/2)=UDD(NI)
      H(((J+8)*(J+9))/2)=UDD(NI)
      H(((J+9)*(J+10))/2)=UDD(NI)
  56  CONTINUE
      DO 57 I=1,NLin
  57  F(I)=H(I)
      JA=NFIRST(2)
      JB=NLAST(2)
      JC=NMIDLE(2)
      IA=NFIRST(1)
      IB=NLAST(1)
      IC=NMIDLE(1)
      JT=JB*(JB+1)/2
      J=2
      I=1
      NJ=IAN(2)
      NI=IAN(1)
      CALL H1ELEC(NI,NJ,XI(1,1),XI(1,2),SEMI,SHMAT)
      J1=0
      DO 1 J=JA,JB
      JJ=J*(J-1)/2
      J1=J1+1
      I1=0
      DO 1 I=IA,IB
      JJ=JJ+1
      I1=I1+1
      H(JJ)=SHMAT(I1,J1)
 1    F(JJ)=SHMAT(I1,J1)
      CALL ROTAT(NJ,     NI,     XI(1,2),XI(1,1),SEMI,
     $           W,      KR,     E2A,    E1B,    ENUCLR)
C
C  ENUCLR is summed over core-core repulsion integrals
C
      I2=0
      DO 7 I1=IA,IC
      II=I1*(I1-1)/2+IA-1
      DO 7 J1=IA,I1
      II=II+1
      I2=I2+1
      H(II)=H(II)+E1B(I2)
 7    F(II)=F(II)+E1B(I2)
      DO  17 I1=IC+1,IB
      II=(I1*(I1+1))/2
      F(II)=F(II)+E1B(1)
 17   H(II)=H(II)+E1B(1)
c
      I2=0
      DO 8 I1=JA,JC
      II=I1*(I1-1)/2+JA-1
      DO 8 J1=JA,I1
      II=II+1
      I2=I2+1
      H(II)=H(II)+E2A(I2)
 8    F(II)=F(II)+E2A(I2)
      DO 18 I1=JC+1,JB
      II=(I1*(I1+1))/2
      F(II)=F(II)+E2A(1)
 18   H(II)=H(II)+E2A(1)
      CALL FOCK2(2,      SEMI,   NFIRST, NMIDLE, NLAST,
     $           IFACT,  I1FACT, W,      P,      PA,   F)
      EE=HELECT(NLAST(2),PA,H,F)
c
      IF(UHF) THEN
        DO 9 I=1,NLin
 9      F(I)=H(I)
        CALL FOCK2(2,      SEMI,   NFIRST, NMIDLE, NLAST,
     $             IFACT,  I1FACT, W,      P,      PB,  F)
        EE=EE+HELECT(NLAST(2),PB,H,F)
      ELSE
        EE=EE+EE
      ENDIF
cc      WRITE(6,'('' EE AND ENUCLR IN DHC'',2F13.6)')EE, ENUCLR
      DENER=EE+ENUCLR
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE DENSIT(NBas, NOcc, C, CC, P)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
C
C
C  Computes the density matrix from the MOs
C  (lowest NOcc MOs assumed to be occupied)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions (atomic orbitals)
C  NOcc    -  number of occupied MOs
C  C       -  eigenvectors
C  CC      -  scratch storage for transposed MOs
C
C  on exit
C
C  P       -  density matrix
C
C
      DIMENSION P(*), C(NBas,NOcc), CC(NOcc,NBas)
C
C
C  Transpose MO coefficients
C
      DO 10 I=1,NBas
      DO 10 J=1,NOcc
      CC(J,I) = C(I,J)
 10   CONTINUE
C
C  form Density
C
      IJ = 0
      DO 30 I=1,NBas
      DO 30 J=1,I
      IJ = IJ+1
      VAL = 0.0d0
      DO 20 K=1,NOcc
      VAL = VAL + CC(K,I)*CC(K,J)
 20   CONTINUE
      P(IJ) = VAL
 30   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE DIAT(NI,NJ,XI,XJ,DI)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
C
C
C  Calculates the diatomic overlap integrals between atoms
C  centered at XI and XJ
C
C  ARGUMENTS
C
C  NI      -  atomic number of first atom
C  NJ      -  atomic number of second atom
C  XI      -  Cartesian coordinates of first atom
C  XJ      -  Cartesian coordinates of second atom
C  DI      -  diatomic overlap, in a 9 x 9 matrix, order is
C              1   2   3   4   5            6     7       8     9
C              S   PX  PY  PZ  D(X**2-Y**2) D(XZ) D(Z**2) D(YZ) D(XY
C
C  Limitations
C  -----------
C  NI/NJ must be less than 54
C  Exponents are assumed to be in common EXPONT
C
C
      DIMENSION XI(3),XJ(3)
      INTEGER PQ1,A,PQ2,B,PQA,PQB,AA,BB,PVAL1,PVAL2,YETA
      LOGICAL FIRST
      COMMON /EXPONT/ EMUS(54),EMUP(54),EMUD(54)
      Dimension DI(9,9),S(3,3,3),UL1(3),UL2(3),C(3,5,5),NPQ(54)
C
      DATA NPQ/1,0, 2,2,2,2,2,2,2,0, 3,3,3,3,3,3,3,0, 4,4,4,4,4,4,4,4,
     +4,4,4,4,4,4,4,4,4,0, 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/
      DATA FIRST /.TRUE./
C
C
C  initialize
C
      X1=XI(1)
      X2=XJ(1)
      Y1=XI(2)
      Y2=XJ(2)
      Z1=XI(3)
      Z2=XJ(3)
      PQ1=NPQ(NI)
      PQA=PQ1
      PQ2=NPQ(NJ)
      PQB=PQ2
      UL1(1)=EMUS(NI)
      UL2(1)=EMUS(NJ)
      UL1(2)=EMUP(NI)
      UL2(2)=EMUP(NJ)
      UL1(3)=EMUD(NI)
      UL2(3)=EMUD(NJ)
      CALL ZeroIT(DI,81)
c
      CALL COE(X1,Y1,Z1,X2,Y2,Z2,PQ1,PQ2,C,R)
      If(PQ1.EQ.0.OR.PQ2.EQ.0) RETURN
      If(PQ1.GT.3) GO TO 351
      A=PQ1-1
      GO TO 352
 351  A=2
 352  CONTINUE
      If(PQ2.GT.3) GO TO 353
      B=PQ2-1
      GO TO 354
 353  B=2
 354  CONTINUE
      IA=A+1
      IB=B+1
c
      IF(PQ1.LT.3.AND.PQ2.LT.3) THEN
        CALL DIAT2(NI,EMUS(NI),EMUP(NI),R,NJ,EMUS(NJ),EMUP(NJ),S)
cc      WRITE(6,'(3I4,F12.6)')(((I,J,K,S(I,J,K),I=1,2),J=1,2),K=1,2)
      ELSE
        CALL ZeroIT(S,27)
        DO 363 I=1,IA
        IB=B+1
        DO 363 J=1,IB
        IF (A.LT.B) GO TO 357
        NEWK=B
        GO TO 358
 357    NEWK=A
 358    CONTINUE
        NK1=NEWK+1
        DO 363 K=1,NK1
        IF(K.GT.I.OR.K.GT.J) GOTO 363
        U1=UL1(I)
        U2=UL2(J)
        YETA=PQA+PQB+3
        ISS=I
        JSS=J
        KSS=K
        S(I,J,K)=SS(PQA,PQB,ISS,JSS,KSS,U1,U2,R,YETA,FIRST)
cc      WRITE(6,'(3I4,F12.6)')I,J,K,S(I,J,K)
 363    CONTINUE
      ENDIF
      DO 369 I=1,IA
      DO 369 J=1,IB
      If(I.EQ.1) GO TO 10
      If(I.EQ.2) GO TO 11
      KMIN=1
      KMAX=5
      GO TO 12
 11   KMIN=2
      KMAX=4
      GO TO 12
 10   KMIN=3
      KMAX=3
 12   DO 369 K=KMIN,KMAX
      If(J.EQ.1) GO TO 13
      If(J.EQ.2) GO TO 14
      LMIN=1
      LMAX=5
      GO TO 15
 14   LMIN=2
      LMAX=4
      GO TO 15
 13   LMIN=3
      LMAX=3
 15   DO 369 L=LMIN,LMAX
      IF (J.EQ.2) GO TO 364
      AA=1
      GO TO 365
 364  AA=-1
 365  CONTINUE
      IF (J.GT.2) GO TO 366
      BB=1
      GO TO 367
 366  BB=-1
 367  CONTINUE
      IVAL=I
      KVAL=K
      JVAL=J
      LVAL=L
      CALL VAL(IVAL,KVAL,PVAL1)
      CALL VAL(JVAL,LVAL,PVAL2)
      DI((PVAL1+1),(PVAL2+1))=S(I,J,1)*C(I,K,3)*C(J,L,3)*AA+(C(I,K,4)*C(
     1J,L,4)+C(I,K,2)*C(J,L,2))*BB*S(I,J,2)+(C(I,K,5)*C(J,L,5)+C(I,K,1)*
     2C(J,L,1))*S(I,J,3)
 369  CONTINUE
C
      RETURN
      END
c ........................................................................
c
      SUBROUTINE COE(X1,Y1,Z1,X2,Y2,Z2,PQ1,PQ2,C,R)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      INTEGER PQ1,PQ2,PQ,CO
      DIMENSION C(3,5,5)
      XY=(X2-X1)**2+(Y2-Y1)**2
      R=SQRT(XY+(Z2-Z1)**2)
      XY=SQRT(XY)
      IF (XY.EQ.0.D0) GO TO 800
      CA=(X2-X1)/XY
      CB=(Z2-Z1)/R
      SA=(Y2-Y1)/XY
      SB=XY/R
      GO TO 804
 800  IF (Z2-Z1) 801,802,803
 801  CA=-1.D0
      CB=-1.D0
      SA=0.D0
      SB=0.D0
      GO TO 804
 802  CA=0.D0
      CB=0.D0
      SA=0.D0
      SB=0.D0
      GO TO 804
 803  CA=1.D0
      CB=1.D0
      SA=0.D0
      SB=0.D0
 804  CONTINUE
      CO=0
      CALL ZeroIT(C,75)
      IF (PQ1.GT.PQ2) GO TO 806
      PQ=PQ2
      GO TO 807
 806  PQ=PQ1
 807  CONTINUE
      C(1,3,3)=1.D0
      IF (PQ.LT.2) GO TO 808
      C(2,4,4)=CA*CB
      C(2,4,3)=CA*SB
      C(2,4,2)=-SA
      C(2,3,4)=-SB
      C(2,3,3)=CB
      C(2,3,2)=0.D0
      C(2,2,4)=SA*CB
      C(2,2,3)=SA*SB
      C(2,2,2)=CA
      IF (PQ.LT.3) GO TO 808
      C2A=2*CA*CA-1.D0
      C2B=2*CB*CB-1.D0
      S2A=2*SA*CA
      S2B=2*SB*CB
      C(3,5,5)=C2A*CB*CB+0.5D0*C2A*SB*SB
      C(3,5,4)=0.5D0*C2A*S2B
      C(3,5,3)=0.8660254037841D0*C2A*SB*SB
      C(3,5,2)=-S2A*SB
      C(3,5,1)=-S2A*CB
      C(3,4,5)=-0.5D0*CA*S2B
      C(3,4,4)=CA*C2B
      C(3,4,3)=0.8660254037841D0*CA*S2B
      C(3,4,2)=-SA*CB
      C(3,4,1)=SA*SB
      C(3,3,5)=0.5773502691894D0*SB*SB*1.5D0
      C(3,3,4)=-0.8660254037841D0*S2B
      C(3,3,3)=CB*CB-0.5D0*SB*SB
      C(3,2,5)=-0.5D0*SA*S2B
      C(3,2,4)=SA*C2B
      C(3,2,3)=0.8660254037841D0*SA*S2B
      C(3,2,2)=CA*CB
      C(3,2,1)=-CA*SB
      C(3,1,5)=S2A*CB*CB+0.5D0*S2A*SB*SB
      C(3,1,4)=0.5D0*S2A*S2B
      C(3,1,3)=0.8660254037841D0*S2A*SB*SB
      C(3,1,2)=C2A*SB
      C(3,1,1)=C2A*CB
 808  CONTINUE
C     S(I,J,K)=SS(PQA,PQB,ISS,JSS,KSS,U1,U2,R,YETA,FIRST)
      RETURN
      END
c .....................................................................
c
      SUBROUTINE BFN(X,BF)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms the "B" integrals for the overlap
C
      DIMENSION BF(13)
      Dimension FACT(17)
      DATA FACT/1.D0,2.D0,6.D0,24.D0,120.D0,720.D0,5040.D0,40320.D0,
     $362880.D0,3628800.D0,39916800.D0,479001600.D0,6227020800.D0,
     $8.71782912D10,1.307674368D12,2.092278989D13,3.556874281D14/
C
      K=12
      IO=0
      ABSX = ABS(X)
      If(ABSX.GT.3.D00) GO TO 40
      If(ABSX.LE.2.D00) GO TO 10
      LAST=15
      GO TO 60
   10 If(ABSX.LE.1.D00) GO TO 20
      LAST=12
      GO TO 60
   20 If(ABSX.LE.0.5D00) GO TO 30
      LAST=7
      GO TO 60
   30 If(ABSX.LE.1.D-6) GO TO 90
      LAST=6
      GO TO 60
   40 EXPX=EXP(X)
      EXPMX=1.D00/EXPX
      BF(1)=(EXPX-EXPMX)/X
      DO 50 I=1,K
   50 BF(I+1)=(DBLE(I)*BF(I)+(-1.D00)**I*EXPX-EXPMX)/X
      GO TO 110
   60 DO 80 I=IO,K
      Y=0.0D00
      DO 70 M=IO,LAST
      XF=1.0D00
      If(M.NE.0) XF=FACT(M)
   70 Y=Y+(-X)**M*DBLE(2*MOD(M+I+1,2))/(XF*DBLE(M+I+1))
   80 BF(I+1)=Y
      GO TO 110
   90 DO 100 I=IO,K
  100 BF(I+1)=DBLE(2*MOD(I+1,2))/DBLE(I+1)
  110 CONTINUE
C
      RETURN
      END
c .....................................................................
c
      DOUBLE PRECISION FUNCTION SS(NA,NB,LA,LB,M,UC,UD,R1,YETA,FIRST)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      LOGICAL FIRST
      INTEGER A,C,D,PP,B,Q,YETA
      DIMENSION FA(14),BI(13,13),AFF(3,3,3),AF(20),BF(20)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,cc,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,RGas
c
      DATA AFF/27*0.0d0/
      R=R1
      UA=UC
      UB=UD
      If(UA.GT.0.0d0) GO TO 88
      SA=0.D0
      GO TO 99
  88  If(UB.GT.0.0d0) GO TO 299
      SA=0.0d0
      GO TO 99
 299  CONTINUE
      R=R*ANTOAU
      ER=R
      GO TO 304
 300  FA(1)=1.0d0
      DO 301 I=1,13
      FA(I+1)=FA(I)*I
 301  CONTINUE
      FIRST=.FALSE.
      DO 302 I=1,13
      BI(I,1)=1.0d0
      BI(I,I)=1.0d0
 302  CONTINUE
      DO 303 I=1,12
      I1=I-1
      DO 303 J=1,I1
      BI(I+1,J+1)=BI(I,J+1)+BI(I,J)
 303  CONTINUE
      AFF(1,1,1)=1.0d0
      AFF(2,1,1)=1.0d0
      AFF(2,2,1)=1.0d0
      AFF(3,1,1)=1.5d0
      AFF(3,2,1)=1.73205d0
      AFF(3,3,1)=1.224745d0
      AFF(3,1,3)=-0.5d0
      GO TO 305
 304  If(FIRST) GO TO 300
 305  CONTINUE
      P=(UA+UB)*ER*0.5d0
      BA=(UA-UB)*ER*0.5d0
      EX=EXP(BA)
      QUO=1.0d0/P
      AF(1)=QUO*EXP(-P)
      NANB=NA+NB
      DO 306 N=1,19
      AF(N+1)=N*QUO*AF(N)+AF(1)
 306  CONTINUE
      NANB1=NANB+1
      CALL BFN(BA,BF)
      SUM=0.0d0
      LAM1=LA-M+1
      LBM1=LB-M+1
      DO 311 I=1,LAM1,2
      DO 311 J=1,LBM1,2
      A=NA+I-LA
      B=NB+J-LB
      C=LA-I-M+1
      D=LB-J-M+1
      SUM1=0.0d0
      IA=A+1
      IB=B+1
      IC=C+1
      ID=D+1
      AB=A+B-1
      DO 310 K1=1,IA
      DO 310 K2=1,IB
      DO 310 K3=1,IC
      DO 310 K4=1,ID
      DO 310 K5=1,M
      DO 310 K6=1,M
      Q=AB-K1-K2+K3+K4+2*K5
      PP=K1+K2+K3+K4+2*K6-5
      JX=M+K2+K4+K5+K6-5
      IX=JX/2
      SUM1=SUM1+BI(IA,K1)*BI(IB,K2)*BI(IC,K3)*BI(ID,K4)*BI(M,K5)*BI(
     1M,K6)*2.D0*(IX*2-JX+0.5D0)*AF(Q)*BF(PP)
 310  CONTINUE
      SUM=SUM+SUM1*AFF(LA,M,I)*AFF(LB,M,J)
 311  CONTINUE
      X=R
      DO 312 I=1,NA
      X=X*R*UA
 312  CONTINUE
      DO 313 I=1,NB
      X=X*R*UB
 313  CONTINUE
      SA=SUM*X*SQRT(UA*UB/(FA(NA+NA+1)*FA(NB+NB+1))*((LA+LA-1
     1)*(LB+LB-1)))/(2.D0**M)
 99   CONTINUE
      SS=SA
C
      RETURN
      END
c ........................................................................
c
      SUBROUTINE VAL(L,M,PVAL)
      INTEGER L,M,PVAL
C
      If(L.EQ.1) Then
        PVAL=0
      Else If(L.EQ.2) Then
        If(M.EQ.2) Then
          PVAL=2
        Else If(M.EQ.3) Then
          PVAL=3
        Else If(M.EQ.4) Then
          PVAL=1
        EndIf
      Else If(L.EQ.3) Then
        PVAL=12-L-M
      EndIf
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE DIAT2(NA,ESA,EPA,R12,NB,ESB,EPB,S)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Computes overlap between atomic orbitals for pairs of atoms
C  Can handle 1S, 2S, 3S, 2P and 3P orbitals
C
      DIMENSION S(3,3,3)
      COMMON /SETC/ A(7),B(7),SA,SB,FACTOR,ISP,IPS
      Dimension IROW(17),IOV(3,3)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,RGas
c
      DATA IROW /1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3/
      DATA IOV(1,1),IOV(2,1),IOV(3,1) /1, 2, 3/
      DATA IOV(1,2),IOV(2,2),IOV(3,2) /2, 4, 5/
      DATA IOV(1,3),IOV(2,3),IOV(3,3) /3, 5, 6/
C
C  numbering corresponds to bond type
C     IOV=1  First row - First row elements
C        =2  First row - Second row
C        =3  First row - Third row
C        =4  Second row - Second row
C        =5  Second row - Third row
C        =6  Third row - Third row
C
C  Assign bond type
C
      IR=IROW(NA)
      JR=IROW(NB)
      II=IOV(IR,JR)
      CALL ZeroIT(S,27)
      RAB=R12*ANTOAU
      GO TO (20,30,40,50,60,70), II
C
C ------------------------------------------------------------------
C  The ordering of the elements within S is:
C  S(1,1,1)=(S(B)/S(A))
C  S(1,2,1)=(P-SIGMA(B)/S(A))
C  S(2,1,1)=(S(B)/P-SIGMA(A))
C  S(2,2,1)=(P-SIGMA(B)/P-SIGMA(A))
C  S(2,2,2)=(P-PI(B)/P-PI(A))
C ------------------------------------------------------------------
C
C  First row - First row overlap
C
   20 CALL SET(ESA,ESB,NA,NB,RAB,II)
      S(1,1,1)=.25D00*SQRT((SA*SB*RAB*RAB)**3)*(A(3)*B(1)-B(3)*A(1))
      RETURN
C
C  First row - Second row overlap
C
   30 CALL SET(ESA,ESB,NA,NB,RAB,II)
      W=SQRT((SA**3)*(SB**5))*(RAB**4)*0.125D00
      S(1,1,1) = SQRT(1.D00/3.D00)
      S(1,1,1)=W*S(1,1,1)*(A(4)*B(1)-B(4)*A(1)+A(3)*B(2)-B(3)*A(2))
      If(NA.GT.1) CALL SET(EPA,ESB,NA,NB,RAB,II)
      If(NB.GT.1) CALL SET(ESA,EPB,NA,NB,RAB,II)
      W=SQRT((SA**3)*(SB**5))*(RAB**4)*0.125D00
      S(ISP,IPS,1)=W*(A(3)*B(1)-B(3)*A(1)+A(4)*B(2)-B(4)*A(2))
      RETURN
C
C  First row - Third row overlap
C
   40 CALL SET(ESA,ESB,NA,NB,RAB,II)
      W=SQRT((SA**3)*(SB**7)/7.5D00)*(RAB**5)*0.0625D00
      SROOT3 = SQRT(3.D00)
      S(1,1,1)=W*(A(5)*B(1)-B(5)*A(1)+
     $          2.D00*(A(4)*B(2)-B(4)*A(2)))/SROOT3
      If(NA.GT.1) CALL SET(EPA,ESB,NA,NB,RAB,II)
      If(NB.GT.1) CALL SET(ESA,EPB,NA,NB,RAB,II)
      W=SQRT((SA**3)*(SB**7)/7.5D00)*(RAB**5)*0.0625D00
      S(ISP,IPS,1)=W*(A(4)*(B(1)+B(3))-B(4)*(A(1)+A(3))+
     $              B(2)*(A(3)+A(5))-A(2)*(B(3)+B(5)))
      RETURN
C
C  Second row - Second row overlap
C
   50 CALL SET(ESA,ESB,NA,NB,RAB,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      RT3=1.D00/SQRT(3.D00)
      S(1,1,1)=W*(A(5)*B(1)+B(5)*A(1)-2.0D00*A(3)*B(3))/3.0D00
      CALL SET(ESA,EPB,NA,NB,RAB,II)
      If(NA.GT.NB) CALL SET(EPA,ESB,NA,NB,RAB,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      D=A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
      E=B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
      S(ISP,IPS,1)=W*RT3*(D+E)
      CALL SET(EPA,ESB,NA,NB,RAB,II)
      If(NA.GT.NB) CALL SET(ESA,EPB,NA,NB,RAB,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      D=A(4)*(B(1)-B(3))-A(2)*(B(3)-B(5))
      E=B(4)*(A(1)-A(3))-B(2)*(A(3)-A(5))
      S(IPS,ISP,1)=-W*RT3*(E-D)
      CALL SET(EPA,EPB,NA,NB,RAB,II)
      W=SQRT((SA*SB)**5)*(RAB**5)*0.0625D00
      S(2,2,1)=-W*(B(3)*(A(5)+A(1))-A(3)*(B(5)+B(1)))
      HD = .5D00
      S(2,2,2)=HD*W*(A(5)*(B(1)-B(3))-B(5)*(A(1)-A(3))
     $           -A(3)*B(1)+B(3)*A(1))
      RETURN
C
C  Second row - Third row overlap
C
   60 CALL SET(ESA,ESB,NA,NB,RAB,II)
      W=SQRT((SA**5)*(SB**7)/7.5D00)*(RAB**6)*0.03125D00
      RT3 = 1.D00/SQRT(3.D00)
      TD = 2.D00
      S(1,1,1)=W*(A(6)*B(1)+A(5)*B(2)-TD*(A(4)*B(3)+
     $          A(3)*B(4))+A(2)*B(5)+A(1)*B(6))/3.0d0
      CALL SET(ESA,EPB,NA,NB,RAB,II)
      If(NA.GT.NB) CALL SET(EPA,ESB,NA,NB,RAB,II)
      W=SQRT((SA**5)*(SB**7)/7.5D00)*(RAB**6)*0.03125D00
      TD = 2.D00
      S(ISP,IPS,1)=W*RT3*(A(6)*B(2)+A(5)*B(1)-TD*(A(4)*B(4)+A(3)*B(3))
     $              +A(2)*B(6)+A(1)*B(5))
      CALL SET(EPA,ESB,NA,NB,RAB,II)
      If(NA.GT.NB) CALL SET(ESA,EPB,NA,NB,RAB,II)
      W=SQRT((SA**5)*SB**7/7.5D00)*(RAB**6)*0.03125D00
      TD = 2.D00
      S(IPS,ISP,1)=-W*RT3*(A(5)*(TD*B(3)-B(1))-B(5)*(TD*A(3)-A(1))
     $               -A(2)*(B(6)-TD*B(4))+B(2)*(A(6)-TD*A(4)))
      CALL SET(EPA,EPB,NA,NB,RAB,II)
      W=SQRT((SA**5)*SB**7/7.5D00)*(RAB**6)*0.03125D00
      S(2,2,1)=-W*(B(4)*(A(1)+A(5))-A(4)*(B(1)+B(5))
     $            +B(3)*(A(2)+A(6))-A(3)*(B(2)+B(6)))
      HD = .5D00
      S(2,2,2)=HD*W*(A(6)*(B(1)-B(3))-B(6)*(A(1)-
     $           A(3))+A(5)*(B(2)-B(4))-B(5)
     $           *(A(2)-A(4))-A(4)*B(1)+B(4)*A(1)-A(3)*B(2)+B(3)*A(2))
      RETURN
C
C  Third row - Third row overlap
C
   70 CALL SET(ESA,ESB,NA,NB,RAB,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      RT3 = 1.D00/SQRT(3.D00)
      S(1,1,1)=W*(A(7)*B(1)-3.D00*(A(5)*B(3)-A(3)*B(5))-A(1)*B(7))/3.D00
      CALL SET(ESA,EPB,NA,NB,RAB,II)
      If(NA.GT.NB) CALL SET(EPA,ESB,NA,NB,RAB,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      D=A(6)*(B(1)-B(3))-2.D00*A(4)*(B(3)-B(5))+A(2)*(B(5)-B(7))
      E=B(6)*(A(1)-A(3))-2.D00*B(4)*(A(3)-A(5))+B(2)*(A(5)-A(7))
      S(ISP,IPS,1)=W*RT3*(D-E)
      CALL SET(EPA,ESB,NA,NB,RAB,II)
      If(NA.GT.NB) CALL SET(ESA,EPB,NA,NB,RAB,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      D=A(6)*(B(1)-B(3))-2.D00*A(4)*(B(3)-B(5))+A(2)*(B(5)-B(7))
      E=B(6)*(A(1)-A(3))-2.D00*B(4)*(A(3)-A(5))+B(2)*(A(5)-A(7))
      S(IPS,ISP,1)=-W*RT3*(-D-E)
      CALL SET(EPA,EPB,NA,NB,RAB,II)
      W=SQRT((SA*SB*RAB*RAB)**7)/480.D00
      TD = 2.D00
      S(2,2,1)=-W*(A(3)*(B(7)+TD*B(3))-A(5)*(B(1)+
     $           TD*B(5))-B(5)*A(1)+A(7)*B(3))
      HD = .5D00
      S(2,2,2)=HD*W*(A(7)*(B(1)-B(3))+B(7)*(A(1)-
     $          A(3))+A(5)*(B(5)-B(3)-B(1))
     $          +B(5)*(A(5)-A(3)-A(1))+2.D00*A(3)*B(3))
C
      RETURN
      END
c ......................................................................
c
      SUBROUTINE SET(S1,S2,NA,NB,RAB,II)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /SETC/ A(7),B(7),SA,SB,FACTOR,ISP,IPS
C
      IF(NA.GT.NB) THEN
        ISP = 2
        IPS = 1
        SA = S2
        SB = S1
      ELSE
        ISP = 1
        IPS = 2
        SA = S1
        SB = S2
      ENDIF
      J=II+2
      If(II.GT.3) J=J-1
      ALPHA=0.5d0*RAB*(SA+SB)
      BETA=0.5d0*RAB*(SB-SA)
      JCALL=J-1
      CALL AINTGS(ALPHA,JCALL)
      CALL BINTGS(BETA,JCALL)
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE AINTGS(X,K)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /SETC/ A(7),B(7),SDUM(3),IDUM(2)
C
C  Forms the "A" integrals for the overlap calculation
C
      C=EXP(-X)
      A(1)=C/X
      DO 10 I=1,K
      A(I+1)=(A(I)*DBLE(I)+C)/X
   10 CONTINUE
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE BINTGS(X,K)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /SETC/ A(7),B(7),SDUM(3),IDUM(2)
      DIMENSION FACT(17)
C
C  Forms the "B" integrals for the overlap calculation
C
      DATA FACT/1.D0,2.D0,6.D0,24.D0,120.D0,720.D0,5040.D0,40320.D0,
     $362880.D0,3628800.D0,39916800.D0,479001600.D0,6227020800.D0,
     $8.71782912D10,1.307674368D12,2.092278989D13,3.556874281D14/
C
      IO = 0
      ABSX = ABS(X)
      IF(AbsX.gt.300.0d0) go to 105
      If(ABSX.GT.3.D00) GO TO 40
      If(ABSX.LE.2.D00) GO TO 10
      If(K.LE.10) GO TO 40
      LAST=15
      GO TO 60
   10 If(ABSX.LE.1.D00) GO TO 20
      If(K.LE.7) GO TO 40
      LAST=12
      GO TO 60
   20 If(ABSX.LE.0.5D00) GO TO 30
      If(K.LE.5) GO TO 40
      LAST=7
      GO TO 60
   30 If(ABSX.LE.1.D-6) GOTO 90
      LAST=6
      GO TO 60
   40 EXPX=EXP(X)
      EXPMX=1.D00/EXPX
      B(1)=(EXPX-EXPMX)/X
      DO 50 I=1,K
   50 B(I+1)=(DBLE(I)*B(I)+(-1.D00)**I*EXPX-EXPMX)/X
      GO TO 110
   60 DO 80 I=IO,K
      Y=0.0D00
      DO 70 M=IO,LAST
      XF=1.0D00
      If(M.NE.0) XF=FACT(M)
   70 Y=Y+(-X)**M*DBLE(2*MOD(M+I+1,2))/(XF*DBLE(M+I+1))
   80 B(I+1)=Y
      GO TO 110
   90 DO 100 I=IO,K
  100 B(I+1)=DBLE(2*MOD(I+1,2))/DBLE(I+1)
      GO TO 110
  105 DO I=IO,K
       B(I+1)=0.0d0
      end do
  110 CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE FOCK1(NAtoms, IAN,    NFIRST, NMIDLE, NLAST,
     $                 PTOT,   PA,     PB,     F)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
C
C  Compute the remaining 1-centre contributions to the Fock matrix
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  NFIRST
C  NMIDLE
C  NLAST
C  PTOT    -  total density matrix
C  PA      -  alpha density matrix
C  PB      -  beta density matrix
C
C  on exit
C
C  F       -  partial Fock matrix
C
C
      DIMENSION F(*), PTOT(*), PA(*), PB(*)
      DIMENSION IAN(NAtoms),NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms)
      COMMON /GAUSS / FN1(54),FN2(54)
      COMMON /TWOELE/ GSS(54),GSP(54),GPP(54),GP2(54),HSP(54),
     $                GSD(54),GPD(54),GDD(54)
c ....................................................
      DIMENSION QTOT(NAtoms)
c ....................................................
C
      PARAMETER (Zero=0.0d0)
C
      CALL CHRGE(NAtoms, NFIRST, NLAST,  PTOT,   QTOT)
c
      DO 100 II=1,NAtoms
       IA=NFIRST(II)
       IB=NMIDLE(II)
       IC=NLAST(II)
       NI=IAN(II)
       DTPOP=Zero
       DAPOP=Zero
       PTPOP=Zero
       PAPOP=Zero
       GOTO (11,21,21,21,31,31,31,31,31)IC-IA+1
 31     DTPOP=PTOT((IC*(IC+1))/2)+PTOT(((IC-1)*(IC))/2)
     +       +PTOT(((IC-2)*(IC-1))/2)+PTOT(((IC-3)*(IC-2))/2)
     +       +PTOT(((IC-4)*(IC-3))/2)
        DAPOP=PA((IC*(IC+1))/2)+PA(((IC-1)*(IC))/2)
     +       +PA(((IC-2)*(IC-1))/2)+PA(((IC-3)*(IC-2))/2)
     +       +PA(((IC-4)*(IC-3))/2)
 21     PTPOP=PTOT((IB*(IB+1))/2)+PTOT(((IB-1)*(IB))/2)
     +       +PTOT(((IB-2)*(IB-1))/2)
        PAPOP=PA((IB*(IB+1))/2)+PA(((IB-1)*(IB))/2)
     +       +PA(((IB-2)*(IB-1))/2)
 11     DBPOP=DTPOP-DAPOP
        PBPOP=PTPOP-PAPOP
      IF(NI.EQ.1)THEN
        SUM=Zero
      ELSE
        SUM2=Zero
        SUM1=Zero
        DO 111 I=IA,IB
        DO 112 J=IA,I-1
 112    SUM1=SUM1+PTOT(J+(I*(I-1))/2)**2
 111    SUM2=SUM2+PTOT((I*(I+1))/2)**2
        SUM=SUM1+SUM1+SUM2
        SUM=SQRT(SUM)-QTOT(II)*0.5D0
      ENDIF
cc      WRITE(6,'('' ATOM'',I3,'' CORRECTION'',F12.5)') II,SUM
      SUM=SUM*FN1(NI)
C
C  F(S,S)
C
         KA=(IA*(IA+1))/2
         F(KA)=F(KA)+PB(KA)*GSS(NI)+PTPOP*GSP(NI)
     +         -PAPOP*HSP(NI) + DTPOP*GSD(NI)
     #         + SUM  ! modification to account for dipolar bonds!
c
         If(NI.LT.3) GO TO 100
         IPLUS=IA+1
         L=KA
         DO 80 J=IPLUS,IB
            M=L+IA
            L=L+J
C
C  F(P,P)
C
            F(L)=F(L)+PTOT(KA)*GSP(NI)-PA(KA)*HSP(NI)+
     1      PB(L)*GPP(NI)+(PTPOP-PTOT(L))*GP2(NI)
     2      -0.5D0*(PAPOP-PA(L))*(GPP(NI)-GP2(NI))
     3      +DTPOP*GPD(NI)
     #         + SUM  ! modification to account for dipolar bonds!
C
C  F(S,P)
C
   80    F(M)=F(M)+2.D0*PTOT(M)*HSP(NI)-PA(M)*(HSP(NI)+GSP(NI))
C
C  F(P,P*)
C
      IMINUS=IB-1
      DO 90 J=IPLUS,IMINUS
       ICC=J+1
      DO 90 L=ICC,IB
       M=(L*(L-1))/2+J
  90   F(M)=F(M)+PTOT(M)*(GPP(NI)-GP2(NI))
     +      -0.5D0*PA(M)*(GPP(NI)+GP2(NI))
      DO 95 J=IB+1,IC
      M=(J*(J+1))/2
  95  F(M)=F(M)+PTOT(KA)*GSD(NI)
     +         +PTPOP*GPD(NI)
     +         +DTPOP*GDD(NI)
  100 CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE FOCK2(NAtoms, SEMI,   NFIRST, NMIDLE, NLAST,
     $                 IFACT,  I1FACT, W,      PTOT,   P,  F)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms the 2-electron 2-centre repulsion part of the Fock Matrix
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C  NFIRST
C  NMIDLE
C  NLAST
C  IFACT   -  triangular numbers   I*(I-1)/2
C  I1FACT  -    ditto              I*(I+1)/2
C  W       -  two-electron matrix
C  PTOT    -  total density matrix
C  P       -  alpha- or beta-density matrix
C
C  on exit
C
C  F       -  partial Fock matrix
C
C
      DIMENSION F(*), PTOT(*), P(*), W(*), NFIRST(*), NMIDLE(*),
     +          NLAST(*)
      DIMENSION IFACT(*),I1FACT(*)
c ......................................................
      DIMENSION SPPOP(NAtoms), DPOP(NAtoms)
c ......................................................
C
      Character*20 SEMI
C
C
      IF(SEMI(1:5).NE.'mindo') THEN
cc
       KK=0
       DO 61 II=1,NAtoms
       IA=NFIRST(II)
       IB=NLAST(II)
       IC=NMIDLE(II)
       SUM=0.D0
       DO 74 I=IA,IC
  74   SUM=SUM+PTOT(I1FACT(I))
       SPPOP(II)=SUM
       SUM=0.D0
       DO 75 I=IC+1,IB
  75   SUM=SUM+PTOT(I1FACT(I))
       DPOP(II)=SUM
       IMINUS=II-1
       DO 69 JJ=1,IMINUS
        JA=NFIRST(JJ)
        JB=NLAST(JJ)
        JC=NMIDLE(JJ)
        DREP=W(KK+1)
        DO 60 I=IA,IC
         KA=IFACT(I)
        DO 60 J=IA,I
         KB=IFACT(J)
         IJ=KA+J
         AA=2.0D00
         IF (I.EQ.J) AA=1.0D00
        DO 60 K=JA,JC
         KC=IFACT(K)
         IK=KA+K
         JK=KB+K
        DO 60 L=JA,K
         IL=KA+L
         JL=KB+L
         KL=KC+L
         BB=2.0D00
         If (K.EQ.L) BB=1.0D00
         KK=KK+1
         A=W(KK)
C
C  A is the repulsion integral (IJ|KL) where orbitals I and J are
C  on atom II and K and L are on atom JJ.
C  AA and BB are correction factors (apparently taking into
C  account equivalence under interchange of indices)
C  IJ is the location of the matrix element between orbitals I & J
C
C  THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
C  MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY,
C  BUT IT WORKS
C
         F(IJ)=F(IJ)+BB*A*PTOT(KL)
         F(KL)=F(KL)+AA*A*PTOT(IJ)
         A=A*AA*BB*0.25D0
         F(IK)=F(IK)-A*P(JL)
         F(IL)=F(IL)-A*P(JK)
         F(JK)=F(JK)-A*P(IL)
         F(JL)=F(JL)-A*P(IK)
   60    CONTINUE
C
C D-orbital correction
C
       DO 62 I=IC+1,IB
       KA=IFACT(I)
       DO 62 J=JA,JB
       IJ=KA+J
C
C  Atom J (S, P, and D (if present)) exchange with atom I (D only)
C
  62   F(IJ)=F(IJ)-0.5D0*DREP*P(IJ)
       DO 63 I=IA,IC
       KA=IFACT(I)
       DO 63 J=JC+1,JB
       IJ=KA+J
C
C  Atom J (D (if present)) exchange with atom I (S and P only)
C
  63   F(IJ)=F(IJ)-0.5D0*DREP*P(IJ)
C
C  THE COULOMB REPULSION TERMS
C  ---------------------------
C
C  Atom J (S, P and D shells) with atom I (D shell)
C
       DO 68 J=JA,JB
       J2=I1FACT(J)
  68   F(J2)=F(J2)+DREP*DPOP(II)
C
C  Atom J (D shell) with Atom I (S and P shells)
C
       DO 67 J=JC+1,JB
       J2=I1FACT(J)
  67   F(J2)=F(J2)+DREP*SPPOP(II)
C
C  Atom I (S, P and D shells) with atom J (D shell)
C
       DO 66 I=IA,IB
       I2=I1FACT(I)
  66   F(I2)=F(I2)+DREP*DPOP(JJ)
C
C  Atom I (D shell) with atom J (S and P shells)
C
       DO 65 I=IC+1,IB
       I2=I1FACT(I)
  65   F(I2)=F(I2)+DREP*SPPOP(JJ)
c
  69   CONTINUE
  61   CONTINUE
cc
      ELSE
cc
       KR=0
       DO 210 II=1,NAtoms
       IA=NFIRST(II)
       IB=NLAST(II)
       IM1=II-1
       DO 220 JJ=1,IM1
       KR=KR+1
       ELREP=W(KR)
       JA=NFIRST(JJ)
       JB=NLAST(JJ)
       DO 240 I=IA,IB
       KA=IFACT(I)
       KK=KA+I
       DO 240 K=JA,JB
       LL=I1FACT(K)
       IK=KA+K
       F(KK)=F(KK)+PTOT(LL)*ELREP
       F(LL)=F(LL)+PTOT(KK)*ELREP
 240   F(IK)=F(IK)-P(IK)*ELREP
 220   CONTINUE
 210   CONTINUE
      ENDIF
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE H1ELEC(NI,NJ,XI,XJ,SEMI,SMAT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine forms the one-electron matrix between two atoms
C
C  ARGUMENTS
C
C  NI      -  atomic number of first atom
C  NJ      -  atomic number of second atom
C  XI      -  Cartesian coordinates of first atom
C  XJ      -  Cartesian coordinates of second atom
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C
C  on exit
C
C  SMAT    -  matrix of one-electron interactions
C
      DIMENSION XI(3),XJ(3),SMAT(9,9), BI(9), BJ(9)
      Character*20 SEMI
C
      COMMON /BETAS / BETAS(54),BETAP(54),BETAD(54)
      COMMON /BETA3 / BETA3(153)
      COMMON /VSIPS / VS(54),VP(54),VD(54)
      COMMON /NATORB/ NATORB(54)
C
      PARAMETER (Half=0.5d0)
C
C
      CALL DIAT(NI,NJ,XI,XJ,SMAT)
c
      IF(SEMI(1:5).EQ.'mindo') THEN
       II=MAX(NI,NJ)
       NBOND=(II*(II-1))/2+NI+NJ-II
       BI(1)=BETA3(NBOND)*VS(NI)
       BI(2)=BETA3(NBOND)*VP(NI)
       BI(3)=BI(2)
       BI(4)=BI(2)
       BJ(1)=BETA3(NBOND)*VS(NJ)
       BJ(2)=BETA3(NBOND)*VP(NJ)
       BJ(3)=BJ(2)
       BJ(4)=BJ(2)
      ELSE
       BI(1)=BETAS(NI)*Half
       BI(2)=BETAP(NI)*Half
       BI(3)=BI(2)
       BI(4)=BI(2)
       BI(5)=BETAD(NI)*Half
       BI(6)=BI(5)
       BI(7)=BI(5)
       BI(8)=BI(5)
       BI(9)=BI(5)
       BJ(1)=BETAS(NJ)*Half
       BJ(2)=BETAP(NJ)*Half
       BJ(3)=BJ(2)
       BJ(4)=BJ(2)
       BJ(5)=BETAD(NJ)*Half
       BJ(6)=BJ(5)
       BJ(7)=BJ(5)
       BJ(8)=BJ(5)
       BJ(9)=BJ(5)
      ENDIF
c
      DO 10 J=1,NATORB(NJ)
      DO 10 I=1,NATORB(NI)
      SMAT(I,J)=SMAT(I,J)*(BI(I)+BJ(J))
 10   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE HCORE(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $                 NFIRST, NMIDLE, NLAST,  USPD,   H,
     $                 W,      ENuclr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates the one-electron matrix and the two-electron integrals
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  list of atomic numbers
C  XC      -  Cartesian coordinates
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C  IPRNT   -  print flag
C  NFIRST
C  NMIDLE
C  NLAST
C  USPD    -  one-electron diagonal terms
C
C  on exit
C
C  H       -  one-electron matrix
C  W       -  two-electron integrals
C  ENuclr  -  nuclear repulsion energy
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms)
      DIMENSION NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms)
      DIMENSION USPD(*),H(*),W(*)
      Dimension E1B(10),E2A(10),DI(9,9)
      Character*20 SEMI
C
      PARAMETER (Zero=0.0d0)
C
C
      ENuclr = Zero
      KR = 1
C
C  Loop over all atoms
C
      DO 60 I=1,NAtoms
      IA = NFIRST(I)
      IB = NLAST(I)
      IC = NMIDLE(I)
      NI = IAN(I)
C
C  fill diagonals and off-diagonals of same atom
C
      DO 10 I1=IA,IB
      I2=I1*(I1-1)/2+IA-1
      DO 9 J1=IA,I1
      I2=I2+1
      H(I2)=Zero
 9    CONTINUE
      H(I2) = USPD(I1)
 10   CONTINUE
c
      DO 50 J=1,I-1
      JA = NFIRST(J)
      JB = NLAST(J)
      JC = NMIDLE(J)
      NJ = IAN(J)
      CALL H1ELEC(NI,NJ,XC(1,I),XC(1,J),SEMI,DI)
C
C  fill atom-other atom one-electron matrix
C
      I2=0
      DO 20 I1=IA,IB
      II=I1*(I1-1)/2+JA-1
      I2=I2+1
      J2=0
      DO 19 J1=JA,JB
      II=II+1
      J2=J2+1
      H(II)=DI(I2,J2)
 19   CONTINUE
 20   CONTINUE
C
C  Calculate
C    the two-electron integrals - W
C    the electron-nuclear terms - E1B and E2A
C    the nuclear-nuclear term   - ENuc
C
      CALL ROTAT(NI,     NJ,     XC(1,I),XC(1,J),SEMI,
     $           W(KR),  KR,     E1B,    E2A,    ENuc)
c
      ENuclr = ENuclr + ENuc
C
C  add the electron-nuclear attraction term for atom I
C
      I2=0
      DO 30 I1=IA,IC
      II=I1*(I1-1)/2+IA-1
      DO 29 J1=IA,I1
      II=II+1
      I2=I2+1
      H(II)=H(II)+E1B(I2)
 29   CONTINUE
 30   CONTINUE
      DO 31 I1=IC+1,IB
      II=(I1*(I1+1))/2
      H(II)=H(II)+E1B(1)
 31   CONTINUE
C
C  add the electron-nuclear attraction term for atom J
C
      I2=0
      DO 40 I1=JA,JC
      II=I1*(I1-1)/2+JA-1
      DO 39 J1=JA,I1
      II=II+1
      I2=I2+1
      H(II)=H(II)+E2A(I2)
 39   CONTINUE
 40   CONTINUE
      DO 41 I1=JC+1,JB
      II=(I1*(I1+1))/2
      H(II)=H(II)+E2A(1)
 41   CONTINUE
c
 50   CONTINUE
 60   CONTINUE
c
      If(IPRNT.GT.4) Then
       WRITE(6,'(//10X,''WHOLE OF TWO-ELECTRON MATRIX IN HCORE''/)')
       WRITE(6,111)(W(I),I=1,KR-1)
 111   FORMAT(8F9.4)
      EndIf
C
      RETURN
      END
c =======================================================================
c
      DOUBLE PRECISION FUNCTION HELECT(NBas,P,H,F)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculates the electronic energy of the system in EV
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  P       -  density matrix, packed, lower triangle
C  H       -  one-electron matrix,  ditto
C  F       -  two-electron matrix,  ditto
C
      DIMENSION P(*), H(*), F(*)
C
      ED=0.0d0
      EE=0.0d0
      K=0
      NN=NBas+1
      DO 20 I=2,NN
      K=K+1
      JJ=I-1
      ED=ED+P(K)*(H(K)+F(K))
      If(I.EQ.NN) GO TO 20
      DO 10 J=1,JJ
      K=K+1
   10 EE=EE+P(K)*(H(K)+F(K))
   20 CONTINUE
      EE=EE + 0.5d0*ED
      HELECT=EE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE ITER(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $                NFIRST, NMIDLE, NLAST,  IFACT,  I1FACT,
     $                NBas,   PSPD,   NAlpha, NBeta,  MDiis,
     $                H,      W,      P,      F,      C,
     $                EIGS,   PA,     POLD,   BA,     FB,
     $                CBETA,  EIGB,   PB,     PBOLD,  BB,
     $                ATHEAT, ENuclr, inp,    Z,      EE,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main driving routine for SCF calculation
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C  IPRNT   -  print flag
C  NFIRST
C  NMIDLE
C  NLAST
C  IFACT   -  triangular numbers   I*(I-1)/2
C  I1FACT  -   ditto               I*(I+1)/2
C  PSPD    -  two-electron diagonal terms (?)
C  NBas    -  number of basis functions
C  NAlpha  -  number of occupied alpha or closed-shell orbitals
C  NBeta   -  number of occupied beta orbitals
C  MDiis   -  maximum size of DIIS subspace
C  H       -  one-electron matrix
C  W       -  two-electron matrix
C  P       -  total density matrix
C  F       -  alpha or closed-shell Fock matrix
C  C       -    ditto molecular orbitals
C  EIGS    -    ditto orbital energies
C  PA      -    ditto density matrix
C  POLD    -  old copies of alpha density matrix
C  BA      -  alpha or closed-shell DIIS error matrix
C  FB      -  beta Fock matrix
C  CBETA   -    ditto molecular orbitals
C  EIGB    -    ditto orbital energies
C  PB      -    ditto density matrix
C  PBOLD   -  old copies of beta density matrix
C  BB      -  beta DIIS error matrix (?)
C  Z       -  scratch storage [at least NBas*(NBas+1) ]
C  ATHEAT  -  heat of atomization
C  ENuclr  -  nuclear repulsion energy
C  inp     -  unit number for input file
C             (SCF options are read directly in this routine)
C             (if -1 then diabled - SCF guess)
C
C  on exit
C
C  EE      -  electronic energy
C  IErr    -  error flag for successful termination
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),PSPD(NBas)
      DIMENSION NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms),
     $          IFACT(NBas),I1FACT(NBas)
      DIMENSION H(*),W(*),P(NBas*(NBas+1)/2)
      DIMENSION F(*),C(NBas*NBas),EIGS(NBas),PA(*),
     $          FB(*),CBETA(NBas*NBas),EIGB(NBas),PB(*),
     $          POLD(*),PBOLD(*),
     $          BA(MDiis+1,MDiis+1),BB(MDiis+1,MDiis+1)
      DIMENSION Z(*)
      Character*20 SEMI
      Character*256 jobname,MOS,MOB
      LOGICAL guess,IOMIX,READY,OKNEWD,FRST,FRSTB,UHF,NEWDG,OKPULY
C
C ----------------------------------------------------------------------
C THE MAIN ARRAYS USED IN ITER ARE:
C            P      ONLY EVER CONTAINS THE TOTAL DENSITY MATRIX
C            PA     ONLY EVER CONTAINS THE ALPHA DENSITY MATRIX
C            PB     ONLY EVER CONTAINS THE BETA DENSITY MATRIX
C            C      ONLY EVER CONTAINS THE EIGENVECTORS (MOs)
C            H      ONLY EVER CONTAINS THE ONE-ELECTRON MATRIX
C            F      STARTS OFF CONTAINING THE ONE-ELECTRON MATRIX,
C                   AND IS USED TO HOLD THE FOCK MATRIX
C            W      ONLY EVER CONTAINS THE TWO-ELECTRON MATRIX
C
C  PRINCIPAL REFERENCES:
C
C  Level Shift
C  -----------
C  "Unconditional Convergence in SCF Theory: A General Level Shift Technique"
C   R.Carbo, J.A.Hernandez and F.Sanz,  Chem.Phys.Lett.  47 (1977) 581
C
C  DIIS
C  ----
C  "Improved SCF Convergence Acceleration"
C   P.Pulay,  J.Comput.Chem.  3 (1982) 556
C
C  Pseudodiagonalization
C  ---------------------
C  "Fast Semiempirical Calculations"
C   J.J.P.Stewart, P.Csaszar and P.Pulay,  J.Comput.Chem.  3 (1982) 227
C
C
      Parameter (nopt=10)
      Dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      Character*256 chopv(nopt)
      Character*4 opnames(nopt)
      Data opnames / 'semi','nogu','ethr','dthr','lvsh','diis',
     $               'iter','pseu','uhfs','prin'/
c
      Data ioptyp /21,0,11,11,11,11,1,11,0,1/
c
      CHARACTER ABPRT(3)*5
      DATA ABPRT/'     ','ALPHA',' BETA'/
      Data JUnit/15/, KUnit/16/         ! unit number for direct access files
c
      Common /job/jobname,lenJ
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
      PARAMETER (EVKCAL=23.0605423d0)
C
c      EVKCAL=RGETRVAL('evol')/RGETRVAL('kcal/mol')
C
C  initialize
C
      guess = .TRUE.      ! whether to use previously converged MOs
      EThr = 1.0d-9       ! convergence criterion on energy
      DThr = 1.0d-5       ! convergence criterion on density
      BShift = Zero       ! level shift
      Diis = 0.1d0        ! Diis switch on (maximum change in density)
      ItrMax = 200        ! maximum no. of SCF cycles
C
C  INPUT OPTIONS
C  -------------
      IF(inp.GT.0) THEN
        backspace inp         ! already read line in <presemi>
        call izeroit(iopv,3*nopt)
        call zeroit(ropv,3*nopt)
        Call readop1(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
c
        If(ifound(2).GT.0) guess = .FALSE.
        If(ifound(3).GT.0) EThr = ropv(1,3)
        If(ifound(4).GT.0) DThr = ropv(1,4)
        If(ifound(5).GT.0) BShift = ropv(1,5)
        If(ifound(6).GT.0) Diis = ropv(1,6)
        If(ifound(7).GT.0) ItrMax = iopv(1,7)
      ENDIF
c
      IErr = -1         ! assume the worst
      NITER=0
      EOLD=Zero
      READY=.FALSE.
      IOMIX=.FALSE.
      OKPULY=.TRUE.
      UHF=NBeta.GT.0
      NLin = (NBas*(NBas+1))/2
c
      IALP = 1
      JALP = 1
      IBET = 1
      JBET = 1
C
C  pointers for scratch
C
      K1 = 1
      K2 = K1 + NLin
      K3 = K1 + NBas**2
C
C  set up initial density matrix
C
      NEWDG = .FALSE.
      PL = One
      OKNEWD=ABS(BShift) .LT. 0.001D0
      JErr = 0
C
C ========================================================================
C -- for semiempirical proper, see if MOs from previous cycle available
      IF(guess) THEN
        call tstchval('mos-file',iyes)
        If(iyes.eq.1) Then
          call getchval('mos-file',MOS)
          If(NBeta.GT.0) call getchval('mob-file',MOB)
        Else
          MOS = jobname(1:lenJ)//'.mos'
          MOB = jobname(1:lenJ)//'.mob'
        EndIf
        call rmblan2(MOS,256,lenM)
        itype = 2
        call ReadMOS(NBas,C,jnk,.False.,lenM,MOS,itype,JErr)
        If(JErr.EQ.0) Then
          CALL DENSIT(NBas, NAlpha, C, Z(K1), PA)
          If(NBeta.GT.0) Then
            call rmblan2(MOB,256,lenM)
            call ReadMOS(NBas,CBETA,jnk,.False.,lenM,MOB,itype,JErr)
            If(JErr.EQ.0) Then
              CALL DENSIT(NBas, NBeta, CBETA, Z(K1), PB)
              DO 10 I=1,NLin
 10           P(I) = PA(I)+PB(I)
            EndIf
          Else
            DO 11 I=1,NLin
 11         P(I) = PA(I)+PA(I)
          EndIf
        EndIf
      ENDIF
C =========================================================================
C
      IF(.NOT.guess.OR.JErr.NE.0) THEN
        W1 = DBLE(NAlpha)/DBLE(NAlpha+NBeta)
        W2 = DBLE(NBeta)/DBLE(NAlpha+NBeta)
        IF(NBeta.EQ.0) THEN
         DO 12 I=1,NLin
         P(I) = Zero
         PA(I) = Zero
 12      CONTINUE
         DO 13 I=1,NBas
         J = (I*(I+1))/2
         P(J) = PSPD(I)
         PA(J) = P(J)*Half
 13      CONTINUE
        ELSE
         DO 14 I=1,NLin
         P(I) = Zero
         PA(I) = Zero
         PB(I) = Zero
 14      CONTINUE
         DO 15 I=1,NBas
         J = (I*(I+1))/2
         P(J) = PSPD(I)
         PA(J) = P(J)*W1
         PB(J) = P(J)*W2
 15      CONTINUE
        ENDIF
      ENDIF
cc
      CALL CpyVEC(NLin,PA,POLD)
      If(NBeta.GT.0) CALL CpyVEC(NLin,PB,PBOLD)
c
      If(UHF.AND.NAlpha.EQ.NBeta) IOMIX=.TRUE.
      If(.NOT.UHF) PLB=Zero
      IF(DThr.LT.EThr) DThr=EThr
c
      If(IPRNT.GT.5) Then
        WRITE(6,'(//10X,''ONE-ELECTRON MATRIX AT ENTRANCE TO ITER'')')
        Do I=1,NBas
        I1 = (I*(I-1))/2 + 1
        I2 = (I*(I+1))/2
        WRITE(6,1111) (H(J),J=I1,I2)
 1111   format(1X,6f12.6)
        EndDo
        WRITE(6,*) ' Initial Density matrix is:'
        Do I=1,NBas
        I1 = (I*(I-1))/2 + 1
        I2 = (I*(I+1))/2
        WRITE(6,1111) (P(J),J=I1,I2)
        EndDo
      ENDIF
c
      FRST = .TRUE.
      FRSTB = .TRUE.
      SHIFT = Zero
**********************************************************************
*                                                                    *
*                START THE SCF LOOP HERE                             *
*                                                                    *
**********************************************************************
   40 NITER=NITER+1
C
      IF(BShift.NE.Zero) THEN
        L=0
        SHIFT=BShift*(NITER+One)**(-1.5D0)
        DO 49 I=1,NBas
        DO 53 J=1,I
        L=L+1
  53    F(L)=H(L)+SHIFT*PA(L)
  49    F(L)=F(L)-SHIFT
      ELSE
        DO 50 I=1,NLin
   50   F(I)=H(I)
      ENDIF
C
C  construct the Fock matrix
C
  41  IF(UHF) THEN
        CALL FOCK2(NAtoms, SEMI,   NFIRST, NMIDLE, NLAST,
     $             IFACT,  I1FACT, W,      P,      PA,  F)
        CALL FOCK1(NAtoms, IAN,    NFIRST, NMIDLE, NLAST,
     $             P,      PA,     PB,     F)
        IF(SHIFT.NE.Zero) THEN
          L=0
          DO 58 I=1,NBas
          DO 59 J=1,I
          L=L+1
  59      FB(L)=H(L)+SHIFT*PB(L)
  58      FB(L)=FB(L)-SHIFT
        ELSE
          DO 57 I=1,NLin
  57      FB(I)=H(I)
        ENDIF
        CALL FOCK2(NAtoms, SEMI,   NFIRST, NMIDLE, NLAST,
     $             IFACT,  I1FACT, W,      P,      PB,  FB)
        CALL FOCK1(NAtoms, IAN,    NFIRST, NMIDLE, NLAST,
     $             P,      PB,     PA,     FB)
      ELSE
        CALL FOCK2(NAtoms, SEMI,   NFIRST, NMIDLE, NLAST,
     $             IFACT,  I1FACT, W,      P,      PA,  F)
        CALL FOCK1(NAtoms, IAN,    NFIRST, NMIDLE, NLAST,
     $             P,      PA,     PA,     F)
      ENDIF
c
      If(IPRNT.GT.5) Then
        WRITE(6,308) NITER
  308   FORMAT('   FOCK MATRIX ON ITERATION',I4)
        Do I=1,NBas
        I1 = (I*(I-1))/2 + 1
        I2 = (I*(I+1))/2
        WRITE(6,1111) (F(J),J=I1,I2)
        EndDo
        If(UHF) Then
          WRITE(6,309) NITER
  309     FORMAT('  BETA FOCK MATRIX ON ITERATION',I4)
          Do I=1,NBas
          I1 = (I*(I-1))/2 + 1
          I2 = (I*(I+1))/2
          WRITE(6,1111) (FB(J),J=I1,I2)
          EndDo
        EndIf
      EndIf
C
C  Test PL for SCF convergence
C  Exit if OK so that density matrices are those that give rise
C  to the Fock matrix, and not vice versa
C
      EE=HELECT(NBas,PA,H,F)
      IF(UHF) THEN
        SUM2=HELECT(NBas,PB,H,FB)
        EE=EE+SUM2
      ELSE
        EE=EE+EE
      ENDIF
      ESCF=(EE+ENuclr)*EVKCAL+ATHEAT
      DELTAE=ABS(ESCF-EOLD)
      IF(PL.LT.DThr.AND.
     +   ABS(ESCF-EOLD)/EVKCAL.LT.EThr .AND. READY) THEN
        If(ABS(SHIFT) .LT. 1.D-10) GOTO 63
        SHIFT=Zero
        DO 42 I=1,NLin
  42    F(I)=H(I)
        GO TO 41
      ENDIF
      READY=(ABS(ESCF-EOLD)/EVKCAL.LT.EThr*10.D0)
c
      If(IPRNT.GT.2) Then
        WRITE(6,1000) NITER,ESCF,ESCF-EOLD,MAX(PL,PLB)
 1000   FORMAT(' Iteration:',I4,' E: ',F15.7,' DeltaE: ',F14.7,
     $         ' DeltaD: ',D10.4)
      EndIf
c
      EOLD = ESCF
C
C  Diagonalize the Fock matrix
C
      IF(NEWDG) THEN
        If(OKPULY)
     +  CALL PULAY(NBas,   NLin,   PA,     F,      Z(K1),
     +             Z(K2),  BA,     JALP,   IALP,   MDiis,
     +             JUnit,  FRST,   IPRNT,  PL)
        CALL PseudoDIAG(NBas,NAlpha,EIGS,F,Z(K3),Z(K1),C)
      ELSE
        CALL EXPAND(NBas,F,C)
        CALL DIAGMAT(C,NBas,Z(K1),Z(K3),EIGS,JErr)
        If(JErr.NE.0) call nerror(1,'iter',
     $                            ' Diagonalization Failure!',0,0)
      ENDIF
c
      If(IPRNT.GT.3) Then
        J=1
        If(UHF)J=2
        WRITE(6,311) ABPRT(J),NITER,(EIGS(I),I=1,NBas)
 311    FORMAT(3X,A5,'  EIGENVALUES ON ITERATION',I3,/10(6F13.6,/))
      EndIf
      If(IPRNT.GT.4) Then
        J=1
        If(UHF)J=2
        WRITE(6,'(//3X,A5,
     + '' EIGENVECTORS ON ITERATION'',I3)') ABPRT(J),NITER
        CALL PrntMAT(NBas,NBas,NBas,C)
      EndIf
C
C *** COMPUTE THE BOND-ORDERS      ! what is this comment?  (JB)
C
      IF(UHF) THEN
        CALL DENSIT(NBas, NAlpha, C, Z(K1), PA)
        If(.NOT. (NEWDG.AND.OKPULY))
     +   CALL CNVG(NBas, PA, POLD, Z(K1), NITER, PL)
        IF(NEWDG.AND.OKPULY) THEN
          CALL PULAY(NBas,   NLin,   PB,     FB,     Z(K1),
     +               Z(K2),  BB,     JBET,   IBET,   MDiis,
     +               KUnit,  FRSTB,  IPRNT,  PLB)
          CALL PseudoDIAG(NBas,NBeta,EIGB,FB,Z(K3),Z(K1),CBETA)
        ELSE
          CALL EXPAND(NBas,FB,CBETA)
          CALL DIAGMAT(CBETA,NBas,Z(K1),Z(K3),EIGB,JErr)
          If(JErr.NE.0) call nerror(1,'iter',
     $                            ' Diagonalization Failure!',0,0)
        ENDIF
c
        If(IPRNT.GT.3) Then
          J=3
          WRITE(6,311) ABPRT(J),NITER,(EIGB(I),I=1,NBas)
        EndIf
        If(IPRNT.GT.4) Then
          J=3
          WRITE(6,'(//3X,A5,
     +   '' EIGENVECTORS ON ITERATION'',I3)') ABPRT(J),NITER
          CALL PrntMAT(NBas,NBas,NBas,CBETA)
        EndIf
c
        CALL DENSIT(NBas, NBeta, CBETA, Z(K1), PB)
        IF( .NOT. (NEWDG.AND.OKPULY))
     +       CALL CNVG(NBas, PB, PBOLD, Z(K1), NITER, PLB)
        DO 60 I=1,NLin
  60    P(I)=PA(I)+PB(I)
      ELSE
        CALL DENSIT(NBas, NAlpha, C, Z(K1), PA)
        DO 51 I=1,NLin
  51    P(I)=PA(I)+PA(I)
        IF(.NOT.(NEWDG.AND.OKPULY)) THEN
            CALL CNVG(NBas, P, POLD, Z(K1), NITER, PL)
        ENDIF
      ENDIF
c
      If(IPRNT.GT.5) Then
        WRITE(6,'('' DENSITY MATRIX ON ITERATION'',I4)') NITER
        Do I=1,NBas
        I1 = (I*(I-1))/2 + 1
        I2 = (I*(I+1))/2
        WRITE(6,1111) (P(J),J=I1,I2)
        EndDo
      EndIf
c
      OKNEWD = (PL.LT.EThr .OR. OKNEWD)
      NEWDG = (PL.LT.Diis .AND. OKNEWD .OR. NEWDG)
      If(PL.LT.Diis*0.3333D0) OKNEWD=.TRUE.
C
      IF(NITER.GE.ItrMax.AND.IPRNT.GT.0) THEN
        WRITE (6,260)
        WRITE (6,270) DELTAE,PL
  260   FORMAT (//2X,'***ERROR*** SCF Did Not Converge',/)
  270   FORMAT (/,5X,' DeltaE = ',E12.4,5X,' Delta Density = ',D10.4,/)
        RETURN
      ENDIF
      GO TO 40
**********************************************************************
*                                                                    *
*                    SELF-CONSISTENCY ACHIEVED                       *
*                      END THE SCF LOOP HERE                         *
*                  CALCULATE THE ELECTRONIC ENERGY                   *
*                                                                    *
**********************************************************************
C
  63  CONTINUE
      EE=HELECT(NBas,PA,H,F)
      IF(UHF) THEN
        EE=EE+HELECT(NBas,PB,H,FB)
      ELSE
        EE=EE+EE
      ENDIF
      IF( ABS(SHIFT) .GT. 1.D-5 ) THEN
C
C  Put F and FB into POLD in order not to destroy them
C  and do exact diagonalizations
C
        DO 72 I=1,NLin
   72   POLD(I)=F(I)
        CALL EXPAND(NBas,POLD,C)
        CALL DIAGMAT(C,NBas,Z(K1),Z(K3),EIGS,JErr)
        If(JErr.NE.0) call nerror(1,'iter',
     $                            ' Diagonalization Failure!',0,0)
        IF(UHF) THEN
          DO 65 I=1,NLin
   65     POLD(I)=FB(I)
          CALL EXPAND(NBas,POLD,CBETA)
          CALL DIAGMAT(CBETA,NBas,Z(K1),Z(K3),EIGS,JErr)
          If(JErr.NE.0) call nerror(1,'iter',
     $                            ' Diagonalization Failure!',0,0)
        ENDIF
      ENDIF
c
      IErr = 0
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE MAMULT(N,A,B,C,ONE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Multiplies matrix A by matrix B and puts the result in C
C  (matrices are packed, lower triangle)
C
      DIMENSION A(*),B(*),C(*)
C
      L=0
      DO 1 I=1,N
      II=((I-1)*I)/2
      DO 1 J=1,I
      JJ=((J-1)*J)/2
      L=L+1
      SUM=0.D0
      DO 2 K=1,J
   2  SUM=SUM+A(II+K)*B(JJ+K)
      DO 3 K=J+1,I
   3  SUM=SUM+A(II+K)*B(((K-1)*K)/2+J)
      DO 4 K=I+1,N
      KK=(K*(K-1))/2
   4  SUM=SUM+A(KK+I)*B(KK+J)
   1  C(L)=SUM+ONE*C(L)
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE MOLDAT(NAtoms, NBas,   IAN,    XC,     ICharg,
     $                  SEMI,   NFIRST, NMIDLE, NLAST,  USPD,
     $                  PSPD,   ATHEAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Sets up semiempirical common blocks and other data
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NBas    -  number of basis functions
C  IAN     -  list of atomic numbers
C  XC      -  Cartesian coordinates
C  ICharg  -  charge on system
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C
C  on exit
C
C  NFIRST
C  NMIDLE
C  NLAST
C  USPD    -  one-electron diagonal terms
C  PSPD    -  two-electron diagonal terms
C  ATHEAT  -  heat of atomization?
C
C
C -- common blocks
      COMMON /NATORB/ NATORB(54)
      COMMON /CORES / CORE(54)
      COMMON /BETAS / BETAS(54),BETAP(54),BETAD(54)
      COMMON /VSIPS / VS(54),VP(54),VD(54)
      COMMON /ONELEC/ USS(54),UPP(54),UDD(54)
      COMMON /MULTIP/ DD(54),QQ(54),AM(54),AD(54),AQ(54)
      COMMON /TWOELE/ GSS(54),GSP(54),GPP(54),GP2(54),HSP(54),
     1                GSD(54),GPD(54),GDD(54)
      COMMON /ALPHA / ALP(54)
      COMMON /GAUSS / FN1(54),FN2(54)
      COMMON /MNDO  / USSM(54),UPPM(54),UDDM(54),ZSM(54),ZPM(54),
     1                ZDM(54),BETASM(54),BETAPM(54),BETADM(54),
     2                ALPM(54),EISOLM(54),DDM(54),QQM(54),AMM(54),
     3                ADM(54),AQM(54),GSSM(54),GSPM(54),GPPM(54),
     4                GP2M(54),HSPM(54),POLVOM(54)
      COMMON /PM3 /  USSPM3(54),UPPPM3(54),UDDPM3(54),ZSPM3(54),
     $      ZPPM3(54),ZDPM3(54),BETASP(54),BETAPP(54),BETADP(54),
     $      ALPPM3(54),EISOLP(54),DDPM3(54),QQPM3(54),AMPM3(54),
     $      ADPM3(54),AQPM3(54),GSSPM3(54),GSPPM3(54),GPPPM3(54),
     $      GP2PM3(54),HSPPM3(54),POLVOP(54)
      COMMON /IDEAP / GUESP1(54,10),GUESP2(54,10),GUESP3(54,10)
     $       /IDEAS / GUESS1(54,10),GUESS2(54,10),GUESS3(54,10)
      COMMON /EXPONT/ ZS(54),ZP(54),ZD(54)
      COMMON /ATOMIC/ EISOL(54),EHEAT(54)
c -- common blocks for MINDO/3
      COMMON /ONELE3 /  USS3(18),UPP3(18)
      COMMON /ATOMI3 /  EISOL3(18),EHEAT3(18)
      COMMON /EXPON3 /  ZS3(18),ZP3(18)
c
      DIMENSION IAN(NAtoms),XC(3,NAtoms)
      DIMENSION NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms)
      Dimension USPD(NBas),PSPD(NBas)
      Character*20 SEMI
C
      Parameter (EVKCAL=23.0605423d0)
C
c      EVKCAL=RGETRVAL('evol')/RGETRVAL('kcal/mol')
C
C  determine what parameters to use
C  (Default set up appears to be AM1)
C
      IF(SEMI(1:3).EQ.'pm3') THEN
C
C  switch in PM3 parameters
C
        DO 10 I=1,54
        POLVOM(I)=POLVOP(I)
        ZS(I)=ZSPM3(I)
        ZP(I)=ZPPM3(I)
        ZD(I)=ZDPM3(I)
        USS(I)=USSPM3(I)
        UPP(I)=UPPPM3(I)
        UDD(I)=UDDPM3(I)
        BETAS(I)=BETASP(I)
        BETAP(I)=BETAPP(I)
        BETAD(I)=BETADP(I)
        ALP(I)=ALPPM3(I)
        EISOL(I)=EISOLP(I)
        DD(I)=DDPM3(I)
        QQ(I)=QQPM3(I)
        AM(I)=AMPM3(I)
        AD(I)=ADPM3(I)
        AQ(I)=AQPM3(I)
        GSS(I)=GSSPM3(I)
        GPP(I)=GPPPM3(I)
        GSP(I)=GSPPM3(I)
        GP2(I)=GP2PM3(I)
        HSP(I)=HSPPM3(I)
        GUESS1(I,1) = GUESP1(I,1)
        GUESS2(I,1) = GUESP2(I,1)
        GUESS3(I,1) = GUESP3(I,1)
        GUESS1(I,2) = GUESP1(I,2)
        GUESS2(I,2) = GUESP2(I,2)
        GUESS3(I,2) = GUESP3(I,2)
        DO 9 J=3,10
        GUESS1(I,J) = 0.0d0
        GUESS2(I,J) = 0.0d0
        GUESS3(I,J) = 0.0d0
  9     CONTINUE
 10     CONTINUE
cc
      ELSE IF(SEMI(1:4).EQ.'mndo') THEN
C
C  switch in MNDO parameters
C
        DO 20 I=1,54
        FN1(I)=0.D0
        ZS(I)=ZSM(I)
        ZP(I)=ZPM(I)
        ZD(I)=ZDM(I)
        USS(I)=USSM(I)
        UPP(I)=UPPM(I)
        UDD(I)=UDDM(I)
        BETAS(I)=BETASM(I)
        BETAP(I)=BETAPM(I)
        BETAD(I)=BETADM(I)
        ALP(I)=ALPM(I)
        EISOL(I)=EISOLM(I)
        GSS(I)=GSSM(I)
        GSP(I)=GSPM(I)
        GPP(I)=GPPM(I)
        GP2(I)=GP2M(I)
        HSP(I)=HSPM(I)
        DD(I)=DDM(I)
        QQ(I)=QQM(I)
        AM(I)=AMM(I)
        AD(I)=ADM(I)
        AQ(I)=AQM(I)
 20     CONTINUE
cc
      ELSE IF(SEMI(1:5).EQ.'mindo') THEN
C
C  switch in MINDO parameters
C
        DO 30 I=1,17
        USS(I)=USS3(I)
        UPP(I)=UPP3(I)
        EISOL(I)=EISOL3(I)
        EHEAT(I)=EHEAT3(I)
        GSS(I)=GSSM(I)
        GSP(I)=GSPM(I)
        GPP(I)=GPPM(I)
        GP2(I)=GP2M(I)
        HSP(I)=HSPM(I)
        ZS(I)=ZS3(I)
        ZP(I)=ZP3(I)
 30     CONTINUE
cc
      ENDIF
C
C  initialize
C
      ATHEAT=0.0d0
      EAT=0.0d0
      IA=1
      IB=0
      DO 190 II=1,NAtoms
        NFIRST(II)=IA
        NI=IAN(II)
        ATHEAT=ATHEAT+EHEAT(NI)
        EAT   =EAT   +EISOL(NI)
        IB=IA+NATORB(NI)-1
        NMIDLE(II)=IB
        IF(NATORB(NI).EQ.9) NMIDLE(II)=IA+3
        NLAST(II)=IB
        USPD(IA)=USS(NI)
        If(IA.EQ.IB) GOTO 183
        K=IA+1
        K1=IA+3
        DO 181 J=K,K1
        USPD(J)=UPP(NI)
 181    CONTINUE
 182    If(K1.EQ.IB) GO TO 183
        K=K1+1
        DO 184 J=K,IB
 184    USPD(J)=UDD(NI)
 183    CONTINUE
        IA=IB+1
 190  CONTINUE
C
      ATHEAT=ATHEAT-EAT*EVKCAL
C
      YY=FLOAT(ICharg)/FLOAT(NBas)
      L=0
      DO 191 I=1,NAtoms
      NI=IAN(I)
      XX=1.D0
      IF(NI.GT.2) XX=0.25D0
      W=CORE(NI)*XX-YY
      IA=NFIRST(I)
      IC=NMIDLE(I)
      IB=NLAST(I)
      DO 360 J=IA,IC
      L=L+1
  360 PSPD(L)=W
      DO 361 J=IC+1,IB
      L=L+1
  361 PSPD(L)=0.D0
  191 CONTINUE
C
cc      WRITE(6,3)(USPD(I),I=1,NBas)
cc   3  FORMAT('   ONE-ELECTRON DIAGONAL TERMS',/,10(/,10F8.3))
cc      WRITE(6,5)(PSPD(I),I=1,NBas)
cc   5  FORMAT('   INITIAL P FOR ALL ATOMIC ORBITALS',/,10(/,10F8.3))
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE PseudoDIAG(NBas,   NOcc,   EIG,    FAO,    WS,
     $                      FMO,    VECTOR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Fast pseudodiagonalization routine
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  NOcc    -  number of occupied MOs
C  EIG     -  eigenvalues from an exact diagonalization
C  FAO     -  Fock matrix, lower triangle, packed
C  WS      -  scratch storage
C  FMO     -  storage for occupied/virtual block of Fock matrix
C  VECTOR  -  on input - old eigenvectors
C             on exit  - new eigenvectors
C
C
C  This is a pseudo-diagonalization routine in that the vectors that
C  are generated are more nearly able to block-diagonalize the Fock
C  matrix over MOs than the starting vectors. Note that:
C  (1) No eigenvectors are generated - the secular determinant is
C      not diagonalized, only the occupied-virtual intersection.
C  (2) Many small elements in the secular determinant are ignored
C      as being too small compared with the largest element.
C  (3) When elements are eliminated by rotation, elements created
C      are ignored.
C
C  All the approximations used here become valid at self-consistency
C
C  reference
C  ---------
C  "Fast Semiempirical Calculations"
C   J.J.P.Stewart, P.Csaszar and P.Pulay  J.Comp.Chem. 3 (1982) 227
C
C
      DIMENSION FAO(*),EIG(NBas),VECTOR(NBas,*)
      Dimension WS(NBas),FMO((NBas-NOcc)*NOcc)
C
C
C  First construct the occupied-virtual block
C
      TINY=0.D0
      LUMO=NOcc+1
      IJ=0
      DO 6 I=LUMO,NBas
      KK=0
      DO 3 J=1,NBas
      SUM=0.D0
      DO 1 K=1,J
      KK=KK+1
   1  SUM=SUM+FAO(KK)*VECTOR(K,I)
      If(J.EQ.NBas) GO TO 3
c
      J1=J+1
      K2=KK
      DO 2 K=J1,NBas
      K2=K2+K-1
   2  SUM=SUM+FAO(K2)*VECTOR(K,I)
   3  WS(J)=SUM
      DO 5 J=1,NOcc
      IJ=IJ+1
      SUM=0.D0
      DO 4 K=1,NBas
   4  SUM=SUM+WS(K)*VECTOR(K,J)
      If(TINY.LT.ABS(SUM)) TINY=ABS(SUM)
   5  FMO(IJ)=SUM
   6  CONTINUE
      TINY=0.04D0*TINY
C
C  now do a crude 2 x 2 rotation to "eliminate" significant elements
C
      IJ=0
      DO 10 I=LUMO,NBas
      DO 9 J=1,NOcc
      IJ=IJ+1
      If(ABS(FMO(IJ)).LT.TINY) GO TO 9
C
C  begin 2 x 2 rotations
C
      A=EIG(J)
      B=EIG(I)
      C=FMO(IJ)
      D=A-B
      E=SIGN(SQRT(4.D0*C*C+D*D),D)
c -- be careful here, potential for round-off error   ! JB Feb 2004
c    at convergence, looks like D and E are the same, hence D/E --> 1
c    ALPHA --> 1, BETA --> 0. Round-off error may make (1 - D/E) very
c    small, but NEGATIVE, hence the SQRT blows up
      bchk = MIN(1.0d0,D/E)
      ALPHA=SQRT(0.5D0*(1.D0+bchk))
      BETA=-SIGN(SQRT(0.5D0*(1.D0-bchk)),C)
C
C  rotation of pseudo-eigenvectors
C
      DO 8 M=1,NBas
      A=VECTOR(M,J)
      B=VECTOR(M,I)
      VECTOR(M,J)=ALPHA*A+BETA*B
      VECTOR(M,I)=ALPHA*B-BETA*A
    8 CONTINUE
    9 CONTINUE
   10 CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE PULAY(NBas,   NDim,   P,      F,      FPPF,
     $                 FOCK,   EMAT,   LFOCK,  NFOCK,  MDiis,
     $                 IUnit,  START,  IPRNT,  PL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  DIIS SCF convergence
C  see   P. Pulay, J.Comp.Chem. 3 (1982) 556
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  NDim    -  dimension of Fock matrix [NBas*(NBas+1)/2]
C  P       -  density matrix, packed, lower triangle
C  F       -  Fock matrix, packed, lower triangle
C             on exit contains "best" Fock matrix
C  FPPF    -  work storage of size NDim
C  FOCK    -   ditto
C  EMAT    -  Diis error matrix
C  LFOCK   -  circular counter for location in direct access file
C  NFOCK   -  current size of DIIS subspace
C             (increase from 1 to MDiis then stays constant)
C  MDiis   -  maximum allowed size of DIIS subspace
C  IUnit   -  unit number of direct access file
C  START   -  logical flag to start DIIS
C  IPRNT   -  print flag
C  PL      -  F*P-P*F  (measure of self-consistency)
C
C
      DIMENSION F(NDim),P(NDim),FPPF(NDim),FOCK(NDim)
      DIMENSION EMAT(MDiis+1,MDiis+1), COEFFS(MDiis+1)
      Dimension EVEC((MDiis+1)**2),LL(MDiis+1),MM(MDiis+1)
      LOGICAL START,DEBUG
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
      DATA tol/1.0d-8/
C
      DEBUG=IPRNT.GT.4
      IF(START) THEN
        NFOCK=1
        LFOCK=1
        START=.FALSE.
      ELSE
        If(NFOCK.LT.MDiis)  NFOCK=NFOCK+1
        IF(LFOCK.NE.MDiis) THEN
          LFOCK=LFOCK+1
        ELSE
          LFOCK=1
        ENDIF
      ENDIF
C
C  First form /FOCK*DENSITY-DENSITY*FOCK/ and store in FPPF
C
      CALL MAMULT(NBas,P,F,FPPF,0.D0)
      CALL MAMULT(NBas,F,P,FPPF,-1.D0)
C
C  FPPF now contains the result of FP-PF
C  Save current fock matrix and FPPF to direct access file
C
      WRITE(IUnit,rec=LFOCK) F,FPPF
C
      NFOCK1=NFOCK+1
      DO 10 I=1,NFOCK
      EMAT(NFOCK1,I) = -One
      EMAT(I,NFOCK1) = -One
C
C  read back FPPF (into FOCK)
C
      READ(IUnit,rec=I) F,FOCK
c
      EMAT(LFOCK,I) = SProd(NDim,FOCK,FPPF)
      EMAT(I,LFOCK) = EMAT(LFOCK,I)
 10   CONTINUE
c
      PL=EMAT(LFOCK,LFOCK)/NDim
      EMAT(NFOCK1,NFOCK1) = Zero
c -- not sure why scaling of EMAT (below) is done
c -- perhaps to prevent determinant of matrix becoming too small?
      CONST = One/EMAT(LFOCK,LFOCK)
      DO 14 I=1,NFOCK
      DO 14 J=1,NFOCK
 14   EMAT(I,J) = EMAT(I,J)*CONST
c
      IF(DEBUG) THEN
        WRITE(6,'('' EMAT'')')
        DO 66 I=1,NFOCK1
 66     WRITE(6,'(6E13.6)') (EMAT(J,I),J=1,NFOCK1)
      ENDIF
c
      L=0
      DO 20 I=1,NFOCK1
      DO 20 J=1,NFOCK1
      L=L+1
  20  EVEC(L) = EMAT(I,J)
c -- undo previous scaling
      CONST = One/CONST
      DO 15 I=1,NFOCK
      DO 15 J=1,NFOCK
  15  EMAT(I,J) = EMAT(I,J)*CONST
C********************************************************************
C   The matrix EMAT should have the form
C
C      :<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0:
C      :<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0:
C      :<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0:
C      :<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0:
C      :     .            .      ...     . :
C      :   -1.0         -1.0     ...    0. :
C
C   where <E(I)*E(J)> is the scalar product of FP-PF for iteration I
C   times FP-PF for iteration J
C
C********************************************************************
      CALL OSINV(EVEC,NFOCK1,D,tol,LL,MM)
c
      If(IPRNT.GT.3)
     $   WRITE(6,'(''  Determinant of DIIS error matrix:'',E16.5)') D
C
C  If determinant greater than -1.0d-6, restart DIIS
C
      If(D.GT.-1.0d-6) Then
        START=.TRUE.
        RETURN
      EndIf
c
      If(NFOCK.LT.2) RETURN
      IL=NFOCK*NFOCK1
      DO 29 I=1,NFOCK
  29  COEFFS(I) = -EVEC(I+IL)
c
      IF(DEBUG) THEN
        WRITE(6,'('' EVEC'')')
        WRITE(6,'(6F12.6)')(COEFFS(I),I=1,NFOCK)
        WRITE(6,'(''    LAGRANGIAN MULTIPLIER (ERROR) =''
     +           ,F13.6)') EVEC(NFOCK1*NFOCK1)
      ENDIF
C
C  now form Fock matrix from linear combination of
C  previous Fock matrices
C
      CALL ZeroIT(F,NDim)
      DO 40 J=1,NFOCK
C
C  read back F (into FOCK)
C
      READ(IUnit,rec=J) FOCK,FPPF
c
      CVal = COEFFS(J)
      DO 30 I=1,NDim
      F(I) = F(I) + CVal*FOCK(I)
 30   CONTINUE
 40   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE REPP(NI,NJ,RIJ,RI,CORE)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
C
C
C  Calculates the two-electron repulsion integrals and the
C  nuclear attraction integrals
C
C  ARGUMENTS
C
C  NI      -  atomic number of first atom
C  NJ      -  atomic number of second atom
C  RIJ     -  interatomic distance
C
C  on exit
C
C  RI      -  array of two-electron repulsion integrals
C  CORE    -  4 x 2 array of electron-core attraction integrals
C
C  common blocks
C  -------------
C  DD      -  dipole charge separations
C  QQ      -  quadrupole separations
C  ADD     -  two-electron, one-centre integrals
C  TORE    -  nuclear charges
C
C
      DIMENSION RI(22),CORE(4,2)
      COMMON /MULTIP/ DD(54),QQ(54),ADD(54,3)
      COMMON /CORES/ TORE(54)
      COMMON /NATORB/ NATORB(54)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,RGas
C
C  initialize
C
      R=RIJ*ANTOAU
      PP=2.0D00
      P2=4.0D00
      P3=8.0D00
      P4=16.0D00
C
C  The two-centre repulsion integrals (over local coordinates)
C  are stored as follows (p-sigma=0; p-pi=p and p*):
C     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
C     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
C     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
C     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
C     (PP/P*P*)=21,   (P*P/P*P)=22.
C
C  Nuclear attraction integrals CORE(KL/IJ) are stored as:
C     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
C  where IJ=1 if the orbital is centered on atom I
C           2         "            "             J
C
      CALL ZeroIT(RI,22)
      CALL ZeroIT(CORE,8)
C
C  define charge separations
C
      DA=DD(NI)
      DB=DD(NJ)
      QA=QQ(NI)
      QB=QQ(NJ)
      TD = 2.D00
      OD = 1.D00
      FD = 4.D00
C
C  Hydrogen - Hydrogen
C
      AEE=0.25D00*(OD/ADD(NI,1)+OD/ADD(NJ,1))**2
      EE=OD/SQRT(R**2+AEE)
      RI(1)=EE*27.21D00
      CORE(1,1)=TORE(NJ)*RI(1)
      CORE(1,2)=TORE(NI)*RI(1)
      If(NATORB(NI).LT.3.AND.NATORB(NJ).LT.3) RETURN
      If(NATORB(NI).LT.3) GO TO 30
C
C  Heavy atom - Hydrogen
C
      ADE=0.25D00*(OD/ADD(NI,2)+OD/ADD(NJ,1))**2
      AQE=0.25D00*(OD/ADD(NI,3)+OD/ADD(NJ,1))**2
      DZE=-OD/SQRT((R+DA)**2+ADE)+OD/SQRT((R-DA)**2+ADE)
      QZZE=OD/SQRT((R-TD*QA)**2+AQE)-TD/SQRT(R**2+AQE)+OD/
     $         SQRT((R+TD*QA)**2+AQE)
      QXXE=TD/SQRT(R**2+FD*QA**2+AQE)-TD/SQRT(R**2+AQE)
      DZE=DZE/PP
      QXXE=QXXE/P2
      QZZE=QZZE/P2
      RI(2)=-DZE
      RI(3)=EE+QZZE
      RI(4)=EE+QXXE
      If(NATORB(NJ).LT.3) GO TO 40
C
C  Hydrogen - Heavy atom
C
   30 CONTINUE
      AED=0.25D00*(OD/ADD(NI,1)+OD/ADD(NJ,2))**2
      AEQ=0.25D00*(OD/ADD(NI,1)+OD/ADD(NJ,3))**2
      EDZ=-OD/SQRT((R-DB)**2+AED)+OD/SQRT((R+DB)**2+AED)
      EQZZ=OD/SQRT((R-TD*QB)**2+AEQ)-TD/SQRT(R**2+AEQ)+OD/
     $         SQRT((R+TD*QB)**2+AEQ)
      EQXX=TD/SQRT(R**2+FD*QB**2+AEQ)-TD/SQRT(R**2+AEQ)
      EDZ=EDZ/PP
      EQXX=EQXX/P2
      EQZZ=EQZZ/P2
      RI(5)=-EDZ
      RI(11)=EE+EQZZ
      RI(12)=EE+EQXX
      If(NATORB(NI).LT.3) GO TO 40
C
C  Heavy atom - Heavy atom
C  CAUTION. ADD REPLACES ADD(1,1) IN /MULTIP/ AND MUST BE RESET.
C
      ADD(1,1)=0.25D00*(OD/ADD(NI,2)+OD/ADD(NJ,2))**2
      ADQ=0.25D00*(OD/ADD(NI,2)+OD/ADD(NJ,3))**2
      AQD=0.25D00*(OD/ADD(NI,3)+OD/ADD(NJ,2))**2
      AQQ=0.25D00*(OD/ADD(NI,3)+OD/ADD(NJ,3))**2
      DXDX=TD/SQRT(R**2+(DA-DB)**2+ADD(1,1))
     $         -TD/SQRT(R**2+(DA+DB)**2+ADD(1,1))
      DZDZ=OD/SQRT((R+DA-DB)**2+ADD(1,1))
     $         +OD/SQRT((R-DA+DB)**2+ADD(1,1))-OD/SQRT((
     $         R-DA-DB)**2+ADD(1,1))-OD/SQRT((R+DA+DB)**2+ADD(1,1))
      DZQXX=-TD/SQRT((R+DA)**2+FD*QB**2+ADQ)+TD/SQRT((R-DA)**2
     $         +FD*QB**2+
     $         ADQ)+TD/SQRT((R+DA)**2+ADQ)-TD/SQRT((R-DA)**2+ADQ)
      QXXDZ=-TD/SQRT((R-DB)**2+FD*QA**2+AQD)+TD/SQRT((R+DB)**2
     $         +FD*QA**2+
     $         AQD)+TD/SQRT((R-DB)**2+AQD)-TD/SQRT((R+DB)**2+AQD)
      DZQZZ=-OD/SQRT((R+DA-TD*QB)**2+ADQ)+OD/SQRT((R-DA-TD*
     $   QB)**2+ADQ)-OD/SQRT((R+DA+TD*QB)**2+ADQ)+OD/SQRT((R-DA+TD*QB)
     $   **2+ADQ)-TD/SQRT((R-DA)**2+ADQ)+TD/SQRT((R+DA)**2+ADQ)
      QZZDZ=-OD/SQRT((R+TD*QA-DB)**2+AQD)+OD/SQRT((R+TD*QA+
     $   DB)**2+AQD)-OD/SQRT((R-TD*QA-DB)**2+AQD)+OD/SQRT((R-2.D
     $   00*QA+DB)**2+AQD)+TD/SQRT((R-DB)**2+AQD)-TD/SQRT((R+DB)**2
     $   +AQD)
      QXXQXX=TD/SQRT(R**2+FD*(QA-QB)**2+AQQ)+TD/SQRT(R**2+FD*(QA+QB)**2+
     $AQQ)-FD/SQRT(R**2+FD*QA**2+AQQ)-FD/SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT
     $(R**2+AQQ)
      QXXQYY=FD/SQRT(R**2+FD*QA**2+FD*QB**2+AQQ)-FD/SQRT(R**2+FD*QA**2+A
     $        QQ)-FD/SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
      QXXQZZ=TD/SQRT((R-TD*QB)**2+FD*QA**2+AQQ)+TD/SQRT((R+TD*QB)**2+FD*
     $QA**2+AQQ)-TD/SQRT((R-TD*QB)**2+AQQ)-TD/SQRT((R+TD*QB)**2+AQQ)-FD/
     $SQRT(R**2+FD*QA**2+AQQ)+FD/SQRT(R**2+AQQ)
      QZZQXX=TD/SQRT((R+TD*QA)**2+FD*QB**2+AQQ)+TD/SQRT((R-TD*QA)**2+FD*
     $QB**2+AQQ)-TD/SQRT((R+TD*QA)**2+AQQ)-TD/SQRT((R-TD*QA)**2+AQQ)-FD/
     $SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
      QZZQZZ=OD/SQRT((R+TD*QA-TD*QB)**2+AQQ)+OD/SQRT((R+TD*QA+TD*QB)**2+
     $AQQ)+OD/SQRT((R-TD*QA-TD*QB)**2+AQQ)+OD/SQRT((R-TD*QA+TD*QB)**2+AQ
     $Q)-TD/SQRT((R-TD*QA)**2+AQQ)-TD/SQRT((R+TD*QA)**2+AQQ)-TD/SQRT((R-
     $TD*QB)**2+AQQ)-TD/SQRT((R+TD*QB)**2+AQQ)+FD/SQRT(R**2+AQQ)
      DXQXZ=-TD/SQRT((R-QB)**2+(DA-QB)**2+ADQ)+TD/SQRT((R+QB)**2+(DA-QB)
     $**2+ADQ)+TD/SQRT((R-QB)**2+(DA+QB)**2+ADQ)-TD/SQRT((R+QB)**2+(DA+Q
     $B)**2+ADQ)
      QXZDX=-TD/SQRT((R+QA)**2+(QA-DB)**2+AQD)+TD/SQRT((R-QA)**2+(QA-DB)
     $**2+AQD)+TD/SQRT((R+QA)**2+(QA+DB)**2+AQD)-TD/SQRT((R-QA)**2+(QA+D
     $B)**2+AQD)
      QXYQXY=FD/SQRT(R**2+TD*(QA-QB)**2+AQQ)+FD/SQRT(R**2+TD*(QA+QB)**2+
     $AQQ)-8.D00/SQRT(R**2+TD*(QA**2+QB**2)+AQQ)
      QXZQXZ=TD/SQRT((R+QA-QB)**2+(QA-QB)**2+AQQ)-TD/SQRT((R+QA+QB)**2+(
     $QA-QB)**2+AQQ)-TD/SQRT((R-QA-QB)**2+(QA-QB)**2+AQQ)+TD/SQRT((R-QA+
     $QB)**2+(QA-QB)**2+AQQ)-TD/SQRT((R+QA-QB)**2+(QA+QB)**2+AQQ)+TD/SQR
     $T((R+QA+QB)**2+(QA+QB)**2+AQQ)+TD/SQRT((R-QA-QB)**2+(QA+QB)**2+AQQ
     $)-TD/SQRT((R-QA+QB)**2+(QA+QB)**2+AQQ)
c
      DXDX=DXDX/P2
      DZDZ=DZDZ/P2
      DZQXX=DZQXX/P3
      QXXDZ=QXXDZ/P3
      DZQZZ=DZQZZ/P3
      QZZDZ=QZZDZ/P3
      DXQXZ=DXQXZ/P3
      QXZDX=QXZDX/P3
      QXXQXX=QXXQXX/P4
      QXXQYY=QXXQYY/P4
      QXXQZZ=QXXQZZ/P4
      QZZQXX=QZZQXX/P4
      QZZQZZ=QZZQZZ/P4
      QXZQXZ=QXZQXZ/P4
      QXYQXY=QXYQXY/P4
      RI(6)=DZDZ
      RI(7)=DXDX
      RI(8)=-EDZ-QZZDZ
      RI(9)=-EDZ-QXXDZ
      RI(10)=-QXZDX
      RI(13)=-DZE-DZQZZ
      RI(14)=-DZE-DZQXX
      RI(15)=-DXQXZ
      RI(16)=EE+EQZZ+QZZE+QZZQZZ
      RI(17)=EE+EQZZ+QXXE+QXXQZZ
      RI(18)=EE+EQXX+QZZE+QZZQXX
      RI(19)=EE+EQXX+QXXE+QXXQXX
      RI(20)=QXZQXZ
      RI(21)=EE+EQXX+QXXE+QXXQYY
      RI(22)=0.5D0*(QXXQXX-QXXQYY)
      ADD(1,1)=ADD(1,2)
   40 CONTINUE
C
C  convert into EV
C
      DO 50 I=2,22
   50 RI(I)=RI(I)*27.21D00
C
C  calculate core-electron attraction
C
      CORE(2,1)=TORE(NJ)*RI(2)
      CORE(3,1)=TORE(NJ)*RI(3)
      CORE(4,1)=TORE(NJ)*RI(4)
      CORE(2,2)=TORE(NI)*RI(5)
      CORE(3,2)=TORE(NI)*RI(11)
      CORE(4,2)=TORE(NI)*RI(12)
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE ROTAT(NI,     NJ,     XI,     XJ,     SEMI,
     $                 W,      KR,     E1B,    E2A,    ENuc)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine computes the repulsion and nuclear attraction
C  integrals over molecular frame coordinates. Integrals over
C  local frame coordinates are first evaluated by subroutine REPP
C  and stored in RI as follows (p-sigma=0; p-pi=p and p*):
C     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
C     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
C     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
C     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
C     (PP/P*P*)=21,   (P*P/P*P)=22.
C
C
C  ARGUMENTS
C
C  NI      -  atomic number of first atom
C  NJ      -  atomic number of second atom
C  XI      -  Cartesian coordinates of first atom
C  XJ      -  Cartesian coordinates of second atom
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO, AM1 or PM3
C
C  on exit
C
C  W       -  two-electron repulsion integrals
C  KR      -  counter for number of two-electron integrals
C  E1B     -  electron-nuclear attraction integrals
C  E2A     -   ditto
C  ENuc    -  nuclear-nuclear repulsion term
C
C
      SAVE
      COMMON /NATORB/ NATORB(54)
      COMMON /TWOEL3/ F03(54)
      COMMON /ALPHA3/ ALP3(153)
      COMMON /ALPHA / ALP(54)
      COMMON /CORES / TORE(54)
      COMMON /IDEAS / FN1(54,10),FN2(54,10),FN3(54,10)
c
      DIMENSION XI(3),XJ(3),W(100),E1B(10),E2A(10)
      Dimension X(3),Y(3),Z(3),RI(22),CORE(4,2)
      Character*20 SEMI
      LOGICAL SI,SK,AM1
c --  the following is a dubious way of ensuring equivalences
      COMMON /ROTDU1/ CSS1,CSP1,CPPS1,CPPP1,CSS2,CSP2,CPPS2,CPPP2
      COMMON /ROTDU2/ X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3
      EQUIVALENCE (CORE(1,1),CSS1),(X(1),X1),(Y(1),Y1),(Z(1),Z1)
c -----------------------------------------------------------
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
C  initialize
C
      RIJ=Zero
      DO 15 I=1,3
      X(I)=XI(I)-XJ(I)
  15  RIJ=RIJ+X(I)**2
c
      IF(SEMI(1:5).EQ.'mindo') THEN
       SUM=14.399D0/SQRT(RIJ+(7.1995D0/F03(NI)+7.1995D0/F03(NJ))**2)
       W(1)=SUM
       KR=KR+1
       L=0
       DO 210 I=1,4
       DO 220 J=1,I
       L=L+1
       E1B(L)=Zero
  220  E2A(L)=Zero
       E1B(L)=-SUM*TORE(NJ)
  210  E2A(L)=-SUM*TORE(NI)
       II=MAX(NI,NJ)
       NBOND=(II*(II-1))/2+NI+NJ-II
       RIJ=SQRT(RIJ)
       IF(NBOND.GT.153) THEN
        SCALE=Zero
       ELSE
        If(NBOND.EQ.22 .OR. NBOND.EQ.29) Then
          SCALE = ALP3(NBOND)*EXP(-RIJ)
        Else
          SCALE = EXP(-ALP3(NBOND)*RIJ)
        EndIf
       ENDIF
       IF(ABS(TORE(NI)).GT.20.AND.ABS(TORE(NJ)).GT.20) THEN
        ENUC = Zero
       ELSE
        ENUC = TORE(NI)*TORE(NJ)*SUM +
     +         ABS(TORE(NI)*TORE(NJ)*(14.399D0/RIJ-SUM)*SCALE)
       ENDIF
       RETURN
cc
      ELSE
cc
       AM1 = SEMI(1:3).EQ.'am1'.OR.SEMI(1:3).EQ.'pm3'
       RIJ=SQRT(RIJ)
       CALL REPP(NI,NJ,RIJ,RI,CORE)
       GAM=RI(1)
C
C  The repulsion integrals over the molecular frame are stored in the
C  order in which they will later be used, i.e., (IJ|KL) where
C  J.LE.I and L.LE.K, with L varying most rapidly and I least rapidly
C
       A=One/RIJ
       DO 11 I=1,3
  11   X(I)=X(I)*A
       Z(3)=Zero
       IF(ABS(X(3)).GT.0.99999999D0) THEN
         X(3)=SIGN(1.D0,X(3))
         GO TO 12
       ENDIF
       Z(3)=SQRT(One-X(3)**2)
       A=One/Z(3)
       Y(1)=-A*X(2)*SIGN(One,X(1))
       Y(2)=ABS(A*X(1))
       Y(3)=Zero
       Z(1)=-A*X(1)*X(3)
       Z(2)=-A*X(2)*X(3)
       GO TO 13
  12   Y(1)=Zero
       Y(2)=One
       Y(3)=Zero
       Z(1)=One
       Z(2)=Zero
  13   CONTINUE
       IB=MIN(NATORB(NI),4)
       JB=MIN(NATORB(NJ),4)
       KI=0
       DO 130 I=1,IB
         SI=I.EQ.1
         II=I-1
       DO 130 J=1,I
         JJ=J-1
         IJ=0
         IF (JJ.EQ.0) IJ=-1
         IF (SI) IJ=+1
       DO 130 K=1,JB
         KK=K-1
         SK=KK.GT.0
       DO 130 L=1,K
         KI=KI+1
         If(SK) GO TO 50
C
C  integral (IJ|KL) is of type (IJ|SS)
C
         IF (IJ) 30,40,20
C
C  (SS|SS)
C
   20    W(KI)=RI(1)
         GO TO 131
C
C  (PS|SS)
C
   30    W(KI)=RI(2)*X(II)
         GO TO 131
C
C  (PP|SS)
C
   40    W(KI)=RI(3)*X(II)*X(JJ)+RI(4)*(Y(II)*Y(JJ)+Z(II)*Z(JJ))
         GO TO 131
   50    LL=L-1
         IF (LL.GT.0) GO TO 90
C
C  integral (IJ|KL) is of type (IJ|PS)
C
         IF (IJ) 70,80,60
C
C  (SS|PS)
C
   60    W(KI)=RI(5)*X(KK)
         GO TO 131
C
C  (PS|PS)
C
   70    W(KI)=RI(6)*X(II)*X(KK)+RI(7)*(Y(II)*Y(KK)+Z(II)*Z(KK))
         GO TO 131
C
C  (PP|PS)
C
   80    W(KI)=X(KK)*(RI(8)*X(II)*X(JJ)+RI(9)*(Y(II)*Y(JJ)+Z(II)*Z(JJ)))
     1   +RI(10)*(X(II)*(Y(JJ)*Y(KK)+Z(JJ)*Z(KK))+X(JJ)*(Y(II)*Y(KK)+Z(I
     2   I)*Z(KK)))
         GO TO 131
C
C  integral (IJ|KL) is of type (IJ|PP)
C
   90    IF (IJ) 110,120,101
C
C  (SS|PP)
C
  101    W(KI)=RI(11)*X(KK)*X(LL)+RI(12)*(Y(KK)*Y(LL)+Z(KK)*Z(LL))
         GO TO 131
C
C  (PS|PP)
C
  110    W(KI)=X(II)*(RI(13)*X(KK)*X(LL)+RI(14)*(Y(KK)*Y(LL)+Z(KK)*Z(LL)
     1   ))+RI(15)*(Y(II)*(Y(KK)*X(LL)+Y(LL)*X(KK))+Z(II)*(Z(KK)*X(LL)+Z
     2   (LL)*X(KK)))
         GO TO 131
C
C  (PP|PP)
C
  120    W(KI)=(RI(16)*X(II)*X(JJ)+RI(17)*(Y(II)*Y(JJ)+Z(II)*Z(JJ)))*X(K
     1   K)*X(LL)+RI(18)*X(II)*X(JJ)*(Y(KK)*Y(LL)+Z(KK)*Z(LL))+RI(19)*(Y
     2   (II)*Y(JJ)*Y(KK)*Y(LL)+Z(II)*Z(JJ)*Z(KK)*Z(LL))+RI(20)*(X(II)*(
     3   X(KK)*(Y(JJ)*Y(LL)+Z(JJ)*Z(LL))+X(LL)*(Y(JJ)*Y(KK)+Z(JJ)*Z(KK))
     4   )+X(JJ)*(X(KK)*(Y(II)*Y(LL)+Z(II)*Z(LL))+X(LL)*(Y(II)*Y(KK)+Z(I
     5   I)*Z(KK))))+RI(21)*(Y(II)*Y(JJ)*Z(KK)*Z(LL)+Z(II)*Z(JJ)*Y(KK)*Y
     6   (LL))+RI(22)*(Y(II)*Z(JJ)+Z(II)*Y(JJ))*(Y(KK)*Z(LL)+Z(KK)*Y(LL)
     7   )
  131  CONTINUE
  130  CONTINUE
  150  CONTINUE
       E1B(1)=-CSS1
       IF(NI.GT.2) THEN
        E1B(2) = -CSP1 *X1
        E1B(3) = -CPPS1*X1**2-CPPP1*(Y1**2+Z1**2)
        E1B(4) = -CSP1 *X2
        E1B(5) = -CPPS1*X1*X2-CPPP1*(Y1*Y2+Z1*Z2)
        E1B(6) = -CPPS1*X2*X2-CPPP1*(Y2*Y2+Z2*Z2)
        E1B(7) = -CSP1 *X3
        E1B(8) = -CPPS1*X1*X3-CPPP1*(Y1*Y3+Z1*Z3)
        E1B(9) = -CPPS1*X2*X3-CPPP1*(Y2*Y3+Z2*Z3)
        E1B(10)= -CPPS1*X3*X3-CPPP1*(Y3*Y3+Z3*Z3)
       ENDIF
       E2A(1)=-CSS2
       IF(NJ.GT.2) THEN
        E2A(2) = -CSP2 *X1
        E2A(3) = -CPPS2*X1**2-CPPP2*(Y1**2+Z1**2)
        E2A(4) = -CSP2 *X2
        E2A(5) = -CPPS2*X1*X2-CPPP2*(Y1*Y2+Z1*Z2)
        E2A(6) = -CPPS2*X2*X2-CPPP2*(Y2*Y2+Z2*Z2)
        E2A(7) = -CSP2 *X3
        E2A(8) = -CPPS2*X1*X3-CPPP2*(Y1*Y3+Z1*Z3)
        E2A(9) = -CPPS2*X2*X3-CPPP2*(Y2*Y3+Z2*Z3)
        E2A(10)= -CPPS2*X3*X3-CPPP2*(Y3*Y3+Z3*Z3)
       ENDIF
c
       If(ABS(TORE(NI)).GT.20.AND.ABS(TORE(NJ)).GT.20) Then
C -- sparkle-sparkle interaction (??)
        ENUC=Zero
        RETURN
       EndIf
c
       SCALE = EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
       NT=NI+NJ
       IF(NT.EQ.8.OR.NT.EQ.9) THEN
       If(NI.EQ.7.OR.NI.EQ.8) SCALE=SCALE+(RIJ-One)*EXP(-ALP(NI)*RIJ)
       If(NJ.EQ.7.OR.NJ.EQ.8) SCALE=SCALE+(RIJ-One)*EXP(-ALP(NJ)*RIJ)
       ENDIF
       ENUC = TORE(NI)*TORE(NJ)*GAM
       SCALE = ABS(SCALE*ENUC)
       IF(AM1) THEN
         DO 160 IG=1,10
         If(ABS(FN1(NI,IG)).GT.Zero)
     +      SCALE=SCALE +TORE(NI)*TORE(NJ)/RIJ*
     +      FN1(NI,IG)*EXP(-FN2(NI,IG)*(RIJ-FN3(NI,IG))**2)
         If(ABS(FN1(NJ,IG)).GT.Zero)
     +      SCALE=SCALE +TORE(NI)*TORE(NJ)/RIJ*
     +      FN1(NJ,IG)*EXP(-FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2)
  160    CONTINUE
       ENDIF
c
       ENUC = ENUC + SCALE
C
C *** NOW ROTATE THE NUCLEAR ATTRACTION INTEGRALS.
C *** THE STORAGE OF THE NUCLEAR ATTRACTION INTEGRALS  CORE(KL/IJ) IS
C     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
C
       KR=KR+KI
      ENDIF
C
      RETURN
      END
c =======================================================================
c
      DOUBLE PRECISION FUNCTION SPCG(NAtoms, IAN, SEMI, NFIRST, NLAST,
     $                               C1,     C2,  C3,   C4,     W)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  SPCG calculates the repulsion between electron 1 in MOs
C  C1 and C2 and electron 2 in MOs C3 and C4 for the valence
C  SP shell for MNDO or MINDO/3
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  SEMI    -  semiempirical Hamiltonian
C             currently either MINDO, MNDO or AM1
C  NFIRST
C  NLAST
C  C1      -  first MO of electron 1
C  C2      -  second MO of electron 1
C  C3      -  first MO of electron 2
C  C4      -  second MO of electron 2
C  W       -  two-electron matrix
C
C  on exit
C
C  SPCG    -  <C1(1)*C2(1)|C3(2)*C4(2)>
C
C
      DIMENSION C1(*),C2(*),C3(*),C4(*),W(*)
      COMMON /TWOELE/ GSS(54),GSP(54),GPP(54),GP2(54),HSP(54),
     $                GSD(54),GPD(54),GDD(54)
      Character*20 SEMI
C
C
      IF(SEMI(1:5).NE.'mindo') THEN
cc
       SPCG=0.0d0
       KK=0
       DO 9 II=1,NAtoms
       IA=NFIRST(II)
       IB=NLAST(II)
       IMINUS=II-1
       DO 10 JJ=1,IMINUS
       JA=NFIRST(JJ)
       JB=NLAST(JJ)
       DO 11 I=IA,IB
       DO 11 J=IA,I
       DO 11 K=JA,JB
       DO 11 L=JA,K
       KK=KK+1
       WINT=W(KK)
       SPCG=SPCG+WINT*(C1(I)*C2(J)*C3(K)*C4(L)
     .       + C1(K)*C2(L)*C3(I)*C4(J))
       If(I.NE.J) SPCG=SPCG+WINT*(C1(J)*C2(I)*C3(K)*C4(L)
     .                  + C1(K)*C2(L)*C3(J)*C4(I))
       If(K.NE.L) SPCG=SPCG+WINT*(C1(I)*C2(J)*C3(L)*C4(K)
     .                  + C1(L)*C2(K)*C3(I)*C4(J))
       If((I.NE.J).AND.(K.NE.L))SPCG=SPCG+WINT*(C1(J)*C2(I)*C3(L)*C4(K)
     .                  + C1(L)*C2(K)*C3(J)*C4(I))
   11  CONTINUE
   10  CONTINUE
    9  CONTINUE
       GOTO 301
cc
      ELSE
cc
       SPCG=0.0d0
       KR=0
       DO 210 II=1,NAtoms
       IA=NFIRST(II)
       IB=NLAST(II)
       IM1=II-1
       DO 220 JJ=1,IM1
       KR=KR+1
       ELREP=W(KR)
       JA=NFIRST(JJ)
       JB=NLAST(JJ)
       DO 240 I=IA,IB
       DO 240 K=JA,JB
  240  SPCG=SPCG+ELREP*(C1(I)*C2(I)*C3(K)*C4(K)+
     +                  C1(K)*C2(K)*C3(I)*C4(I))
  220  CONTINUE
  210  CONTINUE
      ENDIF
c
  301 CONTINUE
      ATEMP=SPCG
      IS1=0
      DO 3 I1=1,NAtoms
      IS1=IS1+1
      IZN=IAN(I1)
C
C  <SS|SS>
C
      SPCG=SPCG+C1(IS1)*C2(IS1)*C3(IS1)*C4(IS1)*GSS(IZN)
      IF(IZN.LT.3) GO TO 3
      IS=IS1
      IS1=IS1+1
       IX=IS1
      IS1=IS1+1
       IY=IS1
      IS1=IS1+1
       IZ=IS1
      SPCG=SPCG+GPP(IZN)*
     1                       (
C
C  <PP|PP> for P = X,Y,Z
C
     2   C1(IX)*C2(IX)*C3(IX)*C4(IX)+
     2   C1(IY)*C2(IY)*C3(IY)*C4(IY)+
     3   C1(IZ)*C2(IZ)*C3(IZ)*C4(IZ)
     4                       )
      SPCG=SPCG+GSP(IZN)*
     1                       (
C
C  <SS|PP>+<PP|SS> for P = X,Y,Z
C
     1   C1(IS)*C2(IS)*C3(IX)*C4(IX)+
     2   C1(IS)*C2(IS)*C3(IY)*C4(IY)+
     3   C1(IS)*C2(IS)*C3(IZ)*C4(IZ)+
     4   C1(IX)*C2(IX)*C3(IS)*C4(IS)+
     5   C1(IY)*C2(IY)*C3(IS)*C4(IS)+
     6   C1(IZ)*C2(IZ)*C3(IS)*C4(IS)
     7                       )
      SPCG=SPCG+GP2(IZN)*
     1                       (
C
C (PP|P,P,)+(P,P,|PP) for P.NE.P,=X,Y,Z
C
     1   C1(IX)*C2(IX)*C3(IY)*C4(IY)+
     2   C1(IX)*C2(IX)*C3(IZ)*C4(IZ)+
     3   C1(IY)*C2(IY)*C3(IZ)*C4(IZ)+
     4   C1(IY)*C2(IY)*C3(IX)*C4(IX)+
     5   C1(IZ)*C2(IZ)*C3(IX)*C4(IX)+
     6   C1(IZ)*C2(IZ)*C3(IY)*C4(IY)
     7                       )
      TEMP1=HSP(IZN)
      DO 4 J1=IX,IZ
      SPCG=SPCG+TEMP1*
     1                       (
C
C  (SP|SP)+(SP|PS)+(PS|SP)+(PS|PS) for P = X,Y,Z
C
     2   C1(IS)*C2(J1)*C3(J1)*C4(IS)+
     2   C1(IS)*C2(J1)*C3(IS)*C4(J1)+
     3   C1(J1)*C2(IS)*C3(IS)*C4(J1)+
     4   C1(J1)*C2(IS)*C3(J1)*C4(IS)
     5                       )
    4 CONTINUE
      TEMP1=0.5D0*(GPP(IZN)-GP2(IZN))
      DO 5 J1=IX,IZ
      DO 6 K1=IX,IZ
      IF(J1.EQ.K1) GO TO 6
      SPCG=SPCG+TEMP1*
     1                       (
C
C  (PP,|PP,)+(PP,|P,P)+(P,P|PP,)+(P,P|P,P) for P.NE.P,=X,Y,Z
C
     2   C1(J1)*C2(K1)*C3(J1)*C4(K1)+
     3   C1(J1)*C2(K1)*C3(K1)*C4(J1)
     4                       )
    6 CONTINUE
    5 CONTINUE
    3 CONTINUE
C
      RETURN
      END
c =======================================================================
c
      DOUBLE PRECISION FUNCTION DIPOLE
     $                         (NAtoms, IAN,    XC,     SEMI,   NFIRST,
     $                          P,      Q,      DIP,    DIPVEC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculates the dipole moment
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  NFIRST  -  start of atom orbital counters
C  XC      -  atomic coordinates
C  P       -  total density matrix
C  Q       -  total atomic charges (nuclear + electronic)
C  DIP     -  working storage
C
C  on exit
C
C  DIPVEC  -  dipole moment vector
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),NFIRST(NAtoms)
      DIMENSION P(*),Q(NAtoms),DIP(4,3),DIPVEC(3)
      Character*20 SEMI
      COMMON /MULTIP/ DD(54),QQ(54),AM(54),AD(54),AQ(54)
C
C ...................................................................
C  In the ZDO approximation, only two terms are retained in the
C  calculation of dipole moments:
C    1. The point charge term (independent of parametrization)
C    2. The one-center hybridization term, which arises from
C       matrix elements of the form <NS|R|NP>. This term is a
C       function of the Slater exponents (ZS,ZP) and is thus
C       dependent on parametrization.
C
C  The hybridization factors (HYF) used in this subroutine were
C  calculated from the following formulae:
C  For second row elements  <2S|R|2P>
C     HYF(I)= 469.56193322*(SQRT(((ZS(I)**5)*(ZP(I)**5)))/
C           ((ZS(I) + ZP(I))**6))
C  For third row elements  <3S|R|3P>
C     HYF(I)=2629.54682607*(SQRT(((ZS(I)**7)*(ZP(I)**7)))/
C           ((ZS(I) + ZP(I))**8))
C  For fourth row elements and higher
C     HYF(I)=2*(2.5416)*DD(I)
C  where DD(I) is the charge separation in atomic units.
C
C  References
C  ----------
C  J.A.Pople and D.L.Beveridge, Approximate MO Theory
C
C  S.P.McGlynn et al.  Applied Quantum Chemistry
C
C
      DIMENSION HYF(54,3)
C     MNDO DATA
      DATA   HYF(1,1)     / 0.0D00/
      DATA   HYF(3,1)/10.44577751D00/
C     MINDO/3 DATA
      DATA   HYF(1,2) /0.0D0     /
      DATA   HYF(5,2) /6.520587D0/
      DATA   HYF(6,2) /4.253676D0/
      DATA   HYF(7,2) /2.947501D0/
      DATA   HYF(8,2) /2.139793D0/
      DATA   HYF(9,2) /2.225419D0/
      DATA   HYF(14,2)/6.663059D0/
      DATA   HYF(15,2)/5.657623D0/
      DATA   HYF(16,2)/6.345552D0/
      DATA   HYF(17,2)/2.522964D0/
C     AM1 DATA
      DATA   HYF(1,3) /0.0D00       /
C
      PARAMETER (Debye=0.39342658d0)
C
c       Debye=rgetrval('deby')
C
C  modification (May 1990)
C
      DO 5 I=4,53
      HYF(I,1) = 5.0832d0*DD(I)
      HYF(I,3) = 5.0832d0*DD(I)
   5  CONTINUE
c
      If(SEMI(1:3).EQ.'am1') Then
        IType = 3
      Else If(SEMI(1:4).EQ.'mndo') Then
        IType = 2
      Else
        IType = 1
      EndIf
c
      CALL ZeroIT(DIP,12)
      DO 20 I=1,NAtoms
      NI = IAN(I)
      IA = NFIRST(I)
      DO 20 J=1,3
      K = ((IA+J)*(IA+J-1))/2+IA
      DIP(J,2) = DIP(J,2)-HYF(NI,ITYPE)*P(K)
   20 DIP(J,1) = DIP(J,1)+4.803d0*Q(I)*XC(J,I)
      DO 30 J=1,3
   30 DIP(J,3) = DIP(J,2)+DIP(J,1)
      DO 40 J=1,3
   40 DIP(4,J) = SQRT(DIP(1,J)**2+DIP(2,J)**2+DIP(3,J)**2)
c
      DIPVEC(1)=DIP(1,3)*Debye
      DIPVEC(2)=DIP(2,3)*Debye
      DIPVEC(3)=DIP(3,3)*Debye
      DIPOLE = DIP(4,3)*Debye
C
      RETURN
      END
