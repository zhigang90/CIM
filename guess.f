c ==================================================================
c  SCF GUESS MODULE                         JB    July  1997
c                                      revised    Nov   2000
c                             revised for ECPs    Sep   2004
c                 revised for atomic densities    April 2007
c                 added multiple orbital swaps    June  2008
c ==================================================================
c
      subroutine preguess(inp,gcard)
      implicit real*8(a-h,o-z)
      character*256 chopval,guess
c
c  reads the GUESS line in the input file and writes options
c  (if any) to the <control> file
c
c  ARGUMENTS
c
c  inp     -  unit number of input file
c  gcard   -  logical flag
c              set to .true. if a GUESS card has been found since the last SCF
c
c
      parameter (nopt=14)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character*20 cdum
      character*256 jobname
      logical found,gcard
c
      parameter (IUnit=1)
c
      data options/'gues','swap','swab','mix ','angl','uhfs',
     $             'file','prin','swa1','swa2','swa3','swb1',
     $             'swb2','swb3'/
      data ioptyp/21,2,2,2,11,0,21,1,2,2,2,2,2,2/
c
c
c -- make sure a basis has been defined
      call tstival('ncf',iyes)
      ncf = 0
      if(iyes.eq.1) ncf=igetival('ncf')
      if(iyes.eq.0.or.ncf.le.0) Call nerror(0,'GUESS module',
     $  'No Basis Set has been Defined - Check your Input',0,0)
c ........................................................
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      If(gcard) Then
        call readop1(inp,    nopt,   options,ioptyp, iopval,
     $               ropval, chopval,ifound)
      Else
c -- simulate a GUESS card
        call izeroit(ifound,nopt)
        ifound(1)=1
        chopval(1)='    '
c -- if an MO file exists, set GUESS=READ
        call getchval('jobname',jobname)
        call rmblan2(jobname,256,lenJ)
        inquire(file=jobname(1:lenJ)//'.mos',exist=found)
        if(found) chopval(1)='read'
      EndIf
c
c -- guess option
      if(ifound(1).eq.1) then
        if(chopval(1)(1:1).eq.' ') then
          guess = 'default'
        else
          guess = chopval(1)
        endif
      endif
c put this in the depository
      call setchval('guess',guess)
c
c -- default mixing angle
      if(ifound(5).eq.1) then
        angle = ropval(1,5)
      else
        angle = 30.0d0     ! default rotation angle for mixing
      endif
c
c -- print flag
      if(ifound(8).eq.1) then
        IPRNT = iopval(1,8)
      else
        IPRNT = 1          ! default print flag
      endif
c
c -- if found swap/mix options, but no orbitals given
c -- then we will swap/mix HOMO and LUMO
c
      if(ifound(2).eq.1.and.iopval(1,2).eq.0) iopval(1,2) = -1
      if(ifound(3).eq.1.and.iopval(1,3).eq.0) iopval(1,3) = -1
      if(ifound(4).eq.1.and.iopval(1,4).eq.0) iopval(1,4) = -1
c
c -- multiple swap options - check for sensible input
c   (must have pairs of orbitals with first greater than second)
c
      do i=9,14
      if(ifound(i).eq.1.and.( (iopval(1,i).le.0.and.iopval(2,i).le.0)
     $            .or.(iopval(1,i).ge.iopval(2,i)) ) )
     $  Call nerror(17,'GUESS module',
     $  'Orbitals Wrong for Swap - Check your Input',0,0)
      end do
c
c -- file option - if the name of the MOS file is not <jobname>.mos
c    (and <jobname>.mob ). The file name for the beta MOS must be the
c    same as for the alpha MOS except the last character of the
c    extension is replaced by a "b" (lower case b)
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
      if(ifound(7).eq.1) then
        call setchval('mosfname',chopval(7))
        call rmblan2(chopval(7),256,lenM)
        call CopyFile(chopval(7)(1:lenM)//'.basis',
     $        lenM+6,jobname(1:lenJ)//'.basis2',lenM+7)
      else
        call setchval('mosfname',jobname)
      end if
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call wrguess(IUnit,guess,ifound(6),iopval,angle)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c -- if guess=read make a copy of the basis file
c -- if not already done so
c
      If(guess(9:12).eq.'read') Then
        inquire(file=jobname(1:lenJ)//'.basis2',exist=found)
        if(.not.found) Call CopyFile(jobname(1:lenJ)//'.basis',
     $          lenJ+6,jobname(1:lenJ)//'.basis2',lenJ+7)
      EndIf
c
      return
      end
c  =======================================================================
c
      SUBROUTINE GUESS(NMem,Z,natom)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** SELF-CONSISTENT FIELD (SCF) GUESS PROGRAM **
C
C  This Program generates initial guess MOs for either a
C  closed-shell SCF or open-shell UHF wavefunction.
C
C  There are 4 options:
C  (1)  Use MOs from previous calculation
C  (2)  Semiempirical guess
C        (i) PM3  (ii) AM1  (iii) MNDO  (iv) MINDO/3
C  (3)  Extended Huckel guess
C  (4)  Core guess
C
C  If the first three guesses fail, the standard default (core)
C  guess is used, and the initial MOs are generated in the SCF
C  module itself.
C
C  Alternatively the density matrix can be generated directly
C  by superposition of atomic densities.
C  ..........................................................
C
C  GUESS itself is simply a "wrapper" which reads the
C  <control> file to determine the size of the current
C  system and the guess option requested and, based on
C  this information, allocates the necessary memory.
C  It then calls GUESSMAIN which is itself another "wrapper"
C  which completes job input and is responsible for all
C  other I/O. GUESSMAIN calls <scfguess> which is the main
C  driving routine. More memory may be requested later by
C  other routines called by <scfguess>, but the bulk of the
C  storage space required is allocated here.
C  .............................................................
C
C
      CHARACTER*256 jobname
      CHARACTER*20 GUES,wvfnc,SEMI,cdum
      LOGICAL Symflag
C
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(natom)
c ..................................................
      DIMENSION Z(NMem)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
      wvfnc = 'ab initio'               ! set up default
      Symflag = .FALSE.                 ! temporary
C
C  Read from the <control> file
C    number of atoms
C    type of SCF guess requested
C    number of alpha/closed-shell orbitals
C    number of beta orbitals
C    wavefunction type (i.e., ab initio or semiempirical)
C    if we have a UHF-singlet
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
      call rdcntrl(IUnit,6,'$guess',3,idum,dum,GUES)
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      call fdcntrl(IUnit,9,'$wavefunc',IEnd)
      If(IEnd.EQ.0) call rdcntrl(IUnit,9,'$wavefunc',3,idum,dum,wvfnc)
      If(wvfnc(1:7).EQ.'Semiemp') SEMI = wvfnc(9:13)
      call fdcntrl(IUnit,12,'$uhf_singlet',IEnd)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c   Use existing guess "as is"
      if(GUES(1:3).eq.'use') RETURN
c
c ............................................................
c -- for the MO guess, dummy atoms are ignored
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
          Call nerror(1,'GUESS module',
     $         'UHF Singlet requested and System is Not a Singlet!',
     $          0,0)
        EndIf
      EndIf
c
c .............................................................
c -- Special. If GUESS=CORE use old default guess in <scfmain>
c
      If(GUES(1:4).EQ.'core') RETURN
c .............................................................
C
C  Read from the <basis> file
C    number of basis functions, primitives and shells
C    number of pseudopotentials
C
      NPsp1 = 0
      NPsp2 = 0
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$nbasis',1,NBas1,dum,cdum)
      call rdcntrl(IUnit,7,'$nshell',1,NShl1,dum,cdum)
      call rdcntrl(IUnit,6,'$nprim',1,NPrm1,dum,cdum)
      call fdcntrl(IUnit,7,'$npseud',IEnd)
      If(IEnd.EQ.0)
     $  call rdcntrl(IUnit,7,'$npseud',1,NPsp1,dum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c .............................................................
c -- Special. If GUESS=ATOM do NOT generate MOs - use density
c
      If(GUES(1:4).EQ.'atom') Then
        CALL ATOM_Density(NAtoms, NBas1,  NShl1,  NPrm1,  NPsp1,
     $                    NBeta,  AtSymb, NMem,   Z,      IErr)
        If(IErr.NE.0) Call nerror(9,'GUESS module',
     $                  'Atomic Density guess failed',0,0)
        If(NBeta.GT.0) Then
          OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $          FORM='FORMATTED',STATUS='OLD')
          Call wrcntrl(IUnit,6,'$nbeta',1,NBeta,rdum,cdum)
          CLOSE (UNIT=IUnit,STATUS='KEEP')
        EndIf
        RETURN
      EndIf
c .............................................................
C
C  Determine these quantities for the guess basis
C
      IF(GUES(1:4).EQ.'read'.AND.wvfnc(1:7).NE.'Semiemp') THEN
C
C  Read from the <basis2> file
C    number of basis functions, primitives and shells
C    number of pseudopotentials
C
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis2',
     $        FORM='FORMATTED',STATUS='OLD',ERR=95)
        call rdcntrl(IUnit,7,'$nbasis',1,NBas2,dum,cdum)
        call rdcntrl(IUnit,7,'$nshell',1,NShl2,dum,cdum)
        call rdcntrl(IUnit,6,'$nprim',1,NPrm2,dum,cdum)
        call fdcntrl(IUnit,7,'$npseud',IEnd)
        If(IEnd.EQ.0)
     $    call rdcntrl(IUnit,7,'$npseud',1,NPsp2,dum,cdum)
        CLOSE(UNIT=IUnit,STATUS='KEEP')
c
        NBas = 0
        N2Elec = 0
        NLin = 0
      ELSE
C
C  determine values for a minimal basis
C
        IXC = 1
        ICHG = IXC + 3*NAtoms
        IXM = ICHG + NAtoms
        IAN = IXM + NAtoms
c
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $        FORM='FORMATTED',STATUS='OLD',ERR=95)
        Call RdCoordF(IUnit,  NAtoms, AtSymb, Z(IXC), -1,
     $                jnk,    Z(ICHG),Z(IXM))
        CLOSE (UNIT=IUnit,STATUS='KEEP')
cc
        CALL GetAtNo(NAtoms,AtSymb,Z(IAN))
        CALL GetMinBas(NAtoms,Z(IAN),NBas2,NShl2,NPrm2)
C
C  determine values for a semiempirical guess
C
        CALL GetNHeavy(NAtoms,Z(IAN),NHeavy)
C
C  Get the number of ghost atoms
C
        call GetNGhosts(NAtoms,Z(IAN),Z(ICHG),NHGhosts,NLGhosts)
c
        NLight = NAtoms-NHeavy-NLGhosts
        NHeavy = NHeavy-NHGhosts
        NBas = 4*NHeavy + NLight
        N2Elec = 50*NHeavy*(NHeavy-1) + 10*NHeavy*NLight
     $           + (NLight*(NLight-1))/2
        NLin = (NBas*(NBas+1))/2
      ENDIF
C
C  determine amount of "scratch" memory
C  (this depends on SCF guess)
C
      MDiis = 6                      ! maximum size of Diis subspace
      NSemi = 3*NAtoms + 4*NBas + N2Elec + 7*NLin + NBas**2
     $          + (MDiis+1)**2
      If(NBeta.GT.0) NSemi = NSemi + NBas + 3*NLin + NBas**2
     $                          + (MDiis+1)**2
      NHuck = 2*NBas2**2 + NBas2
      NSemi = MAX(NSemi,127008)
      NHuck = MAX(NHuck,127008)
      NRead = MAX(127008,3*NBas1*(NBas1+1))
c
      NScr = NBas1*NBas1 + NBas2*NBas2 + NBas1*NBas2 + NBas2
c
      NScr = MAX(NSemi,NHuck,NScr)
      If(GUES(1:4).EQ.'read') NScr = NRead
C
C
C  Now get the memory
C
      IMem= 9*NAtoms + 12*(NShl1+NShl2) + 13*(NPrm1+NPrm2) + NBas1*NBas2
     $                + 2*NBas1**2 + NBas2**2 + NBas2 + 5 + NScr
      If(NBeta.GT.0) IMem = IMem + NBas1**2 + NBas2**2 + NBas2
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,5,'GUESS')
      CALL ZeroIT(Z,IMem)               ! clear the memory
C
C
C  Allocate memory pointers
C
      IAN   = iptr                      !  atomic numbers
      IXC   = IAN   + NAtoms            !  geometry
      ICHG  = IXC   + 3*NAtoms          !  atomic charges
      IXM   = ICHG  + NAtoms            !  atomic masses
      IBAS1 = IXM   + NAtoms            !  basis function data
      INX1  = IBAS1 + 13*NPrm1          !  basis indexing array
      IPSP1 = INX1  + 12*NShl1          !  pseudopotential array
      IBAS2 = IPSP1 + NAtoms            !  guess basis function data
      INX2  = IBAS2 + 13*NPrm2          !  guess basis indexing array
      IPSP2 = INX2  + 12*NShl2          !  guess pseudopotential array
      IS11  = IPSP2 + NAtoms            !  basis/basis overlap
      IS12  = IS11  + NBas1*NBas1       !  basis/guess-basis overlap
c --  S12 used as scratch array in Huckel, add 5 more locations --
      ICA1  = IS12  + NBas1*NBas2 + 5   !  actual basis alpha MO coeffs
      ICA2  = ICA1  + NBas1*NBas1       !  guess basis alpha MO coeffs
      IEA   = ICA2  + NBas2*NBas2       !  guess orbital energies
      IEnd  = IEA   + NBas2
c
      IF(NBeta.GT.0) THEN
       ICB1 = IEnd                      !  actual basis beta MO coeffs
       ICB2 = ICB1 + NBas1*NBas1        !  guess basis beta MO coeffs
       IEB  = ICB2 + NBas2*NBas2        !  guess beta orbital energies
       IEnd = IEB  + NBas2
      ELSE
       ICB1 = 1
       ICB2 = 1
       IEB  = 1
      ENDIF
C
C  storage for symmetry
C
      IUQ = IEnd                      !  list of symmetry unique atoms
      IEnd = IUQ + NAtoms
c
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
C
C  general scratch storage
C
      IScr = IEnd
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,5,'GUESS')
c
c -- set memory mark for TEXAS memory manager
      call getmem(IMem,lastx)
C
C
C  ----------------------------------------------------------------------
C
      CALL GUESSMAIN(NAtoms,  GUES,    AtSymb,  Z(IAN),  Z(IXC),
     $               Z(ICHG), Z(IXM),  NBas1,   NShl1,   NPsp1,
     $               Z(IBAS1),Z(INX1), Z(IPSP1),NBas2,   NShl2,
     $               NPsp2,   Z(IBAS2),Z(INX2), Z(IPSP2),Z(IS11),
     $               Z(IS12), NAlpha,  NBeta,   Z(ICA1), Z(ICA2),
     $               Z(IEA),  Z(ICB1), Z(ICB2), Z(IEB),  NBas,
     $               N2Elec,  MDiis,   NTrans,  Z(IUQ),  Z(ITN),
     $               Z(INQ),  wvfnc,   SEMI,    NScr,    Z(IScr), IErr)
C
C  ----------------------------------------------------------------------
C
C  free memory used in this routine
C
        call retmem(1)
C
C  Exit procedure
C
      RETURN
C  ...............................................
C -- Error Section
C
 95   CONTINUE
      Call nerror(2,'GUESS module',
     $     'GUESS=READ Specified and old <basis> File Does Not Exist',
     $      0,0)
c
      END
c  =======================================================================
c
      SUBROUTINE GUESSMAIN(NAtoms, GUESS,  AtSymb, IAN,    XC,
     $                     XCharg, XMass,  NBas1,  NShl1,  NPsp1,
     $                     BASDAT1,INX1,   IPSP1,  NBas2,  NShl2,
     $                     NPsp2,  BASDAT2,INX2,   IPSP2,  S11,
     $                     S12,    NAlpha, NBeta,  CMO1,   CMO2,
     $                     EA,     CBeta1, CBeta2, EB,     NBas,
     $                     N2Elec, MDiis,  NTrans, IUQ,    TRANS,
     $                     NEqATM, wvfnc,  SEMI,   NMem,   Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for SCF guess program
C  GUESSMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),XCharg(NAtoms),XMass(NAtoms),
     $          BASDAT1(13,*),INX1(12,*),BASDAT2(13,*),INX2(12,*),
     $          IPSP1(NAtoms),IPSP2(NAtoms)
      DIMENSION S11(NBas1,NBas1),S12(NBas1,NBas2),CMO1(NBas1,NBas1),
     $          CMO2(NBas2,NBas2),EA(NBas2),CBeta1(NBas1,NBas1),
     $          CBeta2(NBas2,NBas2),EB(NBas2)
      DIMENSION IUQ(NAtoms),TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      CHARACTER*8 AtSymb(NAtoms),GUESS
      Character*256 jobname,scrdir,MOS,MOB
      Character*20 wvfnc,SEMI,cdum
      INTEGER swap1,swap2,swab1,swab2
      INTEGER MSwap(2,6)
      Logical EPrnt
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
C  initialize
C
      IOut=ioutfil('iout')
      swap1 = 0
      swap2 = 0
      swab1 = 0
      swab2 = 0
      mix1  = 0
      mix2  = 0
c -- multiple swaps (actually up to 3) now a possibity
      call IZeroIT(MSwap,12)
      EPrnt = .True.
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
C     total molecular charge
C     print flag
C     whether orbitals are to be swapped
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$charge',2,idum,Charg,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
C
C  are orbitals to be swapped?
C
      call fdcntrl(IUnit,5,'$swap',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1000) swap1,swap2
      EndIf
c
      call fdcntrl(IUnit,10,'$beta_swap',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1100) swab1,swab2
      EndIf
C
C  are alpha/beta orbitals to be mixed?
C
      call fdcntrl(IUnit,4,'$mix',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1200) mix1,mix2,angle
      EndIf
C
C  check for multiple swaps
C
      call fdcntrl(IUnit,5,'$swp1',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1000) MSwap(1,1),MSwap(2,1)
      EndIf
      call fdcntrl(IUnit,5,'$swp2',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1000) MSwap(1,2),MSwap(2,2)
      EndIf
      call fdcntrl(IUnit,5,'$swp3',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1000) MSwap(1,3),MSwap(2,3)
      EndIf
c
      call fdcntrl(IUnit,10,'$beta_swp1',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1100) MSwap(1,4),MSwap(2,4)
      EndIf
      call fdcntrl(IUnit,10,'$beta_swp2',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1100) MSwap(1,5),MSwap(2,5)
      EndIf
      call fdcntrl(IUnit,10,'$beta_swp3',IEnd)
      If(IEnd.EQ.0) Then
        backspace IUnit
        READ(IUnit,1100) MSwap(1,6),MSwap(2,6)
      EndIf
c
      CLOSE(UNIT=IUnit,STATUS='KEEP')
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
     $             Z(1+NAtoms),BASDAT1)
C  ...................................................................
C  Are there pseudopotentials?
C
      If(NPsp1.GT.0) Then
        CALL IZeroIT(IPSP1,NAtoms)
        call rdcntrl(IUnit,7,'$npseud',1,NPsp1,rdum,cdum)
        Do I=1,NPsp1
        READ(IUnit,*) NAtm,Izcore
        IPSP1(NAtm) = Izcore
        EndDO
      EndIf
C  ...................................................................

      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  sort the basis (by shell)
C
      call SortBAS1(NAtoms,NShl1,Z(1+NAtoms),INX1)
cc      write(6,*) ' INX1 array after <SortBAS1> is:'
cc      do j=1,NShl1
cc      write(6,*) (inx1(i,j),i=1,12)
cc      enddo
C
C  Get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
C  ...................................................................
C  If GUESS=READ option requested, read old basis and old MOS
C  from existing files
C
      IF(GUESS(1:4).EQ.'read') THEN
C
C  MOs will be read from file
C  Either ab initio Or semiempirical
C
        IF(wvfnc(1:7).NE.'Semiemp') THEN
C
C  Read the old basis set data
C
          OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis2',
     $          FORM='FORMATTED',STATUS='OLD')
          call rdbasis(IUnit,  NAtoms, AtSymb, XC,     Z(1),
     $                 Z(1+NAtoms),BASDAT2)
C  ...................................................................
C  Are there pseudopotentials?
C
          If(NPsp2.GT.0) Then
            CALL IZeroIT(IPSP2,NAtoms)
            call rdcntrl(IUnit,7,'$npseud',1,NPsp2,rdum,cdum)
            Do I=1,NPsp2
            READ(IUnit,*) NAtm,Izcore
            IPSP2(NAtm) = Izcore
            EndDO
C
C  If pseudopotentials do not match, stop
C
            CALL ChkPseud(NAtoms,IPSP1,IPSP2,IErr)
            If(IErr.NE.0) Call nerror(3,'GUESS nodule',
     $      'Mismatch in number of ECPs in current and guess basis',0,0)
          EndIf
C  ...................................................................
          CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  sort the basis (by shell)
C
          call SortBAS1(NAtoms,NShl2,Z(1+NAtoms),INX2)
        ENDIF
C
C  now read the old MOs (from the existing binary file)
C
        NBasR = NBas2
        If(wvfnc(1:7).EQ.'Semiemp') NBasR = NBas
c
        itype = 0
        call ReadMOS(NBasR,  CMO2,   jnk,    .False.,lenM,
     $               MOS(1:lenM),itype,IErr)
        If(NBeta.GT.0) Then
          call ReadMOS(NBasR,  CBeta2, jnk,    .False.,lenM,
     $                 MOB(1:lenM),itype,IErr)
          If(IErr.NE.0) Then
           call ReadMOS(NBasR,  CBeta2, jnk,    .False.,lenM,
     $                  MOS(1:lenM),itype,IErr)
           call message('**WARNING** in GUESS module',
     $                  '  Reading Closed-Shell MOs as Beta',0,0)
          EndIf
        EndIf
c
        If(IErr.NE.0) Call nerror(4,'GUESS module',
     $     'GUESS=READ Specified and old MOS File Does Not Exist',0,0)
cc
        IF(wvfnc(1:7).EQ.'Semiemp') THEN
C
C  assign memory pointers
C
          INF = 1
          IC  = INF + NBas2
          IEnd = IC + NBas2**2 - 1
          CALL MemCHK(NMem,IEnd,5,'GUESS')
C
C  remove ghost atoms
C
          call RemoveGhosts(NAtoms,AtSymb,IAN,XCharg,XC)
C
C  set up appropriate semiempirical exponents and core basis
C
          CALL SemiEXPONT(SEMI)
          NBas2 = NBas
          CALL SemiCORE(NAtoms, IAN,    XC,     NBas2,   NShl2,
     $                  BASDAT2,INX2,   IPSP2,  Z(INF))
C
C  get right number of alpha and beta electrons for semiempirical
C
          NAlphS = NAlpha
          NBetS = NBeta
c
          CALL ChkSEMI(SEMI,NAtoms,IAN,NAlphS,NBetS,IErr)
C
C  add core orbitals to semiempirical MOS
C
          nmo = NAlphS
          If(MAX(swap2,swab2,mix2).gt.0.OR.IPRNT.GE.2) nmo=NBas
c
          CALL CpyVEC(NBas*nmo,CMO2,Z(IC))
          CALL AddCORE(NBas,NBas2,nmo,Z(IC),Z(INF),CMO2)
          If(NBeta.GT.0) Then
            CALL CpyVEC(NBas*nmo,CBeta2,Z(IC))
            CALL AddCORE(NBas,NBas2,nmo,Z(IC),Z(INF),CBeta2)
          EndIf
        ENDIF
      ENDIF
C  ...................................................................
C
C  Make the Guess!
C
C  ------------------------------------------------------------
      CALL SCFGUESS(NAtoms, GUESS,  IPRNT,  AtSymb, IAN,
     $              XC,     XCharg, NBas1,  NShl1,  NPsp1,
     $              BASDAT1,INX1,   IPSP1,  NBas2,  NShl2,
     $              NPsp2,  BASDAT2,INX2,   IPSP2,  S11,
     $              S12,    NAlpha, NBeta,  CMO1,   CMO2,
     $              EA,     CBeta1, CBeta2, EB,     Charg,
     $              NBas,   N2Elec, MDiis,  swap1,  swap2,
     $              swab1,  swab2,  mix1,   mix2,   angle,
     $              MSwap,  wvfnc,  NMem,   Z,      Z,   IErr)
C  ------------------------------------------------------------
C
C
      If(IErr.EQ.0) Then
C
C  write guess MOs to binary <MOS> file
C  First restore the proper name of the MOS file
C
       MOS=jobname(1:lenJ)//'.mos'
       lenM=lenJ+4
       MOB=MOS
       MOB(lenM:lenM)='b'
       call setchval('mos-file',MOS)
       call setchval('mob-file',MOB)
C
C    Remove the basis2 file
C
       open(IUnit,file=jobname(1:lenJ)//'.basis2')
       close(IUnit,STATUS='DELETE')
c
       nwmo = MIN(NAlpha+20,NBas2)     ! # MOs to write to <MOS> file
       Call WriteMOS(NBas1,  nwmo,   CMO1,   jnk,    .False.,
     $               lenM,   MOS,    1)
       If(NBeta.GT.0) Call WriteMOS(NBas1, nwmo, CBeta1, jnk, .False.,
     $                              lenM,   MOB,    1)
C
C  if orbitals were swapped/mixed, remove keyword from <control> file
C
       OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $       FORM='FORMATTED',STATUS='OLD')
       If(swap1.NE.0) Call wrcntrl(IUnit,5,'$swap',-1,idum,rdum,cdum)
       If(swab1.NE.0)
     $   Call wrcntrl(IUnit,10,'$beta_swap',-1,idum,rdum,cdum)
       If(mix1.NE.0) Call wrcntrl(IUnit,4,'$mix',-1,idum,rdum,cdum)
       If(NBeta.GT.0) Call wrcntrl(IUnit,6,'$nbeta',1,NBeta,rdum,cdum)
       CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  make a copy of the <basis> file in case of future GUESS=READ
C
       Call CopyFile(jobname(1:lenJ)//'.basis',lenJ+6,
     $               jobname(1:lenJ)//'.basis2',lenJ+7)
C
C  finally print out symmetries of guess MOs
C  ----------------------------------------------------------------
       If(IPRNT.GE.2.OR.(IPRNT.EQ.1.AND.GUESS(1:4).NE.'read')) Then
         NVirt = 0
         If(IPRNT.GT.1) Then
           NVirt = MIN(10,NBas2-NAlpha)
           If(wvfnc(1:7).EQ.'Semiemp') NVirt = MIN(10,NBas-NAlpha)
         EndIf
         If(GUESS(1:4).EQ.'read') EPrnt=.False.
         Call PrintGuessSym(NBas1,  NAlpha, NBeta,  NVirt,  EPrnt,
     $                      swap1,  swap2,  swab1,  swab2,  EA,
     $                      CMO1,   EB,     CBeta1, S11,    Z)
       EndIf
C  ----------------------------------------------------------------
      EndIf
C
      RETURN
c
 1000 FORMAT(7X,I4,1X,I4)
 1100 FORMAT(12X,I4,1X,I4)
 1200 FORMAT(7X,I4,1X,I4,2X,F8.4)
c
      END
c  =======================================================================
c
      SUBROUTINE SCFGUESS(NAtoms, GUESS,  IPRNT,  AtSymb, IAN,
     $                    XC,     XCharg, NBas1,  NShl1,  NPsp1,
     $                    BASDAT1,INX1,   IPSP1,  NBas2,  NShl2,
     $                    NPsp2,  BASDAT2,INX2,   IPSP2,  S11,
     $                    S12,    NAlpha, NBeta,  CMO1,   CMO2,
     $                    EA,     CBeta1, CBeta2, EB,     Charg,
     $                    NBas,   N2Elec, MDiis,  swap1,  swap2,
     $                    swab1,  swab2,  mix1,   mix2,   angle,
     $                    MSwap,  wvfnc,  NMem,   Z,      SZ,  IErr)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main Driver for SCF guess
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  GUESS   -  type of SCF guess requested
C  IPRNT   -  print flag
C  AtSymb  -  atomic symbols
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C  XCharg  -  atomic charges
C  NBas1   -  number of contracted basis functions in actual basis
C  NShl1   -  number of contracted shells
C  NPsp1   -  number of ECPs in basis
C  BASDAT1 -  basis function data (Texas format)
C  INX1    -    ditto
C  IPSP1   -  number of ECPs on each atomic center
C  NBas2   -  number of contracted basis functions in guess basis
C  NShl2   -  number of contracted shells
C  NPsp2   -  number of ECPs in guess basis
C  BASDAT2 -  basis function data (Texas format)
C  INX2    -    ditto
C  IPSP2   -  number of ECPs on each atomic center in guess basis
C  S11     -  overlap matrix actual basis
C  S12     -  overlap matrix actual/guess basis
C  NAlpha  -  number of occupied closed-shell/alpha MOs
C  NBeta   -  number of occupied beta spin MOs
C  CMO1    -  storage for actual guess alpha MOs
C  CMO2    -  storage for initial guess basis alpha MOs
C  EA      -  storage for guess alpha orbital energies
C  CBeta1  -  storage for actual guess beta MOs
C  CBeta2  -  storage for initial guess basis beta MOs
C  EB      -  storage for guess beta orbital energies
C  Charg   -  total charge on system
C  NBas    -  number of basis functions for semiempirical guess
C  N2Elec  -  number of semiempirical 2-electron repulsion integrals
C  MDiis   -  maximum size of semiempirical DIIS subspace
C  swap1   -  if non-zero swaps occupied closed-shell/alpha
C  swap2   -   MO swap1 with virtual MO swap2
C  swab1   -  if non-zero swaps occupied beta spin
C  swab2   -   MO swab1 with virtual MO swab2
C  mix1    -  if non-zero mixes alpha spin MO mix1 and
C  mix2    -   beta spin MO mix2
C  angle   -  mixing angle
C  MSwap   -  Array for multiple swaps (up to 3)
C  wvfnc   -  wavefunction type (semiempirical or ab initio)
C  NMem    -  scratch memory available
C  Z       -  scratch storage (at least 127008 double words)
C  SZ      -  scratch storage for binary write
C             (may use same storage as Z - see calling routine)
C  IErr    -  Error flag on exit
C               0 - success; MOs generated
C               anything else - failure
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),BASDAT1(13,*),INX1(12,*),
     $          BASDAT2(13,*),INX2(12,*),XCharg(NAtoms),
     $          IPSP1(NAtoms),IPSP2(NAtoms)
      DIMENSION S11(NBas1,NBas1),S12(NBas1,NBas2),CMO1(NBas1,NBas1),
     $          CMO2(NBas2,NBas2),EA(NBas2),CBeta1(NBas1,NBas1),
     $          CBeta2(NBas2,NBas2),EB(NBas2)
      INTEGER swap1,swap2,swab1,swab2
      INTEGER MSwap(2,6)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER GUESS*8,wvfnc*20,SEMI*20
      CHARACTER*256 scrf,jobname
      Logical LHuck
C
      DIMENSION Z(NMem),SZ(NBas1,NBas1)
c ......................................................
c     common /big/ bl(100)     ! from Texas depository
c ......................................................
C
      Common /job/jobname,lenJ
C
      Parameter (PI=3.14159 26535 89793d0)
      Parameter (EVTOAU=1.0d0/27.21139613d0)
C
cc      Data IOut/6/                      ! temporary definition of IO
      IOut=ioutfil('iout')
C
C  Save original NAlpha/NBeta
C
      NAlpha0 = NAlpha
      NBeta0  = NBeta
      LHuck = .False.
c
      If(NAlpha.EQ.NBeta) WRITE(IOut,1000)
C
C ........................................................
C  Sort out the orbital swapping
C  We either have specific user-input orbitals to swap
C  OR swap1/swab1 is -1, in which case we will swap HOMO
C  and LUMO
C
      If(swap1.EQ.-1) Then
       swap1 = NAlpha0
       swap2 = swap1+1
      EndIf
c
      If(swab1.EQ.-1) Then
       swab1 = NBeta0
       swab2 = swab1+1
      EndIf
C
C  Sort out the orbital mixing
C  This is similar to orbital swapping, above
C
      If(mix1.EQ.-1) Then
       mix1 = NAlpha0
       mix2 = NBeta0+1
      EndIf
C
C .......................................................
C  If there are pseudopotentials, determine how many
C  orbitals are affected
C
      If(NPsp1.GT.0) Then
       nsp = 0
       Do I=1,NAtoms
       nsp = nsp + IPSP1(I)
       EndDO
       nsp = nsp/2       ! number of core orbitals
      EndIf
C
C -------------------------------
C  CHECKPOINT GUESS
C -------------------------------
C  Data already read-in for this guess
C
      If(GUESS(1:4).EQ.'read') Then
       WRITE(IOut,1050)
       GO TO 60
      EndIf
C
C ........................................................
C -- remove ghost atoms
      call RemoveGhosts(NAtoms,AtSymb,IAN,XCharg,XC)
C ........................................................
C
C -------------------------------
C  SEMIEMPIRICAL GUESS
C -------------------------------
C
      IF(GUESS(1:3).EQ.'pm3'.OR.GUESS(1:3).EQ.'am1'.OR.
     $   GUESS(1:4).EQ.'mndo'.OR.GUESS(1:5).EQ.'mindo'.OR.
     $   GUESS(1:7).EQ.'default') THEN
c
        SEMI = GUESS
        If(GUESS(1:7).EQ.'default') SEMI = 'pm3'
C
C  Semiempirical guess
C  Can we do it?
C  (not available for all atoms)
C
        CALL ChkSEMI(SEMI,NAtoms,IAN,NAlpha,NBeta,IErr)
c
        If(IErr.NE.0) Then
         WRITE(IOut,1100)
         If(GUESS(1:7).NE.'default') Then
          RETURN
         Else
          NAlpha = NAlpha0       ! restore original value
          NBeta = NBeta0
          GO TO 40
         EndIf
        EndIf
C
C  If pseudopotential on any atom is so deep that it eats into the
C  valence orbitals (i.e., if ECP includes valence S-shell) then
C  do not do semiempirical guess
C
        If(NPsp1.GT.0) Then
         CALL ChkECP(NAtoms,IAN,IPSP1,IErr)
         If(IErr.NE.0) Then
          WRITE(IOut,1150)
          If(GUESS(1:7).NE.'default') Then
           RETURN
          Else
           NAlpha = NAlpha0       ! restore original value
           NBeta = NBeta0
           GO TO 40
          EndIf
         Else
C
C  need to correct NAlpha and NBeta to account for ECP
C
          NAlpha = NAlpha + nsp
          NBeta = NBeta + nsp
          If(NBeta0.EQ.0) NBeta = 0
         EndIf
        EndIf
c
        NLin = (NBas*(NBas+1))/2
        NScr = 2*NLin
C
C  allocate scratch pointers
C
        INF   = 1                       !  NFIRST array
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
        IPA   = IC   + NBas**2          !    ditto  density matrix
        IPAO  = IPA  + NLin             !    ditto  old density matrix
        IBA   = IPAO + NLin             !    ditto  DIIS error matrix
        IEnd  = IBA  + (MDiis+1)**2
        If(NBeta.GT.0) Then
         IFB  = IEnd                    !  beta Fock matrix
         ICB  = IFB  + NLin             !    ditto  MOs
         IPB  = ICB  + NBas**2          !    ditto  density matrix
         IPBO = IPB  + NLin             !    ditto  old density matrix
         IBB  = IPBO + NLin             !    ditto  DIIS error matrix
         IEnd = IBB  + (MDiis+1)**2
        Else
         IFB  = 1
         ICB  = 1
         IPB  = 1
         IPBO = 1
         IBB  = 1
        EndIf
        IScr = IEnd
        IEnd = IScr + NScr - 1
c
        IErr = NMem - IEnd
        If(IErr.LT.0) CALL MemERR(8*NMem,8,'SCFGUESS')
c
        CALL SemiGUESS(NAtoms, SEMI,   AtSymb, IAN,    XC,
     $                 Charg,  NBas,   NAlpha, NBeta,  N2Elec,
     $                 MDiis,  IPRNT,  Z(INF), Z(INM), Z(INL),
     $                 Z(IFAC),Z(IFC1),Z(IUSP),Z(IPSP),Z(IW),
     $                 Z(IH),  Z(IP),  Z(IF),  Z(IC),  EA,
     $                 Z(IPA), Z(IPAO),Z(IBA), Z(IFB), Z(ICB),
     $                 EB,     Z(IPB), Z(IPBO),Z(IBB), NScr,
     $                 Z(IScr),IErr)
c
        IF(IErr.NE.0) THEN
         WRITE(IOut,1200)
         If(GUESS(1:7).NE.'default') Then
          RETURN
         Else
          GO TO 40
         EndIf
        ELSE
C
C  on successful exit from <SemiGUESS>
C  alpha/closed-shell MOs are in Z(IC)   orbital energies in EA
C    beta             MOs are in Z(ICB)  orbital energies in EB
C
         If(IPRNT.GT.3) Then
          WRITE(IOut,1300)
          CALL PrntMAT(NAlpha,NBas,NBas,Z(IC))
c
          If(NBeta.GT.0) Then
            WRITE(IOut,1310)
            CALL PrntMAT(NBeta,NBas,NBas,Z(ICB))
          EndIf
         EndIF
C
C  need to add Core orbitals to valence MOs and prepare basis function data
C
         NBas2 = NBas
         CALL SemiCORE(NAtoms, IAN,    XC,     NBas2,  NShl2,
     $                 BASDAT2,INX2,   IPSP1,  Z(INF))
c
         If(IPRNT.GT.0) Then
c -- need to pad orbital energy arrays, as no orbital energies for Core MOs
           NCore = NBas2-NBas
           DO I=NBas,1,-1
           EA(I+NCore) = EA(I)*EVTOAU
           EndDO
           CALL ZeroIT(EA,NCore)
           If(NBeta.GT.0) Then
            DO I=NBas,1,-1
            EB(I+NCore) = EB(I)*EVTOAU
            EndDO
            CALL ZeroIT(EB,NCore)
           EndIf
         EndIf
c
         WRITE(IOut,1350) SEMI
         GO TO 50        ! success
        ENDIF
c
      ENDIF
c
 40   CONTINUE
C
C -------------------------------
C  HUCKEL GUESS
C -------------------------------
C
      IF(GUESS(1:4).EQ.'huck'.OR.GUESS(1:4).EQ.'defa') THEN
C
C  Extended Huckel guess (includes core orbitals)
C  Can we do it?
C  (currently available for all elements through krypton)
C
        CALL ChkHuckel(NAtoms,IAN,NBas,IErr)
c
        If(IErr.NE.0) Then
         WRITE(IOut,1400)
         RETURN
        EndIf
C
C  need to correct NBas to account for ECP
C
        If(NPsp1.GT.0) NBas = NBas - nsp
C
C  allocate scratch pointers
C  (S11 and S12 used as scratch in <HUCKEL>)
C
        ih = 1
        ich = ih + MAX(NBas*NBas,13*9*NAtoms)
        IEnd = ich + NBas - 1
        CALL MemCHK(NMem,IEnd,8,'SCFGUESS')
c
        CALL HUCKEL(NAtoms, NBas,   NAlpha, NBeta,  IPRNT,
     $              IAN,    XC,     NShl2,  BASDAT2,INX2,
     $              IPSP1,  Z(ih),  S11,    S12,    EA,
     $              Z,      CMO2,   IErr)
c
        IF(IErr.NE.0) THEN
         WRITE(IOut,1500)
         RETURN
        ELSE
C
C  on successful exit from <HUCKEL> MOs in CMO2
C
         If(IPRNT.GT.3) Then
          WRITE(IOut,1600)
          CALL PrntMAT(NAlpha,NBas,NBas,CMO2)
         EndIF
c
         NBas2 = NBas
c
         WRITE(IOut,1650)
         If(NBeta.GT.0) CALL CpyVEC(NBas,EA,EB)
         LHuck = .True.  ! set Huckel flag
         GO TO 60        ! success
        ENDIF
c
      ENDIF
C
C
 50   CONTINUE
C
C  At this point we have semiempirical valence MO guess
C  need to add core orbitals
C  Note: May need MORE than NAlpha MOs if any orbitals are swapped
C
      nmo = NAlpha
      multi_swap = MSwap(1,1) + MSwap(2,1) + MSwap(1,2) + MSwap(2,2)
     $           + MSwap(1,3) + MSwap(2,3) + MSwap(1,4) + MSwap(2,4)
     $           + MSwap(1,5) + MSwap(2,5) + MSwap(1,6) + MSwap(2,6)
      If(MAX(swap2,swab2,mix2,multi_swap).gt.0.OR.IPRNT.GE.2) nmo=NBas
c
      CALL AddCORE(NBas,NBas2,nmo,Z(IC),Z(INF),CMO2)
c
      If(IPRNT.GT.3) Then
       WRITE(IOut,1700)
       CALL PrntMAT(NAlpha0,NBas2,NBas2,CMO2)
      EndIF
c
      IF(NBeta.GT.0) THEN
C
C  add core orbitals for Beta MOs
C
        CALL AddCORE(NBas,NBas2,nmo,Z(ICB),Z(INF),CBeta2)
c
        If(IPRNT.GT.3) Then
         WRITE(IOut,1710)
         CALL PrntMAT(NBeta0,NBas2,NBas2,CBeta2)
        EndIF
      ENDIF
C
C
 60   CONTINUE
C
C  At this point we have initial guess MOs in the "guess basis"
C  (either semiempirical, Huckel or read in)
C  Have to convert these MOs from the guess basis to the actual
C  basis for the current calculation
C
C  normalize the basis sets
C  ** WARNING **  <normaliz> should only be called ONCE
C                 Already called if successful HUCKEL guess
C
      CALL normaliz(NShl1,INX1,BASDAT1)
      If(.NOT.LHuck) CALL normaliz(NShl2,INX2,BASDAT2)
C
C  reallocate scratch storage
C
      NShl = MAX(NShl1,NShl2)
      NBS  = MAX(NBas1,NBas2)
      i1 = 1
      i2 = i1 + 12*NShl
      i3 = i2 + NShl
      IEnd = i3 + NBS - 1
      CALL MemCHK(NMem,IEnd,8,'SCFGUESS')
C
C  reorder the actual basis set (to Wolinski special ordering)
C
      CALL reorder(NShl1,INX1,Z(i3),Z(i2),Z(i1),IAN)
      If(GUESS(1:4).EQ.'read'.AND.wvfnc(1:7).NE.'Semiemp')
     $  CALL reorder(NShl2,INX2,Z(i3),Z(i2),Z(i1),IAN)
C
C  form overlap between actual and guess basis - S12
C
      CALL inton2(0,      NAtoms, S12,    INX1,   INX2,
     $            0,      0,      BASDAT1,BASDAT2,XC,
     $            IAN,    NShl1,  NShl2,  NBas1,  NBas2,
     $            Z(i1))
cc      write(6,*) ' S12 overlap matrix:'
cc      call prntmat(nbas2,nbas1,nbas2,s12)
C
C  form overlap matrix in actual basis - S11
C
      CALL inton2(0,      NAtoms, S11,    INX1,   INX1,
     $            0,      0,      BASDAT1,BASDAT1,XC,
     $            IAN,    NShl1,  NShl1,  NBas1,  NBas1,
     $            Z(i1))
cc      write(6,*) ' INX1 array after reordering is:'
cc      do j=1,NShl1
cc      write(6,*) (inx1(i,j),i=1,12)
cc      enddo
cc      write(6,*) ' S11 overlap matrix:'
cc      call prntmat(nbas1,nbas1,nbas1,s11)
C
C  reallocate scratch storage
C  Note: May need MORE than NOcc MOs if any orbitals are swapped
C
       nmoA = MAX(NAlpha0,swap2,mix2,MSwap(2,1),MSwap(2,2),MSwap(2,3))
       nmoB = MAX(NBeta0,swab2,mix2,MSwap(2,4),MSwap(2,5),MSwap(2,6))
       nmo = MAX(nmoA,nmoB)
       If(nmo.GT.NBas2) Call nerror(5,'GUESS module',
     $     'More MOs requested than basis functions in guess basis',
     $      0,0)
c
       If(IPRNT.GE.2) Then
         nmo=NBas2
         nmoA=NBas2
         nmoB=NBas2
       EndIf
C
      i1 = 1
      i2 = i1 + NBas1*nmo
      i3 = i2 + nmo*nmo
      i4 = i3 + NBas1*NBas1
      IEnd = i4 + nmo - 1
      CALL MemCHK(NMem,IEnd,8,'SCFGUESS')
C
C  now form actual guess MOs
C  see Pulay memorandum "Starting Orbitals for SCF wavefunctions"
C    C1 = S11**-1*S12*C2*[C2(t)*S12(t)*S11**-1*S12*C2]**-1/2
C  (see also  D.Cremer & J.Gauss,  J.Comp.Chem. 7 (1986) 3)
C
C
C 1. ALPHA/CLOSED-SHELL MOS
C =========================
C
C ----------------------------------------------------------------
c -- Be Careful!  Both S11 and S12 are overwritten in <GuessMOS>
c -- S11 may be needed for symmetry printing later; S12 will be
c -- needed if there are beta MOs to transform
      If(NBeta.GT.0.OR.IPRNT.GE.1) Then
        call getchval('scrf',scrf)
        call rmblan(scrf,256,lnscr)
        open(unit=20,file=scrf(1:lnscr)//'_tmp.20',
     $      form='unformatted')
        write(20) S11,S12
      EndIf
c ----------------------------------------------------------------
      CALL GuessMOS(nmoA,   NBas1,  NBas2,  S11,    S12,
     $              CMO2,   Z(i1),  Z(i2),  Z(i3),  S12,
     $              Z(i4),  CMO1,   .True.)
C
C
C 2. BETA MOS
C ===========
C
      IF(NBeta0.GT.0) THEN
C
C  Transform beta MOs
C  If Huckel guess, simply copy alpha MOs into Beta
C
        If(LHuck) Then
         CALL CpyVEC(nmo*NBas1,CMO1,CBeta1)
        Else
c -- restore S11 and S12
         rewind 20
         read (20) S11,S12
         CALL GuessMOS(nmoB,   NBas1,  NBas2,  S11,    S12,
     $                 CBeta2, Z(i1),  Z(i2),  Z(i3),  S12,
     $                 Z(i4),  CBeta1, .False.)
        EndIf
      ENDIF
C --------------------------------------------------------------
C
C  alpha/closed-shell guess MOs are now in CMO1
C  beta guess MOs are in CBeta1
C  are any orbitals to be swapped?
C
      IF(swap1.GT.0.AND.swap2.GT.swap1) THEN
        DO 70 I=1,NBas1
        tmp = CMO1(I,swap1)
        CMO1(I,swap1) = CMO1(I,swap2)
        CMO1(I,swap2) = tmp
 70     CONTINUE
        WRITE(IOut,1800) swap1,swap2
      ENDIF
c
      IF(swab1.GT.0.AND.swab2.GT.swab1) THEN
        DO 80 I=1,NBas1
        tmp = CBeta1(I,swab1)
        CBeta1(I,swab1) = CBeta1(I,swab2)
        CBeta1(I,swab2) = tmp
 80     CONTINUE
        WRITE(IOut,1850) swab1,swab2
      ENDIF
C
C  are any orbitals to be mixed?
C
      IF(mix1.GT.0.AND.mix2.GT.0) THEN
        rangle = angle*PI/180.0d0
        smix1=COS(rangle)
        smix2=SIN(rangle)
        DO 90 I=1,NBas1
        tmp  = smix1*CMO1(I,mix1)+smix2*CBeta1(I,mix2)
        tmp2 = smix1*CBeta1(I,mix1)-smix2*CMO1(I,mix2)
        CMO1(I,mix1)   = tmp
        CBeta1(I,mix1) = tmp2
 90     CONTINUE
        WRITE(IOut,1900) mix1,mix2,angle
      ENDIF
C
C multiple swaps
C
      DO L=1,3
      IF(MSwap(1,L).GT.0.AND.MSwap(2,L).GT.MSwap(1,L)) Then
        nswap1 = MSwap(1,L)
        nswap2 = MSwap(2,L)
        DO 91 I=1,NBas1
        tmp = CMO1(I,nswap1)
        CMO1(I,nswap1) = CMO1(I,nswap2)
        CMO1(I,nswap2) = tmp
 91     CONTINUE
        WRITE(IOut,1800) nswap1,nswap2
      ENDIF
      EndDO
c
      DO L=4,6
      IF(MSwap(1,L).GT.0.AND.MSwap(2,L).GT.MSwap(1,L)) Then
        nswap1 = MSwap(1,L)
        nswap2 = MSwap(2,L)
        DO 92 I=1,NBas1
        tmp = CBeta1(I,nswap1)
        CBeta1(I,nswap1) = CBeta1(I,nswap2)
        CBeta1(I,nswap2) = tmp
 92     CONTINUE
        WRITE(IOut,1850) nswap1,nswap2
      ENDIF
      EndDO
c
      If(IPRNT.GT.3) Then
       WRITE(IOut,1950)
       CALL PrntMAT(NAlpha0,NBas1,NBas1,CMO1)
       If(NBeta0.GT.0) Then
        WRITE(IOut,1960)
        CALL PrntMAT(NBeta0,NBas1,NBas1,CBeta1)
       EndIf
      EndIF
C
C  restore original NAlpha/NBeta
C
      NAlpha = NAlpha0
      NBeta  = NBeta0
c
      If(NBeta.GT.0.OR.IPRNT.GE.1) Then
c -- restore S11 and delete temporary file
        Rewind 20
        Read (20)S11,S12
        Close (Unit=20,STATUS='DELETE')
      EndIf
c
      IErr = 0          ! Success! guess MOs generated
c
 1000 FORMAT(' System is a UHF Singlet')
 1050 FORMAT(' SCF Guess from Previous Calculation')
 1100 FORMAT(' Semiempirical Guess not Possible for this System')
 1150 FORMAT(' ECP too deep for Semiempirical Guess')
 1200 FORMAT(' Semiempirical Guess Failed')
 1300 FORMAT(/,' Converged Semiempirical Molecular Orbitals:')
 1310 FORMAT(/,' Converged Semiempirical Beta Molecular Orbitals:')
 1350 FORMAT(' Semiempirical Guess   Method: ',A5)
 1400 FORMAT(' Huckel Guess not Possible for this System')
 1500 FORMAT(' Huckel Guess Failed')
 1600 FORMAT(/,' Huckel Molecular Orbitals')
 1650 FORMAT(' Extended Huckel Guess')
 1700 FORMAT(/,' Expanded Valence Guess Molecular Orbitals:')
 1710 FORMAT(/,' Expanded Valence Guess Beta Molecular Orbitals:')
 1800 FORMAT(' Alpha/Closed-Shell Orbitals ',I4,' and ',I4,' Swapped')
 1850 FORMAT(' Beta Orbitals ',I4,' and ',I4,' Swapped')
 1900 FORMAT(' Alpha Orbital ',I4,' and Beta Orbital ',I4,' Mixed',/,
     $       '   Rotation Angle ',F9.4,' degrees')
 1950 FORMAT(/,' Final Alpha/Closed-Shell Guess MOs:')
 1960 FORMAT(/,' Final Beta Guess MOs:')
c
      END
c  =======================================================================
c
      SUBROUTINE AddCORE(NBas,NBas2,NOcc,CMO,IEx,CMO2)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  expands semiempirical/Huckel vectors by adding additional zeroes
C  corresponding to added core orbitals
C
C  ARGUMENTS
C
C  NBas    -  number of semiempirical (valence) basis functions
C  NBas2   -  number of basis functions including core
C  NOcc    -  number of occupied valence MOs
C  CMO     -  semiempirical/Huckel valence MOs
C  IEx     -  Expansion vector
C             (contains 1's where additional core orbitals needed)
C  CMO2    -  on exit contains expanded MOs
C
C
      DIMENSION CMO(NBas,NOcc),IEx(NBas2),CMO2(NBas2,*)
c
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
C  first add the extra core MOs
C
      IT = 0
      DO 10 I=1,NBas2
      If(IEx(I).EQ.1) Then
C
C  add unit vector
C
        IT = IT+1
        CALL ZeroIT(CMO2(1,IT),NBas2)
        CMO2(I,IT) = One
      EndIf
 10   CONTINUE
C
C  now expand out the valence MOs
C
      JT = 0
      DO 30 I=1,NBas2
      IF(IEx(I).EQ.1) THEN
        DO 20 J=1,NOcc
        CMO2(I,IT+J) = Zero
 20     CONTINUE
      ELSE
        JT = JT+1
        DO 25 J=1,NOcc
        CMO2(I,IT+J) = CMO(JT,J)
 25     CONTINUE
      ENDIF
 30   CONTINUE
C
cc      write(*,*) 'Original MOs'
cc      do i=1,nocc
cc        write(*,'(i3,5f13.5,/,(3x,5f13.5))') i, (CMO(k,i),k=1,NBas)
cc      end do
cc      write(*,*) 'Expanded MOs'
cc      do i=1,nocc+it
cc        write(*,'(i3,5f13.5,/,(3x,5f13.5))') i, (CMO2(k,i),k=1,NBas2)
cc      end do
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE ChkECP(NAtoms,IAN,IPSP,IErr)
      IMPLICIT INTEGER(A-Z)
C
C  Checks to see if pseudopotential is so deep that it encompasses
C  the valence S orbital - if so we will not do a semiempirical guess
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IPSP    -  number of ECPs per atom
C  IErr    -  Error flag on exit
C               0 - ECP does not encroach upon natural valence space
C              -1 - encroachment for 1 or more atoms
C
C
      DIMENSION IAN(NAtoms),IPSP(NAtoms)
c
      DIMENSION CORE(54)
C
      DATA CORE/0,0, 2,2,2,2,2,2,2,2, 10,10,10,10,10,10,10,10,
     $          18,18,18,18,18,18,18,18,18,18,18,18, 30,30,30,30,30,30,
     $          36,36,36,36,36,36,36,36,36,36,36,36, 48,48,48,48,48,48/
C
C
      IErr = -1
c
      DO 10 I=1,NAtoms
      If(IPSP(I).GT.CORE(IAN(I))) RETURN
 10   CONTINUE
C
      IErr = 0
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE ChkPseud(NAtoms,IPSP1,IPSP2,IErr)
      IMPLICIT INTEGER(A-Z)
C
C  Checks to see if pseudopotentials on guess and current basis match
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IPSP1   -  number of ECPs per atom for current basis
C  IPSP2   -  number of ECPs per atom for guess   basis
C  IErr    -  Error flag on exit
C               0 - exact ECP match
C              -1 - mismatch in ECPs for 1 or more atoms
C
      Dimension IPSP1(NAtoms),IPSP2(NAtoms)
C
      IErr = -1
c
      Do 10 I=1,NAtoms
      If(IPSP1(I).NE.IPSP2(I)) RETURN
 10   CONTINUE
C
      IErr = 0
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE GetMinBas(NAtoms,IAN,NBasis,NShell,NPrim)
      IMPLICIT INTEGER(A-Z)
C
C  sets basis set parameters for a minimal basis
C
C  ARGUMENTS
C
C  NAtoms  -  number of (real) atoms
C  IAN     -  list of atomic numbers
C
C  on exit
C
C  NBasis  -  number of basis functions
C  NShell  -  number of shells
C  NPrim   -  number of primitive shells
C
C
      DIMENSION IAN(NAtoms)
C
C
C  Initialize
C
      NBasis = 0
      NShell = 0
      NPrim = 0
C
C  Loop over atoms
C
      DO 10 I=1,NAtoms
      II = IAN(I)
      If(II.GT.54) Call nerror(13,'GUESS module',
     $      'Sorry! No Guess Available for atoms above Xe',0,0)
c
c 1S Shell
      If(II.LE.2) THEN
        NBasis = NBasis + 1
        NShell = NShell + 1
c 2S Shell
cc      Else If(II.LE.4) Then
cc        NBasis = NBasis + 2
cc        NShell = NShell + 2
c 2P Shell
      Else If(II.LE.10) Then
        NBasis = NBasis + 5
        NShell = NShell + 3
c 3S Shell
cc      Else If(II.LE.12) Then
cc        NBasis = NBasis + 6
cc        NShell = NShell + 4
c 3P Shell
      Else If(II.LE.18) Then
        NBasis = NBasis + 9
        NShell = NShell + 5
c 4S Shell
cc      Else If(II.LE.20) Then
cc        NBasis = NBasis + 10
cc        NShell = NShell + 6
c 3D Shell
      Else If(II.LE.29) Then     ! Zn excluded
        NBasis = NBasis + 15
        NShell = NShell + 7
c 4P Shell
      Else If(II.LE.36) Then
        NBasis = NBasis + 18
        NShell = NShell + 8
c 5S Shell
cc      Else If(II.LE.38) Then
cc        NBasis = NBasis + 19
cc        NShell = NShell + 9
c 4D Shell
      Else If(II.LE.47) Then     ! Cd excluded
        NBasis = NBasis + 24
        NShell = NShell + 10
c 5P Shell
      Else If(II.LE.54) Then
        NBasis = NBasis + 27
        NShell = NShell + 11
      EndIf
 10   CONTINUE
c
      NPrim = 3*NShell     ! assume 3 primitives per shell
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE GuessMOS(NOcc,   NBas1,  NBas2,  S11,    S12,
     $                    CMO2,   SM1,    SM2,    SM3,    SS,
     $                    V,      CMO1,   invert)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms MOs in basis 1 from MOs in basis 2 and overlap matrix
C  see Pulay memorandum "Starting Orbitals for SCF wavefunctions"
C    C1 = S11**-1*S12*C2*[C2(t)*S12(t)*S11**-1*S12*C2]**-1/2
C
C Instead:
c First solve S11 C' = S11*C2 (using fast and well scaling LAPACK)
c Then get C1=C'*[C2(t)*S12(t)*C']**-1/2
c A lot faster than inversion using slow routines!
c
C  ARGUMENTS
C
C  NOcc    -  number of occupied MOs
C  NBas1   -  number of basis functions in basis 1
C  NBas2   -  number of basis functions in basis 2
C  IPRNT   -  print flag
C  S11     -  overlap matrix, basis 1
C  S12     -  overlap matrix, basis 1/basis 2
C  CMO2    -  MO coefficients, basis 2
C  SM1     -  work array
C  SM2     -   ditto
C  SM3     -   ditto
C  SS      -   ditto  (can share storage with S12 - see call)
C  V       -  scratch vector
C  CMO1    -  on exit contains MO coefficients, basis 1
C  invert  -  Logical flag indicating whether on not to invert
C             S11 matrix (if transforming beta MOs, already inverted)
C
C
      DIMENSION S11(NBas1,NBas1),S12(NBas1,NBas2),CMO2(NBas2,NOcc),
     $          SM1(NBas1,NOcc),SM2(NOcc,NOcc),SM3(NBas1,NBas1),
     $          SS(NOcc,NOcc),V(NOcc),CMO1(NBas1,NOcc)
      LOGICAL invert
      integer*4 i4err
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
c ................................................................
c     IF(invert) THEN
C
C  invert S11
C
c       CALL CpyVEC(NBas1*NBas1,S11,SM3)
c       CALL INVMAT(SM3,NBas1,SM1,CMO1,S11,IErr)
c
c       If(IErr.NE.0) Then
c         Call nerror(6,'GUESS module',
c    $       'Unable to Invert Overlap Matrix!',0,0)
c       EndIf
c
c     ENDIF
c ................................................................
C
C  form and save S12*CMO2  (in SM1)
C
      Call DGEMM('N',    'N',    NBas1,  NOcc,   NBas2,
     $            One,    S12,   NBas1,  CMO2,   NBas2,
     $            Zero,   SM1,   NBas1)
c
c form C' (in CMO1)
c
      CALL CpyVEC(NBas1*NOcc,SM1,CMO1)
      call DPOSV('U',NBas1,NOcc,S11,NBas1,CMO1,NBas1,i4err)
      ierr=int(i4err)
      If(IErr.NE.0) Then
         Call nerror(7,'GUESS module',
     $       'Unable to Convert Coeff Matrix!',IErr,0)
      EndIf
C
C  form SM1(t) * C'  (in SM2)
C
      Call DGEMM('T',    'N',    NOcc,   NOcc,   NBas1,
     $            One,    SM1,   NBas1,  CMO1,   NBas1,
     $            Zero,   SM2,   NOcc)

C
C  diagonalize SM2
C
      CALL DIAGMAT(SM2,NOcc,SM3,SS,V,IErr)
c
      If(IErr.NE.0) Then
        Call nerror(8,'GUESS module',
     $     'Unable to Diagonalize Composite Overlap!',0,0)
      EndIf
C
C  now form SM2**-1/2
C
      DO 65 I=1,NOcc
      V(I) = One/SQRT(V(I))
 65   CONTINUE
c
      DO 80 I=1,NOcc
      DO 80 J=1,I
      Val = Zero
      DO 70 K=1,NOcc
      Val = Val + SM2(I,K)*SM2(J,K)*V(K)
 70   CONTINUE
      SS(I,J) = Val
      SS(J,I) = Val
 80   CONTINUE
C
C  get final MOs
C  CMO1 = C' * SS
C
      CALL CpyVEC(NBas1*NOcc,CMO1,SM1)
      Call DGEMM('N',    'N',    NBas1,  NOcc,   NOcc,
     $            One,    SM1,   NBas1,  SS,     NOcc,
     $            Zero,   CMO1,  NBas1)
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE SemiGUESS(NAtoms, SEMI,   AtSymb, IAN,    XC,
     $                     Charg,  NBas,   NAlpha, NBeta,  N2Elec,
     $                     MDiis,  IPRNT,  NFIRST, NMIDLE, NLAST,
     $                     IFACT,  I1FACT, USPD,   PSPD,   W,
     $                     H,      P,      F,      C,      EIGA,
     $                     PA,     POLD,   BA,     FB,     CB,
     $                     EIGB,   PB,     PBOLD,  BB,     NMem,
     $                     Z,      IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Semiempirical Guess
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  SEMI    -  which semiempirical method
C  AtSymb  -  atomic symbols
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C  Charg   -  charge on system
C  NBas    -  number of basis functions (on input only upper bound)
C  NAlpha  -  number of occupied alpha/closed shell orbitals
C  NBeta   -  number of occupied beta orbitals
C  N2Elec  -  number of 2-electron repulsion integrals
C  MDiis   -  maximum size of DIIS subspace
C  IPRNT   -  print flag
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
C  NMem    -  available scratch memory
C  Z       -  scratch storage
C  IErr    -  error flag
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms)
      DIMENSION NFIRST(NAtoms),NMIDLE(NAtoms),NLAST(NAtoms),
     $          IFACT(NBas),I1FACT(NBas)
      DIMENSION USPD(NBas),PSPD(NBas),W(N2Elec),C(NBas,NBas),
     $          H(*),P(*),F(*),PA(*),FB(*),PB(*),POLD(*),PBOLD(*),
     $          EIGA(NBas),EIGB(NBas),CB(NBas,NBas),BA(*),BB(*)
      CHARACTER*8 AtSymb(NAtoms)
      Character*20 SEMI
      DIMENSION Z(NMem)
c
      Character*256 scrf,jobname
      Common /job/jobname,lenJ
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,clight,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
c
      Data JUnit/15/, KUnit/16/         ! unit number for direct access files
C
C
      AUTOAN = 1.0d0/ANTOAU
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
C
C  convert input Cartesians to angstroms
C
      CALL VScal(3*NAtoms,AUTOAN,XC)
C
C  set up semiempirical common blocks
C
      ICharg = NINT(Charg)    ! WARNING - this COULD be a problem!   JB
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
      inp = -1          ! disable read from input file for SCF parameters
c
      CALL ITER(NAtoms, IAN,    XC,     SEMI,   IPRNT,
     $          NFIRST, NMIDLE, NLAST,  IFACT,  I1FACT,
     $          NBas,   PSPD,   NAlpha, NBeta,  MDiis,
     $          H,      W,      P,      F,      C,
     $          EIGA,   PA,     POLD,   BA,     FB,
     $          CB,     EIGB,   PB,     PBOLD,  BB,
     $          ATHEAT, ENuclr, inp,    Z,      EE,  IErr)
C
C  Close DIIS direct access files
C
      CLOSE (UNIT=JUnit,STATUS='DELETE')
      If(NBeta.GT.0) CLOSE (UNIT=KUnit,STATUS='DELETE')
C
C  convert Cartesians back to atomic units
C
      CALL VScal(3*NAtoms,ANTOAU,XC)
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE SemiCORE(NAtoms, IAN,    XC,     NBas,   ncs,
     $                    BASDAT, INX,    IPSP,   IEx)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Determines number of basis functions in expanded semiempirical basis set
C  (including Core MOs) and prepares basis set data
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C             (needed in BASDAT - semiredundant?)
C  NBas    -  on input number of semiempirical basis functions
C             on exit number of expanded semiempirical basis functions
C  ncs     -  on exit number of contracted shells
C  BASDAT  -  on exit contains basis function data (Texas format)
C  INX     -    ditto
C  IPSP    -  number of ECPs per atom
C  IEx     -  on exit contains expansion vector
C             (i.e. where to add extra functions)
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),BASDAT(13,*),INX(12,*),
     $          IPSP(NAtoms),IEx(*)
      Dimension XSTO1S(3,34),XSTO2S(3,26),XSTO2P(3,26),XSTO3S(3,18),
     $          XSTO3P(3,18),XSTO3D(3,16),XSTO4S(3,9),XSTO4P(3,9),
     $          XSTO4D(3,7),dsto(3,9)
c
      COMMON /EXPONT/ ZS(54),ZP(54),ZD(54)       ! semiempirical common
c
c  sto-3g expansion to convert Mindo "slater" MOs to gaussians
c  values from:  Pople et al.  J.Chem.Phys.  51 (1969) 2657
      DATA  e11/0.109818d0/,  e12/0.405771d0/,  e13/2.22766d0/,
     $     d1s1/0.444635d0/, d1s2/0.535328d0/, d1s3/0.154329d0/,
     $      e21/0.0751386d0/, e22/0.231031d0/,  e23/0.994203d0/,
     $     d2s1/0.700115d0/, d2s2/0.399513d0/, d2s3/-0.0999672d0/,
     $     d2p1/0.391957d0/, d2p2/0.607684d0/, d2p3/0.155916d0/
c
c  actual sto-3g exponents for 1S (core) orbital
      DATA XSTO1S / 16.1195755d0,  2.9362006d0,  0.7946505d0,   ! Li
     $              30.1678715d0,  5.4951153d0,  1.4871927d0,   ! Be
     $              48.7911148d0,  8.8873625d0,  2.4052670d0,   ! B
     $              71.6168365d0, 13.0450964d0,  3.5305121d0,   ! C
     $              99.1061710d0, 18.0523129d0,  4.8856602d0,   ! N
     $             130.7093200d0, 23.8088665d0,  6.4436083d0,   ! O
     $             166.6791230d0, 30.3608112d0,  8.2168207d0,   ! F
     $             207.0156100d0, 37.7081510d0, 10.2052970d0,   ! Ne
     $             250.7724300d0, 45.6785110d0, 12.3623880d0,   ! Na
     $             299.2374000d0, 54.5064700d0, 14.7515800d0,   ! Mg
     $             351.4214770d0, 64.0118607d0, 17.3241076d0,   ! Al
     $             407.7975510d0, 74.2808330d0, 20.1032923d0,   ! Si
     $             468.3656380d0, 85.3133856d0, 23.0891316d0,   ! P
     $             533.1257360d0, 97.1095183d0, 26.2816254d0,   ! S
     $             601.3456140d0,109.5358540d0, 29.6446769d0,   ! Cl
     $             674.4465180d0,122.8512750d0, 33.2483494d0,   ! Ar
     $             771.5103680d0,140.5315770d0, 38.0333290d0,   ! K
     $             854.0324950d0,155.5630850d0, 42.1014418d0,   ! Ca
c -------------------- transition metals missing ------------------------
     $             1929.432300d0,351.4485020d0, 95.1156802d0,   ! Zn
     $             2061.424530d0,375.4910520d0, 101.622532d0,   ! Ga
     $             2196.384230d0,400.0741290d0, 108.275673d0,   ! Ge
     $             2337.065670d0,425.6994300d0, 115.210879d0,   ! As
     $             2480.626160d0,451.8490970d0, 122.288592d0,   ! Se
     $             2629.997470d0,479.0573220d0, 129.651607d0,   ! Br
     $             2782.160050d0,506.7739270d0, 137.152802d0,   ! Kr
     $             2938.601530d0,535.2699370d0, 144.864934d0,   ! Rb
     $             3100.983950d0,564.8480980d0, 152.869939d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $             4950.260610d0,901.6960380d0, 244.035319d0,   ! Cd
     $             5158.223360d0,939.5767090d0, 254.287357d0,   ! In
     $             5370.465000d0,978.2367850d0, 264.750333d0,   ! Sn
     $             5586.985540d0,1017.676260d0, 275.424247d0,   ! Sb
     $             5810.060070d0,1058.309560d0, 286.421257d0,   ! Te
     $             6035.182040d0,1099.315810d0, 297.519200d0,   ! I
     $             6264.582900d0,1141.101460d0, 308.828082d0 /  ! Xe
c
c  actual sto-3g exponents for 2S (core) orbital
      DATA XSTO2S / 12.0401930d0,  2.7978819d0,  0.9099580d0,   ! Na
     $              15.1218200d0,  3.5139870d0,  1.1428570d0,   ! Mg
     $              18.8993962d0,  4.3918132d0,  1.4283540d0,   ! Al
     $              23.1936561d0,  5.3897069d0,  1.75289995d0,  ! Si
     $              28.0326396d0,  6.5141826d0,  2.11861435d0,  ! P
     $              33.3297517d0,  7.7451175d0,  2.5189526d0,   ! S
     $              38.9604189d0,  9.0535635d0,  2.9444998d0,   ! Cl
     $              45.1642439d0, 10.4951990d0,  3.41336445d0,  ! Ar
     $              52.4020398d0, 12.1771071d0,  3.96037317d0,  ! K
     $              59.5602994d0, 13.8405327d0,  4.50137080d0,  ! Ca
c -------------------- transition metals missing ------------------------
     $              155.841676d0, 36.2142539d0, 11.7779993d0,   ! Zn
     $              167.761868d0, 38.9842503d0, 12.6788881d0,   ! Ga
     $              180.389038d0, 41.9185330d0, 13.6332079d0,   ! Ge
     $              193.197054d0, 44.8948404d0, 14.6011955d0,   ! As
     $              206.157934d0, 47.9065882d0, 15.5807401d0,   ! Se
     $              219.835026d0, 51.0849322d0, 16.6144055d0,   ! Br
     $              233.951412d0, 54.3652768d0, 17.6812753d0,   ! Kr
     $              248.507037d0, 57.7476910d0, 18.7813410d0,   ! Rb
     $              263.501901d0, 61.2321749d0, 19.9146037d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              433.447056d0, 100.723602d0, 32.7585061d0,   ! Cd
     $              452.331445d0, 105.111920d0, 34.1857262d0,   ! In
     $              472.051661d0, 109.694466d0, 35.6761153d0,   ! Sn
     $              492.192623d0, 114.374784d0, 37.1983032d0,   ! Sb
     $              512.754331d0, 119.152875d0, 38.7522896d0,   ! Te
     $              533.736787d0, 124.028738d0, 40.3380748d0,   ! I
     $              555.139989d0, 129.002374d0, 41.9556585d0 /  ! Xe
c
c  actual sto-3g exponents for 2P (core) orbital
      DATA XSTO2P / 12.0401930d0,  2.7978819d0,  0.9099580d0,   ! Na
     $              15.1218200d0,  3.5139870d0,  1.1428570d0,   ! Mg
     $              18.8993962d0,  4.3918132d0,  1.4283540d0,   ! Al
     $              23.1936561d0,  5.3897069d0,  1.75289995d0,  ! Si
     $              28.0326396d0,  6.5141826d0,  2.11861435d0,  ! P
     $              33.3297517d0,  7.7451175d0,  2.5189526d0,   ! S
     $              38.9604189d0,  9.0535635d0,  2.9444998d0,   ! Cl
     $              45.1642439d0, 10.4951990d0,  3.41336445d0,  ! Ar
     $              52.4020398d0, 12.1771071d0,  3.96037317d0,  ! K
     $              59.5602994d0, 13.8405327d0,  4.50137080d0,  ! Ca
c -------------------- transition metals missing ------------------------
     $              155.841676d0, 36.2142539d0, 11.7779993d0,   ! Zn
     $              167.761868d0, 38.9842503d0, 12.6788881d0,   ! Ga
     $              180.389038d0, 41.9185330d0, 13.6332079d0,   ! Ge
     $              193.197054d0, 44.8948404d0, 14.6011955d0,   ! As
     $              206.157934d0, 47.9065882d0, 15.5807401d0,   ! Se
     $              219.835026d0, 51.0849322d0, 16.6144055d0,   ! Br
     $              233.951412d0, 54.3652768d0, 17.6812753d0,   ! Kr
     $              248.507037d0, 57.7476910d0, 18.7813410d0,   ! Rb
     $              263.501901d0, 61.2321749d0, 19.9146037d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              433.447056d0, 100.723602d0, 32.7585061d0,   ! Cd
     $              452.331445d0, 105.111920d0, 34.1857262d0,   ! In
     $              472.051661d0, 109.694466d0, 35.6761153d0,   ! Sn
     $              492.192623d0, 114.374784d0, 37.1983032d0,   ! Sb
     $              512.754331d0, 119.152875d0, 38.7522896d0,   ! Te
     $              533.736787d0, 124.028738d0, 40.3380748d0,   ! I
     $              555.139989d0, 129.002374d0, 41.9556585d0 /  ! Xe
c
c  actual sto-3g exponents for 3S (core) orbital
      DATA XSTO3S / 3.65158398d0, 1.01878266d0, 0.39874463d0,   ! K
     $              4.37470626d0, 1.22053194d0, 0.47770793d0,   ! Ca
c -------------------- transition metals missing ------------------------
     $              12.2815274d0, 3.74625733d0, 1.44542254d0,   ! Zn
     $              12.6150552d0, 3.84799393d0, 1.48467568d0,   ! Ga
     $              14.1966562d0, 4.33043264d0, 1.67081554d0,   ! Ge
     $              15.8716358d0, 4.84135482d0, 1.86794520d0,   ! As
     $              17.6399760d0, 5.38074398d0, 2.07606597d0,   ! Se
     $              19.5017311d0, 5.94864958d0, 2.29517394d0,   ! Br
     $              21.4568467d0, 6.54502216d0, 2.52527302d0,   ! Kr
     $              23.5053410d0, 7.16987820d0, 2.76636191d0,   ! Rb
     $              25.5788669d0, 7.80236971d0, 3.01039679d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              52.5927382d0, 16.0424289d0, 6.18969063d0,   ! Cd
     $              55.9753401d0, 17.0742281d0, 6.58779234d0,   ! In
     $              59.1513510d0, 18.0430107d0, 6.96158016d0,   ! Sn
     $              62.5217334d0, 19.0710827d0, 7.35824375d0,   ! Sb
     $              65.9854944d0, 20.1276381d0, 7.76589716d0,   ! Te
     $              69.5426339d0, 21.2126768d0, 8.18454038d0,   ! I
     $              73.0776598d0, 22.2909702d0, 8.60058103d0 /  ! Xe
c
c  actual sto-3g exponents for 3P (core) orbital
      DATA XSTO3P / 3.65158398d0, 1.01878266d0, 0.39874463d0,   ! K
     $              4.37470626d0, 1.22053194d0, 0.47770793d0,   ! Ca
c -------------------- transition metals missing ------------------------
     $              12.2815274d0, 3.74625733d0, 1.44542254d0,   ! Zn
     $              12.6150552d0, 3.84799393d0, 1.48467568d0,   ! Ga
     $              14.1966562d0, 4.33043264d0, 1.67081554d0,   ! Ge
     $              15.8716358d0, 4.84135482d0, 1.86794520d0,   ! As
     $              17.6399760d0, 5.38074398d0, 2.07606597d0,   ! Se
     $              19.5017311d0, 5.94864958d0, 2.29517394d0,   ! Br
     $              21.4568467d0, 6.54502216d0, 2.52527302d0,   ! Kr
     $              23.5053410d0, 7.16987820d0, 2.76636191d0,   ! Rb
     $              25.5788669d0, 7.80236971d0, 3.01039679d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              52.5927382d0, 16.0424289d0, 6.18969063d0,   ! Cd
     $              55.9753401d0, 17.0742281d0, 6.58779234d0,   ! In
     $              59.1513510d0, 18.0430107d0, 6.96158016d0,   ! Sn
     $              62.5217334d0, 19.0710827d0, 7.35824375d0,   ! Sb
     $              65.9854944d0, 20.1276381d0, 7.76589716d0,   ! Te
     $              69.5426339d0, 21.2126768d0, 8.18454038d0,   ! I
     $              73.0776598d0, 22.2909702d0, 8.60058103d0 /  ! Xe
c
c  actual sto-3g exponents for 3D (core) orbital
      DATA XSTO3D / 10.9473708d0, 3.33929702d0, 1.28840460d0,   ! Zn
     $              12.6150552d0, 3.84799393d0, 1.48467568d0,   ! Ga
     $              14.1966562d0, 4.33043264d0, 1.67081554d0,   ! Ge
     $              15.8716358d0, 4.84135482d0, 1.86794520d0,   ! As
     $              17.6399760d0, 5.38074398d0, 2.07606597d0,   ! Se
     $              19.5017311d0, 5.94864958d0, 2.29517394d0,   ! Br
     $              21.4568467d0, 6.54502216d0, 2.52527302d0,   ! Kr
     $              23.5053410d0, 7.16987820d0, 2.76636191d0,   ! Rb
     $              25.5788669d0, 7.80236971d0, 3.01039679d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              52.5927382d0, 16.0424289d0, 6.18969063d0,   ! Cd
     $              55.9753401d0, 17.0742281d0, 6.58779234d0,   ! In
     $              59.1513510d0, 18.0430107d0, 6.96158016d0,   ! Sn
     $              62.5217334d0, 19.0710827d0, 7.35824375d0,   ! Sb
     $              65.9854944d0, 20.1276381d0, 7.76589716d0,   ! Te
     $              69.5426339d0, 21.2126768d0, 8.18454038d0,   ! I
     $              73.0776598d0, 22.2909702d0, 8.60058103d0 /  ! Xe
c
c  actual sto-3g exponents for 4S (core) orbital
      DATA XSTO4S / 2.24779682d0, 0.82957839d0, 0.36635057d0,   ! Rb
     $              2.46103240d0, 0.90827573d0, 0.40110414d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              5.67485388d0, 2.20975776d0, 0.97274043d0,   ! Cd
     $              5.04855104d0, 1.96587878d0, 0.865384344d0,  ! In
     $              5.58314058d0, 2.17404509d0, 0.957019631d0,  ! Sn
     $              6.12069540d0, 2.38336606d0, 1.04916320d0,   ! Sb
     $              6.70795939d0, 2.61204352d0, 1.14982754d0,   ! Te
     $              7.29599388d0, 2.84102101d0, 1.25062396d0,   ! I
     $              7.90873119d0, 3.07961764d0, 1.35565474d0 /  ! Xe
c
c  actual sto-3g exponents for 4P (core) orbital
      DATA XSTO4P / 2.24779682d0, 0.82957839d0, 0.36635057d0,   ! Rb
     $              2.46103240d0, 0.90827573d0, 0.40110414d0,   ! Sr
c -------------------- transition metals missing ------------------------
     $              5.67485388d0, 2.20975776d0, 0.97274043d0,   ! Cd
     $              5.04855104d0, 1.96587878d0, 0.865384344d0,  ! In
     $              5.58314058d0, 2.17404509d0, 0.957019631d0,  ! Sn
     $              6.12069540d0, 2.38336606d0, 1.04916320d0,   ! Sb
     $              6.70795939d0, 2.61204352d0, 1.14982754d0,   ! Te
     $              7.29599388d0, 2.84102101d0, 1.25062396d0,   ! I
     $              7.90873119d0, 3.07961764d0, 1.35565474d0 /  ! Xe
c
c  actual sto-3g exponents for 4D (core) orbital
      DATA XSTO4D / 5.67485388d0, 2.20975776d0, 0.97274043d0,   ! Cd
     $              5.04855104d0, 1.96587878d0, 0.865384344d0,  ! In
     $              5.58314058d0, 2.17404509d0, 0.957019631d0,  ! Sn
     $              6.12069540d0, 2.38336606d0, 1.04916320d0,   ! Sb
     $              6.70795939d0, 2.61204352d0, 1.14982754d0,   ! Te
     $              7.29599388d0, 2.84102101d0, 1.25062396d0,   ! I
     $              7.90873119d0, 3.07961764d0, 1.35565474d0 /  ! Xe
c
      DATA dsto / 0.15432897d0,  0.53532815d0,  0.44463453d0,   ! 1S
     $           -0.09996723d0,  0.39951283d0,  0.70011547d0,   ! 2S
     $            0.15591627d0,  0.60768372d0,  0.39195739d0,   ! 2P
     $           -0.21962037d0,  0.22559543d0,  0.90039843d0,   ! 3S
     $            0.01058760d0,  0.59516701d0,  0.46200101d0,   ! 3P
     $            0.21976795d0,  0.65554736d0,  0.28657326d0,   ! 3D
     $           -0.30884412d0,  0.01960641d0,  1.13103444d0,   ! 4S
     $           -0.12154686d0,  0.57152276d0,  0.54989495d0,   ! 4P
     $            0.12506621d0,  0.66867856d0,  0.30524682d0 /  ! 4D
C
C
      NShl = 0
      ncs = 0
      nst = 0
      ncf = 0
c
      DO 10 I=1,NAtoms
      II = IAN(I)
      X = XC(1,I)
      Y = XC(2,I)
      Z = XC(3,I)
      JPSP = IPSP(I)     ! no. of core electrons in ECP
c
      IF(II.EQ.1) THEN
        NShl = NShl+1
        BASDAT(1,NShl) = ZS(1)*e11
        BASDAT(2,NShl) = d1s1
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = ZS(1)*e12
        BASDAT(2,NShl) = d1s2
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = ZS(1)*e13
        BASDAT(2,NShl) = d1s3
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
c
        ncs = ncs+1
        INX(1,ncs) = nst
        INX(2,ncs) = I
        INX(3,ncs) = 1
        nst = nst+3
        INX(5,ncs) = nst
        INX(11,ncs) = ncf
        ncf = ncf+1
        INX(10,ncs) = ncf
        INX(12,ncs) = 1
        IEx(ncf) = 0
cc
      ELSE
C
C  use STO-3G basis for this atom for the core
C
C  1.  1S orbital
C  --------------
C
        If(II.LE.20) Then
          I2 = II-2
        Else If(II.LE.38) Then
          I2 = II-11
        Else
          I2 = II-20
        EndIf
c
        If(JPSP.EQ.0) Then
        NShl = NShl+1
        BASDAT(1,NShl) = XSTO1S(1,I2)
        BASDAT(2,NShl) = dsto(1,1)
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = XSTO1S(2,I2)
        BASDAT(2,NShl) = dsto(2,1)
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = XSTO1S(3,I2)
        BASDAT(2,NShl) = dsto(3,1)
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
c
        ncs = ncs+1
        INX(1,ncs) = nst
        INX(2,ncs) = I
        INX(3,ncs) = 1
        nst = nst+3
        INX(5,ncs) = nst
        INX(11,ncs) = ncf
        ncf = ncf+1
        INX(10,ncs) = ncf
        INX(12,ncs) = 1
        IEx(ncf) = 1
c
        NBas = NBas+1
        EndIf
cc
        IF(II.GT.10) THEN
C
C  2.  2S orbital
C  --------------
C
          If(II.LE.20) Then
            I2 = II-10
          Else If(II.LE.38) Then
            I2 = II-19
          Else
            I2 = II-28
          EndIf
c
          If(JPSP.LE.2) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO2S(1,I2)
          BASDAT(2,NShl) = dsto(1,2)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO2S(2,I2)
          BASDAT(2,NShl) = dsto(2,2)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO2S(3,I2)
          BASDAT(2,NShl) = dsto(3,2)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 1
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+1
          INX(10,ncs) = ncf
          INX(12,ncs) = 1
          IEx(ncf) = 1
c
          NBas = NBas+1
          EndIf
C
C  3.  2P orbital
C  --------------
C
          If(JPSP.LE.4) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO2P(1,I2)
          BASDAT(2,NShl) = dsto(1,3)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO2P(2,I2)
          BASDAT(2,NShl) = dsto(2,3)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO2P(3,I2)
          BASDAT(2,NShl) = dsto(3,3)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 3
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+3
          INX(10,ncs) = ncf
          INX(12,ncs) = 2
          IEx(ncf) = 1
          IEx(ncf-1) = 1
          IEx(ncf-2) = 1
c
          NBas = NBas+3
          EndIf
cc
        ENDIF
cc
        IF(II.GT.18) THEN
C
C  4.  3S orbital
C  --------------
C
          If(II.LE.20) Then
            I2 = II-18
          Else If(II.LE.38) Then
            I2 = II-27
          Else
            I2 = II-36
          EndIf
c
          If(JPSP.LE.10) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3S(1,I2)
          BASDAT(2,NShl) = dsto(1,4)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3S(2,I2)
          BASDAT(2,NShl) = dsto(2,4)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3S(3,I2)
          BASDAT(2,NShl) = dsto(3,4)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 1
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+1
          INX(10,ncs) = ncf
          INX(12,ncs) = 1
          IEx(ncf) = 1
c
          NBas = NBas+1
          EndIf
C
C  5.  3P orbital
C  --------------
C
          If(JPSP.LE.12) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3P(1,I2)
          BASDAT(2,NShl) = dsto(1,5)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3P(2,I2)
          BASDAT(2,NShl) = dsto(2,5)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3P(3,I2)
          BASDAT(2,NShl) = dsto(3,5)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 3
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+3
          INX(10,ncs) = ncf
          INX(12,ncs) = 2
          IEx(ncf) = 1
          IEx(ncf-1) = 1
          IEx(ncf-2) = 1
c
          NBas = NBas+3
          EndIf
cc
        ENDIF
cc
        IF(II.GT.29) THEN
C
C  6.  3D orbital
C  --------------
C
          If(II.LE.38) Then
            I2 = II-29
          Else
            I2 = II-38
          EndIf
c
          If(JPSP.LE.18) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3D(1,I2)
          BASDAT(2,NShl) = dsto(1,6)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3D(2,I2)
          BASDAT(2,NShl) = dsto(2,6)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO3D(3,I2)
          BASDAT(2,NShl) = dsto(3,6)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 5
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+5
          INX(10,ncs) = ncf
          INX(12,ncs) = 4
          IEx(ncf) = 1
          IEx(ncf-1) = 1
          IEx(ncf-2) = 1
          IEx(ncf-3) = 1
          IEx(ncf-4) = 1
c
          NBas = NBas+5
          EndIf
cc
        ENDIF
cc
        IF(II.GT.36) THEN
C
C  4.  4S orbital
C  --------------
C
          If(II.LE.38) Then
            I2 = II-36
          Else
            I2 = II-45
          EndIf
c
          If(JPSP.LE.28) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4S(1,I2)
          BASDAT(2,NShl) = dsto(1,7)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4S(2,I2)
          BASDAT(2,NShl) = dsto(2,7)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4S(3,I2)
          BASDAT(2,NShl) = dsto(3,7)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 1
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+1
          INX(10,ncs) = ncf
          INX(12,ncs) = 1
          IEx(ncf) = 1
c
          NBas = NBas+1
          EndIf
C
C  5.  4P orbital
C  --------------
C
          If(JPSP.LE.30) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4P(1,I2)
          BASDAT(2,NShl) = dsto(1,8)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4P(2,I2)
          BASDAT(2,NShl) = dsto(2,8)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4P(3,I2)
          BASDAT(2,NShl) = dsto(3,8)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 3
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+3
          INX(10,ncs) = ncf
          INX(12,ncs) = 2
          IEx(ncf) = 1
          IEx(ncf-1) = 1
          IEx(ncf-2) = 1
c
          NBas = NBas+3
          EndIf
cc
        ENDIF
cc
        IF(II.GT.47) THEN
C
C  6.  4D orbital
C  --------------
C
          I2 = II-47
c
          If(JPSP.LE.36) Then
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4D(1,I2)
          BASDAT(2,NShl) = dsto(1,9)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4D(2,I2)
          BASDAT(2,NShl) = dsto(2,9)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
          NShl = NShl+1
          BASDAT(1,NShl) = XSTO4D(3,I2)
          BASDAT(2,NShl) = dsto(3,9)
          BASDAT(11,NShl) = X
          BASDAT(12,NShl) = Y
          BASDAT(13,NShl) = Z
c
          ncs = ncs+1
          INX(1,ncs) = nst
          INX(2,ncs) = I
          INX(3,ncs) = 5
          nst = nst+3
          INX(5,ncs) = nst
          INX(11,ncs) = ncf
          ncf = ncf+5
          INX(10,ncs) = ncf
          INX(12,ncs) = 4
          IEx(ncf) = 1
          IEx(ncf-1) = 1
          IEx(ncf-2) = 1
          IEx(ncf-3) = 1
          IEx(ncf-4) = 1
c
          NBas = NBas+5
          EndIf
cc
        ENDIF
C
C  now add valence
C
        NShl = NShl+1
        BASDAT(1,NShl) = ZS(II)*e21
        BASDAT(2,NShl) = d2s1
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = ZS(II)*e22
        BASDAT(2,NShl) = d2s2
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = ZS(II)*e23
        BASDAT(2,NShl) = d2s3
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
c
        ncs = ncs+1
        INX(1,ncs) = nst
        INX(2,ncs) = I
        INX(3,ncs) = 1
        nst = nst+3
        INX(5,ncs) = nst
        INX(11,ncs) = ncf
        ncf = ncf+1
        INX(10,ncs) = ncf
        INX(12,ncs) = 1
        IEx(ncf) = 0
c
        NShl = NShl+1
        BASDAT(1,NShl) = ZP(II)*e21
        BASDAT(2,NShl) = d2p1
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = ZP(II)*e22
        BASDAT(2,NShl) = d2p2
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
        NShl = NShl+1
        BASDAT(1,NShl) = ZP(II)*e23
        BASDAT(2,NShl) = d2p3
        BASDAT(11,NShl) = X
        BASDAT(12,NShl) = Y
        BASDAT(13,NShl) = Z
c
        ncs = ncs+1
        INX(1,ncs) = nst
        INX(2,ncs) = I
        INX(3,ncs) = 3
        nst = nst+3
        INX(5,ncs) = nst
        INX(11,ncs) = ncf
        ncf = ncf+3
        INX(10,ncs) = ncf
        INX(12,ncs) = 2
        IEx(ncf) = 0
        IEx(ncf-1) = 0
        IEx(ncf-2) = 0
c
      ENDIF
 10   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE PrintGuessSym(NBasis, NAlpha, NBeta,  NVirt,  EPrnt,
     $                         swap1,  swap2,  swab1,  swab2,  EA,
     $                         CMO,    EB,     CBeta,  S,      SS)
      use memory
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Prints out symmetry species of initial guess MOs by interfacing
C  with Texas routine <PrintOrbSym>
C
C  ARGUMENTS
C
C  NBasis  -  number of basis functions
C  NAlpha  -  number of alpha/closed-shell occupied orbitals
C  NBeta   -  number of beta occupied orbitals
C  NVirt   -  maximum number of virtual MOs to be printed
C  EPrnt   -  logical flag for printing orbital energies
C  swap1   -  if non-zero swaps occupied closed-shell/alpha
C  swap2   -   MO swap1 with virtual MO swap2
C  swab1   -  if non-zero swaps occupied beta spin
C  swab2   -   MO swab1 with virtual MO swab2
C  EA      -  alpha semiempirical/Huckel guess orbital energies
C  CMO     -  alpha/closed-shell MO coefficients
C  EB      -  beta semiempirical/Huckel guess orbital energies
C  CBeta   -  beta MO coefficients
C  S       -  overlap matrix
C  SS      -  storage for overlap as lower triangle
C
C
      DIMENSION CMO(NBasis,*),CBeta(NBasis,*),EA(*),EB(*),
     $          S(NBasis,NBasis),SS(*)
      INTEGER swap1,swap2,swab1,swab2
      Logical EPrnt
c     common /intbl/maxsh,inx(100)
c
c
c -- get data from old Texas depository
      call getival('nsym',nsym)
      call getival('SymFunPr',ifp)
      call getival('iout',IOut)
c
c -- convert overlap matrix to lower triangle
      ij = 0
      do i=1,nbasis
      do j=1,i
      ij = ij+1
      ss(ij) = s(j,i)
      enddo
      enddo
c
c -- now call orbital symmetry printer
      If(nsym.gt.0) Then
        If(NBeta.GT.0) WRITE(IOut,1000)
        If(swap1.GT.0) Then
         ETemp = EA(swap1)
         EA(swap1) = EA(swap2)
         EA(swap2) = ETemp
        EndIf
        Call PrintOrbSym(nsym,   NBasis, NAlpha, NVirt,  CMO,
     $                   SS,   bl(ifp), EPrnt,  EA)
        If(swap1.GT.0) Then
         ETemp = EA(swap1)
         EA(swap1) = EA(swap2)
         EA(swap2) = ETemp
        EndIf
        If(NBeta.GT.0) Then
          WRITE(IOut,1100)
          If(swab1.GT.0) Then
           ETemp = EB(swab1)
           EB(swab1) = EB(swab2)
           EB(swab2) = ETemp
          EndIf
          Call PrintOrbSym(nsym,   NBasis, NBeta, NVirt,  CBeta,
     $                     SS,   bl(ifp), EPrnt, EB)
          If(swab1.GT.0) Then
           ETemp = EB(swab1)
           EB(swab1) = EB(swab2)
           EB(swab2) = ETemp
          EndIf
        EndIf
      EndIf
c
      RETURN
c
 1000 FORMAT(/,'  ALPHA SPIN')
 1100 FORMAT(/,'  BETA SPIN')
c
      END
c  =======================================================================
c
      SUBROUTINE SemiEXPONT(SEMI)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  sets up basis exponents for semiempirical orbitals
C
      COMMON /MNDO / USSM(54),UPPM(54),UDDM(54),ZSM(54),ZPM(54),
     1               ZDM(54),BETASM(54),BETAPM(54),BETADM(54),
     2               ALPM(54),EISOLM(54),DDM(54),QQM(54),AMM(54),
     3               ADM(54),AQM(54),GSSM(54),GSPM(54),GPPM(54),
     4               GP2M(54),HSPM(54),POLVOM(54)
      COMMON /PM3 /  USSPM3(54),UPPPM3(54),UDDPM3(54),ZSPM3(54),
     $      ZPPM3(54),ZDPM3(54),BETASP(54),BETAPP(54),BETADP(54),
     $      ALPPM3(54),EISOLP(54),DDPM3(54),QQPM3(54),AMPM3(54),
     $      ADPM3(54),AQPM3(54),GSSPM3(54),GSPPM3(54),GPPPM3(54),
     $      GP2PM3(54),HSPPM3(54),POLVOP(54)
      COMMON /EXPON3 /  ZS3(18),ZP3(18)
      COMMON /EXPONT/ ZS(54),ZP(54),ZD(54)
c
      Character*20 SEMI
C
C  determine what parameters to use
C  (Default set up appears to be AM1)
C
      IF(SEMI(1:3).EQ.'pm3') THEN
        DO 10 I=1,54
        ZS(I)=ZSPM3(I)
        ZP(I)=ZPPM3(I)
        ZD(I)=ZDPM3(I)
 10     CONTINUE
      ELSE IF(SEMI(1:4).EQ.'mndo') THEN
        DO 20 I=1,54
        ZS(I)=ZSM(I)
        ZP(I)=ZPM(I)
        ZD(I)=ZDM(I)
 20     CONTINUE
      ELSE IF(SEMI(1:5).EQ.'mindo') THEN
        DO 30 I=1,17
        ZS(I)=ZS3(I)
        ZP(I)=ZP3(I)
 30     CONTINUE
      ENDIF
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE SortBAS1(NAtoms,NShell,ILST,INX)
      IMPLICIT INTEGER(A-Z)
C
C  Forms INX array from precursor ILST array
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NShell  -  total number of contracted shells
C  ILST    -  array containing data pertaining to each shell
C              (1,J) - shell type of Jth shell
C              (2,J) - # primitive shells in Jth shell
C              (3,J) - atom associated with Jth shell
C              (4,J) - # general contractions
C  INX     -  on exit contains shell data sorted for PQS integral routines
C             (pointers refer to locations in BASDAT array)
C              (1,J) - pointer to start of primitive contraction
C              (2,J) - pointer to end of primitive contraction
C              (3,J) - shell size (see ISIZE array, below)
C              (4,J) - number of general contraction
C                      (0 in segmented scheme, i.e., "standard" basis sets)
C              (5,J) - last primitive in contraction
C          (6,J)-(8,J) - used in general contraction scheme
C              (9,J) - unused
C             (10,J) - ending index of contracted shell
C             (11,J) - starting index of contracted shell
C             (12,J) - function type (1=S,2=P,3=L etc..)
C
C  ----------------------------------------------------------------------
C   The input ILST array should have the basis functions ordered
C   per shell.  SEE <RdBASIS>.
C  ----------------------------------------------------------------------
C
C
      DIMENSION ILST(4,NShell),INX(12,NShell)
      Dimension ISIZE(13)
C
c                   S    P    L    D   D6    F  F10  G15  H21  I28
      Data ISIZE /  1,   3,   4,   5,   6,   7,  10,  15,  21,  28 ,
C                   G    H    I
     $              9,  11,  13 /
C
C
      CALL IZeroIT(INX,12*NShell)
C
C  first get running total of primitives in ILST(2,J)
C
      Do J=2,NShell
      ILST(2,J) = ILST(2,J) + ILST(2,J-1)
      EndDo
C
C  put shells as they come directly into INX array
C
      ncs = 0                ! shell counter
c
      DO 20 J=1,NShell
      IAtom = ILST(3,J)
      IType = ILST(1,J)
      ngr = ILST(4,J)
        ncs = ncs+1
        If(ncs.EQ.1) Then
         INX(1,ncs) = 0
         INX(10,ncs) = (1+ngr)*ISIZE(IType)
         INX(11,ncs) = 0
        Else
         INX(1,ncs) =  ILST(2,J-1)
         INX(10,ncs) = INX(10,ncs-1) + ISIZE(IType)*(1+ngr)
         INX(11,ncs) = INX(10,ncs-1)
        EndIf
        INX(2,ncs) = IAtom
        INX(3,ncs) = ISIZE(IType)
        INX(4,ncs) = ngr
        INX(5,ncs) = ILST(2,J)
        INX(12,ncs) = IType
 20   CONTINUE
C
C  If the INX array produced by this routine is input to <reorder>
C  then it should produce the special Wolinski ordering for efficient
C  integral evaluation
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE SortBAS2(ncs,INX0,iswap,iorder,INX)
      IMPLICIT INTEGER(A-Z)
C
C  Sort basis according to atom order and provide indexing array
C  relating basis function order in orginal and new order
C
C  ARGUMENTS
C
C  ncs     -  total number of contracted shells
C  INX0    -  work array for copy of INX array
C  iswap   -  work array for intermediate shell reordering
C  iorder  -  relates original basis function order to new order
C  INX     -  on entry shell data sorted in original Wolinski order
C             on exit contains shell data sorted per atom order
C             (pointers refer to locations in BASDAT array)
C             SEE <SortBAS1>
C
      DIMENSION INX0(12,ncs),iswap(ncs),iorder(*),INX(12,ncs)
C
C  ** NOTE: the main reorder of INX was previously added by PP to
C           <reorder> now included here in separate subroutine
C
C
C  first make copy of original INX array and initialize iswap
C
      Call ICpyVEC(12*ncs,INX,INX0)
      Do i=1,ncs
      iswap(i) = i
      EndDo
C
C  now reorder INX per atom
C
 50   continue
      iexch=0
      do i=1,ncs-1
      ics0=i
      ics1=i+1
      ina0=inx(2,ics0)
      ina1=inx(2,ics1)
      if(ina0.gt.ina1) then
        iexch=1
        do k=1,12
        itemp=inx(k,ics0)
        inx(k,ics0)=inx(k,ics1)
        inx(k,ics1)=itemp
        end do
        itemp=iswap(ics0)
        iswap(ics0)=iswap(ics1)
        iswap(ics1)=itemp
      end if
      end do
      icf=0
      do i=1,ncs
      inx(11,i)=icf
      inx(10,i)=icf+inx(3,i)*(1+inx(4,i))
      icf=inx(10,i)
      end do
      if(iexch.eq.1) go to 50
C
C  At this point INX array is sorted per ATOM
C  e.g.   same order as original coordinates
C
c  now create the array iorder
      do 800 ics=1,ncs
        ibeg=inx(11,ics)+1
        iend=inx(10,ics)
c  these are the indices of the functions belonging to the present
c  contraction. this contraction used to be iswap(ics) in the
c  original order.
        icso=iswap(ics)
        ibego=inx0(11,icso)+1
        iendo=inx0(10,icso)
        idiff=(iend-ibeg)-(iendo-ibego)
        if (idiff.ne.0) then
          Call nerror(1,'Call from <SortBAS2>',
     $          'Shell Size Wrong When Reordering Per Atom',0,0)
        end if
        ii=0
        do 850 i=ibego,iendo
          iorder(i)=ibeg+ii
          ii=ii+1
 850    continue
 800  continue
C
      RETURN
      END
