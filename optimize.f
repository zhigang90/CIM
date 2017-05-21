c ==============================================================
c  GEOMETRY OPTIMIZATION MODULE          JB   October 1999
c ==============================================================
c
      subroutine preopt(inp,optcyc,natoms)
      implicit real*8(a-h,o-z)
      character*256 chopval,jobname
c
c  reads the OPTIMIZE line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=26)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      integer optcyc
      character*4 options(nopt)
      character*80 Char
      logical found,Tors,internal
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(natoms)
      DIMENSION XC(3,natoms),IC(natoms,natoms),IAN(natoms)
c ..................................................
      Data nprim/0/, thrbnd/0.85d0/, IType/0/, MaxC/8/
      COMMON /IO/ IOut,ICond                       ! standard/short output
c
      parameter (IUnit=1)
c
      data options/'coor','rege','type','mode','gdii','dmax','gtol',
     $             'dtol','etol','stol','upda','proj','tran','optc',
     $             'ctol','line','scal','back','cuto','hcnv','noto',
     $             'hess','prin','file','qmmm','opti'/
      data ioptyp/21,21,21, 1, 1,11,11,11,11,11,21,21, 0, 1,11,11,11,21,
     $            11, 0, 0,21, 1,21, 0, 0/
c
      Common /job/jobname,lenJ
c
c
c -- get standard and short output
c
       IOut=ioutfil('iout')
       ICond=ioutfil('icond')
       call setival('isumscf',1)      ! switch off SCF summary print
c
      internal=.TRUE.
      Tors=.TRUE.
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c ...........................................................
      IF(optcyc.gt.0) THEN         ! only interpret options first time
c
c -- skip over any additional options directly in input file
 10    CONTINUE
       READ(inp,900) Char
       If(Char(1:1).NE.'$') Then
        BACKSPACE inp
        RETURN
       EndIf
c
c -- if reach here, must be extra embedded options
c -- continue reading lines until find a terminating $
 20    CONTINUE
       READ(inp,900,END=95,ERR=95) Char
       If(Char(1:1).EQ.'$') GO TO 10
       GO TO 20
      ENDIF
c ...........................................................
c
c
c -- open the <opt> file
c -- first check if constraints or additional options need to
c -- be included from a separate file
c
      If(ifound(24).eq.1) Then
       call rmblan(chopval(24),256,Len)
       Call CopyFile(chopval(24)(1:Len),Len,jobname(1:lenJ)//'.opt2',
     $               lenJ+5)
      EndIf
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='UNKNOWN')
c
      WRITE(IUnit,'(a)') '$optimize'
c
c -- optimization coordinates
      IF(ifound(1).eq.1) THEN
        call lowerca2(chopval(1),7)
        If(chopval(1)(1:4).eq.'cart') Then
          IVal = 0
          internal = .FALSE.
        Else If(chopval(1)(1:3).eq.'int'.or.
     $          chopval(1)(1:5).eq.'deloc') Then
          IVal = 1
          IType = 0
        Else If(chopval(1)(1:7).eq.'surface') Then
          IVal = 1
          IType = 1
        Else If(chopval(1)(1:7).eq.'cluster') Then
          IVal = 1
          IType = 2
        Else If(chopval(1)(1:4).eq.'zmat') Then
          IVal = 2
          internal = .FALSE.
        EndIf
        WRITE(IUnit,1000) IVal
        WRITE(IUnit,1010) IType
      ENDIF
c
c --  use of torsions when generating natural internals
      If(ifound(21).eq.1) Then
        Tors = .FALSE.
        WRITE(IUnit,'(a)') '$notors'
      EndIf
c
c -- regeneration of active space with delocalized internals
      IF(ifound(2).eq.1) THEN
        call lowerca2(chopval(2),3)
        If(chopval(2)(1:3).eq.'all') Then
          WRITE(IUnit,'(a)') '$regenerate  2'
        Else
          WRITE(IUnit,'(a)') '$regenerate  1'
        EndIf
      ENDIF
c
c -- stationary point sought (minimum or TS)
      IF(ifound(3).eq.1) THEN
        call lowerca2(chopval(3),3)
        If(chopval(3)(1:2).eq.'ts')
     $     WRITE(IUnit,'(a)') '$transition_state'
      ENDIF
c
c -- mode following during TS search
      If(ifound(4).eq.1) WRITE(IUnit,1100) iopval(1,4)
c
c -- gdiis
      IF(ifound(5).eq.1) THEN
        idiis = iopval(1,5)
        If(idiis.eq.0) idiis = -1
        WRITE(IUnit,1200) idiis
      ENDIF
c
c -- maximum allowed stepsize
      If(ifound(6).eq.1) WRITE(IUnit,1300) ropval(1,6)
c
c -- gradient convergence criterion
      If(ifound(7).eq.1) WRITE(IUnit,1400) '$gtol',ropval(1,7)
c
c -- displacement convergence criterion
      If(ifound(8).eq.1) WRITE(IUnit,1400) '$dtol',ropval(1,8)
c
c -- energy convergence criterion
      If(ifound(9).eq.1) WRITE(IUnit,1500) ropval(1,9)
c
c -- rms gradient criterion for steepest descent step
      If(ifound(10).eq.1) WRITE(IUnit,1600) ropval(1,10)
c
c -- hessian update
      IF(ifound(11).eq.1) THEN
        call lowerca2(chopval(11),9)
        If(chopval(11)(1:2).eq.'no') Then
          IVal = 0
        Else If(chopval(11)(1:2).eq.'ms') Then
          IVal = 1
        Else If(chopval(11)(1:6).eq.'powell') Then
          IVal = 2
        Else If(chopval(11)(1:6).eq.'bofill') Then
          IVal = 3
        Else If(chopval(11)(1:9).eq.'bfgs-safe') Then
          IVal = 5
        Else If(chopval(11)(1:4).eq.'bfgs') Then
          IVal = 4
        EndIf
        WRITE(IUnit,1700) IVal
      ENDIF
c
c -- projection of translations/rotations from Cartesian Hessian
      IF(ifound(12).eq.1) THEN
        call lowerca2(chopval(12),7)
        If(chopval(12)(1:2).eq.'no') Then
          IVal = 0
        Else If(chopval(12)(1:7).eq.'partial') Then
          IVal = 1
        Else If(chopval(11)(1:4).eq.'full') Then
          IVal = 2
        EndIf
        WRITE(IUnit,1800) IVal
      ENDIF
c
c -- inclusion of gradient term in Hessian transformation
      If(ifound(13).eq.1) WRITE(IUnit,'(a)') '$gradterm'
c
c -- number of optimization cycles
      If(ifound(14).eq.1) WRITE(IUnit,1900) iopval(1,14)
c
c -- criterion for satisfying constraints
      If(ifound(15).eq.1) WRITE(IUnit,1400) '$ctol',ropval(1,15)
c
c -- cutoff for near-linear bond angle (default 165 degrees)
      If(ifound(16).eq.1) WRITE(IUnit,2000) ropval(1,16)
c
c -- scale factor for inverse-distance coordinates
c
      If(ifound(17).eq.1) WRITE(IUnit,2100) ropval(1,17)
c
c -- backtransformation option for delocalized internal coordinates
c
      IF(ifound(18).eq.1) THEN
        call lowerca2(chopval(18),4)
        IVal = 0
        If(chopval(18)(1:4).eq.'full') IVal=1
        WRITE(IUnit,2200) IVal
      ENDIF
c
c -- distance cutoff for bonding in surface/cluster optimization
      If(ifound(19).eq.1) WRITE(IUnit,2300) ropval(1,19)
c
c -- convert to Cartesian Hessian every cycle
      If(ifound(20).eq.1) WRITE(IUnit,'(a)') '$hcnvrt'
c
c -- starting hessian option
      IF(ifound(22).eq.1) THEN
        call lowerca2(chopval(22),4)
        If(chopval(22)(1:4).eq.'unit') Then
          IVal = -2
        Else If(chopval(22)(1:4).eq.'defa') Then
          IVal = -3
        EndIf
        WRITE(IUnit,2400) IVal
      ENDIF
c
c -- print flag (for optimization section only)
      If(ifound(23).eq.1) WRITE(IUnit,2500) iopval(1,23)
c
c -- qmmmm flag
      If(ifound(25).eq.1) WRITE(IUnit,'(a)') '$qmmm'
c
c
c -- now read input file to see if any additional lines of data
c    follow the OPTImize command line
c
 30   CONTINUE
      READ(inp,900) Char
      If(Char(1:1).NE.'$') Then
       BACKSPACE inp
       GO TO 50
      EndIf
c
c -- write line directly into <opt> file
      call rmblan(Char,80,len)
      write(IUnit,'(a)') Char(1:len)
c
c -- if reach here, must be extra embedded options
c -- continue reading lines until find a terminating $
c
 40   CONTINUE
      READ(inp,900,END=95,ERR=95) Char
      call rmblan(Char,80,len)
      write(IUnit,'(a)') Char(1:len)
      If(Char(1:1).EQ.'$') GO TO 30
      GO TO 40
c
 50   CONTINUE
c
c
c -- finished options
c -- append any previous optimization data
c
      inquire(file=jobname(1:lenJ)//'.opt2',exist=found)
      If(found) Then
        OPEN (Unit=40,FILE=jobname(1:lenJ)//'.opt2',
     $        FORM='FORMATTED',STATUS='OLD')
        Call AppendFile(40,IUnit,-1)
        CLOSE (Unit=40,STATUS='DELETE')
      EndIf
c
      WRITE(IUnit,'(a)') '$end'
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c -------------------------------------------------------------------
c -- determine in advance how many primitive internals
      IF(internal) THEN
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.coord',
     $        FORM='FORMATTED',STATUS='OLD')
        CALL RdCoord(40,natoms,AtSymb,XC,-1,idum)
        CLOSE (UNIT=40,STATUS='KEEP')
        If(chopval(1)(1:7).NE.'cluster') Then
          CALL GetAtNo(natoms,AtSymb,IAN)
          CALL IZeroIT(IC,natoms*natoms)
c -- check for additional user-defined connectivity
          chopval(24) = jobname(1:lenJ)//'.opt'
          CALL ChkCNNCT(chopval(24),NAtoms,MaxC,IC)
          CALL CONNECTM(natoms,IAN,XC,thrbnd,IC,IErr)
          If(IErr.NE.0) call nerror(1,'preopt',
     $              '**ERROR FORMING CONNECTIVITY MATRIX**',0,0)
          CALL GetNPrim(natoms,Tors,IC,nprim)
c -- allow for additional constraints, replacement of 180 degree
c -- bends etc.. which may INCREASE number of primitives
          nprim = MAX(nprim+2*natoms,2*nprim)
        Else
c -- estimate number of primitives based on distance and cutoff
          CutOff = 5.0d0
          if(ifound(19).eq.1) CutOff=ropval(1,19)
          CALL EstNPrim(natoms,XC,CutOff,nprim)
c -- allow for additional constraints, change in underlying space, replacement
c -- of 180 degree bends etc.. which may INCREASE number of primitives
          nprim = nprim+5*natoms
          write(IOut,*) ' Estimated maximum number of primitives is:',
     $                    nprim
        EndIf
      ENDIF
c
      call setival('nprim',nprim)     ! store in depository
c -------------------------------------------------------------------
c
      RETURN
c
c .....................................................
c  Error Section
c
 95   CONTINUE
      WRITE(IOut,3000)
      CALL OptExit(9)
c
  900 Format(A80)
 1000 FORMAT('$optcoord  ',I2)
 1010 FORMAT('$coordtype  ',I2)
 1100 FORMAT('$mode  ',I4)
 1200 FORMAT('$gdiis  ',I4)
 1300 FORMAT('$stepsize  ',F9.6)
 1400 FORMAT(A5,2X,F9.6)
 1500 FORMAT('$etol  ',F10.8)
 1600 FORMAT('$stol  ',F6.3)
 1700 FORMAT('$update  ',I2)
 1800 FORMAT('$project  ',I2)
 1900 FORMAT('$optcycle  ',I4)
 2000 FORMAT('$linear  ',F7.3)
 2100 FORMAT('$bskal  ',F7.3)
 2200 FORMAT('$backtrans  ',I2)
 2300 FORMAT('$cutoff  ',F9.6)
 2400 FORMAT('$hessian  ',I2)
 2500 FORMAT('$print  ',I2)
 3000 FORMAT(/,2X,'***ERROR*** Problems Reading Optimization Options',/,
     $            '   Check Your Input File')
c
      end
c  =======================================================================
c
      SUBROUTINE OPTGEOM(NMem,Z,Cnvgd)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** GEOMETRY OPTIMIZATION PROGRAM **
C
C  This Program will locate minima and transition states
C  using either Cartesian coordinates, Z-matrix internal
C  coordinates or a set of delocalized internal coordinates
C  which are generated automatically.
C  It also handles fixed distance, bond angle and dihedral
C  angle constraints in either Cartesian or (if appropriate)
C  internal coordinates.
C  ..........................................................
C
C  OPTGEOM itself is simply a "wrapper" which reads the
C  various input files to determine the size of the current
C  system and the optimization options requested and, based
C  on this information, allocates the necessary memory.
C  It then calls OPTMAIN which is itself another "wrapper"
C  which completes job input and is responsible for all
C  other I/O. OPTMAIN calls <optimize> which is the main
C  driving routine.
C
C  FILES
C
C  Data is transferred to and from OPTGEOM via the following
C  files:
C
C  <sym>      -  number of atoms, symmetry data
C  <coord>    -  current geometry (Cartesian coordinates)
C  <zmat>     -  current geometry (Z-matrix)
C  <control>  -  current energy
C  <grad>     -  current gradient
C  <hess>     -  current Hessian (optional)
C  <opt>      -  optimization options (optional)
C  <optchk>   -  information pertinent to next optimization cycle
C                (produced by OPTGEOM itself)
C
C  At the start of a new job, only the <control> and <coord>/<zmat>
C  files are mandatory. Optimization can proceed without Hessian
C  (second derivative) information, but if this is available it
C  really should be utilized. For TS searches and all optimizations
C  in Cartesian coordinates a good initial Hessian is almost
C  mandatory. Often a simple mechanics Hessian is perfectly adequate.
C  The <opt> file is available to request non-standard optimization
C  options or to impose geometrical constraints; for "standard"
C  minimizations the default parameters are normally sufficient
C  and the <opt> file is not needed.
C  The <optchk> file is used to communicate data between
C  optimization cycles and for restarts.
C  .............................................................
C
C
      CHARACTER cdum*20,jobname*256
      LOGICAL Cnvgd,internal,zmat,scan,force
C
      DIMENSION Z(NMem)
C
      PARAMETER (MaxL=8)         !  maximum # of real atoms used to
c                                   define position of a dummy atom + 1
      PARAMETER (thrbnd=0.85d0)  !  distance ratio for bonding
c
      Data IUnit/1/              !  unit number for checkpoint I/O
c
      COMMON /IO/ IOut,ICond     ! standard/short output
      Common /job/jobname,lenJ
C
C
C --------------------------------------------------------------
C  ** CHECK LICENSE **
C
cc      CALL ChkLicense
C --------------------------------------------------------------
C
c -- initialize cpu & elapsed time
      call secund(tcs)
      call elapsec(tes)
C
C  Read from the <control> file
C    total number of atomic centres, including dummies
C    total number of molecules
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      call rdcntrl(IUnit,5,'$nmol',1,NMol,rdum,cdum)
c -- check for potential scan
      CALL RdSCAN0(IUnit,zmat,npcs,force,scan)
c -- if we have an external force scan, do nothing
      If(force) scan=.False.
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Read from the <sym> file
C    number of atoms
C    number of symmetry operations
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,rdum,cdum)
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Is there a Z-matrix?
C
      CALL ScanZMAT(NZ,IErr)
C
C  Read from the <opt> file those optimization parameters that
C  affect memory requirements
C
C  These are:
C
C  IOptC    -  optimization in Cartesian, Z-matrix or natural
C              internal coordinates
C  IGen     -  whether to generate new delocalized internals each cycle
C  IType    -  internal coordinate type (normal/surface/cluster)
C  NFix     -  number of "fixed" atoms
C  NDum     -  number of dummy atoms
C  NCons    -  total number of geometric constraints to be imposed
C  NComp    -  number of composite constraints
C  NPComp   -  total number of primitives in composite constraints
C  NCnnct   -  number of lines of atomic connectivity data
C  NDrive   -  number of coordinates to drive
C  MaxDiis  -  maximum allowed size of GDIIS subspace
C
      CALL RdOPT0(IOptC,  IGen,   IType,  NFix,   NDum,
     $            NCons,  NComp,  NPComp, NCnnct, NDrive,
     $            MaxDiis)
c
c
      If(MaxDiis.EQ.-1) MaxDiis = 4          !  maximum allowed default value
c
      internal = Abs(IOptC).EQ.1.OR.IOptC.EQ.4
c
      NCNTR = MIN(2*NAtom,NAtom+10)          !  maximum # of atomic centres
      NCNTR = MAX(NCNTR,NAtoms+3*NDum+NCons)
      NAT3 = 3*NAtoms
      NDim = 3*NCNTR                         !  maximum possible dimension
      MCons = NCons
      MPComp = NPComp
      MComp = NComp
      If(scan) Then
        MCons = MCons+1
        If(npcs.GT.0) Then
          MComp = MComp+1
          MPComp = MPComp+npcs
        EndIf
      EndIf
c
      If(internal) Then
c .....................................................................
c -- get number of primitives from Texas depository
        call getival('nprim',NIC)
c .....................................................................
        NZ = NAtoms               !  for back-transformation
c -- possible modification of NIC
cc        If(IType.EQ.1) Then
          If(NFix.GT.0) MCons = NCons + NIC   ! allow for surface constraints
          If(IType.EQ.1) NIC = NIC + 99       ! allow for connecting bonds
cc        EndIf
      EndIf
C ...............................................................
C
C  Scratch Storage
C  allocate here ALL scratch storage
C
      NScr = MAX(3*NZ,3*NAT3*NAT3)
      mem0 = 4*NDim*NDim + NDim*MCons + NDim + NAT3
     $            + MAX(NAT3*NAT3,MaxDiis+1)
      If(MaxDiis.GT.1) mem0 = mem0 + (MaxDiis+1)**2 +
     $                          NAT3*MaxDiis + MAX(NAT3,MaxDiis+1)
c
      IF(internal) THEN
        imax = MAX(NAT3,NIC)
        mem1 = NAT3**2 + 2*NAT3*NIC + 2*NIC**2 + 2*NIC + imax*NIC
cc        mem1 = NAT3**2 + 26*NIC + (NIC*(NIC+1))/2
        mem2 = NAT3*NIC + 3*NAT3**2
        IMem = MAX(mem0,mem1,mem2)
c
        If(IGen.GT.0) IMem = IMem + NAT3*NIC
        If(IGen.GT.1) IMem = IMem + 7*NIC
c
        NScr = NScr + IMem
      ELSE IF(Abs(IOptC).EQ.2) THEN
        NScr = NScr + mem0 + 4*(9*NZ*NZ) + 6*NZ
      ELSE
        NScr = NScr + mem0
      ENDIF
c .....................................................................
C
C
C  Now get the memory
C
      MDim = MAX(NDim,MCons)
c
      IMem = 8*NDim + MDim + NAT3**2 + 2*NAtoms + 13*NZ
     $         + 9*NTrans + NAtoms*NTrans + NScr
      If(Abs(IOptC).GT.0) IMem = IMem + 3*NDim + MDim + NDim*NDim
      If(internal) Then
         IMem = IMem + 8*NIC + NAtoms**2
         If(IGen.GE.0)  IMem = IMem + NDim*NIC
         If(IGen.EQ.-1) IMem = IMem + NIC + NAT3 + 24*NIC
      EndIf
      If(IGen.GT.0) IMem = IMem + NIC*NIC
      If(MCons.GT.0) IMem = IMem + 10*MCons + 7*MPComp + MComp + 3
      If(NFix+NDum.GT.0) IMem = IMem + NDim
      If(NDum.GT.0) IMem = IMem + NDum*MaxL
      If(NDrive.GT.0) IMem = IMem + 7*NDrive
      If(NTrans.GT.1) IMem = IMem + 9*NTrans + NAtoms*NTrans
      If(MaxDiis.GT.0) IMem = IMem + 2*NAT3*MaxDiis
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,7,'OPTGEOM')
C
C  Allocate memory pointers
C
C  (a) storage for cartesian coordinates
C
      IXC = iptr                 !  current coordinates
      IXO = IXC + NDim           !  previous coordinates
      IMOL = IXO + NDim          !  pointer to start/end of molecules
      IXCG = IMOL + NAtoms       !  atomic charges
      IXCM = IXCG + NDim         !  atomic masses
      IGC = IXCM + NDim          !  current gradient
      IGO = IGC + NDim           !  previous gradient
      IHS = IGO + NDim           !  Hessian matrix
      IDD = IHS + NAT3*NAT3      !  displacement vector
      IVM = IDD + NDim           !  eigenvector for mode following
      ICNM = IVM + NDim          !  general integer storage
      IEnd = ICNM + MDim
C
C  (b) storage for internal coordinates
C
      IF(Abs(IOptC).GT.0) THEN
c
       IXT = IEnd                !  current geometry
       IGT = IXT + MDim          !  current gradient
       IOT = IGT + NDim          !  previous gradient
       IHT = IOT + NDim          !  Hessian matrix
       IDT = IHT + NDim*NDim     !  displacement vector
       IEnd = IDT + NDim
      ELSE
       IXT = 1
       IGT = 1
       IOT = 1
       IHT = 1
       IDT = 1
c
      ENDIF
C
C  (c) storage for delocalized internal coordinates
C
      IF(internal) THEN
c
       iccn = IEnd                   !  atomic connectivity matrix
       ityp = iccn + NAtoms**2       !  type of primitive internal
       ilst = ityp + NIC             !  list of atoms in primitives
       icof = ilst + 4*NIC           !  primitive weights
       isav = icof + NIC             !  values of primitive torsions
       ixpm = isav + NIC             !  primitive values
       IEnd = ixpm + NIC
c
       If(IGen.EQ.-1) Then
        inp1 = IEnd                  ! # of primitives for each NIC
        int1 = inp1 + NIC            ! integer*2 indexing array
        iut1 = int1 + 12*NIC         ! non-zero NIC components
        IEnd = iut1 + 12*NIC
        iut  = 1
       Else
        iut  = IEnd
        IEnd = iut + NDim*NIC        ! transformation vectors
        inp1 = 1
        int1 = 1
        iut1 = 1
       EndIf
c
       If(IGen.GT.0) Then
        ihpm = IEnd                  !  Hessian in primitive internals
        IEnd = ihpm + NIC*NIC
       Else
        ihpm = 1
       EndIf
c
      ELSE
       iccn = 1
       ityp = 1
       ilst = 1
       icof = 1
       isav = 1
       iut  = 1
       ixpm = 1
       ihpm = 1
       inp1 = 1
       int1 = 1
       iut1 = 1
c
      ENDIF
C
C  (d) storage for Z-matrix coordinates
C
      IF(NZ.GT.0.OR.internal) THEN
c
       IZO = IEnd                !  Z-matrix parameter values
       IZS = IZO + 3*NZ          !  Z-matrix connectivity data
       IXZ = IZS + 4*NZ          !  Z-matrix coordinates
       IZG = IXZ + 3*NZ          !  parameter optimization array
       IEnd = IZG + 3*NZ
      ELSE
       IZO = 1
       IZS = 1
       IXZ = 1
       IZG = 1
c
      ENDIF
C
C  (e) storage for constraints
C
      IF(MCons.GT.0) THEN
c
       ICT = IEnd                !  type of constraint
       IRC = ICT + MCons         !  constraint value
       ICN = IRC + MCons         !  list of atoms involved in each constraint
       ICM = ICN + 4*MCons       !  compound constraint (internals)
       IEos = ICM + MCons        !  storage for 3 energies
       IGos = IEos + 3           !  storage for 3 sets of Lagrange multipliers
       ICmp = IGos + 3*MCons     !  # of primitives in each composite constraint
       IPcm = ICmp + MComp       !   ditto  list of primitives
       IPWt = IPcm + 5*MPComp    !   ditto  weight of each primitive
       IPNM = IPWt + MPComp      !  mapping of constraints to primitives
       IEnd = IPNM + MPComp
      ELSE
       ICT = 1
       IRC = 1
       ICN = 1
       ICM = 1
       IEos = 1
       IGos = 1
       ICmp = 1
       IPcm = 1
       IPWt = 1
       IPNM = 1
c
      ENDIF
C
C  (f) storage for fixed coordinates
C
      IF(NFix+NDum.GT.0) THEN
c
       IFX = IEnd                !  list of fixed/active coordinates
       IEnd = IFX + NDim
      ELSE
       IFX = 1
c
      ENDIF
C
C  (g) storage for dummy atom dependencies
C
      IF(NDum.GT.0) THEN
c
       IDM = IEnd
       IEnd = IDM + NDum*MaxL
      ELSE
       IDM = 1
c
      ENDIF
C
C  (h) storage for coordinate driving
C
      IF(NDrive.GT.0) THEN
c
       IDRT = IEnd
       IDRV = IDRT + NDrive
       IFRV = IDRV + 4*NDrive
       ILRV = IFRV + NDrive
       IEnd = ILRV + NDrive
      ELSE
       IDRT = 1
       IDRV = 1
       IFRV = 1
       ILRV = 1
c
      ENDIF
C
C  (i) storage for symmetry
C
      IUQ = IEnd                 !  list of symmetry-unique atoms
      ITN = IUQ + NAtoms         !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*NTrans       !  list of atomic equivalences
      IEnd = INQ + NATOMS*NTrans
C
C  (j) storage for GDIIS
C
      IF(MaxDiis.GT.0) THEN
c
       IXS = IEnd                !  storage of MaxDiis previous geometries
       IGS = IXS + NAT3*MaxDiis  !  storage of MaxDiis previous gradients
       IEnd = IGS + NAT3*MaxDiis
      ELSE
       IXS = 1
       IGS = 1
c
      ENDIF
C
C  (j) general scratch storage
C
      IScr = IEnd
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,7,'OPTGEOM')
C
C
c memory status
c -- assign memory for high water mark (old TEXAS)
      call getmem(IEnd,lastx)
C
C  ----------------------------------------------------------------------
C
      CALL OPTMAIN(NDim,    NAtoms,  NIC,     NZ,      MCons,
     $             NCons,   MComp,   NComp,   MPComp,  NPComp,
     $             NDum,    MaxL,    NFix,    NCnnct,  NDrive,
     $             NTrans,  MaxDiis, NMol,    Z(IXC),  Z(IMOL),
     $             Z(IXO),  Z(IXCG), Z(IXCM), Z(IGC),  Z(IGO),
     $             Z(IHS),  Z(IDD),  Z(IVM),  Z(IFX),  Z(ICT),
     $             Z(IRC),  Z(ICN),  Z(ICM),  Z(ICNM), Z(IEos),
     $             Z(IGos), Z(ICmp), Z(IPcm), Z(IPWt), Z(IPNM),
     $             Z(IDM),  Z(IDRT), Z(IDRV), Z(IFRV), Z(ILRV),
     $             Z(IUQ),  Z(ITN),  Z(INQ),  Z(iccn), Z(ityp),
     $             Z(ilst), Z(icof), Z(isav), Z(iut),  Z(inp1),
     $             Z(int1), Z(iut1), Z(ixpm), Z(ihpm), Z(IZO),
     $             Z(IZS),  Z(IXZ),  Z(IZG),  Z(IXT),  Z(IGT),
     $             Z(IOT),  Z(IHT),  Z(IDT),  Z(IXS),  Z(IGS),
     $             npcs,    NScr,    Z(IScr), force,   Cnvgd)
C
C  ----------------------------------------------------------------------
C
      call retmem(1)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
      write(IOut,1100) IEnd,mxmem,memtot
 1100 format(/,' Optimize memory status:',/,
     $         ' memory needed=',i15,' high water=',i15,/,
     $         ' total available memory=',i15)
c
c  ----------------------------------------------------------------------
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE OPTMAIN(NDim,   NAtoms, NIC,    NZ,     MCons,
     $                   NCons,  MComp,  NComp,  MPComp, NPComp,
     $                   NDum,   MaxL,   NFix,   NCnnct, NDrive,
     $                   NTrans, MDiis,  NMol,   XC,     IMOL,
     $                   XOld,   XCharg, XMass,  GC,     GCOld,
     $                   HESS,   D,      VMODE,  IFIX,   ICTYP,
     $                   RCON,   ICON,   LCON,   ICNUM,  EOS,
     $                   GOS,    ICOMP,  IPComp, PCWght, IPCNUM,
     $                   IFUNC,  IDRTYP, IDRIVE, FDRIVE, LDRIVE,
     $                   IUNQ,   TRANS,  NEqATM, ICNNCT, ktyp,
     $                   klist,  Coeff,  SavTOR, UT,     NP1,
     $                   INT1,   UT1,    XPrim,  HPRIM,  GEO,
     $                   IGEO,   XZ,     IG,     XINT,   GINT,
     $                   GOld,   HINT,   DINT,   XSTR,   GSTR,
     $                   npcs,   IMem,   Z,      force,  Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for geometry optimization
C  OPTMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(NDim),XOld(NDim),XCharg(NAtoms),XMass(NAtoms),GC(NDim),
     $       GCOld(NDim),HESS(9*NAtoms*NAtoms),D(NDim),VMODE(NDim)
      DIMENSION ktyp(NIC),klist(4,NIC),Coeff(NIC),SavTOR(NIC),
     $          UT(3*NAtoms*NIC),XPRIM(NIC),HPRIM(NIC*NIC)
      DIMENSION UT1(12*NIC)
      INTEGER NP1(NIC),INT1(12*NIC)
      DIMENSION GEO(NZ,3),IGEO(NZ,4),XZ(3,NZ),IG(3*NZ)
      REAL*8 XINT(NDim),GINT(NDim),GOld(NDim),HINT(NDim,NDim),
     $       DINT(NDim)
      REAL*8 XSTR(3*NAtoms*MDiis),GSTR(3*NAtoms*MDiis)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans),
     $          IUNQ(NAtoms),IMOL(NAtoms)
      DIMENSION ICTYP(MCons),ICON(4,MCons),LCON(MCons),ICNUM(MCons),
     $          ICOMP(MComp),IPComp(5,MPComp),PCWght(MPComp),
     $          IPCNUM(MPComp),IFix(NDim),IFunc(NDum,MaxL),
     $          ICNNCT(NAtoms,NAtoms)
      DIMENSION IDRTYP(NDrive),IDRIVE(4,NDrive),FDRIVE(NDrive),
     $          LDRIVE(NDrive)
      REAL*8 RCON(MCons),EOS(3),GOS(3*MCons),RM(3,3)
      INTEGER IOP(30)
c ............................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NDim),ZSymb(NZ),VARNAM(3*NZ)
      CHARACTER*8 CZ(14*NZ+3*NAtoms)
c ............................................................
      CHARACTER GROUP*4,cdum*20,jobname*256,Coord*4
      LOGICAL Symflag,TSflag,oscil,Steep,Cnvgd,zmat,found,force
C
      DIMENSION Z(IMem)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      COMMON /IO/ IOut,ICond                     ! standard/short output
      Common /job/jobname,lenJ
C
      Data IUnit/1/
C
      PARAMETER (MaxCON=9999999)                 ! maximum # of constraints
      PARAMETER (MaxC=8)                         ! maximum # of bonds per atom
C
      EQUIVALENCE (IOptC,   IOP(1))
      EQUIVALENCE (IGen,    IOP(2))
      EQUIVALENCE (ICons,   IOP(3))
      EQUIVALENCE (MaxDiis, IOP(7))
      EQUIVALENCE (IHess,   IOP(14))
      EQUIVALENCE (IUpDat,  IOP(15))
      EQUIVALENCE (MaxCYC,  IOP(18))
      EQUIVALENCE (LThrsh,  IOP(20))
      EQUIVALENCE (IQmm,    IOP(23))
      EQUIVALENCE (IType,   IOP(28))
      EQUIVALENCE (IPRNT,   IOP(29))
C
C
      AUTOAN = 1.0d0/ANTOAU
      TORAD = PI/180.0d0
      NCNTR = NAtoms+NDum                        ! total # of atomic centres
      NAT3 = 3*NAtoms
      NCon = 0
      Cnvgd = .False.
C
C ...............................................................
C  Read from the <sym> file
C    rotation matrix
C    point group symbol
C    number of degrees of freedom
C    number of symmetry-unique atoms
C    symmetry operations
C    symmetry-equivalent atoms array
C  This information read even if system is C1
C
      Symflag = NTrans.GT.1
      CALL RdSYM(Symflag,NAtoms, RM,     GROUP,  NTrans,
     $           NDEG,   NQ,     IUNQ,   TRANS,  NEqATM)
C
C  Read from <control> file
C   total number of atomic centres (including dummies)
C   current energy (if available)
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,rdum,cdum)
      call fdcntrl(IUnit,7,'$energy',idum)
      If(idum.EQ.0) call rdcntrl(IUnit,7,'$energy',2,idum,EC,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  read initial coordinates from <coord> file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL RdCoordF(IUnit,  NAtoms, AtSymb, XOld,   NMol,
     $              IMOL,   XCharg, XMass)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Attempt to read <optchk> file
C  If no <optchk> file is available, this is the first
C  step of a new job
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.optchk',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
c
      READ(40) NCycle,NAtoms,NDEG,NS,NDum,NFix,NCons,NCon,NComp,
     $         NPComp,NDrive,AtSymb,XC,XOld,XCharg,XMass,NDiis,IOP,EOld,
     $         GC,GCOld,HESS,D,VMODE,CZ,(ICNUM(I),I=1,MAX(NDim,MCons))
      If(MDiis.GT.1) READ(40) XSTR,GSTR
      If(NFix.NE.0) READ(40) IFix
      If(NDum.NE.0) READ(40) IFunc
      If(NCons.NE.0) READ(40) ICTYP,RCON,ICON,LCON
      If(NComp.NE.0) READ(40) ICOMP,IPComp,IPCNUM,PCWght
      If(NDrive.GT.0) READ(40) IDRTYP,IDRIVE,FDRIVE,LDRIVE
      If(Abs(ICons).EQ.1) READ(40) NCount,EOS,GOS,IOS,oscil
      If(IOptC.NE.0) READ(40) XINT,GINT,HINT,DINT
      If(Abs(IOptC).EQ.1.OR.IOptC.EQ.4) Then
         READ(40) NPrim,IOrd,ktyp,klist,Coeff,SavTOR
         If(IGen.EQ.-1) Then
            READ(40) NCmp,NP1,INT1,UT1
         Else
            READ(40) UT
         EndIf
      EndIf
      If(IGen.GT.0) READ(40) IHPrim,HPRIM
      If(Abs(IOptC).EQ.2.OR.Abs(IOptC).EQ.1.OR.IOptC.EQ.4)
     $   READ(40) NVar,ZSymb,GEO,IGEO,XZ,IG,VARNAM
c
      CLOSE (UNIT=40,STATUS='KEEP')
      GO TO 30
C
 10   CONTINUE
C
C  Nothing on file
C  Initialize
C
      CALL ZeroIT(XC,NDim)
      CALL CpyVEC(NAT3,XOld,XC)
C
C  Check for data on OPT file
C
      CALL RdOPT1(IOP)
C
      IF(Abs(IOptC).EQ.2) THEN
C
C  Read Z-Matrix data
C  allocate scratch memory pointers
C
       ismb = 1
       is1 = ismb + NZ
       is2 = is1 + 7*NZ
       IEnd = is2 + 6*NZ - 1
c
       CALL RdZMAT(NZ,  CZ(ismb), CZ(is1),CZ(is2),Z(1),
     $             GEO,   IGEO,   IG,     VARNAM, NVar)
C
C  get atomic symbols from Z-matrix symbols
C  and assign dummy atoms
C
       CALL GetAtSym(NZ,CZ(ismb),ZSymb)
      ENDIF
c
      Symflag = IOP(4).EQ.1
      TSflag = IOP(5).EQ.1
c
      NCycle = -1
      NCount = 0
      IOS = 0
      oscil = .False.
      NDiis = 0
      EOld = 0.0d0
      NPrim = 0
      RMSG = 1.0d0
C
C  Read in any constraint data
C
      If(NFix+NDum.GT.0) CALL IZeroIT(IFix,NDim)
      If(NFix.GT.0.AND.IType.NE.1) CALL RdFIX(NAtoms,NFix,IFix)
      If(NFix.GT.0.AND.IType.EQ.1)
     $   CALL RdSurf(NAtoms,NMol,IMOL,IPRNT,NFix,IFix)
      IF(NDum.GT.0) THEN
       CALL RdDUM(NDum,MaxL,IFunc)
C
C  dummy atoms must be fixed, so make sure they are
C
       NFix = NFix + 3*NDum
       DO 20 I=NAT3+1,3*NCNTR
       IFix(I) = 1
 20    CONTINUE
C
C  give dummy atoms an "atomic symbol"
C
       DO 25 I=NAtoms+1,NCNTR
       AtSymb(I) = 'Du      '
 25    CONTINUE
      ENDIF
c
      If(NCons.GT.0) CALL RdCON(MaxCON, ICons,  NCons,  NComp,  NPComp,
     $                          IOptC,  NAtoms, XC,     ICTYP,  RCON,
     $                          ICON,   ICOMP,  IPComp, PCWght)
C
C  read in any coordinates to drive
C
      If(NDrive.GT.0) CALL RdDRIVE(NDrive, IDRTYP, IDRIVE, FDRIVE)
C
C ---------------------------------------------------------------------
C  Check for existence of <scan> file (from potential scan)
C  (Don't do this if external force scan as there is nothing to do)
C
      inquire(file=jobname(1:lenJ)//'.scan',exist=found)
      IF(found.AND..NOT.force) THEN
C
C  optimized scan only available in delocalized internals
C
       If(Abs(IOptC).NE.1) Then
         WRITE(IOut,1000)
         CALL OptExit(9)
       EndIf
C
C  force internal coordinates only
C
       IOptC = 1
C
C  read scanned primitive from <control> file
C
       ipcs = 1
       iwts = ipcs + 5*npcs
c
       OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $       FORM='FORMATTED',STATUS='OLD')
       call rdscan(IUnit,zmat,Coord,I,J,K,L,npcs,
     $             Z(ipcs),Z(iwts),R1,R2,Step,force)
       CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  add scanned primitive to existing constraints (if any)
C
       If(.NOT.zmat) Then
         CALL AddSCAN(MaxCON, NCons,  NAtoms, XC,     Coord,
     $                I,      J,      K,      L,      NComp,
     $                NPComp, ICOMP,  IPComp, PCWght, npcs,
     $                Z(ipcs),Z(iwts),ICTYP,  RCON,   ICON)
       Else
C
C  Z-matrix scan + optimization in delocalized internals
C  Read Z-Matrix data
C  allocate scratch memory pointers
C
         ismb = 1
         is1 = ismb + NZ
         is2 = is1 + 7*NZ
         IEnd = is2 + 6*NZ - 1
c
         CALL RdZMAT(NZ,  CZ(ismb), CZ(is1),CZ(is2),Z(1),
     $               GEO,   IGEO,   IG,     VARNAM, NVar)
C
C **WARNING**
C   Z-matrix used in optimized potential scan CANNOT have any dummy
C   atoms under the current implementation
C
         If(NZ.NE.NAtoms) Then
           WRITE(IOut,1500)
           CALL OptExit(9)
         EndIf
C
         CALL AddSCANZ(MaxCON, NCons,  NAtoms, XC,     Coord,
     $                 VARNAM, IG,     IGEO,   ICTYP,  RCON,
     $                 ICON)
       EndIf
      ENDIF
C ---------------------------------------------------------------------
C
C  Read in any atomic connectivity data
C
      IF(Abs(IOptC).EQ.1) THEN
       CALL IZeroIT(ICNNCT,NAtoms*NAtoms)
       If(NCnnct.GT.0) CALL RdCnnct(NAtoms,NCnnct,MaxC,IPRNT,Z,ICNNCT)
      ENDIF
C
C  ..................................................................
C    Check for errors and incompatible options
C
      CALL ChkOPT(NCNTR,  NAtoms, NDum,   IOptC,  IType,
     $            ICons,  NCons,  ICTYP,  ICON,   RCON,
     $            LThrsh, NFix,   IFix,   MaxL,   IFunc,
     $            MaxDiis,TSflag, Symflag)
C  ..................................................................
C
C  Check that any imposed constraints do not break symmetry
C
      If(NFix.GT.0) CALL SymFIX(NAtoms,IFix,NTrans,NEqATM)
      If(NCons.GT.0) CALL SymCON(NAtoms,NCons,ICTYP,RCON,ICON,
     $                           NTrans,NEqATM)
      If(NDrive.GT.0) CALL SymDRIVE(NAtoms,NDrive,IDRTYP,FDRIVE,IDRIVE,
     $                              NTrans,NEqATM)
C
 30   CONTINUE
C
C  If job is a QM/MM optimization, need to adjust number of molecules
C
      If(IQmm.EQ.1) Then
        Mol2 = IMOL(2)
        DO 40 I=2,NMol
        IMOL(I) = IMOL(I+1)
 40     CONTINUE
        NMol = NMol-1
      EndIf
C
C  If coordinates AND primitives are being regenerated
C  clear the connectivity matrix
C
      If(IGen.EQ.2.AND.NCycle.NE.-1) Then
        CALL IZeroIT(ICNNCT,NAtoms*NAtoms)
        If(NCnnct.GT.0) CALL RdCnnct(NAtoms,NCnnct,MaxC,IPRNT,Z,ICNNCT)
      EndIf
C
C  Now see if there's a Hessian we can use
C  ** NOTE:  If IUpDAT=0 it is assumed that a calculated Hessian
C            is available every cycle
C
      IF(NCycle.EQ.-1.OR.IUpDAT.EQ.0) THEN
       If(IHess.GT.-2)
     $  CALL RdHESS(jobname(1:lenJ)//'.hess',lenJ+5,NAT3,IPRNT,
     $              HESS,IHess)
       If(IHess.LE.-1) CALL SetDiagMat(NAT3,1.0d0,HESS)
       If(IHess.EQ.-3) IHess = -1
      ENDIF
C
C  Read gradient from <grad> file
C
      If(NCycle.GT.-1) Then
C
C  First check for existence of <grad> file
C  If file does NOT exist, it is assumed that this is a restart
C  and we need to repeat the previous cycle
C
       inquire(file=jobname(1:lenJ)//'.grad',exist=found)
       If(.NOT.found) Then
        WRITE(IOut,1800)
        RETURN
       EndIf
       CALL RdGRAD(NAtoms,GC,'dele')
       CALL VScal(NAT3,-1.0d0,GC)            ! file contains forces
      EndIf
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  set rotation matrix to unity for other than Z-matrix optimization
       If(Abs(IOptC).NE.2) CALL SetDiagMat(3,1.0d0,RM)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C  Now optimize!
C
C ---------------------------------------------------------------
      IF(IGen.EQ.-1) THEN
cc
        CALL OPTIMIZE(NCycle, NAtoms, AtSymb, XC,     XOld,
     $                NS,     NMol,   IMOL,   NDum,   MaxL,
     $                IFunc,  NFix,   IFix,   NCons,  ICTYP,
     $                RCON,   ICON,   LCON,   ICNUM,  NCon,
     $                NComp,  NPComp, ICOMP,  IPComp, PCWght,
     $                IPCNUM, NDrive, IDRTYP, IDRIVE, FDRIVE,
     $                LDRIVE, RM,     GROUP,  NDEG,   NQ,
     $                NTrans, TRANS,  NEqATM, NDiis,  IOP,
     $                EC,     EOld,   GC,     GCOld,  HESS,
     $                D,      VMODE,  NZ,     NVar,   ZSymb,
     $                GEO,    IGEO,   XZ,     IG,     VARNAM,
     $                NIC,    NPrim,  ICNNCT, ktyp,   klist,
     $                Coeff,  SavTOR, NCmp,   NP1,    INT1,
     $                UT1,    XPrim,  HPRIM,  IOrd,   IHPrim,
     $                XINT,   GINT,   GOld,   HINT,   DINT,
     $                XSTR,   GSTR,   NCount, EOS,    GOS,
     $                IOS,    oscil,  Steep,  IMem,   Z,
     $                RMSG,   Cnvgd)
cc
      ELSE
cc
        CALL OPTIMIZE(NCycle, NAtoms, AtSymb, XC,     XOld,
     $                NS,     NMol,   IMOL,   NDum,   MaxL,
     $                IFunc,  NFix,   IFIX,   NCons,  ICTYP,
     $                RCON,   ICON,   LCON,   ICNUM,  NCon,
     $                NComp,  NPComp, ICOMP,  IPComp, PCWght,
     $                IPCNUM, NDrive, IDRTYP, IDRIVE, FDRIVE,
     $                LDRIVE, RM,     GROUP,  NDEG,   NQ,
     $                NTrans, TRANS,  NEqATM, NDiis,  IOP,
     $                EC,     EOld,   GC,     GCOld,  HESS,
     $                D,      VMODE,  NZ,     NVar,   ZSymb,
     $                GEO,    IGEO,   XZ,     IG,     VARNAM,
     $                NIC,    NPrim,  ICNNCT, ktyp,   klist,
     $                Coeff,  SavTOR, NCmp,   NP1,    INT1,
     $                UT,     XPrim,  HPRIM,  IOrd,   IHPrim,
     $                XINT,   GINT,   GOld,   HINT,   DINT,
     $                XSTR,   GSTR,   NCount, EOS,    GOS,
     $                IOS,    oscil,  Steep,  IMem,   Z,
     $                RMSG,   Cnvgd)
cc
      ENDIF
C --------------------------------------------------------------
      IF(Cnvgd) THEN
c -- write Cartesian Hessian
       If(IHess.EQ.1) CALL WrHESS(jobname(1:lenJ)//'.hess',lenJ+5,
     $                            1,NAT3,HESS)
c -- summary printout of fixed atoms
       If(NFix.GT.3*NDum) Then
        WRITE(ICond,*)
        CALL PrntFIX(ICond,NAtoms,IFix,cz,Z)
       EndIf
c -- summary printout of constraints
       If(NCons.GT.0) Then
        WRITE(ICond,*)
        CALL PrntCON(ICond,  NCNTR,  NCons,  XC,     ICTYP,
     $               ICON,   RCON,   ICOMP,  IPCOMP, PCWght)
       EndIf
c -- delete <opt> and <optchk> files
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',FORM='FORMATTED')
       CLOSE (UNIT=40,STATUS='DELETE')
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.optchk',FORM='UNFORMATTED')
       CLOSE (UNIT=40,STATUS='DELETE')
c -- may be force field optimization, so delete <ffchk> files
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.ffchk',FORM='UNFORMATTED')
       CLOSE (UNIT=40,STATUS='DELETE')
c -- if Z-matrix potential scan need to update Z matrix values
       inquire(file=jobname(1:lenJ)//'.scan',exist=found)
       If(found) Then
c -- read scanned coordinate type from <control> file
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $        FORM='FORMATTED',STATUS='OLD')
        call rdscan0(IUnit,zmat,jnk,force,found)
        CLOSE (UNIT=IUnit,STATUS='KEEP')
        If(zmat) CALL UpDateZMAT(NAtoms, XC,   CZ(1),  GEO,    IGEO,
     $                           IG,     VARNAM)
       EndIf
       RETURN
      ENDIF
C -------------------------------------------------------------------------
C
C
C  ** NEW STEP NEEDED **
C  invoke updating
C
      IF(NCycle.GT.0) THEN
       If(IUpDAT.NE.0.AND..NOT.Steep) IHess = 1
       If(IOptC.GT.2) IHess = 2
      ENDIF
C
C  Write <optchk> file
C
 90   CONTINUE
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.optchk',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
c
      WRITE(40) NCycle,NAtoms,NDEG,NS,NDum,NFix,NCons,NCon,NComp,
     $         NPComp,NDrive,AtSymb,XC,XOld,XCharg,XMass,NDiis,IOP,EOld,
     $          GC,GCOLD,HESS,D,VMODE,CZ,(ICNUM(I),I=1,MAX(NDim,MCons))
      If(MDiis.GT.1) WRITE(40) XSTR,GSTR
      If(NFix.NE.0) WRITE(40) IFix
      If(NDum.NE.0) WRITE(40) IFunc
      If(NCons.NE.0) WRITE(40) ICTYP,RCON,ICON,LCON
      If(NComp.NE.0) WRITE(40) ICOMP,IPComp,IPCNUM,PCWght
      If(NDrive.GT.0) WRITE(40) IDRTYP,IDRIVE,FDRIVE,LDRIVE
      If(Abs(ICons).EQ.1) WRITE(40) NCount,EOS,GOS,IOS,oscil
      If(IOptC.NE.0) WRITE(40) XINT,GINT,HINT,DINT
      If(Abs(IOptC).EQ.1.OR.IOptC.EQ.4) Then
         WRITE(40) NPrim,IOrd,ktyp,klist,Coeff,SavTOR
         If(IGen.EQ.-1) Then
            WRITE(40) NCmp,NP1,INT1,UT1
         Else
            WRITE(40) UT
         EndIf
      EndIf
      If(IGen.GT.0) WRITE(40) IHPrim,HPRIM
      If(Abs(IOptC).EQ.2.OR.Abs(IOptC).EQ.1.OR.IOptC.EQ.4)
     $   WRITE(40) NVar,ZSymb,GEO,IGEO,XZ,IG,VARNAM
c
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  If job is a QM/MM optimization, need to restore number of molecules
C
      If(IQmm.EQ.1) Then
        NMol = NMol+1
        DO 41 I=2,NMol
        IMOL(I+1) = IMOL(I)
 41     CONTINUE
        IMOL(2) = Mol2
      EndIf
C
C  Write necessary data to <coord>/<zmat> file
C
      IF(Abs(IOptC).EQ.2) THEN
       CALL WrZMAT(NZ,GEO,IG,VARNAM,.true.)          ! <coord> file DELETED
       CALL WrHESS(jobname(1:lenJ)//'.hint',lenJ+5,2,NVar,HINT)
      ELSE
c -- check if there are dummy atoms whose coordinates need to be written
       If(NAtom.GT.NAtoms) Then
         OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $         FORM='FORMATTED',STATUS='OLD')
         CALL RdCoordF(IUnit,  NAtom,  AtSymb, Z,      NMol,
     $                 IMOL,   XCharg, XMass)
         CLOSE (UNIT=IUnit,STATUS='KEEP')
         CALL CpyVEC(3*(NAtom-NAtoms),Z(3*NAtoms+1),XC(3*NAtoms+1))
       EndIF
       CALL WrCoord(NAtom,  AtSymb, XC,     NMol,   IMOL,
     $              XCharg, XMass)
       If(IHess.EQ.1) CALL WrHESS(jobname(1:lenJ)//'.hess',lenJ+5,
     $                            1,NAT3,HESS)
      ENDIF
C
C  Write RMS gradient to texas depository
C
      call setrval('rmsg',RMSG)
c
c  ---------------------------------------------------------------------
c  write to the <control> file infomation about hessian quality
c    ihessq=-1 (negative) poor, crude obtained from geom.opt
c    ihessq=+1 (positive) good, obtained from Hessian calc.
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      ihessq=-1
      If(IUpDat.EQ.0) ihessq=+1    ! assume no update means good Hessian
      call  wrcntrl(IUnit,9,'$hessqual',1,ihessq,Rdum,Cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c  ----------------------------------------------------------------------
C
C  ...............................................
C     Maximum Optimization Cycles Reached
C
      IF(NCycle.EQ.MaxCYC) THEN
       WRITE(IOut,2000)
       WRITE(ICond,2000)
       Cnvgd = .TRUE.
      ENDIF
C  ...............................................
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Optimized Potential Scan Unavailable',
     $            ' in Z-Matrix coordinates')
 1500 FORMAT(/,2X,'***ERROR*** Optimized Z-Matrix Scan Unavailable',
     $            ' if Z Matrix has Dummy Atoms')
 1800 FORMAT('**WARNING** Gradient file NOT found - Assuming Job is',
     $       ' a restart')
 2000 FORMAT(//,' *******************************************',/,
     $          ' **  MAXIMUM OPTIMIZATION CYCLES REACHED  **',/,
     $          ' *******************************************')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTIMIZE(NCycle, NAtoms, AtSymb, XC,     XOld,
     $                    NS,     NMol,   IMOL,   NDum,   MaxL,
     $                    IFunc,  NFix,   IFIX,   NCons,  ICTYP,
     $                    RCON,   ICON,   LCON,   ICNUM,  NCon,
     $                    NComp,  NPComp, ICOMP,  IPComp, PCWght,
     $                    IPCNUM, NDrive, IDRTYP, IDRIVE, FDRIVE,
     $                    LDRIVE, RM,     GROUP,  NDEG,   NQ,
     $                    NTrans, TRANS,  NEqATM, NDiis,  IOP,
     $                    EC,     EOld,   GC,     GCOld,  HESS,
     $                    D,      VMODE,  NZ,     NVar,   ZSymb,
     $                    GEO,    IGEO,   XZ,     IG,     VARNAM,
     $                    NIC,    NPrim,  ICNNCT, ktyp,   klist,
     $                    Coeff,  SavTOR, NCmp,   NP1,    INT1,
     $                    UT,     XPrim,  HPRIM,  IOrd,   IHPrim,
     $                    XINT,   GINT,   GOld,   HINT,   DINT,
     $                    XSTR,   GSTR,   NCount, EOS,    GOS,
     $                    IOS,    oscil,  Steep,  IMem,   Z,
     $                    RMSG,   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ----------------------------------------------------------
C  MAIN DRIVING ROUTINE FOR GEOMETRY OPTIMIZATION
C  ----------------------------------------------------------
C  This routine will locate minima and transition states
C  using either Cartesian coordinates, Z-Matrix internal
C  coordinates or a set of delocalized internal coordinates
C  which are generated automatically.
C
C  It also handles a variety of constraints (distance, bond
C  angle, out-of-plane bend, torsion) in delocalized internal
C  coordinates or Cartesians. Constraints DO NOT need to be
C  satisfied in the starting geometry
C  -----------------------------------------------------------
C
C
C  ARGUMENTS
C
C  NCycle   -  cycle number
C              i.e. number of times this routine has been called
C              in this job step
C              initial value should be -1 (for coordinate generation)
C              incremented by 1 on exit
C  NAtoms   -  number of real atoms
C  AtSymb   -  atomic symbols (used for nice printout)
C  XC       -  current Cartesian coordinates
C              contains new coordinates on exit
C  XOld     -  Cartesian coordinates on previous cycle
C  NS       -  actual dimension of optimization space
C              (on first entry should be set to zero)
C  NMol     -  number of molecules (for, e.g., cluster optimizations)
C  IMOL     -  pointers to start/end of molecules in XC array
C  ..................................................................
C  NDum     -  number of dummy atoms
C              NOTE: The ONLY use of dummy atoms is to help specify
C                    constraints.
C              All dummy atoms are defined with reference to a list of
C              real atoms and dummy atom coordinates will be generated
C              from the coordinates of the real atoms in its defining
C              list. There are two types of dummy atom:
C                1.  Those positioned at the arithmetic mean of the
C                    real atoms in the defining list
C                2.  Those positioned a unit distance along the normal
C                    to a plane defined by 3 atoms, centred on the
C                    middle atom of the 3
C  MaxL     -  maximum number of real atoms involved in definition
C              of dummy atom position
C  IFunc    -  list for each dummy atom of the atom type
C              and the real atoms defining its position
C
C              if there are no dummy atoms the above array is
C              not needed and can be a dummy argument
C  ..................................................................
C  NFix     -  number of fixed Cartesian coordinates
C  IFIX     -  list of fixed/active coordinates
C               0 - coordinate active
C               1 - coordinate inactive (will be "fixed")
C
C              if no atomic coordinates are fixed the above
C              array is not needed and can be a dummy argument
C  ..................................................................
C  NCons    -  total number of constraints
C  ICTYP    -  constraint type
C               1 - fixed distance               stre
C               2 - fixed bond angle             bend
C               3 - fixed out-of-plane bend      outp
C               4 - fixed dihedral angle         tors
C               5 - fixed coplanar bend          linc
C               6 - fixed perpendicular bend     linp
C               9 - composite constraint
C  RCON     -  constraint value (in atomic units)
C  ICON     -  atoms involved in constraint
C                IC1-IC2           distance constraint
C                IC1-IC2-IC3       bond angle constraint
C                IC1-IC2-IC3-IC4   all other constraints
C  LCON     -  which constrains are "active" (independent)
C              LCON(I) = 0  constraint is active
C              LCON(I) = 1  constraint should be eliminated
C  ICNUM    -  indicates which primitive in list corresponds
C              to the desired constraints
C  NCon     -  on exit contains actual number of constraints
C              (some may be equivalent by symmetry or otherwise
C               not independent)
C  NComp    -  number of composite constraints
C  NPComp   -  number of primitives in composite constraints
C  ICOMP    -  number of primitives in each composite constraint
C  IPComp   -  constraint type and atoms involved in constraint
C               IPComp(1,I) - constraint type (same definition as ICTYP array)
C               IPComp(2,I) to IPComp(5,I) - atoms in constraint
C  PCWght   -  weight of each primitive in composite constraint
C             (if no value given, assumed to be unity)
C  IPCNUM   -  indicates which primitive in list corresponds
C              to the desired composite constraints
C
C              if constraints are not being imposed none of the
C              above arrays are needed and can be dummy arguments
C  ..................................................................
C  NDrive   -  number of coordinates to be driven
C  IDRTYP   -  type of each primitive to drive
C               1 - distance               stre
C               2 - bond angle             bend
C               3 - out-of-plane bend      outp
C               4 - dihedral angle         tors
C  IDRIVE   -  definition of each primitive (i.e., list of atoms)
C  FDRIVE   -  force (in internal coordinates) to be applied
C  LDRIVE   -  primitive number, i.e., which primitive we are
C              referring to in list of all primitives
C
C              if there is no coordinate driving none of the
C              above arrays are needed and can be dummy arguments
C  ..................................................................
C  RM       -  rotation matrix
C              i.e. whether any axes reorientation occured
C              during symmetry checking
C  GROUP    -  molecular point group
C  NDEG     -  number of internal degrees of freedom
C  NQ       -  number of symmetry unique atoms
C  NTrans   -  number of symmetry operations
C  TRANS    -  symmetry operations as 3x3 transformation matrices
C  NEqATM   -  list of atomic equivalences under symmetry operations
C
C              if symmetry is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ..................................................................
C  NDiis    -  current size of GDIIS subspace
C              may be incremented/decremented on exit
C  ..................................................................
C
C  IOP
C
C  Array containing all user definable options
C  These are:
C
C  IOP(1)       IOptC      COORDINATE TYPE
C    0 - optimize in Cartesian coordinates
C    1 - generate and optimize in internal coordinates
C        if this fails abort
C   -1 - generate and optimize in internal coordinates
C        if this fails during any stage of the optimization
C        switch to Cartesians and continue
C    2 - optimize in Z-matrix coordinates
C        if this fails abort
C   -2 - optimize in Z-matrix coordinates
C        if this fails during any stage of the optimization
C        switch to Cartesians and continue
C    3 - large-molecule optimization in Cartesian coordinates
C    4 - large-molecule optimization in internals
C
C  IOP(2)       IGen        GENERATION OF DELOCALIZED INTERNALS
C   -1 - generate a set of natural internal coordinates
C    0 - generate a set of delocalized internal coordinates
C        and use throughout the optimization
C    1 - generate a new set of non-redundant internals
C        from the same primitive space
C    2 - generate a new set of underlying primitives
C        and a new set of non-redundant internals
C
C  IOP(3)       ICons       CONSTRAINT TYPE FOR CARTESIAN COORDINATES
C    0 - full optimization (no constraints)
C    1 - constrained optimization using Lagrange multipliers
C    2 - constrained optimization using Penalty functions
C   -1 - attempt a constrained optimization using
C        Lagrange multipliers; if this fails switch
C        to Penalty functions. At convergence, tidy up
C        using Lagrange multipliers if this failed originally
C   -2 - constrained optimization using Penalty functions
C        tidy up at convergence using Lagrange multipliers
C    3 - constrained optimization using delocalized internals
C
C  IOP(4)       Symflag      USE OF SYMMETRY
C    0 - do not use molecular symmetry
C    1 - make use of point group symmetry
C
C  IOP(5)       TSflag       TYPE OF STATIONARY POINT SOUGHT
C    0 - minimum
C    1 - transition state
C
C  IOP(6)       mode         HESSIAN MODE FOLLOWED DURING TS SEARCH
C    0 - mode following switched off (follow lowest mode)
C    N - maximize along mode N
C
C  IOP(7)       MaxDiis      MAXIMUM SIZE OF SUBSPACE FOR GDIIS
C    0 - do not use GDIIS
C   -1 - default size   MIN(NDEG,NATOMS,4)
C    N - size specified by user (DO NOT set N too large)
C
C  IOP(8)       TolG         CONVERGENCE ON MAXIMUM GRADIENT COMPONENT
C  IOP(9)       TolD         CONVERGENCE ON MAXIMUM ATOMIC DISPLACEMENT
C  IOP(10)      TolE         CONVERGENCE ON ENERGY CHANGE
C        To converge TolG MUST be satisfied (default 0.0003) and
C        ONE of TolD (default 0.0003) or TolE (default 0.000001)
C
C  IOP(11)      STol         GRADIENT TOLERANCE FOR STEEPEST DESCENT STEP
C        If the RMS gradient for an unconstrained optimization is
C        greater that STol, a steepest descent step will be taken
C
C  IOP(12)      DMAX         MAXIMUM ALLOWED STEPSIZE
C
C  IOP(14)      IHess        HESSIAN STATUS
C    0 - have "exact" or initial Cartesian Hessian
C        use as is for Cartesian; transform if internals
C    1 - have Hessian from previous step   need to update
C   -1 - set up default (diagonal) Hessian
C   -2 - use unit Hessian
C
C  IOP(15)      IUpDat       HESSIAN UPDATE
C    0 - do not update the Hessian (!?)
C   -1 - use default update
C    1 - Murtagh-Sargent update
C    2 - Powell update
C    3 - Powell/Murtagh-Sargent update (default for TS)
C    4 - BFGS update (default for minimization)
C    5 - BFGS with safeguards to ensure retention of
C        positive definiteness (default for GDIIS)
C
C  IOP(16)      IProj        CARTESIAN HESSIAN PROJECTION
C    0 - do not project; leave Hessian "as is"
C    1 - project out translations
C    2 - project out translations and rotations
C
C  IOP(17)      IDB          dB/dcart TERM IN INTERNAL COORDINATE HESSIAN
C    0 - do not include this term in Hessian transformation
C        (for approximate Hessian & retention of Cartesian
C         Hessian eigenvalue structure)
C    1 - include it
C
C  IOP(18)      MaxCYC       MAXIMUM NUMBER OF OPTIMIZATION CYCLES
C
C  IOP(19)      CTol         TOLERANCE FOR SATISFYING CONSTRAINT
C        If all imposed constraints are satisfied to better than CTol
C        (default 1.0d-6) then can switch to internal coordinates
C
C  IOP(20)      PThrsh       THRESHOLD FOR NEAR-LINEAR BOND ANGLE
C        (default is 165 degrees)
C
C  IOP(21)      BSkal        SCALE FACTOR FOR INVERSE-DISTANCE COORDINATES
C        (default is 1 - no scaling)

C  IOP(22)      ITors        USE OF TORSIONS IN DELOCALIZED INTERNALS
C    0 - use torsions throughout
C    1 - do not use torsions for FIRST molecule only
C        (in surface optimizations, the first "molecule" is the surface)
C
C  IOP(25)      IBack        BACK-TRANSFORMATION TYPE
C    0 - try Z-matrix iterative back-transformation     O(N)
C    1 - normal full iterative back-transformation      O(N**3)
C
C  IOP(26)      CutOff       "BONDING" DISTANCE THRESHOLD FOR VDW
C        (default is 3 Angstroms for surface optimization
C                    5 Angstroms for cluster optimization)
C
C  IOP(27)      HCnvrt       CONVERSION OF INTERNAL HESSIAN TO CARTESIANS
C    0 - convert Hessian only at convergence
C    1 - do conversion every optimization cycle
C
C  IOP(28)      IType        FLAG FOR OPTIMIZATION TYPE (INTERNAL COORDINATES)
C    0 - standard molecular optimization
C    1 - surface adsorption/reaction using delocalized internals
C    2 - cluster coordinates (delocalized/inverse-distance)
C    3 - distance-only coordinates
C
C  IOP(29)      IPRNT        AMOUNT OF PRINT OUT
C    0 - NO printout (except for error messages)
C    1 - summary and warning printout only
C    2 - "standard" printout
C    3 - slightly more printout (including gradient)
C    4 - heavier printout (including Hessian)
C    5 - heavier still (includes iterative printout)
C    6 - very heavy (including internal generation)
C    7 - debug printout
C
C  ......................................................................
C
C
C  EC       -  current energy
C  EOld     -  previous energy
C  GC       -  current Cartesian gradient
C  GCOld    -  previous Cartesian gradient
C  HESS     -  Hessian matrix in Cartesian coordinates
C              (the full 3*NATOMS by 3*NATOMS matrix)
C              The status of the Hessian is controlled by IHess
C  D        -  previous displacement vector (step)
C  VMODE    -  Actual mode (vector) followed on previous cycle
C              on exit will contain mode followed on current cycle
C              determined by criterion of maximum overlap with
C              the previous mode (TS search only)
C  ..................................................................
C  NZ       -  number of atomic centres in Z-Matrix
C  NVar     -  number of variables in Z-Matrix
C  ZSymb    -  Z-Matrix symbols (used for nice printout)
C  GEO      -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO     -  Z-matrix connectivity
C  XZ       -  Cartesian coordinates of all centres in
C              Z-Matrix orientation
C  IG       -  array determining what to do with Z-matrix parameter
C                  0 - optimize it
C                  J - assign same value as previous (Jth) variable
C                 -J - assign same value, opposite sign
C               1000 - fixed
C  VARNAM   -  names of all "variables" (including fixed)
C
C              if Z-Matrix optimization is not being used none of
C              the above arrays are needed and can be dummy arguments
C
C  ** NOTE     Some of the Z-matrix arrays are now used for the
C              back-transformation with delocalized internals
C  ..................................................................
C  NIC      -  estimate of the maximum number of primitive internal
C              coordinates (stretches, bends and torsions) in system
C  NPrim    -  on entry size of primitive space on previous cycle
C              on exit size of primitive space on current cycle
C  ICNNCT   -  atomic connectivity matrix
C  intcor   -  actual number of internal coordinates
C  ktyp     -  integer array containing internal coordinate type
C  klist    -  list of atoms involved in each primitive
C  Coeff    -  coefficients of primitives in (normalized)
C              compound internal coordinate
C  SavTOR   -  array for storing primitive torsions
C              (possible sign changes near limiting values)
C  NCmp     -  total number of non-zero natural internal components
C  NP1      -  number of primitives in each NIC
C  INT1     -  indices for all non-zero components per NIC
C  UT       -  non-zero NIC components OR full set of
C              delocalized internal coordinates
C  XPrim    -  values of primitive internals
C  HPRIM    -  Hessian in primitive internal coordinates
C  IOrd     -  flag indicating change of order between original
C              atom order and order in Z-matrix back-transformation
C               0 - no change;  1 - change
C  IHPrim   -  default diagonal elemewnt in primitive Hessian
C              (used in cluster optimization for new inverse stretch)
C               0 - use default
C               1 - use unity
C
C  ** NOTE  on first entry all quantities from intcor - IHPrim
C           are unset. They will be specified if a full set of
C           internal coordinates is successfully generated
C
C  XINT     -  array to hold internal coordinates
C  GINT     -  gradient in internal coordinates
C              (on first entry empty; on subsequent entries
C               contains gradient from previous cycle)
C  GOld     -  internal gradient on previous cycle
C  HINT     -  Hessian in internal coordinates
C  DINT     -  step in internal coordinates
C
C              if the optimization is in Cartesian coordinates the
C              above arrays are not needed and can be dummy arguments
C  .................................................................
C  XSTR     -  Array for storing all geometries in GDIIS subspace
C  GSTR     -  Array for storing all gradients in GDIIS subspace
C
C              if GDIIS is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ................................................................
C  NCount   -  counter for oscillation checking
C  EOS      -  storage for current and two previous energies
C  GOS      -  storage for current and two previous Lagrange gradients
C  IOS      -  number of times oscillatory behaviour observed
C  oscil    -  Logical flag indicating oscillation status
C           -  on entry  .true.  - oscillation on previous cycle
C                        .false. - no oscillation
C              on exit   .true.  - oscillation on this cycle
C                        .false. - no oscillation
C
C              if constraints are not being imposed none of the
C              above are needed and can be dummy arguments
C  ..................................................................
C  Steep    -  Logical flag indicating if steepest descent step taken
C  IMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  Cnvgd    -  Logical flag indicating convergence
C              on exit   .true.  - convergence achieved
C                        .false. - another step needed
C
C
C  References
C  ----------
C
C  "An Algorithm for the Location of Transition States"
C   J.Baker  J.Comp.Chem.  7 (1986) 385
C
C  "Geometry Optimization by Direct Inversion in the Iterative Subspace"
C   P.Csaszar and P.Pulay  J.Mol.Struct.  114 (1984) 31
C
C  "Geometry Optimization in Cartesian Coordinates:
C   The End of the Z-Matrix?"
C   J.Baker and W.J.Hehre  J.Comp.Chem.  12 (1991) 606
C
C  "Geometry Optimization in Cartesian Coordinates:
C   Constrained Optimization"
C   J.Baker  J.Comp.Chem.  13 (1992) 240
C
C  "Techniques for Geometry Optimization: A Comparison
C   of Cartesian and Natural Internal Coordinates"
C   J.Baker  J.Comp.Chem.  14 (1993) 1085
C
C  "Constrained Optimization in Cartesian Coordinates"
C   J.Baker and D.Bergeron  J.Comp.Chem.  14 (1993) 1339
C
C  "Geometry Optimization in Redundant Internal Coordinates"
C   P.Pulay and G.Fogarasi  J.Chem.Phys.  96 (1992) 2856
C
C  "The Generation and Use of Delocalized Internal Coordinates
C   in Geometry Optimization"
C   J.Baker, A.Kessi and B.Delley  J.Chem.Phys.  105 (1996) 192
C
C  "Constrained Optimization in Delocalized Internal Coordinates"
C   J.Baker  J.Comp.Chem.  18 (1997) 1079
C
C  "Geometry Optimization in Delocalized Internal Coordinates:
C   An Efficient Quadratically Scaling Algorithm for Large Molecules"
C   J.Baker, D.Kinghorn and P.Pulay  J.Chem.Phys.  110 (1999) 4986
C
C  "Efficient Geometry Optimization of Molecular Clusters"
C   J.Baker and P.Pulay  J.Comp.Chem.  21 (2000) 69
C
C  -------------------------------------------------------------------
C
C
      REAL*8 XC(*),XOld(*),RM(3,3),GC(*),D(*),GCOld(*),
     $       HESS(3*NAtoms,3*NAtoms),VMODE(*),XINT(NDEG),
     $       GINT(NDEG),GOld(NDEG),HINT(NDEG*NDEG),DINT(NDEG)
      DIMENSION GEO(NZ,3),IGEO(NZ,4),IG(3*NZ)
      CHARACTER*8 ZSymb(NZ),VARNAM(3*NZ)
      DIMENSION ktyp(NIC),klist(4,NIC),Coeff(NIC),SavTOR(NIC),
     $          UT(*),XPrim(NIC),HPRIM(NIC*NIC)
      INTEGER NP1(NIC),INT1(12*NIC),IMOL(NMol+1)
      REAL*8 XSTR(3*NAtoms,*),GSTR(3*NAtoms,*)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION ICTYP(NCons),RCON(NCons),ICON(4,NCons),
     $          ICNUM(NCons),LCON(NCons),ICOMP(NComp),
     $          IPComp(5,NPComp),PCWght(NPComp),IPCNUM(NPComp)
      DIMENSION IDRTYP(NDrive),IDRIVE(4,NDrive),FDRIVE(NDrive),
     $          LDRIVE(NDrive)
      DIMENSION IFIX(*),IFunc(NDum,MaxL),ICNNCT(NAtoms,NAtoms)
      DIMENSION IOP(30),EOS(3),GOS(NCons,3)
      CHARACTER*8 AtSymb(NAtoms+NDum)
      CHARACTER*4 GROUP,group1
      LOGICAL Symflag,TSflag,rotate,oscil,Steep,Cnvgd
      Logical Change,Full,invrt,HCnvrt
C
      DIMENSION Z(IMem)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*1 xyzch(3*NAtoms)
c ..................................................
C
      COMMON /IO/ IOut,ICond                       ! standard/short output
C
      PARAMETER (small=1.0d-8,thrsh=1.0d-5)
      PARAMETER (TORAD=3.14159265358979323844d0/180.0d0)
      PARAMETER (One=1.0d0)
c
      DATA Change/.FALSE./, HCnvrt/.FALSE./
C
C
C  Check options
C  ------------------------------
      IOptC = IOP(1)
      IGen = IOP(2)
      ICons = IOP(3)
      Symflag = IOP(4).EQ.1
      TSflag = IOP(5).EQ.1
      mode = IOP(6)
      MaxDiis = IOP(7)
      TolG = DFloat(IOP(8))*1.0d-6
      TolD = DFloat(IOP(9))*1.0d-6
      TolE = DFloat(IOP(10))*1.0d-8
      STol = DFloat(IOP(11))*1.0d-3
      DMAX = DFloat(IOP(12))*1.0d-3
      IHess = IOP(14)
      IUpDat = IOP(15)
      IProj = IOP(16)
      IDB = IOP(17)
      CTol = DFloat(IOP(19))*1.0d-6
      PThrsh = DFloat(IOP(20))*TORAD
      BSkal = IOP(21)*1.0d-3
      ITors = IOP(22)
      Full = IOP(25).EQ.1.OR.IOptC.EQ.0
      If(Full) IOrd = -1
      CutOff = DFloat(IOP(26))*1.0d-6
      HCnvrt = IOP(27).EQ.1
      IType = IOP(28)
      IPRNT = IOP(29)
C  --------------------------------
C
      NCNTR = NAtoms+NDum           !  total # of atomic centres
      NAT3 = 3*NAtoms
C
C  Diis option
C
      IF(MaxDiis.LT.0) THEN
       MaxDiis = MIN(NATOMS,NDEG,4)
       If(MaxDiis.LT.3) MaxDiis = 3
      ENDIF
C
C  Hessian update
C
      IF(IUpDat.LT.0) THEN
       IUpDat = 4
       If(TSflag) IUpDat = 3
       If(MaxDiis.GT.0) IUpDat = 5
      ENDIF
C
C  Use of Hessian inverse for large molecule optimization
C
      IF(IOptC.GT.2) THEN
       invrt = IHess.NE.-2
       If(IOptC.EQ.3) invrt = IHess.EQ.0
      ENDIF
C
C
C  allocate scratch pointers
C
      IV1 = 1
      IV2 = IV1  + NAT3*NAT3
      IV3 = IV2  + NAT3*NAT3
      IEnd = IV3 + NAT3*NAT3
      JMem = IMem - IEnd       ! remaining scratch memory
C
C  ....................................................................
C    INTERNAL COORDINATE GENERATION
C
      IF( (NCycle.EQ.-1.OR.(IGen.GT.0.AND.NCycle.NE.0) ) .AND.
     $    (Abs(IOptC).EQ.1.OR.IOptC.EQ.4) ) THEN
c
       call secund(t1)
       NPrim1 = NPrim
       LGen = IGen
       If(NCycle.EQ.-1.AND.IGen.NE.-1) LGen=2
c
       If(IPRNT.GT.1) Then
        If(NCycle.GT.0) Then
         WRITE(IOut,901)
        Else
         WRITE(IOut,900)
        EndIf
       EndIf
c
       IB   = 1
       INB  = IB  + 12*NIC
       IEnd = INB + 12*NIC
C
C  If coordinates are being regenerated, allocate memory for
C  duplicate arrays
C
       IF(IGen.GT.0) THEN
        IUT = IEnd
        IKT = IUT + NDEG*NPrim
        IKL = IKT + NPrim
        INDX = IKL + 4*NPrim
        IEnd = INDX + 2*NPrim
C
C  make copy of previous arrays
C
        CALL CpyVEC(NDEG*NPrim,UT,Z(IUT))     ! save previous
        If(IGen.EQ.2) Then
         CALL ICpyVEC(NPrim,ktyp,Z(IKT))
         CALL ICpyVEC(4*NPrim,klist,Z(IKL))
        EndIf
       ENDIF
c
       CALL MemCHK(IMem,IEnd,8,'OPTIMIZE')
       JMem = IMem - IEnd       ! remaining scratch memory
c
       CALL MakeINTC(NAtoms, AtSymb, XC,     NMol,   IMOL,
     $               GROUP,  NDEG,   NCons,  ICTYP,  RCON,
     $               ICON,   ICNUM,  NPComp, IPComp, IPCNUM,
     $               NCon,   NFix,   IFIX,   NDrive, IDRTYP,
     $               IDRIVE, FDRIVE, LDRIVE, NQ,     NTrans,
     $               TRANS,  NEqATM, IPRNT,  NIC,    NPrim,
     $               NPrim0, ICNNCT, LGen,   IType,  ITors,
     $               CutOff, BSkal,  PThrsh, ktyp,   klist,
     $               Z(IB),  Z(INB), Coeff,  SavTOR, NCmp,
     $               NP1,    INT1,   UT,     XPrim,  IGEO,
     $               GEO,    IG,     IOrd,   JMem,   Z(IEnd),
     $               IErr)
C
C  sort internals w.r.t. old internals if newly generated primitives
C
       If(IType.GE.2.AND.IGen.EQ.2.AND.NCycle.GT.0.AND.IErr.EQ.0) Then
          CALL SortINT(NPrim0, NPrim1, Z(IKT), Z(IKL), NPrim,
     $                 ktyp,   klist,  Change, MPrim,  Z(INDX))
c -------------------------------------------
c -- IMPORTANT - total number of primitives must NOT exceed maximum allowed
          If(NPrim.GT.NIC) Then
           WRITE(IOut,2001)
           CALL OptEXit(9)
          EndIf
c ------------------------------------------
       EndIf
       call secund(t2)
       t1 = t2-t1
cc       write(6,*) ' Total time generating internal coordinates is ',t1
C
C  If internal coordinate generation failed we may still
C  switch to Cartesian coordinates
C
       IF(IErr.NE.0) THEN
        IF(IOptC.GT.0.OR.(TSflag.AND.mode.GT.0)) THEN
         WRITE(IOut,2100)
         CALL OptExit(9)
        ELSE
         IOP(1) = 0
         IOptC = 0                        ! switch to cartesians
         NDiis = 0
         If(IPRNT.GT.0) WRITE(IOut,2200)
        ENDIF
       ENDIF
c
      ENDIF
C
C  ....................................................................
C
C  increment NCycle
C
      NCycle = NCycle+1
c
      IF(NCycle.EQ.0) THEN
       If(IPRNT.GT.1) Then
        If(IOptC.EQ.0) WRITE(IOut,910)
        If(Abs(IOptC).EQ.1) Then
         If(IGen.GE.0) Then
          If(IType.LE.1) WRITE(IOut,920)
          If(IType.EQ.2) WRITE(IOut,925)
         EndIf
         If(IGen.EQ.-1) WRITE(IOut,930)
        EndIf
        If(Abs(IOptC).EQ.2) WRITE(IOut,940)
       EndIf
       RETURN
      ENDIF
C
C  ...................................................................
C    ROTATION SECTION - CHECK FOR AXIS REORIENTATION
C
      CALL ChkROT(RM,IPRNT,thrsh,rotate)
C
C  If rotation did occur, need to transform the old gradient, the
C  displacement and the Hessian into the new coordinate frame
C
      IF(rotate.AND.NZ.EQ.0) THEN
       CALL RotVEC(NAtoms,RM,GCOld)
       CALL RotVEC(NAtoms,RM,D)
       CALL RotHES(NAtoms,RM,Z(IV1),HESS)
       If(IOptC.EQ.0.AND.NDiis.GT.0) THEN
        DO 10 I=1,NDiis
        CALL RotVEC(NAtoms,RM,XSTR(1,I))
        CALL RotVEC(NAtoms,RM,GSTR(1,I))
 10     CONTINUE
       EndIf
      ENDIF
C
C  ....................................................................
C    SYMMETRY SECTION
C    make sure the geometry and gradient conform fully
C    to the molecular symmetry
C
      CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IV1),
     $            small,  XC)
      CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IV1),
     $            small,  GC)
C
      IF(IHess.EQ.0) THEN
       CALL CpyVEC(NAT3*NAT3,HESS,Z(IV1))
       CALL SymHES(NAtoms, NTrans, NEqATM, TRANS,  Z(IV1),
     $             Z(IV2), Z(IV3), small,  HESS)
      ENDIF
C
C  ....................................................................
C    DUMMY ATOM SECTION
C    generate Cartesian coordinates for dummy atoms from the
C    real atoms in their defining list
C
      If(NDum.GT.0) CALL PutDUM(NCycle, NAtoms, NDum,   NCons,  XC,
     $                          MaxL,   IFunc,  ICTYP,  RCON,   ICON)
C
C  ....................................................................
C    FIXED CARTESIAN COORDINATE SECTION
C    if there are any fixed coordinates, make sure all relevant
C    vector components are zero
C
      If(IOptC.EQ.0.AND.NFix-3*NDum.GT.0)
     $              CALL ZeroFIX(NAT3,IFix,GC,D,HESS)
C
C  ....................................................................
C    STEEPEST DESCENT SECTION
C    for an unconstrained minimization, if the RMS Cartesian gradient
C    is greater than STol then we will take a simple steepest descent
C    step and NOT update the Hessian
C
      Steep = .FALSE.
      IF(ICons.EQ.0) THEN
       RMSG = SProd(NAT3,GC,GC)
       RMSG = SQRT(RMSG/DFloat(NAtoms))
       Steep = RMSG .GT. STol .AND. .NOT. TSflag
      ENDIF
C
C  ....................................................................
C    HESSIAN SECTION
C
      IF(IHess.EQ.1.AND.IOptC.LE.0.AND..NOT.Steep.AND.Full) THEN
C
C  Update the Cartesian Hessian
C  This is done even if the optimization is performed in internal
C  coordinates in case:
C    (i)   we switch to Cartesians (due to problems with internals)
C    (ii)  we switch to a different set of internals
C
       If(IPRNT.GT.1) WRITE(IOut,1000)
c
       CALL UpdHES(NAT3,   IUpDat, IPRNT,  D,      GC,
     $             GCOld,  Z(IV1), Z(IV2), HESS)
c
      ENDIF
C
C  ....................................................................
C    DISPLAY SECTION
C
      IF(IPRNT.GT.1) THEN
       If(IOptC.EQ.0) Then
        If(NCons.EQ.0) WRITE(IOut,1100)
        If(NCons.NE.0) WRITE(IOut,1200)
       Else If(Abs(IOptC).EQ.1) Then
        If(NCons.EQ.0) Then
         If(IGen.GE.0.AND.IType.EQ.0) WRITE(IOut,1300)
         If(IGen.GE.0.AND.IType.EQ.1) WRITE(IOut,1310)
         If(IGen.GE.0.AND.IType.EQ.2) WRITE(IOut,1320)
         If(IGen.EQ.-1) WRITE(IOut,1330)
        Else
         If(IType.EQ.0) WRITE(IOut,1400)
         If(IType.EQ.1) WRITE(IOut,1410)
         If(IType.EQ.2) WRITE(IOut,1420)
        EndIf
       Else If(Abs(IOptC).EQ.2) Then
        WRITE(IOut,1500)
       Else If(IOptC.EQ.3) Then
        WRITE(IOut,1550)
       Else
        WRITE(IOut,1575)
       EndIf
       If(TSflag) Then
        WRITE(IOut,1600)
       Else
        WRITE(IOut,1700)
       EndIf
       WRITE(IOut,1800) NCycle
       If(Abs(IOptC).EQ.2) Call PrntZMAT(IOut,NZ,ZSymb,GEO,IGEO)
       If(Abs(IOptC).NE.2) THEN
        CALL PrntCAR(IOut,0,NCNTR,AtSymb,XC)
        group1=GROUP
        call uppercase(group1,1)
        WRITE(IOut,1900) group1,NDEG
       EndIf
       WRITE(IOut,2000) EC
       If(NFix.GT.3*NDum) CALL PrntFIX(IOut,NAtoms,IFix,xyzch,Z(IV2))
       If(NCons.GT.0)
     $    CALL PrntCON(IOut,   NCNTR,  NCons,  XC,     ICTYP,
     $                 ICON,   RCON,   ICOMP,  IPCOMP, PCWght)
       If(IPRNT.GT.2) CALL PrntGRD(IOut,NATOMS,AtSymb,GC)
       If(NDrive.GT.0) CALL PrntDRIVE(IOut,NDrive,IDRTYP,IDRIVE,FDRIVE)
      ENDIF
C
C  ...................................................................
C    INTERNAL COORDINATE TRANSFORMATION SECTION
C
      IF(Abs(IOptC).EQ.1.OR.IOptC.EQ.4) THEN
C
C  make a copy of the current internal gradient
C  (assuming there is one)
C
       CALL CpyVEC(NDEG,GINT,GOld)
c
       call secund(t1)
       CALL GetINTL(NCycle, NAtoms, AtSymb, XC,     NMol,
     $              IMOL,   GROUP,  NDEG,   NCons,  ICTYP,
     $              RCON,   ICON,   LCON,   ICNum,  NCon,
     $              NComp,  NPComp, ICOMP,  PCWght, IPCNUM,
     $              NDrive, LDRIVE, FDRIVE, NTrans, TRANS,
     $              NEqATM, Steep,  IPRNT,  IHess,  IUpDat,
     $              IDB,    GC,     HESS,   NPrim,  NPrim1,
     $              MPrim,  IGen,   BSkal,  ktyp,   klist,
     $              Coeff,  SavTOR, NS,     NCmp,   NP1,
     $              INT1,   UT,     XPrim,  HPRIM,  Z(IUT),
     $              Z(INDX),Change, IOrd,   IHPrim, GOld,
     $              DINT,   XINT,   GINT,   HINT,   JMem,
     $              Z(IEnd),IErr)
c
       call secund(t2)
       t1 = t2-t1
cc       write(6,*) ' Total time transforming to internals is ',t1
C
C  Set back-transformation flag
C  (use Z-matrix if possible; otherwise full)
C
       If(IOrd.EQ.-1.OR.NCon.GT.0.OR.NTrans.GT.1.OR.
     $    IType.GT.0) Full=.TRUE.
c
       IF(IErr.NE.0) THEN
        IF(IOptC.GT.0.OR.(TSflag.AND.mode.GT.0).OR..NOT.Full.OR.
     $     NDrive.GT.0) THEN
         WRITE(IOut,2100)
         CALL OptExit(9)
        ELSE
         IOP(1) = 0
         IOptC = 0                        ! switch to cartesians
         NDiis = 0
         If(IPRNT.GT.0) WRITE(IOut,2200)
        ENDIF
       ENDIF
c
       If(NCycle.EQ.1.AND.NCon.GT.0) CALL ZeroIT(XINT(NDEG+1),NCons)
cc
      ENDIF
C
C  ...................................................................
C    OPTIMIZATION SECTION
C
C  save current Cartesian coordinates
C
      CALL CpyVEC(NAT3,XC,XOld)
C
C
 100  CONTINUE
      IF(IOptC.EQ.0) THEN
cc
       IF(ICons.EQ.0) THEN
C
C  ** UNCONSTRAINED/FIXED ATOM OPTIMIZATION IN CARTESIAN COORDINATES **
C
        CALL OPTCART(NCycle, NAtoms, AtSymb, XC,     NFix,
     $               IFix,   Symflag,GROUP,  NDEG,   NTrans,
     $               TRANS,  NEqATM, TSflag, MaxDiis,NDiis,
     $               TolG,   TolD,   TolE,   DMAX,   mode,
     $               IHess,  IProj,  IPRNT,  EC,     EOld,
     $               GC,     HESS,   D,      XSTR,   GSTR,
     $               VMODE,  Steep,  JMem,   Z(IEnd),RMSG,
     $               Cnvgd)
c
       ELSE IF(Abs(ICons).EQ.1) THEN
C
C  ** CONSTRAINED OPTIMIZATION IN CARTESIANS USING LAGRANGE MULTIPLIERS **
C
        CALL OPTCON(NCycle, NAtoms, NDum,   NCons,  AtSymb,
     $              XC,     ICTYP,  RCON,   ICON,   NFix,
     $              IFix,   MaxL,   IFunc,  Symflag,GROUP,
     $              NDEG,   NTrans, TRANS,  NEqATM, TSflag,
     $              TolG,   TolD,   TolE,   DMAX,   mode,
     $              IProj,  IPRNT,  EC,     EOld,   GC,
     $              HESS,   D,      VMODE,  NCount, EOS,
     $              GOS,    IOS,    oscil,  JMem,   Z(IEnd),
     $              RMSG,   Cnvgd,  IErr)
c
        IF(IErr.NE.0) THEN
         IF(ICons.EQ.1.OR.TSflag) THEN
          WRITE(IOut,2300)
          CALL OptExit(9)
         ENDIF
c
         IOP(3) = -2
         ICons = -2                      ! switch to penalty functions
         IHess = 0
         If(IPRNT.GT.0) WRITE(IOut,2400)
c
         GO TO 100
        ENDIF
c
       ELSE IF(Abs(ICons).EQ.2) THEN
C
C  ** CONSTRAINED OPTIMIZATION IN CARTESIANS USING PENALTY FUNCTIONS **
C
C  make copy of current gradient in case we need to tidy
C  up the converged geometry later
C  (the gradient is modified in <OPTPEN>)
C
        If(ICons.EQ.-2) CALL CpyVEC(NAT3,GC,VMODE)
c
        CALL OPTPEN(NCycle, NAtoms, NDum,   NCons,  AtSymb,
     $              XC,     ICTYP,  RCON,   ICON,   NFix,
     $              IFix,   MaxL,   IFunc,  Symflag,GROUP,
     $              NDEG,   NTrans, TRANS,  NEqATM, MaxDiis,
     $              NDiis,  TolG,   TolD,   TolE,   DMAX,
     $              IProj,  IPRNT,  EC,     EOld,   GC,
     $              HESS,   D,      XSTR,   GSTR,   JMem,
     $              Z(IEnd),RMSG,   Cnvgd)
c
         IF(Cnvgd.AND.ICons.EQ.-2) THEN
          If(IPRNT.GT.0) WRITE(IOut,2500)
          IOP(3) = 1
          ICons = 1                       ! tidy up using Lagrange multipliers
          CALL SetDiagMat(NAT3,One,HESS)  ! and a unit Hessian
          CALL CpyVEC(NAT3,VMODE,GC)
          CALL ZeroIT(XC(3*(NAtoms+NDum)+1),NCons)
          NCount = 0
          IOS = 0
          GO TO 100
         ENDIF
c
       ENDIF
C  ---------------------------------------------------------------
cc
      ELSE IF(Abs(IOptC).EQ.1) THEN
cc
C  ** OPTIMIZATION IN DELOCALIZED INTERNAL COORDINATES **
C
C  First check constraints to see which routine to call
C
       CALL ConINT(NAtoms, NCons,  XC,     ICTYP,  RCON,
     $             ICON,   ICOMP,  IPComp, PCWght, LCON,
     $             CTol,   IPRNT,  Z(IV1), Z(IV2), IAlg)
c
       IF(IAlg.EQ.0) THEN
C
C  If there are constraints, all are now satisfied
C  need to take only active part of Hessian
C
        CALL CpyMAT(NS,NDEG,NS,HINT,Z(IV3))
        CALL ZeroIT(DINT(NS+1),NCon)
        CALL OPTINT(NCycle, NS,     TSflag, MaxDiis,NDiis,
     $              TolG,   TolD,   TolE,   DMAX,   IHess,
     $              IUpDat, mode,   IPRNT,  EC,     EOld,
     $              XINT,   GINT,   GOld,   Z(IV3), DINT,
     $              XSTR,   GSTR,   VMODE,  Steep,  JMem,
     $              Z(IEnd),RMSG,   Cnvgd)
       ELSE
        CALL OPTCONINT(NCycle, NDEG,   NS,     NCon,   TSflag,
     $                 TolG,   TolD,   TolE,   DMAX,   mode,
     $                 IPRNT,  EC,     EOld,   Z(IV1), Z(IV2),
     $                 XINT,   GINT,   HINT,   DINT,   VMODE,
     $                 JMem,   Z(IEnd),RMSG,   Cnvgd,  IErr)
c
        IF(IErr.NE.0) THEN
          WRITE(6,2300)
          CALL OptExit(9)
         ENDIF
c
       ENDIF
c
       IF(.NOT.Cnvgd) THEN
C
C  generate Cartesian coordinates from new internals
C
        CALL CpyVEC(NDEG,DINT,Z(IV1))
cc
        IF(Full) THEN
c -- full O(N**3) iterative back transformation
         call secund(t1)
         CALL CpyVEC(NAT3,XC,D)
         CALL GetCART(NAtoms, AtSymb, NDEG,   NPrim,  IGen,
     $                XINT,   Z(IV1), ktyp,   klist,  Coeff,
     $                SavTOR, UT,     NTrans, TRANS,  NEqATM,
     $                BSkal,  IPRNT,  Scal,   XC,     JMem,
     $                Z(IEnd),IErr)
         If(Scal.NE.One) CALL CpyVEC(NDEG,Z(IV1),DINT)
         call secund(t2)
         t1 = t2-t1
cc         write(6,*) ' CPU time for full back-transformation: ',t1
c
         IF(IErr.NE.0) THEN
          IF(IOptC.EQ.1.OR.(TSflag.AND.mode.GT.0)) THEN
           WRITE(IOut,2100)
           CALL OptExit(9)
          ENDIF
c
          IOP(1) = 0
          IOptC = 0                        ! switch to cartesians
          NDiis = 0
          CALL CpyVEC(NAT3,D,XC)
          If(IPRNT.GT.0) WRITE(IOut,2200)
c
          GO TO 100
         ENDIF
c
         If(IPRNT.GT.1) CALL PrntCAR(IOut,0,NCNTR,AtSymb,XC)
C
C  form Cartesian displacement vector for Hessian update
C  on next cycle
C
         CALL GetD(NAT3,IPRNT,XC,D)
cc
        ELSE IF(IType.EQ.3) THEN
c -- distance-only coordinates
         call secund(t1)
         CALL GetCARTR(NAtoms, AtSymb, NDEG,   NPrim,  XINT,
     $                 Z(IV1), ktyp,   klist,  Z(IV2), SavTOR,
     $                 UT,     NTrans, NEqATM, BSkal,  IPRNT,
     $                 XPRIM,  Scal,   XC,     JMem,  Z(IEnd), IErr)
         If(Scal.NE.One) CALL CpyVEC(NDEG,Z(IV1),DINT)
         call secund(t2)
         t1 = t2-t1
cc         write(6,*) ' CPU time for distance back-transformation: ',t1
        ELSE
c -- O(N) Z-matrix back transformation
         call secund(t1)
         CALL GetCARTZ(NAtoms, AtSymb, NDEG,   NPrim,  IGen,
     $                 XINT,   Z(IV1), ktyp,   klist,  Z(IV2),
     $                 SavTOR, UT,     NCmp,   NP1,    INT1,
     $                 UT,     NTrans, NEqATM, IPRNT,  XPRIM,
     $                 IGEO,   GEO,    IG,     IOrd,   XC,
     $                 JMem,   Z(IEnd),IErr)
         call secund(t2)
         t1 = t2-t1
cc         write(6,*) ' CPU time for Z-matrix back-transformation: ',t1
        ENDIF
       ENDIF
C
C  ..................................................................
C  Possible Hessian transformation (Internal --> Cartesian)
C
       If(Cnvgd.OR.HCnvrt) Then
        IB = 1
        IBP = IB + NPrim*NAT3
        IEnd = IBP + NDEG*NAT3
        CALL MemCHK(IMem,IEnd,8,'OPTIMIZE')
c
        CALL HSSCART(NAtoms, NPrim,  NDEG,   IDB,    IPRNT,
     $               XC,     GINT,   HINT,   ktyp,   klist,
     $               Coeff,  SavTOR, BSkal,  UT,     Z(ib),
     $               Z(ibp), Z(ib),  Z(ibp), Z(ib),  HESS)
c
       EndIf
c  ..................................................................
cc
      ELSE IF(Abs(IOptC).EQ.2) THEN
cc
C  ** OPTIMIZATION IN Z-MATRIX COORDINATES **
C
C  prepare for optimization
C  convert Cartesian quantities to Z-Matrix internals
C
       NIC = MAX(1,3*NZ-6)
C
C  first make a copy of the current internal gradient
C  (assuming there is one)
C
       CALL CpyVEC(NVar,GINT,GOld)
c
       CALL GetZINT(NZ,     NAtoms, NIC,    NVar,   IPRNT,
     $              IHess,  ZSymb,  GEO,    IGEO,   IG,
     $              IDB,    GC,     HESS,   RM,     XZ,
     $              Z(IV1), XINT,   GINT,   HINT,   JMem,
     $              Z(IEnd) )
c
       CALL OPTINT(NCycle, NVar,   TSflag, MaxDiis,NDiis,
     $             TolG,   TolD,   TolE,   DMAX,   IHess,
     $             IUpDat, mode,   IPRNT,  EC,     EOld,
     $             XINT,   GINT,   GOld,   HINT,   DINT,
     $             XSTR,   GSTR,   VMODE,  Steep,  JMem,
     $             Z(IEnd),RMSG,   Cnvgd)
C
C  convert new geometry (internal coordinates on XINT)
C  back to Z-Matrix internals
C
       CALL BackToZ(NZ,     NIC,    NVar,   IG,     XINT,
     $              Z(IV1), GEO)
cc
      ELSE
cc
C  ** BFGS OPTIMIZATION FOR LARGE SYSTEMS **
C
       If(IOptC.EQ.3) Then
C
C  Cartesian coordinates
C
       call secund(t1)
        CALL OPTBFGS(NCycle, NAT3,   TolG,   TolD,   TolE,
     $               DMAX,   invrt,  IPRNT,  EC,     EOld,
     $               XC,     GC,     GCOld,  HESS,   Z(IV1),
     $               Z(IV2), D,      RMSG,   Cnvgd)
       call secund(t2)
       t1 = t2-t1
cc       write(6,*) ' CPU time for Cartesian Optimization Step: ',t1
       Else
C
C  delocalized internals
C
       call secund(t1)
        CALL OPTBFGS(NCycle, NDEG,   TolG,   TolD,   TolE,
     $               DMAX,   invrt,  IPRNT,  EC,     EOld,
     $               XINT,   GINT,   GOld,   HINT,   Z(IV1),
     $               Z(IV2), DINT,   RMSG,   Cnvgd)
       call secund(t2)
       t1 = t2-t1
cc       write(6,*) ' CPU time for Optimization Step: ',t1
c
        IF(.NOT.Cnvgd) THEN
C
C  generate Cartesian coordinates from new internals
C
         IF(Full) THEN
c -- full O(N**3) iterative back transformation
          call secund(t1)
          CALL CpyVEC(NAT3,XC,D)
          CALL CpyVEC(NDEG,DINT,Z(IV1))
          CALL GetCART(NAtoms, AtSymb, NDEG,   NPrim,  IGen,
     $                 XINT,   Z(IV1), ktyp,   klist,  Coeff,
     $                 SavTOR, UT,     NTrans, TRANS,  NEqATM,
     $                 BSkal,  IPRNT,  Scal,   XC,     JMem,
     $                 Z(IEnd),IErr)
          If(Scal.NE.One) CALL CpyVEC(NDEG,Z(IV1),DINT)
          call secund(t2)
          t1 = t2-t1
cc          write(6,*) ' CPU time for full back-transformation: ',t1
c
          IF(IErr.NE.0) THEN
           IF(IOptC.EQ.1.OR.(TSflag.AND.mode.GT.0)) THEN
            WRITE(IOut,2100)
            CALL OptExit(9)
           ENDIF
c
           IOP(1) = 3
           IOptC = 3                        ! switch to cartesians
           NDiis = 0
           CALL CpyVEC(NAT3,D,XC)
           If(IPRNT.GT.0) WRITE(IOut,2200)
c
           GO TO 100
          ENDIF
        ELSE
c -- O(N) Z-matrix back transformation
          call secund(t1)
          CALL CpyVEC(NDEG,DINT,Z(IV1))
          CALL GetCARTZ(NAtoms, AtSymb, NDEG,   NPrim,  IGen,
     $                  XINT,   Z(IV1), ktyp,   klist,  Z(IV2),
     $                  SavTOR, UT,     NCmp,   NP1,    INT1,
     $                  UT,     NTrans, NEqATM, IPRNT,  XPRIM,
     $                  IGEO,   GEO,    IG,     IOrd,   XC,
     $                  JMem,   Z(IEnd),IErr)
          call secund(t2)
          t1 = t2-t1
cc          write(6,*) ' CPU time for Z-matrix back-transformation: ',t1
         ENDIF
        ENDIF
       EndIf
cc
      ENDIF
C
C  make a copy of the Cartesian gradient for next Hessian update
C
      CALL CpyVEC(NAT3,GC,GCOld)
      EOld = EC                  ! save current energy
      If(TSflag) IOP(6) = mode   ! save mode being followed
c
      RETURN
c
  900 FORMAT(/,'** GENERATION OF INTERNAL COORDINATES **')
  901 FORMAT(/,'** REGENERATION OF INTERNAL COORDINATES **')
  910 FORMAT(/,' OPTIMIZATION WILL USE CARTESIAN COORDINATES')
  920 FORMAT(/,' OPTIMIZATION WILL USE DELOCALIZED INTERNAL',
     $         ' COORDINATES')
  925 FORMAT(/,' OPTIMIZATION WILL USE DELOCALIZED CLUSTER',
     $         ' COORDINATES')
  930 FORMAT(/,' OPTIMIZATION WILL USE REDUNDANT INTERNAL',
     $         ' COORDINATES')
  940 FORMAT(/,' OPTIMIZATION WILL USE Z-MATRIX COORDINATES')
 1000 FORMAT(/,' Cartesian Hessian Update')
 1100 FORMAT(//,'** GEOMETRY OPTIMIZATION IN CARTESIAN COORDINATES **')
 1200 FORMAT(//,'** CONSTRAINED OPTIMIZATION IN CARTESIAN',
     $            ' COORDINATES **')
 1300 FORMAT(//,'** GEOMETRY OPTIMIZATION IN DELOCALIZED INTERNAL',
     $            ' COORDINATES **')
 1310 FORMAT(//,'** SURFACE OPTIMIZATION IN DELOCALIZED INTERNAL',
     $            ' COORDINATES **')
 1320 FORMAT(//,'** GEOMETRY OPTIMIZATION IN DELOCALIZED CLUSTER',
     $            ' COORDINATES **')
 1330 FORMAT(//,'** GEOMETRY OPTIMIZATION IN REDUNDANT INTERNAL',
     $            ' COORDINATES **')
 1400 FORMAT(//,'** CONSTRAINED OPTIMIZATION IN DELOCALIZED INTERNAL',
     $            ' COORDINATES **')
 1410 FORMAT(//,'** CONSTRAINED SURFACE OPTIMIZATION IN DELOCALIZED',
     $            ' INTERNAL COORDINATES **')
 1420 FORMAT(//,'** CONSTRAINED OPTIMIZATION IN DELOCALIZED CLUSTER',
     $            ' COORDINATES **')
 1500 FORMAT(//,'** GEOMETRY OPTIMIZATION IN Z-MATRIX COORDINATES **')
 1550 FORMAT(//,'** LARGE MOLECULE OPTIMIZATION IN CARTESIAN',
     $            ' COORDINATES **')
 1575 FORMAT(//,'** LARGE MOLECULE OPTIMIZATION IN DELOCALIZED',
     $            ' INTERNAL COORDINATES **')
 1600 FORMAT('   Searching for a Transition State')
 1700 FORMAT('   Searching for a Minimum')
 1800 FORMAT(/,'   Optimization Cycle: ',I3)
 1900 FORMAT('   Point Group: ',A4,'  Number of degrees of',
     $         ' freedom: ',I3,/)
 2000 FORMAT(/,'   Energy is ',F16.9,/)
 2001 FORMAT(/,2X,'***ERROR*** Too Many Primitives in Cluster',
     $            ' Optimization',/,
     $         14X,'Suggest Restart with Current Geometry')
 2100 FORMAT(/,2X,'***ERROR*** Problems with Internal Coordinates',/,
     $         14X,'Terminating Optimization')
 2200 FORMAT('**WARNING** Problems with Internal Coordinates',/,
     $       '  Switching to Cartesian Coordinates')
 2300 FORMAT(/,2X,'***ERROR*** Problems with Lagrangian - Terminating',
     $              ' Optimization')
 2400 FORMAT('**WARNING** Problems with Lagrangian - Switching to',
     $         ' Penalty functions')
 2500 FORMAT(/,'** Tidying up geometry using Lagrange multiplier',
     $           ' algorithm **')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTBFGS(NCycle, N,      TolG,   TolD,   TolE,
     $                   DMAX,   Invrt,  IPRNT,  EC,     EOld,
     $                   X,      g,      gold,   HINV,   V1,
     $                   V2,     d,      RMSG,   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Simple BFGS minimizer
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C              i.e. number of times this routine has been called
C              in this job step
C  N       -  number of variables to be optimized
C             (3*NAtoms for Cartesian optimization)
C  TolG    -  convergence tolerance on maximum gradient component
C  TolD    -  convergence tolerance on maximum displacement
C  TolE    -  convergence tolerance on energy
C              (if energy change from previous cycle is less than
C               TolE and TolG is satisfied, stop regardless of TolD;
C               this prevents unnecessary steps for floppy molecules)
C  DMAX    -  maximum allowed stepsize
C  Invrt   -  Logical flag for inverting initial Hessian
C             (otherwise unit matrix assumed)
C  IPRNT   -  print flag
C  EC      -  current energy
C  EOld    -  previous energy
C  X       -  current geometry
C  g       -  current gradient
C  gold    -  previous gradient
C  HINV    -  current estimate of inverse Hessian
C  V1      -  work array
C  V2      -   ditto
C  d       -  on exit contains new step
C  RMSG    -  on exit contains root mean square gradient
C  Cnvgd   -  Logical flag indicating convergence
C             on exit   .true.  - convergence achieved
C                       .false. - another step needed
C
C
      DIMENSION X(N),g(N),gold(N),HINV(N,N),d(N)
      Dimension V1(N),V2(N),vv(n,n)
      LOGICAL Invrt,Cnvgd
      COMMON /IO/ IOut,ICond
C
      PARAMETER (Zero=0.0d0,One=1.0d0,TollZero=1.0d-8)
      Parameter (NEG=0,NegReq=0)
C
C
C  is this the first cycle?
C  if invrt flag true, invert starting Hessian
C
      IF(NCycle.EQ.1) THEN
C
C  invert initial guess Hessian
C
        If(Invrt) Then
          CALL CpyVEC(N*N,HINV,VV)
          CALL INVMAT(VV,N,V1,V2,HINV,IErr)
          If(IErr.LT.0) call nerror(1,'optbfgs',
     $              '** ERROR INVERTING HESSIAN **',0,0)
        EndIf
c
      ELSE
C
C  update the inverse Hessian using BFGS formula
C
        DV = Zero
        Do I=1,N
        V1(I) = g(I)-gold(I)
        DV = DV + V1(I)*d(I)
        EndDO
C
C  if DV is negative, retention of positive definiteness is not
C  guaranteed. Print a warning, skip Hessian update and take
C  steepest descent step instead
C
        If(DV.LT.TollZero) Then
         WRITE(IOut,1200)
         CALL CpyVEC(N,g,d)
         GO TO 50
        EndIf
C
C  form HINV*V1
C
        CALL MATVEC(N,HINV,V1,V2)
        DT = One + SProd(N,V1,V2)/DV
C
C  do the update
C
        Do I=1,N
        Do J=1,I
        Tmp = DT*d(I)*d(J) - d(I)*V2(J) - V2(I)*d(J)
        HINV(I,J) = HINV(I,J) + Tmp/DV
        HINV(J,I) = HINV(I,J)
        EndDo
        EndDo
c
        If(IPRNT.GT.2) WRITE(IOut,1100)
cc
      ENDIF
C
C  work out the quasi Newton step
C
      CALL MATVEC(N,HINV,g,d)
C
 50   CONTINUE
C
C  D should be -D
C
      Do I=1,N
      d(I) = -d(I)
      EndDo
C
C  we have a new step in D
C  check the stepsize
C
      CALL ChkD(N,DMAX,IPRNT,d)
C
C  check for convergence and take the next step
C
      CALL CnvINT(N,      X,      d,      g,      EC,
     $            EOld,   IPRNT,  TolG,   TolD,   TolE,
     $            NEG,    NegReq, Cnvgd)
C
C  .............................................................
C    Standard Printout
C
C  form the RMS displacement
C
       RMSG = SProd(N,g,g)
       RMSG = SQRT(RMSG/DFloat(N))
       RMSD = SProd(N,d,d)
       RMSD = SQRT(RMSD/DFloat(N))
C
       If(IPRNT.EQ.1) WRITE(IOut,2000) NCycle,EC,RMSG,RMSD
       WRITE(ICond,2000) NCycle,EC,RMSG,RMSD       ! summary printout
       call f_lush(ICond)
C  ............................................................
C
      If(.NOT.Cnvgd) CALL AddVEC(N,X,d,X)
C
      RETURN
c
 1100 FORMAT(' Inverse Hessian updated using BFGS update')
 1200 FORMAT('**WARNING** Hereditary positive definiteness',
     $       ' endangered',/,
     $       ' Hessian Update Skipped this cycle',/,
     $       ' Taking Steepest Descent Step')
 2000 FORMAT(/,' ** Cycle ',I4,'  Energy ',F16.9,'   RMSG ',F8.5,
     $         '   RMSD ',F8.5,' **')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTCART(NCycle, NAtoms, AtSymb, XC,     NFix,
     $                   IFix,   Symflag,GROUP,  NDEG,   NTrans,
     $                   TRANS,  NEqATM, TSflag, MaxDiis,NDiis,
     $                   TolG,   TolD,   TolE,   DMAX,   mode,
     $                   IHess,  IProj,  IPRNT,  EC,     EOld,
     $                   GC,     HESS,   D,      XSTR,   GSTR,
     $                   VMODE,  Steep,  IMem,   Z,      RMSG,
     $                   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  MAIN DRIVING ROUTINE FOR GEOMETRY OPTIMIZATION IN
C  CARTESIAN COORDINATES
C
C
C  ARGUMENTS
C
C  NCycle   -  cycle number
C              i.e. number of times this routine has been called
C              in this job step
C  NAtoms   -  number of atoms
C  AtSymb   -  atomic symbols (used for nice printout)
C  XC       -  current cartesian coordinates
C              contains new coordinates on exit
C  NFix     -  number of fixed Cartesian coordinates
C  IFix     -  list of fixed/active coordinates
C               0 - coordinate active
C               1 - coordinate inactive (will be "fixed")
C  Symflag  -  Logical flag   .true.  - use symmetry
C                             .false. - do not use symmetry
C  GROUP    -  molecular point group
C  NDEG     -  number of internal degrees of freedom
C  NTrans   -  number of symmetry operations
C  ..................................................................
C  TRANS    -  symmetry operations as 3x3 transformation matrices
C  NEqATM   -  list of atomic equivalences under symmetry operations
C
C              if symmetry is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ..................................................................
C  TSflag   -  Logical flag   .true.  - TS search
C                             .false. - minimization
C  MaxDiis  -  maximum size of subspace for GDIIS
C               0 - do not use GDIIS
C               N - size specified by user (DO NOT set N too large)
C  NDiis    -  current size of GDIIS subspace
C              may be incremented/decremented on exit
C  TolG     -  convergence tolerance on maximum gradient component
C  TolD     -  convergence tolerance on maximum atomic displacement
C  TolE     -  convergence tolerance on energy
C              (if energy change from previous cycle is less than
C               TolE and TolG is satisfied, stop regardless of TolD;
C               this prevents unnecessary steps for floppy molecules)
C  DMAX     -  maximum allowed stepsize
C  mode     -  Hessian mode being followed during a TS search
C              if mode=0 mode following is switched off and
C              maximization will occur along the lowest mode
C  IHess    -  Hessian status flag
C               0 - new Hessian; symmetrize
C               1 - Hessian has been updated; use "as is"
C  IProj    -  Hessian projection flag
C               0 - do not project; leave Hessian "as is"
C               1 - project out translations
C               2 - project out translations and rotations
C  IPRNT    -  flag for controlling printout
C               0 - NO printout (except for error messages)
C               1 - summary and warning printout only
C               2 - "standard" printout
C               3 - heavier printout (including Hessian)
C               4 - very heavy (including some debug) printout
C  EC       -  current energy
C  EOld     -  previous energy
C  GC       -  current gradient
C  HESS     -  Hessian (Force Constant) matrix
C  D        -  displacement vector (step)
C  ................................................................
C  XSTR     -  Array for storing all geometries in GDIIS subspace
C  GSTR     -  Array for storing all gradients in GDIIS subspace
C
C              if GDIIS is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ................................................................
C  VMODE    -  Actual mode (vector) followed on previous cycle
C              on exit will contain mode followed on current cycle
C              determined by criterion of maximum overlap with
C              the previous mode
C  Steep    -  Logical flag   .true.  - take steepest descent step
C                             .false. - determine step as usual
C  ..................................................................
C  IMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  RMSG     -  on exit contains root mean square gradient
C  Cnvgd    -  Logical flag indicating convergence
C              on exit   .true.  - convergence achieved
C                        .false. - another step needed
C
C
      REAL*8 XC(3,NAtoms),GC(3,NAtoms),HESS(3*NAtoms,3*NAtoms),
     $       D(3*NAtoms),VMODE(3*NAtoms)
      REAL*8 XSTR(3*NAtoms,MaxDiis),GSTR(3*NAtoms,MaxDiis)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION IFix(3*NAtoms)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 GROUP
      LOGICAL Symflag,TSflag,Steep,Cnvgd,Diisflag
      COMMON /IO/ IOut,ICond                       ! standard/short output
C
      DIMENSION Z(IMem)
C
      PARAMETER (small=1.0d-8,jnk=0)
C
C
      NAT3 = 3*NAtoms
      NS = NAT3 - NFix        ! actual size of optimization space
C
C  Check options and set up defaults
C
      EigMIN = 0.0001d0
      EigMAX = 25.0d0
      TolDiis = 0.1d0
c
      NegReq = 0
      If(TSflag) NegReq = 1
C
C  allocate scratch pointers
C
      IU = 1
      IEig = IU + NAT3*NAT3
      IGC = IEig + NS
      IM1 = IGC + NAT3
      IM2 = IM1 + NAT3*NAT3
      IEnd = IM2 + MAX(NAT3*NAT3,MaxDiis+1)
      CALL MemCHK(IMem,IEnd,7,'OPTCART')
C
C
C  START THE OPTIMIZATION PROPER
C  ** NOTE:  HESSIAN UPDATING ALREADY DONE IN CALLING ROUTINE **
C
C  If there are no fixed coordinates, project translations
C  and rotations out of the Hessian
C
      IF(NFix.EQ.0) THEN
       IF(.NOT.Steep) THEN
        CALL CpyVEC(NAT3*NAT3,HESS,Z(IU))
        If(IProj.EQ.0) Then
         If(IPRNT.GT.1) WRITE(IOut,1000)
        Else
         CALL ProjTR(NAtoms, XC,     IProj,  IPRNT,  Z(IM1),
     $               Z(IM2),Z(IU))
        EndIf
       ENDIF
C
C  make a copy of the gradient
C
       CALL CpyVEC(NAT3,GC,Z(IGC))
      ELSE
C
C  remove fixed coordinates from parameter space
C
       CALL CpyVEC(NAT3*NAT3,HESS,Z(IM1))
       CALL RmFIX(NAT3,   IFix,   GC,     Z(IM1), Steep,
     $            NS,     Z(IGC), Z(IU))
      ENDIF
C
C  form the RMS gradient
C
      RMSG = SProd(NS,Z(IGC),Z(IGC))
      RMSG = SQRT(RMSG/DFloat(NS))
C
C  ..................................................................
C    STEEPEST DESCENT
C
      IF(Steep) THEN
       DO 10 I=1,NS
       D(I) = -Z(IGC+I-1)
 10    CONTINUE
       If(IPRNT.GT.1) WRITE(IOut,1100)
       GO TO 95
      ENDIF
C  ..................................................................
C
      IF(IPRNT.GT.3) THEN
       WRITE(IOut,1200)
       CALL PrntMAT(NS,NS,NS,Z(IU))
      ENDIF
C
C  Diagonalize
C
      CALL DIAGMAT(Z(IU),NS,Z(IM1),Z(IM2),Z(IEig),IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1300)
       CALL OptExit(9)
      ENDIF
C
C  check for and remove eigenvectors corresponding to zero
C  eigenvalues and possible symmetry-redundant modes
C
      thrsh = small                      ! zero-overlap threshold
      If(IHess.EQ.0) thrsh = DMIN1(small,RMSG*1.0d-4)
c
      CALL ChkHES(NS,     jnk,    IPRNT,  Z(IU),  Z(IEig),
     $            Symflag,Z(IGC), thrsh,  jnk,    NC)
C
C  on exit, NC is the number of Hessian modes that will be
C  used to form the next step.
C
C  check eigenvalue magnitudes are between allowed values
C
      CALL ChkMAG(NC,EigMIN,EigMAX,Z(IEig))
C
C  find the number of negative Hessian eigenvalues
C
      CALL FndNEG(NC,Z(IEig),NEG)
C
C  calculate the next step
C  this can be done using either GDIIS or the
C  standard EF algorithm
C
C  see if GDIIS can be used (minimization only)
C
      Diisflag = MaxDiis.GT.1.AND.RMSG.LT.TolDiis.AND.
     $           NEG.EQ.0.AND..NOT.TSflag
C
      IF(Diisflag) THEN
C
C  calculate scratch pointers for GDIIS
C
       IB = IEnd
       IErr = IB + (MaxDiis+1)**2
       IDS = IErr + NAT3*MaxDiis
       IEnd = IDS + MAX(NAT3,MaxDiis+1)
       CALL MemCHK(IMem,IEnd,7,'OPTCART')
c
       CALL GDIIS(NC,     NAT3,   MaxDiis,NDiis,  IPRNT,
     $            Z(IU),  Z(IEig),Z(IM1), XC,     GC,
     $            XSTR,   GSTR,   Z(IB),  Z(IErr),Z(IErr),
     $            Z(IM2), Z(IDS), D)
c
       If(NDiis.LE.1) Diisflag = .FALSE.
c
      ENDIF
C
      IF(.NOT.Diisflag) THEN
c
       CALL OPTEF(NCycle, NC,     NS,     jnk,    NEG,
     $            NegReq, mode,   IPRNT,  Z(IU),  Z(IEig),
     $            Z(IGC), Z(IM1), VMODE,  D,      IErr)
c
       If(IErr.NE.0) CALL OptExit(9)
c
      ENDIF
C
 95   CONTINUE
C
C  we have a new step in D
C  check the stepsize
C
      CALL ChkD(NS,DMAX,IPRNT,D)
C
C  restore fixed coordinates
C
      IF(NFix.GT.0) THEN
       CALL CpyVEC(NS,D,Z(IM1))
       CALL CpyVEC(NS,Z(IGC),Z(IM2))
       CALL RstFIX(NAtoms, 0,      NAT3,   IFix,   Z(IM1),
     $             Z(IM2), D,      Z(IGC))
      ENDIF
C
C  check for convergence and take the next step
C
      CALL CnvCART(NAtoms, AtSymb, XC,     D,      Z(IGC),
     $             EC,     EOld,   IPRNT,  TolG,   TolD,
     $             TolE,   NEG,    NegReq, Cnvgd)
C
C  .............................................................
C    Standard Printout
C
C  form the RMS displacement
C
       RMSD = SProd(NAT3,D,D)
       RMSD = SQRT(RMSD/DFloat(NC))
C
       If(IPRNT.EQ.1) WRITE(IOut,2000) NCycle,EC,RMSG,RMSD
       WRITE(ICond,2000) NCycle,EC,RMSG,RMSD       ! summary output
       call f_lush(ICond)
C  .............................................................
C
      IF(.NOT.Cnvgd) THEN
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $             small,  D)
       CALL AddVEC(NAT3,XC,D,XC)
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $             small,  XC)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,' Projection of Translations/Rotations Skipped',
     $         ' by Request')
 1100 FORMAT(/,' ** Taking Steepest Descent Step **')
 1200 FORMAT(/,'   Hessian Matrix after Projection')
 1300 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize Hessian Matrix')
 2000 FORMAT(/,' ** Cycle ',I4,'  Energy ',F16.9,'   RMSG ',F8.5,
     $         '   RMSD ',F8.5,' **')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTCON(NCycle, NAtoms, NDum,   NCons,  AtSymb,
     $                  XC,     ICTYP,  RCON,   IC,     NFix,
     $                  IFix,   MaxL,   IFunc,  Symflag,GROUP,
     $                  NDEG,   NTrans, TRANS,  NEqATM, TSflag,
     $                  TolG,   TolD,   TolE,   DMAX,   mode,
     $                  IProj,  IPRNT,  EC,     EOld,   GC,
     $                  HESS,   D,      VMODE,  NCount, EOS,
     $                  GOS,    IOS,    oscil,  IMem,   Z,
     $                  RMSG,   Cnvgd,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  MAIN DRIVING ROUTINE FOR CONSTRAINED GEOMETRY OPTIMIZATION
C  IN CARTESIAN COORDINATES
C
C
C  ARGUMENTS
C
C  NCycle   -  cycle number
C              i.e. number of times this routine has been called
C              in this job step
C  NAtoms   -  number of real atoms
C  NDum     -  number of dummy atoms
C  NCons    -  number of constraints
C  AtSymb   -  atomic symbols (used for nice printout)
C  XC       -  current cartesian coordinates
C              contains new coordinates on exit
C              ** XC also incorporates the Lagrange multipliers **
C  ICTYP    -  constraint type
C               1 - fixed distance
C               2 - fixed bond angle
C               3 - fixed dihedral angle
C  RCON     -  constraint value (in atomic units)
C  IC       -  atoms involved in constraint
C                IC1-IC2           distance constraint
C                IC1-IC2-IC3       angle constraint
C                IC1-IC2-IC3-IC4   dihedral constraint
C  ..................................................................
C  NFix     -  number of fixed Cartesian coordinates
C  IFix     -  list of fixed/active coordinates
C               0 - coordinate active
C               1 - coordinate inactive (will be "fixed")
C
C              if no atomic coordinates are fixed the above
C              array is not needed and can be a dummy argument
C  ..................................................................
C  MaxL     -  maximum number of real atoms involved in definition
C              of dummy atom position
C  IFunc    -  list for each dummy atom of the real atoms (if any)
C              defining its position
C
C              if there are no dummy atoms the above array is
C              not needed and can be a dummy argument
C  ...................................................................
C  Symflag  -  Logical flag   .true.  - use symmetry
C                             .false. - do not use symmetry
C  GROUP    -  molecular point group
C  NDEG     -  number of internal degrees of freedom
C  NTrans   -  number of symmetry operations
C  ..................................................................
C  TRANS    -  symmetry operations as 3x3 transformation matrices
C  NEqATM   -  list of atomic equivalences under symmetry operations
C
C              if symmetry is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ..................................................................
C  TSflag   -  Logical flag   .true.  - TS search
C                             .false. - minimization
C  TolG     -  convergence tolerance on maximum gradient component
C  TolD     -  convergence tolerance on maximum atomic displacement
C  TolE     -  convergence tolerance on energy
C              (if energy change from previous cycle is less than
C               TolE and TolG is satisfied, stop regardless of TolD;
C               this prevents unnecessary steps for floppy molecules)
C  DMAX     -  maximum allowed stepsize
C  mode     -  Hessian mode being followed during a TS search
C              if mode=0 mode following is switched off and
C              maximization will occur along the lowest mode
C  IProj    -  Hessian projection flag
C               0 - do not project; leave Hessian "as is"
C               1 - project out translations
C               2 - project out translations and rotations
C  IPRNT    -  flag for controlling printout
C               0 - NO printout (except for error messages)
C               1 - summary and warning printout only
C               2 - "standard" printout
C               3 - heavier printout (including Hessian)
C               4 - very heavy (including some debug) printout
C  EC       -  current energy
C  EOld     -  previous energy
C  GC       -  current gradient
C              ** GC also incorporates the Lagrange multipliers **
C  HESS     -  Hessian (Force Constant) matrix
C  D        -  displacement vector (step)
C              ** D also incorporates the Lagrange multipliers **
C  VMODE    -  Actual mode (vector) followed on previous cycle
C              on exit will contain mode followed on current cycle
C              determined by criterion of maximum overlap with
C              the previous mode
C              ** VMODE also incorporates the Lagrange multipliers **
C  NCount   -  counter for oscillation checking
C  EOS      -  storage for current and two previous energies
C  GOS      -  storage for current and two previous Lagrange gradients
C  IOS      -  number of times oscillatory behaviour observed
C  oscil    -  Logical flag indicating oscillation status
C           -  on entry  .true.  - oscillation on previous cycle
C                        .false. - no oscillation
C              on exit   .true.  - oscillation on this cycle
C                        .false. - no oscillation
C  ..................................................................
C  IMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  RMSG     -  on exit contains root mean square gradient
C  Cnvgd    -  Logical flag indicating convergence
C              on exit   .true.  - convergence achieved
C                        .false. - another step needed
C  IErr     -  error flag    0 - successful step
C                           -1 - something went wrong
C
C
      REAL*8 XC(3*(NAtoms+NDum)+NCons),GC(3,NAtoms),
     $       HESS(3*NAtoms,3*NAtoms),VMODE(3*NAtoms+NCons),
     $       D(3*(NAtoms+NDum)+NCons)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION ICTYP(NCons),IC(4,NCons),IFix(3*(NAtoms+NDum)+NCons),
     $          IFunc(NDum,MaxL)
      REAL*8 RCON(NCons),EOS(3),GOS(NCons,3)
      CHARACTER*8 AtSymb(NAtoms+NDum)
      CHARACTER*4 GROUP
      LOGICAL Symflag,TSflag,oscil,Cnvgd
      COMMON /IO/ IOut,ICond                       ! standard/short output
C
      DIMENSION Z(IMem)
C
      PARAMETER (One=1.0d0)
      PARAMETER (small=1.0d-8,MaxIOS=5)
C
C
      If(IPRNT.GT.1) WRITE(IOut,1000)
c
      NCNTR = NAtoms + NDum   ! total # of atomic centres
      NAT3 = 3*NAtoms
      NCTR3 = 3*NCNTR
      NDim = NCTR3 + NCons    ! maximum possible dimension
      NS = NDim - NFix        ! actual size of optimization space
C
C  Check options and set up defaults
C
      EigMIN = 0.0001d0
      EigMAX = 25.0d0
      scale = One
C
C  allocate scratch storage
C
      IHC = 1
      IU = IHC + NDim*NDim
      IEig = IU + NDim*NDim
      IGC = IEig + NS
      IGCC = IGC + NDim
      IM1 = IGCC + NDim*NCons
      IM2 = IM1 + NDim*NDim
      IM3 = IM2 + NDim*NDim
      IEnd = IM3 + NAT3*NAT3
      CALL MemCHK(IMem,IEnd,6,'OPTCON')
C
C
C  START THE CONSTRAINED OPTIMIZATION PROPER
C  ** NOTE:  HESSIAN UPDATING ALREADY DONE IN CALLING ROUTINE **
C
C  take a copy of the gradient and expand out to
C  include dummy atoms
C
      CALL CpyVEC(NAT3,GC,Z(IGC))
      If(NDum.GT.0) CALL ZeroIT(Z(IGC+NAT3),3*NDum)
C
C  modify the gradient to incorporate the constraints
C
      CALL ConGRAD(NCNTR,  NCons,  1,      ICTYP,  RCON,
     $             IC,     XC,     Z(IGC), Z(IGCC) )
C
C  Check for Oscillatory behaviour
C
      CALL ChkOSCIL(NCount, NDim,   NCons,  ICTYP,  IC,
     $              IPRNT,  TolG,   EC,     XC,     Z(IM1),
     $              D,      Z(IGC), EOS,    GOS,    IOS,
     $              oscil )
      If(oscil) RETURN
c
 10   CONTINUE
C
C  modify Hessian to incorporate constraint second derivatives
C  and dummy atoms
C
      CALL CpyMAT(NAT3,NAT3,NCTR3,HESS,Z(IU))
      CALL ConHESS(NCNTR,  NCons,  1,      ICTYP,  IC,
     $             XC,     Z(IGCC),Z(IM1), Z(IM2), Z(IU) )
C
C  make sure the modified gradient and Hessian conform fully
C  to the molecular symmetry
C
      CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $            small,  Z(IGC))
      CALL CpyMAT(NAT3,NCTR3,NAT3,Z(IU),Z(IHC))
      CALL SymHES(NAtoms, NTrans, NEqATM, TRANS,  Z(IHC),
     $            Z(IM1), Z(IM2), small,  Z(IM3))
C
C  If there are no fixed coordinates, project translations
C  and rotations out of the Hessian
C
      IF(NFix-3*NDum.EQ.0) THEN
       If(IProj.EQ.0) Then
        If(IPRNT.GT.1) WRITE(IOut,1100)
       Else
        CALL ProjTR(NAtoms, XC,     IProj,  IPRNT,  Z(IM1),
     $              Z(IM2), Z(IM3))
       EndIf
      ENDIF
      CALL CpyMAT(NAT3,NAT3,NCTR3,Z(IM3),Z(IU))
C
C  modify Hessian to incorporate constraint normals
C
      CALL LgngeHESS(NCTR3,NCons,Z(IU),Z(IGCC),Z(IHC))
c
      IF(NFix.GT.0) THEN
C
C  remove fixed coordinates from parameter space
C
       CALL CpyVEC(NDim,Z(IGC),Z(IM2))
       CALL CpyVEC(NDim*NDim,Z(IHC),Z(IM1))
       CALL RmFIX(NDim,   IFix,   Z(IM2), Z(IM1), .false.,
     $            NS,     Z(IGC), Z(IHC))
      ENDIF
c
      IF(IPRNT.GT.3) THEN
       If(NFix-3*NDum.EQ.0.AND.IProj.GT.0) Then
        WRITE(IOut,1200)
       Else
        WRITE(IOut,1300)
       EndIf
       CALL PrntMAT(NS,NS,NS,Z(IHC))
      ENDIF
C
C  Diagonalize
C
      CALL DIAGMAT(Z(IHC),NS,Z(IM1),Z(IM2),Z(IEig),IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1400)
       CALL OptExit(9)
      ENDIF
C
C  check for and remove eigenvectors corresponding to zero
C  eigenvalues and possible symmetry-redundant modes
C
      CALL ChkHES(NS,     NCons,  IPRNT,  Z(IHC), Z(IEig),
     $            Symflag,Z(IGC), small,  NCon,   NC)
C
C  on exit, NC is the number of Hessian modes that will be
C  used to form the next step; NCon is the number of surviving
C  constraints (some may have been removed by symmetry)
C
C  check eigenvalue magnitudes are between allowed values
C
      CALL ChkMAG(NC,EigMIN,EigMAX,Z(IEig))
C
C  find the number of negative Hessian eigenvalues
C
      CALL FndNEG(NC,Z(IEig),NEG)
C
C  make sure the constraint modes are the first vectors stored
C
      If(NEG.GT.NCon)
     $   CALL ConORD(NCTR3-NFix, NCon,  NS, NEG,    Z(IEig),
     $               Z(IHC), IPRNT, Z(IM1), Z(IM2), Z(IU),
     $               Z(IM3))
C
C  set NegReq
C
      NegReq = NCon
      If(TSflag) NegReq = NCon+1
C
C  form the RMS gradient
C
      RMSG = SProd(NS,Z(IGC),Z(IGC))
      RMSG = SQRT(RMSG/DFloat(NC))
C
C  calculate the next step
C  this is done using a modied version of the EF algorithm
C
      CALL OPTEF(NCycle, NC,     NS,     NCon,   NEG,
     $           NegReq, mode,   IPRNT,  Z(IHC), Z(IEig),
     $           Z(IGC), Z(IM1), VMODE,  D,      IErr )
c
      IF(IErr.EQ.-1) THEN
       CALL ChkSCAL(NAT3,IPRNT,Z(IEig),HESS,scale)
       If(scale.EQ.One) GO TO 95
       If(scale.NE.One) GO TO 10
      ENDIF
C
C  we have a new step in D
C  check the stepsize
C
      CALL ChkD(NS,DMAX,IPRNT,D)
C
C  restore fixed coordinates
C
      IF(NFix.GT.0) THEN
       CALL CpyVEC(NS,D,Z(IM1))
       CALL CpyVEC(NS,Z(IGC),Z(IM2))
       CALL RstFIX(NAtoms, NDum,   NDim,   IFix,   Z(IM1),
     $             Z(IM2), D,      Z(IGC) )
      ENDIF
C
C  check for convergence and take the next step
C
      CALL CnvCON(NAtoms, NCNTR,  NCon,   AtSymb, XC,
     $            D,      Z(IGC), EC,     EOld,   IPRNT,
     $            TolG,   TolD,   TolE,   NEG,    NegReq,
     $            Cnvgd )
C
C  .............................................................
C    Standard Printout
C
C  form the RMS displacement
C
       RMSD = SProd(NDim,D,D)
       RMSD = SQRT(RMSD/DFloat(NC))
C
       If(IPRNT.EQ.1) WRITE(IOut,2000) NCycle,EC,RMSG,RMSD
       WRITE(ICond,2000) NCycle,EC,RMSG,RMSD       ! summary printout
       call f_lush(ICond)
C  ............................................................
C
C  If persistent oscillation then quit
C
      IF(IOS.GT.MaxIOS.AND..NOT.Cnvgd) THEN
       WRITE(IOut,1500)
       Cnvgd = .TRUE.
      ENDIF
C
      IF(.NOT.Cnvgd) THEN
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $             small,  D)
       CALL AddVEC(NDim,XC,D,XC)
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $             small,  XC)
      ENDIF
C
C
 95   CONTINUE
      RETURN
c
 1000 FORMAT(/,'   Using Lagrange Multiplier Algorithm',/)
 1100 FORMAT(/,' Projection of Translations/Rotations Skipped',
     $         ' by Request')
 1200 FORMAT(/,'   Hessian Matrix after Projection')
 1300 FORMAT(/,'   Hessian Matrix (unprojected)')
 1400 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize Hessian Matrix')
 1500 FORMAT(/,' ** PERSISTENT OSCILLATORY BEHAVIOUR **',/,
     $     '    This is the best structure you are likely to get.',/,
     $     '    If unsatisfactory, try using the Penalty function',/,
     $     '    algorithm and then tidy up using Lagrange multipliers')
 2000 FORMAT(/,' ** Cycle ',I4,'  Energy ',F16.9,'   RMSG ',F8.5,
     $         '   RMSD ',F8.5,' **')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTCONINT(NCycle, NDEG,   NS,     NCons,  TSflag,
     $                     TolG,   TolD,   TolE,   DMAX,   mode,
     $                     IPRNT,  EC,     EOld,   ICHK,   RCHK,
     $                     XINT,   GC,     HESS,   D,      VMODE,
     $                     IMem,   Z,      RMSG,   Cnvgd,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  MAIN DRIVING ROUTINE FOR CONSTRAINED GEOMETRY OPTIMIZATION IN
C  DELOCALIZED INTERNAL COORDINATES FOR UNSATISFIED CONSTRAINTS
C
C
C  ARGUMENTS
C
C  NCycle  -  cycle number
C             i.e. number of times this routine has been called
C             in this job step
C  NDEG    -  total number of degrees of freedom
C  NS      -  number of unconstrained degrees of freedom
C  NCons   -  number of constraints
C  TSflag  -  Logical flag   .true.  - TS search
C                            .false. - minimization
C  TolG    -  convergence tolerance on maximum gradient component
C  TolD    -  convergence tolerance on maximum displacement
C  TolE    -  convergence tolerance on energy
C             (if energy change from previous cycle is less than
C              TolE and TolG is satisfied, stop regardless of TolD;
C              this prevents unnecessary steps for floppy molecules)
C  DMAX    -  maximum allowed stepsize
C  mode    -  Hessian mode being followed during a TS search
C             if mode=0 mode following is switched off and
C             maximization will occur along the lowest mode
C  IPRNT   -  flag for controlling printout
C              1 - summary and warning printout only
C              2 - "standard" printout
C              3 - slightly more printout (including gradient)
C              4 - heavier printout (including Hessian)
C              5 - heavier still (includes iterative printout)
C  EC      -  current energy
C  EOld    -  previous energy
C  ICHK    -  array indicating which constrains are satisfied
C               ICHK(I) = 0 - ith constraint is satisfied
C               ICHK(I) = 1 - ith constraint is not satisfied
C  RCHK    -  values for constraint differences
C             (current value minus desired value)
C  XINT    -  current internal coordinates
C             contains new coordinates on exit
C             ** XINT also incorporates the Lagrange multipliers **
C  GC      -  current gradient (internal)
C  HESS    -  Hessian (Force Constant) matrix
C  D       -  previous displacement vector (step)
C  VMODE   -  Actual mode (vector) followed on previous cycle
C             on exit will contain mode followed on current cycle
C             determined by criterion of maximum overlap with
C             the previous mode
C  ..................................................................
C  IMem    -  Amount of available scratch space
C  Z       -  Scratch array
C  ..................................................................
C  RMSG    -  on exit contains root mean square gradient
C  Cnvgd   -  Logical flag indicating convergence
C             on exit   .true.  - convergence achieved
C                       .false. - another step needed
C
C
      REAL*8 XINT(NDEG+NCons),GC(NDEG),D(NDEG+NCons),VMODE(NDEG+NCons),
     $       HESS(NDEG,NDEG)
      DIMENSION ICHK(NCons),RCHK(NCons)
      LOGICAL TSflag,Cnvgd
      COMMON /IO/ IOut,ICond                       ! standard/short output
C
      DIMENSION Z(IMem)
C
      PARAMETER (One=1.0d0)
      PARAMETER (small=1.0d-12,jnk=0)
C
C
      If(IPRNT.GT.1) WRITE(IOut,1000)
C
C  determine how many constraints are satisfied
C
      NCon = 0
      DO 10 I=1,NCons
      NCon = NCon + ICHK(I)
 10   CONTINUE
c
      NDim = NS + 2*NCon              ! full dimension of space
C
C  Check options and set up defaults
C
      EigMIN = 0.0001d0
      EigMAX = 25.0d0
      scale = One
C
C  allocate scratch pointers
C
      IU = 1
      IEig = IU + NDim*NDim
      IGC = IEig + NDim
      IM1 = IGC + NDim
      IM2 = IM1 + NDim*NDim
      IEnd = IM2 + NDim*NDim
      CALL MemCHK(IMem,IEnd,6,'OPTINT')
C
C
C  START THE OPTIMIZATION PROPER
C  ** NOTE:  HESSIAN UPDATING ALREADY DONE IN CALLING ROUTINE **
C
C  modify the gradient to incorporate the constraints
C
      CALL IntConGRAD(NS,     NCons,  NCon,   ICHK,   RCHK,
     $                XINT,   GC,     Z(IGC))
C
C  form the RMS gradient
C
      RMSG = SProd(NDim,Z(IGC),Z(IGC))
      RMSG = SQRT(RMSG/DFloat(NDim))
c
 50   CONTINUE
C
C  modify Hessian to incorporate constraint normals (unit vectors)
C
      CALL IntLgngeHESS(NS,NCons,NCon,ICHK,HESS,Z(IU))
c
      If(IPRNT.GT.3) THEN
       WRITE(IOut,1200)
       CALL PrntMAT(NDim,NDim,NDim,Z(IU))
      EndIf
C
C  Diagonalize
C
      CALL DIAGMAT(Z(IU),NDim,Z(IM1),Z(IM2),Z(IEig),IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1300)
       CALL OptExit(9)
      ENDIF
C
C  check for and remove eigenvectors corresponding to zero
C  eigenvalues and possible symmetry-redundant modes
C
      CALL ChkHES(NDim,   jnk,    IPRNT,  Z(IU),  Z(IEig),
     $            .true., Z(IGC), small,  jnk,    NC)
C
C  on exit, NC is the number of Hessian modes that will be
C  used to form the next step.
C
C  check eigenvalue magnitudes are between allowed values
C
      CALL ChkMAG(NC,EigMIN,EigMAX,Z(IEig))
C
C  find the number of negative Hessian eigenvalues
C
      CALL FndNEG(NC,Z(IEig),NEG)
C
C  set NegReq
C
      NegReq = NCon
      If(TSflag) NegReq = NCon+1
C
C  calculate the next step
C  this is done using a modified version of the EF algorithm
C
      CALL OPTEF(NCycle, NC,     NDim,   NCon,   NEG,
     $           NegReq, mode,   IPRNT,  Z(IU),  Z(IEig),
     $           Z(IGC), Z(IM1), VMODE,  D,      IErr)
c
      IF(IErr.EQ.-1) THEN
       CALL ChkSCAL(NS+NCon,IPRNT,Z(IEig),HESS,scale)
       If(scale.EQ.One) GO TO 95
       If(scale.NE.One) GO TO 50
      ENDIF
C
C  we have a new step in D
C  check the stepsize
C
      CALL ChkD(NS+NCon,DMAX,IPRNT,D)
C
C  check for convergence and take the next step
C
      CALL IntCnvCON(NS,     NCons,  NCon,   ICHK,   XINT,
     $               D,      Z(IGC), EC,     EOld,   IPRNT,
     $               TolG,   TolD,   TolE,   NEG,    NegReq,
     $               Cnvgd)
C
C  .............................................................
C    Standard Printout
C
C  form the RMS displacement
C
       RMSD = SProd(NS+NCon,D,D)
       RMSD = SQRT(RMSD/DFloat(NC))
C
       If(IPRNT.EQ.1) WRITE(IOut,2000) NCycle,EC,RMSG,RMSD
       WRITE(ICond,2000) NCycle,EC,RMSG,RMSD       ! summary printout
       call f_lush(ICond)
C  ............................................................
C
C  now expand D to NDEG+NCons
C
      CALL CpyVEC(NDim,D,Z(IM1))
      CALL ExpandV(NS,NCons,NCon,ICHK,Z(IM1),D)
c
      If(.NOT.Cnvgd) CALL AddVEC(NDEG+NCons,XINT,D,XINT)
C
C
 95   CONTINUE
      RETURN
c
 1000 FORMAT(/,'   Using Lagrange Multiplier Algorithm',/)
 1200 FORMAT(/,'   Hessian Matrix')
 1300 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize Hessian Matrix')
 2000 FORMAT(/,' ** Cycle ',I4,'  Energy ',F16.9,'   RMSG ',F8.5,
     $         '   RMSD ',F8.5,' **')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTINT(NCycle, NDEG,   TSflag, MaxDiis,NDiis,
     $                  TolG,   TolD,   TolE,   DMAX,   IHess,
     $                  IUpDat, mode,   IPRNT,  EC,     EOld,
     $                  XINT,   GC,     GOld,   HESS,   D,
     $                  XSTR,   GSTR,   VMODE,  Steep,  IMem,
     $                  Z,      RMSG,   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  MAIN DRIVING ROUTINE FOR GEOMETRY OPTIMIZATION IN INTERNAL
C  COORDINATES  (both delocalized internals and Z-Matrix)
C
C
C  ARGUMENTS
C
C  NCycle   -  cycle number
C              i.e. number of times this routine has been called
C              in this job step
C  NDEG     -  number of variables to be optimized
C  TSflag   -  Logical flag   .true.  - TS search
C                             .false. - minimization
C  MaxDiis  -  maximum size of subspace for GDIIS
C               0 - do not use GDIIS
C               N - size specified by user (DO NOT set N too large)
C  NDiis    -  current size of GDIIS subspace
C              may be incremented/decremented on exit
C  TolG     -  convergence tolerance on maximum gradient component
C  TolD     -  convergence tolerance on maximum displacement
C  TolE     -  convergence tolerance on energy
C              (if energy change from previous cycle is less than
C               TolE and TolG is satisfied, stop regardless of TolD;
C               this prevents unnecessary steps for floppy molecules)
C  DMAX     -  maximum allowed stepsize
C  IHess    -  Hessian status flag
C               0 - have "exact" or initial Hessian   use as is
C               1 - have Hessian from previous step   need to update
C  IUpDat   -  Hessian update flag
C               0 - do not update the Hessian (!?)
C              -1 - use default update
C               1 - Powell update (default for TS search)
C               2 - BFGS update (default for minimization)
C               3 - BFGS with safeguards to ensure retention of
C                   positive definiteness (default for GDIIS)
C  mode     -  Hessian mode being followed during a TS search
C              if mode=0 mode following is switched off and
C              maximization will occur along the lowest mode
C  IPRNT    -  flag for controlling printout
C               0 - NO printout (except for error messages)
C               1 - summary and warning printout only
C               2 - "standard" printout
C               3 - heavier printout (including Hessian)
C               4 - very heavy (including some debug) printout
C  EC       -  current energy
C  EOld     -  previous energy
C  XINT     -  current internal coordinates
C              (will be updated on exit)
C  GC       -  current gradient (internal)
C  GOld     -  previous gradient (internal)
C  HESS     -  Hessian (Force Constant) matrix
C              The status of the Hessian is controlled by IHess
C  D        -  previous displacement vector (step)
C  ................................................................
C  XSTR     -  Array for storing all geometries in GDIIS subspace
C  GSTR     -  Array for storing all gradients in GDIIS subspace
C
C              if GDIIS is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ................................................................
C  VMODE    -  Actual mode (vector) followed on previous cycle
C              on exit will contain mode followed on current cycle
C              determined by criterion of maximum overlap with
C              the previous mode
C  Steep    -  Logical flag   .true.  - take steepest descent step
C                             .false. - determine step as usual
C  ..................................................................
C  IMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ..................................................................
C  RMSG    -  on exit contains root mean square gradient
C  Cnvgd    -  Logical flag indicating convergence
C              on exit   .true.  - convergence achieved
C                        .false. - another step needed
C
C
      REAL*8 XINT(NDEG),GC(NDEG),GOld(NDEG),VMODE(NDEG),
     $       HESS(NDEG,NDEG),D(NDEG)
      REAL*8 XSTR(NDEG,MaxDiis),GSTR(NDEG,MaxDiis)
      LOGICAL TSflag,Steep,Cnvgd,Diisflag
      COMMON /IO/ IOut,ICond                       ! standard/short output
C
      DIMENSION Z(IMem)
C
      PARAMETER (small=1.0d-8,jnk=0)
C
C
C  Check options and set up defaults
C
      EigMIN = 0.0001d0
      EigMAX = 25.0d0
      TolDiis = 0.1d0
c
      NegReq = 0
      If(TSflag) NegReq = 1
C
C  form the RMS gradient
C
      RMSG = SProd(NDEG,GC,GC)
      RMSG = SQRT(RMSG/DFloat(NDEG))
C
C  ..................................................................
C    STEEPEST DESCENT
C
      IF(Steep) THEN
       DO 10 I=1,NDEG
       D(I) = -GC(I)
 10    CONTINUE
       If(IPRNT.GT.1) WRITE(IOut,1000)
       NC = NDEG      !  for proper calculation of RMSD
       GO TO 95
      ENDIF
C  ..................................................................
C
C  allocate scratch pointers
C
      IU = 1
      IEig = IU + NDEG*NDEG
      IM1 = IEig + NDEG
      IM2 = IM1 + NDEG*NDEG
      IEnd = IM2 + MAX(NDEG*NDEG,MaxDiis+1)
      CALL MemCHK(IMem,IEnd,6,'OPTINT')
C
C
C  START THE OPTIMIZATION PROPER
C  First update the Hessian if necessary
C
      If(IHess.EQ.1) CALL UpdHES(NDEG,   IUpDat, IPRNT,  D,      GC,
     $                           GOld,   Z(IM1), Z(IM2), HESS)
c
      CALL CpyVEC(NDEG*NDEG,HESS,Z(IU))
C
      IF(IPRNT.GT.3) THEN
       WRITE(IOut,1100)
       CALL PrntMAT(NDEG,NDEG,NDEG,Z(IU))
      ENDIF
C
C  Diagonalize
C
      CALL DIAGMAT(Z(IU),NDEG,Z(IM1),Z(IM2),Z(IEig),IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1200)
       CALL OptExit(9)
      ENDIF
C
C  check for and remove eigenvectors corresponding to zero
C  eigenvalues and possible symmetry-redundant modes
C
      thrsh = small                      ! zero-overlap threshold
      If(IHess.EQ.0) thrsh = DMIN1(small,RMSG*1.0d-4)
c
      CALL ChkHES(NDEG,   jnk,    IPRNT,  Z(IU),  Z(IEig),
     $            .true., GC,     thrsh,  jnk,    NC)
C
C  on exit, NC is the number of Hessian modes that will be
C  used to form the next step.
C
C  check eigenvalue magnitudes are between allowed values
C
      CALL ChkMAG(NC,EigMIN,EigMAX,Z(IEig))
C
C  find the number of negative Hessian eigenvalues
C
      CALL FndNEG(NC,Z(IEig),NEG)
C
C  calculate the next step
C  this can be done using either GDIIS or the
C  standard EF algorithm
C
C  see if GDIIS can be used (minimization only)
C
      Diisflag = MaxDiis.GT.1.AND.RMSG.LT.TolDiis.AND.
     $           NEG.EQ.0.AND..NOT.TSflag
C
      IF(Diisflag) THEN
C
C  calculate scratch pointers for GDIIS
C
       IB = IEnd
       IErr = IB + (MaxDiis+1)**2
       IDS = IErr + NDEG*MaxDiis + 1
       IEnd = IDS + MAX(NDEG,MaxDiis+1)
       CALL MemCHK(IMem,IEnd,6,'OPTINT')
c
       CALL GDIIS(NC,     NDEG,   MaxDiis,NDiis,  IPRNT,
     $            Z(IU),  Z(IEig),Z(IM1), XINT,   GC,
     $            XSTR,   GSTR,   Z(IB),  Z(IErr),Z(IErr),
     $            Z(IM2), Z(IDS), D)
c
       If(NDiis.LE.1) Diisflag = .FALSE.
c
      ENDIF
C
      IF(.NOT.Diisflag) THEN
c
       CALL OPTEF(NCycle, NC,     NDEG,   jnk,    NEG,
     $            NegReq, mode,   IPRNT,  Z(IU),  Z(IEig),
     $            GC,     Z(IM1), VMODE,  D,      IErr)
c
       If(IErr.NE.0) CALL OptExit(9)
c
      ENDIF
C
 95   CONTINUE
C
C  we have a new step in D
C  check the stepsize
C
      CALL ChkD(NDEG,DMAX,IPRNT,D)
C
C  check for convergence and take the next step
C
      CALL CnvINT(NDEG,   XINT,   D,      GC,     EC,
     $            EOld,   IPRNT,  TolG,   TolD,   TolE,
     $            NEG,    NegReq, Cnvgd)
C
C  .............................................................
C    Standard Printout
C
C  form the RMS displacement
C
       RMSD = SProd(NDEG,D,D)
       RMSD = SQRT(RMSD/DFloat(NC))
C
       If(IPRNT.EQ.1) WRITE(IOut,2000) NCycle,EC,RMSG,RMSD
       WRITE(ICond,2000) NCycle,EC,RMSG,RMSD       ! summary printout
       call f_lush(ICond)
C  ............................................................
C
      If(.NOT.Cnvgd) CALL AddVEC(NDEG,XINT,D,XINT)
C
      RETURN
c
 1000 FORMAT(/,' ** Taking Steepest Descent Step **')
 1100 FORMAT(/,'   Hessian Matrix')
 1200 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize Hessian Matrix')
 2000 FORMAT(/,' ** Cycle ',I4,'  Energy ',F16.9,'   RMSG ',F8.5,
     $         '   RMSD ',F8.5,' **')
c
      END
c  =======================================================================
c
      SUBROUTINE OPTPEN(NCycle, NAtoms, NDum,   NCons,  AtSymb,
     $                  XC,     ICTYP,  RCON,   IC,     NFix,
     $                  IFix,   MaxL,   IFunc,  Symflag,GROUP,
     $                  NDEG,   NTrans, TRANS,  NEqATM, MaxDiis,
     $                  NDiis,  TolG,   TolD,   TolE,   DMAX,
     $                  IProj,  IPRNT,  EC,     EOld,   GC,
     $                  HESS,   D,      XSTR,   GSTR,   IMem,
     $                  Z,      RMSG,   Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  MAIN DRIVING ROUTINE FOR CONSTRAINED MINIMIZATION IN
C  CARTESIAN COORDINATES USING PENALTY FUNCTIONS
C
C
C  ARGUMENTS
C
C  NCycle   -  cycle number
C              i.e. number of times this routine has been called
C              in this job step
C  NAtoms   -  number of real atoms
C  NDum     -  number of dummy atoms
C  NCons    -  number of constraints
C  AtSymb   -  atomic symbols (used for nice printout)
C  XC       -  current cartesian coordinates
C              contains new coordinates on exit
C  ICTYP    -  constraint type
C               1 - fixed distance
C               2 - fixed bond angle
C               3 - fixed dihedral angle
C  RCON     -  constraint value (in atomic units)
C  IC       -  atoms involved in constraint
C                IC1-IC2           distance constraint
C                IC1-IC2-IC3       angle constraint
C                IC1-IC2-IC3-IC4   dihedral constraint
C  ..................................................................
C  NFix     -  number of fixed Cartesian coordinates
C  IFix     -  list of fixed/active coordinates
C               0 - coordinate active
C               1 - coordinate inactive (will be "fixed")
C
C              if no atomic coordinates are fixed the above
C              array is not needed and can be a dummy argument
C  ..................................................................
C  MaxL     -  maximum number of real atoms involved in definition
C              of dummy atom position
C  IFunc    -  list for each dummy atom of the real atoms (if any)
C              defining its position
C
C              if there are no dummy atoms the above array is
C              not needed and can be a dummy argument
C  ...................................................................
C  Symflag  -  Logical flag   .true.  - use symmetry
C                             .false. - do not use symmetry
C  GROUP    -  molecular point group
C  NDEG     -  number of internal degrees of freedom
C  NTrans   -  number of symmetry operations
C  ..................................................................
C  TRANS    -  symmetry operations as 3x3 transformation matrices
C  NEqATM   -  list of atomic equivalences under symmetry operations
C
C              if symmetry is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ..................................................................
C  MaxDiis  -  maximum size of subspace for GDIIS
C               0 - do not use GDIIS
C               N - size specified by user (DO NOT set N too large)
C  NDiis    -  current size of GDIIS subspace
C              may be incremented/decremented on exit
C  TolG     -  convergence tolerance on maximum gradient component
C  TolD     -  convergence tolerance on maximum atomic displacement
C  TolE     -  convergence tolerance on energy
C              (if energy change from previous cycle is less than
C               TolE and TolG is satisfied, stop regardless of TolD;
C               this prevents unnecessary steps for floppy molecules)
C  DMAX     -  maximum allowed stepsize
C  IProj    -  Hessian projection flag
C               0 - do not project; leave Hessian "as is"
C               1 - project out translations
C               2 - project out translations and rotations
C  IPRNT    -  flag for controlling printout
C               0 - NO printout (except for error messages)
C               1 - summary and warning printout only
C               2 - "standard" printout
C               3 - heavier printout (including Hessian)
C               4 - very heavy (including some debug) printout
C  EC       -  current energy
C  EOld     -  previous energy
C  GC       -  current gradient
C  HESS     -  Hessian (Force Constant) matrix
C  D        -  displacement vector (step)
C  ................................................................
C  XSTR     -  Array for storing all geometries in GDIIS subspace
C  GSTR     -  Array for storing all gradients in GDIIS subspace
C
C              if GDIIS is not being used neither of the above
C              arrays are needed and can be dummy arguments
C  ................................................................
C  IMem     -  Amount of available scratch space
C  Z        -  Scratch array
C  ................................................................
C  RMSG     -  on exit contains root mean square gradient
C  Cnvgd    -  Logical flag indicating convergence
C              on exit   .true.  - convergence achieved
C                        .false. - another step needed
C
C
      REAL*8 XC(3*(NAtoms+NDum)+NCons),GC(3,NAtoms),
     $       HESS(3*NAtoms,3*NAtoms),D(3*(NAtoms+NDum))
      REAL*8 XSTR(3*NAtoms,MaxDiis),GSTR(3*NAtoms,MaxDiis)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION IFix(3*(NAtoms+NDum)),IFunc(NDum,MaxL)
      DIMENSION ICTYP(NCons),RCON(NCons),IC(4,NCons)
      CHARACTER*8 AtSymb(NAtoms+NDum)
      CHARACTER*4 GROUP
      LOGICAL Symflag,Cnvgd,Diisflag
      COMMON /IO/ IOut,ICond                       ! standard/short output
C
      DIMENSION Z(IMem)
C
      PARAMETER (small=1.0d-8,jnk=0)
C
C
      If(IPRNT.GT.1) WRITE(IOut,1000)
c
      NCNTR = NAtoms+NDum      ! total # of atomic centres
      NAT3 = 3*NAtoms
      NCTR3 = 3*NCNTR
      NS = NCTR3 - NFix        ! actual size of optimization space
C
C  Check options and set up defaults
C
      EigMIN = 0.0001d0
      EigMAX = 25.0d0
      TolDiis = 0.1d0
c
      NegReq = 0
      mode = 0
C
C  allocate scratch pointers
C
      IHC = 1
      IU = IHC + NAT3*NAT3
      IEig = IU + NCTR3*NCTR3
      IGC = IEig + NS
      IGCC = IGC + NCTR3
      IM1 = IGCC + NCTR3*NCons
      IM2 = IM1 + NAT3*NAT3
      IEnd = IM2 + MAX(NAT3*NAT3,MaxDiis+1)
      CALL MemCHK(IMem,IEnd,6,'OPTPEN')
C
C
C  START THE CONSTRAINED OPTIMIZATION PROPER
C  ** NOTE:  HESSIAN UPDATING ALREADY DONE IN CALLING ROUTINE **
C
C  take a copy of the gradient and expand out to
C  include dummy atoms
C
      CALL CpyVEC(NAT3,GC,Z(IGC))
      If(NDum.GT.0) CALL ZeroIT(Z(IGC+NAT3),3*NDum)
C
C  modify the gradient to incorporate the constraints
C
      CALL ConGRAD(NCNTR,  NCons,  2,      ICTYP,  RCON,
     $             IC,     XC,     Z(IGC), Z(IGCC) )
C
C  modify Hessian to incorporate constraint second derivatives
C  and dummy atoms
C
      CALL CpyMAT(NAT3,NAT3,NCTR3,HESS,Z(IU))
      CALL ConHESS(NCNTR,  NCons,  2,      ICTYP,  IC,
     $             XC,     Z(IGCC),Z(IM1), Z(IM2), Z(IU) )
C
C  make sure the modified gradient and Hessian conform fully
C  to the molecular symmetry
C
      CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $            small,  Z(IGC))
      CALL CpyMAT(NAT3,NCTR3,NAT3,Z(IU),Z(IHC))
      CALL SymHES(NAtoms, NTrans, NEqATM, TRANS,  Z(IHC),
     $            Z(IM1), Z(IM2), small,  Z(IU))
C
C  If there are no fixed coordinates, project translations
C  and rotations out of the Hessian
C
      IF(NFix-3*NDum.EQ.0) THEN
       If(IProj.EQ.0) Then
        If(IPRNT.GT.1) WRITE(IOut,1100)
       Else
        CALL ProjTR(NAtoms, XC,     IProj,  IPRNT,  Z(IM1),
     $              Z(IM2), Z(IU))
       EndIF
      ENDIF
c
      IF(NFix.GT.0) THEN
C
C  remove fixed coordinates from parameter space
C
       CALL CpyVEC(NAT3,Z(IGC),Z(IM2))
       CALL CpyVEC(NAT3*NAT3,Z(IU),Z(IM1))
       CALL RmFIX(NAT3,   IFix,   Z(IM2), Z(IM1), .false.,
     $            NS,     Z(IGC), Z(IU))
      ENDIF
C
C  form the RMS gradient
C
      RMSG = SProd(NAT3,Z(IGC),Z(IGC))
      RMSG = SQRT(RMSG/DFloat(NS))
c
      IF(IPRNT.GT.3) THEN
       If(NFix.EQ.0.AND.IProj.GT.0) Then
        WRITE(IOut,1200)
       Else
        WRITE(IOut,1300)
       EndIf
       CALL PrntMAT(NS,NS,NS,Z(IU))
      ENDIF
C
C  Diagonalize
C
      CALL DIAGMAT(Z(IU),NS,Z(IM1),Z(IM2),Z(IEig),IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(IOut,1400)
       CALL OptExit(9)
      ENDIF
C
C  check for and remove eigenvectors corresponding to zero
C  eigenvalues and possible symmetry-redundant modes
C
      CALL ChkHES(NS,     jnk,    IPRNT,  Z(IU),  Z(IEig),
     $            Symflag,Z(IGC), small,  jnk,    NC)
C
C  on exit, NC is the number of Hessian modes that will be
C  used to form the next step.
C
C  check eigenvalue magnitudes are between allowed values
C
      CALL ChkMAG(NC,EigMIN,EigMAX,Z(IEig))
C
C  find the number of negative Hessian eigenvalues
C
      CALL FndNEG(NC,Z(IEig),NEG)
C
C  calculate the next step
C  this can be done using either GDIIS or the
C  standard EF algorithm
C
C  see if GDIIS can be used (minimization only)
C
      Diisflag = MaxDiis.GT.1.AND.RMSG.LT.TolDiis.AND.
     $           NEG.EQ.0
C
      IF(Diisflag) THEN
C
C  calculate scratch pointers for GDIIS
C
       IB = IEnd
       IErr = IB + (MaxDiis+1)**2
       IDS = IErr + NAT3*MaxDiis
       IEnd = IDS + MAX(NAT3,MaxDiis+1)
       CALL MemCHK(IMem,IEnd,6,'OPTPEN')
c
       CALL GDIIS(NC,     NAT3,   MaxDiis,NDiis,  IPRNT,
     $            Z(IU),  Z(IEig),Z(IM1), XC,     Z(IGC),
     $            XSTR,   GSTR,   Z(IB),  Z(IErr),Z(IErr),
     $            Z(IM2), Z(IDS), D)
c
       If(NDiis.LE.1) Diisflag = .FALSE.
c
      ENDIF
C
      IF(.NOT.Diisflag) THEN
c
       CALL OPTEF(NCycle, NC,     NS,     jnk,    NEG,
     $            NegReq, mode,   IPRNT,  Z(IU),  Z(IEig),
     $            Z(IGC), Z(IM1), jnk,    D,      IErr)
c
      ENDIF
C
C  we have a new step in D
C  check the stepsize
C
      CALL ChkD(NS,DMAX,IPRNT,D)
C
C  restore fixed coordinates
C
      IF(NFix.GT.0) THEN
       CALL CpyVEC(NS,D,Z(IM1))
       CALL CpyVEC(NS,Z(IGC),Z(IM2))
       CALL RstFIX(NAtoms, NDum,   NCTR3,  IFix,   Z(IM1),
     $             Z(IM2), D,      Z(IGC))
      ENDIF
C
C  check for convergence and take the next step
C
      CALL CnvCART(NAtoms, AtSymb, XC,     D,      Z(IGC),
     $             EC,     EOld,   IPRNT,  TolG,   TolD,
     $             TolE,   NEG,    NegReq, Cnvgd)
C
C  .............................................................
C    Standard Printout
C
C  form the RMS displacement
C
       RMSD = SProd(NAT3,D,D)
       RMSD = SQRT(RMSD/DFloat(NC))
C
       If(IPRNT.EQ.1) WRITE(IOut,2000) NCycle,EC,RMSG,RMSD
       WRITE(ICond,2000) NCycle,EC,RMSG,RMSD       ! summary printout
       call f_lush(ICond)
C  ............................................................
C
      IF(.NOT.Cnvgd) THEN
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $             small,  D)
       CALL AddVEC(NCTR3,XC,D,XC)
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  Z(IM1),
     $             small,  XC)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,'   Using Penalty Function Algorithm',/)
 1100 FORMAT(/,' Projection of Translations/Rotations Skipped',
     $         ' by Request')
 1200 FORMAT(/,'   Hessian Matrix after Projection')
 1300 FORMAT(/,'   Hessian Matrix (unprojected)')
 1400 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize Hessian Matrix')
 2000 FORMAT(/,' ** Cycle ',I4,'  Energy ',F16.9,'   RMSG ',F8.5,
     $         '   RMSD ',F8.5,' **')
c
      END
