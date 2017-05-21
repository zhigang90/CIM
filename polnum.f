c ==================================================================
c  NUMERICAL POLARIZABILITY DERIVATIVES      JB   May 2001
c ==================================================================
c
      subroutine prepol(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval
c
c  reads the POLARIZE line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=5)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character cdum*20
      character*256 jobname
      logical isthere
c
      parameter (IUnit=1)
c
      data options/'nump','fiel','dipd','pold','prin'/
      data ioptyp/0,11,0,0,1/
c
c -- deal with jobname header
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
      call setival('isumscf',1)      ! switch off SCF summary print
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c ...........................................................
c -- check for existence of <polchk> file
c -- if exists, ensure best integral threshold
c -- (skip if semiempirical)         JB 9 April 99
c
      inquire(file=jobname(1:lenJ)//'.polchk',exist=isthere)
      If(isthere) Then
        call tstchval('semi',iprnt)
        if(iprnt.eq.0) then
          call getrval('ithr',thresh)
          call setrval('ith1',thresh)
        endif
      EndIf
c ...........................................................
c
      if(ifound(2).eq.1) then
        field = ropval(1,2)
      else
        field = 0.005d0         ! default field (au)
      endif
c
c -- dipole derivatives?
      if(ifound(3).eq.1) then
        itype = 1               ! dipole derivatives
      else
        itype = 0               ! polarizability only
      endif
c
c -- polarizability derivatives?
      if(ifound(4).eq.1) then
        itype = itype+10        ! polarizability derivatives
      endif
c
c -- print flag
      if(ifound(5).eq.1) then
        IPRNT = iopval(1,5)
      else
        IPRNT = 2               ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call WrPOL(IUnit,field,itype)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c  =======================================================================
c
      SUBROUTINE POLNUM(NMem,Z,Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** NUMERICAL POLARIZABILITY DERIVATIVES PROGRAM **
C
C  This Program determines the Cartesian polarizability derivatives
C  by finite-difference on the (analytic) gradient in appropriate
C  electric fields
C
C  At the same time we can also obtain the polarizabilities and
C  the dipole moment derivatives
C
C  polarizability  =  dipole moment differentiated once w.r.t.field
C  dipole derivatives = gradient differentiated once w.r.t.field
C  polarizability derivatives = gradient differentiated twice w.r.t.field
C
C  FILES
C
C  Data is transferred to and from POLNUM via the following
C  files:
C
C  <sym>      -  number of atoms, symmetry data
C  <coord>    -  current geometry (Cartesian coordinates)
C  <control>  -  program options/energy
C  <grad>     -  current gradient
C  <deriv>    -  final Cartesian polarizability/dipole derivatives
C  <polchk>   -  information pertinent to next finite difference
C                step  (produced by POLNUM itself)
C
C
C  **IMPORTANT**
C    It is assumed that the dipole and gradient (if required) at the
C    initial geometry are already available on the <control> file on
C    first entry to this program module
C  ............................................................
C
      DIMENSION Z(NMem)
      Character jobname*256,cdum*20
      Logical Cnvgd
C
      PARAMETER (MSymOP=120)     !  maximum number of symmetry operations
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    total number of atomic centres, including dummies
C    total number of molecules
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtom,dum,cdum)
      call rdcntrl(IUnit,5,'$nmol',1,NMol,dum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Read from the <sym> file
C    number of atoms
C    number of symmetry operations
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,dum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C --------------------------------------------------------------
C
C  Now get the memory
C
      NAt3 = 3*NAtoms
      IMem = 24*NAt3 + 5*NAtom + NAtoms + NAtoms*MSymOP
     $               + 9*MSymOP + NMol+1 + 32
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,6,'POLNUM')
      CALL ZeroIT(Z,IMem)      ! clear the memory
C
C  Allocate memory pointers
C
      IXC = iptr                    !  current geometry
      IML = IXC + 3*NAtom           !  pointer to molecules
      IXCG = IML +  NMol+1          !  atomic charges
      IXCM = IXCG + NAtom           !  atomic masses
      IUQ = IXCM + NAtom            !  list of symmetry unique atoms
      ITN = IUQ + NAtoms            !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*MSymOP          !  list of atomic equivalences
      IDP = INQ + NAtoms*MSymOP     !  dipole moment
      IG0 = IDP + 3                 !  reference (initial) gradient
      IDM = IG0 + NAt3              !  various finite-difference dipoles
      IGM = IDM + 6*3               !  various finite-difference gradients
      IPO = IGM + 12*NAt3           !  polarizability
      IPD = IPO + 9                 !  polarizability derivatives
      IDD = IPD + 6*NAt3            !  dipole moment derivatives
      IGC = IDD + 3*NAt3            !  current gradient
      IWK = IGC + NAt3              !  work array
      IEnd = IWK + NAt3 - 1
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,6,'POLNUM')
C
C  ----------------------------------------------------------------------
C
      CALL POLMAIN(NAtoms,  NAtom,   MSymOP,  Z(IXC),  NMol,
     $             Z(IML),  Z(IXCG), Z(IXCM), Z(IUQ),  Z(IDP),
     $             Z(IG0),  Z(IDM),  Z(IGM),  Z(IPO),  Z(IPD),
     $             Z(IDD),  Z(IGC),  Z(ITN),  Z(INQ),  Z(IWK),
     $             Cnvgd)
C
C  ----------------------------------------------------------------------
C
C  free memory used in this routine
C
cc      CALL ffree(Z(iptr))
C
C  Exit procedure
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE POLMAIN(NAtoms, NAtom,  MSymOP, XC,     NMol,
     $                   IMOL,   XCharg, XMass,  IUNQ,   Dip,
     $                   G0,     DipM,   GM,     Pol,    PolD,
     $                   DipD,   GC,     TRANS,  NEqATM, V,
     $                   Finish)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for numerical polarizability derivatives
C  POLMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(3,NAtom),XCharg(NAtom),XMass(NAtom),V(3*NAtoms),
     $       DipM(3,6),G0(3*NAtoms),GM(3*NAtoms,12),GC(3*NAtoms)
      INTEGER IMOL(NMol+1),IUNQ(NAtoms)
      DIMENSION TRANS(3,3,MSymOP),NEqATM(NAtoms,MSymOP),RM(3,3)
      DIMENSION EField(3),EField0(3),Dip0(3),Dip(3),DipD(3*NAtoms,3),
     $          Pol(3,3),PolD(3*NAtoms,6)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtom)
c ..................................................
      CHARACTER GROUP*4,cdum*20,jobname*256
      LOGICAL Symflag,Finish,isthere,DipolD,PolarD,Dipol
      LOGICAL Eflag,SCFFlag,Dflag,Sflag
C
      PARAMETER (Zero=0.0d0)
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
c
      DATA Eflag/.False./, SCFflag/.False./,
     $     Dflag/.False./, Sflag/.False./
c
      Common /job/jobname,lenJ
C
C
C  initialize
C
      IOut=ioutfil('iout')
      NAt3 = 3*NAtoms
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
      CALL RdSYM(.true., NAtoms, RM,     GROUP,  NTrans,
     $           NDEG,   NQ,     IUNQ,   TRANS,  NEqATM)
C
C  Open the <control> file

      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
C
C ................................................................
C  After a full Polarizability calculation, the energy, dipole moment
C  and <S**2> (if UHF) written to the <control> file will be
C  those for the last applied electric field and NOT those for
C  the original geometry. So - on first entry - read these values
C  (if extant) and save them, restoring them back to the <control>
C  file when the calculation is complete
C
      inquire(file=jobname(1:lenJ)//'.polchk',exist=isthere)
      If(.not.isthere) Then
        call fdcntrl(IUnit,7,'$energy',IEnd)
        If(IEnd.EQ.0) Then
          call rdcntrl(IUnit,7,'$energy',2,idum,E0,cdum)
          Eflag = .True.
        EndIf
        call fdcntrl(IUnit,5,'$escf',IEnd)
        If(IEnd.EQ.0) Then
          call rdcntrl(IUnit,5,'$escf',2,idum,EScf,cdum)
          SCFflag = .True.
        EndIf
        Call RdDIP(IUnit,Dip0,DFlag)
        Call RdS2(IUnit,S2,XMult,Sflag)
      EndIf
C .................................................................
C
C  Read from the <control> file
C    print flag
C    symmetry threshold
C    polarizability data
C    dipole moment
C
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      call rdcntrl(IUnit,14,'$sym_threshold',2,idum,thrsym,cdum)
      call RdPOL(IUnit,field,IType)
      call RdDIP(IUnit,Dip,DFlag)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      Call RdCoordF(IUnit,  NAtom,  AtSymb, XC,    NMol,
     $              IMOL,   XCharg, XMass)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
      DipolD = MOD(IType,10).EQ.1 ! dipole derivatives requested
      PolarD = IType.GE.10        ! polarizability derivatives requested
      Dipol = DipolD.OR.PolarD
C
C  Use of symmetry during finite-difference steps
C
      Symflag = thrsym.gt.Zero
c
      If(IPRNT.GT.0) Then
        WRITE(IOut,1000)
        WRITE(IOut,1100) field
      EndIf
C
C  Attempt to read <polchk> file
C  If no <polchk> file exists, this is the first step
C  of a new job
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.polchk',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
c
      READ(40) Eflag,SCFflag,Dflag,Sflag,E0,EScf,S2,XMult,Dip0
      READ(40) IField,EField,EField0,G0,GM,DipM,Pol,DipD,PolD,
     $         NTrans,TRANS,NEqATM,NDEG,RM,GROUP
c
      CLOSE (UNIT=40,STATUS='KEEP')
C
      GO TO 20
 10   CONTINUE
C
C  .............................................................
C  Nothing on file
C  Initialize
C
      If(IPRNT.GT.0) Then
        WRITE(IOut,1200)
        If(PolarD) WRITE(IOut,1300)
      EndIf
c
      CALL ZeroIT(DipD,NAt3*3)           ! zero dipole derivatives
      CALL ZeroIT(PolD,NAt3*6)           ! zero polarizability derivatives
C
C  initialize field indicator
C
      IField = 0
C
C  ...............................................................
C
 20   CONTINUE
C
C  If derivatives, read gradient from <grad> file
C
      If(PolarD.OR.(Dipol.AND.IField.GT.0)) Then
        CALL RdGRAD(NAtoms,GC,'save')
        CALL VScal(NAt3,-1.0d0,GC)            ! file contains forces
      EndIf
C
 30   CONTINUE
C
C  Do the finite difference Step
C
C  ----------------------------------------------------------------
C
      CALL NUMPOL(NAtoms, NTrans, IPRNT,  field,  SymFlag,
     $            GROUP,  AtSymb, XC,     Dipol,  PolarD,
     $            Dip,    DipM,   G0,     GC,     GM,
     $            V,      TRANS,  NEqATM, IField, EField,
     $            EField0,Pol,    PolD,   DipD,   Finish,
     $            IStatus)
C
C  ----------------------------------------------------------------
C
C  Write <polchk> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.polchk',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
      If(.NOT.Finish) Then
        WRITE(40) Eflag,SCFflag,Dflag,Sflag,E0,EScf,S2,XMult,Dip0
        WRITE(40) IField,EField,EField0,G0,GM,DipM,Pol,DipD,PolD,
     $            NTrans,TRANS,NEqATM,NDEG,RM,GROUP
        CLOSE (UNIT=40,STATUS='KEEP')
      Else
        CLOSE (UNIT=40,STATUS='DELETE')
      EndIf
c
      IF(Finish) THEN
c -- write dipole/polarizability derivatives to <deriv> file
        CALL WrDeriv(NAt3,DipD,PolD,DipolD,PolarD)
c -- delete <polchk> file
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.polchk',
     $        FORM='UNFORMATTED')
        CLOSE (UNIT=40,STATUS='DELETE')
c -- remove $field and $noorient flags and restore old data to <control> file
c .....................................................................
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $        FORM='FORMATTED',STATUS='OLD')
        call wrcntrl(IUnit,9,'$noorient',-1,idum,rdum,cdum)
        call wrcntrl(IUnit,6,'$field',-1,idum,rdum,cdum)
        call setival('field',0)    ! turn off field in depository
        If(Eflag) Then
         call wrcntrl(IUnit,7,'$energy',2,idum,E0,cdum)
        Else
         call fdcntrl(IUnit,7,'$energy',IEnd)
         If(IEnd.EQ.0) call wrcntrl(IUnit,7,'$energy',-1,idum,rdum,cdum)
        EndIf
        If(SCFflag) Then
         call wrcntrl(IUnit,5,'$escf',2,idum,EScf,cdum)
        Else
         call fdcntrl(IUnit,5,'$escf',IEnd)
         If(IEnd.EQ.0) call wrcntrl(IUnit,5,'$escf',-1,idum,rdum,cdum)
        EndIf
        If(Dflag) Then
         Call WrDIP(IUnit,Dip0)
        Else
         call fdcntrl(IUnit,7,'$dipole',IEnd)
         If(IEnd.EQ.0) call wrcntrl(IUnit,7,'$dipole',-1,idum,rdum,cdum)
        EndIf
        If(Sflag) Then
         Call WrS2(IUnit,S2,XMult)
        Else
         call fdcntrl(IUnit,5,'$<s2>',IEnd)
         If(IEnd.EQ.0) call wrcntrl(IUnit,5,'$<s2>',-1,idum,rdum,cdum)
        EndIf
        CLOSE (UNIT=IUnit,STATUS='KEEP')
c .....................................................................
        If(IPRNT.GT.1) WRITE(IOut,2000)
      ENDIF
C
      RETURN
c
  900 Format(9X,I4,2X,I4)
 1000 FORMAT(/,' CALCULATING POLARIZABILITY BY CENTRAL DIFFERENCE ON',
     $         ' DIPOLE MOMENT IN FIELD')
 1100 FORMAT(/,' Field strength: ',F10.6,' au')
 1200 FORMAT(/,' Dipole moment derivatives calculated by',
     $         ' finite-difference on gradient')
 1300 FORMAT(' Polarizability derivatives calculated by',
     $       ' double finite-difference on gradient')
 2000 FORMAT(//,' ***********************************************',/,
     $          ' **  APPARENTLY SUCCESSFUL POLAR CALCULATION  **',/,
     $          ' ***********************************************',//)
c
      END
c ========================================================================
c
      SUBROUTINE NUMPOL(NAtoms, NTrans, IPRNT,  field,  SymFlag,
     $                  GROUP,  AtSymb, XC,     DipolD, PolarD,
     $                  Dip,    DipM,   G0,     GC,     GM,
     $                  V,      TRANS,  NEqATM, IField, EField,
     $                  EField0,Pol,    PolD,   DipD,   Finish,
     $                  IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates polarizability and optionally its derivatives by finite-
C  difference in the presence of an external electric field
C
C  polarizability  =  dipole moment differentiated once w.r.t.field
C  dipole derivatives = gradient differentiated once w.r.t.field
C  polarizability derivatives = gradient differentiated twice w.r.t.field
C
C  ** Makes good - but not optimal - use of symmetry **
C
C  ARGUMENTS
C
C  NAtoms  -  number of real atoms
C  NTrans  -  number of symmetry operations
C  IPRNT   -  flag for controlling printout
C               0 - NO printout (except for error messages)
C               1 - summary and warning printout only
C               2 - "standard" printout
C               3 - more printout (including gradient)
C               4 - debug printout
C  field   -  field strength (atomic units)
C  Symflag -  Logical flag   .true.  - use symmetry
C                            .false. - do not use symmetry
C  GROUP   -  molecular point group
C  AtSymb  -  atomic symbols (used for nice printout)
C  XC      -  Cartesian coordinates
C  DipolD  -  logical flag for calculation of dipole derivatives
C  PolarD  -  logical flag for calculation of polarizability derivatives
C  Dip     -  current dipole moment
C  DipM    -  storage for dipole moments in fields
C  G0      -  original gradient (in zero field)
C  GC      -  current gradient
C  GM      -  storage for gradients in fields
C  V       -  general work array (3*NAtoms)
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  IField  -  field indicator
C               1  -  field in +X direction
C               2  -  field in -X direction
C               3  -  field in +Y direction
C               4  -  field in -Y direction
C               5  -  field in +Z direction
C               6  -  field in -Z direction
C               7  -  field in +X & +Y directions
C               8  -  field in -X & -Y directions
C               9  -  field in +X & +Z directions
C              10  -  field in -X & -Z directions
C              11  -  field in +Y & +Z directions
C              12  -  field in -Y & -Z directions
C  EField  -  field components (X,Y,Z)
C  EField0 -  initial field components (if any)
C
C  on final exit
C
C  Pol     -  polarizability
C  PolD    -  polarizability derivatives
C  DipD    -  dipole derivatives
C  Finish  -  Logical flag indicating status of calculation
C              .true.  -  all requested quantities calculated
C              .false. -  need to take another finite difference step
C  IStatus -  Exit status flag
C              1 - successful completion of this step
C              9 - something went wrong
C             -1 - special indicating first entry only
C
C
      REAL*8 XC(3,NAtoms),G0(3*NAtoms),GC(3*NAtoms),
     $       DipM(3,6),GM(3*NAtoms,12),V(3*NAtoms)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans)
      DIMENSION EField(3),EField0(3),Dip(3),DipD(3*NAtoms,3),
     $          Pol(3,3),PolD(3*NAtoms,6)
      CHARACTER*8 AtSymb(NAtoms+2)
      CHARACTER*4 GROUP
      Character*2 ch2
      LOGICAL Symflag,DipolD,PolarD,Finish
C
      PARAMETER (Zero=0.0d0,One=1.0d0,Two=2.0d0)
C
      DATA thrsh/1.0d-6/       ! zero threshold
C
C
C  Until all finite-difference calculations have been completed
C  all we need do is store the current dipole/gradient and prepare
C  the point charges for the next calculation
C
      NAt3 = 3*NAtoms
      Finish = .FALSE.
      IOut=ioutfil('iout')
      ICon=ioutfil('icon')
c
 20   CONTINUE
C
C
      IF(IField.EQ.6) THEN
C
C  have all dipole moments
C  calculate polarizability and dipole derivatives
C
        CALL CpyVEC(3,Dip,DipM(1,IField))
        If(DipolD) CALL CpyVEC(NAt3,GC,GM(1,IField))
        Call CalcPol(NAtoms,DipM,DipolD,GM,field,Pol,DipD)
        If(.NOT.PolarD) Then
          Finish = .TRUE.
          GO TO 50
        EndIf
      ENDIF
c
      IF(IField.EQ.0) THEN
C
C  initial step
C  save the initial gradient (if there is one) and set up
C  the first (+X) field
C
        CALL GetField0(EField0)    ! check for any existing Field
c
        If(PolarD) Call CpyVEC(NAt3,GC,G0)
        EField(1) = EField0(1) + field
        EField(2) = EField0(2)
        EField(3) = EField0(3)
        CALL SetField(EField)
c
        ch2 = '+X'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.1) THEN
cc
        CALL CpyVEC(3,Dip,DipM(1,IField))
        If(DipolD) CALL CpyVEC(NAt3,GC,GM(1,IField))
c
        If(Symflag) Then
          CALL EquivFLD(NAtoms, NTrans, 1,      NEqATM, TRANS,
     $                  DipolD, V,      Dip,    GC,     IField)
          If(IField.NE.1) Then
           If(IPRNT.GT.1) WRITE(IOut,1010)
           GO TO 20
          EndIf
        EndIf
c
        EField(1) = EField0(1) - field
        EField(2) = EField0(2)
        EField(3) = EField0(3)
        CALL SetField(Efield)
c
        ch2 = '-X'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.2) THEN
cc
        CALL CpyVEC(3,Dip,DipM(1,IField))
        If(DipolD) CALL CpyVEC(NAt3,GC,GM(1,IField))
        EField(1) = EField0(1)
        EField(2) = EField0(2) + field
        EField(3) = EField0(3)
        CALL SetField(EField)
c
        ch2 = '+Y'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.3) THEN
cc
        CALL CpyVEC(3,Dip,DipM(1,IField))
        If(DipolD) CALL CpyVEC(NAt3,GC,GM(1,IField))
c
        If(Symflag) Then
          CALL EquivFLD(NAtoms, NTrans, 2,      NEqATM, TRANS,
     $                  DipolD, V,      Dip,    GC,     IField)
          If(IField.NE.3) Then
           If(IPRNT.GT.1) WRITE(IOut,1010)
           GO TO 20
          EndIf
        EndIf
c
        EField(1) = EField0(1)
        EField(2) = EField0(2) - field
        EField(3) = EField0(3)
        CALL SetField(EField)
c
        ch2 = '-Y'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.4) THEN
cc
        CALL CpyVEC(3,Dip,DipM(1,IField))
        If(DipolD) CALL CpyVEC(NAt3,GC,GM(1,IField))
        EField(1) = EField0(1)
        EField(2) = EField0(2)
        EField(3) = EField0(3) + field
        CALL SetField(EField)
c
        ch2 = '+Z'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.5) THEN
cc
        CALL CpyVEC(3,Dip,DipM(1,IField))
        If(DipolD) CALL CpyVEC(NAt3,GC,GM(1,IField))
c
        If(Symflag) Then
          CALL EquivFLD(NAtoms, NTrans, 3,      NEqATM, TRANS,
     $                  DipolD, V,      Dip,    GC,     IField)
          If(IField.NE.5) Then
           If(IPRNT.GT.1) WRITE(IOut,1010)
           GO TO 20
          EndIf
        EndIf
c
        EField(1) = EField0(1)
        EField(2) = EField0(2)
        EField(3) = EField0(3) - field
        CALL SetField(EField)
c
        ch2 = '-Z'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.6) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
        EField(1) = EField0(1) + field
        EField(2) = EField0(2) + field
        EField(3) = EField0(3)
        CALL SetField(EField)
c
        ch2 = '+X'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
        ch2 = '+Y'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.7) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
c
        If(Symflag) Then
          CALL EquivFLD(NAtoms, NTrans, 4,      NEqATM, TRANS,
     $                  DipolD, V,      Dip,    GC,     IField)
          If(IField.NE.7) Then
           If(IPRNT.GT.1) WRITE(IOut,1010)
           GO TO 20
          EndIf
        EndIf
c
        EField(1) = EField0(1) - field
        EField(2) = EField0(2) - field
        EField(3) = EField0(3)
        CALL SetField(EField)
c
        ch2 = '-X'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
        ch2 = '-Y'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.8) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
        EField(1) = EField0(1) + field
        EField(2) = EField0(2)
        EField(3) = EField0(3) + field
        CALL SetField(EField)
c
        ch2 = '+X'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
        ch2 = '+Z'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.9) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
c
        If(Symflag) Then
          CALL EquivFLD(NAtoms, NTrans, 5,      NEqATM, TRANS,
     $                  DipolD, V,      Dip,    GC,     IField)
          If(IField.NE.9) Then
           If(IPRNT.GT.1) WRITE(IOut,1010)
           GO TO 20
          EndIf
        EndIf
c
        EField(1) = EField0(1) - field
        EField(2) = EField0(2)
        EField(3) = EField0(3) - field
        CALL SetField(EField)
c
        ch2 = '-X'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
        ch2 = '-Z'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.10) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
        EField(1) = EField0(1)
        EField(2) = EField0(2) + field
        EField(3) = EField0(3) + field
        CALL SetField(EField)
c
        ch2 = '+Y'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
        ch2 = '+Z'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.11) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
c
        If(Symflag) Then
          CALL EquivFLD(NAtoms, NTrans, 6,      NEqATM, TRANS,
     $                  DipolD, V,      Dip,    GC,     IField)
          If(IField.NE.11) Then
           If(IPRNT.GT.1) WRITE(IOut,1010)
           GO TO 20
          EndIf
        EndIf
c
        EField(1) = EField0(1)
        EField(2) = EField0(2) - field
        EField(3) = EField0(3) - field
        CALL SetField(EField)
c
        ch2 = '-Y'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
        ch2 = '-Z'
        If(IPRNT.GT.1) WRITE(IOut,1000) ch2
cc
      ELSE IF(IField.EQ.12) THEN
cc
        CALL CpyVEC(NAt3,GC,GM(1,IField))
cc
      ENDIF
C
C
      IField = IField+1
c
 50   CONTINUE
cc      write(6,*) ' dipole moment is:'
cc      call prntmat(1,3,1,dip)
cc      write(6,*) ' DipM matrix is:'
cc      call prntmat(6,3,6,dipm)
cc      write(6,*) ' GM matrix is:'
cc      call prntmat(12,nat3,12,gm)
c
      IF(IField.GT.12) Then
C
C  have all gradients
C  calculate polarizability derivatives
C
cc        write(6,*) ' About to call <CalcPolD>  G0 is:'
cc        call prntmat(1,nat3,1,g0)
        CALL CalcPolD(NAtoms,G0,GM,field,PolD)
        Finish = .TRUE.
      ENDIF
C
C
C -- final print out
C
      IF(Finish) THEN
C
C  symmetrize
C
        If(NTrans.GT.1.AND.NAtoms.GT.2) Then
c
          CALL SymPOL(NTrans,TRANS,GM,DipM,thrsh,Pol)
          If(DipolD) Then
            CALL SymDIP(NAtoms, NTrans, NEqATM, TRANS,  GC,
     $                  GM,     thrsh,  DipD)
          EndIf
          If(PolarD) Then
cc            write(6,*) ' Polarizability derivatives BEFORE symmetry'
cc            call PrntPolD(IOut,NAtoms,AtSymb,PolD)
            CALL SymPOLD(NAtoms, NTrans, NEqATM, TRANS,  DipM,
     $                   GM,     thrsh,  PolD)
cc            write(6,*) ' Polarizability derivatives AFTER symmetry'
cc            call PrntPolD(IOut,NAtoms,AtSymb,PolD)
          EndIf
        EndIf
c
        If(IPRNT.GT.1) Then
          WRITE(IOut,1100)
          WRITE(ICon,1200) field
          WRITE(IOut,1300)
          WRITE(ICon,1300)
          WRITE(IOut,1400) Pol(1,1),Pol(2,1),Pol(2,2),Pol(3,1),
     $                     Pol(3,2),Pol(3,3)
          WRITE(ICon,1400) Pol(1,1),Pol(2,1),Pol(2,2),Pol(3,1),
     $                     Pol(3,2),Pol(3,3)
          If(DipolD) Then
            CALL PrntDipD(IOut,NAtoms,AtSymb,DipD)
            CALL PrntDipD(ICon,NAtoms,AtSymb,DipD)
          EndIf
          If(PolarD) Then
            CALL PrntPolD(IOut,NAtoms,AtSymb,PolD)
            CALL PrntPolD(ICon,NAtoms,AtSymb,PolD)
          EndIf
        EndIf
cc
      ENDIF
C
      RETURN
c
 1000 FORMAT(' Setting External Field in the ',A2,' Direction')
 1010 FORMAT(' ** -Ve Field Equivalent by Symmetry - Skipping **')
 1100 FORMAT(' Finite-difference calculations completed')
 1200 FORMAT(/,' Polarizabilities calculated by finite-difference',
     $         '   field strength: ',F9.6,' au')
 1300 FORMAT(/,' <== POLARIZABILITY TENSOR  (atomic units) ==>')
 1400 FORMAT(/,6X,'XX',10X,'XY',10X,'YY',10X,'XZ',10X,'YZ',10X,'ZZ',/,
     $         6F12.6)
c
      END
c ========================================================================
c
      SUBROUTINE EquivFLD(NAtoms, NTrans, IC,     NEqATM, TRANS,
     $                    DipolD, GM,     Dip,    GC,     IField)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks to see if a -Ve field is formerly equivalent to the
C  corresponding +Ve field, and - if so - generates the step-down
C  dipole and gradient from the previous step-up quantities
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  IC      -  location of field axis
C             (1=X, 2=Y, 3=Z, 4=XY, 5=XZ, 6=YZ)
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  DipolD  -  logical flag for use of gradient
C  GM      -  work array for gradient transformation
C  Dip     -  step-up dipole
C             on exit may contain step-down dipole
C  GC      -  step-up gradient
C             on exit may contain step-down gradient
C  IField  -  current value of finite-difference step indicator
C             on exit may be incremented by 1
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans)
      DIMENSION Dip(3),GC(3,NAtoms),GM(3,NAtoms)
      Logical DipolD
C
      PARAMETER (One=1.0d0)
C
C
C  Loop over symmetry operations and see if we can find reflection
C  planes or c2 axes
C
      IF(IC.EQ.1) THEN
C
C  look for R(YZ) or C2(Y)/C2(Z) or CI
C
        DO 10 IOP=1,NTrans
        If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.One.AND.
     $     TRANS(3,3,IOP).EQ.One) Then
          dip(1) = -dip(1)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(1) = -dip(1)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.One) Then
          dip(1) = -dip(1)
          dip(2) = -dip(2)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(1) = -dip(1)
          dip(2) = -dip(2)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
        EndIf
 10     CONTINUE
cc
      ELSE IF(IC.EQ.2) THEN
C
C  look for R(XZ) or C2(X)/C2(Z) or CI
C
        DO 20 IOP=1,NTrans
        If(TRANS(1,1,IOP).EQ.One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $     TRANS(3,3,IOP).EQ.One) Then
          dip(2) = -dip(2)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(2) = -dip(2)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.One) Then
          dip(1) = -dip(1)
          dip(2) = -dip(2)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(1) = -dip(1)
          dip(2) = -dip(2)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
        EndIf
 20     CONTINUE
cc
      ELSE IF(IC.EQ.3) THEN
C
C  look for R(XY) or C2(X)/C2(Y) or CI
C
        DO 30 IOP=1,NTrans
        If(TRANS(1,1,IOP).EQ.One.AND.Trans(2,2,IOP).EQ.One.AND.
     $     TRANS(3,3,IOP).EQ.-One) Then
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(2) = -dip(2)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(1) = -dip(1)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          dip(1) = -dip(1)
          dip(2) = -dip(2)
          dip(3) = -dip(3)
          If(DipolD)
     $      CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                    GM,    GC)
          GO TO 90
        EndIf
 30     CONTINUE
cc
      ELSE IF(IC.EQ.4) THEN
C
C  look for C2(Z) or CI
C
        DO 40 IOP=1,NTrans
        If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $     TRANS(3,3,IOP).EQ.One) Then
          CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                  GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                  GM,    GC)
          GO TO 90
        EndIf
 40     CONTINUE
cc
      ELSE IF(IC.EQ.5) THEN
C
C  look for C2(Y) or CI
C
        DO 50 IOP=1,NTrans
        If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.One.AND.
     $     TRANS(3,3,IOP).EQ.-One) Then
          CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                  GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                  GM,    GC)
          GO TO 90
        EndIf
 50     CONTINUE
cc
      ELSE IF(IC.EQ.6) THEN
C
C  look for C2(X) or CI
C
        DO 60 IOP=1,NTrans
        If(TRANS(1,1,IOP).EQ.One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $     TRANS(3,3,IOP).EQ.-One) Then
          CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                  GM,    GC)
          GO TO 90
cc
        Else If(TRANS(1,1,IOP).EQ.-One.AND.Trans(2,2,IOP).EQ.-One.AND.
     $          TRANS(3,3,IOP).EQ.-One) Then
          CALL SymGradF(NAtoms,NEqATM(1,IOP),TRANS(1,1,IOP),
     $                  GM,    GC)
          GO TO 90
        EndIf
 60     CONTINUE
cc
      ENDIF
C
      RETURN
c
 90   CONTINUE
      IField = IField+1
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE SymGradF(NAtoms, NEqATM, TRANS,  GM,     GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Sets up symmetry-equivalent step-down gradient from previous
C  step-up gradient     **ABELIAN SYMMETRY ONLY**
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalencies for this symmetry operator
C  TRANS   -  symmetry operation
C  GM      -  copy of current gradient
C  GC      -  on exit contains step-down gradient
C
C
      DIMENSION NEqATM(NAtoms),TRANS(3,3),GM(3,NAtoms),GC(3,NAtoms)
C
C
      CALL CpyVEC(3*NAtoms,GC,GM)
      DO 10 IAtm=1,NAtoms
      JAtm = NEqATM(IAtm)
      GC(1,JAtm) = TRANS(1,1)*GM(1,IAtm)
      GC(2,JAtm) = TRANS(2,2)*GM(2,IAtm)
      GC(3,JAtm) = TRANS(3,3)*GM(3,IAtm)
 10   CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE GetField0(EField0)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Checks for presence of an already existing field
C
      Dimension EField0(3)
      Parameter (Zero=0.0d0)
c
      Call getival('field',ifield)
      If(ifield.eq.1) Then
        call getrval('fieldX',EField0(1))
        call getrval('fieldY',EField0(2))
        call getrval('fieldZ',EField0(3))
      Else
        EField0(1) = Zero
        EField0(2) = Zero
        EField0(3) = Zero
      EndIf
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE SetField(EField)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  writes current electric field components to both the <control>
C  file and the texas depository
C **WARNING** There is also a common block!!!
C
      DIMENSION EField(3)
      Character*256 jobname
      common /fieldx/ xfield,xfgrad,field(9)
c
      Data IUnit/1/
c
      Common /job/jobname,lenJ
C
C  write to control file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      CALL WrField(IUnit,EField)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  write to depository
C
      call setival('field',1)
      call setrval('fieldX',EField(1))
      call setrval('fieldY',EField(2))
      call setrval('fieldZ',EField(3))
C
C  put in the common block (my GOD!!!)
      field(1) = EField(1)
      field(2) = EField(2)
      field(3) = EField(3)
      xfield = 1.1d0
      xfgrad = 0.0d0
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE CalcPol(NAtoms, DipM,   DipolD, GM,     field,
     $                   Pol,    DipD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  computes the polarizability and dipole moment derivatives from the
C  various dipole moments/gradients in different electric fields
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  DipM    -  dipole moments
C  DipolD  -  logical flag for computing dipole derivatives
C  GM      -  gradients
C  field   -  strength of applied field
C
C  on exit
C
C  Pol     -  polarizabilities
C  DipD    -  dipole moment derivatives
C
C
      DIMENSION DipM(3,6),GM(3*NAtoms,6),Pol(3,3),DipD(3*NAtoms,6)
      Logical DipolD
C
C
      twoF = field + field
      fourF = twoF + twoF
C
C  compute polarizabilities
C
      Pol(1,1) = ( DipM(1,1) - DipM(1,2) )/twoF      ! XX term
      Pol(2,2) = ( DipM(2,3) - DipM(2,4) )/twoF      ! YY
      Pol(3,3) = ( DipM(3,5) - DipM(3,6) )/twoF      ! ZZ
c
      Pol(2,1) = ( DipM(2,1) - DipM(2,2) +           ! XY
     $             DipM(1,3) - DipM(1,4) )/fourF
      Pol(3,1) = ( DipM(3,1) - DipM(3,2) +           ! XZ
     $             DipM(1,5) - DipM(1,6) )/fourF
      Pol(3,2) = ( DipM(2,5) - DipM(2,6) +           ! YZ
     $             DipM(3,3) - DipM(3,4) )/fourF
c
      Pol(1,2) = Pol(2,1)
      Pol(1,3) = Pol(3,1)
      Pol(2,3) = Pol(3,2)
C
C  now compute dipole derivatives (if requested)
C
      If(DipolD) Then
        DO 10 I=1,3*NAtoms
        DipD(I,1) = -( GM(I,1) - GM(I,2) )/twoF
        DipD(I,2) = -( GM(I,3) - GM(I,4) )/twoF
        DipD(I,3) = -( GM(I,5) - GM(I,6) )/twoF
 10     CONTINUE
      EndIf
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE CalcPolD(NAtoms,G0,GM,field,PolD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  computes the polarizability derivatives from the various
C  gradients in different electric fields
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  G0      -  initial gradient (no field)
C  GM      -  gradients in various fields
C  field   -  strength of applied field
C
C  on exit
C
C  PolD    -  polarizability derivatives
C
C
      DIMENSION G0(3*NAtoms),GM(3*NAtoms,12),PolD(3*NAtoms,6)
C
C
      field2 = field**2
      Twofield2 = field2 + field2
c
      DO 10 I=1,3*NAtoms
      G0(I) = G0(I) + G0(I)
      PolD(I,1) = -( GM(I,1) + GM(I,2) - G0(I) )/field2
      PolD(I,3) = -( GM(I,3) + GM(I,4) - G0(I) )/field2
      PolD(I,6) = -( GM(I,5) + GM(I,6) - G0(I) )/field2
c
      PolD(I,2) = -( GM(I,7) + GM(I,8) + G0(I) - GM(I,1) -
     $               GM(I,2) - GM(I,3) - GM(I,4) )/Twofield2
      PolD(I,4) = -( GM(I,9) + GM(I,10) + G0(I) - GM(I,1) -
     $               GM(I,2) - GM(I,5) - GM(I,6) )/Twofield2
      PolD(I,5) = -( GM(I,11) + GM(I,12) + G0(I) - GM(I,3) -
     $               GM(I,4) - GM(I,5) - GM(I,6) )/Twofield2
 10   CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE PrntDipD(IOut,NAtoms,AtSymb,DipD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  prints the dipole derivatives
C
C  ARGUMENTS
C
C  IOut    -  unit number to print to
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  DipD    -  dipole moment derivatives
C
C
      DIMENSION DipD(3,NAtoms,3)
      CHARACTER*8 AtSymb(NAtoms)
C
c
      WRITE(IOut,1000)
c
      DO 10 IAtm=1,NAtoms
      WRITE(IOut,1100) IAtm,AtSymb(IAtm),'X',
     $                 DipD(1,IAtm,1),DipD(1,IAtm,2),DipD(1,IAtm,3)
      WRITE(IOut,1200) 'Y',DipD(2,IAtm,1),DipD(2,IAtm,2),DipD(2,IAtm,3)
      WRITE(IOut,1200) 'Z',DipD(3,IAtm,1),DipD(3,IAtm,2),DipD(3,IAtm,3)
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,' <== DIPOLE MOMENT DERIVATIVES  (atomic units) ==>',/,
     $       /,20X,'X',10X,'Y',10X,'Z')
 1100 FORMAT(I4,1X,A4,1X,A2,1X,3(F10.5,1X))
 1200 FORMAT(10X,A2,1X,3(F10.5,1X))
c
      END
c ========================================================================
c
      SUBROUTINE PrntPolD(IOut,NAtoms,AtSymb,PolD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  prints the polarizability derivatives
C
C  ARGUMENTS
C
C  IOut    -  unit number to print to
C  NAtoms  -  number of atoms
C  AtSymb  -  atomic symbols
C  PolD    -  polarizability derivatives
C
C
      DIMENSION PolD(3,NAtoms,6)
      CHARACTER*8 AtSymb(NAtoms)
C
      WRITE(IOut,1000)
c
      DO 10 IAtm=1,NAtoms
      WRITE(IOut,1100) IAtm,AtSymb(IAtm),'X',
     $                 (PolD(1,IAtm,I),I=1,6)
      WRITE(IOut,1200) 'Y',(PolD(2,IAtm,I),I=1,6)
      WRITE(IOut,1200) 'Z',(PolD(3,IAtm,I),I=1,6)
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,' <== POLARIZABILITY DERIVATIVES  (atomic units) ==>',/,
     $       /,19X,'XX',9X,'XY',9X,'YY',9X,'XZ',9X,'YZ',9X,'ZZ')
 1100 FORMAT(I4,1X,A4,1X,A2,1X,6(F10.4,1X))
 1200 FORMAT(10X,A2,1X,6(F10.4,1X))
c
      END
c ========================================================================
c
      SUBROUTINE SymPOL(NTrans,TRANS,V,W,thrsh,Pol)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Symmetrizes the polarizability tensor according to
C  all operations of the molecular point group
C  Additionally zeros all elements below thrsh
C
C  ARGUMENTS
C
C  NTrans  -  number of symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  V       -  scratch space (3*3)
C  W       -  scratch space (3*3)
C  thrsh   -  threshold below which elements will be set to zero
C  Pol     -  on exit contains symmetrized polarizabilities
C
      DIMENSION TRANS(3,3,NTrans),V(3,3),W(3,3),Pol(3,3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      If(NTrans.EQ.1) GO TO 75
c
      CALL CpyVEC(9,Pol,V)
      DO 50 IOP=2,NTrans
C
C  form TRANS * Pol * TRANS
C
      CALL ZeroIT(W,9)
      DO 20 I=1,3
      DO 19 J=1,3
      DO 18 K=1,3
      W(I,J) = W(I,J) + TRANS(I,K,IOP)*Pol(K,J)
 18   CONTINUE
 19   CONTINUE
 20   CONTINUE
c
      DO 30 I=1,3
      DO 29 J=1,3
      DO 28 K=1,3
      V(I,J) = V(I,J) + W(I,K)*TRANS(J,K,IOP)
 28   CONTINUE
 29   CONTINUE
 30   CONTINUE
c
 50   CONTINUE
C
C  scale the final polarizability tensor
C
      skal = One/FLOAT(NTrans)
      DO 60 I=1,3
      Pol(I,1) = skal*V(I,1)
      Pol(I,2) = skal*V(I,2)
      Pol(I,3) = skal*V(I,3)
 60   CONTINUE
c
 75   CONTINUE
C
C  zero out elements below thrsh
C
      DO 80 I=1,3
      If(Abs(Pol(I,1)).LT.thrsh) Pol(I,1) = Zero
      If(Abs(Pol(I,2)).LT.thrsh) Pol(I,2) = Zero
      If(Abs(Pol(I,3)).LT.thrsh) Pol(I,3) = Zero
 80   CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE SymPOLD(NAtoms, NTrans, NEqATM, TRANS,  P,
     $                   V,      thrsh,  PolD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Symmetrizes the Cartesian polarizability derivative tensor according
C  to all operations of the molecular point group
C  Additionally zeros all elements below thrsh
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  P       -  scratch space (3*6)
C  V       -  scratch space
C  thrsh   -  threshold below which elements will be set to zero
C  PolD    -  on entry contains the unsymmetrized polarizability derivatives
C             on exit contains symmetrized polarizability derivatives
C             (defined as PolD(i,atom,j,k)=dp(j,k)/dX(i,atom)
C              where i,j,k=1,3  atom=1,NAtoms)
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),
     $          P(3,6),PolD(3*NAtoms,6),V(3*NAtoms,6)
      Dimension PP(3,3,3),PP1(3,3,3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      If(NTrans.EQ.1) GO TO 85
c
      CALL CpyVEC(18*NAtoms,PolD,V)
      DO 70 IAtm=1,NAtoms
      II = 3*(IAtm-1)
c
      DO 60 IOP=2,NTrans
      KAtm = NEqATM(IAtm,IOP)
      KK = 3*(KAtm-1)
C
C  Form TRANS * DipD * TRANS(t)
C
      CALL ZeroIT(P,18)
      DO 20 I=1,3
      DO 19 J=1,6
      DO 18 K=1,3
      IT = II+K
      P(I,J) = P(I,J) + TRANS(I,K,IOP)*PolD(IT,J)
 18   CONTINUE
 19   CONTINUE
 20   CONTINUE
C
C  expand polarizability derivatives to full matrix
C
      DO I=1,3
      IT = 0
      DO J=1,3
      DO K=1,J
      IT = IT+1
      PP(I,J,K) = P(I,IT)
      PP(I,K,J) = P(I,IT)
      EndDO
      EndDO
      EndDO
c
      CALL ZeroIT(PP1,27)
C
C  transform the last polarizability component index
C
      DO 30 I=1,3
      DO 29 J=1,3
      DO 28 K=1,3
      DO 27 L=1,3
      PP1(I,J,K) = PP1(I,J,K) + TRANS(K,L,IOP)*PP(I,J,L)
 27   CONTINUE
 28   CONTINUE
 29   CONTINUE
 30   CONTINUE
C
C  transform the first polarizability component index
C
      CALL ZeroIT(PP,27)
c
      DO 40 I=1,3
      DO 39 J=1,3
      DO 38 K=1,3
      DO 37 L=1,3
      PP(I,K,J) = PP(I,K,J) + PP1(I,L,J)*TRANS(K,L,IOP)
 37   CONTINUE
 38   CONTINUE
 39   CONTINUE
 40   CONTINUE
C
C   contract polarizability derivatives and accumulate
C
      DO I=1,3
      IT = 0
      DO J=1,3
      DO K=1,J
      IT = IT+1
      P(I,IT) = PP(I,J,K)
      EndDO
      EndDO
      EndDO
c
      DO 50 I=1,3
      IT = KK+I
      DO 49 J=1,6
      V(IT,J) = V(IT,J) + P(I,J)
 49   CONTINUE
 50   CONTINUE
c
 60   CONTINUE
 70   CONTINUE
C
C  scale the final partial derivative vector
C
      skal = One/FLOAT(NTrans)
      DO 80 I=1,3*NAtoms
      PolD(I,1) = skal*V(I,1)
      PolD(I,2) = skal*V(I,2)
      PolD(I,3) = skal*V(I,3)
      PolD(I,4) = skal*V(I,4)
      PolD(I,5) = skal*V(I,5)
      PolD(I,6) = skal*V(I,6)
 80   CONTINUE
c
 85   CONTINUE
C
C  zero out elements below thrsh
C
      DO 90 I=1,3*NAtoms
      If(Abs(PolD(I,1)).LT.thrsh) PolD(I,1) = Zero
      If(Abs(PolD(I,2)).LT.thrsh) PolD(I,2) = Zero
      If(Abs(PolD(I,3)).LT.thrsh) PolD(I,3) = Zero
      If(Abs(PolD(I,4)).LT.thrsh) PolD(I,4) = Zero
      If(Abs(PolD(I,5)).LT.thrsh) PolD(I,5) = Zero
      If(Abs(PolD(I,6)).LT.thrsh) PolD(I,6) = Zero
 90   CONTINUE
C
      RETURN
      END
