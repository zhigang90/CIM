c ==================================================================
c  NUMERICAL HESSIAN MODULE            JB   October 1997
c ==================================================================
c
      subroutine prenumH(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval
c
c  reads the NUMHESS line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=4)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      character*4 options(nopt)
      character cdum*20
      character*256 jobname
      logical found,isthere
c
      parameter (IUnit=1)
c
      data options/'numh','fdst','file','prin'/
      data ioptyp/0,11,21,1/
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
c -- check for existence of <hesschk> file
c -- if exists, ensure best integral threshold
c -- (skip if semiempirical)         JB 9 April 99
c
      inquire(file=jobname(1:lenJ)//'.hesschk',exist=isthere)
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
        delta = ropval(1,2)
      else
        delta = 0.02d0     ! default finite difference step size
      endif
c
c -- finite-difference atom list from file?
c -- (copy to generic file)
c
      If(ifound(3).eq.1) Then
       call rmblan(chopval(3),256,Len)
       Call CopyFile(chopval(3)(1:Len),Len,jobname(1:lenJ)//'.hss',
     $               lenJ+4)
      EndIf
c
c -- print flag
      if(ifound(4).eq.1) then
        IPRNT = iopval(1,4)
      else
        IPRNT = 2          ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call wrcntrl(IUnit,7,'$fdstep',2,idum,delta,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c  =======================================================================
c
      SUBROUTINE HESSNUM(NMem,Z,Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** NUMERICAL HESSIAN PROGRAM **
C
C  This Program determines the Hessian (2nd derivative)
C  matrix in Cartesian coordinates by finite-difference
C  on the (analytical) gradient
C
C  FILES
C
C  Data is transferred to and from HESSNUM via the following
C  files:
C
C  <sym>      -  number of atoms, symmetry data
C  <coord>    -  current geometry (Cartesian coordinates)
C  <control>  -  program options
C  <grad>     -  current gradient
C  <hess>     -  final Cartesian Hessian matrix
C  <hesschk>  -  information pertinent to next finite difference
C                step  (produced by HESSNUM itself)
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
C
C  Now get the memory
C
      NAT3 = 3*NAtoms
      IMem = 4*NAT3**2 + 5*NAT3 + 8*NAtom + 2*NAtoms + NAtoms*MSymOP
     $                 + 9*MSymOP + NMol+1
c
cc      CALL falloc(Z(0),8*IMem,iptr,IErr)
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,7,'HESSNUM')
      CALL ZeroIT(Z,IMem)      ! clear the memory
C
C  Allocate memory pointers
C
      IXC = iptr                    !  current geometry
      IXO = IXC + 3*NAtom           !  input geometry
      IML = IXO + 3*NAtom           !  pointer to molecules
      IXCG = IML +  NMol+1          !  atomic charges
      IXCM = IXCG + NAtom           !  atomic masses
      IUQ = IXCM + NAtom            !  list of symmetry unique atoms
      ISY = IUQ + NAtoms            !   ditto as ones and zeros
      ITN = ISY + NAtoms            !  symmetry operations as 3x3 matrices
      INQ = ITN + 9*MSymOP          !  list of atomic equivalences
      IgP = INQ + NAtoms*MSymOP     !  step-up gradient
      IgM = IgP + NAT3              !  step-down gradient
      IHS = IgM + NAT3              !  Hessian matrix
      IHO = IHS + NAT3**2           !  old Hessian matrix
      IDD = IHO + NAT3**2           !  dipole moment derivatives
      IP  = IDD + 3*NAT3            !  work array
      IR  = IP  + NAT3**2           !  another work array
      IEnd = IR + NAT3**2 - 1
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,7,'HESSNUM')
C
C  ----------------------------------------------------------------------
C
      CALL NumHMAIN(NAtoms,  NAtom,   MSymOP,  Z(IXC),  Z(IXO),
     $              NMol,    Z(IML),  Z(IXCG), Z(IXCM), Z(IUQ),
     $              Z(IgP),  Z(IgM),  Z(IHS),  Z(IHO),  Z(IDD),
     $              Z(IP),   Z(IR),   Z(ISY),  Z(ITN),  Z(INQ),
     $              Cnvgd)
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
      SUBROUTINE NumHMAIN(NAtoms, NAtom,  MSymOP, XC,     XOld,
     $                    NMol,   IMOL,   XCharg, XMass,  IUNQ,
     $                    gP,     gM,     HESS,   HOld,   DipD,
     $                    P,      R,      ISYM,   TRANS,  NEqATM,
     $                    Finish)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for numerical Hessian
C  HESSMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(3,NAtom),XOld(3,NAtom),XCharg(NAtom),XMass(NAtom),
     $       gP(3*NAtoms),gM(3*NAtoms),HESS(3*NAtoms,3*NAtoms),
     $       HOld(3*NAtoms,3*NAtoms),P(9*NAtoms*NAtoms),
     $       R(9*NAtoms*NAtoms),RM(3,3)
      INTEGER IMOL(NMol+1),IUNQ(NAtoms),ISYM(NAtoms)
      DIMENSION TRANS(3,3,MSymOP),NEqATM(NAtoms,MSymOP)
      DIMENSION Dip(3),DipD(3*NAtoms,3),Dip0(3)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtom)
c ..................................................
      CHARACTER GROUP*4,cdum*20,jobname*256
      LOGICAL Dipole,Symflag,Finish,isthere
      LOGICAL Eflag,Dflag,Sflag
C
      PARAMETER (Zero=0.0d0)
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
      PARAMETER (ICond=8)        !  unit number for summary file
c
      DATA Eflag/.False./, Dflag/.False./, Sflag/.False./
c
      Common /job/jobname,lenJ
C
C
      NAT3 = 3*NAtoms
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
C  After a full Hessian calculation, the energy, dipole moment
C  and <S**2> (if UHF) written to the <control> file will be
C  those for the last finite-difference displacement and NOT
C  those for the original geometry. So - on first entry - read
C  these values (if extant) and save them, restoring them back
C  to the <control> file when the Hessian is calculated
C
      inquire(file=jobname(1:lenJ)//'.hesschk',exist=isthere)
      If(.not.isthere) Then
        call fdcntrl(IUnit,7,'$energy',IEnd)
        If(IEnd.EQ.0) Then
          call rdcntrl(IUnit,7,'$energy',2,idum,E0,cdum)
          Eflag = .True.
        EndIf
        Call RdDIP(IUnit,Dip0,DFlag)
        Call RdS2(IUnit,S2,XMult,Sflag)
      EndIf
C .................................................................
C
C  Read from the <control> file
C    print flag
C    finite-difference step size
C    symmetry threshold
C    current energy
C    dipole moment (if available)
C    <S**2> and multiplicity (UHF; if available)
C
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      call rdcntrl(IUnit,7,'$fdstep',2,idum,delta,cdum)
      call rdcntrl(IUnit,14,'$sym_threshold',2,idum,thrsym,cdum)
c .................................................................
c -- if <hesschk> exists, assume successful SCF + gradient
      inquire(file=jobname(1:lenJ)//'.hesschk',exist=isthere)
      If(isthere) Then
        call rdcntrl(IUnit,7,'$energy',2,idum,EC,cdum)
        Call RdDIP(IUnit,Dip,Dipole)
      EndIf
c .................................................................
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      Call RdCoordF(IUnit,  NAtom,  AtSymb, XOld,   NMol,
     $              IMOL,   XCharg, XMass)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
      IStatus = 1
C
C  Use of symmetry during finite-difference steps
C
      Symflag = thrsym.gt.Zero
C
C  Attempt to read <hesschk> file
C  If no <hesschk> file exists, this is the first step
C  of a new job
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.hesschk',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
c
      READ(40) Eflag,Dflag,Sflag,E0,S2,XMult,Dip0
      READ(40) IQ,IC,IS,NQ,IUNQ,ISYM,XOld,AtSymb,gP,DipD,HESS,
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
      If(IPRNT.GT.0) WRITE(ICond,1000) delta
c
      CALL ZeroIT(HESS,NAT3*NAT3)        ! zero Hessian
      CALL ZeroIT(DipD,NAT3*3)           ! zero dipole derivatives
C
C  first check if there is a list of user-defined atomic centres
C  for partial finite-difference Hessian
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.hss',
     $      FORM='FORMATTED',STATUS='OLD',ERR=96)
C
C  data found
C  first line in file should just be a title
C
      If(IPRNT.GT.1) WRITE(6,1500)
c
      READ(40,*)
C
C  read in list of atoms (maximum 10 per line)
C  (overwrite IUNQ array)
C
      CALL IZeroIT(IUNQ,NAtoms)
      IT=1
c
 12   CONTINUE
      READ (40,*,END=94) (IUNQ(I),I=IT,MIN(IT+9,NAtoms))
c
 94   DO 15 I=IT,MIN(IT+9,NAtoms)
      If(IUNQ(I).EQ.0) Then
       NQ = I-1
       GO TO 95
      EndIf
 15   CONTINUE
c
      IT = IT+10
      NQ = IT-1
      GO TO 12
c
 95   CONTINUE
      CLOSE (UNIT=40,STATUS='DELETE')
c
 96   CONTINUE
C
C  generate ISYM array
C
      CALL IZeroIT(ISYM,NAtoms)
      DO 25 I=1,NQ
      II = IUNQ(I)
      ISYM(II) = 1
 25   CONTINUE
C
C  Only the symmetry-unique atoms will be moved
C  Set up step indicators
C
      IQ = 1                   !  which unique atom we are moving
      IC = 1                   !  X, Y or Z
      IS = 1                   !  step up or step down
C
C  ** Special for diatomics only **
C  (assumed to lie entirely along Z-axis)
C
      If(NAtoms.EQ.2) IC = 3
C
C  special for first entry (no gradient, just displace)
C
      IStatus = -1
      GO TO 30
C
C  ...............................................................
C
 20   CONTINUE
C
C  Read gradient from <grad> file
C
      CALL RdGRAD(NATOMS,gM,'save')
      CALL VScal(NAT3,-1.0d0,gM)            ! file contains forces
C
 30   CONTINUE
C
C  Do the finite difference Step
C
C  ----------------------------------------------------------------
C
      CALL NUMHESS(NAtoms, NTrans, IPRNT,  delta,  Symflag,
     $             GROUP,  AtSymb, NQ,     IUNQ,   XC,
     $             XOld,   EC,     gP,     gM,     HESS,
     $             HOld,   Dipole, Dip,    DipD,   P,
     $             R,      ISYM,   TRANS,  NEqATM, RM,
     $             IQ,     IC,     IS,     Finish, IStatus)
C
C  ----------------------------------------------------------------
C
C  Write <hesschk> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.hesschk',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
      If(.NOT.Finish) Then
        WRITE(40) Eflag,Dflag,Sflag,E0,S2,XMult,Dip0
        WRITE(40) IQ,IC,IS,NQ,IUNQ,ISYM,XOld,AtSymb,gP,DipD,HESS,
     $            NTrans,TRANS,NEqATM,NDEG,RM,GROUP
        CLOSE (UNIT=40,STATUS='KEEP')
      Else
        CLOSE (UNIT=40,STATUS='DELETE')
      EndIf
C
C  write new coordinates to <coord> file
C
      If(NAtom.GT.NAtoms)
     $   CALL CpyVEC(3*(NAtom-NAtoms),XOld(1,NAtoms+1),XC(1,NAtoms+1))
      Call WrCoord(NAtom,  AtSymb, XC,     NMol,   IMOL,
     $             XCharg, XMass)
c
      IF(Finish) THEN
c -- write final Hessian to <jobname>.hes file
        CALL WrHESS(jobname(1:lenJ)//'.hess',lenJ+5,0,NAT3,HESS)
c -- write dipole derivatives to <deriv> file
        CALL WrDeriv(NAT3,DipD,jnk,.true.,.false.)
c -- delete <hesschk> file
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.hesschk',
     $        FORM='UNFORMATTED')
        CLOSE (UNIT=40,STATUS='DELETE')
c -- remove $noorient flag and restore old data to <control> file
c .....................................................................
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $        FORM='FORMATTED',STATUS='OLD')
        call wrcntrl(IUnit,9,'$noorient',-1,idum,rdum,cdum)
        If(Eflag) Then
         call wrcntrl(IUnit,7,'$energy',2,idum,E0,cdum)
        Else
         call fdcntrl(IUnit,7,'$energy',IEnd)
         If(IEnd.EQ.0) call wrcntrl(IUnit,7,'$energy',-1,idum,rdum,cdum)
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
c
c write to the <control> file an info about hessian quality
c
c    ihessq=-1 (negative) poor, crude obtained from geom.opt
c    ihessq=+1 (positive) good, obtained from Hessian calc.
c
      ihessq=+1
      call  wrcntrl(IUnit,9,'$hessqual',1,ihessq,rdum,cdum)
c .....................................................................
        CLOSE (UNIT=IUnit,STATUS='KEEP')
        If(IPRNT.GT.1) WRITE(6,2000)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,'  Hessian calculated by central differences',
     $         '   stepsize: ',F7.4,' au',/)
 1500 FORMAT(' List of Finite-Difference Atomic Centres Read In')
 2000 FORMAT(//,' *************************************************',/,
     $          ' **  APPARENTLY SUCCESSFUL HESSIAN CALCULATION  **',/,
     $          ' *************************************************',//)
c
      END
c ========================================================================
c
      SUBROUTINE NUMHESS(NAtoms, NTrans, IPRNT,  delta,  Symflag,
     $                   GROUP,  AtSymb, NQ,     IUNQ,   XC,
     $                   XOld,   EC,     gP,     gM,     HESS,
     $                   HOld,   Dipole, Dip,    DipD,   P,
     $                   R,      ISYM,   TRANS,  NEqATM, RM,
     $                   IQ,     IC,     IS,     Finish, IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates Hessian by Finite Difference on Gradient
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
C               3 - same as 2
C               4 - more printout (including gradient)
C               5 - heavier printout (including Hessian)
C               6 - debug printout
C  delta   -  finite-difference step size
C  Symflag -  Logical flag   .true.  - use symmetry
C                            .false. - do not use symmetry
C  GROUP   -  molecular point group
C  AtSymb  -  atomic symbols (used for nice printout)
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of symmetry unique atoms in ascending order
C  XC      -  current Cartesian coordinates (displaced)
C             contains new coordinates on exit
C  XOld    -  original Cartesian coordinates
C  EC      -  current energy
C  gP      -  gradient for forward displacement
C  gM      -  on entry contains current gradient
C             (also gradient for backward displacement)
C  HESS    -  Hessian matrix (partial and final)
C  HOld    -  scratch storage used for symmetrizing Hessian
C  Dipole  -  Logical flag   .true.  - dipole moment available
C                            .false. - no dipole moment
C  Dip     -  dipole moment
C  DipD    -  dipole derivatives (partial and final)
C  P       -  scratch storage used during Hessian projection
C  R       -    ditto
C  ISYM    -  list of all atoms
C              ISYM(I) = 1 if atom I symmetry-unique
C              ISYM(I) = 0 otherwise
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  RM      -  rotation matrix (redundant under this implementation)
C  IQ      -  step indicator for which unique atom we are moving
C  IC      -  which coordinate we are moving
C              1 - X;  2 - Y;  3 - Z
C  IS      -  step up or step down
C              1 - step up
C             -1 - step down
C  Finish  -  Logical flag indicating status of Hessian generation
C              .true.  -  Hessian generation completed
C              .false. -  need to take another finite difference step
C  IStatus -  Exit status flag
C              1 - successful completion of this step
C              9 - something went wrong
C             -1 - special indicating first entry only
C
C
      REAL*8 XC(3,NAtoms),XOld(3,NAtoms),gP(3*NAtoms),gM(3*NAtoms),
     $       HESS(3*NAtoms,3*NAtoms),HOld(3*NAtoms,3*NAtoms),
     $       P(9*NAtoms*NAtoms),R(9*NAtoms*NAtoms)
      INTEGER IUNQ(NAtoms),ISYM(NAtoms)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans),RM(3,3)
      DIMENSION Dip(3),DipD(3*NAtoms,3)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 GROUP
      LOGICAL Symflag,Dipole,Finish,Equiv,rotate
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
C
      DATA thrsh/1.0d-6/       ! zero threshold
C
C .................................................................
      If(IPRNT.GT.0) WRITE(6,1000)
      If(IPRNT.GT.2) CALL PrntCAR(6,0,NATOMS,AtSymb,XOld)
      If(IPRNT.GT.1) Then
       WRITE(6,1100) delta
       WRITE(6,1200)
       WRITE(6,1201) (IUNQ(I),I=1,NQ)
      EndIf
c .................................................................
c
      NAT3 = 3*NAtoms
      Equiv = .FALSE.
      Finish = .FALSE.
c
      IATOM = IUNQ(IQ)
      ICol = 3*(IATOM-1) + IC
c
      If(IStatus.EQ.-1) THEN
       IStatus = 1
       If(IPRNT.GT.1) WRITE(6,1300) IATOM,IC
       GO TO 95
      EndIf
C
 20   CONTINUE
C
C
      IF(IS.EQ.1) THEN
cc
cc  step up
cc
C  Put gradient into gP and dipole into DipD
C
       CALL CpyVEC(NAT3,gM,gP)
       If(Dipole) THEN
        DipD(ICol,1) = Dip(1)
        DipD(ICol,2) = Dip(2)
        DipD(ICol,3) = Dip(3)
       EndIF
c
       If(IPRNT.GT.3) Then                      ! debug printout
        write(6,*) ' Energy is: ',EC
        write(6,*) ' Step-Up Gradient is:'
        do i=1,natoms
        ii=3*(i-1)
        write(6,*) gP(ii+1),gP(ii+2),gP(ii+3)
        enddo
       EndIf
C
C  ....................................................................
C    ROTATION SECTION - CHECK FOR AXIS REORIENTATION
C
cc      CALL ChkROT(RM,IPRNT,thrsh,rotate)     ! disable this
      rotate = .false.
C
C  If rotation did occur, need to transform the gradient
C  into the old coordinate system
C
      IF(rotate) THEN
       CALL RotVEC(NAtoms,RM,gP)
       If(IPRNT.GT.3) Then                      ! debug printout
        write(6,*) ' Step-Up Gradient in original axis system is:'
        do i=1,natoms
        ii=3*(i-1)
        write(6,*) gP(ii+1),gP(ii+2),gP(ii+3)
        enddo
       EndIf
      ENDIF
C
C  ....................................................................
C
C  prepare for step down
C
       IS = -1
C
C  check if we need to take step down
C
       IF(Symflag.AND.XOld(IC,IATOM).EQ.Zero) THEN
        CALL EquivGRD(NAtoms,NTrans,IC,XOld(1,IATOM),NEqATM,TRANS,
     $                Dipole,gP,gM,Dip,Equiv)
        IF(Equiv) THEN
         If(IPRNT.GT.1) WRITE(6,1400)
         GO TO 20
        ENDIF
       ENDIF
c
       If(IPRNT.GT.1) WRITE(6,1500) IATOM,IC
cc
      ELSE
cc
cc step down
cc
       If(IPRNT.GT.3) Then                      ! debug printout
        write(6,*) ' Energy is: ',EC
        write(6,*) ' Step-Down Gradient is:'
        do i=1,natoms
        ii=3*(i-1)
        write(6,*) gM(ii+1),gM(ii+2),gM(ii+3)
        enddo
       EndIf
C
C  ....................................................................
C    ROTATION SECTION - CHECK FOR AXIS REORIENTATION
C
cc      CALL ChkROT(RM,IPRNT,thrsh,rotate)     ! disable this
      rotate = .false.
C
C  If rotation did occur, need to transform the gradient
C  into the old coordinate system
C
      IF(rotate.AND..NOT.Equiv) THEN
       CALL RotVEC(NAtoms,RM,gM)
       If(IPRNT.GT.3) Then                      ! debug printout
        write(6,*) ' Step-Down Gradient in original axis system is:'
        do i=1,natoms
        ii=3*(i-1)
        write(6,*) gM(ii+1),gM(ii+2),gM(ii+3)
        enddo
       EndIf
      ENDIF
C
C  ....................................................................
C
C
C  Construct this row/column of Hessian
C
       delta2 = delta + delta
c
       DO 30 J=1,NAT3
       HESS(ICol,J) = ( gP(J) - gM(J) )/delta2
 30    CONTINUE
c
       If(IPRNT.GT.4) Then                      ! debug printout
        write(6,*) ' Hessian is:'
        call prntmat(nat3,nat3,nat3,Hess)
       EndIf
C
C  Construct dipole derivative
C
       If(Dipole) Then
        DipD(ICol,1) = (DipD(ICol,1) - Dip(1))/delta2
        DipD(ICol,2) = (DipD(ICol,2) - Dip(2))/delta2
        DipD(ICol,3) = (DipD(ICol,3) - Dip(3))/delta2
       EndIf
C
C  prepare for next step up
C
       IS = 1
       IC = IC+1
       If(IC.GT.3) Then
        IC = 1
        IQ = IQ+1
       EndIF
C
C  Check for termination
C
       If(IQ.GT.NQ.OR.NAtoms.EQ.2) GO TO 98
c
       If(IPRNT.GT.1) Then
        IATOM = IUNQ(IQ)
        WRITE(6,1300) IATOM,IC
       EndIf
cc
      ENDIF
C
 95   CONTINUE
C
C  Take new step
C
      IATOM = IUNQ(IQ)
      CALL CpyVEC(NAT3,XOld,XC)
      XC(IC,IATOM) = XC(IC,IATOM) + IS*delta
C
      RETURN
C
C  ................................................................
C  ** COMPLETED FINITE-DIFFERENCE CALCULATION **
C
 98   CONTINUE
C
C  Construct full Hessian from Symmetry-unique rows
C
      IF(NTrans.GT.1.AND.NAtoms.GT.2) THEN
C
C  Copy Rows into Columns
C
       DO 50 IQ=1,NQ
       JATOM = IUNQ(IQ)
       JJ = 3*(JATOM-1)
       DO 49 IAtm=1,NAtoms
       II = 3*(IAtm-1)
       IF(ISYM(IAtm).EQ.0) THEN
        DO 42 I=1,3
        IT = II+I
        DO 41 J=1,3
        JT = JJ+J
        HESS(IT,JT) = HESS(JT,IT)
 41     CONTINUE
 42     CONTINUE
       ENDIF
c
 49    CONTINUE
 50    CONTINUE
C
C  Now generate all Hessian blocks from
C  symmetry-unique blocks
C
       CALL FulHES(NAtoms, NTrans, ISYM,   NEqATM, TRANS,
     $             P,      R,      HESS)
c
       If(IPRNT.GT.4) Then                      ! debug printout
        write(6,*) ' Full Hessian generated by symmetry is:'
        call prntmat(nat3,nat3,nat3,Hess)
       EndIf
C
C  now symmetrize
C
       CALL CpyVEC(NAT3*NAT3,HESS,HOLD)
       CALL SymHES(NAtoms, NTrans, NEqATM, TRANS,  HOLD,
     $             P,      R,      thrsh,  HESS)
C
C  generate the full set of dipole derivatives if available
C
       IF(Dipole) THEN
        CALL FulDIP(NAtoms, NTrans, ISYM,   NEqATM, TRANS,
     $              P,      R,      DipD)
c
        If(IPRNT.GT.4) Then                    ! debug printout
         write(6,*) ' Full dipole derivs generated by symmetry is:'
         call prntmat(3,nat3,3,dipd)
        EndIf
C
C  now symmetrize
C
        CALL SymDIP(NAtoms, NTrans, NEqATM, TRANS,  P,
     $              R,      thrsh,  DipD)
       ENDIF
cc
      ELSE
cc
C  No symmetry
C  Simply ensure overall matrix is symmetric
C
       FAC = Half
       If(NAtoms.EQ.2) FAC = One
       DO 60 I=2,NAT3
       DO 60 J=1,I-1
       HESS(I,J) = FAC*(HESS(I,J)+HESS(J,I))
       HESS(J,I) = HESS(I,J)
 60    CONTINUE
       If(NATOMS.EQ.2) Then
        HESS(6,6) = HESS(3,3)
        DipD(6,3) = -DipD(3,3)    ! fix dipole derivatives
       EndIf
cc
      ENDIF
C
C  finally add Unit diagonal elements if we have only
C  calculated a partial Hessian
C
      DO 70 I=1,NAT3
      If(HESS(I,I).EQ.Zero) HESS(I,I) = One
 70   CONTINUE
c
      IF(IPRNT.GT.2) THEN
       WRITE(6,1600)
       CALL PrntMAT(NAT3,NAT3,NAT3,HESS)
       If(Dipole) Then
        WRITE(6,1700)
        CALL PrntMAT(3,NAT3,3,DIPD)
       EndIf
      ENDIF
C
C  restore original coordinates
C
      CALL CpyVEC(NAT3,XOld,XC)
      Finish = .TRUE.
      RETURN
c
 1000 FORMAT(/' CALCULATING HESSIAN BY CENTRAL DIFFERENCE ON GRADIENT')
 1100 FORMAT(/,' Finite Difference Step Size: ',F9.6)
 1200 FORMAT(' Finite-Difference Atoms')
 1201 FORMAT(10I5)
 1300 FORMAT(' ** Atom No: ',I3,' Step No: ',I2,' Step Up **')
 1400 FORMAT(' ** Step Down Equivalent to Step Up - Skipping **')
 1500 FORMAT(' ** Atom No: ',I3,' Step No: ',I2,' Step Down **')
 1600 FORMAT(' Final Hessian Matrix is:')
 1700 FORMAT(' Cartesian Dipole Derivatives')
c
      END
c ========================================================================
c
      SUBROUTINE EquivGRD(NATOMS, NTrans, IC,     XC,     NEqATM,
     $                    TRANS,  Dipole, gP,     gM,     Dip,
     $                    Equiv)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks to see if a step down finite-difference step is
C  formerly equivalent to a step up, and - if so - generates
C  the step-down gradient from the previous step up one
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  NTrans  -  number of symmetry operations
C  IC      -  location of reflection plane
C             (1=X, 2=Y, 3=Z)
C  XC      -  current coordinates of atom being moved
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  Dipole  -  Logical flag   .true.  - dipole moment available
C                            .false. - no dipole moment
C  gP      -  step up gradient
C  gM      -  on exit may contain generated step down gradient
C  Dip     -  on entry contains dipole moment vector
C             on exit may contain generated step down dipole
C  Equiv   -  Logical flag
C              .true.  - step down gradient generated
C              .false. - unable to generate step down gradient
C
C
      REAL*8 gP(3,NATOMS),gM(3,NATOMS),XC(3),Dip(3)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NATOMS,NTrans)
      Dimension D(3),P(3)
      LOGICAL Dipole,Equiv,CmpVEC
C
      PARAMETER (thrsh=1.0d-5)
      PARAMETER (delta=0.005d0)
C
C
C  Find which symmetry operation converts
C   XC + delta  ---> XC - delta
C
      CALL CpyVEC(3,XC,D)
      D(IC) = D(IC) + delta
c
      DO 10 IOP=1,NTrans
      P(1) = TRANS(1,1,IOP)*D(1) + TRANS(1,2,IOP)*D(2)
     $           + TRANS(1,3,IOP)*D(3)
      P(2) = TRANS(2,1,IOP)*D(1) + TRANS(2,2,IOP)*D(2)
     $           + TRANS(2,3,IOP)*D(3)
      P(3) = TRANS(3,1,IOP)*D(1) + TRANS(3,2,IOP)*D(2)
     $           + TRANS(3,3,IOP)*D(3)
      P(IC) = -P(IC)
      If(CmpVEC(3,D,P,thrsh)) GO TO 20
 10   CONTINUE
C
C  ** NOTE: Should never get here! **
C  Unable to find suitable operator
C  Step down needs to be taken
C
      RETURN
 20   CONTINUE
C
C  At this point IOP is the operation we need
C  See which atoms are equivalent under IOP and
C  form step down gradient accordingly
C
      DO 30 IAtm=1,NATOMS
      JAtm = NEqATM(IAtm,IOP)
      gM(1,JAtm) =   TRANS(1,1,IOP)*gP(1,IAtm) +
     $               TRANS(1,2,IOP)*gP(2,IAtm) +
     $               TRANS(1,3,IOP)*gP(3,IAtm)
      gM(2,JAtm) =   TRANS(2,1,IOP)*gP(1,IAtm) +
     $               TRANS(2,2,IOP)*gP(2,IAtm) +
     $               TRANS(2,3,IOP)*gP(3,IAtm)
      gM(3,JAtm) =   TRANS(3,1,IOP)*gP(1,IAtm) +
     $               TRANS(3,2,IOP)*gP(2,IAtm) +
     $               TRANS(3,3,IOP)*gP(3,IAtm)
 30   CONTINUE
C
C  form step down dipole if available
C
      IF(Dipole) THEN
       CALL CpyVEC(3,Dip,D)
       Dip(1) = TRANS(1,1,IOP)*D(1) +
     $          TRANS(1,2,IOP)*D(2) +
     $          TRANS(1,3,IOP)*D(3)
       Dip(2) = TRANS(2,1,IOP)*D(1) +
     $          TRANS(2,2,IOP)*D(2) +
     $          TRANS(2,3,IOP)*D(3)
       Dip(3) = TRANS(3,1,IOP)*D(1) +
     $          TRANS(3,2,IOP)*D(2) +
     $          TRANS(3,3,IOP)*D(3)
      ENDIF
cc      write(6,*) ' generated dipole moment is:'
cc      write(6,*) dip
C
      Equiv = .TRUE.
      RETURN
c
 2000 FORMAT(/,2X,'***ERROR*** Unable to find Reflection Operation')
c
      END
c ========================================================================
c
      SUBROUTINE FulDIP(NATOMS, NTrans, ISYM,   NEqATM, TRANS,
     $                  IStr,   P,      DipD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates a full dipole derivative matrix from a partial matrix formed
C  by finite-difference on symmetry unique atoms only
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  NTrans  -  number of symmetry operations
C  ISYM    -  indicates which atoms are symmetry-unique
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  IStr    -  space for duplicate copy of ISYM array
C  P       -  scratch space (3*3)
C  DipD    -  Cartesian dipole derivatives
C
      REAL*8 DipD(3*NATOMS,3),TRANS(3,3,NTrans),P(3,3)
      INTEGER NEqATM(NATOMS,NTrans),ISYM(NATOMS),IStr(NATOMS)
C
C
      DO 10 I=1,NATOMS
      IStr(I) = ISYM(I)
 10   CONTINUE
C
C  Loop over unique atoms
C
      DO 100 IAtm=1,NATOMS
c
      If(ISYM(IAtm).EQ.0) GO TO 100
      II = 3*(IAtm-1)
C
C  Loop over symmetry operations
C
      DO 90 IOP=1,NTrans
C
C  find operation that converts symmetry-unique atom
C  IAtm into non-unique atom KAtm
C
      KAtm = NEqATM(IAtm,IOP)
      If(IStr(KAtm).EQ.1) GO TO 90
C
C  At this point operation IOP will convert unique atom
C  IAtm into non-unique atom KAtm
C
      KK = 3*(KAtm-1)
cc      write(6,*) ' Atom ',IAtm,' ---> Atom ',KAtm,' by operation ',iop
cc      write(6,*) ' II:',ii,' KK:',kk
      IStr(KAtm) = 1
C
C  .......................................................................
C
C  The block in row KAtm can be obtained by transforming
C  the block in row IAtm
C
      CALL ZeroIT(P,9)
      DO 50 I=1,3
      DO 49 J=1,3
      DO 48 K=1,3
      IT = II+K
      P(I,J) = P(I,J) + TRANS(I,K,IOP)*DipD(IT,J)
 48   CONTINUE
 49   CONTINUE
 50   CONTINUE
c
      DO 60 I=1,3
      IT = KK+I
      DO 59 J=1,3
      DO 58 K=1,3
      DipD(IT,J) = DipD(IT,J) + P(I,K)*TRANS(J,K,IOP)
 58   CONTINUE
 59   CONTINUE
 60   CONTINUE
cc
 90   CONTINUE
 100  CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE FulHES(NATOMS, NTrans, ISYM,   NEqATM, TRANS,
     $                  IStr,   P,      HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates a full Hessian matrix from a partial matrix formed
C  by finite-difference on symmetry unique atoms only
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  NTrans  -  number of symmetry operations
C  ISYM    -  indicates which atoms are symmetry-unique
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  IStr    -  space for duplicate copy of ISYM array
C  P       -  scratch space (3*3)
C  HESS    -  Hessian matrix
C
      REAL*8 HESS(3*NATOMS,3*NATOMS),TRANS(3,3,NTrans),P(3,3)
      INTEGER NEqATM(NATOMS,NTrans),ISYM(NATOMS),IStr(NATOMS)
C
C
      DO 10 I=1,NATOMS
      IStr(I) = ISYM(I)
 10   CONTINUE
C
C  Loop over unique atoms
C  ** Columns **
C
      DO 100 IAtm=1,NATOMS
c
      If(ISYM(IAtm).EQ.0) GO TO 100
      II = 3*(IAtm-1)
C
C  Loop over symmetry operations
C
      DO 90 IOP=1,NTrans
C
C  find operation that converts symmetry-unique atom
C  IAtm into non-unique atom KAtm
C
      KAtm = NEqATM(IAtm,IOP)
      If(IStr(KAtm).EQ.1) GO TO 90
C
C  At this point operation IOP will convert unique atom
C  IAtm into non-unique atom KAtm
C
      KK = 3*(KAtm-1)
cc      write(6,*) ' Atom ',IAtm,' ---> Atom ',KAtm,' by operation ',iop
cc      write(6,*) ' II:',ii,' KK:',kk
      IStr(KAtm) = 1
C
C  .......................................................................
C
C  The blocks in row KAtm can be obtained by transforming
C  appropriate blocks in row IAtm
C
C  Loop over non-unique atoms
C  ** Rows **
C
      DO 80 JAtm=1,NATOMS
c
      JJ = 3*(JAtm-1)
C
C  Determine the atom LAtm that JAtm is converted to
C  under operation IOP
C
      LAtm = NEqATM(JAtm,IOP)
      If(ISYM(LAtm).EQ.1) GO TO 80      ! already got this block
      LL = 3*(LAtm-1)
cc      write(6,*) ' Atom ',JAtm,' ---> Atom ',LAtm,' by operation ',iop
cc      write(6,*) ' JJ:',jj,' LL:',ll
C
C  Now transform the symmetry-unique Hessian Matrix block
C  HESS(IAtm,JAtm) into HESS(KAtm,LAtm)
C
cc      write(6,*) ' FORMING HESSIAN BLOCK ',katm,latm
      CALL ZeroIT(P,9)
      DO 50 I=1,3
      DO 49 J=1,3
      JT = JJ+J
      DO 48 K=1,3
      IT = II+K
      P(I,J) = P(I,J) + TRANS(I,K,IOP)*HESS(IT,JT)
 48   CONTINUE
 49   CONTINUE
 50   CONTINUE
c
      DO 60 I=1,3
      IT = KK+I
      DO 59 J=1,3
      JT = LL+J
      DO 58 K=1,3
      HESS(IT,JT) = HESS(IT,JT) + P(I,K)*TRANS(J,K,IOP)
 58   CONTINUE
 59   CONTINUE
 60   CONTINUE
cc
 80   CONTINUE
cc
 90   CONTINUE
 100  CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE SymDIP(NAtoms, NTrans, NEqATM, TRANS,  P,
     $                  V,      thrsh,  DipD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Symmetrizes the Cartesian dipole derivative vector according to
C  all operations of the molecular point group
C  Additionally zeros all elements below thrsh
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  P       -  scratch space (3*3)
C  V       -  scratch space
C  thrsh   -  threshold below which elements will be set to zero
C  DipD    -  on input contains the unsymmetrized dipole derivatives
C             on exit contains symmetrized dipole derivatives
C             (this array is equivalent to DipD(3,NAtoms,3) where the first
C             subscript is the atomic coordinate (x,y,z) and the last is
C             the dipole direction (x,y,z))
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),
     $          P(3,3),DipD(3*NAtoms,3),V(3*NAtoms,3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      If(NTrans.EQ.1) GO TO 75
c
      CALL CpyVEC(9*NAtoms,DipD,V)
      DO 10 IAtm=1,NAtoms
      II = 3*(IAtm-1)
c
      DO 20 IOP=2,NTrans
      KAtm = NEqATM(IAtm,IOP)
      KK = 3*(KAtm-1)
C
C  Form TRANS * DipD * TRANS(t)
C
      CALL ZeroIT(P,9)
      DO 50 I=1,3
      DO 49 J=1,3
      DO 48 K=1,3
      IT = II+K
      P(I,J) = P(I,J) + TRANS(I,K,IOP)*DipD(IT,J)
 48   CONTINUE
 49   CONTINUE
 50   CONTINUE
c
      DO 60 I=1,3
      IT = KK+I
      DO 59 J=1,3
      DO 58 K=1,3
      V(IT,J) = V(IT,J) + P(I,K)*TRANS(J,K,IOP)
 58   CONTINUE
 59   CONTINUE
 60   CONTINUE
c
 20   CONTINUE
 10   CONTINUE
C
C  scale the final partial derivative vector
C
      skal = One/FLOAT(NTrans)
      DO 70 I=1,3*NAtoms
      DipD(I,1) = skal*V(I,1)
      DipD(I,2) = skal*V(I,2)
      DipD(I,3) = skal*V(I,3)
 70   CONTINUE
c
 75   CONTINUE
C
C  zero out elements below thrsh
C
      DO 80 I=1,3*NAtoms
      If(Abs(DipD(I,1)).LT.thrsh) DipD(I,1) = Zero
      If(Abs(DipD(I,2)).LT.thrsh) DipD(I,2) = Zero
      If(Abs(DipD(I,3)).LT.thrsh) DipD(I,3) = Zero
 80   CONTINUE
C
      RETURN
      END
