c ==================================================================
c  NUMERICAL GRADIENT MODULE            JB   November 2009
c ==================================================================
c
      subroutine prenumG(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval
c
c  reads the NUMGRAD line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=3)
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
      data options/'numg','fdst','prin'/
      data ioptyp/0,11,1/
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
c -- check for existence of <gradchk> file
c -- if exists, ensure best integral threshold
c -- (skip if semiempirical)
c
      inquire(file=jobname(1:lenJ)//'.gradchk',exist=isthere)
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
c -- print flag
      if(ifound(3).eq.1) then
        IPRNT = iopval(1,3)
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
      SUBROUTINE GRADNUM(NMem,Z,Cnvgd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** NUMERICAL GRADIENT PROGRAM **
C
C  This Program determines the gradient (1st derivative)
C  vector in Cartesian coordinates by finite-difference
C  on the Energy
C
C  FILES
C
C  Data is transferred to and from GRADNUM via the following
C  files:
C
C  <sym>      -  number of atoms, symmetry data
C  <coord>    -  current geometry (Cartesian coordinates)
C  <control>  -  program options and current energy
C  <grad>     -  final Cartesian gradient
C  <gradchk>  -  information pertinent to next finite difference
C                step  (produced by GRADNUM itself)
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
      IMem = 2*NAT3 + 8*NAtom + 2*NAtoms + NAtoms*MSymOP
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
      IGC = INQ + NAtoms*MSymOP     !  gradient vector
      IV  = IGC + NAT3              !  work array
      IEnd = IV + NAT3 - 1
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,7,'GRADNUM')
C
C  ----------------------------------------------------------------------
C
      CALL NumGMAIN(NAtoms,  NAtom,   MSymOP,  Z(IXC),  Z(IXO),
     $              NMol,    Z(IML),  Z(IXCG), Z(IXCM), Z(IUQ),
     $              Z(ISY),  Z(ITN),  Z(INQ),  Z(IGC),  Z(IV),
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
      SUBROUTINE NumGMAIN(NAtoms, NAtom,  MSymOP, XC,     XOld,
     $                    NMol,   IMOL,   XCharg, XMass,  IUNQ,
     $                    ISYM,   TRANS,  NEqATM, GC,     V,
     $                    Finish)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for numerical gradient
C  GRADMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(3,NAtom),XOld(3,NAtom),XCharg(NAtom),XMass(NAtom),
     $       GC(3,NAtom),V(3,NAtom),RM(3,3)
      INTEGER IMOL(NMol+1),IUNQ(NAtoms),ISYM(NAtoms)
      DIMENSION TRANS(3,3,MSymOP),NEqATM(NAtoms,MSymOP)
      DIMENSION Dip(3),Dip0(3)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtom)
c ..................................................
      CHARACTER GROUP*4,cdum*20,jobname*256
      LOGICAL Symflag,Finish,isthere
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
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
C
C ................................................................
C  After a full gradient calculation, the energy, dipole moment
C  and <S**2> (if UHF) written to the <control> file will be
C  those for the last finite-difference displacement and NOT
C  those for the original geometry. So - on first entry - read
C  these values (if extant) and save them, restoring them back
C  to the <control> file when gradient calculation is complete
C
      inquire(file=jobname(1:lenJ)//'.gradchk',exist=isthere)
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
c -- if <gradchk> exists, assume successful SCF
      inquire(file=jobname(1:lenJ)//'.gradchk',exist=isthere)
      If(isthere) Then
        call rdcntrl(IUnit,7,'$energy',2,idum,EC,cdum)
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
C  Attempt to read <gradchk> file
C  If no <gradchk> file exists, this is the first step
C  of a new job
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.gradchk',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=10)
c
      READ(40) Eflag,Dflag,Sflag,E0,S2,XMult,Dip0
      READ(40) IQ,IC,IS,NQ,IUNQ,ISYM,XOld,AtSymb,EP,GC,
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
      CALL ZeroIT(GC,NAT3)        ! zero gradient
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
C  special for first entry (no energy, just displace)
C
      IStatus = -1
C
C  ...............................................................
C
 20   CONTINUE
C
C  Do the finite difference Step
C
C  ----------------------------------------------------------------
C
      CALL NUMGRAD(NAtoms, NTrans, IPRNT,  delta,  Symflag,
     $             GROUP,  AtSymb, NQ,     IUNQ,   XC,
     $             XOld,   EC,     EP,     GC,     V,
     $             ISYM,   TRANS,  NEqATM, RM,     IQ,
     $             IC,     IS,     Finish, IStatus)
C
C  ----------------------------------------------------------------
C
C  Write <gradchk> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.gradchk',
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
      If(.NOT.Finish) Then
        WRITE(40) Eflag,Dflag,Sflag,E0,S2,XMult,Dip0
        WRITE(40) IQ,IC,IS,NQ,IUNQ,ISYM,XOld,AtSymb,EP,GC,
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
c -- write final gradient to <jobname>.grad file
        CALL VScal(NAT3,-1.0d0,GC)      ! we need forces on the file
        CALL WrGRAD(NAtoms,GC)
c -- delete <gradchk> file
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.gradchk',
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
c .....................................................................
        CLOSE (UNIT=IUnit,STATUS='KEEP')
        If(IPRNT.GT.1) WRITE(6,2000)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,'  Gradient calculated by central differences',
     $         '   stepsize: ',F7.4,' au',/)
 2000 FORMAT(//,' **************************************************',/,
     $          ' **  APPARENTLY SUCCESSFUL GRADIENT CALCULATION  **',/,
     $         ' **************************************************',//)
c
      END
c ========================================================================
c
      SUBROUTINE NUMGRAD(NAtoms, NTrans, IPRNT,  delta,  Symflag,
     $                   GROUP,  AtSymb, NQ,     IUNQ,   XC,
     $                   XOld,   EC,     EP,     GC,     V,
     $                   ISYM,   TRANS,  NEqATM, RM,     IQ,
     $                   IC,     IS,     Finish, IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates Gradient by Finite Difference on Energy
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
C  EP      -  step-up energy
C  GC      -  gradient vector (partial and final)
C  V       -  scratch array
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
C  Finish  -  Logical flag indicating status of gradient generation
C              .true.  -  gradient generation completed
C              .false. -  need to take another finite difference step
C  IStatus -  Exit status flag
C              1 - successful completion of this step
C              9 - something went wrong
C             -1 - special indicating first entry only
C
C
      REAL*8 XC(3,NAtoms),XOld(3,NAtoms),GC(3,NAtoms),V(3*NAtoms)
      INTEGER IUNQ(NAtoms),ISYM(NAtoms)
      DIMENSION TRANS(3,3,NTrans),NEqATM(NAtoms,NTrans),RM(3,3)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 GROUP
      LOGICAL Symflag,Dipole,Finish,Equiv
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
 10   CONTINUE
      IATOM = IUNQ(IQ)
c
      IF(IStatus.EQ.-1) THEN
C
C  We are ready to start
C  First thing to check is the geometry
C  If any coordinate is exactly zero, it is assumed the gradient is zero also
C  So skip that step
C
        If(XOld(IC,IATOM).EQ.Zero) Then
          If(IPRNT.GT.1) WRITE(6,1400) IATOM,IC
          IC = IC+1
          IS = 1
          If(IC.GT.3) Then
            IC = 1
            IQ = IQ+1
          EndIF
          GO TO 10
        EndIf
        IStatus = 1
        If(IPRNT.GT.1) WRITE(6,1300) IATOM,IC
        GO TO 95
      ENDIF
C
C  If we got here, then we are ready to start the finite-difference
c ......................................................................
      IF(IS.EQ.1) THEN
cc
cc  step up
cc
C  copy current energy into EP
C
       EP = EC
C
C  prepare for step down
C
       IS = -1
c
       If(IPRNT.GT.1) WRITE(6,1500) IATOM,IC
cc
      ELSE
cc
cc step down
cc
C
C  Construct current gradient component
C
       delta2 = delta + delta
c
       GC(IC,IAtom) = (EP-EC)/delta2
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
 20    CONTINUE
       IATOM = IUNQ(IQ)
C
C  Check if step should be skipped
C
       If(XOld(IC,IATOM).EQ.Zero) Then
         If(IPRNT.GT.1) WRITE(6,1400) IATOM,IC
         IC = IC+1
         IS = 1
         If(IC.GT.3) Then
           IC = 1
           IQ = IQ+1
         EndIF
C
C  Check for termination
C
         If(IQ.GT.NQ.OR.NAtoms.EQ.2) GO TO 98
         GO TO 20
       EndIf
c
       If(IPRNT.GT.1) WRITE(6,1300) IATOM,IC
cc
      ENDIF
c ......................................................................
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
      IF(NTrans.GT.1.AND.NAtoms.GT.2) THEN
C
C  Construct full gradient from Symmetry-unique elements
C
       CALL FulGrad(NAtoms, NTrans, ISYM,   NEqATM, TRANS,
     $              V,      GC)
c
       If(IPRNT.GT.4) Then                      ! debug printout
        write(6,*) ' Full gradient generated by symmetry is:'
        call prntmat(natoms,3,natoms,GC)
       EndIf
C
C  now symmetrize
C
       CALL SymVEC(NAtoms, NTrans, NEqATM, TRANS,  V,
     $             thrsh,  GC)
cc
      ELSE
cc
C  No symmetry
C  do nothing
cc
      ENDIF
C
      IF(IPRNT.GT.2) THEN
       WRITE(6,1600)
       CALL PrntMAT(NAtoms,3,NAtoms,GC)
      ENDIF
C
C  restore original coordinates
C
      CALL CpyVEC(NAT3,XOld,XC)
      Finish = .TRUE.
      RETURN
c
 1000 FORMAT(/' CALCULATING GRADIENT BY CENTRAL DIFFERENCE ON ENERGY')
 1100 FORMAT(/,' Finite Difference Step Size: ',F9.6)
 1200 FORMAT(' Finite-Difference Atoms')
 1201 FORMAT(10I5)
 1300 FORMAT(' ** Atom No: ',I3,' Step No: ',I2,' Step Up **')
 1400 FORMAT(' ** Atom No: ',I3,' Step No: ',I2,' Skipped **')
 1500 FORMAT(' ** Atom No: ',I3,' Step No: ',I2,' Step Down **')
 1600 FORMAT(' Final Gradient is:')
c
      END
c ========================================================================
c
      SUBROUTINE FulGRAD(NAtoms, NTrans, ISYM,   NEqATM, TRANS,
     $                   IStr,   GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates a full gradient vector from a partial vector formed
C  by finite-difference on symmetry unique atoms only
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  ISYM    -  indicates which atoms are symmetry-unique
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  IStr    -  space for duplicate copy of ISYM array
C  GC      -  gradient vector
C
      REAL*8 GC(3,NAtoms),TRANS(3,3,NTrans)
      INTEGER NEqATM(NAtoms,NTrans),ISYM(NAtoms),IStr(NAtoms)
C
C
      DO 10 I=1,NAtoms
      IStr(I) = ISYM(I)
 10   CONTINUE
C
C  Loop over unique atoms
C
      DO 100 IAtm=1,NAtoms
c
      If(ISYM(IAtm).EQ.0) GO TO 100
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
cc      write(6,*) ' Atom ',IAtm,' ---> Atom ',KAtm,' by operation ',iop
      IStr(KAtm) = 1
C
C  .......................................................................
C
C  The gradient for KAtm can be obtained by transforming
C  the gradient for IAtm
C
      DO 50 I=1,3
      DO 49 J=1,3
      GC(I,KAtm) = GC(I,KAtm) + TRANS(I,J,IOP)*GC(J,IAtm)
 49   CONTINUE
 50   CONTINUE
c
 90   CONTINUE
 100  CONTINUE
C
      RETURN
      END
