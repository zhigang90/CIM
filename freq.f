c ==================================================================
c  VIBRATIONAL FREQUENCY & THERMODYNAMICS MODULE     JB   Nov 1997
c ==================================================================
c
      subroutine prefreq(inp)
      implicit real*8(a-h,o-z)
      character*256 chopval,jobname
c
      Common /job/jobname,lenJ
c
c  reads the FREQ line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=6)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      real*8 delta,thresh
      character*4 options(nopt)
      character cdum*20
      logical found
c
      parameter (IUnit=1)
c
      data options/'freq','temp','pres','lowf','mass','prin'/
      data ioptyp/0,13,13,11,21,1/
c-----------------------------------------------------------
      call getival('iout',iout)
      write(iout,5001)
 5001 format(/72('='))
      write(iout,*)
     *  '        Vibrational Frequency & Thermodynamics Module '
      write(iout,*)' '
c
      call f_lush(iout)
c-----------------------------------------------------------
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c -- temperature for thermodynamic analysis (in degree K)
      if(ifound(2).eq.1) then
        Temp = ropval(1,2)
        Tend = ropval(2,2)
        Tstep = ropval(3,2)
c -- error check
        If(Temp.LT.0.0d0.OR.Tend.LT.0.0d0)
     $    CALL nerror(3,'FREQUENCY module',
     $   'Thermodynamic Analysis meaningless below Zero Kelvin!',0,0)
        If(Tend.NE.0.0d0.AND.Temp.GT.Tend)
     $    CALL nerror(4,'FREQUENCY module',
     $   'Start Temperature must be LESS than end Temperature',0,0)
        If(Tstep.NE.0.0d0.AND.Tstep.GT.(Tend-Temp))
     $    CALL nerror(5,'FREQUENCY module',
     $   'Temperature step greater than Temperature range!',0,0)
        If(Tend.NE.0.0d0.AND.Tstep.EQ.0.0d0)
     $    Tstep = Tend-Temp
      else
        Temp = 298.18d0
        Tend = 0.0d0
        Tstep = 0.0d0
      endif
c -- pressure for thermodynamic analysis (in atm)
      if(ifound(3).eq.1) then
        Pres = ropval(1,3)
        Pend = ropval(2,3)
        Pstep = ropval(3,3)
c -- error check
        If(Pres.LT.0.0d0.OR.Pend.LT.0.0d0)
     $    CALL nerror(6,'FREQUENCY module',
     $   'Thermodynamic Analysis meaningless below Zero Pressure!',0,0)
        If(Pend.NE.0.0d0.AND.Pres.GT.Pend)
     $    CALL nerror(7,'FREQUENCY module',
     $   'Start Pressure must be LESS than end Pressure',0,0)
        If(Pstep.NE.0.0d0.AND.Pstep.GT.(Pend-Pres))
     $    CALL nerror(8,'FREQUENCY module',
     $   'Pressure step greater than Pressure range!',0,0)
        If(Pend.NE.0.0d0.AND.Pstep.EQ.0.0d0)
     $    Pstep = Pend-Pres
      else
        Pres = 1.0d0
        Pend = 0.0d0
        Pstep = 0.0d0
      endif
c -- lowest frequency (i.e., smallest acceptable frequency in wavenumbers)
c    values less than this in magnitude will be removed as "zero" during print-out
      if(ifound(4).eq.1) then
        FLow = ropval(1,4)
      else
        FLow = 10.0d0
      endif
c -- mass weighting (isotope abundance average or most common isotope)
c    keywords are:  MASS=abun(dant)  MASS=aver(age)
      IMass=0
      if(ifound(5).eq.1) then
        if(chopval(5)(1:4).EQ.'abun') IMass=1
      endif
c -- print flag
      if(ifound(6).eq.1) then
        IPRNT = iopval(1,6)
      else
        IPRNT = 3          ! default print flag
      endif
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,5,'$temp',2,idum,Temp,cdum)
      call wrcntrl(IUnit,5,'$tend',2,idum,Tend,cdum)
      call wrcntrl(IUnit,6,'$tstep',2,idum,Tstep,cdum)
      call wrcntrl(IUnit,5,'$pres',2,idum,Pres,cdum)
      call wrcntrl(IUnit,5,'$pend',2,idum,Pend,cdum)
      call wrcntrl(IUnit,6,'$pstep',2,idum,Pstep,cdum)
      call wrcntrl(IUnit,5,'$lowF',2,idum,FLow,cdum)
      call wrcntrl(IUnit,5,'$mass',1,IMass,rdum,cdum)
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      return
      end
c  =======================================================================
c
      SUBROUTINE FREQ(NMem,Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ..........................................................
C  ** VIBRATIONAL FREQUENCY AND THERMODYNAMICS PROGRAM **
C
C  This programs reads a Hessian in Cartesian Coordinates
C  and carries out a frequency analysis (+ some thermodynamics)
C
C  FILES
C
C  Data is transferred to and from FREQ via the following files:
C
C  <sym>      -  number of atoms, symmetry data
C  <coord>    -  current geometry (Cartesian coordinates)
C  <control>  -  program options
C  <hess>     -  Cartesian Hessian matrix
C  <deriv>    -  dipole & polarizability derivatives
C  <aat>      -  atomic axial tensors
C  ..........................................................
C
      DIMENSION Z(NMem)
      Character*20 cdum
      Character*256 jobname
C
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
      PARAMETER (MaxIR=16)       !  no. of IRs in "largest" group
c
      Common /job/jobname,lenJ
C
C
C  Read from the <sym> file
C    number of atoms
C    number of symmetry operations
C    point group symbol
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,dum,cdum)
      call rdcntrl(IUnit,6,'$group',3,idum,dum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C ...................................................................
C ** WARNING **
C    If the point group is c*v or d*h it will be defaulted to
C    c6v/d6h respectively. The standard defaults are c2v/d2h.
C    Increase the number of symmetry operations to reflect this
C
      If(cdum(1:3).EQ.'c*v') NTrans=12
      If(cdum(1:3).EQ.'d*h') NTrans=24
C ...................................................................
C
C
C  Now get the memory
C
      NAT3 = 3*NAtoms
c
      NScr = 27 + 3*NTrans + MaxIR + (NTrans+3)*MaxIR*MaxIR
      NPM  = MAX(NAT3*NAT3,48)       ! not enough memory for D*h diatomics!
c
      IMem = 2*NAT3*NAT3 + NPM + 21*NAT3 + NAtoms + 9*NTrans
     $         + 4*NAtoms*NTrans + NTrans*MaxIR + NScr
c
cc      CALL falloc(Z(0),8*IMem,iptr,IErr)
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,4,'FREQ')
C
C  Allocate memory pointers
C
      IXC = iptr                 !  geometry
      IHS = IXC + NAT3           !  Hessian matrix
      IDD = IHS + NAT3*NAT3      !  dipole derivatives
      IPD = IDD + 3*NAT3         !  polarizability derivatives
      INQ = IPD + 6*NAT3         !  list of atomic equivalences
      ITN = INQ + NAtoms*NTrans  !  symmetry operations as 3x3 matrices
      ICH = ITN + 9*NTrans       !  character table
      ICN = ICH + NTrans*MaxIR   !  representation matrix
      IPM = ICN + NTrans*NAT3    !  scratch matrix
      ISV = IPM + NPM            !  scratch vector
      IAM = ISV + 4*NAT3         !  atomic masses
      IAA = IAM + NAtoms         !  atomic axial tensors
      IVB = IAA + 3*NAT3         !  vibrational frequencies
      IVM = IVB + NAT3           !  normal modes
      IDN = IVM + NAT3*NAT3      !  dipole derivatives over normal modes
      IEnd = IDN + 3*NAT3
C
C  general scratch storage
C
      IScr = IEnd
      IEnd = IScr + NScr
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,4,'FREQ')
C
C  ----------------------------------------------------------------------
C
      CALL FREQMAIN(NAtoms,  NTrans,  MaxIR,   Z(IXC),  Z(IHS),
     $              Z(IDD),  Z(IPD),  Z(INQ),  Z(ITN),  Z(ICH),
     $              Z(ICN),  Z(IPM),  Z(ISV),  Z(IAM),  Z(IVB),
     $              Z(IVM),  Z(IAA),  Z(IDN),  NScr,    Z(IScr))
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
c ======================================================================
c
      SUBROUTINE FREQMAIN(NAtoms, MTrans, MaxIR,  XC,     HESS,
     $                    DipD,   PolD,   NEqATM, TRANS,  CHARAC,
     $                    ChNorm, P,      SVec,   AtMASS, Vib,
     $                    RMode,  AAT,    SNDip,  NMem,   Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for frequency analysis
C  FREQMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      REAL*8 XC(3,NATOMS),HESS(3*NATOMS,3*NATOMS),Vib(3*NATOMS),
     $       P(3*NATOMS,3*NATOMS),RMode(3*NATOMS,3*NATOMS),
     $       AtMASS(NATOMS),SVec(3*NATOMS,4),DipD(3*NATOMS,3),
     $       PolD(3*NATOMS,6),CHARAC(MTrans,MaxIR),
     $       ChNorm(MTrans,3*NATOMS),AAT(3*NAtoms,3),
     $       SNDip(3,3*NATOMS)
      DIMENSION TRANS(3,3,MTrans),NEqATM(NATOMS,MTrans),RM(3,3)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 ITyp(MaxIR)
      CHARACTER*3 Label(3*NAtoms,3)
c ..................................................
      CHARACTER GROUP*4,cdum*20
      character*256 jobname
      LOGICAL Dipole,Polar,AtomAx
C
      DIMENSION Z(NMem)
C
      PARAMETER (IUnit=1)                  !  unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
      NAT3 = 3*NATOMS
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
     $           NDEG,   NQ,     Z,      TRANS,  NEqATM)
C
C  Read from the <control> file
C    temperature & pressure data, lowest "non-zero" frequency,
C    print flag and hessian quality flag :
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,5,'$temp',2,idum,Temp,cdum)
      call rdcntrl(IUnit,5,'$tend',2,idum,Tend,cdum)
      call rdcntrl(IUnit,6,'$tstep',2,idum,Tstep,cdum)
      call rdcntrl(IUnit,5,'$pres',2,idum,Pres,cdum)
      call rdcntrl(IUnit,5,'$pend',2,idum,Pend,cdum)
      call rdcntrl(IUnit,6,'$pstep',2,idum,Pstep,cdum)
      call rdcntrl(IUnit,5,'$lowF',2,idum,FLow,cdum)
      call rdcntrl(IUnit,5,'$mass',1,IMass,dum,cdum)
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)  !   print flag
c
c
      call fdcntrl(IUnit,9,'$hessqual',iendh)
      if(iendh.eq.0) then
         call rdcntrl(IUnit,9,'$hessqual',1,IHESSQ,dum,cdum) ! hess quality flag
      else
         ihessq=0
      endif
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      if(ihessq.le.0) then
         write(6,*)'                    *** WARNING *** '
         if(ihessq.lt.0) then
            write(6,*)
     *      '       Frequencies are calculated from a crude Hessian'
            write(6,*)
     *      '            obtained during a geometry optimiaztion     '
         else
            write(6,*)
     *      '   Frequencies are calculated from unknown quality Hessian'
            write(6,*)
     *      '            obtained probably in previous steps           '
         endif
         write(6,*)'   '
      endif
c
      write(6,*)
      If(Tend.NE.0.0d0) WRITE(6,1000) Temp,Tend,Tstep
      If(Pend.NE.0.0d0) WRITE(6,1100) Pres,Pend,Pstep
      If(IMass.EQ.1) Then
        WRITE(6,1200)
      Else
        WRITE(6,1201)
      EndIf
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call RdCoordF(IUnit,NAtoms,AtSymb,XC,-1,jnk,Z,AtMASS)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read in Cartesian Hessian matrix
C
      CALL RdHESS(jobname(1:lenJ)//'.hess',lenJ+5,NAT3,IPRNT,
     $            HESS,IHess)
c
      If(IHess.NE.0) Then
        Call nerror(1,'FREQUENCY module',
     $                'No <hess> file found',0,0)
      EndIf
C
C  See if there are any dipole derivatives
C
      CALL RdDeriv(NAT3,DipD,PolD,Dipole,Polar)
C
C  See if there are any atomic axial tensors
C
      CALL RdAAT(NAT3,AAT,AtomAx)
C
C  Now do frequency analysis!
C
      CALL VibFREQ(NAtoms, NTrans, MTrans, MaxIR,  GROUP,
     $             AtSymb, XC,     HESS,   Dipole, DipD,
     $             Polar,  PolD,   IPRNT,  NEqATM, TRANS,
     $             CHARAC, ITyp,   ChNorm, P,      SVec,
     $             Label,  AtMASS, Vib,    RMode,  AAT,
     $             AtomAx, SNDip,  Temp,   Tend,   Tstep,   
     $             Pres,   Pend,   Pstep,  FLow,   IMass,
     $             NMem,   Z)
C
      RETURN
c
 1000 FORMAT(' Thermodynamic properties will be computed at the',
     $       ' following Temperatures(K):',/,
     $       ' Start Temperature: ',F10.5,' End Temperature: ',F10.5,
     $       ' Step: ',F9.5)
 1100 FORMAT(' Thermodynamic properties will be computed at the',
     $       ' following Pressures(atm):',/,
     $       ' Start Pressure: ',F10.5,' End Pressure: ',F10.5,
     $       ' Step: ',F9.5)
 1200 FORMAT(' Most Abundant Isotope Atomic Masses will be Used')
 1201 FORMAT(' Isotope-averaged Atomic Masses will be Used')
c
      END
c ======================================================================
c
      SUBROUTINE VibFREQ(NAtoms, NTrans, MTrans, MaxIR,  GROUP,
     $                   AtSymb, XC,     HESS,   Dipole, DipD,
     $                   Polar,  PolD,   IPRNT,  NEqATM, TRANS,
     $                   CHARAC, ITyp,   ChNorm, P,      SVec,
     $                   Label,  AtMASS, Vib,    RMode,  AAT,
     $                   AtomAx, SNDip,  Temp,   Tend,   Tstep,
     $                   Pres,   Pend,   Pstep,  FLow,   IMass,
     $                   NMem,   Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines vibrational frequencies from Hessian matrix
C  and calls <THERMODYN> to calculate thermodynamic properties
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  MTrans  -  number of symmetry operations assuming linear groups
C             taken to be c6v/d6h, respectively
C             NOTE - MTrans is the actual dimension of arrays
C  MaxIR   -  maximum number of different irreps possible in
C             "largest" point group
C  GROUP   -  molecular point group
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C  HESS    -  Hessian matrix in Cartesian coordinates
C  Dipole  -  Logical flag   .true.  - dipole derivatives available
C                            .false. - no dipole derivatives
C  DipD    -  dipole derivatives
C             (in 3 columns containing all X Y Z derivatives)
C  Polar   -  Logical flag   .true.  - polarizability derivatives available
C                            .false. - no polarizability derivatives
C  PolD    -  polarizability derivatives
C             (in 6 columns containing all XX XY YY XZ YZ ZZ derivatives)
C  IPRNT   -  controls level of printout
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  CHARAC  -  storage for character table
C  ITyp    -  storage for irreducible representation labels
C  ChNorm  -  storage for representation matrix
C  P       -  storage (9*NATOMS*NATOMS)
C  SVec    -  storage (12*NATOMS)
C  Label   -  character storage for symmetry labels
C  AAT     -  atomic axial tensors, if any
C  AtomAx  -  Logical flag   .true.  - atomic axial tensors available
C                            .false. - no atomic axial tensors
C  SNDip   -  storage for dipole moment derivatives over normal modes
C  Temp    -  start temperature for thermodynamic analysis (K)
C  Tend    -  end temperature for thermodynamic analysis (may be zero)
C  Tstep   -  temperature step (may be zero)
C  Pres    -  start pressure for thermodynamic analysis (atm)
C  Pend    -  end pressure for thermodynamic analysis (may be zero)
C  Pstep   -  pressure step (may be zero)
C  FLow    -  lowest "non-zero" frequency  (in cm**-1)
C             (any frequency lower in magnitude than this value will be
C              removed during print-out)
C  IMass   -  0 - use isotope-averaged atomic masses
C             1 - use most abundant isotope atomic masses
C
C  on exit
C
C  AtMASS  -  atomic masses
C  Vib     -  vibrational frequencies
C  RMode   -  normal modes
C
C  NMem    -  amount of scratch storage
C  Z       -  scratch storage
C
      REAL*8 XC(3,NATOMS),HESS(3*NATOMS,3*NATOMS),Vib(3*NATOMS),
     $       P(3*NATOMS,3*NATOMS),RMode(3*NATOMS,3*NATOMS),
     $       AtMASS(NATOMS),SVec(3*NATOMS,4),DipD(3*NATOMS,3),
     $       PolD(3*NATOMS,6),CHARAC(MTrans,MaxIR),
     $       ChNorm(MTrans,3*NATOMS),AAT(3*NATOMS,3),SNDip(3,3*NATOMS)
      DIMENSION TRANS(3,3,MTrans),NEqATM(NATOMS,MTrans)
      CHARACTER*8 AtSymb(NATOMS)
      CHARACTER*3 Label(3*NAtoms,3)
      CHARACTER*4 GROUP,ITyp(MaxIR)
      LOGICAL Dipole,Polar,AtomAx
C
      DIMENSION Z(NMem)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
      PARAMETER (ToCM=5.89141d-07,ToCM2=3.94562d0)
      PARAMETER (Vthrsh=1.0d-4)
C
C
      NAT3 = 3*NATOMS
C
C  set the value of TolZero
C
      TolZero = ToCM*(FLow/ToCM2)**2
C
C  get default atomic masses
C
cc      CALL GETATNO(NATOMS,AtSymb,SVec)
cc      CALL DefMASS(NATOMS,SVec,AtMASS)
      DO 5 I=1,NATOMS
      Call nugrep(AtSymb(I)(1:2),INum)
      If(IMass.EQ.1) Then
        CALL AbunMASS(1,INum,AtMASS(I))
      Else 
        If(AtMASS(I).EQ.Zero) CALL DefMASS(1,INum,AtMASS(I))
      EndIf
 5    CONTINUE
C
C  Mass-Weight the Hessian
C
      CALL HessWT(NATOMS,AtMASS,HESS)
c
      IF(IPRNT.GT.4) THEN
       WRITE(6,1000)
       CALL PrntMAT(NAT3,NAT3,NAT3,HESS)
      ENDIF
C
C  Project out translations/rotations from Mass-Weighted Hessian
C
      CALL ProjTRM(NATOMS,IPRNT,XC,AtMASS,P,RMode,HESS)
c
      IF(IPRNT.GT.4) THEN
       WRITE(6,1100)
       CALL PrntMAT(NAT3,NAT3,NAT3,HESS)
      ENDIF
C
C  Diagonalize
C
      CALL CpyVEC(NAT3*NAT3,HESS,RMode)
      CALL DIAGMAT(RMode,NAT3,P,SVec,Vib,IErr)
c
      IF(IPRNT.GT.2) THEN
       WRITE(6,1200)
       WRITE(6,1300) (SIGN(SQRT(ABS(Vib(I))),Vib(I)),I=1,NAT3)
      ENDIF
C
C  Remove Zero frequencies
C
      CALL ChkHES(NAT3,0,0,RMode,Vib,.false.,Jnk,TolZero,NCon,NVib)
C
C  Find the number of negative frequencies
C
      CALL FndNEG(NVib,Vib,NEG)
C
C  Convert frequencies to cm**-1
C  Mass Weight and normalize the normal modes
C  Calculated IR intensities (using unnormalized modes)
C
      DO 10 I=1,NVib
      Vib(I) = SIGN(SQRT(Abs(Vib(I))/ToCM),Vib(I))*ToCM2
      CALL MassWT(-1,NATOMS,AtMASS,RMode(1,I))
      If(Dipole) CALL IRTNSE(NAT3,RMode(1,I),DipD,SNDip(1,I),SVec(I,1))
      If(Polar) CALL RAMAN(NAT3,RMode(1,I),PolD,SVec(I,2),SVec(I,3))
      If(AtomAx) CALL VCD(NAT3,RMode(1,I),DipD,AAT,SVec(I,4))
      VNRM = One/SQRT(SProd(NAT3,RMode(1,I),RMode(1,I)))
      CALL VScal(NAT3,VNRM,RMode(1,I))
 10   CONTINUE
C
C  Set factor for relative intensities
C
      IF(Dipole) THEN
       Vibmx = Zero
       DO 20 I=1,NVib
       Vibmx = MAX(Vibmx,SVec(I,1))
 20    CONTINUE
       If(Vibmx.GT.Vthrsh) Then
        Vibmx = 100.0d0/Vibmx
       Else
        Vibmx = Zero
       EndIf
      ENDIF
C
C ..................................................................
C  Get character table and IR symmetry labels
C
      CALL IrrepLBL(MaxIR,  MTrans, GROUP,  TRANS,  NMem,
     $              Z,      NSym,   CHARAC, ITyp)
C
C  Symmetry operations have been regenerated, possibly in a
C  new order.  Regenerate equivalent atoms array
C
      CALL EQVATM(NATOMS, XC,     MTrans, TRANS,  TolZero,
     $            NEqATM, IErr)
C
C  Find symmetry labels
C
      CALL FrqSYML(NAtoms, NVib,   NSym,   MTrans, IPRNT,
     $             TRANS,  NEqATM, CHARAC, ITyp,   RMode,
     $             Vib,    ChNorm, P(1,1), P(1,4), Label(1,1),
     $             Label(1,2), Label(1,3))
C ..................................................................
C
C  Print Out
C  ---------
C
      IOut=igetival('iout')
      ICond=igetival('icond')
      CALL NormOUT(IOut,   NAtoms, AtSymb,    NVib,      RMode,
     $             Vib,    Dipole, SVec(1,1), Polar,     SVec(1,2),
     $           SVec(1,3),AtomAx, SVec(1,4), Label(1,1),Label(1,2),
     $           Label(1,3),SNDip,1)
c
      CALL NormCOND(ICond, NVib,      Vib,       Dipole,  SVec(1,1),
     $              Polar, SVec(1,2), SVec(1,3), AtomAx,  SVec(1,4),
     $              Label(1,1), Label(1,2), Label(1,3))
C
C  Thermodynamic Analysis
C
      CALL THERMODYN(NAtoms, GROUP,  XC,     AtSymb, AtMASS,
     $               Temp,   Tend,   Tstep,  Pres,   Pend,
     $               Pstep,  NEG,    NVib,   Vib,    P)
C
      RETURN
c
 1000 FORMAT(/,' Mass-Weighted Hessian Matrix:')
 1100 FORMAT(/,' Projected Mass-Weighted Hessian Matrix:')
 1200 FORMAT(/,' Vibrational Frequencies in atomic units')
 1300 FORMAT(1X,6F12.6)
c
      END
c ======================================================================
c
      SUBROUTINE COM(NATOMS,ATMASS,XC,X,Y,Z)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transform into centre-of-mass coordinate system
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  ATMASS  -  atomic masses
C  XC      -  on input  original coordinates
C             on exit   centre-of-mass coordinates
C  X       -  on exit contains X centre-of-mass in original coordinate frame
C  Y       -  on exit contains Y centre-of-mass in original coordinate frame
C  Z       -  on exit contains Z centre-of-mass in original coordinate frame
C
      REAL*8 ATMASS(NATOMS),XC(3,NATOMS)
C
      PARAMETER (Zero=0.0d0)
C
C
      TOTMAS = Zero
      X = Zero
      Y = Zero
      Z = Zero
c
      DO 10 IAtm=1,NATOMS
      ami = AtMASS(IAtm)
      TOTMAS = TOTMAS + ami
      X = X + ami*XC(1,IAtm)
      Y = Y + ami*XC(2,IAtm)
      Z = Z + ami*XC(3,IAtm)
 10   CONTINUE
c
      X = X/TOTMAS
      Y = Y/TOTMAS
      Z = Z/TOTMAS
c
      DO 20 IAtm=1,NATOMS
      XC(1,IAtm) = XC(1,IAtm) - X
      XC(2,IAtm) = XC(2,IAtm) - Y
      XC(3,IAtm) = XC(3,IAtm) - Z
 20   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE DefMASS(NATOMS,IAN,AtMASS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Assign default atomic masses using isotope-averaged atomic masses
C
      DIMENSION IAN(NATOMS),AtMASS(NATOMS)
      parameter (nelem=92)
      dimension amass(nelem)
c
c --- define atomic weights for the first <nelem> elements ---
c --- data available on-line from NIST Standard Reference Database 144
c --- Atomic Weights and Isotopic Compositions   July 2010
c
c                         H - Zr
      data (amass(i), i = 1 , 40) /
     1  1.00794d0,   4.002602d0,  6.941d0,     9.012182d0, 10.811d0,
     2 12.0107d0,   14.0067d0,   15.9994d0,   18.998403d0, 20.1797d0,
     3 22.989769d0, 24.3050d0,   26.981539d0, 28.0855d0,   30.973762d0,
     4 32.065d0,    35.453d0,    39.948d0,    39.0983d0,   40.078d0,
     5 44.955912d0, 47.867d0,    50.9415d0,   51.9961d0,   54.938045d0,
     6 55.845d0,    58.933195d0, 58.6934d0,   63.546d0,    65.38d0,
     7 69.723d0,    72.64d0,     74.921597d0, 78.96d0,     79.904d0,
     8 83.798d0,    85.4678d0,   87.62d0,     88.905848d0, 91.224d0 /
c                         Nb - Hg
      data (amass(i), i = 41 , 80) /
     1  92.906378d0, 95.96d0,     98.000d0,   101.07d0,    102.905504d0,
     2 106.42d0,    107.8682d0,  112.411d0,   114.818d0,   118.710d0,
     3 121.760d0,   127.60d0,    126.904473d0,131.293d0,   132.905452d0,
     4 137.327d0,   138.90547d0, 140.116d0,   140.907653d0,144.242d0,
     5 145.000d0,   150.36d0,    151.964d0,   157.25d0,    158.925347d0,
     6 162.500d0,   164.930322d0,167.259d0,   168.934213d0,173.054d0,
     7 174.9668d0,  178.49d0,    180.94788d0, 183.84d0,    186.207d0,
     8 190.23d0,    192.217d0,   195.084d0,   196.966569d0,200.59d0 /
c                         Tl - U
      data (amass(i), i = 81 , 92) /
     1 204.3833d0,  207.2d0,     208.980399d0,209.000d0,   210.000d0,
     2 222.000d0,   223.000d0,   226.000d0,   227.000d0,   232.038055d0,
     3 231.035884d0,238.02891d0 /
c
c
      DO 10 IATM=1,NATOMS
      IATNO = IAN(IATM)
      IF(IATNO.LT.1.OR.IATNO.GT.nelem) THEN
        Call nerror(2,'GEOMETRY/FREQUENCY module',
     $    'Unknown Atom  Cannot assign atomic mass',IATM,0)
      ENDIF
      AtMASS(IATM) = amass(IATNO)
 10   CONTINUE
c
      RETURN
      END
c ======================================================================
c
      SUBROUTINE AbunMASS(NATOMS,IAN,AtMASS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Assign atomic masses using the most adundant isotope
C
      DIMENSION IAN(NATOMS),AtMASS(NATOMS)
      parameter (nelem=92)
      dimension amass(nelem)
c
c --- define atomic weights for the first <nelem> elements ---
c --- data available on-line from NIST Standard Reference Database 144
c --- Atomic Weights and Isotopic Compositions   July 2010
c
c                         H - Zr
      data (amass(i), i = 1 , 40) /
     1  1.007825d0,  4.002603d0,  7.016005d0,  9.012182d0, 11.009305d0,
     2 12.000000d0, 14.003074d0, 15.994915d0, 18.998403d0, 19.992440d0,
     3 22.989769d0, 23.985042d0, 26.981539d0, 27.976926d0, 30.973762d0,
     4 31.972071d0, 34.968853d0, 39.962383d0, 38.963707d0, 39.962591d0,
     5 44.955912d0, 47.947946d0, 50.943960d0, 51.940508d0, 54.938045d0,
     6 55.934937d0, 58.933195d0, 57.935343d0, 62.929598d0, 63.929142d0,
     7 68.925574d0, 73.921178d0, 74.921597d0, 79.916521d0, 78.918337d0,
     8 83.911507d0, 84.911790d0, 87.905612d0, 88.905848d0, 89.904704d0 /
c                         Nb - Hg
      data (amass(i), i = 41 , 80) /
     1  92.906378d0, 97.905408d0, 97.907216d0,101.904349d0,102.905504d0,
     2 105.903486d0,106.905097d0,113.903359d0,114.903878d0,119.902195d0,
     3 120.903816d0,129.906224d0,126.904473d0,131.904154d0,132.905452d0,
     4 137.905247d0,138.906353d0,139.905439d0,140.907653d0,141.907723d0,
     5 144.912749d0,151.919732d0,152.921230d0,157.924104d0,158.925347d0,
     6 163.929175d0,164.930322d0,165.930293d0,168.934213d0,173.938862d0,
     7 174.940772d0,179.946550d0,180.947996d0,183.950931d0,186.955753d0,
     8 191.961481d0,192.962926d0,194.964791d0,196.966569d0,201.970643d0/
c                         Tl - U
      data (amass(i), i = 81 , 92) /
     1 204.974428d0,207.976652d0,208.980399d0,208.982430d0,208.987148d0,
     2 222.017578d0,223.019736d0,226.025410d0,227.027752d0,232.038055d0,
     3 231.035884d0,238.050788d0 /
c
c
      DO 10 IATM=1,NATOMS
      IATNO = IAN(IATM)
      IF(IATNO.LT.1.OR.IATNO.GT.nelem) THEN
        Call nerror(2,'GEOMETRY/FREQUENCY module',
     $    'Unknown Atom  Cannot assign atomic mass',IATM,0)
      ENDIF
      AtMASS(IATM) = amass(IATNO)
 10   CONTINUE
c
      RETURN
      END
c ======================================================================
c
      SUBROUTINE FrqSYML(NATOMS, NVib,   NSym,   NTrans, IPRNT,
     $                   TRANS,  NEqATM, CHARAC, ITyp,   UMode,
     $                   Vib,    ChNorm, V1,     V2,     SLabel,
     $                   ILabel, RLabel)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines symmetry labels for the normal modes
C  and assigns IR and Raman activity
C
C  ARGUMENTS
C
C  NATOMS  -  number of real atoms
C  NVib    -  number of modes
C  NSym    -  number of different irreducible representations
C  NTrans  -  number of symmetry operations
C  IPRNT   -  flag for controlling printout
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  CHARAC  -  character table
C  ITyp    -  irreducible representation labels
C  UMode   -  normal modes
C  Vib     -  vibrational frequencies
C  ChNorm  -  scratch storage for representation matrix
C  V1      -  vector scratch storage
C  V2      -   ditto
C
C  on exit
C
C  SLabel  -  symmetry labels for each mode
C  ILabel  -  IR activity
C               YES - mode is IR active
C               NO  - mode is IR inactive
C  RLabel  -  Raman activity
C               YES - mode is Raman active
C               NO  - mode is Raman inactive
C
C
      DIMENSION TRANS(9,NTrans),NEqATM(NATOMS,NTrans),
     $          CHARAC(NTrans,NSym),UMode(3*NATOMS,NVib),Vib(NVib)
      REAL*8 ChNorm(NTrans,NVib),V1(*),V2(*)
      CHARACTER*4 ITyp(NSym)
      CHARACTER*3 SLabel(NVib)
      CHARACTER*3 ILabel(NVib),RLabel(NVib)
      CHARACTER unknown*3
C
      PARAMETER (Zero=0.0d0,TollZero=2.0d-1)
C
C
      NAT3 = 3*NATOMS
C
C ---------------------------------------------------------------------
C      get symmetry label and determine IR/Raman activity
C ---------------------------------------------------------------------

      unknown = '?'//'?'//'?'         ! to avoid pre-processor error on IBM
C
C  Loop over symmetry operations
C
      DO 40 IOP=1,NTrans
C
C  Loop over normal modes
C
      DO 30 IVec=1,NVib
c
      CALL CpyVEC(NAT3,UMode(1,IVec),V1)
C
C  generate new vector y = TRANS(x)
C
      DO 10 IAtm=1,NATOMS
      II = 3*(IAtm-1) + 1
      JAtm = NEqATM(IAtm,IOP)
      JJ = 3*(JAtm-1) + 1
      call mult3(V2(JJ),TRANS(1,IOP),V1(II),1)
 10   CONTINUE
c
c  calculate character of representation matrix
c  corresponding to vector x
c
      ss = Zero
      DO 20 K=1,NAT3
 20   ss = ss + V2(K)*V1(K)
c
      ChNorm(IOP,IVec) = ss
c
 30   CONTINUE
c
 40   CONTINUE
cc      write(6,*) ' Characters of normal modes'
cc      do i=1,nvib
cc      write(6,2000) (chnorm(j,i),j=1,NTrans)
cc      enddo
cc 2000 format(8F8.5)
C
C
C  now determine the symmetry labels
C  look for non-degenerate IRs first, then 2, 3,
C  4 and 5-fold degenerate
C
      DO 50 IVec=1,NVib
 50   SLabel(IVec) = unknown
      ifound = 0
      NDegen = 0
c
 100  CONTINUE
      NDegen = NDegen + 1
      IDegen = NDegen - 1
c
      IVec = 1
 150  CONTINUE
      If(SLabel(IVec).NE.unknown) GO TO 200
c
      DO 60 ISym=1,NSym
      DO 59 IOP=1,NTrans
      Char = Zero
      DO 58 I=0,IDegen
      Char = Char + ChNorm(IOP,IVec+I)
 58   CONTINUE
      If(Abs(Char-CHARAC(IOP,ISym)).GT.TollZero) GO TO 60
 59   CONTINUE
C
C  found a match
C
      DO 55 I=0,IDegen
      SLabel(IVec+I) = ITyp(ISym)(1:3)
      DO 54 IOP=1,NTrans
      ChNorm(IOP,IVec+I) = CHARAC(IOP,ISym)
 54   CONTINUE
 55   CONTINUE
      IVec = IVec + IDegen
      ifound = ifound + NDegen
      GO TO 200
 60   CONTINUE
c
 200  IVec = IVec + 1
      If(IVec.LE.NVib-IDegen) GO TO 150
C
C  has every mode been assigned?
C
      If(ifound.LT.NVib.AND.NDegen.LT.5) GO TO 100
C
C  if we've under or over assigned modes
C  looks like problems
C
      If(ifound.NE.NVib.AND.IPRNT.GT.0) WRITE(6,1000)
C
C ................................................................
C
C  final analysis : which fundamentals may be
C  IR- or RAMAN-active ?
C
C  form dipole/quadrapole symmetry factors
C
      DO 70 I=1,NTrans
      V1(I) = TRANS(1,I) + TRANS(5,I) + TRANS(9,I)
      V2(I) = TRANS(1,I)*TRANS(1,I) + TRANS(5,I)*TRANS(5,I) +
     $        TRANS(9,I)*TRANS(9,I) +
     $        TRANS(2,I)*TRANS(4,I) + TRANS(1,I)*TRANS(5,I) +
     $        TRANS(3,I)*TRANS(7,I) + TRANS(1,I)*TRANS(9,I) +
     $        TRANS(6,I)*TRANS(8,I) + TRANS(5,I)*TRANS(9,I)
 70   CONTINUE
c
      DO 80 IVec=1,NVib
c
      If(SLabel(IVec).EQ.unknown) Then
       ILabel(IVec) = unknown
       RLabel(IVec) = unknown
       GO TO 80
      EndIf
C
C  --- IR
C
      ss = Zero
      DO 81 IOP=1,NTrans
 81   ss = ss + V1(IOP)*ChNorm(IOP,IVec)
      nir = NINT(ss/NTrans)
      ILabel(IVec) = 'NO '
      If(nir.GT.0) ILabel(IVec) = 'YES'
cc      If(Abs(ss).GT.TollZero) ILabel(IVec) = 'YES'
C
C  --- RAMAN
C
      ss = Zero
      DO 82 IOP=1,NTrans
 82   ss = ss + V2(IOP)*ChNorm(IOP,IVec)
      nraman = NINT(ss/NTrans)
      RLabel(IVec) = 'NO '
      If(nraman.GT.0) RLabel(IVec) = 'YES'
cc      If(Abs(ss).GT.TollZero) RLabel(IVec) = 'YES'
c
 80   CONTINUE
C
      RETURN
c
 1000 FORMAT('**WARNING** Incomplete Symmetry Analysis of Normal',
     $        ' Modes')
c
      END
c ======================================================================
c
      SUBROUTINE HessWT(NATOMS,AtMASS,HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Mass-Weight the Force Constant Matrix
C
      REAL*8 AtMASS(NATOMS),HESS(3*NATOMS,3*NATOMS)
C
      PARAMETER (One=1.0d0)
C
      DO 20 IATM=1,NATOMS
      AMASI = One/SQRT(AtMASS(IATM))
      II = 3*(IATM-1)
      DO 20 I=1,3
      DO 10 JATM=1,NATOMS
      AMAS = AMASI/SQRT(AtMASS(JATM))
      JJ = 3*(JATM-1)
      DO 10 J=1,3
      HESS(II+I,JJ+J) = HESS(II+I,JJ+J)*AMAS
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE IrrepLBL(MaxIR,  NTrans, GROUP,  TRANS,  NMem,
     $                    Z,      NSym,   CHARAC, ITyp)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determine character table and irreducible representation labels
C  based on point group symbol
C
C  ARGUMENTS
C
C  MaxIR   -  maximum number of different IRs possible in
C             "largest" point group
C  NTrans  -  number of symmetry operations
C  GROUP   -  point group symbol
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  NMem    -  available scratch memory (double words)
C  Z       -  scratch storage
C
C  on exit
C
C  NSym    -  number of different IRs
C  CHARAC  -  character table
C  ITyp    -  IR labels
C
C
      DIMENSION TRANS(3,3,NTrans),CHARAC(NTrans,MaxIR)
      CHARACTER*4 ITyp(MaxIR),GROUP,Sflies
      CHARACTER*1 csf,cla
C
      DIMENSION Z(NMem)
C
      Parameter (thrsym=1.0d-5)
C
C
C  if no symmetry, set up for c1 and exit
C
cc      IF(GROUP.EQ.'c1  ') THEN
cc        NSym = 1
cc        CHARAC(1,1) = 1.0d0
cc        ITyp(1) = 'a  '
cc        RETURN
cc      ENDIF
C
C  set up scratch pointers
C
      igen = 1
      invt = igen + 3*9
      imi  = invt + NTrans
      imj  = imi  + NTrans
      idim = imj  + NTrans
      i1   = idim + MaxIR
      i2   = i1   + 3*MaxIR*MaxIR
      IEnd = i2   + NTrans*MaxIR*MaxIR - 1
      CALL MemCHK(NMem,IEnd,8,'IrrepLBL')
C
C
C  generate symmetry information
C
C  --------------------------------------------------
C  ** WARNING ** Take defaults for linear molecules
C  --------------------------------------------------
      Sflies = GROUP
      If(Sflies.EQ.'d*h ') Sflies = 'd6h '
      If(Sflies.EQ.'c*v ') Sflies = 'c6v '
C
C -- set up generators of group and group transformation matrices
      nprt = 0
      CALL grpsmb(Sflies,csf,cla,nn,nprt)
      CALL getgen(csf,nn,cla,Z(igen),ngen,IErr,nprt)
      CALL groupL(Z(igen),ngen,   TRANS,  nt,     NTrans,
     $            thrsym, Z(invt),Z(imi), Z(imj), nprt,
     $            IErr)
C
C -- derive number (NSym), type (ITyp), and dimension (idim) of real
C -- group representations
      call getrep(csf,    nn,     cla,    NSym,   ITyp,
     $            Z(idim),MaxIR,  nsymsq, nsumid)
C
C -- find generating matrices of all irreducible representations
      call repgen(csf,    nn,     cla,    Z(igen),ngen,
     $            nsymsq, NSym,   ITyp,   Z(idim),Z(i1))
C
C -- calculate all representation matrices
      call grep(NTrans, Z(imi), Z(imj), NSym,    nsymsq,
     $          Z(idim),CHARAC, Z(i2),  Z(i1))
c
cc      call prntchar(NTrans,NSym,CHARAC,ITyp)
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE IRTNSE(NAT3, UMode, DipD, SNDip, SIR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transform Cartesian dipole derivatives to normal mode
C  and get IR intensity
C
C  ARGUMENTS
C
C  NAT3    -  3 x number of atoms
C  UMode   -  normal mode
C  DipD    -  dipole derivatives
C
C  on exit
C
C  SNDip   -  dipole derivatives over normal modes
C  SIR     -  IR intensity
C
C
      REAL*8 UMode(NAT3),DipD(NAT3,3),SNDip(3)
C
      PARAMETER (Zero=0.0d0)
      PARAMETER (FCnvrt=975.12d0)  ! convert au -->  km/mol
C
C
C  Transform dipole derivatives
C
      dd1 = Zero
      dd2 = Zero
      dd3 = Zero
      DO 10 I=1,NAT3
      dd1 = dd1 + DipD(I,1)*UMode(I)
      dd2 = dd2 + DipD(I,2)*UMode(I)
      dd3 = dd3 + DipD(I,3)*UMode(I)
 10   CONTINUE
C
C  now get intensities
C
      SIR = (dd1*dd1 + dd2*dd2 + dd3*dd3)*FCnvrt
c
      SNDip(1) = dd1
      SNDip(2) = dd2
      SNDip(3) = dd3
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE MassWT(IMW,NATOMS,AtMASS,V)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transforms a vector (either the geometry or the gradient)
C  into/from mass-weighted coordinate frame
C
C  ARGUMENTS
C
C  IMW     -  flag indicating how to transform
C              1 - scale by (AtMASS)**1/2
C             -1 - scale by (AtMASS)**-1/2
C  NATOMS  -  number of atoms
C  AtMASS  -  array of atomic masses
C  V       -  vector to be transformed
C
C
      REAL*8 AtMASS(NATOMS),V(3,NATOMS)
C
      PARAMETER (One=1.0d0)
C
C
      IF(IMW.EQ.1) THEN
cc
C  transform coordinates
C
       DO 10 IATM=1,NATOMS
       AMS = SQRT(AtMASS(IATM))
       DO 10 J=1,3
       V(J,IATM) = AMS*V(J,IATM)
 10    CONTINUE
cc
      ELSE IF(IMW.EQ.-1) THEN
cc
C  transform gradient/back-transform coordinates
C
       DO 20 IATM=1,NATOMS
       AMS = One/SQRT(AtMASS(IATM))
       DO 20 J=1,3
       V(J,IATM) = AMS*V(J,IATM)
 20    CONTINUE
cc
      ENDIF
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE MomINRT(NATOMS, XC,     AtMASS, TOTMAS, SCR1,
     $                   SCR2,   T,      PMom,   ITOP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculates principal moments of inertia
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  XC      -  Cartesian coordinates
C  AtMASS  -  atomic masses
C  TOTMAS  -  total mass
C  SCR1    -  general scratch storage
C  SCR2    -   ditto
C  T       -  scratch storage for inertia tensor
C  PMom    -  on exit contains principal moments of inertia
C  ITOP    -  on exit contains number of different moments
C
C
      REAL*8 XC(3,NATOMS),AtMASS(NATOMS),SCR1(9),SCR2(9),
     $       T(9),PMom(3)
C
      PARAMETER (Zero=0.0d0,thrsh=1.0d-5)
C
C
C  initialize
C
      CALL ZeroIT(T,9)
      CX = Zero
      CY = Zero
      CZ = Zero
C
C  find centre of mass
C
      DO 10 IAtm=1,NATOMS
      ami = AtMASS(IAtm)
      CX = CX + ami*XC(1,IAtm)
      CY = CY + ami*XC(2,IAtm)
      CZ = CZ + ami*XC(3,IAtm)
 10   CONTINUE
c
      CX = CX/TOTMAS
      CY = CY/TOTMAS
      CZ = CZ/TOTMAS
C
C  form inertia tensor
C
      DO 20 IAtm=1,NATOMS
      X = XC(1,IAtm) - CX
      Y = XC(2,IAtm) - CY
      Z = XC(3,IAtm) - CZ
      ami = AtMASS(IAtm)
      T(1) = T(1) + ami*(Y*Y + Z*Z)
      T(5) = T(5) + ami*(X*X + Z*Z)
      T(9) = T(9) + ami*(X*X + Y*Y)
      T(2) = T(2) - ami*X*Y
      T(3) = T(3) - ami*X*Z
      T(6) = T(6) - ami*Y*Z
 20   CONTINUE
      T(4) = T(2)
      T(7) = T(3)
      T(8) = T(6)
C
C  diagonalize T
C
      CALL DIAGMAT(T,3,SCR1,SCR2,PMom,IErr)
c
      IF(IErr.NE.0) THEN
       WRITE(6,1000)
       CALL OptEXIT(9)
      ENDIF
C
C  determine the number of different principal moments
C
      ITOP = 3
      If(Abs(PMom(1)-PMom(2)).LT.thrsh) ITOP = ITOP-1
      If(Abs(PMom(2)-PMom(3)).LT.thrsh) ITOP = ITOP-1
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize inertia tensor')
c
      END
c ======================================================================
c
      SUBROUTINE NormOUT(IOut,   NAtoms, AtSymb, NVib,   UMode,
     $                   Vib,    Dipole, SIR,    Polar,  SRA,
     $                   SDR,    AtomAx, SVCD,   SLabel, ILabel,
     $                   RLabel, SNDip,  IPRNT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Neat printout of frequencies and normal modes
C
C  ARGUMENTS
C
C  IOut    -  output unit for printout
C  NATOMS  -  number of atoms
C  AtSymb  -  atomic symbols
C  NVib    -  number of normal modes
C  UMode   -  normal modes
C  Vib     -  vibrational frequencies (in cm**-1)
C  Dipole  -  Logical flag   .true.  - dipole derivatives available
C                            .false. - no dipole derivatives
C  SIR     -  IR intensities (in km/mol)
C  Polar   -  Logical flag   .true.  - polarizability derivatives available
C                            .false. - no polarizability derivatives
C  SRA     -  Raman intensities
C  SDR     -  Raman depolarization ratios
C  AtomAx  -  Logical flag   .true.  - atomic axial tensors available
C                            .false. - no atomic axial tensors
C  SVCD    -  VCD rotational strengths
C  SLabel  -  symmetry labels for each mode
C  ILabel  -  IR activity
C               YES - mode is IR active
C               NO  - mode is IR inactive
C  RLabel  -  Raman activity
C               YES - mode is Raman active
C               NO  - mode is Raman inactive
C  SNDip   -  dipole derivatives over normal modes
C  IPRNT   -  print flag for printing of normal modes
C               >0 - print; otherwise do not print
C
C
      REAL*8 Vib(NVib),UMode(3,NAtoms,NVib),SIR(NVib),
     $       SRA(NVib),SDR(NVib),SVCD(NVib),SNDip(3,NVib)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*3 SLabel(NVib)
      CHARACTER*3 ILabel(NVib),RLabel(NVib)
      LOGICAL Dipole,Polar,AtomAx
C
      PARAMETER (maxcol=3)
C
C
      If(IPRNT.GT.0) WRITE(IOut,1000)
      If(IPRNT.LE.0) WRITE(IOut,1001)
c
      NT = NVib/maxcol
      If(NT.EQ.0) GO TO 30
c
      DO 20 I=1,NT
      Imin = (I-1)*maxcol + 1
      Imax = I*maxcol
      WRITE(IOut,1060) (L,L=Imin,Imax)
      WRITE(IOut,1070) (SLabel(L),L=Imin,Imax)
      WRITE(IOut,1100) (Vib(L),L=Imin,Imax)
      WRITE(IOut,1080) (ILabel(L),L=Imin,Imax)
      If(Dipole) Then
        WRITE(IOut,1085) (SIR(L),L=Imin,Imax)
        WRITE(IOut,1086) 'x',(SNDip(1,L),L=Imin,Imax)
        WRITE(IOut,1086) 'y',(SNDip(2,L),L=Imin,Imax)
        WRITE(IOut,1086) 'z',(SNDip(3,L),L=Imin,Imax)
      EndIf
      WRITE(IOut,1090) (RLabel(L),L=Imin,Imax)
      If(Polar) WRITE(IOut,1095) (SRA(L),L=Imin,Imax)
      If(Polar) WRITE(IOut,1096) (SDR(L),L=Imin,Imax)
      If(AtomAx) WRITE(IOut,1099) (SVCD(L),L=Imin,Imax)
      IF(IPRNT.GT.0) THEN
        WRITE(IOut,1097)
        DO 10 K=1,NATOMS
        WRITE(IOut,1200) AtSymb(K),((UMode(J,K,L),J=1,3),L=Imin,IMax)
 10     CONTINUE
      ENDIF
 20   CONTINUE
c
 30   CONTINUE
      NS = NT*maxcol
      NLeft = NVib - NS
      If(NLeft.EQ.0) RETURN
c
      WRITE(IOut,1060) (L,L=NS+1,NVib)
      WRITE(IOut,1070) (SLabel(L),L=NS+1,NVib)
      WRITE(IOut,1100) (Vib(L),L=NS+1,NVib)
      WRITE(IOut,1080) (ILabel(L),L=NS+1,NVib)
      If(Dipole) Then
        WRITE(IOut,1085) (SIR(L),L=NS+1,NVib)
        WRITE(IOut,1086) 'x',(SNDip(1,L),L=NS+1,NVib)
        WRITE(IOut,1086) 'y',(SNDip(2,L),L=NS+1,NVib)
        WRITE(IOut,1086) 'z',(SNDip(3,L),L=NS+1,NVib)
      EndIf
      WRITE(IOut,1090) (RLabel(L),L=NS+1,NVib)
      If(Polar) WRITE(IOut,1095) (SRA(L),L=NS+1,NVib)
      If(Polar) WRITE(IOut,1096) (SDR(L),L=NS+1,NVib)
      If(AtomAx) WRITE(IOut,1099) (SVCD(L),L=NS+1,NVib)
      IF(IPRNT.GT.0) THEN
        WRITE(IOut,1098)           ! warning  only works for maxcol=3  cludge
        DO 40 K=1,NATOMS
        WRITE(IOut,1200) AtSymb(K),((UMode(J,K,L),J=1,3),L=NS+1,NVib)
 40     CONTINUE
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,' ** VIBRATIONAL FREQUENCIES (CM**-1) AND',
     $            ' NORMAL MODES **')
 1001 FORMAT(/,' ** VIBRATIONAL FREQUENCIES (CM**-1) **')
 1060 FORMAT(/,' Label:   ' 10X,I3,10X,2(10X,I3,10X))
 1070 FORMAT(' Symmetry:',11X,A3,10X,2(10X,A3,10X))
 1080 FORMAT(' IR Active:',10X,A3,10X,2(10X,A3,10X))
 1085 FORMAT(' IR Inten:',5X,F9.3,8X,2(6X,F9.3,8X))
 1086 FORMAT(' dmu',A1,'/dQ:',6X,F9.5,8X,2(6X,F9.5,8X))
 1090 FORMAT(' Raman Active:',7X,A3,10X,2(10X,A3,10X))
 1095 FORMAT(' Raman Inten:',2X,F9.3,8X,2(6X,F9.3,8X))
 1096 FORMAT(' Depolar:',8X,F7.3,8X,2(8X,F7.3,8X))
 1097 FORMAT(15X,3('X',6X,'Y',6X,'Z',8X))
 1098 FORMAT(15X,'X',6X,'Y',6X,'Z')
 1099 FORMAT(' Rot.Strength:',1X,F9.3,8X,2(6X,F9.3,8X))
 1100 FORMAT(' Frequency:',5X,F8.2,7X,2(8X,F8.2,7X))
 1200 FORMAT(1X,A8,3(1X,3F7.3,1X))
c
      END
c ======================================================================
c
      SUBROUTINE NormCOND(IOut,   NVib,   Vib,    Dipole, SIR,
     $                    Polar,  SRA,    SDR,    AtomAx, SVCD,
     $                    SLabel, ILabel, RLabel)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Condensed printout of frequencies and intensities
C
C  ARGUMENTS
C
C  IOut    -  output unit for printout
C  NATOMS  -  number of atoms
C  NVib    -  number of normal modes
C  Vib     -  vibrational frequencies (in cm**-1)
C  Dipole  -  Logical flag   .true.  - dipole derivatives available
C                            .false. - no dipole derivatives
C  SIR     -  IR intensities (in km/mol)
C  Polar   -  Logical flag   .true.  - polarizability derivatives available
C                            .false. - no polarizability derivatives
C  SRA     -  Raman intensities
C  SDR     -  Raman depolarization ratios
C  AtomAx  -  Logical flag   .true.  - atomic axial tensors available
C                            .false. - no atomic axial tensors
C  SVCD    -  VCD rotational strengths
C  SLabel  -  symmetry labels for each mode
C  ILabel  -  IR activity
C               YES - mode is IR active
C               NO  - mode is IR inactive
C  RLabel  -  Raman activity
C               YES - mode is Raman active
C               NO  - mode is Raman inactive
C
      REAL*8 Vib(NVib),SIR(NVib),SRA(NVib),SDR(NVib),SVCD(NVib)
      CHARACTER*3 SLabel(NVib),ILabel(NVib),RLabel(NVib)
      LOGICAL Dipole,Polar,AtomAx
      Character*42 tpr
C
      PARAMETER (maxcol=3)
C
C
      write(iout,*) 'Vibrational frequencies and intensities'
      write(iout,*)
      write(iout,*)
     1' No. Symm. Freq.     IR  Raman  IR int. Raman int.  Depol.',
     2'   Rot.Str.'
c
      do i=NVib,1,-1
        call blankit(tpr,42)
        if(Dipole) write(tpr(1:),'(f10.3)') sir(i)
        if(Polar) write(tpr(11:),'(f11.3,f8.4)') sra(i),sdr(i)
        if(AtomAx) write(tpr(31:),'(f11.3)') svcd(i)
        write(IOut,1030) i,slabel(i),vib(i),ilabel(i),rlabel(i),tpr
cc        if(Dipole.and.Polar) then
cc          write(IOut,1010) i,slabel(i),vib(i),ilabel(i),rlabel(i),
cc     1      sir(i), sra(i), sdr(i)
cc        else if(Dipole.and..not.Polar) then
cc          write(IOut,1010) i,slabel(i),vib(i),ilabel(i),rlabel(i),
cc     1      sir(i)
cc        else if(.not.Dipole.and.Polar) then
cc          write(IOut,1020) i,slabel(i),vib(i),ilabel(i),rlabel(i),
cc     1      sir(i), sra(i), sdr(i)
cc        else
cc          write(IOut,1010) i,slabel(i),vib(i),ilabel(i),rlabel(i)
cc        end if
      end do
cc 1010 format(i4,3x,a3,f9.2,2x,a3,2x,a3,2f10.3,f8.4)
cc 1020 format(i4,3x,a3,f9.2,2x,a3,2x,a3,10x,f10.3,f8.4)
 1030 format(i4,3x,a3,f9.2,2x,a3,2x,a3,a)
c
      END
c ======================================================================
c
      SUBROUTINE OrderP(GROUP,nn)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER csf,cla,zn(0:9),GROUP*4
      data zn /'0','1','2','3','4','5','6','7','8','9'/
C
C  nn is the symmetry number of the principal axis
C
      nn = 0
      csf = GROUP(1:1)
c
      DO 5 I=1,9
      If(GROUP(2:2).EQ.zn(I)) nn=I
    5 CONTINUE
c
      IF(nn.EQ.0) THEN
cc
        cla = GROUP(2:2)
        If( csf.EQ.'c'.AND.(cla.EQ.'s'.OR.
     $                    (nn.EQ.0.AND.cla.EQ.'h')) ) Then
         nn=1
        Else If(csf.EQ.'c'.AND.cla.EQ.'i') Then
         nn=2
        EndIf
cc
      ELSE
cc
        np=-1
        DO 10 I=0,9
        If(GROUP(3:3).EQ.zn(i)) np=I
   10   CONTINUE
        If(np.GE.0) nn = 10*nn + np
c
      ENDIF
c
      RETURN
      END
c ======================================================================
c
      SUBROUTINE PrntCHAR(NTrans,NSym,CHARAC,LType)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CHARAC(NTrans,NSym)
      CHARACTER*4 LType(NSym)
c
      WRITE(6,*) ' Irrep Labels Are:'
      WRITE(6,1000) (LType(I),I=1,NSym)
      WRITE(6,*) ' Character Table is:'
      do i=1,nsym
      write(6,2000) (charac(j,i),j=1,NTrans)
      enddo
c
      return
c
 1000 format(16A5)
 2000 format(16F5.2)
      end
c ======================================================================
c
      SUBROUTINE ProjTRM(NATOMS, IPRNT,  XC,     ATMASS,  P,
     $                   TRVec,  HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Projects out from the Mass-Weighted Hessian matrix in Cartesian
C  coordinates vectors corresponding to translations and
C  infinitesimal rotations
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  IPRNT   -  print flag
C  XC      -  Cartesian coordinates
C  ATMASS  -  atomic masses
C  P       -  scratch space for projection matrix
C  TRVec   -  scratch space for TR vectors
C  HESS    -  on entry contains Cartesian Hessian matrix
C             on exit contains projected Hessian
C
      REAL*8 XC(3,NATOMS),ATMASS(NATOMS),P(3*NATOMS,*),
     $       TRVec(3*NATOMS,*),HESS(3*NATOMS,3*NATOMS)
      REAL*8 T(9),PMom(3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
C  Transform to centre-of-mass coordinate system
C
      CALL CoM(NATOMS,ATMASS,XC,CX,CY,CZ)
C
C  Find the principal moments & rotation generators
C
      CALL ZeroIT(T,9)
      DO 20 I=1,NATOMS
      X = XC(1,I)
      Y = XC(2,I)
      Z = XC(3,I)
      ami = ATMASS(I)
      T(1) = T(1) + ami*(Y*Y + Z*Z)
      T(5) = T(5) + ami*(X*X + Z*Z)
      T(9) = T(9) + ami*(X*X + Y*Y)
      T(2) = T(2) - ami*X*Y
      T(3) = T(3) - ami*X*Z
      T(6) = T(6) - ami*Y*Z
 20   CONTINUE
      T(4) = T(2)
      T(7) = T(3)
      T(8) = T(6)
C
C  Diagonalize T
C
      CALL DIAGMAT(T,3,TRVec,P,PMom,IErr)    ! TRVec & P used as scratch
c
      IF(IErr.NE.0) THEN
       WRITE(6,1000)
       CALL OptEXIT(9)
      ENDIF
C
C  Set up Orthogonal coordinate vectors for translation and
C  rotation about principal axes of inertia
C
      NAT3 = 3*NATOMS
cc      CALL ZeroIT(TRVec,NAT3*6)
      CALL ZeroIT(TRVec,NAT3*NAT3)   ! zero out more due to bug in <DGemm>
      CALL FormTRM(NATOMS,ATMASS,XC,T,TRVec)
c
      IF(IPRNT.GT.5) THEN
       WRITE(6,1100)
       CALL PrntMAT(6,NAT3,6,TRVec)
      ENDIF
C
C  Now form the Projection Matrix
C
      CALL ZeroIT(P,NAT3*NAT3)
c
      DO 30 K=1,6
      DO 30 J=1,NAT3
      DO 30 I=1,NAT3
      P(I,J) = P(I,J) - TRVec(I,K)*TRVec(J,K)
 30   CONTINUE
      DO 40 I=1,NAT3
      P(I,I) = One + P(I,I)
 40   CONTINUE
c
      IF(IPRNT.GT.5) THEN
       WRITE(6,1200)
       CALL PrntMat(NAT3,NAT3,NAT3,P)
      ENDIF
C
C  Project out the translations/rotations from Hessian
C     HESS = P * HESS * P(t)
C
      CALL DGemm('N',    'N',    NAT3,   NAT3,   NAT3,
     $            One,    HESS,  NAT3,   P,      NAT3,
     $            Zero,   TRVec, NAT3)
      CALL DGemm('N',    'N',    NAT3,   NAT3,   NAT3,
     $            One,    P,     NAT3,   TRVec,  NAT3,
     $            Zero,   HESS,  NAT3)
cc      CALL MatAB(NAT3,NAT3,NAT3,HESS,P,TRVec,0)
cc      CALL MatAB(NAT3,NAT3,NAT3,P,TRVec,HESS,0)
c
      If(IPRNT.GT.1) WRITE(6,1300)
C
C  Restore original coordinates
C
      DO 50 I=1,NATOMS
      XC(1,I) = XC(1,I) + CX
      XC(2,I) = XC(2,I) + CY
      XC(3,I) = XC(3,I) + CZ
 50   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unable to Diagonalize inertia tensor')
 1100 FORMAT(/,' Vectors for Translations and Rotations')
 1200 FORMAT(/,' The Projection Matrix')
 1300 FORMAT(/,' Translations and Rotations Projected Out of Hessian')
c
      END

*Deck formtrm
      SUBROUTINE FormTRM(NATOMS,ATMASS,XC,T,V)
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8 ATMASS(NATOMS),XC(3,NATOMS),T(9),V(3,NATOMS,6)
C
      PARAMETER (One=1.0d0,TollZERO=1.0d-8)
C
C
C  This routine generates vectors corresponding to translations
C  and infinitesimal rotations given the coordinates (in centre
C  of mass frame) and the eigenvectors of the inertia tensor
C
C
      NAT3 = 3*NATOMS
c
      DO 10 I=1,NATOMS
      X = XC(1,I)
      Y = XC(2,I)
      Z = XC(3,I)
      ami = SQRT(ATMASS(I))
      CX = ami*(X*T(1) + Y*T(2) + Z*T(3))
      CY = ami*(X*T(4) + Y*T(5) + Z*T(6))
      CZ = ami*(X*T(7) + Y*T(8) + Z*T(9))
      V(1,I,1) = ami
      V(2,I,2) = ami
      V(3,I,3) = ami
      V(1,I,4) = CY*T(7) - CZ*T(4)
      V(2,I,4) = CY*T(8) - CZ*T(5)
      V(3,I,4) = CY*T(9) - CZ*T(6)
      V(1,I,5) = CZ*T(1) - CX*T(7)
      V(2,I,5) = CZ*T(2) - CX*T(8)
      V(3,I,5) = CZ*T(3) - CX*T(9)
      V(1,I,6) = CX*T(4) - CY*T(1)
      V(2,I,6) = CX*T(5) - CY*T(2)
      V(3,I,6) = CX*T(6) - CY*T(3)
 10   CONTINUE
c
      DO 20 I=1,6
      skal = SProd(NAT3,V(1,1,I),V(1,1,I))
      IF(skal.GT.TollZERO) THEN
       skal = One/SQRT(skal)
       CALL VScal(NAT3,skal,V(1,1,I))
      ENDIF
 20   CONTINUE
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE RAMAN(NAT3, UMode, PolD, SRA, SDR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transform Cartesian polarizability derivatives to normal mode
C  and get Raman intensities and depolarization ratios
C
C  ARGUMENTS
C
C  NAT3    -  3 x number of atoms
C  UMode   -  normal mode
C  PolD    -  polarizability derivatives
C  SRA     -  Raman intensity
C  SDR     -  dipolarization ratio
C
C
      REAL*8 UMode(NAT3),PolD(NAT3,6)
C
      PARAMETER (Zero=0.0d0,PTol=1.0d-5)
      PARAMETER (FCnvrt=0.078415936d0)   ! au -->  angstroms**4/amu
C
C
C  Transform polarizability derivatives
C
      VXX = Zero
      VXY = Zero
      VXZ = Zero
      VYY = Zero
      VYZ = Zero
      VZZ = Zero
      DO 10 I=1,NAT3
      VXX = VXX + PolD(I,1)*UMode(I)
      VXY = VXY + PolD(I,2)*UMode(I)
      VYY = VYY + PolD(I,3)*UMode(I)
      VXZ = VXZ + PolD(I,4)*UMode(I)
      VYZ = VYZ + PolD(I,5)*UMode(I)
      VZZ = VZZ + PolD(I,6)*UMode(I)
 10   CONTINUE
C
C  now get intensities
C
      alpha = (VXX+VYY+VZZ)/3.0d0
      gamm2 = 0.5d0*( (VXX-VYY)**2 + (VXX-VZZ)**2 + (VYY-VZZ)**2
     $             +   6.0d0*(VXY*VXY + VXZ*VXZ + VYZ*VYZ) )
c
      RInt = 45.0d0*alpha*alpha + 7.0d0*gamm2
c
      perp = 3.0d0*gamm2
      para = RInt - perp
c
      SDR = Zero
      If(para.GT.PTol) SDR = perp/para
      SRA = RInt*FCnvrt
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE SymNUM(GROUP, ISYM)
      INTEGER ISYM,nn
      CHARACTER GROUP*4
C
C  get rotational symmetry number
C
      IF(GROUP.EQ.'c*v ') THEN
       ISYM = 1
       RETURN
      ELSE IF(GROUP.EQ.'d*h ') THEN
       ISYM = 2
       RETURN
      ELSE
C
C  get order of principal axis
C
       CALL OrderP(GROUP,nn)
c
       IF(GROUP(1:1).EQ.'c') THEN
        ISYM = nn
       ELSE IF(GROUP(1:1).EQ.'d') THEN
        ISYM = 2*nn
       ELSE IF(GROUP(1:1).EQ.'s') THEN
        ISYM = nn
       ELSE IF(GROUP(1:1).EQ.'t') THEN
        ISYM = 12
       ELSE IF(GROUP(1:1).EQ.'o') THEN
        ISYM = 24
       ELSE IF(GROUP(1:1).EQ.'i') THEN
        ISYM = 60
       ELSE
        ISYM = 1
       ENDIF
      ENDIF
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE THERMODYN(NAtoms, GROUP,  XC,     AtSymb, AtMASS,
     $                     Temp,   Tend,   Tstep,  Pres,   Pend,
     $                     Pstep,  NEG,    NVib,   Vib,    SCR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Determines classical thermodynamic properties following
C  a vibrational analysis
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  GROUP   -  molecular point group
C  XC      -  Cartesian coordinates
C  AtSymb  -  atomic symbols
C  AtMASS  -  atomic masses
C  Temp    -  start temperature for thermodynamic analysis (K)
C  Tend    -  end temperature for thermodynamic analysis (may be zero)
C  Tstep   -  temperature step (may be zero)
C  Pres    -  start pressure for thermodynamic analysis (atm)
C  Pend    -  end pressure for thermodynamic analysis (may be zero)
C  Pstep    - pressure step (may be zero)
C  NEG     -  number of negative frequencies
C  NVib    -  total number of frequencies
C  Vib     -  frequencies (cm**-1)
C  SCR     -  general scratch storage
C
C
      REAL*8 XC(3,NAtoms),AtMASS(NAtoms),Vib(NVib)
      REAL*8 SCR(30)
      CHARACTER*8 AtSymb(NAtoms)
      CHARACTER*4 GROUP
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0,THalf=1.5d0)
      PARAMETER (ATMTOP=1.01317d+5)
C
C
      call getival('iout',IOut)
      call getival('icon',ICon)
c
      AUTOAN = One/ANTOAU
      h = 2.0d0*PI*hbar
C
C  Calculate ZPVE (in kcal/mol)
C
      ZPVE = Zero
      DO 10 I=NEG+1,NVib
      ZPVE = ZPVE + Vib(I)
 10   CONTINUE
      ZPVE = Half*h*c*avogad*ZPVE/caljou
c
      T = Temp
      P = Pres
 90   CONTINUE
c
      WRITE(IOut,1000) T,P
      WRITE(ICon,1000) T,P
c
      WRITE(IOut,1050) NEG
      WRITE(IOut,1100) ZPVE
c
      Pr = P*ATMTOP
      RT = R*T/caljou
C
C  Calculate vibrational contribution to the enthalpy and entropy
C
      hovrkt = h/(boltz*T)
c
      hvib = Zero
      svib = Zero
      DO 20 I=NEG+1,NVib
      FRQSI = Vib(I)*c
      UI = FRQSI*hovrkt
      EXPUI = EXP(UI)
      hvib = hvib + FRQSI/(EXPUI-One)
      svib = svib + UI/(EXPUI-One) - LOG(One - One/EXPUI)
 20   CONTINUE
c
      svib = 1000d0*svib*R/caljou
      hvib = hvib*h*avogad/caljou + ZPVE
C
C  determine total mass
C
      TOTMAS = Zero
      DO 30 I=1,NAtoms
      TOTMAS = TOTMAS + AtMASS(I)
 30   CONTINUE
C
C  Calculate translational contribution to the enthalpy and entropy
C
      strn = R*( LOG( ((2.0d0*PI*(TOTMAS*amu/h)*(boltz/h)*T)**THalf)
     $                  * (boltz*T/Pr) ) + 5.0d0*Half)
      strn = 1000d0*strn/caljou
      htrn = THalf*R*T/caljou
C
C  Calculate rotational contribution to the enthalpy and entropy
C
      If(GROUP.EQ.'c*v '.OR.GROUP.EQ.'d*h ') Then
       hrot = R*T/caljou
      Else
       hrot = THalf*R*T/caljou
      EndIf
c
      CALL MomINRT(NAtoms, XC,     AtMASS, TOTMAS, SCR(22),
     $             SCR(13),SCR(4), SCR(1), ITOP)
c
      VCONST = 8.0d0*PI*PI*(boltz/h)*T*(amu/h)*AUTOAN*AUTOAN*1.0d-20
      va = VCONST*SCR(1)
      vb = VCONST*SCR(2)
      vc = VCONST*SCR(3)
c
      CALL SymNUM(GROUP,ISYM)
c
      If(GROUP.EQ.'c*v '.OR.GROUP.EQ.'d*h ') THEN
       srot = LOG( vc/Float(ISYM) )
       srot = srot + One
      Else
       srot = LOG( SQRT(PI*va*vb*vc)/Float(ISYM) )
       srot = srot + THalf
      EndIf
      srot = 1000d0*R*srot/caljou
c
      htot = htrn + hrot + hvib + RT
      stot = strn + srot + svib
C
C  .........................................................
C
c -- standard output
      DO 40 I=1,NATOMS
      WRITE(IOut,1200) I,AtSymb(I),AtMASS(I)
 40   CONTINUE
      WRITE(IOut,1300) TOTMAS
      WRITE(IOut,1400)
      WRITE(IOut,1410) SCR(1),SCR(2),SCR(3)
      WRITE(IOut,1420) 'X',SCR(4),SCR(7),SCR(10)
      WRITE(IOut,1420) 'Y',SCR(5),SCR(8),SCR(11)
      WRITE(IOut,1420) 'Z',SCR(6),SCR(9),SCR(12)
      WRITE(IOut,1500) ISYM
      If(ITOP.EQ.1) WRITE(IOut,1510)
      If(ITOP.EQ.2) WRITE(IOut,1520)
      If(ITOP.EQ.3) WRITE(IOut,1530)
      WRITE(IOut,1600) htrn
      WRITE(IOut,1610) hrot
      WRITE(IOut,1620) hvib
      WRITE(IOut,1625) RT
      WRITE(IOut,1630) strn
      WRITE(IOut,1640) srot
      WRITE(IOut,1650) svib
      WRITE(IOut,1660) htot
      WRITE(IOut,1670) stot
c
c -- summary output
      WRITE(ICon,1050) NEG
      WRITE(ICon,1100) ZPVE
      DO 41 I=1,NATOMS
      WRITE(ICon,1200) I,AtSymb(I),AtMASS(I)
 41   CONTINUE
      WRITE(ICon,1300) TOTMAS
      WRITE(ICon,1400)
      WRITE(ICon,1410) SCR(1),SCR(2),SCR(3)
      WRITE(ICon,1420) 'X',SCR(4),SCR(7),SCR(10)
      WRITE(ICon,1420) 'Y',SCR(5),SCR(8),SCR(11)
      WRITE(ICon,1420) 'Z',SCR(6),SCR(9),SCR(12)
      WRITE(ICon,1500) ISYM
      If(ITOP.EQ.1) WRITE(ICon,1510)
      If(ITOP.EQ.2) WRITE(ICon,1520)
      If(ITOP.EQ.3) WRITE(ICon,1530)
      WRITE(ICon,1600) htrn
      WRITE(ICon,1610) hrot
      WRITE(ICon,1620) hvib
      WRITE(ICon,1625) RT
      WRITE(ICon,1630) strn
      WRITE(ICon,1640) srot
      WRITE(ICon,1650) svib
      WRITE(ICon,1660) htot
      WRITE(ICon,1670) stot
c
c -- any more?
      If(P.GE.Pend) GO TO 95
      P = P + Pstep
      GO TO 90
c
 95   CONTINUE
      P = Pres
      If(T.GE.Tend) RETURN
      T = T + Tstep
      GO TO 90
c
 1000 FORMAT(/,' STANDARD THERMODYNAMIC QUANTITIES AT ',F9.3,' K ',
     $         ' AND ',F9.3,' ATM',/)
 1050 FORMAT(2X,' This Molecule has ',I2,' Imaginary Frequencies')
 1100 FORMAT(2X,' Zero point vibrational energy: ',F12.3,
     $            ' kcal/mol',/)
 1200 FORMAT(2X,' Atom ',I4, '  Element ',A2,'  Has Mass ',F12.6)
 1300 FORMAT(2X,' Molecular Mass: ',F12.6,' amu')
 1400 FORMAT(2X,' Principal axes and moments of inertia in',
     $          ' atomic units:',/,
     $       2X,'                          1           2           3')
 1410 FORMAT(2X,'  Eigenvalues --   ',3F12.5)
 1420 FORMAT(10X,A1,10X,3F12.5)
 1500 FORMAT(2X,' Rotational Symmetry Number is ',I3)
 1510 FORMAT(2X,' The Molecule is a Spherical Top')
 1520 FORMAT(2X,' The Molecule is a Symmetric Top')
 1530 FORMAT(2X,' The Molecule is an Asymmetric Top')
 1600 FORMAT(2X,' Translational Enthalpy: ',F12.3,' kcal/mol')
 1610 FORMAT(2X,' Rotational Enthalpy:    ',F12.3,' kcal/mol')
 1620 FORMAT(2X,' Vibrational Enthalpy:   ',F12.3,' kcal/mol')
 1625 FORMAT(2X,' gas constant (RT):      ',F12.3,' kcal/mol')
 1630 FORMAT(2X,' Translational Entropy:  ',F12.3,'  cal/mol.K')
 1640 FORMAT(2X,' Rotational Entropy:     ',F12.3,'  cal/mol.K')
 1650 FORMAT(2X,' Vibrational Entropy:    ',F12.3,'  cal/mol.K',/)
 1660 FORMAT(2X,' Total Enthalpy:         ',F12.3,' kcal/mol')
 1670 FORMAT(2X,' Total Entropy:          ',F12.3,'  cal/mol.K')
c
      END
c ======================================================================
c
      SUBROUTINE VCD(NAT3, UMode, DipD, AAT, SVCD)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transform Cartesian dipole derivatives to normal mode
C  and get IR intensity
C
C  ARGUMENTS
C
C  NAT3    -  3 x number of atoms
C  UMode   -  normal mode
C  DipD    -  dipole derivatives
C  SIR     -  IR intensity
C
C
      REAL*8 UMode(NAT3),DipD(NAT3,3),AAT(nat3,3)
C
      PARAMETER (Zero=0.0d0)
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
C      
C Convert au -->  10**-44 esu**2 cm**2
C
      fcnvrt=1.0d59*hbar**3/(c*amu*xme)
C      
C e**2 a0 hbar / (c amu ) is the factor to the SI unit C**2 m**2
C
C  Transform derivatives
C
      dd1 = Zero
      dd2 = Zero
      dd3 = Zero
      dmd1 = Zero
      dmd2 = Zero
      dmd3 = Zero
      DO 10 I=1,NAT3
      dd1 = dd1 + DipD(I,1)*UMode(I)
      dd2 = dd2 + DipD(I,2)*UMode(I)
      dd3 = dd3 + DipD(I,3)*UMode(I)
      dmd1 = dmd1 + AAT(I,1)*UMode(I)
      dmd2 = dmd2 + AAT(I,2)*UMode(I)
      dmd3 = dmd3 + AAT(I,3)*UMode(I)
 10   CONTINUE
C
C  now get intensities
C
      SVCD = (dd1*dmd1 + dd2*dmd2 + dd3*dmd3)*FCnvrt
C
      RETURN
      END
