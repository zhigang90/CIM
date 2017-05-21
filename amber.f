c ==================================================================
c  AMBER FORCE FIELD           JB   July 2008/January 2012
c ==================================================================
c
      SUBROUTINE AMBER_MAIN(NAtoms, IAN,    IType,  NBend,  NTors,
     $                      XC,     IC,     IAT,    RD,     RLEN,
     $                      B2,     QCH,    CP,     IB1,    IB2,
     $                      IB3,    RANG,   B3,     IT1,    IT2,
     $                      IT3,    IT4,    TORS,   B4,     NT,
     $                      GU,     GC,     HESS,   ffcyc)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Preliminary Driver for AMBER molecular mechanics force field
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  IType   -  calculation type
C               0 - calculate energy + gradient
C               1 -  ditto + Hessian
C              10 - preoptimize + energy + gradient
C              11 -  ditto + Hessian
C  NBend   -  maximum number of bends allowed (dimensioning)
C  NTors   -  maximum number of proper/improper torsions allowed
C  XC      -  Cartesian coordinates
C  IC      -  atomic connectivity matrix
C  IAT     -  AMBER atom type
C  RD      -  distance matrix
C  RLEN    -  AMBER equilibrium bond lengths
C  B2      -  AMBER stretching force constants
C  QCH     -  AMBER atomic charges
C             ** IMPORTANT ** These are defined in the peptide fragments
C  CP      -  AMBER electrostatic factors
C  IB1     -  1st atom in bend
C  IB2     -  2nd atom in bend
C  IB3     -  3rd atom in bend
C  RANG    -  AMBER equilibrium bond angles
C  B3      -  AMBER bending force constants
C  IT1     -  1st atom in proper/improper torsion
C  IT2     -  2nd atom in proper/improper torsion
C  IT3     -  3rd atom in proper/improper torsion
C  IT4     -  4th atom in proper/improper torsion
C  TORS    -  AMBER equilibrium torsions
C  B4      -  AMBER torsional force constants
C  NT      -  integer torsion order
C  GU      -  step-up gradient for finite difference Hessian
C  GC      -  Cartesian gradient (step-down gradient)
C  HESS       Cartesian Hessian matrix
C  ffcyc   -  integer; controls force field initialization
C
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),IC(NAtoms,NAtoms),
     $          IAT(NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms),
     $          B2(NAtoms,NAtoms),IB1(NBend),IB2(NBend),IB3(NBend),
     $          RANG(NBend),B3(NBend),IT1(NTors),IT2(NTors),
     $          IT3(NTors),IT4(NTors),NT(NTors),B4(NTors),
     $          TORS(NTors),CP(NAtoms,NAtoms),QCH(NAtoms)
      DIMENSION GU(3*NAtoms),GC(3,NAtoms),HESS(9*NAtoms*NAtoms)
      COMMON/NUMINT/NumB,NumT,NumP
      INTEGER ffcyc
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
c
      Character jobname*256,cdum*20,wvfnc*20
      Common /job/jobname,lenJ
c
c      PARAMETER (CVRT=0.0015936D0)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
C
      cvrt=rgetrval('kcal/mol')
C
      wvfnc = 'Force Field'
c
C  Read from the <control> file
C     print flag
C     VdW cutoff
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,6,'$print',1,IPRNT,dum,cdum)
      call rdcntrl(IUnit,11,'$vdw_cutoff',2,idum,CutOff,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XC,-1,idum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Attempt to read <ffchk> file
C  If no <ffchk> file is available, this is the first
C  step of a new job
C
C  WARNING: May be QM/MM Job
C  In this case, we will have TWO <ffchk> files, one for
C  the full system and one for the model system
C
C  Check for existence of <qmmm> file
C
      Iqmmm = 0
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.qmmm',
     $      FORM='UNFORMATTED',STATUS='OLD',ERR=90)
      Iqmmm = 1
      READ(40) IQM
      CLOSE (UNIT=40,STATUS='KEEP')
 90   CONTINUE
c
      If(Iqmmm.EQ.1.AND.IQM.EQ.1) Then
        JUnit = 41
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.ffchk1',
     $        FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      Else
        JUnit = 40
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.ffchk',
     $        FORM='UNFORMATTED',STATUS='OLD',ERR=10)
      EndIf
c
      READ(JUnit) IType,IAN,IC,IAT,RD,RLEN,B2,QCH,CP,IB1,IB2,IB3,RANG,
     $            B3,IT1,IT2,IT3,IT4,TORS,B4,NT,NumB,NumT,NumP
      CLOSE (UNIT=JUnit,STATUS='KEEP')
      GO TO 20
C
C  Get atomic numbers from atomic symbols
C
 10   CONTINUE
      CALL GetAtNo(NAtoms,AtSymb,IAN)
c
 20   CONTINUE
C
C  Call Force Field
C  ----------------
      CALL AMBER_Field(NAtoms, IAN,    IType,  NBend,  NTors,
     $                 IPRNT,  XC,     IC,     IAT,    RD,
     $                 RLEN,   B2,     QCH,    CP,     IB1,
     $                 IB2,    IB3,    RANG,   B3,     IT1,
     $                 IT2,    IT3,    IT4,    TORS,   B4,
     $                 NT,     GU,     E,      GC,     HESS,
     $                 CutOff, ffcyc,  IStatus)
C
C
C  Write <ffchk> file
C
      If(JUnit.EQ.40) Then
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.ffchk',
     $        FORM='UNFORMATTED',STATUS='UNKNOWN')
      Else
        OPEN (UNIT=41,FILE=jobname(1:lenJ)//'.ffchk1',
     $        FORM='UNFORMATTED',STATUS='UNKNOWN')
      EndIf
c
      WRITE(JUnit) IType,IAN,IC,IAT,RD,RLEN,B2,QCH,CP,IB1,IB2,IB3,RANG,
     $             B3,IT1,IT2,IT3,IT4,TORS,B4,NT,NumB,NumT,NumP
      CLOSE (UNIT=JUnit,STATUS='KEEP')
C
C  convert energy to atomic units
C  (gradient and Hessian already converted)
C
      E = E*CVRT
c
c ==========================================================
c -- write final energy to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,0,rdum,wvfnc)
      Call wrcntrl(IUnit,7,'$energy',2,0,E,cdum)
c
      If(IType.EQ.1.OR.IType.EQ.11) then
c        write to the <control> file an info about hessian quality
         ihessq=+1
         call  wrcntrl(IUnit,9,'$hessqual',1,ihessq,rdum,cdum)
      Endif
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c -- write forces to <grad> file
      CALL WrGRAD(NAtoms,GC)
c -- write Hessian (if calculated) to <hess> file
      If(IType.EQ.1.OR.IType.EQ.11)
     $   CALL WrHESS(jobname(1:lenJ)//'.hess',lenJ+5,4,3*NAtoms,HESS)
c ==========================================================
C
      RETURN
      END
c  =======================================================================
c
*Deck amber-field
      SUBROUTINE AMBER_Field(NAtoms, IAN,    IType,  NBend,  NTors,
     $                       IPRNT,  XC,     IC,     IAT,    RD,
     $                       RLEN,   B2,     QCH,    CP,     IB1,
     $                       IB2,    IB3,    RANG,   B3,     IT1,
     $                       IT2,    IT3,    IT4,    TORS,   B4,
     $                       NT,     GU,     E,      G,      HESS,
     $                       CutOff, ffcyc,  IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main Driver for AMBER molecular mechanics force field
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  IType   -  calculation type
C               0 - calculate energy + gradient
C               1 -  ditto + Hessian
C              10 - preoptimize + energy + gradient
C              11 -  ditto + Hessian
C  NBend   -  maximum number of bends allowed (dimensioning)
C  NTors   -  maximum number of torsions/out-of-plane bends allowed
C  IPRNT   -  print flag
C  XC      -  Cartesian coordinates
C  IC      -  atomic connectivity matrix
C  IAT     -  AMBER atom type
C  RD      -  distance matrix
C  RLEN    -  AMBER equilibrium bond lengths
C  B2      -  AMBER stretching force constants
C  QCH     -  AMBER atomic charges
C             ** IMPORTANT ** These are defined in the peptide fragments
C  CP      -  AMBER electrostatic factor
C  IB1     -  1st atom in bend
C  IB2     -  2nd atom in bend
C  IB3     -  3rd atom in bend
C  RANG    -  AMBER equilibrium bond angles
C  B3      -  AMBER bending force constants
C  IT1     -  1st atom in proper/impropertorsion
C  IT2     -  2nd atom in proper/impropertorsion
C  IT3     -  3rd atom in proper/impropertorsion
C  IT4     -  4th atom in proper/impropertorsion
C  TORS    -  AMBER equilibrium torsions
C  B4      -  AMBER torsional force constants
C  NT      -  integer torsion order
C  GU      -  gradient for finite difference Hessian
C  E       -  molecular mechanics energy
C  G       -  Cartesian gradient
C  HESS       Cartesian Hessian matrix
C  CutOff  -  distance cutoff to limit Van der Waals interactions
C  ffcyc   -  integer; controls force field initialization
C  IStatus -  status on exit
C              0 - successful calculation
C             -1 - something went wrong
C
C
      DIMENSION IAN(NAtoms),XC(3*NAtoms),IC(NAtoms,NAtoms),
     $          IAT(NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms),
     $          B2(NAtoms,NAtoms),IB1(NBend),IB2(NBend),IB3(NBend),
     $          RANG(NBend),B3(NBend),IT1(NTors),IT2(NTors),
     $          IT3(NTors),IT4(NTors),NT(NTors),B4(NTors),
     $          TORS(NTors),CP(NAtoms,NAtoms),QCH(NAToms)
      DIMENSION GU(3*NAtoms),G(3*NAtoms),HESS(3*NAtoms,3*NAtoms)
      CHARACTER jobname*256
      INTEGER ffcyc
      LOGICAL ReadPQB
      SAVE ReadPQB
c
      PARAMETER (ZERO=0.0D0,delta=0.001D0,HALF=0.5D0)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
c
      PARAMETER (CVRT=0.0015936D0)
      DATA ReadPQB /.False./
c
      Common /job/jobname,lenJ
C
C
C  THIS ROUTINE CAN CALCULATE THE ENERGY, GRADIENT AND HESSIAN
C  of a molecular system using classical model potential
C  functions taken from the AMBER force field.
C  The Gradient is calculated analytically by differentiating
C  the energy with respect to Cartesian displacements and the
C  Hessian is estimated by central differences on the gradient.
C
C  ..............................................................
C     The AMBER force field was coded from the original literature
C
C     W.D.Cornell, P.Cieplak, C.I.Bayly, I.R.Gould, K.M.Merz,
C     D.M.Ferguson, D.C.Spellmeyer, T.Fox, J.W.Caldwell and P.A.Kollman
C     J.Am.Chem.Soc. 117 (1995) 5179 
C ...............................................................
C
C
      AUTOAN = 1.0d0/ANTOAU
c
      IStatus=0
      NAT3 = 3*NAtoms
      CALL VScal(NAT3,AUTOAN,XC)
c
      IF(ffcyc.EQ.0) THEN
C
C  Assign atom connectivity/bond-order data
C  NEW: May be a <pqb> file which will contain both the connectivity
C       and the force field symbols; if so - read it
C
       OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.pqb',
     $       FORM='FORMATTED',STATUS='OLD',ERR=10)
C
C  IF QM/MM JOB
C  We are going to read the <pqb> file twice; the first time for
C  the full system, the second time for the model system.
C  In the latter case, only a subset of the data is required
C  and we have to add the link atom data
C
       If(ReadPQB) Then
c -- read total number of atoms and number of link atoms
        call getival('nqmmm',NQmmm)
        call getival('nlink',NLink)
        Nqm = NQmmm - NAtoms - NLink
        CALL RdPQBFF(40,'Amber',NAtoms,Nqm,IAT,IC)
        CALL AddLink('Amber',NAtoms,NLink,XC,IAT,IC)
       Else
        CALL RdPQBFF(40,'Amber',NAtoms,NAtoms,IAT,IC)
        REWIND 40
        CALL RdQCH(40,NAtoms,QCH)      ! read AMBER atomic charges
        ReadPQB = .True.
       EndIf
c
       CLOSE (UNIT=40,STATUS='KEEP')
       GO TO 15
C
 10    CONTINUE
C
C  If we get here there is NO <pqb> file
C  Must have a <pqb> file for AMBER to define the atomic charges
C
       Call nerror(11,'AMBER Force Field',
     $   'Must Have a <pqb> file for AMBER Force Field',0,0)
C
C  Fill the AMBER interaction and parameter arrays
C
 15    CONTINUE
       CALL FillAMBER(NAtoms, NBend,  NTors,  IPRNT,  XC,
     $                IAT,    IC,     RD,     RLEN,   B2,
     $                QCH,    CP,     IB1,    IB2,    IB3,
     $                RANG,   B3,     IT1,    IT2,    IT3,
     $                IT4,    TORS,   B4,     NT,     IErr)
c
       If(IErr.NE.0) CALL nerror(13,'AMBER Force Field',
     %                'One or more bonding parameters undefined',0,0)
c
       If(IType.GE.10)
     $    CALL AMBER_PREAMBLE(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $                        IC,     RD,     RLEN,   B2,     CP,
     $                        IB1,    IB2,    IB3,    RANG,   B3,
     $                        IT1,    IT2,    IT3,    IT4,    TORS,
     $                        B4,     NT,     G,      GU)
c
      ENDIF
c
      ffcyc = ffcyc + 1
c
      IF(IType.EQ.1.OR.IType.EQ.11) THEN
C
C  form the AMBER Hessian by central difference on the gradient
C
       If(IPRNT.GT.1)
     $    WRITE(6,*) '  Hessian Calculated Using AMBER Force Field'
c
       CALL ZeroIT(HESS,NAT3*NAT3)
c
       DO 20 I=1,NAT3
C
C  step up
C
       XC(I) = XC(I) + delta
c
       CALL EGAMBER(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $              IC,     RD,     RLEN,   B2,     CP,
     $              IB1,    IB2,    IB3,    RANG,   B3,
     $              IT1,    IT2,    IT3,    IT4,    TORS,
     $              B4,     NT,     EU,     GU)
C
C  step down
C
       XC(I) = XC(I) - delta - delta
c
       CALL EGAMBER(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $              IC,     RD,     RLEN,   B2,     CP,
     $              IB1,    IB2,    IB3,    RANG,   B3,
     $              IT1,    IT2,    IT3,    IT4,    TORS,
     $              B4,     NT,     ED,     G)
C
C  restore XC(I) and form row of Hessian
C
       XC(I) = XC(I) + delta
c
       DO 21 J=I,NAT3
       HESS(I,J) = (GU(J)-G(J))/(delta+delta)
 21    CONTINUE
c
 20    CONTINUE
C
C  estimation complete
C  symmetrize the numerical Hessian
C
       DO 30 I=2,NAT3
       DO 30 J=1,I-1
       HESS(I,J) = HALF*(HESS(I,J)+HESS(J,I))
       HESS(J,I) = HESS(I,J)
 30    CONTINUE
c
       If(IPRNT.GT.1) WRITE(6,*) '  Hessian Calculation Complete'
C
C scale the Hessian to get into au
C
       CALL VScal(NAT3*NAT3,CVRT/(ANTOAU**2),HESS)
      ENDIF
c
      IF(IType.GE.0) THEN
C
C  calculate energy and gradient
C
       CALL EGAMBER(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $              IC,     RD,     RLEN,   B2,     CP,
     $              IB1,    IB2,    IB3,    RANG,   B3,
     $              IT1,    IT2,    IT3,    IT4,    TORS,
     $              B4,     NT,     E,      G)
       CALL VScal(NAT3,-CVRT/ANTOAU,G)
      ENDIF
c
 99   CALL VScal(NAT3,ANTOAU,XC)
C
C  -- ERROR HANDLING --
      If(IStatus.LT.0) Call nerror(5,'Force Field module',
     $  'Unable To Fully Assign Atomic Connectivity/Bond Orders',0,0)
C
      RETURN
      END
*Deck eamber
      SUBROUTINE EAMBER(NAtoms, CutOff, IPRNT,  C,      IAT,
     $                  IC,     RD,     RLEN,   B2,     CP,
     $                  IB1,    IB2,    IB3,    RANG,   B3,
     $                  IT1,    IT2,    IT3,    IT4,    TORS,
     $                  B4,     NT,     E)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(3,NAtoms),RD(NAtoms,NAtoms),CP(NAtoms,NAtoms)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      REAL*8 RANG(*),B3(*),B4(*),TORS(*)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),NT(*)
      DIMENSION IB1(*),IB2(*),IB3(*),IT1(*),IT2(*),IT3(*),IT4(*)
      DIMENSION RVdW(40)
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (FOUR=4.0D0,TWELVE=12.0D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
      PARAMETER (TOLLZERO=1.0D-8)
c
      DATA RVDW/1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080,
     $          1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.8240,
     $          1.8240, 1.8240, 1.8240, 1.8240, 1.8240, 1.8750, 1.7683,
     $          1.7210, 1.6837, 1.6612, 1.6612, 2.0000, 2.0000, 2.1000,
     $          0.6000, 0.6000, 0.6000, 0.6000, 1.4590, 1.4870, 1.3870,
     $          1.2870, 1.1870, 1.1000, 1.4090, 1.3590/
C
C
C   THIS ROUTINE CALCULATES THE ENERGY ONLY
C   of a molecule using the AMBER forcefield
C   Energy contributions coded are:
C
C   1. Bond Stretching
C      ES = CS*(R-R0)**2
C   (R0 is an equilibrium bond length and CS a scaling parameter
C   for each bond)
C
C   2. Non-Bonded VdW interaction
C      (1-3 interactions excluded; 1-4 interactions scaled)
C      EV = CV*[ (RV/R)**12 - 2*(RV/R)**6 ]
C   (RV is the sum of the vdw radii of the two atoms and CV is a
C   scaling parameter taken as the square root of the product of
C   the vdw well depth in kcal/mol)
C
C   3. Electrostatic term
C      EE = qi*qj/(e*Rij)
C   (qi and qj are atomic charges on electronegative atoms; e is a scale factor;
C    these factors are collected as CP)
C
C   4. Bending
C      EB = CB*(Th-Th0)**2
C   (Th0 is an equilibrium bond angle and CB a scaling parameter
C   for each bond angle)
C
C   5. Torsion
C      ET = Vn/2 [1 + COS(n*Th - Th0)]
C   (Vn is a torsion potential, n is the periodicity of the torsion and
C   Th0 is an equilibrium torsion)
C
C
      ES = ZERO
      EV = ZERO
      EE = ZERO
      EB = ZERO
      ET = ZERO
C
C  Form the distance matrix
C
      RD(1,1) = ZERO
      DO 10 I=2,NAtoms
      RD(I,I) = ZERO
      DO 10 J=1,I-1
      RD(J,I) = SQRT( (C(1,I)-C(1,J))**2 + (C(2,I)-C(2,J))**2
     $              + (C(3,I)-C(3,J))**2 )
      RD(I,J) = RD(J,I)
 10   CONTINUE
C
C  Go through each interaction type and calculate its
C  contribution to the energy
C
      DO 20 I=2,NAtoms
      I1=IAT(I)
      DO 20 J=1,I-1
      J1=IAT(J)
c
      R = RD(I,J)
c .......................................
      If(R.GT.CutOff) GO TO 20
c .......................................
C
C  check connectivity
C
      IF(IC(I,J).NE.0) THEN
C
C  1. Bond Stretching
C     Atoms I and J are directly bonded
C
        ES = ES + B2(I,J)*(R-RLEN(I,J))**2
c
      ELSE
C
C  2. Van der Waals interaction
C     Atoms I and J are not directly bonded
C     (For 1-3 interactions B2(I,J) has been set to zero in
C      routine FillB2; likewise 1-4 interactions have been scaled)
C
        RV = (RVdW(I1)+RVdW(J1))/R
c
        EV = EV + B2(I,J)*( (RV**12) - TWO*(RV**6) )
c
      ENDIF
C
C  3. Electrostatic term
C     (take another look at this - loop only over charged atoms)
C
      EE = CP(I,J)/R
c
 20   CONTINUE
C
C
C  4. Bending
C     ** Units appear to be radians **
C
      If(NAtoms.LT.3) GO TO 95
c
      DO 30 M=1,NumB
      I=IB1(M)
      J=IB2(M)
      K=IB3(M)
c
      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
c
      CosTh = (RIJ*RIJ+RJK*RJK-RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosTh).GT.ONE) CosTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)
c
      EB = EB + B3(M)*(Th-RANG(M))**2
c
 30   CONTINUE
C
C
C  5. Torsions (proper & improper)
C     ** Units appear to be radians **
C
      If(NAtoms.LT.4) GO TO 95
c
      DO 40 M=1,NumT+NumP
      I=IT1(M)
      J=IT2(M)
      K=IT3(M)
      L=IT4(M)
      n=NT(M)

      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
      RJL=RD(J,L)
      RKL=RD(K,L)
c
      XIJ =C(1,I) - C(1,J)
      XJK =C(1,J) - C(1,K)
      XKL =C(1,K) - C(1,L)
      YIJ =C(2,I) - C(2,J)
      YJK =C(2,J) - C(2,K)
      YKL =C(2,K) - C(2,L)
      ZIJ =C(3,I) - C(3,J)
      ZJK =C(3,J) - C(3,K)
      ZKL =C(3,K) - C(3,L)
c
      R1(1) = -XIJ
      R1(2) = -YIJ
      R1(3) = -ZIJ
      R2(1) = -XJK
      R2(2) = -YJK
      R2(3) = -ZJK
      R3(1) = -XKL
      R3(2) = -YKL
      R3(3) = -ZKL
      CALL VECMUL(R2,R1,R4)
      CALL VECMUL(R3,R2,R5)
c
      CosIJK = (RIJ*RIJ+RJK*RJK-RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosIJK).GT.ONE) CosIJK=SIGN(ONE,CosIJK)
      SinIJK = SQRT(ONE - CosIJK*CosIJK)
      CosJKL = (RJK*RJK+RKL*RKL-RJL*RJL)/(TWO*RJK*RKL)
      If(Abs(CosJKL).GT.ONE) CosJKL=SIGN(ONE,CosJKL)
      SinJKL = SQRT(ONE - CosJKL*CosJKL)
C
C  ...........................................................
C     **  WARNING  **
C  If any three atoms are linear there are problems with the
C  related dihedral angle. Set the angle to zero and skip
C  derivative evaluation
C
      If(SinIJK.LT.TOLLZERO.OR.SinJKL.LT.TOLLZERO) GO TO 40
C  ...........................................................
C
      CosDih = SPROD(3,R4,R5)/(RIJ*RJK*RJK*RKL*SinIJK*SinJKL)
      If(Abs(CosDih).GT.ONE) CosDih=SIGN(ONE,CosDih)
      Dih = ACOS(CosDih)
cc      write(6,*) ' Torsion is ',Dih*degree
c
      ET = ET + B4(M)*(ONE + COS(n*Dih-TORS(M)))
c
 40   CONTINUE
c
 95   CONTINUE
C
C  Form the total molecular energy
C
      E = ES + EV + EE + EB + ET
c
      if(iprnt.gt.2) write(6,1000) E,ES,EB,ET,EE,EV
      RETURN
c
 1000 FORMAT(/,' TOTAL ENERGY       ',f12.8,/,
     $         ' Stretching         ',f12.8,/,
     $         ' Bending            ',f12.8,/,
     $         ' Torsion            ',f12.8,/,
     $         ' Electrostatic      ',f12.8,/,
     $         ' VdW energy         ',f12.8)
c
      END
*Deck egamber
      SUBROUTINE EGAMBER(NAtoms, CutOff, IPRNT,  C,      IAT,
     $                   IC,     RD,     RLEN,   B2,     CP,
     $                   IB1,    IB2,    IB3,    RANG,   B3,
     $                   IT1,    IT2,    IT3,    IT4,    TORS,
     $                   B4,     NT,     E,      G)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(3,NAtoms),RD(NAtoms,NAtoms),CP(NAtoms,NAtoms)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      REAL*8 RANG(*),B3(*),B4(*),TORS(*),G(3,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),NT(*)
      DIMENSION IB1(*),IB2(*),IB3(*),IT1(*),IT2(*),IT3(*),IT4(*)
      DIMENSION RVdW(40)
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (FOUR=4.0D0,TWELVE=12.0D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
      PARAMETER (TOLLZERO=1.0D-8)
c
      DATA RVDW/1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080,
     $          1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.9080, 1.8240,
     $          1.8240, 1.8240, 1.8240, 1.8240, 1.8240, 1.8750, 1.7683,
     $          1.7210, 1.6837, 1.6612, 1.6612, 2.0000, 2.0000, 2.1000,
     $          0.6000, 0.6000, 0.6000, 0.6000, 1.4590, 1.4870, 1.3870,
     $          1.2870, 1.1870, 1.1000, 1.4090, 1.3590/
C
C
C   THIS ROUTINE CALCULATES THE ENERGY AND CARTESIAN GRADIENT
C   of a molecule using the AMBER forcefield
C   Energy contributions coded are:
C
C   1. Bond Stretching
C      ES = CS*(R-R0)**2
C   (R0 is an equilibrium bond length and CS a scaling parameter
C   for each bond)
C
C   2. Non-Bonded VdW interaction
C      (1-3 interactions excluded; 1-4 interactions scaled)
C      EV = CV*[ (RV/R)**12 - 2*(RV/R)**6 ]
C   (RV is the sum of the vdw radii of the two atoms and CV is a
C   scaling parameter taken as the square root of the product of
C   the vdw well depth in kcal/mol)
C
C   3. Electrostatic term
C      EE = qi*qj/(e*Rij)
C   (qi and qj are atomic charges on electronegative atoms; e is a scale factor;
C    these are collected as CP)
C
C   4. Bending
C      EB = CB*(Th-Th0)**2
C   (Th0 is an equilibrium bond angle and CB a scaling parameter
C   for each bond angle)
C
C   5. Torsion
C      ET = Vn/2 [1 + COS(n*Th - Th0)]
C   (Vn is a torsion potential, n is the periodicity of the torsion and
C   Th0 is an equilibrium torsion)
C
C
      ES = ZERO
      EV = ZERO
      EE = ZERO
      EB = ZERO
      ET = ZERO
      CALL ZeroIT(G,3*NAtoms)
C
C  Form the distance matrix
C
      DO 10 J=1,NAtoms
      DO 10 I=1,NAtoms
      RD(I,J) = SQRT( (C(1,I)-C(1,J))**2 + (C(2,I)-C(2,J))**2
     $              + (C(3,I)-C(3,J))**2 )
 10   CONTINUE
C
C  Go through each interaction type and calculate its
C  contribution to the energy and gradient
C
      DO 20 I=2,NAtoms
      I1=IAT(I)
      DO 20 J=1,I-1
      J1=IAT(J)
c
      R = RD(I,J)
c .......................................
      If(R.GT.CutOff) GO TO 20
c .......................................
c
      X12 = C(1,I)-C(1,J)
      Y12 = C(2,I)-C(2,J)
      Z12 = C(3,I)-C(3,J)
C
C  check connectivity
C
      IF(IC(I,J).NE.0) THEN
C
C  1. Bond Stretching
C     Atoms I and J are directly bonded
C
        ES = ES + B2(I,J)*(R-RLEN(I,J))**2
c
        DCS = B2(I,J)*(R-RLEN(I,J))/R
        G(1,I) = G(1,I) + DCS*X12
        G(1,J) = G(1,J) - DCS*X12
        G(2,I) = G(2,I) + DCS*Y12
        G(2,J) = G(2,J) - DCS*Y12
        G(3,I) = G(3,I) + DCS*Z12
        G(3,J) = G(3,J) - DCS*Z12
c
      ELSE
C
C  2. Van der Waals interaction
C     Atoms I and J are not directly bonded
C     (For 1-3 interactions B2(I,J) has been set to zero in
C      routine FillB2; likewise 1-4 interactions have been scaled)
C
        RV = (RVdW(I1)+RVdW(J1))/R
c
        EV = EV + B2(I,J)*( (RV**12) - TWO*(RV**6) )
c
        DCV = TWELVE*B2(I,J)*( -(RV**12) + (RV**6) )/(R*R)
        G(1,I) = G(1,I) + DCV*X12
        G(1,J) = G(1,J) - DCV*X12
        G(2,I) = G(2,I) + DCV*Y12
        G(2,J) = G(2,J) - DCV*Y12
        G(3,I) = G(3,I) + DCV*Z12
        G(3,J) = G(3,J) - DCV*Z12
c
      ENDIF
C
C  3. Electrostatic term
C     (take another look at this - loop only over charged atoms)
C
      EE = EE + CP(I,J)/R
c
      DCE = -EE/R
      G(1,I) = G(1,I) + DCE*X12
      G(1,J) = G(1,J) - DCE*X12
      G(2,I) = G(2,I) + DCE*Y12
      G(2,J) = G(2,J) - DCE*Y12
      G(3,I) = G(3,I) + DCE*Z12
      G(3,J) = G(3,J) - DCE*Z12
c
 20   CONTINUE
cc      write(6,*) ' At 20 CONTINUE'
cc      write(6,*) ' ES is ',es
cc      write(6,*) ' EV is ',ev
cc      write(6,*) ' EE is ',ee
cc      write(6,*) ' gradient is:'
cc      call prntmat(natoms,3,natoms,g)
C
C
C  4. Bending
C     ** Units appear to be radians **
C
      If(NAtoms.LT.3) GO TO 95
c
      DO 30 M=1,NumB
      I=IB1(M)
      J=IB2(M)
      K=IB3(M)
c
      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
c
      CosTh = (RIJ*RIJ+RJK*RJK-RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosTh).GT.ONE) CosTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)
c
      EB = EB + B3(M)*(Th-RANG(M))**2
c
      SinTh = SQRT(ONE - CosTh*CosTh)
C
C  ............................................................
C  If the three atoms are linear there are problems with the
C  angle bend derivative. Skip derivative evaluation
C
      If(SinTh.LT.TOLLZERO) GO TO 30
C  ............................................................
C
      XIJ = C(1,I)-C(1,J)
      XJK = C(1,J)-C(1,K)
      YIJ = C(2,I)-C(2,J)
      YJK = C(2,J)-C(2,K)
      ZIJ = C(3,I)-C(3,J)
      ZJK = C(3,J)-C(3,K)
c
      D1 = TWO*RIJ*RJK
      D2 = RIJ*RIJ + RJK*RJK - RIK*RIK
cc      DCB = HALF*B3(M)*(Th-RANG(M))*(-FOUR*DEGREE/SinTh)/(D1*D1)
      DCB = B3(M)*(Th-RANG(M))*(-TWO/SinTh)/(D1*D1)
      RJKIJ = RJK/RIJ
      RIJJK = RIJ/RJK
c
      G(1,I) = G(1,I) -  DCB*( D1*XJK + D2*XIJ*RJKIJ )
      G(1,J) = G(1,J) +  DCB*( D1*(XJK-XIJ) -
     $                         D2*(XJK*RIJJK - XIJ*RJKIJ) )
      G(1,K) = G(1,K) +  DCB*( D1*XIJ + D2*XJK*RIJJK )
      G(2,I) = G(2,I) -  DCB*( D1*YJK + D2*YIJ*RJKIJ )
      G(2,J) = G(2,J) +  DCB*( D1*(YJK-YIJ) -
     $                         D2*(YJK*RIJJK - YIJ*RJKIJ) )
      G(2,K) = G(2,K) +  DCB*( D1*YIJ + D2*YJK*RIJJK )
      G(3,I) = G(3,I) -  DCB*( D1*ZJK + D2*ZIJ*RJKIJ )
      G(3,J) = G(3,J) +  DCB*( D1*(ZJK-ZIJ) -
     $                         D2*(ZJK*RIJJK - ZIJ*RJKIJ) )
      G(3,K) = G(3,K) +  DCB*( D1*ZIJ + D2*ZJK*RIJJK )
c
 30   CONTINUE
cc      write(6,*) ' At 30 CONTINUE'
cc      write(6,*) ' EB is ',eb
cc      write(6,*) ' gradient is:'
cc      call prntmat(natoms,3,natoms,g)
C
C
C  5. Torsions (proper & improper)
C     ** Units appear to be radians **
C
      If(NAtoms.LT.4) GO TO 95
c
      DO 40 M=1,NumT+NUmP
      I=IT1(M)
      J=IT2(M)
      K=IT3(M)
      L=IT4(M)
      n=NT(M)

      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
      RJL=RD(J,L)
      RKL=RD(K,L)
c
      XIJ =C(1,I) - C(1,J)
      XIK =C(1,I) - C(1,K)
      XJK =C(1,J) - C(1,K)
      XJL =C(1,J) - C(1,L)
      XKL =C(1,K) - C(1,L)
      YIJ =C(2,I) - C(2,J)
      YIK =C(2,I) - C(2,K)
      YJK =C(2,J) - C(2,K)
      YJL =C(2,J) - C(2,L)
      YKL =C(2,K) - C(2,L)
      ZIJ =C(3,I) - C(3,J)
      ZIK =C(3,I) - C(3,K)
      ZJK =C(3,J) - C(3,K)
      ZJL =C(3,J) - C(3,L)
      ZKL =C(3,K) - C(3,L)
c
      R1(1) = -XIJ
      R1(2) = -YIJ
      R1(3) = -ZIJ
      R2(1) = -XJK
      R2(2) = -YJK
      R2(3) = -ZJK
      R3(1) = -XKL
      R3(2) = -YKL
      R3(3) = -ZKL
      CALL VECMUL(R2,R1,R4)
      CALL VECMUL(R3,R2,R5)
c
      CosIJK = (RIJ*RIJ+RJK*RJK-RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosIJK).GT.ONE) CosIJK=SIGN(ONE,CosIJK)
      SinIJK = SQRT(ONE - CosIJK*CosIJK)
      CosJKL = (RJK*RJK+RKL*RKL-RJL*RJL)/(TWO*RJK*RKL)
      If(Abs(CosJKL).GT.ONE) CosJKL=SIGN(ONE,CosJKL)
      SinJKL = SQRT(ONE - CosJKL*CosJKL)
C
C  ...........................................................
C     **  WARNING  **
C  If any three atoms are linear there are problems with the
C  related dihedral angle. Set the angle to zero and skip
C  derivative evaluation
C
      If(SinIJK.LT.TOLLZERO.OR.SinJKL.LT.TOLLZERO) GO TO 40
C  ...........................................................
C
      CosDih = SPROD(3,R4,R5)/(RIJ*RJK*RJK*RKL*SinIJK*SinJKL)
      If(Abs(CosDih).GT.ONE) CosDih=SIGN(ONE,CosDih)
      Dih = ACOS(CosDih)
c
      ET = ET + B4(M)*(ONE + COS(n*Dih-TORS(M)))
C
C  Set terms for torsional derivatives
C
C  First set the derivatives of each term
C  1. The numerator
C      d/dx ( SProd(3,R4,R5) )
C
      DT1XI = -YJK*(XKL*YJK-XJK*YKL) + ZJK*(-XKL*ZJK+XJK*ZKL)
      DT1XJ = -YKL*(XJK*YIJ-XIJ*YJK) + ZKL*(-XJK*ZIJ+XIJ*ZJK)
     $        -ZIK*(XJK*ZKL-XKL*ZJK) + (YIJ+YJK)*(XKL*YJK-XJK*YKL)
      DT1XK = -YIJ*(XKL*YJK-XJK*YKL) + ZIJ*(-XKL*ZJK+XJK*ZKL)
     $        -ZJL*(XIJ*ZJK-XJK*ZIJ) + (YJK+YKL)*(XJK*YIJ-XIJ*YJK)
      DT1XL = -YJK*(XJK*YIJ-XIJ*YJK) + ZJK*(-XJK*ZIJ+XIJ*ZJK)
c
      DT1YI =  XJK*(XKL*YJK-XJK*YKL) - ZJK*(YKL*ZJK-YJK*ZKL)
      DT1YJ =  XKL*(XJK*YIJ-XIJ*YJK) - ZKL*(YJK*ZIJ-YIJ*ZJK)
     $        -XIK*(XKL*YJK-XJK*YKL) + (ZIJ+ZJK)*(YKL*ZJK-YJK*ZKL)
      DT1YK =  XIJ*(XKL*YJK-XJK*YKL) - ZIJ*(YKL*ZJK-YJK*ZKL)
     $        -XJL*(XJK*YIJ-XIJ*YJK) + (ZJK+ZKL)*(YJK*ZIJ-YIJ*ZJK)
      DT1YL =  XJK*(XJK*YIJ-XIJ*YJK) - ZJK*(YJK*ZIJ-YIJ*ZJK)
c
      DT1ZI = -XJK*(-XKL*ZJK+XJK*ZKL) + YJK*(YKL*ZJK-YJK*ZKL)
      DT1ZJ = -XKL*(-XJK*ZIJ+XIJ*ZJK) + YKL*(YJK*ZIJ-YIJ*ZJK)
     $        -YIK*(YKL*ZJK-YJK*ZKL)  + (XIJ+XJK)*(XJK*ZKL-XKL*ZJK)
      DT1ZK = -XIJ*(-XKL*ZJK+XJK*ZKL) + YIJ*(YKL*ZJK-YJK*ZKL)
     $        -YJL*(YJK*ZIJ-YIJ*ZJK)  + (XJK+XKL)*(XIJ*ZJK-XJK*ZIJ)
      DT1ZL = -XJK*(-XJK*ZIJ+XIJ*ZJK) + YJK*(YJK*ZIJ-YIJ*ZJK)
C
C  2. The Denominator
C      d/dx ( RIJ*RJK*RJK*RKL*SinIJK*SinJKL )
C
C  This is more complicated.
C  first consider the derivatives of each term of the product
C    (a)  RIJ
C    (b)  RJK*RJK
C    (c)  RKL
C    (d)  SinIJK
C    (e)  SinJKL
C
      DaXI =  XIJ/RIJ
      DaXJ = -DaXI
      DaYI =  YIJ/RIJ
      DaYJ = -DaYI
      DaZI =  ZIJ/RIJ
      DaZJ = -DaZI
c
      DbXJ =  XJK+XJK
      DbXK = -DbXJ
      DbYJ =  YJK+YJK
      DbYK = -DbYJ
      DbZJ =  ZJK+ZJK
      DbZK = -DbZJ
c
      DcXK =  XKL/RKL
      DcXL = -DcXK
      DcYK =  YKL/RKL
      DcYL = -DcYK
      DcZK =  ZKL/RKL
      DcZL = -DcZK
c
      RIJK =  CosIJK/(RIJ*RJK)
      RIJ2 =  CosIJK*CosIJK/(RIJ*RIJ)
      RJK2 =  CosIJK*CosIJK/(RJK*RJK)
c
      DdXI =  (-(XIJ-XIK)*RIJK + XIJ*RIJ2 )/SinIJK
      DdXJ =  ( (XIJ-XJK)*RIJK - XIJ*RIJ2 + XJK*RJK2 )/SinIJK
      DdXK =  (-(XIK-XJK)*RIJK - XJK*RJK2 )/SinIJK
c
      DdYI =  (-(YIJ-YIK)*RIJK + YIJ*RIJ2 )/SinIJK
      DdYJ =  ( (YIJ-YJK)*RIJK - YIJ*RIJ2 + YJK*RJK2 )/SinIJK
      DdYK =  (-(YIK-YJK)*RIJK - YJK*RJK2 )/SinIJK
c
      DdZI =  (-(ZIJ-ZIK)*RIJK + ZIJ*RIJ2 )/SinIJK
      DdZJ =  ( (ZIJ-ZJK)*RIJK - ZIJ*RIJ2 + ZJK*RJK2 )/SinIJK
      DdZK =  (-(ZIK-ZJK)*RIJK - ZJK*RJK2 )/SinIJK
c
      RJKL =  CosJKL/(RJK*RKL)
      RJK2 =  CosJKL*CosJKL/(RJK*RJK)
      RKL2 =  CosJKL*CosJKL/(RKL*RKL)
c
      DeXJ =  (-(XJK-XJL)*RJKL + XJK*RJK2 )/SinJKL
      DeXK =  ( (XJK-XKL)*RJKL - XJK*RJK2 + XKL*RKL2 )/SinJKL
      DeXL =  (-(XJL-XKL)*RJKL - XKL*RKL2 )/SinJKL
c
      DeYJ =  (-(YJK-YJL)*RJKL + YJK*RJK2 )/SinJKL
      DeYK =  ( (YJK-YKL)*RJKL - YJK*RJK2 + YKL*RKL2 )/SinJKL
      DeYL =  (-(YJL-YKL)*RJKL - YKL*RKL2 )/SinJKL
c
      DeZJ =  (-(ZJK-ZJL)*RJKL + ZJK*RJK2 )/SinJKL
      DeZK =  ( (ZJK-ZKL)*RJKL - ZJK*RJK2 + ZKL*RKL2 )/SinJKL
      DeZL =  (-(ZJL-ZKL)*RJKL - ZKL*RKL2 )/SinJKL
c
      RA = RJK*RJK*RKL*SinJKL
      DT2XI = RA*(RIJ*DdXI + SinIJK*DaXI)
      DT2YI = RA*(RIJ*DdYI + SinIJK*DaYI)
      DT2ZI = RA*(RIJ*DdZI + SinIJK*DaZI)
c
      DT2XK = RIJ*(RA*DdXK + SinIJK*(RKL*RJK*RJK*DeXK
     $              + SinJKL*(RKL*DbXK + RJK*RJK*DcXK)))
      DT2YK = RIJ*(RA*DdYK + SinIJK*(RKL*RJK*RJK*DeYK
     $              + SinJKL*(RKL*DbYK + RJK*RJK*DcYK)))
      DT2ZK = RIJ*(RA*DdZK + SinIJK*(RKL*RJK*RJK*DeZK
     $              + SinJKL*(RKL*DbZK + RJK*RJK*DcZK)))
c
      RA = RIJ*RJK*RJK*SinIJK
      DT2XL = RA*(RKL*DeXL + SinJKL*DcXL)
      DT2YL = RA*(RKL*DeYL + SinJKL*DcYL)
      DT2ZL = RA*(RKL*DeZL + SinJKL*DcZL)
c
      DT2XJ = RKL*(RA*DeXJ + SinJKL*(RIJ*RJK*RJK*DdXJ
     $              + SinIJK*(RIJ*DbXJ + RJK*RJK*DaXJ)))
      DT2YJ = RKL*(RA*DeYJ + SinJKL*(RIJ*RJK*RJK*DdYJ
     $              + SinIJK*(RIJ*DbYJ + RJK*RJK*DaYJ)))
      DT2ZJ = RKL*(RA*DeZJ + SinJKL*(RIJ*RJK*RJK*DdZJ
     $              + SinIJK*(RIJ*DbZJ + RJK*RJK*DaZJ)))
C
C  Now construct the final derivative
C
      A1 = RIJ*RJK*RJK*RKL*SinIJK*SinJKL
      A2 = SProd(3,R4,R5)
c
c -- set DCT
      DCT = ZERO
      DC1 = COS(TORS(M))
      DC2 = SIN(TORS(M))
      CosDih = COS(Dih)
      IF(n.EQ.2) THEN
        DCT = FOUR*CosDih*DC1 +
     $         (TWO*SIN(Dih) - ONE/TAN(Dih))*DC2
      ELSE If(n.EQ.3) THEN
        DCT = (TWELVE*CosDih**2 - 3.0d0)*DC1 +
     $   ( 8.0d0*CosDih*Sin(Dih) - 
     $    (FOUR*CosDih**2 - ONE)/TAN(Dih) )*DC2
      ELSE If(n.EQ.1) THEN
        DCT = COS(TORS(M)) - SIN(TORS(M))/TAN(Dih)
      ELSE IF(n.EQ.4) THEN
        DCT = 8.0d0*(FOUR*CosDih**3 - TWO*CosDih)*DC1 +
     $   ( FOUR*CosDih*(ONE - TWO*CosDih**2)/TAN(Dih) +
     $     FOUR*SIN(Dih)*(6.0d0*CosDih**2 - ONE) )*DC2
      ELSE
c -- should never get here!
        call nerror(1,'egamber','execution error',0,0)
      ENDIF
      DCT = B4(M)*DCT/(A1*A1)
c
      G(1,I) = G(1,I) + DCT*(A1*DT1XI - A2*DT2XI)
      G(1,J) = G(1,J) + DCT*(A1*DT1XJ - A2*DT2XJ)
      G(1,K) = G(1,K) + DCT*(A1*DT1XK - A2*DT2XK)
      G(1,L) = G(1,L) + DCT*(A1*DT1XL - A2*DT2XL)
      G(2,I) = G(2,I) + DCT*(A1*DT1YI - A2*DT2YI)
      G(2,J) = G(2,J) + DCT*(A1*DT1YJ - A2*DT2YJ)
      G(2,K) = G(2,K) + DCT*(A1*DT1YK - A2*DT2YK)
      G(2,L) = G(2,L) + DCT*(A1*DT1YL - A2*DT2YL)
      G(3,I) = G(3,I) + DCT*(A1*DT1ZI - A2*DT2ZI)
      G(3,J) = G(3,J) + DCT*(A1*DT1ZJ - A2*DT2ZJ)
      G(3,K) = G(3,K) + DCT*(A1*DT1ZK - A2*DT2ZK)
      G(3,L) = G(3,L) + DCT*(A1*DT1ZL - A2*DT2ZL)
c
 40   CONTINUE
      write(6,*) ' At 40 CONTINUE'
      write(6,*) ' ET is ',et
      write(6,*) ' gradient is:'
      call prntmat(natoms,3,natoms,g)
C
 95   CONTINUE
C
C  Form the total molecular energy
C
      E = ES + EV + EE + EB + ET
c
      if(iprnt.gt.2) write(6,1000) E,ES,EB,ET,EE,EV
      RETURN
c
 1000 FORMAT(/,' TOTAL ENERGY       ',f12.8,/,
     $         ' Stretching         ',f12.8,/,
     $         ' Bending            ',f12.8,/,
     $         ' Torsion            ',f12.8,/,
     $         ' Electrostatic      ',f12.8,/
     $         ' VdW energy         ',f12.8)
c
      END
*Deck preamble
      SUBROUTINE AMBER_PREAMBLE(NAtoms, CutOff, IPRNT,  C,      IAT,
     $                          IC,     RD,     RLEN,   B2,     CP,
     $                          IB1,    IB2,    IB3,    RANG,   B3,
     $                          IT1,    IT2,    IT3,    IT4,    TORS,
     $                          IB4,    B4,     NT,     G,      CT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(3*NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms)
      REAL*8 B2(NAtoms,NAtoms),RANG(*),B3(*),B4(*),G(3*NAtoms),
     $       CT(3*NAtoms),CP(NAtoms,NAtoms),TORS(*)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),NT(*)
      DIMENSION IB1(*),IB2(*),IB3(*),IT1(*),IT2(*),IT3(*),IT4(*)
      PARAMETER (DMAX=0.2D0)
C
C  This routine does a steepest descent step from the initial
C  geometry together with a crude line search to lower the
C  energy before commencing the optimization proper.
C
C
      NAT3 = 3*NAtoms
      If(IPRNT.GT.0) WRITE(6,*) ' Starting with Steepest Descent'
c
      CALL EGAMBER(NAtoms, CutOff, IPRNT,  C,      IAT,
     $             IC,     RD,     RLEN,   B2,     CP,
     $             IB1,    IB2,    IB3,    RANG,   B3,
     $             IT1,    IT2,    IT3,    IT4,    TORS,
     $             B4,     NT,     E0,     G)
C
C find the gradient norm
C
      GNORM = SQRT(SProd(NAT3,G,G))
C
C  determine the scale factor
C  if GNORM < DMAX take step as predicted
C  else reduce stepsize to DMAX
C
      If(GNORM.GT.DMAX) Then
       SKAL = DMAX/GNORM
       DO 10 I=1,NAT3
       G(I) =  SKAL*G(I)
 10    CONTINUE
      EndIf
c
 100  CONTINUE
C
C  save a copy of the current coordinates
C
      CALL CpyVEC(NAT3,C,CT)
C
C  take the step
C
      DO 20 I=1,NAT3
      C(I) = C(I) - G(I)
 20   CONTINUE
C
C  calculate a new energy
C
      CALL EAMBER(NAtoms, CutOff, IPRNT,  C,      IAT,
     $            IC,     RD,     RLEN,   B2,     CP,
     $            IB1,    IB2,    IB3,    RANG,   B3,
     $            IT1,    IT2,    IT3,    IT4,    TORS,
     $            B4,     NT,     E1)
C
C  if the energy goes down, take another steepest descent step
C  otherwise revert to previous geometry and start Hessian
C  optimization
C
      IF(E1.LT.E0) THEN
       E0 = E1
       GO TO 100
      ELSE
       CALL CpyVEC(NAT3,CT,C)
      ENDIF
c
      If(IPRNT.GT.0) WRITE(6,*) ' End of Steepest Descent'
c
      RETURN
      END
*Deck filldata
      SUBROUTINE FillAMBER(NAtoms, NBend,  NTors,  IPRNT,  C,
     $                     IAT,    IC,     RD,     RLEN,   B2,
     $                     QCH,    CP,     IB1,    IB2,    IB3,
     $                     RANG,   B3,     IT1,    IT2,    IT3,
     $                     IT4,    TORS,   B4,     NT,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(3,NAtoms),IAT(NAtoms),RD(NAtoms,NAtoms)
      dimension rlen(natoms,natoms),b2(natoms,natoms),qch(natoms),
     $          cp(natoms,natoms)
      dimension ib1(*),ib2(*),ib3(*),rang(*),b3(*)
      dimension it1(*),it2(*),it3(*),it4(*),tors(*),b4(*),NT(*)
      common/numint/numb,numt,nump
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
C
C
C  This routine fills the various AMBER parameter arrays from the known atom
C  types, which are (taken from the original AMBER paper):
C
C    1     CT      any carbon sp3
C    2     C       any carbonyl sp2 carbon
C    3     CA      any aromatic sp2 carbon (and C-e of arganine)
C    4     CM      any double-bonded sp2 carbon
C    5     CC      sp2 aromatic in 5-membered ring with 1 substituent + next to N
C    6     CV      sp2 aromatic in 5-membered ring next to C and lone-pair N
C    7     CW      sp2 aromatic in 5-membered ring next to C and N-H
C    8     CR      sp2 aromatic in 5-membered ring next to two N
C    9     CB      sp2 aromatic at junction of 5- and 6-membered rings
C   10     C*      sp2 aromatic in 5-membered ring next to two C
C   11     CN      sp2 junction between 5- and 6-membered rings bonded to C-H & N-H
C   12     CK      sp2 carbon in 5-membered aromatic ring between N and N-R
C   13     CQ      sp2 carbon in 6-membered ring between lone-pair Ns
C   14     N       sp2 nitrogen in amides
C   15     NA      sp2 nitrogen in aromatic rings with H attached
C   16     NB      sp2 nitrogen in 5-membered ring with lone pair
C   17     NC      sp2 nitrogen in 6-membered ring with lone pair
C   18     N*      sp2 nitrogen in 5-membered ring with carbon substituent
C   19     N2      sp2 nitrogen of aromatic amines (and guanidinium ions)
C   20     N3      sp3 nitrogen
C   21     OW      sp3 oxygen in TIP3P water
C   22     OH      sp3 oxygen in alcohols, tyrosine and protonated carboxylic acids
C   23     OS      sp3 oxygen in ethers
C   24     O       sp2 oxygen in amides
C   25     O2      sp2 oxygen in anionic acids
C   26     S       sulphur in methionine and cysteine
C   27     SH      sulphur with hydrogen attached ??
C   28     P       phosphorus in phosphates
C   29     H       hydrogen attached to N
C   30     HW      hydrogen in TIP3P water
C   31     HO      hydrogen in alcohols and acids
C   32     HS      hydrogen attached to S
C   33     HA      hydrogen attached to aromatic C
C   34     HC      hydrogen attached to aliphatic C with no electron-withdrawing group
C   35     H1      hydrogen attached to aliphatic C with 1 electron-withdrawing group
C   36     H2      hydrogen attached to aliphatic C with 2 electron-withdrawing groups
C   37     H3      hydrogen attached to aliphatic C with 3 electron-withdrawing groups
C   38     HP      hydrogen attached to C directly bonded to formally +ve atoms
C   39     H4      hydrogen attached to aromatic C with 1 electronegative neighbor
C   40     H5      hydrogen attached to aromatic C with 2 electronegative neighbors
C   41     F       fluorine attached to sp3 carbon
C
C ** IMPORTANT **
C    All Angles are stored in radians (NOT degrees)
C
C
C Form the distance matrix
C
      DO 10 J=1,NAtoms
      DO 10 I=1,NAtoms
      RD(I,J) = SQRT( (C(1,I)-C(1,J))**2 + (C(2,I)-C(2,J))**2
     $              + (C(3,I)-C(3,J))**2 )
10    CONTINUE
c
      CALL AMBERB2(NAtoms,IAT,IC,RD,RLEN,B2,QCH,CP,IErr)
      If(IErr.NE.0) call nerror(13,'AMBER Force Field',
     $    'One or more bonding parameters undefined',0,0)
      CALL AMBERB3(NAtoms,NBend,IAT,IC,IB1,IB2,IB3,RANG,B3,IErr)
      CALL AMBERB4(NAtoms,NTors,IAT,IC,IT1,IT2,IT3,IT4,TORS,B4,NT,IErr)
C
C  convert the angles & torsions into radians
C
      CALL VScal(numb,1.0d0/DEGREE,RANG)
      CALL VScal(numt+nump,1.0d0/DEGREE,TORS)
c
cc      If(IPRNT.gt.3) THEN
       zero=0.0d0
       write(6,*) ' AMBER Parameter Arrays are:'
       write(6,*) ' Distance Parameters'
       do 11 i=2,natoms
       do 11 j=1,i-1
       if(b2(i,j).ne.zero) write(6,1111) i,j,rlen(i,j),b2(i,j),cp(i,j)
 1111  format(1x,2i4,3f12.6)
 11    continue
       if(numb.gt.0) then
        write(6,*) ' Angle Parameters ',numb
        do 22 i=1,numb
        ang = rang(i)*degree
        write(6,2222) ib1(i),ib2(i),ib3(i),ang,b3(i)
 2222   format(1x,3i4,2f12.6)
 22     continue
       endif
       if(numt.gt.0) then
        write(6,*) ' Torsion Parameters ',numt
        do 33 i=1,numt
        dih = tors(i)*degree
        write(6,3333) it1(i),it2(i),it3(i),it4(i),dih,b4(i),nt(i)
 3333   format(1x,4i4,2f12.6,i6)
 33     continue
       endif
       if(nump.gt.0) then
        write(6,*) ' Improper Torsion Parameters ',nump
        do 44 i=numt+1,numt+nump
        dih = tors(i)*degree
        write(6,3333) it1(i),it2(i),it3(i),it4(i),dih,b4(i),nt(i)
 44     continue
       endif
cc      EndIf
      RETURN
      END
*Deck fillb2
      SUBROUTINE AMBERB2(NAtoms, IAT,    IC,     RD,     RLEN,
     *                   B2,     QCH,    CP,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),QCH(NAtoms),
     $          CP(NAtoms,NAtoms)
      DIMENSION IATI(82),IATJ(82),R0(82),CIJ(82),EH(40)
      LOGICAL Donor,Accept
      PARAMETER (NumDAT=82,ES=1.0d0/1.2d0)
C
c                 CT       C        CA       CM       CC       CV
      DATA EH/0.1094d0,0.0860d0,0.0860d0,0.0860d0,0.0860d0,0.0860d0,
c                 CW       CR       CB       C*       CN       CK
     $        0.0860d0,0.0860d0,0.0860d0,0.0860d0,0.0860d0,0.0860d0,
c                 CQ       N        NA       NB       NC       N*
     $        0.0860d0,0.1700d0,0.1700d0,0.1700d0,0.1700d0,0.1700d0,
c                 N2       N3       OW       OH       OS       O
     $        0.1700d0,0.1700d0,0.1520d0,0.2104d0,0.1700d0,0.2100d0,
c                 O2       S        SH       P        H        HW
     $        0.2100d0,0.2500d0,0.2500d0,0.2000d0,0.0157d0,0.0157d0,
c                 HO       HS       HA       HC       H1       H2
     $        0.0157d0,0.0157d0,0.0150d0,0.0157d0,0.0157d0,0.0157d0,
c                 H3       HP       H4       H5
     $        0.0157d0,0.0157d0,0.0150d0,0.0150d0/
c
      DATA IATI/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,10,10,10,10, 3, 3, 3,
     $           3, 3, 3, 3, 3, 3, 3, 9, 9, 9, 9, 9, 5, 5, 5, 5, 5,12,
     $          12,12, 4, 4, 4, 4, 4, 4,11,13,13, 8, 8, 8, 1, 1, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 6, 7, 7,29,29,29,
     $          29,29,31,31,32,25,22,23,21,26/
c
      DATA IATJ/ 3, 9, 4, 1,14,18,15,17,24,25,22, 9, 1, 7,34, 3, 9, 4,
     $          11, 1,39,33,19,15,17, 9,11,18,16,17, 1, 6, 7,15,16,40,
     $          18,16, 4, 1,39,40,33,18,15,40,17,40,15,16, 1,41,35,36,
     $          37,34,38,14,18,19,20,22,23,26,27,39,16,39,15,14,18,19,
     $          20,15,22,23,27,28,28,28,30,26/
c -- ** WARNING **  JACS paper has bonding data for fluorine (here assumed 41)
c                   but this is not one of the atom types **
c
      DATA R0/1.409d0,1.419d0,1.444d0,1.522d0,1.335d0,1.383d0,1.388d0,
     $        1.358d0,1.229d0,1.250d0,1.364d0,1.459d0,1.495d0,1.352d0,
     $        1.080d0,1.400d0,1.404d0,1.433d0,1.400d0,1.510d0,1.080d0,
     $        1.080d0,1.340d0,1.381d0,1.339d0,1.370d0,1.419d0,1.374d0,
     $        1.391d0,1.354d0,1.504d0,1.375d0,1.371d0,1.385d0,1.394d0,
     $        1.080d0,1.371d0,1.304d0,1.350d0,1.510d0,1.080d0,1.080d0,
     $        1.080d0,1.365d0,1.380d0,1.080d0,1.324d0,1.080d0,1.343d0,
     $        1.335d0,1.526d0,1.380d0,1.090d0,1.090d0,1.090d0,1.090d0,
     $        1.090d0,1.449d0,1.475d0,1.463d0,1.471d0,1.410d0,1.410d0,
     $        1.810d0,1.810d0,1.080d0,1.394d0,1.080d0,1.381d0,1.010d0,
     $        1.010d0,1.010d0,1.010d0,1.010d0,0.960d0,0.960d0,1.336d0,
     $        1.480d0,1.610d0,1.610d0,0.9572d0,2.038d0/
c
      DATA CIJ/469.0d0,447.0d0,410.0d0,317.0d0,490.0d0,424.0d0,418.0d0,
     $         457.0d0,570.0d0,656.0d0,450.0d0,388.0d0,317.0d0,546.0d0,
     $         367.0d0,469.0d0,469.0d0,427.0d0,469.0d0,317.0d0,367.0d0,
     $         367.0d0,481.0d0,427.0d0,483.0d0,520.0d0,447.0d0,436.0d0,
     $         414.0d0,461.0d0,317.0d0,512.0d0,518.0d0,422.0d0,410.0d0,
     $         367.0d0,440.0d0,529.0d0,549.0d0,317.0d0,367.0d0,367.0d0,
     $         367.0d0,448.0d0,428.0d0,367.0d0,502.0d0,367.0d0,477.0d0,
     $         488.0d0,310.0d0,367.0d0,340.0d0,340.0d0,340.0d0,340.0d0,
     $         340.0d0,337.0d0,337.0d0,337.0d0,367.0d0,320.0d0,320.0d0,
     $         227.0d0,237.0d0,367.0d0,410.0d0,367.0d0,427.0d0,434.0d0,
     $         434.0d0,434.0d0,434.0d0,434.0d0,553.0d0,553.0d0,274.0d0,
     $         525.0d0,230.0d0,230.0d0,553.0d0,166.0d0/
C
C
C  This Subroutine fills the RLEN and B2 arrays for bond stretches
C  by going through each data column looking for a match
C  It also sets B2 for any VdW interactions and CP for electrostatics
C
C
      IErr = 0
      CALL ZeroIT(RLEN,NAtoms**2)
      CALL ZeroIT(B2,NAtoms**2)
      CALL ZeroIT(CP,NAtoms**2)
c
      DO 10 I=2,NAtoms
      I1=IAT(I)
      DO 10 J=1,I-1
      J1=IAT(J)
C
C  check connectivity
C
      IF(IC(I,J).NE.0) THEN
C
C  see if there is data for the bonded pair I1,J1
C
       DO 11 K=1,NumDAT
       IF(I1.EQ.IATI(K)) THEN
C
C  match found for I1
C  check if everything else fits
C
        If(J1.EQ.IATJ(K).AND.IC(I,J).GT.0) GO TO 95
       ENDIF
c
       IF(J1.EQ.IATI(K)) THEN
C
C  match found for J1
C  check it everything else fits
C
        If(I1.EQ.IATJ(K).AND.IC(I,J).GT.0) GO TO 95
       ENDIF
 11    CONTINUE
C
C  no match found
C
       WRITE(6,*) ' No match found for bond: ',I,' - ',J
       IErr = -1
       GO TO 10
c
95     CONTINUE
       RLEN(I,J) = R0(K)
       B2(I,J) = CIJ(K)
c
      ELSE
C
C  Van der Waals and electrostatic interactions
C  remove 1,3-interactions
C  reduce 1,4-interactions by half (vdw) or 1.0/1.2 (electrostatic)
C
       DO 20 K=1,NAtoms
       IF(IC(I,K).NE.0) THEN
c -- I is bonded to K
        If(IC(J,K).NE.0) GO TO 10      ! remove 1,3-interaction
        DO 19 L=1,NAtoms
        IF(IC(L,J).NE.0) THEN
c -- J is bonded to L
c -- If L is bonded to K we have a 1,4 interaction   
         If(IC(K,L).NE.0) Then
          B2(I,J) = 0.5d0*SQRT(EH(I1)*EH(J1))
          CP(I,J) = ES*QCH(I)*QCH(J)
          GO TO 10
         EndIf
        ENDIF
 19     CONTINUE
       ENDIF
 20    CONTINUE
C
C set interaction terms
C
       B2(I,J) = SQRT(EH(I1)*EH(J1))
       CP(I,J) = QCH(I)*QCH(J)
      ENDIF
c
10    CONTINUE
c
      RETURN
      END
*Deck fillb3
      SUBROUTINE AMBERB3(NAtoms, NBend,  IAT,    IC,     IB1,
     $                   IB2,    IB3,    RANG,   B3,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RANG(*),B3(*)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IB1(*),IB2(*),IB3(*)
      COMMON/NUMINT/NumB,NumT,NumP
C
C  This Subroutine calculates the number of bends. storing the
C  atoms involved in the arrays IB1, IB2 and IB3, and fills
C  the RANG and B3 arrays by going through the relevant
C  AMBER data columns looking for a match
C
      NumB=0
      If(NAtoms.LT.3) RETURN
      DO 10 I=3,NAtoms
      II=IAT(I)
      DO 10 J=2,I-1
      JJ=IAT(J)
      DO 10 K=1,J-1
      KK=IAT(K)
C
C  determine if I, J and K are connected
C  options are:
C
C  (a) IJ-JK connected
C  (b) IK-KJ connected
C  (c) JI-IK connected
C
      If(IC(I,J).NE.0.AND.IC(J,K).NE.0)
     $   CALL AMMATCH(NAtoms,I,J,K,II,JJ,KK,RD,IB1,IB2,IB3,RANG,B3)
      If(IC(I,K).NE.0.AND.IC(K,J).NE.0)
     $   CALL AMMATCH(NAtoms,I,K,J,II,KK,JJ,RD,IB1,IB2,IB3,RANG,B3)
      If(IC(J,I).NE.0.AND.IC(I,K).NE.0)
     $   CALL AMMATCH(NAtoms,J,I,K,JJ,II,KK,RD,IB1,IB2,IB3,RANG,B3)
c
c -- check number of bends does not exceed maximum allowed
      If(NumB.GT.NBend) Call nerror(2,'Force Field module',
     $           'More Bends in System than Dimensioned for',0,0)
c
 10   CONTINUE
c
      RETURN
      END
*Deck ammatch
      SUBROUTINE AMMATCH(NAtoms,I,J,K,I1,J1,K1,RD,IB1,IB2,IB3,RANG,B3)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RD(NAtoms,NAtoms),RANG(*),B3(*)
      DIMENSION IB1(*),IB2(*),IB3(*)
      DIMENSION IATI(187),IATJ(187),IATK(187),Theta(187),CIJK(187)
      DIMENSION ISTART(40),INUM(40)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ONE=1.0D0,TWO=2.0D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
c
      DATA IATI/ 2, 2, 2, 2, 2, 2,10,10, 3, 3, 5, 5, 4, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,41,41,35,35,
     $          35,35,35,35,35,35,36,36,36,34,38,38,18, 3, 3,
     $           9, 9, 4, 4, 1, 1, 1,14,18,18,18,15,17,24,25,
     $           2, 2, 3, 3, 3, 3, 3, 3, 9, 9, 9, 9, 4, 4,11,
     $          19,19,19,15, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4,
     $          39, 1, 1, 1, 1, 6, 7, 7, 5, 5,39,10,10, 5, 5,
     $          39,40,40,15,15, 2, 2,10,10, 3, 3, 3, 9, 9, 9,
     $          18, 9, 9, 1, 3, 3, 9,40,40,18,40,17, 2, 2, 1,
     $           1,29, 2, 2, 2, 3, 5, 5,11,11, 8, 8, 7, 9, 5,
     $           8, 2, 3, 3, 9, 2, 2, 2, 9, 9, 9,12,12, 4, 4,
     $           3, 3, 1,29, 1,29,30, 2, 1,31, 1, 1,28, 1, 1,
     $           1,32,25,25,25,22,23/
c -- ** WARNING **  JACS paper has bonding data for fluorine (here assumed 41)
c                   but this is not one of the atom types **
c
      DATA IATJ/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
     $           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
     $           3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     $           4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7,
     $           7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
     $           9,10,10,10,11,11,11,12,12,12,13,13,14,14,14,
     $          14,14,15,15,15,15,15,15,15,15,15,15,15,16,16,
     $          16,17,17,17,17,18,18,18,18,18,18,18,18,18,18,
     $          19,19,19,19,20,20,21,22,22,22,23,23,23,26,26,
     $          27,27,28,28,28,28,28/
c
      DATA IATK/ 1,35,34,38,14,20, 1,34, 1,34, 1,34,34, 1,35,
     $          36,34,38,14,18,19,20,22,23,26,27,41,35,35,14,
     $          18,19,22,23,26,27,36,18,23,34,38,20,23, 3,22,
     $          15,24,15,24,14,24,25,24,15,17,24,24,24,24,25,
     $           3,33, 3, 9,11, 1,39,33,39,33,19,17,19,17,33,
     $          19,15,17,17, 4, 1,39,33, 4,39,33, 1,39,33,18,
     $          18, 6, 7,15,16,15,15,16,39,16,16,39,15,39,15,
     $          15,15,16,15,16, 9,16, 3,11, 9,11,16,18,16,17,
     $          17, 1, 7, 7, 9,15,15,18,16,16,17,17, 1,29, 1,
     $          29,29, 2, 3,29,29, 8,29, 7,29, 7,29,29,12, 8,
     $           6, 3, 9,13,13, 4, 1,29,12, 1,29, 1,29, 1,29,
     $           1,29,29,29,29,29,30,31,31,28, 1,28,28, 1,26,
     $          32,32,25,22,23,23,23/
c -- ** WARNING **  JACS paper has bonding data for fluorine (here assumed 41)
c                   but this is not one of the atom types **
c
      DATA Theta/111.10d0,109.50d0,109.50d0,109.50d0,110.10d0,111.20d0,
     $           115.60d0,109.50d0,114.00d0,109.50d0,113.10d0,109.50d0,
     $           109.50d0,109.50d0,109.50d0,109.50d0,109.50d0,109.50d0,
     $           109.70d0,109.50d0,111.20d0,111.20d0,109.50d0,109.50d0,
     $           114.70d0,108.60d0,109.10d0,109.50d0,109.50d0,109.50d0,
     $           109.50d0,109.50d0,109.50d0,109.50d0,109.50d0,109.50d0,
     $           109.50d0,109.50d0,109.50d0,109.50d0,109.50d0,109.50d0,
     $           109.50d0,120.00d0,120.00d0,111.30d0,128.80d0,114.10d0,
     $           125.30d0,116.60d0,120.40d0,117.00d0,122.90d0,115.40d0,
     $           118.60d0,120.90d0,120.60d0,122.50d0,126.00d0,126.00d0,
     $           120.00d0,120.00d0,120.00d0,120.00d0,120.00d0,120.00d0,
     $           120.00d0,120.00d0,120.00d0,120.00d0,123.50d0,117.30d0,
     $           120.10d0,121.50d0,120.00d0,120.00d0,116.00d0,119.30d0,
     $           123.30d0,120.70d0,119.70d0,119.70d0,119.70d0,117.00d0,
     $           123.30d0,123.30d0,119.70d0,119.70d0,119.70d0,121.20d0,
     $           119.10d0,120.00d0,120.00d0,120.00d0,120.00d0,120.00d0,
     $           120.00d0,120.00d0,120.00d0,120.00d0,120.00d0,120.00d0,
     $           108.70d0,120.00d0,120.00d0,120.00d0,120.00d0,120.00d0,
     $           120.00d0,120.00d0,119.20d0,130.00d0,134.90d0,108.80d0,
     $           117.30d0,116.20d0,132.40d0,106.20d0,110.40d0,127.70d0,
     $           126.20d0,128.60d0,106.40d0,125.00d0,122.70d0,132.80d0,
     $           104.40d0,123.05d0,123.05d0,113.90d0,115.45d0,129.10d0,
     $           121.90d0,120.00d0,118.00d0,118.04d0,120.00d0,126.40d0,
     $           125.20d0,116.80d0,118.00d0,120.00d0,120.00d0,111.60d0,
     $           123.10d0,120.00d0,120.00d0,120.00d0,103.80d0,117.00d0,
     $           117.00d0,120.50d0,112.20d0,118.60d0,111.00d0,121.60d0,
     $           117.60d0,119.20d0,105.40d0,125.80d0,125.80d0,128.80d0,
     $           128.80d0,121.20d0,121.20d0,123.20d0,120.00d0,118.40d0,
     $           120.00d0,109.50d0,109.50d0,104.52d0,113.00d0,108.50d0,
     $           108.50d0,109.50d0,120.50d0,120.50d0, 98.90d0,103.70d0,
     $            96.00d0, 92.07d0,119.90d0,108.23d0,108.23d0,102.60d0,
     $           102.60d0/
c
      DATA CIJK/ 63.0d0, 50.0d0, 50.0d0, 50.0d0, 63.0d0, 80.0d0, 63.0d0,
     $           50.0d0, 63.0d0, 50.0d0, 63.0d0, 50.0d0, 50.0d0, 40.0d0,
     $           50.0d0, 50.0d0, 50.0d0, 50.0d0, 80.0d0, 50.0d0, 80.0d0,
     $           80.0d0, 50.0d0, 50.0d0, 50.0d0, 50.0d0, 77.0d0, 35.0d0,
     $           35.0d0, 50.0d0, 50.0d0, 50.0d0, 50.0d0, 50.0d0, 50.0d0,
     $           50.0d0, 35.0d0, 50.0d0, 50.0d0, 35.0d0, 35.0d0, 50.0d0,
     $           50.0d0, 63.0d0, 70.0d0, 70.0d0, 80.0d0, 70.0d0, 80.0d0,
     $           70.0d0, 80.0d0, 70.0d0, 80.0d0, 70.0d0, 70.0d0, 80.0d0,
     $           80.0d0, 80.0d0, 80.0d0, 80.0d0, 63.0d0, 35.0d0, 63.0d0,
     $           63.0d0, 63.0d0, 70.0d0, 35.0d0, 35.0d0, 35.0d0, 35.0d0,
     $           70.0d0, 70.0d0, 70.0d0, 70.0d0, 35.0d0, 70.0d0, 70.0d0,
     $           70.0d0, 70.0d0, 63.0d0, 70.0d0, 35.0d0, 35.0d0, 63.0d0,
     $           35.0d0, 35.0d0, 70.0d0, 35.0d0, 35.0d0, 70.0d0, 35.0d0,
     $           70.0d0, 70.0d0, 70.0d0, 70.0d0, 70.0d0, 70.0d0, 70.0d0,
     $           35.0d0, 70.0d0, 35.0d0, 35.0d0, 70.0d0, 35.0d0, 70.0d0,
     $           35.0d0, 35.0d0, 35.0d0, 70.0d0, 70.0d0, 63.0d0, 70.0d0,
     $           63.0d0, 63.0d0, 63.0d0, 63.0d0, 70.0d0, 70.0d0, 70.0d0,
     $           70.0d0, 70.0d0, 70.0d0, 63.0d0, 70.0d0, 63.0d0, 70.0d0,
     $           70.0d0, 35.0d0, 35.0d0, 70.0d0, 35.0d0, 70.0d0, 50.0d0,
     $           30.0d0, 50.0d0, 30.0d0, 35.0d0, 70.0d0, 70.0d0, 30.0d0,
     $           30.0d0, 70.0d0, 30.0d0, 70.0d0, 30.0d0, 70.0d0, 30.0d0,
     $           30.0d0, 70.0d0, 70.0d0, 70.0d0, 70.0d0, 70.0d0, 70.0d0,
     $           70.0d0, 70.0d0, 70.0d0, 30.0d0, 70.0d0, 70.0d0, 30.0d0,
     $           70.0d0, 30.0d0, 70.0d0, 30.0d0, 50.0d0, 35.0d0, 35.0d0,
     $           35.0d0, 50.0d0, 35.0d0,100.0d0, 35.0d0, 55.0d0, 45.0d0,
     $           60.0d0,100.0d0,100.0d0, 62.0d0, 68.0d0, 43.0d0, 35.0d0,
     $          140.0d0, 45.0d0,100.0d0, 45.0d0, 45.0d0/
c
      DATA ISTART/  1, 44, 61, 80, 92, 99,102,107,111,122,
     $            125,128,131,133,138,149,152,156,166,170,
     $            172,173,176,  0,  0,179,181,183,  0,  0,
     $              0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
      DATA INUM/43,17,19,12, 7, 3, 5, 4,11, 3, 3, 3, 2, 5,11,
     $           3, 4,10, 4, 2, 1, 3, 3, 0, 0, 2, 2, 5, 0, 0,
     $           0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
C
C
      NumB=NumB+1
      IB1(NumB)=I
      IB2(NumB)=J
      IB3(NumB)=K
C
C Check the central atom
C
      IF(INUM(J1).NE.0) THEN
C
C match found for J1
C see if anything else fits
C
       L1 = ISTART(J1)
       L2 = L1+INUM(J1)-1
       DO 10 L=L1,L2
       If( (I1.EQ.IATI(L).AND.K1.EQ.IATK(L)).OR.
     $     (I1.EQ.IATK(L).AND.K1.EQ.IATI(L)) ) GO TO 95
 10    CONTINUE
C
C  no explicit fit found
C
       WRITE(6,*) ' No match found for angle: ',I,' - ',J,' - ',K
       call nerror(13,'AMBER Force Field',
     $    'One or more angle parameters undefined',0,0)
c
      ENDIF
c
 95   CONTINUE
      RANG(NumB) = Theta(L)
      B3(NumB) = CIJK(L)
c
      RETURN
      END
*Deck fillb4
      SUBROUTINE AMBERB4(NAtoms, NTors,  IAT,    IC,     IT1,
     $                   IT2,    IT3,    IT4,    TORS,   B4,
     $                   NT,     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),TORS(*),B4(*),NT(*)
      DIMENSION IOP(6)
      COMMON/NUMINT/NumB,NumT,NumP
C
C  This Subroutine calculates the number of proper and
C  improper torsions, storing the atoms involved in the
C  arrays IT1, IT2, IT3 and IT4, and fills the TORS, NT and B4
C  arrays by going through the relevant AMBER data
C  columns looking for a match
C
      NumT=0
      NumP=0
      If(NAtoms.LT.4) RETURN
C
C -- proper torsions
C
      DO 10 I=4,NAtoms
      I1=IAT(I)
      DO 10 J=3,I-1
      J1=IAT(J)
      DO 10 K=2,J-1
      K1=IAT(K)
      DO 10 L=1,K-1
      L1=IAT(L)
C
C  determine if I, J, K and L form a dihedral angle
C  options are:
C
C  (a)  KI-IJ-JL  connected
C  (b)  KJ-JI-IL  connected
C  (c)  JI-IK-KL  connected
C  (d)  JK-KI-IL  connected
C  (e)  JI-IL-LK  connected
C  (t)  JL-LI-IK  connected
C  (g)  IJ-JK-KL  connected
C  (h)  IK-KJ-JL  connected
C  (i)  IJ JL-LK  connected
C  (j)  IL-LJ-JK  connected
C  (k)  IK-KL-LJ  connected
C  (1)  IL-LK-KJ  connected
C
      IF(IC(I,J).NE.0) THEN
C
C  I-J connected
C
      If(IC(K,I).NE.0.AND.IC(J,L).NE.0)
     $ CALL DMAMBER(NAtoms,K,I,J,L,K1,I1,J1,L1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(K,J).NE.0.AND.IC(I,L).NE.0)
     $ CALL DMAMBER(NAtoms,K,J,I,L,K1,J1,I1,L1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(I,K).NE.0.AND.IC(K,L).NE.0)
     $ CALL DMAMBER(NAtoms,J,I,K,L,J1,I1,K1,L1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(I,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMAMBER(NAtoms,J,I,L,K,J1,I1,L1,K1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(J,K).NE.0.AND.IC(K,L).NE.0)
     $ CALL DMAMBER(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(J,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMAMBER(NAtoms,I,J,L,K,I1,J1,L1,K1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      ENDIF
c
      IF(IC(J,K).NE.0) THEN
C
C  J-K connected
C
      If(IC(K,I).NE.0.AND.IC(I,L).NE.0)
     $ CALL DMAMBER(NAtoms,J,K,I,L,J1,K1,I1,L1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(I,K).NE.0.AND.IC(J,L).NE.0)
     $ CALL DMAMBER(NAtoms,I,K,J,L,I1,K1,J1,L1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(I,L).NE.0.AND.IC(L,J).NE.0)
     $ CALL DMAMBER(NAtoms,I,L,J,K,I1,L1,J1,K1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(I,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMAMBER(NAtoms,I,L,K,J,I1,L1,K1,J1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      ENDIF
c
      IF(IC(I,K).NE.0) THEN
C
C  I-K connected
C
      If(IC(J,L).NE.0.AND.IC(L,I).NE.0)
     $ CALL DMAMBER(NAtoms,J,L,I,K,J1,L1,I1,K1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      If(IC(K,L).NE.0.AND.IC(L,J).NE.0)
     $ CALL DMAMBER(NAtoms,I,K,L,J,I1,K1,L1,J1,IC,IT1,IT2,IT3,IT4,
     $              TORS,NT,B4)
      ENDIF
c
c -- check number of torsions does not exceed maximum allowed
      If(NumT.GT.NTors) Call nerror(3,'Force Field module',
     $           'More Torsions in System than Dimensioned for',0,0)
c
 10   CONTINUE
C
C
C -- improper torsions
C
      DO 20 I=4,NAtoms
      I1=IAT(I)
      DO 20 J=3,I-1
      J1=IAT(J)
      DO 20 K=2,J-1
      K1=IAT(K)
      DO 20 L=1,K-1
      L1=IAT(L)
C
C  determine if I, J, K and L form an improper torsion
C  options are:
C
C  (a)  I central atom  IJ, IK, IL connected
C  (b)  J central atom  JI, JK, JL connected
C  (c)  K central atom  KI, KJ, KL connected
C  (d)  L central atom  LI, LJ, LK connected
C
      IF(IC(I,J).NE.0.AND.IC(I,K).NE.0.AND.IC(I,L).NE.0) THEN
       CALL IMPROPER(NAtoms,J,K,I,L,J1,K1,I1,L1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,K,L,I,J,K1,L1,I1,J1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,L,J,I,K,L1,J1,I1,K1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
      ENDIF
      IF(IC(J,I).NE.0.AND.IC(J,K).NE.0.AND.IC(J,L).NE.0) THEN
       CALL IMPROPER(NAtoms,I,K,J,L,I1,K1,J1,L1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,K,L,J,I,K1,L1,J1,I1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,L,I,J,K,L1,I1,J1,K1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
      ENDIF
      IF(IC(K,I).NE.0.AND.IC(K,J).NE.0.AND.IC(K,L).NE.0) THEN
       CALL IMPROPER(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,J,L,K,I,J1,L1,K1,I1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,L,I,K,J,L1,I1,K1,J1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
      ENDIF
      IF(IC(L,I).NE.0.AND.IC(L,J).NE.0.AND.IC(L,K).NE.0) THEN
       CALL IMPROPER(NAtoms,I,J,L,K,I1,J1,L1,K1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,J,K,L,I,J1,K1,L1,I1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
       CALL IMPROPER(NAtoms,K,I,L,J,K1,I1,L1,J1,IC,IT1,IT2,IT3,IT4,
     $               TORS,NT,B4)
      ENDIF
c
c -- check number of improper torsions does not exceed maximum allowed
      If(NumT+NumP.GT.NTors) Call nerror(4,'Force Field module',
     $  'More Improper Torsions in System than Dimensioned for',0,0)
c
 20   CONTINUE
C
      RETURN
      END
*Deck dmamber
      SUBROUTINE DMAMBER(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,
     $                   IT3,IT4,TORS,NT,B4)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),TORS(*),NT(*),B4(*)
      DIMENSION IATI(88),IATJ(88),IATK(88),IATL(88)
      DIMENSION THETA(88),IBT(88),CIJKL(88),IS(88)
      DIMENSION ISTART(2),INUM(2)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (TolZero=1.0d-6)
c
      DATA IATI/ 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,14,14,14,14,
     $          22,22,23,23,23,23,23,23,23,26, 1, 1,29,29,22,22,23,23,
     $          52*0/
c
      DATA IATJ/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,26,26,14,14,28,28,28,28,
     $           2, 2, 2, 2, 2, 2, 2, 2, 2,10,10,10, 3, 3, 3, 3, 3, 3,
     $           3, 3, 9, 9, 9, 9, 9, 5, 5, 5, 5, 5,12,12, 4, 4, 4,11,
     $          13, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 7,22,23/
c
      DATA IATK/14,14,14,14, 2, 2, 2, 2,14,14,14,14,23,23, 2, 2, 2, 2,
     $           1, 1, 1, 1, 1, 1,18,18,18,18,26,26, 2, 2,23,23,23,23,
     $           3, 9, 4, 1,14,18,15,17,22, 9, 1, 7, 3, 9, 4,11, 1,19,
     $          15,17, 9,11,18,16,17, 1, 6, 7,15,16,18,16, 4, 1,18,15,
     $          17,15,16, 1,14,18,19,20,22,23,26,27,16,15,28,28/
c
      DATA IATL/ 2, 2, 2, 2,14,14,14,14, 2, 2, 2, 2, 1, 1,14,14,14,14,
     $          22,22,22,22,23,23,12,12, 4, 4, 1, 1,24,24, 1, 1, 1, 1,
     $          52*0/
c
      DATA IBT/  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $           4, 4, 4, 4, 4, 4, 4, 2, 2, 4, 6, 4, 4, 4, 4, 4, 6, 4,
     $           4, 2, 4, 4, 4, 2, 2, 6, 4, 4, 4, 2, 4, 2, 4, 6, 4, 4,
     $           2, 4, 2, 9, 6, 6, 6, 9, 3, 3, 3, 3, 2, 4, 3, 3/
c
      DATA THETA/  0.0d0,180.0d0,180.0d0,180.0d0, 0.0d0, 0.0d0, 0.0d0,
     $           180.0d0,180.0d0,180.0d0,180.0d0, 0.0d0, 0.0d0,180.0d0,
     $           180.0d0,  0.0d0,180.0d0,180.0d0, 0.0d0, 0.0d0, 0.0d0,
     $             0.0d0,  0.0d0,  0.0d0,180.0d0,0.0d0,180.0d0, 0.0d0,
     $             0.0d0,  0.0d0,180.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0,
     $             0.0d0,180.0d0,180.0d0,180.0d0,0.0d0,180.0d0,180.0d0,
     $           180.0d0,180.0d0,180.0d0,180.0d0,0.0d0,180.0d0,180.0d0,
     $           180.0d0,180.0d0,180.0d0,0.0d0,180.0d0,180.0d0,180.0d0,
     $           180.0d0,180.0d0,180.0d0,180.0d0,180.0d0,0.0d0,180.0d0,
     $           180.0d0,180.0d0,180.0d0,180.0d0,180.0d0,180.0d0,0.0d0,
     $           180.0d0,180.0d0,180.0d0,180.0d0,180.0d0,0.0d0,  0.0d0,
     $             0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,0.0d0,  0.0d0,
     $           180.0d0,180.0d0,  0.0d0,  0.0d0/
c
      DATA CIJKL/ 0.0d0, 0.0d0, 0.2d0, 0.0d0, 0.1d0, 0.0d0, 0.07d0,
     $            0.0d0, 0.5d0, 0.15d0,0.0d0, 0.53d0,0.383d0, 0.1d0,
     $            0.4d0, 0.0d0, 1.35d0,0.75d0,0.144d0,1.0d0,0.144d0,
     $            1.0d0,0.144d0,1.0d0, 0.5d0, 2.5d0, 0.5d0, 2.5d0,
     $            0.6d0, 3.5d0, 2.5d0, 2.0d0, 0.25d0, 1.2d0,0.25d0,
     $            1.2d0,14.5d0,12.0d0, 8.7d0, 0.0d0,10.0d0, 5.8d0,
     $            5.4d0, 8.0d0, 1.8d0, 6.7d0, 0.0d0,26.1d0,14.5d0,
     $           14.0d0,10.2d0,14.5d0, 0.0d0, 9.6d0, 6.0d0, 9.6d0,
     $           21.8d0,12.0d0, 6.6d0, 5.1d0, 8.3d0, 0.0d0,20.6d0,
     $           21.5d0, 5.6d0, 4.8d0, 6.8d0,20.0d0,26.6d0, 0.0d0,
     $            7.4d0, 6.1d0,13.6d0, 9.3d0,10.0d0, 1.4d0, 0.0d0,
     $            0.0d0, 0.0d0, 1.4d0, 0.5d0,1.15d0, 1.0d0, 0.75d0,
     $            4.8d0, 6.0d0, 0.75d0,0.75d0/

      DATA   IS/ 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1,
     $           3, 2, 3, 2, 3, 2, 2, 1, 2, 1, 3, 2, 2, 1, 3, 2, 3, 2,
     $           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2,
     $           2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3/
c
      DATA ISTART/1,37/, INUM/36,52/
C
C
      NumT=NumT+1
      IT1(NumT)=I
      IT2(NumT)=J
      IT3(NumT)=K
      IT4(NumT)=L
C
C  1. check for explicit fits
C
      M1=ISTART(1)
      M2=M1+INUM(1)-1
      DO 10 M=M1,M2
      IF(J1.EQ.IATJ(M).AND.K1.EQ.IATK(M)) THEN
C
C  match found for central atoms
C  check if everything else fits
C
       If(I1.EQ.IATI(M).AND.L1.EQ.IATL(M).AND.
     $    IC(J,K).EQ.IBT(M)) GO TO 95
      ENDIF
c
      IF(K1.EQ.IATJ(M).AND.J1.EQ.IATK(M)) THEN
C
C  match found for central atoms
C  check if everything else fits
C
       If(L1.EQ.IATI(M).AND.I1.EQ.IATL(M).AND.
     $    IC(K,J).EQ.IBT(M)) Go TO 95
      ENDIF
 10   CONTINUE
C
C  2. No single dummy fits in AMBER
C     check for double dummy fits
C
      M1=ISTART(2)
      M2=M1+INUM(2)-1
      DO 30 M=M1,M2
      If(J1.EQ.IATJ(M).AND.K1.EQ.IATK(M)) GO TO 95
      If(K1.EQ.IATJ(M).AND.J1.EQ.IATK(M)) GO TO 95
 30   CONTINUE
C
C  At this point nothing fits
C
      WRITE(6,*) ' No match found for torsion: ',I,' - ',J,' - ',K,
     $           ' - ',L
      call nerror(13,'AMBER Force Field',
     $    'One or more torsion parameters undefined',0,0)
c
 95   CONTINUE
C
C  throw out any torsions with zero potential
C
      If(CIJKL(M).GT.TolZero) Then
        TORS(NumT) = THETA(M)
        B4(NumT) = CIJKL(M)/IBT(M)
        NT(NumT)= IS(M)
      Else
        NumT = NumT-1
      EndIf
c
      RETURN
      END
*Deck improper
      SUBROUTINE IMPROPER(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,
     $                    IT3,IT4,TORS,NT,B4)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),TORS(*),NT(*),B4(*)
      DIMENSION IATI(30),IATJ(30),IATK(30),IATL(30)
      DIMENSION CIJKL(30)
      DIMENSION ISTART(3),INUM(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      DATA IATI/ 3, 3, 9,12, 4, 4, 1, 7,17,15,15,15,16,17*0/
c
      DATA IATJ/ 3, 3,17, 9, 2, 2, 4, 9, 4, 6, 7,17, 7, 1,19,25,14*0/
c
      DATA IATK/ 2, 3, 3,18, 4,18, 4,10, 3, 5, 5, 3, 5,14, 3, 2, 2, 3,
     $           3, 3,12, 4, 4,13, 8, 6, 7,14,19,15/
c
      DATA IATL/22, 1,19, 1, 1, 1, 2, 1,19, 1, 1,19, 1, 1,19,25,24,39,
     $          40,33,40,39,33,40,40,39,39,29,29,29/
c
      DATA CIJKL/ 1.1d0, 1.1d0, 1.1d0, 1.0d0, 1.1d0, 1.0d0, 1.1d0,
     $            1.1d0, 1.1d0, 1.1d0, 1.1d0, 1.1d0, 1.1d0, 1.0d0,
     $           10.5d0,10.5d0,10.5d0, 1.1d0, 1.1d0, 1.1d0, 1.1d0,
     $            1.1d0, 1.1d0, 1.1d0, 1.1d0, 1.1d0, 1.1d0, 1.0d0,
     $            1.0d0, 1.0d0/
c
      DATA ISTART/1,14,17/, INUM/13,3,14/
C
C
C  The bond angle, Theta, for improper torsions is always 180 degrees
C  and the periodicity is always 2
C
      NumP=NumP+1
      Num=NumP+NumT
      IT1(Num)=I
      IT2(Num)=J
      IT3(Num)=K
      IT4(Num)=L
C
C  1. check for explicit fits
C
      M1=ISTART(1)
      M2=M1+INUM(1)-1
      DO 10 M=M1,M2
      IF(K1.EQ.IATK(M)) THEN
C
C  match found for central atom
C  check if everything else fits
C
       If(I1.EQ.IATI(M).AND.J1.EQ.IATJ(M).AND.
     $    L1.EQ.IATL(M)) GO TO 95
      ENDIF
c
 10   CONTINUE
C
C  2. check for single dummy fits
C
      M1=ISTART(2)
      M2=M1+INUM(2)-1
      DO 20 M=M1,M2
      IF(K1.EQ.IATK(M)) THEN
C
C  match found for central atom
C  check if everything else fits
C
       If(J1.EQ.IATJ(M).AND.L1.EQ.IATL(M)) GO TO 95
      ENDIF
 20   CONTINUE
C
C  3. check for double dummy fits
C
      M1=ISTART(3)
      M2=M1+INUM(3)-1
      DO 30 M=M1,M2
      If(K1.EQ.IATK(M).AND.L1.EQ.IATL(M)) GO TO 95
 30   CONTINUE
C
C  At this point nothing fits
C  We are going to ignore this coordinate
C
      NumP=NumP-1
      RETURN
c
 95   CONTINUE
      TORS(Num) = 180.0d0
      B4(Num) = CIJKL(M)
      NT(Num)= 2
c
      RETURN
      END
*Deck rdqch
      SUBROUTINE RdQCH(IUnit, NAToms, QCH)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Read the AMBER atomic charges from the <pqb> file
C
C  ARGUMENTS
C
C  IUnit   -  unit number of <pqb> file
C  NAtoms  -  number of atoms
C  QCH     -  on exit contains AMBER Atomic charges
C
C  NOTE: The format of the <PQB> file is not yet determined, but is assumed
C  to be    Cartesian coordinates    Force Field atomic symbols    Charges
C             X   Y   Z    atomic   Sybyl   UFF   AMBER   atomic charges
C
      DIMENSION QCH(NAtoms)
      CHARACTER*20 Char
C
      WRITE(6,*) ' AMBER atomic charges read from <pqb> file'
C
C  Skip the first two lines
C
      READ(IUnit,910) Char
      READ(IUnit,910) Char
c
      DO IAtm=1,NAtoms
      READ(IUnit,920) QCH(IAtm)
      EndDO
C
      RETURN
c
 910  Format(A20)
 920  Format(75X,F7.4)
c
      END
