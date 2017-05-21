c ==================================================================
c  SYBYL_5.2 FORCE FIELD           JB   Jan 2000
c  (based on old code originally written in 1992)
c ==================================================================
c
      SUBROUTINE SYBYL_MAIN(NAtoms, IAN,    IType,  NBend,  NTors,
     $                      XC,     IC,     IAT,    RD,     RLEN,
     $                      B2,     IB1,    IB2,    IB3,    RANG,
     $                      B3,     IT1,    IT2,    IT3,    IT4,
     $                      IS4,    B4,     GU,     GC,     HESS,
     $                      ffcyc)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Preliminary Driver for Sybyl v.5.2 molecular mechanics force field
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
C  XC      -  Cartesian coordinates
C  IC      -  atomic connectivity matrix
C  IAT     -  Sybyl atom type
C  RD      -  distance matrix
C  RLEN    -  Sybyl equilibrium bond lengths
C  B2      -  Sybyl stretching force constants
C  IB1     -  1st atom in bend
C  IB2     -  2nd atom in bend
C  IB3     -  3rd atom in bend
C  RANG    -  Sybyl equilibrium bond angles
C  B3      -  Sybyl bending force constants
C  IT1     -  1st atom in torsion/out-of-plane bend
C  IT2     -  2nd atom in torsion/out-of-plane bend
C  IT3     -  3rd atom in torsion/out-of-plane bend
C  IT4     -  4th atom in torsion/out-of-plane bend
C  IS4     -  Sybyl sign factor for torsion/out-of-plane bend
C  B4      -  Sybyl torsional force constants
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
     $          IT3(NTors),IT4(NTors),IS4(NTors),B4(NTors)
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
C  In this case, we will havr TWO <ffchk> files, one for
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
      READ(JUnit) IType,IAN,IC,IAT,RD,RLEN,B2,IB1,IB2,IB3,RANG,B3,
     $            IT1,IT2,IT3,IT4,IS4,B4,NumB,NumT,NumP
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
      CALL ForceField(NAtoms, IAN,    IType,  NBend,  NTors,
     $                IPRNT,  XC,     IC,     IAT,    RD,
     $                RLEN,   B2,     IB1,    IB2,    IB3,
     $                RANG,   B3,     IT1,    IT2,    IT3,
     $                IT4,    IS4,    B4,     GU,     E,
     $                GC,     HESS,   CutOff, ffcyc,  IStatus)
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
      WRITE(JUnit) IType,IAN,IC,IAT,RD,RLEN,B2,IB1,IB2,IB3,RANG,B3,
     $             IT1,IT2,IT3,IT4,IS4,B4,NumB,NumT,NumP
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
*Deck forcefield
      SUBROUTINE ForceField(NAtoms, IAN,    IType,  NBend,  NTors,
     $                      IPRNT,  XC,     IC,     IAT,    RD,
     $                      RLEN,   B2,     IB1,    IB2,    IB3,
     $                      RANG,   B3,     IT1,    IT2,    IT3,
     $                      IT4,    IS4,    B4,     GU,     E,
     $                      G,      HESS,   CutOff, ffcyc,  IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main Driver for Sybyl v.5.2 molecular mechanics force field
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
C  IAT     -  Sybyl atom type
C  RD      -  distance matrix
C  RLEN    -  Sybyl equilibrium bond lengths
C  B2      -  Sybyl stretching force constants
C  IB1     -  1st atom in bend
C  IB2     -  2nd atom in bend
C  IB3     -  3rd atom in bend
C  RANG    -  Sybyl equilibrium bond angles
C  B3      -  Sybyl bending force constants
C  IT1     -  1st atom in torsion/out-of-plane bend
C  IT2     -  2nd atom in torsion/out-of-plane bend
C  IT3     -  3rd atom in torsion/out-of-plane bend
C  IT4     -  4th atom in torsion/out-of-plane bend
C  IS4     -  Sybyl sign factor for torsion/out-of-plane bend
C  B4      -  Sybyl torsional force constants
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
     $          IT3(NTors),IT4(NTors),IS4(NTors),B4(NTors)
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
C  functions taken from the TRIPOS force field (lone pairs
C  are neglected?). The Gradient is calculated analytically
C  by differentiating the energy with respect to Cartesian
C  displacements and the Hessian is estimated by central
C  differences on the gradient.
C
C  ..............................................................
C     ***  IMPORTANT  ***
C     The TRIPOS force field was originally coded from the SYBYL
C     Theory Manual (see also J.Comp.Chem. 10 (1989) 982).
C     Comparisons with the version provided by TRIPOS themselves
C     in ALCHEMY II indicate that the force field is incorrectly
C     formulated in the manual - there is a missing factor of
C     1/2 in most of the energy terms.
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
        CALL RdPQBFF(40,'Sybyl',NAtoms,Nqm,IAT,IC)
        CALL AddLink('UFF',NAtoms,NLink,XC,IAT,IC)
       Else
        CALL RdPQBFF(40,'Sybyl',NAtoms,NAtoms,IAT,IC)
        ReadPQB = .True.
       EndIf
c
       CLOSE (UNIT=40,STATUS='KEEP')
       GO TO 15
C
 10    CONTINUE
       CALL CONNECTF(NAtoms, IPRNT,  IAN,    XC,     RD,
     $               IT1,    IT2,    IAT,    IC,     IStatus)
       If(IStatus.LT.0) GO TO 99
C
C  Fill the TRIPOS interaction and parameter arrays
C
 15    CONTINUE
       CALL FillData(NAtoms, NBend,  NTors,  IPRNT,  XC,
     $               IAT,    IC,     RD,     RLEN,   B2,
     $               IB1,    IB2,    IB3,    RANG,   B3,
     $               IT1,    IT2,    IT3,    IT4,    IS4, B4)
c
       If(IType.GE.10)
     $    CALL PREAMBLE(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $                  IC,     RD,     RLEN,   B2,     IB1,
     $                  IB2,    IB3,    RANG,   B3,     IT1,
     $                  IT2,    IT3,    IT4,    IS4,    B4,
     $                  G,      GU)
c
      ENDIF
c
      ffcyc = ffcyc + 1
c
      IF(IType.EQ.1.OR.IType.EQ.11) THEN
C
C  form the TRIPOS Hessian by central difference on the gradient
C
       If(IPRNT.GT.1)
     $    WRITE(6,*) '  Hessian Calculated Using Sybyl-5.2 Force Field'
c
       CALL ZeroIT(HESS,NAT3*NAT3)
c
       DO 20 I=1,NAT3
C
C  step up
C
       XC(I) = XC(I) + delta
c
       CALL EGFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $           IC,     RD,     RLEN,   B2,     IB1,
     $           IB2,    IB3,    RANG,   B3,     IT1,
     $           IT2,    IT3,    IT4,    IS4,    B4,
     $           EU,     GU)
C
C  step down
C
       XC(I) = XC(I) - delta - delta
c
       CALL EGFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $           IC,     RD,     RLEN,   B2,     IB1,
     $           IB2,    IB3,    RANG,   B3,     IT1,
     $           IT2,    IT3,    IT4,    IS4,    B4,
     $           ED,     G)
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
       CALL EGFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $           IC,     RD,     RLEN,   B2,     IB1,
     $           IB2,    IB3,    RANG,   B3,     IT1,
     $           IT2,    IT3,    IT4,    IS4,    B4,
     $           E,      G)
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
*Deck eff
      SUBROUTINE EFF(NAtoms, CutOff, IPRNT,  C,      IAT,
     $               IC,     RD,     RLEN,   B2,     IB1,
     $               IB2,    IB3,    RANG,   B3,     IT1,
     $               IT2,    IT3,    IT4,    IS4,    B4,
     $               E)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(3,NAtoms),RD(NAtoms,NAtoms)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      REAL*8 RANG(*),B3(*),B4(*)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),IS4(*)
      DIMENSION IB1(*),IB2(*),IB3(*),IT1(*),IT2(*),IT3(*),IT4(*)
      DIMENSION RVdW(31)
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (FOUR=4.0D0,TWELVE=12.0D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
      PARAMETER (TOLLZERO=1.0D-8)
c
      DATA RVdW/1.7, 1.7, 1.7, 1.7, 1.55,1.55,1.55,1.52,1.52,1.8,
     $          1.55,1.8, 1.5, 1.85,1.75,1.47,1.98,1.8, 1.55,0.0,
     $          1.2, 1.2, 1.2, 1.2, 1.2, 0.0, 1.2, 1.55,1.7, 1.7,1.55/
C
C
C   THIS ROUTINE CALCULATES THE ENERGY ONLY
C   of a molecule using the SYBYL v.5.2 forcefield
C   Energy contributions coded are:
C
C   1. Bond Stretching
C      ES = CS*(R-R0)**2
C   (R0 is an equilibrium bond length and CS a scaling parameter
C   for each bond; if no parameters are known R0 is taken as the
C   initial bond length, R, and Cs = 600)
C
C   2. Non-Bonded VdW interaction (1,3 and H-bond excluded)
C      EV = CV*( 1/(R/RV)**12 - 2/(R/RV)**6)
C   (RV is the sum of the vdw radii of the two atoms and CV is a
C   scaling parameter taken as the square root of the product of
C   the hardness parameters for each atom)
C
C   3. Bending
C      EB = CB*(Th-Th0)**2
C   (Th0 is an equilibrium bond angle and CB a scaling parameter
C   for each bond angle; if no parameters are known Th0 is taken
C   as the original angle and CB = 0.02)
C
C   4. Torsion
C      ET = CT*( 1 + s/|s|*Cos(|s|*Dih)
C   (s and CT are scaling parameters for each torsion; if no
C   parameters are known for the two bonded atoms use s=3 and CT = 0.2)
C
C   5. Out-of-Plane Bend
C      EP = CP*d**2
C   (d is the distance from the atom involved to the plane defined by
C   its three attached atoms and CP is a scaling parameter)
C
C
      ES = ZERO
      EV = ZERO
      EB = ZERO
      ET = ZERO
      EP = ZERO
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
C     (For 1,3-interactions and H-bonding B2(I,J) has been
C      set to zero in routine FillB2)
C
        RV = (RVdW(I1)+RVdW(J1))/R
c
        EV = EV + B2(I,J)*( (RV**12) - TWO*(RV**6) )
c
      ENDIF
 20   CONTINUE
C
C
C  3. Bending
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
      Th = ACOS(CosTh)*DEGREE
c
      EB = EB + B3(M)*(Th-RANG(M))**2
c
 30   CONTINUE
C
C
C  4. Torsions
C
      If(NAtoms.LT.4) GO TO 95
c
      DO 40 M=1,NumT
      I=IT1(M)
      J=IT2(M)
      K=IT3(M)
      L=IT4(M)

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
c
      S = DFLOAT(IS4(M))
      Dih = Abs(S)*Dih
c
      ET = ET + B4(M)*(ONE + COS(Dih)*S/Abs(S))
c
 40   CONTINUE
C
C
C  5. Out-of-Plane bends
C
      DO 50 M=NumT+1,NumT+NumP
      I=IT1(M)
      J=IT2(M)
      K=IT3(M)
      L=IT4(M)
c
      XIJ = C(1,I) - C(1,J)
      XJK = C(1,J) - C(1,K)
      XJL = C(1,J) - C(1,L)
      YIJ = C(2,I) - C(2,J)
      YJK = C(2,J) - C(2,K)
      YJL = C(2,J) - C(2,L)
      ZIJ = C(3,I) - C(3,J)
      ZJK = C(3,J) - C(3,K)
      ZJL = C(3,J) - C(3,L)

      R1(1) = -XJK
      R1(2) = -YJK
      R1(3) = -ZJK
      R2(1) = -XJL
      R2(2) = -YJL
      R2(3) = -ZJL
      R3(1) = -XIJ
      R3(2) = -YIJ
      R3(3) = -ZIJ
      CALL VECMUL(R1,R2,R4)
c
      dn = SProd(3,R4,R3)
      If(dn.LT.TOLLZERO) GO TO 50
c
      dm = SQRT(SProd(3,R4,R4))
      d =  -dn/dm
c
      EP = EP + B4(M)*d*d
c
 50   CONTINUE
c
 95   CONTINUE
C
C  Form the total molecular energy
C
      ES = HALF*ES
      EB = HALF*EB
      ET = HALF*ET
      EP = HALF*EP
c
      E = ES + EV + EB + ET + EP
c
      if(iprnt.gt.2) write(6,1000) E,ES,EB,EP,ET,EV
      RETURN
c
 1000 FORMAT(/,' TOTAL ENERGY       ',f12.8,/,
     $         ' Stretching         ',f12.8,/,
     $         ' Bending            ',f12.8,/,
     $         ' Out-of-Plane Bend  ',f12.8,/,
     $         ' Torsion            ',f12.8,/,
     $         ' VdW energy         ',f12.8)
c
      END
*Deck egff
      SUBROUTINE EGFF(NAtoms, CutOff, IPRNT,  C,      IAT,
     $                IC,     RD,     RLEN,   B2,     IB1,
     $                IB2,    IB3,    RANG,   B3,     IT1,
     $                IT2,    IT3,    IT4,    IS4,    B4,
     $                E,      G)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(3,NAtoms),RD(NAtoms,NAtoms)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      REAL*8 RANG(*),B3(*),B4(*),G(3,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),IS4(*)
      DIMENSION IB1(*),IB2(*),IB3(*),IT1(*),IT2(*),IT3(*),IT4(*)
      DIMENSION RVdW(31)
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (FOUR=4.0D0,TWELVE=12.0D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
      PARAMETER (TOLLZERO=1.0D-8)
c
      DATA RVdW/1.7, 1.7, 1.7, 1.7, 1.55,1.55,1.55,1.52,1.52,1.8,
     $          1.55,1.8, 1.5, 1.85,1.75,1.47,1.98,1.8, 1.55,0.0,
     $          1.2, 1.2, 1.2, 1.2, 1.2, 0.0, 1.2, 1.55,1.7, 1.7,1.55/
C
C
C   THIS ROUTINE CALCULATES THE ENERGY AND CARTESIAN GRADIENT
C   of a molecule using the SYBYL v.5.2 forcefield
C   Energy contributions coded are:
C
C   1. Bond Stretching
C      ES = CS*(R-R0)**2
C   (R0 is an equilibrium bond length and CS a scaling parameter
C   for each bond; if no parameters are known R0 is taken as the
C   initial bond length, R, and Cs = 600)
C
C   2. Non-Bonded VdW interaction (1,3 and H-bond excluded)
C      EV = CV*( 1/(R/RV)**12 - 2/(R/RV)**6)
C   (RV is the sum of the vdw radii of the two atoms and CV is a
C   scaling parameter taken as the square root of the product of
C   the hardness parameters for each atom)
C
C   3. Bending
C      EB = CB*(Th-Th0)**2
C   (Th0 is an equilibrium bond angle and CB a scaling parameter
C   for each bond angle; if no parameters are known Th0 is taken
C   as the original angle and CB = 0.02)
C
C   4. Torsion
C      ET = CT*( 1 + s/|s|*Cos(|s|*Dih)
C   (s and CT are scaling parameters for each torsion; if no
C   parameters are known for the two bonded atoms use s=3 and CT = 0.2)
C
C   5. Out-of-Plane Bend
C      EP = CP*d**2
C   (d is the distance from the atom involved to the plane defined by
C   its three attached atoms and CP is a scaling parameter)
C
C
      ES = ZERO
      EV = ZERO
      EB = ZERO
      ET = ZERO
      EP = ZERO
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
C     (For 1,3-interactions and H-bonding B2(I,J) has been
C      set to zero in routine FillB2)
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
 20   CONTINUE
C
C
C  3. Bending
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
      Th = ACOS(CosTh)*DEGREE
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
      DCB = HALF*B3(M)*(Th-RANG(M))*(-FOUR*DEGREE/SinTh)/(D1*D1)
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
C
C
C  4. Torsions
C
      If(NAtoms.LT.4) GO TO 95
c
      DO 40 M=1,NumT
      I=IT1(M)
      J=IT2(M)
      K=IT3(M)
      L=IT4(M)

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
      S = DFLOAT(IS4(M))
      Dih = Abs(S)*Dih
c
      ET = ET + B4(M)*(ONE + COS(Dih)*S/Abs(S))
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
      n = IAbs(IS4(M))
      IF(n.EQ.1) THEN
       ST = ONE
      ELSE IF(n.EQ.2) THEN
       ST = TWO*CosDih
      ELSE IF(n.EQ.3) THEN
       ST = (TWO*CosDih)**2 - ONE
      ENDIF
      DCT = HALF*B4(M)*S*ST/(A1*A1)
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
C
C
C  5. Out-of-Plane bends
C
      DO 50 M=NumT+1,NumT+NumP
      I=IT1(M)
      J=IT2(M)
      K=IT3(M)
      L=IT4(M)
c
      XIJ = C(1,I) - C(1,J)
      XJK = C(1,J) - C(1,K)
      XJL = C(1,J) - C(1,L)
      YIJ = C(2,I) - C(2,J)
      YJK = C(2,J) - C(2,K)
      YJL = C(2,J) - C(2,L)
      ZIJ = C(3,I) - C(3,J)
      ZJK = C(3,J) - C(3,K)
      ZJL = C(3,J) - C(3,L)

      R1(1) = -XJK
      R1(2) = -YJK
      R1(3) = -ZJK
      R2(1) = -XJL
      R2(2) = -YJL
      R2(3) = -ZJL
      R3(1) = -XIJ
      R3(2) = -YIJ
      R3(3) = -ZIJ
      CALL VECMUL(R1,R2,R4)
c
      dn = SProd(3,R4,R3)
      If(dn.LT.TOLLZERO) GO TO 50
c
      dm = SQRT(SProd(3,R4,R4))
      d =  -dn/dm
c
      EP = EP + B4(M)*d*d
C
C  Set terms for out-of-plane bend derivatives
C
C  First set the derivatives of each term
C  1. The numerator
C      d/dx ( SProd(3,R4,R3) )
C
      DT1XI = YJL*ZJK - YJK*ZJL
      DT1XJ = -DT1XI + ZIJ*(YJK-YJL) - YIJ*(ZJK-ZJL)
      DT1XK = YJL*ZIJ - YIJ*ZJL
      DT1XL = YIJ*ZJK - YJK*ZIJ
c
      DT1YI = XJK*ZJL - XJL*ZJK
      DT1YJ = -DT1YI + XIJ*(ZJK-ZJL) - ZIJ*(XJK-XJL)
      DT1YK = XIJ*ZJL - XJL*ZIJ
      DT1YL = XJK*ZIJ - XIJ*ZJK
c
      DT1ZI = XJL*YJK - XJK*YJL
      DT1ZJ = -DT1ZI + YIJ*(XJK-XJL) - XIJ*(YJK-YJL)
      DT1ZK = XJL*YIJ - XIJ*YJL
      DT1ZL = XIJ*YJK - XJK*YIJ
C
C  2. The Denominator
C      d/dx ( SQRT(SProd(3,R4,R4)) )
C
      RA = (XJK*YJL-XJL*YJK)
      RB = (XJK*ZJL-XJL*ZJK)
      RC = (YJK*ZJL-YJL*ZJK)
c
      DT2XK = -YJL*RA - ZJL*RB
      DT2XL =  YJK*RA + ZJK*RB
      DT2XJ = -(YJK-YJL)*RA - (ZJK-ZJL)*RB
c
      DT2YK =  XJL*RA - ZJL*RC
      DT2YL = -XJK*RA + ZJK*RC
      DT2YJ = (XJK-XJL)*RA - (ZJK-ZJL)*RC
c
      DT2ZK =  XJL*RB + YJL*RC
      DT2ZL = -XJK*RB - YJK*RC
      DT2ZJ = (XJK-XJL)*RB + (YJK-YJL)*RC
C
C  Now construct the final derivative
C
      DCP = -B4(M)*d/(dm*dm)
c
      G(1,I) = G(1,I) + DCP*dm*DT1XI
      G(1,J) = G(1,J) + DCP*(dm*DT1XJ + d*DT2XJ)
      G(1,K) = G(1,K) + DCP*(dm*DT1XK + d*DT2XK)
      G(1,L) = G(1,L) + DCP*(dm*DT1XL + d*DT2XL)
      G(2,I) = G(2,I) + DCP*dm*DT1YI
      G(2,J) = G(2,J) + DCP*(dm*DT1YJ + d*DT2YJ)
      G(2,K) = G(2,K) + DCP*(dm*DT1YK + d*DT2YK)
      G(2,L) = G(2,L) + DCP*(dm*DT1YL + d*DT2YL)
      G(3,I) = G(3,I) + DCP*dm*DT1ZI
      G(3,J) = G(3,J) + DCP*(dm*DT1ZJ + d*DT2ZJ)
      G(3,K) = G(3,K) + DCP*(dm*DT1ZK + d*DT2ZK)
      G(3,L) = G(3,L) + DCP*(dm*DT1ZL + d*DT2ZL)
c
 50   CONTINUE
c
 95   CONTINUE
C
C  Form the total molecular energy
C
      ES = HALF*ES
      EB = HALF*EB
      ET = HALF*ET
      EP = HALF*EP
c
      E = ES + EV + EB + ET + EP
c
      if(iprnt.gt.2) write(6,1000) E,ES,EB,EP,ET,EV
      RETURN
c
 1000 FORMAT(/,' TOTAL ENERGY       ',f12.8,/,
     $         ' Stretching         ',f12.8,/,
     $         ' Bending            ',f12.8,/,
     $         ' Out-of-Plane Bend  ',f12.8,/,
     $         ' Torsion            ',f12.8,/,
     $         ' VdW energy         ',f12.8)
c
      END
*Deck preamble
      SUBROUTINE PREAMBLE(NAtoms, CutOff, IPRNT,  C,      IAT,
     $                    IC,     RD,     RLEN,   B2,     IB1,
     $                    IB2,    IB3,    RANG,   B3,     IT1,
     $                    IT2,    IT3,    IT4,    IS4,    B4,
     $                    G,      CT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 C(3*NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms)
      REAL*8 B2(NAtoms,NAtoms),RANG(*),B3(*),B4(*),G(3*NAtoms),
     $       CT(3*NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms),IS4(*)
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
      CALL EGFF(NAtoms, CutOff, IPRNT,  C,      IAT,
     $          IC,     RD,     RLEN,   B2,     IB1,
     $          IB2,    IB3,    RANG,   B3,     IT1,
     $          IT2,    IT3,    IT4,    IS4,    B4,
     $          E0,     G)
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
      CALL EFF(NAtoms, CutOff, IPRNT,  C,      IAT,
     $         IC,     RD,     RLEN,   B2,     IB1,
     $         IB2,    IB3,    RANG,   B3,     IT1,
     $         IT2,    IT3,    IT4,    IS4,    B4,
     $         E1)
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
*Deck connectf
      SUBROUTINE CONNECTF(NAtoms, IPRNT,  IAN,    C,      RD,
     $                    NBOND,  IVAL,   IAT,    IC,     IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL ERROR
C
      REAL*8 C(3,NAtoms),RD(NAtoms,NAtoms)
      DIMENSION IAN(NAtoms),IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION NBOND(NAtoms),IVAL(NAtoms)
      DIMENSION ISYBYL(54)
      Character*256 jobname
c
      Common /job/jobname,lenJ
      Parameter (IUnit=1)
C
C  ISYBYL - Likely SYBYL atom type for each atom
C
C                    H He Li Be  B  C  N  O F  Ne
      DATA ISYBYL/  13, 0,24, 0, 0, 1, 5, 8,16, 0,
C
C                         Na Mg Al Si  P  S Cl Ar
     *                    21, 0,25,27,12,10,15, 0,
C
C                          K Ca
     *                    22,23,
C
C                   Sc Ti  V Cr Mn Fe Co Ni Cu Zn
     *               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
C
C                               Ga Gs As Se Br Kr
     *                           0, 0, 0, 0,14, 0,
C
C                         Rb Sr
     *                     0, 0,
C
C                    Y Zr Nb Mo To Ru Rh Pd Ag Cd
     *               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
C
C                                In Sn Sb Te I Xe
     *                           0, 0, 0, 0,17, 0 /
C
C
C  Assign/read in the atom types and connectivity for the SYBYL force field
C  If there is any user-defined bonding data to be read in, it should be
C  in the file <jobname.fld>
C
C  The standard SYBYL atom types and their numerical values
C  are listed below (taken from the SYBYL Theory Manual
C  version 5.2 February 1989)
C
C    1     C.3     carbon sp3
C    2     C.2     carbon sp2
C    3     C.ar    carbon aromatic
C    4     C.1     carbon sp
C    5     N.3     nitrogen sp3
C    6     N.2     nitrogen sp2
C    7     N.1     nitrogen sp
C    8     0.3     oxygen sp3
C    9     0.2     oxygen sp2
C    10    S.3     sulphur sp3
C    11    N.ar    nitrogen aromatic
C    12    P.3     phosphorus sp3
C    13    H       hydrogen
C    14    Br      bromine
C    15    Cl      chlorine
C    16    F       fluorine
C    17    I       iodine
C    18    S.2     sulphur sp2
C    19    N.pl3   nitrogen trigonal planar
C    20    LP      lone pair                     REDUNDANT
C    21    Na      sodium
C    22    K       potassium
C    23    Ca      calcium
C    24    Li      lithium
C    25    Al      aluminium
C    26    Du      dummy                         REDUNDANT
C    27    Si      silicon
C    28    N.am    nitrogen amide
C    29    S.0     sulphoxide sulphur
C    30    S.02    sulphone sulphur
C    31    N.4     nitrogen sp3 +ve charge
C    32            unknown atom type
C
C  Standard bond types are:
C  1 = single bond; 2 = double bond; 3 = triple bond;
C  4 = amide bond; 5 = aromatic bond
C
C
C  initialize
C
      DO 5 J=1,NAtoms
      DO 5 I=1,NAtoms
      IC(I,J)=0
 5    CONTINUE
C
C  Assign preliminary bonding using model builder
C
      CALL BOND(NAtoms, IPRNT,  IAN,    C,      RD,
     $          NBOND,  IVAL,   IAT,    IC,     IStatus)
C
C  Is there any user-defined bonding data?
C  This will overwrite bonding produced by <BOND>
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.fld',
     $      FORM='FORMATTED',STATUS='OLD',ERR=96)
C
C  data found
C  first line in file should just be a title
C
      If(IPRNT.GT.1) WRITE(6,1000)
c
      READ(IUnit,*)
C
C  read in lines of connectivity data
C     atom I    atom J      Sybyl bond type
C
 10   CONTINUE
      READ (IUnit,*,END=95) I,J,IB
c
      IC(I,J)=IB
      IC(J,I)=IB
      GO TO 10
c
 95   CONTINUE
      CLOSE (UNIT=IUnit,STATUS='DELETE')
c
 96   CONTINUE
C
C
C  Connectivity matrix complete
C  determine SYBYL atom types
C
      DO 20 IA=1,NAtoms
C
C  check atoms that can have multiple SYBYL types
C
      IF(IAN(IA).EQ.6) THEN
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.EQ.1) IAT(IA) = 1
       If(LIC.EQ.2) IAT(IA) = 2
       If(LIC.EQ.3) IAT(IA) = 4
       If(LIC.EQ.5) IAT(IA) = 3
      ELSE IF(IAN(IA).EQ.7) THEN
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.EQ.1) IAT(IA) = 5
       If(LIC.EQ.2) IAT(IA) = 6
       If(LIC.EQ.3) IAT(IA) = 7
       If(LIC.EQ.5) IAT(IA) = 11
      ELSE IF(IAN(IA).EQ.8) THEN
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.EQ.1) IAT(IA) = 8
       If(LIC.EQ.2) IAT(IA) = 9
      ELSE IF(IAN(IA).EQ.16) THEN
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.EQ.1) IAT(IA) = 10
       If(LIC.EQ.2) IAT(IA) = 18
      ELSE
       IAT(IA) = ISYBYL(IAN(IA))
      ENDIF
 20   CONTINUE
C
C Print Out atom types and connectivities selected
C and check for errors
C
      If(IPRNT.GT.2) WRITE(6,1100)
c
      ERROR=.FALSE.
      DO 30 IA=1,NAtoms
      If(IPRNT.GT.2) WRITE(6,1200) IA,IAT(IA)
      If(IAT(IA).LT.0) ERROR=.TRUE.
 30   CONTINUE
c
      If(IPRNT.GT.2) WRITE(6,1300)
      DO 40 IA=2,NAtoms
      DO 40 JA=1,IA-1
      IF(IC(IA,JA).NE.0) THEN
       If(IPRNT.GT.2) WRITE(6,1400) IA,JA,IC(IA,JA)
       If(IC(IA,JA).LT.0) ERROR=.TRUE.
      ENDIF
 40   CONTINUE
c
      IF(ERROR) THEN
       WRITE(6,1500)
       ISTATUS=-2
       RETURN
      ENDIF
C
      RETURN
C
 1000 FORMAT(' Sybyl Bonding Data Read In')
 1100 FORMAT(/,2X,' Sybyl Atom types',//,
     $         3X,' ATOM  TYPE',/,
     $         3X,' ----  ----')
 1200 FORMAT (4X,I3,2X,I4)
 1300 FORMAT (/,2X,'BONDS info',//,
     $          3X,' Atom1  Atom2  Bond Type',/,
     $          3X,' -----  -----  ---------')
 1400 FORMAT(4X,I4,4X,I4,4X,I4)
 1500 FORMAT(/,2X,'***ERROR*** Unable to determine connectivities.',
     $            ' Data must be read in')
c
      END
*Deck bond
      SUBROUTINE BOND(NAtoms, IPRNT,  IAN,    CATM,   RD,
     $                NBOND,  IVAL,   IAT,    IC,     IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL ERROR,DONE
C
      DIMENSION IAN(NAtoms),CATM(3,NAtoms),RD(NAtoms,NAtoms)
      DIMENSION NBOND(NAtoms),IVAL(NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      REAL*8 NND1,NND2
C
      DIMENSION RADII(54),NVAL(54)
C
C                    H       He
      DATA RADII/  0.6D0,  0.5D0,
C
C                    Li      Be      B       C       N       O
     *             1.25D0, 1.00D0, 0.90D0, 0.85D0, 0.75D0, 0.75D0,
C
C                                                    F       Ne
     *                                             0.75D0, 0.75D0,
C
C                    Na      Mg      Al      Si      P       S
     *             1.60D0, 1.40D0, 1.30D0, 1.20D0, 1.15D0, 1.10D0,
C
C                                                    Cl      Ar
     *                                             1.05D0, 1.00D0,
C
C                    K       Ca
     *             2.00D0, 1.75D0,
C
C                    Sc      Ti      V       Cr      Mn      Fe
     *             1.60D0, 1.50D0, 1.40D0, 1.40D0, 1.40D0, 1.40D0,
C
C                    Co      Ni      Cu      Zn
     *             1.40D0, 1.40D0, 1.38D0, 1.35D0,
C
C                    Ga      Ge      As      Se      Br      Kr
     *             1.33D0, 1.30D0, 1.28D0, 1.26D0, 1.24D0, 1.22D0,
C
C                    Rb      Sr
     *             2.20D0, 1.95D0,
C
C                    Y       Zr      Nb      Mo      Tc      Ru
     *             1.75D0, 1.60D0, 1.45D0, 1.45D0, 1.45D0, 1.45D0,
C
C                    Rh      Pd      Ag      Cd
     *             1.45D0, 1.45D0, 1.50D0, 1.50D0,
C
C                    In      Sn      Sb      Te      I       Xe
     *             1.55D0, 1.45D0, 1.42D0, 1.40D0, 1.38D0, 1.35D0 /
C
C  NVAL - maximum valency of each atom (normal)
C
C                 H  He Li Be B  C  N  O  F  Ne
      DATA NVAL/  1, 0, 1, 2, 3, 4, 3, 2, 1, 0,
C
C                       Na Mg Al Si P  S  Cl Ar
     *                  1, 2, 3, 4, 3, 2, 1, 0,
C
C                       K  Ca
     *                  1, 2,
C
C                 Sc Ti V  Cr Mn Fe Co Ni Cu Zn
     *            3, 4, 5, 6, 5, 4, 3, 2, 1, 0,
C
C                             Ga Ge As Se Br Kr
     *                        3, 4, 3, 2, 1, 0,
C
C                       Rb Sr
     *                  1, 2,
C
C                 Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd
     *            3, 4, 5, 6, 5, 4, 3, 2, 1, 0,
C
C                             In Sn Sb Te I  Xe
     *                        3, 4, 3, 2, 1, 0 /
C
      DATA CCD1/1.30d0/, CCD2/1.38d0/, CCT1/1.17d0/, CCT2/1.22d0/,
     *     CCA1/1.38d0/, CCA2/1.42d0/, CND1/1.19d0/, CND2/1.33d0/,
     *     CNA1/1.33d0/, CNA2/1.38d0/, NND1/1.22d0/, NND2/1.28d0/
C
C
C  Find which atoms are likely to be bonded.
C
      DO 5 IA=1,NAtoms
      NBOND(IA) = 0
      IVAL(IA) = NVAL(IAN(IA))
 5    CONTINUE
c
      DO 10 IA=2,NAtoms
      DO 10 JA=1,IA-1
      R = SQRT( (CATM(1,IA) - CATM(1,JA))**2 +
     $          (CATM(2,IA) - CATM(2,JA))**2 +
     $          (CATM(3,IA) - CATM(3,JA))**2 )
      RD(IA,JA)  = R
      RD(JA,IA)  = R
      BL = RADII(IAN(IA)) + RADII(IAN(JA))
      IF(R.LE.BL) THEN
       IC(IA,JA) = -1
       IC(JA,IA) = -1
       NBOND(IA) = NBOND(IA) + 1
       NBOND(JA) = NBOND(JA) + 1
      ENDIF
 10   CONTINUE
C
C  Check that atomic valency has not been exceeded or is zero
C
      ERROR=.FALSE.
      DO 15 IA=1,NAtoms
      IF(NBOND(IA).EQ.0) THEN
       WRITE(6,1000) IA
       ERROR=.TRUE.
      ELSE IF(NBOND(IA).GT.IVAL(IA)) THEN
       WRITE(6,1100) IA
       ERROR=.TRUE.
      ENDIF
 15   CONTINUE
      IF(ERROR) THEN
       WRITE(6,1200)
       ISTATUS=-1
       RETURN
      ENDIF
C
C  Go through and find definite single links only
C  i.e., atoms which are only bonded to one other atom
C
 20   CONTINUE
      DONE=.TRUE.
c
      DO 25 IA=2,NAtoms
      DO 25 JA=1,IA-1
      IF(IC(IA,JA).LT.0) THEN
       IF(NBOND(IA).EQ.1) THEN
        DONE=.FALSE.
        IC(IA,JA) = IVAL(IA)
        IC(JA,IA) = IVAL(IA)
        IVAL(JA) = IVAL(JA) - IVAL(IA)
        IVAL(IA) = 0
        NBOND(IA) = 0
        NBOND(JA) = NBOND(JA) - 1
       ELSE IF(NBOND(JA).EQ.1) THEN
        DONE=.FALSE.
        IC(IA,JA) = IVAL(JA)
        IC(JA,IA) = IVAL(JA)
        IVAL(IA) = IVAL(IA) - IVAL(JA)
        IVAL(JA) = 0
        NBOND(JA) = 0
        NBOND(IA) = NBOND(IA) - 1
       ENDIF
      ENDIF
 25   CONTINUE
C
C  The process of finding single links can be iterated since
C  when definite connectivities are known, this reduces the
C  options on remaining bonds which may themselves be found
C  on another pass
C
      If(.NOT.DONE) GO TO 20
C
C  Now check for cyclic systems which may not have been
C  detected by the above code
C
      DO 30 IA=2,NAtoms
      IF(IVAL(IA).EQ.NBOND(IA).AND.NBOND(IA).NE.0) THEN
       DO 35 JA=1,IA-1
       IF(IC(IA,JA).LT.0) THEN
        IC(IA,JA) = 1
        IC(JA,IA) = 1
        IVAL(IA) = IVAL(IA) - 1
        IVAL(JA) = IVAL(JA) - 1
        NBOND(IA) = NBOND(IA) - 1
        NBOND(JA) = NBOND(JA) - 1
       ENDIF
 35    CONTINUE
      ENDIF
30    CONTINUE
C
C  Check for completion
C
      DONE=.TRUE.
      DO 40 IA=1,NAtoms
      If(NBOND(IA).NE.0) DONE=.FALSE.
 40   CONTINUE

      If(DONE) RETURN
C
C
C  At this point there are one or more indeterminate
C  multiple bonds
C  ............................................
C     TEMPORARY FOR BONDS INVOLVING C,N ONLY
C  ............................................
C
C  Attempt to arbitrate multiple bonds by comparing atomic
C  distances with standard lengths
C
      DO 50 IA=2,NAtoms
      DO 50 JA=1,IA-1
      IF(IC(IA,JA).LT.0) THEN
       ITYPE=0
       If(IAN(IA).EQ.6.AND.IAN(JA).EQ.6) ITYPE=1
       If((IAN(IA).EQ.6.AND.IAN(JA).EQ.7).OR.
     $    (IAN(IA).EQ.7.AND.IAN(JA).EQ.6)) ITYPE=2
       If(IAN(IA).EQ.7.AND.IAN(JA).EQ.7) ITYPE=3
       IF(ITYPE.GT.0) THEN
        R = RD(IA,JA)
        IF(ITYPE.EQ.1) THEN
         IF(R.GE.CCD1.AND.R.LT.CCD2) THEN
          IC(IA,JA) = 2
          IC(JA,IA) = 2
         ELSE IF(R.GE.CCT1.AND.R.LT.CCT2) THEN
          IC(IA,JA) = 3
          IC(JA,IA) = 3
         ELSE IF(R.GE.CCA1.AND.R.LT.CCA2) THEN
          IC(IA,JA) = 5
          IC(JA,IA) = 5
         ELSE
          IC(IA,JA) = 1
          IC(JA,IA) = 1
         ENDIF
        ELSE IF(ITYPE.EQ.2) THEN
         IF(R.GE.CND1.AND.R.LT.CND2) THEN
          IC(IA,JA) = 2
          IC(JA,IA) = 2
         ELSE IF(R.GE.CNA1.AND.R.LT.CNA2) THEN
          IC(IA,JA) = 5
          IC(JA,IA) = 5
         ELSE
          IC(IA,JA) = 1
          IC(JA,IA) = 1
         ENDIF
        ELSE IF(ITYPE.EQ.3) THEN
         IF(R.GE.NND1.AND.R.LT.NND2) THEN
          IC(IA,JA) = 2
          IC(JA,IA) = 2
         ELSE
          IC(IA,JA) = 1
          IC(JA,IA) = 1
         ENDIF
        ENDIF
       ENDIF
      ENDIF
 50   CONTINUE
C
C  .......................
C  Final aromaticity check
C  .......................
C
      CALL ChkAromatic(NAtoms,IAN,IC)
C  ................................................
C
      RETURN
c
 1000 FORMAT(' Atom ',I2,' does not appear to be bonded to anything ')
 1100 FORMAT(' Valency of atom ',I2,' exceeds standard valency')
 1200 FORMAT(/,2X,'***ERROR*** Unable to determine connectivities.',
     $            ' Data must be read in')
c
      END
*Deck filldata
      SUBROUTINE FillData(NAtoms, NBend,  NTors,  IPRNT,  C,
     $                    IAT,    IC,     RD,     RLEN,   B2,
     $                    IB1,    IB2,    IB3,    RANG,   B3,
     $                    IT1,    IT2,    IT3,    IT4,    IS4, B4)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(3,NAtoms),IAT(NAtoms),RD(NAtoms,NAtoms)
      dimension rlen(natoms,natoms),b2(natoms,natoms)
      dimension ib1(*),ib2(*),ib3(*),rang(*),b3(*)
      dimension it1(*),it2(*),it3(*),it4(*),is4(*),b4(*)
      common/numint/numb,numt,nump
C
C  First set all unknown Sybyl types to 32
C
      DO 5 I=1,NAtoms
      If(IAT(I).LT.0) IAT(I)=32
5     CONTINUE
C
C Form the distance matrix
C
      DO 10 J=1,NAtoms
      DO 10 I=1,NAtoms
      RD(I,J) = SQRT( (C(1,I)-C(1,J))**2 + (C(2,I)-C(2,J))**2
     $              + (C(3,I)-C(3,J))**2 )
10    CONTINUE
c
      CALL FillB2(NAtoms,IAT,IC,RD,RLEN,B2)
      CALL FillB3(NAtoms,NBend,IAT,IC,RD,IB1,IB2,IB3,RANG,B3)
      CALL FillB4(NAtoms,NTors,IAT,IC,IT1,IT2,IT3,IT4,IS4,B4)
c
      If(IPRNT.gt.3) THEN
       zero=0.0d0
       write(6,*) ' TRIPOS Parameter Arrays are:'
       write(6,*) ' Distance Parameters'
       do 11 i=2,natoms
       do 11 j=1,i-1
       if(b2(i,j).ne.zero) write(6,1111) i,j,rlen(i,j),b2(i,j)
 1111  format(1x,2i4,2f12.6)
 11    continue
       if(numb.gt.0) then
        write(6,*) ' Angle Parameters'
        do 22 i=1,numb
        write(6,2222) ib1(i),ib2(i),ib3(i),rang(i),b3(i)
 2222   format(1x,3i4,2f12.6)
 22     continue
       endif
       if(numt.gt.0) then
        write(6,*) ' Torsion Parameters'
        do 33 i=1,numt
        write(6,3333) it1(i),it2(i),it3(i),it4(i),b4(i),is4(i)
 3333   format(1x,4i4,f12.6,i6)
 33     continue
       endif
       if(nump.gt.0) then
        write(6,*) ' Out-of-Plane Bend Parameters'
        do 44 i=numt+1,numt+nump
        write(6,4444) it1(i),it2(i),it3(i),it4(i),b4(i)
 4444   format(1x,i4,2x,3i4,f12.6)
 44     continue
       endif
      EndIf
      RETURN
      END
*Deck fillb2
      SUBROUTINE FillB2(NAtoms,IAT,IC,RD,RLEN,B2)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IATI(112),IATJ(112),IBT(112),R0(112),CIJ(112),EH(32)
      LOGICAL Donor,Accept
      PARAMETER (NumDAT=112,B2D=600.0d0)
C
      DATA EH/0.107d0,0.107d0,0.107d0,0.107d0,0.095d0,0.095d0,0.095d0,
     $        0.116d0,0.116d0,0.314d0,0.095d0,0.314d0,0.042d0,0.434d0,
     $        0.314d0,0.109d0,0.623d0,0.314d0,0.095d0,0.000d0,0.400d0,
     $        0.400d0,0.600d0,0.400d0,0.042d0,0.000d0,0.042d0,0.095d0,
     $        0.314d0,0.314d0,0.095d0,1.000d0/
c
      DATA IATI/ 4, 4,14, 4, 4, 2, 2, 4, 2, 1,14, 4, 2, 1, 3, 3, 2, 1,
     $           3, 2, 1, 3, 0, 4, 2, 1, 3, 3, 4, 4, 4, 2, 2, 1, 3, 6,
     $           6, 2, 1, 3,13, 2, 1, 3, 2, 2, 1, 3,13, 6,28, 3,11, 2,
     $           1, 3,13, 6, 2,28,19, 2, 1, 3,13, 6,19, 8, 1, 9, 9, 8,
     $           2, 1, 3, 2, 1, 3, 5,31, 9,10, 2, 1, 9, 8, 1, 9, 8,27,
     $          27,27,27,27,27, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 8, 8,
     $          10,10,19,19/
c
      DATA IATJ/ 4, 4, 2, 2, 2, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3,15,15,
     $          15,16,16,16,13,13,13,13,13,17, 7, 6, 6, 6, 6, 6, 6, 6,
     $           6, 5, 5, 5, 5,31,31,31,28,28,28,28,28,28,28,11,11,19,
     $          19,19,19,19, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8,12,12,12,12,
     $          18,18,18,10,10,10,10,10,10,10,29,29,29,29,30,30,30, 1,
     $           2, 3, 4, 8,27, 5, 8,10,19,28, 5, 6, 8,19,28,10,10,28,
     $          19,28,19,28/
c
      DATA  IBT/ 3, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1,
     $           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 2, 2, 1, 1, 1, 2,
     $           1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 5, 5, 1,
     $           1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
     $           2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1,
     $           22*1/
c
      DATA R0/1.204d0,1.380d0,1.890d0,1.440d0,1.440d0,1.335d0,1.470d0,
     $        1.458d0,1.501d0,1.540d0,1.850d0,1.440d0,1.510d0,1.525d0,
     $        1.395d0,1.480d0,1.750d0,1.767d0,1.750d0,1.330d0,1.360d0,
     $        1.330d0,1.008d0,1.056d0,1.089d0,1.100d0,1.084d0,2.050d0,
     $        1.158d0,1.330d0,1.330d0,1.270d0,1.444d0,1.440d0,1.346d0,
     $        1.346d0,1.418d0,1.330d0,1.470d0,1.410d0,1.080d0,1.330d0,
     $        1.470d0,1.410d0,1.345d0,1.345d0,1.450d0,1.416d0,1.000d0,
     $        1.440d0,1.450d0,1.346d0,1.330d0,1.300d0,1.450d0,1.350d0,
     $        1.030d0,1.350d0,1.220d0,1.240d0,1.210d0,1.330d0,1.430d0,
     $        1.390d0,0.950d0,1.405d0,1.400d0,1.480d0,1.830d0,1.490d0,
     $        1.490d0,1.600d0,1.710d0,1.800d0,1.740d0,1.780d0,1.817d0,
     $        1.770d0,1.625d0,1.625d0,1.450d0,2.030d0,1.710d0,1.800d0,
     $        1.450d0,1.500d0,1.800d0,1.450d0,1.500d0,1.840d0,1.840d0,
     $        1.840d0,1.840d0,1.620d0,2.290d0,1.480d0,1.390d0,1.820d0,
     $        1.480d0,1.480d0,1.460d0,1.470d0,1.407d0,1.470d0,1.470d0,
     $        1.670d0,1.540d0,1.407d0,1.670d0,1.670d0,1.480d0,1.480d0/
c
      DATA CIJ/1400.0d0,700.0d0,500.0d0,1340.0d0,1340.0d0,1340.0d0,
     $          700.0d0,640.0d0,639.0d0,633.6d0,500.0d0,1340.0d0,
     $          1340.0d0,640.0d0,1400.0d0,1000.0d0,520.0d0,600.0d0,
     $          513.36d0,1200.0d0,600.0d0,500.0d0,700.0d0,700.0d0,
     $          692.0d0,662.4d0,692.0d0,490.0d0,1600.0d0,1300.0d0,
     $          1300.0d0,1305.94d0,1300.0d0,760.2d0,1305.94d0,1305.94d0,
     $          1300.0d0,1300.0d0,760.0d0,720.0d0,692.0d0,1300.0d0,
     $          760.0d0,720.0d0,870.1d0,870.1d0,677.6d0,1090.08d0,
     $          700.0d0,667.6d0,744.48d0,1305.94d0,1400.0d0,1200.0d0,
     $          676.0d0,1306.0d0,692.0d0,1305.94d0,1555.2d0,1120.0d0,
     $          680.0d0,699.84d0,618.9d0,700.0d0,1007.5d0,1200.0d0,
     $          620.0d0,1172.16d0,407.6d0,1400.0d0,1400.0d0,800.0d0,
     $          400.0d0,381.6d0,700.0d0,360.0d0,381.6d0,360.0d0,
     $          360.0d0,360.0d0,600.0d0,600.0d0,360.0d0,381.6d0,
     $          600.0d0,600.0d0,381.6d0,600.0d0,600.0d0,23*600.0d0/
C
C  ** WARNING! **  All parameters for Silicon (SYBYL atom type 27)
C                  and atoms following made up
C
C
C  This Subroutine fills the RLEN and B2 arrays for bond stretches
C  by going through each data column looking for a match
C  It also sets B2 for any VdW interactions
C
C
      CALL ZeroIT(RLEN,NAtoms*NAtoms)
      CALL ZeroIT(B2,NAtoms*NAtoms)
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
        If(J1.EQ.IATJ(K).AND.IC(I,J).EQ.IBT(K)) GO TO 95
       ENDIF
c
       IF(J1.EQ.IATI(K)) THEN
C
C  match found for J1
C  check it everything else fits
C
        If(I1.EQ.IATJ(K).AND.IC(I,J).EQ.IBT(K)) GO TO 95
       ENDIF
 11    CONTINUE
C
C  no match found
C  check for standard bond to hydrogen, else use default
C
       IF(IATI(K).EQ.13.OR.IATJ(K).EQ.13) THEN
        K=23
        GO TO 95
       ENDIF
c
       RLEN(I,J) = RD(I,J)
       B2(I,J) = B2D
       GO TO 10

95     CONTINUE
       RLEN(I,J) = R0(K)
       B2(I,J) = CIJ(K)
c
      ELSE
C
C  Van der Waals interaction
C  check for 1,3 and H-bonds
C
       DO 20 K=1,NAtoms
       K1=IAT(K)
       Donor = (K1.EQ.5.OR.K1.EQ.6.OR.K1.EQ.8.OR.K1.EQ.16.OR.
     $          K1.EQ.19.OR.K1.EQ.28.OR.K1.EQ.31)
       IF(IC(I,K).NE.0) THEN
        If(IC(J,K).NE.0) GO TO 10
        Accept = (J1.EQ.5.OR.J1.EQ.6.OR.J1.EQ.7.OR.J1.EQ.8.OR.
     $            J1.EQ.9.OR.J1.EQ.11.OR.J1.EQ.16)
        If(I1.EQ.13.AND.Donor.AND.Accept) GO TO 10
       ELSE IF(IC(J,K).NE.0) THEN
        Accept = (I1.EQ.5.OR.I1.EQ.6.OR.I1.EQ.7.OR.I1.EQ.8.OR.
     $            I1.EQ.9.OR.I1.EQ.11.OR.I1.EQ.16)
        If(J1.EQ.13.AND.Donor.AND.Accept) GO TO 10
       ENDIF
 20    CONTINUE
C
C set VdW interaction term
C
       B2(I,J) = SQRT(EH(I1)*EH(J1))
      ENDIF
c
10    CONTINUE
c
      RETURN
      END
*Deck fillb3
      SUBROUTINE FillB3(NAtoms,NBend,IAT,IC,RD,IB1,IB2,IB3,RANG,B3)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RD(NAtoms,NAtoms),RANG(*),B3(*)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IB1(*),IB2(*),IB3(*)
      COMMON/NUMINT/NumB,NumT,NumP
C
C  This Subroutine calculates the number of bends. storing the
C  atoms involved in the arrays IB1, IB2 and IB3, and fills
C  the RANG and B3 arrays by going through the relevant
C  TRIPOS data columns looking for a match
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
     $   CALL AMATCH(NAtoms,I,J,K,II,JJ,KK,RD,IB1,IB2,IB3,RANG,B3)
      If(IC(I,K).NE.0.AND.IC(K,J).NE.0)
     $   CALL AMATCH(NAtoms,I,K,J,II,KK,JJ,RD,IB1,IB2,IB3,RANG,B3)
      If(IC(J,I).NE.0.AND.IC(I,K).NE.0)
     $   CALL AMATCH(NAtoms,J,I,K,JJ,II,KK,RD,IB1,IB2,IB3,RANG,B3)
c
c -- check number of bends does not exceed maximum allowed
      If(NumB.GT.NBend) Call nerror(2,'Force Field module',
     $           'More Bends in System than Dimensioned for',0,0)
c
 10   CONTINUE
c
      RETURN
      END
*Deck amatch
      SUBROUTINE AMATCH(NAtoms,I,J,K,I1,J1,K1,RD,IB1,IB2,IB3,RANG,B3)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RD(NAtoms,NAtoms),RANG(*),B3(*)
      DIMENSION IB1(*),IB2(*),IB3(*)
      DIMENSION IATI(202),IATJ(202),IATK(202),Theta(202),CIJK(202)
      DIMENSION ISTART(32),INUM(32)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ONE=1.0D0,TWO=2.0D0,B3D=0.02d0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
c
      DATA IATI/ 0, 4, 2, 1, 3, 7, 0,14,14, 2, 4, 2, 1, 4, 2,
     $           1, 3, 2, 3,15, 0, 4, 2, 1, 3, 2, 1, 6, 5, 2,
     $           1, 3, 6,28, 2, 1, 6, 4, 2, 1, 3, 5,28,19, 2,
     $           1, 3,28, 9, 6,28, 9, 0, 2, 4, 2, 1, 2, 1, 3,
     $           1,15, 3,16, 0, 2,13, 1, 2, 1, 3, 2, 1, 3,13,
     $           6,28,16, 2, 1, 3, 2, 1, 3,28, 8, 1, 2, 1,28,
     $           2, 1, 3,28, 8, 0,14, 2, 1, 3, 3, 3, 3, 1, 3,
     $           6, 3, 3, 6,28, 2, 1, 3,28, 1, 3, 6,28,19, 1,
     $           3, 3, 3, 0, 0, 4, 2, 2, 2, 1, 3, 2, 1, 3, 2,
     $           2, 2, 0, 2, 1, 1, 3, 1, 0, 1, 0, 2, 2, 1, 2,
     $           1, 3, 2, 1, 2, 1, 3, 2, 1, 3, 2, 1, 3, 9, 0,
     $           3, 0, 2, 2, 1, 3, 2, 3, 2, 9, 0, 2, 2, 1, 2,
     $           1, 3, 2, 1, 1, 0, 9,26, 9, 8, 0, 1, 3, 0, 2,
     $           1, 3, 1, 0, 0, 9, 0/
c
      DATA IATJ/6*4,46*2,43*1,28*3,7,13*6,6*5,2*31,19*28,2*11,
     $          9*19,10*8,5*12,3*18,5*10,29,30,30,27/
c
      DATA IATK/ 0, 2, 7, 7, 7, 8, 0,14, 2, 2, 1, 1, 1, 3, 3,
     $           3, 3,15,15,15,13, 6, 6, 6, 6, 5, 5, 5, 5,28,
     $          28,28,28,28,19,19,19, 9, 9, 9, 9, 9, 9, 9, 8,
     $           8, 8, 8, 8,10,10,10, 0, 2, 1, 1, 1, 3, 3, 3,
     $          15,15,16,16,13,13,13, 6, 5, 5, 5,28,28,28,28,
     $          28,28,11,19,19,19, 8, 8, 8, 8, 8,12,18,18,18,
     $          10,10,10,10,10, 0, 3, 3, 3, 3,15,16,17, 6, 6,
     $           6, 5,28,28,28,11,11,11,11,19,19,19,19,19, 8,
     $           8,18,10, 0, 0, 2, 2, 1, 3, 3, 3, 6, 6, 6,28,
     $          19, 8, 0, 1, 1, 3, 3,18, 0, 1, 0, 2, 1, 1, 3,
     $           3, 3,13,13, 6, 6, 6,28,28,28, 9, 9, 9, 9, 0,
     $           3, 0, 2, 3, 3, 3, 6, 6, 9, 9, 0, 2, 1, 1, 3,
     $           3, 3, 6, 8,12, 0, 9, 8, 8, 8, 0, 5, 5, 0, 1,
     $           1, 3,10, 9, 9, 9, 0/
c
      DATA Theta/6*180.0d0,3*120.0d0,121.7d0,120.0d0,121.0d0,116.4d0,
     $           6*120.0d0,122.0d0,120.0d0,123.0d0,120.0d0,118.0d0,
     $           120.0d0,120.0d0,118.0d0,121.8d0,116.4d0,120.0d0,
     $           117.0d0,120.0d0,123.0d0,120.0d0,120.0d0,117.0d0,
     $           123.0d0,5*120.0d0,123.0d0,123.0d0,120.0d0,114.0d0,
     $           120.0d0,110.5d0,120.0d0,125.6d0,111.5d0,125.0d0,
     $           5*109.5d0,109.47d0,4*109.5d0,110.0d0,109.5d0,
     $           109.5d0,110.0d0,8*109.5d0,110.0d0,11*109.5d0,
     $           112.0d0,3*109.5d0,3*107.8d0,109.5d0,107.8d0,
     $           13*120.0d0,118.0d0,4*120.0d0,118.0d0,9*120.0d0,
     $           180.0d0,120.0d0,120.0d0,123.0d0,110.0d0,123.0d0,
     $           110.0d0,120.0d0,112.0d0,118.0d0,118.0d0,120.0d0,
     $           120.0d0,105.0d0,109.5d0,110.0d0,109.5d0,118.0d0,
     $           118.0d0,3*109.5d0,120.0d0,120.0d0,118.0d0,122.0d0,
     $           120.0d0,118.0d0,120.0d0,119.0d0,117.0d0,120.0d0,
     $           120.0d0,109.5d0,16*120.0d0,112.0d0,127.0d0,109.5d0,
     $           110.0d0,109.5d0,109.5d0,3*110.0d0,108.5d0,103.9d0,
     $           120.0d0,5*109.5d0,110.0d0,111.0d0,111.0d0,97.0d0,
     $           94.3d0,98.0d0,97.5d0,102.9d0,107.0d0,107.0d0,
     $           118.0d0,109.5d0/
c
      DATA CIJK/0.040d0,0.040d0,0.040d0,0.040d0,0.040d0,0.040d0,0.024d0,
     $          0.020d0,0.036d0,0.018d0,0.024d0,0.024d0,0.046d0,0.024d0,
     $          0.026d0,0.024d0,0.024d0,0.036d0,0.036d0,0.030d0,0.012d0,
     $          0.070d0,0.024d0,0.020d0,0.040d0,0.024d0,0.040d0,0.030d0,
     $          0.030d0,0.024d0,0.020d0,0.040d0,0.070d0,0.030d0,0.024d0,
     $          0.020d0,0.070d0,0.060d0,0.026d0,0.026d0,0.026d0,0.026d0,
     $          0.030d0,0.030d0,0.072d0,0.030d0,0.030d0,0.014d0,0.030d0,
     $          0.028d0,0.030d0,0.016d0,0.020d0,0.018d0,0.024d0,0.018d0,
     $          0.024d0,0.018d0,0.024d0,0.018d0,0.020d0,0.020d0,0.024d0,
     $          0.040d0,0.016d0,0.016d0,0.024d0,0.018d0,0.018d0,0.024d0,
     $          0.018d0,0.022d0,0.018d0,0.020d0,0.020d0,0.020d0,0.040d0,
     $          0.040d0,0.018d0,0.020d0,0.020d0,0.022d0,0.022d0,0.018d0,
     $          0.020d0,0.020d0,0.014d0,0.018d0,0.018d0,0.040d0,0.018d0,
     $          0.018d0,0.018d0,0.024d0,0.020d0,0.024d0,0.036d0,0.024d0,
     $          0.024d0,0.024d0,0.036d0,0.036d0,0.036d0,0.040d0,0.040d0,
     $          0.040d0,0.062d0,0.062d0,0.040d0,0.030d0,0.040d0,0.040d0,
     $          0.024d0,0.040d0,0.040d0,0.040d0,0.040d0,0.040d0,0.040d0,
     $          0.040d0,0.062d0,0.062d0,0.062d0,0.080d0,0.040d0,0.040d0,
     $          0.080d0,0.082d0,0.080d0,0.082d0,0.040d0,0.044d0,0.040d0,
     $          0.040d0,0.044d0,0.044d0,0.044d0,0.040d0,0.040d0,0.018d0,
     $          0.040d0,0.040d0,0.040d0,0.010d0,0.018d0,0.020d0,0.018d0,
     $          0.044d0,0.040d0,0.052d0,0.044d0,0.044d0,0.016d0,0.020d0,
     $          0.018d0,0.024d0,0.044d0,0.018d0,0.024d0,0.052d0,0.024d0,
     $          0.020d0,0.024d0,0.020d0,0.020d0,0.040d0,0.040d0,0.040d0,
     $          0.040d0,0.040d0,0.040d0,0.018d0,0.018d0,0.040d0,0.060d0,
     $          0.020d0,0.020d0,0.044d0,0.044d0,0.020d0,0.020d0,0.020d0,
     $          0.044d0,0.094d0,0.010d0,0.020d0,0.020d0,0.014d0,0.020d0,
     $          0.020d0,0.040d0,0.040d0,0.040d0,0.020d0,0.022d0,0.020d0,
     $          0.062d0,0.060d0,0.040d0,0.040d0,0.040d0,0.020d0/
c
      DATA ISTART/53,7,96,1,138,125,124,176,0,194,165,186,0,0,0,0,0,
     $            191,167,8*0,146,199,200,144,0/
      DATA INUM/43,46,28,6,6,13,1,10,0,5,2,5,5*0,3,9,8*0,19,1,2,2,0/
C
C
C  ** WARNING! **  All parameters for Silicon made up
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
C  check for wild card fit
C
       IF(I1.EQ.13.OR.K1.EQ.13) THEN
        If(J1.EQ.2) Then
         L=21
         GO TO 95
        Else If(J1.EQ.1) Then
         L=65
         GO TO 95
        EndIf
       ELSE IF(I1.EQ.9.OR.K1.EQ.9) THEN
        If(J1.EQ.29) Then
         L=199
         GO TO 95
        Else If(J1.EQ.30) Then
         L=200
         GO TO 95
        EndIf
       ENDIF
c
       L=L1
       If(L.LT.199) GO TO 95
c
      ENDIF
C
C At this point nothing fits
C
      CosTh = (RD(I,J)*RD(I,J)+RD(J,K)*RD(J,K)-RD(I,K)*RD(I,K))/
     $           (TWO*RD(I,J)*RD(J,K))
      If(Abs(CosTh).GT.ONE) COSTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)*DEGREE
c
      RANG(NumB) = Th
      B3(NumB) = B3D
      RETURN
c
 95   CONTINUE
      RANG(NumB) = Theta(L)
      B3(NumB) = CIJK(L)
c
      RETURN
      END
*Deck fillb4
      SUBROUTINE FillB4(NAtoms,NTors,IAT,IC,IT1,IT2,IT3,IT4,IS4,B4)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),IS4(*),B4(*)
      DIMENSION IOP(6)
      COMMON/NUMINT/NumB,NumT,NumP
      DATA CP1/480.0d0/, CP2/120.0d0/
C
C  This Subroutine calculates the number of torsions and
C  out-ot-plane bends, storing the atoms involved in the
C  arrays IT1, IT2, IT3 and IT4, and fills the IS4 and B4
C  arrays by going through the relevant TRIPOS data
C  columns looking for a match
C
      NumT=0
      NumP=0
      If(NAtoms.LT.4) RETURN
c
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
     $ CALL DMATCH(NAtoms,K,I,J,L,K1,I1,J1,L1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(K,J).NE.0.AND.IC(I,L).NE.0)
     $ CALL DMATCH(NAtoms,K,J,I,L,K1,J1,I1,L1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(I,K).NE.0.AND.IC(K,L).NE.0)
     $ CALL DMATCH(NAtoms,J,I,K,L,J1,I1,K1,L1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(I,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMATCH(NAtoms,J,I,L,K,J1,I1,L1,K1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(J,K).NE.0.AND.IC(K,L).NE.0)
     $ CALL DMATCH(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(J,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMATCH(NAtoms,I,J,L,K,I1,J1,L1,K1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      ENDIF
c
      IF(IC(J,K).NE.0) THEN
C
C  J-K connected
C
      If(IC(K,I).NE.0.AND.IC(I,L).NE.0)
     $ CALL DMATCH(NAtoms,J,K,I,L,J1,K1,I1,L1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(I,K).NE.0.AND.IC(J,L).NE.0)
     $ CALL DMATCH(NAtoms,I,K,J,L,I1,K1,J1,L1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(I,L).NE.0.AND.IC(L,J).NE.0)
     $ CALL DMATCH(NAtoms,I,L,J,K,I1,L1,J1,K1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(I,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMATCH(NAtoms,I,L,K,J,I1,L1,K1,J1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      ENDIF
c
      IF(IC(I,K).NE.0) THEN
C
C  I-K connected
C
      If(IC(J,L).NE.0.AND.IC(L,I).NE.0)
     $ CALL DMATCH(NAtoms,J,L,I,K,J1,L1,I1,K1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      If(IC(K,L).NE.0.AND.IC(L,J).NE.0)
     $ CALL DMATCH(NAtoms,I,K,L,J,I1,K1,L1,J1,IC,IT1,IT2,IT3,IT4,IS4,B4)
      ENDIF
c
c -- check number of torsions does not exceed maximum allowed
      If(NumT.GT.NTors) Call nerror(3,'Force Field module',
     $           'More Torsions in System than Dimensioned for',0,0)
c
 10   CONTINUE
C
C
      DO 20 I=1,NAtoms
      I1=IAT(I)
C
C  go through each atom and see it forms an out-of-plane bend
C
      IF(I1.EQ.2.OR.I1.EQ.3) THEN
       CP=CP1
       GO TO 21
      ELSE IF(I1.EQ. 6.OR.I1.EQ.11.OR.
     $        I1.EQ.19.OR.I1.EQ.28) THEN
       CP=CP2
       GO TO 21
      ELSE
       GO TO 20
      ENDIF
c
 21   CONTINUE
C
C  atom type is suitable for an out-of-plane-bend
C  check if there are connected atoms
C
      NC=0
      DO 22 J=1,NAtoms
      IF(IC(I,J).NE.0) THEN
       NC=NC+1
       IOP(NC)=J
      ENDIF
 22   CONTINUE
c
      IF(NC.EQ.3) THEN
       NumP=NumP+1
       Num=NumP+NumT
       IT1(Num)=I
       IT2(Num)=IOP(1)
       IT3(Num)=IOP(2)
       IT4(Num)=IOP(3)
       B4(Num)=CP
      ENDIF
c
c -- check number of outp does not exceed maximum allowed
      If(NumT+NumP.GT.NTors) Call nerror(4,'Force Field module',
     $  'More Out-of-Plane Bends in System than Dimensioned for',0,0)
c
 20   CONTINUE

      RETURN
      END
*Deck dmatch
      SUBROUTINE DMATCH(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,
     $                  IT3,IT4,IS4,B4)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),IS4(*),B4(*)
      DIMENSION IATI(76),IATJ(76),IATK(76),IATL(76)
      DIMENSION IBT(76),CIJKL(76),IS(76)
      DIMENSION ISTART(3),INUM(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ISD=3,B4D=0.2)
c
      DATA IATI/ 9, 2, 2, 1, 2, 2, 1, 1, 1,13,13,13, 2, 2, 1,
     $          13,60*0/
c
      DATA IATJ/ 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $           2, 2, 2, 2, 1, 4, 4, 4, 4, 2, 2, 4, 2, 1, 4,
     $           2, 1, 3, 3, 4, 4, 2, 2, 1, 3, 6, 6, 2, 1, 3,
     $           5, 2, 2, 1, 3, 6, 5,28, 3, 2, 1, 3, 6,19, 2,
     $           1, 3, 6, 5, 2, 1, 3, 8, 2, 1, 3, 5, 2, 1, 3,10/
c
      DATA IATK/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $           1, 1, 1, 1, 1, 4, 4, 2, 2, 2, 2, 1, 1, 1, 3,
     $           3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5,
     $           5,28,28,28,28,28,28,28,11,19,19,19,19,19, 8,
     $           8, 8, 8, 8,12,12,12,12,18,18,18,18,10,10,10,10/
c
      DATA IATL/ 1, 2, 1, 1, 2,13, 2, 1,13, 2, 1,13, 1, 0, 0,
     $           0, 2, 1,13,13,56*0/
c
      DATA IBT/20*1,3,1,1,2,2,7*1,5,1,1,2,2,3*1,2,5*1,4,6*1,
     $            5,14*1,2,7*1/
c
      DATA CIJKL/0.7d0,0.04d0,0.126d0,0.5d0,0.126d0,0.273d0,0.126d0,
     $           0.126d0,0.274d0,0.274d0,0.274d0,0.274d0,0.126d0,
     $           0.126d0,0.126d0,0.274d0,0.126d0,0.126d0,0.274d0,
     $           0.32d0,0.0d0,0.0d0,0.0d0,0.0d0,12.5d0,1.424d0,0.0d0,
     $           0.12d0,0.2d0,0.0d0,1.6d0,0.12d0,2.0d0,0.6d0,0.0d0,
     $           0.0d0,12.0d0,12.0d0,0.4d0,1.6d0,1.6d0,1.6d0,0.12d0,
     $           0.2d0,0.12d0,0.2d0,6.46d0,6.46d0,0.2d0,1.6d0,1.6d0,
     $           0.12d0,1.6d0,1.6d0,12.0d0,0.4d0,1.6d0,1.6d0,1.6d0,
     $           5.8d0,1.2d0,1.2d0,1.0d0,0.2d0,1.0d0,0.4d0,1.0d0,0.4d0,
     $           1.0d0,0.4d0,1.0d0,0.4d0,1.0d0,0.4d0,1.0d0,4.0d0/
c
      DATA   IS/-3, 3, 3, 3,-3,-3, 3, 3, 3, 3, 3, 3,-3,-3, 3,
     $           3, 3, 3, 3, 3, 1, 1, 1, 1,-2,-2, 1,-3, 3, 1,
     $          -2,-3,-2,-2, 1, 1,-2,-2,-3,-2,-2,-2,-3, 3,-3,
     $           3,-2,-2, 3,-2,-2,-3,-2,-2,-2,-3,-2,-2,-2,-2,
     $           3,-2, 2, 3,-2, 3, 3, 3,-2, 3, 3, 3,-2, 3, 3, 3/
c
      DATA ISTART/1,14,21/, INUM/13,7,56/
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
C  2. check for single dummy fits
C
      M1=ISTART(2)
      M2=M1+INUM(2)-1
      DO 20 M=M1,M2
      IF(J1.EQ.IATJ(M).AND.K1.EQ.IATK(M)) THEN
C
C  match found for central atoms
C  check if everything else fits
C
       If( ((IATI(M).EQ.0.AND.L1.EQ.IATL(M)).OR.(I1.EQ.IATI(M).AND.
     $       IATL(M).EQ.0)).AND.IC(J,K).EQ.IBT(M) ) GO TO 95
      ENDIF
c
      IF(K1.EQ.IATJ(M).AND.J1.EQ.IATK(M)) THEN
C
C  match found for central atoms
C  check if everything else fits
C
       If( ((IATI(M).EQ.0.AND.I1.EQ.IATL(M)).OR.(L1.EQ.IATI(M).AND.
     $       IATL(M).EQ.0)).AND.IC(J,K).EQ.IBT(M) ) GO TO 95
      ENDIF
 20   CONTINUE
C
C  3. check for double dummy fits
C
      M1=ISTART(3)
      M2=M1+INUM(3)-1
      DO 30 M=M1,M2
      If(J1.EQ.IATJ(M).AND.K1.EQ.IATK(M)) GO TO 95
      If(K1.EQ.IATJ(M).AND.J1.EQ.IATK(M)) GO TO 95
 30   CONTINUE
C
C  At this point nothing fits
C
      B4(NumT) = B4D
      IS4(NumT)= ISD
      RETURN
c
 95   CONTINUE
      B4(NumT) = CIJKL(M)
      IS4(NumT)= IS(M)
c
      RETURN
      END
