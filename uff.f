c ==================================================================
c  UNIVERSAL FORCE FIELD         JB   May 2005
c ==================================================================
c
      SUBROUTINE UFF_MAIN(NAtoms, IAN,    IType,  NBend,  NTors,
     $                    NInv,   XC,     IC,     IAT,    RD,
     $                    RLEN,   B2,     IB1,    IB2,    IB3,
     $                    C0,     C1,     C2,     B3,     NB,
     $                    IT1,    IT2,    IT3,    IT4,    CosD,
     $                    B4,     NT,     IP1,    IP2,    IP3,
     $                    IP4,    S0,     S1,     S2,     B5,
     $                    GU,     GC,     HESS,   ffcyc)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Preliminary Driver for universal force field
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
C  NTors   -  maximum number of torsions
C  NInv    -  maximum number of inversions
C  XC      -  Cartesian coordinates
C  IC      -  atomic connectivity matrix
C  IAT     -  UFF atom type
C  RD      -  distance matrix
C  RLEN    -  UFF equilibrium/van der Waals bond lengths
C  B2      -  UFF stretching/van der Waals force constants
C  IB1     -  1st atom in bend
C  IB2     -  2nd atom in bend
C  IB3     -  3rd atom in bend
C  C0      -  1st bending coefficient
C  C1      -  2nd bending coefficient
C  C2      -  3rd bending coefficient
C  B3      -  UFF bending force constants
C  NB      -  integer bend order
C  IT1     -  1st atom in torsion/out-of-plane bend
C  IT2     -  2nd atom in torsion/out-of-plane bend
C  IT3     -  3rd atom in torsion/out-of-plane bend
C  IT4     -  4th atom in torsion/out-of-plane bend
C  CosD    -  Cosine factor for torsion
C  B4      -  UFF torsional force constants
C  NT      -  integer torsion order
C  IP1     -  1st atom in inversion/out-of-plane bend
C  IP2     -  2nd atom in inversion/out-of-plane bend
C  IP3     -  3rd atom in inversion/out-of-plane bend
C  IP4     -  4th atom in inversion/out-of-plane bend
C  S0      -  1st inversion coefficient
C  S1      -  2nd inversion coefficient
C  S2      -  3rd inversion coefficient
C  B5      -  UFF inversion force constants
C  GU      -  step-up gradient for finite difference Hessian
C  GC      -  Cartesian gradient (step-down gradient)
C  HESS    -  Cartesian Hessian matrix
C  ffcyc   -  integer; controls force field initialization
C
C
      DIMENSION IAN(NAtoms),XC(3*NAtoms),IC(NAtoms,NAtoms),
     $          IAT(NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms),
     $          B2(NAtoms,NAtoms),IB1(NBend),IB2(NBend),IB3(NBend),
     $          C0(NBend),C1(NBend),C2(NBend),B3(NBend),NB(NBend),
     $          IT1(NTors),IT2(NTors),IT3(NTors),IT4(NTors),
     $          CosD(NTors),B4(NTors),NT(NTors),
     $          IP1(NInv),IP2(NInv),IP3(NInv),IP4(NInv),
     $          S0(NInv),S1(NInv),S2(NInv),B5(NInv)
      DIMENSION GU(3*NAtoms),GC(3*NAtoms),HESS(3*NAtoms,3*NAtoms)
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
      READ(JUnit) IType,IAN,IC,IAT,RD,RLEN,B2,IB1,IB2,IB3,C0,C1,C2,B3,
     $            NB,IT1,IT2,IT3,IT4,CosD,B4,NT,IP1,IP2,IP3,IP4,S0,S1,
     $            S2,B5,NumB,NumT,NumP
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
      CALL UFF_Field(NAtoms, IAN,    IType,  NBend,  NTors,
     $               NInv,   CutOff, IPRNT,  XC,     IC,
     $               IAT,    RD,     RLEN,   B2,     IB1,
     $               IB2,    IB3,    C0,     C1,     C2,
     $               B3,     NB,     IT1,    IT2,    IT3,
     $               IT4,    CosD,   B4,     NT,     IP1,
     $               IP2,    IP3,    IP4,    S0,     S1,
     $               S2,     B5,     GU,     E,      GC,
     $               HESS,   ffcyc,  IStatus)
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
      WRITE(JUnit) IType,IAN,IC,IAT,RD,RLEN,B2,IB1,IB2,IB3,C0,C1,C2,B3,
     $             NB,IT1,IT2,IT3,IT4,CosD,B4,NT,IP1,IP2,IP3,IP4,S0,S1,
     $             S2,B5,NumB,NumT,NumP
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
c --  write to the <control> file info about hessian quality
        ihessq=+1
        call wrcntrl(IUnit,9,'$hessqual',1,ihessq,rdum,cdum)
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
*Deck uff-field
      SUBROUTINE UFF_Field(NAtoms, IAN,    IType,  NBend,  NTors,
     $                     NInv,   CutOff, IPRNT,  XC,     IC,
     $                     IAT,    RD,     RLEN,   B2,     IB1,
     $                     IB2,    IB3,    C0,     C1,     C2,
     $                     B3,     NB,     IT1,    IT2,    IT3,
     $                     IT4,    CosD,   B4,     NT,     IP1,
     $                     IP2,    IP3,    IP4,    S0,     S1,
     $                     S2,     B5,     GU,     E,      G,
     $                     HESS,   ffcyc,  IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main Driver for UFF force field
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
C  NTors   -  maximum number of torsions
C  NInv    -  maximum number of inversions
C  CutOff  -  distance cutoff to limit Van der Waals interactions
C  IPRNT   -  debug print flag
C  XC      -  Cartesian coordinates
C  IC      -  atomic connectivity matrix
C  IAT     -  UFF atom type
C  RD      -  distance matrix
C  RLEN    -  UFF equilibrium/vam der Waals bond lengths
C  B2      -  UFF stretching/van der Waals force constants
C  IB1     -  1st atom in bend
C  IB2     -  2nd atom in bend
C  IB3     -  3rd atom in bend
C  C0      -  1st bending coefficient
C  C1      -  2nd bending coefficient
C  C2      -  3rd bending coefficient
C  B3      -  UFF bending force constants
C  NB      -  integer bend order
C  IT1     -  1st atom in torsion/out-of-plane bend
C  IT2     -  2nd atom in torsion/out-of-plane bend
C  IT3     -  3rd atom in torsion/out-of-plane bend
C  IT4     -  4th atom in torsion/out-of-plane bend
C  CosD    -  Cosine factor for torsion
C  B4      -  UFF torsional force constants
C  NT      -  integer torsion order
C  IP1     -  1st atom in inversion/out-of-plane bend
C  IP2     -  2nd atom in inversion/out-of-plane bend
C  IP3     -  3rd atom in inversion/out-of-plane bend
C  IP4     -  4th atom in inversion/out-of-plane bend
C  S0      -  1st inversion coefficient
C  S1      -  2nd inversion coefficient
C  S2      -  3rd inversion coefficient
C  B5      -  UFF inversion force constants
C  GU      -  gradient for finite-difference Hessian
C  E       -  UFF energy
C  G       -  Cartesian gradient
C  HESS    -  Cartesian Hessian matrix
C  ffcyc   -  integer; controls force field initialization
C  IStatus -  status on exit
C              0 - successful calculation
C             -1 - something went wrong
C
C
      DIMENSION IAN(NAtoms),XC(3*NAtoms),IC(NAtoms,NAtoms),
     $          IAT(NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms),
     $          B2(NAtoms,NAtoms),IB1(NBend),IB2(NBend),IB3(NBend),
     $          C0(NBend),C1(NBend),C2(NBend),B3(NBend),IB(NBend),
     $          IT1(NTors),IT2(NTors),IT3(NTors),IT4(NTors),
     $          CosD(NTors),B4(NTors),NT(NTors),
     $          IP1(NInv),IP2(NInv),IP3(NInv),IP4(NInv),
     $          S0(NInv),S1(NInv),S2(NInv),B5(NInv)
      DIMENSION GU(3*NAtoms),G(3*NAtoms),HESS(3*NAtoms,3*NAtoms)
      CHARACTER jobname*256
      INTEGER ffcyc
      LOGICAL ReadPQB
      SAVE ReadPQB
c
      PARAMETER (ZERO=0.0D0,delta=0.001D0,HALF=0.5D0)
      PARAMETER (ANTOAU=1.889725991249198d0,CVRT=0.0015936D0)
      DATA ReadPQB /.False./
c
      Common /job/jobname,lenJ
C
C
C  THIS ROUTINE CAN CALCULATE THE ENERGY, GRADIENT AND HESSIAN
C  of a molecular system using classical model potential
C  functions taken from the UFF force field.
C  The Gradient is calculated analytically by differentiating
C  the energy with respect to Cartesian displacements and the
C  Hessian is estimated by central differences on the gradient.
C
C  ..............................................................
C     The UFF force field was coded from the original literature
C
C     A.K.Rappe, C.J.Casewit, K.S.Colwell, W.A.Goddard III and W.M.Skiff
C     J.Am.Chem.Soc. 114 (1992) 10024
C
C     Thanks also to Marcus G. Martin, Sandia Laboratories, for
C     supplying unpublished electronegativity parameters from his
C     Towhee molecular mechanics package and for useful discussions
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
        CALL RdPQBFF(40,'UFF',NAtoms,Nqm,IAT,IC)
        CALL AddLink('UFF',NAtoms,NLink,XC,IAT,IC)
       Else
        CALL RdPQBFF(40,'UFF',NAtoms,NAtoms,IAT,IC)
        ReadPQB = .True.
       EndIf
c
       If(IPRNT.GT.2) Then
         write(6,1100)
         do i=1,NAtoms
         write(6,1200) i,IAT(i)
         enddo
         write(6,1300)
         do i=2,NAtoms
         do j=1,i-1
         if(IC(i,j).NE.0) write(6,1400) i,j,IC(i,j)
         enddo
         enddo
       EndIf
c
       CLOSE (UNIT=40,STATUS='KEEP')
       GO TO 15
c
 10    CONTINUE
       CALL UFF_CONNECTF(NAtoms, IPRNT,  IAN,    XC,     RD,
     $                   IT1,    IT2,    IAT,    IC,     IStatus)
       If(IStatus.LT.0) GO TO 99
C
C  Fill the UFF interaction and parameter arrays
C
 15    CONTINUE
       CALL FillUFF(NAtoms, NBend,  NTors,  NInv,   IPRNT,
     $              IAT,    IC,     RLEN,   B2,     IB1,
     $              IB2,    IB3,    C0,     C1,     C2,
     $              B3,     NB,     IT1,    IT2,    IT3,
     $              IT4,    CosD,   B4,     NT,     IP1,
     $              IP2,    IP3,    IP4,    S0,     S1,
     $              S2,     B5)
c
       If(IType.GE.10)
     $    CALL UFF_PREAMBLE(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $                      IC,     RD,     RLEN,   B2,     IB1,
     $                      IB2,    IB3,    C0,     C1,     C2,
     $                      B3,     NB,     IT1,    IT2,    IT3,
     $                      IT4,    CosD,   B4,     NT,     IP1,
     $                      IP2,    IP3,    IP4,    S0,     S1,
     $                      S2,     B5,     G,      GU)
c
      ENDIF
c
      ffcyc = ffcyc + 1
c
      IF(IType.EQ.1.OR.Itype.EQ.11) THEN
C
C  form the UFF  Hessian by central difference on the gradient
C
       If(IPRNT.GT.1)
     $    WRITE(6,*) '  Hessian Calculated Using UFF Force Field'
c
       CALL ZeroIT(HESS,NAT3*NAT3)
c
       DO 20 I=1,NAT3
C
C  step up
C
       XC(I) = XC(I) + delta
c
       CALL EGUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $            IC,     RD,     RLEN,   B2,     IB1,
     $            IB2,    IB3,    C0,     C1,     C2,
     $            B3,     NB,     IT1,    IT2,    IT3,
     $            IT4,    CosD,   B4,     NT,     IP1,
     $            IP2,    IP3,    IP4,    S0,     S1,
     $            S2,     B5,     EU,     GU)
C
C  step down
C
       XC(I) = XC(I) - delta - delta
c
       CALL EGUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $            IC,     RD,     RLEN,   B2,     IB1,
     $            IB2,    IB3,    C0,     C1,     C2,
     $            B3,     NB,     IT1,    IT2,    IT3,
     $            IT4,    CosD,   B4,     NT,     IP1,
     $            IP2,    IP3,    IP4,    S0,     S1,
     $            S2,     B5,     ED,     G)
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
cc       write(6,*) ' Hessian is:'
cc       call prntmat(nat3,nat3,nat3,hess)
      ENDIF
c
      IF(IType.GE.0) THEN
C
C  calculate energy and gradient
C
       CALL EGUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $            IC,     RD,     RLEN,   B2,     IB1,
     $            IB2,    IB3,    C0,     C1,     C2,
     $            B3,     NB,     IT1,    IT2,    IT3,
     $            IT4,    CosD,   B4,     NT,     IP1,
     $            IP2,    IP3,    IP4,    S0,     S1,
     $            S2,     B5,     E,      G)
       CALL VScal(NAT3,-CVRT/ANTOAU,G)
c
       If(IPRNT.GT.4) Then
         write(6,*) ' Force Field Cartesian Gradient (au) is:'
         do i=1,natoms
         ii = 3*(i-1)
         write(6,1234) i,g(ii+1),g(ii+2),g(ii+3)
 1234    format(1X,i4,2X,3(F15.7))
         enddo
       EndIf
c
      ENDIF
c
 99   CALL VScal(NAT3,ANTOAU,XC)
C
C  -- ERROR HANDLING --
      If(IStatus.LT.0) Call nerror(5,'Force Field module',
     $  'Unable To Fully Assign Atomic Connectivity/Bond Orders',0,0)
C
      RETURN
c
 1100 FORMAT(/,2X,' UFF Atom types',//,
     $         3X,' ATOM  TYPE',/,
     $         3X,' ----  ----')
 1200 FORMAT (4X,I3,2X,I4)
 1300 FORMAT (/,2X,'BONDS info',//,
     $          3X,' Atom1  Atom2  Bond Type',/,
     $          3X,' -----  -----  ---------')
 1400 FORMAT(4X,I4,4X,I4,4X,I4)
c
      END
*Deck euff
      SUBROUTINE EUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $                IC,     RD,     RLEN,   B2,     IB1,
     $                IB2,    IB3,    C0,     C1,     C2,
     $                B3,     NB,     IT1,    IT2,    IT3,
     $                IT4,    CosD,   B4,     NT,     IP1,
     $                IP2,    IP3,    IP4,    S0,     S1,
     $                S2,     B5,     E)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XC(3,NAtoms),RD(NAtoms,NAtoms)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      REAL*8 B3(*),B4(*),B5(*)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IB1(*),IB2(*),IB3(*),C0(*),C1(*),C2(*),NB(*)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),CosD(*),NT(*)
      DIMENSION IP1(*),IP2(*),IP3(*),IP4(*),S0(*),S1(*),S2(*)
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (THREE=3.0d0,FOUR=4.0D0,TWELVE=12.0D0)
      PARAMETER (HalfPI=0.5d0*3.14159265357979D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
      PARAMETER (TOLLZERO=1.0D-8)
C
C
C   THIS ROUTINE CALCULATES THE ENERGY
C   of a molecule using the UFF forcefield
C   Energy contributions coded are:
C
C   1. Bond Stretching
C      ES = 1/2 Kij(R-Rij)**2
C   R is the current interatomic distance in angstroms.
C   Rij is the sum of standard radii for each atom, plus a
C   bond-order correction, plus an electronegativity correction
C
C      Rij = Ri + Rj - Rbo - Ren
C
C   where
C
C      Rbo = 0.1332 (Ri+Rj) ln(n)   n = bond order
C   (note: C-N amide bond order is 1.41; bond order in aromatic rings is 1.5)
C
C      Ren = Ri*Rj (SQRT(Xi) - SQRT(Xj))**2 / (Xi*Ri + Xj*Rj)
C   (Xi is the GMP electronegativity of atom i)
C
C   The stretching force constants Kij are atom-based and are obtained from a
C   generalization of Badger's rules. They are given by
C
C      Kij = 664.12 (Zi*Zj)/(Rij)**3
C
C   where Zi, Zj are effective atomic charges
C
C   2. Non-Bonded VdW interaction (1,3 excluded)
C      EV = Dij[ -2(Xij/X)**6 + (Xij/X)**12 ]
C   Dij is the well depth and Xij is the van der Waals bond length
C   X is the actual distance between atoms i and j
C
C   Dij = SQRT(Di*Dj)       Di,Dj are for individual atoms
C   Xij = SQRT(Xi*Xj)       Xi,Xj are for individual atoms
C
C   3. Bending
C    (a) general nonlinear bend
C      EB = Kijk(C0+C1*Cos(Th)+C2*Cos(2*Th))
C   Th is the current bond angle and Kijk is the bending force
C   constant between atoms I,J,K. The three expansion coefficients
C   are given by
C
C      C2 = 1/(4(Sin(Th0))**2)
C      C1 = -4 C2 Cos(Th0)
C      C0 = C2 (2(Cos(Th0))**2 +1)
C   where Th0 is the idealized equilibrium angle at the central atom (j)
C
C    (b) linear, trigonal-planar, square-planar and octahedral
C      EB = (Kijk/n**2) (1-cos(n*Th))
C   (linear: n=1; trigonal-planar: n=3; square-planar/octahedral n=4)
C
C   The bending force constants in both cases are given by
C
C      Kijk = (664.12/Rij*Rjk) (Zi*Zk/Rik**5) (Rij*Rjk)
C                                [3*Rij*Rjk(1-(Cos(Th0))**2) - Rik Cos(Th0)]
C
C   4. Torsion
C      ET = 1/2 V ( 1 - Cos(n*Dih0) Cos(n*Dih) )
C   V=Kijkl is the torsional barrier and Dih0 is the idealized dihedral angle
C   Specific general cases include
C
C    (a) j=sp3 hybridized center; k=sp3 hybridized center
C        where  n=3 and Dih0=180 (or 60) and
C
C      V = SQRT(Vj*Vk)   with   Vi being an individual main group atom value
C                               (non main-group elements are assigned V=0)
C
C   The torsional terms for pairs of sp3 hybridized group 6 central atoms are exceptions
C   Here  Vj=2 for oxygen and Vj=6.8 for the remaining group elements
C   with n=2 and Dih0=90
C
C    (b) j=sp2 hybridized center; k=sp3 hybridized center
C        where  n=6 and Dih0=0 and V=1.0
C
C   For a single bond involving an sp3 hybridized group 6 central atom
C   and an sp2 atom of another column, V is defined as in (c), below,
C   with n=2 and Dih0=90
C
C   For a single bond where the sp2 hydridized center in (b) is connected
C   to another sp2 hybridized center (e.g., propene) then
C   n=3 and Dih0=180 and V=2.0
C
C    (c) j=sp2 hybridized center; k=sp2 hybridized center
C        where n=2 and Dih0 = 180 (or 60??) and
C
C      V = 5 SQRT(Uj*Uk) (1 + 4.18*ln(BO))  with  BO the bond order between i and j
C      (Uj constants take values 2, 1.25, 0.7, 0.2 and 0.1 for the second through
C      the sixth period/first through the fifth rows)
C
C   5. Inversion
C      EP = Kijkl[C0 + C1 Cos(Yijkl) + C2 Cos(2*Yijkl)]
C   This term applies to exactly 3 atoms (j,k,l) bonded to a central atom (i)
C   Yijkl is the angle between the il axis and the ijk plane
C
C   For C.2 and C.R sp2 atom types   C0=1; C1=-1; c2=0
C   If carbon is bonded to O.2  Kijkl=50  otherwise  Kijkl=6
C
C   Central atoms for which as inversion term is defined are: C,N,P,As,Sb,Bi
C   The inversion force constants (Kijkl) for all other atoms are set to zero
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
      RD(I,J) = SQRT( (XC(1,I)-XC(1,J))**2 + (XC(2,I)-XC(2,J))**2
     $              + (XC(3,I)-XC(3,J))**2 )
      RD(J,I) = RD(I,J)
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
C     (For 1,3-interactions B2(I,J) has been
C      set to zero in routine UFFB2)
C
        RV = RLEN(I,J)/R
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
      n=NB(M)
c
      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
c
      CosTh = (RIJ*RIJ+RJK*RJK-RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosTh).GT.ONE) CosTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)
cc      write(6,*) ' Bend is ',Th*degree
c
      If(n.EQ.0) Then
        EB = EB + B3(M)*(C0(M) + C1(M)*CosTh + C2(M)*DCOS(TWO*Th))
      Else If(n.EQ.1) Then
c -- linear
        EB = EB + B3(M)*(ONE + CosTh)
      Else
        EB = EB + (B3(M)/DFloat(n**2))*(ONE - DCOS(n*Th))
      EndIf
c
      SinTh = SQRT(ONE - CosTh**2)
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
      n=NT(M)

      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
      RJL=RD(J,L)
      RKL=RD(K,L)
c
      XIJ =XC(1,I) - XC(1,J)
      XIK =XC(1,I) - XC(1,K)
      XJK =XC(1,J) - XC(1,K)
      XJL =XC(1,J) - XC(1,L)
      XKL =XC(1,K) - XC(1,L)
      YIJ =XC(2,I) - XC(2,J)
      YIK =XC(2,I) - XC(2,K)
      YJK =XC(2,J) - XC(2,K)
      YJL =XC(2,J) - XC(2,L)
      YKL =XC(2,K) - XC(2,L)
      ZIJ =XC(3,I) - XC(3,J)
      ZIK =XC(3,I) - XC(3,K)
      ZJK =XC(3,J) - XC(3,K)
      ZJL =XC(3,J) - XC(3,L)
      ZKL =XC(3,K) - XC(3,L)
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
      ET = ET + B4(M)*(ONE - CosD(M)*DCOS(n*Dih))
c
 40   CONTINUE
C
C
C  5. Inversions/Out-of-Plane bends
C
      DO 50 M=1,NumP
      I=IP1(M)
      J=IP2(M)
      K=IP3(M)
      L=IP4(M)
c
      XIJ = XC(1,I) - XC(1,J)
      XIK = XC(1,I) - XC(1,K)
      XIL = XC(1,I) - XC(1,L)
      YIJ = XC(2,I) - XC(2,J)
      YIK = XC(2,I) - XC(2,K)
      YIL = XC(2,I) - XC(2,L)
      ZIJ = XC(3,I) - XC(3,J)
      ZIK = XC(3,I) - XC(3,K)
      ZIL = XC(3,I) - XC(3,L)

      R1(1) = -XIJ
      R1(2) = -YIJ
      R1(3) = -ZIJ
      R2(1) = -XIK
      R2(2) = -YIK
      R2(3) = -ZIK
      R3(1) = -XIL
      R3(2) = -YIL
      R3(3) = -ZIL
      CALL VECMUL(R1,R2,R4)
c
      dn = SProd(3,R4,R3)
      dm = SQRT(SProd(3,R4,R4))*SQRT(SProd(3,R3,R3))
      CosTh = dn/dm
      If(Abs(CosTh).GT.ONE) CosTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)
c -- this is the angle IL makes with the normal
c -- the angle we want is 90 degrees minus this
      Th = HalfPI-Th
cc      write(6,*) ' Out-of-plane bend is ',Th*degree
c
      EP = EP + B5(M)*(S0(M) + S1(M)*DCOS(Th) + S2(M)*DCOS(TWO*Th))
c
 50   CONTINUE
c
 95   CONTINUE
C
C  Form the total molecular energy
C
      ES = HALF*ES
c
      E = ES + EV + EB + ET + EP
c
      if(iprnt.gt.2) write(6,1000) E,ES,EB,EP,ET,EV
      RETURN
c
 1000 FORMAT(/,' TOTAL ENERGY       ',f15.8,/,
     $         ' Stretching         ',f15.8,/,
     $         ' Bending            ',f15.8,/,
     $         ' Out-of-Plane Bend  ',f15.8,/,
     $         ' Torsion            ',f15.8,/,
     $         ' VdW energy         ',f15.8)
c
      END
*Deck eguff
      SUBROUTINE EGUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $                 IC,     RD,     RLEN,   B2,     IB1,
     $                 IB2,    IB3,    C0,     C1,     C2,
     $                 B3,     NB,     IT1,    IT2,    IT3,
     $                 IT4,    CosD,   B4,     NT,     IP1,
     $                 IP2,    IP3,    IP4,    S0,     S1,
     $                 S2,     B5,     E,      G)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XC(3,NAtoms),RD(NAtoms,NAtoms)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      REAL*8 B3(*),B4(*),B5(*),G(3,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IB1(*),IB2(*),IB3(*),C0(*),C1(*),C2(*),NB(*)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),CosD(*),NT(*)
      DIMENSION IP1(*),IP2(*),IP3(*),IP4(*),S0(*),S1(*),S2(*)
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (THREE=3.0d0,FOUR=4.0D0,TWELVE=12.0D0)
      PARAMETER (HalfPI=0.5d0*3.14159265357979D0)
      PARAMETER (DEGREE=180.0D0/3.14159265357979D0)
      PARAMETER (TOLLZERO=1.0D-8)
C
C
C   THIS ROUTINE CALCULATES THE ENERGY AND CARTESIAN GRADIENT
C   of a molecule using the UFF forcefield
C   Energy contributions coded are:
C
C   1. Bond Stretching
C      ES = 1/2 Kij(R-Rij)**2
C   R is the current interatomic distance in angstroms.
C   Rij is the sum of standard radii for each atom, plus a
C   bond-order correction, plus an electronegativity correction
C
C      Rij = Ri + Rj - Rbo - Ren
C
C   where
C
C      Rbo = 0.1332 (Ri+Rj) ln(n)   n = bond order
C   (note: C-N amide bond order is 1.41; bond order in aromatic rings is 1.5)
C
C      Ren = Ri*Rj (SQRT(Xi) - SQRT(Xj))**2 / (Xi*Ri + Xj*Rj)
C   (Xi is the GMP electronegativity of atom i)
C
C   The stretching force constants Kij are atom-based and are obtained from a
C   generalization of Badger's rules. They are given by
C
C      Kij = 664.12 (Zi*Zj)/(Rij)**3
C
C   where Zi, Zj are effective atomic charges
C
C   2. Non-Bonded VdW interaction (1,3 excluded)
C      EV = Dij[ -2(Xij/X)**6 + (Xij/X)**12 ]
C   Dij is the well depth and Xij is the van der Waals bond length
C   X is the actual distance between atoms i and j
C
C   Dij = SQRT(Di*Dj)       Di,Dj are for individual atoms
C   Xij = SQRT(Xi*Xj)       Xi,Xj are for individual atoms
C
C   3. Bending
C    (a) general nonlinear bend
C      EB = Kijk(C0+C1*Cos(Th)+C2*Cos(2*Th))
C   Th is the current bond angle and Kijk is the bending force
C   constant between atoms I,J,K. The three expansion coefficients
C   are given by
C
C      C2 = 1/(4(Sin(Th0))**2)
C      C1 = -4 C2 Cos(Th0)
C      C0 = C2 (2(Cos(Th0))**2 +1)
C   where Th0 is the idealized equilibrium angle at the central atom (j)
C
C    (b) linear, trigonal-planar, square-planar and octahedral
C      EB = (Kijk/n**2) (1-cos(n*Th))
C   (linear: n=1; trigonal-planar: n=3; square-planar/octahedral n=4)
C
C   The bending force constants in both cases are given by
C
C      Kijk = (664.12/Rij*Rjk) (Zi*Zk/Rik**5) (Rij*Rjk)
C                                [3*Rij*Rjk(1-(Cos(Th0))**2) - Rik Cos(Th0)]
C
C   4. Torsion
C      ET = 1/2 V ( 1 - Cos(n*Dih0) Cos(n*Dih) )
C   V=Kijkl is the torsional barrier and Dih0 is the idealized dihedral angle
C   Specific general cases include
C
C    (a) j=sp3 hybridized center; k=sp3 hybridized center
C        where  n=3 and Dih0=180 (or 60) and
C
C      V = SQRT(Vj*Vk)   with   Vi being an individual main group atom value
C                               (non main-group elements are assigned V=0)
C
C   The torsional terms for pairs of sp3 hybridized group 6 central atoms are exceptions
C   Here  Vj=2 for oxygen and Vj=6.8 for the remaining group elements
C   with n=2 and Dih0=90
C
C    (b) j=sp2 hybridized center; k=sp3 hybridized center
C        where  n=6 and Dih0=0 and V=1.0
C
C   For a single bond involving an sp3 hybridized group 6 central atom
C   and an sp2 atom of another column, V is defined as in (c), below,
C   with n=2 and Dih0=90
C
C   For a single bond where the sp2 hydridized center in (b) is connected
C   to another sp2 hybridized center (e.g., propene) then
C   n=3 and Dih0=180 and V=2.0
C
C    (c) j=sp2 hybridized center; k=sp2 hybridized center
C        where n=2 and Dih0 = 180 (or 60??) and
C
C      V = 5 SQRT(Uj*Uk) (1 + 4.18*ln(BO))  with  BO the bond order between i and j
C      (Uj constants take values 2, 1.25, 0.7, 0.2 and 0.1 for the second through
C      the sixth period/first through the fifth rows)
C
C   5. Inversion
C      EP = Kijkl[C0 + C1 Cos(Yijkl) + C2 Cos(2*Yijkl)]
C   This term applies to exactly 3 atoms (j,k,l) bonded to a central atom (i)
C   Yijkl is the angle between the il axis and the ijk plane
C
C   For C.2 and C.R sp2 atom types   C0=1; C1=-1; c2=0
C   If carbon is bonded to O.2  Kijkl=50  otherwise  Kijkl=6
C
C   Central atoms for which as inversion term is defined are: C,N,P,As,Sb,Bi
C   The inversion force constants (Kijkl) for all other atoms are set to zero
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
      RD(1,1) = ZERO
      DO 10 I=2,NAtoms
      RD(I,I) = ZERO
      DO 10 J=1,I-1
      RD(I,J) = SQRT( (XC(1,I)-XC(1,J))**2 + (XC(2,I)-XC(2,J))**2
     $              + (XC(3,I)-XC(3,J))**2 )
      RD(J,I) = RD(I,J)
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
      X12 = XC(1,I)-XC(1,J)
      Y12 = XC(2,I)-XC(2,J)
      Z12 = XC(3,I)-XC(3,J)
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
C     (For 1,3-interactions B2(I,J) has been
C      set to zero in routine UFFB2)
C
        RV = RLEN(I,J)/R
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
      n=NB(M)
c
      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
c
      CosTh = (RIJ*RIJ+RJK*RJK-RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosTh).GT.ONE) CosTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)
cc      write(6,*) ' Bend',m,' is ',Th*degree
c
      If(n.EQ.0) Then
        EB = EB + B3(M)*(C0(M) + C1(M)*CosTh + C2(M)*DCOS(TWO*Th))
      Else If(n.EQ.1) Then
c -- linear
        EB = EB + B3(M)*(ONE + CosTh)
      Else
        EB = EB + (B3(M)/DFloat(n**2))*(ONE - DCOS(n*Th))
      EndIf
c
      SinTh = SQRT(ONE - CosTh**2)
C
C  ............................................................
C  If the three atoms are linear there are problems with the
C  angle bend derivative. Skip derivative evaluation
C
      If(n.NE.1.AND.SinTh.LT.TOLLZERO) GO TO 30
C  ............................................................
C
C  Set terms for bend derivatives
C  For the derivative, we consider the bending energy expression
C  as a function of CosTh -- dE/dX = dE/dCosTh x dCosTh/dX
C
      XIJ = XC(1,I)-XC(1,J)
      XJK = XC(1,J)-XC(1,K)
      YIJ = XC(2,I)-XC(2,J)
      YJK = XC(2,J)-XC(2,K)
      ZIJ = XC(3,I)-XC(3,J)
      ZJK = XC(3,J)-XC(3,K)
c
      D1 = TWO*RIJ*RJK
      D2 = RIJ*RIJ + RJK*RJK - RIK*RIK
c
      RIJJK = RIJ/RJK
      RJKIJ = RJK/RIJ
c
c -- set DCB
      If(n.EQ.0) Then
        DCB = B3(M)*(C1(M)+FOUR*C2(M)*CosTh)*TWO/(D1**2)
      Else If(n.EQ.4) Then
        DCB = -B3(M)*(TWO*CosTh**3 - CosTh)*TWO/(D1**2)
      Else If(n.EQ.3) Then
        DCB = -(B3(M)/THREE)*(FOUR*CosTh**2 - ONE)*TWO/(D1**2)
      Else If(n.EQ.1) Then
        DCB = B3(M)*TWO/(D1**2)
      Else
c -- should never get here!
        call nerror(1,'eguff','execution error',0,0)
      EndIf
c
      G(1,I) = G(1,I) - DCB*( D1*XJK + D2*XIJ*RJKIJ )
      G(1,J) = G(1,J) + DCB*( D1*(XJK-XIJ) -
     $                        D2*(XJK*RIJJK - XIJ*RJKIJ) )
      G(1,K) = G(1,K) + DCB*( D1*XIJ + D2*XJK*RIJJK )
      G(2,I) = G(2,I) - DCB*( D1*YJK + D2*YIJ*RJKIJ )
      G(2,J) = G(2,J) + DCB*( D1*(YJK-YIJ) -
     $                        D2*(YJK*RIJJK - YIJ*RJKIJ) )
      G(2,K) = G(2,K) + DCB*( D1*YIJ + D2*YJK*RIJJK )
      G(3,I) = G(3,I) - DCB*( D1*ZJK + D2*ZIJ*RJKIJ )
      G(3,J) = G(3,J) + DCB*( D1*(ZJK-ZIJ) -
     $                        D2*(ZJK*RIJJK - ZIJ*RJKIJ) )
      G(3,K) = G(3,K) + DCB*( D1*ZIJ + D2*ZJK*RIJJK )
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
      n=NT(M)
c
      RIJ=RD(I,J)
      RIK=RD(I,K)
      RJK=RD(J,K)
      RJL=RD(J,L)
      RKL=RD(K,L)
c
      XIJ =XC(1,I) - XC(1,J)
      XIK =XC(1,I) - XC(1,K)
      XJK =XC(1,J) - XC(1,K)
      XJL =XC(1,J) - XC(1,L)
      XKL =XC(1,K) - XC(1,L)
      YIJ =XC(2,I) - XC(2,J)
      YIK =XC(2,I) - XC(2,K)
      YJK =XC(2,J) - XC(2,K)
      YJL =XC(2,J) - XC(2,L)
      YKL =XC(2,K) - XC(2,L)
      ZIJ =XC(3,I) - XC(3,J)
      ZIK =XC(3,I) - XC(3,K)
      ZJK =XC(3,J) - XC(3,K)
      ZJL =XC(3,J) - XC(3,L)
      ZKL =XC(3,K) - XC(3,L)
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
cc      write(6,*) ' Torsion',m,' is ',Dih*degree
c
      ET = ET + B4(M)*(ONE - CosD(M)*DCOS(n*Dih))
C
C  Set terms for torsional derivatives
C  For the derivative, we consider the torsion energy expression
C  as a function of CosDih -- dE/dX = dE/dCosDih x dCosDih/dX
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
      If(n.EQ.3) THEN
       DCT = TWELVE*CosDih**2 - THREE
      Else If(n.EQ.6) THEN
       DCT = TWELVE*CosDih*
     $         ( 16.0d0*CosDih**2*(CosDih**2 - ONE) + THREE )
      Else If(n.EQ.2) THEN
       DCT = FOUR*CosDih
      Else
c -- should never get here!
        call nerror(1,'eguff','execution error',0,0)
      EndIf
      DCT = -B4(M)*CosD(M)*DCT/(A1*A1)
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
C  5. Inversions/Out-of-Plane bends
C
      DO 50 M=1,NumP
      I=IP1(M)
      J=IP2(M)
      K=IP3(M)
      L=IP4(M)
c
      XIJ = XC(1,I) - XC(1,J)
      XIK = XC(1,I) - XC(1,K)
      XIL = XC(1,I) - XC(1,L)
      YIJ = XC(2,I) - XC(2,J)
      YIK = XC(2,I) - XC(2,K)
      YIL = XC(2,I) - XC(2,L)
      ZIJ = XC(3,I) - XC(3,J)
      ZIK = XC(3,I) - XC(3,K)
      ZIL = XC(3,I) - XC(3,L)

      R1(1) = -XIJ
      R1(2) = -YIJ
      R1(3) = -ZIJ
      R2(1) = -XIK
      R2(2) = -YIK
      R2(3) = -ZIK
      R3(1) = -XIL
      R3(2) = -YIL
      R3(3) = -ZIL
      CALL VECMUL(R1,R2,R4)
c
      dn = SProd(3,R4,R3)
      dm1 = SQRT(SProd(3,R4,R4))
      dm2 = SQRT(SProd(3,R3,R3))
      dm = dm1*dm2
      CosTh = dn/dm
      If(Abs(CosTh).GT.ONE) CosTh=SIGN(ONE,CosTh)
      Th = ACOS(CosTh)
c -- this is the angle IL makes with the normal
c -- the angle we want is 90 degrees minus this
      Th = HalfPI-Th
cc      write(6,*) ' Out-of-plane bend ',m,' is ',Th*degree
C
C  ...........................................................
C     **  WARNING  **
C  If the inversion angle is close to zero there could be
C  problems with the derivative. Set the angle to zero and
C  skip derivative evaluation
C
      If(ABS(Th).LT.TOLLZERO) GO TO 50
C  ...........................................................
c
c -- in terms of CosTh, the energy expression becomes
c
      EP = EP + B5(M)*(S0(M) + S1(M)*SQRT(ONE - CosTh**2) -
     $                         S2(M)*(TWO*CosTh**2 - ONE) )
C
C  Set terms for out-of-plane bend derivatives
C
C  First set the derivatives of each term
C  1. The numerator
C      d/dx ( SProd(3,R4,R3) )
C
      DT1XI = -YIJ*ZIK + ZIJ*YIK + YIL*(ZIK-ZIJ) - ZIL*(YIK-YIJ)
      DT1XJ = -YIL*ZIK + ZIL*YIK
      DT1XK =  YIL*ZIJ - ZIL*YIJ
      DT1XL =  YIJ*ZIK - ZIJ*YIK
c
      DT1YI = -XIL*(ZIK-ZIJ) + XIJ*ZIK - ZIJ*XIK - ZIL*(XIJ-XIK)
      DT1YJ =  XIL*ZIK - ZIL*XIK
      DT1YK = -XIL*ZIJ + ZIL*XIJ
      DT1YL = -XIJ*ZIK + ZIJ*XIK
c
      DT1ZI = -XIL*(YIJ-YIK) + YIL*(XIJ-XIK) - XIJ*YIK + YIJ*XIK
      DT1ZJ = -XIL*YIK + YIL*XIK
      DT1ZK =  XIL*YIJ - YIL*XIJ
      DT1ZL =  XIJ*YIK - YIJ*XIK
C
C  2. The Denominator
C      d/dx ( SQRT(SProd(3,R4,R4)) )
C
      RA = (YIJ*ZIK-ZIJ*YIK)
      RB = (ZIJ*XIK-XIJ*ZIK)
      RC = (XIJ*YIK-YIJ*XIK)
c
      DT2XI = (ZIJ-ZIK)*RB + (YIK-YIJ)*RC
      DT2XJ =  ZIK*RB - YIK*RC
      DT2XK = -ZIJ*RB + YIJ*RC
c
      DT2YI = (ZIK-ZIJ)*RA + (XIJ-XIK)*RC
      DT2YJ = -ZIK*RA + XIK*RC
      DT2YK =  ZIJ*RA - XIJ*RC
c
      DT2ZI = (YIJ-YIK)*RA + (XIK-XIJ)*RB
      DT2ZJ =  YIK*RA - XIK*RB
      DT2ZK = -YIJ*RA + XIJ*RB
C
C  Now construct the final derivative
C
      DCP = -B5(M)*(S1(M)*CosTh/SQRT(ONE-CosTh**2) + FOUR*S2(M)*CosTh)
      DCP = DCP/dm**2
c
      dnm1 = dn*dm1/dm2
      dnm2 = dn*dm2/dm1
c
      G(1,I) = G(1,I) + DCP*(dm*DT1XI - dnm2*DT2XI - dnm1*XIL)
      G(1,J) = G(1,J) + DCP*(dm*DT1XJ - dnm2*DT2XJ)
      G(1,K) = G(1,K) + DCP*(dm*DT1XK - dnm2*DT2XK)
      G(1,L) = G(1,L) + DCP*(dm*DT1XL + dnm1*XIL)
      G(2,I) = G(2,I) + DCP*(dm*DT1YI - dnm2*DT2YI - dnm1*YIL)
      G(2,J) = G(2,J) + DCP*(dm*DT1YJ - dnm2*DT2YJ)
      G(2,K) = G(2,K) + DCP*(dm*DT1YK - dnm2*DT2YK)
      G(2,L) = G(2,L) + DCP*(dm*DT1YL + dnm1*YIL)
      G(3,I) = G(3,I) + DCP*(dm*DT1ZI - dnm2*DT2ZI - dnm1*ZIL)
      G(3,J) = G(3,J) + DCP*(dm*DT1ZJ - dnm2*DT2ZJ)
      G(3,K) = G(3,K) + DCP*(dm*DT1ZK - dnm2*DT2ZK)
      G(3,L) = G(3,L) + DCP*(dm*DT1ZL + dnm1*ZIL)
c
 50   CONTINUE
c
 95   CONTINUE
C
C  Form the total molecular energy
C
      ES = HALF*ES
c
      E = ES + EV + EB + ET + EP
c
      if(iprnt.gt.2) write(6,1000) E,ES,EB,EP,ET,EV
      RETURN
c
 1000 FORMAT(/,' TOTAL ENERGY       ',f15.8,/,
     $         ' Stretching         ',f15.8,/,
     $         ' Bending            ',f15.8,/,
     $         ' Out-of-Plane Bend  ',f15.8,/,
     $         ' Torsion            ',f15.8,/,
     $         ' VdW energy         ',f15.8)
c
      END
*Deck filluff
      SUBROUTINE FillUFF(NAtoms, NBend,  NTors,  NInv,   IPRNT,
     $                   IAT,    IC,     RLEN,   B2,     IB1,
     $                   IB2,    IB3,    C0,     C1,     C2,
     $                   B3,     NB,     IT1,    IT2,    IT3,
     $                   IT4,    CosD,   B4,     NT,     IP1,
     $                   IP2,    IP3,    IP4,    S0,     S1,
     $                   S2,     B5)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAT(NAtoms)
      dimension rlen(natoms,natoms),b2(natoms,natoms)
      dimension ib1(*),ib2(*),ib3(*),c0(*),c1(*),c2(*),b3(*),nb(*)
      dimension it1(*),it2(*),it3(*),it4(*),CosD(*),b4(*),nt(*)
      dimension ip1(*),ip2(*),ip3(*),ip4(*),s0(*),s1(*),s2(*),b5(*)
      common/numint/numb,numt,nump
c
c
c -- fill the bond and vdw arrays
      CALL UFFB2(NAtoms,IAT,IC,RLEN,B2)
c
c -- fill the angle bend arrays
      CALL UFFB3(NAtoms, NBend,  IAT,    IC,     RLEN,
     $           IB1,    IB2,    IB3,    C0,     C1,
     $           C2,     B3,     NB)
c
c -- fill the torsion arrays
      CALL UFFB4(NAtoms, NTors,  IAT,    IC,     IT1,
     $           IT2,    IT3,    IT4,    CosD,   B4,
     $           NT)
c
c -- fill the inversion arrays
      CALL UFFB5(NAtoms, NInv,   IAT,    IC,     IP1,
     $           IP2,    IP3,    IP4,    S0,     S1,
     $           S2,     B5)
c
      If(IPRNT.gt.3) THEN
       zero=0.0d0
       write(6,*) ' UFF Parameter Arrays are:'
       write(6,*) ' Distance Parameters'
       do 11 i=2,natoms
       do 11 j=1,i-1
       if(b2(i,j).ne.zero) write(6,1111) i,j,rlen(i,j),b2(i,j)
 1111  format(1x,2i4,2f12.6)
 11    continue
       if(numb.gt.0) then
        write(6,*) ' Angle Parameters'
        write(6,*)
     $'   I   J   K       C0          C1          C2          B3     NB'
        do 22 i=1,numb
        write(6,2222) ib1(i),ib2(i),ib3(i),c0(i),c1(i),c2(i),b3(i),nb(i)
 2222   format(1x,3i4,4f12.6,i4)
 22     continue
       endif
       if(numt.gt.0) then
        write(6,*) ' Torsion Parameters'
        write(6,*)
     $ '   I   J   K   L     CosD        B4       NT'
        do 33 i=1,numt
        write(6,3333) it1(i),it2(i),it3(i),it4(i),cosd(i),b4(i),nt(i)
 3333   format(1x,4i4,2f12.6,i4)
 33     continue
       endif
       if(nump.gt.0) then
        write(6,*) ' Inversion/Out-of-Plane Bend Parameters'
        write(6,*)
     $ '   I   J   K   L      S0          S1          S2          B5'
        do 44 i=1,nump
        write(6,4444)ip1(i),ip2(i),ip3(i),ip4(i),s0(i),s1(i),s2(i),b5(i)
 4444   format(1x,4i4,4f12.6)
 44     continue
       endif
      EndIf
      RETURN
      END
*Deck uffb2
      SUBROUTINE UFFB2(NAtoms,IAT,IC,RLEN,B2)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RLEN(NAtoms,NAtoms),B2(NAtoms,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      Dimension RR(127),          ! bonding atomic radii
     $          XX(127),          ! vdW atomic radii
     $          GMP(127),         ! GMP electronegativity
     $          DD(127)           ! atomic vdW energy
      COMMON /EFFECTIVE/ ZZ(127)  ! effective atomic charge
      COMMON /BONDORD/ BT(127)    ! atomic bond orders
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA RR /0.354d0,0.460d0,0.849d0,1.336d0,1.074d0,0.848d0,0.838d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         0.757d0,0.729d0,0.732d0,0.706d0,0.700d0,0.699d0,0.685d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         0.656d0,0.658d0,0.528d0,0.680d0,0.634d0,0.639d0,0.668d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $         0.920d0,1.539d0,1.421d0,1.244d0,1.117d0,1.101d0,1.056d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         1.056d0,1.064d0,1.049d0,1.027d0,1.077d0,0.854d0,1.044d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         1.032d0,1.953d0,1.761d0,1.513d0,1.412d0,1.412d0,1.402d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         1.345d0,1.382d0,1.270d0,1.335d0,1.241d0,1.164d0,1.302d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         1.193d0,1.260d0,1.197d0,1.211d0,1.190d0,1.192d0,1.147d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         2.260d0,2.052d0,1.698d0,1.564d0,1.473d0,1.467d0,1.484d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         1.322d0,1.478d0,1.332d0,1.338d0,1.386d0,1.403d0,1.459d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         1.398d0,1.407d0,1.386d0,1.382d0,1.267d0,2.570d0,2.277d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         1.943d0,1.841d0,1.823d0,1.816d0,1.801d0,1.780d0,1.771d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         1.735d0,1.732d0,1.710d0,1.696d0,1.673d0,1.660d0,1.637d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         1.671d0,1.611d0,1.511d0,1.392d0,1.526d0,1.380d0,1.372d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         1.314d0,1.372d0,1.371d0,1.364d0,1.262d0,1.340d0,1.518d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         1.459d0,1.512d0,1.50d0, 1.545d0,1.420d0,2.880d0,2.512d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         1.983d0,1.721d0,1.711d0,1.684d0,1.666d0,1.657d0,1.660d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         1.801d0,1.761d0,1.750d0,1.724d0,1.712d0,1.689d0,1.679d0,
c                Lw6+3
     $         1.698d0/
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA XX /2.886d0,2.886d0,2.362d0,2.451d0,2.745d0,4.083d0,4.083d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         3.851d0,3.851d0,3.851d0,3.851d0,3.660d0,3.660d0,3.660d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         3.660d0,3.500d0,3.500d0,3.500d0,3.500d0,3.500d0,3.364d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $         3.243d0,2.983d0,3.021d0,4.499d0,4.295d0,4.147d0,4.147d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         4.147d0,4.035d0,4.035d0,4.035d0,4.035d0,4.035d0,3.947d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         3.868d0,3.812d0,3.399d0,3.295d0,3.175d0,3.175d0,3.144d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         3.023d0,2.961d0,2.912d0,2.912d0,2.872d0,2.834d0,3.495d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         2.763d0,4.383d0,4.280d0,4.230d0,4.205d0,4.189d0,4.141d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         4.114d0,3.641d0,3.345d0,3.124d0,3.165d0,3.052d0,3.052d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         2.998d0,2.963d0,2.929d0,2.899d0,3.148d0,2.848d0,4.463d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         4.392d0,4.420d0,4.470d0,4.50d0, 4.404d0,4.517d0,3.703d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         3.522d0,3.556d0,3.606d0,3.575d0,3.547d0,3.520d0,3.493d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         3.368d0,3.451d0,3.428d0,3.409d0,3.391d0,3.374d0,3.355d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         3.640d0,3.141d0,3.170d0,3.069d0,3.069d0,3.069d0,2.954d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         2.954d0,3.120d0,2.840d0,2.754d0,3.293d0,2.705d0,4.347d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         4.297d0,4.370d0,4.709d0,4.750d0,4.765d0,4.90d0, 3.677d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         3.478d0,3.396d0,3.424d0,3.395d0,3.424d0,3.424d0,3.381d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         3.326d0,3.339d0,3.313d0,3.299d0,3.286d0,3.274d0,3.248d0,
c                Lw6+3
     $         3.236d0/
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA ZZ /0.712d0,0.712d0,0.098d0,1.026d0,1.565d0,1.755d0,1.755d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         1.912d0,1.912d0,1.912d0,1.912d0,2.544d0,2.544d0,2.544d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         2.544d0,2.300d0,2.300d0,2.300d0,2.300d0,2.300d0,1.735d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $         0.194d0,1.081d0,1.787d0,1.792d0,2.323d0,2.863d0,2.863d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         2.863d0,2.703d0,2.703d0,2.703d0,2.703d0,2.703d0,2.348d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         0.300d0,1.165d0,2.141d0,2.592d0,2.659d0,2.659d0,2.679d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         2.463d0,2.43d0, 2.43d0, 2.43d0, 2.43d0, 2.43d0, 1.756d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         1.308d0,1.821d0,2.789d0,2.864d0,2.764d0,2.519d0,0.452d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         1.592d0,2.449d0,3.257d0,3.667d0,3.618d0,3.40d0, 3.40d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         3.40d0, 3.40d0, 3.508d0,3.21d0, 1.956d0,1.65d0, 2.07d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         2.961d0,2.704d0,2.882d0,2.65d0, 0.556d0,1.573d0,2.727d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         3.30d0, 3.30d0, 3.30d0, 3.30d0, 3.30d0, 3.30d0, 3.30d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         3.30d0, 3.30d0, 3.30d0, 3.416d0,3.30d0, 3.30d0, 2.618d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         3.271d0,3.921d0,4.075d0,3.70d0, 3.70d0, 3.70d0, 3.70d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         3.70d0, 3.70d0, 3.731d0,3.382d0,2.625d0,1.75d0, 2.068d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         2.846d0,2.470d0,2.33d0, 2.24d0, 0.583d0,1.847d0,2.92d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         3.90d0, 4.202d0,3.90d0, 3.90d0, 3.90d0, 3.90d0, 3.90d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         3.90d0, 3.90d0, 3.90d0, 3.90d0, 3.90d0, 3.90d0, 3.90d0,
c                Lw6+3
     $         3.90d0 /
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA GMP/4.528d0,4.528d0,9.66d0, 3.006d0,4.877d0,5.11d0, 5.11d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         5.343d0,5.343d0,5.343d0,5.343d0,6.899d0,6.899d0,6.899d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         6.899d0,8.741d0,8.741d0,8.741d0,8.741d0,8.741d0,10.874d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $        11.04d0, 2.843d0,3.951d0,4.06d0, 4.168d0,5.463d0,5.463d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         5.463d0,6.928d0,6.928d0,6.928d0,6.928d0,6.928d0,8.564d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         9.465d0,2.421d0,3.231d0,3.395d0,3.47d0, 3.47d0, 3.65d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         3.415d0,3.325d0,3.76d0, 3.76d0, 4.105d0,4.465d0,3.729d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         5.106d0,3.641d0,4.051d0,5.188d0,6.428d0,7.790d0,8.505d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         2.331d0,3.024d0,3.83d0, 3.40d0, 3.55d0, 3.465d0,3.465d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         3.29d0, 3.575d0,3.975d0,4.32d0, 4.436d0,5.034d0,3.506d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         3.987d0,4.899d0,5.816d0,6.822d0,7.595d0,2.183d0,2.814d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         2.8355, 2.774d0,2.858d0,2.8685, 2.881d0,2.9115, 2.8785,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         3.1665, 3.018d0,3.0555, 3.127d0,3.1865, 3.2514, 3.2889,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         2.9629, 3.70d0, 5.10d0, 4.63d0, 4.63d0, 4.63d0, 3.96d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         3.96d0, 5.14d0, 5.00d0, 4.79d0, 4.894d0,6.27d0, 3.20d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         3.90d0, 4.69d0, 4.21d0, 4.75d0, 5.37d0, 2.00d0, 2.843d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         2.835d0,3.175d0,2.985d0,3.341d0,3.549d0,3.243d0,2.9895,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         2.8315, 3.1935, 3.197d0,3.333d0,3.40d0, 3.47d0, 3.475d0,
c                Lw6+3
     $         3.500d0/
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA DD /0.044d0,0.044d0,0.056d0,0.025d0,0.085d0,0.180d0,0.180d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         0.105d0,0.105d0,0.105d0,0.105d0,0.069d0,0.069d0,0.069d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         0.069d0,0.060d0,0.060d0,0.060d0,0.060d0,0.060d0,0.050d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $         0.042d0,0.030d0,0.111d0,0.505d0,0.402d0,0.305d0,0.305d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         0.305d0,0.274d0,0.274d0,0.274d0,0.274d0,0.274d0,0.227d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         0.185d0,0.035d0,0.238d0,0.019d0,0.017d0,0.017d0,0.016d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         0.015d0,0.013d0,0.013d0,0.013d0,0.014d0,0.015d0,0.005d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         0.124d0,0.415d0,0.379d0,0.309d0,0.291d0,0.251d0,0.220d0,
c                Rb      Sr6+2   Y.3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         0.04d0, 0.235d0,0.072d0,0.069d0,0.059d0,0.056d0,0.056d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         0.048d0,0.056d0,0.053d0,0.048d0,0.036d0,0.228d0,0.599d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         0.567d0,0.449d0,0.398d0,0.339d0,0.332d0,0.045d0,0.364d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         0.017d0,0.013d0,0.010d0,0.010d0,0.009d0,0.008d0,0.008d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         0.009d0,0.007d0,0.007d0,0.007d0,0.007d0,0.006d0,0.228d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         0.041d0,0.072d0,0.081d0,0.067d0,0.067d0,0.067d0,0.066d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         0.066d0,0.037d0,0.073d0,0.080d0,0.039d0,0.385d0,0.680d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         0.663d0,0.518d0,0.325d0,0.284d0,0.248d0,0.050d0,0.404d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         0.033d0,0.026d0,0.022d0,0.022d0,0.019d0,0.016d0,0.014d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         0.013d0,0.013d0,0.013d0,0.012d0,0.012d0,0.011d0,0.011d0,
c                Lw6+3
     $         0.011d0/
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA BT /1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  2.0d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         1.0d0,  1.5d0,  2.0d0,  3.0d0,  1.0d0,  1.5d0,  2.0d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         3.0d0,  1.0d0,  1.0d0,  1.5d0,  2.0d0,  3.0d0,  1.0d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.5d0,  2.0d0,  1.0d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  3.0d0,  1.0d0,  1.0d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  3.0d0,  1.0d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,  1.0d0,
c                Lw6+3
     $         1.0d0/
C
C
      CALL ZeroIT(RLEN,NAtoms*NAtoms)
      CALL ZeroIT(B2,NAtoms*NAtoms)
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
C  set up data for the bonded pair I1,J1
C
       RI = RR(I1)
       RJ = RR(J1)
       XI = GMP(I1)
       XJ = GMP(J1)
       REN = RI*RJ*(SQRT(XI) - SQRT(XJ))**2 / (XI*RI+XJ*RJ)
c
c -- set up bond order
c    assumed to be 1 unless determined otherwise
       BI = BT(I1)
       BJ = BT(J1)
       If(BI.EQ.BJ) Then
        bndrdr = BI
        If(IC(I,J).EQ.4) Then
         bndrdr = 1.41d0       ! special hack for amide bond
        EndIf
       Else
        bndrdr = 1.0d0
       EndIf
c
       RBO = 0.1332d0*(RI+RJ)*LOG(bndrdr)
C
C  now sum up
C
       RLEN(I,J) = RI + RJ - RBO - REN
       B2(I,J) = 664.12d0*(ZZ(I1)*ZZ(J1)/RLEN(I,J)**3)
cc
      ELSE
C
C  Van der Waals interaction
C  do not set for 1,3 bonds
C
       DO 20 K=1,NAtoms
       If(IC(I,K).NE.0.AND.IC(J,K).NE.0) GO TO 19
 20    CONTINUE
C
C set VdW interaction term
C
       B2(I,J) = SQRT(DD(I1)*DD(J1))
 19    CONTINUE
c -- always set distance term as may be needed for bend
       RLEN(I,J) = SQRT(XX(I1)*XX(J1))
      ENDIF
c
      RLEN(J,I) = RLEN(I,J)
      B2(J,I) = B2(I,J)
10    CONTINUE
C
      RETURN
      END
*Deck uffb3
      SUBROUTINE UFFB3(NAtoms, NBend,  IAT,    IC,     RLEN,
     $                 IB1,    IB2,    IB3,    C0,     C1,
     $                 C2,     B3,     NB)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RLEN(NAtoms,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IB1(*),IB2(*),IB3(*),C0(*),C1(*),C2(*),B3(*),NB(*)
      COMMON/NUMINT/NumB,NumT,NumP
C
C  This Subroutine calculates the number of bends. storing the
C  atoms involved in the arrays IB1, IB2 and IB3, and fills
C  the various parameter arrays, C0, C1, C2, B3 and NB
C
      NumB=0
      If(NAtoms.LT.3) RETURN
c
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
     $   CALL AMAUFF(NAtoms,I,J,K,II,JJ,KK,RLEN,IB1,IB2,IB3,
     $               C0,C1,C2,B3,NB)
      If(IC(I,K).NE.0.AND.IC(K,J).NE.0)
     $   CALL AMAUFF(NAtoms,I,K,J,II,KK,JJ,RLEN,IB1,IB2,IB3,
     $               C0,C1,C2,B3,NB)
      If(IC(J,I).NE.0.AND.IC(I,K).NE.0)
     $   CALL AMAUFF(NAtoms,J,I,K,JJ,II,KK,RLEN,IB1,IB2,IB3,
     $               C0,C1,C2,B3,NB)
c
c -- check number of bends does not exceed maximum allowed
      If(NumB.GT.NBend) Call nerror(2,'Force Field module',
     $           'More Bends in System than Dimensioned for',0,0)
c
 10   CONTINUE
c
      RETURN
      END
*Deck amauff
      SUBROUTINE AMAUFF(NAtoms,I,J,K,I1,J1,K1,RLEN,IB1,IB2,IB3,
     $                  C0,C1,C2,B3,NB)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 RLEN(NAtoms,NAtoms)
      DIMENSION IAT(NAtoms)
      DIMENSION IB1(*),IB2(*),IB3(*),C0(*),C1(*),C2(*),B3(*),NB(*)
      Dimension TH(127)           ! standard central-atom bond angles
      Dimension NN(127)           ! coordination parameter
c          n=1 linear; n=3 trigonal planar; n=4 square-planar/octahedral
      COMMON /EFFECTIVE/ ZZ(127)  ! effective atomic charge
      COMMON/NUMINT/NumB,NumT,NumP
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA TH /180.0d0, 83.5d0, 90.0d0,180.0d0,109.47, 109.47, 120.0d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         109.47, 120.0d0,120.0d0,180.0d0,106.7d0,120.0d0,111.2d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         180.0d0,104.51, 146.0d0,110.0d0,120.0d0,180.0d0,180.0d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $          90.0d0,180.0d0,109.47, 109.47, 109.47,  93.8d0,109.47,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         109.47,  92.1d0,103.2d0,109.47,  92.2d0,120.0d0,180.0d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $          90.0d0,180.0d0, 90.0d0,109.47, 109.47,  90.0d0,109.47,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $          90.0d0, 90.0d0,109.47,  90.0d0, 90.0d0, 90.0d0,109.47,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         109.47, 109.47, 109.47,  92.1d0, 90.6d0,180.0d0, 90.0d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         180.0d0, 90.0d0,109.47, 109.47, 109.47,  90.0d0,109.47,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $          90.0d0, 90.0d0, 90.0d0, 90.0d0,180.0d0,109.47, 109.47,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         109.47,  91.6d0, 90.23, 180.0d0, 90.0d0,180.0d0, 90.0d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         109.47,  90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $          90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $          90.0d0,109.47, 109.47,  90.0d0,109.47, 109.47,  90.0d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         109.47,  90.0d0, 90.0d0, 90.0d0, 90.0d0,180.0d0,120.0d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         109.47,  90.0d0, 90.0d0,180.0d0, 90.0d0,180.0d0, 90.0d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $          90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $          90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0, 90.0d0,
c                Lw6+3
     $          90.0d0/
c
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA NN /  0,      0,      4,      0,      0,      0,      3,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $           0,      3,      3,      1,      0,      3,      3,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $           1,      0,      0,      3,      3,      1,      0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $           4,      0,      0,      0,      0,      0,      0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $           0,      0,      0,      0,      3,      3,      0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $           4,      0,      4,      0,      0,      4,      0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $           4,      4,      0,      4,      4,      4,      0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $           0,      0,      0,      0,      0,      0,      4,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $           0,      4,      0,      0,      0,      4,      0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $           4,      4,      4,      4,      1,      0,      0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $           0,      0,      0,      0,      4,      0,      4,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $           0,      4,      4,      4,      4,      4,      4,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $           4,      4,      4,      4,      4,      4,      4,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $           4,      0,      0,      4,      0,      0,      4,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $           0,      4,      4,      4,      4,      1,      0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $           0,      0,      0,      0,      4,      0,      4,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $           4,      4,      4,      4,      4,      4,      4,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $           4,      4,      4,      4,      4,      4,      4,
c                Lw6+3
     $           4 /
c
      PARAMETER (One=1.0D0,Two=2.0d0,Three=3.0d0,Four=4.0d0)
      PARAMETER (DEGTORAD=3.14159265357979D0/180.0D0)
C
C
      NumB=NumB+1
      IB1(NumB)=I
      IB2(NumB)=J
      IB3(NumB)=K
C
C  Set the force constant and NB
C
      Theta = TH(J1)*DEGTORAD
      NB(NumB) = NN(J1)
      B3(NumB) = 664.12d0*(ZZ(I1)*ZZ(K1)/RLEN(I,K)**5)*
     $            (Three*RLEN(I,J)*RLEN(J,K)*(One - DCOS(Theta)**2) -
     $             RLEN(I,K)**2 * DCOS(Theta))
C
C Check the central atom
C
      If(NN(J1).EQ.0) Then
C
C  general non-linear bend
C
       C2(NumB) = One/(Four*DSIN(Theta)**2)
       C1(NumB) = -Four*C2(NumB)*DCOS(Theta)
       C0(NumB) = C2(NumB)*(Two*DCOS(Theta)**2 + One)
      EndIf
C
      RETURN
      END
*Deck uffb4
      SUBROUTINE UFFB4(NAtoms, NTors,  IAT,    IC,     IT1,
     $                  IT2,    IT3,    IT4,    CosD,   B4,
     $                  NT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),CosD(*),B4(*),NT(*)
      COMMON/NUMINT/NumB,NumT,NumP
c .........................................................
c   dynamically allocated temporary storage
      Dimension KIC(NAtoms,NAtoms)
c .........................................................
C
C  This Subroutine calculates the number of torsions, storing the
C  atoms involved in the arrays IT1, IT2, IT3 and IT4, and fills
C  the various parameter arrays, CosD, B4 and NB
C
      NumT=0
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
     $ CALL DMAUFF(NAtoms,K,I,J,L,K1,I1,J1,L1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(K,J).NE.0.AND.IC(I,L).NE.0)
     $ CALL DMAUFF(NAtoms,K,J,I,L,K1,J1,I1,L1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(I,K).NE.0.AND.IC(K,L).NE.0)
     $ CALL DMAUFF(NAtoms,J,I,K,L,J1,I1,K1,L1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(I,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMAUFF(NAtoms,J,I,L,K,J1,I1,L1,K1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(J,K).NE.0.AND.IC(K,L).NE.0)
     $ CALL DMAUFF(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(J,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMAUFF(NAtoms,I,J,L,K,I1,J1,L1,K1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      ENDIF
c
      IF(IC(J,K).NE.0) THEN
C
C  J-K connected
C
      If(IC(K,I).NE.0.AND.IC(I,L).NE.0)
     $ CALL DMAUFF(NAtoms,J,K,I,L,J1,K1,I1,L1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(I,K).NE.0.AND.IC(J,L).NE.0)
     $ CALL DMAUFF(NAtoms,I,K,J,L,I1,K1,J1,L1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(I,L).NE.0.AND.IC(L,J).NE.0)
     $ CALL DMAUFF(NAtoms,I,L,J,K,I1,L1,J1,K1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(I,L).NE.0.AND.IC(L,K).NE.0)
     $ CALL DMAUFF(NAtoms,I,L,K,J,I1,L1,K1,J1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      ENDIF
c
      IF(IC(I,K).NE.0) THEN
C
C  I-K connected
C
      If(IC(J,L).NE.0.AND.IC(L,I).NE.0)
     $ CALL DMAUFF(NAtoms,J,L,I,K,J1,L1,I1,K1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      If(IC(K,L).NE.0.AND.IC(L,J).NE.0)
     $ CALL DMAUFF(NAtoms,I,K,L,J,I1,K1,L1,J1,IC,IT1,IT2,IT3,IT4,
     $             CosD,B4,NT)
      ENDIF
c
c -- check number of torsions does not exceed maximum allowed
      If(NumT.GT.NTors) Call nerror(3,'Force Field module',
     $           'More Torsions in System than Dimensioned for',0,0)
c
 10   CONTINUE
C
C  torsion force constants have to be divided by the number of torsions
C  around a given central bond
C
      CALL IZeroIT(KIC,NAtoms**2)
      DO 20 M=1,NumT
C
C  count the number for each J-K pair
C
      J = IT2(M)
      K = IT3(M)
      KIC(J,K) = KIC(J,K)+1
      KIC(K,J) = KIC(J,K)
 20   CONTINUE
c
      DO 30 M=1,NumT
C
C  divide the force constant by the number of torsions
C
      J = IT2(M)
      K = IT3(M)
      B4(M) = B4(M)/DBLE(KIC(J,K))
 30   CONTINUE
C
      RETURN
      END
*Deck dmauff
      SUBROUTINE DMAUFF(NAtoms,I,J,K,L,I1,J1,K1,L1,IC,IT1,IT2,IT3,IT4,
     $                  CosD,B4,NT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IC(NAtoms,NAtoms)
      DIMENSION IT1(*),IT2(*),IT3(*),IT4(*),CosD(*),B4(*),NT(*)
      Dimension VV(127)           ! atomic torsion barriers
      COMMON/BONDORD/BT(127)
      COMMON/NUMINT/NumB,NumT,NumP
c
c -- ** negative values imply sp2 hybridization **
c                H_      H_b     He4+4   Li      Be3+2   B_3     B_2
      DATA VV /0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0, -2.0d0,
c                C_3     C_R     C_2     C_1     N_3     N_R     N_2
     $         2.119d0,-2.0d0,-2.0d0,  0.0d0,  0.450d0,-2.0d0,-2.0d0,
c                N_1     O_3     O_3_z   O_R     O_2     O_1     F_
     $         0.0d0,  0.018d0,0.018d0,-2.0d0,-2.0d0,  0.0d0,  0.0d0,
c                Ne4+4   Na      Mg3+2   Al3     Si3     P_3+3   P_3+5
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.225d0,2.400d0,2.400d0,
c                P_3+q   S_3+2   S_3+4   S_3+6   S_R     S_2     Cl
     $         2.400d0,0.484d0,0.484d0,0.484d0,-1.25d0,-1.25d0,  0.0d0,
c                Ar4+4   K_      Ca6+2   Sc3+3   Ti3+4   Ti6+4   V_3+5
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Cr6+3   Mn6+2   Fe3+2   Fe6+2   Co6+3   Ni4+2   Cu3+1
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Zn3+2   Ga3+3   Ge3     As3+3   Se3+2   Br      Kr4+4
     $         0.0d0,  0.0d0,  0.701d0,1.5d0,  0.335d0,0.0d0,  0.0d0,
c                Rb      Sr6+2   Y_3+3   Zr3+4   Nb3+5   Mo6+6   Mo3+6
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Tc6+5   Ru6+2   Rh6+3   Pd4+2   Ag1+1   Cd3+2   In3+3
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Sn3     Sb3+3   Te3+2   I_      Xe4+4   Cs      Ba6+2
     $         0.199d0,1.1d0,  0.3d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                La3+3   Ce6+3   Pr6+3   Nd6+3   Pm6+3   Sm6+3   Eu6+3
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Gd6+3   Tb6+3   Dy6+3   Ho6+3   Er6+3   Tm6+3   Yb6+3
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Lu6+3   Hf3+4   Ta3+5   W_6+6   W_3+4   W_3+6   Re6+5
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Re3+7   Os6+6   Ir6+3   Pt4+2   Au4+3   Hg1+2   Tl3+3
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Pb3     Bi3+3   Po3+2   At      Rn4+4   Fr      Ra6+2
     $         0.1d0,  1.0d0,  0.3d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Ac6+3   Th6+4   Pa6+4   U_6+4   Np6+4   Pu6+4   Am6+4
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Cm6+3   Bk6+3   Cf6+3   Es6+3   Fm6+3   Md6+3   No6+3
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
c                Lw6+3
     $         0.0d0/
c
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
C
C
C  For a torsion to be acceptable, there must be values
C  for BOTH central atoms
C
      If(VV(J1).EQ.Zero.OR.VV(K1).EQ.Zero) RETURN
c
      NumT=NumT+1
      IT1(NumT)=I
      IT2(NumT)=J
      IT3(NumT)=K
      IT4(NumT)=L
C
C  Check hybridization
C
      IF(VV(J1).GT.Zero.AND.VV(K1).GT.Zero) THEN
C
C (a) sp3-sp3
C
       B4(NumT) = Half*SQRT(VV(J1)*VV(K1))
       CosD(NumT) = -One
       NT(NumT) = 3
c
      ELSE IF(VV(J1)*VV(K1).LT.Zero) THEN
C
C (b) sp3-sp2
C
       B4(NumT) = Half
       CosD(NumT) = One
       NT(NumT) = 6
c
      ELSE
C
C (c) sp2-sp2
C
c -- set up bond order
c    assumed to be 1 unless determined otherwise
       BJ = BT(J1)
       BK = BT(K1)
       If(BJ.EQ.BK) Then
        bndrdr = BJ
        If(IC(J,K).EQ.4) Then
         bndrdr = 1.41d0       ! special hack for amide bond
        EndIf
       Else
        bndrdr = One
       EndIf
c
       B4(NumT) = 2.5d0*SQRT(VV(J1)*VV(K1))*(One + 4.18d0*LOG(bndrdr))
       CosD(NumT) = One
       NT(NumT) = 2
      ENDIF
C
      RETURN
      END
*Deck uffb5
      SUBROUTINE UFFB5(NAtoms, NInv,   IAT,    IC,     IP1,
     $                 IP2,    IP3,    IP4,    S0,     S1,
     $                 S2,     B5)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION IP1(*),IP2(*),IP3(*),IP4(*),S0(*),S1(*),S2(*),B5(*)
      Dimension IOP(6)
      Logical carbon
      COMMON/NUMINT/NumB,NumT,NumP
c
      PARAMETER (DEGTORAD=3.14159265357979D0/180.0D0)
C
C  This Subroutine calculates the number of out-of-plane bends, storing the
C  atoms involved in the arrays IT1, IT2, IT3 and IT4, and fills the various
C  parameter arrays, S0, S1, S2 and B5
C
      NumP=0
      If(NAtoms.LT.4) RETURN
c
      DO 20 I=1,NAtoms
      carbon = .FALSE.
      I1 = IAT(I)
C
C  go through each atom and see it forms an out-of-plane bend
C  only certain atoms are allowed an inversion term and they
C  must be bonded to exactly three other atoms
C  these are: C_2, C_R, N_3, N_2, N_R, P_3, As3, Sb3 and Bi3
C
      IF(I1.EQ.9.OR.I1.EQ.10) THEN
c -- carbon
       C2 =  0.0d0
       C1 = -1.0d0
       C0 =  1.0d0
       FC =  6.0d0/3.0d0
       carbon = .TRUE.
      ELSE IF(I1.EQ.12.OR.I1.EQ.13.OR.I1.EQ.14) THEN
c -- nitrogen
       C2 =  0.0d0
       C1 = -1.0d0
       C0 =  1.0d0
       FC =  6.0d0/3.0d0
      ELSE IF(I1.EQ.27.OR.I1.EQ.28.OR.I1.EQ.29) THEN
c -- phosphorus
       phi = 84.4339d0*DEGTORAD      ! position of minimum energy
       C2 = 1.0d0                    ! arbitrary
       C1 = -4.0d0*DCOS(phi)
       C0 = -C1*DCOS(phi) + C2*DCOS(2.0d0*phi)
       FC = (22.0d0/3.0d0)/(C0+C1+C2)
      ELSE IF(I1.EQ.53) THEN
c -- arsenic
       phi = 86.9735d0*DEGTORAD      ! position of minimum energy
       C2 = 1.0d0                    ! arbitrary
       C1 = -4.0d0*DCOS(phi)
       C0 = -C1*DCOS(phi) + C2*DCOS(2.0d0*phi)
       FC = (22.0d0/3.0d0)/(C0+C1+C2)
      ELSE IF(I1.EQ.72) THEN
c -- antimony
       phi = 87.7047d0*DEGTORAD      ! position of minimum energy
       C2 = 1.0d0                    ! arbitrary
       C1 = -4.0d0*DCOS(phi)
       C0 = -C1*DCOS(phi) + C2*DCOS(2.0d0*phi)
       FC = (22.0d0/3.0d0)/(C0+C1+C2)
      ELSE IF(I1.EQ.107) THEN
c -- bismuth
       phi = 90.0d0*DEGTORAD         ! position of minimum energy
       C2 = 1.0d0                    ! arbitrary
       C1 = -4.0d0*DCOS(phi)
       C0 = -C1*DCOS(phi) + C2*DCOS(2.0d0*phi)
       FC = (22.0d0/3.0d0)/(C0+C1+C2)
      ELSE
c -- atom has no inversion defined
       GO TO 20
      ENDIF

C  At this point the atom could potentially have an inversion
C  Check that there are three connected atoms
C
      nc = 0
      DO 10 J=1,NAtoms
      If(IC(I,J).NE.0) Then
       nc = nc+1
       IOP(nc) = J
      EndIf
 10   CONTINUE
c
      IF(nc.EQ.3) THEN
C
C  there is an inversion
C
       NumP = NumP+1
       IP1(NumP) = I
       IP2(NumP) = IOP(1)
       IP3(NumP) = IOP(2)
       IP4(NumP) = IOP(3)
C
C  If the central atom is carbon and any one of the other atoms
C  is oxygen (O_2), reset the force constant to 50 kcal/mol
C
       If(carbon) Then
        If(IAT(IOP(1)).EQ.19.OR.IAT(IOP(2)).EQ.19.OR.
     $     IAT(IOP(3)).EQ.19) FC = 50.0d0/3.0d0
       EndIf
c
       S0(NumP) = C0
       S1(NumP) = C1
       S2(NumP) = C2
       B5(NumP) = FC
C
C  now cycle through all 3 atoms
C
       NumP = NumP+1
       IP1(NumP) = I
       IP2(NumP) = IOP(2)
       IP3(NumP) = IOP(3)
       IP4(NumP) = IOP(1)
c
       S0(NumP) = C0
       S1(NumP) = C1
       S2(NumP) = C2
       B5(NumP) = FC
c
       NumP = NumP+1
       IP1(NumP) = I
       IP2(NumP) = IOP(3)
       IP3(NumP) = IOP(1)
       IP4(NumP) = IOP(2)
c
       S0(NumP) = C0
       S1(NumP) = C1
       S2(NumP) = C2
       B5(NumP) = FC
      ENDIF
c
c -- check number of inversions does not exceed maximum allowed
      If(NumP.GT.NInv) Call nerror(4,'Force Field module',
     $  'More Out-of-Plane Bends in System than Dimensioned for',0,0)
c
 20   CONTINUE
C
      RETURN
      END
*Deck getnprimf
      SUBROUTINE GetNPRIMF(NAtoms,IC,NBend,NTors,NInv)
      IMPLICIT INTEGER(A-Z)
C
C  Estimates maximum number of bends and torsions (including
C  out-of-plane bends) in the current molecule based on the
C  atomic connectivity
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IC      -  connectivity matrix
C  NBend   -  estimate of maximum number of bends
C  NTors   -  estimate of maximum number of torsions/oop bends
C
      DIMENSION IC(NAtoms,NAtoms)
C
C
C -- BENDS --
C
      NBend = 0
c
      DO 10 I=3,NAtoms
      DO 10 J=2,I-1
      DO 10 K=1,J-1
c
      If(IC(I,J).NE.0.AND.IC(J,K).NE.0) NBend = NBend+1
      If(IC(J,K).NE.0.AND.IC(K,I).NE.0) NBend = NBend+1
      If(IC(K,I).NE.0.AND.IC(I,J).NE.0) NBend = NBend+1
c
 10   CONTINUE
C
C -- OUT-OF-PLANE BENDS (INVERTIONS)
C
      NInv = 0
c
      DO 20 I=4,NAtoms
      DO 20 J=3,I-1
      DO 20 K=2,J-1
      DO 20 L=1,K-1
c
      If(IC(I,J).NE.0.AND.IC(I,K).NE.0.AND.IC(I,L).NE.0) NInv = NInv+3
      If(IC(J,I).NE.0.AND.IC(J,K).NE.0.AND.IC(J,L).NE.0) NInv = NInv+3
      If(IC(K,I).NE.0.AND.IC(K,J).NE.0.AND.IC(K,L).NE.0) NInv = NInv+3
      If(IC(L,I).NE.0.AND.IC(L,J).NE.0.AND.IC(L,K).NE.0) NInv = NInv+3
c
 20   CONTINUE
C
C -- TORSIONS --
C
      NTors = 0
c
      DO 32 I=2,NAtoms
      DO 32 J=1,I-1
      IF(IC(I,J).NE.0) THEN
C
C  consider I-J as middle 2 atoms in proper torsion
C
       DO 31 K=1,NAtoms
       IF(IC(I,K).NE.0.AND.K.NE.J) THEN
C
C  have K-I-J
C
        DO 30 L=1,NAtoms
        If(IC(J,L).NE.0.AND.I.NE.L.AND.K.NE.L) NTors = NTors+1
 30     CONTINUE
       ENDIF
 31    CONTINUE
      ENDIF
c
 32   CONTINUE
C
      RETURN
      END
*Deck uff-preamble
      SUBROUTINE UFF_PREAMBLE(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $                        IC,     RD,     RLEN,   B2,     IB1,
     $                        IB2,    IB3,    C0,     C1,     C2,
     $                        B3,     NB,     IT1,    IT2,    IT3,
     $                        IT4,    CosD,   B4,     NT,     IP1,
     $                        IP2,    IP3,    IP4,    S0,     S1,
     $                        S2,     B5,     G,      CT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XC(3*NAtoms),RD(NAtoms,NAtoms),RLEN(NAtoms,NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION B2(NAtoms,NAtoms),IB1(*),IB2(*),IB3(*),C0(*),C1(*),
     $          C2(*),B3(*),NB(*),IT1(*),IT2(*),IT3(*),IT4(*),
     $          CosD(*),B4(*),NT(*),IP1(*),IP2(*),IP3(*),IP4(*),
     $          s0(*),S1(*),S2(*),B5(*)
      REAL*8 G(3*NAtoms),CT(3*NAtoms)
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
      CALL EGUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $           IC,     RD,     RLEN,   B2,     IB1,
     $           IB2,    IB3,    C0,     C1,     C2,
     $           B3,     NB,     IT1,    IT2,    IT3,
     $           IT4,    CosD,   B4,     NT,     IP1,
     $           IP2,    IP3,    IP4,    S0,     S1,
     $           S2,     B5,     E0,     G)
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
      CALL CpyVEC(NAT3,XC,CT)
C
C  take the step
C
      DO 20 I=1,NAT3
      XC(I) = XC(I) - G(I)
 20   CONTINUE
C
C  calculate a new energy
C
      CALL EUFF(NAtoms, CutOff, IPRNT,  XC,     IAT,
     $          IC,     RD,     RLEN,   B2,     IB1,
     $          IB2,    IB3,    C0,     C1,     C2,
     $          B3,     NB,     IT1,    IT2,    IT3,
     $          IT4,    CosD,   B4,     NT,     IP1,
     $          IP2,    IP3,    IP4,    S0,     S1,
     $          S2,     B5,     E1)
C
C  if the energy goes down, take another steepest descent step
C  otherwise revert to previous geometry and start Hessian
C  optimization
C
      IF(E1.LT.E0) THEN
       E0 = E1
       GO TO 100
      ELSE
       CALL CpyVEC(NAT3,CT,XC)
      ENDIF
c
      If(IPRNT.GT.0) WRITE(6,*) ' End of Steepest Descent'
c
      RETURN
      END
*Deck connectf
      SUBROUTINE UFF_CONNECTF(NAtoms, IPRNT,  IAN,    XC,     RD,
     $                        NBOND,  IVAL,   IAT,    IC,     IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL ERROR
C
      REAL*8 XC(3,NAtoms),RD(NAtoms,NAtoms)
      DIMENSION IAN(NAtoms),IAT(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION NBOND(NAtoms),IVAL(NAtoms)
      DIMENSION IUFF(103)
      Character*256 jobname
c
      Common /job/jobname,lenJ
      Parameter (IUnit=1)
C
C  IUFF - Likely UFF atom type for each atom
C
C                    H He Li Be  B  C  N  O F  Ne
      DATA IUFF/     1, 3, 4, 5, 6, 8,12,16,21,22,
C
C                         Na Mg Al Si  P  S Cl Ar
     *                    23,24,25,26,27,30,35,36,
C
C                          K Ca
     *                    37,38,
C
C                   Sc Ti  V Cr Mn Fe Co Ni Cu Zn
     *              39,40,42,43,44,45,47,48,49,50,
C
C                               Ga Ge As Se Br Kr
     *                          51,52,53,54,55,56,
C
C                         Rb Sr
     *                    57,58,
C
C                    Y Zr Nb Mo To Ru Rh Pd Ag Cd
     *              59,60,61,62,64,65,66,67,68,69,
C
C                               In Sn Sb Te  I Xe
     *                          70,71,72,73,74,75,
C
C                         Cs Ba
     *                    76,77,
C
C                   La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
     *              78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
C
C                   Hf Ta  W Re  Os  Ir  Pt  Au  Hg
     *              93,94,95,98,100,101,102,103,104,
C
C                               Tl  Pb  Bi  Po  At  Rn
     *                         105,106,107,108,109,110,
C
C                          Fr Ra
     *                    111,112,
C
C                    Ac  Th  Pa   U  Np  Pu  Am  Cm  Bk  Cf  Es  Fm
     *              113,114,115,116,117,118,119,120,121,122,123,123,
C
C                    Md  No  Lr
     *              125,126,127 /
C
C
C  Assign/read in the atom types and connectivity for the UFF force field
C  If there is any user-defined bonding data to be read in, it should be
C  in the file <jobname.fld>
C
C  The standard UFF atom types and their numerical values
C  are listed below
C
C    1     H_      hydrogen (linear)
C    2     H_b     hydrogen bridging (e.g., in B2H6)
C    3     He4     helium square-planar
C    4     Li      lithium (linear)
C    5     Be3     beryllium sp3
C    6     B_3     boron sp3
C    7     B_2     boron sp2
C    8     C_3     carbon sp3
C    9     C_R     carbon aromatic
C   10     C_2     carbon sp2
C   11     C_1     carbon linear
C   12     N_3     nitrogen sp3
C   13     N_R     nitrogen aromatic
C   14     N_2     nitrogen sp2
C   15     N_1     nitrogen linear
C   16     O_3     oxygen sp3
C   17     O_3_z   oxygen sp3 (in zeolites)
C   18     O_R     oxygen aromatic
C   19     O_2     oxygen sp2
C   20     O_1     oxygen linear
C   21     F       fluorine (linear)
C   22     Ne4     neon square-planar
C   23     Na      sodium (linear)
C   24     Mg3     magnesium sp3
C   25     Al3     aluminium sp3
C   26     Si3     silicon sp3
C   27     P_3+3   phosphorus(3) tetrahedral
C   28     P_3+5   phosphorus(5) tetrahedral
C   29     P_3+q   phosphorus(3) tetrahedral, 4-coordinate
C   30     S_3+2   sulphur(2) tetrahedral
C   31     S_3+4   sulphur(4) tetrahedral
C   32     S_3+6   sulphur(6) tetrahedral
C   33     S_R     sulphur aromatic
C   34     S_2     sulphur trigonal planar
C   35     Cl      chlorine (linear)
C   36     Ar4     argon square-planar
C   37     K_      potassium (linear)
C   38     Ca6     calcium octahedral
C   39     Sc3     scandium(3) tetrahedral
C   40     Ti3     titanium(4) tetrahedral
C   41     Ti6     titanium(4) octahedral
C   42     V_3     vanadium(5) tetrahedral
C   43     Cr6     chromium(3) octahedral
C   44     Mn6     manganese(2) octahedral
C   45     Fe3     iron(2) tetrahedral
C   46     Fe6     iron(2) octahedral
C   47     Co6     cobalt(3) octahedral
C   48     Ni4     nickel(2) square-planar
C   49     Cu3     copper(1) tetrahedral
C   50     Zn3     zinc(2) tetrahedral
C   51     Ga3     gallium(3) tetrahedral
C   52     Ge3     germanium tetrahedral
C   53     As3     arsenic(3) tetrahedral
C   54     Se3     selenium(2) tetrahedral
C   55     Br      bromine (linear)
C   56     Kr4     krypton square-planar
C   57     Rb      rubidium (linear)
C   58     Sr6     strontium octahedral
C   59     Y_3     yttrium(3) tetrahedral
C   60     Zr3     zirconium(3) tetrahedral
C   61     Nb3     niobium(5) tetrahedral
C   62     Mo6     molybdenum(6) octahedral
C   63     Mo3     molybdenum(6) tetrahedral
C   64     Tc6     technetium(5) octahedral
C   65     Ru6     ruthenium(2) octahedral
C   66     Rh6     rhodium(3) octahedral
C   67     Pd4     palladium(2) square-planar
C   68     Ag1     silver(1) linear
C   69     Cd3     cadmium(2) tetrahedral
C   70     In3     indium(3) tetrahedral
C   71     Sn3     tin tetrahedral
C   72     Sb3     antimony(3) tetrahedral
C   73     Te3     tellurium(2) tetrahedral
C   74     I       iodine (linear)
C   75     Xe4     xenon square-planar
C   76     Cs      caesium (linear)
C   77     Ba6     barium octahedral
C   78     La3     lanthanum(3) tetrahedral
C   79     Ce6     cerium(3) octahedral
C   80     Pr6     praseodymium(3) octahedral
C   81     Nd6     neodymium(3) octahedral
C   82     Pm6     promethium(3) octahedral
C   83     Sm6     samarium(3) octahedral
C   84     Eu6     europium(3) octahedral
C   85     Gd6     gadolinium(3) octahedral
C   86     Tb6     terbium(3) octahedral
C   87     Dy6     dysprosium(3) octahedral
C   88     Ho6     holmium(3) octahedral
C   89     Er6     erbium(3) octahedral
C   90     Tm6     thulium(3) octahedral
C   91     Yb6     ytterbium(3) octahedral
C   92     Lu6     lutetium(3) octahedral
C   93     Hf3     hafnium(4) tetrahedral
C   94     Ta3     tantalum(5) tetrahedral
C   95     W_6     tungsten(6) octahedral
C   96     W_3     tungsten(4) tetrahedral
C   97     W_3     tungsten(6) tetrahedral
C   98     Re6     rhenium(5) octahedral
C   99     Re3     rhenium(7) tetrahedral
C  100     Os6     osmium(6) octahedral
C  101     Ir6     iridium(3) octahedral
C  102     Pt4     platinum(2) square-planar
C  103     Au4     gold(3) square-planar
C  104     Hg1     mercury(2) linear
C  105     Tl3     thallium(3) tetrahedral
C  106     Pb3     lead tetrahedral
C  107     Bi3     bismuth(3) tetrahedral
C  108     Po3     polonium(2) tetrahedral
C  109     At      astatine (linear)
C  110     Rn4     radon square-planar
C  111     Fr      francium (linear)
C  112     Ra6     radium(2) octahedral
C  113     Ac6     actinium(3) octahedral
C  114     Th6     thorium(4) octahedral
C  115     Pa6     protactinium(4) octahedral
C  116     U_6     uranium(4) octahedral
C  117     Np6     neptunium(4) octahedral
C  118     Pu6     plutonium(4) octahedral
C  119     Am6     americium(4) octahedral
C  120     Cm6     curium(3) octahedral
C  121     Bk6     berkelium(3) octahedral
C  122     Cf6     californium(3) octahedral
C  123     Es6     einsteinium(3) octahedral
C  124     Fm6     fermium(3) octahedral
C  125     Md6     mendelevium(3) octahedral
C  126     No6     nobelium(3) octahedral
C  127     Lr6     lawrencium(3) octahedral
C
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
      CALL BONDF(NAtoms, IPRNT,  IAN,    XC,     RD,
     $           NBOND,  IVAL,   IAT,    IC,     IStatus)
C
C  Is there any user-defined bonding data?
C  This will overwrite bonding produced by <BONDF>
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
C     atom I    atom J      bond type
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
C  determine UFF atom types
C
      DO 20 IA=1,NAtoms
      IAA = IAN(IA)
C
C  check atoms that can have multiple UFF types
C
      IF(IAA.EQ.1) THEN
c -- hydrogen
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(NIC.EQ.1) Then
         IAT(IA) = 1
       Else
         IAT(IA) = 2     ! bridging H
       EndIf
      ELSE IF(IAA.EQ.5) THEN
c -- boron
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(NIC.EQ.3) Then
         IAT(IA) = 7
       Else
         IAT(IA) = 6
       EndIf
      ELSE IF(IAA.EQ.6) THEN
c -- carbon
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.EQ.1) Then
         IAT(IA) = 8
       Else If(LIC.EQ.2) Then
         IAT(IA) = 10
       Else If(LIC.EQ.3) Then
         IAT(IA) = 11
       Else
         IAT(IA) = 9     ! aromatic
       EndIf
      ELSE IF(IAA.EQ.7) THEN
c -- nitrogen
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.LE.1) Then
         IAT(IA) = 12
       Else If(LIC.EQ.2) Then
         IAT(IA) = 14
       Else If(LIC.EQ.3) Then
         IAT(IA) = 15
       Else
         IAT(IA) = 13    ! aromatic
       EndIf
      ELSE IF(IAA.EQ.8) THEN
c -- oxygen
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.LE.1) Then
         IAT(IA) = 16
cc ** need to check for zeolite **
       Else If(LIC.EQ.2) Then
         IAT(IA) = 19
       Else If(LIC.EQ.3) Then
         IAT(IA) = 20
       Else
         IAT(IA) = 18    ! aromatic
       EndIf
      ELSE IF(IAA.EQ.15) THEN
c -- phosphorus
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(NIC.EQ.4) Then
         IAT(IA) = 29
       Else If(MIC.GT.3) Then
         IAT(IA) = 28
       Else
         IAT(IA) = 27
       EndIf
      ELSE IF(IAA.EQ.16) THEN
c -- sulphur
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(LIC.EQ.2.AND.NIC.EQ.1) Then
         IAT(IA) = 34
       Else If(LIC.EQ.5) Then
         IAT(IA) = 33    ! aromatic
       Else If(MIC.EQ.6) Then
         IAT(IA) = 32
       Else If(MIC.EQ.4) Then
         IAT(IA) = 31
       Else
         IAT(IA) = 30
       EndIf
      ELSE IF(IAA.EQ.22) THEN
c -- titanium
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(MIC.GT.4) Then
         IAT(IA) = 41
       Else
         IAT(IA) = 40
       EndIf
      ELSE IF(IAA.EQ.26) THEN
c -- iron
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(MIC.GT.4) Then
         IAT(IA) = 46
       Else
         IAT(IA) = 45
       EndIf
      ELSE IF(IAA.EQ.42) THEN
c -- molybdenum
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(MIC.GT.4) Then
         IAT(IA) = 63
       Else
         IAT(IA) = 62
       EndIf
      ELSE IF(IAA.EQ.74) THEN
c -- tungsten
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(NIC.GT.4) Then
         IAT(IA) = 95
       Else If(MIC.EQ.4) Then
         IAT(IA) = 96
       Else
         IAT(IA) = 97
       EndIf
      ELSE IF(IAA.EQ.75) THEN
c -- rhenium
       CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
       If(MIC.GT.4) Then
         IAT(IA) = 98
       Else
         IAT(IA) = 99
       EndIf
      ELSE
c -- may not be correct, but currently only one option
        IAT(IA) = IUFF(IAA)
      ENDIF
 20   CONTINUE
C
C  Print Out atom types and connectivities selected
C  and check for errors
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
c -- set all unknown bond types to single bonds
       If(IC(IA,JA).LT.0) Then
        If(IPRNT.GT.0) WRITE(6,1450) IA,JA
        IC(IA,JA) = 1
        IC(JA,IA) = 1
       EndIf
       If(IPRNT.GT.2) WRITE(6,1400) IA,JA,IC(IA,JA)
cc       If(IC(IA,JA).LT.0) ERROR=.TRUE.
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
 1000 FORMAT(' UFF Bonding Data Read In')
 1100 FORMAT(/,2X,' UFF Atom types',//,
     $         3X,' ATOM  TYPE',/,
     $         3X,' ----  ----')
 1200 FORMAT (4X,I3,2X,I4)
 1300 FORMAT (/,2X,'BONDS info',//,
     $          3X,' Atom1  Atom2  Bond Type',/,
     $          3X,' -----  -----  ---------')
 1400 FORMAT(4X,I4,4X,I4,4X,I4)
 1450 FORMAT('**WARNING** single bond assumed between atoms: ',2I4)
 1500 FORMAT(/,2X,'***ERROR*** Unable to determine connectivities.',
     $            ' Data must be read in')
c
      END
*Deck bondf
      SUBROUTINE BONDF(NAtoms, IPRNT,  IAN,    XC,     RD,
     $                 NBOND,  IVAL,   IAT,    IC,     IStatus)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL ERROR,DONE,amide
C
      DIMENSION IAN(NAtoms),XC(3,NAtoms),RD(NAtoms,NAtoms)
      DIMENSION NBOND(NAtoms),IVAL(NAtoms)
      DIMENSION IAT(NAtoms),IC(NAtoms,NAtoms)
      REAL*8 NND1,NND2
C
      DIMENSION RADII(103),NVAL(103)
C
C                   H     He    Li    Be    B     C     N     O
      DATA RADII/  0.60, 0.90, 1.40, 1.10, 0.90, 0.85, 0.75, 0.75,
C
C                   F     Ne    Na    Mg    Al    Si    P     S
     *             0.75, 0.95, 1.60, 1.45, 1.30, 1.20, 1.15, 1.10,
C
C                   Cl    Ar    K     Ca    Sc    Ti    V     Cr
     *             1.05, 1.05, 2.00, 1.80, 1.60, 1.50, 1.45, 1.40,
C
C                   Mn    Fe    Co    Ni    Cu    Zn    Ga    Ge
     *             1.40, 1.40, 1.35, 1.35, 1.40, 1.30, 1.35, 1.30,
C
C                   As    Se    Br    Kr    Rb    Sr    Y     Zr
     *             1.30, 1.25, 1.25, 1.22, 2.30, 2.20, 1.75, 1.65,
C
C                   Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
     *             1.55, 1.55, 1.45, 1.55, 1.40, 1.40, 1.45, 1.50,
C
C                   In    Sn    Sb    Te    I     Xe    Cs    Ba
     *             1.55, 1.45, 1.45, 1.42, 1.40, 1.35, 2.70, 2.45,
C
C                   La    Ce    Pr    Nd    Pm    Sm    Eu    Gd
     *             2.00, 1.90, 1.90, 1.85, 1.85, 1.85, 1.80, 1.78,
C
C                   Tb    Dy    Ho    Er    Tm    Yb    Lu    Hf
     *             1.78, 1.75, 1.75, 1.72, 1.72, 1.70, 1.72, 1.70,
C
C                   Ta    W     Re    Os    Ir    Pt    Au    Hg
     *             1.60, 1.60, 1.45, 1.45, 1.45, 1.45, 1.35, 1.40,
C
C                   Tl    Pb    Bi    Po    At    Rn    Fr    Ra
     *             1.60, 1.55, 1.58, 1.58, 1.60, 1.50, 3.00, 2.60,
C
C                   Ac    Th    Pa    U     Np    Pu    Am    Cm
     *             2.00, 1.80, 1.78, 1.75, 1.70, 1.70, 1.70, 1.90,
C
C                   Bk    Cf    Es    Fm    Md    No    Lr
     *             1.85, 1.80, 1.90, 1.78, 1.75, 1.75, 1.75 /
C
C
C  NVAL - maximum valency of each atom (normal)
C         (The number of "bonds" can exceed this)
C
C                H  He Li Be B  C  N  O  F  Ne Na Mg Al Si P  S  Cl Ar
      DATA NVAL/ 1, 4, 1, 2, 3, 4, 3, 2, 1, 4, 1, 2, 3, 4, 3, 2, 1, 4,
C
C                K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
     *           1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1, 4,
C
C                Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
     *           1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1, 4,
C
C                Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf
     *           1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
C
C                Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th
     *           5, 6, 5, 4, 3, 2, 1, 2, 3, 2, 3, 2, 1, 4, 1, 2, 3, 4,
C
C                Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr
     *           4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3 /
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
      R = SQRT( (XC(1,IA) - XC(1,JA))**2 +
     $          (XC(2,IA) - XC(2,JA))**2 +
     $          (XC(3,IA) - XC(3,JA))**2 )
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
C  Check that atomic valency is not zero, i.e., atom is not bonded
C  At same time check certain atoms to see if default valency
C  has been exceeded. (Currently N, P and S.)
C
      ERROR=.FALSE.
      DO 15 IA=1,NAtoms
      IAtm = IAN(IA)
      IF(NBOND(IA).EQ.0) THEN
       WRITE(6,1000) IA
       ERROR=.TRUE.
      ELSE IF(NBOND(IA).GT.IVAL(IA)) THEN
       If(IAtm.EQ.7.OR.IAtm.EQ.15) Then
        IVAL(IA) = 5
       Else If(IAtm.EQ.16) Then
        IVAL(IA) = 4
        If(NBOND(IA).GT.4) IVAL(IA) = 6
       EndIf
      ENDIF
 15   CONTINUE
c
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
        If(IVAL(JA).LT.0) IVAL(JA) = 0
        IVAL(IA) = 0
        NBOND(IA) = 0
        NBOND(JA) = NBOND(JA) - 1
       ELSE IF(NBOND(JA).EQ.1) THEN
        DONE=.FALSE.
        IC(IA,JA) = IVAL(JA)
        IC(JA,IA) = IVAL(JA)
        IVAL(IA) = IVAL(IA) - IVAL(JA)
        If(IVAL(IA).LT.0) IVAL(IA) = 0
        IVAL(JA) = 0
        NBOND(JA) = 0
        NBOND(IA) = NBOND(IA) - 1
       ENDIF
      ENDIF
 25   CONTINUE
cc      write(6,*) ' IC matrix after Loop 25 is:'
cc      do i=1,natoms
cc      write(6,*) ' number of bonds to atom:',i,' is:',nbond(i)
cc      write(6,*) ' IVal to atom:',i,' is:',ival(i)
cc      do j=1,natoms
cc      if(IC(i,j).ne.0) write(6,*) i,j,ic(i,j)
cc      enddo
cc      enddo
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
 30   CONTINUE
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
C
C  .....................
C  Check for amide bonds
C  .....................
C
      DO 60 IA=1,NAtoms
      IF(IAN(IA).EQ.7) THEN
        CALL GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
        If(NIC.NE.3) GO TO 60
c -- found a nitrogen atom attached to 3 other atoms
c -- is it attached to a carbon?
        DO 59 JA=1,NAtoms
        If(IC(IA,JA).NE.0.AND.IAN(JA).EQ.6) Then
c -- we have found a carbon
c -- is it attached to a double-bond oxygen?
          amide=.FALSE.
          DO 58 KA=1,NAtoms
          If(IC(JA,KA).EQ.2.AND.IAN(KA).EQ.8) Then
            amide=.TRUE.
            exit
          EndIf
 58       CONTINUE
          If(amide) Then
            IC(IA,JA) = 4
            IC(JA,IA) = 4
          EndIf
        EndIf
 59     CONTINUE
      ENDIF
 60   CONTINUE
C
      RETURN
c
 1000 FORMAT(' Atom ',I2,' does not appear to be bonded to anything ')
 1200 FORMAT(/,2X,'***ERROR*** Unable to determine connectivities.',
     $            ' Data must be read in')
c
      END
