c ==================================================================
c  FORCE FIELD MODULE            JB   Jan 2000/May 2005/
c                                     July 2008/Jan 2012
c  This routine is called for all force fields to determine the
c  force field type and call the corresponding driving routine
c ==================================================================
c
      subroutine preff(inp,ffcyc)
      implicit real*8(a-h,o-z)
      character*256 chopval,jobname,ff
c
c  reads the FFLD line in the input file and writes options
c  (if any) to the <control> file
c
      parameter (nopt=6)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
      integer ffcyc
      character*4 options(nopt)
      character cdum*20
c
      parameter (IUnit=1)
c
      Common /job/jobname,lenJ
c
      data options/'ffld','hess','preo','cuto','file','prin'/
      data ioptyp/21,0,0,11,21,1/
c
c
      call izeroit(iopval,3*nopt)
      call zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c ..............................................................
      If(ffcyc.gt.0) RETURN        ! only read first time
c ..............................................................
c
c -- what type of force field?
      if(ifound(1).eq.1) then
        if(chopval(1)(1:1).eq.' ') then
          ff = 'sybyl_5.2'
        else
          ff = chopval(1)
        endif
      endif
c put this in the depository
      call lowercas(ff,20)
      call setchval('ffld',ff)
c
c -- set calculation type
c     0 - just calculate E+G
c     1 - E+G + Hessian
c    10 - preoptimization + E+G
c    11 - preoptimization + E+G + Hessian
c
      IType = 0
c
c -- hessian calculation?
      if(ifound(2).eq.1) IType = IType+1
c
c -- preoptimization?
      if(ifound(3).eq.1) IType = IType+10
c
c -- cutoff for Van der Waals interactions
      CutOff = 10.0d0      ! default 10 A
      if(ifound(4).eq.1) CutOff = ropval(1,4)
c
c  -- additional force field parameters from file?
c  -- (copy to generic file)
c
      If(ifound(5).eq.1) Then
       call rmblan(chopval(5),256,Len)
       Call CopyFile(chopval(5)(1:Len),Len,jobname(1:lenJ)//'.fld',
     $               lenJ+4)
      EndIf
c
c -- print flag
      if(ifound(6).eq.1) then
        IPRNT = iopval(1,6)
      else
        IPRNT = 1          ! default print flag
      endif
c
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call RmBlank(256,ff,Len)
      call lowerca2(ff,Len)
      call wrcntrl(IUnit,6,'$basis',3,idum,dum,ff(1:Len))
      call wrcntrl(IUnit,7,'$fftype',1,IType,rdum,cdum)
      call wrcntrl(IUnit,11,'$vdw_cutoff',2,idum,CutOff,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
c
      return
      end
c  =======================================================================
c
      SUBROUTINE FFIELD(ffcyc,NMem,Z)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ...............................................................
C  ** FORCE FIELD PROGRAM **
C
C  This Program calculates both the energy and gradient and
C  (optionally) the Hessian for molecular mechanics force fields
C
C  Current supported force fields are:
C  1.  Sybyl v.5.2
C  2.  Universal Force Field (UFF)
C  3.  AMBER
C  4.  MMFF94 (under development)
C  ................................................................
C
C
      CHARACTER jobname*256,cdum*20,ff*20
      INTEGER ffcyc
      Logical ScrewYou
C
      DIMENSION Z(NMem)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
c
      Common /job/jobname,lenJ
C
C
C  Read from the <control> file
C    number of atoms
C    force field
C    calculation type
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
      call rdcntrl(IUnit,6,'$basis',3,idum,dum,ff)
      call rdcntrl(IUnit,7,'$fftype',1,IType,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c ............................................................
c -- for force field, dummy atoms are ignored
c
      NAtoms = NAtoms-Ndum1-Ndum2
c ............................................................
c
C  (over)estimate the number of bends and torsions/out-of-plane bends
C  [take a simple function of the number of atoms]
C
      NBend = 6*NAtoms
      NTors = 18*NAtoms
      NInv  = 12*NAtoms
C
C  ............................................................
C  Eventually need to overhaul this, getting better estimate
C  and removing NAtoms**2 dependence on distance parameters
C  ............................................................
C
C  Determine the memory for the appropriate force field
C
      ScrewYou = .False.
 50   CONTINUE
      If(ff(1:9).EQ.'sybyl_5.2') Then
        IMem = 11*NAtoms + 4*NAtoms**2 + 5*NBend + 6*NTors
      Else If(ff(1:3).EQ.'uff') Then
        IMem = 11*NAtoms + 4*NAtoms**2 + 8*NBend + 7*NTors + 8*NInv
      Else If(ff(1:5).EQ.'amber') Then
        IMem = 12*NAtoms + 5*NAtoms**2 + 5*NBend + 7*NTors
      Else
c -- looks like a QM/MM Job,just done QM part
c    get force field type from depository
        If(ScrewYou)
     $    call nerror(1,'Force Field module',
     $       'Unrecognizable Force Field - Check Input File',0,0)
        call getchval('ffld',ff)
        ScrewYou = .True.
        GO TO 50
      EndIf
C
C  Do we need a Hessian matrix?
C
      If(IType.EQ.1.OR.IType.EQ.11) IMem = IMem + 9*NAtoms**2
c
      iptr = 1
      IErr = NMem - IMem
      If(IErr.LT.0) CALL MemERR(8*IMem,6,'FFIELD')
      CALL ZeroIT(Z,IMem)             ! clear the memory
C
C
C  Now allocate the memory pointers
C  (1) ubiquitous storage
C
      IAN   = iptr                    !  atomic numbers
      IXC   = IAN  + NAtoms           !  geometry
      IGC   = IXC  + 3*NAtoms         !  gradient
      IAT   = IGC  + 3*NAtoms         !  forcefield atom types
      ICC   = IAT  + NAtoms           !  connectivity matrix
      IRD   = ICC  + NAtoms**2        !  distance matrix
      IRL   = IRD  + NAtoms**2        !  forcefield equilibrium distances
      IBB2  = IRL  + NAtoms**2        !  distance force constants
C
C  (2) force field storage
C
      IF(ff(1:9).EQ.'sybyl_5.2') THEN
cc
      IB1   = IBB2 + NAtoms**2        !  first atom in bend
      IB2   = IB1  + NBend            !  second atom in bend
      IB3   = IB2  + NBend            !  third atom in bend
      IANG  = IB3  + NBend            !  Sybyl equilibrium angles
      IBB3  = IANG + NBend            !  Sybyl bend force constants
      IT1   = IBB3 + NBend            !  first atom in torsion/outp
      IT2   = IT1  + NTors            !  second atom in torsion/outp
      IT3   = IT2  + NTors            !  third atom in torsion/outp
      IT4   = IT3  + NTors            !  fourth atom in torsion/outp
      IS4   = IT4  + NTors            !  Sybyl sign factor for torsion/outp
      IBB4  = IS4  + NTors            !  Sybyl torsion/outp force constants
      IEnd  = IBB4 + NTors
cc
      ELSE IF(ff(1:3).EQ.'uff') THEN
cc
      IB1   = IBB2 + NAtoms**2        !  first atom in bend
      IB2   = IB1  + NBend            !  second atom in bend
      IB3   = IB2  + NBend            !  third atom in bend
      IC0   = IB3  + NBend            !  first coefficient in bend energy term
      IC1   = IC0  + NBend            !  second coefficient in bend energy term
      IC2   = IC1  + NBend            !  third coefficient in bend energy term
      IBB3  = IC2  + NBend            !  UFF bend force constants
      INB   = IBB3 + NBend            !  integer bend order
      IT1   = INB  + NBend            !  first atom in torsion
      IT2   = IT1  + NTors            !  second atom in torsion
      IT3   = IT2  + NTors            !  third atom in torsion
      IT4   = IT3  + NTors            !  fourth atom in torsion
      ICOS  = IT4  + NTors            !  Cosine factor for torsion
      IBB4  = ICOS + NTors            !  UFF torsion force constants
      ITT   = IBB4 + NTors            !  integer torsion order
      IP1   = ITT  + NTors            !  first atom in inversion
      IP2   = IP1  + NInv             !  second atom in inversion
      IP3   = IP2  + NInv             !  third atom in inversion
      IP4   = IP3  + NInv             !  fourth atom in inversion
      IS0   = IP4  + NInv             !  first coefficient in inversion energy term
      IS1   = IS0  + NInv             !  second coefficient in inversion energy term
      IS2   = IS1  + NInv             !  third coefficient in inversion energy term
      IBB5  = IS2  + NInv             !  UFF inversion force constant
      IEnd  = IBB5 + NInv
cc
      ELSE IF(ff(1:5).EQ.'amber') THEN
cc
      IQC   = IBB2 + NAtoms**2        !  AMBER atomic charges
      ICP   = IQC  + NAtoms           !  AMBER electrostatic factors
      IB1   = ICP  + NAtoms**2        !  first atom in bend
      IB2   = IB1  + NBend            !  second atom in bend
      IB3   = IB2  + NBend            !  third atom in bend
      IANG  = IB3  + NBend            !  AMBER equilibrium angles
      IBB3  = IANG + NBend            !  AMBER bend force constants
      IT1   = IBB3 + NBend            !  first atom in proper/improper torsion
      IT2   = IT1  + NTors            !  second atom in proper/improper torsion
      IT3   = IT2  + NTors            !  third atom in proper/improper torsion
      IT4   = IT3  + NTors            !  fourth atom in proper/improper torsion
      ITOR  = IT4  + NTors            !  AMBER equilibrium torsion
      ITT   = ITOR + NTors            !  integer torsion order
      IBB4  = ITT  + NTors            !  AMBER torsion force constants
      IEnd  = IBB4 + NTors
      ENDIF
C
C  (3) Hessian storage
C
      IGU   = IEnd                    !  spare copy of gradient
      IEnd  = IGU  + 3*NAtoms
C
C  Do we need a Hessian matrix?
C
      If(IType.EQ.1.OR.IType.EQ.11) Then
       IHS   = IEnd                   !  Hessian matrix
       IEnd  = IHS  + (3*NAtoms)**2
      Else
       IHS = 1
      EndIf
C
C  Check memory storage not exceeded
C
      IEnd = IEnd - iptr
      CALL MemCHK(IMem,IEnd,6,'FFIELD')
c
c -- set memory mark for TEXAS memory manager
      call getmem(IMem,lastx)
C
C
C  ----------------------------------------------------------------------
C
      IF(ff(1:9).EQ.'sybyl_5.2') THEN
        CALL SYBYL_MAIN(NAtoms,  Z(IAN),  IType,   NBend,   NTors,
     $                  Z(IXC),  Z(ICC),  Z(IAT),  Z(IRD),  Z(IRL),
     $                  Z(IBB2), Z(IB1),  Z(IB2),  Z(IB3),  Z(IANG),
     $                  Z(IBB3), Z(IT1),  Z(IT2),  Z(IT3),  Z(IT4),
     $                  Z(IS4),  Z(IBB4), Z(IGU),  Z(IGC),  Z(IHS),
     $                  ffcyc)
      ELSE IF(ff(1:3).EQ.'uff') THEN
        CALL UFF_MAIN(NAtoms,  Z(IAN),  IType,   NBend,   NTors,
     $                NInv,    Z(IXC),  Z(ICC),  Z(IAT),  Z(IRD),
     $                Z(IRL),  Z(IBB2), Z(IB1),  Z(IB2),  Z(IB3),
     $                Z(IC0),  Z(IC1),  Z(IC2),  Z(IBB3), Z(INB),
     $                Z(IT1),  Z(IT2),  Z(IT3),  Z(IT4),  Z(ICOS),
     $                Z(IBB4), Z(ITT),  Z(IP1),  Z(IP2),  Z(IP3),
     $                Z(IP4),  Z(IS0),  Z(IS1),  Z(IS2),  Z(IBB5),
     $                Z(IGU),  Z(IGC),  Z(IHS),  ffcyc)
      ELSE IF(ff(1:5).EQ.'amber') THEN
        CALL AMBER_MAIN(NAtoms,  Z(IAN),  IType,   NBend,   NTors,
     $                  Z(IXC),  Z(ICC),  Z(IAT),  Z(IRD),  Z(IRL),
     $                  Z(IBB2), Z(IQC),  Z(ICP),  Z(IB1),  Z(IB2),
     $                  Z(IB3),  Z(IANG), Z(IBB3), Z(IT1),  Z(IT2),
     $                  Z(IT3),  Z(IT4),  Z(ITOR), Z(IBB4), Z(ITT),
     $                  Z(IGU),  Z(IGC),  Z(IHS),  ffcyc)
      ENDIF
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
      END
c  =======================================================================
c
*Deck vecmul
      SUBROUTINE VECMUL(R1,R2,R3)
      REAL*8 R1(3),R2(3),R3(3)
c
      R3(1) =  R1(2)*R2(3) - R1(3)*R2(2)
      R3(2) = -R1(1)*R2(3) + R1(3)*R2(1)
      R3(3) =  R1(1)*R2(2) - R1(2)*R2(1)
c
      RETURN
      END
c  =======================================================================
c
*Deck chkaromatic
      SUBROUTINE ChkAromatic(NAtoms,IAN,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This subroutine checks total bond orders for potentially
C  aromatic atoms (carbon and nitrogen) to make sure bonds
C  have not been incorrectly assigned as double bonds
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  IC      -  connectivity matrix
C
      DIMENSION IAN(NAtoms),IC(NAtoms,NAtoms)
      Dimension BondOrder(5),IBond(4),JBond(4)
      Logical Change
c
      Data BondOrder/1.0d0, 2.0d0, 3.0d0, 1.0d0, 1.5d0/
      Parameter (half=0.5d0,zero=0.0d0)
C
C
      DO 50 IA=1,NAtoms
c
      IAtNo = IAN(IA)
      If(IAtNo.EQ.6) Then
        Val = 4.5d0
      Else If(IAtNo.EQ.7) Then
        Val = 3.0d0
      Else
        GO TO 50
      EndIf
c
      Nbond = 0
      bo = zero
      DO JA=1,NAtoms
      If(IC(JA,IA).GT.0) Then
        bo = bo + BondOrder(IC(JA,IA))
        Nbond = Nbond+1
        If(NBond.GT.4) GO TO 50    ! too many bonds - all bets off for this atom
        IBond(Nbond) = JA          ! atom bonded to IA
        JBond(Nbond) = IC(JA,IA)   ! bond type
      EndIf
      EndDO
c
 10   CONTINUE
c
c -- at this point we have either C or N and the total valency
c
      IF( (IAtNo.EQ.6.OR.IAtNo.EQ.7) .AND. bo.GT.Val) THEN
c
c -- find a double bond and change it to aromatic if possible
        Change = .False.
        DO JA=1,Nbond
        If(JBond(JA).EQ.2) Then
          JJ = IAN(IBond(JA))
          If(JJ.EQ.6.OR.JJ.EQ.7) Then
c -- change bond to aromatic
            IC(IA,IBond(JA)) = 5
            IC(IBond(JA),IA) = 5
            bo = bo-half
            Change = .True.
          EndIf
        EndIf
        EndDO
c
c -- go back and check again
        If(Change) GO TO 10
c
      ENDIF
c
c -- for carbon, bond order of 4.5 only possible if we have
c    three bonds, and all three are aromatic
c    change any remaining double bonds to aromatic
      If(IAtNo.EQ.6.AND.bo.GT.4.4) Then
        DO JA=1,Nbond
        If(JBond(JA).EQ.2) Then
          IC(IA,IBond(JA)) = 5
          IC(IBond(JA),IA) = 5
        EndIf
        EndDO
      EndIf
c
 50   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
*Deck getlic
      SUBROUTINE GetLIC(NAtoms,IC,IA,LIC,NIC,MIC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Obtains various bonding parameters for the requested atom
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IC      -  connectivity matrix
C  IA      -  atom for which data is required
C
C  on exit
C
C  LIC     -  maximum bond type
C             standard bond types are: 1 = single bond; 2 = double bond;
C               3 = triple bond; 4 = amide bond; 5 = aromatic bond
C  NIC     -  total number of bonds
C  MIC     -  sum of formal bond orders
C             formal bond orders are: single bond = 1; double bond = 2;
C               triple bond = 3; amide bond = 1; aromatic bond = 1.5
C
C
      DIMENSION IC(NAtoms,NAtoms)
      Dimension BondOrder(5)
      Data BondOrder/1.0d0, 2.0d0, 3.0d0, 1.0d0, 1.5d0/
C
      LIC = 0
      NIC = 0
      bond = 0.0d0
c
      DO 10 JA=1,NAtoms
      IJ = IC(JA,IA)
      If(IJ.GT.0) Then
        NIC = NIC+1
        bond =  bond + BondOrder(IJ)
        If(IJ.GT.LIC) LIC=IC(JA,IA)
      EndIf
 10   CONTINUE
c
      MIC = NINT(bond)
c
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE RdPQBFF(IUnit,FFType,NAtoms,Nqm,IAT,IC)
      IMPLICIT INTEGER(A-Z)
C
C  Reads the <pqb> file, extracts the force field atomic symbols,
C  converts them to atom types and reads atomic connectivity data
C
C  ARGUMENTS
C
C  IUnit   -  Unit number of <pqb> file (should already be opened)
C  FFType  -  name of force field
C             currently either Sybyl_5.2, UFF or AMBER
C  NAtoms  -  number of atoms in current system
C  Nqm     -  number of atoms in QM part (in QMMM job)
C              (if full system, Nqm should be set to NAtoms)
C  IAT     -  on exit force field atom types
C  IC      -  on exit connectivity matrix with bonding types
C
      DIMENSION IAT(NAToms),IC(NAtoms,NAtoms)
      CHARACTER FFType*9, Char*20
      CHARACTER*5 SYBYL(32),UFF(127),AMBER(40),ffS
c
      DATA SYBYL/'C.3  ','C.2  ','C.ar ','C.1  ','N.3  ','N.2  ',
     $           'N.1  ','O.3  ','O.2  ','S.3  ','N.ar ','P.3  ',
     $           'H    ','Br   ','Cl   ','F    ','I    ','S.2  ',
     $           'N.pl3','LP   ','Na   ','K    ','Ca   ','Li   ',
     $           'Al   ','Du   ','Si   ','N_am ','S.O  ','S.O2 ',
     $           'N.4  ','     '/
c
      DATA UFF/'H_   ','H_b  ','He4+4','Li   ','Be3+2','B_3  ','B_2  ',
     $         'C_3  ','C_R  ','C_2  ','C_1  ','N_3  ','N_R  ','N_2  ',
     $         'N_1  ','O_3  ','O_3_z','O_R  ','O_2  ','O_1  ','F_   ',
     $         'Ne4+4','Na   ','Mg3+2','Al3  ','Si3  ','P_3+3','P_3+5',
     $         'P_3+q','S_3+2','S_3+4','S_3+6','S_R  ','S_2  ','Cl   ',
     $         'Ar4+4','K_   ','Ca6+2','Sc3+3','Ti3+4','Ti6+4','V_3+5',
     $         'Cr6+3','Mn6+2','Fe3+2','Fe6+2','Co6+3','Ni4+2','Cu3+1',
     $         'Zn3+2','Ga3+3','Ge3  ','As3+3','Se3+2','Br   ','Kr4+4',
     $         'Rb   ','Sr6+2','Y_3+3','Zr3+4','Nb3+5','Mo6+6','Mo3+6',
     $         'Tc6+5','Ru6+2','Rh6+3','Pd4+2','Ag1+1','Cd3+2','In3+3',
     $         'Sn3  ','Sb3+3','Te3+2','I_   ','Xe4+4','Cs   ','Ba6+2',
     $         'La3+3','Ce6+3','Pr6+3','Nd6+3','Pm6+3','Sm6+3','Eu6+3',
     $         'Gd6+3','Tb6+3','Dy6+3','Ho6+3','Er6+3','Tm6+3','Yb6+3',
     $         'Lu6+3','Hf3+4','Ta3+5','W_6+6','W_3+4','W_3+6','Re6+5',
     $         'Re3+7','Os6+6','Ir6+3','Pt4+2','Au4+3','Hg1+2','Tl3+3',
     $         'Pb3  ','Bi3+3','Po3+2','At   ','Rn4+4','Fr   ','Ra6+2',
     $         'Ac6+3','Th6+4','Pa6+4','U_6+4','Np6+4','Pu6+4','Am6+4',
     $         'Cm6+3','Bk6+3','Cf6+3','Es6+3','Fm6+3','Md6+3','No6+3',
     $         'Lw6+3'/
C
      DATA AMBER/'CT   ','C    ','CA   ','CM   ','CC   ','CV   ',
     $           'CW   ','CR   ','CB   ','C*   ','CN   ','CK   ',
     $           'CQ   ','N    ','NA   ','NB   ','NC   ','N*   ',
     $           'N2   ','N3   ','OW   ','OH   ','OS   ','O    ',
     $           'O2   ','S    ','SH   ','P    ','H    ','HW   ',
     $           'HO   ','HS   ','HA   ','HC   ','H1   ','H2   ',
     $           'H3   ','HP   ','H4   ','H5   '/
C
C
C  Sort out the force field
C
      WRITE(6,*) ' Force Field data read from <pqb> file'
c
      NEl = 32                            ! default is Sybyl
      If(FFType(1:3).EQ.'UFF') NEl=127
      If(FFType(1:5).EQ.'Amber') NEl=40
C
C  Skip the first two lines
C
      READ(IUnit,910) Char
      READ(IUnit,910) Char
c
      DO 30 IAtm=1,Nqm
      If(NEl.EQ.32) Then
        READ(IUnit,900) ffS
      Else If(NEl.EQ.127) Then
        READ(IUnit,901) ffS
      Else
        READ(IUnit,902) ffS
      EndIf
c
      DO 10 I=1,NEl
      If(NEl.EQ.32) Then
        If(ffS.EQ.SYBYL(I)) GO TO 20
      Else If(NEl.EQ.127) Then
        If(ffS.EQ.UFF(I)) GO TO 20
      Else
        If(ffS.EQ.AMBER(I)) GO TO 20
      EndIf
 10   CONTINUE
      I = NEl
c
 20   CONTINUE
      IAT(IAtm) = I
 30   CONTINUE
C
C  Now read the connectivity data
C
      CALL IZeroIT(IC,NAtoms**2)
c
 35   CONTINUE
      READ(IUnit,910,END=95) Char
      If(Char(1:4).NE.'$con') GO TO 35
c
 40   CONTINUE
      READ(IUnit,910,END=96) Char
      If(Char(1:4).EQ.'$end') RETURN
c
      READ(Char(1:20),*) I,J,IB
      If(I.GT.Nqm.OR.J.GT.Nqm) GO TO 40
      IC(I,J) = IB
      IC(J,I) = IB
      GO TO 40
c
C -- ERROR SECTION
 95   CONTINUE
      CALL nerror(1,'<pqb> file',
     $         'Problem reading atom connectivity data',0,0)
c
 96   CONTINUE
      CALL nerror(2,'<pqb> file','Reached End of File in Error',0,0)
c
  900 Format(56X,A5)
  901 Format(63X,A5)
  902 Format(70X,A5)
  910 Format(A20)
c
      END
c  =======================================================================
c
      SUBROUTINE AddLINK(FFType,NAtoms,NLink,XC,IAT,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Adds Link Hydrogen atoms to the atom types and connectivity
C  for the QM part in a QM/MM calculation
C
C  ARGUMENTS
C
C  FFType  -  name of force field
C             currently either Sybyl_5.2 or UFF
C  NAtoms  -  number of atoms in model system
C  NLink   -  number of link atoms
C              (NOTE: number of atoms in QM part is  NAtoms-NLink
C  XC      -  Cartesian coordinates of model system
C  IAT     -  on exit appended force field atom types
C  IC      -  on exit appended connectivity matrix
C
      DIMENSION XC(3,NAtoms),IAT(NAToms),IC(NAtoms,NAtoms)
      CHARACTER*9 FFType
C
      PARAMETER (nH_Syb=13,nH_UFF=1)
C  
C
C  First add the Hydrogens
C
      Nqm = NAtoms - NLink         ! no. of QM atoms
c
      DO 10 I=Nqm+1,NAtoms
      If(FFType(1:3).EQ.'UFF') Then
        IAT(I) = nH_UFF
      Else
        IAT(I) = nH_Syb
      EndIf
 10   CONTINUE
C
C  Now determine which atoms the link atoms are connected to
C  (use closest distance)
C
      DO 30 IAtm=Nqm+1,NAtoms
c
      DSmall = 999999.0d0
      DO 20 I=1,Nqm
      Dist2 =   (XC(1,IAtm)-XC(1,I))**2 + (XC(2,IAtm)-XC(2,I))**2
     $        + (XC(3,IAtm)-XC(3,I))**2
      If(Dist2.LT.DSmall) Then
        DSmall = Dist2
        LL = I
      EndIf
 20   CONTINUE
C
C  Link H atom IAtm is connected to main system atom LL
C
      IC(IAtm,LL) = 1
      IC(LL,IAtm) = 1
c
 30   CONTINUE
C
      RETURN
      END
