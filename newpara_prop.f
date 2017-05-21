c ==================================================================
c  PARALLEL PROPERTY MODULE            JB   March 2001
c ==================================================================
c
      SUBROUTINE Para_PROPMAIN(NAtoms, IAN,    QA,     XC,     NBas,
     $                         NShell, nsh,    BASDAT, INX,    factor,
     $                         r0f,    LMax,   NRad,   NAng,   ExpMIN,
     $                         BL,     AIJ,    DISTN,  RDIST,  radii,
     $                         radwght,XXA,    WTA,    NAlpha, NBeta,
     $                         rhf,    Ylm,    CMO,    CA,     CMOB,
     $                         CB,     NTrans, IUNQ,   TRANS,  NEqATM,
     $                         SIGNA,  SIGNB,  XCharg, XSpin,  XEFG,
     $                         NScr,   Z,      IErr)

      use memory, block => bl
      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Secondary "wrapper" for property module
C  PROPMAIN is responsible for all file I/O and tidying
C  up prior to exit
C
C
      DIMENSION IAN(NAtoms),QA(NAtoms),XC(3,NAtoms),BASDAT(13,*),
     $          INX(12,*),ExpMIN(NShell),BL(*),DISTN(NAtoms),
     $          AIJ(NAtoms,NAtoms),RDIST(NAtoms,NAtoms),
     $          radii(*),radwght(*),XXA(3,1130),WTA(1130),
     $          Ylm(0:LMax,-LMax:LMax),CMO(NAlpha,NBas),
     $          CA(NAlpha,0:LMax,-LMax:LMax),CMOB(NBeta,NBas),
     $          CB(NAlpha,0:LMax,-LMax:LMax)
      DIMENSION IUNQ(NAtoms),TRANS(9,NAtoms),NEqATM(NAtoms,NTrans)
      INTEGER   SIGNA(NAlpha,7),SIGNB(NAlpha,7)
      DIMENSION XCharg(NAtoms,2),XSpin(NAtoms,2),XEFG(9,NAtoms)
c ..................................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
c ..................................................
      DIMENSION RM(3,3)
      CHARACTER GROUP*4,cdum*20
      Character*256 jobname,MOS,MOB
      Logical Symflag,abelian,rhf,spin,efg,full
c  ...................................................................
c     common /intbl/maxsh,ixx(100)     ! old texas integer depository
      common /symm/nsym,nsy(7)         ! old texas symmetry data
c ....................................................................
      DIMENSION Z(NScr)
c
      PARAMETER (Zero=0.0d0,One=1.0d0,Three=3.0d0)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      Data IradQ/0/
c
      Common /job/jobname,lenJ
C
C
C  initialize
C
      IOut=ioutfil('iout')
      ICon=ioutfil('icon')
      thrsh = 1.0d-12
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
      if(lenM.le.76) MOS(lenM+1:lenM+4)='.mos'
      lenM=lenM+4
      MOB=MOS
      MOB(lenM:lenM)='b'
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
C  determine if system has abelian symmetry
C
      abelian = GROUP.EQ.'cs  '.OR.GROUP.EQ.'ci  '.OR.
     $          GROUP.EQ.'c2  '.OR.GROUP.EQ.'c2v '.OR.
     $          GROUP.EQ.'c2h '.OR.GROUP.EQ.'d2  '.OR.
     $          GROUP.EQ.'d2h '
cc      abelian = .true.
C
C  Read from the <control> file
C     print flag
C     which properties to be computed
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      call rdcntrl(IUnit,5,'$prop',3,idum,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  which properties?
C
      spin = cdum(1:4).eq.'spin'
      efg = cdum(1:3).eq.'efg'.OR.cdum(5:7).eq.'efg'
C
C  print header
C
      WRITE(IOut,1000)
      WRITE(ICon,1000)
c
      If(IPrnt.gt.1) Then
       write(iout,*) ' Rassolov/Chipman radial factor: ',r0f
       write(iout,*) ' maximum angular order (LMax):    ',LMax
       write(iout,*) ' radial grid factor:             ',factor
       write(iout,*) ' properties to be computed:  ',cdum(1:4),
     $               '  ',cdum(5:7)
      EndIf
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XC,-1,jnk)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C  Read the basis set data
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdbasis(IUnit,  NAtoms, AtSymb, XC,     Z(1),
     $             Z(1+NAtoms),BASDAT)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C
C ....................................................................
C  Be Careful with basis order!
C  We will use Wolinski order for property module
C  Wolinski ordering (for integral evaluation and hence for the MOs)
C  is PER SHELL (e.g. all S functions together regardless of atom)
C
C  The following routines are relevant here
C   SortBAS1  -  sorts basis per shell (but NOT full Wolinski)
C   reorder   -  orders basis from SortBAS1 into Wolinski
C ....................................................................

      CALL SortBAS1(NAtoms,NShell,Z(1+NAtoms),INX)
      CALL normaliz(NShell,INX,BASDAT)
C
C  get atomic numbers from atomic symbols
C
      CALL GetAtNo(NAtoms,AtSymb,IAN)
C
C  get actual atomic charges (there may, e.g., be ghost atoms)
C
      CALL GetAtChrg(NAtoms,QA)
C
C  need to reorder basis functions to Wolinski order
C
      i1 = 1
      i2 = i1 + NBas
      i3 = i2 + NShell
      i4 = i3 + 12*NShell
      IEnd = i4 + NBas*NBas - 1
      CALL MemCHK(NScr,IEnd,8,'PROPMAIN')
c
      Call reorder(NShell,INX,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinski order
C
C ================================================================
C  get location of abelian symmetry data in old texas depository
C
      call getival('ngener',ngen)         ! number of generators
      nupr = 1                            ! not set if symmetry not used
      ifp = 1                             ! ditto
      If(nsym.gt.0) Then
       call getival('SymNuPr1',nupr)      ! equivalent atoms
       call getival('SymFunPr',ifp)       ! equivalent basis functions
      EndIf
C ================================================================
C
C  now read the old MOs (from the existing binary file)
C  and determine how MOs transform over the abelian symmetry
C  operations (i.e., whether they change sign or not)
C
      itype = 1
      CALL ReadMOS(NBas,Z(i4),jnk,.False.,lenM,MOS(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(1,'Properties module',
     $   'MOS File Does Not Exist',0,0)
c
      If(abelian)
     $ CALL FindMOSym(NBas,   NAlpha, NSym,   nsy,    block(ifp),
     $                Z(i4),  Z(i1),  Z(i2),  SIGNA,  full)
C
C  now transpose MOs
C
      call trspmo(Z(i4),NBas,CMO,NAlpha)
C
C  repeat for beta MOs (if unrestricted wavefunction)
C
      If(NBeta.GT.0) Then
        CALL ReadMOS(NBas,Z(i4),jnk,.False.,lenM,MOB(1:lenM),itype,IErr)
        If(IErr.NE.0) Call nerror(1,'Properties module',
     $     'MOB File Does Not Exist',0,0)
c
        If(abelian.AND.full)
     $   CALL FindMOSym(NBas,   NBeta,  NSym,   nsy,  block(ifp),
     $                  Z(i4),  Z(i1),  Z(i2),  SIGNB,  full)
        call trspmo(Z(i4),NBas,CMOB,NBeta)
      EndIf
C
C
C  Precalculate data for the numerical grid
c  ----------------------------------------
c -- initialize slaves and send over preliminary data
      If(nslv.GT.NQ) Then
c -- inefficiency warning
        call message('**WARNING** in prop',
     $  'Not possible to use slaves efficiently  nslv>NQ:',nslv,NQ)
      EndIf
c
c -- check if the slaves are OK
      call f_lush(6)
      call para_check
c
c -- tell the slaves that Property calculations are next
      call para_bcast(TxDoProp,TxJobType)
c
c -- send initial data to slaves
      call para_initsend
      call para_pack_int(NAlpha,1)
      call para_pack_int(NBeta,1)
      call para_pack_int(rhf,1)
      call para_pack_int(spin,1)
      call para_pack_int(efg,1)
      call para_pack_int(NRad,1)
      call para_pack_int(NAng,1)
      call para_pack_int(LMax,1)
      call para_pack_real(factor,1)
      call para_pack_real(thrsh,1)
c
      call para_bcast_pack(TxPropInit)
c
c -- send the remaining basic data necessary (geometry,symmetry etc...)
      call para_send_info
c
c -- send MO data to slaves
      call para_bcast_real(CMO,NAlpha*NBas,TxPropDat)
      If(NBeta.GT.0) call para_bcast_real(CMOB,NBeta*NBas,TxPropDat)
C
C  switch off symmetry for charge/spin density calculation for
C  non-abelian point group and if MO symmetry analysis incomplete
C
      MSym = 0
      If(abelian.AND.full) MSym=NSym
C
C
C  Now calculate the charge and spin densities
C
C  ------------------------------------------------------------
      IF(spin) THEN
C
C  Charge/Spin density
C  -------------------
C
        Call ZeroIT(XCharg,2*NAtoms)
        If(NBeta.GT.0) CALL ZeroIT(XSpin,2*NAtoms)
c
c -- check if the slaves are OK
        call para_check
c
c -- send charge/spin data to slaves
        call para_initsend
        call para_pack_int(MSym,1)
        call para_pack_real(r0f,1)
        call para_pack_real(SIGNA,NAlpha*7)
        If(NBeta.GT.0) call para_pack_real(SIGNB,NBeta*7)
        call para_bcast_pack(TxSpinDat)
c
c -- loop over atomic centres, sending each atom to a different slave
        call para_distr_atoms(iunq,nq)
c
c -- now sum up charge and spin densities from each slave
        ms = 2*NAtoms  ! size of arrays
        call para_reduce(XCharg,ms,TxProp1)
        If(NBeta.GT.0) call para_reduce(XSpin,ms,TxProp2)
c
c -- That's it!
c
        ENDIF
C
C  print out results
C  -----------------
C
        IF(NTrans.GT.1) THEN
          DO 200 IAtm=1,NQ
          ICntr = IUNQ(IAtm)
C
C  generate all symmetry-related centres from current
C  symmetry-unique centre
C
          DO 190 IOP=1,NTrans
          JAtm = NEqATM(ICntr,IOP)
          If(JAtm.NE.ICntr) Then
            XCharg(JAtm,1) = XCharg(ICntr,1)
            XCharg(JAtm,2) = XCharg(ICntr,2)
            If(.NOT.rhf) Then
              XSpin(JAtm,1) = XSpin(ICntr,1)
              XSpin(JAtm,2) = XSpin(ICntr,2)
            EndIf
          EndIf
 190      CONTINUE
 200      CONTINUE
C
C  now print out the results
C
        WRITE(IOut,1100)
        WRITE(ICon,1100)
c
        If(rhf) Then
          WRITE(IOut,1200)
          WRITE(ICon,1200)
        Else
          WRITE(IOut,1300)
          WRITE(ICon,1300)
        EndIf
c
        DO 300 IAtm=1,NAtoms
        If(rhf) Then
          WRITE(IOut,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2)
          WRITE(ICon,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2)
        Else
          WRITE(IOut,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2),
     $                   XSpin(IAtm,1),XSpin(IAtm,2)
          WRITE(ICon,1400) IAtm,AtSymb(IAtm),
     $                   XCharg(IAtm,1),XCharg(IAtm,2),
     $                   XSpin(IAtm,1),XSpin(IAtm,2)
        EndIf
 300    CONTINUE
cc
      ENDIF
C
      IF(efg) THEN
C
C  Electric Field Gradient
C  -----------------------
C
        CALL ZeroIT(XEFG,9*NAtoms)
c
c *******************************************
c -- Symmetry not working properly
c -- switch off until fixed
        MSym = NSym
        call TempSYM(NAtoms,NSym,NGen,NQ,IUNQ)
c *******************************************
C
c -- check if the slaves are OK
        call para_check
c
c -- send EFG data to slaves
        call para_initsend
        call para_pack_int(IradQ,1)
        call para_pack_int(NQ,1)
        call para_pack_int(IUNQ,NAtoms)
        call para_bcast_pack(TxEFGDat)
c
c -- loop over atomic centres, sending each atom to a different slave
        call para_distr_atoms(iunq,nq)
c
c -- now sum up EFG from each slave
        ms = 9*NAtoms  ! size of arrays
        call para_reduce(XEFG,ms,TxProp3)
c
c -- That's it!
c
C
C  print out results
C  -----------------
C
        IF(NSym.GT.0) THEN
          DO 500 IAtm=1,NQ
          ICntr = IUNQ(IAtm)  ! WARNING - This is WRONG  needs to be abelian
C
C  scale EFG tensor by symmetry factor
C
          CALL RedEFG(NAtoms,ICntr,NGen,bl(nupr),XEFG(1,ICntr))
C
C  generate all symmetry-related centres from current
C  symmetry-unique centre
C
          CALL SymEFG(NAtoms,ICntr,NGen,nsy,bl(nupr),XEFG)
c
 500      CONTINUE
        ENDIF
c
c *******************************************
c -- restore real NSym
        NSym = MSym
c *******************************************
C
C  now print out the results
C
        WRITE(IOut,1500)
        WRITE(ICon,1500)
c
cc        DO 600 IAtm=1,NAtoms
cc        XEFG(2,IAtm) = Three*XEFG(2,IAtm)
cc        XEFG(4,IAtm) = Three*XEFG(4,IAtm)
cc        XEFG(5,IAtm) = Three*XEFG(5,IAtm)
c
cc        WRITE(IOut,1600) IAtm,AtSymb(IAtm),
cc     $                   XEFG(7,IAtm),XEFG(8,IAtm),XEFG(9,IAtm),
cc        WRITE(ICon,1600) IAtm,AtSymb(IAtm),
cc     $                   XEFG(7,IAtm),XEFG(8,IAtm),XEFG(9,IAtm),
cc     $                   XEFG(2,IAtm),XEFG(4,IAtm),XEFG(5,IAtm)
cc 600    CONTINUE
C
C  Now do Tensor representation
C
        DO 650 IAtm=1,NAtoms
        XEFG(2,IAtm) = Three*XEFG(2,IAtm)
        XEFG(4,IAtm) = Three*XEFG(4,IAtm)
        XEFG(5,IAtm) = Three*XEFG(5,IAtm)
c
        WRITE(IOut,1600) IAtm,AtSymb(IAtm),
     $                   XEFG(1,IAtm),XEFG(3,IAtm),XEFG(6,IAtm),
     $                   XEFG(2,IAtm),XEFG(4,IAtm),XEFG(5,IAtm)
        WRITE(ICon,1600) IAtm,AtSymb(IAtm),
     $                   XEFG(1,IAtm),XEFG(3,IAtm),XEFG(6,IAtm),
     $                   XEFG(2,IAtm),XEFG(4,IAtm),XEFG(5,IAtm)
 650    CONTINUE
C
C  diagonalize tensor for each atom and print out eigenvalues
C
        WRITE(IOut,1800)
        WRITE(ICon,1800)
c
        I1 = 4
        I2 = I1 + 9
        I3 = I2 + 9
c
        DO 700 IAtm=1,NAtoms
        CALL EXPAND(3,XEFG(1,IAtm),Z(I1))
        CALL DiagMAT(Z(I1),3,Z(I2),Z(I3),Z,IErr)
c
        If(IErr.NE.0) Then
          WRITE(IOut,2000) IAtm
          GO TO 700
        EndIf
c
        WRITE(IOut,1600) IAtm,AtSymb(IATM),Z(1),Z(2),Z(3)
        WRITE(ICon,1600) IAtm,AtSymb(IATM),Z(1),Z(2),Z(3)
 700    CONTINUE
cc
      ENDIF
C  ------------------------------------------------------------
C
      Call para_next(0)
C
      RETURN
c
 1000 FORMAT(/,'                          ** NUCLEAR PROPERTIES **')
 1100 FORMAT(/,' <== CHARGE AND SPIN DENSITY AT THE NUCLEUS',
     $          '  (atomic units)  ==>')
 1200 FORMAT(/,'                    CHARGE  DENSITY',/,
     $         '    Atom     Delta Function  Rassolov/Chipman')
 1300 FORMAT(/,20X,'CHARGE  DENSITY',20X,'SPIN  DENSITY',/,
     $         '  Atom     Delta Function  Rassolov/Chipman',
     $         '  Delta Function  Rassolov/Chipman')
 1400 FORMAT(I4,1X,A4,2X,4(F12.5,5X))
 1500 FORMAT(/,' ELECTRIC FIELD GRADIENT COMPONENTS in atomic units',
     $         ' (Tensor Representation)',
     $        /,'  Atom       3XX-RR     3YY-RR     3ZZ-RR',
     $          '      3XY        3XZ        3YZ')
 1600 FORMAT(I4,1X,A4,1X,6(F10.5,1X))
 1800 FORMAT(/,' INVARIANT EIGENVALUES OF THE EFG TENSOR',
     $       /,'  Atom',15X,'Eigenvalues')
 2000 FORMAT(2X,'***ERROR*** Unable to diagonalize EFG tensor',
     $          ' for atom: ',I4)
c
      END
