C***********************************************************************
C                    PQS VERSION OF NBO
C                (FOR PQS98  by JON BAKER)
C
C       Jon Bohmann  UW MADISON   1/96
C       Jon Bohmann  at Univ. of Ark.  8/97
C  DRIVER ROUTINES:
C
C      SUBROUTINE RUNNBO
C      SUBROUTINE FEAOIN(CORE,ICORE,MEMORY1,NBOOPT)
C      NO DELSCF
C      NO CHKNBO(T)
C
C***********************************************************************
C real code can be written for the following  SR's in the future, but
C their presence is not critical for basic NBO analysis
C***********************************************************************
       SUBROUTINE DELSCF
C***********************************************************************
       RETURN
       END
C***********************************************************************
       SUBROUTINE CHKNBO(Idum)
C***********************************************************************
       RETURN
       END
C***********************************************************************
       SUBROUTINE RUNNBO(inp)
C***********************************************************************
C 12-Jan-96  JAB First TEXAS version
C 19-Jan-96  Fixed erroneous parameters in calls to SR GETMEM
C  6-Aug-97  Changes for Texas 97
C             and SR RETMEM
C------------------------------

      use memory

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 CINPUT,jobname,nbostr
      INTEGER NLINES,N,MEMORY1
      LOGICAL license, chklicense
      DIMENSION NBOOPT(10)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +           LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +           LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
C     TEXAS COMMON BLOCKS
c     common /mmanag/lcore
c     common /big/bl(1000)
c
      call getival('lcore',lcore)
      CALL GETIVAL('iout',ICOND)
C Input file for TEXAS is unit INP.
C "Short" output is unit ICOND.
      LFNIN=INP
      LFNPR=ICOND
      LTEMP = 55
c
c -- NBO INPUT SECTION --
c    Two methods
c  (1) Read input file until find ENDNBO (indicating end of NBO input)
c  (2) If ENDNBO not found, read to end of file
c    In either case, the NBO section of the input is written to a
c    temporary file (file LFNNBO = 55)
c
c -- read in (again) initial NBO line
      READ(LFNIN,*)
C
C Open unit LTEMP as scratch file
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
      OPEN(UNIT=LTEMP,FILE=jobname(1:lenJ)//'.nbo',FORM ='FORMATTED')
c
C NLINES counts the number of records (cards) that we copy from input
C to unit LTEMP
      NLINES=0
C Read input lines from unit INP = 30 and write them to unit LTEMP
C
   20 CONTINUE
      READ(INP,1000,END=30) CINPUT
 1000 FORMAT(A80)
      If(CINPUT(1:6).EQ.'ENDNBO') GO TO 40
c
      WRITE(LTEMP,1000) CINPUT
      NLINES = NLINES+1
      GO TO 20
c
 30   CONTINUE
C Now, backspace input file so that TEXAS will use remaining records
      DO I=1,NLINES+1
        BACKSPACE(INP)
      END DO
c
 40   CONTINUE
c
c ............................................................
c  LICENSE CHECKING      ! JB    Oct 2006
c  The NBO module is now optional as we have to pay Frank
c  Weinhold a license fee and we can no longer absorb this
c  for software-only sales
c ............................................................
c
      call getchval('nbostr',nbostr)  ! NBO license
      license = chklicense(nbostr)
c
      If(.NOT.license) Then
        WRITE(6,900)
  900   FORMAT(/,' Sorry! You do NOT have a license for NBO',/,
     $           ' ** Contact PQS to enable this option **',/)
        CLOSE (UNIT=LTEMP,STATUS='DELETE')
        RETURN
      EndIf
c
      LFNIN=LTEMP           ! change NBO input definition
C End of keylist input file creation
      WRITE(LTEMP,*)'$NBO $END'
      REWIND(LTEMP)
C
C
C Set NBO options.
C
      NBOOPT(1)  =  0
C NBOOPT(1)=1, No Keylist
C NBOOPT(1)=0 Make NBO look for a keylist
      NBOOPT(2)  =  0
      NBOOPT(3)  =  0
      NBOOPT(4)  =  0
      NBOOPT(5)  =  0
      NBOOPT(6)  =  0
      NBOOPT(7)  =  0
      NBOOPT(8)  =  0
      NBOOPT(9)  =  0
      NBOOPT(10) =  6
C NBOOPT(10)=6   Make NBO act like Gamess version
C
C For now, reserve a constant amount of memory for NBO from the
C TEXAS program memory array.
        call getmem(0,lastx)
        call retmem(1)
C  Perform the NPA/NBO/NLMO analyses.
      CALL NBO(BL(LASTX),LCORE-LASTX,NBOOPT)
C      IF(NBOOPT(10).LT.0) RETURN
C [Here's an example of how energetic analysis might be performed
C in the future:]
C Perform the energetic analysis.
C      NBOOPT(7) = 0
C   10 NBOOPT(1) = 2
C      CALL NBOEAN(bl(LAST+N),MEMORY1,NBOOPT,IDONE)
C      IF(IDONE.NE.0) GOTO 20
C     CALL DELSCF
C     NBOOPT(1) = 3
C     CALL NBOEAN(CORE(ICUR+1),MEMORY1,NBOOPT,IDONE)
C     GOTO 10
C  FREE THE LAST RESERVED MEMORY BLOCK
      CLOSE(UNIT=LTEMP,STATUS='DELETE')
      RETURN
      END

C***********************************************************************
      SUBROUTINE FEAOIN(CORE,ICORE,MEMORY1,NBOOPT)

      use memory

c     COMMON/intbl/MAXSH,INX(1000)

      CALL GETIVAL('ictr',ICTR)
      CALL FEAOTX(CORE,ICORE,MEMORY1,NBOOPT,bl(ICTR))
      RETURN
      END

C***********************************************************************
      SUBROUTINE FEAOTX(CORE,ICORE,MEMORY1,NBOOPT,INX)
C***********************************************************************
C 12-Jan-96  JAB FIRST TEXAS VERSION
C  6-Feb-96  Fixed label ordering of cartesian F--functions
C  7-Aug-97  Name changed to FEAOTX to work with TEXAS storage
C------------------------------

      use memory

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 TITLE,jobname
      DIMENSION CORE(MEMORY1),ICORE(MEMORY1),NBOOPT(10)
      Logical found,toobig
C ----------------------------------------------------------------------
C
C This routine fetches AO basis set information from the TEXAS
C binary files and COMMON blocks
C and stores it in the NBO common blocks and direct access file
C (DAF) for use by the NBO analysis.
C
C  NBO Common blocks
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,KOPT,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBAO/LCTR(MAXBAS),LANG(MAXBAS)
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)

C    TEXAS COMMON BLOCKS
c      COMMON /BIG/BL(1000)
       DIMENSION INX(12,*)
C The following DIMENSION LABELS statement
C allocates memory for angular momentum labels of contracted functions.
C When spherical f functions or g functions are added, the size of the
C labels array must be changed accordingly
C  LABELS array contains NBO labels for the atomic orbitals
C  the ordering is specific to texas93
C In the future, this array should be expanded to provide for
C spherical F functions and the cartesian G functions, once their
C order in TEXAS is known.

      DIMENSION LABELS(25)
      DATA LABELS /
C
C           s
C          ---
     +      1,
C
C           px    py    pz
C          ---   ---   ---
     +      101,  102,  103,
C
C          dxx   dyy   dzz   dxy   dxz   dyz
C          ---   ---   ---   ---   ---   ---
     +      201,  204,  206,  202,  203,  205,
C        Ordering for TEXAS spherical D functions:
     +     255,  254,  251,  252,  253,
C
C        Ordering for TEXAS cartesian F functions:
     +     301,  302,  303,  304,  305,  306,  307,  308,  309,310/
C
C*****
       DATA ZERO/0.0d00/
      character*8 atsymb         ! these are needed for setting
      real*8 xname               ! the atomic numbers
      Equivalence (xname,atsymb) ! correctly
C
C Get the job title of the calculation
C To get the calculation title right, write the title to a scratch file,
C and read it back in using the format for NBO.
C First, open scratch file

      OPEN(UNIT=56,FORM='FORMATTED',STATUS='SCRATCH')
C Write title in TEXAS format
C  First check if a title exists
      CALL TSTCHVAL('titl',iyes)
      if(iyes.eq.0) then
        TITLE='PQS CALCULATION'
      else
        CALL GETCHVAL('titl',TITLE)
      END IF
      WRITE(56,'(A)') TITLE
      REWIND(56)
C Read title in NBO format
      READ (56,1100) (CORE(I),I=1,8)
 1100 FORMAT(8(A8))
      CLOSE (56,STATUS='DELETE')
      NFILE = 2
C Write calculation title to NBO daf
      CALL NBWRIT(CORE,8,NFILE)
C    NATOMS      Number of atomic centers
C    NDIM        Dimension of matrices (overlap and density)
C
C    NBAS        Number of basis functions (.le.NDIM)

      call getival('na',NATOMS)
      call getival('ndum1',Ndum1)
      call getival('ndum2',Ndum2)
      call getival('ncf',NBAS)
      NDIM  = NBAS
c
c -- do not do NBO analysis on charged dummy atoms
      NATOM = NATOMS-Ndum1-Ndum2      ! number of real atoms
      NATOMS = NATOMS-Ndum1           ! real atoms + uncharged dummies
c
c -- NBO currently has a 200 atom and 2000 basis function limit
c -- Check that these limits are not exceeded
      toobig = .False.
      If(NATOMS.GT.MAXATM) Then
       WRITE(6,1200) NATOMS,MAXATM
       toobig = .True.
      EndIf
      If(NBAS.GT.MAXBAS) Then
       WRITE(6,1210) NBAS,MAXBAS
       toobig = .True.
      EndIf
 1200 FORMAT(/,2X,'***ERROR*** Too Many atoms for NBO analysis',/,
     $          '     System has ',I4,' atoms;  Maximum Allowed is ',I4)
 1210 FORMAT(/,2X,'***ERROR*** Too Many basis functions for',
     $            ' NBO analysis',/,'   System has ',I5,
     $            ' functions;  Maximum Allowed is ',I5)
c
      If(toobig) Call nerror(1,'NBO Analysis','  System too big',0,0)
C
C Inititialize various NBO options depending upon the wavefunction
C type and basis set type.
C For TEXAS, now use only HF SCF wavefunction
C Support for UHF and MCSCF
C should be added in the future. (JAB)
C  Obtain the following information:
C
C    ROHF        =.TRUE. If RHF open shell wavefunction
C                =.FALSE. otherwise
C
C    UHF         =.TRUE. If UHF wavefunction
C                =.FALSE. otherwise
C
C    AUHF        =.TRUE. If spin-annihilated UHF wavefunction
C                =.FALSE. otherwise
C
C    CI          =.TRUE. If CI wavefunction
C                =.FALSE. otherwise
C
C    OPEN        =.TRUE. If open shell wavefunction
C                =.FALSE. otherwise
C
C    COMPLX      =.TRUE. If complex wavefunction
C                =.FALSE. otherwise
C                (Note: The program is not capable of handling this.)
C
C    IPSEUD      Set to one if pseudopotentials are used.
C
C    IWCUBF      This pertains only basis sets with F functions.
C        Note: as of Jan 96 neither "STANDARD" or "CUBIC" pure F
C              functions seems appropriate for TEXAS 93.
C                If cartesian F functions are input, set IWCUBF to:
C                    0,  if these are to be transformed to the
C                        standard set of pure F functions
C                    1,  if these are to be transformed to the
C                        cubic set of pure F functions
C
C                If pure F functions are input, set to IWCUBF to:
C                    0,  if these are standard F functions
C                    1,  if these are cubic F functions
      COMPLX = .FALSE.
      IWCUBF =  0
      ORTHO  = .FALSE.
C TEXAS uses the charge--bond order type density matrix
C SO we must have IWDM=1
      IWDM   =  1
      ROHF   = .FALSE.
      UHF    = .FALSE.
      CI     = .FALSE.
      OPEN   = .FALSE.
      MCSCF  = .FALSE.
      AUHF   = .FALSE.
C
      MUNIT = 0
      IF(NBOOPT(2).EQ.0) WRITE(LFNPR,1300)
 1300 FORMAT(/1X,'Analyzing the SCF density')
C
C  Store NATOMS, NDIM, NBAS, MUNIT, wavefunction flags, ISWEAN:
C
      ICORE(1)  = NATOMS
      ICORE(2)  = NDIM
      ICORE(3)  = NBAS
      ICORE(4)  = MUNIT
      ICORE(5)  = 0
      If(ROHF)  ICORE(5)  = 1
      ICORE(6)  = 0
      If(UHF)   ICORE(6)  = 1
      ICORE(7)  = 0
      If(CI)    ICORE(7)  = 1
      ICORE(8)  = 0
      If(OPEN)  ICORE(8)  = 1
      ICORE(9)  = 0
      If(MCSCF) ICORE(9)  = 1
      ICORE(10) = 0
      If(AUHF)  ICORE(10) = 1
      ICORE(11) = 0
      If(ORTHO) ICORE(11) = 1
      ICORE(12) = 1
      NFILE = 3
      CALL NBWRIT(ICORE,12,NFILE)
C
C  Store the total energy on the NBO DAF:
C    Energy stored in  ETOT in TEXAS common/ener/
C When DELSCF exists, CORE(1) below should be the deletion energy
C and CORE(2) should be the total energy.
C? Where is ETOT in TEXAS 97  ?
      CORE(1) = ZERO ! ETOT
      CORE(2) = ZERO ! ETOT
      NFILE = 8
      CALL NBWRIT(CORE,2,NFILE)
C Prepare
C      IATNO(I), I=1,NATOMS
C                     List of atomic numbers
C
C      LCTR(I),I=1,NBAS
C                     List of atomic centers of the basis functions
C                     (LCTR(3)=2 if basis function 3 is on atom 2
C
C      LANG(I) I=1,NBAS
C                     List of angular symmetry information for the
C                     basis functions
C
C Get charge (atomic number) and x,y,z coordinates of atoms.
C Nuclear information is stored in TEXAS linear array bl(inuc)
C with the charge,x-coordinate, y--coordinate, z-coordinate and name
C (N=   card from input) for  each nucleus stored consecutively.
C? J=INUC is Charge of first atom
C? J=INUC+1 is X coordinate of first atom,etc.
      CALL GETRVAL('angs',ANG)
      CALL GETIVAL('inuc',INUC)
c     PRINT*,'INUC=', INUC
c     PRINT*,'ICTR=', ICTR
      J = INUC
      I = 0
c
c  set pseudopotential flag
c
      call getival('npsp',npsp)
      if( npsp .gt. 0 ) ipseud = 1
c
c -- loop over REAL atoms first
      DO 110 IAT = 1,NATOM
C  Store atomic numbers in
C  IATNO and nuclear charges in IZNUC.
        xname = bl(j+4)
        call nugrep( atsymb(1:2), iatn )
        IATNO(IAT)=iatn
        IZNUC(IAT)=BL(J)
c -- get X,Y,Z coordinates
        CORE(I+1) = BL(J+1)/ANG
        CORE(I+2) = BL(J+2)/ANG
        CORE(I+3) = BL(J+3)/ANG
        I=I+3
        J=J+5
  110 CONTINUE
c
c -- now loop over uncharged dummies (if any) IGNORING charged dummies
      If(Ndum2.GT.0) Then
        J = J+5*Ndum1
        DO 111 IAT=NATOM+1,NATOMS
        IATNO(IAT)=0
        IZNUC(IAT)=0
c -- get X,Y,Z coordinates
        CORE(I+1) = BL(J+1)/ANG
        CORE(I+2) = BL(J+2)/ANG
        CORE(I+3) = BL(J+3)/ANG
        I=I+3
        J=J+5
  111   CONTINUE
      EndIf
c
      NFILE = 9
      CALL NBWRIT(CORE,3*NATOMS,NFILE)
C Figure out LCTR and LANG arrays
C Get number of contracted shells in basis set
      CALL GETIVAL('ncs',ncs)
C ncs is the number of "shells" in Texas which MAY include
C generally contracted shells, each one of which increases
C the number of REAL shells (for NBO) by 1
      II = 0
      MM = 0
C loop over number of contracted shells (ncs)
C at the same time, determine the number of NBO shells (NSHELL)
C and NBO primitive exponents (NEXP)
      NSHELL = 0
      NEXP = 0
      DO 30 K = 1,ncs
      ngr = inx(4,k)
      NSHELL = NSHELL + ngr + 1
      NEXP = NEXP + (INX(5,K)-INX(1,K))*(ngr+1)
cc        IF (INX(4,K).NE.0) THEN
cc          WRITE(LFNPR,*)'Generalized contractions cannot be used by NBO'
cc          STOP
cc        ENDIF
      do 29 igr=0,ngr
C find atom on which this shell is centered
        IATOM = INX(2,K)
c       PRINT*, 'IATOM=', IATOM
C find shell size (gives angular momentum type of basis functions)
        MAX  =  INX(3,K)
c       PRINT*, 'MAX =' , MAX
C   each shell is composed of MAX contracted basis functions
C     shell size = 1   s
C     shell size = 3   px,py,pz
C     shell size = 4   s, px ,py ,pz
C     shell size = 6   cartesian D's
C     shell size = 5   pure D's
C     shell size = 7   pure F's
C     shell size = 10  cartesian F's
C     shell size = 15  cartesian G's
C In the following, JJ is the starting index to use for the label array
C JJ is determined by the value of MAX
C S functions
        IF(MAX.EQ.1) THEN
          JJ = 1
C L type functions (S+P)
        ELSE IF (MAX.EQ.4) THEN
          JJ = 1
C  P functions
        ELSE IF (MAX.EQ.3) THEN
          JJ = 2
C pure D
        ELSE IF (MAX.EQ.5) THEN
          JJ = 11
C cartesian D
        ELSE IF (MAX.EQ.6) THEN
          JJ = 5
C pure F
        ELSE IF (MAX.EQ.7) THEN
C When a set of LABELS array values for the pure F functions
C are found, a value of JJ may be determined and the following
C stop statement may be removed.
          call nerror(1,'FEAOIN',
     $     ' Pure F functions not yet supported.',0,0)
C cartesian F's
        ELSE IF (MAX.EQ.10) THEN
          JJ = 16
C cartesian G's
C When a set of LABELS array values for the cartesian G functions
C are found, a value of JJ may be determined and the following stop
C statement may be removed.
        ELSE IF (MAX.EQ.15) THEN
          call nerror(1,'FEAOIN',
     $     'G functions not yet supported.',0,0)
        ELSE
          call nerror(1,'FEAOIN','erroneus max value',0,0)
        ENDIF
C Assign angular momentum code to each contracted basis function

        DO 20 J = 1,MAX
          II = II + 1
          LCTR(II) = IATOM
          LANG(II) = LABELS(JJ)
          JJ = JJ + 1
   20   CONTINUE
c
   29 continue
   30 CONTINUE
C
C  Store IATNO, IZNUC, LCTR, and LANG on NBO DAF:
C
      II = 0
      DO 70 I = 1,NATOMS
        II = II + 1
        ICORE(II) = IATNO(I)
   70 CONTINUE
      DO 80 I = 1,NATOMS
        II = II + 1
        ICORE(II) = IZNUC(I)
   80 CONTINUE
      DO 90 I = 1,NBAS
        II = II + 1
        ICORE(II) = LCTR(I)
   90 CONTINUE
      DO 95 I = 1,NBAS
        II = II + 1
        ICORE(II) = LANG(I)
   95 CONTINUE
      NFILE = 4
      CALL NBWRIT(ICORE,2*NATOMS+2*NBAS,NFILE)

C -------------------------------------------------------------
C Get information for DAF from TEXAS binary files.
      NP1 = 1
      NP4 = 4
      L2 = NDIM*(NDIM+1)/2
C
C  first get overlap, h0 and KE matrices (always needed)
C
      CALL REA (CORE,L2,NP1,'s matrix')
      NFILE = 10
      CALL NBWRIT(CORE,L2,NFILE)
c
      CALL REA(CORE,L2,NP1,'kinetic ')
      NFILE = 18
      CALL NBWRIT(CORE,L2,NFILE)
C
C subtract ke from h0 to get nuclear attraction integrals
C
      I1 = 1
      I2 = I1 + L2
      CALL REA(CORE(I2),L2,NP1,'h0-mtx ')
C
C write out nuclear attraction integrals
C
      DO I = 1,L2
        CORE(I2+I-1) = CORE(I2+I-1) - CORE(I1+I-1)
      ENDDO
      NFILE = 19
      CALL NBWRIT(CORE(I2),L2,NFILE)
C
C
C  Closed/Open-Shell tested for by existence of nonzero <s2>
C  in depository
C
      call getrval('<s2>',s2)
c
      IF(s2.NE.ZERO) THEN
C
C  Open Shell
C
        UHF  = .TRUE.
        OPEN = .TRUE.
c
c -- density matrices
        CALL REA(CORE,L2,NP4,'dena_uhf')
        NFILE = 20
        CALL NBWRIT(CORE,L2,NFILE)
        CALL REA(CORE,L2,NP4,'denb_uhf')
        NFILE = 21
        CALL NBWRIT(CORE,L2,NFILE)
c
c -- Fock matrices
        IF(IWFOCK.NE.0) THEN
          CALL REA(CORE,L2,NP4,'foca_uhf')
          NFILE = 30
          CALL NBWRIT(CORE,L2,NFILE)
          CALL REA(CORE,L2,NP4,'focb_uhf')
          NFILE = 31
          CALL NBWRIT(CORE,L2,NFILE)
        ENDIF
c
c -- Molecular Orbitals
        L3 = NDIM*NDIM
        CALL REA(CORE,L3,NP4,'evea_uhf')
        NFILE = 40
        CALL NBWRIT(CORE,L3,NFILE)
        CALL REA(CORE,L3,NP4,'eveb_uhf')
        NFILE = 41
        CALL NBWRIT(CORE,L3,NFILE)
cc
      ELSE
C
C  Closed Shell
C
c -- density matrix
        CALL REA(CORE,L2,NP4,'den0_rhf')
        NFILE = 20
        CALL NBWRIT(CORE,L2,NFILE)
c
c -- Fock matrix
        IF(IWFOCK.NE.0) THEN
          CALL REA(CORE,L2,NP4,'fock_rhf')
          NFILE = 30
          CALL NBWRIT(CORE,L2,NFILE)
        END IF
c
c -- Molecular Orbitals
        L3 = NDIM*NDIM
        CALL REA(CORE,L3,NP4,'evec_rhf')
        NFILE = 40
        CALL NBWRIT(CORE,L3,NFILE)
      ENDIF
C
C  Store the dipole integrals on the NBODAF:
C
C X Dipole
      CALL REA(CORE,L2,NP1,'aoX    ')
      DO I = 1,L2
        CORE(I)=CORE(I)/ANG
      END DO
      NFILE = 50
      CALL NBWRIT(CORE,L2,NFILE)
C Y Dipole
      CALL REA(CORE,L2,NP1,'aoY    ')
      DO  I = 1,L2
        CORE(I)=CORE(I)/ANG
      END DO
      NFILE = 51
      CALL NBWRIT(CORE,L2,NFILE)
C Z Dipole
      CALL REA(CORE,L2,NP1,'aoZ    ')
      DO  I = 1,L2
        CORE(I)=CORE(I)/ANG
      END DO
      NFILE = 52
      CALL NBWRIT(CORE,L2,NFILE)

ccc      CALL GETIVAL ('nsh',nsh)
      CALL GETIVAL ('ibas',IBAS)
c
      LEN=2+3*NSHELL+6*NEXP
C
C Later we will need NEXP  more memory for sorting gaussians by
C S,P,D,F,G type.  We will start after the LEN locations needed
C here (we add +1 later)
C
      DO 160 I = 1,LEN
        CORE(I) = ZERO
  160 CONTINUE
      ICORE(1) = NSHELL
      ICORE(2) = NEXP
C
C  -- the number of components in the Ith shell:
C
      II = 2
      DO 170 I = 1,ncs
       do igr=0,inx(4,i)
        II = II + 1
        ICORE(II) = INX(3,I)
       enddo
  170 CONTINUE

C
C  NPRIM   -- the number of gaussian primitives in the Ith shell:
C
      DO 180 I = 1,ncs
       do igr=0,inx(4,i)
        II = II + 1
        ICORE(II) = INX(5,i)-INX(1,i)
       enddo
  180 CONTINUE

C
C  NPTR    -- pointer for the Ith shell into the gaussian parameters,
C             EXP, CS, CP, etc.:
C
      ngroff=inx(1,1)     ! offset due to expansion of contracted shells
      DO 190 I = 1,ncs
       ngr = inx(4,i)
       nump = INX(5,i)-INX(1,i)
       do igr=0,ngr
        II = II + 1
        ICORE(II) = ngroff + 1
        ngroff=ngroff+nump
       enddo
  190 CONTINUE
C
C  get exponents and contraction coefficients from BASDAT and
C  load them into CORE for NBO analysis
C
      j=0
      DO 200 I = 1,ncs
      IType = inx(12,i)             ! function type
      nump = INX(5,i)-INX(1,i)      ! number of primitives in function
      ngr = inx(4,i)                ! number of general contraction
      j = 13*inx(1,i)
      jstart = j
      do igr=0,ngr
        j = jstart
        do k=1,nump
          II = II + 1
          CORE(II) = bl(ibas+j)     ! load exponents
          kk = ibas+j+1+igr         ! index into BASDAT for coeff.
          If(IType.EQ.1) Then
            CORE(II+NEXP) = bl(kk)
          Else If(IType.EQ.2) Then
            CORE(II+2*NEXP) = bl(kk)
          Else If(IType.EQ.3) Then
            CORE(II+NEXP) = bl(kk)
            CORE(II+2*NEXP) = bl(kk+1)
          Else If(IType.EQ.4.OR.IType.EQ.5) Then
            CORE(II+3*NEXP) = bl(kk)
          Else If(IType.EQ.6.OR.IType.EQ.7) Then
            CORE(II+4*NEXP) = bl(kk)
          Else If(IType.EQ.8) Then
            CORE(II+5*NEXP) = bl(kk)
          EndIf
          j=j+13
        enddo
      enddo
  200 CONTINUE
C
C  write to NBO direct access file
C
      NFILE = 5
      CALL NBWRIT(CORE,LEN,NFILE)
      RETURN
      END


C***********************************************************************UTILTY
C                                                                       UTILTY
C  DEFAULT UTILITY ROUTINES (INCLUDE SOME OPERATING SYSTEM DEPENDENCE): UTILTY
C                                                                       UTILTY
C      SUBROUTINE DAFILE(NAME,LEN,INBO,ISINGL,LENREC,ERROR)             UTILTY
C      SUBROUTINE DACLOS(INBO,DEL)                                      UTILTY
C      SUBROUTINE DEBYTE(I,IBYTE)                                       UTILTY
C      SUBROUTINE SQFILE(NEW,ERROR)                                     UTILTY
C                                                                       UTILTY
C***********************************************************************UTILTY
      SUBROUTINE DAFILE(NAME,LEN,INBO,ISINGL,LENREC,ERROR)              UTILTY
C***********************************************************************UTILTY
C 22-Jan-93  EDG  New subroutine                                        UTILTY
C------------------------------                                         UTILTY

      use kinds

      IMPLICIT REAL*8 (A-H,O-Z)                                         UTILTY
      LOGICAL ERROR                                                     UTILTY
      CHARACTER*256 NAME                                                UTILTY
      CHARACTER*3 EXT                                                   UTILTY
C                                                                       UTILTY
C  The following parameter is the length of the IDXNBO array in         UTILTY
C  COMMON/NBONAV/:                                                      UTILTY
C                                                                       UTILTY
      PARAMETER(LENGTH = 256)                                           UTILTY
C                                                                       UTILTY
      PARAMETER (MAXFIL = 40)                                           UTILTY
      COMMON/NBNAME/FILENM,NFILE,IFILE(MAXFIL)                          UTILTY
      CHARACTER*256 FILENM                                              UTILTY
C                                                                       UTILTY
      ISINGL = 0                                                        UTILTY
      LENREC = 0                                                        UTILTY
C                                                                       UTILTY
C  Are we working on a 32 (ISINGL=2) or 64 (ISINGL=1) bit machine?      UTILTY
C                                                                       UTILTY
      ISINGL = intsize
*      ISINGL = 1                                                        UNICOS
*     ISINGL = 2                                                        UNIX
*      ISINGL = 2                                                        VMS
C                                                                       UTILTY
C  Determine an appropriate record length for the NBO DAF:              UTILTY
C                                                                       UTILTY
      LENREC = 8 / intsize * LENGTH
*     LENREC = 8 * LENGTH                                               UNICOS
*     LENREC = 4 * LENGTH                                               UNIX
*      LENREC = LENGTH                                                   VMS
C                                                                       UTILTY
C  Filename for the DAF:                                                UTILTY
C                                                                       UTILTY
      NAME = FILENM                                                     UTILTY
      LEN = MIN(LENNB(NAME),252)                                        UTILTY
      IF(INBO.LT.10) THEN                                               UTILTY
        L = 1                                                           UTILTY
        WRITE(EXT,'(I1)') INBO                                          UTILTY
      ELSE IF(INBO.LT.100) THEN                                         UTILTY
        L = 2                                                           UTILTY
        WRITE(EXT,'(I2)') INBO                                          UTILTY
      ELSE                                                              UTILTY
        L = 3                                                           UTILTY
        WRITE(EXT,'(I3)') INBO                                          UTILTY
      END IF                                                            UTILTY
      NAME = NAME(1:LEN) // '.' // EXT(1:L)                             UTILTY
      LEN = LENNB(NAME)                                                 UTILTY
C                                                                       UTILTY
      ERROR = .FALSE.                                                   UTILTY
      IF(ISINGL.EQ.0) ERROR = .TRUE.                                    UTILTY
      IF(LENREC.EQ.0) ERROR = .TRUE.                                    UTILTY
      RETURN                                                            UTILTY
      END                                                               UTILTY
C***********************************************************************UTILTY
      SUBROUTINE DACLOS(INBO,DEL)                                       UTILTY
C***********************************************************************UTILTY
C 14-Feb-93  EDG  New subroutine                                        UTILTY
C------------------------------                                         UTILTY
      IMPLICIT REAL*8 (A-H,O-Z)                                         UTILTY
      LOGICAL DEL                                                       UTILTY
C                                                                       UTILTY
      IF(DEL) THEN                                                      UTILTY
        CLOSE(UNIT=INBO, STATUS='DELETE')                               UTILTY
      ELSE                                                              UTILTY
        CLOSE(UNIT=INBO, STATUS='KEEP')                                 UTILTY
      END IF                                                            UTILTY
      RETURN                                                            UTILTY
      END                                                               UTILTY
C***********************************************************************UTILTY
      SUBROUTINE DEBYTE(I,IBYTE)                                        UTILTY
C***********************************************************************UTILTY
      IMPLICIT REAL*8 (A-H,O-Z)                                         UTILTY
      DIMENSION IBYTE(4),KB(4)                                          UTILTY
C                                                                       UTILTY
      SAVE KB,KPAD,KSW,KTMP                                             UTILTY
C                                                                       UTILTY
      DATA KSW/0/                                                       UTILTY
      DATA KTMP/4HABCD/                                                 UTILTY
C                                                                       UTILTY
C  Extract four Hollerith characters from I, store in IBTYE:            UTILTY
C                                                                       UTILTY
C  If this is the first time that this routine is called, determine     UTILTY
C  in which bytes of an integer word the Hollerith characters reside:   UTILTY
C                                                                       UTILTY
      IF(KSW.EQ.0) THEN                                                 UTILTY
        KSW   = 1                                                       UTILTY
        DO 10 K = 1,4                                                   UTILTY
          KB(K) = 0                                                     UTILTY
   10   CONTINUE                                                        UTILTY
        KBYTE = 0                                                       UTILTY
   20   KBYTE = KBYTE + 1                                               UTILTY
          IF(KBYTE.GT.8)
     $    call nerror(1,'DEBYTE','limited to integer*8',0,0)
          KTEST = MOD(KTMP,256)                                         UTILTY
          IF(KTEST.EQ.65) KB(1) = KBYTE                                 UTILTY
          IF(KTEST.EQ.66) KB(2) = KBYTE                                 UTILTY
          IF(KTEST.EQ.67) KB(3) = KBYTE                                 UTILTY
          IF(KTEST.EQ.68) KB(4) = KBYTE                                 UTILTY
          KTMP = KTMP/256                                               UTILTY
        IF(KTMP.NE.0) GOTO 20                                           UTILTY
        DO 30 K = 1,4                                                   UTILTY
          IF(KB(K).EQ.0)                                                UTILTY
     $    call nerror(1,'DEBYTE','error',0,0)
   30   CONTINUE                                                        UTILTY
C                                                                       UTILTY
C  Determine the bit padding:                                           UTILTY
C                                                                       UTILTY
        KPAD = 0                                                        UTILTY
        KMLT = 1                                                        UTILTY
        DO 40 K = 1,KBYTE                                               UTILTY
          IF(K.NE.KB(1)) KPAD = KPAD + 32 * KMLT                        UTILTY
          IF(K.NE.KBYTE) KMLT = KMLT * 256                              UTILTY
   40   CONTINUE                                                        UTILTY
C                                                                       UTILTY
        DO 60 K = 1,4                                                   UTILTY
          KMAX  = KB(K) - 1                                             UTILTY
          KB(K) = 1                                                     UTILTY
          DO 50 L = 1,KMAX                                              UTILTY
            KB(K) = KB(K) * 256                                         UTILTY
   50     CONTINUE                                                      UTILTY
   60   CONTINUE                                                        UTILTY
      END IF                                                            UTILTY
C                                                                       UTILTY
C  Extract four Hollerith characters from I:                            UTILTY
C                                                                       UTILTY
      DO 100 K = 1,4                                                    UTILTY
        IBYTE(K) = MOD(I/KB(K),256)*KB(1) + KPAD                        UTILTY
  100 CONTINUE                                                          UTILTY
      RETURN                                                            UTILTY
      END                                                               UTILTY
C***********************************************************************UTILTY
      SUBROUTINE SQFILE(NEW,ERROR)                                      UTILTY
C***********************************************************************UTILTY
C 10-Feb-93  EDG  New subroutine                                        UTILTY
C------------------------------                                         UTILTY
      IMPLICIT REAL*8 (A-H,O-Z)                                         UTILTY
      LOGICAL NEW,ERROR,THERE                                           UTILTY
      CHARACTER*3 CNT,EXT                                               UTILTY
      CHARACTER*256 TEMP                                                UTILTY
C                                                                       UTILTY
      PARAMETER (MAXFIL = 40)                                           UTILTY
C
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,KOPT,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBNAME/FILENM,NFILE,IFILE(MAXFIL)                          UTILTY
      CHARACTER*256 FILENM                                              UTILTY
C                                                                       UTILTY
      DATA IREAD/4HREAD/                                                UTILTY
C                                                                       UTILTY
C  Select an alternate filename if this one is not acceptable:          UTILTY
C                                                                       UTILTY
      IO = IOINQR(IWPNAO)                                               UTILTY
      JO = IOINQR(IWTNAO)                                               UTILTY
      KO = IOINQR(IWTNAB)                                               UTILTY
      LO = IOINQR(IWNBBP)                                               UTILTY
      IF(NEW.AND.IO.NE.IREAD.AND.JO.NE.IREAD.AND.                       UTILTY
     +           KO.NE.IREAD.AND.LO.NE.IREAD) THEN                      UTILTY
        TEMP = FILENM                                                   UTILTY
        LEN = MIN(LENNB(TEMP),252)                                      UTILTY
        DO 20 I = 0,999                                                 UTILTY
          IF(I.EQ.0) THEN                                               UTILTY
            L1 = 0                                                      UTILTY
          ELSE IF(I.LT.10) THEN                                         UTILTY
            L1 = 1                                                      UTILTY
            WRITE(CNT,'(I1)') I                                         UTILTY
          ELSE IF(I.LT.100) THEN                                        UTILTY
            L1 = 2                                                      UTILTY
            WRITE(CNT,'(I2)') I                                         UTILTY
          ELSE                                                          UTILTY
            L1 = 3                                                      UTILTY
            WRITE(CNT,'(I3)') I                                         UTILTY
          END IF                                                        UTILTY
          IF(L1.GT.0) THEN                                              UTILTY
            TEMP = TEMP(1:LEN) // CNT(1:L1)                             UTILTY
          ELSE                                                          UTILTY
            TEMP = TEMP(1:LEN)                                          UTILTY
          END IF                                                        UTILTY
          LENGTH = MIN(LENNB(TEMP),252)                                 UTILTY
C                                                                       UTILTY
C  First check the DAF:                                                 UTILTY
C                                                                       UTILTY
          IF(ABS(LFNDAF).LT.10) THEN                                    UTILTY
            L2 = 1                                                      UTILTY
            WRITE(EXT,'(I1)') ABS(LFNDAF)                               UTILTY
          ELSE IF(ABS(LFNDAF).LT.100) THEN                              UTILTY
            L2 = 2                                                      UTILTY
            WRITE(EXT,'(I2)') ABS(LFNDAF)                               UTILTY
          ELSE IF(ABS(LFNDAF).LT.1000) THEN                             UTILTY
            L2 = 3                                                      UTILTY
            WRITE(EXT,'(I3)') ABS(LFNDAF)                               UTILTY
          END IF                                                        UTILTY
          TEMP = TEMP(1:LENGTH) // '.' // EXT(1:L2)                     UTILTY
          INQUIRE(FILE=TEMP, EXIST=THERE)                               UTILTY
          IF(THERE) GO TO 20                                            UTILTY
C                                                                       UTILTY
C  Now check the rest:                                                  UTILTY
C                                                                       UTILTY
          DO 10 J = 1,NFILE                                             UTILTY
            IF(ABS(IFILE(J)).LT.10) THEN                                UTILTY
              L2 = 1                                                    UTILTY
              WRITE(EXT,'(I1)') ABS(IFILE(J))                           UTILTY
            ELSE IF(ABS(IFILE(J)).LT.100) THEN                          UTILTY
              L2 = 2                                                    UTILTY
              WRITE(EXT,'(I2)') ABS(IFILE(J))                           UTILTY
            ELSE IF(ABS(IFILE(J)).LT.1000) THEN                         UTILTY
              L2 = 3                                                    UTILTY
              WRITE(EXT,'(I3)') ABS(IFILE(J))                           UTILTY
            END IF                                                      UTILTY
            TEMP = TEMP(1:LENGTH) // '.' // EXT(1:L2)                   UTILTY
            INQUIRE(FILE=TEMP, EXIST=THERE)                             UTILTY
            IF(THERE) GO TO 20                                          UTILTY
   10     CONTINUE                                                      UTILTY
          GO TO 30                                                      UTILTY
   20   CONTINUE                                                        UTILTY
        WRITE(LFNPR,900)                                                UTILTY
        ERROR = .TRUE.                                                  UTILTY
        RETURN                                                          UTILTY
C                                                                       UTILTY
C  This is a good one!!  If the filename has changed, write a warning:  UTILTY
C                                                                       UTILTY
   30   CONTINUE                                                        UTILTY
        IF(FILENM(1:LENGTH).NE.TEMP(1:LENGTH)) THEN                     UTILTY
          FILENM(1:LENGTH) = TEMP(1:LENGTH)                             UTILTY
          WRITE(LFNPR,910) FILENM(1:52)                                 UTILTY
        END IF                                                          UTILTY
      END IF                                                            UTILTY
C                                                                       UTILTY
C  Open external files:                                                 UTILTY
C                                                                       UTILTY
      TEMP = FILENM                                                     UTILTY
      LENGTH = LENNB(TEMP)                                              UTILTY
      DO 40 J = 1,NFILE                                                 UTILTY
        IF(ABS(IFILE(J)).LT.10) THEN                                    UTILTY
          L2 = 1                                                        UTILTY
          WRITE(EXT,'(I1)') ABS(IFILE(J))                               UTILTY
        ELSE IF(ABS(IFILE(J)).LT.100) THEN                              UTILTY
          L2 = 2                                                        UTILTY
          WRITE(EXT,'(I2)') ABS(IFILE(J))                               UTILTY
        ELSE IF(ABS(IFILE(J)).LT.1000) THEN                             UTILTY
          L2 = 3                                                        UTILTY
          WRITE(EXT,'(I3)') ABS(IFILE(J))                               UTILTY
        END IF                                                          UTILTY
        TEMP = TEMP(1:LENGTH) // '.' // EXT(1:L2)                       UTILTY
        LEN2 = LENNB(TEMP)                                              UTILTY
        IF(IFILE(J).GT.0) THEN                                          UTILTY
          OPEN(UNIT=IFILE(J), FILE=TEMP(1:LEN2), STATUS='NEW')          UTILTY
        ELSE                                                            UTILTY
          OPEN(UNIT=ABS(IFILE(J)), FILE=TEMP(1:LEN2), STATUS='OLD',     UTILTY
     +          ERR=50)                                                 UTILTY
        END IF                                                          UTILTY
   40 CONTINUE                                                          UTILTY
      ERROR = .FALSE.                                                   UTILTY
      RETURN                                                            UTILTY
C                                                                       UTILTY
   50 WRITE(LFNPR,920) TEMP                                             UTILTY
      ERROR = .TRUE.                                                    UTILTY
      RETURN                                                            UTILTY
C                                                                       UTILTY
  900 FORMAT(/1X,'THE SEARCH FOR AN ACCEPTABLE FILENAME HAS FAILED.')   UTILTY
  910 FORMAT(/1X,'FILENAME CHANGED TO ',A52)                            UTILTY
  920 FORMAT(/1X,'ROUTINE SQFILE COULD NOT FIND FILE ',A)               UTILTY
      END                                                               UTILTY
C***********************************************************************UTILTY
C*****************************************************************************
C End of NBO 5.0, Natural Bond Orbital Analysis Programs
C*****************************************************************************
C NBO 5.0 -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1994 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C*****************************************************************************
