c ==================================================================
c  GEOMETRY OPTIMIZATION ROUTINES R-Z          JB   October 1999
c ==================================================================
c
      SUBROUTINE RdCnnct(NAtoms,NCnnct,MaxC,IPRNT,IAC,IC)
      IMPLICIT INTEGER(A-Z)
C
C
C  Attempt to read additional connectivity data from <opt> file.
C  This data will be used in addition to automatically generated
C  connectivity data based on distance criteria
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NCnnct  -  number of lines of connectivity data
C  MaxC    -  maximum atomic connectivity
C  IPRNT   -  print flag
C  IAC     -  scratch array
C  IC      -  on exit contains connectivity matrix
C
      INTEGER IC(NAtoms,NAtoms),IAC(MaxC)
      CHARACTER*80 CHAR
      CHARACTER*256 jobname
c
      Common /job/jobname,lenJ
c
      IOut=ioutfil('iout')
C
C  open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $opt section
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  locate start of connectivity section
C
 20   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:8).NE.'$connect') GO TO 20
C
C  connectivity data found
C
      If(IPRNT.GT.1) WRITE(IOut,1000)
C
C  EXPECTED FORMAT
C   atom number     list of connected atoms
C
      DO 40 I=1,NCnnct
c
      CALL IZeroIT(IAC,MaxC)
      READ(40,910,ERR=97) IAtom,(IAC(L),L=1,MaxC)
c
      DO 30 J=1,MaxC
      JJ = IAC(J)
      If(JJ.EQ.0) GO TO 40
c
      If(IAtom.GT.NAtoms.OR.JJ.GT.NAtoms) Then
       Call nerror(37,'OPTIMIZE module',
     $ 'Connectivity Data (Line # given) Contains too Large Entry',
     $  I,NAtoms)
      EndIF
c
      IC(IAtom,JJ) = 1
      IC(JJ,IAtom) = 1
 30   CONTINUE
c
 40   CONTINUE
cc
      CLOSE (UNIT=40,STATUS='KEEP')
cc
C
C  All data read
C  check that no atom is connected to itself
C
      DO 50 I=1,NAtoms
      IF(IC(I,I).NE.0) THEN
       Call nerror(38,'OPTIMIZE module',
     $  'Atom (# given) is Connected to itself!  Check input data',
     $   I,0)
      EndIF
 50   CONTINUE
C
C  If we get here, we're OK
C
      RETURN
C
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no OPT file found
      Call nerror(39,'OPTIMIZE module',
     $  'There Should be Connectivity Data but no <opt> file found!',
     $   0,0)
c
 96   CONTINUE          ! no connectivity data in OPT file
      Call nerror(40,'OPTIMIZE module',
     $ 'There Should be Connectivity Data but none found in <opt> file!'
     $ ,0,0)
c
 97   CONTINUE          ! error reading connectivity data
      Call nerror(41,'OPTIMIZE module',
     $      'Problem Reading Atomic Connectivity',0,0)
c
  900 Format(A80)
  910 Format(I4,2X,8I4)
 1000 FORMAT(' Additional Atomic Connectivity Read from <opt> file')
c
      END
c  =======================================================================
c
      SUBROUTINE RdCON(MaxCON, ICons,  NCons,  NComp,  NPComp,
     $                 IOptC,  NAtoms, XC,     ICTYP,  RCON,
     $                 IC,     ICOMP,  IPComp, PCWght)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Read user defined constraints from <opt> file and set up
C  default constraint flag if none given
C
C  ARGUMENTS
C
C  MaxCON  -  maximum number of constraints allowed
C  ICons   -  constraint flag for Cartesian optimization
C              0 - full optimization (no constraints)
C              1 - constrained optimization using Lagrange multipliers
C              2 - constrained optimization using Penalty functions
C             -1 - attempt a constrained optimization using
C                  Lagrange multipliers; if this fails switch
C                  to Penalty functions. At convergence, tidy up
C                  using Lagrange multipliers if this failed originally
C             -2 - constrained optimization using Penalty functions
C                  tidy up at convergence using Lagrange multipliers
C  NCons   -  total number of constraints
C  NComp   -  number of composite constraints
C  NPComp  -  number of primitives in composite constraints
C  IOptC   -  coordinate type for optimization
C              0 - Cartesian coordinates
C              1 - delocalized internal coordinates
C              2 - Z-matrix coordinates
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C
C  on exit ...
C
C  ICTYP   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane-bend
C              4 - fixed dihedral angle
C              5 - fixed linear coplanar bend
C              6 - fixed linear perpendicular bend
C              9 - composite constraint
C  RCON    -  constraint values
C             read in from OPT file in angstroms/degrees
C             on exit converted to atomic units
C             If no value given, constraint will be taken as value in
C             starting geometry
C  IC      -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   all other constraints
C  ICOMP   -  number of primitives in each composite constraint
C  IPComp  -  constraint type and atoms involved in constraint
C               IPComp(1,I) - constraint type (same definition as ICTYP array)
C               IPComp(2,I) to IPComp(5,I) - atoms in constraint
C  PCWght  -  weight of each primitive in composite constraint
C             (if no value given, assumed to be unity)
C
      REAL*8 RCon(NCons)
      DIMENSION ICTYP(NCons),IC(4,NCons)
      DIMENSION ICOMP(NComp),IPComp(5,NPComp),PCWght(NPComp)
      Logical found
c
      CHARACTER CHAR*80,type*4
      Character*256 jobname
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      Common /job/jobname,lenJ
C
      PARAMETER (TORAD=3.14159 26535 89793d0/180.0d0)
      PARAMETER (Zero=0.0d0,One=1.0d0,AngMAX=180.000001d0*TORAD)
cc
cc  If any dihedral angle is constrained to be within DTOL1
cc  of 0 or 180 the default is to use Penalty functions and
cc  attempt to tidy up using Lagrange multipliers
cc  If any dihedral angle is constrained to be within DTOL2
cc  of 0 or 180 use Penalty functions only (do not tidy up)
cc
      PARAMETER (DTOL1=20.0d0*TORAD)
      PARAMETER (DTOL2=10.0d0*TORAD)
C
C
C  check number of constraints
C
      If(NCons.GT.MaxCON) GO TO 94
      CALL IZeroIT(IC,4*NCons)
      If(NComp.GT.0) Then
       ncmp = 0        ! # composite constraints counter
       nt = 0          ! # primitives counter
       CALL IZeroIT(IPComp,5*NPComp)
      EndIf
C
C  open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $optimize section
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  locate start of constraints section
C
 20   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:11).NE.'$constraint') GO TO 20
C
C  constraints located
C
      DO 30 I=1,NCons
c
      READ(40,900) CHAR
C
C  how many "parameters" are on this card?
C
      CALL NumFIELD(CHAR,80,NParam)
c
      type = CHAR(1:4)
c
      IF(type.EQ.'stre') THEN
cc
C  distance constraint
C
       ICTYP(I) = 1
c
       IF(NParam.EQ.4) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),RCON(I)
        If(RCON(I).LT.Zero) GO TO 97
        RCON(I) = RCON(I)*ANTOAU
       ELSE IF(NParam.EQ.3) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I)
        CALL StreGRAD(NAtoms,IC(1,I),IC(2,I),XC,RCON(I),.false.,jnk)
       ELSE
        Call nerror(50,'OPTIMIZE module',
     $       'Problems with input for Distance Constraint',0,0)
       ENDIF
cc
      ELSE IF(type.EQ.'bend') THEN
cc
C  bond angle constraint
C
       ICTYP(I) = 2
c
       IF(NParam.EQ.5) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),RCON(I)
        RCON(I) = RCON(I)*TORAD
        If(RCON(I).LT.Zero.OR.RCON(I).GT.AngMAX) GO TO 98
       ELSE IF(NParam.EQ.4) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I)
        CALL AngGRAD(NAtoms,IC(1,I),IC(2,I),IC(3,I),XC,RCON(I),
     $               .false.,jnk)
       ELSE
        Call nerror(51,'OPTIMIZE module',
     $       'Problems with input for Angle Constraint',0,0)
       ENDIF
cc
      ELSE IF(type.EQ.'outp') THEN
cc
C  out-of-plane bend constraint
C
       ICTYP(I) = 3
c
       IF(NParam.EQ.6) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I),RCON(I)
        RCON(I) = RCON(I)*TORAD
        If(Abs(RCON(I)).GT.AngMAX) GO TO 99
       ELSE IF(NParam.EQ.5) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I)
        CALL OutpGRAD(NAtoms,IC(1,I),IC(2,I),IC(3,I),IC(4,I),XC,
     $                RCON(I),.false.,jnk)
       ELSE
        Call nerror(52,'OPTIMIZE module',
     $       'Problems with input for Out-of-Plane Bend Constraint',0,0)
       ENDIF
cc
      ELSE IF(type.EQ.'tors') THEN
cc
C  dihedral angle constraint
C
       ICTYP(I) = 4
c
       IF(NParam.EQ.6) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I),RCON(I)
        RCON(I) = RCON(I)*TORAD
        If(Abs(RCON(I)).GT.AngMAX) GO TO 100
       ELSE IF(NParam.EQ.5) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I)
        CALL DihGRAD(NAtoms,IC(1,I),IC(2,I),IC(3,I),IC(4,I),XC,
     $               RCON(I),.false.,jnk)
       ELSE
        Call nerror(53,'OPTIMIZE module',
     $       'Problems with input for Dihedral Angle Constraint',0,0)
       ENDIF
cc
      ELSE IF(type.EQ.'linc') THEN
cc
C  linear coplanar bend
C
       ICTYP(I) = 5
c
       IF(NParam.EQ.6) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I),RCON(I)
        RCON(I) = RCON(I)*TORAD
        If(Abs(RCON(I)).GT.AngMAX) GO TO 101
       ELSE IF(NParam.EQ.5) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I)
        CALL LincGRAD(NAtoms,IC(1,I),IC(2,I),IC(3,I),IC(4,I),XC,
     $                RCON(I),.false.,jnk)
       ELSE
        Call nerror(54,'OPTIMIZE module',
     $       'Problems with input for Colinear Bend Constraint',0,0)
       ENDIF
cc
      ELSE IF(type.EQ.'linp') THEN
cc
C  linear perpendicular bend
C
       ICTYP(I) = 6
c
       If(NParam.EQ.6) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I),RCON(I)
        RCON(I) = RCON(I)*TORAD
        If(Abs(RCON(I)).GT.AngMAX) GO TO 101
       ELSE IF(NParam.EQ.5) THEN
        READ(CHAR(5:80),*) IC(1,I),IC(2,I),IC(3,I),IC(4,I)
        CALL LinpGRAD(NAtoms,IC(1,I),IC(2,I),IC(3,I),IC(4,I),XC,
     $                RCON(I),.false.,jnk)
       ELSE
        Call nerror(55,'OPTIMIZE module',
     $       'Problems with input for Coplanar Bend Constraint',0,0)
       ENDIF
cc
      ELSE IF(type.EQ.'#com') THEN
cc
C  composite constraint
C  **WARNING** currently does NOT include linear bends
C
       ICTYP(I) = 9
       ncmp = ncmp+1   ! # composite contraint counter
       val = 0.0d0     ! value of constraint
       found = .False.
C
C  need to read in this constraint and determine the
C  number of primitives involved
C  data is of the form:  type   atoms   weight(if any)  value (if any)
C  NOTE: If you want to give a value, you must give a weight
C
 25    CONTINUE
       READ(40,900) CHAR
C
C  how many "parameters" are on this card?
C
       CALL NumFIELD(CHAR,80,NParam)
c
       type = CHAR(1:4)
c
       IF(type.EQ.'#end') THEN
cc
C  end of composite constraint definition
C
        ICOMP(ncmp) = nt
        GO TO 30
cc
       ELSE IF(type.EQ.'stre') THEN
cc
C  distance component
C
        nt = nt+1
        IPComp(1,nt) = 1
c
        IF(NParam.EQ.5) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),PCWght(nt),Value
         found = .True.
        ELSE IF(NParam.EQ.4) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),PCWght(nt)
        ELSE IF(NParam.EQ.3) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt)
         PCWght(nt) = One
        ELSE
         Call nerror(50,'OPTIMIZE module',
     $    'Problems with input for Distance Composite Constraint',0,0)
        ENDIF
c
        CALL StreGRAD(NAtoms,IPComp(2,nt),IPComp(3,nt),XC,va,
     $                .false.,jnk)
        val = val+va*PCWght(nt)
cc
       ELSE IF(type.EQ.'bend') THEN
cc
C  bond angle component
C
        nt = nt+1
        IPComp(1,nt) = 2
c
        IF(NParam.EQ.6) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      PCWght(nt),Value
         found = .True.
        ELSE IF(NParam.EQ.5) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      PCWght(nt)
        ELSE IF(NParam.EQ.4) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt)
         PCWght(nt) = One
        ELSE
         Call nerror(51,'OPTIMIZE module',
     $    'Problems with input for Angle Composite Constraint',0,0)
        ENDIF
c
        CALL AngGRAD(NAtoms,IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),XC,
     $               va,.false.,jnk)
        val = val+va*PCWght(nt)
cc
       ELSE IF(type.EQ.'outp') THEN
cc
C  out-of-plane bend component
C
        nt = nt+1
        IPComp(1,nt) = 3
c
        IF(NParam.EQ.7) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      IPComp(5,nt),PCWght(nt),Value
         found = .True.
        ELSE IF(NParam.EQ.6) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      IPComp(5,nt),PCWght(nt)
        ELSE IF(NParam.EQ.5) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      IPComp(5,nt)
         PCWght(nt) = One
        ELSE
         Call nerror(52,'OPTIMIZE module',
     $    'Problems with input for Out-of-Plane Bend Composite',
     $    ' Constraint',0,0)
        ENDIF
c
        CALL OutpGRAD(NAtoms,IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                IPComp(5,nt),XC,va,.false.,jnk)
        val = val+va*PCWght(nt)
cc
       ELSE IF(type.EQ.'tors') THEN
cc
C  torsion component
C
        nt = nt+1
        IPComp(1,nt) = 4
c
        IF(NParam.EQ.7) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      IPComp(5,nt),PCWght(nt),Value
         found = .True.
        ELSE IF(NParam.EQ.6) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      IPComp(5,nt),PCWght(nt)
        ELSE IF(NParam.EQ.5) THEN
         READ(CHAR(5:80),*) IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $                      IPComp(5,nt)
         PCWght(nt) = One
        ELSE
         Call nerror(53,'OPTIMIZE module',
     $    'Problems with input for Torsion Composite Constraint',0,0)
        ENDIF
c
        CALL DihGRAD(NAtoms,IPComp(2,nt),IPComp(3,nt),IPComp(4,nt),
     $               IPComp(5,nt),XC,va,.false.,jnk)
        val = val+va*PCWght(nt)
cc
       ELSE
cc
C  unknown/unacceptable composite constraint type
C
        GO TO 103
cc
       ENDIF
cc
       RCON(I) = val      ! value for composite constraint
       If(found) RCON(I) = Value
       GO TO 25
cc
      ELSE
cc
C  unknown constraint type
C
       GO TO 102
cc
      ENDIF
 30   CONTINUE
C
C  If there are composite constraints, check that we have
C  the correct number and composition
C
      If(NComp.GT.0) Then
       If(ncmp.NE.NComp.AND.nt.NE.NPComp) GO TO 104
C
C  reassign ICOMP array
C
       If(NComp.GT.1) Then
        DO I=NCmp,2,-1
        ICOMP(I) = ICOMP(I) - ICOMP(I-1)
        EndDO
       EndIf
      EndIf
C
C  If ICons has not been set, select default
C
      IF(ICons.EQ.0) THEN
       ICons = -1
       DO 40 I=1,NCons
       IF(ICTYP(I).EQ.4.AND.IOptC.EQ.0) THEN
        Ang = Abs(RCON(I))
        If(Ang.LT.DTOL1.OR.(AngMax-Ang).LT.DTOL1) ICons = -2
        If(Ang.LT.DTOL2.OR.(AngMax-Ang).LT.DTOL2) ICons =  2
       ENDIF
 40    CONTINUE
      ENDIF
C
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C  ..............................................
C    ERROR SECTION
C
 94   CONTINUE          ! too many constraints
      Call nerror(28,'OPTIMIZE module',
     $  'Trying to Impose Too Many Constraints.  Maximum Allowed is',
     $   MaxCON,0)
c
 95   CONTINUE          ! no OPT file found
      Call nerror(29,'OPTIMIZE module',
     $  'There Should be Constraints but no <opt> file found!',NCons,0)
c
 96   CONTINUE          ! no constraints found
      Call nerror(30,'OPTIMIZE module',
     $  'There Should be Constraints but none found in <opt> file!',
     $   NCons,0)
c
 97   CONTINUE          ! zero or negative distance constraint requested
      Call nerror(31,'OPTIMIZE module',
     $  'Trying to Impose Zero or Negative Distance Constraint',0,0)
c
 98   CONTINUE          ! bond angle constraint out of range
      Call nerror(32,'OPTIMIZE module',
     $  'Bond Angle Constraint must lie between 0 and 180 degrees',0,0)
c
 99   CONTINUE          ! out-of-plane bend constraint out of range
      Call nerror(33,'OPTIMIZE module',
     $'Out-of-Plane Angle Constraint must lie between 0 and 180 degrees'
     $  ,0,0)
c
 100  CONTINUE          ! dihedral angle constraint out of range
      Call nerror(34,'OPTIMIZE module',
     $ 'Dihedral Angle Constraint must lie between -180 and 180 degrees'
     $  ,0,0)
c
 101  CONTINUE          ! colinear bends out of range
      Call nerror(35,'OPTIMIZE module',
     $ 'Colinear Bend Constraint must lie between -180 and 180 degrees'
     $  ,0,0)
c
 102  CONTINUE          ! unknown constraint type
      CHAR = 'Unknown Constraint Type: '//type
      Call nerror(36,'OPTIMIZE module',CHAR,0,0)
c
 103  CONTINUE
      CHAR = 'Unknown Composite Constraint Type: '//type
      Call nerror(37,'OPTIMIZE module',CHAR,0,0)
c
 104  CONTINUE          ! error in composite constraints
      Call nerror(38,'OPTIMIZE module',
     $ 'Error in number or composition of Composite Constraints',0,0)
c
  900 Format(A80)
c
      END
c  =======================================================================
c
      SUBROUTINE RdDRIVE(NDrive, IDRTYP, IDRIVE, FDRIVE)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads <opt> file for user defined coordinate driving
C
C  ARGUMENTS
C
C  NDrive  -  number of primitives to drive
C  IDRTYP  -  type of each primitive to be driven
C               1 - distance               stre
C               2 - bond angle             bend
C               3 - out-of-plane bend      outp
C               4 - dihedral angle         tors
C  IDRIVE  -  definitive of each primitive to orive
C  FDRIVE  -  force (in internal coordinates) to be applied
C
C
      DIMENSION IDRTYP(NDRive),IDRIVE(4,NDrive),FDRIVE(NDrive)
      CHARACTER CHAR*80,type*4
      CHARACTER*256 jobname
c
      Common /job/jobname,lenJ
C
C
C  open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $optimize section
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  locate start of coordinate driving section
C
 20   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:6).NE.'$drive') GO TO 20
C
C  section located
C
      DO 30 I=1,NDrive
c
      READ(40,900) CHAR
      type = CHAR(1:4)
c
      IF(type.EQ.'stre') THEN
       IDRTYP(I) = 1
       READ(CHAR(5:80),*) IDRIVE(1,I),IDRIVE(2,I),FDRIVE(I)
      ELSE IF(type.EQ.'bend') THEN
       IDRTYP(I) = 2
       READ(CHAR(5:80),*) IDRIVE(1,I),IDRIVE(2,I),IDRIVE(3,I),FDRIVE(I)
      ELSE IF(type.EQ.'outp') THEN
       IDRTYP(I) = 3
       READ(CHAR(5:80),*) IDRIVE(1,I),IDRIVE(2,I),IDRIVE(3,I),
     $                    IDRIVE(4,I),FDRIVE(I)
      ELSE IF(type.EQ.'tors') THEN
       IDRTYP(I) = 4
       READ(CHAR(5:80),*) IDRIVE(1,I),IDRIVE(2,I),IDRIVE(3,I),
     $                    IDRIVE(4,I),FDRIVE(I)
      ELSE
       Call nerror(73,'OPTIMIZE module',
     $      'Unknown primitive type for Coordinate Driving',0,0)
      ENDIF
c
 30   CONTINUE
C
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no OPT file found
      Call nerror(74,'OPTIMIZE module',
     $  'There Should be Coordinate Driving but no <opt> file found!',
     $   NDrive,0)
c
 96   CONTINUE          ! no driving coordinates found
      Call nerror(75,'OPTIMIZE module',
     $  'There Should be Coordinates to Drive but none found!',
     $   NDrive,0)
c
  900 Format(A80)
c
      END
c  =======================================================================
c
      SUBROUTINE RdDUM(NDum,MaxL,IFunc)
      IMPLICIT INTEGER(A-Z)
C
C
C  Reads <opt> file for user defined dummy atoms
C
C  ARGUMENTS
C
C  NDum    -  number of dummy atoms
C  MaxL    -  maximum number of real atoms involved in
C             definition of dummy atom position + 1
C  IFunc   -  list for each dummy atom of the atom type
C             and the real atoms defining its position
C
C  All dummy atoms are defined with reference to a list of
C  real atoms and dummy atom coordinates will be generated
C  from the coordinates of the real atoms in its defining
C  list. There are two types of dummy atom:
C
C    1.  Positioned at the arithmetic mean of the up to
C        MaxL-1 real atoms in the defining list
C    2.  Positioned a unit distance along the normal to
C        a plane defined by 3 atoms, centred on the middle
C        atom of the 3
C
C
      DIMENSION IFunc(NDum,MaxL)
      CHARACTER*80 CHAR
      CHARACTER*256 jobname
c
      Common /job/jobname,lenJ
C
C
C  open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $optimize section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  locate start of dummy atom section
C
 20   CONTINUE
      READ(40,900,END=97) CHAR
      If(CHAR(1:6).NE.'$dummy') GO TO 20
C
C  dummy atoms located
C
      CALL IZeroIT(IFunc,NDum*MaxL)
c
      DO 30 I=1,NDum
      READ(40,910) IAT,IFunc(I,1),(IFunc(I,J),J=2,MaxL)
 30   CONTINUE
C
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no OPT file found
      Call nerror(25,'OPTIMIZE module',
     $  'There Should be Dummy Atoms but no <opt> file found!',NDum,0)
c
 96   CONTINUE          ! no $optimize section found
      Call nerror(26,'OPTIMIZE module',
     $      'No $optimize section found on <opt> file!',0,0)
c
 97   CONTINUE          ! no dummy atoms found
      Call nerror(27,'OPTIMIZE module',
     $  'There Should be Dummy Atoms but none found in <opt> file!',
     $   NDum,0)
c
  900 FORMAT(A80)
  910 FORMAT(I4,2X,I4,2X,7I4)
c
      END
c  =======================================================================
c
      SUBROUTINE RdFIX(NAtoms,NFix,IFix)
      IMPLICIT INTEGER(A-Z)
C
C
C  Reads <opt> file for user defined fixed atomic
C  (Cartesian) coordinates
C  NOTE: "Fixed" coordinates will not be included in the
C        optimization space and will not move EXCEPT for
C        reorientation to centre of mass axis system
C
C  ARGUMENTS
C
C  NAtoms  -  number of real atoms
C  NFix    -  on input  number of data lines to read
C             on exit   actual number of fixed coordinates
C  IFix    -  on exit contains list of fixed/active coordinates
C              0 - coordinate active
C              1 - coordinate will be fixed during optimization
C  ** IMPORTANT:  IFix must be zeroed prior to calling RdFIX **
C
      DIMENSION IFix(3*NAtoms)
      CHARACTER CHAR*80,jobname*256
c
      Common /job/jobname,lenJ
C
C
C  zero coordinate counter
C
      inf = 0
C
C  open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $optimize section
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  locate start of fixed coordinate section
C
 20   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:4).NE.'$fix') GO TO 20
C
C  fixed coordinates located
C
      DO 30 I=1,NFix
c
      READ(40,910) IAT,CHAR
C
C  the character string CHAR indicates which coordinates
C  (X and/or Y and/or Z) are to be fixed for atom IAT
C  set corresponding IFix entries
C
      If(IAT.LT.1.OR.IAT.GT.NAtoms) Then
       Call nerror(22,'OPTIMIZATION module',
     $   'Trying to Fix Coordinates of Non-Existent Atom!',IAT,NAtoms)
      EndIf
c
      IAT = 3*(IAT-1)
c
      IF(CHAR(1:1).EQ.'X'.OR.CHAR(2:2).EQ.'X'.OR.CHAR(3:3).EQ.'X') THEN
       II = IAT+1
       IFix(II) = 1
       inf = inf+1
      ENDIF
c
      IF(CHAR(1:1).EQ.'Y'.OR.CHAR(2:2).EQ.'Y'.OR.CHAR(3:3).EQ.'Y') THEN
       II = IAT+2
       IFix(II) = 1
       inf = inf+1
      ENDIF
c
      IF(CHAR(1:1).EQ.'Z'.OR.CHAR(2:2).EQ.'Z'.OR.CHAR(3:3).EQ.'Z') THEN
       II = IAT+3
       IFix(II) = 1
       inf = inf+1
      ENDIF
c
 30   CONTINUE
C
C  set actual number of fixed coordinates
C
      NFix = inf
c
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE           ! no OPT file found
       Call nerror(23,'OPTIMIZATION module',
     $   'There Should be Fixed Atoms but no <opt> file found!',NFix,0)
C
 96   CONTINUE           ! no fixed coordinates found
       Call nerror(24,'OPTIMIZATION module',
     $   'There Should be Fixed Atoms but none found in <opt> file!',
     $    NFix,0)
c
  900 FORMAT(A8)
  910 FORMAT(I4,2X,A3)
c
      END
c  =======================================================================
c
      SUBROUTINE RdHESS(HFile,Len,NDim,IPRNT,HESS,IHflag)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Attempt to read Hessian matrix from file <HFile>
C
C  ARGUMENTS
C
C  HFile   -  File name
C  Len     -  Length of character string HFile
C  NDim    -  Dimension of Hessian
C  IPRNT   -  Print flag
C  HESS    -  Hessian matrix
C  IHflag  -  Hessian found flag
C              0 - Hessian read successfully
C             -1 - no Hessian found
C
C
      DIMENSION HESS(NDim,NDim)
      CHARACTER*(*) HFile
      CHARACTER*80 CHAR
C
      IOut=ioutfil('iout')
C
      IHflag = -1
C
C  open Hessian file
C
      OPEN (UNIT=40,FILE=HFile,FORM='FORMATTED',STATUS='OLD',ERR=95)
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:5).NE.'$hess') GO TO 10
c
      READ(40,910) NDim1
c
      If(NDim1.NE.NDim) Then
       CHAR = 'Hessian dimension on file <'//HFile//
     $        '> incompatible with current system'
       Call nerror(42,'OPTIMIZE module',CHAR,0,0)
      EndIf
C
C  now read the Hessian
C
      DO 20 I=1,NDim
      READ(40,*) (HESS(I,J),J=1,I)
      DO 19 J=1,I
      HESS(J,I) = HESS(I,J)
 19   CONTINUE
 20   CONTINUE
C
C  check for end of data
C
      READ(40,900) CHAR
      If(CHAR(1:4).NE.'$end') Then
       CHAR = 'No End-of-File $end found on Hessian file'//HFile
       Call nerror(43,'OPTIMIZE module',CHAR,0,0)
      EndIf
c
      IHflag = 0
      If(IPRNT.GT.1) WRITE(IOut,*) 'Hessian read from file ',HFile
c
      CLOSE(UNIT=40,STATUS='KEEP')
      RETURN
C
 95   CONTINUE
      RETURN
C  ..............................................
C    ERROR SECTION
C
 96   CONTINUE
      CHAR = 'No $hess marker found on Hessian file'//HFile
      Call nerror(44,'OPTIMIZE module',CHAR,0,0)
c
 900  Format(A8)
 910  Format(11X,I5)
 1100 FORMAT(/,' Hessian read from file ',A55)
c
      END
c  =======================================================================
c
      SUBROUTINE RdHPRIM(NDim,HESS,NCol,ICNUM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Routine to read Primitive Hessian and constraint data from
C  generic HPRIM file
C
C  ARGUMENTS
C
C  NDim    -  Dimension of Hessian
C  HESS    -  Hessian matrix
C  NCol    -  number of fixed primitives
C  ICNUM   -  list of fixed primitives
C             (-1 appended to end)
C
C
      REAL*8 HESS(NDim,NDim)
      INTEGER ICNUM(*)
      CHARACTER*80 CHAR
      CHARACTER*256 jobname
C
      Common /job/jobname,lenJ
c
      IOut=ioutfil('iout')
C
C  Open Hessian file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.hprim',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  first look for primitive Hessian
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:9).NE.'$hessprim') GO TO 10
C
C  Primitive Hessian found
C  Read and check the Hessian dimension
C
      READ(40,910) NDim1
c
      IF(NDim1.NE.NDim) THEN
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
C  now read the Hessian
C
      DO 20 I=1,NDim
      READ(40,*) (HESS(I,J),J=1,I)
      DO 19 J=1,I
      HESS(J,I) = HESS(I,J)
 19   CONTINUE
 20   CONTINUE
C
C  check for end of data
C
      READ(40,900) CHAR
      If(CHAR(1:4).NE.'$end') GO TO 96
C
C  read the fixed primitive data
C
      READ(40,920) NCol
      READ(40,930) (ICNUM(I),I=1,NCol)
C
C  append -1  (to signify end of data)
C
      ICNUM(NCol+1) = -1
c
      READ(40,900) Char
c
      If(Char(1:4).NE.'$end') GO TO 95
c
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C  .................................................
C  ERROR SECTION
C
 95   CONTINUE
      WRITE(IOut,1100)
      CALL OptExit(9)
c
 96   CONTINUE
      WRITE(IOut,1200)
      CALL OptExit(9)
c
 900  Format(A80)
 910  Format(11X,I5)
 920  Format(18X,I5)
 930  Format(10I5)
 1000 FORMAT(/,2X,'***ERROR*** Primitive Hessian on HPRIM file',
     $            ' Incompatible with Current System')
 1100 FORMAT(/,2X,'***ERROR*** Unable to find HPRIM file in <RdHPRIM>')
 1200 FORMAT(/,2X,'***ERROR*** Problems Reading Primitive Hessian Data',
     $              ' in <RdHPRIM>')
c
      END
c  =======================================================================
c
      SUBROUTINE RdOPT0(IOptC,  IGen,   IType,  NFix,   NDum,
     $                  NCons,  NComp,  NPComp, NCnnct, NDrive,
     $                  MaxDiis)
      IMPLICIT INTEGER(A-Z)
C
C
C  Reads from the OPT file those optimization parameters
C  that affect memory requirements
C
C  ARGUMENTS
C
C  IOptC    -  coordinate type flag
C    0 - optimize in Cartesian coordinates
C    1 - generate and optimize in internal coordinates
C        if this fails abort
C   -1 - generate and optimize in internal coordinates
C        if this fails during any stage of the optimization
C        switch to Cartesians and continue
C    2 - optimize in Z-Matrix internal coordinates
C        if this fails abort
C   -2 - optimize in Z-Matrix internal coordinates
C        if this fails during any stage of the optimization
C        switch to Cartesians and continue
C    3 - large-molecule optimization in Cartesian coordinates
C    4 - large-molecule optimization in delocalized internals
C  IGen     -  generation of delocalized internal coordinates
C   -1 - generate a set of natural internal coordinates
C        and use throughout the optimization
C    0 - generate a set of non-redundant internals
C        and use throughout the optimization
C    1 - generate a new set of non-redundant internals
C        from the same primitive space every cycle
C    2 - generate a new set of underlying primitives
C        and a new set of non-redundant internals every cycle
C  IType    -  special flag for surface/cluster optimization
C    0 - standard molecular optimization
C    1 - optimization of molecule reacting/adsorbed on model surface
C    2 - optimization of molecular clusters
C  NFix     -  number of "fixed" atoms or number of lines
C              of surface constraint data
C  NDum     -  number of dummy atoms
C  NCons    -  total number of geometric constraints imposed
C  NComp    -  number of composite constraints
C  NPComp   -  total number of primitives in composite constraints
C  NCnnct   -  number of lines of atom connectivity data
C  NDrive   -  number of coordinates to drive
C  MaxDiis  -  maximum size of subspace for GDIIS
C    0 - do not use GDIIS
C   -1 - default size (depends on system; typically 3 or 4)
C    N - size specified by user (DO NOT set N too large)

C
      CHARACTER CHAR*80
      Character*256 jobname
c
      Common /job/jobname,lenJ
C
C
C  Set up default parameters
C
      IOptC = -1
      IGen = 0
      IType = 0
      NFix = 0
      NDum = 0
      NCons = 0
      NComp = 0
      NPComp = 0
      NCnnct = 0
      NDrive = 0
      MaxDiis = 0
C
C  now attempt to open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  find $optimize section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:4).NE.'$opt') GO TO 10
C
C  read a line at a time and attempt to decypher the keywords
C
 200  CONTINUE
      READ(40,900,END=97) CHAR
cc
      If(CHAR(1:4).EQ.'$end') GO TO 94      ! end of $optimize section
cc
      IF(CHAR(1:8).EQ.'$optcoor') THEN
cc
       READ(CHAR(11:80),*) IOptC
       If(Abs(IOptC).GT.4) Then
        Call nerror(1,'OPTIMIZE module',
     $        'Unallowed Coordinate Type for Optimization',0,0)
       EndIf
cc
      ELSE IF(CHAR(1:8).EQ.'$regener') THEN
cc
       READ(CHAR(13:80),*) IGen
       If(IGen.LT.0.AND.IGen.GT.2) Then
        Call nerror(2,'OPTIMIZE module',
     $   'Unallowed $regenerate option.  Allowed Values: 0, 1 or 2',
     $    0,0)
       EndIf
cc
      ELSE IF(CHAR(1:6).EQ.'$gdiis') THEN
cc
       READ(CHAR(8:80),*) MaxDiis
       If(MaxDiis.LT.-1) Then
        Call nerror(3,'OPTIMIZE module',
     $   'Subspace Size for $gdiis Must be a Positive Integer',0,0)
       EndIf
cc
      ELSE IF(CHAR(1:8).EQ.'$coordty') THEN
cc
       READ(CHAR(11:80),*) IType
       If(IType.LT.0.AND.IType.GT.2) Then
        Call nerror(2,'OPTIMIZE module',
     $   'Unallowed $coordtype option.  Allowed Values: 0, 1 or 2',
     $    0,0)
       EndIf
cc
      ELSE IF(CHAR(1:4).EQ.'$fix') THEN
cc
C  there are fixed cartesian coordinates
C  determine how many atoms are involved
C
 20    CONTINUE
       READ(40,900,END=98) CHAR
       If(CHAR(1:7).EQ.'$endfix') GO TO 200
c
       NFix = NFix+1
       GO TO 20
cc
      ELSE IF(CHAR(1:6).EQ.'$dummy') THEN
cc
C  there are dummy atoms
C  determine how many
C
 30    CONTINUE
       READ(40,900,END=99) CHAR
       If(CHAR(1:7).EQ.'$enddum') GO TO 200
c
       NDum = NDum+1
       GO TO 30
cc
      ELSE IF(CHAR(1:8).EQ.'$constra') THEN
cc
C  there are constraints
C  determine how many
C  note that we can now include composite constraints
C  i.e., constraints involving more than one primitive
C
 40    CONTINUE
       READ(40,900,END=100) CHAR
       If(CHAR(1:8).EQ.'$endcons') GO TO 200
c
       If(CHAR(1:8).EQ.'#composi') Then
C
C  we have a composite constraint
C  locate the ending string
C
        NComp = NComp+1
        NPComp0 = NPComp
 45     CONTINUE
        READ(40,900,END=103) CHAR
        If(CHAR(1:8).EQ.'#endcomp') GO TO 46
        NPComp = NPComp+1
        GO TO 45
c
 46     CONTINUE
        If(NPComp.EQ.NPComp0) Then
c -- a composite constraint with NO primitives
         Call nerror(8,'OPTIMIZE module',
     $    'There is a Composite Constraint with NO primitives',0,0)
        EndIf
       EndIf
c
       NCons = NCons+1
       GO TO 40
cc
      ELSE IF(CHAR(1:8).EQ.'$connect') THEN
cc
C  there is atom connectivity data
C  determine how many lines
C
 50    CONTINUE
       READ(40,900,END=101) CHAR
       If(CHAR(1:8).EQ.'$endconn') GO TO 200
c
       NCnnct = NCnnct+1
       GO TO 50
cc
      ELSE IF(CHAR(1:8).EQ.'$surface') THEN
cc
C  there is surface constraint data
C  determine how many lines
C
 60    CONTINUE
       READ(40,900,END=102) CHAR
       If(CHAR(1:8).EQ.'$endsurf') GO TO 200
c
       NFix = NFix+1
       GO TO 60
cc
      ELSE IF(CHAR(1:6).EQ.'$drive') THEN
cc
C  there are coordinates to drive
C  determine how many
C
 70    CONTINUE
       READ(40,900,END=104) CHAR
       If(CHAR(1:8).EQ.'$enddriv') GO TO 200
c
       NDrive = NDrive+1
       GO TO 70
cc
      ENDIF
      GO TO 200
c
 94   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no <opt> file found
      RETURN
c
 96   CONTINUE          ! no $optimize section found
      Call nerror(4,'OPTIMIZE module',
     $      'No $optimize keyword found on <opt> file!',0,0)
c
 97   CONTINUE          ! no $end section marker found
      Call nerror(5,'OPTIMIZE module',
     $      'No $end section marker found in <opt> file!',0,0)
c
 98   CONTINUE          ! no $endfix marker for fixed atoms
      Call nerror(6,'OPTIMIZE module',
     $      'No $endfix marker found for Fixed Atoms',0,0)
c
 99   CONTINUE          ! no $enddummy marker for dummy atoms
      Call nerror(7,'OPTIMIZE module',
     $      'No $enddummy marker found for Dummy Atoms',0,0)
c
 100  CONTINUE          ! no $endconstraint marker for constraints
      Call nerror(8,'OPTIMIZE module',
     $      'No $endconstraint marker found for Constraints',0,0)
c
 101  CONTINUE          ! no $endconnect marker for connectivity data
      Call nerror(9,'OPTIMIZE module',
     $ 'No $endconnect marker found for Atom Connectivity Data',0,0)
c
 102  CONTINUE
      Call nerror(10,'OPTIMIZE module',
     $ 'No $endsurface marker found for Surface Constraint Data',0,0)
c
 103  CONTINUE
      Call nerror(11,'OPTIMIZE module',
     $ 'No $endcomposite marker found for Composite Constraint',0,0)
c
 104  CONTINUE
      Call nerror(12,'OPTIMIZE module',
     $ 'No $enddrive marker found for Coordinate Driving',0,0)
c
  900 Format(A80)
c
      END
c =======================================================================
c
      SUBROUTINE RdOPT1(IOP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sets up default parameters for optimization and reads
C  OPT file for user defined options
C
C  ARGUMENTS
C
C  IOP      -  integer array containing all user options
C
C  These are:
C
C  IOP(1)       IOptC      COORDINATE TYPE
C    0 - optimize in Cartesian coordinates
C    1 - generate and optimize in internal coordinates
C        if this fails abort
C   -1 - generate and optimize in internal coordinates
C        if this fails during any stage of the optimization
C        switch to Cartesians and continue
C    2 - optimize in Z-matrix coordinates
C        if this fails abort
C   -2 - optimize in Z-matrix coordinates
C        if this fails during any stage of the optimization
C        switch to Cartesians and continue
C    3 - large-molecule optimization in Cartesian coordinates
C    4 - large-molecule optimization in delocalized internals
C
C  IOP(2)       IGen        GENERATION OF DELOCALIZED INTERNALS
C   -1 - generate a set of redundant primitive internal coordinates
C        and use throughout the optimization (NOT IMPLEMENTED)
C    0 - generate a set of non-redundant internals
C        and use throughout the optimization
C    1 - generate a new set of non-redundant internals
C        from the same primitive space every cycle
C    2 - generate a new set of underlying primitives
C        and a new set of non-redundant internals every cycle
C
C  IOP(3)       ICons       CONSTRAINT TYPE FOR CARTESIAN COORDINATES
C    0 - full optimization (no constraints)
C    1 - constrained optimization using Lagrange multipliers
C    2 - constrained optimization using Penalty functions
C   -1 - attempt a constrained optimization using
C        Lagrange multipliers; if this fails switch
C        to Penalty functions. At convergence, tidy up
C        using Lagrange multipliers if this failed originally
C   -2 - constrained optimization using Penalty functions
C        tidy up at convergence using Lagrange multipliers
C
C  IOP(4)       Symflag      USE OF SYMMETRY
C    0 - do not use molecular symmetry
C    1 - make use of point group symmetry
C
C  IOP(5)       TSflag       TYPE OF STATIONARY POINT SOUGHT
C    0 - minimum
C    1 - transition state
C
C  IOP(6)       mode         HESSIAN MODE FOLLOWED DURING TS SEARCH
C    0 - mode following switched off (follow lowest mode)
C    N - maximize along mode N
C
C  IOP(7)       MaxDiis      MAXIMUM SIZE OF SUBSPACE FOR GDIIS
C    0 - do not use GDIIS
C   -1 - default size   MIN(NDEG,NATOMS,4)
C    N - size specified by user (DO NOT set N too large)
C
C  IOP(8)       TolG         CONVERGENCE ON MAXIMUM GRADIENT COMPONENT
C  IOP(9)       TolD         CONVERGENCE ON MAXIMUM ATOMIC DISPLACEMENT
C  IOP(10)      TolE         CONVERGENCE ON ENERGY CHANGE
C        To converge TolG MUST be satisfied (default 0.0003) and
C        ONE of TolD (default 0.0003) or TolE (default 0.000001)
C
C  IOP(11)      STol         GRADIENT TOLERANCE FOR STEEPEST DESCENT STEP
C        If the RMS gradient for an unconstrained optimization is
C        greater than STol, a steepest descent step will be taken
C
C  IOP(12)      DMAX         MAXIMUM ALLOWED STEPSIZE
C
C  IOP(14)      IHess        HESSIAN STATUS
C    0 - have "exact" or initial Cartesian Hessian
C        use as is for Cartesian; transform if internals
C    1 - have Hessian from previous step   need to update
C   -1 - set up default (diagonal) Hessian
C   -2 - use unit matrix as starting Hessian
C
C  IOP(15)      IUpDat       HESSIAN UPDATE
C    0 - do not update the Hessian (!?)
C   -1 - use default update
C    1 - Murtagh-Sargent update
C    2 - Powell update
C    3 - Powell/Murtagh-Sargent update (default for TS)
C    4 - BFGS update (default for minimization)
C    5 - BFGS with safeguards to ensure retention of
C        positive definiteness (default for GDIIS)
C
C  IOP(16)      IProj        CARTESIAN HESSIAN PROJECTION
C    0 - do not project; leave Hessian "as is"
C    1 - project out translations
C    2 - project out translations and rotations
C
C  IOP(17)      IDB          dB/dcart TERM IN INTERNAL COORDINATE HESSIAN
C    0 - do not include this term in Hessian transformation
C        (for approximate Hessian & retention of Cartesian
C         Hessian eigenvalue structure)
C    1 - include it
C
C  IOP(18)      MaxCYC       MAXIMUM NUMBER OF OPTIMIZATION CYCLES
C
C  IOP(19)      CTol         TOLERANCE FOR SATISFYING CONSTRAINT
C        If all imposed constraints are satisfied to better than CTol
C        (default 1.0d-6) then can switch to internal coordinates
C
C  IOP(20)      Linear       THRESHOLD FOR NEAR-LINEAR BOND ANGLE
C        (default is 165 degrees)
C
C  IOP(21)      BSkal        SCALE FACTOR FOR INVERSE-DISTANCE COORDINATES
C        (default is 1 - no scaling)
C
C  IOP(22)      Tors         USE OF TORSIONS IN DELOCALIZED INTERNALS
C    0 - use torsions (default)
C    1 - do not use torsions
C
C  IOP(23)      IQmm         QM/MM OPTIMIZATION
C    0 - normal optimization
C    1 - QM/MM optimization (need to reset number of molecules)
C
C  IOP(25)      IBack        CONTROLS INTERNAL COORDINATE BACK-TRANSFORMATION
C        (default is 1 - use full back-transformation)
C
C  IOP(26)      CutOff       "BONDING" DISTANCE THRESHOLD FOR VDW
C        (default is 5 Angstroms for cluster optimizations
C              and   3 Angstroms for surface optimizations)
C
C  IOP(27)      HCnvrt       CONVERSION OF INTERNAL HESSIAN TO CARTESIANS
C    0 - convert Hessian only at convergence
C    1 - do conversion every optimization cycle
C
C  IOP(28)      IType        SPECIAL FLAG FOR SURFACE/CLUSTER OPTIMIZATION
C    0 - standard ab initio/semiempirical optimization
C    1 - surface optimization
C    2 - cluster optimization (inverse-power distance coordinates)
C
C  IOP(29)      IPRNT        AMOUNT OF PRINT OUT
C    0 - NO printout (except for error messages)
C    1 - summary and warning printout only
C    2 - "standard" printout
C    3 - slightly more printout (including gradient)
C    4 - heavier printout (including Hessian)
C    5 - heavier still (includes iterative printout)
C    6 - very heavy (including internal generation)
C    7 - debug printout
C
      CHARACTER CHAR*80,jobname*256
      INTEGER IOP(30)
C
      PARAMETER (Zero=0.0d0)
      Data CutOff /0.0d0/
c
      Common /job/jobname,lenJ
C
C
C  Set up default parameters
C
      CALL IZeroIT(IOP,30)
      IOP(1) = -1
      IOP(4) = 1
      IOP(8) = 300
      IOP(9) = 300
      IOP(10) = 100
      IOP(11) = 300
      IOP(12) = 300
      IOP(14) = -1
      IOP(15) = -1
      IOP(16) = 2
      IOP(18) = 50
      IOP(19) = 1
      IOP(20) = 165
      IOP(21) = 1000
      IOP(25) = 1
      IOP(29) = 2
C
C  now attempt to open <opt> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  find $optimize section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  read a line at a time and attempt to decypher the keywords
C
 100  CONTINUE
      READ(40,900,END=96) CHAR
cc
      If(CHAR(1:5).EQ.'$end ') GO TO 96      ! end of $optimize section
cc
      IF(CHAR(1:8).EQ.'$optcoor') THEN
cc
       READ(CHAR(11:80),*) IOptC
       If(Abs(IOptC).GT.2) Then
        Call nerror(1,'OPTIMIZE module',
     $        'Unallowed Coordinate Type for Optimization',0,0)
       EndIf
       IOP(1) = IOptC
cc
      ELSE IF(CHAR(1:8).EQ.'$contype') THEN
cc
       READ(CHAR(10:80),*) ICons
       If(Abs(ICons).GT.2) THEN
        Call nerror(10,'OPTIMIZE module',
     $      'Unallowed $contype option.  Allowed Values: -2 to 2',0,0)
       EndIf
       IOP(3) = ICons
cc
      ELSE IF(CHAR(1:8).EQ.'$transit') THEN
cc
       IOP(5) = 1
cc
      ELSE IF(CHAR(1:8).EQ.'$project') THEN
cc
       READ(CHAR(10:80),*) IProj
       If(IProj.LT.0.OR.IProj.GT.2) Then
        Call nerror(11,'OPTIMIZE module',
     $      'Unallowed $project option.  Allowed Values: 0, 1 or 2',0,0)
       EndIf
       IOP(16) = IProj
cc
      ELSE IF(CHAR(1:6).EQ.'$gdiis') THEN
cc
       READ(CHAR(8:80),*) MaxDiis
       If(MaxDiis.LT.-1) Then
        Call nerror(3,'OPTIMIZE module',
     $   'Subspace Size for $gdiis Must be a Positive Integer',0,0)
       EndIf
       IOP(7) = MaxDiis
cc
      ELSE IF(CHAR(1:5).EQ.'$gtol') THEN
cc
 1800 FORMAT(/,2X,'***ERROR*** DMAX Must be Greater than Zero')
       READ(CHAR(7:80),*) TolG
       If(TolG.LT.1.0d-6) Then
        Call nerror(12,'OPTIMIZE module',
     $    'Gradient Convergence $gtol must not be less than 10**-6',0,0)
       EndIf
       IOP(8) = NINT(TolG*10**6)
cc
      ELSE IF(CHAR(1:5).EQ.'$dtol') THEN
cc
       READ(CHAR(7:80),*) TolD
       If(TolD.LT.1.0d-6) Then
        Call nerror(13,'OPTIMIZE module',
     $    'Maximum Stepsize $dtol must not be less than 10**-6',0,0)
       EndIf
       IOP(9) = NINT(TolD*10**6)
cc
      ELSE IF(CHAR(1:5).EQ.'$etol') THEN
cc
       READ(CHAR(7:80),*) TolE
       If(TolE.LT.1.0d-8) Then
        Call nerror(14,'OPTIMIZE module',
     $    'Energy Convergence $etol must not be less than 10**-8',0,0)
       EndIf
       IOP(10) = NINT(TolE*10**8)
cc
      ELSE IF(CHAR(1:5).EQ.'$stol') THEN
cc
       READ(CHAR(7:80),*) STol
       If(STol.LT.1.0d-3) Then
        Call nerror(15,'OPTIMIZE module',
     $  'RMS Gradient for Steepest Descent ($stol) must be > 10**-3',
     $   0,0)
       EndIf
       IOP(11) = NINT(STol*10**3)
cc
      ELSE IF(CHAR(1:8).EQ.'$stepsiz') THEN
cc
       READ(CHAR(11:80),*) DMAX
       If(DMAX.LT.1.0d-3) THEN
        Call nerror(16,'OPTIMIZE module',
     $    'Maximum Allowed Stepsize must not be less than 10**-3',0,0)
       EndIf
       IOP(12) = NINT(DMAX*10**3)
cc
      ELSE IF(CHAR(1:5).EQ.'$mode') THEN
cc
       READ(CHAR(7:80),*) mode
       If(mode.LT.0) Then
        Call nerror(17,'OPTIMIZE module',
     $   'Mode Following $mode must be a Positive Integer',0,0)
       EndIf
       IOP(6) = mode
cc
      ELSE IF(CHAR(1:8).EQ.'$hessian') THEN
cc
       READ(CHAR(10:80),*) IHess
       If(IHess.GT.0.OR.IHess.LT.-3) Then
        Call nerror(18,'OPTIMIZE module',
     $   'Unallowed $hessian option.  Allowed Values: -3 to 0',0,0)
       EndIf
       IOP(14) = IHess
cc
      ELSE IF(CHAR(1:7).EQ.'$update') THEN
cc
       READ(CHAR(9:80),*) IUpDat
       If(IUpDat.GT.5.OR.IUpDat.LT.-1) Then
        Call nerror(18,'OPTIMIZE module',
     $      'Unallowed $update option.  Allowed Values: -1 to 5',0,0)
       EndIf
       IOP(15) = IUpDat
cc
      ELSE IF(CHAR(1:8).EQ.'$optcycl') THEN
cc
       READ(CHAR(11:80),*) MaxCYC
       If(MaxCYC.LT.1) Then
        Call nerror(18,'OPTIMIZE module',
     $   'Maximum Optimization Cycles must be a Positive Integer',0,0)
       EndIf
       IOP(18) = MaxCYC
cc
      ELSE IF(CHAR(1:8).EQ.'$regener') THEN
cc
       READ(CHAR(13:80),*) IGen
       If(IGen.LT.0.OR.IGen.GT.2) Then
        Call nerror(2,'OPTIMIZE module',
     $   'Unallowed $regenerate option.  Allowed Values: 0, 1 or 2',0,0)
       EndIf
       IOP(2) = IGen
cc
      ELSE IF(CHAR(1:8).EQ.'$gradter') THEN
cc
       IOP(17) = 1
cc
      ELSE IF(CHAR(1:5).EQ.'$ctol') THEN
cc
       READ(CHAR(7:80),*) CTol
       If(CTol.LT.1.0d-6) Then
        Call nerror(19,'OPTIMIZE module',
     $    'Constraint Criterion $ctol must not be less than 10**-6',0,0)
       ENDIF
       IOP(19) = NINT(CTol*10**6)
cc
      ELSE IF(CHAR(1:7).EQ.'$linear') THEN
cc
       READ(CHAR(9:80),*) Linear
       If(Linear.LT.150.OR.Linear.GT.180) Then
        Call nerror(20,'OPTIMIZE module',
     $   'Near-Linear Bend $linear must be between 150 and 180 degrees',
     $    0,0)
       EndIf
       IOP(20) = Linear
cc
      ELSE IF(CHAR(1:6).EQ.'$bskal') THEN
cc
       READ(CHAR(9:80),*) BSkal
       If(BSkal.LE.Zero) Then
        Call nerror(21,'OPTIMIZE module',
     $   'Scaling Factor for Inverse-Distances must be positive',0,0)
       EndIf
       IOP(21) = NINT(BSkal*10**3)
cc
      ELSE IF(CHAR(1:7).EQ.'$notors') THEN
cc
       IOP(22) = 1
cc
      ELSE IF(CHAR(1:5).EQ.'$qmmm') THEN
cc
       IOP(23) = 1
cc
      ELSE IF(CHAR(1:8).EQ.'$backtra') THEN
cc
       READ(CHAR(12:80),*) IBack
       If(IBack.LT.0.AND.IBack.GT.1) Then
        Call nerror(2,'OPTIMIZE module',
     $   'Unallowed $backtrans option.  Allowed Values: 0 or 1',
     $    0,0)
       EndIf
       IOP(25) = IBack
cc
      ELSE IF(CHAR(1:7).EQ.'$cutoff') THEN
cc
       READ(CHAR(9:80),*) CutOff
       If(CutOff.LE.Zero) Then
        Call nerror(22,'OPTIMIZE module',
     $   'Bonding Distance Cutoff must be positive',0,0)
       EndIf
       IOP(26) = NINT(CutOff*10**6)
cc
      ELSE IF(CHAR(1:7).EQ.'$hcnvrt') THEN
cc
       IOP(27) = 1
cc
      ELSE IF(CHAR(1:8).EQ.'$coordty') THEN
cc
       READ(CHAR(12:80),*) IType
       If(IType.LT.0.AND.IType.GT.2) Then
        Call nerror(2,'OPTIMIZE module',
     $   'Unallowed $coordtype option.  Allowed Values: 0, 1 or 2',
     $    0,0)
       EndIf
       IOP(28) = IType
cc
      ELSE IF(CHAR(1:6).EQ.'$print') THEN
cc
       READ(CHAR(8:80),*) IPRNT
       IOP(29) = IPRNT
cc
      ELSE IF(CHAR(1:4).EQ.'$fix') THEN
cc
 20    CONTINUE
       READ(40,900) CHAR
       If(CHAR(1:7).EQ.'$endfix') GO TO 100
       GO TO 20
cc
      ELSE IF(CHAR(1:6).EQ.'$dummy') THEN
cc
 30    CONTINUE
       READ(40,900) CHAR
       If(CHAR(1:7).EQ.'$enddum') GO TO 100
       GO TO 30
cc
      ELSE IF(CHAR(1:8).EQ.'$constra') THEN
cc
 40    CONTINUE
       READ(40,900) CHAR
       If(CHAR(1:8).EQ.'$endcons') GO TO 100
       GO TO 40
cc
      ELSE IF(CHAR(1:8).EQ.'$connect') THEN
cc
 50    CONTINUE
       READ(40,900) CHAR
       If(CHAR(1:8).EQ.'$endconn') GO TO 100
       GO TO 50
cc
      ELSE IF(CHAR(1:8).EQ.'$surface') THEN
cc
 60    CONTINUE
       READ(40,900) CHAR
       If(CHAR(1:8).EQ.'$endsurf') GO TO 100
       GO TO 60
cc
      ELSE IF(CHAR(1:6).EQ.'$drive') THEN
cc
 70    CONTINUE
       READ(40,900) CHAR
       If(CHAR(1:8).EQ.'$enddriv') GO TO 100
       GO TO 70
cc
      ELSE
cc
       call rmblan(CHAR,80,lchar)
       If(lchar.eq.0) GO TO 100
       CHAR = 'Unknown keyword in <opt> file: '//CHAR(1:8)
       Call nerror(50,'OPTIMIZE module',CHAR,0,0)
cc
      ENDIF
      GO TO 100
C
 95   CONTINUE
      RETURN         ! no OPT section found
C
 96   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  check CutOff defaults for surface/cluster optimization
C
      If(IType.EQ.1.AND.CutOff.EQ.Zero) IOP(26) = 3000000
      If(IType.EQ.2.AND.CutOff.EQ.Zero) IOP(26) = 5000000
C
      RETURN
c
  900 Format(A80)
  901 Format(A20)
c
      END
c =======================================================================
c
      SUBROUTINE RdSurf(NAtoms,NMol,IMOL,IPRNT,NFix,IFIX)
      IMPLICIT INTEGER(A-Z)
C
C
C  Attempt to read surface constraint data from OPT File.
C  This will be a list of surface atoms which will be fixed
C  during the optimization.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NMol    -  number of molecules (should be two in this implementation)
C  IMOL    -  pointers to start/end of molecules in XC array
C  IPRNT   -  print flag
C  NFix    -  on input number of lines of fixed atoms
C             on exit number of fixed atoms
C  IFIX    -  list of fixed atoms
C             (for compatibility with Cartesian optimization,
C              each atom with have X, Y and Z components fixed)
C
      PARAMETER (MaxC=10)
C
      INTEGER IMOL(NMol+1),IFix(3,NAtoms),IAC(MaxC)
      CHARACTER CHAR*80,jobname*256
c
      Common /job/jobname,lenJ
C
C
      IOut = ioutfil('iout')
c
      CALL IZeroIT(IFIX,3*NAtoms)
      MAtoms = IMol(2)    ! number of atoms in surface
C
C  open OPT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.opt',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $opt section
C
 10   CONTINUE
      READ(40,900,END=95) CHAR
      IF(CHAR(1:4).NE.'$opt') GO TO 10
C
C  locate start of surface section
C
 20   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:8).NE.'$surface') GO TO 20
C
C  surface data found
C
      If(IPRNT.GT.1) WRITE(IOut,1000)
C
C  EXPECTED FORMAT
C   FIXED ALL
C   FIXED LIST   (followed by list of fixed atoms, 10 per line)
C
      READ(40,900) CHAR
      IF(CHAR(1:9).EQ.'fixed all') THEN
       NFix = MAtoms
       DO 30 I=1,MAtoms
       IFix(1,I) = 1
       IFix(2,I) = 1
       IFix(3,I) = 1
 30    CONTINUE
cc
      ELSE
C
C  read NFix lines
C
       BACKSPACE 40
c
       DO 60 I=1,NFix
c
       CALL IZeroIT(IAC,MaxC)
       READ(40,910,ERR=97) (IAC(L),L=1,MaxC)
C
C  determine how many fixed atoms there really are
C
       DO 50 J=1,MaxC
       JJ = IAC(J)
       If(JJ.EQ.0) GO TO 60
c
       IF(JJ.GT.MAtoms) THEN
        WRITE(IOut,1100) JJ,MAtoms
        CALL OptExit(9)
       EndIF
c
       IFIX(1,JJ) = 1
       IFIX(2,JJ) = 1
       IFIX(3,JJ) = 1
 50    CONTINUE
c
 60    CONTINUE
cc
      ENDIF
C
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
C
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no OPT file found
      WRITE(IOut,1300)
      CALL OptExit(9)
c
 96   CONTINUE          ! no surface data in OPT file
      WRITE(IOut,1400)
      CALL OptExit(9)
c
 97   CONTINUE          ! error reading surface data
      WRITE(IOut,1500)
      CALL OptExit(9)
c
  900 Format(A80)
  910 Format(10X,10I4)
 1000 FORMAT(' Surface Constraint Data Read from OPT file')
 1100 FORMAT(/,2X,'***ERROR*** Fixed Atom ',I4,' of surface data',
     $            ' is greater than ',I3,/,'   the total number of',
     $            ' surface atoms. Check your input data')
 1300 FORMAT(/,2X,'***ERROR*** There should be surface data',
     $            ' but no OPT file found!')
 1400 FORMAT(/,2X,'***ERROR*** There should be surface data',
     $            ' but none found in OPT file!')
 1500 FORMAT(/,2X,'***ERROR*** Problem Reading Surface Data',/,
     $         2X,'   Check your OPT file')
c
      END
c =======================================================================
c
      SUBROUTINE RedTOPG(NIC,NVar,IG,GINT,ITMP,GVAR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reduce the full internal coordinate gradient to the set of
C  parameters to be optimized (Z-Matrix variables)
C  (allowing for possible over-constrained optimization)
C
C  ARGUMENTS
C
C  NIC     -  total number of internal coordinates
C  NVar    -  number of Z-Matrix variables
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  GINT    -  full vector of internal coordinate gradients
C  ITMP    -  integer scratch array
C  GVAR    -  on exit contains Z-matrix gradient
C
C
      DIMENSION IG(NIC),GINT(NIC),ITMP(NIC),GVAR(NVar)
C
C
      CALL IZeroIT(ITMP,NIC)
      NV = 0
c
      DO 10 I=1,NIC
      IT = IAbs(IG(I))
      IF(IT.EQ.0) THEN
       NV = NV+1
       ITMP(I) = NV
       GVAR(NV) = GINT(I)
      ELSE IF(IT.NE.1000) THEN
       JT = ITMP(IT)
       GVAR(JT) = GVAR(JT) + GINT(I)*ISIGN(1,IG(I))
      ENDIF
 10   CONTINUE
C
C  Check on number of variables
C
      IF(NV.NE.NVar) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Wrong Number of Variables in <RedTOPG>')
c
      END
c =======================================================================
c
      SUBROUTINE RedTOPX(NIC,NVar,IG,ZINT,XINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reduce the full internal coordinate set to the set of
C  parameters to be optimized (Z-Matrix variables)
C
C  ARGUMENTS
C
C  NIC     -  total number of internal coordinates
C  NVar    -  number of Z-Matrix variables
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  ZINT    -  full vector of internal coordinates
C  XINT    -  on exit contains Z-matrix variables
C
C
      DIMENSION IG(NIC),ZINT(NIC),XINT(NVar)
C
C
      NV = 0
c
      DO 10 I=1,NIC
      IF(IG(I).EQ.0) THEN
       NV = NV+1
       XINT(NV) = ZINT(I)
      ENDIF
 10   CONTINUE
C
C  Check on number of variables
C
      IF(NV.NE.NVar) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Wrong Number of Variables in <RedTOPX>')
c
      END
c =======================================================================
c
      SUBROUTINE ReOrderZ(NZ,INDX,XZ,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reorders Cartesian coordinates from Z-matrix
C  back-transformation to original atom order
C
C  ARGUMENTS
C
C  NZ      -  number of centers
C  INDX    -  index reordering vector
C  XZ      -  initial Cartesian coordinates
C  XC      -  reordered Cartesians
C
C
      DIMENSION INDX(NZ),XZ(3,NZ),XC(3,NZ)
C
C
      DO 10 I=1,NZ
      II = INDX(I)
      XC(1,II) = XZ(1,I)
      XC(2,II) = XZ(2,I)
      XC(3,II) = XZ(3,I)
 10   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE rij2cart( n,  nrij,  rijindx, rijlst, M,
     $                     EVals, EVecs, work, Iwork,IFail,
     $                     xyzlst )
      IMPLICIT NONE
c	..Scalar Arguments..
      INTEGER n, nrij
c	..Array Arguments

      INTEGER rijindx(4,nrij)
      DOUBLE PRECISION rijlst(nrij), xyzlst(3,n)
c =====================================================================
c
c Donald B. Kinghorn
c University of Arkansas
c May 14 1998
c
c Last Modified
c =============
c May 21 1998 DBK
c
c =====================================================================
c Purpose
c =======
c
c Given a list of distances rijlst and an index array rijindx
c compute a set of cartesian coordinates that are consistent
c with the distances rij
c
c The program computes the metric tensor M where:
c
c     mij = xi.xj = ri0 rj0 cos(a) = 1/2 (ri0^2 + rj0^2 - rij^2)
c and
c     ri0 = 1/n Sum(k->n)[rik^2] - 1/n^2 Sum(j<k)[rjk^2]
c
c [di0 = the distance from the centroid of the structure]
c The eigen decomposition of M is then found M = W E W'
c The cartesian coordinates are then given using the 3 largest
c eigenvalues and corresponding vectors
c
c     x(k,i) = e(k)^1/2 w(i,k)
c
c You really should have all n(n-1)/2 distances dij to get a
c meaningful result. However the method is relatively insensitive to
c small errors in the distances and will give a "best fit" to the
c cartesians.
c
c Reference
c =========
c See Aszodi and Taylor "Computers Chem." Vol.21, No.1, pp13-23 (1997)
c and references therin for a discussion of the method
c
c Arguments
c =========
c
c INPUT:
c
c  n        -  (int) number of atoms [ dim of M ]
c
c  nrij     -  (int) number of distances rij
c                    NOTE: this should be n(n-1)/2 for best results
c
c  rijindx  -  (int 4xnrij Array) index set for distances
c				i j 0 0  for distances
c				[ the last 2 columns of rijindx are not used ]
c
c  rijlst   -  (dbl nrij Array)  values for the distances
c
c OUTPUT:
c
c  xyzlst   -  (dbl 3xn Array)  cartesion coordinates
c
c =====================================================================
c	..
c     ..Parameters..
      DOUBLE PRECISION    zero, onehalf
      PARAMETER           (zero = 0.0d0,
     &                     onehalf = 0.5d0)

c     ..
c     ..Local Scalars..
      INTEGER  i, j, k
      DOUBLE PRECISION    sum1, sum2, lx, ly, lz
c     ..
c     ..Local Arrays..
      DOUBLE PRECISION    M(n,n)
c     ..
c     ..Logicals


c     ..
c     ..Local scalars and Arrays needed for the eigen code..
      INTEGER             INFO, IWORK(5*n), nEV,
     &                    LWORK, IFAIL(n), IL, IU

      DOUBLE PRECISION    WORK(32*n), ABSTOL,
     &                    VL, VU, EVALS(n), EVECS(n,n)

c     ..
c     .. External Subroutines ..
      EXTERNAL            DSYEVX


c     ..
c     ..Executable Statements..
c
c**********************************************************************
c
c	Build the metric matrix M
c
c**********************************************************************

c     first fill off diag of M with squared distances rij^2
      DO k=1,nrij
          i = rijindx(1,k)
          j = rijindx(2,k)
          M(i,j) = rijlst(k)*rijlst(k)
          M(j,i) = M(i,j)
      END DO

c     now put the squared distance of each point from the centroid
c     on the diagonal of M

      sum1 = zero
      DO j=1,n-1
          DO i=j+1,n
              sum1 = sum1 + M(i,j)
          END DO
      END DO
      sum1 = sum1/(n*n)

      DO i=1,n
          sum2 = zero
          M(i,i) = zero
          DO j=1,n
              sum2 = sum2 + M(j,i)
          END DO
          M(i,i) = sum2/n - sum1
      END DO

c     mij = xi.xj = 1/2 (ri0^2 + rj0^2 - rij^2)
      DO j=1,n-1
          DO i=j+1,n
              M(i,j) = onehalf * ( M(i,i) + M(j,j) - M(i,j) )
              M(j,i) = M(i,j)
          END DO
      END DO

c**********************************************************************
c
c Find the largest 3 eigenvalues and vectors of M
c
c**********************************************************************

c     using calls to LAPACK routines for eigenvalues and vectors

c     set ABSTOL to something reasonable
c     (machine safe min gives most accurate result)
      ABSTOL = 2.226d-308
      IL = n-2
      IU = n
c      VL = -1.0d0
c      VU = 10.0d0
      LWORK = 32*n

      CALL DSYEVX( 'V','I', 'L', n, M, n, VL, VU, IL, IU, ABSTOL, nEV,
     &             EVALS, EVECS, n, WORK, LWORK, IWORK, IFAIL, INFO )

c**********************************************************************
c
c Generate the cartesian coordinates and load them into xyzlst
c
c**********************************************************************

      lz = DSQRT(EVALS(1))
      ly = DSQRT(EVALS(2))
      lx = DSQRT(EVALS(3))

      DO j=1,n
          xyzlst(1,j) = lx * EVECS(j,3)
          xyzlst(2,j) = ly * EVECS(j,2)
          xyzlst(3,j) = lz * EVECS(j,1)
      END DO
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE RmFIX(N1,     IFIX,   GOld,   HOld,   Gflag,
     $                 N2,     GNew,   HNew)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Removes fixed Cartesian coordinates from the parameter
C  set prior to calculating the optimization step
C
C  ARGUMENTS
C
C  N1      -  dimension of full space
C  IFix    -  array indicating optimization status for each coordinate
C              0 - optimize
C              1 - fixed (remove)
C  GOld    -  full cartesian gradient
C  HOld    -  full cartesian Hessian
C  Gflag   -  logical flag   .true. - only handle gradient
C                           .false. - gradient and Hessian
C
C  on exit ...
C
C  N2      -  dimension after fixed coordinates removed
C  GNew    -  gradient with fixed coordinates removed
C  HNew    -  Hessian with fixed coordinates removed
C
      DIMENSION IFix(N1),GOld(N1),HOld(N1,N1),
     $          GNew(N2),HNew(N2,N2)
      LOGICAL Gflag
C
C
C  remove fixed coordinates from gradient
C
      IT = 0
      JT = 0
      DO 10 I=1,N1
      IF(IFix(I).EQ.0) THEN
       IT = IT+1
       JT = JT+1
       GNew(IT) = GOld(JT)
      ELSE IF(IFix(I).EQ.1) THEN
       JT = JT+1
      ENDIF
 10   CONTINUE
c
      If(Gflag) RETURN
C
C  remove fixed coordinates from Hessian
C
      IT = 0
      JT = 0
      DO 30 I=1,N1
      IF(IFix(I).EQ.0) THEN
       IT = IT+1
       JT = JT+1
       KT = 0
       LT = 0
       DO 20 J=1,I
       IF(IFix(J).EQ.0) THEN
        KT = KT+1
        LT = LT+1
        HNew(IT,KT) = HOld(JT,LT)
        HNew(KT,IT) = HNew(IT,KT)
       ELSE IF(IFix(J).EQ.1) THEN
        LT = LT+1
       ENDIF
 20    CONTINUE
      ELSE IF(IFix(I).EQ.1) THEN
       JT = JT+1
      ENDIF
 30   CONTINUE
C
C  Check that correct number of elements set
C
      IF(IT.NE.N2.AND.JT.NE.N1) THEN
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       CALL OptExit(9)
      ENDIF
c
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Incorrect number of Fixed',
     $            ' Coordinates in <RmFIX>')
c
      END
c =======================================================================
c
      SUBROUTINE RotHES(NAtoms,RM,SCR,HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8 RM(3,3),SCR(3,3),HESS(3*NAtoms,3*NAtoms)
C
C  Transforms the Hessian in Cartesian coordinates according
C  to the rotation matrix, RM
C  This is done one (3x3) block at a time
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  RM      -  rotation matrix
C  SCR     -  scratch storage
C  HESS    -  Hessian matrix
C             on exit contains transformed (rotated) Hessian
C
C
      DO 30 IAtm=1,NAtoms
      II = (IAtm-1)*3
      DO 30 JAtm=1,IAtm
      JJ = (JAtm-1)*3
C
C  transform the (IAtM,JAtM) block
C
      DO 10 I=1,3
      DO 10 J=1,3
      JT = JJ+J
      SCR(I,J) = RM(I,1)*HESS(II+1,JT) + RM(I,2)*HESS(II+2,JT)
     $                 + RM(I,3)*HESS(II+3,JT)
 10   CONTINUE
c
      DO 20 I=1,3
      IT = II+I
      DO 20 J=1,3
      JT = JJ+J
      HESS(IT,JT) = SCR(I,1)*RM(J,1) + SCR(I,2)*RM(J,2)
     $                     + SCR(I,3)*RM(J,3)
 20   CONTINUE
C
 30   CONTINUE
C
C  Form the full Hessian
C  (since we only transformed the lower "triangle")
C
      DO 50 IAtm=2,NAtoms
      II = (IAtm-1)*3
      DO 50 JAtm=1,IAtm-1
      JJ = (JAtm-1)*3
c
      DO 40 I=1,3
      IT = II+I
      DO 40 J=1,3
      JT = JJ+J
      HESS(JT,IT) = HESS(IT,JT)
 40   CONTINUE
 50   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE RotVEC(NAtoms,RM,V)
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8 RM(3,3),V(3,NAtoms)
C
C  Transforms a "coordinate" vector V according to
C  the rotation matrix, RM
C  On exit V contains the transformed vector
C
C
      DO  10 I=1,NAtoms
      V1 = RM(1,1)*V(1,I) + RM(1,2)*V(2,I) + RM(1,3)*V(3,I)
      V2 = RM(2,1)*V(1,I) + RM(2,2)*V(2,I) + RM(2,3)*V(3,I)
      V3 = RM(3,1)*V(1,I) + RM(3,2)*V(2,I) + RM(3,3)*V(3,I)
      V(1,I) = V1
      V(2,I) = V2
      V(3,I) = V3
 10   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE RstFIX(NAtoms, NDum,   NDim,   IFix,   DOld,
     $                  GOld,   DNew,   GNew)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Restores fixed Cartesian coordinates to the displacement and
C  gradient vectors and transforms to centre of mass coordinates
C
C  ARGUMENTS
C
C  NAtoms  -  number of real atoms
C  NDum    -  number of dummy atoms
C  NDim    -  total dimension (including constraints)
C  IFix    -  array indicating optimization status for each coordinate
C              0 - optimize
C              1 - fixed (remove)
C  DOld    -  displacement vector with fixed coordinates removed
C  GOld    -  gradient vector with fixed coordinates removed
C  DNew    -  full cartesian displacement vector
C  GNew    -  full cartesian gradient vector
C
      DIMENSION IFix(NDim),GOld(*),GNew(NDim),DOld(*),DNew(NDim)
C
      PARAMETER (Zero=0.0d0)
C
C
      IT = 0
      JT = 0
      DO 10 I=1,NDim
      IF(IFix(I).EQ.0) THEN
       IT = IT+1
       JT = JT+1
       DNew(IT) = DOld(JT)
       GNew(IT) = GOld(JT)
      ELSE IF(IFix(I).EQ.1) THEN
       IT = IT+1
       DNew(IT) = Zero
       GNew(IT) = Zero
      ENDIF
 10   CONTINUE
C
C  now put displacement vector in CMS frame
C  (dummy atoms have zero mass)
C
      CALL CMS(NAtoms,X,Y,Z,DNew)
c
      DO 20 IATM=NAtoms+1,NAtoms+NDum
      II = 3*(IATM-1)
      DNew(II+1) = DNew(II+1) - X
      DNew(II+2) = DNew(II+2) - Y
      DNew(II+3) = DNew(II+3) - Z
 20   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE S2Bend(I,J,K,L,XC,SS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Determines "sign" (direction) of bond I-J with respect
C  to plane J-K-L for two angle Z-matrix definition
C  If value is too small, atoms are nearly planar and
C  this choice of bends should be rejected
C
C  ARGUMENT
C
C  I       -  first atom in Z matrix
C  J       -  second atom in Z matrix (bonded to I)
C  K       -  third atom in Z matrix (making angle I-J-K)
C  L       -  fourth atom in Z matrix (making angle I-J-L)
C  XC      -  Cartesian coordinates
C  SS      -  on exit sign and magnitude of vector product JKxJL
C
C
      DIMENSION XC(3,*)
      Dimension U(3),V(3),W(3)
C
C  define vectors J-K and J-L
C
      U(1) = XC(1,K) - XC(1,J)
      U(2) = XC(2,K) - XC(2,J)
      U(3) = XC(3,K) - XC(3,J)
      V(1) = XC(1,L) - XC(1,J)
      V(2) = XC(2,L) - XC(2,J)
      V(3) = XC(3,L) - XC(3,J)
      W(1) = U(2)*V(3)-U(3)*V(2)
      W(2) = U(3)*V(1)-U(1)*V(3)
      W(3) = U(1)*V(2)-U(2)*V(1)
      SS = SProd(3,W,W)
C
C  take dot product with J-I
C
      U(1) = XC(1,I) - XC(1,J)
      U(2) = XC(2,I) - XC(2,J)
      U(3) = XC(3,I) - XC(3,J)
      SS = SProd(3,U,W)/SS
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE S_BTRAN(NAtoms, XC,     NPrim,  ktyp,   klist,
     $                   SavTOR, makeq,  makeb,  itor,   IPRNT,
     $                   NDEG,   NCmp,   NP1,    INT1,   UT,
     $                   IBM,    BM,     nb,     NB1,    INB1,
     $                   B,      NB2,    INB2,   BT,     XPrim,  IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Direct construction of B-Matrix in natural internal coordinates
C  [primitive B matrix is NOT explicitly formed - the non-zero elements
C   of each column (maximum 12) multiply non-zero elements of UT
C   directly to form B(UT)]     SPARSE MATRIX ALGORITHM
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  NPrim   -  total number of primitive internal coordinates
C  ktyp    -  integer array containing internal coordinate type
C  klist   -  list of atoms involved in each primitive
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  makeq   -  integer flag for determining current value of
C             primitive internal coordinate
C               0 - NO;  1 - YES
C  makeb   -  integer flag for calculating B-Matrix column
C               0 - NO;  1 - YES
C  itor    -  integer flag for saving/checking values of primitive
C             torsions (this is done due to possible convergence
C             problems in iterative generation of new Cartesians)
C              0 - no action
C              1 - save initial primitive torsions in SavTOR
C              2 - check current torsions against initial
C                  (may need to change sign of angles > PI)
C  IPRNT   -  flag controlling print out
C  NDEG    -  number of active delocalized internal coordinates
C  NCmp    -  total number of non-zero natural internal components
C  NP1     -  number of primitives contributing to each NIC
C  INT1    -  indices for all non-zero primitives for each NIC
C  UT      -  non-zero NIC components ordered per NIC
C  ................................................................
C  IBM     -  work array for non-zero B-matrix indices
C  BM      -  work array for non-zero B-matrix entries
C             ** NOTE - Assumed that each B-matrix row contains
C                       at most 52 non-zero entries
C  ................................................................
C
C  on exit
C
C  nb      -  total number of non-zero B-matrix elements
C  NB1     -  number of non-zero entries for each row
C  INB1    -  indices for all non-zero B-matrix rows
C  B       -  non-zero B-matrix elements ordered per row
C  NB2     -  number of non-zero entries for each column
C  INB2    -  indices for all non-zero B-matrix columns
C  BT      -  non-zero B-matrix elements ordered per column
C  XPrim   -  primitive internal coordinate values
C  IErr    -  error flag   0 - success
C                         -1 - something went wrong
C
C
      DIMENSION XC(3,NAtoms),ktyp(NPrim),klist(4,NPrim),SavTOR(NPrim),
     $          UT(NCmp),XPrim(NPrim)
      DIMENSION B(12*NDEG),BT(12*NDEG),BM(3*NAtoms,52)
      INTEGER NP1(NDEG),INT1(NCmp),IBM(3*NAtoms,52)
      INTEGER NB1(3*NAtoms),INB1(12*NDEG),NB2(NDEG),INB2(12*NDEG)
      Dimension INDX(12),BVal(12),kval(7)
c -------------------------------------------
      dimension BB(NDEG,3*Natoms)
c -------------------------------------------
C
      PARAMETER (Zero=0.0d0,TolZero=1.0d-9)
C
      Data kval/ 6, 9, 12, 12, 12, 12, 6 /
C
C
C  initialize
C
      nc = 0
      nbb = 0
C
C  Loop over all natural internals
C
      DO 60 IC=1,NDEG
      nb = 0
C
C  Loop over all primitives in this NIC
C
      DO 20 I=1,NP1(IC)
      nc = nc+1
      IP = INT1(nc)
      UVal = UT(nc)
C
C  construct primitive B-Matrix column for coordinate IC
C
      CALL BMAT2(IP,     NAtoms, XC,     NPrim,  ktyp,
     $           klist,  SavTOR, makeq,  makeb,  itor,
     $           BSkal,  IPRNT,  INDX,   BVal,   QQ,  IErr)
c
      If(IERr.NE.0) CALL OptExit(9)
      If(makeq.eq.1) XPRIM(IP) = QQ
C
C  on exit from <BMAT2>
C     INDX - index of non-zero entries to B-Matrix column
C     BVal - non-zero B Matrix values  (maximum 12)
C
      nval = kval(ktyp(IP))
      DO 10 J=1,nval
      JJ = INDX(J)
      BB(IC,JJ) = BB(IC,JJ) + UVal*BVal(J)
      nb = nb + 1                    ! increment counter for row JJ of B
      INB1(nb) = JJ                  ! store column (Cart) index
      B(nb) = UVal*BVal(J)           ! store column (Cart) value
 10   CONTINUE
 20   CONTINUE
C
C  finished internal coordinate IC
C  Sum up any duplications among the 3*NAtoms contributions to each
C  B matrix row
C
      DO I=1,3*NAtoms
      NB1(I) = 0                ! NB1 used as scratch
      EndDO
      DO 40 J=1,nb              ! number of non-zero contributions to row IC
      JJ = INB1(J)              ! index of contribution
      IF(NB1(JJ).NE.0) THEN
C
C  index is duplicate - sum up
C
       JJJ = NB1(JJ)            ! original position of index
       INB1(J) = 0              ! mark duplicate index for removal
       B(JJJ) = B(JJJ) + B(J)
      ELSE
C
C  store index
C
       NB1(JJ) = J
      ENDIF
 40   CONTINUE
C
C  now remove duplications
C
      ns = nbb
      DO 45 J=1,nb
      IF(INB1(J).NE.0.AND.Abs(B(J)).GT.TolZero) THEN
       nbb = nbb+1
       INB2(nbb) = INB1(J)
       BT(nbb) = B(J)
      ENDIF
 45   CONTINUE
      NB2(IC) = nbb-ns          ! # non-zero B-matrix elements for row IC
cc
 60   CONTINUE
C
C  At this point we have non-zero B-matrix elements stored for
C  each natural internal coordinate, i.e., B-transpose
C  Resort to form non-zero elements for each Cartesian component
C
      Do I=1,3*NAtoms
      NB1(I) = 0
      EndDo
c
      nb = 0
      DO 70 IC=1,NDEG
      DO 70 I=1,NB2(IC)
      nb = nb+1
      II = INB2(nb)
      NB1(II) = NB1(II) + 1
      IBM(II,NB1(II)) = IC
      BM(II,NB1(II)) = BT(nb)
 70   CONTINUE
c
      nb = 0
      DO 80 I=1,3*NAtoms
      DO 80 J=1,NB1(I)
      nb = nb+1
      INB1(nb) = IBM(I,J)
      B(nb) = BM(I,J)
 80   CONTINUE
c
      If(makeb.EQ.1.AND.IPRNT.GT.6) Then
       write(6,*) ' Print Out of non-zero elements from full B'
       nbb = 0
       do i=1,ndeg
       write(6,*) ' non-zero entries for NIC:',i
       do j=1,3*natoms
       if(abs(bb(i,j)).gt.tolzero) then
        write(6,*) j,bb(i,j)
        nbb = nbb+1
       endif
       enddo
       enddo
       write(6,*) ' TOTAL NUMBER OF NON-ZERO ELEMENTS IS:',nbb
       write(6,*) ' number of non-zero entries per row:'
       do i=1,ndeg
       nb=0
       do j=1,3*natoms
       if(abs(bb(i,j)).gt.tolzero) then
        nb = nb+1
        write(6,*) i,j,bb(i,j)
       endif
       enddo
       write(6,'(1x,i6,2x,i6)') i,nb
       enddo
      EndIf
C
      RETURN
c
 1000 FORMAT(/,1x,I9,' Non-Zero B-Matrix Elements',/)
c
      END
c =======================================================================
c
      SUBROUTINE S_CONJUGATE(NAT3,   N,      Nb,     NB1,    INB1,
     $                       B,      NB2,    INB2,   BT,     C,
     $                       D,      X,      thrsh,  MaxIT,  IPRNT,
     $                       PV,     RV,     YV,     ZV,     ST,
     $                       IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine solves the linear equation set  (B*B(t))X = C by the
C  method of conjugate gradients
C  Adapted from PP routine <CONGRAD> used in sparse
C  matrix transformation
C  ** SPARSE MATRIX ALGORITHM **
C
C  ARGUMENTS
C
C  NAT3    -  dimension of B-matrix (3*NAtoms)
C  N       -  number of active natural internal coordinates
C  Nb      -  total number of non-zero B-matrix elements
C  NB1     -  number of non-zero entries for each row
C  INB1    -  indices for all non-zero B-matrix rows
C  B       -  non-zero B-matrix elements ordered per row
C  NB2     -  number of non-zero entries for each column
C  INB2    -  indices for all non-zero B-matrix columns
C  BT      -  non-zero B-matrix elements ordered per column
C  C       -  right-hand side vector  (B*gcart)
C  D       -  inverse diagonals of A matrix
C  X       -  on entry initial guess for solution
C             on exit  solution (if converged)
C  thrsh   -  convergence criterion
C  MaxIT   -  maximum number of cycles allowed
C  IPRNT   -  print flag
C  PV      -  scratch vector
C  RV      -   ditto
C  YV      -   ditto
C  ZV      -   ditto
C  ST      -   ditto
C  IErr    -  error flag    0 - converged
C                          -1 - failed
C
C
      DIMENSION B(Nb),BT(Nb),C(N),D(N),X(N),
     $          PV(N),RV(N),YV(N),ZV(N),ST(NAT3)
      INTEGER NB1(NAT3),INB1(Nb),NB2(N),INB2(Nb)
      Dimension cc(2)
C
      PARAMETER (Zero=0.0d0,thrs=1.0d-30)
C
C
C  initialize circular counters
C
      IOut = ioutfil('iout')
      IErr = 0
c
      i0 = 1
      i1 = 2
C
C  set YV = BT*(B*X(0))
C
      np = 0
      DO 10 I=1,NAT3
      ST(I) = Zero
      DO 9 J=1,NB1(I)
      np = np+1
      ST(I) = ST(I) + B(np)*X(INB1(np))
  9   CONTINUE
 10   CONTINUE
c
      np = 0
      DO 12 I=1,N
      YV(I) = Zero
      DO 11 J=1,NB2(I)
      np = np+1
      YV(I) = YV(I) + BT(np)*ST(INB2(np))
 11   CONTINUE
 12   CONTINUE
C
C  set RV = C - BT*B*X(0)
C      ZV = D*RV
C      PV = ZV
C      cc(1) = RV(t)*ZV
C
      Sum = Zero
      DO 20 I=1,N
      RV(I) = C(I) - YV(I)
      ZV(I) = D(I)*RV(I)
      PV(I) = ZV(I)
      Sum = Sum + RV(I)*ZV(I)
 20   CONTINUE
C
C  if sum is very small already means initial gradient
C  is essentially exact
C
      If(Abs(Sum).LT.thrs) RETURN
C
C  start of iterative loop
C
      cc(i0) = Sum
      k = 0
c
 100  CONTINUE
      k = k+1
C
C  set YV = BT*(B*PV(k-1))
C
      np = 0
      DO 30 I=1,NAT3
      ST(I) = Zero
      DO 29 J=1,NB1(I)
      np = np+1
      ST(I) = ST(I) + B(np)*PV(INB1(np))
 29   CONTINUE
 30   CONTINUE
c
      np = 0
      DO 32 I=1,N
      YV(I) = Zero
      DO 31 J=1,NB2(I)
      np = np+1
      YV(I) = YV(I) + BT(np)*ST(INB2(np))
 31   CONTINUE
 32   CONTINUE
c
      Sum = SProd(N,PV,YV)
      bb = cc(i0)/Sum
C
C  X(k)  = X(k-1)  + bb*PV(k-1)
C  RV(k) = RV(k-1) - bb*YV
C
      DO 40 I=1,N
      X(I) = X(I) + bb*PV(I)
      RV(I) = RV(I) - bb*YV(I)
 40   CONTINUE
C
C  check for convergence
C  Absolute value of RV < thrsh
C
      DMax = Zero
      QSQR = Zero
C
      DO 50 I=1,N
      DX = RV(I)
      If(Abs(DX).GT.DMax) DMax = Abs(DX)
      QSQR = QSQR + DX**2
 50   CONTINUE
C
      QSQR = SQRT(QSQR/N)
c
      If(IPRNT.GT.2) WRITE(IOut,1200) k,DMax,QSQR
C
C  converged?
C
cc      If(DMax.LT.thrsh) RETURN
      If(DMax.LT.thrsh) Then
       write(IOut,'('' Gradient converged in'',i4,'' cycles'')')k
       return
      endif
c
      If(k.GT.MaxIT) Then
        WRITE(IOut,1300)
        IErr = -1
        RETURN
      EndIf
C
C  prepare for next cycle
C
      DO 60 I=1,N
      ZV(I) = RV(I)*D(I)
 60   CONTINUE
c
      Sum = SProd(N,ZV,RV)
      cc(i1) = Sum
      bb = cc(i1)/cc(i0)
C
C  PV(k) = ZV + bb*PV(k-1)
C
      DO 70 I=1,N
      PV(I) = ZV(I) + bb*PV(I)
 70   CONTINUE
c
      ii = i1
      i1 = i0
      i0 = ii
      GO TO 100
c
 1200 FORMAT(5X,'Cycle: ',I3,'  Maximum deviation: ',F12.8,
     $          '  RMS deviation: ',F12.8)
 1300 FORMAT(/,2X,'***ERROR*** Exceeded allowed number of iterative',
     $            ' cycles in <S_CONJUGATE>')
c
      END
c =======================================================================
c
      SUBROUTINE S_DefHES(NDEG,   intcor, NCon,   ktyp,   klist,
     $                    XC,     IAN,    FC,     NCmp,   NP1,
     $                    INT1,   UT,     HINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate Hessian by forming UT(t)*H(prim)*UT
C  where H(prim) is a "Hessian" diagonal matrix in the
C  primitive space with appropriate force constants
C  given for each primitive type
C  ** SPARSE MATRIX ALGORITHM **
C
C  ARGUMENTS
C
C  NDEG    -  number of degrees of freedom
C             (and number of natural internal coordinates)
C  intcor  -  number of primitive internals
C  NCon    -  number of constraints
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C  XC      -  Cartesian coordinates
C  IAN     -  list of atomic numbers
C  FC      -  space for primitive force constants
C  NCmp    -  total number of non-zero natural internal components
C  NP1     -  number of primitives in each NIC
C  INT1    -  indices for all non-zero components per NIC
C  UT      -  non-zero components ordered per NIC
C  HINT    -  on exit contains scaled-guess Hessian/Inverse Hessian
C
C
      REAL*8 XC(3,*),FC(intcor),UT(NCmp),HINT(NDEG,NDEG)
      DIMENSION ktyp(intcor),klist(4,intcor),IAN(*)
      INTEGER NP1(NDEG),INT1(NCmp)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
C  form primitive force constant array
C
      DO 10 I=1,intcor
      FC(I) = ForceCnst(ktyp(I),klist(1,I),klist(2,I),klist(3,I),
     $                  klist(4,I),NCon,IAN,XC)
 10   CONTINUE
C
C  now form the Hessian
C  this is simply UT(t)*UT scaled by force constants
C  ** LOOKS LIKE THIS NEEDS TO BE SORTED OUT **
C
      DO 30 I=1,NDEG
      DO 30 J=1,I
      Val = Zero
      DO 20 K=1,NP1(I)
cc      Val = Val + UT(K,I)*UT(K,J)*FC(K)
 20   CONTINUE
      HINT(I,J) = Val
      HINT(J,I) = Val
 30   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE S_GrdINT2(NAT3,   N,      gcart,  IPRNT,  Nb,
     $                     NB1,    INB1,   B,      NB2,    INB2,
     $                     BT,     Bg,     D,      Z,      gint,
     $                     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Iterative transformation of Cartesian gradient to gradient
C  over delocalized internal coordinates
C  Solve
C     Gint(k+1)  =  D**-1 { B*Gcart + (D - B*Bt)*Gint(k) }
C  where D is the diagonal of B*Bt
C  ** SPARSE MATRIX ALGORITHM **
C
C  ARGUMENTS
C
C  NAT3    -  3*number of atoms
C  N       -  size of optimization space
C  gcart   -  gradient in Cartesian coordinates
C  IPRNT   -  print flag
C  Nb      -  total number of non-zero B-matrix elements
C  NB1     -  number of non-zero entries for each row
C  INB1    -  indices for all non-zero B-matrix rows
C  B       -  non-zero B-matrix elements ordered per row
C  NB2     -  number of non-zero entries for each column
C  INB2    -  indices for all non-zero B-matrix columns
C  BT      -  non-zero B-matrix elements ordered per column
C  Bg      -  storage for B*gcart
C  D       -  storage for diagonal elements of B*BT
C  Z       -  general scratch storage (at least 4*N + NAT3)
C  gint    -  on exit contains gradient in delocalized internals
C  IErr    -  error flag    0 - new gradient found
C                          -1 - unable to converge to new gradient
C
C
      DIMENSION gcart(NAT3),B(Nb),BT(Nb),Bg(N),D(N),
     $          Z(*),gint(N)
      INTEGER NB1(NAT3),INB1(Nb),NB2(N),INB2(Nb)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=5.0d-9,MaxIT=999)
C
C
      IOut = ioutfil('iout')
C
C  form BT*gcart
C
      np = 0
      DO 10 I=1,N
      Bg(I) = Zero
      DO 9 J=1,NB2(I)
      np = np+1
      Bg(I) = Bg(I) + BT(np)*gcart(INB2(np))
  9   CONTINUE
 10   CONTINUE
C
C  form and invert diagonal elements of B*B(t)
C
      np = 0
      DO 30 I=1,N
      Val = Zero
      DO 20 K=1,NB2(I)
      np = np+1
      Val = Val + BT(np)**2
 20   CONTINUE
      D(I) = One/Val
      gint(I) = Bg(I)*D(I)       ! first estimate of gradient
 30   CONTINUE
c
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
C  now call <CONJUGATE> to solve linear equations
C
      i1 = 1
      i2 = i1 + N
      i3 = i2 + N
      i4 = i3 + N
      i5 = i4 + N
c
      CALL S_CONJUGATE(NAT3,   N,      Nb,     NB1,    INB1,
     $                 B,      NB2,    INB2,   BT,     Bg,
     $                 D,      gint,   thrsh,  MaxIT,  IPRNT,
     $                 Z(i1),  Z(i2),  Z(i3),  Z(i4),  Z(i5),
     $                 IErr)
C
      RETURN
c
 1000 FORMAT(/,5X,'Iterative generation of Internal Gradient')
c
      END
c =======================================================================
c
      SUBROUTINE S_GrdINT3(NAT3,   N,      gcart,  IPRNT,  Nb,
     $                     NB1,    INB1,   B,      NB2,    INB2,
     $                     BT,     Y,      D,      Z,      gint,
     $                     IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Iterative transformation of Cartesian gradient to gradient
C  over delocalized internal coordinates
C  Solve
C     Y(k+1)  =  D**-1 { Gcart + (D - Bt*B)*Y(k) }
C  where D is the diagonal of Bt*B
C  Then
C     Gint = B*Y
C  ** SPARSE MATRIX ALGORITHM **
C
C  ARGUMENTS
C
C  NAT3    -  3*number of atoms
C  N       -  size of optimization space
C  gcart   -  gradient in Cartesian coordinates
C  IPRNT   -  print flag
C  Nb      -  total number of non-zero B-matrix elements
C  NB1     -  number of non-zero entries for each row
C  INB1    -  indices for all non-zero B-matrix rows
C  B       -  non-zero B-matrix elements ordered per row
C  NB2     -  number of non-zero entries for each column
C  INB2    -  indices for all non-zero B-matrix columns
C  BT      -  non-zero B-matrix elements ordered per column
C  Y       -  storage for Y-vector
C  D       -  storage for diagonal elements of BT*B
C  Z       -  general scratch storage (at least 4*NAT3 + N)
C  gint    -  on exit contains gradient in delocalized internals
C  IErr    -  error flag    0 - new gradient found
C                          -1 - unable to converge to new gradient
C
C
      DIMENSION gcart(NAT3),B(Nb),BT(Nb),Y(NAT3),D(NAT3),
     $          Z(*),gint(N)
      INTEGER NB1(NAT3),INB1(Nb),NB2(N),INB2(Nb)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=5.0d-9,MaxIT=999)
C
C
      IOut = ioutfil('iout')
C
C  form and invert diagonal elements of B(t)*B
C
      np = 0
      DO 30 I=1,NAT3
      Val = Zero
      DO 20 K=1,NB1(I)
      np = np+1
      Val = Val + B(np)**2
 20   CONTINUE
      D(I) = One/Val
      Y(I) = gcart(I)*D(I)       ! first estimate of Y
 30   CONTINUE
c
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
C  now call <CONJUGATE> to solve linear equations
C
      i1 = 1
      i2 = i1 + NAT3
      i3 = i2 + NAT3
      i4 = i3 + NAT3
      i5 = i4 + NAT3
c
      CALL S_CONJUGATE(N,      NAT3,   Nb,     NB2,    INB2,
     $                 BT,     NB1,    INB1,   B,      gcart,
     $                 D,      Y,      thrsh,  MaxIT,  IPRNT,
     $                 Z(i1),  Z(i2),  Z(i3),  Z(i4),  Z(i5),
     $                 IErr)

C
C  form internal gradient from Y
C
      np = 0
      DO 40 I=1,N
      gint(I) = Zero
      DO 39 J=1,NB2(I)
      np = np+1
      gint(I) = gint(I) + BT(np)*Y(INB2(np))
 39   CONTINUE
 40   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,5X,'Iterative generation of Internal Gradient')
c
      END
c =======================================================================
c
      SUBROUTINE S_TranINT(NPrim,  NDEG,   XPrim,  NCmp,   NP1,
     $                     INT1,   UT1,    XINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transforms primitive internal coordinates to compound
C  natural internal set
C  ** SPARSE MATRIX VERSION **
C
      REAL*8 XPrim(NPrim),UT1(NCmp),XINT(NDEG)
      INTEGER NP1(NDEG),INT1(NCmp)
C
      PARAMETER (Zero=0.0d0)
C
      np = 0
      DO 20 I=1,NDEG
      XINT(I) = Zero
      DO 10 J=1,NP1(I)
      np = np+1
      XINT(I) = XINT(I) + UT1(np)*XPrim(INT1(np))
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SCHMIDT(NDEG,   intcor, NCons,  VC,     thrsh,
     $                   IPRNT,  VM,     NS,     LCON,   UT,
     $                   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Modify the non-redundant set of natural internals
C  to account for any imposed constraints
C
C  What is done is to Schmidt-orthogonalize the current set of active
C  internal coordinates to the set of constraint vectors.
C
C  ARGUMENTS
C
C  NDEG    -  number of active internals
C             (number of degrees of freedom)
C  intcor  -  number of primitive internals
C  NCons   -  number of constrained (fixed) primitives
C             on exit actual number of constraints
C  VC      -  initial set of constraint vectors
C             (one vector for each fixed primitive)
C  thrsh   -  zero norm threshold for rejecting internal coordinate
C  IPRNT   -  print flag
C  VM      -  scratch array for vectors
C  NS      -  number of active internals left after
C             Schmidt-Orthogonalization
C  LCON    -  on exit indicates which constraints are "active" (independent)
C              LCON(I) = 0  constraint is active
C              LCON(I) = 1  constraint should be eliminated
C  UT      -  on entry contains current set of natural internals
C             on exit contains Schmidt-Orthogonalized set
C  IErr    -  error flag
C
C
      REAL*8 VM(intcor,NCons+NDEG),UT(intcor,NCons+NDEG),
     $       VC(intcor,NCons)
      DIMENSION LCON(NCons)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      IOut = ioutfil('iout')
c
      IErr = -1
      CALL IZeroIT(LCON,NCons)
c
      If(IPRNT.GT.2) WRITE(IOut,1000)
C
C  normalize all the constraint vectors
C  (at the same time eliminate zero constraint vectors)
C
      II = 0
      DO 40 I=1,NCons
      II = II+1
      snorm = SProd(intcor,VM(1,I),VM(1,I))
      IF(snorm.LT.thrsh) THEN
       If(IPRNT.GT.2) WRITE(IOut,1100) I
       LCON(I) = 1
       II = II-1
      ELSE
       snorm = One/SQRT(snorm)
       CALL VScal(intcor,snorm,VM(1,I))
       If(II.NE.I) Then
        CALL CpyVEC(intcor,VM(1,I),VM(1,II))
        CALL CpyVEC(intcor,VC(1,I),VC(1,II))
       EndIf
      ENDIF
 40   CONTINUE
      NCon = II
c
      If(IPRNT.GT.5) Then
       WRITE(IOut,1200)
       CALL PrntMAT(NCon,intcor,NCon,VM)
      EndIf
C
C  Constraints may not be independent
C  see if any can be eliminated by Schmidt-Orthogonalization
C
      IT = 1
      II = 2
      DO 70 I=2,NCons
      If(LCON(I).EQ.1) GO TO 70
      IT = IT+1
      DO 60 J=II-1,1,-1
c
      cf = SProd(intcor,VM(1,J),VM(1,IT))
      DO 50 K=1,intcor
      VM(K,IT) = VM(K,IT) - cf*VM(K,J)
 50   CONTINUE
c
      snorm = SProd(intcor,VM(1,IT),VM(1,IT))
c
      If(snorm.LT.thrsh) Then
       LCON(I) = 1
       GO TO 70
      EndIf
c
 60   CONTINUE
C
C  normalize current constraint vector
C
      snorm = One/SQRT(snorm)
      CALL VScal(intcor,snorm,VM(1,IT))
      If(II.NE.IT) Then
       CALL CpyVEC(intcor,VM(1,IT),VM(1,II))
       CALL CpyVEC(intcor,VC(1,IT),VC(1,II))
      EndIf
      II = II+1
c
 70   CONTINUE
C
C  number of constraints surviving is II-1
C
      NCons = NCon
      IF(II-1.NE.NCons) THEN
       If(IPRNT.GT.2) WRITE(IOut,1300) NCons-II+1
       NCons = II-1
      ENDIF
c
      If(IPRNT.GT.5) Then
       WRITE(IOut,1400)
       CALL PrntMAT(NCons,intcor,NCons,VM)
      EndIf
C
C  Now Schmidt-Orthogonalize the active internals
C  to the remaining constraint vectors
C
      II = NCons+1
      DO 100 I=1,NDEG
      CALL CpyVEC(intcor,UT(1,I),VM(1,II))
c
      DO 90 J=II-1,1,-1
c
      cf = SProd(intcor,VM(1,J),VM(1,II))
      DO 80 K=1,intcor
      VM(K,II) = VM(K,II) - cf*VM(K,J)
 80   CONTINUE
c
      snorm = SProd(intcor,VM(1,II),VM(1,II))
c
      IF(snorm.LT.thrsh) THEN
       If(IPRNT.GT.2) WRITE(IOut,1500) I
cc       II = II-1     ! WRONG  (JB)
       GO TO 100
      ENDIF
c
 90   CONTINUE
C
C  normalize current internal coordinate
C
      snorm = One/SQRT(snorm)
      CALL VScal(intcor,snorm,VM(1,II))
      II = II+1
      If(II.GT.NDEG) Then
       II = II-1
       Exit
      EndIf
c
 100  CONTINUE
C
C  copy orthogonalized active internals back into UT
C
      NS = II-NCons
c
      If(IPRNT.GT.2) WRITE(IOut,1600) NCons,NS
c
      DO 110 I=1,NS
      CALL CpyVEC(intcor,VM(1,I+NCons),UT(1,I))
 110  CONTINUE
C
C  now restore constraint vectors
C  (these are needed in order to span the full space
C   for the back-transformation)
C
      DO 120 I=1,NCons
      II = I+NS
      CALL CpyVEC(intcor,VC(1,I),UT(1,II))
 120  CONTINUE
C
C  possibly see what we've got
C
      IF(IPRNT.GT.4) THEN
       WRITE(IOut,1700) NS,NCons
       CALL PrntMAT(NS+NCons,intcor,NS+NCons,UT)
      EndIF
C
C  set all small components to zero
C
      DO 150 I=1,NDEG
      DO 150 J=1,intcor
      If(Abs(UT(J,I)).LT.thrsh) UT(J,I) = Zero
 150  CONTINUE
C
C  check for error
C
      IF(NS+NCons.NE.NDEG) THEN
       WRITE(IOut,1800)
      ELSE
       IErr = 0
      ENDIF
cc      write(6,2222)
cc 2222 format(/,' $$$$$ Turning Off ALL Constraints $$$$',/)
cc      ns = ndeg
cc      ncon = 0
C
      RETURN
c
 1000 FORMAT(/,' Imposing Constraints by Schmidt Orthogonalization')
 1100 FORMAT(' Eliminating Constraint ',I5,'  This has a Zero',
     $       ' Constraint Vector')
 1200 FORMAT(' Normalized Projected Constraint Vectors')
 1300 FORMAT(' Imposed Constraints are not Independent.  Eliminating ',
     $         I5,' Constraints')
 1400 FORMAT(' Final Set of Constraint Vectors')
 1500 FORMAT('**WARNING** Internal Coordinate ',I5,' has Norm Below',
     $        ' Threshold - Eliminating')
 1600 FORMAT(' Eliminated ',I5,' Coordinates',/,' Number of Internal',
     $       ' Coordinates Left in Active Space is ',I5)
 1700 FORMAT(/,'   Schmidt-Orthogonalized Set of ',I5,' Active and ',I5,
     $         ' Constraint Vectors')
 1800 FORMAT(/,2X,'***ERROR*** Optimization Space has Incorrect',
     $            ' Dimension',/,
     $            '   after Schmidt-Orthogonalization of Constraints')
c
      END
c =======================================================================
c
      SUBROUTINE SortBEND(int0,klist,nc,mc)
      IMPLICIT INTEGER(A-Z)
C
C  Checks whether primitive in current natural internal coordinate
C  has not already been found in a previous coordinate
C
C  ARGUMENTS
C
C  int0    -  start of bends in klist array
C  klist   -  list of atoms involved in each primitive bend
C  nc      -  current location of bend in klist array (end)
C             if bend exists, will be decremented by 1 on exit
C  mc      -  final location of bend in klist array
C
C
      DIMENSION klist(4,nc)
C
C
C  check from start of bends to end whether current bend exists
C
      DO 10 I=int0,nc-1
      IF(klist(3,I).EQ.klist(3,nc)) THEN
C
C  central atom of bend is identical
C  check ends
C
       If( (klist(1,I).EQ.klist(1,nc).AND.
     $      klist(2,I).EQ.klist(2,nc)) .OR.
     $     (klist(1,I).EQ.klist(2,nc).AND.
     $      klist(2,I).EQ.klist(1,nc)) ) Then
C
C  bend is identical
C
         mc = I
         nc = nc-1
         RETURN
       EndIf
      ENDIF
 10   CONTINUE
C
C  If we gat here, didn't find bend
C
      mc = nc
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SortINT(NPrim0, NPrim1, ktyp1,  klist1, NPrim2,
     $                   ktyp2,  klist2, Change, MPrim,  INDX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Compares the underlying set of primitive internal coordinates as
C  generated on two different optimization cycles
C  ** CURRENTLY HANDLES INVERSE DISTANCE COORDINATES ONLY **
C
C  ARGUMENTS
C
C  NPrim0  -  number of intramolecular primitives
C             (assumed to remain the same)
C  NPrim1  -  number of primitive internals in first coordinate set
C  ktyp1   -  primitive internal type for first set
C  klist1  -  list of atoms involved in each primitive
C  NPrim2  -  number of primitive internals in second coordinate set
C  ktyp2   -  primitive internal type for second set
C  klist2  -  list of atoms involved in each primitive
C
C  on exit
C
C  Change  -  logical flag indicating whether primitive space has changed
C  MPrim   -  total number of entries in INDX array
C  INDX    -  on exit contains array of differences between two sets
C             (specifically what to do to first set to get second)
C                0 - primitive common to both sets
C               -1 - primitive needs to be deleted from first set
C               +1 - primitive needs to be added to first set
C
C
      DIMENSION ktyp1(NPrim1),klist1(4,NPrim1),INDX(*),
     $          ktyp2(NPrim2),klist2(4,NPrim2)
      LOGICAL Change
C
C
      IOut = ioutfil('iout')
      Change = .FALSE.
C
C  loop over two different sets of inverse-power distances
C  ** IMPORTANT **
C  It is assumed that all other primitives remain the SAME
C
      I1 = NPrim0
      I2 = NPrim0
      IT = NPrim0
C
C  set first NPrim0 entries of INDX array to zero
C
      DO 5 I=1,NPrim0
      INDX(I) = 0
  5   CONTINUE
C
C  now check the intermolecular inverse-power distance primitives
C
 10   CONTINUE
      I1 = I1+1
      I2 = I2+1
      IT = IT+1
cc
      IF(klist1(1,I1).EQ.klist2(1,I2)) THEN
        If(klist1(2,I1).LT.klist2(2,I2)) Then
         INDX(IT) = -1
         I2 = I2-1
        Else If(klist1(2,I1).GT.klist2(2,I2)) Then
         INDX(IT) = +1
         I1 = I1-1
        Else
         INDX(IT) = 0
        EndIf
      ELSE
        If(klist1(1,I1).LT.klist2(1,I2)) Then
         INDX(IT) = -1
         I2 = I2-1
        Else
         INDX(IT) = +1
         I1 = I1-1
        EndIf
      ENDIF
cc
      If(I1.LT.NPrim1.AND.I2.LT.NPrim2) GO TO 10
C
C  check for last entry
C
      If(I2.EQ.NPrim2) Then
       NLeft = NPrim1 - I1
      Else
       NLeft = I2 - NPrim2
      EndIf
c
      MPrim = IT + Abs(NLeft)
c
      IF(NLeft.GT.0) THEN
       DO 20 I=IT+1,IT+NLeft
       INDX(I) = -1
 20    CONTINUE
      ELSE IF(NLeft.LT.0) THEN
       DO 30 I=IT+1,IT-NLeft
       INDX(I) = 1
 30    CONTINUE
      ENDIF
C
C  Check dimension is correct
C  NPrim1 + changes should equal NPrim2
C
      ITot = 0
      DO 40 I=1,MPrim
      If(INDX(I).NE.0) Change = .TRUE.
      ITot = ITot + INDX(I)
 40   CONTINUE
c
      If( (NPrim1+ITot).NE.NPrim2 ) Then
       WRITE(IOut,1000)
       CALL OptExit(9)
      EndIf
C
      If(Change) WRITE(IOut,1100)
      If(.NOT.Change) WRITE(IOut,1200)
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Problem with Dimension of Primitive',
     $              ' Space')
 1100 FORMAT(' Change in Underlying Primitive Space')
 1200 FORMAT(' Underlying Primitive Space Unchanged')
c
      END
c =======================================================================
c
      SUBROUTINE SortOUTP(int0,klist,nc,mc)
      IMPLICIT INTEGER(A-Z)
C
C  Checks whether out-of-plane bend in current natural internal coordinate
C  has not already been found in a previous coordinate
C
C  ARGUMENTS
C
C  int0    -  start of out-of-plane bends in klist array
C  klist   -  list of atoms involved in each primitive
C  nc      -  current location of primitive in klist array (end)
C             if out-of-plane bend exists, will be decremented by 1 on exit
C  mc      -  final location of out-of-plane bend in klist array
C
C
      DIMENSION klist(4,nc)
C
C
C  check from start of list to end whether current out-of-plane bend exists
C
      DO 10 I=int0,nc-1
      IF( (klist(1,I).EQ.klist(1,nc) .AND.
     $     klist(2,I).EQ.klist(2,nc) .AND.
     $     klist(3,I).EQ.klist(3,nc) .AND.
     $     klist(4,I).EQ.klist(4,nc))  .OR.
     $    (klist(1,I).EQ.klist(1,nc) .AND.
     $     klist(2,I).EQ.klist(3,nc) .AND.
     $     klist(3,I).EQ.klist(2,nc) .AND.
     $     klist(4,I).EQ.klist(4,nc)) ) Then
C
C  out-of-plane bend is identical
C
        mc = I
        nc = nc-1
        RETURN
      ENDIF
 10   CONTINUE
C
C  If we get here, didn't find this out-of-plane bend
C
      mc = nc
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SortPRIM(NPrim,  ktyp,   klist)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads and sorts out primitives generated from routine <INTC>
C
C  ARGUMENTS
C
C  NPrim   -  on entry maximum size of primitive space
C             on exit number of primitives found
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C  klist   -  list of atoms involved in each primitive internal
C
C
      DIMENSION ktyp(NPrim),klist(4,NPrim)
c
      Character*4 Type(6),typ
      Data Type / 'STRE','BEND','OUTP','TORS','LIN1','LIN2' /
C
C
C  initialize
C
      CALL IZeroIT(klist,4*NPrim)
      IType = 0
      mc = 0          ! number of primitives
      int0 = 1        ! positioned at first stretch
C
      OPEN (UNIT=7,FILE='NIC',FORM='FORMATTED',STATUS='OLD')
C
C  Loop over all 6 allowed coordinate types
C
      DO 40 IType=1,6
C
C  Read information from unit 7 until hit EOF
C
 30   CONTINUE
      READ(7,910,END=96) typ,IC,i1,i2,i3,i4
  910 Format(20X,A4,I4,2X,4(I4,5X))
c
      IF(typ.EQ.Type(IType)) THEN
c
       IF(typ.EQ.'STRE') THEN
C
C  stretches
C
        mc = mc+1
        ktyp(mc) = 1
        klist(1,mc) = i1
        klist(2,mc) = i2
C
C  check that we don't already have this primitive
C
        CALL SortSTRE(int0,klist,mc,lc)
c
       ELSE IF(typ.EQ.'BEND') THEN
C
C  bends
C
        mc = mc+1
        ktyp(mc) = 2
        klist(1,mc) = i1
        klist(2,mc) = i2
        klist(3,mc) = i3
C
C  check that we don't already have this primitive
C
        CALL SortBEND(int0,klist,mc,lc)
c
       ELSE
C
C  all other coordinate types
C
        mc = mc+1
        ktyp(mc) = IType
        klist(1,mc) = i1
        klist(2,mc) = i2
        klist(3,mc) = i3
        klist(4,mc) = i4
C
C  check that we don't already have this primitive
C
        If(typ.EQ.'TORS') Then
         CALL SortTORS(int0,klist,mc,lc)
        Else If(typ.EQ.'OUT ') Then
         CALL SortOUTP(int0,klist,mc,lc)
        Else
         lc = mc    ! WARNING - assume only a few linear coordinates
        EndIf
c
       ENDIF
      ENDIF
c
      GO TO 30
 96   CONTINUE
C
C  reached end of file
C
      If(IType.EQ.6) GO TO 50
      int0 = mc+1        ! marker to beginning of next primitive type
c
      REWIND 7
 40   CONTINUE
 50   CONTINUE
C
      CLOSE (UNIT=7,STATUS='KEEP')
C
C  finished sorting coordinates
C  number of primitives found is mc
C
      NPrim = mc
c
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SortSTRE(int0,klist,nc,mc)
      IMPLICIT INTEGER(A-Z)
C
C  Checks whether stretch in current natural internal coordinate
C  has not already been found in a delocalized ring coordinate
C
C  ARGUMENTS
C
C  int0    -  start of stretches in klist array
C  klist   -  list of atoms involved in each primitive stretch
C  nc      -  current location of stretch in klist array (end)
C             if stretch exists, will be decremented by 1 on exit
C  mc      -  final location of stretch in klist array
C
C
      DIMENSION klist(4,nc)
C
C
C  check from start of stretch to end whether current stretch exists
C
      DO 10 I=int0,nc-1
      If( (klist(1,I).EQ.klist(1,nc).AND.
     $     klist(2,I).EQ.klist(2,nc)) .OR.
     $    (klist(1,I).EQ.klist(2,nc).AND.
     $     klist(2,I).EQ.klist(1,nc)) ) Then
C
C  stretch is identical
C
        mc = I
        nc = nc-1
        RETURN
      EndIf
 10   CONTINUE
C
C  If we gat here, didn't find stretch
C
      mc = nc
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SortTORS(int0,klist,nc,mc)
      IMPLICIT INTEGER(A-Z)
C
C  Checks whether torsion in current natural internal coordinate
C  has not already been found in a previous coordinate
C
C  ARGUMENTS
C
C  int0    -  start of torsions in klist array
C  klist   -  list of atoms involved in each primitive torsion
C  nc      -  current location of torsion in klist array (end)
C             if torsion exists, will be decremented by 1 on exit
C  mc      -  final location of torsion in klist array
C
C
      DIMENSION klist(4,nc)
C
C
C  check from start of torsions to end whether current torsion exists
C
      DO 10 I=int0,nc-1
      IF( (klist(1,I).EQ.klist(1,nc) .AND.
     $     klist(2,I).EQ.klist(2,nc) .AND.
     $     klist(3,I).EQ.klist(3,nc) .AND.
     $     klist(4,I).EQ.klist(4,nc))  .OR.
     $    (klist(1,I).EQ.klist(4,nc) .AND.
     $     klist(2,I).EQ.klist(3,nc) .AND.
     $     klist(3,I).EQ.klist(2,nc) .AND.
     $     klist(4,I).EQ.klist(1,nc)) ) Then
C
C  torsion is identical
C
        mc = I
        nc = nc-1
        RETURN
      ENDIF
 10   CONTINUE
C
C  If we gat here, didn't find torsion
C
      mc = nc
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SurfCNNCT(NAtoms, NMol,   IMOL,   XC,     NIC,
     $                     Cutoff, IPRNT,  IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates additional connectivity between adsorbed molecule
C  and the surface
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  NMol    -  number of molecules (should be two in this implementation)
C  IMOL    -  pointers to start/end of molecules in XC array
C  XC      -  Cartesian coordinates
C  NIC     -  maximum number of primitive internal coordinates
C             that can be generated
C             ** NOTE:  This is typically far more than 3*NAtoms-6
C                       as many redundancies are found that are
C                       dealt with later
C  CutOff  -  distance cutoff for bonding between adsorbate and surface
C              (in Angstroms**2)
C             atoms less than a distance cutoff from one another
C             are considered as bonded
C  IPRNT   -  print flag
C  IC      -  connectivity matrix
C             (may be modified on exit)
C
      DIMENSION XC(3,NAtoms),IMOL(NMol+1)
      DIMENSION IC(NAtoms,NAtoms)
C
      PARAMETER (ANTOAU=1.889725991249198d0)
      PARAMETER (One=1.0d0,CutMax=10.0d0)
C
C
      IOut = ioutfil('iout')
c
      CutOff0 = CutOff                 ! save initial cutoff
      nc = 0                           ! number of connections
c
 500  CONTINUE
c
      If(IPRNT.GT.1) Then
       WRITE(IOut,1000)
       WRITE(IOut,1100) CutOff
      EndIf
C
C  convert cutoff to distance squared in au
C
      CutOf2 = (CutOff*ANTOAU)**2
C
C  generate additional connectivity between adsorbate and surface
C
c -- surface
      IStrt = IMOL(1)+1
      IEnd  = IMOL(2)
c -- adsorbate
      JStrt = IMOL(2)+1
      JEnd  = IMOL(3)
C
C  Loop over atoms in the two regions
C
      DO 20 I=IStrt,IEnd
      DO 10 J=JStrt,JEnd
C
C  get interatomic distance I-J
C
      X = XC(1,I) - XC(1,J)
      Y = XC(2,I) - XC(2,J)
      Z = XC(3,I) - XC(3,J)
      R2 = X*X + Y*Y + Z*Z
C
C  atoms "bonded" if distance less than cutoff
C
      IF(R2.LT.CutOf2) THEN
       IC(I,J) = 1
       IC(J,I) = 1
       nc = nc+1
      ENDIF
 10   CONTINUE
 20   CONTINUE
C
C  check
C  do we have any connectivity?
C
      IF(nc.EQ.0) THEN
       CutOff = CutOff + One
       If(IPRNT.GT.1) WRITE(IOut,1500)
       If(CutOff.GT.CutMax) Then
        WRITE(IOut,2000) CutMax
        CALL OptExit(9)
       EndIf
       GO TO 500
      ENDIF
c
 1000 FORMAT(' Generating Connectivity Between Adsorbate and Surface')
 1100 FORMAT(' Cutoff for bonding is ',F10.6,' Angstroms')
 1500 FORMAT('**WARNING**  No Connectivity - Increasing Cutoff')
 2000 FORMAT(/,2X,'***ERROR*** Adsorbate More Than ',F8.4,' Angstroms ',
     $            ' From Surface - NO BONDING!')
c
      END
c =======================================================================
c
      SUBROUTINE SurfINT(NAtoms, NMol,   IMOL,   XC,     PThrsh,
     $                   NIC,    Cutoff, IPRNT,  IC,     intcor,
     $                   ktyp,   klist)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates linking primitive internals between adsorbed
C  molecule and the surface
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atoms
C  NMol    -  number of molecules (should be two in this implementation)
C  IMOL    -  pointers to start/end of molecules in XC array
C  XC      -  Cartesian coordinates
C  PThrsh  -  smallest value (in radians) allowed for near-linear angle
C  NIC     -  maximum number of primitive internal coordinates
C             that can be generated
C             ** NOTE:  This is typically far more than 3*NAtoms-6
C                       as many redundancies are found that are
C                       dealt with later
C  CutOff  -  distance cutoff for bonding between adsorbate and surface
C              (in Angstroms**2)
C             atoms less than a distance cutoff from one another
C             are considered as bonded
C  IPRNT   -  print flag
C  IC      -  connectivity matrix
C
C  on exit
C
C  intcor  -  number of primitive internals found
C  ktyp    -  integer array containing internal coordinate type
C             ** NOTE:  These types are:
C                       1 - stretch
C                       2 - bend
C                       3 - out-of-plane bend
C                       4 - torsion
C                       5 - linear coplanar bend
C                       6 - linear perpendicular bend
C                       7 - inverse power stretch
C  klist   -  list of atoms involved in each primitive
C
C
      DIMENSION XC(3,NAtoms),IMOL(NMol+1)
      DIMENSION IC(NAtoms,NAtoms),ktyp(NIC),klist(4,NIC)
C
      PARAMETER (ANTOAU=1.889725991249198d0)
      PARAMETER (One=1.0d0,CutMax=10.0d0)
C
C
      IOut = ioutfil('iout')
c
      CutOff0 = CutOff                 ! save initial cutoff
      intcor0 = intcor                 !  ditto
      nc = 0                           ! number of connections
c
 500  CONTINUE
c
      If(IPRNT.GT.1) Then
       WRITE(IOut,1000)
       WRITE(IOut,1100) CutOff
      EndIf
C
C  convert cutoff to distance squared in au
C
      CutOf2 = (CutOff*ANTOAU)**2
C
C  generate linking primitives between adsorbate and surface
C
c -- surface
      IStrt = IMOL(1)+1
      IEnd  = IMOL(2)
c -- adsorbate
      JStrt = IMOL(2)+1
      JEnd  = IMOL(3)
C
C  Loop over atoms in the two regions
C
      DO 20 I=IStrt,IEnd
      DO 10 J=JStrt,JEnd
C
C  get interatomic distance I-J
C
      X = XC(1,I) - XC(1,J)
      Y = XC(2,I) - XC(2,J)
      Z = XC(3,I) - XC(3,J)
      R2 = X*X + Y*Y + Z*Z
C
C  atoms "bonded" if distance less than cutoff
C
      IF(R2.LT.CutOf2) THEN
       IC(I,J) = 1
       IC(J,I) = 1
       nc = nc+1
      ENDIF
 10   CONTINUE
 20   CONTINUE
C
C  check
C  do we have any connectivity?
C
      IF(nc.EQ.0) THEN
       CutOff = CutOff + One
       If(IPRNT.GT.1) WRITE(IOut,1500)
       If(CutOff.GT.CutMax) Then
        WRITE(6,2000) CutMax
        CALL OptExit(9)
       EndIf
       GO TO 500
      ENDIF
C
C  assign connecting primitives
C
      DO 60 I=IStrt,IEnd
      DO 50 J=JStrt,JEnd
      IF(IC(I,J).EQ.1) THEN
C
C  atoms I and J are connected
C  add the stretch I-J
C
       intcor = intcor+1
       ktyp(intcor) = 1
       klist(1,intcor) = I
       klist(2,intcor) = J
C
C  add bends for whatever atoms in the molecule J is connected to
C
       DO 40 K=JStrt,JEnd
       IF(IC(J,K).EQ.1) THEN
C
C  do not add if bend is linear
C
        CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
        If(Th.GT.PThrsh) GO TO 40
c
        intcor = intcor+1
        ktyp(intcor) = 2
        klist(1,intcor) = I
        klist(2,intcor) = K
        klist(3,intcor) = J
C
C  now add torsions involving atoms J,K and the surface
C
        DO 30 L=IStrt,IEnd
        IF(IC(I,L).EQ.1) THEN
C
C  do not add if L-I-J is linear
C
         CALL AngGRAD(NAtoms,L,I,J,XC,Th,.false.,jnk)
         If(Th.GT.PThrsh) GO TO 30
c
         intcor = intcor+1
         ktyp(intcor) = 4
         klist(1,intcor) = L
         klist(2,intcor) = I
         klist(3,intcor) = J
         klist(4,intcor) = K
        ENDIF
 30     CONTINUE
       ENDIF
 40    CONTINUE
C
      ENDIF
 50   CONTINUE
 60   CONTINUE
C
C  now look for additional bends involving 2 surface atoms
C  ** Beware of Duplication **
C
      intcor1 = intcor+1
c
      DO 90 J=JStrt,JEnd
      DO 80 I=IStrt,IEnd
      IF(IC(I,J).EQ.1) THEN
C
C  look for surface bend
C
       DO 70 K=IStrt,IEnd
       IF(IC(I,K).EQ.1) THEN
C
C  do not add if bend is linear
C
        CALL AngGRAD(NAtoms,J,I,K,XC,Th,.false.,jnk)
        If(Th.GT.PThrsh) GO TO 70
c
        intcor = intcor+1
        ktyp(intcor) = 2
        klist(1,intcor) = J
        klist(2,intcor) = K
        klist(3,intcor) = I
C
C  check for duplication
C
        DO L=intcor1,intcor-1
        IF(klist(3,L).EQ.I.AND.
     $    (klist(1,L).EQ.J.AND.klist(2,L).EQ.K) .OR.
     $    (klist(1,L).EQ.K.AND.klist(2,L).EQ.J) ) THEN
         intcor = intcor-1
         GO TO 70
        ENDIF
        ENDDO
       ENDIF
 70    CONTINUE
C
      ENDIF
 80   CONTINUE
 90   CONTINUE
c
      If(IPRNT.GT.4) Then
       WRITE(IOut,1200)
       Do I=intcor0+1,intcor
       If(klist(3,I).EQ.0) Then
        WRITE(IOut,1300) klist(1,I),klist(2,I)
       Else If(klist(4,I).EQ.0) Then
        WRITE(IOut,1310) klist(1,I),klist(3,I),klist(2,I)
       Else
        WRITE(IOut,1320) klist(1,I),klist(2,I),klist(3,I),klist(4,I)
       EndIf
       EndDo
      EndIf
c
      If(IPRNT.GT.2) Then
       WRITE(IOut,1400) intcor-intcor0
       WRITE(IOut,1450) intcor
      EndIf
c
      CutOff = CutOff0             ! restore original cutoff
      RETURN
c
 1000 FORMAT(' Generating Primitives Connecting Adsorbate and Surface')
 1100 FORMAT(' Cutoff for bonding is ',F10.6,' Angstroms')
 1200 FORMAT(' Additional Linking Primitives are:')
 1300 FORMAT(1X,' Stretch:',5X,2I5)
 1310 FORMAT(1X,' Bend:   ',5X,3I5)
 1320 FORMAT(1X,' Torsion:',5X,4I5)
 1400 FORMAT(' Added ',I4,' Linking Primitives')
 1450 FORMAT(' There are now ',I5,' Primitives')
 1500 FORMAT('**WARNING**  No Connectivity - Increasing Cutoff')
 2000 FORMAT(/,2X,'***ERROR*** Adsorbate More Than ',F8.4,' Angstroms ',
     $            ' From Surface - NO BONDING!')
c
      END
c =======================================================================
c
      SUBROUTINE SymCON(NAtoms, NCons,  ICTYP,  RCON,   IC,
     $                  NTrans, NEqATM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks whether the imposed constraints will preserve the
C  molecular symmetry.  If not, prints out the additional
C  constraints required and Exits
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NCons   -  number of constraints
C  ICTYP   -  constraint type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane-bend
C              4 - fixed dihedral angle
C              5 - fixed linear coplanar bend
C              6 - fixed linear perpendicular bend
C              9 - composite constraint
C  RCON    -  value of constraint
C  IC      -  atoms involved in constraint
C               IC1-IC2           distance constraint
C               IC1-IC2-IC3       bond angle constraint
C               IC1-IC2-IC3-IC4   dihedral constraint
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C
      DIMENSION ICTYP(NCons),RCON(NCons),IC(4,NCons)
      DIMENSION NEqATM(NAtoms,NTrans)
C
C
      If(NTrans.EQ.1) RETURN     ! there is no symmetry to conserve
c
      IOut = ioutfil('iout')
      IErr = 0
C
      DO 80 IQ=1,NCons
C
C  Jump out for composite constraints
C
      If(ICTYP(IQ).EQ.9) GO TO 80
c
      IF(ICTYP(IQ).EQ.1) THEN
cc
C  distance constraint (I-J)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       VAL = RCON(IQ)
C
C  check for equivalent atoms
C  these must also be constrained for symmetry
C  to be maintained
C
       DO 10 IOP=2,NTrans
       IF(I.GT.NAtoms) THEN
        IT = I
       ELSE
        IT = NEqATM(I,IOP)
       ENDIF
       IF(J.GT.NAtoms) THEN
        JT = J
       ELSE
        JT = NEqATM(J,IOP)
       ENDIF
       If( (IT.EQ.I.AND.JT.EQ.J).OR.
     $     (JT.EQ.I.AND.IT.EQ.J) ) GO TO 10
C
C  check to see if this constraint is present
C
       DO 9 JQ=1,NCons
       IF(ICTYP(JQ).EQ.1) THEN
        IS = IC(1,JQ)
        JS = IC(2,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS).OR.(JT.EQ.IS.AND.IT.EQ.JS))
     $          .AND.RCON(JQ).EQ.VAL ) GO TO 10
       ENDIF
 9     CONTINUE
C
C  if we reach here the symmetry-required constraint is absent
C
       IErr = -1
       WRITE(IOut,1000) IT,JT,I,J
C
 10    CONTINUE
cc
      ELSE IF(ICTYP(IQ).EQ.2) THEN
cc
C  angle constraint (I-J-K)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       VAL = RCON(IQ)
C
C  check for equivalent atoms
C  these must also be constrained for symmetry
C  to be maintained
C
       DO 20 IOP=2,NTrans
       IF(I.GT.NAtoms) THEN
        IT = I
       ELSE
        IT = NEqATM(I,IOP)
       ENDIF
       IF(J.GT.NAtoms) THEN
        JT = J
       ELSE
        JT = NEqATM(J,IOP)
       ENDIF
       IF(K.GT.NAtoms) THEN
        KT = K
       ELSE
        KT = NEqATM(K,IOP)
       ENDIF
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K).OR.
     $     (KT.EQ.I.AND.JT.EQ.J.AND.IT.EQ.K) ) GO TO 20
C
C  check to see if this constraint is present
C
       DO 19 JQ=1,NCons
       IF(ICTYP(JQ).EQ.2) THEN
        IS = IC(1,JQ)
        JS = IC(2,JQ)
        KS = IC(3,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS.AND.KT.EQ.KS).OR.
     $       (KT.EQ.IS.AND.JT.EQ.JS.AND.IT.EQ.KS))
     $          .AND.RCON(JQ).EQ.VAL ) GO TO 20
       ENDIF
 19    CONTINUE
C
C  if we reach here the symmetry-required constraint is absent
C
       IErr = -1
       WRITE(IOut,1100) IT,JT,KT,I,J,K
C
 20    CONTINUE
cc
      ELSE IF(ICTYP(IQ).EQ.3) THEN
cc
C  out-of-plane bend constraint (I-J-K-L)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       L = IC(4,IQ)
       VAL = RCON(IQ)
C
C  check for equivalent atoms
C  these must also be constrained for symmetry
C  to be maintained
C  ** WARNING: Dummy atoms ignored **
C
       If(I.GT.NAtoms.OR.J.GT.NAtoms.OR.K.GT.NAtoms.OR.L.GT.NAtoms) Then
        WRITE(IOut,1150) I,J,K,L
        GO TO 80
       EndIf
cc
       DO 30 IOP=2,NTrans
        IT = NEqATM(I,IOP)
        JT = NEqATM(J,IOP)
        KT = NEqATM(K,IOP)
        LT = NEqATM(L,IOP)
cc
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K.AND.LT.EQ.L).OR.
     $     (IT.EQ.I.AND.KT.EQ.J.AND.JT.EQ.K.AND.LT.EQ.L) ) GO TO 30
C
C  check to see if this constraint is present
C
       DO 29 JQ=1,NCons
       IF(ICTYP(JQ).EQ.3) THEN
        IS = IC(1,JQ)
        JS = IC(2,JQ)
        KS = IC(3,JQ)
        LS = IC(4,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS.AND.KT.EQ.KS.AND.LT.EQ.LS).OR.
     $       (IT.EQ.IS.AND.KT.EQ.JS.AND.JT.EQ.KS.AND.LT.EQ.LS))
     $          .AND.Abs(RCON(JQ)).EQ.Abs(VAL) ) GO TO 30
cc note: temporary fix since oop bends have signature and original
cc       test may erroneously flag angles with same magnitude but
cc       opposite sign (however test no longer foolproof)
       ENDIF
 29    CONTINUE
C
C  if we reach here the symmetry-required constraint is absent
C
       IErr = -1
       WRITE(IOut,1200) IT,JT,KT,LT,I,J,K,L
C
 30    CONTINUE
cc
      ELSE IF(ICTYP(IQ).EQ.4) THEN
cc
C  dihedral constraint (I-J-K-L)
C
       I = IC(1,IQ)
       J = IC(2,IQ)
       K = IC(3,IQ)
       L = IC(4,IQ)
       VAL = RCON(IQ)
C
C  check for equivalent atoms
C  these must also be constrained for symmetry
C  to be maintained
C  ** WARNING: Dummy atoms ignored **
C
       If(I.GT.NAtoms.OR.J.GT.NAtoms.OR.K.GT.NAtoms.OR.L.GT.NAtoms) Then
        WRITE(IOut,1250) I,J,K,L
        GO TO 80
       EndIf
cc
       DO 40 IOP=2,NTrans
        IT = NEqATM(I,IOP)
        JT = NEqATM(J,IOP)
        KT = NEqATM(K,IOP)
        LT = NEqATM(L,IOP)
cc
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K.AND.LT.EQ.L).OR.
     $     (LT.EQ.I.AND.KT.EQ.J.AND.JT.EQ.K.AND.IT.EQ.L) ) GO TO 40
C
C  check to see if this constraint is present
C
       DO 39 JQ=1,NCons
       IF(ICTYP(JQ).EQ.4) THEN
        IS = IC(1,JQ)
        JS = IC(2,JQ)
        KS = IC(3,JQ)
        LS = IC(4,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS.AND.KT.EQ.KS.AND.LT.EQ.LS).OR.
     $       (LT.EQ.IS.AND.KT.EQ.JS.AND.JT.EQ.KS.AND.IT.EQ.LS))
     $          .AND.Abs(RCON(JQ)).EQ.Abs(VAL) ) GO TO 40
cc note: temporary fix since dihedrals have signature and original
cc       test may erroneously flag angles with same magnitude but
cc       opposite sign (however test no longer foolproof)
       ENDIF
 39    CONTINUE
C
C  if we reach here the symmetry-required constraint is absent
C
       IErr = -1
       WRITE(IOut,1300) IT,JT,KT,LT,I,J,K,L
C
 40    CONTINUE
cc
      ENDIF
c
 80   CONTINUE
C
      IF(IErr.NE.0) THEN
       WRITE(IOut,1400)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(' Distances ',2I4,8X,' and ',2I4,8X,' Should be Identical')
 1100 FORMAT(' Angles    ',3I4,4X,' and ',3I4,4X,' Should be Identical')
 1150 FORMAT('**WARNING** OOP bend Constraint ',4I4,' Involves Dummy',
     $   ' Atoms',/'  Check Carefully for Symmetry-Related Constraints')
 1200 FORMAT(' OOP bends ',4I4,   ' and ',4I4,   ' Should be Identical')
 1250 FORMAT('**WARNING** Dihedral Constraint ',4I4,' Involves Dummy',
     $   ' Atoms',/'  Check Carefully for Symmetry-Related Constraints')
 1300 FORMAT(' Dihedrals ',4I4,   ' and ',4I4,   ' Should be Identical')
 1400 FORMAT(/,2X,'***ERROR*** Imposed Constraints Insufficient to',
     $              ' Preserve Symmetry')
c
      END
c =======================================================================
c
      SUBROUTINE SymDRIVE(NAtoms, NDrive, IDRTYP,  FDRIVE, IDRIVE,
     $                    NTrans, NEqATM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks whether the imposed coordinate driving will preserve the
C  molecular symmetry.  If not, prints out the additional driving
C  required and Exits
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NDrive  -  number of primitives to drive
C  IDRTYP  -  coordinate type
C              1 - fixed distance
C              2 - fixed bond angle
C              3 - fixed out-of-plane-bend
C              4 - fixed dihedral angle
C  FDRIVE  -  driving force for each coordinate
C  IDRIVE  -  atoms involved in primitives being driven
C               ID1-ID2           distance constraint
C               ID1-ID2-ID3       bond angle constraint
C               ID1-ID2-ID3-ID4   out-of-plane bend or dihedral constraint
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C
      DIMENSION IDRTYP(NDrive),FDRIVE(NDrive),IDRIVE(4,NDrive)
      DIMENSION NEqATM(NAtoms,NTrans)
C
C
      If(NTrans.EQ.1) RETURN     ! there is no symmetry to conserve
c
      IOut = ioutfil('iout')
      IErr = 0
C
      DO 80 IQ=1,NDrive
c
      IF(IDRTYP(IQ).EQ.1) THEN
cc
C  driving a distance (I-J)
C
       I = IDRIVE(1,IQ)
       J = IDRIVE(2,IQ)
       VAL = FDRIVE(IQ)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C
       DO 10 IOP=2,NTrans
       IF(I.GT.NAtoms) THEN
        IT = I
       ELSE
        IT = NEqATM(I,IOP)
       ENDIF
       IF(J.GT.NAtoms) THEN
        JT = J
       ELSE
        JT = NEqATM(J,IOP)
       ENDIF
       If( (IT.EQ.I.AND.JT.EQ.J).OR.
     $     (JT.EQ.I.AND.IT.EQ.J) ) GO TO 10
C
C  check to see if this coordinate is being driven
C
       DO 9 JQ=1,NDrive
       IF(IDRTYP(JQ).EQ.1) THEN
        IS = IDRIVE(1,JQ)
        JS = IDRIVE(2,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS).OR.(JT.EQ.IS.AND.IT.EQ.JS))
     $          .AND.FDRIVE(JQ).EQ.VAL ) GO TO 10
       ENDIF
 9     CONTINUE
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1000) IT,JT,I,J
C
 10    CONTINUE
cc
      ELSE IF(IDRTYP(IQ).EQ.2) THEN
cc
C  driving an angle (I-J-K)
C
       I = IDRIVE(1,IQ)
       J = IDRIVE(2,IQ)
       K = IDRIVE(3,IQ)
       VAL = FDRIVE(IQ)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C
       DO 20 IOP=2,NTrans
       IF(I.GT.NAtoms) THEN
        IT = I
       ELSE
        IT = NEqATM(I,IOP)
       ENDIF
       IF(J.GT.NAtoms) THEN
        JT = J
       ELSE
        JT = NEqATM(J,IOP)
       ENDIF
       IF(K.GT.NAtoms) THEN
        KT = K
       ELSE
        KT = NEqATM(K,IOP)
       ENDIF
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K).OR.
     $     (KT.EQ.I.AND.JT.EQ.J.AND.IT.EQ.K) ) GO TO 20
C
C  check to see if this coordinate is being driven
C
       DO 19 JQ=1,NDRive
       IF(IDRTYP(JQ).EQ.2) THEN
        IS = IDRIVE(1,JQ)
        JS = IDRIVE(2,JQ)
        KS = IDRIVE(3,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS.AND.KT.EQ.KS).OR.
     $       (KT.EQ.IS.AND.JT.EQ.JS.AND.IT.EQ.KS))
     $          .AND.FDRIVE(JQ).EQ.VAL ) GO TO 20
       ENDIF
 19    CONTINUE
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1100) IT,JT,KT,I,J,K
C
 20    CONTINUE
cc
      ELSE IF(IDRTYP(IQ).EQ.3) THEN
cc
C  driving an out-of-plane bend (I-J-K-L)
C
       I = IDRIVE(1,IQ)
       J = IDRIVE(2,IQ)
       K = IDRIVE(3,IQ)
       L = IDRIVE(4,IQ)
       VAL = FDRIVE(IQ)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C  ** WARNING: Dummy atoms ignored **
C
       If(I.GT.NAtoms.OR.J.GT.NAtoms.OR.K.GT.NAtoms.OR.L.GT.NAtoms) Then
        WRITE(IOut,1150) I,J,K,L
        GO TO 80
       EndIf
cc
       DO 30 IOP=2,NTrans
        IT = NEqATM(I,IOP)
        JT = NEqATM(J,IOP)
        KT = NEqATM(K,IOP)
        LT = NEqATM(L,IOP)
cc
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K.AND.LT.EQ.L).OR.
     $     (IT.EQ.I.AND.KT.EQ.J.AND.JT.EQ.K.AND.LT.EQ.L) ) GO TO 30
C
C  check to see if this coordinate is being driven
C
       DO 29 JQ=1,NDrive
       IF(IDRTYP(JQ).EQ.3) THEN
        IS = IDRIVE(1,JQ)
        JS = IDRIVE(2,JQ)
        KS = IDRIVE(3,JQ)
        LS = IDRIVE(4,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS.AND.KT.EQ.KS.AND.LT.EQ.LS).OR.
     $       (IT.EQ.IS.AND.KT.EQ.JS.AND.JT.EQ.KS.AND.LT.EQ.LS))
     $          .AND.Abs(FDRIVE(JQ)).EQ.Abs(VAL) ) GO TO 30
cc note: temporary fix since oop bends have signature and original
cc       test may erroneously flag angles with same magnitude but
cc       opposite sign (however test no longer foolproof)
       ENDIF
 29    CONTINUE
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1200) IT,JT,KT,LT,I,J,K,L
C
 30    CONTINUE
cc
      ELSE IF(IDRTYP(IQ).EQ.4) THEN
cc
C  dihedral constraint (I-J-K-L)
C
       I = IDRIVE(1,IQ)
       J = IDRIVE(2,IQ)
       K = IDRIVE(3,IQ)
       L = IDRIVE(4,IQ)
       VAL = FDRIVE(IQ)
C
C  check for equivalent atoms
C  these must also be driven for symmetry
C  to be maintained
C  ** WARNING: Dummy atoms ignored **
C
       If(I.GT.NAtoms.OR.J.GT.NAtoms.OR.K.GT.NAtoms.OR.L.GT.NAtoms) Then
        WRITE(IOut,1250) I,J,K,L
        GO TO 80
       EndIf
cc
       DO 40 IOP=2,NTrans
        IT = NEqATM(I,IOP)
        JT = NEqATM(J,IOP)
        KT = NEqATM(K,IOP)
        LT = NEqATM(L,IOP)
cc
       If( (IT.EQ.I.AND.JT.EQ.J.AND.KT.EQ.K.AND.LT.EQ.L).OR.
     $     (LT.EQ.I.AND.KT.EQ.J.AND.JT.EQ.K.AND.IT.EQ.L) ) GO TO 40
C
C  check to see if this coordinate is being driven
C
       DO 39 JQ=1,NDrive
       IF(IDRTYP(JQ).EQ.4) THEN
        IS = IDRIVE(1,JQ)
        JS = IDRIVE(2,JQ)
        KS = IDRIVE(3,JQ)
        LS = IDRIVE(4,JQ)
        If( ((IT.EQ.IS.AND.JT.EQ.JS.AND.KT.EQ.KS.AND.LT.EQ.LS).OR.
     $       (LT.EQ.IS.AND.KT.EQ.JS.AND.JT.EQ.KS.AND.IT.EQ.LS))
     $          .AND.Abs(FDRIVE(JQ)).EQ.Abs(VAL) ) GO TO 40
cc note: temporary fix since dihedrals have signature and original
cc       test may erroneously flag angles with same magnitude but
cc       opposite sign (however test no longer foolproof)
       ENDIF
 39    CONTINUE
C
C  if we reach here the symmetry-required coordinate is absent
C
       IErr = -1
       WRITE(IOut,1300) IT,JT,KT,LT,I,J,K,L
C
 40    CONTINUE
cc
      ENDIF
c
 80   CONTINUE
C
      IF(IErr.NE.0) THEN
       WRITE(IOut,1400)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(' Distances ',2I4,8X,' and ',2I4,8X,' Should be Identical')
 1100 FORMAT(' Angles    ',3I4,4X,' and ',3I4,4X,' Should be Identical')
 1150 FORMAT('**WARNING** OOP bend Driving ',4I4,' Involves Dummy',
     $   ' Atoms',/'  Check Carefully for Symmetry-Related Driving')
 1200 FORMAT(' OOP bends ',4I4,   ' and ',4I4,   ' Should be Identical')
 1250 FORMAT('**WARNING** Dihedral Driving ',4I4,' Involves Dummy',
     $   ' Atoms',/'  Check Carefully for Symmetry-Related Driving')
 1300 FORMAT(' Dihedrals ',4I4,   ' and ',4I4,   ' Should be Identical')
 1400 FORMAT(/,2X,'***ERROR*** Coordinate Driving Insufficient to',
     $              ' Preserve Symmetry')
c
      END
c =======================================================================
c
      SUBROUTINE SymFIX(NAtoms,IFix,NTrans,NEqATM)
      IMPLICIT INTEGER(A-Z)
C
C
C  Checks whether fixed atom constraints will preserve the
C  molecular symmetry.  If not, prints out the additional
C  constraints required and Exits
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IFix    -  list of fixed/active coordinates
C              0 - coordinate active
C              1 - coordinate will be fixed during optimization
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C
C
      DIMENSION IFix(3,NAtoms),NEqATM(NAtoms,NTrans)
      CHARACTER*1 CFix(3)
C
      DATA CFix/'X','Y','Z'/
C
C
      IOut = ioutfil('iout')
c
      If(NTrans.EQ.1) RETURN     ! there is no symmetry to conserve
      IErr = 0
C
      DO 20 I=1,NAtoms
      DO 20 J=1,3
      IF(IFix(J,I).EQ.1) THEN
C
C  we have a fixed coordinate
C  check for equivalent atoms
C  these must also be fixed for symmetry
C  to be maintained
C
       DO 10 IOP=2,NTrans
       IT = NEqATM(I,IOP)
       IF(IT.NE.I.AND.IFix(J,IT).NE.1) THEN
        IErr = -1
        WRITE(IOut,1000) IT,CFix(J)
       ENDIF
 10    CONTINUE
c
      ENDIF
 20   CONTINUE
C
      IF(IErr.NE.0) THEN
       WRITE(IOut,1100)
       CALL OptExit(9)
      ENDIF
C
      RETURN
c
 1000 FORMAT(' To Preserve Symmetry  Atom ',I3,2X,A1,' Coordinate',
     $       ' Should be Fixed')
 1100 FORMAT(/,2X,'***ERROR*** Fixing Only the Requested Atomic',
     $              ' Coordinates will Break Symmetry')
c
      END
c =======================================================================
c
      SUBROUTINE SymHES(NAtoms, NTrans, NEqATM, TRANS,  HOld,
     $                  R,      U,      thrsh,  HNew )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Symmetrizes a Cartesian Hessian matrix according to all
C  operations of the molecular point group
C  Additionally zeros all elements below thrsh
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  HOld    -  matrix (Hessian) to be symmetrized
C  R       -  scratch space (for transformation matrix)
C  U       -  scratch space
C  thrsh   -  threshold below which elements will be set to zero
C  HNew    -  on exit contains symmetrized matrix
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),
     $          HOld(3*NAtoms,3*NAtoms),R(3*NAtoms,3*NAtoms),
     $          U(3*NAtoms,3*NAtoms),HNew(3*NAtoms,3*NAtoms)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
C
C
      NAT3 = 3*NAtoms
      CALL CpyVEC(NAT3*NAT3,HOld,HNew)
      If(NTrans.EQ.1) GO TO 75
c
      DO 70 IOP=2,NTrans
C
C  Construct the R matrix for each symmetry operation
C
      CALL ZeroIT(R,NAT3*NAT3)
c
      DO 20 IAtm=1,NAtoms
      II = (IAtm-1)*3
      JAtm = NEqATM(IAtm,IOP)
      JJ = (JAtm-1)*3
c
      DO 10 K=1,3
      KK = II+K
      DO 10 L=1,3
      LL = JJ+L
      R(KK,LL) = TRANS(L,K,IOP)
 10   CONTINUE
c
 20   CONTINUE
C
C  Form R * HESS * R(t)
C
      CALL ZeroIT(U,NAT3*NAT3)
c
      DO 40 J=1,NAT3
      DO 40 K=1,NAT3
      VAL = HOld(K,J)
      DO 30 I=1,NAT3
      U(I,J) = U(I,J) + R(I,K)*VAL
 30   CONTINUE
 40   CONTINUE
c
      DO 60 J=1,NAT3
      DO 60 K=1,NAT3
      VAL = R(J,K)
      DO 50 I=1,NAT3
      HNew(I,J) = HNew(I,J) + U(I,K)*VAL
 50   CONTINUE
 60   CONTINUE
c
 70   CONTINUE
C
C  scale the final matrix
C
      skal = One/DFloat(NTrans)
      CALL Vscal(NAT3*NAT3,skal,HNew)
C
 75   CONTINUE
C
C  zero out elements below thrsh
C
      DO 80 I=1,NAT3
      DO 80 J=1,I
      If(Abs(HNew(I,J)).LT.thrsh) HNew(I,J) = Zero
      HNew(J,I) = HNew(I,J)
 80   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE SymNIC(intcor, NVib,   NDEG,   NAT3,   thrsh,
     $                  IPRNT,  UT,     GINT,   XINT,   IHess,
     $                  BINV,   IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Remove symmetry-redundant composite internals
C
C  ARGUMENTS
C
C  intcor  -  total number of primitive internals generated
C  NVib    -  number of non-redundant natural internals
C             (in general case this should be 3*NATOMS-6)
C  NDEG    -  number of degrees of freedom accounting
C             for all symmetry
C  NAT3    -  3 * number of atoms (leading dimension of BINV)
C  thrsh   -  threshold for elimination of natural internals
C             (if gradient below thrsh - eliminate)
C  IPRNT   -  print flag
C  UT      -  transformation matrix
C             (vectors of linear coefficients relating
C              primitive and natural internals)
C  GINT    -  gradient over NVib natural internals
C  XINT    -  current internal coordinate values
C  IHess   -  flag controlling Hessian transformation
C             if Hessian needs to be transformed we need to
C             eliminate symmetry-redundant columns of BINV
C  BINV    -  inverse B-Matrix
C  IErr    -  error flag
C              0 - success; eliminated appropriate number or
C                  natural internals due to symmetry
C             -1 - failure; wrong number eliminated
C
      REAL*8 GINT(NVib),XINT(NVib),UT(intcor,NVib)
      REAL*8 BINV(NAT3,NVib)
C
C
      IOut = ioutfil('iout')
C
C  eliminate all natural internals with gradient below thrsh
C
      IT = 0
      DO 10 I=1,NVib
      Val = Abs(GINT(I))
      IF(Val.GT.thrsh) THEN
       IT = IT + 1
       CALL CpyVEC(intcor,UT(1,I),UT(1,IT))
       If(IHess.EQ.0) CALL CpyVEC(NAT3,BINV(1,I),BINV(1,IT))
       XINT(IT) = XINT(I)
       GINT(IT) = GINT(I)
      ENDIF
 10   CONTINUE
c
      If(IPRNT.GT.2) WRITE(IOut,1000) NVib-IT
c
      IF(IT.NE.NDEG) THEN
       WRITE(IOut,1100) NDEG,IT
cc       If(IT.LT.NDEG) IErr = -1
       NDEG = IT
       RETURN
      ENDIF
c
      If(IPRNT.GT.4) THEN
       WRITE(IOut,1200)
       CALL PrntMAT(NDEG,intcor,NDEG,UT)
      EndIF
C
      RETURN
c
 1000 FORMAT(' Eliminated ',I5,' Coordinates due to Symmetry')
 1100 FORMAT('**WARNING** Wrong Number of Internal Coordinates',/,
     $       '  Should be: ',I3,'  Found: ',I3)
 1200 FORMAT(/,'   Final Set of Delocalized Internal Coordinates')
c
      END
c =======================================================================
c
      SUBROUTINE SymPRIM(NPrim,IG,XPrim)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Sets symmetry-related primitives to same values
C  **WARNING**  Actually this is not implemented properly as
C    values are averaged pairwise as they are found instead
C    of finding them all and then averaging
C    However, the iterative procedure should take care of this
C
C
C  ARGUMENTS
C
C  NPrim   -  number of primitive internals
C  IG      -  primitive symmetry equivalences
C                 0 - parent primitive
C                 J - has same value as previous (Jth) primitive
C                -J - has same value, opposite sign
C  XPrim   -  primitive values
C
C
      DIMENSION IG(NPrim),XPrim(NPrim)
C
      PARAMETER (Half=0.5d0)
C
C
      DO 10 I=1,NPrim
      II = IG(I)
      If(II.GT.0) Then
       Val = Half*(XPrim(I) + XPrim(II))
       XPrim(I) = Val
       XPrim(II) = Val
      Else If(II.LT.0) Then
       Val = Half*(XPrim(I) - XPrim(-II))
       XPrim(I) =  Val
       XPrim(-II) = -Val
      EndIf
 10   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE TidyUP(NDEG,   intcor, thrsh,  IPRNT,  UTT,
     $                  Coeff,  ktyp,   klist,  SavTOR, UT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Remove all internals from the primitive space with
C  weights less than thrsh
C
C  ARGUMENTS
C
C  NDEG    -  size of active space
C  intcor  -  on input current size of primitive space
C             on output new size of space
C  thrsh   -  weight threshold for eliminating primitive
C  IPRNT   -  print flag
C  UTT     -  scratch array for transpose of UT
C  Coeff   -  weighting of primitives in active space
C  ktyp    -   array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C  klist   -  list of atoms involved in each primitive internal
C  SavTOR  -  array for storing primitive torsions
C             (possible sign changes near limiting values)
C  UT      -  transformation matrix
C              i.e. which linear combination of primitive internals
C                   make up each compound (natural) internal coordinate
C
C
      REAL*8 Coeff(intcor),SavTOR(intcor),UT(intcor*NDEG),
     $       UTT(NDEG,intcor)
      DIMENSION ktyp(intcor),klist(4,intcor)
      DIMENSION ityp(6)
C
      DATA ityp/6*0/
C
C
C  first check if we actually have anything to eliminate
C
      DO 10 I=1,intcor
      If(Coeff(I).LT.thrsh) GO TO 20
 10   CONTINUE
C
C  at this point every primitive should be kept
C  so why are we still here?
C
      RETURN
c
 20   CONTINUE
C
C  initialize
C
      II = 0
C
C now check weights for elimination
C
      DO 40 I=1,intcor
      IF(Coeff(I).GE.thrsh) THEN
C
C  keep primitive coordinate
C
       II = II+1
       IF(I.NE.II) THEN
        ktyp(II) = ktyp(I)
        klist(1,II) = klist(1,I)
        klist(2,II) = klist(2,I)
        klist(3,II) = klist(3,I)
        klist(4,II) = klist(4,I)
        Coeff(II) = Coeff(I)
        SavTOR(II) = SavTOR(I)
       ENDIF
C
C  take care with UT as leading dimension is changing
C  copy into UTT array as a transpose
C
       DO 30 J=1,NDEG
       IJ = (J-1)*intcor + I
       UTT(J,II) = UT(IJ)
 30    CONTINUE
cc
      ELSE
C
C  eliminate it
C  (keep tabs for print out)
C
       IT = ktyp(I)
       ityp(IT) = ityp(IT) + 1
cc
      ENDIF
 40   CONTINUE
C
C  reset number of primitives
C
      intcor = II
C
C  restore revised internal coordinates
C
      DO 60 I=1,NDEG
      II = (I-1)*intcor
      DO 50 J=1,intcor
      UT(II+J) = UTT(I,J)
 50   CONTINUE
 60   CONTINUE
c
      If(IPRNT.GT.2) Then
       IOut = ioutfil('iout')
       WRITE(IOut,1000)
       If(ityp(1).GT.0) WRITE(IOut,1100) ityp(1)
       If(ityp(2).GT.0) WRITE(IOut,1200) ityp(2)
       If(ityp(3).GT.0) WRITE(IOut,1300) ityp(3)
       If(ityp(4).GT.0) WRITE(IOut,1400) ityp(4)
       If(ityp(5).GT.0) WRITE(IOut,1500) ityp(5)
       If(ityp(6).GT.0) WRITE(IOut,1600) ityp(6)
       WRITE(IOut,1700) intcor
      EndIf
c
      RETURN
c
 1000 FORMAT(/,' Eliminating Redundant Primitive Internals from Space')
 1100 FORMAT('  Removed ',I6,' Stretches')
 1200 FORMAT('  Removed ',I6,' Bends')
 1300 FORMAT('  Removed ',I6,' Out-of-Plane Bends')
 1400 FORMAT('  Removed ',I6,' Torsions')
 1500 FORMAT('  Removed ',I6,' Colinear bends')
 1600 FORMAT('  Removed ',I6,' Perpendicular bends')
 1700 FORMAT(' There are now ',I6,' Primitive Internals')
c
      END
c =======================================================================
c
      SUBROUTINE TOPOLOGY(NAtoms, NMol,   IMOL,   ITors,  NIC,
     $                    GROUP,  NQ,     IC,     IPRNT,  intcor,
     $                    ktyp,   klist)
      IMPLICIT INTEGER(A-Z)
C
C
C  Generates bond stretches and primitive bends and torsions
C  from the atomic connectivity matrix
C  ** TEMPORARY BRUTE FORCE ALGORITHM **
C  ** NEED TO PUT SAFETY CHECKS IN    **
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NMol    -  number of molecules (for, e.g., cluster optimizations)
C  IMOL    -  pointers to start/end of molecules
C  ITors   -  use of primitive torsions
C               0 - use torsions throughout
C               1 - do not use torsions for FIRST molecule only
C             (in surface optimizations, the first "molecule" is the surface)
C  NIC     -  maximum number of primitive internal coordinates
C             that can be generated
C             ** NOTE:  This is typically far more than 3*NATOMS-6
C                       as many redundancies are found that are
C                       dealt with later
C  GROUP   -  molecular point group
C  NQ      -  number of symmetry unique atoms
C  IC      -  connectivity matrix
C  IPRNT   -  print flag
C
C  on exit
C
C  intcor  -  number of primitive internals found
C  ktyp    -  integer array containing internal coordinate type
C             ** NOTE:  These types are:
C                       1 - stretch
C                       2 - bend
C                       3 - out-of-plane bend
C                       4 - torsion
C  klist   -  list of atoms involved in each primitive
C
C
      DIMENSION IMOL(NAtoms),IC(NAtoms,NAtoms)
      DIMENSION ktyp(NIC),klist(4,NIC)
      CHARACTER*4 GROUP
      LOGICAL Bend,Outp,Tors
C
      DATA Bend/.TRUE./, Outp/.FALSE./
C
C
C  initialize
C
      IOut = ioutfil('iout')
      CALL IZeroIT(klist,4*NIC)
      intcor = 0
c
      If(IPRNT.GT.3) WRITE(IOut,1000)
C
C
C  Loop over all molecules
C  generate primitives for each molecule separately
C
      DO 70 Mol=1,NMol
      IStrt = IMOL(Mol)+1
      IEnd  = IMOL(Mol+1)

C  whether or not to use torsions

      If(ITors.EQ.1.AND.Mol.EQ.1) Then
       Tors = .FALSE.
      Else
       Tors = .TRUE.
      EndIf
C
C  (a) Stretches
C
      int0 = intcor
      DO 10 I=IStrt+1,IEnd
      DO 10 J=IStrt,I-1
      IF(IC(I,J).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 1
       klist(1,intcor) = I
       klist(2,intcor) = J
      ENDIF
 10   CONTINUE
c
      If(IPRNT.GT.4.AND.intcor.GT.int0) Then
       WRITE(IOut,1100)
       Do I=int0+1,intcor
       WRITE(IOut,1200) I,klist(1,I),klist(2,I)
       EndDo
      EndIf
c
      If(IPRNT.GT.3) WRITE(IOut,1300) intcor-int0
c
      If(.NOT.Bend) GO TO 25
C
C  (b) Bends
C
      int0 = intcor
      DO 20 I=IStrt+2,IEnd
      DO 20 J=IStrt+1,I-1
      DO 20 K=IStrt,J-1
c
      IF(IC(I,J).NE.0.AND.IC(J,K).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 2
       klist(1,intcor) = I
       klist(2,intcor) = K
       klist(3,intcor) = J
      ENDIF
c
      IF(IC(J,K).NE.0.AND.IC(K,I).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 2
       klist(1,intcor) = J
       klist(2,intcor) = I
       klist(3,intcor) = K
      ENDIF
c
      IF(IC(K,I).NE.0.AND.IC(I,J).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 2
       klist(1,intcor) = K
       klist(2,intcor) = J
       klist(3,intcor) = I
      ENDIF
c
 20   CONTINUE
c
      If(IPRNT.GT.4.AND.intcor.GT.int0) Then
       WRITE(IOut,1400)
       Do I=int0+1,intcor
       WRITE(IOut,1200) I,klist(1,I),klist(3,I),klist(2,I)
       EndDo
      EndIf
c
      If(IPRNT.GT.3) WRITE(IOut,1500) intcor-int0
c
 25   CONTINUE
C
      IF(.NOT.Outp) THEN            ! usually skip out-of-plane bends
       If(NAtoms.NE.4) GO TO 35
       If( .NOT.(GROUP.EQ.'cs  '.AND.NQ.EQ.4) .AND.
     $     .NOT.(GROUP.EQ.'c2v ') .AND.
     $     .NOT.(GROUP.EQ.'d3h ') ) GO TO 35
C
C -- Special Case   Planar molecule --
C    With 3 atoms attached to a central atom there will
C    be no proper torsions
C    Use out-of-plane-bend instead
C
       I = 0
       J = 0
       K = 0
       L = 0
c
       DO 27 II=1,4
       I = I + IC(1,II)
       J = J + IC(2,II)
       K = K + IC(3,II)
       L = L + IC(4,II)
 27    CONTINUE
c
       If( .NOT.(I.EQ.3.AND.J.EQ.1.AND.K.EQ.1.AND.L.EQ.1) .AND.
     $     .NOT.(I.EQ.1.AND.J.EQ.3.AND.K.EQ.1.AND.L.EQ.1) .AND.
     $     .NOT.(I.EQ.1.AND.J.EQ.1.AND.K.EQ.3.AND.L.EQ.1) .AND.
     $     .NOT.(I.EQ.1.AND.J.EQ.1.AND.K.EQ.1.AND.L.EQ.3) ) GO TO 35
      ENDIF
C
C  (c) Out-of-Plane Bends
C
      int0 = intcor
      DO 30 I=IStrt+3,IEnd
      DO 30 J=IStrt+2,I-1
      DO 30 K=IStrt+1,J-1
      DO 30 L=IStrt,K-1
c
      IF(IC(I,J).NE.0.AND.IC(I,K).NE.0.AND.IC(I,L).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = J
       klist(2,intcor) = K
       klist(3,intcor) = L
       klist(4,intcor) = I
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = K
       klist(2,intcor) = L
       klist(3,intcor) = J
       klist(4,intcor) = I
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = L
       klist(2,intcor) = J
       klist(3,intcor) = K
       klist(4,intcor) = I
      ENDIF
c
      IF(IC(J,I).NE.0.AND.IC(J,K).NE.0.AND.IC(J,L).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = I
       klist(2,intcor) = K
       klist(3,intcor) = L
       klist(4,intcor) = J
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = K
       klist(2,intcor) = L
       klist(3,intcor) = I
       klist(4,intcor) = J
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = L
       klist(2,intcor) = I
       klist(3,intcor) = K
       klist(4,intcor) = J
      ENDIF
c
      IF(IC(K,I).NE.0.AND.IC(K,J).NE.0.AND.IC(K,L).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = I
       klist(2,intcor) = J
       klist(3,intcor) = L
       klist(4,intcor) = K
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = J
       klist(2,intcor) = L
       klist(3,intcor) = I
       klist(4,intcor) = K
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = L
       klist(2,intcor) = I
       klist(3,intcor) = J
       klist(4,intcor) = K
      ENDIF
c
      IF(IC(L,I).NE.0.AND.IC(L,J).NE.0.AND.IC(L,K).NE.0) THEN
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = I
       klist(2,intcor) = J
       klist(3,intcor) = K
       klist(4,intcor) = L
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = J
       klist(2,intcor) = K
       klist(3,intcor) = I
       klist(4,intcor) = L
c
       intcor = intcor+1
       ktyp(intcor) = 3
       klist(1,intcor) = K
       klist(2,intcor) = I
       klist(3,intcor) = J
       klist(4,intcor) = L
      ENDIF
c
 30   CONTINUE
c
      If(IPRNT.GT.4.AND.intcor.GT.int0) Then
       WRITE(IOut,1600)
       Do I=int0+1,intcor
       WRITE(IOut,1200) I,klist(1,I),klist(2,I),klist(3,I),klist(4,I)
       EndDo
      EndIf
c
      If(IPRNT.GT.3) WRITE(IOut,1700) intcor-int0
C
 35   CONTINUE
C
      If(.NOT.Tors) GO TO 45
C
C  (d) Torsions
C
      int0 = intcor
      DO 42 I=IStrt+1,IEnd
      DO 42 J=IStrt,I-1
      IF(IC(I,J).NE.0) THEN
C
C  consider I-J as middle 2 atoms in proper torsion
C
       DO 41 K=IStrt,IEnd
       IF(IC(I,K).NE.0.AND.K.NE.J) THEN
C
C  have K-I-J
C
        DO 40 L=IStrt,IEnd
        IF(IC(J,L).NE.0.AND.I.NE.L.AND.K.NE.L) THEN
         intcor = intcor+1
         ktyp(intcor) = 4
         klist(1,intcor) = K
         klist(2,intcor) = I
         klist(3,intcor) = J
         klist(4,intcor) = L
        ENDIF
 40     CONTINUE
       ENDIF
 41    CONTINUE
      ENDIF
c
 42   CONTINUE
c
      If(IPRNT.GT.4.AND.intcor.GT.int0) Then
       WRITE(IOut,1800)
       Do I=int0+1,intcor
       WRITE(IOut,1200) I,klist(1,I),klist(2,I),klist(3,I),klist(4,I)
       EndDo
      EndIf
c
      If(IPRNT.GT.3) WRITE(IOut,1900) intcor-int0
C
 45   CONTINUE
C
C -- End loop over molecules
C
 70   CONTINUE
c
      If(IPRNT.GT.2) WRITE(IOut,2000) intcor
C
      RETURN
c
 1000 FORMAT(' Generating Primitive Internal Coordinates')
 1100 FORMAT(' Primitive Stretches:')
 1200 FORMAT(5X,I4,2X,4I5)
 1300 FORMAT(' There are ',I6,' Stretches')
 1400 FORMAT(' Primitive Bends:')
 1500 FORMAT(' There are ',I6,' Bends')
 1600 FORMAT(' Primitive Out-of-Plane Bends')
 1700 FORMAT(' There are ',I6,' Out-of-Plane Bends')
 1800 FORMAT(' Primitive Torsions:')
 1900 FORMAT(' There are ',I6,' Torsions')
 2000 FORMAT(' There are ',I6,' Primitive Internals')
c
      END
c =======================================================================
c
      SUBROUTINE Tor2ZN(NAtoms, NPrim,  IGen,   ktyp,   klist,
     $                  XPrim,  XC,     IPRNT,  IVec,   Zindx,
     $                  ZVals,  ZMap,   INDX,   IOrd)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms a Z matrix for the iterative back-transformation
C
C  DELOCALIZED INTERNALS
C  ---------------------
C  Keys off torsions. By construction, all torsions are
C  present in primitive set, so should always work.
C
C
C  NATURAL INTERNALS
C  -----------------
C  Keys off out-of-plane bends. If a suitable out-of-plane bend
C  cannot be found, tries first bend-pairs and then torsions.
C  If no torsions can be found, fails.
C
C .............................................................
C  NOTE:  Formerly this routine keyed off torsions only.
C         This was changed to accomodate natural internals
C         which often do not have all the necessary torsions
C         in the primitive space
C .............................................................
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NPrim   -  number of primitives
C  IGen    -  how to form Z matrix (coordinate type)
C              -1 keys off bend-pairs, out-of-plane bends & torsions
C                 in that order (natural internals)
C              otherwise keys off torsions only (delocalized)
C  ktyp    -  array indicating each primitive internal type
C               currently these are
C                1 - stretch
C                2 - bend
C                3 - out-of-plane bend
C                4 - torsion
C                5 - linear coplanar bend
C                6 - linear perpendicular bend
C                7 - 1/R distance coordinates
C             **WARNING** only types 1-4 are handled here
C  klist   -  list of atoms involved in each primitive internal
C                I  J  0  0      stretch I-J
C                I  K  J  0      bend    I-J-K
C                I  K  J  L      torsion I-J-K-L
C                            or out-of-plane bend
C             **WARNING** primitives are assumed to be in order
C               i.e.  all stretches, then all bends and so on
C  XPrim   -  values of primitive internals
C  XC      -  Cartesian coordinates
C  IPRNT   -  print flag
C  IVec    -  integer work array to keep track of which primitives
C             are already accounted for during Z-matrix construction
C
C  on exit
C
C  Zindx   -  index of atoms in Z-matrix
C             The first column of Zindx is the atom numbers; it
C             is needed for renumbering the atoms into canonical
C             order on exit
C              col 1 2      are bond indices
C              col 1 2 3    give angle indices
C              col 1 2 3 4  give either torsion, 2-bend or
C                            out-of-plane bend indices
C              column 5 is 0 if connectivity is via a standard torsion
C                          1 if connectivity is via two bends
C                          2 if connectivity is via out-of-plane bend
C              on exit col 1 is 1,2,3,...,n and is not needed
C               for the Z-matrix to Cartesian routine
C  ZVals    -  values of parameters in Z matrix
C  ZMap     -  canonical vertex mapping
C               e.g.  1  1
C                     3  2
C                     2  3
C               means that input atom 1 maps to output atom 1
C                          input atom 3 maps to output atom 2
C                          inout atom 2 maps to output atom 3
C               This is needed to get the Z-matrix in canonical form
C  INDX     -  indexing array indicating location in primitive list
C               of those primitives extracted to construct Z-matrix
C  IOrd     -  integer flag
C                -1 - something went wrong
C                 0 - success, no canonical reordering
C                 1 - success, atoms reordered in Z-matrix
C
C
      INTEGER ktyp(Nprim),klist(4,Nprim),IVec(NPrim),
     &	      Zindx(NAtoms,5),ZMap(NAtoms,2),INDX(3*NAtoms)
      REAL*8  XPrim(NPrim),XC(3,NAtoms),ZVals(NAtoms,3)
      Logical found
C
      PARAMETER (Zero=0.0d0,One=1.0d0,thrsh=0.1d0)
C
C
      IOut = ioutfil('iout')
C
C  first find out how many bends, torsions and out-of-plane bends
C  there are and where they start in the primitive list
C
      DO 10 I=1,NPrim
      If(ktyp(I).NE.1) Exit
 10   CONTINUE
      NB1 = I
      DO 11 I=NB1,NPrim
      If(ktyp(I).NE.2) Exit
 11   CONTINUE
      NB2 = I-1
      NP1 = I
      DO 12 I=NB2+1,NPrim
      If(Ktyp(I).NE.3) Exit
 12   CONTINUE
      NP2 = I-1
      NT1 = I
      DO 13 I=NT1,NPrim
      If(Ktyp(I).NE.4) Exit
 13   CONTINUE
      NT2 = I-1
C
C  bends              start at NB1, end at NB2
C  out-of-plane bends start at NP1, end at NP2
C  torsions           start at NT1, end at NT2
C
      IOrd = -1                 ! assume the worst
C
C ----------------------------------------------------------------
C  For delocalized internals, pick the first torsion in the
C  list and use it to construct the first 4 rows of the Z matrix.
C  For natural internals, pick the first bend and construct
C  the first 3 rows of the Z matrix
C ----------------------------------------------------------------
C
      CALL IZeroIT(Zindx,5*NAtoms)
      CALL IZeroIT(IVec,NPrim)
c
      IF(IGen.EQ.-1.OR.NAtoms.LT.4) THEN
C
C  start with first bend
C
       i1 = klist(1,NB1)
       i2 = klist(3,NB1)
       i3 = klist(2,NB1)
c
       Zindx(1,1) = i3
       Zindx(2,1) = i2
       Zindx(2,2) = i3
       Zindx(3,1) = i1
       Zindx(3,2) = i2
       Zindx(3,3) = i3
C
C  set initial entries in IVec
C
       IVec(i1) = 1
       IVec(i2) = 1
       IVec(i3) = 1
C
C  set counter (number of entries in Z matrix so far)
C
       IC = 3
cc
      ELSE
C
C  start with first torsion
C
       i1 = klist(1,NT1)
       i2 = klist(2,NT1)
       i3 = klist(3,NT1)
       i4 = klist(4,NT1)
c
       Zindx(1,1) = i4
       Zindx(2,1) = i3
       Zindx(2,2) = i4
       Zindx(3,1) = i2
       Zindx(3,2) = i3
       Zindx(3,3) = i4
       Zindx(4,1) = i1
       Zindx(4,2) = i2
       Zindx(4,3) = i3
       Zindx(4,4) = i4
C
C  set initial entries in IVec
C
       IVec(i1) = 1
       IVec(i2) = 1
       IVec(i3) = 1
       IVec(i4) = 1
C
C  set counter (number of entries in Z matrix so far)
C
       IC = 4
cc
      ENDIF
C  ---------------------------------------------------------
C
C
C  First look for an out-of-plane bend
C
 20   CONTINUE
      If(IC.EQ.NAtoms) GO TO 95
      If(IGen.NE.-1) GO TO 94      ! go direct to torsions
c
      DO 30 I=NP1,NP2
c
      i1 = klist(1,I)
      i2 = klist(2,I)
      i3 = klist(3,I)
      i4 = klist(4,I)
C
C  in order to be acceptable as the next line in the Z matrix
C  i1 must be a new atom, and the other three indices must
C  already be present in the Z matrix
C  i4 is always the central atom (by definition)
C
      IF(IVec(i1).EQ.0) THEN
       If(IVec(i2)+IVec(i3)+IVec(i4).EQ.3) Then
C
C  new Z-matrix entry found
C
        IC = IC+1
        IVec(i1) = 1
        Zindx(IC,1) = i1
        Zindx(IC,2) = i4
        Zindx(IC,3) = i2
        Zindx(IC,4) = i3
        Zindx(IC,5) = 2
        GO TO 20
       EndIf
      ENDIF
 30   CONTINUE
C
C  If we get here then we didn't find an out-of-plane bend
C  Look for a suitable bend-pair instead
C
      DO 50 I=NB1+1,NB2
c
      found = .False.
      i1 = klist(1,I)
      i2 = klist(3,I)
      i3 = klist(2,I)
C
C  is this bend acceptable?
C  bends must be I-J-K and I-J-L where I is a new atom and
C  the other three indices must be present in the Z matrix
C
      IF(IVec(i1).EQ.0.AND.(IVec(i2)+IVec(i3).EQ.2)) THEN
       found = .True.
      ELSE IF(IVec(i3).EQ.0.AND.(IVec(i1)+IVec(i2).EQ.2)) THEN
       found = .True.
       itemp = i1
       i1 = i3
       i3 = itemp
      ENDIF
c
      IF(found) THEN
C
C  one suitable bend found
C  do we have a partner?
C
       DO 40 J=NB1+1,NB2
       If(J.EQ.I) GO TO 40    ! don't check the bend just found
c
       j1 = klist(1,J)
       j2 = klist(3,J)
       j3 = klist(2,J)
C
C  is this bend a partner?
C
       IF(j2.EQ.i2) THEN
        If(j1.EQ.i1.AND.IVec(j3).EQ.1) Then
C
C  potential new Z-matrix entry found
C  check for near-planarity and get sign (orientation) of bends
C
         CALL S2Bend(i1,i2,i3,j3,XC,SS)
         If(Abs(SS).LT.thrsh) GO TO 40  ! near-planar, reject
c
         IC = IC+1
         IVec(i1) = 1
         Zindx(IC,1) = i1
         Zindx(IC,2) = i2
         Zindx(IC,3) = i3
         Zindx(IC,4) = j3
         Zindx(IC,5) = 1
C
C  make sure we have a right-handed axis system
C
         If(SS.LT.Zero) Then
          Zindx(IC,3) = j3
          Zindx(IC,4) = i3
         EndIf
         GO TO 20
        Else If(j3.EQ.i1.AND.IVec(j1).EQ.1) Then
C
C  potential new Z-matrix entry found
C  check for near-planarity and get sign (orientation) of bends
C
         CALL S2Bend(i1,i2,i3,j1,XC,SS)
         If(Abs(SS).LT.thrsh) GO TO 40  ! near-planar, reject
c
         IC = IC+1
         IVec(i1) = 1
         Zindx(IC,1) = i1
         Zindx(IC,2) = i2
         Zindx(IC,3) = i3
         Zindx(IC,4) = j1
         Zindx(IC,5) = 1
C
C  make sure we have a right-handed axis system
C
         If(SS.LT.Zero) Then
          Zindx(IC,3) = j1
          Zindx(IC,4) = i3
         EndIf
         GO TO 20
        EndIf
       ENDIF
 40    CONTINUE
      ENDIF
 50   CONTINUE
C
C  If we get here then we didn't find a bend-pair
C  Look for a torsion instead
C
 94   CONTINUE
      DO 60 I=NT1,NT2
c
      i1 = klist(1,I)
      i2 = klist(2,I)
      i3 = klist(3,I)
      i4 = klist(4,I)
C
C  in order to be acceptable as the next line in the Z matrix
C  either i1 or i4 must be a new atom, and the other three
C  indices must already be present in the Z matrix
C
      IF(IVec(i1).EQ.0) THEN
       If(IVec(i2)+IVec(i3)+IVec(i4).EQ.3) Then
C
C  new Z-matrix entry found
C
        IC = IC+1
        IVec(i1) = 1
        Zindx(IC,1) = i1
        Zindx(IC,2) = i2
        Zindx(IC,3) = i3
        Zindx(IC,4) = i4
        GO TO 20
       EndIf
      ELSE IF(IVec(i4).EQ.0) THEN
       If(IVec(i1)+IVec(i2)+IVec(i3).EQ.3) Then
C
C  new Z-matrix entry found
C
        IC = IC+1
        IVec(i4) = 1
        Zindx(IC,1) = i4
        Zindx(IC,2) = i3
        Zindx(IC,3) = i2
        Zindx(IC,4) = i1
        GO TO 20
       EndIf
      ENDIF
 60   CONTINUE
C
C  If we get here it means that we haven't found a suitable
C  torsion for the next Z-matrix entry
      WRITE(IOut,1100)
      RETURN
C
C ................................................................
C
 95   CONTINUE
C
C  Successful Z-matrix construction
C  Check that all stretches and simple bends for torsions exist
C  and load Z-matrix values into ZVals
C
      IC = 0
C
C  load bonds
C
      DO 70 I=2,NAtoms
      i1 = Zindx(I,1)
      i2 = Zindx(I,2)
      DO 69 J=1,NB1-1
      IF( (i1.EQ.klist(1,J).AND.i2.EQ.klist(2,J)) .OR.
     $    (i1.EQ.klist(2,J).AND.i2.EQ.klist(1,J)) ) THEN
       IC = IC+1
       ZVals(I,1) = XPRIM(J)
       INDX(IC) = J
       GO TO 70
      ENDIF
 69   CONTINUE
C
C  if we get here then we didn't find the bond
C
      If(IPRNT.GT.1) WRITE(IOut,1200) i1,i2
      RETURN
 70   CONTINUE
C
C  load angles
C
      DO 80 I=3,NAtoms
      i1 = Zindx(I,1)
      i2 = Zindx(I,2)
      i3 = Zindx(I,3)
      DO 79 J=NB1,NB2
      IF(  i2.EQ.klist(3,J) .AND.
     $   ((i1.EQ.klist(1,J).AND.i3.EQ.klist(2,J)) .OR.
     $    (i1.EQ.klist(2,J).AND.i3.EQ.klist(1,J))) ) THEN
       IC = IC+1
       ZVals(I,2) = XPRIM(J)
       INDX(IC) = J
       GO TO 80
      ENDIF
 79   CONTINUE
C
C  if we get here then we didn't find the angle
C
      If(IPRNT.GT.1) WRITE(IOut,1300) i1,i2,i3
      RETURN
 80   CONTINUE
C
C  load torsions, angle pairs and out-of-plane bends
C
      DO 90 I=4,NAtoms
      i1 = Zindx(I,1)
      i2 = Zindx(I,2)
      i3 = Zindx(I,3)
      i4 = Zindx(I,4)
      ii = Zindx(I,5)
c
      IF(ii.EQ.0) THEN
C
C  we have a torsion
C
       DO 84 J=NT1,NT2
       IF( (i1.EQ.klist(1,J).AND.i2.EQ.klist(2,J).AND.
     $      i3.EQ.klist(3,J).AND.i4.EQ.klist(4,J)) .OR.
     $     (i1.EQ.klist(4,J).AND.i2.EQ.klist(3,J).AND.
     $      i3.EQ.klist(2,J).AND.i4.EQ.klist(1,J)) ) THEN
        IC = IC+1
        ZVals(I,3) = XPRIM(J)
        INDX(IC) = J
        GO TO 90
       ENDIF
 84    CONTINUE
cc
      ELSE IF(ii.EQ.1) THEN
C
C  we have a second bend
C
       DO 85 J=NB1,NB2
       IF(  i2.EQ.klist(3,J) .AND.
     $    ((i1.EQ.klist(1,J).AND.i4.EQ.klist(2,J)) .OR.
     $     (i1.EQ.klist(2,J).AND.i4.EQ.klist(1,J))) ) THEN
        IC = IC+1
        ZVals(I,3) = XPRIM(J)
        INDX(IC) = J
        GO TO 90
       ENDIF
 85    CONTINUE
cc
      ELSE
C
C  we have an out-of-plane bend
C
       DO 86 J=NP1,NP2
       IF( (i1.EQ.klist(1,J).AND.i2.EQ.klist(4,J)) .AND.
     $     (i3.EQ.klist(2,J).AND.i4.EQ.klist(3,J)) .OR.
     $     (i3.EQ.klist(3,J).AND.i4.EQ.klist(2,J)) ) THEN
        IC = IC+1
        ZVals(I,3) = XPRIM(J)
        INDX(IC) = J
        GO TO 90
       ENDIF
 86    CONTINUE
cc
      ENDIF
cc
 90   CONTINUE
C
C ................................................................
C
C  The atoms now need to be renumbered to a canonical
C  z-matrix form. The mapping from original to canonical
C  numbering is given in zmap.
C
C
c -- create the mapping
      IOrd = 0
      DO i=1,NAtoms
      Zmap(i,1) = Zindx(i,1)
      Zmap(i,2) = i
      If ( Zmap(i,1) .NE. Zmap(i,2) ) IOrd = 1
      END DO
c
c  now use zmap to put zindx in canonical order
c
      IF ( IOrd.EQ.1 ) THEN
        DO j=1,4
        DO i=1,NAtoms
        ii = Zindx(i,j)
        DO k=1,NAtoms
        If ( ii.EQ.Zmap(k,1) ) Zindx(i,j) = Zmap(k,2)
        END DO
        END DO
        END DO
      END IF
C
C  Zindx is now in canonical order so we can shift the columns
C  to the left so that the first three columns of Zindx contain
C  the Z matrix to Cartesian input
C
      DO j=1,4
      DO i=1,NAtoms
        Zindx(i,j) = Zindx(i,j+1)
      END DO
      END DO
C
      RETURN
c
 1100 FORMAT(' Unable to find suitable torsion in <Tor2ZN>')
 1200 FORMAT(' Stretch ',2I4,' not found in primitive list')
 1300 FORMAT(' Bend ',3I4,' not found in primitive list')
c
      END
c =======================================================================
c
      SUBROUTINE TranBM(NAtoms,intcor,NDEG,UT,B,VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transform B-Matrix according to vectors in UT
C  On exit, VM contains the transformed B-Matrix
C
      REAL*8 UT(intcor,NDEG),B(3*NAtoms,intcor),VM(3*NAtoms,NDEG)
C
      NAT3 = 3*NAtoms
      CALL ZeroIT(VM,NAT3*NDEG)
c
      DO 20 J=1,NDEG
      DO 20 K=1,intcor
      UVal = UT(K,J)
      DO 10 I=1,NAT3
      VM(I,J) = VM(I,J) + UVal*B(I,K)
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE TranINT(intcor,NDEG,ZINT,UT,XINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transforms primitive internal coordinates to compound
C  natural internal set
C
      REAL*8 ZINT(intcor),UT(intcor,NDEG),XINT(NDEG)
C
      DO 10 I=1,NDEG
      XINT(I) = SProd(intcor,UT(1,I),ZINT)
 10   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE TrnHPRIM(NDEG,NPrim,HPRIM,UT,VV,HINT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Transforms primitive Hessian to Hessian over latest set
C  of non-redundant coordinates
C
C  ARGUMENTS
C
C  NDEG    -  size of non-redundant subspace
C  NPrim   -  total number of primitive internal coordinates
C  HPRIM   -  Hessian in primitive internals
C  UT      -  current set of non-redundant internal coordinates
C             (transformation matrix)
C  VV      -  scratch matrix
C  HINT    -  on exit contains transformed Hessian
C
C
      REAL*8 HPRIM(NPrim,NPrim),UT(NPrim,NDEG),VV(NDEG,NPrim),
     $       HINT(NDEG,NDEG)
C
      PARAMETER (Zero=0.0d0)
C
C
C  transform columns
C
      DO 20 I=1,NDEG
      DO 20 J=1,NPrim
      VAL = Zero
      DO 10 K=1,NPrim
      VAL = VAL + UT(K,I)*HPRIM(K,J)
 10   CONTINUE
      VV(I,J) = VAL
 20   CONTINUE
C
C  transform rows
C
      DO 40 I=1,NDEG
      DO 40 J=1,NDEG
      VAL = Zero
      DO 30 K=1,NPrim
      VAL = VAL + VV(I,K)*UT(K,J)
 30   CONTINUE
      HINT(I,J) = VAL
 40   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE UpdHES(N,      IUpDat, IPRNT,  D,      G,
     $                  GOld,   V1,     V2,     HESS )
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine is responsible for updating the Hessian
C  during a geometry optimization
C
C  Several updating procedures are available:
C   (a) The Murtagh-Sargent update
C       A simple rank-one update formula which allows
C       the Hessian eigenvalue structure to change
C       (see Comp.J. 13 (1970) 185)
C   (b) The Powell update
C       This is a flexible, general rank-two update which allows
C       the Hessian eigenvalue structure to change.
C       (see Math.Prog. 1 (1971) 26)
C   (c) A Composite Powell/Murtagh-Sargent update
C       Suggested by Bofill as an improved update in practice
C       (see J.Comp.Chem. 15 (1994) 1)
C       Default update for transition state search
C   (d) The BFGS update
C       This update is more likely to retain positive definiteness.
C       Default update for minimization
C
C
C  ARGUMENTS
C
C  N       -  dimension of Hessian
C  IUpDat  -  update control parameter
C              0 - skip update this cycle
C              1 - Murtagh-Sargent update
C              2 - Powell update
C              3 - Powell/Murtagh-Sargent update
C              4 - BFGS update
C              5 - BFGS update with positive definite check
C                  (skip update if threatened)
C  IPRNT   -  controls level of print out
C  D       -  geometry displacement on previous cycle
C  G       -  current gradient
C  GOld    -  previous gradient
C  V1      -  scratch vector (dimension N)
C  V2      -  scratch vector (dimension N)
C  HESS    -  Hessian        old Hessian on input
C                            new Hessian on exit
C
C  ...............................................................
C  ** NOTE: The update will be skipped if any calculated scalar **
C  **       product in the update formula is less than TollZero **
C  ...............................................................
C
C
      REAL*8 D(N),G(N),GOld(N),V1(N),V2(N),HESS(N,N)
C
      PARAMETER (One=1.0d0,TollZero=1.0d-8)
C
C
      IOut = ioutfil('iout')
C
C  Skip update if requested
C
      IF(IUpDat.EQ.0) THEN
       If(IPRNT.GT.1) WRITE(IOut,1000)
       RETURN
      ENDIF
C
C  Prepare for update
C  Form and save Hess*D in V1
C
      CALL MatVEC(N,HESS,D,V1)
c
      IF(IUpDat.EQ.1) THEN
cc
C  (a) Murtagh-Sargent update
C
       If(IPRNT.GT.1) WRITE(IOut,1100)
c
       DO 10 I=1,N
       V2(I) = G(I) - GOld(I) - V1(I)
 10    CONTINUE
c
       DT = SProd(N,V2,D)
c
       If(Abs(DT).LT.TollZero) Then
        If(IPRNT.GT.1) WRITE(IOut,1200)
        RETURN
       EndIf
c
       DO 20 I=1,N
       DO 20 J=1,I
       HESS(I,J) = HESS(I,J) + V2(I)*V2(J)/DT
       HESS(J,I) = HESS(I,J)
 20    CONTINUE
cc
      ELSE IF(IUpDat.EQ.2) THEN
cc
C  (b) Powell update
C
       If(IPRNT.GT.1) WRITE(IOut,1300)
c
       DO 30 I=1,N
       V2(I) = G(I) - GOld(I) - V1(I)
 30    CONTINUE
c
       DD = SProd(N,D,D)
c
       If(DD.LT.TollZero) Then
        If(IPRNT.GT.1) WRITE(IOut,1200)
        RETURN
       EndIf
c
       DT = SProd(N,V2,D)
       DTDD = DT/DD
c
       DO 40 I=1,N
       DO 40 J=1,I
       TmP = V2(I)*D(J) + D(I)*V2(J) - D(I)*DTDD*D(J)
       HESS(I,J) = HESS(I,J) + TmP/DD
       HESS(J,I) = HESS(I,J)
 40    CONTINUE
cc
      ELSE IF(IUpDat.EQ.3) THEN
cc
C  (c) Powell/Murtagh-Sargent update
C
       If(IPRNT.GT.1) WRITE(IOut,1400)
C
       DO 50 I=1,N
       V2(I) = G(I) - GOld(I) - V1(I)
 50    CONTINUE
c
       DD = SProd(N,D,D)
       DT = SProd(N,V2,D)
c
       If(DD.LT.TollZero.OR.Abs(DT).LT.TollZero) Then
        If(IPRNT.GT.1) WRITE(IOut,1200)
        RETURN
       EndIf
c
       TT = SProd(N,V2,V2)
       DTDD = DT/DD
C
C  find mixing factor
C
       Zp = ONE - (DT*DT)/(DD*TT)
       Zm = ONE - Zp
c
       If(IPRNT.GT.1) WRITE(IOut,1500) Zp,Zm
c
       DO 60 I=1,N
       DO 60 J=1,I
       TmP = V2(I)*D(J) + D(I)*V2(J) - D(I)*DTDD*D(J)
       TmM = V2(I)*V2(J)
       HESS(I,J) = HESS(I,J) + Zp*TmP/DD + Zm*TmM/DT
       HESS(J,I) = HESS(I,J)
 60    CONTINUE
cc
      ELSE
cc
C  (d) BFGS update
C
       If(IPRNT.GT.1) WRITE(IOut,1600)
C
       DO 70 I=1,N
       V2(I) = G(I) - GOld(I)
 70    CONTINUE
c
       DD = SProd(N,V2,D)
C
C  If DD is negative, retention of positive definiteness is not
C  guaranteed. Print a warning and skip update if requested
C
       IF(DD.LT.TollZero) THEN
        If(IPRNT.GT.1) WRITE(IOut,1700)
        IF(IUpDat.EQ.5) THEN
         If(IPRNT.GT.1) WRITE(IOut,1800)
         RETURN
        ENDIF
       ENDIF
c
       DT = SProd(N,V1,D)
c
       If(Abs(DT).LT.TollZero.OR.Abs(DD).LT.TollZero) Then
        If(IPRNT.GT.1) WRITE(IOut,1200)
        RETURN
       EndIf
c
       DO 80 I=1,N
       DO 80 J=1,I
       TmP = (V2(I)*V2(J))/DD - (V1(I)*V1(J))/DT
       HESS(I,J) = HESS(I,J) + TmP
       HESS(J,I) = HESS(I,J)
 80    CONTINUE
cc
      ENDIF
C
      RETURN
c
 1000 FORMAT(' Hessian Update Skipped by Request')
 1100 FORMAT(' Hessian Updated using Murtagh-Sargent Update')
 1200 FORMAT('**WARNING** Small Scalar Product',/,
     $       ' Hessian Update Skipped this cycle')
 1300 FORMAT(' Hessian Updated using Powell Update')
 1400 FORMAT(' Hessian Updated using Powell/Murtagh-Sargent Update')
 1500 FORMAT('  Mixing factors: ',F12.6,' Powell',/
     $       18X,F12.6,' Murtagh-Sargent')
 1600 FORMAT(' Hessian Updated using BFGS Update')
 1700 FORMAT('**WARNING** Hereditary positive definiteness',
     $       ' endangered')
 1800 FORMAT(' Hessian Update Skipped this cycle')
c
      END
c =======================================================================
c
      SUBROUTINE UTPrim(NDEG,   NPrim,  NPrim1, GOld,   DINT,
     $                  GINT,   UT,     UTOld,  GPOld,  DP,
     $                  GP)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Tranforms old gradient, step and new gradient over non-redundant
C  coordinate space to their equivalents over primitive internals
C  prior to updating the primitive Hessian
C
C  ARGUMENTS
C
C  NDEG    -  size of non-redundant subspace
C  NPrim   -  number of internals in current primitive space
C  NPrim1  -  number of internals in previous primitive space
C  GOld    -  old gradient over non-redundant coordinates
C  DINT    -  displacement on previous cycle
C  GINT    -  current gradient over non-redundant coordinates
C  UT      -  current set of non-redundant coordinates
C  UTOld   -  previous set of non-redundant coordinates
C
C  on exit
C
C  GPOld   -  old gradient in primitive coordinates
C  DP      -  displacement in primitive coordinates
C  GP      -  gradient in primitive coordinates
C
C
      REAL*8 GOld(NDEG),DINT(NDEG),GINT(NDEG)
      REAL*8 UT(NPrim,NDEG),UTOld(NPrim1,NDEG)
      REAL*8 GPOld(NPrim1),DP(NPrim1),GP(NPrim)
C
C
C  transform quantities
C
      CALL ZeroIT(GPOld,NPrim1)
      CALL ZeroIT(DP,NPrim1)
      CALL ZeroIT(GP,NPrim)
c
      DO 30 J=1,NDEG
      Val1 = GOld(J)
      Val2 = DINT(J)
      Val3 = GINT(J)
      DO 10 I=1,NPrim1
      GPOld(I) = GPOld(I) + UTOld(I,J)*Val1
      DP(I) = DP(I) + UTOld(I,J)*Val2
 10   CONTINUE
      DO 20 I=1,NPrim
      GP(I) = GP(I) + UT(I,J)*Val3
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE ZeroFIX(NAT3,IFix,GC,D,HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Zeros all gradient, displacement and Hessian components
C  corresponding to fixed coordinates
C
      DIMENSION IFix(NAT3),GC(NAT3),D(NAT3),HESS(NAT3,NAT3)
C
      PARAMETER (Zero=0.0d0)
C
      DO 20 I=1,NAT3
      IF(IFix(I).EQ.1) THEN
       GC(I) = Zero
       D(I) = Zero
       DO 10 J=1,NAT3
       HESS(I,J) = Zero
       HESS(J,I) = Zero
 10    CONTINUE
      ENDIF
 20   CONTINUE
C
      RETURN
      END
