c ==================================================================
c  Z-MATRIX READ/WRITE ROUTINES        JB   Aug 1997
c ==================================================================
c
      SUBROUTINE RdZMAT(NZ,     ZSymb,  ZMTINP, CINTNL, FINTNL,
     $                  GEO,    IGEO,   IG,     VARNAM, NVar)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C  Reads Z-Matrix from control file
C  ** Written initially by Fora Chan **
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres (including dummy atoms)
C             in Z-matrix
C  ZSymb   -  Z-matrix symbols (real & dummy atoms)
C  ZMTINP  -  intermediate scratch storage for the up to 7 possible
C             separate data types per Z-matrix line
C  CINTNL  -  character scratch array for input parameters
C  FINTNL  -  real scratch array for values of input parameters
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO    -  Z-matrix connectivity
C              Atoms are usually defined with respect to previously
C              defined atoms by a stretch, a bend and a torsion;
C              an alternative is a stretch and 2 bends, or a stretch,
C              a bend and an out-of-plane bend. The fourth column
C              of IGEO is used to distinguish these cases
C                0 - bend + torsion
C                1 - 2 bends
C                2 - bend + out-of-plane bend
C      ** WARNING - CURRENTLY SET UP FOR TORSIONS ONLY **
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  VARNAM  -  names of all "variables" (including fixed)
C  NVar    -  number of variables
C
C
      DIMENSION GEO(NZ,3),IGEO(NZ,4),IG(3*NZ),FINTNL(3*NZ)
      CHARACTER*8 ZSymb(NZ),ZMTINP(NZ,7),CINTNL(3*NZ,2),VARNAM(3*NZ)
      CHARACTER BLANK,TAB,BLANK8*8
      CHARACTER CHR*1,ZLINE*80
      character*256 jobname
C
      PARAMETER (BLANK = ' ',TAB = '	',BLANK8 = '        ')
      PARAMETER (MXCHAR=80)
      COMMON /COLLZERO/ JZERO
c
      Common /job/jobname,lenJ
C
C
C  open ZMAT file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.zmat',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  first locate $zmat section
C
 10   CONTINUE
      READ(40,900,END=96) ZLINE
      If(ZLINE(1:5).NE.'$zmat') GO TO 10
C
C  Initialize position of zero in Fortran collating series.
C
      JZERO = ICHAR('0')
C
C  Initialize arrays
C
      DO 20 J = 1,7
      DO 20 I = 1,NZ
      ZMTINP(I,J) = BLANK8
 20   CONTINUE
      CALL IZeroIT(IG,3*NZ)
      CALL IZeroIT(IGEO,4*NZ)
C
C  Read Line of Z-Matrix into ZLINE
C  NZ is number of lines in Z-Matrix
C  NINTNL is number of lines of internal variables.
C
      NINTNL = 0
c
      DO 30 I = 1,NZ
      Read (40,900,END=97) ZLINE
C
C  Assign main body of Z-matrix to ZMTINP
C
      IFIRST = 1
      CALL RmBLNK1(IFIRST,ZLINE,BLANK,BLANK,MXCHAR)
      J = 1
      INEXT = IFIRST + 1
      DO WHILE (INEXT .LT. MXCHAR)
      IF ( (ZLINE(INEXT:INEXT) .EQ. BLANK ) .OR.
     +     (ZLINE(INEXT:INEXT) .EQ. TAB) ) THEN
        ZMTINP(I,J) = ZLINE(IFIRST:INEXT-1)
        J = J + 1
        If (J .GT. 7) GO TO 30
        IFIRST = INEXT + 1
        CALL RmBLNK1(IFIRST,ZLINE,BLANK,BLANK,MXCHAR)
        INEXT = IFIRST + 1
      ELSE
        INEXT = INEXT + 1
      ENDIF
      ENDDO
 30   CONTINUE
C
C  Read and assign variables to CINTRL
C
      Read (40,900) ZLINE
  40  Read (40,900,END=60) ZLINE
      If(ZLINE(1:4).EQ.'$end') GO TO 60
c
      NINTNL = NINTNL+1
      IFIRST = 1
      CALL RmBLNK1(IFIRST,ZLINE,BLANK,BLANK,MXCHAR)
      J = 1
      INEXT = IFIRST + 1
      DO WHILE (INEXT .LT. MXCHAR)
      CHR = ZLINE(INEXT:INEXT)
      IF ((CHR .EQ. BLANK) .OR. (CHR .EQ. TAB) .OR.
     +    (CHR .EQ. '=')) THEN
        CINTNL(NINTNL,J) = ZLINE(IFIRST:INEXT-1)
        J = J + 1
        If (J .GT. 2) GO TO 50
        IFIRST = INEXT + 1
        CALL RmBLNK1(IFIRST,ZLINE,BLANK,'=',MXCHAR)
        INEXT = IFIRST + 1
      ELSE
        INEXT = INEXT + 1
      ENDIF
      ENDDO
 50   CONTINUE
c
      GO TO 40
 60   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  Check number of parameters
C
      If(NINTNL.GT.3*NZ) Then
       Call nerror(1,'Z-Matrix Read',
     $  'Too Many Parameters in Z-Matrix Parameter List',0,0)
      EndIf
C
C  Convert CINTNL(*,2) to FINTNL which is an array of real numbers.
C
      DO 70 I = 1,NINTNL
      CALL CHRTOF(CINTNL(I,2),FLT)
      FINTNL(I) = FLT
 70   CONTINUE
C
C  Assign ZMTINP(*,1) to ZSymb.
C
      DO 80 I = 1,NZ
      ZSymb(I) = ZMTINP(I,1)
 80   CONTINUE
c
      NIC = 0
C
C  Assign ZMTINP(*,2) to IGEO(*,1)
C
      CALL TOIGEO(ZMTINP,2,IGEO,1,ZSymb,NZ)
C
C  Assign ZMTINP(*,3) to GEO(*,1).
C
      CALL TOGEO(ZMTINP, 3,      GEO,    1,      NZ,
     +           NINTNL, CINTNL, FINTNL, IG,     VARNAM,
     +           NIC)
C
C  Assign ZMTINP(*,4) to IGEO(*,2)
C
      CALL TOIGEO(ZMTINP,4,IGEO,2,ZSymb,NZ)
C
C  Assign ZMTINP(*,5) to GEO(*,2).
C
      CALL TOGEO(ZMTINP, 5,      GEO,    2,      NZ,
     +           NINTNL, CINTNL, FINTNL, IG,     VARNAM,
     +           NIC)
C
C  Assign ZMTINP(*,6) to IGEO(*,3)
C
      CALL TOIGEO(ZMTINP,6,IGEO,3,ZSymb,NZ)
C
C  Assign ZMTINP(*,7) to GEO(*,3).
C
      CALL TOGEO(ZMTINP, 7,      GEO,    3,      NZ,
     +           NINTNL, CINTNL, FINTNL, IG,     VARNAM,
     +           NIC)
C
C  ..................................................................
C
C  CHECK Z-MATRIX
C  check for bad connectivity and near-tetrahedral angles
C
      CALL ChkZMAT(NZ,IGEO,GEO)
C  ..................................................................
C
C  Determine the actual number of Z-Matrix variables by
C  looking for zeros in the IG array
C
      NVar = 0
      DO 90 I=1,NIC
      IF(IG(I).EQ.0) NVar=NVar+1
 90   CONTINUE
C
      RETURN
C
C  -------------------------------------------------------------------
C  ERROR SECTION
C
 95   CONTINUE          ! no ZMAT file
      Call nerror(2,'Z-Matrix Read',
     $  'Z-Matrix Read Requested but no <zmat> file found!',0,0)
c
 96   CONTINUE          ! no $zmatrix section found
      Call nerror(3,'Z-Matrix Read',
     $      'No $zmatrix section found on <zmat> file!',0,0)
c
 97   CONTINUE          ! unexpected end of Z-matrix
      Call nerror(4,'Z-Matrix Read',
     $      'Unexpected End of Z-Matrix on <zmat> file',0,0)
c
  900 Format(A80)
c
      END
c  =======================================================================
c
      SUBROUTINE CHRTOF(STRNG,FLT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C  Converts a character string to a floating point number
C
C  ARGUMENTS
C
C  STRNG   -  input character string
C  FLT     -  on exit contains floating point number
C
C
      INTEGER ASIZE
      CHARACTER BLANK,TAB
      CHARACTER*8 STRNG,DEC1,DEC2
      COMMON /COLLZERO/ JZERO
      LOGICAL NEGTVE
C
      PARAMETER (BLANK = ' ',TAB = '	',ASIZE = 8)
C
      iout=igetival('iout')
C
C Data checking.
C
      DO 10 I = 1,ASIZE
      ICH = ICHAR(STRNG(I:I)) - JZERO
      IF ( ((ICH .LT. 0) .OR. (ICH .GT. 9)) .AND.
     +	 (STRNG(I:I) .NE. '.') .AND. (STRNG(I:I) .NE. BLANK)
     +    .AND. (STRNG(I:I) .NE. TAB) .AND. (STRNG(I:I) .NE. '-')
     +    .AND. (STRNG(I:I) .NE. '+') ) THEN
        WRITE(iout,1000) STRNG
        CALL OptExit(9)
      ENDIF
  10  CONTINUE
C
C Check negative numbers.
C
      NEGTVE = .FALSE.
      IF (STRNG(1:1) .EQ. '+') THEN
        IFIRST = 2
      ELSE IF (STRNG(1:1) .EQ. '-') THEN
        IFIRST = 2
        NEGTVE = .TRUE.
      ELSE
        IFIRST = 1
      ENDIF
C	
      DO 20 INEXT = 1,ASIZE
      IF (STRNG(INEXT:INEXT) .EQ. '.') THEN
        DEC1 = STRNG(IFIRST:INEXT-1)
        DEC2 = STRNG(INEXT+1:ASIZE)
        GO TO 21
      else if(strng(inext:inext).eq.' ') then
        dec1=strng(ifirst:inext-1)
        dec2='00000000'
        go to 21
      ENDIF
  20  CONTINUE
  21  CONTINUE
C
      CALL CHRTOI(DEC1,IDEC1,ILEN)
      CALL CHRTOI(DEC2,IDEC2,ILEN)
C
      IEXPNT = ILEN
      IFCTOR = 10 ** IEXPNT
      FDEC = DFLOAT(IDEC2) / DFLOAT(IFCTOR)
      FLT = DFLOAT(IDEC1) + FDEC
C
      If (NEGTVE) FLT = -FLT
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Real Number Expected for Z-Matrix',
     $            ' Variable: ',A8)
c
      END
c  =======================================================================
c
      SUBROUTINE CHRTOI(STRNG,INTGR,ILEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C  Converts a character string to an integer
C
C  ARGUMENTS
C
C  STRNG   -  input character string
C  INTGR   -  on exit contains integer
C  ILEN    -  length of integer
C
      INTEGER ASIZE
      CHARACTER BLANK,TAB
      CHARACTER*8 STRNG
      COMMON /COLLZERO/ JZERO
C
      PARAMETER (BLANK = ' ',TAB = '	',ASIZE = 8)
      iout=igetival('iout')
C
C Data checking.
C
      DO 10 I = 1,ASIZE
      ICH = ICHAR(STRNG(I:I)) - JZERO
      IF ( ((ICH .LT. 0) .OR. (ICH .GT. 9)) .AND.
     +     (STRNG(I:I) .NE. BLANK) .AND. (STRNG(I:I) .NE. TAB) ) THEN
        WRITE(iout,1000) STRNG
        CALL OptExit(9)
      ENDIF
  10  CONTINUE
C
      IFIRST = 1
      DO 20 INEXT = 1,ASIZE
      IF ( (STRNG(INEXT:INEXT) .EQ. BLANK)  .OR.
     +     (STRNG(INEXT:INEXT) .EQ. TAB) ) THEN
        ILEN = LEN(STRNG(IFIRST:INEXT-1))
        GO TO 21
      ENDIF
  20  CONTINUE
  21  CONTINUE
C
      INTGR = 0
      DO 30 I = 1,ILEN
      ICH = ICHAR(STRNG(I:I)) - JZERO
      IEXPNT = ILEN - I
      IFCTOR = 10 ** IEXPNT
      INTGR = INTGR + ICH * IFCTOR
  30  CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Integer Expected on Z-Matrix Line:',
     $       A80)
c
      END
c  =======================================================================
c
      SUBROUTINE ChkZMAT(NZ,IGEO,GEO)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks the final input Z-Matrix for bad connectivity,
C  negative bond lengths or angles and converts
C  near-tetrahedral angles to exact value
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  IGEO    -  Z-matrix connectivity
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C
C
      DIMENSION IGEO(NZ,3),GEO(NZ,3)
      LOGICAL Error
C
      PARAMETER (Tetra=109.4712206344907d0,thrsh=1.0d-3)
      DATA Error/.FALSE./
C
      iout=igetival('iout')
C
C  Check connectivity
C
      DO 10 I=2,NZ
C
C  check for not yet defined centres
C
      If(IGEO(I,1).GE.I.OR.IGEO(I,2).GE.I.OR.IGEO(I,3).GE.I) Then
       Error = .TRUE.
       WRITE(iout,1000) I
      EndIf
C
C  check for multiply connected centres on same line
C
      If( (IGEO(I,1).EQ.IGEO(I,2)) .OR.
     $    (IGEO(I,1).EQ.IGEO(I,3)) .OR.
     $    (IGEO(I,2).EQ.IGEO(I,3).AND.I.GT.2) ) Then
       Error = .TRUE.
       WRITE(iout,1100) I
      EndIf
C
C  check for negative bond lengths or angles
C
      If(GEO(I,1).LT.0) Then
       Error = .TRUE.
       WRITE(iout,1200) I
      EndIf
      If(GEO(I,2).LT.0) Then
       Error = .TRUE.
       WRITE(iout,1300) I
      EndIf
C
 10   CONTINUE
c
      If(Error) CALL OptExit(9)
C
C  check for near-tetrahedral angles
C
      DO 20 I=3,NZ
      Ang = GEO(I,2)
      If(Abs(Ang-Tetra).LT.thrsh) GEO(I,2)=Tetra
 20   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Z-Matrix Line ',I3,' refers to a',
     $            ' not yet defined centre')
 1100 FORMAT(/,2X,'***ERROR*** Z-Matrix Line ',I3,' has multiple',
     $            ' references to the same centre')
 1200 FORMAT(/,2X,'***ERROR*** Z-Matrix Line ',I3,' has a',
     $            ' negative bond length')
 1300 FORMAT(/,2X,'***ERROR*** Z-Matrix Line ',I3,' has a',
     $            ' negative bond angle')
c
      END
c  =======================================================================
c
      SUBROUTINE GMETRY(NZ,GEO,IGEO,XZ)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Computes Cartesian coordinates from Z-matrix
C  **  ADAPTED FROM CODE WRITTEN BY M.J.S.DEWAR **
C  **  initial modification:  Fora Chan         **
C
C ............................................................
C  Originally this routine used bond lengths, bond angles
C  and torsions only in the Z matrix. Subsequently extended
C  to include pairs of angles and out-of-plane bends
C ............................................................
C
C  ARGUMENTS
C
C  NZ      -  number of centres (including dummy atoms)
C  GEO     -  internal coordinate values
C  IGEO    -  Z-matrix connectivity data
C  XZ      -  on exit contains Cartesian coordinates
C
C  NOTE:   1st  centre at origin
C          2nd  centre along Z-Axis
C          3rd  centre in XZ plane
C
C
      DIMENSION GEO(NZ,3),IGEO(NZ,4),XZ(3,NZ)
      Dimension U(3),V(3),W(3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0)
      PARAMETER (TolZero=1.0d-16,thrsh1=0.1d0,thrsh4=1.0d-4)
      Logical jubel
C
C
      iout = ioutfil('iout')
C
C  initialize
C
      XZ(1,1) = Zero
      XZ(2,1) = Zero
      XZ(3,1) = Zero
      XZ(1,2) = Zero
      XZ(2,2) = Zero
      XZ(3,2) = GEO(2,1)
c
      If(NZ.EQ.2) GO TO 20
c
      CCOS = COS(GEO(3,2))
      IF(IGEO(3,1).EQ.1)THEN
       XZ(3,3) = XZ(3,1) + GEO(3,1)*CCOS
      ELSE
       XZ(3,3) = XZ(3,2) - GEO(3,1)*CCOS
      ENDIF
      XZ(1,3) = GEO(3,1)*SIN(GEO(3,2))
      XZ(2,3) = Zero
c
      If(NZ.EQ.3) GO TO 20
c
      DO 10 I=4,NZ
      IF(IGEO(I,4).EQ.0) THEN
C
C  atom defined with a torsion
C
       COSA = COS(GEO(I,2))
       MB = IGEO(I,2)
       MC = IGEO(I,1)
       XB = XZ(1,MB)-XZ(1,MC)
       YB = XZ(2,MB)-XZ(2,MC)
       ZB = XZ(3,MB)-XZ(3,MC)
       RBC = XB*XB+YB*YB+ZB*ZB
c
       IF(RBC.LT.TolZero) THEN
C
C  Two atoms are coincident. A fatal error
C
        WRITE(iout,1000) MB,MC
        CALL OptExit(9)
       ELSE
        RBC = One/SQRT(RBC)
       ENDIF
c
       MA = IGEO(I,3)
       XA = XZ(1,MA)-XZ(1,MC)
       YA = XZ(2,MA)-XZ(2,MC)
       ZA = XZ(3,MA)-XZ(3,MC)
C
C  Rotate about the Z-Axis to make YB=0 and XB>0
C  If XYB is too small, first rotate the Y-Axis by 90 degrees
C
       XYB = SQRT(XB*XB+YB*YB)
       K = -1
       IF(XYB.LT.thrsh1) THEN
        XPA = ZA
        ZA = -XA
        XA = XPA
        XPB = ZB
        ZB = -XB
        XB = XPB
        XYB = SQRT(XB*XB+YB*YB)
        K = +1
       ENDIF
C
C  Rotate about the Y-Axis to make ZB vanish
C
       COSTH = XB/XYB
       SINTH = YB/XYB
       XPA = XA*COSTH+YA*SINTH
       YPA = YA*COSTH-XA*SINTH
       SINPH = ZB*RBC
       COSPH = SQRT(ABS(ONE-SINPH*SINPH))
       XQA = XPA*COSPH+ZA*SINPH
       ZQA = ZA*COSPH-XPA*SINPH
C
C  Rotate about the X-Axis to make ZA=0 and YA>0
C
       YZA = SQRT(YPA**2+ZQA**2)
       IF(YZA.GT.thrsh4) THEN
        COSKH = YPA/YZA
        SINKH = ZQA/YZA
       ELSE
C
C  Angle too small to be important
C
        COSKH = One
        SINKH = Zero
       ENDIF
C
C  Coordinates are:  A = (XQA,YZA,0)  B = (RBC,0,0)  C = (0,0,0)
C  None are negative.
C  Coordinates of the current centre in the new frame are:
C
       SINA = SIN(GEO(I,2))
       SIND = -SIN(GEO(I,3))
       COSD = COS(GEO(I,3))
       XD = GEO(I,1)*COSA
       YD = GEO(I,1)*SINA*COSD
       ZD = GEO(I,1)*SINA*SIND
C
C  Now transform back to the original frame
C
       YPD = YD*COSKH-ZD*SINKH
       ZPD = ZD*COSKH+YD*SINKH
       XPD = XD*COSPH-ZPD*SINPH
       ZQD = ZPD*COSPH+XD*SINPH
       XQD = XPD*COSTH-YPD*SINTH
       YQD = YPD*COSTH+XPD*SINTH
c
       IF(K.GT.0) THEN
        XRD = -ZQD
        ZQD = XQD
        XQD = XRD
       ENDIF
c
       XZ(1,I) = XQD+XZ(1,MC)
       XZ(2,I) = YQD+XZ(2,MC)
       XZ(3,I) = ZQD+XZ(3,MC)
cc
      ELSE IF(IGEO(I,4).EQ.1) THEN
C
C  atom defined using 2 bends
C     I J K L:   bends  I-J-K  and I-J-L
C
       J = IGEO(I,1)
       K = IGEO(I,2)
       L = IGEO(I,3)
C
C  define unit vector U along bond J-K
C                     V along bond J-L
C
       U(1) = XZ(1,K) - XZ(1,J)
       U(2) = XZ(2,K) - XZ(2,J)
       U(3) = XZ(3,K) - XZ(3,J)
       V(1) = XZ(1,L) - XZ(1,J)
       V(2) = XZ(2,L) - XZ(2,J)
       V(3) = XZ(3,L) - XZ(3,J)
       CALL NOM(U,3,jubel,thrsh4)
       CALL NOM(V,3,jubel,thrsh4)
       CALL NORMAL(U,V,W)
c
       C = SProd(3,U,V)
       Y = One-C*C
       COSA = COS(GEO(I,2))
       COSB = COS(GEO(I,3))
c
       X = (COSA - COSB*C)/Y
       Y = (COSB - COSA*C)/Y
       Z = One - X**2 - Y**2 - 2.0d0*X*Y*C
       If(Z.LT.Zero) Then
        If(-Z.LT.thrsh4) Then
         Z = Zero
        Else
         WRITE(iout,1100) J
         z = -sqrt(-z)
        EndIf
       Else
        Z = SQRT(Z)
       EndIf
C
C  form direction of vector along bond J-I
C
       W(1) = X*U(1) + Y*V(1) + Z*W(1)
       W(2) = X*U(2) + Y*V(2) + Z*W(2)
       W(3) = X*U(3) + Y*V(3) + Z*W(3)
C
C  get new Cartesian coordinates
C
       BL = GEO(I,1)
       XZ(1,I) = XZ(1,J) + BL*W(1)
       XZ(2,I) = XZ(2,J) + BL*W(2)
       XZ(3,I) = XZ(3,J) + BL*W(3)
cc
      ELSE
C
C  atom defined using out-of-plane bend
C
       J = IGEO(I,2)
       K = IGEO(I,3)
       L = IGEO(I,1)    ! central atom
C
C  define unit vector U along bond L-J
C                     V along bond L-K
C
       U(1) = XZ(1,J) - XZ(1,L)
       U(2) = XZ(2,J) - XZ(2,L)
       U(3) = XZ(3,J) - XZ(3,L)
       V(1) = XZ(1,K) - XZ(1,L)
       V(2) = XZ(2,K) - XZ(2,L)
       V(3) = XZ(3,K) - XZ(3,L)
       CALL NOM(U,3,jubel,thrsh4)
       CALL NOM(V,3,jubel,thrsh4)
       CALL NORMAL(U,V,W)
       CALL NORMAL(U,W,V)
c
       COSP = COS(GEO(I,3))
       SINP = SIN(GEO(I,3))
       COSB = COS(GEO(I,2))/COSP
       SINB = SQRT(One - COSB*COSB)
c
       U(1) = COSB*U(1) + SINB*V(1)
       U(2) = COSB*U(2) + SINB*V(2)
       U(3) = COSB*U(3) + SINB*V(3)
C
C  form direction of vector along bond L-I
C
       W(1) = COSP*U(1) + SINP*W(1)
       W(2) = COSP*U(2) + SINP*W(2)
       W(3) = COSP*U(3) + SINP*W(3)
C
C  get new Cartesian coordinates
C
       BL = GEO(I,1)
       XZ(1,I) = XZ(1,L) + BL*W(1)
       XZ(2,I) = XZ(2,L) + BL*W(2)
       XZ(3,I) = XZ(3,L) + BL*W(3)
cc
      ENDIF
 10   CONTINUE
 20   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Atoms ',I6,' and ',I6,' Are Coincident!')
 1100 FORMAT(/,2X,'***ERROR*** Bad Bend Pair at atom ',I6)
c
      END
c  =======================================================================
c
      SUBROUTINE PrntZMAT(IOut,NZ,ZSymb,GEO,IGEO)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GEO(NZ,3),IGEO(NZ,3)
      CHARACTER*8 ZSymb(NZ)
C
C  Prints out the Z-Matrix (Angstroms and degrees)
C
      WRITE(IOut,1000)
c
      INum=1
      WRITE(IOut,1100) INum,ZSymb(1)
      If(NZ.GT.1) WRITE(IOut,1100) INum+1,ZSymb(2),IGEO(2,1),GEO(2,1)
      If(NZ.GT.2) WRITE(IOut,1100) INum+2,ZSymb(3),IGEO(3,1),GEO(3,1),
     $                                          IGEO(3,2),GEO(3,2)
c
      DO 10 I=4,NZ
      WRITE(IOut,1100) I,ZSymb(I),IGEO(I,1),GEO(I,1),IGEO(I,2),GEO(I,2),
     $                         IGEO(I,3),GEO(I,3)
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,10X,' Z-Matrix  (Angstroms and Degrees)',/)
 1100 FORMAT(2X,I3,2X,A8,1X,I4,1X,F10.6,1X,I4,1X,F10.4,1X,I4,1X,F10.4)
c
      END
c  =======================================================================
c
      SUBROUTINE RmBLNK1(IFIRST,STRNG,RM1,RM2,MXCHAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER BLANK,TAB
      PARAMETER (BLANK = ' ',TAB = '	')
      CHARACTER RM1,RM2,CHR,STRNG*80
C
      DO WHILE (IFIRST .LT. MXCHAR)
       CHR = STRNG(IFIRST:IFIRST)
       IF( (CHR .EQ. RM1) .OR. (CHR .EQ. RM2) .OR.
     $       (CHR .EQ. TAB) ) THEN
        IFIRST = IFIRST + 1
       ELSE
        GO TO 95
       ENDIF
      ENDDO
 95   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE RmDummy(NZ,ZSymb,XZ,NAtoms,AtSymb,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Removes dummy atoms from coordinate list
C  ** NOTE:  This routine assumes that dummy atom symbols
C            start with x (but not xe), q or Du/du
C
C  ARGUMENTS
C
C  NZ      -  total number of centres (including dummy atoms)
C  ZSymb   -  atomic symbols for all centres
C             (Du for dummy atoms)
C  XZ      -  Cartesian coordinates of all centres
C
C  on exit
C
C  NAtoms  -  number of real atoms
C             (on input set to NZ)
C  AtSymb  -  atomic symbols for real atoms only
C  XC      -  Cartesian coordinates of real atoms only
C
C
      REAL*8 XZ(3,NZ),XC(3,NAtoms)
      CHARACTER*8 ZSymb(NZ),AtSymb(NAtoms)
C
C
C  Loop over all centres in the system
C
      IAtm = 0
      DO 10 I=1,NZ
      If( (ZSymb(I)(1:1).EQ.'x'.AND.ZSymb(I)(2:2).NE.'e') .OR.
     $    ((ZSymb(I)(1:1).EQ.'D'.OR.ZSymb(I)(1:1).EQ.'d') .AND.
     $     ZSymb(I)(2:2).EQ.'u') .OR. ZSymb(I)(1:1).EQ.'q') GO TO 10
      IAtm = IAtm+1
      AtSymb(IAtm) = ZSymb(I)
      XC(1,IAtm) = XZ(1,I)
      XC(2,IAtm) = XZ(2,I)
      XC(3,IAtm) = XZ(3,I)
 10   CONTINUE
C
      NAtoms = IAtm
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE ScanZMAT(NZ,IErr)
      IMPLICIT INTEGER(A-Z)
C
C
C  Scans Z-Matrix to determine number of atomic centres
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres (including dummy atoms)
C             in Z-matrix
C  IErr    -  Error flag   0 - success
C                         -1 - no ZMAT file found or error
C
C
      CHARACTER*8 CHAR,BLANK8
      character*256 jobname
      PARAMETER (BLANK8 = '        ')
c
      Common /job/jobname,lenJ
C
C
      NZ = 0
      IErr = -1
C
C  open <zmat> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.zmat',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  first locate the $zmatrix section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      If(CHAR(1:5).NE.'$zmat') GO TO 10
C
C  now read until blank line or $end encountered
C  this is the number of Z-Matrix centres
C
 20   READ(40,900,END=97) CHAR
      If(CHAR.EQ.BLANK8.OR.CHAR(1:4).EQ.'$end') GO TO 94
      NZ = NZ+1
      GO TO 20
c
 94   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
      IErr = 0
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE          ! no ZMAT file
      RETURN
c
 96   CONTINUE          ! no $zmatrix section found
      Call nerror(3,'Z-Matrix Read',
     $      'No $zmatrix section found on <zmat> file!',0,0)
c
 97   CONTINUE          ! no $end section marker found
      Call nerror(5,'Z-Matrix Read',
     $      'No End-of-File marker ($end) found on <zmat> file!',0,0)
c
  900 Format(A8)
c
      END
c  =======================================================================
c
      SUBROUTINE TOGEO(ZMTINP, IZDIM,  GEO,    IGDIM,  NZ,
     $                 NINTNL, CINTNL, FINTNL, IG,     VARNAM,
     $                 NVar)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C  Assigns parameter entries in ZMTINP array to values in GEO
C
C  ARGUMENTS
C
C  ZMTINP  -  intermediate scratch storage for the up to 7 possible
C             separate data types per Z-matrix line
C  IZDIM   -  location of entry in ZMTINP array to be assigned
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGDIM   -  location in GEO array where assigned ZMTINP entry put
C  NZ      -  total number of centres (lines) in Z-Matrix
C  NINTNL  -  number of parameters in Z-Matrix parameter section
C  CINTNL  -  character scratch array for input parameters
C  FINTNL  -  real scratch array for values of input parameters
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  VARNAM  -  names of all "variables" (including fixed)
C  NVar    -  number of variables
C
C
      DIMENSION GEO(NZ,3),IG(3*NZ),FINTNL(3*NZ)
      CHARACTER*8 ZMTINP(NZ,7),VARNAM(3*NZ),CINTNL(3*NZ,2)
      CHARACTER BLANK,TAB,FIX*8,NEGCIN*9
C
      PARAMETER (BLANK = ' ',TAB = '	',FIX = '        ')
      PARAMETER (IFIX = 1000)
      PARAMETER (ZERO=0.0d0)
C
C
      DO 20 I=1,NZ
C
C  Process blank ZMTINP entries
C
      IF(ZMTINP(I,IZDIM)(1:1).EQ.BLANK) THEN
       GEO(I,IGDIM) = ZERO
       GO TO 20
      ENDIF
c
      NVar = NVar + 1
C
C  Process variables for optimization
C
      DO 10 J=1,NINTNL
      NEGCIN = '-'//CINTNL(J,1)
      IF(ZMTINP(I,IZDIM).EQ.CINTNL(J,1)) THEN
       GEO(I,IGDIM) = FINTNL(J)
       VARNAM(NVar) = CINTNL(J,1)
       CALL TOIG(NZ,ZMTINP,IZDIM,I,NVar,IG,VARNAM)
       GO TO 20
      ELSE IF(ZMTINP(I,IZDIM).EQ.NEGCIN) THEN
       GEO(I,IGDIM) = -FINTNL(J)
       VARNAM(NVar) = CINTNL(J,1)
       CALL TOIG(NZ,ZMTINP,IZDIM,I,NVar,IG,VARNAM)
       GO TO 20
      ENDIF
 10   CONTINUE
C
C  Process fixed variables.
C
      CALL CHRTOF(ZMTINP(I,IZDIM),FLT)
      GEO(I,IGDIM) = FLT
      VARNAM(NVar) = FIX
      IG(NVar) = IFIX
C
 20   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE TOIG(NZ,     ZMTINP, IZDIM,  ICUR,   NVar,
     $                IG,     VARNAM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Assigns optimization options to each Z-Matrix parameter
C
C  ARGUMENTS
C
C  NZ      -  total number of centres (lines) in Z-Matrix
C  ZMTINP  -  intermediate scratch storage for the up to 7 possible
C             separate data types per Z-matrix line
C  IZDIM   -  location of entry in ZMTINP array to be assigned
C  ICUR    -  Z-Matrix parameter being examined
C  NVar    -  variables already assigned
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  VARNAM  -  names of all "variables" (including fixed)
C
C
      DIMENSION IG(3*NZ)
      CHARACTER*8 ZMTINP(NZ,7),VARNAM(3*NZ)
      CHARACTER NEGVNM*9
C
      PARAMETER (IOPT=0)
C
C
C  Search for the same internal variables.
C
      DO 10 I=1,NVar-1
      NEGVNM = '-'//VARNAM(I)
      IF(ZMTINP(ICUR,IZDIM).EQ.VARNAM(I)) THEN
       IG(NVar) = I
       GO TO 20
      ELSE IF(ZMTINP(ICUR,IZDIM).EQ.NEGVNM) THEN
       IG(NVar) = -I
       GO TO 20
      ELSE
       IG(NVar) = IOPT
      ENDIF
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE TOIGEO(ZMTINP,IZDIM,IGEO,IIGDIM,ZSymb,NZ)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Assigns connectivity entries in ZMTINP array to values in IGEO
C
C  ARGUMENTS
C
C  ZMTINP  -  intermediate scratch storage for the up to 7 possible
C             separate data types per Z-matrix line
C  IZDIM   -  location of entry in ZMTINP array to be assigned
C  IGEO    -  Z-matrix connectivity
C  IGDIM   -  location in IGEO array where assigned ZMTINP entry put
C  ZSymb   -  Z matrix symbols (real & dummy atoms)
C  NZ      -  total number of centres (lines) in Z-Matrix
C
C
      DIMENSION IGEO(NZ,3)
      CHARACTER*8 ZMTINP(NZ,7),ZSymb(NZ)
      CHARACTER BLANK
      PARAMETER (BLANK = ' ')
C
C
      DO 20 I=1,NZ
      IF(ZMTINP(I,IZDIM)(1:1).EQ.BLANK) THEN
       IGEO(I,IIGDIM) = 0
       GO TO 20
      ENDIF
c
      DO 10 J=1,NZ
      IF(ZMTINP(I,IZDIM).EQ.ZSymb(J)) THEN
       IGEO(I,IIGDIM) = J
       GO TO 20
      ENDIF
 10   CONTINUE
c
      CALL CHRTOI(ZMTINP(I,IZDIM),INTGR,ILEN)
      IGEO(I,IIGDIM) = INTGR
 20   CONTINUE
C
      RETURN
      END
c  =======================================================================
c
      SUBROUTINE UpDateZMAT(NZ,     XC,     ZSymb,  GEO,    IGEO,
     $                      IG,     VARNAM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine updates the values in a Z Matrix directly from
C  the Z-matrix connectivity and the current Cartesian coordinates
C  (Currently used during a potential scan)
C
C  ARGUMENTS
C
C  NZ      -  number of rows in the Z-matrix
C             (should be same as number of real atoms)
C  XC      -  Cartesian coordinates
C  ZSymb   -  Z-matrix atomic symbol
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO    -  Z-matrix connectivity
C              Atoms are usually defined with respect to previously
C              defined atoms by a stretch, a bend and a torsion;
C              an alternative is a stretch and 2 bends, or a stretch,
C              a bend and an out-of-plane bend. The fourth column
C              of IGEO is used to distinguish these cases
C                0 - bend + torsion
C                1 - 2 bends
C                2 - bend + out-of-plane bend
C      ** WARNING - CURRENTLY SET UP FOR TORSIONS ONLY **
C  IG      -  array determining what to do with Z-matrix variable
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C             For a potential scan, all variables should be fixed
C             EXCEPT the scanned one. Here we are going to update all
C             the "fixed" variables and leave the scanned one.
C  VARNAM  -  names of all Z-matrix variables (including "fixed")
C
C
      DIMENSION XC(3,NZ),GEO(NZ,3),IGEO(NZ,4),IG(3*NZ)
      CHARACTER*8 ZSymb(NZ),VARNAM(3*NZ)
      character*80 Char,Char1
      character*256 jobname
      Logical found
C
      Common/ job/ jobname,lenJ
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
C
C  Update ALL Z-matrix parameters using current Cartesian geometry
C
      ToDEG = 180.0d0/PI
      AUTOAN = 1.0d0/ANTOAU
C
C  Assign stretches
C
      DO 10 I=2,NZ
      J = IGEO(I,1)
      CALL StreGRAD(NZ,I,J,XC,Th,.false.,jnk)
      GEO(I,1) = Th*AUTOAN
 10   CONTINUE
C
C  Assign bends
C
      DO 20 I=3,NZ
      J = IGEO(I,1)
      K = IGEO(I,2)
      CALL AngGRAD(NZ,I,J,K,XC,Th,.false.,jnk)
      GEO(I,2) = Th*ToDEG
 20   CONTINUE
C
C  Assign torsions/out-of-plane bends
C
      DO 30 I=4,NZ
      J = IGEO(I,1)
      K = IGEO(I,2)
      L = IGEO(I,3)
      If(IGEO(I,4).EQ.0) Then
        CALL DihGRAD(NZ,I,J,K,L,XC,Th,.false.,jnk)
      Else
        CALL OutpGRAD(NZ,I,J,K,L,XC,Th,.false.,jnk)
      EndIF
      GEO(I,3) = Th*ToDEG
 30   CONTINUE
C
C  now write to Z Matrix
C  we are going to update all the "fixed" variables
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.zmat',
     $      FORM='FORMATTED',STATUS='UNKNOWN')
C
C  write the first two lines
C
      WRITE(40,'(a)') '$zmatrix'
      WRITE(40,900) ZSymb(1)
C
C  second line of Z Matrix
C
      Char(1:8) = ZSymb(2)
      WRITE(Char(9:14),901) IGEO(2,1)
      If(IG(1).EQ.1000) Then
        WRITE(Char(15:28),902) GEO(2,1)
      Else
        Char(15:28) = VARNAM(1)
      EndIf
      WRITE(40,'(a)') Char(1:28)
C
C  third line of Z Matrix
C
      Char(1:8) = ZSymb(3)
      WRITE(Char(9:14),901) IGEO(3,1)
      If(IG(2).EQ.1000) Then
        WRITE(Char(15:28),902) GEO(3,1)
      Else
        Char(15:28) = VARNAM(2)
      EndIf
      WRITE(Char(29:34),901) IGEO(3,2)
      If(IG(NZ).EQ.1000) Then
        WRITE(Char(35:48),902) GEO(3,2)
      Else
        Char(35:48) = VARNAM(NZ)
      EndIf
      WRITE(40,'(a)') Char(1:48)
C
C  all other Z-matrix lines
C
      IT = 2
      DO 40 I=4,NZ
      Char(1:8) = ZSymb(I)
      WRITE(Char(9:14),901) IGEO(I,1)
      IT = IT+1
      If(IG(IT).EQ.1000) Then
        WRITE(Char(15:28),902) GEO(I,1)
      Else
        Char(15:28) = VARNAM(IT)
      EndIf
      WRITE(Char(29:34),901) IGEO(I,2)
      JT = NZ-2+IT
      If(IG(JT).EQ.1000) Then
        WRITE(Char(35:48),902) GEO(I,2)
      Else
        Char(35:48) = VARNAM(JT)
      EndIf
      WRITE(Char(49:54),901) IGEO(I,3)
      JT = 2*NZ-5+IT
      If(IG(JT).EQ.1000) Then
        WRITE(Char(55:68),902) GEO(I,3)
      Else
        Char(55:68) = VARNAM(JT)
      EndIf
      WRITE(40,'(a)') Char(1:68)
 40   CONTINUE
C
C  now write variables section
C
      WRITE(40,*)
c
      DO 50 I=1,NZ-1
      If(IG(I).EQ.0) WRITE(40,910) VARNAM(I),GEO(I+1,1)
 50   CONTINUE
c
      IT = NZ-1
c
      DO 60 I=1,NZ-2
      IT = IT+1
      If(IG(IT).EQ.0) WRITE(40,910) VARNAM(IT),GEO(I+2,2)
 60   CONTINUE
c
      DO 70 I=1,NZ-3
      IT = IT+1
      If(IG(IT).EQ.0) WRITE(40,910) VARNAM(IT),GEO(I+3,3)
 70   CONTINUE
c
      WRITE(40,'(a)') '$end'
C
C  close file and exit
C
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
c
  900 Format(A8)
  901 Format(1X,I4,1X)
  902 Format(1X,F12.6,1X)
  910 Format(A8,1X,F12.6)
c
      END
c  =======================================================================
c
      SUBROUTINE WrZMAT(NZ,GEO,IG,VARNAM,delC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Updates the parameter block in the Z-Matrix file
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres in Z-Matrix
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  VARNAM  -  names of all "variables" (including fixed)
C  delC    -  logical variable, if .true. <coord> file will be deleted
C
C
      DIMENSION GEO(NZ,3),IG(*)
      CHARACTER*8 VARNAM(*),CHAR
      character*256 jobname
      Logical delC
c
      Common /job/jobname,lenJ
C
C
C  Try to open the <zmat> file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.zmat',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  position file at beginning of variable section
C
      DO 10 I=1,NZ+2
      READ(40,900) CHAR
 10   CONTINUE
C
C  run over all internal coordinates and extract and
C  write variables to ZMAT
C
      DO 20 I=1,NZ-1
      If(IG(I).EQ.0) WRITE(40,910) VARNAM(I),GEO(I+1,1)
 20   CONTINUE
c
      IT = NZ-1
c
      DO 30 I=1,NZ-2
      IT = IT+1
      If(IG(IT).EQ.0) WRITE(40,910) VARNAM(IT),GEO(I+2,2)
 30   CONTINUE
c
      DO 40 I=1,NZ-3
      IT = IT+1
      If(IG(IT).EQ.0) WRITE(40,910) VARNAM(IT),GEO(I+3,3)
 40   CONTINUE
c
      WRITE(40,'(a)') '$end'
C
C  close file and exit
C
      CLOSE (UNIT=40,STATUS='KEEP')
C
C  See if <coord> file needs to be deleted
C  (for proper working of GEOMETRY)
C
      If(delC) Then
        OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.coord',STATUS='UNKNOWN')
        CLOSE (UNIT=40,STATUS='DELETE')
      EndIf
C
      RETURN
C
 95   CONTINUE
      Call nerror(1,'Z-Matrix Write',
     $      'No Z-Matrix file found in <WrZMAT>!',0,0)
c
  900 Format(A8)
  910 Format(A8,1X,F12.6)
c
      END
