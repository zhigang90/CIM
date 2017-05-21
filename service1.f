c ==================================================================
c  SERVICE ROUTINES taken from OPTIMIZE module
c ==================================================================
c
c  (A) Matrix Diagonalization
c  --------------------------
c
      SUBROUTINE DIAGMAT(A,N,VM,V,D,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
      integer*4 i4err
C
C
C  Diagonalizes a Real Symmetric matrix A by Householder
C  reduction to Tridiagonal form
C
C  ARGUMENTS
C
C  A     -  input matrix
C           on exit contains eigenvectors
C  N     -  dimension of A
C  VM    -  scratch space for eigenvector ordering (N*N)
C  V     -  scratch space (N)
C  D     -  on exit contains eigenvalues
C
C
      DIMENSION A(N,N),VM(N,N),V(N),D(N)
      Dimension S(10)
C
C
cc      CALL TRED2b(A,N,D,V)
cc      CALL TQLIb(D,V,N,A,IErr)
cc      If(IErr.NE.0) RETURN           ! error exit
C
C  Order the eigenvalues
C
cc      CALL SortEV(N,D,A,VM,V)
C
c --------------------------------------
c -- now using fast, optimized LAPACK routine
c -- If test because need at least 3*N storage
C
      If(N.GT.3) Then
        call dsyev('V','U',N,A,N,D,VM,N*N,i4err)
      Else
        call dsyev('V','U',N,A,N,D,S,10,i4err)
      EndIf
c
c The following is even faster and better in parallel, but a lot
c more scratch is needed for the divide and conquer algorithm
cc      iw1=1+5*N+2*N*INT(log(DBLE(N))/log(2.0D0)+1)+3*N**2
cc      iw2=2+5*N
cc      call getmem(iw1,itmp1)
cc      call getmem(iw2,itmp2)
cc      call dsyevd('V','U',N,A,N,D,bl(itmp1),iw1,bl(itmp2),iw2,IErr)
cc      call retmem(2)
C
      IERR=int(i4err)
      RETURN
      END
c =================================================================
c
      SUBROUTINE TRED2b(A,N,D,E)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Householder Reduction of a real symmetric matrix to
C  tridiagonal form
C  Taken from:
C  "Numerical Recipes: The Art of Scientific Computing"
C  Press, Flannery, Teukolsky and Vetterling
C
C  ARGUMENTS
C
C  A     -  input matrix
C           on exit contains orthogonal matrix effecting
C           the transformation
C  N     -  dimension of A
C  D     -  on exit returns diagonal elements
C           of tridiagonal matrix
C  E     -  on exit returns off-diagonal elements
C           of tridiagonal matrix
C
C  This routine is used along with TQLIb to diagonalize a matrix
C  and find its eigenvectors and eigenvalues
C
C
      DIMENSION A(N,N),D(N),E(N)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0)
C
C
      IF(N.EQ.1) THEN                ! special case
       D(1) = A(1,1)
       A(1,1) = ONE
       E(1) = ZERO
       RETURN
      ENDIF
C
      DO 18 I=N,2,-1
      L = I-1
      H = ZERO
      Scal = ZERO
      IF(L.GT.1) THEN
       DO 11 K=1,L
       Scal = Scal + Abs(A(I,K))
 11    CONTINUE
       IF(Scal.EQ.ZERO) THEN          ! skip transformation
        E(I) = A(I,L)
       ELSE
        DO 12 K=1,L
        A(I,K) = A(I,K)/Scal
        H = H + A(I,K)**2
 12     CONTINUE
        F = A(I,L)
        G = -SIGN(SQRT(H),F)
        E(I) = Scal*G
        H = H - F*G
        A(I,L) = F-G
        F = ZERO
        DO 15 J=1,L
        A(J,I) = A(I,J)/H
        G = ZERO
        DO 13 K=1,J
        G = G + A(J,K)*A(I,K)
 13     CONTINUE
        IF(L.GT.J) THEN
         DO 14 K=J+1,L
         G = G + A(K,J)*A(I,K)
 14      CONTINUE
        ENDIF
        E(J) = G/H
        F = F + E(J)*A(I,J)
 15     CONTINUE
        HH = F/(H+H)
        DO 17 J=1,L
        F = A(I,J)
        G = E(J) - HH*F
        E(J) = G
        DO 16 K=1,J
        A(J,K) = A(J,K) - F*E(K) - G*A(I,K)
 16     CONTINUE
 17     CONTINUE
       ENDIF
      ELSE
       E(I) = A(I,L)
      ENDIF
      D(I) = H
 18   CONTINUE
c
      D(1) = ZERO
      E(1) = ZERO
      DO 23 I=1,N
      L = I-1
      IF(D(I).NE.ZERO) THEN
       DO 21 J=1,L
       G = ZERO
       DO 19 K=1,L
       G = G + A(I,K)*A(K,J)
 19    CONTINUE
       DO 20 K=1,L
       A(K,J) = A(K,J) - G*A(K,I)
 20    CONTINUE
 21    CONTINUE
      ENDIF
      D(I) = A(I,I)
      A(I,I) = ONE
      IF(L.GE.1) THEN
       DO 22 J=1,L
       A(I,J) = ZERO
       A(J,I) = ZERO
 22    CONTINUE
      ENDIF
 23   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE TQLIb(D,E,N,A,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine determines the eigenvectors and eigenvalues
C  of a Real Symmetric Tridiagonal Matrix
C  Taken from:
C  "Numerical Recipes: The Art of Scientific Computing"
C  Press, Flannery, Teukolsky and Vetterling
C
C  ARGUMENTS
C
C  D     -  diagonal elements of tridiagonal matrix
C           on exit contains eigenvalues
C  E     -  off-diagonal elements of tridiagonal matrix
C  N     -  dimension of A
C  A     -  on input should be unit matrix or matrix output
C           by TRED2b if this was used to get tridiagonal form
C           on exit contains eigenvectors
C  IErr  -  error flag    0 - success
C                        -1 - something's gone wrong
C
C  This routine is used along with TRED2b to diagonalize a matrix
C
C
      DIMENSION D(N),E(N),A(N,N)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0)
      PARAMETER (MaxIT=50)
C
C
      If(N.EQ.1) GO TO 95                ! special case
C
      IErr = -1
C
      DO 11 I=2,N
      E(I-1) = E(I)
 11   CONTINUE
      E(N) = ZERO
c
      DO 15 L=1,N
      ITER = 0
c
 1    CONTINUE
      DO 12 M=L,N-1
      DD = Abs(D(M)) + Abs(D(M+1))
      If(Abs(E(M))+DD.EQ.DD) GO TO 2
 12   CONTINUE
      M = N
c
 2    CONTINUE
      IF(M.NE.L) THEN
cc
       If(ITER.EQ.MaxIT) RETURN            ! Error Exit
cc
       ITER = ITER + 1
       G = (D(L+1)-D(L))/(TWO*E(L))
       R = SQRT(ONE + G*G)
       G = D(M) - D(L) + E(L)/(G + SIGN(R,G))
       S = ONE
       C = ONE
       P = ZERO
       DO 14 I=M-1,L,-1
       F = S*E(I)
       B = C*E(I)
       IF(Abs(F).GE.Abs(G)) THEN
        C = G/F
        R = SQRT(ONE + C*C)
        E(I+1) = F*R
        S = ONE/R
        C = C*S
       ELSE
        S = F/G
        R = SQRT(ONE + S*S)
        E(I+1) = G*R
        C = ONE/R
        S = S*C
       ENDIF
       G = D(I+1) - P
       R = (D(I)-G)*S + TWO*C*B
       P = S*R
       D(I+1) = G+P
       G = C*R - B
c
       DO 13 K=1,N
       F = A(K,I+1)
       A(K,I+1) = S*A(K,I) + C*F
       A(K,I) = C*A(K,I) - S*F
 13    CONTINUE
c
 14    CONTINUE
       D(L) = D(L) - P
       E(L) = G
       E(M) = ZERO
       GO TO 1
cc
      ENDIF
 15   CONTINUE
C
 95   CONTINUE
      IErr = 0
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE SortEV(N,EigVal,U,VM,INDX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Sorts all eigenvalues and associated eigenvectors of a
C  matrix into ascending order
C
C  ARGUMENTS
C
C  N       -  dimension of eigensystem
C  EigVal  -  array containing eigenvalues in random order
C  U       -  corresponding eigenvectors in columns
C  VM      -  scratch space (N*N)
C  INDX    -  scratch space for indexing array (N)
C
C
      REAL*8 EigVal(N),U(N,N),VM(N,N)
      INTEGER INDX(N)
C
C
C  first make the index table for ordering the eigenvalues
C
      CALL INDEXA(N,EigVal,INDX)
C
C  order the eigenvalues
C
      CALL CpyVEC(N,EigVal,VM(1,1))
      DO 10 I=1,N
      EigVal(I) = VM(INDX(I),1)
 10   CONTINUE
C
C  order the eigenvectors
C
      CALL CpyVEC(N*N,U,VM)
      DO 20 I=1,N
      CALL CpyVEC(N,VM(1,INDX(I)),U(1,I))
 20   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE INDEXA(N,A,INDX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Indexes an array A of length N
C  on exit the array INDX is such that A(INDX(J)) is in ascending
C  order for J=1,2...N. The initial array A is unchanged.
C
C  This routine combines a number of the sorting algorithms in
C  "Numerical Recipes: The Art of Scientific Computing"
C  Press, Flannery, Teukolsky and Vetterling
C
C  The algorithm chosen depends on the size of N
C
C
      DIMENSION A(N),INDX(N)
C
      PARAMETER (ALN2I=1.0d0/0.69314718d0,small=1.0d-5)
      PARAMETER (NSmall=50,NBig=500)
C
C
C  set up initial index array
C
      DO 10 I=1,N
      INDX(I) = I
 10   CONTINUE
C
      IF(N.LT.NSmall) THEN
cc
C  sort by straight insertion
C
       DO 30 J=2,N
       AT = A(J)
       IND = INDX(J)
       DO 20 I=J-1,1,-1
       If(A(INDX(I)).LE.AT) GO TO 21
       INDX(I+1) = INDX(I)
 20    CONTINUE
       I = 0
 21    INDX(I+1) = IND
 30    CONTINUE
cc
      ELSE IF(N.GE.NSmall.AND.N.LT.NBig) THEN
cc
C  Shell sort
C
       LOGNB2 = INT( LOG(Float(N))*ALN2I + small )
       M = N
c
       DO 60 NN=1,LOGNB2
       M = M/2
       K = N-M
       DO 50 J=1,K
       I = J
 40    CONTINUE
       L = I+M
       IF(A(INDX(L)).LT.A(INDX(I))) THEN
        IND = INDX(I)
        INDX(I) = INDX(L)
        INDX(L) = IND
        I = I-M
        If(I.GE.1) GO TO 40
       ENDIF
 50    CONTINUE
 60    CONTINUE
cc
      ELSE
cc
C  Heap sort
C
       L = N/2 + 1
       IR = N
c
 70    CONTINUE
       IF(L.GT.1) THEN
        L = L-1
        IND = INDX(L)
        AT = A(IND)
       ELSE
        IND = INDX(IR)
        AT = A(IND)
        INDX(IR) = INDX(1)
        IR = IR-1
        IF(IR.EQ.1) THEN
         INDX(1) = IND
         RETURN
        ENDIF
       ENDIF
c
       I = L
       J = L+L
c
 80    CONTINUE
       IF(J.LE.IR) THEN
        IF(J.LT.IR) THEN
         If(A(INDX(J)).LT.A(INDX(J+1))) J = J+1
        ENDIF
        IF(AT.LT.A(INDX(J))) THEN
         INDX(I) = INDX(J)
         I = J
         J = J+J
        ELSE
         J = IR+1
        ENDIF
        GO TO 80
       ENDIF
       INDX(I) = IND
       GO TO 70
cc
      ENDIF
C
      RETURN
      END
c =================================================================
c
c  (B) Matrix Inversion
c  --------------------
c
      SUBROUTINE INVMAT(A,N,VV,INDX,AINV,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Inverts a matrix A by LU decomposition and back substitution
C  **NOTE** on exit A contains the LU decomposition of the
C           original input matrix, which is destroyed
C
C  ARGUMENTS
C
C  A     -  input matrix (destroyed on exit)
C  N     -  dimension of A
C  VV    -  scratch vector (real)
C  INDX  -  scratch vector (integer)
C  AINV  -  on exit contains the matrix inverse
C  IErr  -  error flag
C             0 - matrix successfully inverted
C            -1 - LU decomposition failed
C
C
      DIMENSION A(N,N),VV(N),INDX(N),AINV(N,N)
C
      PARAMETER (ONE=1.0d0)
C
C
      CALL LUDCMP(A,N,VV,INDX,ID,IErr)
      If(IErr.NE.0) RETURN
C
C  find inverse by columns
C
      CALL SetDiagMat(N,ONE,AINV)
      DO 10 I=1,N
      CALL LUBKSB(A,N,INDX,AINV(1,I))
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE LUDCMP(A,N,VV,INDX,ID,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine does LU decomposition of a Square Matrix
C  Taken from:
C  "Numerical Recipes: The Art of Scientific Computing"
C  Press, Flannery, Teukolsky and Vetterling
C
C  ARGUMENTS
C
C  A     -  input matrix
C           on exit contains LU decomposition
C  N     -  dimension of A
C  VV    -  scratch vector to hold implicit row scaling
C  INDX  -  on exit contains row permutation effected by
C           partial pivoting
C  ID    -  on exit +/-1 depending on whether the total
C           number of row interchanges was even or odd
C           (used for sign of matrix determinant)
C  IErr  -  error flag    0 - success
C                        -1 - LU decomposition failed
C
C  This routine is used along with LUBKSB to solve linear
C  equations or to invert a matrix
C
C
      DIMENSION A(N,N),VV(N),INDX(N)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0,Tiny=1.0d-20)
C
C
      IErr = -1
      ID = 1
C
      DO 20 I=1,N
      AAMAX = ZERO
      DO 10 J=1,N
      If(Abs(A(I,J)).GT.AAMAX) AAMAX = Abs(A(I,J))
 10   CONTINUE
      If(AAMAX.LT.Tiny) RETURN         ! matrix is singular
      VV(I) = ONE/AAMAX                ! save the scaling factor
 20   CONTINUE
c
      DO 90 J=1,N
c
      DO 40 I=1,J-1
      SUM = A(I,J)
      DO 30 K=1,I-1
      SUM = SUM - A(I,K)*A(K,J)
 30   CONTINUE
      A(I,J) = SUM
 40   CONTINUE
c
      AAMAX = ZERO                     ! search for largest pivot
      DO 60 I=J,N
      SUM = A(I,J)
      DO 50 K=1,J-1
      SUM = SUM - A(I,K)*A(K,J)
 50   CONTINUE
      A(I,J) = SUM
      DUM = VV(I)*Abs(SUM)
      IF(DUM.GE.AAMAX) THEN
       IMAX = I
       AAMAX = DUM
      ENDIF
 60   CONTINUE
c
      IF(J.NE.IMAX) THEN                ! do we need to interchange rows?
       DO 70 K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
 70    CONTINUE
       ID = -ID                         ! change parity of ID
       VV(IMAX) = VV(J)                 ! also interchange scaling factor
      ENDIF
c
      INDX(J) = IMAX
      If(A(J,J).EQ.ZERO) A(J,J) = Tiny
      IF(J.NE.N) THEN
       DUM = ONE/A(J,J)
       DO 80 I=J+1,N
       A(I,J) = A(I,J)*DUM
 80    CONTINUE
      ENDIF
c
 90   CONTINUE
C
      IErr = 0
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE LUBKSB(A,N,INDX,B)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Solves the set of N linear equations:  A * X = B
C  Taken from:
C  "Numerical Recipes: The Art of Scientific Computing"
C  Press, Flannery, Teukolsky and Vetterling
C
C  ARGUMENTS
C
C  A     -  input matrix
C           **NOTE** NOT A but its LU decomposition
C           (as determined from routine LUDCMP)
C  N     -  dimension of A
C  INDX  -  permutation vector returned by LUDCMP
C  B     -  on input contains rhs vector
C           on exit contains solution vector X
C
C
C  This routine is used along with LUDCMP to solve linear
C  equations or to invert a matrix
C
C
      DIMENSION A(N,N),INDX(N),B(N)
C
      PARAMETER (ZERO=0.0d0)
C
C
      II = 0
      DO 20 I=1,N
      LL = INDX(I)
      SUM = B(LL)
      B(LL) = B(I)
      IF(II.NE.0) THEN
       DO 10 J=II,I-1
       SUM = SUM - A(I,J)*B(J)
 10    CONTINUE
      ELSE IF(SUM.NE.ZERO) THEN
       II = I
      ENDIF
      B(I) = SUM
 20   CONTINUE
c
      DO 40 I=N,1,-1
      SUM = B(I)
      IF(I.LT.N) THEN
       DO 30 J=I+1,N
       SUM = SUM - A(I,J)*B(J)
 30    CONTINUE
      ENDIF
      B(I) = SUM/A(I,I)      ! store a component of the solution vector
 40   CONTINUE
C
      RETURN
      END
c =================================================================
c
c  (C) General Utilities
c  ---------------------
c
      SUBROUTINE SetDiagMat(N,Val,A)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(N,N)
C
C  Set up a diagonal matrix with all entries
C  equal to Val
C
      CALL ZeroIT(A,N*N)          ! use old texas utility
      DO 10 I=1,N
      A(I,I) = Val
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
c
      SUBROUTINE AddVEC(N,A,B,C)
      REAL*8 A(N),B(N),C(N)
C
C  Adds two vectors:   C = A + B
C
      DO 10 I=1,N
      C(I) = A(I) + B(I)
 10   CONTINUE
c
      RETURN
      END
c =================================================================
c
      SUBROUTINE MinusVEC(N,A,B,C)
      REAL*8 A(N),B(N),C(N)
C
C  Subtracts two vectors:   C = A - B
C
      DO 10 I=1,N
      C(I) = A(I) - B(I)
 10   CONTINUE
c
      RETURN
      END
c =================================================================
      LOGICAL FUNCTION CmpVEC(N,A,B,thrsh)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Compares two vectors of length N
C  On exit returns .true. if no component of the
C  two vectors differs by more than thrsh;
C  otherwise returns .false.
C
      REAL*8 A(N),B(N)
C
      CmpVEC = .FALSE.
c
      DO 10 I=1,N
      If(Abs(A(I)-B(I)).GT.thrsh) RETURN
 10   CONTINUE
c
      CmpVEC = .TRUE.
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE CpyVEC(N,A,B)
      REAL*8 A(N),B(N)
C
C  Copies vector A into vector B
C
      DO 10 I=1,N
      B(I) = A(I)
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE ICpyVEC(N,A,B)
      INTEGER A(N),B(N)
C
C  Copies integer array A into vector B
C
      DO 10 I=1,N
      B(I) = A(I)
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE EXPAND(N,A,B)
      REAL*8 A(N*(N+1)/2),B(N,N)
C
C  Expands a lower triangle stored as a linear array in A
C  to a full square symmetric matrix stored in B
C
      IT = 0
      DO 10 I=1,N
      DO 10 J=1,I
      IT = IT+1
      B(I,J) = A(IT)
      B(J,I) = A(IT)
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE MatVEC(N,A,V1,V2)
      REAL*8 A(N,N),V1(N),V2(N),Val
C
C  Square matrix A x vector V1:   V2 = A x V1
C
      CALL ZeroIT(V2,N)
c
      DO 20 J=1,N
      Val = V1(J)
      DO 10 I=1,N
      V2(I) = V2(I) + A(I,J)*Val
 10   CONTINUE
 20   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE GetAtNo(NATOMS,AtSymb,IAN)
      IMPLICIT INTEGER(A-Z)
      Data IOut/6/
C
C  Get Atomic numbers from Atomic symbols
C
      INTEGER IAN(NATOMS)
      CHARACTER*8 AtSymb(NATOMS)
C
      DO 10 IAtm=1,NATOMS
      CALL nugrep(AtSymb(IAtm)(1:2),INum)
      IF(INum.LT.0) THEN
       WRITE(IOut,1000) IAtm,AtSymb(IAtm)(1:2)
       CALL nerror(1,'Subroutine GetAtNo',
     $               ' Unidentified Atomic Symbol',0,0)
      ENDIF
      IAN(IAtm) = INum
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unidentified Atomic Symbol',/,
     $         5X,'Atomic centre: ',I3,' Atomic Symbol: ',A2)
c
      END
c =================================================================
c
      subroutine nugrep(atosym,nuchar)
      implicit real*8 (a-h,o-z)

c ..............................................................
c      deduce nuclear charge from element symbol
c      elements included : 1-94
c         ATOSYM                    NUCHAR
c      symbol for real atom      returns atomic number
c      symbol for dummy atom     returns 0
c      unknown symbol            returns -1
c ..............................................................

      parameter ( nelem = 94 )

      character*2 atosym,nuc(nelem)

      character*26 alpbet,betalp

      data nuc /'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     1          'na','mg','al','si','p ','s ','cl','ar',
     2          'k ','ca','sc','ti','v ','cr','mn','fe','co','ni',
     3          'cu','zn','ga','ge','as','se','br','kr',
     4          'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd',
     5          'ag','cd','in','sn','sb','te','i ','xe',
     6          'cs','ba','la','ce','pr','nd','pm','sm','eu','gd',
     7          'tb','dy','ho','er','tm','yb','lu','hf','ta','w ',
     8          're','os','ir','pt','au','hg','tl','pb','bi','po',
     9          'at','rn','fr','ra','ac','th','pa','u ','np','pu' /

      data alpbet / 'abcdefghijklmnopqrstuvwxyz' /
      data betalp / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /

      i=index(betalp,atosym(1:1))
      if(i.gt.0) atosym(1:1)=alpbet(i:i)
      i=index(betalp,atosym(2:2))
      if(i.gt.0) atosym(2:2)=alpbet(i:i)

      nuchar=-1
      i=0
  100 i=i+1
      if(i.gt.nelem) goto 150
      if(atosym.ne.nuc(i)) goto 100
      nuchar=i
      return

  150 continue
c -- 2 symbols failed  If just first symbol fits, accept

      If(atosym(1:1).eq.'h') Then
       nuchar=1
      Else If(atosym(1:1).eq.'b') Then
       nuchar=5
      Else if(atosym(1:1).eq.'c') Then
       nuchar=6
      Else If(atosym(1:1).eq.'n') Then
       nuchar=7
      Else If(atosym(1:1).eq.'o') Then
       nuchar=8
      Else If(atosym(1:1).eq.'f') Then
       nuchar=9
      Else If(atosym(1:1).eq.'p') Then
       nuchar=15
      Else If(atosym(1:1).eq.'s') Then
       nuchar=16
      Else If(atosym(1:1).eq.'k') Then
       nuchar=19
      Else If(atosym(1:1).eq.'v') Then
       nuchar=23
      Else If(atosym(1:1).eq.'y') Then
       nuchar=39
      Else If(atosym(1:1).eq.'i') Then
       nuchar=53
      Else If(atosym(1:1).eq.'w') Then
       nuchar=74
      Else If(atosym(1:1).eq.'u') Then
       nuchar=92
c -- dummy atom
      Else If( atosym(1:1).eq.'x'.or.atosym(1:1).eq.'q'.or.
     $        (atosym(1:1).eq.'d'.and.atosym(2:2).eq.'u')) Then
       nuchar=0
      EndIf
c
      return
      end
c =================================================================
c
      SUBROUTINE SymbolM(NAtoms, Symb, AtSymb)
      IMPLICIT INTEGER(A-Z)
C
C  This routine massages the user-defined atom symbols
C  removing all numbers from the symbol
C
C  Symbols are typically of the form:
C    atomic symbol - string of numbers - additional characters
C  E.g.  C1A,  H13,  SI33B
C  In general numbers are used to label particular atoms (and are
C  ignored). Any characters NOT part of the real atomic symbol
C  itself are used to distinguish between different "classes" of
C  the same atom for, e.g., assigning different basis sets and
C  determining the molecular symmetry.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  Symb    -  user-defined atom symbol
C  AtSymb  -  on exit massaged atomic symbol
C
C
      CHARACTER*8 Symb(NAtoms),AtSymb(NAtoms),blank8
c
      Data blank8/'        '/
C
C
C  Initialize position of zero in collating series
C
      JZero = ICHAR('0')
C
C  Loop over all atoms
C
      DO 50 I=1,NAtoms
C
C  clear AtSymb
C
      AtSymb(I) = blank8
C
C  first symbol always needed
C
      AtSymb(I)(1:1) = Symb(I)(1:1)
C
C  Loop over the remaining  characters
C
      jt = 1
c
      DO 40 j=2,8
C
C  If character is a blank, exit
C
      If(Symb(I)(j:j).EQ.' ') exit
C
C  If character is a number then remove it
C
      ICh = ICHAR(Symb(I)(j:j)) - JZero
      If(ICh.LT.0.OR.ICh.GT.9) Then
       jt = jt+1
       AtSymb(I)(jt:jt) = Symb(I)(j:j)
      EndIf
c
 40   CONTINUE
C
C  That should be all
C
 50   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE GetSymbM(AtSymb,SymS)
      IMPLICIT INTEGER(A-Z)
C
C  detects any "special symbols" in a given atomic symbol
C
C  ARGUMENTS
C
C  AtSymb  -  input atomic symbol
C  Symb    -  on exit contains any special symbol (maximum 2)
C
C
      CHARACTER AtSymb*8,Symb*8,SymS*2
      Character*1 Special(15)
C
      Data Special/'~','!','@','#','$','%','^','&','*','+',
     $             '=','<','>','?','"'/
C
C
C  first remove any numbers
C
      Call SymbolM(1,AtSymb,Symb)
C
C  now look for special symbols
C  (first symbol cannot be special)
C
      SymS = '  '
      IT = 0
c
      DO I=2,8
      If(Symb(I:I).EQ.' ') Exit
      DO J=1,15
      If(Symb(I:I).EQ.Special(J)) Then
       IT = IT+1
       SymS(IT:IT) = Special(J)
       Exit
      EndIf
      EndDO
      If(IT.EQ.2) RETURN
      EndDO
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE AddSymbM(SymS,Symb)
      IMPLICIT REAL(A-H,O-Z)
C
C  appends a "special symbol" onto the end of a standard
C  (maximum 2 symbol) atomic symbol
C
C  ARGUMENTS
C
C  SymS    -  special symbol
C  Symb    -  standard atomic symbol
C
      CHARACTER SymS*2,Symb*4
C
C  remove all blanks and join symbols
C
      call rmblan(Symb,4,len0)
      Symb = Symb(1:len0)//SymS
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE NumFIELD(Char,Len,NumF)
      IMPLICIT INTEGER(A-Z)
C
C
C  Determines how many fields (strings separated by blanks)
C  are in a given character string
C
C  ARGUMENTS
C
C  Char    -  input character string
C  Len     -  length of charater string
C  NumF    -  on exit number of fields in string
C
C
      Character*1 Char(*)
C
C
C  initialize
C
      NumF = 0
      I = 0
c
 10   CONTINUE
      I = I+1
      If(I.GT.Len) RETURN
      If(Char(I).EQ.' ') GO TO 10
C
C  found non-blank character (start of field)
C
      NumF = NumF + 1
C
C  search for next blank
C
 20   CONTINUE
      I = I+1
      If(I.GT.Len) RETURN
      If(Char(I).NE.' ') GO TO 20
C
C  found blank (end of field)
C
      GO TO 10
c
      END
c =================================================================
c
      subroutine getnucdat(na,x,ian)

      use memory

      implicit real*8 (a-h,o-z)
c  returns the nuclear data (coordinates and atomic numbers)
c  in the arrays x and ian
c  these arrays must be predefined
c  data obtained from depository
      dimension x(3,na),ian(na)
      character*8 atsymb,symb
      equivalence (atsymb,xname)
c     common /big/bl(30000)
      call getival('inuc',inn)
      do 100 i=1,na
        x(1,i)=bl(inn+1)
        x(2,i)=bl(inn+2)
        x(3,i)=bl(inn+3)
        xname=bl(inn+4)
        call nugrep(atsymb(1:2),ian(i))   ! get atomic number from atomic symbol
        inn=inn+5
 100  continue
      return
      end
c =================================================================
c
      subroutine getatchrg(na,qa)

      use memory

      implicit real*8(a-h,o-z)
c  returns the actual atomic charge in the array qa
c  data obtained from depository
      dimension qa(na)
c     common /big/bl(30000)
      call getival('inuc',inn)
      do 100 i=1,na
        qa(i)=bl(inn)
        inn=inn+5
 100  continue
      return
      end
c =================================================================
c
      function arc1 (x,y)
      implicit real*8(a-h,o-z)

c ---------------------------------------------------------------------
c      compute arc1 = arctan(y/x)
c ---------------------------------------------------------------------

      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree

      if (abs(x).lt.1.0d-11) then
        arc1=0.5d0*PI
      else
        arc1=atan(y/x)
        if (x.lt.0.d0) arc1=arc1+PI
      end if

      return
      end
c =================================================================
c
      subroutine cross(a,b,c)
      implicit real*8 (a-h,o-z)
c ---------------------------------------------------------------------
c      vector c = vector product of vectors a and b
c ---------------------------------------------------------------------
      dimension a(3),b(3),c(3)
c
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
c
      return
      end
c =================================================================
c
      subroutine mult3 (a,b,c,n)
      implicit real*8 (a-h,o-z)

c ---------------------------------------------------------------------
c     multiplication of two matrices b,c
c      a = b*c     dim: a(3,n),b(3,3),c(3,n)
c ---------------------------------------------------------------------

      dimension a(3,n),b(3,3),c(3,n)
      do 30 j=1,n
      do 20 i=1,3
      s = 0.0d0
      do 10 k=1,3
  10  s = s + b(i,k)*c(k,j)
  20  a(i,j) = s
  30  continue

      return
      end
c ======================================================================
c
      subroutine normal(u,v,w)
      implicit real*8(a-h,o-z)

c ---------------------------------------------------------------------
c      compute normalized vector cross product w = u x v
c ---------------------------------------------------------------------

      dimension u(3),v(3),w(3)
      logical jubel

      data tol/1.d-8/

      w(1)=u(2)*v(3)-u(3)*v(2)
      w(2)=u(3)*v(1)-u(1)*v(3)
      w(3)=u(1)*v(2)-u(2)*v(1)

      call nom(w,3,jubel,tol)

      return
      end
c =================================================================
c
      subroutine nom(u,ndim,jubel,small)
      implicit real*8 (a-h,o-z)

c ---------------------------------------------------------------------
c     normalize the vector u of dimension ndim
c     logical jubel : .true.  if norm squared is larger than <small>
c                             => vector will be normalized
c                     .false. => vector will be set to zero !
c ---------------------------------------------------------------------

      dimension u(ndim)
      logical jubel

      data one/1.d0/

      x=SProd(ndim,u,u)
      jubel=x.gt.small
      if (jubel) then
        x=one/sqrt(x)
        call vscal(ndim,x,u)
      else
        call ZeroIT(u,ndim)
      endif

      return
      end
c =================================================================
c
      double precision function s2(x)
      implicit real*8 (a-h,o-z)

c ---------------------------------------------------------------------
c      s2(x) = sqrt ( 1 - x*x )
c ---------------------------------------------------------------------

      if (abs(x).gt.1.00001d0) then
       s2=0.d0
       call message('function s2',
     1 ' ARGUMENT OF FUNCTION s2 IS greater than 1.000001',0,0)
      elseif (abs(x).ge.1.d0) then
       s2=0.d0
      else
cc       s2=sqrt(1.d0-x*x)
c -- this is more accurate when x is small      ! JB  Oct 2009
       s2=sqrt( (1.0d0+x)*(1.0d0-x) )
      endif

      return
      end
c =================================================================
c
      DOUBLE PRECISION FUNCTION SProd(N,V1,V2)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     Forms scalar (dot) product of two vectors
C
      REAL*8 V1(N),V2(N)
C
      PARAMETER (ZERO=0.0d0)
c
      SProd = ZERO
      DO 10 I=1,N
      SProd = SProd + V1(I)*V2(I)
 10   CONTINUE
c
      RETURN
      END
c =================================================================
c
      subroutine vecdif(u,r,c1,c2)
      implicit real*8(a-h,o-z)

c ---------------------------------------------------------------------
c     compute the difference u between vectors c1 and c2, u = c1-c2
c     (all vectors of dimension 3) and return the norm of u on r
c ---------------------------------------------------------------------

      dimension u(3),c1(3),c2(3)
      logical jubel

      data tol/1.d-08/

      u(1)=c1(1)-c2(1)
      u(2)=c1(2)-c2(2)
      u(3)=c1(3)-c2(3)
      r=sqrt(SProd(3,u,u))
      call nom(u,3,jubel,tol)

      return
      end
c =================================================================
c
      SUBROUTINE VScal(N,skal,V)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Scales all elements of a vector by skal
C    V = V * skal
C
      REAL*8 V(N)
C
      DO 10 I=1,N
      V(I) = V(I)*skal
 10   CONTINUE
C
      RETURN
      END
c =================================================================
c
      SUBROUTINE PrntCAR(IOut,IUnit,NATOMS,AtSymb,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Prints out the geometry in Cartesian coordinates
C  (in Angstroms)
C
C  ARGUMENTS
C
C  IOut    -  output stream (usually 6 or 8)
C  IUnit   -  Unit flag
C              0 - coordinates in atomic units
C              1 - coordinates in angstroms
C  NATOMS  -  number of atoms
C  AtSymb  -  atomic symbols
C  XC      -  Cartesian coordinates
C
      REAL*8 XC(3,NATOMS)
      CHARACTER*8 AtSymb(NATOMS)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
      PARAMETER (One=1.0d0)
C
C
      WRITE(IOut,1000)
c
      SKAL = One/ANTOAU
      If(IUnit.EQ.1) SKAL = One
c
      DO 10 IAtm=1,NATOMS
      CX = XC(1,IAtm)*SKAL
      CY = XC(2,IAtm)*SKAL
      CZ = XC(3,IAtm)*SKAL
      WRITE(IOut,1100) IAtm,AtSymb(IAtm),CX,CY,CZ
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,22X,' Coordinates (Angstroms)',/,
     $          5X,'ATOM',14X,'X',11X,'Y',11X,'Z')
 1100 FORMAT(2X,I3,2x,A8,3(2X,F10.6))
c
      END
c =================================================================
c
      SUBROUTINE PrntGRD(IOut,NATOMS,AtSymb,GC)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 GC(3,NATOMS)
      CHARACTER*8 AtSymb(NATOMS)
C
C  Prints out the forces in Cartesian coordinates
C  (in atomic units)
C
      WRITE(IOut,1000)
c
      DO 10 IAtm=1,NATOMS
      CX = -GC(1,IAtm)
      CY = -GC(2,IAtm)
      CZ = -GC(3,IAtm)
      WRITE(IOut,1100) IAtm,AtSymb(IAtm),CX,CY,CZ
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,25X,' Cartesian Forces (au)',/,
     $          5X,'ATOM',15X,'X',12X,'Y',12X,'Z')
 1100 FORMAT(2X,I3,2x,A8,3(2X,F11.7))
c
      END
c =================================================================
c
      SUBROUTINE PrntMAT(N,NRow,NCol,A)
      REAL*8 A(NRow,NCol)
C
C  Prints out N columns of an NRow * NCol matrix A
C
      PARAMETER (maxcol=6)
      Data iout/6/
C
      NP = N
      If(NP.GT.NCol) NP = NCol     ! can't print more than NCol
c
      NT = NP/maxcol
      If(NT.EQ.0) GO TO 30
c
      DO 20 I=1,NT
      Imin = (I-1)*maxcol + 1
      Imax = I*maxcol
      write(iout,1000)
      DO 10 J=1,NRow
      write(iout,1100) (A(J,K),K=Imin,Imax)
 10   CONTINUE
 20   CONTINUE
c
 30   CONTINUE
      NS = NT*maxcol
      NLeft = NP - NS
      If(NLeft.EQ.0) RETURN
c
      write(iout,1000)
      DO 40 J=1,NRow
      write(iout,1100) (A(J,K),K=NS+1,NP)
 40   CONTINUE
C
      RETURN
c
 1000 FORMAT(/)
 1100 FORMAT(1X,6F12.6)
c
      END
c
c ========================================================================
c  MEMORY HANDLING AND ERROR CHECKING ROUTINES taken from OPTIMIZE module
c ========================================================================
c
      SUBROUTINE MemCHK(IMem,IEnd,Length,subR)
      CHARACTER*(*) subR
      Data IOut/6/
C
C  Checks that available scratch memory is not exceeded
C  during calculation of pointers
C
C  ARGUMENTS
C
C  IMem    -  number of double words available
C  IEnd    -  number of double words needed
C  Length  -  number of characters in subroutine name
C  subR    -  the subroutine in which the memory is being divided up
C           (used for error printout)
C
      IF(IEnd.GT.IMem) THEN
       WRITE(IOut,1000) subR(1:Length),IMem,IEnd
       CALL OptExit(9)
      ENDIF
c
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Insufficient Scratch Memory in',
     $            ' Subroutine ',A12,/,
     $            '     Available: ',I12,'   Requested: ',I12)
c
      END
c ======================================================================
c
      SUBROUTINE MemERR(NBytes,Length,subR)
      CHARACTER*(*) subR
      Data IOut/6/
C
C  Subroutine "subR" was unable to allocate NBytes of memory
C
C  ARGUMENTS
C
C  NBytes  -  number of bytes requested
C  Length  -  number of characters in subroutine name
C  subR    -  the subroutine in which the request was made
C             (used for error printout)
C
      WRITE(IOut,1000) NBytes,subR(1:Length)
      CALL OptExit(9)
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unable to allocate ',I12,' bytes ',
     $            ' in Subroutine ',A12)
c
      END
c ======================================================================
c
      SUBROUTINE OptExit(NExit)
      INTEGER NExit
C
C  Exit Handler
C  < to be determined>
C
cc      If(NExit.EQ.9) CALL WrExit(NExit)
      CALL Exit(NExit)
C
      RETURN
      END
c ======================================================================
c
      SUBROUTINE WrExit(IStatus)
      IMPLICIT INTEGER(A-Z)
C
C  Writes IStatus flag from optimization to temporary file
C  "converged" to indicate job termination
C  (needed due to problems with status flag on IBM)
C
C
      OPEN (UNIT=40,FILE='converged',FORM='FORMATTED',
     $      STATUS='UNKNOWN')
      WRITE(40,900) IStatus
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
c
  900 Format(I6)
c
      END
c
c =======================================================================
c  A number of routines, originally developed for the OPTIMIZE module
c  but also finding general use elsewhere in the PQS program
c =======================================================================
c
      integer function ioutfil(file)
      character*(*) file
      ioutfil=igetival(file)
      end
c .......................................................................
c
      SUBROUTINE CMS(NAtoms,X,Y,Z,XC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Transform into centre of mass coordinate system
C  (all atoms assumed to have unit mass)
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  X       -  on exit contains X centre-of-mass in original coordinate frame
C  Y       -  on exit contains Y centre-of-mass in original coordinate frame
C  Z       -  on exit contains Z centre-of-mass in original coordinate frame
C  XC      -  on input  original coordinates
C             on exit   cms coordinates
C
      REAL*8 XC(3,NAtoms)
C
      PARAMETER (Zero=0.0d0)
C
C
      X = Zero
      Y = Zero
      Z = Zero
      DO 10 IAtm=1,NAtoms
      X = X + XC(1,IAtm)
      Y = Y + XC(2,IAtm)
      Z = Z + XC(3,IAtm)
 10   CONTINUE
      X = X/DFloat(NAtoms)
      Y = Y/DFloat(NAtoms)
      Z = Z/DFloat(NAtoms)
c
      DO 20 IAtm=1,NAtoms
      XC(1,IAtm) = XC(1,IAtm) - X
      XC(2,IAtm) = XC(2,IAtm) - Y
      XC(3,IAtm) = XC(3,IAtm) - Z
 20   CONTINUE
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE GetAtSym(NZ,ZSymb,AtSymb)
      IMPLICIT INTEGER(A-Z)
C
C
C  Extracts Atomic symbols from symbols given to atoms in
C  the Z-matrix
C
C  ARGUMENTS
C
C  NZ      -  total number of atomic centres
C  ZSymb   -  Z-matrix atomic symbols
C  AtSymb  -  on exit contains real atomic symbols
C             (including dummy symbol for dummy atoms)
C
C
      CHARACTER*8 ZSymb(NZ),AtSymb(NZ),Blank
      CHARACTER atosym*2
C
      DATA Blank/'        '/
C
C  ...........................................................
C  ** RULES **
C  Z-matrix symbols are such that if either the first OR
C  the first two symbols are those for a genuine atom then
C  the centre will be interpreted as such.
C
C  Symbols beginning with the characters 'X' or 'Q' will be
C  interpreted as dummy atoms EXCEPT for 'XE' which will be
C  taken as Xenon.
C  ...........................................................
C
C
C  Initialize position of zero in Fortran collating series.
C
      JZero = ICHAR('0')
c
      DO 10 I=1,NZ
      AtSymb(I) = Blank
C
C  see if the first two symbols make sense
C
      atosym = ZSymb(I)(1:2)
C
C  second symbol MAY be a number, if so set to blank
C
      ICh = ICHAR(atosym(2:2)) - JZero
      If(ICh.GE.0.AND.ICh.LE.9) atosym(2:2) = ' '
c
      CALL nugrep(atosym,INuc)
c
      IF(INuc.EQ.-1) THEN
C
C  try the first symbol only
C
       atosym = ZSymb(I)(1:1)//' '
       CALL nugrep(atosym,INuc)
      ENDIF
C
C  trap dummy atoms
C
      IF(INuc.EQ.0) THEN
        atosym = 'Du'
      ELSE IF(INuc.EQ.-1) THEN
        IOut = ioutfil('iout')
        WRITE(IOut,1000) I
        CALL OptExit(9)
      ENDIF
c
      AtSymb(I)(1:2) = atosym
 10   CONTINUE
C
      RETURN
c
 1000 FORMAT(/,2X,'***ERROR*** Unknown Atomic Symbol in Line ',I5,
     $            ' of Z-Matrix')
c
      END
c .......................................................................
c
      SUBROUTINE SymVEC(NAtoms, NTrans, NEqATM, TRANS,  SCR,
     $                  thrsh,  V)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Symmetrizes a Cartesian vector (typically the coordinates,
C  the gradient or a Hessian eigenvector) according to all
C  operations of the molecular point group
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for symmetry operations
C  SCR     -  scratch space (3*NATOMS)
C  thrsh   -  threshold below which elements will be set to zero
C  V       -  vector to be symmetrized
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),
     $          SCR(3,NAtoms),V(3,NATOMS)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0)
C
C
      If(NTrans.EQ.1) GO TO 25
c
      CALL CpyVEC(3*NAtoms,V,SCR)
      DO 10 IAtm=1,NAtoms
      DO 10 IOP=2,NTrans
      V1 = TRANS(1,1,IOP)*V(1,IAtm) + TRANS(1,2,IOP)*V(2,IAtm)
     $                    + TRANS(1,3,IOP)*V(3,IAtm)
      V2 = TRANS(2,1,IOP)*V(1,IAtm) + TRANS(2,2,IOP)*V(2,IAtm)
     $                    + TRANS(2,3,IOP)*V(3,IAtm)
      V3 = TRANS(3,1,IOP)*V(1,IAtm) + TRANS(3,2,IOP)*V(2,IAtm)
     $                    + TRANS(3,3,IOP)*V(3,IAtm)
      ISYM = NEqATM(IAtm,IOP)
      SCR(1,ISYM) = SCR(1,ISYM) + V1
      SCR(2,ISYM) = SCR(2,ISYM) + V2
      SCR(3,ISYM) = SCR(3,ISYM) + V3
 10   CONTINUE
c
      skal = One/DFloat(NTrans)
      DO 20 IAtm=1,NAtoms
      V(1,IAtm) = SCR(1,IAtm)*skal
      V(2,IAtm) = SCR(2,IAtm)*skal
      V(3,IAtm) = SCR(3,IAtm)*skal
 20   CONTINUE
C
 25   CONTINUE
C
C  zero out elements below thrsh
C
      DO 30 IAtm=1,NAtoms
      DO 30 J=1,3
      If(Abs(V(J,IAtm)).LT.thrsh) V(J,IAtm) = Zero
 30   CONTINUE
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE StreGRAD(NAtoms,I,J,XC,dd,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the stretch I-J
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  first atom in stretch
C  J       -  central atom in stretch
C  XC      -  Cartesian coordinates
C  dd      -  on exit contains stretch distance
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate stretch gradient
C              .false. -  skip gradient calculation
C  G       -  on exit contains stretch gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms)
      LOGICAL grd
C
      PARAMETER (One=1.0d0)
C
C
      XIJ = XC(1,I) - XC(1,J)
      YIJ = XC(2,I) - XC(2,J)
      ZIJ = XC(3,I) - XC(3,J)
      dd = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
c
      If(.NOT.grd) RETURN
c
      DCS = One/dd
      G(1,I) = -DCS*XIJ
      G(1,J) =  DCS*XIJ
      G(2,I) = -DCS*YIJ
      G(2,J) =  DCS*YIJ
      G(3,I) = -DCS*ZIJ
      G(3,J) =  DCS*ZIJ
c
      RETURN
      END
c .......................................................................
c
      SUBROUTINE AngGRAD(NAtoms,I,J,K,XC,Th,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the bond angle I-J-K
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  first atom in bond angle
C  J       -  central atom in bond angle
C  K       -  third atom in bond angle
C  XC      -  Cartesian coordinates
C  Th      -  on exit contains bond angle
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate angle gradient
C              .false. -  skip gradient calculation
C  G       -  on exit contains angle gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms)
      LOGICAL grd
C
      PARAMETER (One=1.0d0,Two=2.0d0,small=1.0d-6)
C
C
      XIJ = XC(1,I) - XC(1,J)
      XIK = XC(1,I) - XC(1,K)
      XJK = XC(1,J) - XC(1,K)
      YIJ = XC(2,I) - XC(2,J)
      YIK = XC(2,I) - XC(2,K)
      YJK = XC(2,J) - XC(2,K)
      ZIJ = XC(3,I) - XC(3,J)
      ZIK = XC(3,I) - XC(3,K)
      ZJK = XC(3,J) - XC(3,K)
      RIJ = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
      RIK = SQRT(XIK*XIK + YIK*YIK + ZIK*ZIK)
      RJK = SQRT(XJK*XJK + YJK*YJK + ZJK*ZJK)
c
      D1 = Two*RIJ*RJK
      D2 = (RIJ*RIJ + RJK*RJK - RIK*RIK)
c
      CosTh = D2/D1
      If(Abs(CosTh).GT.One) CosTh = SIGN(One,CosTh)
      SinTh = SQRT(One - CosTh*CosTh)
      Th = ACOS(CosTh)
c
      If(.NOT.grd) RETURN
C
C  .............................................................
C    ** WARNING  **
C  If the three atoms are linear there are problems with the
C  angle bend derivative.  Skip derivative evaluation
C
      If(SinTh.LT.small) RETURN
C  .............................................................
C
      CALL ZeroIT(G,3*NAtoms)
c
      DCB = (Two/SinTh)/(D1*D1)
      RJKIJ = RJK/RIJ
      RIJJK = RIJ/RJK
c
      G(1,I) = -DCB*( D1*XJK + D2*XIJ*RJKIJ )
      G(1,J) =  DCB*( D1*(XJK-XIJ) - D2*(XJK*RIJJK - XIJ*RJKIJ) )
      G(1,K) =  DCB*( D1*XIJ + D2*XJK*RIJJK )
      G(2,I) = -DCB*( D1*YJK + D2*YIJ*RJKIJ )
      G(2,J) =  DCB*( D1*(YJK-YIJ) - D2*(YJK*RIJJK - YIJ*RJKIJ) )
      G(2,K) =  DCB*( D1*YIJ + D2*YJK*RIJJK )
      G(3,I) = -DCB*( D1*ZJK + D2*ZIJ*RJKIJ )
      G(3,J) =  DCB*( D1*(ZJK-ZIJ) - D2*(ZJK*RIJJK - ZIJ*RJKIJ) )
      G(3,K) =  DCB*( D1*ZIJ + D2*ZJK*RIJJK )
c
      RETURN
      END
c .......................................................................
c
      SUBROUTINE OutpGRAD(NAtoms,I,J,K,L,XC,outp,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the out-of-plane-bend I-J-K-L
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  out-of-plane atom in out-of-plane bend
C  J       -  second atom in out-of-plane bend
C  K       -  third atom in out-of-plane bend
C  L       -  central atom in out-of-plane bend
C  XC      -  Cartesian coordinates
C  outp    -  on exit contains out-of-plane bend
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate outp gradient   NOT IMPLEMENTED
C              .false. -  skip gradient calculation
C  G       -  on exit contains out-of-plane bend gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms)
      LOGICAL grd
C
      DIMENSION R1(3),R2(3),R3(3),R4(3)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,small=1.0d-6)
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
C
C
      outp = Zero
c
      CALL VecDIF(R1,rn,XC(1,I),XC(1,L))
      CALL VecDIF(R2,rn,XC(1,J),XC(1,L))
      CALL VecDIF(R3,rn,XC(1,K),XC(1,L))
      CO = SProd(3,R2,R3)
      SI = S2(CO)
c
c -- three atoms defining plane are linear
      If(SI.LT.small) RETURN
c
      CALL Normal(R2,R3,R4)
      CP = SProd(3,R1,R4)
      SJ = S2(CP)
c
      CX = -One
      If(CP.LT.Zero) CX = One
      outp = -CX*ACOS(SJ)
C
C  As defined here the out-of-plane bend cannot have a magnitude
C  greater than 90 degrees; we want it to go up to 180 degrees.
C  We use the torsion ILJK to decide this; if magnitude of torsion
C  is greater than 90 degrees then so is the magnitude of the
C  out-of-plane bend
C
cc      Call DihGRAD(NAtoms,I,L,J,K,XC,Dih,.False.,Jnk)
cc      If(ABS(Dih).GT.0.5d0*PI) Then
cc        If(outp.GT.ZERO) outp = PI-outp
cc        If(outp.LT.ZERO) outp = -(PI+outp)
cc      EndIf
C  
C
      If(.NOT.grd) RETURN
C
C
C  set terms for out-of-plane bend derivatives
C  ** NOT IMPLEMENTED **
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE DihGRAD(NAtoms,I,J,K,L,XC,Dih,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the dihedral angle I-J-K-L
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  first atom in torsion
C  J       -  second atom in torsion
C  K       -  third atom in torsion
C  L       -  fourth atom in torsion
C  XC      -  Cartesian coordinates
C  Dih     -  on exit contains dihedral angle
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate dihedral gradient
C              .false. -  skip gradient calculation
C  G       -  on exit contains dihedral angle gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms)
      LOGICAL grd
C
      DIMENSION R1(3),R2(3),R3(3),R4(3),R5(3)
C
      PARAMETER (ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0,small=1.0d-6)
C
C
      Dih = ZERO
c
      XIJ = XC(1,I) - XC(1,J)
      XIK = XC(1,I) - XC(1,K)
      XJK = XC(1,J) - XC(1,K)
      XJL = XC(1,J) - XC(1,L)
      XKL = XC(1,K) - XC(1,L)
      YIJ = XC(2,I) - XC(2,J)
      YIK = XC(2,I) - XC(2,K)
      YJK = XC(2,J) - XC(2,K)
      YJL = XC(2,J) - XC(2,L)
      YKL = XC(2,K) - XC(2,L)
      ZIJ = XC(3,I) - XC(3,J)
      ZIK = XC(3,I) - XC(3,K)
      ZJK = XC(3,J) - XC(3,K)
      ZJL = XC(3,J) - XC(3,L)
      ZKL = XC(3,K) - XC(3,L)
      RIJ = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
      RIK = SQRT(XIK*XIK + YIK*YIK + ZIK*ZIK)
      RJK = SQRT(XJK*XJK + YJK*YJK + ZJK*ZJK)
      RJL = SQRT(XJL*XJL + YJL*YJL + ZJL*ZJL)
      RKL = SQRT(XKL*XKL + YKL*YKL + ZKL*ZKL)
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
      CALL Cross(R2,R1,R4)
      CALL Cross(R3,R2,R5)
c
      CosIJK = (RIJ*RIJ + RJK*RJK - RIK*RIK)/(TWO*RIJ*RJK)
      If(Abs(CosIJK).GT.ONE) CosIJK = SIGN(ONE,CosIJK)
      SinIJK = SQRT(ONE - CosIJK*CosIJK)
      CosJKL = (RJK*RJK + RKL*RKL - RJL*RJL)/(TWO*RJK*RKL)
      If(Abs(CosJKL).GT.ONE) CosJKL = SIGN(ONE,CosJKL)
      SinJKL = SQRT(ONE - CosJKL*CosJKL)
C
C  ..............................................................
C    **  WARNING  **
C  If any three atoms are linear there are problems defining
C  the related dihedral angle.  Set the angle to zero and
C  skip derivative evaluation
C
      If(SinIJK.LT.small.OR.SinJKL.LT.small) RETURN
C  ..............................................................
C
      CosDih = SProd(3,R4,R5)/(RIJ*RJK*RJK*RKL*SinIJK*SinJKL)
      If(Abs(CosDih).GT.ONE) CosDih = SIGN(ONE,CosDih)
      Dih = ACOS(CosDih)
C
C  get sign
C
      CALL Cross(R4,R5,R1)
      Sig = SProd(3,R1,R2)
      If(Sig.LT.small) Dih = -Dih
C
      If(.NOT.grd) RETURN
C
C
C  set terms for torsional derivatives
C
      CALL ZeroIT(G,3*NAtoms)
C
C  first set the derivatives of the numerator
C        d/dx ( SProd(3,R4,R5) )
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
     $        -YIK*(YKL*ZJK-YJK*ZKL) + (XIJ+XJK)*(XJK*ZKL-XKL*ZJK)
      DT1ZK = -XIJ*(-XKL*ZJK+XJK*ZKL) + YIJ*(YKL*ZJK-YJK*ZKL)
     $        -YJL*(YJK*ZIJ-YIJ*ZJK) + (XJK+XKL)*(XIJ*ZJK-XJK*ZIJ)
      DT1ZL = -XJK*(-XJK*ZIJ+XIJ*ZJK) + YJK*(YJK*ZIJ-YIJ*ZJK)
C
C  now set the derivatives of the denominator
C        d/dx ( RIJ*RJK*RJK*RKL*SinIJK*SinJKL )
C
C  this is more complicated
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
      RIJK = CosIJK/(RIJ*RJK)
      RIJ2 = CosIJK*CosIJK/(RIJ*RIJ)
      RJK2 = CosIJK*CosIJK/(RJK*RJK)
c
      DdXI = ( -(XIJ-XIK)*RIJK + XIJ*RIJ2 )/SinIJK
      DdXJ = (  (XIJ-XJK)*RIJK - XIJ*RIJ2 + XJK*RJK2 )/SinIJK
      DdXK = ( -(XIK-XJK)*RIJK - XJK*RJK2 )/SinIJK
c
      DdYI = ( -(YIJ-YIK)*RIJK + YIJ*RIJ2 )/SinIJK
      DdYJ = (  (YIJ-YJK)*RIJK - YIJ*RIJ2 + YJK*RJK2 )/SinIJK
      DdYK = ( -(YIK-YJK)*RIJK - YJK*RJK2 )/SinIJK
c
      DdZI = ( -(ZIJ-ZIK)*RIJK + ZIJ*RIJ2 )/SinIJK
      DdZJ = (  (ZIJ-ZJK)*RIJK - ZIJ*RIJ2 + ZJK*RJK2 )/SinIJK
      DdZK = ( -(ZIK-ZJK)*RIJK - ZJK*RJK2 )/SinIJK
c
      RJKL = CosJKL/(RJK*RKL)
      RJK2 = CosJKL*CosJKL/(RJK*RJK)
      RKL2 = CosJKL*CosJKL/(RKL*RKL)
c
      DeXJ = ( -(XJK-XJL)*RJKL + XJK*RJK2 )/SinJKL
      DeXK = (  (XJK-XKL)*RJKL - XJK*RJK2 + XKL*RKL2 )/SinJKL
      DeXL = ( -(XJL-XKL)*RJKL - XKL*RKL2 )/SinJKL
c
      DeYJ = ( -(YJK-YJL)*RJKL + YJK*RJK2 )/SinJKL
      DeYK = (  (YJK-YKL)*RJKL - YJK*RJK2 + YKL*RKL2 )/SinJKL
      DeYL = ( -(YJL-YKL)*RJKL - YKL*RKL2 )/SinJKL
c
      DeZJ = ( -(ZJK-ZJL)*RJKL + ZJK*RJK2 )/SinJKL
      DeZK = (  (ZJK-ZKL)*RJKL - ZJK*RJK2 + ZKL*RKL2 )/SinJKL
      DeZL = ( -(ZJL-ZKL)*RJKL - ZKL*RKL2 )/SinJKL
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
      DCT = SIGN(ONE,Dih)/(A1*A1)
c
      G(1,I) = DCT*(A1*DT1XI - A2*DT2XI)
      G(1,J) = DCT*(A1*DT1XJ - A2*DT2XJ)
      G(1,K) = DCT*(A1*DT1XK - A2*DT2XK)
      G(1,L) = DCT*(A1*DT1XL - A2*DT2XL)
      G(2,I) = DCT*(A1*DT1YI - A2*DT2YI)
      G(2,J) = DCT*(A1*DT1YJ - A2*DT2YJ)
      G(2,K) = DCT*(A1*DT1YK - A2*DT2YK)
      G(2,L) = DCT*(A1*DT1YL - A2*DT2YL)
      G(3,I) = DCT*(A1*DT1ZI - A2*DT2ZI)
      G(3,J) = DCT*(A1*DT1ZJ - A2*DT2ZJ)
      G(3,K) = DCT*(A1*DT1ZK - A2*DT2ZK)
      G(3,L) = DCT*(A1*DT1ZL - A2*DT2ZL)
C
      RETURN
      END
c .......................................................................
c
      SUBROUTINE ReadBinary1(IUnit,N,A)
      IMPLICIT REAL*8(A-H,O-Z)
c
c -- binary read of part/all of a one-dimensional matrix
c
      DIMENSION A(N)
c
      READ(IUnit) A
c
      RETURN
      END
c .......................................................................
c
      SUBROUTINE WriteBinary1(IUnit,N,A)
      IMPLICIT REAL*8(A-H,O-Z)
c
c -- binary write of part/all of a one-dimensional matrix
c
      DIMENSION A(N)
c
      WRITE(IUnit) A
c
      RETURN
      END
c .......................................................................
c
      SUBROUTINE SAVE_POT(NDim,ext,A)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Saves the FTC Coulomb potential over the full FTC grid
C
C  ARGUMENTS
C
C  NDim    -  total dimension of FTC grid
C  ext     -  file extension
C             either "potS" (smooth) or "potM" (mixed)
C  A       -  potential to be saved
C
C
      REAL*8 A(NDim)
      CHARACTER ext*4,jobname*256
      Parameter (IUnit=1)
      Common /job/jobname,lenJ
C
C  Open the file
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.'//ext,
     $      FORM='UNFORMATTED',STATUS='UNKNOWN')
c
      WRITE(IUnit) A
c
      CLOSE (UNIT=IUnit,STATUS='keep')
c
      RETURN
      END
c .......................................................................
c
      subroutine prntbas(ncs,nsh,ncf,nbf,inx,basdat,dens)
      implicit real*8(a-h,o-z)
      dimension inx(12,ncs),basdat(13,nsh)
      dimension dens(*)
c
      write(6,*) ' ncs:',ncs,' nsh:',nsh,' ncf:',ncf,' nbf:',nbf
      write(6,*) ' INX array is:'
      do i=1,12
      write(6,*) (inx(i,j),j=1,ncs)
      enddo
      write(6,*) ' BASDAT array is:'
      call prntmat(nsh,13,nsh,basdat)
      write(6,*) ' density matrix is:'
      ij=0
      do i=1,ncf
      do j=1,i
      ij = ij+1
      write(6,*) i,'  ',j,'  ',dens(ij)
      enddo
      enddo
c
      RETURN
      end
c .......................................................................
c
      subroutine prntprim(nprim,ktyp,klist,xprim)
      implicit real*8(a-h,o-z)
c
c -- prints out current values of primitive internals
c
      dimension ktyp(nprim),klist(4,nprim),xprim(nprim)
      parameter (ToAng=1.0d0/1.88972687777435527243d0)
      parameter (ToDeg=180.0d0/3.14159265358979323844d0)
c
      Do i=1,nprim
      If(ktyp(i).eq.1) then
        write(6,1001) klist(1,i),klist(2,i),xprim(i)*ToAng
      Else If(ktyp(i).eq.2) then
        write(6,1002) klist(1,i),klist(2,i),klist(3,i),
     $                xprim(i)*ToDeg
      Else If(ktyp(i).eq.3) then
        write(6,1003) klist(1,i),klist(2,i),klist(3,i),klist(4,i),
     $                xprim(i)*ToDeg
      Else If(ktyp(i).eq.4) then
        write(6,1004) klist(1,i),klist(2,i),klist(3,i),klist(4,i),
     $                xprim(i)*ToDeg
      EndIf
      EndDo
c
 1001 Format(1X,'Stretch:   ',2I4,F12.6)
 1002 Format(1X,'Bend:      ',3I4,F10.3)
 1003 Format(1X,'OOP Bend:  ',4I4,F10.3)
 1004 Format(1X,'Torsion:   ',4I4,F10.3)
c
      return
      end
