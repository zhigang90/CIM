      SUBROUTINE PreGRID(NAtoms, XNuc,   IAN,    WXXA,   WWTA,
     $                   DISTN,  RDIST,  AIJ,    XXA,    WTA)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  calculates inverse interatomic distance and Becke AIJ
C  arrays prior to grid construction
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XNuc    -  nuclear coordinates
C  IAN     -  charge (atomic numbers)
C  WXXA    -  work array for angular grid
C  WWTA    -  work array for angular weights
C
C  on exit
C
C  DISTN   -  distance to nearest neighbour
C  RDIST   -  inverse interatomic distances
C  AIJ     -  Becke aij coefficients
C  XXA     -  angular quadrature grid points (ALL orders)
C  WTA     -  angular quadrature weights
C             (** NOTE:  grid points are distance presorted **)
C
C
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms)
      dimension IAN(NAtoms),xnuc(3,NAtoms)
      dimension WXXA(3,1130),WWTA(1130)
      dimension XXA(3,1130),WTA(1130)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
C
C
      DO 20 J=1,NAtoms
        call braggslaterradius(IAN(j),rj)
        AIJ(J,J) = Zero
        RDIST(J,J) = Zero
        DO 10 I=1,J-1
C
C  calculation of aij parameter - Becke article, Eq. (A2)
C
        call braggslaterradius(IAN(i),ri)
        thekhi = SQRT(ri/rj)      ! NOTE - Ahlrichs modification
        AIJ(i,j) = (thekhi-One)/(thekhi+One)
        AIJ(i,j) = AIJ(i,j)/(AIJ(i,j)**2 - One)
cc  |aij| > 1/2 is forbidden
        If(Abs(AIJ(i,j)).GT.Half) AIJ(i,j) = SIGN(Half,AIJ(i,j))
c
        thekhi = SQRT(rj/ri)      ! NOTE - Ahlrichs modification
        AIJ(j,i) = (thekhi-One)/(thekhi+One)
        AIJ(j,i) = AIJ(j,i)/(AIJ(j,i)**2 - One)
cc  |aij| > 1/2 is forbidden
        If(Abs(AIJ(j,i)).GT.Half) AIJ(j,i) = SIGN(Half,AIJ(j,i))
C
C  inverse interatomic distance array
C
        RDIST(I,J) = One/SQRT( (XNuc(1,I)-XNuc(1,J))**2 +
     $                         (XNuc(2,I)-XNuc(2,J))**2 +
     $                         (XNuc(3,I)-XNuc(3,J))**2 )
        RDIST(J,I) = RDIST(I,J)
c
 10     CONTINUE
 20   CONTINUE
C
C  for each atom, find distance to nearest neighbour
C
      DO 40 I=1,NAtoms
      Dist = RDIST(I,1)
        DO 30 J=1,NAtoms
        If(RDIST(J,I).GT.Dist) Dist = RDIST(J,I)
 30     CONTINUE
      DISTN(I) = One/Dist
 40   CONTINUE
C
C  precompute generic angular grids and weights
C
      CALL AngGRID(WTA,XXA)
cc      CALL AngGRID(WWTA,WXXA)
cc      write(6,*) ' Angular grid points before sorting'
cc      do i=1,1130
cc      write(6,*) i,wxxa(1,i),wxxa(2,i),wxxa(3,i)
cc      enddo
cc      CALL SortGRID(WWTA,WXXA,WTA,XXA)
cc      write(6,*) ' Angular grid points after sorting'
cc      do i=1,1130
cc      write(6,*) i,xxa(1,i),xxa(2,i),xxa(3,i)
cc      enddo
C
      RETURN
      END
********************************************************************
      SUBROUTINE SortGRID(WWTA,WXXA,WTA,XXA)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Orders angular quadrature grid points so that the sum of the
C  distances between each point and all other points preceeding
C  it in the array XXA is a minimum
C
C  ARGUMENTS
C
C  WWTA    -  angular quadrature weight
C  WXXA    -  angular points on spherical surface
C  on exit
C
C  WTA     -  sorted weights
C  XXA     -  distance sorted points
C
C
      DIMENSION WWTA(1130),WXXA(3,1130),WTA(1130),XXA(3,1130)
      dimension INDX(8),RDIST(434,434),IPNT(434),ILST(434)
c
      data INDX / 1, 15, 41, 91, 201, 395, 697, 1131/
      data big/10000.0d0/, Zero/0.0d0/
C
C
C  sort each of the 7 possible angular shells
C
      DO 50 IOdr=1,7
      I1 = INDX(IOdr)
      I2 = INDX(IOdr+1)-1
ccccc
      II=0
      DO I=I1,I2
      II = II+1
      JJ=0
      DO J=I1,I-1
      JJ = JJ+1
      RDIST(II,JJ) = SQRT( (WXXA(1,I)-WXXA(1,J))**2
     $                   + (WXXA(2,I)-WXXA(2,J))**2
     $                   + (WXXA(3,I)-WXXA(3,J))**2 )
      RDIST(JJ,II) = RDIST(II,JJ)
      ENDDO
      ENDDO
C
      NumPt = I2-I1+1
      CALL IZeroIT(IPNT,NumPt)
c
      IPNT(1) = 1
      ILST(1) = 1
c
c  find point closest to first point
c
      DMin = Big
      DO J=2,NumPt
      Dist = RDIST(J,1)
      If(Dist.LT.DMin) Then
        JPtr = J
        DMin = Dist
      EndIf
      ENDDO
c
c  now have closest point to first point
c  mark as being found
c
      IPNT(JPtr) = 1
      ILST(2) = JPtr
c
c  now find point with sum of distances from these two points smallest
c  then continue with sum of distances from 3 points etc...
c
      numf = 2
c
 20   CONTINUE
      DMin = big
      DO K=1,NumPt
      If(IPNT(K).EQ.0) Then
       Dist = Zero
       DO L=1,numf
       LL = ILST(L)
       Dist = Dist + RDIST(LL,K)
       ENDDO
       If(Dist.LT.DMin) Then
        KPtr = K
        DMin = Dist
       EndIf
      EndIf
      ENDDO
c
      numf = numf+1
      IPNT(KPtr) = 1
      ILST(numf) = KPtr
      If(numf.LT.NumPt) GO TO 20
c
c  at this point should have sorted this shell
c  order the grid points and weights
c
      I1 = I1-1
      DO I=1,NumPt
      II = I1+ILST(I)
      JJ = I1+I
      WTA(JJ) = WWTA(II)
      XXA(1,JJ) = WXXA(1,II)
      XXA(2,JJ) = WXXA(2,II)
      XXA(3,JJ) = WXXA(3,II)
      ENDDO
ccccc
 50   CONTINUE
C
      RETURN
      END
