      subroutine rwinit
c  this routine initializes the  rw system
      IMPLICIT REAL*8 (A-H,O-Z)
      character*8 nam
      common/tape/inpx,inp2x,ioutx,ipunx,iarc,icondx,itestx,nplx(9),nbi
      COMMON /TAPE1/ INP,INP2,IOUT,IPUN,IX,ICOND,ITEST,NENTRY,LTAB,NTAP,
     1 NPL(9),NBL(9),NEN(9),LENTRY(9),NAM(2000),NUM(2000),IREC(2000),
     2 NFIL(2000)
      ltab=2000
      ntap=9
      itest=0
      inp=inpx
      inp2=inp2x
      iout=ioutx
      ipun=ipunx
      icond=icondx
c
      DO 30 I=1,NTAP
         NEN(I)=1
         NPLX(I)=10+i
         npl(i)=nplx(i)
         NBL(I)=1
         LENTRY(I)=0
   30 CONTINUE
      end
c
c
      SUBROUTINE REA (A,L,N1,NAME)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*8 nam,na,na1
      character*(*) name
      DIMENSION A(1)
      COMMON /TAPE1/ INP,INP2,IOUT,IPUN,IX,ICOND,ITEST,NENTRY,LTAB,NTAP,
     1 NPL(9),NBL(9),NEN(9),LENTRY(9),NAM(2000),NUM(2000),IREC(2000),
     2 NFIL(2000)
      N=N1
      IF (N.LT.0) N=-N1
      L2=L
      IF (N.GT.0.AND.N.LT.NTAP) GO TO 10
      CALL nerror(1,'rea','tape no. out of range',n,ntap)
   10 continue
      na=name
      NP=NPL(N)
      IF (N1.GT.0) GO TO 30
      IF (L2.GT.0) GO TO 20
      READ (NP) NA1,L1,L1
      BACKSPACE NP
      L2=L1
      IF (L.NE.L2.AND.L.NE.0) WRITE (IOUT,100) NA,N,L,L2
   20 CALL REAR (NP,NA1,A,L2)
      NBL(N)=NBL(N)+1
      GO TO 60
   30 CALL SUDAT (N,NA,NI)
      NP=NPL(N)
      IF (NI.NE.0) GO TO 40
      CALL nerror(2,'rea','data not found, code='//na,n,l)
   40 NU=NI
      CALL SPUL (N,NU)
      IF (ITEST.LE.0.AND.L.GT.0) GO TO 50
      READ (NP) NA1,L1,L2
      NBL(N)=NBL(N)+1
      CALL SPUL(N,NU)
      IF (L.LE.0) L2=L1
      IF (L.NE.L2.AND.L.NE.0)  then
        WRITE (IOUT,100) NA,N,L,L1
        WRITE (icond,100) NA,N,L,L1
      END IF
      IF (NA1.NE.NA) GO TO 70
   50 CALL REAR (NP,NA1,A,L2)
      NBL(N)=NU+1
   60 IF (ITEST.GE.1) WRITE (IOUT,110) N,L,NEN(N),NBL(N),NA
      IF (NA1.EQ.NA) RETURN
 70   continue
      CALL nerror(1,'wri','tape dir. wrong'//','//na//','//na1,np,nu)
      RETURN
C
  100 FORMAT (/1X,27HRECORD LENGTH ERROR,NAME=  ,I7,7H  UNIT=,I4,15H  LE
     1NGTH SPECD=,I7,15H  LENGTH FOUND=,I7,/)
  110 FORMAT (/1X,25HSUBROUTINE REA,UNIT NO,= ,I3,3X,11HREC.LENGTH=,I6,3
     1X,15HNUMBER OF REC.=,I4,3X,9HPOSITION=,I4,3X,5HNAME=,A8)
C
      END
c
C
      SUBROUTINE WRI (A,L,N1,K,NAME)
      IMPLICIT REAL*8 (A-H,O-Z)
c This is the main I/O write routine. See also the description
c in old TEXAS. It writes a record on file npl(n1), consisting of
c a name (8 char), a length L and an array A(1:L)
c N1 may be negative - in this case no entry is made in the
c file table but it is considered a continuation record group
c The following tables are kept:
c nam(2000): the names of the records
c num(2000): the position of the record in its file
c irec(2000) is the record grup number. It differs from num
c  by this: irec counts the record groups while num points to
c the first record of the group
c ifil(2000): the file number corresponding to record group i
c (this is the internal file number, 1 to 9)
c Note that if a record group is written, all record groups
c which are behind it on the same file are considered destroyed
c See also below
c ARGUMENTS:
c INPUT:
c A(L) a real*8 array
c L its length
c N1 is the internal file number, between 1 and 9,
c The unit numbers default to 11,12,.. 19
c K gives the position of the record:
c    K=0: write AFTER the other records on this file
c    K=positive integer, write at ABSOLUTE POSITION K
c    K=negative, rewind by K records and write there. E.g.
c    k=-1 repllaces the last record on this file
c Name (8 characters)= name under which the record is accessed
c
C     .... NENTRY IS THE HIGHEST ENTRY IN THE TABLES NAM,NUM,IREC
C     .... NBL IS THE CURRENT TAPE POSITION, NEN IS RECORD COUNT
C     .... LENTRY IS THE NUMBER OF RECORD GROUPS AND, OF COURDE, OF ENTR
C     .... IN THE TABLES FOR A PARTICULAR TAPE
C     .... NAM CONTAINS THE RECORD GROUP NAMES, NUM THE FIRST RECORD'S P
C     .... TION FOR EACH RECORD GROUP, IREC IS THE NUMBER OF THE RECORD
      character*(*) name
      character*8 nam,na,na1,bl8
      DIMENSION A(L)
      COMMON /TAPE1/ INP,INP2,IOUT,IPUN,IX,ICOND,ITEST,NENTRY,LTAB,NTAP,
     1 NPL(9),NBL(9),NEN(9),LENTRY(9),NAM(2000),NUM(2000),IREC(2000),
     2 NFIL(2000)
      data bl8/'        '/
c
      N=N1
      IF (N.LT.0) N=-N1
      IF (N.le.0.or.N.gt.NTAP) then
        CALL nerror(1,'wri','tape no. out of range',n,ntap)
      end if
      NA=NAME
      NP=NPL(N)
      M=K
      IF (K.LE.0) M=LENTRY(N)+1
      IF (K.LT.0) M=M+K
      NM=NEN(N)
      IFREE=0
      IF (N1.LT.0) GO TO 50
      DO 30 I=1,NENTRY
         II=nfil(I)
         IF (II.NE.N) GO TO 20
         IF (IREC(I).EQ.M) NM=NUM(I)
         IF (IREC(I).ge.M) then
c remove record groups preceding the current one
           IREC(I)=0
           NAM(I)=bl8
           NUM(I)=0
         end if
   20    IF (IFREE.EQ.0.AND.NAM(I).EQ.bl8) IFREE=I
         IF (NAM(I).eq.NA) then
           NAM(I)(8:8)='X'
cc           WRITE (IOUT,110) NA,N,NUM(I),nam(i)
cc           WRITE (ICOND,110) NA,N,NUM(I),nam(i)
  110 FORMAT (1X,'DUPLICATE DATA,OLD ONE CHANGED TO X IN POS.8 ' ,A8,
     1 5HFILE=,I3,' RECORD=',I4,' NEW NAME=',A8)
         end if
   30 CONTINUE
   40 CONTINUE
      IF (IFREE.EQ.0) IFREE=NENTRY+1
      CALL SPUL (N,NM)
      WRITE (NP) NA,L,L,A
      LENTRY(N)=M
      NBL(N)=NM+1
      NEN(N)=NM+1
      GO TO 60
   50 CALL SPUL (N,NM)
      WRITE (NP) bl8,L,L,A
      NBL(N)=NBL(N)+1
      NEN(N)=NEN(N)+1
C
C     .... CONTINUATION OF A RECORD GROUP -NO TABLE ENTRY IS MADE
C
   60 IF (ITEST.LT.1) GO TO 70
      WRITE (IOUT,120) N,L,NEN(N),NBL(N),NA
   70 IF (IFREE.LE.LTAB) GO TO 80
      CALL nerror(2,'wri','table overflow, wanted & max.',ifree,ltab)
   80 IF (IFREE.EQ.0) GO TO 90
      NAM(IFREE)=NA
      NUM(IFREE)=NM
      IREC(IFREE)=M
      NFIL(IFREE)=N
      IF (IFREE.GT.NENTRY) NENTRY=IFREE
   90 RETURN
C
  100 FORMAT (/1X,21HTAPE NO. OUT OF RANGE,I5)
  120 FORMAT (/1X,25HSUBROUTINE WRI,UNIT NO,= ,I3,3X,11HREC.LENGTH=,I6,3
     1X,15HNUMBER OF REC.=,I4,3X,9HPOSITION=,I4,3X,5HNAME=,A8)
  130 FORMAT (/1X,33HTOO MANY RECORDS, TABLE OVERFLOW ,/)
C
      END
C
      SUBROUTINE REAR (NP,NA1,A,L2)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*8 na1
      DIMENSION A(L2)
      READ (NP) NA1,L,L,A
      RETURN
      END
c
      SUBROUTINE SPUL (N,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*8 nam
      COMMON /TAPE1/ INP,INP2,IOUT,IPUN,IX,ICOND,ITEST,NENTRY,LTAB,NTAP,
     1 NPL(9),NBL(9),NEN(9),LENTRY(9),NAM(2000),NUM(2000),IREC(2000),
     2 NFIL(2000)
      NP=NPL(N)
      M=L-NBL(N)
C     IF (-M.LE.L/2) GO TO 10
      IF(M.GE.0) GO TO 10
      REWIND NP
      M=L-1
   10 IF (M) 20,60,40
   20 M=-M
      DO 30 I=1,M
   30 BACKSPACE NP
      GO TO 60
   40 DO 50 I=1,M
   50 READ (NP)
C
C      READ (NP) IS PROBABLY BETTER - JUST THAT IT DOES NOT WORK ON THE
C      HARRIS MINICOMPUTER - IF POSSIBLE SUBSTITUTE READ(NP)
C
   60 NBL(N)=L
      RETURN
C
      END
c=================================================================
      SUBROUTINE SUDAT (N,NA,NI)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*8 name,na,nam
c sudat returns the block number ni (NOT the position in the nam,num,irec arrays!)
c for name na on unit n
      COMMON /TAPE1/ INP,INP2,IOUT,IPUN,IX,ICOND,ITEST,NENTRY,LTAB,NTAP,
     1 NPL(9),NBL(9),NEN(9),LENTRY(9),NAM(2000),NUM(2000),IREC(2000),
     2 NFIL(2000)
CPP
      dimension ic1(8),ic2(8)
      NI=0
      If (NENTRY.EQ.0) GO TO 30
      DO 10 I=1,NENTRY
c both name and file No. are ok
         IF (na.EQ.nam(i).AND.N.EQ.nfil(I)) then
           NI=I
         end if
   10 CONTINUE
      IF (NI.NE.0) then
        ni=num(ni)
        RETURN
      end if
      DO 20 I=1,NENTRY
c file No. is wrong but data are still on a file
         IF (NA.EQ.NAM(I)) then
           n3=nfil(i)
           NI=I
         end if
   20 CONTINUE
      IF (NI.EQ.0) GO TO 30
      WRITE (IOUT,40) N,NA,N3,nam(i)
      write (icond,40) n,na,n3,nam(i)
   30 IF (ITEST.EQ.1.AND.NI.EQ.0) WRITE (IOUT,50) N,NA
      if(ni.eq.0) return
      ni=num(ni)
      n=n3
      RETURN
C
   40 FORMAT (/,'DATA ON THE WRONG UNIT- SOUGHT',I3,1x,a8,
     1 'FOUND',I3,1x,a8,' PROGRAM PROCEEDS',/)
   50 FORMAT (1X,19HNO SUCH DATA EXIST ,I5,2x,a8)
C
      END
c======================================================================
      subroutine restart(nlogi,nrec)
c  This routine reconstructs the internal record table for
c  INTERNAL logical unit nlogi. Usually, the unit number is
c  nlogi+10, e.g. for INTERNAL unit 1, it is 11
c  nrec is the maximum record number until restore goes - if zero,
c  the whole file is restored.
      implicit real*8 (a-h,o-z)
c      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
c     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
c     20)
      character*8 nam,na1,na,bl8
      common/tape/inpx,inp2x,ioutx,ipunx,iarc,icondx,itestx,nplx(9),nbi
      COMMON /TAPE1/ INP,INP2,IOUT,IPUN,IX,ICOND,ITEST,NENTRY,LTAB,NTAP,
     1 NPL(9),NBL(9),NEN(9),LENTRY(9),NAM(2000),NUM(2000),IREC(2000),
     2 NFIL(2000)
      data bl8/'        '/
         if (nlogi.le.0.or.nlogi.gt.ntap) go to 50
         np=npl(nlogi)
         rewind np
         nbl(nlogi)=1
         nen(nlogi)=1
         nblo=0
         na=bl8
         lentry(nlogi)=0
   10    read (np,end=40) na1
   20    continue
         nblo=nblo+1
         nbl(nlogi)=nbl(nlogi)+1
         nen(nlogi)=nen(nlogi)+1
         if (nblo.eq.nrec) go to 50
         if (na1.ne.bl8) go to 30
c   *** the distinction between na1 and na is only important
c       for record groups or file 2 but this is not restarted any more
         go to 10
 30      call sudat(nlogi,na,ni)
c    *** check if this name occurs earlier
         if (ni.gt.0.and.nam(ni).eq.na1) then
           nam(ni)(8:8)='X'
         endif
         nentry=nentry+1
         lentry(nlogi)=lentry(nlogi)+1
         irec(nentry)=lentry(nlogi)
         nam(nentry)=na1
         num(nentry)=nbl(nlogi)-1
         nfil(nentry)=nlogi
c
         na=na1
         go to 10
   40    nbl(nlogi)=nbl(nlogi)+1
   50 continue
c
      end
