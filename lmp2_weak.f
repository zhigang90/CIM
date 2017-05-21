C================pcksetup===========================================
      subroutine pcksetup(nel,inx,ncore,lbas,idm,
     1                    iprint)

      use memory

        implicit real*8(a-h,o-z)
c       common/big/bl(30000)
        dimension inx(12,*)
        integer*2 lbas(*),idm(*)
C
C        setup for pick taken form scf1.f from texas.
C        modified for PQS by Svein Saebo, September 2000
C        everything in CAPS from original Boughton Pulay program
C
        na=igetival('na')
        ncf=igetival('ncf')
        ncs=igetival('ncs')
        nval=nel-ncore
C--------- Obtain LMAXA=Maximum number of functions on 4 atoms --------C
C          (This is mindlessly crude, a better way must be found.)
      NNN1=0
      NNN2=0
      NNN3=0
      NNN4=0
      DO 255 IATOM=1,NA
         NNNA=0
         DO 251 ISHELL=1,NCS
            NUMAT=INX(2,ISHELL)
C           IF(NUMAT.EQ.IATOM) NNNA=NNNA + INX(3,ISHELL)
C        ss fix for general contraction
        MMM =inx(3,ishell)*(inx(4,ishell)+1)
            IF(NUMAT.EQ.IATOM) NNNA=NNNA + MMM
  251    CONTINUE
         IF(NNNA.GT.NNN1) THEN
           NNN1=NNNA
           IA1=IATOM
         ENDIF
  255 CONTINUE
      DO 265 IATOM=1,NA
         IF(IATOM.EQ.IA1) GOTO 265
         NNNA=0
         DO 261 ISHELL=1,NCS
            NUMAT=INX(2,ISHELL)
        MMM =inx(3,ishell)*(inx(4,ishell)+1)
            IF(NUMAT.EQ.IATOM) NNNA=NNNA + mmm
  261    CONTINUE
         IF(NNNA.GT.NNN2) THEN
           NNN2=NNNA
           IA2=IATOM
         ENDIF
  265 CONTINUE
      DO 275 IATOM=1,NA
         IF(IATOM.EQ.IA1.OR.IATOM.EQ.IA2) GOTO 275
         NNNA=0
         DO 271 ISHELL=1,NCS
            NUMAT=INX(2,ISHELL)
        MMM =inx(3,ishell)*(inx(4,ishell)+1)
            IF(NUMAT.EQ.IATOM) NNNA=NNNA + mmm
  271    CONTINUE
         IF(NNNA.GT.NNN3) THEN
            NNN3=NNNA
            IA3=IATOM
         ENDIF
  275 CONTINUE
      DO 285 IATOM=1,NA
         IF(IATOM.EQ.IA1.OR.IATOM.EQ.IA2.OR.IATOM.EQ.IA3) GOTO  285
         NNNA=0
         DO 281 ISHELL=1,NCS
            NUMAT=INX(2,ISHELL)
        MMM =inx(3,ishell)*(inx(4,ishell)+1)
            IF(NUMAT.EQ.IATOM) NNNA=NNNA + mmm
  281    CONTINUE
         IF(NNNA.GT.NNN4) NNN4=NNNA
  285 CONTINUE
      LMAXA= NNN1 + NNN2 + NNN3 + NNN4
C----------------------------------------------------------------------C
        call mmark
        call getmem(ncf**2,is2)
        call getmem(ncf*nval,ic)
         CALL getmem(LMAXA**2,I5)
         CALL getmem(LMAXA,I6)
         CALL getmem(LMAXA,I7)
         CALL getmem(LMAXA,I8)
         CALL getmem(LMAXA,I9)
         CALL getmem(NCF,I10)
         CALL getint(4*NEL,I11)
         CALL getmem(4*NEL,I12)
         CALL getint(NCF,I13)
         CALL getint(NA,I14)
         CALL getmem(LMAXA,I15)
        np1=1
        np4=4
        iloca=mataddr('occu')
        call matdef('ovls','s',ncf,ncf)
        call matread('ovls',np1,'s matrix')
        iovla=mataddr('ovls')
C        call matprint('loca',6)
      call pick(bl(ic),bl(iloca),bl(iovla),bl(is2),bl(i5),
     1          bl(i6),bl(i7),bl(i8),bl(i9),bl(i10),
     2          BL(I11),BL(I12),bL(I13),BL(I14),NCF,
     3          NEL,LMAXA,NA,np1,NP4,
     4          BL(I15),inx,ncore,nval,lbas,
     5          idm,iprint)
C----------------------------------------------------------------------C
        call retmark
        end
C=========pick=========================================
      SUBROUTINE PICK (C,wf,ovl,S2,SA,
     1                 U,V,W,A,CC,
     2                 IPCK,PCK,MAP,MORPH,NAO,
     3                 NMO,LMAXA,NATOMS,NP1,NP4,
     4                 Z,inx,ncore,nval,lbas,
     5                 idm,iprint)
C----------------------------------------------------------------------C
C        This is the  Boughton, Pulay pick routine modified for
C        PQS by Svein Saebo, September 2000
C        everything in CAPS from original Boughton, Pulay program
C        comments are also from original subroutine
C----------------------------------------------------------------------C
C     Pick the atomic centers for each localized molecular orbital to  C
C   be used in the localized electron correlation treatment.  The      C
C   overlap and wavefunction matrices are mapped to atom-contiguous    C
C   form.  They are therefore incompatible with the rest of the        C
C   molecular calculation and should be used only for the purposes of  C
C   this routine.                                                      C
C     The arrays IPCK and PCK will end up with entries as follows;     C
C     IPCK(1,*)= Number of the first-order atom                        C
C     IPCK(2,*)= Number of the second-order atom                       C
C     IPCK(3,*)= Number of the third-order atom                        C
C     IPCK(4,*)= Number of the fourth-order atom                       C
C      PCK(1,*)= Least-squares residual if first-order atom is used    C
C      PCK(2,*)= Least-squares residual if second-order atom also used C
C      PCK(3,*)= Least-squares residual if three atoms are used        C
C      PCK(4,*)= Least-squares residual if all four atoms are used     C
C                                                                      C
C     This version of PICK truncates the picking process at two atoms  C
C   (for now at least).  All unique permutations of atoms taken two at C
C   a time ARE NOT EXAMINED.  Experience indicates that the same pair  C
C   of atoms is picked if the process commits to the first atom and    C
C   then searches for a complimentary atom to complete the pair.  This C
C   is certainly quicker and simpler, but the nonrigorousness must be  C
C   recognized.  A more satisfactory method may perhaps later be found.C
C----------------------------------------------------------------------C
        implicit real*8(a-h,o-z)
C     DIMENSION Q(5,NATOMS)
      DIMENSION IPCK(4,nval),PCK(4,nval),CC(NAO),MAP(NAO),MORPH(NATOMS)
      DIMENSION C(NAO,nval),S2(NAO,NAO),Z(LMAXA)
      DIMENSION A(LMAXA),U(LMAXA),V(LMAXA),W(LMAXA),SA(LMAXA*LMAXA)
        dimension ovl(*),wf(nao,nval)
        integer*2 lbas(nao,*),idm(*),en
      dimension inx(12,*)
        parameter(zero=0.0d0,one=1.0d0)
C
        ncs=igetival('ncs')
      FITCRIT=0.02D0
C------------------------- CREATE MAPPING VECTOR ----------------------C
        en=1
        na=natoms
      LAG=0
      DO 200 IATOM=1,NA
      NDX=0
      DO 150 ISHELL=1,NCS
         NAT=INX(2,ISHELL)
        NMM =inx(3,ishell)*(inx(4,ishell)+1)
C        NMM=INX(3,ISHELL)
         IF(NAT.NE.IATOM) NDX=NDX + NMM
         IF(NAT.EQ.IATOM) THEN
            DO 100 I=1,NMM
               NDX=NDX+1
               LAG=LAG+1
               MAP(NDX)=LAG
  100       CONTINUE
         ENDIF
  150 CONTINUE
      MORPH(IATOM)=LAG
  200 CONTINUE
C-------------- MAP WAVEFUNCTION TO ATOM-CONTIGUOUS FORM --------------C
      DO 300 I=1,NAO
         IM=MAP(I)
         DO 250 J=1,nval
            C(IM,J)=wf(I,J)
  250    CONTINUE
  300 CONTINUE
C------------- MAP OVERLAP MATRIX TO ATOM-CONTIGUOUS FORM -------------C
      K=0
      DO 400 J=1,NAO
         JM=MAP(J)
         DO 350 I=1,J
            IM=MAP(I)
            K=K+1
            S2(IM,JM)=ovl(k)
            S2(JM,IM)=ovl(k)
  350    CONTINUE
  400 CONTINUE
C----------------------- BEGIN PICKING ALGORITHM ----------------------C
      DO 2000 MO=1,nval
         PCK(1,MO)=one
         PCK(2,MO)=zero
         PCK(3,MO)=zero
         PCK(4,MO)=zero
         IPCK(1,MO)=0
         IPCK(2,MO)=0
         IPCK(3,MO)=0
         IPCK(4,MO)=0
C     ----------- FIND FIRST-ORDER ATOMIC CENTER FOR MO ----------
         CALL DCOPY(NAO,C(1,MO),1,CC(1),1)
         IB=1
         DO 600 IAT=1,NA
            IE=MORPH(IAT)
            L = IE-IB + 1
            CALL DSCAL(L,ZERO,V,1)
            DO 450 J=1,NAO
               CALL DAXPY(L,CC(J),S2(IB,J),1,V,1)
  450       CONTINUE
            CALL DCOPY(L,V,1,Z,1)
            DO 500 J=1,L
               CALL DCOPY(L,S2(IB,IB+J-1),1,SA((J-1)*L+1),1)
  500       CONTINUE
            CALL HOUSEFIT(SA,L,L,V,L,1,U,W)
            RESID=1.0D0 - DDOT(L,U,1,Z,1)
            IF(RESID.LT.PCK(1,MO)) THEN
               PCK(1,MO)=RESID
               IPCK(1,MO)=IAT
               IB1=IB
               LE1=L
            ENDIF
  550        IB=IE+1
  600    CONTINUE
         IF(PCK(1,MO).LT.FITCRIT) GOTO 2000
C     ---------- FIND SECOND-ORDER ATOMIC CENTER FOR MO ----------
         PCK(2,MO)=PCK(1,MO)
         IB=1
         DO 800 IAT=1,NA
            IE=MORPH(IAT)
            IF(IAT.EQ.IPCK(1,MO)) GOTO 750
            L = IE-IB + 1
            LL=LE1 + L
            CALL DSCAL(LL,ZERO,V,1)
            DO 650 J=1,NAO
               CALL DAXPY(LE1,CC(J),S2(IB1,J),1,V,1)
               CALL DAXPY(L,CC(J),S2(IB,J),1,V(LE1+1),1)
  650       CONTINUE
            CALL DCOPY(LL,V,1,Z,1)
            DO 700 J=1,LE1
               CALL DCOPY(LE1,S2(IB1,IB1+J-1),1,SA((J-1)*LL+1),1)
               CALL DCOPY(L,S2(IB,IB1+J-1),1,SA((J-1)*LL+1+LE1),1)
  700       CONTINUE
            LLL=LE1*LL
            DO 720 J=1,L
               CALL DCOPY(LE1,S2(IB1,IB+J-1),1,SA((J-1)*LL+1+LLL),1)
               CALL DCOPY(L,S2(IB,IB+J-1),1,SA((J-1)*LL+1+LE1+LLL),1)
  720       CONTINUE
            CALL HOUSEFIT(SA,LL,LL,V,LL,1,U,W)
            RESID=1.0D0 - DDOT(LL,U,1,Z,1)
            IF(RESID.LT.PCK(2,MO)) THEN
               PCK(2,MO)=RESID
               IPCK(2,MO)=IAT
               IB2=IB
               LE2=L
            ENDIF
  750       IB=IE+1
  800    CONTINUE
         IF(PCK(2,MO).LT.FITCRIT) GOTO 2000
C     ---------------- TRUNCATION OF THIS SUBROUTINE -------------
C         GOTO 2000
C     ------------------------------------------------------------
C     ----------- FIND THIRD-ORDER ATOMIC CENTER FOR MO ----------
         PCK(3,MO)=PCK(2,MO)
         IB=1
         SSS=0.0
         DO 1000 IAT=1,NA
            IE=MORPH(IAT)
            IF(IAT.EQ.IPCK(1,MO).OR.IAT.EQ.IPCK(2,MO)) GOTO 950
            SUM=0.0
            DO 820 J=IB,IE
              DO 810 I=IB,IE
                 SUM=SUM + CC(I)*CC(J)*S2(I,J)
  810         CONTINUE
  820       CONTINUE
            IF(SUM.LT.SSS) GOTO 950
            SSS=SUM
            L = IE-IB + 1
            LL= LE1 + LE2 + L
            CALL DSCAL(LL,ZERO,V,1)
            LELE=LE1 + LE2 + 1
            DO 850 J=1,NAO
               CALL DAXPY(LE1,CC(J),S2(IB1,J),1,V,1)
               CALL DAXPY(LE2,CC(J),S2(IB2,J),1,V(LE1+1),1)
               CALL DAXPY(L,CC(J),S2(IB,J),1,V(LELE),1)
  850       CONTINUE
            CALL DCOPY(LL,V,1,Z,1)
            DO 900 J=1,LE1
               CALL DCOPY(LE1,S2(IB1,IB1+J-1),1,SA((J-1)*LL+1),1)
               CALL DCOPY(LE2,S2(IB2,IB1+J-1),1,SA((J-1)*LL+1+LE1),1)
               CALL DCOPY(L,S2(IB,IB1+J-1),1,SA((J-1)*LL+LELE),1)
  900       CONTINUE
            LLL=LE1*LL
            DO 920 J=1,LE2
               CALL DCOPY(LE1,S2(IB1,IB2+J-1),1,SA((J-1)*LL+1+LLL),1)
             CALL DCOPY(LE2,S2(IB2,IB2+J-1),1,SA((J-1)*LL+1+LE1+LLL),1)
               CALL DCOPY(L,S2(IB,IB2+J-1),1,SA((J-1)*LL+LLL+LELE),1)
  920       CONTINUE
            LLL=LLL+LE2*LL
            DO 930 J=1,L
               CALL DCOPY(LE1,S2(IB1,IB+J-1),1,SA((J-1)*LL+1+LLL),1)
               CALL DCOPY(LE2,S2(IB2,IB+J-1),1,SA((J-1)*LL+1+LE1+LLL),1)
               CALL DCOPY(L,S2(IB,IB+J-1),1,SA((J-1)*LL+LLL+LELE),1)
  930       CONTINUE
            CALL HOUSEFIT(SA,LL,LL,V,LL,1,U,W)
            RESID=1.0D0 - DDOT(LL,U,1,Z,1)
C            IF(RESID.LT.PCK(3,MO)) THEN
               PCK(3,MO)=RESID
               IPCK(3,MO)=IAT
               IB3=IB
               LE3=L
C            ENDIF
  950       IB=IE+1
 1000    CONTINUE
         IF(PCK(3,MO).LT.FITCRIT) GOTO 2000
C     ----------- FIND FOURTH-ORDER ATOMIC CENTER FOR MO ----------
         PCK(4,MO)=PCK(3,MO)
         IB=1
         SSS=0.0
         DO 1400 IAT=1,NA
            IE=MORPH(IAT)
            IF(IAT.EQ.IPCK(1,MO).OR.IAT.EQ.IPCK(2,MO)) GOTO 1350
            IF(IAT.EQ.IPCK(3,MO)) GOTO 1350
            SUM=0.0
            DO 1020 J=IB,IE
              DO 1010 I=IB,IE
                 SUM=SUM + CC(I)*CC(J)*S2(I,J)
 1010         CONTINUE
 1020       CONTINUE
            IF(SUM.LT.SSS) GOTO 1350
            SSS=SUM
            L = IE-IB + 1
            LL= LE1 + LE2 + LE3 + L
            CALL DSCAL(LL,ZERO,V,1)
            LELE=LE1 + LE2 + 1
            LEL3=LELE + LE3
            DO 1100 J=1,NAO
               CALL DAXPY(LE1,CC(J),S2(IB1,J),1,V,1)
               CALL DAXPY(LE2,CC(J),S2(IB2,J),1,V(LE1+1),1)
               CALL DAXPY(LE3,CC(J),S2(IB3,J),1,V(LELE),1)
               CALL DAXPY(L,CC(J),S2(IB,J),1,V(LEL3),1)
 1100       CONTINUE
            CALL DCOPY(LL,V,1,Z,1)
            DO 1150 J=1,LE1
               CALL DCOPY(LE1,S2(IB1,IB1+J-1),1,SA((J-1)*LL+1),1)
               CALL DCOPY(LE2,S2(IB2,IB1+J-1),1,SA((J-1)*LL+1+LE1),1)
               CALL DCOPY(LE3,S2(IB3,IB1+J-1),1,SA((J-1)*LL+LELE),1)
               CALL DCOPY(L,S2(IB,IB1+J-1),1,SA((J-1)*LL+LEL3),1)
 1150       CONTINUE
            LLL=LE1*LL
            DO 1200 J=1,LE2
               CALL DCOPY(LE1,S2(IB1,IB2+J-1),1,SA((J-1)*LL+1+LLL),1)
             CALL DCOPY(LE2,S2(IB2,IB2+J-1),1,SA((J-1)*LL+1+LE1+LLL),1)
               CALL DCOPY(LE3,S2(IB3,IB2+J-1),1,SA((J-1)*LL+LLL+LELE),1)
               CALL DCOPY(L,S2(IB,IB2+J-1),1,SA((J-1)*LL+LLL+LEL3),1)
 1200       CONTINUE
            LLL=LLL+LE2*LL
            DO 1250 J=1,LE3
               CALL DCOPY(LE1,S2(IB1,IB3+J-1),1,SA((J-1)*LL+1+LLL),1)
             CALL DCOPY(LE2,S2(IB2,IB3+J-1),1,SA((J-1)*LL+1+LE1+LLL),1)
               CALL DCOPY(LE3,S2(IB3,IB3+J-1),1,SA((J-1)*LL+LLL+LELE),1)
               CALL DCOPY(L,S2(IB,IB3+J-1),1,SA((J-1)*LL+LLL+LEL3),1)
 1250       CONTINUE
            LLL=LLL+LE3*LL
            DO 1300 J=1,L
               CALL DCOPY(LE1,S2(IB1,IB+J-1),1,SA((J-1)*LL+1+LLL),1)
               CALL DCOPY(LE2,S2(IB2,IB+J-1),1,SA((J-1)*LL+1+LE1+LLL),1)
               CALL DCOPY(LE3,S2(IB3,IB+J-1),1,SA((J-1)*LL+LLL+LELE),1)
               CALL DCOPY(L,S2(IB,IB+J-1),1,SA((J-1)*LL+LLL+LEL3),1)
 1300       CONTINUE
            CALL HOUSEFIT(SA,LL,LL,V,LL,1,U,W)
            RESID=1.0D0 - DDOT(LL,U,1,Z,1)
C            IF(RESID.LT.PCK(4,MO)) THEN
               PCK(4,MO)=RESID
               IPCK(4,MO)=IAT
C            ENDIF
 1350       IB=IE+1
 1400    CONTINUE
 2000 CONTINUE
C------------------------ END PICKING ALGORITHM -----------------------C
C---------------------------- PRINT RESULTS ---------------------------C
      if(iprint.ge.2)write(6,*) 'Local domains for weak pairs:'
      DO 3000 M=1,nval
         I1=IPCK(1,M)
         I2=IPCK(2,M)
         I3=IPCK(3,M)
         I4=IPCK(4,M)
        idd=0
        do iao=1,nao
        lbas(iao,m) =0
        enddo
      do ics=1,ncs
            NAT=INX(2,ics)
            IF(NAT.EQ.I1.OR.NAT.EQ.I2.OR.NAT.EQ.I3.OR.NAT.EQ.I4) THEN
        ifirst=inx(11,ics)+1
        ilast=inx(10,ics)
        do icf =ifirst,ilast
        lbas(icf,m)=en
        idd=idd+1
        enddo
        endif
        enddo
        idm(m)=idd
        iao=0
        do ix=1,nao
        if(lbas(ix,m).eq.en) then
        iao=iao+1
        lbas(iao,m)=ix
        endif
        enddo
 3000 CONTINUE
        iddmax=0
        do ii=1,nval
        idd=idm(ii)
        iddmax=max0(iddmax,idd)
        if(iprint.ge.2) then
      write(6,*) 'Correlated MO no: ',ii,
     1' Atoms contributing: ',(ipck(iw,ii),iw=1,4),
     2' dimension', idm(ii)
      if(iprint.ge.3)write(6,*) (lbas(iw,ii),iw=1,idd)
        endif
        end do
        if(iprint.ge.2)write(6,*)'Max domain size: ',iddmax
      END
C
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE HOUSEFIT(A,M,N,B,K,L,X,W)
C======================================================================C
C     Solve a system of simultaneous equations, AX=B. If B is a matrix C
C   then X, the  matrix of multiple solutions,  must  have  the same   C
C   dimensions as B.                                                   C
C     The system is allowed to be  overdetermined, in  which case you  C
C   will get the least-squares solution(s).                            C
C======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), B(K,L), X(K,L), W(M)
      CALL HOUSETRANS(A,M,N,B,K,L,W)
      CALL BACKSUB(A,M,N,B,K,L,X)
      RETURN
      END
C
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE HOUSETRANS(A,M,N,B,KB,LB,W)
C======================================================================C
C                 Householder transformation.                          C
C======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(M,N) , B(KB,LB) , W(M)
      DO 300 K=1,N
         MX=IDAMAX(M-K+1,A(K,K),1) + K-1
         RMS=0.0
         DO 50 I=K,M
            W(I)=A(I,K)/DABS(A(MX,K))
            RMS=RMS+W(I)*W(I)
   50    CONTINUE
         RMS=DSQRT(RMS)
         BK=1/(RMS*(RMS + DABS(W(K))))
         ZZ=W(K)
         W(K)=W(K) + DSIGN(RMS,ZZ)
         DO 100 J=1,N
            S=DDOT(M-K+1,W(K),1,A(K,J),1)
            S=BK*S
            CALL DAXPY(M-K+1,-S,W(K),1,A(K,J),1)
  100    CONTINUE
         DO 200 J=1,LB
            S=DDOT(M-K+1,W(K),1,B(K,J),1)
            S=BK*S
            CALL DAXPY(M-K+1,-S,W(K),1,B(K,J),1)
  200    CONTINUE
  300 CONTINUE
      RETURN
      END
C
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE BACKSUB(A,M,N,B,KB,LB,X)
C======================================================================C
C                       Execute backsubstitution.                      C
C======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N) , B(KB,LB) , X(KB,LB)
      CALL DCOPY(KB*LB,B,1,X,1)
      DO 300 L=1,LB
         DO 200 J=N,2,-1
            X(J,L)=X(J,L)/A(J,J)
            CALL DAXPY(J-1,-X(J,L),A(1,J),1,X(1,L),1)
  200    CONTINUE
         X(1,L)=X(1,L)/A(1,1)
  300 CONTINUE
      RETURN
      END
C
C===========Weakpair===================================
      subroutine Weakpair(ncf,nval,lbasw,idmw,ncore,
     1                    iprint,nweak,nkunix,listp,mstrong,
     2                    uniq,lfil,wrdw,giga)
C        main routine for calculation of weak pairs using
C        multipole expansion
C        dipoles and quadrapoles only for now
C        the weak pair energy is returned in eweak2
C        called from mp2iter in lmp2_iter.f
C
C        Svein Saebo, Fayetteville, AR Summer and Fall 2000
C
C        lbas,idm,lbasij,idmij are the original local bases used
C        for the strong pairs
C        lbasw,idmw are (small) local domains used for weak pairs
C

      use memory

        implicit real*8(a-h,o-z)
        integer*2 lbasw(ncf,*),idmw(*)
        dimension listp(*)
        integer uniq(*)
c       common/big/bl(30000)
        parameter(zero=0.0d0,one=1.0d0)
C
        call secund(tst1)
        dtmax=0.0d0
        np4=4
        np1=1
        call mmark
C
        call matdef('rinv','s',nval,nval)
        call matread('rinv',np4,'Rijinv')
        irinv=mataddr('rinv')-1
C
C        YYx,YYy,YYZ the x,y and z components of matrix Y(ncf,nval)
C        used for the dipole part, should be saved
C
        call matdef('YYx','r',ncf,nval)
        call matdef('YYy','r',ncf,nval)
        call matdef('YYz','r',ncf,nval)
C
C        keep matrices above
C
        call matdef('prdx','q',ncf,ncf)
        call matdef('prdy','q',ncf,ncf)
        call matdef('prdz','q',ncf,ncf)
        call matdef('dipx','s',ncf,ncf)
        call matdef('dipy','s',ncf,ncf)
        call matdef('dipz','s',ncf,ncf)
C        dipx.. dipole integrals , proj projection matrix
        call matread('dipx',np1,'aoX')
        call matread('dipy',np1,'aoY')
        call matread('dipz',np1,'aoZ')
        call matmmult('dipx','proj','prdx')
        call matmmult('dipy','proj','prdy')
        call matmmult('dipz','proj','prdz')
        call matrem('dipz')
        call matrem('dipy')
        call matrem('dipx')
         call matpose('prdx')
         call matpose('prdy')
         call matpose('prdz')
        call matmmult('prdx','occu','YYx')
        call matmmult('prdy','occu','YYy')
        call matmmult('prdz','occu','YYz')
        call matrem('prdz')
        call matrem('prdy')
        call matrem('prdx')
C
        iyyx=mataddr('YYx')
        iyyy=mataddr('YYy')
        iyyz=mataddr('YYz')
C        
C        Xxx,Xxy..... are the 6 components of matrix X used
C        for the qudrapole part, should be kept
C
        call matdef('Xxx','r',ncf,nval)
        call matdef('Xxy','r',ncf,nval)
        call matdef('Xxz','r',ncf,nval)
        call matdef('Xyy','r',ncf,nval)
        call matdef('Xyz','r',ncf,nval)
        call matdef('Xzz','r',ncf,nval)
C
        call matdef('Rijx','s',nval,nval)
        call matdef('Rijy','s',nval,nval)
        call matdef('Rijz','s',nval,nval)
C        distance matrices between localized orbitals
        irijx=mataddr('Rijx')-1
        irijy=mataddr('Rijy')-1
        irijz=mataddr('Rijz')-1
        call matread('Rijx',np4,'Rijx')
        call matread('Rijy',np4,'Rijy')
        call matread('Rijz',np4,'Rijz')
C        radiii of centroids
        call matdef('rx2','d',nval,nval)
        call matdef('ry2','d',nval,nval)
        call matdef('rz2','d',nval,nval)
C
        ix2a=mataddr('rx2')-1
        iy2a=mataddr('ry2')-1
        iz2a=mataddr('rz2')-1
        call matread('rx2',np4,'r2x')
        call matread('ry2',np4,'r2y')
        call matread('rz2',np4,'r2z')
C        keep the Xxx matrices....and Rx,Ry,Rz
        call matdef('pmxx','q',ncf,ncf)
        call matdef('pmxy','q',ncf,ncf)
        call matdef('pmxz','q',ncf,ncf)
        call matdef('pmyy','q',ncf,ncf)
        call matdef('pmyz','q',ncf,ncf)
        call matdef('pmzz','q',ncf,ncf)
C        smxx for second moments integrals calculated in quadrap
        call matdef('smxx','s',ncf,ncf)
        call matdef('smxy','s',ncf,ncf)
        call matdef('smxz','s',ncf,ncf)
        call matdef('smyy','s',ncf,ncf)
        call matdef('smyz','s',ncf,ncf)
        call matdef('smzz','s',ncf,ncf)
C        matrices for second moments:
         call quadrap
        call matmmult('smxx','proj','pmxx')
        call matmmult('smxy','proj','pmxy')
        call matmmult('smxz','proj','pmxz')
        call matmmult('smyy','proj','pmyy')
        call matmmult('smyz','proj','pmyz')
        call matmmult('smzz','proj','pmzz')
        call matrem('smzz')
        call matrem('smyz')
        call matrem('smyy')
        call matrem('smxz')
        call matrem('smxy')
        call matrem('smxx')
        call matpose('pmxx')
        call matpose('pmxy')
        call matpose('pmxz')
        call matpose('pmyy')
        call matpose('pmyz')
        call matpose('pmzz')
        call matmmult('pmxx','occu','Xxx')
        call matmmult('pmxy','occu','Xxy')
        call matmmult('pmxz','occu','Xxz')
        call matmmult('pmyy','occu','Xyy')
        call matmmult('pmyz','occu','Xyz')
        call matmmult('pmzz','occu','Xzz')
        ixxa=mataddr('Xxx')
        ixya=mataddr('Xxy')
        ixza=mataddr('Xxz')
        iyya=mataddr('Xyy')
        iyza=mataddr('Xyz')
        izza=mataddr('Xzz')
        call matrem('pmzz')
        call matrem('pmyz')
        call matrem('pmyy')
        call matrem('pmxz')
        call matrem('pmxy')
        call matrem('pmxx')
C        loop over weak pairs (flagged with zero)
        call secund(tst2)
        tst=tst2-tst1
        write(6,*) 'weak startup', tst/60.d0
        rmin=100.0d0
        rmax=zero
        rsum=zero
        ncoun=0
        ten=zero
        npar=0
        nwek=0
        do ii=1,nval
        x1=bl(ix2a+ii)
        y1=bl(iy2a+ii)
        z1=bl(iz2a+ii)
        idd=idmw(ii)
        do jj=1,ii
        npar=npar+1
        if(listp(npar).ne.0) cycle
        nwek=nwek+1
        lpar=mstrong+nwek
         iuniq=uniq(lpar)
         if(iuniq.le.0) cycle
        jdd=idmw(jj)
        x2=bl(ix2a+jj)
        y2=bl(iy2a+jj)
        z2=bl(iz2a+jj)
        Rinv=bl(irinv+npar)
        Rij3=Rinv*Rinv*Rinv
        Rij4=Rij3*Rinv
        Rij=one/Rinv
C        for statistics
        rmin=min(rmin,rij)
        rmax=max(rmax,rij)
        rsum=rsum+rij
        ncoun=ncoun+1
C
        Rx=bl(irijx+npar)
        Ry=bl(irijy+npar)
        Rz=bl(irijz+npar)
C     write(6,*)'rinv= ',Rinv,'Rij= ',Rij,RIJ3,Rij4
C     write(6,*)'rx,ry,rz'
C        write(6,66) rx,ry,rz
C        write(6,66) x1,y1,z1
C        write(6,66) x2,y2,z2
C  66 format(1x,3f10.7)
C        everything needed  above
C
        call matdef('kweak','r',idd,jdd)
        kaddr=mataddr('kweak')
C
C        mkkweak contructs Kij for the weak pair ij, result in 'kweak'
C
        call mkkweak(bl(kaddr),bl(iyyx),bl(iyyy),bl(iyyz),ncf,
     1                 ii,jj,lbasw(1,ii),lbasw(1,jj),idd,
     2                 jdd,Rij3,Rij4,Rx,Ry,
     3                 Rz,bl(ixxa),bl(ixya),bl(ixza),bl(iyya),
     4                 bl(iyza),bl(izza),x1,y1,z1,
     5                 x2,y2,z2)
C
        wrdw=wrdw+idd*jdd
        if(wrdw.gt.giga) then
        lfil=lfil+1
        wrdw=idd*jdd
        endif
        nkunit=nkunix+lfil
        call wriarr(bl(kaddr),idd*jdd,nkunit)
C
        call matrem('kweak')
        enddo
        enddo
        ravrg=rsum/float(ncoun)
      write(6,*) 'Information about the weak pairs in Bohrs:'
      write(6,*) 'Minimum distance : ',rmin
      write(6,*) 'Maximum distance : ',rmax
      write(6,*) 'Average distance : ',ravrg
      write(6,*) ncoun,' symmetry unique weak pairs'
        call retmark
        end
C=================mkkweak=====================================
       subroutine mkkweak(RK,YYx,YYy,YYz,ncf,
     1                    ii,jj,lbasi,lbasj,idd,
     2                    jdd,rij3,rij4,Rx,Ry,
     3                    Rz,Xxx,Xxy,Xxz,Xyy,
     4                    Xyz,Xzz,x1,y1,z1,
     2                    x2,y2,z2)
C
C        this subroutine constructs the approximate internal
C        exchange matrix Kij for a single pair using multipoles
C
C        Svein Saebo, Faytteville, AR Summer 2000
C
C        called from weakpair
C        result returned in matrix 'kweak'
C
        implicit real*8(a-h,o-z)
        integer*2 lbasi(*),lbasj(*)
      dimension RK(idd,jdd),yyX(ncf,*),yyY(ncf,*),yyZ(ncf,*)
      dimension Xxx(ncf,*),Xxy(ncf,*),Xxz(ncf,*)
      dimension Xyy(ncf,*),Xyz(ncf,*),Xzz(ncf,*)
C
        do iao=1,idd
        my=lbasi(iao)
        xmi=yyX(my,ii)
        ymi=yyY(my,ii)
        zmi=yyZ(my,ii)
        xxmi=xXX(my,ii)
        xymi=xXY(my,ii)
        xzmi=xXZ(my,ii)
        yymi=xYY(my,ii)
        yzmi=xYZ(my,ii)
        zzmi=xZZ(my,ii)
C        second moments with resepct to orbital 1
        xxmi1=xxmi-x1*xmi*2.0d0
        xymi1=xymi-x1*ymi-y1*xmi
        xzmi1=xzmi-x1*zmi-z1*xmi
        yymi1=yymi-y1*ymi*2.0d0
        yzmi1=yzmi-y1*zmi-z1*ymi
        zzmi1=zzmi-z1*zmi*2.0d0
C
        Raxi=rx*xmi+ry*ymi+rz*zmi
C                
              ar2i=xxmi1+yymi1+zzmi1
C
        do jao=1,jdd
        ny=lbasj(jao)
        xnj=yyX(ny,jj)
        ynj=yyY(ny,jj)
        znj=yyZ(ny,jj)
        xxnj=xXX(ny,jj)
        xynj=xXY(ny,jj)
        xznj=xXZ(ny,jj)
        yynj=xYY(ny,jj)
        yznj=xYZ(ny,jj)
        zznj=xZZ(ny,jj)
C        second moments with resepct to orbital 1
        xxnj2=xxnj-x2*xnj*2.0d0
        xynj2=xynj-x2*ynj-y2*xnj
        xznj2=xznj-x2*znj-z2*xnj
        yynj2=yynj-y2*ynj*2.0d0
        yznj2=yznj-y2*znj-z2*ynj
        zznj2=zznj-z2*znj*2.0d0
C
        Rbxj=rx*xnj+ry*ynj+rz*znj
C
C        first dipole dipole terms:
C
C      (R**-3)*(a|rv|i)*(b|rv|j)
C        - 3*(R**-5)*[Rv.(a|rv|i) Rv.(b|rv|j)]
C        note with a normalized Rv as we have the second term
C        should be:
C        - 3*(R**-3)*[Rv.(a|rv|i) Rv.(b|rv|j)]
C
      aribrj=xmi*xnj+ymi*ynj+zmi*znj
C
C        quadrapole contributions
C
C        3/2 * (a|r2|i) * {Rv . (b|rv|j)}
C        -3/2 * (b|r2|j) * {Rv . (a|rv|i)}
C
C
      br2j=xxnj2+yynj2+zznj2
C
        t23=1.5d0*(ar2i*Rbxj-br2j*Raxi)
C
      Rxxbxj=rx*(xxmi1*xnj+xymi1*ynj+xzmi1*znj) +
     1       ry*(xymi1*xnj+yymi1*ynj+yzmi1*znj) +
     2       rz*(xzmi1*xnj+yzmi1*ynj+zzmi1*znj)
C
      Rxxaxi=rx*(xxnj2*xmi+xynj2*ymi+xznj2*zmi) +
     1       ry*(xynj2*xmi+yynj2*ymi+yznj2*zmi) +
     2       rz*(xznj2*xmi+yznj2*ymi+zznj2*zmi)
C
        t45=3.0d0*(Rxxbxj-Rxxaxi)
C
      RRaxxi=rx*(rx*xxmi1+ry*xymi1+rz*xzmi1) +
     1       ry*(rx*xymi1+ry*yymi1+rz*yzmi1) +
     2       rz*(rx*xzmi1+ry*yzmi1+rz*zzmi1)
C
      RRbxxj=rx*(rx*xxnj2+ry*xynj2+rz*xznj2) +
     1       ry*(rx*xynj2+ry*yynj2+rz*yznj2) +
     2       rz*(rx*xznj2+ry*yznj2+rz*zznj2)
C
        t67=7.5d0*(RRbxxj*raxi-RRaxxi*rbxj)

C        this (above)may look ugly, loops over (x,y,z) components
C        would probably be nicer, but no reason to change once it has
C        been done this way
C
        Rk(iao,jao)=Rij3*(aribrj-3.0d0*Raxi*Rbxj)
     1               +Rij4*(t23+t45+t67)
C
C  first line from dipoles second line from quadrapoles
C
        enddo
        enddo
        end
C========updatw==================================================
        subroutine updatw(ncore,ncf,lbasw,ii,jj,
     1                      idd,jdd,dtmaxw,iorbst,ipps)
C
C        uses the internal exchange matrix stored in 'Kweak'
C
C        Svein Saebo Fayetteville AR, Fall 2000
C
C        called from weakpair
C

      use memory

        implicit real*8(a-h,o-z)
c       common/big/bl(30000)
        integer iorbst(*),ipps(*)
        integer*2 lbasw(*)
C
        call matdef('orbi','q',idd,idd)
        call matdef('eigi','d',idd,idd)
        call matdef('orbj','q',jdd,jdd)
        call matdef('eigj','d',jdd,jdd)
        iorba=mataddr('orbi')
        iega=mataddr('eigi')
        jorba=mataddr('orbj')
        jega=mataddr('eigj')
        ipp=ipps(ii)
        jpp=ipps(jj)
        ipoint=iorbst(ii)
        jpoint=iorbst(jj)
        korbs=mataddr('morbs')
        keigs=mataddr('meigs')
        id2=idd*idd
        jd2=jdd*jdd
        do kk=1,id2
        bl(iorba-1+kk)=bl(korbs+ipoint-1+kk)
        end do
        do ll=1,idd
        bl(iega-1+ll)=bl(keigs+ipp-1+ll)
        end do
        do kk=1,jd2
        bl(jorba-1+kk)=bl(korbs+jpoint-1+kk)
        end do
        do ll=1,jdd
        bl(jega-1+ll)=bl(keigs+jpp-1+ll)
        end do
        call matdef('Rmo','r',idd,jdd)
C     call orbfix(ii,jj,ncf,idd,jdd,lbasw)
C  orbitals now in place continue as before
C        transform residuum to temporary MO basis
C
          call matscal('Rij',2.0d0)
        call traMOw(idd,jdd)
C
C        the Rij in MO basis is now in  'Rmo'
        irad=mataddr('Rmo')
C        calculate deltat in MO basis
        ieoa=mataddr('eval')-1
        EO=-bl(ieoa+ii+ncore)-bl(ieoa+jj+ncore)
      call DeltaTw(bl(irad),bl(irad),idd,jdd,bl(iega),bl(jega),
     1 EO,tmax)
C        transform deltaT back to AO basis
        call traaow(idd,jdd)
C
        dtmaxw=max(dtmaxw,tmax)
        call matrem('Rmo')
        call matrem('eigj')
        call matrem('orbj')
        call matrem('eigi')
        call matrem('orbi')
        end
C=========quadrap===========================================
        subroutine quadrap

      use memory

        implicit real*8(a-h,o-z)
c     common /intbl/maxsh,inx(100)
c       common/big/bl(30000)
        call getival('ictr',ictr)
        call getival('ibas',ibas)
        call getival('inuc',inuc)
        call getival('ncs',ncs)
CSS 5 for xx,xy etc??
        nul=0
        ism=5
C        XX
      call inton(ism,nul,bl(mataddr('smxx')),bl(ictr),1,1,bl(ibas),
     1           bl(inuc),ncs)
C        XY
      call inton(ism,nul,bl(mataddr('smxy')),bl(ictr),2,1,bl(ibas),
     1           bl(inuc),ncs)
C        XZ
      call inton(ism,nul,bl(mataddr('smxz')),bl(ictr),3,1,bl(ibas),
     1           bl(inuc),ncs)
C        XY
      call inton(ism,nul,bl(mataddr('smyy')),bl(ictr),2,2,bl(ibas),
     1           bl(inuc),ncs)
C        YZ
      call inton(ism,nul,bl(mataddr('smyz')),bl(ictr),3,2,bl(ibas),
     1           bl(inuc),ncs)
C        ZZ
      call inton(ism,nul,bl(mataddr('smzz')),bl(ictr),3,3,bl(ibas),
     1           bl(inuc),ncs)
        end
C===============mkweaklt===========================
        subroutine mkweaklt(nweakl,nval,invl,listp)
        implicit integer(a-z)
        dimension invl(nval,*),listp(*)
        logical*1 nweakl(nval,*)
        do ii=1,nval
         do jj=1,ii
        ij=invl(ii,jj)
        ijp=listp(ij)
        if(ijp.eq.0) then
        nweakl(ii,jj) =.false.
        nweakl(jj,ii) =.false.        
        else
        nweakl(ii,jj)=.true.
        nweakl(jj,ii)=.true.
        endif
        enddo
        enddo
        end
C=======mkstrng===================================
        subroutine mkstrng(istron,jcmax,nval,nstron)
        implicit integer(a-z)
        integer*2 istron(jcmax,nval)
        logical*1 nstron(nval,*)
        do ii=1,nval
        iput=1
        do kk=1,nval
        if(nstron(kk,ii)) cycle
        iput=iput+1
        istron(iput,ii)=kk
        enddo
        istron(1,ii)=iput-1
        enddo
        end        
C=========mklbmW========================================
      subroutine mklbmW(lbas,idm,lbmap,ncf,nval,
     1                   lbasij,idmij,listp,iddmx,iprint)
C
C        and the union if the domains of ij
C
C        Svein Saebo Fayetteville, AR fall 2000
C
        implicit real*8(a-h,o-z)
        integer*2 lbas(ncf,*),idm(*),lbmap(ncf,*),null,isum
        integer*2 lbasij(iddmx,*),idmij(*),en,to
        dimension listp(*)
        parameter(null=0,en=1,to=2)
C
        epsx=1.0d-07
        do ii=1,nval
        do iao=1,ncf
        lbmap(iao,ii) = null
        enddo
        enddo
        do ii=1,nval
        idd=idm(ii)
        do ild=1,idd
        iao=lbas(ild,ii)
        lbmap(iao,ii)=  en
        enddo
        enddo
        ij=0
        nweak=0
        do ii=1,nval
        do jj=1,ii
        ij=ij+1
        if(listp(ij).ne.0) cycle
        nweak=nweak+1
        do ik=1,ncf
        isum=lbmap(ik,ii)+lbmap(ik,jj)
        if(isum.eq.to) then
        write(6,*)'warning:', nweak, 'have ao ', ik,' in common'
        call nerror(98,'mklbmw','stopping',0,0)
        endif
        enddo
        enddo
        enddo
        ij=0
        nweak=0
        do ii=1,nval
        idd=idm(ii)
        do jj=1,ii
        ij=ij+1
        if(listp(ij).ne.0) cycle
        nweak=nweak+1
        jdd=idm(jj)
        do im=1,idd
        lbasij(im,nweak)=lbas(im,ii)
        enddo
        do jm=1,jdd
        lbasij(idd+jm,nweak)=lbas(jm,jj)
        enddo
        idds=idd+jdd
        ix=-nweak
      call checLBw(idds,lbasij(1,nweak),ncf,epsx,ix,
     1            iprint,lbas(1,ii),lbas(1,jj),idd,jdd)
        idmij(nweak)=idds
        if(idds.ne.idd+jdd) Then
        call nerror(idds,'checklbw','wrong dimension',idd,jdd)
        endif
        enddo
        enddo
        if(iprint.ge.3) then
        write(6,*) 'symmetrical dim. for ',nweak,' weak pairs'
        write(6,*) (idmij(iw),iw=1,nweak)
        endif
        if(iprint.ge.2) then
        iav=0
        do iw=1,nweak
        iav=iav+idmij(iw)
        enddo
        ave=float(iav)/float(nweak)
        write(6,*)'Average dimension for ',nweak,' weak pairs: ',ave
        endif
        end
C=========mksym================================
        subroutine mksym(A,B,idd,jdd,idds)
        implicit real*8(a-h,o-z)
        dimension A(idds,*),B(idd,jdd)
        do jj=1,jdd
        do ii=1,idd
        A(ii,idd+jj)=B(ii,jj)
        enddo
        enddo
        end
C=========checlbw===============================================
        subroutine ChecLBw(idd,lbas,ncf,epsi,ii,
     1                      iprint,lbi,lbj,idi,jdi)
C        this subroutine checks a local domain for redundant functions
C        idd  dimension of local domain
C        lbas local domain
C        ncf number of basis functions
C        epsi  threshold for removing redundant functions
C
C        Svein Saebo Fayetteville, AR June 1999
C        modified for weak pairs September 2000
C        

      use memory

        implicit real*8(a-h,o-z)
        integer*2 lbas(*),lbi(*),lbj(*)
c       common/big/bl(30000)
C
        Parameter(zero=0.0d0)
        ilop=0
C
        iovl=mataddr('ovl')
  100 continue
         call matdef('Oiis','s',idd,idd)
        call matdef('ttms','d',idd*8,idd*8)
        indim=(idd*5)/2+1
        call matdef('iwork','d',indim,indim)
         ifad=idd/2+1
         call matdef('ifail','d',ifad,ifad)
         call matdef('Oii','q',idd,idd)
         isiia=mataddr('Oii')
         call compab(bl(iovl),bl(isiia),ncf,idd,idd,lbas(1),lbas(1))
         call matdef('orbi','q',idd,idd)
         iorbi=mataddr('orbi')-1
         call matdef('eigi','d',idd,idd)
         ieigi=mataddr('eigi')
         call matcopy('Oii','Oiis')
         iossa=mataddr('Oiis')
        ittms=mataddr('ttms')
        iiwork=mataddr('iwork')
         ifai=mataddr('ifail')
        abstol=2*Dlamch('S')
c        call sdiag2(idd,idd,bl(isiia),bl(ieigi),bl(iorbi+1))
         call dspeVX('V','I','U',idd,bl(iossa),zero,zero,1,1,
     1 abstol,mrem,bl(ieigi),bl(iorbi+1),idd,bl(ittms),
     2 bl(iiwork),bl(ifai),info)
         if(info.ne.0) call nerror(1,'dspeVX','diagonalization failed',
     1 info,ii)
        if(bl(ieigi).lt.epsi)then
        vmax=zero
        do lts=1,idd
        val=abs(bl(iorbi+lts))
        if(val.lt.vmax) cycle
        imx=lts
        vmax=val
        end do
        ilop=ilop+1
      if(imx.eq.idd) then
        call matrem('eigi')
        call matrem('orbi')
        call matrem('Oii')
         call matrem('ifail')
        call matrem('iwork')
        call matrem('ttms')
         call matrem('Oiis')
C        go back and check if you need to remove another function
        idd=idd-1
        jdi=jdi-1
        goto 100
        endif
        do irm=imx,idd-1
        lbas(irm)=lbas(irm+1)
        end do
        idd=idd-1
C        function must also be removed from original unsymmetrical
C        bases
        if(imx.eq.idi) then
        idi=idi-1
        goto 100
        endif
        if(imx.gt.idi) then
        imy=imx-idi
        do irm=imy,jdi-1
        lbj(irm)=lbj(irm+1)
        enddo
        jdi=jdi-1
        else
        do irm=imx,idi-1
        lbi(irm)=lbi(irm+1)
        enddo
        idi=idi-1
        endif
        call matrem('eigi')
        call matrem('orbi')
        call matrem('Oii')
         call matrem('ifail')
        call matrem('iwork')
        call matrem('ttms')
         call matrem('Oiis')
C        go back and check if you need to remove another function
        goto 100
        end if
        call matrem('eigi')
        call matrem('orbi')
        call matrem('Oii')
         call matrem('ifail')
        call matrem('iwork')
        call matrem('ttms')
         call matrem('Oiis')
        if(iprint.ge.3) then
        write(6,*) ilop,' functions removed from domain for ',ii
        write(6,*) 'Dimension: ',idd
        write(6,*) 'AOs: ', (lbas(iwr),iwr=1,idd)
        endif
        end
C=========deltaTw==============================================
        subroutine DeltaTw(R,R2,idd,jdd,egi,egj,EO,tmax)
        implicit real*8(a-h,o-z)
        dimension R(idd,jdd),egi(*),egj(*)
        dimension R2(*)
        parameter(zero=0.0d0)
        tmax=zero
        do jj=1,jdd
        del1=eo+egj(jj)
        do ii=1,idd
         del=del1+egi(ii)
        R(ii,jj)=-R(ii,jj)/del
        end do
        end do
C        find the absolute largest element in deltat
        imaxa=idamax(idd*jdd,R2,1)
        tmax=abs(R2(imaxa))
        end
C=========traMOw================================================
        subroutine tramow(idd,jdd)

      use memory

        implicit real*8(a-h,o-z)
c       common/big/bl(30000)
        call matpose('orbi')
        call matdef('tmp1','r',idd,jdd)
        call matmmult('orbi','Rij','tmp1')
        call matmmult('tmp1','orbj','Rmo')
C         irad=mataddr('Rmo')
C         iwst=mataddr('tmp1')
C         call matcopy('Rmo','tmp1')
C         call regen(bl(irad),bl(iwst),idd,jdd)
C
         call matrem('tmp1')
C         call matrem('tempx')
        end
C=========traAOw==============================================
        subroutine traAOw(idd,jdd)
        implicit real*8(a-h,o-z)
        call matpose('orbi')
        call matpose('orbj')
        call matdef('tmp1','r',idd,jdd)
        call matmmult('orbi','Rmo','tmp1')
        call matmmult('tmp1','orbj','Rij')
        call matrem('tmp1')
        end
C=======moveTw=================================================
        subroutine moveTw(lpar,idij,itst)

      use memory

        implicit real*8(a-h,o-z)
        dimension itst(*)
c       common/big/bl(30000)
C
        itfus=mataddr('tful')+itst(lpar)-1
        idtst=mataddr('Rij')-1
        do imv=1,idij
        bl(itfus+imv)=bl(itfus+imv)+bl(idtst+imv)
        end do
        end
C=======penw======================================================
      subroutine penw(lpar,ii,jj,eij,idd,
     1                jdd,itst)

      use memory

        implicit real*8(a-h,o-z)
        dimension itst(*)
c       common/big/bl(30000)
C
       ijdim=idd*jdd
        call matdef('Tij','r',idd,jdd)
        call gettw(lpar,ijdim,itst)
        itij=mataddr('Tij')
         irij=mataddr('Rij')
        call ptrac2(bl(irij),bl(itij),ijdim,eij)
         eij=4.0d0*eij
        call matrem('Tij')
        end
C=======GetTw===================================================
        subroutine getTw(lpar,idd,itst)

      use memory

        implicit real*8(a-h,o-z)
        dimension itst(*)
c       common/big/bl(30000)
C
        itfus=mataddr('tful')+itst(lpar)-1
        itijs=mataddr('Tij')-1
        do imv=1,idd
        bl(itijs+imv)=bl(itfus+imv)
        end do
        end
C================rijmp2w===========================================
      subroutine RijMP2w(lpar,ii,jj,idm,lbas,
     1                  ncf,nval,FF,nstron,nconv,
     2                  invl,listp,nkunix,itst,ifunpair,
     3                  isorb,uniq,equi,nsym,ilist,
     4                  jlist,trij,lbtab,idmx,idmij,
     5                  lbasij,lpt,ndiag,inndp,idd,
     6                  jdd,iff,jff,tcopl,ttmat,
     7                  redu,idmw,lbasw,listw,npairs,
     8                  coro,krec,ifil,lfil,wrdw,
     9                  giga)
C
C        Calculates MP2 residuum using generator state formalism
C
C        Calculates Second order residuum for one pair.
C        lpar pair number
C        ii, jj orbital indices
C        idm(*) dimension of local domains
C        lbas(ncf,*)  contains the local bases
C        ncf  number of contraced basis functions
C        nval number of valence orbitals
C        FF FOck matrix in localized orbitals
C        ilist and jlist give  i and j for a given pair
C        'X1' used for temporary storage ncf*ncf defined in the
C        calling program
C
C        The final residuum will be in 'Rij' in local dimension
C
C        Svein Saebo, Fayetteville, AR Summer 1999
C
C        called from mp2iter
C

      use memory

        implicit real*8(a-h,o-z)
        dimension ifunpair(7,*),isorb(nval,*)
        integer*4 equi(2,*)
        integer*4 uniq(*)
        integer*2 idmw(*),lbasw(*)
        dimension ilist(*),jlist(*)
        dimension krec(*),ifil(*)
        logical*1 nstron(nval,nval)
        logical nconv,ndiag,redu,weak,coro
        dimension invl(nval,*),listp(*)
        dimension FF(nval,nval)
        integer*2 lbas(*),idm(*),idmij(*),lbasij(*)
        integer*2 lbtab(idmx,nsym,*),listw(nval,*)
        integer*2 inndp(*)
        integer*4 itst(*)
c       common/big/bl(30000)
          parameter(half=0.5d0)
C
        call secund(t1rij)
        ifoca=mataddr('fck')
        iovad=mataddr('ovl')
        ijdim=idd*jdd
        call matdef('Rij','r',idd,jdd)
        irad=mataddr('Rij')
C        first calculate FCS and SCF terms
        call matdef('FCS','r',idd,jdd)
        ifcsa=mataddr('FCS')
        call matdef('Tij','r',idd,jdd)
        call getTw(lpar,ijdim,itst)
C
        wrdw=wrdw+ijdim
        if(wrdw.gt.giga) then
        lfil=lfil+1
        wrdw=ijdim
        endif
        nkunit=nkunix+lfil
        call reaarr(bl(irad),ijdim,nkunit)
C
C        FCijS-term
C    calculate F*Cij --> Y1
        call matdef('FCij','r',idd,jdd)
        call matdef('Fii','q',idd,idd)
        ifiia=mataddr('Fii')
       call compab(bl(ifoca),bl(ifiia),ncf,idd,idd,
     1             lbasw(iff),lbasw(iff))
        call matmmult('Fii','Tij','FCij')
        call matrem('Fii')
C    multiply with S from the right and add to Rij
        call matdef('Sjj','q',jdd,jdd)
        isjja=mataddr('Sjj')
       call compab(bl(iovad),bl(isjja),ncf,jdd,jdd,
     1             lbasw(jff),lbasw(jff))
        call matmmult('FCij','Sjj','FCS')
        call matrem('Sjj')
        call matadd('FCS','Rij')
C        SCijF-terms
C        calculate SCij --> FCij
        call matdef('Sii','q',idd,idd)
        isiia=mataddr('Sii')
       call compab(bl(iovad),bl(isiia),ncf,idd,idd,
     1             lbasw(iff),lbasw(iff))
        call matmmult('Sii','Tij','FCij')
        call matrem('Sii')
C    multiply with F from the right and add to Rij
        call matdef('Fjj','q',jdd,jdd)
        ifjja=mataddr('Fjj')
       call compab(bl(ifoca),bl(ifjja),ncf,jdd,jdd,
     1             lbasw(jff),lbasw(jff))
        call matmmult('FCij','Fjj','FCS')
        call matrem('Fjj')
        call matrem('FCij')
        call matadd('FCS','Rij')
C        extra terms
C    calculate F*Cji --> Y1
        call matdef('Tji','r',jdd,idd)
        call matpose2('Tij','Tji','n')
        call matdef('FCji','q',idd,idd)
        call matdef('Fij','r',idd,jdd)
        call matdef('Sij','r',idd,jdd)
        ifija=mataddr('Fij')
        isija=mataddr('Sij')
      call compab(bl(ifoca),bl(ifija),ncf,idd,jdd,
     1            lbasw(iff),lbasw(jff))
        call matmmult('Fij','Tji','FCji')
C    multiply with S from the right and add to Rij
      call compab(bl(iovad),bl(isija),ncf,idd,jdd,
     2            lbasw(iff),lbasw(jff))
        call matmmult('FCji','Sij','FCS')
      call matscal('FCS',-half)
        call matadd('FCS','Rij')
C        SCjiF-terms
C        calculate SCji --> SCji
        call matmmult('Sij','Tji','FCji')
        call matrem('Sij')
C    multiply with F from the right and add to Rij
        call matmmult('FCji','Fij','FCS')
        call matrem('Fij')
        call matrem('FCji')
      call matscal('FCS',-half)
        call matadd('FCS','Rij')
        call matrem('Tji')
C        extra contributions finished
        call matrem('Tij')
C        now calculate coupling terms
        call secund(ttco1)
        call getint_2(ncf,lbun)
C        symmetry  used, multiply by S ouside loop over k
C        only one subroutine here
      call TcouplW(ii,jj,nval,nstron,idm,
     1                   invl,listp,FF,lbas,ncf,
     2                   itst,ifunpair,isorb,uniq,equi,
     3                   nsym,ilist,jlist,lbtab,idmx,
     4                   idmij,lbasij,lpt,ndiag,inndp,
     5                   idd,jdd,iff,jff,ttmat,
     6                   bl(lbun),lbasw,idmw,listw,npairs,
     7                   coro,krec,ifil)
        call retmem(1)
        call secund(ttco2)
        tcopl=tcopl+ttco2-ttco1
        call matrem('FCS')
C
        if(redu) then
        call penalty(lpar,ncf,idd,iff,lbas,
     1                 itst,lpt,lbasij,ndiag)
        endif
C
        call secund(t2rij)
        trij=trij+t2rij-t1rij
        end
C=========tcouplw===================================================
      subroutine TcouplW(ii,jj,nval,nstron,idm,
     1                   invl,listp,FF,lbas,ncf,
     2                   itst,ifunpair,isorb,uniq,equi,
     3                   nsym,ilist,jlist,lbtab,idmx,
     4                   idmij,lbasij,lpt,ndiag,inndp,
     5                   idd,jdd,iff,jff,ttmat,
     6                   lbasun,lbasw,idmw,listw,npairs,
     7                   coro,krec,ifil)
C        calculates coupling terms for the MP2 residuum.
C        with symmetry
C        with this subroutine the multiplications with S are outside
C        the sum over k.
C        Input:
C        ii,jj orbital indices
C        nval number of valence orbitals
C        nstron logical array telling if a pair is strong or not
C        idm  dimensions of the local domains
C        invl inverse pair list
C        listp pair list
C        FF  Fock matrix in localized orbitals
C        lbas local domains
C        ncf number of contracted basis functions
C        itst pointer to the start of T for a given pair
C        ifunpair symmtry relations between AOs
C        isorb  symmetry relations between MOs
C        uniq gives symmetry uniq pair no for unique pairs
C        (see equip in lmp2_sym.f)
C
C        Svein Saebo, Fayetteville AR June 1999
C        Modified by SS August 200
C                

      use memory

        implicit real*8(a-h,o-z)
        logical*1 nstron(nval,*)
        logical coro
        dimension invl(nval,*),listp(*),ilist(*),jlist(*)
        dimension ifunpair(7,*),isorb(nval,*)
        integer*4 uniq(*),equi(2,*)
        dimension krec(*),ifil(*)
        dimension FF(nval,nval)
        integer*2 lbas(*),idm(*),lbtab(idmx,nsym,*),lbasw(*)
        integer*2 idmij(*),lbasij(*),idmw(*),listw(nval,*)
        integer*2 inndp(*),lbasun(*)
        integer*4 itst(*)
        logical wein
c       common/big/bl(30000)
C
        idmxs=idmx*idmx
        itmxa=mataddr('Tmax')
C
        wein=.false.
C
        ixadr=mataddr('x1')
        ioadr=mataddr('ovl')
      call mkunios(jj,ncf,nval,lbas,idm,
     1             lbasij,idmij,invl,listp,nstron,
     2             inndp,lbasun,idun,equi,lbtab,
     3             idmx,nsym,ilist,jlist,wein,
     4             idmw,listw,lbasw,npairs)
        call zerox1(bl(ixadr),ncf,idun,lbasun)
C        loop over valence orbitals kk
        do kk=1,nval
         if(nstron(kk,jj)) cycle
         kjpa=invl(kk,jj)
        kjpar=listp(kjpa)
        iuniq=uniq(kjpar)
        if(kk.eq.jj)then
        lpx=0
        kdkj=idm(kk)
        kfkj=ncf*(kk-1)+1
        else
        lpx=inndp(kjpar)
        kdkj=idmij(lpx)
        kfkj=ncf*(lpx-1)+1
        endif
        kjdim=kdkj*kdkj
      if(jj.gt.kk) then                ! case 1
C        Case 1:  need pair Tkj but jj>kk so we'll have to
C        get Tjk and transpose it
        call matdef('Tkj','q',kdkj,kdkj)
        call matdef('Tij','q',kdkj,kdkj)
        if(iuniq.gt.0)then
        if(coro) then
        call gett(kjpar,kdkj,itst)
        else
C
        itij=mataddr('Tij')
        nunit=ifil(iuniq)
        irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
        endif
        
C
        else
C        getts does the same as gett (puts T into 'Tij') when kjpar
C        is not one of the symmetry related parents stored
      call getTs2(kjpar,ncf,nval,itst,idm,
     1           lbas,equi,ifunpair,isorb,ilist,
     2           jlist,ipt,isym,idmij,lbasij,
     3           lpx,uniq,coro,idmxs,ifil,
     4           krec)
        endif
        call matpose2('Tij','Tkj','n')
        call matrem('Tij')
      else
C        Case 2                simply get Tkj
        call matdef('Tij','q',kdkj,kdkj)
        if(iuniq.gt.0) then
        if(coro) then
        call gett(kjpar,kdkj,itst)
        else
C
        itij=mataddr('Tij')
        nunit=ifil(iuniq)
        irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
        endif
C
        else
      call getTs2(kjpar,ncf,nval,itst,idm,
     1           lbas,equi,ifunpair,isorb,ilist,
     2           jlist,ipt,isym,idmij,lbasij,
     3           lpx,uniq,coro,idmxs,ifil,
     4           krec)
        endif
        call matredef('Tij','Tkj','q',kdkj,kdkj)
      endif
C        all this just to get Tkj which now is in 'Tkj'
C
        fik=-FF(ii,kk)
        call matscal('Tkj',fik)
        itadr=mataddr('Tkj')
C
        if(ii.eq.jj) call nerror(96,'tcouplw','error',0,0)
C
         if(iuniq.gt.0) then
        if(kk.eq.jj) then
        call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,lbas(kfkj),
     1 lbas(kfkj))
         else
        call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,lbasij(kfkj),
     1 lbasij(kfkj))
        endif
        else
C
        if(lpx.gt.0) ipt=nval+lpx
C
         call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,
     1 lbtab(1,isym,ipt),lbtab(1,isym,ipt))
         endif
C
        call matrem('Tkj')
        end do         !loop over kk
C
C        now multiply by S from both sides in combined local dimension
C
        call matdef('xun','q',idun,idun)
C
        ixun=mataddr('xun')
        call compab(bl(ixadr),bl(ixun),ncf,idun,idun,
     1                lbasun,lbasun)
C
          call atoap(bl(ixun),idun)
C
        call matdef('oxm','r',idd,idun)
        call matdef('ovlfj','r',idun,jdd)
        call matdef('ovlif','r',idd,idun)
        iovifa=mataddr('ovlif')
        iovjfa=mataddr('ovlfj')
        call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbasw(iff),lbasun)
        call compab(bl(ioadr),bl(iovjfa),ncf,idun,jdd,
     1               lbasun,lbasw(jff))
        call secund(ttmat1)
        call matmmult('ovlif','xun','oxm')
        call matrem('ovlif')
        call matmmult('oxm','ovlfj','FCS')
        call secund(ttmat2)
        ttmat=ttmat+ttmat2-ttmat1
        call matrem('ovlfj')
        call matrem('oxm')
        call matrem('xun')
         call matadd('FCS','Rij')
C
      call mkunios(ii,ncf,nval,lbas,idm,
     1             lbasij,idmij,invl,listp,nstron,
     2             inndp,lbasun,idun,equi,lbtab,
     3             idmx,nsym,ilist,jlist,wein,
     4             idmw,listw,lbasw,npairs)
        call zerox1(bl(ixadr),ncf,idun,lbasun)
C left half done  do right part
        do kk=1,nval
        if(nstron(ii,kk)) cycle
        ikpa=invl(ii,kk)
        ikpar=listp(ikpa)
        iuniq=uniq(ikpar)
        if(kk.eq.ii) then
        lpx=0
        idik=idm(ii)
        ifik=ncf*(ii-1)+1
        else
        lpx=inndp(ikpar)
        idik=idmij(lpx)
        ifik=ncf*(lpx-1)+1
        endif
        ikdim=idik*idik
      if(kk.gt.ii) then                !Case 1
        call matdef('Tik','q',idik,idik)
        call matdef('Tij','q',idik,idik)
        if(iuniq.gt.0) then
        if(coro) then
        call gett(ikpar,idik,itst)
        else
C
        itij=mataddr('Tij')
        nunit=ifil(iuniq)
        irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
        endif

        else
      call getTs2(ikpar,ncf,nval,itst,idm,
     1            lbas,equi,ifunpair,isorb,ilist,
     2            jlist,ipt,isym,idmij,lbasij,
     3            lpx,uniq,coro,idmxs,ifil,
     4            krec)
        endif
        call matpose2('Tij','Tik','n')
        call matrem('Tij')
      else                                !Case 2
        call matdef('Tij','q',idik,idik)
        if(iuniq.gt.0) then
        if(coro) then
        call gett(ikpar,idik,itst)
        else
        itij=mataddr('Tij')
        nunit=ifil(iuniq)
        irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
        endif
C
        else
      call getTs2(ikpar,ncf,nval,itst,idm,
     1            lbas,equi,ifunpair,isorb,ilist,
     2            jlist,ipt,isym,idmij,lbasij,
     3            lpx,uniq,coro,idmxs,ifil,
     4            krec)
        endif
        call matredef('Tij','Tik','q',idik,idik)
      endif
C
        fkj=-FF(kk,jj)
        call matscal('Tik',fkj)
        itadr=mataddr('Tik')
C
         if(iuniq.gt.0)then
        if(ii.eq.kk) then
        call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,lbas(ifik),
     1 lbas(ifik))
        else
        call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,lbasij(ifik),
     1 lbasij(ifik))
        endif
C
         else
C
        if(lpx.gt.0) ipt=nval+lpx
C
         call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,
     1 lbtab(1,isym,ipt),lbtab(1,isym,ipt))
         endif
        call matrem('Tik')
        end do                                ! end loop over kk
C
C        now multiply by S from both sides in combined local dimension
C
C
C        first determine the combined local local basis
C
C
        call matdef('xun','q',idun,idun)
        ixun=mataddr('xun')
        call compab(bl(ixadr),bl(ixun),ncf,idun,idun,
     1                lbasun,lbasun)
C
           call atoap(bl(ixun),idun)
C
        call matdef('oxm','r',idd,idun)
        call matdef('ovlfj','r',idun,jdd)
        call matdef('ovlif','r',idd,idun)
        iovifa=mataddr('ovlif')
        iovjfa=mataddr('ovlfj')
        call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbasw(iff),lbasun)
        call compab(bl(ioadr),bl(iovjfa),ncf,idun,jdd,
     1               lbasun,lbasw(jff))
        call secund(ttmat1)
        call matmmult('ovlif','xun','oxm')
        call matrem('ovlif')
        call matmmult('oxm','ovlfj','FCS')
        call secund(ttmat2)
        ttmat=ttmat+ttmat2-ttmat1
        call matrem('ovlfj')
        call matrem('oxm')
        call matrem('xun')
         call matadd('FCS','Rij')
        end
C=========lbnwea===========================================
        subroutine lbnwea(ncf,nval,iprint,uniq,ifunpair,
     1                      nsym,lbasij,idmij,epsx,equi,
     2                      ncda,ttda)
C
C        remove redundant functions
C        for weak pairs
C
C        Svein Saebo, Fayetteville, AR Summer 2000
C
        implicit real*8(a-h,o-z)
      dimension ifunpair(7,*)
        integer*2 lbasij(ncf,*),idmij(*)
        integer uniq(*),equi(2,*)
C
C        check local basis for non-diagonal pairs for redundencies
C
        lpar=0
        do ii=1,nval
        idd=idmij(ii)
      call checLB(idd,lbasij(1,ii),ncf,epsx,-ii,
     1            iprint,ncda,ttda)
        idmij(ii)=idd
        end do
C        final check of local domains if symmetry is used
C     if(nsym.gt.0)
C    1  call symche2(equi,ifunpair,lbasij,idmij,nval,
C    2               ncf,iprint,nstron,inndp)
        end
C=========orbfix==============================================
      subroutine orbfix(ii,jj,ncf,idd,jdd,lbas)
C
C     modifies orbitals for pair ii,jj when idd ne jdd
C

      use memory

      implicit real*8(a-h,o-z)
      integer*2 lbas(*)
c     common/big/ bl(30000)
      Parameter(fakt=1.0d-6)
C
        ifu=ncf*(ii-1)+1
        jfu=ncf*(jj-1)+1
C     transform s to MO basis
C
      ijdif=idd-jdd
      if(ijdif.ge.0) then
       call matdef('tempx','d',jdd,jdd)
      else
       call matdef('tempx','d',idd,idd)
      endif
      ixst=mataddr('tempx')
C
      iovs=mataddr('ovl')
      ieigs=mataddr('eigi')
      jeigs=mataddr('eigj')
      call matdef('tempa','r',idd,jdd)
      iast=mataddr('tempa')
      call compab(bl(iovs),bl(iast),ncf,idd,jdd,lbas(ifu),lbas(jfu))
      call matdef('tempb','r',idd,jdd)
      call matpose('orbi')
      call matmmult('orbi','tempa','tempb')
      call matpose('orbi')
      call matmmult('tempb','orbj','tempa')
C        write(6,*) 'S old MO bas:'
C        call matprint('tempa',6)
C    calculate S(tr)*S (jddxjdd)
      call matrem('tempb')
      call matdef('tempc','q',jdd,jdd)
      call matdef('tempb','r',jdd,idd)
      call matpose2('tempa','tempb','n')
      call matmmult('tempb','tempa','tempc')
      call matrem('tempb')
C     **** add fact*(orbital energy) to diagonal
      icad=mataddr('tempc')
      call matdef('tempi','d',jdd,jdd)
      itist=mataddr('tempi')
      jdia=-jdd-1
      do jj2=1,jdd
        jdia=jdia+jdd+1
        bl(icad+jdia)=bl(icad+jdia)+bl(jeigs+jj2-1)*fakt
      enddo
C********  replace sdiag2 with dspevx **********
        call matdef('fsym','s',jdd,jdd)
        call matcopy('tempc','fsym')
        ifsym=mataddr('fsym')
        call eig1(bl(ifsym),bl(icad),bl(itist),jdd,jdd)
C        call sdiag2(jdd,jdd,bl(icad),bl(itist),bl(icad))
        call matrem('fsym')
C********  replace sdiag2 with dspevx **********
      call matrem('tempi')
      call matdef('tempb','q',jdd,jdd)
      call matmmult('orbj','tempc','tempb')
      call matcopy('tempb','orbj')
      call matrem('tempb')
      call matrem('tempc')
C     j-part finished
C     calculate S*S(tr) (iddxidd)
C     S still in 'tempa'
      call matdef('tempc','q',idd,idd)
      call matdef('tempb','r',jdd,idd)
      call matpose2('tempa','tempb','n')
      call matmmult('tempa','tempb','tempc')
      call matrem('tempb')
C    add fact*EO to diagonal
      icad=mataddr('tempc')
      call matdef('tempi','d',idd,idd)
      itist=mataddr('tempi')
      idia=-idd-1
      do ii2=1,idd
        idia=idia+idd+1
        bl(icad+idia)=bl(icad+idia)+bl(ieigs+ii2-1)*fakt
      end do
C********  replace sdiag2 with dspevx **********
        call matdef('fsym','s',idd,idd)
        call matcopy('tempc','fsym')
        ifsym=mataddr('fsym')
        call eig1(bl(ifsym),bl(icad),bl(itist),idd,idd)
C        call sdiag2(idd,idd,bl(icad),bl(itist),bl(icad))
        call matrem('fsym')
C********  replace sdiag2 with dspevx **********
      call matrem('tempi')
      call matdef('tempb','q',idd,idd)
      call matmmult('orbi','tempc','tempb')
      call matcopy('tempb','orbi')
      call matrem('tempb')
      call matrem('tempc')
      call matrem('tempa')
C     transform S to new MO basis
      call matdef('tempa','r',idd,jdd)
      iast=mataddr('tempa')
      call compab(bl(iovs),bl(iast),ncf,idd,jdd,lbas(ifu),lbas(jfu))
      call matpose('orbi')
      call matdef('tempb','r',idd,jdd)
      call matmmult('orbi','tempa','tempb')
      call matpose('orbi')
      call matmmult('tempb','orbj','tempa')
C        write(6,*) 'S new MO bas:'
C        call matprint('tempa',6)
      call matrem('tempb')
C        put diagonal elements in 'tempx'
  666 continue
      if(ijdif.le.0) then
      idia=(-ijdif-1)*idd -1
      DO ii1=1,idd
      idia=idia+idd+1
      bl(ixst-1+ii1)=bl(iast+idia)
      end do
      else
      idia=ijdif-idd-1
      DO jj1=1,jdd
      idia=idia+idd+1
      bl(ixst-1+jj1)=bl(iast+idia)
      end do
      endif
      call matrem('tempa')
C        matrix tempx kept to be used in subroutine regen
C
C     transfor F to new MO basis, I part first
      ifos=mataddr('fck')
      call matdef('tempa','q',idd,idd)
      iast=mataddr('tempa')
      call compab(bl(ifos),bl(iast),ncf,idd,idd,lbas(ifu),lbas(ifu))
      call matsimtr('tempa','orbi','tempa')
C        write(6,*)'Fii new MO'
C        call matprint('tempa',6)
      idia=-idd-1
      DO II2=1,idd
        idia=idia+idd+1
        diag=bl(iast+idia)
        bl(ieigs+ii2-1)=diag
      end do
C        replace eigenvalues with diagonal elments
      call matrem('tempa')
C     j-part
      call matdef('tempa','q',jdd,jdd)
      iast=mataddr('tempa')
      call compab(bl(ifos),bl(iast),ncf,jdd,jdd,lbas(jfu),lbas(jfu))
      call matsimtr('tempa','orbj','tempa')
C        write(6,*)'Fjj new MO'
C        call matprint('tempa',6)
      jdia=-jdd-1
      DO jj2=1,jdd
        jdia=jdia+jdd+1
        diag=bl(iast+jdia)
        bl(jeigs+jj2-1)=diag
      end do
C        call matprint('eigi',6)
C        call matprint('eigj',6)
      call matrem('tempa')
      END
C=======regen=============================================
      subroutine regen(T,W,idd,jdd)
C

      use memory

      implicit real*8(a-h,o-z)
      dimension T(idd,jdd),W(idd,jdd)
c     common /big/bl(30000)
      Parameter(four=4.0d0,two=2.0d0)
C                
C        write(6,*)'regen: old r:'
C        call matprint('Rmo',6)
C        call matprint('tempx',6)
      ixst=mataddr('tempx')
      istart=1
      jstart=1
      if(idd.gt.jdd) istart=idd-jdd+1
      if(jdd.gt.idd) jstart=jdd-idd+1
      j1=-1
      do jj1=jstart,jdd
      j1=j1+1
      i1=-1
      do ii1=istart,idd
      i1=i1+1
      F1=bl(ixst+i1)*bl(ixst+j1)
C
      FAKT=TWO/(FOUR-F1*F1)
      IIP=II1+JSTART-ISTART
      JJP=JJ1+ISTART-JSTART
      T(II1,JJ1)=FAKT*(TWO*W(II1,JJ1)+F1*W(JJP,IIP))
      enddo
      enddo
C        write(6,*)'After regen:'
C        call matprint('Rmo',6)
      end
