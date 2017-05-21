C This contains the subroutines needed for CIM module
C Since both serial and parallel versions need these subroutines,
C I collect them in this independent file
C NZG_6/28/2016 @UARK

c     ##############################################################
c     ##  subroutine NJ_lower  --  make a string lowercase        ##
c     ##  2004.04.16 by Wei Li; Update 2005.10.17 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_lower(line)
      implicit none
      integer i,k,ich
      character line*(*)

      k=len(line)
      do i=1,k
         ich=ichar(line(i:i))
         if (ich.ge.65.and.ich.le.90) line(i:i)=char(ich+32)
      enddo

      end subroutine NJ_lower


c     ##############################################################
c     ##  subroutine NJ_prtcol -- print column(c1:c2) of mat(n,n) ##
c     ##  2005.12.22 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c
c     format: f(x.y) x>7 sugg 12.5,12.7,14.9
c
      subroutine NJ_prtcol(io,n,mat,c1,c2,fm)
      implicit none
      integer i,j,jj,n,io,c1,c2,n5,nf,nc,x,y,k
      real*8 mat(n,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4

      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y

      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
 101  format('(i7,5',a1,i2,'.',i2,')')
 102  format('(i7,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
 103  format('(3x,5(',i2,'x,i7))')
 104  format('(3x,',i2,'(',i2,'x,i7))')

      do jj=1,n5
         write(io,fmt3) (j,j=c1+(jj-1)*5,c1+jj*5-1)
         write(io,fmt1) (i,(mat(i,j),j=c1+(jj-1)*5,c1+jj*5-1),i=1,n)
C         if (jj.ne.n5.or.nf.ne.0) write(io,*)
      enddo

      if (nf.ne.0) then
         write(io,fmt4)(j,j=c1+n5*5,c2)
         write(io,fmt2) (i,(mat(i,j),j=c1+n5*5,c2),i=1,n)
      endif

      end subroutine NJ_prtcol


      subroutine NJ_prtcol5(io,m,n,mat,ntran,nprt,fm)
      implicit none
      integer i,j,jj,m,n,io,c1,c2,n5,nf,nc,x,y,k,nprt(n),ntran(n)
      real(kind=8) mat(m,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4
      integer,allocatable::idxx(:)
C
      allocate(idxx(n))
C
      c1=1
      c2=0
      do j=1,n
         if (nprt(j).ne.0) then
            c2=c2+1
            idxx(c2)=j
         endif
      enddo
      if (c2.eq.0) then
         write(io,*) 'No inputs for this matrix'
         return
      endif
C
C
      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
C
      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
 101  format('(i7,5',a1,i2,'.',i2,')')
 102  format('(i7,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
 103  format('(3x,5(',i2,'x,i7))')
 104  format('(3x,',i2,'(',i2,'x,i7))')
      do jj=1,n5
         write(io,fmt3) (ntran(idxx(j)),j=c1+(jj-1)*5,c1+jj*5-1)
         write(io,fmt1) (ntran(i),(mat(i,idxx(j)),
     &        j=c1+(jj-1)*5,c1+jj*5-1),i=1,m)
      enddo
C
      if (nf.ne.0) then
       write(io,fmt4)(ntran(idxx(j)),j=c1+n5*5,c2)
         write(io,fmt2) (ntran(i),(mat(i,idxx(j)),j=c1+n5*5,c2),i=1,m)
      endif
      deallocate(idxx)
C
      end
CC

      subroutine INDEX_REORDER(io,JF1,JF,nmo,Sover,S2,S3,co1,ZA)
      implicit none
      integer::io,nmo,JF1,JF,ZA(JF1),i,j,k,l
      real(kind=8)::Sover(JF1,JF1),S2(JF1,JF1),S3(JF,JF),co1(JF1,nmo)
      real(kind=8),allocatable::SMO(:)
C
      allocate(SMO(JF1))
      do j=1,JF1
         do k=1,JF1
            S2(j,k)=Sover(ZA(j),ZA(k))
         enddo
      enddo

      do j=1,JF
         do k=1,JF
            S3(j,k)=S2(j,k)
         enddo
      enddo
C
      do i=1,Nmo
         SMO(:)=co1(:,i)
         do j=1,JF1
            co1(j,i)=SMO(ZA(j))
         enddo
      enddo
C
      deallocate(SMO)    
      return
      end subroutine INDEX_REORDER


C --- Project nmo2 of nmo occ. LMOs from total system onto small area
      subroutine projorb(io,nbs,nbs2,nmo,nmo2,mo,mo2,s,wf1,wf2)
      integer io,nbs,nbs2,nmo,i,j,k,L,info
      real(kind=8) mo(nbs,nmo),mo2(nbs2,nmo2),s(nbs,nbs)
      real(kind=8) wf1(nmo2),wf2(nmo2)
      real(kind=8),allocatable::S4(:,:),L1(:,:),work(:),work2(:),L2(:,:)
      real(kind=8) mo3(nbs2,nmo2),A(nbs2,nmo2),t,tstart,tend
      integer lwork,ipiv(nbs2)
C     
      if (nbs2.gt.nbs .or. nmo2.gt.nmo) call nerror(1,'CIM Module',
     &               "Error to project orbital to a larger domain!",0,0)
      if (io>0) then
         if (nmo.eq.nmo2) then
            write(io,100) nmo,nbs,nbs2
         else
            write(io,105) nmo,nmo2,nbs,nbs2
         endif
      endif
100   format(3x,'[Project MOs] MO=',i4,'  AO:',i4,'-->',i4)
105   format(3x,'[Project MOs] MO:',i4,'-->',i4,'  AO:',i4,'-->',i4)
      
      if (nbs2.eq.nbs) then
         do j=1,nmo2
            mo2(:,j)=mo(:,j)
         enddo
      else
         allocate(L1(nbs2,nbs2),S4(nbs2,nbs2),L2(nbs2,nbs2))
         DO IV=1,nbs2
            DO IU=1,nbs2
               S4(IU,IV)=s(IU,IV)
           END DO
         END DO
         
         do j=1,nmo2
            DO I=1,nbs2
               mo2(i,j)=0.0D0
               DO K=1,nbs
                  mo2(i,j)=mo2(i,j)+MO(K,j)*S(k,i)
               end do
            end do
            L1=S4; L2=S4
            mo3=mo2

C --- Zhigang 2016.3.2 Use MKL subroutine to solve the linear equation
            call dpotrf('L',nbs2,L1,nbs2,info)
            call dpotrs('L',nbs2,1,L1,nbs2,mo2(:,j),nbs2,info)
            if (info/=0) call nerror(1,'CIM Module',
     &               "Fail to solve the linear equation!",0,0)

C --- Zhigang 2016.9.14 Try to use matrix formulation 
C --- J. Boughton, P. Pulay, J. Comp. Chem., 14, 736 (1993)
C            call elapsec(tstart)
C            lwork=64*nbs2
C            allocate(work(lwork),work2(nbs2))
C            call dsytrf('L',nbs2,L2,nbs2,ipiv,work,lwork,info)
C            call dsytri('L',nbs2,L2,nbs2,ipiv,work2,info)
C            deallocate(work,work2)
C            do i=1,nbs2
C               do k=i+1,nbs2
C                  L2(i,k)=L2(k,i)
C               enddo
C            enddo
C            A(:,j)=0.0D0
C            call matmul_mkl(L2,mo3(:,j),A(:,j),nbs2,nbs2,1)
C            call elapsec(tend)
C            t=(tend-tstart)/60.0D0
C            write(6,*) "time for matrix multiply",t
C            
C            do i=1,nbs2
C               if (mo2(i,j)-A(i,j)>1.0D-6) then
C                  write(6,*) "different"
C                  stop
C               endif
C            enddo
C --- The method is OK but it is not faster the the linear equation way

         end do
      endif

      do j=1,nmo2
         call WFUNC(io,nbs,s,nbs2,mo(1,j),mo2(1,j),wf1(j),wf2(j))
      enddo
     
      end subroutine projorb


      subroutine WFUNC(io,n,s,m,mo1,mo2,P1,P2)
      integer i,j,k,l,m,n,mu,nu,io
      real(kind=8) s(n,n),mo1(n),mo2(m),P1,P2,PP
      real(kind=8),allocatable::c(:),d(:)

      allocate(c(n))
      do mu  = 1,n
         if (mu.le.m) then
            c(mu) = mo1(mu)-mo2(mu)
         else
            c(mu) = mo1(mu)
         endif
      enddo
C
      P1 = 0.0D+00
      do mu = 1,m
         do nu = 1,m
            P1 = P1 + c(mu)*c(nu)*s(nu,mu)
         enddo
      enddo
C
      P2 = 0.0D+00
      do mu = 1,n
         do nu = 1,n
            P2 = P2 + c(mu)*c(nu)*s(nu,mu)
         enddo
      enddo

      deallocate(c)
      end subroutine WFUNC
C
C
      subroutine locindx_occ(io,nbs1,mo1,nbs2,mo2,nmo,s,indx,eta,nstar)
      implicit none
      integer io,nbs1,nbs2,nmo,i,j,L,i1,i2,j1,j2,nstar,imax,imin
      real(kind=8) mo1(nbs1,nmo),mo2(nbs2,nmo),s(nbs2,nbs2)
      real(kind=8) p1,pp,over,pmax,pmin,indx(nmo),eta
      character(len=1),allocatable::star(:)
C
      allocate(star(nmo))
C
      do i=1,nmo
         star = ' '
      enddo
C
      nstar=0
      do i=1,nmo
         p1=0.0d0
         do j2=1,nbs2
            do j1=1,nbs2
              p1=p1+mo2(j1,i)*mo2(j2,i)*s(j1,j2) ! overlap
            end do
         end do
         indx(i)=1.0d0-p1
**       indx(i)=dabs(p1-1.0d0)
         if (indx(i).lt.eta) then
            star(i) = '*'
            nstar = nstar + 1
         endif
      end do
C
      pmax=indx(1);imax=1;pmin=indx(1);imin=1
      do i=2,nmo
         if (indx(i).gt.pmax) then
            pmax=indx(i)
            imax=i
         endif
         if (indx(i).lt.pmin) then
            pmin=indx(i)
            imin=i
         endif
      enddo
C
      if (io.gt.0) then
         write(io,100) nmo,nbs1,nbs2
         write(io,110) (i,indx(i),star(i),i=1,nmo)
         write(io,120) nstar,eta
         write(io,130) pmax,imax,pmin,imin
      end if
C
      deallocate(star)
 100  format(1x,'local index of',i5,' MOs with',i5,' ->',i5,' AOs')
 110  format(5(i5,'-',f7.4,a1))   
 120  format(1x,'*',i5,' localization indices little than',f7.4)
 130  format(1x,'Max:',f7.4,' in',i5,' and Min:',f7.4,' in',i5)
      end subroutine locindx_occ

      subroutine Add_Atom(io,nvi,NATOM,link,dis,BA,J0,J1,thre)
      implicit none
      integer::io,nvi,NATOM,indx,i,j,k,l,AA,BB,J0,J1,J00
      integer::SNBsub(nvi),SOBsub(nvi,NATOM),link(NATOM,NATOM)
      integer::iindex(NATOM),BA(NATOM)
      real(kind=8)::dis(NATOM,NATOM),thre

      iindex(:)=0;J00=J0
      do i=J0+1,NATOM
         AA=BA(i)
         do j=1,J0
            BB=BA(j)
            if (dis(AA,BB)<thre) then
               J00=J00+1;iindex(J00)=i
               exit
            endif
         end do
      end do
C
      do i=J0+1,J00
         AA=BA(i)
         BA(i)=BA(iindex(i))
         BA(iindex(i))=AA
      end do 
      J1=J00
      if (io>0) then
         WRITE(io,'(I8,a27,f8.2,a)')
     &             J00-J0,'ATOM ADDED TO AO DOMAIN IN',thre,' Angstrom'
         WRITE(io,'(10I8)')(BA(i),i=J0+1,J00)
         WRITE(io,'(I8,a27,f8.2,a)')
     &             J00-J0,'ATOM ADDED TO AO DOMAIN IN',thre,' Angstrom'
         WRITE(io,*)''
      endif
C
      return
      end subroutine Add_Atom

      subroutine BASIS_REORDER(io,NW,NATOM,ZA,nb,BA,J0,JF)
      implicit none
      integer::NW,io,NATOM,J0,JF,JF1,K1,K2,J,K
      integer::ZA(NW),BA(NATOM),nb(NATOM)
       
      ZA=0;JF=0
      DO J=1,J0
         K1=nb(BA(J))
         if (BA(j).ne.NATOM) then
            K2=nb(BA(J)+1)-1
         else
            K2=NW
         endif
         DO K=K1,K2
            JF=JF+1
            ZA(JF)=K
         end do
      end do
      JF1=JF
      DO J=J0+1,NATOM
         K1=nb(BA(J))
         if (BA(j).ne.NATOM) then
            K2=nb(BA(J)+1)-1
         else
            K2=NW
         endif
         DO K=K1,K2
            JF1=JF1+1
            ZA(JF1)=K
         end do
      end do
      if(JF1/=NW) then
        write(io,*)'THERE MUST STH WRONG WHEN DETER BASIS INDEX',JF1,NW
        stop
      end if
      return
      end subroutine BASIS_REORDER
C
C
      subroutine locindx_red(io,nbs1,mo1,nbs2,mo2,nmo,s,indx,eta,
     &                       nstar,paored)
      implicit none
      integer io,nbs1,nbs2,nmo,i,j,k,L,i1,i2,j1,j2,k1,k2,nstar
      real(kind=8) mo1(nbs1,nmo),mo2(nbs2,nmo),s(nbs2,nbs2)
      real(kind=8) p1,pp,over,pmax,pmin,indx(nmo),eta
      integer::paored(nmo)
C
      nstar=0; paored=1
      do i=1,nmo
         p1=0.0d0
         do j2=1,nbs2
            do j1=1,nbs2
               p1=p1+mo2(j1,i)*mo2(j2,i)*s(j1,j2) ! overlap
            end do
         end do
         indx(i)=1.0d0-p1
         if (indx(i)<=eta) then
            nstar=nstar+1
            paored(i)=0 !paored=0 means not be reduced else will be reduced
         endif
      end do
      end subroutine locindx_red


C***** Construct the Full Atomic AO Domain. 2015.1.14 Zhigang *****
      subroutine Full_AAO(NATOM,nat_ao,J1_aao,BA_AAO,BA_FAAO,nat_faao)
      implicit none
      integer::NATOM,nat_ao,nat_faao,outfaao,j,k,l,in
      integer::BA_AAO(NATOM,nat_ao),BA_FAAO(NATOM),J1_AAO(nat_ao)

      BA_FAAO(:J1_aao(1))=BA_AAO(:J1_aao(1),1)
      nat_faao=J1_aao(1)
      do j=2,nat_ao
         do k=1,J1_aao(j)
            call inarray(BA_AAO(k,j),BA_FAAO(:nat_faao),nat_faao,in)
            if (in==0) then
               nat_faao=nat_faao+1
               BA_FAAO(nat_faao)=BA_AAO(k,j)
            end if
         end do
      end do

      outfaao=nat_faao
      do j=1,NATOM
         call inarray(j,BA_FAAO(:nat_faao),nat_faao,in)
         if (in==0) then
            outfaao=outfaao+1
            BA_FAAO(outfaao)=j
         end if
      end do
      return
      end subroutine Full_AAO


      subroutine MKL_Tfock(io,n,m,Mat,X,Mat2)
      implicit none
      integer io,i,j,k,l,m,n
      real*8 Mat(n,n),X(n,m),Mat2(m,m),PP,P1,P2
      real(kind=8),dimension(:,:),allocatable::Y,XT
      allocate(Y(n,m))

      P1=1.0;P2=0.0
      Y=0.0
      call DGEMM('N','N',n,m,n,P1,Mat,n,X,n,P2,Y,n)
      Mat2=0.0
      call DGEMM('t','N',m,m,n,P1,X,n,Y,n,P2,Mat2,m)
      deallocate(Y)

      return
      end subroutine MKL_Tfock


C -- FOR DIAGONALIZATION. 2003.5.25.
        SUBROUTINE NJ_qr(io,A,N,Q,B,C,L,sort)
        implicit none
        INTEGER N,L,io,sort,i,j,k
        real*8 A(N,N),Q(N,N),B(N),C(N),valu,vect,eps
        eps=1.0d-16
        CALL CSTRQ(A,N,Q,B,C,eps)
        CALL CSSTQ(N,B,C,Q,L,eps)
        if (L.eq.0) call nerror(1,'CIM Module',
     &                          "QR Diagonalization failure!",0,0)
C --- Sort for eigenvalues ---
        if (sort>0) then   ! order eigenvalues from small to large
           do i=1,N-1
              valu=B(i)
              do j=i+1,N
                 if (B(j).lt.valu) then
                    valu=B(j)
                    B(j)=B(i)
                    B(i)=valu
                    do k=1,N
                       vect=Q(k,i)
                       Q(k,i)=Q(k,j)
                       Q(k,j)=vect
                    enddo
                 endif
              enddo
           enddo
           if (io.gt.0) write(io,110) N
        elseif (sort<0) then
           do i=1,N-1
              valu=B(i)
              do j=i+1,N
                 if (B(j).gt.valu) then ! order eigenvalues from large to small
                    valu=B(j)
                    B(j)=B(i)
                    B(i)=valu
                    do k=1,N
                       vect=Q(k,i)
                       Q(k,i)=Q(k,j)
                       Q(k,j)=vect
                    enddo
                 endif
              enddo
           enddo
           if (io.gt.0) write(io,120) N
        else
           if (io.gt.0) write(io,130) N
        endif
        if (io.gt.0) write(io,140) (B(i),i=1,N)
 110  format(1x,'QR diagonalization with',i4,' ascending eigenvalues')
 120  format(1x,'QR diagonalization with',i4,' descending eigenvalues')
 130  format(1x,'QR diagonalization with',i4,' unsorted eigenvalues')
 140  format(5d14.6)
      END
C
        SUBROUTINE CSSTQ(N,B,C,Q,L,EPS)
        real*8 Q(N,N),B(N),C(N)
        real*8 EPS,D,H,P,R,F,E,S,G
        INTEGER L
        C(N)=0.0d0
        D=0.0d0
        F=0.0d0
        DO 50 J=1,N
           IT=0
           H=EPS*(DABS(B(J))+DABS(C(J)))
           IF (H.GT.D) D=H
           M=J-1
10         M=M+1
           IF (M.LE.N) THEN
              IF (DABS(C(M)).GT.D) GOTO 10
           END IF
           IF (M.NE.J) THEN
15            IF (IT.EQ.60) THEN
                 L=0
                 WRITE(*,18)
18               FORMAT(1X,'  FAIL')
                 WRITE(11,*) 'FAIL IN QR'
                 RETURN
              END IF
              IT=IT+1
              G=B(J)
              P=(B(J+1)-G)/(2.0d0*C(J))
              R=DSQRT(P*P+1.0d0)
              IF (P.GE.0.0d0) THEN
                 B(J)=C(J)/(P+R)
              ELSE
                 B(J)=C(J)/(P-R)
              END IF
              H=G-B(J)
              DO 20 I=J+1,N
20               B(I)=B(I)-H
                 F=F+H
                 P=B(M)
            E=1.0d0
            S=0.0d0
            DO 40 I=M-1,J,-1
              G=E*C(I)
              H=E*P
              IF (DABS(P).GE.DABS(C(I))) THEN
                E=C(I)/P
                R=DSQRT(E*E+1.0d0)
                C(I+1)=S*P*R
                S=E/R
                E=1.0d0/R
              ELSE
                E=P/C(I)
                R=DSQRT(E*E+1.0d0)
                C(I+1)=S*C(I)*R
                S=1.0d0/R
                E=E/R
              END IF
              P=E*B(I)-S*G
              B(I+1)=H+S*(E*G+S*B(I))
              DO 30 K=1,N
                H=Q(K,I+1)
                Q(K,I+1)=S*Q(K,I)+E*H
                Q(K,I)=E*Q(K,I)-S*H
30            CONTINUE
40          CONTINUE
            C(J)=S*P
            B(J)=E*P
            IF (DABS(C(J)).GT.D) GOTO 15
          END IF
          B(J)=B(J)+F
50      CONTINUE
        DO 80 I=1,N
          K=I
          P=B(I)
          IF (I+1.LE.N) THEN
            J=I
60          J=J+1
            IF (J.LE.N) THEN
              IF (B(J).LE.P) THEN
                K=J
                P=B(J)
                GOTO 60
              END IF
            END IF
          END IF
          IF (K.NE.I) THEN
            B(K)=B(I)
            B(I)=P
            DO 70 J=1,N
              P=Q(J,I)
              Q(J,I)=Q(J,K)
              Q(J,K)=P
70          CONTINUE
          END IF
80      CONTINUE
        L=1
        RETURN
        END

        SUBROUTINE CSTRQ(A,N,Q,B,C,eps)
        real*8 A(N,N),Q(N,N),B(N),C(N)
        real*8 F,H,G,H2,eps
        DO 10 I=1,N
        DO 10 J=1,N
10      Q(I,J)=A(I,J)
        DO 80 I=N,2,-1
          H=0.0d0
          IF (I.GT.2) THEN
            DO 20 K=1,I-1
20          H=H+Q(I,K)*Q(I,K)
          END IF
          IF (dabs(H).le.eps) THEN
            C(I)=0.0d0
            IF (I.EQ.2) C(I)=Q(I,I-1)
            B(I)=0.0d0
          ELSE
            C(I)=DSQRT(H)
            IF (Q(I,I-1).GT.0.0d0) C(I)=-C(I)
            H=H-Q(I,I-1)*C(I)
            Q(I,I-1)=Q(I,I-1)-C(I)
            F=0.0d0
            DO 50 J=1,I-1
              Q(J,I)=Q(I,J)/H
              G=0.0d0
              DO 30 K=1,J
30            G=G+Q(J,K)*Q(I,K)
              IF (J+1.LE.I-1) THEN
                DO 40 K=J+1,I-1
40              G=G+Q(K,J)*Q(I,K)
              END IF
              C(J)=G/H
              F=F+G*Q(J,I)
50          CONTINUE
            H2=F/(H+H)
            DO 70 J=1,I-1
              F=Q(I,J)
              G=C(J)-H2*F
              C(J)=G
              DO 60 K=1,J
60            Q(J,K)=Q(J,K)-F*C(K)-G*Q(I,K)
70          CONTINUE
            B(I)=H
          END IF
80      CONTINUE
        DO 85 I=1,N-1
85      C(I)=C(I+1)
        C(N)=0.0d0
        B(1)=0.0d0
        DO 130 I=1,N
          IF ((dabs(B(I)).gt.eps).AND.(I-1.GE.1)) THEN
            DO 110 J=1,I-1
              G=0.0d0
              DO 90 K=1,I-1
90            G=G+Q(I,K)*Q(K,J)
              DO 100 K=1,I-1
100           Q(K,J)=Q(K,J)-G*Q(K,I)
110         CONTINUE
          END IF
          B(I)=Q(I,I)
          Q(I,I)=1.0d0
          IF (I-1.GE.1) THEN
            DO 120 J=1,I-1
              Q(I,J)=0.0d0
              Q(J,I)=0.0d0
120         CONTINUE
          END IF
130     CONTINUE
        RETURN
        END

C --- Count the number of values \in [z1,z2] of A(n)
      subroutine numpoint(n,A,z1,z2,num)
      implicit none
      integer n,num,i,j
      real(kind=8) A(n),z1,z2,zzz
C
      num=0
      do i=1,n
         zzz=A(i)
         if (zzz.lt.z1) cycle
         if (zzz.gt.z2) cycle
         num=num+1
      end do
C
      end


c *************************************************************** 
C * Reorder Fock matrix so basis functions are ordered per atom *
C * as opposed to Wolinski special ordering                     *
c * Modify the subroutine ReorderMO                             *
C * NZG_5/4/2016 @UARK                                          *
C ***************************************************************
      SUBROUTINE ReorderFock(NBas,NOcc,INDX,COld,CNew)
      IMPLICIT REAL*8(A-H,O-Z)

C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  NOcc    -  number of occupied MOs (Here nocc==nabs)
C  INDX    -  index array relating atom ordering to Wolinski ordering
C             (output from subroutine <SortBAS2>)
C  COld    -  Full set of Fock matrix in Wolinski order
C  CNew    -  on output contains Fock matrix in atom order
C
C
      DIMENSION INDX(NBas),COld(NBas,nocc),CNew(nbas,nocc)
C
C
      DO I=1,NBas
         II = INDX(I)
         DO J=1,NOcc
            CNew(II,J) = COld(I,J)
         end do
      end do

      cold=cnew
      do i=1,nocc
         ii=indx(i)
         do j=1,nbas
            cnew(j,ii)=cold(j,i)
         end do
      end do
C
      RETURN
      END subroutine ReorderFock


c ********************************************************* 
C * Reorder Fock matrix from atom order to Wolinski order *
c * It is opposite to subroutine ReoderFock               *
C * NZG_6/13/2016 @UARK                                   *
C *********************************************************
      SUBROUTINE ReorderFock2(NBas,NOcc,INDX,COld,CNew)
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION INDX(NBas),COld(NBas,nocc),CNew(nbas,nocc)
C
      DO I=1,NBas
         II = INDX(I)
         DO J=1,NOcc
            CNew(I,J) = COld(II,J)
         end do
      end do

      cold=cnew
      do i=1,nocc
         ii=indx(i)
         do j=1,nbas
            cnew(j,i)=cold(j,ii)
         end do
      end do
C
      RETURN
      END subroutine ReorderFock2


c ********************************************************* 
C * Reorder Fock matrix from atom order to Wolinski order *
c * It is opposite to subroutine ReoderFock               *
c * Store the Fock matrix in a symmetric matrix           *
c * Modify subroutine ReorderFock2                        *
C * NZG_3/30/2017 @UARK                                   *
C *********************************************************
      SUBROUTINE ReorderFock3(NBas,INDX,COld,CNew2)
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION INDX(NBas),COld(NBas,NBas),CNew2(nbas*(nbas+1)/2)
C
      real*8 CNew(nbas,nbas)

      DO I=1,NBas
         II = INDX(I)
         DO J=1,NBas
            CNew(I,J) = COld(II,J)
         end do
      end do

      cold=cnew
      do i=1,nbas
         ii=indx(i)
         do j=1,nbas
            cnew(j,i)=cold(j,ii)
         end do
      end do
C
      call squa2sym(nbas,cnew,cnew2)

      RETURN
      END subroutine ReorderFock3


      subroutine squa2sym(nbas,cnew,cnew2)
      implicit none
      
      integer nbas,i,j,ij
      real*8 cnew(nbas,nbas),cnew2(nbas*(nbas+1)/2)

      ij=0
      do i=1,nbas
         do j=1,i
            ij=ij+1
            cnew2(ij)=cnew(i,j)
         end do
      end do
      
      return
      end subroutine squa2sym


C --- 2005.11.19 --- Transformation matrix for overlap by Canonical Orthogonalization  ---
      subroutine NJ_canorth(io,nbs,nif,Evalu,Evect,Tmatr)
      implicit none
      integer io,nbs,i,j,k,nif
      real*8 Evalu(nif),Evect(nbs,nif),Tmatr(nbs,nif),PP

      Tmatr=0d0
      do j=1,nif; do i=1,nbs
            PP=Evect(i,j)/dsqrt(Evalu(j))
            Tmatr(i,j)=Tmatr(i,j)+PP
      end do; end do

      if (io>0) then
         write(io,*) 'Transformation Matrix by Canonical',
     &               ' Orthogonalization'
         write(io,*)
      endif

      end subroutine NJ_canorth


C ***********************************************
C * Gram-Schmidt orthogonalization for orbitals *
C * NZG_4/18/2017 @UARK                         *
C ***********************************************
      subroutine Schmidt_orth(nbas,norb,C,S)
      implicit none
      integer nbas,norb,i,j
      real*8 C(nbas,norb),S(nbas,nbas),tmp(nbas,1),tmp2(nbas,1)
      real*8 Sij

      do i=2,norb
         do j=1,i-1
            tmp2(:nbas,1)=C(:,j)
            tmp=matmul(S,tmp2)
            Sij=-dot_product(C(1:nbas,i),tmp(:,1))
            C(:,i)=C(:,i)+Sij*C(:,j)
         enddo
         call normorb(nbas,1,C(:,i),S)
      enddo

      return
      end subroutine Schmidt_orth
   

C****** 2004.12.21 trace for a matrix [real*8]******
      function dtrace2(n,mat)
      implicit none
      integer n,i
      real*8 mat(n,n),dtrace2
      dtrace2=0.0d0
      do i=1,n
         dtrace2=dtrace2+mat(i,i)
      enddo
      end function dtrace2



c  ################################################################
c  ##  Make PQS input of clusters                                ##
c  ##  Originally in GAMESS-2004.12.26; Update 2005.12.26 by WL  ##
C  ##  5/1/2016_NZG @UARK                                        ##
c  ################################################################
      subroutine SUBINP(inp,nat,AtSymb,coor,snat,BA,clumtd,mwords,
     &                  icharg,mult,lineardep,calforce)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer inp,nat,snat,BA(snat),i,j,k,l,m,n,m1,m2,clumtd
      real(kind=8) coor(3,nat)
      integer k1,k2,k3,k4,k5,k6,stat,mplevl,lenJ,nfor,lineardep
      integer memnew,icharg,L1,L2,mult,L3,L4,mwords,L5,L6
      character*256 line,jobname,basline
      character charg*100,mulp*100
      character*8 AtSymb(nat)
      character*100 forline(30)
      logical calforce
      data IUnit   /1/ ! unit number for checkpoint I/O
      common /job/ jobname,lenJ
     
C
      rewind(inp)
C
C --- MEMORY
      write(line,*) mwords
      call NJ_trim(line,k1,k2)
      line='%MEM=200'
      write(inp,100) trim(line)
C
C --- GEOMETRY
      write(charg,*) icharg
      call NJ_trim(charg,L1,L2)
      write(mulp,*) mult
      call NJ_trim(mulp,L3,L4)
      line='GEOM=PQS '//'CHAR='//charg(L1:L2)//' MULT='//mulp(L3:L4)//
     &     ' SYMM=0'
      write(inp,100) trim(line)

      do k=1,snat
         L=BA(k)
         write(inp,120) AtSymb(L),(coor(j,L),j=1,3)
      enddo
C
C --- BASIS
      open(unit=iunit,file=jobname(1:lenJ)//'.inp')
      nfor=-1
      do while(.true.)
         read(unit=iunit,fmt="(A50)",iostat=stat) line
         call NJ_trim(line,L5,L6)
         call NJ_lower(line(L5:L6))
         if (line(L5:L5+4)=='basis') basline=line
         if (line(L6-3:L6)=='next') then
            nfor=0
            cycle
         endif
         if (line(L5:L5+2)=='for' .and. nfor>-1) then
            nfor=nfor+1      
            forline(nfor)=line
         endif
         if (line(L5:L5+2)=='cim') then
            write(inp,100) trim(basline)
            do k=1,nfor
               write(inp,100) trim(forline(k))
            enddo
            exit
         end if
      end do
      close(iunit)
C      write(inp,100) 'BASIS NEXT'
C      open(unit=iunit,file=jobname(1:lenJ)//'.basis')
C      do while(.true.)
C         read(iunit,fmt="(A50)",iostat=stat) line
C         if (stat/=0) exit
C         if (line(1:1)=='$') then
C            cycle
C         else
C            write(inp,100) trim(line)
C         end if
C      end do
C
C --- METHOD
      select case(clumtd)
         case(1)
            line='CIMSUB RIMP2'
         case(2)
            line='CIMSUB MP2'
         case(3)
            line='CIMSUB CCD'
         case(4)
            if (lineardep==1) then
               line='CIMSUB CCSD THRE=14'
            else
               line='CIMSUB CCSD'
            endif
         case(5)
            line='CIMSUB CCSD TRIPLES'
      end select
      if (calforce) line=trim(line)//' FORCE'
      write(inp,100) trim(line)

 100  format(a)
 120  format(A8,3f16.8)
C
      end subroutine SUBINP


      subroutine iwrit8(iunit, key, n, M)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      integer M(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(8I10)') (M(k),k=1,n)
      write(iunit,*) '$END'
      end subroutine iwrit8


      subroutine rwrit8(iunit, key, n, A)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      double precision A(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(1P,4E20.12)') (A(k),k=1,n)
      write(iunit,*) '$END'
      end subroutine rwrit8


      subroutine rread8(iunit, key, n, A)
      implicit none
      character key*(*),line*200
      integer iunit,n,L,k,iyes,itry
      double precision A(n)
C
      A=0.0D0
C
      iyes=0
      itry=0
      L=len(key)
 50   do while(iyes.eq.0)
         read(iunit,'(a)',end=100,err=100) line
         if (index(line,key(1:L)).ne.0) then
            iyes=1
            goto 100
         endif
      end do
C
 100  itry=itry+1
      if (iyes.eq.0) then
         if (itry.eq.1) then
            rewind(iunit)
            goto 50
         else
            write(*,*) 'ERROR: can not find '//key(1:L)//' in .txt'
            stop
         endif
      else
         read(iunit,'(1P,4E20.12)') (A(k),k=1,n)
         read(iunit,*,end=999)
      endif

 999  return

      end subroutine rread8


c     ###########################################################
c     ##  subroutine dislink -- atoms distance & link matrix   ##
c     ##  2005.01.09 by Wei Li; Update 2005.10.16 by Wei Li    ##
c     ###########################################################
c
c     2005.04.12 modify ad from 0.2 to 0.1
c     2005.05.07 modify ad from 0.1 to 0.15
c     2016.04.19 modify nuchar to IAN -NZG @UARK
      subroutine dislink(io,nat,IAN,coor,dis,link)
      implicit none
      integer io,nat,link(nat,nat),i,j,k,l,lk(20),IAN(nat)
      real*8 coor(3,nat),dis(nat,nat),x,y,z,dis0,radius,ad
      real*4 bdtp(20)
      parameter (ad=0.168d0) ! 2005.05.22 0.16 --> 0.168
      external radius

      dis=0.0d0; link=0
      do i=1,nat-1
         do j=i+1,nat
            x=coor(1,i)-coor(1,j)
            y=coor(2,i)-coor(2,j)
            z=coor(3,i)-coor(3,j)
            dis(i,j)=dsqrt(x*x+y*y+z*z)
            dis(j,i)=dis(i,j)
            k=IAN(i)
            l=IAN(j)
            dis0=radius(k)+radius(l)+ad
            if (dis(i,j)<dis0) then
               link(i,j)=1
               link(j,i)=1
            endif
         enddo
      enddo

      if (io>0) then
         write(io,*) '*** Distance matrix (Angstroms) ***'
         call NJ_prtsym(io,nat,dis,'f11.6')
         write(io,*)
         write(io,*) '*** Geometry links ***'
         do i=1,nat
            lk=0; bdtp=0.0D0; k=0
            do j=i+1,nat
               if (link(i,j).ne.0) then
                  k=k+1; lk(k)=j
                  bdtp(k)=real(link(i,j))
               endif
            enddo
            write(io,'(i6,2x,10(i6,f4.1))') i,(lk(j),bdtp(j),j=1,k)
         enddo
         write(io,*)
      endif

      end subroutine dislink

c     ##############################################################
c     ##  function radius  --  radius for an atom on order        ##
c     ##  2005.01.09 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c
c     2005.05.07 modify r(S) from 1.01 to 1.11
c
      function radius(n)
      implicit none
      integer n
      real*8 radius,Radi(109)

      data Radi/0.30d0,1.16d0,1.23d0,0.89d0,0.88d0,
     &     0.77d0,0.70d0,0.66d0,0.58d0,0.55d0,
     &     1.40d0,1.36d0,1.25d0,1.17d0,1.10d0,
     &     1.11d0,0.99d0,1.58d0,2.03d0,1.74d0,
     &     1.44d0,1.32d0,1.20d0,1.13d0,1.17d0,
     &     1.16d0,1.16d0,1.15d0,1.17d0,1.25d0,
     &     1.25d0,1.22d0,1.21d0,1.17d0,1.14d0,
     &     1.89d0,2.25d0,1.92d0,1.62d0,1.45d0,
     &     1.34d0,1.29d0,1.23d0,1.24d0,1.25d0,
     &     1.28d0,1.34d0,1.41d0,1.50d0,1.40d0,
     &     1.41d0,1.37d0,1.33d0,2.09d0,2.35d0,
     &     1.98d0,1.69d0,1.65d0,1.65d0,1.64d0,
     &     1.64d0,1.66d0,1.85d0,1.61d0,1.59d0,
     &     1.59d0,1.58d0,1.57d0,1.56d0,1.70d0,
     &     1.56d0,1.44d0,1.34d0,1.30d0,1.28d0,
     &     1.26d0,1.26d0,1.29d0,1.34d0,1.44d0,
     &     1.55d0,1.54d0,1.52d0,1.53d0,1.52d0,
     &     1.53d0,2.45d0,2.02d0,1.70d0,1.63d0,
     &     1.46d0,1.40d0,1.36d0,1.25d0,1.57d0,
     &     1.58d0,1.54d0,1.53d0,1.84d0,1.61d0,
     &     1.50d0,1.49d0,1.38d0,1.36d0,1.26d0,
     &     1.20d0,1.16d0,1.14d0,1.06d0/

      if (n>109) stop 'Nuclear charge >109'
      radius=Radi(n)

      end function radius


C --- 2006.05.21 IDis Calc.
      subroutine NJ_idis(nat,link,IDis)
      implicit none
      integer nat,link(nat,nat),IDis(nat,nat),i,j,k,L,m
C
      IDis=100000
C
      do i=1,nat-1
         do j=i+1,nat
            if (link(i,j)==1) then
               IDis(i,j)=1
               IDis(j,i)=1
            endif
         enddo
      enddo
C
      do m=1,10
         do i=1,nat-1
            do j=i+1,nat
               do k=1,nat
                  if (k==i.or.k==j) cycle
                  L=IDis(i,k)+IDis(j,k)
                  if (L.ne.m) cycle
                  if (L<IDis(i,j)) then
                     IDis(i,j)=L
                     IDis(j,i)=L
                  endif
               enddo
            enddo
         enddo
      enddo
C
      do i=1,nat
         IDis(i,i)=0
      enddo
C
      end subroutine NJ_idis


C **** Construct atomic cluster for each atom 2015.2.4 Zhigang ****
      subroutine Atom_clu(NATOM,dis,atmclu,natmclu)
      implicit none
      integer::NATOM,i,j,thre,n,clu_last
      real(kind=8)::dis(NATOM,NATOM)
      integer::natmclu(NATOM,7)      ! num of atoms of each atomic cluster
      integer::atmclu(NATOM,NATOM,7) ! label of atoms in each atomic clu

      do thre=1,7                    ! thre=1~7 means distance=1.0~7.0
         do i=1,NATOM
            natmclu(i,thre)=1
            atmclu(1,i,thre)=i
            clu_last=NATOM
            do j=1,NATOM
               if (j==i) cycle
               if (dis(j,i)<(dble(thre))) then
                  natmclu(i,thre)=natmclu(i,thre)+1
                  n=natmclu(i,thre)
                  atmclu(n,i,thre)=j
               else
                  atmclu(clu_last,i,thre)=j
                  clu_last=clu_last-1
               end if
            end do
         end do
      end do

      end subroutine Atom_clu

C
C *************************************************************
C *  Calc Mulliken Population and The Atomic Labels for LMOs  *
C *  Modified by NZG_4/25/2016 @UARK                          *
C *************************************************************

      subroutine CalMPop(io,NATOM,nbas,nocc,nmo,SOVER,SMO,nbatm,
     &                   batom,SR,SR1,SOB,atmlevl,molevl)
      implicit none
      integer io,NATOM,nbas,nocc,Nmo,i,j,k,m,n,J2,J4,I2,kk,k1,k2,k3,k4
      integer nbatm(natom),batom(nbas,natom),SOB(Nmo,NATOM),n1,n2
      real(kind=8) SMO(nbas,nmo),SR(NATOM,Nmo),SR1(NATOM,Nmo)
      real(kind=8) SOVER(nbas,nbas) 
      real(kind=8) P1,P2,PP,NPP
      integer atmlevl(NATOM),molevl(nocc)
C
C----------------------------------------------------------------------------
C Calc Mulliken population (2001.04.21)
C SR(j,i)=\sum_{k}^{NW}\sum_{l\in{j}}C_{ki}S_{kl}C_{li}  (j=1,natom; i=1,Nmo)
C----------------------------------------------------------------------------
      write(6,*) "norb for pop",nmo
      DO I=1,NMO
         DO J=1,NATOM
            SOB(I,J)=J
            P1=0.0D0
            DO J2=1,nbatm(J)
               J4=BATOM(J2,J)
               P2=0.0D0
               DO K=1,nbas
                  P2=P2+SMO(K,I)*SOVER(K,J4)
               END DO
               P1=P1+SMO(J4,I)*P2
            END DO
            SR(J,I)=P1
         END DO
      END DO

C --- SR1: unsorted SR: Mulliken population of basis in atoms
      SR1=SR
C   
      if (io>0) WRITE(io,*) '+++ MULLIKEN GROSS POPULATION ON ATOMS +++'
      P2=0.0d0
      DO I=1,NATOM
         P1=0.0d0
         DO J=1,nocc
            P1=P1+SR(I,J)*2.0d0
         enddo
         if (io>0) then
            if (mod(I,3)==0 .or. I==NATOM) then
               WRITE(io,'(i6,f12.7)') I,P1
            else
               write(io,'(i6,f12.7)',advance='no') I,P1
            endif
         end if
         P2=P2+P1
      enddo
      if (io>0) then
         write(io,'('' Mulliken Population Total:'',f9.4)') P2
         write(io,*)
      end if
C
C-------------------------------------------------------------------------
C (SR(J,I),J=1,NATOM) are reordered from large to small for the ith Occ MO
C-------------------------------------------------------------------------
      DO I=1,nocc
         DO J=1,NATOM-1
            P1=dabs(SR(J,I))
            DO K=J+1,NATOM
               IF (dabs(SR(K,I)).GT.P1) THEN
                  P1=SR(K,I)
                  SR(K,I)=SR(J,I)
                  SR(J,I)=P1
                  P1=dabs(P1)
                  NPP=SOB(I,K)
                  SOB(I,K)=SOB(I,J)
                  SOB(I,J)=NPP
               ENDIF
            end do
         end do
      end do
C
      molevl=0
      do I=1,nocc
          k=SOB(I,1)
          k1=molevl(I)
          k2=atmlevl(k)
          if (k2.gt.k1) then
             molevl(I)=k2
          endif
      enddo
C
      if (io>0) then
         WRITE(io,*) 'Sorted Mulliken population of LMOs in atoms'
         call NJ_prtcol3(io,NATOM,NMO,SOB,SR,1,nocc,'f10.5')
      endif
  
      end subroutine CalMPop


      subroutine Fragment_Deter(io,nmo,ncore,FXY,SYS,NS0,NS1,
     &    molevl,thre,thre2,highorlow)
      implicit none
      integer::io,nmo,ncore,highorlow(nmo)
      integer::SYS(nmo-ncore,nmo-ncore),NS0(nmo-ncore),NS1(nmo-ncore)
      real(kind=8)::FXY(nmo-ncore,nmo-ncore),thre,thre2,thre0
      integer::i,j,k,l,tmp(nmo),KK,k1,k2,molevl(nmo),levl_low,levl_high
      integer,allocatable::Ktmp(:)
      character::line*1000
      logical multilevl
C
      multilevl = .false.
      highorlow = 0
      levl_low=100
      levl_high=0
      do i=ncore+1,nmo
         if (levl_low > molevl(i)) levl_low = molevl(i)
         if (levl_high< molevl(i)) levl_high= molevl(i)
      enddo
      if (levl_low/=levl_high) then
         do i=ncore+1,nmo
            if (molevl(i)==levl_high) highorlow(i)=1
         end do
      end if
      if (levl_low.ne.levl_high) then
         multilevl = .true.
         write(io,*) 'Lowest level  =',levl_low
         write(io,*) 'Highest level =',levl_high
      endif
C
      SYS(:,:)=0;NS0(:)=0;NS1(:)=0
      do i=1,nmo-ncore
         SYS(1,i)=i+ncore;NS0(i)=1;NS1(i)=1
      end do

      do i=1,nmo-ncore
         thre0 = thre
         if (multilevl.and.molevl(i+ncore)==levl_high) thre0=thre2
         do j=1,nmo-ncore
            if(j==i) cycle
            if(FXY(j,i)<=thre0) then
               NS1(i)=NS1(i)+1
               SYS(NS1(i),i)=j+ncore
           end if
         end do
C
         tmp(:)=0
         tmp(2:NS1(i))=SYS(2:NS1(i),i)
         do j=2,NS1(i)
            KK=tmp(2);l=2
            do k=2,NS1(i)
               if(tmp(k)>KK) then
                  KK=tmp(k);l=k
               end if
            end do
            SYS(NS1(i)-j+2,i)=KK;tmp(l)=0
         end do
      end do

      if (io>0) then
         write(io,'(1x,a)')'----- Original SYS -----'
         write(io,'(1x,a8,a8,a8)')'SYS','NUM','INDEX'
         do i=1,nmo-ncore
            allocate(Ktmp(NS1(i)))
            Ktmp(1:NS1(i))=SYS(1:NS1(i),i)
            call NJ_prtlab(line,NS1(i),Ktmp)
            deallocate(Ktmp)
            call NJ_trim(line,k1,k2)
            write(io,'(1x,I8,I8,1x,a)')i,NS1(i),line(k1:k2)
        end do
        write(io,'(1x,a)')'------------------------'
      end if

      return
      end subroutine Fragment_Deter


C *********************************************************************
C * Eliminate small clusters that are totally included in larger ones *
C * Always do not merge clusters with different levels                *
C * Modified by NZG_4/26/2016 @UARK                                   *
C *********************************************************************
      subroutine Reduce_Frag(io,nmo,ncore,Frag,NF1,NF2,nsys,molevl,
     &                       act_orb)
      implicit none
      integer::io,nmo,ncore,nsys
      integer::Frag(nmo-ncore,nmo-ncore),NF1(nmo-ncore),NF2(nmo-ncore)
      integer::i,j,k,NUW
      integer,allocatable::Sub1(:),Sub2(:),Sub3(:)
      integer,allocatable::SS1(:,:),NSS1(:),NSS2(:)
      integer::max2,L,mtmp,n,nfrg,isys,nfg,k1,k2,k3,k4,KK,molevl(nmo)
      integer:: L1,L2
      real(kind=8),allocatable::num1(:)
      integer,allocatable::Ktmp(:)
      character::line1*1000,line2*1000
      integer levl_low, levl_high,mrgsub
      logical multilevl,act_orb(nmo)
C
      NUW=nmo-ncore
C
      allocate(Sub1(NUW),Sub2(NUW),Sub3(NUW))
C   -----  Reduce the Redandunt SYS  ----
      do i=1,NUW-1
         do j=1+i,NUW
            if(NF2(j)>=NF2(i)) then
               n=NF2(i);NF2(i)=NF2(j);NF2(j)=n
               n=NF1(i);NF1(i)=NF1(j);NF1(j)=n
               Sub1=Frag(:,i);Frag(:,i)=Frag(:,j);Frag(:,j)=Sub1
            end if
         end do
      end do

      nsys=0
      do i=1,NUW-1
         if (NF2(i)==0) cycle
         do j=i+1,NUW
            if (NF2(j)==0) cycle
            L1 =  Frag(1,i)
            L2 =  Frag(1,j)
            if (molevl(L1).ne.molevl(L2)) cycle
            if (act_orb(L1).ne.act_orb(L2)) cycle
            Sub1(1:NUW)=Frag(1:NUW,i)
            Sub2(1:NUW)=Frag(1:NUW,j)
            call NJ_findover(0,NUW,Sub1,Sub2,Sub3,L)
            if (L==NF2(j)) then
               nsys=nsys+1
               l=0;NF1(i)=NF1(i)+1;Frag(NF1(i),i)=Frag(1,j)
C
               do k=NF2(i),NF1(i),-1
                  if(Sub1(k)==Frag(1,j)) cycle
                  l=l+1
                  Frag(NF2(i)-l+1,i)=Sub1(k)
               end do
               if ((NF2(i)-NF1(i)-l)/=0) then
                  write(io,*)'ERROR IN REDUCING FRAG',NF2(i),NF1(i),l
                  write(io,*)'i',i,(Sub1(k),k=1,NF2(i))
                  write(io,*)'j',j,(Sub2(k),k=1,NF2(j))
                  flush(io)
                  stop
                end if
                Frag(:,j)=0;NF2(j)=0;NF1(j)=0
                Sub1(:)=Frag(:,i)
            end if
         end do
      end do

      do i=1,NUW
         if(NF2(i)==0) cycle
         Sub1(:)=0
         Sub1(1:NF1(i))=Frag(1:NF1(i),i)
         do j=1,NF1(i)
            KK=Sub1(1)
            l=1
            do k=1,NF1(i)
               if (Sub1(k)>KK) then
                  KK=Sub1(k)
                  l=k
               end if
            end do
            Frag(NF1(i)-j+1,i)=KK
            Sub1(l)=0
         end do
      end do

      nsys=NUW-nsys

      nfrg=0
      allocate(SS1(NUW,nsys),NSS1(nsys),NSS2(nsys))
      SS1(:,:)=0;NSS1(:)=0
      do i=1,NUW
         if(NF2(i)==0) cycle
         nfrg=nfrg+1
         SS1(:,nfrg)=Frag(:,i);NSS1(nfrg)=NF1(i);NSS2(nfrg)=NF2(i)
      end do
C
      Frag(:,:)=0;NF1(:)=0;NF2(:)=0
      do i=1,nfrg
         Frag(:,i)=SS1(:,i);NF1(i)=NSS1(i);NF2(i)=NSS2(i)
      end do
      
      if (io>0) then
         if (nfrg==1) then
            write(io,*) "Warning: Only one cluster generated!"
            return
         end if
         write(io,'(1x,a)')'----------   Reduced cluster   ---------'
         write(io,'(1x,a3,a4,1x,a30,1x,a30)')
     &                         'CLU','NUM','Cent MOs','Buff MOs'
         do i=1,nsys
           allocate(Ktmp(NSS2(i)))
           Ktmp(1:NSS2(i))=SS1(1:NSS2(i),i)
           call NJ_prtlab(line1,NSS1(i),Ktmp(1:NSS1(i)))
           call NJ_prtlab(line2,NSS2(i)-NSS1(i),Ktmp(NSS1(i)+1:NSS2(i)))
           deallocate(Ktmp)
           call NJ_trim(line1,k1,k2)
           call NJ_trim(line2,k3,k4)
           write(io,'(1x,I3,I4,1x,a30,1x,a30)')
     &                               i,NSS2(i),line1(k1:k2),line2(k3:k4)

         end do
         write(io,'(1x,a)')'---------------------------------------'
      end if

      return
      end subroutine Reduce_Frag
C 
C 
      subroutine NJ_findover(io,nfg,Sub1,Sub2,Sub3,num)
      implicit none
      integer io,nfg,Sub1(nfg),Sub2(nfg),Sub3(nfg),num,i,j,k1,k2
      character line*100
C
      Sub3=0; num=0
      do i=1,nfg
         k1=Sub1(i)
         if (k1==0) exit
         do j=1,nfg
            k2=Sub2(j)
            if (k2==0) exit
            if (k2==k1) then
               num=num+1
               Sub3(num)=k1
            endif
         enddo
      enddo
C
      if (io>0.and.num>0) then
         call NJ_prtlab(line,nfg,Sub3)
         write(io,'(a)') trim(line)
      endif
      end subroutine NJ_findover


c     ##############################################################
c     ##  subroutine NJ_prtsym  --  print symmetric mat(n,n)      ##
c     ##  2005.12.22 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c
c     format: f(x.y) x>7 sugg 12.5,12.7,14.9
c
      subroutine NJ_prtsym(io,n,mat,fm)
      implicit none
      integer i,j,jj,n,io,n5,nf,nc,x,y,ini,ifi,k
      real*8 mat(n,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4

      n5=n/5; nf=mod(n,5)
      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y

      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
 101  format('(i7,5',a1,i2,'.',i2,')')
 102  format('(i7,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
 103  format('(3x,5(',i2,'x,i7))')
 104  format('(3x,',i2,'(',i2,'x,i7))')

      do jj=1,n5
         ini=1+(jj-1)*5
         write(io,fmt3) (j,j=ini,jj*5)
         do k=1+(jj-1)*5,n
            ifi=min(jj*5,k)
            write(io,fmt1) k,(mat(k,j),j=ini,ifi)
         enddo
      enddo

      if (nf.ne.0) then
         ini=n-nf+1
         write(io,fmt4)(j,j=ini,n)
         do k=ini,n
            write(io,fmt2) k,(mat(k,j),j=ini,k)
         enddo
      endif
      call flush(io)

      end subroutine NJ_prtsym


c     ##############################################################
c     ##  subroutine NJ_prtlab  --  write Ktmp(kk) into line      ##
c     ##  2005.01.09 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c
c     Format: "1,3,5-8,2*0,9"
c
      subroutine NJ_prtlab(line,kk,Ktmp)
      integer Ktmp(kk),kkk
      parameter (linemax=80)
      character ch,ch2,line*(*),line1*2000,line2*2000
      parameter (ch=',',ch2='-')

      line=' '; ini=1; fc=1
      if (kk==0) return

      write(line1,*) Ktmp(1)
      call NJ_trim(line1,k1,k2)
      line(ini:ini+k2-k1)=line1(k1:k2)
      ini=ini+k2-k1+1
      nz=0
C
      kkk=min(2,kk)
C
      if (Ktmp(1)==0.and.Ktmp(kkk)==0) then
         fc=0; nz=1; ini=1
      endif

      do 110 i=2,kk
         write(line1,*) Ktmp(i)
         call NJ_trim(line1,k1,k2)
         if (Ktmp(i)-Ktmp(i-1)==1) then
            if (i<kk) kkk=Ktmp(i+1)-Ktmp(i)
            if (i==kk.or.kkk.ne.1) then
               line(ini:ini+k2-k1+1)=ch2//line1(k1:k2)
               ini=ini+k2-k1+2
            endif
         elseif (Ktmp(i)==0) then
            nz=nz+1
            if (i==kk) then
               if (nz==1) then
                  line(ini:ini+1)=ch//'0'; ini=ini+2
               else
                  write(line2,*) nz; call NJ_trim(line2,k3,k4)
                  if (fc==1) then  ! 2005.01.09 add
                     line(ini:ini+k4-k3+3)=ch//line2(k3:k4)//'*0'
                     ini=ini+k4-k3+4; exit
                  else
                     line(ini:ini+k4-k3+2)=line2(k3:k4)//'*0'
                     ini=ini+k4-k3+3; fc=1; exit
                  endif
               endif
            else
               if (Ktmp(i+1).ne.0.and.nz==1) then
                  line(ini:ini+1)=ch//'0'; ini=ini+2; nz=0
               elseif (Ktmp(i+1).ne.0.and.nz.ne.1) then
                  write(line2,*) nz; call NJ_trim(line2,k3,k4)
                  if (fc==1) then  ! 2005.01.09 add
                     line(ini:ini+k4-k3+3)=ch//line2(k3:k4)//'*0'
                     ini=ini+k4-k3+4; nz=0
                  else
                     line(ini:ini+k4-k3+2)=line2(k3:k4)//'*0'
                     ini=ini+k4-k3+3; nz=0; fc=1
                  endif
               endif
            endif
         else
            line(ini:ini+k2-k1+1)=ch//line1(k1:k2)
            ini=ini+k2-k1+2
         endif
 110  enddo

      end subroutine NJ_prtlab

C
c     ##############################################################
c     ##  subroutine NJ_prtcol -- print column(c1:c2) of mat(m,n) ##
c     ##  2005.12.22 by Wei Li; Update 2005.12.25 by Wei Li       ##
c     ##############################################################
c
c     format: f(x.y) x>7 sugg 12.5,12.7,14.9
c
      subroutine NJ_prtcol2(io,m,n,mat,c1,c2,fm)
      implicit none
      integer i,j,jj,m,n,io,c1,c2,n5,nf,nc,x,y,k
      real(kind=8) mat(m,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4
C
      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
C
      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
 101  format('(i7,5',a1,i2,'.',i2,')')
 102  format('(i7,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
 103  format('(3x,5(',i2,'x,i7))')
 104  format('(3x,',i2,'(',i2,'x,i7))')
C
      do jj=1,n5
         write(io,fmt3) (j,j=c1+(jj-1)*5,c1+jj*5-1)
         write(io,fmt1) (i,(mat(i,j),j=c1+(jj-1)*5,c1+jj*5-1),i=1,m)
      enddo
C
      if (nf.ne.0) then
         write(io,fmt4)(j,j=c1+n5*5,c2)
         write(io,fmt2) (i,(mat(i,j),j=c1+n5*5,c2),i=1,m)
      endif
C
      end subroutine NJ_prtcol2
C
C
c     ##############################################################
c     ##  subroutine NJ_prtcol -- print column(c1:c2) of mat(m,n) ##
c     ##  2005.12.22 by Wei Li; Update 2005.12.25 by Wei Li       ##
c     ##############################################################
c
c     format: f(x.y) x>7 sugg 12.5,12.7,14.9
c
      subroutine NJ_prtcol3(io,m,n,mat0,mat,c1,c2,fm)
      implicit none
      integer i,j,jj,m,n,io,c1,c2,n5,nf,nc,x,y,k,mat0(n,m)
      real(kind=8) mat(m,n)
      character fm*(*),ch,fm2*10
      character*40 fmt1,fmt2,fmt3,fmt4
C
      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
C
      write(fmt1,101) ch,x,y; write(fmt2,102) nf,nf,ch,x,y
 101  format('(i4,1x,5i4,5',a1,i2,'.',i2,')')
 102  format('(i4,1x,',i2,'i4,',i2,a1,i2,'.',i2,')')
      write(fmt3,103) x-7; write(fmt4,104) nf,nf,x-7
 103  format('(5x,5i4,5(',i2,'x,i7))')
 104  format('(5x,',i2,'i4,',i2,'(',i2,'x,i7))')
C
      do jj=1,n5
         write(io,fmt3) (j,j=c1+(jj-1)*5,c1+jj*5-1),
     &                  (j,j=c1+(jj-1)*5,c1+jj*5-1)
         write(io,fmt1)
     &  (i,(mat0(j,i),j=c1+(jj-1)*5,c1+jj*5-1),
     &      (mat(i,j),j=c1+(jj-1)*5,c1+jj*5-1),i=1,m)
         if (jj.ne.n5.or.nf.ne.0) write(io,'(1x,74(''-''))')
      enddo
C
      if (nf.ne.0) then
         write(io,fmt4)(j,j=c1+n5*5,c2),(j,j=c1+n5*5,c2)
         write(io,fmt2) (i,(mat0(j,i),j=c1+n5*5,c2),
     &                  (mat(i,j),j=c1+n5*5,c2),i=1,m)
      endif
      call flush(io)
C
      end subroutine NJ_prtcol3

C
c     ##############################################################
c     ##  subroutine NJ_trim  --  move blank of two sides         ##
c     ##  2005.01.07 by Wei Li; Update 2005.11.01 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_trim(line,k1,k2)
      implicit none
      integer k1,k2,i,j
      character line*(*)

      j=len(line)
      do i=1,j
         if (line(i:i).ne.' ') then
            k1=i; exit
         endif
         if (i==j) then
            k1=1; k2=1
            return
         endif
      enddo

      do i=j,1,-1
         if (line(i:i).ne.' ') then
            k2=i; exit
         endif
      enddo

      end subroutine NJ_trim


C --- Count the number of frozen core orbitals automatically
      subroutine FRZORB(io,NATOM,IAN,nfocc)
      implicit double precision(A-H,O-Z)
      dimension IAN(NATOM)
C
      nfocc=0
      do i=1,NATOM
         if (IAN(i).gt.86) then
            nfocc=nfocc+43
         elseif (IAN(i).gt.54) then
            nfocc=nfocc+27
         elseif (IAN(i).gt.36) then
            nfocc=nfocc+18
         elseif (IAN(i).gt.18) then
            nfocc=nfocc+9
         elseif (IAN(i).gt.10) then
            nfocc=nfocc+5
         elseif (IAN(i).gt.2) then
            nfocc=nfocc+1
         endif
      enddo
C
      end subroutine FRZORB
C
C --- Normalization of a set of orbital mo(nbs,nmo)
      subroutine normorb(nbs,nmo,mo,s)
      implicit none
      integer nbs,nmo,i,j,k,L,m,n
      real(kind=8) mo(nbs,nmo),s(nbs,nbs),P1,PP
C
      do i=1,nmo
         call normfact(nbs,mo(1,i),s,P1)
         do j=1,nbs
            mo(j,i)=mo(j,i)/P1
         enddo
      enddo
C
      end subroutine normorb
C

C --- Determine the normalization factor of one orbital c(nbs), vv returned
C     For normali,ization use c=c/p
      subroutine normfact(nbs,ac,s,p)
      implicit none
      integer nbs,k1,k2
      real(kind=8) ac(nbs),s(nbs,nbs),p
C
      p=0.0d0
      do k1=1,nbs
         do k2=1,nbs
            p=p+ac(k1)*ac(k2)*s(k2,k1)
         end do
      end do
      if (p.le.0.0d0) stop 'Error in normvalu(): negotive module'
      p=dsqrt(p)
      end subroutine normfact


      subroutine NJ_denmat(io,nbs,nif,noc,cmo,dm)
      implicit none
      integer io,nbs,nif,noc,i,j,k
      real*8 cmo(nbs,nif),dm(nbs,nbs),PP
      real*8,allocatable::cmo1(:,:)
      allocate(cmo1(noc,nbs))
      cmo1=transpose(cmo(:,1:noc))
      dm=0.0d0
      do i=1,nbs; do j=1,nbs
         PP=0.0d0
         do k=1,noc
            PP=PP+cmo1(k,j)*cmo1(k,i)
         enddo
         dm(j,i)=2.0d0*PP
      enddo; enddo
      deallocate(cmo1)

      if (io>0) then
         write(io,*) '*** Density matrix ***'
         call NJ_prtsym(io,nbs,dm,'d14.6')
         write(io,*)
      endif
      return
      end


      subroutine NJ_denmat_sym(nbs,nif,noc,cmo,dm)
      implicit none
      integer nbs,nif,noc,i,j,k,ij
      real*8 cmo(nbs,nif),dm(nbs*(nbs+1)/2),PP

      dm=0.0d0; ij=0
      do i=1,nbs
         do j=1,i
            ij=ij+1
            PP=0.0d0
            do k=1,noc
               PP=PP+cmo(i,k)*cmo(j,k)
            enddo
            dm(ij)=2.0d0*PP
         enddo
      enddo

      return
      end subroutine NJ_denmat_sym


      subroutine inarray(element,array,num,init)
      implicit none
      integer::array(num),element,num,i,init
     
      do i=1,num
         if (element==array(i)) then
            init=1
            return
         else
            init=0  
         end if
      end do
      return
      end subroutine inarray


      subroutine dunit(n,a)
      real(kind=8) a(n,n)
      do i=1,n
         do j=1,n
            if (j.eq.i) then
               a(j,i)=1.0d0
            else
               a(j,i)=0.0d0
            endif
         enddo
      enddo
      end subroutine dunit


C --- 2008.02.21 --- Project out the mo1(nbs,nmo1) from mo2(nbs,nmo2)
      subroutine pjotorb(nbs,nmo1,nmo2,mo1,mo2,s)
      implicit none
      integer io,nbs,nmo1,nmo2,i,j,k,L,m,n
      real(kind=8) mo1(nbs,nmo1),mo2(nbs,nmo2),s(nbs,nbs)
      real(kind=8),allocatable::p(:,:)
      real(kind=8)::P1,P2
C
      allocate(p(nbs,nbs))
C
      call HF_DenMat(nbs,nmo1,mo1,p)
      P1=-1.0D0; P2=1.0D0
      call DGEMM('N','N',nbs,nbs,nbs,P1,p,nbs,s,nbs,P2,mo2,nbs)
C
      deallocate(p)
C
      end subroutine pjotorb

      subroutine HF_DenMat(nbs,noc,smo,dm)
      implicit none
      integer io,nbs,noc
      real*8 smo(nbs,noc),dm(nbs,nbs),P1,P2
C
      P1=1.0D0; p2=0.0D0
      call DGEMM('N','t',nbs,nbs,noc,P1,smo,nbs,smo,nbs,P2,dm,nbs)
C
      return
      end subroutine HF_DenMat


      subroutine matcopy_cim(a,b,m,n)
      implicit none
      integer::m,n
      real*8 a(m,n),b(m,n)

      b=a
      end subroutine matcopy_cim



C ****************************************************
C * Split the string by comma and return each number *
C * NZG_11/20/2016 @NJU                              *
C ****************************************************
      subroutine get_act_atm(act_str,act_atm,NAtom)
      implicit none
      
      integer i,j,k,NAtom,length,actnum,asciinum,totlen
      character*256 act_str
      logical act_atm(NAtom)

      length=0
      totlen=len_trim(act_str)
      do i=1,totlen
         asciinum=ichar(act_str(i:i))
         if (asciinum>=48 .and. asciinum<=57) then
            length=length+1
         else
            select case(length)
            case(4)
               actnum=(ichar(act_str(i-1:i-1))-48)+
     &                (ichar(act_str(i-2:i-2))-48)*10+
     &                (ichar(act_str(i-3:i-3))-48)*100+
     &                (ichar(act_str(i-4:i-4))-48)*1000
            case(3)
               actnum=(ichar(act_str(i-1:i-1))-48)+
     &                (ichar(act_str(i-2:i-2))-48)*10+
     &                (ichar(act_str(i-3:i-3))-48)*100
            case(2)
               actnum=(ichar(act_str(i-1:i-1))-48)+
     &                (ichar(act_str(i-2:i-2))-48)*10
            case(1)
               actnum=ichar(act_str(i-1:i-1))-48
            case(0)
               write(6,*) "Please Check the Active Space Input!"
               stop
            case default
               write(6,*) "Atomic Number Over 9999!"
               stop
            end select
            act_atm(actnum)=.true.
            length=0
         endif
      enddo

      i=totlen+1
      select case(length)
      case(4)
         actnum=(ichar(act_str(i-1:i-1))-48)+
     &          (ichar(act_str(i-2:i-2))-48)*10+
     &          (ichar(act_str(i-3:i-3))-48)*100+
     &          (ichar(act_str(i-4:i-4))-48)*1000
      case(3)
         actnum=(ichar(act_str(i-1:i-1))-48)+
     &          (ichar(act_str(i-2:i-2))-48)*10+
     &          (ichar(act_str(i-3:i-3))-48)*100
      case(2)
         actnum=(ichar(act_str(i-1:i-1))-48)+
     &          (ichar(act_str(i-2:i-2))-48)*10
      case(1)
         actnum=ichar(act_str(i-1:i-1))-48
      case(0)
         return
      case default
         write(6,*) "Atomic Number Over 9999!"
         stop
      end select
      act_atm(actnum)=.true.
      
      end subroutine get_act_atm


C *********************************************************************
C * Calculate the overlap of MOs to check if they are orthonormalized *
C * Zhigang_4/9/2017 @UARK                                            *
C *********************************************************************

      subroutine MO_over(Sij,C,S,nmo,nbas)
      implicit none
      
      integer nmo,nbas,i,j,k,l
      real*8 Sij(nmo,nmo),C(nbas,nmo),S(nbas,nbas)

      Sij=0.0D0
      do j=1,nmo
         do i=1,nmo
            do k=1,nbas
               do l=1,nbas
                  Sij(i,j)=Sij(i,j)+C(k,i)*S(k,l)*C(l,j)
               enddo
            enddo
         enddo
      enddo

      return
      end subroutine MO_over


      subroutine checkiden(m,n,a,b,same)
      implicit none
      
      integer m,n,i,j,same
      real*8 a(m,n),b(m,n)

      do i=1,n
         do j=1,m
            if (abs(a(i,j)-b(i,j))>1.0d06) then
               same=0
               return
            endif
         enddo
      enddo
      same=1
      return
      end subroutine checkiden
