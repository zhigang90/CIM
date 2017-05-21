      subroutine formgmat(nq,nek,b,ibcontr,gmat)
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              OUTPUT
c              gmat(*) G matrix (BBt) in square form; note that Bt is stored in the program for some reason
      implicit real*8 (a-h,o-z)
      parameter(one=1.0d0,eps=1.0d-9)
      dimension b(54,nq),ibcontr(20,nq),gmat(*)
c       write(*,*) 'ibcontr'
c      do nn=1,6
c         write(*,60) (ibcontr(nn,kq),kq=1,nq)
c 60   format((5i5,2x,5i5))
c      end do
      ij=0
      do i=1,nq
        natomi=ibcontr(2,i)
        do j=1,nq
          sum=0.0d0
          natomj=ibcontr(2,j)
          do njj=1,natomj
c  This is the true number of the atom. However, here we keep only the
c  non-zero elements
            jatom=ibcontr(njj+2,j)
            do nii=1,natomi
              iatom=ibcontr(nii+2,i)
              if(iatom.eq.jatom) then
                i3=3*nii
                j3=3*njj
c      write(*,5) i,j,ibcontr(nii+2,i),jatom,(b(i3-k,i),b(j3-k,j),k=-2,0)
    5 format(2i3,2x,2i3,2x,6f12.6)
                sum=sum+b(i3-2,i)*b(j3-2,j)+b(i3-1,i)*b(j3-1,j)+b(i3,i)*
     &                  b(j3,j)  
              end if
            end do
          end do
          ij=ij+1
          gmat(ij)=sum
c          print *, 'GMAT, i,j=',i,j,sum
        end do
      end do 
c      call prntmat(nq,nq,nq,gmat)
      end
c======================================================================
      subroutine geninv (nq,c,d,tol,u,    x,y,ninv,proj)
C     GENERALIZED INVERSE ROUTINE 
C     DIAGONALIZE, INVERT IF THE DIAGONAL IS NOT TOO SMALL,
C     INVERT BACK
c     nq=dimension of the square matrix c
c     c= matrix to get the gen,. inverse. It is replaced by its gen. inverse 
c      nq=dimension of it
c     d is the determinant - no function here
c     tol a small number, if an eigenvalue of c is less than tol 
c          then this root is not inverted
c     u holds the eigenvectors of C
c     x and y are two temporary vectors nq long. x will hold the eigenvalues 
c     ninv is a result: the number of non-zero (not small) eigenvalues
c     ginv gives the projector to the (numerically) non-null eigenspace of c
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension c(nq,nq),u(nq,nq),x(*),y(*),proj(nq,nq)
      include "common/big"
c
      call sdiag2(nq,nq,c,x,u)
c
      det =1.0d0
      ninv=0  
c  eigenvalues in ascending order
      ifirst=0
      do 100 i=1,nq
         if(abs(x(i)).gt.tol) then
           if(ifirst.eq.0) ifirst=i
           det=det*x(i)
           x(i)=1.0d0/x(i)
           ninv=ninv+1
         else
           x(i)=0.0d0
         end if  
      d=det   
 100  continue
c      write(*,*) 'ifirst, Eigenvalues'
c      write(*,150) ifirst,(x(k),k=1,nq) 
 150  format(i5,2x,5f12.6,/,(7x,5f12.6))
c       write(*,*) 'Eigenvectors of G'
c       call prntmatn(6,nq,nq,nq,u,'G eigenv.')
c       call flush6()
c       call outsqr(u,nq,nq,nq,'Eigenvectors of G')
      do 400 i=1,nq
        do 300 j=1,nq
          sum=0.0d0
          sum1=0.0d0
          do 200 k=ifirst,nq
            sum=sum+u(i,k)*u(j,k)*x(k)
            sum1=sum1+u(i,k)*u(j,k)
 200        continue
           c(i,j)=sum
           proj(i,j)=sum1
 300     continue
 400     continue
c
c        write(*,*) 'generalized inverse of G in geninv'
c        call prntmatn(6,nq,nq,nq,c,'G-')
c        write(*,*) 'projector to the nonzero subspace of G'
c        call prntmatn(6,nq,nq,nq,proj,'P')
c       call flush6()
c       call outsqr(c,nq,nq,nq,'generalized inverse of G')
c       call flush6()
      end
c======================================================================
      subroutine ginvxbf(nq,ginv,bf,phi)
c   This routine simply performs the matrix-vector product G-*[Bf]
c   where G- = the generalized inverse of G, [Bf]=B*f, B= B mtx, f= Cart. forces
c   It should be done, of course, by LAPACK routines
c   The result is the internal forces in phi
      implicit real*8 (a-h,o-z)
      dimension ginv(nq,nq),bf(nq),phi(nq)
      do i=1,nq
        sum=0.0d0
        do j=1,nq
          sum=sum+ginv(i,j)*bf(j)
        end do
        phi(i)=sum
      end do
      end
c======================================================================
      subroutine matxdiag(n,a,d)
c   Forms the product A*D where D is diagonal
c   The result is returned in A
      implicit real*8 (a-h,o-z)
      dimension a(n,n),d(n)
      do j=1,n
        do i=1,n
          a(i,j)=a(i,j)*d(j)
        end do
      end do
      end
c======================================================================
      subroutine transfdfc(nq,d,hess)
      implicit real*8 (a-h,o-z)
c  This routine simply pots the elements of D, interpreted as a diagonal matrix,
c  in the diagonal of Hess
c
c  Parameters INTENT(IN)
c  nq = dimension of the square matrices
c  d  = the diagonal matrix
c  hess=the force constant matrix. It will be a diagonal mtx hess9i,i)=d(i)
c
      dimension d(nq),hess(nq,nq)
      call zeroit(hess,nq**2)
      do i=1,nq
        hess(i,i)=d(i)
      end do
      end
c======================================================================

      SUBROUTINE PrntMATN(iout,N,NRow,NCol,A,title)
      REAL*8 A(NRow,NCol)
      character*(*) title
C
C  Prints out N columns of an NRow * NCol matrix A
C
      PARAMETER (maxcol=10)
C
      write(iout,*) title
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
 1100 FORMAT(1X,10F12.6)
c
      END     
c
c ========================================================================

