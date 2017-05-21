c  fixed dimension in <SDIAG2> now dynamically allocated
c  JB   Oct 18 2000
c ==================================================================
c May 6, 98 PP: I have modified subroutine nerror to print out
c   the memory usage etc.
c ==================================================================
c  fixed dimension (formerly 1501) increased to 3000 in <SDIAG2>
c  JB   July 25 1997
c ==================================================================
c
c     this module contains the home-made service routines
c     like matrix manipulations etc.
c     some of these duplicate blas routines and should be replaced
c     by them later
c.... now come the homebrewn
c..  replaced by blas routine (same calling sequence)
c     subroutine add1 (a,con,b,nn)
c     implicit real*8 (a-h,o-z)
c     dimension a(1), b(1)
c     n=nn
c     do 10 i=1,n
c  10 b(i)=b(i)+con*a(i)
c     return
c
c     end
c     subroutine tfer (a,b,n)
c     implicit real*8 (a-h,o-z)
c     dimension a(1), b(1)
c
c           transfer a to b
c
c     nn=n
c     do 10 i=1,nn
c  10 b(i)=a(i)
c     return
c
c     end
c=============================================================
      subroutine add(a,b,n)
      implicit double precision (a-h,o-z)
      dimension a(n),b(n)
c...  adds a to b, result in b
c%include '/sys/ins/base.ins.ftn'
c%include '/sys/ins/vec.ins.ftn'
c      call vec_$dadd_vector(a,b,n,b)
      do 100 i=1,n
        b(i)=b(i)+a(i)
 100  continue
      end
c==============================================================
      subroutine add1(a,con,b,n)
      implicit double precision (a-h,o-z)
      dimension a(n),b(n)
      call daxpy(n,con,a,1,b,1)
      end
c==============================================================
      subroutine tfer (a,b,n)
      implicit double precision (a-h,o-z)
      dimension a(n),b(n)
      call dcopy(n,a,1,b,1)
      end
c==============================================================
      subroutine mxmtr(a,b,c,m,n,p)
c...  this is not a true blas routine- it is our own creation
c...  this subroutine performs the matrix multiply c=a(trspse)*b
c     both the storage dimension and the real dimension of the
c     matrices is given by a(p,m),b(p,n),c(m,n)
c     it uses an efficient blocking algorithm
c     note that in general it is more efficient to perform a(tr)b
c     than a*b. transpose the matrix before if you have to
      implicit real*8 (a-h,o-z)
      integer p
      dimension  a(p,m), b(p,n),c(m,n)
      parameter (ibl=16)
c --  performs a(tr)*b=c
      do 500 i=1,m,ibl
        iend=min0(i+ibl-1,m)
        do 500 j=1,n,ibl
          jend=min0(j+ibl-1,n)
          do 500 i1=i,iend
          do 500 j1=j,jend
          sum=0.0d0
          sum=ddot(p,a(1,i1),1,b(1,j1),1)
          c(i1,j1)=sum
 500  continue
      return
      end
c==============================================================
      subroutine sdiag2 (m,n,a,d,x)
      implicit real*8 (a-h,o-z)
c
c
c      computation of all eigenvalues and eigenvectors of a real
c      symmetric matrix by the method of qr transformations.
c      if the euclidean norm of the rows varies   s t r o n g l y
c      most accurate results may be obtained by permuting rows and
c      columns to give an arrangement with increasing norms of rows.
c
c      two machine constants must be adjusted appropriately,
c      eps = minimum of all x such that 1+x is greater than 1 on the
c      e     computer,
c      tol = inf / eps  with inf = minimum of all positive x represen-
c            table within the computer.
c
c      input
c
c      (m)   corresponding value of the actual
c            dimension statement a(m,m), d(m), x(m,m),
c      (n)   not larger than (m), order of the matrix,
c      (a)   the matrix to be diagonalized, its lower triangle has to
c            be given as  ((a(i,j), j=1,i), i=1,n),
c
c      output
c
c      (d)   components d(1), ..., d(n) hold the computed eigenvalues
c            in ascending sequence. the remaining components of (d) are
c            unchanged,
c      (x)   the computed eigenvector corresponding to the j-th eigen-
c            value is stored as column (x(i,j), i=1,n). the eigenvectors
c            are normalized and orthogonal to working accuracy. the
c            remaining entries of (x) are unchanged.
c
c      array (a) is unaltered. however, the actual parameters
c      corresponding to (a) and (x)  may be identical, ''overwriting''
c      the eigenvectors on (a).
c
c      leibniz-rechenzentrum, munich 1965
c
c
      dimension a(m,m), d(m), x(m,m)
      dimension e(m)    ! dynamically allocated
c
c     correct adjustment for ieee 64  bit floating  point numbers`
c       is about (not too sensitive)
      eps=1.7764d-15
      tol=2.782d-307
cc      eps=1.0d-14
cc      tol=1.0d- 200
c
      if (n.eq.1) go to 350
      do 10 i=1,n
      do 10 j=1,i
   10 x(i,j)=a(i,j)
c
c     householder's reduction
c     simulation of loop do 150 i=n,2,(-1)
c
      do 150 ni=2,n
         ii=n+2-ni
c
c     fake loop for recursive address calculation
c
      do 150 i=ii,ii
         l=i-2
         h=0.0d0
         g=x(i,i-1)
         if (l) 140,140,20
   20    do 30 k=1,l
   30    h=h+x(i,k)**2
         s=h+g*g
         if (s.ge.tol) go to 40
         h=0.0d0
         go to 140
   40    if (h) 140,140,50
   50    l=l+1
         f=g
         g=sqrt(s)
         if (f) 70,70,60
   60    g=-g
   70    h=s-f*g
         x(i,i-1)=f-g
         f=0.0d0
c
         do 110 j=1,l
            x(j,i)=x(i,j)/h
            s=0.0d0
            s=ddot(j,x(j,1),m,x(i,1),m)
            j1=j+1
            if (j1.gt.l) go to 100
            j1l=l-j1+1
            s=s+ddot(j1l,x(j1,j),1,x(i,j1),m)
  100       e(j)=s/h
  110    f=f+s*x(j,i)
c
         f=f/(h+h)
         call daxpy(l,-f,x(i,1),m,e,1)
c
         do 130 j=1,l
            f=x(i,j)
            s=e(j)
         call daxpy(j,-f,e,1,x(j,1),m)
         call daxpy(j,-s,x(i,1),m,x(j,1),m)
  130    continue
c
  140    d(i)=h
  150 e(i-1)=g
c
c     accumulation of transformation matrices
c
      d(1)=x(1,1)
      x(1,1)=1.0d0
      do 200 i=2,n
         l=i-1
         if (d(i)) 190,190,160
  160    do 180 j=1,l
            s=0.0d0
            s=ddot(l,x(i,1),m,x(1,j),1)
       call daxpy(l,-s,x(1,i),1,x(1,j),1)
  180  continue
  190    d(i)=x(i,i)
         x(i,i)=1.0d0
      do 200 j=1,l
         x(i,j)=0.0d0
  200 x(j,i)=0.0d0
c
c     diagonalization of the tridiagonal matrix
c
      b=0.0d0
      f=0.0d0
      e(n)=0.0d0
c
      do 310 l=1,n
         h=eps*(abs(d(l))+abs(e(l)))
         if (h.gt.b) b=h
c
c     test for splitting
c
         do 210 j=l,n
            if (abs(e(j)).le.b) go to 220
  210    continue
c
c     test for convergence
c
  220    if (j.eq.l) go to 310
c
c     shift from upper 2*2 minor
c
  230    p=(d(l+1)-d(l))*0.5/e(l)
         r=sqrt(p*p+1.0)
         if (p) 240,250,250
  240    p=p-r
         go to 260
  250    p=p+r
  260    h=d(l)-e(l)/p
         do 270 i=l,n
  270    d(i)=d(i)-h
         f=f+h
c
c     qr transformation
c
         p=d(j)
         c=1.0
         s=0.0
c
c     simulation of loop do 330 i=j-1,l,(-1)
c
         j1=j-1
         do 300 ni=l,j1
            ii=l+j1-ni
c
c     fake loop for recursive address calculation
c
         do 300 i=ii,ii
            g=c*e(i)
            h=c*p
c
c     protection against underflow of exponents
c
            if (abs(p).lt.abs(e(i))) go to 280
            c=e(i)/p
            r=sqrt(c*c+1.0)
            e(i+1)=s*p*r
            s=c/r
            c=1.0/r
            go to 290
  280       c=p/e(i)
            r=sqrt(c*c+1.0)
            e(i+1)=s*e(i)*r
            s=1.0/r
            c=c/r
  290       p=c*d(i)-s*g
            d(i+1)=h+s*(c*g+s*d(i))
         call drot(n,x(1,i+1),1,x(1,i),1,c,s)
  300    continue
c
         e(l)=s*p
         d(l)=c*p
         if (abs(e(l)).gt.b) go to 230
c
c     convergence
c
  310 d(l)=d(l)+f
c
c     ordering of eigenvalues
c
      ni=n-1
      do 340 i=1,ni
         k=i
         p=d(i)
         j1=i+1
         do 320 j=j1,n
            if (d(j).ge.p) go to 320
            k=j
            p=d(j)
  320    continue
         if (k.eq.i) go to 340
         d(k)=d(i)
         d(i)=p
         do 330 j=1,n
            p=x(j,i)
            x(j,i)=x(j,k)
  330    x(j,k)=p
  340 continue
      go to 360
c
c     special treatment of case n = 1
c
  350 d(1)=a(1,1)
      x(1,1)=1.0
  360 return
c
      end
c==================================================================
      subroutine zeroit (a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
c     call vec_$dzero (a,n)
      do  100 i=1,n
        a(i)=0.0d0
 100  continue
      return
      end
c=============================================================
      subroutine izeroit (a,n)
      implicit integer(a-z)
      dimension a(n)
      do 100 i=1,n
       a(i) = 0
 100  continue
      return
      end
c==============================================================
      subroutine mult(tmat,con,num)
      implicit real*8 (a-h,o-z)
      dimension tmat(num)
c
c     this subroutine multiplies an array by a constant
c     this is  jmb's invention; it is the same as dscal
      call  dscal(num,con,tmat,1)
c
      end
c==============================================================
      subroutine mave(amat,vold,vnew,num)
      implicit real*8 (a-h,o-z)
      dimension amat(*),vold(num),vnew(num)
c
c     this subroutine multiplies a triangular matrix by a column vector
c     and returns a new column vector.  vnew can't be same place as vold
c
      n=num
      do 13 i=1,n
      ij=i*(i-1)/2
      sum=0.0
         do 20 j=1,n
         ij=ij+1
         if(j.gt.i) ij=ij+j-2
         sum=sum+amat(ij)*vold(j)
   20 continue
      vnew(i)=sum
   13 continue
      return
c
      end
c==============================================================
      subroutine tri(a,b,m)
      implicit real*8 (a-h,o-z)
      dimension a(*), b(*)
c
c     a symmetric quadratic matrix a is changed to triangular form b
c
      mm=m
      ij=0
      do 200 i=1,mm
         iad=i*mm-mm
         do 100 j=1,i
            ij=ij+1
            b(ij)=a(iad+j)
  100 continue
  200 continue
      return
c
      end
c==============================================================
      subroutine quad(a,b,c,m)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(m,m)
c
c     make a quadratic symmetric or antisymmetric matrix b from
c     triangular matrix a
c
      c1=c
      con=abs(c1)
      ij=0
      do 10 i=1,m
         do 20 j=1,i
            ij=ij+1
            b(i,j)=c1*a(ij)
            b(j,i)=con*a(ij)
   20 continue
      b(i,i)=b(i,i)*(c1+con)/2
   10 continue
      return
c
      end
c==============================================================
      subroutine vecmat(amat,istcol,icol,vecin,vecout,idime)
      implicit real*8 (a-h,o-z)
      dimension amat(1), vecin(1), vecout(1)
c
c     multiplies amat(t) times a vector starting from stcol to stcol+ind
c     looks like josep bofill's creation
      ist=istcol
      ic=icol
      id=idime
      indcol=ist+ic-1
      do 10 i=ist,indcol
         nad=i*id-id
         sum=0.0
         do 20 j=1,id
            sum=sum+amat(j+nad)*vecin(j)
   20 continue
         vecout(i)=sum
   10 continue
      return
c
      end
c==============================================================
      subroutine matmat(a,ias,ia,b,ibs,ib,d,id)
      implicit real*8 (a-h,o-z)
      dimension a(1),b(1),d(1)
c
c     multiplies a(t) x b with common dimension id, result in d
c
      ja=ia
      jas=ias
      jb=ib
      jbs=ibs
      jae=ja+jas-1
      jbe=jb+jbs-1
      idim=id
      do 30 i=jas,jae
         iad=i*idim-idim
         do 20 j=jbs,jbe
            jad=j*idim-idim
            sum=0.0d0
            do 10 k=1,idim
               sum=sum+a(k+iad)*b(k+jad)
   10 continue
            d(jad+i)=sum
   20 continue
   30 continue
      return
c
      end
c==============================================================
      subroutine eig (b,a,d,n)

      use memory

      implicit real*8 (a-h,o-z)
c     common /big/bl(1000)
      dimension a(n,n), b(1), d(1)
cc      ij=0
cc      do 10 i=1,n
cc      do 10 j=1,i
cc         ij=ij+1
cc         a(i,j)=b(ij)
cc   10 continue
cc      call sdiag2 (n,n,a,d,a)
c
      call getrval('zero',zero)
      call getmem(8*n,iadrtmp)
      call getmem(5*n,iwork)
      call getmem(n,iadrfail)
      call dspevx('V','A','U',n,b,zero,zero,1,n,
     1            0.0d0,m,d,a,n,bl(iadrtmp),
     2            bl(iwork),bl(iadrfail),info)
      call retmem(3)
c this is the diagonalizer used in mindo, LAPACK is a lot faster
c dspevd is the fastest and best in parallel, but uses a lot of scratch memory
c dspevx is also fast and does not need much scratch
c in the guess usually there is plenty of memory
cc      iw1=1+5*n+2*n*int(log(dble(n))/log(2.0d0)+1)+2*n**2
cc      iw2=2+5*n
cc      call getmem(iw1,itmp1)
cc      call getmem(iw2,itmp2)
cc      call dspevd('V','U',n,b,d,a,n,bl(itmp1),iw1,bl(itmp2),iw2,IErr)
cc      call retmem(2)
       end
c==============================================================
      subroutine eig1(b,a,d,n,neig)

      use memory

      implicit real*8 (a-h,o-z)
c  This routine calculates neig eigenvalues of the matrix b
c  n is the order of the matrix
c  b contains the matrix, in packed upper triangular form
c  a holds, on exit, the eigenvalues
c     common /big/bl(1000)
      dimension a(n,n), b(1), d(1)
      call getrval('zero',zero)
      call getmem(8*n,iadrtmp)
      call getmem(5*n,iwork)
      call getmem(n,iadrfail)
      abstol=2*Dlamch('S')
      call dspevx('V','I','U',n,b,zero,zero,1,neig,abstol,
     1            m,d,a,n,bl(iadrtmp),bl(iwork),bl(iadrfail),info)
      call retmem(3)
      end
c=============================================================================
      subroutine spur (A,B,n,s)
c  This subroutine calculates the trace of the product of two symmetrical
c  matrices, stored in triangular form (upper triangle column-wise)
c  Arguments
c  INTENT(IN)
c  A, B = two symmatrical matrices, stored as n*(n+1)/2 long arrays
c         indexing: (ij)=i*(i-1)/2+j where i>=j
c  n    = dimension
c  INTENT(OUT)
c  s = Trace(A*B)
c  Spur is German for Trace
c
      implicit real*8 (a-h,o-z)
      parameter(two=2.0d0)
      dimension A(*), B(*)
      ntri=n*(n+1)/2
      s=ddot(ntri,A,1,B,1)*two
      ii=0
      do i=1,n
        ii=ii+i
        s=s-A(ii)*B(ii)
      end do
      end
c      s=zero
c      ij=0
c      do 10 i=1,ncf
c      do 10 j=1,i
c         ij=ij+1
c         s=s+a(ij)*b(ij)
c         if (i.eq.j) s=s-a(ij)*b(ij)*half
c   10 continue
c      s=s*two
c      return
cc
c      end
c==============================================================
      subroutine trans (a,ncf)
      implicit real*8 (a-h,o-z)
      dimension a(ncf,ncf)
      do 10 i=1,ncf
      do 10 j=1,i
         s=a(i,j)
         a(i,j)=a(j,i)
         a(j,i)=s
   10 continue
      return
c
      end
c======================================================================
c linking routine between mamu and mxma which translates existing calls
c of mamu directly into calls of mxma.
c richard g.a. bone, march, 1991.
c main interpretative problem is that mamu performs b*a matrix multiplication
c whereas mxma performs a*b.
c mamu works with either square or triangular matrices according to
c arguments k,l,m, whereas mxma requires squares only.
c n.b. in pulay's group mxma may not take the same array address twice.
c
      subroutine mamu (b,a,c,k,l,m,n,w)

      use memory

      implicit double precision (a-h,o-z)
      dimension a(1),b(1),c(1),w(1)
c     common /big/bl(1)
      nsq = n**2
c
      if((k.ne.0).and.(k.ne.1)) goto 9998
      if((l.ne.0).and.(l.ne.1)) goto 9998
      if((m.ne.0).and.(m.ne.1)) goto 9998
c
      call getmem(nsq,ia)
      call getmem(nsq,ic)
c
      if (l.eq.0) call squr(a,bl(ia),n)
      if (l.eq.1) call tfer(a,bl(ia),nsq)
c
      if (k.eq.0) then
       call getmem(nsq,ib)
       call squr(b,bl(ib),n)
       call mxma(bl(ia),1,n,bl(ib),1,n,bl(ic),1,n,n,n,n)
       call retmem(1)
      endif
      if (k.eq.1) then
cc       ir=ir-1
       call mxma(bl(ia),1,n,b,1,n,bl(ic),1,n,n,n,n)
      endif
c
      if (m.eq.0) call tri(bl(ic),c,n)
      if (m.eq.1) call tfer(bl(ic),c,nsq)
c
c     call retmem(ir)
c
      call retmem(2)
      return
 9998  call message('mamu',
     1  'incorrect use of mamu: integer arguments must be 1 or 0',
     2  k,l)
      return
      end
c
c======================================================================
      subroutine squr(t,sq,n)
      implicit real*8(a-h,o-z)
c     triangle to square
      dimension sq(1),t(1)
      ij=0
      ii=0
      do 20 i=1,n
      jj=0
c$dir no_recurrence
cdir$ ivdep
      do 10 j=1,i
      ij=ij+1
      sq(ii+j)=t(ij)
      sq(jj+i)=t(ij)
10    jj=jj+n
20    ii=ii+n
      return
      end
c==============================================================
      subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
c     this code should run at about 34 mflops for matrices
c     of dimensions between 50 and 500 on a 25 mhz rs6000
      implicit real*8 (a-h,o-z)
      dimension  a(*),b(*),r(*)
      parameter (nb=60)
c     written based on tips by  ron bell, ibm united kingdom, ltd.
c     seeibm document no. gg24-3611
c     note that this code, in spite of all efforts, does not reach
c     the numerical performance (43 mflops/s) claimed in the above
c     publication. it is assumed now that the performance figures
c     of bell refer to the inner loop only and do not reflect
c     the overhead.  the maximum performance of this code is about
c     35 mflops/s
      dimension rr(nb,nb),aa(nb,nb),bb(nb,nb)
c     performs r=a*b; the actual dimensions are naxma,maxmb,naxmb
c     the addressing function is: r(i,j)=r((i-1)*mcolr+(j-1)*mrowr+1)
c     a(i,k)=a((i-1)*mcola+(k-1)*mrowa+1)
c     b(k,j)=b((k-1)*mcolb+(j-1)*mrowb+1)
c     good example of obscure addressing in fortran
c     please do not write code like this. unfortunately, this piece of
c     code is widely used on the crays
c     the buffer size has been roughly optimized for the apollo dn10k
      ir=1
      do 3 j=1,nrow
      irr=ir
       do 2 i=1,ncol
c       ncol is the number of rows. the dutch mind is perplexing
         r(irr)=0.0d0
         irr=irr+mcolr
 2     continue
       ir=ir+mrowr
 3     continue
      do 1400 ii=1,ncol,nb
        ie=min0(ii+nb-1,ncol)
        ie1=ie-ii+1
        ie2=3*(ie1/3)
c       write(*,*) ie,ie1
          do 1300 jj=1,nrow,nb
       je=min0(jj+nb-1,nrow)
       je1=je-jj+1
            je2=3*(je1/3)
c      make sure that ie2 and je2 are divisible by 3 for the
c      loop unrolling
       ir=(ii-1)*mcolr+(jj-1)*mrowr+1
       do 200 i=1,ie1
                  irr=ir
              do 100 j=1,je1
c               rr(i,j)=r(irr)
                rr(i,j)=0.0d0
                irr=irr+mrowr
 100          continue
                  ir=ir+mcolr
 200        continue
        do 1000 kk=1,nlink,nb
         ke=min0(kk+nb-1,nlink)
         ke1=ke-kk+1
         ia=(ii-1)*mcola+(kk-1)*mrowa+1
         do 400 i=1,ie1
                iaa=ia
                do 300 k=1,ke1
                  aa(k,i)=a(iaa)
                  iaa=iaa+mrowa
 300            continue
                ia=ia+mcola
 400          continue
          ib=(kk-1)*mcolb+(jj-1)*mrowb+1
           do 600 j=1,je1
             ibb=ib
             do 500 k=1,ke1
               bb(k,j)=b(ibb)
               ibb=ibb+mcolb
 500              continue
                  ib=ib+mrowb
 600            continue
            do 900 i=1,ie2,3
               do 800 j=1,je2,3
               s00=rr(i,j)
               s10=rr(i+1,j)
                    s20=rr(i+2,j)
               s01=rr(i,j+1)
               s11=rr(i+1,j+1)
                    s21=rr(i+2,j+1)
                    s02=rr(i,j+2)
                    s12=rr(i+1,j+2)
                    s22=rr(i+2,j+2)
               do 700 k=1,ke1
        s00=s00+aa(k,i)*bb(k,j)
        s01=s01+aa(k,i)*bb(k,j+1)
                      s02=s02+aa(k,i)*bb(k,j+2)
        s10=s10+aa(k,i+1)*bb(k,j)
        s11=s11+aa(k,i+1)*bb(k,j+1)
                      s12=s12+aa(k,i+1)*bb(k,j+2)
                      s20=s20+aa(k,i+2)*bb(k,j)
                      s21=s21+aa(k,i+2)*bb(k,j+1)
                      s22=s22+aa(k,i+2)*bb(k,j+2)
 700                continue
                    rr(i,j)=s00
                      rr(i+1,j)=s10
                    rr(i+2,j)=s20
              rr(i,j+1)=s01
               rr(i+1,j+1)=s11
                    rr(i+2,j+1)=s21
                    rr(i,j+2)=s02
                    rr(i+1,j+2)=s12
                    rr(i+2,j+2)=s22
 800              continue
                  do 820 j=je2+1,je1
                    s00=rr(i,j)
                    s10=rr(i+1,j)
                    s20=rr(i+2,j)
                    do 810 k=1,ke1
        s00=s00+aa(k,i)*bb(k,j)
        s10=s10+aa(k,i+1)*bb(k,j)
                      s20=s20+aa(k,i+2)*bb(k,j)
 810             continue
                    rr(i,j)=s00
                      rr(i+1,j)=s10
                    rr(i+2,j)=s20
 820             continue
 900             continue
                 do 950 i=ie2+1,ie1
                  do 940 j=1,je1
                   do 930 k=1,ke1
                    rr(i,j)=rr(i,j)+aa(k,i)*bb(k,j)
 930               continue
 940              continue
 950             continue
 1000           continue
                    ir=(ii-1)*mcolr+(jj-1)*mrowr+1
           do 1200 i=1,ie1
             irr=ir
             do 1100 j=1,je1
               r(irr)=r(irr)+rr(i,j)
                    irr=irr+mrowr
 1100             continue
                 ir=ir+mcolr
 1200            continue
 1300           continue
 1400          continue
              end
c==============================================================
      subroutine nerror(noer,routine,message,n1,n2)

      use memory
      use messages

      real*8 tcpu,telap
      character*(*) routine
      character*(*) message
      character*256 errm1,errm2
      character*26 datestr,timestr
      parameter (sixty=60.0d0)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest
      write(errm1,2000)noer,routine(1:len_trim(routine))
2000  format('** Error no. ',i0,' in ',a,' **')
      If(n1.NE.0.OR.n2.NE.0) Then
       write(errm2,2001),message(1:len_trim(message)),n1,n2
      Else
       write(errm2,2002),message(1:len_trim(message))
      EndIf
2001  format(' ',a,' parameters:',i0,1x,i0)
2002  format(' ',a)
      call setchval('errm1',errm1)
      call setchval('errm2',errm2)
c
      call date1(datestr)
      call chartime(timestr)
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
      call secund(tcpu)
      call elapsec(telap)
      wbuf=''
      call addbuf(errm1)
      write(iout,'(a)')errm1
      call addbuf(errm2)
      write(iout,'(a)')errm2
      call display_buf(iout,'PQS Fatal Error:')
      write(iout,3500)
      write(iout,3400) nreq,nmark,lastadr-ioffset,mxmem,memtot
 3400 format(/,'Memory status:',/,
     *' request number=',i4,' memory marks=',i3,
     *' last used address=',i9,/,
     *' high water=',i9,' total available memory=',i9)
      write(iout,3700) tcpu/sixty,telap/sixty
      write(iout,4000) timestr
      write(iout,4500)
 3500 format(72('='))
 3700 format(/'Total CPU and elapsed time =',2f10.2,' min')
 4000 format(/'Termination on ',a24)
 4500 format(72('=')/)
      call para_stop
      write(iout,'(a)') '*** PQS Error Termination ***'
cc    call abort()
      STOP '*** PQS Error Termination ***'
      end
c==============================================================
      subroutine message(messag1,messag2,n1,n2)
      character*(*) messag1
      character*(*) messag2
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest
      iout = igetival('iout')
      write(iout,'(a)') trim(messag1)
      If(n1.NE.0.OR.n2.NE.0) Then
       write(iout,'(a,1x,i0,1x,i0)') trim(messag2),n1,n2
      Else
       write(iout,'(a)') trim(messag2)
      EndIf
      return
      end
c==============================================================
      subroutine osinv (a,n,d,tol,l,m)
      implicit real*8 (a-h,o-z)
c
c     parameters:  a - input matrix , destroyed in computation and repla
c                      by resultant inverse (must be a general matrix)
c                  n - order of matrix a
c                  d - resultant determinant
c            l and m - work vectors of length n
c                tol - if pivot element is less than this parameter the
c                      matrix is taken for singular (usually = 1.0e-8)
c     a determinant of zero indicates that the matrix is singular
c
      dimension a(1), m(1), l(1)
      d=1.0d0
      nk=-n
      do 180 k=1,n
         nk=nk+n
         l(k)=k
         m(k)=k
         kk=nk+k
         biga=a(kk)
         do 20 j=k,n
            iz=n*(j-1)
         do 20 i=k,n
            ij=iz+i
c
c     10 follows
c
            if (abs(biga)-abs(a(ij))) 10,20,20
   10       biga=a(ij)
            l(k)=i
            m(k)=j
   20    continue
         j=l(k)
         if (j-k) 50,50,30
   30    ki=k-n
         do 40 i=1,n
            ki=ki+n
            holo=-a(ki)
            ji=ki-k+j
            a(ki)=a(ji)
   40    a(ji)=holo
   50    i=m(k)
         if (i-k) 80,80,60
   60    jp=n*(i-1)
         do 70 j=1,n
            jk=nk+j
            ji=jp+j
            holo=-a(jk)
            a(jk)=a(ji)
   70    a(ji)=holo
   80    if (abs(biga)-tol) 90,100,100
   90    d=0.0d0
         return
  100    do 120 i=1,n
            if (i-k) 110,120,110
  110       ik=nk+i
            a(ik)=a(ik)/(-biga)
  120    continue
         do 150 i=1,n
            ik=nk+i
            ij=i-n
         do 150 j=1,n
            ij=ij+n
            if (i-k) 130,150,130
  130       if (j-k) 140,150,140
  140       kj=ij-i+k
            a(ij)=a(ik)*a(kj)+a(ij)
  150    continue
         kj=k-n
         do 170 j=1,n
            kj=kj+n
            if (j-k) 160,170,160
  160       a(kj)=a(kj)/biga
  170    continue
         d=d*biga
         a(kk)=1.0d0/biga
  180 continue
      k=n
  190 k=k-1
      if (k) 260,260,200
  200 i=l(k)
      if (i-k) 230,230,210
  210 jq=n*(k-1)
      jr=n*(i-1)
      do 220 j=1,n
         jk=jq+j
         holo=a(jk)
         ji=jr+j
         a(jk)=-a(ji)
  220 a(ji)=holo
  230 j=m(k)
      if (j-k) 190,190,240
  240 ki=k-n
      do 250 i=1,n
         ki=ki+n
         holo=a(ki)
         ji=ki+j-k
         a(ki)=-a(ji)
  250 a(ji)=holo
      go to 190
  260 return
      end
c==============================================================
      subroutine absmin(n,a,i,xmin)
c  this routine returns the lowest absolute value in the
c array a, from a(1) to a(n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      i=1
      xmin=abs(a(1))
      do 100 j=2,n
        if (xmin.gt.abs(a(j))) then
           xmin=abs(a(j))
           i=j
        end if
 100  continue
      end
c==============================================================
      subroutine absmax(n,a,i,xmax)
c  this routine returns the highest absolute value in the
c array a, from a(1) to a(n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
c blas      
      i=idamax(n,a,1)
      xmax=abs(a(i))
c     i=1
c     xmax=abs(a(1))
c     do 100 j=2,n
c       if (xmax.lt.abs(a(j))) then
c          xmax=abs(a(j))
c          i=j
c       end if
c100  continue
      end
c==============================================================
      subroutine absmaxdif(n,a,b,i,xmax)
c  this routine returns the highest absolute value in the
c array a, from a(1) to a(n)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      i=1
      xmax=abs(a(1)-b(1))
      do 100 j=2,n
        if (xmax.lt.abs(a(j)-b(j))) then
           xmax=abs(a(j)-b(j))
           i=j
        end if
 100  continue
      end
c==============================================================
      subroutine mxmb(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
c     this code should run at about 34 Mflops for matrices
c     of dimensions between 50 and 500 on a 25 MHz RS6000
      implicit real*8 (a-h,o-z)
      dimension  a(*),b(*),r(*)
      parameter (nb=60)
      dimension rr(nb,nb),aa(nb,nb),bb(nb,nb)
c     Written based on tips by  Ron Bell, IBM United Kingdom, Ltd.
c     SeeIBM Document No. GG24-3611
c     Note that this code, in spite of all efforts, does not reach
c     the numerical performance (43 Mflops/s) claimed in the above
c     publication. It is assumed now that the performance figures
c     of Bell refer to the inner loop only and do not reflect
c     the overhead.  The maximum performance of this code is about
c     35 Mflops/s
c     performs r=a*b; the actual dimensions are naxma,maxmb,naxmb
c     the addressing function is: r(i,j)=r((i-1)*mcolr+(j-1)*mrowr+1)
c     a(i,k)=a((i-1)*mcola+(k-1)*mrowa+1)
c     b(k,j)=b((k-1)*mcolb+(j-1)*mrowb+1)
c     It calculates  the actual matrix r(ncol,nrow) from a(ncol,nlink)
c      and b(nlink,nrow)
c    For normal matrices, i.e. the usual case, mcola, mcolb and mcolr
c    are 1, and mrowa=ma, mrowb=mb and mrowr=mr if the
c    DECLARATION of these matrices is a(ma,*), b(mb,*), r(mr,*)
c    For transposed matrices, interchange ncola and nrowa. Now
c    ncola=ma, nrowa=1, etc.
c     Good example of obscure addressing in FORTRAN
c     Please do not write code like this. Unfortunately, this piece of
c     code is widely used on the Crays
c     the buffer size has been roughly optimized for the apollo dn10k
      ir=1
      do 1400 ii=1,ncol,nb
        ie=min0(ii+nb-1,ncol)
        ie1=ie-ii+1
        ie2=3*(ie1/3)
c       write(*,*) ie,ie1
          do 1300 jj=1,nrow,nb
       je=min0(jj+nb-1,nrow)
       je1=je-jj+1
            je2=3*(je1/3)
c      make sure that ie2 and je2 are divisible by 3 for the
c      loop unrolling
       ir=(ii-1)*mcolr+(jj-1)*mrowr+1
       do 200 i=1,ie1
                  irr=ir
              do 100 j=1,je1
c               rr(i,j)=r(irr)
                rr(i,j)=0.0d0
                irr=irr+mrowr
 100          continue
                  ir=ir+mcolr
 200        continue
        do 1000 kk=1,nlink,nb
         ke=min0(kk+nb-1,nlink)
         ke1=ke-kk+1
         ia=(ii-1)*mcola+(kk-1)*mrowa+1
         do 400 i=1,ie1
                iaa=ia
                do 300 k=1,ke1
                  aa(k,i)=a(iaa)
                  iaa=iaa+mrowa
 300            continue
                ia=ia+mcola
 400          continue
          ib=(kk-1)*mcolb+(jj-1)*mrowb+1
           do 600 j=1,je1
             ibb=ib
             do 500 k=1,ke1
               bb(k,j)=b(ibb)
               ibb=ibb+mcolb
 500              continue
                  ib=ib+mrowb
 600            continue
            do 900 i=1,ie2,3
               do 800 j=1,je2,3
               s00=rr(i,j)
               s10=rr(i+1,j)
                    s20=rr(i+2,j)
               s01=rr(i,j+1)
               s11=rr(i+1,j+1)
                    s21=rr(i+2,j+1)
                    s02=rr(i,j+2)
                    s12=rr(i+1,j+2)
                    s22=rr(i+2,j+2)
               do 700 k=1,ke1
        s00=s00+aa(k,i)*bb(k,j)
        s01=s01+aa(k,i)*bb(k,j+1)
                      s02=s02+aa(k,i)*bb(k,j+2)
        s10=s10+aa(k,i+1)*bb(k,j)
        s11=s11+aa(k,i+1)*bb(k,j+1)
                      s12=s12+aa(k,i+1)*bb(k,j+2)
                      s20=s20+aa(k,i+2)*bb(k,j)
                      s21=s21+aa(k,i+2)*bb(k,j+1)
                      s22=s22+aa(k,i+2)*bb(k,j+2)
 700                continue
                    rr(i,j)=s00
                      rr(i+1,j)=s10
                    rr(i+2,j)=s20
              rr(i,j+1)=s01
               rr(i+1,j+1)=s11
                    rr(i+2,j+1)=s21
                    rr(i,j+2)=s02
                    rr(i+1,j+2)=s12
                    rr(i+2,j+2)=s22
 800              continue
                  do 820 j=je2+1,je1
                    s00=rr(i,j)
                    s10=rr(i+1,j)
                    s20=rr(i+2,j)
                    do 810 k=1,ke1
        s00=s00+aa(k,i)*bb(k,j)
        s10=s10+aa(k,i+1)*bb(k,j)
                      s20=s20+aa(k,i+2)*bb(k,j)
 810             continue
                    rr(i,j)=s00
                      rr(i+1,j)=s10
                    rr(i+2,j)=s20
 820             continue
 900             continue
                 do 950 i=ie2+1,ie1
                  do 940 j=1,je1
                   do 930 k=1,ke1
                    rr(i,j)=rr(i,j)+aa(k,i)*bb(k,j)
 930               continue
 940              continue
 950             continue
 1000           continue
                    ir=(ii-1)*mcolr+(jj-1)*mrowr+1
           do 1200 i=1,ie1
             irr=ir
             do 1100 j=1,je1
               r(irr)=r(irr)+rr(i,j)
                    irr=irr+mrowr
 1100             continue
                 ir=ir+mcolr
 1200            continue
 1300           continue
 1400          continue
              end
c=======================================================
      subroutine trspmo(amat,m,bmat,n)
c This is matrix transpose A(T)--> b  (blas replacement)
      implicit real*8 (a-h,o-z)
      dimension amat(m,n), bmat(n,m)
c
      do 100 i=1,m
      do 100 j=1,n
      bmat(j,i)=amat(i,j)
  100 continue
      end
c========================================================
c
      subroutine stringcat(fname,l1,ext,l2,fname1,j)
c     this routine removes leading and trailing blanks from
c     fname and concatenates it with extn.
c     the result is in fname1.
c     Parameters:
c     Input:
c     fname: irst part of the file name
c     l1: the maximum length of both fname and the result, fname1
c     ext: the extension (characters). Note that the dot (.) is not
c     automatically assumed
c     l2: the maximum length of ext
c     Output:
c     fname1: the result
c     j: the real length of fname1, omitting trailing blanks
c
      implicit none
      character*(*) fname, fname1
      character*(*) ext
      integer l1,l2,j,lenr,lenf
      fname1=''
      fname1=adjustl(fname)
      lenf=len_trim(fname1)
      lenr=l1-lenf
      if(lenr.ge.l2)then
        fname1(lenf+1:lenf+l2)=ext(1:l2)
      else if(lenr.gt.0)then
         fname1(lenf+1:lenf+lenr)=ext(1:lenr)
      endif
      j=len_trim(fname1)
      end
c==============================================================
      subroutine druma (d,ncf,ifi,txt)
      implicit real*8 (a-h,o-z)
      character*8 txt
      dimension d(1)
c
      write (ifi,20) txt
c----------------------------------------------------------------------c
c do not print anything if it is too long :
      if(ncf.gt.200) then
          write(ifi,*)
     *  ' it is not printed out for NCF=',ncf,' > 200 '
         return
      endif
c----------------------------------------------------------------------c
      n=ncf
      ja=1
      do 10 i=1,n
         je=ja+i-1
         write (ifi,30) i,(d(j),j=ja,je)
         ja=je+1
   10 continue
      return
c
   20 format (30x,3h***,a8,3h***)
   30 format (1x,i4,2x,10f12.7,/,(7x,10f12.7))
cc 30 format (1x,i4,2x,10f10.4,/,(7x,10f10.4))
c
      end
c=======================================================================
      REAL*8 FUNCTION ARCOS (X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (X.GE.1.0D0) GO TO 10
      IF (X.LE.-1.0D0) GO TO 20
      X1=SQRT(1.0D0-X**2)
      IF (ABS(X).LT.1.0D-11) GO TO 30
      S=ATAN(X1/X)
      IF (X.LT.0.0) S=S+3.1415926536D0
      ARCOS=S
      RETURN
   10 ARCOS=0.0D0
      RETURN
   20 ARCOS=3.1415926536D0
      RETURN
   30 ARCOS=1.5707963268D0
      RETURN
C
      END
c====================================================================
      SUBROUTINE VEKTOR (U,R,I,J,XNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3), Xnuc(3,*)
C
C       BILDET DEN NORMIERTEN ENTFERNUNGSVEKTOR VOM KERN J NACH KERN I
C        UND DIE ENTFERNUNG R
C
      r=0.0d0
      do l=1,3
      u(l)=xnuc(l,i)-xnuc(l,j)
      r=r+u(l)**2
      end do
      call nom1(u)
      end
c=========================================================================
      SUBROUTINE NOM1(U)
c This was originally NOM but JB has a different nom
c it normalizes a vector
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3)
      X=1.0d0/SQRT(SCALAR(U,U))
      DO 10 I=1,3
      U(I)=U(I)*X
   10 CONTINUE
      RETURN
      END
C========================================================================
c      REAL*8 FUNCTION S2 (X)
C See in service1
c      IMPLICIT REAL*8 (A-H,O-Z)
c      S2=SQRT(1.0D0-X**2)
c      END
C========================================================================
      REAL*8 FUNCTION SCALAR (U,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3), V(3)
      SCALAR=0.0D0
      DO 10 I=1,3
      SCALAR=SCALAR+U(I)*V(I)
   10 CONTINUE
      END
C
c      SUBROUTINE NORMAL (U,V,W)
c See in service1
c      IMPLICIT REAL*8 (A-H,O-Z)
c      DIMENSION U(3), V(3), W(3)
Cc
Cc     99999...  W WIRD EIN SENKRECHT AUF DIE EBENE(U,V) STEHENDER EINHEI
Cc      TOR
Cc
c      W(1)=U(2)*V(3)-U(3)*V(2)
c      W(2)=U(3)*V(1)-U(1)*V(3)
c      W(3)=U(1)*V(2)-U(2)*V(1)
c      CALL NOM (W)
c      RETURN
C
c      END
cGG===================================================================
c      subroutine DetectNaN(a,n,NaNFound)
cc  This is for the PC - the Digital compiler can check for NaN
c      real*8 a
c      logical IsNaN
c      dimension a(n)
c      NaNFound=0
c      do i=1,n
c        if(IsNaN(a(i))) then
c          NaNFound=1
c        end if
c      end do
c      end
c========================================================================
      subroutine get_ij_half(ij, i, j)
c
c extracts indices i & j from a common index ij=i*(i-1)/2 + j
c
      implicit none
      integer ij, i, j
      intrinsic sqrt, float
c
      i = sqrt(float(ij + ij))
      j = ij - i*(i-1)/2
      if (i .lt. j) then
         i = i + 1
         j = ij - i*(i-1)/2
      endif
c
      end
c-----------------------------------------------------------------
      subroutine get_ij_full(ij,nj, i, j)
c
c     extracts indices i and j a from common index ij=(i-1)*nj + j
c
      implicit none
      integer ij,nj, i, j
c
      i=ij/nj+1
      j=ij-(i-1)*nj
c
      if(j.eq.0) then
         i=i-1
         j=ij-(i-1)*nj
      endif
c
      end
c-----------------------------------------------------------------
      subroutine symfunc(reuse)

      use memory 

      implicit real*8 (a-h,o-z)
c  this routine reserves real & integer memory for the symmetry determining
c  routine sympar. The latter determines symmetry pairs of contracted shells
c  and individual functions.
c  INPUT:  reuse - a logical flag indicating whether or not the
c          basis/shell symmetry pair data can OVERWRITE existing data
c
c  Uses common /intbl/ (the INX array) and common/symm/
c     common /big/bl(100)
c     common /intbl/maxsh,inx(100)
      common /symm/nsym,nsy(7)
      logical reuse
c
c  ** reserve space for the symmetry equivalence array
c   get job parameters
      call getival('ncs',ncs)
      call getival('ncf',ncf)
      call getival('nsh',nsh)
      call getival('inuc',inuc)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      call getival('iout',iout)
c ......................................
c -- be careful with memory pointers here
      call tstival('SymFunPr',iexist)
      If(.not.reuse.or.iexist.eq.0) Then
        call getint(7*ncf,ifp)
        call getint(7*ncs,ifp1)
        call setival('SymFunPr',ifp)
        call setival('SymShPr',ifp1)
      Else
        call getival('SymFunPr',ifp)
        call getival('SymShPr',ifp1)
c -- need to put current coordinates into basdat array
        call getival('na  ',na)
        call GetCoord(na,nsh,ncs,bl(ictr),bl(inuc),bl(ibas))
      EndIf
c ......................................
c  bl(ifp) is the starting address of an integer array, ifup(7,ncf)
c  in sympar, which holds the symmetry pairs of the contracted functions
      call izeroit(bl(ifp1),7*ncs)
c  bl(ifp1) is the starting address of an integer array, ishp(7,ncs)
c   in sympar, which holds the symmetry pairs of the contracted shells
c
      call getrval('acc',acc)
      eps=acc*1000.0d0
c     establish symmetry pairs
      call getmem(13*nsh,ibas2)
      call sympar(ncs,nsh,ncf,nsym,nsy,bl(inuc),bl(ibas),bl(ibas2),
     1            bl(ictr),bl(ifp),bl(ifp1),eps,iout)
      call retmem(1)
      end
c-----------------------------------------------------------------
      subroutine symmoff
      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /memors/ nsymx,ijshp,isymm
      common /symmet/ rnsym
      common /symm/nsymm,msy(7)
c
c -- save current value of nsym
      call setival('nsymC',nsym)
c
cc      write(6,*)' symmetry off : nsym=',nsym
c
      nsym=0
      nsymx=0
      rnsym=1.0d0
      nsymm=0
      call setival('nsym',0)
c
      end
c-----------------------------------------------------------------
      subroutine symmon
      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /memors/ nsymx,ijshp,isymm
      common /symmet/ rnsym
      common /symm/nsymm,msy(7)
      data one /1.d0/
c
c -- restore original value of nsym
      call getival('nsymC',nsymC)
c
      nsym = nsymC
      nsymx=nsymC
      rnsym=one/(one+dble(nsym))
      nsymm=nsymC
      call setival('nsym',nsym)
c
cc      write(6,*)' symmetry on : nsym=',nsym
c
      end
