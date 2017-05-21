c-----------------------------------------------------------------------
      subroutine opt_mol(nq,ngeo,x,g,e,fc,cq,enext,dismax,icmx,hq)
c-----------------------------------------------------------------------
c     Optimization in natural internal (BMAT) coordinates.
c     This routine interfaces the Molpro Optimization methods to BMAT.
c     F.Eckert Jan 96 - Jun 96, latest additions Jan 99
c
      implicit double precision (a-h,o-z)
c      character*16 opttyp
      include "common/cgeom"
      include "common/tapes"
      include "common/coptg"
      include "common/czmat"
      include "common/big"
      common/cgadr_opt/ ix,ig,ie,ih
      dimension x(nq,ngeo),g(nq,ngeo),e(ngeo),fc(nq),cq(nq),hq(nq,nq)
c
      write(iout,'(1X,''BMAT using opt '',A/)') opttyp(1:lenstr(opttyp))
c
      nsave=n
      nptsave=npt
      n=nq
      npt=ngeo
c     build hessian
      ibase=icorr(0)
      if (.not.modelhes .and. .not.hstart .and. zinitcart.le.0) then
        call fzero(hq,n*n)
        do i=1,n
         if (dabs(fc(i)).gt.1.d-6) hq(i,i)=1.0d0/fc(i)
        end do
      end if
      if (iprhes.gt.0) call outsqr(hq,n,n,n,'In BMAT: Hessian')
c     scale everything to [a.u]
      call skal0(x,g,hq,n,npt)
c      call dscal_X(n*npt,-1.0d0,g,1)
      if (iprhes.gt.0)
     > call outsqr(hq,n,n,n,'In BMAT: Hessian after scaling')
c...  updating of the hessian
      if (ibfgs.gt.0.and.maxupd.ne.0) then
       call opt_bfgs(x,g,e,hq)
      else
       call gethes(x,g,e,hq)
      end if
c...  update for transition states
      if (iroot.ne.1.and.ihflag.ne.1) then
       call tsupdate(q(ix),q(ig),q(ih))
      end if
c     call optimizer
      ixnew=icorr(n)
      ignew=icorr(n)
      idone=0
      mexclu=0
      if (opttyp(1:2).eq.'AH') call opt_ah(x,g,e,hq,q(ixnew),
     >    q(ignew),enext,conv,convl)
      if (opttyp(1:2).eq.'AK') call opt_ak(x,g,e,hq,q(ixnew),
     >    q(ignew),enext,conv,convl)
      if (opttyp(1:2).eq.'SR') then
        write (6,'(/'' SR not available, making RF step instead'')')
        call opt_rf(x,g,e,hq,q(ixnew),q(ignew),enext)
      end if
      if (opttyp(1:3).eq.'QSD') call opt_qsd(x,g,e,
     >    hq,q(ixnew),q(ignew),enext,convl)
      if (opttyp(1:2).eq.'DI') call opt_ds(x,g,e,hq,q(ixnew),
     >    q(ignew),enext,mexclu,conv,convl)
      if (opttyp(1:2).eq.'RF') call opt_rf(x,g,e,hq,q(ixnew),
     >    q(ignew),enext)
c     scale things back
      call skal1(x,g,hq,q(ixnew),n,npt)
c      call dscal_X(n*npt,-1.0d0,g,1)
c     bmat needs dx not xnew
      call subvec(cq,q(ixnew),x(1,npt),n)
c     get largest displacement
      call absmax(n,cq,icmx,dismax)
      n=nsave
      npt=nptsave
c
      call corlsr(ibase)
      return
      end
c-----------------------------------------------------------------------
      subroutine skal0(x,g,h,n,npt)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(n,npt),g(n,npt),h(n,n)
c
c     scale coordinates, gradients and hessian for Pulay's optimizer
c
      ang=get_scalar('ANGSTROM')
      ajoule=0.22936757d0
c
      call dscal_X(n*npt,ang,x,1)
      call dscal_X(n*npt,-1.0d0*ajoule/ang,g,1)
      call dscal_X(n*n,ajoule/(ang*ang),h,1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine skal1(x,g,h,xnew,n,npt)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(n,npt),g(n,npt),h(n,n),xnew(n)
c
c     scale back coordinates, gradients and hessian
c
      ang=get_scalar('ANGSTROM')
      ajoule=0.22936757d0
c
      call dscal_X(n*npt,1.0d0/ang,x,1)
      call dscal_X(n*npt,-1.0d0*ang/ajoule,g,1)
      call dscal_X(n*n,ang*ang/ajoule,h,1)
      call dscal_X(n,1.0d0/ang,xnew,1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hes_bmat(hx,hq,BB,B,ibcontr,nq,nek)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c This routine transforma a Cartesian Hessian hx to (possibly  redundant)
c internal coordinates
c BB is storage for the full (sparse) B matrix
c B  contains the non-zero elements of B
c ibcontr gives the encoding (the atoms) to expand B to BB
c nq= number of interbnal coordinates
c nek=3*natoms
      dimension hx(nek,nek),hq(nq,nq),bb(nq,nek),
     >          b(54,nq),ibcontr(20,nq)
      include "common/big"
c
c     Transform cartesian hessian to internal BMAT hessian
c     b and ibcontr hold encoded B-matrix
c
      ang=get_scalar('ANGSTROM')
      ajoule=0.22936757d0
      ibase=icorr(0)
c
c     explicitly construct B-matrix
c
      call fzero(bb,nq*nek)
      do i=1,nq
       natom=ibcontr(2,i)
       k3=0
       do k=1,natom
        k3=k3+3
        iatom=ibcontr(k+2,i)
        iat3=iatom*3
        bb(i,iat3-2)=b(k3-2,i)
        bb(i,iat3-1)=b(k3-1,i)
        bb(i,iat3)=b(k3,i)
       end do
      end do
c      call outsqr(BB,nq,nq,nek,'B-matrix')
c
c     scale model hessian to Angstrom
c
c      call outsqr(Hx,nek,nek,nek,'Cart. Hessian before scaling')
      call dscal_X(nek*nek,ang*ang/ajoule,hx,1)
c      call outsqr(Hx,nek,nek,nek,'Cart. Hessian after scaling')
c
c     calculate the generalized inverse of B: A=Bt(BBt)^-1
c
      ig=icorr(nq*nq)
      call mxma(bb,1,nq,bb,nq,1,q(ig),1,nq,nq,nek,nq)
c      call outsqr(q(iG),nq,nq,nq,'G')
      indx=icorr(nq*nq)
      iwork1=icorr(nq*nq)
      iwork2=icorr(nq*nq)
      iv1=icorr(nek)
      irtps=icorr(nq)
      ivec=icorr(nq*nq)
      call fmove(q(ig),q(ivec),nq*nq)
      call diag2(nq,nq,q(irtps),q(ivec))
      do i=1,nq
        if(dabs(q(irtps+i-1)).gt.1.d-6) then
          q(irtps+i-1)=1.0d0/q(irtps+i-1)
        else
          q(irtps+i-1)=0.0d0
        end if
      end do
      iginv=icorr(nq*nq)
      call fzero(q(iginv),nq*nq)
      do i=1,nq
       q(iginv+nq*(i-1)+i-1)=q(irtps+i-1)
      end do
      call mxma(q(ivec),1,nq,q(iginv),1,nq,q(iwork1),1,nq,nq,nq,nq)
      call mxma(q(iwork1),1,nq,q(ivec),nq,1,q(ig),1,nq,nq,nq,nq)
c      call matinv(q(ig),nq,nq,q(indx),q(iwork1),q(iwork2),q(iv1))
c
      ia=icorr(nq*nek)
      call mxma(bb,nq,1,q(ig),1,nq,q(ia),1,nek,nek,nq,nq)
c      call outsqr(q(iA),nek,nek,nq,'Inverse(B) = A')
c     test
c      icheck=icorr(nq*nq)
c      call mxma(BB,1,nq,q(iA),1,nek,q(icheck),1,nq,nq,nek,nq)
c      call outsqr(q(icheck),nq,nq,nq,' B*A should be 1')
c      icheck1=icorr(nek*nek)
c      call mxma(q(iA),1,nek,BB,1,nq,q(icheck1),1,nek,nek,nq,nek)
c      call outsqr(q(icheck1),nek,nek,nek,' A*B should be ?')
c
c     Hq=AtHxA
c
      ie1=icorr(nq*nek)
      call mxma(q(ia),nek,1,hx,1,nek,q(ie1),1,nq,nq,nek,nek)
      call mxma(q(ie1),1,nq,q(ia),1,nek,hq,1,nq,nq,nek,nq)
c      call outsqr(Hq,nq,nq,nq,'Hessian in internal coordinates')
c
c  Scale the Hessian back
      call dscal_X(nek*nek,ajoule/(ang*ang),hx,1)
      call corlsr(ibase)
      return
      end
c-----------------------------------------------------------------------
      subroutine relaxq(nq,fi,hq,cq,edecr,dismax,icmx,scfac)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
C     THIS ROUTINE CALCULATES THE NEW INTERNAL COORDINATES
C     AND THE PROJECTED ENERGY LOWERINGS
C     NQ (INPUT) NUMBER OF INTERNAL COORDINATES
C     FI (INPUT) INTERNAL FORCES
C     hc (INPUT) 
C     cq  (OUTPUT) CHANGE IN THE INTERNAL COORDINATES
c     edecr INTENT(OUT) energy lowering
C     DISMAX (OUTPUT) MAXIMUM CHANGE IN THE INTERNAL COORDINATES
C     ICMX (OUTPUT) COORDINATE OF THE LARGEST DISPLACEMENT
C     SCFAC  SCALE FACTOR FOR THE DISPLACEMENTS, NORMALLY 1

      dimension fi(nq),hq(nq,nq),cq(nq)
      include "common/cgeom"
      include "common/coptg"
      include "common/big"
      ibase=icorr(0)
      ihinv=icorr(nq*nq)
      indx=icorr(nq*nq)
      iwork1=icorr(nq*nq)
      iwork2=icorr(nq*nq)
      iv1=icorr(nq)
c
      call fmove(hq,q(ihinv),nq*nq)
      if(iprhes.ge.0)
     >call outsqr(q(iHinv),nq,nq,nq,'Hessian [ang^2/ajoule]')
c      call matinv(q(ihinv),nq,nq,q(indx),q(iwork1),q(iwork2),q(iv1))
      call hesinv(nq,q(iHinv),q(indx),q(iv1),q(iwork1),q(iwork2))
      if(iprhes.ge.0)
     >call outsqr(q(iHinv),nq,nq,nq,'Inverse Hessian')
      call mxva(q(ihinv),1,nq,fi,1,cq,1,nq,nq)
      if (scfac.ne.1.0d0) call dscal_X(nq,scfac,cq,1)
      edecr=-0.5d0*ddot_X(nq,cq,1,fi,1)
      call absmax(nq,cq,icmx,dismax)
c
      call corlsr(ibase)
      return
      end
c-----------------------------------------------------------------------
      subroutine relaxqred(nq,   fi,   hq,    proj,tol,
     &                     cq,   edecr,dismax,icmx,scfac)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
C     THIS ROUTINE CALCULATES THE NEW INTERNAL COORDINATES
C     AND THE PROJECTED ENERGY LOWERINGS
C     NQ (INPUT) NUMBER OF INTERNAL COORDINATES
C     FI (INPUT) INTERNAL FORCES
C     hq (INPUT) The Hessian (force constant) matrix in internal coordinates
c     proj  Projector to the non-redundant subspace of the int. coordinates
c     tol: the tolerance for zero eigenvalue, say 10**-6
C     cq  (OUTPUT) CHANGE IN THE INTERNAL COORDINATES
c     edecr INTENT(OUT) energy lowering
C     DISMAX (OUTPUT) MAXIMUM CHANGE IN THE INTERNAL COORDINATES
C     ICMX (OUTPUT) COORDINATE OF THE LARGEST DISPLACEMENT
C     SCFAC  SCALE FACTOR FOR THE DISPLACEMENTS, NORMALLY 1

      dimension fi(nq),hq(nq,nq),proj(nq,nq),cq(nq)
      include "common/cgeom"
      include "common/coptg"
      include "common/big"
      ibase=icorr(0)
      ihinv=icorr(nq*nq)
      indx=icorr(nq*nq)
      iwork1=icorr(nq*nq)
      iwork2=icorr(nq*nq)
      iv1=icorr(nq)
      iv2=icorr(nq)
c
      one=1.0d0
      zero=0.0d0
cpp
c  Try just to use a diagonal F.C. matrix on fi and then project the distortion
       do i=1,nq
         q(iv1+i-1)=fi(i)/hq(i,i)
       end do 
       call dgemv('n',nq,nq,one,proj,nq,q(iv1),1,zero,cq,1)
       goto 1000
cpp
      call fmove(hq,q(ihinv),nq*nq)
c  post- and premultiply H by P
      call dgemm_x('N','N',nq,nq,nq,
     &             one,q(ihinv),nq,proj,nq,
     &             zero,q(iwork1),nq)
      call dgemm_x('N','N',nq,nq,nq,
     &             one,proj,nq,q(iwork1),nq,
     &             zero,q(ihinv),nq)
c  Form the generalized inverse of H
      call geninv(nq,q(ihinv),det,tol,q(iwork1),  
     &            q(iv1),q(iv2),ninv,q(memproj))
c  post- and premultiply H- by P
      call dgemm_x('N','N',nq,nq,nq,
     &             one,q(ihinv),nq,proj,nq,
     &             zero,q(iwork1),nq)
      call dgemm_x('N','N',nq,nq,nq,
     &             one,proj,nq,q(iwork1),nq,
     &             zero,q(ihinv),nq)
c
cpp
      write(*,*) 'Projected force constant matrix'
      call prntmat(nq,nq,nq,q(ihinv))
      if(iprhes.ge.0)
     >call outsqr(q(iHinv),nq,nq,nq,'Hessian [ang^2/ajoule]')
      if(iprhes.ge.0)
     >call outsqr(q(iHinv),nq,nq,nq,'Inverse Hessian')
c      call mxva(q(ihinv),1,nq,fi,1,cq,1,nq,nq)
      call dgemv_x('n',nq,nq,one,q(ihinv),nq,fi,1,zero,cq,1)
cpp
 1000 continue
      if (scfac.ne.1.0d0) call dscal_X(nq,scfac,cq,1)
      edecr=-0.5d0*ddot_X(nq,cq,1,fi,1)
      call absmax(nq,cq,icmx,dismax)
c
      call corlsr(ibase)
      return
      end
c-----------------------------------------------------------------------
      subroutine hesinv(nq,h,vec,root,hinv,work)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h(nq,nq),vec(nq,nq),root(nq),hinv(nq,nq),work(nq,nq)
c
c     calculate the general inverse of a nq*nq matrix
c
      call fmove(h,vec,nq*nq)
      call diag2(nq,nq,root,vec)
      do i=1,nq
        if(dabs(root(i)).gt.1.d-6)
     >  root(i)=1.0d0/root(i)
      end do
      call fzero(hinv,nq*nq)
      do i=1,nq
       hinv(i,i)=root(i)
      end do
      call mxma(vec,1,nq,hinv,1,nq,work,1,nq,nq,nq,nq)
      call mxma(work,1,nq,vec,nq,1,h,1,nq,nq,nq,nq)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine osinv (a,n,d,tol,l,m)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
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
            if (dabs(biga)-dabs(a(ij))) 10,20,20
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
   80    if (dabs(biga)-tol) 90,100,100
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
c--------------------------------------------------------------
      subroutine add(a,b,n)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      do 100 i=1,n
        b(i)=b(i)+a(i)
 100  continue
      end
      subroutine add1(a,con,b,n)
c  Adds con times a to b, result in b
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      call daxpy_X(n,con,a,1,b,1)
      end
c--------------------------------------------------------------
      subroutine tfer (a,b,n)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      call dcopy_X(n,a,1,b,1)
      end
c--------------------------------------------------------------
      subroutine zeroit (a,n)
      implicit real*8 (a-h,o-z)
      dimension a(*)
c     call vec_$dzero (a,n)
      do  100 i=1,n
        a(i)=0.0d0
 100  continue
      end
c--------------------------------------------------------------
      subroutine mult(tmat,con,num)
      implicit real*8 (a-h,o-z)
      dimension tmat(1)
c
c     this subroutine multiplies an array by a constant
c     this is  jmb's invention; it is the same as dscal
      call  dscal_X(num,con,tmat,1)
c
      end
c--------------------------------------------------------------
      subroutine inner(x,y,sum,ne)
c   WHO ON EARTH WROTE THIS - THIS DUPLICATES SCALAR
      implicit real*8 (a-h,o-z)
      dimension x(1),y(1)
c     multiplies the elements in two arrays with the same index and adds
c     them similar to inner product
      sum=0.0d0
      m=ne
      do 30 l=1,m
      sum=sum+x(l)*y(l)
   30 continue
      return
      end
c--------------------------------------------------------------
c
      subroutine absmax(n,x,imax,xmax)
c
c     finds the maximum element xmax with index imax in
c     vector x(n) - maybe there is a similar blas routine
c
      implicit double precision (a-h,o-z)
      dimension x(n)
c
      imax = 1
      xmax = dabs(x(1))
      do 30 i = 2,n
         if(xmax.lt.dabs(x(i))) then
           xmax = dabs(x(i))
           imax = i
         end if
  30  continue
      return
      end
