      subroutine intforc (nq,   nek, b,    ibcontr, f,
     1                    fi,   qq,  iout, icond,   gmat,
     2                    ginv, mxf, fmax1)
c    nq (input) number of coordinates
c    nek (input) 3 times the number of atoms (# of cartesiians)
c    b(54,nq) (input) b matrix stored compressed
c    ibcontr(20,nq) coding of b
c    lq,mq dimensions of b (in the main program)
c    f(1..nek) (input) cartesian forces
c    fi(1..nq) (output) forces
c    qq (input) internal coordinates
c    iout,icond (input)  output files
c    gmat (output): contains the G matrix on output
c    ginv (output): contains the projector GG(-) where G(-)=the generalized inverse of G
c    mxf: coordinate of the maximum internal force
c    fmax1: maximum internal force
      implicit real*8 (a-h,o-z)
      include "common/big"
      external bbt
      dimension b(54,nq),ibcontr(20,nq),f(nek),fi(nq),qq(nq)
      dimension gmat(*),ginv(*)
c
      ibase=icorr(0)
      ic=icorr(nq)
c
cpp
c      write(*,*) 'inp,iout,icond=',inp,iout,icond
c      write(*,*) 'Entered intforc, GMAT='
c      call prntmat(nq,nq,nq,gmat)
c      call outsqr(gmat,nq,nq,nq,'G matrix in intforc')
c      call outsqr(ginv,nq,nq,nq,'G inverse in intforc')
      call bxv(nq,nek,b,ibcontr,f,q(ic))
c      call bdiag(nq,nek,b,ibcontr,q(idm))
c      call congrad(nq,nek,b,ibcontr,q(ix),q(iy),q(idm),q(ic),
c     1    q(it),q(ip),q(ir),q(iz),bbt)
c      call tfer(q(ix),fi,nq)
c      call twocol (iout,nq,qq,fi,'Internal coordinates and forces')
c
c  Multiply Bf (in q(ic)) by the generalized inverse of G
      call ginvxbf(nq,ginv,q(ic),fi)
      call twocol (icond,nq,qq,fi,'Internal coordinates and forces')
cpp
      call twocol (6,nq,qq,fi,'Internal coordinates and forces')
      call absmax(nq,fi,mxf,fmax1)
      fmax1=sign(fmax1,fi(mxf))
c      write(iout,400) mxf,fmax1
      if(icond.ne.0) write(icond,400) mxf,fmax1
 400  format(/' Maximum force on coord. ',i4, ' = ',f12.7)
 310  continue
c
c      call retmem(8)
c
      call corlsr(ibase)
      return
      end
c
      subroutine dispri (idis,iout,icond,inp2,ipu,na,ia,symb,xy)
c     this routine prints and punches the distorted cartesians
      implicit real*8 (a-h,o-z)
      dimension ia(na),xy(3,na)
      character*2 symb(na)
      include "common/cmpp"
c     dimension symb(na)
c      write (iout,1290)
      if (icond.ne.0) write (icond,1290)
 1290 format (/' New cartesian coordinates')
      if (ipu.ne.0.and.iprocs.eq.0) then
      write(ipu,1300) idis
 1300 format(' Displacement no.=',i4)
      end if
      i3=0
      if(iprocs.eq.0) rewind inp2
      do 740 i=1,na
         xiat=DBLE(ia(i))
c         write (iout,1310) i,symb(i),xiat,xy(1,i),xy(2,i),xy(3,i)
         if (icond.ne.0)
     1   write (icond,1310) i,symb(i),xiat,xy(1,i),xy(2,i),xy(3,i)
 1310 format (1x,i3,1x,'n=',a4,4x,f6.2,4x,3f10.6,f10.6)
         if(iprocs.eq.0) write (inp2,1320) symb(i),xiat,xy(1,i),
     1                   xy(2,i),xy(3,i)
 1320 format ('n=',a4,4x,f6.2,4x,3f10.6,f10.6)
         if (ipu.eq.0.or.iprocs.ne.0) go to 740
         write (ipu,1320) symb(i),xiat,xy(1,i),xy(2,i),xy(3,i)
  740 continue
      if(iprocs.eq.0) rewind inp2
c
      end
c
      subroutine readdiag(f,nq,nfix,ifix,inp,iout,icond)
c     this routine reads a diagonal force constant matrix, inverts them, 
c   and substitutes zero for the inverse if there are fixed coordinates
      implicit real*8 (a-h,o-z)
      real*8 l,m
      character*128 line
c     parameters:
c     f(nq) - output (force constants, inverted in program)
c     nq - input, number of int. coordinates'
c     nfix  input - number of coordinates to be fixed
c     ifix  input - coordinates to be fixed
c      inp  input - main input file
c     iout and icond input - output files
      dimension f(nq),ifix(*)
      include "common/cmpp"
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
cpp
      write(*,1180)
      if(icond.ne.0) write (icond,1180)
 1180 format (/1x,'Force constants'/)
      do il=1,nq,8
      j1=il
      j2=min(nq,j1+7)
      if(iprocs.eq.0) read(inp,'(a)') line
      if (nprocs.gt.1) call global_brdcst_str(line)
      read (line,1170) (f(j),j=j1,j2)
      end do
 1170 format (8f10.7)
      if(icond.ne.0) write (icond,1220) (f(j),j=1,nq)
cpp
      write (*,1220) (f(j),j=1,nq)
 1220 format (1x,10f10.5,/,(1x,10f10.5))
      do 340 i=1,nq
        f(i)=one/f(i)
  340 continue
c     end of reading the force constants
c
c     fixed coordinates - set the corresponding rows and columns
c     of the inverse force constant matrix zero
      if (nfix.ne.0) then
      do 390 k=1,nfix
         ii=ifix(k)
         f(i)=zero
  390 continue
      end if
c     end of readfc
      end
c
      subroutine relaxd(nq,fi,fc,cq,edecr,dismax,icmx,scfac)
      implicit real*8 (a-h,o-z)
c     this routine calculates the new internal coordinates
c     and the projected energy lowerings
c     in this version,. it can only use diagonal force constants for
c       relaxation
c     nq (input) number of internal coordinates
c     fi (input) internal forces
c     fc (input) inverse force constant matrix
c     cq  (output) change in the internal coordinates
c     edecr (input) decrease of the total energy in this step
c     dismax (output) maximum change in the internal coordinates
c     icmx (output) coordinate of the largest displacement
c     scfac  scale factor for the displacements, normally 1

      dimension fi(nq),fc(nq),cq(nq)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1  i,acc
      edecr=zero
      do 480 i=1,nq
         cq(i)=fc(i)*fi(i)
         edecr=edecr-cq(i)*fi(i)
 480  continue
      edecr=edecr*0.5d0
      write(6,*) 'Estimated energy decrease in aJ=',edecr
      if(edecr.gt.zero) scfac=-one
      if(scfac.ne.one) call mult(cq,scfac,nq)
      call absmax(nq,cq,icmx,dismax)
      end
c
      subroutine fixtor (should,found)
c      fix torsional angles if they turn around fully
      implicit real*8 (a-h,o-z)
      logical correct
      dimension xnfo(8)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      data xnfo/1.0d0,2.0d0,3.0d0,4.0d0,6.0d0,9.0d0,4.5d0,1.5d0/
      deviat=should-found
c       deviat=qq(i)+ccin(i)-qq1(i)
c
c      write(*,*) ' FIXTOR: should,found,diff',should,found,deviat
c
c pp
c        the purpose of the following code is to identify cases when a
c        torsional angle changes suddenly by 2*pi
c        because of the 1/n normalization, the change may be 2*pi,
c        2*pi/2, /3, /4, /6 or /9 or /4.5 or /1.5
         adevi=abs(deviat)
         if(adevi.gt.half) then
            do  1285 nf=1,8
              xnfold=xnfo(nf)
c             this is the number of torsional angles linearly combined
c             it can be 1,2,3,4,6 or 9
              offset=two*pi/xnfold
              upper=1.1d0*offset
              xlower=0.9d0*offset
              if(adevi.lt.upper.and.adevi.gt.xlower) then
                if(deviat.gt.zero) then
c                  qq1(i)=qq1(i)+offset
                   found=found+offset
                   go to 1286
                else
c                  qq1(i)=qq1(i)-offset
                   found=found-offset
                   go to 1286
                endif
              endif
 1285       continue
         endif
 1286    continue
c
      end
c
      subroutine precon_distortnew(x,dm1,y,nq)
c  calculates y=dm1*x
      implicit real*8 (a-h,o-z)
      dimension x(nq),dm1(nq),y(nq)
      do i=1,nq
        y(i)=dm1(i)*x(i)
      end do
      end
c
      subroutine congrad(nq,nek,b,ibcontr,x,y,dm1,c,t,p,r,z,oper)
c   this routine solves iteratively the equation Ax=c
c   inek,b,ibcontr,t are strictly arguments for "oper"
c   A is the operation = oper can be bbt or btb;
c   bbt forms y= (BBtranspose)x; t is a temporary variable, nek long
c   Btb forms y= Btranspose*B*x ; t is a  temp. variable
c   Algorithm describe e.g. by A. Greenbaum, C. Li and H. Z. Chao,
c   "Practical Iterative Methods for Large-Scale Computations",
c   ed.  D.L. Boley, D.G. Truhlar, Y. Saad, R.E. Wyatt, L.E. Collins
c   North-Holland, Amsterdam 1989.
c
c   The righ-hand side is c, and the final solution is x.
c   dm1 will hold the inverse diagonal of the matrix A=BBt
c
c   For general implementation, replace oneiter and its arguments
c   and also bdiag. The latter forms and inverts diag(A)**-1
c   precon_distortnew multiplies a vector with diag(A)**-1
c   y,dm1,t,p,r,z are just workspaces
c
c   Parameters: INPUT
c   nq: number of internal coordinates
c   nek: number of Cartesians (3*na)
c   b: B matrix
c   ibcontr: contraction info for the compactly stored B matrix
c   c: right hand side of the equation
c   dm1: inverse diagonal of the matrix (in this case BBt)
c      OUTPUT
c   x: solution of the equation
c   oper: the external subroutine
      implicit real*8 (a-h,o-z)
      parameter(onem=-1.0d0)
      dimension b(54,nq),ibcontr(20,nq),c(nq),t(nek),dm1(nq)
      dimension x(nq),y(nq),z(nq),p(nq),r(nq),aa(2),a(2)
      external oper
c  set up circular counter
      i0=1
      i1=2
      damp=1.0d-7
c  set the first colums of x to D**-1*c   x0=DM*c
      call precon_distortnew(c,dm1,x,nq)
c   y=A*x0
      call oper(nq,nek,b,ibcontr,x,y,t)
c   put in a little damping
      call add1(x,damp,y,nq)
c   r0=c-Ax0
      call tfer(c,r,nq)
      call add1(y,onem,r,nq)
c   z=DM*r0
      call precon_distortnew(r,dm1,z,nq)
c   p0=z
      call tfer(z,p,nq)
c   aa0=z(t)r0
      call inner(z,r,sum,nq)
      aa(1)=sum
c  return address for iteration
      iter=0
100   continue
      iter=iter+1
c   y=Ap(k-1)
        call oper(nq,nek,b,ibcontr,p,y,t)
c   put in a little damping
      call add1(p,damp,y,nq)
c   a(k-1)=aa(k-1)/p(k-1)(t)y
        call inner(p,y,sum,nq)
        a(i0)=aa(i0)/sum
c   x(k)=x(k-1)+a(k-1)p(k-1)
        call add1(p,a(i0),x,nq)
c   r(k)=r(k-1)-a(k-1)*y
        s=-a(i0)
        call add1(y,s,r,nq)
c   check for convergence
        call absmax(nq,r,imx,absmx)
        toler=1.0d-9
        if(absmx.lt.toler) go to 1000
c   z=DM*r(k)
        call precon_distortnew(r,dm1,z,nq)
c   aa(k)=z(t)*r(k);   bb=aa(k)/aa(k-1)
        call inner(z,r,sum,nq)
        aa(i1)=sum
        bb=aa(i1)/aa(i0)
c   p(k)=z+bb*p(k-1)
        call mult(p,bb,nq)
        call add(z,p,nq)
        ii=i1
        i1=i0
        i0=ii
      go to 100
 1000 continue
c     call tfer(x(1,i1),x(1,1),nq)
c      print *, 'congrad iter=',iter
      end
c=======================================================================
      subroutine distort (nek,nq,na,qc,qq,  qq1,xa,ibcode,ibcontr,b,
     1  shftpt,xy)
c-----------------------------------------------------------------------
c   this routine distorts the molecule along given internal
c   coordinates and gives the Cartesians as a result
c
c   INPUT
c   nek=3*na, (input) number of cartesians
c   na (input) number of nuclei
c   nq (input) number of internal coordinates
c   qc (input) change in the internal coordinates
c   qq (input) original internal coordinates
c   qq1 (workspace, contains the final coordinates on entry)
c   xa(3,na)(input) cartesians in angstrom units
c   ibcode(6,*): coding info for the B matrix
c   ibcontr(20,nq): contraction for the B matrix
c   b (workspace) b matrix
c   shftpt (input) see machb
c   OUTPUT
c   xy(3,na)  (output) distorted cartesians
      implicit real*8 (a-h,o-z)
      external btb,bbt
      dimension qc(nq),qq(nq),qq1(nq),xa(3,na),xy(3,na)
      dimension ibcontr(20,nq),b(54,nq)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /tapeneu/ inp,inp2,iout,ipun,ix,icond,inpf3,inpf4
      include "common/big"
c      common bl(1000)
c   this is the maximum internal coordinate change -otherwise
c   use scaling
      parameter (coormax=0.10d0)
c      parameter (cartmax=0.25d0)
c qq1 is the current set of internal coordinates, qq is the original set
c   qq1 will be replaced by the current error
      toler=1.0d-9
c   define temporary storage
      ibase=icorr(0)
      ix=icorr(nek)
      iy=icorr(nek)
      idm1=icorr(nek)
      it=icorr(nek)
      ip=icorr(nek)
      ir=icorr(nek)
      iz=icorr(nek)
      iu=icorr(nek)
c  in the redundant case, first project the internal coordinate
c  changes
c  THIS DOES NOT GO - it projects OK but that changes the other
c  directions if there is redundancy. I do not understand!
c     call bdiag(nq,nek,b,ibcontr,bl(idm1))
c     call bbt(nq,nek,b,ibcontr,qc,bl(iu),t)
c     call congrad(nq,nek,b,ibcontr,qc,bl(iy),bl(idm1),bl(iu),
c    1  bl(it),bl(ip),bl(ir),bl(iz),bbt)
c  end projecting delta q
c  this should be done only if the coordinate system is redundant
c  calculate the diagonal of the matrix BtB
      call bdiag1(nq,nek,b,ibcontr,q(idm1))
      call tfer (qq,qq1,nq)
c qc holds the desired set:
      call add(qq,qc,nq)
c xy is the current set of Cartesians
      call tfer (xa,xy,nek)
c   set the first value of the change
      call tfer(qc,q(ix),nq)
c
c     this is the return point for iteration
c
c      call outvec(xy,nek,'old cartesians')
      onem=-one
      itmax=40
      do 900 irep=1,itmax
c  calculate the error in the internal coordinates, put it in qq1
      call mult(qq1,onem,nq)
      call add(qc,qq1,nq)
c  do not calculate actual - should values
c  rather, look at the change
      call absmax(nq,q(ix),imax,delmax)
      if(delmax.lt. toler) go to 1000
      call absmax(nq,qq1,imax,delmax)
      if(delmax.gt.coormax) then
c   scale back the change
        scale=coormax/delmax
        call mult(qq1,scale,nq)
      else
        scale=one
      end if
c  delta X = Bt * (BBt)**-1 * Delta Q  to first order
c  first evaluate  x=[(BBt)**-1]Delta Q
c  Interface   subroutine congrad(nq,nek,b,ibcontr,x,y,dm1,c,t,p,r,z,op)
c   alllocate memory for the temporary stuff
c     call congrad(nq,nek,b,ibcontr,q(ix),q(iy),q(idm1),qq1,
c    1  q(it),q(ip),q(ir),q(iz))
c   bl(ix) holds (BBt)**-1*qq1
c interface subroutine btxv(nq,nek,b,ibcontr,v,w): gives w=Bt*v
c     call btxv(nq,nek,b,ibcontr,q(ix),q(iy))
c  bl(iy) now holds the first-order change in the Cartesians
c     only for checking - calc. B*q(iy)
c     call bxv(nq,nek,b,ibcontr,q(iy),q(it))
c     call add1(qq1,onem,q(it),nq)
c     call absmax(nq,q(it),ii,dmax)
c     print *, 'devi=', ii,dmax
c   end of testing
c
c   new method: x is the conj. grad solution of Bt*B*x=Bt*q
      call btxv(nq,nek,b,ibcontr,qq1,q(iu))
      call congrad(nek,nq,b,ibcontr,q(ix),q(iy),q(idm1),q(iu),
     1 q(it),q(ip),q(ir),q(iz),btb)
c   re-introduce the max. Cartesian idea - it slows things down
c     call absmax(nek,q(iy),ii,cmax)
c     if(cmax.gt.cartmax) then
c       scale=cartmax/cmax
c       call mult(q(iy),scale,nek)
c     end if
      call add(q(ix),xy,nek)
c      call outvec(xy,nek,'new cartesians')
c  calculate the "should" value of the internal coordinates
c  this is  the original value plus the scaled change
      sc1=one-one/scale
      call tfer(qc,q(ix),nq)
      call add1(qq1,sc1,q(ix),nq)
c     this call to the b matrix routine only serves to calculate the
c    coordinates
c   Interface:
c  subroutine machbnew(na,xa,nq,qonly,shftpt,ibcode,ibcontr,bmat,qq)
c   save the old values of the int. coord. in bl(ix) with negative sign
c???????????????????????????????
c      call tfer(qc,q(ix),nq)
c????????????????????????????????
      call add1(qq1,onem,q(ix),nq)
      call machbnew(na,xy,nq,.true.,shftpt,ibcode,ibcontr,b,qq1)
      do 720 i=1,nq
        call fixtor(q(ix-1+i),qq1(i))
  720 continue
c   calculate change
      call add1(qq1,onem,q(ix),nq)
c      call outvec(qq1,nq,'qnew')
c      call outvec(q(ix),nq,'delta q')
 900  continue
      if(icond.ne.0) write (icond,*) 'Subroutine distort: ',
     >            'Back transformation to cartesian coordinates:'
      if(icond.ne.0) write (icond,*)
     >                 '?Caution: within 40 steps no convergence'
c      call error('error detected','distort')
 1000 continue
c
c restore qc
      if(icond.ne.0) write(icond,*)
      if(icond.ne.0) write(icond,*) 'Subroutine distort: ',
     >               'Back transformation to cartesian coordinates:'
      if(icond.ne.0) write(icond,*) 'Convergence in ',irep,' steps'
      call add1(qq,onem,qc,nq)
      call corlsr(ibase)
      return
      end
c=======================================================================
      subroutine distortnew (nek,nq,na,qc,qq,
     &                       qq1,xa,ibcode,ibcontr,b,
     &                       ginv,shftpt,xy)
c-----------------------------------------------------------------------
c   this routine distorts the molecule along given internal
c   coordinates and gives the Cartesians as a result
c
c                  Delta X = B(t)*[G-]*Delta q
c
c   The internal coordinate changes shoul;d have been projected at this
c   point but only in first order
c
c   INPUT
c   nek=3*na, (input) number of cartesians
c   na (input) number of nuclei
c   nq (input) number of internal coordinates
c   qc (input) change in the internal coordinates, scaled if necessary
c   qq (input) original internal coordinates
c   qq1 (workspace, contains the final coordinates on entry)
c   xa(3,na)(input) cartesians in angstrom units
c   ibcode(6,*): coding info for the B matrix
c   ibcontr(20,nq): contraction for the B matrix
c   b (workspace) b matrix
c   ginv: (generalized) inverse of the G matrix
c   shftpt (input) see machb
c   OUTPUT
c   xy(3,na)  (output) distorted cartesians  iIN ATOMIC UNITS
      implicit real*8 (a-h,o-z)
      dimension qc(nq),qq(nq),qq1(nq),xa(3,na),xy(3,na),ginv(nq,nq)
      dimension ibcontr(20,nq),b(54,nq)
c      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
c     1i,acc
      common /tapeneu/ inp,inp2,iout,ipun,ix,icond,inpf3,inpf4
      include "common/big"
c   this is the maximum internal coordinate change -otherwise use scaling
      parameter (coormax=0.20d0,zero=0.0d0,half=0.5d0,one=1.0d0,
     &           pi=3.14159 26535 89793 24d0)
c      parameter (cartmax=0.25d0)
c qq1 is the desired set of internal coordinates, qq is the original set
      shftpt=half*pi
      call absmax(nq,qc,imax,chmax)
      if(chmax.gt.coormax) then
        scale=coormax/chmax
      write(*,*) 'Max. change, max. allowed change, scale factor',
     &            chmax,coormax,scale
        call mult(qc,scale,nq)
      end if
c   define temporary storage
      ibase=icorr(0)
      memvec=icorr(nq)
      ix=icorr(nek)
cpp
c      write(*,*) 'Changes in the internal coordinates'
c      write(*,'(5f12.6)') (qc(k),k=1,nq)
c
CPP
c      call prntmatn(6,nq,nq,nq,ginv,'G inverse on entry to distortnew')
      toler=1.0d-8
      onem=-one
      itmax=40
      call dcopy_x(nq,qc,1,qq1,1)
      call dcopy_x(nek,xa,1,xy,1)
c      write(*,*) 'Original cartesians'
c      do i=1,na
c        write(*,200) i,xy(1,i),xy(2,i),xy(3,i)
  200   format(i5,3f12.6)
c      end do
c
      call machbnew(na,xy,nq,.false.,shftpt,ibcode,ibcontr,b,q(memvec))
cpp
c      write(*,'("int. coord.",/,(4(i5,f12.6)))')(k,q(memvec+k-1),k=1,nq)
      residuum=100.0d0
      do 900 irep=1,itmax
      call dgemv_x('n',nq,nq,one,ginv,nq,  qq1,1,zero,q(memvec),1)
cpp
c      write(*,'("qq1",(5f12.6))') (qq1(k),k=1,nq)
c      write(*,*) 'G-*Delta(q)'
c      write(*,*) (q(memvec+kk-1),kk=1,nq)
c      write(*,*) 'B matrix on entry in btxv'
c      do i=1,nq
c         write(*,950) i,(b(kk,i),kk=1,12) 
c  950  format(i5, 12f9.5)
c      end do
c      write(*,*) 'ibcontr on entry in btxv'
c      do i=1,nq
c        natms=ibcontr(2,i)
c        write(*,'(i4,2x,18i4)') natms,(ibcontr(2+k,i),k=1,natms)
c      end do
cpp
      call btxv(nq,nek,b,ibcontr,q(memvec),q(ix))
cpp
c         write(*,*) 'Delta Cartesians' 
c        write(*,'(3f12.6)') (q(ix+k-1),k=1,nek)
cpp
      call add(q(ix),xy,nek)
cpp
c      write(*,*) 'New Cartesians'
c      do i=1,na
c        write(*,200) i,xy(1,i),xy(2,i),xy(3,i)
c      end do
c  Calculate new internal coordinates       
c   Interface:
c   call: machbnew(na,xa,nq,qonly,shftpt,ibcode,ibcontr,bmat,qq)
      call machbnew(na,xy,nq,.true.,shftpt,ibcode,ibcontr,b,q(memvec))
cpp
c      write(*,*) 'New internal coordinates, iter=',irep
c      write(*,250) (kk,q(memvec+kk-1),kk=1,nq)
c      write(*,*) 'Original internal coordinates, iter=',irep
c      write(*,250) (kk,qq(kk),kk=1,nq)
c  Calculate new changes
      call dcopy_x(nq,qq,1,qq1,1)
      call add(qc,qq1,nq)
      call add1(q(memvec),-one,qq1,nq)
cpp
c      write(*,*) 'Required changes in internal coordinates, iter=',irep
c      write(*,250) (kk,qq1(kk),kk=1,nq)
 250  format(i5,f14.6)

c  qq1 should have the error in the internal coordinates
c  There may be a problem with redundant coordinates because the residuum
c  cannot be exactly zero. However, in the final stages, when the geom.
c  distortion is small, the remaining residuum should go to zero
      call absmax(nq,qq1,imax,xmax)
      if(xmax.le.toler) exit
c      write(*,*) 'Previous and current residuum',residuum,xmax
c  Exit if the residuum increases
      if(residuum-xmax.gt.toler) then
        residuum=xmax
      else
        exit
      end if
c      call outvec(qq1,nq,'qnew')
c      call outvec(q(ix),nq,'delta q')
       if(irep.eq.itmax) then
         write(*,*) 'In distort within',itmax,' steps no convergence'
       end if
 900  continue
cpp
c      write(*,*) 'Converged internal coordinates, final residue=',
c     &            residuum
c      do i=1,nq
c        write(*,'(i5,f12.6)') i,q(memvec+i-1)
c      end do
cpp 
      call corlsr(ibase)
c  Quit if the residuum is too big
      if(xmax.gt.0.2d0) then
        write(*,*) 'Residuum is too big', residuum
        call fehler
      end if
c      toang=one/get_scalar('TOANG')
c      call dscal_x(nek,toang,xy,1)
cpp
c      write(*,*) 'New Cartesians in Bohrs'
c      do i=1,na
c        write(*,200) i,xy(1,i),xy(2,i),xy(3,i)
c      end do
      end
c=======================================================================
      subroutine distortnew1(nek,nq,na,qc,qq,
     &                       qq1,xa,ibcode,ibcontr,b,
     &                       ginv,xy)
c-----------------------------------------------------------------------
c   this routine distorts the molecule along given internal
c   coordinates and gives the Cartesians as a result
c
c                  Delta X = B(t)*[G-]*Delta q
c
c   The internal coordinate changes shoul;d have been projected at this
c   point but only in first order
c
c   INTENT(IN)
c   nek=3*na, (input) number of cartesians
c   nq (input) number of internal coordinates
c   na (input) number of nuclei
c   qc (input) change in the internal coordinates, scaled if necessary
c   qq (input) original internal coordinates
c   qq1 (workspace, contains the final coordinates on entry)
c   xa(3,na)(input) cartesians in angstrom units
c   ibcode(6,*): coding info for the B matrix
c   ibcontr(20,nq): contraction for the B matrix
c   b (workspace) b matrix
c   ginv: (generalized) inverse of the G matrix
c
c  COMMONS:
c   inpf3  Unit number of the file containing the internal cordinate definitions and values
c   inpf4  Unit number of the second file containing the internal cordinate definitions and
c   values. inpf3 and inpf4 are swapped. The previuos file is used to correct the primitive
c   torsion values by 2*Pi in intcvalues
c
c   INTENT(OUT)
c   xy(3,na)  (output) distorted cartesians  iIN ATOMIC UNITS
      implicit real*8 (a-h,o-z)
      dimension qc(nq),qq(nq),qq1(nq),xa(3,na),xy(3,na),ginv(nq,nq)
      dimension ibcontr(20,nq),b(54,nq)
c      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
c     1i,acc
      common /tapeneu/ inp,inp2,iout,ipun,ix,icond,inpf3,inpf4
      include "common/big"
c   this is the maximum internal coordinate change -otherwise use scaling
      parameter (coormax=0.10d0,zero=0.0d0,one=1.0d0,
     &           pi=3.14159 26535 89793 24d0)
c      parameter (cartmax=0.15d0)
c qq1 is the desired set of internal coordinates, qq is the original set
      shftpt=0.5d0*pi
      call absmax(nq,qc,imax,chmax)
      if(chmax.gt.coormax) then
        scale=coormax/chmax
      write(*,*) 'Max. change, max. allowed change, scale factor',
     &            chmax,coormax,scale
        call mult(qc,scale,nq)
      end if
c   define temporary storage
      ibase=icorr(0)
      memvec=icorr(nq)
cpp
      memvec1=icorr(nq)
      ix=icorr(nek)
c
CPP
c      call prntmatn(6,nq,nq,nq,ginv,'G inverse on entry to distortnew')
      toler=1.0d-8
      onem=-one
      itmax=40
      call dcopy_x(nq,qc,1,qq1,1)
      call dcopy_x(nek,xa,1,xy,1)
c      write(*,*) 'Original cartesians'
c      do i=1,na
c        write(*,200) i,xy(1,i),xy(2,i),xy(3,i)
  200   format(i5,3f12.6)
c      end do
c
      call machbnew(na,xy,nq,.false.,shftpt,ibcode,ibcontr,b,q(memvec1))
c  Calculate the new internal coordinates using intcvalues
       
      call intcvalues(na,xa,shftpt,nq,q(memvec))
cpp  compare the internal coordinates obtained by BMAT, intcvalues and BMAT again
c      write(*,*) 'Composite internal coordinates: BMAT,intcvalues'
c      do k=1,nq
c        write(*,'(i4,3F12.6)') 
c     &          k,q(memvec1+k-1),q(memvec+k-1)
c      end do
      residuum=100.0d0
      do 900 irep=1,itmax
      call dgemv_x('n',nq,nq,one,ginv,nq,  qq1,1,zero,q(memvec),1)
cpp
c      write(*,*) 'G-*Delta(q)'
c      write(*,*) (q(memvec+kk-1),kk=1,nq)
cpp
c      write(*,*) 'B matrix on entry in btxv'
c      do i=1,nq
c         write(*,950) i,(b(kk,i),kk=1,12) 
c  950  format(i5, 12f9.5)
c      end do
c      write(*,*) 'ibcontr on entry in btxv'
c      do i=1,nq
c        natms=ibcontr(2,i)
c        write(*,'(i4,2x,18i4)') natms,(ibcontr(2+k,i),k=1,natms)
c      end do
cpp
      call btxv(nq,nek,b,ibcontr,q(memvec),q(ix))
cpp
c        write(*,'(3f12.6)') (q(ix+k-1),k=1,nek)
cpp
      call add(q(ix),xy,nek)
cpp
c      write(*,*) 'New Cartesians'
c      do i=1,na
c        write(*,200) i,xy(1,i),xy(2,i),xy(3,i)
c      end do
c  Calculate new internal coordinates       
c   Interface:
c   call: machbnew(na,xa,nq,qonly,shftpt,ibcode,ibcontr,bmat,qq)
c      call machbnew(na,xy,nq,.true.,shftpt,ibcode,ibcontr,b,q(memvec1))
      call intcvalues(na,xy,shftpt,nq,q(memvec))
c      write(*,*) 'Composite internal coordinates: BMAT,intcvalues'
c      do k=1,nq
c        write(*,'(i4,3F12.6)') 
c     &          k,q(memvec1+k-1),q(memvec+k-1)
c      end do
cpp
c      write(*,*) 'New internal coordinates, iter=',irep
c      write(*,250) (kk,q(memvec+kk-1),kk=1,nq)
c      write(*,*) 'Original internal coordinates, iter=',irep
c      write(*,250) (kk,qq(kk),kk=1,nq)
c  Calculate new changes
      call dcopy_x(nq,qq,1,qq1,1)
      call add(qc,qq1,nq)
      call add1(q(memvec),-one,qq1,nq)
cpp
c      write(*,*) '    Original   , change,      current  , iter=',irep
c      do k=1,nq
c        write(*,'(i4,3f12.6)') k,qq(k),qc(k),q(memvec+k-1)
c      end do
cpp

c  qq1 should have the error in the internal coordinates
c  There may be a problem with redundant coordinates because the residuum
c  cannot be exactly zero. However, in the final stages, when the geom.
c  distortion is small, the remaining residuum should go to zero
      call absmax(nq,qq1,imax,xmax)
c      write(*,*) 'Previous and current residuum',residuum,xmax
      if(xmax.le.toler) exit
c  Exit if the residuum increases
      if(residuum-xmax.gt.toler) then
        residuum=xmax
      else
        exit
      end if
c      call outvec(qq1,nq,'qnew')
c      call outvec(q(ix),nq,'delta q')
       if(irep.eq.itmax) then
         write(*,*) 'In distort within',itmax,' steps no convergence'
       end if
 900  continue
cpp
c      write(*,*) 'Converged internal coordinates, final residue=',
c     &            residuum
c      do i=1,nq
c        write(*,'(i5,f12.6)') i,q(memvec+i-1)
c      end do
cpp 
      call corlsr(ibase)
c  Quit if the residuum is too big
      if(xnax.gt.0.2d0) then
        write(*,*) 'Residuum is too big', xmax
        call fehler
      end if
      toang=one/get_scalar('TOANG')
      call dscal_x(nek,toang,xy,1)
cpp
c      write(*,*) 'New Cartesians in Bohrs'
c      do i=1,na
c        write(*,200) i,xy(1,i),xy(2,i),xy(3,i)
c      end do
      end
c=======================================================================
      subroutine btb(nek,nq,b,ibcontr,x,y,t)
c   calculates y=bt*b*x
      implicit real*8 (a-h,o-z)
      dimension x(nek),y(nek),t(nq)
      call bxv(nq,nek,b,ibcontr,x,t)
      call btxv(nq,nek,b,ibcontr,t,y)
      end
c=======================================================================
      subroutine btxvnew(nq,nek,b,ibcontr,x,y)
      implicit real*8 (a-h,o-z)
      dimension b(54,nq),ibcontr(20,nq),x(nq),y(nek)
      do ii=1,nq
        natms=ibcontr(2,ii)
        do iat=1,natms
          kat=ibcontr(2+iat,ii)
          do k=-2,0
            y(kat*3+k)=y(kat*3+k)+b(iat*3+k,ii)*x(ii)
          end do
        end do
      end do
c
      end
c=======================================================================
