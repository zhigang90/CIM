c======================================================================
      subroutine f1shift_xyz(fock1,ncf,ntri,lind,dens0,over0,over1,
     *                       ds1d,bl,xlvsh)
c-----------------------------------------------------------------
c Calculates
c          -1/4 * xlvsh * S0*( D0*S1*D0 )*S0
c needed for F1const with level shift and adds it to Fock1
c
c It does it for three directions at once (x,y,z)
c-----------------------------------------------------------------
c INPUT :
c
c fock1    - constant part of F1=H1+G(D0,g1)
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c dens0    - unpreturbed density
c over0    - unpreturbed overlap
c over1    - Ist-order overlap matrix
c ds1d     - scratch for d0*S1*d0 matrix
c bl       - storage for everything
c xlvsh    - level shift
c
c Output :
c
c fock1    - constant part of F1=H1+G(D0,g1) - 0.5*xlvsh*S0*[0.5*D0*S1*D0]*S0
c
c-----------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension bl(*)
      dimension fock1(ntri,3)
      dimension dens0(*)
      dimension over0(*)
      dimension over1(ntri,3)
      dimension ds1d(ntri,3)
c--------------------------------------------------------------
c      return
c--------------------------------------------------------------
c - 0.5 D0*S1*D0
c
      call getmem(ncf*ncf*3,lwork3)
c
      call zeroit(ds1d,ntri*3)
      call ds1d_xyz(bl(lwork3),ncf,ntri,lind,dens0, over1,ds1d)
c---------------------------------------------------------------
c        write(6,*)' Fock1 on entry  '
c        call drumh(Fock1(1,1),ncf, 6  ,'Fock1-x ')
c        call drumh(Fock1(1,2),ncf, 6  ,'Fock1-y ')
c        call drumh(Fock1(1,3),ncf, 6  ,'Fock1-z ')
c---------------------------------------------------------------
c we need to add -1/4*xlvsh*S( D S1 D ) S  to  Fock1=H1 + G(D0,g1)
c change sign of -1/2DS1D to +
c
      factor=-1.d0
      call dscal(ntri*3,factor,ds1d,1)
c
c - 0.5 S0*( +0.5*D0*S1*D0 )*S0* xlvsh  and add it to f1const :
c
      call  sd1s_xyz(bl(lwork3),ncf,ntri,lind,over0,ds1d,xlvsh,fock1)
c---------------------------------------------------------------
c        write(6,*)' Fock1 on exit   '
c        call drumh(Fock1(1,1),ncf, 6  ,'Fock1-x ')
c        call drumh(Fock1(1,2),ncf, 6  ,'Fock1-y ')
c        call drumh(Fock1(1,3),ncf, 6  ,'Fock1-z ')
c---------------------------------------------------------------
c
      call retmem(1)
c--------------------------------------------------------------
      end
c======================================================================
      subroutine ds1d_xyz(work,ncf,ntri,lind,den,ove,hfc)
      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension den(*),ove(*), hfc(*)
      dimension work(3,ncf,ncf)
      data zero,half /0.d0, 0.5d0/
c
      ntr2=ntri*2
c---------------------------------------------
c  calculate  the -0.5*d0*s1*d0  matrix
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               kj=jj+k
               if(k.gt.j) kj=kk+j
               sux=sux+ove(ik)*den(kj)
               suy=suy+ove(ntri+ik)*den(kj)
               suz=suz+ove(ntr2+ik)*den(kj)
            enddo
            work(1,i,j)=sux
            work(2,i,j)=suy
            work(3,i,j)=suz
         enddo
      enddo
c
      ij=0
      do i=1,ncf
         ii=lind(i)
         do j=1,i
            ij=ij+1
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               sux=sux+den(ik)*work(1,k,j)
               suy=suy+den(ik)*work(2,k,j)
               suz=suz+den(ik)*work(3,k,j)
            enddo
c
            hfc(ij)     =hfc(ij)      - half*sux
            hfc(ntri+ij)=hfc(ntri+ij) - half*suy
            hfc(ntr2+ij)=hfc(ntr2+ij) - half*suz
c
         enddo
      enddo
c
      end
c======================================================================
c it calculates S0*D1*S0 as needed for level shifted cphf :
c
      subroutine sd1s_xyz(work,ncf,ntri,lind,over0,dens1,xlvsh,hfc)
      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension over0(*),dens1(*), hfc(*)
      dimension work(3,ncf,ncf)
      data zero,half /0.d0, 0.5d0/
c---------------------------------------------
      ntr2=ntri*2
c---------------------------------------------
      halfx= half*xlvsh
c---------------------------------------------
c  calculate  the  0.5*s0*d1*s0  matrix
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               kj=jj+k
               if(k.gt.j) kj=kk+j
               sux=sux+dens1(ik)*over0(kj)
               suy=suy+dens1(ntri+ik)*over0(kj)
               suz=suz+dens1(ntr2+ik)*over0(kj)
            enddo
            work(1,i,j)=sux
            work(2,i,j)=suy
            work(3,i,j)=suz
         enddo
      enddo
c
      ij=0
      do i=1,ncf
         ii=lind(i)
         do j=1,i
            ij=ij+1
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               sux=sux+over0(ik)*work(1,k,j)
               suy=suy+over0(ik)*work(2,k,j)
               suz=suz+over0(ik)*work(3,k,j)
            enddo
c
            hfc(ij)     =hfc(ij)      -halfx*sux
            hfc(ntri+ij)=hfc(ntri+ij) -halfx*suy
            hfc(ntr2+ij)=hfc(ntr2+ij) -halfx*suz
c
         enddo
      enddo
c
      end
c======================================================================
      subroutine cphf_enx(r3,thrs,ntri,mgo,liter,cpuit,elait,lend,
     *                    errmax,residx,residy,residz,oscylates)
      implicit real*8 (a-h,o-z)
      logical oscylates
      dimension r3(*)     ! last residium
      data zero /0.d0/
c
      call getival('iout',iout)
c
      deltx=zero
      delty=zero
      deltz=zero
      do l=1,ntri
         deltax=abs(r3(l))
         deltay=abs(r3(l+ntri))
         deltaz=abs(r3(l+ntri*2))
         if(deltax.gt.deltx) deltx=deltax
         if(deltay.gt.delty) delty=deltay
         if(deltaz.gt.deltz) deltz=deltaz
      enddo
c
c allow mixing directions ONLY if loose integ.thres.
c and it is not converged
c
c     mix_dir=1    ! positive value means mix directions
c     if(mgo.eq.1) then
c        if(deltx.le.thrs) mix_dir=-1
c        if(delty.le.thrs) mix_dir=-1
c        if(deltz.le.thrs) mix_dir=-1
c     else
c        mix_dir=-1
c     endif
c
c     call setival('mix_dir',mix_dir)
c
      oscylates=.false.
      if(deltx.gt.residx .and. deltx.gt.thrs) oscylates=.true.
      if(delty.gt.residy .and. delty.gt.thrs) oscylates=.true.
      if(deltz.gt.residz .and. deltz.gt.thrs) oscylates=.true.
c
      lend=1
      errmax=max(deltx,delty,deltz)
      if(errmax.gt.thrs) lend=0
c
      cput=cpuit/60.d0
      elat=elait/60.d0
c
      peltx=max(thrs,deltx)
      pelty=max(thrs,delty)
      peltz=max(thrs,deltz)
c
      write(iout,420) liter,peltx,pelty,peltz,cput,elat,oscylates
      call f_lush(iout)
c 420 format(i3,3x,3(1x,1pe11.4),1x,2(f8.2,1x),l5)
  420 format(i3,3x,3(1x,e11.4),1x,2(f8.2,1x),l5)
c
      residx=deltx
      residy=delty
      residz=deltz
c
      end
c======================================================================
      subroutine d1const_xyz(rhf,   ncf,   ntri,  nocc,  lind,
     *                       vec,   val,   d0,    f1,    s1,
     *                       d1c,   ivcd,  pve,   bl)
c-----------------------------------------------------------------
c Calculates the constant part of the D1 matrix for :
c
c  (1) electric field perturbation (dipole polarizability)
c  (2) nuclear displacement perturbation (analytical hessian)
c
c It does it for three directions at once (x,y,z)
c-----------------------------------------------------------------
c INPUT :
c
c rhf      - rhf or uhf flag
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c nocc     - number of occupied orbitals
c lind     - diagonals of i*(i-1)/2
c vec      - eigen vectors
c val      - eigen values
c d0       - unpreturbed density
c f1       - Ist-order fock matrix
c s1       - Ist-order overlap matrix
c ivcd     - integer flag for VCD  0 - NO
c bl       - storage for everything
c
c OUTPUT :
c
c d1c      - constant part of D1 :
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eo*S1)Cv*[CoCv+ + CvCo+]
c pve      - perturbed coefficient matrix for VCD
c-----------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      dimension lind(*)
      dimension bl(*)
      dimension d0(*), vec(*), val(*)
      dimension s1(ntri,3), f1(ntri,3)
      dimension d1c(ntri,3), pve(ncf*nocc,3)          ! output
c--------------------------------------------------------------
c    Calculate the constant part of D10
c   ----------------------------------------
c 1. contributions from [ F1(d0,g1) - ei S1 ]
c 2. contributions from the 0.5 D0*S1*D0
c
c  S1 is not zero only for field dependent basis set
c
c The total constant part of the D1 stored in d1con(*)
c
c--------------------------------------------------------------
cc      call getrval('xlvsh',xlvsh)      ! disabled    JB
cc - NOTE: If this is reenabled, it should be passed as an
cc   argument and NOT obtained from the depository - this
cc   routine is called on the slaves where depository
cc   stuff does NOT work    <-----   IMPORTANT
      xlvsh = 0.0d0
c--------------------------------------------------------------
      fact1=2.0d0     ! for dens1_part1 like occupancy number
      fact2=0.5d0     ! for dS1d
      if(.not.rhf) then
        fact1=1.0d0
        fact2=1.0d0
      endif
c--------------------------------------------------------------
c 1. calculate contributions from [ F1(d0,g1) - ei S1 ]
c
c one direction at the time :
c
      call getmem(ncf**2,lw1)
      nvirt=ncf-nocc
      memvirt = ncf*MAX(nvirt,nocc)
      call getmem(memvirt,lw2)
      call getmem(ncf*nocc,lw3)
c
      do icr=1,3
        call dens1_1dir2n(fact1,ncf, nocc, xlvsh, f1(1,icr),
     *                    s1(1,icr),vec, val, d1c(1,icr), bl(lw1),
     *                    bl(lw2),bl(lw3))
c       call printd1('new',ncf,d1c(1,icr))
c part of the pert. coeff. matrix related to this part is in lw2
        If(ivcd.NE.0) call dcopy(ncf*nocc,bl(lw2),1,pve(1,icr),1)
      enddo
c
      call retmem(3)
c
c  2. calculate -fact2*D0*S1*D0 and add it to the constant part of D10.
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     call getmem(ncf*ncf,lw3)
c
c     do icr=1,3
c     call ds1d_1dir(fact2,bl(lw3),ncf,lind,d0,s1(1,icr),d1c(1,icr))
c     enddo
c
c     call retmem(1)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      call getmem(ncf*ncf,lw0)
      call getmem(ncf*ncf,lw1)
      call getmem(ncf*ncf,lw2)
c
      do icr=1,3
        call ds1d_1m(fact2,bl(lw0),bl(lw1),bl(lw2),ncf,
     *               d0,s1(1,icr),d1c(1,icr))
c other part of the pert. coeff. matrix
        if(ivcd.NE.0)
     *     call dgemm('n','n',ncf,nocc,ncf,-0.25d0,bl(lw2),ncf,
     *                vec,ncf,1.0d0,pve(1,icr),ncf)
      enddo
c
      call retmem(3)
c---------------------------------------------------------------
cc      write(6,*)' D1const  exit   '
cc      call drumh(d1c(1,1),ncf, 6  ,'D1con-x ')
cc      call drumh(d1c(1,2),ncf, 6  ,'D1con-y ')
cc      call drumh(d1c(1,3),ncf, 6  ,'D1con-z ')
c--------------------------------------------------------------
      end
c======================================================================
c new ds1d_1dir routine named ds1d_1m (matrix) :
c======================================================================
ccc   subroutine ds1d_1m(w0,w1,w2,ncf,d0,s1,d1c)
      subroutine ds1d_1m(fact,w0,w1,w2,ncf,d0,s1,d1c)
      implicit real*8 (a-h,o-z)
c
      dimension d0(*),s1(*)
      dimension w0(ncf,ncf),w1(ncf,ncf),w2(ncf,ncf)
      dimension d1c(*)         ! inp/out
c
      data zero,half,one /0.d0, 0.5d0, 1.0d0/
c
c  expand the dens and overlap1 matrices to quadratic
c
      call quad(d0,w0, one,ncf)
      call quad(s1,w1, one,ncf)
c  calculate W2=D0*S1
      call dgemm('n','n',ncf,ncf,ncf,
     1           one,w0 , ncf, w1, ncf,
     2           zero, w2, ncf)
c  calculate w1=0.5*W2*D0=0.5*D0*S1*D0   RHF
c  calculate w1=1.0*W2*D0=1.0*D0*S1*D0   UHF
      call dgemm('n','n',ncf,ncf,ncf,
     1           fact,w2,ncf,w0,ncf,
ccc  1           half,w2,ncf,w0,ncf,
     2           zero, w1, ncf)
c
      ij=0
      do i=1,ncf
         do j=1,i
            ij=ij+1
            d1c(ij)=d1c(ij)-w1(j,i)
         enddo
      enddo
c
      end
c======================================================================
      subroutine ds1d_1dir(fact,work,ncf,lind,den,ove,hfc)
c
c one-dimensional D S1 D constructor
c
      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension den(*),ove(*),hfc(*)      ! all ntri dimension
      dimension work(ncf,ncf)
      parameter (zero=0.0d0)
c-----------------------------------------------------------------
c  calculate  the -0.5*d0*s1*d0  matrix       for RHF
c   or
c  calculate  the -1.0*d0*s1*d0  matrix       for UHF
c
         do i=1,ncf
            ii=lind(i)
            do j=1,ncf
               jj=lind(j)
               sux=zero
               do k=1,ncf
                  kk=lind(k)
                  ik=ii+k
                  if(k.gt.i) ik=kk+i
                  kj=jj+k
                  if(k.gt.j) kj=kk+j
                  sux=sux+ove(ik)*den(kj)
               enddo
               work(i,j)=sux
            enddo
         enddo
c
         ij=0
         do i=1,ncf
            ii=lind(i)
            do j=1,i
               ij=ij+1
               sux=zero
               do k=1,ncf
                  kk=lind(k)
                  ik=ii+k
                  if(k.gt.i) ik=kk+i
                  sux=sux+den(ik)*work(k,j)
               enddo
               hfc(ij)=hfc(ij) - fact*sux
            enddo
         enddo
c
      end
c====================================================================
      subroutine dens1_1dir2(fact ,ncf,ntri,nocc,lind,
     *                       xlvsh,foc,ove ,vec ,val ,
     *                       dn1  ,w1 ,w2  )
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,two=2.d0)
      dimension lind(*)
      dimension foc(*),ove(*),vec(*),val(*),dn1(*)
      dimension w1(*),w2(*)
c-----------------------------------------------------------------------
c one direction at the time
c-----------------------------------------------------------------------
c This routine calculates a part of the D1 matrix :
c
c D1 = SUM(jocc, kvir) { Cj+*G(D1,g0)*Ck/(ej-ek) *[ Cj*Ck+ + Ck*Cj+] }
c
c Projection of F(D0,g1)-ej*S1 called from dens1_const
c-----------------------------------------------------------------------
c  fact  - occupancy : 2 for RHF, 1 for UHF
c  xlvsh - level shift
c  foc - Fock matrix Ist-order
c  ove - overlap matrix
c  vec - eigenvectors
c  val - eigenvalues
c  w1,w2 -working arrays
c
c  Output :
c  dn1 - density Ist-order (part of it)
c---------------------------------------------------
      call zeroit(dn1,ntri)
c----------------------------
      do 900 jor=1,nocc
         ej=val(jor)
         je=jor*ncf
         ja=je-ncf
         do i=1,ncf
            ii=lind(i)
            sux=zero
            do j=1,ncf
               ij=ii+j
               if(j.gt.i) ij=lind(j)+i
               sux=sux+vec(ja+j)*(foc(ij)-ej*ove(ij))
            enddo
            w1(i)     =sux
         enddo
c        -----------------------------------------
c        at this point we have the Cj+*F vector stored
c        in w1 for a current j-orbital (occupied)
c        -----------------------------------------
         call zeroit(w2,ncf)
c        -----------------------------------------
         do 1000 kor=nocc+1,ncf
            ek=val(kor)
            ejk=fact/(ej-ek-xlvsh)
            ke=kor*ncf
            ka=ke-ncf
c---
            sux=ddot(ncf,w1(1)     ,1,vec(ka+1),1)
c        -----------------------------------------
c        scalars for each pair of jorb & korb ( occ & virt)
c
            sjkx=+sux*ejk
c
c---
c           do 1150 i=1,ncf
c              w2(i)     = w2(i)     +sjkx*vec(ka+i)
c1150       continue
c---
            call daxpy(ncf,sjkx,vec(ka+1),1,w2(1),1)
c
c           in w2() virtuals rescaled by Cj+*A*Ck/(ej-ek)
c
 1000    continue
c
         do i=1,ncf
            ii=lind(i)
            vecjai=vec(ja+i)
            do j=1,ncf
               cvx=vecjai*w2(j)
               ij=ii+j
               if(j.gt.i) ij=lind(j)+i
               dn1(ij)      =dn1(ij)       + cvx
            enddo
         enddo
  900 continue
c
c diagonal elements need to be multiplied by 2
c
      do i=1,ncf
         ii=lind(i)+i
         dn1(ii)=2.0d0*dn1(ii)
      enddo
c
      end
c======================================================================
      subroutine dens1_1dir1(fact ,ncf,ntri,nocc,lind,
     *                       xlvsh,foc,vec ,val ,dn1 ,
     *                       w1 ,w2  )
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,two=2.d0)
      dimension lind(*)
      dimension foc(*),vec(*),val(*),dn1(*)
      dimension w1(*),w2(*)
c-----------------------------------------------------------------------
c one direction at the time
c-----------------------------------------------------------------------
c This routine calculates a part of the D1 matrix :
c
c D1 = SUM(jocc, kvir) { Cj+*G(D1,g0)*Ck/(ej-ek) *[ Cj*Ck+ + Ck*Cj+] }
c
c projection of F(D1,g0)=G(D1,g0) called from cphf_solver
c-----------------------------------------------------------------------
c  fact  - occupancy : 2 for RHF, 1 for UHF
c  xlvsh - level shift
c  foc - Fock matrix Ist-order
c  vec - eigenvectors
c  val - eigenvalues
c  w1,w2 -working arrays
c
c  Output :
c  dn1 -  part of the Ist-order Density
c---------------------------------------------------
      call zeroit(dn1,ntri)
c----------------------------
      do 900 jor=1,nocc
         ej=val(jor)
         je=jor*ncf
         ja=je-ncf
         do i=1,ncf
            ii=lind(i)
            sux=zero
            do j=1,ncf
               ij=ii+j
               if(j.gt.i) ij=lind(j)+i
               sux=sux + vec(ja+j)*foc(ij)
            enddo
            w1(i)     =sux
         enddo
c        -----------------------------------------
c        at this point we have the Cj+*F vector stored
c        in w1 for a current j-orbital (occupied)
c        -----------------------------------------
         call zeroit(w2,ncf)
c        -----------------------------------------
         do 1000 kor=nocc+1,ncf
            ek=val(kor)
            ejk=fact/(ej-ek-xlvsh)
            ke=kor*ncf
            ka=ke-ncf
c---
            sux=ddot(ncf,w1(1)     ,1,vec(ka+1),1)
c        -----------------------------------------
c        scalars for each pair of jorb & korb ( occ & virt)
c
            sjkx=+sux*ejk
c
c---
c           do 1150 i=1,ncf
c              w2(i)     = w2(i)     +sjkx*vec(ka+i)
c1150       continue
c---
            call daxpy(ncf,sjkx,vec(ka+1),1,w2(1),1)
c
c           in w2() virtuals rescaled by Cj+*G(D1,g0)*Ck/(ej-ek)
c
 1000    continue
c
         do i=1,ncf
            ii=lind(i)
            vecjai=vec(ja+i)
            do j=1,ncf
               cvx=vecjai*w2(j)
               ij=ii+j
               if(j.gt.i) ij=lind(j)+i
               dn1(ij)      =dn1(ij)       + cvx
            enddo
         enddo
  900 continue
c
c diagonal elements need to be multiplied by 2
c
      do i=1,ncf
         ii=lind(i)+i
         dn1(ii)      =2.0d0*dn1(ii)
      enddo
c
c      call printd1('old',ncf,dn1)
      end
c======================================================================
      subroutine dens1_1dir1n(fact ,ncf,nocc, xlv,  fock,
     1                       cvec ,eval ,d1, w1,   w2)
c-----------------------------------------------------------------------
c one direction at the time
c-----------------------------------------------------------------------
c This routine calculates a part of the first-order density matrix :
c
c D1 = SUM(jocc,avir){Cj(T)*F(D1,g0)*Ca/(ej-ea-xlv)*[Cj*Ca(T)+Ca*Cj(T)]}
c
c  (T) = transpose
c-----------------------------------------------------------------------
c  INTENT(IN)
c  fact    = occupancy : 2 for RHF, 1 for UHF
c  ncf     = number of contracted basis functions
c  nocc    = number of occupied orbitals
c  xlv     = level shift  (see below)
c  fock    = first-order Fock matrix ,F(D1,g0) above, in triangular form
c  cvec    = MO coefficients, C, the first nocc columns are occupied
c  eval    = orbital energies e above
c  INTENT(OUT)
c  d1      = the above part of the first-order density
c  STORAGE
c  w1      = a working array ncf**2 long
c  w2      = a working array ncf*nvirt long  (nvirt=ncf-nocc)
c---------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.0d0,two=2.d0)
      dimension fock(*),cvec(ncf,ncf),eval(*),d1(*)
      dimension w1(*),w2(*)
c
      nvirt=ncf-nocc
c  expand the Fock matrix to quadratic
      call quad(fock,w1,one,ncf)
c  W2=Cocc(T)*F
      call dgemm('t','n',nocc,ncf,ncf,
     1           one, cvec, ncf, w1, ncf,
     2           zero, w2, nocc)
c  W1=Cocc(T)*F*Cvirt=W2*Cvirt  nocc x nvirt matrix
      call dgemm('n','n', nocc, nvirt,ncf,
     1            one, w2, nocc, cvec(1,nocc+1), ncf,
     2            zero,w1, nocc)
c  Now scale W1=Cocc(T)FCvirt with 1.0/(e(i)-e(a)-xlv)
      call scalebydenom(nocc,nvirt,eval,eval(nocc+1),xlv,w1)
c  Calculate Cvirt*W1(T). Multiply with the factor(occupancy)
c  Result is the perturbed coeff mx in W2!      
      call dgemm('n','t', ncf, nocc, nvirt,
     1            one, cvec(1,nocc+1), ncf, w1, nocc,
     2            zero,w2, ncf)
c Calculate W2*Cocc(T). Factor moved to this call by GM
      call dgemm('n','t',ncf,ncf,nocc,
     1            fact,w2,ncf,cvec,ncf,
     2            zero, w1,ncf)
c  Result, in quadratic form, is in W1.
c  Add transpose to W1 and transfer it to the triangular array d1
      call symtrian(ncf,w1,d1)
c      call printd1('new',ncf,d1)
      end
c======================================================================
      subroutine dens1_1dir2n(fact ,ncf,nocc, xlv,  fock,
     1                       smat, cvec ,eval ,d1, w1,
     2                       w2,   w3)
c-----------------------------------------------------------------------
c one direction at the time
c-----------------------------------------------------------------------
c This routine calculates the constant part of the first-order density matrix :
c
c D1 = SUM(jocc,avir){Cj(T)*[F(D0,g1)-ej*S]*Ca/(ej-ea-xlv)*[Cj*Ca(T)+Ca*Cj(T)]}
c
c  It is calculated here as
c
c D1(p,q)=SUM(j,a)[C(T)*F-Eocc*C(T)S]*Cvirt](j,a)/(ej-ea-xlv)*[Cpj*Cqa+Cpa*Cqj]
c
c  (T) = transpose
c-----------------------------------------------------------------------
c  INTENT(IN)
c  fact    = occupancy : 2 for RHF, 1 for UHF
c  ncf     = number of contracted basis functions
c  nocc    = number of occupied orbitals
c  xlv     = level shift  (should not be used here, set it to zero)
c  fock    = first-order Fock matrix ,F(D0,g1) above, in triangular form
c  smat    = first-order AO overlap matrix in triangular form
c  cvec    = MO coefficients, C, the first nocc columns are occupied
c  eval    = orbital energies e above
c  INTENT(OUT)
c  d1      = the above part of the first-order density
c  STORAGE
c  w1      = a working array ncf**2 long
c  w2      = a working array ncf*nvirt long
c  w3      = a working array ncf*nocc long
c---------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.0d0,two=2.d0)
      dimension fock(*),cvec(ncf,ncf),eval(*),d1(*)
      dimension w1(ncf,ncf),w2(*),w3(nocc,ncf)
c
      nvirt=ncf-nocc
c      ivirtst=nocc*ncf+1
c  expand the Fock matrix to quadratic
      call quad(fock,w1,one,ncf)
c  W2=Cocc(T)*F
      call dgemm('t','n',nocc,ncf,ncf,
     1           one, cvec, ncf, w1, ncf,
     2           zero, w2, nocc)
c  expand the overlap matrix to quadratic
      call quad(smat,w1,one,ncf)
c  W3=Cocc(T)*S
      call dgemm('t','n',nocc,ncf,ncf,
     1           one, cvec, ncf, w1, ncf,
     2           zero, w3, nocc)
c  Multiply the rows of W3 by the occupied orbital energies and subtract
c  them from W2
      call RowMultiply(nocc,ncf,eval,w3,w2)
c  Build W1=[Cocc(T)*F-Eocc*Cocc(T)*S]Cvirt=W2*Cvirt
      call dgemm('n','n', nocc, nvirt,ncf,
     1            one, w2, nocc, cvec(1,nocc+1), ncf,
     2            zero,w1, nocc)
c  Now scale W1=Cocc(T)FCvirt with 1.0/(e(i)-e(a)-xlv)
      call scalebydenom(nocc,nvirt,eval,eval(nocc+1),xlv,w1)
c  Calculate Cvirt*W1(T). Multiply with the factor (occupancy)
c  Result is the perturbed coeff mx      
      call dgemm('n','t', ncf, nocc, nvirt,
     1            one, cvec(1,nocc+1), ncf, w1, nocc,
     2            zero,w2, ncf)
c Calculate W2*Cocc(T). Factor moved to this call by GM
      call dgemm('n','t',ncf,ncf,nocc,
     1            fact,w2,ncf,cvec,ncf,
     2            zero, w1,ncf)
c  Result, in quadratic form, is in W1.
c  Add transpose to W1 and transfer it to the triangular array d1
      call symtrian(ncf,w1,d1)
c      call printd1('new',ncf,d1)
      end
c======================================================================
      subroutine scalebydenom(nocc,nvirt,eocc,evirt,xlv,f)
c  This routine divides element (i,a) of the matrix F by
c  (eocc(i)-eocc(a)-xlv); F(i,a)=F(i,a)/(eocc(i)-eocc(a)-xlv)
c
c  Arguments:
c  INTENT(IN)
c  nocc     = number of occupied orbitals, number of rows of F
c  nvirt    = number of virtual orbitals, number of columns of F
c  eocc     = occupied orbital energies
c  evirt    = virtual orbital energies
c  INTENT(INOUT)
c  f        = Fock matrix (occupied x virtual part in MO basis)
      implicit real*8 (a-h,o-z)
      integer a
      dimension f(nocc,nvirt),eocc(nocc),evirt(nvirt)
      do a=1,nvirt
        xx=evirt(a)+xlv
        do i=1,nocc
          yy=eocc(i)-xx
          f(i,a)=f(i,a)/yy
        end do
      end do
      end
c======================================================================
      subroutine symtrian(n,a,b)
c  This routine adds A+A(T) to the symmetrical matrix B stored as
c  the upper triangle row-wise.
c  Arguments:
c  INTENT(IN)
c  n - dimension of the square matrix A
c  A(n,n) - square matrix
c  INTENT(OUT)
c  B(1:n*(n+1)/2): triangular matrix, it i,j element is A(i,j)+A(j,i)
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(*)
      ij=0
      do i=1,n
        do j=1,i
          ij=ij+1
          b(ij)=a(i,j)+a(j,i)
        end do
      end do
      end
c======================================================================
      subroutine RowMultiply(nocc,ncf,eocc,a,b)
c  B=B-Eocc*A
c  A,B= nocc x ncf matrices
c  Eocc is diagonal (orb. energies)
c  Arguments:
c  INTENT(IN)
c  nocc   = number of rows of A,B, Eocc  (number of occupied orbitals)
c  ncf    = number of columns of A,B (number of AOs)
c  Eocc   = occupied orbital energies (vector, diagonal matrix)
c  A      = nocc x ncf
c  INTENT(INOUT)
c  B      = nocc x ncf
      implicit real*8 (a-h,o-z)
      dimension eocc(nocc),a(nocc,ncf),b(nocc,ncf)
      do k=1,ncf
        do i=1,nocc
           b(i,k)=b(i,k)-eocc(i)*a(i,k)
        end do
      end do
      end
c======================================================================
      subroutine printd1(text,n,d1)
      implicit real*8 (a-h,o-z)
      dimension d1(*)
      character text*3
      write(6,*) text
      ib=0
      do i=1,n
        ie=ib+i
        write(6,100) (d1(k),k=ib+1,ie)
        write(6,*)
        ib=ie
      end do
  100 format((5f12.6))
      end
c======================================================================
c use only for hessian cphf
      subroutine d1const_ds1d(rhf,ncf,ntri,nocc,lind,
     *                       vec,val,d0  ,f1  ,s1  ,
     *                       d1c,bl,ds1d )
c-----------------------------------------------------------------
c Calculates the constant part of the D1 matrix for :
c
c  (1) electric field perturbation (dipole polarizability)
c  (2) nuclear displacement perturbation (analytical hessian)
c
c It does it for three directions at once (x,y,z)
c-----------------------------------------------------------------
c INPUT :
c
c rhf      - rhf or uhf flag
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c nocc     - number of occupied orbitals
c lind     - diagonals of i*(i-1)/2
c vec      - eigen vectors
c val      - eigen values
c d0       - unpreturbed density
c f1       - Ist-order fock matrix
c s1       - Ist-order overlap matrix
c bl       - storage for everything
c OUTPUT :
c
c d1c      - constant part of D1 :
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eo*S1)Cv*[CoCv+ + CvCo+]
c ds1d     - -1/2DS1D
c-----------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      dimension lind(*)
      dimension bl(*)
      dimension d0(*), vec(*), val(*)
      dimension s1(ntri,3), f1(ntri,3)
      dimension d1c(ntri,3)                      ! output
      dimension ds1d(ntri,3)                     ! output
c--------------------------------------------------------------
c    Calculate the constant part of D10
c   ----------------------------------------
c 1. contributions from [ F1(d0,g1) - ei S1 ]
c 2. contributions from the 0.5 D0*S1*D0
c
c  S1 is not zero only for field dependent basis set
c
c The total constant part of the D1 stored in d1con(*)
c
c--------------------------------------------------------------
cc      call getrval('xlvsh',xlvsh)      ! disabled    JB
cc - NOTE: If this is reenabled, it should be passed as an
cc   argument and NOT obtained from the depository - this
cc   routine is called on the slaves where depository
cc   stuff does NOT work    <-----   IMPORTANT
      xlvsh = 0.0d0
c--------------------------------------------------------------
      fact1=2.0d0     ! for dens1_part1 like occupancy number
      fact2=0.5d0     ! for dS1d
      if(.not.rhf) then
        fact1=1.0d0
        fact2=1.0d0
      endif
c--------------------------------------------------------------
c 1. calculate contributions from [ F1(d0,g1) - ei S1 ]
c
c one direction at the time :
c
      call getmem(ncf**2,lw1)
      nvirt=ncf-nocc
      memvirt = ncf*MAX(nvirt,nocc)
      call getmem(memvirt,lw2)
      call getmem(ncf*nocc,lw3)
c
      do icr=1,3
         call dens1_1dir2n(fact1,ncf, nocc, xlvsh, f1(1,icr),
     *                    s1(1,icr),vec, val, d1c(1,icr), bl(lw1),
     *                    bl(lw2),bl(lw3))
c      call printd1('new',ncf,d1c(1,icr))
      enddo
c
      call retmem(3)
c
c  2. calculate -fact2*D0*S1*D0 and add it to the constant part of D10.
c
      call zeroit(ds1d,3*ntri)
      call getmem(ncf*ncf,lw3)
      do icr=1,3
      call ds1d_1dir(fact2,bl(lw3),ncf,lind,d0,s1(1,icr),ds1d(1,icr))
      enddo
c
c make final d1const :
c
      call add2vec(3*ntri,d1c,ds1d,d1c)
c
      call retmem(1)
c
c---------------------------------------------------------------
cc      write(6,*)' D1const  exit   '
cc      call drumh(d1c(1,1),ncf, 6  ,'D1con-x ')
cc      call drumh(d1c(1,2),ncf, 6  ,'D1con-y ')
cc      call drumh(d1c(1,3),ncf, 6  ,'D1con-z ')
c--------------------------------------------------------------
      end
c======================================================================
