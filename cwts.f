c======================================================================
c
c  MM 2003
c
c  this file contains several numerical quadrature subroutines
c  for computing the dft contributions to the constant part of
c  derivative fock matrices and hessian, taking into account
c  the weight derivatives.
c
c  there are three sets of subroutines, for different cases
c  depending on the amount of memory available for the
c  calculation:
c
c  Case 1: all the derivative fock matrices can be computed
c          in one pass.
c
c          RHF case -->  cwopfh
c          UHF case -->  cwopfhu
c
c              These subroutines perform quadrature for all
c              the components of the fock derivative, plush
c              quadrature for the hessian, without applying
c              translational invarince. The translational
c              invariance can be applied just once for each
c              atom, at the end of the quadrature.
c
c  Case 2: multiple passes, component icntr has to be computed
c          in current pass.
c
c          RHF case --> cwmpfh
c          UHF case --> cwmpfhu
c
c              These subroutines perform quadrature for the
c              requested component of the fock derivative,
c              component icntr is computed applying translational
c              invariance at each grid pont. Quadrature for the
c              Hessian is also performed, without translational
c              invariance. Translational invariance for the
c              hessian can be applied separately, at the end of
c              the quadrature.
c
c  Case 3: multiple passes, component icntr is not computed
c          in current pass.
c
c          RHF case --> cwmpf
c          UHF case --> cwmpfu
c
c              These subroutines perform quadrature for the
c              requested components of the fock derivative only.
c              (use only if component icntr is not computed in
c               current pass). No need for translational invariance
c
c =====================================================================
      subroutine cwopfh(dft,    npp,    nbas,   nbf,    natoms,
     $                  nbatm,  thrsh,  da,     dm,     gden,
     $                  wght,   pra,    prara,  prarb,  pga,
     $                  pgc,    praga,  pragb,  pragc,  pgaga,
     $                  pgagb,  pgagc,  pgcgc,  vao,    vaox,
     $                  vaoxx,  vaoxxx, inb,    vm,     indx,
     $                  denx,   denxx,  gdx,    gdxx,    gx,
     $                  gxx,     sv,     sw,    icntr,   gwt,
     $                  hwt,    fda,    hess,   td1,     tg1,
     $                  tsw,    tqf,    tqh)
c     implicit real*8 (a-h,o-z)
      implicit none             ! trying to avoid typing errors
c
c
c  carries out numerical integration and accumulates contribution
c  from current grid point into derivative fock matrices and hessian
c  taking into account the weight derivatives.
c
c  All components of the derivative fock matrices are computed in
c  one pass. Note that the icntr components of derivative fock and
c  hessian are not computed by this subroutine, as they can be
c  obtained by translational invariance.
c
c  ** closed shell **
c
c  arguments
c
c  dft     -  method flag (note: all methods include slater exchange)
c              1 - 3 - local correlation
c             >3 - nonlocal
c  npp     -  number of contributing (non-zero) grid points this batch
c  nbas    -  total number of basis functions
c  nbf     -  indexing array to location of "non-zero" basis functions
c  natoms  -  number of atoms
c  nbatm   -  basis functions --> atom index
c  thrsh   -  threshold for neglect of contribution
c  da      -  closed-shell density matrix (lower triangle)
c  dm      -  maximum density matrix element per column
c  gden    -  density gradient at grid points (non local only, dft > 3)
c  wght    -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density (dft > 3)
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
c  vao     -  "non-zero" basis function values at grid points
c  potx    -  gradient functional derivative (dft > 3 only)
c  potxx   -  gradient functional second derivative (dft > 3)
c  potxd   -  density-gradient mixed second derivative (dft > 3)
c  vaox    -  basis function 1st derivatives at grid points
c  vaoxx   -  basis function 2nd derivatives at grid points
c  vaoxxx  -  basis function 3rd derivatives at grid points (dft > 3)
c  inb     -  indexing array for non-zero entries to vao
c  vm      -  array containing maximum magnitude ao per grid point
c  indx    -  index into contributing columns of vao
c  denx    -  scratch storage for density gradient per grid point
c  denxx   -    ditto for density hessian per grid point
c  gdx     -  ditto for atomic gradient of density gradient (dft > 3)
c  gdxx    -  ditto for atomic hessian of density gradient (dft > 3)
c  gx      -  ditto for atomic gradient of gradient invariant (dft > 3)
c  gxx     -  ditto for atomic hessian of gradient invariant (dft > 3)
c  sv      -  storage for coefficient of vao(i)*vao(j) in quadrature
c          -  formula (dft > 3)
c  sw      -  storage for coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  icntr   -  current atomic center
c  gwt     -  gradient of quadrature weight
c  hwt     -  hessian of quadrature weight
c
c  on exit
c
c  fda     -  contribution to derivative fock matrices
c  hess    -  contribution to hessian matrix
c
c
      integer npp,nbas,natoms,icntr
      real*8 wght(*),gwt(3,natoms,*),hwt(3,natoms,3,natoms,*)
      real*8 vao(*),vaox(3,*),vaoxx(6,*),vaoxxx(10,*)
      real*8 vm(*),da(nbas*(nbas+1)/2),dm(nbas),gden(3,*)
      real*8 pra(npp),pga(npp),pgc(npp),prara(npp)
      real*8 prarb(npp),praga(npp),pragb(npp),pragc(npp)
      real*8 pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
      real*8 denx(3,natoms),denxx(3,natoms,3,natoms)
      real*8 gdx(3,3,natoms),gdxx(3,3,natoms,3,natoms)
      real*8 gx(3,natoms),gxx(3,natoms,3,natoms)
      real*8 sv(3,natoms),sw(3,3,natoms)
      integer dft,nbf(*),inb(*),indx(npp),nbatm(nbas)
      real*8 fda(3,natoms,*),hess(3,natoms,3,natoms)
      real*8 thrsh
      real*8 td1,tg1,tsw,tqf,tqh
c
      real*8 zero,two,three,four,six,twelve,half
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0,half=0.5d0)
      real*8 epsi,alot
      parameter (epsi=1.0d-15,alot=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,doval
      integer nat3,ip,ipp,i,j,ii,jj,ij,it
      integer ia,ja,iatm,jatm,iixx,katm,k1,isv,isw
      real*8 ra,ra2,ga,gc,rara,rara2,raga,ragb,ragc,prap
      real*8 gaga,gagb,gagc,gcgc,vmx,vmx2,wg,prg,prg2,pgg,pgg2,pg,pg2
      real*8 abra,abrara,thrx,thrx1,thrxx,valt,abdaij,abdaijm,abval
      real*8 valmm,valmm1,valmm2,valxx,valxy,valyy,valxz,valyz,valzz
      real*8 abvj,abvvj,valj,valjx,valjy,valjz,abijt,abij
      real*8 thrsh1,thtest
      real*8 t1,t2,t3,t4,t5,t6
      real*8 gwtmax,dmx,valx,valy,valz,daij,xyc,xzc,yzc
      real*8 dmaxyz,vmax,val,vali,vdxc,valm,vdjx,xc
      real*8 hvalx,hvaly,hvalz,dx,dy,dz,sxx,sxy,sxz
      real*8 vd1,vd2,vd3,svmax,swmax,tswmax
      real*8 vmax1,vmax2,vmax3
      real*8 valix,valiy,valiz,xvi,xxvx,xxvy,xxvz,valtest,valmxj
      real*8 valmx,valmx1,xvj,valm1,val1,smaxmx,val2,val3,val4
      real*8 vij,vijx,vijy,vijz,gij,hdx,hdy,hdz,hgx,hgy,hgz
      real*8 daijt,smax,dmax,potp,potxp,potp2,potxp2
      real*8 sxxi,sxyi,sxzi,dpmax2
      real*8 valt1,valt2
c
c
      nat3 = 3*natoms
c
      if(dft.le.3) then
c
c  local only
c
      do 50 ip=1,npp
      ipp = indx(ip)
      wg  = wght(ipp)
      prap=pra(ip)
      ra = prap*wg
      rara = prara(ip)*wg
      vmx = vm(ipp)
c
c   initializations for weight derivatives
c
      gwtmax=zero
      do ia=1,natoms
        if(ia.ne.icntr)then
          gwt(1,ia,ipp)=gwt(1,ia,ipp)*prap
          gwt(2,ia,ipp)=gwt(2,ia,ipp)*prap
          gwt(3,ia,ipp)=gwt(3,ia,ipp)*prap
          gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +    abs(gwt(3,ia,ipp)))
        endif
      enddo
      call zeroit(denx,nat3)
      call zeroit(denxx,nat3**2)
c
c    form gradient and hessian density at current grid point.
c
c    note:
c    for the closed shell case, the array da contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx and denxx will contain half the gradient and hessian
c    of the total closed shell density
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrara=abs(rara)
      dodenx=abrara.gt.epsi
      thrx=alot
      if(dodenx)thrx=thrsh/abrara
      thrx1=thrx/vmx2
      abra=abs(ra)
      dodenxx=abra.gt.epsi
      thrxx=alot
      if(dodenxx)thrxx=thrsh/abra
      if(.not.(dodenx.or.dodenxx))goto 25 ! should never happen!
c
c   now compute density derivatives
c
      do 20 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        dmx = dm(ii)
        if(dmx*vmx**2.lt.thrsh) go to 20
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        if(iatm.eq.icntr)goto 20
        valx = vaox(1,i)
        valy = vaox(2,i)
        valz = vaox(3,i)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx*vmx
        dox=(dodenx.and.(valt.gt.thrx1.or.valt**2.gt.thrx))
        doxx=(dodenxx.and.dmx*vmx2.gt.thrxx)
        if(.not.(dox.or.doxx)) goto 20
        do 18 j=nbf(ipp)+1,i
          jj = inb(j)
          ij = it + jj
          jatm = nbatm(jj)
          daij = da(ij)*vao(j)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denx(1,iatm) = denx(1,iatm) - daij*valx
              denx(2,iatm) = denx(2,iatm) - daij*valy
              denx(3,iatm) = denx(3,iatm) - daij*valz
          endif
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
          if(doxx.and.abdaijm.gt.thrxx)then
              denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*vaoxx(1,i)
              xyc=daij*vaoxx(2,i)
              denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
              denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
              denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*vaoxx(3,i)
              xzc=daij*vaoxx(4,i)
              denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
              denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
              yzc=daij*vaoxx(5,i)
              denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
              denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
              denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*vaoxx(6,i)
          endif
          if(jatm.lt.iatm) go to 18
          if(jatm.eq.icntr)goto 18
          abdaijm=abs(da(ij))*valm*vmx
          if(doxx.and.abdaijm.gt.thrxx)then
c
c    (b) iatm with jatm
c
            daij = da(ij)*vaox(1,j)
              denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
              denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
              denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
              daij = da(ij)*vaox(2,j)
              denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
              denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
              denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
              daij = da(ij)*vaox(3,j)
              denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
              denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
              denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
          endif
 18     continue
        do 19 j=i+1,nbf(ipp+1)
          jj = inb(j)
          ij = (jj*(jj-1))/2 + ii
          jatm = nbatm(jj)
          daij = da(ij)*vao(j)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denx(1,iatm) = denx(1,iatm) - daij*valx
              denx(2,iatm) = denx(2,iatm) - daij*valy
              denx(3,iatm) = denx(3,iatm) - daij*valz
          endif
          if(doxx.and.abdaijm.gt.thrxx)then
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
              denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*vaoxx(1,i)
              xyc=daij*vaoxx(2,i)
              denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
              denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
              denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*vaoxx(3,i)
              xzc=daij*vaoxx(4,i)
              denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
              denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
              yzc=daij*vaoxx(5,i)
              denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
              denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
              denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*vaoxx(6,i)
          endif
c
c    (b) iatm with jatm
c
          if(jatm.lt.iatm) go to 19
          if(jatm.eq.icntr)goto 19
          abdaijm=abs(da(ij))*valm*vmx
          if(doxx.and.abdaijm.gt.thrxx)then
              daij = da(ij)*vaox(1,j)
              denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
              denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
              denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
              daij = da(ij)*vaox(2,j)
              denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
              denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
              denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
              daij = da(ij)*vaox(3,j)
              denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
              denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
              denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
          endif
 19     continue
 20   continue
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c    numerical quadrature for the derivative fock matrix.
c
c -- get maximum absolute value of 1st-order density at this grid point
c
      call absmax(nat3,denx,iixx,dmaxyz)
c
c -- global threshold testing
c
      vmax = max(abra*vmx2,abrara*dmaxyz*vmx2)
      if(vmax.lt.thrsh) go to 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      do 40 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        vali = vao(i)
        val = vali*ra
        vdxc = vali*rara
        valx = vaox(1,i)*ra
        valy = vaox(2,i)*ra
        valz = vaox(3,i)*ra
        abval= abs(val)
        doval= abval.gt.thrsh1
        valm = max(abs(valx),abs(valy),abs(valz))
        valt = max(abval,valm,abs(vdxc)*dmaxyz)
        if(valt.gt.thrsh1) then
          do 30 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            jatm = nbatm(jj)
            valj = vao(j)
            if(abs(valj)*valm.gt.thrsh)then
              if(iatm.ne.icntr) then
                fda(1,iatm,ij) = fda(1,iatm,ij) - valx*valj
                fda(2,iatm,ij) = fda(2,iatm,ij) - valy*valj
                fda(3,iatm,ij) = fda(3,iatm,ij) - valz*valj
              endif
            endif
            if(doval)then
              if(jatm.ne.icntr) then
                fda(1,jatm,ij) = fda(1,jatm,ij) - val*vaox(1,j)
                fda(2,jatm,ij) = fda(2,jatm,ij) - val*vaox(2,j)
                fda(3,jatm,ij) = fda(3,jatm,ij) - val*vaox(3,j)
              endif
            endif
c
c -- contribution of functional derivative
c
            vdjx = vdxc*valj
            if(abs(vdjx)*dmaxyz.gt.thrsh) then
              do katm=1,natoms
                fda(1,katm,ij)=fda(1,katm,ij) + vdjx*denx(1,katm)
                fda(2,katm,ij)=fda(2,katm,ij) + vdjx*denx(2,katm)
                fda(3,katm,ij)=fda(3,katm,ij) + vdjx*denx(3,katm)
              enddo
            endif
c
c -- weight derivatives contribution
c
            xc=vali*valj
            do ia=1,natoms
              if(ia.ne.icntr)then
                fda(1,ia,ij)=fda(1,ia,ij)+xc*gwt(1,ia,ipp)
                fda(2,ia,ij)=fda(2,ia,ij)+xc*gwt(2,ia,ipp)
                fda(3,ia,ij)=fda(3,ia,ij)+xc*gwt(3,ia,ipp)
              endif
            enddo
 30       continue
        endif
 40   continue
c
c    numerical quadrature for the hessian matrix
c
 45   continue
      call secund(t3)
      tqf=tqf+t3-t2
c
c -- direct contribution to hessian matrix.
c    a factor two is applied
c
      rara2=two*rara
      ra2=two*ra
      do iatm=1,natoms
        hvalx=rara2*denx(1,iatm)
        hvaly=rara2*denx(2,iatm)
        hvalz=rara2*denx(3,iatm)
        do jatm=iatm,natoms
          hess(1,iatm,1,jatm) = hess(1,iatm,1,jatm) +
     $     denx(1,jatm)*hvalx + ra2*denxx(1,iatm,1,jatm)
          hess(2,iatm,1,jatm) = hess(2,iatm,1,jatm) +
     $     denx(1,jatm)*hvaly + ra2*denxx(2,iatm,1,jatm)
          hess(3,iatm,1,jatm) = hess(3,iatm,1,jatm) +
     $     denx(1,jatm)*hvalz + ra2*denxx(3,iatm,1,jatm)
          hess(1,iatm,2,jatm) = hess(1,iatm,2,jatm) +
     $     denx(2,jatm)*hvalx + ra2*denxx(1,iatm,2,jatm)
          hess(2,iatm,2,jatm) = hess(2,iatm,2,jatm) +
     $     denx(2,jatm)*hvaly + ra2*denxx(2,iatm,2,jatm)
          hess(3,iatm,2,jatm) = hess(3,iatm,2,jatm) +
     $     denx(2,jatm)*hvalz + ra2*denxx(3,iatm,2,jatm)
          hess(1,iatm,3,jatm) = hess(1,iatm,3,jatm) +
     $     denx(3,jatm)*hvalx + ra2*denxx(1,iatm,3,jatm)
          hess(2,iatm,3,jatm) = hess(2,iatm,3,jatm) +
     $     denx(3,jatm)*hvaly + ra2*denxx(2,iatm,3,jatm)
          hess(3,iatm,3,jatm) = hess(3,iatm,3,jatm) +
     $     denx(3,jatm)*hvalz + ra2*denxx(3,iatm,3,jatm)
        enddo
      enddo
c
c  --  add weight derivatives contribution
c
      do ia=1,natoms
      if(ia.ne.icntr)then
        do ja=ia,natoms
        if(ja.ne.icntr)then
        hess(1,ia,1,ja)=hess(1,ia,1,ja) + two*(gwt(1,ia,ipp)*denx(1,ja)
     $   + gwt(1,ja,ipp)*denx(1,ia)) + hwt(1,ia,1,ja,ipp)
        hess(2,ia,1,ja)=hess(2,ia,1,ja) + two*(gwt(2,ia,ipp)*denx(1,ja)
     $   + gwt(1,ja,ipp)*denx(2,ia)) + hwt(2,ia,1,ja,ipp)
        hess(3,ia,1,ja)=hess(3,ia,1,ja) + two*(gwt(3,ia,ipp)*denx(1,ja)
     $   + gwt(1,ja,ipp)*denx(3,ia)) + hwt(3,ia,1,ja,ipp)
        hess(1,ia,2,ja)=hess(1,ia,2,ja) + two*(gwt(1,ia,ipp)*denx(2,ja)
     $   + gwt(2,ja,ipp)*denx(1,ia)) + hwt(1,ia,2,ja,ipp)
        hess(2,ia,2,ja)=hess(2,ia,2,ja) + two*(gwt(2,ia,ipp)*denx(2,ja)
     $   + gwt(2,ja,ipp)*denx(2,ia)) + hwt(2,ia,2,ja,ipp)
        hess(3,ia,2,ja)=hess(3,ia,2,ja) + two*(gwt(3,ia,ipp)*denx(2,ja)
     $   + gwt(2,ja,ipp)*denx(3,ia)) + hwt(3,ia,2,ja,ipp)
        hess(1,ia,3,ja)=hess(1,ia,3,ja) + two*(gwt(1,ia,ipp)*denx(3,ja)
     $   + gwt(3,ja,ipp)*denx(1,ia)) + hwt(1,ia,3,ja,ipp)
        hess(2,ia,3,ja)=hess(2,ia,3,ja) + two*(gwt(2,ia,ipp)*denx(3,ja)
     $   + gwt(3,ja,ipp)*denx(2,ia)) + hwt(2,ia,3,ja,ipp)
        hess(3,ia,3,ja)=hess(3,ia,3,ja) + two*(gwt(3,ia,ipp)*denx(3,ja)
     $   + gwt(3,ja,ipp)*denx(3,ia)) + hwt(3,ia,3,ja,ipp)
        endif
        enddo
      endif
      enddo
      call secund(t4)
      tqh=tqh+t4-t3
c
c   local dft ends here
c
 50   continue
cc
      else
cc
c
c   non-local dft
c
      do 200 ip=1,npp
        ipp = indx(ip)
        wg = wght(ipp)
        vmx = vm(ipp)
        ra = pra(ip)*wg
        ga = pga(ip)*wg
        gc = pgc(ip)*wg
        rara = (prara(ip)+prarb(ip))*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  some sums of the potentials that will be used later
c
        prg=raga+ragb+ragc
        pgg=gaga+gagb+two*gagc+half*gcgc
        pg=ga+half*gc
c
c  potential combinations for weight derivatives contributions
c
        potp = pra(ip)
        potxp = pga(ip)+half*pgc(ip)
        potp2 = potp+potp
        potxp2 = potxp+potxp
c
c  density gradient at current point
c
        dx=gden(1,ipp)
        dy=gden(2,ipp)
        dz=gden(3,ipp)
c
c  zero out derivatives of densities and gradients
c
        call zeroit(denx,nat3)
        call zeroit(denxx,nat3**2)
        call zeroit(gdx,3*nat3)
        call zeroit(gdxx,3*nat3**2)
        call zeroit(gx,nat3)
        call zeroit(gxx,nat3**2)
c
c   initializations for weight derivatives
c
c   unlike the local case above, we do not multiply
c   the weight gradient by the potential at this stage
c
        gwtmax=zero
        do ia=1,natoms
          if(ia.ne.icntr)then
            gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +        abs(gwt(3,ia,ipp)))
          endif
        enddo
c
c    form the atomic gradient and hessian  of density
c    and density gradient at current grid point
c
c    note:
c    for the closed shell case, the array da contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx and denxx will contain half the atomic gradient and hessian
c    of the total closed shell density. likevise, gdx and gdxx will
c    contain half the atomic gradient and hessian of the total closed
c    shell density gradient
c
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c    in addition, we just check whether the second derivatives
c    of density and density gradient have to be computed at all.
c
        call secund(t1)
        abra=abs(ra)
        dodenxx=abra.gt.epsi
        do 220 i=nbf(ipp)+1,nbf(ipp+1)
          ii = inb(i)
          dmx = dm(ii)       ! max. element of first order density
          if(dmx.gt.epsi)then
             thtest=thrsh/(vmx*dmx)
          else
            goto 220
          endif
          it = (ii*(ii-1))/2
          iatm = nbatm(ii)
          if(iatm.eq.icntr) goto 220
          valx = vaox(1,i)
          valy = vaox(2,i)
          valz = vaox(3,i)
          valm = max(abs(valx),abs(valy),abs(valz))
          valxx = vaoxx(1,i)
          valxy = vaoxx(2,i)
          valyy = vaoxx(3,i)
          valxz = vaoxx(4,i)
          valyz = vaoxx(5,i)
          valzz = vaoxx(6,i)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          do 218 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            valj=vao(j)
            jatm = nbatm(jj)
            abvj=abs(valj)
            valjx=vaox(1,j)
            valjy=vaox(2,j)
            valjz=vaox(3,j)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            daijt=da(ij)
            abijt=abs(daijt)
            daij = daijt*valj
            abij=abs(daij)
c
c -- atomic gradient of density
c
              if(abij*valm.gt.thrsh)then
                denx(1,iatm) = denx(1,iatm) - daij*valx
                denx(2,iatm) = denx(2,iatm) - daij*valy
                denx(3,iatm) = denx(3,iatm) - daij*valz
              endif
c
c -- atomic gradient of density gradient
c
             if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              gdx(1,1,iatm)=gdx(1,1,iatm)-daijt*(valj*valxx+valx*valjx)
              gdx(1,2,iatm)=gdx(1,2,iatm)-daijt*(valj*valxy+valy*valjx)
              gdx(1,3,iatm)=gdx(1,3,iatm)-daijt*(valj*valxz+valz*valjx)
              gdx(2,1,iatm)=gdx(2,1,iatm)-daijt*(valj*valxy+valx*valjy)
              gdx(2,2,iatm)=gdx(2,2,iatm)-daijt*(valj*valyy+valy*valjy)
              gdx(2,3,iatm)=gdx(2,3,iatm)-daijt*(valj*valyz+valz*valjy)
              gdx(3,1,iatm)=gdx(3,1,iatm)-daijt*(valj*valxz+valx*valjz)
              gdx(3,2,iatm)=gdx(3,2,iatm)-daijt*(valj*valyz+valy*valjz)
              gdx(3,3,iatm)=gdx(3,3,iatm)-daijt*(valj*valzz+valz*valjz)
             endif
c
c -- (a) one center terms: iatm with iatm
c
c
c -- atomic hessian of density
c
              if(dodenxx.and.abij*valmm.gt.thrsh)then
                denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*valxx
                xyc=daij*valxy
                denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
                denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
                denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*valyy
                xzc=daij*valxz
                denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
                denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
                yzc=daij*valyz
                denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
                denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
                denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*valzz
              endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
                gdxx(1,1,iatm,1,iatm)=gdxx(1,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(1,i)*valj+valxx*valjx)
                gdxx(1,1,iatm,2,iatm)=gdxx(1,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,1,iatm)=gdxx(1,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,2,iatm)=gdxx(1,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valyy*valjx)
                gdxx(1,1,iatm,3,iatm)=gdxx(1,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,3,iatm,1,iatm)=gdxx(1,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,2,iatm,3,iatm)=gdxx(1,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,2,iatm)=gdxx(1,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,3,iatm)=gdxx(1,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valzz*valjx)
c
                gdxx(2,1,iatm,1,iatm)=gdxx(2,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxx*valjy)
                gdxx(2,1,iatm,2,iatm)=gdxx(2,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,1,iatm)=gdxx(2,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,2,iatm)=gdxx(2,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(4,i)*valj+valyy*valjy)
                gdxx(2,1,iatm,3,iatm)=gdxx(2,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,3,iatm,1,iatm)=gdxx(2,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,2,iatm,3,iatm)=gdxx(2,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,2,iatm)=gdxx(2,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,3,iatm)=gdxx(2,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valzz*valjy)

                gdxx(3,1,iatm,1,iatm)=gdxx(3,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxx*valjz)
                gdxx(3,1,iatm,2,iatm)=gdxx(3,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,1,iatm)=gdxx(3,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,2,iatm)=gdxx(3,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyy*valjz)
                gdxx(3,1,iatm,3,iatm)=gdxx(3,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,3,iatm,1,iatm)=gdxx(3,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,2,iatm,3,iatm)=gdxx(3,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,2,iatm)=gdxx(3,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,3,iatm)=gdxx(3,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(10,i)*valj+valzz*valjz)
              endif
c
c -- (b) two center terms: iatm with jatm
            if(jatm.lt.iatm) go to 218
            if(jatm.eq.icntr)goto 218
c
c -- atomic hessian of density
c
              if(dodenxx.and.abijt*abvvj*valm.gt.thrsh)then
                daij = daijt*valjx
                denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
                denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
                denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
                daij = daijt*valjy
                denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
                denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
                denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
                daij = daijt*valjz
                denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
                denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
                denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
              endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
                gdxx(1,1,iatm,1,jatm)=gdxx(1,1,iatm,1,jatm)+
     $                      daijt * (valxx*valjx+valx*vaoxx(1,j))
                gdxx(1,1,iatm,2,jatm)=gdxx(1,1,iatm,2,jatm)+
     $                      daijt * (valxx*valjy+valx*vaoxx(2,j))
                gdxx(1,2,iatm,1,jatm)=gdxx(1,2,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valy*vaoxx(1,j))
                gdxx(1,2,iatm,2,jatm)=gdxx(1,2,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valy*vaoxx(2,j))
                gdxx(1,1,iatm,3,jatm)=gdxx(1,1,iatm,3,jatm)+
     $                      daijt * (valxx*valjz+valx*vaoxx(4,j))
                gdxx(1,3,iatm,1,jatm)=gdxx(1,3,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valz*vaoxx(1,j))
                gdxx(1,2,iatm,3,jatm)=gdxx(1,2,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valy*vaoxx(4,j))
                gdxx(1,3,iatm,2,jatm)=gdxx(1,3,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valz*vaoxx(2,j))
                gdxx(1,3,iatm,3,jatm)=gdxx(1,3,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valz*vaoxx(4,j))
c
                gdxx(2,1,iatm,1,jatm)=gdxx(2,1,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valx*vaoxx(2,j))
                gdxx(2,1,iatm,2,jatm)=gdxx(2,1,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valx*vaoxx(3,j))
                gdxx(2,2,iatm,1,jatm)=gdxx(2,2,iatm,1,jatm)+
     $                      daijt * (valyy*valjx+valy*vaoxx(2,j))
                gdxx(2,2,iatm,2,jatm)=gdxx(2,2,iatm,2,jatm)+
     $                      daijt * (valyy*valjy+valy*vaoxx(3,j))
                gdxx(2,1,iatm,3,jatm)=gdxx(2,1,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valx*vaoxx(5,j))
                gdxx(2,3,iatm,1,jatm)=gdxx(2,3,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valz*vaoxx(2,j))
                gdxx(2,2,iatm,3,jatm)=gdxx(2,2,iatm,3,jatm)+
     $                      daijt * (valyy*valjz+valy*vaoxx(5,j))
                gdxx(2,3,iatm,2,jatm)=gdxx(2,3,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valz*vaoxx(3,j))
                gdxx(2,3,iatm,3,jatm)=gdxx(2,3,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valz*vaoxx(5,j))
c
                gdxx(3,1,iatm,1,jatm)=gdxx(3,1,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valx*vaoxx(4,j))
                gdxx(3,1,iatm,2,jatm)=gdxx(3,1,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valx*vaoxx(5,j))
                gdxx(3,2,iatm,1,jatm)=gdxx(3,2,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valy*vaoxx(4,j))
                gdxx(3,2,iatm,2,jatm)=gdxx(3,2,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valy*vaoxx(5,j))
                gdxx(3,1,iatm,3,jatm)=gdxx(3,1,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valx*vaoxx(6,j))
                gdxx(3,3,iatm,1,jatm)=gdxx(3,3,iatm,1,jatm)+
     $                      daijt * (valzz*valjx+valz*vaoxx(4,j))
                gdxx(3,2,iatm,3,jatm)=gdxx(3,2,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valy*vaoxx(6,j))
                gdxx(3,3,iatm,2,jatm)=gdxx(3,3,iatm,2,jatm)+
     $                      daijt * (valzz*valjy+valz*vaoxx(5,j))
                gdxx(3,3,iatm,3,jatm)=gdxx(3,3,iatm,3,jatm)+
     $                      daijt * (valzz*valjz+valz*vaoxx(6,j))
              endif
 218      continue
          do 219 j=i+1,nbf(ipp+1)
            jj = inb(j)
            ij = (jj*(jj-1))/2 + ii
            jatm = nbatm(jj)
            valj=vao(j)
            abvj=abs(valj)
            valjx=vaox(1,j)
            valjy=vaox(2,j)
            valjz=vaox(3,j)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            daijt= da(ij)
            abijt=abs(daijt)
            daij = daijt*valj
            abij=abs(daij)
c
c -- atomic gradient of density
c
              if(abij*valm.gt.thrsh)then
                denx(1,iatm) = denx(1,iatm) - daij*valx
                denx(2,iatm) = denx(2,iatm) - daij*valy
                denx(3,iatm) = denx(3,iatm) - daij*valz
              endif
c
c -- atomic gradient of density gradient
c
             if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              gdx(1,1,iatm)=gdx(1,1,iatm)-daijt*(valj*valxx+valx*valjx)
              gdx(1,2,iatm)=gdx(1,2,iatm)-daijt*(valj*valxy+valy*valjx)
              gdx(1,3,iatm)=gdx(1,3,iatm)-daijt*(valj*valxz+valz*valjx)
              gdx(2,1,iatm)=gdx(2,1,iatm)-daijt*(valj*valxy+valx*valjy)
              gdx(2,2,iatm)=gdx(2,2,iatm)-daijt*(valj*valyy+valy*valjy)
              gdx(2,3,iatm)=gdx(2,3,iatm)-daijt*(valj*valyz+valz*valjy)
              gdx(3,1,iatm)=gdx(3,1,iatm)-daijt*(valj*valxz+valx*valjz)
              gdx(3,2,iatm)=gdx(3,2,iatm)-daijt*(valj*valyz+valy*valjz)
              gdx(3,3,iatm)=gdx(3,3,iatm)-daijt*(valj*valzz+valz*valjz)
             endif
c
c -- (a) one center terms: iatm with iatm
c
c
c -- atomic hessian of density
c
            if(dodenxx.and.abij*valmm.gt.thrsh)then
              denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*valxx
              xyc=daij*valxy
              denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
              denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
              denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*valyy
              xzc=daij*valxz
              denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
              denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
              yzc=daij*valyz
              denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
              denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
              denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*valzz
            endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
                gdxx(1,1,iatm,1,iatm)=gdxx(1,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(1,i)*valj+valxx*valjx)
                gdxx(1,1,iatm,2,iatm)=gdxx(1,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,1,iatm)=gdxx(1,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,2,iatm)=gdxx(1,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valyy*valjx)
                gdxx(1,1,iatm,3,iatm)=gdxx(1,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,3,iatm,1,iatm)=gdxx(1,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,2,iatm,3,iatm)=gdxx(1,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,2,iatm)=gdxx(1,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,3,iatm)=gdxx(1,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valzz*valjx)
c
                gdxx(2,1,iatm,1,iatm)=gdxx(2,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxx*valjy)
                gdxx(2,1,iatm,2,iatm)=gdxx(2,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,1,iatm)=gdxx(2,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,2,iatm)=gdxx(2,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(4,i)*valj+valyy*valjy)
                gdxx(2,1,iatm,3,iatm)=gdxx(2,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,3,iatm,1,iatm)=gdxx(2,3,iatm,1,iatm)+
     $                 daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,2,iatm,3,iatm)=gdxx(2,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,2,iatm)=gdxx(2,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,3,iatm)=gdxx(2,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valzz*valjy)
c
                gdxx(3,1,iatm,1,iatm)=gdxx(3,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxx*valjz)
                gdxx(3,1,iatm,2,iatm)=gdxx(3,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,1,iatm)=gdxx(3,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,2,iatm)=gdxx(3,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyy*valjz)
                gdxx(3,1,iatm,3,iatm)=gdxx(3,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,3,iatm,1,iatm)=gdxx(3,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,2,iatm,3,iatm)=gdxx(3,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,2,iatm)=gdxx(3,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,3,iatm)=gdxx(3,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(10,i)*valj+valzz*valjz)
              endif
c
c -- (b) two center terms: iatm with jatm
c
            if(jatm.lt.iatm) go to 219
            if(jatm.eq.icntr)goto 219
c
c -- atomic hessian of density
c
              if(dodenxx.and.abijt*abvvj*valm.gt.thrsh)then
                daij = daijt*valjx
                denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
                denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
                denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
                daij = daijt*valjy
                denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
                denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
                denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
                daij = daijt*valjz
                denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
                denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
                denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
              endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
                gdxx(1,1,iatm,1,jatm)=gdxx(1,1,iatm,1,jatm)+
     $                      daijt * (valxx*valjx+valx*vaoxx(1,j))
                gdxx(1,1,iatm,2,jatm)=gdxx(1,1,iatm,2,jatm)+
     $                      daijt * (valxx*valjy+valx*vaoxx(2,j))
                gdxx(1,2,iatm,1,jatm)=gdxx(1,2,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valy*vaoxx(1,j))
                gdxx(1,2,iatm,2,jatm)=gdxx(1,2,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valy*vaoxx(2,j))
                gdxx(1,1,iatm,3,jatm)=gdxx(1,1,iatm,3,jatm)+
     $                      daijt * (valxx*valjz+valx*vaoxx(4,j))
                gdxx(1,3,iatm,1,jatm)=gdxx(1,3,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valz*vaoxx(1,j))
                gdxx(1,2,iatm,3,jatm)=gdxx(1,2,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valy*vaoxx(4,j))
                gdxx(1,3,iatm,2,jatm)=gdxx(1,3,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valz*vaoxx(2,j))
                gdxx(1,3,iatm,3,jatm)=gdxx(1,3,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valz*vaoxx(4,j))
c
                gdxx(2,1,iatm,1,jatm)=gdxx(2,1,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valx*vaoxx(2,j))
                gdxx(2,1,iatm,2,jatm)=gdxx(2,1,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valx*vaoxx(3,j))
                gdxx(2,2,iatm,1,jatm)=gdxx(2,2,iatm,1,jatm)+
     $                      daijt * (valyy*valjx+valy*vaoxx(2,j))
                gdxx(2,2,iatm,2,jatm)=gdxx(2,2,iatm,2,jatm)+
     $                      daijt * (valyy*valjy+valy*vaoxx(3,j))
                gdxx(2,1,iatm,3,jatm)=gdxx(2,1,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valx*vaoxx(5,j))
                gdxx(2,3,iatm,1,jatm)=gdxx(2,3,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valz*vaoxx(2,j))
                gdxx(2,2,iatm,3,jatm)=gdxx(2,2,iatm,3,jatm)+
     $                      daijt * (valyy*valjz+valy*vaoxx(5,j))
                gdxx(2,3,iatm,2,jatm)=gdxx(2,3,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valz*vaoxx(3,j))
                gdxx(2,3,iatm,3,jatm)=gdxx(2,3,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valz*vaoxx(5,j))
c
                gdxx(3,1,iatm,1,jatm)=gdxx(3,1,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valx*vaoxx(4,j))
                gdxx(3,1,iatm,2,jatm)=gdxx(3,1,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valx*vaoxx(5,j))
                gdxx(3,2,iatm,1,jatm)=gdxx(3,2,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valy*vaoxx(4,j))
                gdxx(3,2,iatm,2,jatm)=gdxx(3,2,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valy*vaoxx(5,j))
                gdxx(3,1,iatm,3,jatm)=gdxx(3,1,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valx*vaoxx(6,j))
                gdxx(3,3,iatm,1,jatm)=gdxx(3,3,iatm,1,jatm)+
     $                      daijt * (valzz*valjx+valz*vaoxx(4,j))
                gdxx(3,2,iatm,3,jatm)=gdxx(3,2,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valy*vaoxx(6,j))
                gdxx(3,3,iatm,2,jatm)=gdxx(3,3,iatm,2,jatm)+
     $                      daijt * (valzz*valjy+valz*vaoxx(5,j))
                gdxx(3,3,iatm,3,jatm)=gdxx(3,3,iatm,3,jatm)+
     $                      daijt * (valzz*valjz+valz*vaoxx(6,j))
              endif
 219      continue
 220    continue
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and hessian of density gradient invariant
c
        do iatm=1,natoms
          gx(1,iatm)=two*(dx*gdx(1,1,iatm)+
     $                   dy*gdx(2,1,iatm)+dz*gdx(3,1,iatm))
          gx(2,iatm)=two*(dx*gdx(1,2,iatm)+
     $                   dy*gdx(2,2,iatm)+dz*gdx(3,2,iatm))
          gx(3,iatm)=two*(dx*gdx(1,3,iatm)+
     $                   dy*gdx(2,3,iatm)+dz*gdx(3,3,iatm))
        enddo
c
        do iatm=1,natoms
          do jatm=iatm,natoms
            gxx(1,iatm,1,jatm)=two*(
     $        dx*gdxx(1,1,iatm,1,jatm)+gdx(1,1,iatm)*gdx(1,1,jatm)+
     $        dy*gdxx(2,1,iatm,1,jatm)+gdx(2,1,iatm)*gdx(2,1,jatm)+
     $        dz*gdxx(3,1,iatm,1,jatm)+gdx(3,1,iatm)*gdx(3,1,jatm))
            gxx(1,iatm,2,jatm)=two*(
     $        dx*gdxx(1,1,iatm,2,jatm)+gdx(1,1,iatm)*gdx(1,2,jatm)+
     $        dy*gdxx(2,1,iatm,2,jatm)+gdx(2,1,iatm)*gdx(2,2,jatm)+
     $        dz*gdxx(3,1,iatm,2,jatm)+gdx(3,1,iatm)*gdx(3,2,jatm))
            gxx(2,iatm,1,jatm)=two*(
     $        dx*gdxx(1,2,iatm,1,jatm)+gdx(1,2,iatm)*gdx(1,1,jatm)+
     $        dy*gdxx(2,2,iatm,1,jatm)+gdx(2,2,iatm)*gdx(2,1,jatm)+
     $        dz*gdxx(3,2,iatm,1,jatm)+gdx(3,2,iatm)*gdx(3,1,jatm))
            gxx(2,iatm,2,jatm)=two*(
     $        dx*gdxx(1,2,iatm,2,jatm)+gdx(1,2,iatm)*gdx(1,2,jatm)+
     $        dy*gdxx(2,2,iatm,2,jatm)+gdx(2,2,iatm)*gdx(2,2,jatm)+
     $        dz*gdxx(3,2,iatm,2,jatm)+gdx(3,2,iatm)*gdx(3,2,jatm))
            gxx(1,iatm,3,jatm)=two*(
     $        dx*gdxx(1,1,iatm,3,jatm)+gdx(1,1,iatm)*gdx(1,3,jatm)+
     $        dy*gdxx(2,1,iatm,3,jatm)+gdx(2,1,iatm)*gdx(2,3,jatm)+
     $        dz*gdxx(3,1,iatm,3,jatm)+gdx(3,1,iatm)*gdx(3,3,jatm))
            gxx(3,iatm,1,jatm)=two*(
     $        dx*gdxx(1,3,iatm,1,jatm)+gdx(1,3,iatm)*gdx(1,1,jatm)+
     $        dy*gdxx(2,3,iatm,1,jatm)+gdx(2,3,iatm)*gdx(2,1,jatm)+
     $        dz*gdxx(3,3,iatm,1,jatm)+gdx(3,3,iatm)*gdx(3,1,jatm))
            gxx(2,iatm,3,jatm)=two*(
     $        dx*gdxx(1,2,iatm,3,jatm)+gdx(1,2,iatm)*gdx(1,3,jatm)+
     $        dy*gdxx(2,2,iatm,3,jatm)+gdx(2,2,iatm)*gdx(2,3,jatm)+
     $        dz*gdxx(3,2,iatm,3,jatm)+gdx(3,2,iatm)*gdx(3,3,jatm))
            gxx(3,iatm,2,jatm)=two*(
     $        dx*gdxx(1,3,iatm,2,jatm)+gdx(1,3,iatm)*gdx(1,2,jatm)+
     $        dy*gdxx(2,3,iatm,2,jatm)+gdx(2,3,iatm)*gdx(2,2,jatm)+
     $        dz*gdxx(3,3,iatm,2,jatm)+gdx(3,3,iatm)*gdx(3,2,jatm))
            gxx(3,iatm,3,jatm)=two*(
     $        dx*gdxx(1,3,iatm,3,jatm)+gdx(1,3,iatm)*gdx(1,3,jatm)+
     $        dy*gdxx(2,3,iatm,3,jatm)+gdx(2,3,iatm)*gdx(2,3,jatm)+
     $        dz*gdxx(3,3,iatm,3,jatm)+gdx(3,3,iatm)*gdx(3,3,jatm))
          enddo
        enddo
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the fock matrix
c    (i.e., quantities v, w and x at page 7436 of johnson and frisch).
c
        prg2=two*prg
        pgg2=two*pgg
        pg2=two*pg
        sxx=pg2*dx
        sxy=pg2*dy
        sxz=pg2*dz
        do iatm=1,natoms
          vd1=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          vd2=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          vd3=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          sv(1,iatm)=rara*denx(1,iatm)+prg*gx(1,iatm)
          sv(2,iatm)=rara*denx(2,iatm)+prg*gx(2,iatm)
          sv(3,iatm)=rara*denx(3,iatm)+prg*gx(3,iatm)
          sw(1,1,iatm)=pg2*gdx(1,1,iatm)+dx*vd1
          sw(2,1,iatm)=pg2*gdx(2,1,iatm)+dy*vd1
          sw(3,1,iatm)=pg2*gdx(3,1,iatm)+dz*vd1
          sw(1,2,iatm)=pg2*gdx(1,2,iatm)+dx*vd2
          sw(2,2,iatm)=pg2*gdx(2,2,iatm)+dy*vd2
          sw(3,2,iatm)=pg2*gdx(3,2,iatm)+dz*vd2
          sw(1,3,iatm)=pg2*gdx(1,3,iatm)+dx*vd3
          sw(2,3,iatm)=pg2*gdx(2,3,iatm)+dy*vd3
          sw(3,3,iatm)=pg2*gdx(3,3,iatm)+dz*vd3
        enddo
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    fock matrix.
c
c -- get the maximum absolute value of coefficients v and w
c
        call absmax(nat3,sv,isv,svmax)
        call absmax(3*nat3,sw,isw,swmax)
        tswmax=three*swmax
c -- global threshold testing
        vmx2=vmx*vmx
        smax = abs(sxx+sxy+sxz)
        dmax = abs(dx+dy+dz)
        vmax3 = (svmax+tswmax+tswmax)*vmx2
        vmax1 = (abra+two*smax)*vmx2
        dpmax2 = potxp2*dmax*(vmx2+vmx2)
        vmax2 = gwtmax*(potp*vmx2+dpmax2)
        vmax=max(vmax1,vmax2,vmax3)
        if(vmax.lt.thrsh)  go to 245
c
c -- numerical quadrature  for derivative fock matrix
c
        do 240 i=nbf(ipp)+1,nbf(ipp+1)
          vali = vao(i)
          valix = vaox(1,i)
          valiy = vaox(2,i)
          valiz = vaox(3,i)
          xvi=sxx*valix+sxy*valiy+sxz*valiz+ra*vali
          sxxi=sxx*vali
          sxyi=sxy*vali
          sxzi=sxz*vali
          xxvx=sxx*vaoxx(1,i)+sxy*vaoxx(2,i)+sxz*vaoxx(4,i)+ra*valix
          xxvy=sxx*vaoxx(2,i)+sxy*vaoxx(3,i)+sxz*vaoxx(5,i)+ra*valiy
          xxvz=sxx*vaoxx(4,i)+sxy*vaoxx(5,i)+sxz*vaoxx(6,i)+ra*valiz
          valm=abs(vali)
          valm1=abs(xvi)
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1=max(abs(xxvx),abs(xxvy),abs(xxvz))
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          val4=gwtmax*(potp*valm*vmx+dpmax2)
          valt1=max(val1,val3)
          valt2=max(val3,val4)
          valtest=max(valt1,valt2)
          if(valtest.gt.thrsh)then
            ii = inb(i)
            it = (ii*(ii-1))/2
            iatm = nbatm(ii)
            do 230 j=nbf(ipp)+1,i
              jj = inb(j)
              ij = it + jj
              jatm = nbatm(jj)
              valj = vao(j)
              valjx = vaox(1,j)
              valjy = vaox(2,j)
              valjz = vaox(3,j)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              xvj=sxx*valjx+sxy*valjy+sxz*valjz
              if(valmx1*abs(valj)+valmx*abs(xvj).gt.thrsh)then
              if(iatm.ne.icntr) then
                fda(1,iatm,ij) = fda(1,iatm,ij) - xxvx*valj - valix*xvj
                fda(2,iatm,ij) = fda(2,iatm,ij) - xxvy*valj - valiy*xvj
                fda(3,iatm,ij) = fda(3,iatm,ij) - xxvz*valj - valiz*xvj
              endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
              if(jatm.ne.icntr) then
                fda(1,jatm,ij) = fda(1,jatm,ij) - xvi*valjx -
     $          vaoxx(1,j)*sxxi - vaoxx(2,j)*sxyi - vaoxx(4,j)*sxzi
                fda(2,jatm,ij) = fda(2,jatm,ij) - xvi*valjy -
     $          vaoxx(2,j)*sxxi - vaoxx(3,j)*sxyi - vaoxx(5,j)*sxzi
                fda(3,jatm,ij) = fda(3,jatm,ij) - xvi*valjz -
     $          vaoxx(4,j)*sxxi - vaoxx(5,j)*sxyi - vaoxx(6,j)*sxzi
              endif
              endif
              vij=vali*valj
              vijx=valix*valj+valjx*vali
              vijy=valiy*valj+valjy*vali
              vijz=valiz*valj+valjz*vali
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).gt.thrsh)then
                do katm=1,natoms
                  fda(1,katm,ij) = fda(1,katm,ij) + (sv(1,katm)*vij+
     $            sw(1,1,katm)*vijx+sw(2,1,katm)*vijy+sw(3,1,katm)*vijz)
                  fda(2,katm,ij) = fda(2,katm,ij) + (sv(2,katm)*vij+
     $            sw(1,2,katm)*vijx+sw(2,2,katm)*vijy+sw(3,2,katm)*vijz)
                  fda(3,katm,ij) = fda(3,katm,ij) + (sv(3,katm)*vij+
     $            sw(1,3,katm)*vijx+sw(2,3,katm)*vijy+sw(3,3,katm)*vijz)
                enddo
              endif
c
c -- weight derivatives contribution
c
              gij=potp*vij+potxp2*(dx*vijx+dy*vijy+dz*vijz)
              if(gwtmax*abs(gij).gt.thrsh)then
                do ia=1,natoms
                  if(ia.ne.icntr)then
                    fda(1,ia,ij)=fda(1,ia,ij)+gwt(1,ia,ipp)*gij
                    fda(2,ia,ij)=fda(2,ia,ij)+gwt(2,ia,ipp)*gij
                    fda(3,ia,ij)=fda(3,ia,ij)+gwt(3,ia,ipp)*gij
                  endif
                enddo
              endif
 230        continue
          endif
 240    continue
c
c   numerical quadrature for direct contribution to hessian matrix
c
 245    continue
        call secund(t5)
        tqf=tqf+t5-t4
c
c -- direct contribution to hessian matrix.
c    a factor two is applied
c
        ra2=two*ra
        rara2=two*rara
        do iatm=1,natoms
          hdx=rara2*denx(1,iatm)+prg2*gx(1,iatm)
          hdy=rara2*denx(2,iatm)+prg2*gx(2,iatm)
          hdz=rara2*denx(3,iatm)+prg2*gx(3,iatm)
          hgx=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          hgy=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          hgz=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          do jatm=iatm,natoms
       hess(1,iatm,1,jatm)=hess(1,iatm,1,jatm)+hdx*denx(1,jatm)+
     $ hgx*gx(1,jatm)+ra2*denxx(1,iatm,1,jatm)+pg2*gxx(1,iatm,1,jatm)
       hess(2,iatm,1,jatm)=hess(2,iatm,1,jatm)+hdy*denx(1,jatm)+
     $ hgy*gx(1,jatm)+ra2*denxx(2,iatm,1,jatm)+pg2*gxx(2,iatm,1,jatm)
       hess(3,iatm,1,jatm)=hess(3,iatm,1,jatm)+hdz*denx(1,jatm)+
     $ hgz*gx(1,jatm)+ra2*denxx(3,iatm,1,jatm)+pg2*gxx(3,iatm,1,jatm)
       hess(1,iatm,2,jatm)=hess(1,iatm,2,jatm)+hdx*denx(2,jatm)+
     $ hgx*gx(2,jatm)+ra2*denxx(1,iatm,2,jatm)+pg2*gxx(1,iatm,2,jatm)
       hess(2,iatm,2,jatm)=hess(2,iatm,2,jatm)+hdy*denx(2,jatm)+
     $ hgy*gx(2,jatm)+ra2*denxx(2,iatm,2,jatm)+pg2*gxx(2,iatm,2,jatm)
       hess(3,iatm,2,jatm)=hess(3,iatm,2,jatm)+hdz*denx(2,jatm)+
     $ hgz*gx(2,jatm)+ra2*denxx(3,iatm,2,jatm)+pg2*gxx(3,iatm,2,jatm)
       hess(1,iatm,3,jatm)=hess(1,iatm,3,jatm)+hdx*denx(3,jatm)+
     $ hgx*gx(3,jatm)+ra2*denxx(1,iatm,3,jatm)+pg2*gxx(1,iatm,3,jatm)
       hess(2,iatm,3,jatm)=hess(2,iatm,3,jatm)+hdy*denx(3,jatm)+
     $ hgy*gx(3,jatm)+ra2*denxx(2,iatm,3,jatm)+pg2*gxx(2,iatm,3,jatm)
       hess(3,iatm,3,jatm)=hess(3,iatm,3,jatm)+hdz*denx(3,jatm)+
     $ hgz*gx(3,jatm)+ra2*denxx(3,iatm,3,jatm)+pg2*gxx(3,iatm,3,jatm)
          enddo
        enddo
c
c  --  add weight derivatives contribution.
c      the sv array is used to store a partial sum
c
        do ia=1,natoms
          sv(1,ia)=potp2*denx(1,ia)+potxp2*gx(1,ia)
          sv(2,ia)=potp2*denx(2,ia)+potxp2*gx(2,ia)
          sv(3,ia)=potp2*denx(3,ia)+potxp2*gx(3,ia)
        enddo
        do ia=1,natoms
        if(ia.ne.icntr)then
          do ja=ia,natoms
          if(ja.ne.icntr)then
            hess(1,ia,1,ja)=hess(1,ia,1,ja)+gwt(1,ia,ipp)*sv(1,ja)
     $               +gwt(1,ja,ipp)*sv(1,ia)+hwt(1,ia,1,ja,ipp)
            hess(2,ia,1,ja)=hess(2,ia,1,ja)+gwt(2,ia,ipp)*sv(1,ja)
     $               +gwt(1,ja,ipp)*sv(2,ia)+hwt(2,ia,1,ja,ipp)
            hess(3,ia,1,ja)=hess(3,ia,1,ja)+gwt(3,ia,ipp)*sv(1,ja)
     $               +gwt(1,ja,ipp)*sv(3,ia)+hwt(3,ia,1,ja,ipp)
            hess(1,ia,2,ja)=hess(1,ia,2,ja)+gwt(1,ia,ipp)*sv(2,ja)
     $               +gwt(2,ja,ipp)*sv(1,ia)+hwt(1,ia,2,ja,ipp)
            hess(2,ia,2,ja)=hess(2,ia,2,ja)+gwt(2,ia,ipp)*sv(2,ja)
     $               +gwt(2,ja,ipp)*sv(2,ia)+hwt(2,ia,2,ja,ipp)
            hess(3,ia,2,ja)=hess(3,ia,2,ja)+gwt(3,ia,ipp)*sv(2,ja)
     $               +gwt(2,ja,ipp)*sv(3,ia)+hwt(3,ia,2,ja,ipp)
            hess(1,ia,3,ja)=hess(1,ia,3,ja)+gwt(1,ia,ipp)*sv(3,ja)
     $               +gwt(3,ja,ipp)*sv(1,ia)+hwt(1,ia,3,ja,ipp)
            hess(2,ia,3,ja)=hess(2,ia,3,ja)+gwt(2,ia,ipp)*sv(3,ja)
     $               +gwt(3,ja,ipp)*sv(2,ia)+hwt(2,ia,3,ja,ipp)
            hess(3,ia,3,ja)=hess(3,ia,3,ja)+gwt(3,ia,ipp)*sv(3,ja)
     $               +gwt(3,ja,ipp)*sv(3,ia)+hwt(3,ia,3,ja,ipp)
          endif
          enddo
        endif
        enddo
        call secund(t6)
        tqh=tqh+t6-t5
 200  continue
      endif
c
      return
      end
c =====================================================================
      subroutine cwopfhu(dft,    npp,    nbas,   nbf,    natoms,
     $                   nbatm,  thrsh,  da,     db,     dm,
     $                   gdena,  gdenb,  wght,   pra,    prb,
     $                   prara,  prbrb,  prarb,  pga,    pgb,
     $                   pgc,    praga,  pragb,  pragc,  prbga,
     $                   prbgb,  prbgc,  pgaga,  pgagb,  pgagc,
     $                   pgbgb,  pgbgc,  pgcgc,  vao,    vaox,
     $                   vaoxx,  vaoxxx, inb,    vm,     indx,
     $                   denxa,  denxb,  denxxa, denxxb, gdxa,
     $                   gdxb,   gdxxa,  gdxxb,  gxaa,   gxbb,
     $                   gxab,   gxxaa,  gxxbb,  gxxab,  sva,
     $                   svb,    swa,    swb,    sga,    sgb,
     $                   sgc,    icntr,  gwt,    hwt,    fda,
     $                   fdb,    hess,   td1,    tg1,    tsw,
     $                   tqf,    tqh)
c     implicit real*8 (a-h,o-z)
      implicit none             ! trying to avoid typing errors
c
c
c  carries out numerical integration and accumulates contribution
c  from current grid point into derivative fock matrices and hessian
c  taking into account the weight derivatives.
c
c  All components of the derivative fock matrices are computed in
c  one pass. Note that the icntr components of derivative fock and
c  hessian are not computed by this subroutine, as they can be
c  obtained by translational invariance.
c
c  ** open shell **
c
c  arguments
c
c  dft     -  method flag (note: all methods include slater exchange)
c              1 - 3 - local correlation
c             >3 - nonlocal
c  npp     -  number of contributing (non-zero) grid points this batch
c  nbas    -  total number of basis functions
c  nbf     -  indexing array to location of "non-zero" basis functions
c  natoms  -  number of atoms
c  nbatm   -  basis functions --> atom index
c  thrsh   -  threshold for neglect of contribution
c  da      -  alpha density matrix (lower triangle)
c  db      -  beta density matrix (lower triangle)
c  dm      -  maximum density matrix element per column
c  gdena   -  alpha density gradient at grid points (only dft > 3)
c  gdena   -  beta density gradient at grid points (dft > 3)
c  wght    -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prb   -  Functional derivative w.r.t. beta density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgb   -  Funct. deriv. w.r.t. beta gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  prbga -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha grad.
C           (dft > 3)
C  prbgb -  Funct. 2nd. deriv. w.r.t. beta dens. and  grad. (dft > 3)
C  prbgc -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgbgb -  Funct. 2nd. deriv. w.r.t. beta grad. (dft > 3)
C  pgbgc -  Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
c  vao     -  "non-zero" basis function values at grid points
c  vaox    -  basis function 1st derivatives at grid points
c  vaoxx   -  basis function 2nd derivatives at grid points
c  vaoxxx  -  basis function 3rd derivatives at grid points (dft > 3)
c  inb     -  indexing array for non-zero entries to vao
c  vm      -  array containing maximum magnitude ao per grid point
c  indx    -  index into contributing columns of vao
c  denxa   -  scratch storage for alpah density gradient per grid point
c  denxb   -  scratch storage for beta density gradient per grid point
c  denxxa  -  ditto for alpha density hessian per grid point
c  denxxb  -  ditto for beta density hessian per grid point
c  gdxa    -  ditto for alpha atomic gradient of density gradient
c             (dft > 3)
c  gdxb    -  ditto for beta atomic gradient of density gradient
c             (dft > 3)
c  gdxxa   -  ditto for atomic hessian of alpha density gradient
c             (dft > 3)
c  gdxxb   -  ditto for atomic hessian of beta density gradient
c             (dft > 3)
c  gxaa    -  ditto for atomic gradient of alpha gradient invariant
c             (dft > 3)
c  gxbb    -  ditto for atomic gradient of beta gradient invariant
c             (dft > 3)
c  gxbb    -  ditto for atomic gradient of alpha beta gradient
c             invariant (dft > 3)
c  gxxaa   -  ditto for atomic hessian of alpha gradient invariant
c             (dft > 3)
c  gxxbb   -  ditto for atomic hessian of beta gradient invariant
c             (dft > 3)
c  gxxbb   -  ditto for atomic hessian of alpha beta gradient invariant
c             (dft > 3)
c  sva     -  storage for alpha coefficient of vao(i)*vao(j) in
c             quadrature formula
c  svb     -  storage for beta coefficient of vao(i)*vao(j) in
c             quadrature formula
c  swa     -  storage for alpha coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  swb     -  storage for beta coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  sga     -  storage for partial sum for hessian quadrature (dft > 3)
c  sgb     -  storage for partial sum for hessian quadrature (dft > 3)
c  sgc     -  storage for partial sum for hessian quadrature (dft > 3)
c  icntr   -  current atomic center
c  gwt     -  gradient of quadrature weight
c  hwt     -  hessian of quadrature weight
c
c  on exit
c
c  fda     -  contribution to alpha derivative fock matrices
c  fdb     -  contribution to beta derivative fock matrices
c  hess    -  contribution to hessian matrix
c
c
c
c
      integer npp,nbas,natoms,icntr
      real*8 wght(*),gwt(3,natoms,*),hwt(3,natoms,3,natoms,*)
      real*8 vao(*),vaox(3,*),vaoxx(6,*),vaoxxx(10,*)
      real*8 vm(*),da(nbas*(nbas+1)/2),dm(nbas),gdena(3,*)
      real*8 db(nbas*(nbas+1)/2),gdenb(3,*)
      real*8 pra(npp),pga(npp),pgc(npp),prara(npp)
      real*8 prb(npp),pgb(npp),prbrb(npp)
      real*8 prarb(npp),praga(npp),pragb(npp),pragc(npp)
      real*8 prbga(npp),prbgb(npp),prbgc(npp)
      real*8 pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
      real*8 pgbgb(npp),pgbgc(npp)
      real*8 denxa(3,natoms),denxxa(3,natoms,3,natoms)
      real*8 denxb(3,natoms),denxxb(3,natoms,3,natoms)
      real*8 gdxa(3,3,natoms),gdxxa(3,3,natoms,3,natoms)
      real*8 gdxb(3,3,natoms),gdxxb(3,3,natoms,3,natoms)
      real*8 gxaa(3,natoms),gxxaa(3,natoms,3,natoms)
      real*8 gxbb(3,natoms),gxxbb(3,natoms,3,natoms)
      real*8 gxab(3,natoms),gxxab(3,natoms,3,natoms)
      real*8 sva(3,natoms),svb(3,natoms)
      real*8 swa(3,3,natoms),swb(3,3,natoms)
      real*8 sga(3,natoms),sgb(3,natoms),sgc(3,natoms)
      integer dft,nbf(*),inb(*),indx(npp),nbatm(nbas)
      real*8 fda(3,natoms,*),hess(3,natoms,3,natoms)
      real*8 fdb(3,natoms,*)
      real*8 thrsh
      real*8 td1,tg1,tsw,tqf,tqh
c
      real*8 zero,two,three,four,six,twelve,half
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0,half=0.5d0)
      real*8 epsi,alot
      parameter (epsi=1.0d-15,alot=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,doval
      integer nat3,ip,ipp,i,j,ii,jj,ij,it
      integer ia,ja,iatm,jatm,katm,iixx,k1
      real*8 t1,t2,t3,t4,t5,t6
      real*8 wg,vmx,vmx2,prap,prbp,gwtmax,pgap,pgbp,pgcp
      real*8 ra,rb,rara,rbrb,rarb,abrr,abr
      real*8 ga,gb,gc,raga,ragb,ragc,rbga,rbgb,rbgc
      real*8 gaga,gagb,gagc,gbgb,gbgc,gcgc
      real*8 dax,day,daz,dbx,dby,dbz
      real*8 dax2,day2,daz2,dbx2,dby2,dbz2
      real*8 ga2,gb2,pgap2,pgbp2
      real*8 thrx,thrxx,thrx1,thrsh1
      real*8 dmx2,daij,dbij,daij2,dbij2,abdij,abdijm
      real*8 valx,valy,valz,valm,valt,abval
      real*8 vali,vala,valxa,valya,valza
      real*8 valj,valb,valxb,valyb,valzb
      real*8 valij,valija,valijb
      real*8 xyc,xzc,yzc
      real*8 svamax,svbmax,svmax,vmax,abgmax
      real*8 swamax,swbmax,swmax,tswmax,vmax1,vmax3
      real*8 valxx,valxy,valxz,valyy,valyz,valzz
      real*8 valmm,valm1,valmm1,valmm2
      real*8 valmx,valmx1a,valmx1b,valmx1,val1,val2,val3,valtest
      real*8 thtest,valmxj,xvja,xvjb,abxvj
      real*8 valjx,valjy,valjz,abvj,abvvj,abijt,abij
      real*8 valix,valiy,valiz,xvia,xvib
      real*8 sxxa,sxya,sxza,sxxb,sxyb,sxzb
      real*8 sxxia,sxyia,sxzia,sxxib,sxyib,sxzib,smax,smaxmx
      real*8 xxvxa,xxvya,xxvza,xxvxb,xxvyb,xxvzb
      real*8 vij,vijx,vijy,vijz
      real*8 dxija,dxijb,gija,gijb
c
c
      nat3 = 3*natoms
c
      if(dft.le.3) then
c
c  local only
c
      do 50 ip=1,npp
      ipp = indx(ip)
      wg  = wght(ipp)
      prap=pra(ip)
      prbp=prb(ip)
      ra = prap*wg
      rb = prbp*wg
      rara = prara(ip)*wg
      rbrb = prbrb(ip)*wg
      rarb = prarb(ip)*wg
      vmx = vm(ipp)
c
c   initializations for weight derivatives
c
      gwtmax=zero
      do ia=1,natoms
        if(ia.ne.icntr)then
          gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +    abs(gwt(3,ia,ipp)))
        endif
      enddo
      call zeroit(denxa,nat3)
      call zeroit(denxb,nat3)
      call zeroit(denxxa,nat3**2)
      call zeroit(denxxb,nat3**2)
c
c    form gradient and hessian density at current grid point.
c
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrr=max(abs(rara),abs(rbrb))+abs(rarb)
      dodenx=abrr.gt.epsi
      thrx=alot
      if(dodenx)thrx=thrsh/abrr
      thrx1=thrx/vmx2
      abr=max(abs(ra),abs(rb))
      dodenxx=abr.gt.epsi
      thrxx=alot
      if(dodenxx)thrxx=thrsh/abr
      if(.not.(dodenx.or.dodenxx))goto 25 ! should never happen!
c
c   now compute density derivatives
c
      do 20 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        dmx2 = dm(ii)
        if(dmx2*vmx2.lt.thrsh) go to 20
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        if(iatm.eq.icntr)goto 20
        valx = vaox(1,i)
        valy = vaox(2,i)
        valz = vaox(3,i)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx2*vmx
        dox=(dodenx.and.(two*valt.gt.thrx1.or.valt**2.gt.thrx))
        doxx=(dodenxx.and.dmx2*vmx2.gt.thrxx)
        if(.not.(dox.or.doxx)) goto 20
        do 18 j=nbf(ipp)+1,i
          jj = inb(j)
          ij = it + jj
          jatm = nbatm(jj)
          daij2=two*da(ij)
          dbij2=two*db(ij)
          daij = daij2*vao(j)
          dbij = dbij2*vao(j)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denxa(1,iatm) = denxa(1,iatm) - daij*valx
              denxa(2,iatm) = denxa(2,iatm) - daij*valy
              denxa(3,iatm) = denxa(3,iatm) - daij*valz
c
              denxb(1,iatm) = denxb(1,iatm) - dbij*valx
              denxb(2,iatm) = denxb(2,iatm) - dbij*valy
              denxb(3,iatm) = denxb(3,iatm) - dbij*valz
          endif
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
          if(doxx.and.abdijm.gt.thrxx)then
             denxxa(1,iatm,1,iatm)=denxxa(1,iatm,1,iatm)+daij*vaoxx(1,i)
             xyc=daij*vaoxx(2,i)
             denxxa(1,iatm,2,iatm)=denxxa(1,iatm,2,iatm)+xyc
             denxxa(2,iatm,1,iatm)=denxxa(2,iatm,1,iatm)+xyc
             denxxa(2,iatm,2,iatm)=denxxa(2,iatm,2,iatm)+daij*vaoxx(3,i)
             xzc=daij*vaoxx(4,i)
             denxxa(1,iatm,3,iatm)=denxxa(1,iatm,3,iatm)+xzc
             denxxa(3,iatm,1,iatm)=denxxa(3,iatm,1,iatm)+xzc
             yzc=daij*vaoxx(5,i)
             denxxa(2,iatm,3,iatm)=denxxa(2,iatm,3,iatm)+yzc
             denxxa(3,iatm,2,iatm)=denxxa(3,iatm,2,iatm)+yzc
             denxxa(3,iatm,3,iatm)=denxxa(3,iatm,3,iatm)+daij*vaoxx(6,i)
c
             denxxb(1,iatm,1,iatm)=denxxb(1,iatm,1,iatm)+dbij*vaoxx(1,i)
             xyc=dbij*vaoxx(2,i)
             denxxb(1,iatm,2,iatm)=denxxb(1,iatm,2,iatm)+xyc
             denxxb(2,iatm,1,iatm)=denxxb(2,iatm,1,iatm)+xyc
             denxxb(2,iatm,2,iatm)=denxxb(2,iatm,2,iatm)+dbij*vaoxx(3,i)
             xzc=dbij*vaoxx(4,i)
             denxxb(1,iatm,3,iatm)=denxxb(1,iatm,3,iatm)+xzc
             denxxb(3,iatm,1,iatm)=denxxb(3,iatm,1,iatm)+xzc
             yzc=dbij*vaoxx(5,i)
             denxxb(2,iatm,3,iatm)=denxxb(2,iatm,3,iatm)+yzc
             denxxb(3,iatm,2,iatm)=denxxb(3,iatm,2,iatm)+yzc
             denxxb(3,iatm,3,iatm)=denxxb(3,iatm,3,iatm)+dbij*vaoxx(6,i)
          endif
          if(jatm.lt.iatm) go to 18
          if(jatm.eq.icntr) go to 18
          abdijm=max(abs(da(ij)),abs(db(ij)))*valm*vmx
          if(doxx.and.abdijm.gt.thrxx)then
c
c    (b) iatm with jatm
c
              daij = daij2*vaox(1,j)
              denxxa(1,iatm,1,jatm)=denxxa(1,iatm,1,jatm)+daij*valx
              denxxa(2,iatm,1,jatm)=denxxa(2,iatm,1,jatm)+daij*valy
              denxxa(3,iatm,1,jatm)=denxxa(3,iatm,1,jatm)+daij*valz
              daij = daij2*vaox(2,j)
              denxxa(1,iatm,2,jatm)=denxxa(1,iatm,2,jatm)+daij*valx
              denxxa(2,iatm,2,jatm)=denxxa(2,iatm,2,jatm)+daij*valy
              denxxa(3,iatm,2,jatm)=denxxa(3,iatm,2,jatm)+daij*valz
              daij = daij2*vaox(3,j)
              denxxa(1,iatm,3,jatm)=denxxa(1,iatm,3,jatm)+daij*valx
              denxxa(2,iatm,3,jatm)=denxxa(2,iatm,3,jatm)+daij*valy
              denxxa(3,iatm,3,jatm)=denxxa(3,iatm,3,jatm)+daij*valz
c
              dbij = dbij2*vaox(1,j)
              denxxb(1,iatm,1,jatm)=denxxb(1,iatm,1,jatm)+dbij*valx
              denxxb(2,iatm,1,jatm)=denxxb(2,iatm,1,jatm)+dbij*valy
              denxxb(3,iatm,1,jatm)=denxxb(3,iatm,1,jatm)+dbij*valz
              dbij = dbij2*vaox(2,j)
              denxxb(1,iatm,2,jatm)=denxxb(1,iatm,2,jatm)+dbij*valx
              denxxb(2,iatm,2,jatm)=denxxb(2,iatm,2,jatm)+dbij*valy
              denxxb(3,iatm,2,jatm)=denxxb(3,iatm,2,jatm)+dbij*valz
              dbij = dbij2*vaox(3,j)
              denxxb(1,iatm,3,jatm)=denxxb(1,iatm,3,jatm)+dbij*valx
              denxxb(2,iatm,3,jatm)=denxxb(2,iatm,3,jatm)+dbij*valy
              denxxb(3,iatm,3,jatm)=denxxb(3,iatm,3,jatm)+dbij*valz
          endif
 18     continue
        do 19 j=i+1,nbf(ipp+1)
          jj = inb(j)
          ij = (jj*(jj-1))/2 + ii
          jatm = nbatm(jj)
          daij2 = two*da(ij)
          dbij2 = two*db(ij)
          daij = daij2*vao(j)
          dbij = dbij2*vao(j)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denxa(1,iatm) = denxa(1,iatm) - daij*valx
              denxa(2,iatm) = denxa(2,iatm) - daij*valy
              denxa(3,iatm) = denxa(3,iatm) - daij*valz
c
              denxb(1,iatm) = denxb(1,iatm) - dbij*valx
              denxb(2,iatm) = denxb(2,iatm) - dbij*valy
              denxb(3,iatm) = denxb(3,iatm) - dbij*valz
          endif
          if(doxx.and.abdijm.gt.thrxx)then
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
             denxxa(1,iatm,1,iatm)=denxxa(1,iatm,1,iatm)+daij*vaoxx(1,i)
             xyc=daij*vaoxx(2,i)
             denxxa(1,iatm,2,iatm)=denxxa(1,iatm,2,iatm)+xyc
             denxxa(2,iatm,1,iatm)=denxxa(2,iatm,1,iatm)+xyc
             denxxa(2,iatm,2,iatm)=denxxa(2,iatm,2,iatm)+daij*vaoxx(3,i)
             xzc=daij*vaoxx(4,i)
             denxxa(1,iatm,3,iatm)=denxxa(1,iatm,3,iatm)+xzc
             denxxa(3,iatm,1,iatm)=denxxa(3,iatm,1,iatm)+xzc
             yzc=daij*vaoxx(5,i)
             denxxa(2,iatm,3,iatm)=denxxa(2,iatm,3,iatm)+yzc
             denxxa(3,iatm,2,iatm)=denxxa(3,iatm,2,iatm)+yzc
             denxxa(3,iatm,3,iatm)=denxxa(3,iatm,3,iatm)+daij*vaoxx(6,i)
c
             denxxb(1,iatm,1,iatm)=denxxb(1,iatm,1,iatm)+dbij*vaoxx(1,i)
             xyc=dbij*vaoxx(2,i)
             denxxb(1,iatm,2,iatm)=denxxb(1,iatm,2,iatm)+xyc
             denxxb(2,iatm,1,iatm)=denxxb(2,iatm,1,iatm)+xyc
             denxxb(2,iatm,2,iatm)=denxxb(2,iatm,2,iatm)+dbij*vaoxx(3,i)
             xzc=dbij*vaoxx(4,i)
             denxxb(1,iatm,3,iatm)=denxxb(1,iatm,3,iatm)+xzc
             denxxb(3,iatm,1,iatm)=denxxb(3,iatm,1,iatm)+xzc
             yzc=dbij*vaoxx(5,i)
             denxxb(2,iatm,3,iatm)=denxxb(2,iatm,3,iatm)+yzc
             denxxb(3,iatm,2,iatm)=denxxb(3,iatm,2,iatm)+yzc
             denxxb(3,iatm,3,iatm)=denxxb(3,iatm,3,iatm)+dbij*vaoxx(6,i)
          endif
c
c    (b) iatm with jatm
c
          if(jatm.lt.iatm) go to 19
          if(jatm.eq.icntr) go to 19
          abdijm=max(abs(da(ij)),abs(db(ij)))*valm*vmx
          if(doxx.and.abdijm.gt.thrxx)then
              daij = daij2*vaox(1,j)
              denxxa(1,iatm,1,jatm)=denxxa(1,iatm,1,jatm)+daij*valx
              denxxa(2,iatm,1,jatm)=denxxa(2,iatm,1,jatm)+daij*valy
              denxxa(3,iatm,1,jatm)=denxxa(3,iatm,1,jatm)+daij*valz
              daij = daij2*vaox(2,j)
              denxxa(1,iatm,2,jatm)=denxxa(1,iatm,2,jatm)+daij*valx
              denxxa(2,iatm,2,jatm)=denxxa(2,iatm,2,jatm)+daij*valy
              denxxa(3,iatm,2,jatm)=denxxa(3,iatm,2,jatm)+daij*valz
              daij = daij2*vaox(3,j)
              denxxa(1,iatm,3,jatm)=denxxa(1,iatm,3,jatm)+daij*valx
              denxxa(2,iatm,3,jatm)=denxxa(2,iatm,3,jatm)+daij*valy
              denxxa(3,iatm,3,jatm)=denxxa(3,iatm,3,jatm)+daij*valz
c
              dbij = dbij2*vaox(1,j)
              denxxb(1,iatm,1,jatm)=denxxb(1,iatm,1,jatm)+dbij*valx
              denxxb(2,iatm,1,jatm)=denxxb(2,iatm,1,jatm)+dbij*valy
              denxxb(3,iatm,1,jatm)=denxxb(3,iatm,1,jatm)+dbij*valz
              dbij = dbij2*vaox(2,j)
              denxxb(1,iatm,2,jatm)=denxxb(1,iatm,2,jatm)+dbij*valx
              denxxb(2,iatm,2,jatm)=denxxb(2,iatm,2,jatm)+dbij*valy
              denxxb(3,iatm,2,jatm)=denxxb(3,iatm,2,jatm)+dbij*valz
              dbij = dbij2*vaox(3,j)
              denxxb(1,iatm,3,jatm)=denxxb(1,iatm,3,jatm)+dbij*valx
              denxxb(2,iatm,3,jatm)=denxxb(2,iatm,3,jatm)+dbij*valy
              denxxb(3,iatm,3,jatm)=denxxb(3,iatm,3,jatm)+dbij*valz
          endif
 19     continue
 20   continue
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c -- now form the coefficients that multiply the basis functions
c    in the expression for the Fock Matrix
c    (i.e., quantities V at page 7436 of Johnson and Frisch).
c
        DO IAtm=1,NAtoms
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)
          svb(1,iatm)=rbrb*denxb(1,iatm)+rarb*denxa(1,iatm)
          svb(2,iatm)=rbrb*denxb(2,iatm)+rarb*denxa(2,iatm)
          svb(3,iatm)=rbrb*denxb(3,iatm)+rarb*denxa(3,iatm)
        EndDO
c
c    numerical quadrature for the derivative fock matrix.
c
c -- get maximum absolute value of sva and svb
c
      Call absmax(NAt3,sva,iixx,svamax)
      Call absmax(NAt3,svb,iixx,svbmax)
      svmax=max(svamax,svbmax)
c
c -- global threshold testing
c
      VMax = Max(Abr*VMx2,Abrr*svmax*VMx2)
      if(vmax.lt.thrsh) go to 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      do 40 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        vali = vao(i)
        vala = vali*ra
        valb = vali*rb
        valxa = vaox(1,i)*ra
        valya = vaox(2,i)*ra
        valza = vaox(3,i)*ra
        valxb = vaox(1,i)*rb
        valyb = vaox(2,i)*rb
        valzb = vaox(3,i)*rb
        abval= max(abs(vala),abs(valb))
        doval= abval.gt.thrsh1
        valm = max(abs(vaox(1,i)),abs(vaox(2,i)),abs(vaox(3,i)))*abr
        Valt = MAX(abval,valm,Abs(vali)*svmax)
        if(valt.gt.thrsh1) then
          do 30 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            jatm = nbatm(jj)
            valj = vao(j)
            if(abs(valj)*valm.gt.thrsh)then
              if(iatm.ne.icntr) then
                fda(1,iatm,ij) = fda(1,iatm,ij) - valxa*valj
                fda(2,iatm,ij) = fda(2,iatm,ij) - valya*valj
                fda(3,iatm,ij) = fda(3,iatm,ij) - valza*valj
                fdb(1,iatm,ij) = fdb(1,iatm,ij) - valxb*valj
                fdb(2,iatm,ij) = fdb(2,iatm,ij) - valyb*valj
                fdb(3,iatm,ij) = fdb(3,iatm,ij) - valzb*valj
              endif
            endif
            if(doval)then
              if(jatm.ne.icntr) then
                fda(1,jatm,ij) = fda(1,jatm,ij) - vala*vaox(1,j)
                fda(2,jatm,ij) = fda(2,jatm,ij) - vala*vaox(2,j)
                fda(3,jatm,ij) = fda(3,jatm,ij) - vala*vaox(3,j)
                fdb(1,jatm,ij) = fdb(1,jatm,ij) - valb*vaox(1,j)
                fdb(2,jatm,ij) = fdb(2,jatm,ij) - valb*vaox(2,j)
                fdb(3,jatm,ij) = fdb(3,jatm,ij) - valb*vaox(3,j)
              endif
            endif
c
c -- contribution of functional derivative
c
            valij = vali*valj
            if(abs(valij)*svmax.gt.thrsh) then
              do katm=1,natoms
                fda(1,katm,ij)=fda(1,katm,ij) + valij*sva(1,katm)
                fda(2,katm,ij)=fda(2,katm,ij) + valij*sva(2,katm)
                fda(3,katm,ij)=fda(3,katm,ij) + valij*sva(3,katm)
                fdb(1,katm,ij)=fdb(1,katm,ij) + valij*svb(1,katm)
                fdb(2,katm,ij)=fdb(2,katm,ij) + valij*svb(2,katm)
                fdb(3,katm,ij)=fdb(3,katm,ij) + valij*svb(3,katm)
              enddo
            endif
c
c -- weight derivatives contribution
c
            valija=valij*prap
            valijb=valij*prbp
            do ia=1,natoms
              if(ia.ne.icntr)then
                fda(1,ia,ij)=fda(1,ia,ij)+valija*gwt(1,ia,ipp)
                fda(2,ia,ij)=fda(2,ia,ij)+valija*gwt(2,ia,ipp)
                fda(3,ia,ij)=fda(3,ia,ij)+valija*gwt(3,ia,ipp)
                fdb(1,ia,ij)=fdb(1,ia,ij)+valijb*gwt(1,ia,ipp)
                fdb(2,ia,ij)=fdb(2,ia,ij)+valijb*gwt(2,ia,ipp)
                fdb(3,ia,ij)=fdb(3,ia,ij)+valijb*gwt(3,ia,ipp)
              endif
            enddo
 30       continue
        endif
 40   continue
c
c    numerical quadrature for the hessian matrix
c
 45   continue
      call secund(t3)
      tqf=tqf+t3-t2
c
c -- direct contribution to hessian matrix.
c
      do iatm=1,natoms
        do jatm=iatm,natoms
          hess(1,iatm,1,jatm) = hess(1,iatm,1,jatm) +
     $     sva(1,iatm)*denxa(1,jatm) + svb(1,iatm)*denxb(1,jatm)+
     $     ra*denxxa(1,iatm,1,jatm) + rb*denxxb(1,iatm,1,jatm)
          hess(2,iatm,1,jatm) = hess(2,iatm,1,jatm) +
     $     sva(2,iatm)*denxa(1,jatm) + svb(2,iatm)*denxb(1,jatm)+
     $     ra*denxxa(2,iatm,1,jatm) + rb*denxxb(2,iatm,1,jatm)
          hess(3,iatm,1,jatm) = hess(3,iatm,1,jatm) +
     $     sva(3,iatm)*denxa(1,jatm) + svb(3,iatm)*denxb(1,jatm)+
     $     ra*denxxa(3,iatm,1,jatm) + rb*denxxb(3,iatm,1,jatm)
          hess(1,iatm,2,jatm) = hess(1,iatm,2,jatm) +
     $     sva(1,iatm)*denxa(2,jatm) + svb(1,iatm)*denxb(2,jatm)+
     $     ra*denxxa(1,iatm,2,jatm) + rb*denxxb(1,iatm,2,jatm)
          hess(2,iatm,2,jatm) = hess(2,iatm,2,jAtm) +
     $     sva(2,iatm)*denxa(2,jatm) + svb(2,iatm)*denxb(2,jatm)+
     $     ra*denxxa(2,iatm,2,jatm) + rb*denxxb(2,iatm,2,jatm)
          hess(3,iatm,2,jatm) = hess(3,iatm,2,jatm) +
     $     sva(3,iatm)*denxa(2,jatm) + svb(3,iatm)*denxb(2,jatm)+
     $     ra*denxxa(3,iatm,2,jatm) + rb*denxxb(3,iatm,2,jatm)
          hess(1,iatm,3,jatm) = hess(1,iatm,3,jatm) +
     $     sva(1,iatm)*denxa(3,jatm) + svb(1,iatm)*denxb(3,jatm)+
     $     ra*denxxa(1,iatm,3,jatm) + rb*denxxb(1,iatm,3,jatm)
          hess(2,iatm,3,jatm) = hess(2,iatm,3,jatm) +
     $     sva(2,iatm)*denxa(3,jatm) + svb(2,iatm)*denxb(3,jatm)+
     $     ra*denxxa(2,iatm,3,jatm) + rb*denxxb(2,iatm,3,jatm)
          hess(3,iatm,3,jatm) = hess(3,iatm,3,jatm) +
     $     sva(3,iatm)*denxa(3,jatm) + svb(3,iatm)*denxb(3,jatm)+
     $     ra*denxxa(3,iatm,3,jatm) + rb*denxxb(3,iatm,3,jatm)
        EndDo
      EndDo
c
c  --  add weight derivatives contribution
c
c      the sva array is used to store a partial sum
c
      do ia=1,natoms
        sva(1,ia)=prap*denxa(1,ia)+prbp*denxb(1,ia)
        sva(2,ia)=prap*denxa(2,ia)+prbp*denxb(2,ia)
        sva(3,ia)=prap*denxa(3,ia)+prbp*denxb(3,ia)
      enddo
      do ia=1,natoms
      if(ia.ne.icntr)then
        do ja=ia,natoms
        if(ja.ne.icntr)then
          hess(1,ia,1,ja)=hess(1,ia,1,ja)+gwt(1,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(1,ia)+hwt(1,ia,1,ja,ipp)
          hess(2,ia,1,ja)=hess(2,ia,1,ja)+gwt(2,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(2,ia)+hwt(2,ia,1,ja,ipp)
          hess(3,ia,1,ja)=hess(3,ia,1,ja)+gwt(3,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(3,ia)+hwt(3,ia,1,ja,ipp)
          hess(1,ia,2,ja)=hess(1,ia,2,ja)+gwt(1,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(1,ia)+hwt(1,ia,2,ja,ipp)
          hess(2,ia,2,ja)=hess(2,ia,2,ja)+gwt(2,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(2,ia)+hwt(2,ia,2,ja,ipp)
          hess(3,ia,2,ja)=hess(3,ia,2,ja)+gwt(3,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(3,ia)+hwt(3,ia,2,ja,ipp)
          hess(1,ia,3,ja)=hess(1,ia,3,ja)+gwt(1,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(1,ia)+hwt(1,ia,3,ja,ipp)
          hess(2,ia,3,ja)=hess(2,ia,3,ja)+gwt(2,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(2,ia)+hwt(2,ia,3,ja,ipp)
          hess(3,ia,3,ja)=hess(3,ia,3,ja)+gwt(3,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(3,ia)+hwt(3,ia,3,ja,ipp)
        endif
        enddo
      endif
      enddo
      call secund(t4)
      tqh=tqh+t4-t3
c
c   local dft ends here
c
 50   continue
cc
      else
cc
c
c   non-local dft
c
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        prap=pra(ip)
        prbp=prb(ip)
        ra = prap*wg
        rb = prbp*wg
        pgap = pga(ip)
        pgbp = pgb(ip)
        pgcp = pgc(ip)
        pgap2=pgap+pgap
        pgbp2=pgbp+pgbp
        ga = pgap*wg
        gb = pgbp*wg
        gc = pgcp*wg
        rara = prara(ip)*wg
        rbrb = prbrb(ip)*wg
        rarb = prarb(ip)*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        rbga = prbga(ip)*wg
        rbgb = prbgb(ip)*wg
        rbgc = prbgc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gbgb = pgbgb(ip)*wg
        gbgc = pgbgc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  density gradient at current point
c
        DAX=GDENA(1,IPP)
        DAY=GDENA(2,IPP)
        DAZ=GDENA(3,IPP)
        DBX=GDENB(1,IPP)
        DBY=GDENB(2,IPP)
        DBZ=GDENB(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenXA,NAt3)
        CALL ZeroIT(DenXXA,NAt3**2)
        CALL ZeroIT(GDXA,3*NAt3)
        CALL ZeroIT(GDXXA,3*NAt3**2)
        CALL ZeroIT(DenXB,NAt3)
        CALL ZeroIT(DenXXB,NAt3**2)
        CALL ZeroIT(GDXB,3*NAt3)
        CALL ZeroIT(GDXXB,3*NAt3**2)
c   initializations for weight derivatives
c
        gwtmax=zero
        do ia=1,natoms
          if(ia.ne.icntr)then
            gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +        abs(gwt(3,ia,ipp)))
          endif
        enddo
c
c    form the atomic gradient and Hessian  of density
c    and density gradient at current grid point
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c    in addition, we check whether the second derivatives of the
c    density gradient should be computed at all.
c
        abgmax=max(abs(ga),abs(gb),abs(gc))
        dogradxx=abgmax.gt.epsi
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          dmx2 = DM(II)*two   ! max. element of first order density
          if(dmx2.gt.epsi)then
             thtest=thrsh/(vmx*dmx2)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          if(iatm.eq.icntr)goto 220
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2=two*DA(IJ)
            DBIJ2=two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
              DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
              DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
              DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
              DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
              DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VALZZ
c
              DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VALXX
              xyc=DBIJ*VALXY
              DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
              DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
              DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VALYY
              xzc=DBIJ*VALXZ
              DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
              DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
              yzc=DBIJ*VALYZ
              DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
              DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
              DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
c  alpha
              GDXXA(1,1,IAtm,1,IAtm)=GDXXA(1,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXA(1,1,IAtm,2,IAtm)=GDXXA(1,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,1,IAtm)=GDXXA(1,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,2,IAtm)=GDXXA(1,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXA(1,1,IAtm,3,IAtm)=GDXXA(1,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,3,IAtm,1,IAtm)=GDXXA(1,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,2,IAtm,3,IAtm)=GDXXA(1,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,2,IAtm)=GDXXA(1,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,3,IAtm)=GDXXA(1,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXA(2,1,IAtm,1,IAtm)=GDXXA(2,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXA(2,1,IAtm,2,IAtm)=GDXXA(2,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,1,IAtm)=GDXXA(2,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,2,IAtm)=GDXXA(2,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXA(2,1,IAtm,3,IAtm)=GDXXA(2,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,3,IAtm,1,IAtm)=GDXXA(2,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,2,IAtm,3,IAtm)=GDXXA(2,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,2,IAtm)=GDXXA(2,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,3,IAtm)=GDXXA(2,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXA(3,1,IAtm,1,IAtm)=GDXXA(3,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXA(3,1,IAtm,2,IAtm)=GDXXA(3,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,1,IAtm)=GDXXA(3,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,2,IAtm)=GDXXA(3,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXA(3,1,IAtm,3,IAtm)=GDXXA(3,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,3,IAtm,1,IAtm)=GDXXA(3,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,2,IAtm,3,IAtm)=GDXXA(3,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,2,IAtm)=GDXXA(3,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,3,IAtm)=GDXXA(3,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
c  beta
              GDXXB(1,1,IAtm,1,IAtm)=GDXXB(1,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXB(1,1,IAtm,2,IAtm)=GDXXB(1,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,1,IAtm)=GDXXB(1,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,2,IAtm)=GDXXB(1,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXB(1,1,IAtm,3,IAtm)=GDXXB(1,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,3,IAtm,1,IAtm)=GDXXB(1,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,2,IAtm,3,IAtm)=GDXXB(1,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,2,IAtm)=GDXXB(1,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,3,IAtm)=GDXXB(1,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXB(2,1,IAtm,1,IAtm)=GDXXB(2,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXB(2,1,IAtm,2,IAtm)=GDXXB(2,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,1,IAtm)=GDXXB(2,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,2,IAtm)=GDXXB(2,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXB(2,1,IAtm,3,IAtm)=GDXXB(2,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,3,IAtm,1,IAtm)=GDXXB(2,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,2,IAtm,3,IAtm)=GDXXB(2,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,2,IAtm)=GDXXB(2,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,3,IAtm)=GDXXB(2,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXB(3,1,IAtm,1,IAtm)=GDXXB(3,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXB(3,1,IAtm,2,IAtm)=GDXXB(3,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,1,IAtm)=GDXXB(3,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,2,IAtm)=GDXXB(3,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXB(3,1,IAtm,3,IAtm)=GDXXB(3,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,3,IAtm,1,IAtm)=GDXXB(3,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,2,IAtm,3,IAtm)=GDXXB(3,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,2,IAtm)=GDXXB(3,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,3,IAtm)=GDXXB(3,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
c
            If(JAtm.LT.IAtm) GO TO 218
            If(JAtm.eq.icntr) GO TO 218
c
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJ2*VALJX
              DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJY
              DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJZ
              DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
              DBIJ = DBIJ2*VALJX
              DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJY
              DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJZ
              DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
c    alpha
              GDXXA(1,1,IAtm,1,JAtm)=GDXXA(1,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXA(1,1,IAtm,2,JAtm)=GDXXA(1,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXA(1,2,IAtm,1,JAtm)=GDXXA(1,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXA(1,2,IAtm,2,JAtm)=GDXXA(1,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXA(1,1,IAtm,3,JAtm)=GDXXA(1,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXA(1,3,IAtm,1,JAtm)=GDXXA(1,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXA(1,2,IAtm,3,JAtm)=GDXXA(1,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXA(1,3,IAtm,2,JAtm)=GDXXA(1,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXA(1,3,IAtm,3,JAtm)=GDXXA(1,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXA(2,1,IAtm,1,JAtm)=GDXXA(2,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXA(2,1,IAtm,2,JAtm)=GDXXA(2,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXA(2,2,IAtm,1,JAtm)=GDXXA(2,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXA(2,2,IAtm,2,JAtm)=GDXXA(2,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXA(2,1,IAtm,3,JAtm)=GDXXA(2,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXA(2,3,IAtm,1,JAtm)=GDXXA(2,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXA(2,2,IAtm,3,JAtm)=GDXXA(2,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXA(2,3,IAtm,2,JAtm)=GDXXA(2,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXA(2,3,IAtm,3,JAtm)=GDXXA(2,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXA(3,1,IAtm,1,JAtm)=GDXXA(3,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXA(3,1,IAtm,2,JAtm)=GDXXA(3,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXA(3,2,IAtm,1,JAtm)=GDXXA(3,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXA(3,2,IAtm,2,JAtm)=GDXXA(3,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXA(3,1,IAtm,3,JAtm)=GDXXA(3,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXA(3,3,IAtm,1,JAtm)=GDXXA(3,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXA(3,2,IAtm,3,JAtm)=GDXXA(3,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXA(3,3,IAtm,2,JAtm)=GDXXA(3,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXA(3,3,IAtm,3,JAtm)=GDXXA(3,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
c    beta
              GDXXB(1,1,IAtm,1,JAtm)=GDXXB(1,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXB(1,1,IAtm,2,JAtm)=GDXXB(1,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXB(1,2,IAtm,1,JAtm)=GDXXB(1,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXB(1,2,IAtm,2,JAtm)=GDXXB(1,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXB(1,1,IAtm,3,JAtm)=GDXXB(1,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXB(1,3,IAtm,1,JAtm)=GDXXB(1,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXB(1,2,IAtm,3,JAtm)=GDXXB(1,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXB(1,3,IAtm,2,JAtm)=GDXXB(1,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXB(1,3,IAtm,3,JAtm)=GDXXB(1,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXB(2,1,IAtm,1,JAtm)=GDXXB(2,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXB(2,1,IAtm,2,JAtm)=GDXXB(2,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXB(2,2,IAtm,1,JAtm)=GDXXB(2,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXB(2,2,IAtm,2,JAtm)=GDXXB(2,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXB(2,1,IAtm,3,JAtm)=GDXXB(2,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXB(2,3,IAtm,1,JAtm)=GDXXB(2,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXB(2,2,IAtm,3,JAtm)=GDXXB(2,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXB(2,3,IAtm,2,JAtm)=GDXXB(2,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXB(2,3,IAtm,3,JAtm)=GDXXB(2,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXB(3,1,IAtm,1,JAtm)=GDXXB(3,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXB(3,1,IAtm,2,JAtm)=GDXXB(3,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXB(3,2,IAtm,1,JAtm)=GDXXB(3,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXB(3,2,IAtm,2,JAtm)=GDXXB(3,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXB(3,1,IAtm,3,JAtm)=GDXXB(3,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXB(3,3,IAtm,1,JAtm)=GDXXB(3,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXB(3,2,IAtm,3,JAtm)=GDXXB(3,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXB(3,3,IAtm,2,JAtm)=GDXXB(3,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXB(3,3,IAtm,3,JAtm)=GDXXB(3,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2= two*DA(IJ)
            DBIJ2= two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
              DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
              DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
              DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
              DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
              DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VALZZ
c
              DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VALXX
              xyc=DBIJ*VALXY
              DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
              DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
              DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VALYY
              xzc=DBIJ*VALXZ
              DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
              DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
              yzc=DBIJ*VALYZ
              DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
              DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
              DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
c    alpha
              GDXXA(1,1,IAtm,1,IAtm)=GDXXA(1,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXA(1,1,IAtm,2,IAtm)=GDXXA(1,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,1,IAtm)=GDXXA(1,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,2,IAtm)=GDXXA(1,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXA(1,1,IAtm,3,IAtm)=GDXXA(1,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,3,IAtm,1,IAtm)=GDXXA(1,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,2,IAtm,3,IAtm)=GDXXA(1,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,2,IAtm)=GDXXA(1,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,3,IAtm)=GDXXA(1,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXA(2,1,IAtm,1,IAtm)=GDXXA(2,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXA(2,1,IAtm,2,IAtm)=GDXXA(2,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,1,IAtm)=GDXXA(2,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,2,IAtm)=GDXXA(2,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXA(2,1,IAtm,3,IAtm)=GDXXA(2,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,3,IAtm,1,IAtm)=GDXXA(2,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,2,IAtm,3,IAtm)=GDXXA(2,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,2,IAtm)=GDXXA(2,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,3,IAtm)=GDXXA(2,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXA(3,1,IAtm,1,IAtm)=GDXXA(3,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXA(3,1,IAtm,2,IAtm)=GDXXA(3,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,1,IAtm)=GDXXA(3,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,2,IAtm)=GDXXA(3,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXA(3,1,IAtm,3,IAtm)=GDXXA(3,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,3,IAtm,1,IAtm)=GDXXA(3,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,2,IAtm,3,IAtm)=GDXXA(3,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,2,IAtm)=GDXXA(3,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,3,IAtm)=GDXXA(3,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
c    beta
              GDXXB(1,1,IAtm,1,IAtm)=GDXXB(1,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXB(1,1,IAtm,2,IAtm)=GDXXB(1,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,1,IAtm)=GDXXB(1,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,2,IAtm)=GDXXB(1,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXB(1,1,IAtm,3,IAtm)=GDXXB(1,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,3,IAtm,1,IAtm)=GDXXB(1,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,2,IAtm,3,IAtm)=GDXXB(1,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,2,IAtm)=GDXXB(1,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,3,IAtm)=GDXXB(1,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXB(2,1,IAtm,1,IAtm)=GDXXB(2,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXB(2,1,IAtm,2,IAtm)=GDXXB(2,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,1,IAtm)=GDXXB(2,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,2,IAtm)=GDXXB(2,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXB(2,1,IAtm,3,IAtm)=GDXXB(2,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,3,IAtm,1,IAtm)=GDXXB(2,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,2,IAtm,3,IAtm)=GDXXB(2,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,2,IAtm)=GDXXB(2,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,3,IAtm)=GDXXB(2,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXB(3,1,IAtm,1,IAtm)=GDXXB(3,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXB(3,1,IAtm,2,IAtm)=GDXXB(3,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,1,IAtm)=GDXXB(3,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,2,IAtm)=GDXXB(3,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXB(3,1,IAtm,3,IAtm)=GDXXB(3,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,3,IAtm,1,IAtm)=GDXXB(3,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,2,IAtm,3,IAtm)=GDXXB(3,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,2,IAtm)=GDXXB(3,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,3,IAtm)=GDXXB(3,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
            If(JAtm.LT.IAtm) GO TO 219
            If(JAtm.EQ.ICntr) GO TO 219
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJ2*VALJX
              DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJY
              DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJZ
              DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
              DBIJ = DBIJ2*VALJX
              DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJY
              DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJZ
              DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
c     alpha
              GDXXA(1,1,IAtm,1,JAtm)=GDXXA(1,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXA(1,1,IAtm,2,JAtm)=GDXXA(1,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXA(1,2,IAtm,1,JAtm)=GDXXA(1,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXA(1,2,IAtm,2,JAtm)=GDXXA(1,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXA(1,1,IAtm,3,JAtm)=GDXXA(1,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXA(1,3,IAtm,1,JAtm)=GDXXA(1,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXA(1,2,IAtm,3,JAtm)=GDXXA(1,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXA(1,3,IAtm,2,JAtm)=GDXXA(1,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXA(1,3,IAtm,3,JAtm)=GDXXA(1,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXA(2,1,IAtm,1,JAtm)=GDXXA(2,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXA(2,1,IAtm,2,JAtm)=GDXXA(2,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXA(2,2,IAtm,1,JAtm)=GDXXA(2,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXA(2,2,IAtm,2,JAtm)=GDXXA(2,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXA(2,1,IAtm,3,JAtm)=GDXXA(2,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXA(2,3,IAtm,1,JAtm)=GDXXA(2,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXA(2,2,IAtm,3,JAtm)=GDXXA(2,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXA(2,3,IAtm,2,JAtm)=GDXXA(2,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXA(2,3,IAtm,3,JAtm)=GDXXA(2,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXA(3,1,IAtm,1,JAtm)=GDXXA(3,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXA(3,1,IAtm,2,JAtm)=GDXXA(3,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXA(3,2,IAtm,1,JAtm)=GDXXA(3,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXA(3,2,IAtm,2,JAtm)=GDXXA(3,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXA(3,1,IAtm,3,JAtm)=GDXXA(3,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXA(3,3,IAtm,1,JAtm)=GDXXA(3,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXA(3,2,IAtm,3,JAtm)=GDXXA(3,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXA(3,3,IAtm,2,JAtm)=GDXXA(3,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXA(3,3,IAtm,3,JAtm)=GDXXA(3,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
c     beta
              GDXXB(1,1,IAtm,1,JAtm)=GDXXB(1,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXB(1,1,IAtm,2,JAtm)=GDXXB(1,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXB(1,2,IAtm,1,JAtm)=GDXXB(1,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXB(1,2,IAtm,2,JAtm)=GDXXB(1,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXB(1,1,IAtm,3,JAtm)=GDXXB(1,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXB(1,3,IAtm,1,JAtm)=GDXXB(1,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXB(1,2,IAtm,3,JAtm)=GDXXB(1,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXB(1,3,IAtm,2,JAtm)=GDXXB(1,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXB(1,3,IAtm,3,JAtm)=GDXXB(1,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXB(2,1,IAtm,1,JAtm)=GDXXB(2,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXB(2,1,IAtm,2,JAtm)=GDXXB(2,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXB(2,2,IAtm,1,JAtm)=GDXXB(2,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXB(2,2,IAtm,2,JAtm)=GDXXB(2,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXB(2,1,IAtm,3,JAtm)=GDXXB(2,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXB(2,3,IAtm,1,JAtm)=GDXXB(2,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXB(2,2,IAtm,3,JAtm)=GDXXB(2,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXB(2,3,IAtm,2,JAtm)=GDXXB(2,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXB(2,3,IAtm,3,JAtm)=GDXXB(2,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXB(3,1,IAtm,1,JAtm)=GDXXB(3,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXB(3,1,IAtm,2,JAtm)=GDXXB(3,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXB(3,2,IAtm,1,JAtm)=GDXXB(3,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXB(3,2,IAtm,2,JAtm)=GDXXB(3,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXB(3,1,IAtm,3,JAtm)=GDXXB(3,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXB(3,3,IAtm,1,JAtm)=GDXXB(3,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXB(3,2,IAtm,3,JAtm)=GDXXB(3,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXB(3,3,IAtm,2,JAtm)=GDXXB(3,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXB(3,3,IAtm,3,JAtm)=GDXXB(3,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and Hessian of density gradient invariant
c
c
        DO IAtm=1,NAtoms
c   alpha alpha
          GXAA(1,IAtm)=two*(DAX*GDXA(1,1,IAtm)+
     $                   DAY*GDXA(2,1,IAtm)+DAZ*GDXA(3,1,IAtm))
          GXAA(2,IAtm)=two*(DAX*GDXA(1,2,IAtm)+
     $                   DAY*GDXA(2,2,IAtm)+DAZ*GDXA(3,2,IAtm))
          GXAA(3,IAtm)=two*(DAX*GDXA(1,3,IAtm)+
     $                   DAY*GDXA(2,3,IAtm)+DAZ*GDXA(3,3,IAtm))
c   beta beta
          GXBB(1,IAtm)=two*(DBX*GDXB(1,1,IAtm)+
     $                   DBY*GDXB(2,1,IAtm)+DBZ*GDXB(3,1,IAtm))
          GXBB(2,IAtm)=two*(DBX*GDXB(1,2,IAtm)+
     $                   DBY*GDXB(2,2,IAtm)+DBZ*GDXB(3,2,IAtm))
          GXBB(3,IAtm)=two*(DBX*GDXB(1,3,IAtm)+
     $                   DBY*GDXB(2,3,IAtm)+DBZ*GDXB(3,3,IAtm))
c   alpha beta
          GXAB(1,IAtm)=DAX*GDXB(1,1,IAtm)+DBX*GDXA(1,1,IAtm)+
     $                 DAY*GDXB(2,1,IAtm)+DBY*GDXA(2,1,IAtm)+
     $                 DAZ*GDXB(3,1,IAtm)+DBZ*GDXA(3,1,IAtm)
          GXAB(2,IAtm)=DAX*GDXB(1,2,IAtm)+DBX*GDXA(1,2,IAtm)+
     $                 DAY*GDXB(2,2,IAtm)+DBY*GDXA(2,2,IAtm)+
     $                 DAZ*GDXB(3,2,IAtm)+DBZ*GDXA(3,2,IAtm)
          GXAB(3,IAtm)=DAX*GDXB(1,3,IAtm)+DBX*GDXA(1,3,IAtm)+
     $                 DAY*GDXB(2,3,IAtm)+DBY*GDXA(2,3,IAtm)+
     $                 DAZ*GDXB(3,3,IAtm)+DBZ*GDXA(3,3,IAtm)
        EndDO
c
c  compute second derivatives of gradient invariants
c  only of needed
c
        if(dogradxx)then
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
c     alpha alpha
            GXXAA(1,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,1,JAtm)+GDXA(1,1,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,1,JAtm)+GDXA(2,1,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,1,JAtm)+GDXA(3,1,IAtm)*GDXA(3,1,JAtm))
            GXXAA(1,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,2,JAtm)+GDXA(1,1,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,2,JAtm)+GDXA(2,1,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,2,JAtm)+GDXA(3,1,IAtm)*GDXA(3,2,JAtm))
            GXXAA(2,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,1,JAtm)+GDXA(1,2,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,1,JAtm)+GDXA(2,2,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,1,JAtm)+GDXA(3,2,IAtm)*GDXA(3,1,JAtm))
            GXXAA(2,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,2,JAtm)+GDXA(1,2,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,2,JAtm)+GDXA(2,2,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,2,JAtm)+GDXA(3,2,IAtm)*GDXA(3,2,JAtm))
            GXXAA(1,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,3,JAtm)+GDXA(1,1,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,3,JAtm)+GDXA(2,1,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,3,JAtm)+GDXA(3,1,IAtm)*GDXA(3,3,JAtm))
            GXXAA(3,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,1,JAtm)+GDXA(1,3,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,1,JAtm)+GDXA(2,3,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,1,JAtm)+GDXA(3,3,IAtm)*GDXA(3,1,JAtm))
            GXXAA(2,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,3,JAtm)+GDXA(1,2,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,3,JAtm)+GDXA(2,2,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,3,JAtm)+GDXA(3,2,IAtm)*GDXA(3,3,JAtm))
            GXXAA(3,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,2,JAtm)+GDXA(1,3,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,2,JAtm)+GDXA(2,3,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,2,JAtm)+GDXA(3,3,IAtm)*GDXA(3,2,JAtm))
            GXXAA(3,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,3,JAtm)+GDXA(1,3,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,3,JAtm)+GDXA(2,3,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,3,JAtm)+GDXA(3,3,IAtm)*GDXA(3,3,JAtm))
c     beta beta
            GXXBB(1,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,1,JAtm)+GDXB(1,1,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,1,JAtm)+GDXB(2,1,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,1,JAtm)+GDXB(3,1,IAtm)*GDXB(3,1,JAtm))
            GXXBB(1,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,2,JAtm)+GDXB(1,1,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,2,JAtm)+GDXB(2,1,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,2,JAtm)+GDXB(3,1,IAtm)*GDXB(3,2,JAtm))
            GXXBB(2,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,1,JAtm)+GDXB(1,2,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,1,JAtm)+GDXB(2,2,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,1,JAtm)+GDXB(3,2,IAtm)*GDXB(3,1,JAtm))
            GXXBB(2,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,2,JAtm)+GDXB(1,2,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,2,JAtm)+GDXB(2,2,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,2,JAtm)+GDXB(3,2,IAtm)*GDXB(3,2,JAtm))
            GXXBB(1,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,3,JAtm)+GDXB(1,1,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,3,JAtm)+GDXB(2,1,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,3,JAtm)+GDXB(3,1,IAtm)*GDXB(3,3,JAtm))
            GXXBB(3,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,1,JAtm)+GDXB(1,3,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,1,JAtm)+GDXB(2,3,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,1,JAtm)+GDXB(3,3,IAtm)*GDXB(3,1,JAtm))
            GXXBB(2,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,3,JAtm)+GDXB(1,2,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,3,JAtm)+GDXB(2,2,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,3,JAtm)+GDXB(3,2,IAtm)*GDXB(3,3,JAtm))
            GXXBB(3,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,2,JAtm)+GDXB(1,3,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,2,JAtm)+GDXB(2,3,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,2,JAtm)+GDXB(3,3,IAtm)*GDXB(3,2,JAtm))
            GXXBB(3,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,3,JAtm)+GDXB(1,3,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,3,JAtm)+GDXB(2,3,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,3,JAtm)+GDXB(3,3,IAtm)*GDXB(3,3,JAtm))
          EndDO
        EndDO
c
c     alpha beta is done separately, as only few functionals
c     actually need this one
c
        if(abs(gc).gt.epsi)then
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
            GXXAB(1,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,1,JAtm)+GDXA(1,1,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,1,JAtm)+GDXB(1,1,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,1,JAtm)+GDXA(2,1,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,1,JAtm)+GDXB(2,1,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,1,JAtm)+GDXA(3,1,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,1,JAtm)+GDXB(3,1,IAtm)*GDXA(3,1,JAtm)
            GXXAB(2,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,1,JAtm)+GDXA(1,2,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,1,JAtm)+GDXB(1,2,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,1,JAtm)+GDXA(2,2,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,1,JAtm)+GDXB(2,2,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,1,JAtm)+GDXA(3,2,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,1,JAtm)+GDXB(3,2,IAtm)*GDXA(3,1,JAtm)
            GXXAB(3,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,1,JAtm)+GDXA(1,3,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,1,JAtm)+GDXB(1,3,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,1,JAtm)+GDXA(2,3,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,1,JAtm)+GDXB(2,3,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,1,JAtm)+GDXA(3,3,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,1,JAtm)+GDXB(3,3,IAtm)*GDXA(3,1,JAtm)
            GXXAB(1,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,2,JAtm)+GDXA(1,1,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,2,JAtm)+GDXB(1,1,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,2,JAtm)+GDXA(2,1,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,2,JAtm)+GDXB(2,1,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,2,JAtm)+GDXA(3,1,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,2,JAtm)+GDXB(3,1,IAtm)*GDXA(3,2,JAtm)
            GXXAB(2,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,2,JAtm)+GDXA(1,2,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,2,JAtm)+GDXB(1,2,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,2,JAtm)+GDXA(2,2,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,2,JAtm)+GDXB(2,2,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,2,JAtm)+GDXA(3,2,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,2,JAtm)+GDXB(3,2,IAtm)*GDXA(3,2,JAtm)
            GXXAB(3,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,2,JAtm)+GDXA(1,3,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,2,JAtm)+GDXB(1,3,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,2,JAtm)+GDXA(2,3,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,2,JAtm)+GDXB(2,3,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,2,JAtm)+GDXA(3,3,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,2,JAtm)+GDXB(3,3,IAtm)*GDXA(3,2,JAtm)
            GXXAB(1,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,3,JAtm)+GDXA(1,1,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,3,JAtm)+GDXB(1,1,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,3,JAtm)+GDXA(2,1,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,3,JAtm)+GDXB(2,1,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,3,JAtm)+GDXA(3,1,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,3,JAtm)+GDXB(3,1,IAtm)*GDXA(3,3,JAtm)
            GXXAB(2,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,3,JAtm)+GDXA(1,2,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,3,JAtm)+GDXB(1,2,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,3,JAtm)+GDXA(2,2,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,3,JAtm)+GDXB(2,2,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,3,JAtm)+GDXA(3,2,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,3,JAtm)+GDXB(3,2,IAtm)*GDXA(3,3,JAtm)
            GXXAB(3,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,3,JAtm)+GDXA(1,3,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,3,JAtm)+GDXB(1,3,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,3,JAtm)+GDXA(2,3,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,3,JAtm)+GDXB(2,3,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,3,JAtm)+GDXA(3,3,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,3,JAtm)+GDXB(3,3,IAtm)*GDXA(3,3,JAtm)
          EndDO
        EndDO
        else
          CALL ZeroIT(GXXAB,NAt3**2)
        endif
        else
          CALL ZeroIT(GXXAA,NAt3**2)
          CALL ZeroIT(GXXBB,NAt3**2)
          CALL ZeroIT(GXXAB,NAt3**2)
        endif
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
C    and also compute partial sums used in Hessian quadrature
c
        ga2= ga+ga
        gb2= gb+gb
        DAX2=DAX+DAX
        DAY2=DAY+DAY
        DAZ2=DAZ+DAZ
        DBX2=DBX+DBX
        DBY2=DBY+DBY
        DBZ2=DBZ+DBZ
c
c     coeff. X
c
        SXXA=ga2*DAX+gc*DBX
        SXYA=ga2*DAY+gc*DBY
        SXZA=ga2*DAZ+gc*DBZ
        SXXB=gb2*DBX+gc*DAX
        SXYB=gb2*DBY+gc*DAY
        SXZB=gb2*DBZ+gc*DAZ
c
        DO IAtm=1,NAtoms
c
c     coeff. V
c
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)+
     $        raga*gxaa(1,iatm)+ragb*gxbb(1,iatm)+ragc*gxab(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)+
     $        raga*gxaa(2,iatm)+ragb*gxbb(2,iatm)+ragc*gxab(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)+
     $        raga*gxaa(3,iatm)+ragb*gxbb(3,iatm)+ragc*gxab(3,iatm)
          svb(1,iatm)=rarb*denxa(1,iatm)+rbrb*denxb(1,iatm)+
     $        rbga*gxaa(1,iatm)+rbgb*gxbb(1,iatm)+rbgc*gxab(1,iatm)
          svb(2,iatm)=rarb*denxa(2,iatm)+rbrb*denxb(2,iatm)+
     $        rbga*gxaa(2,iatm)+rbgb*gxbb(2,iatm)+rbgc*gxab(2,iatm)
          svb(3,iatm)=rarb*denxa(3,iatm)+rbrb*denxb(3,iatm)+
     $        rbga*gxaa(3,iatm)+rbgb*gxbb(3,iatm)+rbgc*gxab(3,iatm)
c
c  partial sums for Hessian
c
          sga(1,iatm)=raga*denxa(1,iatm)+rbga*denxb(1,iatm)+
     $        gaga*gxaa(1,iatm)+gagb*gxbb(1,iatm)+gagc*gxab(1,iatm)
          sga(2,iatm)=raga*denxa(2,iatm)+rbga*denxb(2,iatm)+
     $        gaga*gxaa(2,iatm)+gagb*gxbb(2,iatm)+gagc*gxab(2,iatm)
          sga(3,iatm)=raga*denxa(3,iatm)+rbga*denxb(3,iatm)+
     $        gaga*gxaa(3,iatm)+gagb*gxbb(3,iatm)+gagc*gxab(3,iatm)
c
          sgb(1,iatm)=ragb*denxa(1,iatm)+rbgb*denxb(1,iatm)+
     $        gagb*gxaa(1,iatm)+gbgb*gxbb(1,iatm)+gbgc*gxab(1,iatm)
          sgb(2,iatm)=ragb*denxa(2,iatm)+rbgb*denxb(2,iatm)+
     $        gagb*gxaa(2,iatm)+gbgb*gxbb(2,iatm)+gbgc*gxab(2,iatm)
          sgb(3,iatm)=ragb*denxa(3,iatm)+rbgb*denxb(3,iatm)+
     $        gagb*gxaa(3,iatm)+gbgb*gxbb(3,iatm)+gbgc*gxab(3,iatm)
c
          sgc(1,iatm)=ragc*denxa(1,iatm)+rbgc*denxb(1,iatm)+
     $        gagc*gxaa(1,iatm)+gbgc*gxbb(1,iatm)+gcgc*gxab(1,iatm)
          sgc(2,iatm)=ragc*denxa(2,iatm)+rbgc*denxb(2,iatm)+
     $        gagc*gxaa(2,iatm)+gbgc*gxbb(2,iatm)+gcgc*gxab(2,iatm)
          sgc(3,iatm)=ragc*denxa(3,iatm)+rbgc*denxb(3,iatm)+
     $        gagc*gxaa(3,iatm)+gbgc*gxbb(3,iatm)+gcgc*gxab(3,iatm)
c
c     coeff. W
c
          swa(1,1,iatm)=ga2*gdxa(1,1,iatm)+gc*gdxb(1,1,Iatm)+
     $        DAX2*sga(1,iatm)+DBX*sgc(1,iatm)
          swa(1,2,iatm)=ga2*gdxa(1,2,iatm)+gc*gdxb(1,2,Iatm)+
     $        DAX2*sga(2,iatm)+DBX*sgc(2,iatm)
          swa(1,3,iatm)=ga2*gdxa(1,3,iatm)+gc*gdxb(1,3,Iatm)+
     $        DAX2*sga(3,iatm)+DBX*sgc(3,iatm)
          swa(2,1,iatm)=ga2*gdxa(2,1,iatm)+gc*gdxb(2,1,Iatm)+
     $        DAY2*sga(1,iatm)+DBY*sgc(1,iatm)
          swa(2,2,iatm)=ga2*gdxa(2,2,iatm)+gc*gdxb(2,2,Iatm)+
     $        DAY2*sga(2,iatm)+DBY*sgc(2,iatm)
          swa(2,3,iatm)=ga2*gdxa(2,3,iatm)+gc*gdxb(2,3,Iatm)+
     $        DAY2*sga(3,iatm)+DBY*sgc(3,iatm)
          swa(3,1,iatm)=ga2*gdxa(3,1,iatm)+gc*gdxb(3,1,Iatm)+
     $        DAZ2*sga(1,iatm)+DBZ*sgc(1,iatm)
          swa(3,2,iatm)=ga2*gdxa(3,2,iatm)+gc*gdxb(3,2,Iatm)+
     $        DAZ2*sga(2,iatm)+DBZ*sgc(2,iatm)
          swa(3,3,iatm)=ga2*gdxa(3,3,iatm)+gc*gdxb(3,3,Iatm)+
     $        DAZ2*sga(3,iatm)+DBZ*sgc(3,iatm)
c
          swb(1,1,iatm)=gb2*gdxb(1,1,iatm)+gc*gdxa(1,1,Iatm)+
     $        DBX2*sgb(1,iatm)+DAX*sgc(1,iatm)
          swb(1,2,iatm)=gb2*gdxb(1,2,iatm)+gc*gdxa(1,2,Iatm)+
     $        DBX2*sgb(2,iatm)+DAX*sgc(2,iatm)
          swb(1,3,iatm)=gb2*gdxb(1,3,iatm)+gc*gdxa(1,3,Iatm)+
     $        DBX2*sgb(3,iatm)+DAX*sgc(3,iatm)
          swb(2,1,iatm)=gb2*gdxb(2,1,iatm)+gc*gdxa(2,1,Iatm)+
     $        DBY2*sgb(1,iatm)+DAY*sgc(1,iatm)
          swb(2,2,iatm)=gb2*gdxb(2,2,iatm)+gc*gdxa(2,2,Iatm)+
     $        DBY2*sgb(2,iatm)+DAY*sgc(2,iatm)
          swb(2,3,iatm)=gb2*gdxb(2,3,iatm)+gc*gdxa(2,3,Iatm)+
     $        DBY2*sgb(3,iatm)+DAY*sgc(3,iatm)
          swb(3,1,iatm)=gb2*gdxb(3,1,iatm)+gc*gdxa(3,1,Iatm)+
     $        DBZ2*sgb(1,iatm)+DAZ*sgc(1,iatm)
          swb(3,2,iatm)=gb2*gdxb(3,2,iatm)+gc*gdxa(3,2,Iatm)+
     $        DBZ2*sgb(2,iatm)+DAZ*sgc(2,iatm)
          swb(3,3,iatm)=gb2*gdxb(3,3,iatm)+gc*gdxa(3,3,Iatm)+
     $        DBZ2*sgb(3,iatm)+DAZ*sgc(3,iatm)
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    fock matrix.
c
c -- get the maximum absolute value of coefficients v and w
c
        Call absmax(NAt3,sva,iixx,svamax)
        Call absmax(NAt3,svb,iixx,svbmax)
        svmax=max(svamax,svbmax)
        Call absmax(3*NAt3,swa,iixx,swamax)
        Call absmax(3*NAt3,swb,iixx,swbmax)
        swmax=max(swamax,swbmax)
        tswmax=Three*swmax
c
c -- global threshold testing
c
        abr=max(abs(ra),abs(rb))
        vmx2=vmx*vmx
        SMax = max(Abs(SXXA+SXYA+SXZA),abs(SXXB+SXYB+SXZB))
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abr+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVIA=SXXA*ValIX+SXYA*ValIY+SXZA*ValIZ+ra*ValI
          XVIB=SXXB*ValIX+SXYB*ValIY+SXZB*ValIZ+rb*ValI
          SXXIA=SXXA*ValI
          SXYIA=SXYA*ValI
          SXZIA=SXZA*ValI
          SXXIB=SXXB*ValI
          SXYIB=SXYB*ValI
          SXZIB=SXZB*ValI
          XXVXA=SXXA*VAOXX(1,I)+SXYA*VAOXX(2,I)+SXZA*VAOXX(4,I)+ra*VALIX
          XXVYA=SXXA*VAOXX(2,I)+SXYA*VAOXX(3,I)+SXZA*VAOXX(5,I)+ra*VALIY
          XXVZA=SXXA*VAOXX(4,I)+SXYA*VAOXX(5,I)+SXZA*VAOXX(6,I)+ra*VALIZ
          XXVXB=SXXB*VAOXX(1,I)+SXYB*VAOXX(2,I)+SXZB*VAOXX(4,I)+rb*VALIX
          XXVYB=SXXB*VAOXX(2,I)+SXYB*VAOXX(3,I)+SXZB*VAOXX(5,I)+rb*VALIY
          XXVZB=SXXB*VAOXX(4,I)+SXYB*VAOXX(5,I)+SXZB*VAOXX(6,I)+rb*VALIZ
          valm=abs(vali)
          valm1=max(abs(xvia),abs(xvib))
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1a=max(abs(xxvxa),abs(xxvya),abs(xxvza))
          valmx1b=max(abs(xxvxb),abs(xxvyb),abs(xxvzb))
          valmx1=max(valmx1a,valmx1b)
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJA=SXXA*ValJX+SXYA*ValJY+SXZA*ValJZ
              XVJB=SXXB*ValJX+SXYB*ValJY+SXZB*ValJZ
              abxvj=max(abs(xvja),abs(xvjb))
              if(valmx1*abs(valj)+valmx*abxvj.gt.thrsh)then
              if(iatm.ne.icntr) then
                fda(1,iatm,ij)=fda(1,iatm,ij)-xxvxa*valj-valix*xvja
                fda(2,iatm,ij)=fda(2,iatm,ij)-xxvya*valj-valiy*xvja
                fda(3,iatm,ij)=fda(3,iatm,ij)-xxvza*valj-valiz*xvja
                fdb(1,iatm,ij)=fdb(1,iatm,ij)-xxvxb*valj-valix*xvjb
                fdb(2,iatm,ij)=fdb(2,iatm,ij)-xxvyb*valj-valiy*xvjb
                fdb(3,iatm,ij)=fdb(3,iatm,ij)-xxvzb*valj-valiz*xvjb
              endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
              if(jatm.ne.icntr) then
                fda(1,jAtm,ij) = fda(1,jatm,ij) - xvia*Valjx -
     $          vaoxx(1,j)*sxxia - vaoxx(2,j)*sxyia - vaoxx(4,j)*sxzia
                fda(2,jatm,ij) = fda(2,jatm,ij) - xvia*valjy -
     $          vaoxx(2,j)*sxxia - vaoxx(3,j)*sxyia - vaoxx(5,j)*sxzia
                fda(3,jatm,ij) = fda(3,jAtm,ij) - xvia*valjz -
     $          vaoxx(4,j)*sxxia - vaoxx(5,j)*sxyia - vaoxx(6,j)*sxzia
c
                fdb(1,jatm,ij) = fdb(1,jatm,ij) - xvib*valjx -
     $          vaoxx(1,j)*sxxib - vaoxx(2,j)*sxyib - vaoxx(4,j)*sxzib
                fdb(2,jatm,ij) = fdb(2,jatm,ij) - xvib*valjy -
     $          vaoxx(2,j)*sxxib - vaoxx(3,j)*sxyib - vaoxx(5,j)*sxzib
                fdb(3,jatm,ij) = fdb(3,jatm,ij) - xvib*valjz -
     $          vaoxx(4,j)*sxxib - vaoxx(5,j)*sxyib - vaoxx(6,j)*sxzib
              endif
              endif
              vij=vali*valj
              vijx=valix*valj+valjx*vali
              vijy=valiy*valj+valjy*vali
              vijz=valiz*valj+valjz*vali
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                do katm=1,natoms
                  fda(1,katm,ij) = fda(1,katm,ij) + (sva(1,katm)*vij+
     $                           swa(1,1,katm)*vijx+swa(2,1,katm)*vijy+
     $                           swa(3,1,katm)*vijz)
                  fda(2,katm,ij) = fda(2,katm,ij) + (sva(2,katm)*vij+
     $                           swa(1,2,katm)*vijx+swa(2,2,katm)*vijy+
     $                           swa(3,2,katm)*vijz)
                  fda(3,katm,ij) = fda(3,katm,ij) + (sva(3,katm)*vij+
     $                           swa(1,3,katm)*vijx+swa(2,3,katm)*vijy+
     $                           swa(3,3,katm)*vijz)
c
                  fdb(1,katm,ij) = fdb(1,katm,ij) + (svb(1,katm)*vij+
     $                           swb(1,1,katm)*vijx+swb(2,1,katm)*vijy+
     $                           swb(3,1,katm)*vijz)
                  fdb(2,katm,ij) = fdb(2,katm,ij) + (svb(2,katm)*vij+
     $                           swb(1,2,katm)*vijx+swb(2,2,katm)*vijy+
     $                           swb(3,2,katm)*vijz)
                  fdb(3,katm,ij) = fdb(3,katm,ij) + (svb(3,katm)*vij+
     $                           swb(1,3,katm)*vijx+swb(2,3,katm)*vijy+
     $                           swb(3,3,katm)*vijz)
                enddo
              endif
c
c -- weight derivatives contribution
c
              dxija=dax*vijx+day*vijy+daz*vijz
              dxijb=dbx*vijx+dby*vijy+dbz*vijz
              gija=prap*vij+pgap2*dxija+pgcp*dxijb
              gijb=prbp*vij+pgbp2*dxijb+pgcp*dxija
              if(gwtmax*max(abs(gija),abs(gijb)).gt.thrsh)then
                do ia=1,natoms
                  if(ia.ne.icntr)then
                    fda(1,ia,ij)=fda(1,ia,ij)+gwt(1,ia,ipp)*gija
                    fda(2,ia,ij)=fda(2,ia,ij)+gwt(2,ia,ipp)*gija
                    fda(3,ia,ij)=fda(3,ia,ij)+gwt(3,ia,ipp)*gija
                    fdb(1,ia,ij)=fdb(1,ia,ij)+gwt(1,ia,ipp)*gijb
                    fdb(2,ia,ij)=fdb(2,ia,ij)+gwt(2,ia,ipp)*gijb
                    fdb(3,ia,ij)=fdb(3,ia,ij)+gwt(3,ia,ipp)*gijb
                  endif
                enddo
              endif
 230        continue
          endif
 240    continue
c
c   numerical quadrature for direct contribution to hessian matrix
c
 245    continue
        call secund(t5)
        tqf=tqf+t5-t4
c
c -- direct contribution to hessian matrix.
c
        Do IAtm=1,NAtoms
          Do JAtm=IAtm,NAtoms
            hess(1,IAtm,1,JAtm) = hess(1,IAtm,1,JAtm) +
     $       sva(1,IAtm)*DenXA(1,JAtm) + svb(1,IAtm)*DenXB(1,JAtm) +
     $       sga(1,IAtm)*gxaa(1,JAtm)  + sgb(1,IAtm)*gxbb(1,JAtm) +
     $       sgc(1,IAtm)*gxab(1,JAtm)  + ra*DenXXA(1,IAtm,1,JAtm) +
     $       rb*DenXXB(1,IAtm,1,JAtm)  + ga*gxxaa(1,IAtm,1,JAtm) +
     $       gb*gxxbb(1,IAtm,1,JAtm)   + gc*gxxab(1,IAtm,1,JAtm)
            hess(2,IAtm,1,JAtm) = hess(2,IAtm,1,JAtm) +
     $       sva(2,IAtm)*DenXA(1,JAtm) + svb(2,IAtm)*DenXB(1,JAtm) +
     $       sga(2,IAtm)*gxaa(1,JAtm)  + sgb(2,IAtm)*gxbb(1,JAtm) +
     $       sgc(2,IAtm)*gxab(1,JAtm)  + ra*DenXXA(2,IAtm,1,JAtm) +
     $       rb*DenXXB(2,IAtm,1,JAtm)  + ga*gxxaa(2,IAtm,1,JAtm) +
     $       gb*gxxbb(2,IAtm,1,JAtm)   + gc*gxxab(2,IAtm,1,JAtm)
            hess(3,IAtm,1,JAtm) = hess(3,IAtm,1,JAtm) +
     $       sva(3,IAtm)*DenXA(1,JAtm) + svb(3,IAtm)*DenXB(1,JAtm) +
     $       sga(3,IAtm)*gxaa(1,JAtm)  + sgb(3,IAtm)*gxbb(1,JAtm) +
     $       sgc(3,IAtm)*gxab(1,JAtm)  + ra*DenXXA(3,IAtm,1,JAtm) +
     $       rb*DenXXB(3,IAtm,1,JAtm)  + ga*gxxaa(3,IAtm,1,JAtm) +
     $       gb*gxxbb(3,IAtm,1,JAtm)   + gc*gxxab(3,IAtm,1,JAtm)
            hess(1,IAtm,2,JAtm) = hess(1,IAtm,2,JAtm) +
     $       sva(1,IAtm)*DenXA(2,JAtm) + svb(1,IAtm)*DenXB(2,JAtm) +
     $       sga(1,IAtm)*gxaa(2,JAtm)  + sgb(1,IAtm)*gxbb(2,JAtm) +
     $       sgc(1,IAtm)*gxab(2,JAtm)  + ra*DenXXA(1,IAtm,2,JAtm) +
     $       rb*DenXXB(1,IAtm,2,JAtm)  + ga*gxxaa(1,IAtm,2,JAtm) +
     $       gb*gxxbb(1,IAtm,2,JAtm)   + gc*gxxab(1,IAtm,2,JAtm)
            hess(2,IAtm,2,JAtm) = hess(2,IAtm,2,JAtm) +
     $       sva(2,IAtm)*DenXA(2,JAtm) + svb(2,IAtm)*DenXB(2,JAtm) +
     $       sga(2,IAtm)*gxaa(2,JAtm)  + sgb(2,IAtm)*gxbb(2,JAtm) +
     $       sgc(2,IAtm)*gxab(2,JAtm)  + ra*DenXXA(2,IAtm,2,JAtm) +
     $       rb*DenXXB(2,IAtm,2,JAtm)  + ga*gxxaa(2,IAtm,2,JAtm) +
     $       gb*gxxbb(2,IAtm,2,JAtm)   + gc*gxxab(2,IAtm,2,JAtm)
            hess(3,IAtm,2,JAtm) = hess(3,IAtm,2,JAtm) +
     $       sva(3,IAtm)*DenXA(2,JAtm) + svb(3,IAtm)*DenXB(2,JAtm) +
     $       sga(3,IAtm)*gxaa(2,JAtm)  + sgb(3,IAtm)*gxbb(2,JAtm) +
     $       sgc(3,IAtm)*gxab(2,JAtm)  + ra*DenXXA(3,IAtm,2,JAtm) +
     $       rb*DenXXB(3,IAtm,2,JAtm)  + ga*gxxaa(3,IAtm,2,JAtm) +
     $       gb*gxxbb(3,IAtm,2,JAtm)   + gc*gxxab(3,IAtm,2,JAtm)
            hess(1,IAtm,3,JAtm) = hess(1,IAtm,3,JAtm) +
     $       sva(1,IAtm)*DenXA(3,JAtm) + svb(1,IAtm)*DenXB(3,JAtm) +
     $       sga(1,IAtm)*gxaa(3,JAtm)  + sgb(1,IAtm)*gxbb(3,JAtm) +
     $       sgc(1,IAtm)*gxab(3,JAtm)  + ra*DenXXA(1,IAtm,3,JAtm) +
     $       rb*DenXXB(1,IAtm,3,JAtm)  + ga*gxxaa(1,IAtm,3,JAtm) +
     $       gb*gxxbb(1,IAtm,3,JAtm)   + gc*gxxab(1,IAtm,3,JAtm)
            hess(2,IAtm,3,JAtm) = hess(2,IAtm,3,JAtm) +
     $       sva(2,IAtm)*DenXA(3,JAtm) + svb(2,IAtm)*DenXB(3,JAtm) +
     $       sga(2,IAtm)*gxaa(3,JAtm)  + sgb(2,IAtm)*gxbb(3,JAtm) +
     $       sgc(2,IAtm)*gxab(3,JAtm)  + ra*DenXXA(2,IAtm,3,JAtm) +
     $       rb*DenXXB(2,IAtm,3,JAtm)  + ga*gxxaa(2,IAtm,3,JAtm) +
     $       gb*gxxbb(2,IAtm,3,JAtm)   + gc*gxxab(2,IAtm,3,JAtm)
            hess(3,IAtm,3,JAtm) = hess(3,IAtm,3,JAtm) +
     $       sva(3,IAtm)*DenXA(3,JAtm) + svb(3,IAtm)*DenXB(3,JAtm) +
     $       sga(3,IAtm)*gxaa(3,JAtm)  + sgb(3,IAtm)*gxbb(3,JAtm) +
     $       sgc(3,IAtm)*gxab(3,JAtm)  + ra*DenXXA(3,IAtm,3,JAtm) +
     $       rb*DenXXB(3,IAtm,3,JAtm)  + ga*gxxaa(3,IAtm,3,JAtm) +
     $       gb*gxxbb(3,IAtm,3,JAtm)   + gc*gxxab(3,IAtm,3,JAtm)
          EndDo
        EndDo
c
c  --  add weight derivatives contribution.
c      the sva array is used to store a partial sum
c
        do ia=1,natoms
          sva(1,ia)=prap*denxa(1,ia)+prbp*denxb(1,ia)+pgap*gxaa(1,ia)+
     $               pgbp*gxbb(1,ia)+pgcp*gxab(1,ia)
          sva(2,ia)=prap*denxa(2,ia)+prbp*denxb(2,ia)+pgap*gxaa(2,ia)+
     $               pgbp*gxbb(2,ia)+pgcp*gxab(2,ia)
          sva(3,ia)=prap*denxa(3,ia)+prbp*denxb(3,ia)+pgap*gxaa(3,ia)+
     $               pgbp*gxbb(3,ia)+pgcp*gxab(3,ia)
        enddo
        do ia=1,natoms
        if(ia.ne.icntr)then
          do ja=ia,natoms
          if(ja.ne.icntr)then
          hess(1,ia,1,ja)=hess(1,ia,1,ja)+gwt(1,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(1,ia)+hwt(1,ia,1,ja,ipp)
          hess(2,ia,1,ja)=hess(2,ia,1,ja)+gwt(2,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(2,ia)+hwt(2,ia,1,ja,ipp)
          hess(3,ia,1,ja)=hess(3,ia,1,ja)+gwt(3,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(3,ia)+hwt(3,ia,1,ja,ipp)
          hess(1,ia,2,ja)=hess(1,ia,2,ja)+gwt(1,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(1,ia)+hwt(1,ia,2,ja,ipp)
          hess(2,ia,2,ja)=hess(2,ia,2,ja)+gwt(2,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(2,ia)+hwt(2,ia,2,ja,ipp)
          hess(3,ia,2,ja)=hess(3,ia,2,ja)+gwt(3,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(3,ia)+hwt(3,ia,2,ja,ipp)
          hess(1,ia,3,ja)=hess(1,ia,3,ja)+gwt(1,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(1,ia)+hwt(1,ia,3,ja,ipp)
          hess(2,ia,3,ja)=hess(2,ia,3,ja)+gwt(2,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(2,ia)+hwt(2,ia,3,ja,ipp)
          hess(3,ia,3,ja)=hess(3,ia,3,ja)+gwt(3,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(3,ia)+hwt(3,ia,3,ja,ipp)
          endif
          enddo
        endif
        enddo
        call secund(t6)
        tqh=tqh+t6-t5
 200    continue
      endif
c
      return
      end
c =====================================================================
      subroutine cwmpfh(dft,    npp,    nbas,   nbf,    natoms,
     $                  natonce,nb,     ne,     nbatm,  thrsh,
     $                  da,     dm,     gden,   wght,   pra,
     $                  prara,  prarb,  pga,    pgc,    praga,
     $                  pragb,  pragc,  pgaga,  pgagb,  pgagc,
     $                  pgcgc,  vao,    vaox,   vaoxx,  vaoxxx,
     $                  inb,    vm,     indx,   denx,   denxx,
     $                  gdx,    gdxx,   gx,     gxx,    sv,
     $                  sw,     icntr,  gwt,    hwt,    trp,
     $                  fda,    hess,   td1,    tg1,    tsw,
     $                  tqf,    tqh)
c     implicit real*8 (a-h,o-z)
      implicit none             ! trying to avoid typing errors
c
c
c  carries out numerical integration and accumulates contribution
c  from current grid point into derivative fock matrices and hessian
c  taking into account the weight derivatives.
c
c  natonce components of the derivative fock matrices are computed in
c  one pass (components nb to ne, with nb <= icntr <= ne). The
c  icntr component of the derivative fock matrix is computed
c  by applying translational invariance at each grid point.
c  The Hessian components involving icntr are not computed, as
c  they can be obtained by translational invariance at the end
c  of the calculation
c
c  ** closed shell **
c
c  arguments
c
c  dft     -  method flag (note: all methods include slater exchange)
c              1 - 3 - local correlation
c             >3 - nonlocal
c  npp     -  number of contributing (non-zero) grid points this batch
c  nbas    -  total number of basis functions
c  nbf     -  indexing array to location of "non-zero" basis functions
c  natoms  -  number of atoms
c  natonce -  number of atomic derivatives to compute in this pass
c  nb,ne   -  first and last deivative to compute (natonce=ne-nb+1)
c  nbatm   -  basis functions --> atom index
c  thrsh   -  threshold for neglect of contribution
c  da      -  closed-shell density matrix (lower triangle)
c  dm      -  maximum density matrix element per column
c  gden    -  density gradient at grid points (non local only, dft > 3)
c  wght    -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density (dft > 3)
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
c  vao     -  "non-zero" basis function values at grid points
c  potx    -  gradient functional derivative (dft > 3 only)
c  potxx   -  gradient functional second derivative (dft > 3)
c  potxd   -  density-gradient mixed second derivative (dft > 3)
c  vaox    -  basis function 1st derivatives at grid points
c  vaoxx   -  basis function 2nd derivatives at grid points
c  vaoxxx  -  basis function 3rd derivatives at grid points (dft > 3)
c  inb     -  indexing array for non-zero entries to vao
c  vm      -  array containing maximum magnitude ao per grid point
c  indx    -  index into contributing columns of vao
c  denx    -  scratch storage for density gradient per grid point
c  denxx   -    ditto for density hessian per grid point
c  gdx     -  ditto for atomic gradient of density gradient (dft > 3)
c  gdxx    -  ditto for atomic hessian of density gradient (dft > 3)
c  gx      -  ditto for atomic gradient of gradient invariant (dft > 3)
c  gxx     -  ditto for atomic hessian of gradient invariant (dft > 3)
c  sv      -  storage for coefficient of vao(i)*vao(j) in quadrature
c          -  formula (dft > 3)
c  sw      -  storage for coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  icntr   -  current atomic center
c  gwt     -  gradient of quadrature weight
c  hwt     -  hessian of quadrature weight
c  trp     -  scratch area for translational invariance
c
c  on exit
c
c  fda     -  contribution to derivative fock matrices
c  hess    -  contribution to hessian matrix
c
c
      integer npp,nbas,natoms,icntr,natonce,nb,ne
      real*8 wght(*),gwt(3,natoms,*),hwt(3,natoms,3,natoms,*)
      real*8 vao(*),vaox(3,*),vaoxx(6,*),vaoxxx(10,*)
      real*8 vm(*),da(nbas*(nbas+1)/2),dm(nbas),gden(3,*)
      real*8 pra(npp),pga(npp),pgc(npp),prara(npp)
      real*8 prarb(npp),praga(npp),pragb(npp),pragc(npp)
      real*8 pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
      real*8 denx(3,natoms),denxx(3,natoms,3,natoms)
      real*8 gdx(3,3,natoms),gdxx(3,3,natoms,3,natoms)
      real*8 gx(3,natoms),gxx(3,natoms,3,natoms)
      real*8 sv(3,natoms),sw(3,3,natoms)
      integer dft,nbf(*),inb(*),indx(npp),nbatm(nbas)
      real*8 fda(3,nb:ne,*),hess(3,natoms,3,natoms)
      real*8 trp(3,natoms)
      real*8 thrsh
      real*8 td1,tg1,tsw,tqf,tqh
c
      real*8 zero,two,three,four,six,twelve,half
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0,half=0.5d0)
      real*8 epsi,alot
      parameter (epsi=1.0d-15,alot=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,doval
      integer nat3,ip,ipp,i,j,ii,jj,ij,it
      integer ia,ja,iatm,jatm,iixx,katm,k1,isv,isw
      real*8 ra,ra2,ga,gc,rara,rara2,raga,ragb,ragc,prap
      real*8 gaga,gagb,gagc,gcgc,vmx,vmx2,wg,prg,prg2,pgg,pgg2,pg,pg2
      real*8 abra,abrara,thrx,thrx1,thrxx,valt,abdaij,abdaijm,abval
      real*8 valmm,valmm1,valmm2,valxx,valxy,valyy,valxz,valyz,valzz
      real*8 abvj,abvvj,valj,valjx,valjy,valjz,abijt,abij
      real*8 thrsh1,thtest
      real*8 t1,t2,t3,t4,t5,t6
      real*8 gwtmax,dmx,valx,valy,valz,daij,xyc,xzc,yzc
      real*8 dmaxyz,vmax,val,vali,vdxc,valm,vdjx,xc
      real*8 hvalx,hvaly,hvalz,dx,dy,dz,sxx,sxy,sxz
      real*8 vd1,vd2,vd3,svmax,swmax,tswmax
      real*8 vmax1,vmax2,vmax3
      real*8 valix,valiy,valiz,xvi,xxvx,xxvy,xxvz,valtest,valmxj
      real*8 valmx,valmx1,xvj,valm1,val1,smaxmx,val2,val3,val4
      real*8 vij,vijx,vijy,vijz,gij,hdx,hdy,hdz,hgx,hgy,hgz
      real*8 daijt,smax,dmax,potp,potxp,potp2,potxp2
      real*8 sxxi,sxyi,sxzi,dpmax2
      real*8 valt1,valt2
c
c consistency check
c
      if(icntr.lt.nb.or.icntr.gt.ne)then
        call nerror(1,'cwmphf','icntr is wrong',icntr,nb)
      endif
      nat3 = 3*natoms
c
      if(dft.le.3) then
c
c  local only
c
      do 50 ip=1,npp
      ipp = indx(ip)
      wg  = wght(ipp)
      prap=pra(ip)
      ra = prap*wg
      rara = prara(ip)*wg
      vmx = vm(ipp)
c
c   initializations for weight derivatives
c
      gwtmax=zero
      do ia=1,natoms
        if(ia.ne.icntr)then
          gwt(1,ia,ipp)=gwt(1,ia,ipp)*prap
          gwt(2,ia,ipp)=gwt(2,ia,ipp)*prap
          gwt(3,ia,ipp)=gwt(3,ia,ipp)*prap
          gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +    abs(gwt(3,ia,ipp)))
        endif
      enddo
      call zeroit(denx,nat3)
      call zeroit(denxx,nat3**2)
c
c    form gradient and hessian density at current grid point.
c
c    note:
c    for the closed shell case, the array da contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx and denxx will contain half the gradient and hessian
c    of the total closed shell density
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrara=abs(rara)
      dodenx=abrara.gt.epsi
      thrx=alot
      if(dodenx)thrx=thrsh/abrara
      thrx1=thrx/vmx2
      abra=abs(ra)
      dodenxx=abra.gt.epsi
      thrxx=alot
      if(dodenxx)thrxx=thrsh/abra
      if(.not.(dodenx.or.dodenxx))goto 25 ! should never happen!
c
c   now compute density derivatives
c
      do 20 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        dmx = dm(ii)
        if(dmx*vmx**2.lt.thrsh) go to 20
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        if(iatm.eq.icntr)goto 20
        valx = vaox(1,i)
        valy = vaox(2,i)
        valz = vaox(3,i)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx*vmx
        dox=(dodenx.and.(valt.gt.thrx1.or.valt**2.gt.thrx))
        doxx=(dodenxx.and.dmx*vmx2.gt.thrxx)
        if(.not.(dox.or.doxx)) goto 20
        do 18 j=nbf(ipp)+1,i
          jj = inb(j)
          ij = it + jj
          jatm = nbatm(jj)
          daij = da(ij)*vao(j)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denx(1,iatm) = denx(1,iatm) - daij*valx
              denx(2,iatm) = denx(2,iatm) - daij*valy
              denx(3,iatm) = denx(3,iatm) - daij*valz
          endif
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
          if(doxx.and.abdaijm.gt.thrxx)then
              denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*vaoxx(1,i)
              xyc=daij*vaoxx(2,i)
              denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
              denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
              denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*vaoxx(3,i)
              xzc=daij*vaoxx(4,i)
              denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
              denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
              yzc=daij*vaoxx(5,i)
              denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
              denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
              denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*vaoxx(6,i)
          endif
          if(jatm.lt.iatm) go to 18
          if(jatm.eq.icntr)goto 18
          abdaijm=abs(da(ij))*valm*vmx
          if(doxx.and.abdaijm.gt.thrxx)then
c
c    (b) iatm with jatm
c
            daij = da(ij)*vaox(1,j)
              denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
              denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
              denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
              daij = da(ij)*vaox(2,j)
              denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
              denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
              denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
              daij = da(ij)*vaox(3,j)
              denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
              denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
              denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
          endif
 18     continue
        do 19 j=i+1,nbf(ipp+1)
          jj = inb(j)
          ij = (jj*(jj-1))/2 + ii
          jatm = nbatm(jj)
          daij = da(ij)*vao(j)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denx(1,iatm) = denx(1,iatm) - daij*valx
              denx(2,iatm) = denx(2,iatm) - daij*valy
              denx(3,iatm) = denx(3,iatm) - daij*valz
          endif
          if(doxx.and.abdaijm.gt.thrxx)then
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
              denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*vaoxx(1,i)
              xyc=daij*vaoxx(2,i)
              denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
              denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
              denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*vaoxx(3,i)
              xzc=daij*vaoxx(4,i)
              denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
              denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
              yzc=daij*vaoxx(5,i)
              denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
              denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
              denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*vaoxx(6,i)
          endif
c
c    (b) iatm with jatm
c
          if(jatm.lt.iatm) go to 19
          if(jatm.eq.icntr)goto 19
          abdaijm=abs(da(ij))*valm*vmx
          if(doxx.and.abdaijm.gt.thrxx)then
              daij = da(ij)*vaox(1,j)
              denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
              denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
              denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
              daij = da(ij)*vaox(2,j)
              denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
              denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
              denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
              daij = da(ij)*vaox(3,j)
              denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
              denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
              denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
          endif
 19     continue
 20   continue
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c    numerical quadrature for the derivative fock matrix.
c
c -- get maximum absolute value of 1st-order density at this grid point
c
      call absmax(nat3,denx,iixx,dmaxyz)
c
c -- global threshold testing
c
      vmax = max(abra*vmx2,abrara*dmaxyz*vmx2)
      if(vmax.lt.thrsh) go to 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      do 40 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        vali = vao(i)
        val = vali*ra
        vdxc = vali*rara
        valx = vaox(1,i)*ra
        valy = vaox(2,i)*ra
        valz = vaox(3,i)*ra
        abval= abs(val)
        doval= abval.gt.thrsh1
        valm = max(abs(valx),abs(valy),abs(valz))
        valt = max(abval,valm,abs(vdxc)*dmaxyz)
        if(valt.gt.thrsh1) then
          do 30 j=nbf(ipp)+1,i
            call zeroit(trp,nat3)
            jj = inb(j)
            ij = it + jj
            jatm = nbatm(jj)
            valj = vao(j)
            if(abs(valj)*valm.gt.thrsh)then
              if(iatm.ne.icntr) then
                trp(1,iatm) = trp(1,iatm) - valx*valj
                trp(2,iatm) = trp(2,iatm) - valy*valj
                trp(3,iatm) = trp(3,iatm) - valz*valj
              endif
            endif
            if(doval)then
              if(jatm.ne.icntr) then
                trp(1,jatm) = trp(1,jatm) - val*vaox(1,j)
                trp(2,jatm) = trp(2,jatm) - val*vaox(2,j)
                trp(3,jatm) = trp(3,jatm) - val*vaox(3,j)
              endif
            endif
c
c -- contribution of functional derivative
c
            vdjx = vdxc*valj
            if(abs(vdjx)*dmaxyz.gt.thrsh) then
              do katm=1,natoms
                trp(1,katm)=trp(1,katm) + vdjx*denx(1,katm)
                trp(2,katm)=trp(2,katm) + vdjx*denx(2,katm)
                trp(3,katm)=trp(3,katm) + vdjx*denx(3,katm)
              enddo
            endif
c
c -- weight derivatives contribution
c
            xc=vali*valj
            do ia=1,natoms
              if(ia.ne.icntr)then
                trp(1,ia)=trp(1,ia)+xc*gwt(1,ia,ipp)
                trp(2,ia)=trp(2,ia)+xc*gwt(2,ia,ipp)
                trp(3,ia)=trp(3,ia)+xc*gwt(3,ia,ipp)
              endif
            enddo
c
c -- apply translational invariance and accumulate contributions
c    into the fda array
c
              do ia=1,natoms
                if(ia.ne.icntr)then
                  trp(1,icntr)=trp(1,icntr)-trp(1,ia)
                  trp(2,icntr)=trp(2,icntr)-trp(2,ia)
                  trp(3,icntr)=trp(3,icntr)-trp(3,ia)
                endif
              enddo
              do ia=nb,ne
                fda(1,ia,ij)=fda(1,ia,ij)+trp(1,ia)
                fda(2,ia,ij)=fda(2,ia,ij)+trp(2,ia)
                fda(3,ia,ij)=fda(3,ia,ij)+trp(3,ia)
            EndDO
 30       continue
        endif
 40   continue
c
c    numerical quadrature for the hessian matrix
c
 45   continue
      call secund(t3)
      tqf=tqf+t3-t2
c
c -- direct contribution to hessian matrix.
c    a factor two is applied
c
      rara2=two*rara
      ra2=two*ra
      do iatm=1,natoms
        hvalx=rara2*denx(1,iatm)
        hvaly=rara2*denx(2,iatm)
        hvalz=rara2*denx(3,iatm)
        do jatm=iatm,natoms
          hess(1,iatm,1,jatm) = hess(1,iatm,1,jatm) +
     $     denx(1,jatm)*hvalx + ra2*denxx(1,iatm,1,jatm)
          hess(2,iatm,1,jatm) = hess(2,iatm,1,jatm) +
     $     denx(1,jatm)*hvaly + ra2*denxx(2,iatm,1,jatm)
          hess(3,iatm,1,jatm) = hess(3,iatm,1,jatm) +
     $     denx(1,jatm)*hvalz + ra2*denxx(3,iatm,1,jatm)
          hess(1,iatm,2,jatm) = hess(1,iatm,2,jatm) +
     $     denx(2,jatm)*hvalx + ra2*denxx(1,iatm,2,jatm)
          hess(2,iatm,2,jatm) = hess(2,iatm,2,jatm) +
     $     denx(2,jatm)*hvaly + ra2*denxx(2,iatm,2,jatm)
          hess(3,iatm,2,jatm) = hess(3,iatm,2,jatm) +
     $     denx(2,jatm)*hvalz + ra2*denxx(3,iatm,2,jatm)
          hess(1,iatm,3,jatm) = hess(1,iatm,3,jatm) +
     $     denx(3,jatm)*hvalx + ra2*denxx(1,iatm,3,jatm)
          hess(2,iatm,3,jatm) = hess(2,iatm,3,jatm) +
     $     denx(3,jatm)*hvaly + ra2*denxx(2,iatm,3,jatm)
          hess(3,iatm,3,jatm) = hess(3,iatm,3,jatm) +
     $     denx(3,jatm)*hvalz + ra2*denxx(3,iatm,3,jatm)
        enddo
      enddo
c
c  --  add weight derivatives contribution
c
      do ia=1,natoms
      if(ia.ne.icntr)then
        do ja=ia,natoms
        if(ja.ne.icntr)then
        hess(1,ia,1,ja)=hess(1,ia,1,ja) + two*(gwt(1,ia,ipp)*denx(1,ja)
     $   + gwt(1,ja,ipp)*denx(1,ia)) + hwt(1,ia,1,ja,ipp)
        hess(2,ia,1,ja)=hess(2,ia,1,ja) + two*(gwt(2,ia,ipp)*denx(1,ja)
     $   + gwt(1,ja,ipp)*denx(2,ia)) + hwt(2,ia,1,ja,ipp)
        hess(3,ia,1,ja)=hess(3,ia,1,ja) + two*(gwt(3,ia,ipp)*denx(1,ja)
     $   + gwt(1,ja,ipp)*denx(3,ia)) + hwt(3,ia,1,ja,ipp)
        hess(1,ia,2,ja)=hess(1,ia,2,ja) + two*(gwt(1,ia,ipp)*denx(2,ja)
     $   + gwt(2,ja,ipp)*denx(1,ia)) + hwt(1,ia,2,ja,ipp)
        hess(2,ia,2,ja)=hess(2,ia,2,ja) + two*(gwt(2,ia,ipp)*denx(2,ja)
     $   + gwt(2,ja,ipp)*denx(2,ia)) + hwt(2,ia,2,ja,ipp)
        hess(3,ia,2,ja)=hess(3,ia,2,ja) + two*(gwt(3,ia,ipp)*denx(2,ja)
     $   + gwt(2,ja,ipp)*denx(3,ia)) + hwt(3,ia,2,ja,ipp)
        hess(1,ia,3,ja)=hess(1,ia,3,ja) + two*(gwt(1,ia,ipp)*denx(3,ja)
     $   + gwt(3,ja,ipp)*denx(1,ia)) + hwt(1,ia,3,ja,ipp)
        hess(2,ia,3,ja)=hess(2,ia,3,ja) + two*(gwt(2,ia,ipp)*denx(3,ja)
     $   + gwt(3,ja,ipp)*denx(2,ia)) + hwt(2,ia,3,ja,ipp)
        hess(3,ia,3,ja)=hess(3,ia,3,ja) + two*(gwt(3,ia,ipp)*denx(3,ja)
     $   + gwt(3,ja,ipp)*denx(3,ia)) + hwt(3,ia,3,ja,ipp)
        endif
        enddo
      endif
      enddo
      call secund(t4)
      tqh=tqh+t4-t3
c
c   local dft ends here
c
 50   continue
cc
      else
cc
c
c   non-local dft
c
      do 200 ip=1,npp
        ipp = indx(ip)
        wg = wght(ipp)
        vmx = vm(ipp)
        ra = pra(ip)*wg
        ga = pga(ip)*wg
        gc = pgc(ip)*wg
        rara = (prara(ip)+prarb(ip))*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  some sums of the potentials that will be used later
c
        prg=raga+ragb+ragc
        pgg=gaga+gagb+two*gagc+half*gcgc
        pg=ga+half*gc
c
c  potential combinations for weight derivatives contributions
c
        potp = pra(ip)
        potxp = pga(ip)+half*pgc(ip)
        potp2 = potp+potp
        potxp2 = potxp+potxp
c
c  density gradient at current point
c
        dx=gden(1,ipp)
        dy=gden(2,ipp)
        dz=gden(3,ipp)
c
c  zero out derivatives of densities and gradients
c
        call zeroit(denx,nat3)
        call zeroit(denxx,nat3**2)
        call zeroit(gdx,3*nat3)
        call zeroit(gdxx,3*nat3**2)
        call zeroit(gx,nat3)
        call zeroit(gxx,nat3**2)
c
c   initializations for weight derivatives
c
c   unlike the local case above, we do not multiply
c   the weight gradient by the potential at this stage
c
        gwtmax=zero
        do ia=1,natoms
          if(ia.ne.icntr)then
            gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +        abs(gwt(3,ia,ipp)))
          endif
        enddo
c
c    form the atomic gradient and hessian  of density
c    and density gradient at current grid point
c
c    note:
c    for the closed shell case, the array da contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx and denxx will contain half the atomic gradient and hessian
c    of the total closed shell density. likevise, gdx and gdxx will
c    contain half the atomic gradient and hessian of the total closed
c    shell density gradient
c
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c    in addition, we just check whether the second derivatives
c    of density and density gradient have to be computed at all.
c
        call secund(t1)
        abra=abs(ra)
        dodenxx=abra.gt.epsi
        do 220 i=nbf(ipp)+1,nbf(ipp+1)
          ii = inb(i)
          dmx = dm(ii)       ! max. element of first order density
          if(dmx.gt.epsi)then
             thtest=thrsh/(vmx*dmx)
          else
            goto 220
          endif
          it = (ii*(ii-1))/2
          iatm = nbatm(ii)
          if(iatm.eq.icntr) goto 220
          valx = vaox(1,i)
          valy = vaox(2,i)
          valz = vaox(3,i)
          valm = max(abs(valx),abs(valy),abs(valz))
          valxx = vaoxx(1,i)
          valxy = vaoxx(2,i)
          valyy = vaoxx(3,i)
          valxz = vaoxx(4,i)
          valyz = vaoxx(5,i)
          valzz = vaoxx(6,i)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          do 218 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            valj=vao(j)
            jatm = nbatm(jj)
            abvj=abs(valj)
            valjx=vaox(1,j)
            valjy=vaox(2,j)
            valjz=vaox(3,j)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            daijt=da(ij)
            abijt=abs(daijt)
            daij = daijt*valj
            abij=abs(daij)
c
c -- atomic gradient of density
c
              if(abij*valm.gt.thrsh)then
                denx(1,iatm) = denx(1,iatm) - daij*valx
                denx(2,iatm) = denx(2,iatm) - daij*valy
                denx(3,iatm) = denx(3,iatm) - daij*valz
              endif
c
c -- atomic gradient of density gradient
c
             if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              gdx(1,1,iatm)=gdx(1,1,iatm)-daijt*(valj*valxx+valx*valjx)
              gdx(1,2,iatm)=gdx(1,2,iatm)-daijt*(valj*valxy+valy*valjx)
              gdx(1,3,iatm)=gdx(1,3,iatm)-daijt*(valj*valxz+valz*valjx)
              gdx(2,1,iatm)=gdx(2,1,iatm)-daijt*(valj*valxy+valx*valjy)
              gdx(2,2,iatm)=gdx(2,2,iatm)-daijt*(valj*valyy+valy*valjy)
              gdx(2,3,iatm)=gdx(2,3,iatm)-daijt*(valj*valyz+valz*valjy)
              gdx(3,1,iatm)=gdx(3,1,iatm)-daijt*(valj*valxz+valx*valjz)
              gdx(3,2,iatm)=gdx(3,2,iatm)-daijt*(valj*valyz+valy*valjz)
              gdx(3,3,iatm)=gdx(3,3,iatm)-daijt*(valj*valzz+valz*valjz)
             endif
c
c -- (a) one center terms: iatm with iatm
c
c
c -- atomic hessian of density
c
              if(dodenxx.and.abij*valmm.gt.thrsh)then
                denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*valxx
                xyc=daij*valxy
                denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
                denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
                denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*valyy
                xzc=daij*valxz
                denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
                denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
                yzc=daij*valyz
                denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
                denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
                denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*valzz
              endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
                gdxx(1,1,iatm,1,iatm)=gdxx(1,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(1,i)*valj+valxx*valjx)
                gdxx(1,1,iatm,2,iatm)=gdxx(1,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,1,iatm)=gdxx(1,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,2,iatm)=gdxx(1,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valyy*valjx)
                gdxx(1,1,iatm,3,iatm)=gdxx(1,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,3,iatm,1,iatm)=gdxx(1,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,2,iatm,3,iatm)=gdxx(1,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,2,iatm)=gdxx(1,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,3,iatm)=gdxx(1,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valzz*valjx)
c
                gdxx(2,1,iatm,1,iatm)=gdxx(2,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxx*valjy)
                gdxx(2,1,iatm,2,iatm)=gdxx(2,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,1,iatm)=gdxx(2,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,2,iatm)=gdxx(2,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(4,i)*valj+valyy*valjy)
                gdxx(2,1,iatm,3,iatm)=gdxx(2,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,3,iatm,1,iatm)=gdxx(2,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,2,iatm,3,iatm)=gdxx(2,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,2,iatm)=gdxx(2,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,3,iatm)=gdxx(2,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valzz*valjy)

                gdxx(3,1,iatm,1,iatm)=gdxx(3,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxx*valjz)
                gdxx(3,1,iatm,2,iatm)=gdxx(3,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,1,iatm)=gdxx(3,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,2,iatm)=gdxx(3,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyy*valjz)
                gdxx(3,1,iatm,3,iatm)=gdxx(3,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,3,iatm,1,iatm)=gdxx(3,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,2,iatm,3,iatm)=gdxx(3,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,2,iatm)=gdxx(3,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,3,iatm)=gdxx(3,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(10,i)*valj+valzz*valjz)
              endif
c
c -- (b) two center terms: iatm with jatm
            if(jatm.lt.iatm) go to 218
            if(jatm.eq.icntr)goto 218
c
c -- atomic hessian of density
c
              if(dodenxx.and.abijt*abvvj*valm.gt.thrsh)then
                daij = daijt*valjx
                denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
                denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
                denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
                daij = daijt*valjy
                denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
                denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
                denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
                daij = daijt*valjz
                denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
                denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
                denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
              endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
                gdxx(1,1,iatm,1,jatm)=gdxx(1,1,iatm,1,jatm)+
     $                      daijt * (valxx*valjx+valx*vaoxx(1,j))
                gdxx(1,1,iatm,2,jatm)=gdxx(1,1,iatm,2,jatm)+
     $                      daijt * (valxx*valjy+valx*vaoxx(2,j))
                gdxx(1,2,iatm,1,jatm)=gdxx(1,2,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valy*vaoxx(1,j))
                gdxx(1,2,iatm,2,jatm)=gdxx(1,2,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valy*vaoxx(2,j))
                gdxx(1,1,iatm,3,jatm)=gdxx(1,1,iatm,3,jatm)+
     $                      daijt * (valxx*valjz+valx*vaoxx(4,j))
                gdxx(1,3,iatm,1,jatm)=gdxx(1,3,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valz*vaoxx(1,j))
                gdxx(1,2,iatm,3,jatm)=gdxx(1,2,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valy*vaoxx(4,j))
                gdxx(1,3,iatm,2,jatm)=gdxx(1,3,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valz*vaoxx(2,j))
                gdxx(1,3,iatm,3,jatm)=gdxx(1,3,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valz*vaoxx(4,j))
c
                gdxx(2,1,iatm,1,jatm)=gdxx(2,1,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valx*vaoxx(2,j))
                gdxx(2,1,iatm,2,jatm)=gdxx(2,1,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valx*vaoxx(3,j))
                gdxx(2,2,iatm,1,jatm)=gdxx(2,2,iatm,1,jatm)+
     $                      daijt * (valyy*valjx+valy*vaoxx(2,j))
                gdxx(2,2,iatm,2,jatm)=gdxx(2,2,iatm,2,jatm)+
     $                      daijt * (valyy*valjy+valy*vaoxx(3,j))
                gdxx(2,1,iatm,3,jatm)=gdxx(2,1,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valx*vaoxx(5,j))
                gdxx(2,3,iatm,1,jatm)=gdxx(2,3,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valz*vaoxx(2,j))
                gdxx(2,2,iatm,3,jatm)=gdxx(2,2,iatm,3,jatm)+
     $                      daijt * (valyy*valjz+valy*vaoxx(5,j))
                gdxx(2,3,iatm,2,jatm)=gdxx(2,3,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valz*vaoxx(3,j))
                gdxx(2,3,iatm,3,jatm)=gdxx(2,3,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valz*vaoxx(5,j))
c
                gdxx(3,1,iatm,1,jatm)=gdxx(3,1,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valx*vaoxx(4,j))
                gdxx(3,1,iatm,2,jatm)=gdxx(3,1,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valx*vaoxx(5,j))
                gdxx(3,2,iatm,1,jatm)=gdxx(3,2,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valy*vaoxx(4,j))
                gdxx(3,2,iatm,2,jatm)=gdxx(3,2,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valy*vaoxx(5,j))
                gdxx(3,1,iatm,3,jatm)=gdxx(3,1,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valx*vaoxx(6,j))
                gdxx(3,3,iatm,1,jatm)=gdxx(3,3,iatm,1,jatm)+
     $                      daijt * (valzz*valjx+valz*vaoxx(4,j))
                gdxx(3,2,iatm,3,jatm)=gdxx(3,2,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valy*vaoxx(6,j))
                gdxx(3,3,iatm,2,jatm)=gdxx(3,3,iatm,2,jatm)+
     $                      daijt * (valzz*valjy+valz*vaoxx(5,j))
                gdxx(3,3,iatm,3,jatm)=gdxx(3,3,iatm,3,jatm)+
     $                      daijt * (valzz*valjz+valz*vaoxx(6,j))
              endif
 218      continue
          do 219 j=i+1,nbf(ipp+1)
            jj = inb(j)
            ij = (jj*(jj-1))/2 + ii
            jatm = nbatm(jj)
            valj=vao(j)
            abvj=abs(valj)
            valjx=vaox(1,j)
            valjy=vaox(2,j)
            valjz=vaox(3,j)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            daijt= da(ij)
            abijt=abs(daijt)
            daij = daijt*valj
            abij=abs(daij)
c
c -- atomic gradient of density
c
              if(abij*valm.gt.thrsh)then
                denx(1,iatm) = denx(1,iatm) - daij*valx
                denx(2,iatm) = denx(2,iatm) - daij*valy
                denx(3,iatm) = denx(3,iatm) - daij*valz
              endif
c
c -- atomic gradient of density gradient
c
             if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              gdx(1,1,iatm)=gdx(1,1,iatm)-daijt*(valj*valxx+valx*valjx)
              gdx(1,2,iatm)=gdx(1,2,iatm)-daijt*(valj*valxy+valy*valjx)
              gdx(1,3,iatm)=gdx(1,3,iatm)-daijt*(valj*valxz+valz*valjx)
              gdx(2,1,iatm)=gdx(2,1,iatm)-daijt*(valj*valxy+valx*valjy)
              gdx(2,2,iatm)=gdx(2,2,iatm)-daijt*(valj*valyy+valy*valjy)
              gdx(2,3,iatm)=gdx(2,3,iatm)-daijt*(valj*valyz+valz*valjy)
              gdx(3,1,iatm)=gdx(3,1,iatm)-daijt*(valj*valxz+valx*valjz)
              gdx(3,2,iatm)=gdx(3,2,iatm)-daijt*(valj*valyz+valy*valjz)
              gdx(3,3,iatm)=gdx(3,3,iatm)-daijt*(valj*valzz+valz*valjz)
             endif
c
c -- (a) one center terms: iatm with iatm
c
c
c -- atomic hessian of density
c
            if(dodenxx.and.abij*valmm.gt.thrsh)then
              denxx(1,iatm,1,iatm)=denxx(1,iatm,1,iatm)+daij*valxx
              xyc=daij*valxy
              denxx(1,iatm,2,iatm)=denxx(1,iatm,2,iatm)+xyc
              denxx(2,iatm,1,iatm)=denxx(2,iatm,1,iatm)+xyc
              denxx(2,iatm,2,iatm)=denxx(2,iatm,2,iatm)+daij*valyy
              xzc=daij*valxz
              denxx(1,iatm,3,iatm)=denxx(1,iatm,3,iatm)+xzc
              denxx(3,iatm,1,iatm)=denxx(3,iatm,1,iatm)+xzc
              yzc=daij*valyz
              denxx(2,iatm,3,iatm)=denxx(2,iatm,3,iatm)+yzc
              denxx(3,iatm,2,iatm)=denxx(3,iatm,2,iatm)+yzc
              denxx(3,iatm,3,iatm)=denxx(3,iatm,3,iatm)+daij*valzz
            endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
                gdxx(1,1,iatm,1,iatm)=gdxx(1,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(1,i)*valj+valxx*valjx)
                gdxx(1,1,iatm,2,iatm)=gdxx(1,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,1,iatm)=gdxx(1,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxy*valjx)
                gdxx(1,2,iatm,2,iatm)=gdxx(1,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valyy*valjx)
                gdxx(1,1,iatm,3,iatm)=gdxx(1,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,3,iatm,1,iatm)=gdxx(1,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxz*valjx)
                gdxx(1,2,iatm,3,iatm)=gdxx(1,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,2,iatm)=gdxx(1,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valyz*valjx)
                gdxx(1,3,iatm,3,iatm)=gdxx(1,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valzz*valjx)
c
                gdxx(2,1,iatm,1,iatm)=gdxx(2,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(2,i)*valj+valxx*valjy)
                gdxx(2,1,iatm,2,iatm)=gdxx(2,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,1,iatm)=gdxx(2,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(3,i)*valj+valxy*valjy)
                gdxx(2,2,iatm,2,iatm)=gdxx(2,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(4,i)*valj+valyy*valjy)
                gdxx(2,1,iatm,3,iatm)=gdxx(2,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,3,iatm,1,iatm)=gdxx(2,3,iatm,1,iatm)+
     $                 daijt * (vaoxxx(6,i)*valj+valxz*valjy)
                gdxx(2,2,iatm,3,iatm)=gdxx(2,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,2,iatm)=gdxx(2,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyz*valjy)
                gdxx(2,3,iatm,3,iatm)=gdxx(2,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valzz*valjy)
c
                gdxx(3,1,iatm,1,iatm)=gdxx(3,1,iatm,1,iatm)+
     $                   daijt * (vaoxxx(5,i)*valj+valxx*valjz)
                gdxx(3,1,iatm,2,iatm)=gdxx(3,1,iatm,2,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,1,iatm)=gdxx(3,2,iatm,1,iatm)+
     $                   daijt * (vaoxxx(6,i)*valj+valxy*valjz)
                gdxx(3,2,iatm,2,iatm)=gdxx(3,2,iatm,2,iatm)+
     $                   daijt * (vaoxxx(7,i)*valj+valyy*valjz)
                gdxx(3,1,iatm,3,iatm)=gdxx(3,1,iatm,3,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,3,iatm,1,iatm)=gdxx(3,3,iatm,1,iatm)+
     $                   daijt * (vaoxxx(8,i)*valj+valxz*valjz)
                gdxx(3,2,iatm,3,iatm)=gdxx(3,2,iatm,3,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,2,iatm)=gdxx(3,3,iatm,2,iatm)+
     $                   daijt * (vaoxxx(9,i)*valj+valyz*valjz)
                gdxx(3,3,iatm,3,iatm)=gdxx(3,3,iatm,3,iatm)+
     $                   daijt * (vaoxxx(10,i)*valj+valzz*valjz)
              endif
c
c -- (b) two center terms: iatm with jatm
c
            if(jatm.lt.iatm) go to 219
            if(jatm.eq.icntr)goto 219
c
c -- atomic hessian of density
c
              if(dodenxx.and.abijt*abvvj*valm.gt.thrsh)then
                daij = daijt*valjx
                denxx(1,iatm,1,jatm)=denxx(1,iatm,1,jatm)+daij*valx
                denxx(2,iatm,1,jatm)=denxx(2,iatm,1,jatm)+daij*valy
                denxx(3,iatm,1,jatm)=denxx(3,iatm,1,jatm)+daij*valz
                daij = daijt*valjy
                denxx(1,iatm,2,jatm)=denxx(1,iatm,2,jatm)+daij*valx
                denxx(2,iatm,2,jatm)=denxx(2,iatm,2,jatm)+daij*valy
                denxx(3,iatm,2,jatm)=denxx(3,iatm,2,jatm)+daij*valz
                daij = daijt*valjz
                denxx(1,iatm,3,jatm)=denxx(1,iatm,3,jatm)+daij*valx
                denxx(2,iatm,3,jatm)=denxx(2,iatm,3,jatm)+daij*valy
                denxx(3,iatm,3,jatm)=denxx(3,iatm,3,jatm)+daij*valz
              endif
c
c -- atomic hessian of density gradient
c
              if(abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
                gdxx(1,1,iatm,1,jatm)=gdxx(1,1,iatm,1,jatm)+
     $                      daijt * (valxx*valjx+valx*vaoxx(1,j))
                gdxx(1,1,iatm,2,jatm)=gdxx(1,1,iatm,2,jatm)+
     $                      daijt * (valxx*valjy+valx*vaoxx(2,j))
                gdxx(1,2,iatm,1,jatm)=gdxx(1,2,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valy*vaoxx(1,j))
                gdxx(1,2,iatm,2,jatm)=gdxx(1,2,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valy*vaoxx(2,j))
                gdxx(1,1,iatm,3,jatm)=gdxx(1,1,iatm,3,jatm)+
     $                      daijt * (valxx*valjz+valx*vaoxx(4,j))
                gdxx(1,3,iatm,1,jatm)=gdxx(1,3,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valz*vaoxx(1,j))
                gdxx(1,2,iatm,3,jatm)=gdxx(1,2,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valy*vaoxx(4,j))
                gdxx(1,3,iatm,2,jatm)=gdxx(1,3,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valz*vaoxx(2,j))
                gdxx(1,3,iatm,3,jatm)=gdxx(1,3,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valz*vaoxx(4,j))
c
                gdxx(2,1,iatm,1,jatm)=gdxx(2,1,iatm,1,jatm)+
     $                      daijt * (valxy*valjx+valx*vaoxx(2,j))
                gdxx(2,1,iatm,2,jatm)=gdxx(2,1,iatm,2,jatm)+
     $                      daijt * (valxy*valjy+valx*vaoxx(3,j))
                gdxx(2,2,iatm,1,jatm)=gdxx(2,2,iatm,1,jatm)+
     $                      daijt * (valyy*valjx+valy*vaoxx(2,j))
                gdxx(2,2,iatm,2,jatm)=gdxx(2,2,iatm,2,jatm)+
     $                      daijt * (valyy*valjy+valy*vaoxx(3,j))
                gdxx(2,1,iatm,3,jatm)=gdxx(2,1,iatm,3,jatm)+
     $                      daijt * (valxy*valjz+valx*vaoxx(5,j))
                gdxx(2,3,iatm,1,jatm)=gdxx(2,3,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valz*vaoxx(2,j))
                gdxx(2,2,iatm,3,jatm)=gdxx(2,2,iatm,3,jatm)+
     $                      daijt * (valyy*valjz+valy*vaoxx(5,j))
                gdxx(2,3,iatm,2,jatm)=gdxx(2,3,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valz*vaoxx(3,j))
                gdxx(2,3,iatm,3,jatm)=gdxx(2,3,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valz*vaoxx(5,j))
c
                gdxx(3,1,iatm,1,jatm)=gdxx(3,1,iatm,1,jatm)+
     $                      daijt * (valxz*valjx+valx*vaoxx(4,j))
                gdxx(3,1,iatm,2,jatm)=gdxx(3,1,iatm,2,jatm)+
     $                      daijt * (valxz*valjy+valx*vaoxx(5,j))
                gdxx(3,2,iatm,1,jatm)=gdxx(3,2,iatm,1,jatm)+
     $                      daijt * (valyz*valjx+valy*vaoxx(4,j))
                gdxx(3,2,iatm,2,jatm)=gdxx(3,2,iatm,2,jatm)+
     $                      daijt * (valyz*valjy+valy*vaoxx(5,j))
                gdxx(3,1,iatm,3,jatm)=gdxx(3,1,iatm,3,jatm)+
     $                      daijt * (valxz*valjz+valx*vaoxx(6,j))
                gdxx(3,3,iatm,1,jatm)=gdxx(3,3,iatm,1,jatm)+
     $                      daijt * (valzz*valjx+valz*vaoxx(4,j))
                gdxx(3,2,iatm,3,jatm)=gdxx(3,2,iatm,3,jatm)+
     $                      daijt * (valyz*valjz+valy*vaoxx(6,j))
                gdxx(3,3,iatm,2,jatm)=gdxx(3,3,iatm,2,jatm)+
     $                      daijt * (valzz*valjy+valz*vaoxx(5,j))
                gdxx(3,3,iatm,3,jatm)=gdxx(3,3,iatm,3,jatm)+
     $                      daijt * (valzz*valjz+valz*vaoxx(6,j))
              endif
 219      continue
 220    continue
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and hessian of density gradient invariant
c
        do iatm=1,natoms
          gx(1,iatm)=two*(dx*gdx(1,1,iatm)+
     $                   dy*gdx(2,1,iatm)+dz*gdx(3,1,iatm))
          gx(2,iatm)=two*(dx*gdx(1,2,iatm)+
     $                   dy*gdx(2,2,iatm)+dz*gdx(3,2,iatm))
          gx(3,iatm)=two*(dx*gdx(1,3,iatm)+
     $                   dy*gdx(2,3,iatm)+dz*gdx(3,3,iatm))
        enddo
c
        do iatm=1,natoms
          do jatm=iatm,natoms
            gxx(1,iatm,1,jatm)=two*(
     $        dx*gdxx(1,1,iatm,1,jatm)+gdx(1,1,iatm)*gdx(1,1,jatm)+
     $        dy*gdxx(2,1,iatm,1,jatm)+gdx(2,1,iatm)*gdx(2,1,jatm)+
     $        dz*gdxx(3,1,iatm,1,jatm)+gdx(3,1,iatm)*gdx(3,1,jatm))
            gxx(1,iatm,2,jatm)=two*(
     $        dx*gdxx(1,1,iatm,2,jatm)+gdx(1,1,iatm)*gdx(1,2,jatm)+
     $        dy*gdxx(2,1,iatm,2,jatm)+gdx(2,1,iatm)*gdx(2,2,jatm)+
     $        dz*gdxx(3,1,iatm,2,jatm)+gdx(3,1,iatm)*gdx(3,2,jatm))
            gxx(2,iatm,1,jatm)=two*(
     $        dx*gdxx(1,2,iatm,1,jatm)+gdx(1,2,iatm)*gdx(1,1,jatm)+
     $        dy*gdxx(2,2,iatm,1,jatm)+gdx(2,2,iatm)*gdx(2,1,jatm)+
     $        dz*gdxx(3,2,iatm,1,jatm)+gdx(3,2,iatm)*gdx(3,1,jatm))
            gxx(2,iatm,2,jatm)=two*(
     $        dx*gdxx(1,2,iatm,2,jatm)+gdx(1,2,iatm)*gdx(1,2,jatm)+
     $        dy*gdxx(2,2,iatm,2,jatm)+gdx(2,2,iatm)*gdx(2,2,jatm)+
     $        dz*gdxx(3,2,iatm,2,jatm)+gdx(3,2,iatm)*gdx(3,2,jatm))
            gxx(1,iatm,3,jatm)=two*(
     $        dx*gdxx(1,1,iatm,3,jatm)+gdx(1,1,iatm)*gdx(1,3,jatm)+
     $        dy*gdxx(2,1,iatm,3,jatm)+gdx(2,1,iatm)*gdx(2,3,jatm)+
     $        dz*gdxx(3,1,iatm,3,jatm)+gdx(3,1,iatm)*gdx(3,3,jatm))
            gxx(3,iatm,1,jatm)=two*(
     $        dx*gdxx(1,3,iatm,1,jatm)+gdx(1,3,iatm)*gdx(1,1,jatm)+
     $        dy*gdxx(2,3,iatm,1,jatm)+gdx(2,3,iatm)*gdx(2,1,jatm)+
     $        dz*gdxx(3,3,iatm,1,jatm)+gdx(3,3,iatm)*gdx(3,1,jatm))
            gxx(2,iatm,3,jatm)=two*(
     $        dx*gdxx(1,2,iatm,3,jatm)+gdx(1,2,iatm)*gdx(1,3,jatm)+
     $        dy*gdxx(2,2,iatm,3,jatm)+gdx(2,2,iatm)*gdx(2,3,jatm)+
     $        dz*gdxx(3,2,iatm,3,jatm)+gdx(3,2,iatm)*gdx(3,3,jatm))
            gxx(3,iatm,2,jatm)=two*(
     $        dx*gdxx(1,3,iatm,2,jatm)+gdx(1,3,iatm)*gdx(1,2,jatm)+
     $        dy*gdxx(2,3,iatm,2,jatm)+gdx(2,3,iatm)*gdx(2,2,jatm)+
     $        dz*gdxx(3,3,iatm,2,jatm)+gdx(3,3,iatm)*gdx(3,2,jatm))
            gxx(3,iatm,3,jatm)=two*(
     $        dx*gdxx(1,3,iatm,3,jatm)+gdx(1,3,iatm)*gdx(1,3,jatm)+
     $        dy*gdxx(2,3,iatm,3,jatm)+gdx(2,3,iatm)*gdx(2,3,jatm)+
     $        dz*gdxx(3,3,iatm,3,jatm)+gdx(3,3,iatm)*gdx(3,3,jatm))
          enddo
        enddo
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the fock matrix
c    (i.e., quantities v, w and x at page 7436 of johnson and frisch).
c
        prg2=two*prg
        pgg2=two*pgg
        pg2=two*pg
        sxx=pg2*dx
        sxy=pg2*dy
        sxz=pg2*dz
        do iatm=1,natoms
          vd1=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          vd2=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          vd3=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          sv(1,iatm)=rara*denx(1,iatm)+prg*gx(1,iatm)
          sv(2,iatm)=rara*denx(2,iatm)+prg*gx(2,iatm)
          sv(3,iatm)=rara*denx(3,iatm)+prg*gx(3,iatm)
          sw(1,1,iatm)=pg2*gdx(1,1,iatm)+dx*vd1
          sw(2,1,iatm)=pg2*gdx(2,1,iatm)+dy*vd1
          sw(3,1,iatm)=pg2*gdx(3,1,iatm)+dz*vd1
          sw(1,2,iatm)=pg2*gdx(1,2,iatm)+dx*vd2
          sw(2,2,iatm)=pg2*gdx(2,2,iatm)+dy*vd2
          sw(3,2,iatm)=pg2*gdx(3,2,iatm)+dz*vd2
          sw(1,3,iatm)=pg2*gdx(1,3,iatm)+dx*vd3
          sw(2,3,iatm)=pg2*gdx(2,3,iatm)+dy*vd3
          sw(3,3,iatm)=pg2*gdx(3,3,iatm)+dz*vd3
        enddo
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    fock matrix.
c
c -- get the maximum absolute value of coefficients v and w
c
        call absmax(nat3,sv,isv,svmax)
        call absmax(3*nat3,sw,isw,swmax)
        tswmax=three*swmax
c -- global threshold testing
        vmx2=vmx*vmx
        smax = abs(sxx+sxy+sxz)
        dmax = abs(dx+dy+dz)
        vmax3 = (svmax+tswmax+tswmax)*vmx2
        vmax1 = (abra+two*smax)*vmx2
        dpmax2 = potxp2*dmax*(vmx2+vmx2)
        vmax2 = gwtmax*(potp*vmx2+dpmax2)
        vmax=max(vmax1,vmax2,vmax3)
        if(vmax.lt.thrsh)  go to 245
c
c -- numerical quadrature  for derivative fock matrix
c
        do 240 i=nbf(ipp)+1,nbf(ipp+1)
          vali = vao(i)
          valix = vaox(1,i)
          valiy = vaox(2,i)
          valiz = vaox(3,i)
          xvi=sxx*valix+sxy*valiy+sxz*valiz+ra*vali
          sxxi=sxx*vali
          sxyi=sxy*vali
          sxzi=sxz*vali
          xxvx=sxx*vaoxx(1,i)+sxy*vaoxx(2,i)+sxz*vaoxx(4,i)+ra*valix
          xxvy=sxx*vaoxx(2,i)+sxy*vaoxx(3,i)+sxz*vaoxx(5,i)+ra*valiy
          xxvz=sxx*vaoxx(4,i)+sxy*vaoxx(5,i)+sxz*vaoxx(6,i)+ra*valiz
          valm=abs(vali)
          valm1=abs(xvi)
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1=max(abs(xxvx),abs(xxvy),abs(xxvz))
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          val4=gwtmax*(potp*valm*vmx+dpmax2)
          valt1=max(val1,val3)
          valt2=max(val3,val4)
          valtest=max(valt1,valt2)
          if(valtest.gt.thrsh)then
            ii = inb(i)
            it = (ii*(ii-1))/2
            iatm = nbatm(ii)
            do 230 j=nbf(ipp)+1,i
              call zeroit(trp,nat3)
              jj = inb(j)
              ij = it + jj
              jatm = nbatm(jj)
              valj = vao(j)
              valjx = vaox(1,j)
              valjy = vaox(2,j)
              valjz = vaox(3,j)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              xvj=sxx*valjx+sxy*valjy+sxz*valjz
              if(valmx1*abs(valj)+valmx*abs(xvj).gt.thrsh)then
              if(iatm.ne.icntr) then
                trp(1,iatm) = trp(1,iatm) - xxvx*valj - valix*xvj
                trp(2,iatm) = trp(2,iatm) - xxvy*valj - valiy*xvj
                trp(3,iatm) = trp(3,iatm) - xxvz*valj - valiz*xvj
              endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
              if(jatm.ne.icntr) then
                trp(1,jatm) = trp(1,jatm) - xvi*valjx -
     $          vaoxx(1,j)*sxxi - vaoxx(2,j)*sxyi - vaoxx(4,j)*sxzi
                trp(2,jatm) = trp(2,jatm) - xvi*valjy -
     $          vaoxx(2,j)*sxxi - vaoxx(3,j)*sxyi - vaoxx(5,j)*sxzi
                trp(3,jatm) = trp(3,jatm) - xvi*valjz -
     $          vaoxx(4,j)*sxxi - vaoxx(5,j)*sxyi - vaoxx(6,j)*sxzi
              endif
              endif
              vij=vali*valj
              vijx=valix*valj+valjx*vali
              vijy=valiy*valj+valjy*vali
              vijz=valiz*valj+valjz*vali
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).gt.thrsh)then
                do katm=1,natoms
                  trp(1,katm) = trp(1,katm) + (sv(1,katm)*vij+
     $            sw(1,1,katm)*vijx+sw(2,1,katm)*vijy+sw(3,1,katm)*vijz)
                  trp(2,katm) = trp(2,katm) + (sv(2,katm)*vij+
     $            sw(1,2,katm)*vijx+sw(2,2,katm)*vijy+sw(3,2,katm)*vijz)
                  trp(3,katm) = trp(3,katm) + (sv(3,katm)*vij+
     $            sw(1,3,katm)*vijx+sw(2,3,katm)*vijy+sw(3,3,katm)*vijz)
                enddo
              endif
c
c -- weight derivatives contribution
c
              gij=potp*vij+potxp2*(dx*vijx+dy*vijy+dz*vijz)
              if(gwtmax*abs(gij).gt.thrsh)then
                do ia=1,natoms
                  if(ia.ne.icntr)then
                    trp(1,ia)=trp(1,ia)+gwt(1,ia,ipp)*gij
                    trp(2,ia)=trp(2,ia)+gwt(2,ia,ipp)*gij
                    trp(3,ia)=trp(3,ia)+gwt(3,ia,ipp)*gij
                  endif
                enddo
              endif
c
c -- apply translational invariance and accumulate contributions
c    into the fda array
c
              do ia=1,natoms
                if(ia.ne.icntr)then
                  trp(1,icntr)=trp(1,icntr)-trp(1,ia)
                  trp(2,icntr)=trp(2,icntr)-trp(2,ia)
                  trp(3,icntr)=trp(3,icntr)-trp(3,ia)
                endif
              enddo
              do ia=nb,ne
                fda(1,ia,ij)=fda(1,ia,ij)+trp(1,ia)
                fda(2,ia,ij)=fda(2,ia,ij)+trp(2,ia)
                fda(3,ia,ij)=fda(3,ia,ij)+trp(3,ia)
            EndDO
 230        continue
          endif
 240    continue
c
c   numerical quadrature for direct contribution to hessian matrix
c
 245    continue
        call secund(t5)
        tqf=tqf+t5-t4
c
c -- direct contribution to hessian matrix.
c    a factor two is applied
c
        ra2=two*ra
        rara2=two*rara
        do iatm=1,natoms
          hdx=rara2*denx(1,iatm)+prg2*gx(1,iatm)
          hdy=rara2*denx(2,iatm)+prg2*gx(2,iatm)
          hdz=rara2*denx(3,iatm)+prg2*gx(3,iatm)
          hgx=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          hgy=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          hgz=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          do jatm=iatm,natoms
       hess(1,iatm,1,jatm)=hess(1,iatm,1,jatm)+hdx*denx(1,jatm)+
     $ hgx*gx(1,jatm)+ra2*denxx(1,iatm,1,jatm)+pg2*gxx(1,iatm,1,jatm)
       hess(2,iatm,1,jatm)=hess(2,iatm,1,jatm)+hdy*denx(1,jatm)+
     $ hgy*gx(1,jatm)+ra2*denxx(2,iatm,1,jatm)+pg2*gxx(2,iatm,1,jatm)
       hess(3,iatm,1,jatm)=hess(3,iatm,1,jatm)+hdz*denx(1,jatm)+
     $ hgz*gx(1,jatm)+ra2*denxx(3,iatm,1,jatm)+pg2*gxx(3,iatm,1,jatm)
       hess(1,iatm,2,jatm)=hess(1,iatm,2,jatm)+hdx*denx(2,jatm)+
     $ hgx*gx(2,jatm)+ra2*denxx(1,iatm,2,jatm)+pg2*gxx(1,iatm,2,jatm)
       hess(2,iatm,2,jatm)=hess(2,iatm,2,jatm)+hdy*denx(2,jatm)+
     $ hgy*gx(2,jatm)+ra2*denxx(2,iatm,2,jatm)+pg2*gxx(2,iatm,2,jatm)
       hess(3,iatm,2,jatm)=hess(3,iatm,2,jatm)+hdz*denx(2,jatm)+
     $ hgz*gx(2,jatm)+ra2*denxx(3,iatm,2,jatm)+pg2*gxx(3,iatm,2,jatm)
       hess(1,iatm,3,jatm)=hess(1,iatm,3,jatm)+hdx*denx(3,jatm)+
     $ hgx*gx(3,jatm)+ra2*denxx(1,iatm,3,jatm)+pg2*gxx(1,iatm,3,jatm)
       hess(2,iatm,3,jatm)=hess(2,iatm,3,jatm)+hdy*denx(3,jatm)+
     $ hgy*gx(3,jatm)+ra2*denxx(2,iatm,3,jatm)+pg2*gxx(2,iatm,3,jatm)
       hess(3,iatm,3,jatm)=hess(3,iatm,3,jatm)+hdz*denx(3,jatm)+
     $ hgz*gx(3,jatm)+ra2*denxx(3,iatm,3,jatm)+pg2*gxx(3,iatm,3,jatm)
          enddo
        enddo
c
c  --  add weight derivatives contribution.
c      the sv array is used to store a partial sum
c
        do ia=1,natoms
          sv(1,ia)=potp2*denx(1,ia)+potxp2*gx(1,ia)
          sv(2,ia)=potp2*denx(2,ia)+potxp2*gx(2,ia)
          sv(3,ia)=potp2*denx(3,ia)+potxp2*gx(3,ia)
        enddo
        do ia=1,natoms
        if(ia.ne.icntr)then
          do ja=ia,natoms
          if(ja.ne.icntr)then
            hess(1,ia,1,ja)=hess(1,ia,1,ja)+gwt(1,ia,ipp)*sv(1,ja)
     $               +gwt(1,ja,ipp)*sv(1,ia)+hwt(1,ia,1,ja,ipp)
            hess(2,ia,1,ja)=hess(2,ia,1,ja)+gwt(2,ia,ipp)*sv(1,ja)
     $               +gwt(1,ja,ipp)*sv(2,ia)+hwt(2,ia,1,ja,ipp)
            hess(3,ia,1,ja)=hess(3,ia,1,ja)+gwt(3,ia,ipp)*sv(1,ja)
     $               +gwt(1,ja,ipp)*sv(3,ia)+hwt(3,ia,1,ja,ipp)
            hess(1,ia,2,ja)=hess(1,ia,2,ja)+gwt(1,ia,ipp)*sv(2,ja)
     $               +gwt(2,ja,ipp)*sv(1,ia)+hwt(1,ia,2,ja,ipp)
            hess(2,ia,2,ja)=hess(2,ia,2,ja)+gwt(2,ia,ipp)*sv(2,ja)
     $               +gwt(2,ja,ipp)*sv(2,ia)+hwt(2,ia,2,ja,ipp)
            hess(3,ia,2,ja)=hess(3,ia,2,ja)+gwt(3,ia,ipp)*sv(2,ja)
     $               +gwt(2,ja,ipp)*sv(3,ia)+hwt(3,ia,2,ja,ipp)
            hess(1,ia,3,ja)=hess(1,ia,3,ja)+gwt(1,ia,ipp)*sv(3,ja)
     $               +gwt(3,ja,ipp)*sv(1,ia)+hwt(1,ia,3,ja,ipp)
            hess(2,ia,3,ja)=hess(2,ia,3,ja)+gwt(2,ia,ipp)*sv(3,ja)
     $               +gwt(3,ja,ipp)*sv(2,ia)+hwt(2,ia,3,ja,ipp)
            hess(3,ia,3,ja)=hess(3,ia,3,ja)+gwt(3,ia,ipp)*sv(3,ja)
     $               +gwt(3,ja,ipp)*sv(3,ia)+hwt(3,ia,3,ja,ipp)
          endif
          enddo
        endif
        enddo
        call secund(t6)
        tqh=tqh+t6-t5
 200  continue
      endif
c
      return
      end
c =====================================================================
      subroutine cwmpfhu(dft,    npp,    nbas,   nbf,    natoms,
     $                   natonce, nb,    ne,     nbatm,  thrsh,
     $                   da,     db,     dm,     gdena,  gdenb,
     $                   wght,   pra,    prb,    prara,  prbrb,
     $                   prarb,  pga,    pgb,    pgc,    praga,
     $                   pragb,  pragc,  prbga,  prbgb,  prbgc,
     $                   pgaga,  pgagb,  pgagc,  pgbgb,  pgbgc,
     $                   pgcgc,  vao,    vaox,   vaoxx,  vaoxxx,
     $                   inb,    vm,     indx,   denxa,  denxb,
     $                   denxxa, denxxb, gdxa,   gdxb,   gdxxa,
     $                   gdxxb,  gxaa,   gxbb,   gxab,   gxxaa,
     $                   gxxbb,  gxxab,  sva,    svb,    swa,
     $                   swb,    sga,    sgb,    sgc,    icntr,
     $                   gwt,    hwt,    trpa,   trpb,   fda,
     $                   fdb,    hess,   td1,    tg1,    tsw,
     $                   tqf,    tqh)
c     implicit real*8 (a-h,o-z)
      implicit none             ! trying to avoid typing errors
c
c
c  carries out numerical integration and accumulates contribution
c  from current grid point into derivative fock matrices and hessian
c  taking into account the weight derivatives.
c
c  natonce components of the derivative fock matrices are computed in
c  one pass (components nb to ne, with nb <= icntr <= ne). The
c  icntr component of the derivative fock matrix is computed
c  by applying translational invariance at each grid point.
c  The Hessian components involving icntr are not computed, as
c  they can be obtained by translational invariance at the end
c  of the calculation
c
c  ** open shell **
c
c  arguments
c
c  dft     -  method flag (note: all methods include slater exchange)
c              1 - 3 - local correlation
c             >3 - nonlocal
c  npp     -  number of contributing (non-zero) grid points this batch
c  nbas    -  total number of basis functions
c  nbf     -  indexing array to location of "non-zero" basis functions
c  natoms  -  number of atoms
c  natonce -  number of atomic derivatives to compute in this pass
c  nb,ne   -  first and last deivative to compute (natonce=ne-nb+1)
c  nbatm   -  basis functions --> atom index
c  thrsh   -  threshold for neglect of contribution
c  da      -  alpha density matrix (lower triangle)
c  db      -  beta density matrix (lower triangle)
c  dm      -  maximum density matrix element per column
c  gdena   -  alpha density gradient at grid points (only dft > 3)
c  gdena   -  beta density gradient at grid points (dft > 3)
c  wght    -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prb   -  Functional derivative w.r.t. beta density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgb   -  Funct. deriv. w.r.t. beta gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  prbga -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha grad.
C           (dft > 3)
C  prbgb -  Funct. 2nd. deriv. w.r.t. beta dens. and  grad. (dft > 3)
C  prbgc -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgbgb -  Funct. 2nd. deriv. w.r.t. beta grad. (dft > 3)
C  pgbgc -  Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
c  vao     -  "non-zero" basis function values at grid points
c  vaox    -  basis function 1st derivatives at grid points
c  vaoxx   -  basis function 2nd derivatives at grid points
c  vaoxxx  -  basis function 3rd derivatives at grid points (dft > 3)
c  inb     -  indexing array for non-zero entries to vao
c  vm      -  array containing maximum magnitude ao per grid point
c  indx    -  index into contributing columns of vao
c  denxa   -  scratch storage for alpah density gradient per grid point
c  denxb   -  scratch storage for beta density gradient per grid point
c  denxxa  -  ditto for alpha density hessian per grid point
c  denxxb  -  ditto for beta density hessian per grid point
c  gdxa    -  ditto for alpha atomic gradient of density gradient
c             (dft > 3)
c  gdxb    -  ditto for beta atomic gradient of density gradient
c             (dft > 3)
c  gdxxa   -  ditto for atomic hessian of alpha density gradient
c             (dft > 3)
c  gdxxb   -  ditto for atomic hessian of beta density gradient
c             (dft > 3)
c  gxaa    -  ditto for atomic gradient of alpha gradient invariant
c             (dft > 3)
c  gxbb    -  ditto for atomic gradient of beta gradient invariant
c             (dft > 3)
c  gxbb    -  ditto for atomic gradient of alpha beta gradient
c             invariant (dft > 3)
c  gxxaa   -  ditto for atomic hessian of alpha gradient invariant
c             (dft > 3)
c  gxxbb   -  ditto for atomic hessian of beta gradient invariant
c             (dft > 3)
c  gxxbb   -  ditto for atomic hessian of alpha beta gradient invariant
c             (dft > 3)
c  sva     -  storage for alpha coefficient of vao(i)*vao(j) in
c             quadrature formula
c  svb     -  storage for beta coefficient of vao(i)*vao(j) in
c             quadrature formula
c  swa     -  storage for alpha coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  swb     -  storage for beta coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  sga     -  storage for partial sum for hessian quadrature (dft > 3)
c  sgb     -  storage for partial sum for hessian quadrature (dft > 3)
c  sgc     -  storage for partial sum for hessian quadrature (dft > 3)
c  icntr   -  current atomic center
c  gwt     -  gradient of quadrature weight
c  hwt     -  hessian of quadrature weight
c  trpa    -  scratch area for translational invariance, alpha
c  trpb    -  scratch area for translational invariance, beta
c
c  on exit
c
c  fda     -  contribution to alpha derivative fock matrices
c  fdb     -  contribution to beta derivative fock matrices
c  hess    -  contribution to hessian matrix
c
c
c
c
      integer npp,nbas,natoms,icntr,natonce,nb,ne
      real*8 wght(*),gwt(3,natoms,*),hwt(3,natoms,3,natoms,*)
      real*8 vao(*),vaox(3,*),vaoxx(6,*),vaoxxx(10,*)
      real*8 vm(*),da(nbas*(nbas+1)/2),dm(nbas),gdena(3,*)
      real*8 db(nbas*(nbas+1)/2),gdenb(3,*)
      real*8 pra(npp),pga(npp),pgc(npp),prara(npp)
      real*8 prb(npp),pgb(npp),prbrb(npp)
      real*8 prarb(npp),praga(npp),pragb(npp),pragc(npp)
      real*8 prbga(npp),prbgb(npp),prbgc(npp)
      real*8 pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
      real*8 pgbgb(npp),pgbgc(npp)
      real*8 denxa(3,natoms),denxxa(3,natoms,3,natoms)
      real*8 denxb(3,natoms),denxxb(3,natoms,3,natoms)
      real*8 gdxa(3,3,natoms),gdxxa(3,3,natoms,3,natoms)
      real*8 gdxb(3,3,natoms),gdxxb(3,3,natoms,3,natoms)
      real*8 gxaa(3,natoms),gxxaa(3,natoms,3,natoms)
      real*8 gxbb(3,natoms),gxxbb(3,natoms,3,natoms)
      real*8 gxab(3,natoms),gxxab(3,natoms,3,natoms)
      real*8 sva(3,natoms),svb(3,natoms)
      real*8 swa(3,3,natoms),swb(3,3,natoms)
      real*8 sga(3,natoms),sgb(3,natoms),sgc(3,natoms)
      integer dft,nbf(*),inb(*),indx(npp),nbatm(nbas)
      real*8 fda(3,nb:ne,*),hess(3,natoms,3,natoms)
      real*8 fdb(3,nb:ne,*)
      real*8 trpa(3,natoms),trpb(3,natoms)
      real*8 thrsh
      real*8 td1,tg1,tsw,tqf,tqh
c
      real*8 zero,two,three,four,six,twelve,half
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0,half=0.5d0)
      real*8 epsi,alot
      parameter (epsi=1.0d-15,alot=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,doval
      integer nat3,ip,ipp,i,j,ii,jj,ij,it
      integer ia,ja,iatm,jatm,katm,iixx,k1
      real*8 t1,t2,t3,t4,t5,t6
      real*8 wg,vmx,vmx2,prap,prbp,gwtmax,pgap,pgbp,pgcp
      real*8 ra,rb,rara,rbrb,rarb,abrr,abr
      real*8 ga,gb,gc,raga,ragb,ragc,rbga,rbgb,rbgc
      real*8 gaga,gagb,gagc,gbgb,gbgc,gcgc
      real*8 dax,day,daz,dbx,dby,dbz
      real*8 dax2,day2,daz2,dbx2,dby2,dbz2
      real*8 ga2,gb2,pgap2,pgbp2
      real*8 thrx,thrxx,thrx1,thrsh1
      real*8 dmx2,daij,dbij,daij2,dbij2,abdij,abdijm
      real*8 valx,valy,valz,valm,valt,abval
      real*8 vali,vala,valxa,valya,valza
      real*8 valj,valb,valxb,valyb,valzb
      real*8 valij,valija,valijb
      real*8 xyc,xzc,yzc
      real*8 svamax,svbmax,svmax,vmax,abgmax
      real*8 swamax,swbmax,swmax,tswmax,vmax1,vmax3
      real*8 valxx,valxy,valxz,valyy,valyz,valzz
      real*8 valmm,valm1,valmm1,valmm2
      real*8 valmx,valmx1a,valmx1b,valmx1,val1,val2,val3,valtest
      real*8 thtest,valmxj,xvja,xvjb,abxvj
      real*8 valjx,valjy,valjz,abvj,abvvj,abijt,abij
      real*8 valix,valiy,valiz,xvia,xvib
      real*8 sxxa,sxya,sxza,sxxb,sxyb,sxzb
      real*8 sxxia,sxyia,sxzia,sxxib,sxyib,sxzib,smax,smaxmx
      real*8 xxvxa,xxvya,xxvza,xxvxb,xxvyb,xxvzb
      real*8 vij,vijx,vijy,vijz
      real*8 dxija,dxijb,gija,gijb
c
c consistency check
c
      if(icntr.lt.nb.or.icntr.gt.ne)then
        call nerror(1,'cwmphfu','icntr is wrong',icntr,nb)
      endif
c
      nat3 = 3*natoms
c
      if(dft.le.3) then
c
c  local only
c
      do 50 ip=1,npp
      ipp = indx(ip)
      wg  = wght(ipp)
      prap=pra(ip)
      prbp=prb(ip)
      ra = prap*wg
      rb = prbp*wg
      rara = prara(ip)*wg
      rbrb = prbrb(ip)*wg
      rarb = prarb(ip)*wg
      vmx = vm(ipp)
c
c   initializations for weight derivatives
c
      gwtmax=zero
      do ia=1,natoms
        if(ia.ne.icntr)then
          gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +    abs(gwt(3,ia,ipp)))
        endif
      enddo
      call zeroit(denxa,nat3)
      call zeroit(denxb,nat3)
      call zeroit(denxxa,nat3**2)
      call zeroit(denxxb,nat3**2)
c
c    form gradient and hessian density at current grid point.
c
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrr=max(abs(rara),abs(rbrb))+abs(rarb)
      dodenx=abrr.gt.epsi
      thrx=alot
      if(dodenx)thrx=thrsh/abrr
      thrx1=thrx/vmx2
      abr=max(abs(ra),abs(rb))
      dodenxx=abr.gt.epsi
      thrxx=alot
      if(dodenxx)thrxx=thrsh/abr
      if(.not.(dodenx.or.dodenxx))goto 25 ! should never happen!
c
c   now compute density derivatives
c
      do 20 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        dmx2 = dm(ii)
        if(dmx2*vmx2.lt.thrsh) go to 20
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        if(iatm.eq.icntr)goto 20
        valx = vaox(1,i)
        valy = vaox(2,i)
        valz = vaox(3,i)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx2*vmx
        dox=(dodenx.and.(two*valt.gt.thrx1.or.valt**2.gt.thrx))
        doxx=(dodenxx.and.dmx2*vmx2.gt.thrxx)
        if(.not.(dox.or.doxx)) goto 20
        do 18 j=nbf(ipp)+1,i
          jj = inb(j)
          ij = it + jj
          jatm = nbatm(jj)
          daij2=two*da(ij)
          dbij2=two*db(ij)
          daij = daij2*vao(j)
          dbij = dbij2*vao(j)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denxa(1,iatm) = denxa(1,iatm) - daij*valx
              denxa(2,iatm) = denxa(2,iatm) - daij*valy
              denxa(3,iatm) = denxa(3,iatm) - daij*valz
c
              denxb(1,iatm) = denxb(1,iatm) - dbij*valx
              denxb(2,iatm) = denxb(2,iatm) - dbij*valy
              denxb(3,iatm) = denxb(3,iatm) - dbij*valz
          endif
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
          if(doxx.and.abdijm.gt.thrxx)then
             denxxa(1,iatm,1,iatm)=denxxa(1,iatm,1,iatm)+daij*vaoxx(1,i)
             xyc=daij*vaoxx(2,i)
             denxxa(1,iatm,2,iatm)=denxxa(1,iatm,2,iatm)+xyc
             denxxa(2,iatm,1,iatm)=denxxa(2,iatm,1,iatm)+xyc
             denxxa(2,iatm,2,iatm)=denxxa(2,iatm,2,iatm)+daij*vaoxx(3,i)
             xzc=daij*vaoxx(4,i)
             denxxa(1,iatm,3,iatm)=denxxa(1,iatm,3,iatm)+xzc
             denxxa(3,iatm,1,iatm)=denxxa(3,iatm,1,iatm)+xzc
             yzc=daij*vaoxx(5,i)
             denxxa(2,iatm,3,iatm)=denxxa(2,iatm,3,iatm)+yzc
             denxxa(3,iatm,2,iatm)=denxxa(3,iatm,2,iatm)+yzc
             denxxa(3,iatm,3,iatm)=denxxa(3,iatm,3,iatm)+daij*vaoxx(6,i)
c
             denxxb(1,iatm,1,iatm)=denxxb(1,iatm,1,iatm)+dbij*vaoxx(1,i)
             xyc=dbij*vaoxx(2,i)
             denxxb(1,iatm,2,iatm)=denxxb(1,iatm,2,iatm)+xyc
             denxxb(2,iatm,1,iatm)=denxxb(2,iatm,1,iatm)+xyc
             denxxb(2,iatm,2,iatm)=denxxb(2,iatm,2,iatm)+dbij*vaoxx(3,i)
             xzc=dbij*vaoxx(4,i)
             denxxb(1,iatm,3,iatm)=denxxb(1,iatm,3,iatm)+xzc
             denxxb(3,iatm,1,iatm)=denxxb(3,iatm,1,iatm)+xzc
             yzc=dbij*vaoxx(5,i)
             denxxb(2,iatm,3,iatm)=denxxb(2,iatm,3,iatm)+yzc
             denxxb(3,iatm,2,iatm)=denxxb(3,iatm,2,iatm)+yzc
             denxxb(3,iatm,3,iatm)=denxxb(3,iatm,3,iatm)+dbij*vaoxx(6,i)
          endif
          if(jatm.lt.iatm) go to 18
          if(jatm.eq.icntr) go to 18
          abdijm=max(abs(da(ij)),abs(db(ij)))*valm*vmx
          if(doxx.and.abdijm.gt.thrxx)then
c
c    (b) iatm with jatm
c
              daij = daij2*vaox(1,j)
              denxxa(1,iatm,1,jatm)=denxxa(1,iatm,1,jatm)+daij*valx
              denxxa(2,iatm,1,jatm)=denxxa(2,iatm,1,jatm)+daij*valy
              denxxa(3,iatm,1,jatm)=denxxa(3,iatm,1,jatm)+daij*valz
              daij = daij2*vaox(2,j)
              denxxa(1,iatm,2,jatm)=denxxa(1,iatm,2,jatm)+daij*valx
              denxxa(2,iatm,2,jatm)=denxxa(2,iatm,2,jatm)+daij*valy
              denxxa(3,iatm,2,jatm)=denxxa(3,iatm,2,jatm)+daij*valz
              daij = daij2*vaox(3,j)
              denxxa(1,iatm,3,jatm)=denxxa(1,iatm,3,jatm)+daij*valx
              denxxa(2,iatm,3,jatm)=denxxa(2,iatm,3,jatm)+daij*valy
              denxxa(3,iatm,3,jatm)=denxxa(3,iatm,3,jatm)+daij*valz
c
              dbij = dbij2*vaox(1,j)
              denxxb(1,iatm,1,jatm)=denxxb(1,iatm,1,jatm)+dbij*valx
              denxxb(2,iatm,1,jatm)=denxxb(2,iatm,1,jatm)+dbij*valy
              denxxb(3,iatm,1,jatm)=denxxb(3,iatm,1,jatm)+dbij*valz
              dbij = dbij2*vaox(2,j)
              denxxb(1,iatm,2,jatm)=denxxb(1,iatm,2,jatm)+dbij*valx
              denxxb(2,iatm,2,jatm)=denxxb(2,iatm,2,jatm)+dbij*valy
              denxxb(3,iatm,2,jatm)=denxxb(3,iatm,2,jatm)+dbij*valz
              dbij = dbij2*vaox(3,j)
              denxxb(1,iatm,3,jatm)=denxxb(1,iatm,3,jatm)+dbij*valx
              denxxb(2,iatm,3,jatm)=denxxb(2,iatm,3,jatm)+dbij*valy
              denxxb(3,iatm,3,jatm)=denxxb(3,iatm,3,jatm)+dbij*valz
          endif
 18     continue
        do 19 j=i+1,nbf(ipp+1)
          jj = inb(j)
          ij = (jj*(jj-1))/2 + ii
          jatm = nbatm(jj)
          daij2 = two*da(ij)
          dbij2 = two*db(ij)
          daij = daij2*vao(j)
          dbij = dbij2*vao(j)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
          if(dox)then
c
c -- atomic gradient of density
c
              denxa(1,iatm) = denxa(1,iatm) - daij*valx
              denxa(2,iatm) = denxa(2,iatm) - daij*valy
              denxa(3,iatm) = denxa(3,iatm) - daij*valz
c
              denxb(1,iatm) = denxb(1,iatm) - dbij*valx
              denxb(2,iatm) = denxb(2,iatm) - dbij*valy
              denxb(3,iatm) = denxb(3,iatm) - dbij*valz
          endif
          if(doxx.and.abdijm.gt.thrxx)then
c
c -- atomic hessian of density
c -- (a) iatm with iatm
c
             denxxa(1,iatm,1,iatm)=denxxa(1,iatm,1,iatm)+daij*vaoxx(1,i)
             xyc=daij*vaoxx(2,i)
             denxxa(1,iatm,2,iatm)=denxxa(1,iatm,2,iatm)+xyc
             denxxa(2,iatm,1,iatm)=denxxa(2,iatm,1,iatm)+xyc
             denxxa(2,iatm,2,iatm)=denxxa(2,iatm,2,iatm)+daij*vaoxx(3,i)
             xzc=daij*vaoxx(4,i)
             denxxa(1,iatm,3,iatm)=denxxa(1,iatm,3,iatm)+xzc
             denxxa(3,iatm,1,iatm)=denxxa(3,iatm,1,iatm)+xzc
             yzc=daij*vaoxx(5,i)
             denxxa(2,iatm,3,iatm)=denxxa(2,iatm,3,iatm)+yzc
             denxxa(3,iatm,2,iatm)=denxxa(3,iatm,2,iatm)+yzc
             denxxa(3,iatm,3,iatm)=denxxa(3,iatm,3,iatm)+daij*vaoxx(6,i)
c
             denxxb(1,iatm,1,iatm)=denxxb(1,iatm,1,iatm)+dbij*vaoxx(1,i)
             xyc=dbij*vaoxx(2,i)
             denxxb(1,iatm,2,iatm)=denxxb(1,iatm,2,iatm)+xyc
             denxxb(2,iatm,1,iatm)=denxxb(2,iatm,1,iatm)+xyc
             denxxb(2,iatm,2,iatm)=denxxb(2,iatm,2,iatm)+dbij*vaoxx(3,i)
             xzc=dbij*vaoxx(4,i)
             denxxb(1,iatm,3,iatm)=denxxb(1,iatm,3,iatm)+xzc
             denxxb(3,iatm,1,iatm)=denxxb(3,iatm,1,iatm)+xzc
             yzc=dbij*vaoxx(5,i)
             denxxb(2,iatm,3,iatm)=denxxb(2,iatm,3,iatm)+yzc
             denxxb(3,iatm,2,iatm)=denxxb(3,iatm,2,iatm)+yzc
             denxxb(3,iatm,3,iatm)=denxxb(3,iatm,3,iatm)+dbij*vaoxx(6,i)
          endif
c
c    (b) iatm with jatm
c
          if(jatm.lt.iatm) go to 19
          if(jatm.eq.icntr) go to 19
          abdijm=max(abs(da(ij)),abs(db(ij)))*valm*vmx
          if(doxx.and.abdijm.gt.thrxx)then
              daij = daij2*vaox(1,j)
              denxxa(1,iatm,1,jatm)=denxxa(1,iatm,1,jatm)+daij*valx
              denxxa(2,iatm,1,jatm)=denxxa(2,iatm,1,jatm)+daij*valy
              denxxa(3,iatm,1,jatm)=denxxa(3,iatm,1,jatm)+daij*valz
              daij = daij2*vaox(2,j)
              denxxa(1,iatm,2,jatm)=denxxa(1,iatm,2,jatm)+daij*valx
              denxxa(2,iatm,2,jatm)=denxxa(2,iatm,2,jatm)+daij*valy
              denxxa(3,iatm,2,jatm)=denxxa(3,iatm,2,jatm)+daij*valz
              daij = daij2*vaox(3,j)
              denxxa(1,iatm,3,jatm)=denxxa(1,iatm,3,jatm)+daij*valx
              denxxa(2,iatm,3,jatm)=denxxa(2,iatm,3,jatm)+daij*valy
              denxxa(3,iatm,3,jatm)=denxxa(3,iatm,3,jatm)+daij*valz
c
              dbij = dbij2*vaox(1,j)
              denxxb(1,iatm,1,jatm)=denxxb(1,iatm,1,jatm)+dbij*valx
              denxxb(2,iatm,1,jatm)=denxxb(2,iatm,1,jatm)+dbij*valy
              denxxb(3,iatm,1,jatm)=denxxb(3,iatm,1,jatm)+dbij*valz
              dbij = dbij2*vaox(2,j)
              denxxb(1,iatm,2,jatm)=denxxb(1,iatm,2,jatm)+dbij*valx
              denxxb(2,iatm,2,jatm)=denxxb(2,iatm,2,jatm)+dbij*valy
              denxxb(3,iatm,2,jatm)=denxxb(3,iatm,2,jatm)+dbij*valz
              dbij = dbij2*vaox(3,j)
              denxxb(1,iatm,3,jatm)=denxxb(1,iatm,3,jatm)+dbij*valx
              denxxb(2,iatm,3,jatm)=denxxb(2,iatm,3,jatm)+dbij*valy
              denxxb(3,iatm,3,jatm)=denxxb(3,iatm,3,jatm)+dbij*valz
          endif
 19     continue
 20   continue
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c -- now form the coefficients that multiply the basis functions
c    in the expression for the Fock Matrix
c    (i.e., quantities V at page 7436 of Johnson and Frisch).
c
        DO IAtm=1,NAtoms
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)
          svb(1,iatm)=rbrb*denxb(1,iatm)+rarb*denxa(1,iatm)
          svb(2,iatm)=rbrb*denxb(2,iatm)+rarb*denxa(2,iatm)
          svb(3,iatm)=rbrb*denxb(3,iatm)+rarb*denxa(3,iatm)
        EndDO
c
c    numerical quadrature for the derivative fock matrix.
c
c -- get maximum absolute value of sva and svb
c
      Call absmax(NAt3,sva,iixx,svamax)
      Call absmax(NAt3,svb,iixx,svbmax)
      svmax=max(svamax,svbmax)
c
c -- global threshold testing
c
      VMax = Max(Abr*VMx2,Abrr*svmax*VMx2)
      if(vmax.lt.thrsh) go to 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      do 40 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        vali = vao(i)
        vala = vali*ra
        valb = vali*rb
        valxa = vaox(1,i)*ra
        valya = vaox(2,i)*ra
        valza = vaox(3,i)*ra
        valxb = vaox(1,i)*rb
        valyb = vaox(2,i)*rb
        valzb = vaox(3,i)*rb
        abval= max(abs(vala),abs(valb))
        doval= abval.gt.thrsh1
        valm = max(abs(vaox(1,i)),abs(vaox(2,i)),abs(vaox(3,i)))*abr
        Valt = MAX(abval,valm,Abs(vali)*svmax)
        if(valt.gt.thrsh1) then
          do 30 j=nbf(ipp)+1,i
            call zeroit(trpa,nat3)
            call zeroit(trpb,nat3)
            jj = inb(j)
            ij = it + jj
            jatm = nbatm(jj)
            valj = vao(j)
            if(abs(valj)*valm.gt.thrsh)then
              if(iatm.ne.icntr) then
                trpa(1,iatm) = trpa(1,iatm) - valxa*valj
                trpa(2,iatm) = trpa(2,iatm) - valya*valj
                trpa(3,iatm) = trpa(3,iatm) - valza*valj
                trpb(1,iatm) = trpb(1,iatm) - valxb*valj
                trpb(2,iatm) = trpb(2,iatm) - valyb*valj
                trpb(3,iatm) = trpb(3,iatm) - valzb*valj
              endif
            endif
            if(doval)then
              if(jatm.ne.icntr) then
                trpa(1,jatm) = trpa(1,jatm) - vala*vaox(1,j)
                trpa(2,jatm) = trpa(2,jatm) - vala*vaox(2,j)
                trpa(3,jatm) = trpa(3,jatm) - vala*vaox(3,j)
                trpb(1,jatm) = trpb(1,jatm) - valb*vaox(1,j)
                trpb(2,jatm) = trpb(2,jatm) - valb*vaox(2,j)
                trpb(3,jatm) = trpb(3,jatm) - valb*vaox(3,j)
              endif
            endif
c
c -- contribution of functional derivative
c
            valij = vali*valj
            if(abs(valij)*svmax.gt.thrsh) then
              do katm=1,natoms
                trpa(1,katm)=trpa(1,katm) + valij*sva(1,katm)
                trpa(2,katm)=trpa(2,katm) + valij*sva(2,katm)
                trpa(3,katm)=trpa(3,katm) + valij*sva(3,katm)
                trpb(1,katm)=trpb(1,katm) + valij*svb(1,katm)
                trpb(2,katm)=trpb(2,katm) + valij*svb(2,katm)
                trpb(3,katm)=trpb(3,katm) + valij*svb(3,katm)
              enddo
            endif
c
c -- weight derivatives contribution
c
            valija=valij*prap
            valijb=valij*prbp
            do ia=1,natoms
              if(ia.ne.icntr)then
                trpa(1,ia)=trpa(1,ia)+valija*gwt(1,ia,ipp)
                trpa(2,ia)=trpa(2,ia)+valija*gwt(2,ia,ipp)
                trpa(3,ia)=trpa(3,ia)+valija*gwt(3,ia,ipp)
                trpb(1,ia)=trpb(1,ia)+valijb*gwt(1,ia,ipp)
                trpb(2,ia)=trpb(2,ia)+valijb*gwt(2,ia,ipp)
                trpb(3,ia)=trpb(3,ia)+valijb*gwt(3,ia,ipp)
              endif
            enddo
c
c -- apply translational invariance and accumulate contributions
c    into the fda array
c
            do ia=1,natoms
              if(ia.ne.icntr)then
                trpa(1,icntr)=trpa(1,icntr)-trpa(1,ia)
                trpa(2,icntr)=trpa(2,icntr)-trpa(2,ia)
                trpa(3,icntr)=trpa(3,icntr)-trpa(3,ia)
                trpb(1,icntr)=trpb(1,icntr)-trpb(1,ia)
                trpb(2,icntr)=trpb(2,icntr)-trpb(2,ia)
                trpb(3,icntr)=trpb(3,icntr)-trpb(3,ia)
              endif
            enddo
            do ia=nb,ne
              fda(1,ia,ij)=fda(1,ia,ij)+trpa(1,ia)
              fda(2,ia,ij)=fda(2,ia,ij)+trpa(2,ia)
              fda(3,ia,ij)=fda(3,ia,ij)+trpa(3,ia)
              fdb(1,ia,ij)=fdb(1,ia,ij)+trpb(1,ia)
              fdb(2,ia,ij)=fdb(2,ia,ij)+trpb(2,ia)
              fdb(3,ia,ij)=fdb(3,ia,ij)+trpb(3,ia)
            enddo
 30       continue
        endif
 40   continue
c
c    numerical quadrature for the hessian matrix
c
 45   continue
      call secund(t3)
      tqf=tqf+t3-t2
c
c -- direct contribution to hessian matrix.
c
      do iatm=1,natoms
        do jatm=iatm,natoms
          hess(1,iatm,1,jatm) = hess(1,iatm,1,jatm) +
     $     sva(1,iatm)*denxa(1,jatm) + svb(1,iatm)*denxb(1,jatm)+
     $     ra*denxxa(1,iatm,1,jatm) + rb*denxxb(1,iatm,1,jatm)
          hess(2,iatm,1,jatm) = hess(2,iatm,1,jatm) +
     $     sva(2,iatm)*denxa(1,jatm) + svb(2,iatm)*denxb(1,jatm)+
     $     ra*denxxa(2,iatm,1,jatm) + rb*denxxb(2,iatm,1,jatm)
          hess(3,iatm,1,jatm) = hess(3,iatm,1,jatm) +
     $     sva(3,iatm)*denxa(1,jatm) + svb(3,iatm)*denxb(1,jatm)+
     $     ra*denxxa(3,iatm,1,jatm) + rb*denxxb(3,iatm,1,jatm)
          hess(1,iatm,2,jatm) = hess(1,iatm,2,jatm) +
     $     sva(1,iatm)*denxa(2,jatm) + svb(1,iatm)*denxb(2,jatm)+
     $     ra*denxxa(1,iatm,2,jatm) + rb*denxxb(1,iatm,2,jatm)
          hess(2,iatm,2,jatm) = hess(2,iatm,2,jAtm) +
     $     sva(2,iatm)*denxa(2,jatm) + svb(2,iatm)*denxb(2,jatm)+
     $     ra*denxxa(2,iatm,2,jatm) + rb*denxxb(2,iatm,2,jatm)
          hess(3,iatm,2,jatm) = hess(3,iatm,2,jatm) +
     $     sva(3,iatm)*denxa(2,jatm) + svb(3,iatm)*denxb(2,jatm)+
     $     ra*denxxa(3,iatm,2,jatm) + rb*denxxb(3,iatm,2,jatm)
          hess(1,iatm,3,jatm) = hess(1,iatm,3,jatm) +
     $     sva(1,iatm)*denxa(3,jatm) + svb(1,iatm)*denxb(3,jatm)+
     $     ra*denxxa(1,iatm,3,jatm) + rb*denxxb(1,iatm,3,jatm)
          hess(2,iatm,3,jatm) = hess(2,iatm,3,jatm) +
     $     sva(2,iatm)*denxa(3,jatm) + svb(2,iatm)*denxb(3,jatm)+
     $     ra*denxxa(2,iatm,3,jatm) + rb*denxxb(2,iatm,3,jatm)
          hess(3,iatm,3,jatm) = hess(3,iatm,3,jatm) +
     $     sva(3,iatm)*denxa(3,jatm) + svb(3,iatm)*denxb(3,jatm)+
     $     ra*denxxa(3,iatm,3,jatm) + rb*denxxb(3,iatm,3,jatm)
        EndDo
      EndDo
c
c  --  add weight derivatives contribution
c
c      the sva array is used to store a partial sum
c
      do ia=1,natoms
        sva(1,ia)=prap*denxa(1,ia)+prbp*denxb(1,ia)
        sva(2,ia)=prap*denxa(2,ia)+prbp*denxb(2,ia)
        sva(3,ia)=prap*denxa(3,ia)+prbp*denxb(3,ia)
      enddo
      do ia=1,natoms
      if(ia.ne.icntr)then
        do ja=ia,natoms
        if(ja.ne.icntr)then
          hess(1,ia,1,ja)=hess(1,ia,1,ja)+gwt(1,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(1,ia)+hwt(1,ia,1,ja,ipp)
          hess(2,ia,1,ja)=hess(2,ia,1,ja)+gwt(2,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(2,ia)+hwt(2,ia,1,ja,ipp)
          hess(3,ia,1,ja)=hess(3,ia,1,ja)+gwt(3,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(3,ia)+hwt(3,ia,1,ja,ipp)
          hess(1,ia,2,ja)=hess(1,ia,2,ja)+gwt(1,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(1,ia)+hwt(1,ia,2,ja,ipp)
          hess(2,ia,2,ja)=hess(2,ia,2,ja)+gwt(2,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(2,ia)+hwt(2,ia,2,ja,ipp)
          hess(3,ia,2,ja)=hess(3,ia,2,ja)+gwt(3,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(3,ia)+hwt(3,ia,2,ja,ipp)
          hess(1,ia,3,ja)=hess(1,ia,3,ja)+gwt(1,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(1,ia)+hwt(1,ia,3,ja,ipp)
          hess(2,ia,3,ja)=hess(2,ia,3,ja)+gwt(2,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(2,ia)+hwt(2,ia,3,ja,ipp)
          hess(3,ia,3,ja)=hess(3,ia,3,ja)+gwt(3,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(3,ia)+hwt(3,ia,3,ja,ipp)
        endif
        enddo
      endif
      enddo
      call secund(t4)
      tqh=tqh+t4-t3
c
c   local dft ends here
c
 50   continue
cc
      else
cc
c
c   non-local dft
c
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        prap=pra(ip)
        prbp=prb(ip)
        ra = prap*wg
        rb = prbp*wg
        pgap = pga(ip)
        pgbp = pgb(ip)
        pgcp = pgc(ip)
        pgap2=pgap+pgap
        pgbp2=pgbp+pgbp
        ga = pgap*wg
        gb = pgbp*wg
        gc = pgcp*wg
        rara = prara(ip)*wg
        rbrb = prbrb(ip)*wg
        rarb = prarb(ip)*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        rbga = prbga(ip)*wg
        rbgb = prbgb(ip)*wg
        rbgc = prbgc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gbgb = pgbgb(ip)*wg
        gbgc = pgbgc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  density gradient at current point
c
        DAX=GDENA(1,IPP)
        DAY=GDENA(2,IPP)
        DAZ=GDENA(3,IPP)
        DBX=GDENB(1,IPP)
        DBY=GDENB(2,IPP)
        DBZ=GDENB(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenXA,NAt3)
        CALL ZeroIT(DenXXA,NAt3**2)
        CALL ZeroIT(GDXA,3*NAt3)
        CALL ZeroIT(GDXXA,3*NAt3**2)
        CALL ZeroIT(DenXB,NAt3)
        CALL ZeroIT(DenXXB,NAt3**2)
        CALL ZeroIT(GDXB,3*NAt3)
        CALL ZeroIT(GDXXB,3*NAt3**2)
c   initializations for weight derivatives
c
        gwtmax=zero
        do ia=1,natoms
          if(ia.ne.icntr)then
            gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +        abs(gwt(3,ia,ipp)))
          endif
        enddo
c
c    form the atomic gradient and Hessian  of density
c    and density gradient at current grid point
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c    in addition, we check whether the second derivatives of the
c    density gradient should be computed at all.
c
        abgmax=max(abs(ga),abs(gb),abs(gc))
        dogradxx=abgmax.gt.epsi
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          dmx2 = DM(II)*two   ! max. element of first order density
          if(dmx2.gt.epsi)then
             thtest=thrsh/(vmx*dmx2)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          if(iatm.eq.icntr)goto 220
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2=two*DA(IJ)
            DBIJ2=two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
              DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
              DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
              DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
              DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
              DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VALZZ
c
              DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VALXX
              xyc=DBIJ*VALXY
              DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
              DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
              DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VALYY
              xzc=DBIJ*VALXZ
              DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
              DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
              yzc=DBIJ*VALYZ
              DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
              DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
              DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
c  alpha
              GDXXA(1,1,IAtm,1,IAtm)=GDXXA(1,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXA(1,1,IAtm,2,IAtm)=GDXXA(1,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,1,IAtm)=GDXXA(1,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,2,IAtm)=GDXXA(1,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXA(1,1,IAtm,3,IAtm)=GDXXA(1,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,3,IAtm,1,IAtm)=GDXXA(1,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,2,IAtm,3,IAtm)=GDXXA(1,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,2,IAtm)=GDXXA(1,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,3,IAtm)=GDXXA(1,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXA(2,1,IAtm,1,IAtm)=GDXXA(2,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXA(2,1,IAtm,2,IAtm)=GDXXA(2,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,1,IAtm)=GDXXA(2,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,2,IAtm)=GDXXA(2,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXA(2,1,IAtm,3,IAtm)=GDXXA(2,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,3,IAtm,1,IAtm)=GDXXA(2,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,2,IAtm,3,IAtm)=GDXXA(2,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,2,IAtm)=GDXXA(2,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,3,IAtm)=GDXXA(2,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXA(3,1,IAtm,1,IAtm)=GDXXA(3,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXA(3,1,IAtm,2,IAtm)=GDXXA(3,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,1,IAtm)=GDXXA(3,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,2,IAtm)=GDXXA(3,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXA(3,1,IAtm,3,IAtm)=GDXXA(3,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,3,IAtm,1,IAtm)=GDXXA(3,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,2,IAtm,3,IAtm)=GDXXA(3,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,2,IAtm)=GDXXA(3,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,3,IAtm)=GDXXA(3,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
c  beta
              GDXXB(1,1,IAtm,1,IAtm)=GDXXB(1,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXB(1,1,IAtm,2,IAtm)=GDXXB(1,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,1,IAtm)=GDXXB(1,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,2,IAtm)=GDXXB(1,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXB(1,1,IAtm,3,IAtm)=GDXXB(1,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,3,IAtm,1,IAtm)=GDXXB(1,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,2,IAtm,3,IAtm)=GDXXB(1,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,2,IAtm)=GDXXB(1,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,3,IAtm)=GDXXB(1,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXB(2,1,IAtm,1,IAtm)=GDXXB(2,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXB(2,1,IAtm,2,IAtm)=GDXXB(2,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,1,IAtm)=GDXXB(2,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,2,IAtm)=GDXXB(2,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXB(2,1,IAtm,3,IAtm)=GDXXB(2,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,3,IAtm,1,IAtm)=GDXXB(2,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,2,IAtm,3,IAtm)=GDXXB(2,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,2,IAtm)=GDXXB(2,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,3,IAtm)=GDXXB(2,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXB(3,1,IAtm,1,IAtm)=GDXXB(3,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXB(3,1,IAtm,2,IAtm)=GDXXB(3,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,1,IAtm)=GDXXB(3,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,2,IAtm)=GDXXB(3,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXB(3,1,IAtm,3,IAtm)=GDXXB(3,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,3,IAtm,1,IAtm)=GDXXB(3,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,2,IAtm,3,IAtm)=GDXXB(3,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,2,IAtm)=GDXXB(3,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,3,IAtm)=GDXXB(3,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
c
            If(JAtm.LT.IAtm) GO TO 218
            If(JAtm.eq.icntr) GO TO 218
c
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJ2*VALJX
              DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJY
              DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJZ
              DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
              DBIJ = DBIJ2*VALJX
              DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJY
              DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJZ
              DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
c    alpha
              GDXXA(1,1,IAtm,1,JAtm)=GDXXA(1,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXA(1,1,IAtm,2,JAtm)=GDXXA(1,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXA(1,2,IAtm,1,JAtm)=GDXXA(1,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXA(1,2,IAtm,2,JAtm)=GDXXA(1,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXA(1,1,IAtm,3,JAtm)=GDXXA(1,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXA(1,3,IAtm,1,JAtm)=GDXXA(1,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXA(1,2,IAtm,3,JAtm)=GDXXA(1,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXA(1,3,IAtm,2,JAtm)=GDXXA(1,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXA(1,3,IAtm,3,JAtm)=GDXXA(1,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXA(2,1,IAtm,1,JAtm)=GDXXA(2,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXA(2,1,IAtm,2,JAtm)=GDXXA(2,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXA(2,2,IAtm,1,JAtm)=GDXXA(2,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXA(2,2,IAtm,2,JAtm)=GDXXA(2,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXA(2,1,IAtm,3,JAtm)=GDXXA(2,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXA(2,3,IAtm,1,JAtm)=GDXXA(2,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXA(2,2,IAtm,3,JAtm)=GDXXA(2,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXA(2,3,IAtm,2,JAtm)=GDXXA(2,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXA(2,3,IAtm,3,JAtm)=GDXXA(2,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXA(3,1,IAtm,1,JAtm)=GDXXA(3,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXA(3,1,IAtm,2,JAtm)=GDXXA(3,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXA(3,2,IAtm,1,JAtm)=GDXXA(3,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXA(3,2,IAtm,2,JAtm)=GDXXA(3,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXA(3,1,IAtm,3,JAtm)=GDXXA(3,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXA(3,3,IAtm,1,JAtm)=GDXXA(3,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXA(3,2,IAtm,3,JAtm)=GDXXA(3,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXA(3,3,IAtm,2,JAtm)=GDXXA(3,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXA(3,3,IAtm,3,JAtm)=GDXXA(3,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
c    beta
              GDXXB(1,1,IAtm,1,JAtm)=GDXXB(1,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXB(1,1,IAtm,2,JAtm)=GDXXB(1,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXB(1,2,IAtm,1,JAtm)=GDXXB(1,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXB(1,2,IAtm,2,JAtm)=GDXXB(1,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXB(1,1,IAtm,3,JAtm)=GDXXB(1,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXB(1,3,IAtm,1,JAtm)=GDXXB(1,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXB(1,2,IAtm,3,JAtm)=GDXXB(1,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXB(1,3,IAtm,2,JAtm)=GDXXB(1,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXB(1,3,IAtm,3,JAtm)=GDXXB(1,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXB(2,1,IAtm,1,JAtm)=GDXXB(2,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXB(2,1,IAtm,2,JAtm)=GDXXB(2,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXB(2,2,IAtm,1,JAtm)=GDXXB(2,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXB(2,2,IAtm,2,JAtm)=GDXXB(2,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXB(2,1,IAtm,3,JAtm)=GDXXB(2,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXB(2,3,IAtm,1,JAtm)=GDXXB(2,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXB(2,2,IAtm,3,JAtm)=GDXXB(2,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXB(2,3,IAtm,2,JAtm)=GDXXB(2,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXB(2,3,IAtm,3,JAtm)=GDXXB(2,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXB(3,1,IAtm,1,JAtm)=GDXXB(3,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXB(3,1,IAtm,2,JAtm)=GDXXB(3,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXB(3,2,IAtm,1,JAtm)=GDXXB(3,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXB(3,2,IAtm,2,JAtm)=GDXXB(3,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXB(3,1,IAtm,3,JAtm)=GDXXB(3,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXB(3,3,IAtm,1,JAtm)=GDXXB(3,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXB(3,2,IAtm,3,JAtm)=GDXXB(3,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXB(3,3,IAtm,2,JAtm)=GDXXB(3,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXB(3,3,IAtm,3,JAtm)=GDXXB(3,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2= two*DA(IJ)
            DBIJ2= two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
              DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
              DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
              DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
              DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
              DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VALZZ
c
              DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VALXX
              xyc=DBIJ*VALXY
              DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
              DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
              DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VALYY
              xzc=DBIJ*VALXZ
              DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
              DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
              yzc=DBIJ*VALYZ
              DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
              DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
              DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
c    alpha
              GDXXA(1,1,IAtm,1,IAtm)=GDXXA(1,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXA(1,1,IAtm,2,IAtm)=GDXXA(1,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,1,IAtm)=GDXXA(1,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,2,IAtm)=GDXXA(1,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXA(1,1,IAtm,3,IAtm)=GDXXA(1,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,3,IAtm,1,IAtm)=GDXXA(1,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,2,IAtm,3,IAtm)=GDXXA(1,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,2,IAtm)=GDXXA(1,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,3,IAtm)=GDXXA(1,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXA(2,1,IAtm,1,IAtm)=GDXXA(2,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXA(2,1,IAtm,2,IAtm)=GDXXA(2,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,1,IAtm)=GDXXA(2,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,2,IAtm)=GDXXA(2,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXA(2,1,IAtm,3,IAtm)=GDXXA(2,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,3,IAtm,1,IAtm)=GDXXA(2,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,2,IAtm,3,IAtm)=GDXXA(2,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,2,IAtm)=GDXXA(2,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,3,IAtm)=GDXXA(2,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXA(3,1,IAtm,1,IAtm)=GDXXA(3,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXA(3,1,IAtm,2,IAtm)=GDXXA(3,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,1,IAtm)=GDXXA(3,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,2,IAtm)=GDXXA(3,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXA(3,1,IAtm,3,IAtm)=GDXXA(3,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,3,IAtm,1,IAtm)=GDXXA(3,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,2,IAtm,3,IAtm)=GDXXA(3,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,2,IAtm)=GDXXA(3,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,3,IAtm)=GDXXA(3,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
c    beta
              GDXXB(1,1,IAtm,1,IAtm)=GDXXB(1,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXB(1,1,IAtm,2,IAtm)=GDXXB(1,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,1,IAtm)=GDXXB(1,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,2,IAtm)=GDXXB(1,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXB(1,1,IAtm,3,IAtm)=GDXXB(1,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,3,IAtm,1,IAtm)=GDXXB(1,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,2,IAtm,3,IAtm)=GDXXB(1,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,2,IAtm)=GDXXB(1,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,3,IAtm)=GDXXB(1,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXB(2,1,IAtm,1,IAtm)=GDXXB(2,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXB(2,1,IAtm,2,IAtm)=GDXXB(2,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,1,IAtm)=GDXXB(2,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,2,IAtm)=GDXXB(2,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXB(2,1,IAtm,3,IAtm)=GDXXB(2,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,3,IAtm,1,IAtm)=GDXXB(2,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,2,IAtm,3,IAtm)=GDXXB(2,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,2,IAtm)=GDXXB(2,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,3,IAtm)=GDXXB(2,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXB(3,1,IAtm,1,IAtm)=GDXXB(3,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXB(3,1,IAtm,2,IAtm)=GDXXB(3,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,1,IAtm)=GDXXB(3,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,2,IAtm)=GDXXB(3,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXB(3,1,IAtm,3,IAtm)=GDXXB(3,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,3,IAtm,1,IAtm)=GDXXB(3,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,2,IAtm,3,IAtm)=GDXXB(3,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,2,IAtm)=GDXXB(3,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,3,IAtm)=GDXXB(3,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
            If(JAtm.LT.IAtm) GO TO 219
            If(JAtm.EQ.ICntr) GO TO 219
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJ2*VALJX
              DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJY
              DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJZ
              DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
              DBIJ = DBIJ2*VALJX
              DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJY
              DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJZ
              DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
c     alpha
              GDXXA(1,1,IAtm,1,JAtm)=GDXXA(1,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXA(1,1,IAtm,2,JAtm)=GDXXA(1,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXA(1,2,IAtm,1,JAtm)=GDXXA(1,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXA(1,2,IAtm,2,JAtm)=GDXXA(1,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXA(1,1,IAtm,3,JAtm)=GDXXA(1,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXA(1,3,IAtm,1,JAtm)=GDXXA(1,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXA(1,2,IAtm,3,JAtm)=GDXXA(1,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXA(1,3,IAtm,2,JAtm)=GDXXA(1,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXA(1,3,IAtm,3,JAtm)=GDXXA(1,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXA(2,1,IAtm,1,JAtm)=GDXXA(2,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXA(2,1,IAtm,2,JAtm)=GDXXA(2,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXA(2,2,IAtm,1,JAtm)=GDXXA(2,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXA(2,2,IAtm,2,JAtm)=GDXXA(2,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXA(2,1,IAtm,3,JAtm)=GDXXA(2,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXA(2,3,IAtm,1,JAtm)=GDXXA(2,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXA(2,2,IAtm,3,JAtm)=GDXXA(2,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXA(2,3,IAtm,2,JAtm)=GDXXA(2,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXA(2,3,IAtm,3,JAtm)=GDXXA(2,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXA(3,1,IAtm,1,JAtm)=GDXXA(3,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXA(3,1,IAtm,2,JAtm)=GDXXA(3,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXA(3,2,IAtm,1,JAtm)=GDXXA(3,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXA(3,2,IAtm,2,JAtm)=GDXXA(3,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXA(3,1,IAtm,3,JAtm)=GDXXA(3,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXA(3,3,IAtm,1,JAtm)=GDXXA(3,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXA(3,2,IAtm,3,JAtm)=GDXXA(3,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXA(3,3,IAtm,2,JAtm)=GDXXA(3,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXA(3,3,IAtm,3,JAtm)=GDXXA(3,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
c     beta
              GDXXB(1,1,IAtm,1,JAtm)=GDXXB(1,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXB(1,1,IAtm,2,JAtm)=GDXXB(1,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXB(1,2,IAtm,1,JAtm)=GDXXB(1,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXB(1,2,IAtm,2,JAtm)=GDXXB(1,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXB(1,1,IAtm,3,JAtm)=GDXXB(1,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXB(1,3,IAtm,1,JAtm)=GDXXB(1,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXB(1,2,IAtm,3,JAtm)=GDXXB(1,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXB(1,3,IAtm,2,JAtm)=GDXXB(1,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXB(1,3,IAtm,3,JAtm)=GDXXB(1,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXB(2,1,IAtm,1,JAtm)=GDXXB(2,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXB(2,1,IAtm,2,JAtm)=GDXXB(2,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXB(2,2,IAtm,1,JAtm)=GDXXB(2,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXB(2,2,IAtm,2,JAtm)=GDXXB(2,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXB(2,1,IAtm,3,JAtm)=GDXXB(2,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXB(2,3,IAtm,1,JAtm)=GDXXB(2,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXB(2,2,IAtm,3,JAtm)=GDXXB(2,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXB(2,3,IAtm,2,JAtm)=GDXXB(2,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXB(2,3,IAtm,3,JAtm)=GDXXB(2,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXB(3,1,IAtm,1,JAtm)=GDXXB(3,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXB(3,1,IAtm,2,JAtm)=GDXXB(3,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXB(3,2,IAtm,1,JAtm)=GDXXB(3,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXB(3,2,IAtm,2,JAtm)=GDXXB(3,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXB(3,1,IAtm,3,JAtm)=GDXXB(3,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXB(3,3,IAtm,1,JAtm)=GDXXB(3,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXB(3,2,IAtm,3,JAtm)=GDXXB(3,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXB(3,3,IAtm,2,JAtm)=GDXXB(3,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXB(3,3,IAtm,3,JAtm)=GDXXB(3,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and Hessian of density gradient invariant
c
c
        DO IAtm=1,NAtoms
c   alpha alpha
          GXAA(1,IAtm)=two*(DAX*GDXA(1,1,IAtm)+
     $                   DAY*GDXA(2,1,IAtm)+DAZ*GDXA(3,1,IAtm))
          GXAA(2,IAtm)=two*(DAX*GDXA(1,2,IAtm)+
     $                   DAY*GDXA(2,2,IAtm)+DAZ*GDXA(3,2,IAtm))
          GXAA(3,IAtm)=two*(DAX*GDXA(1,3,IAtm)+
     $                   DAY*GDXA(2,3,IAtm)+DAZ*GDXA(3,3,IAtm))
c   beta beta
          GXBB(1,IAtm)=two*(DBX*GDXB(1,1,IAtm)+
     $                   DBY*GDXB(2,1,IAtm)+DBZ*GDXB(3,1,IAtm))
          GXBB(2,IAtm)=two*(DBX*GDXB(1,2,IAtm)+
     $                   DBY*GDXB(2,2,IAtm)+DBZ*GDXB(3,2,IAtm))
          GXBB(3,IAtm)=two*(DBX*GDXB(1,3,IAtm)+
     $                   DBY*GDXB(2,3,IAtm)+DBZ*GDXB(3,3,IAtm))
c   alpha beta
          GXAB(1,IAtm)=DAX*GDXB(1,1,IAtm)+DBX*GDXA(1,1,IAtm)+
     $                 DAY*GDXB(2,1,IAtm)+DBY*GDXA(2,1,IAtm)+
     $                 DAZ*GDXB(3,1,IAtm)+DBZ*GDXA(3,1,IAtm)
          GXAB(2,IAtm)=DAX*GDXB(1,2,IAtm)+DBX*GDXA(1,2,IAtm)+
     $                 DAY*GDXB(2,2,IAtm)+DBY*GDXA(2,2,IAtm)+
     $                 DAZ*GDXB(3,2,IAtm)+DBZ*GDXA(3,2,IAtm)
          GXAB(3,IAtm)=DAX*GDXB(1,3,IAtm)+DBX*GDXA(1,3,IAtm)+
     $                 DAY*GDXB(2,3,IAtm)+DBY*GDXA(2,3,IAtm)+
     $                 DAZ*GDXB(3,3,IAtm)+DBZ*GDXA(3,3,IAtm)
        EndDO
c
c  compute second derivatives of gradient invariants
c  only of needed
c
        if(dogradxx)then
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
c     alpha alpha
            GXXAA(1,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,1,JAtm)+GDXA(1,1,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,1,JAtm)+GDXA(2,1,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,1,JAtm)+GDXA(3,1,IAtm)*GDXA(3,1,JAtm))
            GXXAA(1,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,2,JAtm)+GDXA(1,1,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,2,JAtm)+GDXA(2,1,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,2,JAtm)+GDXA(3,1,IAtm)*GDXA(3,2,JAtm))
            GXXAA(2,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,1,JAtm)+GDXA(1,2,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,1,JAtm)+GDXA(2,2,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,1,JAtm)+GDXA(3,2,IAtm)*GDXA(3,1,JAtm))
            GXXAA(2,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,2,JAtm)+GDXA(1,2,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,2,JAtm)+GDXA(2,2,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,2,JAtm)+GDXA(3,2,IAtm)*GDXA(3,2,JAtm))
            GXXAA(1,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,3,JAtm)+GDXA(1,1,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,3,JAtm)+GDXA(2,1,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,3,JAtm)+GDXA(3,1,IAtm)*GDXA(3,3,JAtm))
            GXXAA(3,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,1,JAtm)+GDXA(1,3,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,1,JAtm)+GDXA(2,3,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,1,JAtm)+GDXA(3,3,IAtm)*GDXA(3,1,JAtm))
            GXXAA(2,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,3,JAtm)+GDXA(1,2,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,3,JAtm)+GDXA(2,2,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,3,JAtm)+GDXA(3,2,IAtm)*GDXA(3,3,JAtm))
            GXXAA(3,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,2,JAtm)+GDXA(1,3,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,2,JAtm)+GDXA(2,3,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,2,JAtm)+GDXA(3,3,IAtm)*GDXA(3,2,JAtm))
            GXXAA(3,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,3,JAtm)+GDXA(1,3,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,3,JAtm)+GDXA(2,3,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,3,JAtm)+GDXA(3,3,IAtm)*GDXA(3,3,JAtm))
c     beta beta
            GXXBB(1,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,1,JAtm)+GDXB(1,1,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,1,JAtm)+GDXB(2,1,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,1,JAtm)+GDXB(3,1,IAtm)*GDXB(3,1,JAtm))
            GXXBB(1,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,2,JAtm)+GDXB(1,1,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,2,JAtm)+GDXB(2,1,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,2,JAtm)+GDXB(3,1,IAtm)*GDXB(3,2,JAtm))
            GXXBB(2,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,1,JAtm)+GDXB(1,2,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,1,JAtm)+GDXB(2,2,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,1,JAtm)+GDXB(3,2,IAtm)*GDXB(3,1,JAtm))
            GXXBB(2,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,2,JAtm)+GDXB(1,2,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,2,JAtm)+GDXB(2,2,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,2,JAtm)+GDXB(3,2,IAtm)*GDXB(3,2,JAtm))
            GXXBB(1,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,3,JAtm)+GDXB(1,1,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,3,JAtm)+GDXB(2,1,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,3,JAtm)+GDXB(3,1,IAtm)*GDXB(3,3,JAtm))
            GXXBB(3,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,1,JAtm)+GDXB(1,3,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,1,JAtm)+GDXB(2,3,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,1,JAtm)+GDXB(3,3,IAtm)*GDXB(3,1,JAtm))
            GXXBB(2,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,3,JAtm)+GDXB(1,2,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,3,JAtm)+GDXB(2,2,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,3,JAtm)+GDXB(3,2,IAtm)*GDXB(3,3,JAtm))
            GXXBB(3,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,2,JAtm)+GDXB(1,3,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,2,JAtm)+GDXB(2,3,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,2,JAtm)+GDXB(3,3,IAtm)*GDXB(3,2,JAtm))
            GXXBB(3,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,3,JAtm)+GDXB(1,3,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,3,JAtm)+GDXB(2,3,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,3,JAtm)+GDXB(3,3,IAtm)*GDXB(3,3,JAtm))
          EndDO
        EndDO
c
c     alpha beta is done separately, as only few functionals
c     actually need this one
c
        if(abs(gc).gt.epsi)then
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
            GXXAB(1,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,1,JAtm)+GDXA(1,1,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,1,JAtm)+GDXB(1,1,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,1,JAtm)+GDXA(2,1,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,1,JAtm)+GDXB(2,1,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,1,JAtm)+GDXA(3,1,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,1,JAtm)+GDXB(3,1,IAtm)*GDXA(3,1,JAtm)
            GXXAB(2,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,1,JAtm)+GDXA(1,2,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,1,JAtm)+GDXB(1,2,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,1,JAtm)+GDXA(2,2,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,1,JAtm)+GDXB(2,2,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,1,JAtm)+GDXA(3,2,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,1,JAtm)+GDXB(3,2,IAtm)*GDXA(3,1,JAtm)
            GXXAB(3,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,1,JAtm)+GDXA(1,3,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,1,JAtm)+GDXB(1,3,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,1,JAtm)+GDXA(2,3,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,1,JAtm)+GDXB(2,3,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,1,JAtm)+GDXA(3,3,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,1,JAtm)+GDXB(3,3,IAtm)*GDXA(3,1,JAtm)
            GXXAB(1,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,2,JAtm)+GDXA(1,1,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,2,JAtm)+GDXB(1,1,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,2,JAtm)+GDXA(2,1,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,2,JAtm)+GDXB(2,1,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,2,JAtm)+GDXA(3,1,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,2,JAtm)+GDXB(3,1,IAtm)*GDXA(3,2,JAtm)
            GXXAB(2,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,2,JAtm)+GDXA(1,2,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,2,JAtm)+GDXB(1,2,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,2,JAtm)+GDXA(2,2,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,2,JAtm)+GDXB(2,2,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,2,JAtm)+GDXA(3,2,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,2,JAtm)+GDXB(3,2,IAtm)*GDXA(3,2,JAtm)
            GXXAB(3,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,2,JAtm)+GDXA(1,3,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,2,JAtm)+GDXB(1,3,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,2,JAtm)+GDXA(2,3,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,2,JAtm)+GDXB(2,3,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,2,JAtm)+GDXA(3,3,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,2,JAtm)+GDXB(3,3,IAtm)*GDXA(3,2,JAtm)
            GXXAB(1,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,3,JAtm)+GDXA(1,1,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,3,JAtm)+GDXB(1,1,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,3,JAtm)+GDXA(2,1,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,3,JAtm)+GDXB(2,1,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,3,JAtm)+GDXA(3,1,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,3,JAtm)+GDXB(3,1,IAtm)*GDXA(3,3,JAtm)
            GXXAB(2,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,3,JAtm)+GDXA(1,2,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,3,JAtm)+GDXB(1,2,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,3,JAtm)+GDXA(2,2,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,3,JAtm)+GDXB(2,2,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,3,JAtm)+GDXA(3,2,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,3,JAtm)+GDXB(3,2,IAtm)*GDXA(3,3,JAtm)
            GXXAB(3,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,3,JAtm)+GDXA(1,3,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,3,JAtm)+GDXB(1,3,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,3,JAtm)+GDXA(2,3,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,3,JAtm)+GDXB(2,3,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,3,JAtm)+GDXA(3,3,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,3,JAtm)+GDXB(3,3,IAtm)*GDXA(3,3,JAtm)
          EndDO
        EndDO
        else
          CALL ZeroIT(GXXAB,NAt3**2)
        endif
        else
          CALL ZeroIT(GXXAA,NAt3**2)
          CALL ZeroIT(GXXBB,NAt3**2)
          CALL ZeroIT(GXXAB,NAt3**2)
        endif
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
C    and also compute partial sums used in Hessian quadrature
c
        ga2= ga+ga
        gb2= gb+gb
        DAX2=DAX+DAX
        DAY2=DAY+DAY
        DAZ2=DAZ+DAZ
        DBX2=DBX+DBX
        DBY2=DBY+DBY
        DBZ2=DBZ+DBZ
c
c     coeff. X
c
        SXXA=ga2*DAX+gc*DBX
        SXYA=ga2*DAY+gc*DBY
        SXZA=ga2*DAZ+gc*DBZ
        SXXB=gb2*DBX+gc*DAX
        SXYB=gb2*DBY+gc*DAY
        SXZB=gb2*DBZ+gc*DAZ
c
        DO IAtm=1,NAtoms
c
c     coeff. V
c
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)+
     $        raga*gxaa(1,iatm)+ragb*gxbb(1,iatm)+ragc*gxab(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)+
     $        raga*gxaa(2,iatm)+ragb*gxbb(2,iatm)+ragc*gxab(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)+
     $        raga*gxaa(3,iatm)+ragb*gxbb(3,iatm)+ragc*gxab(3,iatm)
          svb(1,iatm)=rarb*denxa(1,iatm)+rbrb*denxb(1,iatm)+
     $        rbga*gxaa(1,iatm)+rbgb*gxbb(1,iatm)+rbgc*gxab(1,iatm)
          svb(2,iatm)=rarb*denxa(2,iatm)+rbrb*denxb(2,iatm)+
     $        rbga*gxaa(2,iatm)+rbgb*gxbb(2,iatm)+rbgc*gxab(2,iatm)
          svb(3,iatm)=rarb*denxa(3,iatm)+rbrb*denxb(3,iatm)+
     $        rbga*gxaa(3,iatm)+rbgb*gxbb(3,iatm)+rbgc*gxab(3,iatm)
c
c  partial sums for Hessian
c
          sga(1,iatm)=raga*denxa(1,iatm)+rbga*denxb(1,iatm)+
     $        gaga*gxaa(1,iatm)+gagb*gxbb(1,iatm)+gagc*gxab(1,iatm)
          sga(2,iatm)=raga*denxa(2,iatm)+rbga*denxb(2,iatm)+
     $        gaga*gxaa(2,iatm)+gagb*gxbb(2,iatm)+gagc*gxab(2,iatm)
          sga(3,iatm)=raga*denxa(3,iatm)+rbga*denxb(3,iatm)+
     $        gaga*gxaa(3,iatm)+gagb*gxbb(3,iatm)+gagc*gxab(3,iatm)
c
          sgb(1,iatm)=ragb*denxa(1,iatm)+rbgb*denxb(1,iatm)+
     $        gagb*gxaa(1,iatm)+gbgb*gxbb(1,iatm)+gbgc*gxab(1,iatm)
          sgb(2,iatm)=ragb*denxa(2,iatm)+rbgb*denxb(2,iatm)+
     $        gagb*gxaa(2,iatm)+gbgb*gxbb(2,iatm)+gbgc*gxab(2,iatm)
          sgb(3,iatm)=ragb*denxa(3,iatm)+rbgb*denxb(3,iatm)+
     $        gagb*gxaa(3,iatm)+gbgb*gxbb(3,iatm)+gbgc*gxab(3,iatm)
c
          sgc(1,iatm)=ragc*denxa(1,iatm)+rbgc*denxb(1,iatm)+
     $        gagc*gxaa(1,iatm)+gbgc*gxbb(1,iatm)+gcgc*gxab(1,iatm)
          sgc(2,iatm)=ragc*denxa(2,iatm)+rbgc*denxb(2,iatm)+
     $        gagc*gxaa(2,iatm)+gbgc*gxbb(2,iatm)+gcgc*gxab(2,iatm)
          sgc(3,iatm)=ragc*denxa(3,iatm)+rbgc*denxb(3,iatm)+
     $        gagc*gxaa(3,iatm)+gbgc*gxbb(3,iatm)+gcgc*gxab(3,iatm)
c
c     coeff. W
c
          swa(1,1,iatm)=ga2*gdxa(1,1,iatm)+gc*gdxb(1,1,Iatm)+
     $        DAX2*sga(1,iatm)+DBX*sgc(1,iatm)
          swa(1,2,iatm)=ga2*gdxa(1,2,iatm)+gc*gdxb(1,2,Iatm)+
     $        DAX2*sga(2,iatm)+DBX*sgc(2,iatm)
          swa(1,3,iatm)=ga2*gdxa(1,3,iatm)+gc*gdxb(1,3,Iatm)+
     $        DAX2*sga(3,iatm)+DBX*sgc(3,iatm)
          swa(2,1,iatm)=ga2*gdxa(2,1,iatm)+gc*gdxb(2,1,Iatm)+
     $        DAY2*sga(1,iatm)+DBY*sgc(1,iatm)
          swa(2,2,iatm)=ga2*gdxa(2,2,iatm)+gc*gdxb(2,2,Iatm)+
     $        DAY2*sga(2,iatm)+DBY*sgc(2,iatm)
          swa(2,3,iatm)=ga2*gdxa(2,3,iatm)+gc*gdxb(2,3,Iatm)+
     $        DAY2*sga(3,iatm)+DBY*sgc(3,iatm)
          swa(3,1,iatm)=ga2*gdxa(3,1,iatm)+gc*gdxb(3,1,Iatm)+
     $        DAZ2*sga(1,iatm)+DBZ*sgc(1,iatm)
          swa(3,2,iatm)=ga2*gdxa(3,2,iatm)+gc*gdxb(3,2,Iatm)+
     $        DAZ2*sga(2,iatm)+DBZ*sgc(2,iatm)
          swa(3,3,iatm)=ga2*gdxa(3,3,iatm)+gc*gdxb(3,3,Iatm)+
     $        DAZ2*sga(3,iatm)+DBZ*sgc(3,iatm)
c
          swb(1,1,iatm)=gb2*gdxb(1,1,iatm)+gc*gdxa(1,1,Iatm)+
     $        DBX2*sgb(1,iatm)+DAX*sgc(1,iatm)
          swb(1,2,iatm)=gb2*gdxb(1,2,iatm)+gc*gdxa(1,2,Iatm)+
     $        DBX2*sgb(2,iatm)+DAX*sgc(2,iatm)
          swb(1,3,iatm)=gb2*gdxb(1,3,iatm)+gc*gdxa(1,3,Iatm)+
     $        DBX2*sgb(3,iatm)+DAX*sgc(3,iatm)
          swb(2,1,iatm)=gb2*gdxb(2,1,iatm)+gc*gdxa(2,1,Iatm)+
     $        DBY2*sgb(1,iatm)+DAY*sgc(1,iatm)
          swb(2,2,iatm)=gb2*gdxb(2,2,iatm)+gc*gdxa(2,2,Iatm)+
     $        DBY2*sgb(2,iatm)+DAY*sgc(2,iatm)
          swb(2,3,iatm)=gb2*gdxb(2,3,iatm)+gc*gdxa(2,3,Iatm)+
     $        DBY2*sgb(3,iatm)+DAY*sgc(3,iatm)
          swb(3,1,iatm)=gb2*gdxb(3,1,iatm)+gc*gdxa(3,1,Iatm)+
     $        DBZ2*sgb(1,iatm)+DAZ*sgc(1,iatm)
          swb(3,2,iatm)=gb2*gdxb(3,2,iatm)+gc*gdxa(3,2,Iatm)+
     $        DBZ2*sgb(2,iatm)+DAZ*sgc(2,iatm)
          swb(3,3,iatm)=gb2*gdxb(3,3,iatm)+gc*gdxa(3,3,Iatm)+
     $        DBZ2*sgb(3,iatm)+DAZ*sgc(3,iatm)
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    fock matrix.
c
c -- get the maximum absolute value of coefficients v and w
c
        Call absmax(NAt3,sva,iixx,svamax)
        Call absmax(NAt3,svb,iixx,svbmax)
        svmax=max(svamax,svbmax)
        Call absmax(3*NAt3,swa,iixx,swamax)
        Call absmax(3*NAt3,swb,iixx,swbmax)
        swmax=max(swamax,swbmax)
        tswmax=Three*swmax
c
c -- global threshold testing
c
        abr=max(abs(ra),abs(rb))
        vmx2=vmx*vmx
        SMax = max(Abs(SXXA+SXYA+SXZA),abs(SXXB+SXYB+SXZB))
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abr+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVIA=SXXA*ValIX+SXYA*ValIY+SXZA*ValIZ+ra*ValI
          XVIB=SXXB*ValIX+SXYB*ValIY+SXZB*ValIZ+rb*ValI
          SXXIA=SXXA*ValI
          SXYIA=SXYA*ValI
          SXZIA=SXZA*ValI
          SXXIB=SXXB*ValI
          SXYIB=SXYB*ValI
          SXZIB=SXZB*ValI
          XXVXA=SXXA*VAOXX(1,I)+SXYA*VAOXX(2,I)+SXZA*VAOXX(4,I)+ra*VALIX
          XXVYA=SXXA*VAOXX(2,I)+SXYA*VAOXX(3,I)+SXZA*VAOXX(5,I)+ra*VALIY
          XXVZA=SXXA*VAOXX(4,I)+SXYA*VAOXX(5,I)+SXZA*VAOXX(6,I)+ra*VALIZ
          XXVXB=SXXB*VAOXX(1,I)+SXYB*VAOXX(2,I)+SXZB*VAOXX(4,I)+rb*VALIX
          XXVYB=SXXB*VAOXX(2,I)+SXYB*VAOXX(3,I)+SXZB*VAOXX(5,I)+rb*VALIY
          XXVZB=SXXB*VAOXX(4,I)+SXYB*VAOXX(5,I)+SXZB*VAOXX(6,I)+rb*VALIZ
          valm=abs(vali)
          valm1=max(abs(xvia),abs(xvib))
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1a=max(abs(xxvxa),abs(xxvya),abs(xxvza))
          valmx1b=max(abs(xxvxb),abs(xxvyb),abs(xxvzb))
          valmx1=max(valmx1a,valmx1b)
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              call zeroit(trpa,nat3)
              call zeroit(trpb,nat3)
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJA=SXXA*ValJX+SXYA*ValJY+SXZA*ValJZ
              XVJB=SXXB*ValJX+SXYB*ValJY+SXZB*ValJZ
              abxvj=max(abs(xvja),abs(xvjb))
              if(valmx1*abs(valj)+valmx*abxvj.gt.thrsh)then
              if(iatm.ne.icntr) then
                trpa(1,iatm)=trpa(1,iatm)-xxvxa*valj-valix*xvja
                trpa(2,iatm)=trpa(2,iatm)-xxvya*valj-valiy*xvja
                trpa(3,iatm)=trpa(3,iatm)-xxvza*valj-valiz*xvja
                trpb(1,iatm)=trpb(1,iatm)-xxvxb*valj-valix*xvjb
                trpb(2,iatm)=trpb(2,iatm)-xxvyb*valj-valiy*xvjb
                trpb(3,iatm)=trpb(3,iatm)-xxvzb*valj-valiz*xvjb
              endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
              if(jatm.ne.icntr) then
                trpa(1,jAtm) = trpa(1,jatm) - xvia*Valjx -
     $          vaoxx(1,j)*sxxia - vaoxx(2,j)*sxyia - vaoxx(4,j)*sxzia
                trpa(2,jatm) = trpa(2,jatm) - xvia*valjy -
     $          vaoxx(2,j)*sxxia - vaoxx(3,j)*sxyia - vaoxx(5,j)*sxzia
                trpa(3,jatm) = trpa(3,jAtm) - xvia*valjz -
     $          vaoxx(4,j)*sxxia - vaoxx(5,j)*sxyia - vaoxx(6,j)*sxzia
c
                trpb(1,jatm) = trpb(1,jatm) - xvib*valjx -
     $          vaoxx(1,j)*sxxib - vaoxx(2,j)*sxyib - vaoxx(4,j)*sxzib
                trpb(2,jatm) = trpb(2,jatm) - xvib*valjy -
     $          vaoxx(2,j)*sxxib - vaoxx(3,j)*sxyib - vaoxx(5,j)*sxzib
                trpb(3,jatm) = trpb(3,jatm) - xvib*valjz -
     $          vaoxx(4,j)*sxxib - vaoxx(5,j)*sxyib - vaoxx(6,j)*sxzib
              endif
              endif
              vij=vali*valj
              vijx=valix*valj+valjx*vali
              vijy=valiy*valj+valjy*vali
              vijz=valiz*valj+valjz*vali
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                do katm=1,natoms
                  trpa(1,katm) = trpa(1,katm) + (sva(1,katm)*vij+
     $                           swa(1,1,katm)*vijx+swa(2,1,katm)*vijy+
     $                           swa(3,1,katm)*vijz)
                  trpa(2,katm) = trpa(2,katm) + (sva(2,katm)*vij+
     $                           swa(1,2,katm)*vijx+swa(2,2,katm)*vijy+
     $                           swa(3,2,katm)*vijz)
                  trpa(3,katm) = trpa(3,katm) + (sva(3,katm)*vij+
     $                           swa(1,3,katm)*vijx+swa(2,3,katm)*vijy+
     $                           swa(3,3,katm)*vijz)
c
                  trpb(1,katm) = trpb(1,katm) + (svb(1,katm)*vij+
     $                           swb(1,1,katm)*vijx+swb(2,1,katm)*vijy+
     $                           swb(3,1,katm)*vijz)
                  trpb(2,katm) = trpb(2,katm) + (svb(2,katm)*vij+
     $                           swb(1,2,katm)*vijx+swb(2,2,katm)*vijy+
     $                           swb(3,2,katm)*vijz)
                  trpb(3,katm) = trpb(3,katm) + (svb(3,katm)*vij+
     $                           swb(1,3,katm)*vijx+swb(2,3,katm)*vijy+
     $                           swb(3,3,katm)*vijz)
                enddo
              endif
c
c -- weight derivatives contribution
c
              dxija=dax*vijx+day*vijy+daz*vijz
              dxijb=dbx*vijx+dby*vijy+dbz*vijz
              gija=prap*vij+pgap2*dxija+pgcp*dxijb
              gijb=prbp*vij+pgbp2*dxijb+pgcp*dxija
              if(gwtmax*max(abs(gija),abs(gijb)).gt.thrsh)then
                do ia=1,natoms
                  if(ia.ne.icntr)then
                    trpa(1,ia)=trpa(1,ia)+gwt(1,ia,ipp)*gija
                    trpa(2,ia)=trpa(2,ia)+gwt(2,ia,ipp)*gija
                    trpa(3,ia)=trpa(3,ia)+gwt(3,ia,ipp)*gija
                    trpb(1,ia)=trpb(1,ia)+gwt(1,ia,ipp)*gijb
                    trpb(2,ia)=trpb(2,ia)+gwt(2,ia,ipp)*gijb
                    trpb(3,ia)=trpb(3,ia)+gwt(3,ia,ipp)*gijb
                  endif
                enddo
              endif
c
c -- apply translational invariance and accumulate contributions
c    into the fda array
c
              do ia=1,natoms
                if(ia.ne.icntr)then
                  trpa(1,icntr)=trpa(1,icntr)-trpa(1,ia)
                  trpa(2,icntr)=trpa(2,icntr)-trpa(2,ia)
                  trpa(3,icntr)=trpa(3,icntr)-trpa(3,ia)
                  trpb(1,icntr)=trpb(1,icntr)-trpb(1,ia)
                  trpb(2,icntr)=trpb(2,icntr)-trpb(2,ia)
                  trpb(3,icntr)=trpb(3,icntr)-trpb(3,ia)
                endif
              enddo
              do ia=nb,ne
                fda(1,ia,ij)=fda(1,ia,ij)+trpa(1,ia)
                fda(2,ia,ij)=fda(2,ia,ij)+trpa(2,ia)
                fda(3,ia,ij)=fda(3,ia,ij)+trpa(3,ia)
                fdb(1,ia,ij)=fdb(1,ia,ij)+trpb(1,ia)
                fdb(2,ia,ij)=fdb(2,ia,ij)+trpb(2,ia)
                fdb(3,ia,ij)=fdb(3,ia,ij)+trpb(3,ia)
              enddo
 230        continue
          endif
 240    continue
c
c   numerical quadrature for direct contribution to hessian matrix
c
 245    continue
        call secund(t5)
        tqf=tqf+t5-t4
c
c -- direct contribution to hessian matrix.
c
        Do IAtm=1,NAtoms
          Do JAtm=IAtm,NAtoms
            hess(1,IAtm,1,JAtm) = hess(1,IAtm,1,JAtm) +
     $       sva(1,IAtm)*DenXA(1,JAtm) + svb(1,IAtm)*DenXB(1,JAtm) +
     $       sga(1,IAtm)*gxaa(1,JAtm)  + sgb(1,IAtm)*gxbb(1,JAtm) +
     $       sgc(1,IAtm)*gxab(1,JAtm)  + ra*DenXXA(1,IAtm,1,JAtm) +
     $       rb*DenXXB(1,IAtm,1,JAtm)  + ga*gxxaa(1,IAtm,1,JAtm) +
     $       gb*gxxbb(1,IAtm,1,JAtm)   + gc*gxxab(1,IAtm,1,JAtm)
            hess(2,IAtm,1,JAtm) = hess(2,IAtm,1,JAtm) +
     $       sva(2,IAtm)*DenXA(1,JAtm) + svb(2,IAtm)*DenXB(1,JAtm) +
     $       sga(2,IAtm)*gxaa(1,JAtm)  + sgb(2,IAtm)*gxbb(1,JAtm) +
     $       sgc(2,IAtm)*gxab(1,JAtm)  + ra*DenXXA(2,IAtm,1,JAtm) +
     $       rb*DenXXB(2,IAtm,1,JAtm)  + ga*gxxaa(2,IAtm,1,JAtm) +
     $       gb*gxxbb(2,IAtm,1,JAtm)   + gc*gxxab(2,IAtm,1,JAtm)
            hess(3,IAtm,1,JAtm) = hess(3,IAtm,1,JAtm) +
     $       sva(3,IAtm)*DenXA(1,JAtm) + svb(3,IAtm)*DenXB(1,JAtm) +
     $       sga(3,IAtm)*gxaa(1,JAtm)  + sgb(3,IAtm)*gxbb(1,JAtm) +
     $       sgc(3,IAtm)*gxab(1,JAtm)  + ra*DenXXA(3,IAtm,1,JAtm) +
     $       rb*DenXXB(3,IAtm,1,JAtm)  + ga*gxxaa(3,IAtm,1,JAtm) +
     $       gb*gxxbb(3,IAtm,1,JAtm)   + gc*gxxab(3,IAtm,1,JAtm)
            hess(1,IAtm,2,JAtm) = hess(1,IAtm,2,JAtm) +
     $       sva(1,IAtm)*DenXA(2,JAtm) + svb(1,IAtm)*DenXB(2,JAtm) +
     $       sga(1,IAtm)*gxaa(2,JAtm)  + sgb(1,IAtm)*gxbb(2,JAtm) +
     $       sgc(1,IAtm)*gxab(2,JAtm)  + ra*DenXXA(1,IAtm,2,JAtm) +
     $       rb*DenXXB(1,IAtm,2,JAtm)  + ga*gxxaa(1,IAtm,2,JAtm) +
     $       gb*gxxbb(1,IAtm,2,JAtm)   + gc*gxxab(1,IAtm,2,JAtm)
            hess(2,IAtm,2,JAtm) = hess(2,IAtm,2,JAtm) +
     $       sva(2,IAtm)*DenXA(2,JAtm) + svb(2,IAtm)*DenXB(2,JAtm) +
     $       sga(2,IAtm)*gxaa(2,JAtm)  + sgb(2,IAtm)*gxbb(2,JAtm) +
     $       sgc(2,IAtm)*gxab(2,JAtm)  + ra*DenXXA(2,IAtm,2,JAtm) +
     $       rb*DenXXB(2,IAtm,2,JAtm)  + ga*gxxaa(2,IAtm,2,JAtm) +
     $       gb*gxxbb(2,IAtm,2,JAtm)   + gc*gxxab(2,IAtm,2,JAtm)
            hess(3,IAtm,2,JAtm) = hess(3,IAtm,2,JAtm) +
     $       sva(3,IAtm)*DenXA(2,JAtm) + svb(3,IAtm)*DenXB(2,JAtm) +
     $       sga(3,IAtm)*gxaa(2,JAtm)  + sgb(3,IAtm)*gxbb(2,JAtm) +
     $       sgc(3,IAtm)*gxab(2,JAtm)  + ra*DenXXA(3,IAtm,2,JAtm) +
     $       rb*DenXXB(3,IAtm,2,JAtm)  + ga*gxxaa(3,IAtm,2,JAtm) +
     $       gb*gxxbb(3,IAtm,2,JAtm)   + gc*gxxab(3,IAtm,2,JAtm)
            hess(1,IAtm,3,JAtm) = hess(1,IAtm,3,JAtm) +
     $       sva(1,IAtm)*DenXA(3,JAtm) + svb(1,IAtm)*DenXB(3,JAtm) +
     $       sga(1,IAtm)*gxaa(3,JAtm)  + sgb(1,IAtm)*gxbb(3,JAtm) +
     $       sgc(1,IAtm)*gxab(3,JAtm)  + ra*DenXXA(1,IAtm,3,JAtm) +
     $       rb*DenXXB(1,IAtm,3,JAtm)  + ga*gxxaa(1,IAtm,3,JAtm) +
     $       gb*gxxbb(1,IAtm,3,JAtm)   + gc*gxxab(1,IAtm,3,JAtm)
            hess(2,IAtm,3,JAtm) = hess(2,IAtm,3,JAtm) +
     $       sva(2,IAtm)*DenXA(3,JAtm) + svb(2,IAtm)*DenXB(3,JAtm) +
     $       sga(2,IAtm)*gxaa(3,JAtm)  + sgb(2,IAtm)*gxbb(3,JAtm) +
     $       sgc(2,IAtm)*gxab(3,JAtm)  + ra*DenXXA(2,IAtm,3,JAtm) +
     $       rb*DenXXB(2,IAtm,3,JAtm)  + ga*gxxaa(2,IAtm,3,JAtm) +
     $       gb*gxxbb(2,IAtm,3,JAtm)   + gc*gxxab(2,IAtm,3,JAtm)
            hess(3,IAtm,3,JAtm) = hess(3,IAtm,3,JAtm) +
     $       sva(3,IAtm)*DenXA(3,JAtm) + svb(3,IAtm)*DenXB(3,JAtm) +
     $       sga(3,IAtm)*gxaa(3,JAtm)  + sgb(3,IAtm)*gxbb(3,JAtm) +
     $       sgc(3,IAtm)*gxab(3,JAtm)  + ra*DenXXA(3,IAtm,3,JAtm) +
     $       rb*DenXXB(3,IAtm,3,JAtm)  + ga*gxxaa(3,IAtm,3,JAtm) +
     $       gb*gxxbb(3,IAtm,3,JAtm)   + gc*gxxab(3,IAtm,3,JAtm)
          EndDo
        EndDo
c
c  --  add weight derivatives contribution.
c      the sva array is used to store a partial sum
c
        do ia=1,natoms
          sva(1,ia)=prap*denxa(1,ia)+prbp*denxb(1,ia)+pgap*gxaa(1,ia)+
     $               pgbp*gxbb(1,ia)+pgcp*gxab(1,ia)
          sva(2,ia)=prap*denxa(2,ia)+prbp*denxb(2,ia)+pgap*gxaa(2,ia)+
     $               pgbp*gxbb(2,ia)+pgcp*gxab(2,ia)
          sva(3,ia)=prap*denxa(3,ia)+prbp*denxb(3,ia)+pgap*gxaa(3,ia)+
     $               pgbp*gxbb(3,ia)+pgcp*gxab(3,ia)
        enddo
        do ia=1,natoms
        if(ia.ne.icntr)then
          do ja=ia,natoms
          if(ja.ne.icntr)then
          hess(1,ia,1,ja)=hess(1,ia,1,ja)+gwt(1,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(1,ia)+hwt(1,ia,1,ja,ipp)
          hess(2,ia,1,ja)=hess(2,ia,1,ja)+gwt(2,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(2,ia)+hwt(2,ia,1,ja,ipp)
          hess(3,ia,1,ja)=hess(3,ia,1,ja)+gwt(3,ia,ipp)*sva(1,ja)
     $             +gwt(1,ja,ipp)*sva(3,ia)+hwt(3,ia,1,ja,ipp)
          hess(1,ia,2,ja)=hess(1,ia,2,ja)+gwt(1,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(1,ia)+hwt(1,ia,2,ja,ipp)
          hess(2,ia,2,ja)=hess(2,ia,2,ja)+gwt(2,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(2,ia)+hwt(2,ia,2,ja,ipp)
          hess(3,ia,2,ja)=hess(3,ia,2,ja)+gwt(3,ia,ipp)*sva(2,ja)
     $             +gwt(2,ja,ipp)*sva(3,ia)+hwt(3,ia,2,ja,ipp)
          hess(1,ia,3,ja)=hess(1,ia,3,ja)+gwt(1,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(1,ia)+hwt(1,ia,3,ja,ipp)
          hess(2,ia,3,ja)=hess(2,ia,3,ja)+gwt(2,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(2,ia)+hwt(2,ia,3,ja,ipp)
          hess(3,ia,3,ja)=hess(3,ia,3,ja)+gwt(3,ia,ipp)*sva(3,ja)
     $             +gwt(3,ja,ipp)*sva(3,ia)+hwt(3,ia,3,ja,ipp)
          endif
          enddo
        endif
        enddo
        call secund(t6)
        tqh=tqh+t6-t5
 200    continue
      endif
c
      return
      end
c =====================================================================
      subroutine cwmpf(dft,    npp,    nbas,   nbf,    natoms,
     $                 natonce,nb,     ne,     nbatm,  thrsh,
     $                 da,     dm,     gden,   wght,   pra,
     $                 prara,  prarb,  pga,    pgc,    praga,
     $                 pragb,  pragc,  pgaga,  pgagb,  pgagc,
     $                 pgcgc,  vao,    vaox,   vaoxx,  inb,
     $                 vm,     indx,   denx,   gdx,    gx,
     $                 sv,     sw,     icntr,  gwt,    fda,
     $                 td1,     tg1,   tsw,    tqf)
c     implicit real*8 (a-h,o-z)
      implicit none             ! trying to avoid typing errors
c
c
c  carries out numerical integration and accumulates contribution
c  from current grid point into derivative fock matrices
c  taking into account the weight derivatives.
c
c  natonce components of the derivative fock matrices are computed in
c  one pass (components nb to ne, with icntr <= nb or icntr >= ne).
c
c  ** closed shell **
c
c  arguments
c
c  dft     -  method flag (note: all methods include slater exchange)
c              1 - 3 - local correlation
c             >3 - nonlocal
c  npp     -  number of contributing (non-zero) grid points this batch
c  nbas    -  total number of basis functions
c  nbf     -  indexing array to location of "non-zero" basis functions
c  natoms  -  number of atoms
c  natonce -  number of atomic derivatives to compute in this pass
c  nb,ne   -  first and last deivative to compute (natonce=ne-nb+1)
c  nbatm   -  basis functions --> atom index
c  thrsh   -  threshold for neglect of contribution
c  da      -  closed-shell density matrix (lower triangle)
c  dm      -  maximum density matrix element per column
c  gden    -  density gradient at grid points (non local only, dft > 3)
c  wght    -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density (dft > 3)
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
c  vao     -  "non-zero" basis function values at grid points
c  potx    -  gradient functional derivative (dft > 3 only)
c  potxx   -  gradient functional second derivative (dft > 3)
c  potxd   -  density-gradient mixed second derivative (dft > 3)
c  vaox    -  basis function 1st derivatives at grid points
c  vaoxx   -  basis function 2nd derivatives at grid points
c  inb     -  indexing array for non-zero entries to vao
c  vm      -  array containing maximum magnitude ao per grid point
c  indx    -  index into contributing columns of vao
c  denx    -  scratch storage for density gradient per grid point
c  gdx     -  ditto for atomic gradient of density gradient (dft > 3)
c  gx      -  ditto for atomic gradient of gradient invariant (dft > 3)
c  sv      -  storage for coefficient of vao(i)*vao(j) in quadrature
c          -  formula (dft > 3)
c  sw      -  storage for coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  icntr   -  current atomic center
c  gwt     -  gradient of quadrature weight
c
c  on exit
c
c  fda     -  contribution to derivative fock matrices
c  hess    -  contribution to hessian matrix
c
c
      integer npp,nbas,natoms,icntr,natonce,nb,ne
      real*8 wght(*),gwt(3,natoms,*)
      real*8 vao(*),vaox(3,*),vaoxx(6,*)
      real*8 vm(*),da(nbas*(nbas+1)/2),dm(nbas),gden(3,*)
      real*8 pra(npp),pga(npp),pgc(npp),prara(npp)
      real*8 prarb(npp),praga(npp),pragb(npp),pragc(npp)
      real*8 pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
      real*8 denx(3,natoms)
      real*8 gdx(3,3,natoms)
      real*8 gx(3,natoms)
      real*8 sv(3,natoms),sw(3,3,natoms)
      integer dft,nbf(*),inb(*),indx(npp),nbatm(nbas)
      real*8 fda(3,nb:ne,*)
      real*8 thrsh
      real*8 td1,tg1,tsw,tqf
c
      real*8 zero,two,three,four,six,twelve,half
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0,half=0.5d0)
      real*8 epsi,alot
      parameter (epsi=1.0d-15,alot=1.0d17)
      logical  dodenx,dogradxx,dox,doval
      integer nat3,ip,ipp,i,j,ii,jj,ij,it
      integer ia,ja,iatm,jatm,iixx,katm,k1,isv,isw
      real*8 ra,ra2,ga,gc,rara,rara2,raga,ragb,ragc,prap
      real*8 gaga,gagb,gagc,gcgc,vmx,vmx2,wg,prg,prg2,pgg,pgg2,pg,pg2
      real*8 abra,abrara,thrx,thrx1,thrxx,valt,abdaij,abdaijm,abval
      real*8 valmm,valmm1,valmm2,valxx,valxy,valyy,valxz,valyz,valzz
      real*8 abvj,abvvj,valj,valjx,valjy,valjz,abijt,abij
      real*8 thrsh1,thtest
      real*8 t1,t2,t3,t4,t5,t6
      real*8 gwtmax,dmx,valx,valy,valz,daij,xyc,xzc,yzc
      real*8 dmaxyz,vmax,val,vali,vdxc,valm,vdjx,xc
      real*8 hvalx,hvaly,hvalz,dx,dy,dz,sxx,sxy,sxz
      real*8 vd1,vd2,vd3,svmax,swmax,tswmax
      real*8 vmax1,vmax2,vmax3
      real*8 valix,valiy,valiz,xvi,xxvx,xxvy,xxvz,valtest,valmxj
      real*8 valmx,valmx1,xvj,valm1,val1,smaxmx,val2,val3,val4
      real*8 vij,vijx,vijy,vijz,gij,hdx,hdy,hdz,hgx,hgy,hgz
      real*8 daijt,smax,dmax,potp,potxp,potp2,potxp2
      real*8 sxxi,sxyi,sxzi,dpmax2
      real*8 valt1,valt2
c
c consistency check
c
      if(icntr.ge.nb.and.icntr.le.ne)then
        call nerror(1,'cwmpf','icntr is wrong',icntr,nb)
      endif
c
      nat3 = 3*natoms
c
      if(dft.le.3) then
c
c  local only
c
      do 50 ip=1,npp
      ipp = indx(ip)
      wg  = wght(ipp)
      prap=pra(ip)
      ra = prap*wg
      rara = prara(ip)*wg
      vmx = vm(ipp)
c
c   initializations for weight derivatives
c
      gwtmax=zero
      do ia=1,natoms
        if(ia.ne.icntr)then
          gwt(1,ia,ipp)=gwt(1,ia,ipp)*prap
          gwt(2,ia,ipp)=gwt(2,ia,ipp)*prap
          gwt(3,ia,ipp)=gwt(3,ia,ipp)*prap
          gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +    abs(gwt(3,ia,ipp)))
        endif
      enddo
      call zeroit(denx,nat3)
c
c    form gradient density at current grid point.
c    only components nb through ne are formed
c
c    note:
c    for the closed shell case, the array da contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient densities is omitted here, so
c    denx will contain half the gradient and hessian
c    of the total closed shell density
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrara=abs(rara)
      dodenx=abrara.gt.epsi
      thrx=alot
      if(dodenx)thrx=thrsh/abrara
      thrx1=thrx/vmx2
      abra=abs(ra)
      if(.not.dodenx)goto 25 ! should never happen!
c
c   now compute density derivatives
c
      do 20 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        dmx = dm(ii)
        if(dmx*vmx**2.lt.thrsh) go to 20
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        if(iatm.lt.nb.or.iatm.gt.ne)goto 20
        valx = vaox(1,i)
        valy = vaox(2,i)
        valz = vaox(3,i)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx*vmx
        dox=(dodenx.and.(valt.gt.thrx1.or.valt**2.gt.thrx))
        if(.not.dox) goto 20
        do 18 j=nbf(ipp)+1,i
          jj = inb(j)
          ij = it + jj
          jatm = nbatm(jj)
          daij = da(ij)*vao(j)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
c
c -- atomic gradient of density
c
          denx(1,iatm) = denx(1,iatm) - daij*valx
          denx(2,iatm) = denx(2,iatm) - daij*valy
          denx(3,iatm) = denx(3,iatm) - daij*valz
 18     continue
        do 19 j=i+1,nbf(ipp+1)
          jj = inb(j)
          ij = (jj*(jj-1))/2 + ii
          jatm = nbatm(jj)
          daij = da(ij)*vao(j)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
c
c -- atomic gradient of density
c
          denx(1,iatm) = denx(1,iatm) - daij*valx
          denx(2,iatm) = denx(2,iatm) - daij*valy
          denx(3,iatm) = denx(3,iatm) - daij*valz
 19     continue
 20   continue
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c    numerical quadrature for the derivative fock matrix.
c
c -- get maximum absolute value of 1st-order density at this grid point
c
      call absmax(nat3,denx,iixx,dmaxyz)
c
c -- global threshold testing
c
      vmax = max(abra*vmx2,abrara*dmaxyz*vmx2)
      if(vmax.lt.thrsh) go to 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      do 40 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        vali = vao(i)
        val = vali*ra
        vdxc = vali*rara
        valx = vaox(1,i)*ra
        valy = vaox(2,i)*ra
        valz = vaox(3,i)*ra
        abval= abs(val)
        doval= abval.gt.thrsh1
        valm = max(abs(valx),abs(valy),abs(valz))
        valt = max(abval,valm,abs(vdxc)*dmaxyz)
        if(valt.gt.thrsh1) then
          do 30 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            jatm = nbatm(jj)
            valj = vao(j)
            if(abs(valj)*valm.gt.thrsh)then
              if(iatm.ge.nb.and.iatm.le.ne) then
                fda(1,iatm,ij) = fda(1,iatm,ij) - valx*valj
                fda(2,iatm,ij) = fda(2,iatm,ij) - valy*valj
                fda(3,iatm,ij) = fda(3,iatm,ij) - valz*valj
              endif
            endif
            if(doval)then
              if(jatm.ge.nb.and.jatm.le.ne) then
                fda(1,jatm,ij) = fda(1,jatm,ij) - val*vaox(1,j)
                fda(2,jatm,ij) = fda(2,jatm,ij) - val*vaox(2,j)
                fda(3,jatm,ij) = fda(3,jatm,ij) - val*vaox(3,j)
              endif
            endif
c
c -- contribution of functional derivative
c
            vdjx = vdxc*valj
            if(abs(vdjx)*dmaxyz.gt.thrsh) then
              do katm=nb,ne
                fda(1,katm,ij)=fda(1,katm,ij) + vdjx*denx(1,katm)
                fda(2,katm,ij)=fda(2,katm,ij) + vdjx*denx(2,katm)
                fda(3,katm,ij)=fda(3,katm,ij) + vdjx*denx(3,katm)
              enddo
            endif
c
c -- weight derivatives contribution
c
            xc=vali*valj
            do ia=nb,ne
              fda(1,ia,ij)=fda(1,ia,ij)+xc*gwt(1,ia,ipp)
              fda(2,ia,ij)=fda(2,ia,ij)+xc*gwt(2,ia,ipp)
              fda(3,ia,ij)=fda(3,ia,ij)+xc*gwt(3,ia,ipp)
            enddo
 30       continue
        endif
 40   continue
 45   continue
      call secund(t3)
      tqf=tqf+t3-t2
c
c   local dft ends here
c
 50   continue
cc
      else
cc
c
c   non-local dft
c
      do 200 ip=1,npp
        ipp = indx(ip)
        wg = wght(ipp)
        vmx = vm(ipp)
        ra = pra(ip)*wg
        ga = pga(ip)*wg
        gc = pgc(ip)*wg
        rara = (prara(ip)+prarb(ip))*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  some sums of the potentials that will be used later
c
        prg=raga+ragb+ragc
        pgg=gaga+gagb+two*gagc+half*gcgc
        pg=ga+half*gc
c
c  potential combinations for weight derivatives contributions
c
        potp = pra(ip)
        potxp = pga(ip)+half*pgc(ip)
        potp2 = potp+potp
        potxp2 = potxp+potxp
c
c  density gradient at current point
c
        dx=gden(1,ipp)
        dy=gden(2,ipp)
        dz=gden(3,ipp)
c
c  zero out derivatives of densities and gradients
c
        call zeroit(denx,nat3)
        call zeroit(gdx,3*nat3)
        call zeroit(gx,nat3)
c
c   initializations for weight derivatives
c
c   unlike the local case above, we do not multiply
c   the weight gradient by the potential at this stage
c
        gwtmax=zero
        do ia=1,natoms
          if(ia.ne.icntr)then
            gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +        abs(gwt(3,ia,ipp)))
          endif
        enddo
c
c    form the atomic gradient of density
c    and density gradient at current grid point
c    only components nb through ne are formed
c
c    note:
c    for the closed shell case, the array da contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx will contain half the atomic gradient
c    of the total closed shell density. likevise, gdx will
c    contain half the atomic gradient of the total closed
c    shell density gradient
c
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c
        call secund(t1)
        abra=abs(ra)
        do 220 i=nbf(ipp)+1,nbf(ipp+1)
          ii = inb(i)
          dmx = dm(ii)       ! max. element of first order density
          if(dmx.gt.epsi)then
             thtest=thrsh/(vmx*dmx)
          else
            goto 220
          endif
          it = (ii*(ii-1))/2
          iatm = nbatm(ii)
          if(iatm.lt.nb.or.iatm.gt.ne) goto 220
          valx = vaox(1,i)
          valy = vaox(2,i)
          valz = vaox(3,i)
          valm = max(abs(valx),abs(valy),abs(valz))
          valxx = vaoxx(1,i)
          valxy = vaoxx(2,i)
          valyy = vaoxx(3,i)
          valxz = vaoxx(4,i)
          valyz = vaoxx(5,i)
          valzz = vaoxx(6,i)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          do 218 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            valj=vao(j)
            abvj=abs(valj)
            valjx=vaox(1,j)
            valjy=vaox(2,j)
            valjz=vaox(3,j)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            daijt=da(ij)
            abijt=abs(daijt)
            daij = daijt*valj
            abij=abs(daij)
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              denx(1,iatm) = denx(1,iatm) - daij*valx
              denx(2,iatm) = denx(2,iatm) - daij*valy
              denx(3,iatm) = denx(3,iatm) - daij*valz
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              gdx(1,1,iatm)=gdx(1,1,iatm)-daijt*(valj*valxx+valx*valjx)
              gdx(1,2,iatm)=gdx(1,2,iatm)-daijt*(valj*valxy+valy*valjx)
              gdx(1,3,iatm)=gdx(1,3,iatm)-daijt*(valj*valxz+valz*valjx)
              gdx(2,1,iatm)=gdx(2,1,iatm)-daijt*(valj*valxy+valx*valjy)
              gdx(2,2,iatm)=gdx(2,2,iatm)-daijt*(valj*valyy+valy*valjy)
              gdx(2,3,iatm)=gdx(2,3,iatm)-daijt*(valj*valyz+valz*valjy)
              gdx(3,1,iatm)=gdx(3,1,iatm)-daijt*(valj*valxz+valx*valjz)
              gdx(3,2,iatm)=gdx(3,2,iatm)-daijt*(valj*valyz+valy*valjz)
              gdx(3,3,iatm)=gdx(3,3,iatm)-daijt*(valj*valzz+valz*valjz)
            endif
 218      continue
          do 219 j=i+1,nbf(ipp+1)
            jj = inb(j)
            ij = (jj*(jj-1))/2 + ii
            valj=vao(j)
            abvj=abs(valj)
            valjx=vaox(1,j)
            valjy=vaox(2,j)
            valjz=vaox(3,j)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            daijt= da(ij)
            abijt=abs(daijt)
            daij = daijt*valj
            abij=abs(daij)
c
c -- atomic gradient of density
c
              if(abij*valm.gt.thrsh)then
                denx(1,iatm) = denx(1,iatm) - daij*valx
                denx(2,iatm) = denx(2,iatm) - daij*valy
                denx(3,iatm) = denx(3,iatm) - daij*valz
              endif
c
c -- atomic gradient of density gradient
c
             if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              gdx(1,1,iatm)=gdx(1,1,iatm)-daijt*(valj*valxx+valx*valjx)
              gdx(1,2,iatm)=gdx(1,2,iatm)-daijt*(valj*valxy+valy*valjx)
              gdx(1,3,iatm)=gdx(1,3,iatm)-daijt*(valj*valxz+valz*valjx)
              gdx(2,1,iatm)=gdx(2,1,iatm)-daijt*(valj*valxy+valx*valjy)
              gdx(2,2,iatm)=gdx(2,2,iatm)-daijt*(valj*valyy+valy*valjy)
              gdx(2,3,iatm)=gdx(2,3,iatm)-daijt*(valj*valyz+valz*valjy)
              gdx(3,1,iatm)=gdx(3,1,iatm)-daijt*(valj*valxz+valx*valjz)
              gdx(3,2,iatm)=gdx(3,2,iatm)-daijt*(valj*valyz+valy*valjz)
              gdx(3,3,iatm)=gdx(3,3,iatm)-daijt*(valj*valzz+valz*valjz)
             endif
 219      continue
 220    continue
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient of density gradient invariant
c
        do iatm=1,natoms
          gx(1,iatm)=two*(dx*gdx(1,1,iatm)+
     $                   dy*gdx(2,1,iatm)+dz*gdx(3,1,iatm))
          gx(2,iatm)=two*(dx*gdx(1,2,iatm)+
     $                   dy*gdx(2,2,iatm)+dz*gdx(3,2,iatm))
          gx(3,iatm)=two*(dx*gdx(1,3,iatm)+
     $                   dy*gdx(2,3,iatm)+dz*gdx(3,3,iatm))
        enddo
c
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the fock matrix
c    (i.e., quantities v, w and x at page 7436 of johnson and frisch).
c
        prg2=two*prg
        pgg2=two*pgg
        pg2=two*pg
        sxx=pg2*dx
        sxy=pg2*dy
        sxz=pg2*dz
        do iatm=1,natoms
          vd1=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          vd2=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          vd3=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          sv(1,iatm)=rara*denx(1,iatm)+prg*gx(1,iatm)
          sv(2,iatm)=rara*denx(2,iatm)+prg*gx(2,iatm)
          sv(3,iatm)=rara*denx(3,iatm)+prg*gx(3,iatm)
          sw(1,1,iatm)=pg2*gdx(1,1,iatm)+dx*vd1
          sw(2,1,iatm)=pg2*gdx(2,1,iatm)+dy*vd1
          sw(3,1,iatm)=pg2*gdx(3,1,iatm)+dz*vd1
          sw(1,2,iatm)=pg2*gdx(1,2,iatm)+dx*vd2
          sw(2,2,iatm)=pg2*gdx(2,2,iatm)+dy*vd2
          sw(3,2,iatm)=pg2*gdx(3,2,iatm)+dz*vd2
          sw(1,3,iatm)=pg2*gdx(1,3,iatm)+dx*vd3
          sw(2,3,iatm)=pg2*gdx(2,3,iatm)+dy*vd3
          sw(3,3,iatm)=pg2*gdx(3,3,iatm)+dz*vd3
        enddo
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    fock matrix.
c
c -- get the maximum absolute value of coefficients v and w
c
        call absmax(nat3,sv,isv,svmax)
        call absmax(3*nat3,sw,isw,swmax)
        tswmax=three*swmax
c -- global threshold testing
        vmx2=vmx*vmx
        smax = abs(sxx+sxy+sxz)
        dmax = abs(dx+dy+dz)
        vmax3 = (svmax+tswmax+tswmax)*vmx2
        vmax1 = (abra+two*smax)*vmx2
        dpmax2 = potxp2*dmax*(vmx2+vmx2)
        vmax2 = gwtmax*(potp*vmx2+dpmax2)
        vmax=max(vmax1,vmax2,vmax3)
        if(vmax.lt.thrsh)  go to 245
c
c -- numerical quadrature  for derivative fock matrix
c
        do 240 i=nbf(ipp)+1,nbf(ipp+1)
          vali = vao(i)
          valix = vaox(1,i)
          valiy = vaox(2,i)
          valiz = vaox(3,i)
          xvi=sxx*valix+sxy*valiy+sxz*valiz+ra*vali
          sxxi=sxx*vali
          sxyi=sxy*vali
          sxzi=sxz*vali
          xxvx=sxx*vaoxx(1,i)+sxy*vaoxx(2,i)+sxz*vaoxx(4,i)+ra*valix
          xxvy=sxx*vaoxx(2,i)+sxy*vaoxx(3,i)+sxz*vaoxx(5,i)+ra*valiy
          xxvz=sxx*vaoxx(4,i)+sxy*vaoxx(5,i)+sxz*vaoxx(6,i)+ra*valiz
          valm=abs(vali)
          valm1=abs(xvi)
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1=max(abs(xxvx),abs(xxvy),abs(xxvz))
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          val4=gwtmax*(potp*valm*vmx+dpmax2)
          valt1=max(val1,val3)
          valt2=max(val3,val4)
          valtest=max(valt1,valt2)
          if(valtest.gt.thrsh)then
            ii = inb(i)
            it = (ii*(ii-1))/2
            iatm = nbatm(ii)
            do 230 j=nbf(ipp)+1,i
              jj = inb(j)
              ij = it + jj
              jatm = nbatm(jj)
              valj = vao(j)
              valjx = vaox(1,j)
              valjy = vaox(2,j)
              valjz = vaox(3,j)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              xvj=sxx*valjx+sxy*valjy+sxz*valjz
              if(valmx1*abs(valj)+valmx*abs(xvj).gt.thrsh)then
              if(iatm.ge.nb.and.iatm.le.ne) then
                fda(1,iatm,ij) = fda(1,iatm,ij) - xxvx*valj - valix*xvj
                fda(2,iatm,ij) = fda(2,iatm,ij) - xxvy*valj - valiy*xvj
                fda(3,iatm,ij) = fda(3,iatm,ij) - xxvz*valj - valiz*xvj
              endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
              if(jatm.ge.nb.and.jatm.le.ne) then
                fda(1,jatm,ij) = fda(1,jatm,ij) - xvi*valjx -
     $          vaoxx(1,j)*sxxi - vaoxx(2,j)*sxyi - vaoxx(4,j)*sxzi
                fda(2,jatm,ij) = fda(2,jatm,ij) - xvi*valjy -
     $          vaoxx(2,j)*sxxi - vaoxx(3,j)*sxyi - vaoxx(5,j)*sxzi
                fda(3,jatm,ij) = fda(3,jatm,ij) - xvi*valjz -
     $          vaoxx(4,j)*sxxi - vaoxx(5,j)*sxyi - vaoxx(6,j)*sxzi
              endif
              endif
              vij=vali*valj
              vijx=valix*valj+valjx*vali
              vijy=valiy*valj+valjy*vali
              vijz=valiz*valj+valjz*vali
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).gt.thrsh)then
                do katm=nb,ne
                  fda(1,katm,ij) = fda(1,katm,ij) + (sv(1,katm)*vij+
     $            sw(1,1,katm)*vijx+sw(2,1,katm)*vijy+sw(3,1,katm)*vijz)
                  fda(2,katm,ij) = fda(2,katm,ij) + (sv(2,katm)*vij+
     $            sw(1,2,katm)*vijx+sw(2,2,katm)*vijy+sw(3,2,katm)*vijz)
                  fda(3,katm,ij) = fda(3,katm,ij) + (sv(3,katm)*vij+
     $            sw(1,3,katm)*vijx+sw(2,3,katm)*vijy+sw(3,3,katm)*vijz)
                enddo
              endif
c
c -- weight derivatives contribution
c
              gij=potp*vij+potxp2*(dx*vijx+dy*vijy+dz*vijz)
              if(gwtmax*abs(gij).gt.thrsh)then
                do ia=nb,ne
                  fda(1,ia,ij)=fda(1,ia,ij)+gwt(1,ia,ipp)*gij
                  fda(2,ia,ij)=fda(2,ia,ij)+gwt(2,ia,ipp)*gij
                  fda(3,ia,ij)=fda(3,ia,ij)+gwt(3,ia,ipp)*gij
                enddo
              endif
 230        continue
          endif
 240    continue
 245    continue
        call secund(t5)
        tqf=tqf+t5-t4
 200    continue
      endif
c
      return
      end
c =====================================================================
      subroutine cwmpfu(dft,    npp,    nbas,   nbf,    natoms,
     $                   natonce,nb,    ne,     nbatm,  thrsh,
     $                   da,     db,     dm,    gdena,  gdenb,
     $                   wght,   pra,    prb,   prara,  prbrb,
     $                   prarb,  pga,    pgb,   pgc,    praga,
     $                   pragb,  pragc,  prbga, prbgb,  prbgc,
     $                   pgaga,  pgagb,  pgagc, pgbgb,  pgbgc,
     $                   pgcgc,  vao,    vaox,  vaoxx,  inb,
     $                   vm,     indx,   denxa, denxb,  gdxa,
     $                   gdxb,   gxaa,   gxbb,  gxab,   sva,
     $                   svb,    swa,    swb,   sga,    sgb,
     $                   sgc,    icntr,  gwt,   fda,    fdb,
     $                   td1,    tg1,    tsw,   tqf)
c     implicit real*8 (a-h,o-z)
      implicit none             ! trying to avoid typing errors
c
c
c  carries out numerical integration and accumulates contribution
c  from current grid point into derivative fock matrices
c  taking into account the weight derivatives.
c
c  natonce components of the derivative fock matrices are computed in
c  one pass (components nb to ne, with icntr <= nb or icntr >= ne).
c
c  ** open shell **
c
c  arguments
c
c  dft     -  method flag (note: all methods include slater exchange)
c              1 - 3 - local correlation
c             >3 - nonlocal
c  npp     -  number of contributing (non-zero) grid points this batch
c  nbas    -  total number of basis functions
c  nbf     -  indexing array to location of "non-zero" basis functions
c  natoms  -  number of atoms
c  natonce -  number of atomic derivatives to compute in this pass
c  nb,ne   -  first and last deivative to compute (natonce=ne-nb+1)
c  nbatm   -  basis functions --> atom index
c  thrsh   -  threshold for neglect of contribution
c  da      -  alpha density matrix (lower triangle)
c  db      -  beta density matrix (lower triangle)
c  dm      -  maximum density matrix element per column
c  gdena   -  alpha density gradient at grid points (only dft > 3)
c  gdena   -  beta density gradient at grid points (dft > 3)
c  wght    -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prb   -  Functional derivative w.r.t. beta density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgb   -  Funct. deriv. w.r.t. beta gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  prbga -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha grad.
C           (dft > 3)
C  prbgb -  Funct. 2nd. deriv. w.r.t. beta dens. and  grad. (dft > 3)
C  prbgc -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgbgb -  Funct. 2nd. deriv. w.r.t. beta grad. (dft > 3)
C  pgbgc -  Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
c  vao     -  "non-zero" basis function values at grid points
c  vaox    -  basis function 1st derivatives at grid points
c  vaoxx   -  basis function 2nd derivatives at grid points
c  inb     -  indexing array for non-zero entries to vao
c  vm      -  array containing maximum magnitude ao per grid point
c  indx    -  index into contributing columns of vao
c  denxa   -  scratch storage for alpah density gradient per grid point
c  denxb   -  scratch storage for beta density gradient per grid point
c  gdxa    -  ditto for alpha atomic gradient of density gradient
c             (dft > 3)
c  gdxb    -  ditto for beta atomic gradient of density gradient
c             (dft > 3)
c  gxaa    -  ditto for atomic gradient of alpha gradient invariant
c             (dft > 3)
c  gxbb    -  ditto for atomic gradient of beta gradient invariant
c             (dft > 3)
c  gxbb    -  ditto for atomic gradient of alpha beta gradient
c             invariant (dft > 3)
c  sva     -  storage for alpha coefficient of vao(i)*vao(j) in
c             quadrature formula
c  svb     -  storage for beta coefficient of vao(i)*vao(j) in
c             quadrature formula
c  swa     -  storage for alpha coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  swb     -  storage for beta coefficient of
c               (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
c  sga     -  storage for partial sum for hessian quadrature (dft > 3)
c  sgb     -  storage for partial sum for hessian quadrature (dft > 3)
c  sgc     -  storage for partial sum for hessian quadrature (dft > 3)
c  icntr   -  current atomic center
c  gwt     -  gradient of quadrature weight
c
c  on exit
c
c  fda     -  contribution to alpha derivative fock matrices
c  fdb     -  contribution to beta derivative fock matrices
c
c
      integer npp,nbas,natoms,icntr,natonce,nb,ne
      real*8 wght(*),gwt(3,natoms,*)
      real*8 vao(*),vaox(3,*),vaoxx(6,*)
      real*8 vm(*),da(nbas*(nbas+1)/2),dm(nbas),gdena(3,*)
      real*8 db(nbas*(nbas+1)/2),gdenb(3,*)
      real*8 pra(npp),pga(npp),pgc(npp),prara(npp)
      real*8 prb(npp),pgb(npp),prbrb(npp)
      real*8 prarb(npp),praga(npp),pragb(npp),pragc(npp)
      real*8 prbga(npp),prbgb(npp),prbgc(npp)
      real*8 pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
      real*8 pgbgb(npp),pgbgc(npp)
      real*8 denxa(3,natoms)
      real*8 denxb(3,natoms)
      real*8 gdxa(3,3,natoms)
      real*8 gdxb(3,3,natoms)
      real*8 gxaa(3,natoms)
      real*8 gxbb(3,natoms)
      real*8 gxab(3,natoms)
      real*8 sva(3,natoms),svb(3,natoms)
      real*8 swa(3,3,natoms),swb(3,3,natoms)
      real*8 sga(3,natoms),sgb(3,natoms),sgc(3,natoms)
      integer dft,nbf(*),inb(*),indx(npp),nbatm(nbas)
      real*8 fda(3,nb:ne,*)
      real*8 fdb(3,nb:ne,*)
      real*8 thrsh
      real*8 td1,tg1,tsw,tqf,tqh
c
      real*8 zero,two,three,four,six,twelve,half
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0,half=0.5d0)
      real*8 epsi,alot
      parameter (epsi=1.0d-15,alot=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,doval
      integer nat3,ip,ipp,i,j,ii,jj,ij,it
      integer ia,ja,iatm,jatm,katm,iixx,k1
      real*8 t1,t2,t3,t4,t5,t6
      real*8 wg,vmx,vmx2,prap,prbp,gwtmax,pgap,pgbp,pgcp
      real*8 ra,rb,rara,rbrb,rarb,abrr,abr
      real*8 ga,gb,gc,raga,ragb,ragc,rbga,rbgb,rbgc
      real*8 gaga,gagb,gagc,gbgb,gbgc,gcgc
      real*8 dax,day,daz,dbx,dby,dbz
      real*8 dax2,day2,daz2,dbx2,dby2,dbz2
      real*8 ga2,gb2,pgap2,pgbp2
      real*8 thrx,thrxx,thrx1,thrsh1
      real*8 dmx2,daij,dbij,daij2,dbij2,abdij,abdijm
      real*8 valx,valy,valz,valm,valt,abval
      real*8 vali,vala,valxa,valya,valza
      real*8 valj,valb,valxb,valyb,valzb
      real*8 valij,valija,valijb
      real*8 xyc,xzc,yzc
      real*8 svamax,svbmax,svmax,vmax,abgmax
      real*8 swamax,swbmax,swmax,tswmax,vmax1,vmax3
      real*8 valxx,valxy,valxz,valyy,valyz,valzz
      real*8 valmm,valm1,valmm1,valmm2
      real*8 valmx,valmx1a,valmx1b,valmx1,val1,val2,val3,valtest
      real*8 thtest,valmxj,xvja,xvjb,abxvj
      real*8 valjx,valjy,valjz,abvj,abvvj,abijt,abij
      real*8 valix,valiy,valiz,xvia,xvib
      real*8 sxxa,sxya,sxza,sxxb,sxyb,sxzb
      real*8 sxxia,sxyia,sxzia,sxxib,sxyib,sxzib,smax,smaxmx
      real*8 xxvxa,xxvya,xxvza,xxvxb,xxvyb,xxvzb
      real*8 vij,vijx,vijy,vijz
      real*8 dxija,dxijb,gija,gijb
c
c consistency check
c
      if(icntr.ge.nb.and.icntr.le.ne)then
        call nerror(1,'cwmpfu','icntr is wrong',icntr,nb)
      endif
c
      nat3 = 3*natoms
c
      if(dft.le.3) then
c
c  local only
c
      do 50 ip=1,npp
      ipp = indx(ip)
      wg  = wght(ipp)
      prap=pra(ip)
      prbp=prb(ip)
      ra = prap*wg
      rb = prbp*wg
      rara = prara(ip)*wg
      rbrb = prbrb(ip)*wg
      rarb = prarb(ip)*wg
      vmx = vm(ipp)
c
c   initializations for weight derivatives
c
      gwtmax=zero
      do ia=1,natoms
        if(ia.ne.icntr)then
          gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +    abs(gwt(3,ia,ipp)))
        endif
      enddo
      call zeroit(denxa,nat3)
      call zeroit(denxb,nat3)
c
c    form gradient density at current grid point.
c    only components nb through ne are formed
c
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrr=max(abs(rara),abs(rbrb))+abs(rarb)
      dodenx=abrr.gt.epsi
      thrx=alot
      if(dodenx)thrx=thrsh/abrr
      thrx1=thrx/vmx2
      abr=max(abs(ra),abs(rb))
      if(.not.dodenx)goto 25 ! should never happen!
c
c   now compute density derivatives
c
      do 20 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        dmx2 = dm(ii)
        if(dmx2*vmx2.lt.thrsh) go to 20
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        if(iatm.lt.nb.or.iatm.gt.ne)goto 20
        valx = vaox(1,i)
        valy = vaox(2,i)
        valz = vaox(3,i)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx2*vmx
        dox=(dodenx.and.(two*valt.gt.thrx1.or.valt**2.gt.thrx))
        if(.not.dox) goto 20
        do 18 j=nbf(ipp)+1,i
          jj = inb(j)
          ij = it + jj
          daij2=two*da(ij)
          dbij2=two*db(ij)
          daij = daij2*vao(j)
          dbij = dbij2*vao(j)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
c
c -- atomic gradient of density
c
          denxa(1,iatm) = denxa(1,iatm) - daij*valx
          denxa(2,iatm) = denxa(2,iatm) - daij*valy
          denxa(3,iatm) = denxa(3,iatm) - daij*valz
c
          denxb(1,iatm) = denxb(1,iatm) - dbij*valx
          denxb(2,iatm) = denxb(2,iatm) - dbij*valy
          denxb(3,iatm) = denxb(3,iatm) - dbij*valz
 18     continue
        do 19 j=i+1,nbf(ipp+1)
          jj = inb(j)
          ij = (jj*(jj-1))/2 + ii
          daij2 = two*da(ij)
          dbij2 = two*db(ij)
          daij = daij2*vao(j)
          dbij = dbij2*vao(j)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
c
c -- atomic gradient of density
c
          denxa(1,iatm) = denxa(1,iatm) - daij*valx
          denxa(2,iatm) = denxa(2,iatm) - daij*valy
          denxa(3,iatm) = denxa(3,iatm) - daij*valz
c
          denxb(1,iatm) = denxb(1,iatm) - dbij*valx
          denxb(2,iatm) = denxb(2,iatm) - dbij*valy
          denxb(3,iatm) = denxb(3,iatm) - dbij*valz
 19     continue
 20   continue
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c -- now form the coefficients that multiply the basis functions
c    in the expression for the Fock Matrix
c    (i.e., quantities V at page 7436 of Johnson and Frisch).
c
        DO IAtm=1,NAtoms
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)
          svb(1,iatm)=rbrb*denxb(1,iatm)+rarb*denxa(1,iatm)
          svb(2,iatm)=rbrb*denxb(2,iatm)+rarb*denxa(2,iatm)
          svb(3,iatm)=rbrb*denxb(3,iatm)+rarb*denxa(3,iatm)
        EndDO
c
c    numerical quadrature for the derivative fock matrix.
c
c -- get maximum absolute value of sva and svb
c
      Call absmax(NAt3,sva,iixx,svamax)
      Call absmax(NAt3,svb,iixx,svbmax)
      svmax=max(svamax,svbmax)
c
c -- global threshold testing
c
      VMax = Max(Abr*VMx2,Abrr*svmax*VMx2)
      if(vmax.lt.thrsh) go to 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      do 40 i=nbf(ipp)+1,nbf(ipp+1)
        ii = inb(i)
        it = (ii*(ii-1))/2
        iatm = nbatm(ii)
        vali = vao(i)
        vala = vali*ra
        valb = vali*rb
        valxa = vaox(1,i)*ra
        valya = vaox(2,i)*ra
        valza = vaox(3,i)*ra
        valxb = vaox(1,i)*rb
        valyb = vaox(2,i)*rb
        valzb = vaox(3,i)*rb
        abval= max(abs(vala),abs(valb))
        doval= abval.gt.thrsh1
        valm = max(abs(vaox(1,i)),abs(vaox(2,i)),abs(vaox(3,i)))*abr
        Valt = MAX(abval,valm,Abs(vali)*svmax)
        if(valt.gt.thrsh1) then
          do 30 j=nbf(ipp)+1,i
            jj = inb(j)
            ij = it + jj
            jatm = nbatm(jj)
            valj = vao(j)
            if(abs(valj)*valm.gt.thrsh)then
              if(iatm.ge.nb.and.iatm.le.ne) then
                fda(1,iatm,ij) = fda(1,iatm,ij) - valxa*valj
                fda(2,iatm,ij) = fda(2,iatm,ij) - valya*valj
                fda(3,iatm,ij) = fda(3,iatm,ij) - valza*valj
                fdb(1,iatm,ij) = fdb(1,iatm,ij) - valxb*valj
                fdb(2,iatm,ij) = fdb(2,iatm,ij) - valyb*valj
                fdb(3,iatm,ij) = fdb(3,iatm,ij) - valzb*valj
              endif
            endif
            if(doval)then
              if(jatm.ge.nb.and.jatm.le.ne) then
                fda(1,jatm,ij) = fda(1,jatm,ij) - vala*vaox(1,j)
                fda(2,jatm,ij) = fda(2,jatm,ij) - vala*vaox(2,j)
                fda(3,jatm,ij) = fda(3,jatm,ij) - vala*vaox(3,j)
                fdb(1,jatm,ij) = fdb(1,jatm,ij) - valb*vaox(1,j)
                fdb(2,jatm,ij) = fdb(2,jatm,ij) - valb*vaox(2,j)
                fdb(3,jatm,ij) = fdb(3,jatm,ij) - valb*vaox(3,j)
              endif
            endif
c
c -- contribution of functional derivative
c
            valij = vali*valj
            if(abs(valij)*svmax.gt.thrsh) then
              do katm=nb,ne
                fda(1,katm,ij)=fda(1,katm,ij) + valij*sva(1,katm)
                fda(2,katm,ij)=fda(2,katm,ij) + valij*sva(2,katm)
                fda(3,katm,ij)=fda(3,katm,ij) + valij*sva(3,katm)
                fdb(1,katm,ij)=fdb(1,katm,ij) + valij*svb(1,katm)
                fdb(2,katm,ij)=fdb(2,katm,ij) + valij*svb(2,katm)
                fdb(3,katm,ij)=fdb(3,katm,ij) + valij*svb(3,katm)
              enddo
            endif
c
c -- weight derivatives contribution
c
            valija=valij*prap
            valijb=valij*prbp
            do ia=nb,ne
              fda(1,ia,ij)=fda(1,ia,ij)+valija*gwt(1,ia,ipp)
              fda(2,ia,ij)=fda(2,ia,ij)+valija*gwt(2,ia,ipp)
              fda(3,ia,ij)=fda(3,ia,ij)+valija*gwt(3,ia,ipp)
              fdb(1,ia,ij)=fdb(1,ia,ij)+valijb*gwt(1,ia,ipp)
              fdb(2,ia,ij)=fdb(2,ia,ij)+valijb*gwt(2,ia,ipp)
              fdb(3,ia,ij)=fdb(3,ia,ij)+valijb*gwt(3,ia,ipp)
            enddo
 30       continue
        endif
 40   continue
 45   continue
      call secund(t3)
      tqf=tqf+t3-t2
c
c   local dft ends here
c
 50   continue
cc
      else
cc
c
c   non-local dft
c
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        prap=pra(ip)
        prbp=prb(ip)
        ra = prap*wg
        rb = prbp*wg
        pgap = pga(ip)
        pgbp = pgb(ip)
        pgcp = pgc(ip)
        pgap2=pgap+pgap
        pgbp2=pgbp+pgbp
        ga = pgap*wg
        gb = pgbp*wg
        gc = pgcp*wg
        rara = prara(ip)*wg
        rbrb = prbrb(ip)*wg
        rarb = prarb(ip)*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        rbga = prbga(ip)*wg
        rbgb = prbgb(ip)*wg
        rbgc = prbgc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gbgb = pgbgb(ip)*wg
        gbgc = pgbgc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  density gradient at current point
c
        DAX=GDENA(1,IPP)
        DAY=GDENA(2,IPP)
        DAZ=GDENA(3,IPP)
        DBX=GDENB(1,IPP)
        DBY=GDENB(2,IPP)
        DBZ=GDENB(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenXA,NAt3)
        CALL ZeroIT(GDXA,3*NAt3)
        CALL ZeroIT(DenXB,NAt3)
        CALL ZeroIT(GDXB,3*NAt3)
c   initializations for weight derivatives
c
        gwtmax=zero
        do ia=1,natoms
          if(ia.ne.icntr)then
            gwtmax=max(gwtmax,abs(gwt(1,ia,ipp)),abs(gwt(2,ia,ipp)),
     +        abs(gwt(3,ia,ipp)))
          endif
        enddo
c
c    form the atomic gradient of density
c    and density gradient at current grid point
c    only components nb through ne are formed
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c
        abgmax=max(abs(ga),abs(gb),abs(gc))
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          dmx2 = DM(II)*two   ! max. element of first order density
          if(dmx2.gt.epsi)then
             thtest=thrsh/(vmx*dmx2)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          if(iatm.lt.nb.or.iatm.gt.ne)goto 220
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2=two*DA(IJ)
            DBIJ2=two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2= two*DA(IJ)
            DBIJ2= two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and Hessian of density gradient invariant
c
c
        DO IAtm=1,NAtoms
c   alpha alpha
          GXAA(1,IAtm)=two*(DAX*GDXA(1,1,IAtm)+
     $                   DAY*GDXA(2,1,IAtm)+DAZ*GDXA(3,1,IAtm))
          GXAA(2,IAtm)=two*(DAX*GDXA(1,2,IAtm)+
     $                   DAY*GDXA(2,2,IAtm)+DAZ*GDXA(3,2,IAtm))
          GXAA(3,IAtm)=two*(DAX*GDXA(1,3,IAtm)+
     $                   DAY*GDXA(2,3,IAtm)+DAZ*GDXA(3,3,IAtm))
c   beta beta
          GXBB(1,IAtm)=two*(DBX*GDXB(1,1,IAtm)+
     $                   DBY*GDXB(2,1,IAtm)+DBZ*GDXB(3,1,IAtm))
          GXBB(2,IAtm)=two*(DBX*GDXB(1,2,IAtm)+
     $                   DBY*GDXB(2,2,IAtm)+DBZ*GDXB(3,2,IAtm))
          GXBB(3,IAtm)=two*(DBX*GDXB(1,3,IAtm)+
     $                   DBY*GDXB(2,3,IAtm)+DBZ*GDXB(3,3,IAtm))
c   alpha beta
          GXAB(1,IAtm)=DAX*GDXB(1,1,IAtm)+DBX*GDXA(1,1,IAtm)+
     $                 DAY*GDXB(2,1,IAtm)+DBY*GDXA(2,1,IAtm)+
     $                 DAZ*GDXB(3,1,IAtm)+DBZ*GDXA(3,1,IAtm)
          GXAB(2,IAtm)=DAX*GDXB(1,2,IAtm)+DBX*GDXA(1,2,IAtm)+
     $                 DAY*GDXB(2,2,IAtm)+DBY*GDXA(2,2,IAtm)+
     $                 DAZ*GDXB(3,2,IAtm)+DBZ*GDXA(3,2,IAtm)
          GXAB(3,IAtm)=DAX*GDXB(1,3,IAtm)+DBX*GDXA(1,3,IAtm)+
     $                 DAY*GDXB(2,3,IAtm)+DBY*GDXA(2,3,IAtm)+
     $                 DAZ*GDXB(3,3,IAtm)+DBZ*GDXA(3,3,IAtm)
        EndDO
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
C    and also compute partial sums used in Hessian quadrature
c
        ga2= ga+ga
        gb2= gb+gb
        DAX2=DAX+DAX
        DAY2=DAY+DAY
        DAZ2=DAZ+DAZ
        DBX2=DBX+DBX
        DBY2=DBY+DBY
        DBZ2=DBZ+DBZ
c
c     coeff. X
c
        SXXA=ga2*DAX+gc*DBX
        SXYA=ga2*DAY+gc*DBY
        SXZA=ga2*DAZ+gc*DBZ
        SXXB=gb2*DBX+gc*DAX
        SXYB=gb2*DBY+gc*DAY
        SXZB=gb2*DBZ+gc*DAZ
c
        DO IAtm=1,NAtoms
c
c     coeff. V
c
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)+
     $        raga*gxaa(1,iatm)+ragb*gxbb(1,iatm)+ragc*gxab(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)+
     $        raga*gxaa(2,iatm)+ragb*gxbb(2,iatm)+ragc*gxab(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)+
     $        raga*gxaa(3,iatm)+ragb*gxbb(3,iatm)+ragc*gxab(3,iatm)
          svb(1,iatm)=rarb*denxa(1,iatm)+rbrb*denxb(1,iatm)+
     $        rbga*gxaa(1,iatm)+rbgb*gxbb(1,iatm)+rbgc*gxab(1,iatm)
          svb(2,iatm)=rarb*denxa(2,iatm)+rbrb*denxb(2,iatm)+
     $        rbga*gxaa(2,iatm)+rbgb*gxbb(2,iatm)+rbgc*gxab(2,iatm)
          svb(3,iatm)=rarb*denxa(3,iatm)+rbrb*denxb(3,iatm)+
     $        rbga*gxaa(3,iatm)+rbgb*gxbb(3,iatm)+rbgc*gxab(3,iatm)
c
c  partial sums for Hessian
c
          sga(1,iatm)=raga*denxa(1,iatm)+rbga*denxb(1,iatm)+
     $        gaga*gxaa(1,iatm)+gagb*gxbb(1,iatm)+gagc*gxab(1,iatm)
          sga(2,iatm)=raga*denxa(2,iatm)+rbga*denxb(2,iatm)+
     $        gaga*gxaa(2,iatm)+gagb*gxbb(2,iatm)+gagc*gxab(2,iatm)
          sga(3,iatm)=raga*denxa(3,iatm)+rbga*denxb(3,iatm)+
     $        gaga*gxaa(3,iatm)+gagb*gxbb(3,iatm)+gagc*gxab(3,iatm)
c
          sgb(1,iatm)=ragb*denxa(1,iatm)+rbgb*denxb(1,iatm)+
     $        gagb*gxaa(1,iatm)+gbgb*gxbb(1,iatm)+gbgc*gxab(1,iatm)
          sgb(2,iatm)=ragb*denxa(2,iatm)+rbgb*denxb(2,iatm)+
     $        gagb*gxaa(2,iatm)+gbgb*gxbb(2,iatm)+gbgc*gxab(2,iatm)
          sgb(3,iatm)=ragb*denxa(3,iatm)+rbgb*denxb(3,iatm)+
     $        gagb*gxaa(3,iatm)+gbgb*gxbb(3,iatm)+gbgc*gxab(3,iatm)
c
          sgc(1,iatm)=ragc*denxa(1,iatm)+rbgc*denxb(1,iatm)+
     $        gagc*gxaa(1,iatm)+gbgc*gxbb(1,iatm)+gcgc*gxab(1,iatm)
          sgc(2,iatm)=ragc*denxa(2,iatm)+rbgc*denxb(2,iatm)+
     $        gagc*gxaa(2,iatm)+gbgc*gxbb(2,iatm)+gcgc*gxab(2,iatm)
          sgc(3,iatm)=ragc*denxa(3,iatm)+rbgc*denxb(3,iatm)+
     $        gagc*gxaa(3,iatm)+gbgc*gxbb(3,iatm)+gcgc*gxab(3,iatm)
c
c     coeff. W
c
          swa(1,1,iatm)=ga2*gdxa(1,1,iatm)+gc*gdxb(1,1,Iatm)+
     $        DAX2*sga(1,iatm)+DBX*sgc(1,iatm)
          swa(1,2,iatm)=ga2*gdxa(1,2,iatm)+gc*gdxb(1,2,Iatm)+
     $        DAX2*sga(2,iatm)+DBX*sgc(2,iatm)
          swa(1,3,iatm)=ga2*gdxa(1,3,iatm)+gc*gdxb(1,3,Iatm)+
     $        DAX2*sga(3,iatm)+DBX*sgc(3,iatm)
          swa(2,1,iatm)=ga2*gdxa(2,1,iatm)+gc*gdxb(2,1,Iatm)+
     $        DAY2*sga(1,iatm)+DBY*sgc(1,iatm)
          swa(2,2,iatm)=ga2*gdxa(2,2,iatm)+gc*gdxb(2,2,Iatm)+
     $        DAY2*sga(2,iatm)+DBY*sgc(2,iatm)
          swa(2,3,iatm)=ga2*gdxa(2,3,iatm)+gc*gdxb(2,3,Iatm)+
     $        DAY2*sga(3,iatm)+DBY*sgc(3,iatm)
          swa(3,1,iatm)=ga2*gdxa(3,1,iatm)+gc*gdxb(3,1,Iatm)+
     $        DAZ2*sga(1,iatm)+DBZ*sgc(1,iatm)
          swa(3,2,iatm)=ga2*gdxa(3,2,iatm)+gc*gdxb(3,2,Iatm)+
     $        DAZ2*sga(2,iatm)+DBZ*sgc(2,iatm)
          swa(3,3,iatm)=ga2*gdxa(3,3,iatm)+gc*gdxb(3,3,Iatm)+
     $        DAZ2*sga(3,iatm)+DBZ*sgc(3,iatm)
c
          swb(1,1,iatm)=gb2*gdxb(1,1,iatm)+gc*gdxa(1,1,Iatm)+
     $        DBX2*sgb(1,iatm)+DAX*sgc(1,iatm)
          swb(1,2,iatm)=gb2*gdxb(1,2,iatm)+gc*gdxa(1,2,Iatm)+
     $        DBX2*sgb(2,iatm)+DAX*sgc(2,iatm)
          swb(1,3,iatm)=gb2*gdxb(1,3,iatm)+gc*gdxa(1,3,Iatm)+
     $        DBX2*sgb(3,iatm)+DAX*sgc(3,iatm)
          swb(2,1,iatm)=gb2*gdxb(2,1,iatm)+gc*gdxa(2,1,Iatm)+
     $        DBY2*sgb(1,iatm)+DAY*sgc(1,iatm)
          swb(2,2,iatm)=gb2*gdxb(2,2,iatm)+gc*gdxa(2,2,Iatm)+
     $        DBY2*sgb(2,iatm)+DAY*sgc(2,iatm)
          swb(2,3,iatm)=gb2*gdxb(2,3,iatm)+gc*gdxa(2,3,Iatm)+
     $        DBY2*sgb(3,iatm)+DAY*sgc(3,iatm)
          swb(3,1,iatm)=gb2*gdxb(3,1,iatm)+gc*gdxa(3,1,Iatm)+
     $        DBZ2*sgb(1,iatm)+DAZ*sgc(1,iatm)
          swb(3,2,iatm)=gb2*gdxb(3,2,iatm)+gc*gdxa(3,2,Iatm)+
     $        DBZ2*sgb(2,iatm)+DAZ*sgc(2,iatm)
          swb(3,3,iatm)=gb2*gdxb(3,3,iatm)+gc*gdxa(3,3,Iatm)+
     $        DBZ2*sgb(3,iatm)+DAZ*sgc(3,iatm)
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    fock matrix.
c
c -- get the maximum absolute value of coefficients v and w
c
        Call absmax(NAt3,sva,iixx,svamax)
        Call absmax(NAt3,svb,iixx,svbmax)
        svmax=max(svamax,svbmax)
        Call absmax(3*NAt3,swa,iixx,swamax)
        Call absmax(3*NAt3,swb,iixx,swbmax)
        swmax=max(swamax,swbmax)
        tswmax=Three*swmax
c
c -- global threshold testing
c
        abr=max(abs(ra),abs(rb))
        vmx2=vmx*vmx
        SMax = max(Abs(SXXA+SXYA+SXZA),abs(SXXB+SXYB+SXZB))
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abr+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVIA=SXXA*ValIX+SXYA*ValIY+SXZA*ValIZ+ra*ValI
          XVIB=SXXB*ValIX+SXYB*ValIY+SXZB*ValIZ+rb*ValI
          SXXIA=SXXA*ValI
          SXYIA=SXYA*ValI
          SXZIA=SXZA*ValI
          SXXIB=SXXB*ValI
          SXYIB=SXYB*ValI
          SXZIB=SXZB*ValI
          XXVXA=SXXA*VAOXX(1,I)+SXYA*VAOXX(2,I)+SXZA*VAOXX(4,I)+ra*VALIX
          XXVYA=SXXA*VAOXX(2,I)+SXYA*VAOXX(3,I)+SXZA*VAOXX(5,I)+ra*VALIY
          XXVZA=SXXA*VAOXX(4,I)+SXYA*VAOXX(5,I)+SXZA*VAOXX(6,I)+ra*VALIZ
          XXVXB=SXXB*VAOXX(1,I)+SXYB*VAOXX(2,I)+SXZB*VAOXX(4,I)+rb*VALIX
          XXVYB=SXXB*VAOXX(2,I)+SXYB*VAOXX(3,I)+SXZB*VAOXX(5,I)+rb*VALIY
          XXVZB=SXXB*VAOXX(4,I)+SXYB*VAOXX(5,I)+SXZB*VAOXX(6,I)+rb*VALIZ
          valm=abs(vali)
          valm1=max(abs(xvia),abs(xvib))
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1a=max(abs(xxvxa),abs(xxvya),abs(xxvza))
          valmx1b=max(abs(xxvxb),abs(xxvyb),abs(xxvzb))
          valmx1=max(valmx1a,valmx1b)
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJA=SXXA*ValJX+SXYA*ValJY+SXZA*ValJZ
              XVJB=SXXB*ValJX+SXYB*ValJY+SXZB*ValJZ
              abxvj=max(abs(xvja),abs(xvjb))
              if(valmx1*abs(valj)+valmx*abxvj.gt.thrsh)then
              if(iatm.ge.nb.and.iatm.le.ne) then
                fda(1,iatm,ij)=fda(1,iatm,ij)-xxvxa*valj-valix*xvja
                fda(2,iatm,ij)=fda(2,iatm,ij)-xxvya*valj-valiy*xvja
                fda(3,iatm,ij)=fda(3,iatm,ij)-xxvza*valj-valiz*xvja
                fdb(1,iatm,ij)=fdb(1,iatm,ij)-xxvxb*valj-valix*xvjb
                fdb(2,iatm,ij)=fdb(2,iatm,ij)-xxvyb*valj-valiy*xvjb
                fdb(3,iatm,ij)=fdb(3,iatm,ij)-xxvzb*valj-valiz*xvjb
              endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
              if(jatm.ge.nb.and.jatm.le.ne) then
                fda(1,jAtm,ij) = fda(1,jatm,ij) - xvia*Valjx -
     $          vaoxx(1,j)*sxxia - vaoxx(2,j)*sxyia - vaoxx(4,j)*sxzia
                fda(2,jatm,ij) = fda(2,jatm,ij) - xvia*valjy -
     $          vaoxx(2,j)*sxxia - vaoxx(3,j)*sxyia - vaoxx(5,j)*sxzia
                fda(3,jatm,ij) = fda(3,jAtm,ij) - xvia*valjz -
     $          vaoxx(4,j)*sxxia - vaoxx(5,j)*sxyia - vaoxx(6,j)*sxzia
c
                fdb(1,jatm,ij) = fdb(1,jatm,ij) - xvib*valjx -
     $          vaoxx(1,j)*sxxib - vaoxx(2,j)*sxyib - vaoxx(4,j)*sxzib
                fdb(2,jatm,ij) = fdb(2,jatm,ij) - xvib*valjy -
     $          vaoxx(2,j)*sxxib - vaoxx(3,j)*sxyib - vaoxx(5,j)*sxzib
                fdb(3,jatm,ij) = fdb(3,jatm,ij) - xvib*valjz -
     $          vaoxx(4,j)*sxxib - vaoxx(5,j)*sxyib - vaoxx(6,j)*sxzib
              endif
              endif
              vij=vali*valj
              vijx=valix*valj+valjx*vali
              vijy=valiy*valj+valjy*vali
              vijz=valiz*valj+valjz*vali
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                do katm=nb,ne
                  fda(1,katm,ij) = fda(1,katm,ij) + (sva(1,katm)*vij+
     $                           swa(1,1,katm)*vijx+swa(2,1,katm)*vijy+
     $                           swa(3,1,katm)*vijz)
                  fda(2,katm,ij) = fda(2,katm,ij) + (sva(2,katm)*vij+
     $                           swa(1,2,katm)*vijx+swa(2,2,katm)*vijy+
     $                           swa(3,2,katm)*vijz)
                  fda(3,katm,ij) = fda(3,katm,ij) + (sva(3,katm)*vij+
     $                           swa(1,3,katm)*vijx+swa(2,3,katm)*vijy+
     $                           swa(3,3,katm)*vijz)
c
                  fdb(1,katm,ij) = fdb(1,katm,ij) + (svb(1,katm)*vij+
     $                           swb(1,1,katm)*vijx+swb(2,1,katm)*vijy+
     $                           swb(3,1,katm)*vijz)
                  fdb(2,katm,ij) = fdb(2,katm,ij) + (svb(2,katm)*vij+
     $                           swb(1,2,katm)*vijx+swb(2,2,katm)*vijy+
     $                           swb(3,2,katm)*vijz)
                  fdb(3,katm,ij) = fdb(3,katm,ij) + (svb(3,katm)*vij+
     $                           swb(1,3,katm)*vijx+swb(2,3,katm)*vijy+
     $                           swb(3,3,katm)*vijz)
                enddo
              endif
c
c -- weight derivatives contribution
c
              dxija=dax*vijx+day*vijy+daz*vijz
              dxijb=dbx*vijx+dby*vijy+dbz*vijz
              gija=prap*vij+pgap2*dxija+pgcp*dxijb
              gijb=prbp*vij+pgbp2*dxijb+pgcp*dxija
              if(gwtmax*max(abs(gija),abs(gijb)).gt.thrsh)then
                do ia=nb,ne
                  fda(1,ia,ij)=fda(1,ia,ij)+gwt(1,ia,ipp)*gija
                  fda(2,ia,ij)=fda(2,ia,ij)+gwt(2,ia,ipp)*gija
                  fda(3,ia,ij)=fda(3,ia,ij)+gwt(3,ia,ipp)*gija
                  fdb(1,ia,ij)=fdb(1,ia,ij)+gwt(1,ia,ipp)*gijb
                  fdb(2,ia,ij)=fdb(2,ia,ij)+gwt(2,ia,ipp)*gijb
                  fdb(3,ia,ij)=fdb(3,ia,ij)+gwt(3,ia,ipp)*gijb
                enddo
              endif
 230        continue
          endif
 240    continue
 245    continue
        call secund(t5)
        tqf=tqf+t5-t4
 200    continue
      endif
c
      return
      end
