      subroutine loca(nmo,meth,iprscf,emin,fock,coef)

      use memory

      implicit real*8 (a-h,o-z)
c  Arguments:
c  nmo=number of occupied MOs
c  meth=localization method, either 'Boys' or 'Pipek'
c  iprsc= SCF print flag
c  emin= the lowest orbital energy still counted as valence
c   orbital and included in the localization. Default (in rhf5) is -3.0
c  fock=the NAME of the Fock matrix (i.e., alpha or beta spin)
c  coef=the NAME of the SCF coefficient matrix to be localized
      character*(*) meth
      character*4 method
      character*(*) fock,coef
      logical exchanged
c     common /big/bl(10000)
c     common /intbl/maxsh,inx(100)
c   automatically allocated variable
      integer order(nmo)
c
      call mmark
c
      method=meth(1:4)
      if(meth(1:1).eq.'b') meth(1:1)='B'
      if(meth(1:1).eq.'p') meth(1:1)='P'
      call getival('ncf',ncf)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('ncs',ncs)
      call getival('ictr',ictr)
      call getival('na',na)
c determine the first valence orbital
      iev=mataddr('diag')
      ifirst=1
      do i=1,nmo
        if(bl(iev+i-1).lt.emin) then
          ifirst=ifirst+1
        end if
      end do
c no valence obitals if ifirst>nmo
      if(ifirst.gt.nmo) then
        call retmark
        return
      end if
c  write only the orbitals if there is only 1 valence orbital
      if(nmo-ifirst.lt.1) go to 50
      call lowerca2(method,4)
      if(method.eq.'boys') then
        call boys(coef,ifirst,nmo,nrot)
      else if(method.eq.'pipe') then
        call matdef('smat','s',ncf,ncf)
        ism=mataddr('smat')
c   place for the h matrix
c
c  parameters: itype: 0=overlap, 1=h0 mtx, see the rest in the routine
c  na=number of atoms,integrals,contraction info,2 directions for
c  x,y,z-dependent quantities,basis set info,nucl..info, numer of
c  contracted shells
c  determine the overlap and h0 matrices
        call inton(0,na,bl(ism),bl(ictr),0,0,bl(ibas),bl(inuc),ncs)
        call pipek(coef,ifirst,nmo,nrot)
      else
        return
      end if
  50  continue
c
      call getival('iout',iout)
      call matdef('xmo','d',nmo,nmo)
      call matdef('ymo','d',nmo,nmo)
      call matdef('zmo','d',nmo,nmo)
c Calculate Coulson orb. energies (diagonal Fock matrix elements)
      call matdef('eloc','d',nmo,nmo)
      call matsub('occu',coef,1,nmo)
      call matsimtr(fock,'occu','eloc')
c
c  Sort the localized orbitals by increasing orb. energy
CPP
c     call matprint('eloc',6)
c first fill up order (integer!) with 1,2,...,nmo
      do i=1,nmo
        order(i)=i
      end do
   80 continue
      exchanged=.false.
      do i=ifirst+1,nmo 
        call matelem('eloc',i,i,ei)
c        call matsub('icolumn',coef,i,i)
        do j=ifirst,i
          call matelem('eloc',j,j,ej)
          if(ej.gt.ei) then
            exchanged=.true.
            call mateset('eloc',i,i,ej)
            call mateset('eloc',j,j,ei)
            iii=order(i)
            order(i)=order(j)
            order(j)=iii
c            call matsub('jcolumn',coef,j,j)
c            call matcopy('icolumn','tempxyz9')
c            call matcopy('jcolumn','icolumn')
c            call matcopy('tempxyz9','jcolumn')
c            call matrem('jcolumn')
            exit
          end if
        end do
c        call matrem('icolumn')
      end do
      if(exchanged) go to 80 
c
c  Now reorder the eigenvectors
      call matdef('tempmocp','r',ncf,nmo-ifirst+1)
      call matsub('valence',coef,ifirst,nmo)
      call matcopy('valence','tempmocp')
      call matrem('valence')
      do i=ifirst,nmo
        j=order(i)
        j1=j-ifirst+1
c  new(i)=old(order(i))
        call matsub('jcolumn','tempmocp',j1,j1)
        call matsub('icolumn',coef,i,i)
        call matcopy('jcolumn','icolumn') 
        call matrem('icolumn')
        call matrem('jcolumn')
      end do
      call matrem('tempmocp')
CPP
c     call matprint('eloc',6)
c  Write out orbital coefficients and localized "eigenvalues"
      np4=4
      if(coef(1:5).ne.'coefB') then
cc        write(6,*) ' Writing Alpha localized orbitals'
        call matwrite(coef,np4,0,'loca_rhf')
        call matwrite('eloc',np4,0,'eloc_rhf')
       else
cc        write(6,*) ' Writing Beta localized orbitals'
        call matwrite(coef,np4,0,'loca_uhf')
        call matwrite('eloc',np4,0,'eloc_uhf')
      end if
c
      call matsimtr('dipx','occu','xmo')
      call matsimtr('dipy','occu','ymo')
      call matsimtr('dipz','occu','zmo')
      call matdef('rsquare','s',ncf,ncf)
      ir2=mataddr('rsquare')
      call matdef('r2mo','d',nmo,nmo)
      call inton(4,na,bl(ir2),bl(ictr),0,0,bl(ibas),bl(inuc),ncs)
      call matsimtr('rsquare','occu','r2mo')
c  write out these quantities
      if(coef(1:5).ne.'coefB') then
        call matwrite('xmo',np4,0,'xcor_loA')
        call matwrite('ymo',np4,0,'ycor_loA')
        call matwrite('zmo',np4,0,'zcor_loA')
        call matwrite('r2mo',np4,0,'rr_loA')
      else
        call matwrite('xmo',np4,0,'xcor_loB')
        call matwrite('ymo',np4,0,'ycor_loB')
        call matwrite('zmo',np4,0,'zcor_loB')
        call matwrite('r2mo',np4,0,'rr_loB')
      endif
c
      ieloc=mataddr('eloc')
      iwf=mataddr(coef)
      ixmo=mataddr('xmo')
      iymo=mataddr('ymo')
      izmo=mataddr('zmo')
      ir2mo=mataddr('r2mo')
      call getrval('angs',angs)
      if(iprscf.gt.3) call druwf(bl(ieloc),coef,0.0d0,0,0,nmo,iout)
      write(iout,100) meth(1:5),nrot
c      write(icond,100) meth(1:5),nrot
 100  format(
     1   /1x,'The ',a,' localized orbitals after ',i0,' rotations:'
     2  //6x,'  Energy (Eh)',15x,'Centroid (Ang)',13x'Radius (Ang)'/)
      do i=1,nmo
        eavr=bl(ieloc+i-1)
        xc=bl(ixmo+i-1)
        yc=bl(iymo+i-1)
        zc=bl(izmo+i-1)
        rr=sqrt(bl(ir2mo+i-1)-xc**2-yc**2-zc**2)
c       write(icond,200) i,eavr,xc/angs,yc/angs,zc/angs,rr/angs
        write(iout,200)  i,eavr,xc/angs,yc/angs,zc/angs,rr/angs
      end do
 200  format(i6,f13.7,3f13.7,f13.7)
      call matrem('r2mo')
      call matrem('rsquare')
      call matrem('occu')
      call matrem('eloc')
      call matrem('zmo')
      call matrem('ymo')
      call matrem('xmo')
      call matsear('smat',ifound)
      if(ifound.ne.0) call matrem('smat')
      call retmark
      end
c
c----------------------------------------------------------------------
      subroutine boys(coef,ifirst,nmo,nrot)

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) coef
c  this routine performs the Boys localization.
c  ARGUMENTS
c  coef = NAME of the matrix containing the SCF coefficients  
c  ifirst = the first occupied orbital included in the localization
c  nmo = total number of occupied MOs
c  it takes its other arguments from the matrix system:
c  'dipx', 'dipy', 'dipz'= matrix of the X, Y and Z operators in AO basis
c  INTENT(OUT)
c  nrot = number of orbital rotations needed
      parameter (thres=5.0d-8)
c this number is the tangent of pi/8=22.5 degrees
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
c     common /big/bl(10000)
      call getrval('one',one)
      call getrval('zero',zero)
      call matsub('cval',coef,ifirst,nmo)
      mval=nmo+1-ifirst
      call matdef('xmo','s',mval,mval)
      call matdef('ymo','s',mval,mval)
      call matdef('zmo','s',mval,mval)
      ixmo=mataddr('xmo')
      iymo=mataddr('ymo')
      izmo=mataddr('zmo')
      call matsimtr('dipx','cval','xmo')
      call matsimtr('dipy','cval','ymo')
      call matsimtr('dipz','cval','zmo')
      call matdef('loctr','q',mval,mval)
      loctr=mataddr('loctr')
      call matzero('loctr')
      do i=1,mval
        call mateset('loctr',i,i,one)
      end do
c return label for the iteration
      nrot=0
 200  continue
c first establish the max. rotational angle
      rotmax=zero
      ii=0
      do i=2,mval
        ii=ii+i-1
       jj=0
        do j=1,i-1
         jj=jj+j-1
        call locrot(i,j,ii,jj,mval,bl(ixmo),bl(iymo),bl(izmo),s)
c  the routine locrot determines tan(4*alpha)=t4 for the 2x2 rotation
c  between orbitals i and j which minimizes the sum of squared orbital
c  radii
          if(abs(s).gt.rotmax) rotmax=abs(s)
        end do
      end do
c  now do the same thing again but rotate the orbitals
      ii=0
      do i=2,mval
        ii=ii+i-1
        icol=loctr+(i-1)*mval
      jj=0
        do j=1,i-1
        jj=jj+j-1
        jcol=loctr+(j-1)*mval
        call locrot(i,j,ii,jj,mval,bl(ixmo),bl(iymo),bl(izmo),s)
        if(abs(s).gt.0.2d0*rotmax) then
          nrot=nrot+1
            c=sqrt(one-s**2)
          call drot(mval,bl(icol),1,bl(jcol),1,c,-s)
c rotate the X,Y Z matrices, too
             call rotsym(mval,bl(ixmo),i,j,ii,jj,c,s)
             call rotsym(mval,bl(iymo),i,j,ii,jj,c,s)
             call rotsym(mval,bl(izmo),i,j,ii,jj,c,s)
          end if
      end do
      end do
      if(rotmax.gt.thres.and.nrot.lt.1000000) go to 200
      call matmmult('cval','loctr','cval')
      call matrem('loctr')
      call matrem('zmo')
      call matrem('ymo')
      call matrem('xmo')
      call matrem('cval')
      end
c
c----------------------------------------------------------------------
      subroutine locrot(i,j,ii,jj,m,x,y,z,s)
      implicit real*8 (a-h,o-z)
c  Input parameters
c  i,j= rotation indices; i>j
c  ii=i*(i-1)/2; jj=j*(j-1)/2
c  m=dimension (number of orbitals to be localized)
c  x,y and z: the matrix of the x,y and z operators in MO basis
c  Output parameter:
c  t4=tan[4 alpha] where alpha is the rotational angle which maximizes
c  localization using the rotation (i,j)
c  s=result: sine of the rotation angle
      parameter (zero=0.0d0,one=1.0d0,pi2=1.570796327d0,thres=5.0d-8)
c  this number is the tangent of pi/8=22.5 degrees
      data tpi8/0.414213562373095d0/
      dimension x(*),y(*),z(*)
c     ii=i*(i-1)/2
c     jj=j*(j-1)/2
      xii=x(ii+i)
      xjj=x(jj+j)
      xij=x(ii+j)
      yii=y(ii+i)
      yjj=y(jj+j)
      yij=y(ii+j)
      zii=z(ii+i)
      zjj=z(jj+j)
      zij=z(ii+j)
      tt=xij*(xii-xjj)+yij*(yii-yjj)+zij*(zii-zjj)
      den=xij**2+yij**2+zij**2-0.25d0*((xii-xjj)**2+(yii-yjj)**2
     1 +(zii-zjj)**2)
      if(abs(den).lt.1.0d-6) then
        if(tt.gt.zero) then
        t4=-tpi8
        else
        t4=tpi8
        end if
      else
        t4=tt/den
      end if
        at4=abs(t4)
          if(at4.lt.0.20d0) then
c  use the power series expansion Sin[ArcTan[x]/4,{x,0,6}]=x/4-11 x^3/128
c  +431 x^5/8192 + O[x^7]
            s=0.25d0*t4-0.0859375d0*t4**3
            if(den.gt.zero) then
              c=sqrt(one-s**2)
                s=0.7071067812d0*(c+s)
c  this means that we had a mximum instead of a minimum and increase
c  4alpha by pi, i.e. alpha by pi/4
              end if
          else
            a4=atan(t4)
            c4=cos(a4)
            s4=sin(a4)
            s=sin(a4*0.25d0)
            if(c4*den+s4*tt.gt.zero) s=sin(a4*0.25d0+pi2)
          end if
      end
c
c----------------------------------------------------------------------
      subroutine rotsym(m,x,i,j,ii,jj,c,s)
      implicit real*8 (a-h,o-z)
c  rotate a symmetrical matrix
c   Input:
c   m=dimension
c   x: matrix to rotate, stored as upper triangle columnwise
c   i,j: columns & rows to rotate
c   ii=i*(i-1)/2; jj=j*(j-1)/2
c   c=cosine of the angle
c   s=sine of the angle
c     rotation: |i'>=c|i>-s|j>; |j'>=s|i>+c|j>
c   Output:
c   x=rotated matrix

      parameter(two=2.0d0)
      dimension x(*)
      cs=c*s
      cs2=two*cs
      xii=x(ii+i)
      xjj=x(jj+j)
      xij=x(ii+j)
      i1=ii+1
      j1=jj+1
      do k=1,m
        z=x(i1)*c-x(j1)*s
      zz=x(i1)*s+x(j1)*c
      x(i1)=z
      x(j1)=zz
      i1=i1+1
      j1=j1+1
      if(k.ge.i)i1=i1+k-1
      if(k.ge.j)j1=j1+k-1
      end do
      x(ii+i)=c**2*xii+s**2*xjj-cs2*xij
      x(jj+j)=s**2*xii+c**2*xjj+cs2*xij
      x(ii+j)=(xii-xjj)*cs+xij*(c**2-s**2)
      end
c
c----------------------------------------------------------------------
      subroutine pipek(coef,ifirst,nmo,nrot)

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) coef
c  this routine performs the Pipek-Mezey localization.
c  ARGUMENTS
c  INTENT(IN)
c  coef= NAME of the matrix containing the SCF coefficients
c  ifirst = first occupied orbital included in the localization
c  nmo=number of occupied MOs
c  it takes its other arguments from the matrix system:
c  'dipx', 'dipy', 'dipz'= matrix of the X, Y and Z operators in AO basis
c  INTENT(OUT)
c  nrot = number of orbital rotations
c
      parameter (thres=5.0d-8)
c this number is the tangent of pi/8=22.5 degrees
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
c     common /big/bl(10000)
c     common /intbl/maxsh,inx(100)
c number of valence orbitals
      mval=nmo-ifirst+1
      call getrval('one',one)
      call getrval('zero',zero)
      call getival('ncf',ncf)
      call getival('ictr',ictr)
      call getival('na',na)
      call matsub('cval',coef,ifirst,nmo)
      icoef=mataddr('cval')
      call matdef('SxC','r',ncf,mval)
      isxc=mataddr('SxC')
      call matmmult('smat','cval','SxC')
c reserve memory for qia, qja, qija
      call getmem(na,iqi)
      call getmem(na,iqj)
      call getmem(na,iqij)
      nrot=0
c return label for the iteration
 200  continue
c  establish the max. rotational angle
      rotmax=zero
      ii=0
      imval=0
      do i=2,mval
      ii=ii+i-1
      imval=imval+ncf
c       imval=(i-1)*ncf
      icol=icoef+imval
      icol1=isxc+imval
      jj=0
      jmval=-ncf
      do j=1,i-1
        jj=jj+j-1
          jmval=jmval+ncf
c         jmval=(j-1)*ncf
        jcol=icoef+jmval
        jcol1=isxc+jmval
        call piprot(i,j,mval,ncf,na,bl(icoef),bl(isxc),s,
     1     bl(ictr),bl(iqi),bl(iqj),bl(iqij))
        if(abs(s).gt.rotmax)rotmax=abs(s)
        if(abs(s).gt.0.2d0*rotmax) then
          nrot=nrot+1
          c=sqrt(one-s**2)
c  rotate the matrices C and SC
          call drot(ncf,bl(icol),1,bl(jcol),1,c,s)
          call drot(ncf,bl(icol1),1,bl(jcol1),1,c,s)
        end if
      end do
      end do
      if(rotmax.gt.thres.and.nrot.lt.1000000) go to 200
      call retmem(3)
      call matrem('SxC')
      call matrem('cval')
      end
c
c----------------------------------------------------------------------
      subroutine piprot(i,j,mval,ncf,na,cmat,sc,s,inx,qi,qj,qij)
c This routine determines the sine of the rotational angle for the Pipek-Mezey
c  localization. cmat is the SCF coefficient matrix, SC =  SC, S is the
c  overlap matrix, s is the sine of the rotational angle
      implicit real*8 (a-h,o-z)
c  Input parameters
c  i,j= rotation indices; i>j
c  mval=dimension (number of orbitals to be localized)
c  ncf=number of contracted functions
c  ncs=number of contracted shells
c  na=number of nuclei
c  cmat=SCF coefficients; sc=SC; S=overlap matrix
c  qia,qja,qija: Storage
c  Output parameter:
c  s=Sin[ArcTan[4 alpha]/4] where alpha is the rotational angle which maximizes
c  localization using the rotation (i,j)
      parameter (zero=0.0d0,one=1.0d0,pi2=1.570796327d0,thres=5.0d-8)
c  this number is the tangent of pi/8=22.5 degrees
      data tpi8/0.414213562373095d0/
      dimension cmat(ncf,mval),sc(ncf,mval),inx(12,*),qi(na),qj(na),
     1 qij(na)
c     ii=i*(i-1)/2
c     jj=j*(j-1)/2
c loop over the AOs; determine to which atom do they belong
c  zero out qi, qj, qij
c  number of contracted shells
      call getival('ncs',ncs)
c very small number
      call getrval('acc',acc)
      call zeroit(qi,na)
      call zeroit(qj,na)
      call zeroit(qij,na)
      icf=0
c cycle through contracted (possibly GENERALLY CONTRACTED) shells
      do ics=1,ncs
c number of general contractions-1 (zero for segmented contraction)
        ngc=inx(4,ics)
c atom associated with this contr. shell
      nat=inx(2,ics)
      qia=zero
      qja=zero
      qija=zero
c length of the shell: 1 for s, 3 for p, 4 for sp, 4 for d, 6 for d6...
      length=inx(3,ics)
        do igc=1,ngc+1
          do ifun=1,length
           icf=icf+1
           qia=qia+cmat(icf,i)*sc(icf,i)
           qja=qja+cmat(icf,j)*sc(icf,j)
           qija=qija+cmat(icf,i)*sc(icf,j)+cmat(icf,j)*sc(icf,i)
60      format(3i3,4f9.6,2x,3f9.6)
50      format('i,j,icf,c(icf,i),c(icf,j),sc(icf,i),sc(icf,j),qia,qja,
     1  qija')
        end do
        end do
      qi(nat)=qi(nat)+qia
      qj(nat)=qj(nat)+qja
c  this qij is the double of the quantity in the Boughton paper
      qij(nat)=qij(nat)+qija
      end do
      a=zero
      b=zero
      do nat=1,na
        a=a+(qi(nat)-qj(nat))**2-qij(nat)**2
        b=b+(qi(nat)-qj(nat))*qij(nat)
      end do
c  form Mathematica:
c the function to be maximized is const+1/4(a*Cos(4*phi)+2*b*Sin[4*phi])
c Its derivative is -a*sin[4*phi]+2*b*Cos[4*phi]  = 0
c Tan[4*phi]=2b/a
c The second drivative at phi=0 is -4*a. It should be negative for
c a maxmum. If it is positive, then we need the next solution of Tan:
c increase 4*phi by Pi or Phi by Pi/4
c  test for minimum or maximum
      if(abs(a).lt.acc) then
        if(abs(b).gt.acc) then
          t4=tpi8
        else
          t4=zero
      end if
      else
          t4=(b+b)/a
      end if
        at4=abs(t4)
          if(at4.lt.0.20d0) then
c  use the power series expansion Sin[ArcTan[x]/4,{x,0,6}]=x/4-11 x^3/128
c  +431 x^5/8192 + O[x^7]
            s=0.25d0*t4-0.0859375d0*t4**3
          else
            a4=atan(t4)
            c4=cos(a4)
            s4=sin(a4)
            s=sin(a4*0.25d0)
          end if
         if(a.lt.zero) then
         c=sqrt(one-s**2)
           s=0.7071067812d0*(c+s)
c  this means that we had a minimum instead of a maximum and increase
c  4alpha by pi, i.e. alpha by pi/4
       end if
      end
