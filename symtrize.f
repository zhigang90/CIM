      subroutine symtrize(na,xnuc,iprnt,symthr,nsym,nsy)
c  this routine finds the molecular symmetry or approximate
c  symmetry, and symmetrizes the nuclear coordinates if the
c  symmetry is accurate to within symthr. If symthr is zero,
c  it returns with nsym=0
c  parameters:
c  INPUT:
c  na=number of nuclei
c  xnuc(5,na)= nuclear data: q,x,y,z,name for each atom
c  iprnt= print flag
c  symthr= threshold
c  OUTPUT:
c  nsym= number of symmetry elements found
c  nsy(k)= the nature of the k-th symmetry element. The symmetrey
c  elements are numbered according to the coordinates which change sign:
c  1=x=sigma(yz); 2=y=sigma(xz); 3=z=sigma(xy); 4=xy=C2(z); 5=xz=C2(y);
c  6=yz=C2(x); 7=xyz=i

      use memory

      implicit real*8 (a-h,o-z)
      logical nopair
      character*2 csyop(7)
      character*3 groups(8),group,group1
c ..................................................
c -- automatic allocation of arrays in F90
      real*8 xyz(3,na),qq(na)
c ..................................................
      parameter (nsyop=7)
      dimension xnuc(5,na),nsy(nsyop),ipair(nsyop),devmax(nsyop)
      dimension xx(3)
c  mirror gives the symmetry behavior of x,y,z; ns12 is the group
c  multiplication table
      dimension mirror(3,nsyop),ns12(nsyop,nsyop)
c     common /intbl/ maxsh,inx(100)
c
      data mirror(1,1),mirror(2,1),mirror(3,1),mirror(1,2),mirror(2,2),m
     1irror(3,2),mirror(1,3),mirror(2,3),mirror(3,3),mirror(1,4),mirror(
     22,4),mirror(3,4),mirror(1,5),mirror(2,5),mirror(3,5),mirror(1,6),m
     3irror(2,6),mirror(3,6),mirror(1,7),mirror(2,7),mirror(3,7)/-1,1,1,
     41,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
      data ns12 /0,0,0,0,0,0,0, 4,0,0,0,0,0,0, 5,6,0,0,0,0,0,
     1  2,1,7,0,0,0,0, 3,7,1,6,0,0,0, 7,3,2,5,4,0,0, 6,5,4,3,2,1,0/
      data csyop /'x ','y ','z ','xy','xz','yz','i '/
      data zero/0.0d0/
      data groups/'d2h','d2 ','c2v','c2h','c2 ','cs ','ci ','c1 '/
c  long output file
      call getival('iout',iout)
c      call getival('icond',icond)
      nsym=0
c if symthr=0 then symmetry is supressed
      if(symthr.eq.zero) then
        go to 1000
      end if
c
      do i=1,na
        xyz(1,i)=xnuc(2,i)
        xyz(2,i)=xnuc(3,i)
        xyz(3,i)=xnuc(4,i)
        qq(i)=xnuc(1,i)
      end do
c this is a pre-symmetrizer - finds planes more efficiently
c than the main symmetrizer
      call center_of_mass(na,qq,xyz)
      call sym_rot(na,qq,xyz,symthr)
      do i=1,na
        xnuc(2,i)=xyz(1,i)
        xnuc(3,i)=xyz(2,i)
        xnuc(4,i)=xyz(3,i)
      end do
c
   50 continue
      epsmax=zero
      npairt=0
c      this is set true if there is an atom with no pair under any symm.
c  devmax is the maximum deviation from symmetry. It is used to pinpoint
c  almost symmetrical atoms in a checkrun
        call zeroit(devmax,nsyop)
        do 150 k=1,nsyop
           nopair=.false.
           npair=0
           do 140 i=1,na
              do 110 l=1,3
                 xx(l)=xnuc(l+1,i)
                 if (mirror(l,k).lt.0) xx(l)=-xnuc(l+1,i)
  110         continue
              ipairk=0
         epsmin=1000.0d0
              do 130 j=1,na
                 rr=zero
                 do 120 l=1,3
                    rr=rr+abs(xx(l)-xnuc(l+1,j))
  120            continue
                 if (rr.lt.symthr) then
                   ipairk=j
                   npair=npair+1
                   npairt=npairt+1
                   if(rr.gt.epsmax) epsmax=rr
                 end if
            if (rr.lt.epsmin) epsmin=rr
  130         continue
              if(epsmin.gt.devmax(k)) devmax(k)=epsmin
c  devmax(k) = max(i) min(j) [R(i)-j]
c  for exact symmetry, this is zero
              if(ipairk.eq.0) nopair=.true.
  140      continue
           if(.not.nopair) then
                  nsym=nsym+1
         nsy(nsym)=k
               end if
           if(nopair.and.devmax(k).lt.0.2d0) then
        write(iout,200) csyop(k),devmax(k)
 200         format(' approximate symmetry ',a2,2x,f12.6)
           end if
  150    continue
c   make sure that the symmetries form a group
      do 160 k=1,nsym
        ns1=nsy(k)
        do 160 l=k+1,nsym
          ns2=nsy(l)
          ns3=ns12(ns1,ns2)
          if(ns3.ne.0) then
            do 170 m=1,nsym
              if(nsy(m).eq.ns3) go to 160
 170        continue
            nsym=nsym+1
            nsy(nsym)=ns3
          end if
 160  continue
ccc
      If(iprnt.gt.2) Then
        write(iout,900) (csyop(nsy(k)),k=1,nsym)
 900    format(' atomic symmetry pairs',/,13x,7(a2,3x))
      EndIf
ccc
c     reserve integer memory for the atomic symmetry pairs
      call getint(na*nsym,nupair)
      call setival('SymNuPr',nupair)
      call setival('SymNuPr1',nupair)     ! JB  19 Mar 99  NEEDS CHECKING ***
c  the atomic nuclear pairs are stored as an array
c  ((nupair(i,k), i=1,na), k=1,nsym)
c   subscript for the array nupair
        do 250 i=1,na
           npair=0
           do 240 k=1,nsym
              kk=nsy(k)
              do 210 l=1,3
                 xx(l)=xnuc(l+1,i)
                 if (mirror(l,kk).lt.0) xx(l)=-xnuc(l+1,i)
  210         continue
              ipair(k)=0
              do 230 j=1,na
                 rr=zero
                 do 220 l=1,3
                    rr=rr+abs(xx(l)-xnuc(l+1,j))
  220            continue
                 if (rr.lt.symthr) then
                   ipair(k)=j
                   npair=npair+1
                   npairt=npairt+1
                 end if
  230         continue
  240      continue
c          now force-symmetrize the coordinates of atoms which are
c          pair to the present one
              do 242 l=1,3
                 xx(l)=xnuc(l+1,i)/dble(npair+1)
  242         continue
           do 245 k=1,nsym
             kk=nsy(k)
             j=ipair(k)
             if(j.eq.0) go  to 245
             do 244 l=1,3
               xx(l)=xx(l)+xnuc(l+1,j)*mirror(l,kk)/dble(npair+1)
 244         continue
 245       continue
           do 248 k=1,nsym
              kk=nsy(k)
              j=ipair(k)
              if(j.eq.0) go to 248
              do 247 l=1,3
                xnuc(l+1,i)=xx(l)
                xnuc(l+1,j)=xx(l)*mirror(l,kk)
 247          continue
c   set nuclear pair info - needed in gradients etc. for symmetry,
c   also in NMR
c            inup=nupair+i+(k-1)*na-1
c            bl(inup)=ipair(k)
             call int_to_bl( bl(nupair), i+(k-1)*na, ipair(k) )
 248       continue
c
 249       continue
ccc
           If(iprnt.gt.2) Then
           write (iout,910) i,(ipair(k),k=1,nsym)
  910      format (1x,i3,6x,7i5)
           EndIf
ccc
  250   continue
 1000   continue
      call setival('nsym',nsym)
      call getint(7,iadr)
      call izeroit(bl(iadr), 7 )
      do 1100 i=1,7
c      inx(iadr+i-1)=0
c      if(i.le.nsym) inx(iadr+i-1)=nsy(i)
       if(i.le.nsym) call int_to_bl( bl(iadr), i, nsy(i) )
 1100 continue
c
c---------------------------------------------------------
c the symmetry elements are in common /symm/ but they are also
c stored as options (nsym and the starting address of nsy)
c
      call setival('nsyo',iadr)
c---------------------------------------------------------
c determine the number of group generators :
c
      ngener=nsym
      if(nsym.gt.0) then
         if(nsym.eq.3) ngener=2
         if(nsym.eq.7) ngener=3
      endif
c
c stored ngener as a option
c
      call setival('ngener',ngener)
c---------------------------------------------------------
      if(nsym.gt.0) then
c         write(icond,1200)symthr,(csyop(nsy(k)),k=1,nsym)
         write(iout,1200)symthr,(csyop(nsy(k)),k=1,nsym)
      else
c         write(icond,1300)  symthr
      endif
c
 1200 format('Symmetry found at threshold ',e12.5,1x,
     1 ': ',7(a3,1x))
 1300 format('No symmetry found at threshold', 1x,e12.5)
      if(nsym.eq.7) then
c  D2h
        group=groups(1)
      else if(nsym.eq.3) then
        if(nsy(1).le.3.and.nsy(2).le.3) then
c generators are two symmetry planes - C2v
          group=groups(3)
        else if(nsy(1).le.3.or.nsy(2).le.3) then
c  generators are one plane and one C2 - C2h
          group=groups(4)
        else if(nsy(1).gt.3.and.nsy(2).gt.3.and.nsy(2).lt.7) then
c  generators are two C2 axes - D2
          group=groups(2)
        end if
      else if(nsym.eq.1) then
        if(nsy(1).le.3) then
c one plane only - Cs
          group=groups(6)
        else if(nsy(1).lt.7) then
c one C2 axis
          group=groups(5)
        else if(nsy(1).eq.7) then
          group=groups(7)
        end if
      else
         group=groups(8)
      end if
      if(iprnt.gt.0) then
        group1=group
        call uppercase(group1,1)
        write(iout,1400) group1
      end if
c      if(iprnt.gt.1) write(icond,1400) group
 1400 format('Largest Abelian subgroup of the mol. point group is ',
     1 a3)
c  store the pont group information
      call setchval('Ptgroup',group)
c---------------------------------------------------------
      end
c============================================================
      subroutine sym_rot(na,q,xnuc,tol)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,
     1 pi=3.14159 26535 89793d0)
      logical found,testplane
      dimension q(na),xnuc(3,na),plane(3),zz(3),vv(3)
c Try first the coordinate planes
      nfound=0
      do ll=1,3
        plane(1)=zero
        plane(2)=zero
        plane(3)=zero
        plane(ll)=one
        if(testplane(na,q,xnuc,tol,plane)) then
          nfound=nfound+1
        end if
      end do
      if(nfound.ge.2) then
        return
      end if
      call symplane(na,q,xnuc,tol,1,plane,found)
c      print *, found,plane
      if(.not.found) return
c xx is first the z axis
      zz(1)=zero
      zz(2)=zero
      zz(3)=one
      cc=scalar(zz,plane)
      zpl=abs(arcos(cc))
c      print *, 'zpl',zpl
      if(abs(cc).lt.0.995d0) then
      call normal(plane,zz,vv)
c      print *, 'vv', vv
c plane is the normal to the mirror plane, vv is perpendicular to the z axis
c and the normal to the plane
c3456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
c xx is the x unit vector
c angle between the x axis and the normal
      if(abs(vv(2)).gt.1.0d-9) then
        if(abs(vv(1)).gt.1.0d-9) then
            angle=atan(vv(2)/vv(1))
        else
          angle=pi*half
        end if
      else
        angle=zero
      end if
c      print *, 'angle', angle
c first rotate the whole thing around the z axis by -angle
c      print *, 'before rotation'
c      call printnuc(na,q,xnuc)
      call rotmol(na,xnuc,3,-angle)
c      print *, 'after first rot.'
c      call printnuc(na,q,xnuc)
c rotate the symm. plane as well
      call rotmol(1,plane,3,-angle)
      if(plane(3).lt.zero) then
        plane(1)=-plane(1)
        plane(2)=-plane(2)
        plane(3)=-plane(3)
      end if
c      print *, 'plane',plane
c next rotate around x to get the plane to the z axis
        cc=scalar(zz,plane)
        zpl=arcos(cc)
      if(plane(2).lt.zero) zpl=-zpl
        call rotmol(na,xnuc,1,zpl)
        call rotmol(1,plane,1,zpl)
c      print *, 'plane2',plane
      end if
c      print *, 'after second rot.'
c      call printnuc(na,q,xnuc)
      call symmetrizz(na,xnuc,tol,3)
c      print *, 'after symmetrization'
c      call printnuc(na,q,xnuc)
c look for a symmetry plane perpendicular to the z axis
      call symplane(na,q,xnuc,tol,2,plane,found)
c      print *, 'second trial', found,plane
      if(.not.found) return
c eliminate the z component of the vector
      plane(3)=zero
      call nom1(plane)
      if(abs(plane(2)).gt.1.0d-9) then
        angle=atan(plane(2)/plane(1))
      else
        angle=zero
      end if
c      print *, 'angle', angle
      call rotmol(na,xnuc,3,-angle)
c      call printnuc(na,q,xnuc)
      call symmetrizz(na,xnuc,tol,1)
      call symmetrizz(na,xnuc,tol,2)
c after second symmetrization step
c      call printnuc(na,q,xnuc)
      end
c=======================================================================================
      subroutine symmetrizz(na,xnuc,tol,iaxis)
      implicit real*8 (a-h,o-z)
      parameter(half=0.5d0)
      dimension xnuc(3,na),xx(3),yy(3)
      do i=1,na
      if(iaxis.eq.3) then
        xx(1)=xnuc(1,i)
        xx(2)=xnuc(2,i)
        xx(3)=-xnuc(3,i)
      else if(iaxis.eq.1) then
        xx(1)=-xnuc(1,i)
        xx(2)=xnuc(2,i)
        xx(3)=xnuc(3,i)
      else if(iaxis.eq.2) then
        xx(1)=xnuc(1,i)
        xx(2)=-xnuc(2,i)
        xx(3)=xnuc(3,i)
      end if
      rr=sqrt(xx(1)**2+xx(2)**2+xx(3)**2)
      tol1=rr*tol
      if(tol.gt.tol1) tol1=tol
c       print *,'tol1',tol1,xx
      do j=1,i
        yy(1)=xnuc(1,j)
        yy(2)=xnuc(2,j)
        yy(3)=xnuc(3,j)
c         print *,'dist.',i,j,distance(xx,yy)
        if(distance(xx,yy).lt.tol1) then
c          print *, 'averaging ', i,j
          x=(xx(1)+yy(1))*half
          y=(xx(2)+yy(2))*half
          z=(xx(3)+yy(3))*half
          xnuc(1,i)=x
          xnuc(1,j)=x
          xnuc(2,i)=y
          xnuc(2,j)=y
          xnuc(3,i)=-z
          xnuc(3,j)=z
          if(iaxis.eq.1) then
            xnuc(1,i)=-x
            xnuc(3,i)=z
          else if(iaxis.eq.2) then
            xnuc(2,i)=-y
            xnuc(3,i)=z
          end if
          go to 500
        end if
      end do
  500   continue
      end do
      end
c=============================================================================
      subroutine rotmol(na,xnuc,iaxis,angle)
      implicit real*8 (a-h,o-z)
      dimension xnuc(3,na)
c rotates a molecule around an axis (x=1,y=2,z=3) by angle counterclockwise
      cc=cos(angle)
      ss=sin(angle)
c      print *, 'cc,ss',cc,ss
      do i=1,na
        if(iaxis.eq.1) then
          yy=cc*xnuc(2,i)-ss*xnuc(3,i)
          xnuc(3,i)=ss*xnuc(2,i)+cc*xnuc(3,i)
          xnuc(2,i)=yy
        else if (iaxis.eq.3) then
c          print *,i,xnuc(1,i),xnuc(2,i),xnuc(3,i)
          xx=cc*xnuc(1,i)-ss*xnuc(2,i)
          xnuc(2,i)=ss*xnuc(1,i)+cc*xnuc(2,i)
          xnuc(1,i)=xx
c          print *,i,xnuc(1,i),xnuc(2,i),xnuc(3,i)
        else if (iaxis.eq.2) then
          zz=cc*xnuc(3,i)-ss*xnuc(1,i)
          xnuc(1,i)=ss*xnuc(3,i)+cc*xnuc(1,i)
          xnuc(3,i)=zz
        end if
      end do
      end
c========================================================================================
c      subroutine printnuc(na,q,xnuc)
c      implicit real*8 (a-h,o-z)
c      dimension q(na),xnuc(3,na)
c      do i=1,na
c      print 100, q(i),(xnuc(l,i),l=1,3)
c      end do
c  100 format('N=',8x,f10.2,3f10.6)
c      end
c===========================================================================
      subroutine center_of_mass(na,q,xnuc)
      implicit real*8 (a-h,o-z)
      dimension xnuc(3,na),q(na),com(3)
      do l=1,3
      com(l)=0.0d0
      end do
      qq=0.0d0
      do i=1,na
      do l=1,3
        com(l)=com(l)+abs(q(i))*xnuc(l,i)
      end do
      qq=qq+abs(q(i))
      end do
      do l=1,3
      com(l)=com(l)/qq
      end do
      do i=1,na
      do l=1,3
        xnuc(l,i)=xnuc(l,i)-com(l)
      end do
      end do
      end
c==========================================================================
      subroutine symplane(na,q,xnuc,tol,iwhich,plane,found)
      implicit real*8 (a-h,o-z)
c this program will find symmetry planes in a molecule
c parameters:
c  INPUT:
c     na=number of nuclei
c     q(na): nuclear charges
c     xnuc(3,na)=coordinates (already shifted to the center of charge
c     tol: tolerance for approximate symmetry, a dimensionless number (relative)
c     iwhich= integer, which plane do we want (1 for the first, 2 for the second)
c  OUTPUT:
c     plane: the cooordinates of the unit vector normal to the plane
c     found= logical, TRUE if a plane was found
c
      logical found,testplane
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,
     1 pi=3.14159 26535 89793d0)
      dimension xnuc(3,na),q(na),plane(3),x0(3),x1(3),xx(3)
      do i=1,na
      r0=zero
      do l=1,3
        x0(l)=xnuc(l,i)
        r0=r0+x0(l)**2
      end do
c        print *, 'i,xi,ri ',i,(xnuc(l,i),l=1,3), r0
      do 1000 j=1,i-1
c cycle through atom pairs
        if(abs(q(i)-q(j)).gt.1.0d-8) go to 1000
        r1=zero
        do l=1,3
          x1(l)=xnuc(l,j)
          r1=r1+x1(l)**2
        end do
        rmax=max(r0,r1)
        rmax=sqrt(rmax)
c          print *, 'atoms & max. dist. ', i,j,rmax
c  exit if the two atoms are at different distances
        if(abs(r0-r1).gt.rmax*tol) go to 1000
c find their difference
        do l=1,3
          plane(l)=x0(l)-x1(l)
        end do
c normalized difference vector i-j (from j to i) in plane
        call nom1(plane)
c test if all atoms obey the reflection symmetry to this plane
          found=testplane(na,q,xnuc,tol,plane)
          if(.not.found) go to 1000
c      print *,'plane found',i,j,plane
          if(iwhich.eq.1) then
            return
          end if
          if(iwhich.eq.2) then
c         offaxis=sqrt(plane(1)**2+plane(2)**2)
            offaxis=abs(plane(3))
c try to find a symmetry plane perpendicular to z
            if(offaxis.le.tol) then
              return
            else
             found=.false.
            end if
          end if
 1000   continue
c j loop ends
      end do
c  i loop ends
      found=.false.
      end
C=================================================================================
      real*8 function distance(u,v)
      implicit real*8 (a-h,o-z)
      dimension u(3),v(3)
      s=0.0d0
      do l=1,3
        s=s+(u(l)-v(l))**2
      end do
      distance=sqrt(s)
      end
c===========================================================================
      logical function testplane(na,q,xnuc,tol,plane)
      implicit real*8 (a-h,o-z)
c This function returns .TRUE. if the plane is indeed a
c symmetry plane
c Arguments:
c INPUT:
c  na=number of atoms
c  xnuc(3,na)=Cartesian coordinates
c  q(na)=nuclear charges (double)
c  plane(3): vector normal to the plane
c
      dimension xnuc(3,na),q(na),plane(3),xx(3)
      do k=1,na
        rr=0.0d0
        do l=1,3
          xx(l)=xnuc(l,k)
          rr=rr+xx(l)**2
        end do
        rr=sqrt(rr)
        tol1=tol*rr
        if(tol.gt.tol1) tol1=tol
        ss=scalar(xx,plane)
        do l=1,3
          xx(l)=xx(l)-2.0d0*ss*plane(l)
        end do
c  xx is the reflected image vector of k
c      print *,'reflected image of ',k,xx
        qk=q(k)
        do m=1,na
          qm=q(m)
c compare only identical atoms
          if(abs(qk-qm).lt.1.0d-8) then
            dev=distance(xx,xnuc(1,m))
c           print *,'dist. between ',k,m,'is ',dev
            if(dev.lt.tol1) go to 500
          end if
        end do
c m loop ends
c  no atom pair was found to atom k under this operation
c exits this iteration of the j loop
       go to 1000
  500     continue
c                   print *,i,j,k,m
c  if it got so far, then atom k has a symmetry partner under this plane (ij)
c atom m is a symmetry pair of atom k
      end do
c k loop end
      testplane=.true.
      return
 1000 continue
      testplane=.false.
      end
