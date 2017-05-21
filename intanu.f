c----------------------------------------------------------
c  This differ from source/intanu.f only in the subroutine
c  INTAMX which is used for the NMR/CPHF method.
c     Here the intamx rout. calculates :
c
c     <mi|ni(x-1)>,<mi|ni(x+1)>
c     <mi|ni(y-1)>,<mi|ni(y+1)>  in order to get <mi|dni/dx,y,z>
c     <mi|ni(z-1)>,<mi|ni(z+1)>        in the inc123.f
c
c    instead of
c
c    <mi(x-1) | ni>, <mi(x+1) | ni>
c    <mi(y-1) | ni>, <mi(y+1) | ni>
c    <mi(z-1) | ni>, <mi(z+1) | ni>
c
c    and then the <dmi/dx,y,z|ni> was formed in the inc123
c=============================================================
c  there are the following subroutines :
c
c        intash
c        intamp
c        intamx
c        nuch01
c        fcfh01
c        frmh01
c=============================================================
      subroutine intash(ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,l3,xint,kl)
c  see comments in intar
      implicit real*8 (a-h,o-z)
      dimension xa(3),xb(3),xint(900), xfc(63),yfc(63),zfc(63), idg(8)
      dimension nfu(9),mulx(88),muly(88),mulz(88), la(3), lb(3), l3(3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension len(8)
      dimension l1r(3),l2r(3)
      data len/1,3,4,6,10,15,21,28/
      data idg/0,1,1,2,3,4,5,6/
      data nfu/0,1,4,8,14,24,39,60,88/
      data mulx/0, 1,0,0, 0,1,0,0, 2,0,0,1,1,0, 3,2,2,1,1,1,0,0,0,0,
     1          4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, 5,4,4,3,3,3,2,2,2,2,
     2          1,1,1,1,1,0,0,0,0,0,0,  6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,
     3          1,1,1,1,1,1,0,0,0,0,0,0,0/
      data muly/0, 0,1,0, 0,0,1,0, 0,2,0,1,0,1, 0,1,0,2,1,0,3,2,1,0,
     1          0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, 0,1,0,2,1,0,3,2,1,0,
     2          4,3,2,1,0,5,4,3,2,1,0,  0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     3          5,4,3,2,1,0,6,5,4,3,2,1,0/
      data mulz/0, 0,0,1, 0,0,0,1, 0,0,2,0,1,1, 0,0,1,0,1,2,0,1,2,3,
     1          0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, 0,0,1,0,1,2,0,1,2,3,
     2          0,1,2,3,4,0,1,2,3,4,5,  0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     3          0,1,2,3,4,5,0,1,2,3,4,5,6/
      data dfac/2.309401076758d0/
c
c this is 4/sqrt(3)
c
      twoa=two*a
      ln=len(ityp)*len(jtyp)
      iadg=idg(ityp)
      ibdg=idg(jtyp)
      do 300 k1=1,4
      loc=(k1-1)*ln
         do 30 i=1,3
         la(i)=0
      l1r(i)=iadg+la(i)+1
  30  l2r(i)=ibdg+lb(i)+1
 40   continue
      if(k1.eq.4) go to 41
         l1r(k1)=l1r(k1)+1
 41   continue
         call xtor(l1r(1)-1,l2r(1)-1,l3(1),a,b,xa(1),xb(1),xfc)
         call xtor(l1r(2)-1,l2r(2)-1,l3(2),a,b,xa(2),xb(2),yfc)
         call xtor(l1r(3)-1,l2r(3)-1,l3(3),a,b,xa(3),xb(3),zfc)
      ij=0
      ifu=nfu(ityp)+1
      ifu2=nfu(ityp+1)
      jfu=nfu(jtyp)+1
      jfu2=nfu(jtyp+1)
      do 20 i=ifu,ifu2
         if(k1.ne.4)  la(k1)=+1
         ixf=mulx(i)+la(1)+1
         ixff=(ixf-1)*l2r(1)
         iyf=muly(i)+la(2)+1
         iyff=(iyf-1)*l2r(2)
         izf=mulz(i)+la(3)+1
         izff=(izf-1)*l2r(3)
      do 20 j=jfu,jfu2
         jxf=mulx(j)+lb(1)+1
         jyf=muly(j)+lb(2)+1
         jzf=mulz(j)+lb(3)+1
         ij=ij+1
         if (ixf.le.0.or.iyf.le.0.or.izf.le.0.or.jxf.le.0.or.jyf.le.0.or
     1   .jzf.le.0) go to 10
         xint(ij+loc)=xfc(ixff+jxf)*yfc(iyff+jyf)*zfc(izff+jzf)*s0
         go to 20
   10    xint(ij+loc)=zero
   20   continue
  300   continue
      return
c
      end
c=============================================================
      subroutine intamp(ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,l3,xint,kl,
     *ll6)
      implicit real*8 (a-h,o-z)
      dimension xa(3),xb(3), xfc(63),yfc(63),zfc(63), idg(8)
      dimension nfu(9),mulx(88),muly(88),mulz(88), la(3), lb(3), l3(3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension len(8)
      dimension xint(ll6)
      dimension l1r(3),l2r(3)
      dimension lxyz(6,3)
      data lxyz/ 0, 0,-1,+1,+1,-1,
     *          +1,-1, 0, 0,-1,+1,
     *          -1,+1,+1,-1, 0, 0/
      data len/1,3,4,6,10,15,21,28/
      data idg/0,1,1,2,3,4,5,6/
      data nfu/0,1,4,8,14,24,39,60,88/
      data mulx/0, 1,0,0, 0,1,0,0, 2,0,0,1,1,0, 3,2,2,1,1,1,0,0,0,0,
     1          4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, 5,4,4,3,3,3,2,2,2,2,
     2          1,1,1,1,1,0,0,0,0,0,0,  6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,
     3          1,1,1,1,1,1,0,0,0,0,0,0,0/
      data muly/0, 0,1,0, 0,0,1,0, 0,2,0,1,0,1, 0,1,0,2,1,0,3,2,1,0,
     1          0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, 0,1,0,2,1,0,3,2,1,0,
     2          4,3,2,1,0,5,4,3,2,1,0,  0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     3          5,4,3,2,1,0,6,5,4,3,2,1,0/
      data mulz/0, 0,0,1, 0,0,0,1, 0,0,2,0,1,1, 0,0,1,0,1,2,0,1,2,3,
     1          0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, 0,0,1,0,1,2,0,1,2,3,
     2          0,1,2,3,4,0,1,2,3,4,5,  0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     3          0,1,2,3,4,5,0,1,2,3,4,5,6/
      data dfac/2.309401076758d0/
c
c this is 4/sqrt(3)
c
       twoa=two*a
       ln=len(ityp)*len(jtyp)
       iadg=idg(ityp)
       ibdg=idg(jtyp)
                       do 300 k1=1,6
                       loc=(k1-1)*ln
         do 30 i=1,3
         la(i)=0
         lb(i)=lxyz(k1,i)
         l1r(i)=iadg+la(i)+1
  30     l2r(i)=ibdg+lb(i)+1
         call xtor(l1r(1)-1,l2r(1)-1,l3(1),a,b,xa(1),xb(1),xfc)
         call xtor(l1r(2)-1,l2r(2)-1,l3(2),a,b,xa(2),xb(2),yfc)
         call xtor(l1r(3)-1,l2r(3)-1,l3(3),a,b,xa(3),xb(3),zfc)
      ij=0
      ifu=nfu(ityp)+1
      ifu2=nfu(ityp+1)
      jfu=nfu(jtyp)+1
      jfu2=nfu(jtyp+1)
      do 20 i=ifu,ifu2
         ixf=mulx(i)+la(1)+1
         ixff=(ixf-1)*l2r(1)
         iyf=muly(i)+la(2)+1
         iyff=(iyf-1)*l2r(2)
         izf=mulz(i)+la(3)+1
         izff=(izf-1)*l2r(3)
      do 20 j=jfu,jfu2
         jxf=mulx(j)+lb(1)+1
         jyf=muly(j)+lb(2)+1
         jzf=mulz(j)+lb(3)+1
         ij=ij+1
         if (ixf.le.0.or.iyf.le.0.or.izf.le.0.or.jxf.le.0.or.jyf.le.0.or
     1   .jzf.le.0) go to 10
         xint(ij+loc)=xfc(ixff+jxf)*yfc(iyff+jyf)*zfc(izff+jzf)*s0
         go to 20
   10    xint(ij+loc)=zero
   20   continue
  300   continue
      return
      end
c=============================================================
      subroutine intamx(ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,l3,xint,kl,
     *ll6)
      implicit real*8 (a-h,o-z)
      dimension xa(3),xb(3), xfc(63),yfc(63),zfc(63), idg(8)
      dimension nfu(9),mulx(88),muly(88),mulz(88), la(3), lb(3), l3(3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension len(8)
      dimension xint(ll6)
      dimension l1r(3),l2r(3)
      dimension lxyz(6,3)
      data lxyz/-1,+1, 0, 0, 0, 0,
     *           0, 0,-1,+1, 0, 0,
     *           0, 0, 0, 0,-1,+1/
      data len/1,3,4,6,10,15,21,28/
      data idg/0,1,1,2,3,4,5,6/
      data nfu/0,1,4,8,14,24,39,60,88/
      data mulx/0, 1,0,0, 0,1,0,0, 2,0,0,1,1,0, 3,2,2,1,1,1,0,0,0,0,
     1          4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, 5,4,4,3,3,3,2,2,2,2,
     2          1,1,1,1,1,0,0,0,0,0,0,  6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,
     3          1,1,1,1,1,1,0,0,0,0,0,0,0/
      data muly/0, 0,1,0, 0,0,1,0, 0,2,0,1,0,1, 0,1,0,2,1,0,3,2,1,0,
     1          0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, 0,1,0,2,1,0,3,2,1,0,
     2          4,3,2,1,0,5,4,3,2,1,0,  0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     3          5,4,3,2,1,0,6,5,4,3,2,1,0/
      data mulz/0, 0,0,1, 0,0,0,1, 0,0,2,0,1,1, 0,0,1,0,1,2,0,1,2,3,
     1          0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, 0,0,1,0,1,2,0,1,2,3,
     2          0,1,2,3,4,0,1,2,3,4,5,  0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     3          0,1,2,3,4,5,0,1,2,3,4,5,6/
      data dfac/2.309401076758d0/
c
c this is 4/sqrt(3)
c
c      twoa=two*a
       ln=len(ityp)*len(jtyp)
       iadg=idg(ityp)
       ibdg=idg(jtyp)
                       do 300 k1=1,6
                       loc=(k1-1)*ln
         do 30 i=1,3
         la(i)=0
         lb(i)=lxyz(k1,i)
         l1r(i)=iadg+la(i)+1
  30     l2r(i)=ibdg+lb(i)+1
         call xtor(l1r(1)-1,l2r(1)-1,l3(1),a,b,xa(1),xb(1),xfc)
         call xtor(l1r(2)-1,l2r(2)-1,l3(2),a,b,xa(2),xb(2),yfc)
         call xtor(l1r(3)-1,l2r(3)-1,l3(3),a,b,xa(3),xb(3),zfc)
      ij=0
      ifu=nfu(ityp)+1
      ifu2=nfu(ityp+1)
      jfu=nfu(jtyp)+1
      jfu2=nfu(jtyp+1)
      do 20 i=ifu,ifu2
         ixf=mulx(i)+la(1)+1
         ixff=(ixf-1)*l2r(1)
         iyf=muly(i)+la(2)+1
         iyff=(iyf-1)*l2r(2)
         izf=mulz(i)+la(3)+1
         izff=(izf-1)*l2r(3)
      do 20 j=jfu,jfu2
         jxf=mulx(j)+lb(1)+1
         jyf=muly(j)+lb(2)+1
         jzf=mulz(j)+lb(3)+1
         ij=ij+1
         if (ixf.le.0.or.iyf.le.0.or.izf.le.0.or.jxf.le.0.or.jyf.le.0.or
     1   .jzf.le.0) go to 10
         xint(ij+loc)=xfc(ixff+jxf)*yfc(iyff+jyf)*zfc(izff+jzf)*s0
         go to 20
   10    xint(ij+loc)=zero
   20   continue
  300   continue
      return
      end
c=============================================================
      subroutine nuch01 (n1,n2,a,b,sab,ax,ay,az,bx,by,bz,cx,cy,cz,la,lb,
     1                   kk,xint,ln)
c-----------------------------------------------------------------------
c
c     this is for total 2 pluses
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
      dimension ndeg(8), la(3), lb(3), xint(225)
      data ndeg/0,1,1,2,3,4,5,6/
c-----------------------------------------------------------------------
c     calculate number of roots and xrys for rys
c
ccccc nroots=(ndeg(n1)+ndeg(n2)+3)/2
c     nroots=(ndeg(n1)+ndeg(n2)+4)/2
c     if(nroots.eq.0) nroots=1
c-----------------------------------------------------------------------
      p=a+b
      p1=1.d0/p
      px=(a*ax+b*bx)*p1
      py=(a*ay+b*by)*p1
      pz=(a*az+b*bz)*p1
      rpc2=(px-cx)**2+(py-cy)**2+(pz-cz)**2
c-----------------------------------------------------------------------
c     xrys=p*rpc2
c
c     calculate roots
c
c     if(nroots.le.3) call rt123
c     if(nroots.eq.4) call root4
c-----------------------------------------------------------------------
c     calculate two dimentional integrals
c
      call fcfh01(n1,n2,sab,p,px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz,la,lb)
c-----------------------------------------------------------------------
c     assemble core attraction integrals
c
c---> call frmh01(n1,n2,xint)
      call frmca (n1,n2,xint)
c-----------------------------------------------------------------------
c
      end
c=============================================================
      subroutine fcfh01 (n1,n2,sab,p,px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz
     1,la,lb)
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c
c     calculation of the two-dimensional integrals, ix, iy, and iz
c     for the evaluation of the core attraction integrals
c
c-----------------------------------------------------------------------
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /xyzfc/ xfc(343),yfc(343),zfc(343)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      common /hermit/ h(55)
      common /wermit/ w(55)
      dimension itm4(7), ndege(8)
      dimension minh(8), maxh(8)
      dimension la(3),lb(3)
c
c-----------------------------------------------------------------------
      dimension xaa(3,36,6),xbb(3,36,6)
c
      dimension xyza(3,36),xyzb(3,36)
c-----------------------------------------------------------------------
      data itm4/0,7,14,21,28,35,42/
c     this is a relic and simulates a 3-index array xfc(i2,i1,iroot)
      data ndege/1,2,2,3,4,5,6,7/
      data minh/1,2,4,7,11,16,22,29/,maxh/1,3,6,10,15,21,28,36/
c-----------------------------------------------------------------------
c
      ppx=p*px
      ppy=p*py
      ppz=p*pz
      const=two*sqrt(p/pi)*sab
      ie1=ndege(n1)
      ie2=ndege(n2)
      mxh1=(ie1+ie2+2)/2
      mxh1=maxh(mxh1)
c
      indz=0
      do 120 iroot=1,nroots
         u2=p*rysr(iroot)
         rw=rysw(iroot)*const
         ppu2=p+u2
         ppu21=one/ppu2
         sqrt1=sqrt(ppu2)
         sqrt11=one/sqrt1
         ttx=(ppx+u2*cx)*ppu21
         tty=(ppy+u2*cy)*ppu21
         ttz=(ppz+u2*cz)*ppu21
         do 20 i=1,mxh1
            hh1=h(i)*sqrt11
            x1=hh1+ttx
            y1=hh1+tty
            z1=hh1+ttz
ckw99       if (ie1.eq.0) go to 10
c
            xyza(1,i)=x1-ax
            xyza(2,i)=y1-ay
            xyza(3,i)=z1-az
c
             do 200 ki=1,3
                xaa(ki,i,1)=one
                if(la(ki).eq. 2) xaa(ki,i,1)=xyza(ki,i)*xyza(ki,i)
                if(la(ki).eq. 1) xaa(ki,i,1)=xyza(ki,i)
ckw99           if(la(ki).eq. 0) xaa(ki,i,1)=one
                if(la(ki).eq.-1) xaa(ki,i,1)=zero
  200        continue
            if (ie1.eq.1) go to 10
c---
             do 201 ki=1,3
             xaa(ki,i,2)=xyza(ki,i)
             if(la(ki).eq.-1) then
                xaa(ki,i,2)=one
             else
                xaa(ki,i,2)=xyza(ki,i)*xaa(ki,i,1)
             endif
  201        continue
c
ckw99
ckw99       if (ie1.eq.2) go to 10
ckw99
             do iaa=3,ie1
                do ki=1,3
                   xaa(ki,i,iaa)=xyza(ki,i)*xaa(ki,i,iaa-1)
                enddo
             enddo
ckw99
c---
   10       continue
c---
ckw99       if (ie2.eq.0) go to 20
c
               xyzb(1,i)=x1-bx
               xyzb(2,i)=y1-by
               xyzb(3,i)=z1-bz
c
             do 206 ki=1,3
                xbb(ki,i,1)=one
                if(lb(ki).eq. 2) xbb(ki,i,1)=xyzb(ki,i)*xyzb(ki,i)
                if(lb(ki).eq. 1) xbb(ki,i,1)=xyzb(ki,i)
ckw99           if(lb(ki).eq. 0) xbb(ki,i,1)=one
                if(lb(ki).eq.-1) xbb(ki,i,1)=zero
  206        continue
            if (ie2.eq.1) go to 20
c---
             do 207 ki=1,3
             xbb(ki,i,2)=xyzb(ki,i)
             if(lb(ki).eq.-1) then
                xbb(ki,i,2)=one
             else
                xbb(ki,i,2)=xyzb(ki,i)*xbb(ki,i,1)
             endif
  207        continue
ckw99
ckw99       if (ie2.eq.2) go to 20
ckw99
             do ibb=3,ie2
                do ki=1,3
                   xbb(ki,i,ibb)=xyzb(ki,i)*xbb(ki,i,ibb-1)
                enddo
             enddo
c---
   20    continue
c
         do 110 i1=1,ie1
            ind1=itm4(i1)
         do 110 i2=1,ie2
            ind=ind1+i2
            mpts=(i1+i2+2)/2
            minh1=minh(mpts)
            maxh1=maxh(mpts)
c
            fcx=zero
            fcy=zero
            fcz=zero
            do 100 i=minh1,maxh1
ckw99
               f1x=xaa(1,i,i1)*xbb(1,i,i2)
               f1y=xaa(2,i,i1)*xbb(2,i,i2)
               f1z=xaa(3,i,i1)*xbb(3,i,i2)
c
               fcx=fcx+f1x*w(i)
               fcy=fcy+f1y*w(i)
               fcz=fcz+f1z*w(i)
c
  100       continue
            xfc(ind+indz)=fcx
            yfc(ind+indz)=fcy
            zfc(ind+indz)=fcz*rw
  110    continue
c
  120 indz=indz+49
c
      end
