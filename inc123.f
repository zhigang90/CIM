ckwolinski Feb. , 1994
c in the incsh2 subroutine some terms are commented out since
c it has been found that they cancel each other
c    I(mi,ni) - II(mi,ni) and elements of the H11N
c
c=============================================================
c   there are subroutines : incsh1, incsh2, incsh3
c=============================================================
      subroutine incsh1 (ityp,jtyp,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,
     *                   xint,lxin)

      use memory

      implicit real*8 (a-h,o-z)
      dimension xa(3),xb(3),la(3),lb(3),ll(3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
ckwol common /tape/ inp,inp2,iout
      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /gauge/ gauger(3)

      ! dimensions and data up to I-type functions

                             !28*28*4
      dimension xint(lxin),yint(3136)
                             !28*28*6
      dimension ximp(4704)

      dimension ndeg(8)
      data ndeg/0,1,1,2,3,4,5,6/
      dimension len(8)
      data len/1,3,4,6,10,15,21,28/
      dimension  nfu(9)
      data nfu/0,1,4,8,14,24,39,60,88/
      dimension mulx(88), muly(88), mulz(88)
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

      twoa=two*a
      twob=two*b
      ln=len(ityp)*len(jtyp)
      ln4=ln*4
      ia=nfu(ityp)+1
      ie=nfu(ityp+1)
      ja=nfu(jtyp)+1
      je=nfu(jtyp+1)
      do 20 i=1,3
         la(i)=0
         lb(i)=0
         ll(i)=0
   20 continue
      go to(30,40,500) key
c*************
   30 continue
      do 10 i=1,ln4
         xint(i)=zero
   10 continue
c
c this is for an external electric field
c and electric field gradient (k1=1,...9)
      if(k1.ne.0) then
         if(k1.le.3) then
              ll(k1)=1
         else
            if(k1.eq.4) ll(1)=2
            if(k1.eq.6) ll(2)=2
            if(k1.eq.9) ll(3)=2
            if(k1.eq.5) then
               ll(1)=1
               ll(2)=1
            endif
            if(k1.eq.7) then
               ll(1)=1
               ll(3)=1
            endif
            if(k1.eq.8) then
               ll(2)=1
               ll(3)=1
            endif
         endif
      endif
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,xint,k1)
      go to 240
c
  40  continue
      do 41 l=1,ln4
  41  yint(l)=zero
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      ij=0
      do 50 i=ia,ie
      do 50 j=ja,je
         ij=ij+1
         ft=float(2*(mulx(j)+muly(j)+mulz(j))+3)*b
         xint(ij)=yint(ij)*ft
         xint(ij+ln)=yint(ij+ln)*ft
         xint(ij+2*ln)=yint(ij+2*ln)*ft
         xint(ij+3*ln)=yint(ij+3*ln)*ft
   50 continue
      lb(1)=2
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      ft=-two*b**2
      do 60 i=1,ln4
   60 xint(i)=xint(i)+yint(i)*ft
      lb(1)=0
      lb(2)=2
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      do 70 i=1,ln4
   70 xint(i)=xint(i)+yint(i)*ft
      lb(2)=0
      lb(3)=2
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      lb(3)=0
      do 80 i=1,ln4
   80 xint(i)=xint(i)+yint(i)*ft
      if (jtyp.lt.4) go to 240
      lb(3)=0
      lb(1)=-2
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      ij=0
      do 90 i=ia,ie
      do 90 j=ja,je
         ij=ij+1
         ft=float((mulx(j)*(mulx(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
         xint(ij+ln)=xint(ij+ln)-yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
         xint(ij+3*ln)=xint(ij+3*ln)-yint(ij+3*ln)*ft
   90 continue
      lb(1)=0
      lb(2)=-2
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      ij=0
      do 100 i=ia,ie
      do 100 j=ja,je
         ij=ij+1
         ft=float((muly(j)*(muly(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
         xint(ij+ln)=xint(ij+ln)-yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
         xint(ij+3*ln)=xint(ij+3*ln)-yint(ij+3*ln)*ft
  100 continue
      lb(2)=0
      lb(3)=-2
      call intash (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,yint,k1)
      lb(3)=0
      ij=0
      do 110 i=ia,ie
      do 110 j=ja,je
         ij=ij+1
         ft=float((mulz(j)*(mulz(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
         xint(ij+ln)=xint(ij+ln)-yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
         xint(ij+3*ln)=xint(ij+3*ln)-yint(ij+3*ln)*ft
  110 continue
      go to 240
c*************
 500  continue
      ll3=ln*3
      ll6=ln*6
      do 502 l=1,ll6
 502  ximp(l)=zero
      call intamp (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,ximp,k1,ll6)
      ij=0
      do 31 i=ia,ie
      do 31 j=ja,je
      fxj=float(mulx(j))
      fyj=float(muly(j))
      fzj=float(mulz(j))
      ij=ij+1
      xint(ij)=fzj*ximp(ij)-fyj*ximp(ij+ln)
      xint(ij+ln)=fxj*ximp(ij+2*ln)-fzj*ximp(ij+3*ln)
      xint(ij+2*ln)=fyj*ximp(ij+4*ln)-fxj*ximp(ij+5*ln)
  31  continue
      if(nogiao.ne.0) then
         do 503 l=1,ll6
 503     ximp(l)=zero
c
       call intamx (ityp,jtyp,a,b,sqa,sqb,s0,xa,xb,la,lb,ll,ximp,k1,ll6)
c
         xbg=xb(1)-gauger(1)
         ybg=xb(2)-gauger(2)
         zbg=xb(3)-gauger(3)
c
         ij=0
         do 32 i=ia,ie
         do 32 j=ja,je
         fxj=float(mulx(j))
         fyj=float(muly(j))
         fzj=float(mulz(j))
         ij=ij+1
c
         adx= fxj*ximp(ij)     -twob*ximp(ij+ln)
         ady= fyj*ximp(ij+2*ln)-twob*ximp(ij+3*ln)
         adz= fzj*ximp(ij+4*ln)-twob*ximp(ij+5*ln)
c------
c        xint(ij)=xint(ij)          +   xb(2)*adz - xb(3)*ady
c        xint(ij+ln)=xint(ij+ln)    +   xb(3)*adx - xb(1)*adz
c        xint(ij+2*ln)=xint(ij+2*ln)+   xb(1)*ady - xb(2)*adx
c-----
         xint(ij)=xint(ij)          +   ybg*adz - zbg*ady
         xint(ij+ln)=xint(ij+ln)    +   zbg*adx - xbg*adz
         xint(ij+2*ln)=xint(ij+2*ln)+   xbg*ady - ybg*adx
  32     continue
      endif
c*************
  240 continue
      return
      end
c=============================================================
      subroutine incsh2(ityp,jtyp,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,
     *                  xint,   xra   ,lxin,
     *                  rp01,rp10,rpmm,rpmp,rppm,
     *                  rpp0,rm01,rm10,xin0, tab_mm,tab_pmm,tab_mp,
     *                  tab_pmp,tab_ppm )
c------------------------------------------------------------------
c This routine calculates three types of 1el integrals over the GIAO
c basis set :
c  (1)          (mi1|h01n|ni0) + (mi0|h01n|ni1)
c                              +
c  (2)                  (mi0|h11n|ni0)
c
c These two are needed for the diamagnetic part of the shielding tensor
c
c and
c                     (mi0 | h01n | ni0)
c
c needed for the paramagnetic part of the shielding tensor
c-----------------------------------------------------------------------
c This routine was changed on February 28, 1999 by Tomasz Janowski
c in order to
c (1) eliminate fixed dimension arrays
c (2) speed up calculations by reducing number of calls of nuch01
c     (instead, the shifting of ang.mom. has been used where possible).
c-----------------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
      dimension xa(3),xb(3),la(3),lb(3),xra(3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /tape/ inp,inp2,iout,it  ,jt,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /gauge/ gauger(3)
c
      dimension ikb(9),ikc(9),iks(9),ikt(9)
      dimension kai(9),kbi(9)
      data ikb/3,3,3,1,1,1,2,2,2/
      data ikc/2,2,2,3,3,3,1,1,1/
      data iks/2,3,1,2,3,1,2,3,1/
      data ikt/3,1,2,3,1,2,3,1,2/
      data kai/1,1,1,2,2,2,3,3,3/
      data kbi/1,2,3,1,2,3,1,2,3/
      dimension xint(lxin)
c----------------------------------------------------------------------
c maximum dimension needed here is 27*leni*lenj .
c For ii it is 27*28*28= 21168
c----------------------------------------------------------------------
      dimension rp01(*),rp10(*),rpmm(*),rpmp(*),rppm(*),
     *          rpp0(*),rm01(*),rm10(*),xin0(*)
c----------------------------------------------------------------------
c New arrays for temporary storage of some intermediate integrals
c
      dimension tab_mm(*),tab_pmm(*),tab_mp(*),tab_pmp(*),tab_ppm(*)
c----------------------------------------------------------------------
      dimension idk(3,3)
      data idk/1,0,0,
     1         0,1,0,
     2         0,0,1/

      ! dimensions and data up to I-type functions

      dimension ndeg(8)
      data ndeg/0,1,1,2,3,4,5,6/
      dimension len(8)
      data len/1,3,4,6,10,15,21,28/
      dimension  nfu(9)
      data nfu/0,1,4,8,14,24,39,60,88/

      dimension mulx(88), muly(88), mulz(88)
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

      twoa=two*a
      twob=two*b
      ln=len(ityp)*len(jtyp)
      ln4=ln*4
      ia=nfu(ityp)+1
      ie=nfu(ityp+1)
      ja=nfu(jtyp)+1
      je=nfu(jtyp+1)
c---------------------------------------------------------------------
c  for given pair of primitive shells and for given reference atom
c     the number of roots , xrys and roots  for rys polinominal
c                     are calculated
c---------------------------------------------------------------------
c
      nroots=(ndeg(ityp)+ndeg(jtyp)+4)/2
c
      ax=xa(1)
      ay=xa(2)
      az=xa(3)
c
      bx=xb(1)
      by=xb(2)
      bz=xb(3)
c
      cx=xra(1)
      cy=xra(2)
      cz=xra(3)
c
      p=a+b
      p1=one/p
      px=(a*ax+b*bx)*p1
      py=(a*ay+b*by)*p1
      pz=(a*az+b*bz)*p1
      rpc2=(px-cx)**2+(py-cy)**2+(pz-cz)**2
      xrys=p*rpc2
c
c     calculate roots
c
      if (nroots.le.3) call rt123
      if (nroots.eq.4) call root4
      if (nroots.eq.5) call root5
      if (nroots.gt.5) call root6

c*********************************************************************
c   {mi(c-1)|r1n**-1|ni(s-1)} integrals : c,s=x,y,z
c
c               in tab_mm(9*ln)
c*****************                             ***********************
      do 20 i=1,3
         la(i)=0
         lb(i)=0
 20   continue
      if(ityp.eq.1.or.jtyp.eq.1) then
         do 510 l=1,9*ln
 510     tab_mm(l)=zero
      else
         do 500 ic=1,3
            do 500 is=1,3
               ics=(ic-1)*3+is
               loc=(ics-1)*ln
               la(ic)=-1
               lb(is)=-1
               call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),
     *                     xb(2),xb(3),xra(1),xra(2),xra(3),la,lb,k1,
     *                     tab_mm(1+loc),ln)
               la(ic)=0
               lb(is)=0
 500     continue
      endif
c
c-------------------------------------------------------------------
c       {mi(t-1)|r1n**-1|ni(s+1)}  t,s=x,yz
c                 in tab_mp(9*ln)
c
c---------------                             -----------------------
c
      if(ityp.eq.1) then
         do 403 l=1,9*ln
 403     tab_mp(l)=zero
      else
         do 401 i=1,3
           do 401 j=1,3
c*
              ij=(i-1)*3 + j
              loc=(ij-1)*ln
              la(i)=-1
              lb(j)=+1
              call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),
     *                    xb(2),xb(3),xra(1),xra(2),xra(3),la,lb,k1,
     *                    tab_mp(1+loc),ln)
              la(i)=0
              lb(j)=0
 401     continue
      endif
c
c-------------------------------------------------------------------
c       {mi(t+1)|r1n**-1|ni(s-1)}  t,s=x,yz
c                 in tab_mp(9*ln-18*ln)
c
c---------------                             -----------------------
c
      ln9=9*ln
      ln91=ln9+1
      if(jtyp.eq.1) then
         do 404 l=ln91,18*ln
 404     tab_mp(l)=zero
      else
         do 400 i=1,3
            do 400 j=1,3
               loc=ln9 + ( (i-1)*3 + (j-1) )*ln
               la(i)=+1
               lb(j)=-1
               call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),
     *                     xb(2),xb(3),xra(1),xra(2),xra(3),la,lb,k1,
     *                     tab_mp(1+loc),ln)
               la(i)=0
               lb(j)=0
c*
 400     continue
      endif
c*
c----------------------------------------------------------------------
c   begining   for   i(a,r)(mi,ni) *
c----------------------------------------------------------------------
c   {mi0|r1n**-1|ni(t-1)} integrals ;t=x,y,z  /3*ln integrals/
c
c               in rm01(3*ln)     /1,2/
c****************                                   *******************
      do 24 i=1,3
         la(i)=0
         lb(i)=0
 24   continue
c
      if(jtyp.eq.1) then
         do 101 l=1,3*ln
            rm01(l)=zero
 101     continue
      else
         do 100 i=1,3
            loc=(i-1)*ln
            lb(i)=-1
            call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),
     *                  xb(3),xra(1),xra(2),xra(3),la,lb,k1,rm01(1+loc),
     *                  ln)
            lb(i)=0
 100     continue
      endif
c*
c**********************************************************************
c
c   {mi0|r1n**-1|ni(t+1)} integrals ;t=x,y,z  /3*ln integrals/
c
c               in rp01(3*ln)      /3,4/
c***********************                   ****************************
      do 110 i=1,3
         loc=(i-1)*ln
         lb(i)=+1
         call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),
     *               xb(3),xra(1),xra(2),xra(3),la,lb,k1,rp01(1+loc),ln)
         lb(i)=0
 110  continue
c*
c**********************************************************************
c
c   {mi(c+1,s-1)|r1n**-1|ni(t-1)} integrals : c,s,t=x,y,z ,t.ne.s
c
c               in tab_pmm(27*ln)      /5,6,9,10/
c*******************                            ***********************
      if(jtyp.eq.1 .or. ityp.eq.1) then
         do 121 l=1,27*ln
            tab_pmm(l)=zero
 121     continue
      else
         do 120 ic=1,3
            do 120 is=1,3
               do 120 it=1,3
                  if(it.eq.is) go to 120
                  icst= (ic-1)*9 + (is-1)*3 + it
                  loc=(icst-1)*ln
                  la(ic)=+1
                  la(is)=-1
                  lb(it)=-1
                  if(ic.eq.is) then
                     do 122 l=1,ln
                        tab_pmm(l+loc)=rm01( (it-1)*ln+l)
 122                 continue
                  else
                     call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),
     *                           xb(1),xb(2),xb(3),xra(1),xra(2),xra(3),
     *                           la,lb,k1,tab_pmm(1+loc),ln)
                  endif
                  la(ic)=0
                  la(is)=0
                  lb(it)=0
 120     continue
      endif
c*********************************************************************
c
c   {mi(c+1,s-1)|r1n**-1|ni(t+1)} integrals : c,s,t=x,y,z ,t.ne.s
c
c               in tab_pmp(27*ln)      /7,8,11,12/
c*****************                             ***********************
      if (ityp .eq. 0) then
         do 222 l=1, 27*ln
            tab_pmp(l)=0.0d0
 222     continue
      else
         do 130 ic=1,3
            do 130 is=1,3
               do 130 it=1,3
                  if(it.eq.is) go to 130
                  icst= (ic-1)*9 + (is-1)*3 + it
                  loc=(icst-1)*ln
                  la(ic)=+1
                  la(is)=-1
                  lb(it)=+1
                  if(ic.eq.is) then
                     do 131 l=1,ln
                        tab_pmp(l+loc)=rp01( (it-1)*ln+l )
 131                 continue
                  else
                     call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),
     *                           xb(1),xb(2),xb(3),xra(1),xra(2),
     *                           xra(3),la,lb,k1,tab_pmp(1+loc),ln)
                  endif
                  la(ic)=0
                  la(is)=0
                  lb(it)=0
 130     continue
      endif
c*********************************************************************
c
c   {mi(c+1,s+1)|r1n**-1|ni(t-1)} integrals : c,s,t=x,y,z ,t.ne.s
c
c               in tab_ppm(27*ln)     /13,14,15,16/
c*****************                             ***********************
      if(jtyp.eq.1) then
         do 141 l=1,27*ln
            tab_ppm(l)=zero
 141     continue
      else
         do 140 ic=1,3
            do 140 is=1,3
               do 140 it=1,3
                  if(it.eq.is) go to 140
                  icst= (ic-1)*9 + (is-1)*3 + it
                  loc=(icst-1)*ln
                  la(ic)=+1
                  la(is)=+1
                  if(ic.eq.is) la(ic)=2
                  lb(it)=-1
                  call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),
     *                        xb(2),xb(3),xra(1),xra(2),xra(3),la,lb,k1,
     *                        tab_ppm(1+loc),ln)
                  la(ic)=0
                  la(is)=0
                  lb(it)=0
 140     continue
      endif
c*********************************************************************
c
c   {mi(c+1,s+1)|r1n**-1|ni0)} integrals : c,s=x,y,z ,c.ge.s
c
c               in rpp0(6*ln)       /17,18,19,20/
c*****************                             ***********************
      do 150 ic=1,3
         ii=ic*(ic-1)/2
         do 150 is=1,ic
            ics=ii+is
            loc=(ics-1)*ln
            la(ic)=+1
            la(is)=+1
            if(ic.eq.is) la(ic)=2
            call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),
     *                  xb(3),xra(1),xra(2),xra(3),la,lb,k1,rpp0(1+loc),
     *                  ln)
            la(ic)=0
            la(is)=0
 150  continue
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do 750 icase=1,9
         lcas=(icase-1)*ln
         kb=ikb(icase)
         kc=ikc(icase)
         ks=iks(icase)
         kt=ikt(icase)
c
         bmi=xa(kb)
         cmi=xa(kc)
c
         bmig=bmi-gauger(kb)
         cmig=cmi-gauger(kc)
c
         smi=xa(ks)
         sni=xb(ks)
         tmi=xa(kt)
         tni=xb(kt)
c
         fcs=float( idk(kc,ks) )
         fbs=float( idk(kb,ks) )
         fct=float( idk(kc,kt) )
         fbt=float( idk(kb,kt) )
c
         kt1=kt-1
         ks1=ks-1
         kc1=kc-1
         kb1=kb-1
c
         loc1=kt1*ln
         loc2=ks1*ln
         loc3=(kc1*9 + ks1*3 + kt1)*ln
         loc4=(kc1*9 + kt1*3 + ks1)*ln
         loc5=(kb1*9 + ks1*3 + kt1)*ln
         loc6=(kb1*9 + kt1*3 + ks1)*ln
         loc7=(kc1*kc/2+ks1)*ln
         loc8=(kc1*kc/2+kt1)*ln
         loc9=(kb1*kb/2+ks1)*ln
         loc0=(kb1*kb/2+kt1)*ln
         if(ks.gt.kc) loc7=(ks1*ks/2+kc1)*ln
         if(kt.gt.kc) loc8=(kt1*kt/2+kc1)*ln
         if(ks.gt.kb) loc9=(ks1*ks/2+kb1)*ln
         if(kt.gt.kb) loc0=(kt1*kt/2+kb1)*ln
         ij=0
         do 71 i=ia,ie
            if(ks.eq.1)  fmis=float(mulx(i))
            if(ks.eq.2)  fmis=float(muly(i))
            if(ks.eq.3)  fmis=float(mulz(i))
c
            if(kt.eq.1)  fmit=float(mulx(i))
            if(kt.eq.2)  fmit=float(muly(i))
            if(kt.eq.3)  fmit=float(mulz(i))
            do 71 j=ja,je
               if(ks.eq.1)  fnis=float(mulx(j))
               if(ks.eq.2)  fnis=float(muly(j))
               if(ks.eq.3)  fnis=float(mulz(j))
c
               if(kt.eq.1)  fnit=float(mulx(j))
               if(kt.eq.2)  fnit=float(muly(j))
               if(kt.eq.3)  fnit=float(mulz(j))
               ij=ij+1
c
      xint(ij+lcas)=
c------it is canceled out by corresponding part of H11N-----C
c        only if the GIAO is used
c    1              +(-bmi*fcs + cmi*fbs)*fnit* rm01(loc1+ij)
c    1              -(-bmi*fct + cmi*fbt)*fnis* rm01(loc2+ij)
c    2    -(-bmi*fcs + cmi*fbs )*twob* rp01(loc1+ij)
c    2    +(-bmi*fct + cmi*fbt )*twob* rp01(loc2+ij)
c---
     1              +(-bmig*fcs + cmig*fbs)*fnit* rm01(loc1+ij)
     1              -(-bmig*fct + cmig*fbt)*fnis* rm01(loc2+ij)
     2    -(-bmig*fcs + cmig*fbs )*twob* rp01(loc1+ij)
     2    +(-bmig*fct + cmig*fbt )*twob* rp01(loc2+ij)
c------it is canceled out by corresponding part of H11N-----C
     1    -bmi*fmis*fnit* tab_pmm(loc3+ij)
     1    +bmi*fmit*fnis* tab_pmm(loc4+ij)
     1    +cmi*fmis*fnit* tab_pmm(loc5+ij)
     1    -cmi*fmit*fnis* tab_pmm(loc6+ij)
     2    +bmi*twob*( fmis*tab_pmp(loc3+ij)
     2              - fmit*tab_pmp(loc4+ij) )
     2    -cmi*twob*( fmis*tab_pmp(loc5+ij)
     2               -fmit*tab_pmp(loc6+ij) )
     2    +bmi*twoa*( fnit*tab_ppm(loc3+ij)
     2              - fnis*tab_ppm(loc4+ij) )
     2    -cmi*twoa*( fnit*tab_ppm(loc5+ij)
     2              - fnis*tab_ppm(loc6+ij) )
     3  -twoa*twob*bmi*( (tmi-tni)*rpp0(loc7+ij)
     3                  -(smi-sni)*rpp0(loc8+ij) )
     3  +twoa*twob*cmi*( (tmi-tni)*rpp0(loc9+ij)
     3                  -(smi-sni)*rpp0(loc0+ij) )
c
  71     continue
 750  continue
c**********************************************************************
c   end for  i(a,r)(mi,ni) *
c**********************************************************************
c**********************************************************************
c   begining   for   ii(a,r)(mi,ni) *
c**********************************************************************
      do 34 i=1,3
         la(i)=0
         lb(i)=0
 34   continue
c*
      call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),xb(3),
     *            xra(1),xra(2),xra(3),la,lb,k1,xin0(1),ln)
c
c**********************************************************************
c
c   {mi(t-1)|r1n**-1|ni0)} integrals ;t=x,y,z  /3*ln integrals/
c
c               in rm10(3*ln)     /1,2/
c****************                                   *******************
      if(ityp.eq.1) then
         do 301 l=1,3*ln
            rm10(l)=zero
 301     continue
      else
         do 300 i=1,3
            loc=(i-1)*ln
            la(i)=-1
            call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),
     *                  xb(3),xra(1),xra(2),xra(3),la,lb,k1,rm10(1+loc),
     *                  ln)
            la(i)=0
 300     continue
      endif
c
c**********************************************************************
c
c   {mi(t+1)|r1n**-1|ni0} integrals ;t=x,y,z  /3*ln integrals/
c
c               in rp10(3*ln)      /3,4/
c***********************                   ****************************
c
      do 310 i=1,3
         loc=(i-1)*ln
         do 311 l=1,ln
            rp10(l+loc)=rp01(l+loc)-(xa(i)-xb(i))*xin0(l)
 311     continue
 310  continue
c**********************************************************************
c  CHANGED on February 28, 1999 by Tomasz Janowski.
c
c   {mi(t-1)|r1n**-1|ni(c+1,s-1)} integrals : c,s,t=x,y,z ,t.ne.s
c
c               in rpmm(27*ln)      /5,6,9,10/
c*******************                            ***********************
      if(ityp.eq.1 .or. jtyp.eq.1) then
         do 321 l=1,27*ln
            rpmm(l)=zero
 321     continue
      else
         do 320 ic=1,3
            do 320 is=1,3
               do 320 it=1,3
               if(it.eq.is) go to 320
               ipoz1=(it-1+3*(is-1)+9*(ic-1))*ln
               ipoz2=(is-1+3*(it-1)+9*(ic-1))*ln
               ipoz3=(is-1+3*(it-1))*ln
               do l=1, ln
                  rpmm(ipoz1+l)=tab_pmm(ipoz2+l)+(xa(ic)-xb(ic))*
     *                          tab_mm(ipoz3+l)
               enddo
 320     continue
      endif
c*********************************************************************
c  CHANGED on February 28, 1999 by Tomasz Janowski.
c
c   {mi(t+1)|r1n**-1|ni(c+1,s-1)} integrals : c,s,t=x,y,z ,t.ne.s
c
c               in rpmp(27*ln)      /7,8,11,12/
c*****************                             ***********************
      if (jtyp .eq. 1) then
         do 333 l=1, 27*ln
            rpmp(l)=0.0d0
 333     continue
      else
         do 330 ic=1,3
            do 330 is=1,3
               do 330 it=1,3
                  if(it.eq.is) go to 330
                  ipoz1=((it-1)+3*(is-1)+9*(ic-1))*ln
                  ipoz2=((is-1)+3*(it-1)+9*(ic-1))*ln
                  ipoz3=((is-1)+3*(it-1)+9)*ln
                  do l=1, ln
                     rpmp(ipoz1+l)=tab_ppm(ipoz2+l)+(xa(ic)-xb(ic))*
     *                             tab_mp(ipoz3+l)
                  enddo
 330     continue
      endif
c*********************************************************************
c  CHANGED on February 28, 1999 by Tomasz Janowski.
c
c   {mi(t-1)|r1n**-1|ni(c+1,s+1)} integrals : c,s,t=x,y,z ,t.ne.s
c
c               in rppm(27*ln)     /13,14,15,16/
c*****************                             ***********************
      if(ityp.eq.1) then
         do 341 l=1,27*ln
              rppm(l)=zero
 341     continue
      else
         do 340 ic=1,3
            do 340 is=1,3
               do 340 it=1,3
                  if(it.eq.is) go to 340
                  ipoz1=((it-1)+3*(is-1)+9*(ic-1))*ln
                  ipoz2=((is-1)+3*(it-1)+9*(ic-1))*ln
                  ipoz3=((is-1)+3*(it-1))*ln
                  do l=1, ln
                     rppm(ipoz1+l)=tab_pmp(ipoz2+l)+(xa(ic)-xb(ic))*
     *                             tab_mp(ipoz3+l)
                  enddo
 340     continue
      endif
c*********************************************************************
c
c   {mi0|r1n**-1|ni(c+1,s+1)} integrals : c,s=x,y,z  , c.ge.s
c
c               in rpp0(6*ln)       /17,18,19,20/
c*****************                             ***********************
      do 350 ic=1,3
         ii=ic*(ic-1)/2
         lic=(ic-1)*ln
         do 350 is=1,ic
            lis=(is-1)*ln
            ics=ii+is
            loc=(ics-1)*ln
            do 351 l=1,ln
               rpp0(l+loc)=rpp0(l+loc)+(xa(is)-xb(is))*rp10(lic+l) +
     *                     (xa(ic)-xb(ic))*rp01(l+lis)
 351        continue
 350  continue
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do 950 icase=1,9
c
         lcas=(icase-1)*ln
c
         kb=ikb(icase)
         kc=ikc(icase)
         ks=iks(icase)
         kt=ikt(icase)
c
         bni=xb(kb)
         cni=xb(kc)
c
         bnig=bni-gauger(kb)
         cnig=cni-gauger(kc)
c
         smi=xa(ks)
         sni=xb(ks)
         tmi=xa(kt)
         tni=xb(kt)
c
c
         fcs=float( idk(kc,ks) )
         fbs=float( idk(kb,ks) )
         fct=float( idk(kc,kt) )
         fbt=float( idk(kb,kt) )
c
         kt1=kt-1
         ks1=ks-1
         kc1=kc-1
         kb1=kb-1
c
         loc1=kt1*ln
         loc2=ks1*ln
         loc3=(kc1*9 + ks1*3 + kt1)*ln
         loc4=(kc1*9 + kt1*3 + ks1)*ln
         loc5=(kb1*9 + ks1*3 + kt1)*ln
         loc6=(kb1*9 + kt1*3 + ks1)*ln
         loc7=(kc1*kc/2+ks1)*ln
         loc8=(kc1*kc/2+kt1)*ln
         loc9=(kb1*kb/2+ks1)*ln
         loc0=(kb1*kb/2+kt1)*ln
         if(ks.gt.kc) loc7=(ks1*ks/2+kc1)*ln
         if(kt.gt.kc) loc8=(kt1*kt/2+kc1)*ln
         if(ks.gt.kb) loc9=(ks1*ks/2+kb1)*ln
         if(kt.gt.kb) loc0=(kt1*kt/2+kb1)*ln
c
         ij=0
         do 91 i=ia,ie
            if(ks.eq.1)  fmis=float(mulx(i))
            if(ks.eq.2)  fmis=float(muly(i))
            if(ks.eq.3)  fmis=float(mulz(i))
c
            if(kt.eq.1)  fmit=float(mulx(i))
            if(kt.eq.2)  fmit=float(muly(i))
            if(kt.eq.3)  fmit=float(mulz(i))
            do 91 j=ja,je
               if(ks.eq.1)  fnis=float(mulx(j))
               if(ks.eq.2)  fnis=float(muly(j))
               if(ks.eq.3)  fnis=float(mulz(j))
c
               if(kt.eq.1)  fnit=float(mulx(j))
               if(kt.eq.2)  fnit=float(muly(j))
               if(kt.eq.3)  fnit=float(mulz(j))
               ij=ij+1
      xint(ij+lcas)= xint(ij+lcas)
c------it is canceled out by corresponding part of H11N-----C
c        only if the GIAO is used
c    1             + (-bni*fcs + cni*fbs)*fmit* rm10(loc1+ij)
c    1              -(-bni*fct + cni*fbt)*fmis* rm10(loc2+ij)
c    2    -(-bni*fcs + cni*fbs )*twoa* rp10(loc1+ij)
c    2    +(-bni*fct + cni*fbt )*twoa* rp10(loc2+ij)
c--
     1             + (-bnig*fcs + cnig*fbs)*fmit* rm10(loc1+ij)
     1              -(-bnig*fct + cnig*fbt)*fmis* rm10(loc2+ij)
     2    -(-bnig*fcs + cnig*fbs )*twoa* rp10(loc1+ij)
     2    +(-bnig*fct + cnig*fbt )*twoa* rp10(loc2+ij)
c------it is canceled out by corresponding part of H11N-----C
     1    -bni*fnis*fmit* rpmm(loc3+ij)
     1    +bni*fnit*fmis* rpmm(loc4+ij)
     1    +cni*fnis*fmit* rpmm(loc5+ij)
     1    -cni*fnit*fmis* rpmm(loc6+ij)
      xint(ij+lcas)= xint(ij+lcas)
     2    +bni*twoa*( fnis*rpmp(loc3+ij)
     2              - fnit*rpmp(loc4+ij))
     2    -cni*twoa*( fnis*rpmp(loc5+ij)
     2               -fnit*rpmp(loc6+ij))
     2    +bni*twob*( fmit*rppm(loc3+ij)
     2              - fmis*rppm(loc4+ij) )
     2    -cni*twob*( fmit*rppm(loc5+ij)
     2              - fmis*rppm(loc6+ij) )
     3  +twoa*twob*bni*( (tmi-tni)*rpp0(loc7+ij)
     3                  -(smi-sni)*rpp0(loc8+ij) )
     3  -twoa*twob*cni*( (tmi-tni)*rpp0(loc9+ij)
     3                  -(smi-sni)*rpp0(loc0+ij) )
c
  91     continue
 950  continue
c***************************
c   end for  ii(a,r)(mi,ni) *
c***************************
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c                       (mi0|h11|ni0)
c*************
  600 continue
c
c  the cases  :    1  for  xx   ka=x  kb=x    ami=xa(kb)  ani=xb(kb)
c                  2       xy   ka=x  kb=y    ami=xa(kb)  ani=xb(kb)
c                  3       xz   ka=x  kb=z    ami=xa(kb)  ani=xb(kb)
c
c                  4       yx   ka=y  kb=x    ami=xa(kb)  ani=xb(kb)
c                  5       yy   ka=y  kb=y    ami=xa(kb)  ani=xb(kb)
c                  6       yz   ka=y  kb=z    ami=xa(kb)  ani=xb(kb)
c
c                  7       zx   ka=z  kb=x    ami=xa(kb)  ani=xb(kb)
c                  8       zy   ka=z  kb=y    ami=xa(kb)  ani=xb(kb)
c                  9       zz   ka=z  kb=z    ami=xa(kb)  ani=xb(kb)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 410 ic=1,3
         ii=ic*(ic-1)/2
         lic=(ic-1)*ln
         do 410 is=1,3
            jj=is*(is-1)/2
            lis=(is-1)*ln
            ics=ii+is
            if(is.gt.ic) ics=jj+ic
            loc=(ics-1)*ln
            locn=( (ic-1)*3+is-1 )*ln + 18*ln
            do 411 l=1,ln
               tab_mp(l+locn)=rpp0(l+loc) - (xa(ic)-xb(ic))*rp01(l+lis)
 411        continue
 410  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 650 icase=1,9
c
         lcas=(icase-1)*ln
c*
         ka=kai(icase)
         kb=kbi(icase)
         ami=xa(ka)
         ani=xb(ka)
c
ckw*****if(nogiao.ne.0) then
         ami=ami-gauger(ka)
         ani=ani-gauger(ka)
ckw*****endif
c***********************************************************
         ij=0
         do 61 i=ia,ie
            if(kb.eq.1)  fmib=float(mulx(i))
            if(kb.eq.2)  fmib=float(muly(i))
            if(kb.eq.3)  fmib=float(mulz(i))
            do 61 j=ja,je
               if(kb.eq.1)  fnib=float(mulx(j))
               if(kb.eq.2)  fnib=float(muly(j))
               if(kb.eq.3)  fnib=float(mulz(j))
               ij=ij+1
               rppm(ij+lcas)= zero
c------it is canceled out by corresponding part of H01N-----C
c        only if the GIAO is used
     *             +    ani*fmib*rm10( (kb-1)*ln+ij )
     *             +    ami*fnib*rm01( (kb-1)*ln+ij )
     1             -    ani*twoa*rp10( (kb-1)*ln+ij )
     1             -    ami*twob*rp01( (kb-1)*ln+ij )
c------it is canceled out by corresponding part of H01N-----C
     2  + float(idk(ka,kb))*xin0(ij)
     2  + fmib*tab_mp( ((kb-1)*3+ka-1)*ln  +ij )
     2  + fnib*tab_mp( ((ka-1)*3+kb-1)*ln  +ij + 9*ln)
     3   -twoa*tab_mp( ((kb-1)*3+ka-1)*ln+ij + 18*ln)
     3   -twob*tab_mp( ((ka-1)*3+kb-1)*ln+ij + 18*ln)
  61     continue
c......  beginning of change
 650  continue
c
c changed by k.wolinski old  h11n,xy->h11n,yx etc. for off-diag
c               order of perturbations in h11n was assumed
c 01-03-1991    as h-first and mi-second, while the true
c               order is exactly opposite.!!!
c         the difference affects only off-diagonal elements
c                 of diamagnetic part of shielding.
c
      if(nogiao.eq.0) then
         do 62 l=1,ln
            xint(l)=xint(l) + rppm(l+4*ln) + rppm(l+8*ln)
            xint(l+ln)=xint(l+ln) - rppm(l+3*ln)
            xint(l+2*ln)=xint(l+2*ln) - rppm(l+6*ln)
            xint(l+3*ln)=xint(l+3*ln) - rppm(l+ln)
            xint(l+4*ln)=xint(l+4*ln) + rppm(l) + rppm(l+8*ln)
            xint(l+5*ln)=xint(l+5*ln) - rppm(l+7*ln)
            xint(l+6*ln)=xint(l+6*ln) - rppm(l+2*ln)
            xint(l+7*ln)=xint(l+7*ln) - rppm(l+5*ln)
            xint(l+8*ln)=xint(l+8*ln) + rppm(l) + rppm(l+4*ln)
  62     continue
c
      else
c
         do 63 l=1,ln
            xint(l)=        + rppm(l+4*ln) + rppm(l+8*ln)
            xint(l+ln)=           - rppm(l+3*ln)
            xint(l+2*ln)=             - rppm(l+6*ln)
            xint(l+3*ln)=             - rppm(l+ln)
            xint(l+4*ln)=             + rppm(l) + rppm(l+8*ln)
            xint(l+5*ln)=             - rppm(l+7*ln)
            xint(l+6*ln)=             - rppm(l+2*ln)
            xint(l+7*ln)=             - rppm(l+5*ln)
            xint(l+8*ln)=             + rppm(l) + rppm(l+4*ln)
  63     continue
      endif
c**********************************************************************
c                        begining for
c
c                     (mi0 | h01n | ni0)
c
c**********************************************************************
      do 550 icase=1,3
         lcas=(icase-1)*ln
         lcah=lcas + 9*ln
c        if(icase.eq.1) then
            kc=2
            ks=3
c        endif
         if(icase.eq.2) then
            kc=3
            ks=1
         endif
         if(icase.eq.3) then
            kc=1
            ks=2
         endif
         cmi=xa(kc)
         cni=xb(kc)
         smi=xa(ks)
         sni=xb(ks)
c*
         loc1=(kc-1)*ln
         loc2=(ks-1)*ln
         loc3=( (kc-1)*3+ks-1)*ln
         loc4=( (ks-1)*3+kc-1)*ln
c*
         ij=0
         do 51 i=ia,ie
            fmic=float(mulx(i))
            if(kc.eq.2)  fmic=float(muly(i))
            if(kc.eq.3)  fmic=float(mulz(i))
c
            fmis=float(mulx(i))
            if(ks.eq.2)  fmis=float(muly(i))
            if(ks.eq.3)  fmis=float(mulz(i))
            do 51 j=ja,je
               ij=ij+1
               fnic=float(mulx(j))
               if(kc.eq.2)  fnic=float(muly(j))
               if(kc.eq.3)  fnic=float(mulz(j))
c
               fnis=float(mulx(j))
               if(ks.eq.2)  fnis=float(muly(j))
               if(ks.eq.3)  fnis=float(mulz(j))
c*
      xint(lcah+ij)=fmic*fnis*tab_mm(loc3+ij)-fmis*fnic*tab_mm(loc4+ij)
     *  +twob*( -fmic*tab_mp(loc3+ij) + fmis*tab_mp(loc4+ij) )
     *  +twoa*( -fnis*tab_mp(loc3+ij+9*ln) + fnic*tab_mp(loc4+ij+9*ln))
     *  +twoa*twob*( (smi-sni)*rp10(loc1+ij)-(cmi-cni)*rp10(loc2+ij) )
c*
  51     continue
 550  continue
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
      end
c=============================================================
      subroutine incsh3(ityp,jtyp,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,
     *                  xint,xra,lxin)

      use memory

      implicit real*8 (a-h,o-z)
      dimension xa(3),xb(3),xra(3)
      dimension la(3),lb(3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /tape/ inp,inp2,iout,it  ,jt,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
c----------------------------------------------------------------------
c This routine calculates integrals needed for 1st-order el-nuc.attraction
c          V10 mi,ni = (mi1|V|ni0) + (mi0|V|ni1)
c over the GIAO basis functions.
c----------------------------------------------------------------------
c This routine was changed by KW in February, 1999  Final integrals are
c formed on the contracted shell level in the intsh3 routine.
c----------------------------------------------------------------------
c output integrals :
c
      dimension xint(lxin)    ! 4*len
c
c xint(1-ln)    contains  {mi(x+1)|r1n**-1|ni0} integrals
c xint(ln-2ln)  contains  {mi(y+1)|r1n**-1|ni0} integrals
c xint(2ln-3ln) contains  {mi(z+1)|r1n**-1|ni0} integrals
c xint(3ln-4ln) contains  {mi0    |r1n**-1|ni0} integrals
c----------------------------------------------------------------------
      dimension ndeg(8)
      data ndeg/0,1,1,2,3,4,5,6/
      dimension len(8)
      data len/1,3,4,6,10,15,21,28/
      dimension  nfu(9)
      data nfu/0,1,4,8,14,24,39,60,88/
c----------------------------------------------------------------------
      twoa=two*a
      twob=two*b
      ln=len(ityp)*len(jtyp)
      ia=nfu(ityp)+1
      ie=nfu(ityp+1)
      ja=nfu(jtyp)+1
      je=nfu(jtyp+1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  for a given pair of primitive shells and for a given reference atom
c
c     the number of roots , xrys and roots  for rys polinominal
c
c                     are calculated
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nroots=(ndeg(ityp)+ndeg(jtyp)+4)/2
c
      ax=xa(1)
      ay=xa(2)
      az=xa(3)
c
      bx=xb(1)
      by=xb(2)
      bz=xb(3)
c
      cx=xra(1)
      cy=xra(2)
      cz=xra(3)
c
      p=a+b
      p1=one/p
      px=(a*ax+b*bx)*p1
      py=(a*ay+b*by)*p1
      pz=(a*az+b*bz)*p1
      rpc2=(px-cx)**2+(py-cy)**2+(pz-cz)**2
      xrys=p*rpc2
c
c     calculate roots
c
      if (nroots.le.3) call rt123
      if (nroots.eq.4) call root4
      if (nroots.eq.5) call root5
      if (nroots.gt.5) call root6
c**********************************************************************
          do 34 i=1,3
          la(i)=0
 34       lb(i)=0
c
c calculate {mi0|r1n**-1|ni0} integrals ; put in xint(3*ln+1 - 4*ln)
c
      call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),xb(3),
     *            xra(1),xra(2),xra(3),la,lb,k1,xint(1+3*ln),ln)
c
c**********************************************************************
c
c   {mi(t+1)|r1n**-1|ni0} integrals ;t=x,y,z  /3*ln integrals/
c
c               in xint(3*ln)      /3,4/
c***********************                   ****************************
c
          do 110 i=1,3
          loc=(i-1)*ln
          la(i)=+1
      call nuch01(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),xb(3),
     *            xra(1),xra(2),xra(3),la,lb,k1,xint(1+loc),ln)
          la(i)=0
 110      continue
c
      end
c=============================================================
