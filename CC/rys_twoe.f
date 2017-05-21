      subroutine twoe(sij,pij,basdat,indxx,indyy,
     1                indzz,fact,nonz,buf)
      implicit real*8 (a-h,o-z)
      dimension sij(*),pij(*),basdat(13,*),indxx(*),indyy(*),indzz(*),
     1 buf(*)
      dimension buf1(50625)
c
c     calculation of the two electron integrals
c
      logical nonz,first,firstc
      common /two1/ iityp,jjtyp,kktyp,lltyp,
     1 ityp,jtyp,ktyp,ltyp,ncount,ncnt,nord,maxl,iwy,ilen1,
     1 jlen1,klen1,len1,ia,ie,ja,je,ka,ke,la,le,lx,eps,eps1,ca3,
     2 ca2,ca1,ca0
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension xp(3),xq(3),xa(3),xb(3),xc(3),xd(3)
      dimension xl(3)
      equivalence (xp(1),px),(xp(2),py), (xp(3),pz), (xq(1),qx),
     1 (xq(2),qy), (xq(3),qz),
     3 (xa(1),ax), (xa(2),ay), (xa(3),az), (xb(1),bx),
     4 (xb(2),by), (xb(3),bz), (xc(1),cx), (xc(2),cy),
     5 (xc(3),cz), (xd(1),dx), (xd(2),dy), (xd(3),dz)
c      common bl(1),buf(1296),sij(1200),indxx(1296),indyy(1296),
c     1 indzz(1296),buf1(1296)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(16)
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots,
     1             nroot2,ie1,ie2,ie3,ie4,ie12,ie34
      save /rys/
      dimension itm275(5),itm55(11),itm11(5),itm5(11)
      dimension gx(11,11),gy(11,11),gz(11,11),xf1(605),yf1(605),zf1(605)
      dimension xfc(1375),yfc(1375),zfc(1375)
c     *** xf1,yf1,zf1 corresponds to the 3-index fortran array xf1(11,5,*)
c     *** xfc,yfc,zfc corresponds to xfc(5,11,5,5)
c     *** these dimensions are good up to and including g functions
c  this equivalence is for printing purposes only
      dimension ndege(6)
      data itm275/0,275,550,825,1100/
      data itm55/0,55,110,165,220,275,330,385,440,495,550/
      data itm11/0,11,22,33,44/
      data itm5/0,5,10,15,20,25,30,35,40,45,50/
c     *** this is 2*pi**(5/2)
      data con1 /34.98683665265d0/
      data ndege/1,2,2,3,4,5/
        logical exists
c
c      inquire(file='testout',exist=exists)
c       if(.not.exists) open(61,file='testout',status='unknown')
      firstc=.true.
      nonz=.false.
c
c      print *,'types',ityp,jtyp,ktyp,ltyp
c       print *,'ka,ke,la,le,ia,ie,ja,je',ka,ke,la,le,ia,ie,ja,je
c
c     inner loop over individual integrals
      do 390 ks=ka,ke
        cc=basdat(1,ks)
          csk=basdat(2,ks)
          cpk=basdat(3,ks)
        coefk=csk
        if(ktyp.eq.3) then
            coefk=cpk
          csk=csk/cpk
          fr3=csk
        endif
        do 300 l=1,3
300         xc(l)=basdat(10+l,ks)
cpp
c       print *, 'cc,csk,cpk,xc',cc,csk,cpk,xc
c
        do 390 ls=la,le
          csl=basdat(2,ls)
            cpl=basdat(3,ls)
            dd=basdat(1,ls)
          qq=cc+dd
            qqr=one/qq
          coefl=csl
          if(ltyp.eq.3) then
              fr4=csl/cpl
            coefl=cpl
            end if
          coefl=coefk*coefl
          r=zero
          do 310 l=1,3
            x1=xc(l)
            xd(l)=basdat(10+l,ls)
            x2=xd(l)
            xq(l)=(cc*x1+dd*x2)*qqr
            r=r+(x1-x2)**2
  310     continue
c       write(61,*) 'dd,csl,cpl,xd=',dd,csl,cpl,xd
          e=cc*dd*qqr
          scd=sqrt(sqrt(cc*dd))**3*exp(-e*r)*fact
          ji=0
          ij3=-3   !index of the p array
c
          do 381 is=ia,ie
             aa=basdat(1,is)
               csi=basdat(2,is)
               cpi=basdat(3,is)
             coefi=csi
             if (ityp.eq.3) then
               coefi=cpi
               csi=csi/cpi
               fr1=csi
             endif
             coefi=coefi*coefl
             do  320 l=1,3
 320           xa(l)=basdat(10+l,is)
             do 380 js=ja,je
               ji=ji+1
               ij3=ij3+3
               bb=basdat(1,js)
                 csj=basdat(2,js)
                 cpj=basdat(3,js)
               coefj=csj
               if(jtyp.eq.3) then
                   coefj=cpj
                 csj=csj/cpj
                 fr2=csj
               endif
               do 330 l=1,3
 330             xb(l)=basdat(10+l,js)
cpp
c                print *, 'overlaps=',scd,sij(ji)
                 sabcd=scd*sij(ji)
                 pp=aa+bb
                 do 360 l=1,3
                   xp(l)=pij(ij3+l)
  360            continue
cpp
c      write(61,7777) ia,ie,ja,je,ka,ke,la,le,nord,aa,bb,cc,dd,pp,qq,
c     1 scd,sij(ji),xa,xb,xc,xd,xp,xq,coefk,coefl,coefi,coefj,
c     2 e,r
 7777 format(' test',9i3,8f10.6,/,1x,12f10.6,/,1x,12f10.6)
c
c     inner loop over individual integrals
                 pq=pp*qq
                 ppq=pp+qq
                 rho=pq/ppq
c
c
c     calculate number of roots and xrys for rys
c
                 rr1=((px-qx)**2+(py-qy)**2+(pz-qz)**2)
                 xrys=rr1*rho
                 if (rr1.lt.acc) rr1=acc
                 facmax=min(one/rr1,rho*1.273239545d0)
                 if (sabcd*sqrt(facmax).lt.eps) go to 380
c     If these is at least one contraction calculated, set nonz to true
                 nonz=.true.
  364            sabcd=sabcd*coefi*coefj
                 const=con1*sabcd/(pq*sqrt(ppq))
                 consts=const*0.88622692545275801d0
c *** o.88622.. is sqrt(pi)/2   con1=2pi**(5/2)
c
c     calculate roots
c
                 first=.true.
c  For (ss|ss) and (ps|ss), nroots=1. For (ss|ss), nroot2=2; for (ps|ss), nroot2=3
                 if (nroots.gt.1) go to 190
                 if(nroot2.gt.2) go to 100
                 call ft0
                 buf1(1)=consts*f0
c       write(61,*) 'consts,f0,buf1,rr1,xrys',consts,f0,buf1(1),rr1,xrys
                 go to 203
c *** the above was a special code for (ss ss) integrals
100              continue
c *** special code for (ps ss), (sp ss) (ss ps) and (ss sp)
                 call ft0
                 intct=0
                 if (ityp.gt.1) then
                   den=aa+bb
                   do 110 l=1,3
110                xl(l)=bb*(xa(l)-xb(l))*f0
                   cs=csi
                   if(ityp.eq.3) intct=1
                 else if (jtyp.gt.1) then
                   den=aa+bb
                   do 102 l=1,3
102                xl(l)=aa*(xb(l)-xa(l))*f0
                   cs=csj
                   if(jtyp.eq.3) intct=1
                 else if (ktyp.gt.1) then
                   den=cc+dd
                   do 103 l=1,3
103                xl(l)=dd*(xc(l)-xd(l))*f0
                   cs=csk
                   f1=-f1
                   if(ktyp.eq.3) intct=1
                 else
                   den=cc+dd
                   do 104 l=1,3
104                xl(l)=cc*(xd(l)-xc(l))*f0
                   cs=fr4
                   f1=-f1
                   if (ltyp.eq.3) intct=1
                 endif
                   denr=one/den
                 buf1(1)=consts*f0*cs
                 do 105 l=1,3
                 buf1(intct+l)=-consts*(xl(l)+rho*f1*(xp(l)-xq(l)))*denr
c                 buf1(intct+l)=-consts*(xl(l)+rho*f1*(xp(l)-xq(l)))/den
105              continue
                 go to 203
190              continue
                 if (nroots.le.3) call rt123
                 if (nroots.eq.4) call root4
                 if (nroots.eq.5) call root5
                 if (nroots.gt.5) call root6
                   do 65 iroot=1,nroots
                   rw=rysw(iroot)*const
c
c
c     calculate two dimensional integrals
c
                   b00=rysr(iroot)/((one+rysr(iroot))*ppq)
                   b00q=b00*qq
                   b00p=b00*pp
c                   b10=one/(pp+pp)-b00q/pp
                   b10=half*(one-b00q)/pp
c                   b01=one/(qq+qq)-b00p/qq
                   b01=half*(one-b00p)*qqr
                   b00=b00*half
                   c0x=px-ax+(qx-px)*b00q
                   cx0=qx-cx+(px-qx)*b00p
                   c0y=py-ay+(qy-py)*b00q
                   cy0=qy-cy+(py-qy)*b00p
                   c0z=pz-az+(qz-pz)*b00q
                   cz0=qz-cz+(pz-qz)*b00p
                   gx(1,1)=one
                   gy(1,1)=one
                   gz(1,1)=rw
                   gx(2,1)=c0x
                   gx(1,2)=cx0
                   gx(2,2)=b00+cx0*c0x
                   gy(2,1)=c0y
                   gy(1,2)=cy0
                   gy(2,2)=b00+cy0*c0y
                   gz(2,1)=c0z*rw
                   gz(1,2)=cz0*rw
                   gz(2,2)=(b00+c0z*cz0)*rw
                   b10t=zero
                   b00t=b00
cpp
c      write(61,1116) b00,b10,b01
cpp
1116  format(' after b00 b10 b01', 3f10.6)
                   do 10 ijd=3,ie12
                     b10t=b10t+b10
                     b00t=b00t+b00
                     gx(ijd,1)=b10t*gx(ijd-2,1)+c0x*gx(ijd-1,1)
                     gy(ijd,1)=b10t*gy(ijd-2,1)+c0y*gy(ijd-1,1)
                     gz(ijd,1)=b10t*gz(ijd-2,1)+c0z*gz(ijd-1,1)
                     if(ie34.eq.1) go to 10
                     gx(ijd,2)=b00t*gx(ijd-1,1)+cx0*gx(ijd,1)
                     gy(ijd,2)=b00t*gy(ijd-1,1)+cy0*gy(ijd,1)
                     gz(ijd,2)=b00t*gz(ijd-1,1)+cz0*gz(ijd,1)
10                 continue
                   b01t=zero
                   b00t=b00
                   do 20 kld=3,ie34
                     b01t=b01t+b01
                     b00t=b00t+b00
                     gx(1,kld)=b01t*gx(1,kld-2)+cx0*gx(1,kld-1)
                     gy(1,kld)=b01t*gy(1,kld-2)+cy0*gy(1,kld-1)
                     gz(1,kld)=b01t*gz(1,kld-2)+cz0*gz(1,kld-1)
                     if(ie12.eq.1) go to 20
                     gx(2,kld)=b00t*gx(1,kld-1)+c0x*gx(1,kld)
                     gy(2,kld)=b00t*gy(1,kld-1)+c0y*gy(1,kld)
                     gz(2,kld)=b00t*gz(1,kld-1)+c0z*gz(1,kld)
20                 continue
                   b10t=zero
                   do 25 ijd=3,ie12
                     b10t=b10t+b10
                     b00t=b00
                     do 25 kld=3,ie34
                     b00t=b00t+b00
                    gx(ijd,kld)=b10t*gx(ijd-2,kld)+b00t*gx(ijd-1,kld-1)+
     1                 c0x*gx(ijd-1,kld)
                    gy(ijd,kld)=b10t*gy(ijd-2,kld)+b00t*gy(ijd-1,kld-1)+
     1                 c0y*gy(ijd-1,kld)
                    gz(ijd,kld)=b10t*gz(ijd-2,kld)+b00t*gz(ijd-1,kld-1)+
     1                 c0z*gz(ijd-1,kld)
25                 continue
                   do 30 i=1,ie12
                     ii1=itm55(i)
                     do 30 k=1,ie34
                       xf1(ii1+k)=gx(i,k)
                       yf1(ii1+k)=gy(i,k)
                       zf1(ii1+k)=gz(i,k)
30                 continue
                   ie11=ie12
                   do 40 j=2,ie2
                     ie11=ie11-1
                     do 40 i=1,ie11
                       ij=itm55(i)+itm11(j)
c     *** (1,j,i) index
c     *** (i,j-1,i+1) index  is ij+44
c     *** 1,j-1,i) index  is ij-11
                       do 40 k=1,ie34
                         xf1(ij+k)=xf1(ij+44+k)+(ax-bx)*xf1(ij-11+k)
                         yf1(ij+k)=yf1(ij+44+k)+(ay-by)*yf1(ij-11+k)
                         zf1(ij+k)=zf1(ij+44+k)+(az-bz)*zf1(ij-11+k)
c***  xf1(k,j,i)=xf1(k,j-1,i+1)+(ax-bx)*xf1(k,j-1,i)
c
40                 continue
                   l=1
                   do 50 i=1,ie1
                     ii1=itm275(i)
                     ii2=itm55(i)
                     do 50 j=1,ie2
                       ij1=ii1+itm55(j)
                       ij2=ii2+itm11(j)
                       do 50 k=1,ie34
                         xfc(ij1+itm5(k)+1)=xf1(ij2+k)
                         yfc(ij1+itm5(k)+1)=yf1(ij2+k)
                         zfc(ij1+itm5(k)+1)=zf1(ij2+k)
c
c***  xfc(1,k,j,i)=xf1(k,j,i)
c
cpp
c      write(61,*) '2dint',i,j,k,l,xfc(ij1+itm5(k)+1),yfc(ij1+itm5(k)+1),
c     1  zfc(ij1+itm5(k)+1),rw
cpp     
50                 continue
                   do 60 i=1,ie1
                     ii4=itm275(i)
                     do 60 j=1,ie2
                       ij4=ii4+itm55(j)
                       ie33=ie34
                       do 60 l=2,ie4
                         ijl4=ij4+l
                         ie33=ie33-1
                         do 60 k=1,ie33
                           ijk4=ijl4+itm5(k)
                           xfc(ijk4)=xfc(ijk4+4)+(cx-dx)*xfc(ijk4-1)
                           yfc(ijk4)=yfc(ijk4+4)+(cy-dy)*yfc(ijk4-1)
                           zfc(ijk4)=zfc(ijk4+4)+(cz-dz)*zfc(ijk4-1)
c
c***  xfc(l,k,j,i)=xfc(l-1,k+1,j,i)+(cx-dx)*xfc(l-1,k,j,i)
c
c      write(61,9999) i,j,k,l,xfc(ijk4),yfc(ijk4),zfc(ijk4),rw
9999  format('factors',4i5,3f10.6,f12.8)
60                 continue
c
c     assembly of the two electron integrals
c
c     offset the indentation by 4 places to the left
               if (first) then
                 do 75 ict1=1,ncount
                   buf1(ict1)=xfc(indxx(ict1))*yfc(indyy(ict1))*
     1                         zfc(indzz(ict1))
75               continue
                 first=.false.
               else
               do 80 ict=1,ncount
                 buf1(ict)=buf1(ict)+xfc(indxx(ict))*yfc(indyy(ict))*
     1           zfc(indzz(ict))
80             continue
             endif
65         continue
c  offset back by 4 places to the right
               if(iityp.eq.3) then
                 nco=0
                 do 120 jxj=1,jlen1
                   do 120 kxk=1,klen1
                     do 120 lxl=1,len1
                       nco=nco+1
                       buf1(nco)=buf1(nco)*fr1
120              continue
               endif
               if(jjtyp.eq.3) then
                 nco=0
                 jinc=(jlen1-1)*klen1*len1
                 do 150 i=1,ilen1
                   do 140 k=1,klen1
                     do 140 l=1,len1
                       nco=nco+1
                       buf1(nco)=buf1(nco)*fr2
140                continue
                   nco=nco+jinc
150              continue
               endif
               if(kktyp.eq.3) then
               nco=0
               kinc=(klen1-1)*len1
               do 180 i=1,ilen1
                 do 180 j=1,jlen1
                   do 170 l=1,len1
                     nco=nco+1
                     buf1(nco)=buf1(nco)*fr3
170                continue
                   nco=nco+kinc
180            continue
             endif
             if(lltyp.eq.3) then
             nco=1
             do 200 i=1,ilen1
               do 200 j=1,jlen1
                 do 200 k=1,klen1
                   buf1(nco)=buf1(nco)*fr4
                   nco=nco+len1
200          continue
           endif
203        if (firstc) then
             do 205 icx=1,ncount
               buf(icx)=buf1(icx)
205          continue
             firstc=.false.
           else
             do 210 icx=1,ncount
               buf(icx)=buf(icx)+buf1(icx)
210          continue
           endif
           go to 380
220        if (firstc) then
             do 230 ixc=1,ncount
230          buf(ixc)=zero
           endif
380      continue
381    continue
390   continue
c      write(61,4444) (buf(i),i=1,ncount)
4444  format(1x,' buf',(5f13.8))
c
      end
