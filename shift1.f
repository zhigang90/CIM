c===================================================================
c Feb.22,96, K.Wolinski:
c
c This file contains now the old int123.f file(intsh1,2,3 routines)
c===================================================================
      subroutine  shift1(bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /aux/ natom,nucnr(200)
      common /cux/ lden,lfoc,lforc,love,lh01, lhfc,lval,lvec,lrem
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c external field info :
      common /fieldx/ xfield,xfgrad,elfiel(9)
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /gauge/ gauger(3)
c-----
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------
      call mmark
c------------------------------------------------
      ifield=xfield
      ifgrad=xfgrad
c------------------------------------------------
c memory reservation
c
      call memres1(bl)
c------------------------------------------------
c set up the gauge-origin for both the CPHF and GIAO
c
       call rgauge(bl(inuc),ngauge,gauger)
c------------------------------------------------
c check-run only or lack of scf convergance :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c------------------------------------------------
      ntri=ncf*(ncf+1)/2
c------------------------------------------------
      call secund(t3b)
      call elapsec(et3b)
c------------------------------------------------
c check if COSMO flag is on and get corresponding data
c
      call getCosmoData(na,bl,bl(inuc))
c
      call getival('cosmo',icosmo)
      if(icosmo.eq.0) then
        inucur=inuc
        natcur=na
      else
        call getival('inucosmo',inucosmo)
        call getival('c_nps',nps)       ! no. of surface segments
        inucur=inucosmo
        natcur=na+nps
      endif
ctest
c       inucur=inuc
c       natcur=na
c------------------------------------------------
      call intsh3 (natcur,bl,inx,last,nprint,iout,bl(ibas),
     *             bl(inucur),ncs,ncf,ntri,bl(lfoc))
c
c     call intsh3 (na,bl,inx,last,nprint,iout,bl(ibas),
c    *             bl(inuc),ncs,ncf,ntri,lfoc)
c------------------------------------------------
      call secund(t1b)
      call elapsec(et1b)
c------------------------------------------------
c change here - removed ifield, ifgrad, elfiel
c
      call intsh1(na,bl,inx,last,nprint,ifield,ifgrad,elfiel,
     *             iout,bl(ibas),bl(inuc),ncs,ncf,ntri,
     *             bl(love),bl(lfoc))
c
c  pseudopotential contribution
c
      call getival('npsp',npsp)
      if(npsp.ne.0.and.nogiao.eq.0)then
        call pspFnmr(na,npsp,ncf,ncs,inx,bl(inuc),bl(ibas),
     $               ntri,bl(lfoc))
        if(nprint.gt.2) then
          write(iout,*)'After psp contribution'
          call druma(bl(lfoc),ncf,iout,'fockmatx')
          call druma(bl(lfoc+ntri),ncf,iout,'fockmaty')
          call druma(bl(lfoc+2*ntri),ncf,iout,'fockmatz')
        endif
      endif
c------------------------------------------------
      call secund(t1end)
      call elapsec(et1end)
c------------------------------------------------
c2002 write H1 and S1 on a disk : unit 61 & 62
c
      irec=1
      ntri3=ntri*3
      call save1mat(61,irec,ntri3,bl(lfoc))
      call save1mat(62,irec,ntri3,bl(love))
c------------------------------------------------
c print timing
c
      call timesh(t3b,t1b,t1end,et3b,et1b,et1end,iout,icond)
c------------------------------------------------
      call retmark
c
      end
c===================================================================
      subroutine memres1(bl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /cux/ lden,lfoc,lforc,love,lh01, lhfc,lval,lvec,lrem
c--
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c----------
      dimension bl(*)
c------------------------------------------------
c lforc is the starting address for the nuclear shielding tensor
c
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
      call getmem(ntri3,lfoc )
      call getmem(ntri3,love )
c------------------------------------------------
c zero out : Forc, Fock, Overlap, H01
c
      call zeroit(bl(lfoc),ntri3)
      call zeroit(bl(love),ntri3)
c------------------------------------------------
c update value of last  :
c
      call getmem(1,last)
      call retmem(1)
c------------------------------------------------
c     if(last+ntri.gt.lcore) then
c       ioffset=igetival('ioffset')
c       write (iout,300) last+ntri-ioffset,lcore-ioffset
c     call nerror(1,'Shift1',
c    *            'Memory: needed and available ',
c    *             last+ntri-ioffset,lcore-ioffset)
c     endif
c 300 format (/1x,'common bl too small for shift1 run, required =',i10,
c    *3x,' available =',i10/)
c------------------------------------------------
      end
c===================================================================
      subroutine timesh(t3b,t1b,t1end,et3b,et1b,et1end,iout,icond)
      implicit real*8 (a-h,o-z)
      sixty=60.0d0
      tvvv =(t1b-t3b)/sixty
      etvvv=(et1b-et3b)/sixty
      tsth =(t1end-t1b)/sixty
      etsth=(et1end-et1b)/sixty
      total=(t1end-t3b)/sixty
      etotal=(et1end-et3b)/sixty
c
      write(iout,400) total,etotal
c2002 write(iout,100) tsth,etsth
c2002 write(iout,200) tvvv,etvvv
      call f_lush(iout)
c
c     write(iout,500)
c
  100 format(24x,'s10, t10, h10 =',f9.2,11x,f9.2,' min')
  200 format(24x,'v10           =',f9.2,11x,f9.2,' min')
  400 format(/'Master CPU time for GIAO 1e integrals = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
  500 format(58('-'))
c
      end
c===================================================================
      subroutine rgauge(datnuc,ngauge,gauger)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
      dimension gauger(3)
      data zero /0.d0/
c
      if(ngauge.eq.0) then
         gauger(1)=zero
         gauger(2)=zero
         gauger(3)=zero
      else
         gauger(1)=datnuc(2,ngauge)
         gauger(2)=datnuc(3,ngauge)
         gauger(3)=datnuc(4,ngauge)
      endif
c
      end
c======================================================================
      subroutine intsh1(nna,bl,inx,last,iprint,ifield,ifgrad,xfld,
     *                   iout,datbas,datnuc,ncs,ncf,ntri,ove1,foc1)
      implicit real*8 (a-h,o-z)
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension bl(*)
      dimension inx(12,*)
c             !  28*28*4  28*28*3
      dimension s4(3136),s3(2352),xra(3)
      dimension xfld(9)
      dimension ove1(ntri,3),foc1(ntri,3)
c--------------------------------
c remember value of last
      lenter=last
c--------------------------------
      xra(1)=zero
      xra(2)=zero
      xra(3)=zero
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         xi=datnuc(2,iat)
         yi=datnuc(3,iat)
         zi=datnuc(4,iat)
         ngci=inx(4,i)
         do 55 igc=0,ngci
c
         jfu=0
         do 50 j=1,i
            len2=inx(3,j)
            jat=inx(2,j)
            xj=datnuc(2,jat)
            yj=datnuc(3,jat)
            zj=datnuc(4,jat)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
c
             xx=xi-xj
             yy=yi-yj
             zz=zi-zj
             xij= yi*zj-zi*yj
             yij=-xi*zj+zi*xj
             zij= xi*yj-yi*xj
             len=len1*len2
c----
        do 45 jgc=0,ngcjx
c----
             do 10 l=1,6*len
  10         bl(last+l)=zero
c
             ls4=4*len
             ls3=3*len
             llx=last+3*len
c***************************************
c
       key=1
c
c      (i10 /j00) + (i00 /j10)  matrix
c
c           multiplicator = +i/2c
c***************************************
c giao
      if(nogiao.eq.0) then
        if(iat.eq.jat) go to 12
          do 11 l=1,ls4
   11     s4(l)=zero
          ld=4
c
          call onesh(i,j,igc,jgc,datbas,bl,s4,inx,key,nna,kk,k2,
     *               xra,ls4,ld)
c
          do 120 l=1,len
          ll=last+l
          bl(ll)=      -zz*s4(l+len) + yy*s4(l+2*len) + xij*s4(l+3*len)
          bl(ll+len)=  +zz*s4(l)      -xx*s4(l+2*len) + yij*s4(l+3*len)
          bl(ll+2*len)=-yy*s4(l)        +xx*s4(l+len) + zij*s4(l+3*len)
 120      continue
   12   continue
      endif
c*******************************
c
        key=3
c
c       (i00 / h10 /j00) matrix
c
c           multiplicator = -i/2c
c*******************************
        do 13 l=1,ls3
   13   s3(l)=zero
        ld=3
c-
        call onesh(i,j,igc,jgc,datbas,bl,s3,inx,key,nna,kk,k2,
     *             xra,ls3,ld)
c-
        do 110 l=1,len
        ll=llx+l
        bl(ll)=bl(ll) - s3(l)
        bl(ll+len)=bl(ll+len) - s3(l+len)
        bl(ll+2*len)=bl(ll+2*len) - s3(l+2*len)
 110    continue
c
c**********************************************
c
       key=2
c
c      (i10 / t /j00) + (i00 / t /j10)  matrix
c
c           multiplicator = +i/2c
c**********************************************
c giao
      if(nogiao.eq.0) then
        if(iat.eq.jat) go to 15
c
         ld=4
         do 14 l=1,ls4
   14    s4(l)=zero
c----
         call onesh(i,j,igc,jgc,datbas,bl,s4,inx,key,nna,kk,k2,
     *              xra,ls4,ld)
c----
         do 150 l=1,len
         ll=llx+l
         bl(ll)=bl(ll) -zz*s4(l+len) + yy*s4(l+2*len) + xij*s4(l+3*len)
         bl(ll+len)=bl(ll+len) +zz*s4(l)-xx*s4(l+2*len)+yij*s4(l+3*len)
         bl(ll+2*len)=bl(ll+2*len)-yy*s4(l)+xx*s4(l+len)+zij*s4(l+3*len)
 150     continue
c
c***************************************
c for an external electric field
c
c      (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)  matrix
c
c and for electric field gradient
c
c      (i10 /xy/ j00) + (i00 /xy/ j10)  matrix
c
c where xy=xx,yx,yy,zx,zy,zz
c
c           multiplicator = +i/2c
c***************************************
c for an external electric field
        if(ifield.eq.1) then
           key=1
           do 16 k1=1,3
           field=xfld(k1)
           if(field.ne.zero) then
             xxf=xx*field
             yyf=yy*field
             zzf=zz*field
             xijf=xij*field
             yijf=yij*field
             zijf=zij*field
               do 17 l=1,ls4
               s4(l)=zero
   17          continue
c------
              call onesh(i,j,igc,jgc,datbas,bl,s4,inx,key,nna,k1,k2,
     *                   xra,ls4,ld)
c------
               do 160 l=1,len
               ll=llx+l
               bl(ll)=bl(ll)
     *                -zzf*s4(l+len)+yyf*s4(l+2*len)+xijf*s4(l+3*len)
               bl(ll+len)=bl(ll+len)
     *                +zzf*s4(l)-xxf*s4(l+2*len)+yijf*s4(l+3*len)
               bl(ll+2*len)=bl(ll+2*len)
     *                -yyf*s4(l)+xxf*s4(l+len)+zijf*s4(l+3*len)
 160           continue
           endif
   16      continue
        endif
c
c for an external electric field gradient
        if(ifgrad.eq.1) then
           key=1
           do 18 k1=4,9
           field=xfld(k1)
           if(field.ne.zero) then
             xxf=xx*field
             yyf=yy*field
             zzf=zz*field
             xijf=xij*field
             yijf=yij*field
             zijf=zij*field
               do 19 l=1,ls4
               s4(l)=zero
   19          continue
c
               call onesh(i,j,igc,jgc,datbas,bl,s4,inx,key,nna,k1,k2,
     *                    xra,ls4,ld)
c
               do 180 l=1,len
               ll=llx+l
               bl(ll)=bl(ll)
     *                -zzf*s4(l+len)+yyf*s4(l+2*len)+xijf*s4(l+3*len)
               bl(ll+len)=bl(ll+len)
     *                +zzf*s4(l)-xxf*s4(l+2*len)+yijf*s4(l+3*len)
               bl(ll+2*len)=bl(ll+2*len)
     *                -yyf*s4(l)+xxf*s4(l+len)+zijf*s4(l+3*len)
  180          continue
           endif
   18      continue
        endif
c
   15   continue
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  now common multiplicator for integrals from  key=3 + key=2          c
c                                                                      c
c                        is equal  +i/2c                               c
c                                                                      c
c  these integrals give contributions to the 1-order fock matrix      c
c  which is used in coupled hartree-fock equations. from this i-order  c
c  density matrix has the factor equal +i/2c .                         c
c
c***********************************************************************
            iij=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
               ii=iff*(iff-1)/2
            do 40 j1=1,len2
               jff=jff+1
               iij=iij+1
               iijy=iij+len
               iijz=iijy+len
               if (jff.gt.iff) go to 40
               ij=ii+jff
            ove1(ij,1)=ove1(ij,1) + bl(last+iij)
            ove1(ij,2)=ove1(ij,2) + bl(last+iijy)
            ove1(ij,3)=ove1(ij,3) + bl(last+iijz)
c***
            foc1(ij,1)=foc1(ij,1) + bl(last+3*len+iij)
            foc1(ij,2)=foc1(ij,2) + bl(last+3*len+iijy)
            foc1(ij,3)=foc1(ij,3) + bl(last+3*len+iijz)
c***
   40       continue
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c*
c--------------------------------
c return value of last
      last=lenter
c--------------------------------
c***********************************************************************
c*
      if(iprint.gt.2) then
c*
      call druma(ove1(1,1),ncf,iout,'overlapx')
      call druma(ove1(1,2),ncf,iout,'overlapy')
      call druma(ove1(1,3),ncf,iout,'overlapz')
        write(iout,444)
c
        write(iout,446)
c
      call druma(foc1(1,1),ncf,iout,'fockmatx')
      call druma(foc1(1,2),ncf,iout,'fockmaty')
      call druma(foc1(1,3),ncf,iout,'fockmatz')
c*
      endif
c*
c***********************************************************************
   70 return
c
 444  format('*******************************')
 446  format(/'***  constant one-el. part of the fock matrix  ***'/)
c
      end
c======================================================================
      subroutine intsh3(nna,bl,inx,last,iprint,iout,
     *                  datbas,datnuc,ncs,ncf,ntri,foc1)
c*********************************************************************
c This routine calculates 1-el integral derivatives over the GIAO's
c                      with 1/r1n operator.
c
c                   (i10 / v /j00) + (i00 / v /j10)
c*********************************************************************
      implicit real*8 (a-h,o-z)
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension bl(1)
      dimension foc1(ntri,3)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,1)
      dimension s4(3136),xra(3)
c--------------------------------------------------------------------
c
c s4(4*28*28) (up to i|i ) used for (mi(x+1) |r-Rn|-1 ni0)
c                                   (mi(y+1) |r-Rn|-1 ni0)
c                                   (mi(z+1) |r-Rn|-1 ni0)
c                                   (mi0     |r-Rn|-1 ni0)
c--------------------------------------------------------------------
      if(nogiao.ne.0) return
c--------------------------------------------------------------------
c remember value of last
      lenter=last
c--------------------------------------------------------------------
c
      kk=0
      k2=0
c
      key=5
c
c
      ifu=0
      do 60 i=1,ncs
         iat=inx(2,i)
         xi=datnuc(2,iat)
         yi=datnuc(3,iat)
         zi=datnuc(4,iat)
         len1=inx(3,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
         jfu=0
         do 50 j=1,i
            jat=inx(2,j)
            xj=datnuc(2,jat)
            yj=datnuc(3,jat)
            zj=datnuc(4,jat)
            len2=inx(3,j)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
            len=len1*len2
            ls4= 4*len
            lastv=last
c
             xx=xi-xj
             yy=yi-yj
             zz=zi-zj
             xij= yi*zj-zi*yj
             yij=-xi*zj+zi*xj
             zij= xi*yj-yi*xj
c
         do 45 jgc=0,ngcjx
c
            do 12 l=1,ls4
  12        bl(lastv+l)=zero
c
            ld=4
c
c          loop over atoms
c
           do 250 nra=1,nna
              zza=datnuc(1,nra)
              xra(1)=datnuc(2,nra)
              xra(2)=datnuc(3,nra)
              xra(3)=datnuc(4,nra)
c---
           call onesh(i,j,igc,jgc,datbas,bl,s4,inx,key,nna,kk,k2,
     *                xra,ls4,ld)
c---
              do 169 l=1,ls4
 169          s4(l)=-zza*s4(l)
c
ckw99
c v     x,y,z
              do 188 l=1,len
                 lvx=lastv+l
                 lvy=lvx+len
                 lvz=lvy+len
          bl(lvx)=bl(lvx) -zz*s4(l+len)+yy*s4(l+2*len)+xij*s4(l+3*len)
          bl(lvy)=bl(lvy) +zz*s4(l)    -xx*s4(l+2*len)+yij*s4(l+3*len)
          bl(lvz)=bl(lvz) -yy*s4(l)    +xx*s4(l+len)  +zij*s4(l+3*len)
 188          continue
ckw99
 250       continue
c***********************************************************************
            iijxx=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
               ii=iff*(iff-1)/2
            do 40 j1=1,len2
               jff=jff+1
c
               iijxx=iijxx+1
               iijxy=iijxx+len
               iijxz=iijxy+len
c
         if(jff.gt.iff) go to 40
               ij=ii+jff
c
            foc1(ij,1)=foc1(ij,1) + bl(lastv+iijxx)
            foc1(ij,2)=foc1(ij,2) + bl(lastv+iijxy)
            foc1(ij,3)=foc1(ij,3) + bl(lastv+iijxz)
c
   40       continue
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c--------------------------------
c return value of last
      last=lenter
c--------------------------------
      end
c======================================================================
      subroutine onesh(ics,jcs,igc,jgc,datbas, bl,sb,inx,key,nn,k1,k2,
     *                 xra,lb9,ldim)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /tape/ inp,inp2,iout
      common /forcdbl/ thre1,thre2,tchf
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension  xra(3), xa(3), xb(3)
      dimension bl(*)
      dimension inx(12,1)
      dimension datbas(13,*)
c
c               ! 28*28*12
      dimension  rb(9408),sb(lb9)
c
      data twopi/0.6366197723675d0/
c-----------------------------------------------------------------------
c
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
c
      ilen=inx(3,ics)
      jlen=inx(3,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
      if (ityp.eq.6) ilen1=10
      if (jtyp.eq.6) jlen1=10
c
      if(ityp.eq.11) ilen1=15
      if(ityp.eq.12) ilen1=21
      if(ityp.eq.13) ilen1=28
c
      if(jtyp.eq.11) jlen1=15
      if(jtyp.eq.12) jlen1=21
      if(jtyp.eq.13) jlen1=28
c
      len1=ilen1*jlen1
c---
      ityp1=ityp
      jtyp1=jtyp
      if(ityp.ge.5) ityp1=ityp-1
      if(jtyp.ge.5) jtyp1=jtyp-1
      if(ityp.ge.7) ityp1=ityp-2
      if(jtyp.ge.7) jtyp1=jtyp-2
      if(ityp.eq.11) ityp1=6
      if(ityp.eq.12) ityp1=7
      if(ityp.eq.13) ityp1=8
      if(jtyp.eq.11) jtyp1=6
      if(jtyp.eq.12) jtyp1=7
      if(jtyp.eq.13) jtyp1=8
c-----------------------------------------------------------------------
      lxin=ldim*len1
      do i=1,lxin
         sb(i)=zero
      enddo
c-----------------------------------------------------------------------
c dynamical memory allocations for key=4:
c
      if(key.eq.4) then
         call getmem( 3*len1,irp01)
         call getmem( 3*len1,irp10)
         call getmem(27*len1,irpmm)
         call getmem(27*len1,irpmp)
         call getmem(27*len1,irppm)
         call getmem( 9*len1,irpp0)
         call getmem( 3*len1,irm01)
         call getmem( 3*len1,irm10)
         call getmem( 3*len1,ixin0)
         call getmem( 9*len1,i_mm )
         call getmem(27*len1,i_pmm)
         call getmem(27*len1,i_mp )
         call getmem(27*len1,i_pmp)
         call getmem(27*len1,i_ppm)
      endif
c-----------------------------------------------------------------------
c beg. and end of contraction :
c
      ia=inx(1,ics)+1
      ja=inx(1,jcs)+1
      ie=inx(5,ics)
      je=inx(5,jcs)
c-----------------------
c
      do 600 i=ia,ie
         a=datbas(1,i)
         sqa=sqrt(a)
         csa=datbas(igc+2,i)
c    cpa is needed only for the l type
         cpa=datbas(3,i)
         do 200 l=1,3
         xa(l)=datbas(l+10,i)
  200    continue
c-------
         do 500 j=ja,je
            b=datbas(1,j)
            sqb=sqrt(b)
            csb=datbas(jgc+2,j)
            cpb=datbas(3,j)
c-------
            apb=a+b
            r=zero
            e=a*b/apb
            do 300 l=1,3
               xb(l)=datbas(l+10,j)
               r=r+(xa(l)-xb(l))**2
  300       continue
c-----
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
c-----
            if(abs(s0).lt.thre1 ) go to 500
c-----------------------------------------------------------------------
       if(key.eq.1.or.key.eq.2.or.key.eq.3) then
         call incsh1(ityp1,jtyp1,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,rb,lxin)
       endif
c-----------------------------------------------------------------------
c      if(key.eq.4) then
c         call incsh2(ityp1,jtyp1,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,rb,
c    *               xra,lxin)
c      endif
c-----------------------------------------------------------------------
       if(key.eq.4) then
         call incsh2(ityp1,jtyp1,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,rb,
     *               xra,lxin,
     *               bl(irp01),bl(irp10),bl(irpmm),bl(irpmp),bl(irppm),
     *               bl(irpp0),bl(irm01),bl(irm10),bl(ixin0),
     *               bl(i_mm),bl(i_pmm),bl(i_mp),bl(i_pmp),bl(i_ppm) )
       endif
c-----------------------------------------------------------------------
       if(key.eq.5) then
         call incsh3(ityp1,jtyp1,key,k1,k2,a,b,sqa,sqb,s0,xa,xb,rb,
     *               xra,lxin)
       endif
c-----------------------------------------------------------------------
c
          ij=0
          do 410 ld=1,ldim
            do 400 i1=1,ilen1
               coefi=csa
               if (ityp.eq.3.and.i1.gt.1) coefi=cpa
            do 400 j1=1,jlen1
               coefj=csb
               if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
               ij=ij+1
               sb(ij)=sb(ij)+rb(ij)*coefi*coefj
  400          continue
  410       continue
  500    continue
  600 continue
c-----------------------------------------------------------------------
c release memory allocated for key=4
c
      if(key.eq.4) call retmem(14)
c-----------------------------------------------------------------------
c transformation d6->d5, f10->f7
c
      incre=len1
      do 710 ld=1,ldim
      iadd=1+(ld-1)*incre
c
        if(ityp.eq.4) then
           call dtran1a(rb,sb(iadd),jlen1)
        end if
        if(ityp.eq.6) then
           call ftran1a(rb,sb(iadd),jlen1)
        endif
        if(ityp.eq.11) then
           call gtran1a(rb,sb(iadd),jlen1)
        end if
        if(ityp.eq.12) then
           call htran1a(rb,sb(iadd),jlen1)
        end if
        if(ityp.eq.13) then
c           call itran1a(rb,sb(iadd),jlen1)
        end if
c
c at this point i-functions are transformed so the length is ilen
c
        if(jtyp.eq.4) then
           call dtran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.6) then
           call ftran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.11) then
           call gtran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.12) then
           call htran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.11) then
c           call itran2a(rb,sb(iadd),ilen)
        end if
  710 continue
c
c end of transformation
c
c move to appropriate location
c
      incre=len1
      ij=0
      do 720 ld=1,ldim
      iadd=(ld-1)*incre
      do 730 ij1=1,len
      ij=ij+1
      sb(ij)=sb(iadd+ij1)
  730 continue
  720 continue
c-----------------------------
      end
c====================================================================
      subroutine get_refat(datnuc,natoms,nrat)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
c
c natoms - total number of atoms in molecule
c
c output : nrat - heaviest atom (reference for D10 accuracy)
c
      qq=0.d0
      do 100 nat=1,natoms
         q =datnuc(1,nat)
         if(q.gt.qq) then
            qq=q
            nrat=nat
         endif
  100 continue
      end
c====================================================================
      subroutine getCosmoData(na,bl,datnuc)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,nX,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      dimension bl(*)
      dimension datnuc(5,na)
c
c    check if the COSMO flag has been defined, otherwise
c    define it now.
c
      call tstival('cosmo',icosmo)
      If(icosmo.EQ.0) Then
        call setival('cosmo',0)
      Else
        call getival('cosmo',icosmo)
      EndIf
c-------------------------------------------------------------------
      IF(icosmo.ne.0)then
c
c   obtain COSMO parameters from depository
c
        call getival('ndum',ndum)
        natom = na-ndum                 ! no. of real atoms
        call getrval('c_fepsi',fepsi)
        call getival('c_nps',nps)       ! no. of surface segments
c-------------------------------------------------------------------
c allocate memory for copy of datnuc with cosmo data in
c 1. segments coordinates
c 2. segments charges rescaled by epsi
c
c "new" datnuc will be dimensioned as datnuc(5,na+nps)
c
        call getmem(5*(na+nps),inucosmo)
        call setival('inucosmo',inucosmo)
c-------------------------------------------------------------------
        last=last+5*(na*nps)
c-------------------------------------------------------------------
c   allocate memory for COSMO forces
c
        call mmark
        call getmem(3*natom,icoxyz)
        call getmem(3*nps,icosurf)
        call getmem(nps,iar)
        call getmem(nps,iqcos)
        call getint(nps,iiatsp)
        call getint(natom,icharge)
c-------------------------------------------------------------------
c
c  read arrays from COSMO data file
c
        call cosmo_rdata(bl(icoxyz),bl(icosurf),bl(iar),bl(iqcos),
     $                   bl(iiatsp),bl(icharge),natom,nps)
c
c-------------------------------------------------------------------
c Make copy of datnuc with cosmo in:
c
        call makeCosDat(datnuc,na,nps,bl(icosurf),bl(iqcos),fepsi,
     *                   bl(inucosmo))
c-------------------------------------------------------------------
        call retmark
      endif
c
c--------------------------------------------------------------------
      end
c====================================================================
      subroutine  makeCosDat(datnuc,na,nps,cosurf,qcos,fepsi,datnucos)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,na)
      dimension cosurf(3,nps),qcos(nps)
      dimension datnucos(5,na+nps)
c
c rescale cosmo segment charges by epsi
c
      call dscal(nps,fepsi,qcos,1)
c
      call zeroit(datnucos,5*(na+nps))
      do i=1,na
         datnucos(1,i)=datnuc(1,i)
         datnucos(2,i)=datnuc(2,i)
         datnucos(3,i)=datnuc(3,i)
         datnucos(4,i)=datnuc(4,i)
      enddo
c
      ip=na
      do i=1,nps
         ip=ip+1
         datnucos(1,ip)=qcos(i)
         datnucos(2,ip)=cosurf(1,i)
         datnucos(3,ip)=cosurf(2,i)
         datnucos(4,ip)=cosurf(3,i)
      enddo
c
      end
c====================================================================
