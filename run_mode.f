      subroutine run_mode(bl,inx,ncf,minpr,maxpr, xminpr_bl,xmaxpr_bl,
     *                    lsemi,scftype,ireturn)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical in_core,on_disk
      character*11 scftype
      common /memmax/ ispblx, maxme1,iforwhat
      common /outfile/ ioutput
      common /lindvec/ lind,idensp
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1b/ nbl2,nbloks
      common /memors/ nsym ,ijshp,isymm
      common /symmet/ rnsym
c
      common /runmode/ sparsitx,min_price,max_price
c----------------------------------------------------------------
c  input :
c
c  ireturn=1 if called for Schwarz integrals otherwise 0
c----------------------------------------------------------------
c data for stored integrals in-core
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c----------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension intnu(15)
      dimension xntnu(15)
      dimension xint_blk(15),cost_blk(15)
c
c integrals from blocks more expensive than ...
c----------------------------------------------------------------
c save it in a common
      min_price=minpr
      max_price=maxpr
c----------------------------------------------------------------
c estimate a sparsity of the integral's file :
c
      call getmem(ncf*(ncf+1)/2,iovlp)
      call sparsity(bl,inx,bl(lind),rnsym,bl(iovlp),sparsitx )
      call retmem(1)
c
      spar_sym=rnsym
      sparsitx=sparsitx*spar_sym
c
      if(sparsitx.gt.1.0d0) sparsitx=1.0d0
c----------------------------------------------------------
c       write(ioutput,*)
c    * ' nsym,rnsym=',nsym,rnsym,' sparsity=',sparsitx,
c    * ' spar_sym=',spar_sym
c----------------------------------------------------------
      if(ireturn.eq.1) return
c----------------------------------------------------------
c calculete and print integral's statistics :
c number of 2e-int. and  price disribution
c
      ncost_int=ncost
      ncost_blk=ncost + nbloks +1
c
      call statint(inx,bl(iisd),bl(jjsd),bl(ijbld),bl(npard),nbl2,
     *             bl(ncost_int),xntnu,intnu,minpr,maxpr,
     *             bl(ncost_blk),xint_blk,cost_blk,xminpr_bl,xmaxpr_bl,
     *             q_size,sparsitx)
c
c integral price is in bl(ncost) , block price is in bl(ncost+nbloks+1)
c Block's price is not used now .
c----------------------------------------------------------
c check for integral storage ONLY for (1) scf (2) cphf (3) FTC
c
      If(lsemi.NE.0. AND.
     *  (iforwhat.eq.1.or.iforwhat.eq.2.or.iforwhat.eq.11)) Then
      Else
         scftype='full-direct'
         incorex=0
         return
      EndIf
c----------------------------------------------------------
c check how many integrals can be stored in-core & on-disk
c
      call check_storage(q_size,xntnu(1),
     *                   xint_core,xint_disk,in_core,on_disk)
c
c output :
c xint_core & xint_disk - number of int. that can be stored
c                         in-core & on-disk (estimated)
c logical in_core,on_disk
c----------------------------------------------------------
c setup to zero the following addresses for stored integrals:
c
       incorex=0
       ilab=0
       jlab=0
       klab=0
       llab=0
c
       isiz=0
       jsiz=0
       ksiz=0
       lsiz=0
c
c----------------------------------------------------------
c
      scftype='full-direct'
      if(in_core .or. on_disk) scftype='semi-direct'
      if(scftype.eq.'full-direct') return
c----------------------------------------------------------
c allocate memory for two arrays :
c (1) number of integrals in each super-block
c (2) number of quartets in each super-block
c
c Both above are the theoretical numbers (without symm. or neglect)
c They are used only for stored integrals
c
      call getmem(nbloks,integbl)
      call getmem(nbloks,iqrtsbl)
c
      call setival('integbl',integbl)
      call setival('iqrtsbl',iqrtsbl)
c----------------------------------------------------------
c find out what the price limit for stored inregrals is :
c
c     call find_blk_price(xint_blk,cost_blk,xint_core,xint_disk,
c    *                xminpr_bl,xmaxpr_bl,xprice_core,xprice_disk)
c
      call find_int_price(xntnu,intnu,xint_core,xint_disk,
     *                minpr,maxpr,iprice_core,iprice_disk)
c
c output :
c iprice_core (integ. more expensive than iprice_core - in-core)
c iprice_disk (integ. more expensive than iprice_disk - on-disk)
c----------------------------------------------------------
c decide which integrals will be stored in-core or on-disk and
c which will be recalculated
c
      if(in_core) then
c        change  prices of stored integrals to negative values :
c Total available storage in double words ON ALL nodes
         call getival('incore',incore)
         call getival('nslv', nproc)
         if(nproc.eq.0) nproc=nproc+1
c Memory needed for in-core stored intgrals is
         xmem_core=dble((incore-4)*nproc)
         call set_int_price(nbl2,bl(npard),bl(nsupb),
     *                  bl(ijbld),bl(iisd),bl(jjsd),inx,
     *                  'in_core', xmem_core, iprice_core,
     *                  sparsitx,
     *                  bl(ncost_int),bl(integbl),bl(iqrtsbl),
     *        xmem_2_store,xint_2_store,xqrt_2_store,nblk_2_store)
c
c        allocate memory to handle stored integrals :
         call set_coredat(xmem_2_store,xint_2_store,xqrt_2_store,
     *                    nblk_2_store,
     *                    incorex,ilab,jlab,klab,llab,
     *                            isiz,jsiz,ksiz,lsiz,iqrt)
         on_disk=.false.
      endif
c
      if(on_disk) then
         call getival('indisk',indisk) ! in MB
         xmem_disk=indisk
         xmem_disk=xmem_disk*1000000/8
         call set_int_price(nbl2,bl(npard),bl(nsupb),
     *                  bl(ijbld),bl(iisd),bl(jjsd),inx,
     *                  'on_disk', xmem_disk, iprice_disk,
     *                  sparsitx,
     *                  bl(ncost_int),bl(integbl),bl(iqrtsbl),
     *        xmem_2_store,xint_2_store,xqrt_2_store,nblk_2_store)
c
         intmem=int( 0.9d0*dble(indisk) ) 
         call setival('intmem',intmem)    ! in MWords
      endif
c
c output from the set_price() call  : for each super-block :
c  ncost(ikbl) is now : negative for stored integrals
c                       possitive for recalculated integrlas
c                       (absolute value is the original price)
c----------------------------------------------------------
cc      write(ioutput,*)' nsym,rnsym=',nsym,rnsym,' sparsity=',sparsitx
c----------------------------------------------------------
      end
c================================================================
      subroutine statint(inx,iis,jjs,ijbl,npar,nbl2,
     *             ncost_int,xntnu   ,intnu   ,minpr    ,maxpr,
     *             xcost_blk,xint_blk,cost_blk,xminpr_bl,xmaxpr_bl,
     *             q_size, sparsity)
      implicit real*8 (a-h,o-z)
      common /memmax/ ispblx, maxme1,iforwhat
      common /outfile/ ioutput
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*), npar(*)
      dimension ncost_int(*)
      dimension xcost_blk(*)
      dimension intnu(15)
      dimension xntnu(15)
      dimension xint_blk(15), cost_blk(15)
      data half /0.5d0/
c-------------------------------------------
c number of integrals more expensive than :
c    0.0%  of maximum integral-price (all integrals)
c    0.1%  of maximum integral-price
c    0.3%  of maximum integral-price
c    0.5%  of maximum integral-price
c    0.7%  of maximum integral-price
c    1.0%  of maximum integral-price
c    1.5%  of maximum integral-price
c    3.0%  of maximum integral-price
c    5.0%  of maximum integral-price
c    7.5%  of maximum integral-price
c    10.%  of maximum integral-price
c    15.%  of maximum integral-price
c    25.%  of maximum integral-price
c    50.%  of maximum integral-price
c    75.%  of maximum integral-price
c
c these prices in intnu() and number of integ. in xntnu()
c--------------------------
c number of integrals FROM big-blocks more expensive than
c    0.0%  of maximum block-price (all integrals)
c    0.1%  of maximum block-price
c    0.3%  of maximum block-price
c    0.5%  of maximum block-price
c    0.7%  of maximum block-price
c    1.0%  of maximum block-price
c    1.5%  of maximum block-price
c    3.0%  of maximum block-price
c    5.0%  of maximum block-price
c    7.5%  of maximum block-price
c    10.%  of maximum block-price
c    15.%  of maximum block-price
c    25.%  of maximum block-price
c    50.%  of maximum block-price
c    75.%  of maximum block-price
c
c these prices in cost_blk() and number of integ. in xint_blk()
c--------------------------
      do 15 ii=1,15
      xntnu(ii)=0.d0
      xint_blk(ii)=0.d0
   15 continue
c
      xprice=dble(maxpr)
      x000=0.d0
      x001=0.001d0*xprice
      x003=0.003d0*xprice
      x005=0.005d0*xprice
      x007=0.007d0*xprice
      x010=0.010d0*xprice
      x015=0.015d0*xprice
      x030=0.030d0*xprice
      x050=0.050d0*xprice
      x075=0.075d0*xprice
      x100=0.10d0*xprice
      x150=0.15d0*xprice
      x250=0.25d0*xprice
      x500=0.50d0*xprice
      x750=0.75d0*xprice
c
      intnu(1)=int(x000)
      intnu(2)=int(x001)
      intnu(3)=int(x003)
      intnu(4)=int(x005)
      intnu(5)=int(x007)
      intnu(6)=int(x010)
      intnu(7)=int(x015)
      intnu(8)=int(x030)
      intnu(9)=int(x050)
      intnu(10)=int(x075)
      intnu(11)=int(x100)
      intnu(12)=int(x150)
      intnu(13)=int(x250)
      intnu(14)=int(x500)
      intnu(15)=int(x750)
c-------------------------
      xprice=xmaxpr_bl
      x000=0.d0
      x001=0.001d0*xprice
      x003=0.003d0*xprice
      x005=0.005d0*xprice
      x007=0.007d0*xprice
      x010=0.010d0*xprice
      x015=0.015d0*xprice
      x030=0.030d0*xprice
      x050=0.050d0*xprice
      x075=0.075d0*xprice
      x100=0.10d0*xprice
      x150=0.15d0*xprice
      x250=0.25d0*xprice
      x500=0.50d0*xprice
      x750=0.75d0*xprice
c
      cost_blk(1)=x000
      cost_blk(2)=x001
      cost_blk(3)=x003
      cost_blk(4)=x005
      cost_blk(5)=x007
      cost_blk(6)=x010
      cost_blk(7)=x015
      cost_blk(8)=x030
      cost_blk(9)=x050
      cost_blk(10)=x075
      cost_blk(11)=x100
      cost_blk(12)=x150
      cost_blk(13)=x250
      cost_blk(14)=x500
      cost_blk(15)=x750
c-------------------------
c for average size of a quartet up to (ll|ss) & (ll|ll)
c (needed to estimate number of integrals that can be hold in-core)
c
      xint_part1=0.d0
      xqrt_part1=0.d0
      xint_part2=0.d0
      xqrt_part2=0.d0
c
      xint_all=0.d0
      xqrt_all=0.d0
c-------------------------
c contracted integrals counter
c
      ikbl=0
      do 100 ibl=1,nbl2
      nparij=npar(ibl)
            ijcs1=ijbl(ibl,1)
            ics1=iis(ijcs1)
            jcs1=jjs(ijcs1)
c
            ijs1=inx(3,ics1)*inx(3,jcs1)
            ijg1=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
c-
ccc   call parinteg(ibl,nparij,ijbl,nbl2,iis,jjs,ijnte,inte,ijxint)
c-
         do 100 kbl=1,ibl
         ikbl=ikbl+1
         nparkl=npar(kbl)
            klcs1=ijbl(kbl,1)
            kcs1=iis(klcs1)
            lcs1=jjs(klcs1)
c
            kls1=inx(3,kcs1)*inx(3,lcs1)
            klg1=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
c-----------------
            ijkls1=ijs1*kls1
            ijklg1=ijg1*klg1
c-----------------
            if(ibl.eq.kbl) then
               ijklq=nparij*(nparij+1)/2
            else
               ijklq=nparij*nparkl
            endif
c
c  quartets & integrals in one super-block :
c
            ijklq=ijklq*ijklg1
c
            qijkl=dble(ijklq)
c
c include sparsity
c
            qijkl=qijkl*sparsity
c
            xinteg1=qijkl*dble(ijkls1)
c
c  total number of quartets & integrals:
c
            xint_all=xint_all + xinteg1
            xqrt_all=xqrt_all + qijkl
c-----------------
         nprice=ncost_int(ikbl)
         xprice=xcost_blk(ikbl)
c-----------------
c calculate average size of a quartet (for in-core storage)
c (average over few first blocks)
c
         if(ijkls1.le. 16) then   ! only up to (ll|ss)
           xint_part1=xint_part1 + xinteg1
           xqrt_part1=xqrt_part1 + qijkl
         endif
         if(ijkls1.le.256) then   ! only up to (ll|ll)
           xint_part2=xint_part2 + xinteg1
           xqrt_part2=xqrt_part2 + qijkl
         endif
c-----------------
c
c count number of integrals more expensive than ..
c
         do ii=1,15
          if(nprice.gt.intnu(ii)   ) xntnu(ii)=xntnu(ii)+xinteg1
          if(xprice.gt.cost_blk(ii)) xint_blk(ii)=xint_blk(ii)+xinteg1
         enddo
c
  100 continue
c-----------------
c average size of a quartet calculated over blocks up to (ll|ls)
c
      q_size_part1=xint_part1/xqrt_part1
      q_size_part2=xint_part2/xqrt_part2
      q_size_full=xint_all /xqrt_all
c
      q_size=(q_size_part1 + q_size_part2)*0.5d0
c-----------------
c     write(ioutput,*)' total number of quartets =',xqrt_all
c     write(ioutput,*)' total number of integrals=',xint_all
c     write(ioutput,*)' avergae size of quartet : '
c     write(ioutput,*)'  overall=',int(q_size_full),
c    *                   ' part1=',int(q_size_part1),
c    *                   ' part2=',int(q_size_part2),
c    *                   '  used=',int(q_size)
c-----------------
c
      call statint1(xntnu,intnu,minpr,maxpr,
     *              xint_blk,cost_blk,xminpr_bl,xmaxpr_bl,
     *              sparsity)
c
      end
c================================================
      subroutine parinteg(ibl,nparij,ijbl,nbl2,iis,jjs,
     *                    ijnt,int, ijxint)
      implicit real*8 (a-h,o-z)
      dimension iis(*),jjs(*),ijbl(nbl2,*)
c-
       ijeq=0
       do 1001 ijpar=1,nparij
        ijcs=ijbl(ibl,ijpar)
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        if(jcs.eq.ics) ijeq=ijeq+1
 1001  continue
       ijxint=(nparij-ijeq)*ijnt + ijeq*int*(int+1)/2
c-
      end
c================================================
      subroutine statint1(xntnu,intnu,minpr,maxpr,
     *                    xint_blk,cost_blk,xminpr_bl,xmaxpr_bl,
     *                    sparsity)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      dimension xntnu(15), intnu(15)
      dimension xint_blk(15), cost_blk(15)
c-----------------------------------
c print statistics : number of 2e-int. and its price disribution
c
c--
ccc   xinteg=xntnu(1)*sparsity
      xinteg=xntnu(1)            !    sparsity already included
      x1024=1024.d0
      xmega=1.d0/(x1024*x1024)
      if(xinteg.gt.0.d0) xint1=1.d0/xinteg
c--
      write(ioutput,500)
  500 format(72('.'))
      write(ioutput,501)
c 501 format(' Integral statistics (w/o symmetry and neglect)'/)
c     write(ioutput,502) xntnu(1) ,minpr,maxpr,xminpr_bl,xmaxpr_bl
  501 format(' Integral statistics                           '/)
      write(ioutput,502) xinteg   ,minpr,maxpr,xminpr_bl,xmaxpr_bl
  502 format('  total number of integrals to calc.=',f18.0/
     *   '  computational cost/integral : min.=',i3,' , max=',i6/
     *   '  computational cost/big-block: min.=',f15.0,' max=',f15.0)
      write(ioutput,500)
c...........................
c     write(ioutput,503) sparsity
c 503 format(' disk storage for integral file with sparsity=',f5.3/)
c     part=0.d0
c
c     do 20 ii=1,15
c       if(xntnu(ii).gt.0.d0) then
c          xintx=xntnu(ii)*sparsity
c          xbyte=xintx*8.d0*xmega
c          if(xinteg.gt.0.d0) part=xintx*xint1
c          part=part*100.d0
c          write(ioutput,504) xintx,part,intnu(ii),xbyte
c include sparsity :
ccc        xntnu(ii)=xintx
c       endif
c  20 continue
c
c     write(ioutput,500)
c 504 format
c    *(f17.0,' (',f5.1,'%) int. stored for ipay=',i6,' take ',
c    *                                               f10.2,'MB')
c...........................
c---
      end
c================================================================
      subroutine sparsity(bl,inx,lind,rnsym,overlap,sparx)
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /outfile/ ioutput
      dimension bl(*)
      dimension inx(12,*)
      dimension lind(*)
      dimension overlap(*)
c---------------------------
      thres=sqrt(eps)
c---------------------------
ccccc call inton (0,na,overlap,inx,0,0,bl(ibas),bl(inuc),ncs)
ccccc call over_prt(overlap,ncf,'full')
      call intons(0,na,overlap,inx,0,0,bl(ibas),bl(inuc),ncs)
ccccc call over_prt(overlap,ncf,'ss00')
c---------------------------
      xint0=0.d0
      xint1=0.d0
         do 150 icf=1,ncf
         do 150 jcf=1,icf
         ij=lind(icf)+jcf
         sij=overlap(ij)
         xint0=xint0+1.d0
         if(abs(sij).ge.thres) xint1=xint1+1.d0
  150    continue
c
c---------------------------
c
      spar1=xint1/xint0
      const=1.8d0              !  assumption
      spar2=const*spar1*spar1
      sparx=spar2
cccc  sparx=rnsym*spar2
      if(sparx.gt.1.d0) sparx=1.d0
c
c     write(ioutput,800) rnsym,xint0,xint1,xint1/xint0,sparx
c 800 format('rnsym=',f8.5,' all=',f10.0,' non-0=',f10.0/
c    *       'n-0/all=',f10.5,' sparsity=',f10.5)
      write(ioutput,800) sparx
  800 format(/' estimated sparsity of integrals =',f5.2/)
c
      end
c==========================================================
      subroutine over_prt(over,ncf,text)
      implicit real*8 (a-h,o-z)
      character*4 text
      dimension over(*)
c
      write(8,*) 'overlap ',text
c
      do 200 i=1,ncf
      ii=i*(i-1)/2
      write(8,88)(over(ii+j),j=1,i)
  200 continue
   88 format(13(f9.4,1x))
      end
c==========================================================
      subroutine intons (itype,nna,oneint,inx,kk1,kk2,basdat,datnuc,ncs)
      implicit real*8 (a-h,o-z)
c
c     itype (in) is 0 for overlap,
c
      dimension inx(12,*), s(784)
      dimension oneint(*),basdat(13,*),datnuc(5,*)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      k1=kk1
      k2=kk2
c------------------------------------------------------------------
      key=itype
      if (itype.le.1) key=itype+1
c     total number of 1-el. integrals
      iza=0
      ifu=0
c     cycle through the contracted shells and calculate the integrals
      do 60 i=1,ncs
c       number of gen. contractions
        ngci=inx(4,i)
        iat=inx(2,i)
        do 55 igc=0,ngci
         jfu=0
         len1=inx(3,i)
         do 50 j=1,i
           ngcj=inx(4,j)
           jat=inx(2,j)
           ngcjx=ngcj
           if(j.eq.i) ngcjx=igc
           do 45 jgc=0,ngcjx
            len2=inx(3,j)
            len=len1*len2
            call onels(i,j,igc,jgc,.true.,basdat,basdat,datnuc,s,
     1      inx,key,nna,k1,k2)
            iij=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
               ii=iff*(iff-1)/2
            do 40 j1=1,len2
               jff=jff+1
               iij=iij+1
               if (jff.le.iff) then
                   ij=ii+jff
                   oneint(ij)=s(iij)
                   iza=iza+1
               end if
   40       continue
            jfu=jfu+len2
c        end gen. contr. loop j
   45     continue
   50    continue
         ifu=ifu+len1
c        end gen. contr. loop i
   55  continue
   60 continue
ctest
cc      if(iza.gt.0) then
cc        if (iprint.eq.1) then
cc          write (6,90) (oneint(i),i=1,iza)
cc   90     format (2x,5f14.9)
cc        end if
cc      end if
c
      end
      subroutine onels(ics,jcs,igc,jgc,same,basdat,basdat2,datnuc,s,
     1  inx,key,nn,k1,k2)
      implicit real*8 (a-h,o-z)
      logical same
      dimension basdat(13,*),basdat2(13,*),datnuc(5,*)
      dimension s(784), xa(3), xb(3), inx(12,1), xint(784), p(3)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      data sqtw/0.2886751345905d0/
c
c this is 1/sqrt(12)
c
      data twopi/0.6366197723675d0/
c
c key=1 for overlap,
c
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
      ilen=inx(3,ics)
      jlen=inx(3,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
c
      if(ityp.eq.6) ilen1=10
      if(jtyp.eq.6) jlen1=10
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
c
c for d functions, 6 spaces are needed prior to the transformation
c
      do 10 i=1,len1
   10 s(i)=zero
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      ja=inx(1,jcs)+1
      je=inx(5,jcs)
      do 60 i=ia,ie
         a=basdat(1,i)
         csa=basdat(igc+2,i)
c    cpa is needed only for the l type
         cpa=basdat(3,i)
         sqa=sqrt(a)
         do 20 l=1,3
           xa(l)=basdat(l+10,i)
   20    continue
         do 50 j=ja,je
            if (same) then
              b=basdat(1,j)
              csb=basdat(jgc+2,j)
              cpb=basdat(3,j)
              do 24 l=1,3
               xb(l)=basdat(l+10,j)
  24          continue
            else
              b=basdat2(1,j)
              csb=basdat2(jgc+2,j)
              cpb=basdat2(3,j)
              do 26 l=1,3
               xb(l)=basdat2(l+10,j)
  26          continue
            end if
            sqb=sqrt(b)
            apb=a+b
            r=zero
            e=a*b/apb
            do 30 l=1,3
               p(l)=(a*xa(l)+b*xb(l))/apb
               r=r+(xa(l)-xb(l))**2
   30       continue
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
c
c     write(iout,1000) ics,jcs,ityp,jtyp,i,j,ilen1,jlen1,csa,cpa,csb,cpb
c
      ityp1=ityp
      jtyp1=jtyp
c     correspondance between ityp and ityp1:
c     ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28)
c     ityp1 1    2    3    4    4     5    5      6      7       8
c     ityp 11(g)  12(h)  13(i)
c     ityp1 6      7      8
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
c
c     *** no special type is needed for d6 at this point
c
            ij=0
            do 40 i1=1,ilen1
               coefi=csa
               if (ityp.eq.3.and.i1.gt.1) coefi=cpa
            do 40 j1=1,jlen1
               coefj=csb
               if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
               ij=ij+1
               s(ij)=s(ij)+    s0  *coefi*coefj
   40       continue
   50    continue
   60 continue
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
c     transform the d functions
        if (ityp.eq.4) then
           call dtran1a(xint,s,jlen1)
        end if
c     transform the f functions
        if (ityp.eq.6) then
           call ftran1a(xint,s,jlen1)
        end if
c     transform the g functions
        if(ityp.eq.11) then
           call gtran1a(xint,s,jlen1)
        end if
        if(ityp.eq.12) then
           call htran1a(xint,s,jlen1)
        end if
        if(ityp.eq.13) then
c           call itran1a(xint,s,jlen1)
        end if
c
c at this point the first (i) functions have been transformed so the length is ilen
c
c     transform the d functions
        if(jtyp.eq.4) then
           call dtran2a(xint,s,ilen)
        end if
c transform the f functions to 7 independent ones
        if(jtyp.eq.6) then
          call ftran2a (xint,s,ilen)
        end if
c  transform the g functions
        if(jtyp.eq.11) then
          call gtran2a(xint,s,ilen)
        end if
c  transform the h functions
        if(jtyp.eq.12) then
          call htran2a(xint,s,ilen)
        end if
c  transform the g functions
        if(jtyp.eq.13) then
c          call itran2a(xint,s,ilen)
        end if
czzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
c
      end
c================================================================
      subroutine check_storage(q_size,x100,
     *                         xint_core,xint_disk,in_core,on_disk)
      implicit real*8 (a-h,o-z)
      logical in_core,on_disk
      common /outfile/ ioutput
c---------------------------------------------------------------
c input :
c q_size - estimated size of a quartet (from first 15% of blocks)
c x100 - total number of integrals (with sparsity ?)
c
c output :
c xint_core - number of integrals that can be stored in-core
c xint_disk - number of integrals that can be stored on-disk
c
c xint_core & xint_disk might be 0 if it is below 10% .
c---------------------------------------------------------------
      i41=ioutput
c---------------------------------------------------------------
c get storage info : available storage in double words :
c
c  This is on each SINGLE node
c
      call getival('incore',incore)
      call getival('indisk',indisk)  ! in MW(ords)
c
      xdiske=dble(indisk)
      xncore=dble(incore)
c---------------------------------------------------------------
c At least 2M must be given for in-core storage & 100M for on-disk
c
c  incore is given in  Words
c  indisk is given in MBytes
c
      xint_core=0.d0
      xint_disk=0.d0
c
      on_disk=.true.
      in_core=.true.
      if(incore.lt. 2000000) then
c        write(6,*)' In-core storage requires at least 2MW of memory'
         in_core=.false.
      endif
      if(indisk.lt.     100) then
c        write(6,*)' On-disk storage requires at least 100MB of disk'
         on_disk=.false.
      endif
c
      if( (.not.in_core) .and. (.not.on_disk) ) RETURN
c---------------------------------------------------------------
c Total available storage in double words ON ALL nodes
c
      call getival('nslv', nproc)
c
      if(nproc.eq.0) nproc=nproc+1
c
      xncore=xncore*nproc
      xdiske=xdiske*nproc
c---------------------------------------------------------------
c Memory needed for in-core stored intgrals is
c (1) for integrals themsleves ( in R*8 )
c (2) for label's info: 4 starting labels ( in I*2 )
c
c  Note : 4 shell-sizes are stored for each superblock
c
c Thus, for xint to store we need : xint + xqrt*4/4 double words
c
c Assuming average quartet size q_size (integrals/quartet)
c gives :
c
c M = xint + xint/q_size = xint*[1 + 1/q_size]
c
c Hence
c
c  xint = M/(1 + 1/q_size)
c
c Integrals from low ang.mom. quartets will be more likely stored
c which means short quartets (ssss),(psss),(ppss) etc.
c The average  q_size is etimated in statint().
c
c It can also be determined upon number of quartets :
c
c M = xint + quart = quart*q_size + quart
c
c M = quart*(q_size +1)
c
c Number of quartest that can be stored in M is :
c
c quart=M/(q_size + 1)
c
c Then for integrals we have   M - quart
c---------------------------------------------------------------
c In-core storage : integrals R*8 and eight integers I*2 for label info
c
ccc   factor=1.d0/(q_size + 1.d0)           ! for quartets
c
      factor=q_size/(q_size + 1.d0)  ! for integrals
c
      xint_core=xncore*factor
c
c On-disk storage
c
      xint_disk=xdiske*1 000 000
c
c---------------------------------------------------------------
c calculate percentage of integrals that can be stored
c
c     if(xint_disk.lt.x100) then
c        x_on_disk=xint_disk/x100
c     else
c        x_on_disk=1.d0
c     endif
c
c     if(xint_core.lt.x100) then
c        x_in_core=xint_core/x100
c     else
c        x_in_core=1.d0
c     endif
c
c     x_in_core=x_in_core*100.d0
c     x_on_disk=x_on_disk*100.d0
c
c---------------------------------------------------------------
C On Jon's request 10% condition is not in effect anymore
c
c if it is not possible to store at least 10% then do not bother
c with storage (full-direct mode)
c
c     x010=0.10d0*x100      ! 10% of all integrals
c
c     on_disk=.true.
c     in_core=.true.
c     if(xint_disk.lt.x010) on_disk=.false.
c     if(xint_core.lt.x010) in_core=.false.
c---------------------------------------------------------------
c     write(i41,100) 'in core',xncore,xint_core,x_in_core
c     write(i41,100) 'on disk',xdiske,xint_disk,x_on_disk
c
c 100 format(' number of integrals that can be stored ',a8,
c    * '(', f15.0 ,') is ',f15.0,' (',f5.1,' %)' )
c---------------------------------------------------------------
      end
c==============================================================
      subroutine find_int_price(xntnu,intnu,xint_core,xint_disk,
     *                          minpr,maxpr, iprice_core,iprice_disk)
c---------------------------------------------------------------
c input:
c xntnu(15) -number of int. more expensive than
c intnu(15)   corresponding prices
c
c xint_core - number of integrals that can be stored in-core
c xint_disk - number of integrals that can be stored on-disk
c
c minpr - minimum price (per integral)
c maxpr - maximum price (per integral)
c
c output :
c
c iprice_core - integrals more expensive than IPRICE_CORE
c               will be stored in-core
c iprice_disk - integrals more expensive than IPRICE_DISK
c               will be stored on disk
c
c NOTE : only in-core OR on-disk storage (not both)
c
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      dimension xntnu(15)
      dimension intnu(15)
c
c--
      do ii=1,15
         xinteg=xntnu(ii)
         xprice=dble(intnu(ii))
         write(ioutput,5021) xprice,xinteg
      enddo
 5021 format('  number of integ. more expensive than ',f8.0,1x,f15.0)
c--
c
c for in-core storage
c
      iprice_core=maxpr
      iicore=0
      do ii=1,15
         xinteg=xntnu(ii)
         if(xinteg.lt.xint_core) then
            iprice_core=intnu(ii)
            iicore=ii
            go to 100
         endif
      enddo
 100  continue
c
      if(iicore.gt.1) iprice_core=intnu(iicore-1)
c
c independently for disk-storage
c
      iprice_disk=maxpr
      iidisk=0
      do ii=1,15
         xinteg=xntnu(ii)
         if(xinteg.lt.xint_disk) then
            iprice_disk=intnu(ii)
            iidisk=ii
            go to 200
         endif
      enddo
 200  continue
c
      write(ioutput,*)'------------------------------------------------'
      write(ioutput,*)'store in-core more expensive than :',iprice_core
      write(ioutput,*)'store on-disk more expensive than :',iprice_disk
      write(ioutput,*)'------------------------------------------------'
c
      end
c
c==============================================================
      subroutine find_blk_price(xint_blk,cost_blk,xint_core,xint_disk,
     *                 xminpr_bl,xmaxpr_bl,xprice_core,xprice_disk)
c---------------------------------------------------------------
c input:
c xint_blk(15) -number of int. more expensive than
c cost_blk(15)   corresponding prices
c
c xint_core - number of integrals that can be stored in-core
c xint_disk - number of integrals that can be stored on-disk
c
c xminpr_bl - minimum price (per block )
c xmaxpr_bl - maximum price (per block )
c
c output :
c
c xprice_core - integrals from blocks more expensive than XPRICE_CORE
c               will be stored in-core
c xprice_disk - integrals from blocks more expensive than XPRICE_DISK
c               will be stored on disk
c
c NOTE : only in-core OR on-disk storage (not both)
c
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      dimension xint_blk(15)
      dimension cost_blk(15)
c
c--
      do ii=1,15
         xinteg=xint_blk(ii)
         xprice=cost_blk(ii)
         write(ioutput,5021) xprice,xinteg
      enddo
 5021 format('  number of integ. from blocks more expensive than ',
     *       f15.0,1x,f15.0)
c--
c
c for in-core storage
c
      xprice_core=xmaxpr_bl
      do ii=1,15
         xinteg=xint_blk(ii)
         if(xinteg.lt.xint_core) then
            xprice_core=cost_blk(ii)
            go to 100
         endif
      enddo
 100  continue
c
c independently for disk-storage
c
      xprice_disk=xmaxpr_bl
      do ii=1,15
         xinteg=xint_blk(ii)
         if(xinteg.lt.xint_disk) then
            xprice_disk=cost_blk(ii)
            go to 200
         endif
      enddo
 200  continue
c
      write(ioutput,*)'------------------------------------------------'
      write(ioutput,*)'store in-core more expensive than :',xprice_core
      write(ioutput,*)'store on-disk more expensive than :',xprice_disk
      write(ioutput,*)'------------------------------------------------'
c
      end
c
c==============================================================
      subroutine set_coredat(xmem_in_core,xint_in_core,xqrt_in_core,
     *                       nblk_in_core,
     *                       incorex,ilab,jlab,klab,llab,
     *                               isiz,jsiz,ksiz,lsiz,iqrt)

      use memory

      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
c-------------------------------------------------------------
c input :
c xmem_in_core - memory needed for in-core storage (as estimated)
c xint_in_core - number of integrals to be stored in-core
c xqrt_in_core - number of quartests to be stored in-core
c nblk_in_core - number of super-blocks to be stored in-core
c
c NOTE : all above concern ALL NODES - on each single node
c                                      1/nproc * above
c
c output :
c
c incorex - address in bl for the core() array
c ilab-llab - addresses of corresponding label arrays
c isiz-lsiz - addresses of corresponding size's arrays
c-------------------------------------------------------------
      i41=ioutput
c-------------------------------------------------------------
      call getival('incore',incore)      ! memmory
      call getival('incorex',incore_add)    ! address in bl
c-------------------------------------------------------------
c Number of processors :
c
      call getival('nslv', nproc)
      if(nproc.eq.0) nproc=1
c
c
c precdicted memory, no of integ. & number of quartets on ONE CPU:
c
ccc   xmem_in_cor1=xmem_in_core/nproc
c
      xint_in_cor1=xint_in_core/nproc+nproc+1
      xqrt_in_cor1=xqrt_in_core/nproc+nproc+1
c-------------------------------------------------------------
c     write(6,*)'in-core=', xint_in_core,xqrt_in_core,nblk_in_core
c-------------------------------------------------------------
c Availabe in-core space split in two parts :
c (1) for label's info : address = incore_add
c (2) for integrals themselvs address=incore_add + quartet number
c
c Addressing for label's info (kept as I*2):
c
c ilab(nqrt)
c jlab(nqrt)
c klab(nqrt)
c llab(nqrt)
c
      ispac4= int(xqrt_in_cor1)/4 +1
c
      ilab=incore_add
      jlab=ilab + ispac4
      klab=jlab + ispac4
      llab=klab + ispac4
c
c      write(6,*)' memory: labels=',4*ispac4
c
c integral's address
c
      incorex=llab + ispac4
c
c............................
c put sizes out of in-core area
c
c Allocate memory for shell's sizes and number of quartets
c for each super-block that will be stored in-core
c ( sizes I*2 , number of qrt in I*4 )
c
      isize=nblk_in_core
      nqrts=nblk_in_core
c
c      write(6,*)' memory: sizse =',4*isize,' nqrts=',nqrts
c
      call getint_2(isize,isiz)
      call getint_2(isize,jsiz)
      call getint_2(isize,ksiz)
      call getint_2(isize,lsiz)
c
      call getint_4(nqrts,iqrt)
c-------------------------------------------------------------
c memory for labels
c
      labmem=4*ispac4+4
c
c memory for integrals
c
      intmem=incore - labmem
ctest
c        write(6,*)' prepint: incore=',incore,labmem
c        write(6,*)' xqrt: xqrt_in_core=',xqrt_in_core
c        write(6,*)'       xqrt_in_cor1=',xqrt_in_cor1
c        write(6,*)'       ispac4      =',ispac4
c        if(intmem.le.0) stop
ctest
c
      call setival('intmem',intmem)
c-------------------------------------------------------------
c actual number of integrals to store in-core :
c
c On one node
c
      integ1=int(xint_in_cor1)
      nqrte1=int(xqrt_in_cor1)
c
      memor1= integ1 + 4*ispac4+4
c
c On all nodes
c
      integx=integ1*nproc
      nqrtex=nqrte1*nproc
      nmemory=memor1*nproc
c
      mem_in_core=int(xmem_in_core)
c     if(nmemory.ge.incore*nproc-nproc) then
c        write(6,*)' not enough in-core memory to store '
c        write(6,*) integx,' integrals from ',nqrtex ,' quartets'
c        write(6,*) ' Memory avail.=',incore*nproc,'; needed=',nmemory
c        write(6,*) ' nproc',nproc,'; incore=',incore,')'
c        stop
c     endif
c
c     write(i41,*)'Memory for in-core needed on all nodes =',mem_in_core
c     write(i41,*)'   expected number of stored integrals =',integx
c     write(i41,*)'   expected number of stored quartets  =',nqrtex
c
      end
c==============================================================
      subroutine set_int_price(nbl2,npar,nsupb,ijbl,iis,jjs,inx,
     *                     storage,xmem_store, iprice_store,
     *                     sparsity,
     *                     ncost,blinteg,bliqrts,
     *      xmem_2_store,xint_2_store,xqrt_2_store,nblk_2_store)
      implicit real*8 (a-h,o-z)
      character*7 storage        ! in_core or on_disk
c--------------------------------------------------------------
c this routine determines which integrals will be stored
c in-core or on-disk
c
c input : storage - 'in-core' or 'on-disk'
c         xmem_store - available memory to store integrals & label's info
c         iprice_store - price above which int. will be stored
c         sparsity - predicted sparsity of integral file
c
c output: ncost(iklb) - negative price for stored integrals
c                       possitive price for recalculated
c         blinteg(ikbl)- number of integrals in a block (w/o sparsity)
c         bliqrts(ikbl)- number of quartets  in a block (w/o sparsity)
c         xmem_2_store - total memory taht WILL be used for in-core storage
c         xint_2_store - number of integrals that WILL be stored
c         xqrt_2_store - number of integrals that WILL be stored
c         nblk_2_store - number of super blocks that WILL be stored
c--------------------------------------------------------------
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /outfile/ ioutput
c
      dimension npar(nbl2),ijbl(nbl2,*)
      dimension iis(*),jjs(*),inx(12,*)
      dimension ncost(*),nsupb(*)
      dimension blinteg(*), bliqrts(*)
c---------------------------------------------------------------
cNote : sparsity as it was estimated, is used for integrals only
c for quartets (for safety) 1.5* sparsity (but not >1)
c
      q_sparsity=sparsity*1.5d0
      if(q_sparsity.gt. 1.0d0) q_sparsity=1.0d0
      if(q_sparsity.lt. 0.5d0) q_sparsity=0.5d0
c
c However, if quatets are short (e.g. sto-3g basis set) then
c sparsity for quartets should be sharper
c
      q_length=dble(ncf)/dble(ncs)
      if(q_length .le. 2.5d0 .and. ncf.le.100)  q_sparsity=1.0d0
c---------------------------------------------------------------
c
c storage & recalculations
c
      int0=0        ! block counter
      int2=0        ! block counter
      xnt0=0.d0     ! re-calc.integral counter
      xnt2=0.d0     ! stored integral counter
      xtot=0.d0     ! total number of integ.
c
      qtot=0.d0     ! total number of quartets
      qnt0=0.d0     ! re-calc.quartets counter
      qnt2=0.d0     ! stored  quartets counter
c
      nblk=0.d0     ! total number of small blocks
      nsm0=0.d0     ! re-calc.small blocks
      nsm2=0.d0     ! stored  small blocks
c
      xmem_used=0.d0
c
      ikbl=0
      do ibl=1,nbl2
         ijpar=npar(ibl)
         ijcs1=ijbl(ibl,1)
         ics1=iis(ijcs1)
         jcs1=jjs(ijcs1)
         ijs1=inx(3,ics1)*inx(3,jcs1)
         ijg1=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
         do kbl=1,ibl
            klpar=npar(kbl)
            klcs1=ijbl(kbl,1)
            kcs1=iis(klcs1)
            lcs1=jjs(klcs1)
            kls1=inx(3,kcs1)*inx(3,lcs1)
            klg1=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
c
            ikbl=ikbl+1
            iprice=ncost(ikbl)
            ismall=nsupb(ikbl)
c
            if(ibl.eq.kbl) then
               ijklq=ijpar*(ijpar+1)/2
            else
               ijklq=ijpar*klpar
            endif
c
            ijklq=ijklq*ijg1*klg1
            q=dble(ijklq)
c
            xinteg=q*dble(ijs1*kls1)
c
c setup the blinteg(ikbl) and bliqrts(ikbl) arrays :
c
            blinteg(ikbl)=xinteg
            bliqrts(ikbl)=q
c
c include sparsity for integrals:
c
            x=xinteg*sparsity
c
c include q_sparsity for quartets :
c
            q=q*q_sparsity
c
            xtot=xtot+x
            qtot=qtot+q
c
            nblk=nblk+ismall
c memory
            xmem = x + q
c
            if(iprice.gt.iprice_store) then
               if(xmem_used+xmem.LT.xmem_store) then
                  xmem_used=xmem_used+xmem
                  xnt2=xnt2+x
                  int2=int2+1
                  ncost(ikbl)=-iprice
                  qnt2=qnt2+q
                  nsm2=nsm2+ismall
               else
                  xnt0=xnt0+x
                  int0=int0+1
cccccc            ncost(ikbl)= no change
                  qnt0=qnt0+q
                  nsm0=nsm0+ismall
               endif
            else
                  xnt0=xnt0+x
                  int0=int0+1
cccccc            ncost(ikbl)= no change
                  qnt0=qnt0+q
                  nsm0=nsm0+ismall
            endif
c
c..         write(ioutput,*)
c..  *     ' block=',ikbl,' prices=',iprice,ncost(ikbl)
c
         enddo
      enddo
c-------------------------------------------------
      xnt2_1=xnt2
      qnt2_1=qnt2
      xmem_1=xmem_used
c-------------------------------------------------
c try to put in-core more integrals from remaining blocks
c
      ikbl=0
      do ibl=1,nbl2
         ijpar=npar(ibl)
         ijcs1=ijbl(ibl,1)
         ics1=iis(ijcs1)
         jcs1=jjs(ijcs1)
         ijs1=inx(3,ics1)*inx(3,jcs1)
         ijg1=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
         do kbl=1,ibl
            klpar=npar(kbl)
            klcs1=ijbl(kbl,1)
            kcs1=iis(klcs1)
            lcs1=jjs(klcs1)
            kls1=inx(3,kcs1)*inx(3,lcs1)
            klg1=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
            ikbl=ikbl+1
            iprice=ncost(ikbl)
            if(iprice.gt.0) then
               ismall=nsupb(ikbl)
c
               if(ibl.eq.kbl) then
                  ijklq=ijpar*(ijpar+1)/2
               else
                  ijklq=ijpar*klpar
               endif
c
               ijklq=ijklq*ijg1*klg1
c
               q=dble(ijklq)
c
               xinteg=q*dble(ijs1*kls1)
c
c include sparsity for integrals :
c
               x=xinteg*sparsity
c
c include q_sparsity for quartets :
c
               q=q*q_sparsity
c
c memory
               xmem = x + q
c
               if(xmem_used+xmem.LT.xmem_store) then
                  xmem_used=xmem_used+xmem
                  xnt2=xnt2+x
                  int2=int2+1
                  ncost(ikbl)=-iprice
                  qnt2=qnt2+q
                  nsm2=nsm2+ismall
c
                  xnt0=xnt0-x
                  int0=int0-1
                  qnt0=qnt0-q
                  nsm0=nsm0-ismall
               endif
            endif
         enddo
      enddo
c-------------------------------------------------
      nbloks=nbl2*(nbl2+1)/2
c-------------------------------------------------
      xmem_2_store=xmem_used
      xint_2_store=xnt2
      xqrt_2_store=qnt2
ckw   nblk_2_store=int2
      int2=(int2*15)/10 ! increase by 1.5
      nblk_2_store=int2
      if(nblk_2_store.gt.nbloks) nblk_2_store=nbloks
c-------------------------------------------------
c     if(xnt2.gt.xnt2_1) then
c        write(ioutput,*)' more mem.in-core=',xmem_1,xmem_used
c        write(ioutput,*)' totalmem.in-core=',xmem_store
c        write(ioutput,*)' more int.in-core=',xnt2_1,xnt2
c        write(ioutput,*)' more qrt.in-core=',qnt2_1,qnt2
c     endif
c-------------------------------------------------
c
      write(ioutput,199)
  199 format(/'   Prediction for integral storage '/)
c     write(ioutput,200) nbloks,storage,  int2,int0
      write(ioutput,200) nbloks,storage,  nblk_2_store,int0
  200 format(
     *' Total number of big blocks =',i15/
     *'  ',a7, ' stored big blocks =',i15/
     *'    recalculated big blocks =',i15/)
c
      write(ioutput,210) nblk ,nsm2, nsm0
  210 format(
     *' Total number of small blks =',i15/
     *'          stored small blks =',i15/
     *'    recalculated small blks =',i15/)
c
      xtot1=1.d0/xtot
      x2p=xnt2*xtot1*100.d0
      x0p=xnt0*xtot1*100.d0
c
      write(ioutput,220) xtot,xnt2,xnt0
  220 format(
     *' Total number of integrals  =',f15.0/
     *'          stored integrals  =',f15.0/
     *'    recalculated integrals  =',f15.0/)
c
      write(ioutput,230) qtot,qnt2,qnt0
  230 format(
     *' Total number of quartests  =',f15.0/
     *'          stored quartests  =',f15.0/
     *'    recalculated quartests  =',f15.0/)
c
ctest
c       do ii=1,nbloks
c          xx=blinteg(ii)
c          qq=bliqrts(ii)
c          write(6,*)'block=',ii,' no int=',xx,' no qrts=',qq
c       enddo
ctest
c
      end
c=======================================================================
      subroutine update_run_mode(bl,inx,ncf,lsemi,scftype)

      use memory, block => bl

c
c called only for FTC to store some integrals in core
c
c
      implicit real*8 (a-h,o-z)
      logical in_core,on_disk
      character*11 scftype
      character*11 scftypex
      character*11 run_type
c
      common /runtype/ scftypex
c
      common /memmax/ ispblx, maxme1,iforwhat
      common /outfile/ ioutput
      common /lindvec/ lind,idensp
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1b/ nbl2,nbloks
      common /memors/ nsym ,ijshp,isymm
      common /symmet/ rnsym
c
      common /runmode/ sparsitx,minpr,maxpr
c----------------------------------------------------------------
c data for stored integrlas in-core
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c----------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension intnu(15)
      dimension xntnu(15)
      dimension xint_blk(15),cost_blk(15)
c----------------------------------------------------------------
c integrals from blocks more expensive than ...
c----------------------------------------------------------------
c     run_type=scftype
c----------------------------------------------------------------
c change back prices from negative to possitive 
c
      call change_price(nbloks,bl(ncost))
c----------------------------------------------------------
c calculete and print integral's statistics :
c number of 2e-int. and  price disribution
c
      ncost_int=ncost
c
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc1)
c
      call statint_ftc(inx,bl(iisd),bl(jjsd),bl(ijbld),bl(npard),nbl2,
     *                 bl(ncost_int),xntnu,intnu,minpr,maxpr,
     *                 q_size,sparsitx, bl(iftc1))
c
c integral price is in bl(ncost) 
c----------------------------------------------------------
c check for integral storage ONLY for (1) scf (2) cphf (3) FTC
c
      if(lsemi.ne.0 .and.
     $  (iforwhat.eq.1.or.iforwhat.eq.2.or.iforwhat.eq.11)) then
      else
         scftype ='full-direct'
         scftypex='full-direct'
         incorex=0
         return
      endif
c----------------------------------------------------------
c check how many integrals can be stored in-core & on-disk
c
      call check_storage(q_size,xntnu(1),
     *                   xint_core,xint_disk,in_core,on_disk)
c
c output :
c xint_core & xint_disk - number of int. that can be stored
c                         in-core & on-disk (estimated)
c logical in_core,on_disk
c
c----------------------------------------------------------
c setup to zero the following addresses for stored integrals:
c
       incorex=0
       ilab=0
       jlab=0
       klab=0
       llab=0
c
       isiz=0
       jsiz=0
       ksiz=0
       lsiz=0
c
c----------------------------------------------------------
c
      scftype='full-direct'
      if(in_core .or. on_disk) scftype='semi-direct'
c
      scftypex=scftype
      if(scftype.eq.'full-direct') return
c----------------------------------------------------------
c allocate memory for two arrays :
c (1) number of integrals in each super-block
c (2) number of quartets in each super-block
c
c Both above are the theoretical numbers (without symm. or neglect)
c They are used only for stored integrals
c
c      write(6,*)' run_type=',run_type,' scftype=',scftype
c
c     if(run_type.eq.'full-direct') then
         call getmem(nbloks,integbl)
         call getmem(nbloks,iqrtsbl)
c
         call setival('integbl',integbl)
         call setival('iqrtsbl',iqrtsbl)
c     else
c        call getival('integbl',integbl)
c        call getival('iqrtsbl',iqrtsbl)
c     endif
c----------------------------------------------------------
c find out what the price limit for stored inregrals is :
c
      call find_int_price(xntnu,intnu,xint_core,xint_disk,
     *                minpr,maxpr,iprice_core,iprice_disk)
c
c output :
c iprice_core (integ. more expensive than iprice_core - in-core)
c iprice_disk (integ. more expensive than iprice_disk - on-disk)
c----------------------------------------------------------
c decide which integrals will be stored in-core or on-disk and
c which will be recalculated
c
      if(in_core) then
c        change  prices of stored integrals to negative values :
c Total available storage in double words ON ALL nodes
         call getival('incore',incore)
         call getival('nslv', nproc)
         if(nproc.eq.0) nproc=nproc+1
c Memory needed for in-core stored intgrals is
         xmem_core=dble((incore-4)*nproc)
         call set_int_price_ftc(nbl2,bl(npard),bl(nsupb),
     *                  bl(ijbld),bl(iisd),bl(jjsd),inx,
     *                  'in_core', xmem_core, iprice_core,
     *                  sparsitx,
     *                  bl(ncost_int),bl(integbl),bl(iqrtsbl),
     *                  bl(iftc1),
     *        xmem_2_store,xint_2_store,xqrt_2_store,nblk_2_store)
c
c        allocate memory to handle stored integrals :
         call set_coredat(xmem_2_store,xint_2_store,xqrt_2_store,
     *                    nblk_2_store,
     *                    incorex,ilab,jlab,klab,llab,
     *                            isiz,jsiz,ksiz,lsiz,iqrt)
         on_disk=.false.
      endif
c
      if(on_disk) then
         call getival('indisk',indisk) ! in MB
         xmem_disk=indisk
         xmem_disk=xmem_disk*1000000/8
         call set_int_price_ftc(nbl2,bl(npard),bl(nsupb),
     *                  bl(ijbld),bl(iisd),bl(jjsd),inx,
     *                  'on_disk', xmem_disk, iprice_disk,
     *                  sparsitx,
     *                  bl(ncost_int),bl(integbl),bl(iqrtsbl),
     *                  bl(iftc1),
     *        xmem_2_store,xint_2_store,xqrt_2_store,nblk_2_store)
c
         intmem=int( 0.9d0*dble(indisk) )
         call setival('intmem',intmem)
ccccc    call set_diskdat()
      endif
c
c output from the set_price() call  : for each super-block :
c  ncost(ikbl) is now : negative for stored integrals
c                       possitive for recalculated integrlas
c                       (absolute value is the original price)
c----------------------------------------------------------
cc      write(ioutput,*)' nsym,rnsym=',nsym,rnsym,' sparsity=',sparsitx
c----------------------------------------------------------
      end
c================================================================
c================================================================
      subroutine statint_ftc(inx,iis,jjs,ijbl,npar,nbl2,
     *             ncost_int,xntnu   ,intnu   ,minpr    ,maxpr,
     *             q_size, sparsity, iftc1)
      implicit real*8 (a-h,o-z)
      common /memmax/ ispblx, maxme1,iforwhat
      common /outfile/ ioutput
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*), npar(*)
      dimension ncost_int(*)
      dimension intnu(15)
      dimension xntnu(15)
      dimension iftc1(*)
      data half /0.5d0/
c-------------------------------------------
c number of integrals more expensive than :
c    0.0%  of maximum integral-price (all integrals)
c    0.1%  of maximum integral-price
c    0.3%  of maximum integral-price
c    0.5%  of maximum integral-price
c    0.7%  of maximum integral-price
c    1.0%  of maximum integral-price
c    1.5%  of maximum integral-price
c    3.0%  of maximum integral-price
c    5.0%  of maximum integral-price
c    7.5%  of maximum integral-price
c    10.%  of maximum integral-price
c    15.%  of maximum integral-price
c    25.%  of maximum integral-price
c    50.%  of maximum integral-price
c    75.%  of maximum integral-price
c
c these prices in intnu() and number of integ. in xntnu()
c--------------------------
      do 15 ii=1,15
      xntnu(ii)=0.d0
   15 continue
c
      xprice=dble(maxpr)
      x000=0.d0
      x001=0.001d0*xprice
      x003=0.003d0*xprice
      x005=0.005d0*xprice
      x007=0.007d0*xprice
      x010=0.010d0*xprice
      x015=0.015d0*xprice
      x030=0.030d0*xprice
      x050=0.050d0*xprice
      x075=0.075d0*xprice
      x100=0.10d0*xprice
      x150=0.15d0*xprice
      x250=0.25d0*xprice
      x500=0.50d0*xprice
      x750=0.75d0*xprice
c
      intnu(1)=int(x000)
      intnu(2)=int(x001)
      intnu(3)=int(x003)
      intnu(4)=int(x005)
      intnu(5)=int(x007)
      intnu(6)=int(x010)
      intnu(7)=int(x015)
      intnu(8)=int(x030)
      intnu(9)=int(x050)
      intnu(10)=int(x075)
      intnu(11)=int(x100)
      intnu(12)=int(x150)
      intnu(13)=int(x250)
      intnu(14)=int(x500)
      intnu(15)=int(x750)
c-------------------------
      xprice=xmaxpr_bl
      x000=0.d0
      x001=0.001d0*xprice
      x003=0.003d0*xprice
      x005=0.005d0*xprice
      x007=0.007d0*xprice
      x010=0.010d0*xprice
      x015=0.015d0*xprice
      x030=0.030d0*xprice
      x050=0.050d0*xprice
      x075=0.075d0*xprice
      x100=0.10d0*xprice
      x150=0.15d0*xprice
      x250=0.25d0*xprice
      x500=0.50d0*xprice
      x750=0.75d0*xprice
c
c-------------------------
c for average size of a quartet up to (ll|ss) & (ll|ll)
c (needed to estimate number of integrals that can be hold in-core)
c
      xint_part1=0.d0
      xqrt_part1=0.d0
      xint_part2=0.d0
      xqrt_part2=0.d0
c
      xint_all=0.d0
      xqrt_all=0.d0
c-------------------------
c contracted integrals counter
c
      ikbl=0
      do 100 ibl=1,nbl2
      nparij=npar(ibl)
            ijcs1=ijbl(ibl,1)
            ics1=iis(ijcs1)
            jcs1=jjs(ijcs1)
c
            ijs1=inx(3,ics1)*inx(3,jcs1)
            ijg1=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
c
         do 100 kbl=1,ibl
         ikbl=ikbl+1
         nparkl=npar(kbl)
            klcs1=ijbl(kbl,1)
            kcs1=iis(klcs1)
            lcs1=jjs(klcs1)
c
            kls1=inx(3,kcs1)*inx(3,lcs1)
            klg1=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
c-----------------
            ijkls1=ijs1*kls1
            ijklg1=ijg1*klg1
c-----------------
            ijklq =0
            do ijpar=1,nparij
               ijcs=ijbl(ibl,ijpar)
               ics=iis(ijcs)
               jcs=jjs(ijcs)
               ijftc=iftc1(ics)+iftc1(jcs)
               if(ibl.eq.kbl) then
                  nparklx=ijpar
               else
                  nparklx=nparkl 
               endif
               do klpar=1,nparklx
                  klcs=ijbl(kbl,klpar)
                  kcs=iis(klcs)
                  lcs=jjs(klcs)
                  klftc=iftc1(kcs)+iftc1(lcs)
                  if(ijftc+klftc.gt.1) ijklq=ijklq+1
               enddo
            enddo
c-----------------
c
c  quartets & integrals in one super-block :
c
            ijklq=ijklq*ijklg1
c
            qijkl=dble(ijklq)
c
c include sparsity
c
            qijkl=qijkl*sparsity
c
            xinteg1=qijkl*dble(ijkls1)
c
c  total number of quartets & integrals:
c
            xint_all=xint_all + xinteg1
            xqrt_all=xqrt_all + qijkl
c-----------------
         nprice=ncost_int(ikbl)
c-----------------
c calculate average size of a quartet (for in-core storage)
c (average over few first blocks)
c
         if(ijkls1.le. 16) then   ! only up to (ll|ss)
           xint_part1=xint_part1 + xinteg1
           xqrt_part1=xqrt_part1 + qijkl
         endif
         if(ijkls1.le.256) then   ! only up to (ll|ll)
           xint_part2=xint_part2 + xinteg1
           xqrt_part2=xqrt_part2 + qijkl
         endif
c-----------------
c
c count number of integrals more expensive than ..
c
         do ii=1,15
          if(nprice.gt.intnu(ii)   ) xntnu(ii)=xntnu(ii)+xinteg1
         enddo
c
  100 continue
c-----------------
c average size of a quartet calculated over blocks up to (ll|ls)
c
      q_size_part1=xint_part1/xqrt_part1
      q_size_part2=xint_part2/xqrt_part2
      q_size_full=xint_all /xqrt_all
c
      q_size=(q_size_part1 + q_size_part2)*0.5d0
c-----------------
c     write(ioutput,*)' total number of quartets =',xqrt_all
c     write(ioutput,*)' total number of integrals=',xint_all
c     write(ioutput,*)' avergae size of quartet : '
c     write(ioutput,*)'  overall=',int(q_size_full),
c    *                   ' part1=',int(q_size_part1),
c    *                   ' part2=',int(q_size_part2),
c    *                   '  used=',int(q_size)
c-----------------
c
      end
c================================================
      subroutine change_price(nbloks,ncost)
      dimension ncost(*)
c
      do ikbl=1,nbloks
        iprice=ncost(ikbl)
        if(iprice.lt.0) iprice=-iprice
        ncost(ikbl)=iprice
      enddo
      end
c==============================================================
      subroutine set_int_price_ftc(nbl2,npar,nsupb,ijbl,iis,jjs,inx,
     *                     storage,xmem_store, iprice_store,
     *                     sparsity, ncost,blinteg,bliqrts,
     *                     iftc1,
     *      xmem_2_store,xint_2_store,xqrt_2_store,nblk_2_store)
      implicit real*8 (a-h,o-z)
      character*7 storage        ! in_core or on_disk
c--------------------------------------------------------------
c this routine determines which integrals will be stored
c in-core or on-disk
c
c input : storage - 'in-core' or 'on-disk'
c         xmem_store - available memory to store integrals & label's info
c         iprice_store - price above which int. will be stored
c         sparsity - predicted sparsity of integral file
c
c output: ncost(iklb) - negative price for stored integrals
c                       possitive price for recalculated
c         blinteg(ikbl)- number of integrals in a block (w/o sparsity)
c         bliqrts(ikbl)- number of quartets  in a block (w/o sparsity)
c         xmem_2_store - total memory taht WILL be used for in-core storage
c         xint_2_store - number of integrals that WILL be stored
c         xqrt_2_store - number of integrals that WILL be stored
c         nblk_2_store - number of super blocks that WILL be stored
c--------------------------------------------------------------
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /outfile/ ioutput
c
      dimension npar(nbl2),ijbl(nbl2,*)
      dimension iis(*),jjs(*),inx(12,*)
      dimension ncost(*),nsupb(*)
      dimension blinteg(*), bliqrts(*)
      dimension iftc1(*)
c---------------------------------------------------------------
cNote : sparsity as it was estimated, is used for integrals only
c for quartets (for safety) 1.5* sparsity (but not >1)
c
      q_sparsity=sparsity*1.5d0
      if(q_sparsity.gt. 1.0d0) q_sparsity=1.0d0
      if(q_sparsity.lt. 0.5d0) q_sparsity=0.5d0
c
c However, if quatets are short (e.g. sto-3g basis set) then
c sparsity for quartets should be sharper
c
      q_length=dble(ncf)/dble(ncs)
      if(q_length .le. 2.5d0 .and. ncf.le.100)  q_sparsity=1.0d0
c---------------------------------------------------------------
c     write(6,*)' storage=',storage,' memory 4 st.=', xmem_store,
c    *  '  store more exp than=', iprice_store
c---------------------------------------------------------------
c
c storage & recalculations
c
      int0=0        ! block counter of recalc blocks
      int2=0        ! block counter of stored blocks
      int22=0       ! block counter of stored blocks
c
      xnt0=0.d0     ! re-calc.integral counter
      xnt2=0.d0     ! stored integral counter
      xtot=0.d0     ! total number of integ.
c
      qtot=0.d0     ! total number of quartets
      qnt0=0.d0     ! re-calc.quartets counter
      qnt2=0.d0     ! stored  quartets counter
c
      nblk=0.d0     ! total number of small blocks
      nsm0=0.d0     ! re-calc.small blocks
      nsm2=0.d0     ! stored  small blocks
c
      xmem_used=0.d0
c
      ikbl=0
      ikbl1=0
      do ibl=1,nbl2
         ijpar=npar(ibl)
         ijcs1=ijbl(ibl,1)
         ics1=iis(ijcs1)
         jcs1=jjs(ijcs1)
         ijs1=inx(3,ics1)*inx(3,jcs1)
         ijg1=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
         do kbl=1,ibl
            klpar=npar(kbl)
            klcs1=ijbl(kbl,1)
            kcs1=iis(klcs1)
            lcs1=jjs(klcs1)
            kls1=inx(3,kcs1)*inx(3,lcs1)
            klg1=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
c
            ikbl=ikbl+1
            iprice=ncost(ikbl)
            ismall=nsupb(ikbl)
c
            ijklq=0
            do ijp=1,ijpar
               ijcs=ijbl(ibl,ijp)
               ics=iis(ijcs)
               jcs=jjs(ijcs)
               ijftc=iftc1(ics)+iftc1(jcs)
               klparx=klpar
               if(ibl.eq.kbl) klparx=ijp
               do klp=1,klparx
                  klcs=ijbl(kbl,klp)
                  kcs=iis(klcs)
                  lcs=jjs(klcs)
                  klftc=iftc1(kcs)+iftc1(lcs)
                  if(ijftc+klftc.gt.1) ijklq=ijklq+1
               enddo
            enddo
            if(ijklq.gt.0) then
               ikbl1=ikbl1+1
               nblk=nblk+ismall
            endif
c
c
            ijklq=ijklq*ijg1*klg1
            q=dble(ijklq)
c
            xinteg=q*dble(ijs1*kls1)
c
c setup the blinteg(ikbl) and bliqrts(ikbl) arrays :
c
            blinteg(ikbl)=xinteg
            bliqrts(ikbl)=q
c
c include sparsity for integrals:
c
c2005       x=xinteg*sparsity
            x=xinteg
c
c include q_sparsity for quartets :
c
c2005       q=q*q_sparsity
c
            xtot=xtot+x
            qtot=qtot+q
c
c memory
            xmem = x + q
c
            if(iprice.gt.iprice_store) then
               if(xmem_used+xmem.LT.xmem_store) then
                  xmem_used=xmem_used+xmem
                  xnt2=xnt2+x
                  if(ijklq.gt.0) int2=int2+1
                  int22=int22+1
                  ncost(ikbl)=-iprice
                  qnt2=qnt2+q
                  if(ijklq.gt.0) nsm2=nsm2+ismall
               else
                  xnt0=xnt0+x
                  if(ijklq.gt.0) int0=int0+1
cccccc            ncost(ikbl)= no change
                  qnt0=qnt0+q
                  if(ijklq.gt.0) nsm0=nsm0+ismall
               endif
            else
                  xnt0=xnt0+x
                  int0=int0+1
cccccc            ncost(ikbl)= no change
                  qnt0=qnt0+q
                  if(ijklq.gt.0) nsm0=nsm0+ismall
            endif
c
c..         write(ioutput,*)
c..  *     ' block=',ikbl,' prices=',iprice,ncost(ikbl)
c
         enddo
      enddo
c-------------------------------------------------
c-------------------------------------------------
ctest only
c
c     ilt0=0
c     igt0=0
c     do isbl=1,nbl2*(nbl2+1)/2
c         iprice=ncost(ikbl)
c         if(iprice.lt.0) then
c            ilt0=ilt0+1
c         else
c            igt0=igt0+1
c         endif
c     enddo
c
c     write(6,*)'1 blocks with negative price=',ilt0
c     write(6,*)'1 blocks with positive price=',igt0
c-------------------------------------------------
      if(int22.eq.nbl2*(nbl2+1)/2) go to 9876
c-------------------------------------------------
c try to store more integrals from remaining blocks
c
      ikbl=0
      do ibl=1,nbl2
         ijpar=npar(ibl)
         ijcs1=ijbl(ibl,1)
         ics1=iis(ijcs1)
         jcs1=jjs(ijcs1)
         ijs1=inx(3,ics1)*inx(3,jcs1)
         ijg1=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
         do kbl=1,ibl
            klpar=npar(kbl)
            klcs1=ijbl(kbl,1)
            kcs1=iis(klcs1)
            lcs1=jjs(klcs1)
            kls1=inx(3,kcs1)*inx(3,lcs1)
            klg1=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
            ikbl=ikbl+1
            iprice=ncost(ikbl)
            if(iprice.gt.0) then
               ismall=nsupb(ikbl)
c
               ijklq=0
               do ijp=1,ijpar
                  ijcs=ijbl(ibl,ijp)
                  ics=iis(ijcs)
                  jcs=jjs(ijcs)
                  ijftc=iftc1(ics)+iftc1(jcs)
                  klparx=klpar
                  if(ibl.eq.kbl) klparx=ijp
                  do klp=1,klparx
                     klcs=ijbl(kbl,klp)
                     kcs=iis(klcs)
                     lcs=jjs(klcs)
                     klftc=iftc1(kcs)+iftc1(lcs)
                     if(ijftc+klftc.gt.1) ijklq=ijklq+1
                  enddo
               enddo
c
               ijklq=ijklq*ijg1*klg1
c
               q=dble(ijklq)
c
               xinteg=q*dble(ijs1*kls1)
c
c include sparsity for integrals :
c
c2005          x=xinteg*sparsity
               x=xinteg
c
c include q_sparsity for quartets :
c
c2005          q=q*q_sparsity
c
c memory
               xmem = x + q
c
               if(xmem_used+xmem.LT.xmem_store) then
                  xmem_used=xmem_used+xmem
                  xnt2=xnt2+x
                  if(ijklq.gt.0) int2=int2+1
                  int22=int22+1
                  ncost(ikbl)=-iprice
                  qnt2=qnt2+q
                  if(ijklq.gt.0) nsm2=nsm2+ismall
c
                  xnt0=xnt0-x
                  if(ijklq.gt.0) int0=int0-1
                  qnt0=qnt0-q
                  if(ijklq.gt.0) nsm0=nsm0-ismall
               endif
            endif
         enddo
      enddo
c-------------------------------------------------
 9876 continue
c-------------------------------------------------
c2005 nbloks=nbl2*(nbl2+1)/2
      nbloks=ikbl1   ! non-FTC blocks
c-------------------------------------------------
      xmem_2_store=xmem_used
      xint_2_store=xnt2
      xqrt_2_store=qnt2
      int2=(int2*15)/10 ! increase by 1.5
      nblk_2_store=int2
      if(nblk_2_store.gt.nbloks) nblk_2_store=nbloks
c-------------------------------------------------
ctest only
c
c     ilt0=0
c     igt0=0
c     do isbl=1,nbl2*(nbl2+1)/2
c         iprice=ncost(ikbl)
c         if(iprice.lt.0) then
c            ilt0=ilt0+1
c         else
c            igt0=igt0+1
c         endif
c     enddo
c
c     write(6,*)'2 blocks with negative price=',ilt0
c     write(6,*)'2 blocks with positive price=',igt0
c     call f_lush(6)
c-------------------------------------------------
      write(ioutput,199)
  199 format(/'Prediction for non-FTC integral storage '/)
      write(ioutput,200) nbloks,storage,  int2,int0
  200 format(
     *' Total number of big blocks =',i15/
     *'  ',a7, ' stored big blocks =',i15/
     *'    recalculated big blocks =',i15/)
c
      write(ioutput,210) nblk ,nsm2, nsm0
  210 format(
     *' Total number of small blks =',i15/
     *'          stored small blks =',i15/
     *'    recalculated small blks =',i15/)
c
      xtot1=1.d0/xtot
      x2p=xnt2*xtot1*100.d0
      x0p=xnt0*xtot1*100.d0
c
      write(ioutput,220) xtot,xnt2,xnt0
  220 format(
     *' Total number of integrals  =',f15.0/
     *'          stored integrals  =',f15.0/
     *'    recalculated integrals  =',f15.0/)
c
      write(ioutput,230) qtot,qnt2,qnt0
  230 format(
     *' Total number of quartests  =',f15.0/
     *'          stored quartests  =',f15.0/
     *'    recalculated quartests  =',f15.0/)
c
ctest
c       do ii=1,nbloks
c          xx=blinteg(ii)
c          qq=bliqrts(ii)
c          write(6,*)'block=',ii,' no int=',xx,' no qrts=',qq
c       enddo
ctest
c
      end
