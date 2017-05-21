c=====================================================================
      subroutine int_cphf(nblocks,bl,inx,ntri,thres1,
     *                    dscreen,dens,fock,labels,mywork,
     *                    igran)
c---------------------------------------------------------------------
c This routine is called from coupled-perturbed HF for all three modes :
c           non -, semi-, and full-direct.
c Three Fock matrices are constructed F(DX,g0),F(DY,g0) and F(DZ,g0)
c using stored integrals first (if any) and then re-calculated integrals
c
c For non - and semi-direct mode the stored integrals are used first
c and corresponding part of Fock matrix is built .Then, for semi-direct
c mode remaining integrals are calculated (with where='fock') and
c building of the Fock matrix is finished.
c
c For recalculated integrals value of where='fock' is set up HERE
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      logical rescale
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NOT EXISTING super block
c---------------------------------------------------------------------
c known from int_store   ;
      common /datstore/ thres_stored,isto,isecond,ito
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
c
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),dscreen(*),dens(ntri,3),fock(ntri,3)
      dimension labels(*)
      call secund(time0)
c----------------------------------------------------------------
c
c  if it is a full-direct run, reset common datstore
c  and skip storage retrieval
c
      if(scftype.eq.'full-direct')then
        isto=0
        ito=0
        isecond=0
        goto 9999
      endif
c
c  use stored integrals only once per iteration
c
      if(isto.eq.0) go to 9999
c----------------------------------------------------------------
c First construct a part of fock matrix with stored integrals:
c (those integrals were calculated by calling int_store routine)
c Do it only once each iteration (isto=1 means Use stored!)
c
c       return threshold used for stored integrals:
c
        thres0=thres_stored
c make sure stored integrals are only used once each cycle
        isto=0
c
c get number of stored big-blocks and address of stored quartets array
c
c       first check where integrals are stored (in-core or on-disk):
c
        if(incorex.ne.0) then
c
c          in-core storage :
c
           call fock_core_cphf(
     $        bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $        bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $        bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     bl(lind),
     $        dens,        fock,       ntri)
        else
c
c          on-disk storage :
c
           rewind(160)
           rewind(161)
           rewind(162)
           rewind(163)
           rewind(164)
           rewind(165)
           rewind(166)
           rewind(167)
           nblcount=0
           rescale=.false.
 9876      continue
           call read_33(160,iblstore,iqstore,nsplit) 
           if(iblstore.eq.0) go to 9999
           call getint(nsplit*2,ifromto)
           call getint(nsplit*5,nofinteg)
           call read_44(165,bl(ifromto),bl(nofinteg),nsplit)
           nblcount=nblcount+iblstore
           if(nblcount.EQ.nblstore) rescale=.true.
           call getint_2(4*iblstore , ijklsiz) ! I*2
           call getint_2(  iblstore , nquarts) ! I*2
           call getint_2(4*iqstore  , ijkllab) ! I*2
           call getival('intbu72',integbu72)
           call getival('intbu62',integbu62)
           call getival('intbuf2',integbuf2)
           call getival('intbuf4',integbuf4)
           call getival('intbuf8',integbuf8)
           call getmem(integbuf8,ixint8)
           call getint_4(integbuf4,ixint4)
           call getint_2(integbuf2,ixint2)
           call getint_2(integbu62,ixin62)
           call getint_2(integbu72,ixin72)
           call fock_disk_cphf(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   dens,       fock,
     $        ncs,         ntri)
           call retmem(10)
           if(.not.rescale) go to 9876
        endif
ccccc   call dens_fock_prt(dens,fock,'part_stored')
c
        if(scftype.eq.'non -direct') then
           call secund(time1)
c          call term_info(thres1, 0.d0,time1-time0,where)
C???       return
           go to 8888
        endif
c
 9999 continue
c
c  check if we actually have to compute any integrals,
c  because it might happen that a slave executes this routine
c  even with mywork=-1, in order to retrieve the stored integrals
c  (this is to avoid a race condition MM 08/29/2006 )
c
      if(mywork.lt.0)return
c
c----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c re-calculate integrals
c
      where='fock'
ccc   screen='cphf'
      screen='cph2'
c
c set inetgral timings to zero as well as neglect stuff
c
      call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c total number of super-blocks  nsupblks=nbl2*(nbl2+1)/2
c
      nsupblks=nblocks
      if(mywork.eq.0) then
         istart=1
         istop=nsupblks
      else
         istart=mywork
         istop=mywork+igran-1
      end if
c----------------------------------------------------------------
      do isupbl=istart,istop
c
c        get integral's price
c
         call get_price(bl(ncost),isupbl,iprice)
c blocks were stored:
c  - if they have negative price and they are before the last stored
c  - if they have negative price and all negative blocks were stored
c  - if all negative blocks were stored and our block is before the last
         if((iprice.LE.0.and.isupbl.le.ito).or.(iprice.le.0.and.
     *   isecond.eq.1).or.(isecond.eq.1.and.isupbl.le.ito)) go to 1111
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,dscreen,where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
c
c check if requested super-block was in a list :
c
        if(stopnow) go to 1111
c
c integrals arrived in bl(ibuffz); check if they are :
c
c  nintez is the size of a given quartet.It is set up
c  to zero if there is no integrals
c
        if(nintez.gt.0) then
           call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
           call fock_cphf(bl(ibuffz),dens,fock,bl(lind),ntri,
     *                    nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
        endif
c
c----------------------------------------------------------------
c
        if(moreint) go to 11
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
 8888 continue
c
c re-scale fock matrices :
c
c     call rescale_fock(fock,ntri,bl(lind),thres1)
c
c   moved to the para_cphf routine .
c
      call zerout_fock_diag(fock,ntri,bl(lind))
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
      call secund(time2)
      call term_info(thres1,time2-time1,time1-time0,where)
c----------------------------------------------------------------
c
      end
c=====================================================================
c Three fock-builder routines follows :
c
c   1. using recalculated integrals from buffer : fock_cphf
c   2. using stored integrals (stored in-core   : fock_core_cphf
c   3. using stored integrals (stored on-disk   : fock_disk_cphf
c
c=====================================================================
      subroutine fock_cphf(buf,dens,fock,lind,ntri,
     *                     nbls,ngcd,lnijkl,labels,length,lgenct )
c----------------------------------------------------------------
c Fock builder (using integrals in buffer)
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension dens(ntri,3),fock(ntri,3)
      dimension lind(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
      do 100 ijklp=1,nbls
      ngcq=lgenct(ijklp)
        do 150 iqu=1,ngcq
        icff=labels(1,iqu,ijklp)
        jcff=labels(2,iqu,ijklp)
        kcff=labels(3,iqu,ijklp)
        lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
          integ=0
          do 200 iii=1,ilen
          icf=icff+iii
          ii=lind(icf)
          do 200 jjj=1,jlen
          jcf=jcff+jjj
          jj=lind(jcf)
          do 200 kkk=1,klen
          kcf=kcff+kkk
          kk=lind(kcf)
c ik
             if(icf.ge.kcf) then
                ikf=ii+kcf
                fik=1.d0
             else
                ikf=kk+icf
                fik=-1.d0
             endif
                xik=dens(ikf,1)
                yik=dens(ikf,2)
                zik=dens(ikf,3)
c jk
             if(jcf.ge.kcf) then
                jkf=jj+kcf
                fjk=1.d0
             else
                jkf=kk+jcf
                fjk=-1.d0
             endif
                xjk=dens(jkf,1)
                yjk=dens(jkf,2)
                zjk=dens(jkf,3)
c
          do 200 lll=1,llen
          lcf=lcff+lll
          ll=lind(lcf)
c il
             if(icf.ge.lcf) then
                ilf=ii+lcf
                fil=1.d0
             else
                ilf=ll+icf
                fil=-1.d0
             endif
                xil=dens(ilf,1)
                yil=dens(ilf,2)
                zil=dens(ilf,3)
c jl
             if(jcf.ge.lcf) then
                jlf=jj+lcf
                fjl=1.d0
             else
                jlf=ll+jcf
                fjl=-1.d0
             endif
                xjl=dens(jlf,1)
                yjl=dens(jlf,2)
                zjl=dens(jlf,3)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
ctest
c             write(6,1234) icf,jcf,kcf,lcf,xint0*1.d-10
c1234 format(4i3,1x,2(f12.8,1x))
c-------------------------------------------------------
c Coulomb part is zero for anti-symmetric perturbation
c Only exchange part is present
c-------------------------------------------------------
c         magnitude checking again
          if(abs(xint0).gt.0.5d0) then
ctest
c             write(6,1234) icf,jcf,kcf,lcf,xint0*1.d-10
c1234 format(4i3,1x,2(f12.8,1x))
c
              xiljk=xint0*fil*fjk
              xjlik=xint0*fjl*fik
c
c il,jl, ik,jk
c
                fock(ilf,1)=fock(ilf,1) + xjk*xiljk
                fock(jlf,1)=fock(jlf,1) + xik*xjlik
                fock(ikf,1)=fock(ikf,1) + xjl*xjlik
                fock(jkf,1)=fock(jkf,1) + xil*xiljk
c
                fock(ilf,2)=fock(ilf,2) + yjk*xiljk
                fock(jlf,2)=fock(jlf,2) + yik*xjlik
                fock(ikf,2)=fock(ikf,2) + yjl*xjlik
                fock(jkf,2)=fock(jkf,2) + yil*xiljk
c
                fock(ilf,3)=fock(ilf,3) + zjk*xiljk
                fock(jlf,3)=fock(jlf,3) + zik*xjlik
                fock(ikf,3)=fock(ikf,3) + zjl*xjlik
                fock(jkf,3)=fock(jkf,3) + zil*xiljk
c
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c==============================================================
      subroutine fock_core_cphf(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        dens,        fock,       ntri)
c----------------------------------------------------------------
c input :
c core(*) - array with stored integrals
c i-llab()- starting indecies for integrals
c i-lsiz()- shell's sizes for each super-block
c iqrt()  - number of stored quartets for each super-block
c integrals - number of stored integrals
c nqstore   - number of stored quartets
c nblstore  - number of stored big-blocks
c----------------------------------------------------------------
c Fock builder (using in-core integrals )
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dens(ntri,3),fock(ntri,3)
      dimension lind(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      intx=0
      ijklp=0
      do isbl=1,nblstore
         ilen=isiz(isbl)
         jlen=jsiz(isbl)
         klen=ksiz(isbl)
         llen=lsiz(isbl)
c
         nqrt1=iqrt(isbl)
c
         do iqrt1=1,nqrt1
            ijklp=ijklp+1
            icff=ilab(ijklp)
            jcff=jlab(ijklp)
            kcff=klab(ijklp)
            lcff=llab(ijklp)
c
         do iii=1,ilen
            icf=icff+iii
            do jjj=1,jlen
               jcf=jcff+jjj
               do kkk=1,klen
                  kcf=kcff+kkk
                  do lll=1,llen
                     intx=intx+1
                     xint0=core(intx)
                     lcf=lcff+lll
c-----------------------------------------------------------
                     ii=lind(icf)
                     jj=lind(jcf)
                     kk=lind(kcf)
                     ll=lind(lcf)
c
      call get_dxyz(icf,kcf,ii,kk,dens,ntri, ikf,xik,yik,zik)
      call get_dxyz(jcf,kcf,jj,kk,dens,ntri, jkf,xjk,yjk,zjk)
      call get_dxyz(icf,lcf,ii,ll,dens,ntri, ilf,xil,yil,zil)
      call get_dxyz(jcf,lcf,jj,ll,dens,ntri, jlf,xjl,yjl,zjl)
c
c-----------------------------------------------------------
c  ***  il  ***
c
              if( icf.gt.lcf ) then
                fock(ilf,1)=fock(ilf,1) + xjk*xint0
                fock(ilf,2)=fock(ilf,2) + yjk*xint0
                fock(ilf,3)=fock(ilf,3) + zjk*xint0
              endif
              if( lcf.gt.icf ) then
                fock(ilf,1)=fock(ilf,1) - xjk*xint0
                fock(ilf,2)=fock(ilf,2) - yjk*xint0
                fock(ilf,3)=fock(ilf,3) - zjk*xint0
              endif
c
c  ***  jl  ***
c
              if( jcf.gt.lcf ) then
                fock(jlf,1)=fock(jlf,1) + xik*xint0
                fock(jlf,2)=fock(jlf,2) + yik*xint0
                fock(jlf,3)=fock(jlf,3) + zik*xint0
              endif
              if( lcf.gt.jcf ) then
                fock(jlf,1)=fock(jlf,1) - xik*xint0
                fock(jlf,2)=fock(jlf,2) - yik*xint0
                fock(jlf,3)=fock(jlf,3) - zik*xint0
              endif
c
c  ***  ik  ***
c
              if( icf.gt.kcf ) then
                fock(ikf,1)=fock(ikf,1) + xjl*xint0
                fock(ikf,2)=fock(ikf,2) + yjl*xint0
                fock(ikf,3)=fock(ikf,3) + zjl*xint0
              endif
              if( kcf.gt.icf ) then
                fock(ikf,1)=fock(ikf,1) - xjl*xint0
                fock(ikf,2)=fock(ikf,2) - yjl*xint0
                fock(ikf,3)=fock(ikf,3) - zjl*xint0
              endif
c
c  ***  jk  ***
c
              if( jcf.gt.kcf ) then
                fock(jkf,1)=fock(jkf,1) + xil*xint0
                fock(jkf,2)=fock(jkf,2) + yil*xint0
                fock(jkf,3)=fock(jkf,3) + zil*xint0
              endif
              if( kcf.gt.jcf ) then
                fock(jkf,1)=fock(jkf,1) - xil*xint0
                fock(jkf,2)=fock(jkf,2) - yil*xint0
                fock(jkf,3)=fock(jkf,3) - zil*xint0
              endif
c-----------------------------------------------------------
                  enddo         !   do lll=1,llen
               enddo            !   do kkk=1,klen
            enddo               !   do jjj=1,jlen
         enddo                  !   do iii=1,ilen
         enddo                  !   do iqrt1=1,nqrt1
      enddo                     !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
         resc=thres0/thres1
         call dscal(ntri*3,resc,fock,1)
c        do 200 ij=1,ntri
c           fock(ij,1)=fock(ij,1)*resc
c           fock(ij,2)=fock(ij,2)*resc
c           fock(ij,3)=fock(ij,3)*resc
c 200    continue
      endif
c
      end
c=====================================================================
      subroutine fock_disk_cphf(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dens,       fock,
     $        ncs,         ntri)
c----------------------------------------------------------------
c input :
c core(*) - array with stored integrals
c i-llab()- starting indecies for integrals
c i-lsiz()- shell's sizes for each super-block
c nquarts()  - number of stored quartets for each super-block
c integrals - number of stored integrals
c nqstore   - number of stored quartets
c nblstore  - number of stored big-blocks
c----------------------------------------------------------------
c Fock builder (using on-disk integrals )
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical rescale
      logical doint72,doint62
c----------------------------------------------------------------
      common /screening/ dens_max_el 
c----------------------------------------------------------------
      integer*2 ijkllab(4,iqstore) !   ilab,jlab,klab,llab
      integer*2 ijklsiz(4,iblstore)!   isiz(isbl),jsiz,ksiz,lsiz  
      integer*2 nquarts(iblstore)     !   no of quartets in each block
c----------------------------------------------------------------
      integer*2 xinte72(*)
      integer*2 xinte62(*)
      integer*2 xinteg2(*)
      integer*4 xinteg4(*)
      dimension xinteg8(*)
      dimension dens(ntri,3),fock(ntri,3)
c
      dimension lind(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      data half /0.5d0/
c----------------------------------------------------------------
c integrals on-disk have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
c units for stored stuff :
c
      iunit0=160 ! for no of stored bl & qr
      iunit1=161 ! sizes
      iunit2=162 ! quarts
      iunit3=163 ! labels
      iunit4=164 ! integrals in  real*8 & i*4 & i*2
      iunit6=166 ! integrals in  i*2
      iunit7=167 ! integrals in  i*2
c----------------------------------------------------------------
c maximum real value of integrals stored in xinteg2() 
c
      xmax7=0.5d-7
      xmax6=0.5d-6
c 
c----------------------------------------------------------------
c kept as xmax7/thres0 ! e.g. 0.5*10**-7/10**-10= 500
c
c if actual threshold thres1 is bigger then or equal to 10**-7
c then we do NOT need these small integrals. They are needed 
c only if thres1 < 10**-7. 
c Integrlas have to be 
c
c         x > 0.5*thres1
c Including maximum density d gives 
c      
c      d* x > 0.5*thres1
c
c Thus we have :
c
c   d* 0.5 *10**-7 > d*x > 0.5*thres1
c
c which finally yields
c
c   d* 10**-7 > thres1
c
c----------------------------------------------------------------
c Thresholds ratio is
c
      t01=thres0/thres1
c----------------------------------------------------------------
      dens_max=4.0d0*sqrt(dens_max_el)
      dens_max=min(dens_max,1.0d0)
c
      doint72=.false.
      if(dens_max*1.0d-7.GT.thres1) doint72=.true.
      doint62=.false.
      if(dens_max*1.0d-6.GT.thres1) doint62=.true.
c----------------------------------------------------------------
c     write(6,666) dens_max,thres1,doint72,doint62
c666  format('   dens_max=',f10.6,' thers1=',1pe8.1' doint=',2(L2,1x))
c----------------------------------------------------------------
           dmax1=1.0d0
           dmax=dmax1*t01
c----------------------------------------------------------------
c
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
      ijklp=0
      do isplit=1,nsplit
      isbl_beg=ifromto(1,isplit)
      isbl_end=ifromto(2,isplit)
      integi72=nofinteg(1,isplit)
      integi62=nofinteg(2,isplit)
      integin2=nofinteg(3,isplit)
      integin4=nofinteg(4,isplit)
      integin8=nofinteg(5,isplit)
c
      if(doint72) then
      if(integi72.gt.0) then
         call read_int2(iunit6,xinte72,integi72)
      endif
c
      endif
      if(doint62) then
      if(integi62.gt.0) then
         call read_int2(iunit7,xinte62,integi62)
      endif
      endif
c
      if(integin2.gt.0) then
         call read_int2(iunit4,xinteg2,integin2)
      endif
c
      if(integin4.gt.0) then
         call read_int4(iunit4,xinteg4,integin4)
      endif
c
      if(integin8.gt.0) then
         call read_int8(iunit4,xinteg8,integin8)
      endif
c
      intx8=0
      intx4=0
      intx2=0
      int62=0
      int72=0
      do isbl=isbl_beg,isbl_end
         ilen=ijklsiz(1,isbl)
         jlen=ijklsiz(2,isbl)
         klen=ijklsiz(3,isbl)
         llen=ijklsiz(4,isbl)
         nqrt=nquarts(isbl)
         isize=ilen*jlen*klen*llen
         do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
c
            icase=5
            if(integin8.gt.0) then
               if(icff.lt.0) then
                  icase=1
                  icff=-icff
               endif
            endif
c
            if(integin4.gt.0.and.icase.eq.5) then
               if(jcff.lt.0) then
                  icase=2
                  jcff=-jcff
               endif
            endif
c
            if(integin2.gt.0.and.icase.eq.5) then
               if(kcff.lt.0) then
                  icase=3
                  kcff=-kcff
               endif
            endif
c
            if(integi62.gt.0.and.icase.eq.5) then
               if(lcff.lt.0) then
                  icase=4
                  lcff=-lcff
               endif
            endif
c--------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
c
              if(integin8.gt.intx8 .and. icase.eq.1) then
                 call make_fock8_cphf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dens,        fock,       lind,       intx8,
     $             xinteg8,     ntri)
              endif
c
              if(integin4.gt.intx4 .and. icase.eq.2) then
                 call make_fock4_cphf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dens,        fock,       lind,       intx4,
     $             xinteg4,     ntri)
              endif
c
              if(integin2.gt.intx2 .and. icase.eq.3) then
                 call make_fock2_cphf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dens,        fock,       lind,       intx2,
     $             xinteg2,     ntri)
              endif
c
              if(doint62) then
               if(integi62.gt.int62 .and. icase.eq.4) then
                 if(dmax1*1.0d-6.gt.thres1) then
                   call make_fock2_cphf(
     $               ilen,        jlen,       klen,       llen,
     $               icff,        jcff,       kcff,       lcff,
     $               dens,        fock,       lind,       int62,
     $               xinte62,     ntri)
                 endif
               endif
              endif
c
              if(doint72) then
               if(integi72.gt.int72 .and. icase.eq.5) then
                 if(dmax1*1.0d-7.gt.thres1) then
                   call make_fock2_cphf(
     $               ilen,        jlen,       klen,       llen,
     $               icff,        jcff,       kcff,       lcff,
     $               dens,        fock,       lind,       int72,
     $               xinte72,     ntri)
                 endif
               endif
              endif
c
         enddo                  !   do iqrt=1,nqrt
      enddo                     !   do isbl=1,nblstore
      enddo                     !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
         resc=thres0/thres1
         call dscal(ntri*3,resc,fock,1)
      endif
c
      end
c=====================================================================
      subroutine get_dxyz(icf,kcf,ii,kk,dens,ntri, ikf,xik,yik,zik)
      implicit real*8 (a-h,o-z)
      dimension dens(ntri,3)
c-- output : ikf,xik,yik,zik
c
          if(icf.ge.kcf) then
            ikf=ii+kcf
            xik=dens(ikf,1)
            yik=dens(ikf,2)
            zik=dens(ikf,3)
          else
            ikf=kk+icf
            xik=-dens(ikf,1)
            yik=-dens(ikf,2)
            zik=-dens(ikf,3)
          endif
c
      end
c=====================================================================
      subroutine select_dens(dens,ntri, select)
      implicit real*8 (a-h,o-z)
      dimension dens(ntri,3),select(ntri)
c
      do 100 i=1,ntri
      x=abs( dens(i,1) )
      y=abs( dens(i,2) )
      z=abs( dens(i,3) )
c
      select(i)=max(x,y,z)
c
  100 continue
c
      end
c=====================================================================
      subroutine rescale_fock(fock,ntri,lind,thres1)
      implicit real*8 (a-h,o-z)
      dimension fock(ntri,3)
      dimension lind(*)
      data half,two /0.5d0,2.d0/
c
      call getival('ncf',ncf)
c
      factor=thres1*half
c
      do i=1,ncf
         ii=lind(i)+i
         fock(ii,1)=fock(ii,1)*two
         fock(ii,2)=fock(ii,2)*two
         fock(ii,3)=fock(ii,3)*two
      enddo
c
      do ij=1,ntri
         fock(ij,1)=fock(ij,1)*factor
         fock(ij,2)=fock(ij,2)*factor
         fock(ij,3)=fock(ij,3)*factor
      enddo
c
      end
c=====================================================================
      subroutine zerout_fock_diag(fock,ntri,lind)
      implicit real*8 (a-h,o-z)
      dimension fock(ntri,3)
      dimension lind(*)
      data zero /0.0d0/
c
      call getival('ncf',ncf)
c
      do i=1,ncf
         ii=lind(i)+i
         fock(ii,1)=zero
         fock(ii,2)=zero
         fock(ii,3)=zero
      enddo
c
      end
c=====================================================================
      subroutine make_fock8_cphf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dens,        fock,       lind,       intx8,
     $        xinteg8,    ntri)
      implicit real*8 (a-h,o-z)
      dimension dens(ntri,3),fock(ntri,3)
      dimension xinteg8(*)
      dimension lind(*)
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
             do lcf=lcff+1,lcff+llen
               intx8=intx8+1
               xint0=xinteg8(intx8)
               ii=lind(icf)
               jj=lind(jcf)
               kk=lind(kcf)
               ll=lind(lcf)
c
               call get_dxyz(icf,kcf,ii,kk,dens,ntri, ikf,xik,yik,zik)
               call get_dxyz(jcf,kcf,jj,kk,dens,ntri, jkf,xjk,yjk,zjk)
               call get_dxyz(icf,lcf,ii,ll,dens,ntri, ilf,xil,yil,zil)
               call get_dxyz(jcf,lcf,jj,ll,dens,ntri, jlf,xjl,yjl,zjl)
c
c             ***  il  ***
c
              if( icf.gt.lcf ) then
                fock(ilf,1)=fock(ilf,1) + xjk*xint0
                fock(ilf,2)=fock(ilf,2) + yjk*xint0
                fock(ilf,3)=fock(ilf,3) + zjk*xint0
              endif
              if( lcf.gt.icf ) then
                fock(ilf,1)=fock(ilf,1) - xjk*xint0
                fock(ilf,2)=fock(ilf,2) - yjk*xint0
                fock(ilf,3)=fock(ilf,3) - zjk*xint0
              endif
c
c             ***  jl  ***
c
              if( jcf.gt.lcf ) then
                fock(jlf,1)=fock(jlf,1) + xik*xint0
                fock(jlf,2)=fock(jlf,2) + yik*xint0
                fock(jlf,3)=fock(jlf,3) + zik*xint0
              endif
              if( lcf.gt.jcf ) then
                fock(jlf,1)=fock(jlf,1) - xik*xint0
                fock(jlf,2)=fock(jlf,2) - yik*xint0
                fock(jlf,3)=fock(jlf,3) - zik*xint0
              endif
c
c             ***  ik  ***
c
              if( icf.gt.kcf ) then
                fock(ikf,1)=fock(ikf,1) + xjl*xint0
                fock(ikf,2)=fock(ikf,2) + yjl*xint0
                fock(ikf,3)=fock(ikf,3) + zjl*xint0
              endif
              if( kcf.gt.icf ) then
                fock(ikf,1)=fock(ikf,1) - xjl*xint0
                fock(ikf,2)=fock(ikf,2) - yjl*xint0
                fock(ikf,3)=fock(ikf,3) - zjl*xint0
              endif
c
c             ***  jk  ***
c
              if( jcf.gt.kcf ) then
                fock(jkf,1)=fock(jkf,1) + xil*xint0
                fock(jkf,2)=fock(jkf,2) + yil*xint0
                fock(jkf,3)=fock(jkf,3) + zil*xint0
              endif
              if( kcf.gt.jcf ) then
                fock(jkf,1)=fock(jkf,1) - xil*xint0
                fock(jkf,2)=fock(jkf,2) - yil*xint0
                fock(jkf,3)=fock(jkf,3) - zil*xint0
              endif

c
            enddo         !   over lll=1,llen
          enddo            !   over kkk=1,klen
        enddo               !   over jjj=1,jlen
      enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_fock4_cphf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dens,        fock,       lind,       intx4,
     $        xinteg4,    ntri)
      implicit real*8 (a-h,o-z)
      dimension dens(ntri,3),fock(ntri,3)
      integer*4 xinteg4(*)
      dimension lind(*)
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
             do lcf=lcff+1,lcff+llen
               intx4=intx4+1
               xint0=xinteg4(intx4)
               ii=lind(icf)
               jj=lind(jcf)
               kk=lind(kcf)
               ll=lind(lcf)
c
               call get_dxyz(icf,kcf,ii,kk,dens,ntri, ikf,xik,yik,zik)
               call get_dxyz(jcf,kcf,jj,kk,dens,ntri, jkf,xjk,yjk,zjk)
               call get_dxyz(icf,lcf,ii,ll,dens,ntri, ilf,xil,yil,zil)
               call get_dxyz(jcf,lcf,jj,ll,dens,ntri, jlf,xjl,yjl,zjl)
c
c             ***  il  ***
c
              if( icf.gt.lcf ) then
                fock(ilf,1)=fock(ilf,1) + xjk*xint0
                fock(ilf,2)=fock(ilf,2) + yjk*xint0
                fock(ilf,3)=fock(ilf,3) + zjk*xint0
              endif
              if( lcf.gt.icf ) then
                fock(ilf,1)=fock(ilf,1) - xjk*xint0
                fock(ilf,2)=fock(ilf,2) - yjk*xint0
                fock(ilf,3)=fock(ilf,3) - zjk*xint0
              endif
c
c             ***  jl  ***
c
              if( jcf.gt.lcf ) then
                fock(jlf,1)=fock(jlf,1) + xik*xint0
                fock(jlf,2)=fock(jlf,2) + yik*xint0
                fock(jlf,3)=fock(jlf,3) + zik*xint0
              endif
              if( lcf.gt.jcf ) then
                fock(jlf,1)=fock(jlf,1) - xik*xint0
                fock(jlf,2)=fock(jlf,2) - yik*xint0
                fock(jlf,3)=fock(jlf,3) - zik*xint0
              endif
c
c             ***  ik  ***
c
              if( icf.gt.kcf ) then
                fock(ikf,1)=fock(ikf,1) + xjl*xint0
                fock(ikf,2)=fock(ikf,2) + yjl*xint0
                fock(ikf,3)=fock(ikf,3) + zjl*xint0
              endif
              if( kcf.gt.icf ) then
                fock(ikf,1)=fock(ikf,1) - xjl*xint0
                fock(ikf,2)=fock(ikf,2) - yjl*xint0
                fock(ikf,3)=fock(ikf,3) - zjl*xint0
              endif
c
c             ***  jk  ***
c
              if( jcf.gt.kcf ) then
                fock(jkf,1)=fock(jkf,1) + xil*xint0
                fock(jkf,2)=fock(jkf,2) + yil*xint0
                fock(jkf,3)=fock(jkf,3) + zil*xint0
              endif
              if( kcf.gt.jcf ) then
                fock(jkf,1)=fock(jkf,1) - xil*xint0
                fock(jkf,2)=fock(jkf,2) - yil*xint0
                fock(jkf,3)=fock(jkf,3) - zil*xint0
              endif

c
            enddo         !   over lll=1,llen
          enddo            !   over kkk=1,klen
        enddo               !   over jjj=1,jlen
      enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_fock2_cphf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dens,        fock,       lind,       intx2,
     $        xinteg2,    ntri)
      implicit real*8 (a-h,o-z)
      dimension dens(ntri,3),fock(ntri,3)
      integer*2 xinteg2(*)
      dimension lind(*)
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
             do lcf=lcff+1,lcff+llen
               intx2=intx2+1
               xint0=xinteg2(intx2)
               ii=lind(icf)
               jj=lind(jcf)
               kk=lind(kcf)
               ll=lind(lcf)
c
               call get_dxyz(icf,kcf,ii,kk,dens,ntri, ikf,xik,yik,zik)
               call get_dxyz(jcf,kcf,jj,kk,dens,ntri, jkf,xjk,yjk,zjk)
               call get_dxyz(icf,lcf,ii,ll,dens,ntri, ilf,xil,yil,zil)
               call get_dxyz(jcf,lcf,jj,ll,dens,ntri, jlf,xjl,yjl,zjl)
c
c             ***  il  ***
c
              if( icf.gt.lcf ) then
                fock(ilf,1)=fock(ilf,1) + xjk*xint0
                fock(ilf,2)=fock(ilf,2) + yjk*xint0
                fock(ilf,3)=fock(ilf,3) + zjk*xint0
              endif
              if( lcf.gt.icf ) then
                fock(ilf,1)=fock(ilf,1) - xjk*xint0
                fock(ilf,2)=fock(ilf,2) - yjk*xint0
                fock(ilf,3)=fock(ilf,3) - zjk*xint0
              endif
c
c             ***  jl  ***
c
              if( jcf.gt.lcf ) then
                fock(jlf,1)=fock(jlf,1) + xik*xint0
                fock(jlf,2)=fock(jlf,2) + yik*xint0
                fock(jlf,3)=fock(jlf,3) + zik*xint0
              endif
              if( lcf.gt.jcf ) then
                fock(jlf,1)=fock(jlf,1) - xik*xint0
                fock(jlf,2)=fock(jlf,2) - yik*xint0
                fock(jlf,3)=fock(jlf,3) - zik*xint0
              endif
c
c             ***  ik  ***
c
              if( icf.gt.kcf ) then
                fock(ikf,1)=fock(ikf,1) + xjl*xint0
                fock(ikf,2)=fock(ikf,2) + yjl*xint0
                fock(ikf,3)=fock(ikf,3) + zjl*xint0
              endif
              if( kcf.gt.icf ) then
                fock(ikf,1)=fock(ikf,1) - xjl*xint0
                fock(ikf,2)=fock(ikf,2) - yjl*xint0
                fock(ikf,3)=fock(ikf,3) - zjl*xint0
              endif
c
c             ***  jk  ***
c
              if( jcf.gt.kcf ) then
                fock(jkf,1)=fock(jkf,1) + xil*xint0
                fock(jkf,2)=fock(jkf,2) + yil*xint0
                fock(jkf,3)=fock(jkf,3) + zil*xint0
              endif
              if( kcf.gt.jcf ) then
                fock(jkf,1)=fock(jkf,1) - xil*xint0
                fock(jkf,2)=fock(jkf,2) - yil*xint0
                fock(jkf,3)=fock(jkf,3) - zil*xint0
              endif

c
            enddo         !   over lll=1,llen
          enddo            !   over kkk=1,klen
        enddo               !   over jjj=1,jlen
      enddo                  !   over iii=1,ilen
c
      end
