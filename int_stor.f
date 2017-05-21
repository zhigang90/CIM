c=====================================================================
      subroutine int_store(ifrom,ito,bl,inx,dens,thres,labels,
     1                     isecond)
c---------------------------------------------------------------------
c Calculates integrals which will be stored :
c
c 1. in-core  (if incorex>0)   where='core'
c 2. on-disk  (if incorex=0)   where='disk'
c
c for non - and semi-direct modes.
c
c Input :
c ifrom - starting superblock
c ito   - closing superblock
c bl(*)   - calculation area
c inx(12,*) - basis set info
c dens      - fake density matrix
c thres     - integral threshold
c labels(*) - area for integral's labels
c isecond   - switch shows if we have already stored all integrals with
c             negative price
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
      integer*8 xntegrals , integ62,integ72,integx2,integx4,integx8
      common /datadisk/xntegrals,integ62,integ72,integx2,integx4,integx8
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c---------------------------------------------------------------------
      common /datadis2/nrec72,nrec62,nrecd2,nrecd4,nrecd8,npasses
c---------------------------------------------------------------------
      common /runtype/ scftype,where
c---------------------------------------------------------------------
      parameter (msplit=10000)
c---------------------------------------------------------------------
      dimension inx(12,*)
      dimension bl(*),dens(*)
      dimension labels(*)
c-----------------------------------------------------------------
      dimension ifromto(2,msplit),nofinteg(5,msplit)
c-----------------------------------------------------------------
      data iisto62_max,iisto72_max,iistor2_max /0,0,0/
      save iisto72_max,iistor2_max
      data iistor4_max,iistor8_max /0,0/
      save iistor4_max,iistor8_max
c-----------------------------------------------------------------
      data storet,storee /0.d0,0.d0/
      save storet,storee
c-----------------------------------------------------------------
      data writcpu1,writcpu2,writcpu3,writcpu4 /0.d0,0.d0,0.d0,0.d0/
      data writela1,writela2,writela3,writela4 /0.d0,0.d0,0.d0,0.d0/
c
      save writcpu1,writcpu2,writcpu3,writcpu4
      save writela1,writela2,writela3,writela4
c-----------------------------------------------------------------
ccc   call getival('iout',iout)
      iout=6
c-----------------------------------------------------------------
      call secund(time1)
      call elapsec(etim1)
c-----------------------------------------------------------------
      where ='fock'
      screen='fock'
c----------------------------------------------------------------
c set integral timings to zero as well as neglect stuff
c
c     call init_info(where)
c
c----------------------------------------------------------------
      call getival('integbl',integbl) ! pointer to max. blocksizes
      call getival('iqrtsbl',iqrtsbl) ! pointer to max. quartetsizes
c----------------------------------------------------------------
        call tstival('ftc0',iyes)
        if(iyes.eq.0) call setival('ftc0',0)
        call getival('ftc0',iftc0)
        If(iftc0.NE.0) Then
           call getival('nftcbl4',nftcbl4)  ! pointer
        endif
c----------------------------------------------------------------
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
      if(incorex.eq.0) then
c
cccc    if(npasses.eq.0) call open4store
c
        call getival('intpack',ipack) ! no of int to store
        call getival('qrtpack',iqpac) ! no of qrt to store
        call getival('blkpack',ibpac) ! no of blk to store
c
c          write(6,*)' IBPAC=',ibpac
c
        IF(IBPAC.EQ.0) go to 9876
           npasses=npasses+1
c
c
c        write(6,*)' Pass no=',npasses,' packs=',ipack,iqpac,ibpac
c
c check memory
c
c       find maximum size of a super-block; 
c       that is minimum for the integral's buffer
c
        call getival('maxsint',maxsize)
c       call setival('maxsqrt',maxsq)
c
        call getmem(0,last0)
        call retmem(1)
        call getival('lcore',lcore)
c
        memoryT=lcore-last0        ! Total free memory 
        memory23T=(2*memoryT)/3    ! 2/3 of free memory
c                                   keep 1/3 for calcul.
c
c       memory for labels   :
c
        needed4lab=(5*ibpac/2+1)/2+2  ! sizes+quarts
     *              +(2*iqpac+1)/2+1  ! starting labels
c
c       memory for integrals:
c
        needed4int=(3*maxsize+18)/4 + maxsize/2+1 + maxsize
c
        neededT=needed4lab+needed4int
c
c       write(6,*)' max=',maxsize,' mT=',memoryT
c
        if(neededT.ge.memory23T) then
           write(iout,*)'not enough memory for stored integrals+labels'
           call nerror(1,'int_stor',
     *     'not enough memory to stored integrals+labels',
     *                  neededT,memory23T)
        endif
c
c       write(6,*)'memory: Total=',memoryT,memory23T,' needed=',neededT
c
c
cccc    files opened in para_stor (open4store)
c
        iunit0=160
        iunit1=161
        iunit2=162
        iunit3=163
        iunit4=164   ! big integrals
        iunit5=165
        iunit6=166   ! small integrals 10^-7>X
        iunit7=167   ! small integrals 10^-6>X>10^-7
c
        call getint_2(4*ibpac , ijklsiz) ! I*2
        call getint_2(  ibpac , nquarts) ! I*2
        call getint_2(4*iqpac , ijkllab) ! I*2
c
        call getmem(0,last1)
        call retmem(1)
c
c
        memoryI=last1-last0
        memoryR=lcore-last1  ! free memory 
        memory23R=(2*memoryR)/3    ! 2/3 of free memory
c
c
c       There will be 5 types of integrals to store :
c                                      distribution:
c
c       10^-7 > x            stored as I*2     50%
c       10^-6 > x > 10^-7    stored as I*2     10%
c       2^15-1> x > 10^-6    stored as I*2     10%
c       2^31-1> x > 2^15-1   stored as I*4     25%
c               x > 2^31-1   stored as R*8      5%
c   
        m100=memory23R/100
c
        integbuf8=max( maxsize , 5*m100)
        integbuf4=max( maxsize ,25*m100)
        integbuf2=max( maxsize ,10*m100)
        integbu62=max( maxsize ,10*m100)
        integbu72=max( maxsize ,50*m100)
c
c
        intbuflim= 5000000      ! limit for integral buf size
c
        if(intbuflim.lt.maxsize) intbuflim=maxsize
c
c
        if(integbuf8.gt.intbuflim) integbuf8=intbuflim
        if(integbuf4.gt.intbuflim) integbuf4=intbuflim
        if(integbuf2.gt.intbuflim) integbuf2=intbuflim
        if(integbu62.gt.intbuflim) integbu62=intbuflim
        if(integbu72.gt.intbuflim) integbu72=intbuflim
c
c       check memory again
c
        neededT=(integbu72+integbu62+integbuf2+6)/4
     *            + integbuf4/2+1 + integbuf8
c
c       write(6,*)'memor2: Total=',memoryR,memory23R,' needed=',neededT
c       call f_lush(6)
c
        if(neededT.ge.memory23R) then
          write(iout,*)'not enough memory for stored integrals buffer'
          call nerror(2,'int_stor','not enough memory for int.buffer',
     *                 neededt,memory23R)
        endif
c
ctest
c       write(91,*)' INT_STOR:NPASS =',npasses,' Ibpac=',ibpac
c       write(91,*)' memoryI=',memoryI,' memoryR=',memoryR 
c       write(91,*)
c    *            ' intbu72=',integbu72,
c    *            ' intbu62=',integbu62,
c    *            ' intbuf2=',integbuf2,
c    *            ' intbuf4=',integbuf4,
c    *            ' intbuf8=',integbuf8,
c    *            '  minbuff=',maxsize
c       call f_lush(91)
ctest
c
        call getint_2(integbu72,ixin72)  ! for small integ
        call getint_2(integbu62,ixin62)  ! for small integ
        call getint_2(integbuf2,ixint2)
        call getint_4(integbuf4    ,ixint4)
        call getmem(   integbuf8    ,ixint8)
c
        call getmem(0,last2)
        call retmem(1)
c 
c       write(6,*)' Total memory allocated : alloc1=',last1-last0,
c    *            ' alloc2=',last2-last1,' Total=',last2-last0
c
      endif
c----------------------------------------------------------------
c loop over super-blocks :
c
      iistor8=0    ! count no of stored INTEGRALSin this call of int_stor
      iistor4=0    ! count no of stored INTEGRALSin this call of int_stor
      iistor2=0    ! count no of stored INTEGRALSin this call of int_stor
      iisto62=0    ! count no of stored INTEGRALSin this call of int_stor
      iisto72=0    ! count no of stored INTEGRALSin this call of int_stor
c
      iblstore=0   ! count no of stored BLOCKS   in this call of int_stor
      iqstore=0    ! count no of stored QUARTETS in this call of int_stor
c
      isplit=0
      ifromto(1,1)=1
c
      DO ISUP=IFROM,ITO
c
       if(iftc0.EQ.0) then
         isupbl=isup
       else
         call get_block4(isup,bl(nftcbl4),isupbl)
       endif
c
        call get_price(bl(ncost),isupbl,iprice) !  get integral's price
c
c
c       do not store:
c              - if price is positive and we are in the first round
c              - if price is negative and we are in the second round
c
c
        if(isecond.eq.0)then
          if(iprice.gt.0) goto 1111
        else
          if(iprice.lt.0) goto 1111
        endif
c
        call get_supbl_size(isupbl,bl(integbl),bl(iqrtsbl),isize,iqsiz)
c
        if(incorex.eq.0) then
          if( (iisto72+isize.gt.integbu72).or.
     *        (iisto62+isize.gt.integbu62).or.
     *        (iistor2+isize.gt.integbuf2).or.
     *        (iistor4+isize.gt.integbuf4).or.
     *        (iistor8+isize.gt.integbuf8)     )then
            isplit=isplit+1
            ifromto(2,isplit)=iblstore
            ifromto(1,isplit+1)=iblstore+1
            nofinteg(1,isplit)=iisto72
            nofinteg(2,isplit)=iisto62
            nofinteg(3,isplit)=iistor2
            nofinteg(4,isplit)=iistor4
            nofinteg(5,isplit)=iistor8
            iisto72_max=max(iisto72,iisto72_max)
            iisto62_max=max(iisto62,iisto62_max)
            iistor2_max=max(iistor2,iistor2_max)
            iistor4_max=max(iistor4,iistor4_max)
            iistor8_max=max(iistor8,iistor8_max)
            call secund(w1) 
            call elapsec(e1)
            if(iisto72.gt.0) then 
              call write_int2(iunit6,bl(ixin72),iisto72)
              nrec72=nrec72+1
              integ72=integ72+iisto72
c             write(6,*)' record72=',nrec72,' length=',iisto72
c             call f_lush(6)
            endif
            if(iisto62.gt.0) then 
              call write_int2(iunit7,bl(ixin62),iisto62)
              nrec62=nrec62+1
              integ62=integ62+iisto62
c             write(6,*)' record62=',nrec62,' length=',iisto62
c             call f_lush(6)
            endif
            if(iistor2.gt.0) then 
              call write_int2(iunit4,bl(ixint2),iistor2)
              nrecd2=nrecd2+1
              integx2=integx2+iistor2
c             write(6,*)' records2=',nrecd2,' length=',iistor2
c             call f_lush(6)
            endif
            if(iistor4.gt.0) then 
              call write_int4(iunit4,bl(ixint4),iistor4)
              nrecd4=nrecd4+1
              integx4=integx4+iistor4
c             write(6,*)' records4=',nrecd4,' length=',iistor4
c             call f_lush(6)
            endif
            if(iistor8.gt.0) then 
              call write_int8(iunit4,bl(ixint8),iistor8)
              nrecd8=nrecd8+1
              integx8=integx8+iistor8
c             write(6,*)' records8=',nrecd8,' length=',iistor8
c             call f_lush(6)
            endif
            call secund(w2)
            call elapsec(e2)
            writcpu4=writcpu4+w2-w1
            writela4=writela4+e2-e1
            iisto72=0
            iisto62=0
            iistor2=0
            iistor4=0
            iistor8=0
          endif
        endif
c
        iblstore=iblstore+1 
        nblstore=nblstore+1
c
c       remember number of quartets stored up to this block
c
        nqstore_before=nqstore
c
   11   continue
c
        call calcint2(isupbl,bl,inx,thres,dens,where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
c
c check if requested super-block was in a list :
c
        if(stopnow) go to 1111
c
c integrals arrived in bl(ibuffz); check if they are :
c nintez is the size of a given quartet.It is set up
c to zero if there are no integrals
c
        if(nintez.gt.0) then
          call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
          if(incorex.gt.0) then
            call core_store(
     *        bl(ibuffz),  bl(incorex), bl(ilab),    bl(jlab),
     *        bl(klab),    bl(llab),    bl(isiz),    bl(jsiz),
     *        bl(ksiz),    bl(lsiz),    bl(iqrt),    nqstore_before,
     *        integrals,   nqstore,     nblstore,    nblsiz,
     *        ngctoz,      nintez,      labels(lab1),labels(lab2),
     *        labels(lab3))
          else
            call disk_store(
     *        bl(ixin62),     bl(ixin72),    bl(ixint2),    bl(ixint4),
     *        bl(ixint8),     thres,         bl(ibuffz),    bl(ijklsiz),
     *        bl(nquarts),    bl(ijkllab),   iblstore,      iqstore,
     *        iisto62,        iisto72,       iistor2,       iistor4,
     *        iistor8,        nqstore_before,nblstore,      nqstore,
     *        xntegrals,      nblsiz,        ngctoz,        nintez,
     *        labels(lab1),   labels(lab2),  labels(lab3))
c
          endif
        endif
c---------------------------------------------------------------
        if(moreint) go to 11
c---------------------------------------------------------------
        if(nqstore.eq.nqstore_before) then
           nblstore=nblstore-1
           iblstore=iblstore-1
        endif
c---------------------------------------------------------------
c
 1111   continue
c
        if(incorex.eq.0) then
          if(isup.eq.ito) then
            isplit=isplit+1
            ifromto(2,isplit)=iblstore
            nofinteg(1,isplit)=iisto72
            nofinteg(2,isplit)=iisto62
            nofinteg(3,isplit)=iistor2
            nofinteg(4,isplit)=iistor4
            nofinteg(5,isplit)=iistor8
            iisto72_max=max(iisto72,iisto72_max)
            iisto62_max=max(iisto62,iisto62_max)
            iistor2_max=max(iistor2,iistor2_max)
            iistor4_max=max(iistor4,iistor4_max)
            iistor8_max=max(iistor8,iistor8_max)
            call secund(w1)   
            call elapsec(e1)
            if(iisto72.gt.0) then 
              call write_int2(iunit6,bl(ixin72),iisto72)
              nrec72=nrec72+1
              integ72=integ72+iisto72
c             write(6,*)' record72=',nrec72,' length=',iisto72
c             call f_lush(6)
            endif
            if(iisto62.gt.0) then 
              call write_int2(iunit7,bl(ixin62),iisto62)
              nrec62=nrec62+1
              integ62=integ62+iisto62
c             write(6,*)' record62=',nrec62,' length=',iisto62
c             call f_lush(6)
            endif
            if(iistor2.gt.0) then 
              call write_int2(iunit4,bl(ixint2),iistor2)
              nrecd2=nrecd2+1
              integx2=integx2+iistor2
c             write(6,*)' records2=',nrecd2,' length=',iistor2
c             call f_lush(6)
            endif
            if(iistor4.gt.0) then 
              call write_int4(iunit4,bl(ixint4),iistor4)
              nrecd4=nrecd4+1
              integx4=integx4+iistor4
c             write(6,*)' records4=',nrecd4,' length=',iistor4
c             call f_lush(6)
            endif
            if(iistor8.gt.0) then 
              call write_int8(iunit4,bl(ixint8),iistor8)
              nrecd8=nrecd8+1
              integx8=integx8+iistor8
c             write(6,*)' records8=',nrecd8,' length=',iistor8
c             call f_lush(6)
            endif
            call secund(w2)
            call elapsec(e2)
            writcpu4=writcpu4+w2-w1
            writela4=writela4+e2-e1
          endif
        endif
c
      ENDDO
c----------------------------------------------------------------
      nsplit=isplit
c----------------------------------------------------------------
c     write(6,*)'INT_S:secnd=',isecond,'Pass=',npasses,' Nsplit=',nsplit
c     write(6,*)'stored here :iblstore=',iblstore,
c    *                   ' iqstore=',iqstore
c     if(incorex.eq.0) then
c        write(6,*)'stored total:nblstore=',nblstore,
c    *                   ' nqstore=',nqstore,' integ=',xntegrals
c     else
c        write(6,*)'stored total:nblstore=',nblstore,
c    *                   ' nqstore=',nqstore,' integ=',integrals
c     endif
c
c     do isplit=1,nsplit
c        write(6,*)
c    * 'isplit=',isplit,' from-to=',ifromto(1,isplit),ifromto(2,isplit),
c    * ' integ=',nofinteg(isplit)
c     enddo
c----------------------------------------------------------------
      if(incorex.eq.0) then
c       write on disk ijklsiz, nquarts and ijkllab
c
        if(iblstore.gt.0) then
          call secund(w1)   
          call elapsec(e1)
          call write_33(iunit0,iblstore,iqstore,nsplit)
          call secund(w2)
          call elapsec(e2)
          writcpu1=writcpu1+w2-w1
          writela1=writela1+e2-e1
          call write_44(iunit5,ifromto,nofinteg,nsplit)
          call secund(w3)
          call elapsec(e3)
          writcpu2=writcpu2+w3-w2
          writela2=writela2+e3-e2
          call write_int2(iunit1,bl(ijklsiz),4*iblstore)
          call write_int2(iunit2,bl(nquarts),iblstore)
          call write_int2(iunit3,bl(ijkllab),4*iqstore)
          call secund(w4)
          call elapsec(e4)
          writcpu3=writcpu3+w4-w3
          writela3=writela3+e4-e3
        endif
c
        call retmem(8)
c
c
        call setival('intbu72',iisto72_max) ! needed in int_fock
        call setival('intbu62',iisto62_max) ! needed in int_fock
        call setival('intbuf2',iistor2_max) ! needed in int_fock
        call setival('intbuf4',iistor4_max) ! needed in int_fock
        call setival('intbuf8',iistor8_max) ! needed in int_fock
c
        write(91,880) iisto72_max,iisto62_max,iistor2_max,
     *                iistor4_max,iistor8_max
 880  format(
     *' Maximum record length: 72=',i8,' 62=',i8,
     *                       ' I2=',i8,' I4=',i8,' R8=',i8)
c
      endif
c----------------------------------------------------------------
      call secund(time2)
      call elapsec(etim2)
      storet=storet+(time2-time1)
      storee=storee+(etim2-etim1)
c----------------------------------------------------------------
ccc   IF(IBPAC.NE.0) return
c----------------------------------------------------------------
 9876 continue
c----------------------------------------------------------------
      write(91,*)'-----------------------------------------------------'
      write(91,*)' In the int_store routine : scf type is ',scftype
      write(91,*)'  '
c
c timing 
      if(incorex.eq.0) then
         write(91,883) npasses,xntegrals,
     *        nrec62+nrec72+nrecd2+nrecd4+nrecd8
         write(91,882) integ72,nrec72,
     *                 integ62,nrec62,
     *                 integx2,nrecd2,
     *                 integx4,nrecd4,
     *                 integx8,nrecd8
           write(91,881) writcpu1/60.d0, writela1/60.d0,
     *                   writcpu2/60.d0, writela2/60.d0,
     *                   writcpu3/60.d0, writela3/60.d0,
     *                   writcpu4/60.d0, writela4/60.d0 
c
  881   format('  Timie for writting on disk   1  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/,
     *         '                               2  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/,
     *         '                               3  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/,
     *         '                               4  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/)
c
  883   format
     *('  In ',i4,' passes ',i15,' integrals written on disk in ',i5,
     * ' records'/)
  882   format
     *('  type I*2<10^-7:',i15,' integrals written on disk in ',i5,
     * ' records'/,
     * '  type I*2<10^-6:',i15,' integrals written on disk in ',i5,
     * ' records'/,
     * '  type I*2      :',i15,' integrals written on disk in ',i5,
     * ' records'/,
     * '  type I*4      :',i15,' integrals written on disk in ',i5,
     * ' records'/,
     * '  type R*8      :',i15,' integrals written on disk in ',i5,
     * ' records'/)
c
      endif
c----------------------------------------------------------------
c
        write(91,88) storet/60.d0, storee/60.d0 
   88   format('  Timie for stored integrals    : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min')
      write(91,*)'-----------------------------------------------------'
c-----------------------------------------------------------------
      call f_lush(91)
c-----------------------------------------------------------------
c
      end
c
c=====================================================================
      subroutine core_store(
     *        buf,         core,        ilab,        jlab,
     *        klab,        llab,        isiz,        jsiz,
     *        ksiz,        lsiz,        iqrt,        nqstore_before,
     *        integrals,   nqstore,     nblstore,    nbls,
     *        ngcd,        lnijkl,      labels,      length,
     *        lgenct )
c
c----------------------------------------------------------------
c Input :
c buf() -integrals buffer
c----------------------------------------------------------------
c Store integrals in core :
c ONLY from quartets containing integrals > 0.5*thers
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(*),jlab(*),klab(*),llab(*)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)
      integer*4 iqrt(*)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension core(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)     ! NOT USED anymore
c----------------------------------------------------------------
      data half /0.5d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
c       nqrts=0
        do 100 ijklp=1,nbls
          do 150 iqu=1,ngcd
c-------------------------------------------------------
             intct=0
             xmax=0.d0
             do 200 iii=1,ilen
             do 200 jjj=1,jlen
             do 200 kkk=1,klen
             do 200 lll=1,llen
                intct=intct+1
                xint=buf(ijklp,intct,iqu)
c
c check value of integrals
c
                xabs=abs(xint)
                if(xabs.gt.xmax) xmax=xabs
c------------------------------------------------------
c
                integrals=integrals+1
                core(integrals)=xint
c
c------------------------------------------------------
  200        continue
             if(xmax.gt.half) then
                nqstore=nqstore+1
                ilab(nqstore)=labels(1,iqu,ijklp)
                jlab(nqstore)=labels(2,iqu,ijklp)
                klab(nqstore)=labels(3,iqu,ijklp)
                llab(nqstore)=labels(4,iqu,ijklp)
c               nqrts=nqrts+1
             else
c               remove integrals stored from the quartet below threshold
                integrals=integrals-intct
             endif
  150   continue
  100 continue
c----------------------------------------------------------------
c
      iqrt(nblstore)=nqstore - nqstore_before
      isiz(nblstore)=ilen
      jsiz(nblstore)=jlen
      ksiz(nblstore)=klen
      lsiz(nblstore)=llen
c
c----------------------------------------------------------------
      end
c=====================================================================
      subroutine disk_store(
     *        xinte62,        xinte72,       xinteg2,       xinteg4,  
     *        xinteg8,        thres,         buf,           ijklsiz,
     *        nquarts,        ijkllab,       iblstore,      iqstore,  
     *        iisto62,        iisto72,       iistor2,       iistor4,
     *        iistor8,        nqstore_before,nblstore,      nqstore,
     *        integrals,      nbls,          ngcd,          lnijkl,
     *        labels,         length,        lgenct)
c
c----------------------------------------------------------------
c Input :
c buf() -integrals buffer
c----------------------------------------------------------------
c Store integrals on disk :
c ONLY from quartets containing integrals > 0.5*thers
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      integer*8 integrals         ! counting of integrals stored
c
      integer*2 ijkllab(4,*)      ! iqstore
      integer*2 ijklsiz(4,*)      ! iblstore
      integer*2 nquarts(*)        ! iblstore
c
      integer*2 ilen,jlen,klen,llen,nqrt
c
      integer*2 xinte72(*)
      integer*2 xinte62(*)
      integer*2 xinteg2(*)
      integer*4 xinteg4(*)
      dimension xinteg8(*)
      dimension buf(nbls,lnijkl,ngcd)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)     ! NOT USED anymore
c----------------------------------------------------------------
      data half /0.5d0/
ccc   parameter(xi4_max=2147483647.0d0)  ! 2**(4*8-1)-1=2**31 - 1
ccc   parameter(xi4_max= 214748365.0d0)  ! above/10
      parameter(xi4_max=  21474836.0d0)  ! above/100
ccc   parameter(xi4_max=   2147483.0d0)  ! above/1000
      parameter(xi2_max=     32767.0d0)  ! 2**(2*8-1)-1=2**15 - 1
c----------------------------------------------------------------
c to select integrals < 0.5*10**-7 . 
c
      xmax7=0.5d-7/thres  ! e.g.0.5* 10**-7/10**-10 = 0.5*10**+3
      xmax6=0.5d-6/thres  ! e.g.0.5* 10**-6/10**-10 = 0.5*10**+4
c
c      write(6,*)' int_stor: xmax7=',xmax7
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
      isize=ilen*jlen*klen*llen
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        ibegi72=iisto72
        ibegi62=iisto62
        ibegin2=iistor2
        ibegin4=iistor4
        ibegin8=iistor8
        do 100 ijklp=1,nbls
          do 150 iqu=1,ngcd
c-------------------------------------------------------
             xmax=0.d0
             do intct=1,isize
                xint=buf(ijklp,intct,iqu)
                xabs=abs(xint)
                xmax=max(xmax,xabs)
             enddo
c-------------------------------------------------------
             if(xmax.gt.half) then
                if(xmax.ge.xi4_max) then
                   do intct=1,isize
                      xinteg8(ibegin8+intct)=buf(ijklp,intct,iqu)
                   enddo
                   ibegin8=ibegin8+isize
                   iistor8=iistor8+isize
                   icf1=-1
                   jcf1= 1
                   kcf1= 1
                   lcf1= 1
                else
                   if(xmax.ge.xi2_max) then
                      do intct=1,isize
                         xinteg4(ibegin4+intct)=buf(ijklp,intct,iqu)
                      enddo
                      ibegin4=ibegin4+isize
                      iistor4=iistor4+isize
                      icf1= 1
                      jcf1=-1
                      kcf1= 1
                      lcf1= 1
                   else
                      if(xmax.ge.xmax6) then  ! 1/2*  10**-6 * 10**+10
                         do intct=1,isize
                            xinteg2(ibegin2+intct)=buf(ijklp,intct,iqu)
                         enddo
                         ibegin2=ibegin2+isize
                         iistor2=iistor2+isize
                         icf1= 1
                         jcf1= 1
                         kcf1=-1
                         lcf1= 1
                      else
                         if(xmax.ge.xmax7) then  ! 1/2*  10**-7 * 10**+10
                            do intct=1,isize
                            xinte62(ibegi62+intct)=buf(ijklp,intct,iqu)
                            enddo
                            ibegi62=ibegi62+isize
                            iisto62=iisto62+isize
                            icf1= 1
                            jcf1= 1
                            kcf1= 1
                            lcf1=-1
                         else
                            do intct=1,isize
                            xinte72(ibegi72+intct)=buf(ijklp,intct,iqu)
                            enddo
                            ibegi72=ibegi72+isize
                            iisto72=iisto72+isize
                            icf1= 1
                            jcf1= 1
                            kcf1= 1
                            lcf1= 1
                         endif
                      endif
                   endif
                endif
             endif
c-------------------------------------------------------
             if(xmax.gt.half) then
                nqstore=nqstore+1
                iqstore=iqstore+1
                integrals=integrals+isize
                
                icff=1+labels(1,iqu,ijklp)
                jcff=1+labels(2,iqu,ijklp)
                kcff=1+labels(3,iqu,ijklp)
                lcff=1+labels(4,iqu,ijklp)
c
                if(icf1.eq.-1) icff=-icff
                if(jcf1.eq.-1) jcff=-jcff
                if(kcf1.eq.-1) kcff=-kcff
                if(lcf1.eq.-1) lcff=-lcff
c
                ijkllab(1,iqstore)=icff
                ijkllab(2,iqstore)=jcff
                ijkllab(3,iqstore)=kcff
                ijkllab(4,iqstore)=lcff
             endif
  150   continue
  100 continue
c----------------------------------------------------------------
      nqrt=nqstore - nqstore_before
      if(nqrt.gt.0) then
         ijklsiz(1,iblstore)=ilen
         ijklsiz(2,iblstore)=jlen
         ijklsiz(3,iblstore)=klen
         ijklsiz(4,iblstore)=llen
         nquarts(iblstore)=nqrt
      endif
c----------------------------------------------------------------
      end
c=====================================================================
      subroutine get_price(ncost,isupbl,iprice)
      dimension ncost(*)
c
      iprice=ncost(isupbl)
c
      end
c==============================================================
      subroutine write_int(iunit4,xint,isize)
      implicit real*8 (a-h,o-z)
      dimension xint(isize)
      write(iunit4) xint
      end
c==============================================================
      subroutine write_int8(iunit,xint,isize)
      implicit real*8 (a-h,o-z)
      dimension xint(isize)
cccc       write(6,*) ' integ R*8'
c          do i=1,isize
c             x=xint(i)*1.0d-10
c             write(6,66) i,x
c66   format(i5,3x,f15.8)
cccc       enddo
      write(iunit) xint
      end
c==============================================================
      subroutine write_int4(iunit,xint,isize)
      implicit real*8 (a-h,o-z)
      integer*4 xint(isize)
ccc        write(6,*) ' integ I*4'
c          do i=1,isize
c             x=xint(i)*1.0d-10
c             write(6,66) i,x
c66   format(i5,3x,f15.8)
ccc        enddo
      write(iunit) xint
      end
c==============================================================
      subroutine write_int2(iunit,xint,isize)
      implicit real*8 (a-h,o-z)
      integer*2 xint(isize)
cccc       write(6,*) ' integ I*4'
c          do i=1,isize
c             x=xint(i)*1.0d-10
c             write(6,66) i,x
c66   format(i5,3x,f15.8)
ccc        enddo
      write(iunit) xint
      end
c==============================================================
      subroutine write_33(iunit,iblstore,iqstore,nsplit) 
      integer iblstore,iqstore,nsplit
        write(iunit) iblstore,iqstore,nsplit
      end
c==============================================================
      subroutine write_44(iunit,ifromto,nofinteg,nsplit)
      integer ifromto(2,nsplit),nofinteg(5,nsplit)
c
         write(iunit) ifromto,nofinteg
c
c      write(6,*) ' Nsplit=',nsplit
c      do i=1,nsplit
c      write(6,*)
c    *' isp=',i,' i248=',nofinteg(1,i),nofinteg(2,i),nofinteg(3,i)
c      enddo
c
      end
c==============================================================
