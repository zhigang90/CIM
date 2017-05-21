c=====================================================================
c
      subroutine pstore(nblocks,bl,inx,thres,labels)
c
c=====================================================================
c  handles integral storage based on available space
c  assembles packages that can fit in (based on max. theoretical size)
c  first stores the superblocks with negative prices, then
c  start storing positives
c
c  nblocks, number of blocks
c  thres, final threshold
c  bl,inx,labels storage areas
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      integer*8 xnte, xpack, xim, ximq, xntmem, xndisk
      integer*8 xntegrals,integ62,integ72,integx2,integx4,integx8
      common /datadisk/xntegrals,integ62,integ72,integx2,integx4,integx8
      common /datstore/ thres_stored,isto,isecond,ito
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /outfile/ ioutput
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz,nqstore,iqrt,nblstore
      common /datadis2/nrec72,nrec62,nrecd2,nrecd4,nrecd8,npasses
      common /memor1a/ npard,ncost,nsupb, mxsize
      dimension inx(12,*)
      dimension bl(*)
      dimension labels(*)
c-----------------------------------------------------------------
      parameter(i4_max=2147483647)  ! 2**(4*8-1)-1=2**31 - 1
c-----------------------------------------------------------------
      call init_info('fock')
c
      call secund(time1)
c----------------------------------------------------------------
c save integral threshold used to calculate stored integrals
c (in common /datstore/ )
c
      thres_stored=thres
c
c set up starting parameters
c
      nrec72=0
      nrec62=0
      nrecd2=0
      nrecd4=0
      nrecd8=0
      npasses=0
c
      integ72=0
      integ62=0
      integx2=0
      integx4=0
      integx8=0
      xntegrals=0 ! number of stored integrals (on current slave)
c
      integrals=0 ! number of stored integrals (on current slave)
      inte=0      ! number of integrals altogether
      xnte=0      ! number of integrals altogether
      nqstore=0   ! number of quartets (on curr. slave)
      nqst=0      ! nqstore is summed up here
      nblstore=0  ! number of superblocks stored
      nblst=0     ! nblst is summed up here
      ifrom=1     ! the starting superblock
      isecond=0   ! this is 0 in the first round, 1 in the second
      iwher=0     ! slave number
c
      xc4=ncf*(ncf+1)/2
      xc4=xc4*(xc4+1.0d0)*0.5d0 ! rough estimate of all integrals
      maxq=ncs*(ncs+1)/2
      maxq=maxq*(maxq+1)/2      ! all quartets
c
      call getival('incore',incore)  !  in  Words
      call getival('indisk',indisk)  !  in MBytes
c
      call getival('intmem',intmem)   ! memory for integrals
c
c     write(6,*)' intmem=',intmem,' indisk=',indisk
c
      xntmem=0
      xndisk=0
      if(indisk.gt.0) then
         xntmem=intmem
         xndisk=indisk
         xntmem=xntmem*1000000
         xndisk=xndisk*1000000
      endif
c
c       write(6,*)' xntmem=',xntmem,' xndisk=',xndisk
c
      call getival('integbl',integbl) ! pointer to max. blocksizes
      call getival('iqrtsbl',iqrtsbl) ! pointer to max. quartetsizes
      call getival('nslv',nslv)
      nproc = nslv   ! the number of active slaves
ckw
c        write(6,*) '  ACTIVE SLAVES :'
c        do icpu=1,nslv
c          write(6,*)'     slave =',icpu
c        enddo
ckw
      if(nproc.eq.0) nproc=1
c
c do not do H2+ in parallel!
c
      if(nblocks.lt.nslv) then
         call nerror(1,'pstore ',
     1   'too many slaves for too few blocks',
     2   nslv,nblocks)
      end if
c----------------------------------------------------------------
      if(nslv.eq.0) then
         nbl4cpu=nblocks
      else
         nbl4cpu=nblocks/nslv
      endif
c----------------------------------------------------------------
        call tstival('ftc0',iyes)
        if(iyes.eq.0) call setival('ftc0',0)
        call getival('ftc0',iftc0)
        If(iftc0.NE.0) Then
           call getival('nftcbl4',nftcbl4)  ! pointer
        endif
c----------------------------------------------------------------
c make denspar(ics,jcs)=1.0
c
      call getmem(ncs*ncs,idensp)
      call setup_denschw(ncs,bl(idensp))
      call elapsec(timo)
      if(nslv.gt.0) call getrval('imbal',tlost)    ! JB jan 23 99
c----------------------------------------------------------------
      if(indisk.gt.0) call open4store
c----------------------------------------------------------------
c
c find maximum size of a block
c
      maxsi=0
      maxsq=0
      do isupbl=ifrom,nblocks
c
       if(iftc0.EQ.0) then
         isup=isupbl
       else
         call get_block4(isupbl,bl(nftcbl4),isup)
       endif
c
         call get_supbl_size(isup,bl(integbl),bl(iqrtsbl),isize,iqsiz)
         maxsi=max(maxsi,isize)
         maxsq=max(maxsq,iqsiz)
      enddo
c
c     write(6,*)' Max Integ=',maxsi,' Max Quart=',maxsq
c
      call setival('maxsint',maxsi)
c     call setival('maxsqrt',maxsq) not used MM 08/22/06
c
c----------------------------------------------------------------
 111  continue
c
c find out where the next batch will be stored and
c find out how much of the memory is used up (in integrals, common datacore)
c
      call get_storage(iwher)
c
c     write(6,*)' Iwher=',iwher,'   nblstore=',nblstore,nblst
c
      xpack = 0  !the theoretical max. size of the integralpack
      ipack = 0  !the theoretical max. size of the integralpack
      iqpac = 0  !the theoretical max. quartets in the integralpack
      ibpac = 0  ! no of stored blocks
c
c
      do isupbl=ifrom,nblocks
c
       if(iftc0.EQ.0) then
         isup=isupbl
       else
         call get_block4(isupbl,bl(nftcbl4),isup)
       endif
c
         call get_price(bl(ncost),isup,iprice)
c
c     write(6,*)' pstore: bl=',isupbl,' price=',iprice
c
c       attempt to store: - if price is negative and in the first round
c                         - if price is positive and in the second round
c
         if((iprice.lt.0.and.isecond.eq.0).or.
     $      (iprice.gt.0.and.isecond.ne.0))then
c
         call get_supbl_size(isup,bl(integbl),bl(iqrtsbl),isize,iqsiz)
c
            if(incore.gt.0) then
               ipack=ipack+isize
               iqpac=iqpac+iqsiz
            else
               xpack=xpack+isize
               iqpac=iqpac+iqsiz
            endif
            ibpac=ibpac+1
c
c if the integralpack would not fit, exit the loop
c lim is the empty storage (intmem-integrals) most of the time
c xc4 and nproc is used to make sure that slaves get proportional loads
c actually one slave can not take more than one half ?? of its share
c xc4 is only important for extremely small storage needs
c
c this xc4 is a nonsense 
c
c
           if(incore.ne.0) then
              lim =       intmem-integrals
c             limq=incore-intmem-nqstore
              limq=incorex-ilab-nqstore
cccc          if(ipack.ge.lim.or.iqpac.ge.limq) goto 113
              if(ipack.ge.lim.or.iqpac.ge.limq) then
                 ipack=ipack-isize
                 iqpac=iqpac-iqsiz
                 ibpac=ibpac-1
c                write(6,*)'  will jump to 113 : block=',isupbl
                 goto 113
              endif
              if(ibpac.gt.nbl4cpu) then
c                write(6,*)' exit because of iBpac '
                 ipack=ipack-isize
                 iqpac=iqpac-iqsiz
                 ibpac=ibpac-1
                 goto 113
              endif
           else
c             xim =       xntmem-xntegrals/2
c             ximq=xndisk-xntmem-nqstore/2
c
              xim=xntmem - (integ72+integ62+integx2)*2
     *                   -  integx4*4-integx8*8       
              ximq=xndisk-xntmem-nqstore*4
c
c             write(6,*)' int:',integ72,integ62,integx2,integx4,integx8
c    *                    , xntegrals,' limits=',xim,ximq,' bl=',isupbl
c             write(6,*)' stored : nqstore=',nqstore,' xnt=',xntegrals
c
c
c
c             write(6,*)' bl=',isupbl,' sizes=',isize,iqsiz
c             write(6,*)'    packs=',xpack,iqpac,ibpac
c             write(6,*)'    limis=',xim,ximq 
c
ccccc
              if(ibpac.gt.nbl4cpu) then
c                write(6,*)' exit because of iBpac '
                 xpack=xpack-isize
                 ipack=xpack
                 iqpac=iqpac-iqsiz
                 ibpac=ibpac-1
                 goto 113
              endif
ccccc         if(iqpac.ge.10*maxsq) then
              if(iqpac.ge.   maxsi) then
c                write(6,*)' exit because of iqpac '
                 xpack=xpack-isize
                 ipack=xpack
                 iqpac=iqpac-iqsiz
                 ibpac=ibpac-1
                 goto 113
              endif
              if(xpack.ge.i4_max) then
c                write(6,*)' exit because of i4_max '
                 xpack=xpack-isize
                 ipack=xpack
                 iqpac=iqpac-iqsiz
                 ibpac=ibpac-1
                 goto 113
              endif
              if(4.0d0*xpack.ge.xim.or.4*iqpac.ge.ximq) then
c                write(6,*)' exit because of xlim,q '
                 xpack=xpack-isize
                 ipack=xpack
                 iqpac=iqpac-iqsiz
                 ibpac=ibpac-1
                 goto 113
              endif
           endif
c
         endif
c
c     write(6,*) 'nbl isup inte pack',nblocks,isupbl,integrals,ipack
      enddo
c
 113  continue
c
      ito=isupbl-1
c
c if the integral pack is empty this processor is finished
c
      if(ito.lt.ifrom) then
         call elapsec(timb)
         tlost=tlost+(nproc-nslv)*(timb-timo)
         timo=timb
         nproc=nproc-1
         inte=inte+integrals
         xnte=xnte+xntegrals
         nqst=nqst+nqstore
         nblst=nblst+nblstore
c
         write(ioutput,210) iwher
 210     format('  Slave code: ',i9)
      endif
c----------------------------------------------------------------
c     write(6,*)' pstore 1: second=',isecond,' from-to=',ifrom,ito,
c    *    ' nproc=',nproc
c     write(6,*)' int pack=',ipack,' qrt pack=',iqpac,' blk pack=',ibpac
c
      call setival('intpack',ipack)
      call setival('qrtpack',iqpac)
      call setival('blkpack',ibpac)
c     call setival('nproces',nproc) not used MM 08/22/06
c----------------------------------------------------------------
      if(ito.GE.ifrom) then
         call para_store(ifrom,ito,bl,inx,bl(idensp),thres,labels,
     *                   iwher,isecond)
      endif
c     write(6,*)'2stored: iwher=',iwher,' incorex=',incorex,
c    $          ' nqstore=',nqstore,' xnt=',xntegrals
c----------------------------------------------------------------
c
c the next pack starts from where we left
c
      ifrom=ito+1
c
c if we have finished storing the superblocks with negative price
c we start over from the beginning and store positives
      if(ito.eq.nblocks.and.isecond.eq.0) then
         isecond=1
         ifrom=1
      endif
ckw1
c     write(6,*)'pstor2: nproc=',nproc,' sec=',isecond,' ifrom=',ifrom
ckw1
c-----------------------------------------------------------------
c check if we are finished
c
      if(nproc.gt.0) goto 111
c
c       give no more work to all slaves
c
        do iwher=1,nslv
           call para_store(ifrom,ito,bl,inx,bl(idensp),thres,labels,
     *                     iwher,isecond)
        enddo
c-----------------------------------------------------------------
      call setrval('imbal',tlost)
      call send_limits(isecond,ito)
c----------------------------------------------------------------
      if(incore.gt.0) then
         write(ioutput,*)' in-core stored integrals=',inte
      else
         write(ioutput,*)' on-disk stored integrals=',xnte
      endif
      write(ioutput,*)'                 quartets=',nqst
      write(ioutput,*)'             super-blocks=',nblst
      if(nqst.ne.0) then
        if(incore.gt.0) then
        write(ioutput,*)'     average quartet-size=',inte/nqst
        else
        write(ioutput,*)'     average quartet-size=',xnte/nqst
        endif
      end if

c calculate memory used to store inregrals & label's info :
c
      call getival('nslv',nslv)
      if(nslv.eq.0) nslv=1
      if(incore.gt.0) then
         memory_avai= incore*nslv
      else
         memory_avai= indisk*nslv
      endif

      meavai= intmem*nslv
c
      if(incore.gt.0) then
        memory_used= inte + nqst
        write(ioutput,*)' in-core memory available=',memory_avai
        write(ioutput,*)' in-core memory for int. =',meavai,
     *                  ' used=',inte
        write(ioutput,*)' in-core memory for qrt. =',memory_avai-meavai,
     *                ' used=',nqst
        write(ioutput,*)' in-core memory used     =',memory_used
      else
        memory_usei= ((integ72+integ62+integx2)*2
     *                +integx4*4+integx8*8 )/1000000
        memory_useq= ( nqst*4 )/1000000
c
        memory_used=memory_usei+memory_useq

      write(ioutput,*)' on-disk memory available=',memory_avai,'(*10^6)'
      write(ioutput,*)' on-disk memory for int. =',meavai,
     *                  ' used=',memory_usei,' (*10^6) '
      write(ioutput,*)' on-disk memory for qrt. =',memory_avai-meavai,
     *                  ' used=',memory_useq,' (*10^6) '
      write(ioutput,*)' on-disk memory used   =',memory_used,' (*10^6)'
      endif
c----------------------------------------------------------------
c release memory allocated for idensp
c
      call retmem(1)
c----------------------------------------------------------------
c print integ. timings
c
      call secund(time2)
c
      if(incore.gt.0) then
      write(ioutput,66) (time2-time1)/60.d0
      else
      write(ioutput,67) (time2-time1)/60.d0
      endif
  66  format('  in-core stored integrals calculated in :',f6.2,' min')
  67  format('  on-disk stored integrals calculated in :',f6.2,' min')
      call f_lush(ioutput)
c
      call term_info(thres,time2-time1,time2-time1,'fock')
c
      end
c
c =====================================================================
c
      subroutine get_supbl_size(isupbl,blinteg,blqrts,isize,iqsiz)
c
c =====================================================================
      implicit real*8 (a-h,o-z)
      dimension blinteg(*),blqrts(*)
c find out size of particular superblock
c
c isupbl - number of superblock
c blinteg(n of superblocks) how many integrals in a superblock,
c blqrts() how man quartets in a superblock (max.)
c isize - (output) size of particular sup.bl.
c iqsiz - (out) quartets in this supbl.
c
      isize=int(blinteg(isupbl))
      iqsiz=int(blqrts(isupbl))
      end
c =====================================================================
c
      subroutine open4store
      implicit real*8(a-h,o-z)
      character*256 scrf,filename
c----------------------------------------------------
c open sequential access files for stored integrals (FTC)
c units 160-167
c----------------------------------------------------
      nfile = 160
c
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,len1)
      len = len1 + 6
c
      filename = scrf(1:len1)//'.numbq'
      open (unit=nfile+0,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.sizes'
      open (unit=nfile+1,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.quart'
      open (unit=nfile+2,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.label'
      open (unit=nfile+3,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.integ'
      open (unit=nfile+4,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.split'
      open (unit=nfile+5,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.int72'
      open (unit=nfile+6,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      filename = scrf(1:len1)//'.int62'
      open (unit=nfile+7,file=filename(1:len),
     *      form='unformatted',status='unknown')
c
      end
c =====================================================================
      subroutine clos4int()
c
           close(160,status='delete')
           close(161,status='delete')
           close(162,status='delete')
           close(163,status='delete')
           close(164,status='delete')
           close(165,status='delete')
           close(166,status='delete')
           close(167,status='delete')
c
      end
