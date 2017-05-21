      subroutine int_fock(idft,ax,nblocks,nfock,rhf,ncf,inx,thres1,
     *                    dens,fock,fockB,dn,DenB,labels,mywork,igran)
c---------------------------------------------------------------------
c kw Sep,2000 : "density" for prescreening of integr. DENS
c                is NOT used anymore in case of UHF
c---------------------------------------------------------------------
c  KW
c This routine is called from SCF for all three modes :
c           non -, semi-, and full-direct.
c The closed-shell Fock matrix is constructed using stored integrals
c first (if any) and then re-calculated integrals. For re-calculated
c integrals variable where is set up to value of 'fock' .
c
c For non - and semi-direct mode the stored integrals are used first
c and corresponding part of Fock matrix is built .Then, for semi-direct
c mode remaining integrals are calculated (with where='fock') and
c building of the Fock matrix is finished.
c
c For recalculated integrals value of where='fock' is set up HERE
c
c PP :
c This routine calls two-el. int. block by block
c and constructs the closed-shell Fock matrix .
c
c INPUT:
c  idft      -  dft flag:  0 - no dft
c  ax        -  factor to multiply exchange contribution for ACM
c  nblocks   -  number of superblocks, not strictly needed
c  nfock     -  number of Fock matrices to contruct
c  rhf       -  logical flag, .true. if closed shell RHF
c  inx       -  common inx (contraction info)
c  thres1    -  threshold for integral*density
c  dens      -  density matrix, used to determine if a given
c               shell quartet is to be calculated
c               for open shell systems, contains alpha+beta densities
c  fock      -  alpha/closed-shell Fock matrix/matrices
c  fockB     -  beta Fock matrix
c  dn        -  alpha/closed-shell density matrix/matrices
c  DenB      -  beta density matrix
c
c  NOTE: the idea of storing both the densities and the Fock matrices
c        was abandoned because of difficulties with the parallel
c        implementation
c  labels     = memory area for labels
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NOT EXISTING super block
c---------------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
      logical rhf,moreint,stopnow
      logical rescale
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
c known from int_store   ;
      common /datstore/ thres_stored,isto,isecond,ito
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c
      common /infob/ inuc,ibas,na,nbf,nsh,NCX,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      common /cpu/ intsize,iacc,icache,memreal
      common /counters/nintsum1,nintsum2,nquartsum
      save iarray
      dimension inx(12,*)
      dimension dens(*),fock(nfock,*),dn(nfock,*)
      dimension fockB(*),DenB(*)
      dimension labels(*)
c-----------------------------------------------------------------
      data calcit,calcie /0.d0,0.d0/
      data buildt,builde /0.d0,0.d0/
      data recalt,recale /0.d0,0.d0/
c-----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
c-----------------------------------------------------------------
      parameter (Zero=0.0d0)
c-----------------------------------------------------------------
      readcpu1=0.d0
      readcpu2=0.d0
      readcpu3=0.d0
      readcpu4=0.d0
c
      readela1=0.d0
      readela2=0.d0
      readela3=0.d0
      readela4=0.d0
c-----------------------------------------------------------------
      call secund(time0)
      call elapsec(etim0)
c----------------------------------------------------------------
      where='    '
      call init_info(where)
c----------------------------------------------------------------
      write(91,*)'-----------------------------------------------------'
      write(91,*)' In the int_fock routine : scf type is ',scftype
      write(91,*) '   '
c
c  for full-direct calculations we might need to reset
c  the common datstore
c
      if(scftype.eq.'full-direct')then
        ito=0
        isecond=0
      endif
c
c----------------------------------------------------------------
c make sure the first part is only done once each scf cycle
c in parallel runs this was done for all superblock
      if(isto.eq.0) go to 9999
       isto=0
       nc22=ncf**2
       call getint(nc22,iarray)
       call setiarray(ncf,bl(iarray))
c----------------------------------------------------------------
c transform ordinary density dens(ij) into denspar(ics,jcs)
c
      call getmem(ncs*ncs,idensp)
      call getmem(ncf,map_fs)
      if(rhf) then
         call setup_densp2(inx,ncf,ncs,dens,bl(idensp),
     *                    bl(map_fs))
      else
         ntri=ncf*(ncf+1)/2
         call getmem(ntri,idenab)
         call addvec_abs(ntri,dn,denb,bl(idenab))
         call vscal(ntri, 0.5d0 , bl(idenab) )
         call setup_densp2(inx,ncf,ncs,bl(idenab),bl(idensp),
     *                    bl(map_fs))
         call retmem(1)
      endif
CCCC  call retmem(1)
c----------------------------------------------------------------
c First construct a part of fock matrix with stored integrals:
c (those integrals were calculated by calling int_store routine)
c do it only once in each SCF iteration (isto=1 means Use stored!)
c
      if(scftype.eq.'full-direct') go to 9999
c
c       return threshold used for stored integrals:
c
        thres0=thres_stored
c
c       first check where integrals are stored (in-core or on-disk):
c----------------------------------------------------------------------
        IF(INCOREX.NE.0) THEN
c         in-core storage :
          IF(idft.EQ.0) THEN
           if(rhf) then
            call fock_core(
     $        bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $        bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $        bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     bl(lind),
     $        dens,        fock)
           else
            call uhf_core(
     $        bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $        bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $        bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     bl(lind),
     $        dn,          DenB,       dens,       fock,
     $        fockB)
           endif
          ELSE IF(ax.NE.Zero) THEN
           if(rhf) then
            call acmfock_core(ax,bl(incorex),
     *                     bl(ilab),bl(jlab),bl(klab),bl(llab),
     *                     bl(isiz),bl(jsiz),bl(ksiz),bl(lsiz),bl(iqrt),
     *                     integrals,nqstore,nblstore,
     *                     thres0,thres1, bl(lind),dens,fock)
           else
            call acmuhf_core(ax,bl(incorex),
     *                    bl(ilab),bl(jlab),bl(klab),bl(llab),
     *                    bl(isiz),bl(jsiz),bl(ksiz),bl(lsiz),bl(iqrt),
     *                    integrals,nqstore,nblstore,
     *                    thres0,thres1,
     *                    bl(lind),dn,DenB,dens,fock,fockB)
           endif
          ELSE
           if(rhf) then
            call dftfock_core(bl(incorex),
     *                     bl(ilab),bl(jlab),bl(klab),bl(llab),
     *                     bl(isiz),bl(jsiz),bl(ksiz),bl(lsiz),bl(iqrt),
     *                     integrals,nqstore,nblstore,
     *                     thres0,thres1, bl(lind),dens,fock)
           else
            call dftuhf_core(bl(incorex),
     *                    bl(ilab),bl(jlab),bl(klab),bl(llab),
     *                    bl(isiz),bl(jsiz),bl(ksiz),bl(lsiz),bl(iqrt),
     *                    integrals,nqstore,nblstore,
     *                    thres0,thres1,
     *                    bl(lind),dens,fock,fockB)
           endif
          ENDIF
c----------------------------------------------------------------------
        ELSE
c----------------------------------------------------------------------
c          on-disk storage :
c
c          call open4store
c
           rewind(160)
           rewind(161)
           rewind(162)
           rewind(163)
           rewind(164)
           rewind(165)
           rewind(166)
           rewind(167)
c
c       write(6,*)' INT_FOCK : NBlstore=',nblstore,' NQstore=',nqstore
c       call f_lush(6)
c
           nblcount=0
           rescale=.false.
 9876      continue
c
           call secund(r1)
           call elapsec(e1)
c
           call read_33(160,iblstore,iqstore,nsplit) 
c
c          write(6,*)' iblstore,iqstore,nsplit=',
c    *                 iblstore,iqstore,nsplit
c
           if(iblstore.eq.0) go to 9999
c
           call secund(r2)
           call elapsec(e2)
           readcpu1=readcpu1+(r2-r1)
           readela1=readela1+(e2-e1)
c
           call getint(nsplit*2,ifromto)
           call getint(nsplit*5,nofinteg)
c
           call secund(r1)
           call elapsec(e1)
c
           call read_44(165,bl(ifromto),bl(nofinteg),nsplit)
c
           call secund(r2)
           call elapsec(e2)
           readcpu2=readcpu2+(r2-r1)
           readela2=readela2+(e2-e1)
c
           nblcount=nblcount+iblstore
           if(nblcount.EQ.nblstore) rescale=.true.
c
c       write(6,*)' this pass: IBlstore=',iblstore,' IQstore=',iqstore,
c    *       ' rescale=',rescale,' Nsplit=',nsplit
c
           call getint_2(4*iblstore , ijklsiz) ! I*2
           call getint_2(  iblstore , nquarts) ! I*2
           call getint_2(4*iqstore  , ijkllab) ! I*2
c
           call getival('intbu72',integbu72)
           call getival('intbu62',integbu62)
           call getival('intbuf2',integbuf2)
           call getival('intbuf4',integbuf4)
           call getival('intbuf8',integbuf8)
c
c        write(6,*)' int buff=',integbu72,integbuf2,integbuf4,integbuf8
c
           call getmem(integbuf8,ixint8)
           call getint_4(integbuf4,ixint4)
           call getint_2(integbuf2,ixint2)
           call getint_2(integbu62,ixin62)
           call getint_2(integbu72,ixin72)
c
c----------------Fock-builders with on-disk integrals------------------
c
          IF(idft.EQ.0) THEN
             if(rhf) then
                call  fock_disk(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   dens,       fock,
     $        ncs,         bl(idensp), bl(map_fs))
             else
                call  uhf_disk(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   dn,         denb,
     $        dens,        fock,       fockb,      ncs,
     $        bl(idensp),  bl(map_fs))
             endif
          ELSE IF(ax.NE.Zero) THEN
             if(rhf) then
                call  acmfock_disk(
     $        ax,          ncf,        bl(ixin72), bl(ixin62),
     $        bl(ixint2),  bl(ixint4), bl(ixint8), bl(ijkllab),
     $        bl(ijklsiz), bl(nquarts),iqstore,    iblstore,
     $        rescale,     nsplit,     bl(ifromto),bl(nofinteg),
     $        thres0,      thres1,     bl(lind),   dens,
     $        fock,        ncs,        bl(idensp), bl(map_fs))
             else
                call  acmuhf_disk(
     $        ax,          ncf,         bl(ixin72), bl(ixin62),
     $        bl(ixint2),  bl(ixint4),  bl(ixint8), bl(ijkllab),
     $        bl(ijklsiz), bl(nquarts), iqstore,    iblstore,
     $        rescale,     nsplit,      bl(ifromto),bl(nofinteg),
     $        thres0,      thres1,      bl(lind),   dn,
     $        denb,        dens,        fock,       fockb,
     $        ncs,         bl(idensp),  bl(map_fs))
             endif
          ELSE
             if(rhf) then
                call dftfock_disk(
     $        ncf,         bl(ixin72),  bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8),  bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,     iblstore,   rescale,
     $        nsplit,      bl(ifromto), bl(nofinteg),thres0,
     $        thres1,      bl(lind),    dens,       fock,
     $        ncs,         bl(idensp),  bl(map_fs))
             else
                call dftuhf_disk(
     $        ncf,         bl(ixin72),  bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8),  bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,     iblstore,   rescale,
     $        nsplit,      bl(ifromto), bl(nofinteg),thres0,
     $        thres1,      bl(lind),    dens,       fock,
     $        fockb,       ncs,         bl(idensp), bl(map_fs))
             endif
          ENDIF
c
c----------------Fock-builders with on-disk integrals--end-------------

           call retmem(10)
           if(.not.rescale) go to 9876
        ENDIF
c----------------------------------------------------------------------
c
        if(scftype.eq.'non -direct') then
           call secund(time1)
           call term_info(thres1, 0.d0,time1-time0,where)
           return
        endif
c
        call secund(time1)
        call elapsec(etim1)
c
        buildt=buildt+(time1-time0)
        builde=builde+(etim1-etim0)
c
        if(incorex.ne.0) then
           write(91,87) buildt/60.d0, builde/60.d0 
        else
           write(91,88) buildt/60.d0, builde/60.d0 
c
           write(91,881) readcpu1/60.d0, readela1/60.d0,
     *                   readcpu2/60.d0, readela2/60.d0,
     *                   readcpu3/60.d0, readela3/60.d0,
     *                   readcpu4/60.d0, readela4/60.d0 
  881   format('  Timie for reading fromdisk   1  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/,
     *         '                               2  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/,
     *         '                               3  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/,
     *         '                               4  : cpu=',f10.2,
     *  ' elapsed=',f10.2,' min',/)
c
        endif
   87   format('  Timie for in-core Fock builder: cpu=',f10.2,
     *  ' elapsed=',f10.2,' min')
   88   format('  Timie for on-disk Fock builder: cpu=',f10.2,
     *  ' elapsed=',f10.2,' min')
        write(91,*) ' no of stored       blocks=',nblstore
      write(91,*)'-----------------------------------------------------'
c----------------------------------------------------------------
 9999 continue
c
c  check if we actually have to compute any integrals,
c  because it might happen that a slave executes this routine
c  even with mywork=-1, in order to retrieve the stored integrals
c  (this is to avoid a race condition MM 08/29/2006 )
c
      if(mywork.lt.0)return
c
      call mmark
c----------------------------------------------------------------
      call secund(time1)
      call elapsec(etim1)
c----------------------------------------------------------------
c re-calculate integrals
c
      where='fock'
      screen='fock'
      if(idft.GT.0.AND.ax.EQ.Zero) screen='coul'    ! ** NEW  JB **
c----------------------------------------------------------------
c set integral timings to zero as well as neglect stuff
c
      call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
        call tstival('ftc0',iyes)
        if(iyes.eq.0) call setival('ftc0',0)
        call getival('ftc0',iftc0)
        If(iftc0.NE.0) Then
           call getival('nftcbl4',nftcbl4)  ! pointer
        endif
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
      irecalc=0
cccc  DO isupbl=istart,istop
      DO isup  =istart,istop
c
         if(iftc0.EQ.0) then
           isupbl=isup
         else
           call get_block4(isup,bl(nftcbl4),isupbl)
         endif
c
c        get integral's price
         call get_price(bl(ncost),isupbl,iprice)
c blocks were stored:
c  - if they have negative price and they are before the last stored
c  - if they have negative price and all negative blocks were stored
c  - if all negative blocks were stored and our block is before the last
         if((iprice.LE.0.and.isup.le.ito).or.(iprice.le.0.and.
     *   isecond.eq.1).or.(isecond.eq.1.and.isup.le.ito)) go to 1111
c
         irecalc=irecalc+1
c
   11 continue
c
        call secund(timc1)
        call elapsec(etic1)
        call calcint2(isupbl,bl,inx,thres1,bl(idensp),where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
        call secund(timc2)
        call elapsec(etic2)
       calcit=calcit+timc2-timc1
       calcie=calcie+etic2-etic1
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
           lsh=labels(lab2)*labels(lab2+1)*labels(lab2+2)*labels(lab2+3)
           nintsum2=nintsum2+nblsiz*lsh
           nquartsum=nquartsum+nblsiz
c   these are sums: nintsum2 sum up nblsiz*lsh and nquartsum nblsiz
c   they must be zeroed in the SCF program
         IF(idft.EQ.0) THEN
c  ...........................
c   Standard HF
c  ...........................
           If(lsh.gt.1296) Then
            if(rhf) then
             call fock_bldr(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
            else
             call uhf_bldr(ncf,bl(ibuffz),fock,fockB,dn,DenB,
     $                   dens,bl(iarray),nblsiz,ngctoz,nintez,
     $                   labels(lab1),labels(lab2),labels(lab3))
            endif
           Else
            if(labels(lab2+3).gt.1) then
              if(nfock.eq.1) then
                if(rhf) then
                 call fock_bldr1(ncf,bl(ibuffz),fock,dn,
     *                     bl(iarray),nblsiz,ngctoz,nintez,
     *                     labels(lab1),labels(lab2),labels(lab3))
                else
                 call uhf_bldr1(ncf,bl(ibuffz),fock,fockB,dn,DenB,
     $                        dens,bl(iarray),nblsiz,ngctoz,nintez,
     $                        labels(lab1),labels(lab2),labels(lab3))
                endif
              else
                call fock_bldr2(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              end if
            else
                if(nfock.gt.1) then
                call fock_bldr3(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
                else
                if(rhf) then
                 call fock_bldr4(ncf,bl(ibuffz),fock,dn,
     *                     bl(iarray),nblsiz,ngctoz,nintez,
     *                     labels(lab1),labels(lab2),labels(lab3))
                else
                 call uhf_bldr4(ncf,bl(ibuffz),fock,fockB,dn,DenB,
     $                        dens,bl(iarray),nblsiz,ngctoz,nintez,
     $                        labels(lab1),labels(lab2),labels(lab3))
                endif
              end if
              end if
           EndIf
         ELSE IF(ax.NE.Zero) THEN
c  ...........................
c   Hybrid HF-DFT (ACM)
c  ...........................
           If(lsh.gt.1296) Then
            If(rhf) then
             call acmfock_bldr(ax,nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
            else
             call acmuhf_bldr(ax,ncf,bl(ibuffz),fock,fockB,dn,DenB,
     $                      dens,bl(iarray),nblsiz,ngctoz,nintez,
     $                      labels(lab1),labels(lab2),labels(lab3))
            endif
           Else
            if(labels(lab2+3).gt.1) then
              if(nfock.eq.1) then
               if(rhf) then
                call acmfock_bldr1(ax,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
               else
                call acmuhf_bldr1(ax,ncf,bl(ibuffz),fock,fockB,dn,DenB,
     $                          dens,bl(iarray),nblsiz,ngctoz,nintez,
     $                          labels(lab1),labels(lab2),labels(lab3))
               endif
              else
                call acmfock_bldr2(ax,nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              end if
            else
              if(rhf) then
               call acmfock_bldr3(ax,nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              else
               call acmuhf_bldr3(ax,ncf,bl(ibuffz),fock,fockB,dn,DenB,
     $                    dens,bl(iarray),nblsiz,ngctoz,nintez,
     $                    labels(lab1),labels(lab2),labels(lab3))
              endif
             end if
           EndIf
         ELSE
c  ...........................
c   "Pure" DFT (NO EXCHANGE)
c  ...........................
           If(lsh.gt.1296) Then
            if(rhf) then
             call dftfock_bldr(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
            else
             call dftuhf_bldr(ncf,bl(ibuffz),fock,fockB,dens,
     $                      bl(iarray),nblsiz,ngctoz,nintez,
     $                      labels(lab1),labels(lab2),labels(lab3))
            endif
           Else
             if(labels(lab2+3).gt.1) then
               if(nfock.eq.1) then
                if(rhf) then
                 call dftfock_bldr1(ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
                else
                 call dftuhf_bldr1(ncf,bl(ibuffz),fock,fockB,dens,
     $                           bl(iarray),nblsiz,ngctoz,nintez,
     $                           labels(lab1),labels(lab2),labels(lab3))
                endif
               else
                 call dftfock_bldr2(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
               end if
             else
              if(rhf) then
               call dftfock_bldr3(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              else
               call dftuhf_bldr3(ncf,bl(ibuffz),fock,fockB,dens,
     $                           bl(iarray),nblsiz,ngctoz,nintez,
     $                           labels(lab1),labels(lab2),labels(lab3))
              endif
             end if
           EndIf
         ENDIF
        endif
c
c----------------------------------------------------------------
c
        if(moreint) go to 11
c
 1111 continue
c
      END DO
      call retmark
c----------------------------------------------------------------
ccccc   call dens_fock_prt(dens,fock,'total-fock ')
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
c  DO NOT rescale the fock matrix here - it is made in para_fock
c     call resc_nfock(ncf,nfock,fock,thres1)
c
      call secund(time2)
      call elapsec(etim2)
      call term_info(thres1,time2-time1,time1-time0,where)
c
      recalt=recalt+time2-time1
      recale=recale+etim2-etim1
c
        write(91,89) recalt/60.d0, recale/60.d0 
        write(91,90) calcit/60.d0, calcie/60.d0
   89   format('  Timie for rec.int+Fock builder: cpu=',f10.2,
     *  ' elapsed=',f10.2,' min')
   90   format('  Timie for recalculating inetg.: cpu=',f10.2,
     *  ' elapsed=',f10.2,' min')
        write(91,*) ' no of recalculated blocks=',irecalc
      write(91,*)'-----------------------------------------------------'
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine fock_core(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        dens,        fock)
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
cccc  integer*2 isiz(nqstore),jsiz(nqstore),ksiz(nqstore),lsiz(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      dimension core(*)       !      dimension core(integrals)
      dimension dens(*),fock(*)
      dimension lind(*)
      dimension iix(4)
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      intx=0
      ijklp=0
      do isbl=1,nblstore
         ilen=isiz(isbl)
         jlen=jsiz(isbl)
         klen=ksiz(isbl)
         llen=lsiz(isbl)
c
c     write(6,*)' fock_core: isbl=',isbl,' sizes=',ilen,jlen,klen,llen
c
         nqrt1=iqrt(isbl)
c
c     write(6,*)' fock_core: isbl=',isbl,' sizes=',ilen,jlen,klen,llen,
c    *  ' quarts=',nqrt1
c
         do iqrt1=1,nqrt1
            ijklp=ijklp+1
            icff=ilab(ijklp)
            jcff=jlab(ijklp)
            kcff=klab(ijklp)
            lcff=llab(ijklp)
c
ctest
c        write(6,*)' ijklp=',ijklp
c        write(6,*)'i-l,len=',ilen,jlen,klen,llen
c        write(6,*)'i-l,cff=',icff,jcff,kcff,lcff
c
         do iii=1,ilen
            icf=icff+iii
            do jjj=1,jlen
               jcf=jcff+jjj
               do kkk=1,klen
                  kcf=kcff+kkk
                  do lll=1,llen
                     intx=intx+1
                     xint=core(intx)
                     xin4=4.d0*xint
                     lcf=lcff+lll
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
ctest.......................................................
c               write(6,55) icf,jcf,kcf,lcf,xint*1.0d-10
 55              format(4i3,f15.9,2x,'fock-core')
ctest.......................................................
c
                     ii=lind(iff)
                     jj=lind(jff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
                     ik=ii+kff
                     il=ii+lff
c
                     if(jff.ge.kff) then
                       jk=jj+kff
                     else
                       jk=kk+jff
                     end if
                     if(jff.ge.lff) then
                       jl=jj+lff
                     else
                       ll=lind(lff)
                       jl=ll+jff
                     end if
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
                     fock(ik)=fock(ik)-xint*dens(jl)
                     fock(jl)=fock(jl)-xint*dens(ik)
                     fock(jk)=fock(jk)-xint*dens(il)
                     fock(il)=fock(il)-xint*dens(jk)
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
         enddo                  !   do iqrt1=1,nqrt1
      enddo                     !   do isbl=1,nblstore
c
      if(thres0.ne.thres1) then
         resc=thres0/thres1
         call dscal(ncf*(ncf+1)/2,resc,fock,1)
      endif
c
ctest
          if(integrals.NE.intx) then
             call nerror(1,'fock_core',
     $           'integrals .ne. intx',integrals,intx)
          endif
ctest
c
      end
c=====================================================================
      subroutine fock_disk(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dens,       fock,
     $        ncs,         densp,      map_fs)
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
      dimension dens(*),fock(*)
c
      dimension lind(*)
c----------------------------------------------------------------
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
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
c if screenig density densp() is not used then 
c
           dmax1=1.0d0
           dmax=dmax1*t01
c----------------------------------------------------------------
c
           call secund(r1)
           call elapsec(e1)
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
           call secund(r2)
           call elapsec(e2)
           readcpu3=readcpu3+(r2-r1)
           readela3=readela3+(e2-e1)
c----------------------------------------------------------------
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
c     write(6,*)' integ=',integi72,integin2,integin4,integin8
c     call f_lush(6)
c
           call secund(r1)
           call elapsec(e1)
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
           call secund(r2)
           call elapsec(e2)
           readcpu4=readcpu4+(r2-r1)
           readela4=readela4+(e2-e1)
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
c
c        write(6,*)' bl=',isbl,' size=',isize,' nq=',nqrt
c
         do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
c
c        write(6,*)' bl=',isbl,' size=',isize,' nq=',nqrt
c        write(6,*)'         qr=',iqrt,' labs=',icff,jcff,kcff,lcff
c        call f_lush(6)
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
c--------------------------------------------------------
c
c        write(6,*)' case=',icase,' labs=',icff,jcff,kcff,lcff
c        write(6,*)' case=',icase,' lens=',ilen,jlen,klen,llen
c             call f_lush(6)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c shell :
c          ics=map_fs(icff+1)
c          jcs=map_fs(jcff+1)
c          kcs=map_fs(kcff+1)
c          lcs=map_fs(lcff+1)
c
c screening density
c
c          dij=densp(ics,jcs)
c          dkl=densp(kcs,lcs)
c          dik=densp(ics,kcs)
c          dil=densp(ics,lcs)
c          djk=densp(jcs,kcs)
c          djl=densp(jcs,lcs)
c
c          dmax1=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c          dmax1=1.0d0
c          dmax=dmax1*t01
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         using  dmax*abs(x) > 0.5 is the same as 
c
c                     abs(x) > 0.5/dmax 
c
c         and no multiplications below 
c
              dmax2=half/dmax
c
c
              if(integin8.gt.intx8 .and. icase.eq.1) then
                 call make_fock8(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          intx8,xinteg8)
              endif
c
              if(integin4.gt.intx4 .and. icase.eq.2) then
                 call make_fock4(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          intx4,xinteg4)
              endif
c
              if(integin2.gt.intx2 .and. icase.eq.3) then
                 call make_fock2(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          intx2,xinteg2)
              endif
c
              if(doint62) then
               if(integi62.gt.int62 .and. icase.eq.4) then
                 if(dmax1*1.0d-6.gt.thres1) then
                 call make_fock2(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          int62,xinte62)
                 endif
               endif
              endif
c
              if(doint72) then
               if(integi72.gt.int72 .and. icase.eq.5) then
                 if(dmax1*1.0d-7.gt.thres1) then
                 call make_fock2(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          int72,xinte72)
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
         call dscal(ncf*(ncf+1)/2,resc,fock,1)
      endif
c
      end
c=====================================================================
c=====================================================================
c initiate and terminate routines :
c
      subroutine init_info(where)
      implicit real*8 (a-h,o-z)
      character*4 where
ctest only
c     common /CONUT_JUMPS/
c    * icond_1,icond_2,icond_3,icond_4,icond_5,icond_6,icond_7,icond_8
ctest only
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /stabili/ nstable,nostable
ctest only
c     icond_1=0
c     icond_2=0
c     icond_3=0
c     icond_4=0
c     icond_5=0
c     icond_6=0
c     icond_7=0
c     icond_8=0
ctest only
c
         time0=0.d0
         call timepr(time0)
c
         nstable=0
         nostable=0
c
      if(where.eq.'fock') then
         timreca2=0.d0
      else
         timstor2=0.d0
      endif
c
         lpartot=0
         lpareal=0
c
         ntotal=0
         noijkl=0
         nopres=0
         nohalf=0
c
         nrimtot=0
         nrimret=0
c
      end
c-------------------------
      subroutine term_info(thres,time12,time01,where)
      implicit real*8 (a-h,o-z)
      character*4 where
ctest only
c     common /CONUT_JUMPS/
c    * icond_1,icond_2,icond_3,icond_4,icond_5,icond_6,icond_7,icond_8
ctest only
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                          isiz,jsiz,ksiz,lsiz, nqstore
      common /outfile/ ioutput
      common /stabili/ nstable,nostable
c
      if(where.eq.'    ') then
         write(ioutput,130)
c        only non-direct mode; time for fock building time01
         write(ioutput,145) timprep2+timstor2
         write(ioutput,150) timprep2+timstor2 + time01
         return
      endif
c
c--------------------------------------------
c percent of retained quartets :
c
      if(ntotal.le.0 .or. nrimtot.le.0) RETURN
c
      pcshq=dble(nopres)/dble(ntotal)*100.d0
      ppshq=dble(nrimret)/dble(nrimtot)*100.d0
c
      write(ioutput,130)
      write(ioutput,1001) thres
      write(ioutput,1002) ntotal,nopres,pcshq,nrimtot,nrimret,ppshq
 1001 format('screening of shell quartets with threshold=',1pd15.4)
 1002 format(
     * '  contracted level :  total=',i10,' retained=',i10,' ',f5.1,'%'/
     * '  primitive  level :  total=',i10,' retained=',i10,' ',f5.1,'%')
c--------------------------------------------
c satble & unstable primitive quartets :
c
      pstab=dble(nostable)/dble(nstable+nostable)*100.d0
      write(ioutput,1003) nstable,nostable,pstab
 1003 format(
     * '        numerically stable =',i10,' unstable=',i10,' ',f5.1,'%')
c--------------------------------------------
c
c print subroutine timing :
c
c     call timepr(timblok2 + 0.000001d0)
c
      RETURN
c--------------------------------------------
      if(where.eq.'core') then
         write(ioutput,*)
     *   ' number of integrals stored in/on ',where,' =',integrals
         write(ioutput,135) timstor2
         write(ioutput,136) time12
      endif
      if(where.eq.'disk') then
         write(ioutput,*)
     *   ' number of integrals stored in/on ',where,' =',integrals
         write(ioutput,135) timstor2
         write(ioutput,137) time12
      endif
      if(where.eq.'fock') then
ccccc    write(ioutput,140) timreca2
ccccc    write(ioutput,138) time12
ccccc    write(ioutput,145) timprep2+timstor2+timreca2
         write(ioutput,150) timprep2+timstor2+          time12+time01
      endif
      if(where.eq.'forc') then
c???     write(ioutput,245) timprep2+timstor2+timreca2
c???     write(ioutput,245) timprep2+timstor2
         write(ioutput,250) timprep2+timstor2+          time12+time01
      endif
c
  130 format('======================================================')
  135 format(' Total time for stored two-el.integrals = ',f10.2,' sec' )
  136 format(' Total time for two-el.ineg. + in-core  = ',f10.2,' sec'/)
  137 format(' Total time for two-el.ineg. + on-disk  = ',f10.2,' sec'/)
  138 format(' Total time for two-el.ineg. + in-fock  = ',f10.2,' sec'/)
c
  140 format(' Total time for recalc two-el.integrals = ',f10.2,' sec' )
  145 format(' Total time for two-electron integrals  = ',f10.2,' sec' )
  150 format(' Total time for two-el.ineg. + fockmat  = ',f10.2,' sec'/)
c
  245 format(' Total time for Ist-deriv.two-el.integ. = ',f10.2,' sec' )
  250 format(' Total time for two-el.forces on atoms  = ',f10.2,' sec'/)
c-----------------------------------------------------------------------
c
      i9=ioutput
c97   write(i9,*)'-----------------------------------------------------'
c97   write(i9,*) 'Pre-calculations for pairs of contracted shells have'
c97   write(i9,*) 'been executed ',lpareal,'(out of ',lpartot,') times '
c
      write(i9,*)'-----------------------------------------------------'
      write(i9,*)'Total number of contracted shell quartets ',ntotal
      write(i9,*)'has been reduced by molecular symmetry to ',noijkl
      write(i9,*)'and due to the screening procedure to     ',nopres
      write(i9,*)'-- - - - - - - - - - - - - - - - - - - - - - - - - - '
      write(i9,*)'number of primitive quartets after symmet.',nrimtot
      write(i9,*)'reduced by the screening procedure to     ',nrimret
      write(i9,*)'-----------------------------------------------------'
c
      write(ioutput,130)
c-----------------------------------------------------------------------
      end
c=====================================================================
      subroutine get_lab123(nbls,ngctot,lab1,lab2,lab3)
c
        lab1=1
        lab2=lab1+4*nbls*ngctot
        lab3=lab2+4
c
      end
c=====================================================================
      subroutine dens_fock_prt(dens,fock,text)
      implicit real*8 (a-h,o-z)
      character*11 text
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      dimension dens(*),fock(*)
c
      write(6,*) text
c
      write(6,*)' fock   '
      do 100 i=1,ncf
      ii=i*(i-1)/2
      write(6,88)(1.d-10*fock(ii+j),j=1,i)
  100 continue
c
      write(6,*)' density'
      do 200 i=1,ncf
      ii=i*(i-1)/2
      write(6,88)(dens(ii+j),j=1,i)
  200 continue
   88 format(13(f9.4,1x))
      end
c====================================================================
      subroutine setiarray(ncf,iarray)
      dimension iarray(ncf,ncf)
      ij=0
      do i=1,ncf
        do j=1,i
          ij=ij+1
          iarray(i,j)=ij
          iarray(j,i)=ij
        end do
      end do
      end
c===================================================
      subroutine resc_nfock(ncf,nfock,fock,thres)
      implicit real*8 (a-h,o-z)
      dimension fock(nfock,*)
      data half /0.5d0/
c------------------------------------------------------
c divide the non-diagonal elements of the Fock matrix by 2
c and scale back to the normal values from the large ones
c-------------------------------------------------------
      thinv2=half*thres
      ij=0
      do icf=1,ncf
        do jcf=1,icf
          ij=ij+1
          do ifock=1,nfock
            fock(ifock,ij)=fock(ifock,ij)*thinv2
          if(icf.eq.jcf) fock(ifock,ij)=fock(ifock,ij)+fock(ifock,ij)
          end do
        end do
      end do
c-------------------------------------------------------
      end
c=====================================================================
      subroutine resu_nfock(ncf,fockA,fockB,thres)
      implicit real*8 (a-h,o-z)
      dimension fockA(*),fockB(*)
c------------------------------------------------------
c for open shell need to double all Fock matrix elements
c so, compared to closed-shell, need to double only the
c diagonal elements (see <resc_nfock> above)
c-------------------------------------------------------
      ij=0
      do 20 icf=1,ncf
      do 10 jcf=1,icf
      ij=ij+1
      fockA(ij)=fockA(ij)*thres
      fockB(ij)=fockB(ij)*thres
 10   continue
c ....................................................
c -- WARNING! next 2 lines depend on ij retaining
c --          its final value in loop 10
      fockA(ij) = fockA(ij) + fockA(ij)
      fockB(ij) = fockB(ij) + fockB(ij)
 20   continue
c-------------------------------------------------------
      end
c=====================================================================
      subroutine setup_densp(inx,lind,ncf,ncs,densmat,denspar,map_fs)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /screening/ dens_max_el      ! output
      dimension inx(12,ncs)
      dimension lind(*)
      dimension densmat(*)
      dimension map_fs(ncf)
      dimension denspar(ncs,ncs)          ! output
      data zero /0.d0/
c--------------------------------------------------------------------
c densmat is assumed to be ncf*(ncf+1)/2 in all cases
c--------------------------------------------------------------------
c input :
c inx(12,ncs)       = basis set info
c lind(i)=i*(i-1)/2 = pre-computed in prepint2
c ncf, ncs          = contracted functions and shells
c densmat           = full density matrix
c
c output :
c densp(ics,jcs)    = max absolute elemet of density for each ics,jcs
c                     it is SQUARED for ordinary scf integrals
c dens_max_el       = maximum element overall
c                     save in common /screening/
c
c local
c map_fs(icf)=ics     mapping from contr.funct to contr.shells
c--------------------------------------------------------------------
c
      do 10 ics=1,ncs
      icf_b=inx(11,ics)+1
      icf_e=inx(10,ics)
      do 10 icf=icf_b,icf_e
      map_fs(icf)=ics
   10 continue
c
      call zeroit(denspar(1,1),ncs*ncs)
c
      dens_max_el=zero
      ijcf=0
      do 100 icf=1,ncf
      ics=map_fs(icf)
      do 100 jcf=1,icf
      jcs=map_fs(jcf)
c
      ijcf=ijcf+1
      dijcf=abs(densmat(ijcf))
      dijcs=denspar(ics,jcs)
c
      if(dijcf.gt.dijcs) then
         denspar(ics,jcs)=dijcf
         denspar(jcs,ics)=dijcf
         if(dijcf.gt.dens_max_el) dens_max_el=dijcf
      endif
  100 continue
c
      dens_max_el=dens_max_el*dens_max_el
c
      if(where.eq.'forc' .or. where.eq.'hess') return
c
c square it for scf integrals
c
      do 200 ics=1,ncs
      do 200 jcs=1,ics
      dmx=denspar(ics,jcs)
      dmx2=dmx*dmx
      denspar(ics,jcs)=dmx2
      denspar(jcs,ics)=dmx2
  200 continue
c
      end
c====================================================================
      subroutine setup_densp2(inx,ncf,ncs,densmat,denspar,map_fs)
      implicit real*8 (a-h,o-z)
      common /screening/ dens_max_el      ! output
      dimension inx(12,ncs)
      dimension densmat(*)
      dimension map_fs(ncf)
      dimension denspar(ncs,ncs)          ! output
      data zero /0.d0/
c--------------------------------------------------------------------
c densmat is assumed to be ncf*(ncf+1)/2 in all cases
c--------------------------------------------------------------------
c input :
c inx(12,ncs)       = basis set info
c ncf, ncs          = contracted functions and shells
c densmat           = full density matrix
c
c output :
c densp(ics,jcs)    = max absolute elemet of density for each ics,jcs
c dens_max_el       = maximum element overall
c                     save in common /screening/
c
c local
c map_fs(icf)=ics     mapping from contr.funct to contr.shells
c--------------------------------------------------------------------
c
      do 10 ics=1,ncs
      icf_b=inx(11,ics)+1
      icf_e=inx(10,ics)
      do 10 icf=icf_b,icf_e
      map_fs(icf)=ics
   10 continue
c
      call zeroit(denspar(1,1),ncs*ncs)
c
      dens_max_el=zero
      ijcf=0
      do 100 icf=1,ncf
      ics=map_fs(icf)
      do 100 jcf=1,icf
      jcs=map_fs(jcf)
c
      ijcf=ijcf+1
      dijcf=abs(densmat(ijcf))
      dijcs=denspar(ics,jcs)
c
      if(dijcf.gt.dijcs) then
         denspar(ics,jcs)=dijcf
         denspar(jcs,ics)=dijcf
         if(dijcf.gt.dens_max_el) dens_max_el=dijcf
      endif
  100 continue
c
c do NOT square DENSPAR
c
c Square this guy
c
      dens_max_el=dens_max_el*dens_max_el
c
c
      end
c====================================================================
      subroutine setup_densp3(inx,ncf,ncs,densmat,denspar,map_fs)
      implicit real*8 (a-h,o-z)
      dimension inx(12,ncs)
      dimension densmat(ncf,ncf)
      dimension map_fs(ncf)
      dimension denspar(ncs,ncs)          ! output
c--------------------------------------------------------------------
c input :
c inx(12,ncs)       = basis set info
c ncf, ncs          = contracted functions and shells
c densmat           = full density matrix
c map_fs(icf)=ics     mapping from contr.funct to contr.shells
c
c output :
c densp(ics,jcs)    = max absolute elemet of density for each ics,jcs
c--------------------------------------------------------------------
c
      do 10 ics=1,ncs
      icf_b=inx(11,ics)+1
      icf_e=inx(10,ics)
      do 10 icf=icf_b,icf_e
      map_fs(icf)=ics
   10 continue
c
      call zeroit(denspar,ncs*ncs)
c
      do 100 icf=1,ncf
      ics=map_fs(icf)
      do 100 jcf=1,ncf
      jcs=map_fs(jcf)
c
      dijcf=densmat(icf,jcf)
      dijcs=denspar(ics,jcs)
      if(dijcf.gt.dijcs) denspar(ics,jcs)=dijcf
  100 continue
c
c do NOT square DENSPAR
c
      end
c====================================================================
      subroutine make_fock8(ilen,jlen,klen,llen,
     *                       icff,jcff,kcff,lcff,
     *                       dens,fock,lind,dmax,
     *                       intx8,xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(*),fock(*)
      dimension xinteg8(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx8=intx8+1
                     xint=xinteg8(intx8)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin4=4.d0*xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     jj=lind(jff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
                     ik=ii+kff
                     il=ii+lff
c
                     if(jff.ge.kff) then
                       jk=jj+kff
                     else
                       jk=kk+jff
                     end if
                     if(jff.ge.lff) then
                       jl=jj+lff
                     else
                       ll=lind(lff)
                       jl=ll+jff
                     end if
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
                     fock(ik)=fock(ik)-xint*dens(jl)
                     fock(jl)=fock(jl)-xint*dens(ik)
                     fock(jk)=fock(jk)-xint*dens(il)
                     fock(il)=fock(il)-xint*dens(jk)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_fock4(ilen,jlen,klen,llen,
     *                       icff,jcff,kcff,lcff,
     *                       dens,fock,lind,dmax,
     *                       intx4,xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(*),fock(*)
      integer*4 xinteg4(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx4=intx4+1
                     xint=xinteg4(intx4)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin4=4.d0*xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     jj=lind(jff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
                     ik=ii+kff
                     il=ii+lff
c
                     if(jff.ge.kff) then
                       jk=jj+kff
                     else
                       jk=kk+jff
                     end if
                     if(jff.ge.lff) then
                       jl=jj+lff
                     else
                       ll=lind(lff)
                       jl=ll+jff
                     end if
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
                     fock(ik)=fock(ik)-xint*dens(jl)
                     fock(jl)=fock(jl)-xint*dens(ik)
                     fock(jk)=fock(jk)-xint*dens(il)
                     fock(il)=fock(il)-xint*dens(jk)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_fock2(ilen,jlen,klen,llen,
     *                       icff,jcff,kcff,lcff,
     *                       dens,fock,lind,dmax,
     *                       intx2,xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(*),fock(*)
      integer*2 xinteg2(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx2=intx2+1
                     xint=xinteg2(intx2)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
c
                     xin4=4.d0*xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     jj=lind(jff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
                     ik=ii+kff
                     il=ii+lff
c
                     if(jff.ge.kff) then
                       jk=jj+kff
                     else
                       jk=kk+jff
                     end if
                     if(jff.ge.lff) then
                       jl=jj+lff
                     else
                       ll=lind(lff)
                       jl=ll+jff
                     end if
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
                     fock(ik)=fock(ik)-xint*dens(jl)
                     fock(jl)=fock(jl)-xint*dens(ik)
                     fock(jk)=fock(jk)-xint*dens(il)
                     fock(il)=fock(il)-xint*dens(jk)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine dftfock_disk(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dens,       fock,
     *        ncs,         densp,      map_fs)
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
      dimension dens(*),fock(*)
c
      dimension lind(*)
c----------------------------------------------------------------
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
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
c     write(6,*)
c    * 'dens_max=',dens_max,' thers1=',thres1,' do=',doint72,doint62
c----------------------------------------------------------------
c if screenig density densp() is not used then 
c
           dmax1=1.0d0
           dmax=dmax1*t01
c----------------------------------------------------------------
c
           call secund(r1)
           call elapsec(e1)
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
           call secund(r2)
           call elapsec(e2)
           readcpu3=readcpu3+(r2-r1)
           readela3=readela3+(e2-e1)
c----------------------------------------------------------------
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
c     write(6,*)' integ=',integi72,integi62,integin2,integin4,integin8
c     call f_lush(6)
c
           call secund(r1)
           call elapsec(e1)
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
           call secund(r2)
           call elapsec(e2)
           readcpu4=readcpu4+(r2-r1)
           readela4=readela4+(e2-e1)
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
c
c        write(6,*)' bl=',isbl,' size=',isize,' nq=',nqrt
c
         do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)

c        write(6,*)' bl=',isbl,' size=',isize,' nq=',nqrt
c        write(6,*)'         qr=',iqrt,' labs=',icff,jcff,kcff,lcff
c        call f_lush(6)
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
c--------------------------------------------------------
c
c        write(6,*)' case=',icase,' labs=',icff,jcff,kcff,lcff
c        write(6,*)' case=',icase,' lens=',ilen,jlen,klen,llen
c             call f_lush(6)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c shell :
c          ics=map_fs(icff+1)
c          jcs=map_fs(jcff+1)
c          kcs=map_fs(kcff+1)
c          lcs=map_fs(lcff+1)
c
c screening density
c
c          dij=densp(ics,jcs)
c          dkl=densp(kcs,lcs)
c          dmax1=4.0d0*max(dij,dkl)
c          dmax1=1.0d0
c          dmax=dmax1*t01
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         using  dmax*abs(x) > 0.5 is the same as 
c
c                     abs(x) > 0.5/dmax 
c
c         and no multiplications below 
c
c
              dmax2=half/dmax
c
c
              if(integin8.gt.intx8 .and. icase.eq.1) then
                 call make_dftfock8(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          intx8,xinteg8)
              endif
c
              if(integin4.gt.intx4 .and. icase.eq.2) then
                 call make_dftfock4(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          intx4,xinteg4)
              endif
c
              if(integin2.gt.intx2 .and. icase.eq.3) then
                 call make_dftfock2(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          intx2,xinteg2)
              endif
c
              if(doint62) then
               if(integi62.gt.int62 .and. icase.eq.4) then
                 if(dmax1*1.0d-6.gt.thres1) then
                 call make_dftfock2(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          int62,xinte62)
                 endif
               endif
              endif
c
              if(doint72) then
               if(integi72.gt.int72 .and. icase.eq.5) then
                 if(dmax1*1.0d-7.gt.thres1) then
                 call make_dftfock2(ilen,jlen,klen,llen,
     *                          icff,jcff,kcff,lcff,
     *                          dens,fock,lind,dmax2,
     *                          int72,xinte72)
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
         call dscal(ncf*(ncf+1)/2,resc,fock,1)
      endif
c
      end
c=====================================================================
      subroutine make_dftfock8(ilen,jlen,klen,llen,
     *                       icff,jcff,kcff,lcff,
     *                       dens,fock,lind,dmax,
     *                       intx8,xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(*),fock(*)
      dimension xinteg8(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx8=intx8+1
                     xint=xinteg8(intx8)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin4=4.d0*xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_dftfock4(ilen,jlen,klen,llen,
     *                       icff,jcff,kcff,lcff,
     *                       dens,fock,lind,dmax,
     *                       intx4,xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(*),fock(*)
      integer*4 xinteg4(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx4=intx4+1
                     xint=xinteg4(intx4)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin4=4.d0*xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_dftfock2(ilen,jlen,klen,llen,
     *                       icff,jcff,kcff,lcff,
     *                       dens,fock,lind,dmax,
     *                       intx2,xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(*),fock(*)
      integer*2 xinteg2(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx2=intx2+1
                     xint=xinteg2(intx2)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin4=4.d0*xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
c
                     fock(ij)=fock(ij)+xin4*dens(kl)
                     fock(kl)=fock(kl)+xin4*dens(ij)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine read_int(iunit ,xint,isize)
      implicit real*8 (a-h,o-z)
      dimension xint(isize)
      read(iunit) xint
      end
c=====================================================================
      subroutine read_int8(iunit ,xint,isize)
      implicit real*8 (a-h,o-z)
      dimension xint(isize)
      read(iunit) xint
c          write(6,*) ' after read-in integ R*8'
c          do i=1,isize
c             x=xint(i)*1.0d-10
c             write(6,66) i,x
c          enddo
c66   format(i5,3x,f15.8)
      end
c=====================================================================
      subroutine read_int4(iunit ,xint,isize)
      implicit real*8 (a-h,o-z)
      integer*4 xint(isize)
      read(iunit) xint
c          write(6,*) ' after read-in integ I*4'
c          do i=1,isize
c             x=xint(i)*1.0d-10
c             write(6,66) i,x
c          enddo
c66   format(i5,3x,f15.8)
      end
c=====================================================================
      subroutine read_int2(iunit ,xint,isize)
      implicit real*8 (a-h,o-z)
      integer*2 xint(isize)
      read(iunit) xint
c          write(6,*) ' after read-in integ I*4'
c          do i=1,isize
c             x=xint(i)*1.0d-10
c             write(6,66) i,x
c          enddo
c66   format(i5,3x,f15.8)
      end
c==============================================================
      subroutine read_33(iunit,iblstore,iqstore,nsplit) 
      integer iblstore,iqstore,nsplit
      read(iunit) iblstore,iqstore,nsplit
      end
c==============================================================
      subroutine read_44(iunit,ifromto,nofinteg,nsplit)
      integer ifromto(2,nsplit),nofinteg(5,nsplit)
      read (iunit) ifromto,nofinteg
      end
c=====================================================================
      subroutine get_block4(isup,nftcbl4,isupbl)
      dimension nftcbl4(*)
c
      isupbl=nftcbl4(isup)
c
      end
