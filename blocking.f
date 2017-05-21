C===============================================================
      subroutine blockin1(iforwhat,natoms,ncs,inx,bl,datnuc,datbas,nbl1)
c---------------------------------------------------------------
c blocking procedure for single shells (contracted)
c
c Input
c---------
c iforwhat  - shows type of task (integrals to be calculated)
c natoms    - number of atoms
c ncs       - number of contracted shells
c inx()     - basis set info
c bl(*)     - storage for everything
c datnuc()  - nuclear data
c datbas()  - basis set data
c
c
c Output
c-------
c nbl1      - number of blocks of single shells
c---------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*4 text
c
      common /outfile/ ioutput
      common /route/ iroute
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
c
      dimension bl(*)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c---------------------------------------------------------------
c make blocks of single shells :
c
c allocate memory
c
      call getint(ncs+1,nblock1)
      call getint(ncs  ,nblock1_back)
c
c save these addresses :
c
      call setival('nblock1' ,nblock1)
      call setival('nblock1b',nblock1_back)
c
c999  call blk_shells(iforwhat,ncs,inx,iroute,datnuc,datbas,
c999 *                bl(nblock1),nbl1)
c no limits :
      call blk_shells(ncs,inx,iroute,datnuc,datbas,bl(nblock1),nbl1)
c
      call setival('nbl1' ,nbl1)
c
c output : nblock1(i-1)+1  : first shell in i-block1
c output : nblock1(i)      : last  shell in i-block1
c output : nbl1 - number of blocks of shells
c---------------------------------------------------------------
c make nblock1_back(*) array
c
      call make_nblock1_back(ncs,nbl1,bl(nblock1),bl(nblock1_back))
c
c output : nblock1_back(ics) => block of shells ics belongs to
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_shells(ncs,inx,iroute,datnuc,datbas,nblock1,nbl1)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      logical txs93,txs95, cond1,cond2,cond3,cond4,cond5,cond6
      logical cond0
      dimension nblock1(0:ncs)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c-----------------------------------------------------------------
c calculate number of different shells (i.e. different blocks)
c according to two different criterions (txs93 or txs95)
c-----------------------------------------------------------------
c NO  limits on a number of shells per block
c-----------------------------------------------------------------
      iexch=0
      do 55 ics0=1,ncs-1
         ist0=inx(12,ics0)                 ! type of shell
         iat0=inx(2,ics0)                  ! atom
         nzi0=0                            ! charge
         if(iat0.gt.0) nzi0=datnuc(1,iat0)
         isc0=inx(5,ics0)-inx(1,ics0)      ! contraction lenght
         isp0=inx(1,ics0)+1                ! first primitive
         isp0l=inx(5,ics0)                 ! last  primitive
         exi0=datbas(1,isp0)               ! exponent of first primitive
         exi0l=datbas(1,isp0l)             ! exponent of last  primitive
         coe0=datbas(2,isp0)               ! con.coef of first primitive
         ngci0=inx(4,ics0)                 ! general contraction deep
c
c next shell :
c
         ics1=ics0+1
         ist1=inx(12,ics1)
         iat1=inx(2,ics1)
         nzi1=0
         if(iat1.gt.0) nzi1=datnuc(1,iat1)
         isc1=inx(5,ics1)-inx(1,ics1)
         isp1=inx(1,ics1)+1
         isp1l=inx(5,ics1)                 ! last  primitive
         exi1=datbas(1,isp1)
         exi1l=datbas(1,isp1l)             ! exponent of last  primitive
         coe1=datbas(2,isp1)
         ngci1=inx(4,ics1)
c
c---------------------------------------------
c check against criterion (3 or 6 conditions) :
c
           cond1=.false.
           cond2=.false.
           cond3=.false.
           cond4=.false.
           cond5=.false.
           cond6=.false.
           if(ist0 .eq. ist1 ) cond1=.true.
           if(isc0 .eq. isc1 ) cond2=.true.
           if(ngci0.eq.ngci1 ) cond3=.true.
           if(nzi0 .eq. nzi1 ) cond4=.true.
           if(exi0 .eq. exi1 ) cond5=.true.
           if(coe0 .eq. coe1 ) cond6=.true.
c---------------------------------------------
           if(cond1 .and. iroute.eq.1) then
c -- because of possible numerical instability
              call get_cond0(exi0,exi0l,exi1,exi1l,ist0,cond0)
           endif
c---------------------------------------------
c
           txs93=.false.
           txs95=.false.
           if(cond1.and.cond2.and.cond3) txs93=.true.
           if(cond4.and.cond5.and.cond6.and.txs93) txs95=.true.
c
           if(iroute.eq.1) then
c2002         if(.not.txs93) then
              if( (.not.txs93) .or. (.not.cond0) ) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
              endif
           else
              if(.not.txs95) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
              endif
           endif
c
   55 continue
c
      nblock1(0)=0
      nblock1(iexch+1)=ncs
c
c--------------------------------------------------------
      nbl1=iexch+1              ! number of shlell-blocks
c--------------------------------------------------------
      maxshell=0
      do ibl1=1,nbl1
         ibeg=nblock1(ibl1-1)+1
         iend=nblock1(ibl1)
         ishell=iend-ibeg+1
         maxshell=max(maxshell,ishell)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
      enddo
c--------------------------------------------------------
      write(ioutput,495) nbl1,ncs, maxshell
  495 format( ' Number of Blocks of Contracted Shells         =',i7/
     *        ' total number of contracted shells             =',i7/
     *        ' maximum number of shells/block                =',i7/)
ctest
c     do ibl1=1,nbl1
c        ibeg=nblock1(ibl1-1)+1
c        iend=nblock1(ibl1)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
c        write(91,*)' blk1=',ibl1,' type=',itype,' contr.=',icont,
c    *              ' no of shells=',iend-ibeg+1
c     enddo
ctest
c---------------------------------------------
      end
c===============================================================
c NOT used now
      subroutine blk_shells2(ncs,inx,iroute,datnuc,datbas,nblock1,nbl1)
      implicit real*8 (a-h,o-z)
c999  parameter (maxvalue=48)
      parameter (maxvalue=150)
      common /outfile/ ioutput
      logical txs93,txs95, cond1,cond2,cond3,cond4,cond5,cond6
      dimension nblock1(0:ncs)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c-----------------------------------------------------------------
c calculate number of different shells (i.e. different blocks)
c according to two different criterions (txs93 or txs95)
c-----------------------------------------------------------------
c put limits on a number of shells per block depending on a type
c of shells according to : n_limit*size*contr=150
c-----------------------------------------------------------------
      iexch=0
      ishells=0
      do 55 ics0=1,ncs-1
         ist0=inx(12,ics0)                 ! type of shell
         iat0=inx(2,ics0)                  ! atom
         nzi0=0                            ! charge
         if(iat0.gt.0) nzi0=datnuc(1,iat0)
         isc0=inx(5,ics0)-inx(1,ics0)      ! contraction lenght
         isp0=inx(1,ics0)+1                ! first primitive
         exi0=datbas(1,isp0)               ! exponent of first primitive
         coe0=datbas(2,isp0)               ! con.coef of first primitive
         ngci0=inx(4,ics0)                 ! general contraction deep
c
c limits : n_limit*size*contr=48
c
         isize=inx(3,ics0)
         if(isize.eq.5) isize=6
         if(isize.eq.7) isize=10
         limshell=maxvalue/(isize*isc0)
         if(limshell.le.0) limshell=1
c
c next shell :
c
         ics1=ics0+1
         ist1=inx(12,ics1)
         iat1=inx(2,ics1)
         nzi1=0
         if(iat1.gt.0) nzi1=datnuc(1,iat1)
         isc1=inx(5,ics1)-inx(1,ics1)
         isp1=inx(1,ics1)+1
         exi1=datbas(1,isp1)
         coe1=datbas(2,isp1)
         ngci1=inx(4,ics1)
c
c---------------------------------------------
c check against criterion (3 or 6 conditions) :
c
           cond1=.false.
           cond2=.false.
           cond3=.false.
           cond4=.false.
           cond5=.false.
           cond6=.false.
           if(ist0 .eq. ist1 ) cond1=.true.
           if(isc0 .eq. isc1 ) cond2=.true.
           if(ngci0.eq.ngci1 ) cond3=.true.
           if(nzi0 .eq. nzi1 ) cond4=.true.
           if(exi0 .eq. exi1 ) cond5=.true.
           if(coe0 .eq. coe1 ) cond6=.true.
c
           txs93=.false.
           txs95=.false.
           if(cond1.and.cond2.and.cond3) txs93=.true.
           if(cond4.and.cond5.and.cond6.and.txs93) txs95=.true.
c
           ishells=ishells+1
C
           if(iroute.eq.1) then
              if(.not.txs93 .or. ishells.ge.limshell) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
                 ishells=0
              endif
           else
              if(.not.txs95 .or. ishells.ge.limshell) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
ctest..............................................
c     write(91,*)
c    * ' block1=',iexch,' type=',ist0,' con=',isc0,
c    * ' limit=',limshell,' n-shell=',ishells
ctest..............................................
                 ishells=0
              endif
           endif
c
   55 continue
c
      nblock1(0)=0
      nblock1(iexch+1)=ncs
c
c--------------------------------------------------------
      nbl1=iexch+1              ! number of shlell-blocks
c--------------------------------------------------------
      maxshell=0
      do ibl1=1,nbl1
         ibeg=nblock1(ibl1-1)+1
         iend=nblock1(ibl1)
         ishell=iend-ibeg+1
         maxshell=max(maxshell,ishell)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
      enddo
c--------------------------------------------------------
      write(ioutput,495) nbl1,ncs, maxshell
  495 format( ' Number of Blocks of Contracted Shells         =',i7/
     *        ' total number of contracted shells             =',i7/
     *        ' maximum number of shells/block                =',i7/)
ctest
c     do ibl1=1,nbl1
c        ibeg=nblock1(ibl1-1)+1
c        iend=nblock1(ibl1)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
c        write(91,*)' blk1=',ibl1,' type=',itype,' contr.=',icont,
c    *              ' no of shells=',iend-ibeg+1
c     enddo
ctest
c---------------------------------------------
      end
c===============================================================
c NOT used now
      subroutine blk_shells1
     *       (iforwhat,ncs,inx,iroute,datnuc,datbas,nblock1,nbl1)
      implicit real*8 (a-h,o-z)
c     parameter (maxvalue=36)
      parameter (maxvalue=48)
      common /outfile/ ioutput
      logical txs93,txs95, cond1,cond2,cond3,cond4,cond5,cond6
      dimension nblock1(0:ncs)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c-----------------------------------------------------------------
c calculate number of different shells (i.e. different blocks)
c according to two different criterions (txs93 or txs95)
c-----------------------------------------------------------------
c put limits on a number of shells per block depending on a type
c of shells according to : n_limit*size*contr=36
c-----------------------------------------------------------------
c factors for limit:
c
      iexch=0
      ishells=0
      do 55 ics0=1,ncs-1
         ist0=inx(12,ics0)                 ! type of shell
         iat0=inx(2,ics0)                  ! atom
         nzi0=0                            ! charge
         if(iat0.gt.0) nzi0=datnuc(1,iat0)
         isc0=inx(5,ics0)-inx(1,ics0)      ! contraction lenght
         isp0=inx(1,ics0)+1                ! first primitive
         exi0=datbas(1,isp0)               ! exponent of first primitive
         coe0=datbas(2,isp0)               ! con.coef of first primitive
         ngci0=inx(4,ics0)                 ! general contraction deep
c
c limits : n_limit*size*contr=36
c
         isize=inx(3,ics0)
         if(isize.eq.5) isize=6
         if(isize.eq.7) isize=10
         limshell=maxvalue/(isize*isc0)
         if(limshell.le.0) limshell=1
         if(iforwhat.eq.1) then
c           2-el. ordinary integrals (scf,cphf)
            limshell=10*limshell
         endif
         if(iforwhat.eq.2) then
c           2-el. giao integrals (nmr)
            limshell= 5*limshell
         endif
         if(iforwhat.eq.3) then
c           2-el. integral gradient derivatives (forecs)
            limshell= 3*limshell
         endif
c.....   if(iforwhat.eq.4) then
c           2-el. integral hessian derivatives
c           limshell= 3*limshell
c.....   endif
         if(iforwhat.eq.5) then
c           2-el. ordinary integrals for lmp2
            limshell=min(limshell,12)
         endif
c
c next shell :
c
         ics1=ics0+1
         ist1=inx(12,ics1)
         iat1=inx(2,ics1)
         nzi1=0
         if(iat1.gt.0) nzi1=datnuc(1,iat1)
         isc1=inx(5,ics1)-inx(1,ics1)
         isp1=inx(1,ics1)+1
         exi1=datbas(1,isp1)
         coe1=datbas(2,isp1)
         ngci1=inx(4,ics1)
c
c---------------------------------------------
c check against criterion (3 or 6 conditions) :
c
           cond1=.false.
           cond2=.false.
           cond3=.false.
           cond4=.false.
           cond5=.false.
           cond6=.false.
           if(ist0 .eq. ist1 ) cond1=.true.
           if(isc0 .eq. isc1 ) cond2=.true.
           if(ngci0.eq.ngci1 ) cond3=.true.
           if(nzi0 .eq. nzi1 ) cond4=.true.
           if(exi0 .eq. exi1 ) cond5=.true.
           if(coe0 .eq. coe1 ) cond6=.true.
c
           txs93=.false.
           txs95=.false.
           if(cond1.and.cond2.and.cond3) txs93=.true.
           if(cond4.and.cond5.and.cond6.and.txs93) txs95=.true.
c
           ishells=ishells+1
C
           if(iroute.eq.1) then
              if(.not.txs93 .or. ishells.ge.limshell) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
                 ishells=0
              endif
           else
              if(.not.txs95 .or. ishells.ge.limshell) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
ctest..............................................
c     write(91,*)
c    * ' block1=',iexch,' type=',ist0,' con=',isc0,
c    * ' limit=',limshell,' n-shell=',ishells
ctest..............................................
                ishells=0
              endif
           endif
c
   55 continue
c
      nblock1(0)=0
      nblock1(iexch+1)=ncs
c
c--------------------------------------------------------
      nbl1=iexch+1              ! number of shlell-blocks
c--------------------------------------------------------
      maxshell=0
      do ibl1=1,nbl1
         ibeg=nblock1(ibl1-1)+1
         iend=nblock1(ibl1)
         ishell=iend-ibeg+1
         maxshell=max(maxshell,ishell)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
      enddo
c--------------------------------------------------------
      write(ioutput,495) nbl1,ncs, maxshell
  495 format( ' Number of Blocks of Contracted Shells         =',i7/
     *        ' total number of contracted shells             =',i7/
     *        ' maximum number of shells/block                =',i7/)
ctest
c     do ibl1=1,nbl1
c        ibeg=nblock1(ibl1-1)+1
c        iend=nblock1(ibl1)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
c        write(91,*)' blk1=',ibl1,' type=',itype,' contr.=',icont,
c    *              ' no of shells=',iend-ibeg+1
c     enddo
ctest
c---------------------------------------------
      end
c===============================================================
      subroutine make_nblock1_back(ncs,nbl1,nblock1,nblock1_back)
      dimension nblock1(0:ncs),nblock1_back(ncs)
c construct nblock1_back() array
      do ibl=1,nbl1
         icsb=nblock1(ibl-1)+1
         icse=nblock1(ibl)
         do ics=icsb,icse
            nblock1_back(ics)=ibl
         enddo
      enddo
      end
c===============================================================
      subroutine blockin2(lcore,nprint,ncs,inx,bl,nbl1,
c output :
     *                    nbl2,minpr,maxpr,  xminpr_bl, xmaxpr_bl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c
      common /memmax/ ispblx, maxme1,iforwhat
      common /route/ iroute
c all below are output adresses :
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1c/ map_ij_bl2
      common /memors/ nsym,ijshp,isymm
c up to here
c
c     common /intbl/ifp,inxx(100)
c
      dimension bl(*)
      dimension inx(12,*)
c---------------------------------------------------------------
c blocking procedure for pairs
c---------------------------------------------------------------
c allocate memory for iis(ijcs)--> ics  & jjs(ijcs)-->jcs arrays
c allocate memory for map_ij_bl2(ijcs)-->ijbl2
c
      ncsp=ncs*(ncs+1)/2
c
      call getint(ncsp,iisd)
      call getint(ncsp,jjsd)
      call getint(ncsp,map_ij_bl2)
c
c Store these addresses in the  common /memor1/ iisd,jjsd,ijbld
c                               common /memor1c/ map_ij_bl2
c---------------------------------------------------------------
c Constructe blocks of contracted shell pairs , setup the arrays:
c npar, ijbl and nsupb
c
        call blk_pairs(bl,ncs,inx,bl(iisd),bl(jjsd),nbl1,
     *                nbl2,nsupb,npard,ijbld, mxsize,map_ij_bl2,
     *                lcore,nsym)
c
c  output : arrays' addresses : nsupb, npar, ijbl and  mxsize )
c---------------------------------------------------------------
      if(nsym.gt.0) then
         ncsp=ncs*(ncs+1)/2
         call getint(nsym*ncsp, ijshp)
         call getival('SymShPr',ifp1)
         call parsym(ncs,bl(ifp1),bl(iisd),bl(jjsd),bl(ijshp),nsym)
      endif
c---------------------------------------------------------------
c About prices :
c I need two prices : per 1 integral & per one WHOLE big-block
c
      nbl4=nbl2*(nbl2+1)/2
      call getmem(2*nbl4, ncost)
      call price(inx,bl(iisd),bl(jjsd),bl(ijbld),bl(npard),
     *           bl(nsupb),nbl2,nprint,
     *           bl(ncost),minpr,maxpr,
     *           bl(ncost+nbl4+1),xminpr_bl,xmaxpr_bl)
c
c bl(ncost) - prices per integral
c bl(ncost+nbl4+1) - prices for the whole block
c
c minpr & maxpr - minimum & maximum price for one integral
c xminpr_bl & xmaxpr_bl - minimum & maximum price for a block
c---------------------------------------------------------------
c This is the end of blocking procedure for PAIRS
c
c The following is known at this point :
c
c  pairs ijcs are given by ics=iis(ijcs) and jcs=jjs(ijcs)
c  number of pair-blocks = NBL2
c  number of pairs in the block ibl is NPAR(ibl)
c  which pairs belong to this block : ijcs=ijbl(ibl,1-npar)
c---
c The number of BLOCKS of contracted shell Quartets would be
c NBL4=NBL2*(NBL2+1)/2 if they are not too big. I call
c them Super-blocks of c.s.quartets.The maximum size for
c each of them is known and kept in array MAXSIZE(i=1,NSUPB)
c (bl(mxsize)) (no of quartets). According to this max.size
c super-blocks will be split into a number of smaller blocks
c which is also known (NBLOKS) but is NOT transfered.
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairs(bl,ncs,inx,iis,jjs,nbl1,
     *                     nbl2,nsupb,npard,ijbld,mxsize,map_ij_bl2,
     *                     lcore,nsym)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /memmax/ ispblx, maxme1,iforwhat
      common /route/ iroute
      dimension bl(*)
      dimension inx(12,*)
      dimension iis(*),jjs(*)
c---------------------------------------------------------------
c Construct pairs of contracted shells ( I>=J )
c and setup iis(ijcs) & jjs(ijcs) arrays :
c
      call pairs(ncs,inx,iis,jjs)  ! eliminate it latter on
c---------------------------------------------------------------
c Calculate number of pair-blocks NBL2 & maximum pairs/block MAXPAR
c
      call getival('nblock1' ,nblock1)
      call blk_pairs_calc(nbl1,bl(nblock1),nbl2,maxpar)
c
c output : nbl2
c output : maxpar    - maximum number of pairs per block2
c---------------------------------------------------------------
c number of quartet-blocks (super-blocks) :
c
      nbl4=nbl2*(nbl2+1)/2
c
c---------------------------------------------------------------
c Allocate memory for NPAR and IJBL and for MXSIZE
c
      call getint(nbl2       ,npard)
      call getint(nbl2*maxpar,ijbld)
      call getint(nbl4       ,mxsize)
c
c output : addresses : npard,ijbld & mxsize
c---------------------------------------------------------------
c Constructe blocks of contracted pairs:Setup IJBL(nbl2,*)and  NPAR(*)
c
      call blk_pairs_make(nbl1,bl(nblock1),nbl2,bl(npard),bl(ijbld),
     *                    bl(map_ij_bl2) )
c---------------------------------------------------------------
c  Set up the MXSIZE and NSUPB vectors
c mxsize(isupbl)=> maiximum size of a small block in isupbl big-block
c nsupb (isupbl)=> number of small blocks in isupbl big-block
c
      call getint(nbl4,nsupb)
c
      call blkmemor(bl,nbl2,bl(ijbld),bl(npard),
     *              iis,jjs,inx, lcore,iroute)
c
      call blksizer(bl,bl(nsupb),nbl2,bl(ijbld),bl(npard),
     *              iis,jjs,inx, bl(mxsize),lcore,nsym,ncs,iroute)
c
c---------------------------------------------------------------
      end
c===============================================================
      subroutine pairs(ncs,inx,iis,jjs)
      dimension inx(12,*)
      dimension iis(*),jjs(*)
c
c  constructe pairs of contracted shells with reordered basis set: I>=J
c
      ijcs=0
      do 450 ics=1,ncs
      do 450 jcs=1,ics
        ijcs=ijcs+1
        iis(ijcs)=ics
        jjs(ijcs)=jcs
  450 continue
c
      end
c===============================================================
      subroutine blk_pairs_calc(nbl1,nblock1,nbl2,maxpar)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
c
c input : nbl1 - number of single shell blocks
c
c output: nbl2 - number of shell-pairs  blocks
c output: maxpar - max. number of pairs/block
c
c
c
      ncsp=0
      maxpar=0
      ijbl=0             !    blocks counter (with limit)
      do 100 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 200 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar=0
            do 300 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 400 jcs=jbeg,jenx
                  ncsp=ncsp+1
                  if(ijpar.ge.limpair) then
                     ijbl=ijbl+1
                     ijpar=0
                  endif
                  ijpar=ijpar+1
ccccccccccccccc   npar(ijbl)=ijpar     ccccccccccccccccccccccc
                  if(ijpar.gt.maxpar) maxpar=ijpar
  400          continue
  300       continue
  200    continue
  100 continue
c
      nbl2=ijbl   ! number of pair-blocks with the limit limpair
c
c---------------------------------------------------------------
      call setival('ncspairs',ncsp)
c---------------------------------------------------------------
      write(ioutput,495) nbl2,ncsp,maxpar
  495 format( ' Number of Blocks of Contracted Shell Pairs    =',i7/
     *        ' total number of contracted shell pairs        =',i7/
     *        ' maximum number of shell-pairs/block           =',i7)
c
c     nbl4=nbl2*(nbl2+1)/2
c     ntotq=ncsp*(ncsp+1)/2
c     write(ioutput,496) nbl4,ntotq
c 496 format(/' Minimum number of Blocks of Contracted Quartets =',i6/
c    *        '    total number of contracted quartets =',i15/)
c---------------------------------------------------------------
c
      end
c===============================================================
      subroutine blk_pairs_make(nbl1,nblock1,nbl2,npar,ijbl,map_ij_bl2)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      common /intlim/ limxmem,limblks,limpair
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncsp)
c
c nbl1 - number of single shell blocks
c nbl2 - number of shell-pairs  blocks
c
c
c  constructe pairs of contracted shells
c
      ijblock=0             !    blocks counter
      do 100 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 200 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar=0
            do 300 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 400 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  if(ijpar.ge.limpair) then
                     ijblock=ijblock+1
                     ijpar=0
                  endif
                  ijpar=ijpar+1
                  npar(ijblock)=ijpar
                  ijbl(ijblock,ijpar)=ijcs
                  map_ij_bl2(ijcs)=ijblock
  400          continue
  300       continue
  200    continue
  100 continue
c
      last_big2=ijblock
      call setival('last_big2',ijblock)
      call setival('last_mid2',ijblock)
c---------------------------------------------------------------
c     write(6,*)' nbl2=',nbl2,' ijblock=',ijblock
c
c     ijpar_sum=0
c     do 2500 ibl=1,nbl2
c     ijpar=npar(ibl)
c     write(ioutput,502) ibl,ijpar
c     ijpar_sum=ijpar_sum+ijpar
c2500 continue
c 502 format('Pair-Block =',i4,' : no of shell-pairs =',i4)
c     write(6,*)'total number of pairs=',ijpar_sum
c---------------------------------------------------------------
c
      end
c=======================================================================
      subroutine blksizer(bl,nsupb ,nbl2,ijbl,npar,iis,jjs,inx,
     *                    mxsize,lcore,nsym,ncs,iroute)
c-----------------------------------------------------------
c This routine determines the block-size for two-electron
c integrals calculations.
c-----------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical stable
      common /outfile/ ioutput
      common /intlim/ limxmem,limblks,limpair
      common /cpu/ intsize,iacc,icache,memreal
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
c
      common /logic1/ ndege(1)
      common /logic2/ lenn(1)
      common /logic3/ lensm(1)
c
      dimension bl(*)
      dimension mxsize(*), nsupb(*)
      dimension npar(*),ijbl(nbl2,*)
      dimension inx(12,*),iis(*),jjs(*)
c-----------------------------------------------------------
c  iforwhat shows where the Blockin2 routins are called from :
c    from twoint  with iforwhat=1  (ordinary integrals)
c    from Shift2  with iforwhat=2  (GIAO derivatives  )
c    from Force2  with iforwhat=3  (I-st derivatives  )
c    from ??????  with iforwhat=4  (IIed derivatives  )
c    from ptwoint with iforwhat=5  (mp2 integrals     )
c
c  This is used to send info about memory to the Calcint2
c-----------------------------------------------------------
      call getmem(0,last1)
      call retmem(1)
c-----------------------------------------------------------
c After blksizer returns to blkpair there are two more calls
c of memo1_int ( integer allocation using getmem)
c
      last2=nsym*ncs*(ncs+1)/2  +  nbl2*(nbl2+1)/2
      last2=last2/intsize+1
c
      last12=last1+last2
c
      memaval=lcore-last12
c-----------------------------------------------------------
c  Calculates the total number of blocks of contracted
c  shells quartets NBLOKS (only for info).
c  Sets up the MXSIZE(*) and NSUPB(*) arrays for each
c  Super-block (number of them is NBL2*(NBL2+1)/2 ) .
c-----------------------------------------------------------
      call getival('ibas',ibas)
c-----------------------------------------------------------
      maxibuf=0    ! find the maximum size of a integral buffer
      maxindx=0    ! find the maximum size of labels array
c
      nblsmax=0
      maxme1=0
      nbloks=0     ! counter of all small blocks
      nbl12=0      ! counter of super-blocks
c
      do 100 ibl=1,nbl2
      ijpar=npar(ibl)
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
c
      if(iroute.eq.1) then
         call check_exp(ibl,inx,ijbl,nbl2,ijpar,iis,jjs,bl(ibas),
     *                  boamax,rapbmax)
      endif
c
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
      call get_type1(itype,jtype, itype1,jtype1)
c
      nfij=lenn(itype1)*lenn(jtype1)
c
      lci=inx(5,ics1)-inx(1,ics1)
      lcj=inx(5,jcs1)-inx(1,jcs1)
c
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
      ngcij=(ngci+1)*(ngcj+1)
c
      do 100 kbl=1,ibl
      nbl12=nbl12+1
c
      klpar=npar(kbl)
      klcs1=ijbl(kbl,1)
      kcs1=iis(klcs1)
      lcs1=jjs(klcs1)
c
      if(iroute.eq.1) then
         if(kbl.eq.ibl) then
            docmax=boamax
            rcpdmax=rapbmax
         else
            call check_exp(kbl,inx,ijbl,nbl2,klpar,iis,jjs,bl(ibas),
     *                     docmax,rcpdmax)
         endif
      endif
c
      ktype=inx(12,kcs1)
      ltype=inx(12,lcs1)
      call get_type1(ktype,ltype, ktype1,ltype1)
c
      nfkl=lenn(ktype1)*lenn(ltype1)
      nfijkl=nfij*nfkl
c
      lck=inx(5,kcs1)-inx(1,kcs1)
      lcl=inx(5,lcs1)-inx(1,lcs1)
c
      ngck=inx(4,kcs1)
      ngcl=inx(4,lcs1)
      ngckl=(ngck+1)*(ngcl+1)
c-----------------------------------------------------------
c calculate maximum memory needed for one block :
c
c 1) for ordinary two-el.integrals (ifor=1)
c 2) for GIAO two-el. derivatives  (ifor=2)
c 3) for gradient derivatives      (ifor=3)
c 4) for second derivatives        (ifor=4)
c 5) for mp2/lmp2 integrals        (ifor=5)
c
c 6) for non-FTC  integrals        (ifor=11)
c 7) for non-FTC integrals derivatives (ifor=13)
c
      ifor=iforwhat
      if(ifor.eq.5) ifor=1
      if(ifor.eq.11) ifor=1
      if(ifor.eq.13) ifor=3
c
      IF( iroute.eq.1 ) THEN
         call blksize1(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *                 lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *                 memor2,memor2ij,memor2kl,memor4,ifor)
      ELSE
         call blksize2(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *                 lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *                 memor2,memor2ij,memor2kl,memor4,ifor)
      ENDIF
c
c output:
c           memor2,memor2ij,memor2kl,memor4
c-----------------------------------------------------------
c2000 New stuff : limits according to integral's type
c
      call get_max_am(itype1,jtype1,ktype1,ltype1,mmax,nsij,nskl,nqmax)
      call get_limit(mmax,icache,nfijkl,ifor,ibl,kbl,
     *               itype1,jtype1,ktype1,ltype1,
     *               nqrt_limit)
c
      maxsize=nqrt_limit
c---------------------------------------------------
      call getival('stab',istab)
      if(iroute.eq.1 .and. istab.lt.0) then
         call num_stab(mmax,nsij,nskl,nqmax,
     *                 boamax,rcpdmax,docmax,rapbmax,stable)
c        if(.not.stable) then
c        write(6,*)'non-stable s-blok=',nbl12,' pair-blocks=',ibl,kbl
c        endif
         if(.not.stable) maxsize=1
      endif
c---------------------------------------------------
      jump1234=0
 1234 continue
      ikbl=0       ! counter of blocks belonging to one super-block
      maxmem=0
      maxqrt=0
      mxsize(nbl12)=maxsize
      nquart=ijpar*klpar
      if(ibl.eq.kbl) nquart=ijpar*(ijpar+1)/2
c--------------------------------------------------------------------
      call count_small(intsize,ibl,kbl,maxsize,nquart,ijpar,klpar,
     *                 memor2,memor2ij,memor2kl,memor4,
     *                 maxmem,maxqrt,memin4,ikbl)    ! output
c--------------------------------------------------------------------
c put the upper limit on the maxmem :
c
      if(maxmem.ge.limxmem) then
         jump1234=jump1234+1
c
         if(jump1234.eq.1) then
c           estimate the maximum size that fits into limxmem from :
c           limxmem=mem_con+npij*memor2ij+npkl*memor2kl+nq*memor4
c           assuming npij=npij, npkl=1, nq=npij
c
            mem_con=memor2+memin4
            xmaxsize=dble(limxmem-mem_con-memor2kl)
            xmaxsize=xmaxsize/dble(memor2ij+memor4)
            if(xmaxsize .lt. 1.d0) xmaxsize=1.d0
            maxsize=int(xmaxsize)
         else
            if(maxsize.eq.1) then
               write(6,*)' Memory problem: limxmem too small ',limxmem
               call nerror
     *              (970,'blksizer',' Memory problem ',maxmem,limxmem)
            else
               maxsize=2*maxsize/3
               if(maxsize.eq.0) maxsize=1
            endif
         endif
         go to 1234
      endif
c--------------------------------------------------------------------
c select a super-block with maximum memory requirement for ONE block :
c
      nsupb(nbl12)=ikbl
      nbloks=nbloks+ikbl
      if(maxmem.gt.maxme1) then
         maxme1=maxmem
         ispblx=nbl12
      endif
c--------------------------------------------------------------------
c Select a super-block with maximum block-size for ONE block :
c
      if(maxqrt.gt.nblsmax) then
         nblsmax=maxqrt
         ispblq=nbl12
      endif
c--------------------------------------------------------------------
c Select a block (small) with maximum size for an integral buffer :
c        find maximum size of the integral buffer
c        and maximum size of labels array (4,maxindx)
c
      indxsiz1=maxqrt*ngcij*ngckl
      ibufsize=indxsiz1*nfijkl
      indxsize=indxsiz1*4 +maxqrt +4
      if(ibufsize.gt.maxibuf) maxibuf=ibufsize
      if(indxsize.gt.maxindx) maxindx=indxsize
c--------------------------------------------------------------------
  100 continue
c--------------------------------------------------------------------
c Print out info concerning memory requiments :
c--------------------------------------------------------------------
c get number of contracted shell pairs calculated in blk_pairs_cal..
c
c     ncspairs=ncs*(ncs+1)/2  ! if no schwarz screening on pairs
c
c total number of pairs & quartets if no preliminary Schwarz screening :
      ncspatot=ncs*(ncs+1)/2
c number of pairs & quartets after preliminary Schwarz screening :
      call getival('ncspairs',ncspairs)
c
      xquartota=0.5d0*dble(ncspatot)*dble(ncspatot+1)
      xquartets=0.5d0*dble(ncspairs)*dble(ncspairs+1)
      write(ioutput,*)
      write(ioutput,496) nbl12,nbloks,nbloks/nbl12,xquartota,xquartets,
     *                   int(xquartets/nbloks), nblsmax
  496 format( ' Number of Super-Blocks of Contracted Quartets =',i20/
     *        ' number of Small-Blocks of Contracted Quartets =',i20/
     *        '           small-blocks/super-block on average =',i20/
     *        ' Total number of Contracted Shells    Quartets =',f20.0/
     *        ' number of Contracted Shells Quartets retained =',f20.0/
     *        '    average small-block size (in quartets)     =',i20/
     *        '    maximum small-block size (in quartets)     =',i20)
      write(ioutput,900)
  900 format('=================================================')
c
c     write(ioutput,497) maxibuf, maxindx, maxme1
  497 format(/' Maximum lenght of the integral buffer         =',i20/
     *        ' and corresponding array for labels (*4)       =',i20/
     *        ' maximum memory/block (i.e.minimum in two-el.) =',i20/)
c--------------------------------------------------------------------
      if(maxme1.ge.memaval) then
         ioffset=igetival('ioffset')
         irequ=maxme1+last12-ioffset
         iavail=lcore-ioffset
c
c        write(6,*)' mamme1=',maxme1,' memaval=',memaval
c        write(6,*)' lcore=',lcore,' ioffset=',ioffset
c
         call nerror(925,'blksizer',' Memory problem ',irequ,iavail)
      endif
c--------------------------------------------------------------------
      end
c=======================================================================
      subroutine get_type1(itype,jtype, itype1,jtype1)
c
      itype1=itype
      jtype1=jtype
      if(itype.gt.4) itype1=itype-1
      if(jtype.gt.4) jtype1=jtype-1
      if(itype1.gt.5) itype1=itype1-1
      if(jtype1.gt.5) jtype1=jtype1-1
      if(itype.eq.11) itype1=6
      if(itype.eq.12) itype1=7
      if(itype.eq.13) itype1=8
c
      if(jtype.eq.11) jtype1=6
      if(jtype.eq.12) jtype1=7
      if(jtype.eq.13) jtype1=8
c
      end
c=======================================================================
      subroutine count_small(intsize,ibl,kbl,maxsize,nquart,ijpar,klpar,
     *                       memor2,memor2ij,memor2kl,memor4,
     *                       maxmem,maxqrt,memin4,ikbl)
c
c  calculate the number of small blocks of contracted quartets
c
c     write(6,*)' count_small: ibl,kbl=',ibl,kbl
c     write(6,*)' maxsize=',maxsize,' pairs=',ijpar,klpar,
c    *' nquart=',nquart
c
c     call f_lush(6)
c
      if(ibl.gt.kbl) then
        if(nquart.le.maxsize) then
          ikbl=ikbl+1
          memory=ijpar*memor2ij+klpar*memor2kl+ nquart*memor4
          if(memory.gt.maxmem) maxmem=memory
          if(nquart.gt.maxqrt) maxqrt=nquart
          ijsize=ijpar
          klsize=klpar
        else
          ijsize=ijpar
          klsize=klpar
          ijdev=1
          kldev=1
          ijrem=0
          klrem=0
c case 1
          if( ijpar.gt.maxsize.and.klpar.gt.maxsize) then
             ijsize=maxsize
             ijdev=ijpar/ijsize
             ijrem=mod(ijpar,ijsize)
c
             klsize=1
             kldev=klpar/klsize
             klrem=mod(klpar,klsize)
          endif
c case 2
          if( ijpar.gt.maxsize.and.klpar.le.maxsize) then
             ijsize=maxsize
             ijdev=ijpar/ijsize
             ijrem=mod(ijpar,ijsize)
c
             klsize=1
             kldev=klpar/klsize
             klrem=mod(klpar,klsize)
          endif
c case 3
          if( ijpar.le.maxsize.and.klpar.gt.maxsize) then
             klsize=maxsize/ijpar
             kldev=klpar/klsize
             klrem=mod(klpar,klsize)
          endif
c case 4
          if( ijpar.le.maxsize.and.klpar.le.maxsize) then
CKWOL        if(ijpar.lt.klpar) then
CKWOL          ijsize=maxsize/klpar
CKWOL          ijdev=ijpar/ijsize
CKWOL          ijrem=mod(ijpar,ijsize)
CKWOL        else
               klsize=maxsize/ijpar
               kldev=klpar/klsize
               klrem=mod(klpar,klsize)
CKWOL        endif
          endif
c case 5 (no-blocking: 1 quartet per block )
          if(maxsize.eq.1) then
            ijsize=1
            klsize=1
            ijdev=ijpar
            kldev=klpar
            ijrem=0
            klrem=0
          endif
c
          ikbl=ikbl+ijdev*kldev             ! blocks of ijsize x klsize
          if(ijrem.gt.0) ikbl=ikbl+kldev    ! blocks ijrem x klsize
          if(klrem.gt.0) ikbl=ikbl+ijdev    ! blocks ijsize x klrem
          if(ijrem.gt.0 .and. klrem.gt.0) ikbl=ikbl+1
cmemory
          nprij1=ijsize
          nprkl1=klsize
          nprij2=ijrem
          nprkl2=klsize
          nprij3=ijsize
          nprkl3=klrem
          nprij4=ijrem
          nprkl4=klrem
          nqrt1=ijsize*klsize
          nqrt2=ijrem*klsize
          nqrt3=ijsize*klrem
          nqrt4=ijrem*klrem
c--
          memreq1=nprij1*memor2ij+nprkl1*memor2kl+nqrt1*memor4
          memreq2=nprij2*memor2ij+nprkl2*memor2kl+nqrt2*memor4
          memreq3=nprij3*memor2ij+nprkl3*memor2kl+nqrt3*memor4
          memreq4=nprij4*memor2ij+nprkl4*memor2kl+nqrt4*memor4
c--
          memory=max(memreq1,memreq2,memreq3,memreq4)
          if(memory.gt.maxmem) maxmem=memory
          nqrt1234=max(nqrt1,nqrt2,nqrt3,nqrt4)
          if(nqrt1234.gt.maxqrt) maxqrt=nqrt1234
cmemory
        endif
c
      else
c         diagonal case
c
          if(nquart.le.maxsize) then
              ikbl=ikbl+1
              memory=ijpar*memor2ij+klpar*memor2kl+ nquart*memor4
              if(memory.gt.maxmem) maxmem=memory
              if(nquart.gt.maxqrt) maxqrt=nquart
              ijsize=ijpar
              klsize=klpar
          else
              ijsize=sqrt( dble(maxsize) )
              ndev=ijpar/ijsize
              nrem=mod(ijpar,ijsize)
              if(ndev.ge.1) then
                 nqrtd=ijsize*(ijsize+1)/2
                 memory=((3*ijsize+1)*memor2ij)/2 +nqrtd*memor4
                 if(memory.gt.maxmem) then
                    maxmem=memory
                 endif
                 if(nqrtd .gt.maxqrt) maxqrt=nqrtd
                 do 300 n=1,ndev
                 ikbl=ikbl+1
  300            continue
              endif
              if(ndev.ge.2) then
                 nqrtn=ijsize*ijsize
                 memory=2*ijsize*memor2ij +nqrtn*memor4
                 if(memory.gt.maxmem) then
                    maxmem=memory
                 endif
                 if(nqrtn .gt.maxqrt) maxqrt=nqrtn
                 do 400 n1=2,ndev
                 do 400 n2=1,n1-1
                    ikbl=ikbl+1
  400            continue
              endif
              if(nrem.gt.0) then
                 nqrtr=nrem*ijsize
                 ikbl=ikbl+1
                 memory=(nrem+ijsize)*memor2ij+nqrtr*memor4
                 if(memory.gt.maxmem) maxmem=memory
                 if(nqrtr .gt.maxqrt) maxqrt=nqrtr
                 nqrtrem=nqrtr
                 do 500 n=2,ndev
                 nn=(n-1)*ijsize
                 if(nqrtr+nqrtrem.le.maxsize ) then
                     nqrtrem=nqrtrem+nqrtr
                     memory=(nrem+n*ijsize)*memor2ij+nqrtrem*memor4
                     if(memory.gt.maxmem) maxmem=memory
                     if(nqrtrem.gt.maxqrt) maxqrt=nqrtrem
                 else
                     ikbl=ikbl+1
                     nqrtrem=nqrtr
                 endif
  500            continue
                 ikbl=ikbl+1
                 nqrtl=nrem*(nrem+1)/2
                 memory=((3*nrem+1)*memor2ij)/2+ nqrtl*memor4
                 if(memory.gt.maxmem) maxmem=memory
                 if(nqrtl.gt.maxqrt) maxqrt=nqrtl
              endif
          endif
      endif
c--------------------------------------------------------------------
c MAXMEM is the memory needed for pairs and quartets in a small block.
c In addition to it there is memory allocation in prec2ij & prec2kl
c for a whole super-block this small belongs to (memor2) and latter
c on there will be another memory allocation in blockin4 (memin4)
c
      if(intsize.ne.1) then
        memin4=7*(ikbl/intsize+1) +  maxqrt+2*( maxqrt/intsize +1 )
      else
        memin4=7*ikbl+3*maxqrt
      endif
c
c Thus, the total memory needed to handle integral calculations is :
c
      maxmem=maxmem + memor2 + memin4
c
c--------------------------------------------------------------------
      end
c===============================================================
      subroutine blksize1(ibl,kbl,ijpar,klpar,ityp,jtyp,ktyp,ltyp,
     *                    lcix,lcjx,lckx,lclx, ngcix,ngcjx,ngckx,ngclx,
     *                    memory2,memor2ij,memor2kl,memory4,ifor)
      implicit real*8 (a-h,o-z)
      character*4 where
c----------------------------------------------------------------------
c  Calculates  memory needed in a Super-block for :
c
c   1) prec2ij and prec2kl (called for each super-block)
c                                     - memory2
c
c  for each block in a super-block calculates :
c
c   2) memory needed for ONE pair ij  - memor2ij
c   3) memory needed for ONE pair kl  - memor2kl
c   4) memory needed for ONE quartet  - memory4
c
c----------------------------------------------------------------------
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c----------------------------------------------------------------------
      common /logic4/ nfu(1)
      common /logic1/ ndege(1)
      common /logic2/ len(1)
      common /logic3/ lensm(1)
c----------------------------------------------------------------------
c ifor=1   for ordinary two-el.integrals
c ifor=2   for NMR GIAO integral derivatives
c ifor=3   for gradient integral derivatives (first  derivatives)
c ifor=4   for hessian  integral derivatives (second derivatives)
c----------------------------------------------------------------------
c
      n_times=1
      if(ifor.eq.2) n_times= 6 ! giao Ist derivatives
      if(ifor.eq.3) n_times= 9 ! gradient derivatives
      if(ifor.eq.4) n_times=45 ! hessian  derivatives
c
      where='    '
      if(ifor.eq.2) where='shif'
      if(ifor.eq.3) where='forc'
      if(ifor.eq.4) where='hess'
c----------------------------------------------------------------------
      lci=lcix
      lcj=lcjx
      lck=lckx
      lcl=lclx
      ngci=ngcix
      ngcj=ngcjx
      ngck=ngckx
      ngcl=ngclx
c----------------------------------------------------------------------
c set up commons /obarai/ and /shell/
c
      NQI=NDEGE(ITYP)
      NQJ=NDEGE(JTYP)
      NQK=NDEGE(KTYP)
      NQL=NDEGE(LTYP)
      NSIJ=NQI+NQJ-1
      NSKL=NQK+NQL-1
c
      if(where.eq.'shif'.or. where.eq.'forc') then
        nsij=nsij+1
        nskl=nskl+1
      endif
      if(where.eq.'hess') then
        nsij=nsij+2
        nskl=nskl+2
      endif
C
      MMAX=NSIJ+NSKL-1
      mmax1=mmax-1
C
      LNI=LEN(ITYP)
      LNJ=LEN(JTYP)
      LNK=LEN(KTYP)
      LNL=LEN(LTYP)
      LNIJKL=LNI*LNJ*LNK*LNL
C
      LNIJ=LENSM(NSIJ)
      LNKL=LENSM(NSKL)
C
       NSIJ1=NSIJ+1
       NSKL1=NSKL+1
       NQIJ=NQI
       IF(NQJ.GT.NQI) NQIJ=NQJ
       NQIJ1=NQIJ+1
       NQKL=NQK
       IF(NQL.GT.NQK) NQKL=NQL
       NQKL1=NQKL+1
c
       LSHELLT=0
       IF(ITYP.EQ.3) LSHELLT=LSHELLT+1
       IF(JTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(KTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(LTYP.EQ.3) LSHELLT=LSHELLT+1
c----------------------------------------------------------------------
       lcij=lci*lcj
       lckl=lck*lcl
c----------------------------------------------------------------------
c for new general contraction handling :
c
      ngcij=(ngci+1)*(ngcj+1)
      ngckl=(ngck+1)*(ngcl+1)
      ngcd =ngcij*ngckl
c----------------------------------------------------------------------
c Memory reserved for the WHOLE super-block in prec2ij,prec2kl :
c
      mem2ij=3*ijpar*lcij
      mem2kl=3*klpar*lckl
c
      if(kbl.ne.ibl) then
         memory2=mem2ij+mem2kl
      else
         memory2=mem2ij
      endif
c
c This much memory is allocated for a WHOLE super-block , regaldless
c of its futher splitting into a number of small blocks. Thus,
c this is a CONSTANT part of memory with respect to the small blocks
c
c----------------------------------------------------------------------
c Memory for ONE pair ij, ONE pair kl and ONE quartet IJKL in a S-block
c
      ijpar1=1
      klpar1=1
      nbls1=1
c
      nfumax=nfu(mmax)
      if(nsij.ge.nskl) then
         nfha=nfumax*klpar1*lckl
      else
         nfha=nfumax*ijpar1*lcij
      endif
c----------------------------------------------------------------------
c pair-precalculations:
c
      call in5a(ijpar1,mmax1, memprij)
      call in5b(klpar1,mmax1, memprkl)
c
        if(nsij.ge.nskl) then
           memprkl=memprkl + nfha*3
        else
           memprij=memprij + nfha*3
        endif
c------------------------
      memor2ij=memprij
      memor2kl=memprkl
c
      if(where.eq.'shif') then
c        add memory reserved in memo6 routine :
         memor2ij=memor2ij+3
         memor2kl=memor2kl+3
      endif
c------------------------
c quartet-calculations:
c
      call in5c(nbls1,mmax1,nfumax,mempre4)
c
c------------------------
c quartet-trobsa,assemble:
c
      call in4a(nbls1,memasse,memtrob,where)
c------------------------
c quartet-amshift:
c
      call in4b(nbls1,memamsh,where)
c------------------------
c convertion in erinteg_1 for assembling (iaax,ibbx,iccx) (deriv only)
c
      if(ifor.eq.3 .or. ifor.eq.4) memasse=memasse+3
c
      memasse=memasse+3            ! convert in assemble
      memamsh=memamsh+6*n_times    ! convert in amshift
      if(ngcd.gt.1) memasse=memasse+(ngcij+ngckl+ngcd*lnij*lnkl)
c----------------------------------------------------------------------
c What is needed at once ? Memory for :
c
c    1. Precalc4 + Trobsa+Assemble     proportional to NBLS
c    2. Precalc4 + Assemble+Amshift    proportional to NBLS
c----------------------------------------------------------------------
c
      mem1 = mempre4+memtrob+memasse
      mem2 = mempre4+memasse+memamsh
      memory4 =max(mem1,mem2)
c
c     if(ifor.eq.3) then
c      write(6,*)'   size1: mem1: mempre4+memtrob+memasse=',
c    *               mem1, mempre4,memtrob,memasse
c      write(6,*)'   size1: mem2: mempre4+memasse+memamsh=',
c    *               mem2, mempre4,memasse,memamsh
c     endif
c----------------------------------------------------------------------
c final output is:
c
c     memory2 =  constatnt part of mem. for a small block
c     memor2ij=  mem. for one ij-pair
c     memor2kl=  mem. for one kl-pair
c     memory4 =  mem. for one quartet
c
c formula for memory needed for one small block is as follows:
c
c   memory=memory2 + nparij*memor2ij + nparkl*memor2kl+ nqrt*memor4
c----------------------------------------------------------------------
c
      end
c===============================================================
      subroutine blksize2(ibl,kbl,ijpar,klpar,ityp,jtyp,ktyp,ltyp,
     *                    lcix,lcjx,lckx,lclx, ngcix,ngcjx,ngckx,ngclx,
     *                    memory2,memor2ij,memor2kl,memory4,ifor)
      implicit real*8 (a-h,o-z)
      character*4 where
c----------------------------------------------------------------------
c  Calculates a memory needed in a Super-block for :
c   1) prec2ij and prec2kl (called for each super-block)
c                                     - memory2
c
c  for each block in a super-block calculates :
c
c   2) memory needed for ONE pair ij  - memor2ij
c   3) memory needed for ONE pair kl  - memor2kl
c   4) memory needed for ONE quartet  - memory4
c
c----------------------------------------------------------------------
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c----------------------------------------------------------------------
      common /logic4/ nfu(1)
      common /logic1/ ndege(1)
      common /logic2/ len(1)
      common /logic3/ lensm(1)
c----------------------------------------------------------------------
c ifor=1   for ordinary two-el.integrals
c ifor=2   for NMR GIAO integral derivatives
c ifor=3   for gradient integral derivatives (first  derivatives)
c ifor=4   for hessian  integral derivatives (second derivatives)
c----------------------------------------------------------------------
c
      n_times=1
      if(ifor.eq.2) n_times= 6 ! giao Ist derivatives
      if(ifor.eq.3) n_times= 9 ! gradient derivatives
      if(ifor.eq.4) n_times=45 ! hessian  derivatives
c
      where='    '
      if(ifor.eq.2) where='shif'
      if(ifor.eq.3) where='forc'
      if(ifor.eq.4) where='hess'
c----------------------------------------------------------------------
      lci=lcix
      lcj=lcjx
      lck=lckx
      lcl=lclx
      lcij=lci*lcj
      lckl=lck*lcl
      ngci=ngcix
      ngcj=ngcjx
      ngck=ngckx
      ngcl=ngclx
c----------------------------------------------------------------------
c set up commons /obarai/ and /shell/
c for routines in4a,in4b
c
      NQI=NDEGE(ITYP)
      NQJ=NDEGE(JTYP)
      NQK=NDEGE(KTYP)
      NQL=NDEGE(LTYP)
      NSIJ=NQI+NQJ-1
      NSKL=NQK+NQL-1
c
      if(where.eq.'shif'.or. where.eq.'forc') then
        nsij=nsij+1
        nskl=nskl+1
      endif
      if(where.eq.'hess') then
        nsij=nsij+2
        nskl=nskl+2
      endif
C
      MMAX=NSIJ+NSKL-1
      mmax1=mmax-1
C
      LNI=LEN(ITYP)
      LNJ=LEN(JTYP)
      LNK=LEN(KTYP)
      LNL=LEN(LTYP)
      LNIJKL=LNI*LNJ*LNK*LNL
C
      LNIJ=LENSM(NSIJ)
      LNKL=LENSM(NSKL)
C
       NSIJ1=NSIJ+1
       NSKL1=NSKL+1
       NQIJ=NQI
       IF(NQJ.GT.NQI) NQIJ=NQJ
       NQIJ1=NQIJ+1
       NQKL=NQK
       IF(NQL.GT.NQK) NQKL=NQL
       NQKL1=NQKL+1
c
       LSHELLT=0
       IF(ITYP.EQ.3) LSHELLT=LSHELLT+1
       IF(JTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(KTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(LTYP.EQ.3) LSHELLT=LSHELLT+1
c----------------------------------------------------------------------
c for new general contraction handling :
c
      ngcij=(ngci+1)*(ngcj+1)
      ngckl=(ngck+1)*(ngcl+1)
      ngcd =ngcij*ngckl
c----------------------------------------------------------------------
c Memory reserved for the WHOLE super-block in prec2ij,prec2kl :
c
      mem2ij=(ijpar+2)*lcij
      mem2kl=(klpar+2)*lckl
c
      if(kbl.ne.ibl) then
         memory2=mem2ij+mem2kl
      else
         memory2=mem2ij
      endif
c
c----------------------------------------------------------------------
c Memory for the ONE pair ij,kl and ONE quartet in a Super-block is :
c----------------------------------------------------------------------
c pair-precalculations (as reserved in memo5a,b):
c
c proportional to the number of pairs :
c
      memprij=3+13*lcij
      memprkl=3+13*lckl
c
      if(ifor.eq.3 .or. ifor.eq.4) then
         memprij=memprij+lci+lcj
         memprkl=memprkl+lck
      endif
c
c
c constat part (not proportional to pairs):
c
      memcoij=lci+lcj+lcij + mmax1*lcij
      memcokl=lck+lcl+lckl + mmax1*lckl
      if(ngcd.gt.1) then
        memcoij=memcoij+ngcij*lcij
        memcokl=memcokl+ngckl*lckl
      endif
      memcons=memcoij+memcokl
c------------------------
      memor2ij=memprij
      memor2kl=memprkl
c
      if(where.eq.'shif') then
c      add memory reserved in memo6 routine :
         memor2ij=memor2ij+3
         memor2kl=memor2kl+3
      endif
c------------------------
c quartet-calculations (reserved in memo5c :
c
c proportional to block-szie (nbls):
c
      mempre4=20
      if(ngcd.gt.1) mempre4=20+ngcd+1
c
c constant (not prop. to nbls):
c
      nfumax=nfu(mmax)
      nfha=nfumax*max(lcij,lckl)
      memcons=memcons+(4*(ngcd+1)+3*nfha)
c------------------------
c add constant part of requested memory to memory2
c
      memory2=memory2 + memcons
c------------------------
c quartet-trobsa,assemble:
c
      nbls1=1
c
      call in4a(nbls1,memasse,memtrob,where)
c------------------------
c quartet-amshift:
c
      call in4b(nbls1,memamsh,where)
c------------------------
c convertion in erinteg_2 for assembling (iaax,ibbx,iccx) (deriv only)
c
      if(ifor.eq.3 .or. ifor.eq.4) memasse=memasse+3
c
      memamsh=memamsh+6*n_times   ! convert in amshift
      if(ngcd.gt.1) memasse=memasse+(ngcij+ngckl+ngcd*lnij*lnkl)
c----------------------------------------------------------------------
c What is needed at once ? Memory for :
c
c    1. Precalc4 + Trobsa+Assemble     proportional to NBLS
c    2. Precalc4 + Assemble+Amshift    proportional to NBLS
c----------------------------------------------------------------------
      mem1 = mempre4+memtrob+memasse
      mem2 = mempre4+memasse+memamsh
      memory4 =max(mem1,mem2)
c----------------------------------------------------------------------
c final output is:
c
c     memory2 =  constatnt part of mem. for a small block
c     memor2ij=  mem. for one ij-pair
c     memor2kl=  mem. for one kl-pair
c     memory4 =  mem. for one quartet
c
c formula for memory needed for one small block is as follows:
c
c   memory=memory2 + nparij*memor2ij + nparkl*memor2kl+ nqrt*memor4
c----------------------------------------------------------------------
      end
c===============================================================
c This set of routines named in... is called from
c the blksize1 subroutine in oreder to estimate
c the memory request for a given Super-block
c and use it latter to set up the block-size .
c==============================================================
      subroutine in4a(nbls,memasse,memtrob,where)
      character*4 where
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /logic1/ ndege(1)
      common /logic2/ len(1)
      common /logic3/ lensm(1)
      common /logic4/ nfu(1)
c------------------------------------------
c Memory requested in the Memo4a subroutine
c------------------------------------------
      memtrob=0
      memasse=0
c------------------------------------------
c for trobsa and  for assemble :
c
      mem0=lnij*lnkl
c ------------------------------------------
c       Memory for "assemble"
c ------------------------------------------
      ngcijkl=(ngci+1)*(ngcj+1)*(ngck+1)*(ngcl+1)
      nblsg=nbls*ngcijkl
c
ccccc if(where.ne.'shif' .or. where.ne.'forc') then
c       ----------------------------
        memasse=nblsg*(lnijkl + mem0)
c       ----------------------------
ccccc endif
      if(where.eq.'shif') then
c       ----------------------------
        memasse=nblsg*(7*lnijkl + mem0+6*nfu(nsij)*nfu(nskl)  )
c       ----------------------------
      endif
      if(where.eq.'forc') then
c       ----------------------------
change  memasse=nblsg*(9*lnijkl +4*mem0+10*nfu(nsij)*nfu(nskl)  )
        memasse=nblsg*max(9*lnijkl,4*mem0)
        memasse=memasse + 10*nfu(nsij)*nfu(nskl)
c       ----------------------------
      endif
      if(where.eq.'hess') then
c       ----------------------------
c       belove 54 if grad. are returned to
ccc     memasse=nblsg*max(45*lnijkl,10*mem0)  ! only second returned
        memasse=nblsg*max(54*lnijkl,10*mem0)  ! first & second returned
        memasse=memasse + 55*nfu(nsij)*nfu(nskl)
c       ----------------------------
      endif
c
      if(mmax.le.2) then
         memasse=memasse+2*nbls
         return
      endif
c
        IF(LSHELLT.GT.0) THEN
           mbfkl12=lnij*nfu(nqkl+1)*nbls
           mbfij12=nfu(nqij+1)*lnkl*nbls
c
          if(where.eq.'shif') then
           mbfkl12=lnij*nfu(nqkl1+1)*nbls + 6*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12=nfu(nqij1+1)*lnkl*nbls + 6*nfu(nqij+1)*nfu(nskl)*nbls
          endif
          if(where.eq.'forc') then
           mbfkl12= 4*lnij*nfu(nqkl1+1)*nbls
     *            +10*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12= 4*nfu(nqij1+1)*lnkl*nbls
     *            +10*nfu(nqij+1)*nfu(nskl)*nbls
          endif
c
          if(lshellt.gt.1) then
c           ----------------------
            memasse=memasse+2*(mbfij12+mbfkl12)
c           ----------------------
          else
c           ----------------------
            memasse=memasse+(mbfij12+mbfkl12)
c           ----------------------
          endif
c
        IF( LSHELLT.GT.1 ) THEN
c
          mbf2l =nfu(nqij+1)*nfu(nqkl+1)*nbls
          mbfkl3=lnij*nbls
          mbfij3=lnkl*nbls
          if(where.eq.'shif') then
            mbf2l=nfu(nqij1+1)*nfu(nqkl1+1)*nbls
     *         +6*nfu(nqij +1)*nfu(nqkl +1)*nbls
            mbfkl3=lnij*4*nbls + 6*nfu(nsij)*nbls
            mbfij3=4*lnkl*nbls + 6*nfu(nskl)*nbls
          endif
          if(where.eq.'forc') then
            mbf2l=4*nfu(nqij1+1)*nfu(nqkl1+1)*nbls
     *         +10*nfu(nqij +1)*nfu(nqkl +1)*nbls
            mbfkl3=4*lnij*4*nbls +10*nfu(nsij)*nbls
            mbfij3=4*4*lnkl*nbls +10*nfu(nskl)*nbls
          endif
c
          if(lshellt.gt.2) then
c           ----------------------
            memasse=memasse+4*mbf2l
c           ----------------------
          else
c           ----------------------
            memasse=memasse+2*mbf2l
c           ----------------------
          endif
c           ----------------------
            memasse=memasse+(mbfij3+mbfkl3)
c           ----------------------
c
        IF( LSHELLT.GT.2 ) THEN
c
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=mbf3l0*nbls
          if(where.eq.'shif') then
            mbf3l1=max( nfu(nqij1+1),nfu(nqkl1+1) )
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=4*mbf3l1*nbls + 6*mbf3l0*nbls
          endif
          if(where.eq.'forc') then
            mbf3l1=max( nfu(nqij1+1),nfu(nqkl1+1) )
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=(4*4*mbf3l1 +10*mbf3l0)*nbls
          endif
c
          if(lshellt.gt.3) then
c           ----------------------
            memasse=memasse+4*mbf3l
c           ----------------------
           else
c           ----------------------
            memasse=memasse+2*mbf3l
c           ----------------------
           endif
c
        IF( LSHELLT.GT.3 ) then
c
          i4s =nbls
          if(where.eq.'shif') i4s =16*nbls + 6*nbls
          if(where.eq.'forc') i4s =4*16*nbls +10*nbls
c           ----------------------
            memasse=memasse+i4s
c           ----------------------
        ENDIF
        ENDIF
        ENDIF
        ENDIF
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Memory handling for Obara-Saika-Tracy method
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
c     l11=mmax
c     l12=lensm(mmax)
c     mem1=l11*l12
      mem1=mmax*lensm(mmax)
cc
      mem2=0
      if(nsij.ge.nskl) then
        klstep=0
        do 10 ijstep=mmax,nsij,-1
        klstep=klstep+1
        ijdim=lensm(ijstep)
        kldim=lensm(klstep)
        ijkld=ijdim*kldim
        mem2=mem2+ijkld
   10   continue
      else
        ijstep=0
        do 11 klstep=mmax,nskl,-1
        ijstep=ijstep+1
        ijdim=lensm(ijstep)
        kldim=lensm(klstep)
        ijkld=ijdim*kldim
        mem2=mem2+ijkld
   11   continue
      endif
c           ----------------------
            memtrob=nbls*(mem0+mem1+mem2)
c           ----------------------
c
      end
c
c==============================================================
      subroutine in4b(nbls,memamsh,where)
      character*4 where
c--
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c------------------------------------------
c Memory requested in the Memo4b subroutine
c------------------------------------------
c       Memory for amshift
c
            mwvus=max(lnij,lnkl)*max(nfu(nqj+1),nfu(nql+1))
            mxij=nfu(nqi+1)*nfu(nqij+1)*lnkl
c
            mwij=mwvus
            mwij=mwij*nbls
            mxij=mxij*nbls
        if(where.eq.'shif') then
            mwij=6*mwij
            mxij=6*mxij
        endif
        if(where.eq.'forc') then
            mwij=10*mwij
            mxij=10*mxij
        endif
c           ----------------------
            memamsh=mwij+mxij
c           ----------------------
        IF(LSHELLT.GT.0) THEN
c
            mvus=mwvus
            myz=nfu(nqi+1)*nfu(nqj+1)*nfu(nqkl+1)
            mvus=mvus*nbls
            myz=myz*nbls
        if(where.eq.'shif') then
            mvus=6*mvus
            myz =6*myz
        endif
        if(where.eq.'forc') then
            mvus=10*mvus
            myz =10*myz
        endif
c           ----------------------
            memamsh=memamsh+(mvus+myz)
c           ----------------------
c
        IF( LSHELLT.GT.1 ) THEN
            mbf2l=nfu(nqij+1)*nfu(nqkl+1) *nbls
            if(where.eq.'shif') then
               mbf2l=6*mbf2l
            endif
            if(where.eq.'forc') then
               mbf2l=10*mbf2l
            endif
c           ----------------------
            memamsh=memamsh+(2*mvus+myz)
c           ----------------------
c
          if(lshellt.gt.2) then
c           ----------------------
            memamsh=memamsh+4*mbf2l
c           ----------------------
          else
c           ----------------------
            memamsh=memamsh+mbf2l
c           ----------------------
          endif
c
        IF( LSHELLT.GT.2 ) THEN
c
         mnbls=nbls
         if(where.eq.'shif') mnbls=6*nbls
c
         if(lshellt.gt.3) then
c           ----------------------
            memamsh=memamsh+4*mnbls
c           ----------------------
          else
c           ----------------------
            memamsh=memamsh+2*mnbls
c           ----------------------
          endif
c
        ENDIF
        ENDIF
        ENDIF
c
      end
c==============================================================
      subroutine in5a(npij,mmax1, memory)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory requested in the Memo5a subroutine
c------------------------------------------
c Memory handling for left-hand pairs:
c Total number of calls of Getmem is 12 or 14 (if gen.con.)
c------------------------------------------
      ijpar=npij
c------------------------------------------
       ndi=   ijpar*lci
       ndj=   ijpar*lcj
c     ---------------------------------
      memory=2*(ndi+ndj) + 3*ijpar
c     ---------------------------------
       ndij =ndi*lcj
       ndij3=ndij*3
c     ---------------------
      memory=memory+3*ndij3
c     ---------------------
      memory=memory+2*ndij
c     ---------------------
      memory=memory+ndij3
c     ---------------------
      ndijm=ndij*mmax1
c     ---------------------
      memory=memory+ndijm
c     ---------------------
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
      if(ngcd.gt.1) then
        ndig=ndi*ngci1
        ndjg=ndj*ngcj1
c       -----------------------
        memory=memory+ndig+ndjg
c       -----------------------
      endif
c
      end
c==============================================================
      subroutine in5b(npkl,mmax1, memory)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory requested in the Memo5b subroutine
c------------------------------------------
c Total number of calls of Getmem is 12 or 14 (if gen.con.)
c------------------------------------------
      klpar=npkl
c------------------------------------------
       ndk=   klpar*lck
       ndl=   klpar*lcl
c     ---------------------------------
      memory=2*(ndk+ndl) + 3*klpar
c     ---------------------------------
       ndkl=ndk*lcl
       ndkl3=ndkl*3
c     ---------------------------------
      memory=memory+3*ndkl3
c     ---------------------------------
      memory=memory+2*ndkl
c     ---------------------------------
      memory=memory+ndkl3
c     ---------------------------------
      ndklm=ndkl*mmax1
c     ---------------------------------
      memory=memory+ndklm
c     ---------------------------------
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
      if(ngcd.gt.1) then
        ndkg=ndk*ngck1
        ndlg=ndl*ngcl1
c       -------------------------------
        memory=memory+ndkg+ndlg
c       -------------------------------
      endif
c
      end
c==============================================================
      subroutine in5c(nbls,mmax1,nfumax,memor4)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c memory requested in the Memo5c subroutine
c------------------------------------------
c memory handling
c Total number of calls of Getmem is 24 or 26 (if gen.cont)
c reserve memory for quartets ijkl
c------------------------------------------
      nblsi=nbls
      if(intsize.ne.1) nblsi=nbls/intsize+1
c     ------------------------------------
      memor4=4*nblsi
c     ------------------------
      memor4=memor4 + 4*nbls
c     ------------------------
      nbmx=nbls*mmax1
c     ------------------------
      memor4=memor4 + 2*nbmx
c     ------------------------
      nbls3=nbls*3
c     ------------------------
      memor4=memor4 + 5*nbls3 + 3*nbls
c     ------------------------
      memor4=memor4 + nbls3*nfumax
c     ------------------------
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c     ------------------------
c-->  memor4=memor4 + 4*ngcd
c     ------------------------
      if(ngcd.gt.1) then
c       ----------------------
        memor4=memor4 + nbls*(1+ngcd)
c       ----------------------
      endif
c
      end
c==============================================================
c symmetry :
c
      subroutine parsym(ncs,inx,iis,jjs,ijshp,nsym)
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      dimension inx(7,ncs)
      dimension iis(*),jjs(*),ijshp(nsym,*)
      logical firstd
c
      ncsp=ncs*(ncs+1)/2
c
c     write(8,*)'          pairs from parsym       '
c     write(8,*)' pair  i,j-shells   symm.op.  image i,j-shells '
c
      do 100 ijcs=1,ncsp
         ics=iis(ijcs)
         jcs=jjs(ijcs)
c
c find images of ijcs for each symm. oper.
c
          do 150 ns=1,nsym
             ics0=inx(ns,ics)
             jcs0=inx(ns,jcs)
c canonical order for images
             icsi=max0(ics0,jcs0)
             jcsi=min0(ics0,jcs0)
c find what pair it is
c
             ijshp(ns,ijcs)=icsi*(icsi-1)/2 + jcsi
c
c            write(8,80) ijcs,ics,jcs, ns,ijshp(ns,ijcs),ics0,jcs0
  150     continue
  100 continue
c
c  80 format(i5,3x,2(i3,1x),5x,i3,5x,i5,3x,2(i3,1x))
c
c---------------------------------------------------
c     if( nprint.gt.2) then
c       write(8,*)' PAIR and its IMAGES'
c       do 250 ijcs=1,ncsp
c       write(8,90) ijcs,(ijshp(ns,ijcs),ns=1,nsym)
c 250   continue
c  90   format(i4,7x,3(i4,2x))
c     endif
c---------------------------------------------------
      end
c===============================================================
      subroutine price(inx,iis,jjs,ijbl,npar,nsupb,nbl2,nprint,
     *                 ncost_int,minpr,maxpr,
     *                 xcost_blk,xminpr_bl,xmaxpr_bl)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /logic1/ ndege(1)
      common /logic2/ lenn(1)
      common /logic3/ lensm(1)
      common /logic4/ nfu(1)
c
      dimension nsupb(*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*),npar(*)
      dimension ncost_int(*)
      dimension xcost_blk(*)
c-------------------------
c*
      if(nprint.gt.3) then
         call getival('iout',iout)
         write(iout,505)
  505 format(
     *' Block  MMAX   type     contr    nfun/q      Nqrt   Nfunc',2x,
     * 'Price/i','    Price/blk')
      endif
c
c contracted integrals counter
c
      xmaxpr_bl=0.d0          ! maximum price for the whole block
      xminpr_bl=1.d+20        ! minimum price for the whole block
c
      maxpr=0              ! maximum price for 1 integral in a block
      minpr=1000000        ! minimum price for 1 integral in a block
      maxst=0
      ntotq=0
c
      ikbl=0
      do 100 ibl=1,nbl2
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
c     int=inx(3,ics1)
c     jnt=inx(3,jcs1)
c     ijnt=int*jnt
      itype1=itype
      jtype1=jtype
      if(itype.gt.4) itype1=itype-1
      if(jtype.gt.4) jtype1=jtype-1
      if(itype1.gt.5) itype1=itype1-1
      if(jtype1.gt.5) jtype1=jtype1-1
      if(itype.eq.11) itype1=6
      if(itype.eq.12) itype1=7
      if(jtype.eq.13) jtype1=8
      if(jtype.eq.11) jtype1=6
      if(jtype.eq.12) jtype1=7
      if(itype.eq.13) itype1=8
      icont=inx(5,iis(ijcs1))-inx(1,iis(ijcs1))
      jcont=inx(5,jjs(ijcs1))-inx(1,jjs(ijcs1))
      ijcont=icont*jcont
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
      ngci1=ngci+1
      ngcj1=ngcj+1
      nparij=npar(ibl)
      nfij  =lenn(itype1)*lenn(jtype1)
      ngcij =ngci1*ngcj1
c-
         do 100 kbl=1,ibl
         ikbl=ikbl+1
         klcs1=ijbl(kbl,1)
         kcs1=iis(klcs1)
         lcs1=jjs(klcs1)
         nparkl=npar(kbl)
c
         nquart=nparij*nparkl
         if(kbl.eq.ibl) nquart=nparij*(nparij+1)/2
c
         ntotq=ntotq+nquart
c
c        knt=inx(3,kcs1)
c        lnt=inx(3,lcs1)
c        klnt=knt*lnt
         ktype=inx(12,kcs1)
         ltype=inx(12,lcs1)
         ktype1=ktype
         ltype1=ltype
         if(ktype.gt.4) ktype1=ktype-1
         if(ltype.gt.4) ltype1=ltype-1
         if(ktype1.gt.5) ktype1=ktype1-1
         if(ltype1.gt.5) ltype1=ltype1-1
         if(ktype.eq.11) ktype1=6
         if(ktype.eq.12) ktype1=7
         if(ktype.eq.13) ktype1=8
         if(ltype.eq.11) ltype1=6
         if(ltype.eq.12) ltype1=7
         if(ltype.eq.13) ltype1=8
c
         nfijkl=nfij*lenn(ktype1)*lenn(ltype1)
         nfun=nquart*nfijkl
c
         kcont=inx(5,iis(klcs1))-inx(1,iis(klcs1))
         lcont=inx(5,jjs(klcs1))-inx(1,jjs(klcs1))
         ijklc=ijcont*kcont*lcont
c
c number of general contractions
c
         ngck=inx(4,kcs1)
         ngcl=inx(4,lcs1)
         ngck1=ngck+1
         ngcl1=ngcl+1
cccccc   ngcd=ngci1*ngcj1*ngck1*ngcl1
         ngcd=ngcij*ngck1*ngcl1
c-
         nqi=ndege(itype1)
         nqj=ndege(jtype1)
         nsij=nqi+nqj-1
         lnij=lensm(nsij)
         nqij=nqi
         if(nqj.gt.nqi) nqij=nqj
c
         nqk=ndege(ktype1)
         nql=ndege(ltype1)
         nskl=nqk+nql-1
         lnkl=lensm(nskl)
         nqkl=nqk
         if(nql.gt.nqk) nqkl=nql
c
         mmax=nsij+nskl-1
c
c estimate the price for a block
c
         if ( mmax.le.2 ) then
           nprice=39
           if(mmax.eq.2) nprice=54
           if( max(itype,jtype,ktype,ltype).eq.3) nprice=56
         else
c from prec4neg(36) and XWPQ(15) :
           nprice=51
           nprice=nprice+4*mmax -3
c
c from Obasai :
c
           isum=0
           do 851 im=1,mmax-2
           isum=isum + (mmax-1-im)*( nfu(im+1)-nfu(im) )
  851      continue
           nfact=30-mmax
           if(mmax.ge.9) nfact=18
           nprice=nprice+nfact*isum
c
           nsxx=nskl
           nqxx=nqij
           if(nskl.gt.nsij) then
             nsxx=nsij
             nqxx=nqkl
           endif
c
c from Tracy :
c
           isum=0
           do 852 kp=2,nsxx
           ndix=nqxx-nsxx+kp
           if(ndix.le.0) ndix=1
           isum=isum+(nfu(mmax+1-kp)-nfu( ndix ))*(nfu(kp+1)-nfu(kp))
  852      continue
           nprice=nprice+7*isum
c
c from assemble ;
c
           nassem=(lnkl-nfu(nqkl))*(lnij-nfu(nqij))
c
           nprice=nprice+nassem
c
           if(itype.eq.3) nprice=nprice+2*nassem
           if(jtype.eq.3) nprice=nprice+2*nassem
           if(ktype.eq.3) nprice=nprice+2*nassem
           if(ltype.eq.3) nprice=nprice+2*nassem
c
         endif
c
ccccccc  nprice=nprice*icont*jcont*kcont*lcont
         nprice=nprice*ijklc
c
c from shifting of angular momentum :
c
         if(nqij.eq.nsij .and. nqkl.eq.nskl) then
         else if (mmax.gt.2) then
            nq1=nqj
            nq2=nqi
               if(nqj.gt.nqi) then
                  nq1=nqi
                  nq2=nqj
               endif
            isum=0
            do 853 j=2,nq1
            do 853 i=nsij+1-j,nq2,-1
            isum=isum+(nfu(j+1)-nfu(j))*(nfu(i+1)-nfu(i))
  853       continue
            isum=2*isum
c from tfer
            isum=isum*(nfu(nskl+1)-nfu(nqkl))
c
            nq1=nql
            nq2=nqk
               if(nql.gt.nqk) then
                  nq1=nqk
                  nq2=nql
               endif
            ksum=0
            do 854 l=2,nq1
            do 854 k=nskl+1-l,nq2,-1
            ksum=ksum+(nfu(l+1)-nfu(l))*(nfu(k+1)-nfu(k))
  854       continue
            ksum=2*ksum
c from tfer
            ksum=ksum*(nfu(nqi+1)-nfu(nqi))*(nfu(nqj+1)-nfu(nqj))
c
            nprice=nprice + isum+ksum
c
         endif
c
        nprice=nprice/nfijkl
c
        ngen1=16*ngcd
        ngen2=(ngci1+1)*(ngcj1+1)*(ngck1+1)*(ngcl1+1)
        ngen12=ngen1/ngen2
c
        nprice=nprice*ngen12
        nprice=nprice/ngcd
c
        xprice=dble(nprice)*dble(nfun)
c
        ncost_int(ikbl)=nprice       ! price for 1 integral in a block
        xcost_blk(ikbl)=xprice       ! price for the whole block
c
        if(nprice.gt.maxpr) maxpr=nprice
        if(nprice.lt.minpr) minpr=nprice
c
        if(xprice.gt.xmaxpr_bl) xmaxpr_bl=xprice
        if(xprice.lt.xminpr_bl) xminpr_bl=xprice
c
c  end of price
c
       if(nprint.gt.3) then
         if(mmax.ne.maxst) then
c* for printing purpose only
           write(iout,*)'   '
           maxst=mmax
         endif
         write(iout,506) ikbl,mmax,itype,jtype,ktype,ltype,
     *                icont,jcont,kcont,lcont,nfijkl,nquart,nfun,
     *                nprice,xprice
         if(nsupb(ikbl).gt.1) then
            write(iout,507) nsupb(ikbl)
         endif
       endif
  100 continue
c
  506 format(2x,i4,2x,i2,2x,4i2,2x,4i2,3x,i5,4x,i6,1x,i8,3x,i6,2x,f12.0)
  507 format(5x,' this super-block is split into ',i5,' blocks')
c---------------------------------------------------------------
c
      end
c================================================================
c    subroutines for blocking quartets :
c================================================================
      subroutine blockin4(isbl,ibl,kbl,nbloks,mxsize,
     *                    bl,inx,npard,nbl2)
c---------------------------------------------
c this routine is called for each super-block
c NBLOKS is a number of small blocks belonging
c to the given super-block ISBL .
c---------------------------------------------
      implicit real*8 (a-h,o-z)
      logical firstd
c
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /memor1/ iisd,jjsd,ijbld
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
      common /memor3/ nblok1d
      common /memors/ nsym,ijshp,isymm
c
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------
c  constructe blocks of quartets of contracted shells :
c  set up vectors : nqrt,nibl,nkbl, nijb,nije, nklb,nkle
c
      call memo2(nbloks)  ! res.mem. for 7 arrays (above)
c
      call blockqur(isbl,ibl,kbl,nbloks,maxqrt,nbl2,bl(ijbld),bl(npard),
     *bl(nqrtd),bl(nibld),bl(nkbld),bl(nijbd),bl(nijed),
     *bl(nklbd),bl(nkled),bl(iisd),bl(jjsd),inx,   mxsize )
c
      call memo3(maxqrt)
c
c end of the blocking procedure
c
      end
c===============================================================
      subroutine blockqur(isbl,ibl,kbl,nbloks,maxqrt,nbl2,ijbl,npar,
     * nqrt,nibl,nkbl, nijb,nije,nklb,nkle,iis,jjs,inx,mxsize)
c
cccc  dimension mxsize(*)
      dimension npar(*),ijbl(nbl2,*)
      dimension inx(12,*),iis(*),jjs(*)
c
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension nqrt(*)
c-----------------------------------------------------------
c This is the super-block ISBL  and it is constructed
c from pair-blocks ibl and kbl .
c
c Construct NBLOKS blocks of quartets of contracted shells
c belonging to the super-block ISBL
c-----------------------------------------------------------
c input :
c
c isbl - given super-block to be split in nbloks small parts
c ibl,kbl - pair-blocks making this super-block
c nbloks - number of small blocks isbl will be split into
c          known from blockin2 (pairs blocking)
c nbl2 - number of pair-blocks
c ijbl(ibl,ijp)-> ijcs pair of contr.shells
c
c npar(ibl) - number of pairs in ibl pair-block
c nqrt(ikbl)- number of quartets in the small block IKBL
c
c Output :
c
c maxqrt - maximum number of quartets in one block
c-----------------------------------------------------------
c
      ikbl=0
c
      nbls=mxsize
      if(ibl.gt.kbl) then
         call nondiax(nbls,nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *                nibl,nkbl, nijb,nije,nklb,nkle)
      else
         call diagonx(nbls,nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *                nibl,nkbl, nijb,nije,nklb,nkle)
      endif
c
c It should be  nbloks=ikbl
c
      if(ikbl.ne.nbloks) then
c1999    write(6,*)' In the Super-block no.=',isbl
c1999    write(6,*)' number of blocks of contracted quartets'
c1999    write(6,*)' from Blksizer and Blockqur subroutines'
c1999    write(6,*)'          is not the same             '
c1999    write(6,*)' should be=',nbloks,' it is =',ikbl
c
c       do iix=ikbl+1,nbloks
c          nibl(iix)=0
c          nkbl(iix)=0
c          nijb(iix)=0
c          nije(iix)=0
c          nklb(iix)=0
c          nkle(iix)=0
c          nqrt(iix)=0
c       enddo
c
         call nerror(isbl,'blksizer',
     *               ' wrong number of blocks ',ikbl,nbloks)
      endif
c
      maxqrt=0
      ikblmx=1
      ntotq=0
      do 650 ikbl=1,nbloks
      ntotq=ntotq+nqrt(ikbl)
          if(nqrt(ikbl).gt.maxqrt) then
             maxqrt=nqrt(ikbl)
             ikblmx=ikbl
          endif
  650 continue
c------------------------------------------------------------------
c     write(8,497) isbl,nbloks,ntotq,maxqrt,ikblmx
c 497 format(/'     In the Super-block no.=',i6/
c    *        '  number of Blocks of Contracted Quartets =',i6/
c    *        '  total number of contracted quartets =',i7/
c    *        '  maximum number of quartets =',i6,' in block=',i6/)
c
c------------------------------------------------------------------
ctest
c     maxqrt=0
c     ikblmx=1
c     ntotq=0
c     do 651 ikbl=1,nbloks
c     ntotq=ntotq+nqrt(ikbl)
c         if(nqrt(ikbl).gt.maxqrt) then
c            maxqrt=nqrt(ikbl)
c            ikblmx=ikbl
c         endif
c 651 continue
c     write(8,497) nbloks,ntotq,maxqrt,ikblmx
c------------------------------------------------------------------
      end
c===============================================================
      subroutine nondiax(nbls,nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *nibl,nkbl, nijb,nije,nklb,nkle)
      implicit real*8 (a-h,o-z)
c
      dimension npar(*),ijbl(nbl2,*)
      dimension nqrt(*)
c
      dimension nibl(*), nkbl(*)
      dimension nijb(*),nije(*),nklb(*),nkle(*)
c
      ijpar=npar(ibl)
      klpar=npar(kbl)
      nquart=ijpar*klpar
c
      maxsize=nbls
c
          ijsize=ijpar
          klsize=klpar
          ijdev=1
          kldev=1
          ijrem=0
          klrem=0
c
          if(nquart.le.maxsize) go to 99
c
c case 1
          if( ijpar.gt.maxsize.and.klpar.gt.maxsize) then
c
             ijsize=maxsize
             ijdev=ijpar/ijsize
             ijrem=mod(ijpar,ijsize)
c
             klsize=1
             kldev=klpar/klsize
             klrem=mod(klpar,klsize)
          endif
c case 2
          if( ijpar.gt.maxsize.and.klpar.le.maxsize) then
c
             ijsize=maxsize
             ijdev=ijpar/ijsize
             ijrem=mod(ijpar,ijsize)
c
             klsize=1
             kldev=klpar/klsize
             klrem=mod(klpar,klsize)
          endif
c case 3
          if( ijpar.le.maxsize.and.klpar.gt.maxsize) then
             klsize=maxsize/ijpar
             kldev=klpar/klsize
             klrem=mod(klpar,klsize)
          endif
c case 4
          if( ijpar.le.maxsize.and.klpar.le.maxsize) then
c            if(ijpar.lt.klpar) then
c              ijsize=maxsize/klpar
c              ijdev=ijpar/ijsize
c              ijrem=mod(ijpar,ijsize)
c            else
c              klsize=maxsize/ijpar
c              kldev=klpar/klsize
c              klrem=mod(klpar,klsize)
c            endif
c-----
               klsize=maxsize/ijpar
               kldev=klpar/klsize
               klrem=mod(klpar,klsize)
c-----
          endif
c case 5 (no-blocking: 1 quartet per block )
          if(maxsize.eq.1) then
            ijsize=1
            klsize=1
            ijdev=ijpar
            kldev=klpar
            ijrem=0
            klrem=0
          endif
c
   99 continue
c
          nqrt1=ijsize*klsize
          nqrt2=ijrem*klsize
          nqrt3=ijsize*klrem
          nqrt4=ijrem*klrem
          ijds=ijdev*ijsize
          klds=kldev*klsize
c
          do 100 ij=1,ijdev
          ij1=(ij-1)*ijsize
             do 110 kl=1,kldev
             kl1=(kl-1)*klsize
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ij1+1
             nije(ikbl)=ij1+ijsize
             nklb(ikbl)=kl1+1
             nkle(ikbl)=kl1+klsize
             nqrt(ikbl)=nqrt1
  110        continue
             if(klrem.gt.0) then
                ikbl=ikbl+1
                nibl(ikbl)=ibl
                nkbl(ikbl)=kbl
                nijb(ikbl)=ij1+1
                nije(ikbl)=ij1+ijsize
                nklb(ikbl)=klds+1
                nkle(ikbl)=klpar
                nqrt(ikbl)=nqrt3
             endif
  100     continue
c
           if(ijrem.gt.0) then
             do 200 kl=1,kldev
             kl1=(kl-1)*klsize
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ijds+1
             nije(ikbl)=ijpar
             nklb(ikbl)=kl1+1
             nkle(ikbl)=kl1+klsize
             nqrt(ikbl)=nqrt2
  200        continue
           endif
c
ckw        if(klrem.gt.0) then
ckw          do 300 ij=1,ijdev
ckw          ij1=(ij-1)*ijsize
ckw          ikbl=ikbl+1
ckw          nibl(ikbl)=ibl
ckw          nkbl(ikbl)=kbl
ckw          nijb(ikbl)=ij1+1
ckw          nije(ikbl)=ij1+ijsize
ckw          nklb(ikbl)=klds+1
ckw          nkle(ikbl)=klpar
ckw          nqrt(ikbl)=nqrt3
c 300        continue
ckw        endif
c
           if(nqrt4.gt.0) then
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ijds+1
             nije(ikbl)=ijpar
             nklb(ikbl)=klds+1
             nkle(ikbl)=klpar
             nqrt(ikbl)=nqrt4
           endif
c
      end
c===============================================================
      subroutine diagonx(nbls,nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *nibl,nkbl, nijb,nije,nklb,nkle)
      implicit real*8 (a-h,o-z)
c
      dimension npar(*),ijbl(nbl2,*)
c
      dimension nqrt(*),nibl(*),nkbl(*)
      dimension nijb(*),nije(*),nklb(*),nkle(*)
c
      ijpar=npar(ibl)
c
      nquart=ijpar*(ijpar+1)/2
c
      if(nquart.le.nbls) then
c
         ikbl=ikbl+1
c
         nibl(ikbl)=ibl
         nkbl(ikbl)=kbl
         nijb(ikbl)=1
         nije(ikbl)=ijpar
         nklb(ikbl)=1
         nkle(ikbl)=0
         nqrt(ikbl)=nquart
c
cxxxxxxxxxx
c        ijkl=0
c        do 100 ijp=1,ijpar
c        ijcs=ijbl(ibl,ijp)
c        do 100 klp=1,ijp
c        klcs=ijbl(kbl,klp)
c        ijkl=ijkl+1
c 100    continue
cxxxxxxxxxx
c
      else
c
         ijsize=sqrt( dble(nbls) )
         ndev=ijpar/ijsize
         nrem=mod(ijpar,ijsize)
c
c  ijpar=ndev*ijsize+nrem
c
         nqrtd=ijsize*(ijsize+1)/2
         nqrtn=ijsize*ijsize
c
c  true diagonal blocks
c
         do 200 n=1,ndev
         nn=(n-1)*ijsize
         ikbl=ikbl+1
c
         nibl(ikbl)=ibl
         nkbl(ikbl)=kbl
         nijb(ikbl)=nn+1
         nije(ikbl)=nn+ijsize
         nklb(ikbl)=nn+1
         nkle(ikbl)=0
         nqrt(ikbl)=nqrtd
c
cxxxxxxxxxx
c            ijkl=0
c            do 250 ijp=nn+1,nn+ijsize
c            ijcs=ijbl(ibl,ijp)
c            do 250 klp=nn+1,ijp
c            klcs=ijbl(kbl,klp)
c            ijkl=ijkl+1
c 250        continue
cxxxxxxxxxx
  200    continue
c
c  non-diagonal square blocks ( ndev*(ndev-1)/2 )
c
         do 300 n1=2,ndev
         nn1=(n1-1)*ijsize
         do 300 n2=1,n1-1
         nn2=(n2-1)*ijsize
         ikbl=ikbl+1
c
         nibl(ikbl)=ibl
         nkbl(ikbl)=kbl
         nijb(ikbl)=nn1+1
         nije(ikbl)=nn1+ijsize
         nklb(ikbl)=nn2+1
         nkle(ikbl)=nn2+ijsize
         nqrt(ikbl)=nqrtn
c
cxxxxxxxxxx
c            ijkl=0
c            do 350 ijp=nn1+1,nn1+ijsize
c            ijcs=ijbl(ibl,ijp)
c            do 350 klp=nn2+1,nn2+ijsize
c            klcs=ijbl(kbl,klp)
c            ijkl=ijkl+1
c 350        continue
cxxxxxxxxxx
  300    continue
c
         if(nrem.gt.0) then
c
c  non-diagonal rectangle blocks ( ndev blocks nrem*ijsize )
c
             nqrtr=nrem*ijsize
c
c  for the first
c
             ikbl=ikbl+1
c
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ndev*ijsize+1
             nije(ikbl)=ndev*ijsize+nrem
             nklb(ikbl)=1
             nkle(ikbl)=ijsize
             nqrt(ikbl)=nqrtr
c
cxxxxxxxxxx
c               ijkl=0
c               do 451 ijp=ndev*ijsize+1,ndev*ijsize+nrem
c               ijcs=ijbl(ibl,ijp)
c               do 451 klp=1,ijsize
c               klcs=ijbl(kbl,klp)
c               ijkl=ijkl+1
c 451           continue
cxxxxxxxxxx
c
             nqrtrem=nqrtr
             do 400 n=2,ndev
             nn=(n-1)*ijsize
c         if(nqrtr+nqrt(ikbl).le.nbls ) then
          if(nqrtr+nqrtrem.le.nbls ) then
c            ijkl=nqrt(ikbl)
             nkle(ikbl)=nn+ijsize
             nqrt(ikbl)=nqrtrem+nqrtr
             nqrtrem=nqrtrem+nqrtr
          else
             ikbl=ikbl+1
c            ijkl=0
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ndev*ijsize+1
             nije(ikbl)=ndev*ijsize+nrem
             nklb(ikbl)=nn+1
             nkle(ikbl)=nn+ijsize
             nqrt(ikbl)=nqrtr
             nqrtrem=nqrtr
c
          endif
cxxxxxxx
c               do 450 ijp=ndev*ijsize+1,ndev*ijsize+nrem
c               ijcs=ijbl(ibl,ijp)
c               do 450 klp=nn+1,nn+ijsize
c               klcs=ijbl(kbl,klp)
c               ijkl=ijkl+1
c 450           continue
c            nqrt(ikbl)=ijkl
cxxxxxxx
  400        continue
c
c  last small diagonal block ( nrem*(nrem+1)/2 )
c  it may be a separate block or it may go together with
c  the first true diagonal block, depending on NBLS
c
              ikbl=ikbl+1
c
              nibl(ikbl)=ibl
              nkbl(ikbl)=kbl
              nijb(ikbl)=ndev*ijsize+1
              nije(ikbl)=ndev*ijsize+nrem
              nklb(ikbl)=ndev*ijsize+1
              nkle(ikbl)=0
              nqrt(ikbl)=nrem*(nrem+1)/2
c
              ijkl=0
c
cxxxxxxxxxx
c               do 500 ijp=ndev*ijsize+1,ndev*ijsize+nrem
c               ijcs=ijbl(ibl,ijp)
c               do 500 klp=ndev*ijsize+1,ijp
c               klcs=ijbl(kbl,klp)
c               ijkl=ijkl+1
c 500           continue
cxxxxxxxxxx
         endif
      endif
c
      end
c===============================================================
      subroutine block1(ikbl,nibl,nkbl,nijb,nije,nklb,nkle,ijbl,nbl2,
     *                  nblok1,nqrt)
      implicit real*8 (a-h,o-z)
c
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension ijbl(nbl2,*)
      dimension nblok1(2,*)
      dimension nqrt(*)
c
      ibl=nibl(ikbl)
      kbl=nkbl(ikbl)
      ijbeg=nijb(ikbl)
      ijend=nije(ikbl)
      klbeg=nklb(ikbl)
      klend=nkle(ikbl)
c
c klend can be 0 for "diagonal" blocks
c
      ijkl=0
      if(klend.ne.0) then
         do ijp=ijbeg,ijend
            ijcs=ijbl(ibl,ijp)
            do klp=klbeg,klend
               klcs=ijbl(kbl,klp)
               ijkl=ijkl+1
               nblok1(1,ijkl)=ijcs
               nblok1(2,ijkl)=klcs
            enddo
         enddo
      else
         do ijp=ijbeg,ijend
            ijcs=ijbl(ibl,ijp)
            do klp=klbeg,ijp
               klcs=ijbl(kbl,klp)
               ijkl=ijkl+1
               nblok1(1,ijkl)=ijcs
               nblok1(2,ijkl)=klcs
            enddo
         enddo
      endif
c
      nqrt(ikbl)=ijkl
c
      end
c===============================================================
      subroutine onesym(ibl,kbl,ijbl,nbl2,ijbeg,ijend,klbeg,klend,
     *                  ijshp,nsym,isymm,nblsym )
      dimension ijbl(nbl2,*),ijshp(nsym,*),isymm(*)
      dimension imij(8),imkl(8)
c number 8 is the maximum number of symmetry oper. (nsym)
c
      ijkls=0
      ijkl =0
      do 900 ijp=ijbeg,ijend
         ijcs=ijbl(ibl,ijp)
         klendx=klend
         if(klend.eq.0) klendx=ijp
         do 900 klp=klbeg,klendx
            klcs=ijbl(kbl,klp)
            ijkl=ijkl+1
c
         if(isymm(ijkl).eq.0) go to 900
c canonical order
         ij0=max0(ijcs,klcs)
         kl0=min0(ijcs,klcs)
c find the latest among images of ijcs,klcs for each symm.oper. :
         ijcs1=ijshp(1,ijcs)
         klcs1=ijshp(1,klcs)
         ij1=max0(ijcs1,klcs1)
         kl1=min0(ijcs1,klcs1)
c
         imij(1)=ij1
         imkl(1)=kl1
c
         do 200 ns=2,nsym
            ijcsi=ijshp(ns,ijcs)
            klcsi=ijshp(ns,klcs)
c canonical order
            iji=max0(ijcsi,klcsi)
            kli=min0(ijcsi,klcsi)
c
            imij(ns)=iji
            imkl(ns)=kli
c
            if(iji.gt.ij1) then
               ij1=iji
               kl1=kli
            endif
            if(iji.eq.ij1 .and. kli.gt.kl1) then
               ij1=iji
               kl1=kli
            endif
  200    continue
c check if the image is later or not
         if( ij1.lt.ij0 ) then
            ijkls=ijkls+1
cold        isymm(ijkl)=1
            isymm(ijkl)=-1
         endif
         if( ij1.eq.ij0 .and. kl1.le.kl0) then
            ijkls=ijkls+1
cold        isymm(ijkl)=1
            isymm(ijkl)=-1
         endif
c
         if(isymm(ijkl).eq.-1) then
c           find how many distinct quartets it represents
            call distinct(nsym,imij,imkl,ij0,kl0,ndist)
            isymm(ijkl)=ndist
         else
            isymm(ijkl)=0
         endif
  900 continue
c
      nblsym=ijkls
c
      end
c========================================================
      subroutine distinct(nsym,imij,imkl,ij0,kl0,ndist)
      dimension imij(*),imkl(*)
      dimension ijkl(9)
c
c number 9 comes from nsym+1 where nsym can be up to 8
c
      ij1=imij(1)
      kl1=imkl(1)
c
      if(nsym.eq.1) then
        ndist=1
        if(ij1.ne.ij0 .or. kl1.ne.kl0) ndist=2
        return
      endif
c
      ijkl(1)=ij0*(ij0-1)/2+kl0
      do ns=2,nsym+1
        ij=imij(ns-1)
        kl=imkl(ns-1)
        ijkl(ns)=ij*(ij-1)/2+kl
      enddo
c
      call findif(ijkl,nsym+1,ndist)
c
      end
c=====
      subroutine findif(ijkl,ndim,ndist)
      dimension ijkl(ndim)
c
   10 continue
      iexch=0
      do 20 i=1,ndim-1
      i0=ijkl(i)
      i1=ijkl(i+1)
c
      if(i1.lt.i0) then
        iexch=1
        ijkl(i)=i1
        ijkl(i+1)=i0
      endif
   20 continue
      if(iexch.eq.1) go to 10
c
      iexch=1
      do 30 i=1,ndim-1
      i0=ijkl(i)
      i1=ijkl(i+1)
      if(i0.ne.i1) then
        iexch=iexch+1
      endif
   30 continue
      ndist=iexch
      return
      end
c==========test==============================
      subroutine print_pairs(ibl,kbl,ijbl,nbl2,nparij,nparkl)
      dimension ijbl(nbl2,*)
c
      do ijpar=1,nparij
      ijcs=ijbl(ibl,ijpar)
      write(8,*)' ij-pairs :',ijcs
      enddo
c
      do klpar=1,nparkl
      klcs=ijbl(kbl,klpar)
      write(8,*)' kl-pairs :',klcs
      enddo
c
      end
c================================================================
      subroutine blockin2_again(lcore,nprint,ncs,inx,bl,nbl1,
c output :
     *                    nbl2,minpr,maxpr,  xminpr_bl, xmaxpr_bl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c
      common /memmax/ ispblx, maxme1,iforwhat
      common /route/ iroute
c all below are output adresses :
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1c/ map_ij_bl2
      common /memors/ nsym,ijshp,isymm
c up to here
c
c     common /intbl/ifp,inxx(100)
c
      dimension bl(*)
      dimension inx(12,*)
c---------------------------------------------------------------
c blocking procedure for pairs
c---------------------------------------------------------------
c allocate memory for iis(ijcs)--> ics  & jjs(ijcs)-->jcs arrays
c allocate memory for map_ij_bl2(ijcs)-->ijbl2
c
      ncsp=ncs*(ncs+1)/2
c
      call getint(ncsp,iisd)
      call getint(ncsp,jjsd)
      call getint(ncsp,map_ij_bl2)
c
c Store these addresses in the  common /memor1/ iisd,jjsd,ijbld
c                               common /memor1c/ map_ij_bl2
c---------------------------------------------------------------
c Constructe blocks of contracted shell pairs , setup the arrays:
c npar, ijbl and nsupb
c
        call getival('schwarz',ischwarz)
        call blk_pairs_again(bl,ncs,inx,bl(iisd),bl(jjsd),nbl1,
     *                nbl2,nsupb,npard,ijbld, mxsize,map_ij_bl2,
     *                lcore,nsym,bl(ischwarz))
c
c  output : arrays' addresses : nsupb, npar, ijbl and  mxsize )
c---------------------------------------------------------------
      if(nsym.gt.0) then
         ncsp=ncs*(ncs+1)/2
         call getint(nsym*ncsp, ijshp)
         call getival('SymShPr',ifp1)
         call parsym(ncs,bl(ifp1),bl(iisd),bl(jjsd),bl(ijshp),nsym)
      endif
c---------------------------------------------------------------
c About prices :
c I need two prices : per 1 integral & per one WHOLE big-block
c
      nbl4=nbl2*(nbl2+1)/2
      call getmem(2*nbl4, ncost)
      call price(inx,bl(iisd),bl(jjsd),bl(ijbld),bl(npard),
     *           bl(nsupb),nbl2,nprint,
     *           bl(ncost),minpr,maxpr,
     *           bl(ncost+nbl4+1),xminpr_bl,xmaxpr_bl)
c
c bl(ncost) - prices per integral
c bl(ncost+nbl4+1) - prices for the whole block
c
c minpr & maxpr - minimum & maximum price for one integral
c xminpr_bl & xmaxpr_bl - minimum & maximum price for a block
c---------------------------------------------------------------
c This is the end of blocking procedure for PAIRS
c
c The following is known at this point :
c
c  pairs ijcs are given by ics=iis(ijcs) and jcs=jjs(ijcs)
c  number of pair-blocks = NBL2
c  number of pairs in the block ibl is NPAR(ibl)
c  which pairs belong to this block : ijcs=ijbl(ibl,1-npar)
c---
c The number of BLOCKS of contracted shell Quartets would be
c NBL4=NBL2*(NBL2+1)/2 if they are not too big. I call
c them Super-blocks of c.s.quartets.The maximum size for
c each of them is known and kept in array MAXSIZE(i=1,NSUPB)
c (bl(mxsize)) (no of quartets). According to this max.size
c super-blocks will be split into a number of smaller blocks
c which is also known (NBLOKS) but is NOT transfered.
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairs_again(bl,ncs,inx,iis,jjs,nbl1,
     *                     nbl2,nsupb,npard,ijbld,mxsize,map_ij_bl2,
     *                     lcore,nsym,schwarz)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /memmax/ ispblx, maxme1,iforwhat
      common /route/ iroute
      dimension bl(*)
      dimension inx(12,*)
      dimension iis(*),jjs(*)
      dimension schwarz(ncs,ncs)
c---------------------------------------------------------------
c Construct pairs of contracted shells ( I>=J )
c and setup iis(ijcs) & jjs(ijcs) arrays :
c
      call pairs(ncs,inx,iis,jjs)  ! eliminate it latter on
c---------------------------------------------------------------
c NOTE :
c _cal3 & _mak3 split blocks of pairs in THREE groups (Big,Med.,Small)
c _cal2 & _mak2 split blocks of pairs in TWO   groups (Big,    ,Small)
c _cal1 & _mak1 do not split blocks of pairs in groups(Big,    ,Small)
c---------------------------------------------------------------
c Calculate number of pair-blocks NBL2 & maximum pairs/block MAXPAR
c
      call getival('nblock1' ,nblock1)
c333  call blk_pairs_cal3(nbl1,bl(nblock1),nbl2,maxpar,ncs,schwarz)
c222  call blk_pairs_cal2(nbl1,bl(nblock1),nbl2,maxpar,ncs,schwarz)
      call blk_pairs_cal1(nbl1,bl(nblock1),nbl2,maxpar,ncs,schwarz)
c
c output : nbl2
c output : maxpar    - maximum number of pairs per block2
c---------------------------------------------------------------
c number of quartet-blocks (super-blocks) :
c
      nbl4=nbl2*(nbl2+1)/2
c
c---------------------------------------------------------------
c Allocate memory for NPAR and IJBL and for MXSIZE
c
      call getint(nbl2       ,npard)
      call getint(nbl2*maxpar,ijbld)
      call getint(nbl4       ,mxsize)
c
c output : addresses : npard,ijbld & mxsize
c---------------------------------------------------------------
c Constructe blocks of contracted pairs:Setup IJBL(nbl2,*)and  NPAR(*)
c
c
c333  call blk_pairs_mak3(nbl1,bl(nblock1),nbl2,bl(npard),bl(ijbld),
c333 *                    bl(map_ij_bl2) ,ncs,schwarz)
c222  call blk_pairs_mak2(nbl1,bl(nblock1),nbl2,bl(npard),bl(ijbld),
c222 *                    bl(map_ij_bl2) ,ncs,schwarz)
      call blk_pairs_mak1(nbl1,bl(nblock1),nbl2,bl(npard),bl(ijbld),
     *                    bl(map_ij_bl2) ,ncs,schwarz)
c---------------------------------------------------------------
c reorder pairs in each block according to the schwarz integ. values :
c
      call blk_pairs_reor(ncs,schwarz,nbl2,nbl2,bl(npard),bl(ijbld))
c---------------------------------------------------------------
c  Set up the MXSIZE and NSUPB vectors
c mxsize(isupbl)=> maiximum size of a small block in isupbl big-block
c nsupb (isupbl)=> number of small blocks in isupbl big-block
c
      call getint(nbl4,nsupb)
c
      call blksizer(bl,bl(nsupb),nbl2,bl(ijbld),bl(npard),
     *              iis,jjs,inx, bl(mxsize),lcore,nsym,ncs,iroute)
c
c---------------------------------------------------------------
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc0)
      If(iftc0.NE.0) Then
       call getival('ccdd',iccdd)
       call getint(nbl4,nftcbl4)
cc     call print2cd(ncs,schwarz,nbl2,bl(npard),bl(ijbld),bl(iftc0))
       call ftc_bl4(ncs,nbl2,bl(ijbld),bl(iftc0),iccdd,
     *              bl(nftcbl4),nbl4ftc)
       call retmem(1)
       call getint(nbl4ftc,nftcbl4)
       call setival('nftcbl4',nftcbl4)   ! pointer
       call setival('nbl4ftc',nbl4ftc)   ! no of FTC blocks
      endif
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairs_cal3(nbl1,nblock1,nbl2,maxpar,ncs,schwarz)
c--------------------------------------------------------------
c split pair blocks into THREE groups with :
c
c (1) big    schwarz int. :              (ij|ij) >= sqrt(eps)
c                                                   10**-5
c (2) midium schwarz int. :  sqrt(eps) > (ij|ij) >= eps*sqrt(eps)
c                              10*-5                 10**-10
c (3) small  schwarz int. :  eps*sqrt()> (ij|ij) > eps**2
c                              10**-10             10**-20
c note :
c  all pairs with (ij|ij) < 10**-20 are removed
c--------------------------------------------------------------
c
c
      implicit real*8 (a-h,o-z)
      logical present
      common /neglect/ eps,eps1,epsr,eps8
      common /outfile/ ioutput
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension schwarz(ncs,ncs)
c--------------------------------------------------------------
c input : nbl1 - number of single shell blocks
c
c output: nbl2 - number of shell-pairs  blocks
c output: maxpar - max. number of pairs/block
c--------------------------------------------------------------
      eps20=eps*eps      ! 10**-20 by default
      eps10=eps          ! 10**-10 by default
      eps5=sqrt(eps)     ! 10**-5  by default
      eps15=eps*eps5     ! 10**-15 by default
c
      ncsp_tot=0     ! total number of pairs
      ncsp_ret=0     ! number of retained pairs
      maxpar=0
      ijbl=0                    !    block's counter (with limit)
c--------------------------------------------------------------
c First make pair blocks containing pairs with big (ij|ij) > 10**-5
c--------------------------------------------------------------
c make pair blocks with big schwraz integ. : (ij|ij) >= eps5
      do 101 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 201 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_big=0
            present=.false.
            do 301 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 401 jcs=jbeg,jenx
                  ncsp_tot=ncsp_tot+1
                  x_ij_ij=schwarz(ics,jcs)
c................ if(x_ij_ij.GE.eps20) then
                     if(x_ij_ij.GE.eps5) then
                        present=.true.
                        ncsp_ret=ncsp_ret+1
                        if(ijpar_big.ge.limpair) then
                           ijbl=ijbl+1
                           ijpar_big=0
                        endif
                        ijpar_big=ijpar_big+1
                        maxpar=max(maxpar,ijpar_big)
                     endif
c................ endif
  401          continue
  301       continue
            if(.not.present) ijbl=ijbl-1
  201    continue
  101 continue
c
c--------------------------------------------------------------
c remember number of pair blocks with "big" Schwarz integrals :
c
      ijbl_big=ijbl
c--------------------------------------------------------------
c make pair blocks with midium schwraz integ.: 10-5 > (ij|ij) >= 10-10
      do 102 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 202 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_mid=0
            present=.false.
            do 302 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 402 jcs=jbeg,jenx
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps10) then
                     if(x_ij_ij.LT.eps5) then
                        present=.true.
                        ncsp_ret=ncsp_ret+1
                        if(ijpar_mid.ge.limpair) then
                           ijbl=ijbl+1
                           ijpar_mid=0
                        endif
                        ijpar_mid=ijpar_mid+1
                        maxpar=max(maxpar,ijpar_mid)
                     endif
                  endif
  402          continue
  302       continue
            if(.not.present) ijbl=ijbl-1
  202    continue
  102 continue
c--------------------------------------------------------------
c remember number of pair blocks with "midium" Schwarz integrals :
c
      ijbl_mid=ijbl-ijbl_big
c--------------------------------------------------------------
c make pair blocks with small schwraz integ. : eps10> (ij|ij) >= eps20
      do 103 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 203 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_sml=0
            present=.false.
            do 303 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 403 jcs=jbeg,jenx
                  x_ij_ij=schwarz(ics,jcs)
ccccc?????        if(x_ij_ij.GE.eps20) then
                  if(x_ij_ij.GE.eps15) then
                     if(x_ij_ij.LT.eps10) then
                        present=.true.
                        ncsp_ret=ncsp_ret+1
                        if(ijpar_sml.ge.limpair) then
                           ijbl=ijbl+1
                           ijpar_sml=0
                        endif
                        ijpar_sml=ijpar_sml+1
                        maxpar=max(maxpar,ijpar_sml)
                     endif
                  endif
  403          continue
  303       continue
            if(.not.present) ijbl=ijbl-1
  203    continue
  103 continue
c--------------------------------------------------------------
      nbl2=ijbl   ! number of pair-blocks with the limit limpair
c---------------------------------------------------------------
c save number of pairs which left after schwarz screening :
c
      call setival('ncspairs',ncsp_ret)
c---------------------------------------------------------------
      ijbl_sml=ijbl-ijbl_big-ijbl_mid
      write(ioutput,*  ) '                 '
      write(ioutput,*  ) 'Second Blocking Procedure '
      write(ioutput,495) nbl2,ijbl_big,ijbl_mid,ijbl_sml,
     *                   ncsp_tot,ncsp_ret,maxpar
  495 format( ' Number of Blocks of Contracted Shell Pairs    =',i7/
     *        '    blocks with big    schwarz integrals       =',i7/
     *        '    blocks with midium schwarz integrals       =',i7/
     *        '    blocks with small  schwarz integrals       =',i7/
     *        ' total number of contracted shell pairs        =',i7/
     *        ' number of contracted shell pairs retained     =',i7/
     *        ' maximum number of shell-pairs/block           =',i7)
c---------------------------------------------------------------
c
      end
c===============================================================
      subroutine blk_pairs_mak3(nbl1,nblock1,nbl2,npar,ijbl,map_ij_bl2,
     *                          ncs,schwarz)
      implicit real*8 (a-h,o-z)
      logical present
      common /outfile/ ioutput
      common /intlim/ limxmem,limblks,limpair
      common /neglect/ eps,eps1,epsr,eps8
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncsp)
      dimension schwarz(ncs,ncs)
c---------------------------------------------------------------
c nbl1 - number of single shell blocks
c nbl2 - number of shell-pairs  blocks
c
c constructe pairs of contracted shells
c---------------------------------------------------------------
      eps20=eps*eps      ! 10**-20
      eps10=eps          ! 10**-10
      eps5=sqrt(eps)     ! 10**-5
      eps15=eps*eps5     ! 10**-15 by default
c---------------------------------------------------------------
      ijblock=0             !    blocks counter
c---------------------------------------------------------------
c make pair blocks with big schwraz integ. : (ij|ij) >= eps5
      do 101 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 201 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_big=0
            present=.false.
            do 301 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 401 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  x_ij_ij=schwarz(ics,jcs)
c................ if(x_ij_ij.GE.eps20) then
cok                  if(x_ij_ij.GE.eps) then
                     if(x_ij_ij.GE.eps5) then
                        present=.true.
                        if(ijpar_big.ge.limpair) then
                           ijblock=ijblock+1
                           ijpar_big=0
                        endif
                        ijpar_big=ijpar_big+1
                        npar(ijblock)=ijpar_big
                        ijbl(ijblock,ijpar_big)=ijcs
                        map_ij_bl2(ijcs)=ijblock
                     endif
c................ endif
  401          continue
  301       continue
            if(.not.present) ijblock=ijblock-1
  201    continue
  101 continue
c--------------------------------------------------------------
c remember the last pair block with "big" Schwarz integrals :
c
      last_big2=ijblock
      call setival('last_big2',ijblock)
c---------------------------------------------------------------
c make pair blocks with midium schwraz integ.: 10-5> (ij|ij) >=10-10
      do 102 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 202 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_mid=0
            present=.false.
            do 302 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 402 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps10) then
                     if(x_ij_ij.LT.eps5) then
                        present=.true.
                        if(ijpar_mid.ge.limpair) then
                           ijblock=ijblock+1
                           ijpar_mid=0
                        endif
                        ijpar_mid=ijpar_mid+1
                        npar(ijblock)=ijpar_mid
                        ijbl(ijblock,ijpar_mid)=ijcs
                        map_ij_bl2(ijcs)=ijblock
                     endif
                  endif
  402          continue
  302       continue
            if(.not.present) ijblock=ijblock-1
  202    continue
  102 continue
c--------------------------------------------------------------
c remember the last pair block with "midium" Schwarz integrals :
c
      last_mid2=ijblock
      call setival('last_mid2',ijblock)
c---------------------------------------------------------------
c make pair blocks with small schwraz integ. : eps10> (ij|ij) >= eps20
      do 103 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 203 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_sml=0
            present=.false.
            do 303 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 403 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  x_ij_ij=schwarz(ics,jcs)
cccc??????        if(x_ij_ij.GE.eps20) then
                  if(x_ij_ij.GE.eps15) then
                     if(x_ij_ij.LT.eps10) then
                        present=.true.
                        if(ijpar_sml.ge.limpair) then
                           ijblock=ijblock+1
                           ijpar_sml=0
                        endif
                        ijpar_sml=ijpar_sml+1
                        npar(ijblock)=ijpar_sml
                        ijbl(ijblock,ijpar_sml)=ijcs
                        map_ij_bl2(ijcs)=ijblock
                     endif
                  endif
  403          continue
  303       continue
            if(.not.present) ijblock=ijblock-1
  203    continue
  103 continue
c---------------------------------------------------------------
c calculate number of pairs with small schwarz integrals
c
      ijpar_big=0
      do ibl=1,last_big2
         ijpar=npar(ibl)
         ijpar_big= ijpar_big + ijpar
      enddo
      ijpar_mid=0
      do ibl=last_big2+1,last_mid2
         ijpar=npar(ibl)
         ijpar_mid= ijpar_mid + ijpar
      enddo
      ijpar_sml=0
      do ibl=last_mid2+1,nbl2
         ijpar=npar(ibl)
         ijpar_sml= ijpar_sml + ijpar
      enddo
c
      ijpar_tot=ijpar_big+ijpar_mid+ijpar_sml
c
      write(ioutput,495) ijpar_big,ijpar_mid,ijpar_sml,ijpar_tot
  495 format( ' number of pairs with big    schwarz integrals =',i7/
     *        ' number of pairs with midium schwarz integrals =',i7/
     *        ' number of pairs with small  schwarz integrals =',i7/
     *        ' Total number of contracted shells pairs       =',i7)
c---------------------------------------------------------------
c     write(6,*)' nbl2=',nbl2,' ijblock=',ijblock
c     ijpar_sum=0
c     do 2500 ibl=1,nbl2
c     ijpar=npar(ibl)
c     write(ioutput,502) ibl,ijpar
c     write(ioutput,503) (ijbl(ibl,ij),ij=1,ijpar)
c     ijpar_sum=ijpar_sum+ijpar
c2500 continue
c 502 format('Pair-Block =',i4,' : no of shell-pairs =',i4)
c 503 format('  pairs:',14(i4,1x) )
c     write(ioutput,*)'total number of pairs=',ijpar_sum
c---------------------------------------------------------------
      end
c=======================================================================
      subroutine get_max_am(itype1,jtype1,ktype1,ltype1,mmax,
     *                      nsij,nskl,nqmax)
      common /logic1/ ndege(1)
c
c returns total angular momentum for (ij,kl) quartet
c-
         nqi=ndege(itype1)
         nqj=ndege(jtype1)
         nsij=nqi+nqj-1
c
         nqk=ndege(ktype1)
         nql=ndege(ltype1)
         nskl=nqk+nql-1
c
         mmax=nsij+nskl-1
c
         nqmax= max(nqi,nqj,nqk,nql)
c
      end
c=======================================================================
      subroutine get_limit(mmax,icache,nfijkl,ifor,ibl,kbl,
     *                     itype1,jtype1,ktype1,ltype1,
     *                     nqrt_limit)
      implicit real*8 (a-h,o-z)
      common /intlim/ limxmem,limblks,limpair
      dimension limitqp(14)    ! limits for a given tot.ang.mom.
      dimension limitql(14)    ! limits for a given tot.ang.mom.
c maximum number of functions for a given total angular momentum
      dimension limitfp(14)    ! if there is no l-shells
      dimension limitfl(14)    ! if there are l-shells
c----------------------------------------------------------------------
c Input :
c
c   mmax  - total angular momentum (+1) for (ij|kl)
c  icache - current cache memory size
c  nfijkl - number of functions (integrals) in (ij|kl)
c    ifor - type of task (scf,giao,force,hess=1,2,3,4)
c ibl,kbl - pair blocks
c    itype1,jtype1,ktype1,ltype1 - type of shells
c
c Output :
c nqrt_limit - block size limit
c----------------------------------------------------------------------
c     ang.mom.+1    1   2   3   4   5   6   7   8   9  10 11 12 13 14
c optimized for p-shell:
ccc   data limitqp /1000,500,250, 75, 50, 50, 25, 25, 12, 7, 4, 3, 2, 1/
      data limitqp /1000,500,250,100, 75, 50, 40, 30, 12, 7, 4, 3, 2, 1/
c                 ssss psss ppss ppps pppp dppp ddpp
c----------------------------------------------------------------------
      data limitql /1000,350,200, 75, 55, 40, 30, 25, 12, 7, 4, 3, 2, 1/
c                 ssss lsss llss llls llll dlll ddll
c----------------------------------------------------------------------
c number of functions :
c                   ssss psss ppss ppps pppp dppp ddpp
      data limitfp /   1,   3,   9,  27,  81, 162, 324,
     *               648,1296,2160,3600,6000,10000,15000/
c                   dddp dddd fddd ffdd fffd ffff  gfff
c----------------------------------------------------------------------
c                   ssss lsss llss llls llll dlll ddll
      data limitfl /   1,   4,  16,  64, 256, 384, 576,
     *               864,1296,2160,3600,6000,10000,15000/
c                   dddl dddd fddd ffdd fffd ffff  gfff
c------------------------------------------------------------------
      data min_cache /16384/
c------------------------------------------------------------------
c from tests on aspirin/6-31g (l shells segmented into s,p)
c 1. limits for (ss|ss) do not change scf or force integral timings
c------------------------------------------------------------------
c 1) for ordinary two-el.integrals (ifor=1)
c 2) for GIAO two-el. derivatives  (ifor=2)
c 3) for gradient derivatives      (ifor=3)
c 4) for second derivatives        (ifor=4)
c 5) for mp2/lmp2 integrals        (ifor=5)
c------------------------------------------------------------------
c        N E W    B L O C K    S I Z E    L I M I T S
c------------------------------------------------------------------
c  maximum number of quartets as a function of integral's type
c     fitted to the minimum cache size = 16K = 16384 bytes
c                    min_cache= 16384
c        limit in        integral  angular    no of    limits from
c      no of quartes       type    mom. +1  integrals  icache/nfunc
c------------------------------------------------------------------
c     limitq( 1)=1000   !  ssss      1           1       16384
c     limitq( 2)= 500   !  psss      2           3        5461
c
c     limitq( 3)= 250   !  ppss      3           9        1820
c                          dsss                  6        2731
c
c     limitq( 4)=  75   !  ppps      4          27         606
c                          dpss                 18         910
c                          fsss                 10        1638
c
c     limitq( 5)=  50   !  pppp      5          81         202
c                          dpps                 54         303
c                          ddss                 36         455
c                          fpss                 30         546
c                          gsss                 15        1092
c
c     limitq( 6)=  50   !  dppp      6         162         101
c                          ddps                108         152
c                          fdss                 60         273
c                          gpss                 45         364
c                          hsss                 21         780
c
c     limitq( 7)=  25   !  ddpp      7         324          50
c                          ddds                216          76
c                          fdps                180          91
c                          ffss                100         164
c                          gdss                 90         182
c                          hpss                 63         260
c                          isss                 28         585
c
c     limitq( 8)= 25    !  dddp      8         648          25
c                          fdds                360          46
c                          ffps                300          55
c                          gfss                150         110
c                          hdss                126         130
c                          ipss                 84         195
c
c     limitq( 9)= 12    !  dddd      9        1296          12
c                          fddp               1080          15
c                          ffds                600          27
c                          gfps                450          36
c                          ggss                225          73
c                          hfss                210          78
c                          idss                168          98
c                          jpss                108         152
c                          ksss                 45         364
c
c     limitq(10)=  7    !  fddd     10        2160           7
c                          ffdp               1800           9
c                          fffs               1000          16
c                          gfds                900          18
c                          ggps                675          24
c                          hgss                315          52
c
c     limitq(11)=  4    !  ffdd     11        3600           4
c                          fffp               3000           5
c                          gffs               1500          11
c                          hfds               1260          13
c                          hgps                945          17
c                          hhss                441          37
c                          igss                420          39
c                          jfss                280          59
c                          kdss                270          61
c                          lpss                165          99
c                          msss                 66         248
c
c     limitq(12)=  3    !  fffd     12        6000           3
c     limitq(13)=  2    !  ffff     13       10000           2
c     limitq(14)=  1    !  gfff     14       15000           1
c
c  all blocks with mmax over 13 have block szie = 1
c
c     limitq(15)=  1    !  ggff     15       22500           1
c     limitq(16)=  1    !  gggf     16       33750           1
c     limitq(17)=  1    !  gggg     17       50625           1
c     limitq(18)=  1    !  hggg     18       70875           1
c     limitq(19)=  1    !  hhgg     19       99225           1
c     limitq(20)=  1    !  hhhg     20      138915           1
c     limitq(21)=  1    !  hhhh     21      194481           1
c------------------------------------------------------------------
      if(limblks.NE.0) then
         nqrt_limit=limblks
         return
      endif
c------------------------------------------------------------------
      itask=ifor
c
      if(mmax.ge.14) then
         nqrt_limit=1
         return
      endif
c------------------------------------------------------------------
c check for p- and l-shells
c
      ip=0
      if(itype1.eq.2.or.jtype1.eq.2.or.ktype1.eq.2.or.ltype1.eq.2) ip=1
      il=0
      if(itype1.eq.3.or.jtype1.eq.3.or.ktype1.eq.3.or.ltype1.eq.3) il=1
c
      if(ip.eq.0 .and. il.eq.0) then
c        no l-, no p-shells
         nqrt_limit=limitqp(mmax)
         nfunct=    limitfp(mmax)
      endif
      if(ip.eq.1 .and. il.eq.0) then
         nqrt_limit=limitqp(mmax)
         nfunct=    limitfp(mmax)
      endif
      if(ip.eq.0 .and. il.eq.1) then
         nqrt_limit=limitql(mmax)
         nfunct=    limitfl(mmax)
      endif
      if(ip.eq.1 .and. il.eq.1) then
         nqrt_limit=( limitqp(mmax)+limitql(mmax) )/2
         nfunct=    ( limitfp(mmax)+limitfl(mmax) )/2
      endif
c------------------------------------------------------------------
c cache factor
c
      cachef=dble(icache)/dble(min_cache)
      if(cachef.lt.1.d0) cachef=1.d0
      xqrt_limit=cachef*dble(nqrt_limit)
c------------------------------------------------------------------
c number of functions factor : e.g. bigger blocks for fsss than for ppps
c
      functf=dble(nfunct)/dble(nfijkl)
      if(functf.lt.1.d0) functf=1.d0
      xqrt_limit=functf*xqrt_limit
c------------------------------------------------------------------
c task factor : make smaller blocks for more demending tasks
c
      taskf=1.d0
cgiao if(itask.eq.2) taskf=0.85d0
c     if(itask.eq.3) taskf=0.2500d0
c     if(itask.eq.4) taskf=0.3333d0
cgiao if(itask.eq.2) taskf=0.85d0
      if(itask.eq.3) taskf=0.2500d0
      if(itask.eq.4) taskf=0.1000d0
c.....................................
      xqrt_limit=taskf*xqrt_limit
c------------------------------------------------------------------
      nqrt_limit=int(xqrt_limit)
      if(nqrt_limit.lt.1) nqrt_limit=1
c
      end
c=======================================================================
      subroutine blk_pairs_cal2(nbl1,nblock1,nbl2,maxpar,ncs,schwarz)
c--------------------------------------------------------------
c split pair blocks into TWO  groups with :
c
c (1) big    schwarz int. :              (ij|ij) >= eps
c                                                   10**-10
c (2) small  schwarz int. :        eps > (ij|ij) >= eps**2
c                                10**-10            10**-20
c note :
c  all pairs with (ij|ij) < 10**-20 are removed
c--------------------------------------------------------------
c
c
      implicit real*8 (a-h,o-z)
      logical present
      common /neglect/ eps,eps1,epsr,eps8
      common /outfile/ ioutput
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension schwarz(ncs,ncs)
c--------------------------------------------------------------
c input : nbl1 - number of single shell blocks
c
c output: nbl2 - number of shell-pairs  blocks
c output: maxpar - max. number of pairs/block
c--------------------------------------------------------------
      eps20=eps*eps       ! 10**-20 by default
      eps15=eps*sqrt(eps) ! 10**-20 by default
      eps10=eps           ! 10**-10 by default
c
      ncsp_tot=0     ! total number of pairs
      ncsp_ret=0     ! number of retained pairs
      maxpar=0
      ijbl=0                    !    block's counter (with limit)
c--------------------------------------------------------------
c make pair blocks with big schwraz integ. : (ij|ij) >= eps10
c
      do 101 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 201 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_big=0
            present=.false.
            do 301 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 401 jcs=jbeg,jenx
                  ncsp_tot=ncsp_tot+1
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps10) then
c15               if(x_ij_ij.GE.eps15) then
                     present=.true.
                     ncsp_ret=ncsp_ret+1
                     if(ijpar_big.ge.limpair) then
                        ijbl=ijbl+1
                        ijpar_big=0
                     endif
                     ijpar_big=ijpar_big+1
                     maxpar=max(maxpar,ijpar_big)
                  endif
  401          continue
  301       continue
            if(.not.present) ijbl=ijbl-1
  201    continue
  101 continue
c
c--------------------------------------------------------------
c remember number of pair blocks with "big" Schwarz integrals :
c
      ijbl_big=ijbl
c--------------------------------------------------------------
c make pair blocks with small schwraz integ. : eps10> (ij|ij) >= eps20
      do 103 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 203 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_sml=0
            present=.false.
            do 303 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 403 jcs=jbeg,jenx
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps20) then
                     if(x_ij_ij.LT.eps10) then
c15                  if(x_ij_ij.LT.eps15) then
                        present=.true.
                        ncsp_ret=ncsp_ret+1
                        if(ijpar_sml.ge.limpair) then
                           ijbl=ijbl+1
                           ijpar_sml=0
                        endif
                        ijpar_sml=ijpar_sml+1
                        maxpar=max(maxpar,ijpar_sml)
                     endif
                  endif
  403          continue
  303       continue
            if(.not.present) ijbl=ijbl-1
  203    continue
  103 continue
c--------------------------------------------------------------
      nbl2=ijbl   ! number of pair-blocks with the limit limpair
c---------------------------------------------------------------
c save number of pairs which left after schwarz screening :
c
      call setival('ncspairs',ncsp_ret)
c---------------------------------------------------------------
      ijbl_sml=ijbl-ijbl_big
      write(ioutput,*  ) '                 '
      write(ioutput,*  ) 'Second Blocking Procedure '
      write(ioutput,495) nbl2,ijbl_big,ijbl_sml,
     *                   ncsp_tot,ncsp_ret,maxpar
  495 format( ' Number of Blocks of Contracted Shell Pairs    =',i7/
     *        '    blocks with big    schwarz integrals       =',i7/
     *        '    blocks with small  schwarz integrals       =',i7/
     *        ' total number of contracted shell pairs        =',i7/
     *        ' number of contracted shell pairs retained     =',i7/
     *        ' maximum number of shell-pairs/block           =',i7)
c---------------------------------------------------------------
c
      end
c===============================================================
      subroutine blk_pairs_mak2(nbl1,nblock1,nbl2,npar,ijbl,map_ij_bl2,
     *                          ncs,schwarz)
      implicit real*8 (a-h,o-z)
      logical present
      common /outfile/ ioutput
      common /intlim/ limxmem,limblks,limpair
      common /neglect/ eps,eps1,epsr,eps8
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncsp)
      dimension schwarz(ncs,ncs)
c---------------------------------------------------------------
c nbl1 - number of single shell blocks
c nbl2 - number of shell-pairs  blocks
c
c constructe pairs of contracted shells
c---------------------------------------------------------------
      eps20=eps*eps       ! 10**-20 by default
      eps15=eps*sqrt(eps) ! 10**-20 by default
      eps10=eps           ! 10**-10 by default
c---------------------------------------------------------------
      ijblock=0             !    blocks counter
c---------------------------------------------------------------
c make pair blocks with big schwraz integ. : (ij|ij) >= eps10
c
      do 101 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 201 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_big=0
            present=.false.
            do 301 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 401 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps10) then
c15               if(x_ij_ij.GE.eps15) then
                     present=.true.
                     if(ijpar_big.ge.limpair) then
                        ijblock=ijblock+1
                        ijpar_big=0
                     endif
                     ijpar_big=ijpar_big+1
                     npar(ijblock)=ijpar_big
                     ijbl(ijblock,ijpar_big)=ijcs
                     map_ij_bl2(ijcs)=ijblock
                  endif
  401          continue
  301       continue
            if(.not.present) ijblock=ijblock-1
  201    continue
  101 continue
c--------------------------------------------------------------
c remember the last pair block with "big" Schwarz integrals :
c
      last_big2=ijblock
      call setival('last_big2',ijblock)
      call setival('last_mid2',ijblock)  ! needed in get_limits
c---------------------------------------------------------------
c make pair blocks with small schwraz integ. : eps10> (ij|ij) >= eps20
c
      do 103 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 203 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_sml=0
            present=.false.
            do 303 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 403 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps20) then
                     if(x_ij_ij.LT.eps10) then
c15                  if(x_ij_ij.LT.eps15) then
                        present=.true.
                        if(ijpar_sml.ge.limpair) then
                           ijblock=ijblock+1
                           ijpar_sml=0
                        endif
                        ijpar_sml=ijpar_sml+1
                        npar(ijblock)=ijpar_sml
                        ijbl(ijblock,ijpar_sml)=ijcs
                        map_ij_bl2(ijcs)=ijblock
                     endif
                  endif
  403          continue
  303       continue
            if(.not.present) ijblock=ijblock-1
  203    continue
  103 continue
c---------------------------------------------------------------
c calculate number of pairs with big and small schwarz integrals
c
      ijpar_big=0
      do ibl=1,last_big2
         ijpar=npar(ibl)
         ijpar_big= ijpar_big + ijpar
      enddo
      ijpar_sml=0
      do ibl=last_big2+1,nbl2
         ijpar=npar(ibl)
         ijpar_sml= ijpar_sml + ijpar
      enddo
c
      ijpar_tot=ijpar_big+ijpar_sml
c
      write(ioutput,495) ijpar_big,ijpar_sml,ijpar_tot
  495 format( ' number of pairs with big    schwarz integrals =',i7/
     *        ' number of pairs with small  schwarz integrals =',i7/
     *        ' Total number of contracted shells pairs       =',i7)
c---------------------------------------------------------------
c     write(6,*)' nbl2=',nbl2,' ijblock=',ijblock
c     ijpar_sum=0
c     do 2500 ibl=1,nbl2
c     ijpar=npar(ibl)
c     write(ioutput,502) ibl,ijpar
c     write(ioutput,503) (ijbl(ibl,ij),ij=1,ijpar)
c     ijpar_sum=ijpar_sum+ijpar
c2500 continue
c 502 format('Pair-Block =',i4,' : no of shell-pairs =',i4)
c 503 format('  pairs:',14(i4,1x) )
c     write(ioutput,*)'total number of pairs=',ijpar_sum
c---------------------------------------------------------------
      end
c=======================================================================
      subroutine blk_pairs_cal1(nbl1,nblock1,nbl2,maxpar,ncs,schwarz)
c--------------------------------------------------------------
c Do NOT split pair blocks into groups with big, medium or small
c schwarz integrals. Just remove all pairs with
c                 (ij|ij) < 10**-20
c--------------------------------------------------------------
c
c
      implicit real*8 (a-h,o-z)
      logical present
      common /neglect/ eps,eps1,epsr,eps8
      common /outfile/ ioutput
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension schwarz(ncs,ncs)
c--------------------------------------------------------------
c input : nbl1 - number of single shell blocks
c
c output: nbl2 - number of shell-pairs  blocks
c output: maxpar - max. number of pairs/block
c--------------------------------------------------------------
      eps20=eps*eps        ! 10**-20 by default
ccc   eps15=eps*sqrt(eps)  ! 10**-15 by default
c
      ncsp_tot=0     ! total number of pairs
      ncsp_ret=0     ! number of retained pairs
      maxpar=0
      ijbl=0                    !    block's counter (with limit)
c--------------------------------------------------------------
c make pair blocks ; include only pairs wih : (ij|ij) >= eps20
c
      do 101 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 201 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_big=0
            present=.false.
            do 301 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 401 jcs=jbeg,jenx
                  ncsp_tot=ncsp_tot+1
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps20) then
                     present=.true.
                     ncsp_ret=ncsp_ret+1
                     if(ijpar_big.ge.limpair) then
                        ijbl=ijbl+1
                        ijpar_big=0
                     endif
                     ijpar_big=ijpar_big+1
                     maxpar=max(maxpar,ijpar_big)
                  endif
  401          continue
  301       continue
            if(.not.present) ijbl=ijbl-1
  201    continue
  101 continue
c--------------------------------------------------------------
      nbl2=ijbl   ! number of pair-blocks with the limit limpair
c---------------------------------------------------------------
c save number of pairs which left after schwarz screening :
c
      call setival('ncspairs',ncsp_ret)
c---------------------------------------------------------------
      write(ioutput,*  ) '                 '
      write(ioutput,*  ) 'Second Blocking Procedure '
      write(ioutput,495) nbl2,
     *                   ncsp_tot,ncsp_ret,eps,maxpar
  495 format( ' Number of Blocks of Contracted Shell Pairs    =',i7/
     *        ' total number of contracted shell pairs        =',i7/
     *        ' number of contracted shell pairs retained     =',i7/
     *        ' (those with (ij|ij) >= eps**2, where eps =',1pe11.4,')'/
     *        ' maximum number of shell-pairs/block           =',i7)
c---------------------------------------------------------------
c
      end
c===============================================================
      subroutine blk_pairs_mak1(nbl1,nblock1,nbl2,npar,ijbl,map_ij_bl2,
     *                          ncs,schwarz)
      implicit real*8 (a-h,o-z)
      logical present
      common /outfile/ ioutput
      common /intlim/ limxmem,limblks,limpair
      common /neglect/ eps,eps1,epsr,eps8
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncsp)
      dimension schwarz(ncs,ncs)
c---------------------------------------------------------------
c nbl1 - number of single shell blocks
c nbl2 - number of shell-pairs  blocks
c
c constructe pairs of contracted shells
c---------------------------------------------------------------
      eps20=eps*eps        ! 10**-20
c---------------------------------------------------------------
      ijblock=0             !    blocks counter
c---------------------------------------------------------------
c make pair blocks with schwraz integ. : (ij|ij) >= eps20
c
      do 101 ibl=1,nbl1
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 201 jbl=1,ibl
            ijblock=ijblock+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
c
            ijpar_big=0
            present=.false.
            do 301 ics=ibeg,iend
               iics=ics*(ics-1)/2
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 401 jcs=jbeg,jenx
                  ijcs=iics+jcs
                  x_ij_ij=schwarz(ics,jcs)
                  if(x_ij_ij.GE.eps20) then
                     present=.true.
                     if(ijpar_big.ge.limpair) then
                        ijblock=ijblock+1
                        ijpar_big=0
                     endif
                     ijpar_big=ijpar_big+1
                     npar(ijblock)=ijpar_big
                     ijbl(ijblock,ijpar_big)=ijcs
                     map_ij_bl2(ijcs)=ijblock
                  endif
  401          continue
  301       continue
            if(.not.present) ijblock=ijblock-1
  201    continue
  101 continue
c--------------------------------------------------------------
c remember the last pair block with "big" Schwarz integrals :
c
      last_big2=ijblock
      call setival('last_big2',ijblock)
      call setival('last_mid2',ijblock)  ! needed in get_limits
c---------------------------------------------------------------
c     write(6,*)' nbl2=',nbl2,' ijblock=',ijblock
c     ijpar_sum=0
c     do 2500 ibl=1,nbl2
c     ijpar=npar(ibl)
c     write(ioutput,502) ibl,ijpar
c     write(ioutput,503) (ijbl(ibl,ij),ij=1,ijpar)
c     ijpar_sum=ijpar_sum+ijpar
c2500 continue
c 502 format('Pair-Block =',i4,' : no of shell-pairs =',i4)
c 503 format('  pairs:',14(i4,1x) )
c     write(ioutput,*)'total number of pairs=',ijpar_sum
c---------------------------------------------------------------
      end
c=======================================================================
      subroutine get_cond0(exi0,exi0l, exi1,exi1l,ist0,cond0)
      implicit real*8 (a-h,o-z)
      logical cond0
c
      cond0=.false.
      if(ist0.eq.1) then                  !    s-shell
         emax=100000.000d0
         emin=     0.005d0
      endif
      if(ist0.eq.2 .or. ist0.eq.3) then   !    p,l-shell
         emax= 10000.000d0
         emin=     0.050d0
      endif
      if(ist0.eq.4 .or. ist0.eq.5) then   !    d-shell
         emax=  1000.000d0
         emin=     0.500d0
      endif
      if(ist0.eq.6 .or. ist0.eq.7) then   !    f-shell
         emax=   100.000d0
         emin=     1.000d0
      endif
      if(ist0.gt.7) then                  !    g-shell and higher
         emax=    10.000d0
         emin=     1.000d0
      endif
c
      cond0=.true.
      if(exi0.gt.emax .or. exi1.gt.emax) cond0=.false.
      if(exi0l.lt.emin .or. exi1l.lt.emin) cond0=.false.
c
      end
c=======================================================================
      subroutine check_exp(ibl,inx,ijbl,nbl2,ijpar,iis,jjs,datbas,
     *                     boamax,rapbmax)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
      dimension datbas(13,*)
c
      exp_max_j=0.d0
      exp_min_i=1000000.d0
      exp_min_j=1000000.d0
      do ijpax=1,ijpar
         ijcs=ijbl(ibl,ijpax)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         ips_b=inx(1,ics)+1
         ips_e=inx(5,ics)
         jps_b=inx(1,jcs)+1
         jps_e=inx(5,jcs)
         do ips=ips_b,ips_e
            exp_i=datbas(1,ips)
            exp_min_i=min(exp_min_i,exp_i)
         enddo
         do jps=jps_b,jps_e
            exp_j=datbas(1,jps)
            exp_max_j=max(exp_max_j,exp_j)
            exp_min_j=min(exp_min_j,exp_j)
         enddo
      enddo
c
      apbmin=exp_min_i + exp_min_j
      boamax=exp_max_j/max(1.d0,exp_min_i)
      boamax= 2.d0*boamax
      rapbmax=1.d0/apbmin
c
      end
c=======================================================================
      subroutine num_stab(mmax,nsij,nskl,nqmax,
     *                    boamax,rcpdmax,docmax,rapbmax,stable)
      implicit real*8 (a-h,o-z)
      logical stable
      common /neglect/ eps,eps1,epsr
      lost_allow=15+log10(eps)   ! 15 from double precision accuracy
      if(lost_allow.lt.0) lost_allow=0
c
      stable=.true.
c
      nshifts=min(nsij,nskl)-1
      if( mmax.le.6 ) return                 ! dp|dp = 7
      if( nqmax.le.2) return  ! up to p only
      if(nshifts.le.2) return
c-----------------------------------------------------------------------
      if(nsij.ge.nskl) then
c        shifting in tracy's req. from 1 to 3
         call digit_lost(boamax,rcpdmax,nshifts,lost_digit)
         if(lost_digit.gt.lost_allow) stable=.false.
      endif
c-----------------------------------------------------------------------
      if(nsij.lt.nskl) then
c        shifting in tracy's req. from 3 to 1
         call digit_lost(docmax,rapbmax,nshifts,lost_digit)
         if(lost_digit.gt.lost_allow) stable=.false.
      endif
c-----------------------------------------------------------------------
      end
c====================================================================
      subroutine redefmem(memtot,needed,limxmem)
      common /outfile/ ioutput
c
c     lcore - total memory ( with offset in)
c     mqrt1 - memory needed for ONE (highest) quartet
c     limxmem - limit memory (too small)
c
c     redefine limxmem
      if(needed.ge.memtot) then
         call nerror(912,'redefmem',
     *              'lcore too samll for new limxmem:',
     *               memtot,needed)
      else
         write(ioutput,*)
     *   '  limxmem redefined:  old=',limxmem,'  new=',needed
         call f_lush(6)
         limxmem=needed
      endif
c
      end
c====================================================================
      subroutine blkmemor(bl,nbl2,ijbl,npar,iis,jjs,inx,
     *                    lcore,iroute)
c-----------------------------------------------------------
c This routine checks out if LIMXMEM (memory limit for integrals)
c should be increased
c-----------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      common /intlim/ limxmem,limblks,limpair
      common /memmax/ ispblx, maxme1,iforwhat
c
      common /logic1/ ndege(1)
      common /logic2/ lenn(1)
      common /logic3/ lensm(1)
c
      dimension bl(*)
      dimension npar(*),ijbl(nbl2,*)
      dimension inx(12,*),iis(*),jjs(*)
c-----------------------------------------------------------
c  iforwhat shows where the Blockin2 routins are called from :
c    from twoint  with iforwhat=1  (ordinary integrals)
c    from Shift2  with iforwhat=2  (GIAO derivatives  )
c    from Force2  with iforwhat=3  (I-st derivatives  )
c    from ??????  with iforwhat=4  (IIed derivatives  )
c    from ptwoint with iforwhat=5  (mp2 integrals     )
c-----------------------------------------------------------
c
      ioffset=igetival('ioffset')
      memtot=lcore-ioffset
c
      if(memtot.lt.limxmem) then
        call nerror(911,'blkmemor','lcore smaller than limxmem:',
     *              memtot,limxmem)
      endif
c-----------------------------------------------------------
c remember limxmem
c
      limxnew=limxmem
      nblx=nbl2-1
      if(nblx.eq.0) nblx=1
c-----------------------------------------------------------
 4321 continue    ! return address for redefining LIMXMEM
c-----------------------------------------------------------
c08   do 100 ibl=nbl2,nbl2
      do 100 ibl=  1 ,nbl2
      ijpar=npar(ibl)
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
c
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
      call get_type1(itype,jtype, itype1,jtype1)
c
      nfij=lenn(itype1)*lenn(jtype1)
c
      lci=inx(5,ics1)-inx(1,ics1)
      lcj=inx(5,jcs1)-inx(1,jcs1)
c
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
      ngcij=(ngci+1)*(ngcj+1)
c
c08   do 100 kbl=nblx,nbl2
      do 100 kbl=1   ,ibl
c
      klpar=npar(kbl)
      klcs1=ijbl(kbl,1)
      kcs1=iis(klcs1)
      lcs1=jjs(klcs1)
c
      ktype=inx(12,kcs1)
      ltype=inx(12,lcs1)
      call get_type1(ktype,ltype, ktype1,ltype1)
c
      nfkl=lenn(ktype1)*lenn(ltype1)
      nfijkl=nfij*nfkl
c
      lck=inx(5,kcs1)-inx(1,kcs1)
      lcl=inx(5,lcs1)-inx(1,lcs1)
c
      ngck=inx(4,kcs1)
      ngcl=inx(4,lcs1)
      ngckl=(ngck+1)*(ngcl+1)
c-----------------------------------------------------------
c calculate maximum memory needed for one block :
c
c 1) for ordinary two-el.integrals (ifor=1)
c 2) for GIAO two-el. derivatives  (ifor=2)
c 3) for gradient derivatives      (ifor=3)
c 4) for second derivatives        (ifor=4)
c 5) for mp2/lmp2 integrals        (ifor=5)
c 6) for non-FTC  integrals        (ifor=11)
c 7) for non-FTC integrals derivatives (ifor=13)
c
      ifor=iforwhat
      if(ifor.eq.5) ifor=1
      if(ifor.eq.11) ifor=1
      if(ifor.eq.13) ifor=3
c
      IF( iroute.eq.1 ) THEN
         call blksize1(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *                 lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *                 memor2,memor2ij,memor2kl,memor4,ifor)
      ELSE
         call blksize2(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *                 lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *                 memor2,memor2ij,memor2kl,memor4,ifor)
      ENDIF
c
c output:
c           memor2,memor2ij,memor2kl,memor4
c-----------------------------------------------------------
c check if there is enough memory given to handle the minimum
c i.e. one quartet, one ij- & one kl-pair :
c
      mqrt1=memor2ij+memor2kl+memor4
      mqrt1=mqrt1+memor2  ! memor2 is an amount allocated in prec2ij,kl
c
c
      if(mqrt1.ge.limxnew) then
c       not enough memory (limxmem) for one block
        xneed=dble(mqrt1)*1.15d0 ! add 15% more
        needed=xneed
c      write(6,911) ibl,kbl, mqrt1, limxnew, needed
c      write(6,*)'     types=',itype1,jtype1,ktype1,ltype1
c      write(6,*)' functions=', nfij, nfkl
c 911 format('block2=',2i4,' quad1=',i8,' mem:old new= ',2i10)
        limxnew=needed
        go to 4321
      endif
c-----------------------------------------------------------
  100 continue
c--------------------------------------------------------------------
c    try to increase the limxmem limit but it must be consistent
c    with total memory given (lcore)
c
      if(limxnew.NE.limxmem) then
        call redefmem(memtot,limxnew, limxmem)
      endif
c--------------------------------------------------------------------
      end
c=======================================================================
      subroutine print2cd(ncs,schwarz,nbl2,npar,ijbl,iftc0)
      implicit real*8 (a-h,o-z)
      dimension schwarz(ncs,ncs)
      dimension npar(nbl2)
      dimension ijbl(nbl2,*)
      dimension iftc0(ncs)
c
cc       write(6,*) iftc0
c
      do ibl=1,nbl2
         ijpar=npar(ibl)
cc         write(6,*) '      pair-block no=',ibl,' pairs=',ijpar
cc         write(6,*) ' pair   shells    iftc  schwarz  '
         do ij=1,ijpar
            ijcs=ijbl(ibl,ij)
            call get_ij_half(ijcs, ics, jcs)
            ijcd=iftc0(ics)+iftc0(jcs)
            write(6,66) ijcs,ics,jcs,ijcd,schwarz(ics,jcs)
         enddo
      enddo
c
  66  format(2x,i3,2x,2(i3,1x),2x,i2,2x,1pe12.6)
      end
c=======================================================================
      subroutine ftc_bl4(ncs,nbl2,ijbl,iftc0,iccdd,nftcbl4,nbl4ftc)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      dimension ijbl(nbl2,*)
      dimension iftc0(ncs)
      dimension nftcbl4(*)
c
      nbl4ftc=0
      nbl4   =0
      do ibl=1,nbl2
         ijcs1=ijbl(ibl,1)
         call get_ij_half(ijcs1, ics1, jcs1)
         ijcd1=iftc0(ics1)+iftc0(jcs1)
c
         if(iccdd.eq.0 .and. ijcd1.eq.0 ) then
            nbl4=nbl4+ibl
            go to 9876
         endif
c
         do kbl=1,ibl
            nbl4=nbl4+1
            klcs1=ijbl(kbl,1)
            call get_ij_half(klcs1, kcs1, lcs1)
            klcd1=iftc0(kcs1)+iftc0(lcs1)
            if(iccdd.eq.0 .and. klcd1.eq.0 ) go to 9877
            if(ijcd1+klcd1.GE.2) then
               nbl4ftc= nbl4ftc+1
               nftcbl4(nbl4ftc)=nbl4
            endif
 9877       continue
         enddo
 9876    continue
      enddo
c
      write(ioutput,66) nbl4,nbl4ftc
 66   format(/' Number of blocks : Total=',i10,'  for FTC=',i10,/)
c
      end
c=======================================================================
