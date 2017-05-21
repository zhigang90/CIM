c===============================================================
      subroutine blocking1(iforwhat,datnuc,no_basis,list1,
     *                     ncs,inx,datbas, ncs_b,ncs_e,
     *                     nblock1,nblock1_back,nbl1,maxshell)
c---------------------------------------------------------------
c blocking procedure for single shells (contracted)
c
c Input
c---------
c iforwhat  - shows type of task (integrals to be calculated)
c datnuc()  - nuclear data
c no_basis  - basis set number (from 1 to 4)
c list1(*)  - list of shells (may be not consecutive)
c ncs       - total number of contr.shells in this basis set (no_basis)
c inx()     - basis set info for no_basis basis set
c datbas()  - basis set data for no_basis basis set
c ncs_b     - consider shells from ncs_b to ncs_e only
c ncs_e
c
c Output
c-------
c nblock1() - array showing :
c             nblock1(1,i)  : first shell in i-block1
c             nblock1(2,i)  : last  shell in i-block1
c nblock1_back() - shows single-shell block no a given shell belongs to
c                  nblock1_back(ics)--> iblock1
c nbl1      - number of blocks of single shells
c maxshell  - maximum number of single shells/block
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*4 text
      common /outfile/ ioutput
      dimension list1(*)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
      dimension nblock1(2,ncs),nblock1_back(ncs)
c---------------------------------------------------------------
      call getival('iroute',iroute)
c---------------------------------------------------------------
c     write(ioutput,*)
c    * ' -----------------------------------------------'
c     write(ioutput,*)
c    * ' Blocks of Single Shells for Basis Set No=',no_basis
c     write(ioutput,*)'     number of shells =',ncs
c     write(ioutput,*)'     shell considered = (',ncs_b,' :',ncs_e,')'
c---------------------------------------------------------------
c make blocks of single shells :
c
      call blk_shells_be(iforwhat,ncs, ncs_b,ncs_e,list1,inx,iroute,
     *                   datnuc,datbas, nblock1,nbl1,maxshell)
c---------------------------------------------------------------
c make nblock1_back(*) array
c
      call make_nbback1(ncs,nbl1,nblock1,nblock1_back)
c
c output : nblock1_back(ics) => block of shells ics belongs to
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_shells_be(iforwhat,ncs,ncs_b,ncs_e,list1,inx,
     *                 iroute, datnuc,datbas,nblock1,nbl1,maxshell)
      implicit real*8 (a-h,o-z)
      parameter (maxvalue=48)
      common /outfile/ ioutput
      logical txs93,txs95, cond1,cond2,cond3,cond4,cond5,cond6
      logical cond0
      dimension list1(*)
      dimension nblock1(2,*)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c-----------------------------------------------------------------
c Main routine to makes blocks of single shells.
c
c calculate number of different shells (i.e. different blocks)
c according to two different criterions (txs93 or txs95)
c-----------------------------------------------------------------
c NO limits on a number of shells per block
c-----------------------------------------------------------------
c Input :
c
c iforwhat - values from 1 to 6 (so far) shows the type of
c            the task integrals are needed for ( scf,force,etc)
c ncs      - number of contracted shells in the LIST
c ncs_b, ncs_e  - range of shells of interest (part of a whole
c                 basis set
c list1(*)    - list of shells to be considered here
c inx(12,*  ) - basis set info
c iroute      - 1 or to showing blocking procedure (tx93 or tx95)
c datnuc(5,*) - nuclear data
c datbas(13,*)- basis set data
c
c Output :
c
c nblock1() - array showing :
c             nblock1(1,i)  : first shell in i-block1
c             nblock1(2,i)  : last  shell in i-block1
c nbl1      - number of single-shell blocks
c maxshell  - maximum number of shells per block
c-----------------------------------------------------------------
c
      iexch=0
      ishells=0
      do 55 ics0=ncs_b,ncs_e-1
ccc   do 55 icsx=ncs_b,ncs_e-1   ! shells from the LIST
ccc      ics0=list1(icsx)
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
c next shell in a list :
c
ccccc    ics1=list1(icsx+1)
         ics1=      ics0+1
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
c check against criterion (3 or 5 conditions) :
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
c---------------------------------------------
           if(cond1 .and. iroute.eq.1) then
c             because of possible numerical instability
              call get_cond0(exi0,exi0l, exi1,exi1l,ist0,cond0)
           endif
c---------------------------------------------
c
           ics10=ics1-ics0
c
           ishells=ishells+1
C
           if(iroute.eq.1) then
c2002         if(.not.txs93 .or. ics10.gt.1 ) then
              if( (.not.txs93) .or. ics10.gt.1 .or. (.not.cond0) ) then
c                because of possible numerical instability
                 iexch=iexch+1
                 nblock1(2,iexch)=ics0
                 nblock1(1,iexch+1)=ics1
                 ishells=0
              endif
           else
              if(.not.txs95 .or. ics10.gt.1) then
                 iexch=iexch+1
                 nblock1(2,iexch)=ics0
                 nblock1(1,iexch+1)=ics1
                 ishells=0
              endif
           endif
   55 continue
c
      nbl1=iexch+1              ! number of shlell-blocks
c
cccc  nblock1(1,1)=list1(ncs_b)
cccc  nblock1(2,nbl1)=list1(ncs_e)
      nblock1(1,1)=      ncs_b
      nblock1(2,nbl1)=      ncs_e
c
c--------------------------------------------------------
      maxshell=0
      do ibl1=1,nbl1
         ibeg=nblock1(1,ibl1)
         iend=nblock1(2,ibl1)
         ishell=iend-ibeg+1
         maxshell=max(maxshell,ishell)
      enddo
c--------------------------------------------------------
c     write(ioutput,495) nbl1,ncs_e-ncs_b+1,ncs_b,ncs_e, maxshell
c 495 format( '      Number of Blocks of Contracted Shells =',i5/
c    *        '      number of considerd contr. shells     =',i5/
c    *        '      shell considered = (',i3,' :',i3,')'/
c    *        '      maximum number of shells/block        =',i5/)
ctest
c     do ibl1=1,nbl1
c        ibeg=nblock1(1,ibl1)
c        iend=nblock1(2,ibl1)
c        itype=inx(12,ibeg)
c        icont=inx(5,ibeg)-inx(1,ibeg)   ! contraction lenght
c        write(91,*)' blk1=',ibl1,' type=',itype,' contr.=',icont,
c    *              ' no of shells=',iend-ibeg+1
c        do ics=ibeg,iend
c           write(91,*)'    shells=',ics
c        enddo
c     enddo
ctest
c---------------------------------------------
      end
c=======================================================================
      subroutine blocking2_same(bl,called4,inx,
     *                             ds,
     *                             nbl1_1,nblock1_1,maxsh_1,
     *                             nbl1_2,nblock1_2,maxsh_2)
      implicit real*8 (a-h,o-z)
      character*7 called4
      common /intlim/ limxmem,limblks,limpair
      dimension bl(*)
      dimension ds(*)  ! screening "density" for mp2&mp2d integrals
      dimension inx(12,*)
      dimension nblock1_1(2,*),nblock1_2(2,*)
c---------------------------------------------------------------
c ASSUMING THAT ALL SHELLS BELONG TO THE SAME BASIS SET.
c---------------------------------------------------------------
c
c blocking procedure for pairs from two sets of blocks of single shells:
c
c Input :
c bl(*)    - buffer for everything
c called4  - shows if it is called for left side (IJ| or right side |KL)
c            pairs  ( can be 'ijpairs' or 'klpairs' )
c            This distinction is needed for different memory allocations
c
c Set of shell-blocks no 1 :
c  nbl1_1       - no of shell-blocks in set 1
c  nblock1_1(*) - array showing shells belonging to blocks:
c                 nblock1(1,i)  : first shell in i-block1
c                 nblock1(2,i)  : last  shell in i-block1
c  maxsh_1      - maximum number of shells in one block in set 1
c
c Set of shell-blocks no 2 : as above with ext. _2
c
c---------------------------------------------------------------
c total number of shell pairs :  ncspairs=ncsh_1*ncsh_2
c total number of pair-blocks WITHOUT limit: nbl2 = nbl1_1*nbl1_2
c maximum pairs in one pair-block:  maxpairs=maxsh_1*maxsh_2
c---------------------------------------------------------------
c Constructe blocks of contracted shell pairs; setup the arrays:
c nparx, ijblx
c
      call blk_pairx(bl,ds,called4,inx,limpair, 
     *               nbl1_1,nblock1_1, nbl1_2,nblock1_2)
c
c---------------------------------------------------------------
c This is the end of blocking procedure for PAIRS
c
c The following is known at this point :
c
c  pairs ijcs are given by ics=iis(ijcs) and jcs=jjs(ijcs)
c  number of pair-blocks = NBL2
c  number of pairs in the block ibl is nparx(ibl)
c  which pairs belong to this block : ijcs=ijblx(ibl,1-npar)
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairx(bl,ds,called4,inx,limpair, 
     *                     nbl1_1,nblock1_1,nbl1_2,nblock1_2)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*7 called4
      dimension bl(*)
      dimension ds(*)
      dimension inx(12,*)
      dimension nblock1_1(2,*),nblock1_2(2,*)
c
      data timereor /0.d0/
      save timereor
c---------------------------------------------------------------
      call getival('ibas',ibas)
      call getival('iroute',iroute)
c---------------------------------------------------------------
c Calculate number of pair-blocks NBL2 & maximum pairs/block MAXPAR
c with the limit=limpair . This is very inconvenient but needed to
c know how much memory to allocate. If there was no limit then this
c step would not be needed at all.
c
      call blk_pairs_calx(nbl1_1,nblock1_1, nbl1_2,nblock1_2,
     *                    iroute,bl(ibas),
     *                    inx,limpair,nbl2,maxpar,npartot)
c
c output : nbl2   - number of pair-blocks (with a limit)
c          maxpar - maximum number of pairs per block2
c          npartot- total number of pairs in all blocks
c---------------------------------------------------------------
      call getival('ncs',ncs)
      call getmem(ncs,idshell)
c---------------------------------------------------------------
c Predicted number of pair-blocks :
c
      nbl2_pred=nbl2
c---------------------------------------------------------------
c Allocate memory for NPARX and IJBLX
c
      call getint(nbl2_pred       ,nparx)
      call getint(nbl2_pred*maxpar,ijblx)
      call getint(npartot,map_ij_blx)
c---------------------------------------------------------------
c Constructe blocks of contracted pairs: setup arrays :
c       ijblx(nbl2,*) , nparx(nbl2) and map_ij_blx(*)
c
      call getival('schwarz',ischwarz)
cccc  call getival('ncs',ncs)
c
      call getival('mpres_in_bg',mpres_in_bg)
c
cccc  call getmem(ncs,idshell)
      call blk_pairs_makx(ncs,bl(ischwarz),
     *                    ds,bl(mpres_in_bg),bl(idshell),
     *                    nbl1_1,nblock1_1, nbl1_2,nblock1_2,
     *                    iroute,bl(ibas),
     *                    inx,limpair,nbl2_pred,nbl2,
     *                    bl(nparx),bl(ijblx),bl(map_ij_blx))
cccc  call retmem(1)
c
c nbl2 is a real number of pair-blocks (should be nbl2.LE.nbl2_pred)
c---------------------------------------------------------------
c reorder pairs in each block according to the schwarz integ. values :
c
      call secund(treo1)
      call blk_pairs_reor(ncs,bl(ischwarz),nbl2_pred,nbl2,
     *                    bl(nparx),bl(ijblx))
      call secund(treo2)
      timereor=timereor+treo2-treo1
      call setrval('timereor',timereor)
c---------------------------------------------------------------
c and save addresses and number of pair-blocks:
c
      if(called4.eq.'ijpairs') then
         call setival('nparx',nparx)
         call setival('ijblx',ijblx)
         call setival('mapijblx',map_ij_blx)
         call setival('blocksij',nbl2)
         call setival('blpredij',nbl2_pred)
      else
         call setival('npary',nparx)
         call setival('ijbly',ijblx)
         call setival('mapijbly',map_ij_blx)
         call setival('blockskl',nbl2)
         call setival('blpredkl',nbl2_pred)
      endif
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairs_calx(nbl1_1,nblock1_1, nbl1_2,nblock1_2,
     *                          iroute,datbas,
     *                          inx,limpair, nbl2,maxpar,npartot)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
      dimension nblock1_1(2,*),nblock1_2(2,*)
      dimension inx(12,*)
      dimension datbas(13,*)
c-----------------------------------------------------------------
c This routine estimates number of pair-blocks which is needed
c for memory allocation before pair-blocks are really constructed
c-----------------------------------------------------------------
c Input :
c nbl1_1 - number of single shell blocks in set 1
c nbl1_2 - number of single shell blocks in set 2
c nblock_1(*)  - showing beginning & end of a block of singl-shell
c nblock_2(*)
c limpair  - limit for number of shell-pairs per block
c
c Output :
c
c nbl2 - number of shell-pairs  blocks
c maxpar - maximum number of pairs in one block
c mpartot- total number of pairs
c-----------------------------------------------------------------
c calculate number of contr.shell pair-blocks WITH a limit
c
      ijpar_sum=0
      maxpar=0
      ijblock=0             !    blocks counter
      do 100 ibl=1,nbl1_1
         ibeg=nblock1_1(1,ibl)
         iend=nblock1_1(2,ibl)
         do 200 jbl=1,nbl1_2
            jbeg=nblock1_2(1,jbl)
            jend=nblock1_2(2,jbl)
cccc        call get_lim2(inx,ibeg,jbeg,limpair,limit2)
            call get_lim2(datbas,iroute,inx,ibeg,iend,jbeg,jend,
     *                    limpair,limit2)
            ijblock=ijblock+1
            ijpar=0
            do 300 ics=ibeg,iend
               iics=ics*(ics-1)/2
               do 400 jcs=jbeg,jend
                  if(ics.ge.jcs) then
                     ijcs=iics+jcs
                  else
                     ijcs=jcs*(jcs-1)/2 +ics
                  endif
c???............................................
c2001             if(ijpar.ge.limpair) then
                  if(ijpar.ge.limit2 ) then
                     ijblock=ijblock+1
                     ijpar=0
                  endif
c???............................................
                  ijpar_sum=ijpar_sum + 1
                  ijpar=ijpar+1
                  if(ijpar.gt.maxpar) maxpar=ijpar
  400          continue
  300       continue
  200    continue
  100 continue
c
      nbl2=ijblock
      npartot=ijpar_sum
c---------------------------------------------------------------
c     write(91,*)' blk_pairs_calx routine '
c     write(ioutput,495) nbl2,npartot,maxpar
c 495 format( ' Number of Blocks of Contracted Shell Pairs    =',i7/
c    *        ' total number of contracted shell pairs        =',i7/
c    *        ' maximum number of shell-pairs/block           =',i7)
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairs_makx(ncs,schwarz,
     *                          ds,mpres_in_bg,dshell,
     *                          nbl1_1,nblock1_1, nbl1_2,nblock1_2,
     *                          iroute,datbas,
     *                          inx,limpair,nbl2_pred,nbl2,
     *                          npar,ijbl,map_ij_bl2)
      implicit real*8 (a-h,o-z)
      logical do12,do21,dopa
      common /outfile/ ioutput
      common /neglect/ eps,eps1,epsr,eps8
      dimension nblock1_1(2,*),nblock1_2(2,*)
      dimension inx(12,*)
      dimension npar(nbl2_pred)                    ! output
      dimension ijbl(nbl2_pred,*)                  ! output
      dimension map_ij_bl2(*)                      ! output (ncspairs)
      dimension schwarz(ncs,ncs)
      dimension ds(ncs,ncs), dshell(ncs)
      dimension datbas(13,*)
c-----------------------------------------------------------------
      dimension mpres_in_bg(*)
c-----------------------------------------------------------------
c constructe pair-blocks of contracted shells
c
c-----------------------------------------------------------------
c Input :
c ncs    - number of contracted shells (total)
c schwarz(ncs,ncs) - an array with schwarz integrals
c nbl1_1 - number of single shell blocks in set 1
c nblock_1(*)  - showing beginning & end of a block of singl-shell
c nbl1_2 - number of single shell blocks in set 2
c nblock_2(*)  - showing beginning & end of a block of singl-shell
c limpair  - limit for number of shell-pairs per block
c nbl2_pred - predicted number of pair-blocks (for dimension only)
c
c Output :
c nbl2      - real number of pair-blocks
c npar()    - an array with number of pairs in each pair-block
c ijbl(nbl2_pred,*) - shows pairs belonging to a given pair-blocks
c                     ijbl(ijblock,ijpar)-> ijcs (=ics*(ics-1)/2+jcs)
c map_ij_bl2(*) - mapping from ijcs to ijblock
c
c-----------------------------------------------------------------
c        write(6,*) (mpres_in_bg(i),i=1,6)
c        write(6,*) '------------------------------'
c-----------------------------------------------------------------
      eps2=eps*eps
c---------------------------------------------------------------
c     write(ioutput,*)'-----------------------------------------'
c     write(ioutput,*)' Blocks of Shell-Pairs  (blk_pairs_makx) '
c
      ijpar_sum=0
      IF(NBL1_1.GE.NBL1_2) THEN
         ijblock=0             !    blocks counter
         do ibl=1,nbl1_2
            ibeg=nblock1_1(1,ibl)
            iend=nblock1_1(2,ibl)
            do jbl=1,ibl
               jbeg=nblock1_2(1,jbl)
               jend=nblock1_2(2,jbl)
ccc            call get_lim2(inx,ibeg,jbeg,limpair,limit2)
               call get_lim2(datbas,iroute,inx,ibeg,iend,jbeg,jend,
     *                       limpair,limit2)
               if(ibl.eq.jbl) then
                  ijblock=ijblock+1
                  call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock,ibl,ibeg,iend,jbl,jbeg,jend,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001                          limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
                  go to 100
               endif
c.................................................................
check if we need to do just 2-1 or 1-2 aslo or maybe 2-1 and small part
c
               ibl1=jbl
               i1beg=nblock1_1(1,jbl)
               i1end=nblock1_1(2,jbl)
               ibl2=ibl
               i2beg=ibeg
               i2end=iend
c
               jbl1=jbl
               j1beg=jbeg
               j1end=jend
               jbl2=ibl
               j2beg=nblock1_2(1,ibl)
               j2end=nblock1_2(2,ibl)
c
               call whattodo(i1beg,i1end, i2beg,i2end,
     *                       j1beg,j1end, j2beg,j2end,
     *                       ipa_beg,ipa_end,
     *                       jpa_beg,jpa_end,
     *                       do21,do12,dopa)
c              ...................................................
               if(do21) then
                  ijblock=ijblock+1
                  call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock, ibl2,i2beg,i2end,
     *                                     jbl1,j1beg,j1end,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
               endif
c              ...................................................
               if(do12) then
                  ijblock=ijblock+1
                  call make_1_bl2_ics(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock, ibl1,i1beg,i1end,
     *                                     jbl2,j2beg,j2end,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
               endif
c              ...................................................
               if(dopa) then
                  ijblock=ijblock+1
                  call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock, ibl ,ipa_beg,ipa_end,
     *                                     jbl ,jpa_beg,jpa_end,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
               endif
c.................................................................
c
  100          continue
            enddo
         enddo
c
         do ibl=nbl1_2+1,nbl1_1
            ibeg=nblock1_1(1,ibl)
            iend=nblock1_1(2,ibl)
            do jbl=1,nbl1_2
               jbeg=nblock1_2(1,jbl)
               jend=nblock1_2(2,jbl)
ccccccc        call get_lim2(inx,ibeg,jbeg,limpair,limit2)
               call get_lim2(datbas,iroute,inx,ibeg,iend,jbeg,jend,
     *                       limpair,limit2)
               ijblock=ijblock+1
               call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                         ijblock,ibl,ibeg,iend, jbl,jbeg,jend,
     *                      limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                      limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
               ijpar_sum=ijpar_sum+ npar(ijblock)
            enddo
         enddo
      ELSE
         ijblock=0
         do ibl=1,nbl1_1
            ibeg=nblock1_1(1,ibl)
            iend=nblock1_1(2,ibl)
            do jbl=1,ibl
               jbeg=nblock1_2(1,jbl)
               jend=nblock1_2(2,jbl)
ccccc          call get_lim2(inx,ibeg,jbeg,limpair,limit2)
               call get_lim2(datbas,iroute,inx,ibeg,iend,jbeg,jend,
     *                       limpair,limit2)
               if(ibl.eq.jbl) then
                  ijblock=ijblock+1
                  call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock,ibl,ibeg,iend,jbl,jbeg,jend,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
                  go to 200
               endif
c.................................................................
check if we need to do just 2-1 or 1-2 aslo or maybe 2-1 and small part
c
               ibl1=jbl
               i1beg=nblock1_1(1,jbl)
               i1end=nblock1_1(2,jbl)
               ibl2=ibl
               i2beg=ibeg
               i2end=iend
c
               jbl1=jbl
               j1beg=jbeg
               j1end=jend
               jbl2=ibl
               j2beg=nblock1_2(1,ibl)
               j2end=nblock1_2(2,ibl)
c
               call whattodo(i1beg,i1end, i2beg,i2end,
     *                       j1beg,j1end, j2beg,j2end,
     *                       ipa_beg,ipa_end,
     *                       jpa_beg,jpa_end,
     *                       do21,do12,dopa)
c              ...................................................
               if(do21) then
                  ijblock=ijblock+1
                  call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock, ibl2,i2beg,i2end,
     *                                     jbl1,j1beg,j1end,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
               endif
c              ...................................................
               if(do12) then
                  ijblock=ijblock+1
                  call make_1_bl2_ics(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock, ibl1,i1beg,i1end,
     *                                     jbl2,j2beg,j2end,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
               endif
c              ...................................................
               if(dopa) then
                  ijblock=ijblock+1
                  call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                            ijblock, ibl ,ipa_beg,ipa_end,
     *                                     jbl ,jpa_beg,jpa_end,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
                  ijpar_sum=ijpar_sum+ npar(ijblock)
               endif
c.................................................................
  200          continue
            enddo
         enddo
c
         do ibl=1,nbl1_1
            ibeg=nblock1_1(1,ibl)
            iend=nblock1_1(2,ibl)
            do jbl=nbl1_1+1,nbl1_2
               jbeg=nblock1_2(1,jbl)
               jend=nblock1_2(2,jbl)
ccccc          call get_lim2(inx,ibeg,jbeg,limpair,limit2)
               call get_lim2(datbas,iroute,inx,ibeg,iend,jbeg,jend,
     *                       limpair,limit2)
               ijblock=ijblock+1
               call make_1_bl2_jcs(ncs,schwarz,eps2,
     *                            mpres_in_bg,ds,dshell,
     *                         ijblock,ibl,ibeg,iend, jbl,jbeg,jend,
     *                         limit2 ,nbl2_pred,npar,ijbl,map_ij_bl2)
c2001*                         limpair,nbl2_pred,npar,ijbl,map_ij_bl2)
               ijpar_sum=ijpar_sum+ npar(ijblock)
            enddo
         enddo
      ENDIF
c---------------------------------------------------------------
c
c real number of pair-blocks
c
      nbl2=ijblock
c
      if(nbl2_pred.LT.nbl2) then
         write(ioutput,*)
     *  ' potential problem: more pair-blocks than predicted'
      endif
c---------------------------------------------------------------
c printing
c
c     do ijblock=1,nbl2
c        ijpar=npar(ijblock)
c        write(ioutput,*)
c    *   ' Block of Shell-pairs no=',ijblock,' npairs=',ijpar
c        do ijp=1,ijpar
c           ijcs=ijbl(ijblock,ijp)
c           call get_ij_half(ijcs,ics,jcs)
c           write(ioutput,*) '  shell-pairs :',ijcs,' shells=',ics,jcs
c        enddo
c     enddo
c---------------------------------------------------------------
c     write(ioutput,*)' Total number of pair-blocks=',nbl2
c     write(ioutput,*)'                   predicted=',nbl2_pred
c     write(ioutput,*)' Total number of shell pairs=',ijpar_sum
c---------------------------------------------------------------
      end
c===============================================================
      subroutine make_1_bl2(ncs,schwarz,eps2,
     *                            mpres_in_bg,muse,
     *                        ijblock,ibl,ibeg,iend, jbl,jbeg,jend,
     *                        limpair,nbl2,npar,ijbl,map_ij_bl2)
      implicit real*8 (a-h,o-z)
      logical doit
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncspairs)
      dimension schwarz(ncs,ncs)
c---------------------------------------------------------------
      dimension mpres_in_bg(*)
c---------------------------------------------------------------
c mpres_in_bg is a mpres_in_bg(ics)= 1 or 0
c muse shows how to use this list : if 1 apply to ICS, if 2 to JCS
c---------------------------------------------------------------
check if sigle-shell blocks ibl & jbl overlap :
c
            icase=0
            if(iend.le.jbeg .or. ibeg.ge.jend) then
c                       no overlaping
c              write(91,*)' case 1 : no overlaping '
               icase=1
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif  ! if(iend.le.jbeg .or. ibeg.ge.jend) then
            if(ibeg.eq.jbeg .and. iend.eq.jend) then
c                         i1                i3
c                          ...................
c                          ...................
c                         j1                j3
c
c              write(91,*)' case 2 : i3=j3,i1=j1 the same'
               icase=2
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=ibeg,ics
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
c
            if(iend.gt.jbeg .and. iend.le.jend .and. ibeg.le.jbeg) then
c              i1         i2     i3
c              ...................
c                          ...................
c                         j1     j2         j3
c
c              write(91,*)'case 3 : i3>j1, i3<=j3, i1<=j1'
               icase=3
               ijpar=0
               doit=.false.
               do ics=jbeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,ics
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,jbeg-1
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=jbeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=iend+1,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.ge.jend .and. ibeg.ge.jbeg) then
c                      i1         i2     i3
c                      ...................
c                ...................
c               j1     j2         j3
c
c               write(91,*)' case 4 : i3>j1, i3=>j3, i1=>j1'
               icase=4
               ijpar=0
               doit=.false.
               do ics=ibeg,jend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=ibeg,ics
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,ibeg-1
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=jend+1,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=ibeg,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.gt.jend .and. ibeg.lt.jbeg) then
c                    i  .............
c                    j    .......
c
c               write(91,*)' case 5 : i3>j1, i3>j3, i1<j1'
               icase=5
               ijpar=0
               doit=.false.
               do ics=jbeg,jend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,ics
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,jbeg-1
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=jend+1,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.lt.jend .and. ibeg.gt.jbeg) then
c                    i     .....
c                    j .............
c               write(91,*)' case 6 : i3>j1, i3<j3, i1>j1'
               icase=6
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=ibeg,ics
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=jbeg,ibeg-1
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,iend
                  if(muse.eq.1) ipresent=mpres_in_bg(ics)
                  do jcs=iend+1,jend
                     if(muse.eq.2) ipresent=mpres_in_bg(jcs)
                     x_ij_ij=schwarz(ics,jcs)*dble(ipresent)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
 1234       continue
            if(icase.eq.0) then
               write(91,*)' case=0'
               write(91,*)' ibl,jbl=',ibl,jbl,' i_b-e=',ibeg,iend,
     *                                        ' j_b-e=',jbeg,jend
            endif
c---------------------------------------------------------------
      end
c===============================================================
c two new versions of the above make_1_bl2 routine
c make_1_bl2_Ics where mpres_in_bg operates on the first shell Ics
c make_1_bl2_Jcs where mpres_in_bg operates on the secnd shell Jcs
c===============================================================
c===============================================================
      subroutine make_1_bl2_Ics(ncs,schwarz,eps2,
     *                          mpres_in_bg,ds,dshell,
     *                        ijblock,ibl,ibeg,iend, jbl,jbeg,jend,
     *                        limpair,nbl2,npar,ijbl,map_ij_bl2)
      implicit real*8 (a-h,o-z)
      logical doit
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncspairs)
      dimension schwarz(ncs,ncs)
      dimension ds(ncs,ncs), dshell(ncs)
c---------------------------------------------------------------
c this is the same as make_1_bl2_Jcs but mpres_in_bg operates on
c the first ICS shell
c---------------------------------------------------------------
      dimension mpres_in_bg(*)
c---------------------------------------------------------------
c mpres_in_bg is a mpres_in_bg(ics)= 1 or 0 showing if a shell is present
c---------------------------------------------------------------
         isr1=min(ibeg,jbeg)
         isr2=max(iend,jend)
c---------------------------------------------------------------
         do ics=ibeg,iend
            if(mpres_in_bg(ics).eq.1) then
cccccc        call absmax( ncs       ,ds(  1 ,ics),iiii,dsmax)
c             call absmax(isr2-isr1+1,ds(isr1,ics),iiii,dsmax)
c             dshell(ics)=dsmax
              iiii=idamax(isr2-isr1+1,ds(isr1,ics),1)
              dshell(ics)=ds(iiii+isr1-1,ics)
            else
              dshell(ics)=0.d0
            endif
         enddo
c---------------------------------------------------------------
ccccc             write(6,*)' _Ics is called'
c---------------------------------------------------------------
check if sigle-shell blocks ibl & jbl overlap :
c
            icase=0
            if(iend.le.jbeg .or. ibeg.ge.jend) then
c                       no overlaping
c              write(91,*)' case 1 : no overlaping '
               icase=1
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif  ! if(iend.le.jbeg .or. ibeg.ge.jend) then
            if(ibeg.eq.jbeg .and. iend.eq.jend) then
c                         i1                i3
c                          ...................
c                          ...................
c                         j1                j3
c
c              write(91,*)' case 2 : i3=j3,i1=j1 the same'
               icase=2
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=ibeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
c
            if(iend.gt.jbeg .and. iend.le.jend .and. ibeg.le.jbeg) then
c              i1         i2     i3
c              ...................
c                          ...................
c                         j1     j2         j3
c
c              write(91,*)'case 3 : i3>j1, i3<=j3, i1<=j1'
               icase=3
               ijpar=0
               doit=.false.
               do ics=jbeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=ibeg,jbeg-1
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=jbeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=iend+1,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.ge.jend .and. ibeg.ge.jbeg) then
c                      i1         i2     i3
c                      ...................
c                ...................
c               j1     j2         j3
c
c               write(91,*)' case 4 : i3>j1, i3=>j3, i1=>j1'
               icase=4
               ijpar=0
               doit=.false.
               do ics=ibeg,jend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=ibeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=ibeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,ibeg-1
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=jend+1,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=ibeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.gt.jend .and. ibeg.lt.jbeg) then
c                    i  .............
c                    j    .......
c
c               write(91,*)' case 5 : i3>j1, i3>j3, i1<j1'
               icase=5
               ijpar=0
               doit=.false.
               do ics=jbeg,jend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=ibeg,jbeg-1
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=jend+1,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.lt.jend .and. ibeg.gt.jbeg) then
c                    i     .....
c                    j .............
c               write(91,*)' case 6 : i3>j1, i3<j3, i1>j1'
               icase=6
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=ibeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=ibeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=jbeg,ibeg-1
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               do ics=ibeg,iend
                if(mpres_in_bg(ics).eq.1) then
                  do jcs=iend+1,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(ics)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
                endif
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
 1234       continue
            if(icase.eq.0) then
               write(91,*)' case=0'
               write(91,*)' ibl,jbl=',ibl,jbl,' i_b-e=',ibeg,iend,
     *                                        ' j_b-e=',jbeg,jend
            endif
c---------------------------------------------------------------
      end
c===============================================================
c===============================================================
      subroutine make_1_bl2_Jcs(ncs,schwarz,eps2,
     *                          mpres_in_bg,ds,dshell,
     *                        ijblock,ibl,ibeg,iend, jbl,jbeg,jend,
     *                        limpair,nbl2,npar,ijbl,map_ij_bl2)
      implicit real*8 (a-h,o-z)
      logical doit
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncspairs)
      dimension schwarz(ncs,ncs)
      dimension ds(ncs,ncs), dshell(ncs)
c---------------------------------------------------------------
c this is the same as make_1_bl2_Ics but mpres_in_bg operates on
c the econd JCS shell
c---------------------------------------------------------------
      dimension mpres_in_bg(*)
c---------------------------------------------------------------
c mpres_in_bg is a mpres_in_bg(ics)= 1 or 0
c---------------------------------------------------------------
ccccc             write(6,*)' _Jcs is called'
c---------------------------------------------------------------
         isr1=min(ibeg,jbeg)
         isr2=max(iend,jend)
c---------------------------------------------------------------
         do jcs=jbeg,jend
            if(mpres_in_bg(jcs).eq.1) then
ccccccccc     call absmax(    ncs    ,ds(  1 ,jcs),iiii,dsmax)
c             call absmax(isr2-isr1+1,ds(isr1,jcs),iiii,dsmax)
c             dshell(jcs)=dsmax
              iiii=idamax(isr2-isr1+1,ds(isr1,jcs),1)
              dshell(jcs)=ds(isr1+iiii-1,jcs)
            else
              dshell(jcs)=0.d0
            endif
         enddo
c---------------------------------------------------------------
c          do ics=1,ncs
c              dshell(ics)=1.0d0
c          enddo
c---------------------------------------------------------------
check if sigle-shell blocks ibl & jbl overlap :
c
            icase=0
            if(iend.le.jbeg .or. ibeg.ge.jend) then
c                       no overlaping
c              write(91,*)' case 1 : no overlaping '
               icase=1
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif  ! if(iend.le.jbeg .or. ibeg.ge.jend) then
            if(ibeg.eq.jbeg .and. iend.eq.jend) then
c                         i1                i3
c                          ...................
c                          ...................
c                         j1                j3
c
c              write(91,*)' case 2 : i3=j3,i1=j1 the same'
               icase=2
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                  do jcs=ibeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
c
            if(iend.gt.jbeg .and. iend.le.jend .and. ibeg.le.jbeg) then
c              i1         i2     i3
c              ...................
c                          ...................
c                         j1     j2         j3
c
c              write(91,*)'case 3 : i3>j1, i3<=j3, i1<=j1'
               icase=3
               ijpar=0
               doit=.false.
               do ics=jbeg,iend
                  do jcs=jbeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,jbeg-1
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=jbeg,iend
                  do jcs=iend+1,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.ge.jend .and. ibeg.ge.jbeg) then
c                      i1         i2     i3
c                      ...................
c                ...................
c               j1     j2         j3
c
c               write(91,*)' case 4 : i3>j1, i3=>j3, i1=>j1'
               icase=4
               ijpar=0
               doit=.false.
               do ics=ibeg,jend
                  do jcs=ibeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,iend
                  do jcs=jbeg,ibeg-1
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=jend+1,iend
                  do jcs=ibeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.gt.jend .and. ibeg.lt.jbeg) then
c                    i  .............
c                    j    .......
c
c               write(91,*)' case 5 : i3>j1, i3>j3, i1<j1'
               icase=5
               ijpar=0
               doit=.false.
               do ics=jbeg,jend
                  do jcs=jbeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,jbeg-1
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=jend+1,iend
                  do jcs=jbeg,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
            if(iend.gt.jbeg .and. iend.lt.jend .and. ibeg.gt.jbeg) then
c                    i     .....
c                    j .............
c               write(91,*)' case 6 : i3>j1, i3<j3, i1>j1'
               icase=6
               ijpar=0
               doit=.false.
               do ics=ibeg,iend
                  do jcs=ibeg,ics
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,iend
                  do jcs=jbeg,ibeg-1
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               do ics=ibeg,iend
                  do jcs=iend+1,jend
                     x_ij_ij=schwarz(ics,jcs)*dshell(jcs)
                     if(x_ij_ij.GE.EPS2) then
                        doit=.true.
                        call count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                                  nbl2,npar,ijbl,map_ij_bl2)
                     endif
                  enddo
               enddo
               if(.not.doit) ijblock=ijblock-1
               go to 1234
            endif
 1234       continue
            if(icase.eq.0) then
               write(91,*)' case=0'
               write(91,*)' ibl,jbl=',ibl,jbl,' i_b-e=',ibeg,iend,
     *                                        ' j_b-e=',jbeg,jend
            endif
c---------------------------------------------------------------
      end
c===============================================================
c===============================================================
      subroutine count_subs(ics,jcs,ijblock,ijpar,limpair,
     *                      nbl2,npar,ijbl,map_ij_bl2)
      dimension npar(nbl2)                    ! output
      dimension ijbl(nbl2,*)                  ! output
      dimension map_ij_bl2(*)                 ! output (ncspairs)
c
                     if(ics.ge.jcs) then
                        ijcs=ics*(ics-1)/2 +jcs
                     else
                        ijcs=jcs*(jcs-1)/2 +ics
                     endif
                     if(ijpar.ge.limpair) then
                        ijblock=ijblock+1
                        ijpar=0
                     endif
                     ijpar=ijpar+1
                     npar(ijblock)=ijpar
                     ijbl(ijblock,ijpar)=ijcs
                     map_ij_bl2(ijcs)=ijblock
c
      end
c=======================================================================
      subroutine whattodo(i1beg,i1end, i2beg,i2end,
     *                    j1beg,j1end, j2beg,j2end,
     *                    ipa_beg,ipa_end,
     *                    jpa_beg,jpa_end,
     *                    do21,do12,dopa)
      logical do12,do21,dopa
      logical i1_same_j1, i2_same_j2
      logical i1_diff_j1, i2_diff_j2
c
      do21=.false.
      do12=.false.
      dopa=.false.
c
c compare i1 & j1 blocks as well as i2 & j2
c
      i1_same_j1=.false.      ! means identical
      i2_same_j2=.false.
c
      i1_diff_j1=.false.      ! means compeltly different
      i2_diff_j2=.false.
c
c  if none of above is true then they overlap
c
      if(i1beg.eq.j1beg .and. i1end.eq.j1end) i1_same_j1=.true.
      if(i2beg.eq.j2beg .and. i2end.eq.j2end) i2_same_j2=.true.
c
      if(i1beg.gt.j1end .or. i1end.lt.j1beg) i1_diff_j1=.true.
      if(i2beg.gt.j2end .or. i2end.lt.j2beg) i2_diff_j2=.true.
c
      if(i1_same_j1 .and. i2_same_j2) then
         do21=.true.
         return
      endif
c
      if(i1_diff_j1 .or. i2_diff_j2) then
         do21=.true.
         do12=.true.
         return
      endif
c12
      if(i1_same_j1 .and. (.not.i2_diff_j2) ) then
         i2long=i2end-i2beg+1
         j2long=j2end-j2beg+1
         if(i2long.gt.j2long) then
            do21=.true.
            do12=.false.
            dopa=.false.
         else
            do21=.false.
            do12=.true.
            dopa=.false.
         endif
         return
      endif
c34
      if(i2_same_j2 .and. (.not.i1_diff_j1) ) then
         i1long=i1end-i1beg+1
         j1long=j1end-j1beg+1
         if(i1long.gt.j1long) then
            do21=.false.
            do12=.true.
            dopa=.false.
         else
            do21=.true.
            do12=.false.
            dopa=.false.
         endif
         return
      endif
c56,78
      if( (.not.i1_diff_j1) .and.  (.not.i2_diff_j2) ) then
         i1long=i1end-i1beg+1
         j1long=j1end-j1beg+1
         i2long=i2end-i2beg+1
         j2long=j2end-j2beg+1
         if(i1long.gt.j1long .and. i2long.lt.j2long) then
            do21=.false.
            do12=.true.
            dopa=.false.
         endif
         if(i1long.lt.j1long .and. i2long.gt.j2long) then
            do21=.true.
            do12=.false.
            dopa=.false.
         endif
         if(i1long.gt.j1long .and. i2long.gt.j2long) then
            do21=.true.
            do12=.false.
            dopa=.true.
            ipa_beg=i1beg
            ipa_end=j1beg-1
            jpa_beg=j2beg
            jpa_end=j2end
         endif
         if(i1long.lt.j1long .and. i2long.lt.j2long) then
            do21=.true.
            do12=.false.
            dopa=.true.
            ipa_beg=i1beg
            ipa_end=i1end
            jpa_beg=i2end+1
            jpa_end=j2end
         endif
         return
      endif
c
      end
c=======================================================================
      subroutine blksizer_new(ijbl,nbl2_ijd,nbl2_ij,npar_ij,inx_i,inx_j,
     *                        klbl,nbl2_kld,nbl2_kl,npar_kl,inx_k,inx_l,
     *                        nblock4, nbl4 )
c----------------------------------------------------------------------
c This routine constructs basic data for blocks of quartets
c----------------------------------------------------------------------
c Input  :
c------------
c ijbl(nbl2_ijd,*) - an array showing pairs of cont.shells making
c                    a given ijblok (of pairs):
c                    ijbl(ijblock,ijpar)->ijcs  (=ics*(ics-1)/2+jcs)
c nbl2_ijd        - predicted number of left side pair-blocks IJ
c nbl2_ij         - real      number of left side pair-blocks IJ
c npar_ij(*)      - an array showing number of pairs in
c                   a given pair-block : npar_ij(ijblock)->nparij
c inx_i(12,*)     - basis set data where the ICS shell belong
c inx_j(12,*)     - basis set data where the JCS shell belong
c
c the same for right side pair-blocks KL
c
c klbl(nbl2_kld,*)
c nbl2_kld
c nbl2_kl
c npar_kl(*)
c inx_k(12,*)
c inx_l(12,*)
c
c
c Output :
c------------
c nbl4          - number of quartet-blocks
c nblock4(2,*)  - an array showing pair-blocks making
c                 up a given quartet-block :
c                 nblock4(1,ijklblock)-> ijblock
c                 nblock4(2,ijklblock)-> klblock
c maximum sizes for labels & int.buffer :
c     common /intbuf/ maxibuf,maxindx
c----------------------------------------------------------------------
c NO limits for the size of a block (i.e. number of quartets)
c----------------------------------------------------------------------
c This is good ONLY for two cases :
c (1) all IJ pair-blocks and all KL pair-blocks are THE SAME
c     which is true if
c     (a) full and the same basis set for all four indecies i,j,k,l)
c     (b) the same basis set for all indeces and
c         only one shell for I, all for J, one for K (the same as for I)
c         and all for L
c (2) IJ pair-blocks and all KL pair-blocks are DIFFERENT
c     but the same number of IJ & KL pair-blocks from the same basis set
c     this is true if
c     (a) one shell Ics for I , all for J, one shell Kcs.NE.Ics for K
c         and all for L
c     (b) all shells for I & J, one for K and another one for L
c     (c) another way around
c
c
c In (1a),(1b) and (2b) cases above the number of IJ & KL
c pair-blocks is the same.
c
c----------------------------------------------------------------------
c nbl2_ij & _kl are real numbers of pair-blokcs
c nbl2_ijd & _kld  are predicted ones needed here for dimensions
c----------------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
      logical ij_same_kl, ij1_same_kl1, ij2_same_kl2,ij1_same_kl2
      common /outfile/ ioutput
      common /cpu/ intsize,iacc,icache,memreal
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
c
      common /logic2/ lenn(1)
c
c     dimension npar(*),ijbl(nbl2,*)
c     dimension inx(12,*),iis(*),jjs(*)
c
      dimension inx_i(12,*),inx_j(12,*)
      dimension inx_k(12,*),inx_l(12,*)
      dimension ijbl(nbl2_ijd,*),npar_ij(nbl2_ijd)
      dimension klbl(nbl2_kld,*),npar_kl(nbl2_kld)
c
      dimension nblock4(2,*)
      data nquartot /0/
      save nquartot
c-----------------------------------------------------------
      call getmem(0,last1)
      call retmem(1)
c-----------------------------------------------------------
c     write(ioutput,*)'----------------------------------------'
c     write(ioutput,*)' Blocks of Shell-Quartets (blksizer_new)'
c
c-----------------------------------------------------------
c  Calculates the total number of blocks of contracted
c  shells quartets NBLOKS (an output parameter)
c-----------------------------------------------------------
c
      nquartets=0  ! quartets counter
      maxibuf=0    ! find the maximum size of a integral buffer
      maxindx=0    ! find the maximum size of labels array
      nblsmax=0
      nbloks=0     ! counter of all small blocks
      nbl12=0      ! counter of super-blocks
c-------------------------------------------------------------
c loop twice for parts 1 & 2 :
c
c  if(nbl2_ij.GE.nbl2_kl)              if(nbl2_ij.LT.nbl2_kl)
c
c      nbl2_kl                                nbl2_kl
c     ---------                         -------------------
c     | \     |                         | \     |         |
c     |   \   |                         |   \   |   2     |
c     |  1  \ |                         | 1   \ |         |
c     |-------|                         -------------------
c     |  2    |
c     |       |
c     ---------
c
c------------------------------------------------------------------
c        write(91,*)' blksizer: nbl2_ij &_kl=',nbl2_ij,nbl2_kl
c        call f_lush(91)
c------------------------------------------------------------------
c The simplest case : at least one of NBL2_IJ or NBL2_KL is equal 1
c
      if( nbl2_ij.eq.1 .or. nbl2_kl.eq.1 ) then
         do ibl=1,nbl2_ij
            ijpar=npar_ij(ibl)
            do kbl=1,nbl2_kl
               klpar=npar_kl(kbl)
c              ..........................................
c              check if ibl is the same as kbl :
               call arebl2same(ibl,ijbl,nbl2_ijd,ijpar,
     *                         kbl,klbl,nbl2_kld,klpar,
     *                         ij_same_kl)
c              ..........................................
               nbl12=nbl12+1
               nbloks=nbloks+1
               if(ij_same_kl) then
c......           nquart=ijpar*(ijpar+1)/2
                  nblock4(1,nbl12)=-ibl
                  nblock4(2,nbl12)=-kbl
c                 minus means diagonal block
               else
c......           nquart=ijpar*klpar
                  nblock4(1,nbl12)=ibl
                  nblock4(2,nbl12)=kbl
               endif
c......        nquartets=nquartets + nquart
            enddo
         enddo
         go to 999
      endif
c------------------------------------------------------------------
      if(nbl2_ij.GE.nbl2_kl) then
         iloop1_e=nbl2_kl
c
         iloop2_b=nbl2_kl+1
         iloop2_e=nbl2_ij
         kloop2_b=1
         kloop2_e=nbl2_kl
      else
         iloop1_e=nbl2_ij
c
         iloop2_b=1
         iloop2_e=nbl2_ij
         kloop2_b=nbl2_ij+1
         kloop2_e=nbl2_kl
      endif
c-------------------------------------------------------------
c loop 1 :
c
      do 100 ibl=1,iloop1_e
c        write(91,*)' loop 1 ; ibl=',ibl
c        call f_lush(91)
         ijpar=npar_ij(ibl)
         do 100 kbl=1,ibl
            klpar=npar_kl(kbl)
c           check if ibl is the same as kbl
            call arebl2same(ibl,ijbl,nbl2_ijd,ijpar,
     *                      kbl,klbl,nbl2_kld,klpar,
     *                      ij_same_kl)
c...........................................................
            if(ibl.eq.kbl) then
               nbl12=nbl12+1
               nbloks=nbloks+1
               if(ij_same_kl) then
c......           nquart=ijpar*(ijpar+1)/2
                  nblock4(1,nbl12)=-ibl
                  nblock4(2,nbl12)=-kbl
c                 minus means diagonal block
               else
c......           nquart=ijpar*klpar
                  nblock4(1,nbl12)=ibl
                  nblock4(2,nbl12)=kbl
               endif
c......        nquartets=nquartets + nquart
               go to 888
            endif
c...........................................................
c>>>>>>>>>  if(ibl.GT.kbl) then
c              do we need ij2-kl1 only or ij1-kl2 as well ?
c              compare ij1 with kl1 and ij2 with kl2 :
c
               ibl_1=kbl
               kbl_1=kbl
               ijpar_1=npar_ij(kbl)
               klpar_1=klpar
               call arebl2same(ibl_1,ijbl,nbl2_ijd,ijpar_1,
     *                         kbl_1,klbl,nbl2_kld,klpar_1,
     *                         ij1_same_kl1)
c
               ibl_2=ibl
               kbl_2=ibl
               ijpar_2=ijpar
               klpar_2=npar_kl(ibl)
               call arebl2same(ibl_2,ijbl,nbl2_ijd,ijpar_2,
     *                         kbl_2,klbl,nbl2_kld,klpar_2,
     *                         ij2_same_kl2)
c
c              we always need ij2-kl1 :
c
               nbl12=nbl12+1
               nbloks=nbloks+1
               if(ij_same_kl) then
c......           nquart_21=ijpar*(ijpar+1)/2
                  nblock4(1,nbl12)=-ibl
                  nblock4(2,nbl12)=-kbl
               else
c......           nquart_21=ijpar*klpar
                  nblock4(1,nbl12)=ibl
                  nblock4(2,nbl12)=kbl
               endif
c
c......        nquartets=nquartets + nquart_21
c
               if(ij1_same_kl1 .and. ij2_same_kl2) then
c                 we need only ij2-kl1
               else
c                 we need ij1-kl2 also
                  ikbl=2
c                 check if it is diagonal
                  call arebl2same(ibl_1,ijbl,nbl2_ijd,ijpar_1,
     *                            kbl_2,klbl,nbl2_kld,klpar_2,
     *                            ij1_same_kl2)
                  nbl12=nbl12+1
                  nbloks=nbloks+1
                  if(ij1_same_kl2) then
c......              nquart_12=ijpar_1*(ijpar_1+1)/2
                     nblock4(1,nbl12)=-kbl
                     nblock4(2,nbl12)=-ibl
                  else
c......              nquart_12=ijpar_1*klpar_2
                     nblock4(1,nbl12)=kbl
                     nblock4(2,nbl12)=ibl
                  endif
c......           nquartets=nquartets + nquart_12
               endif
c>>>>>>>>>  endif                   !  (ibl.NE.kbl)
c-----------------------------------------------------------
  888 continue
c
  100 continue
c--------------------------------------------------------------------
c
c  loop 2 :
c
      do 200 ibl=iloop2_b,iloop2_e
         ijpar=npar_ij(ibl)
         do 200 kbl=kloop2_b,kloop2_e
            klpar=npar_kl(kbl)
c           ..........................................
c           check if ibl is the same as kbl :
            call arebl2same(ibl,ijbl,nbl2_ijd,ijpar,
     *                      kbl,klbl,nbl2_kld,klpar,
     *                      ij_same_kl)
c           ..........................................
            nbl12=nbl12+1
            nbloks=nbloks+1
            if(ij_same_kl) then
c......        nquart=ijpar*(ijpar+1)/2
               nblock4(1,nbl12)=-ibl
               nblock4(2,nbl12)=-kbl
            else
c......        nquart=ijpar*klpar
               nblock4(1,nbl12)=ibl
               nblock4(2,nbl12)=kbl
            endif
c......     nquartets=nquartets + nquart
c           ..........................................
c--------------------------------------------------------------------
  200 continue
c--------------------------------------------------------------------
  999 continue
c--------------------------------------------------------------------
ctry it
c
c
      nbl4=nbl12    ! number of quartet-blocks
c
c--------------------------------------------------------------------
c
c Find maximum sizes for labels & integral buffer
c
c (Print out info concerning blocks of contracted shell quartets)
c--------------------------------------------------------------------
c        write(91,*)' blksizer: nbl2_ij &_kl=',nbl2_ij,nbl2_kl
c     write(ioutput,*)'---------- Quartet-Blocks ----------'
c
      isumq=0
      do ibl4=1,nbl4
         ibl=nblock4(1,ibl4)
         kbl=nblock4(2,ibl4)
         if(ibl.le.0) then
cccc        diagonal block
            ibl=-ibl
            kbl=-kbl
            ijpar=npar_ij(ibl)
            nquart=ijpar*(ijpar+1)/2
c           write(ioutput,*)'   Quart-Block no=',ibl4,
c    *      ' ; pair-blocks=',ibl,kbl,' nquarts=',nquart,'  diagonal',
c    *      ' ijpar=',ijpar
         else
            ijpar=npar_ij(ibl)
            klpar=npar_kl(kbl)
            nquart=ijpar*klpar
c           write(ioutput,*)'   Quart-Block no=',ibl4,
c    *      ' ; pair-blocks=',ibl,kbl,' nquarts=',nquart,
c    *      ' ijpar,klpar=',ijpar,klpar
         endif
c
         call get_nij_data(ijbl(ibl,1),inx_i,inx_j,nfij,ngcij)
         call get_nij_data(klbl(kbl,1),inx_k,inx_l,nfkl,ngckl)
c
         isumq=isumq+nquart
         nblsmax=max(nblsmax,nquart)
         labsize=4*nquart*ngcij*ngckl + nquart +4
         maxindx=max(maxindx, labsize)
         maxibuf=max(maxibuf, labsize*nfij*nfkl)
      enddo
c--------------------------------------------------------------------
c     if(isumq.ne.nquartets) then
c        write(6,*)
c    *' WARRNING : total number of contracted quartetsin error'
c        write(6,*) ' estimated=',nquartets,' true=',isumq
c        write(6,*)' blksizer: nbl2_ij &_kl=',nbl2_ij,nbl2_kl
c     endif
c--------------------------------------------------------------------
c printing
c
c     write(ioutput,900)
c..   write(ioutput,*)  ' true number of quartets in all blocks=',isumq
c     write(ioutput,496) nbl12,nbloks,nbloks/nbl12,nquartets,
c    *                   nquartets/nbloks, nblsmax
c 496 format( ' Number of Super-Blocks of Contracted Quartets =',i20/
c    *        ' number of Small-Blocks of Contracted Quartets =',i20/
c    *        '           small-blocks/super-block on average =',i20/
c    *        ' Total number of Contracted Shells    Quartets =',i20/
c    *        '    average small-block size (in quartets)     =',i20/
c    *        '    maximum small-block size (in quartets)     =',i20)
c     write(ioutput,900)
c 900 format('=================================================')
c................................
c
c     write(ioutput,497) maxibuf, maxindx
c 497 format(/' Maximum lenght of the integral buffer         =',i20/
c    *        ' and corresponding array for labels (*4)       =',i20/)
c--------------------------------------------------------------------
c     if(maxme1.ge.memaval) then
c        ioffset=igetival('ioffset')
c        irequ=maxme1+last12-ioffset
c        iavail=lcore-ioffset
c        call nerror(925,'blksizer',' Memory problem ',irequ,iavail)
c     endif
c--------------------------------------------------------------------
      end
c=======================================================================
      subroutine arebl2same(ijblock,ijbl,nbl2_ijd,ijpar,
     *                      klblock,klbl,nbl2_kld,klpar,
     *                      ij_same_kl)
c--------------------------------------------------------------------
c check if two pair-blocks (left- and rght-side ) are the same
c
c input :
c ijblock - left-side IJ block number
c ijbl(nbl2_ijd,*) - see routine blksizer_new for description
c nbl2_ijd         -       ---"------
c ijpar            - number of pairs in ijblock
c klblock - rightside KL block number
c nbl2_kld         -       ---"------
c klpar            - number of pairs in klblock
c
c Output:
c logical ij_same_kl
c--------------------------------------------------------------------
c
      logical ij_same_kl
      dimension ijbl(nbl2_ijd,*)
      dimension klbl(nbl2_kld,*)
c--------------------------------------------------------------------
      ijcs1=ijbl(ijblock,1)       !  first pair
      ijcs2=ijbl(ijblock,ijpar)   !  last pair
c
      klcs1=klbl(klblock,1)       !  first pair
      klcs2=klbl(klblock,klpar)   !  last pair
c--------------------------------------------------------------------
            ij_same_kl=.false.
            if(ijpar.eq.klpar) then
               if(ijcs1.eq.klcs1 .and. ijcs2.eq.klcs2) then
                  ij_same_kl=.true.
               endif
            endif
c-----------------------------------------------------------
      end
c=======================================================================
      subroutine get_nij_data(ijcs1,inx_i,inx_j,nfij,ngcij)
      common /logic2/ lenn(1)
      dimension inx_i(12,*),inx_j(12,*)
c
            call get_ij_half(ijcs1,ics1,jcs1)
            itype=inx_i(12,ics1)
            jtype=inx_j(12,jcs1)
            call get_type1(itype,jtype, itype1,jtype1)
            nfij=lenn(itype1)*lenn(jtype1)
            ngci=inx_i(4,ics1)
            ngcj=inx_j(4,jcs1)
            ngcij=(ngci+1)*(ngcj+1)
c
      end
c=======================================================================
      subroutine blockin4_n(isbl,ibl,kbl,bl)
c---------------------------------------------
c this routine is called for each super-block ISBL
c---------------------------------------------
c  nbloks=1   NO small blocks
c---------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
      common /memors/ nsym,ijshp,isymm
      dimension bl(*)
c------------------------------------------------
c  set up vectors : nqrt,nibl,nkbl, nijb,nije, nklb,nkle
c
      nbloks=1
      call memo2(nbloks)  ! res.mem. for 7 arrays (above)
c
      call getival('nparx',npar_ij)
      call getival('npary',npar_kl)
c
      call getival('nblock4',nblock4)
c
      call blockqur_n(isbl,maxqrt,bl(npar_ij),bl(npar_kl),
     *                bl(nblock4),
     *                bl(nqrtd),bl(nibld),bl(nkbld),bl(nijbd),bl(nijed),
     *                bl(nklbd),bl(nkled) )
c
cccc  call memo3(maxqrt)
      call getint(maxqrt,isymm)
c
      end
c===============================================================
      subroutine blockqur_n(isbl,maxqrt,  npar_ij,npar_kl,
     *                      nblock4,
     *                      nqrt,nibl,nkbl, nijb,nije,nklb,nkle)
      dimension npar_ij(*), npar_kl(*)
      dimension nblock4(2,*)
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension nqrt(*)
c-----------------------------------------------------------
c input :
c
c isbl        - number of a quartet block
c nblock4(2,*)
c npar(ibl) - number of pairs in ibl pair-block
c
c Output :
c
c maxqrt - maximum number of quartets in one block
c arrays :
c nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*),nqrt(*)
c
c-----------------------------------------------------------
      ibl=nblock4(1,isbl)
      kbl=nblock4(2,isbl)
c-----------------------------------------------------------
      ikbl=0
      if(ibl.lt.0) then
c        diagonal block
         ibl=-ibl
         kbl=-kbl
         ijpar=npar_ij(ibl)
         nquart=ijpar*(ijpar+1)/2
         ikbl=ikbl+1
         nibl(ikbl)=ibl
         nkbl(ikbl)=kbl
         nijb(ikbl)=1
         nije(ikbl)=ijpar
         nklb(ikbl)=1
         nkle(ikbl)=0
         nqrt(ikbl)=nquart
      else
         ijpar=npar_ij(ibl)
         klpar=npar_kl(kbl)
         nquart=ijpar*klpar
         ikbl=ikbl+1
         nibl(ikbl)=ibl
         nkbl(ikbl)=kbl
         nijb(ikbl)=1
         nije(ikbl)=ijpar
         nklb(ikbl)=1
         nkle(ikbl)=klpar
         nqrt(ikbl)=nquart
      endif
c
      maxqrt=nquart
c------------------------------------------------------------------
      end
c===============================================================
      subroutine blk_pairs_reor(ncs,schwarz,nbl2_pred,nbl2,npar,ijbl)
      implicit real*8 (a-h,o-z)
      dimension npar(nbl2_pred)                    ! input
      dimension ijbl(nbl2_pred,*)                  ! input/output
      dimension schwarz(ncs,ncs)                   ! input
c-----------------------------------------------------------------
      do ijblock=1,nbl2
         ijpar=npar(ijblock)
  100    continue
       iexch=0
       do ijp1=1,ijpar-1
          ijcs1=ijbl(ijblock,ijp1)
          call get_ij_half(ijcs1,ics,jcs)
          xij1=schwarz(ics,jcs)
          ijp2=ijp1+1
          ijcs2=ijbl(ijblock,ijp2)
          call get_ij_half(ijcs2,ics,jcs)
          xij2=schwarz(ics,jcs)
            if(xij2.gt.xij1) then
             iexch=1
             ijbl(ijblock,ijp1)=ijcs2
             ijbl(ijblock,ijp2)=ijcs1
            endif
         enddo
         if(iexch.eq.1) go to 100
      enddo
c-----------------------------------------------------------------
c     do ijblock=1,nbl2
c        ijpar=npar(ijblock)
c        do ijp=1,ijpar
c           ijcs=ijbl(ijblock,ijp)
c           call get_ij_half(ijcs,ics,jcs)
c           xij=schwarz(ics,jcs)
c           write(91,*)' ijblk=',ijblock,' ijpar=',ijp,' sch=',xij
c        enddo
c     enddo
c-----------------------------------------------------------------
      end
c===============================================================
      subroutine make_nbback1(ncs,nbl1,nblock1,nblock1_back)
      dimension nblock1(2,ncs),nblock1_back(ncs)
c construct nblock1_back() array
      do ibl=1,nbl1
         icsb=nblock1(1,ibl)
         icse=nblock1(2,ibl)
         do ics=icsb,icse
            nblock1_back(ics)=ibl
         enddo
      enddo
      end
c===============================================================
ccccc subroutine get_lim2(inx,ibeg,jbeg,limpair,limit2)
      subroutine get_lim2(datbas,iroute,inx,ibeg,iend,jbeg,jend,
     *                    limpair,limit2)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension datbas(13,*)
      common /memmax/ ispblx, maxme1,iforwhat
c......................................................................
c iforwhat can be 5 (mp2 integrals) or 6 mp2 gradient integrals
c......................................................................
c     ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28)
c
c     put a limit on number of pairs per pair-block
c     if a pair containes g- or higher orbitals then
c     restrict to only 2 or 3 pairs
c     (too high memory in memo4b)
c
c     if this is not a case just take limpair
c......................................................................
c general contraction (0 for segmented basis set)
c     ngci=inx(4,ibeg)
c     ngcj=inx(4,jbeg)
c     ngcij=max(ngci,ngcj)
c......................................................................
      itype=inx(12,ibeg)
      jtype=inx(12,jbeg)
      ijmax=max(itype,jtype)
      ijmin=min(itype,jtype)
c
      limit2=limpair
      if(ijmax.gt.7) then
c        g- and higher orbitals
         if(ijmin.eq.1) limit2=min(limpair,4)
         if(ijmin.eq.2) limit2=min(limpair,3)
         if(ijmin.ge.3) limit2=min(limpair,2)
      endif
      if(iforwhat.eq.6) then ! gradient itegrals
      if(ijmax.ge.6) then
c        f- and higher orbitals
         if(ijmin.eq.1) limit2=min(limpair,4)
         if(ijmin.eq.2) limit2=min(limpair,3)
         if(ijmin.ge.3) limit2=min(limpair,2)
      endif
      endif
c........................................................
c because of possible instability : one pair/block2
c this is important for iroute=1 ONLY
c........................................................
      if(iroute.eq.2) return
      call getival('stab',istab)
      if(istab.gt.0) return
c........................................................
      boamax=0.d0
c???  rapbmax=0.d0
      do ics=ibeg,iend
         isp1=inx(1,ics)+1               ! first primitive
         isp2=inx(5,ics)                 ! last  primitive
         do jcs=jbeg,jend
            jsp1=inx(1,jcs)+1               ! first primitive
            jsp2=inx(5,jcs)                 ! last  primitive
c
            do isp=isp1,isp2
               aa=datbas(1,isp)
               aa1=1.d0/aa
               do jsp=jsp1,jsp2
                  bb=datbas(1,jsp)
                  boamax=max(bb*aa1,boamax)
c???              rapb1=1.d0/(aa+bb)
c???              rapbmax=max(rapbmax,rapb1)
               enddo
            enddo
         enddo
      enddo
c
      boamax=boamax+boamax
c     10**4.5=31622.8
      if(boamax.gt.31622.8d0 .and. ijmax.ge.4) limit2=1 ! d and higher
c
      end
