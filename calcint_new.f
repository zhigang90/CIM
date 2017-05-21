c===================================================================
      subroutine calcint2_new(
     *                    isupb_r,bl,inx,thres, denspar, wherex ,labels,
     *                    ibuffz, nblsiz,nintez,ngctoz, moreint,stopnow)
      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 wherex
      character*4 schw4
c----------------------------------------------------------------
c Input parameters :
c
c isupb_r - requested super-block , ONE small block of it will be
c           calculated at once
c bl() - array for everything
c inx(12,ncs) - basis set and nuclear data info
c thres  - integral threshold
c denspar-density matrix over contracted shells (densp(ics,jcs)
c wherex -(character*4) shows the destination of integrals
c          wherex=where together with price tells the program which
c          blocks of integrals should be calculated
c
c Output parameters :
c
c labels  - reutrning labes info
c ibuffz  - address of integ.buffer (ibuf)
c nblsiz  - block-size(nbls)
c nintez  - number of integrals in each quartet
c           shows if there are any integrals from this small block
c           if nintez=0 no integrals
c ngctoz  - total contraction (general) length (ngcd)
c
c moreint - shows if all small blocks from a given super-block
c           have been processed or not
c stopnow - tells if the requested super-block exists or not
c----------------------------------------------------------------
      common /interface/ibuffx,nblsix,nintex,ngctox
c----------------------------------------------------------------
c  Values of where :
c
c     where.eq.'disk' put integrals on a disk
c     where.eq.'fock' put integrals to the Fock matrix
c     where.eq.'core' put integrals to the core storage
c
c     where.eq.'shif' NMR/GIAO two-el.int.der. to the Fock(D0,G1)
c     where.eq.'forc' for gradient integralS (first geom. derivatives)
c     where.eq.'hess' for hessian integrals (second geom. derivatives)
c--
c     where.eq.'sch1' Schwarz integrals
c     where.eq.'lmp2' integrals for LMP2 (for given ICS,KCS)
c
c----------------------------------------------------------------
      common /schwarz/ schw4    ! used only in this file & in precalc.f
      common /runtype/ scftype,where
      common /memor1/ iisd,jjsd,ijbld
      common /memor1b/ nbl2,nbloks
c
      dimension bl(*)
      dimension inx(12,*)
      dimension denspar(*)
      dimension labels(*)
c----------------------------------------------------------------
      schw4=wherex
c----------------------------------------------------------------
      where=wherex
c----------------------------------------------------------------
ctest
c     write(6,*)' from calcint_new : block=',isupb_r,' where=',where
c     call f_lush(6)
c----------------------------------------------------------------
CX    if(where.NE.'spec') then
CX       write(6,*)' can not do this task=',where
CX       stop ' can not do this task'
CX    endif
c----------------------------------------------------------------
      if(where.eq.'sch1') where='fock'
      if(where.eq.'sch2') where='shif'
      if(where.eq.'sch3') where='forc'
      if(where.eq.'sch4') where='hess'
      if(where.eq.'lmp2') where='fock'
      if(where.eq.'spec') where='fock'
c----------------------------------------------------------------
      nintex=0
c----------------------------------------------------------------
      stopnow=.false.
      if(isupb_r.le.0) then
         call nerror(1,'calcint2',
     *                 ' negative block number; can not do it',
     *                   isupb_r,0)
      endif
c     if(isupb_r.gt.nbloks) then
c        write(6,*)' stopnow=T, isupb_r=',isupb_r,' nbloks=',nbloks
c        stopnow=.true.
c        return
c     endif
c----------------------------------------------------------------
c Set up the integral threshold :
c
      call setup_thres(thres)
c----------------------------------------------------------------
c calculations of two-electron integrals in blocks
c Nbl2   - number of blocks of pairs
c
        call blockint1_n(bl,nbl2,inx,bl(iisd),bl(jjsd),bl(ijbld),
     *                   denspar, labels, isupb_r, moreint)
c
      ibuffz= ibuffx
      nblsiz= nblsix
      nintez= nintex
      ngctoz= ngctox
c----------------------------------------------------------------
      end
c===================================================================
      subroutine blockint1_n(bl,nbl2,inx,iis,jjs,ijbl,
     *                       denspar, labels, isupb_r, moreint)
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint
      character*11 scftype
      character*4 where
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file & in precalc.f
c----------------------------------------------------------------
      common /superbl/ isupb,ibl,kbl
      common /ilepar/ lpartot,lpareal
      common /runtype/ scftype,where
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /lindvec/ lind,idensp
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension denspar(*)
      dimension labels(*)
c----------------------------------------------------------------
      data ncall_2 /0/
c----------------------------------------------------------------
      save ncall_2, isupb_last, ibl_last, kbl_last, last0
c     save iblok_r
c----------------------------------------------------------------
c for check-point :
      if(.not.moreint) then
         call getmem(1,last0)
         call retmem(1)
      endif
c----------------------------------------------------------------
      if(ncall_2.eq.0) then
         isupb_last=0
         ibl_last  =0
         kbl_last  =0
      endif
c
      if(isupb_r.lt.isupb_last) then
c        new task (e.g. new scf cycle or...)
         isupb_last=0
         ibl_last  =0
         kbl_last  =0
      endif
c
      ncall_2=ncall_2+1
c------------------------------------------------------
c calculate integrals from one requested super-block isupb_r :
c
c1999
      nblok_r=1
      ikbl_r =1
c1999
c-----------------------------------
c find which pair-blokcs make this super-block :
c
      call getival('nblock4',nblock4)
      call get_pair_blocks(bl(nblock4),isupb_r,ibl_r,kbl_r)
c----------------------------------------------------------------------
      IF(   where.eq.'fock'                    .or.
     *      where.eq.'shif'                    .or.
     *      where.eq.'forc'                    .or.
     *      where.eq.'hess'                  ) THEN
c----------------------------------------------------------------------
       call prec2ij_new(ibl_r,bl)
       call prec2kl_new(kbl_r,bl)
c----------------------------------------------------------------------
c CALCULATE INTEGRALS :
c
c          write(6,*)' B_B=',isupb_r,' p-bl=',ibl_r,kbl_r
c          call f_lush(6)
      call onesuper1_n(isupb_r,ibl_r,kbl_r,bl,inx,
     *               nbl2, iis,jjs,ijbl,where,
     *               denspar,labels)
c----------------------------------------------------------------------
c release memory reserved in prec2ij WHEN a given super-block is done
c
               call retmem(1)
            call retmem(1)
c----------------------------------------------------------------------
c remember what was calculated :
      isupb_last=isupb_r
        ibl_last=ibl_r
        kbl_last=kbl_r
c----------------------------------------------------------------------
c tell a user if there is anymore integrals
c to be calculated from this super-block :
c
         moreint=.false.
c----------------------------------------------------------------------
      ENDIF
c----------------------------------------------------------------------
c check-point:
      if(.not.moreint) then
        call getmem(1,last1)
        call retmem(1)
        if(last1.ne.last0) then
           call nerror(4,'blocint1',
     *                   ' WRONG memory allocations ',last0,last1)
        endif
      endif
c-------------------------------------------------------------------
      end
c===================================================================
      subroutine onesuper1_n(isupb,ibl,kbl, bl,inx,
     *                     nbl2, iis,jjs,ijbl,where,
     *                     denspar,labels)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*4 where
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file & in precalc.f
      common /ijcsfl/ ijblokp,ijprevf,ijprevl,ijtprev,maxprev,ngcprev
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension denspar(*)
      dimension labels(*)
c----------------------------------------------------------------
c get needed parameters from depository :
c
      call getival('inx_1' ,inx_1 )
      call getival('inx_2' ,inx_2 )
c
      call getival('nparx',nparx)
      call getival('ijblx',ijblx)
c     call getival('mapijblx',map_ij_bl2)
c     call getival('blocksij',nbl2_ij)
      call getival('blpredij',nbl2_ijd)
c
c get needed parameters from depository :
c
      call getival('inx_3' ,inx_3 )
      call getival('inx_4' ,inx_4 )
c
      call getival('npary',npary)
      call getival('ijbly',ijbly)
c     call getival('mapijbly',map_kl_bl2)
c     call getival('blockskl',nbl2_kl)
      call getival('blpredkl',nbl2_kld)
c----------------------------------------------------------------
      call mmark
c----------------------------------------------------------------
c Constract blocks of contracted shell qyuartets for a given
c super-block (isupb) ; set up common /memor2/ and /memor3/ .
c
      call blockin4_n(isupb,ibl,kbl,bl)
c----------------------------------------------------------------
c loop over blocks belonging to the super-block isupb:
c NO loop : no small blocks
c
      ijblokp=0
      ijprevf=0
      ijprevl=0
      ijtprev=0
      maxprev=0
      ngcprev=0
c
      ikbl=1
c
         call oneblock_n(isupb,ikbl,bl, bl(nibld),bl(nkbld),
     *                   bl(nijbd),bl(nijed),bl(nklbd),bl(nkled),
     *                   bl(nqrtd),
     *                   iis,jjs,ijbl,nbl2,inx,
     *                   bl(ijblx),nbl2_ijd,bl(inx_1),bl(inx_2),
     *                   bl(ijbly),nbl2_kld,bl(inx_3),bl(inx_4),
     *                   where,
     *                   denspar,labels)
c----------------------------------------------------------------
      call retmark
c----------------------------------------------------------------
      end
c===================================================================
      subroutine oneblock_n(isupb,ikbl,bl,
     *                      nibl,nkbl,nijb,nije,nklb,nkle,
     *                      nqrt,
     *                      iis,jjs,ijbl,nbl2,inx,
     *                      ijbl_12 ,nbl2_ijd,   inx_1,   inx_2,
     *                      ijbl_34 ,nbl2_kld,   inx_3,   inx_4,
     *                      where,
     *                      denspar,labels)
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      logical stable
c
      character*4 where
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file & in precalc.f
      common /interface/ibuffx,nblsix,nintex,ngctox
      common /route/ iroute
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
cxx
      common /logic2/ len(1)
      common /logic4/ nfu(1)
cxx
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
      common /memor5g/ indxr
c
c nmr derivatives :
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
      common /symmet/ rnsym
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension nqrt(*)
c
      dimension ijbl_12(nbl2_ijd,*),inx_1(12,*),inx_2(12,*)
      dimension ijbl_34(nbl2_kld,*),inx_3(12,*),inx_4(12,*)
c
      dimension denspar(*)
      dimension labels(*)
cpp
c------------------------------------------------------
c  number of ij and kl pairs (npij,npkl)
c  (zero for npkl means that the block is diagonal)
      ibl=nibl(ikbl)
      kbl=nkbl(ikbl)
c
      nijbeg=nijb(ikbl)
      nijend=nije(ikbl)
      npij=nijend-nijbeg+1
c
      nklbeg=nklb(ikbl)
      nklend=nkle(ikbl)
      npkl=nklend-nklbeg+1
      if(nklend.eq.0) npkl=0
c1999
      nbls=nqrt(ikbl)
c1999
c
      if(nbls.eq.0) then
         write(6,*)' nbls=0 after blockin4 : go to 100 executed'
         nintex=0
         go to 100
      endif
c-------------------------------------------------------------
c for Schwarz integrals only diagonal blocks :
c
      if(schw4.eq.'sch1') then
         if(npkl.ne.0) then
            nintex=0
            go to 100
         endif
      endif
c-------------------------------------------------------------
c1999    The use of nblok1(2,nbls) array eleminated !!!!
c-------------------------------------------------------------
         ntotal=ntotal+nbls
c
c        make all quartets present :
c
         call isymm1(nbls,bl(isymm))
c
c-------------------------------------------------------------
c  first quartet of shells :
c
c     ijcs1=nblok1(1,1)
c     klcs1=nblok1(2,1)
ckw1999
      ijcs1=ijbl_12(ibl,nijbeg)
      klcs1=ijbl_34(kbl,nklbeg)
c
c     ijcs1=ijbl(ibl,nijbeg)
c     klcs1=ijbl(kbl,nklbeg)
c-------------------------------------------------------------
c  set up the type (common types) and length of shells
c  and information concerning contractions (commons contr,gcont) :
c
      call shells(inx,iis,jjs,ijcs1,klcs1)
c-------------------------------------------------------------
c setup the obarai and shell commons :
c
      call iobara(itype1,jtype1,ktype1,ltype1,where)
c
c-------------------------------------------------------------
      nfumax=nfu(mmax)
      mmax1=mmax-1
      if(mmax1.eq.0) mmax1=1
c
      npklx=npkl
      if(npkl.eq.0) npklx=npij
c
c98 it is needed only for Tracy's recursive and since this
c   recursive was re-formulated it is now always as follows :
c
      nhabcd=npklx
      nfha=nfumax*nhabcd*lckl
c-------------------------------------------------------------
c memory handling
c
c 1: for pairs precalculations (20 quantities)
c (for whole block of contracted quartets and
c  all primitive quartets belonging to this block)
c
c 2: and quartets precalculations (12 quantities)
c (for whole block of contracted quartets and
c        one primitive quartet )
c
c 3: and for individual shells (8 quantities)
c ( aa, bb, cc, dd - exponents )
c ( cis,cjs,cks,cls -contr coef.)
c
      IF( iroute.eq.1 ) THEN
         call memo5a_1(npij ,mmax1)
         call memo5b_1(npklx,mmax1)
         call memo5c_1(nbls,mmax1,npij,npklx,nfha,nfumax)
c
c        51 or 57 calls of getmem
         ncalls=51
         if(ngcd.gt.1) ncalls=57
      ELSE
         call memo5a_2(npij ,mmax1)
         call memo5b_2(npklx,mmax1)
         call memo5c_2(nbls,mmax1,npij,npklx,nfumax)
c
c        48 or 52 calls of getmem
         ncalls=48
         if(ngcd.gt.1) ncalls=52
      ENDIF
c
c-------------------------------------------------------------
c Perform pre-calculations for pairs of contracted sells
c
c-------------------------------------------------------------
      IF( iroute.eq.1 ) THEN
         call precalc2_1_new(bl,mmax,mmax1,nhabcd,nfumax,
     *                       ibl,nijbeg,nijend,npij,
     *                       kbl,nklbeg,nklend,npklx,npkl)
      ELSE
         call precalc2_2_new(bl,mmax,mmax1, nfumax,
     *                       ibl,nijbeg,nijend,npij,
     *                       kbl,nklbeg,nklend,npkl)
      ENDIF
c
c-------------------------------------------------------------
c Use density matrix to neglect integrals :
c
      call getmem(nbls,idnsx)
      ncalls=ncalls+1
c-------------------------------------------------------------
c for NMR derivatives
c
      if(where.eq.'shif') then
         call memo6(npij,npklx)        !   2 calls of getmem
         call getmem(npij ,ijcent)
         call getmem(npklx,klcent)
         call precal2d(bl(inuc),iis,jjs,inx, npij,npklx,npkl,
     *           ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *           bl(ixyab),bl(ixycd),
     *           bl(isymm),bl(ijcent),bl(klcent) )
         call retmem(2)
      endif
c-------------------------------------------------------------
c for GRADIENT derivatives (eliminate one-center quartets):
c
      if(where.eq.'forc' .or. where.eq.'hess') then
         if(schw4.eq.'sch1') then
         else
            call getmem(npij ,ijcent)
            call getmem(npklx,klcent)
            call onecentX(inx_1,inx_2,npij,
     *                    inx_3,inx_4,npklx,npkl,
     *              ibl,kbl,
     *              ijbl_12,nbl2_ijd,ijbl_34,nbl2_kld,
     *              nijbeg,nijend,nklbeg,nklend,
     *              bl(isymm),bl(ijcent),bl(klcent) )
            call retmem(2)
         endif
      endif
c-------------------------------------------------------------
c for Schwarz integrals (eliminate off-diagonal quartets):
c
      if(schw4.eq.'sch1') then
         if(npkl.ne.0)
     $     call nerror(1,'oneblock_n',' off-diagonal block',0,0)
         if(kbl.ne.ibl)
     $     call nerror(1,'oneblock_n',' off-diagonal block',0,0)
         call schw_quarts(iis,jjs,ibl,ijbl,nbl2,nijbeg,nijend,
     *                    bl(isymm) )
         call make_dens1(bl(idnsx),nbls)
      endif
c-------------------------------------------------------------
c neglect some contracted quartets using Schwarz inequality :
c isym(ijkl)=0 or 1 , 0 means that ijkl quartet is neglected.
c
      if(schw4.eq.'sch1') then
      else
         call schw_neg_new(nbls,ncs,iis,jjs,
     *                     ibl,ijbl_12,nbl2_ijd,nijbeg,npij,
     *                     kbl,ijbl_34,nbl2_kld,nklbeg,npklx,npkl,
     *                     denspar,
c output
     *                     bl(isymm), bl(idnsx))
c
c symmetry handling :
c
         if(schw4.NE.'lmp2') then
            if(nsym.gt.0) then
              call onesym(ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *                     bl(ijshp),nsym,bl(isymm),nblsym )
               noijkl=noijkl+nblsym
            else
               noijkl=noijkl+nbls
            endif
         endif
      endif
c-------------------------------------------------------------
c At this moment some of c.s.quartets may not be present due
c to the symmetry or due to neglect according to the CSHNEG .
c Thus, reduce the block size from  nbls ---> nblsp (present)
c (only present c.s.quartets are considered)
c                  and
c   re-define idx1 and idx2 (from index4)
c   and constructe indx(ijkl)=ijklp for use in prec4neg & precspc
c
      call indexp(npij,npklx,npkl,nijbeg,nijend, nklbeg,nklend,
     *            bl(idx1),bl(idx2),bl(isymm),
     *            bl(indxp),bl(indxr),nblsp)
c
c-------------------------------------------------------------
c FROM hereafter the block-size is reduced nbls-->nblsp
c
      if(nblsp.eq.0) then
        nintex=0
        ncalls=ncalls-1
c
c       ncalls-1 comes from the fact that no memory
c       has been allocated for buf in this case :
c
c       write(6,*)' from calcint_new nblsp=0 isupb=',isupb
        go to 110
      endif
c
      nopres=nopres+nblsp
      nblsrem=nbls
      nbls=nblsp
c
c-------------------------------------------------------------
c get max value of max density and max values of estcd for
c each contraction kl (will be uesd in prec4neg & precspec)
c
      call getmem(lcij,iesti2ij)
      call getmem(lckl,iesti2kl)
      call setival('esti2ij',iesti2ij)
      call setival('esti2kl',iesti2kl)
      ncalls=ncalls+2
c
      call get_max4est2(nblsrem,bl(idnsx),npkl,
     *                  npij ,lcij,bl(ieab),bl(iesti2ij),
     *                  npklx,lckl,bl(iecd),bl(iesti2kl))
c
c for use in prec4neg_ & precspec_
c     call setival('esti2ij',iesti2ij)
c     call setival('esti2kl',iesti2kl)
c-------------------------------------------------------------
c find cases when for a given ijpair all klpairs
c yield quartests which are NOT present by symmetry
c This is only when symmetry is present
c
      isymfac=1
      jump=1
      if(nsym.gt.0) then
         call getint(npij,jump)
         ncalls=ncalls+1
         call symm_jump(bl(isymm),npij,npklx,npkl,bl(jump))
c
         call getmem(npij*npklx,isymfac)
         ncalls=ncalls+1
         call symm_fact(bl(isymm),npij,npklx,npkl,rnsym,bl(isymfac))
      endif
c
c save these addresses as parameters
c
      call setival('jump',jump)
      call setival('isymfac',isymfac)
c
c jum(ijpar) and symfac(ijkl) will be used in prec4neg & precspec .
c-------------------------------------------------------------
c calculate special integrals with mmax<=2
c (ssss), (psss),(lsss) etc.
c
      if(mmax.le.2) then
c        allocate memory for list_ij, list_kl (pairs reordering)
c        used in eriteg_ & erintsp_
c
         call getint(npij ,list_ij)
         call getint(npklx,list_kl)
c
         call setival('list_ij',list_ij)
         call setival('list_kl',list_kl)
c
         IF( iroute.eq.1 ) THEN
            call erintsp_1_n(
     *           bl,first,nbls,acc,ikbl,npij,npklx,npkl,idnsx)
         ELSE
            call erintsp_2_n(
     *           bl,first,nbls,acc,ikbl,npij,npklx,npkl,idnsx)
         ENDIF
         call retmem(2)   ! release mem. allocated for list_ij,list_kl
         if( first ) then
             nintex=0
             go to 110
         endif
         go to 120
      endif
c
c-------------------------------------------------------------
c calculate integrals with MMAX > 2
c
      if(nsij.ge.nskl) then
         immax=mmax-2
         kmmax=nskl-2
         lobsa=2
         if(lshelij.gt.0) lobsa=1
      else
         immax=nsij-2
         kmmax=mmax-2
         lobsa=4
         if(lshelkl.gt.0) lobsa=3
      endif
c
c        allocate memory for list_ij, list_kl (pairs reordering)
c        used in eriteg_ & erintsp_
c
         call getint(npij ,list_ij)
         call getint(npklx,list_kl)
c
         call setival('list_ij',list_ij)
         call setival('list_kl',list_kl)
c
      call check_stab(stable,where)
c
      IF( iroute.eq.1 ) THEN
         call erinteg_1_n(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                  lobsa,immax,kmmax, where,idnsx,stable)
      ELSE
         call erinteg_2_n(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                  lobsa,immax,kmmax, where,idnsx,stable)
      ENDIF
c
         call retmem(2)   ! release mem. allocated for list_ij,list_kl
c-------------------------------------------------------------
      if( first ) then
          nintex=0
          go to 110
      endif
c-------------------------------------------------------------
c transformatin d6->d5 , f10->f7
c
      lenall=ilen+jlen+klen+llen
      lenall1=ilen1+jlen1+klen1+llen1
c
      if(lenall1-lenall.ne.0 ) then
c        call second(trans1)
c
c
c......
         mnbls=nbls
         nbuf=ibuf
         if(where.eq.'shif') then
c        --- nmr derivatives ----
           mnbls=6*nbls
           nbuf=ibuf+ngcd*nbls*lnijkl
         endif
         if(where.eq.'forc' .or. where.eq.'hess') then
c        --- 1st-order derivatives for gradient & hessian
           mnbls=9*nbls
           nbuf=ibuf
         endif
c......
c
         incrt=mnbls*lnijkl
c
c gen.cont.
         do 130 iqu=1,ngcd
           jbuf=nbuf+(iqu-1)*incrt
           call transfor(bl,mnbls,jbuf,ibuf2,itype,jtype,ktype,ltype,
     *                   ilen,jlen,klen,llen,ilen1,jlen1,klen1,llen1)
  130    continue
c
         if(where.eq.'hess') then
c        --- 2ed-order derivatives for hessian
           mnbls=45*nbls   ! for secend and first derivatives together
           nbuf=ibuf+9*ngcd*nbls*lnijkl
c
           incrt=mnbls*lnijkl
           do 131 iqu=1,ngcd
           jbuf=nbuf+(iqu-1)*incrt
           call transfor(bl,mnbls,jbuf,ibuf2,itype,jtype,ktype,ltype,
     *                   ilen,jlen,klen,llen,ilen1,jlen1,klen1,llen1)
  131      continue
         endif
      endif
c
c end of transformation
c-------------------------------------------------------------
c
  120 continue
c-------------------------------------------------------------
c At this moment integrals are in the buffer:
c
c (1) if(whwre.eq.'fock'or'disk'or'core' in bl(ibuf)
c (2) if(whwre.eq.'shif'                 in bl(ibuf+nbls*lnijkl*ngcd)
c
c Corresponding labels arrays will be constracted :
c Then, integrals and labels will be returned to the calling
c place and the USER is resposible to put them somewhere :
c Thus, regardless of the value of where, construct labels :
c
        length=1+4*nbls*ngcd
        lgenct=length+4
c
        call destret_new(ibl,ijbl_12,nbl2_ijd,inx_1,inx_2,
     *                   kbl,ijbl_34,nbl2_kld,inx_3,inx_4,
     *                   nijbeg,nijend,nklbeg,nklend,
     *               bl(isymm),
     *               nbls,iis,jjs,
     *               bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,
     *       bl(indxr), labels(1),labels(length),labels(lgenct) )
c
c       setup common /interface/
c       integ. address in bl :
        ibuffx=ibuf
        if(where.eq.'shif') ibuffx=ibuf+nbls*lnijkl*ngcd
c sizes :
        nblsix=nbls
        nintex=lnijkl
        ngctox=ngcd
c
c-------------------------------------------------------------
  110 continue
c-------------------------------------------------------------
c release memory at the end of a given block :
c
c     ncalls=41
c     if(ngcd.gt.1) ncalls=45
c
c nmr deriv:
      if(where.eq.'shif') ncalls=ncalls+2
change  call retmem(ncalls)
        call retmem(ncalls+1)   ! +1 to release ibuf
c-------------------------------------------------------------
  100 continue     ! jump here if 0 quartets by symmetry
c-------------------------------------------------------------
      end
c===================================================================
      subroutine destret_new(ibl,ijbl_12,nbl2_ijd, inx_1,inx_2,
     *                       kbl,ijbl_34,nbl2_kld, inx_3,inx_4,
     *                       nijbeg,nijend,nklbeg,nklend,
     *                   ipres,
     *                   nbls,iis,jjs,
     *                   icfg,jcfg,kcfg,lcfg,ngcd,indxr,
     *                   labels,length,lgenct )
c
c----------------------------------------------------------------
c This routine is called from calcint2 for each block
c (if it has any intregrals)
c----------------------------------------------------------------
c This is for returning itegrals to the calling place .
c This routine constructs only label's arrays
c----------------------------------------------------------------
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
c
      dimension ipres(*)
      dimension iis(*),jjs(*)
      dimension ijbl_12(nbl2_ijd,*),inx_1(12,*),inx_2(12,*)
      dimension ijbl_34(nbl2_kld,*),inx_3(12,*),inx_4(12,*)
c
      dimension iix(4)
c
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension indxr(*)
c------------------------
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      length(1)=ilen
      length(2)=jlen
      length(3)=klen
      length(4)=llen
c----------------------------------------------------------------
ctest
c     write(8,*)' FROM DESTRET pair-blocks no=',ibl,kbl
c----------------------------------------------------------------
c
c  loop over quartets belonging to the block :
c
      ijkl =0
      do 100 ijp=nijbeg,nijend
         ijcs=ijbl_12(ibl,ijp)
         nklendx=nklend
         if(nklend.eq.0) nklendx=ijp
         do 100 klp=nklbeg,nklendx
            klcs=ijbl_34(kbl,klp)
            ijkl=ijkl+1
            if(ipres(ijkl).eq.0) go to 100
c
            ijklp=indxr(ijkl)
c                                        In the quartet IJKL :
c                                        --------------------
            ics=iis(ijcs)               ! contracted shell ICS
            jcs=jjs(ijcs)               !     -"-          JCS
            kcs=iis(klcs)               !     -"-          KCS
            lcs=jjs(klcs)               !     -"-          LCS
c
            icff1=inx_1(11,ics)
            jcff1=inx_2(11,jcs)
            kcff1=inx_3(11,kcs)
            lcff1=inx_4(11,lcs)
c
            if(ngcd.eq.1) then
               ngcq=1
               icfg(1)=icff1
               jcfg(1)=jcff1
               kcfg(1)=kcff1
               lcfg(1)=lcff1
            else
               call indexg_new(icff1,jcff1,kcff1,lcff1,
     *                         ilen,jlen,klen,llen,
     *                         icfg,jcfg,kcfg,lcfg,ngcq)
            endif
c
            lgenct(ijklp)=ngcq
c
            do 150 iqu=1,ngcq
               icff=icfg(iqu)
               jcff=jcfg(iqu)
               kcff=kcfg(iqu)
               lcff=lcfg(iqu)
c
               labels(1,iqu,ijklp)=icff
               labels(2,iqu,ijklp)=jcff
               labels(3,iqu,ijklp)=kcff
               labels(4,iqu,ijklp)=lcff
  150       continue
  100 continue
c
      end
c==============================================================
      subroutine indexg_new(icff,jcff,kcff,lcff,
     *                      ilen,jlen,klen,llen,
     *                      icfg,jcfg,kcfg,lcfg,ngcq)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension iix(100),jjx(100),kkx(100),llx(100)
c------------------------------------------------------------
c dim. 100 should be enough since the max. ge.con is 9, so 81
c is actually the max. for iix,jjx,kkx,llx
c------------------------------------------------------------
c
             ijpg=0
             icf=icff
             do 2041 igc=0,ngci
               jcf=jcff
               do 2042 jgc=0,ngcj
                  ijpg=ijpg+1
                  iix(ijpg)=icf
                  jjx(ijpg)=jcf
                  jcf=jcf+jlen
 2042          continue
               icf=icf+ilen
 2041        continue
c
             klpg=0
             kcf=kcff
             do 2043 kgc=0,ngck
               lcf=lcff
               do 2044 lgc=0,ngcl
                  klpg=klpg+1
                  kkx(klpg)=kcf
                  llx(klpg)=lcf
                  lcf=lcf+llen
 2044          continue
                kcf=kcf+klen
 2043        continue
c
             ijklg=0
             do 2045 ijp1=1,ijpg
             do 2045 klp1=1,klpg
                ijklg=ijklg+1
                icfg(ijklg)=iix(ijp1)
                jcfg(ijklg)=jjx(ijp1)
                kcfg(ijklg)=kkx(klp1)
                lcfg(ijklg)=llx(klp1)
 2045        continue
c
      ngcq=ijklg
c
      end
c===================================================================
c the only difference compare to original ones is that here
c we call prec4neg_2_n instead prec4neg_2 .
c we call precspec_2_n instead precspec_2 .
c===================================================================
      subroutine erinteg_2_n(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                     lobsa,immax,kmmax, where,idnsx,stable)
c
c character variable WHERE is needed for derivatives only
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      logical stable
      character*4 where
      common /route/ iroute
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
ctime
      common /timex/ tconv1,tconv2,ttrobs
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time3/ tzeroi,tspeci
      common /time4/ tderiv
ctime
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
      common /logic4/ nfu(1)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
      common /memor5g/ indxr
c nmr der.
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
      common /symmet/ rnsym
c
      dimension bl(*)
c-----------------------------------------------
      call getival('jump',jump)
      call getival('isymfac',isymfac)
c-----------------------------------------------
      call getival('esti2ij',iesti2ij)
      call getival('esti2kl',iesti2kl)
c
      call getival('list_ij',list_ij)
      call getival('list_kl',list_kl)
c-----------------------------------------------
c memory handling for a block
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(nbls, l11,l12,mem2,igmcnt)
c
c97....
      nfumax=nfu(mmax)
      mmax1=mmax-1
      nfumax1=nfumax
      narray=1
      nqijr=nqij
      nqklr=nqkl
      if(where.eq.'shif') nfumax1=nfu(mmax1)
      if(where.eq.'forc') then
         nfumax1=nfu(mmax1)
         narray=4
         nqij=nqij-1
         nqkl=nqkl-1
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
      if(where.eq.'hess') then
         nfumax1=nfu(mmax1-1)
         narray=10
         nqij=nqij-2
         nqkl=nqkl-2
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
c
c****
      first=.true.
c****
      ibut2=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcijx)
        call getmem(ngckl*nbls,igcklx)
        nblsx=nbls*narray
        call getmem(ngcd*nblsx*lnij*lnkl,ibut2)
      endif
c
c loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
c2002
      if(.not.stable) then
         call set_boamax2(bl(iaa),bl(ibb),l1,l2,boa_max) ! used in xwpq_22
         call setrval('boa_max',boa_max)
      endif
c2002
      lcij=lcij+1
      if(ngcd.gt.1) then
        call gcpairs(nbls,lcij,ngci1,ngcj1,ngcij,
     *               bl(igcij),bl(igcijx))
      endif
         lckl=0
         nblsr=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         lckl=lckl+1
c2002
         if(.not.stable) then
            call set_boamax2(bl(icc),bl(idd),l3,l4,doc_max) ! used in xwpq_22
            call setrval('doc_max',doc_max)
         endif
c2002
         if(ngcd.gt.1) then
           call gcpairs(nbls,lckl,ngck1,ngcl1,ngckl,
     *                  bl(igckl),bl(igcklx))
      endif
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
c      call second(t4b)
c
c2002 : ok for MP2 but not for MP2 gardient
c neglect some primitive blocks & re-order pairs in these which remain
c
c
c     call neg_primitive(bl,nbls,npij,npklx,npkl,lcij,lckl,lc12,lc34,
c    *    bl(iapb),bl(icpd),bl(ieab),bl(iecd),bl(iesti2ij),bl(iesti2kl),
c    *    bl(list_ij),bl(list_kl),nbls1)
c
c     if(nbls1.eq.0) go to 34
c
c2002 : ok for MP2 but not for MP2 gardient
c
      call prec4neg_2_n(nbls,npij,npklx,npkl,lcij,lckl,
     1     lc12,lc34, bl(idx1),bl(idx2),bl(indxr),
     *     bl(ieab),bl(iecd),bl(idnsx),bl(iesti2ij),bl(iesti2kl),
     *     bl(list_ij),bl(list_kl),
     *     bl(jump),bl(isymfac),
     *     bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq), nsym,bl(isymm),
     * bl(irho),bl(irppq),bl(irhoapb),bl(irhocpd),bl(irys),bl(iconst),
     *                          nbls1,bl(indx) )
c
c
c      call second(t4e)
c      tpre4n=tpre4n+(t4e-t4b)
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c------------------------------------------------------
c note :
c from here a given block (nbls) was reduced to (nbls1) ***
c              for given l1,l2,l3,l4 !!!
c                      and
c rho, rppq, rhoapb, rhocpd, rys, and const
c have dimension nbls1 (not  (nbls) )
c------------------------------------------------------
c substitute zero for all needed bufors for quartets
c which do not appear first time after neglect
c
         if(first) then
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
               if(where.eq.'forc' .or. where.eq.'hess') then
c                 it is zero out in assemblx.f
               else
                  call zeroint(bl,nbls,nbls1,nblsnot,
     *                      lnij,lnkl,indx,ngcd,ibut2)
               endif
            endif
         endif
c
c------------------------------------------------------
c integrals of the type (i+j,s|k+l,s) :
c
c
c        call second(txwpqb)
c
        if(stable) then
         call xwpq_2(nbls1,bl(ixwp),bl(ixwq),bl(ip1234),
     1               ijpar1,lc12, klpar1,lc34,
     *               lcij,lckl,npij,npklx,
     *               bl(idx1),bl(idx2),bl(indx),
     *               bl(irppq),bl(ixp),bl(ixq),bl(ixpp),bl(ixqq),
     *               bl(itxab),bl(itxcd),bl(iabcd),
     *               bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb) )
        else
         call xwpq_22(nbls1,bl(ixwp),bl(ixwq),bl(ip1234),
     1               ijpar1,lc12, klpar1,lc34,
     *               lcij,lckl,npij,npklx,
     *               bl(idx1),bl(idx2),bl(indx),
     *               bl(irppq),bl(ixp),bl(ixq),bl(ixpp),bl(ixqq),
     *               bl(itxab),bl(itxcd),bl(iabcd),
     *               bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb) ,
     *                where,
     *               immax,kmmax,lobsa)
        endif
c------------------------------------------------------
c convert:  abnia,cdnia, xpn,xqn, habcd from pair
c to quartet (nbls1) quantities for trobsa
c
c ij
c
c        call second(tc1xb)
c        txwpq=txwpq+(tc1xb-txwpqb)
c
c
      if(nbls1.ne.nbls .or. nblsr.ne.nbls) then
         call conv1x_2(nbls1,mmax1,npij ,lcij, bl(idx1),bl(indx),
     *                 bl(ixpn), bl(ixpnx)  )
                       iabnix=iabnia+(lcij-1)*mmax1
      endif
c kl
         call conv1x_2(nbls1,mmax1,npklx,lckl, bl(idx2),bl(indx),
     *                 bl(ixqn), bl(ixqnx)  )
                       icdnix=icdnia+(lckl-1)*mmax1
c
                      ihabcdx=ihabcd+(lckl-1)*nfumax*3
c
c------------------------------------------------------
c
         call trobsa(bl,nbls1,l11,l12,mem2,immax,kmmax,lobsa)
c
c------------------------------------------------------
c
      if(ngcd.eq.1) then
         call assemblx(first,nbls,nbls1,lnij,lnkl,
     *                 l1,l2,l3,l4,lcij,lckl,npij,npklx)
      else
         call gcquart(nbls,nbls1, bl(indx),
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(igcijx),ngcij, bl(igcklx),ngckl,
     *                bl(indgc),bl(igcoef)  )
         call assemblg(first,nbls,nbls1,lnij,lnkl,
     *                 l1,l2,l3,l4,ngcd, ibut2)
      endif
c
c        call second(tasse)
c        tassem=tassem+(tasse-ttrob)
c
c------------------------------------------------------
c
         nblsr=nbls1
c
   34    continue
   12 continue
c
c------------------------------------------------------
c
c grad derivatives:
c
      nqij=nqijr
      nqkl=nqklr
c
c------------------------------------------------------
c transpose but2(ngcd, nbls*lnij*lnkl) into buf2(nbls*lnij*lnkl,ngcd)
c for Ist-derivatives put it into buf() not buf2()
c
      if(ngcd.gt.1) then
        lda=ngcd*narray
        ldb=nbls*lnij*lnkl
        into2=ibuf2
        if(where.eq.'forc' .or. where.eq.'hess') into2=ibuf
        call trspmo(bl(ibut2),lda, bl(into2),ldb)
        call retmem(3)
      endif
c------------------------------------------------------
      if( first ) then
          call retmem(igmcnt)
          return
      endif
c------------------------------------------------------
c
c release memory from trobsa
c
      call retmem(3)
c
c------------------------------------------------------
c for NMR, GRADIENT & HESSIAN derivatives :
c
      lnijr=lnij
      lnklr=lnkl
c------------------------------------------------------
c reserves memory for amshift
c
      call memo4b(nbls,igmcnt1)
c
c------------------------------------------------------
c for NMR/GIAO derivatives :
c
c
      if(where.eq.'shif') then
c        time for derivatives
c
c        call second(tderb)
c
c        return to the original values of nsij,nskl and mmax :
c
         ijderiv=1
         klderiv=1
         call iobarb(ijderiv,klderiv)
c
c
c construct (nmr) derivatives : (i+j,s|k+l,s)(x,y,z)
c
         call nmrderx(nbls,lnijr,lnklr,npij,npklx,ngcd,
     *                idx1,idx2, ixab,ixcd,ixyab,ixycd)
c
c
c        call second(tdere)
c        tderiv=tderiv+tdere-tderb
c
      endif
c------------------------------------------------------
c for GRADIENT derivatives :
c
      if(where.eq.'forc') then
c
c        call second(tderb)
c
c        return to the original values of nsij,nskl and mmax :
c
         ijderiv=1
         klderiv=1
         call iobarb(ijderiv,klderiv)
c
c construct gradient derivatives : (i+j,s|k+l,s)(x,y,z)
c
         call force_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c
c        call second(tdere)
c        tderiv=tderiv+tdere-tderb
c
      endif
c
c------------------------------------------------------
c for HESSIAN  derivatives :
c
      if(where.eq.'hess') then
c
c        call second(tderb)
c
c        return to the original values of nsij,nskl and mmax :
c
         ijderiv=2
         klderiv=2
         call iobarb(ijderiv,klderiv)
c
c construct hessian derivatives : (i+j,s|k+l,s)(x,y,z)
c
         call hessian_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c
c        call second(tdere)
c        tderiv=tderiv+tdere-tderb
      endif
c
c------------------------------------------------------
c
c  shift angular momentum (a->b, c->d) :
c
c     call second(tamsb)
c
      call amshift(nbls,lnij,lnkl,npij,npklx,ngcd)
c
c     call second(tamse)
c     tamshf=tamshf+(tamse-tamsb)
c
      call retmem(igmcnt+igmcnt1-3)
c
c------------------------------------------------------
      end
c===================================================================
      subroutine erintsp_2_n(bl,first,nbls,acc,ikbl,npij,npklx,npkl,
     *                       idnsx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /logic4/ nfu(1)
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
      common /memor5g/ indxr
      common /memors/ nsym,ijshp,isymm
      common /symmet/ rnsym
c
      dimension bl(*)
c-----------------------------------------------
      call getival('jump',jump)
      call getival('isymfac',isymfac)
c-----------------------------------------------
      call getival('esti2ij',iesti2ij)
      call getival('esti2kl',iesti2kl)
c
      call getival('list_ij',list_ij)
      call getival('list_kl',list_kl)
c-----------------------------------------------
c
c memory handling for a block
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(nbls, l11,l12,mem2,igmcnt)
c
c****
      first=.true.
c****
      ibut=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcijx)
        call getmem(ngckl*nbls,igcklx)
        call getmem(ngcd*nbls*lnijkl,ibut)
      endif
c
c loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
      lcij=lcij+1
      if(ngcd.gt.1) then
        call gcpairs(nbls,lcij,ngci1,ngcj1,ngcij,
     *               bl(igcij),bl(igcijx))
      endif
         lckl=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         lckl=lckl+1
         if(ngcd.gt.1) then
           call gcpairs(nbls,lckl,ngck1,ngcl1,ngckl,
     *                  bl(igckl),bl(igcklx))
      endif
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
c neglect some primitive blocks & re-order pairs in these which remain
c
c2002 : ok for MP2 but not for MP2 gardient
c
c     call neg_primitive(bl,nbls,npij,npklx,npkl,lcij,lckl,lc12,lc34,
c    *    bl(iapb),bl(icpd),bl(ieab),bl(iecd),bl(iesti2ij),bl(iesti2kl),
c    *    bl(list_ij),bl(list_kl),nbls1)
c
c     if(nbls1.eq.0) go to 34
c
c2002 : ok for MP2 but not for MP2 gardient
c
      call precspec_2_n(nbls,npij,npklx,npkl, lcij,lckl,
     1     lc12,lc34,bl(idx1),bl(idx2),bl(indxr),
     *     bl(ieab),bl(iecd),bl(idnsx),bl(iesti2ij),bl(iesti2kl),
     *     bl(list_ij),bl(list_kl),
     *     bl(jump),bl(isymfac),
     *     bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq),bl(i1apb),bl(i1cpd),bl(itxab),bl(itxcd),
     *     nsym,bl(isymm),
     *     bl(irho),bl(irys),bl(iconst),bl(ixwp),bl(ixwq),
     *                            nbls1,bl(indx) )
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c------------------------------------------------------
c  note :
c from here a given block (nbls) was reduced to (nbls1) ***
c              for given l1,l2,l3,l4 !!!
c                      and
c rho, rppq, rhoapb, rhocpd, rys, and const
c have dimension nbls1 (not  (nbls) )
c------------------------------------------------------
c substitute zero for all needed bufors for quartets
c which do not appear first time after neglect
c
         if(first) then
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
               call zeroint(bl,nbls,nbls1,nblsnot,
     *                      lnij,lnkl,indx,ngcd,ibut)
            endif
         endif
c
c------------------------------------------------------
c
c  special code for (ss ss) and (xs ss), (sx ss), (ss xs), (ss sx)
c  integrals ( x= p or l )
c
      if(ngcd.eq.1) then
          call specase_2(bl,first,nbls,nbls1, bl(indx),
     *          npij,npklx,l1,l2,l3,l4,
     *          bl(icis),bl(icjs),bl(icks),bl(icls),
     *          bl(ibuf),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),bl(irr1) )
      else
         call gcquart(nbls,nbls1, bl(indx),
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(igcijx),ngcij, bl(igcklx),ngckl,
     *                bl(indgc),bl(igcoef)  )
          call specasg(bl,first,nbls,nbls1, bl(indx),bl(idx1),bl(idx2),
chang*          bl(ibuf),bl(ibuf2),
     *          bl(ibut),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),
     *          ngcd,bl(indgc),bl(igcoef),lnijkl )
      endif
c
   34    continue
   12 continue
c
c------------------------------------------------------
c
c transpose but into buf
      if(ngcd.gt.1) then
        lda=ngcd
        ldb=nbls*lnijkl
cccccc  call dgetmo(bl(ibut),lda,lda,ldb,bl(ibuf),ldb)
        call trspmo(bl(ibut),lda,        bl(ibuf),ldb)
        call retmem(3)
      endif
c
      call retmem(igmcnt)
c
c------------------------------------------------------
c
      end
c===================================================================
      subroutine get_pair_blocks(nblock4,isupb_r,ibl_r,kbl_r)
      dimension nblock4(2,*)
c
      ibl=nblock4(1,isupb_r)
      kbl=nblock4(2,isupb_r)
      if(ibl.lt.0) ibl=-ibl
      if(kbl.lt.0) kbl=-kbl
c
      ibl_r=ibl
      kbl_r=kbl
c
      end
c===================================================================
c for iroute=1 :
      subroutine erinteg_1_n(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                       lobsa,immax,kmmax, where,idnsx,stable)
c
c character variable WHERE is needed for derivatives only
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      logical stable
      character*4 where
c
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
c
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
ctime
      common /timex/ tconv1,tconv2,ttrobs
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time3/ tzeroi,tspeci
      common /time4/ tderiv
ctime
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /logic4/ nfu(1)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
c new for grad. derivatives:
      common /memor5dd/ iaax,ibbx,iccx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
c nmr der.
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
      common /symmet/ rnsym
c
      dimension bl(*)
c-----------------------------------------------
      call getival('jump',jump)
      call getival('isymfac',isymfac)
c-----------------------------------------------
cxxxx
c
c  memory handling for a block
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(nbls, l11,l12,mem2,igmcnt)
c
      nfumax=nfu(mmax)
      mmax1=mmax-1
c97....
      nfumax=nfu(mmax)
      mmax1=mmax-1
      nfumax1=nfumax
      narray=1
      nqijr=nqij
      nqklr=nqkl
      if(where.eq.'shif') nfumax1=nfu(mmax1)
      if(where.eq.'forc') then
         nfumax1=nfu(mmax1)
         narray=4
         nqij=nqij-1
         nqkl=nqkl-1
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
      if(where.eq.'hess') then
         nfumax1=nfu(mmax1-1)
         narray=10
         nqij=nqij-2
         nqkl=nqkl-2
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
c97....
c
c****
      first=.true.
c****
      ibut2=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcij)
        call getmem(ngckl*nbls,igckl)
        nblsx=nbls*narray
        call getmem(ngcd*nblsx*lnij*lnkl,ibut2)
      endif
c
c97...
c
      if(where.eq.'forc' .or. where.eq.'hess') then
        call getmem(nbls,iaax)
        call getmem(nbls,ibbx)
        call getmem(nbls,iccx)
      endif
c
c97...
c
c loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
c200
      if(.not.stable) then
         call set_boamax1(bl(iaa),bl(ibb),npij,l1,l2,boa_max)
         call setrval('boa_max',boa_max)
      endif
c200
      if(ngcd.gt.1) then
        call gcparij(nbls, bl(idx1),npij,
     *               l1,l2,ngci1,ngcj1,ngcij,
     *               bl(igci),bl(igcj), bl(igcij))
      endif
      lcij=lcij+1
         lckl=0
         nblsr=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
c2002
         if(.not.stable) then
            call set_boamax1(bl(icc),bl(idd),npklx,l3,l4,doc_max)
            call setrval('doc_max',doc_max)
         endif
c2002
         if(ngcd.gt.1) then
           call gcparkl(nbls, bl(idx2),npklx,
     *                  l3,l4,ngck1,ngcl1,ngckl,
     *                  bl(igck),bl(igcl), bl(igckl))
         endif
         lckl=lckl+1
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
ctime
c      call second(t4b)
c
      call prec4neg_1_n(nbls,npij,npklx,npkl,lcij,lckl,
     1     ijpar1,lc12, klpar1,lc34, bl(idx1),bl(idx2),
     2     bl(ieab),bl(iecd),bl(idnsx),
     *     bl(jump),bl(isymfac),
     *     bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq), nsym,bl(isymm),
     *          bl(irppq),bl(irhoapb),bl(irhocpd),bl(irys),bl(iconst),
     *                          nbls1,bl(indx) )
ctime
c      call second(t4e)
c      tpre4n=tpre4n+(t4e-t4b)
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c  note :
c from here a given block (nbls) was reduced to (nbls1)
c              for given l1,l2,l3,l4 !!!
c                      and
c rho, rppq, rhoapb, rhocpd, rys, and const
c  have dimension nbls1 (not  (nbls) )
c------------------------------------------------------
c substitute zero for all needed bufors for quartets
c which do not appear first time after neglect
c
         if(first) then
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
               if(where.eq.'forc' .or. where.eq.'hess') then
c                 it is zero out in assemblx.f
               else
                  call zeroint(bl,nbls,nbls1,nblsnot,
     *                      lnij,lnkl,indx,ngcd,ibut2)
               endif
            endif
         endif
c
c------------------------------------------------------
c integrals of the type (i+j,s|k+l,s) :
c
        if(stable) then
         call xwpq_1(nbls1,bl(ixwp),bl(ixwq),bl(ip1234),
     1                ijpar1,lc12, klpar1,lc34,
     *                lcij,lckl,npij,npklx,
     *                bl(idx1),bl(idx2),bl(indx),
     *                bl(irppq),bl(ixp),bl(ixq),bl(ixpp),bl(ixqq),
     *                bl(itxab),bl(itxcd),bl(iabcd),
     *                bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb) )
        else
         call xwpq_11(nbls1,bl(ixwp),bl(ixwq),bl(ip1234),
     1                ijpar1,lc12, klpar1,lc34,
     *                lcij,lckl,npij,npklx,
     *                bl(idx1),bl(idx2),bl(indx),
     *                bl(irppq),bl(ixp),bl(ixq),bl(ixpp),bl(ixqq),
     *                bl(itxab),bl(itxcd),bl(iabcd),
     *                bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb) ,
     *                where,
     *                immax,kmmax,lobsa)
        endif
c------------------------------------------------------
c convert:  abnia,cdnia, xpn,xqn, habcd from pair
c to quartet (nbls1) quantities for trobsa
c
c ij
ctime
c        call second(tc1xb)
c        txwpq=txwpq+(tc1xb-txwpqb)
c
      if(nbls1.ne.nbls .or. nblsr.ne.nbls) then
         call conv1x_1(nbls1,mmax1,npij ,lcij, bl(idx1),bl(indx),
     *                 bl(iabnia), bl(ixpn), bl(iabnix),bl(ixpnx)  )
         if(where.eq.'forc' .or. where.eq.'hess') then
            call conv1der(nbls1,npij,l1,bl(idx1),bl(indx),
     *                    bl(iaa),bl(iaax) )
            call conv1der(nbls1,npij ,l2,bl(idx1),bl(indx),
     *                    bl(ibb),bl(ibbx) )
         endif
      endif
c kl
         call conv1x_1(nbls1,mmax1,npklx,lckl, bl(idx2),bl(indx),
     *                 bl(icdnia), bl(ixqn), bl(icdnix),bl(ixqnx)  )
         if(where.eq.'forc' .or. where.eq.'hess') then
            call conv1der(nbls1,npklx,l3,bl(idx2),bl(indx),
     *                    bl(icc),bl(iccx) )
         endif
c
         call conv2x(nbls1,nfumax1,npklx,lckl, bl(idx2),bl(indx),
     *               bl(ihabcd),nfumax, bl(ihabcdx)  )
c------------------------------------------------------
c
         call trobsa(bl,nbls1,l11,l12,mem2,immax,kmmax,lobsa)
c
c------------------------------------------------------
      if(ngcd.eq.1) then
         call assemblx(first,nbls,nbls1,lnij,lnkl,
     *                 l1,l2,l3,l4,lcij,lckl,npij,npklx)
      else
         call gcqijkl(nbls,nbls1, bl(indx), npij,npklx,
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(indgc),bl(igcoef),
     *                bl(igcij),ngcij, bl(igckl),ngckl)
         call assemblg(first,nbls,nbls1,lnij,lnkl,
     *                 l1,l2,l3,l4,ngcd, ibut2)
      endif
ctime
c        call second(tasse)
c        tassem=tassem+(tasse-ttrob)
c------------------------------------------------------
c
         nblsr=nbls1
c
   34    continue
   12 continue
c------------------------------------------------------
c grad  &  hess derivatives:
c
      nqij=nqijr
      nqkl=nqklr
c
c release memory allocated for 'forc' or 'hess':
c
      if(where.eq.'forc' .or. where.eq.'hess') call retmem(3)
c
c transpose but2(ngcd, nbls*lnij*lnkl) into buf2(nbls*lnij*lnkl,ngcd)
c for Ist-derivatives put it into buf() not buf2()
c
      if(ngcd.gt.1) then
        lda=ngcd*narray
        ldb=nbls*lnij*lnkl
        into2=ibuf2
        if(where.eq.'forc' .or. where.eq.'hess') into2=ibuf
        call trspmo(bl(ibut2),lda, bl(into2),ldb)
        call retmem(3)
      endif
c
      if( first ) then
          call retmem(igmcnt)
          return
      endif
c
c release memory from trobsa
c
      call retmem(3)
c
c-------------
c------------------------------------------------------
c for NMR, GRADIENT & HESSIAN derivatives :
c
      lnijr=lnij
      lnklr=lnkl
c------------------------------------------------------
c reserves memory for amshift
c
      call memo4b(nbls,igmcnt1)
c
c------------------------------------------------------
c for NMR/GIAO derivatives :
c
c
      if(where.eq.'shif') then
c
c        return to the original values of nsij,nskl and mmax :
c
         ijderiv=1
         klderiv=1
         call iobarb(ijderiv,klderiv)
c
c construct (nmr) derivatives : (i+j,s|k+l,s)(x,y,z)
c
         call nmrderx(nbls,lnijr,lnklr,npij,npklx,ngcd,
     *                idx1,idx2, ixab,ixcd,ixyab,ixycd)
c
      endif
c------------------------------------------------------
c for GRADIENT derivatives :
c
      if(where.eq.'forc') then
c
c        return to the original values of nsij,nskl and mmax :
c
         ijderiv=1
         klderiv=1
         call iobarb(ijderiv,klderiv)
c
c construct gradient derivatives : (i+j,s|k+l,s)(x,y,z)
c
         call force_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c
      endif
c------------------------------------
c for HESSIAN  derivatives :
c
      if(where.eq.'hess') then
c
c        return to the original values of nsij,nskl and mmax :
c
         ijderiv=2
         klderiv=2
         call iobarb(ijderiv,klderiv)
c
c construct hessian derivatives : (i+j,s|k+l,s)(x,y,z)
c
         call hessian_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c
      endif
c------------------------------------------------------
c
c  shift angular momentum (a->b, c->d) :
c
      call amshift(nbls,lnij,lnkl,npij,npklx,ngcd)
c
      call retmem(igmcnt+igmcnt1-3)
c
c------------------------------------------------------
      end
c===================================================================
      subroutine erintsp_1_n(bl,first,nbls,acc,ikbl,npij,npklx,npkl,
     *                       idnsx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c
      logical first
c
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
c
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /logic4/ nfu(1)
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
      common /memors/ nsym,ijshp,isymm
      common /symmet/ rnsym
c
      dimension bl(*)
c-----------------------------------------------
      call getival('jump',jump)
      call getival('isymfac',isymfac)
c-----------------------------------------------
c
c memory handling for a block
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(nbls, l11,l12,mem2,igmcnt)
c
c****
      first=.true.
c****
      ibut=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcij)
        call getmem(ngckl*nbls,igckl)
        call getmem(ngcd*nbls*lnijkl,ibut)
      endif
c
c loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
      if(ngcd.gt.1) then
        call gcparij(nbls, bl(idx1),npij,
     *               l1,l2,ngci1,ngcj1,ngcij,
     *               bl(igci),bl(igcj), bl(igcij))
      endif
      lcij=lcij+1
         lckl=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         if(ngcd.gt.1) then
           call gcparkl(nbls, bl(idx2),npklx,
     *                  l3,l4,ngck1,ngcl1,ngckl,
     *                  bl(igck),bl(igcl), bl(igckl))
      endif
         lckl=lckl+1
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
      call precspec_1_n(nbls,npij,npklx,npkl, lcij,lckl,
     1     ijpar1,lc12, klpar1,lc34,bl(idx1),bl(idx2),
     2     bl(ieab),bl(iecd), bl(idnsx),
     *     bl(jump),bl(isymfac),
     *     bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq),bl(i1apb),bl(i1cpd),bl(itxab),bl(itxcd),
     *     nsym,bl(isymm),
     *              bl(irys),bl(iconst),bl(ixwp),bl(ixwq),
     *                            nbls1,bl(indx) )
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c note :
c from here a given block (nbls) was reduced to (nbls1) ***
c              for given l1,l2,l3,l4 !!!
c                      and
c rho, rppq, rhoapb, rhocpd, rys, and const
c have dimension nbls1 (not  (nbls) )
c------------------------------------------------------
c substitute zero for all needed buffors for quartets
c which do not appear first time after neglect
c
         if(first) then
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
               call zeroint(bl,nbls,nbls1,nblsnot,
     *                      lnij,lnkl,indx,ngcd,ibut)
            endif
         endif
c
c------------------------------------------------------
c    special code for (ss ss) and (xs ss), (sx ss), (ss xs), (ss sx)
c          itegrals ( x= p or l )
c
      if(ngcd.eq.1) then
         call specase_1(bl,first,nbls,nbls1, bl(indx),bl(idx1),bl(idx2),
     *          npij,npklx,l1,l2,l3,l4,
     *          bl(icis),bl(icjs),bl(icks),bl(icls),
     *          bl(ibuf),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),bl(irr1) )
      else
         call gcqijkl(nbls,nbls1, bl(indx), npij,npklx,
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(indgc),bl(igcoef),
     *                bl(igcij),ngcij, bl(igckl),ngckl)
          call specasg(bl,first,nbls,nbls1, bl(indx),bl(idx1),bl(idx2),
chang*          bl(ibuf),bl(ibuf2),
     *          bl(ibut),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),
     *          ngcd,bl(indgc),bl(igcoef),lnijkl )
      endif
c
   34    continue
   12 continue
c------------------------------------------------------
c
c transpose but into buf
      if(ngcd.gt.1) then
         lda=ngcd
         ldb=nbls*lnijkl
cccccc   call dgetmo(bl(ibut),lda,lda,ldb,bl(ibuf),ldb)
         call trspmo(bl(ibut),lda,        bl(ibuf),ldb)
         call retmem(3)
      endif
c
      call retmem(igmcnt)
c
c------------------------------------------------------
      end
c===================================================================
      subroutine get_max4est2(nbls,densmax,ndiag,
     *                        npij,lcij,eab,esti2ij,
     *                        npkl,lckl,ecd,esti2kl)
      implicit real*8 (a-h,o-z)
      dimension densmax(nbls)
      dimension eab(npij,lcij)
      dimension ecd(npkl,lckl)
c output
      dimension esti2ij(lcij)
      dimension esti2kl(lckl)
c
      call getrval('dmx_over',ddmax)
c
      do ij=1,lcij
         esti2max=0.d0
         do ijp=1,npij
            esti2max=max(eab(ijp,ij),esti2max)
         enddo
         esti2ij(ij)=      esti2max
      enddo
c
      if(ndiag.ne.0) then
         do kl=1,lckl
            esti2max=0.d0
            do klp=1,npkl
               esti2max=max(ecd(klp,kl),esti2max)
            enddo
            esti2kl(kl)=ddmax*esti2max
         enddo
      else
         do kl=1,lckl
            esti2kl(kl)=ddmax*esti2ij(kl)
         enddo
      endif
c
      end
c====================================================================
      subroutine neg_primitive(bl,nbls,npij,npkl,ndiag,ij,kl,lc12,lc34,
     *                         apb,cpd,eab,ecd,esti2ij,esti2kl,
     *                         list_ij,list_kl,nbls1)
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr,eps8
      dimension bl(*)
      dimension apb(lc12), cpd(lc34)
      dimension eab(npij,lc12),ecd(npkl,lc34)
      dimension esti2ij(lc12), esti2kl(lc34)
      dimension list_ij(*),list_kl(*)
      data ij_last /0/
      save ij_last
c
      abpcd1=apb(ij)+cpd(kl)
      abpcdr=1.d0/abpcd1
      estimx=esti2ij(ij)*esti2kl(kl)*abpcdr
c
      if(estimx.lt.epsr) then
         nbls1=0
c......  ij_last=0
         if(kl.eq.lc34) ij_last=0
         return
      else
         nbls1=nbls
c        re-order pairs in accord with decreasing estimate value
         if(ij.ne.ij_last) then
            call reorder_primpairs(ij,npij,eab,lc12,list_ij)
         endif
c
         call reorder_primpairs(kl,npkl,ecd,lc34,list_kl)
         ij_last=ij
         if(kl.eq.lc34) ij_last=0
      endif
c
      end
c====================================================================
      subroutine reorder_primpairs(ij,npij,eab,lc12,list_ij)
      implicit real*8 (a-h,o-z)
      dimension eab(npij,lc12)
      dimension list_ij(*)
c
      do ijpar=1,npij
         list_ij(ijpar)=ijpar
      enddo
c
  100 continue
      iexch=0
      do ijpar=1,npij-1
         ijpar1=list_ij(ijpar)
         ijpar2=list_ij(ijpar+1)
         estij1=eab(ijpar1,ij)
         estij2=eab(ijpar2,ij)
         if(estij2.gt.estij1) then
            list_ij(ijpar  )=ijpar2
            list_ij(ijpar+1)=ijpar1
            iexch=1
         endif
      enddo
      if(iexch.eq.1) go to 100
c
      end
c====================================================================
