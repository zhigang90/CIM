c===================================================================
c OCT/NOV 97 KW. For gradient there is a call to the ONECENTR routine
c              which eliminates one-center quartets from calculations.
c              Lines 640-650.
c===================================================================
c SEP/24/96 KW. eleminated a double loop over pair-blocks from
c the BLOCKINT1 subroutine. This loop was executed to find out
c which ij- and kl-pair blocks made a given requested big block.
c Unfortunatly, this pice of code was executed every time the
c integral program was called i.e. for each small block.
c For big calculations (gramicidin /1504 b.f and almost 200 000
c big blocks ) it would take a LOT of time.
c It is gone now (see lines 204-271).
c===================================================================
c SEP/20/96 KW. density matrix here is NOT ordinary density. It is now
c denspar(icshell,jcshell) quadratic matrix, which used to be setup
c here. This setup is performed now in corresponding int_NAME routines
c (name=fock, store, giao, cphf).
c===================================================================
c SEP/20/96 KW. the DGETMO (library) routine is replaced by TRSPMO
c subroutine (put into the service.f file).(Transpose A(n,m)->B(m,n)
c===================================================================
c March 26,96 K.W. : call to wriut has been commented out.
c===================================================================
      subroutine calcint2(isupb_r,bl,inx,thres, denspar, wherex ,labels,
     *                    ibuffz, nblsiz,nintez,ngctoz, moreint,stopnow)
      implicit real*8 (a-h,o-z)
      logical firstd
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
c     where.eq.'sch1' Schwarz integrals for SCF
c     where.eq.'sch2' Schwarz integrals for NMR
c     where.eq.'sch3' Schwarz integrals for GRADIENT
c     where.eq.'sch4' Schwarz integrals for HESSIAN
c
c     last sch2,sch3,sch4   - NOT used now (yet)
c
c----------------------------------------------------------------
      common /schwarz/ schw4    ! used only in this file
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
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
c     write(6,*)' from calcint : block=',isupb_r,' where=',where
c----------------------------------------------------------------
      if(where.eq.'sch1') where='fock'
      if(where.eq.'sch2') where='shif'
      if(where.eq.'sch3') where='forc'
      if(where.eq.'sch4') where='hess'
c----------------------------------------------------------------
      nintex=0
c----------------------------------------------------------------
      stopnow=.false.
      if(isupb_r.le.0) then
         call nerror(1,'calcint2',
     *                 ' negative block number; can not do it',
     *                   isupb_r,0)
      endif
      if(isupb_r.gt.nbloks) then
         stopnow=.true.
         return
      endif
c----------------------------------------------------------------
c Set up the integral threshold :
c
      call setup_thres(thres)
c----------------------------------------------------------------
c calculations of two-electron integrals in blocks
c Nbl2   - number of blocks of pairs
c Nbloks - number of super-blocks of quartets (from Prepint2)
c
        call blockint1(bl,nbl2,inx,bl(iisd),bl(jjsd),bl(ijbld),
     *          npard,bl(ncost),bl(mxsize),bl(nsupb),denspar,
     *          labels, isupb_r, moreint)
c
      ibuffz= ibuffx
      nblsiz= nblsix
      nintez= nintex
      ngctoz= ngctox
c
c----------------------------------------------------------------
      end
c===================================================================
      subroutine blockint1(bl,nbl2,inx,iis,jjs,ijbl,
     *                     npard, ncost,mxsize,nsupb,denspar,
     *                     labels, isupb_r, moreint)
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstd
      logical moreint
      character*11 scftype
      character*4 where
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file
c----------------------------------------------------------------
      common /superbl/ isupb,ibl,kbl
      common /ilepar/ lpartot,lpareal
      common /runtype/ scftype,where
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /lindvec/ lind,idensp
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension ncost(*),mxsize(*),nsupb(*)
      dimension denspar(*)
      dimension labels(*)
      dimension ibf2(3070)
c----------------------------------------------------------------
      save last0
      save iblok_r
c----------------------------------------------------------------
c for check-point :
      if(.not.moreint) then
         call getmem(1,last0)
         call retmem(1)
      endif
c----------------------------------------------------------------
cpp
      ibls=3070
c     counter for the second buffer
      isto=0
cpp
c----------------------------------------------------------------
      if(.not.moreint) then
         iblok_r=1
      else
         iblok_r=iblok_r+1
      endif
c------------------------------------------------------
c calculate integrals from one small block (iblok_r)
c of requested super-block isupb_r :
c
c nblokx=nsupb(isupb_r) is known from prepint2 ( blockin2)
c iprice=ncost(isupb_r) is known from prepint2 ( blockin2)
c
      nblokx=nsupb(isupb_r)  ! number of small blocks in this big-block
      nblok_r=nblokx
      ikbl_r =iblok_r
c
c pair-blocks of this super-block are :
c
      call get_ij_half(isupb_r, ibl_r, kbl_r)
c-------------------------------------------------
c
      iprice=ncost(isupb_r)
c
      IF(  (where.eq.'disk' .and. iprice.eq.1) .or.
     *     (where.eq.'core' .and. iprice.eq.1) .or.
     *      where.eq.'fock'                    .or.
     *      where.eq.'shif'                    .or.
     *      where.eq.'forc'                    .or.
     *      where.eq.'hess'                  ) THEN
c----------------------------------------------------------------------
c check out possibility of reusing previously calculated pair-data
c
      if(ikbl_r.eq.1 ) then    ! new super-block
            call prec2ij(ibl_r, bl,inx,bl(npard),nbl2, iis,jjs,ijbl)
         if(kbl_r.ne.ibl_r) then
            call prec2kl(kbl_r, bl,inx,bl(npard),nbl2, iis,jjs,ijbl)
         endif
      endif
c----------------------------------------------------------------------
c CALCULATE INTEGRALS :
c
      maxsize=mxsize(isupb_r)
      call onesuper1(isupb_r,ibl_r,kbl_r,nblok_r,maxsize,bl,inx,npard,
     *               nbl2, iis,jjs,ijbl,isto,ibls,ibf2,where,
     *               denspar,labels,ikbl_r)
c----------------------------------------------------------------------
c release memory reserved in prec2ij  WHEN a given super-block is finished
c
         if(ikbl_r.eq.nblok_r) then
            if(kbl_r.ne.ibl_r) then
               call retmem(1)
            endif
            call retmem(1)
         endif
c----------------------------------------------------------------------
c tell a user if there is anymore integrals to be calculated in this super-block
c
      if( ikbl_r.eq.nblok_r ) then
         moreint=.false.
      else
         moreint=.true.
      endif
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
      subroutine onesuper1(isupb,ibl,kbl, nblokx,maxsize,bl,inx,npard,
     *                     nbl2, iis,jjs,ijbl,isto,ibls,ibf2,where,
     *                     denspar,labels,iblokx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*4 where
      common /ijcsfl/ ijblokp,ijprevf,ijprevl,ijtprev,maxprev,ngcprev
c
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
      common /memor3/ nblok1d
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension denspar(*)
      dimension labels(*)
      dimension ibf2(3070)
c----------------------------------------------------------------
      if(iblokx.eq.1) then
      call mmark
      endif
c----------------------------------------------------------------
c Constract blocks of contracted shell qyuartets for a given
c super-block (isupb) ; set up common /memor2/ and /memor3/ .
c
      if(iblokx.eq.1) then
         call blockin4(isupb,ibl,kbl, nblokx,maxsize,
     *                 bl,inx,npard,nbl2)
      endif
c
c----------------------------------------------------------------
c loop over blocks belonging to the super-block isupb:
c
      if(iblokx.eq.1) then
         ijblokp=0
         ijprevf=0
         ijprevl=0
         ijtprev=0
         maxprev=0
         ngcprev=0
      endif
c
      ikbl=iblokx
c
cnoloop  do 200 ikbl=1,nblokx
ccccc    write(6,*)'from onesuper no=',isupb,'  block=',ikbl
           call oneblock(isupb,ikbl,bl, bl(nibld),bl(nkbld),
     *                   bl(nijbd),bl(nijed),bl(nklbd),bl(nkled),
     *                   bl(nblok1d),bl(nqrtd),
     *                   iis,jjs,ijbl,nbl2,inx,
     *                   isto,ibls,ibf2,where,
cccc *                   densmat,labels)
     *                   denspar,labels)
  200    continue
c----------------------------------------------------------------
      if(iblokx.eq.nblokx) then
      call retmark
      endif
c----------------------------------------------------------------
      end
c===================================================================
      subroutine oneblock(isupb,ikbl,bl,nibl,nkbl, nijb,nije,nklb,nkle,
     *                    nblok1,nqrt,
     *                    iis,jjs,ijbl,nbl2,inx,
     *                    isto,ibls,ibf2,where,
     *                    denspar,labels)
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first,firstd
      logical skip_sb
      logical stable
      character*4 where
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file
c
      common /interface/ibuffx,nblsix,nintex,ngctox
c
      common /route/ iroute
      common /howmany/ ntotal,noijkl,nopres,nohalf,nrimtot,nrimret
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
ctime
      common /timex/ tconv1,tconv2,ttrobs
      common /time0/ tprec2
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time2/ terint,ttrans
      common /time3/ tzeroi,tspeci
ctime
c*
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
cxx
      common /logic2/ len(1)
      common /logic4/ nfu(1)
      common /lindvec/ lind,idensp
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
      common /esti4kl/ iesti2kl
c
c nmr derivatives :
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
      common /symmet/ rnsym
c
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension nblok1(2,*),nqrt(*)
c
      dimension denspar(*)
      dimension labels(*)
cpp
      dimension ibf2(3070)
cpp
c------------------------------------------------------
      if(icheck.gt.0) return
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
c-------------------------------------------------------------
      ntotal=ntotal+nbls
c------------------------------------------------------
c for Schwarz integrals only diagonal blocks :
c
c2000
      if(schw4.eq.'sch1') then
         if(npkl.ne.0) then
            nintex=0
            go to 100
         endif
      endif
c-------------------------------------------------------------
c  first quartet of shells :
c
c     ijcs1=nblok1(1,1)
c     klcs1=nblok1(2,1)
ckw1999
      ijcs1=ijbl(ibl,nijbeg)
      klcs1=ijbl(kbl,nklbeg)
c-------------------------------------------------------------
c oct. 2000
c check if the whole small block can be skip because of schwarz int.
c get maximum schwarz integrals for this small block
c
      if(schw4.NE.'sch1') then
         call getival('schwarz',ischwarz)
         call smblock_neg(bl(ischwarz),ncs, ijcs1,klcs1,skip_sb)
         if(skip_sb) then
cccc        write(6,*)' bb=',isupb,' sb=',ikbl,' skip'
            nintex=0
            go to 100
         endif
      endif
c-------------------------------------------------------------
c     make all quartets present :
c
      call isymm1(nbls,bl(isymm))
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
c98 it is needed only for Tracy's recursive and since
c   this recursive was re-formulated it is now always
c   as follows :
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
         call precalc2_1(isupb,bl,mmax,mmax1,nhabcd,nfumax,nbl2,nbls,
     *                 inx,iis,jjs,ijbl,nblok1,
     *                 ibl,nijbeg,nijend,npij,
     *                 kbl,nklbeg,nklend,npklx,npkl)
      ELSE
         call precalc2_2(isupb,bl,mmax,mmax1, nfumax, nbl2,nbls,
     *                   inx,iis,jjs,ijbl,nblok1,
     *                   ibl,nijbeg,nijend,npij,
     *                   kbl,nklbeg,nklend,npklx,npkl)
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
            call onecentr(iis,jjs,inx, npij,npklx,npkl,
     *              ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *              bl(isymm),bl(ijcent),bl(klcent) )
            call retmem(2)
         endif
      endif
c-------------------------------------------------------------
c for Schwarz integrals (eliminate off-diagonal quartets):
c
      if(schw4.eq.'sch1') then
         if(npkl.ne.0) 
     $     call nerror(1,'oneblock',' off-diagonal block',0,0)
         if(kbl.ne.ibl)
     $     call nerror(1,'oneblock',' off-diagonal block',0,0)
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
         call schw_neg(nbls,ncs, npij,npklx,npkl, nijbeg,nklbeg,
     *                     nbl2,ibl,kbl,ijbl,iis,jjs, denspar,
     *                     bl(isymm), bl(idnsx),nblsp)    ! output
c2000
         if(nblsp.eq.0) then
           nintex=0
           ncalls=ncalls-1
           go to 110
         endif
c2000
c
c symmetry handling :
c
         if(nsym.gt.0) then
           call onesym(ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *                  bl(ijshp),nsym,bl(isymm),nblsym )
            noijkl=noijkl+nblsym
         else
            noijkl=noijkl+nbls
         endif
      endif
c-------------------------------------------------------------
c At this moment some of c.s.quartets may not be present due
c to the symmetry or due to neglect according to the CSHNEG .
c Thus, reduce the block size from  nbls ---> nblsp (present)
c (only present c.s.quartets are considered)
c                  and
c   define idx1 and idx2
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
      call getmem(lckl,iesti2kl)
      ncalls=ncalls+1
c
      call get_max4kl(nblsrem,npklx,bl(idnsx),lckl,bl(iecd),
     *                bl(iesti2kl))
c
c iesti2kl saved in common /esti4kl/
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
         IF( iroute.eq.1 ) THEN
            call erintsp_1(
     *           bl,first,nbls,acc,ikbl,npij,npklx,npkl,idnsx)
         ELSE
            call erintsp_2(
     *           bl,first,nbls,acc,ikbl,npij,npklx,npkl,idnsx)
         ENDIF
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
c-------------------------------------------------------------
      call check_stab(stable,where)
c
      IF( iroute.eq.1 ) THEN
         call erinteg_1(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                  lobsa,immax,kmmax, where,idnsx,stable)
      ELSE
         call erinteg_2(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                  lobsa,immax,kmmax, where,idnsx,stable)
      ENDIF
c-------------------------------------------------------------
      if( first ) then
          nintex=0
          go to 110
      endif
c-------------------------------------------------------------
c transformation d6->d5 , f10->f7,  g15->g9, h21->h11, i28->i13
cc
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
c (1) if(where.eq.'fock'or'disk'or'core' in bl(ibuf)
c (2) if(where.eq.'shif'                 in bl(ibuf+nbls*lnijkl*ngcd)
c
c Corresponding labels arrays will be constracted :
c Then, integrals and labels will be returned to the calling
c place and the USER is resposible to put them somewhere :
c Thus, regardless of the value of where, construct labels :
c
        length=1+4*nbls*ngcd
        lgenct=length+4
c
        call destret(ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *               bl(isymm),
     *               ikbl,nbls,iis,jjs,ncs,inx,
     *               bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,
     *       bl(indxr), labels(1),labels(length),labels(lgenct) )
c
cctest
      ijklen=ilen*jlen*klen*llen
      icrt=ijklen*nbls
c      write(*,*) 'after destret, ijklen, nbls=',lnijkl,nbls
c      write(*,*) 'itype,jtype,ktype,ltype ', itype,jtype,ktype,ltype
c      write(*,*)'ilen,jlen,klen,llen ',ilen,jlen,klen,llen
c      write(*,*)'labels'
c      do ibl=1,nbls
c        write(*,*) ((labels(4*ibl-4+k)+1),k=1,4)
c      end do
c      if(incrt.le.1000) then
c        write(*,'(5(e10.3,2x))')(bl(ibuf-1+k),k=1,incrt)
c      else
c        write(*,*) 'Only the first 1000 integrals are printed'
c        write(*,'(5(e10.3,2x))')(bl(ibuf-1+k),k=1,1000)
c      end if
cc gen.cont.
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
c for iroute=1 :
      subroutine erinteg_1(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
     *                     lobsa,immax,kmmax, where,idnsx,stable)
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
      common /memor1/ iisd,jjsd,ijbld
      common /memor3/ nblok1d
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
      call prec4neg_1(nbls,npij,npklx,npkl,lcij,lckl,
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
      subroutine erintsp_1(bl,first,nbls,acc,ikbl,npij,npklx,npkl,
     *                     idnsx)

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
      common /memor1/ iisd,jjsd,ijbld
      common /memor3/ nblok1d
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
     *               bl(igci),bl(igcj), bl(igcij) )
      endif
      lcij=lcij+1
         lckl=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         if(ngcd.gt.1) then
           call gcparkl(nbls, bl(idx2),npklx,
     *                  l3,l4,ngck1,ngcl1,ngckl,
     *                  bl(igck),bl(igcl), bl(igckl) )
      endif
         lckl=lckl+1
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
      call precspec_1(nbls,npij,npklx,npkl, lcij,lckl,
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
c for iroute=2 :
c
c
      subroutine erinteg_2(bl,first,nbls,acc, ikbl,npij,npklx,npkl,
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
      common /memor1/ iisd,jjsd,ijbld
      common /memor3/ nblok1d
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
      common /esti4kl/ iesti2kl
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
c2002
      if(.not.stable) then
         call set_boamax2(bl(icc),bl(idd),l3,l4,doc_max) ! used in xwpq_22
         call setrval('doc_max',doc_max)
      endif
c2002
         lckl=lckl+1
         if(ngcd.gt.1) then
           call gcpairs(nbls,lckl,ngck1,ngcl1,ngckl,
     *                  bl(igckl),bl(igcklx))
         endif
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
c      call second(t4b)
c
      call prec4neg_2(nbls,npij,npklx,npkl,lcij,lckl,
     1     lc12,lc34, bl(idx1),bl(idx2),bl(indxr),
     *     bl(ieab),bl(iecd),bl(idnsx),bl(iesti2kl),
     *     bl(jump),bl(isymfac),
     *     bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq), nsym,bl(isymm),
     * bl(irho),bl(irppq),bl(irhoapb),bl(irhocpd),bl(irys),bl(iconst),
     *                          nbls1,bl(indx) )
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
c               it is zero out in assemblx.f now
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
      subroutine erintsp_2(bl,first,nbls,acc,ikbl,npij,npklx,npkl,idnsx)

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
      common /memor1/ iisd,jjsd,ijbld
      common /memor3/ nblok1d
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
      common /esti4kl/ iesti2kl
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
      call precspec_2(
     $    nbls,          npij,          npklx,         npkl, 
     $    lcij,          lckl,          lc12,          lc34,
     $    bl(idx1),      bl(idx2),      bl(indxr),     bl(ieab),
     $    bl(iecd),      bl(idnsx),     bl(iesti2kl),  bl(jump),
     $    bl(isymfac),   bl(iapb),      bl(icpd),      bl(icij),
     $    bl(ickl),      bl(ixp),       bl(ixq),       bl(i1apb),
     $    bl(i1cpd),     bl(itxab),     bl(itxcd),     nsym,
     $    bl(isymm),     bl(irho),      bl(irys),      bl(iconst),
     $    bl(ixwp),      bl(ixwq),      nbls1,         bl(indx) )
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
c------------------------------------------------------
      end
c===================================================================
      subroutine transfor(bl,nbls,ibuf,ibuf2,itype,jtype,ktype,ltype,
     *              ilen,jlen,klen,llen,ilen1,jlen1,klen1,llen1)
      implicit real*8 (a-h,o-z)
c
      dimension bl(*)
c
        ilenx=ilen1
        jlenx=jlen1
        klenx=klen1
        llenx=llen1
c
       if(itype.eq.4) then
         call dtrans1(nbls,bl(ibuf),bl(ibuf2),jlenx*klenx*llenx)
       endif
       if(itype.eq.6) then
         call ftrans1(nbls,bl(ibuf),bl(ibuf2),jlenx,klenx,llenx)
       endif
       if(itype.gt.10) then
         if(itype.eq.11) then
           call gtrans1(nbls,bl(ibuf),bl(ibuf2),jlenx*klenx*llenx)
         end if
         if(itype.eq.12) then
           call htrans1(nbls,bl(ibuf),bl(ibuf2),jlenx*klenx*llenx)
         end if
         if(itype.eq.13) then
c           call gtrans1(nbls,bl(ibuf),bl(ibuf2),jlenx*klenx*llenx)
         end if
       end if
c
       ilenx=ilen
c
       if(jtype.eq.4) then
         call dtrans2(nbls,bl(ibuf),bl(ibuf2),ilenx,klenx*llenx)
       endif
       if(jtype.eq.6) then
         call ftrans2(nbls,bl(ibuf),bl(ibuf2),ilenx,klenx,llenx)
       end if
       if(jtype.gt.10) then
         if(jtype.eq.11) then
           call gtrans2(nbls,bl(ibuf),bl(ibuf2),ilenx,klenx*llenx)
         end if
         if(jtype.eq.12) then
           call htrans2(nbls,bl(ibuf),bl(ibuf2),ilenx,klenx*llenx)
         end if
         if(jtype.eq.13) then
c           call itrans2(nbls,bl(ibuf),bl(ibuf2),ilenx,klenx*llenx)
         end if
       end if
c
       jlenx=jlen
c
       if(ktype.eq.4) then
         call dtrans3(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx,llenx)
       endif
       if(ktype.eq.6) then
         call ftrans3(nbls,bl(ibuf),bl(ibuf2),ilenx,jlenx,llenx)
       end if
       if(ktype.gt.10) then
         if(ktype.eq.11) then
           call gtrans3(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx,llenx)
         endif
         if(ktype.eq.12) then
           call htrans3(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx,llenx)
         endif
         if(ktype.eq.13) then
c           call itrans3(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx,llenx)
         endif
       endif
c
       klenx=klen
c
       if(ltype.eq.4) then
         call dtrans4(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx*klenx)
       endif
       if(ltype.eq.6) then
         call ftrans4(nbls,bl(ibuf),bl(ibuf2),ilenx,jlenx,klenx)
       end if
       if(ltype.gt.10) then
         if(ltype.eq.11) then
           call gtrans4(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx*klenx)
         end if
         if(ltype.eq.12) then
           call htrans4(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx*klenx)
         end if
         if(ltype.eq.13) then
c           call itrans4(nbls,bl(ibuf),bl(ibuf2),ilenx*jlenx*klenx)
         end if
       endif
c
      end
c===================================================================
      subroutine dtrans1(nbls,a,buf,jkllen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:jkllen,1:6) integrals for nbls shell-quartets
c  TEMP. STORAGE c  buf(1:nbls,1:jkllen,1:5)
c
c     *** the order of the raw d functions is:
c     xx, yy, zz, xy, xz, yz
c     1   2   3   4   5   6
c     *** the order of the transformed integrals is
c     (1/12)**(1/2)(2zz-xx-yy), (1/2)(xx-yy),xy,xz,yz
c                 1                   2      3  4  5
c
c There had been a lot of unnecessary memory movement. Speedup is
c considerable. Probably other routines should be also transformed.
c Routine was bottleneck in one test. (gm)
c
      implicit real*8 (a-h,o-z)
      parameter (sqtw=0.2886751345905d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls*jkllen,6),buf(nbls*jkllen,5)
      data sqr40/ 6.32455532033676d0 /
c
      lentot=nbls*jkllen
      do jkl=1,lentot
        buf(jkl,1)=sqtw*(two*a(jkl,3)-a(jkl,1)-a(jkl,2))
        buf(jkl,2)=half*(a(jkl,1)-a(jkl,2))
      end do
      call dcopy(lentot*2,buf(1,1),1,a(1,1),1)
      call dcopy(lentot,a(1,4),1,a(1,3),1)  ! done with three calls to dcopy
      call dcopy(lentot,a(1,5),1,a(1,4),1)  ! cuz some vectorized versions might
      call dcopy(lentot,a(1,6),1,a(1,5),1)  ! fail if the two arrays share
                                            ! storage
cc      do jkl=1,jkllen
cc        do ibl=1,nbls
cc        buf(ibl,jkl,1)=sqtw*(two*a(ibl,jkl,3)-a(ibl,jkl,1)-a(ibl,jkl,2))
cc        end do
cc      end do
cc      do jkl=1,jkllen
cc        do ibl=1,nbls
cc          buf(ibl,jkl,2)=half*(a(ibl,jkl,1)-a(ibl,jkl,2))
cc        end do
cc      end do
cc      do jkl=1,jkllen
cc        do ibl=1,nbls
cc          buf(ibl,jkl,3)=a(ibl,jkl,4)
cc        end do
cc      end do
cc      do jkl=1,jkllen
cc        do ibl=1,nbls
cc          buf(ibl,jkl,4)=a(ibl,jkl,5)
cc        end do
cc      end do
cc      do jkl=1,jkllen
cc        do ibl=1,nbls
cc          buf(ibl,jkl,5)=a(ibl,jkl,6)
cc        end do
cc      end do
c
cc      call tfer(buf(1,1,1),a(1,1,1),nbls*jkllen*5)
c
      end
c===================================================================
      subroutine dtrans2(nbls,a,buf,ilen,kllen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:kllen,1:6,i:ilen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:kllen,1:5,1:ilen)
c
c     *** the order of the raw d functions is:
c     xx, yy, zz, xy, xz, yz
c     1   2   3   4   5   6
c     *** the order of the transformed integrals is
c     (1/12)**(1/2)(2zz-xx-yy), (1/2)(xx-yy),xy,xz,yz
c                 1                   2      3  4  5
c
      implicit real*8 (a-h,o-z)
      parameter (sqtw=0.2886751345905d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,kllen,6,ilen),buf(nbls,kllen,5,ilen)
      do i=1,ilen
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,1,i)=
     1            sqtw*(two*a(ibl,kl,3,i)-a(ibl,kl,1,i)-a(ibl,kl,2,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,2,i)=half*(a(ibl,kl,1,i)-a(ibl,kl,2,i))
        end do
      end do
      do kl=1,kllen
       do ibl=1,nbls
        buf(ibl,kl,3,i)=a(ibl,kl,4,i)
       end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,4,i)=a(ibl,kl,5,i)
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,5,i)=a(ibl,kl,6,i)
        end do
      end do
      end do
c
      call tfer(buf(1,1,1,1),a(1,1,1,1),nbls*ilen*kllen*5)
c
      end
c===================================================================
      subroutine dtrans3(nbls,a,buf,ijlen,llen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:llen,1:6,1:ijlen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:llen,1:5,1:ijlen)
c
c     *** the order of the raw d functions is:
c     xx, yy, zz, xy, xz, yz
c     1   2   3   4   5   6
c     *** the order of the transformed integrals is
c     (1/12)**(1/2)(2zz-xx-yy), (1/2)(xx-yy),xy,xz,yz
c                 1                   2      3  4  5
c
      implicit real*8 (a-h,o-z)
      parameter (sqtw=0.2886751345905d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,llen,6,ijlen),buf(nbls,llen,5,ijlen)
      do ij=1,ijlen
      do l=1,llen
        do ibl=1,nbls
        buf(ibl,l,1,ij)=
     1      sqtw*(two*a(ibl,l,3,ij)-a(ibl,l,1,ij)-a(ibl,l,2,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,2,ij)=half*(a(ibl,l,1,ij)-a(ibl,l,2,ij))
        end do
      end do
      do l=1,llen
       do ibl=1,nbls
        buf(ibl,l,3,ij)=a(ibl,l,4,ij)
       end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,4,ij)=a(ibl,l,5,ij)
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,5,ij)=a(ibl,l,6,ij)
        end do
      end do
      end do
c
      call tfer(buf(1,1,1,1),a(1,1,1,1),nbls*ijlen*llen*5)
c
      end
c===================================================================
      subroutine dtrans4(nbls,a,buf,ijklen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:6,1:ijklen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:5,1:ijklen)
c
c     *** the order of the raw d functions is:
c     xx, yy, zz, xy, xz, yz
c     1   2   3   4   5   6
c     *** the order of the transformed integrals is
c     (1/12)**(1/2)(2zz-xx-yy), (1/2)(xx-yy),xy,xz,yz
c                 1                   2      3  4  5
c
      implicit real*8 (a-h,o-z)
      parameter (sqtw=0.2886751345905d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,6,ijklen),buf(nbls,5,ijklen)
      data sqr40/ 6.32455532033676d0 /
      do ijk=1,ijklen
        do ibl=1,nbls
        buf(ibl,1,ijk)=sqtw*(two*a(ibl,3,ijk)-a(ibl,1,ijk)-a(ibl,2,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,2,ijk)=half*(a(ibl,1,ijk)-a(ibl,2,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,3,ijk)=a(ibl,4,ijk)
        end do
        do ibl=1,nbls
          buf(ibl,4,ijk)=a(ibl,5,ijk)
        end do
        do ibl=1,nbls
          buf(ibl,5,ijk)=a(ibl,6,ijk)
        end do
      end do
c
      call tfer(buf(1,1,1),a(1,1,1),nbls*ijklen*5)
c
      end
c===================================================================
      subroutine ftrans1(nbls,a,buf,jlen,klen,llen)
c
c     ***   <fj,kl>  ***
c
c     *** the order of the raw f functions xxx,xxy,xxz,xyy,xyz,
c     *** xzz,yyy,yyz,yzz,zzz
c     *** the order of the transformed integrals (5xxy-rry), (5xxz-rrz),
c     *** (5yyx-rrx), (5yyz-rrz), (5zzx-rrx), (5zzy-rry), xyz
c     *** these are not orthogonal and are not normalized now
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,*),buf(nbls,*)
      dimension in0(10)
      data sqr40/ 6.32455532033676d0 /
c
 1000 continue
      istep3=llen
      istep2=istep3*klen
      istep1=istep2*jlen
c
      ncount=istep1*7
c
      in0(1)=0
      do 10 i=2,10
      in0(i)=in0(i-1)+istep1
   10 continue
c
      do 100 j=1,jlen
      jj0=(j-1)*istep2
      do 100 k=1,klen
      kk0=(k-1)*istep3
      jk0=jj0+kk0
      do 100 l=1,llen
      ind0=jk0+l
      ind7=ind0
      i10=in0(1)+ind0
      i20=in0(2)+ind0
      i30=in0(3)+ind0
      i40=in0(4)+ind0
      i50=in0(5)+ind0
      i60=in0(6)+ind0
      i70=in0(7)+ind0
      i80=in0(8)+ind0
      i90=in0(9)+ind0
      i00=in0(10)+ind0
c
      i17=in0(1) +ind7
      i27=in0(2) +ind7
      i37=in0(3) +ind7
      i47=in0(4) +ind7
      i57=in0(5) +ind7
      i67=in0(6) +ind7
      i77=in0(7) +ind7
cxxx
      do 100 n=1,nbls
c
      buf(n,i17)=four*a(n,i20)-a(n,i70)-a(n,i90)
      buf(n,i27)=four*a(n,i30)-a(n,i80)-a(n,i00)
      buf(n,i37)=four*a(n,i40)-a(n,i10)-a(n,i60)
      buf(n,i47)=four*a(n,i80)-a(n,i30)-a(n,i00)
      buf(n,i57)=four*a(n,i60)-a(n,i10)-a(n,i40)
      buf(n,i67)=four*a(n,i90)-a(n,i20)-a(n,i70)
      buf(n,i77)=     a(n,i50)*sqr40
  100 continue
c
      call tfer(buf,a,ncount*nbls)
c      do 110 i=1,ncount
c      do 110 n=1,nbls
c  110 a(n,i)=buf(n,i)
c      return
      end
c===================================================================
      subroutine ftrans2(nbls,a,buf,ilen,klen,llen)
c
c     ***   <if,kl>  ***
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,*),buf(nbls,*)
      dimension in0(10)
      data sqr40/ 6.32455532033676d0 /
c
 1000 continue
c
      istep3=llen
      istep2=istep3*klen
      istep10=istep2*10
      istep17=istep2*7
c
      ncount=istep17*ilen
cxxxxxxxxxxxxxxxxxxxxxxxxxxx
      in0(1)=0
      do 10 i=2,10
      in0(i)=in0(i-1)+istep2
   10 continue
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
      do 100 i=1,ilen
      ii0=(i-1)*istep10
      ii7=(i-1)*istep17
      do 100 k=1,klen
      kk=(k-1)*istep3
      ik0=ii0+kk
      ik7=ii7+kk
      do 100 l=1,llen
      ind0=ik0+l
      ind7=ik7+l
c
      i10=in0(1) +ind0
      i20=in0(2) +ind0
      i30=in0(3) +ind0
      i40=in0(4) +ind0
      i50=in0(5) +ind0
      i60=in0(6) +ind0
      i70=in0(7) +ind0
      i80=in0(8) +ind0
      i90=in0(9) +ind0
      i00=in0(10) +ind0
c
      i17=in0(1) +ind7
      i27=in0(2) +ind7
      i37=in0(3) +ind7
      i47=in0(4) +ind7
      i57=in0(5) +ind7
      i67=in0(6) +ind7
      i77=in0(7) +ind7
cxxx
      do 100 n=1,nbls
c
      buf(n,i17)=four*a(n,i20)-a(n,i70)-a(n,i90)
      buf(n,i27)=four*a(n,i30)-a(n,i80)-a(n,i00)
      buf(n,i37)=four*a(n,i40)-a(n,i10)-a(n,i60)
      buf(n,i47)=four*a(n,i80)-a(n,i30)-a(n,i00)
      buf(n,i57)=four*a(n,i60)-a(n,i10)-a(n,i40)
      buf(n,i67)=four*a(n,i90)-a(n,i20)-a(n,i70)
      buf(n,i77)=     a(n,i50)*sqr40
  100 continue
c
      call tfer(buf,a,ncount*nbls)
c
c      do 110 i=1,ncount
c      do 110 n=1,nbls
c  110 a(n,i)=buf(n,i)
c      return
      end
c===================================================================
      subroutine ftrans3(nbls,a,buf,ilen,jlen,llen)
c
c     ***   <ij,fl>  ***
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,*),buf(nbls,*)
      dimension in0(10)
      data sqr40/ 6.32455532033676d0 /
c
 1000 continue
c     istep4=1
      istep3=llen
      istep20=istep3*10
      istep27=istep3*7
      istep10=istep20*jlen
      istep17=istep27*jlen
c
      ncount=istep17*ilen
c
      in0(1)=0
      do 10 i=2,10
      in0(i)=in0(i-1)+istep3
   10 continue
c
      do 100 i=1,ilen
      ii0=(i-1)*istep10
      ii7=(i-1)*istep17
      do 100 j=1,jlen
      jj0=(j-1)*istep20
      jj7=(j-1)*istep27
      ij0=ii0+jj0
      ij7=ii7+jj7
      do 100 l=1,llen
      ind0=ij0+l
      ind7=ij7+l
c
      i10=in0(1) +ind0
      i20=in0(2) +ind0
      i30=in0(3) +ind0
      i40=in0(4) +ind0
      i50=in0(5) +ind0
      i60=in0(6) +ind0
      i70=in0(7) +ind0
      i80=in0(8) +ind0
      i90=in0(9) +ind0
      i00=in0(10) +ind0
c
      i17=in0(1) +ind7
      i27=in0(2) +ind7
      i37=in0(3) +ind7
      i47=in0(4) +ind7
      i57=in0(5) +ind7
      i67=in0(6) +ind7
      i77=in0(7) +ind7
cxxx
      do 100 n=1,nbls
c
      buf(n,i17)=four*a(n,i20)-a(n,i70)-a(n,i90)
      buf(n,i27)=four*a(n,i30)-a(n,i80)-a(n,i00)
      buf(n,i37)=four*a(n,i40)-a(n,i10)-a(n,i60)
      buf(n,i47)=four*a(n,i80)-a(n,i30)-a(n,i00)
      buf(n,i57)=four*a(n,i60)-a(n,i10)-a(n,i40)
      buf(n,i67)=four*a(n,i90)-a(n,i20)-a(n,i70)
      buf(n,i77)=     a(n,i50)*sqr40
  100 continue
c
c
      call tfer(buf,a,ncount*nbls)
c      do 110 i=1,ncount
c      do 110 n=1,nbls
c  110 a(n,i)=buf(n,i)
c      return
      end
c===================================================================
      subroutine ftrans4(nbls,a,buf,ilen,jlen,klen)
c
c     ***   <ij,kf>  ***
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,*),buf(nbls,*)
      dimension in0(10)
      data sqr40/ 6.32455532033676d0 /
c
c
 1000 continue
      istep4=1
      istep30=10
      istep37=7
      istep20=istep30*klen
      istep27=istep37*klen
      istep10=istep20*jlen
      istep17=istep27*jlen
c
      ncount=istep17*ilen
c
c     in0(1)=0
      in0(1)=istep4
      do 10 i=2,10
      in0(i)=in0(i-1)+istep4
   10 continue
c
c
      do 100 i=1,ilen
      ii0=(i-1)*istep10
      ii7=(i-1)*istep17
      do 100 j=1,jlen
      jj0=(j-1)*istep20
      jj7=(j-1)*istep27
      ij0=ii0+jj0
      ij7=ii7+jj7
      do 100 k=1,klen
      kk0=(k-1)*istep30
      kk7=(k-1)*istep37
c
      ind0=ij0+kk0
      ind7=ij7+kk7
c
      i10=in0(1) +ind0
      i20=in0(2) +ind0
      i30=in0(3) +ind0
      i40=in0(4) +ind0
      i50=in0(5) +ind0
      i60=in0(6) +ind0
      i70=in0(7) +ind0
      i80=in0(8) +ind0
      i90=in0(9) +ind0
      i00=in0(10) +ind0
c
      i17=in0(1) +ind7
      i27=in0(2) +ind7
      i37=in0(3) +ind7
      i47=in0(4) +ind7
      i57=in0(5) +ind7
      i67=in0(6) +ind7
      i77=in0(7) +ind7
cxxx
      do 100 n=1,nbls
c
      buf(n,i17)=four*a(n,i20)-a(n,i70)-a(n,i90)
      buf(n,i27)=four*a(n,i30)-a(n,i80)-a(n,i00)
      buf(n,i37)=four*a(n,i40)-a(n,i10)-a(n,i60)
      buf(n,i47)=four*a(n,i80)-a(n,i30)-a(n,i00)
      buf(n,i57)=four*a(n,i60)-a(n,i10)-a(n,i40)
      buf(n,i67)=four*a(n,i90)-a(n,i20)-a(n,i70)
      buf(n,i77)=     a(n,i50)*sqr40
  100 continue
c
c
      call tfer(buf,a,ncount*nbls)
c
c      do 110 i=1,ncount
c      do 110 n=1,nbls
c  110 a(n,i)=buf(n,i)
c      return
      end
c===================================================================
      subroutine gtrans1(nbls,a,buf,jkllen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:jkllen,1:15) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:jkllen,1:9)
c
c     *** the order of the raw g functions is:
c     xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
c       1     2     3     4     5    6      7    8      9    10
c     yyyy,  yyyz, yyzz, yzzz,  zzzz
c      11     12    13    14     15
c     *** the order of the transformed integrals is
c     xxxx+yyyy-6xxyy, xxxx+zzzz-6xxzz, yyyy+zzzz-6yyzz
c            1                        2                        3
c     (xxxy-3xyzz), (xyyy-3xyzz), (xxxz-3xyyz), (xzzz-3xyyz).
c          4             5             6             7
c     (yyyz-3xxyz), (yzzz-3xxyz)
c          8             9
c
      implicit real*8 (a-h,o-z)
      parameter (six=6.0d0,yn=3.2659863237109d0,xn=1.15470053837925d0)
c   3.265986.. is sqrt(32/3); 1.1547005 is sqrt(4/3)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,jkllen,15),buf(nbls,jkllen,9)
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,1)=
     1    xn*(a(ibl,jkl,1)+a(ibl,jkl,11)-six*a(ibl,jkl,4))
ctest
c        write(*,*) ibl, a(ibl,jkl,1),a(ibl,jkl,11),a(ibl,jkl,4),
c    1  'buf=',buf(ibl,jkl,1)
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,2)=
     1    xn*(a(ibl,jkl,1)+a(ibl,jkl,15)-six*a(ibl,jkl,6))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,3)=
     1    xn*(a(ibl,jkl,11)+a(ibl,jkl,15)-six*a(ibl,jkl,13))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,4)=yn*(a(ibl,jkl,2)-three*a(ibl,jkl,9))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,5)=yn*(a(ibl,jkl,7)-three*a(ibl,jkl,9))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,6)=yn*(a(ibl,jkl,3)-three*a(ibl,jkl,8))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,7)=yn*(a(ibl,jkl,10)-three*a(ibl,jkl,8))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,8)=yn*(a(ibl,jkl,12)-three*a(ibl,jkl,5))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,9)=yn*(a(ibl,jkl,14)-three*a(ibl,jkl,5))
        end do
      end do
c
      call tfer(buf(1,1,1),a(1,1,1),nbls*jkllen*9)
c
      end
c===================================================================
      subroutine gtrans2(nbls,a,buf,ilen,kllen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:kllen,1:15,i:ilen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:kllen,1:9,1:ilen)
c
c     *** the order of the raw g functions is:
c     xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
c       1     2     3     4     5    6      7    8      9    10
c     yyyy,  yyyz, yyzz, yzzz,  zzzz
c      11     12    13    14     15
c     *** the order of the transformed integrals is
c     xxxx+yyyy-6xxyy, xxxx+zzzz-6xxzz, yyyy+zzzz-6yyzz
c            1                        2                        3
c     (xxxy-3xyzz), (xyyy-3xyzz), (xxxz-3xyyz), (xzzz-3xyyz).
c          4             5             6             7
c     (yyyz-3xxyz), (yzzz-3xxyz)
c          8             9
c
      implicit real*8 (a-h,o-z)
      parameter (six=6.0d0,yn=3.2659863237109d0,xn=1.15470053837925d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,kllen,15,ilen),buf(nbls,kllen,9,ilen)
      do i=1,ilen
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,1,i)=
     1        xn*(a(ibl,kl,1,i)+a(ibl,kl,11,i)-six*a(ibl,kl,4,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,2,i)=
     1    xn*(a(ibl,kl,1,i)+a(ibl,kl,15,i)-six*a(ibl,kl,6,i))
        end do
      end do
      do kl=1,kllen
       do ibl=1,nbls
        buf(ibl,kl,3,i)=
     1      xn*(a(ibl,kl,11,i)+a(ibl,kl,15,i)-six*a(ibl,kl,13,i))
       end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,4,i)=yn*(a(ibl,kl,2,i)-three*a(ibl,kl,9,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,5,i)=yn*(a(ibl,kl,7,i)-three*a(ibl,kl,9,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,6,i)=yn*(a(ibl,kl,3,i)-three*a(ibl,kl,8,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,7,i)=yn*(a(ibl,kl,10,i)-three*a(ibl,kl,8,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,8,i)=yn*(a(ibl,kl,12,i)-three*a(ibl,kl,5,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,9,i)=yn*(a(ibl,kl,14,i)-three*a(ibl,kl,5,i))
        end do
      end do
      end do
c
      call tfer(buf(1,1,1,1),a(1,1,1,1),nbls*ilen*kllen*9)
c
      end
c===================================================================
      subroutine gtrans3(nbls,a,buf,ijlen,llen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:llen,1:15,1:ijlen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:llen,1:9,1:ijlen)
c
c     *** the order of the raw g functions is:
c     xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
c       1     2     3     4     5    6      7    8      9    10
c     yyyy,  yyyz, yyzz, yzzz,  zzzz
c      11     12    13    14     15
c     *** the order of the transformed integrals is
c     xxxx+yyyy-6xxyy, xxxx+zzzz-6xxzz, yyyy+zzzz-6yyzz
c            1                        2                        3
c     (xxxy-3xyzz), (xyyy-3xyzz), (xxxz-3xyyz), (xzzz-3xyyz).
c          4             5             6             7
c     (yyyz-3xxyz), (yzzz-3xxyz)
c          8             9
c
      implicit real*8 (a-h,o-z)
      parameter (six=6.0d0,yn=3.2659863237109d0,xn=1.15470053837925d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,llen,15,ijlen),buf(nbls,llen,9,ijlen)
      do ij=1,ijlen
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,1,ij)=
     1    xn*(a(ibl,l,1,ij)+a(ibl,l,11,ij)-six*a(ibl,l,4,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,2,ij)=
     1    xn*(a(ibl,l,1,ij)+a(ibl,l,15,ij)-six*a(ibl,l,6,ij))
        end do
      end do
      do l=1,llen
       do ibl=1,nbls
        buf(ibl,l,3,ij)=
     1      xn*(a(ibl,l,11,ij)+a(ibl,l,15,ij)-six*a(ibl,l,13,ij))
       end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,4,ij)=yn*(a(ibl,l,2,ij)-three*a(ibl,l,9,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,5,ij)=yn*(a(ibl,l,7,ij)-three*a(ibl,l,9,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,6,ij)=yn*(a(ibl,l,3,ij)-three*a(ibl,l,8,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,7,ij)=yn*(a(ibl,l,10,ij)-three*a(ibl,l,8,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,8,ij)=yn*(a(ibl,l,12,ij)-three*a(ibl,l,5,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,9,ij)=yn*(a(ibl,l,14,ij)-three*a(ibl,l,5,ij))
        end do
      end do
      end do
c
      call tfer(buf(1,1,1,1),a(1,1,1,1),nbls*ijlen*llen*9)
c
      end
c===================================================================
      subroutine gtrans4(nbls,a,buf,ijklen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:15,1:ijklen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:9,1:ijklen)
c
c     *** the order of the raw g functions is:
c     xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
c       1     2     3     4     5    6      7    8      9    10
c     yyyy,  yyyz, yyzz, yzzz,  zzzz
c      11     12    13    14     15
c     *** the order of the transformed integrals is
c     xxxx+yyyy-6xxyy, xxxx+zzzz-6xxzz, yyyy+zzzz-6yyzz
c            1                        2                        3
c     (xxxy-3xyzz), (xyyy-3xyzz), (xxxz-3xyyz), (xzzz-3xyyz).
c          4             5             6             7
c     (yyyz-3xxyz), (yzzz-3xxyz)
c          8             9
c
      implicit real*8 (a-h,o-z)
      parameter (six=6.0d0,yn=3.2659863237109d0,xn=1.15470053837925d0)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension a(nbls,15,ijklen),buf(nbls,9,ijklen)
      do ijk=1,ijklen
        do ibl=1,nbls
          buf(ibl,1,ijk)=
     1    xn*(a(ibl,1,ijk)+a(ibl,11,ijk)-six*a(ibl,4,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,2,ijk)=
     1    xn*(a(ibl,1,ijk)+a(ibl,15,ijk)-six*a(ibl,6,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,3,ijk)=
     1    xn*(a(ibl,11,ijk)+a(ibl,15,ijk)-six*a(ibl,13,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,4,ijk)=yn*(a(ibl,2,ijk)-three*a(ibl,9,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,5,ijk)=yn*(a(ibl,7,ijk)-three*a(ibl,9,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,6,ijk)=yn*(a(ibl,3,ijk)-three*a(ibl,8,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,7,ijk)=yn*(a(ibl,10,ijk)-three*a(ibl,8,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,8,ijk)=yn*(a(ibl,12,ijk)-three*a(ibl,5,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,9,ijk)=yn*(a(ibl,14,ijk)-three*a(ibl,5,ijk))
        end do
      end do
c
      call tfer(buf(1,1,1),a(1,1,1),nbls*ijklen*9)
c
      end
c===================================================================
      subroutine htrans1(nbls,a,buf,jkllen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:jkllen,1:21) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:jkllen,1:11)
c
c     *** the order of the raw h functions is:
c     x5,    x4y,   x4z,   x3y2,  x3yz,  x3z2,  x2y3,  x2y2z, x2yz2, x2z3,
c     1       2      3      4      5      6      7       8      9     10
c     xy4,   xy3z,  xy3z2, xyz3,  xz4.   y5,    y4z,   y3z2.  y2z3.  yz4,    z5
c     11      12     13     14     15    16     17      18     19    20      21
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
c
      implicit real*8 (a-h,o-z)
      PARAMETER (two=2.0d0,three=3.0d0, five=5.0d0,six=6.0d0,ten=10.0d0)
      parameter (xn=0.730296743340221d0,yn=2.3094010767585d0,
     1    zn=9.23760430703401d0,sn=16.0d0/3.0d0)
c  xn=4/sqrt(30)  yn=4/sqrt(3)  zn=16/sqrt(3)   sn=16/3
      dimension a(nbls,jkllen,21),buf(nbls,jkllen,11)
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,1)=
     1    xn*(a(ibl,jkl,1)-ten*a(ibl,jkl,4)+five*a(ibl,jkl,11))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,2)=
     1    xn*(a(ibl,jkl,1)-ten*a(ibl,jkl,6)+five*a(ibl,jkl,15))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,3)=
     1    xn*(a(ibl,jkl,16)-ten*a(ibl,jkl,7)+five*a(ibl,jkl,2))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,4)=
     1    xn*(a(ibl,jkl,16)-ten*a(ibl,jkl,18)+five*a(ibl,jkl,20))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,5)=
     1    xn*(a(ibl,jkl,21)-ten*a(ibl,jkl,10)+five*a(ibl,jkl,3))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,6)=
     1    xn*(a(ibl,jkl,21)-ten*a(ibl,jkl,19)+five*a(ibl,jkl,17))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,7)=
     1    yn*(a(ibl,jkl,3)+a(ibl,jkl,17)-six*a(ibl,jkl,8))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,8)=
     1    yn*(a(ibl,jkl,2)+a(ibl,jkl,20)-six*a(ibl,jkl,9))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,9)=
     1    yn*(a(ibl,jkl,11)+a(ibl,jkl,15)-six*a(ibl,jkl,13))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,10)=zn*(a(ibl,jkl,5)-a(ibl,jkl,12))
        end do
      end do
      do jkl=1,jkllen
        do ibl=1,nbls
          buf(ibl,jkl,11)=
     1    sn*(a(ibl,jkl,5)-two*a(ibl,jkl,14)+a(ibl,jkl,12))
        end do
      end do
c
      call tfer(buf(1,1,1),a(1,1,1),nbls*jkllen*11)
c
      end
c===================================================================
      subroutine htrans2(nbls,a,buf,ilen,kllen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:kllen,1:21,i:ilen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:kllen,1:11,1:ilen)
c
c     *** the order of the raw h functions is:
c     x5,    x4y,   x4z,   x3y2,  x3yz,  x3z2,  x2y3,  x2y2z, x2yz2, x2z3,
c     1       2      3      4      5      6      7       8      9     10
c     xy4,   xy3z,  xy3z2, xyz3,  xz4.   y5,    y4z,   y3z2.  y2z3.  yz4,    z5
c     11      12     13     14     15    16     17      18     19    20      21
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
c
      implicit real*8 (a-h,o-z)
      PARAMETER (two=2.0d0,three=3.0d0, five=5.0d0,six=6.0d0,ten=10.0d0)
      parameter (xn=0.730296743340221d0,yn=2.3094010767585d0,
     1   zn=9.23760430703401d0,sn=16.0d0/3.0d0)
c  xn=4/sqrt(30)  yn=4/sqrt(3)  zn=16/sqrt(3)   sn=16/3
      dimension a(nbls,kllen,21,ilen),buf(nbls,kllen,11,ilen)
      do i=1,ilen
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,1,i)=
     1      xn*(a(ibl,kl,1,i)-ten*a(ibl,kl,4,i)+five*a(ibl,kl,11,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,2,i)=
     1    xn*(a(ibl,kl,1,i)-ten*a(ibl,kl,6,i)+five*a(ibl,kl,15,i))
        end do
      end do
      do kl=1,kllen
       do ibl=1,nbls
        buf(ibl,kl,3,i)=
     1   xn*(a(ibl,kl,16,i)-ten*a(ibl,kl,7,i)+five*a(ibl,kl,2,i))
       end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,4,i)=
     1    xn*(a(ibl,kl,16,i)-ten*a(ibl,kl,18,i)+five*a(ibl,kl,20,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,5,i)=
     1    xn*(a(ibl,kl,21,i)-ten*a(ibl,kl,10,i)+five*a(ibl,kl,3,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,6,i)=
     1    xn*(a(ibl,kl,21,i)-ten*a(ibl,kl,19,i)+five*a(ibl,kl,17,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,7,i)=
     1    yn*(a(ibl,kl,3,i)+a(ibl,kl,17,i)-six*a(ibl,kl,8,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,8,i)=
     1    yn*(a(ibl,kl,2,i)+a(ibl,kl,20,i)-six*a(ibl,kl,9,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,9,i)=
     1    yn*(a(ibl,kl,11,i)+a(ibl,kl,15,i)-six*a(ibl,kl,13,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,10,i)=
     1    zn*(a(ibl,kl,5,i)-a(ibl,kl,12,i))
        end do
      end do
      do kl=1,kllen
        do ibl=1,nbls
          buf(ibl,kl,11,i)=
     1    sn*(a(ibl,kl,5,i)-two*a(ibl,kl,14,i)+a(ibl,kl,12,i))
        end do
      end do
      end do
c
      call tfer(buf(1,1,1,1),a(1,1,1,1),nbls*ilen*kllen*11)
c
      end
c===================================================================
      subroutine htrans3(nbls,a,buf,ijlen,llen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:llen,1:21,1:ijlen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:llen,1:11,1:ijlen)
c
c     *** the order of the raw h functions is:
c     x5,    x4y,   x4z,   x3y2,  x3yz,  x3z2,  x2y3,  x2y2z, x2yz2, x2z3,
c     1       2      3      4      5      6      7       8      9     10
c     xy4,   xy3z,  xy3z2, xyz3,  xz4.   y5,    y4z,   y3z2.  y2z3.  yz4,    z5
c     11      12     13     14     15    16     17      18     19    20      21
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
c
      implicit real*8 (a-h,o-z)
      PARAMETER (two=2.0d0,three=3.0d0, five=5.0d0,six=6.0d0,ten=10.0d0)
      parameter (xn=0.730296743340221d0,yn=2.3094010767585d0,
     1    zn=9.23760430703401d0,sn=16.0d0/3.0d0)
c  xn=4/sqrt(30)  yn=4/sqrt(3)  zn=16/sqrt(3)   sn=16/3
      dimension a(nbls,llen,21,ijlen),buf(nbls,llen,11,ijlen)
      do ij=1,ijlen
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,1,ij)=
     1    xn*(a(ibl,l,1,ij)-ten*a(ibl,l,4,ij)+five*a(ibl,l,11,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,2,ij)=
     1    xn*(a(ibl,l,1,ij)-ten*a(ibl,l,6,ij)+five*a(ibl,l,15,ij))
        end do
      end do
      do l=1,llen
       do ibl=1,nbls
        buf(ibl,l,3,ij)=
     1    xn*(a(ibl,l,16,ij)-ten*a(ibl,l,7,ij)+five*a(ibl,l,2,ij))
       end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,4,ij)=
     1    xn*(a(ibl,l,16,ij)-ten*a(ibl,l,18,ij)+five*a(ibl,l,20,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,5,ij)=
     1      xn*(a(ibl,l,21,ij)-ten*a(ibl,l,10,ij)+five*a(ibl,l,3,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,6,ij)=
     1    xn*(a(ibl,l,21,ij)-ten*a(ibl,l,19,ij)+five*a(ibl,l,17,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,7,ij)=
     1    yn*(a(ibl,l,3,ij)+a(ibl,l,17,ij)-six*a(ibl,l,8,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,8,ij)=
     1    yn*(a(ibl,l,2,ij)+a(ibl,l,20,ij)-six*a(ibl,l,9,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,9,ij)=
     1    yn*(a(ibl,l,11,ij)+a(ibl,l,15,ij)-six*a(ibl,l,13,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,10,ij)=
     1    zn*(a(ibl,l,5,ij)-a(ibl,l,12,ij))
        end do
      end do
      do l=1,llen
        do ibl=1,nbls
          buf(ibl,l,11,ij)=
     1    sn*(a(ibl,l,5,ij)-two*a(ibl,l,14,ij)+a(ibl,l,12,ij))
        end do
      end do
      end do
c
      call tfer(buf(1,1,1,1),a(1,1,1,1),nbls*ijlen*llen*11)
c
      end
c===================================================================
      subroutine htrans4(nbls,a,buf,ijklen)
c  Arguments:
c  INTENT(IN)
c  nbls=blocksize
c  jlen,klen,llen: shell sizes (e.g. d=5) of the shells J,K,L in (IJ|KL);
c  jklen=jlen*klen*llen
c  INTENT(INOUT)
c  a(1:nbls,1:21,1:ijklen) integrals for nbls shell-quartets
c  TEMP. STORAGE
c  buf(1:nbls,1:11,1:ijklen)
c
c     *** the order of the raw h functions is:
c     x5,    x4y,   x4z,   x3y2,  x3yz,  x3z2,  x2y3,  x2y2z, x2yz2, x2z3,
c     1       2      3      4      5      6      7       8      9     10
c     xy4,   xy3z,  xy3z2, xyz3,  xz4.   y5,    y4z,   y3z2.  y2z3.  yz4,    z5
c     11      12     13     14     15    16     17      18     19    20      21
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
c
      implicit real*8 (a-h,o-z)
      PARAMETER (two=2.0d0,three=3.0d0, five=5.0d0,six=6.0d0,ten=10.0d0)
      parameter (xn=0.730296743340221d0,yn=2.3094010767585d0,
     1    zn=9.23760430703401d0,sn=16.0d0/3.0d0)
c  xn=4/sqrt(30)  yn=4/sqrt(3)  zn=16/sqrt(3)   sn=16/3
      dimension a(nbls,21,ijklen),buf(nbls,11,ijklen)
      do ijk=1,ijklen
        do ibl=1,nbls
          buf(ibl,1,ijk)=
     1    xn*(a(ibl,1,ijk)-ten*a(ibl,4,ijk)+five*a(ibl,11,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,2,ijk)=
     1    xn*(a(ibl,1,ijk)-ten*a(ibl,6,ijk)+five*a(ibl,15,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,3,ijk)=
     1    xn*(a(ibl,16,ijk)-ten*a(ibl,7,ijk)+five*a(ibl,2,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,4,ijk)=
     1    xn*(a(ibl,16,ijk)-ten*a(ibl,18,ijk)+five*a(ibl,20,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,5,ijk)=
     1    xn*(a(ibl,21,ijk)-ten*a(ibl,10,ijk)+five*a(ibl,3,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,6,ijk)=
     1    xn*(a(ibl,21,ijk)-ten*a(ibl,19,ijk)+five*a(ibl,17,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,7,ijk)=
     1    yn*(a(ibl,3,ijk)+a(ibl,17,ijk)-six*a(ibl,8,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,8,ijk)=
     1    yn*(a(ibl,2,ijk)+a(ibl,20,ijk)-six*a(ibl,9,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,9,ijk)=
     1    yn*(a(ibl,11,ijk)+a(ibl,15,ijk)-six*a(ibl,13,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,10,ijk)=
     1    zn*(a(ibl,5,ijk)-a(ibl,12,ijk))
        end do
        do ibl=1,nbls
          buf(ibl,11,ijk)=
     1    sn*(a(ibl,5,ijk)-two*a(ibl,14,ijk)+a(ibl,12,ijk))
        end do
      end do
      call tfer(buf(1,1,1),a(1,1,1),nbls*ijklen*11)
c
      end
c===================================================================
      subroutine conv1x_1(nbls1,mmax1,npij,lcij, idx1,indx,
     *                    abnia,xpn,abnix,xpnx )
      implicit real*8 (a-h,o-z)
c
      dimension idx1(*),indx(*)
      dimension xpn(npij,3,*)
      dimension abnia(npij,mmax1,*)
c
      dimension xpnx(nbls1,3)
      dimension abnix(nbls1,mmax1)
c
      do 10 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        xpnx(i,1)=xpn(ijpar,1,lcij)
        xpnx(i,2)=xpn(ijpar,2,lcij)
        xpnx(i,3)=xpn(ijpar,3,lcij)
   10 continue
c
      do 20 m=1,mmax1
      do 20 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        abnix(i,m)=abnia(ijpar,m,lcij)
   20 continue
c
      end
c===================================================================
      subroutine conv1x_2(nbls1,mmax1,npij,lcij, idx1,indx,
     *                    xpn,xpnx )
      implicit real*8 (a-h,o-z)
      dimension idx1(*),indx(*)
      dimension xpn(npij,3,*)
      dimension xpnx(nbls1,3)
c
      do 10 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        xpnx(i,1)=xpn(ijpar,1,lcij)
        xpnx(i,2)=xpn(ijpar,2,lcij)
        xpnx(i,3)=xpn(ijpar,3,lcij)
   10 continue
c
      end
c===================================================================
      subroutine conv2x(nbls1,nfumax1,npkl,lckl, idx2,indx,
     *                  habcd,nfumax, habcdx  )
      implicit real*8 (a-h,o-z)
      dimension idx2(*),indx(*)
      dimension habcd(npkl,3,nfumax,*)
      dimension habcdx(nbls1,3,nfumax)
c
         do 32 ifu=1,nfumax1
         do 32 i=1,nbls1
         ijkl=indx(i)
         klpar=idx2(ijkl)
           habcdx(i,1,ifu)=habcd(klpar,1,ifu,lckl)
           habcdx(i,2,ifu)=habcd(klpar,2,ifu,lckl)
           habcdx(i,3,ifu)=habcd(klpar,3,ifu,lckl)
   32    continue
c
      end
c===================================================================
      subroutine conv1der(nbls1,npij,lci,idx1,indx, aa, aax)
      implicit real*8 (a-h,o-z)
      dimension idx1(*),indx(*)
      dimension aa(npij,*)
c output :
      dimension aax(nbls1)
c
      do 10 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        aax(i)=aa(ijpar,lci)*2.0d0
   10 continue
c
      end
c===================================================================
      subroutine isymm1(nbls,isymm)
      dimension isymm(*)
      do 10 ii=1,nbls
      isymm(ii)=1
   10 continue
      return
      end
c===================================================================
c Used when iroute=2 (new) :
c==========================
      subroutine gcpairs(nbls,lcij,ngci1,ngcj1,ngcij, gcij, gcijx)
      implicit real*8 (a-h,o-z)
c------------------------------------------------------
      dimension gcij(ngci1,ngcj1,*)
      dimension gcijx(ngcij,nbls)
c------------------------------------------------------
c      This is called from Erinteg and Erintsp
c             From contraction loops
c
c          FOR GENERAL CONTRACTED SHELLS
c
c------------------------------------------------------
c
      do ijkl=1,nbls
         ijpg=0
         do igc=1,ngci1
            do jgc=1,ngcj1
               ijpg=ijpg+1
               gcijx(ijpg,ijkl)=gcij(igc,jgc,lcij)
            enddo
         enddo
      enddo
c
      end
c====================================================================
      subroutine gcquart(nbls,nbls1, index,
     *                   ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                   gcij,ngcij,  gckl,ngckl,
ccc              output :
     *                   indgc,gcoef)
      implicit real*8 (a-h,o-z)
c------------------------------------------------------
      dimension index(*)
      dimension indgc(nbls)
      dimension gcoef(ngcd,nbls), gcij(ngcij,nbls),gckl(ngckl,nbls)
c------------------------------------------------------
c      This is called from Erinteg and Erintsp
c             From contraction loops
c
c         FOR GENERAL CONTRACTED SHELLS
c------------------------------------------------------
      ijpg=ngci1*ngcj1
      klpg=ngck1*ngcl1
c
      do 204 i=1,nbls1
         ijkl=index(i)
         ijklg=0
         do ijp1=1,ijpg
            gcoefij=gcij(ijp1,ijkl)
            do klp1=1,klpg
               ijklg=ijklg+1
               gcoef(ijklg,ijkl)=gcoefij*gckl(klp1,ijkl)
            enddo
         enddo
         indgc(ijkl)=ijklg
  204 continue
c
      end
c====================================================================
c Used when iroute=1 (old) :
c==========================
      subroutine gcparij(nbls, indxij,npij,
     *                   ii,jj,ngci1,ngcj1,ngcij,
     *                   gci,gcj,gcij)
      implicit real*8 (a-h,o-z)
c--------------------------------------------------------
      dimension indxij(*)
      dimension gci(npij,ngci1,*),gcj(npij,ngcj1,*)
      dimension gcij(ngcij,nbls)
c--------------------------------------------------------
c      This is called from Erinteg and Erintsp
c             From contraction loops
c
c       FOR GENERAL CONTRACTED SHELLS
c--------------------------------------------------------
c
      do 204 ijkl=1,nbls
      ijpar=indxij(ijkl)
             ijpg=0
             do 2041 igc=1,ngci1
             coefi=gci(ijpar,igc,ii)
             do 2041 jgc=1,ngcj1
             coefj=gcj(ijpar,jgc,jj)*coefi
cccccc       if(jcs.eq.ics .and. jgc.eq.igc) coefj=coefj*0.5d0
             ijpg=ijpg+1
             gcij(ijpg,ijkl)=coefj
 2041        continue
  204 continue
c
      end
c====================================================================
      subroutine gcparkl(nbls,indxkl,npkl,
     *                   kk,ll,ngck1,ngcl1,ngckl,
     *                   gck,gcl,gckl)
      implicit real*8 (a-h,o-z)
c******************************************************
      dimension indxkl(*)
      dimension gck(npkl,ngck1,*),gcl(npkl,ngcl1,*)
c
      dimension gckl(ngckl,nbls)
c------------------------------------------------------
c
      do 204 ijkl=1,nbls
      klpar=indxkl(ijkl)
             klpg=0
             do 2042 kgc=1,ngck1
             coefk=gck(klpar,kgc,kk)
             do 2042 lgc=1,ngcl1
             coefl=gcl(klpar,lgc,ll)*coefk
ccccccc      if(lcs.eq.kcs .and. lgc.eq.kgc) coefl=coefl*0.5d0
             klpg=klpg+1
             gckl(klpg,ijkl)=coefl
 2042        continue
c
  204 continue
c
      end
c====================================================================
      subroutine gcqijkl(nbls,nbls1, index,npij,npkl,
     *                   ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                   indgc,gcoef,
     *                   gcij,ngcij, gckl,ngckl)
      implicit real*8 (a-h,o-z)
c------------------------------------------------------
      dimension index(*)
      dimension indgc(nbls)
      dimension gcoef(ngcd,nbls)
      dimension gcij(ngcij,nbls),gckl(ngckl,nbls)
c------------------------------------------------------
c      This is called from Erinteg and Erintsp
c             From contraction loops
c
c   FOR GENERAL CONTRACTED SHELLS
c
      ijpg=ngci1*ngcj1
      klpg=ngck1*ngcl1
c
      do 204 i=1,nbls1
         ijkl=index(i)
         ijklg=0
         do 2043 ijp1=1,ijpg
             gcoefij=gcij(ijp1,ijkl)
             do 2043 klp1=1,klpg
             ijklg=ijklg+1
             gcoef(ijklg,ijkl)=gcoefij*gckl(klp1,ijkl)
 2043    continue
         indgc(ijkl)=ijklg
  204 continue
c
      end
c====================================================================
      subroutine get_max4kl(nbls,npkl,densmax,lckl,ecd,esti2kl)
      implicit real*8 (a-h,o-z)
      dimension densmax(nbls)
      dimension ecd(npkl,lckl)
c output
      dimension esti2kl(lckl)
c
ckw   call absmax(nbls,densmax,ii,ddmax)
c
cccc  call getrval('dmx_over',dmx_overall)
      call getrval('dmx_over',ddmax)
c
      do 100 kl=1,lckl
ccc   call absmax(npkl,ecd(1,kl),ii,esti2max)
      esti2max=0.d0
      do klp=1,npkl
         esti2max=max(ecd(klp,kl),esti2max)
      enddo
      esti2kl(kl)=ddmax*esti2max
  100 continue
c
      end
c====================================================================
      subroutine symm_jump(isymm,npij,npklx,npkl,jump)
      logical jump(npij)
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file
      dimension isymm(*)
c
c not for schwarz :
c
      if(schw4.eq.'sch1') then
         do ijpar=1,npij
            jump(ijpar)=.false.
         enddo
         return
      endif
c
      ijkl=0
      do ijpar=1,npij
         npklend=npklx
         if(npkl.eq.0) npklend=ijpar
c
         jump(ijpar)=.true.
         do klpar=1,npklend
            ijkl=ijkl+1
            if(isymm(ijkl).gt.0) then
               jump(ijpar)=.false.
               ijkl=ijkl+(npklend-klpar)
               go to 100
            endif
         enddo
  100 continue
      enddo
c
      end
c====================================================================
      subroutine symm_fact(isymm,npij,npkl,ndiag,rnsym,symfac)
      implicit real*8 (a-h,o-z)
      character*4 schw4
      common /schwarz/ schw4    ! used only in this file
      dimension isymm(*)
      dimension symfac(*)
c
      if(ndiag.eq.0) then
         ndim=npij*(npij+1)/2
      else
         ndim=npij*npkl
      endif
c
c for Schwarz integrals symmetry is not used but symfac()=1 must be setup
c
      if(schw4.eq.'sch1') then
         do ijkl=1,ndim
            ijklsm=isymm(ijkl)
            if(ijklsm.gt.0) symfac(ijkl)=1.d0
         enddo
         RETURN
      endif
c
c for everything else symfact() is use when symmetry is specified :
c
         do ijkl=1,ndim
            ijklsm=isymm(ijkl)
            if(ijklsm.gt.0) symfac(ijkl)=rnsym*dble(ijklsm)
         enddo
c---------------------------------------------------------
      end
c====================================================================
      subroutine schw_quarts(iis,jjs,ibl,ijbl,nbl2,nijbeg,nijend,
     *                       ipres )
c
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension ipres(*)
c
c It is called for diagonal blocks (ibl=kbl)
c eliminate off-diagonal contracted quartets
c
      kbl=ibl
c
      ijkl=0
      do ijp=nijbeg,nijend
        ijcs=ijbl(ibl,ijp)
        ics=iis(ijcs)
        jcs=jjs(ijcs)
           do klp=nijbeg,ijp
              klcs=ijbl(kbl,klp)
              kcs=iis(klcs)
              lcs=jjs(klcs)
              ijkl=ijkl+1
              if(kcs.eq.ics .and. lcs.eq.jcs) then
ccccccc          ipres(ijkl)=0
              else
                 ipres(ijkl)=0
              endif
           enddo
      enddo
c
      end
c====================================================================
      subroutine make_dens1(densmax,nbls)
c called only for Schwarz integrals
c
      implicit real*8 (a-h,o-z)
      dimension densmax(nbls)
      data xmil /1.d6/
c
      do ijkl=1,nbls
         densmax(ijkl)=xmil
      enddo
c
      call setrval('dmx_over',  xmil     )
c
      end
c====================================================================
      subroutine set_boamax2(aa,bb,l1,l2,boa_max) ! stability in xwpq_22
      implicit real*8 (a-h,o-z)
      dimension aa(*),bb(*)  ! contraction length
c
      a_exp=aa(l1)
      b_exp=bb(l2)
c
      boa_max=2.d0*b_exp/max(1.d0,a_exp)
c
cccc  call setrval('boa_max',boa)
c
      end
c====================================================================
      subroutine set_boamax1(aa,bb,npij,l1,l2,boa_max)
      implicit real*8 (a-h,o-z)
      dimension aa(npij,*),bb(npij,*)  ! contraction length
c
      boa_max=0.d0
      do ijpar=1,npij
         a_exp=aa(ijpar,l1)
         b_exp=bb(ijpar,l2)
         boa=2.d0*b_exp/max(1.d0,a_exp)
         boa_max=max(boa_max,boa)
      enddo
c
cccc  call setrval('boa_max',boa_max)
c
      end
c====================================================================
      subroutine smblock_neg(schwarz,ncs, ijcs1,klcs1,skip_sb)
      implicit real*8 (a-h,o-z)
      logical skip_sb
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension schwarz(ncs,ncs)
c
c first shell-pairs in left- & right-side pair-blocks
c
      call get_ij_half(ijcs1, ics, jcs)
      call get_ij_half(klcs1, kcs, lcs)
c
      eps2=eps*eps
c
      x_ij_ij_max=schwarz(ics,jcs)
      x_kl_kl_max=schwarz(kcs,lcs)
c
      call setrval('x_ij_ij',x_ij_ij_max)
      call setrval('x_kl_kl',x_kl_kl_max)
c
      estim_max=dens_max_el * x_ij_ij_max * x_kl_kl_max
c
      skip_sb=.false.
      if(estim_max.LT.eps2 ) skip_sb=.true.
c
      end
c====================================================================
