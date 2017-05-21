c===============================================================
c
C July 4,97 KW: Changes concerning symmetry have been made in TWO
C               subroutines : precspec_2 and prec4neg_2
c----------------------------------------------------------------
c All routines of the type name_1 are used when   IROUTE=1
c All routines of the type name_2 are used when   IROUTE=2
c----------------------------------------------------------------
c New overlap :
c
c Sab=(a*b)**3/4 * (a+b)**-1 * exp( -ab/(a+b) * Rab**2)
c
c====================================================================
      subroutine prec2ij(ibl, bl,inx,npar,nbl2, iis,jjs,ijbl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstd
      common /route/ iroute
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primij/ iabprim, ijdim ,ijpar1
      common /time0/ tprec2
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*), npar(*)
c---------------------------------------
      if(icheck.gt.0) then
         call getmem(0,l0)
         return
      endif
c---------------------------------------
      call secund(tprec2b)
c---------------------------------------
      nparij=npar(ibl)
c---------------------------------------
      call dimenij(ibl,inx,nparij,nbl2, iis,jjs,ijbl,ijdim,ijcont)
c---------------------------------------
      ijpar1=nparij
      IF( iroute.eq.1 ) THEN
         call getmem(3*ijdim,iabprim)
         call ab_prim_1(ibl,nparij,ijbl,nbl2,
     *                  bl(inuc),bl(ibas),inx,iis,jjs,
     *                  bl(iabprim),ijcont)
      ELSE
         call getmem(2*ijcont+ijdim,iabprim)
         iapb =iabprim
         i1apb=iabprim+ijcont
         isab =iabprim+ijcont*2
c
         call ab_prim_2(ibl,nparij,ijbl,nbl2,
     *                  bl(inuc),bl(ibas),inx,iis,jjs,
     *                  bl(iapb),bl(i1apb),bl(isab),ijcont)
c200
      ENDIF
c---------------------------------------
      call secund(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c=======================
      subroutine prec2kl(ibl, bl,inx,npar,nbl2, iis,jjs,ijbl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstd
      common /route/ iroute
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primkl/ kabprim, kldim ,klpar1
      common /time0/ tprec2
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*), npar(*)
c---------------------------------------
      if(icheck.gt.0) then
         call getmem(0,l0)
         return
      endif
c---------------------------------------
      call secund(tprec2b)
c---------------------------------------
      nparkl=npar(ibl)
c---------------------------------------
      call dimenij(ibl,inx,nparkl,nbl2, iis,jjs,ijbl,kldim,klcont)
c---------------------------------------
      klpar1=nparkl
      IF( iroute.eq.1 ) THEN
         call getmem(3*kldim,kabprim)
         call ab_prim_1(ibl,nparkl,ijbl,nbl2,
     *                  bl(inuc),bl(ibas),inx,iis,jjs,
     *                  bl(kabprim),klcont)
c
      ELSE
         call getmem(2*klcont+kldim,kabprim)
         icpd =kabprim
         i1cpd=kabprim+klcont
         iscd =kabprim+klcont*2
c
         call ab_prim_2(ibl,nparkl,ijbl,nbl2,
     *                  bl(inuc),bl(ibas),inx,iis,jjs,
     *                  bl(icpd),bl(i1cpd),bl(iscd),klcont)
      ENDIF
c---------------------------------------
      call secund(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c=======================
      subroutine dimenij(ibl,inx,nparij,nbl2, iis,jjs,ijbl,ijdim,ijcont)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
c
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
      icont=inx(5,ics1)-inx(1,ics1)
      jcont=inx(5,jcs1)-inx(1,jcs1)
      ijcont=icont*jcont
      ijdim=nparij*ijcont
      end
c====================================================================
      subroutine indexp(npij,npkl,npklx,nijbeg,nijend, nklbeg,nklend,
     *                  indxij,indxkl,ipres,indxp,indxr,nblsp)
      dimension indxij(*),indxkl(*),ipres(*),indxp(*),indxr(*)
c----------------------------------------------------------------
c  setup relation between PRESENT c.s.quartets and pairs
c            ijklp and ijpar=1,..., and klpar=1,...
c----------------------------------------------------------------
c
      ijklp=0
      ijkl=0
      ijpar=0
      do ijp=nijbeg,nijend
         ijpar=ijpar+1
         klpar=0
         nklendx=nklend
         if(npklx.eq.0) nklendx=ijp
         do klp=nklbeg,nklendx
            klpar=klpar+1
            ijkl=ijkl+1
            indxr(ijkl )=0
            if(ipres(ijkl).ne.0) then
               ijklp=ijklp+1
               indxij(ijklp)=ijpar
               indxkl(ijklp)=klpar
                indxp(ijklp)=ijkl
                indxr(ijkl )=ijklp
            endif
         enddo
      enddo
c
      nblsp=ijklp
c
c         write(8,88) (indxij(ii),ii=1,nblsp)
c         write(8,88) (indxkl(ii),ii=1,nblsp)
c 88      format('from index4 ;indxij,kl=',5i3)
c         write(8,89) (indxp(ii),ii=1,nblsp)
c 89      format('from index4 ;indxp    =',5i3)
c
      end
c====================================================================
      subroutine specasg(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                   buf,buf1, const,rysx,xpqr,txxr,
     *                   ngcd,indgc,gcoef,ijkln)
c--------------------------------------------------------------------
c Works with Transposed BUF array
c--------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
c--------------------------------------------------------------------
c
      dimension bl(*)
      dimension index(*),indxij(*),indxkl(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension const(*),rysx(*)
      dimension buf(ngcd,nbls,ijkln),buf1(nbls1,*)
      dimension indgc(nbls)
      dimension gcoef(ngcd,nbls)
c
c--------------------------------------------------------------------
c     FOR GENERAL CONTRACTED SHELLS
c  this subroutine constitues the special code for
c  two types of integrals over nbls quartets of primitive shells
c  1. (ss|ss)
c  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p
c  these integrals are also contracted here.
c
c  input
c  ------
c  all precalculated values for whole block :
c
c  const(nbls) - contains consts=pi3*sabcd/(pq*sqrt(ppq)) for all int.
c  rysx(nbls) - contains  xrys i.e. arguments for fm,ft0 routines
c  xp,xp      - geometry for p,q
c
c  output
c  ------
c  buf(ngcd,nbls,ijkln) - contains final integrals
c--------------------------------------------------------------------
c
c memory for f00,f11 :
c
      call getmem(nbls1,if00)
      call getmem(nbls1,if11)
c
      if00=if00-1
      if11=if11-1
c
      do 100 i=1,nbls1
      xrys=rysx(i)
      call ft0
      bl(if00+i)=f0
      bl(if11+i)=f1
  100 continue
c
c  special code for (ss ss) integrals
c
      if(mmax.eq.1) then
          do 2031 i=1,nbls1
          buf1(i,1)=const(i)*bl(if00+i)
 2031     continue
      go to 203
      endif
c
c  special code for (ps ss), (sp ss) (ss ps) and (ss sp)
c
cxxxx if(mmax.eq.2) then
      if (ityp.gt.1) then
        do 101 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 101    continue
      else if (jtyp.gt.1) then
        do 102 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(2,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 102    continue
      else if (ktyp.gt.1) then
        do 103 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 103    continue
      else
        do 104 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 104    continue
      endif
c
        do 106 i=1,nbls1
          buf1(i,1)=-xpqr(1,i)*const(i)
          buf1(i,2)=-xpqr(2,i)*const(i)
          buf1(i,3)=-xpqr(3,i)*const(i)
106     continue
c
c
  203 continue
c
      if(first) then
c
           do 204 icx=1,lnijkl
           do 204 i=1,nbls1
           xint=buf1(i,icx)
             ijkl=index(i)
ckw          ngcq=indgc(ijkl)
             do 2041 iqu=1,ngcd
             buf(iqu,ijkl,icx)=xint*gcoef(iqu,ijkl)
 2041        continue
  204    continue
c
         first=.false.
      else
c
           do 205 icx=1,lnijkl
           do 205 i=1,nbls1
           xint=buf1(i,icx)
             ijkl=index(i)
ckw          ngcq=indgc(ijkl)
             do 2051 iqu=1,ngcd
             buf(iqu,ijkl,icx)=buf(iqu,ijkl,icx)+xint*gcoef(iqu,ijkl)
 2051        continue
  205    continue
c
      endif
c
c release memory
c
      call retmem(2)
c
      return
      end
c====================================================================
c1999 subroutine reuseij(ibl,nbls,mmax,nblok1,reuse)
      subroutine reuseij(ibl,mmax,ijcurr_f,ijcurr_l,reuse)
      character*3 reuse
      logical cond1,cond2
      common /ijcsfl/ ijblokp,ijprevf,ijprevl,ijtprev,maxprev,ngcprev
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c1999 dimension nblok1(2,*)
c-----------------------------------------------------------
      reuse='no '
c---------------------------------------------------------------
c new in order to remember left-hand pairs IJ in a block ikbl
c keep it in common /ijcsfl/ and use it
c
c1999 ijcurrf=nblok1(1,1)
c1999 ijcurrl=nblok1(1,nbls)
      ijcurrf=ijcurr_f
      ijcurrl=ijcurr_l
      ijtcurr=0
      if(mmax.eq.2) then
        if(itype1.gt.1 .or. jtype1.gt.1) ijtcurr=2
      endif
      maxcurr=mmax
c     ngccurr=ngcd
      ngccurr=1000*ngci1+100*ngcj1+10*ngck1+ngcl1
c----
      if(ibl.ne.ijblokp) go to 1000
c----------------------------------------
      cond1=.false.
      cond2=.false.
c-----
        if(ijcurrf.eq.ijprevf .and. ijcurrl.eq.ijprevl) cond1=.true.
ctest.............
c       if(ngccurr.eq.ngcprev. or.
c    *  (ngccurr.gt.1.and.ngcprev.gt.1)) cond2=.true.
ctest.............
        if(ngccurr.eq.ngcprev ) cond2=.true.
c----------------------------------------
c  check if the current block has the same left-hand pairs IJ
c  like the previous one. If it has DO NOT precalculate
c  corresponding IJ-pair quantities. Instead, re-use old ones
c  ( the memory location is the same )
c
c ijtprev and ijtcurr may have value of 0 or 2. They are needed only
c for the mmax=2 and give information about re-scaling (2) or not (0)
c precalculated TXAB
c
      if(cond1 .and. cond2) then
c------>  if(maxprev.ne.2) go to 9999
          if(maxcurr.eq.maxprev) go to 9999
          if(ijtprev.eq.2 .and. ijtcurr.ne.2) then
          else
            go to 9999
          endif
      endif
          go to 8888
 9999 continue
      reuse='yes'
 8888 continue
c
c----------------------------------------
cnew in order to remember left-hand pairs IJ in a block ikbl
c keep it in common /ijcsfl/ and use it in precalc...
c
 1000 continue
      ijblokp=ibl
      ijprevf=ijcurrf
      ijprevl=ijcurrl
      ijtprev=ijtcurr
      maxprev=maxcurr
      ngcprev=ngccurr
c----------------------------------------
      end
c================================================================
      subroutine precdiag
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /route/ iroute
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c---------------------------------------------------------------
c* for a diagonal case : npklx=0
c
c* since this is for a diagonal block
c* pairs IJ and KL are the same
c
      ixcd=ixab
      ixq =ixp
      ixqn=ixpn
      ixqq=ixpp
c
      icpd=iapb
      i1cpd=i1apb
      ickl=icij
      ifkl=ifij
      iscd=isab
c
      icc=iaa
      idd=ibb
c
      icks=icis
      icls=icjs
c
      iecd=ieab
c
      itxcd=itxab
c
      icdnia=iabnia
c-----------------------------------
      IF( iroute.eq.1 ) THEN
        igck=igci
        igcl=igcj
      ELSE
        igckl=igcij
      ENDIF
c-----------------------------------
      end
c===============================================================
c======= Duplicated routines for different values of IROUTE ====
c============================================================
c for iroute=1 :
c============================================================
      subroutine ab_prim_1(ibl,nparij,ijbl,nbl2,
     *                   datnuc,datbas,inx,iis,jjs,
     *                   abprim,ijcont)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension abprim(nparij,ijcont,3)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c---------------------------------------------------------------
      do 100 ijpar=1,nparij
      ijcs=ijbl(ibl,ijpar)
      ics=iis(ijcs)
      jcs=jjs(ijcs)
c
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      iat=inx(2,ics)
        ja=inx(1,jcs)+1
        je=inx(5,jcs)
        jat=inx(2,jcs)
        xab=datnuc(2,iat)-datnuc(2,jat)
        yab=datnuc(3,iat)-datnuc(3,jat)
        zab=datnuc(4,iat)-datnuc(4,jat)
        rr=xab*xab+yab*yab+zab*zab
        ij=0
           do 200 is=ia,ie
           aa=datbas(1,is)
c----------------------------------
             do 200 js=ja,je
             bb=datbas(1,js)
             ij=ij+1
c----------------------------------
c
             axb=aa*bb
             apb=aa+bb
             apb1=one/apb
             e=axb*apb1
c98
             err=e*rr
             if(err.gt.32) then
                exp_err=0.d0         ! exp(-32)=1.27*10**-14
             else
                exp_err=exp(-err)
             endif
c98
             abprim(ijpar,ij,1)=apb
             abprim(ijpar,ij,2)=apb1
c98          abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp(-e*rr)
             abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp_err
c----------------------------------
  200     continue
  100 continue
      end
c============================================================
      subroutine precalc2_1(isupb,bl,mmax,mmax1,nhabcd,nfumax,nbl2,nbls,
     *                       inx,iis,jjs,ijbl,nblok1,
     *                       ibl,nijbeg,nijend,npij,
     *                       kbl,nklbeg,nklend,npklx,npkl)
      implicit real*8 (a-h,o-z)
      character*3 reuse
c-----------------------------------------------------------
      common /neglect/ eps,eps1,epsr
c
      common /begin/ ijbegin,klbegin
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /ilepar/ lpartot,lpareal
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nblok1(2,*)
c---------------------------------------------------------------
      par268=eps1*8.d0
c---------------------------------------------------------------
      ijcurr_f=ijbl(ibl,nijbeg)
      ijcurr_l=ijbl(ibl,nijend)
      call reuseij(ibl,mmax,ijcurr_f,ijcurr_l,reuse)
c-----------------------------------------------------------
      ijbegin=nijbeg
      klbegin=nklbeg
c-----------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c--------
      lpartot=lpartot+1
      if(reuse.eq.'no ') then
          lpareal=lpareal+1
          call precal2a_1(bl(ibas),bl(inuc),iis,jjs,inx, npij,
     *                  ibl,ijbl,nbl2,nijbeg,nijend,
     *                  bl(iabprim      ),ijpar1 ,
     *                  lcij,bl(iaa),bl(ibb),bl(ieab),bl(icis),bl(icjs),
     *                  bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *                  bl(igci),bl(igcj),ngci1,ngcj1,'left ',par268)
c----------
c for abnia
          if(mmax.gt.2) then
             call precal2b_1(mmax1,lcij,npij, bl(i1apb),
     *                       ijpar1,ijbegin, bl(iabnia))
          endif
      endif
c----------------------------------------
c precalculations for pairs kl
c
      if(kbl.eq.ibl) then
          kabprim=iabprim
          kldim  =ijdim
          klpar1 =ijpar1
      endif
c
      if(npkl.ne.0) then
         call precal2x_1(kabprim,kldim,icpd,i1cpd,iscd)
c
         lpartot=lpartot+1
         lpareal=lpareal+1
         call precal2a_1(bl(ibas),bl(inuc),iis,jjs,inx, npkl,
     *                  kbl,ijbl,nbl2,nklbeg,nklend,
     *                  bl(kabprim      ),klpar1 ,
     *                  lckl,bl(icc),bl(idd),bl(iecd),bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igck),bl(igcl),ngck1,ngcl1,'right',par268)
      endif
c----------------------------------------
c for cdnia
c
         if(npkl.ne.0) then
            if(mmax.gt.2) then
               call precal2b_1(mmax1,lckl,npklx,bl(i1cpd),
     *                         klpar1,klbegin, bl(icdnia))
            endif
         else
            call precdiag
         endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c_1(npklx,bl(i1cpd),klpar1,lckl,
     *                   bl(ihabcd),nhabcd,nfumax )
      endif
c---------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_1(datbas,datnuc,iis,jjs,inx,npij,
     * ibl,ijbl,nbl2,nijbeg,nijend,
     * abprim,ijpar1,lcij, aaa,bbb,estab,cis,cjs,
     * xab,xparij,coefij,factij,txab, gci,gcj,ngci1,ngcj1,which,par268)
c
      implicit real*8 (a-h,o-z)
      character*5 which
      COMMON /types/itype,jtype,ktype,ltype,ityp,jtyp,ktyp,ltyp
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /neglect/ eps,eps1,epsr
c
      dimension abprim(ijpar1,lcij,3)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension aaa(npij,*),bbb(npij,*), cis(npij,*),cjs(npij,*)
      dimension xab(npij,3), xparij(npij,3,lcij,3)
      dimension estab(npij,lcij)
      dimension coefij(npij,lcij), factij(npij,lcij)
      dimension xa(3),xb(3),txab(npij,3,*)
      dimension gci(npij,ngci1,*),gcj(npij,ngcj1,*)
c---------------------------------------------------------------
c     par268=eps1*8.d0
c---------------------------------------------------------------
c 2002 : stability issue : find max 2b/a and max(1/(a+b)
c 2002 : stability issue : find max 2d/c and max(1/(c+d)
c---------------------------------------------------------------
c
      if(which.eq.'left ') then
          itypp=ityp
          jtypp=jtyp
          ngcii=ngci
          ngcjj=ngcj
          nqii=nqi
          nqjj=nqj
          lcii=lci
          lcjj=lcj
      else
          itypp=ktyp
          jtypp=ltyp
          ngcii=ngck
          ngcjj=ngcl
          nqii=nqk
          nqjj=nql
          lcii=lck
          lcjj=lcl
      endif
c
c-------------------------------------------------------------
c for primitive integral neglect : to prevent from throwing away
c small primitive contributions which might sum up to a big value
c
c98   contrij=dble(lcii*lcjj)
c98   contrij=sqrt(contrij)
c-------------------------------------------------------------
c for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c---------------
c
c precalculations for the pairs IJ :
c
c2002 instablity issue :
      boamax=0.d0
      rapbmax=0.d0
c
      ijpar=0
      do 100 ijp=nijbeg,nijend
      ijpar=ijpar+1
        ijcs=ijbl(ibl,ijp)
c
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        fact1=one
        if(ics.eq.jcs) fact1=fact1*half
c
c      starting contr
        ia=inx(1,ics)+1
        ja=inx(1,jcs)+1
c      last contr
        ie=inx(5,ics)
        je=inx(5,jcs)
c
c       number of general contr.
c       ngci=inx(4,ics)
c       ngcj=inx(4,jcs)
c in the common block contr
c
        iatom=inx(2,ics)
        jatom=inx(2,jcs)
c
c97
c
c center23 (atom number 0 is the one with s-orbital/zero-exp.)
c
        if(iatom.eq.0) iatom=jatom
        if(jatom.eq.0) jatom=iatom
c
        do 150 i=1,3
          if(iatom.eq.0 .and. jatom.eq.0) then
            xa(i)=zero
            xb(i)=zero
          else
            xa(i)=datnuc(1+i,iatom)
            xb(i)=datnuc(1+i,jatom)
          endif
          xab(ijpar,i)=xa(i)-xb(i)
  150   continue
c
c97
c98
ccc      factab_max=0.d0
c98
         ji=0
         is1=0
         do 200 is=ia,ie
            is1=is1+1
            aa=datbas(1,is)
c2002
            aor1=1.d0/max(1.d0,aa)
            aaa(ijpar,is1)=aa
               if(ngctot.eq.0) then
                  csi=datbas(2,is)
                  cpi=datbas(3,is)
                  coefi=csi
ckw99n
                  est_i=max( abs(csi),abs(cpi) )
ckw99n
                  if(itypp.eq.3) then
                     coefi=cpi
                     csi=csi/cpi
                     facti=csi
                  endif
                  coefi=coefi*fact1
ckw99
ccccc????         est_i=coefi
ckw99
                  if(which.eq.'right') coefi=coefi*par268
                  cis(ijpar,is1)=csi
               else
c gen.contr. shell is somewhere
ckw99
                  est_i=one
ckw99
                  do 210 ig=0,ngcii
                  gci(ijpar,ig+1,is1)=datbas(ig+2,is)
  210             continue
               endif
c
         js1=0
         do 200 js=ja,je
            js1=js1+1
            ji=ji+1
            bb=datbas(1,js)
c2002
            boamax=max(boamax,bb*aor1)      ! b/max(1,a)
            bbb(ijpar,js1)=bb
               if(ngctot.eq.0) then
                  csj=datbas(2,js)
                  cpj=datbas(3,js)
                  coefj=csj
                  if(jtypp.eq.3) coefj=cpj
                  coefij(ijpar,ji )=coefi*coefj
ckw99n
                  est_j=max( abs(csj),abs(cpj) )
c
                  estab(ijpar,ji)=est_i*est_j
ckw99n
ckw99
cvc?????          estab(ijpar,ji)=est_i*coefj
ckw99
c
                  if(jtypp.eq.3) then
                     csj=csj/cpj
                     factj=csj
                     if(itypp.eq.3) then
                       factij(ijpar,ji  )=facti*factj
                     endif
                  endif
                     cjs(ijpar,js1)=csj
               else
c gen.contr.
                  coefij(ijpar,ji )=fact1
ckw99
                  estab(ijpar,ji)=coefij(ijpar,ji)
ckw99
                  if(which.eq.'right') coefij(ijpar,ji)=fact1*par268
c
                  do 220 jg=0,ngcjj
                  gcj(ijpar,jg+1,js1)=datbas(jg+2,js)
  220             continue
               endif
c--------------------------------
ckw99
c              factab=abs(estab(ijpar,ji))
c              if(factab.gt.factab_max) factab_max=factab
c--------------------------------
c
            rapb=abprim(ijp,ji,2)
            sab =abprim(ijp,ji,3)
c
            aa1=aa*abprim(ijp,ji,2)
            bb1=bb*abprim(ijp,ji,2)
c
            coefij(ijpar,ji)=coefij(ijpar,ji)     *sab
ckw99
            estab(ijpar,ji)=estab(ijpar,ji)*sab
ckw99
c2002
            rapbmax=max(rapbmax,rapb)
c
            xpn_max=1.d0
            do 230 l=1,3
            xparij(ijpar,l,ji,1)=aa1*xa(l) + bb1*xb(l)  ! xp(ijpar,l,ji
                xxl=xa(l)
                if(nqii.lt.nqjj) xxl=xb(l)
            xparij(ijpar,l,ji,2)=xparij(ijpar,l,ji,1)-xxl  ! xpn
            xparij(ijpar,l,ji,3)=aa*xa(l)+bb*xb(l)         ! xpp
ckw99
c for neglect :
            xpn_abs=abs( xparij(ijpar,l,ji,2) )
            if(xpn_abs.gt.xpn_max) xpn_max=xpn_abs
ckw99
c center23:
            if(aa.le.zero .or. bb.le.zero) then
               xparij(ijpar,l,ji,2)=zero
            endif
  230       continue
c------------------------------------------------------------------
ckw99
            estab(ijpar,ji)=estab(ijpar,ji)*xpn_max
            estab(ijpar,ji)=estab(ijpar,ji)*estab(ijpar,ji)
ckw99
c------------------------------------------------------------------
  200    continue       ! end of the loop over contractions
c
c
c.............................................................
c        factab_norm=1.d0/factab_max
c        factab_norm= factab_norm * factab_norm
c        do ijprim=1,lcii*lcjj
c           estab(ijpar,ijprim)=estab(ijpar,ijprim)*factab_norm
c        enddo
c.............................................................
c------------------------------------------------------------------
  100 continue
c2002
      boamax=2.d0*boamax
      if(which.eq.'left ') then
         call setrval('boamax',boamax)
         call setrval('rapbmax',rapbmax)
         call setrval('docmax',boamax)
         call setrval('rcpdmax',rapbmax)
      else
         call setrval('docmax',boamax)
         call setrval('rcpdmax',rapbmax)
      endif
c2002
c
ctxab
c
      if(nqii.ge.nqjj) then
         ijs1=0
         do 151 is1=1,lcii
         do 151 js1=1,lcjj
         ijs1=ijs1+1
         do 151 ijpar=1,npij
            txab(ijpar,1,ijs1)=-bbb(ijpar,js1)*xab(ijpar,1)
            txab(ijpar,2,ijs1)=-bbb(ijpar,js1)*xab(ijpar,2)
            txab(ijpar,3,ijs1)=-bbb(ijpar,js1)*xab(ijpar,3)
c center23:
            aa=aaa(ijpar,is1)
            if(aa.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  151    continue
      else
         ijs1=0
         do 152 is1=1,lcii
         do 152 js1=1,lcjj
         ijs1=ijs1+1
         do 152 ijpar=1,npij
            txab(ijpar,1,ijs1)= aaa(ijpar,is1)*xab(ijpar,1)
            txab(ijpar,2,ijs1)= aaa(ijpar,is1)*xab(ijpar,2)
            txab(ijpar,3,ijs1)= aaa(ijpar,is1)*xab(ijpar,3)
c center23:
            bb=bbb(ijpar,js1)
            if(bb.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  152    continue
      endif
c
c for the case with mmax=2 - special cases
c
      if ( mmax.eq.2 ) then
         if (itypp.gt.1 .or. jtypp.gt.1 ) then
           do 153 ijs1=1,lcij
           do 153 ijpar=1,npij
            rapb1=abprim(nijbeg-1+ijpar,ijs1,2)
c
            txab(ijpar,1,ijs1)=txab(ijpar,1,ijs1)*rapb1
            txab(ijpar,2,ijs1)=txab(ijpar,2,ijs1)*rapb1
            txab(ijpar,3,ijs1)=txab(ijpar,3,ijs1)*rapb1
  153      continue
         endif
      endif
c--------------
      end
c====================================================================
      subroutine precal2b_1(mmax1,lcij,npij, rapb,ijpar1,ijbegin, abnia)
c-------------------------------------------------------------------
c  Description :
c
c  This subroutine calculates 1 quantity for contracted
c  shell pairs which belong to a given block of contracted
c  of contracted quartets.
c    This subroutine is called (from Blockint) only for blocks
c  with total angular momentum MMAX > 2 ,where MMAX is the sum
c  of angular momentum of four contracted shells (+1) /for exemple :
c  (ss ss):mmax=1, (pp,ps):mmax=4, (dd,ll):mmax=7, (ff,ff):mmax=13 /
c
c  INPUT
c  ------
c  MMAX1 = MMAX-1 , where MMAX is the sum of angular momentum
c  NPIJ  - number of pairs <IJ in a given block of contracted quartets
c  RAPB(ij,lcij) - 1/(a+b) : reversed sum of exponents of each
c                     pair of primitive shells from IJth pair of
c                     contracted shells.
c  OUTPUT
c  -------
c  For each pair of primitive shells ij (kl) from IJ (KL) contracted
c  pair belonging to a given block of quartets with total ang.mom. MMAX
c
c  1. ABNIA(IJ,L,ij) -   L*( 0.5/(a+b) )  with L=1,2,...MMAX-1
c
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---> dimension rapb(npij,*), abnia(npij,mmax1,*)
      dimension rapb(ijpar1,*), abnia(npij,mmax1,*)
      data half /0.5d0/
c--------------------------
      ijstart=ijbegin-1
c--------------------------
c
      do 150 ij=1,lcij
         do 200 ijpar=1,npij
            apb2=half*rapb(ijstart+ijpar,ij)
            abnia(ijpar,1,ij)=apb2
               do 250 i=2,mmax1
               abnia(ijpar,i,ij)=abnia(ijpar,i-1,ij)+apb2
  250          continue
  200    continue
  150  continue
c
      end
c====================================================================
      subroutine precal2c_1(npkl,rcpd,klpar1,lckl,habcd,nhabcd,nfumax)
c---------------------------------------------------------------------
c This subroutine is called for blocks with total angular momentum
c MMAX > 2 (MMAX is the sum of angular momentum of four c.s.+1
c for exemple :(ss ss):mmax=1, (pp,ps):mmax=4, (dd,ll):mmax=7 etc.)
c
c INPUT
c---------
c
c NPKL  - number of pairs KL> in a given block of contracted quartets
c RCPD(kl,lckl) - 1/(c+d)
c
c OUTPUT
c  ---------
c
c HABCD(KL,lx,ifu,kl)
c
c where IFU denotes number of elementary function from 1 up to the
c total number of functions corresponding to the MMAX-1 /for example
c for mmax=3 it is from 1 to 10 - s,x,y,z,xx,yy,zz,xy,xz,yz /
c The second index lx stays for x,y or z and is used to find the power
c of an elementary function in these directions from matrix HNIA which
c is constat and is set up in BLOCK DATA logobsa. For example :
c
c habcd(klpar,1=x,ifu,kl)=hnia(1,ifu)*rcpd(klpar,kl)
c
c HABCD is used only in TRACY's recursive in subroutines TRACIJ, TRACKL
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /begin/ ijbegin,klbegin
      common /logicd/ hnia(3,1)
c
      dimension rcpd(klpar1,*)
      dimension habcd(nhabcd,3,nfumax,*)
c-----------------------
      klstart=klbegin-1
c-----------------------
      do kl=1,lckl
         do ifu=1,nfumax
            hx=hnia(1,ifu)
            hy=hnia(2,ifu)
            hz=hnia(3,ifu)
            do klpar=1,npkl
               habcd(klpar,1,ifu,kl)=hx*rcpd(klstart+klpar,kl)
               habcd(klpar,2,ifu,kl)=hy*rcpd(klstart+klpar,kl)
               habcd(klpar,3,ifu,kl)=hz*rcpd(klstart+klpar,kl)
            enddo
         enddo
      enddo
c
      end
c====================================================================
      subroutine prec4neg_1(nbls,npij,npkl,ndiag,ij,kl,
     1     ijpar1,lc12, klpar1,lc34,indxij,indxkl,
     2     estab,estcd,densmax,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, nsym,isymm,
coutput
     *     rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c----------------------------------------------------------------
c*
      implicit real*8 (a-h,o-z)
      logical jump(*)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
cflops
      common /flops/ iflop(20)
      common /neglect/ eps,eps1,epsr
      common /begin/ ijbegin,klbegin
c
      dimension indxij(*),indxkl(*),index(*)
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension  densmax(*)
c-------
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension        rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c symmetry
      dimension isymm(*)
      dimension symfac(*)
c
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
c---------------------------------------------------------------
c for test print only :pi256=(pi/256)**-1
c     pi256=81.48733086d0
c---------------------------------------------------------------
      par268=eps1*8.d0
      permut=half
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
      sqrpold=one
      rpold=one
      abpcd1r=zero
      abpcdr=zero
c-----------------------------------------
c
      ijkl=0
      ijkl1=0
      ijklp=0
      do 100 ijpar=1,npij
         coef1=coefij(ijpar,ij)
         npx=npkl
         if(ndiag.eq.0) then
            npx=ijpar
            coef1=coef1*par268
         endif
         if(nsym.gt.0) then
            if(jump(ijpar)) then
               ijkl=ijkl+npx
               go to 100
            endif
         endif
ckw99
         apb1=apb(ijstart+ijpar,ij)
         esti1=estab (ijpar,ij)
         xp1= xp(ijpar,1,ij)
         xp2= xp(ijpar,2,ij)
         xp3= xp(ijpar,3,ij)
ckw99
      do 150 klpar=1,npx
      ijkl=ijkl+1
      ijklsm=isymm(ijkl)
      if(ijklsm.eq.0) go to 150
ckw99    if(nsym.gt.0)  symfac=rnsym*dble(ijklsm)
         ijklp=ijklp+1
         cpd1=cpd(klstart+klpar,kl)
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
         abpcd1=apb1+cpd1
         abxcd=apb1*cpd1
c
      if(abpcd1.ne.abpcd1r) then
         abpcd1r=abpcd1
         abpcdr=one/abpcd1r
      endif
         rho1=abxcd*abpcdr
c
         estim=esti1*esti2
         estim=estim*abpcdr
         estim=estim*densmax(ijkl)
c
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            ijkl1=ijkl1+1
c------->   index(ijkl1)=ijkl
            index(ijkl1)=ijklp
c
            rppq(ijkl1)=abpcdr
            rhoapb(ijkl1)=abpcdr*cpd1
            rhocpd(ijkl1)=abpcdr*apb1
c
            x1= xp1 - xq(klpar,1,kl)
            x2= xp2 - xq(klpar,2,kl)
            x3= xp3 - xq(klpar,3,kl)
c
            rr2=x1*x1 + x2*x2 + x3*x3
            rys(ijkl1)=rr2*rho1
c
              if(abpcdr.ne.rpold) then
                rpold=abpcdr
                sqrpold=sqrt(rpold)
              endif
c
            const(ijkl1)=coef12*sqrpold
            if(ndiag.eq.0.and.ijpar.eq.klpar) then
               const(ijkl1)=const(ijkl1)*permut
            endif
ckw99       if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac
            if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac(ijkl)
         endif
 150  continue
 100  continue
c
      nbls1=ijkl1
c
cflops
cxxxxx iflop(6)=iflop(6)+3*nbls+2*nbls1
c
      end
c====================================================================
      subroutine precspec_1(nbls,npij,npkl,ndiag, ij,kl,
     1     ijpar1,lc12, klpar1,lc34, indxij,indxkl,
     2     estab,estcd,densmax,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, rapb,rcpd,txab,txcd,
     *                   nsym,isymm,
c output
     *     rys,const,xpqr,txxr,nbls1,index)
c---------------------------------------------------------------
c*
      implicit real*8 (a-h,o-z)
      logical jump(*)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
cflops
      common /flops/ iflop(20)
      common /neglect/ eps,eps1,epsr
      common /begin/ ijbegin,klbegin
c
      dimension indxij(*),indxkl(*),index(*)
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension rapb(ijpar1,*),rcpd(klpar1,*)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension  densmax(*)
c-------
      dimension txab(npij,3,*),txcd(npkl,3,*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension        rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
c
      dimension isymm(*)
      dimension  symfac(*)
c
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
c---------------------------------------------------------------
c for test print only :pi256=(pi/256)**-1
c     pi256=81.48733086d0
c---------------------------------------------------------------
      par268=eps1*8.d0
      permut=half
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
      sqrpold=one
      rpold=one
c
      abpcd1r=zero
      abpcdr=zero
c
      ijkl=0
      ijkl1=0
      ijklp=0
      do 100 ijpar=1,npij
         coef1=coefij(ijpar,ij)
         npx=npkl
         if(ndiag.eq.0) then
            npx=ijpar
            coef1=coef1*par268
         endif
         if(nsym.gt.0) then
            if(jump(ijpar)) then
               ijkl=ijkl+npx
               go to 100
            endif
         endif
ckw99
         apb1=apb(ijstart+ijpar,ij)
         esti1=estab (ijpar,ij)
         xp1= xp(ijpar,1,ij)
         xp2= xp(ijpar,2,ij)
         xp3= xp(ijpar,3,ij)
ckw99
      do 150 klpar=1,npx
      ijkl=ijkl+1
c
      ijklsm=isymm(ijkl)
      if(ijklsm.eq.0) go to 150
ckw99    if(nsym.gt.0) symfac=rnsym*dble(ijklsm)
c-----------
         ijklp=ijklp+1
c
         cpd1=cpd(klstart+klpar,kl)
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
         abpcd1=apb1+cpd1
         abxcd=apb1*cpd1
cnew
c
      if(abpcd1.ne.abpcd1r) then
         abpcd1r=abpcd1
         abpcdr=one/abpcd1r
      endif
         rho1=abxcd*abpcdr
c
         estim=esti1*esti2
         estim=estim*abpcdr
         estim=estim*densmax(ijkl)
cTTTT
c        write(8,66) esti1,esti2,estim
c        write(8,67) estim*pi256*eps1    , eps
c        if(estim.gt.epsr) then
c           write(8,*) ' integrals retained '
c        else
c           write(8,*) ' integrals neglected'
c        endif
c 66  format(' estij=',e15.5,' estkl=',e15.5,' estim=',e15.5)
c 67  format(' est. uesd=',e15.5,' esp=',e15.5)
cTTTT
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            ijkl1=ijkl1+1
c---->      index(ijkl1)=ijkl
            index(ijkl1)=ijklp
c
            x1= xp1 - xq(klpar,1,kl)
            x2= xp2 - xq(klpar,2,kl)
            x3= xp3 - xq(klpar,3,kl)
c
            xpqr(1,ijkl1)=x1*rho1
            xpqr(2,ijkl1)=x2*rho1
            xpqr(3,ijkl1)=x3*rho1
c
            rr2=x1*x1 + x2*x2 + x3*x3
            rys(ijkl1)=rr2*rho1
c
              if(abpcdr.ne.rpold) then
                rpold=abpcdr
                sqrpold=sqrt(rpold)
              endif
c
            const(ijkl1)=coef12*sqrpold
            if(ndiag.eq.0.and.ijpar.eq.klpar) then
               const(ijkl1)=const(ijkl1)*permut
            endif
ckw99       if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac
            if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac(ijkl)
c
         endif
 150  continue
 100  continue
c
      nbls1=ijkl1
      if(nbls1.eq.0) return
c
cnew
      if (ityp.gt.1 .or. jtyp.gt.1 ) then
        do 210 i=1,nbls1
        ijkl=index(i)
        ijpar=indxij(ijkl)
c
        xpqr(1,i)=xpqr(1,i)*rapb(ijstart+ijpar,ij)
        xpqr(2,i)=xpqr(2,i)*rapb(ijstart+ijpar,ij)
        xpqr(3,i)=xpqr(3,i)*rapb(ijstart+ijpar,ij)
c
        txxr(1,i)=txab(ijpar,1,ij)
        txxr(2,i)=txab(ijpar,2,ij)
        txxr(3,i)=txab(ijpar,3,ij)
  210   continue
      endif
c
      if (ktyp.gt.1 .or. ltyp.gt.1 ) then
        do 220 i=1,nbls1
        ijkl=index(i)
        klpar=indxkl(ijkl)
        xpqr(1,i)=xpqr(1,i)*rcpd(klstart+klpar,kl)
        xpqr(2,i)=xpqr(2,i)*rcpd(klstart+klpar,kl)
        xpqr(3,i)=xpqr(3,i)*rcpd(klstart+klpar,kl)
c
        txxr(1,i)=txcd(klpar,1,kl)
        txxr(2,i)=txcd(klpar,2,kl)
        txxr(3,i)=txcd(klpar,3,kl)
  220   continue
      endif
cflops
cxxxxx iflop(6)=iflop(6)+3*nbls+2*nbls1
c
      end
c====================================================================
      subroutine xwpq_1(nbls1,xwp,xwq,p1234, ijpar1,lc12, klpar1,lc34,
     *                  lcij,lckl,npij,npkl, indxij,indxkl,index,
     *                  rppq,xp,xq,xpp,xqq,
     *                  txab,txcd,abcd,apb,rcpd,cpd,rapb)
c-----------------------------------------------------------------
c  For two primitive pairs lcij and lckl which contribute
c  to all NBLS1 contracted quartets IJKL in a given block
c  4 quatities are calculated for each quartet (xwp,xwq,p1234,abcd)
c
c  INPUT :
c--------
c
c  NBLS1 - blocksize
c  lcij, lckl - pairs of primitive shells
c  NPIJ,NPKL  - number of contracted pairs
c  INDEX(ijkl1) - for ijkl1 quartet in Reduced block shows
c                 the ijkl  quartet in Original block
c  INDXIJ(ijkl) - shows the contracted pair IJPAR for ijkl quartet
c  INDXKL(ijkl) - shows the contracted pair KLPAR for ijkl quartet
c  RPPQ(ijkl1)  = 1/(a+b+c+d)         where a,b,c,d, are exponents
c  XP(ijpar,3,lcij) - coordinates of P
c  XQ(klpar,3,lckl) - coordinates of Q
c  XPP(ijpar,3,lcij)- coordinates of ( a*A + b*B)
c  XQQ(klpar,3,lckl)- coordinates of ( c*C + d*D)
c  TXAB(ijpar,3,  )  - coordinates XAB re-scaled by exponents a or b
c  TXCD(klpar,3,  )  - coordinates XAB re-scaled by exponents c or d
c  APB(ijpar,lcij) =a+b
c  CPD(klpar,lckl) =c+d
c  RAPB(ijpar,lcij)=1/(a+b)
c  RCPD(klpar,lckl)=1/(c+d)
c
c  OUTPUT
c--------
c  XWP(ijkl1) - coordinates of W-P,  where W=(XPP+XQQ)/(a+b+c+d)
c  XWQ(ijkl1) - coordinates of W-Q
c  P1234(ijkl1)=(txab+txcd)*(rcpd OR rapb)
c  ABCD(ijkl1)= apb*rcpd   OR  cpd*rapb
c-----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /begin/ ijbegin,klbegin
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /stabili/ nstable,nostable
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
c
      dimension indxij(*),indxkl(*),index(*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension xpp(npij,3,*),xqq(npkl,3,*)
      dimension rppq(*)
      dimension txab(npij,3,*),txcd(npkl,3,*)
cnowy dimension apb(npij,*),rapb(npij,*),cpd(npkl,*),rcpd(npkl,*)
      dimension apb(ijpar1,*),rapb(ijpar1,*)
      dimension cpd(klpar1,*),rcpd(klpar1,*)
cnowy
      dimension abcd(nbls1)
c---------------------------------------------------------------
      nstable=nstable + nbls1
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
C
      do 100 i=1,nbls1
         ijkl=index(i)
         ijpar=indxij(ijkl)
         klpar=indxkl(ijkl)
         rppq1=rppq(i)
c
         xwl=( xpp(ijpar,1,lcij) + xqq(klpar,1,lckl) )*rppq1
         ywl=( xpp(ijpar,2,lcij) + xqq(klpar,2,lckl) )*rppq1
         zwl=( xpp(ijpar,3,lcij) + xqq(klpar,3,lckl) )*rppq1
c
         xwp(i,1)=xwl-xp(ijpar,1,lcij)
         xwp(i,2)=ywl-xp(ijpar,2,lcij)
         xwp(i,3)=zwl-xp(ijpar,3,lcij)
c
         xwq(i,1)=xwl-xq(klpar,1,lckl)
         xwq(i,2)=ywl-xq(klpar,2,lckl)
         xwq(i,3)=zwl-xq(klpar,3,lckl)
c
c        for tracij_1, & trackl_1
c
         rcpd1=rcpd(klstart+klpar,lckl)
         abcd(i)=apb(ijstart+ijpar,lcij)*rcpd1
c
         p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
         p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
         p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
  100 continue
c
c the abcd(i) parameter used in Tracy's recursive is defined diferently
c for two shifting directions :
c
c for shifting from 1 to 3    abcd=(a+b)/(c+d)     tracij_1 & _2
c for shifting from 3 to 1    abcd=(c+d)/(a+b)     trackl_1 & _2
c
c The shifting direction is determined NOW by maximum efficiency i.e.
c depends on the angular mom. on left- & right side :
c     if( nsij.ge.nskl) then shift from 1 to 3
c     if( nsij.lt.nskl) then shift from 3 to 1
c However, the shifting direction might be choisen upon different
c criterion. For instance, numerical stability in Tracy's recursive.
c
c Thus, for now the shifting direction is determined by efficiency
c only so we may have both directions of shifting and we need both
c possible abcd :
c The case below is very rare so it is not expensive to do abcd again
c
      if( nsij.lt.nskl ) then
         do i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)
            klpar=indxkl(ijkl)
            abcd(i)=rapb(ijstart+ijpar,lcij)*cpd(klstart+klpar,lckl)
         enddo
      endif
c
c
      end
c====================================================================
      subroutine specase_1(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                   npij,npkl,ii,jj,kk,ll,
     *                   cis,cjs,cks,cls,
     *                   buf,buf1, const,rysx,xpqr,txxr,concoe)
c****

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
c******************************************************
c
      dimension bl(*)
      dimension index(*),indxij(*),indxkl(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension cis(npij,*),cjs(npij,*),cks(npkl,*),cls(npkl,*)
      dimension const(*),rysx(*),concoe(*)
      dimension buf(nbls,*),buf1(nbls1,4)
c
c***************************************************************
c**  this subroutine constitues the special code for
c**  two types of integrals over nbls quartets of primitive shells
c**  1. (ss|ss)
c**  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p or l
c**  these integrals are also contracted here.
c**
c**  this routine is called from the twoe subroutine.
c**
c**
c**  input
c**  ------
c**  all precalculated values for whole block :
c**
c**  const(nbls) - contains consts=pi3*sabcd/(pq*sqrt(ppq)) for all int.
c**  rysx(nbls) - contains  xrys i.e. arguments for fm,ft0 routines
c**  xp,xp      - geometry for p,q
c**
c**  output
c**  ------
c**  buf(nbls,*) -contains final contracted integrals
c***************************************************************
c
c* memory for f00,f11 :
c
      call getmem(nbls1,if00)
      call getmem(nbls1,if11)
c
      if00=if00-1
      if11=if11-1
c
c*******
c
c*
      do 100 i=1,nbls1
      xrys=rysx(i)
      call ft0
      bl(if00+i)=f0
      bl(if11+i)=f1
  100 continue
c
c *** special code for (ss ss) integrals
c
      if(mmax.eq.1) then
          do 2031 i=1,nbls1
          buf1(i,1)=const(i)*bl(if00+i)
 2031     continue
      go to 203
      endif
c
c *** special code for (ps ss), (sp ss) (ss ps) and (ss sp)
c
cxxxx if(mmax.eq.2) then
      intct=0
      if (ityp.gt.1) then
        do 101 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 101    continue
      else if (jtyp.gt.1) then
        do 102 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(2,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 102    continue
      else if (ktyp.gt.1) then
        do 103 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 103    continue
      else
        do 104 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 104    continue
      endif
c*************
        if( max(ityp,jtyp,ktyp,ltyp).eq.3) then
          if(ityp.eq.3) then
            do 1051 i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)
            concoe(i)=const(i)*cis(ijpar,ii)
 1051       continue
          endif
          if(jtyp.eq.3) then
            do 1052 i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)
            concoe(i)=const(i)*cjs(ijpar,jj)
 1052       continue
          endif
          if(ktyp.eq.3) then
            do 1053 i=1,nbls1
            ijkl=index(i)
            klpar=indxkl(ijkl)
            concoe(i)=const(i)*cks(klpar,kk)
 1053       continue
          endif
          if(ltyp.eq.3) then
            do 1054 i=1,nbls1
            ijkl=index(i)
            klpar=indxkl(ijkl)
            concoe(i)=const(i)*cls(klpar,ll)
 1054       continue
          endif
c*
           do 105 i=1,nbls1
           buf1(i,1)=concoe(i)*bl(if00+i)
  105      continue
           intct=1
        endif
c*************
        do 106 i=1,nbls1
          buf1(i,intct+1)=-xpqr(1,i)*const(i)
          buf1(i,intct+2)=-xpqr(2,i)*const(i)
          buf1(i,intct+3)=-xpqr(3,i)*const(i)
106     continue
c
c**********************************************************
c
  203 continue
      if(first) then
         do 204 icx=1,lnijkl
         do 204 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf1(i,icx)
  204    continue
         first=.false.
      else
         do 205 icx=1,lnijkl
         do 205 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf(ijkl,icx)+buf1(i,icx)
  205    continue
      endif
c
c
c release memory
c
      call retmem(2)
c
      return
      end
c====================================================================
c====================================================================
      subroutine precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c-------------------------------------------------------------
c* precalculations for the pairs IJ :
c
      iapb =iabprim
      i1apb=iabprim   + ijdim
      isab =iabprim   + ijdim*2
c--------------
      end
c============================================================
c used with iroute=2 :
c============================================================
      subroutine ab_prim_2(ibl,nparij,ijbl,nbl2,
     *                     datnuc,datbas,inx,iis,jjs,
     *                     apb,rapb,sab,ijcont)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension apb(ijcont),rapb(ijcont),sab(nparij,ijcont)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension sqrtx(400),eexx(400)
c---------------------------------------------------------------
c for the first pair :
c
      ijcs1=ijbl(ibl,  1  )
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
c
      ia=inx(1,ics1)+1
      ie=inx(5,ics1)
      ja=inx(1,jcs1)+1
      je=inx(5,jcs1)
c
      ij=0
      do is=ia,ie
         aa=datbas(1,is)
         do js=ja,je
            bb=datbas(1,js)
            ij=ij+1
            axb=aa*bb
            apb(ij)=aa+bb
            rapb(ij)=one/apb(ij)
            e=axb*rapb(ij)
            eexx(ij)=e
            sqrtx(ij)=rapb(ij)*sqrt(sqrt(axb))**3
         enddo
      enddo
c---------------------------------------------------------------
      do 100 ijpar=1,nparij
      ijcs=ijbl(ibl,ijpar)
      ics=iis(ijcs)
      jcs=jjs(ijcs)
c
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      iat=inx(2,ics)
        ja=inx(1,jcs)+1
        je=inx(5,jcs)
        jat=inx(2,jcs)
        xab=datnuc(2,iat)-datnuc(2,jat)
        yab=datnuc(3,iat)-datnuc(3,jat)
        zab=datnuc(4,iat)-datnuc(4,jat)
        rr=xab*xab+yab*yab+zab*zab
           ij=0
           do 200 is=ia,ie
           do 200 js=ja,je
             ij=ij+1
c
             e=eexx(ij)
             err=e*rr
             if(err.gt.32.d0) then
c......         exp_err=0.d0         ! exp(-32)=1.27*10**-14
c......         sab(ijpar,ij)=rapb(ij)*sqrt3*exp_err
                sab(ijpar,ij)=0.d0
             else
                exp_err=exp(-err)
c.............  sab(ijpar,ij)=rapb(ij)*sqrt3*exp_err
                sqrt3=sqrtx(ij)    ! already multiplied by rapb(ij)
                sab(ijpar,ij)=         sqrt3*exp_err
             endif
c98
c----------------------------------
  200      continue
  100 continue
      end
c============================================================
      subroutine precalc2_2(isupb,bl,mmax,mmax1,       nfumax,nbl2,nbls,
     *                       inx,iis,jjs,ijbl,nblok1,
     *                       ibl,nijbeg,nijend,npij,
     *                       kbl,nklbeg,nklend,npklx,npkl)
      implicit real*8 (a-h,o-z)
      character*3 reuse
c-----------------------------------------------------------
      common /neglect/ eps,eps1,epsr
c
      common /begin/ ijbegin,klbegin
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /ilepar/ lpartot,lpareal
c
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nblok1(2,*)
c---------------------------------------------------------------
      par268=eps1*8.d0
c---------------------------------------------------------------
      ijcurr_f=ijbl(ibl,nijbeg)
      ijcurr_l=ijbl(ibl,nijend)
      call reuseij(ibl,mmax,ijcurr_f,ijcurr_l,reuse)
c-----------------------------------------------------------
      ijbegin=nijbeg
      klbegin=nklbeg
c-----------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_2(iabprim,lcij,iapb,i1apb,isab)
c--------
      lpartot=lpartot+1
      if(reuse.eq.'no ') then
          lpareal=lpareal+1
          call precal2a_2(bl(ibas),bl(inuc),iis,jjs,inx, npij,
     *                  ibl,ijbl,nbl2,nijbeg,nijend,
     *                  bl(iabprim+lcij),bl(iabprim+2*lcij),ijpar1,
     *                  lcij,
     *                  bl(iaa),bl(ibb),
     *                  bl(ieab), bl(icis),bl(icjs),
     *                  bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *                  bl(igcij),ngci1,ngcj1,'left ',par268)
c
c----------
c for abnia
          if(mmax.gt.2) then
             call precal2b_2(mmax1,lcij, bl(i1apb),bl(iabnia))
          endif
      endif
c----------------------------------------
c precalculations for pairs kl
c
      if(kbl.eq.ibl) then
          kabprim=iabprim
          kldim  =ijdim
          klpar1 =ijpar1
          lckl=lcij
      endif
c
      if(npkl.ne.0) then
         call precal2x_2(kabprim,lckl,icpd,i1cpd,iscd)
c
         lpartot=lpartot+1
         lpareal=lpareal+1
c--
         call precal2a_2(bl(ibas),bl(inuc),iis,jjs,inx, npkl,
     *                  kbl,ijbl,nbl2,nklbeg,nklend,
     *                  bl(kabprim+lckl),bl(kabprim+2*lckl),klpar1,
     *                  lckl,
c not an error, icc should be twice :
c2000*                  bl(icc),bl(icc),
     *                  bl(icc),bl(idd),
     *                  bl(iecd), bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igckl),ngck1,ngcl1,'right',par268)
c
c--
      endif
c----------------------------------------
c for cdnia
         if(npkl.ne.0) then
            if(mmax.gt.2) then
               call precal2b_2(mmax1,lckl, bl(i1cpd),bl(icdnia))
            endif
         else
            call precdiag
         endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c_2(lckl,bl(i1cpd),bl(ihabcd),nfumax )
      endif
c---------------------------------------------------------------
      end
c====================================================================
c  Description  of the Precal2b subroutine :
c
c  This subroutine calculates 1 quantity for contracted
c  shell pairs which belong to a given block of contracted
c  of contracted quartets.
c    This subroutine is called (from Blockint) only for blocks
c  with total angular momentum MMAX > 2 ,where MMAX is the sum
c  of angular momentum of four contracted shells (+1) /for exemple :
c  (ss ss):mmax=1, (pp,ps):mmax=4, (dd,ll):mmax=7, (ff,ff):mmax=13 /
c
c  INPUT
c  ------
c  MMAX1 = MMAX-1 , where MMAX is the sum of angular momentum
c  NPIJ  - number of pairs <IJ in a given block of contracted quartets
c  RAPB(ij,lcij) - 1/(a+b) : reversed sum of exponents of each
c                     pair of primitive shells from IJth pair of
c                     contracted shells.
c  OUTPUT
c  -------
c  For each pair of primitive shells ij (kl) from IJ (KL) contracted
c  pair belonging to a given block of quartets with total ang.mom. MMAX
c
c  1. ABNIA(IJ,L,ij) -   L*( 0.5/(a+b) )  with L=1,2,...MMAX-1
c
c-------------------------------------------------------------------
c
      subroutine precal2b_2(mmax1,lcij, rapb, abnia)
      implicit real*8 (a-h,o-z)
cccc  dimension rapb(ijpar1,*), abnia(npij,mmax1,*)
      dimension rapb(*), abnia(mmax1,*)
      data half /0.5d0/
c--------------------------
c     ijstart=ijbegin-1
c--------------------------
      do 150 ij=1,lcij
            apb2=half*rapb(ij)
            abnia(1,ij)=apb2
               do 250 i=2,mmax1
               abnia(i,ij)=abnia(i-1,ij)+apb2
  250          continue
  150 continue
c
      end
c====================================================================
      subroutine precal2c_2(lckl,rcpd,habcd,nfumax)
      implicit real*8 (a-h,o-z)
      common /logicd/ hnia(3,1)
      dimension rcpd(*), habcd(3,nfumax,*)
c-----------------------
c Do it only for one pair (the first one)
c
      do kl=1,lckl
         rcpdkl=rcpd(kl)
         do ifu=1,nfumax
            hx=hnia(1,ifu)
            hy=hnia(2,ifu)
            hz=hnia(3,ifu)
            habcd(1,ifu,kl)=hx*rcpdkl
            habcd(2,ifu,kl)=hy*rcpdkl
            habcd(3,ifu,kl)=hz*rcpdkl
         enddo
      enddo
c
      end
c====================================================================
      subroutine prec4neg_2(nbls,npij,npkl,ndiag,ij,kl,
     1     lc12,lc34,indxij,indxkl,indxr,
     *     estab,estcd,densmax, esti2kl,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, nsym,isymm,
coutput
     *     rho,rppq,rhoapb,rhocpd,rys,const, nbls1,index)
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      logical jump(*)
      common /neglect/ eps,eps1,epsr
c
      dimension indxij(*),indxkl(*),indxr(*),index(*)
ccccc dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension apb(lc12), cpd(lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension                    esti2kl(lc34)
      dimension  densmax(nbls)
      dimension  symfac(*)
c-------
      dimension xp(npij,3,*),xq(npkl,3,*)
ccc   dimension rho(*),rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c now scalars:  rho,rppq,rhoapb,rhocpd  - rho not used at all
      dimension rys(*),const(*)
c symmetry
      dimension isymm(*)
c
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
c---------------------------------------------------------------
      permut=half
c---------------------------------------------------------------
      apb1=apb(ij)
      cpd1=cpd(kl)
c
      abpcd1=apb1+cpd1
      abxcd=apb1*cpd1
      abpcd1r=abpcd1
      abpcdr=one/abpcd1r
      rho1=abxcd*abpcdr
c
      abpcdrcpd1=abpcdr*cpd1
      abpcdrapb1=abpcdr*apb1
      sqrpold=sqrt(abpcdr)
c
      rppq  =abpcdr
      rhoapb=abpcdrcpd1
      rhocpd=abpcdrapb1
c-------------------------------------------
      if(ndiag.eq.0) sqrpold=sqrpold*eps1*8.d0
c-------------------------------------------
      ijkl=0
      ijkl1=0
c-------------------------------------------
c2002 ddmax1=esti2kl(kl)*abpcdr
      ddmax1=esti2kl(kl)
c
      IF(NSYM.EQ.0) THEN
c
c .......no symmetry...................
c
         do 100 ijpar=1,npij
            npklx=npkl
            if(ndiag.eq.0) npklx=ijpar
            esti1=estab(ijpar,ij)*abpcdr
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 100
            endif
            coef1=coefij(ijpar,ij)*sqrpold
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do  50 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 50
               estim=esti1*estcd(klpar,kl)*densmax(ijkl)
               if(estim.gt.epsr) then
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
c
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
c
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
c
                  const(ijkl1)=coef1*coefkl(klpar,kl)
                  if(ndiag.eq.0.and.ijpar.eq.klpar) then
                     const(ijkl1)=const(ijkl1)*permut
                  endif
               endif
  50        continue
 100     continue
c
      ELSE
c
c ..........symmetry...................
c
         do 200 ijpar=1,npij
            npklx=npkl
            if(ndiag.eq.0) npklx=ijpar
               if(jump(ijpar)) then
                  ijkl=ijkl+npklx
                  go to 200
               endif
            esti1=estab(ijpar,ij)*abpcdr
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 200
            endif
            coef1=coefij(ijpar,ij)*sqrpold
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 150 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 150
               estim=esti1*estcd(klpar,kl)*densmax(ijkl)
               if(estim.gt.epsr) then
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
c
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
c
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
c
                  const(ijkl1)=coef1*coefkl(klpar,kl)*symfac(ijkl)
                  if(ndiag.eq.0.and.ijpar.eq.klpar) then
                     const(ijkl1)=const(ijkl1)*permut
                  endif
               endif
 150        continue
 200     continue
c
      ENDIF
c
      nbls1=ijkl1
c
      end
c====================================================================
      subroutine precspec_2(
     $    nbls,          npij,          npkl,          ndiag,
     $    ij,            kl,            lc12,          lc34,
     $    indxij,        indxkl,        indxr,         estab,
     $    estcd,         densmax,       esti2kl,       jump,
     $    symfac,        apb,           cpd,           coefij,
     $    coefkl,        xp,            xq,            rapb,
     $    rcpd,          txab,          txcd,          nsym,
     $    isymm,         rho,           rys,           const,
     $    xpqr,          txxr,          nbls1,         index)
c--------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical jump(*)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
c
      dimension indxij(*),indxkl(*),indxr(*),index(*)
c---------------
      dimension apb(lc12),cpd(lc34)
      dimension rapb(lc12),rcpd(lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension                    esti2kl(lc34)
      dimension  densmax(nbls)
      dimension  symfac(*)
c-------
      dimension txab(npij,3,*),txcd(npkl,3,*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
c
      dimension isymm(*)
c
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
c---------------------------------------------------------------
ccc   par268=eps1*8.d0
      permut=half
c---------------------------------------------------------------
      apb1=apb(ij)
      cpd1=cpd(kl)
c
      abpcd1=apb1+cpd1
      abxcd=apb1*cpd1
      abpcd1r=abpcd1
      abpcdr=one/abpcd1r
      rho1=abxcd*abpcdr
      rho2=rho1
      if(ityp.gt.1.or.jtyp.gt.1) rho2=rho2*rapb(ij)
      if(ktyp.gt.1.or.ltyp.gt.1) rho2=rho2*rcpd(kl)
c
      abpcdrcpd1=abpcdr*cpd1
      abpcdrapb1=abpcdr*apb1
      sqrpold=sqrt(abpcdr)
c-------------------------------------------
      if(ndiag.eq.0) sqrpold=sqrpold*eps1*8.d0
c-------------------------------------------
c2002 ddmax1=esti2kl(kl)*abpcdr
      ddmax1=esti2kl(kl)
c
      ijkl=0
      ijkl1=0
      IF(NSYM.EQ.0) THEN
c.........no symmetry...........
         do 100 ijpar=1,npij
            npklx=npkl
            if(ndiag.eq.0) npklx=ijpar
            esti1=estab(ijpar,ij)*abpcdr
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 100
            endif
            coef1=coefij(ijpar,ij)*sqrpold
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 50 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 50
               estim=esti1*estcd(klpar,kl)*densmax(ijkl)
               if(estim.gt.epsr) then
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
c
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
c
                  xpqr(1,ijkl1)=x1*rho2
                  xpqr(2,ijkl1)=x2*rho2
                  xpqr(3,ijkl1)=x3*rho2
c
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
c
                  const(ijkl1)=coef1*coefkl(klpar,kl)
                  if(ndiag.eq.0.and.ijpar.eq.klpar) then
                     const(ijkl1)=const(ijkl1)*permut
                  endif
                endif
 50         continue
 100     continue
      ELSE
c.........symmetry...........
         do 200 ijpar=1,npij
            npklx=npkl
            if(ndiag.eq.0) npklx=ijpar
            if(jump(ijpar)) then
               ijkl=ijkl+npklx
               go to 200
            endif
            esti1=estab(ijpar,ij)*abpcdr
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 200
            endif
            coef1=coefij(ijpar,ij)*sqrpold
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 150 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 150
               estim=esti1*estcd(klpar,kl)*densmax(ijkl)
               if(estim.gt.epsr) then
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
c
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
c
                  xpqr(1,ijkl1)=x1*rho2
                  xpqr(2,ijkl1)=x2*rho2
                  xpqr(3,ijkl1)=x3*rho2
c
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
c
                  const(ijkl1)=coef1*coefkl(klpar,kl)*symfac(ijkl)
                  if(ndiag.eq.0.and.ijpar.eq.klpar) then
                     const(ijkl1)=const(ijkl1)*permut
                  endif
                endif
 150        continue
 200     continue
      ENDIF
c
      nbls1=ijkl1
      if(nbls1.eq.0) return
c
      if (ityp.gt.1 .or. jtyp.gt.1 ) then
        do 210 i=1,nbls1
        ijkl=index(i)
        ijpar=indxij(ijkl)
        txxr(1,i)=txab(ijpar,1,ij)
        txxr(2,i)=txab(ijpar,2,ij)
        txxr(3,i)=txab(ijpar,3,ij)
  210   continue
      endif
c
      if (ktyp.gt.1 .or. ltyp.gt.1 ) then
        do 220 i=1,nbls1
        ijkl=index(i)
        klpar=indxkl(ijkl)
        txxr(1,i)=txcd(klpar,1,kl)
        txxr(2,i)=txcd(klpar,2,kl)
        txxr(3,i)=txcd(klpar,3,kl)
  220   continue
      endif
cflops
cxxxxx iflop(6)=iflop(6)+3*nbls+2*nbls1
c
      end
c====================================================================
      subroutine xwpq_2(nbls1,xwp,xwq,p1234, ijpar1,lc12, klpar1,lc34,
     *                  lcij,lckl,npij,npkl, indxij,indxkl,index,
     *                  rppq,xp,xq,xpp,xqq,
     *                  txab,txcd,abcd,apb,rcpd,cpd,rapb)
c------------------------------------------------------------------
c  For two primitive pairs lcij and lckl which contribute
c  to all NBLS1 contracted quartets IJKL in a given block
c  4 quatities are calculated for each quartet (xwp,xwq,p1234,abcd)
c
c  INPUT :
c--------
c
c  NBLS1 - blocksize
c  lcij, lckl - pairs of primitive shells
c  NPIJ,NPKL  - number of contracted pairs
c  INDEX(ijkl1) - for ijkl1 quartet in Reduced block shows
c                 the ijkl  quartet in Original block
c  INDXIJ(ijkl) - shows the contracted pair IJPAR for ijkl quartet
c  INDXKL(ijkl) - shows the contracted pair KLPAR for ijkl quartet
c  RPPQ(ijkl1)  = 1/(a+b+c+d)         where a,b,c,d, are exponents
c  XP(ijpar,3,lcij) - coordinates of P
c  XQ(klpar,3,lckl) - coordinates of Q
c  XPP(ijpar,3,lcij)- coordinates of ( a*A + b*B)
c  XQQ(klpar,3,lckl)- coordinates of ( c*C + d*D)
c  TXAB(ijpar,3,  )  - coordinates XAB re-scaled by exponents a or b
c  TXCD(klpar,3,  )  - coordinates XAB re-scaled by exponents c or d
c  APB(ijpar,lcij) =a+b
c  CPD(klpar,lckl) =c+d
c  RAPB(ijpar,lcij)=1/(a+b)
c  RCPD(klpar,lckl)=1/(c+d)
c
c  OUTPUT
c--------
c
c  XWP(ijkl1) - coordinates of W-P,  where W=(XPP+XQQ)/(a+b+c+d)
c  XWQ(ijkl1) - coordinates of W-Q
c  P1234(ijkl1)=(txab+txcd)*rcpd
c  ABCD(ijkl1)= apb*rcpd
c
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     common /begin/ ijbegin,klbegin
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /stabili/ nstable,nostable
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
      dimension indxij(*),indxkl(*),index(*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension xpp(npij,3,*),xqq(npkl,3,*)
      dimension txab(npij,3,*),txcd(npkl,3,*)
cnowy dimension apb(npij,*),rapb(npij,*),cpd(npkl,*),rcpd(npkl,*)
c     dimension rppq(*)
c     dimension apb(ijpar1,*),rapb(ijpar1,*)
c     dimension cpd(klpar1,*),rcpd(klpar1,*)
c     dimension abcd(nbls1)
cnowy
      dimension apb(lc12),rapb(lc12)
      dimension cpd(lc34),rcpd(lc34)
c---------------------------------------------------------------
      nstable=nstable + nbls1
c---------------------------------------------------------------
c     ijstart=ijbegin-1
c     klstart=klbegin-1
c---------------------------------------------------------------
      rppq1=rppq
      rcpd1=rcpd(lckl)
c
      if( nsij.ge.nskl ) then
         abcd=apb(lcij)*rcpd1
      else
         abcd=rapb(lcij)*cpd(lckl)
      endif
c
      do 100 i=1,nbls1
         ijkl=index(i)
         ijpar=indxij(ijkl)
         klpar=indxkl(ijkl)
c
         xwl=( xpp(ijpar,1,lcij) + xqq(klpar,1,lckl) )*rppq1
         ywl=( xpp(ijpar,2,lcij) + xqq(klpar,2,lckl) )*rppq1
         zwl=( xpp(ijpar,3,lcij) + xqq(klpar,3,lckl) )*rppq1
c
         xwp(i,1)=xwl-xp(ijpar,1,lcij)
         xwp(i,2)=ywl-xp(ijpar,2,lcij)
         xwp(i,3)=zwl-xp(ijpar,3,lcij)
c
         xwq(i,1)=xwl-xq(klpar,1,lckl)
         xwq(i,2)=ywl-xq(klpar,2,lckl)
         xwq(i,3)=zwl-xq(klpar,3,lckl)
c
c        for tracij_2 & trackl_2
c
         p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
         p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
         p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
c
  100 continue
c
      end
c====================================================================
      subroutine specase_2(bl,first,nbls,nbls1, index,
     *                     npij,npkl,ii,jj,kk,ll,
     *                     cis,cjs,cks,cls,
     *                     buf,buf1, const,rysx,xpqr,txxr,concoe)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      dimension bl(*)
      dimension index(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension cis(*),cjs(*),cks(*),cls(*)
      dimension const(*),rysx(*),concoe(*)
      dimension buf(nbls,*),buf1(nbls1,4)
c---------------------------------------------------------------
c**  this subroutine constitues the special code for
c**  two types of integrals over nbls quartets of primitive shells
c**  1. (ss|ss)
c**  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p or l
c**  these integrals are also contracted here.
c**
c**  this routine is called from the twoe subroutine.
c**
c**
c**  input
c**  ------
c**  all precalculated values for whole block :
c**
c**  const(nbls) - contains consts=pi3*sabcd/(pq*sqrt(ppq)) for all int.
c**  rysx(nbls) - contains  xrys i.e. arguments for fm,ft0 routines
c**  xp,xp      - geometry for p,q
c**
c**  output
c**  ------
c**  buf(nbls,*) -contains final contracted integrals
c---------------------------------------------------------------
c
c* memory for f00,f11 :
c
      call getmem(nbls1,if00)
      call getmem(nbls1,if11)
c
      if00=if00-1
      if11=if11-1
c
c*******
c
c*
      do 100 i=1,nbls1
      xrys=rysx(i)
      call ft0
      bl(if00+i)=f0
      bl(if11+i)=f1
  100 continue
c
c *** special code for (ss ss) integrals
c
      if(mmax.eq.1) then
          do 2031 i=1,nbls1
          buf1(i,1)=const(i)*bl(if00+i)
 2031     continue
      go to 203
      endif
c
c *** special code for (ps ss), (sp ss) (ss ps) and (ss sp)
c
cxxxx if(mmax.eq.2) then
      intct=0
      if (ityp.gt.1) then
        do 101 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 101    continue
      else if (jtyp.gt.1) then
        do 102 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(2,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 102    continue
      else if (ktyp.gt.1) then
        do 103 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 103    continue
      else
        do 104 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 104    continue
      endif
c------------------- NEW ------------------------
        if( max(ityp,jtyp,ktyp,ltyp).eq.3) then
c
          if(ityp.eq.3) then
            coeff=cis(ii)
            do 1051 i=1,nbls1
            concoe(i)=const(i)*coeff
 1051       continue
          endif
          if(jtyp.eq.3) then
            coeff=cjs(jj)
            do 1052 i=1,nbls1
            concoe(i)=const(i)*coeff
 1052       continue
          endif
          if(ktyp.eq.3) then
            coeff=cks(kk)
            do 1053 i=1,nbls1
            concoe(i)=const(i)*coeff
 1053       continue
          endif
          if(ltyp.eq.3) then
            coeff=cls(ll)
            do 1054 i=1,nbls1
            concoe(i)=const(i)*coeff
 1054       continue
          endif
           do 105 i=1,nbls1
           buf1(i,1)=concoe(i)*bl(if00+i)
  105      continue
           intct=1
        endif
c------------------- NEW end --------------------
        do 106 i=1,nbls1
          buf1(i,intct+1)=-xpqr(1,i)*const(i)
          buf1(i,intct+2)=-xpqr(2,i)*const(i)
          buf1(i,intct+3)=-xpqr(3,i)*const(i)
106     continue
c------------------------------------------------
  203 continue
      if(first) then
         do 204 icx=1,lnijkl
         do 204 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf1(i,icx)
  204    continue
         first=.false.
      else
         do 205 icx=1,lnijkl
         do 205 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf(ijkl,icx)+buf1(i,icx)
  205    continue
      endif
c
c release memory
c
      call retmem(2)
c
      end
c====================================================================
      subroutine precal2x_2(iabprim,lcij, iapb,i1apb,isab)
c-------------------------------------------------------------
c* precalculations for the pairs IJ :
c
c     iapb =iabprim
c     i1apb=iabprim   + ijdim
c     isab =iabprim   + ijdim*2
      iapb =iabprim
      i1apb=iabprim+lcij
      isab =iabprim+2*lcij
c--------------
      end
c====================================================================
      subroutine precal2a_2(datbas,datnuc,iis,jjs,inx,npij,
     *                    ibl,ijbl,nbl2,nijbeg,nijend,
     *                    rapb,sab,ijpar1, lcij,
     *                    aaa,bbb,
     *                    estab, cis,cjs,
     *                    xab,xparij,coefij,factij,txab,
     *                    gcij,ngci1,ngcj1,which,par268 )
c
      implicit real*8 (a-h,o-z)
      character*5 which
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      COMMON /types/itype,jtype,ktype,ltype,ityp,jtyp,ktyp,ltyp
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension rapb(lcij),sab(ijpar1,lcij)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension cis(*),cjs(*)
      dimension xab(npij,3), xparij(npij,3,lcij,3)
      dimension estab (npij,lcij)
      dimension coefij(npij,lcij), factij(lcij)
      dimension xa(3),xb(3),txab(npij,3,*)
      dimension gcij(ngci1,ngcj1,lcij)
      dimension aaa(     *),bbb(     *)   ! dimen= lci,lcj -contr.
c local arrays (dimensioned as contraction max length) :
      dimension est_ijx(30,30),coef_ijx(30,30)
c---------------------------------------------------------------
      if(which.eq.'left ') then
          itypp=ityp
          jtypp=jtyp
          ngcii=ngci
          ngcjj=ngcj
          nqii=nqi
          nqjj=nqj
          lcii=lci
          lcjj=lcj
      else
          itypp=ktyp
          jtypp=ltyp
          ngcii=ngck
          ngcjj=ngcl
          nqii=nqk
          nqjj=nql
          lcii=lck
          lcjj=lcl
      endif
c
c-------------------------------------------------------------
c for primitive integral neglect : to prevent from throwing away
c small primitive contributions which might sum up to a big value
c
c98   contrij=dble(lcii*lcjj)
c98   contrij=sqrt(contrij)
c-------------------------------------------------------------
c for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c-------------------------------------------------------------
c precalculate all possible quantities for the first pair only
c
      ijp1=nijbeg
      call first_pair(nbl2,ibl,ijbl,inx,ijp1,itypp,jtypp,
     *                datbas,par268,which,ngctot,
     *                aaa ,cis,bbb ,cjs,
     *                est_ijx,coef_ijx,factij,
     *                gcij,ngcii,ngcjj)
c-------------------------------------------------------------
c precalculations for the pairs IJ :
c
c2002 instablity issue :
      boamax=0.d0
      rapbmax=0.d0
c
      ijpar=0
      do 100 ijp=nijbeg,nijend
      ijpar=ijpar+1
        ijcs=ijbl(ibl,ijp)
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        fact1=one
        if(ics.eq.jcs) fact1=fact1*half
c      starting contr
        ia=inx(1,ics)+1
        ja=inx(1,jcs)+1
c      last contr
        ie=inx(5,ics)
        je=inx(5,jcs)
c.......................................................
        iatom=inx(2,ics)
        jatom=inx(2,jcs)
c
c center23 (atom number 0 is the one with s-orbital/zero-exp.)
c
        if(iatom.eq.0) iatom=jatom
        if(jatom.eq.0) jatom=iatom
c
          do 150 i=1,3
          if(iatom.eq.0 .and. jatom.eq.0) then
            xa(i)=zero
            xb(i)=zero
          else
            xa(i)=datnuc(1+i,iatom)
            xb(i)=datnuc(1+i,jatom)
          endif
          xab(ijpar,i)=xa(i)-xb(i)
  150     continue
c.......................................................
c
         ji=0
         is1=0
         do 200 is=ia,ie
            is1=is1+1
            aa=datbas(1,is)
c2002
            aor1=1.d0/max(1.d0,aa)
c
         js1=0
         do 200 js=ja,je
            js1=js1+1
            ji=ji+1
            bb=datbas(1,js)
c2002
            boamax=max(boamax,bb*aor1)      ! b/max(1,a)
c---------------------------------------------
c              if(ngctot.eq.0) then
c                 coefij(ijpar,ji)=coef_ijx(is1,js1)*fact1
c                 estab(ijpar,ji)=est_ijx(is1,js1)
c              else
c                 gen.contr.
c                 coefij(ijpar,ji )=fact1
c                 estab(ijpar,ji)= one      ! should be =fact1
c                 if(which.eq.'right') coefij(ijpar,ji)=fact1*par268
c              endif
c----------------------------------------------------
               coefij(ijpar,ji)=coef_ijx(is1,js1)*fact1
               estab(ijpar,ji)=est_ijx(is1,js1)
c----------------------------------------------------
            rapb1=rapb(ji)
            sab1 =sab(ijp,ji)
            aa1=aa*rapb1
            bb1=bb*rapb1
c
c-overlap   coefij(ijpar,ji)=coefij(ijpar,ji)*rapb1*sab1
            coefij(ijpar,ji)=coefij(ijpar,ji)*sab1
            estab(ijpar,ji)=estab(ijpar,ji)*sab1
c
            xpn_max=1.d0
            do 230 l=1,3
               xparij(ijpar,l,ji,1)=aa1*xa(l)+bb1*xb(l) ! xp(ijpar,l,ji
               xxl=xa(l)
               if(nqii.lt.nqjj) xxl=xb(l)
               xparij(ijpar,l,ji,2)=xparij(ijpar,l,ji,1)-xxl  ! xpn
               xparij(ijpar,l,ji,3)=aa*xa(l)+bb*xb(l)         ! xpp
c center23:
               if(aa.le.zero .or. bb.le.zero) then
                  xparij(ijpar,l,ji,2)=zero
               endif
ckw99
c for neglect : include factor (P-A) which appears in the O-S recursive
c               and multiplies (ss/ss)(0) integral
c
               xpn_abs=abs( xparij(ijpar,l,ji,2) )
               if(xpn_abs.gt.xpn_max) xpn_max=xpn_abs
c
  230       continue
c------------------------------------------------------------------
            estab(ijpar,ji)=estab(ijpar,ji)*xpn_max
c
c square it for neglect :
            estab(ijpar,ji)=estab(ijpar,ji)*estab(ijpar,ji)
c2002
            rapbmax=max(rapbmax,rapb1)
c------------------------------------------------------------------
  200    continue        ! end of the loop over contractions
c------------------------------------------------------------------
  100 continue        ! end of the loop over contrcted pairs
c------------------------------------------------------------------
c2002
      boamax=2.d0*boamax
      if(which.eq.'left ') then
         call setrval('boamax',boamax)
         call setrval('rapbmax',rapbmax)
         call setrval('docmax',boamax)
         call setrval('rcpdmax',rapbmax)
      else
         call setrval('docmax',boamax)
         call setrval('rcpdmax',rapbmax)
      endif
c2002
ctxab
c
      if(nqii.ge.nqjj) then
         ijs1=0
         do 151 is1=1,lcii
         aa1=aaa(is1)
         do 151 js1=1,lcjj
         bb1=bbb(js1)
         ijs1=ijs1+1
         do 151 ijpar=1,npij
            txab(ijpar,1,ijs1)=-bb1*xab(ijpar,1)
            txab(ijpar,2,ijs1)=-bb1*xab(ijpar,2)
            txab(ijpar,3,ijs1)=-bb1*xab(ijpar,3)
c center23:
            if(aa1.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  151    continue
      else
         ijs1=0
         do 152 is1=1,lcii
         aa1=aaa(is1)
         do 152 js1=1,lcjj
         bb1=bbb(js1)
         ijs1=ijs1+1
         do 152 ijpar=1,npij
            txab(ijpar,1,ijs1)= aa1*xab(ijpar,1)
            txab(ijpar,2,ijs1)= aa1*xab(ijpar,2)
            txab(ijpar,3,ijs1)= aa1*xab(ijpar,3)
c center23:
            if(bb1.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  152    continue
      endif
c
c for the case with mmax=2 - special cases
c
      if ( mmax.eq.2 ) then
         if (itypp.gt.1 .or. jtypp.gt.1 ) then
           do 153 ijs1=1,lcij
           rapb1=rapb(ijs1)
           do 153 ijpar=1,npij
            txab(ijpar,1,ijs1)=txab(ijpar,1,ijs1)*rapb1
            txab(ijpar,2,ijs1)=txab(ijpar,2,ijs1)*rapb1
            txab(ijpar,3,ijs1)=txab(ijpar,3,ijs1)*rapb1
  153      continue
         endif
      endif
c-------------------------------------------------
c2000
      if(where.eq.'forc'.or. where.eq.'hess') then
         do is1=1,lcii
            aaa(is1)=2.d0*aaa(is1)
         enddo
         if(which.eq.'left ') then
            do js1=1,lcjj
               bbb(js1)=2.d0*bbb(js1)
            enddo
         endif
      endif
c-------------------------------------------------
      end
c====================================================================
      subroutine first_pair(nbl2,ibl,ijbl,inx,ijp1,itypp,jtypp,
     *                      datbas,par268,which,ngctot,
     *                      aexp,cis,bexp,cjs,
     *                      est_ijx,coef_ijx,factij,
     *                      gcij,ngcii,ngcjj)
      implicit real*8 (a-h,o-z)
      character*5 which
      dimension datbas(13,*)
ccc   dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
      dimension inx(12,*),ijbl(nbl2,*)
c
      dimension aexp(*),cis(*)
      dimension bexp(*),cjs(*)
c
      dimension est_ijx(30,30),coef_ijx(30,30)
      dimension factij(*)
      dimension gcij(ngcii+1,ngcjj+1,*)
c-----------------------------------------------------------------------
c output :
c     aexp(is1)=aa
c     bexp(js1)=bb
c     cis(is1)=csi
c     cjs(js1)=csj
c     factij(ji)=csi*csj
c     coef_ijx(is1,js1)=coefi*coefj
c     est_ijx(is1,js1)=max(abs(csi),abs(cpi))*max(abs(csj),abs(cpj))
c-----------------------------------------------------------------------
       ijcs=ijbl(ibl,ijp1)
ccc    ics=iis(ijcs)
ccc    jcs=jjs(ijcs)
       call get_ij_half(ijcs, ics, jcs)
c
c      starting contr
       ia=inx(1,ics)+1
       ja=inx(1,jcs)+1
c      last contr
       ie=inx(5,ics)
       je=inx(5,jcs)
c
       is1=0
       do is=ia,ie
          is1=is1+1
          aa=datbas(1,is)
          aexp(is1)=aa
       enddo
       js1=0
       do js=ja,je
          js1=js1+1
          bb=datbas(1,js)
          bexp(js1)=bb
       enddo
c
       ji=0
       is1=0
       do is=ia,ie
          is1=is1+1
          if(ngctot.eq.0) then
             csi=datbas(2,is)
             cpi=datbas(3,is)
             coefi=csi
             est_i=max( abs(csi),abs(cpi) )
             if(itypp.eq.3) then
                coefi=cpi
                csi=csi/cpi
             endif
             if(which.eq.'right') coefi=coefi*par268
             cis(is1)=csi
          endif
c
          js1=0
          do js=ja,je
             ji=ji+1
             js1=js1+1
             if(ngctot.eq.0) then
                csj=datbas(2,js)
                cpj=datbas(3,js)
                coefj=csj
                est_j=max( abs(csj),abs(cpj) )
                est_ijx(is1,js1)=est_i*est_j
                if(jtypp.eq.3) then
                   coefj=cpj
                   csj=csj/cpj
                   if(itypp.eq.3) factij(ji)=csi*csj
                endif
                coef_ijx(is1,js1)=coefi*coefj
                cjs(js1)=csj
             else
cccc            general contraction
                est_ijx(is1,js1)=1.d0
                coef_ijx(is1,js1)=1.d0
                if(which.eq.'right') coef_ijx(is1,js1)=par268
                do ig=0,ngcii
                   gci=datbas(ig+2,is)
                   do jg=0,ngcjj
                     gcj=gci*datbas(jg+2,js)
                     gcij(ig+1,jg+1,ji)=gcj
                   enddo
                enddo
             endif
          enddo
       enddo
c
c-----------------------------------------------------------------------
      end
c=======================================================================
      subroutine xwpq_11(nbls1,xwp,xwq,p1234, ijpar1,lc12, klpar1,lc34,
     *                  lcij,lckl,npij,npkl, indxij,indxkl,index,
     *                  rppq,xp,xq,xpp,xqq,
     *                  txab,txcd,abcd,apb,rcpd,cpd,rapb,
     *                  where,
     *                  immax,kmmax,lobsa)
c-----------------------------------------------------------------
c  For two primitive pairs lcij and lckl which contribute
c  to all NBLS1 contracted quartets IJKL in a given block
c  4 quatities are calculated for each quartet (xwp,xwq,p1234,abcd)
c
c  INPUT :
c--------
c
c  NBLS1 - blocksize
c  lcij, lckl - pairs of primitive shells
c  NPIJ,NPKL  - number of contracted pairs
c  INDEX(ijkl1) - for ijkl1 quartet in Reduced block shows
c                 the ijkl  quartet in Original block
c  INDXIJ(ijkl) - shows the contracted pair IJPAR for ijkl quartet
c  INDXKL(ijkl) - shows the contracted pair KLPAR for ijkl quartet
c  RPPQ(ijkl1)  = 1/(a+b+c+d)         where a,b,c,d, are exponents
c  XP(ijpar,3,lcij) - coordinates of P
c  XQ(klpar,3,lckl) - coordinates of Q
c  XPP(ijpar,3,lcij)- coordinates of ( a*A + b*B)
c  XQQ(klpar,3,lckl)- coordinates of ( c*C + d*D)
c  TXAB(ijpar,3,  )  - coordinates XAB re-scaled by exponents a or b
c  TXCD(klpar,3,  )  - coordinates XAB re-scaled by exponents c or d
c  APB(ijpar,lcij) =a+b
c  CPD(klpar,lckl) =c+d
c  RAPB(ijpar,lcij)=1/(a+b)
c  RCPD(klpar,lckl)=1/(c+d)
c
c  OUTPUT
c--------
c  XWP(ijkl1) - coordinates of W-P,  where W=(XPP+XQQ)/(a+b+c+d)
c  XWQ(ijkl1) - coordinates of W-Q
c  P1234(ijkl1)=(txab+txcd)*(rcpd OR rapb)
c  ABCD(ijkl1)= apb*rcpd   OR  cpd*rapb
c-----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*4 where
      logical stable
      common /begin/ ijbegin,klbegin
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common /neglect/ eps,eps1,epsr
c
      common /stabili/ nstable,nostable
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
c
      dimension indxij(*),indxkl(*),index(*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension xpp(npij,3,*),xqq(npkl,3,*)
      dimension rppq(*)
      dimension txab(npij,3,*),txcd(npkl,3,*)
cnowy dimension apb(npij,*),rapb(npij,*),cpd(npkl,*),rcpd(npkl,*)
      dimension apb(ijpar1,*),rapb(ijpar1,*)
      dimension cpd(klpar1,*),rcpd(klpar1,*)
cnowy
      dimension abcd(nbls1)
c---------------------------------------------------------------
      lost_allow=15+log10(eps)   ! 15 from double precision accuracy
      if(lost_allow.lt.0) lost_allow=0
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
      rapb1_max=0.d0
      rcpd1_max=0.d0
      do 100 i=1,nbls1
         ijkl=index(i)
         ijpar=indxij(ijkl)
         klpar=indxkl(ijkl)
         rppq1=rppq(i)
c
         xwl=( xpp(ijpar,1,lcij) + xqq(klpar,1,lckl) )*rppq1
         ywl=( xpp(ijpar,2,lcij) + xqq(klpar,2,lckl) )*rppq1
         zwl=( xpp(ijpar,3,lcij) + xqq(klpar,3,lckl) )*rppq1
c
         xwp(i,1)=xwl-xp(ijpar,1,lcij)
         xwp(i,2)=ywl-xp(ijpar,2,lcij)
         xwp(i,3)=zwl-xp(ijpar,3,lcij)
c
         xwq(i,1)=xwl-xq(klpar,1,lckl)
         xwq(i,2)=ywl-xq(klpar,2,lckl)
         xwq(i,3)=zwl-xq(klpar,3,lckl)
c
c        for tracij_1, & trackl_1
c
         rcpd1=rcpd(klstart+klpar,lckl)
         abcd(i)=apb(ijstart+ijpar,lcij)*rcpd1
c
         p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
         p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
         p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
c
         rapb1=rapb(ijstart+ijpar,lcij)
         rapb1_max=max(rapb1_max,rapb1)
         rcpd1_max=max(rcpd1_max,rcpd1)
  100 continue
c
c-----------------------------------------------------------
c the abcd(i) parameter used in Tracy's recursive is defined diferently
c for two shifting directions :
c
c for shifting from 1 to 3    abcd=(a+b)/(c+d)     tracij_1 & _2
c for shifting from 3 to 1    abcd=(c+d)/(a+b)     trackl_1 & _2
c
c-----------------------------------------------------------
c numerical stability problem in Tracy's reqursieve:
c-----------------------------------------------------------
      nshifts=min(nsij,nskl)-1
      if(where.eq.'shif' .or. where.eq.'forc') then
         nshifts=nshifts-1
      endif
      if(where.eq.'hess') then
         nshifts=nshifts-2
      endif
c-----------------------------------------------------------
      if(nsij.ge.nskl) then
c        shifting from 1 to 3 : check if stable
c
         call getrval('boa_max',boa_max) ! setup when xwpq_11 is called
         call digit_lost(boa_max,rcpd1_max,nshifts,lost_digit)
         stable=.true.
         if(lost_digit.gt.lost_allow) stable=.false.
         if(stable) then
            immax=mmax-2
            kmmax=nskl-2
            lobsa=2
            if(lshelij.gt.0) lobsa=1
            nstable=nstable + nbls1
         else
            immax=nsij-2
            kmmax=mmax-2
            lobsa=4
            if(lshelkl.gt.0) lobsa=3
            do i=1,nbls1
               abcd(i)=1.d0/abcd(i)
            enddo
            nostable=nostable + nbls1
         endif
      endif
c-----------------------------------------------------------
      if(nsij.lt.nskl) then
c        shifting from 3 to 1 : check if stable
c
         call getrval('doc_max',doc_max) ! setup when xwpq_11 is called
         call digit_lost(doc_max,rapb1_max,nshifts,lost_digit)
         stable=.true.
         if(lost_digit.gt.lost_allow) stable=.false.
         if(stable) then
            immax=nsij-2
            kmmax=mmax-2
            lobsa=4
            if(lshelkl.gt.0) lobsa=3
            do i=1,nbls1
               abcd(i)=1.d0/abcd(i)
            enddo
            nstable=nstable + nbls1
         else
            immax=mmax-2
            kmmax=nskl-2
            lobsa=2
            if(lshelij.gt.0) lobsa=1
            nostable=nostable + nbls1
         endif
      endif
c-----------------------------------------------------------
      end
c====================================================================
      subroutine xwpq_22(nbls1,xwp,xwq,p1234, ijpar1,lc12, klpar1,lc34,
     *                  lcij,lckl,npij,npkl, indxij,indxkl,index,
     *                  rppq,xp,xq,xpp,xqq,
     *                  txab,txcd,abcd,apb,rcpd,cpd,rapb,
     *                  where,
     *                  immax,kmmax,lobsa)
c------------------------------------------------------------------
c  For two primitive pairs lcij and lckl which contribute
c  to all NBLS1 contracted quartets IJKL in a given block
c  4 quatities are calculated for each quartet (xwp,xwq,p1234,abcd)
c
c  INPUT :
c--------
c
c  NBLS1 - blocksize
c  lcij, lckl - pairs of primitive shells
c  NPIJ,NPKL  - number of contracted pairs
c  INDEX(ijkl1) - for ijkl1 quartet in Reduced block shows
c                 the ijkl  quartet in Original block
c  INDXIJ(ijkl) - shows the contracted pair IJPAR for ijkl quartet
c  INDXKL(ijkl) - shows the contracted pair KLPAR for ijkl quartet
c  RPPQ(ijkl1)  = 1/(a+b+c+d)         where a,b,c,d, are exponents
c  XP(ijpar,3,lcij) - coordinates of P
c  XQ(klpar,3,lckl) - coordinates of Q
c  XPP(ijpar,3,lcij)- coordinates of ( a*A + b*B)
c  XQQ(klpar,3,lckl)- coordinates of ( c*C + d*D)
c  TXAB(ijpar,3,  )  - coordinates XAB re-scaled by exponents a or b
c  TXCD(klpar,3,  )  - coordinates XAB re-scaled by exponents c or d
c  APB(ijpar,lcij) =a+b
c  CPD(klpar,lckl) =c+d
c  RAPB(ijpar,lcij)=1/(a+b)
c  RCPD(klpar,lckl)=1/(c+d)
c
c  OUTPUT
c--------
c
c  XWP(ijkl1) - coordinates of W-P,  where W=(XPP+XQQ)/(a+b+c+d)
c  XWQ(ijkl1) - coordinates of W-Q
c  P1234(ijkl1)=(txab+txcd)*rcpd
c  ABCD(ijkl1)= apb*rcpd
c
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*4 where
      logical stable
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common /neglect/ eps,eps1,epsr
c
      common /stabili/ nstable,nostable
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
      dimension indxij(*),indxkl(*),index(*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension xpp(npij,3,*),xqq(npkl,3,*)
      dimension txab(npij,3,*),txcd(npkl,3,*)
cnowy
      dimension apb(lc12),rapb(lc12)
      dimension cpd(lc34),rcpd(lc34)
c---------------------------------------------------------------
      lost_allow=15+log10(eps)   ! 15 from double precision accuracy
      if(lost_allow.lt.0) lost_allow=0
c---------------------------------------------------------------
      rppq1=rppq
      rcpd1=rcpd(lckl)
c
      abcd=apb(lcij)*rcpd1
c
ckw2000 now always calculate abcd as (a+b)/c+d)
c
      do 100 i=1,nbls1
         ijkl=index(i)
         ijpar=indxij(ijkl)
         klpar=indxkl(ijkl)
c
         xwl=( xpp(ijpar,1,lcij) + xqq(klpar,1,lckl) )*rppq1
         ywl=( xpp(ijpar,2,lcij) + xqq(klpar,2,lckl) )*rppq1
         zwl=( xpp(ijpar,3,lcij) + xqq(klpar,3,lckl) )*rppq1
c
         xwp(i,1)=xwl-xp(ijpar,1,lcij)
         xwp(i,2)=ywl-xp(ijpar,2,lcij)
         xwp(i,3)=zwl-xp(ijpar,3,lcij)
c
         xwq(i,1)=xwl-xq(klpar,1,lckl)
         xwq(i,2)=ywl-xq(klpar,2,lckl)
         xwq(i,3)=zwl-xq(klpar,3,lckl)
c
c        for tracij_2 & trackl_2
c
         p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
         p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
         p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
c
  100 continue
c
c-----------------------------------------------------------
c numerical stability problem in Tracy's reqursieve:
c-----------------------------------------------------------
      nshifts=min(nsij,nskl)-1
      if(where.eq.'shif' .or. where.eq.'forc') then
         nshifts=nshifts-1
      endif
      if(where.eq.'hess') then
         nshifts=nshifts-2
      endif
c-----------------------------------------------------------
      if(nsij.ge.nskl) then
c        shifting from 1 to 3 : check if stable
c
         call getrval('boa_max',boa_max) ! setup when xwpq_11 is called
         call digit_lost(boa_max,rcpd1,nshifts,lost_digit)
         stable=.true.
         if(lost_digit.gt.lost_allow) stable=.false.
         if(stable) then
            immax=mmax-2
            kmmax=nskl-2
            lobsa=2
            if(lshelij.gt.0) lobsa=1
            nstable=nstable + nbls1
         else
            immax=nsij-2
            kmmax=mmax-2
            lobsa=4
            if(lshelkl.gt.0) lobsa=3
            abcd=1.d0/abcd
            nostable=nostable + nbls1
         endif
      endif
c-----------------------------------------------------------
      if(nsij.lt.nskl) then
c        shifting from 3 to 1 : check if stable
c
         call getrval('doc_max',doc_max) ! setup when xwpq_11 is called
         rapb1=rapb(lcij)
         call digit_lost(doc_max,rapb1,nshifts,lost_digit)
         stable=.true.
         if(lost_digit.gt.lost_allow) stable=.false.
         if(stable) then
            immax=nsij-2
            kmmax=mmax-2
            lobsa=4
            if(lshelkl.gt.0) lobsa=3
            abcd=1.d0/abcd
            nstable=nstable + nbls1
         else
            immax=mmax-2
            kmmax=nskl-2
            lobsa=2
            if(lshelij.gt.0) lobsa=1
            nostable=nostable + nbls1
         endif
      endif
c-----------------------------------------------------------
      end
c====================================================================
      subroutine check_stab(stable,where)
      implicit real*8 (a-h,o-z)
      character*4 where
      logical stable
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common /neglect/ eps,eps1,epsr
c-----------------------------------------------------------
      stable=.true.
c-----------------------------------------------------------
      call getival('stab',istab)
c     write(6,*)'instability status : istab=',istab
c
      if(istab.eq.1) return
c-----------------------------------------------------------
      maxangmom=mmax
      nshifts=min(nsij,nskl)-1
      if(where.eq.'shif'.or.where.eq.'forc') then
         maxangmom=mmax-2
         nshifts=nshifts-1
      endif
      if(where.eq.'hess') then
         maxangmom=mmax-4
         nshifts=nshifts-2
      endif
c
C2004 if( maxangmom.le.6 ) return            ! dp|dp = 7
ccccc if( maxangmom.le.5 ) return            ! dp|pp = 6
      if( maxangmom.le.4 ) return            ! pp|pp = 5
      if( max(nqi,nqj,nqk,nql).le.2) return  ! up to p only
      if(nshifts.le.1) return
c-----------------------------------------------------------
      xlogeps=log10(eps)
      lost_allow=15+xlogeps   ! 15 from double precision accuracy
      if(where.eq.'shif'.or.where.eq.'forc') lost_allow=14+xlogeps
      if(where.eq.'hess') lost_allow=13+xlogeps
      if(lost_allow.lt.0) lost_allow=0
c-----------------------------------------------------------
      call getrval('boamax',boamax)
      call getrval('rapbmax',rapbmax)
      call getrval('docmax',docmax)
      call getrval('rcpdmax',rcpdmax)
c-----------------------------------------------------------------------
c numerical stability problem in Tracy's reqursieve IF going from 1 to 3
c
      if(nsij.ge.nskl) then
c        shifting in tracy's req. from 1 to 3
         call digit_lost(boamax,rcpdmax,nshifts,lost_digit)
         if(lost_digit.gt.lost_allow) stable=.false.
      endif
c-----------------------------------------------------------------------
c numerical stability problem in Tracy's reqursieve IF going from 3 to 1
c
      if(nsij.lt.nskl) then
c        shifting in tracy's req. from 3 to 1
         call digit_lost(docmax,rapbmax,nshifts,lost_digit)
         if(lost_digit.gt.lost_allow) stable=.false.
      endif
c-----------------------------------------------------------------------
      end
c====================================================================
      subroutine digit_lost(boa_max,rcpd1_max,nshifts,lost_digit)
      implicit real*8 (a-h,o-z)
c
         lost_digit=0
         unstab_max=boa_max*rcpd1_max
         if(unstab_max.gt.10.d0 ) then
            xost_digit  = log10(unstab_max)   ! lost in one shifting
            xost_digit  = nshifts*xost_digit  ! lost in nskl-1 shifts
            lost_digit  = nint( xost_digit )
         endif
cccc     write(6,*)'2b/a =',boa_max,' 1/(c+d)=',rcpd1_max,' ',unstab_max
c
      end
c====================================================================
