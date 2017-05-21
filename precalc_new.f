c====================================================================
c===========                               ==========================
c=========== NEW ROUTINES for NEW BLOCKING ==========================
c===========                               ==========================
c====================================================================
      subroutine prec2ij_new(ibl, bl)

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
c---------------------------------------
c1999
c get needed parameters from depository :
c
      call getival('inx_1' ,inx_1 )
      call getival('ibas_1',ibas_1)
      call getival('inx_2' ,inx_2 )
      call getival('ibas_2',ibas_2)
c
      call getival('nparx',nparx)
      call getival('ijblx',ijblx)
      call getival('mapijblx',map_ij_blx)
ctemp call getival('blocksij',nbl2)
      call getival('blocksij',nbl2_ij)
      call getival('blpredij',nbl2_pred)
c
      nbl2_dim=nbl2_pred
c---------------------------------------
      if(icheck.gt.0) then
         call getmem(0,l0)
         return
      endif
c---------------------------------------
      call secund(tprec2b)
c---------------------------------------
      call dimenij_new(ibl,bl(inx_1),bl(inx_2),
     *                 bl(nparx),bl(ijblx),nbl2_dim,
     *                 nparij,ijdim,ijcont)    ! output
c---------------------------------------
      ijpar1=nparij
      IF( iroute.eq.1 ) THEN
         call getmem(3*ijdim,iabprim)
         call ab_prim_1_new(ibl,bl(inuc),nparij,bl(ijblx),nbl2_dim,
     *                      bl(ibas_1),bl(inx_1),
     *                      bl(ibas_2),bl(inx_2),
     *                      bl(iabprim),ijcont )
      ELSE
         call getmem(2*ijcont+ijdim,iabprim)
         iapb =iabprim
         i1apb=iabprim+ijcont
         isab =iabprim+ijcont*2
c
         call ab_prim_2_new(ibl,bl(inuc),nparij,bl(ijblx),nbl2_dim,
     *                      bl(ibas_1),bl(inx_1),
     *                      bl(ibas_2),bl(inx_2),
     *                      bl(iapb),bl(i1apb),bl(isab),ijcont)
      ENDIF
c---------------------------------------
      call secund(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c=======================
      subroutine prec2kl_new(ibl, bl)

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
c---------------------------------------
c1999
c get needed parameters from depository :
c
      call getival('inx_3' ,inx_3 )
      call getival('ibas_3',ibas_3)
      call getival('inx_4' ,inx_4 )
      call getival('ibas_4',ibas_4)
c
      call getival('npary',npary)
      call getival('ijbly',ijbly)
      call getival('mapijbly',map_ij_bly)
ctemp call getival('blockskl',nbl2)
      call getival('blockskl',nbl2_kl)
      call getival('blpredkl',nbl2_pred)
c
      nbl2_dim=nbl2_pred
c---------------------------------------
      if(icheck.gt.0) then
         call getmem(0,l0)
         return
      endif
c---------------------------------------
      call secund(tprec2b)
c---------------------------------------
      call dimenij_new(ibl,bl(inx_3),bl(inx_4),
     *                 bl(npary),bl(ijbly),nbl2_dim,
     *                 nparkl,kldim,klcont)    ! output
c---------------------------------------
      klpar1=nparkl
      IF( iroute.eq.1 ) THEN
         call getmem(3*kldim,kabprim)
         call ab_prim_1_new(ibl,bl(inuc),nparkl,bl(ijbly),nbl2_dim,
     *                      bl(ibas_3),bl(inx_3),
     *                      bl(ibas_4),bl(inx_4),
     *                      bl(kabprim),klcont)
c
      ELSE
         call getmem(2*klcont+kldim,kabprim)
         icpd =kabprim
         i1cpd=kabprim+klcont
         iscd =kabprim+klcont*2
c
         call ab_prim_2_new(ibl,bl(inuc),nparkl,bl(ijbly),nbl2_dim,
     *                      bl(ibas_3),bl(inx_3),
     *                      bl(ibas_4),bl(inx_4),
     *                      bl(icpd),bl(i1cpd),bl(iscd),klcont)
      ENDIF
c---------------------------------------
      call secund(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c============================================================
      subroutine precalc2_1_new(bl,mmax,mmax1,nhabcd,nfumax,
     *                       ibl,nijbeg,nijend,npij,
     *                       kbl,nklbeg,nklend,npklx,npkl)
      implicit real*8 (a-h,o-z)
      character*3 reuse
c-----------------------------------------------------------
      common /neglect/ eps,eps1,epsr,eps8
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
      dimension bl(*)
c---------------------------------------
c1999 get needed parameters from depository :
c
c for left side pairs :
c
      call getival('inx_1' ,inx_1 )
      call getival('ibas_1',ibas_1)
      call getival('inx_2' ,inx_2 )
      call getival('ibas_2',ibas_2)
      call getival('ijblx',ijblx)
      call getival('blocksij',nbl2_ij)
      call getival('blpredij',nbl2_dimij) ! for dim only
c
c for right side pairs :
c
      call getival('inx_3' ,inx_3 )
      call getival('ibas_3',ibas_3)
      call getival('inx_4' ,inx_4 )
      call getival('ibas_4',ibas_4)
      call getival('ijbly',ijbly)
      call getival('blockskl',nbl2_kl)
      call getival('blpredkl',nbl2_dimkl) ! for dim only
c---------------------------------------------------------------
      par268=eps8
c---------------------------------------------------------------
         reuse='no '
c---------------------------------------------------------------
      ijbegin=nijbeg
      klbegin=nklbeg
c---------------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c--------
      lpartot=lpartot+1
      if(reuse.eq.'no ') then
          lpareal=lpareal+1
          call precal2a_1_new(ibl,npij,nijbeg,nijend,bl(inuc),
     *                        bl(ijblx),nbl2_dimij,
     *                  bl(ibas_1),bl(inx_1),bl(ibas_2),bl(inx_2),
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
      if(npkl.eq.0) then
          kabprim=iabprim
          kldim  =ijdim
          klpar1 =ijpar1
      endif
c
         if(npkl.ne.0) then
            call precal2x_1(kabprim,kldim,icpd,i1cpd,iscd)
            lpartot=lpartot+1
            lpareal=lpareal+1
            call precal2a_1_new(kbl,npkl,nklbeg,nklend,bl(inuc),
     *                          bl(ijbly),nbl2_dimkl,
     *                  bl(ibas_3),bl(inx_3),bl(ibas_4),bl(inx_4),
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
      subroutine precalc2_2_new(bl,mmax,mmax1,       nfumax,
     *                          ibl,nijbeg,nijend,npij,
     *                          kbl,nklbeg,nklend,npkl)
      implicit real*8 (a-h,o-z)
      character*3 reuse
c-----------------------------------------------------------
      common /neglect/ eps,eps1,epsr,eps8
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
      dimension bl(*)
c---------------------------------------
c1999 get needed parameters from depository :
c
c for left side pairs :
c
      call getival('inx_1' ,inx_1 )
      call getival('ibas_1',ibas_1)
      call getival('inx_2' ,inx_2 )
      call getival('ibas_2',ibas_2)
      call getival('ijblx',ijblx)
      call getival('blocksij',nbl2_ij)
      call getival('blpredij',nbl2_dimij) ! for dim only
c
c for right side pairs :
c
      call getival('inx_3' ,inx_3 )
      call getival('ibas_3',ibas_3)
      call getival('inx_4' ,inx_4 )
      call getival('ibas_4',ibas_4)
      call getival('ijbly',ijbly)
      call getival('blockskl',nbl2_kl)
      call getival('blpredkl',nbl2_dimkl) ! for dim only
c---------------------------------------
c---------------------------------------------------------------
      par268=eps8
c---------------------------------------------------------------
         reuse='no '
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
            call precal2a_2_new(ibl,npij,nijbeg,nijend,bl(inuc),
     *                          bl(ijblx),nbl2_dimij,
     *                  bl(ibas_1),bl(inx_1),bl(ibas_2),bl(inx_2),
     *                  bl(iabprim+lcij),bl(iabprim+2*lcij),ijpar1,lcij,
     *                  bl(iaa),bl(ibb), bl(ieab), bl(icis),bl(icjs),
     *                  bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *                  bl(igcij),ngci1,ngcj1,'left ',par268)
c----------
c for abnia
          if(mmax.gt.2) then
             call precal2b_2(mmax1,lcij, bl(i1apb),bl(iabnia))
          endif
      endif
c----------------------------------------
c precalculations for pairs kl
c
      if( npkl.eq.0) then
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
            call precal2a_2_new(kbl,npkl,nklbeg,nklend,bl(inuc),
     *                          bl(ijbly),nbl2_dimkl,
     *                  bl(ibas_3),bl(inx_3),bl(ibas_4),bl(inx_4),
     *                  bl(kabprim+lckl),bl(kabprim+2*lckl),klpar1,lckl,
c not an error, icc should be twice :
c2002*  not anymore     bl(icc),bl(icc), bl(iecd), bl(icks),bl(icls),
     *                  bl(icc),bl(idd), bl(iecd), bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igckl),ngck1,ngcl1,'right',par268)
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
      subroutine dimenij_new(ibl,inx_1,inx_2,npar,ijbl,nbl2_dim,
     *                       nparij,ijdim,ijcont)
      dimension inx_1(12,*),inx_2(12,*)
      dimension npar(nbl2_dim),ijbl(nbl2_dim,*)
c
      nparij=npar(ibl)
      ijcs1=ijbl(ibl,1)
      call get_ij_half(ijcs1,ics1,jcs1)
      icont=inx_1(5,ics1)-inx_1(1,ics1)
      jcont=inx_2(5,jcs1)-inx_2(1,jcs1)
      ijcont=icont*jcont
      ijdim=nparij*ijcont
      end
c====================================================================
      subroutine ab_prim_2_new(ibl,datnuc,nparij,ijbl, nbl2_dim,
     *                         datbas_1,inx_1,
     *                         datbas_2,inx_2,
     *                         apb,rapb,sab,ijcont)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
      dimension datbas_1(13,*),datbas_2(13,*)
      dimension inx_1(12,*), inx_2(12,*), ijbl(nbl2_dim,*)
      dimension apb(ijcont),rapb(ijcont),sab(nparij,ijcont)
c
      dimension sqrtx(400),eexx(400)    ! local
      data zero,one /0.d0, 1.d0/
c---------------------------------------------------------------
c for the first pair :
c
      ijcs1=ijbl(ibl,  1  )
      call get_ij_half(ijcs1,ics1,jcs1)
      ia=inx_1(1,ics1)+1
      ie=inx_1(5,ics1)
      ja=inx_2(1,jcs1)+1
      je=inx_2(5,jcs1)
c
      ij=0
      do is=ia,ie
         aa=datbas_1(1,is)
         do js=ja,je
            bb=datbas_2(1,js)
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
      call get_ij_half(ijcs,ics,jcs)
      ia=inx_1(1,ics)+1
      ie=inx_1(5,ics)
      iat=inx_1(2,ics)
c
      ja=inx_2(1,jcs)+1
      je=inx_2(5,jcs)
      jat=inx_2(2,jcs)
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
c====================================================================
      subroutine first_pair_new(ijcs,datbas_1,inx_1, datbas_2,inx_2,
     *                          itypp,jtypp, par268,which,ngctot,
     *                          aexp,cis,bexp,cjs,
     *                          est_ijx,coef_ijx,factij,
     *                          gcij,ngcii,ngcjj)
      implicit real*8 (a-h,o-z)
      character*5 which
      dimension datbas_1(13,*), datbas_2(13,*)
      dimension inx_1(12,*),inx_2(12,*)
c
      dimension aexp(*),cis(*)
      dimension bexp(*),cjs(*)
c
      dimension est_ijx(30,30),coef_ijx(30,30)
      dimension factij(*)
      dimension gcij(ngcii+1,ngcjj+1,*)
c-----------------------------------------------------------------------
c input :
c ijcs - first shell-pair of interest
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
c
       call get_ij_half(ijcs, ics, jcs)
c
       ia=inx_1(1,ics)+1         ! beginning of the contraction
       ja=inx_2(1,jcs)+1
       ie=inx_1(5,ics)           ! end      of the contraction
       je=inx_2(5,jcs)
c
       is1=0
       do is=ia,ie
          is1=is1+1
          aa=datbas_1(1,is)
          aexp(is1)=aa
       enddo
       js1=0
       do js=ja,je
          js1=js1+1
          bb=datbas_2(1,js)
          bexp(js1)=bb
       enddo
c
       ji=0
       is1=0
       do is=ia,ie
          is1=is1+1
          if(ngctot.eq.0) then
             csi=datbas_1(2,is)
             cpi=datbas_1(3,is)
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
                csj=datbas_2(2,js)
                cpj=datbas_2(3,js)
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
                   gci=datbas_1(ig+2,is)
                   do jg=0,ngcjj
                     gcj=gci*datbas_2(jg+2,js)
                     gcij(ig+1,jg+1,ji)=gcj
                   enddo
                enddo
             endif
          enddo
       enddo
c
c-----------------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_2_new(ibl,npij,nijbeg,nijend,datnuc,
     *                          ijbl,nbl2_dim,
     *                          datbas_1,inx_1,datbas_2,inx_2,
     *                    rapb,sab,ijpar1, lcij,
     *                    aaa,bbb,    estab, cis,cjs,
     *                    xab,xparij,coefij,factij,txab,
     *                    gcij,ngci1,ngcj1,which,par268 )
c
      implicit real*8 (a-h,o-z)
      character*5 which
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
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
      dimension datbas_1(13,*),datbas_2(13,*), datnuc(5,*)
      dimension inx_1(12,*),inx_2(12,*),  ijbl(nbl2_dim,*)
c
      dimension cis(*),cjs(*)
      dimension xab(npij,3), xparij(npij,3,lcij,3)
      dimension estab (npij,lcij)
      dimension coefij(npij,lcij), factij(lcij)
      dimension xa(3),xb(3),txab(npij,3,*)
      dimension gcij(ngci1,ngcj1,lcij)
c local arrays (dimensioned as contraction max length) :
      dimension est_ijx(30,30),coef_ijx(30,30)
c
c
c only for where='forc' or 'hess'
      dimension aaa(     *),bbb(     *)   ! dimen= lci,lcj -contr.
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
c for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c-------------------------------------------------------------
c precalculate all possible quantities for the first pair only
c
      ijcs_first=ijbl(ibl,nijbeg)
      call first_pair_new(ijcs_first,datbas_1,inx_1, datbas_2,inx_2,
     *                    itypp,jtypp, par268,which,ngctot,
     *                    aaa ,cis,bbb ,cjs,
     *                    est_ijx,coef_ijx,factij,
     *                    gcij,ngcii,ngcjj)
c-------------------------------------------------------------
c precalculations for the pairs IJ :
c
      boamax=0.d0
      rapbmax=0.d0
      ijpar=0
      do 100 ijp=nijbeg,nijend
      ijpar=ijpar+1
        ijcs=ijbl(ibl,ijp)
        call get_ij_half(ijcs,ics,jcs)
        fact1=one
        if(ics.eq.jcs) fact1=fact1*half
        ia=inx_1(1,ics)+1      ! starting contr
        ja=inx_2(1,jcs)+1
        ie=inx_1(5,ics)        ! last contr
        je=inx_2(5,jcs)
c.......................................................
        iatom=inx_1(2,ics)
        jatom=inx_2(2,jcs)
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
            aa=datbas_1(1,is)
            aor1=1.d0/max(1.d0,aa)
c
         js1=0
         do 200 js=ja,je
            js1=js1+1
            ji=ji+1
            bb=datbas_2(1,js)
            boamax=max(boamax,bb*aor1)      ! b/max(1,a)
c
cnopermut   coefij(ijpar,ji)=coef_ijx(is1,js1)*fact1
            coefij(ijpar,ji)=coef_ijx(is1,js1)
            estab(ijpar,ji)=est_ijx(is1,js1)
c
            rapb1=rapb(ji)
            sab1 =sab(ijp,ji)
            aa1=aa*rapb1
            bb1=bb*rapb1
            rapbmax=max(rapbmax,rapb1)
c
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
c------------------------------------------------------------------
  200    continue        ! end of the loop over primitives
  100 continue           ! end of the loop over contrcted pairs
c------------------------------------------------------------------
c2002
      boamax=boamax+boamax
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
         aa1= aaa(is1)
         do 151 js1=1,lcjj
         bb1= bbb(js1)
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
         aa1= aaa(is1)
         do 152 js1=1,lcjj
         bb1= bbb(js1)
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
      subroutine ab_prim_1_new(ibl,datnuc,nparij,ijbl, nbl2_dim,
     *                         datbas_1,inx_1,
     *                         datbas_2,inx_2,
     *                         abprim,ijcont)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datnuc(5,*)
      dimension datbas_1(13,*),datbas_2(13,*)
      dimension inx_1(12,*), inx_2(12,*), ijbl(nbl2_dim,*)
      dimension abprim(nparij,ijcont,3)
c---------------------------------------------------------------
      do 100 ijpar=1,nparij
         ijcs=ijbl(ibl,ijpar)
         call get_ij_half(ijcs,ics,jcs)
         ia=inx_1(1,ics)+1
         ie=inx_1(5,ics)
         iat=inx_1(2,ics)
c
         ja=inx_2(1,jcs)+1
         je=inx_2(5,jcs)
         jat=inx_2(2,jcs)
c
         xab=datnuc(2,iat)-datnuc(2,jat)
         yab=datnuc(3,iat)-datnuc(3,jat)
         zab=datnuc(4,iat)-datnuc(4,jat)
         rr=xab*xab+yab*yab+zab*zab
         ij=0
         do 200 is=ia,ie
            aa=datbas_1(1,is)
            do 200 js=ja,je
               bb=datbas_2(1,js)
               ij=ij+1
               axb=aa*bb
               apb=aa+bb
               apb1=one/apb
               e=axb*apb1
               err=e*rr
               if(err.gt.32) then
                  exp_err=0.d0         ! exp(-32)=1.27*10**-14
               else
                  exp_err=exp(-err)
               endif
c98
               abprim(ijpar,ij,1)=apb
               abprim(ijpar,ij,2)=apb1
c98            abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp(-e*rr)
               abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp_err
  200       continue
  100 continue
c---------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_1_new(ibl,npij,nijbeg,nijend,datnuc,
     *                          ijbl,nbl2_dim,
     *                          datbas_1,inx_1,datbas_2,inx_2,
     *                          abprim,ijpar1,lcij,
     *                          aaa,bbb,estab,cis,cjs,
     *                          xab,xparij,coefij,factij,txab,
     *                          gci,gcj,ngci1,ngcj1,which,par268)
c
      implicit real*8 (a-h,o-z)
      character*5 which
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
      dimension abprim(ijpar1,lcij,3)
c
      dimension datbas_1(13,*),datbas_2(13,*), datnuc(5,*)
      dimension inx_1(12,*),inx_2(12,*),  ijbl(nbl2_dim,*)
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
c for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c---------------
c precalculations for the pairs IJ :
c
      boamax=0.d0
      rapbmax=0.d0
      ijpar=0
      do 100 ijp=nijbeg,nijend
         ijpar=ijpar+1
         ijcs=ijbl(ibl,ijp)
         call get_ij_half(ijcs,ics,jcs)
c
cnoper   fact1=one
cnoper   if(ics.eq.jcs) fact1=fact1*half
c
c       starting contr
         ia=inx_1(1,ics)+1
         ja=inx_2(1,jcs)+1
c       last contr
         ie=inx_1(5,ics)
         je=inx_2(5,jcs)
c
c       number of general contr.
c        ngci=inx_1(4,ics)
c        ngcj=inx_2(4,jcs)
c in the common block contr
c
         iatom=inx_1(2,ics)
         jatom=inx_2(2,jcs)
c
c97
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
  150    continue
c
         ji=0
         is1=0
         do 200 is=ia,ie
            is1=is1+1
            aa=datbas_1(1,is)
            aor1=1.d0/max(1.d0,aa)
            aaa(ijpar,is1)=aa
            if(ngctot.eq.0) then
               csi=datbas_1(2,is)
               cpi=datbas_1(3,is)
               coefi=csi
               est_i=max( abs(csi),abs(cpi) )
               if(itypp.eq.3) then
                  coefi=cpi
                  csi=csi/cpi
                  facti=csi
               endif
cnopermut      coefi=coefi*fact1
               if(which.eq.'right') coefi=coefi*par268
               cis(ijpar,is1)=csi
            else
c              gen.contr. shell is somewhere
               est_i=one
               do 210 ig=0,ngcii
                  gci(ijpar,ig+1,is1)=datbas_1(ig+2,is)
  210          continue
            endif
c
         js1=0
         do 200 js=ja,je
            js1=js1+1
            ji=ji+1
            bb=datbas_2(1,js)
            boamax=max(boamax,bb*aor1)      ! b/max(1,a)
            bbb(ijpar,js1)=bb
            if(ngctot.eq.0) then
               csj=datbas_2(2,js)
               cpj=datbas_2(3,js)
               coefj=csj
               if(jtypp.eq.3) coefj=cpj
               coefij(ijpar,ji )=coefi*coefj
ckw99n
               est_j=max( abs(csj),abs(cpj) )
               estab(ijpar,ji)=est_i*est_j
ckw99
c
               if(jtypp.eq.3) then
                  csj=csj/cpj
                  factj=csj
                  if(itypp.eq.3) factij(ijpar,ji  )=facti*factj
               endif
               cjs(ijpar,js1)=csj
            else
c              gen.contr.
cnopermut      coefij(ijpar,ji)=fact1
               coefij(ijpar,ji)=1.d0
               estab(ijpar,ji)=1.d0
cnopermut      if(which.eq.'right') coefij(ijpar,ji)=fact1*par268
               if(which.eq.'right') coefij(ijpar,ji)=par268
c
               do 220 jg=0,ngcjj
                  gcj(ijpar,jg+1,js1)=datbas_2(jg+2,js)
  220          continue
            endif
c--------------------------------
            rapb=abprim(ijp,ji,2)
            sab =abprim(ijp,ji,3)
c
            aa1=aa*abprim(ijp,ji,2)
            bb1=bb*abprim(ijp,ji,2)
c
            coefij(ijpar,ji)=coefij(ijpar,ji)     *sab
            estab(ijpar,ji)=estab(ijpar,ji)*sab
c
            xpn_max=1.d0
            do 230 l=1,3
               xparij(ijpar,l,ji,1)=aa1*xa(l)+bb1*xb(l)  !xp(ijpar,l,ji
               xxl=xa(l)
               if(nqii.lt.nqjj) xxl=xb(l)
               xparij(ijpar,l,ji,2)=xparij(ijpar,l,ji,1)-xxl  ! xpn
               xparij(ijpar,l,ji,3)=aa*xa(l)+bb*xb(l)         ! xpp
c              for neglect :
               xpn_abs=abs( xparij(ijpar,l,ji,2) )
               if(xpn_abs.gt.xpn_max) xpn_max=xpn_abs
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
            rapbmax=max(rapbmax,rapb)
c------------------------------------------------------------------
  200    continue       ! end of the loop over contractions
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
      subroutine prec4neg_2_n(nbls,npij,npkl,ndiag,ij,kl,
     1     lc12,lc34,indxij,indxkl,indxr,
     *     estab,estcd,densmax, esti2ij,esti2kl,
     *     list_ij,list_kl,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, nsym,isymm,
coutput
     *     rho,rppq,rhoapb,rhocpd,rys,const, nbls1,index)
c-------------------------------------------------------------------
c Input :
c   esti2ij(ij) - maximum estim for ij-prim.pair
c   esti2kl(kl) - maximum estim for kl-prim.pair * max of densmax
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical jump(*)
      common /neglect/ eps,eps1,epsr,eps8
c
      dimension indxij(*),indxkl(*),indxr(*),index(*)
      dimension apb(lc12), cpd(lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension  esti2ij(lc12), esti2kl(lc34)
      dimension  list_ij(npij),list_kl(npkl)
      dimension  densmax(nbls)
      dimension  symfac(*)
c-------
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension rys(*),const(*)
      dimension isymm(*)
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
c---------------------------------------------------------------
c     permut=half
c---------------------------------------------------------------
      apb1=apb(ij)
      cpd1=cpd(kl)
c
      abpcd1=apb1+cpd1
      abxcd=apb1*cpd1
      abpcdr=one/abpcd1
c-------------------------------------------
      ddmax1=esti2kl(kl)*abpcdr
c-------------------------------------------
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
      IF(ndiag.eq.0) then
         sqrpold=sqrpold*eps8
         ijkl=0
         ijkl1=0
         IF(NSYM.EQ.0) THEN
         do 100 ijpar=1,npij
            npklx=ijpar
            esti1=estab(ijpar,ij)
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 100
            endif
            coef1=coefij(ijpar,ij)
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do  50 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 50
               esti2=estcd(klpar,kl)
               estim=esti1*esti2*abpcdr*densmax(ijkl)
               if(estim.gt.epsr) then
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
                  coef2=coefkl(klpar,kl)
                  const(ijkl1)=coef1*coef2*sqrpold
cnopermut         if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut            const(ijkl1)=const(ijkl1)*permut
cnopermut         endif
               endif
  50        continue
 100     continue
c
         ELSE
c ..........symmetry...................
c
         do 200 ijpar=1,npij
            npklx=ijpar
            if(jump(ijpar)) then
               ijkl=ijkl+npklx
               go to 200
            endif
            esti1=estab(ijpar,ij)
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 200
            endif
            coef1=coefij(ijpar,ij)
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 150 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 150
               esti2=estcd(klpar,kl)
               estim=esti1*esti2*abpcdr*densmax(ijkl)
               if(estim.gt.epsr) then
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
                  coef2=coefkl(klpar,kl)
                  const(ijkl1)=coef1*coef2*sqrpold
cnopermut         if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut            const(ijkl1)=const(ijkl1)*permut
cnopermut         endif
c-------->        const(ijkl1)=const(ijkl1)*symfac
                  const(ijkl1)=const(ijkl1)*symfac(ijkl)
               endif
 150        continue
 200     continue
         ENDIF
         nbls1=ijkl1
      ENDIF
c
      IF(ndiag.ne.0) then
         ijkl=0
         ijkl1=0
         IF(NSYM.EQ.0) THEN
         do 300 ijpar=1,npij
ckw      do 300 ijp  =1,npij
ckw         ijpar=list_ij(ijp)
            ijij=(ijpar-1)*npkl
            npklx=npkl
            esti1=estab(ijpar,ij)
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 300
ckw            EXIT
            endif
            coef1=coefij(ijpar,ij)
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 350 klpar=1,npklx
ckw         do 350 klp  =1,npklx
ckw            klpar=list_kl(klp)
c............. ijkl=ijkl+1
               ijkl=ijij+klpar
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 350
               esti2=estcd(klpar,kl)
               estim=esti1*esti2*abpcdr*densmax(ijkl)
               if(estim.ge.epsr) then
ckw            if(estim.lt.epsr) EXIT
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
                  coef2=coefkl(klpar,kl)
                  const(ijkl1)=coef1*coef2*sqrpold
cnopermut         if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut            const(ijkl1)=const(ijkl1)*permut
cnopermut         endif
               endif
 350        continue
 300     continue
         ENDIF
         nbls1=ijkl1
      ENDIF
c
      end
c====================================================================
      subroutine precspec_2_n(nbls,npij,npkl,ndiag, ij,kl,
     1     lc12,lc34,indxij,indxkl,indxr,
cccc *     estab,estcd,densmax, esti2kl,
     *     estab,estcd,densmax, esti2ij,esti2kl,
     *     list_ij,list_kl,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, rapb,rcpd,txab,txcd,
     *     nsym,isymm,
     *     rho,rys,const,xpqr,txxr, nbls1,index)
c--------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical jump(*)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr,eps8
c
      dimension indxij(*),indxkl(*),indxr(*),index(*)
c---------------
      dimension apb(lc12),cpd(lc34)
      dimension rapb(lc12),rcpd(lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
ccccc dimension                    esti2kl(lc34)
      dimension  esti2ij(lc12), esti2kl(lc34)
      dimension  list_ij(npij),list_kl(npkl)
      dimension  densmax(nbls)
      dimension  symfac(*)
c-------
      dimension txab(npij,3,*),txcd(npkl,3,*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension isymm(*)
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
c---------------------------------------------------------------
ccc   par268=eps1*8.d0
ccc   permut=half
c---------------------------------------------------------------
      apb1=apb(ij)
      cpd1=cpd(kl)
c
      abpcd1=apb1+cpd1
      abxcd=apb1*cpd1
      abpcdr=one/abpcd1
c-------------------------------------------
      ddmax1=esti2kl(kl)*abpcdr
c-------------------------------------------
      rho1=abxcd*abpcdr
      rho2=rho1
      if(ityp.gt.1.or.jtyp.gt.1) rho2=rho2*rapb(ij)
      if(ktyp.gt.1.or.ltyp.gt.1) rho2=rho2*rcpd(kl)
c
      abpcdrcpd1=abpcdr*cpd1
      abpcdrapb1=abpcdr*apb1
      sqrpold=sqrt(abpcdr)
c-------------------------------------------
      IF(ndiag.eq.0) then
         sqrpold=sqrpold*eps8
         ijkl=0
         ijkl1=0
         IF(NSYM.EQ.0) THEN
         do 100 ijpar=1,npij
            npklx=ijpar
            esti1=estab(ijpar,ij)
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 100
            endif
            coef1=coefij(ijpar,ij)
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 50 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 50
               esti2=estcd(klpar,kl)
               estim=esti1*esti2*abpcdr*densmax(ijkl)
               if(estim.gt.epsr) then
                  coef2=coefkl(klpar,kl)
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
                  xpqr(1,ijkl1)=x1*rho2
                  xpqr(2,ijkl1)=x2*rho2
                  xpqr(3,ijkl1)=x3*rho2
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
                  const(ijkl1)=coef1*coef2*sqrpold
cnopermut         if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut            const(ijkl1)=const(ijkl1)*permut
cnopermut         endif
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
            esti1=estab(ijpar,ij)
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 200
            endif
            coef1=coefij(ijpar,ij)
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 150 klpar=1,npklx
               ijkl=ijkl+1
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 150
               esti2=estcd(klpar,kl)
               estim=esti1*esti2*abpcdr*densmax(ijkl)
               if(estim.gt.epsr) then
                  coef2=coefkl(klpar,kl)
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
                  xpqr(1,ijkl1)=x1*rho2
                  xpqr(2,ijkl1)=x2*rho2
                  xpqr(3,ijkl1)=x3*rho2
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
                  const(ijkl1)=coef1*coef2*sqrpold
cnopermut         if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut            const(ijkl1)=const(ijkl1)*permut
cnopermut         endif
                  const(ijkl1)=const(ijkl1)*symfac(ijkl)
                endif
 150        continue
 200     continue
         ENDIF
         nbls1=ijkl1
c......................................
      ENDIF
      IF(ndiag.NE.0) then
         ijkl=0
         ijkl1=0
         IF(NSYM.EQ.0) THEN
         do 300 ijpar=1,npij
ckw      do 300 ijp  =1,npij
ckw         ijpar=list_ij(ijp)
            ijij=(ijpar-1)*npkl
            npklx=npkl
            esti1=estab(ijpar,ij)
            if(esti1*ddmax1.lt.epsr) then
               ijkl=ijkl+npklx
               go to 300
ckw            EXIT
            endif
            coef1=coefij(ijpar,ij)
            xp1= xp(ijpar,1,ij)
            xp2= xp(ijpar,2,ij)
            xp3= xp(ijpar,3,ij)
            do 350 klpar=1,npklx
ckw         do 350 klp  =1,npklx
c............  ijkl=ijkl+1
ckw            klpar=list_kl(klp)
               ijkl=ijij+klpar
               ijklsm=isymm(ijkl)
               if(ijklsm.eq.0) go to 350
               esti2=estcd(klpar,kl)
               estim=esti1*esti2*abpcdr*densmax(ijkl)
               if(estim.ge.epsr) then
ckw            if(estim.lt.epsr) EXIT
                  coef2=coefkl(klpar,kl)
                  ijkl1=ijkl1+1
                  index(ijkl1)=indxr(ijkl)
                  x1= xp1 - xq(klpar,1,kl)
                  x2= xp2 - xq(klpar,2,kl)
                  x3= xp3 - xq(klpar,3,kl)
                  xpqr(1,ijkl1)=x1*rho2
                  xpqr(2,ijkl1)=x2*rho2
                  xpqr(3,ijkl1)=x3*rho2
                  rr2=x1*x1 + x2*x2 + x3*x3
                  rys(ijkl1)=rr2*rho1
                  const(ijkl1)=coef1*coef2*sqrpold
cnopermut         if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut            const(ijkl1)=const(ijkl1)*permut
cnopermut         endif
               endif
 350        continue
 300     continue
         ENDIF
         nbls1=ijkl1
      ENDIF
c--------------------------------------------------------------------
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
c
      end
c====================================================================
      subroutine prec4neg_1_n(nbls,npij,npkl,ndiag,ij,kl,
     1     ijpar1,lc12, klpar1,lc34,indxij,indxkl,
     2     estab,estcd,densmax,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, nsym,isymm,
coutput
     *     rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical jump(*)
c
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr,eps8
      common /begin/ ijbegin,klbegin
c
      dimension indxij(*),indxkl(*),index(*)
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension  densmax(nbls)
c-------
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension        rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c symmetry
      dimension isymm(*)
      dimension symfac(*)
c
      data zero,half,one /0.d0 , 0.5d0 , 1.d0 /
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
            coef1=coef1*eps8
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
c------->      index(ijkl1)=ijkl
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
cnopermut      if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut         const(ijkl1)=const(ijkl1)*permut
cnopermut      endif
               if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac(ijkl)
            endif
 150     continue
 100  continue
c
      nbls1=ijkl1
c
      end
c====================================================================
      subroutine precspec_1_n(nbls,npij,npkl,ndiag, ij,kl,
     1     ijpar1,lc12, klpar1,lc34, indxij,indxkl,
     2     estab,estcd,densmax,
     *     jump,symfac,
     *     apb,cpd,coefij,coefkl,xp,xq, rapb,rcpd,txab,txcd,
     *                   nsym,isymm,
c output
     *     rys,const,xpqr,txxr,nbls1,index)
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical jump(*)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr,eps8
      common /begin/ ijbegin,klbegin
c
      dimension indxij(*),indxkl(*),index(*)
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension rapb(ijpar1,*),rcpd(klpar1,*)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension  densmax(nbls)
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
            coef1=coef1*eps8
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
            estim=esti1*esti2*abpcdr*densmax(ijkl)
            if(estim.gt.epsr) then
               coef12=coef1*coef2
               ijkl1=ijkl1+1
c---->         index(ijkl1)=ijkl
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
cnopermut      if(ndiag.eq.0.and.ijpar.eq.klpar) then
cnopermut         const(ijkl1)=const(ijkl1)*permut
cnopermut      endif
               if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac(ijkl)
c
            endif
 150     continue
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
           xpqr(1,i)=xpqr(1,i)*rapb(ijstart+ijpar,ij)
           xpqr(2,i)=xpqr(2,i)*rapb(ijstart+ijpar,ij)
           xpqr(3,i)=xpqr(3,i)*rapb(ijstart+ijpar,ij)
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
           txxr(1,i)=txcd(klpar,1,kl)
           txxr(2,i)=txcd(klpar,2,kl)
           txxr(3,i)=txcd(klpar,3,kl)
  220   continue
      endif
c
      end
c====================================================================
