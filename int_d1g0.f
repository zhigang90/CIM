c=====================================================================
      subroutine int_d1g0(idft,ax,nblocks,bl,inx,ntri,thres1,
     *                        dens,fock,labels,mywork,igran)
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
cxxx  common /datadisk/ ?????????????????????????????????????
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
c     common /intbl/ifpp,inxx(100)
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),dens(ntri,3),fock(ntri,3)
      dimension labels(*)
c----------------------------------------------------------------
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
c          in-core storage :
           call fock_core_cphf(bl(incorex),
     *                   bl(ilab),bl(jlab),bl(klab),bl(llab),
     *                   bl(isiz),bl(jsiz),bl(ksiz),bl(lsiz),bl(iqrt),
     *                         integrals,nqstore,nblstore,
     *                         thres0,thres1, bl(lind),dens,fock,ntri)
        else
c          on-disk storage :
c          call fock_disk_cphf( ?????????????????????????????????
c    *                          thres0,thres1, bl(lind),dens,fock)
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
c----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c re-calculate integrals
c
      where='fock'
      screen='fock'
c
c----------------------------------------------------------------
c allocate memory for denspar(ics,jcs) (screening density)
c
      call getmem(ncs*ncs,idensp)
      call getmem(ncf,map_fs)      ! mapping from funct. to shells
c----------------------------------------------------------------
c there are 3 density matrices; select one (max) for calcint2 call
      call getmem(ntri,lselect)
c
      call select_dens(dens,ntri, bl(lselect))
c----------------------------------------------------------------
c transform selected density bl(lselect) dens(ij) into denspar(ics,jcs)
c
      call setup_densp2(inx,ncf,ncs,bl(lselect),bl(idensp),
     *                 bl(map_fs))
c----------------------------------------------------------------
      call retmem(1)    ! lselect
c----------------------------------------------------------------
c set inetgral timings to zero as well as neglect stuff
c
      call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c symmetry stuff :
c
      call getival('nsym',nsym)
      if(nsym.gt.0) call getival('SymFunPr',ifp)
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
c
      do isupbl=istart,istop
c
c        get integral's price
c
         call get_price(bl(ncost),isupbl,iprice)
c
c blocks were stored:
c  - if they have negative price and they are before the last stored
c  - if they have negative price and all negative blocks were stored
c  - if all negative blocks were stored and our block is before the last
c
         if((iprice.LE.0.and.isupbl.le.ito).or.
     *      (iprice.le.0.and. isecond.eq.1).or.
     *      (isecond.eq.1.and.isupbl.le.ito)) go to 1111
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,bl(idensp),where,
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
           call fock_d1g0_nos(idft,ax,
     *                         bl(ibuffz),dens,fock,bl(lind),ntri,
     *                         bl(idensp),bl(map_fs),ncs,
     *                         nblsiz,ngctoz,nintez,thres1,
     *                         labels(lab1),labels(lab2),labels(lab3))
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
c release memory allocated for
c
      call retmem(1)    ! map_fs
      call retmem(1)    ! idensp
c----------------------------------------------------------------
 8888 continue
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
c   2. using stored integrals (stored in-corer  : fock_core_cphf
c   3. using stored integrals (stored on-diskr  : fock_disk_cphf
c
c=====================================================================
      subroutine fock_d1g0_nos(idft,ax,buf,dens,fock,lind,ntri,
     *                         densp,map_fs,ncs,
     *                  nbls,ngcd,lnijkl,thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for polarizability CPHF )
c
c NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension dens(ntri,3),fock(ntri,3)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
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
cccc    do 150 iqu=1,ngcd
           icff=labels(1,iqu,ijklp)
           jcff=labels(2,iqu,ijklp)
           kcff=labels(3,iqu,ijklp)
           lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c shell :
           ics=map_fs(icff+1)
           jcs=map_fs(jcff+1)
           kcs=map_fs(kcff+1)
           lcs=map_fs(lcff+1)
c
c screening density (squared):
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
           dik=densp(ics,kcs)
           dil=densp(ics,lcs)
           djk=densp(jcs,kcs)
           djl=densp(jcs,lcs)
c
           dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c
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
       do 200 lll=1,llen
          lcf=lcff+lll
          ll=lind(lcf)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
c         magnitude checking again
c
          if(dmax*abs(xint0).gt.0.5d0) then
c
             xint4=xint0*4.d0
             if(idft.gt.0) xint0=xint0*ax
c.......................................................
c ....................Coulomb  part ....................
                if(icf.ge.jcf) then
                  ijf = ii+jcf
                else
                  ijf = jj+icf
                endif
                if(kcf.ge.lcf) then
                  klf = kk+lcf
                else
                  klf = ll+kcf
                endif
c ....................Exchange part ....................
                if(icf.ge.kcf) then
                  ikf = ii+kcf
                else
                  ikf = kk+icf
                endif
                if(icf.ge.lcf) then
                  ilf = ii+lcf
                else
                  ilf = ll+icf
                endif
                if(jcf.ge.kcf) then
                  jkf = jj+kcf
                else
                  jkf = kk+jcf
                endif
                if(jcf.ge.lcf) then
                  jlf = jj+lcf
                else
                  jlf = ll+jcf
                endif
c..........................................................
c.....................Coulomb  part ....................
c  ***  ij  ***
                fock(ijf,1)=fock(ijf,1) + dens(klf,1)*xint4
                fock(ijf,2)=fock(ijf,2) + dens(klf,2)*xint4
                fock(ijf,3)=fock(ijf,3) + dens(klf,3)*xint4
c  ***  kl  ***
                fock(klf,1)=fock(klf,1) + dens(ijf,1)*xint4
                fock(klf,2)=fock(klf,2) + dens(ijf,2)*xint4
                fock(klf,3)=fock(klf,3) + dens(ijf,3)*xint4
c
             if(idft.eq.0 .or. (idft.gt.0.and.ax.ne.0.d0)) then
c ....................Exchange part ....................
c  ***  il  ***
                fock(ilf,1)=fock(ilf,1) - dens(jkf,1)*xint0
                fock(ilf,2)=fock(ilf,2) - dens(jkf,2)*xint0
                fock(ilf,3)=fock(ilf,3) - dens(jkf,3)*xint0
c  ***  jl  ***
                fock(jlf,1)=fock(jlf,1) - dens(ikf,1)*xint0
                fock(jlf,2)=fock(jlf,2) - dens(ikf,2)*xint0
                fock(jlf,3)=fock(jlf,3) - dens(ikf,3)*xint0
c  ***  ik  ***
                fock(ikf,1)=fock(ikf,1) - dens(jlf,1)*xint0
                fock(ikf,2)=fock(ikf,2) - dens(jlf,2)*xint0
                fock(ikf,3)=fock(ikf,3) - dens(jlf,3)*xint0
c  ***  jk  ***
                fock(jkf,1)=fock(jkf,1) - dens(ilf,1)*xint0
                fock(jkf,2)=fock(jkf,2) - dens(ilf,2)*xint0
                fock(jkf,3)=fock(jkf,3) - dens(ilf,3)*xint0
c..........................................................
             endif
c
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
c
c Not used anymore, 
c Fder is symetriezed from the skeleton in fdersymm_gen
c
      subroutine fock_d1g0_sym(idft,ax,nsym,ifp,
     *                         buf,dens,fock,lind,ntri,
     *                         densp,map_fs,ncs,
     *                    nbls,ngcd,lnijkl,thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c SYMMETRY
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension dens(ntri,3),fock(ntri,3)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
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
ccc     do 150 iqu=1,ngcd
           icff=labels(1,iqu,ijklp)
           jcff=labels(2,iqu,ijklp)
           kcff=labels(3,iqu,ijklp)
           lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c shell :
           ics=map_fs(icff+1)
           jcs=map_fs(jcff+1)
           kcs=map_fs(kcff+1)
           lcs=map_fs(lcff+1)
c
c screening density (squared):
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
           dik=densp(ics,kcs)
           dil=densp(ics,lcs)
           djk=densp(jcs,kcs)
           djl=densp(jcs,lcs)
c
           dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
          integ=0
       do 200 iii=1,ilen
          ic1=icff+iii
          ii1=lind(ic1)
       do 200 jjj=1,jlen
          jc1=jcff+jjj
          jj1=lind(jc1)
       do 200 kkk=1,klen
          kc1=kcff+kkk
          kk1=lind(kc1)
       do 200 lll=1,llen
          lc1=lcff+lll
          ll1=lind(lc1)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
c-------------------------------------------------------
          if(dmax*abs(xint0).gt.0.5d0) then
c
             xint4=xint0*4.d0
             if(idft.gt.0) xint0=xint0*ax
c
c..............................................................
c ....................Coulomb  part ....................
                if(ic1.ge.jc1) then
                  ijf = ii1+jc1
                else
                  ijf = jj1+ic1
                endif
                if(kc1.ge.lc1) then
                  klf = kk1+lc1
                else
                  klf = ll1+kc1
                endif
c ....................Exchange part ....................
                if(ic1.ge.kc1) then
                  ikf = ii1+kc1
                else
                  ikf = kk1+ic1
                endif
                if(ic1.ge.lc1) then
                  ilf = ii1+lc1
                else
                  ilf = ll1+ic1
                endif
                if(jc1.ge.kc1) then
                  jkf = jj1+kc1
                else
                  jkf = kk1+jc1
                endif
                if(jc1.ge.lc1) then
                  jlf = jj1+lc1
                else
                  jlf = ll1+jc1
                endif
c..........................................................
c.....................Coulomb  part ....................
c  ***  ij  ***
                fock(ijf,1)=fock(ijf,1) + dens(klf,1)*xint4
                fock(ijf,2)=fock(ijf,2) + dens(klf,2)*xint4
                fock(ijf,3)=fock(ijf,3) + dens(klf,3)*xint4
c  ***  kl  ***
                fock(klf,1)=fock(klf,1) + dens(ijf,1)*xint4
                fock(klf,2)=fock(klf,2) + dens(ijf,2)*xint4
                fock(klf,3)=fock(klf,3) + dens(ijf,3)*xint4
c
c ....................Exchange part ....................
c  ***  il  ***
                fock(ilf,1)=fock(ilf,1) - dens(jkf,1)*xint0
                fock(ilf,2)=fock(ilf,2) - dens(jkf,2)*xint0
                fock(ilf,3)=fock(ilf,3) - dens(jkf,3)*xint0
c  ***  jl  ***
                fock(jlf,1)=fock(jlf,1) - dens(ikf,1)*xint0
                fock(jlf,2)=fock(jlf,2) - dens(ikf,2)*xint0
                fock(jlf,3)=fock(jlf,3) - dens(ikf,3)*xint0
c  ***  ik  ***
                fock(ikf,1)=fock(ikf,1) - dens(jlf,1)*xint0
                fock(ikf,2)=fock(ikf,2) - dens(jlf,2)*xint0
                fock(ikf,3)=fock(ikf,3) - dens(jlf,3)*xint0
c  ***  jk  ***
                fock(jkf,1)=fock(jkf,1) - dens(ilf,1)*xint0
                fock(jkf,2)=fock(jkf,2) - dens(ilf,2)*xint0
                fock(jkf,3)=fock(jkf,3) - dens(ilf,3)*xint0
c..............................................................
c
             do ns=1,nsym
c
                xcoul=xint4
                xexch=xint0
c
                icf=ifp(ns,ic1)
                jcf=ifp(ns,jc1)
                kcf=ifp(ns,kc1)
                lcf=ifp(ns,lc1)
c
                if(icf.lt.0) then
c                  write(6,*)' icf < 0 =',icf
                   icf=-icf
                   xcoul=-xcoul
                   xexch=-xexch
                endif
                if(jcf.lt.0) then
c                  write(6,*)' jcf < 0 =',jcf
                   jcf=-jcf
                   xcoul=-xcoul
                   xexch=-xexch
                endif
                if(kcf.lt.0) then
c                  write(6,*)' kcf < 0 =',kcf
                   kcf=-kcf
                   xcoul=-xcoul
                   xexch=-xexch
                endif
                if(lcf.lt.0) then
c                  write(6,*)' lcf < 0 =',lcf
                   lcf=-lcf
                   xcoul=-xcoul
                   xexch=-xexch
                endif
c
                 ii=lind(icf)
                 jj=lind(jcf)
                 kk=lind(kcf)
                 ll=lind(lcf)
c..............................................................
c ....................Coulomb  part ....................
                if(icf.ge.jcf) then
                  ijf = ii+jcf
                else
                  ijf = jj+icf
                endif
                if(kcf.ge.lcf) then
                  klf = kk+lcf
                else
                  klf = ll+kcf
                endif
c ....................Exchange part ....................
                if(icf.ge.kcf) then
                  ikf = ii+kcf
                else
                  ikf = kk+icf
                endif
                if(icf.ge.lcf) then
                  ilf = ii+lcf
                else
                  ilf = ll+icf
                endif
                if(jcf.ge.kcf) then
                  jkf = jj+kcf
                else
                  jkf = kk+jcf
                endif
                if(jcf.ge.lcf) then
                  jlf = jj+lcf
                else
                  jlf = ll+jcf
                endif
c..........................................................
c.....................Coulomb  part ....................
c  ***  ij  ***
                fock(ijf,1)=fock(ijf,1) + dens(klf,1)*xcoul
                fock(ijf,2)=fock(ijf,2) + dens(klf,2)*xcoul
                fock(ijf,3)=fock(ijf,3) + dens(klf,3)*xcoul
c  ***  kl  ***
                fock(klf,1)=fock(klf,1) + dens(ijf,1)*xcoul
                fock(klf,2)=fock(klf,2) + dens(ijf,2)*xcoul
                fock(klf,3)=fock(klf,3) + dens(ijf,3)*xcoul
c
c ....................Exchange part ....................
c  ***  il  ***
                fock(ilf,1)=fock(ilf,1) - dens(jkf,1)*xexch
                fock(ilf,2)=fock(ilf,2) - dens(jkf,2)*xexch
                fock(ilf,3)=fock(ilf,3) - dens(jkf,3)*xexch
c  ***  jl  ***
                fock(jlf,1)=fock(jlf,1) - dens(ikf,1)*xexch
                fock(jlf,2)=fock(jlf,2) - dens(ikf,2)*xexch
                fock(jlf,3)=fock(jlf,3) - dens(ikf,3)*xexch
c  ***  ik  ***
                fock(ikf,1)=fock(ikf,1) - dens(jlf,1)*xexch
                fock(ikf,2)=fock(ikf,2) - dens(jlf,2)*xexch
                fock(ikf,3)=fock(ikf,3) - dens(jlf,3)*xexch
c  ***  jk  ***
                fock(jkf,1)=fock(jkf,1) - dens(ilf,1)*xexch
                fock(jkf,2)=fock(jkf,2) - dens(ilf,2)*xexch
                fock(jkf,3)=fock(jkf,3) - dens(ilf,3)*xexch
c..........................................................
c..............................................................
             enddo !  ns=1,nsym
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
