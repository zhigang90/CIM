      subroutine int_d1g0nat(
     $        natend,      natonce,    idft,       ax,
     $        nblocks,     bl,         inx,        ntri,
     $        thres1,      map_fs,     denspar,    dens,
     $        fock,        labels,     mywork,     igran)
c---------------------------------------------------------------------
c for NATONCE atoms at one : cphf for hessian
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
      logical rescale
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
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*)
      dimension denspar(*),dens(*),fock(*)
      dimension natend(natonce)
      dimension labels(*),map_fs(*)
c
      Parameter (Zero=0.0d0)
c----------------------------------------------------------------
      call secund(time0)
c----------------------------------------------------------------
c symmetry stuff :
c
      call getival('nsym',nsym)
      if(nsym.gt.0) call getival('SymFunPr',ifp)
c----------------------------------------------------------------
c
c  if it is a full-direct run, reset common datstore
c  and skip storage retrieval
c
      if(scftype.eq.'full-direct')then
        isto=0
        ito=0
        isecond=0
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
          If(nsym.eq.0) Then
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nd1g0_core_nos(
     $          bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $          bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $          bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $          nblstore,    thres0,     thres1,     bl(lind),
     $          dens,        fock,       ntri,       natonce,
     $          idft,        ax,         denspar,    ncs,
     $          map_fs)
            else
              call nd1g0_core_nosC(
     $          bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $          bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $          bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $          nblstore,    thres0,     thres1,     bl(lind),
     $          dens,        fock,       ntri,       natonce,
     $          denspar,     ncs,        map_fs)
            endif
          else
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nd1g0_core_sym(
     $          nsym,        bl(ifp),    bl(incorex),bl(ilab),
     $          bl(jlab),    bl(klab),   bl(llab),   bl(isiz),
     $          bl(jsiz),    bl(ksiz),   bl(lsiz),   bl(iqrt),
     $          integrals,   nqstore,    nblstore,   thres0,
     $          thres1,      bl(lind),   dens,       fock,
     $          ntri,        natonce,    idft,       ax,
     $          denspar,     ncs,        map_fs)
             else
              call nd1g0_core_symC(
     $          nsym,        bl(ifp),    bl(incorex),bl(ilab),
     $          bl(jlab),    bl(klab),   bl(llab),   bl(isiz),
     $          bl(jsiz),    bl(ksiz),   bl(lsiz),   bl(iqrt),
     $          integrals,   nqstore,    nblstore,   thres0,
     $          thres1,      bl(lind),   dens,       fock,
     $          ntri,        natonce,    denspar,    ncs,
     $          map_fs)
            endif
          endif
        else
c
c          on-disk storage :
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
           nblcount=0
           rescale=.false.
 9876      continue
c
           call secund(r1)
           call elapsec(e1)
c
           call read_33(160,iblstore,iqstore,nsplit) 
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
           call getmem(integbuf8,ixint8)
           call getint_4(integbuf4,ixint4)
           call getint_2(integbuf2,ixint2)
           call getint_2(integbu62,ixin62)
           call getint_2(integbu72,ixin72)
c
c----------------Fock-builders with on-disk integrals------------------
          If(nsym.eq.0) Then
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
                call  nd1g0_disk_nos(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   dens,       fock,
     $        ntri,        natonce,    idft,       ax,
     $        ncs,         denspar,    map_fs)
            else
                call  nd1g0_disk_nosC(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   dens,       fock,
     $        ntri,        natonce,    ncs,        denspar,
     $        map_fs)
            endif
          else
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
                call  nd1g0_disk_sym(
     $        nsym,        bl(ifp),    ncf,        bl(ixin72),
     $        bl(ixin62),  bl(ixint2), bl(ixint4), bl(ixint8),
     $        bl(ijkllab), bl(ijklsiz),bl(nquarts),iqstore,
     $        iblstore,    rescale,    nsplit,     bl(ifromto),
     $        bl(nofinteg),thres0,     thres1,     bl(lind),
     $        dens,        fock,       ntri,       natonce,
     $        idft,        ax,         ncs,        denspar,
     $        map_fs)
            else
                call  nd1g0_disk_symC(
     $        nsym,        bl(ifp),    ncf,        bl(ixin72),
     $        bl(ixin62),  bl(ixint2), bl(ixint4), bl(ixint8),
     $        bl(ijkllab), bl(ijklsiz),bl(nquarts),iqstore,
     $        iblstore,    rescale,    nsplit,     bl(ifromto),
     $        bl(nofinteg),thres0,     thres1,     bl(lind),
     $        dens,        fock,       ntri,       natonce,
     $        ncs,         denspar,    map_fs)
            endif
          endif
          call retmem(10)
          if(.not.rescale) go to 9876
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
c
c  check if we actually have to compute any integrals,
c  because it might happen that a slave executes this routine
c  even with mywork=-1, in order to retrieve the stored integrals
c  (this is to avoid a race condition MM 08/29/2006 )
c
      if(mywork.lt.0)return
c----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c re-calculate integrals
c
      where='fock'
c
c set integral timings to zero as well as neglect stuff
c
      call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
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
        call calcint2(isupbl,bl,inx,thres1,denspar,where,
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
        IF(nintez.gt.0) THEN
          call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
          If(nsym.eq.0) Then
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nfock_d1g0_nos(
     $        bl(ibuffz),  idft,        ax,          dens,
     $        fock,        bl(lind),    ntri,        natend,
     $        natonce,     denspar,     map_fs,      ncs,
     $        nblsiz,      ngctoz,      nintez,      thres1,
     $        labels(lab1),labels(lab2),labels(lab3))
            Else
c -- "Pure" DFT - no exchange contribution
              call nfock_d1g0_nosC(bl(ibuffz),dens,fock,
     *                             bl(lind),ntri,natend,natonce,
     *                             denspar,map_fs,ncs,
     *                             nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            EndIf
          Else
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nfock_d1g0_sym(nsym,bl(ifp),bl(ibuffz),
     *                            idft,ax,dens,fock,
     *                            bl(lind),ntri,natend,natonce,
     *                            denspar,map_fs,ncs,
     *                            nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            Else
c -- "Pure" DFT - no exchange contribution
              call nfock_d1g0_symC(nsym,bl(ifp),bl(ibuffz),
     *                             dens,fock,
     *                             bl(lind),ntri,natend,natonce,
     *                             denspar,map_fs,ncs,
     *                             nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            EndIf
          EndIf
        ENDIF
c
c----------------------------------------------------------------
c
        if(moreint) go to 11
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
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
      subroutine select_dnat(natend,natonce,dens,ntri, select)
      implicit real*8 (a-h,o-z)
cccc  dimension dens(natonce,ntri,3),select(ntri)
      dimension dens(3,natonce,ntri),select(ntri)
      dimension natend(natonce)
c
      do i=1,ntri
         dmax=0.d0
         do iat=1,natonce
            if(natend(iat).eq.0) then   ! not converaged yet
               x=abs( dens(1,iat,i) )
               y=abs( dens(2,iat,i) )
               z=abs( dens(3,iat,i) )
               dmax=max(x,y,z,dmax)
            endif
         enddo
         select(i)=dmax
      enddo
c
      end
c=====================================================================
      subroutine nfock_d1g0_nos(
     $        buf,         idft,        ax,          dens,
     $        fock,        lind,        ntri,        natend,
     $        nat,         densp,       map_fs,      ncs,
     $        nbls,        ngcd,        lnijkl,      thres1,
     $        labels,      length,      lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
ccccc dimension dens(nat,ntri,3),fock(nat,ntri,3)
ccccc dimension dens(3,nat,ntri),fock(3,nat,ntri)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
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
c screening density
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
           dik=densp(ics,kcs)
           dil=densp(ics,lcs)
           djk=densp(jcs,kcs)
           djl=densp(jcs,lcs)
c
           dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
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
c         output : ijf & density elements xij,yij,zij
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
c..........................................................
c            write(6,1234) icf,jcf,kcf,lcf,xint0*thres1
c1234 format(4i3,1x,2(f12.8,1x))
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
c ....................Coulomb  part ....................
                 call daxpy(nat3,xint4,dens(1,klf),1,fock(1,ijf),1)
                 call daxpy(nat3,xint4,dens(1,ijf),1,fock(1,klf),1)
c ....................Exchange part ....................
                 call daxpy(nat3,-xint0,dens(1,jkf),1,fock(1,ilf),1)
                 call daxpy(nat3,-xint0,dens(1,ikf),1,fock(1,jlf),1)
                 call daxpy(nat3,-xint0,dens(1,jlf),1,fock(1,ikf),1)
                 call daxpy(nat3,-xint0,dens(1,ilf),1,fock(1,jkf),1)
c..........................................................
c
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine nfock_d1g0_sym(
     $        nsym,        ifp,         buf,         idft,
     $        ax,          dens,        fock,        lind,
     $        ntri,        natend,      nat,         densp,
     $        map_fs,      ncs,         nbls,        ngcd,
     $        lnijkl,      thres1,      labels,      length,
     $        lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c SYMMETRY
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      dimension buf(nbls,lnijkl,ngcd)
cccc  dimension dens(nat,ntri,3),fock(nat,ntri,3)
cccc  dimension dens(3,nat,ntri),fock(3,nat,ntri)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
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
c screening density
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
           dik=densp(ics,kcs)
           dil=densp(ics,lcs)
           djk=densp(jcs,kcs)
           djl=densp(jcs,lcs)
c
           dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
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
c
          IF(dmax*abs(xint0).gt.0.5d0) THEN
c
            xint4=xint0*4.d0
            If(idft.gt.0) xint0=xint0*ax
c-------------------------------------------------------
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
c-------------------------------------------------------
c
c ....................Coulomb  part ....................
                 call daxpy(nat3,xint4,dens(1,klf),1,fock(1,ijf),1)
                 call daxpy(nat3,xint4,dens(1,ijf),1,fock(1,klf),1)
c ....................Exchange part ....................
                 call daxpy(nat3,-xint0,dens(1,jkf),1,fock(1,ilf),1)
                 call daxpy(nat3,-xint0,dens(1,ikf),1,fock(1,jlf),1)
                 call daxpy(nat3,-xint0,dens(1,jlf),1,fock(1,ikf),1)
                 call daxpy(nat3,-xint0,dens(1,ilf),1,fock(1,jkf),1)
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
c            write(6,*)' icf < 0 =',icf
             icf=-icf
             xcoul=-xcoul
             xexch=-xexch
          endif
          if(jcf.lt.0) then
c            write(6,*)' jcf < 0 =',jcf
             jcf=-jcf
             xcoul=-xcoul
             xexch=-xexch
          endif
          if(kcf.lt.0) then
c            write(6,*)' kcf < 0 =',kcf
             kcf=-kcf
             xcoul=-xcoul
             xexch=-xexch
          endif
          if(lcf.lt.0) then
c            write(6,*)' lcf < 0 =',lcf
             lcf=-lcf
             xcoul=-xcoul
             xexch=-xexch
          endif
c
          ii=lind(icf)
          jj=lind(jcf)
          kk=lind(kcf)
          ll=lind(lcf)
c
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
c
c ....................Coulomb  part ....................
                 call daxpy(nat3,xcoul,dens(1,klf),1,fock(1,ijf),1)
                 call daxpy(nat3,xcoul,dens(1,ijf),1,fock(1,klf),1)
c
c ....................Exchange part ....................
                 call daxpy(nat3,-xexch,dens(1,jkf),1,fock(1,ilf),1)
                 call daxpy(nat3,-xexch,dens(1,ikf),1,fock(1,jlf),1)
                 call daxpy(nat3,-xexch,dens(1,jlf),1,fock(1,ikf),1)
                 call daxpy(nat3,-xexch,dens(1,ilf),1,fock(1,jkf),1)
c..............................................................
       enddo   !     ns=1,nsym
       ENDIF
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine nfock_d1g0_nosC(buf,dens,fock,
     *                           lind,ntri,natend,nat,densp,
     *                           map_fs,ncs,nbls,ngcd,lnijkl,
     *                           thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c NO Symmetry   ** COULOMB ONLY **
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
ccccc dimension dens(nat,ntri,3),fock(nat,ntri,3)
ccccc dimension dens(3,nat,ntri),fock(3,nat,ntri)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
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
c screening density
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
c
           dmax=4.0d0*max(dij,dkl)
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
c         output : ijf & density elements xij,yij,zij
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
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fock(icoordiat,ijf)=fock(icoordiat,ijf) +
     $                           dens(icoordiat,klf)*xint4
       enddo
       do icoordiat=1,nat3
             fock(icoordiat,klf)=fock(icoordiat,klf) +
     $                           dens(icoordiat,ijf)*xint4
       enddo
c
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine nfock_d1g0_symC(nsym,ifp,buf,dens,fock,
     *                           lind,ntri,natend,nat,densp,
     *                           map_fs,ncs,nbls,ngcd,lnijkl,
     *                           thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c SYMMETRY   ** COULOMB ONLY **
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      dimension buf(nbls,lnijkl,ngcd)
cccc  dimension dens(nat,ntri,3),fock(nat,ntri,3)
cccc  dimension dens(3,nat,ntri),fock(3,nat,ntri)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
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
c screening density
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
c
           dmax=4.0d0*max(dij,dkl)
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
c
          if(dmax*abs(xint0).gt.0.5d0) then
c
             xint4=xint0*4.d0
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
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fock(icoordiat,ijf)=fock(icoordiat,ijf) +
     $                           dens(icoordiat,klf)*xint4
       enddo
       do icoordiat=1,nat3
             fock(icoordiat,klf)=fock(icoordiat,klf) +
     $                           dens(icoordiat,ijf)*xint4
       enddo
c..............................................................
c
       do ns=1,nsym
c
          xcoul=xint4
c
          icf=ifp(ns,ic1)
          jcf=ifp(ns,jc1)
          kcf=ifp(ns,kc1)
          lcf=ifp(ns,lc1)
c
          if(icf.lt.0) then
c            write(6,*)' icf < 0 =',icf
             icf=-icf
             xcoul=-xcoul
          endif
          if(jcf.lt.0) then
c            write(6,*)' jcf < 0 =',jcf
             jcf=-jcf
             xcoul=-xcoul
          endif
          if(kcf.lt.0) then
c            write(6,*)' kcf < 0 =',kcf
             kcf=-kcf
             xcoul=-xcoul
          endif
          if(lcf.lt.0) then
c            write(6,*)' lcf < 0 =',lcf
             lcf=-lcf
             xcoul=-xcoul
          endif
c
          ii=lind(icf)
          jj=lind(jcf)
          kk=lind(kcf)
          ll=lind(lcf)
c
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
c
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fock(icoordiat,ijf)=fock(icoordiat,ijf) +
     $                           dens(icoordiat,klf)*xcoul
       enddo
       do icoordiat=1,nat3
             fock(icoordiat,klf)=fock(icoordiat,klf) +
     $                           dens(icoordiat,ijf)*xcoul
       enddo
c..............................................................
       enddo   !     ns=1,nsym
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c==============================================================
      subroutine nd1g0_core_nos(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        dens,        fock,       ntri,       nat,
     $        idft,        ax,         densp,      ncs,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
              ij=lind(max0(icf,jcf))+min0(icf,jcf)
              do kcf=kcff+1,kcff+klen
                ik=lind(max0(icf,kcf))+min0(icf,kcf)
                jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
                do lcf=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint4=xint0*4.d0
                    if(idft.gt.0) xint0=xint0*ax

                    il=lind(max0(icf,lcf))+min0(icf,lcf)
                    jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                    kl=lind(max0(kcf,lcf))+min0(kcf,lcf)

c ......................Coulomb  part ....................
                    call daxpy(nat3,xint4,dens(1,kl),1,fock(1,ij),1)
                    call daxpy(nat3,xint4,dens(1,ij),1,fock(1,kl),1)
c ......................Exchange part ....................
                    call daxpy(nat3,-xint0,dens(1,jl),1,fock(1,ik),1)
                    call daxpy(nat3,-xint0,dens(1,ik),1,fock(1,jl),1)
                    call daxpy(nat3,-xint0,dens(1,il),1,fock(1,jk),1)
                    call daxpy(nat3,-xint0,dens(1,jk),1,fock(1,il),1)
c..........................................................
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c==============================================================
      subroutine nd1g0_core_sym(
     $        nsym,        ifp,        core,        ilab,
     $        jlab,       klab,        llab,        isiz,
     $        jsiz,       ksiz,        lsiz,        iqrt,
     $        integrals,  nqstore,     nblstore,    thres0,
     $        thres1,     lind,        dens,        fock,
     $        ntri,       nat,         idft,        ax,
     $        densp,      ncs,         map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do ic1=icff+1,icff+ilen
            do jc1=jcff+1,jcff+jlen
              ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
              do kc1=kcff+1,kcff+klen
                ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
                jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
                do lc1=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint4=xint0*4.d0
                    if(idft.gt.0) xint0=xint0*ax

                    il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                    jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                    kl=lind(max0(kc1,lc1))+min0(kc1,lc1)

c ......................Coulomb  part ....................
                    call daxpy(nat3,xint4,dens(1,kl),1,fock(1,ij),1)
                    call daxpy(nat3,xint4,dens(1,ij),1,fock(1,kl),1)
c ......................Exchange part ....................
                    call daxpy(nat3,-xint0,dens(1,jl),1,fock(1,ik),1)
                    call daxpy(nat3,-xint0,dens(1,ik),1,fock(1,jl),1)
                    call daxpy(nat3,-xint0,dens(1,il),1,fock(1,jk),1)
                    call daxpy(nat3,-xint0,dens(1,jk),1,fock(1,il),1)
c..........................................................
                    do ns=1,nsym
                      xco=xint4
                      xex=xint0
                      icf=ifp(ns,ic1)
                      jcf=ifp(ns,jc1)
                      kcf=ifp(ns,kc1)
                      lcf=ifp(ns,lc1)
                      if(icf.lt.0) then
                        icf=-icf
                        xco=-xco
                        xex=-xex
                      endif
                      if(jcf.lt.0) then
                        jcf=-jcf
                        xco=-xco
                        xex=-xex
                      endif
                      if(kcf.lt.0) then
                        kcf=-kcf
                        xco=-xco
                        xex=-xex
                      endif
                      if(lcf.lt.0) then
                        lcf=-lcf
                        xco=-xco
                        xex=-xex
                      endif
                      lic=lind(icf)
                      ljc=lind(jcf)
                      lkc=lind(kcf)
                      llc=lind(lcf)
                      ijc=max0(lic,ljc)+min0(icf,jcf)
                      ikc=max0(lic,lkc)+min0(icf,kcf)
                      ilc=max0(lic,llc)+min0(icf,lcf)
                      jkc=max0(ljc,lkc)+min0(jcf,kcf)
                      jlc=max0(ljc,llc)+min0(jcf,lcf)
                      klc=max0(lkc,llc)+min0(kcf,lcf)
c ......................Coulomb  part ....................
                      call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                      call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c ......................Exchange part ....................
                      call daxpy(nat3,-xex,dens(1,jlc),1,fock(1,ikc),1)
                      call daxpy(nat3,-xex,dens(1,ikc),1,fock(1,jlc),1)
                      call daxpy(nat3,-xex,dens(1,ilc),1,fock(1,jkc),1)
                      call daxpy(nat3,-xex,dens(1,jkc),1,fock(1,ilc),1)
c..........................................................
                    enddo
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c==============================================================
      subroutine nd1g0_core_nosC(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        dens,        fock,       ntri,       nat,
     $        densp,      ncs,         map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, NO Symmetry, Coulomb term only
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
              ij=lind(max0(icf,jcf))+min0(icf,jcf)
              do kcf=kcff+1,kcff+klen
                do lcf=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint4=xint0*4.d0

                    kl=lind(max0(kcf,lcf))+min0(kcf,lcf)

c ......................Coulomb  part ....................
                    call daxpy(nat3,xint4,dens(1,kl),1,fock(1,ij),1)
                    call daxpy(nat3,xint4,dens(1,ij),1,fock(1,kl),1)
c..........................................................
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c==============================================================
      subroutine nd1g0_core_symC(
     $        nsym,        ifp,        core,        ilab,
     $        jlab,       klab,        llab,        isiz,
     $        jsiz,       ksiz,        lsiz,        iqrt,
     $        integrals,  nqstore,     nblstore,    thres0,
     $        thres1,     lind,        dens,        fock,
     $        ntri,       nat,         densp,       ncs,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, Symmetry, Coulomb term only
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do ic1=icff+1,icff+ilen
            do jc1=jcff+1,jcff+jlen
              ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
              do kc1=kcff+1,kcff+klen
                do lc1=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint4=xint0*4.d0

                    kl=lind(max0(kc1,lc1))+min0(kc1,lc1)

c ......................Coulomb  part ....................
                    call daxpy(nat3,xint4,dens(1,kl),1,fock(1,ij),1)
                    call daxpy(nat3,xint4,dens(1,ij),1,fock(1,kl),1)
c..........................................................
                    do ns=1,nsym
                      xco=xint4
                      icf=ifp(ns,ic1)
                      jcf=ifp(ns,jc1)
                      kcf=ifp(ns,kc1)
                      lcf=ifp(ns,lc1)
                      if(icf.lt.0) then
                        icf=-icf
                        xco=-xco
                      endif
                      if(jcf.lt.0) then
                        jcf=-jcf
                        xco=-xco
                      endif
                      if(kcf.lt.0) then
                        kcf=-kcf
                        xco=-xco
                      endif
                      if(lcf.lt.0) then
                        lcf=-lcf
                        xco=-xco
                      endif
                      ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                      klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ......................Coulomb  part ....................
                      call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                      call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c..........................................................
                    enddo
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_nos(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dens,       fock,
     $        ntri,        nat,        idft,       ax,
     $        ncs,         densp,      map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry
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
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
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
      nat3=nat*3
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
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
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
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_nos(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         idft,        ax,         dmax,      intx8,
     $         xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_nos(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         idft,        ax,         dmax,      intx4,
     $         xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_nos(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         idft,        ax,         dmax,      intx2,
     $         xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_nos(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dens,        fock,       lind,      nat3,
     $             idft,        ax,         dmax,      int62,
     $             xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_nos(
     $              ilen,        jlen,       klen,      llen,
     $              icff,        jcff,       kcff,      lcff,
     $              dens,        fock,       lind,      nat3,
     $              idft,        ax,         dmax,      int72,
     $              xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_nosC(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dens,       fock,
     $        ntri,        nat,        ncs,        densp,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry, Coulomb term only
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
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
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
      nat3=nat*3
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
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
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
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_nosC(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         dmax,        intx8,      xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_nosC(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         dmax,        intx4,      xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_nosC(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         dmax,        intx2,      xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_nosC(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dens,        fock,       lind,      nat3,
     $             dmax,        int62,      xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_nosC(
     $              ilen,        jlen,       klen,      llen,
     $              icff,        jcff,       kcff,      lcff,
     $              dens,        fock,       lind,      nat3,
     $              dmax,        int72,      xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_sym(
     $        nsym,        ifp,        ncf,         xinte72,
     $        xinte62,     xinteg2,    xinteg4,     xinteg8,
     $        ijkllab,     ijklsiz,    nquarts,     iqstore,
     $        iblstore,    rescale,    nsplit,      ifromto,
     $        nofinteg,    thres0,     thres1,      lind,
     $        dens,        fock,       ntri,        nat,
     $        idft,        ax,         ncs,         densp,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry
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
      dimension ifp(7,*)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
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
      nat3=nat*3
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
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
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
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_sym(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         idft,        ax,         dmax,      nsym,
     $         ifp,         intx8,      xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_sym(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         idft,        ax,         dmax,      nsym,
     $         ifp,         intx4,      xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_sym(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         idft,        ax,         dmax,      nsym,
     $         ifp,         intx2,      xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_sym(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dens,        fock,       lind,      nat3,
     $             idft,        ax,         dmax,      nsym,
     $             ifp,         int62,      xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_sym(
     $              ilen,        jlen,       klen,      llen,
     $              icff,        jcff,       kcff,      lcff,
     $              dens,        fock,       lind,      nat3,
     $              idft,        ax,         dmax,      nsym,
     $              ifp,         int72,      xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_symC(
     $        nsym,        ifp,        ncf,         xinte72,
     $        xinte62,     xinteg2,    xinteg4,     xinteg8,
     $        ijkllab,     ijklsiz,    nquarts,     iqstore,
     $        iblstore,    rescale,    nsplit,      ifromto,
     $        nofinteg,    thres0,     thres1,      lind,
     $        dens,        fock,       ntri,        nat,
     $        ncs,         densp,      map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry
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
      dimension ifp(7,*)
      dimension dens(3*nat,ntri),fock(3*nat,ntri)
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
      nat3=nat*3
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
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
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
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_symC(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         dmax,        nsym,       ifp,       intx8,
     $         xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_symC(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         dmax,        nsym,       ifp,       intx4,
     $         xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_symC(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dens,        fock,       lind,      nat3,
     $         dmax,        nsym,       ifp,       intx2,
     $         xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_symC(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dens,        fock,       lind,      nat3,
     $             dmax,        nsym,       ifp,       int62,
     $             xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_symC(
     $              ilen,        jlen,       klen,      llen,
     $              icff,        jcff,       kcff,      lcff,
     $              dens,        fock,       lind,      nat3,
     $              dmax,        nsym,       ifp,       int72,
     $              xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,fock,1)
      endif
c
      end
c====================================================================
      subroutine make_nfock8_nos(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        idft,        ax,         dmax,      intx8,
     $        xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            ik=lind(max0(icf,kcf))+min0(icf,kcf)
            jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
            do lcf=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(icf,lcf))+min0(icf,lcf)
                jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dens(1,jl),1,fock(1,ik),1)
                call daxpy(nat3,-xint,dens(1,ik),1,fock(1,jl),1)
                call daxpy(nat3,-xint,dens(1,il),1,fock(1,jk),1)
                call daxpy(nat3,-xint,dens(1,jk),1,fock(1,il),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock8_sym(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        idft,        ax,         dmax,      nsym,
     $        ifp,         intx8,      xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
            jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
            do lc1=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dens(1,jl),1,fock(1,ik),1)
                call daxpy(nat3,-xint,dens(1,ik),1,fock(1,jl),1)
                call daxpy(nat3,-xint,dens(1,il),1,fock(1,jk),1)
                call daxpy(nat3,-xint,dens(1,jk),1,fock(1,il),1)
c......................................................
                do ns=1,nsym
                  xco=xin4
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  lic=lind(icf)
                  ljc=lind(jcf)
                  lkc=lind(kcf)
                  llc=lind(lcf)
                  ijc=max0(lic,ljc)+min0(icf,jcf)
                  ikc=max0(lic,lkc)+min0(icf,kcf)
                  ilc=max0(lic,llc)+min0(icf,lcf)
                  jkc=max0(ljc,lkc)+min0(jcf,kcf)
                  jlc=max0(ljc,llc)+min0(jcf,lcf)
                  klc=max0(lkc,llc)+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                  call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c ..................Exchange part ....................
                  call daxpy(nat3,-xex,dens(1,jlc),1,fock(1,ikc),1)
                  call daxpy(nat3,-xex,dens(1,ikc),1,fock(1,jlc),1)
                  call daxpy(nat3,-xex,dens(1,ilc),1,fock(1,jkc),1)
                  call daxpy(nat3,-xex,dens(1,jkc),1,fock(1,ilc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock8_nosC(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        dmax,        intx8,      xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_nfock8_symC(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        dmax,      nsym,         ifp,       intx8,
     $        xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            do lc1=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c......................................................
                do ns=1,nsym
                  xco=xin4
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                  klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                  call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock4_nos(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        idft,        ax,         dmax,      intx4,
     $        xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            ik=lind(max0(icf,kcf))+min0(icf,kcf)
            jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
            do lcf=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(icf,lcf))+min0(icf,lcf)
                jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dens(1,jl),1,fock(1,ik),1)
                call daxpy(nat3,-xint,dens(1,ik),1,fock(1,jl),1)
                call daxpy(nat3,-xint,dens(1,il),1,fock(1,jk),1)
                call daxpy(nat3,-xint,dens(1,jk),1,fock(1,il),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock4_sym(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        idft,        ax,         dmax,      nsym,
     $        ifp,         intx4,      xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
            jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
            do lc1=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dens(1,jl),1,fock(1,ik),1)
                call daxpy(nat3,-xint,dens(1,ik),1,fock(1,jl),1)
                call daxpy(nat3,-xint,dens(1,il),1,fock(1,jk),1)
                call daxpy(nat3,-xint,dens(1,jk),1,fock(1,il),1)
c......................................................
                do ns=1,nsym
                  xco=xin4
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  lic=lind(icf)
                  ljc=lind(jcf)
                  lkc=lind(kcf)
                  llc=lind(lcf)
                  ijc=max0(lic,ljc)+min0(icf,jcf)
                  ikc=max0(lic,lkc)+min0(icf,kcf)
                  ilc=max0(lic,llc)+min0(icf,lcf)
                  jkc=max0(ljc,lkc)+min0(jcf,kcf)
                  jlc=max0(ljc,llc)+min0(jcf,lcf)
                  klc=max0(lkc,llc)+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                  call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c ..................Exchange part ....................
                  call daxpy(nat3,-xex,dens(1,jlc),1,fock(1,ikc),1)
                  call daxpy(nat3,-xex,dens(1,ikc),1,fock(1,jlc),1)
                  call daxpy(nat3,-xex,dens(1,ilc),1,fock(1,jkc),1)
                  call daxpy(nat3,-xex,dens(1,jkc),1,fock(1,ilc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock4_nosC(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        dmax,        intx4,      xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_nfock4_symC(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        dmax,        nsym,       ifp,       intx4,
     $        xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            do lc1=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c......................................................
                do ns=1,nsym
                  xco=xin4
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                  klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                  call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock2_nos(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        idft,        ax,         dmax,      intx2,
     $        xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            ik=lind(max0(icf,kcf))+min0(icf,kcf)
            jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
            do lcf=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(icf,lcf))+min0(icf,lcf)
                jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dens(1,jl),1,fock(1,ik),1)
                call daxpy(nat3,-xint,dens(1,ik),1,fock(1,jl),1)
                call daxpy(nat3,-xint,dens(1,il),1,fock(1,jk),1)
                call daxpy(nat3,-xint,dens(1,jk),1,fock(1,il),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock2_sym(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        idft,        ax,         dmax,      nsym,
     $        ifp,         intx2,      xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
            jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
            do lc1=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dens(1,jl),1,fock(1,ik),1)
                call daxpy(nat3,-xint,dens(1,ik),1,fock(1,jl),1)
                call daxpy(nat3,-xint,dens(1,il),1,fock(1,jk),1)
                call daxpy(nat3,-xint,dens(1,jk),1,fock(1,il),1)
c......................................................
                do ns=1,nsym
                  xco=xin4
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  lic=lind(icf)
                  ljc=lind(jcf)
                  lkc=lind(kcf)
                  llc=lind(lcf)
                  ijc=max0(lic,ljc)+min0(icf,jcf)
                  ikc=max0(lic,lkc)+min0(icf,kcf)
                  ilc=max0(lic,llc)+min0(icf,lcf)
                  jkc=max0(ljc,lkc)+min0(jcf,kcf)
                  jlc=max0(ljc,llc)+min0(jcf,lcf)
                  klc=max0(lkc,llc)+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                  call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c ..................Exchange part ....................
                  call daxpy(nat3,-xex,dens(1,jlc),1,fock(1,ikc),1)
                  call daxpy(nat3,-xex,dens(1,ikc),1,fock(1,jlc),1)
                  call daxpy(nat3,-xex,dens(1,ilc),1,fock(1,jkc),1)
                  call daxpy(nat3,-xex,dens(1,jkc),1,fock(1,ilc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock2_nosC(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        dmax,        intx2,      xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_nfock2_symC(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dens,        fock,       lind,      nat3,
     $        dmax,        nsym,       ifp,       intx2,
     $        xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(nat3,*),fock(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            do lc1=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin4=4.d0*xint
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                call daxpy(nat3,xin4,dens(1,kl),1,fock(1,ij),1)
                call daxpy(nat3,xin4,dens(1,ij),1,fock(1,kl),1)
c......................................................
                do ns=1,nsym
                  xco=xin4
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                  klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  call daxpy(nat3,xco,dens(1,klc),1,fock(1,ijc),1)
                  call daxpy(nat3,xco,dens(1,ijc),1,fock(1,klc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
