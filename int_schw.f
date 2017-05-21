c=====================================================================
      subroutine int_schwarz(bl,inx,schwarz,labels,schx)
c---------------------------------------------------------------------
c This routine calls two-el. program for Schwarz integrals only (ij|ij)
c and put them into schwarz(*)
c
c INPUT:
c  bl         = common bl (free memory)
c  inx        = common inx (contraction info)
c  labels     = memory area for labels
c
c OUTPUT :
c  schwarz    = an array holding (ij|ij) integrals
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 schx
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NOT EXISTING super block
c---------------------------------------------------------------------
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      common /cpu/ intsize,iacc,icache,memreal
c???  common /counters/nintsum1,nintsum2,nquartsum
      common /memor1b/ nbl2,nbloks
      dimension inx(12,*)
      dimension bl(*),schwarz(*)
      dimension labels(*)
c----------------------------------------------------------------
c for screening in cshneg : it is not important for SCHW. int.
c but for consistency it is included here
c
      screen='fock'
c
c screening will be done like for fock-bulider
c----------------------------------------------------------------
      call zeroit(schwarz,ncs*ncs)
c----------------------------------------------------------------
c calculate Schwarz integrals
c
c     schx= sch1  for scf  Schwarz integrals
c     schx= sch2  for giao Schwarz integrals
c     schx= sch3  for forc Schwarz integrals
c     schx= sch4  for hess Schwarz integrals
c
c----------------------------------------------------------------
c setup dimension of the integral buffer :
c
      ndim=1
      if(schx.eq.'sch2') ndim=6
      if(schx.eq.'sch3') ndim=9
      if(schx.eq.'sch4') ndim=45
c----------------------------------------------------------------
c test :
cccc       write(6,*)' int_schwarz : schx=',schx,'  ncs=',ncs
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c set up denspar(ics,jcs)=1.0
c
      call getmem(ncs*ncs,idensp)
      call setup_denschw(ncs,bl(idensp))
c----------------------------------------------------------------
c make mapping array (icf-->ics) :
c
      call getint(ncf,map_cf_cs)
      call make_map_cf_cs(ncs,ncf,inx,bl(map_cf_cs))
c----------------------------------------------------------------
      thres1=1.0d-20
c----------------------------------------------------------------
c Do ONLY diagonal blocks (IJbl | IJbl)
c
      DO ijbl2=1,nbl2
        isupbl=ijbl2*(ijbl2+1)/2
c
   11 continue
c
        where=schx      ! must be set again because it
c                        is changed to 'fock' etc. in calcint
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
           lsh=labels(lab2)*labels(lab2+1)*labels(lab2+2)*labels(lab2+3)
c???       nintsum2=nintsum2+nblsiz*lsh
c???       nquartsum=nquartsum+nblsiz
c
c put Schwarz integrals into schwartz array :
c
           call schw_bldr(ncs,ncf,bl(ibuffz),ndim,nblsiz,ngctoz,nintez,
     *                    bl(map_cf_cs),
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    schwarz,thres1)
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
c re-scale schwarz integrals to the real values (*thres1)
c
      call schw_resc(schwarz,ncs,thres1)
c----------------------------------------------------------------
c release memory allocated for idensp and map_cf_cs
c
      call retmem(2)
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine setup_denschw(ncs,denspar)
      implicit real*8 (a-h,o-z)
      common /screening/ dens_max_el      ! output
      dimension denspar(ncs,ncs)          ! output
      data one /1.d0/
c--------------------------------------------------------------------
c output :
c densp(ics,jcs)    = one
c dens_max_el       = maximum element overall
c                     save in common /screening/
c--------------------------------------------------------------------
c
      dens_max_el=one
c
      do 200 ics=1,ncs
      do 200 jcs=1,ics
      denspar(ics,jcs)=one
      denspar(jcs,ics)=one
  200 continue
c
      end
c=====================================================================
      subroutine schw_bldr(ncs,ncf,buf,ndim,nbls,ngcd,lnijkl,map_cf_cs,
     *                     labels,length,lgenct, schwarz,thres1)
c----------------------------------------------------------------
c Makes an array with schwarz integrals (ics,jcs | ics,jcs) (shells)
c----------------------------------------------------------------
c  input parameters
c  map_cf_cs  = mapping from contr. functions to contracted shells
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
cccc  dimension buf(nbls,lnijkl,ngcd)
      dimension buf(ndim,nbls,lnijkl,ngcd)
      dimension schwarz(ncs,ncs)        ! output
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension map_cf_cs(ncf)
c----------------------------------------------------------------
c2002 permutation factors OFF :
      permut1=2.d0     ! because (ij|kl) with ij=kl
cccc  if(ngcd.gt.1) permut1=1.d0
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
      ninteg=ilen*jlen*klen*llen
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
              icff=labels(1, 1 ,ijklp)
              jcff=labels(2, 1 ,ijklp)
              kcff=labels(3, 1 ,ijklp)
              lcff=labels(4, 1 ,ijklp)
c
              ics=map_cf_cs(icff+1)
              jcs=map_cf_cs(jcff+1)
              kcs=map_cf_cs(kcff+1)
              lcs=map_cf_cs(lcff+1)
              permut=permut1
              if(ics.eq.jcs) permut=permut1*4.d0
c-------------------------------------------------------
c             if(ics.ne.kcs) stop 'ics.ne.kcs'
c             if(jcs.ne.lcs) stop 'jcs.ne.lcs'
c-------------------------------------------------------
           xint_max=0.d0
           do 150 iqu=1,ngcd
              do 200 intct=1,ninteg
                 do 250 ider=1,ndim
                    xint=buf(ider,ijklp,intct,iqu)
                    xint=abs(xint)
                    xint_max=max(xint_max,xint)
c------------------------------------------------------
  250            continue
  200         continue
  150      continue
c2002 permutation factors OFF :
           xint_max=xint_max*permut
           schwarz(ics,jcs)=xint_max
           schwarz(jcs,ics)=xint_max
  100   continue
c
      end
c==============================================================
      subroutine make_map_cf_cs(ncs,ncf,inx,map_cf_cs)
      dimension map_cf_cs(ncf)
      dimension inx(12,ncs)
c
      do ics=1,ncs
         icf_b=inx(11,ics)+1
         icf_e=inx(10,ics)
         do icf=icf_b,icf_e
            map_cf_cs(icf)=ics
         enddo
      enddo
c
      end
c==============================================================
      subroutine schw_resc(schwarz,ncs,thres1)
      implicit real*8 (a-h,o-z)
      dimension schwarz(ncs,ncs)
c
c     write(6,*)' Schwarz integrals'
c     write(6,*)' thres1=',thres1
      do ics=1,ncs
         do jcs=1,ncs
            schwarz(ics,jcs)=schwarz(ics,jcs)*thres1
c           write(6,66) ics,jcs,schwarz(ics,jcs)
         enddo
      enddo
  66  format(2i3,f15.9)
c
      end
c==============================================================
