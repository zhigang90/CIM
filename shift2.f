c=======================================================================
C July 07,97 KW : symmetrization of the Fock matrices for NMR chemical
C                 shifts has been added. The new routine is named
c                 FOCKSYMM_NMR(ngener,ncf,fock,ifp,nsymop)
c                 and is located at the end of this file.
c
C                 This routine is usd for both the giao and cphf
c                 contributions to the derivative Fock matrices
c                 F(D0,g1) (g1=giao derivatives) and
c                 F(D1,g0) (g0=ordinary two-el. integrals)
c
c                 focksymm_nmr is caled from
c                 (1) shift2    (line 176) to symmetrize F(D0,g1)
c                 (2) para_cphf (line  94) to symmetrize F(D1,g0)
c=======================================================================
c            +-------------------------------------+
c            +                                     +
c            +   NMR chemical shielding tensor     +
c            +                                     +
c            + two-electron integral derivatives   +
c            +                                     +
c            +     with the GIAO basis set         +
c            +                                     +
c            +         Krzysztof Wolinski          +
c            +           December, 1995            +
c            +-------------------------------------+
c=======================================================================
      subroutine shift2(idft,ax,bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstd
      character*11 scftype
      character*4 where
c------------------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
c
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c---
      common /neglect/ eps,eps1,epsr
c------------------------------------------------
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c------------------------------------------------
c NMR options :
c
      common /forcdbl/ thre1,thre2,tchf
      common /nmrint/ ichf,ncache,iprint,nogiao,noacce,ngauge,ntrans
      common /nmrpri/ maxprice
c---
      common /cux/ lden,lfoc,lforc,love,lh01,lhfc,lval,lvec,lrem
c------------------------------------------------
c     common /intbl/ifpp,inxx(100)
c------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension xintxx(9)
      data half,one /0.5d0,1.d0/
c----------------------------------------------------------------------
c remember the common block /neglect/ eps,eps1,epsr which has been used
c in the INTEG for SCF. It is needed because the integral threshold may
c be different for zero-order and first-order GIAO integrals. Previous
c values have to be returned for zero-order integrals in the CPHF
c procedure.
c
      reps=eps
      reps1=eps1
      repsr=epsr
      icacher=icache
c----------------------------------------------------------------------
c prepare the integral system (for both giao and cphf integrals) :
c (like twoint95 for ordinary integrals)
c     call twoint95(lcore,nsym, maxprice,
c    1              ncache,limxmem,limblks,limpair,iroute,
c    2              icheck,iprint, thresh,istat,iforwhat,
c    3              nblocks,maxbuffer,maxlabels,
c    4              scftype,xintxx)
c Full list of parameters for twoint95 (13) is reduced here to only
c FOUR . The remaning parameters are taken from commons set up for
c the SCF . Hopefully, it will be always the case for post-scf calc.
c I want to use this interface for all post-scf applications (nmr,grad
c and second derivatives). Four input parameters will be taken from
c appropriate commons set up when a given application is called. Here
c four input parameters taken from NMR option commons :
c
c
      thresh=thre2
      iforwhat=2
c
c ::::::::::::::::::::::::::::::::::::::::::
c -- parallel?
      call getival('nslv',nslv)
      If(nslv.GT.0) call para_JobInit(iforwhat)
c ::::::::::::::::::::::::::::::::::::::::::
c
      call post_scf(bl,inx,ncache,iprint,thresh,
     *              iforwhat,1,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
      if(scftype.ne.'full-direct')then
        call getival('incore',incore)
        if(incore.gt.0)then
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with core storage'
        else
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with disk storage'
        endif
      endif
c
c The previous call sets up the local integral system
c the following initializes the slaves if we are running
c in parallel. Many input parameters are taken from commons
c or not sent (see comments in ptwoint.f)
c
      call para_post_scf(ncache,iforwhat,nblocks)

c output (see scf) : last line parameters
c----------------------------------------------------------------------
c if it is not the GIAO method :
c     if(nogiao.ne.0) return
c----------------------------------------------------------------------
c reserve memory for the malkin corrections (kept during CPHF I hope)
c
      call getival('malk',malk)
      if (malk.ne.0) then
         call getival('nmo',nmo)
         call getmem(nmo*(ncf-nmo),imalk)
         call zeroit(bl(imalk),nmo*(ncf-nmo))
         last=last+nmo*(ncf-nmo)
      else
         imalk=0  ! it will not be used!!!
      endif
      call setival('imalk',imalk)
c----------------------------------------------------------------------
c Calculate the GIAO two-el.int.deriv. and contract them with density
c Put output to the  fockX,fockY,fockZ matrices
c
c GIAO integrals will NEVER be re-calculated.
c----------------------------------------------------------------------
c reserve memory for label's arrays :
c
      call getmem(maxlabels,labels)
      last=last+maxlabels
c----------------------------------------------------------------------
c if it is not the GIAO method :
c release memory reserved for labels
      if(nogiao.ne.0) then
         call retmem(1)
         last=last-maxlabels
         return
      endif
c----------------------------------------------------------------------
c allocate memory for the G(D,g1) part of the 1st-order Fock matrices X,Y,Z
c
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
      call getmem(ntri3,lfxyz)
      last=last+ntri3
c
      call zeroit(bl(lfxyz),ntri3)
c----------------------------------------------------------------------
c density is already in lden ( chshift) 
c
      call getival('lden',lden)
c----------------------------------------------------------------------
c read-in unperturbed density :
c
c     call getmem(ntri3,lden)
c
c      np4=4
c      call sudat(np4,'den0_rhf',ni)
c      if(ni.gt.0) then
c        call rea (bl(lden),ntri,np4,'den0_rhf')
c      else
c        call restart(np4,0)
c        call restart(np1,0)
c        call sudat(np4,'den0_rhf',ni)
c        if (ni.gt.0) then
c          call rea (bl(lden),ntri,np4,'den0_rhf')
c        else
c          call nerror(2,'shift1 in memres1',
c    1    'Program cannot find the 0th-oder density on <jobname>.14',
c    2     ni,ni)
c        endif
c      endif
c----------------------------------------------------------------------
c update last
      call getmem(1,last)
      call retmem(1)
c----------------------------------------------------------------------
      call secund(tgiao1)
      call elapsec(etgi1)
c----------------------------------------------------------------------
c screening density :
c transform ordinary density dens(ij) into denspar(ics,jcs)
c
      call getmem(2*ncs*ncs,idensp)
c
      call getmem(ncf,map_fs)
c
      call setup_densp2(inx,ncf,ncs,bl(lden),bl(idensp),bl(map_fs))
c
      call retmem(1)  !  map_fs
c
c read in second screening density
c
      call read1mat(69, 1 ,ncs*ncs,bl(idensp+ncs*ncs))
c
c----------------------------------------------------------------------
c calculate two-electron integral derivatives for NMR/GIAO
c (like call int_fock )
c
      call para_giao(idft,ax,nblocks,bl,inx,ntri,thresh,
     *              bl(idensp),bl(lden),bl(lfxyz),bl(labels))
c----------------------------------------------------------------------
      call retmem(1)    !  idensp
c----------------------------------------------------------------------
      call secund(tgiao2)
      call elapsec(etgi2)
      total=(tgiao2-tgiao1)/60.0d0
      elaps=(etgi2-etgi1)/60.0d0
c     write(iout,500)
      write(iout,400) total,elaps
      call f_lush(iout)
c     write(icond,400) total
c     write(iout,500)
  400 format('Master CPU time for GIAO 2e integrals = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
  500 format(58('-'))
c----------------------------------------------------------------------
c re-scale these parts of fock matrices :
c
      call dscal(ntri3,thresh,bl(lfxyz),1)
c----------------------------------------------------------------------
c calculates the GIAO exchange matrix derivatives
c and the malkin corrections
c
      if(idft.ne.0) then
         call getival('ndum',ndum)
         call dft_giao(idft,malk,lfxyz,na-ndum,ncf,ibas,iprint)
c----------------------------------------------------------------------
         call secund(tgiao3)
         call elapsec(etgi3)
         total=(tgiao3-tgiao2)/60.0d0
         elaps=(etgi3-etgi2)/60.0d0
c        write(iout,500)
         write(iout,600) total,elaps
         call f_lush(iout)
c        write(iout,500)
  600 format('Master CPU time for GIAO XC integrals = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
      endif
c----------------------------------------------------------------------
c test
c     write(8,*)'From Shift2 : f(d0,g1)  '
c     lfx=lfxyz
c     lfy=lfx+ntri
c     lfz=lfy+ntri
c     call drumx (bl(lfx),ncf,8   ,'fock   X' )
c     call drumx (bl(lfy),ncf,8   ,'fock   Y' )
c     call drumx (bl(lfz),ncf,8   ,'fock   Z' )
c----------------------------------------------------------------------
c Symetrize these parts of the Fock matrices for exac symmetry :
c
      call getival('nsym',nsym)
      if(nsym.gt.0) then
         call getival('ngener',ngener)
         call getival('nsyo',nsyo)
         call getival('SymFunPr',ifp)
         call focksymm_nmr(ngener,ncf,bl(lfxyz),bl(ifp),bl(nsyo))
      endif
c----------------------------------------------------------------------
c read-in 1e-part of Fock derivative matrices :
c
      call getmem(ntri3,lfoc)
c
      irec=1
      call read1mat(61,irec,ntri3,bl(lfoc))
c----------------------------------------------------------------------
c and add them to the Fock (lfoc)
c
      call daxpy(ntri3,one,bl(lfxyz),1,bl(lfoc),1)
c----------------------------------------------------------------------
c save fock matrices :
c
      call save1mat(61,irec,ntri3,bl(lfoc))
c----------------------------------------------------------------------
      if(Malk.ne.0) then
          if(nsym.gt.0) then    ! do symmetry scaling
             write(iout,1000)
             call VScal(nmo*(ncf-nmo),DFloat(nsym+1),bl(imalk))
          endif
          if(iprint.ge.1) then
              write(iout,*) ' Malkin correction matrix'
              call prntmat(nmo,ncf-nmo,nmo,bl(imalk))
          endif
      endif
c
ctest
c     write(8,*)'From Shift2 : Full focks'
c     call drumx (bl(lfxyz+1),ncf,8   ,'fock   X' )
c     call drumx (bl(lfxyz+ntri),ncf,8   ,'fock   Y' )
c     call drumx (bl(lfxyz+ntri*2),ncf,8   ,'fock   Z' )
c----------------------------------------------------------------------
c release memory reserved for labels, lfxyz and lfoc
c
      call retmem(3)
c
c   keep "malk" allocation
c----------------------------------------------------------------------
c update last
      call getmem(1,last)
      call retmem(1)
c----------------------------------------------------------------------
c  Return the SCF values of eps,eps1,epsr
c
      eps=reps
      eps1=reps1
      epsr=repsr
      icache=icacher
c----------------------------------------------------------------------
c
 1000 Format('**WARNING**  With degenerate Molecular Orbitals turn',/,
     $       '  off symmetry with the Malkin correction')
c
      end
c===================================================================
      subroutine drumx (d,ncf,ifi,txt)
      implicit real*8 (a-h,o-z)
      character*8 txt
      dimension d(1)
      write (ifi,20) txt
      n=ncf
      ja=1
      do 10 i=1,n
         je=ja+i-1
         write (ifi,30) i,(d(j),j=ja,je)
         ja=je+1
   10 continue
      return
c
   20 format (30x,3h***,a8,3h***)
   30 format (1x,i4,2x,10f12.5,/,(7x,10f12.5))
c
      end
c===================================================================
c NMR/GIAO two-el.int.deriv. - precalculation's routine
c
      subroutine precal2d(datnuc,iis,jjs,inx,npij,npkl,npklx,
     *           ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *           xyab,xycd, ipres,ijcent,klcent)
      implicit real*8 (a-h,o-z)
c
      dimension datnuc(5,*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
      dimension xyab(npij,3),xycd(npkl,3)
      dimension ijcent(npij),klcent(npkl)
      dimension ipres(*)
c
c***************
c
c* precalculations for the pairs IJ :
c
      ijpar=0
      do 100 ijp=nijbeg,nijend
      ijpar=ijpar+1
        ijcs=ijbl(ibl,ijp)
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        iatom=inx(2,ics)
        jatom=inx(2,jcs)
          ijcent(ijpar)=2
          if(iatom.eq.jatom) ijcent(ijpar)=1
        xa=datnuc(2,iatom)
        ya=datnuc(3,iatom)
        za=datnuc(4,iatom)
        xb=datnuc(2,jatom)
        yb=datnuc(3,jatom)
        zb=datnuc(4,jatom)
        xyab(ijpar,3)= xa*yb-ya*xb
        xyab(ijpar,2)=-xa*zb+za*xb
        xyab(ijpar,1)= ya*zb-za*yb
  100 continue
c
c
c* precalculations for the pairs KL :
c
      IF(npklx.ne.0) then
c
        klpar=0
        do 200 klp=nklbeg,nklend
        klpar=klpar+1
          klcs=ijbl(kbl,klp)
          kcs=iis(klcs)
          lcs=jjs(klcs)
          katom=inx(2,kcs)
          latom=inx(2,lcs)
            klcent(klpar)=2
            if(katom.eq.latom) klcent(klpar)=1
          xc=datnuc(2,katom)
          yc=datnuc(3,katom)
          zc=datnuc(4,katom)
          xd=datnuc(2,latom)
          yd=datnuc(3,latom)
          zd=datnuc(4,latom)
          xycd(klpar,3)= xc*yd-yc*xd
          xycd(klpar,2)=-xc*zd+zc*xd
          xycd(klpar,1)= yc*zd-zc*yd
  200   continue
c
      ELSE
c* for a diagonal case :
c
c* since this is for a diagonal block
c* pairs IJ and KL are the same
c
        do 300 ijpar=1,npij
           klcent(ijpar)=ijcent(ijpar)
           do 310 i=1,3
           xycd(ijpar,i)=xyab(ijpar,i)
  310      continue
  300   continue
c
      ENDIF
c
c eliminate contracted quartets of the type (AA|BB) (centers)
c
      ijkl=0
      do 400 ijpar=1,npij
      ijc=ijcent(ijpar)
      npklend=npkl
      if(npklx.eq.0) npklend=ijpar
      do 400 klpar=1,npklend
      klc=klcent(klpar)
      ijkl=ijkl+1
      if(ijc.eq.1 .and. klc.eq.1) ipres(ijkl)=0
  400 continue
c
      end
c===================================================================
c
c  NMR/GIAO two-el.integ. derivatives routines :
c----
c  the nmrderx routine is called from Calcint2 when WHERE='shif'
c----
      subroutine nmrderx(nbls,lnijr,lnklr,npij,npkl,ngcd,
     *                   idx1,idx2, ixab,ixcd,ixyab,ixycd)
c

      use memory

      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
c     common /big/ bl(1)
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
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
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
c dimensions for assembling :
      common /dimasse/ lqijr,lqklr,lqmxr,lij3,lkl3,l3l,lsss
c----------------------------------------------------------
c     lqijr=nfu(nqij1+1)
c     lqklr=nfu(nqkl1+1)
c     lqmxr=lqijr
c     if(lqklr.gt.lqijr) lqmxr=lqklr
c
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=lqij
      if(lqkl.gt.lqij) lqmx=lqkl
c----------------------------------------------------------
      call getmem(3*nbls,ixabq)
      call getmem(3*nbls,ixcdq)
      call getmem(3*nbls,ixyabq)
      call getmem(3*nbls,ixycdq)
c----------------------------------------------------------
      call conv24x(nbls,npij,npkl,bl(idx1),bl(idx2),
     *            bl(ixab),bl(ixcd),bl(ixyab),bl(ixycd),
     *            bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq) )
c----------------------------------------------------------
      ibeg =ibuf2
      incr1=ngcd*nbls*lnijr*lnklr
      ider=ibeg+incr1
c
      ijbex=nfu(nqij)+1
      klbex=nfu(nqkl)+1
      ijenx=lnij
      klenx=lnkl
c
      call giao_der(ngcd,nbls,bl(ibeg),lnijr,lnklr,lnij,lnkl,
     *             ijbex,ijenx, klbex,klenx,
     *             bl(ider),bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c
      ibuf2=ibuf2+incr1
c----------------------------------------------------------
      if(lshellt.eq.0) go to 100
c----------------------------------------------------------
c
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-   --- for bfij1 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lqij
         klenx=lnkl
         if(nqij.eq.nsij) ijenx=1
c
         ibeg =ibfij1
         incr1=nbls*lqijr*lnklr
         ider=ibeg+incr1
c
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lnklr,lqij,lnkl,
     *                ijbex,ijenx, klbex,klenx,
     *               bl(ider),bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfij1=ibfij1+incr1
      endif
c----------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-   --- for bfij2 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lqij
         klenx=lnkl
         if(nqij.eq.nsij) ijenx=1  ! ???????????
c
         ibeg =ibfij2
         incr1=nbls*lqijr*lnklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lnklr,lqij,lnkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfij2=ibfij2+incr1
      endif
c----------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-   --- for bfkl1 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lnij
         klenx=lqkl
         if(nqkl.eq.nskl) klenx=1
c
         ibeg =ibfkl1
         incr1=nbls*lnijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lnijr,lqklr,lnij,lqkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfkl1=ibfkl1+incr1
      endif
c----------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-   --- for bfkl2 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lnij
         klenx=lqkl
         if(nqkl.eq.nskl) klenx=1  ! ???????????
c
         ibeg =ibfkl2
         incr1=nbls*lnijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lnijr,lqklr,lnij,lqkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfkl2=ibfkl2+incr1
      endif
c----------
      if(lshellt.eq.1) go to 100
c----------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
      if(lshelij.eq.3) then
c-   --- for bfij3 (nbls,4,lnklr) ; 4 is for nmr only
c
         ijbex=1
         klbex=klbeg
         ijenx=1
         klenx=lnkl
c
         ibeg =ibfij3
         incr1=4*nbls*lnklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,lnklr,1   ,lnkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfij3=ibfij3+incr1
      endif
c----------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.3) then
c-   --- for bfkl3 (nbls,lnijr,4) ; 4 is for nmr only
c
         ijbex=ijbeg
         klbex=1
         ijenx=lnij
         klenx=1
c
         ibeg =ibfkl3
         incr1=4*nbls*lnijr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lnijr,4    ,lnij,1   ,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfkl3=ibfkl3+incr1
      endif
c----------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(1).eq.1) then
c-   --- for bf2l1 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lqij
         klenx=lqkl
         if(nqij.eq.nsij) ijenx=1
         if(nqkl.eq.nskl) klenx=1
c
         ibeg =ibf2l1
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l1=ibf2l1+incr1
      endif
c----------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(2).eq.1) then
c-   --- for bf2l2 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lqij
         klenx=lqkl
         if(nqij.eq.nsij) ijenx=1
         if(nqkl.eq.nskl) klenx=1
c
         ibeg =ibf2l2
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l2=ibf2l2+incr1
      endif
c----------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(3).eq.1) then
c-   --- for bf2l3 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lqij
         klenx=lqkl
         if(nqij.eq.nsij) ijenx=1
         if(nqkl.eq.nskl) klenx=1
c
         ibeg =ibf2l3
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l3=ibf2l3+incr1
      endif
c----------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(4).eq.1) then
c-   --- for bf2l4 ---
c
         ijbex=ijbeg
         klbex=klbeg
         ijenx=lqij
         klenx=lqkl
         if(nqij.eq.nsij) ijenx=1
         if(nqkl.eq.nskl) klenx=1
c
         ibeg =ibf2l4
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l4=ibf2l4+incr1
      endif
c----------
      if(lshellt.eq.2) go to 100
c----------
c     if(lcase.eq.12.or.lcase.eq.16) then
      if(lcas3(1).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) -first
c
cccccc   ijbex=ijbeg
         ijbex=1
         klbex=klbeg
         ijenx=1
         klenx=lqkl
c
         ibeg =ibf3l1
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,lqmxr,1   ,lqmx,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l1=ibf3l1+incr1
      endif
c----------
c     if(lcase.eq.13.or.lcase.eq.16) then
      if(lcas3(2).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) - second
c
cccccc   ijbex=ijbeg
         ijbex=1
         klbex=klbeg
         ijenx=1
         klenx=lqkl
c
         ibeg =ibf3l2
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,lqmxr,1   ,lqmx,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l2=ibf3l2+incr1
      endif
c----------
c     if(lcase.eq.14.or.lcase.eq.16) then
      if(lcas3(3).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - third
c
         ijbex=ijbeg
         klbex=1
ccccc    klbex=klbeg
         ijenx=lqij
         klenx=1
c
         ibeg =ibf3l3
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqmxr,4    ,lqmx,1   ,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l3=ibf3l3+incr1
      endif
c----------
c     if(lcase.eq.15.or.lcase.eq.16) then
      if(lcas3(4).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - fourth
c
         ijbex=ijbeg
         klbex=1
ccccc    klbex=klbeg
         ijenx=lqij
         klenx=1
c
         ibeg =ibf3l4
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqmxr,4    ,lqmx,1   ,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l4=ibf3l4+incr1
      endif
c----------
      if(lcase.eq.16) then
c-   --- for ssss(nbls)  ---
c
         ijbex=1
         klbex=1
         ijenx=1
         klenx=1
c
         ibeg =issss
         incr1=nbls*16
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,4    ,1   ,1   ,
     *                ijbex,ijenx, klbex,klenx,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         issss=issss+incr1
      endif
c
  100 continue
c...................
      call retmem(4)
      end
c=================================================================
      subroutine conv24x(nbls,npij,npkl,idx1,idx2 ,
     *                  xab ,xcd, xyab, xycd ,
     *                  xabq,xcdq,xyabq,xycdq )
      implicit real*8 (a-h,o-z)
c
      dimension idx1(nbls),idx2(nbls)
      dimension xab(npij,3) ,xcd(npkl,3) ,xyab(npij,3) ,xycd(npkl,3)
      dimension xabq(nbls,3),xcdq(nbls,3),xyabq(nbls,3),xycdq(nbls,3)
c
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
      klpar=idx2(ijkl)
        do 150 i=1,3
        xabq(ijkl,i)=xab(ijpar,i)
        xcdq(ijkl,i)=xcd(klpar,i)
        xyabq(ijkl,i)=xyab(ijpar,i)
        xycdq(ijkl,i)=xycd(klpar,i)
  150   continue
  100 continue
c
      end
c=================================================================
      subroutine giao_der(ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                   ijbex,ijenx, klbex,klenx,
     *                   deriv,  xab,xcd, xyab,xycd)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
      dimension buf2(nbls,lnijr,lnklr,ngcd)
      dimension deriv(6,nbls,lnij,lnkl,ngcd)
c--------------------------------------------
c---> dimension typ1x(nbls,lnij,lnkl,ngcd),
c--->*          typ2x(nbls,lnij,lnkl,ngcd)
c---> dimension typ1y(nbls,lnij,lnkl,ngcd),
c--->*          typ2y(nbls,lnij,lnkl,ngcd)
c---> dimension typ1z(nbls,lnij,lnkl,ngcd),
c--->*          typ2z(nbls,lnij,lnkl,ngcd)
c--------------------------------------------
      dimension xab(nbls,3),xcd(nbls,3),xyab(nbls,3),xycd(nbls,3)
c
c-----
c
      do 200 kl=klbex,klenx
      klpx=npxyz(1,kl)
      klpy=npxyz(2,kl)
      klpz=npxyz(3,kl)
      do 200 ij=ijbex,ijenx
      ijpx=npxyz(1,ij)
      ijpy=npxyz(2,ij)
      ijpz=npxyz(3,ij)
c
      do 225 iqu=1,ngcd
c
        do 250 ijkl=1,nbls
c---------------
c--x deriv.
        abzy= -xab(ijkl,3)*buf2(ijkl,ijpy,kl,iqu)
     *        +xab(ijkl,2)*buf2(ijkl,ijpz,kl,iqu)
     *       +xyab(ijkl,1)*buf2(ijkl,ij,kl,iqu)
        cdzy= -xcd(ijkl,3)*buf2(ijkl,ij,klpy,iqu)
     *        +xcd(ijkl,2)*buf2(ijkl,ij,klpz,iqu)
     *       +xycd(ijkl,1)*buf2(ijkl,ij,kl,iqu)
c
cccccc  typ1x(ijkl,ij,kl,iqu)= abzy+cdzy
cccccc  typ2x(ijkl,ij,kl,iqu)=-abzy+cdzy
        deriv(1,ijkl,ij,kl,iqu)= abzy+cdzy
        deriv(2,ijkl,ij,kl,iqu)=-abzy+cdzy
c
c
c--y deriv.
        abzx= +xab(ijkl,3)*buf2(ijkl,ijpx,kl,iqu)
     *        -xab(ijkl,1)*buf2(ijkl,ijpz,kl,iqu)
     *       +xyab(ijkl,2)*buf2(ijkl,ij,kl,iqu)
        cdzx= +xcd(ijkl,3)*buf2(ijkl,ij,klpx,iqu)
     *        -xcd(ijkl,1)*buf2(ijkl,ij,klpz,iqu)
     *       +xycd(ijkl,2)*buf2(ijkl,ij,kl,iqu)
c
cccccc  typ1y(ijkl,ij,kl,iqu)= abzx+cdzx
cccccc  typ2y(ijkl,ij,kl,iqu)=-abzx+cdzx
        deriv(3,ijkl,ij,kl,iqu)= abzx+cdzx
        deriv(4,ijkl,ij,kl,iqu)=-abzx+cdzx
c
c--z deriv.
        abyx= -xab(ijkl,2)*buf2(ijkl,ijpx,kl,iqu)
     *        +xab(ijkl,1)*buf2(ijkl,ijpy,kl,iqu)
     *       +xyab(ijkl,3)*buf2(ijkl,ij,kl,iqu)
        cdyx= -xcd(ijkl,2)*buf2(ijkl,ij,klpx,iqu)
     *        +xcd(ijkl,1)*buf2(ijkl,ij,klpy,iqu)
     *       +xycd(ijkl,3)*buf2(ijkl,ij,kl,iqu)
c
cccccc  typ1z(ijkl,ij,kl,iqu)= abyx+cdyx
cccccc  typ2z(ijkl,ij,kl,iqu)=-abyx+cdyx
        deriv(5,ijkl,ij,kl,iqu)= abyx+cdyx
        deriv(6,ijkl,ij,kl,iqu)=-abyx+cdyx
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
c*******************************************************************
c
c     symmetrization of the first-order NMR fock matrices
c
c
c===================================================================
cccc  subroutine FockSymm_SCF(ngener,ncf,nfock,fock,ifp)
      subroutine focksymm_nmr(ngener,ncf,fock,ifp,nsymop)
c---------------------------------------------------------------------
c
c    symmetrize F(D0,g1) (with the giao integral derivatives)
c                and
c    symmetrize F(D1,g0) (with the ordinary zeroth-order integrals)
c
c input : ngener= number of group generators
c         ncf   = number of basis functions
c         fock  = 3 Fock matrices (X,Y,Z perturbed)
c         ifp   = symmetry images of a given function for each sym.op.
c         nsymop= symmetry operations
c
c output: fock  = symmetrized Fock matrices
c
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension fock(*)
      dimension nsymop(7),ifp(7,ncf)
c local stuff:
      dimension mirror(3,7),negx(7),negy(7),negz(7)
      data mirror/1,-1,-1, -1,1,-1, -1,-1,1, -1,-1,1, -1,1,-1,
     *            1,-1,-1, 1,1,1/
      data half,one /0.5d0,1.d0/
c-----------------------------------------
      if (ncf.eq.1) return
c-----------------------------------------
      ntri=ncf*(ncf+1)/2
      lrix=0
      lriy=lrix+ntri
      lriz=lriy+ntri
c-----------------------------------------
      do 10 ns=1,ngener
         nop=nsymop(ns)
         negx(ns)=mirror(1,nop)
         negy(ns)=mirror(2,nop)
         negz(ns)=mirror(3,nop)
   10 continue
c-----------------------------------------
      do 350 ns=1,ngener
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if(i1.lt.0) then
               i1=-i1
               fct=-fct
            endif
            j1=ifp(ns,j)
            if (j1.lt.0) then
              j1=-j1
              fct=-fct
            endif
            ij1=i1*(i1-1)/2 +j1
            if (j1.gt.i1) then
              ij1=j1*(j1-1)/2 +i1
              fct=-fct
            endif
            ffx=fock(lrix+ij)+fct*fock(lrix+ij1)*negx(ns)
            ffy=fock(lriy+ij)+fct*fock(lriy+ij1)*negy(ns)
            ffz=fock(lriz+ij)+fct*fock(lriz+ij1)*negz(ns)
            if (ij.gt.ij1) then
               ffx=ffx*half
               ffy=ffy*half
               ffz=ffz*half
            endif
            fock(lrix+ij)=ffx
            fock(lrix+ij1)=fct*ffx*negx(ns)
            fock(lriy+ij)=ffy
            fock(lriy+ij1)=fct*ffy*negy(ns)
            fock(lriz+ij)=ffz
            fock(lriz+ij1)=fct*ffz*negz(ns)
  340    continue
  350 continue
c
      end
