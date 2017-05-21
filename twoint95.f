      subroutine twoint95(lcorex,nsymx, maxpricex,ncachex,iroutex,
     *                    icheckx,iprintx, threshx,ioutputx,iforwhatx,
     *                    lsemi,nblocks,maxbuffer,maxlabels,
     *                    scftypex,xintxx)
c-------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
      logical firstd
      character*11 scftype
      character*11 scftypex
      character*4 where
c-------------------------------------------------
c INPUT PARAMETERS (16):
c
c 1.  lcorex    - size of the bl array (common /big/ )
c 2.  nsymx     - number of symmetry operations
c 3.  maxpricex - maximum price a user is willing to pay
c                 for re-calculations of integrals
c                 (higher value more recalculated integrals
c                  i.e. more direct run )
c
c Next five following parameters are responsible for blocking
c structure :
c
c 4.  ncachex   - a factor used to increase cache memory
c                 in order to allow bigger blocks
c 5.  iroutex   - shows which of two different blocking strategy
c                 to use : (a) block contains shells shering
c                              1.type, 2.contraction, 3. gen.cont.
c                          (b) three above plus
c                              4.atom
c 6.  icheckx   - check-run if >0
c 7.  iprintx   - print level
c 8.  threshx   - integral threshold
c
c 9.  ioutputx  - file to write info
c 10. iforwhatx - shows the type of integrlas to be calculated;
c                 1) for ordinary two-el.integrals (iforwhat=1)
c                 2) for GIAO two-el. derivatives  (iforwhat=2)
c                 3) for gradient derivatives      (iforwhat=3)
c                 4) for second derivatives        (iforwhat=4)
c                 5) for MP2 integrals             (iforwhat=5)
c 11. lsemi     - integer flag to use semi-direct if available
c                  0 - do not use; otherwise use
c ......................................................................
c
c OUTPUT PARAMETERS :
c
c nblocks    - total number of blocks of contracted shell quartets
c              (super-blocks)
c maxbuffer  - maximum size for integral buffer , calc in blksizer
c maxlabels  - maximum size for labels info     , calc in blksizer
c
c scftypex   - scf mode to be run
c xintxx(9)  contains number of integrals more expensive than :
c    0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75%, 99% of maximum price
c-------------------------------------------------
c In addition the following commons are used for
c an external communications :
c
      common /cpu/ intsize,iacc,icache,memreal
c     common /big/ bl(1)
c     common /intbl/ maxsh,ifp,inx(1)
c NOTE :
c different /intbl/ maxsh,inx(1) has been used in the newest TX95 !!!!
c-------------------------------------------------
c computer dependent variables:
c      intsize - ratio double precision to integer
c      icache  - cache memory size
c      memreal - real memory
c These parameters are set in DRIVER
c
c     common /cpu/ intsize,iacc,icache,memreal
c-------------------------------------------------
c internal commons of the two-el. program :
c
c-------------------------------------------------
c options for integral calculations (from readin).
c This common block is used ONLY here and in chshift.f
c
      common /intgop/ ncache,maxprice,iprint,iblock
c-------------------------------------------------
c options for integral calculations passed from INTGOP
c to the two-el. integ. program . This common block is
c used in : prepint2.f, blocking.f, precalc.f, calcint2.f
c           and shift2.f
c
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
c-------------------------------------------------
      common /intlim/ limxmem,limblks,limpair
c------------------------------------------------
c this shows which blocking strategy is going to be used:
c
      common /route/ iroute
c-------------------------------------------------
      common /ener/ title(16),etot,enuc,qx(3),qtot,time,tlast,par(30)
      common /tape/ inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c new common for ioutput :
      common /outfile/ ioutput
c-------------------------------------------------
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c-------------------------------------------------
      dimension xintxx(9)
c-------------------------------------------------
c
c prepare for two-electron integrals calculations
c
c-------------------------------------------------
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('ictr',ictr)
c-------------------------------------------------
c transfer all input parameters into local commons :
c
      lcore   =lcorex
      nsym    =nsymx
      maxprice=maxpricex
c responsible for blocking structure :
      ncache  =ncachex
c  now these are defied in intcal
c      limxmem =limxmex
c      limblks =limblkx
c      limpair =limpaix
      iroute  =iroutex
c
c check-run and print level :
      icheck  =icheckx
      iprint  =iprintx
ctest
c     icheck  =1
c     iprint  =1
ctest
      nprint  =iprint
c output file :
      ioutput =ioutputx
c integral's threshold :
      eps     =threshx
      par(1)  =eps
c
c-------------------------------
c not used :logical firstd and :
c
      ndirect=0
      nbeg=0
      nend=0
c-------------------------------
c not needed but still present :
c
      iblock=0
      iblok=iblock
c-------------------------------------------------
c  ICACHE is used to determine the block-size i.e. the
c  number of of contracted shell quartets in a block .
c  Since that, it might be useful to increase this
c  parameter in some cases to have bigger blocks.
c  For this reason the "ncache" parameter is introduced.
c
      icacher=icache
      icache=icache*ncache
c
c-------------------------------------------------------
c
      scftype='    '
      iforwhat=iforwhatx
      call prepint2(bl,eps,inuc,ibas,na,nbf,nsh,ncf,ncs,bl(ictr),
     *              lcore,nsym,lsemi,scftype)
c
      scftypex=scftype
      icache=icacher
c-------------------------------------------------------
c for LMP2 (MP2) integrals change iforwhatx to 1
c
c     if(iforwhatx.eq.5) then
c        iforwhat =1
c        iforwhatx=1
c     endif
c-------------------------------------------------------
c these two are setup in blocking procedure executed by
c prepint2 :
      maxbuffer=maxibuf
      maxlabels=maxindx
      nblocks=nbl2*(nbl2+1)/2
c-------------------------------------------------------
      iyes=0
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc0)
      if(iftc0.NE.0) then
         call getival('nftcbl4',nftcbl4)   ! pointer
         call getival('nbl4ftc',nbl4ftc)   ! no of FTC blocks
         nblocks=nbl4ftc
      endif
c-------------------------------------------------------
      end
