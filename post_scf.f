      subroutine post_scf(bl,inx,ncachex,iprintx, threshx,
     *                    iforwhatx,lsemi,nblocks,maxbuffer,maxlabels,
     *                    scftypex,xintxx)
c-------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical firstd
      character*11 scftype
      character*11 scftypex
      character*4 where
c-------------------------------------------------
c INPUT PARAMETERS (4 13):
c
c
c Next five following parameters are responsible for blocking
c structure :
c
c 4.  ncachex   - a factor used to increase cache memory
c                 in order to allow bigger blocks
c 10. iprintx   - print level
c 11. threshx   - integral threshold
c
c 13. iforwhatx - shows the type of integrlas to be calculated;
c                 1) for ordinary two-el.integrals (iforwhat=1)
c                 2) for GIAO two-el. derivatives  (iforwhat=2)
c                 3) for gradient derivatives      (iforwhat=3)
c                 4) for second derivatives        (iforwhat=4)
c                    (used only in blocking routines)
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
      common /cpu/ intsize,iacc,icache,memreal
      common /intgop/ ncache,maxprice,iprint,iblock
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
c------------------------------------------------
C???  common /ener/ title(16),etot,enuc,qx(3),qtot,time,tlast,par(30)
      common /tape/ inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c-------------------------------------------------
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c-------------------------------------------------
      dimension xintxx(9)
      dimension bl(*)
      dimension inx(12,*)
c-------------------------------------------------
c
c prepare for two-electron integrals calculations
c
c-------------------------------------------------
c transfer all input parameters into local commons :
c
c responsible for blocking structure :
      ncache  =ncachex
c
c check-run and print level :
      iprint  =iprintx
ctest
c     icheck  =1
c     iprint  =1
ctest
      nprint  =iprint
      eps     =threshx
C???  par(1)  =eps
c
c-------------------------------
      icache=icache*ncache
c-------------------------------------------------------
c
      scftype='    '
      iforwhat=iforwhatx
      call prepint2(bl,eps,inuc,ibas,na,nbf,nsh,ncf,ncs,inx,
     *              lcore,nsym,lsemi,scftype)
c
c     if(iforwhat.eq.1) then
c        write(6,100) scftype
c        write(8,100) scftype
c     endif
c     if(iforwhat.eq.2) then
c        write(6,200) scftype
c        write(8,200) scftype
c     endif
c 100 format('===  SCF will be performed in the ',a11,' fashion  ===')
c 200 format('=== CPHF will be performed in the ',a11,' fashion  ===')
c
       scftypex=scftype
c-------------------------------------------------------
c these two are setup in blocking procedure executed by
c prepint2 :
      maxbuffer=maxibuf
      maxlabels=maxindx
      nblocks=nbl2*(nbl2+1)/2
c-------------------------------------------------------
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc0)
      if(iftc0.NE.0) then
         call getival('nftcbl4',nftcbl4)   ! pointer
         call getival('nbl4ftc',nbl4ftc)   ! no of FTC blocks
         nblocks=nbl4ftc
      endif
c-------------------------------------------------------
c-------------------------------------------------------
      end
