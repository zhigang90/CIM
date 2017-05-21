      subroutine ptwoint1(nfock,  threshx,threshy,iforwhat,idft,
     *                   ax,     nrad,   nang,   lrad,   lang,
     *                   Iradq,  NBatch, rhf,    NAlpha, NBeta,
     *                   wherex, xintxx, nblocks,maxbuffer,maxlabels)
      implicit real*8(a-h,o-z)
c
c Input :
c
c  1. nfock    -  number of Fock matrices calculated simultaneously
c  2. threshx  -  integral threshold
c  3. threshy  -  second, sharper integral threshold
c  4. iforwhat -  shows the type of integrals to be calculated;
c                  1) for ordinary two-el.integrals (iforwhat=1)
c                  2) for GIAO two-el. derivatives  (iforwhat=2)
c                  3) for gradient derivatives      (iforwhat=3)
c                  4) for second derivatives        (iforwhat=4)
c                  5) for MP2 integrals             (iforwhat=5)
c                     it will suppress symmetry and it will be 1
c                     in the further calls
c                     (used only in blocking routines)
c ......................................................................
c Following arguments all need to be communicated to slaves for DFT
c
c  5. idft    -  dft flag
c  6. ax      -  amount of HF exchange included in DFT functional
c  7. nrad    -  maximum number of radial grid points per atom
c  8. nang    -  maximum number of angular grid points per shell
c  9. lrad    -  number of radial grid points for numerical grid
c             (if zero, will be determined in grid generation routines)
c 10. lang    -  angular quadrature order for numerical grid
c             (if zero, will be determined in grid generation routines)
c 11. IradQ   -  radial quadrature type
c 12. NBatch  -  maximum number of grid points treated in one batch
c 13. rhf     -  logical flag for closed-shell RHF
c 14. NAlpha  -  number of occupied alpha/closed-shell MOs
C 15. NBeta   -  number of occupied beta MOs
c ......................................................................
c
c Output :
c
c 16. wherex   - scftype (apparently not used)
c 17. xintxx(9)- contains number of integrals more expensive than :
c               0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75%, 99% of maximum price
c 18. nblocks  - total number of blocks of contracted shell quartets
c                 (super-blocks)
c 19. maxbuffer - maximum size for integral buffer , calc in blksizer
c 20. maxlabels - maximum size for labels info     , calc in blksizer
c---------------------------------------------------------------------
      character*(*) wherex
      Logical rhf
      dimension xintxx(9)
c---------------------------------------------------------------------
c get parameters for twoint95 & pvm :
c
      call getival('lcore',lcorex)
      call getival('nsym',nsymx)
      call getival('maxprice',maxpricex)
      call getival('ncache',ncachex)
      call getival('iroute',iroutex)
      call getival('chec',icheckx)
      call getival('iprn',iprintx)
      call getival('istat',istat)
c
c   suppress symmetry in the MP2 integrals
c
      iforwhatx=iforwhat
ckw   if(iforwhat.eq.5) nsymx=0
c---------------------------------------------------------------------
c     Initialize local integral system.
c
      call twoint95(lcorex,nsymx,maxpricex,ncachex,iroutex,
     *              icheckx,iprintx,threshx,istat,iforwhatx,
     *              0,nblocks,maxbuffer,maxlabels,
     *                 wherex,xintxx)
c
c For LMP2 integ. iforwhatx has value of 1 on return
c---------------------------------------------------------------------
c
c Start and initialize slaves
c Arguments are the data to be sent to the slaves + nblocks
c
c The following are not sent -  let me tell why
c lcore    - slaves have their own memory size
c maxprice - only the original value is used now, it is
c            not used at all in post_scf calcs.
c            will be sent when functional
c iroute   - always takes the value setup in setdata, in
c            post_scf keeps the scf value
c icheck   - checkruns and non-zero printlevels do not make sense
c iprint   - in slaves, since the output is sent mixed to the master
c            slaves will always use the value from setdata
c            if needed one has to poke into the program anyway
c istat    - it's always the value from setdata
c (nsym is sent with the basic data)
c
c
c the following call is to an empty toutine, so we might as well
c comment it out MM
c
c     call para_twoint1(iforwhat,ncachex,threshx,threshy,nfock,
c    $                 nblocks, idft,   ax,     nrad,   nang,
c    $                 lrad,    lang,   IradQ,  NBatch, rhf,
c    $                 NAlpha,  NBeta)
c
      end
c=====================================================================
c     subroutine para_twoint1(iforwhat,ncachex,threshx,threshy,
c    *                       nfock,nblocks,idft,ax,nrad,nang,
c    *                       lrad,lang,IradQ,NBatch,rhf,NAlpha,NBeta)
c     implicit real*8 (a-h,o-z)
c     end
c======================================================================
