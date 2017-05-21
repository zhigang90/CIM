      subroutine ptwoint(nfock,  threshx,threshy,iforwhat,idft,
     *                   ax,     nrad,   nang,   lrad,   lang,
     *                   Iradq,  NBatch, rhf,    NAlpha, NBeta,
     *                   lsemi,  scftype,xintxx, nblocks,maxbuffer,
     *                   maxlabels)
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
c                 11) for classical FTC integrals   (iforwhat=11)
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
C 16. lsemi   -  integer flag to use semi-direct if available
C                 0 - do not use; otherwise use
c ......................................................................
c
c Output :
c
c 17. scftype  - scf type - either "full-direct" or "semi-direct"
c 18. xintxx(9)- contains number of integrals more expensive than :
c               0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75%, 99% of maximum price
c 19. nblocks  - total number of blocks of contracted shell quartets
c                 (super-blocks)
c 20. maxbuffer - maximum size for integral buffer , calc in blksizer
c 21. maxlabels - maximum size for labels info     , calc in blksizer
c---------------------------------------------------------------------
      character*(*) scftype
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
     *              lsemi,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
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
      call para_twoint(iforwhat,ncachex,threshx,threshy,nfock,
     $                 lsemi,   nblocks, idft,   ax,    nrad,
     $                 nang,    lrad,    lang,   IradQ, NBatch,
     $                 rhf,     NAlpha,  NBeta)
c
      end
c=====================================================================
c
      subroutine concat(namei,proc)
c     assembles a file name from namei and index.
c     the file name is returned
c     in file name and the length of the file name is returned in ile.
      character*256 namei,proc


c     determine the last character in name
      k=1
      do while(namei(k:k).ne.' ')
         k=k+1
      end do
c
c     find the end of the proc string
c
      n=0
      do while(proc(n+1:n+1).ne.' ')
         n=n+1
      end do
c
c     concatenate the two strings
c
      namei(k:k+n)=proc(1:n)
c
      end
