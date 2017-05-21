c==================================================================
      subroutine schw_neg(nbls,ncs, npij,npkl,npkl0, nijbeg,nklbeg,
     *                    nbl2,ibl,kbl,ijbl,iis,jjs, denspar,
c output
     *                    ipres, densmax,nblsp )
c---------------------------------------------------------------
c The SCREEN variable is used here , not WHERE
c---------------------------------------------------------------
c This routine is called from Calcint2 after pair-precalculations
c are completed. It neglects some of contracted shell quartets
c according to the Schwarz inequality (ij|kl)**2 <(ij|ij)*(kl|kl)
c Ipres(i) is 0 or 1 and shows if i-quartet can be negl.
c---------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
      character*4 screen
      common /screen_type/ screen
c     common /big/ bl(1)
      common /screening/ dens_max_el
      common /neglect/ eps,eps1,epsr
c---------------------------------------------------------------
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(*)   ! it is ncs*ncs or 2*ncs*ncs
      dimension densmax(nbls)
      dimension ipres(nbls)
      data zero /0.d0/
c---------------------------------------------------------------
      call getival('schwarz',ischwarz)
c---------------------------------------------------------------
c find  maximum of (ij|ij) & (kl|kl)
c this is already found in smblock_neg (calcint.f)
c
      call getrval('x_ij_ij',x_ij_ij_max)
      call getrval('x_kl_kl',x_kl_kl_max)
c---------------------------------------------------------------
c calculate an upper bond for the Estimator
c
      if(screen.eq.'fock') then
         call schw_neq_fock(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
      if(screen.eq.'coul') then
         call schw_neq_coul(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
      if(screen.eq.'forc') then
         call schw_neq_forc(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
      if(screen.eq.'shif') then    ! giao G(D0,g1)    OLD
         call schw_neq_shif(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
      if(screen.eq.'giao') then    ! giao G(D0,g1)    NEW
         call schw_neq_giao(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
      if(screen.eq.'cphf') then     ! NMR cphf    OLD
         call schw_neq_cphf(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
cnew
      if(screen.eq.'cph2') then     ! NMR cphf with f with 2-particle screening
         call schw_neq_cph2(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
cnew
      if(screen.eq.'2pde') then     ! HESS cphf with 2-particle screening
c                                   ! but using two densities
         call schw_neq_2pde(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
      if(screen.eq.'2pd1') then     ! HESS cphf with 2-particle screening
c                                   ! but using one density
         call schw_neq_2pd1(nbls,ncs,npij,npkl,npkl0,
     *              nijbeg,nklbeg, nbl2,ibl,kbl,ijbl,iis,jjs,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres,nblsp)
      endif
c---------------------------------------------------------------
c on return the Ipres(ijkl) has a value of 0 if it is to be neglect.
c and max density (squared) densmax() for each quartet .
c---------------------------------------------------------------
c
      end
c==================================================================
      subroutine schw_finder(ncs,npij,ijbegin,nbl2,ibl,ijbl,iis,jjs,
     *                       schwarz, x_ij_ij_max)
c
c this routine finds maximum element of Shwarz integrals
c
      implicit real*8 (a-h,o-z)
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension schwarz(ncs,ncs)
      data zero /0.d0/
c
c pairs have been reordered according to Schwarz so:
c
      ijpar=ijbegin
      ijcs=ijbl(ibl,ijpar)
      ics=iis(ijcs)
      jcs=jjs(ijcs)
      x_ij_ij_max=schwarz(ics,jcs)
c
c2000
c     x_ij_ij_max=zero
c     do ijc=1,npij
c        ijpar=ijbegin-1+ijc
c        ijcs=ijbl(ibl,ijpar)
c        ics=iis(ijcs)
c        jcs=jjs(ijcs)
c        xschw=schwarz(ics,jcs)
c        if(xschw.gt.x_ij_ij_max) x_ij_ij_max=xschw
c        endif
c     enddo
c2000
c
      end
c==================================================================
      subroutine get_dijkl(screen,dij,dik,dil,djk,djl,dkl,dmx)
c
c for everything except forces & hessian the denspar is already squared
c since we do screening in a quadratic form then we need factors=16
c
      implicit real*8 (a-h,o-z)
      character*4 screen
c
c for scf,schwarz:
c
      if(screen.eq.'fock') then
         dij_16=dij*16.d0      ! 16 is 4**2
         dkl_16=dkl*16.d0
         dmx=max(dij_16,dik,dil,djk,djl,dkl_16)
c
c for pure DFT (i.e., coulomb only)
c
      else if(screen.eq.'coul') then
         dmx = max(dij*16.0d0,dkl*16.0d0)    ! ** NEW  JB **
c
c for gradient :
c
      else if(screen.eq.'forc') then
         dmx=4.d0*dij*dkl + dik*djl + dil*djk
         dmx=dmx*dmx                         ! square for gradient
c
c for giao :
c
      else if(screen.eq.'shif') then
         dik_h=dik*0.5d0
         dil_h=dil*0.5d0
         djk_h=djk*0.5d0
         djl_h=djl*0.5d0
         dmx=max(dij,dik_h,dil_h,djk_h,djl_h,dkl)
c
c for nmr cphf :
c
      else if(screen.eq.'cphf') then
         dmx=max(dik,dil,djk,djl)
      endif
c
c for hessian cphf : when G(D1,g0) is calculated
      if(screen.eq.'2pde') then
         dijkl=sqrt(dij*dkl)
         dikjl=sqrt(dik*djl)
         diljk=sqrt(dil*djk)
         dmx=4.d0*dijkl + dikjl + diljk
         dmx=dmx*dmx
      endif
c
c for hessian :
c
      if(screen.eq.'hess') then
          write(6,*)' screen=hess  in use '
           STOP ' stopped '
         dij_16=dij*16.d0      ! 16 is 4**2
         dkl_16=dkl*16.d0
         dmx=16.d0*dij*dkl + dik*djl + dil*djk
      endif
c
      end
c==================================================================
      subroutine schw_neq_fock(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'fock' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
cccc  dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         dij=denspar(ics,jcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  dkl=denspar(kcs,lcs)
                  dik=denspar(ics,kcs)
                  dil=denspar(ics,lcs)
                  djk=denspar(jcs,kcs)
                  djl=denspar(jcs,lcs)
c
                  call get_dmx_fock(dij,dik,dil,djk,djl,dkl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_fock(dij,dik,dil,djk,djl,dkl,dmx)
      implicit real*8 (a-h,o-z)
c
c for scf,schwarz:
c
         dmx=max(4.d0*dij,dik,dil,djk,djl,4.d0*dkl)
         dmx=dmx*dmx
c
      end
c==================================================================
      subroutine schw_neq_coul(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'coul' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         dij=denspar(ics,jcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  dkl=denspar(kcs,lcs)
c
                  call get_dmx_coul(dij,dkl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_coul(dij,dkl,dmx)
c
c for everything except forces & hessian the dens is already squared
c since we do screening in a quadratic form then we need factors=16
c
      implicit real*8 (a-h,o-z)
c
         dmx = 4.0d0*max(dij,dkl)    ! ** NEW  JB **
         dmx=dmx*dmx
c
      end
c==================================================================
      subroutine schw_neq_forc(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'forc' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         dij=denspar(ics,jcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  dkl=denspar(kcs,lcs)
                  dik=denspar(ics,kcs)
                  dil=denspar(ics,lcs)
                  djk=denspar(jcs,kcs)
                  djl=denspar(jcs,lcs)
c
                  call get_dmx_forc(dij,dik,dil,djk,djl,dkl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_forc(dij,dik,dil,djk,djl,dkl,dmx)
c
c for everything except forces & hessian the dens is already squared
c since we do screening in a quadratic form then we need factors=16
c
      implicit real*8 (a-h,o-z)
c
         dmx=4.d0*dij*dkl + dik*djl + dil*djk
         dmx=dmx*dmx                         ! square for gradient
c
      end
c==================================================================
      subroutine schw_neq_shif(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'shif' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         dij=denspar(ics,jcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  dkl=denspar(kcs,lcs)
                  dik=denspar(ics,kcs)
                  dil=denspar(ics,lcs)
                  djk=denspar(jcs,kcs)
                  djl=denspar(jcs,lcs)
c
                  call get_dmx_shif(dij,dik,dil,djk,djl,dkl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_shif(dij,dik,dil,djk,djl,dkl,dmx)
c
c for everything except forces & hessian the dens is already squared
c since we do screening in a quadratic form then we need factors=16
c
      implicit real*8 (a-h,o-z)
c
         dmx1=max(dij,dkl)
         dmx2=0.50d0*max(dik,dil,djk,djl)
         dmx=max(dmx1,dmx2)
         dmx=dmx*dmx
c
      end
c==================================================================
      subroutine schw_neq_cphf(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'cphf' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  dik=denspar(ics,kcs)
                  dil=denspar(ics,lcs)
                  djk=denspar(jcs,kcs)
                  djl=denspar(jcs,lcs)
c
                  call get_dmx_cphf(dik,dil,djk,djl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_cphf(dik,dil,djk,djl,dmx)
c
c for everything except forces & hessian the dens is already squared
c since we do screening in a quadratic form then we need factors=16
c
      implicit real*8 (a-h,o-z)
c
      dmx=max(dik,dil,djk,djl)
      dmx=dmx*dmx
c
      end
c==================================================================
      subroutine schw_neq_2pde(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for '2pde' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs,2), schwarz(ncs,ncs)
c
c  denspar(ncs,ncs,1) is the current D1 (delta)
c  denspar(ncs,ncs,2) is the second screening density
c
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         d1ij=denspar(ics,jcs,1)
         d2ij=denspar(ics,jcs,2)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  d1kl=denspar(kcs,lcs,1)
                  d1ik=denspar(ics,kcs,1)
                  d1il=denspar(ics,lcs,1)
                  d1jk=denspar(jcs,kcs,1)
                  d1jl=denspar(jcs,lcs,1)
c
                  d2kl=denspar(kcs,lcs,2)
                  d2ik=denspar(ics,kcs,2)
                  d2il=denspar(ics,lcs,2)
                  d2jk=denspar(jcs,kcs,2)
                  d2jl=denspar(jcs,lcs,2)
c
                  call get_dmx_2pde(d1ij,d1ik,d1il,d1jk,d1jl,d1kl,
     *                              d2ij,d2ik,d2il,d2jk,d2jl,d2kl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_2pde(d1ij,d1ik,d1il,d1jk,d1jl,d1kl,
     *                        d2ij,d2ik,d2il,d2jk,d2jl,d2kl,dmx)
      implicit real*8 (a-h,o-z)
c for hessian cphf : when G(D1,g0) is calculated and for G(D0,g1)
c
         dijkl12=d1ij*d2kl
         dikjl12=d1ik*d2jl
         diljk12=d1il*d2jk
c
         dijkl21=d2ij*d1kl
         dikjl21=d2ik*d1jl
         diljk21=d2il*d1jk
c
         dmx12=4.d0*dijkl12 + dikjl12 + diljk12
         dmx21=4.d0*dijkl21 + dikjl21 + diljk21
c
         dmx=max(dmx12,dmx21)
c
         dmx=dmx*dmx
c
      end
c==================================================================
c==================================================================
      subroutine schw_neq_cph2(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'cphf' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs,2), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  d1ik=denspar(ics,kcs,1)
                  d1il=denspar(ics,lcs,1)
                  d1jk=denspar(jcs,kcs,1)
                  d1jl=denspar(jcs,lcs,1)
c
                  d2ik=denspar(ics,kcs,2)
                  d2il=denspar(ics,lcs,2)
                  d2jk=denspar(jcs,kcs,2)
                  d2jl=denspar(jcs,lcs,2)
c
                  call get_dmx_cph2(d1ik,d1il,d1jk,d1jl,
     *                              d2ik,d2il,d2jk,d2jl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_cph2(d1ik,d1il,d1jk,d1jl,
     *                        d2ik,d2il,d2jk,d2jl,dmx)
      implicit real*8 (a-h,o-z)
c
         dikjl12=d1ik*d2jl
         diljk12=d1il*d2jk
c
         dikjl21=d2ik*d1jl
         diljk21=d2il*d1jk
c
         dmx12= dikjl12 + diljk12
         dmx21= dikjl21 + diljk21
         dmx=max(dmx12,dmx21)
         dmx=dmx*dmx
c
      end
c==================================================================
      subroutine schw_neq_giao(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c screen   : for 'giao' type screening
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs,2), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         d1ij=denspar(ics,jcs,1)
         d2ij=denspar(ics,jcs,2)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  d1kl=denspar(kcs,lcs,1)
                  d1ik=denspar(ics,kcs,1)
                  d1il=denspar(ics,lcs,1)
                  d1jk=denspar(jcs,kcs,1)
                  d1jl=denspar(jcs,lcs,1)
c
                  d2kl=denspar(kcs,lcs,2)
                  d2ik=denspar(ics,kcs,2)
                  d2il=denspar(ics,lcs,2)
                  d2jk=denspar(jcs,kcs,2)
                  d2jl=denspar(jcs,lcs,2)
c
                  call get_dmx_giao(d1ij,d1ik,d1il,d1jk,d1jl,d1kl,
     *                              d2ij,d2ik,d2il,d2jk,d2jl,d2kl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_giao(d1ij,d1ik,d1il,d1jk,d1jl,d1kl,
     *                        d2ij,d2ik,d2il,d2jk,d2jl,d2kl,dmx)
c
      implicit real*8 (a-h,o-z)
c
c        dik_h=dik*0.5d0
c        dil_h=dil*0.5d0
c        djk_h=djk*0.5d0
c        djl_h=djl*0.5d0
c        dmx=max(dij,dik_h,dil_h,djk_h,djl_h,dkl)
c
c        dmx1=max(dij,dkl)
c        dmx2=0.50d0*max(dik,dil,djk,djl)
c        dmx=max(dmx1,dmx2)
c        dmx=dmx*dmx
c
         dmx1=4.d0*d1ij*d2kl + d1ik*d2jl + d1il*d2jk
         dmx2=4.d0*d2ij*d1kl + d2ik*d1jl + d2il*d1jk
c
         dmx=max(dmx1,dmx2)
         dmx=0.5d0*dmx
         dmx=dmx*dmx
c
      end
c==================================================================
      subroutine schw_neq_2pd1(nbls,ncs,npij,npkl,npkl0,
     *                  nijbeg,nklbeg,nbl2,ibl,kbl,ijbl,iis,jjs,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres,nblsp)
c------------------------------------------------------------------
c INPUT :
c nbls     : block size (number of quartets)
c npij     : number of IJcs pairs in this block of quartets
c npkl     : number of KLcs pairs in this block of quartets
c npkl0    : shows if IJcs & KLcs pairs are the same (npkl0=0)
c screen   : shows type of screening
c denspar  : density over contracted shells
c schwarz  : schwarz integrals over contracted shells
c x_ij_ij_max : maximum of (ij|ij) in the current block
c x_kl_kl_max : maximum of (kl|kl) in the current block
c
c OUTPUT
c
c ipres(nbls)  : is =0 if a given quartet is to be neglected
c densmax(nbls): maximum density(squared) element for prescreening
c nblsp : is =0 if the whole block is skipped
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*),ijbl(nbl2,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
c
      dens1max=dens_max_el
c
c------------------------------------------------------------------
      nblsp=nbls
c------------------------------------------------------------------
ccc   dmx_overall=0.d0
      dmx_overall=dens_max_el
c
      do ijc=1,npij
         ijpar=nijbeg-1+ijc
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         dij=denspar(ics,jcs)
c
         if(npkl0.eq.0)then
            npkl_e=ijc
            ijij=ijc*(ijc-1)/2
         else
            npkl_e=npkl
            ijij=(ijc-1)*npkl
         endif
c
         x_ij_ij=schwarz(ics,jcs)
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           .......................................................
c           all next ijpairs will have x_ij_ij even smaller so
c           neglect all quatrtets from here to the end and exit
c
            do iqrt=ijij+1,nbls
               ipres(iqrt)=0
            enddo
            nblsp=nblsp-(nbls-ijij)
c
            go to 7777
c           .......................................................
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=ijbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) then
                  nblsp=nblsp-1
                  go to 9999
c                 this can be true for giao,gard NOT for scf
               endif
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
c                 .............................
c                 all next klpairs will have even smaller x_kl_kl so
c                 neglect quartets from here to the end of kl-block
c                 and exit kl loop
                  do iqrt=klc,npkl_e
                     ipres(ijij+iqrt)=0
                  enddo
                  nblsp=nblsp-(npkl_e-klc+1)
                  go to 8888
c                 .............................
               else
                  dkl=denspar(kcs,lcs)
                  dik=denspar(ics,kcs)
                  dil=denspar(ics,lcs)
                  djk=denspar(jcs,kcs)
                  djl=denspar(jcs,lcs)
c
                  call get_dmx_2pd1(dij,dik,dil,djk,djl,dkl,dmx)
c
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                     nblsp=nblsp-1
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
c
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
         endif                      ! (estim_max.LT.eps2) ij
 8888    continue
      enddo                         ! ijc=1,npij
c------------------------------------------------------------------
 7777 continue
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine get_dmx_2pd1(dij,dik,dil,djk,djl,dkl,dmx)
      implicit real*8 (a-h,o-z)
c
c for hessian cphf : when G(D1,g0) is calculated
c
      dijkl=dij*dkl
      dikjl=dik*djl
      diljk=dil*djk
      dmx=4.d0*dijkl + dikjl + diljk
      dmx=dmx*dmx
c
      end
c==================================================================
