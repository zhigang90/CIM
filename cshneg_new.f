c==================================================================
      subroutine schw_neg_new(nbls,ncs,iis,jjs,
     *                        ibl,ijbl,nbl2_ijd,nijbeg,npij,
     *                        kbl,klbl,nbl2_kld,nklbeg,npkl,npkl0,
     *                        denspar,
     *                        ipres, densmax )   !   output
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
      common /neglect/ eps,eps1,epsr,eps8
c---------------------------------------------------------------
      dimension iis(*),jjs(*)
      dimension ijbl(nbl2_ijd,*)
      dimension klbl(nbl2_kld,*)
      dimension denspar(ncs,ncs)
      dimension densmax(nbls)
      dimension ipres(nbls)
      data zero /0.d0/
c---------------------------------------------------------------
      call getival('schwarz',ischwarz)
c---------------------------------------------------------------
c find  maximum of (ij|ij) & (kl|kl) for this block
c
      call schw_finder(ncs,npij,nijbeg,nbl2_ijd,ibl,ijbl,iis,jjs,
     *                 bl(ischwarz),x_ij_ij_max )
      if(npkl0.ne.0) then
         call schw_finder(ncs,npkl,nklbeg,nbl2_kld,kbl,klbl,iis,jjs,
     *                    bl(ischwarz),x_kl_kl_max )
      else
         x_kl_kl_max=x_ij_ij_max
      endif
c---------------------------------------------------------------
c calculate an upper bond for the Estimator
c
      call schw_neq_new(nbls,ncs,iis,jjs,screen,
     *                  ibl,ijbl,nbl2_ijd,nijbeg,npij,
     *                  kbl,klbl,nbl2_kld,nklbeg,npkl,npkl0,
     *              denspar,bl(ischwarz),x_ij_ij_max,x_kl_kl_max,
     *              densmax,ipres)
c
c on return the Ipres(ijkl) has a value of 0 if it is to be neglect.
c and max density (squared) densmax() for each quartet .
c---------------------------------------------------------------
c
      end
c==================================================================
      subroutine schw_neq_new(nbls,ncs,iis,jjs,screen,
     *                        ibl,ijbl,nbl2_ijd,nijbeg,npij,
     *                        kbl,klbl,nbl2_kld,nklbeg,npkl,npkl0,
     *                  denspar,schwarz, x_ij_ij_max, x_kl_kl_max,
     *                  densmax,ipres)
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
c OUTPUT :
c
c ipres(nbls)  : is =0 if a given quartet is to neglected
c densmax(nbls): maximum density(squared) element for prescreening
c------------------------------------------------------------------
c densmax(ijkl) will be used in prec4neg & precspec
c               ONLY if ipres(ijkl)>0
c------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*4 screen
c---------------------------------------------------------------
      common /neglect/ eps,eps1,epsr,eps8
      common /screening/ dens_max_el                 ! input
      dimension iis(*),jjs(*)
      dimension ijbl(nbl2_ijd,*)
      dimension klbl(nbl2_kld,*)
      dimension denspar(ncs,ncs), schwarz(ncs,ncs)
      dimension densmax(nbls)                        ! output
      dimension ipres(nbls)                          ! output
      data zero,four /0.d0 , 4.d0/
c------------------------------------------------------------------
ccc   write(6,*)' cshneg : screen=',screen
c------------------------------------------------------------------
      eps2=eps*eps
c------------------------------------------------------------------
      dens1max=dens_max_el
c------------------------------------------------------------------
c Preliminary selection of quartets for further calculations :
c
      estim_max=dens1max * x_ij_ij_max * x_kl_kl_max
      if(estim_max.LT.eps2 ) then
c        write(91,*)' the whole block skipped'
c        all quartets to be neglected
         do ijkl=1,nbls
            ipres(ijkl)=0
         enddo
         RETURN
      endif
c------------------------------------------------------------------
      dmx_overall=0.d0
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
c
         dens1mx_ij=dens1max * x_ij_ij
         estim_max=dens1mx_ij * x_kl_kl_max
         if(estim_max.LT.eps2) then
c           neglect quatrtets from this ijcs pair and ALL klcs pairs
            do ijc1=ijc,npij
               if(npkl0.eq.0)then
                  npkl_e1=ijc1
                  ijij1=ijc1*(ijc1-1)/2
               else
                  npkl_e1=npkl
                  ijij1=(ijc1-1)*npkl
               endif
               do iqrt=1,npkl_e1
                  ipres(ijij1+iqrt)=0
               enddo
            enddo
            go to 7777      ! exit ij loop
         else
            do klc=1,npkl_e
               klpar=nklbeg-1+klc
               klcs=klbl(kbl,klpar)
               kcs=iis(klcs)
               lcs=jjs(klcs)
               ijkl=ijij+klc
               if(ipres(ijkl).eq.0) go to 9999  ! already 0  ????
c
               x_kl_kl=schwarz(kcs,lcs)
               estim_max=dens1mx_ij * x_kl_kl
               if(estim_max.LT.eps2) then
                  do iqrt=klc,npkl_e
                      ipres(ijij+iqrt)=0
                  enddo
                  go to 8888      ! exit kl loop
               else
                  call get_maxscreen(denspar,ncs,ics,jcs,kcs,lcs,screen,
     *                               dmx)
                  estim_max=dmx*x_ij_ij * x_kl_kl
                  if(estim_max.LT.eps2) then
                     ipres(ijkl)=0
                  else
                     densmax(ijkl)=dmx
                     if(dmx.gt.dmx_overall) dmx_overall=dmx
                  endif
               endif                ! (estim_max.LT.eps2) kl
 9999          continue
            enddo                   ! klc=1,npkl_e
 8888       continue
         endif                      ! (estim_max.LT.eps2) ij
      enddo                         ! ijc=1,npij
 7777 continue
c
c------------------------------------------------------------------
c save dmx_overall
c
      call setrval('dmx_over',dmx_overall)
c
c will be used in get_max4kl
c------------------------------------------------------------------
      end
c==================================================================
      subroutine dmax_find_4b(isupbl,bl,denspar,ncs,screen)
c-------------------------------------------------------
c input :
c isupbl - number of a current super-block
c bl(*)  - scratch buffer
c denspar(ncs,ncs) - density matrix over cont.shells
c ncs  - number of cont.shells
c
c output : dens_max kept in common /screening/
c-------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      character*4 screen
      common /screening/ dens_max_el      ! output
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1b/ nbl2,nbloks
      dimension bl(*)
      dimension denspar(ncs,ncs)
      common /infob/ inucx,ibasx,nax,nbfx,nshx,ncfx,ncsx
c----------------------------------------------------------------
      natoms=nax
      if(natoms.le.50) then
         dens_max_el=1.d0
         return
      endif
c----------------------------------------------------------------
      if(screen.eq.'mp2d') then
         call getrval('dmx_mp2d',dsmax)
         dens_max_el=dsmax
cccc     write(6,*)' from cshneg: block=',isupbl,' densmax=',dsmax
         return
      endif
c----------------------------------------------------------------
c get pair-blocks making up this quartet-block :
c
      call getival('nblock4',nblock4)
      call get_pair_blocks(bl(nblock4),isupbl,ibl,kbl)
c----------------------------------------------------------------
c get data concerning pair-blocks :
c
         call getival('nparx',npar_ij)
         call getival('ijblx',ijbl)
c....... call getival('mapijblx',map_ij_bl2_ij)
         call getival('blocksij',nbl2_ij)
         call getival('blpredij',nbl2_ijd)
c
         call getival('npary',npar_kl)
         call getival('ijbly',klbl)
c....... call getival('mapijbly',map_ij_bl2_kl)
         call getival('blockskl',nbl2_kl)
         call getival('blpredkl',nbl2_kld)
c----------------------------------------------------------------
c for these two pair-blocks find shell ranges :
c
      call get_shell_range_n(bl(iisd),bl(jjsd),
     *                       ibl,bl(ijbl),nbl2_ijd,bl(npar_ij),
     *                       kbl,bl(klbl),nbl2_kld,bl(npar_kl),
     *                       ics_b,ics_e, jcs_b,jcs_e,
     *                       kcs_b,kcs_e, lcs_b,lcs_e)
c
c find maximum element of density matrix within above shell range :
c
      if(screen.eq.'mp2 ') then
         call dmax_finder_n(ncs,denspar,'djl_max',
     *                      ics_b,ics_e, jcs_b,jcs_e,
     *                      kcs_b,kcs_e, lcs_b,lcs_e,
     *                      dij,dik,dil,djk,djl,dkl)
      else
         call dmax_finder_n(ncs,denspar,'allijkl',
     *                      ics_b,ics_e, jcs_b,jcs_e,
     *                      kcs_b,kcs_e, lcs_b,lcs_e,
     *                      dij,dik,dil,djk,djl,dkl)
      endif
c
      if(screen.eq.'forc') then
c        dij,dik,dil,djk,djl,dkl are NOT yet squared
c
         dmx=4.d0*dij*dkl + dik*djl + dil*djk
         dmx=dmx*dmx       ! square it for gradient
cc
      else if(screen.eq.'fock') then
c        dij,dik,dil,djk,djl,dkl are already squared
c
         dij_16=dij*16.d0      ! 16 is 4**2
         dkl_16=dkl*16.d0
         dmx=max(dij_16,dik,dil,djk,djl,dkl_16)
cc
      else if(screen.eq.'coul') then
         dmx = max(dij*16.0d0,dkl*16.0d0)    ! ** NEW  JB **
cc
      else if(screen.eq.'mp2 ') then
         dmx=djl
cc
      else if(screen.eq.'shif') then
         dik_h=dik*0.5d0
         dil_h=dil*0.5d0
         djk_h=djk*0.5d0
         djl_h=djl*0.5d0
         dmx=max(dij,dik_h,dil_h,djk_h,djl_h,dkl)
cc
      else if(screen.eq.'cphf') then
         dmx=max(dik,dil,djk,djl)
c
c for hessian (not ready yet) :
c
      else if(screen.eq.'hess') then
         call nerror(1,'dmax_find_4b','not ready for hessian',0,0)
         dmx=max(dij,dik,dil,djk,djl,dkl)
      endif
c
c finally:
c
      dens_max_el=dmx
c
      end
c=======================================================
      subroutine get_shell_range_n(iis,jjs,
     *                             ibl,ijbl,nbl2_ijd,npar_ij,
     *                             kbl,klbl,nbl2_kld,npar_kl,
     *                             ics_b,ics_e, jcs_b,jcs_e,
     *                             kcs_b,kcs_e, lcs_b,lcs_e)
      dimension iis(*),jjs(*)
      dimension ijbl(nbl2_ijd,*), npar_ij(*)
      dimension klbl(nbl2_kld,*), npar_kl(*)
c
      nparij=npar_ij(ibl)
      nijbeg=1
      call shel_finder_n(nparij,nijbeg,nbl2_ijd,ibl,ijbl,iis,jjs,
     *                   ics_b,ics_e, jcs_b,jcs_e)
      nparkl=npar_kl(kbl)
      nklbeg=1
      call shel_finder_n(nparkl,nklbeg,nbl2_kld,kbl,klbl,iis,jjs,
     *                   kcs_b,kcs_e, lcs_b,lcs_e)
c
      end
c==================================================================
      subroutine shel_finder_n(npij,ijbegin,nbl2,ibl,ijbl,iis,jjs,
     *                         ics_b,ics_e, jcs_b,jcs_e)
c
c this routine finds range of ics & jcs-shells in a pair-block .
c
      dimension iis(*),jjs(*),ijbl(nbl2,*)
c
      ijpar_b=ijbegin
      ijpar_e=ijbegin-1+npij
c
      call minimax_ij(ibl,ijbl,nbl2,ijpar_b,ijpar_e,iis,jjs,
     *                ics_b,ics_e, jcs_b,jcs_e)
c
      end
c==================================================================
c this routine finds range of ics & jcs-shells in a pair-block .
c
      subroutine minimax_ij(ibl,ijbl,nbl2,ijpar_b,ijpar_e,iis,jjs,
     *                      ics_b,ics_e, jcs_b,jcs_e)
      dimension iis(*),jjs(*),ijbl(nbl2,*)
c
      ics_min= 100 000
      ics_max=-100 000
      jcs_min= 100 000
      jcs_max=-100 000
      do ijpar=ijpar_b,ijpar_e
         ijcs=ijbl(ibl,ijpar)
         ics=iis(ijcs)
         jcs=jjs(ijcs)
         if(ics.gt.ics_max) ics_max=ics
         if(ics.lt.ics_min) ics_min=ics
         if(jcs.gt.jcs_max) jcs_max=jcs
         if(jcs.lt.jcs_min) jcs_min=jcs
      enddo
c
      ics_b=ics_min
      ics_e=ics_max
      jcs_b=jcs_min
      jcs_e=jcs_max
c
      end
c=======================================================
      subroutine dmax_finder_n(ncs,denspar,what2find,
     *                         ics_b,ics_e, jcs_b,jcs_e,
     *                         kcs_b,kcs_e, lcs_b,lcs_e,
     *             dij_max,dik_max,dil_max,djk_max,djl_max,dkl_max)
      implicit real*8 (a-h,o-z)
      character*7 what2find
c it can be : 'allijkl' ,
c             'dij_max', 'dik_max', 'dil_max',
c             'djk_max', 'djk_max', 'dkl_max',
c
      dimension denspar(ncs,ncs)
      data zero /0.d0/
c
      dij_max=zero
      dik_max=zero
      dil_max=zero
      djk_max=zero
      djl_max=zero
      dkl_max=zero
c
      IF(WHAT2FIND.EQ.'allijkl') THEN
c          ij,ik,il :
         do ics=ics_b,ics_e
            do jcs=jcs_b,jcs_e
               dij=denspar(ics,jcs)
               if(dij.gt.dij_max) dij_max=dij
            enddo
            do kcs=kcs_b,kcs_e
               dik=denspar(ics,kcs)
               if(dik.gt.dik_max) dik_max=dik
            enddo
            do lcs=lcs_b,lcs_e
               dil=denspar(ics,lcs)
               if(dil.gt.dil_max) dil_max=dil
            enddo
         enddo
c          jk,jl :
         do jcs=jcs_b,jcs_e
            do kcs=kcs_b,kcs_e
               djk=denspar(jcs,kcs)
               if(djk.gt.djk_max) djk_max=djk
            enddo
            do lcs=lcs_b,lcs_e
               djl=denspar(jcs,lcs)
               if(djl.gt.djl_max) djl_max=djl
            enddo
         enddo
c
c          kl :
         do kcs=kcs_b,kcs_e
            do lcs=lcs_b,lcs_e
               dkl=denspar(kcs,lcs)
               if(dkl.gt.dkl_max) dkl_max=dkl
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'allijkl') THEN
c
      IF(WHAT2FIND.EQ.'dij_max') THEN
c        ij :
         do ics=ics_b,ics_e
            do jcs=jcs_b,jcs_e
               dij=denspar(ics,jcs)
               if(dij.gt.dij_max) dij_max=dij
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'dij_max') THEN
c
      IF(WHAT2FIND.EQ.'dik_max') THEN
c        ik :
         do ics=ics_b,ics_e
            do kcs=kcs_b,kcs_e
               dik=denspar(ics,kcs)
               if(dik.gt.dik_max) dik_max=dik
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'dik_max') THEN
      IF(WHAT2FIND.EQ.'dil_max') THEN
c        il :
         do ics=ics_b,ics_e
            do lcs=lcs_b,lcs_e
               dil=denspar(ics,lcs)
               if(dil.gt.dil_max) dil_max=dil
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'dil_max') THEN
      IF(WHAT2FIND.EQ.'djk_max') THEN
c          jk :
         do jcs=jcs_b,jcs_e
            do kcs=kcs_b,kcs_e
               djk=denspar(jcs,kcs)
               if(djk.gt.djk_max) djk_max=djk
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'djk_max') THEN
      IF(WHAT2FIND.EQ.'djl_max') THEN
c          jl :
         do jcs=jcs_b,jcs_e
            do lcs=lcs_b,lcs_e
               djl=denspar(jcs,lcs)
               if(djl.gt.djl_max) djl_max=djl
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'djl_max') THEN
      IF(WHAT2FIND.EQ.'dkl_max') THEN
c          kl :
         do kcs=kcs_b,kcs_e
            do lcs=lcs_b,lcs_e
               dkl=denspar(kcs,lcs)
               if(dkl.gt.dkl_max) dkl_max=dkl
            enddo
         enddo
         return
      ENDIF          !       IF(WHAT2FIND.EQ.'dkl_max') THEN
c
      end
c==================================================================
      subroutine get_maxscreen(denspar,ncs,ics,jcs,kcs,lcs,screen,dmx)
      implicit real*8 (a-h,o-z)
      character*4 screen
      common /mp2shells/ ics_mp2, kcs_mp2 ! used here and in calcint2.f
      dimension denspar(ncs,ncs)
c....................................................................
c for everything except forces the denspar is already squared
c since we do screening in a quadratic form then we need factors=16
c....................................................................
c
c for scf:
c
      if(screen.eq.'fock') then
         dij=denspar(ics,jcs)
         dkl=denspar(kcs,lcs)
         dik=denspar(ics,kcs)
         dil=denspar(ics,lcs)
         djk=denspar(jcs,kcs)
         djl=denspar(jcs,lcs)
         dij_16=dij*16.d0      ! 16 is 4**2
         dkl_16=dkl*16.d0
         dmx=max(dij_16,dik,dil,djk,djl,dkl_16)
         return
      endif
c
c for gradient :
c
      if(screen.eq.'forc') then
         dij=denspar(ics,jcs)
         dkl=denspar(kcs,lcs)
         dik=denspar(ics,kcs)
         dil=denspar(ics,lcs)
         djk=denspar(jcs,kcs)
         djl=denspar(jcs,lcs)
         dmx=4.d0*dij*dkl + dik*djl + dil*djk
         dmx=dmx*dmx                        ! squar for gradient
         return
      endif
c
c for pure DFT (coulomb only)
c
      if(screen.eq.'coul') then
         dij = denspar(ics,jcs)
         dkl = denspar(kcs,lcs)
         dmx = max(dij,dkl)*16.0d0
         return
      endif
c
c for giao :
c
      if(screen.eq.'shif') then
         dij=denspar(ics,jcs)
         dkl=denspar(kcs,lcs)
         dik=denspar(ics,kcs)
         dil=denspar(ics,lcs)
         djk=denspar(jcs,kcs)
         djl=denspar(jcs,lcs)
         dik_h=dik*0.5d0
         dil_h=dil*0.5d0
         djk_h=djk*0.5d0
         djl_h=djl*0.5d0
         dmx=max(dij,dik_h,dil_h,djk_h,djl_h,dkl)
         return
      endif
c
c for nmr cphf :
c
      if(screen.eq.'cphf') then
         dik=denspar(ics,kcs)
         dil=denspar(ics,lcs)
         djk=denspar(jcs,kcs)
         djl=denspar(jcs,lcs)
         dmx=max(dik,dil,djk,djl)
         return
      endif
c
c for hessian (not used) :
c
      if(screen.eq.'hess') then
        call nerror(1,'get_maxscreen','no hess type screening',0,0)
         dij=denspar(ics,jcs)
         dkl=denspar(kcs,lcs)
         dik=denspar(ics,kcs)
         dil=denspar(ics,lcs)
         djk=denspar(jcs,kcs)
         djl=denspar(jcs,lcs)
         dmx=max(dij,dik,dil,djk,djl,dkl)
         return
      endif
c
c
c screening for MP2 integrals and MP2 derivatives :
c
      if(screen.eq.'mp2 '.or.screen.eq.'lmp2'.or.screen.eq.'mp2d') then
c        we need only JL element but since shells could be switched
c        we have to take a right element :
c        ...........................................................
                  if( ((ics.EQ.ics_mp2).and.(kcs.EQ.kcs_mp2)) .or.
     *                ((kcs.EQ.ics_mp2).and.(ics.EQ.kcs_mp2)) ) then
ccccc                dmx=denspar(jcs,lcs)
                     dmx=denspar(lcs,jcs)
                  endif
                  if( ((ics.EQ.ics_mp2).and.(lcs.EQ.kcs_mp2)) .or.
     *                ((lcs.EQ.ics_mp2).and.(ics.EQ.kcs_mp2)) ) then
ccccc                dmx=denspar(jcs,kcs)
                     dmx=denspar(kcs,jcs)
                  endif
                  if( ((jcs.EQ.ics_mp2).and.(kcs.EQ.kcs_mp2)) .or.
     *                ((kcs.EQ.ics_mp2).and.(jcs.EQ.kcs_mp2)) ) then
ccccc                dmx=denspar(ics,lcs)
                     dmx=denspar(lcs,ics)
                  endif
                  if( ((jcs.EQ.ics_mp2).and.(lcs.EQ.kcs_mp2)) .or.
     *                ((lcs.EQ.ics_mp2).and.(jcs.EQ.kcs_mp2)) ) then
ccccc                dmx=denspar(ics,kcs)
                     dmx=denspar(kcs,ics)
                  endif
c        ...........................................................
         return
      endif
c
      end
