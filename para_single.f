      subroutine para_driver
      call driver
      end
c======================================================================
      subroutine para_stop
      end
c======================================================================
      subroutine para_mem(lcore,incore,indisk)
      end
c======================================================================
      subroutine para_next(mgo)
      end
c======================================================================
      subroutine para_done(mgo)
      end
c======================================================================
      subroutine para_jobinit(iforwhat)
      end
c======================================================================
      subroutine para_oneint(itype,  nna,    oneint, inx,    kk1,
     $                       kk2,    basdat, datnuc, ncs)
      implicit real*8 (a-h,o-z)
      call inton(itype,  nna,    oneint, inx,    kk1,
     $           kk2,    basdat, datnuc, ncs)
      end
c======================================================================
      subroutine para_twoint(iforwhat,ncachex,threshx,threshy,nfock,
     *                       lsemi,nblocks,idft,ax,nrad,nang,
     *                       lrad,lang,IradQ,NBatch,rhf,NAlpha,NBeta)
      implicit real*8 (a-h,o-z)
      end
c======================================================================
      subroutine get_storage(iwher)
      end
c======================================================================
      subroutine para_store(ifrom,ito,bl,inx,dens,thres,labels,
     *                      iwher,isecond)
      call int_store(ifrom,ito,bl,inx,dens,thres,labels,
     1                isecond)
      end
c======================================================================
      subroutine send_limits(isecond,ito)
      end
c======================================================================
      subroutine para_fock(idft,ax,nblocks,nfock,rhf,ncf,bl,inx,
     *                     thres,dens,fock,fockB,bdn,DenB,
     *                     labels,mywork,igran,mgo,iforwhat)
      implicit real*8 (a-h,o-z)
      call int_fock(idft,ax,nblocks,nfock,rhf,ncf,inx,thres,
     *              dens,fock,fockB,bdn,DenB,labels,mywork,igran)
      end
c======================================================================
      subroutine para_preinit_FTC(
     &       expaccl,   dist2mp,   isharpgrd, isharpness, ncspl,
     &       ncfpl,    iicsplsize, icssharps, iicsdiff,   ilistsd,
     &       iicspltype,iplbasdat, imultipolecut,nomultipole)
      implicit real*8(a-h,o-z)
      integer expaccl
      call preinit_FTC(
     &       expaccl,   dist2mp,   isharpgrd, isharpness, ncspl,
     &       ncfpl,    iicsplsize, icssharps, iicsdiff,   ilistsd,
     &       iicspltype,iplbasdat, imultipolecut,nomultipole)
      end
c======================================================================
      subroutine para_FTC_INITIAL(
     &       griddens,  irangeacc, irangeacc2, expaccl,   ncspl,
     &       ncfpl,     isharpgrd, isharpness, iicsplsize,icssharps,
     &       iicspltype,iplbasdat, Lxmin,      Lxmax,     Lymin,
     &       Lymax,     Lzmin,     Lzmax,      Lxo,       Lyo,
     &       Lzo,       Lxe,       Lye,        Lze,       PLDmax,
     &       npwx,      npwy,      npwz,       npwxe,     npwye,
     &       npwze,     iinit,   igridranges,igridranges2,maxovs,
     &       isharpovs, insharpovs, ii4sharps, icorecut,
     &    iisharpgridrange1,iisharpgridrange2, griddensf)
      implicit real*8(a-h,o-z)
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
      integer expaccl
      call FTC_INITIAL(
     &       griddens,  irangeacc, irangeacc2, expaccl,   ncspl,
     &       ncfpl,     isharpgrd, isharpness, iicsplsize,icssharps,
     &       iicspltype,iplbasdat, Lxmin,      Lxmax,     Lymin,
     &       Lymax,     Lzmin,     Lzmax,      Lxo,       Lyo,
     &       Lzo,       Lxe,       Lye,        Lze,       PLDmax,
     &       npwx,      npwy,      npwz,       npwxe,     npwye,
     &       npwze,     iinit,   igridranges,igridranges2,maxovs,
     &       isharpovs, insharpovs, ii4sharps, icorecut,
     &    iisharpgridrange1,iisharpgridrange2,griddensf)
      end
c======================================================================
      subroutine para_FTC_ENERGY(
     &     griddens,   rLxmin,     rLymin,     rLzmin,     rLxo,
     &     rLyo,       rLzo,       rLxe,       rLye,       rLze,
     &     PLDmax,     npwx,       npwy,       npwz,       npwxe,
     &     npwye,      npwze,      isharpness, iinit,      ncspl,
     &     igridranges,igridranges2,icssharps, iicsdiff,   iicspltype,
     &     maxovs,     isharpovs,  insharpovs, ii4sharps,  icorecut,
     &iisharpgridrange1,iisharpgridrange2,denA,denB,       fockA,
     &     fockB,      thint,      dist2mp,  imultipolecut,nomultipole,
     &     conver)
      implicit real*8(a-h,o-z)
      call FTC_ENERGY(
     &     griddens,   rLxmin,     rLymin,     rLzmin,     rLxo,
     &     rLyo,       rLzo,       rLxe,       rLye,       rLze,
     &     PLDmax,     npwx,       npwy,       npwz,       npwxe,
     &     npwye,      npwze,      isharpness, iinit,      ncspl,
     &     igridranges,igridranges2,icssharps, iicsdiff,   iicspltype,
     &     maxovs,     isharpovs,  insharpovs, ii4sharps,  icorecut,
     &iisharpgridrange1,iisharpgridrange2,denA,denB,       fockA,
     &     fockB,      thint,      dist2mp,  imultipolecut,nomultipole,
     &     conver)
      return
      end
c======================================================================
      subroutine para_post_scf(ncache,iforwhat,nblocks)
      end
c======================================================================
      subroutine para_grad(idft,ax,rhf,nblocks,bl,inx,ntri,
     *                     thres,dens,denB,atforces,labels)
      implicit real*8 (a-h,o-z)
      call int_grad(idft,ax,rhf,nblocks,bl,inx,ntri,thres,
     *              dens,denB,atforces,labels,0,0)
      end
c======================================================================
      subroutine para_FTC_INIT_FORCES(
     $       isharpgrd, isharpness, ncspl,    ncfpl,     iicsplsize,
     $       icssharps, iicsdiff,   ilistsd,  iicspltype, iplbasdat,
     $       rLxmin,    rLxmax,     rLymin,   rLymax,    rLzmin,
     $       rLzmax,    rLx,        rLy,      rLz,       PLDmax,
     $       npwx,      npwy,       npwz,   igridranges,igridranges2,
     $       maxovs,   isharpovs, insharpovs, ii4sharps, icorecut,
     $     iisharpgridrange1, iisharpgridrange2, griddens)
      implicit real*8 (a-h,o-z)
      call FTC_INIT_FORCES(
     $       isharpgrd, isharpness, ncspl,    ncfpl,     iicsplsize,
     $       icssharps, iicsdiff,   ilistsd,  iicspltype, iplbasdat,
     $       rLxmin,    rLxmax,     rLymin,   rLymax,    rLzmin,
     $       rLzmax,    rLx,        rLy,      rLz,       PLDmax,
     $       npwx,      npwy,       npwz,   igridranges,igridranges2,
     $       maxovs,   isharpovs, insharpovs, ii4sharps, icorecut,
     $     iisharpgridrange1, iisharpgridrange2, griddens)
      return
      end
c======================================================================
      subroutine para_FTC_RHF_FORCES(
     $    natom,      ncs,        ncf,        ibas,       ictr,
     $    griddens,   Lxmin,      Lxmax,      Lymin,      Lymax,
     $    Lzmin,      Lzmax,      Lxo,        Lyo,        Lzo,
     $    Dmax,       npwx,       npwy,       npwz,       isharpgrd,
     $    isharpness, ncspl,   igridranges, igridranges2, ncfpl,
     $    iicsplsize, icssharps,  iicsdiff,   idens,      iyregions,
     $    npwyregdim, iicspltype, maxovs,     isharpovs,  insharpovs,
     $    ii4sharps,  iplbasdat,  icorecut,   iro4,       iro5,
     $  iisharpgridrange1, iisharpgridrange2)
      implicit real*8 (a-h,o-z)
      call FTC_RHF_FORCES(
     $    natom,      ncs,        ncf,        ibas,       ictr,
     $    griddens,   Lxmin,      Lxmax,      Lymin,      Lymax,
     $    Lzmin,      Lzmax,      Lxo,        Lyo,        Lzo,
     $    Dmax,       npwx,       npwy,       npwz,       isharpgrd,
     $    isharpness, ncspl,   igridranges, igridranges2, ncfpl,
     $    iicsplsize, icssharps,  iicsdiff,   idens,      iyregions,
     $    npwyregdim, iicspltype, maxovs,     isharpovs,  insharpovs,
     $    ii4sharps,  iplbasdat,  icorecut,   iro4,       iro5,
     $  iisharpgridrange1, iisharpgridrange2)
      return
      end
c======================================================================
      subroutine para_HessInit(ncache,iforwhat,nblocks)
      end
c======================================================================
c for hessian with first derivatives calculated ONLY    :
c
      subroutine para_d0g1(idft,ax,rhf,nblocks,bl,inx,ntri,thref1,
     *                     natb,nate,listreal,do_allat,
     *                     densp, ncenter,
     *                     dens,denB,fder,fderB,labels)
      implicit real*8 (a-h,o-z)
c
      call int_d0g1(idft,ax,rhf,nblocks,bl,inx,ntri,thref1,
     *                     natb,nate,listreal,do_allat,
     *                     densp, ncenter,
     *              dens,denB,fder,fderB,labels,0,0)
      end
c======================================================================
c for hessian with second derivatives calculated ONLY    :
c
      subroutine para_d0g2(idft,ax,rhf,nblocks,bl,inx,ntri,threg2,
     *                     densp, ncenter,
     *                     dens,denB,hessian,labels)
      implicit real*8 (a-h,o-z)
      call int_d0g2(idft,ax,rhf,nblocks,bl,inx,ntri,threg2,
     *                     densp, ncenter,
     *              dens,denB,hessian,labels,0,0)
      end
c======================================================================
      subroutine para_giao(idft,ax,nblocks,bl,inx,ntri,thres,dscreen,
     *                     dens,fock,labels)
      implicit real*8 (a-h,o-z)
      call int_giao(idft,ax,nblocks,bl,inx,ntri,thres,dscreen,dens,
     *              fock,labels,0,0)
      end
c======================================================================
c for "nmr" CPHF
      subroutine para_cphf(nblocks,bl,inx,ntri,thres,dscreen,dens,
     *                     fock,labels,myjob,igran,mgo)
      implicit real*8 (a-h,o-z)
      call int_cphf(nblocks,bl,inx,ntri,thres,dscreen,dens,fock,
     *              labels,myjob,igran)
      end
c======================================================================
c for "polarizability CPHF
      subroutine para_cphf_xyz(idft,ax,nblocks,bl,inx,ntri,thres,dens,
     *                         fock,labels)
      implicit real*8 (a-h,o-z)
      call int_d1g0(idft,ax,nblocks,bl,inx,ntri,thres,dens,fock,
     *              labels,0,0)
      end
c======================================================================
c for constant part of first-order density and FDS1 in CPHF  (RHF Hessian)
      subroutine para_d1const(natoms, natb,   nate,   listunq,ncf,
     $                        ntri,   nocc,   lind,   dens0,  s1,
     $                        f1,     vec,    val,    bl,     fds1,
     $                        work1)
      use memory, block => bl
      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(*),dens0(*),s1(*),f1(*),vec(*),val(*),
     $          bl(*),fds1(ncf,ncf,3),work1(*)
      data nfile60/60/, nfile61/61/, nfile62/62/, nfile65/65/
      data  xlvsh/0.0d0/
c
c
      natonce=nate-natb+1
      ncf2 = ncf**2
      iat1=0
c
c -- for VCD reserve memory for constant part of perturbed
c -- wavefunction coefficients
      call tstival('vcd',ivcd)
      If(ivcd.NE.0) Then
        call getmem(ncf*nocc*3,lpve)
        ncoefile = 79    ! unit for the pert. coeff. file
      EndIf
c
      do iat=natb,nate
        iat1=iat1+1
        nat=listunq(iat)
c
c  read in fock1 and over1 matrices from a disk
c
        call read1mat(nfile61,nat,ntri*3,f1)
        call read1mat(nfile62,nat,ntri*3,s1)
c
c ...................................................................
c  calculate FDS1 and S1DF matrices and F1-(FDS1+S1DF)
c  use work1 for FD(ncf,ncf)
c
        lrec=1
        call read1mat(nfile60,lrec,ncf2,work1)
c
        call calcF1FDS1(ncf,    ntri,   work1,  s1,  f1,
     *                  fds1)
c
        irec=(iat1-1)*3 + lrec
        call save1mat(nfile60,irec+1,ncf2,fds1(1,1,1) )
        call save1mat(nfile60,irec+2,ncf2,fds1(1,1,2) )
        call save1mat(nfile60,irec+3,ncf2,fds1(1,1,3) )
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(xlvsh.ne.0.d0) then
c              when level shift is to be used then we need an extra term
c              in D1const : -0.25*xlvsh*Ci+*[S0(D0*S1*D0)S0]*Ca .Calculate
c              first the  : -0.25*xlvsh*    [S0(D0*S1*D0)S0]   term and
c              add it to the F1=H1 + G(D0,g1). Then it will be projected
c              in dens1_xyz to the requested form :
c              Use work1 for temporary storage :
c
               call getival('lsmat',ls0)
               call f1shift_xyz(f1,ncf,ntri,lind,dens0,bl(ls0),s1,
     *                          work1, bl,xlvsh)
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  calculate constant part of 1st-order density matrices
c  use work1 to calculate d1const and write it on a disk
c  when calculating D1const save 1/2DS1D in fds1 and put
c  also on a disk after D1const
c
        call d1const_xyz(.true.,   ncf,    ntri,   nocc,   lind,
     *                   vec,      val,    dens0,  f1,     s1,
     *                   work1,    ivcd,  bl(lpve), bl )
c
c  write d1const on a disk
        call save1mat(nfile65,iat1,ntri*3,work1)  ! D1const
        If(ivcd.NE.0)
     *    call save1mat(ncoefile,iat1,ncf*nocc*3,bl(lpve))  !  D1const
      enddo
c
      If(ivcd.NE.0) call retmem(1)
c
      end
c======================================================================
c for constant part of first-order density and FDS1 in CPHF  (UHF Hessian)
c
      subroutine para_d1const_uhf(natoms, natb,   nate,   listunq,ncf,
     $                            ntri,   nalpha, nbeta,  lind,  dens0A,
     $                            dens0B, s1,     f1,     vecA,   vecB,
     $                            valA,   valB,   bl,     fds1,   fds1b,
     $                            work1)
      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(*),dens0A(*),dens0B(*),s1(*),f1(*),
     $          vecA(*),vecB(*),valA(*),valB(*),bl(*),
     $          fds1(ncf,ncf,3),fds1B(ncf,ncf,3),work1(*)
      data nfile60/60/, nfile61/61/, nfile62/62/, nfile65/65/
      data nfile71/71/, nfile75/75/
      data  xlvsh/0.0d0/
c
c
      natonce = nate-natb+1
      ncf2 = ncf**2
      iat1=0
c
      do iat=natb,nate
        iat1=iat1+1
        nat=listunq(iat)
c
c  read in first order overlap matrix for current atom
c
        call read1mat(nfile62,nat,ntri*3,s1)
c
c  alpha  component
c  ----------------
c
c  read in constant part of first order alpha fock matrix
c  for current atom
c
        call read1mat(nfile61,nat,ntri*3,f1)
c
c ...................................................................
c  calculate FDS1 and S1DF matrices and F1-(FDS1+S1DF)
c  use work1 for FD(ncf,ncf)
c
        lrec=1
        call read1mat(nfile60,lrec,ncf2,work1)
        call dscal(ncf2, 2.d0 ,work1 ,1)
c
        call calcF1FDS1(ncf,    ntri,   work1,  s1,  f1,
     *                  fds1)
        irec=(iat1-1)*3 + lrec
        call save1mat(nfile60,irec+1,ncf2,fds1(1,1,1) )
        call save1mat(nfile60,irec+2,ncf2,fds1(1,1,2) )
        call save1mat(nfile60,irec+3,ncf2,fds1(1,1,3) )
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(xlvsh.ne.0.d0) then
c              when level shift is to be used then we need an extra term
c              in D1const : -0.25*xlvsh*Ci+*[S0(D0*S1*D0)S0]*Ca .Calculate
c              first the  : -0.25*xlvsh*    [S0(D0*S1*D0)S0]   term and
c              add it to the F1=H1 + G(D0,g1). Then it will be projected
c              in dens1_xyz to the requested form :
c              Use work1 for temporary storage :
c
               call getival('lsmat',ls0)
               call f1shift_xyz(f1,ncf,ntri,lind,dens0A,bl(ls0),s1,
     *                          work1, bl,xlvsh)
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  calculate constant part of 1st-order alpha density mat.
c  use work1 to calculate d1const and write it on a disk
c  work1 does not need to be initialized
c
        call d1const_xyz(.false., ncf,    ntri,  nalpha, lind,
     *                   vecA,   valA,   dens0A, f1,     s1,
     *                   work1,   0,      bl,    bl) ! no VCD here
c
c  write alpha d1const on a disk
        call save1mat(nfile65,iat1,ntri*3,work1)  ! D1const
c
c  beta  component:
c  ----------------
c
c  read in constant part of first order beta fock matrix
c  for current atom
c
        call read1mat(nfile71,nat,ntri*3,f1)
c
c ...................................................................
c  calculate FDS1 and S1DF matrices and F1-(FDS1+S1DF)
c  use work1 for FD(ncf,ncf)
        lrec=natonce*3+2
        call read1mat(nfile60,lrec,ncf2,work1)
        call dscal(ncf2, 2.d0 ,work1 ,1)
        call calcF1FDS1(ncf,    ntri,   work1,  s1,  f1,
     *                  fds1)
        irec=(iat1-1)*3 + lrec
        call save1mat(nfile60,irec+1,ncf2,fds1(1,1,1) )
        call save1mat(nfile60,irec+2,ncf2,fds1(1,1,2) )
        call save1mat(nfile60,irec+3,ncf2,fds1(1,1,3) )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(xlvsh.ne.0.d0) then
c              when level shift is to be used then we need an extra term
c              in D1const : -0.25*xlvsh*Ci+*[S0(D0*S1*D0)S0]*Ca .Calculate
c              first the  : -0.25*xlvsh*    [S0(D0*S1*D0)S0]   term and
c              add it to the F1=H1 + G(D0,g1). Then it will be projected
c              in dens1_xyz to the requested form :
c              Use work1 for temporary storage :
c
               call f1shift_xyz(f1,ncf,ntri,lind,dens0B,bl(ls0),s1,
     *                          work1, bl,xlvsh)
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  calculate constant part of 1st-order beta density mat.
c  use work1 to calculate d1const and write it on a disk
c  work1 does not need to be initialized
c
        call d1const_xyz(.false., ncf,    ntri,  nbeta,  lind,
     *                   vecB,   valB,   dens0B, f1,     s1,
     *                   work1,   0,      bl,    bl) ! no VCD here
c
c  write beta d1const on a disk
        call save1mat(nfile75,iat1,ntri*3,work1)  ! D1const
      enddo
c
      return
      end
c======================================================================
c for weighted first-order density in CPHF (RHF Hessian)
      subroutine para_wdens(natb,   nate,   listunq,ncf,    ntri,
     $                      lind,   dens0,  fd0,    f1,     gmat1,
     $                      d1,     work,   r1)
      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(*),dens0(*),fd0(*),work(ncf,ncf,3),
     $          f1(*),gmat1(ntri,3,*),d1(ntri,3,*),r1(ntri*3)
      data nfile61/61/, nfile63/63/, nfile64/64/

      iat1=0
      do iat=natb,nate
        iat1=iat1+1
        nat=listunq(iat)
c -- read in fock1 matrix from disk
        call read1mat(nfile61,nat,ntri*3,f1)
c
c  calculate W1 1st-order weighted density
c
        call wdens1(.true., ncf,    ntri,   lind,   dens0,
     *              fd0,    f1, gmat1(1,1,iat1),d1(1,1,iat1),
     *              work,   r1)
c  on return r1 contains final 1st-order weighted density
c
c  write density & weighted density to disk
c
        call save1mat(nfile63,nat,ntri*3,d1(1,1,iat1))
        call save1mat(nfile64,nat,ntri*3,r1)
      enddo
      return
      end
c======================================================================
c for weighted first-order density in CPHF (UHF Hessian)
      subroutine para_wdensu_xyz(natb,   nate,   listunq,ncf,    ntri,
     $                           lind,   dens0A, dens0B, fd0A,   fd0B,
     $                           f1,     gmat1A, gmat1B, d1A,    d1B,
     $                           work,   s1,     r1)
      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(*),dens0A(*),dens0B(*),fd0A(*),
     $          fd0B(*),f1(ntri,3),gmat1A(ntri,3,*),gmat1B(ntri,3,*),
     $          d1A(ntri,3,*),d1B(ntri,3,*),work(ncf,ncf,3),
     $          s1(ntri*3),r1(ntri*3)
      data nfile61/61/, nfile63/63/, nfile64/64/,
     $     nfile71/71/, nfile73/73/

      iat1=0
      do iat=natb,nate
        iat1=iat1+1
        nat=listunq(iat)
c -- read in alpha fock1 matrix from disk
        call read1mat(nfile61,nat,ntri*3,f1)
c
c  calculate W1 alpha 1st-order weighted density
c
        call wdens1(.false.,ncf,    ntri,   lind,   dens0A,
     *              fd0A,   f1,  gmat1A(1,1,iat1),d1A(1,1,iat1),
     *              work,   r1)
c  on return r1 contains final alpha 1st-order weighted density
c
c -- read in beta fock1 matrix from disk
        call read1mat(nfile71,nat,ntri*3,f1)
c
c  calculate W1 beta 1st-order weighted density
c
        call wdens1(.false.,ncf,    ntri,   lind,   dens0B,
     *              fd0B,   f1,  gmat1B(1,1,iat1),d1B(1,1,iat1),
     *              work,   s1)
c  on return s1 contains final beta 1st-order weighted density
c
c -- sum alpha and beta components
        call AddVec(ntri*3,r1,s1,r1)
c
c  write density & weighted density to disk
c
        call save1mat(nfile63,nat,ntri*3,d1A(1,1,iat1))
        call save1mat(nfile73,nat,ntri*3,d1B(1,1,iat1))
        call save1mat(nfile64,nat,ntri*3,r1)
      enddo
      return
      end
c======================================================================
      subroutine para_dft(dft,    NAtoms, XNuc,   IAN,    NSym,
     $                    NGen,   ISYM,   NEqATM, NQ,     IUNQ,
     $                    IPRNT,  NRad,   NAng,   IradQ,  factor,
     $                    NBatch, lgrid,  lsemi,  DISTN,  AIJ,
     $                    RDIST,  XXA,    WTA,    XGRID,  WGHT,
     $                    thrsh,  NBas,   NShell, BASDAT, INX,
     $                    BL,     ExpMIN, rhf,    DA,     DB,
     $                    DM,     NSLAVE, LSlave, LFINI,  NScr,
     $                    Z,      XCA,    XCB,    EXC,    EL)
      implicit real*8(a-h,o-z)
c
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),IUNQ(NQ),ISYM(NGen),
     $          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     $          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NBas),
     $          XCA(NBas*(NBas+1)/2),XCB(NBas*(NBas+1)/2)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     $          XXA(3,1130),WTA(1130)
      Dimension NSLAVE(NQ),LFINI(NQ)
      DIMENSION Z(NScr)
      INTEGER dft
      Logical rhf,LSlave(NQ)
c
c
c -- SINGLE PROCESSOR MODE
c
      DO 100 IAtom=1,NQ
      If(IAtom.EQ.1) Then
        IEntry = 1
      Else If(IAtom.EQ.NQ) Then
        IEntry = -1
      EndIf
      ICntr = IUNQ(IAtom)
        If(rhf) Then
          CALL DFTFockC(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                  NSym,   NGen,   ISYM,   NEqATM, IPRNT,
     $                  NRad,   NAng,   IradQ,  factor, NBatch,
     $                  DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                  XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                  BASDAT, INX,    BL,     ExpMIN, DA,
     $                  DM,     NScr,   Z,      XCA,    EXC,
     $                  EL,     lgrid,  lsemi,  IEntry)
        Else
          CALL DFTFockU(dft,    ICntr,  NAtoms, XNuc,   IAN,
     $                  NSym,   NGen,   ISYM,   NEqATM, IPRNT,
     $                  NRad,   NAng,   IradQ,  factor, NBatch,
     $                  DISTN,  AIJ,    RDIST,  XXA,    WTA,
     $                  XGRID,  WGHT,   thrsh,  NBas,   NShell,
     $                  BASDAT, INX,    BL,     ExpMIN, DA,
     $                  DB,     DM,     NScr,   Z,      XCA,
     $                  XCB,    EXC,    EL,     lgrid,  lsemi,
     $                  IEntry)
        EndIf
 100  CONTINUE
c
      end
c======================================================================
      subroutine para_dftg(dft,    NAtoms, XNuc,   IAN,    nsym,
     $                     ngen,   ISYM,   NEqATM, NEqBAS, NQ,
     $                     IUNQ,   IPRNT,  lrad,   lang,   IradQ,
     $                     factor, NBatch, DISTN,  AIJ,    RDIST,
     $                     XXA,    WTA,    XGRID,  WGHT,   thrsh,
     $                     NBas,   NShell, BASDAT, INX,    BL,
     $                     ExpMIN, NBAtm,  IdWt,   rhf,    DA,
     $                     DB,     NScr,   Z,      GX,     GY,
     $                     GZ,     GXB,    GYB,    GZB,    GWT,
     $                     GC,     EL)
      end
c======================================================================
      SUBROUTINE para_dfth(dft,    NAtoms, ipass,  npass,  Nb,
     $                     Ne,     XNuc,   IAN,    nsym,   ngen,
     $                     ISYM,   NEqATM, NQ,     IUNQ,   IPRNT,
     $                     lrad,   lang,   IradQ,  factor, NBatch,
     $                     thrsh,  NBas,   NShell, BASDAT, INX,
     $                     NBAtm,  IdWt,   rhf,    DA,     DB,
     $                     DM,     FDA,    FDB,    HESS,   EL)
      end
c======================================================================
      subroutine para_dftcphf(dft,    NAtoms, NatOnce,NatDo,  XNuc,
     $                        IAN,    NQ,     IUNQ,   XSym,   lrad,
     $                        lang,   IradQ,  factor, NBatch, thrsh,
     $                        NBas,   NShell, BASDAT, INX,    IdWt,
     $                        rhf,    DA,     DB,     DEN1,   DEN1B,
     $                        DM1,    FDA,    FDB,    lgrid)
      end
c======================================================================
      subroutine para_dftgiao(dft,    NAtoms, XNuc,   IAN,    NSym,
     *                        NGen,   ISYM,   NEqATM, NQ,     IUNQ,
     *                        IPRNT,  NRad,   NAng,   IradQ,  factor,
     *                        NBatch, DISTN,  AIJ,    RDIST,  XXA,
     *                        WTA,    XGRID,  WGHT,   thrsh,  NBas,
     *                        NShell, BASDAT, INX,    BL,     ExpMIN,
     *                        NOcc,   CMO,    VCMO,   Malk,
     *                        NScr,   Z,      XC,     XMalkin)
      implicit real*8(a-h,o-z)
c
      DIMENSION XNuc(3,NAtoms),IAN(NAtoms),IUNQ(NQ),ISYM(NGen),
     *          XGRID(3,*),WGHT(*),ExpMIN(NShell),BL(*),
     *          NEqATM(NAtoms,NSym),BASDAT(13,*),INX(12,*),
     *          CMO(NOcc,NBas),
     *          XC(3*NBas*(NBas+1)/2),
     $          VCMO(NBas-NOcc,NBas),XMalkin(NBas-NOcc,NOcc)
      DIMENSION DISTN(NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms),
     *          XXA(3,1130),WTA(1130)
      DIMENSION Z(NScr)
      INTEGER dft
c
c
c -- SINGLE PROCESSOR MODE
c
      DO 100 IAtom=1,NQ
      If(IAtom.EQ.1) Then
        IEntry = 1
      Else If(IAtom.EQ.NQ) Then
        IEntry = -1
      EndIf
      ICntr = IUNQ(IAtom)
      CALL DFTNMRC(dft,    ICntr,  NAtoms, XNuc,   IAN,
     *             NSym,   NGen,   ISYM,   NEqATM, IPRNT,
     *             NRad,   NAng,   IradQ,  factor, NBatch,
     *             DISTN,  AIJ,    RDIST,  XXA,    WTA,
     *             XGRID,  WGHT,   thrsh,  NBas,   NShell,
     *             BASDAT, INX,    BL,     ExpMIN, NOcc,
     *             CMO,    VCMO,   Malk,   NScr,   Z,
     *             XC,     XMalkin,IEntry)
 100  CONTINUE
      END
c=======================================================================
      subroutine para_shift_start(bl,lden,ldn1,ntri)
      implicit real*8 (a-h,o-z)
      end
c=======================================================================
      subroutine para_shift_end(shie,nreq)
      implicit real*8 (a-h,o-z)
      end
c=======================================================================
      subroutine para_shi(na,bl,inx,last,nprint,iout,blibas,blinuc,
     *            ncs,ncf,ntri,nreq,den,dn1,shie,hna,nra,nrat,ldn1,lhna)
      implicit real*8 (a-h,o-z)
      call shield(na,bl,inx,last,nprint,iout,blibas,blinuc,
     *            ncs,ncf,ntri,nreq,den,dn1,shie,hna,nra,nrat,ldn1,lhna)
      end
c=======================================================================
      SUBROUTINE Para_PROPMAIN(NAtoms, IAN,    QA,     XC,     NBas,
     $                         NShell, nsh,    BASDAT, INX,    factor,
     $                         r0f,    LMax,   NRad,   NAng,   ExpMIN,
     $                         BL,     AIJ,    DISTN,  RDIST,  radii,
     $                         radwght,XXA,    WTA,    NAlpha, NBeta,
     $                         rhf,    Ylm,    CMO,    CA,     CMOB,
     $                         CB,     NTrans, IUNQ,   TRANS,  NEqATM,
     $                         SIGNA,  SIGNB,  XCharg, XSpin,  XEFG,
     $                         NScr,   Z,      IErr)
      implicit real*8(a-h,o-z)
      CALL PROPMAIN(NAtoms, IAN,    QA,     XC,     NBas,
     $              NShell, nsh,    BASDAT, INX,    factor,
     $              r0f,    LMax,   NRad,   NAng,   ExpMIN,
     $              BL,     AIJ,    DISTN,  RDIST,  radii,
     $              radwght,XXA,    WTA,    NAlpha, NBeta,
     $              rhf,    Ylm,    CMO,    CA,     CMOB,
     $              CB,     NTrans, IUNQ,   TRANS,  NEqATM,
     $              SIGNA,  SIGNB,  XCharg, XSpin,  XEFG,
     $              NScr,   Z,      IErr)
      end
c=======================================================================
      subroutine para_rmp2(nmo)
      implicit real*8(a-h,o-z)
      call rmp2(nmo)
      end
c=======================================================================
      subroutine para_ump2(nalpha,nbeta)
      implicit real*8(a-h,o-z)
      call ump2(nalpha,nbeta)
      end
c=======================================================================
      subroutine para_cosmo_surfrep(nps)
      implicit real*8(a-h,o-z)
      end
c=======================================================================
      subroutine para_cosmo_pot(den,phin,phi,nps,ntri)
      implicit real*8(a-h,o-z)
      end
c=======================================================================
      subroutine para_cosmo_h0(h0cos,qcos,fepsi,nps,ntri)
      implicit real*8(a-h,o-z)
      end
c======================================================================
      subroutine para_cosmo_data(cosurf,xyz,natom,nps,npspher)
      implicit real*8 (a-h,o-z)
      end
c======================================================================
      subroutine para_cosmo_forc(den,cfor,qcos,cosurf,iatsp,
     $                           na,nps,ntri)
      implicit real*8 (a-h,o-z)
      end
c======================================================================
      subroutine para_printhelp(ichan)
c
c  prints a help message (serial version)
c
      implicit none
      integer ichan
      write(ichan,100)
 100  format(/' USAGE: pqs.x [-v,-h] [-i] molecule-name'/)
      end
c======================================================================
cc      subroutine fafOpenm(filename,ifileID)
cc      character*(*) filename
cc      end
c======================================================================
      subroutine fafCreatem(filename,ifileID)
      character*(*) filename
      end
c======================================================================
      subroutine fafCreates(filename,ifileID)
      character*(*) filename
      end
c======================================================================
      subroutine fafwrite(array,length,num,filename,n,ir)
      end
c======================================================================
      subroutine fafread(array,length,num,filename,n,ir)
      end
c======================================================================
      subroutine fafClosem(ifileID,keep,info)
      end
c======================================================================
      subroutine fafCloses(ifileID,keep,info)
      end
c======================================================================
      subroutine fafrequest
      end
c======================================================================
      subroutine fafnbread
      end
c======================================================================
      subroutine afwritebin(ndisk,lbin,ibin4,ibin1,indxbin,irecord,ij,
     1                      npairs,afname,islvid)
      end
c======================================================================
      subroutine afwritebingr(ndisk,lbin,ibin4,ibin1,indxbin,irecord,ij,
     1                      npairs,afname,islvid,isize)
      end
c======================================================================
cc      subroutine para_coulomb
cc      call coulomb
cc      end
c======================================================================
      subroutine par_diis_matrix()
      end
c======================================================================
      subroutine par_new_coeff()
      end
c======================================================================
      subroutine diis_slave_next()
      end
c======================================================================
      subroutine coupled_cluster()
      call coulomb
      end
c======================================================================
      subroutine array_files()
      end
c======================================================================
      subroutine para_cimsubgen(nclu)
      integer nclu
      call cimsubgen(nclu)
      end
c======================================================================
      subroutine para_cimsub
      call cimsub
      end
c======================================================================
