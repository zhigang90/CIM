      subroutine uhf_bldr(ncf,    buf,    fockA,  fockB,  DenA,
     *                    DenB,   DenAB,  iarray, nbls,   ngcd,
     *                    lnijkl, labels, length, lgenct)
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c  input parameters
c  ncf    -  number of contracted basis functions
c  buf    -  integrals buffer
c  fockA  -  alpha Fock matrix
c  fockB  -  beta Fock matrix
c  DenA   -  alpha density matrix
c  DenB   -  beta density matrix
c  DenAB  -  alpha+beta density matrix
c  iarray -  the triangular numbers, iarray(i,j)=i*(i-1)/2+j
c  nbls   -  integ. block-size , number of contracted shell quartets
c----------------------------------------------------------------
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fockA(*),fockB(*),DenA(*),DenB(*),DenAB(*)
      dimension iarray(ncf,ncf)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      parameter (half=0.5d0)
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
            icff=labels(1,iqu,ijklp)
            jcff=labels(2,iqu,ijklp)
            kcff=labels(3,iqu,ijklp)
            lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
            intct=0
             icf=icff
             do 200 iii=1,ilen
             icf=icf+1
             jcf=jcff
             do 200 jjj=1,jlen
             jcf=jcf+1
             iijf=iarray(jcf,icf)
             kcf=kcff
             do 200 kkk=1,klen
             kcf=kcf+1
             iikf=iarray(kcf,icf)
             jjkf=iarray(kcf,jcf)
             lcf=lcff
             do 200 lll=1,llen
             lcf=lcf+1
             iilf=iarray(lcf,icf)
             jjlf=iarray(lcf,jcf)
             kklf=iarray(lcf,kcf)
             intct=intct+1
             xint=buf(ijklp,intct,iqu)
ctest
c               write(8,55) icf,jcf,kcf,lcf,xint*1.0d-10
c    *
c55              format(4i3,3f15.9)
ctest
c------------------------------------------------------
             if(abs(xint).gt.half) then
                xin2=xint+xint
c
          fockA(iijf)=fockA(iijf)+xin2*DenAB(kklf)
          fockA(kklf)=fockA(kklf)+xin2*DenAB(iijf)
          fockA(iikf)=fockA(iikf)-xint*DenA(jjlf)
          fockA(jjlf)=fockA(jjlf)-xint*DenA(iikf)
          fockA(jjkf)=fockA(jjkf)-xint*DenA(iilf)
          fockA(iilf)=fockA(iilf)-xint*DenA(jjkf)
c
          fockB(iijf)=fockB(iijf)+xin2*DenAB(kklf)
          fockB(kklf)=fockB(kklf)+xin2*DenAB(iijf)
          fockB(iikf)=fockB(iikf)-xint*DenB(jjlf)
          fockB(jjlf)=fockB(jjlf)-xint*DenB(iikf)
          fockB(jjkf)=fockB(jjkf)-xint*DenB(iilf)
          fockB(iilf)=fockB(iilf)-xint*DenB(jjkf)
c
             endif
c------------------------------------------------------
  200        continue
  150   continue
  100 continue
c
      end
c==============================================================
      subroutine uhf_bldr1(ncf,    buf,    fockA,  fockB,  DenA,
     *                     DenB,   DenAB,  iarray, nbls,   ngcd,
     *                     lnijkl, labels, length, lgenct)
      implicit real*8 (a-h,o-z)
c
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fockA(*),fockB(*),DenA(*),DenB(*),DenAB(*)
      dimension iarray(ncf,ncf)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      equivalence(iix(1),icf),(iix(2),jcf),(iix(3),kcf),(iix(4),lcf)
      dimension ishell(1296),jshell(1296),kshell(1296),lshell(1296)
      parameter (half=0.5d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
      len=ilen*jlen*klen*llen
c---------------------------------------------------------------
c  establish a simple loop for the short inner loops
      if(len.le.1296) then
        intct=0
        do 350 iii=1,ilen
          do 350 jjj=1,jlen
            do 350 kkk=1,klen
              do 350 lll=1,llen
                intct=intct+1
                ishell(intct)=iii
                jshell(intct)=jjj
                kshell(intct)=kkk
                lshell(intct)=lll
 350     continue
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 400 ijklp=1,nbls
        ngcq=lgenct(ijklp)
          do 450 iqu=1,ngcq
            icff=labels(1,iqu,ijklp)
            jcff=labels(2,iqu,ijklp)
            kcff=labels(3,iqu,ijklp)
            lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
             do 500 intct=1,len
             xint=buf(ijklp,intct,iqu)
             if(abs(xint).gt.half) then
ctest
c               write(8,55) icf,jcf,kcf,lcf,xint*1.0d-10
c    *
c55              format(4i3,3f15.9)
ctest
c------------------------------------------------------
              icf=icff+ishell(intct)
              jcf=jcff+jshell(intct)
              kcf=kcff+kshell(intct)
              lcf=lcff+lshell(intct)
                xin2=xint+xint
c
                iijf=iarray(jcf,icf)
                kklf=iarray(lcf,kcf)
                iikf=iarray(kcf,icf)
                jjlf=iarray(lcf,jcf)
                iilf=iarray(lcf,icf)
                jjkf=iarray(kcf,jcf)
c
          fockA(iijf)=fockA(iijf)+xin2*DenAB(kklf)
          fockA(kklf)=fockA(kklf)+xin2*DenAB(iijf)
          fockA(iikf)=fockA(iikf)-xint*DenA(jjlf)
          fockA(jjlf)=fockA(jjlf)-xint*DenA(iikf)
          fockA(jjkf)=fockA(jjkf)-xint*DenA(iilf)
          fockA(iilf)=fockA(iilf)-xint*DenA(jjkf)
c
          fockB(iijf)=fockB(iijf)+xin2*DenAB(kklf)
          fockB(kklf)=fockB(kklf)+xin2*DenAB(iijf)
          fockB(iikf)=fockB(iikf)-xint*DenB(jjlf)
          fockB(jjlf)=fockB(jjlf)-xint*DenB(iikf)
          fockB(jjkf)=fockB(jjkf)-xint*DenB(iilf)
          fockB(iilf)=fockB(iilf)-xint*DenB(jjkf)
c------------------------------------------------------
              end if
  500        continue
  450   continue
  400 continue
      else
        print *,'The program should not have gotten here!!!'
      end if
c
      end
c==============================================================
      subroutine uhf_bldr4(ncf,    buf,    fockA,  fockB,  DenA,
     *                     DenB,   DenAB,  iarray, nbls,   ngcd,
     *                     lnijkl, labels, length, lgenct)
      implicit real*8 (a-h,o-z)
c
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fockA(*),fockB(*),DenA(*),DenB(*),DenAB(*)
      dimension iarray(ncf,ncf)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      equivalence(iix(1),icf),(iix(2),jcf),(iix(3),kcf),(iix(4),lcf)
      dimension ishell(1296),jshell(1296),kshell(1296),lshell(1296)
      parameter (half=0.5d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      len=ilen*jlen*klen
c---------------------------------------------------------------
c  establish a simple loop for the short inner loops
      if(len.le.1296) then
        intct=0
        do 50 iii=1,ilen
          do 50 jjj=1,jlen
            do 50 kkk=1,klen
                intct=intct+1
                ishell(intct)=iii
                jshell(intct)=jjj
                kshell(intct)=kkk
 50     continue
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
        ngcq=lgenct(ijklp)
          do 150 iqu=1,ngcq
            icff=labels(1,iqu,ijklp)
            jcff=labels(2,iqu,ijklp)
            kcff=labels(3,iqu,ijklp)
            lcf=labels(4,iqu,ijklp)+1
c-------------------------------------------------------
             do 200 intct=1,len
             xint=buf(ijklp,intct,iqu)
             if(abs(xint).gt.half) then
ctest
c               write(8,55) icf,jcf,kcf,lcf,xint*1.0d-10
c    *
c55              format(4i3,3f15.9)
ctest
c------------------------------------------------------
              icf=icff+ishell(intct)
              jcf=jcff+jshell(intct)
              kcf=kcff+kshell(intct)
c
                iijf=iarray(jcf,icf)
                kklf=iarray(lcf,kcf)
                iikf=iarray(kcf,icf)
                iilf=iarray(lcf,icf)
                jjkf=iarray(kcf,jcf)
                jjlf=iarray(lcf,jcf)
c
                xin2=xint+xint
c
          fockA(iijf)=fockA(iijf)+xin2*DenAB(kklf)
          fockA(kklf)=fockA(kklf)+xin2*DenAB(iijf)
          fockA(iikf)=fockA(iikf)-xint*DenA(jjlf)
          fockA(jjlf)=fockA(jjlf)-xint*DenA(iikf)
          fockA(jjkf)=fockA(jjkf)-xint*DenA(iilf)
          fockA(iilf)=fockA(iilf)-xint*DenA(jjkf)
c
          fockB(iijf)=fockB(iijf)+xin2*DenAB(kklf)
          fockB(kklf)=fockB(kklf)+xin2*DenAB(iijf)
          fockB(iikf)=fockB(iikf)-xint*DenB(jjlf)
          fockB(jjlf)=fockB(jjlf)-xint*DenB(iikf)
          fockB(jjkf)=fockB(jjkf)-xint*DenB(iilf)
          fockB(iilf)=fockB(iilf)-xint*DenB(jjkf)
c------------------------------------------------------
              end if
  200        continue
  150   continue
  100 continue
      else
        print *,'The program should not have gotten here!!!'
      end if
c
      end
c==============================================================
      subroutine uhf_core(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        DenA,        DenB,       DenAB,      fockA,
     $        fockB)
c----------------------------------------------------------------
c Fock builder (using in-core integrals )   UHF
c----------------------------------------------------------------
c input :
c core(*) - array with stored integrals
c i-llab()- starting indecies for integrals
c i-lsiz()- shell's sizes for each super-block
c iqrt()  - number of stored quartets for each super-block
c integrals - number of stored integrals
c nqstore   - number of stored quartets
c nblstore  - number of stored big-blocks
c
c  DenA   -  alpha density matrix
c  DenB   -  beta density matrix
c  DenAB  -  alpha+beta density matrix
c  fockA  -  alpha Fock matrix
c  fockB  -  beta Fock matrix
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
cccc  integer*2 isiz(nqstore),jsiz(nqstore),ksiz(nqstore),lsiz(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      dimension core(*)          ! dimension core(integrals)
      dimension DenA(*),DenB(*),DenAB(*),fockA(*),fockB(*)
      dimension lind(*)
      dimension iix(4)
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
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
c
           do iii=1,ilen
              icf=icff+iii
              do jjj=1,jlen
                 jcf=jcff+jjj
                 do kkk=1,klen
                    kcf=kcff+kkk
                    do lll=1,llen
                       intx=intx+1
                       xint=core(intx)
                       xin2=xint+xint
                       lcf=lcff+lll
canonical   order
                       call descend(icf,jcf,kcf,lcf,iix)
c
                       iff=iix(1)
                       jff=iix(2)
                       kff=iix(3)
                       lff=iix(4)
c
ctest.......................................................
c                 write(6,55) icf,jcf,kcf,lcf,xint*1.0d-10
c55                format(4i3,f15.9,2x,'fock-core')
ctest.......................................................
c
                       ii=lind(iff)
                       jj=lind(jff)
                       kk=lind(kff)
c
                       ij=ii+jff
                       kl=kk+lff
                       ik=ii+kff
                       il=ii+lff
c
                       if(jff.ge.kff) then
                         jk=jj+kff
                       else
                         jk=kk+jff
                       end if
                       if(jff.ge.lff) then
                         jl=jj+lff
                       else
                         ll=lind(lff)
                         jl=ll+jff
                       end if
c
                       fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                       fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                       fockA(ik)=fockA(ik)-xint*DenA(jl)
                       fockA(jl)=fockA(jl)-xint*DenA(ik)
                       fockA(jk)=fockA(jk)-xint*DenA(il)
                       fockA(il)=fockA(il)-xint*DenA(jk)
c
                       fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                       fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                       fockB(ik)=fockB(ik)-xint*DenB(jl)
                       fockB(jl)=fockB(jl)-xint*DenB(ik)
                       fockB(jk)=fockB(jk)-xint*DenB(il)
                       fockB(il)=fockB(il)-xint*DenB(jk)
c
                    enddo         !   over lll=1,llen
                 enddo            !   over kkk=1,klen
              enddo               !   over jjj=1,jlen
           enddo                  !   over iii=1,ilen
         enddo                    !   do iqrt1=1,nqrt1
      enddo                       !   do isbl=1,nblstore
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(ncf*(ncf+1)/2,resc,focka,1)
        call dscal(ncf*(ncf+1)/2,resc,fockb,1)
      endif
c
ctest
          if(integrals.NE.intx) then
             call nerror(1,'fock_core',
     $                   'integrals NE intx',integrals,intx)
          endif
c
ctest
c
      end
c=====================================================================
      subroutine uhf_disk(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dena,       denb,
     $        denab,       focka,      fockb,      ncs,
     $        densp,       map_fs)
c----------------------------------------------------------------
c Fock builder (using on-disk integrals )   UHF
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
      dimension DenA(*),DenB(*),DenAB(*),fockA(*),fockB(*)
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
      dmax1=1.0d0
      dmax=dmax1*t01
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
        endif
c
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
c--------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
c--------------------------------------------------------
            dmax2=half/dmax
c
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_fock8_uhf(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          dena,        denb,       denab,      focka,
     $          fockb,       lind,       dmax2,      intx8,
     $          xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_fock4_uhf(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          dena,        denb,       denab,      focka,
     $          fockb,       lind,       dmax2,      intx4,
     $          xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_fock2_uhf(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          dena,        denb,       denab,      focka,
     $          fockb,       lind,       dmax2,      intx2,
     $          xinteg2)
            endif
c
            if(doint62) then
             if(integi62.gt.int62 .and. icase.eq.4) then
               if(dmax1*1.0d-6.gt.thres1) then
                 call make_fock2_uhf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dena,        denb,       denab,      focka,
     $             fockb,       lind,       dmax2,      int62,
     $             xinte62)
               endif
             endif
            endif
c
            if(doint72) then
             if(integi72.gt.int72 .and. icase.eq.5) then
               if(dmax1*1.0d-7.gt.thres1) then
                 call make_fock2_uhf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dena,        denb,       denab,      focka,
     $             fockb,       lind,       dmax2,      int72,
     $             xinte72)
               endif
             endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                     !   do isbl=1,nblstore
      enddo                     !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(ncf*(ncf+1)/2,resc,focka,1)
        call dscal(ncf*(ncf+1)/2,resc,fockb,1)
      endif
c
      end
c=====================================================================
      subroutine acmuhf_disk(
     $        ax,          ncf,        xinte72,    xinte62,
     $        xinteg2,     xinteg4,    xinteg8,    ijkllab,
     $        ijklsiz,     nquarts,    iqstore,    iblstore,
     $        rescale,     nsplit,     ifromto,    nofinteg,
     $        thres0,      thres1,     lind,       dena,
     $        denb,        denab,      focka,      fockb,
     $        ncs,         densp,      map_fs)
c----------------------------------------------------------------
c Fock builder (using on-disk integrals )   UHF
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
      dimension DenA(*),DenB(*),DenAB(*),fockA(*),fockB(*)
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
      dmax1=1.0d0
      dmax=dmax1*t01
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
        endif
c
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
c--------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
c--------------------------------------------------------
            dmax2=half/dmax
c
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_acmfock8_uhf(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          dena,        denb,       denab,      focka,
     $          fockb,       lind,       dmax2,      intx8,
     $          xinteg8,     ax)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_acmfock4_uhf(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          dena,        denb,       denab,      focka,
     $          fockb,       lind,       dmax2,      intx4,
     $          xinteg4,     ax)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_acmfock2_uhf(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          dena,        denb,       denab,      focka,
     $          fockb,       lind,       dmax2,      intx2,
     $          xinteg2,     ax)
            endif
c
            if(doint62) then
             if(integi62.gt.int62 .and. icase.eq.4) then
               if(dmax1*1.0d-6.gt.thres1) then
                 call make_acmfock2_uhf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dena,        denb,       denab,      focka,
     $             fockb,       lind,       dmax2,      int62,
     $             xinte62,     ax)
               endif
             endif
            endif
c
            if(doint72) then
             if(integi72.gt.int72 .and. icase.eq.5) then
               if(dmax1*1.0d-7.gt.thres1) then
                 call make_acmfock2_uhf(
     $             ilen,        jlen,       klen,       llen,
     $             icff,        jcff,       kcff,       lcff,
     $             dena,        denb,       denab,      focka,
     $             fockb,       lind,       dmax2,      int72,
     $             xinte72,     ax)
               endif
             endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                     !   do isbl=1,nblstore
      enddo                     !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(ncf*(ncf+1)/2,resc,focka,1)
        call dscal(ncf*(ncf+1)/2,resc,fockb,1)
      endif
c
      end
c=====================================================================
      subroutine dftuhf_disk(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       denAB,      fockA,
     $        fockB,       ncs,        densp,      map_fs)
c----------------------------------------------------------------
c input :
c core(*) - array with stored integrals
c i-llab()- starting indecies for integrals
c i-lsiz()- shell's sizes for each super-block
c nquarts()  - number of stored quartets for each super-block
c integrals - number of stored integrals
c nqstore   - number of stored quartets
c nblstore  - number of stored big-blocks
c----------------------------------------------------------------
c Fock builder (using on-disk integrals )
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
      dimension DenAB(*),fockA(*),fockB(*)
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
c if screenig density densp() is not used then 
c
      dmax1=1.0d0
      dmax=dmax1*t01
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
        endif
c
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
c
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
c------------------------------------------------------
            dmax2=half/dmax
c
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_dftuhf8(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          denab,       focka,      fockb,      lind,
     $          dmax,        intx8,      xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_dftuhf4(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          denab,       focka,      fockb,      lind,
     $          dmax,        intx4,      xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_dftuhf2(
     $          ilen,        jlen,       klen,       llen,
     $          icff,        jcff,       kcff,       lcff,
     $          denab,       focka,      fockb,      lind,
     $          dmax,        intx2,      xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_dftuhf2(
     $              ilen,        jlen,       klen,       llen,
     $              icff,        jcff,       kcff,       lcff,
     $              denab,       focka,      fockb,      lind,
     $              dmax,        int62,      xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                  call make_dftuhf2(
     $              ilen,        jlen,       klen,       llen,
     $              icff,        jcff,       kcff,       lcff,
     $              denab,       focka,      fockb,      lind,
     $              dmax,        int72,      xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                  !   do isbl=1,nblstore
      enddo                  !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(ncf*(ncf+1)/2,resc,focka,1)
        call dscal(ncf*(ncf+1)/2,resc,fockb,1)
      endif
c
      end
c====================================================================
      subroutine make_fock8_uhf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dena,        denb,       denab,      focka,
     $        fockb,       lind,       dmax,       intx8,
     $        xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dena(*),denb(*),denab(*),focka(*),fockb(*)
      dimension xinteg8(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
ccccc         if(dmax*abs(xint).gt.half) then
              if(     abs(xint).gt.dmax) then
                xin2=xint+xint
c
c                 canonical order
c
                call descend(icf,jcf,kcf,lcf,iix)
c
                iff=iix(1)
                jff=iix(2)
                kff=iix(3)
                lff=iix(4)
c
                ii=lind(iff)
                jj=lind(jff)
                kk=lind(kff)
c
                ij=ii+jff
                kl=kk+lff
                ik=ii+kff
                il=ii+lff
c
                if(jff.ge.kff) then
                  jk=jj+kff
                else
                  jk=kk+jff
                end if
                if(jff.ge.lff) then
                  jl=jj+lff
                else
                  ll=lind(lff)
                  jl=ll+jff
                end if
c
                fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                fockA(ik)=fockA(ik)-xint*DenA(jl)
                fockA(jl)=fockA(jl)-xint*DenA(ik)
                fockA(jk)=fockA(jk)-xint*DenA(il)
                fockA(il)=fockA(il)-xint*DenA(jk)
c
                fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                fockB(ik)=fockB(ik)-xint*DenB(jl)
                fockB(jl)=fockB(jl)-xint*DenB(ik)
                fockB(jk)=fockB(jk)-xint*DenB(il)
                fockB(il)=fockB(il)-xint*DenB(jk)
c
              endif       !          (abs(xint).gt.half) then
c
            enddo         !   over lll=1,llen
          enddo         !   over kkk=1,klen
        enddo         !   over jjj=1,jlen
      enddo         !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_fock4_uhf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dena,        denb,       denab,      focka,
     $        fockb,       lind,       dmax,       intx4,
     $        xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dena(*),denb(*),denab(*),focka(*),fockb(*)
      integer*4 xinteg4(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
ccccc         if(dmax*abs(xint).gt.half) then
              if(     abs(xint).gt.dmax) then
                xin2=xint+xint
c
c                 canonical order
c
                call descend(icf,jcf,kcf,lcf,iix)
c
                iff=iix(1)
                jff=iix(2)
                kff=iix(3)
                lff=iix(4)
c
                ii=lind(iff)
                jj=lind(jff)
                kk=lind(kff)
c
                ij=ii+jff
                kl=kk+lff
                ik=ii+kff
                il=ii+lff
c
                if(jff.ge.kff) then
                  jk=jj+kff
                else
                  jk=kk+jff
                end if
                if(jff.ge.lff) then
                  jl=jj+lff
                else
                  ll=lind(lff)
                  jl=ll+jff
                end if
c
                fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                fockA(ik)=fockA(ik)-xint*DenA(jl)
                fockA(jl)=fockA(jl)-xint*DenA(ik)
                fockA(jk)=fockA(jk)-xint*DenA(il)
                fockA(il)=fockA(il)-xint*DenA(jk)
c
                fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                fockB(ik)=fockB(ik)-xint*DenB(jl)
                fockB(jl)=fockB(jl)-xint*DenB(ik)
                fockB(jk)=fockB(jk)-xint*DenB(il)
                fockB(il)=fockB(il)-xint*DenB(jk)
c
              endif       !          (abs(xint).gt.half) then
c
            enddo         !   over lll=1,llen
          enddo         !   over kkk=1,klen
        enddo         !   over jjj=1,jlen
      enddo         !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_fock2_uhf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dena,        denb,       denab,      focka,
     $        fockb,       lind,       dmax,       intx2,
     $        xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dena(*),denb(*),denab(*),focka(*),fockb(*)
      integer*2 xinteg2(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
ccccc         if(dmax*abs(xint).gt.half) then
              if(     abs(xint).gt.dmax) then
                xin2=xint+xint
c
c                 canonical order
c
                call descend(icf,jcf,kcf,lcf,iix)
c
                iff=iix(1)
                jff=iix(2)
                kff=iix(3)
                lff=iix(4)
c
                ii=lind(iff)
                jj=lind(jff)
                kk=lind(kff)
c
                ij=ii+jff
                kl=kk+lff
                ik=ii+kff
                il=ii+lff
c
                if(jff.ge.kff) then
                  jk=jj+kff
                else
                  jk=kk+jff
                end if
                if(jff.ge.lff) then
                  jl=jj+lff
                else
                  ll=lind(lff)
                  jl=ll+jff
                end if
c
                fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                fockA(ik)=fockA(ik)-xint*DenA(jl)
                fockA(jl)=fockA(jl)-xint*DenA(ik)
                fockA(jk)=fockA(jk)-xint*DenA(il)
                fockA(il)=fockA(il)-xint*DenA(jk)
c
                fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                fockB(ik)=fockB(ik)-xint*DenB(jl)
                fockB(jl)=fockB(jl)-xint*DenB(ik)
                fockB(jk)=fockB(jk)-xint*DenB(il)
                fockB(il)=fockB(il)-xint*DenB(jk)
c
              endif       !          (abs(xint).gt.half) then
c
            enddo         !   over lll=1,llen
          enddo         !   over kkk=1,klen
        enddo         !   over jjj=1,jlen
      enddo         !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_acmfock8_uhf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dena,        denb,       denab,      focka,
     $        fockb,       lind,       dmax,       intx8,
     $        xinteg8,     ax)
      implicit real*8 (a-h,o-z)
      dimension dena(*),denb(*),denab(*),focka(*),fockb(*)
      dimension xinteg8(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
ccccc         if(dmax*abs(xint).gt.half) then
              if(     abs(xint).gt.dmax) then
                xin2=xint+xint
                xint=ax*xint
c
c                 canonical order
c
                call descend(icf,jcf,kcf,lcf,iix)
c
                iff=iix(1)
                jff=iix(2)
                kff=iix(3)
                lff=iix(4)
c
                ii=lind(iff)
                jj=lind(jff)
                kk=lind(kff)
c
                ij=ii+jff
                kl=kk+lff
                ik=ii+kff
                il=ii+lff
c
                if(jff.ge.kff) then
                  jk=jj+kff
                else
                  jk=kk+jff
                end if
                if(jff.ge.lff) then
                  jl=jj+lff
                else
                  ll=lind(lff)
                  jl=ll+jff
                end if
c
                fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                fockA(ik)=fockA(ik)-xint*DenA(jl)
                fockA(jl)=fockA(jl)-xint*DenA(ik)
                fockA(jk)=fockA(jk)-xint*DenA(il)
                fockA(il)=fockA(il)-xint*DenA(jk)
c
                fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                fockB(ik)=fockB(ik)-xint*DenB(jl)
                fockB(jl)=fockB(jl)-xint*DenB(ik)
                fockB(jk)=fockB(jk)-xint*DenB(il)
                fockB(il)=fockB(il)-xint*DenB(jk)
c
              endif       !          (abs(xint).gt.half) then
c
            enddo         !   over lll=1,llen
          enddo         !   over kkk=1,klen
        enddo         !   over jjj=1,jlen
      enddo         !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_acmfock4_uhf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dena,        denb,       denab,      focka,
     $        fockb,       lind,       dmax,       intx4,
     $        xinteg4,     ax)
      implicit real*8 (a-h,o-z)
      dimension dena(*),denb(*),denab(*),focka(*),fockb(*)
      integer*4 xinteg4(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
ccccc         if(dmax*abs(xint).gt.half) then
              if(     abs(xint).gt.dmax) then
                xin2=xint+xint
                xint=ax*xint
c
c                 canonical order
c
                call descend(icf,jcf,kcf,lcf,iix)
c
                iff=iix(1)
                jff=iix(2)
                kff=iix(3)
                lff=iix(4)
c
                ii=lind(iff)
                jj=lind(jff)
                kk=lind(kff)
c
                ij=ii+jff
                kl=kk+lff
                ik=ii+kff
                il=ii+lff
c
                if(jff.ge.kff) then
                  jk=jj+kff
                else
                  jk=kk+jff
                end if
                if(jff.ge.lff) then
                  jl=jj+lff
                else
                  ll=lind(lff)
                  jl=ll+jff
                end if
c
                fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                fockA(ik)=fockA(ik)-xint*DenA(jl)
                fockA(jl)=fockA(jl)-xint*DenA(ik)
                fockA(jk)=fockA(jk)-xint*DenA(il)
                fockA(il)=fockA(il)-xint*DenA(jk)
c
                fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                fockB(ik)=fockB(ik)-xint*DenB(jl)
                fockB(jl)=fockB(jl)-xint*DenB(ik)
                fockB(jk)=fockB(jk)-xint*DenB(il)
                fockB(il)=fockB(il)-xint*DenB(jk)
c
              endif       !          (abs(xint).gt.half) then
c
            enddo         !   over lll=1,llen
          enddo         !   over kkk=1,klen
        enddo         !   over jjj=1,jlen
      enddo         !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_acmfock2_uhf(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dena,        denb,       denab,      focka,
     $        fockb,       lind,       dmax,       intx2,
     $        xinteg2,     ax)
      implicit real*8 (a-h,o-z)
      dimension dena(*),denb(*),denab(*),focka(*),fockb(*)
      integer*2 xinteg2(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
ccccc         if(dmax*abs(xint).gt.half) then
              if(     abs(xint).gt.dmax) then
                xin2=xint+xint
                xint=ax*xint
c
c                 canonical order
c
                call descend(icf,jcf,kcf,lcf,iix)
c
                iff=iix(1)
                jff=iix(2)
                kff=iix(3)
                lff=iix(4)
c
                ii=lind(iff)
                jj=lind(jff)
                kk=lind(kff)
c
                ij=ii+jff
                kl=kk+lff
                ik=ii+kff
                il=ii+lff
c
                if(jff.ge.kff) then
                  jk=jj+kff
                else
                  jk=kk+jff
                end if
                if(jff.ge.lff) then
                  jl=jj+lff
                else
                  ll=lind(lff)
                  jl=ll+jff
                end if
c
                fockA(ij)=fockA(ij)+xin2*DenAB(kl)
                fockA(kl)=fockA(kl)+xin2*DenAB(ij)
                fockA(ik)=fockA(ik)-xint*DenA(jl)
                fockA(jl)=fockA(jl)-xint*DenA(ik)
                fockA(jk)=fockA(jk)-xint*DenA(il)
                fockA(il)=fockA(il)-xint*DenA(jk)
c
                fockB(ij)=fockB(ij)+xin2*DenAB(kl)
                fockB(kl)=fockB(kl)+xin2*DenAB(ij)
                fockB(ik)=fockB(ik)-xint*DenB(jl)
                fockB(jl)=fockB(jl)-xint*DenB(ik)
                fockB(jk)=fockB(jk)-xint*DenB(il)
                fockB(il)=fockB(il)-xint*DenB(jk)
c
              endif       !          (abs(xint).gt.half) then
c
            enddo         !   over lll=1,llen
          enddo         !   over kkk=1,klen
        enddo         !   over jjj=1,jlen
      enddo         !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_dftuhf8(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dens,        focka,      fockb,      lind,
     $        dmax,        intx8,      xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dens(*),focka(*),fockb(*)
      dimension xinteg8(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx8=intx8+1
                     xint=xinteg8(intx8)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin2=xint+xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
c
                     fockA(ij)=fockA(ij)+xin2*dens(kl)
                     fockA(kl)=fockA(kl)+xin2*dens(ij)
c
                     fockB(ij)=fockB(ij)+xin2*dens(kl)
                     fockB(kl)=fockB(kl)+xin2*dens(ij)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_dftuhf4(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dens,        focka,      fockb,      lind,
     $        dmax,        intx4,      xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dens(*),focka(*),fockb(*)
      integer*4 xinteg4(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx4=intx4+1
                     xint=xinteg4(intx4)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin2=xint+xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
c
                     fockA(ij)=fockA(ij)+xin2*dens(kl)
                     fockA(kl)=fockA(kl)+xin2*dens(ij)
c
                     fockB(ij)=fockB(ij)+xin2*dens(kl)
                     fockB(kl)=fockB(kl)+xin2*dens(ij)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
c=====================================================================
      subroutine make_dftuhf2(
     $        ilen,        jlen,       klen,       llen,
     $        icff,        jcff,       kcff,       lcff,
     $        dens,        focka,      fockb,      lind,
     $        dmax,        intx2,      xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dens(*),focka(*),fockb(*)
      integer*2 xinteg2(*)
      dimension iix(4),lind(*)
      data half /0.5d0/
c
c      here dmax is 0.5/dmax and dmax*abs(x) > half is replaced by
c                                     abs(x) > dmax 
c
c
         do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
               do kcf=kcff+1,kcff+klen
                  do lcf=lcff+1,lcff+llen
                     intx2=intx2+1
                     xint=xinteg2(intx2)
ccccc             if(dmax*abs(xint).gt.half) then
                  if(     abs(xint).gt.dmax) then
                     xin2=xint+xint
canonical order
                     call descend(icf,jcf,kcf,lcf,iix)
c
                     iff=iix(1)
                     jff=iix(2)
                     kff=iix(3)
                     lff=iix(4)
c
                     ii=lind(iff)
                     kk=lind(kff)
c
                     ij=ii+jff
                     kl=kk+lff
c
                     fockA(ij)=fockA(ij)+xin2*dens(kl)
                     fockA(kl)=fockA(kl)+xin2*dens(ij)
c
                     fockB(ij)=fockB(ij)+xin2*dens(kl)
                     fockB(kl)=fockB(kl)+xin2*dens(ij)
c
                  endif         !          (abs(xint).gt.half) then
c
                  enddo         !   over lll=1,llen
               enddo            !   over kkk=1,klen
            enddo               !   over jjj=1,jlen
         enddo                  !   over iii=1,ilen
c
      end
