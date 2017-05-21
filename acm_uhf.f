      subroutine acmuhf_core(ax, core,ilab,jlab,klab, llab,
     *                                isiz,jsiz,ksiz,lsiz,iqrt,
     *                                integrals,nqstore,nblstore,
     *                                thres0,thres1, lind,
     *                       DenA,   DenB,   DenAB,  fockA,  fockB)
c----------------------------------------------------------------
c input :
c core(*) - array with stored integrals
c i-llab()- starting indecies for integrals
c i-lsiz()- shell's sizes for each super-block
c iqrt()  - number of stored quartets for each super-block
c integrals - number of stored integrals
c nqstore   - number of stored quartets
c nblstore  - number of stored big-blocks
c----------------------------------------------------------------
c ACM Fock builder (using in-core integrals )
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
cccc  integer*2 isiz(nqstore),jsiz(nqstore),ksiz(nqstore),lsiz(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      dimension core(*)         ! dimension core(integrals)
      dimension DenA(*),DenB(*),DenAB(*),fockA(*),fockB(*)
      dimension lind(*)
      dimension iix(4)
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
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
c    ................................................
                       xint = ax*xint          ! acm
c    ................................................
                       lcf=lcff+lll
canonical   order
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
                    enddo         !   over lll=1,llen
                 enddo            !   over kkk=1,klen
              enddo               !   over jjj=1,jlen
           enddo                  !   over iii=1,ilen
         enddo                    !   do iqrt1=1,nqrt1
      enddo                       !   do isbl=1,nblstore
c
      if(thres0.ne.thres1) then
         resc=thres0/thres1
         ij=0
         do 200 i=1,ncf
         do 200 j=1,i
            ij=ij+1
            fockA(ij)=fockA(ij)*resc
            fockB(ij)=fockB(ij)*resc
  200    continue
      endif
c
      end
c=====================================================================
      subroutine acmuhf_bldr(ax,     ncf,    buf,    fockA,  fockB,
     *                       DenA,   DenB,   DenAB,  iarray, nbls,
     *                       ngcd,   lnijkl, labels, length, lgenct)
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------
c ACM Fock builder
c----------------------------------------------------------------
c  input parameters
c  ax     -  factor to multiply exchange term
c  ncf    -  number of contracted functions
c  buf    -  intgrals buffer
c  fockA  -  alpha Fock matrix
c  fockB  -  beta Fock matrix
c  DenA   -  alpha density matrix
c  DenB   -  beta density matrix
c  DenAB  -  alpha+beta density matrix
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
c  .................................................
                xint=ax*xint              ! acm
c  .................................................
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
      subroutine acmuhf_bldr1(ax,     ncf,    buf,    fockA,  fockB,
     *                        DenA,   DenB,   DenAB,  iarray, nbls,
     *                        ngcd,   lnijkl, labels, length, lgenct)
      implicit real*8 (a-h,o-z)
c
c----------------------------------------------------------------
c ACM Fock builder
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
c
                xin2=xint+xint
c  .................................................
                xint=ax*xint              ! acm
c  .................................................
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
      subroutine acmuhf_bldr3(ax,     ncf,    buf,    fockA,  fockB,
     *                        DenA,   DenB,   DenAB,  iarray, nbls,
     *                        ngcd,   lnijkl, labels, length, lgenct)
      implicit real*8 (a-h,o-z)
c
c----------------------------------------------------------------
c ACM Fock builder
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
c  .................................................
                xint=ax*xint              ! acm
c  .................................................
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
