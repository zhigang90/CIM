      subroutine fock_bldr(nfock,ncf,buf,fock,dn,
     *                     iarray, nbls,ngcd,lnijkl,
     *                     labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c  input parameters
c  nfock= number of Fock matrices
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  fock = Fock matrices
c  dn   = density matrices
c  iarray = the triangular numbers, iarray(i,j)=i*(i-1)/2+j
c  nbls = integ. block-size , number of contracted shell quartets
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fock(nfock,*),dn(nfock,*)
      dimension iarray(ncf,ncf)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      call getrval('four',four)
      call getrval('half',half)
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
                xin4=four*xint
c  Make separate loop if nfock=1
           if(nfock.eq.1) then
          fock(1,iijf)=fock(1,iijf)+xin4*dn(1,kklf)
          fock(1,kklf)=fock(1,kklf)+xin4*dn(1,iijf)
          fock(1,iikf)=fock(1,iikf)-xint*dn(1,jjlf)
          fock(1,jjlf)=fock(1,jjlf)-xint*dn(1,iikf)
          fock(1,jjkf)=fock(1,jjkf)-xint*dn(1,iilf)
          fock(1,iilf)=fock(1,iilf)-xint*dn(1,jjkf)
	       else
               do i=1,nfock
          fock(i,iijf)=fock(i,iijf)+xin4*dn(i,kklf)
          fock(i,kklf)=fock(i,kklf)+xin4*dn(i,iijf)
          fock(i,iikf)=fock(i,iikf)-xint*dn(i,jjlf)
          fock(i,jjlf)=fock(i,jjlf)-xint*dn(i,iikf)
          fock(i,jjkf)=fock(i,jjkf)-xint*dn(i,iilf)
          fock(i,iilf)=fock(i,iilf)-xint*dn(i,jjkf)
               end do
             end if
             endif
c------------------------------------------------------
  200        continue
  150   continue
  100 continue
c
      end
c==============================================================
      subroutine fock_bldr2(nfock,ncf,buf,fock,dn,
     *                    iarray,
     *                    nbls,ngcd,lnijkl,
     *                    labels,length,lgenct )
c
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fock(nfock,*),dn(nfock,*)
      dimension  iarray(ncf,ncf)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      equivalence(iix(1),icf),(iix(2),jcf),(iix(3),kcf),(iix(4),lcf)
      dimension ishell(1296),jshell(1296),kshell(1296),lshell(1296)
      data four,half/4.0d0,0.5d0/
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
        do 50 iii=1,ilen
          do 50 jjj=1,jlen
            do 50 kkk=1,klen
              do 50 lll=1,llen
                intct=intct+1
                ishell(intct)=iii
                jshell(intct)=jjj
                kshell(intct)=kkk
                lshell(intct)=lll
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
            lcff=labels(4,iqu,ijklp)
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
              lcf=lcff+lshell(intct)
c
                iijf=iarray(jcf,icf)
                kklf=iarray(lcf,kcf)
                iikf=iarray(kcf,icf)
                iilf=iarray(lcf,icf)
                jjkf=iarray(kcf,jcf)
                jjlf=iarray(lcf,jcf)
c
                xin4=four*xint
c
               do i=1,nfock
          fock(i,iijf)=fock(i,iijf)+xin4*dn(i,kklf)
          fock(i,kklf)=fock(i,kklf)+xin4*dn(i,iijf)
          fock(i,iikf)=fock(i,iikf)-xint*dn(i,jjlf)
          fock(i,jjlf)=fock(i,jjlf)-xint*dn(i,iikf)
          fock(i,jjkf)=fock(i,jjkf)-xint*dn(i,iilf)
          fock(i,iilf)=fock(i,iilf)-xint*dn(i,jjkf)
               end do
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
      subroutine fock_bldr1(ncf,buf,fock,dn,
     *                      iarray,nbls,ngcd,lnijkl,
     *                      labels,length,lgenct )
c
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fock(*),dn(*)
      dimension  iarray(ncf,ncf)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      equivalence(iix(1),icf),(iix(2),jcf),(iix(3),kcf),(iix(4),lcf)
      dimension ishell(1296),jshell(1296),kshell(1296),lshell(1296)
      data four,half/4.0d0,0.5d0/
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
                xin4=four*xint
c
                iijf=iarray(jcf,icf)
                kklf=iarray(lcf,kcf)
          fock(iijf)=fock(iijf)+xin4*dn(kklf)
          fock(kklf)=fock(kklf)+xin4*dn(iijf)
                iikf=iarray(kcf,icf)
                jjlf=iarray(lcf,jcf)
          fock(iikf)=fock(iikf)-xint*dn(jjlf)
          fock(jjlf)=fock(jjlf)-xint*dn(iikf)
                iilf=iarray(lcf,icf)
                jjkf=iarray(kcf,jcf)
          fock(jjkf)=fock(jjkf)-xint*dn(iilf)
          fock(iilf)=fock(iilf)-xint*dn(jjkf)
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
      subroutine fock_bldr3(nfock,ncf,buf,fock,dn,
     *                    iarray,
     *                    nbls,ngcd,lnijkl,
     *                    labels,length,lgenct )
c
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fock(nfock,*),dn(nfock,*)
      dimension  iarray(ncf,ncf)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      equivalence(iix(1),icf),(iix(2),jcf),(iix(3),kcf),(iix(4),lcf)
      dimension ishell(1296),jshell(1296),kshell(1296)
      data four,half/4.0d0,0.5d0/
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
                xin4=four*xint
c
               do i=1,nfock
          fock(i,iijf)=fock(i,iijf)+xin4*dn(i,kklf)
          fock(i,kklf)=fock(i,kklf)+xin4*dn(i,iijf)
          fock(i,iikf)=fock(i,iikf)-xint*dn(i,jjlf)
          fock(i,jjlf)=fock(i,jjlf)-xint*dn(i,iikf)
          fock(i,jjkf)=fock(i,jjkf)-xint*dn(i,iilf)
          fock(i,iilf)=fock(i,iilf)-xint*dn(i,jjkf)
               end do
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
      subroutine fock_bldr4(ncf,buf,fock,dn,
     *                      iarray,nbls,ngcd,lnijkl,
     *                      labels,length,lgenct )
c
c----------------------------------------------------------------
c Fock builder
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
      dimension buf(nbls,lnijkl,ngcd)
      dimension fock(*),dn(*)
      dimension  iarray(ncf,ncf)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      equivalence(iix(1),icf),(iix(2),jcf),(iix(3),kcf),(iix(4),lcf)
      dimension ishell(1296),jshell(1296),kshell(1296)
      data four,half/4.0d0,0.5d0/
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
                xin4=four*xint
c
          fock(iijf)=fock(iijf)+xin4*dn(kklf)
          fock(kklf)=fock(kklf)+xin4*dn(iijf)
          fock(iikf)=fock(iikf)-xint*dn(jjlf)
          fock(jjlf)=fock(jjlf)-xint*dn(iikf)
          fock(jjkf)=fock(jjkf)-xint*dn(iilf)
          fock(iilf)=fock(iilf)-xint*dn(jjkf)
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
