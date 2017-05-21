c=====================================================================
      subroutine int_giao(idft,ax,nblocks,bl,inx,ntri,thres1,dscreen,
     *                    dens,fock,labels,mywork,igran)
c---------------------------------------------------------------------
c This routine is called from CHSIFT for all three modes :
c           non -, semi-, and full-direct.
c It calculates the GIAO inegrals for NMR chemical shift.
c The part of the first-order Fock matrix is constructed using GIAO
c two-electron integral derivatives . These integral derivatives are
c calculated only once and then used to construct the fock10 matrix.
c
c value of where='shif' is set up HERE
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NON EXISTING super block
c---------------------------------------------------------------------
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),dscreen(*),dens(*),fock(*)
      dimension labels(*)
c-----------------------------------------------------------------
      time0 = 0.0d0
      call secund(time1)
c----------------------------------------------------------------
c calculate the GIAO integrals
c
      where='shif'
ccc   screen='shif'
      screen='giao'
c----------------------------------------------------------------
c set inetgral timings to zero as well as neglect stuff
c
        call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c total number of super-blocks  nsupblks=nbl2*(nbl2+1)/2
c
      nsupblks=nblocks
      if(mywork.eq.0) then
         istart=1
         istop=nsupblks
      else
         istart=mywork
         istop=mywork+igran-1
      end if
c----------------------------------------------------------------
      do isupbl=istart,istop
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,dscreen,where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
c
c check if requested super-block was in a list :
c
        if(stopnow) go to 1111
c
c integrals arrived in bl(ibuffz); check if they are :
c
c  nintez is the size of a given quartet.It is set up
c  to zero if there is no integrals
c
        if(nintez.gt.0) then
           call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
ckw        call get_ranges(nblsiz,ngctoz,nintez,
ckw  *                     labels(lab1),labels(lab2),labels(lab3),
ckw  *                     ir1,ir2, jr1,jr2, kr1,kr2, lr1,lr2)
           call fock_giao(idft,ax,bl(ibuffz),dens,fock,bl(lind),ntri,
     *                    nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
        endif
c
        if(moreint) then
           go to 11
        endif
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
      call zerout_fock_diag(fock,ntri,bl(lind))
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
      call secund(time2)
      call term_info(thres1,time2-time1,time1-time0,where)
c
c----------------------------------------------------------------
c
      end
c=====================================================================
c Fock-GIAO-builder routine F(D0,G1) :
c
      subroutine fock_giao(idft,ax,buf,dens,fock,lind,ntri,
     *                     nbls,ngcd,lnijkl, labels,length,lgenct )
c----------------------------------------------------------------
c F(D0,G1) builder (using GIAO integrals from buffer)
c----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
      dimension buf(6,nbls,lnijkl,ngcd)
      dimension dens(ntri),fock(ntri,3)
      dimension lind(*)
      dimension iix(4)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      logical doexch
c
      parameter (Zero=0.0d0)
c----------------------------------------------------------------
c needed for scaling exchange contribution in dft
c
      halfax=0.5d0
      if (idft.gt.0) then
         halfax=halfax*ax
         doexch=.false.
      endif
      if(idft.eq.0.or.ax.NE.Zero) doexch=.true.
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
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             ii=lind(icf)
             do 200 jjj=1,jlen
             jcf=jcff+jjj
             jj=lind(jcf)
               if(icf.ge.jcf) then
                 ijf=ii+jcf
                 fij=1.d0
               else
                 ijf=jj+icf
                 fij=-1.d0
               endif
             do 200 kkk=1,klen
             kcf=kcff+kkk
             kk=lind(kcf)
               if(icf.ge.kcf) then
                  ikf=ii+kcf
                  fik=1.d0
               else
                  ikf=kk+icf
                  fik=-1.d0
               endif
               if(jcf.ge.kcf) then
                  jkf=jj+kcf
                  fjk=1.d0
               else
                  jkf=kk+jcf
                  fjk=-1.d0
               endif
             do 200 lll=1,llen
             lcf=lcff+lll
             ll=lind(lcf)
               if(icf.ge.lcf) then
                  ilf=ii+lcf
                  fil=-1.d0
               else
                  ilf=ll+icf
                  fil=1.d0
               endif
               if(jcf.ge.lcf) then
                  jlf=jj+lcf
                  fjl=-1.d0
               else
                  jlf=ll+jcf
                  fjl=1.d0
               endif
               if(kcf.ge.lcf) then
                  klf=kk+lcf
                  fkl=1.d0
               else
                  klf=ll+kcf
                  fkl=-1.d0
               endif
c
             integ=integ+1
c------------------------------------------
             xder1=buf(1,ijklp,integ,iqu)
             yder1=buf(3,ijklp,integ,iqu)
             zder1=buf(5,ijklp,integ,iqu)
             xder2=buf(2,ijklp,integ,iqu)
             yder2=buf(4,ijklp,integ,iqu)
             zder2=buf(6,ijklp,integ,iqu)
c--------------------------------------------------------
c            print integrals
c              write(6,1234) icf,jcf,kcf,lcf,xder1*1.d-10,xder2*1.d-10
c              write(6,1235)                 yder1*1.d-10,yder2*1.d-10
c              write(6,1235)                 zder1*1.d-10,zder2*1.d-10
c
c 1234         format(4i3,1x,2(f12.8,1x))
c 1235         format(13x,   2(f12.8,1x))
c--------------------------------------------------------
               dklij=fkl*dens(ijf)
               dijkl=fij*dens(klf)
c--------------------------------------------------------
c
c Coulomb part :
c
c  ***  ij, kl  ***
                fock(ijf,1)=fock(ijf,1)+dijkl*(xder1-xder2)
                fock(klf,1)=fock(klf,1)+dklij*(xder1+xder2)
c
                fock(ijf,2)=fock(ijf,2)+dijkl*(yder1-yder2)
                fock(klf,2)=fock(klf,2)+dklij*(yder1+yder2)
c
                fock(ijf,3)=fock(ijf,3)+dijkl*(zder1-zder2)
                fock(klf,3)=fock(klf,3)+dklij*(zder1+zder2)
c
c--------------
c Exchange part: not needed for pure dft calculations
c
              IF(doexch) THEN
c
                diljk=fil*dens(jkf)*halfax
                djlik=fjl*dens(ikf)*halfax
                dikjl=fik*dens(jlf)*halfax
                djkil=fjk*dens(ilf)*halfax
c
c  ***  il,jl, ik,jk  ***
c
                fock(ilf,1)=fock(ilf,1) + diljk*xder1
                fock(jlf,1)=fock(jlf,1) + djlik*xder2
                fock(ikf,1)=fock(ikf,1) + dikjl*xder2
                fock(jkf,1)=fock(jkf,1) + djkil*xder1
c
                fock(ilf,2)=fock(ilf,2) + diljk*yder1
                fock(jlf,2)=fock(jlf,2) + djlik*yder2
                fock(ikf,2)=fock(ikf,2) + dikjl*yder2
                fock(jkf,2)=fock(jkf,2) + djkil*yder1
c
                fock(ilf,3)=fock(ilf,3) + diljk*zder1
                fock(jlf,3)=fock(jlf,3) + djlik*zder2
                fock(ikf,3)=fock(ikf,3) + dikjl*zder2
                fock(jkf,3)=fock(jkf,3) + djkil*zder1
c
              ENDIF
c--------------
  250        continue
c--------------
  200        continue
c------------------------------------------------------
  150   continue
  100 continue
c------------------------------------------------------
c
      end
c==============================================================
      subroutine get_ranges(nbls,ngcd,lnijkl,labels,length,lgenct,
     *                      ir1,ir2, jr1,jr2, kr1,kr2, lr1,lr2)
c----------------------------------------------------------------
c Output : range of indeces :
c ir1,ir2, jr1,jr2, kr1,kr2, lr1,lr2
c----------------------------------------------------------------
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        ir1=10 000
        jr1=10 000
        kr1=10 000
        lr1=10 000
        ir2=0
        jr2=0
        kr2=0
        lr2=0
c
        do 100 ijklp=1,nbls
        ngcq=lgenct(ijklp)
          do 150 iqu=1,ngcq
          icff=labels(1,iqu,ijklp)
          jcff=labels(2,iqu,ijklp)
          kcff=labels(3,iqu,ijklp)
          lcff=labels(4,iqu,ijklp)
c
          icf1=icff+1
          icf2=icff+ilen
          jcf1=jcff+1
          jcf2=jcff+jlen
          kcf1=kcff+1
          kcf2=kcff+klen
          lcf1=lcff+1
          lcf2=lcff+llen
c
            if(icf1.lt.ir1) ir1=icf1
            if(jcf1.lt.jr1) jr1=jcf1
            if(kcf1.lt.kr1) kr1=kcf1
            if(lcf1.lt.lr1) lr1=lcf1
c
            if(icf2.gt.ir2) ir2=icf2
            if(jcf2.gt.jr2) jr2=jcf2
            if(kcf2.gt.kr2) kr2=kcf2
            if(lcf2.gt.lr2) lr2=lcf2
c
  150     continue
  100   continue
c------------------------------------------------------
      end
c====================================================================
