c=====================================================================
      subroutine int_d0g2(idft,ax,rhf,nblocks,bl,inx,ntri,threg2,
     *                    densp , ncenter,
     *                    dens,denB,hessian,labels,mywork,igran)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c (2) the second deriv.integ. together with the zero-order density
c     give direct contributions to the final hessian Tr{ D0*G(D0,g2) }
c
c---------------------------------------------------------------------
c This routine is called from para_d0g2 (called from hess2e2)
c for all three modes : non -, semi-, and full-direct.
c 2nd-derivative Integral g2 are calculated only once
c and contracted with two-particle density giving contrib.
c to the final hessian
c
c value of where='hess' is set up HERE
c---------------------------------------------------------------------
c INPUT :
c  idft              -  dft flag:  0 - no dft
c  ax                - factor to multiply exchange contribution for ACM
c  rhf               - logical flag for closed-shell
c  nblocks           - number of blocks of contracted shell quartets
c  bl                - scratch for everything
c  inx(12,ncs)       - basis set & nuclear data
c  ntri              - ncf*(ncf+1)/2
c  threg2            - integral threshold for G(D0,g2)
c  densp(ncs,ncs)    - screening density
c  ncenter(ncf)      - center of each basis function
c  dens(ntri)        - alpha/closed-shell density matrix
c  denB(ntri)        - beta density matrix
c  labels            - labels buffer
c  mywork            - for parallel runs
c  igran             - for parallel runs
C OUTPUT:
c hessian(3,natoms,3,natoms) - hessian matrix
c---------------------------------------------------------------------
      logical rhf,moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
      common /runtype/ scftype,where
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NON EXISTING super block
c---------------------------------------------------------------------
c     common /intbl/ifpp,inxx(100)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /lindvec/ lind,idensp
      dimension inx(12,*)
      dimension densp(*)
      dimension ncenter(*)
      dimension bl(*),dens(*),denB(*)
      dimension labels(*)
      dimension hessian(*)
c-----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c get print flag
c
      call getival('printh',nprint)
c----------------------------------------------------------------
c for symmetrization of two-el. part of hessian for exact symm.
c
      call getival('nsym',nsym)
c----------------------------------------------------------------
      where='hess'
c----------------------------------------------------------------
c set integral timings to zero as well as neglect stuff
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
c loop over blocks of contracted shell quartets
c
         call do_intder2(idft,ax,rhf,istart,istop,
ccc  *                    bl,inx,threg2,bl(idensp),labels,
     *                    bl,inx,threg2,    densp ,labels,
     *                    dens,denB,hessian,
     *                    bl(lind),ntri,   ncenter ,na,nsym )
c----------------------------------------------------------------
c Print the hessian matrix. Tt contains only Tr{ d0*G(d0,g2) } contr.
c
      if(nprint.ge.2) then
         write(6,*)' HESSIAN with Tr{ D0*G(D0,g2) } contributions'
         call hess_print(hessian,na,threg2,
     *                   ' contr. no=1 : Tr{D0*G(D0,g2) } ')
      endif
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
      call secund(time2)
      call term_info(threg2,time2-time1,time1-time1,where)
c
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine make_78_from_45(der2,work)
      implicit real*8 (a-h,o-z)
      dimension der2(45),work(78)
c
c block aa:
       axax=der2(1)
       axay=der2(2)
       axaz=der2(3)
       ayay=der2(4)
       ayaz=der2(5)
       azaz=der2(6)
c
c block ab:
       axbx=der2(7)
       axby=der2(8)
       axbz=der2(9)
       aybx=der2(10)
       ayby=der2(11)
       aybz=der2(12)
       azbx=der2(13)
       azby=der2(14)
       azbz=der2(15)
c block ac:
       axcx=der2(16)
       axcy=der2(17)
       axcz=der2(18)
       aycx=der2(19)
       aycy=der2(20)
       aycz=der2(21)
       azcx=der2(22)
       azcy=der2(23)
       azcz=der2(24)
c block bb:
       bxbx=der2(25)
       bxby=der2(26)
       bxbz=der2(27)
       byby=der2(28)
       bybz=der2(29)
       bzbz=der2(30)
c block bc:
       bxcx=der2(31)
       bxcy=der2(32)
       bxcz=der2(33)
       bycx=der2(34)
       bycy=der2(35)
       bycz=der2(36)
       bzcx=der2(37)
       bzcy=der2(38)
       bzcz=der2(39)
c block cc:
       cxcx=der2(40)
       cxcy=der2(41)
       cxcz=der2(42)
       cycy=der2(43)
       cycz=der2(44)
       czcz=der2(45)
c block ad: from transl. inv.
       axdx=-(axax+axbx+axcx)
       axdy=-(axay+axby+axcy)
       axdz=-(axaz+axbz+axcz)
       aydx=-(aXaY+aybx+aycx)
       aydy=-(ayay+ayby+aycy)
       aydz=-(ayaz+aybz+aycz)
       azdx=-(aXaZ+azbx+azcx)
       azdy=-(aYaZ+azby+azcy)
       azdz=-(azaz+azbz+azcz)
c block bd: from transl. inv.
       bxdx=-(AxBx+bxbx+bxcx)
       bxdy=-(AYBX+bxby+bxcy)
       bxdz=-(AZBX+bxbz+bxcz)
       bydx=-(AXBY+bXbY+bycx)
       bydy=-(AYBY+byby+bycy)
       bydz=-(AZBY+bybz+bycz)
       bzdx=-(AXBZ+bXbZ+bzcx)
       bzdy=-(AyBz+bybz+bzcy)
       bzdz=-(AZBZ+bzbz+bzcz)
c block cd: from transl. inv.
       cxdx=-(AXCX+BXCX+cxcx)
       cxdy=-(AYCX+BYCX+cxcy)
       cxdz=-(AZCX+BZCX+cxcz)
       cydx=-(AXCY+BXCY+cXcY)
       cydy=-(AYCY+BYCY+cycy)
       cydz=-(AZCY+BZCY+cycz)
       czdx=-(AXCZ+BxCZ+cXcZ)
       czdy=-(AYCZ+BYCZ+cYcZ)
       czdz=-(AZCZ+BZCZ+czcz)
c block dd:
       dxdx=-(AXDX+BXDX+CXDX)
cccc   dxdy=-(dxay+dxby+dxcy)
       dxdy=-(AYDX+BYDX+CYDX)
cccc   dxdz=-(dxaz+dxbz+dxcz)
       dxdz=-(AZDX+BZDX+CZDX)
cccc   dydy=-(dyay+dyby+dycy)
       dydy=-(AYDY+BYDY+CYDY)
ccccc  dydz=-(dyaz+dybz+dycz)
       dydz=-(AZDY+BZDY+CZDY)
cccc   dzdz=-(dzaz+dzbz+dzcz)
       dzdz=-(AZDZ+BZDZ+CZDZ)
c------------------------------------------------------------
c construct all 10 blocks of sec.der. (output) from 6 blocks:
c
c          AA AB AC AD                AA AB AC
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c------------------------------------------------------------
c first 24 derivatives are in right order :
c
c blocks AA,AB,AC :
       do m=1,24
          work(m)=der2(m)
       enddo
c
c block AD :
       work(25)=axdx    !  =-(axax+axbx+axcx)
       work(26)=axdy    !  =-(axay+axby+axcy)
       work(27)=axdz    !  =-(axaz+axbz+axcz)
       work(28)=aydx    !  =-(aXaY+aybx+aycx)
       work(29)=aydy    !  =-(ayay+ayby+aycy)
       work(30)=aydz    !  =-(ayaz+aybz+aycz)
       work(31)=azdx    !  =-(aXaZ+azbx+azcx)
       work(32)=azdy    !  =-(aYaZ+azby+azcy)
       work(33)=azdz    !  =-(azaz+azbz+azcz)
c
c block BB :
       work(34)=bxbx    !  =der2(25)
       work(35)=bxby    !  =der2(26)
       work(36)=bxbz    !  =der2(27)
       work(37)=byby    !  =der2(28)
       work(38)=bybz    !  =der2(29)
       work(39)=bzbz    !  =der2(30)
c
c block BC :
       work(40)=bxcx    !  =der2(31)
       work(41)=bxcy    !  =der2(32)
       work(42)=bxcz    !  =der2(33)
       work(43)=bycx    !  =der2(34)
       work(44)=bycy    !  =der2(35)
       work(45)=bycz    !  =der2(36)
       work(46)=bzcx    !  =der2(37)
       work(47)=bzcy    !  =der2(38)
       work(48)=bzcz    !  =der2(39)
c
c block BD :
       work(49)=bxdx    !  =-(AxBx+bxbx+bxcx)
       work(50)=bxdy    !  =-(AYBX+bxby+bxcy)
       work(51)=bxdz    !  =-(AZBX+bxbz+bxcz)
       work(52)=bydx    !  =-(AXBY+bXbY+bycx)
       work(53)=bydy    !  =-(AYBY+byby+bycy)
       work(54)=bydz    !  =-(AZBY+bybz+bycz)
       work(55)=bzdx    !  =-(AXBZ+bXbZ+bzcx)
       work(56)=bzdy    !  =-(AyBz+bybz+bzcy)
       work(57)=bzdz    !  =-(AZBZ+bzbz+bzcz)
c
c block CC :
       work(58)=cxcx    !  =der2(40)
       work(59)=cxcy    !  =der2(41)
       work(60)=cxcz    !  =der2(42)
       work(61)=cycy    !  =der2(43)
       work(62)=cycz    !  =der2(44)
       work(63)=czcz    !  =der2(45)
c
c block CD :
       work(64)=cxdx    !  =-(AXCX+BXCX+cxcx)
       work(65)=cxdy    !  =-(AYCX+BYCX+cxcy)
       work(66)=cxdz    !  =-(AZCX+BZCX+cxcz)
       work(67)=cydx    !  =-(AXCY+BXCY+cXcY)
       work(68)=cydy    !  =-(AYCY+BYCY+cycy)
       work(69)=cydz    !  =-(AZCY+BZCY+cycz)
       work(70)=czdx    !  =-(AXCZ+BxCZ+cXcZ)
       work(71)=czdy    !  =-(AYCZ+BYCZ+cYcZ)
       work(72)=czdz    !  =-(AZCZ+BZCZ+czcz)
c
c block DD :
       work(73)=dxdx    !  =-(AXDX+BXDX+CXDX)
       work(74)=dxdy    !  =-(AYDX+BYDX+CYDX)
       work(75)=dxdz    !  =-(AZDX+BZDX+CZDX)
       work(76)=dydy    !  =-(AYDY+BYDY+CYDY)
       work(77)=dydz    !  =-(AZDY+BZDY+CZDY)
       work(78)=dzdz    !  =-(AZDZ+BZDZ+CZDZ)
c--------------------------------------------------------------------
      end
c==============================================================
      subroutine prtn_der2(buf, lind,ntri,nbls,ngcd,lnijkl,
     *                     ncenter,labels,length,lgenct )
c----------------------------------------------------------------
c prints out second derivative two-electron integrals
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral derivatives buffer
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(45,nbls,lnijkl,ngcd)
      dimension lind(*) , ncenter(*)
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension work(78)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from hess_mtrx () :NBLS=',nbls,' ngcd=',ngcd
      write(6,*)'    --------------------------------------'
      write(6,*)'    Second derivatives of 2-el. integrals '
      write(6,*)'    --------------------------------------'
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
          do 150 iqu=1,ngcd
          icff=labels(1,iqu,ijklp)
          jcff=labels(2,iqu,ijklp)
          kcff=labels(3,iqu,ijklp)
          lcff=labels(4,iqu,ijklp)
c
          iat=ncenter(icff+1)
          jat=ncenter(jcff+1)
          kat=ncenter(kcff+1)
          lat=ncenter(lcff+1)
c........
c how many different centers :
c
           nc=1
           if(jat.ne.iat) nc=nc+1
           if(kat.ne.iat .and. kat.ne.jat) nc=nc+1
           if(lat.ne.iat .and. lat.ne.jat .and. lat.ne.kat) nc=nc+1
c
cccccc    write(6,*)' Ncent=',nc,' are : ',iat,jat,kat,lat,' qrt=',ijklp
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             ii=lind(icf)
                do 250 jjj=1,jlen
                jcf=jcff+jjj
                jj=lind(jcf)
                ijf=ii+jcf
                if(jcf.gt.icf) ijf=jj+icf
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      ilf=ii+lcf
                      jlf=jj+lcf
                      klf=kk+lcf
                      if(lcf.gt.icf) ilf=ll+icf
                      if(lcf.gt.jcf) jlf=ll+jcf
                      if(lcf.gt.kcf) klf=ll+kcf
c
                      integ=integ+1
c--------------------------------------------------------
        call make_78_from_45(buf(1,ijklp,integ,iqu),work)
        fact=1.d-10
        call dscal(78,fact,work,1)
        call print_pnl2(nc,iat,jat,kat,lat, work, icf,jcf,kcf,lcf)
c
c--------------------------------------------------------
c in work 78 derivatives ordered as :
c
c   all 10 blocks of sec.der. made out of 6 blocks:
c
c          AA AB AC AD                AA AB AC
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
  150   continue
  100 continue
      end
c=====================================================================
      subroutine print_pnl2(ncent,iat,jat,kat,lat,
     *                      der2, icf,jcf,kcf,lcf)
c PRINT ONLY
      implicit real*8 (a-h,o-z)
      dimension iix(4)
      dimension der2(78)
c--------------------------------------------------------------------
c     write(6,*)
c    *'Atoms =',iat,jat,kat,lat,'  (',ncent,' different centers)'
c     write(6,*)' '
c
      write(6,60) iat,jat,kat,lat,ncent
  60  format('Atoms =',4(i3,1x),'   (',i2,'  different centers)'/)
c--------------------------------------------------------------------
canonical order :
c
                   call descend(icf,jcf,kcf,lcf,iix)
c--------------------------------------------------------------------
c block AA:
                   write(6,61)'d2/dAidAj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=1,6)
  61  format(a16,4i2,1x,3(f12.6,2x)/39x,2(f12.6,2x)/53x,1(f12.6,2x))
c
c block AB:
                   write(6,62)'d2/dAidBj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=7,15)
  62  format(a16,4i2,1x,3(f12.6,2x)/25x,3(f12.6,2x)/25x,3(f12.6,2x))
c
c block AC:
                   write(6,62)'d2/dAidCj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=16,24)
c
c block AD:
                   write(6,62)'d2/dAidDj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=25,33)
c
c block BB:
                   write(6,61)'d2/dBidBj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=34,39)
c
c block BC:
                   write(6,62)'d2/dBidCj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=40,48)
c
c block BD:
                   write(6,62)'d2/dBidDj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=49,57)
c
c block CC:
                   write(6,61)'d2/dCidCj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=58,63)
c
c block CD:
                   write(6,62)'d2/dCidDj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=64,72)
c
c block DD:
                   write(6,61)'d2/dDidDj: ijkl=', icf,jcf,kcf,lcf,
     *                 (der2(ii),ii=73,78)
c--------------------------------------------------------------------
c
      end
c==============================================================
      subroutine hess_mtrx(buf, dens, hess, lind, thres1,
     *                     ntri,nbls, ngcd, lnijkl,
     *                     ncenter,labels, length, lgenct,
     *                     na, nsym, idft, ax )
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral derivatives buffer
c  dens()   - density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c  na       - number of atoms
c  idft     - dft flag (0 - no DFT)
c  ax       - factor to multiply exchange contribution (for hybrid DFT)
c
c OUTPUT :
c hess       () - hessian matrix
c----------------------------------------------------------------
c ordering of 45 second derivative integrals :
c
c    der2_AxAx=buf(1,nbls,lnijkl,ngcd),
c    der2_AxAy=buf(2,nbls,lnijkl,ngcd),
c    der2_AxAz=buf(3,nbls,lnijkl,ngcd),
c    der2_AyAy=buf(4,nbls,lnijkl,ngcd),
c    der2_AyAz=buf(5,nbls,lnijkl,ngcd),
c    der2_AzAz=buf(6,nbls,lnijkl,ngcd),
c
c    der2_AxBx= buf(7,nbls,lnijkl,ngcd),
c    der2_AxBy= buf(8,nbls,lnijkl,ngcd),
c    der2_AxBz= buf(9,nbls,lnijkl,ngcd),
c    der2_AyBx=buf(10,nbls,lnijkl,ngcd),
c    der2_AyBy=buf(11,nbls,lnijkl,ngcd),
c    der2_AyBz=buf(12,nbls,lnijkl,ngcd),
c    der2_AzBx=buf(13,nbls,lnijkl,ngcd),
c    der2_AzBy=buf(14,nbls,lnijkl,ngcd),
c    der2_AzBz=buf(15,nbls,lnijkl,ngcd),
c
c    der2_AxCx=buf(16,nbls,lnijkl,ngcd),
c    der2_AxCy=buf(17,nbls,lnijkl,ngcd),
c    der2_AxCz=buf(18,nbls,lnijkl,ngcd),
c    der2_AyCx=buf(19,nbls,lnijkl,ngcd),
c    der2_AyCy=buf(20,nbls,lnijkl,ngcd),
c    der2_AyCz=buf(21,nbls,lnijkl,ngcd),
c    der2_AzCx=buf(22,nbls,lnijkl,ngcd),
c    der2_AzCy=buf(23,nbls,lnijkl,ngcd),
c    der2_AzCz=buf(24,nbls,lnijkl,ngcd),
c
c    der2_BxBx=buf(25,nbls,lnijkl,ngcd),
c    der2_BxBy=buf(26,nbls,lnijkl,ngcd),
c    der2_BxBz=buf(27,nbls,lnijkl,ngcd),
c    der2_ByBy=buf(28,nbls,lnijkl,ngcd),
c    der2_ByBz=buf(29,nbls,lnijkl,ngcd),
c    der2_BzBz=buf(30,nbls,lnijkl,ngcd),
c
c    der2_BxCx=buf(31,nbls,lnijkl,ngcd),
c    der2_BxCy=buf(32,nbls,lnijkl,ngcd),
c    der2_BxCz=buf(33,nbls,lnijkl,ngcd),
c    der2_ByCx=buf(34,nbls,lnijkl,ngcd),
c    der2_ByCy=buf(35,nbls,lnijkl,ngcd),
c    der2_ByCz=buf(36,nbls,lnijkl,ngcd),
c    der2_BzCx=buf(37,nbls,lnijkl,ngcd),
c    der2_BzCy=buf(38,nbls,lnijkl,ngcd),
c    der2_BzCz=buf(39,nbls,lnijkl,ngcd),
c
c    der2_CxCx=buf(40,nbls,lnijkl,ngcd),
c    der2_CxCy=buf(41,nbls,lnijkl,ngcd),
c    der2_CxCz=buf(42,nbls,lnijkl,ngcd),
c    der2_CyCy=buf(43,nbls,lnijkl,ngcd),
c    der2_CyCz=buf(44,nbls,lnijkl,ngcd),
c    der2_CzCz=buf(45,nbls,lnijkl,ngcd),
c----------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(45,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension hess(3,na, 3,na)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension wor1(78)
      dimension work(45)
      dimension x(3,4,3,4)
      data zero /0.d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from hess_mtrx () :NBLS=',nbls,' ngcd=',ngcd
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
          do 150 iqu=1,ngcd
          icff=labels(1,iqu,ijklp)
          jcff=labels(2,iqu,ijklp)
          kcff=labels(3,iqu,ijklp)
          lcff=labels(4,iqu,ijklp)
c
          iat=ncenter(icff+1)
          jat=ncenter(jcff+1)
          kat=ncenter(kcff+1)
          lat=ncenter(lcff+1)
c
cccccc    call zeroit(x,144)  ! 144=12*12
          call zeroit(work,45)
cccccc    call zeroit(wor1,78)
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             ii=lind(icf)
                do 250 jjj=1,jlen
                jcf=jcff+jjj
                jj=lind(jcf)
                ijf=ii+jcf
                if(jcf.gt.icf) ijf=jj+icf
                dij4=dens(ijf)*4.d0
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                   dik=dens(ikf)
                   djk=dens(jkf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      ilf=ii+lcf
                      jlf=jj+lcf
                      klf=kk+lcf
                      if(lcf.gt.icf) ilf=ll+icf
                      if(lcf.gt.jcf) jlf=ll+jcf
                      if(lcf.gt.kcf) klf=ll+kcf
c
                      dil=dens(ilf)
                      djl=dens(jlf)
                      dkl=dens(klf)
c
                      If(idft.eq.0) Then
                        dijkl=dij4*dkl - dik*djl -dil*djk
                      Else
                        dijkl=dij4*dkl - ax*(dik*djl +dil*djk)
                      EndIf
c
                      integ=integ+1
c--------------------------------------------------------
          if(abs(dijkl).lt.thres1) go to 350
c--------------------------------------------------------
c in work 78 derivatives ordered as :
c
c   all 10 blocks of sec.der. made out of 6 blocks:
c
c          AA AB AC AD                AA AB AC
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c--------------------------------------------------------
        call daxpy(45,dijkl,buf(1,ijklp,integ,iqu),1,work,1)
c
c       do iel=1,45
c          work(iel)=work(iel)+dijkl*buf(iel,ijklp,integ,iqu)
c       enddo
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
        call make_78_from_45(work,wor1)
c--------------------------------------------------------
c     call x_aa(x,1,wor1,1)      !   Block AA: wor1(1-6)
c
      x(1,1 , 1,1 )= wor1(1)   !  AxAx
      x(1,1 , 2,1 )= wor1(2)   !  AxAy
      x(1,1 , 3,1 )= wor1(3)   !  AxAz
      x(2,1 , 2,1 )= wor1(4)   !  AyAy
      x(2,1 , 3,1 )= wor1(5)   !  AyAz
      x(3,1 , 3,1 )= wor1(6)   !  AzAz
c
c--------------------------------------------------------
c     call x_aa(x,2,wor1,34)     !   block BB: wor1(34-39)
c
      x(1,2 , 1,2 )= wor1(34)   !  BxBx
      x(1,2 , 2,2 )= wor1(35)   !  BxBy
      x(1,2 , 3,2 )= wor1(36)   !  BxBz
      x(2,2 , 2,2 )= wor1(37)   !  ByBy
      x(2,2 , 3,2 )= wor1(38)   !  ByBz
      x(3,2 , 3,2 )= wor1(39)   !  BzBz
c
c--------------------------------------------------------
c     call x_aa(x,3,wor1,58)     !   block CC: wor1(58-63)
c
      x(1,3 , 1,3 )= wor1(58)   !  CxCx
      x(1,3 , 2,3 )= wor1(59)   !  CxCy
      x(1,3 , 3,3 )= wor1(60)   !  CxCz
      x(2,3 , 2,3 )= wor1(61)   !  CyCy
      x(2,3 , 3,3 )= wor1(62)   !  CyCz
      x(3,3 , 3,3 )= wor1(63)   !  CzCz

c--------------------------------------------------------
c     call x_aa(x,4,wor1,73)     !   block DD: wor1(73-78)
c
      x(1,4 , 1,4 )= wor1(73)   !  DxDx
      x(1,4 , 2,4 )= wor1(74)   !  DxDy
      x(1,4 , 3,4 )= wor1(75)   !  DxDz
      x(2,4 , 2,4 )= wor1(76)   !  DyDy
      x(2,4 , 3,4 )= wor1(77)   !  DyDz
      x(3,4 , 3,4 )= wor1(78)   !  DzDz
c
c--------------------------------------------------------
c     call x_ab(x,1,2, wor1,7)   !   block AB: wor1(7-15)
c
      x(1,1 , 1,2 )= wor1(7)   !  AxBx
      x(1,1 , 2,2 )= wor1(8)   !  AxBy
      x(1,1 , 3,2 )= wor1(9)   !  AxBz
      x(2,1 , 1,2 )= wor1(10)  !  AyBx
      x(2,1 , 2,2 )= wor1(11)  !  AyBy
      x(2,1 , 3,2 )= wor1(12)  !  AyBz
      x(3,1 , 1,2 )= wor1(13)  !  AzBx
      x(3,1 , 2,2 )= wor1(14)  !  AzBy
      x(3,1 , 3,2 )= wor1(15)  !  AzBz
c
      x(1,2 , 1,1 )= wor1(7)   !  AxBx
      x(2,2 , 1,1 )= wor1(8)   !  AxBy
      x(3,2 , 1,1 )= wor1(9)   !  AxBz
      x(1,2 , 2,1 )= wor1(10)  !  AyBx
      x(2,2 , 2,1 )= wor1(11)  !  AyBy
      x(3,2 , 2,1 )= wor1(12)  !  AyBz
      x(1,2 , 3,1 )= wor1(13)  !  AzBx
      x(2,2 , 3,1 )= wor1(14)  !  AzBy
      x(3,2 , 3,1 )= wor1(15)  !  AzBz
c--------------------------------------------------------
c     call x_ab(x,1,3, wor1,16)  !   block AC: wor1(16-24)
c
      x(1,1 , 1,3 )= wor1(16)  !  AxCx
      x(1,1 , 2,3 )= wor1(17)  !  AxCy
      x(1,1 , 3,3 )= wor1(18)  !  AxCz
      x(2,1 , 1,3 )= wor1(19)  !  AyCx
      x(2,1 , 2,3 )= wor1(20)  !  AyCy
      x(2,1 , 3,3 )= wor1(21)  !  AyCz
      x(3,1 , 1,3 )= wor1(22)  !  AzCx
      x(3,1 , 2,3 )= wor1(23)  !  AzCy
      x(3,1 , 3,3 )= wor1(24)  !  AzCz
c
      x(1,3 , 1,1 )= wor1(16)  !  AxCx
      x(2,3 , 1,1 )= wor1(17)  !  AxCy
      x(3,3 , 1,1 )= wor1(18)  !  AxCz
      x(1,3 , 2,1 )= wor1(19)  !  AyCx
      x(2,3 , 2,1 )= wor1(20)  !  AyCy
      x(3,3 , 2,1 )= wor1(21)  !  AyCz
      x(1,3 , 3,1 )= wor1(22)  !  AzCx
      x(2,3 , 3,1 )= wor1(23)  !  AzCy
      x(3,3 , 3,1 )= wor1(24)  !  AzCz
c
c--------------------------------------------------------
c     call x_ab(x,1,4, wor1,25)  !   block AD: wor1(25-33)
c
      x(1,1 , 1,4 )= wor1(25)  !  AxDx
      x(1,1 , 2,4 )= wor1(26)  !  AxDy
      x(1,1 , 3,4 )= wor1(27)  !  AxDz
      x(2,1 , 1,4 )= wor1(28)  !  AyDx
      x(2,1 , 2,4 )= wor1(29)  !  AyDy
      x(2,1 , 3,4 )= wor1(30)  !  AyDz
      x(3,1 , 1,4 )= wor1(31)  !  AzDx
      x(3,1 , 2,4 )= wor1(32)  !  AzDy
      x(3,1 , 3,4 )= wor1(33)  !  AzDz
c
      x(1,4 , 1,1 )= wor1(25)  !  AxDx
      x(2,4 , 1,1 )= wor1(26)  !  AxDy
      x(3,4 , 1,1 )= wor1(27)  !  AxDz
      x(1,4 , 2,1 )= wor1(28)  !  AyDx
      x(2,4 , 2,1 )= wor1(29)  !  AyDy
      x(3,4 , 2,1 )= wor1(30)  !  AyDz
      x(1,4 , 3,1 )= wor1(31)  !  AzDx
      x(2,4 , 3,1 )= wor1(32)  !  AzDy
      x(3,4 , 3,1 )= wor1(33)  !  AzDz
c
c--------------------------------------------------------
c     call x_ab(x,2,3, wor1,40)  !   block BC: wor1(40-48)
c
      x(1,2 , 1,3 )= wor1(40)  !  BxCx
      x(1,2 , 2,3 )= wor1(41)  !  BxCy
      x(1,2 , 3,3 )= wor1(42)  !  BxCz
      x(2,2 , 1,3 )= wor1(43)  !  ByCx
      x(2,2 , 2,3 )= wor1(44)  !  ByCy
      x(2,2 , 3,3 )= wor1(45)  !  ByCz
      x(3,2 , 1,3 )= wor1(46)  !  BzCx
      x(3,2 , 2,3 )= wor1(47)  !  BzCy
      x(3,2 , 3,3 )= wor1(48)  !  BzCz
c
      x(1,3 , 1,2 )= wor1(40)  !  BxCx
      x(2,3 , 1,2 )= wor1(41)  !  BxCy
      x(3,3 , 1,2 )= wor1(42)  !  BxCz
      x(1,3 , 2,2 )= wor1(43)  !  ByCx
      x(2,3 , 2,2 )= wor1(44)  !  ByCy
      x(3,3 , 2,2 )= wor1(45)  !  ByCz
      x(1,3 , 3,2 )= wor1(46)  !  BzCx
      x(2,3 , 3,2 )= wor1(47)  !  BzCy
      x(3,3 , 3,2 )= wor1(48)  !  BzCz
c
c--------------------------------------------------------
c     call x_ab(x,2,4, wor1,49)  !   block BD: wor1(49-57)
c
      x(1,2 , 1,4 )= wor1(49)  !  BxDx
      x(1,2 , 2,4 )= wor1(50)  !  BxDy
      x(1,2 , 3,4 )= wor1(51)  !  BxDz
      x(2,2 , 1,4 )= wor1(52)  !  ByDx
      x(2,2 , 2,4 )= wor1(53)  !  ByDy
      x(2,2 , 3,4 )= wor1(54)  !  ByDz
      x(3,2 , 1,4 )= wor1(55)  !  BzDx
      x(3,2 , 2,4 )= wor1(56)  !  BzDy
      x(3,2 , 3,4 )= wor1(57)  !  BzDz
c
      x(1,4 , 1,2 )= wor1(49)  !  BxDx
      x(2,4 , 1,2 )= wor1(50)  !  BxDy
      x(3,4 , 1,2 )= wor1(51)  !  BxDz
      x(1,4 , 2,2 )= wor1(52)  !  ByDx
      x(2,4 , 2,2 )= wor1(53)  !  ByDy
      x(3,4 , 2,2 )= wor1(54)  !  ByDz
      x(1,4 , 3,2 )= wor1(55)  !  BzDx
      x(2,4 , 3,2 )= wor1(56)  !  BzDy
      x(3,4 , 3,2 )= wor1(57)  !  BzDz
c
c--------------------------------------------------------
c     call x_ab(x,3,4, wor1,64)  !   block CD: wor1(64-72)
c
      x(1,3 , 1,4 )= wor1(64)  !  CxDx
      x(1,3 , 2,4 )= wor1(65)  !  CxDy
      x(1,3 , 3,4 )= wor1(66)  !  CxDz
      x(2,3 , 1,4 )= wor1(67)  !  CyDx
      x(2,3 , 2,4 )= wor1(68)  !  CyDy
      x(2,3 , 3,4 )= wor1(69)  !  CyDz
      x(3,3 , 1,4 )= wor1(70)  !  CzDx
      x(3,3 , 2,4 )= wor1(71)  !  CzDy
      x(3,3 , 3,4 )= wor1(72)  !  CzDz
c
      x(1,4 , 1,3 )= wor1(64)  !  CxDx
      x(2,4 , 1,3 )= wor1(65)  !  CxDy
      x(3,4 , 1,3 )= wor1(66)  !  CxDz
      x(1,4 , 2,3 )= wor1(67)  !  CyDx
      x(2,4 , 2,3 )= wor1(68)  !  CyDy
      x(3,4 , 2,3 )= wor1(69)  !  CyDz
      x(1,4 , 3,3 )= wor1(70)  !  CzDx
      x(2,4 , 3,3 )= wor1(71)  !  CzDy
      x(3,4 , 3,3 )= wor1(72)  !  CzDz
c
c--------------------------------------------------------
c--------------------------------------------------------
c Atom-Diagonal Blocks : AA, BB, CC, DD
c
      call hess_aa(hess,na,iat, x,1 )
      call hess_aa(hess,na,jat, x,2 )
      call hess_aa(hess,na,kat, x,3 )
      call hess_aa(hess,na,lat, x,4 )
c
c Atom-Off-Diagonal Blocks :
c
      call hess_ab(hess,na,iat,jat, x,1,2)  ! block AB
      call hess_ab(hess,na,iat,kat, x,1,3)  ! block AC
      call hess_ab(hess,na,iat,lat, x,1,4)  ! block AD
c
      call hess_ab(hess,na,jat,kat, x,2,3)  ! block BC
      call hess_ab(hess,na,jat,lat, x,2,4)  ! block BD
c
      call hess_ab(hess,na,kat,lat, x,3,4)  ! block CD
c--------------------------------------------------------
      if(nsym.gt.0) then
         call getival('SymNuPr1',nupair)
         call hessian_symm(na,nsym,bl(nupair),iat,jat,kat,lat,
     *                     x,hess )
      endif
c--------------------------------------------------------
  150   continue
  100 continue
c
      end
cccccc
      subroutine print_crap(natom,nsym,nupair)
      implicit integer(a-z)
      integer nupair(natom,nsym)
      write(6,*) ' natom:',natom,' nsym:',nsym
      write(6,*) ' nupair array for real atoms is:'
      do i=1,natom
      write(6,*) (nupair(i,j),j=1,nsym)
      enddo
      call f_lush(6)
      return
      end
cccccc
c=====================================================================
      subroutine hess_mtrxC(buf, dens, hess, lind, thres1,
     *                      ntri,nbls, ngcd, lnijkl,
     *                      ncenter, labels, length, lgenct,
     *                      na, nsym )
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral derivatives buffer
c  dens()   - density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c OUTPUT :
c hess       () - hessian matrix
c----------------------------------------------------------------
c
c  ** COULOMB ONLY **
c
c ordering of 45 second derivative integrals :
c
c    der2_AxAx=buf(1,nbls,lnijkl,ngcd),
c    der2_AxAy=buf(2,nbls,lnijkl,ngcd),
c    der2_AxAz=buf(3,nbls,lnijkl,ngcd),
c    der2_AyAy=buf(4,nbls,lnijkl,ngcd),
c    der2_AyAz=buf(5,nbls,lnijkl,ngcd),
c    der2_AzAz=buf(6,nbls,lnijkl,ngcd),
c
c    der2_AxBx= buf(7,nbls,lnijkl,ngcd),
c    der2_AxBy= buf(8,nbls,lnijkl,ngcd),
c    der2_AxBz= buf(9,nbls,lnijkl,ngcd),
c    der2_AyBx=buf(10,nbls,lnijkl,ngcd),
c    der2_AyBy=buf(11,nbls,lnijkl,ngcd),
c    der2_AyBz=buf(12,nbls,lnijkl,ngcd),
c    der2_AzBx=buf(13,nbls,lnijkl,ngcd),
c    der2_AzBy=buf(14,nbls,lnijkl,ngcd),
c    der2_AzBz=buf(15,nbls,lnijkl,ngcd),
c
c    der2_AxCx=buf(16,nbls,lnijkl,ngcd),
c    der2_AxCy=buf(17,nbls,lnijkl,ngcd),
c    der2_AxCz=buf(18,nbls,lnijkl,ngcd),
c    der2_AyCx=buf(19,nbls,lnijkl,ngcd),
c    der2_AyCy=buf(20,nbls,lnijkl,ngcd),
c    der2_AyCz=buf(21,nbls,lnijkl,ngcd),
c    der2_AzCx=buf(22,nbls,lnijkl,ngcd),
c    der2_AzCy=buf(23,nbls,lnijkl,ngcd),
c    der2_AzCz=buf(24,nbls,lnijkl,ngcd),
c
c    der2_BxBx=buf(25,nbls,lnijkl,ngcd),
c    der2_BxBy=buf(26,nbls,lnijkl,ngcd),
c    der2_BxBz=buf(27,nbls,lnijkl,ngcd),
c    der2_ByBy=buf(28,nbls,lnijkl,ngcd),
c    der2_ByBz=buf(29,nbls,lnijkl,ngcd),
c    der2_BzBz=buf(30,nbls,lnijkl,ngcd),
c
c    der2_BxCx=buf(31,nbls,lnijkl,ngcd),
c    der2_BxCy=buf(32,nbls,lnijkl,ngcd),
c    der2_BxCz=buf(33,nbls,lnijkl,ngcd),
c    der2_ByCx=buf(34,nbls,lnijkl,ngcd),
c    der2_ByCy=buf(35,nbls,lnijkl,ngcd),
c    der2_ByCz=buf(36,nbls,lnijkl,ngcd),
c    der2_BzCx=buf(37,nbls,lnijkl,ngcd),
c    der2_BzCy=buf(38,nbls,lnijkl,ngcd),
c    der2_BzCz=buf(39,nbls,lnijkl,ngcd),
c
c    der2_CxCx=buf(40,nbls,lnijkl,ngcd),
c    der2_CxCy=buf(41,nbls,lnijkl,ngcd),
c    der2_CxCz=buf(42,nbls,lnijkl,ngcd),
c    der2_CyCy=buf(43,nbls,lnijkl,ngcd),
c    der2_CyCz=buf(44,nbls,lnijkl,ngcd),
c    der2_CzCz=buf(45,nbls,lnijkl,ngcd),
c----------------------------------------------------------------
      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(45,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension hess(3,na, 3,na)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension work(78)
      dimension x(3,4,3,4)
      data zero /0.d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
          do 150 iqu=1,ngcd
          icff=labels(1,iqu,ijklp)
          jcff=labels(2,iqu,ijklp)
          kcff=labels(3,iqu,ijklp)
          lcff=labels(4,iqu,ijklp)
c
          iat=ncenter(icff+1)
          jat=ncenter(jcff+1)
          kat=ncenter(kcff+1)
          lat=ncenter(lcff+1)
c
          call zeroit(x,144)  ! 144=12*12
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
             integ=0
             do 200 iii=1,ilen
             icf=icff+iii
             ii=lind(icf)
                do 250 jjj=1,jlen
                jcf=jcff+jjj
                jj=lind(jcf)
                ijf=ii+jcf
                if(jcf.gt.icf) ijf=jj+icf
                dij4=dens(ijf)*4.d0
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      klf=kk+lcf
                      if(lcf.gt.kcf) klf=ll+kcf
c
                      dkl=dens(klf)
c
                      dijkl=dij4*dkl
c
                      integ=integ+1
c--------------------------------------------------------
          if(abs(dijkl).lt.thres1) go to 350
c--------------------------------------------------------
        call dscal(45,dijkl,buf(1,ijklp,integ,iqu),1)
c
        call make_78_from_45(buf(1,ijklp,integ,iqu),work)
c--------------------------------------------------------
c in work 78 derivatives ordered as :
c
c   all 10 blocks of sec.der. made out of 6 blocks:
c
c          AA AB AC AD                AA AB AC
c             BB BC BD      from         BB BC
c                CC CD                      CC
c                   DD
c      1-6, 7-15,16-24,25-33         1-6, 7-15,16-24
c          34-39,40-48,49-57             25-30,31-39
c                58-63,64-72                   40-45
c                      73-78
c--------------------------------------------------------
c instead of x_aa :
c
c     call x_aa(x,1,work,1)      !   Block AA: work(1-6)
c
      x(1,1 , 1,1 )=x(1,1 , 1,1 ) + work(1)   !  AxAx
      x(1,1 , 2,1 )=x(1,1 , 2,1 ) + work(2)   !  AxAy
      x(1,1 , 3,1 )=x(1,1 , 3,1 ) + work(3)   !  AxAz
      x(2,1 , 2,1 )=x(2,1 , 2,1 ) + work(4)   !  AyAy
      x(2,1 , 3,1 )=x(2,1 , 3,1 ) + work(5)   !  AyAz
      x(3,1 , 3,1 )=x(3,1 , 3,1 ) + work(6)   !  AzAz
c
c--------------------------------------------------------
c     call x_aa(x,2,work,34)     !   block BB: work(34-39)
c
      x(1,2 , 1,2 )=x(1,2 , 1,2 ) + work(34)   !  BxBx
      x(1,2 , 2,2 )=x(1,2 , 2,2 ) + work(35)   !  BxBy
      x(1,2 , 3,2 )=x(1,2 , 3,2 ) + work(36)   !  BxBz
      x(2,2 , 2,2 )=x(2,2 , 2,2 ) + work(37)   !  ByBy
      x(2,2 , 3,2 )=x(2,2 , 3,2 ) + work(38)   !  ByBz
      x(3,2 , 3,2 )=x(3,2 , 3,2 ) + work(39)   !  BzBz
c
c--------------------------------------------------------
c     call x_aa(x,3,work,58)     !   block CC: work(58-63)
c
      x(1,3 , 1,3 )=x(1,3 , 1,3 ) + work(58)   !  CxCx
      x(1,3 , 2,3 )=x(1,3 , 2,3 ) + work(59)   !  CxCy
      x(1,3 , 3,3 )=x(1,3 , 3,3 ) + work(60)   !  CxCz
      x(2,3 , 2,3 )=x(2,3 , 2,3 ) + work(61)   !  CyCy
      x(2,3 , 3,3 )=x(2,3 , 3,3 ) + work(62)   !  CyCz
      x(3,3 , 3,3 )=x(3,3 , 3,3 ) + work(63)   !  CzCz

c--------------------------------------------------------
c     call x_aa(x,4,work,73)     !   block DD: work(73-78)
c
      x(1,4 , 1,4 )=x(1,4 , 1,4 ) + work(73)   !  DxDx
      x(1,4 , 2,4 )=x(1,4 , 2,4 ) + work(74)   !  DxDy
      x(1,4 , 3,4 )=x(1,4 , 3,4 ) + work(75)   !  DxDz
      x(2,4 , 2,4 )=x(2,4 , 2,4 ) + work(76)   !  DyDy
      x(2,4 , 3,4 )=x(2,4 , 3,4 ) + work(77)   !  DyDz
      x(3,4 , 3,4 )=x(3,4 , 3,4 ) + work(78)   !  DzDz
c
c--------------------------------------------------------
c     call x_ab(x,1,2, work,7)   !   block AB: work(7-15)
c
      x(1,1 , 1,2 )=x(1,1 , 1,2 ) + work(7)   !  AxBx
      x(1,1 , 2,2 )=x(1,1 , 2,2 ) + work(8)   !  AxBy
      x(1,1 , 3,2 )=x(1,1 , 3,2 ) + work(9)   !  AxBz
      x(2,1 , 1,2 )=x(2,1 , 1,2 ) + work(10)  !  AyBx
      x(2,1 , 2,2 )=x(2,1 , 2,2 ) + work(11)  !  AyBy
      x(2,1 , 3,2 )=x(2,1 , 3,2 ) + work(12)  !  AyBz
      x(3,1 , 1,2 )=x(3,1 , 1,2 ) + work(13)  !  AzBx
      x(3,1 , 2,2 )=x(3,1 , 2,2 ) + work(14)  !  AzBy
      x(3,1 , 3,2 )=x(3,1 , 3,2 ) + work(15)  !  AzBz
c
      x(1,2 , 1,1 )=x(1,2 , 1,1 ) + work(7)   !  AxBx
      x(2,2 , 1,1 )=x(2,2 , 1,1 ) + work(8)   !  AxBy
      x(3,2 , 1,1 )=x(3,2 , 1,1 ) + work(9)   !  AxBz
      x(1,2 , 2,1 )=x(1,2 , 2,1 ) + work(10)  !  AyBx
      x(2,2 , 2,1 )=x(2,2 , 2,1 ) + work(11)  !  AyBy
      x(3,2 , 2,1 )=x(3,2 , 2,1 ) + work(12)  !  AyBz
      x(1,2 , 3,1 )=x(1,2 , 3,1 ) + work(13)  !  AzBx
      x(2,2 , 3,1 )=x(2,2 , 3,1 ) + work(14)  !  AzBy
      x(3,2 , 3,1 )=x(3,2 , 3,1 ) + work(15)  !  AzBz
c--------------------------------------------------------
c     call x_ab(x,1,3, work,16)  !   block AC: work(16-24)
c
      x(1,1 , 1,3 )=x(1,1 , 1,3 ) + work(16)  !  AxCx
      x(1,1 , 2,3 )=x(1,1 , 2,3 ) + work(17)  !  AxCy
      x(1,1 , 3,3 )=x(1,1 , 3,3 ) + work(18)  !  AxCz
      x(2,1 , 1,3 )=x(2,1 , 1,3 ) + work(19)  !  AyCx
      x(2,1 , 2,3 )=x(2,1 , 2,3 ) + work(20)  !  AyCy
      x(2,1 , 3,3 )=x(2,1 , 3,3 ) + work(21)  !  AyCz
      x(3,1 , 1,3 )=x(3,1 , 1,3 ) + work(22)  !  AzCx
      x(3,1 , 2,3 )=x(3,1 , 2,3 ) + work(23)  !  AzCy
      x(3,1 , 3,3 )=x(3,1 , 3,3 ) + work(24)  !  AzCz
c
      x(1,3 , 1,1 )=x(1,3 , 1,1 ) + work(16)  !  AxCx
      x(2,3 , 1,1 )=x(2,3 , 1,1 ) + work(17)  !  AxCy
      x(3,3 , 1,1 )=x(3,3 , 1,1 ) + work(18)  !  AxCz
      x(1,3 , 2,1 )=x(1,3 , 2,1 ) + work(19)  !  AyCx
      x(2,3 , 2,1 )=x(2,3 , 2,1 ) + work(20)  !  AyCy
      x(3,3 , 2,1 )=x(3,3 , 2,1 ) + work(21)  !  AyCz
      x(1,3 , 3,1 )=x(1,3 , 3,1 ) + work(22)  !  AzCx
      x(2,3 , 3,1 )=x(2,3 , 3,1 ) + work(23)  !  AzCy
      x(3,3 , 3,1 )=x(3,3 , 3,1 ) + work(24)  !  AzCz
c
c--------------------------------------------------------
c     call x_ab(x,1,4, work,25)  !   block AD: work(25-33)
c
      x(1,1 , 1,4 )=x(1,1 , 1,4 ) + work(25)  !  AxDx
      x(1,1 , 2,4 )=x(1,1 , 2,4 ) + work(26)  !  AxDy
      x(1,1 , 3,4 )=x(1,1 , 3,4 ) + work(27)  !  AxDz
      x(2,1 , 1,4 )=x(2,1 , 1,4 ) + work(28)  !  AyDx
      x(2,1 , 2,4 )=x(2,1 , 2,4 ) + work(29)  !  AyDy
      x(2,1 , 3,4 )=x(2,1 , 3,4 ) + work(30)  !  AyDz
      x(3,1 , 1,4 )=x(3,1 , 1,4 ) + work(31)  !  AzDx
      x(3,1 , 2,4 )=x(3,1 , 2,4 ) + work(32)  !  AzDy
      x(3,1 , 3,4 )=x(3,1 , 3,4 ) + work(33)  !  AzDz
c
      x(1,4 , 1,1 )=x(1,4 , 1,1 ) + work(25)  !  AxDx
      x(2,4 , 1,1 )=x(2,4 , 1,1 ) + work(26)  !  AxDy
      x(3,4 , 1,1 )=x(3,4 , 1,1 ) + work(27)  !  AxDz
      x(1,4 , 2,1 )=x(1,4 , 2,1 ) + work(28)  !  AyDx
      x(2,4 , 2,1 )=x(2,4 , 2,1 ) + work(29)  !  AyDy
      x(3,4 , 2,1 )=x(3,4 , 2,1 ) + work(30)  !  AyDz
      x(1,4 , 3,1 )=x(1,4 , 3,1 ) + work(31)  !  AzDx
      x(2,4 , 3,1 )=x(2,4 , 3,1 ) + work(32)  !  AzDy
      x(3,4 , 3,1 )=x(3,4 , 3,1 ) + work(33)  !  AzDz
c
c--------------------------------------------------------
c     call x_ab(x,2,3, work,40)  !   block BC: work(40-48)
c
      x(1,2 , 1,3 )=x(1,2 , 1,3 ) + work(40)  !  BxCx
      x(1,2 , 2,3 )=x(1,2 , 2,3 ) + work(41)  !  BxCy
      x(1,2 , 3,3 )=x(1,2 , 3,3 ) + work(42)  !  BxCz
      x(2,2 , 1,3 )=x(2,2 , 1,3 ) + work(43)  !  ByCx
      x(2,2 , 2,3 )=x(2,2 , 2,3 ) + work(44)  !  ByCy
      x(2,2 , 3,3 )=x(2,2 , 3,3 ) + work(45)  !  ByCz
      x(3,2 , 1,3 )=x(3,2 , 1,3 ) + work(46)  !  BzCx
      x(3,2 , 2,3 )=x(3,2 , 2,3 ) + work(47)  !  BzCy
      x(3,2 , 3,3 )=x(3,2 , 3,3 ) + work(48)  !  BzCz
c
      x(1,3 , 1,2 )=x(1,3 , 1,2 ) + work(40)  !  BxCx
      x(2,3 , 1,2 )=x(2,3 , 1,2 ) + work(41)  !  BxCy
      x(3,3 , 1,2 )=x(3,3 , 1,2 ) + work(42)  !  BxCz
      x(1,3 , 2,2 )=x(1,3 , 2,2 ) + work(43)  !  ByCx
      x(2,3 , 2,2 )=x(2,3 , 2,2 ) + work(44)  !  ByCy
      x(3,3 , 2,2 )=x(3,3 , 2,2 ) + work(45)  !  ByCz
      x(1,3 , 3,2 )=x(1,3 , 3,2 ) + work(46)  !  BzCx
      x(2,3 , 3,2 )=x(2,3 , 3,2 ) + work(47)  !  BzCy
      x(3,3 , 3,2 )=x(3,3 , 3,2 ) + work(48)  !  BzCz
c
c--------------------------------------------------------
c     call x_ab(x,2,4, work,49)  !   block BD: work(49-57)
c
      x(1,2 , 1,4 )=x(1,2 , 1,4 ) + work(49)  !  BxDx
      x(1,2 , 2,4 )=x(1,2 , 2,4 ) + work(50)  !  BxDy
      x(1,2 , 3,4 )=x(1,2 , 3,4 ) + work(51)  !  BxDz
      x(2,2 , 1,4 )=x(2,2 , 1,4 ) + work(52)  !  ByDx
      x(2,2 , 2,4 )=x(2,2 , 2,4 ) + work(53)  !  ByDy
      x(2,2 , 3,4 )=x(2,2 , 3,4 ) + work(54)  !  ByDz
      x(3,2 , 1,4 )=x(3,2 , 1,4 ) + work(55)  !  BzDx
      x(3,2 , 2,4 )=x(3,2 , 2,4 ) + work(56)  !  BzDy
      x(3,2 , 3,4 )=x(3,2 , 3,4 ) + work(57)  !  BzDz
c
      x(1,4 , 1,2 )=x(1,4 , 1,2 ) + work(49)  !  BxDx
      x(2,4 , 1,2 )=x(2,4 , 1,2 ) + work(50)  !  BxDy
      x(3,4 , 1,2 )=x(3,4 , 1,2 ) + work(51)  !  BxDz
      x(1,4 , 2,2 )=x(1,4 , 2,2 ) + work(52)  !  ByDx
      x(2,4 , 2,2 )=x(2,4 , 2,2 ) + work(53)  !  ByDy
      x(3,4 , 2,2 )=x(3,4 , 2,2 ) + work(54)  !  ByDz
      x(1,4 , 3,2 )=x(1,4 , 3,2 ) + work(55)  !  BzDx
      x(2,4 , 3,2 )=x(2,4 , 3,2 ) + work(56)  !  BzDy
      x(3,4 , 3,2 )=x(3,4 , 3,2 ) + work(57)  !  BzDz
c
c--------------------------------------------------------
c     call x_ab(x,3,4, work,64)  !   block CD: work(64-72)
c
      x(1,3 , 1,4 )=x(1,3 , 1,4 ) + work(64)  !  CxDx
      x(1,3 , 2,4 )=x(1,3 , 2,4 ) + work(65)  !  CxDy
      x(1,3 , 3,4 )=x(1,3 , 3,4 ) + work(66)  !  CxDz
      x(2,3 , 1,4 )=x(2,3 , 1,4 ) + work(67)  !  CyDx
      x(2,3 , 2,4 )=x(2,3 , 2,4 ) + work(68)  !  CyDy
      x(2,3 , 3,4 )=x(2,3 , 3,4 ) + work(69)  !  CyDz
      x(3,3 , 1,4 )=x(3,3 , 1,4 ) + work(70)  !  CzDx
      x(3,3 , 2,4 )=x(3,3 , 2,4 ) + work(71)  !  CzDy
      x(3,3 , 3,4 )=x(3,3 , 3,4 ) + work(72)  !  CzDz
c
      x(1,4 , 1,3 )=x(1,4 , 1,3 ) + work(64)  !  CxDx
      x(2,4 , 1,3 )=x(2,4 , 1,3 ) + work(65)  !  CxDy
      x(3,4 , 1,3 )=x(3,4 , 1,3 ) + work(66)  !  CxDz
      x(1,4 , 2,3 )=x(1,4 , 2,3 ) + work(67)  !  CyDx
      x(2,4 , 2,3 )=x(2,4 , 2,3 ) + work(68)  !  CyDy
      x(3,4 , 2,3 )=x(3,4 , 2,3 ) + work(69)  !  CyDz
      x(1,4 , 3,3 )=x(1,4 , 3,3 ) + work(70)  !  CzDx
      x(2,4 , 3,3 )=x(2,4 , 3,3 ) + work(71)  !  CzDy
      x(3,4 , 3,3 )=x(3,4 , 3,3 ) + work(72)  !  CzDz
c
c--------------------------------------------------------
ccccc call x_aa(x,1,work,1)      !   Block AA: work(1-6)
c     call x_ab(x,1,2, work,7)   !   block AB: work(7-15)
c     call x_ab(x,1,3, work,16)  !   block AC: work(16-24)
c     call x_ab(x,1,4, work,25)  !   block AD: work(25-33)
ccccc call x_aa(x,2,work,34)     !   block BB: work(34-39)
c     call x_ab(x,2,3, work,40)  !   block BC: work(40-48)
c     call x_ab(x,2,4, work,49)  !   block BD: work(49-57)
ccccc call x_aa(x,3,work,58)     !   block CC: work(58-63)
c     call x_ab(x,3,4, work,64)  !   block CD: work(64-72)
ccccc call x_aa(x,4,work,73)     !   block DD: work(73-78)
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
c Atom-Diagonal Blocks : AA, BB, CC, DD
c
      call hess_aa(hess,na,iat, x,1 )
      call hess_aa(hess,na,jat, x,2 )
      call hess_aa(hess,na,kat, x,3 )
      call hess_aa(hess,na,lat, x,4 )
c
c Atom-Off-Diagonal Blocks :
c
      call hess_ab(hess,na,iat,jat, x,1,2)  ! block AB
      call hess_ab(hess,na,iat,kat, x,1,3)  ! block AC
      call hess_ab(hess,na,iat,lat, x,1,4)  ! block AD
c
      call hess_ab(hess,na,jat,kat, x,2,3)  ! block BC
      call hess_ab(hess,na,jat,lat, x,2,4)  ! block BD
c
      call hess_ab(hess,na,kat,lat, x,3,4)  ! block CD
c--------------------------------------------------------
      if(nsym.gt.0) then
         call getival('SymNuPr1',nupair)
         call hessian_symm(na,nsym,bl(nupair),iat,jat,kat,lat,
     *                     x,hess )
      endif
c--------------------------------------------------------
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine hess_aa(hess,na,iat, x,ia )
c     make atom-diagonal blocks for iat from ia .
c
      implicit real*8 (a-h,o-z)
      dimension hess(3,na, 3,na)
      dimension x(3,4,3,4)
c
      hess(1,iat, 1,iat)=hess(1,iat, 1,iat)+x(1, ia , 1, ia )
      hess(1,iat, 2,iat)=hess(1,iat, 2,iat)+x(1, ia , 2, ia )
      hess(1,iat, 3,iat)=hess(1,iat, 3,iat)+x(1, ia , 3, ia )
      hess(2,iat, 2,iat)=hess(2,iat, 2,iat)+x(2, ia , 2, ia )
      hess(2,iat, 3,iat)=hess(2,iat, 3,iat)+x(2, ia , 3, ia )
      hess(3,iat, 3,iat)=hess(3,iat, 3,iat)+x(3, ia , 3, ia )
c
      end
c=====================================================================
      subroutine hess_ab(hess,na,iat,jat, x,ia,ja)
      implicit real*8 (a-h,o-z)
      dimension hess(3,na, 3,na)
      dimension x(3,4,3,4)
c--------------------------------------------------------
c block AB
c 1 :
      if(iat.lt.jat) then
      hess(1,iat, 1,jat)=hess(1,iat, 1,jat) + x(1,ia,1,ja) !  AxBx
      hess(1,iat, 2,jat)=hess(1,iat, 2,jat) + x(1,ia,2,ja) !  AxBy
      hess(1,iat, 3,jat)=hess(1,iat, 3,jat) + x(1,ia,3,ja) !  AxBz
c
      hess(2,iat, 1,jat)=hess(2,iat, 1,jat) + x(2,ia,1,ja) !  AyBx
      hess(2,iat, 2,jat)=hess(2,iat, 2,jat) + x(2,ia,2,ja) !  AyBy
      hess(2,iat, 3,jat)=hess(2,iat, 3,jat) + x(2,ia,3,ja) !  AyBz
c
      hess(3,iat, 1,jat)=hess(3,iat, 1,jat) + x(3,ia,1,ja) !  AzBx
      hess(3,iat, 2,jat)=hess(3,iat, 2,jat) + x(3,ia,2,ja) !  AzBy
      hess(3,iat, 3,jat)=hess(3,iat, 3,jat) + x(3,ia,3,ja) !  AzBz
      endif
c 2 :
      if(iat.gt.jat) then
      hess(1,jat, 1,iat)=hess(1,jat, 1,iat) + x(1,ja,1,ia) !  AxBx
      hess(2,jat, 1,iat)=hess(2,jat, 1,iat) + x(2,ja,1,ia) !  AxBy
      hess(3,jat, 1,iat)=hess(3,jat, 1,iat) + x(3,ja,1,ia) !  AxBz
c
      hess(1,jat, 2,iat)=hess(1,jat, 2,iat) + x(1,ja,2,ia) !  AyBx
      hess(2,jat, 2,iat)=hess(2,jat, 2,iat) + x(2,ja,2,ia) !  AyBy
      hess(3,jat, 2,iat)=hess(3,jat, 2,iat) + x(3,ja,2,ia) !  AyBz
c
      hess(1,jat, 3,iat)=hess(1,jat, 3,iat) + x(1,ja,3,ia) !  AzBx
      hess(2,jat, 3,iat)=hess(2,jat, 3,iat) + x(2,ja,3,ia) !  AzBy
      hess(3,jat, 3,iat)=hess(3,jat, 3,iat) + x(3,ja,3,ia) !  AzBz
      endif
c
c 3:
      if(iat.eq.jat) then
      hess(1,iat, 1,jat)=hess(1,iat, 1,jat) + x(1,ia,1,ja) !  AxBx
      hess(1,iat, 1,jat)=hess(1,iat, 1,jat) + x(1,ja,1,ia) !  AxBx
c
      hess(1,iat, 2,jat)=hess(1,iat, 2,jat) + x(1,ia,2,ja) !  AxBy
     *                                      + x(1,ja,2,ia) !  AxBy
      hess(1,iat, 3,jat)=hess(1,iat, 3,jat) + x(1,ia,3,ja) !  AxBz
     *                                      + x(1,ja,3,ia) !  AxBz
      hess(2,iat, 2,jat)=hess(2,iat, 2,jat) + x(2,ia,2,ja) !  AyBy
      hess(2,iat, 2,jat)=hess(2,iat, 2,jat) + x(2,ja,2,ia) !  AyBy
c
      hess(2,iat, 3,jat)=hess(2,iat, 3,jat) + x(2,ia,3,ja) !  AyBz
     *                                      + x(2,ja,3,ia) !  AyBz
c
      hess(3,iat, 3,jat)=hess(3,iat, 3,jat) + x(3,ia,3,ja) !  AzBz
      hess(3,iat, 3,jat)=hess(3,iat, 3,jat) + x(3,ja,3,ia) !  AzBz
      endif
c--------------------------------------------------------
      end
c=====================================================================
      subroutine x_aa(x,ia,work,i)
      implicit real*8 (a-h,o-z)
      dimension x(3,4,3,4) , work(78)
c
      i1=i
      i2=i1+1
      i3=i2+1
      i4=i3+1
      i5=i4+1
      i6=i5+1
c
      x(1,ia , 1,ia )=x(1,ia , 1,ia ) + work(i1)   !  AxAx
      x(1,ia , 2,ia )=x(1,ia , 2,ia ) + work(i2)   !  AxAy
      x(1,ia , 3,ia )=x(1,ia , 3,ia ) + work(i3)   !  AxAz
c
      x(2,ia , 2,ia )=x(2,ia , 2,ia ) + work(i4)   !  AyAy
      x(2,ia , 3,ia )=x(2,ia , 3,ia ) + work(i5)   !  AyAz
c
      x(3,ia , 3,ia )=x(3,ia , 3,ia ) + work(i6)   !  AzAz
c
      end
c=====================================================================
      subroutine x_ab(x,i1,i2, work,j)
      implicit real*8 (a-h,o-z)
      dimension x(3,4,3,4) , work(78)
c
      j1=j
      j2=j1+1
      j3=j2+1
      j4=j3+1
      j5=j4+1
      j6=j5+1
      j7=j6+1
      j8=j7+1
      j9=j8+1
c
      x(1,i1 , 1,i2 )=x(1,i1 , 1,i2 ) + work(j1)   !  AxBx
      x(1,i1 , 2,i2 )=x(1,i1 , 2,i2 ) + work(j2)   !  AxBy
      x(1,i1 , 3,i2 )=x(1,i1 , 3,i2 ) + work(j3)   !  AxBz
c
      x(2,i1 , 1,i2 )=x(2,i1 , 1,i2 ) + work(j4)  !  AyBx
      x(2,i1 , 2,i2 )=x(2,i1 , 2,i2 ) + work(j5)  !  AyBy
      x(2,i1 , 3,i2 )=x(2,i1 , 3,i2 ) + work(j6)  !  AyBz
c
      x(3,i1 , 1,i2 )=x(3,i1 , 1,i2 ) + work(j7)  !  AzBx
      x(3,i1 , 2,i2 )=x(3,i1 , 2,i2 ) + work(j8)  !  AzBy
      x(3,i1 , 3,i2 )=x(3,i1 , 3,i2 ) + work(j9)  !  AzBz
c
      x(1,i2 , 1,i1 )=x(1,i2 , 1,i1 ) + work(j1)   !  AxBx
      x(2,i2 , 1,i1 )=x(2,i2 , 1,i1 ) + work(j2)   !  AxBy
      x(3,i2 , 1,i1 )=x(3,i2 , 1,i1 ) + work(j3)   !  AxBz
c
      x(1,i2 , 2,i1 )=x(1,i2 , 2,i1 ) + work(j4)  !  AyBx
      x(2,i2 , 2,i1 )=x(2,i2 , 2,i1 ) + work(j5)  !  AyBy
      x(3,i2 , 2,i1 )=x(3,i2 , 2,i1 ) + work(j6)  !  AyBz
c
      x(1,i2 , 3,i1 )=x(1,i2 , 3,i1 ) + work(j7)  !  AzBx
      x(2,i2 , 3,i1 )=x(2,i2 , 3,i1 ) + work(j8)  !  AzBy
      x(3,i2 , 3,i1 )=x(3,i2 , 3,i1 ) + work(j9)  !  AzBz
c
      end
c=====================================================================
      subroutine hessian_symm(na,nsym,nupair,iat,jat,kat,lat,x,hess)
c--------------------------------------------------------
c symmetrize two-electron part of the hessian matrix :
c
c Input/output : hess
c--------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common /negxyz/ negx(7),negy(7),negz(7)
c
      dimension nupair(na,nsym)           ! (natoms,nsym)
      dimension ngxyz(3,7)
      dimension hess(3,na,3,na)
      dimension x(3,4,3,4)
c
      do ns=1,nsym
         ngxyz(1,ns)=negx(ns)
         ngxyz(2,ns)=negy(ns)
         ngxyz(3,ns)=negz(ns)
      enddo
c
      do ns=1,nsym
         iat1=nupair(iat,ns)
         jat1=nupair(jat,ns)
         kat1=nupair(kat,ns)
         lat1=nupair(lat,ns)
c
         call hess_aax(hess,na,iat1, x,1 ,ngxyz(1,ns) )
         call hess_aax(hess,na,jat1, x,2 ,ngxyz(1,ns) )
         call hess_aax(hess,na,kat1, x,3 ,ngxyz(1,ns) )
         call hess_aax(hess,na,lat1, x,4 ,ngxyz(1,ns) )
c
         call hess_abx(hess,na,iat1,jat1, x,1,2,ngxyz(1,ns)) ! block AB
         call hess_abx(hess,na,iat1,kat1, x,1,3,ngxyz(1,ns)) ! block AC
         call hess_abx(hess,na,iat1,lat1, x,1,4,ngxyz(1,ns)) ! block AD
         call hess_abx(hess,na,jat1,kat1, x,2,3,ngxyz(1,ns)) ! block BC
         call hess_abx(hess,na,jat1,lat1, x,2,4,ngxyz(1,ns)) ! block BD
         call hess_abx(hess,na,kat1,lat1, x,3,4,ngxyz(1,ns)) ! block CD
      enddo
c
      end
c==============================================================
      subroutine hess_aax(hess,na,iat, x,ia ,ng)
      implicit real*8 (a-h,o-z)
      dimension ng(3)
      dimension hess(3,na, 3,na)
      dimension x(3,4,3,4)
c
      hess(1,iat, 1,iat)=hess(1,iat, 1,iat)+x(1,ia , 1,ia)*ng(1)*ng(1)
      hess(1,iat, 2,iat)=hess(1,iat, 2,iat)+x(1,ia , 2,ia)*ng(1)*ng(2)
      hess(1,iat, 3,iat)=hess(1,iat, 3,iat)+x(1,ia , 3,ia)*ng(1)*ng(3)
      hess(2,iat, 2,iat)=hess(2,iat, 2,iat)+x(2,ia , 2,ia)*ng(2)*ng(2)
      hess(2,iat, 3,iat)=hess(2,iat, 3,iat)+x(2,ia , 3,ia)*ng(2)*ng(3)
      hess(3,iat, 3,iat)=hess(3,iat, 3,iat)+x(3,ia , 3,ia)*ng(3)*ng(3)
c
      end
c=====================================================================
      subroutine hess_abx(hess,na,iat,jat, x,ia,ja, ng)
      implicit real*8 (a-h,o-z)
      dimension ng(3)
      dimension hess(3,na, 3,na)
      dimension x(3,4,3,4)
c--------------------------------------------------------
c block AB
c 1 :
      if(iat.lt.jat) then
      hess(1,iat, 1,jat)=hess(1,iat, 1,jat)+x(1,ia,1,ja)*ng(1)*ng(1)
      hess(1,iat, 2,jat)=hess(1,iat, 2,jat)+x(1,ia,2,ja)*ng(1)*ng(2)
      hess(1,iat, 3,jat)=hess(1,iat, 3,jat)+x(1,ia,3,ja)*ng(1)*ng(3)
c
      hess(2,iat, 1,jat)=hess(2,iat, 1,jat)+x(2,ia,1,ja)*ng(2)*ng(1)
      hess(2,iat, 2,jat)=hess(2,iat, 2,jat)+x(2,ia,2,ja)*ng(2)*ng(2)
      hess(2,iat, 3,jat)=hess(2,iat, 3,jat)+x(2,ia,3,ja)*ng(2)*ng(3)
c
      hess(3,iat, 1,jat)=hess(3,iat, 1,jat)+x(3,ia,1,ja)*ng(3)*ng(1)
      hess(3,iat, 2,jat)=hess(3,iat, 2,jat)+x(3,ia,2,ja)*ng(3)*ng(2)
      hess(3,iat, 3,jat)=hess(3,iat, 3,jat)+x(3,ia,3,ja)*ng(3)*ng(3)
      endif
c 2 :
      if(iat.gt.jat) then
      hess(1,jat, 1,iat)=hess(1,jat, 1,iat)+x(1,ja,1,ia)*ng(1)*ng(1)
      hess(2,jat, 1,iat)=hess(2,jat, 1,iat)+x(2,ja,1,ia)*ng(2)*ng(1)
      hess(3,jat, 1,iat)=hess(3,jat, 1,iat)+x(3,ja,1,ia)*ng(3)*ng(1)
c
      hess(1,jat, 2,iat)=hess(1,jat, 2,iat)+x(1,ja,2,ia)*ng(1)*ng(2)
      hess(2,jat, 2,iat)=hess(2,jat, 2,iat)+x(2,ja,2,ia)*ng(2)*ng(2)
      hess(3,jat, 2,iat)=hess(3,jat, 2,iat)+x(3,ja,2,ia)*ng(3)*ng(2)
c
      hess(1,jat, 3,iat)=hess(1,jat, 3,iat)+x(1,ja,3,ia)*ng(1)*ng(3)
      hess(2,jat, 3,iat)=hess(2,jat, 3,iat)+x(2,ja,3,ia)*ng(2)*ng(3)
      hess(3,jat, 3,iat)=hess(3,jat, 3,iat)+x(3,ja,3,ia)*ng(3)*ng(3)
      endif
c
c 3:
      if(iat.eq.jat) then
      hess(1,iat, 1,jat)=hess(1,iat, 1,jat)+x(1,ia,1,ja)*ng(1)*ng(1)
      hess(1,iat, 1,jat)=hess(1,iat, 1,jat)+x(1,ja,1,ia)*ng(1)*ng(1)
c
      hess(1,iat, 2,jat)=hess(1,iat, 2,jat)+x(1,ia,2,ja)*ng(1)*ng(2)
     *                                     +x(1,ja,2,ia)*ng(1)*ng(2)
      hess(1,iat, 3,jat)=hess(1,iat, 3,jat)+x(1,ia,3,ja)*ng(1)*ng(3)
     *                                     +x(1,ja,3,ia)*ng(1)*ng(3)
      hess(2,iat, 2,jat)=hess(2,iat, 2,jat)+x(2,ia,2,ja)*ng(2)*ng(2)
      hess(2,iat, 2,jat)=hess(2,iat, 2,jat)+x(2,ja,2,ia)*ng(2)*ng(2)
c
      hess(2,iat, 3,jat)=hess(2,iat, 3,jat)+x(2,ia,3,ja)*ng(2)*ng(3)
     *                                     +x(2,ja,3,ia)*ng(2)*ng(3)
      hess(3,iat, 3,jat)=hess(3,iat, 3,jat)+x(3,ia,3,ja)*ng(3)*ng(3)
      hess(3,iat, 3,jat)=hess(3,iat, 3,jat)+x(3,ja,3,ia)*ng(3)*ng(3)
      endif
c--------------------------------------------------------
      end
c=====================================================================
      subroutine do_intder2(idft,ax,rhf,istart,istop,
     *                      bl,inx,thres1,densp,labels,
     *                    dens,denB,hessian,lind,ntri,ncenter,na,nsym)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
      logical rhf,moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
      common /runtype/ scftype,where
c
      parameter (zero=0.0d0)
c---------------------------------------------------------------------
      dimension inx(12,*)
      dimension labels(*)
      dimension lind(*)
      dimension ncenter(*)
      dimension bl(*)
      dimension densp(*)
      dimension dens(*)
      dimension denB(*)
      dimension hessian(*)
c----------------------------------------------------------------
      call getival('printh',nprint)
c----------------------------------------------------------------
c calculate the hessian integrals
c
      where='hess'
      screen='forc'     ! screening with 2-particle density
c
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c loop over blocks of contracted shell quartets
c
      DO ISUPBL=ISTART,ISTOP
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,densp,where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
c
c check if requested super-block was in a list :
c
        if(stopnow) go to 1111
c
c integrals arrived in bl(ibuffz); check if they are :
c bl(ibuffz) contains now BOTH the gradient and the hessian
c integral derivatives . Pointers to them are :
c   ider1=ibuffz
c   ider2=ider1 + 9*nbls *lnijkl*ngcd
c   which is    + 9*nblsiz*nintez*ngctoz
c
      ider1=ibuffz
      ider2=ider1 + 9*nblsiz*nintez*ngctoz
c
c  nintez is the size of a given quartet.It is set up
c  to zero if there is no integrals
c
      IF(nintez.gt.0) THEN
        call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
c
c calculate two-el. contribution to hessian and derivative fock matrices :
c
        IF(idft.EQ.0) THEN
c  ...........................
c   Standard HF
c  ...........................
          if(nprint.ge.5) then
c          print 2ed der.integ:
           call prtn_der2(bl(ider2 ),lind,
     *                    ntri,nblsiz,ngctoz,nintez,ncenter,
     *                    labels(lab1),labels(lab2),labels(lab3) )
          endif
c
          If(rhf) Then
           call hess_mtrx(bl(ider2 ),dens,hessian ,lind,thres1,
     *                    ntri,nblsiz,ngctoz,nintez,ncenter,
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    na,nsym,idft,ax)
          Else
           call hess_mtrx_uhf(bl(ider2 ),dens,denB,hessian ,lind,
     $                        thres1,ntri,nblsiz,ngctoz,nintez,
     $                        ncenter,labels(lab1),labels(lab2),
     $                        labels(lab3),na,nsym,idft,ax)
          EndIf
cc
        ELSE IF(ax.NE.Zero) THEN
c  ...........................
c   Hybrid HF-DFT (ACM)
c  ...........................
          If(rhf) Then
           call hess_mtrx(bl(ider2 ),dens,hessian ,lind,thres1,
     *                    ntri,nblsiz,ngctoz,nintez,ncenter,
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    na,nsym,idft,ax)
          Else
           call hess_mtrx_uhf(bl(ider2 ),dens,denB,hessian ,lind,
     $                        thres1,ntri,nblsiz,ngctoz,nintez,
     $                        ncenter,labels(lab1),labels(lab2),
     $                        labels(lab3),na,nsym,idft,ax)
          EndIf
cc
        ELSE
c  ...........................
c   "Pure" DFT (NO EXCHANGE)
c  ...........................
cc
          If(rhf) Then
           call hess_mtrxC(bl(ider2 ),dens,hessian ,lind,thres1,
     *                     ntri,nblsiz,ngctoz,nintez,ncenter,
     *                     labels(lab1),labels(lab2),labels(lab3),
     *                     na,nsym)
          Else
           call hess_mtrxC_uhf(bl(ider2 ),dens,denB,hessian ,lind,
     $                         thres1,ntri,nblsiz,ngctoz,nintez,
     $                         ncenter,labels(lab1),labels(lab2),
     $                         labels(lab3),na,nsym)
          EndIf
        ENDIF
      ENDIF
c
        if(moreint) then
           go to 11
        endif
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
      end
c=====================================================================
