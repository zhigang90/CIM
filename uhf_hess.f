c==============================================================
      subroutine hess_mtrx_uhf(buf,    denA,  denB,  hess,  lind,
     $                         thres1, ntri,  nbls,  ngcd,  lnijkl,
     $                         ncenter,labels,length,lgenct,na,
     $                         nsym,   idft,  ax )
c----------------------------------------------------------------
c
c  this subroutine is derived from hess_mtrx
c
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
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
      dimension denA(ntri),denB(ntri)
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
                dij2=(denA(ijf)+denB(ijf))*2.d0
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                   dikA=denA(ikf)
                   djkA=denA(jkf)
                   dikB=denB(ikf)
                   djkB=denB(jkf)
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
                      dilA=denA(ilf)
                      djlA=denA(jlf)
                      dilB=denB(ilf)
                      djlB=denB(jlf)
                      dkl=denA(klf)+denB(klf)
c
                      If(idft.eq.0) Then
                      dijkl=2.0d0*(dij2*dkl - dikA*djlA - dilA*djkA
     $                                      - dikB*djlB - dilB*djkB)
                      Else
                      dijkl=2.0d0*(dij2*dkl - ax*
     $                  (dikA*djlA + dilA*djkA + dikB*djlB + dilB*djkB))
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
c=====================================================================
      subroutine hess_mtrxC_uhf(buf, denA, denB, hess, lind, thres1,
     *                          ntri,nbls, ngcd, lnijkl,
     *                          ncenter, labels, length, lgenct,
     *                          na, nsym )
c----------------------------------------------------------------
c
c  this subroutine is derived from hess_mtrxC
c
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
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
      dimension denA(ntri),denB(ntri)
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
                dij4=(denA(ijf)+denB(ijf))*4.0d0
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      klf=kk+lcf
                      if(lcf.gt.kcf) klf=ll+kcf
c
                      dkl=denA(klf)+denB(klf)
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
      subroutine fder_uhf(buf,denA,denB,fderA,fderB,lind,
     *                  ntri,nbls,ngcd,lnijkl,ncenter,
     *                  labels,length,lgenct,na,nsym )
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c STORAGE:
c labels    - space for the integral indices
c length(4) - the sizes of the 4 shells making up the shell integral.
c   E.g. if the 5 shells are of type L, P, D6 and S, then length=4,3,6,1
c lgenct    - it does not appear to be used. Somethig to do with general
c   contractions?
c na        - number of atoms
c nsym      - number of non-trivial symmetry operations (0 for C1 sym.)
c
c OUTPUT :
c fderA() - alpha Fock matrix derivatives
c fderB() - beta Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension fderA(3,na,ntri),fderB(3,na,ntri)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      data zero,half,one,two,four /0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from fock_deriv() :NBLS=',nbls,' ngcd=',ngcd
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
                dij2=(denA(ijf)+denB(ijf))*two
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                   dikA=denA(ikf)
                   djkA=denA(jkf)
                   dikB=denB(ikf)
                   djkB=denB(jkf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      ilf=ii+lcf
                      jlf=jj+lcf
                      klf=kk+lcf
                      if(lcf.gt.icf) ilf=ll+icf
                      if(lcf.gt.jcf) jlf=ll+jcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl2=(denA(klf)+denB(klf))*two
                      dilA=denA(ilf)
                      djlA=denA(jlf)
                      dilB=denB(ilf)
                      djlB=denB(jlf)
c
                      integ=integ+1
c------------------------------------------
c  the 12 derivatives
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
                      xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
                      yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
                      zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
c  Coulomb term 1  for atom 1
           fderA(1,iat,ijf)=fderA(1,iat,ijf)+xinta*dkl2
           fderA(2,iat,ijf)=fderA(2,iat,ijf)+yinta*dkl2
           fderA(3,iat,ijf)=fderA(3,iat,ijf)+zinta*dkl2
           fderB(1,iat,ijf)=fderB(1,iat,ijf)+xinta*dkl2
           fderB(2,iat,ijf)=fderB(2,iat,ijf)+yinta*dkl2
           fderB(3,iat,ijf)=fderB(3,iat,ijf)+zinta*dkl2
c  Coulomb term 1  for atom 2
           fderA(1,jat,ijf)=fderA(1,jat,ijf)+xintb*dkl2
           fderA(2,jat,ijf)=fderA(2,jat,ijf)+yintb*dkl2
           fderA(3,jat,ijf)=fderA(3,jat,ijf)+zintb*dkl2
           fderB(1,jat,ijf)=fderB(1,jat,ijf)+xintb*dkl2
           fderB(2,jat,ijf)=fderB(2,jat,ijf)+yintb*dkl2
           fderB(3,jat,ijf)=fderB(3,jat,ijf)+zintb*dkl2
c  Coulomb term 1  for atom 3
           fderA(1,kat,ijf)=fderA(1,kat,ijf)+xintc*dkl2
           fderA(2,kat,ijf)=fderA(2,kat,ijf)+yintc*dkl2
           fderA(3,kat,ijf)=fderA(3,kat,ijf)+zintc*dkl2
           fderB(1,kat,ijf)=fderB(1,kat,ijf)+xintc*dkl2
           fderB(2,kat,ijf)=fderB(2,kat,ijf)+yintc*dkl2
           fderB(3,kat,ijf)=fderB(3,kat,ijf)+zintc*dkl2
c  Coulomb term 1  for atom 4
           fderA(1,lat,ijf)=fderA(1,lat,ijf)+xintd*dkl2
           fderA(2,lat,ijf)=fderA(2,lat,ijf)+yintd*dkl2
           fderA(3,lat,ijf)=fderA(3,lat,ijf)+zintd*dkl2
           fderB(1,lat,ijf)=fderB(1,lat,ijf)+xintd*dkl2
           fderB(2,lat,ijf)=fderB(2,lat,ijf)+yintd*dkl2
           fderB(3,lat,ijf)=fderB(3,lat,ijf)+zintd*dkl2
c  Coulomb term 2  for atom 1
           fderA(1,iat,klf)=fderA(1,iat,klf)+xinta*dij2
           fderA(2,iat,klf)=fderA(2,iat,klf)+yinta*dij2
           fderA(3,iat,klf)=fderA(3,iat,klf)+zinta*dij2
           fderB(1,iat,klf)=fderB(1,iat,klf)+xinta*dij2
           fderB(2,iat,klf)=fderB(2,iat,klf)+yinta*dij2
           fderB(3,iat,klf)=fderB(3,iat,klf)+zinta*dij2
c  Coulomb term 2  for atom 2
           fderA(1,jat,klf)=fderA(1,jat,klf)+xintb*dij2
           fderA(2,jat,klf)=fderA(2,jat,klf)+yintb*dij2
           fderA(3,jat,klf)=fderA(3,jat,klf)+zintb*dij2
           fderB(1,jat,klf)=fderB(1,jat,klf)+xintb*dij2
           fderB(2,jat,klf)=fderB(2,jat,klf)+yintb*dij2
           fderB(3,jat,klf)=fderB(3,jat,klf)+zintb*dij2
c  Coulomb term 2  for atom 3
           fderA(1,kat,klf)=fderA(1,kat,klf)+xintc*dij2
           fderA(2,kat,klf)=fderA(2,kat,klf)+yintc*dij2
           fderA(3,kat,klf)=fderA(3,kat,klf)+zintc*dij2
           fderB(1,kat,klf)=fderB(1,kat,klf)+xintc*dij2
           fderB(2,kat,klf)=fderB(2,kat,klf)+yintc*dij2
           fderB(3,kat,klf)=fderB(3,kat,klf)+zintc*dij2
c  Coulomb term 2  for atom 4
           fderA(1,lat,klf)=fderA(1,lat,klf)+xintd*dij2
           fderA(2,lat,klf)=fderA(2,lat,klf)+yintd*dij2
           fderA(3,lat,klf)=fderA(3,lat,klf)+zintd*dij2
           fderB(1,lat,klf)=fderB(1,lat,klf)+xintd*dij2
           fderB(2,lat,klf)=fderB(2,lat,klf)+yintd*dij2
           fderB(3,lat,klf)=fderB(3,lat,klf)+zintd*dij2
c  Exchange term 1  for atom 1
           fderA(1,iat,ikf)=fderA(1,iat,ikf)-xinta*djlA
           fderA(2,iat,ikf)=fderA(2,iat,ikf)-yinta*djlA
           fderA(3,iat,ikf)=fderA(3,iat,ikf)-zinta*djlA
           fderB(1,iat,ikf)=fderB(1,iat,ikf)-xinta*djlB
           fderB(2,iat,ikf)=fderB(2,iat,ikf)-yinta*djlB
           fderB(3,iat,ikf)=fderB(3,iat,ikf)-zinta*djlB
c  Exchange term 1  for atom 2
           fderA(1,jat,ikf)=fderA(1,jat,ikf)-xintb*djlA
           fderA(2,jat,ikf)=fderA(2,jat,ikf)-yintb*djlA
           fderA(3,jat,ikf)=fderA(3,jat,ikf)-zintb*djlA
           fderB(1,jat,ikf)=fderB(1,jat,ikf)-xintb*djlB
           fderB(2,jat,ikf)=fderB(2,jat,ikf)-yintb*djlB
           fderB(3,jat,ikf)=fderB(3,jat,ikf)-zintb*djlB
c  Exchange term 1  for atom 3
           fderA(1,kat,ikf)=fderA(1,kat,ikf)-xintc*djlA
           fderA(2,kat,ikf)=fderA(2,kat,ikf)-yintc*djlA
           fderA(3,kat,ikf)=fderA(3,kat,ikf)-zintc*djlA
           fderB(1,kat,ikf)=fderB(1,kat,ikf)-xintc*djlB
           fderB(2,kat,ikf)=fderB(2,kat,ikf)-yintc*djlB
           fderB(3,kat,ikf)=fderB(3,kat,ikf)-zintc*djlB
c  Exchange term 1  for atom 4
           fderA(1,lat,ikf)=fderA(1,lat,ikf)-xintd*djlA
           fderA(2,lat,ikf)=fderA(2,lat,ikf)-yintd*djlA
           fderA(3,lat,ikf)=fderA(3,lat,ikf)-zintd*djlA
           fderB(1,lat,ikf)=fderB(1,lat,ikf)-xintd*djlB
           fderB(2,lat,ikf)=fderB(2,lat,ikf)-yintd*djlB
           fderB(3,lat,ikf)=fderB(3,lat,ikf)-zintd*djlB
c  Exchange term 2  for atom 1
           fderA(1,iat,jlf)=fderA(1,iat,jlf)-xinta*dikA
           fderA(2,iat,jlf)=fderA(2,iat,jlf)-yinta*dikA
           fderA(3,iat,jlf)=fderA(3,iat,jlf)-zinta*dikA
           fderB(1,iat,jlf)=fderB(1,iat,jlf)-xinta*dikB
           fderB(2,iat,jlf)=fderB(2,iat,jlf)-yinta*dikB
           fderB(3,iat,jlf)=fderB(3,iat,jlf)-zinta*dikB
c  Exchange term 2  for atom 2
           fderA(1,jat,jlf)=fderA(1,jat,jlf)-xintb*dikA
           fderA(2,jat,jlf)=fderA(2,jat,jlf)-yintb*dikA
           fderA(3,jat,jlf)=fderA(3,jat,jlf)-zintb*dikA
           fderB(1,jat,jlf)=fderB(1,jat,jlf)-xintb*dikB
           fderB(2,jat,jlf)=fderB(2,jat,jlf)-yintb*dikB
           fderB(3,jat,jlf)=fderB(3,jat,jlf)-zintb*dikB
c  Exchange term 2  for atom 3
           fderA(1,kat,jlf)=fderA(1,kat,jlf)-xintc*dikA
           fderA(2,kat,jlf)=fderA(2,kat,jlf)-yintc*dikA
           fderA(3,kat,jlf)=fderA(3,kat,jlf)-zintc*dikA
           fderB(1,kat,jlf)=fderB(1,kat,jlf)-xintc*dikB
           fderB(2,kat,jlf)=fderB(2,kat,jlf)-yintc*dikB
           fderB(3,kat,jlf)=fderB(3,kat,jlf)-zintc*dikB
c  Exchange term 2  for atom 4
           fderA(1,lat,jlf)=fderA(1,lat,jlf)-xintd*dikA
           fderA(2,lat,jlf)=fderA(2,lat,jlf)-yintd*dikA
           fderA(3,lat,jlf)=fderA(3,lat,jlf)-zintd*dikA
           fderB(1,lat,jlf)=fderB(1,lat,jlf)-xintd*dikB
           fderB(2,lat,jlf)=fderB(2,lat,jlf)-yintd*dikB
           fderB(3,lat,jlf)=fderB(3,lat,jlf)-zintd*dikB
c  Exchange term 3  for atom 1
           fderA(1,iat,jkf)=fderA(1,iat,jkf)-xinta*dilA
           fderA(2,iat,jkf)=fderA(2,iat,jkf)-yinta*dilA
           fderA(3,iat,jkf)=fderA(3,iat,jkf)-zinta*dilA
           fderB(1,iat,jkf)=fderB(1,iat,jkf)-xinta*dilB
           fderB(2,iat,jkf)=fderB(2,iat,jkf)-yinta*dilB
           fderB(3,iat,jkf)=fderB(3,iat,jkf)-zinta*dilB
c  Exchange term 3  for atom 2
           fderA(1,jat,jkf)=fderA(1,jat,jkf)-xintb*dilA
           fderA(2,jat,jkf)=fderA(2,jat,jkf)-yintb*dilA
           fderA(3,jat,jkf)=fderA(3,jat,jkf)-zintb*dilA
           fderB(1,jat,jkf)=fderB(1,jat,jkf)-xintb*dilB
           fderB(2,jat,jkf)=fderB(2,jat,jkf)-yintb*dilB
           fderB(3,jat,jkf)=fderB(3,jat,jkf)-zintb*dilB
c  Exchange term 3  for atom 3
           fderA(1,kat,jkf)=fderA(1,kat,jkf)-xintc*dilA
           fderA(2,kat,jkf)=fderA(2,kat,jkf)-yintc*dilA
           fderA(3,kat,jkf)=fderA(3,kat,jkf)-zintc*dilA
           fderB(1,kat,jkf)=fderB(1,kat,jkf)-xintc*dilB
           fderB(2,kat,jkf)=fderB(2,kat,jkf)-yintc*dilB
           fderB(3,kat,jkf)=fderB(3,kat,jkf)-zintc*dilB
c  Exchange term 3  for atom 4
           fderA(1,lat,jkf)=fderA(1,lat,jkf)-xintd*dilA
           fderA(2,lat,jkf)=fderA(2,lat,jkf)-yintd*dilA
           fderA(3,lat,jkf)=fderA(3,lat,jkf)-zintd*dilA
           fderB(1,lat,jkf)=fderB(1,lat,jkf)-xintd*dilB
           fderB(2,lat,jkf)=fderB(2,lat,jkf)-yintd*dilB
           fderB(3,lat,jkf)=fderB(3,lat,jkf)-zintd*dilB
c  Exchange term 4  for atom 1
           fderA(1,iat,ilf)=fderA(1,iat,ilf)-xinta*djkA
           fderA(2,iat,ilf)=fderA(2,iat,ilf)-yinta*djkA
           fderA(3,iat,ilf)=fderA(3,iat,ilf)-zinta*djkA
           fderB(1,iat,ilf)=fderB(1,iat,ilf)-xinta*djkB
           fderB(2,iat,ilf)=fderB(2,iat,ilf)-yinta*djkB
           fderB(3,iat,ilf)=fderB(3,iat,ilf)-zinta*djkB
c  Exchange term 4  for atom 2
           fderA(1,jat,ilf)=fderA(1,jat,ilf)-xintb*djkA
           fderA(2,jat,ilf)=fderA(2,jat,ilf)-yintb*djkA
           fderA(3,jat,ilf)=fderA(3,jat,ilf)-zintb*djkA
           fderB(1,jat,ilf)=fderB(1,jat,ilf)-xintb*djkB
           fderB(2,jat,ilf)=fderB(2,jat,ilf)-yintb*djkB
           fderB(3,jat,ilf)=fderB(3,jat,ilf)-zintb*djkB
c  Exchange term 4  for atom 3
           fderA(1,kat,ilf)=fderA(1,kat,ilf)-xintc*djkA
           fderA(2,kat,ilf)=fderA(2,kat,ilf)-yintc*djkA
           fderA(3,kat,ilf)=fderA(3,kat,ilf)-zintc*djkA
           fderB(1,kat,ilf)=fderB(1,kat,ilf)-xintc*djkB
           fderB(2,kat,ilf)=fderB(2,kat,ilf)-yintc*djkB
           fderB(3,kat,ilf)=fderB(3,kat,ilf)-zintc*djkB
c  Exchange term 4  for atom 4
           fderA(1,lat,ilf)=fderA(1,lat,ilf)-xintd*djkA
           fderA(2,lat,ilf)=fderA(2,lat,ilf)-yintd*djkA
           fderA(3,lat,ilf)=fderA(3,lat,ilf)-zintd*djkA
           fderB(1,lat,ilf)=fderB(1,lat,ilf)-xintd*djkB
           fderB(2,lat,ilf)=fderB(2,lat,ilf)-yintd*djkB
           fderB(3,lat,ilf)=fderB(3,lat,ilf)-zintd*djkB
c
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine fder_uacm(buf,denA,denB,fderA,fderB,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,nsym,ax )
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c STORAGE:
c labels    - space for the integral indices
c length(4) - the sizes of the 4 shells making up the shell integral.
c   E.g. if the 5 shells are of type L, P, D6 and S, then length=4,3,6,1
c lgenct    - it does not appear to be used. Somethig to do with general
c   contractions?
c na        - number of atoms
c nsym      - number of non-trivial symmetry operations (0 for C1 sym.)
c  ax       - factor to multiply exchange contribution (for hybrid DFT)
c
c OUTPUT :
c fderA() - alpha Fock matrix derivatives
c fderB() - beta Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension fderA(3,na,ntri),fderB(3,na,ntri)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      data zero,half,one,two,four /0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from fock_deriv() :NBLS=',nbls,' ngcd=',ngcd
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
                dij2=(denA(ijf)+denB(ijf))*two
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                   dikA=denA(ikf)*ax
                   djkA=denA(jkf)*ax
                   dikB=denB(ikf)*ax
                   djkB=denB(jkf)*ax
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      ilf=ii+lcf
                      jlf=jj+lcf
                      klf=kk+lcf
                      if(lcf.gt.icf) ilf=ll+icf
                      if(lcf.gt.jcf) jlf=ll+jcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl2=(denA(klf)+denB(klf))*two
                      dilA=denA(ilf)*ax
                      djlA=denA(jlf)*ax
                      dilB=denB(ilf)*ax
                      djlB=denB(jlf)*ax
c
                      integ=integ+1
c------------------------------------------
c  the 12 derivatives
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
                      xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
                      yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
                      zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
c  Coulomb term 1  for atom 1
           fderA(1,iat,ijf)=fderA(1,iat,ijf)+xinta*dkl2
           fderA(2,iat,ijf)=fderA(2,iat,ijf)+yinta*dkl2
           fderA(3,iat,ijf)=fderA(3,iat,ijf)+zinta*dkl2
           fderB(1,iat,ijf)=fderB(1,iat,ijf)+xinta*dkl2
           fderB(2,iat,ijf)=fderB(2,iat,ijf)+yinta*dkl2
           fderB(3,iat,ijf)=fderB(3,iat,ijf)+zinta*dkl2
c  Coulomb term 1  for atom 2
           fderA(1,jat,ijf)=fderA(1,jat,ijf)+xintb*dkl2
           fderA(2,jat,ijf)=fderA(2,jat,ijf)+yintb*dkl2
           fderA(3,jat,ijf)=fderA(3,jat,ijf)+zintb*dkl2
           fderB(1,jat,ijf)=fderB(1,jat,ijf)+xintb*dkl2
           fderB(2,jat,ijf)=fderB(2,jat,ijf)+yintb*dkl2
           fderB(3,jat,ijf)=fderB(3,jat,ijf)+zintb*dkl2
c  Coulomb term 1  for atom 3
           fderA(1,kat,ijf)=fderA(1,kat,ijf)+xintc*dkl2
           fderA(2,kat,ijf)=fderA(2,kat,ijf)+yintc*dkl2
           fderA(3,kat,ijf)=fderA(3,kat,ijf)+zintc*dkl2
           fderB(1,kat,ijf)=fderB(1,kat,ijf)+xintc*dkl2
           fderB(2,kat,ijf)=fderB(2,kat,ijf)+yintc*dkl2
           fderB(3,kat,ijf)=fderB(3,kat,ijf)+zintc*dkl2
c  Coulomb term 1  for atom 4
           fderA(1,lat,ijf)=fderA(1,lat,ijf)+xintd*dkl2
           fderA(2,lat,ijf)=fderA(2,lat,ijf)+yintd*dkl2
           fderA(3,lat,ijf)=fderA(3,lat,ijf)+zintd*dkl2
           fderB(1,lat,ijf)=fderB(1,lat,ijf)+xintd*dkl2
           fderB(2,lat,ijf)=fderB(2,lat,ijf)+yintd*dkl2
           fderB(3,lat,ijf)=fderB(3,lat,ijf)+zintd*dkl2
c  Coulomb term 2  for atom 1
           fderA(1,iat,klf)=fderA(1,iat,klf)+xinta*dij2
           fderA(2,iat,klf)=fderA(2,iat,klf)+yinta*dij2
           fderA(3,iat,klf)=fderA(3,iat,klf)+zinta*dij2
           fderB(1,iat,klf)=fderB(1,iat,klf)+xinta*dij2
           fderB(2,iat,klf)=fderB(2,iat,klf)+yinta*dij2
           fderB(3,iat,klf)=fderB(3,iat,klf)+zinta*dij2
c  Coulomb term 2  for atom 2
           fderA(1,jat,klf)=fderA(1,jat,klf)+xintb*dij2
           fderA(2,jat,klf)=fderA(2,jat,klf)+yintb*dij2
           fderA(3,jat,klf)=fderA(3,jat,klf)+zintb*dij2
           fderB(1,jat,klf)=fderB(1,jat,klf)+xintb*dij2
           fderB(2,jat,klf)=fderB(2,jat,klf)+yintb*dij2
           fderB(3,jat,klf)=fderB(3,jat,klf)+zintb*dij2
c  Coulomb term 2  for atom 3
           fderA(1,kat,klf)=fderA(1,kat,klf)+xintc*dij2
           fderA(2,kat,klf)=fderA(2,kat,klf)+yintc*dij2
           fderA(3,kat,klf)=fderA(3,kat,klf)+zintc*dij2
           fderB(1,kat,klf)=fderB(1,kat,klf)+xintc*dij2
           fderB(2,kat,klf)=fderB(2,kat,klf)+yintc*dij2
           fderB(3,kat,klf)=fderB(3,kat,klf)+zintc*dij2
c  Coulomb term 2  for atom 4
           fderA(1,lat,klf)=fderA(1,lat,klf)+xintd*dij2
           fderA(2,lat,klf)=fderA(2,lat,klf)+yintd*dij2
           fderA(3,lat,klf)=fderA(3,lat,klf)+zintd*dij2
           fderB(1,lat,klf)=fderB(1,lat,klf)+xintd*dij2
           fderB(2,lat,klf)=fderB(2,lat,klf)+yintd*dij2
           fderB(3,lat,klf)=fderB(3,lat,klf)+zintd*dij2
c  Exchange term 1  for atom 1
           fderA(1,iat,ikf)=fderA(1,iat,ikf)-xinta*djlA
           fderA(2,iat,ikf)=fderA(2,iat,ikf)-yinta*djlA
           fderA(3,iat,ikf)=fderA(3,iat,ikf)-zinta*djlA
           fderB(1,iat,ikf)=fderB(1,iat,ikf)-xinta*djlB
           fderB(2,iat,ikf)=fderB(2,iat,ikf)-yinta*djlB
           fderB(3,iat,ikf)=fderB(3,iat,ikf)-zinta*djlB
c  Exchange term 1  for atom 2
           fderA(1,jat,ikf)=fderA(1,jat,ikf)-xintb*djlA
           fderA(2,jat,ikf)=fderA(2,jat,ikf)-yintb*djlA
           fderA(3,jat,ikf)=fderA(3,jat,ikf)-zintb*djlA
           fderB(1,jat,ikf)=fderB(1,jat,ikf)-xintb*djlB
           fderB(2,jat,ikf)=fderB(2,jat,ikf)-yintb*djlB
           fderB(3,jat,ikf)=fderB(3,jat,ikf)-zintb*djlB
c  Exchange term 1  for atom 3
           fderA(1,kat,ikf)=fderA(1,kat,ikf)-xintc*djlA
           fderA(2,kat,ikf)=fderA(2,kat,ikf)-yintc*djlA
           fderA(3,kat,ikf)=fderA(3,kat,ikf)-zintc*djlA
           fderB(1,kat,ikf)=fderB(1,kat,ikf)-xintc*djlB
           fderB(2,kat,ikf)=fderB(2,kat,ikf)-yintc*djlB
           fderB(3,kat,ikf)=fderB(3,kat,ikf)-zintc*djlB
c  Exchange term 1  for atom 4
           fderA(1,lat,ikf)=fderA(1,lat,ikf)-xintd*djlA
           fderA(2,lat,ikf)=fderA(2,lat,ikf)-yintd*djlA
           fderA(3,lat,ikf)=fderA(3,lat,ikf)-zintd*djlA
           fderB(1,lat,ikf)=fderB(1,lat,ikf)-xintd*djlB
           fderB(2,lat,ikf)=fderB(2,lat,ikf)-yintd*djlB
           fderB(3,lat,ikf)=fderB(3,lat,ikf)-zintd*djlB
c  Exchange term 2  for atom 1
           fderA(1,iat,jlf)=fderA(1,iat,jlf)-xinta*dikA
           fderA(2,iat,jlf)=fderA(2,iat,jlf)-yinta*dikA
           fderA(3,iat,jlf)=fderA(3,iat,jlf)-zinta*dikA
           fderB(1,iat,jlf)=fderB(1,iat,jlf)-xinta*dikB
           fderB(2,iat,jlf)=fderB(2,iat,jlf)-yinta*dikB
           fderB(3,iat,jlf)=fderB(3,iat,jlf)-zinta*dikB
c  Exchange term 2  for atom 2
           fderA(1,jat,jlf)=fderA(1,jat,jlf)-xintb*dikA
           fderA(2,jat,jlf)=fderA(2,jat,jlf)-yintb*dikA
           fderA(3,jat,jlf)=fderA(3,jat,jlf)-zintb*dikA
           fderB(1,jat,jlf)=fderB(1,jat,jlf)-xintb*dikB
           fderB(2,jat,jlf)=fderB(2,jat,jlf)-yintb*dikB
           fderB(3,jat,jlf)=fderB(3,jat,jlf)-zintb*dikB
c  Exchange term 2  for atom 3
           fderA(1,kat,jlf)=fderA(1,kat,jlf)-xintc*dikA
           fderA(2,kat,jlf)=fderA(2,kat,jlf)-yintc*dikA
           fderA(3,kat,jlf)=fderA(3,kat,jlf)-zintc*dikA
           fderB(1,kat,jlf)=fderB(1,kat,jlf)-xintc*dikB
           fderB(2,kat,jlf)=fderB(2,kat,jlf)-yintc*dikB
           fderB(3,kat,jlf)=fderB(3,kat,jlf)-zintc*dikB
c  Exchange term 2  for atom 4
           fderA(1,lat,jlf)=fderA(1,lat,jlf)-xintd*dikA
           fderA(2,lat,jlf)=fderA(2,lat,jlf)-yintd*dikA
           fderA(3,lat,jlf)=fderA(3,lat,jlf)-zintd*dikA
           fderB(1,lat,jlf)=fderB(1,lat,jlf)-xintd*dikB
           fderB(2,lat,jlf)=fderB(2,lat,jlf)-yintd*dikB
           fderB(3,lat,jlf)=fderB(3,lat,jlf)-zintd*dikB
c  Exchange term 3  for atom 1
           fderA(1,iat,jkf)=fderA(1,iat,jkf)-xinta*dilA
           fderA(2,iat,jkf)=fderA(2,iat,jkf)-yinta*dilA
           fderA(3,iat,jkf)=fderA(3,iat,jkf)-zinta*dilA
           fderB(1,iat,jkf)=fderB(1,iat,jkf)-xinta*dilB
           fderB(2,iat,jkf)=fderB(2,iat,jkf)-yinta*dilB
           fderB(3,iat,jkf)=fderB(3,iat,jkf)-zinta*dilB
c  Exchange term 3  for atom 2
           fderA(1,jat,jkf)=fderA(1,jat,jkf)-xintb*dilA
           fderA(2,jat,jkf)=fderA(2,jat,jkf)-yintb*dilA
           fderA(3,jat,jkf)=fderA(3,jat,jkf)-zintb*dilA
           fderB(1,jat,jkf)=fderB(1,jat,jkf)-xintb*dilB
           fderB(2,jat,jkf)=fderB(2,jat,jkf)-yintb*dilB
           fderB(3,jat,jkf)=fderB(3,jat,jkf)-zintb*dilB
c  Exchange term 3  for atom 3
           fderA(1,kat,jkf)=fderA(1,kat,jkf)-xintc*dilA
           fderA(2,kat,jkf)=fderA(2,kat,jkf)-yintc*dilA
           fderA(3,kat,jkf)=fderA(3,kat,jkf)-zintc*dilA
           fderB(1,kat,jkf)=fderB(1,kat,jkf)-xintc*dilB
           fderB(2,kat,jkf)=fderB(2,kat,jkf)-yintc*dilB
           fderB(3,kat,jkf)=fderB(3,kat,jkf)-zintc*dilB
c  Exchange term 3  for atom 4
           fderA(1,lat,jkf)=fderA(1,lat,jkf)-xintd*dilA
           fderA(2,lat,jkf)=fderA(2,lat,jkf)-yintd*dilA
           fderA(3,lat,jkf)=fderA(3,lat,jkf)-zintd*dilA
           fderB(1,lat,jkf)=fderB(1,lat,jkf)-xintd*dilB
           fderB(2,lat,jkf)=fderB(2,lat,jkf)-yintd*dilB
           fderB(3,lat,jkf)=fderB(3,lat,jkf)-zintd*dilB
c  Exchange term 4  for atom 1
           fderA(1,iat,ilf)=fderA(1,iat,ilf)-xinta*djkA
           fderA(2,iat,ilf)=fderA(2,iat,ilf)-yinta*djkA
           fderA(3,iat,ilf)=fderA(3,iat,ilf)-zinta*djkA
           fderB(1,iat,ilf)=fderB(1,iat,ilf)-xinta*djkB
           fderB(2,iat,ilf)=fderB(2,iat,ilf)-yinta*djkB
           fderB(3,iat,ilf)=fderB(3,iat,ilf)-zinta*djkB
c  Exchange term 4  for atom 2
           fderA(1,jat,ilf)=fderA(1,jat,ilf)-xintb*djkA
           fderA(2,jat,ilf)=fderA(2,jat,ilf)-yintb*djkA
           fderA(3,jat,ilf)=fderA(3,jat,ilf)-zintb*djkA
           fderB(1,jat,ilf)=fderB(1,jat,ilf)-xintb*djkB
           fderB(2,jat,ilf)=fderB(2,jat,ilf)-yintb*djkB
           fderB(3,jat,ilf)=fderB(3,jat,ilf)-zintb*djkB
c  Exchange term 4  for atom 3
           fderA(1,kat,ilf)=fderA(1,kat,ilf)-xintc*djkA
           fderA(2,kat,ilf)=fderA(2,kat,ilf)-yintc*djkA
           fderA(3,kat,ilf)=fderA(3,kat,ilf)-zintc*djkA
           fderB(1,kat,ilf)=fderB(1,kat,ilf)-xintc*djkB
           fderB(2,kat,ilf)=fderB(2,kat,ilf)-yintc*djkB
           fderB(3,kat,ilf)=fderB(3,kat,ilf)-zintc*djkB
c  Exchange term 4  for atom 4
           fderA(1,lat,ilf)=fderA(1,lat,ilf)-xintd*djkA
           fderA(2,lat,ilf)=fderA(2,lat,ilf)-yintd*djkA
           fderA(3,lat,ilf)=fderA(3,lat,ilf)-zintd*djkA
           fderB(1,lat,ilf)=fderB(1,lat,ilf)-xintd*djkB
           fderB(2,lat,ilf)=fderB(2,lat,ilf)-yintd*djkB
           fderB(3,lat,ilf)=fderB(3,lat,ilf)-zintd*djkB
c
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine fder_uhfC(buf,denA,denB,fderA,fderB,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,nsym )
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
c----------------------------------------------------------------
c
c  ** COULOMB ONLY **
c
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c STORAGE:
c labels    - space for the integral indices
c length(4) - the sizes of the 4 shells making up the shell integral.
c   E.g. if the 5 shells are of type L, P, D6 and S, then length=4,3,6,1
c lgenct    - it does not appear to be used. Somethig to do with general
c   contractions?
c na        - number of atoms
c nsym      - number of non-trivial symmetry operations (0 for C1 sym.)
c
c OUTPUT :
c fderA() - alpha Fock matrix derivatives
c fderB() - beta Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension fderA(3,na,ntri),fderB(3,na,ntri)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      data zero,half,one,two,four /0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from fock_deriv() :NBLS=',nbls,' ngcd=',ngcd
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
                dij2=(denA(ijf)+denB(ijf))*two
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      klf=kk+lcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl2=(denA(klf)+denB(klf))*two
c
                      integ=integ+1
c------------------------------------------
c  the 12 derivatives
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
                      xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
                      yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
                      zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
c  Coulomb term 1  for atom 1
           fderA(1,iat,ijf)=fderA(1,iat,ijf)+xinta*dkl2
           fderA(2,iat,ijf)=fderA(2,iat,ijf)+yinta*dkl2
           fderA(3,iat,ijf)=fderA(3,iat,ijf)+zinta*dkl2
           fderB(1,iat,ijf)=fderB(1,iat,ijf)+xinta*dkl2
           fderB(2,iat,ijf)=fderB(2,iat,ijf)+yinta*dkl2
           fderB(3,iat,ijf)=fderB(3,iat,ijf)+zinta*dkl2
c  Coulomb term 1  for atom 2
           fderA(1,jat,ijf)=fderA(1,jat,ijf)+xintb*dkl2
           fderA(2,jat,ijf)=fderA(2,jat,ijf)+yintb*dkl2
           fderA(3,jat,ijf)=fderA(3,jat,ijf)+zintb*dkl2
           fderB(1,jat,ijf)=fderB(1,jat,ijf)+xintb*dkl2
           fderB(2,jat,ijf)=fderB(2,jat,ijf)+yintb*dkl2
           fderB(3,jat,ijf)=fderB(3,jat,ijf)+zintb*dkl2
c  Coulomb term 1  for atom 3
           fderA(1,kat,ijf)=fderA(1,kat,ijf)+xintc*dkl2
           fderA(2,kat,ijf)=fderA(2,kat,ijf)+yintc*dkl2
           fderA(3,kat,ijf)=fderA(3,kat,ijf)+zintc*dkl2
           fderB(1,kat,ijf)=fderB(1,kat,ijf)+xintc*dkl2
           fderB(2,kat,ijf)=fderB(2,kat,ijf)+yintc*dkl2
           fderB(3,kat,ijf)=fderB(3,kat,ijf)+zintc*dkl2
c  Coulomb term 1  for atom 4
           fderA(1,lat,ijf)=fderA(1,lat,ijf)+xintd*dkl2
           fderA(2,lat,ijf)=fderA(2,lat,ijf)+yintd*dkl2
           fderA(3,lat,ijf)=fderA(3,lat,ijf)+zintd*dkl2
           fderB(1,lat,ijf)=fderB(1,lat,ijf)+xintd*dkl2
           fderB(2,lat,ijf)=fderB(2,lat,ijf)+yintd*dkl2
           fderB(3,lat,ijf)=fderB(3,lat,ijf)+zintd*dkl2
c  Coulomb term 2  for atom 1
           fderA(1,iat,klf)=fderA(1,iat,klf)+xinta*dij2
           fderA(2,iat,klf)=fderA(2,iat,klf)+yinta*dij2
           fderA(3,iat,klf)=fderA(3,iat,klf)+zinta*dij2
           fderB(1,iat,klf)=fderB(1,iat,klf)+xinta*dij2
           fderB(2,iat,klf)=fderB(2,iat,klf)+yinta*dij2
           fderB(3,iat,klf)=fderB(3,iat,klf)+zinta*dij2
c  Coulomb term 2  for atom 2
           fderA(1,jat,klf)=fderA(1,jat,klf)+xintb*dij2
           fderA(2,jat,klf)=fderA(2,jat,klf)+yintb*dij2
           fderA(3,jat,klf)=fderA(3,jat,klf)+zintb*dij2
           fderB(1,jat,klf)=fderB(1,jat,klf)+xintb*dij2
           fderB(2,jat,klf)=fderB(2,jat,klf)+yintb*dij2
           fderB(3,jat,klf)=fderB(3,jat,klf)+zintb*dij2
c  Coulomb term 2  for atom 3
           fderA(1,kat,klf)=fderA(1,kat,klf)+xintc*dij2
           fderA(2,kat,klf)=fderA(2,kat,klf)+yintc*dij2
           fderA(3,kat,klf)=fderA(3,kat,klf)+zintc*dij2
           fderB(1,kat,klf)=fderB(1,kat,klf)+xintc*dij2
           fderB(2,kat,klf)=fderB(2,kat,klf)+yintc*dij2
           fderB(3,kat,klf)=fderB(3,kat,klf)+zintc*dij2
c  Coulomb term 2  for atom 4
           fderA(1,lat,klf)=fderA(1,lat,klf)+xintd*dij2
           fderA(2,lat,klf)=fderA(2,lat,klf)+yintd*dij2
           fderA(3,lat,klf)=fderA(3,lat,klf)+zintd*dij2
           fderB(1,lat,klf)=fderB(1,lat,klf)+xintd*dij2
           fderB(2,lat,klf)=fderB(2,lat,klf)+yintd*dij2
           fderB(3,lat,klf)=fderB(3,lat,klf)+zintd*dij2
c
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine fderN_uacm(buf,denA,denB,fderA,fderB,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,nsym,ax,natb,nate)
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c STORAGE:
c labels    - space for the integral indices
c length(4) - the sizes of the 4 shells making up the shell integral.
c   E.g. if the 5 shells are of type L, P, D6 and S, then length=4,3,6,1
c lgenct    - it does not appear to be used. Somethig to do with general
c   contractions?
c na        - number of atoms
c nsym      - number of non-trivial symmetry operations (0 for C1 sym.)
c  ax       - factor to multiply exchange contribution (for hybrid DFT)
c  natb     - fisrt componet to be computed in current pass
c  nate     - last componet to be computed in current pass
c
c OUTPUT :
c fderA() - alpha Fock matrix derivatives
c fderB() - beta Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat,dojat,dokat,dolat
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension fderA(3,natb:nate,ntri),fderB(3,natb:nate,ntri)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      data zero,half,one,two,four /0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from fock_deriv() :NBLS=',nbls,' ngcd=',ngcd
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
c-------------------------------------------------------
          doiat=.false.
          if(iat.ge.natb .and. iat.le.nate) doiat=.true.
          dojat=.false.
          if(jat.ge.natb .and. jat.le.nate) dojat=.true.
          dokat=.false.
          if(kat.ge.natb .and. kat.le.nate) dokat=.true.
          dolat=.false.
          if(lat.ge.natb .and. lat.le.nate) dolat=.true.
c-------------------------------------------------------
      if(.not.doiat .and. .not.dojat .and.
     *   .not.dokat .and. .not.dolat       ) go to 150
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
                dij2=(denA(ijf)+denB(ijf))*two
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                   dikA=denA(ikf)*ax
                   djkA=denA(jkf)*ax
                   dikB=denB(ikf)*ax
                   djkB=denB(jkf)*ax
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      ilf=ii+lcf
                      jlf=jj+lcf
                      klf=kk+lcf
                      if(lcf.gt.icf) ilf=ll+icf
                      if(lcf.gt.jcf) jlf=ll+jcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl2=(denA(klf)+denB(klf))*two
                      dilA=denA(ilf)*ax
                      djlA=denA(jlf)*ax
                      dilB=denB(ilf)*ax
                      djlB=denB(jlf)*ax
c
                      integ=integ+1
c------------------------------------------
c  the 12 derivatives
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
                      xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
                      yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
                      zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
         if(doiat) then
c  Coulomb term 1  for atom 1
           fderA(1,iat,ijf)=fderA(1,iat,ijf)+xinta*dkl2
           fderA(2,iat,ijf)=fderA(2,iat,ijf)+yinta*dkl2
           fderA(3,iat,ijf)=fderA(3,iat,ijf)+zinta*dkl2
           fderB(1,iat,ijf)=fderB(1,iat,ijf)+xinta*dkl2
           fderB(2,iat,ijf)=fderB(2,iat,ijf)+yinta*dkl2
           fderB(3,iat,ijf)=fderB(3,iat,ijf)+zinta*dkl2
c  Coulomb term 2  for atom 1
           fderA(1,iat,klf)=fderA(1,iat,klf)+xinta*dij2
           fderA(2,iat,klf)=fderA(2,iat,klf)+yinta*dij2
           fderA(3,iat,klf)=fderA(3,iat,klf)+zinta*dij2
           fderB(1,iat,klf)=fderB(1,iat,klf)+xinta*dij2
           fderB(2,iat,klf)=fderB(2,iat,klf)+yinta*dij2
           fderB(3,iat,klf)=fderB(3,iat,klf)+zinta*dij2
c  Exchange term 1  for atom 1
           fderA(1,iat,ikf)=fderA(1,iat,ikf)-xinta*djlA
           fderA(2,iat,ikf)=fderA(2,iat,ikf)-yinta*djlA
           fderA(3,iat,ikf)=fderA(3,iat,ikf)-zinta*djlA
           fderB(1,iat,ikf)=fderB(1,iat,ikf)-xinta*djlB
           fderB(2,iat,ikf)=fderB(2,iat,ikf)-yinta*djlB
           fderB(3,iat,ikf)=fderB(3,iat,ikf)-zinta*djlB
c  Exchange term 2  for atom 1
           fderA(1,iat,jlf)=fderA(1,iat,jlf)-xinta*dikA
           fderA(2,iat,jlf)=fderA(2,iat,jlf)-yinta*dikA
           fderA(3,iat,jlf)=fderA(3,iat,jlf)-zinta*dikA
           fderB(1,iat,jlf)=fderB(1,iat,jlf)-xinta*dikB
           fderB(2,iat,jlf)=fderB(2,iat,jlf)-yinta*dikB
           fderB(3,iat,jlf)=fderB(3,iat,jlf)-zinta*dikB
c  Exchange term 3  for atom 1
           fderA(1,iat,jkf)=fderA(1,iat,jkf)-xinta*dilA
           fderA(2,iat,jkf)=fderA(2,iat,jkf)-yinta*dilA
           fderA(3,iat,jkf)=fderA(3,iat,jkf)-zinta*dilA
           fderB(1,iat,jkf)=fderB(1,iat,jkf)-xinta*dilB
           fderB(2,iat,jkf)=fderB(2,iat,jkf)-yinta*dilB
           fderB(3,iat,jkf)=fderB(3,iat,jkf)-zinta*dilB
c  Exchange term 4  for atom 1
           fderA(1,iat,ilf)=fderA(1,iat,ilf)-xinta*djkA
           fderA(2,iat,ilf)=fderA(2,iat,ilf)-yinta*djkA
           fderA(3,iat,ilf)=fderA(3,iat,ilf)-zinta*djkA
           fderB(1,iat,ilf)=fderB(1,iat,ilf)-xinta*djkB
           fderB(2,iat,ilf)=fderB(2,iat,ilf)-yinta*djkB
           fderB(3,iat,ilf)=fderB(3,iat,ilf)-zinta*djkB
         endif
         if(dojat) then
c  Coulomb term 1  for atom 2
           fderA(1,jat,ijf)=fderA(1,jat,ijf)+xintb*dkl2
           fderA(2,jat,ijf)=fderA(2,jat,ijf)+yintb*dkl2
           fderA(3,jat,ijf)=fderA(3,jat,ijf)+zintb*dkl2
           fderB(1,jat,ijf)=fderB(1,jat,ijf)+xintb*dkl2
           fderB(2,jat,ijf)=fderB(2,jat,ijf)+yintb*dkl2
           fderB(3,jat,ijf)=fderB(3,jat,ijf)+zintb*dkl2
c  Coulomb term 2  for atom 2
           fderA(1,jat,klf)=fderA(1,jat,klf)+xintb*dij2
           fderA(2,jat,klf)=fderA(2,jat,klf)+yintb*dij2
           fderA(3,jat,klf)=fderA(3,jat,klf)+zintb*dij2
           fderB(1,jat,klf)=fderB(1,jat,klf)+xintb*dij2
           fderB(2,jat,klf)=fderB(2,jat,klf)+yintb*dij2
           fderB(3,jat,klf)=fderB(3,jat,klf)+zintb*dij2
c  Exchange term 1  for atom 2
           fderA(1,jat,ikf)=fderA(1,jat,ikf)-xintb*djlA
           fderA(2,jat,ikf)=fderA(2,jat,ikf)-yintb*djlA
           fderA(3,jat,ikf)=fderA(3,jat,ikf)-zintb*djlA
           fderB(1,jat,ikf)=fderB(1,jat,ikf)-xintb*djlB
           fderB(2,jat,ikf)=fderB(2,jat,ikf)-yintb*djlB
           fderB(3,jat,ikf)=fderB(3,jat,ikf)-zintb*djlB
c  Exchange term 2  for atom 2
           fderA(1,jat,jlf)=fderA(1,jat,jlf)-xintb*dikA
           fderA(2,jat,jlf)=fderA(2,jat,jlf)-yintb*dikA
           fderA(3,jat,jlf)=fderA(3,jat,jlf)-zintb*dikA
           fderB(1,jat,jlf)=fderB(1,jat,jlf)-xintb*dikB
           fderB(2,jat,jlf)=fderB(2,jat,jlf)-yintb*dikB
           fderB(3,jat,jlf)=fderB(3,jat,jlf)-zintb*dikB
c  Exchange term 3  for atom 2
           fderA(1,jat,jkf)=fderA(1,jat,jkf)-xintb*dilA
           fderA(2,jat,jkf)=fderA(2,jat,jkf)-yintb*dilA
           fderA(3,jat,jkf)=fderA(3,jat,jkf)-zintb*dilA
           fderB(1,jat,jkf)=fderB(1,jat,jkf)-xintb*dilB
           fderB(2,jat,jkf)=fderB(2,jat,jkf)-yintb*dilB
           fderB(3,jat,jkf)=fderB(3,jat,jkf)-zintb*dilB
c  Exchange term 4  for atom 2
           fderA(1,jat,ilf)=fderA(1,jat,ilf)-xintb*djkA
           fderA(2,jat,ilf)=fderA(2,jat,ilf)-yintb*djkA
           fderA(3,jat,ilf)=fderA(3,jat,ilf)-zintb*djkA
           fderB(1,jat,ilf)=fderB(1,jat,ilf)-xintb*djkB
           fderB(2,jat,ilf)=fderB(2,jat,ilf)-yintb*djkB
           fderB(3,jat,ilf)=fderB(3,jat,ilf)-zintb*djkB
         endif
         if(dokat) then
c  Coulomb term 1  for atom 3
           fderA(1,kat,ijf)=fderA(1,kat,ijf)+xintc*dkl2
           fderA(2,kat,ijf)=fderA(2,kat,ijf)+yintc*dkl2
           fderA(3,kat,ijf)=fderA(3,kat,ijf)+zintc*dkl2
           fderB(1,kat,ijf)=fderB(1,kat,ijf)+xintc*dkl2
           fderB(2,kat,ijf)=fderB(2,kat,ijf)+yintc*dkl2
           fderB(3,kat,ijf)=fderB(3,kat,ijf)+zintc*dkl2
c  Coulomb term 2  for atom 3
           fderA(1,kat,klf)=fderA(1,kat,klf)+xintc*dij2
           fderA(2,kat,klf)=fderA(2,kat,klf)+yintc*dij2
           fderA(3,kat,klf)=fderA(3,kat,klf)+zintc*dij2
           fderB(1,kat,klf)=fderB(1,kat,klf)+xintc*dij2
           fderB(2,kat,klf)=fderB(2,kat,klf)+yintc*dij2
           fderB(3,kat,klf)=fderB(3,kat,klf)+zintc*dij2
c  Exchange term 1  for atom 3
           fderA(1,kat,ikf)=fderA(1,kat,ikf)-xintc*djlA
           fderA(2,kat,ikf)=fderA(2,kat,ikf)-yintc*djlA
           fderA(3,kat,ikf)=fderA(3,kat,ikf)-zintc*djlA
           fderB(1,kat,ikf)=fderB(1,kat,ikf)-xintc*djlB
           fderB(2,kat,ikf)=fderB(2,kat,ikf)-yintc*djlB
           fderB(3,kat,ikf)=fderB(3,kat,ikf)-zintc*djlB
c  Exchange term 2  for atom 3
           fderA(1,kat,jlf)=fderA(1,kat,jlf)-xintc*dikA
           fderA(2,kat,jlf)=fderA(2,kat,jlf)-yintc*dikA
           fderA(3,kat,jlf)=fderA(3,kat,jlf)-zintc*dikA
           fderB(1,kat,jlf)=fderB(1,kat,jlf)-xintc*dikB
           fderB(2,kat,jlf)=fderB(2,kat,jlf)-yintc*dikB
           fderB(3,kat,jlf)=fderB(3,kat,jlf)-zintc*dikB
c  Exchange term 3  for atom 3
           fderA(1,kat,jkf)=fderA(1,kat,jkf)-xintc*dilA
           fderA(2,kat,jkf)=fderA(2,kat,jkf)-yintc*dilA
           fderA(3,kat,jkf)=fderA(3,kat,jkf)-zintc*dilA
           fderB(1,kat,jkf)=fderB(1,kat,jkf)-xintc*dilB
           fderB(2,kat,jkf)=fderB(2,kat,jkf)-yintc*dilB
           fderB(3,kat,jkf)=fderB(3,kat,jkf)-zintc*dilB
c  Exchange term 4  for atom 3
           fderA(1,kat,ilf)=fderA(1,kat,ilf)-xintc*djkA
           fderA(2,kat,ilf)=fderA(2,kat,ilf)-yintc*djkA
           fderA(3,kat,ilf)=fderA(3,kat,ilf)-zintc*djkA
           fderB(1,kat,ilf)=fderB(1,kat,ilf)-xintc*djkB
           fderB(2,kat,ilf)=fderB(2,kat,ilf)-yintc*djkB
           fderB(3,kat,ilf)=fderB(3,kat,ilf)-zintc*djkB
         endif
         if(dolat) then
c  Coulomb term 1  for atom 4
           fderA(1,lat,ijf)=fderA(1,lat,ijf)+xintd*dkl2
           fderA(2,lat,ijf)=fderA(2,lat,ijf)+yintd*dkl2
           fderA(3,lat,ijf)=fderA(3,lat,ijf)+zintd*dkl2
           fderB(1,lat,ijf)=fderB(1,lat,ijf)+xintd*dkl2
           fderB(2,lat,ijf)=fderB(2,lat,ijf)+yintd*dkl2
           fderB(3,lat,ijf)=fderB(3,lat,ijf)+zintd*dkl2
c  Coulomb term 2  for atom 4
           fderA(1,lat,klf)=fderA(1,lat,klf)+xintd*dij2
           fderA(2,lat,klf)=fderA(2,lat,klf)+yintd*dij2
           fderA(3,lat,klf)=fderA(3,lat,klf)+zintd*dij2
           fderB(1,lat,klf)=fderB(1,lat,klf)+xintd*dij2
           fderB(2,lat,klf)=fderB(2,lat,klf)+yintd*dij2
           fderB(3,lat,klf)=fderB(3,lat,klf)+zintd*dij2
c  Exchange term 1  for atom 4
           fderA(1,lat,ikf)=fderA(1,lat,ikf)-xintd*djlA
           fderA(2,lat,ikf)=fderA(2,lat,ikf)-yintd*djlA
           fderA(3,lat,ikf)=fderA(3,lat,ikf)-zintd*djlA
           fderB(1,lat,ikf)=fderB(1,lat,ikf)-xintd*djlB
           fderB(2,lat,ikf)=fderB(2,lat,ikf)-yintd*djlB
           fderB(3,lat,ikf)=fderB(3,lat,ikf)-zintd*djlB
c  Exchange term 2  for atom 4
           fderA(1,lat,jlf)=fderA(1,lat,jlf)-xintd*dikA
           fderA(2,lat,jlf)=fderA(2,lat,jlf)-yintd*dikA
           fderA(3,lat,jlf)=fderA(3,lat,jlf)-zintd*dikA
           fderB(1,lat,jlf)=fderB(1,lat,jlf)-xintd*dikB
           fderB(2,lat,jlf)=fderB(2,lat,jlf)-yintd*dikB
           fderB(3,lat,jlf)=fderB(3,lat,jlf)-zintd*dikB
c  Exchange term 3  for atom 4
           fderA(1,lat,jkf)=fderA(1,lat,jkf)-xintd*dilA
           fderA(2,lat,jkf)=fderA(2,lat,jkf)-yintd*dilA
           fderA(3,lat,jkf)=fderA(3,lat,jkf)-zintd*dilA
           fderB(1,lat,jkf)=fderB(1,lat,jkf)-xintd*dilB
           fderB(2,lat,jkf)=fderB(2,lat,jkf)-yintd*dilB
           fderB(3,lat,jkf)=fderB(3,lat,jkf)-zintd*dilB
c  Exchange term 4  for atom 4
           fderA(1,lat,ilf)=fderA(1,lat,ilf)-xintd*djkA
           fderA(2,lat,ilf)=fderA(2,lat,ilf)-yintd*djkA
           fderA(3,lat,ilf)=fderA(3,lat,ilf)-zintd*djkA
           fderB(1,lat,ilf)=fderB(1,lat,ilf)-xintd*djkB
           fderB(2,lat,ilf)=fderB(2,lat,ilf)-yintd*djkB
           fderB(3,lat,ilf)=fderB(3,lat,ilf)-zintd*djkB
         endif
c
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine fderN_uhfC(buf,denA,denB,fderA,fderB,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,nsym,natb,nate)
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
c----------------------------------------------------------------
c
c  ** COULOMB ONLY **
c
c INPUT :
c  buf()    - integral derivatives buffer
c  denA()   - alpha density matrix
c  denB()   - beta density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c STORAGE:
c labels    - space for the integral indices
c length(4) - the sizes of the 4 shells making up the shell integral.
c   E.g. if the 5 shells are of type L, P, D6 and S, then length=4,3,6,1
c lgenct    - it does not appear to be used. Somethig to do with general
c   contractions?
c na        - number of atoms
c nsym      - number of non-trivial symmetry operations (0 for C1 sym.)
c  natb     - fisrt componet to be computed in current pass
c  nate     - last componet to be computed in current pass
c
c OUTPUT :
c fderA() - alpha Fock matrix derivatives
c fderB() - beta Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat,dojat,dokat,dolat
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension fderA(3,natb:nate,ntri),fderB(3,natb:nate,ntri)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      data zero,half,one,two,four /0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from fock_deriv() :NBLS=',nbls,' ngcd=',ngcd
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
c-------------------------------------------------------
          doiat=.false.
          if(iat.ge.natb .and. iat.le.nate) doiat=.true.
          dojat=.false.
          if(jat.ge.natb .and. jat.le.nate) dojat=.true.
          dokat=.false.
          if(kat.ge.natb .and. kat.le.nate) dokat=.true.
          dolat=.false.
          if(lat.ge.natb .and. lat.le.nate) dolat=.true.
c-------------------------------------------------------
      if(.not.doiat .and. .not.dojat .and.
     *   .not.dokat .and. .not.dolat       ) go to 150
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
                dij2=(denA(ijf)+denB(ijf))*two
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      klf=kk+lcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl2=(denA(klf)+denB(klf))*two
c
                      integ=integ+1
c------------------------------------------
c  the 12 derivatives
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
                      xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
                      yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
                      zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
         if(doiat) then
c  Coulomb term 1  for atom 1
           fderA(1,iat,ijf)=fderA(1,iat,ijf)+xinta*dkl2
           fderA(2,iat,ijf)=fderA(2,iat,ijf)+yinta*dkl2
           fderA(3,iat,ijf)=fderA(3,iat,ijf)+zinta*dkl2
           fderB(1,iat,ijf)=fderB(1,iat,ijf)+xinta*dkl2
           fderB(2,iat,ijf)=fderB(2,iat,ijf)+yinta*dkl2
           fderB(3,iat,ijf)=fderB(3,iat,ijf)+zinta*dkl2
c  Coulomb term 2  for atom 1
           fderA(1,iat,klf)=fderA(1,iat,klf)+xinta*dij2
           fderA(2,iat,klf)=fderA(2,iat,klf)+yinta*dij2
           fderA(3,iat,klf)=fderA(3,iat,klf)+zinta*dij2
           fderB(1,iat,klf)=fderB(1,iat,klf)+xinta*dij2
           fderB(2,iat,klf)=fderB(2,iat,klf)+yinta*dij2
           fderB(3,iat,klf)=fderB(3,iat,klf)+zinta*dij2
         endif
         if(dojat) then
c  Coulomb term 1  for atom 2
           fderA(1,jat,ijf)=fderA(1,jat,ijf)+xintb*dkl2
           fderA(2,jat,ijf)=fderA(2,jat,ijf)+yintb*dkl2
           fderA(3,jat,ijf)=fderA(3,jat,ijf)+zintb*dkl2
           fderB(1,jat,ijf)=fderB(1,jat,ijf)+xintb*dkl2
           fderB(2,jat,ijf)=fderB(2,jat,ijf)+yintb*dkl2
           fderB(3,jat,ijf)=fderB(3,jat,ijf)+zintb*dkl2
c  Coulomb term 2  for atom 2
           fderA(1,jat,klf)=fderA(1,jat,klf)+xintb*dij2
           fderA(2,jat,klf)=fderA(2,jat,klf)+yintb*dij2
           fderA(3,jat,klf)=fderA(3,jat,klf)+zintb*dij2
           fderB(1,jat,klf)=fderB(1,jat,klf)+xintb*dij2
           fderB(2,jat,klf)=fderB(2,jat,klf)+yintb*dij2
           fderB(3,jat,klf)=fderB(3,jat,klf)+zintb*dij2
         endif
         if(dokat) then
c  Coulomb term 1  for atom 3
           fderA(1,kat,ijf)=fderA(1,kat,ijf)+xintc*dkl2
           fderA(2,kat,ijf)=fderA(2,kat,ijf)+yintc*dkl2
           fderA(3,kat,ijf)=fderA(3,kat,ijf)+zintc*dkl2
           fderB(1,kat,ijf)=fderB(1,kat,ijf)+xintc*dkl2
           fderB(2,kat,ijf)=fderB(2,kat,ijf)+yintc*dkl2
           fderB(3,kat,ijf)=fderB(3,kat,ijf)+zintc*dkl2
c  Coulomb term 2  for atom 3
           fderA(1,kat,klf)=fderA(1,kat,klf)+xintc*dij2
           fderA(2,kat,klf)=fderA(2,kat,klf)+yintc*dij2
           fderA(3,kat,klf)=fderA(3,kat,klf)+zintc*dij2
           fderB(1,kat,klf)=fderB(1,kat,klf)+xintc*dij2
           fderB(2,kat,klf)=fderB(2,kat,klf)+yintc*dij2
           fderB(3,kat,klf)=fderB(3,kat,klf)+zintc*dij2
         endif
         if(dolat) then
c  Coulomb term 1  for atom 4
           fderA(1,lat,ijf)=fderA(1,lat,ijf)+xintd*dkl2
           fderA(2,lat,ijf)=fderA(2,lat,ijf)+yintd*dkl2
           fderA(3,lat,ijf)=fderA(3,lat,ijf)+zintd*dkl2
           fderB(1,lat,ijf)=fderB(1,lat,ijf)+xintd*dkl2
           fderB(2,lat,ijf)=fderB(2,lat,ijf)+yintd*dkl2
           fderB(3,lat,ijf)=fderB(3,lat,ijf)+zintd*dkl2
c  Coulomb term 2  for atom 4
           fderA(1,lat,klf)=fderA(1,lat,klf)+xintd*dij2
           fderA(2,lat,klf)=fderA(2,lat,klf)+yintd*dij2
           fderA(3,lat,klf)=fderA(3,lat,klf)+zintd*dij2
           fderB(1,lat,klf)=fderB(1,lat,klf)+xintd*dij2
           fderB(2,lat,klf)=fderB(2,lat,klf)+yintd*dij2
           fderB(3,lat,klf)=fderB(3,lat,klf)+zintd*dij2
         endif
c
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
  150   continue
  100 continue
c
      end
