c=====================================================================
      subroutine int_d0g1(idft,ax,rhf,nblocks,bl,inx,ntri,thref1,
     *                    natb,nate,listreal,do_allat,
     *                    densp , ncenter,
     *                    dens,denB,fder,fderB, labels,mywork,igran)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c
c (1) the first deriv. integ. are contracted with the zero-rder density
c     to form 2-el. part of the first-order fock matrices G(D0,g1)
c---------------------------------------------------------------------
c This routine is called from para_d0g1 (called from hess2e1)
c for all three modes : non -, semi-, and full-direct.
c Integral derivatives g1=(ij|kl)(Xn,Yn,Zn) are calculated only once
c and are contracted with corresponding density to form Fock derivative
c matrix
c---------------------------------------------------------------------
c INPUT :
c  idft              -  dft flag:  0 - no dft
c  ax                - factor to multiply exchange contribution for ACM
c  rhf               - logical flag for closed-shell
c  nblocks           - number of blocks of contracted shell quartets
c  bl                - scratch for everything
c  inx(12,ncs)       - basis set & nuclear data
c  ntri              - ncf*(ncf+1)/2
c  thref1            - integral threshold for G(D0,g1)
c  densp             - screening density
c  ncenter           - cenetr of each basis function
c  dens(ntri)        - alpha/closed-shell density matrix
c  denB(ntri)        - beta density matrix
c  labels            - labels buffer
c  mywork            - for parallel runs
c  igran             - for parallel runs
C OUTPUT:
c fder(3,natoms,ntri)- derivative fock matrices
c hessian(3,natoms,3,natoms) - hessian matrix
c---------------------------------------------------------------------
      logical rhf,moreint,stopnow,do_allat
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
      dimension bl(*),dens(*),denB(*)
      dimension densp(*)
      dimension ncenter(*)
      dimension listreal(*)
      dimension labels(*)
      dimension fder(*),fderB(*)
c-----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c get print flag
c
      call getival('printh',nprint)
c
      call getival('nsym',nsym)
c----------------------------------------------------------------
      where='forc'
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
      if(.not.rhf) then
         call do_intder1all(idft,ax,rhf,istart,istop,
     *                      bl,inx,thref1,    densp ,labels,
     *                      dens,denB,fder,fderB,
     *                      natb,nate,
     *                      bl(lind),ntri,   ncenter ,na,nsym )
      endif
cccc  if(nsym.eq.0) then
      if(rhf) then
      if(do_allat) then
         call do_intder1all(idft,ax,rhf,istart,istop,
     *                      bl,inx,thref1,    densp ,labels,
     *                      dens,denB,fder,fderB,
     *                      natb,nate,
     *                      bl(lind),ntri,   ncenter ,na,nsym )
      else
         call getival('SymFunPr',ifp)
         call getival('SymNuPr1',nupair)
         call getival('nsyo',nsyo)
         call do_intder1unq(idft,ax,rhf,istart,istop,
     *                      bl,inx,thref1,    densp ,labels,
     *                      dens,denB,fder,fderB,
     *                      natb,nate,listreal,
     *                      bl(lind),ntri,   ncenter ,na,
     *                      bl(nsyo),bl(ifp),bl(nupair),nsym)
      endif
      endif
c----------------------------------------------------------------
c Print derivative fock matrices
c
      if(nprint.ge.3) then
         if(rhf)then
           write(6,*)' Derivative Fock after G(D0,g1)          '
           call fder_print(fder,na,ntri,thref1)
         else
           write(6,*)' Alpha derivative Fock after G(D0,g1)          '
           call fder_print(fder,na,ntri,thref1)
           write(6,*)' Beta derivative Fock after G(D0,g1)          '
           call fder_print(fderB,na,ntri,thref1)
         endif
      endif
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
      call secund(time2)
      call term_info(thref1,time2-time1,time1-time1,where)
c
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine prtn_der1(buf,
     *                     lind,ntri,nbls,ngcd,lnijkl,
     *                     ncenter,labels,length,lgenct,
     *                     na,nsym )
c----------------------------------------------------------------
c INPUT :
c  buf()    - integral Ist derivatives buffer
c  dens()   - density matrix
c  lind()   - precomputed diagonals (  i*(i-1)/2 )
c  ntri     - ncf*(ncf+1)/2
c  nbls     - size of the integral derivative block (no of quartets)
c  ngcd     - product of general contractions of 4 shells
c  lnijkl   - number of integrals from one quartet
c  ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c OUTPUT :
c atforces() - forces at atoms
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      data zero /0.d0/
c----------------------------------------------------------------
      write(6,*)'    --------------------------------------'
      write(6,*)'    First  derivatives of 2-el. integrals '
      write(6,*)'    --------------------------------------'
      write(6,*)
     *'                  d/dA          d/dB         d/dC         d/dD'
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
c------------------------------------------
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
      x10a=1.d-10*xinta
      x10b=1.d-10*xintb
      x10c=1.d-10*xintc
      x10d=1.d-10*xintd
c
      y10a=1.d-10*yinta
      y10b=1.d-10*yintb
      y10c=1.d-10*yintc
      y10d=1.d-10*yintd
c
      z10a=1.d-10*zinta
      z10b=1.d-10*zintb
      z10c=1.d-10*zintc
      z10d=1.d-10*zintd
c     write(6,66) icf,jcf,kcf,lcf, x10a,x10b,x10c,x10d
c    *                           , y10a,y10b,y10c,y10d
c    *                           , z10a,z10b,z10c,z10d
c 66  format(4(i2,1x),3x,4(f12.7,2x)/15x,4(f12.7,2x)/15x,4(f12.7,2x))
      write(6,66) icf,jcf,kcf,lcf, x10a,x10b,x10c,x10d,' d/dx',
     *                             y10a,y10b,y10c,y10d,' d/dy',
     *                             z10a,z10b,z10c,z10d,' d/dz'
  66  format(4(i2,1x),1x,4(f12.6,1x),1x,a5/
     *               13x,4(f12.6,1x),1x,a5/
     *               13x,4(f12.6,1x),1x,a5)
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine do_intder1all(idft,ax,rhf,istart,istop,
     *                         bl,inx,thres1,densp,labels,
     *                         dens,denB,fder,fderB,
     *                         natb,nate,
     *                         lind,ntri,ncenter,na,nsym)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
      logical rhf,moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
      common /runtype/ scftype,where
c
      parameter (Zero=0.0d0)
c---------------------------------------------------------------------
      dimension inx(12,*)
      dimension labels(*)
      dimension lind(*)
      dimension ncenter(*)
      dimension bl(*)
      dimension densp(*)
      dimension dens(*),denB(*)
      dimension fder(*),fderB(*)
c----------------------------------------------------------------
      call getival('printh',nprint)
c----------------------------------------------------------------
c calculate the gradient integrals for Fder (needed for hessian)
c
c
      where='forc'
      screen='2pde'
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
c
c  nintez is the size of a given quartet.It is set up
c  to zero if there is no integrals
c
      IF(nintez.gt.0) THEN
        call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
c
c calculate two-el. contrib. to derivative fock metrices :
c
        IF(idft.EQ.0) THEN
c  ...........................
c   Standard HF
c  ...........................
          if(nprint.ge.5) then
c          print 1st der.integ:
           call prtn_der1(bl(ibuffz),lind,ntri,nblsiz,ngctoz,
     *	   nintez,ncenter,labels(lab1),labels(lab2),labels(lab3),
     *       na,nsym)
          endif
c
          If(rhf) Then
           call fderNrhf(bl(ibuffz),dens,fder,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,nsym,natb,nate)
          Else
           call fderNuhf(bl(ibuffz),dens,denB,fder,fderB,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,nsym,natb,nate)
          EndIf
cc
        ELSE IF(ax.NE.Zero) THEN
c  ...........................
c   Hybrid HF-DFT (ACM)
c  ...........................
          If(rhf) Then
           call fderN_racm(bl(ibuffz),dens,fder,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,nsy,ax,natb,nate)
          Else
           call fderN_uacm(bl(ibuffz),dens,denB,fder,fderB,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,nsy,ax,natb,nate)
          EndIf
cc
        ELSE
c  ...........................
c   "Pure" DFT (NO EXCHANGE)
c  ...........................
cc
          If(rhf) Then
           call fderN_rhfC(bl(ibuffz),dens,fder,lind,
     *                 ntri,nblsiz,ngctoz,nintez,ncenter,
     *                 labels(lab1),labels(lab2),labels(lab3),
     *                 na,nsym,natb,nate)
          Else
           call fderN_uhfC(bl(ibuffz),dens,denB,fder,fderB,lind,
     *                 ntri,nblsiz,ngctoz,nintez,ncenter,
     *                 labels(lab1),labels(lab2),labels(lab3),
     *                 na,nsym,natb,nate)
          EndIf
        ENDIF
      ENDIF
c
ctest
c          write(6,*)' Derivative Fock after block ',isupbl
c          call fder_print(fder,na,ntri,thres1)
ctest
        if(moreint) go to 11
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
      end
c=====================================================================
      subroutine fderNrhf(buf,dens ,fder,lind,
     *                  ntri,nbls,ngcd,lnijkl,ncenter,
     *                  labels,length,lgenct,na,nsym,natb,nate )
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
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
c fder() - Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat,dojat,dokat,dolat
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
cccc  dimension fder(3,na,ntri)
      dimension fder(3,natb:nate,ntri)
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
                dij4=dens(ijf)*four
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
                      dkl4=dens(klf)*four
                      dil=dens(ilf)
                      djl=dens(jlf)
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
c--------------------------------------------------------
         if(doiat) then
c  Coulomb term 1  for atom 1
           fder(1,iat,ijf)=fder(1,iat,ijf)+xinta*dkl4
           fder(2,iat,ijf)=fder(2,iat,ijf)+yinta*dkl4
           fder(3,iat,ijf)=fder(3,iat,ijf)+zinta*dkl4
c  Coulomb term 2  for atom 1
           fder(1,iat,klf)=fder(1,iat,klf)+xinta*dij4
           fder(2,iat,klf)=fder(2,iat,klf)+yinta*dij4
           fder(3,iat,klf)=fder(3,iat,klf)+zinta*dij4
c  Exchange term 1  for atom 1
           fder(1,iat,ikf)=fder(1,iat,ikf)-xinta*djl
           fder(2,iat,ikf)=fder(2,iat,ikf)-yinta*djl
           fder(3,iat,ikf)=fder(3,iat,ikf)-zinta*djl
c  Exchange term 2  for atom 1
           fder(1,iat,jlf)=fder(1,iat,jlf)-xinta*dik
           fder(2,iat,jlf)=fder(2,iat,jlf)-yinta*dik
           fder(3,iat,jlf)=fder(3,iat,jlf)-zinta*dik
c  Exchange term 3  for atom 1
           fder(1,iat,jkf)=fder(1,iat,jkf)-xinta*dil
           fder(2,iat,jkf)=fder(2,iat,jkf)-yinta*dil
           fder(3,iat,jkf)=fder(3,iat,jkf)-zinta*dil
c  Exchange term 4  for atom 1
           fder(1,iat,ilf)=fder(1,iat,ilf)-xinta*djk
           fder(2,iat,ilf)=fder(2,iat,ilf)-yinta*djk
           fder(3,iat,ilf)=fder(3,iat,ilf)-zinta*djk
         endif   ! (doiat) then
c--------------------------------------------------------
c--------------------------------------------------------
         if(dojat) then
c  Coulomb term 1  for atom 2
           fder(1,jat,ijf)=fder(1,jat,ijf)+xintb*dkl4
           fder(2,jat,ijf)=fder(2,jat,ijf)+yintb*dkl4
           fder(3,jat,ijf)=fder(3,jat,ijf)+zintb*dkl4
c  Coulomb term 2  for atom 2
           fder(1,jat,klf)=fder(1,jat,klf)+xintb*dij4
           fder(2,jat,klf)=fder(2,jat,klf)+yintb*dij4
           fder(3,jat,klf)=fder(3,jat,klf)+zintb*dij4
c  Exchange term 1  for atom 2
           fder(1,jat,ikf)=fder(1,jat,ikf)-xintb*djl
           fder(2,jat,ikf)=fder(2,jat,ikf)-yintb*djl
           fder(3,jat,ikf)=fder(3,jat,ikf)-zintb*djl
c  Exchange term 2  for atom 2
           fder(1,jat,jlf)=fder(1,jat,jlf)-xintb*dik
           fder(2,jat,jlf)=fder(2,jat,jlf)-yintb*dik
           fder(3,jat,jlf)=fder(3,jat,jlf)-zintb*dik
c  Exchange term 3  for atom 2
           fder(1,jat,jkf)=fder(1,jat,jkf)-xintb*dil
           fder(2,jat,jkf)=fder(2,jat,jkf)-yintb*dil
           fder(3,jat,jkf)=fder(3,jat,jkf)-zintb*dil
c  Exchange term 4  for atom 2
           fder(1,jat,ilf)=fder(1,jat,ilf)-xintb*djk
           fder(2,jat,ilf)=fder(2,jat,ilf)-yintb*djk
           fder(3,jat,ilf)=fder(3,jat,ilf)-zintb*djk
         endif  ! (dojat) then
c--------------------------------------------------------
c--------------------------------------------------------
         if(dokat) then
c  Coulomb term 1  for atom 3
           fder(1,kat,ijf)=fder(1,kat,ijf)+xintc*dkl4
           fder(2,kat,ijf)=fder(2,kat,ijf)+yintc*dkl4
           fder(3,kat,ijf)=fder(3,kat,ijf)+zintc*dkl4
c  Coulomb term 2  for atom 3
           fder(1,kat,klf)=fder(1,kat,klf)+xintc*dij4
           fder(2,kat,klf)=fder(2,kat,klf)+yintc*dij4
           fder(3,kat,klf)=fder(3,kat,klf)+zintc*dij4
c  Exchange term 1  for atom 3
           fder(1,kat,ikf)=fder(1,kat,ikf)-xintc*djl
           fder(2,kat,ikf)=fder(2,kat,ikf)-yintc*djl
           fder(3,kat,ikf)=fder(3,kat,ikf)-zintc*djl
c  Exchange term 2  for atom 3
           fder(1,kat,jlf)=fder(1,kat,jlf)-xintc*dik
           fder(2,kat,jlf)=fder(2,kat,jlf)-yintc*dik
           fder(3,kat,jlf)=fder(3,kat,jlf)-zintc*dik
c  Exchange term 3  for atom 3
           fder(1,kat,jkf)=fder(1,kat,jkf)-xintc*dil
           fder(2,kat,jkf)=fder(2,kat,jkf)-yintc*dil
           fder(3,kat,jkf)=fder(3,kat,jkf)-zintc*dil
c  Exchange term 4  for atom 3
           fder(1,kat,ilf)=fder(1,kat,ilf)-xintc*djk
           fder(2,kat,ilf)=fder(2,kat,ilf)-yintc*djk
           fder(3,kat,ilf)=fder(3,kat,ilf)-zintc*djk
         endif   ! (dokat) then
c--------------------------------------------------------
c--------------------------------------------------------
         if(dolat) then
c  Coulomb term 1  for atom 4
           fder(1,lat,ijf)=fder(1,lat,ijf)+xintd*dkl4
           fder(2,lat,ijf)=fder(2,lat,ijf)+yintd*dkl4
           fder(3,lat,ijf)=fder(3,lat,ijf)+zintd*dkl4
c  Coulomb term 2  for atom 4
           fder(1,lat,klf)=fder(1,lat,klf)+xintd*dij4
           fder(2,lat,klf)=fder(2,lat,klf)+yintd*dij4
           fder(3,lat,klf)=fder(3,lat,klf)+zintd*dij4
c  Exchange term 1  for atom 4
           fder(1,lat,ikf)=fder(1,lat,ikf)-xintd*djl
           fder(2,lat,ikf)=fder(2,lat,ikf)-yintd*djl
           fder(3,lat,ikf)=fder(3,lat,ikf)-zintd*djl
c  Exchange term 2  for atom 4
           fder(1,lat,jlf)=fder(1,lat,jlf)-xintd*dik
           fder(2,lat,jlf)=fder(2,lat,jlf)-yintd*dik
           fder(3,lat,jlf)=fder(3,lat,jlf)-zintd*dik
c  Exchange term 3  for atom 4
           fder(1,lat,jkf)=fder(1,lat,jkf)-xintd*dil
           fder(2,lat,jkf)=fder(2,lat,jkf)-yintd*dil
           fder(3,lat,jkf)=fder(3,lat,jkf)-zintd*dil
c  Exchange term 4  for atom 4
           fder(1,lat,ilf)=fder(1,lat,ilf)-xintd*djk
           fder(2,lat,ilf)=fder(2,lat,ilf)-yintd*djk
           fder(3,lat,ilf)=fder(3,lat,ilf)-zintd*djk
         endif  ! (dolat) then
c--------------------------------------------------------
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
      subroutine fderNuhf(buf,denA,denB,fderA,fderB,lind,
     *                  ntri,nbls,ngcd,lnijkl,ncenter,
     *                  labels,length,lgenct,na,nsym,natb,nate )
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
         if(doiat) then
           fderA(1,iat,ijf)=fderA(1,iat,ijf)+xinta*dkl2
           fderA(2,iat,ijf)=fderA(2,iat,ijf)+yinta*dkl2
           fderA(3,iat,ijf)=fderA(3,iat,ijf)+zinta*dkl2
           fderB(1,iat,ijf)=fderB(1,iat,ijf)+xinta*dkl2
           fderB(2,iat,ijf)=fderB(2,iat,ijf)+yinta*dkl2
           fderB(3,iat,ijf)=fderB(3,iat,ijf)+zinta*dkl2
           fderA(1,iat,klf)=fderA(1,iat,klf)+xinta*dij2
           fderA(2,iat,klf)=fderA(2,iat,klf)+yinta*dij2
           fderA(3,iat,klf)=fderA(3,iat,klf)+zinta*dij2
           fderB(1,iat,klf)=fderB(1,iat,klf)+xinta*dij2
           fderB(2,iat,klf)=fderB(2,iat,klf)+yinta*dij2
           fderB(3,iat,klf)=fderB(3,iat,klf)+zinta*dij2
           fderA(1,iat,ikf)=fderA(1,iat,ikf)-xinta*djlA
           fderA(2,iat,ikf)=fderA(2,iat,ikf)-yinta*djlA
           fderA(3,iat,ikf)=fderA(3,iat,ikf)-zinta*djlA
           fderB(1,iat,ikf)=fderB(1,iat,ikf)-xinta*djlB
           fderB(2,iat,ikf)=fderB(2,iat,ikf)-yinta*djlB
           fderB(3,iat,ikf)=fderB(3,iat,ikf)-zinta*djlB
           fderA(1,iat,jlf)=fderA(1,iat,jlf)-xinta*dikA
           fderA(2,iat,jlf)=fderA(2,iat,jlf)-yinta*dikA
           fderA(3,iat,jlf)=fderA(3,iat,jlf)-zinta*dikA
           fderB(1,iat,jlf)=fderB(1,iat,jlf)-xinta*dikB
           fderB(2,iat,jlf)=fderB(2,iat,jlf)-yinta*dikB
           fderB(3,iat,jlf)=fderB(3,iat,jlf)-zinta*dikB
           fderA(1,iat,jkf)=fderA(1,iat,jkf)-xinta*dilA
           fderA(2,iat,jkf)=fderA(2,iat,jkf)-yinta*dilA
           fderA(3,iat,jkf)=fderA(3,iat,jkf)-zinta*dilA
           fderB(1,iat,jkf)=fderB(1,iat,jkf)-xinta*dilB
           fderB(2,iat,jkf)=fderB(2,iat,jkf)-yinta*dilB
           fderB(3,iat,jkf)=fderB(3,iat,jkf)-zinta*dilB
           fderA(1,iat,ilf)=fderA(1,iat,ilf)-xinta*djkA
           fderA(2,iat,ilf)=fderA(2,iat,ilf)-yinta*djkA
           fderA(3,iat,ilf)=fderA(3,iat,ilf)-zinta*djkA
           fderB(1,iat,ilf)=fderB(1,iat,ilf)-xinta*djkB
           fderB(2,iat,ilf)=fderB(2,iat,ilf)-yinta*djkB
           fderB(3,iat,ilf)=fderB(3,iat,ilf)-zinta*djkB
         endif  ! (doiat) then
c--------------------------------------------------------
         if(dojat) then
           fderA(1,jat,ijf)=fderA(1,jat,ijf)+xintb*dkl2
           fderA(2,jat,ijf)=fderA(2,jat,ijf)+yintb*dkl2
           fderA(3,jat,ijf)=fderA(3,jat,ijf)+zintb*dkl2
           fderB(1,jat,ijf)=fderB(1,jat,ijf)+xintb*dkl2
           fderB(2,jat,ijf)=fderB(2,jat,ijf)+yintb*dkl2
           fderB(3,jat,ijf)=fderB(3,jat,ijf)+zintb*dkl2
           fderA(1,jat,klf)=fderA(1,jat,klf)+xintb*dij2
           fderA(2,jat,klf)=fderA(2,jat,klf)+yintb*dij2
           fderA(3,jat,klf)=fderA(3,jat,klf)+zintb*dij2
           fderB(1,jat,klf)=fderB(1,jat,klf)+xintb*dij2
           fderB(2,jat,klf)=fderB(2,jat,klf)+yintb*dij2
           fderB(3,jat,klf)=fderB(3,jat,klf)+zintb*dij2
           fderA(1,jat,ikf)=fderA(1,jat,ikf)-xintb*djlA
           fderA(2,jat,ikf)=fderA(2,jat,ikf)-yintb*djlA
           fderA(3,jat,ikf)=fderA(3,jat,ikf)-zintb*djlA
           fderB(1,jat,ikf)=fderB(1,jat,ikf)-xintb*djlB
           fderB(2,jat,ikf)=fderB(2,jat,ikf)-yintb*djlB
           fderB(3,jat,ikf)=fderB(3,jat,ikf)-zintb*djlB
           fderA(1,jat,jlf)=fderA(1,jat,jlf)-xintb*dikA
           fderA(2,jat,jlf)=fderA(2,jat,jlf)-yintb*dikA
           fderA(3,jat,jlf)=fderA(3,jat,jlf)-zintb*dikA
           fderB(1,jat,jlf)=fderB(1,jat,jlf)-xintb*dikB
           fderB(2,jat,jlf)=fderB(2,jat,jlf)-yintb*dikB
           fderB(3,jat,jlf)=fderB(3,jat,jlf)-zintb*dikB
           fderA(1,jat,jkf)=fderA(1,jat,jkf)-xintb*dilA
           fderA(2,jat,jkf)=fderA(2,jat,jkf)-yintb*dilA
           fderA(3,jat,jkf)=fderA(3,jat,jkf)-zintb*dilA
           fderB(1,jat,jkf)=fderB(1,jat,jkf)-xintb*dilB
           fderB(2,jat,jkf)=fderB(2,jat,jkf)-yintb*dilB
           fderB(3,jat,jkf)=fderB(3,jat,jkf)-zintb*dilB
           fderA(1,jat,ilf)=fderA(1,jat,ilf)-xintb*djkA
           fderA(2,jat,ilf)=fderA(2,jat,ilf)-yintb*djkA
           fderA(3,jat,ilf)=fderA(3,jat,ilf)-zintb*djkA
           fderB(1,jat,ilf)=fderB(1,jat,ilf)-xintb*djkB
           fderB(2,jat,ilf)=fderB(2,jat,ilf)-yintb*djkB
           fderB(3,jat,ilf)=fderB(3,jat,ilf)-zintb*djkB
         endif   ! (dojat) then
c--------------------------------------------------------
         if(dokat) then
           fderA(1,kat,ijf)=fderA(1,kat,ijf)+xintc*dkl2
           fderA(2,kat,ijf)=fderA(2,kat,ijf)+yintc*dkl2
           fderA(3,kat,ijf)=fderA(3,kat,ijf)+zintc*dkl2
           fderB(1,kat,ijf)=fderB(1,kat,ijf)+xintc*dkl2
           fderB(2,kat,ijf)=fderB(2,kat,ijf)+yintc*dkl2
           fderB(3,kat,ijf)=fderB(3,kat,ijf)+zintc*dkl2
           fderA(1,kat,klf)=fderA(1,kat,klf)+xintc*dij2
           fderA(2,kat,klf)=fderA(2,kat,klf)+yintc*dij2
           fderA(3,kat,klf)=fderA(3,kat,klf)+zintc*dij2
           fderB(1,kat,klf)=fderB(1,kat,klf)+xintc*dij2
           fderB(2,kat,klf)=fderB(2,kat,klf)+yintc*dij2
           fderB(3,kat,klf)=fderB(3,kat,klf)+zintc*dij2
           fderA(1,kat,ikf)=fderA(1,kat,ikf)-xintc*djlA
           fderA(2,kat,ikf)=fderA(2,kat,ikf)-yintc*djlA
           fderA(3,kat,ikf)=fderA(3,kat,ikf)-zintc*djlA
           fderB(1,kat,ikf)=fderB(1,kat,ikf)-xintc*djlB
           fderB(2,kat,ikf)=fderB(2,kat,ikf)-yintc*djlB
           fderB(3,kat,ikf)=fderB(3,kat,ikf)-zintc*djlB
           fderA(1,kat,jlf)=fderA(1,kat,jlf)-xintc*dikA
           fderA(2,kat,jlf)=fderA(2,kat,jlf)-yintc*dikA
           fderA(3,kat,jlf)=fderA(3,kat,jlf)-zintc*dikA
           fderB(1,kat,jlf)=fderB(1,kat,jlf)-xintc*dikB
           fderB(2,kat,jlf)=fderB(2,kat,jlf)-yintc*dikB
           fderB(3,kat,jlf)=fderB(3,kat,jlf)-zintc*dikB
           fderA(1,kat,jkf)=fderA(1,kat,jkf)-xintc*dilA
           fderA(2,kat,jkf)=fderA(2,kat,jkf)-yintc*dilA
           fderA(3,kat,jkf)=fderA(3,kat,jkf)-zintc*dilA
           fderB(1,kat,jkf)=fderB(1,kat,jkf)-xintc*dilB
           fderB(2,kat,jkf)=fderB(2,kat,jkf)-yintc*dilB
           fderB(3,kat,jkf)=fderB(3,kat,jkf)-zintc*dilB
           fderA(1,kat,ilf)=fderA(1,kat,ilf)-xintc*djkA
           fderA(2,kat,ilf)=fderA(2,kat,ilf)-yintc*djkA
           fderA(3,kat,ilf)=fderA(3,kat,ilf)-zintc*djkA
           fderB(1,kat,ilf)=fderB(1,kat,ilf)-xintc*djkB
           fderB(2,kat,ilf)=fderB(2,kat,ilf)-yintc*djkB
           fderB(3,kat,ilf)=fderB(3,kat,ilf)-zintc*djkB
         endif   !  (dokat) then
c--------------------------------------------------------
         if(dolat) then
           fderA(1,lat,ijf)=fderA(1,lat,ijf)+xintd*dkl2
           fderA(2,lat,ijf)=fderA(2,lat,ijf)+yintd*dkl2
           fderA(3,lat,ijf)=fderA(3,lat,ijf)+zintd*dkl2
           fderB(1,lat,ijf)=fderB(1,lat,ijf)+xintd*dkl2
           fderB(2,lat,ijf)=fderB(2,lat,ijf)+yintd*dkl2
           fderB(3,lat,ijf)=fderB(3,lat,ijf)+zintd*dkl2
           fderA(1,lat,klf)=fderA(1,lat,klf)+xintd*dij2
           fderA(2,lat,klf)=fderA(2,lat,klf)+yintd*dij2
           fderA(3,lat,klf)=fderA(3,lat,klf)+zintd*dij2
           fderB(1,lat,klf)=fderB(1,lat,klf)+xintd*dij2
           fderB(2,lat,klf)=fderB(2,lat,klf)+yintd*dij2
           fderB(3,lat,klf)=fderB(3,lat,klf)+zintd*dij2
           fderA(1,lat,ikf)=fderA(1,lat,ikf)-xintd*djlA
           fderA(2,lat,ikf)=fderA(2,lat,ikf)-yintd*djlA
           fderA(3,lat,ikf)=fderA(3,lat,ikf)-zintd*djlA
           fderB(1,lat,ikf)=fderB(1,lat,ikf)-xintd*djlB
           fderB(2,lat,ikf)=fderB(2,lat,ikf)-yintd*djlB
           fderB(3,lat,ikf)=fderB(3,lat,ikf)-zintd*djlB
           fderA(1,lat,jlf)=fderA(1,lat,jlf)-xintd*dikA
           fderA(2,lat,jlf)=fderA(2,lat,jlf)-yintd*dikA
           fderA(3,lat,jlf)=fderA(3,lat,jlf)-zintd*dikA
           fderB(1,lat,jlf)=fderB(1,lat,jlf)-xintd*dikB
           fderB(2,lat,jlf)=fderB(2,lat,jlf)-yintd*dikB
           fderB(3,lat,jlf)=fderB(3,lat,jlf)-zintd*dikB
           fderA(1,lat,jkf)=fderA(1,lat,jkf)-xintd*dilA
           fderA(2,lat,jkf)=fderA(2,lat,jkf)-yintd*dilA
           fderA(3,lat,jkf)=fderA(3,lat,jkf)-zintd*dilA
           fderB(1,lat,jkf)=fderB(1,lat,jkf)-xintd*dilB
           fderB(2,lat,jkf)=fderB(2,lat,jkf)-yintd*dilB
           fderB(3,lat,jkf)=fderB(3,lat,jkf)-zintd*dilB
           fderA(1,lat,ilf)=fderA(1,lat,ilf)-xintd*djkA
           fderA(2,lat,ilf)=fderA(2,lat,ilf)-yintd*djkA
           fderA(3,lat,ilf)=fderA(3,lat,ilf)-zintd*djkA
           fderB(1,lat,ilf)=fderB(1,lat,ilf)-xintd*djkB
           fderB(2,lat,ilf)=fderB(2,lat,ilf)-yintd*djkB
           fderB(3,lat,ilf)=fderB(3,lat,ilf)-zintd*djkB
         endif  !  (dolat) then
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
      subroutine fderN_racm(buf,dens ,fder,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,nsym,ax,natb,nate)
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
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
c fder() - Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat,dojat,dokat,dolat
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension fder(3,natb:nate,ntri)
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
                dij4=dens(ijf)*four
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                   ikf=ii+kcf
                   jkf=jj+kcf
                   if(kcf.gt.icf) ikf=kk+icf
                   if(kcf.gt.jcf) jkf=kk+jcf
                   dik=dens(ikf)*ax
                   djk=dens(jkf)*ax
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      ilf=ii+lcf
                      jlf=jj+lcf
                      klf=kk+lcf
                      if(lcf.gt.icf) ilf=ll+icf
                      if(lcf.gt.jcf) jlf=ll+jcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl4=dens(klf)*four
                      dil=dens(ilf)*ax
                      djl=dens(jlf)*ax
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
           fder(1,iat,ijf)=fder(1,iat,ijf)+xinta*dkl4
           fder(2,iat,ijf)=fder(2,iat,ijf)+yinta*dkl4
           fder(3,iat,ijf)=fder(3,iat,ijf)+zinta*dkl4
c  Coulomb term 2  for atom 1
           fder(1,iat,klf)=fder(1,iat,klf)+xinta*dij4
           fder(2,iat,klf)=fder(2,iat,klf)+yinta*dij4
           fder(3,iat,klf)=fder(3,iat,klf)+zinta*dij4
c  Exchange term 1  for atom 1
           fder(1,iat,ikf)=fder(1,iat,ikf)-xinta*djl
           fder(2,iat,ikf)=fder(2,iat,ikf)-yinta*djl
           fder(3,iat,ikf)=fder(3,iat,ikf)-zinta*djl
c  Exchange term 2  for atom 1
           fder(1,iat,jlf)=fder(1,iat,jlf)-xinta*dik
           fder(2,iat,jlf)=fder(2,iat,jlf)-yinta*dik
           fder(3,iat,jlf)=fder(3,iat,jlf)-zinta*dik
c  Exchange term 3  for atom 1
           fder(1,iat,jkf)=fder(1,iat,jkf)-xinta*dil
           fder(2,iat,jkf)=fder(2,iat,jkf)-yinta*dil
           fder(3,iat,jkf)=fder(3,iat,jkf)-zinta*dil
c  Exchange term 4  for atom 1
           fder(1,iat,ilf)=fder(1,iat,ilf)-xinta*djk
           fder(2,iat,ilf)=fder(2,iat,ilf)-yinta*djk
           fder(3,iat,ilf)=fder(3,iat,ilf)-zinta*djk
         endif
         if(dojat) then
c  Coulomb term 1  for atom 2
           fder(1,jat,ijf)=fder(1,jat,ijf)+xintb*dkl4
           fder(2,jat,ijf)=fder(2,jat,ijf)+yintb*dkl4
           fder(3,jat,ijf)=fder(3,jat,ijf)+zintb*dkl4
c  Coulomb term 2  for atom 2
           fder(1,jat,klf)=fder(1,jat,klf)+xintb*dij4
           fder(2,jat,klf)=fder(2,jat,klf)+yintb*dij4
           fder(3,jat,klf)=fder(3,jat,klf)+zintb*dij4
c  Exchange term 1  for atom 2
           fder(1,jat,ikf)=fder(1,jat,ikf)-xintb*djl
           fder(2,jat,ikf)=fder(2,jat,ikf)-yintb*djl
           fder(3,jat,ikf)=fder(3,jat,ikf)-zintb*djl
c  Exchange term 2  for atom 2
           fder(1,jat,jlf)=fder(1,jat,jlf)-xintb*dik
           fder(2,jat,jlf)=fder(2,jat,jlf)-yintb*dik
           fder(3,jat,jlf)=fder(3,jat,jlf)-zintb*dik
c  Exchange term 3  for atom 2
           fder(1,jat,jkf)=fder(1,jat,jkf)-xintb*dil
           fder(2,jat,jkf)=fder(2,jat,jkf)-yintb*dil
           fder(3,jat,jkf)=fder(3,jat,jkf)-zintb*dil
c  Exchange term 4  for atom 2
           fder(1,jat,ilf)=fder(1,jat,ilf)-xintb*djk
           fder(2,jat,ilf)=fder(2,jat,ilf)-yintb*djk
           fder(3,jat,ilf)=fder(3,jat,ilf)-zintb*djk
         endif
         if(dokat) then
c  Coulomb term 1  for atom 3
           fder(1,kat,ijf)=fder(1,kat,ijf)+xintc*dkl4
           fder(2,kat,ijf)=fder(2,kat,ijf)+yintc*dkl4
           fder(3,kat,ijf)=fder(3,kat,ijf)+zintc*dkl4
c  Coulomb term 2  for atom 3
           fder(1,kat,klf)=fder(1,kat,klf)+xintc*dij4
           fder(2,kat,klf)=fder(2,kat,klf)+yintc*dij4
           fder(3,kat,klf)=fder(3,kat,klf)+zintc*dij4
c  Exchange term 1  for atom 3
           fder(1,kat,ikf)=fder(1,kat,ikf)-xintc*djl
           fder(2,kat,ikf)=fder(2,kat,ikf)-yintc*djl
           fder(3,kat,ikf)=fder(3,kat,ikf)-zintc*djl
c  Exchange term 2  for atom 3
           fder(1,kat,jlf)=fder(1,kat,jlf)-xintc*dik
           fder(2,kat,jlf)=fder(2,kat,jlf)-yintc*dik
           fder(3,kat,jlf)=fder(3,kat,jlf)-zintc*dik
c  Exchange term 3  for atom 3
           fder(1,kat,jkf)=fder(1,kat,jkf)-xintc*dil
           fder(2,kat,jkf)=fder(2,kat,jkf)-yintc*dil
           fder(3,kat,jkf)=fder(3,kat,jkf)-zintc*dil
c  Exchange term 4  for atom 3
           fder(1,kat,ilf)=fder(1,kat,ilf)-xintc*djk
           fder(2,kat,ilf)=fder(2,kat,ilf)-yintc*djk
           fder(3,kat,ilf)=fder(3,kat,ilf)-zintc*djk
         endif
         if(dolat) then
c  Coulomb term 1  for atom 4
           fder(1,lat,ijf)=fder(1,lat,ijf)+xintd*dkl4
           fder(2,lat,ijf)=fder(2,lat,ijf)+yintd*dkl4
           fder(3,lat,ijf)=fder(3,lat,ijf)+zintd*dkl4
c  Coulomb term 2  for atom 4
           fder(1,lat,klf)=fder(1,lat,klf)+xintd*dij4
           fder(2,lat,klf)=fder(2,lat,klf)+yintd*dij4
           fder(3,lat,klf)=fder(3,lat,klf)+zintd*dij4
c  Exchange term 1  for atom 4
           fder(1,lat,ikf)=fder(1,lat,ikf)-xintd*djl
           fder(2,lat,ikf)=fder(2,lat,ikf)-yintd*djl
           fder(3,lat,ikf)=fder(3,lat,ikf)-zintd*djl
c  Exchange term 2  for atom 4
           fder(1,lat,jlf)=fder(1,lat,jlf)-xintd*dik
           fder(2,lat,jlf)=fder(2,lat,jlf)-yintd*dik
           fder(3,lat,jlf)=fder(3,lat,jlf)-zintd*dik
c  Exchange term 3  for atom 4
           fder(1,lat,jkf)=fder(1,lat,jkf)-xintd*dil
           fder(2,lat,jkf)=fder(2,lat,jkf)-yintd*dil
           fder(3,lat,jkf)=fder(3,lat,jkf)-zintd*dil
c  Exchange term 4  for atom 4
           fder(1,lat,ilf)=fder(1,lat,ilf)-xintd*djk
           fder(2,lat,ilf)=fder(2,lat,ilf)-yintd*djk
           fder(3,lat,ilf)=fder(3,lat,ilf)-zintd*djk
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
      subroutine fderN_rhfC(buf,dens ,fder,lind,
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
c  dens()   - density matrix
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
c fder() - Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat,dojat,dokat,dolat
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension fder(3,natb:nate,ntri)
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
                dij4=dens(ijf)*four
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      klf=kk+lcf
                      if(lcf.gt.kcf) klf=ll+kcf
                      dkl4=dens(klf)*four
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
           fder(1,iat,ijf)=fder(1,iat,ijf)+xinta*dkl4
           fder(2,iat,ijf)=fder(2,iat,ijf)+yinta*dkl4
           fder(3,iat,ijf)=fder(3,iat,ijf)+zinta*dkl4
c  Coulomb term 2  for atom 1
           fder(1,iat,klf)=fder(1,iat,klf)+xinta*dij4
           fder(2,iat,klf)=fder(2,iat,klf)+yinta*dij4
           fder(3,iat,klf)=fder(3,iat,klf)+zinta*dij4
         endif
         if(dojat) then
c  Coulomb term 1  for atom 2
           fder(1,jat,ijf)=fder(1,jat,ijf)+xintb*dkl4
           fder(2,jat,ijf)=fder(2,jat,ijf)+yintb*dkl4
           fder(3,jat,ijf)=fder(3,jat,ijf)+zintb*dkl4
c  Coulomb term 2  for atom 2
           fder(1,jat,klf)=fder(1,jat,klf)+xintb*dij4
           fder(2,jat,klf)=fder(2,jat,klf)+yintb*dij4
           fder(3,jat,klf)=fder(3,jat,klf)+zintb*dij4
         endif
         if(dokat) then
c  Coulomb term 1  for atom 3
           fder(1,kat,ijf)=fder(1,kat,ijf)+xintc*dkl4
           fder(2,kat,ijf)=fder(2,kat,ijf)+yintc*dkl4
           fder(3,kat,ijf)=fder(3,kat,ijf)+zintc*dkl4
c  Coulomb term 2  for atom 3
           fder(1,kat,klf)=fder(1,kat,klf)+xintc*dij4
           fder(2,kat,klf)=fder(2,kat,klf)+yintc*dij4
           fder(3,kat,klf)=fder(3,kat,klf)+zintc*dij4
         endif
         if(dolat) then
c  Coulomb term 1  for atom 4
           fder(1,lat,ijf)=fder(1,lat,ijf)+xintd*dkl4
           fder(2,lat,ijf)=fder(2,lat,ijf)+yintd*dkl4
           fder(3,lat,ijf)=fder(3,lat,ijf)+zintd*dkl4
c  Coulomb term 2  for atom 4
           fder(1,lat,klf)=fder(1,lat,klf)+xintd*dij4
           fder(2,lat,klf)=fder(2,lat,klf)+yintd*dij4
           fder(3,lat,klf)=fder(3,lat,klf)+zintd*dij4
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
      subroutine do_intder1unq(idft,ax,rhf,istart,istop,
     *                         bl,inx,thres1,densp,labels,
     *                         dens,denB,fder,fderB,
     *                         natb,nate,listreal,
     *                         lind,ntri,ncenter,na,
     *                         nsyo,ifp,nupair,nsym)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
      logical rhf,moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
      common /runtype/ scftype,where
c
      parameter (Zero=0.0d0)
c---------------------------------------------------------------------
      dimension inx(12,*)
      dimension labels(*)
      dimension lind(*)
      dimension ncenter(*)
      dimension listreal(*)
      dimension bl(*)
      dimension densp(*)
      dimension dens(*),denB(*)
      dimension fder(*),fderB(*)
      dimension ifp(7,*),nsyo(7)
      dimension nupair(na,*)           ! (natoms,nsym)
      dimension mirror(3,7)
      dimension ngxyz(3,7)
      data mirror/-1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c----------------------------------------------------------------
      do 10 ns=1,nsym
         nop=nsyo(ns)
         do icr=1,3
            ngxyz(icr,ns)=mirror(icr,nop)
         enddo
   10 continue
c----------------------------------------------------------------
      call getival('printh',nprint)
c----------------------------------------------------------------
c calculate the gradient integrals for Fder (needed for hessian)
c
      where='forc'
      screen='2pde'
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
c
c  nintez is the size of a given quartet.It is set up
c  to zero if there is no integrals
c
      IF(nintez.gt.0) THEN
        call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
c
c calculate two-el. contrib. to derivative fock metrices :
c
        IF(idft.EQ.0) THEN
c  ...........................
c   Standard HF
c  ...........................
          if(nprint.ge.5) then
c          print 1st der.integ:
           call prtn_der1(bl(ibuffz),lind,ntri,nblsiz,ngctoz,
     *	   nintez,ncenter,labels(lab1),labels(lab2),labels(lab3),
     *       na,nsym)
          endif
c
          If(rhf) Then
           call fderSrhf(bl(ibuffz),dens,fder,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,natb,nate,listreal,
     *                nsyo,ifp,nupair,nsym,ngxyz)
          Else
           call fderSuhf(bl(ibuffz),dens,denB,fder,fderB,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,ifp,nsym,natb,nate)
          EndIf
cc
        ELSE IF(ax.NE.Zero) THEN
c  ...........................
c   Hybrid HF-DFT (ACM)
c  ...........................
          If(rhf) Then
           call fderS_racm(bl(ibuffz),dens,fder,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,ifp,nsym,ax,natb,nate)
          Else
           call fderS_uacm(bl(ibuffz),dens,denB,fder,fderB,lind,
     *                ntri,nblsiz,ngctoz,nintez,ncenter,
     *                labels(lab1),labels(lab2),labels(lab3),
     *                na,ifp,nsym,ax,natb,nate)
          EndIf
cc
        ELSE
c  ...........................
c   "Pure" DFT (NO EXCHANGE)
c  ...........................
cc
          If(rhf) Then
           call fderS_rhfC(bl(ibuffz),dens,fder,lind,
     *                 ntri,nblsiz,ngctoz,nintez,ncenter,
     *                 labels(lab1),labels(lab2),labels(lab3),
     *                 na,ifp,nsym,natb,nate)
          Else
           call fderS_uhfC(bl(ibuffz),dens,denB,fder,fderB,lind,
     *                 ntri,nblsiz,ngctoz,nintez,ncenter,
     *                 labels(lab1),labels(lab2),labels(lab3),
     *                 na,ifp,nsym,natb,nate)
          EndIf
        ENDIF
      ENDIF
ctest
c          write(6,*)' Derivative Fock after block ',isupbl
cccc       call fder_print(fder,na,ntri,thres1)
c          call fder_print(fder,2 ,ntri,thres1)
ctest
c
        if(moreint) go to 11
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
      end
c=====================================================================
c=====================================================================
c========================SYMMETRY=====================================
c=====================================================================
c=====================================================================
      subroutine fderSrhf(buf,dens ,fder,lind,
     *                  ntri,nbls,ngcd,lnijkl,ncenter,
     *                  labels,length,lgenct,
     *                  na,natb,nate,listreal,
     *                  nsyo,ifp,nupair,nsym,ngxyz)
c----------------------------------------------------------------
c  Fock matrix derivatives G(D0,g1) for BASIC SCF
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
c fder() - Fock matrix derivatives
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat,dojat,dokat,dolat
      logical doiat1,dojat1,dokat1,dolat1
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
cccc  dimension fder(3,na,ntri)
      dimension fder(3,natb:nate,ntri)
      dimension lind(*) , ncenter(*)
      dimension listreal(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension ifp(7,*)
      dimension nsyo(7)
      dimension nupair(na,*)           ! (natoms,nsym)
      dimension xyzA(3),xyzB(3),xyzC(3),xyzD(3)
      dimension ngxyz(3,7)
      data zero,half,one,two,four /0.0d0,0.5d0,1.0d0,2.0d0,4.0d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from fock_deriv() :NBLS=',nbls,' ngcd=',ngcd
c     write(6,*)' NatB,NatE=',natb,nate
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
c check if they are uniqe atoms
c find what unique atoms these real (unique) atoms are
c
          iau=listreal(iat)
          jau=listreal(jat)
          kau=listreal(kat)
          lau=listreal(lat)
c
c notice: if iat is not a unique atom then iau=0
c-------------------------------------------------------
          doiat=.false.
          if(iau.ge.natb .and. iau.le.nate) doiat=.true.
          dojat=.false.
          if(jau.ge.natb .and. jau.le.nate) dojat=.true.
          dokat=.false.
          if(kau.ge.natb .and. kau.le.nate) dokat=.true.
          dolat=.false.
          if(lau.ge.natb .and. lau.le.nate) dolat=.true.
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
                dij4=dens(ijf)*four
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
                      dkl4=dens(klf)*four
                      dil=dens(ilf)
                      djl=dens(jlf)
c
                      integ=integ+1
c------------------------------------------
c  the 12 derivatives
c                     xinta=buf(1,ijklp,integ,iqu)
c                     xintb=buf(2,ijklp,integ,iqu)
c                     xintc=buf(3,ijklp,integ,iqu)
c                     xintd=-(xinta+xintb+xintc) ! trans. inv.
c                     yinta=buf(4,ijklp,integ,iqu)
c                     yintb=buf(5,ijklp,integ,iqu)
c                     yintc=buf(6,ijklp,integ,iqu)
c                     yintd=-(yinta+yintb+yintc) ! trans. inv.
c                     zinta=buf(7,ijklp,integ,iqu)
c                     zintb=buf(8,ijklp,integ,iqu)
c                     zintc=buf(9,ijklp,integ,iqu)
c                     zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
                      xyzA(1)=buf(1,ijklp,integ,iqu)
                      xyzA(2)=buf(4,ijklp,integ,iqu)
                      xyzA(3)=buf(7,ijklp,integ,iqu)
c
                      xyzB(1)=buf(2,ijklp,integ,iqu)
                      xyzB(2)=buf(5,ijklp,integ,iqu)
                      xyzB(3)=buf(8,ijklp,integ,iqu)
c
                      xyzC(1)=buf(3,ijklp,integ,iqu)
                      xyzC(2)=buf(6,ijklp,integ,iqu)
                      xyzC(3)=buf(9,ijklp,integ,iqu)
c
                      xyzD(1)=-(xyzA(1)+xyzB(1)+xyzC(1))
                      xyzD(2)=-(xyzA(2)+xyzB(2)+xyzC(2))
                      xyzD(3)=-(xyzA(3)+xyzB(3)+xyzC(3))
c--------------------------------------------------------
         if(doiat) then
           do icr=1,3
             fder(icr,iau,ijf)=fder(icr,iau,ijf)+xyzA(icr)*dkl4
             fder(icr,iau,klf)=fder(icr,iau,klf)+xyzA(icr)*dij4
             fder(icr,iau,ikf)=fder(icr,iau,ikf)-xyzA(icr)*djl
             fder(icr,iau,jlf)=fder(icr,iau,jlf)-xyzA(icr)*dik
             fder(icr,iau,jkf)=fder(icr,iau,jkf)-xyzA(icr)*dil
             fder(icr,iau,ilf)=fder(icr,iau,ilf)-xyzA(icr)*djk
           enddo
         endif   ! (doiat) then
c--------------------------------------------------------
         if(dojat) then
           do icr=1,3
             fder(icr,jau,ijf)=fder(icr,jau,ijf)+xyzB(icr)*dkl4
             fder(icr,jau,klf)=fder(icr,jau,klf)+xyzB(icr)*dij4
             fder(icr,jau,ikf)=fder(icr,jau,ikf)-xyzB(icr)*djl
             fder(icr,jau,jlf)=fder(icr,jau,jlf)-xyzB(icr)*dik
             fder(icr,jau,jkf)=fder(icr,jau,jkf)-xyzB(icr)*dil
             fder(icr,jau,ilf)=fder(icr,jau,ilf)-xyzB(icr)*djk
           enddo
         endif  ! (dojat) then
c--------------------------------------------------------
         if(dokat) then
           do icr=1,3
             fder(icr,kau,ijf)=fder(icr,kau,ijf)+xyzC(icr)*dkl4
             fder(icr,kau,klf)=fder(icr,kau,klf)+xyzC(icr)*dij4
             fder(icr,kau,ikf)=fder(icr,kau,ikf)-xyzC(icr)*djl
             fder(icr,kau,jlf)=fder(icr,kau,jlf)-xyzC(icr)*dik
             fder(icr,kau,jkf)=fder(icr,kau,jkf)-xyzC(icr)*dil
             fder(icr,kau,ilf)=fder(icr,kau,ilf)-xyzC(icr)*djk
           enddo
         endif   ! (dokat) then
c--------------------------------------------------------
         if(dolat) then
           do icr=1,3
             fder(icr,lau,ijf)=fder(icr,lau,ijf)+xyzD(icr)*dkl4
             fder(icr,lau,klf)=fder(icr,lau,klf)+xyzD(icr)*dij4
             fder(icr,lau,ikf)=fder(icr,lau,ikf)-xyzD(icr)*djl
             fder(icr,lau,jlf)=fder(icr,lau,jlf)-xyzD(icr)*dik
             fder(icr,lau,jkf)=fder(icr,lau,jkf)-xyzD(icr)*dil
             fder(icr,lau,ilf)=fder(icr,lau,ilf)-xyzD(icr)*djk
           enddo
         endif  ! (dolat) then
c--------------------------------------------------------
CCCCCCC    call fder_print(fder,na,ntri,1.d-10)
c          call fder_print(fder,nate-natb+1,ntri,1.d-10)
c--------------------------------------------------------
c symmetry
c........................................................
       do ns=1,nsym
          iat1=nupair(iat,ns)
          jat1=nupair(jat,ns)
          kat1=nupair(kat,ns)
          lat1=nupair(lat,ns)
c-------------------------------------------------------
c check if they are uniqe atoms
c find what unique atoms these real (unique) atoms are
c notice: if iat1 is not a unique atom then iau1=0
c
          iau1=listreal(iat1)
          jau1=listreal(jat1)
          kau1=listreal(kat1)
          lau1=listreal(lat1)
ccc
      if(iau1.eq.0 .and. jau1.eq.0 .and.
     *   kau1.eq.0 .and. lau1.eq.0      ) go to 7654
c-------------------------------------------------------
          doiat1=.false.
          if(iau1.ge.natb .and. iau1.le.nate) doiat1=.true.
          dojat1=.false.
          if(jau1.ge.natb .and. jau1.le.nate) dojat1=.true.
          dokat1=.false.
          if(kau1.ge.natb .and. kau1.le.nate) dokat1=.true.
          dolat1=.false.
          if(lau1.ge.natb .and. lau1.le.nate) dolat1=.true.
c-------------------------------------------------------
      if(.not.doiat1.and. .not.dojat1.and.
     *   .not.dokat1.and. .not.dolat1      ) go to 7654
c-------------------------------------------------------
          ic1=ifp(ns,icf)
          jc1=ifp(ns,jcf)
          kc1=ifp(ns,kcf)
          lc1=ifp(ns,lcf)
c
          fct=one
          if(ic1.lt.0) then
             ic1=-ic1
             fct=-fct
          endif
          if(jc1.lt.0) then
             jc1=-jc1
             fct=-fct
          endif
          if(kc1.lt.0) then
             kc1=-kc1
             fct=-fct
          endif
          if(lc1.lt.0) then
             lc1=-lc1
             fct=-fct
          endif
c
          ii1=lind(ic1)
          jj1=lind(jc1)
          kk1=lind(kc1)
          ll1=lind(lc1)
c
c ....................Coulomb  part ....................
          if(ic1.ge.jc1) then
            ijf1 = ii1+jc1
          else
            ijf1 = jj1+ic1
          endif
          if(kc1.ge.lc1) then
            klf1 = kk1+lc1
          else
            klf1 = ll1+kc1
          endif
c ....................Exchange part ....................
          if(ic1.ge.kc1) then
            ikf1 = ii1+kc1
          else
            ikf1 = kk1+ic1
          endif
          if(ic1.ge.lc1) then
            ilf1 = ii1+lc1
          else
            ilf1 = ll1+ic1
          endif
          if(jc1.ge.kc1) then
            jkf1 = jj1+kc1
          else
            jkf1 = kk1+jc1
          endif
          if(jc1.ge.lc1) then
            jlf1 = jj1+lc1
          else
            jlf1 = ll1+jc1
          endif
c
          fkl4=fct*dens(klf1)*four
          fij4=fct*dens(ijf1)*four
          fjl =fct*dens(jlf1)
          fik =fct*dens(ikf1)
          fil =fct*dens(ilf1)
          fjk =fct*dens(jkf1)
c-------------------------------------------------------
         if(doiat1) then
           do icr=1,3
             xyzAs=xyzA(icr)*ngxyz(icr,ns)
             fder(icr,iau1,ijf1)=fder(icr,iau1,ijf1)+xyzAs*fkl4
             fder(icr,iau1,klf1)=fder(icr,iau1,klf1)+xyzAs*fij4
             fder(icr,iau1,ikf1)=fder(icr,iau1,ikf1)-xyzAs*fjl
             fder(icr,iau1,jlf1)=fder(icr,iau1,jlf1)-xyzAs*fik
             fder(icr,iau1,jkf1)=fder(icr,iau1,jkf1)-xyzAs*fil
             fder(icr,iau1,ilf1)=fder(icr,iau1,ilf1)-xyzAs*fjk
           enddo
         endif   ! (doiat) then
c--------------------------------------------------------
         if(dojat1) then
           do icr=1,3
             xyzBs=xyzB(icr)*ngxyz(icr,ns)
             fder(icr,jau1,ijf1)=fder(icr,jau1,ijf1)+xyzBs*fkl4
             fder(icr,jau1,klf1)=fder(icr,jau1,klf1)+xyzBs*fij4
             fder(icr,jau1,ikf1)=fder(icr,jau1,ikf1)-xyzBs*fjl
             fder(icr,jau1,jlf1)=fder(icr,jau1,jlf1)-xyzBs*fik
             fder(icr,jau1,jkf1)=fder(icr,jau1,jkf1)-xyzBs*fil
             fder(icr,jau1,ilf1)=fder(icr,jau1,ilf1)-xyzBs*fjk
           enddo
         endif  ! (dojat) then
c
c--------------------------------------------------------
         if(dokat1) then
           do icr=1,3
             xyzCs=xyzC(icr)*ngxyz(icr,ns)
             fder(icr,kau1,ijf1)=fder(icr,kau1,ijf1)+xyzCs*fkl4
             fder(icr,kau1,klf1)=fder(icr,kau1,klf1)+xyzCs*fij4
             fder(icr,kau1,ikf1)=fder(icr,kau1,ikf1)-xyzCs*fjl
             fder(icr,kau1,jlf1)=fder(icr,kau1,jlf1)-xyzCs*fik
             fder(icr,kau1,jkf1)=fder(icr,kau1,jkf1)-xyzCs*fil
             fder(icr,kau1,ilf1)=fder(icr,kau1,ilf1)-xyzCs*fjk
           enddo
         endif   ! (dokat) then
c
c--------------------------------------------------------
         if(dolat1) then
           do icr=1,3
             xyzDs=xyzD(icr)*ngxyz(icr,ns)
             fder(icr,lau1,ijf1)=fder(icr,lau1,ijf1)+xyzDs*fkl4
             fder(icr,lau1,klf1)=fder(icr,lau1,klf1)+xyzDs*fij4
             fder(icr,lau1,ikf1)=fder(icr,lau1,ikf1)-xyzDs*fjl
             fder(icr,lau1,jlf1)=fder(icr,lau1,jlf1)-xyzDs*fik
             fder(icr,lau1,jkf1)=fder(icr,lau1,jkf1)-xyzDs*fil
             fder(icr,lau1,ilf1)=fder(icr,lau1,ilf1)-xyzDs*fjk
           enddo
         endif  ! (dolat) then
c
c--------------------------------------------------------
ccccc      call fder_print(fder,na,ntri,1.d-10)
c          call fder_print(fder,nate-natb+1,ntri,1.d-10)
 7654  continue
       enddo !     ns=1,nsym
c symmetry end
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
c NOT READY YET
      subroutine fderSuhf(buf,denA,denB,fderA,fderB,lind,
     *                  ntri,nbls,ngcd,lnijkl,ncenter,
     *                  labels,length,lgenct,na,ifp,
     *                  nsym,natb,nate)
      end
c=====================================================================
c NOT READY YET
      subroutine fderS_racm(buf,dens ,fder,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,ifp,
     *                     nsym,ax,natb,nate)
      end
c=====================================================================
c NOT READY YET
      subroutine fderS_rhfC(buf,dens ,fder,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,ifp,
     *                     nsym,natb,nate)
      end
c=====================================================================
c NOT READY YET
      subroutine fderS_uacm(buf,dens ,fder,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,ifp,
     *                     nsym,ax,natb,nate)
      end
c=====================================================================
c NOT READY YET
      subroutine fderS_uhfC(buf,dens ,fder,lind,
     *                     ntri,nbls,ngcd,lnijkl,ncenter,
     *                     labels,length,lgenct,na,nsym,natb,nate)
      end
c=====================================================================
