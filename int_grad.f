      subroutine int_grad(idft,ax,rhf,nblocks,bl,inx,ntri,thres1,
     *                    dens,denB,atforces,labels,mywork,igran)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c It calculates the two-electron inegral derivatives for gradient
c and contracts them with corresponding two-particle density matrix
c elements to form FORCES on atoms
c---------------------------------------------------------------------
c This routine is called from PARA_GRAD (called from FORCE2)
c for all three modes : non -, semi-, and full-direct.
c Integral derivatives g1=(ij|kl)(Xn,Yn,Zn) are calculated only once
c and are contracted with corresponding two-particle density
c
c value of where='forc' is set up HERE
c---------------------------------------------------------------------
c INPUT :
c  idft              -  dft flag:  0 - no dft
c  ax                - factor to multiply exchange contribution for ACM
c  rhf               - logical flag for closed-shell
c  nblocks           - number of blocks of contracted shell quartets
c  bl                - scratch for everything
c  inx(12,ncs)       - basis set & nuclear data
c  ntri              - ncf*(ncf+1)/2
c  thres1            - integral threshold
c  dens(ntri)        - alpha/closed-shell density matrix
c  denB(ntri)        - beta density matrix
c  labels            - labels buffer
c  mywork            - for parallel runs
c  igran             - for parallel runs
C OUTPUT:
c atforces(3,natoms) - forces on atoms
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
      dimension bl(*),dens(*),denB(*)
      dimension labels(*)
      dimension atforces(*)
c
      parameter (Zero=0.0d0)
c-----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c calculate the GRADIENT integrals
c
      where='forc'
      screen='forc'
c----------------------------------------------------------------
c for symmetrization of two-el. part of forces/atoms for exact symm.
c
      call getival('nsym',nsym)
      if(nsym.gt.0) then
         call getival('nsyo',nsyo)
         call make_negxyz(nsym,bl(nsyo))
      endif
cc        if(nsym.gt.0) call make_negxyz(nsym,nsyo)
c----------------------------------------------------------------
c set integral timings to zero as well as neglect stuff
c
        call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
        call tstival('ftc0',iyes)
        if(iyes.eq.0) call setival('ftc0',0)
        call getival('ftc0',iftc0)
        If(iftc0.NE.0) Then
           call getival('nftcbl4',nftcbl4)  ! pointer
        endif
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
c transform ordinary density dens(ij) into denspar(ics,jcs)
c
      call getmem(ncs*ncs,idensp)
      call getmem(ncf,map_fs)
      If(rhf) Then
        call setup_densp2(inx,ncf,ncs,dens,bl(idensp),
     *                   bl(map_fs))
      Else
c -- for UHF, allocate extra array for alpha+beta densities
        call getmem(ntri,idAB)
csep    Call AddVEC(ntri,dens,denB,bl(idAB))
c we need rather sum of ABSOLUTE values here for prescreening :
        call addvec_abs(ntri,dens,denB,bl(idAB))
        call vscal(ntri, 0.5d0 , bl(idab) )
        call setup_densp2(inx,ncf,ncs,bl(idAB),bl(idensp),
     *                   bl(map_fs))
        call retmem(1)
C **WARNING ** Could use alpha+beta density in <force_uhf>
      EndIf
      call retmem(1)     ! allocated for map_fs
c----------------------------------------------------------------
c print density :
cc    call druma(dens,ncf, 6  ,'density ')
c----------------------------------------------------------------
c allocate memory for mapping array : ncenter(icf)->iatom
c
      call getint(ncf,ncenter)
c
c integer allocation in bl
c----------------------------------------------------------------
c get centers of basis functions :
c
      call get_center(ncs,ncf,inx,bl(ncenter))
c
c----------------------------------------------------------------
c loop over blocks of contracted shell quartets
c
CCCC  DO ISUPBL=ISTART,ISTOP
      DO isup  =istart,istop
c
         if(iftc0.EQ.0) then
           isupbl=isup
         else
           call get_block4(isup,bl(nftcbl4),isupbl)
         endif
c
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,bl(idensp),where,
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
c calculate forces on atoms :
c
        IF(idft.EQ.0) THEN
c  ...........................
c   Standard HF
c  ...........................
          If(rhf) Then
           call force_ona(bl(ibuffz),dens,atforces,bl(lind),
     *                    ntri,nblsiz,ngctoz,nintez,bl(ncenter),
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    na,nsym)
          Else
           call force_uhf(bl(ibuffz),dens,denB,atforces,bl(lind),
     *                    ntri,nblsiz,ngctoz,nintez,bl(ncenter),
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    na,nsym)
          EndIf
cc
        ELSE IF(ax.NE.Zero) THEN
c  ...........................
c   Hybrid HF-DFT (ACM)
c  ...........................
          If(rhf) Then
           call force_acm(ax,bl(ibuffz),dens,atforces,bl(lind),
     *                    ntri,nblsiz,ngctoz,nintez,bl(ncenter),
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    na,nsym)
          Else
           call force_uacm(ax,bl(ibuffz),dens,denB,atforces,bl(lind),
     *                     ntri,nblsiz,ngctoz,nintez,bl(ncenter),
     *                     labels(lab1),labels(lab2),labels(lab3),
     *                     na,nsym)
          EndIf
cc
        ELSE
c  ...........................
c   "Pure" DFT (NO EXCHANGE)
c  ...........................
cc
          If(rhf) Then
           call force_dft(bl(ibuffz),dens,atforces,bl(lind),
     *                    ntri,nblsiz,ngctoz,nintez,bl(ncenter),
     *                    labels(lab1),labels(lab2),labels(lab3),
     *                    na,nsym)
          Else
           call force_udft(bl(ibuffz),dens,denB,atforces,bl(lind),
     *                     ntri,nblsiz,ngctoz,nintez,bl(ncenter),
     *                     labels(lab1),labels(lab2),labels(lab3),
     *                     na,nsym)
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
c release memory allocated for ncenter(icf)->iatom
      call retmem(1)
c
c release memory allocated for idensp
      call retmem(1)
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
      call secund(time2)
        call term_info(thres1,time2-time1,time1-time1,where)
c
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine get_center(ncs,ncf,inx,ncenter)
      dimension inx(12,*)
      dimension ncenter(ncf)      ! output
c
      icff=0
      do ics=1,ncs
         iat=inx(2,ics)
         len=inx(3,ics)
         lgc=inx(4,ics)
         do igc=0,lgc
            do icf=1,len
               icff=icff+1
               ncenter(icff)=iat
            enddo
         enddo
      enddo
c test
      if(icff.ne.ncf) then
        call nerror(1,'get_center','icff .ne. ncf',icff,ncf)
      endif
c
      end
c=====================================================================
c  Forces on atoms for BASIC SCF
c ------------------------------
c
      subroutine force_ona(buf,dens ,atforces,
     *                     lind,ntri,nbls,ngcd,lnijkl,
     *                     ncenter,labels,length,lgenct,
     *                     na,nsym )
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
c atforces() - forces at atoms
c----------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension atforces(3,*)             ! (3,natoms)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension xyz(3,4)
      dimension vec(9)
      data zero /0.d0/
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' from forces_ona() :NBLS=',nbls,' ngcd=',ngcd
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
cccc      call zeroit(xyz,12)
          call zeroit(vec,9)
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
                      dijkl=dij4*dkl - dik*djl -dil*djk
c
                      integ=integ+1
c-----------------------------------------------------------------------
cslower              call daxpy(9,dijkl,buf(1,ijklp,integ,iqu),1,vec,1)
c
                    do iel=1,9
                       vec(iel)=vec(iel)+dijkl*buf(iel,ijklp,integ,iqu)
                    enddo
c-----------------------------------------------------------------------
c for print only
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
c
c     x10a=1.d-10*xinta
c     x10b=1.d-10*xintb
c     x10c=1.d-10*xintc
c     x10d=-(x10a+x10b+x10c)
c
c     y10a=1.d-10*yinta
c     y10b=1.d-10*yintb
c     y10c=1.d-10*yintc
c     y10d=-(y10a+y10b+y10c)
c
c     z10a=1.d-10*zinta
c     z10b=1.d-10*zintb
c     z10c=1.d-10*zintc
c     z10d=-(z10a+z10b+z10c)
c     write(6,66) icf,jcf,kcf,lcf, x10a,x10b,x10c,x10d
c    *                           , y10a,y10b,y10c,y10d
c    *                           , z10a,z10b,z10c,z10d
c 66  format(4(i2,1x),3x,4(f12.7,2x)/15x,4(f12.7,2x)/15x,4(f12.7,2x))
c-----------------------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
             xyz(1,1)=-vec(1)
             xyz(1,2)=-vec(2)
             xyz(1,3)=-vec(3)
             xyz(2,1)=-vec(4)
             xyz(2,2)=-vec(5)
             xyz(2,3)=-vec(6)
             xyz(3,1)=-vec(7)
             xyz(3,2)=-vec(8)
             xyz(3,3)=-vec(9)
c--------------------------------------------------------
c translational invariance :
c
             xyz(1,4)=-xyz(1,1)-xyz(1,2)-xyz(1,3)
             xyz(2,4)=-xyz(2,1)-xyz(2,2)-xyz(2,3)
             xyz(3,4)=-xyz(3,1)-xyz(3,2)-xyz(3,3)
c
             do icr=1,3
                atforces(icr,iat)=atforces(icr,iat)+xyz(icr,1)
                atforces(icr,jat)=atforces(icr,jat)+xyz(icr,2)
                atforces(icr,kat)=atforces(icr,kat)+xyz(icr,3)
                atforces(icr,lat)=atforces(icr,lat)+xyz(icr,4)
             enddo
             if(nsym.gt.0) then
                call getival('SymNuPr',nupair)
                call atforce_symm(na,nsym,bl(nupair),iat,jat,kat,lat,
     *                            xyz,atforces)
             endif
c--------------------------------------------------------
  150   continue
  100 continue
c--------------------------------------------------------
ctest
c     write(6,*)' atoms &  forces'
c     do iat=1,na
c
c        fx3=1.d-10*atforces(1,iat)
c        fy3=1.d-10*atforces(2,iat)
c        fz3=1.d-10*atforces(3,iat)
c
c        write(6,222) iat, fx3,fy3,fz3
c
c     enddo
c 222 format(1x,i3,2x,3(f9.6,1x) )
ctest
      end
c=====================================================================
c  Forces on atoms for ACM DFT
c ----------------------------
c
      subroutine force_acm(ax,buf,dens ,atforces,
     *                     lind,ntri,nbls,ngcd,lnijkl,
     *                     ncenter,labels,length,lgenct,
     *                     na,nsym )
c----------------------------------------------------------------
c INPUT :
c ax       - factor to multiply exchange term
c buf()    - integral derivatives buffer
c dens()   - density matrix
c lind()   - precomputed diagonals (  i*(i-1)/2 )
c ntri     - ncf*(ncf+1)/2
c nbls     - size of the integral derivative block (no of quartets)
c ngcd     - product of general contractions of 4 shells
c lnijkl   - number of integrals from one quartet
c ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c OUTPUT :
c atforces() - forces at atoms
c----------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension atforces(3,*)             ! (3,natoms)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension xyz(3,4)
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
          call zeroit(xyz,12)
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
c
                      dil=dens(ilf)
                      djl=dens(jlf)
                      dkl=dens(klf)
c
                      dijkl=dij4*dkl - dik*djl - dil*djk
c
                      integ=integ+1
c------------------------------------------
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
cdonotdoit            xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
cdonotdoit            yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
cdonotdoit            zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
                      xyz(1,1)=xyz(1,1)-dijkl*xinta
                      xyz(1,2)=xyz(1,2)-dijkl*xintb
                      xyz(1,3)=xyz(1,3)-dijkl*xintc
cdonotdoit            xyz(1,4)=xyz(1,4)-dijkl*xintd
c
                      xyz(2,1)=xyz(2,1)-dijkl*yinta
                      xyz(2,2)=xyz(2,2)-dijkl*yintb
                      xyz(2,3)=xyz(2,3)-dijkl*yintc
cdonotdoit            xyz(2,4)=xyz(2,4)-dijkl*yintd
c
                      xyz(3,1)=xyz(3,1)-dijkl*zinta
                      xyz(3,2)=xyz(3,2)-dijkl*zintb
                      xyz(3,3)=xyz(3,3)-dijkl*zintc
cdonotdoit            xyz(3,4)=xyz(3,4)-dijkl*zintd
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
c translational invariance :
c
             xyz(1,4)=-xyz(1,1)-xyz(1,2)-xyz(1,3)
             xyz(2,4)=-xyz(2,1)-xyz(2,2)-xyz(2,3)
             xyz(3,4)=-xyz(3,1)-xyz(3,2)-xyz(3,3)
c
             do icr=1,3
                atforces(icr,iat)=atforces(icr,iat)+xyz(icr,1)
                atforces(icr,jat)=atforces(icr,jat)+xyz(icr,2)
                atforces(icr,kat)=atforces(icr,kat)+xyz(icr,3)
                atforces(icr,lat)=atforces(icr,lat)+xyz(icr,4)
             enddo
             if(nsym.gt.0) then
                call getival('SymNuPr',nupair)
                call atforce_symm(na,nsym,bl(nupair),iat,jat,kat,lat,
     *                            xyz,atforces)
             endif
c--------------------------------------------------------
  150   continue
  100 continue
      end
c=====================================================================
c  Forces on atoms for "PURE" DFT
c -------------------------------
c
      subroutine force_dft(buf,dens ,atforces,
     *                     lind,ntri,nbls,ngcd,lnijkl,
     *                     ncenter,labels,length,lgenct,
     *                     na,nsym )
c----------------------------------------------------------------
c INPUT :
c buf()    - integral derivatives buffer
c dens()   - density matrix
c lind()   - precomputed diagonals (  i*(i-1)/2 )
c ntri     - ncf*(ncf+1)/2
c nbls     - size of the integral derivative block (no of quartets)
c ngcd     - product of general contractions of 4 shells
c lnijkl   - number of integrals from one quartet
c ncenter()- mapping functions to atoms ( ncenter(icf)->iatom )
c
c OUTPUT :
c atforces() - forces at atoms
c----------------------------------------------------------------
      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dens(ntri)
      dimension atforces(3,*)             ! (3,natoms)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension xyz(3,4)
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
          call ZeroIT(xyz,12)
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
                      dijkl=dij4*dens(klf)
c
                      integ=integ+1
c------------------------------------------
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
cdonotdoit            xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
cdonotdoit            yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
cdonotdoit            zintd=-(zinta+zintb+zintc) ! trans. inv.
c--------------------------------------------------------
                      xyz(1,1)=xyz(1,1)-dijkl*xinta
                      xyz(1,2)=xyz(1,2)-dijkl*xintb
                      xyz(1,3)=xyz(1,3)-dijkl*xintc
cdonotdoit            xyz(1,4)=xyz(1,4)-dijkl*xintd
c
                      xyz(2,1)=xyz(2,1)-dijkl*yinta
                      xyz(2,2)=xyz(2,2)-dijkl*yintb
                      xyz(2,3)=xyz(2,3)-dijkl*yintc
cdonotdoit            xyz(2,4)=xyz(2,4)-dijkl*yintd
c
                      xyz(3,1)=xyz(3,1)-dijkl*zinta
                      xyz(3,2)=xyz(3,2)-dijkl*zintb
                      xyz(3,3)=xyz(3,3)-dijkl*zintc
cdonotdoit            xyz(3,4)=xyz(3,4)-dijkl*zintd
c--------------------------------------------------------
  350                 continue
  300              continue
  250           continue
  200        continue
c--------------------------------------------------------
c translational invariance :
c
             xyz(1,4)=-xyz(1,1)-xyz(1,2)-xyz(1,3)
             xyz(2,4)=-xyz(2,1)-xyz(2,2)-xyz(2,3)
             xyz(3,4)=-xyz(3,1)-xyz(3,2)-xyz(3,3)
c
             do icr=1,3
                atforces(icr,iat)=atforces(icr,iat)+xyz(icr,1)
                atforces(icr,jat)=atforces(icr,jat)+xyz(icr,2)
                atforces(icr,kat)=atforces(icr,kat)+xyz(icr,3)
                atforces(icr,lat)=atforces(icr,lat)+xyz(icr,4)
             enddo
             if(nsym.gt.0) then
                call getival('SymNuPr',nupair)
                call atforce_symm(na,nsym,bl(nupair),iat,jat,kat,lat,
     *                            xyz,atforces)
             endif
c--------------------------------------------------------
  150   continue
  100 continue
      end
c==============================================================
      subroutine make_negxyz(nsym,nsyo)
      common /negxyz/ negx(7),negy(7),negz(7)
      dimension nsyo(7)
      dimension mirror(3,7)
      data mirror/-1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c
      do 10 ns=1,nsym
         nop=nsyo(ns)
         negx(ns)=mirror(1,nop)
         negy(ns)=mirror(2,nop)
         negz(ns)=mirror(3,nop)
   10 continue
c
      end
c==============================================================
      subroutine atforce_symm(na,nsym,nupair,iat,jat,kat,lat,
     *                        xyz,atforces)
c
c--------------------------------------------------------
c symmetrize two-electron part of forces on atoms:
c
c Input :
c
c         atforces(3,natoms) -to be symmetrized
c Output: atforces(3,natoms) -symmetrized
c
c--------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common /negxyz/ negx(7),negy(7),negz(7)
c
      dimension nupair(na,nsym)           ! (natoms,nsym)
      dimension atforces(3,*)             ! (3,natoms)
      dimension xyz(3,4)
c
      do ns=1,nsym
         ngx=negx(ns)
         ngy=negy(ns)
         ngz=negz(ns)
c
         iat1=nupair(iat,ns)
         jat1=nupair(jat,ns)
         kat1=nupair(kat,ns)
         lat1=nupair(lat,ns)
c
         atforces(1,iat1)=atforces(1,iat1)+ngx*xyz(1,1)
         atforces(1,jat1)=atforces(1,jat1)+ngx*xyz(1,2)
         atforces(1,kat1)=atforces(1,kat1)+ngx*xyz(1,3)
         atforces(1,lat1)=atforces(1,lat1)+ngx*xyz(1,4)
c
         atforces(2,iat1)=atforces(2,iat1)+ngy*xyz(2,1)
         atforces(2,jat1)=atforces(2,jat1)+ngy*xyz(2,2)
         atforces(2,kat1)=atforces(2,kat1)+ngy*xyz(2,3)
         atforces(2,lat1)=atforces(2,lat1)+ngy*xyz(2,4)
c
         atforces(3,iat1)=atforces(3,iat1)+ngz*xyz(3,1)
         atforces(3,jat1)=atforces(3,jat1)+ngz*xyz(3,2)
         atforces(3,kat1)=atforces(3,kat1)+ngz*xyz(3,3)
         atforces(3,lat1)=atforces(3,lat1)+ngz*xyz(3,4)
      enddo
c
      end
c==============================================================
      subroutine addvec_abs(ntri,densA,densB,densC)
      implicit real*8 (a-h,o-z)
      dimension densa(ntri),densb(ntri),densc(ntri)
c
      do i=1,ntri
         densc(i)= abs( densa(i) ) + abs( densb(i) )
      enddo
c
      end
