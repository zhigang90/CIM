c  Forces on atoms for BASIC SCF
c ------------------------------
c
      subroutine force_uhf(buf,denA,denB,atforces,
     *                     lind,ntri,nbls,ngcd,lnijkl,
     *                     ncenter,labels,length,lgenct,
     *                     na,nsym )
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
c OUTPUT :
c  atforces() - forces on atoms
c----------------------------------------------------------------

      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension atforces(3,*)             ! (3,natoms)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension xyz(3,4)
      Parameter (Two=2.0d0)
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
                dij2=(denA(ijf)+denB(ijf))*Two
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
                      dijkl=Two*(dij2*dkl - dikA*djlA - dilA*djkA
     $                                    - dikB*djlB - dilB*djkB)
c
                      integ=integ+1
c-------------------------------------------------------
                      xinta=buf(1,ijklp,integ,iqu)
                      xintb=buf(2,ijklp,integ,iqu)
                      xintc=buf(3,ijklp,integ,iqu)
c                     xintd=-(xinta+xintb+xintc) ! trans. inv.
                      yinta=buf(4,ijklp,integ,iqu)
                      yintb=buf(5,ijklp,integ,iqu)
                      yintc=buf(6,ijklp,integ,iqu)
c                     yintd=-(yinta+yintb+yintc) ! trans. inv.
                      zinta=buf(7,ijklp,integ,iqu)
                      zintb=buf(8,ijklp,integ,iqu)
                      zintc=buf(9,ijklp,integ,iqu)
c                     zintd=-(zinta+zintb+zintc) ! trans. inv.
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
ctest
c     write(6,*)' atoms &  forces'
c     do iat=1,na
c       fx3=1.d-10*atforces(1,iat)
c       fy3=1.d-10*atforces(2,iat)
c       fz3=1.d-10*atforces(3,iat)
c       write(6,222) iat, fx3,fy3,fz3
c     enddo
 222  format(1x,i3,2x,3(f9.6,1x) )
ctest
      end
c=====================================================================
c  Forces on atoms for ACM DFT
c ----------------------------
c
      subroutine force_uacm(ax,buf,denA,denB,atforces,
     *                      lind,ntri,nbls,ngcd,lnijkl,
     *                      ncenter,labels,length,lgenct,
     *                      na,nsym )
c----------------------------------------------------------------
c INPUT :
c  ax       - factor to multiply exchange term
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
c  atforces() - forces on atoms
c----------------------------------------------------------------
      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension atforces(3,*)             ! (3,natoms)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension xyz(3,4)
      Parameter (Two=2.0d0)
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
                dij2=(denA(ijf)+denB(ijf))*Two
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
                      dijkl=Two*(dij2*dkl - ax*(dikA*djlA + dilA*djkA
     $                                         +dikB*djlB + dilB*djkB))
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
      subroutine force_udft(buf,denA,denB,atforces,
     *                      lind,ntri,nbls,ngcd,lnijkl,
     *                      ncenter,labels,length,lgenct,
     *                      na,nsym )
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
c OUTPUT :
c  atforces() - forces on atoms
c----------------------------------------------------------------
      use memory

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
c
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension denA(ntri),denB(ntri)
      dimension atforces(3,*)             ! (3,natoms)
      dimension lind(*) , ncenter(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c local :
      dimension xyz(3,4)
      Parameter (Four=4.0d0)
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
                dij4=(denA(ijf)+denB(ijf))*Four
                do 300 kkk=1,klen
                   kcf=kcff+kkk
                   kk=lind(kcf)
                      do 350 lll=1,llen
                      lcf=lcff+lll
                      ll=lind(lcf)
                      klf=kk+lcf
                      if(lcf.gt.kcf) klf=ll+kcf
c
                      dijkl=dij4*(denA(klf)+denB(klf))
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
