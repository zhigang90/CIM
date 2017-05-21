C APR 1 98 PP  I have removed the thresholds on the overlap etc.
c  integrals. They are insignificant compared to the nuclear
c  attraction force. I ngelect the nucl. attraction force by
c  a different algorithm: S0*abs(zzz)*max(sqa,sqb), where
c  zza is the nucl. charge and sqa,sqb are square roots of the
c exponents. Note that one should use a nad b but this is
c consistent with the energy expression.
c===================================================================
c             One-electron part of the hessian
c
c                    July 2000, K.Wolinski
c
c===================================================================
      subroutine hess1(rhf,bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      Logical rhf
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c     common /intbl/ifpp,inxx(100)
c external field info :
      common /fieldx/ xfield,xfgrad,elfiel(9)
c Hessian  options :
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------------------------
c local
      dimension natodo(2,1000)
c------------------------------------------------------------------
      call matreset
      call matmark
c------------------------------------------------------------------
      call secund(tgrad1)
      call elapsec(elagrad1)
c------------------------------------------------------------------
      ifield=xfield
      ifgrad=xfgrad
c------------------------------------------------------------------
c memory was allocated in hessana : get addresses
c
      call getival('lhess', lhess)
      call getival('ldensi',lden  )
      call getival('ldewsi',ldew  )
      if(.not.rhf)then
        call getival('ldensB',mden  )
      endif
      call getival('npsp',npsp)            ! pseudopotentials
c
c  lhess contains already 2-el. contr. Tr{D0*G(D0,g2)} (rescaled)
c lfock1 contains already 2-el. contr.  F(D0,G(D0,g1)} (rescaled)
c
c------------------------------------------------------------------
c check-run only or lack of scf convergance :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c------------------------------------------------------------------
c Calculate nuclear repulsion part of the hessian
c
      call nucrep_hess(na,bl(inuc),bl(lhess))
c
      if(nprint.ge.2) then
         call hess_print(bl(lhess),na,1.d0,
     *                   '+ contr. no=2 : nuclear repulsion')
      endif
c------------------------------------------------------------------
c---------- fisrt -order derivatives ------------------------------
c------------------------------------------------------------------
c find out for how many atoms Fock-derivative matrices (H1)
c can be calculated at once : -5 tells that it is called from hess1
c----------------------------------------------------------------------
c H1 are constracted for symmetry uniqe atoms only
c
      call getival('natunq',natunq)      ! number of unique atoms
      call getival('listunq',listunq)
      call getival('listrea',listreal)
c
      call fder4nat(-5,rhf,natunq,ntri,ntimes,natonce)
c
c        ntimes=2
c        natonce=3
c---------------------------------------------------------------------
      call getival('iout',iout)
         if(ntimes.eq.1) write(iout,100)
  100    format(/'Hder will be calculated in 1 pass')
         if(ntimes.gt.1) write(iout,101) ntimes
  101    format(/'Hder will be calculated in ',i2,' passes')
         write(iout,*) '  '
      call f_lush(iout)
c---------------------------------------------------------------------
c
      call make_natodo(natunq,natonce,ntimes,natodo)
c
c     natodo(1,itime)-> nat_begining
c     natodo(2,itime)-> nat_end
c----------------------------------------------------------------------
      DO ITIME=1,NTIMES
         natb=natodo(1,itime)
         nate=natodo(2,itime)
         natonce=nate-natb+1
         write(iout,1234) itime,natonce
 1234    format('  Constructing Hder matrices ',I3,' time for ',I3,
     $          ' atoms')
         write(iout,*) '  '
         call f_lush(iout)
c------------------------------------------------------------------
c read Fder matrices from a disk and transpose them to (3,nat,ntri)
c
      call getmem(ntri*3*natonce,lfock1)
c
      call getmem(ntri*3,itmp)
      call trsp_fromd(61,ntri,natb,nate,bl(lfock1),bl(itmp),bl(listunq))
      call retmem(1)
c
      if(.not.rhf) then
          call getmem(ntri*3*natonce,lfock1b)
          call getmem(ntri*3,itmp)
          call trsp_fromd(71,ntri,natb,nate,bl(lfock1b),bl(itmp),
     *                    bl(listunq))
          call retmem(1)
      endif
c------------------------------------------------------------------
c Calculate T1 (kinetic) first derivative integrals
c for derivative fock matrix
c
      if(rhf)then
        call intonFder(natb,nate,inx,ifield,ifgrad,elfiel,
     *              bl(ibas),bl(inuc),ncs,ncf,ntri,
     *              bl(listreal),bl(lfock1) )
      else
        call intonFderu(natb,nate,inx,ifield,ifgrad,elfiel,
     *              bl(ibas),bl(inuc),ncs,ncf,ntri,
     *              bl(listreal),
     *              bl(lfock1),bl(lfock1B) )
      endif
c
      if(nprint.ge.3) then
         if(rhf)then
           write(6,*)' Derivative Fock after 1st-order Kinetic '
           call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
         else
           write(6,*)' Alpha derivative Fock after 1st-order Kinetic '
           call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
           write(6,*)' Beta derivative Fock after 1st-order Kinetic '
           call fder_print(bl(lfock1B),nate-natb+1,ntri,1.d0  )
         endif
      endif
c------------------------------------------------------------------
c Calculate first-order electron-nuclear attraction contributions
c to the derivative Fock matrix
c
      if(rhf)then
        call elenuc_Fder(na,natb,nate,inx,bl(ibas),bl(inuc),ncs,ncf,
     *                   ntri,bl(listreal),bl(lfock1) )
      else
        call elenuc_Fdeu(na,natb,nate,inx,bl(ibas),bl(inuc),ncs,ncf,
     *                    ntri,bl(listreal),bl(lfock1),bl(lfock1B))
      endif
c
      if(nprint.ge.3) then
         if(rhf)then
           write(6,*)' Derivative Fock after 1st-order Ele-Nuc  '
           call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
         else
           write(6,*)' Alpha derivative Fock after 1st-order Ele-Nuc  '
           call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
           write(6,*)' Beta derivative Fock after 1st-order Ele-Nuc  '
           call fder_print(bl(lfock1B),nate-natb+1,ntri,1.d0  )
         endif
      endif
c------------------------------------------------------------------
c Calculate first-order pseudopotential contributions
c to the derivative Fock matrix
c
      if(npsp.ne.0)then
cc        call secund(pspt0)
cc        call elapsec(pspet0)
        if(rhf)lfock1B=lfock1
        call pspFder(rhf,na,natb,nate,npsp,ncf,ncs,inx,bl(inuc),
     $               bl(ibas),bl(listreal),bl(lfock1),bl(lfock1B))
        if(nprint.ge.3) then
          if(rhf)then
            write(6,*)' Derivative Fock after 1st-order Pseudopot'
            call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
          else
            write(6,*)' Alpha derivative Fock after 1st-order Pseudopot'
            call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
            write(6,*)' Beta derivative Fock after 1st-order Pseudopot'
            call fder_print(bl(lfock1B),nate-natb+1,ntri,1.d0  )
          endif
        endif
cc        call secund(pspt1)
cc        call elapsec(pspet1)
cc        psptime=pspt1-pspt0
cc        pspelap=pspet1-pspet0
      endif
c------------------------------------------------------------------
c -- First-order derivative Fock matrices are done.
c -- Need to store them on disk in transposed form on unit 61
c -- transpose fock1(3,natoms,ntri) into fock1t(ntri,3,natoms)
c -- For UHF, beta Fock1 is on unit 71
c
c transpose Fder matrices into (ntri,3,nat) and write them on a disk
c
      call getmem(ntri*3,itmp)
ckw   call trsp_ondsk(61,ntri,natb,nate,bl(lfock1),bl(itmp))
      call trsp_ondsx(61,ntri,natb,nate,bl(lfock1),bl(itmp),bl(listunq))
      call retmem(1)
c
      if(.not.rhf) then
          call getmem(ntri*3,itmp)
ckw       call trsp_ondsk(71,ntri,natb,nate,bl(lfock1b),bl(itmp))
          call trsp_ondsx(71,ntri,natb,nate,bl(lfock1b),bl(itmp),
     *                    bl(listunq))
          call retmem(1)
      endif
c------------------------------------------------------------------
c Calculate S1 (overlap) first derivative integrals
c Use lfock1 allocation :
c
      call zeroit(bl(lfock1),natonce*3*ntri)
      call intonSder(natb,nate,inx,bl(ibas),bl(inuc),ncs,ncf,ntri,
     *               bl(listreal),bl(lfock1) )
c
      if(nprint.ge.3) then
         write(6,*)' 1st-order Overlap          '
         call fder_print(bl(lfock1),nate-natb+1,ntri,1.d0  )
      endif
c
c transpose Sder matrices into (ntri,3,nat) and write them on a disk
c
      call getmem(ntri*3,itmp)
ccc   call trsp_ondsk(62,ntri,natb,nate,bl(lfock1),bl(itmp))
      call trsp_ondsx(62,ntri,natb,nate,bl(lfock1),bl(itmp),bl(listunq))
      call retmem(1)
c------------------------------------------------------------------
      if(.not.rhf) call retmem(1)  !   allocated for lfock1b
      call retmem(1)               !   allocated for lfock1
      ENDDO     !        ITIME=1,NTIMES
c------------------------------------------------------------------
c At this point we have F1 and S1 on a disk for symmetry unique atoms
c We need to generate these for all atoms.
c
c
      call getival('nsym',nsym)
      If( nsym.ne.0) Then
        call getival('SymFunPr',ifp)
        call getival('SymNuPr1',nupair)
        call getival('nsyo',nsyo)
c
        call getmem(6*ntri,lfock1)   ! for f1(ntri,3,2) 2 atoms
        iunit=62     !  overlap1
        call fdersymm_hes2(iunit,0, nsym  ,ncf,ntri,na,bl(lfock1),
     *                    bl(nsyo),bl(ifp),bl(nupair),
     *                    natunq,bl(listunq))
        iunit=61     !  Fcock1
        call fdersymm_hes2(iunit,0, nsym  ,ncf,ntri,na,bl(lfock1),
     *                    bl(nsyo),bl(ifp),bl(nupair),
     *                    natunq,bl(listunq))
        if(.not.rhf) then
          iunit=71
          call fdersymm_hes2(iunit,0, nsym  ,ncf,ntri,na,bl(lfock1),
     *                      bl(nsyo),bl(ifp),bl(nupair),
     *                      natunq,bl(listunq))
        endif
c
        call retmem(1)
      EndIf
c------------------------------------------------------------------
c---------- second-order derivatives ------------------------------
c------------------------------------------------------------------
c Calculate S2 & T2 (overlap & kinetic) second derivative integrals
c and put them to the final hessian matrix
c
      if(rhf)then
        call intonHall(na, inx, ifield, ifgrad, elfiel,
     *              bl(ibas),bl(inuc),ncs, ncf, ntri,
     *              bl(lden), bl(ldew), bl(lhess) )
      else
c
c  compute total (alpha + beta) density and store in ldenab
c
        call getmem(ntri,ldenab)
        call zeroit(bl(ldenab),ntri)
        Call AddVEC(ntri,bl(lden),bl(mden),bl(ldenab))
c
c  now compute contribution to hessian using total density
c
        call intonHall(na, inx, ifield, ifgrad, elfiel,
     *              bl(ibas),bl(inuc),ncs, ncf, ntri,
     *              bl(ldenab), bl(ldew), bl(lhess) )
      endif
c
      if(nprint.ge.2) then
         call hess_print(bl(lhess),na,1.d0,
     *                   '+ contr. no=3 : S2*W0 + T2*D0    ')
      endif
c--------------------------------------------------------------------
c Calculate second-order electron-nuclear attraction contributions
c to the final hessian
c
      if(rhf)then
        call elenuc_Hess(na,inx,bl(ibas),bl(inuc),ncs, ncf,ntri,
     *                   bl(lden),bl(lhess) )
      else
c
c  for the uhf case, we use the total (alpha + beta) density,
c  that has been computed and stored in ldenab (see above)
c
        call elenuc_Hess(na,inx,bl(ibas),bl(inuc),ncs, ncf,ntri,
     *                   bl(ldenab),bl(lhess) )
      endif
c
      if(nprint.ge.2) then
         call hess_print(bl(lhess),na,1.d0,
     *                   '+ contr. no=4 :         V2*D0    ')
      endif
c--------------------------------------------------------------------
c Calculate second-order pseudopotential contributions
c to the final hessian
c
      if(npsp.ne.0)then
cc        call secund(pspt0)
cc        call elapsec(pspet0)
        if(rhf)then
          call psphess(na,npsp,ncf,ncs,inx,bl(inuc),bl(ibas),
     $         bl(lden),bl(lhess))
        else
          call psphess(na,npsp,ncf,ncs,inx,bl(inuc),bl(ibas),
     $         bl(ldenab),bl(lhess))
        endif
        if(nprint.ge.2) then
           call hess_print(bl(lhess),na,1.d0,
     *                     '+ contr. no=5 : Pseudopotentials ')
        endif
cc        call secund(pspt1)
cc        call elapsec(pspet1)
cc        psptime=pspt1-pspt0
cc        pspelap=pspet1-pspet0
      endif
c--------------------------------------------------------------------
*
* PIB: Calculate the D0 contribution to the dipole moment derivatives
*
      Call GetIVal('na',nAtoms)
*---- Allocate memory for the D0*Mx contribution to the dip. mom. der.
      Call getmem(3*3*nAtoms,iDMDer)
*
      if(rhf)then
        Call IntDDer(nAtoms,inx,bl(ibas),bl(inuc),nCS,nCF,nTri,
     &               bl(lDen),bl(iDMDer))
      else
c
c  for the uhf case, we use the total (alpha + beta) density,
c  that has been computed and stored in ldenab (see above)
c
        Call IntDDer(nAtoms,inx,bl(ibas),bl(inuc),nCS,nCF,nTri,
     &               bl(ldenab),bl(iDMDer))
c
        call retmem(1)    ! allocated for ldenab
      endif
*
c------------------------------------------------------------------
      call retmem(1)   ! alloc for iDMDer
c------------------------------------------------------------------
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) total/60.0d0,elaps/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for 1e part of hessian  ='
     *,f8.2,' Elapsed = ',f8.2,' min')
c------------------------------------------------------------------
      call matremark
c------------------------------------------------------------------
      end
c======================================================================
      subroutine intonSder(natb,nate,inx,datbas,datnuc,ncs,ncf,ntri,
     *                     listreal,over1)
c---------------------------------------------------------------------
c This routine calculates Ist derivatives of the following integrals :
c             overlap ONLY
c it is used for analytical hessian
c---------------------------------------------------------------------
c INPUT :
c natb,nate       - calc Sder contributions for atoms <natb,nate>
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c
c over1(3,natb:nate,ntri) - deriv. overlatp matrix
c
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat, dojat
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension listreal(*)
      dimension over1(3,natb:nate,ntri)
clocal :
      dimension s3(3*784) ! 3*(28x28) (i|i)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
c
c check if iat is a symmetry unique atom
c if iat is not a unique atom then iau=0
c
         iau=listreal(iat)
         doiat=.false.
         if(iau.ge.natb .and. iau.le.nate) doiat=.true.
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               jau=listreal(jat)
               dojat=.false.
               if(jau.ge.natb .and. jau.le.nate) dojat=.true.
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
                  if(iat.ne.jat) then
                     do l=1,len3
                        s3(l)=zero
                     enddo
c
c................... overlap integral derivatives
c
                     key=1
                     call onef(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
c
                  endif             ! if(iat.ne.jat) then
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
               IF(IAT.NE.JAT) THEN
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
c.................................................................
                 if(doiat) then
                  over1(icr,iau,ij)=over1(icr,iau,ij) + s3(iij_icr)
                 endif
                 if(dojat) then
                  over1(icr,jau,ij)=over1(icr,jau,ij) - s3(iij_icr)
                 endif
c.................................................................
               enddo
               ENDIF
c
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine intonFder(natb,nate,inx,ifield,ifgrad,xfld,
     *                     datbas,datnuc,ncs,ncf,ntri,
     *                     listreal,fock1)
c---------------------------------------------------------------------
c This routine calculates Ist derivatives of the following integrals :
c
c              kinetic , dipole & quadrupole
c
c last two for an external electric field & field gradient
c
c it is used for analytical hessian
c---------------------------------------------------------------------
c INPUT :
c natb,nate       - calc Fder contributions for atoms <natb,nate>
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c
c fock1(3,natb:nate,ntri) - deriv. fock : kinetc Ist-order contr.
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat, dojat
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension listreal(*)
      dimension fock1(3,natb:nate,ntri)
clocal :
      dimension t3(3*784),f3(3*784) ! 3*(28x28) (i|i)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
c
c check if iat is a symmetry unique atom
c if iat is not a unique atom then iau=0
c
         iau=listreal(iat)
         doiat=.false.
         if(iau.ge.natb .and. iau.le.nate) doiat=.true.
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               jau=listreal(jat)
               dojat=.false.
               if(jau.ge.natb .and. jau.le.nate) dojat=.true.
c
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
                  if(iat.ne.jat) then
                     do l=1,len3
                        t3(l)=zero
                     enddo
c
c................... kinetic energy integral derivatives
c
                     key=2
                     call onef(i,j,igc,jgc,datbas,t3,inx,key,kk,k2)
c
c................... an external electric field & field gradient.....
c
                     key=1
                     if(ifield.eq.1) then
c ......................an external electric field
c                       (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)
c
                        do k1=1,3
                           field=xfld(k1)
                           if(field.ne.zero) then
                              do l=1,len3
                                 f3(l)=zero
                              enddo
c
                              call onef(i,j,igc,jgc,datbas,f3,
     *                                  inx,key,k1,k2)
c
                              do l=1,len3
                                 t3(l)=t3(l)+f3(l)*field
                              enddo
                           endif    ! field.ne.zero
                        enddo       ! k1=1,3
                     endif          ! if(ifield.eq.1) then
c
                     if(ifgrad.eq.1) then
c ......................an external electric field gradient
c                       (i10 /xy/ j00) + (i00 /xy/ j10)
c                       where xy=xx,yx,yy,zx,zy,zz
c
                        do k1=4,9
                           field=xfld(k1)
                           if(field.ne.zero) then
                              do l=1,len3
                                 f3(l)=zero
                              enddo
c
                              call onef(i,j,igc,jgc,datbas,f3,
     *                                  inx,key,k1,k2)
c
                              do l=1,len3
                                 t3(l)=t3(l)+f3(l)*field
                              enddo
                           endif    ! field.ne.zero
                        enddo       ! k1=4,9
                     endif          ! if(ifield.eq.1) then
                  endif             ! if(iat.ne.jat) then
c
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
               IF(IAT.NE.JAT) THEN
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
c.................................................................
                  if(doiat) then
                     fock1(icr,iau,ij)=fock1(icr,iau,ij) + t3(iij_icr)
                  endif
                  if(dojat) then
                     fock1(icr,jau,ij)=fock1(icr,jau,ij) - t3(iij_icr)
                  endif
c.................................................................
               enddo
               ENDIF
c
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine intonFderu(natb,nate,inx,ifield,ifgrad,xfld,
     *                  datbas,datnuc,ncs,ncf,ntri,
     *                  listreal,fock1A,fock1B)
c---------------------------------------------------------------------
c This routine calculates Ist derivatives of the following integrals :
c
c              kinetic , dipole & quadrupole
c
c last two for an external electric field & field gradient
c
c it is used for analytical hessian
c---------------------------------------------------------------------
c INPUT :
c natb,nate       - calc Fder contributions for atoms <natb,nate>
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c
c fock1A(3,natb:nate,ntri) - Alpha deriv. fock : kinetc Ist-order contr.
c fock1B(3,natb:nate,ntri) - Beta deriv. fock : kinetc Ist-order contr.
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat, dojat
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension listreal(*)
      dimension fock1A(3,natb:nate,ntri),fock1B(3,natb:nate,ntri)
clocal :
      dimension t3(3*784),f3(3*784) ! 3*(28x28) (i|i)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
c
c check if iat is a symmetry unique atom
c if iat is not a unique atom then iau=0
c
         iau=listreal(iat)
         doiat=.false.
         if(iau.ge.natb .and. iau.le.nate) doiat=.true.
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               jau=listreal(jat)
               dojat=.false.
               if(jau.ge.natb .and. jau.le.nate) dojat=.true.
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
                  if(iat.ne.jat) then
                     do l=1,len3
                        t3(l)=zero
                     enddo
c
c................... kinetic energy integral derivatives
c
                     key=2
                     call onef(i,j,igc,jgc,datbas,t3,inx,key,kk,k2)
c
c................... an external electric field & field gradient.....
c
                     key=1
                     if(ifield.eq.1) then
c ......................an external electric field
c                       (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)
c
                        do k1=1,3
                           field=xfld(k1)
                           if(field.ne.zero) then
                              do l=1,len3
                                 f3(l)=zero
                              enddo
c
                              call onef(i,j,igc,jgc,datbas,f3,
     *                                  inx,key,k1,k2)
c
                              do l=1,len3
                                 t3(l)=t3(l)+f3(l)*field
                              enddo
                           endif    ! field.ne.zero
                        enddo       ! k1=1,3
                     endif          ! if(ifield.eq.1) then
c
                     if(ifgrad.eq.1) then
c ......................an external electric field gradient
c                       (i10 /xy/ j00) + (i00 /xy/ j10)
c                       where xy=xx,yx,yy,zx,zy,zz
c
                        do k1=4,9
                           field=xfld(k1)
                           if(field.ne.zero) then
                              do l=1,len3
                                 f3(l)=zero
                              enddo
c
                              call onef(i,j,igc,jgc,datbas,f3,
     *                                  inx,key,k1,k2)
c
                              do l=1,len3
                                 t3(l)=t3(l)+f3(l)*field
                              enddo
                           endif    ! field.ne.zero
                        enddo       ! k1=4,9
                     endif          ! if(ifield.eq.1) then
                  endif             ! if(iat.ne.jat) then
c
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
               IF(IAT.NE.JAT) THEN
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
c.................................................................
                if(doiat) then
                  fock1A(icr,iau,ij)=fock1A(icr,iau,ij) + t3(iij_icr)
                  fock1B(icr,iau,ij)=fock1B(icr,iau,ij) + t3(iij_icr)
                endif
                if(dojat) then
                  fock1A(icr,jau,ij)=fock1A(icr,jau,ij) - t3(iij_icr)
                  fock1B(icr,jau,ij)=fock1B(icr,jau,ij) - t3(iij_icr)
                endif
c.................................................................
               enddo
               ENDIF
c
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine elenuc_Fder(natoms,natb,nate,inx,datbas,datnuc,ncs,
     *                       ncf, ntri,listreal,fock1 )
      implicit real*8 (a-h,o-z)
      logical donra
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension listreal(*)
      dimension fock1(3,natb:nate,ntri) ! input/output derivative fock
c---------------------------------------------------------------------
c include charged dummies :
c
      call getival('ndum1',ndum1)
c---------------------------------------------------------------------
ccc   do nra=1,natoms
      do nra=1,natoms+ndum1
c
c        check if nra is a symmetry unique atom
c        if nra is not a unique atom then nru=0
c
        if(nra.le.natoms) then
           nru=listreal(nra)
           donra=.false.
           if(nru.ge.natb .and. nru.le.nate) donra=.true.
        else
           nru=0  ! irrelevant, not used in this case
           donra=.false.  !   if charged dummy
        endif
c
         call intofder7(donra,nra,nru,natb,nate,
     *                  inx,datbas,datnuc, ncs,ncf,
     *                  ntri, listreal,fock1)
      enddo
c
      end
c======================================================================
      subroutine intofder7(donra,nra,nru, natb,nate,
     *                     inx,datbas,datnuc,ncs,ncf,
     *                     ntri,listreal,fock1)
c---------------------------------------------------------------------
c
c This routine calculates Ist derivatives of the following integrals :
c
c             nuclear atraction & Hellmann-Feynman forces
c
c---------------------------------------------------------------------
c INPUT :
c nra             - reference atom
c nru   :         - symmetry unique image of nra
c natb,nate       - calc Fder contributions for atoms <natb,nate>
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c fock1(3,natoms,ntri) - fock derivatives - INPUT/OUTPUT
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat, dojat,donra
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension listreal(*)
      dimension fock1(3,natb:nate,ntri)
clocal :
      dimension sa(3*784),sb(3*784) ! 3*(28x28) (i|i)
      dimension xra(3)
      data zero /0.d0/
c---------------------------------------------------------------------
      zza=datnuc(1,nra)
c
      xra(1)=datnuc(2,nra)
      xra(2)=datnuc(3,nra)
      xra(3)=datnuc(4,nra)
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
c
c check if iat is a symmetry unique atom
c if iat is not a unique atom then iau=0
c
         iau=listreal(iat)
         doiat=.false.
         if(iau.ge.natb .and. iau.le.nate) doiat=.true.
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               jau=listreal(jat)
               dojat=.false.
               if(jau.ge.natb .and. jau.le.nate) dojat=.true.
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
c
c................... nuclear attraction derivatives
c
                     call onef7(i,j,igc,jgc,datbas,sa,sb,inx,xra,zza)
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
c...............................................
                if(donra) then
                  fock1(icr,nru,ij)=fock1(icr,nru,ij)
     *                                 +(sa(iij_icr)+sb(iij_icr))*zza
                endif
                if(doiat) then
                  fock1(icr,iau,ij)=fock1(icr,iau,ij)-sa(iij_icr)*zza
                endif
                if(dojat) then
                  fock1(icr,jau,ij)=fock1(icr,jau,ij)-sb(iij_icr)*zza
                endif
c...............................................
               enddo
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine elenuc_Fdeu(natoms,natb,nate,inx,datbas,datnuc,ncs,ncf,
     *                        ntri,listreal,fock1A,fock1B)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
      logical donra
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension listreal(*)
      dimension fock1A(3,natb:nate,ntri),fock1B(3,natb:nate,ntri)
c---------------------------------------------------------------------
c include charged dummies :
c
      call getival('ndum1',ndum1)
c---------------------------------------------------------------------
ccc   do nra=1,natoms
      do nra=1,natoms+ndum1
c
c        check if nra is a symmetry unique atom
c        if nra is not a unique atom then nru=0
c
        if(nra.le.natoms) then
           nru=listreal(nra)
           donra=.false.
           if(nru.ge.natb .and. nru.le.nate) donra=.true.
        else
           nru=0  ! irrelevant, not used in this case
           donra=.false.  !   if charged dummy
        endif
c
         call intofde7u(donra,nra,nru,natb,nate,
     *                  inx,datbas,datnuc,ncs,ncf,
     *                  ntri,listreal,fock1A,fock1B)
      enddo
c
      end
c======================================================================
      subroutine intofde7u(donra,nra,nru, natb,nate,
     *                     inx,datbas,datnuc,ncs,ncf,
     *                     ntri,listreal,fock1A,fock1B)
c---------------------------------------------------------------------
c
c This routine calculates Ist derivatives of the following integrals :
c
c             nuclear atraction & Hellmann-Feynman forces
c
c---------------------------------------------------------------------
c INPUT :
c nra             - reference atom
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c fock1(3,natoms,ntri) - fock derivatives - INPUT/OUTPUT
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doiat, dojat,donra
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension listreal(*)
      dimension fock1A(3,natb:nate,ntri),fock1B(3,natb:nate,ntri)
clocal :
      dimension sa(3*784),sb(3*784) ! 3*(28x28) (i|i)
      dimension xra(3)
      data zero /0.d0/
c---------------------------------------------------------------------
      zza=datnuc(1,nra)
c
      xra(1)=datnuc(2,nra)
      xra(2)=datnuc(3,nra)
      xra(3)=datnuc(4,nra)
c
c check if nra is a symmetry unique atom
c if nra is not a unique atom then nru=0
c
c     nru=listreal(nra)
c     donra=.false.
c     if(nru.ge.natb .and. nru.le.nate) donra=.true.
c
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         iau=listreal(iat)
         doiat=.false.
         if(iau.ge.natb .and. iau.le.nate) doiat=.true.
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               jau=listreal(jat)
               dojat=.false.
               if(jau.ge.natb .and. jau.le.nate) dojat=.true.
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
c no need to zero out - zeroed out in onef7
c                  do l=1,len3
c                     sa(l)=zero
c                     sb(l)=zero
c                  enddo
c
c................... nuclear attraction derivatives
c
                     call onef7(i,j,igc,jgc,datbas,sa,sb,inx,xra,zza)
c
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
c...............................................
                if(donra) then
                  fock1A(icr,nru,ij)=fock1A(icr,nru,ij)
     *                                 +(sa(iij_icr)+sb(iij_icr))*zza
                  fock1B(icr,nru,ij)=fock1B(icr,nru,ij)
     *                                 +(sa(iij_icr)+sb(iij_icr))*zza
                endif
                if(doiat) then
                  fock1A(icr,iau,ij)=fock1A(icr,iau,ij)-sa(iij_icr)*zza
                  fock1B(icr,iau,ij)=fock1B(icr,iau,ij)-sa(iij_icr)*zza
                endif
                if(dojat) then
                  fock1A(icr,jau,ij)=fock1A(icr,jau,ij)-sb(iij_icr)*zza
                  fock1B(icr,jau,ij)=fock1B(icr,jau,ij)-sb(iij_icr)*zza
                endif
c...............................................
               enddo
   40             continue
                  jfu=jfu+len2
c
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine nucrep_hess(natoms,datnuc,hessian)
      implicit none
      integer natoms
      real*8 datnuc(5,natoms), hessian(3,natoms,3,natoms)
c---------------------------------------------------------------------
c Purpose:
c This routine calculates nuclear repulsion terms of the hessian
c
c upper triangle is calculated
c---------------------------------------------------------------------
c INPUT :
c natoms - number of atoms
c datnuc(5,natoms) - nuclear data
c
c OUTPUT:
c hessian(3,natoms,3,natoms) - components of the second derivative
c                            - of the nuclear potential are stored
c                            - directly in the hessian storage area
c                            - [Upper Triangle ONLY]
c---------------------------------------------------------------------
c
c     ..Parameters..
      real*8     one
      parameter (one=1.0D0)
      real*8     three
      parameter (three=3.0D0)
c     ..Local Scalars..
      logical doA,doB,doAB
      integer ndum1
      integer nra,nrb
      real*8  zza,zzb,xra,yra,zra,xrb,yrb,zrb,xab,yab,zab,
     &        r,rm32,rm52,zazbrm32,zazb3rm52
      real*8  val1,val2,val3,val4,val5,val6
c
c     .. Executable Statements ..
c---------------------------------------------------------------------
c check for charged dummy atoms. If they are present then include them
c
      call getival('ndum1',ndum1)
c
ccc   write(6,*)' no of real atoms @ charged dummies = ',natoms,ndum1
c---------------------------------------------------------------------
c.....loop over a < b=2,natoms ... the nuclear potential terms
c
c **** I assume that the hessian has been zeroed or this is not the
c **** first thing in there
c
cccc  do nrb=2,natoms
      do nrb=2,natoms+ndum1
         !get coords and charge for atom a
         zzb=datnuc(1,nrb)
         xrb=datnuc(2,nrb)
         yrb=datnuc(3,nrb)
         zrb=datnuc(4,nrb)
         doB=.true.
         if(nrb.gt.natoms) doB=.false.
         do nra=1,nrb-1
            ! get coords and charge for atom b
            zza=datnuc(1,nra)
            xra=datnuc(2,nra)
            yra=datnuc(3,nra)
            zra=datnuc(4,nra)
            doA=.true.
            if(nrA.gt.natoms) doA=.false.
c
            doAB=.true.
            if(.not.doA .or. .not.doB) doAB=.false.
c
            ! rab^2=(xa-xb)^2 + (ya-yb)^2 + (za-zb)^2
            xab=xra-xrb
            yab=yra-yrb
            zab=zra-zrb
            r=xab*xab+yab*yab+zab*zab

            ! rab^-3/2
            rm32=r*r*r
            rm32=sqrt(rm32)
            rm32=one/rm32

            ! rab^-5/2
            rm52=rm32/r

            ! Za*Zb/r^3/2 and 3*Za*Zb/r^5/2
            zazbrm32=zza*zzb*rm32
            zazb3rm52=zza*zzb*three*rm52

            !
            ! load the a-a, b-b, a-b blocks of the hessian (21 terms)
            ! but only 6 distinct formulas
            !
c
            val1= zazb3rm52 * xab * xab - zazbrm32
            val2= zazb3rm52 * xab * yab
            val3= zazb3rm52 * xab * zab
            val4= zazb3rm52 * yab * yab - zazbrm32
            val5= zazb3rm52 * yab * zab
            val6= zazb3rm52 * zab * zab - zazbrm32
c
ctest
c           if(.not.doB .and. doA) then
c             write(6,*) ' contributions to A,A A=',nra
c             write(6,66)nra,nrb,val1,val2,val3,val4,val5,val6
c66           format('Atoms=',2(i2,1x),' val=',6(f10.6,1x))
c           endif
ctest
c
            ! xaxa=xbxb=xaxb
            if(doA) hessian(1,nra,1,nra)=hessian(1,nra,1,nra) + val1
            if(doB) hessian(1,nrb,1,nrb)=hessian(1,nrb,1,nrb) + val1
            if(doAB)hessian(1,nra,1,nrb)=hessian(1,nra,1,nrb) - val1
c
            ! xaya=xbyb=-xayb=-yaxb
            if(doA) hessian(1,nra,2,nra)=hessian(1,nra,2,nra) + val2
            if(doB) hessian(1,nrb,2,nrb)=hessian(1,nrb,2,nrb) + val2
            if(doAB)hessian(2,nra,1,nrb)=hessian(2,nra,1,nrb) - val2
            if(doAB)hessian(1,nra,2,nrb)=hessian(1,nra,2,nrb) - val2

            ! xaza=xbzb=-xazb=-zaxb
            if(doA) hessian(1,nra,3,nra)=hessian(1,nra,3,nra) + val3
            if(doB) hessian(1,nrb,3,nrb)=hessian(1,nrb,3,nrb) + val3
            if(doAB)hessian(3,nra,1,nrb)=hessian(3,nra,1,nrb) - val3
            if(doAB)hessian(1,nra,3,nrb)=hessian(1,nra,3,nrb) - val3

            ! yaya=ybyb=-yayb
            if(doA) hessian(2,nra,2,nra)=hessian(2,nra,2,nra) + val4
            if(doB) hessian(2,nrb,2,nrb)=hessian(2,nrb,2,nrb) + val4
            if(doAB)hessian(2,nra,2,nrb)=hessian(2,nra,2,nrb) - val4

            ! yaza=ybzb=-yazb=-zayb
            if(doA) hessian(2,nra,3,nra)=hessian(2,nra,3,nra) + val5
            if(doB) hessian(2,nrb,3,nrb)=hessian(2,nrb,3,nrb) + val5
            if(doAB)hessian(3,nra,2,nrb)=hessian(3,nra,2,nrb) - val5
            if(doAB)hessian(2,nra,3,nrb)=hessian(2,nra,3,nrb) - val5


            ! zaza=zbzb=-zazb
            if(doA) hessian(3,nra,3,nra)=hessian(3,nra,3,nra) + val6
            if(doB) hessian(3,nrb,3,nrb)=hessian(3,nrb,3,nrb) + val6
            if(doAB)hessian(3,nra,3,nrb)=hessian(3,nra,3,nrb) - val6

         enddo! nrb=nra,natoms-1
      enddo   ! nra=2,natoms

c---------------------------------------------------------------------
      end
c
c end of nucrep_hess
c======================================================================
c     subroutine intonHess  over IAT.NE.JAT
c
      subroutine intonHess(natoms,  inx,  ifield,  ifgrad,  xfld,
     2                   datbas,  datnuc,ncs,     ncf,     ntri,
     3                   dn,      dw,    athess)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c This routine calculates 2nd derivatives of the  overlap, kinetic,
c dipole & quadrupole integrals
c
c last two for an external electric field & field gradient
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c                   (1 if present, 0 if not)
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix
c dw(ntri)        - weighted density matrix
c INTENT(OUT)
c athess(3,natoms,3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension dn(ntri),dw(ntri)   ! density & weight density
      dimension athess(3,natoms,3,natoms)
      parameter (mxsh=28,mxsh2=mxsh**2,mxsh6=mxsh2*6)
c  28 is the maximum number of angular momentum components in a shell
c  for the Cartesian i-type shell (i28). A given 1-electron
c  integral has only 6 independent second derivatives, say the upper
c  half of the second derivative matrix with respect to the coordinates
c  of the first center, d2/dAx2, d2/dAxdAy, d2/dAxdAz, d2/dAy2,
c  d2/dAudAz, d2/dAz2. The other second derivatives follow from the
c  symmetry of the derivation with respect to  interchange of the
c  coordinates: d2/dAydAx=d2/dAxdAy. The derivatives with respect to the
c  coordinates of the second center follow from the translational
c  invariance condition: d2/dAidBj=-d2/dAidAj, and d2/dBidBj=d2/dAidAj
clocal :
      dimension s3(mxsh6),t3(mxsh6),f3(mxsh6)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len6=len*6
               do 45 jgc=0,ngcjx
                 if(iat.ne.jat) then
                    call zeroit(s3,len6)
                    call zeroit(t3,len6)
c
c................... overlap integral derivatives
c
                   key=1
                   call onehess(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
c
c................... kinetic energy integral derivatives
c
                   key=2
                   call onehess(i, j,  igc, jgc, datbas,
     1                          t3,inx,key, kk,  k2)
c
c................... an external electric field & field gradient.....
c
                   key=1
                   if(ifield.eq.1) then
c ......................an external electric field
c                       (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)
c
                      do k1=1,3
                         field=xfld(k1)
                         if(field.ne.zero) then
                           call zeroit(f3,len6)
c
                            call onehess(i, j,  igc, jgc, datbas,
     *                                   f3,inx,key, k1,  k2)
c
c   replace with DAXPY
                            do l=1,len6
                               t3(l)=t3(l)+f3(l)*field
                            enddo
                         endif    ! field.ne.zero
                      enddo       ! k1=1,3
                   endif          ! if(ifield.eq.1) then
c
                   if(ifgrad.eq.1) then
c ......................an external electric field gradient
c                       (i10 /xy/ j00) + (i00 /xy/ j10)
c                       where xy=xx,yx,yy,zx,zy,zz
c
                      do k1=4,9
                         field=xfld(k1)
                         if(field.ne.zero) then
                            call zeroit(f3,len6)
c
                            call onehess(i, j,  igc, jgc, datbas,
     *                                   f3,inx,key, k1,  k2)
c
c   replace with DAXPY
                            do l=1,len6
                               t3(l)=t3(l)+f3(l)*field
                            enddo
                         endif    ! field.ne.zero
                      enddo       ! k1=4,9
                   endif          ! if(ifgrad.eq.1) then
                 endif            ! iat.ne.jat
c
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
c                    iijx=iij
c                    iijy=iijx+len
c                    iijz=iijx+len*2
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
                     dnij=dn(ij)
                     dwij=dw(ij)
                     if(iff.ne.jff) then
                        dnij=2.d0*dnij
                        dwij=2.d0*dwij
                     endif
c
ccccccc   IF(IAT.NE.JAT) THEN
             icrij=iij-len
             do icr1=1,3
                do icr2=icr1,3
                   icrij=icrij+len
c.................................................................
c The contribution to the Hessian is -Tr[S"W], where W=2CeC+,
c and S" is the second derivative of the overlap matrix
c s3 is d2S/dAidAj  where i,j=1,2 or 3 (x,y,z), and Ai, Aj are the
c coordinates of the first atom. The derivatives for the second atom
c are obtained from the translational invariance conditions:
c d2S/dAidBj=-d2S/dAidBi; d2S/dBidBj=d2S/dAidAj
c.................................................................
c
          s3dw=t3(icrij)*dnij-s3(icrij)*dwij
c
c
c     write(6,*)
c    *' d2S,Tij/dAxdAy=',iat,icr1,icr2,' D2=', s3(icrij), t3(icrij)
c...................................................................
c
c        aiaj=s3dw
c        bibj=aiaj
c
c        aibj=-aiaj
c        ajbi=-aiaj
c...................................................................
cpp       if(iat.gt.jat) then
cpp         athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+s3dw
cpp         athess(icr1,jat,icr2,iat)=athess(icr1,jat,icr2,iat)-s3dw
cpp         athess(icr2,jat,icr1,iat)=athess(icr2,jat,icr1,iat)-s3dw
cpp         athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,iat)+s3dw
cpp       else
cpp         athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+s3dw
cpp         athess(icr1,iat,icr2,jat)=athess(icr1,iat,icr2,jat)-s3dw
cpp         athess(icr2,iat,icr1,jat)=athess(icr2,iat,icr1,jat)-s3dw
cpp         athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,iat)+s3dw
cpp       end if
c...................................................................
          if(iat.gt.jat) then
            athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+s3dw
            athess(icr1,jat,icr2,iat)=athess(icr1,jat,icr2,iat)-s3dw
            athess(icr2,jat,icr1,iat)=athess(icr2,jat,icr1,iat)-s3dw
            athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+s3dw
          endif
          if(iat.lt.jat) then
            athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+s3dw
            athess(icr1,iat,icr2,jat)=athess(icr1,iat,icr2,jat)-s3dw
            athess(icr2,iat,icr1,jat)=athess(icr2,iat,icr1,jat)-s3dw
            athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+s3dw
          end if
c...................................................................
                enddo
             enddo
ccccccc   ENDIF         !    (IAT.NE.JAT) THEN
c
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine onehess(ics,  jcs,  igc,  jgc,  datbas,
     *                   sb,   inx,  key,  k1,   k2)
      implicit real*8 (a-h,o-z)
c this subroutine calculates the one-electron integral 2nd derivatives
c for a contracted pair of shells, and one particular general
c contraction (cf. onel wich includes the general contraction loop)
c ARGUMENTS
c INTENT(IN)
c ics,jcs are the indices of contracted shells
c igc and jgc are the general contraction indices (starting at 0)
c datbas contains the shell info
c inx is the contraction info
c key gives the kind of matrix element wanted (together with k1 and k2
c the latter can take the values 1,2 or 3 =x,y and z
c key=1 for overlap, 2 for kinetic energy; other functionality
c  (dipole moments etc.) is not currently included; if included
c c key=1 for overlap, 2 for kinetic, 3 for dipole, 4 for r**2, 5 for
c second moments
c The nuclear derivatives (key=6 for potential at nucleus nn,
c key=7 for electric field at nn, 8 for field gradient at nn
c are in a separate routine intof7
c INTENT(OUT)
c sb holds the integral derivatives
c format is sb(6,jlen1,ilen1)
c  where ilen1 is the size of contracted shell ics (3 for P, 5 for D)
c  jlen1 is the size of contracted shell jcs,
c  and the last subscript denotes the derivatives (x,y,z)
c
      common /forcdbl/ thre1,thre2,tchf
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension inx(12,*)
      dimension datbas(13,*)
      dimension sb(*)
clocal: rb is 6*28**2=4704
      parameter (mxsh=28,mxsh2=mxsh**2,mxsh6=mxsh2*6)
      dimension rb(mxsh6)
      dimension xa(3),xb(3)
      data twopi/0.6366197723675d0/  !  this is 1/sqrt(12)
c----------------------------------------------------------------------
      ldim=6
c
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
c
      ilen=inx(3,ics)
      jlen=inx(3,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
      if (ityp.eq.6) ilen1=10
      if (jtyp.eq.6) jlen1=10
c
      if(ityp.eq.11) ilen1=15
      if(ityp.eq.12) ilen1=21
      if(ityp.eq.13) ilen1=28
c
      if(jtyp.eq.11) jlen1=15
      if(jtyp.eq.12) jlen1=21
      if(jtyp.eq.13) jlen1=28
c
      len1=ilen1*jlen1
c
      ityp1=ityp
      jtyp1=jtyp
      if(ityp.ge.5) ityp1=ityp-1
      if(jtyp.ge.5) jtyp1=jtyp-1
      if(ityp.ge.7) ityp1=ityp-2
      if(jtyp.ge.7) jtyp1=jtyp-2
      if(ityp.gt.10) ityp1=ityp-5
      if(jtyp.gt.10) jtyp1=jtyp-5
c
      lxin=ldim*len1
c
c----------------------------------------------------------------------
c beg. and end of contraction :
c
      ia=inx(1,ics)+1
      ja=inx(1,jcs)+1
      ie=inx(5,ics)
      je=inx(5,jcs)
c
      do 11 i=1,lxin
  11  sb(i)=zero
c----------------------------------------------------------------------
c
      do 600 i=ia,ie
         a=datbas(1,i)
         sqa=sqrt(a)
         csa=datbas(igc+2,i)
c    cpa is needed only for the l type
         cpa=datbas(3,i)
         do 200 l=1,3
         xa(l)=datbas(l+10,i)
  200    continue
c
         do 500 j=ja,je
            b=datbas(1,j)
            sqb=sqrt(b)
            csb=datbas(jgc+2,j)
            cpb=datbas(3,j)
c
            apb=a+b
            r=zero
            e=a*b/apb
            do 300 l=1,3
               xb(l)=datbas(l+10,j)
               r=r+(xa(l)-xb(l))**2
  300       continue
c
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
c
c key.eq.1.or.key.eq.2  only
c
            call inconh(ityp1,jtyp1,key,k1,k2,a,b,s0,xa,xb,rb)
c
c
c
          ij=0
          do 410 ld=1,ldim
            do 400 i1=1,ilen1
               coefi=csa
               if (ityp.eq.3.and.i1.gt.1) coefi=cpa
            do 400 j1=1,jlen1
               coefj=csb
               if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
               ij=ij+1
               sb(ij)=sb(ij)+rb(ij)*coefi*coefj
  400          continue
  410       continue
  500    continue
  600 continue
c----------------------------------------------------------------------
c transformation d6->d5, f10->f7
c
      incre=len1
      do 710 ld=1,ldim
      iadd=1+(ld-1)*incre
c
        if(ityp.eq.4) then
           call dtran1a(rb,sb(iadd),jlen1)
        end if
        if(ityp.eq.6) then
           call ftran1a(rb,sb(iadd),jlen1)
        endif
        if(ityp.eq.11) then
           call gtran1a(rb,sb(iadd),jlen1)
        endif
        if(ityp.eq.12) then
           call htran1a(rb,sb(iadd),jlen1)
        endif
        if(ityp.eq.13) then
c           call itran1a(rb,sb(iadd),jlen1)
        endif
c
c at this point d and f functions are transformed so the length is ilen
c
        if(jtyp.eq.4) then
           call dtran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.6) then
           call ftran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.11) then
           call gtran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.12) then
           call htran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.13) then
c           call itran2a(rb,sb(iadd),ilen)
        end if
  710 continue
c
c end of transformation
c
c move to appropriate location
c
      incre=len1
      ij=0
      do 720 ld=1,ldim
      iadd=(ld-1)*incre
      do 730 ij1=1,len
      ij=ij+1
      sb(ij)=sb(iadd+ij1)
  730 continue
  720 continue
c----------------------------------------------------------------------
      end
c====================================================================
      subroutine inconh(ityp,jtyp,key,k1,k2,a,b,s0,xa,xb,xint)
      implicit real*8 (a-h,o-z)
c    this routine constructs a set of one-electron integral second
c    derivatives over two primitive shells
c   ARGUMENTS:
c
c  INTENT(IN):
c  ityp,jtyp: the types of the shells, s=1,p=2,l=3,d=4,f=5,g=6,h=7,i=8
c    note that only 6-component d and 10-component f are used here
c    transformation to d5 and f7 is later
c  key is the type of the integral: 1=overlap, 2=kinetic
c  (3=a dipole component, 4=r**2,, 5=a second moment component)
c  Currently only 1 and 2 are permitted
c  k1 and k2 are Cartesian components. They can range from 1=x to 3=z
c   they will be used for the dipole and 2nd moment integrals
c  a and b are the Gaussian exponents, sqa and sqb their square roots
c  s0 is the overlap integral of the spherical parts of the Gaussians
c  xa(3) and xb(3) are the center coordinates of the Gaussians
c  INTENT(OUT):
c  xint(1:lenj,1:leni,6) gives the integral second derivatives with
c  respect to the coordinates of the first center. there are only
c  six components computed; the rest of the derivatives can be
c  reconstructed from these 6
c   leni and lenj are the sizes (lengths) of the shells, e.g.
c   3 for P, 4 for L, 6 for D6
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      dimension xint(*)
      dimension xa(3),xb(3),len(8),la(3),lb(3),ll(3)
      dimension mulx(88),muly(88),mulz(88),nfu(9)
      parameter (mxsh=28,mxsh2=mxsh**2,mxsh6=mxsh2*6)
      dimension yint(mxsh6)
c
      data len/1,3,4,6,10,15,21,28/
      data nfu/0,1,4,8,14,24,39,60,88/
c
      data mulx/0, 1,0,0, 0,1,0,0, 2,0,0,1,1,0, 3,2,2,1,1,1,0,0,0,0,
     1          4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, 5,4,4,3,3,3,2,2,2,2,
     2          1,1,1,1,1,0,0,0,0,0,0,  6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,
     3          1,1,1,1,1,1,0,0,0,0,0,0,0/
      data muly/0, 0,1,0, 0,0,1,0, 0,2,0,1,0,1, 0,1,0,2,1,0,3,2,1,0,
     1          0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, 0,1,0,2,1,0,3,2,1,0,
     2          4,3,2,1,0,5,4,3,2,1,0,  0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     3          5,4,3,2,1,0,6,5,4,3,2,1,0/
      data mulz/0, 0,0,1, 0,0,0,1, 0,0,2,0,1,1, 0,0,1,0,1,2,0,1,2,3,
     1          0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, 0,0,1,0,1,2,0,1,2,3,
     2          0,1,2,3,4,0,1,2,3,4,5,  0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     3          0,1,2,3,4,5,0,1,2,3,4,5,6/
c
c
      ln=len(ityp)*len(jtyp)
      ia=nfu(ityp)+1
      ie=nfu(ityp+1)
      ja=nfu(jtyp)+1
      je=nfu(jtyp+1)
c  ln6=ln*6, the total number of 2nd derivative components of all
c  integrals in the shell
c
      ln6=6*ln
      do 10 i=1,ln6
         xint(i)=zero
         yint(i)=zero
   10 continue
      do 20 i=1,3
         la(i)=0
         lb(i)=0
         ll(i)=0
   20 continue
   25 go to(30,40) key
c
   30 continue
c
c this is for an external electric field
c and electric field gradient (k1=1,...9)
      if(k1.ne.0) then
         if(k1.le.3) then
              ll(k1)=1
         else
            if(k1.eq.4) ll(1)=2
            if(k1.eq.6) ll(2)=2
            if(k1.eq.9) ll(3)=2
            if(k1.eq.5) then
               ll(1)=1
               ll(2)=1
            endif
            if(k1.eq.7) then
               ll(1)=1
               ll(3)=1
            endif
            if(k1.eq.8) then
               ll(2)=1
               ll(3)=1
            endif
         endif
      endif
c
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      go to 240
c
   40 call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 50 i=ia,ie
      do 50 j=ja,je
         ij=ij+1
         ft=dble(2*(mulx(j)+muly(j)+mulz(j))+3)*b
         xint(ij)     =yint(ij)*ft
         xint(ij+ln)  =yint(ij+ln)*ft
         xint(ij+2*ln)=yint(ij+2*ln)*ft
         xint(ij+3*ln)=yint(ij+3*ln)*ft
         xint(ij+4*ln)=yint(ij+4*ln)*ft
         xint(ij+5*ln)=yint(ij+5*ln)*ft
   50 continue
      lb(1)=2
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ft=-two*b**2
      do 60 i=1,ln6
   60 xint(i)=xint(i)+yint(i)*ft
      lb(1)=0
      lb(2)=2
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      do 70 i=1,ln6
   70 xint(i)=xint(i)+yint(i)*ft
      lb(2)=0
      lb(3)=2
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      lb(3)=0
      do 80 i=1,ln6
   80 xint(i)=xint(i)+yint(i)*ft
      if (jtyp.lt.4) go to 240
      lb(3)=0
      lb(1)=-2
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 90 i=ia,ie
      do 90 j=ja,je
         ij=ij+1
         ft=dble((mulx(j)*(mulx(j)-1))/2)
         xint(ij)     =xint(ij)     -yint(ij)*ft
         xint(ij+ln)  =xint(ij+ln)  -yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
         xint(ij+3*ln)=xint(ij+3*ln)-yint(ij+3*ln)*ft
         xint(ij+4*ln)=xint(ij+4*ln)-yint(ij+4*ln)*ft
         xint(ij+5*ln)=xint(ij+5*ln)-yint(ij+5*ln)*ft
   90 continue
      lb(1)=0
      lb(2)=-2
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 100 i=ia,ie
      do 100 j=ja,je
         ij=ij+1
         ft=dble((muly(j)*(muly(j)-1))/2)
         xint(ij)     =xint(ij)     -yint(ij)*ft
         xint(ij+ln)  =xint(ij+ln)  -yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
         xint(ij+3*ln)=xint(ij+3*ln)-yint(ij+3*ln)*ft
         xint(ij+4*ln)=xint(ij+4*ln)-yint(ij+4*ln)*ft
         xint(ij+5*ln)=xint(ij+5*ln)-yint(ij+5*ln)*ft
  100 continue
      lb(2)=0
      lb(3)=-2
      call intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      lb(3)=0
      ij=0
      do 110 i=ia,ie
      do 110 j=ja,je
         ij=ij+1
         ft=dble((mulz(j)*(mulz(j)-1))/2)
         xint(ij)     =xint(ij)     -yint(ij)*ft
         xint(ij+ln)  =xint(ij+ln)  -yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
         xint(ij+3*ln)=xint(ij+3*ln)-yint(ij+3*ln)*ft
         xint(ij+4*ln)=xint(ij+4*ln)-yint(ij+4*ln)*ft
         xint(ij+5*ln)=xint(ij+5*ln)-yint(ij+5*ln)*ft
  110 continue
c
  240 continue
c
      end
c====================================================================
      subroutine intarh(ityp,jtyp,a,b,s0,xa,xb,la,lb,l3,xint)
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------------
c
c this routine calculates one-electron integral second derivatives over
c  gaussian shells. It is assumed that the 1s (spherical) part of the
c  gaussian is normalized (this is denoted by 'gaussian').
c
c the integrand is
c (x-xa(1))**la(1)*...(z-xb(3))**lb(3)*x**l3(1)...* z**l3(3)*
c  (x-xa(1))**mulx1*...(z-zb(3))**mulz2*gaussian1*gaussian2
c  i.e. it can include, besides the x,y,z factors associated with the
c  angular momentum functions, extra (x-xa),...(z-zb) factors given by
c  la and lb, as well as extra x,y,z factors given by l3
c  NOTE: la and lb are NOT the Cartesian factors associated with the
c  Cartesian Gaussian function, but extra factors
c
c here xa(1), xa(2), xa(3) are the coordinates of center a
c xb(1),(2),(3) is the same for center b. a and b are the
c exponents of the two gaussians
c the arrays mulx, muly and mulz contain the x,y, z factors for
c whole shells of functions. the shell types available are (at this
c level) s, p, l (this is s and p with common exponent), d, f and g.
c these are cartesian functions and the transformation to spherical
c harmonics takes place later. for instance, the array elements
c nfu(ityp)+1 to nfu(ityp+1) refer to a shell of type itype.
c the types are 1-s,2-p,3-l,4-d,5-d,6-g. e.g. elements 9 to 14 of
c mulx,muly and mulz refer to d functions. mulx(9) is 2, showing that
c the first function in the d shell is x**2. muly and mulz(9) are 0.
c mulx(12) is 1, muly(12) is 1, mulz(12) is 0, i.e. the 4th function
c of a d shell (remember, we begin with 9) is a d(xy).
c gaussian exponents. sqa and sqb are the square roots of a
c and b, resp. s0 is the overlap integral between normalized
c 1s gaussians and it is explicitly evaluated in subroutine onel
c the vector la contains the powers of (x-xa), (y-ya),(z-za).
c lb contains (x-xb) etc. l3 contains the powers of x,y,z.
c the results is built in ((xint(j,i),j=1,jsh),i=1,ish)
c where ish,jsh are the shell sizes: 1,3,4,6,10,15,21,28 for
c s,p,l,d,f,g,h,i
c it should be ok up to and including i functions, although it would
c be a good idea to test it for high ang. mom. functions befor using it
c it uses hermite-gaussian quadrature. note that this is not very
c economical but transparent.
c
c The derivative is taken with respect to the first center ("A") only.
c This is sufficient; the other derivatives can be calculated from
c  e.g. d2/dAxdBy=-d2/dAxdAy;   d2/dBxdBy=d2/dAxdAy
c  parameter list:
c input:
c ityp,jtyp: function types, 1,2,3..8 for s,p,l,d,f,g,h,i
c a,b: Gaussian exponents
c s0: overlap between normalized 1s gaussians
c xa(3),xb(3): orbital centers
c la(3),lb(3),l3(3) : extra Cartesian exponents, see above.
c It should not change these but apparently it changes la
c  This appears to cause no problems because only lb is non-zero
c  in any call
c output:
c  xint(jsh,ish,6)
c
c  call the routine xtor
c   in order to be correct for i functions and a fourth-degree
c   operator (say, hexadecapole), it needs the hermitian roots
c   and weights up to 9th degree; we have it up to the 10th,
c    so it could be expanded to j functions...
c  dimensioning
c  arguments
      dimension xa(3),xb(3),xint(6*784),la(3),lb(3),l3(3)
c  data
      dimension idg(8),nfu(9),mul(3,88)
c  local variables
      dimension xfc(63),yfc(63),zfc(63)
c
      data idg/0,1,1,2,3,4,5,6/
c  nfu stores the beginnings and endings of shells in the arrays
c  Shell type ityp goes from nfu(ityp)+1 to nfu(ityp)
      data nfu/0,1,4,8,14,24,39,60,88/
c  mul(k,ifu) gives the power of the x,y,z (k=1..3) factors in the shells
      data mul/0,0,0,                                              ! s
     1  1,0,0, 0,1,0, 0,0,1,                                        ! p
     2  0,0,0, 1,0,0, 0,1,0, 0,0,1,                                 ! l
     3  2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1,                   ! d6
     4  3,0,0, 2,1,0, 2,0,1, 1,2,0, 1,1,1, 1,0,2, 0,3,0, 0,2,1,
     5  0,1,2, 0,0,3,                                               ! f10
     6  4,0,0, 3,1,0, 3,0,1, 2,2,0, 2,1,1, 2,0,2, 1,3,0, 1,2,1,
     7  1,1,2, 1,0,3, 0,4,0, 0,3,1, 0,2,2, 0,1,3, 0,0,4,            ! g15
     8  5,0,0, 4,1,0, 4,0,1, 3,2,0, 3,1,1, 3,0,2, 2,3,0, 2,2,1,
     9  2,1,2, 2,0,3, 1,4,0, 1,3,1, 1,2,2, 1,1,3, 1,0,4, 0,5,0,
     &  0,4,1, 0,3,2, 0,2,3, 0,1,4, 0,0,5,                          ! h21
     1  6,0,0, 5,1,0, 5,0,1, 4,2,0, 4,1,1, 4,0,2, 3,3,0, 3,2,1,
     2  3,1,2, 3,0,3, 2,4,0, 2,3,1, 2,2,2, 2,1,3, 2,0,4, 1,5,0,
     3  1,4,1, 1,3,2, 1,2,3, 1,1,4, 1,0,5, 0,6,0, 0,5,1, 0,4,2,
     4  0,3,3, 0,2,4, 0,1,5, 0,0,6/                                 ! i28
      iadg=idg(ityp)
      ibdg=idg(jtyp)
      l1x=iadg+la(1)+2   ! the highest possible power of (x-Ax)
      l2x=ibdg+lb(1)     ! the highest possible power of (x-Bx)
      l1y=iadg+la(2)+2   ! the highest possible power of (y-Ay)
      l2y=ibdg+lb(2)     ! the highest possible power of (y-By)
      l1z=iadg+la(3)+2   ! the highest possible power of (z-Az)
      l2z=ibdg+lb(3)     ! the highest possible power of (z-Bz)
      istart=nfu(ityp)+1
      ilen=nfu(ityp+1)-istart+1
      jstart=nfu(jtyp)+1
      jlen=nfu(jtyp+1)-jstart+1
c
      call inta2h(ilen,  jlen,   mul(1,istart), mul(1,jstart), la,
     1            lb,    l3,     a,             b,             xa,
     3            xb,    l1x,    l1y,           l1z,           l2x,
     4            l2y,   l2z,    s0,            xfc,           yfc,
     5            zfc,   xint)
      end
c========================================================================
      subroutine inta2h(ilen,  jlen,   mul1,  mul2,   la,
     2                  lb,    l3,     a,     b,      xa,
     3                  xb,    l1x,    l1y,   l1z,    l2x,
     4                  l2y,   l2z,    s0,    xfc,   yfc,
     5                  zfc,   xint)
      implicit real*8 (a-h,o-z)
c This is the working routine to calculate 1-el. overlap, kinetic etc.
c integral derivatives
c x(i,j)=integral (x-xa(1))**[la(1)+mul1(1,i)]*(x-xa(1))**[la(2)+mul1(2,i)]
c ...(z-xb(3))**[lb(3)+mul2(3,j)] *x**l3(1)*y**l3(2)*z**l3(3) *
c  exp[-a(R-XA)**2]*[exp[-b(R-XB)**2] where R=(x,y,z),
c  XA=(xa(1),xa(2),xa(3);
c  i.e. it can include, besides the x,y,z factors associated with the
c  angular momentum functions, extra (x-xa),...(z-zb) factors, as well
c  as extra x,y,z factors
c  the loop over i,j is (xint(j,i),j=1,jlen),i=1,ilen)
c Arguments
c INTENT(IN)
c ilen,jlen         = the number of components in the shell, e.g. 3 for P
c mul1x,mul1y,mul1z = arrays containing the power of (x-Ax),(y-Ay),(z-Az)
c  in the function. Here Ax,Ay,Az are the coordinates of the first atom,
c  see xa,ya,za as arguments
c mul2x,mul2y.mul2z = the same for (x-Bx),(y-By),(z-Bz)
c la(3)             = extra powers of (x-ax), (y-ay), (z-az)
c lb(3)             = extra powers of (x-bx), (y-by), (z-bz)
c l3(3)             = extra powers of x,y,z
c a,b               = exponents of the Gaussians
c xa(3),xb(3)       = Cartesian coordinates of the centers, Ax..Bz
c l1x,l1y,l1z       = maximum possible power of (x-ax), (y-ay), (z-az)
c l2x,l2y,l2z       = same for (x-bx), (y-by), (z-bz)
c s0                = overlap of the spheerical part of the 2 Gaussians
c
c INTENT - STORAGE
c xfc,yfc,zfc
c
c INTENT(OUT)
c xint(jlen,ilen,k)   = one-electron integrals
c  arguments
c
      dimension xa(3),xb(3),xint(jlen,ilen,6),la(3),lb(3),l3(3)
      dimension mul1(3,ilen),mul2(3,jlen)
c
c  local variables
      dimension xfc(0:l2x,0:l1x),yfc(0:l2y,0:l1y),zfc(0:l2z,0:l1z)
      dimension ixf(3)
      dimension inc1x(6),inc1y(6),inc1z(6)
      dimension inc2x(6),inc2y(6),inc2z(6)
      dimension inc3x(6),inc3y(6),inc3z(6)
      dimension inc4x(6),inc4y(6),inc4z(6)
c
      PARAMETER (Zero=0.0d0,one=1.0d0,two=2.0d0)
c  the following arrays correspond to the perturbation:
c   xx, xy, xz, yy, yz, zz
      data inc1x/ 2, 1, 1, 0, 0, 0/,
     *     inc1y/ 0, 1, 0, 2, 1, 0/,
     *     inc1z/ 0, 0, 1, 0, 1 ,2/
c
      data inc2x/ 0, 1, 1, 0, 0, 0/,
     *     inc2y/ 0,-1, 0, 0, 1, 0/,
     *     inc2z/ 0, 0,-1, 0,-1, 0/
c
      data inc3x/ 0,-1,-1, 0, 0, 0/,
     *     inc3y/ 0, 1, 0, 0,-1, 0/,
     *     inc3z/ 0, 0, 1, 0, 1, 0/
c
      data inc4x/-2,-1,-1, 0, 0, 0/,
     *     inc4y/ 0,-1, 0,-2,-1, 0/,
     1     inc4z/ 0, 0,-1, 0,-1,-2/
c  call the routine xtor
c   in order to be correct for i functions and a fourth-degree
c   operator (say, hexadecapole), it needs the hermitian roots
c   and weights up to 9th degree; we have it up to the 10th,
c   so it could be expanded to j functions...
c
c     if(ilen.eq.3.and.jlen.eq.3) then
c          write(6,*)' ilen,jlen=',ilen,jlen,' exp=',a,b
c     endif
c
      a2=two*a
      call xtor (l1x,l2x,l3(1),a,b,xa(1),xb(1),xfc)
      call xtor (l1y,l2y,l3(2),a,b,xa(2),xb(2),yfc)
      call xtor (l1z,l2z,l3(3),a,b,xa(3),xb(3),zfc)
      k12=0
      do 70 k1=1,3                   ! differentiation wrt ax(k1)
        do 60 k2=k1,3                ! differentiation wrt ax(k2)
          k12=k12+1
          if(k1.eq.k2) then
            do i=1,ilen
              ixf(1)=mul1(1,i)+la(1)
              ixf(2)=mul1(2,i)+la(2)
              ixf(3)=mul1(3,i)+la(3)
              xmult1=ixf(k1)*(ixf(k1)-1)
              xmult2=2*ixf(k1)+1
c
              ix1=ixf(1)+inc1x(k12)
              iy1=ixf(2)+inc1y(k12)
              iz1=ixf(3)+inc1z(k12)
c ix2 and ix3 are only needed for nondiagonal perturbations, e.g. xy
              ix4=ixf(1)+inc4x(k12)
              iy4=ixf(2)+inc4y(k12)
              iz4=ixf(3)+inc4z(k12)
c
              do j=1,jlen
                jxf=mul2(1,j)+lb(1)
                jyf=mul2(2,j)+lb(2)
                jzf=mul2(3,j)+lb(3)
               if (jxf.ge.0.and.jyf.ge.0.and.jzf.ge.0) then
                xint(j,i,k12)=
     1          a2**2*xfc(jxf,ix1)*yfc(jyf,iy1)*zfc(jzf,iz1)
c
                if (ixf(1).ge.0.and.ixf(2).ge.0.and.ixf(3).ge.0) then
                  xint(j,i,k12)=xint(j,i,k12) -
     2         a2*xmult2*xfc(jxf,ixf(1))*yfc(jyf,ixf(2))*zfc(jzf,ixf(3))
                end if
c
                if (ix4.ge.0.and.iy4.ge.0.and.iz4.ge.0) then
                  xint(j,i,k12)=xint(j,i,k12) +
     2         xmult1*xfc(jxf,ix4)*yfc(jyf,iy4)*zfc(jzf,iz4)
                end if
c
                xint(j,i,k12)=xint(j,i,k12)*s0
               endif
              end do                            ! j
            end do                              ! i
          else
            do i=1,ilen
              ixf(1)=mul1(1,i)+la(1)
              ixf(2)=mul1(2,i)+la(2)
              ixf(3)=mul1(3,i)+la(3)
ckw           xmult1=ixf(k1)
ckw           xmult2=ixf(k2)
              xmult2=ixf(k1)
              xmult1=ixf(k2)
c
              ix1=ixf(1)+inc1x(k12)
              iy1=ixf(2)+inc1y(k12)
              iz1=ixf(3)+inc1z(k12)
c
              ix2=ixf(1)+inc2x(k12)
              iy2=ixf(2)+inc2y(k12)
              iz2=ixf(3)+inc2z(k12)
c
              ix3=ixf(1)+inc3x(k12)
              iy3=ixf(2)+inc3y(k12)
              iz3=ixf(3)+inc3z(k12)
c
              ix4=ixf(1)+inc4x(k12)
              iy4=ixf(2)+inc4y(k12)
              iz4=ixf(3)+inc4z(k12)
c     if(ilen.eq.3.and.jlen.eq.3) then
c         write(6,*)' ix1,iy1,iz1=', ix1,iy1,iz1
c         write(6,*)' ix2,iy2,iz2=', ix2,iy2,iz2
c         write(6,*)' ix3,iy3,iz3=', ix3,iy3,iz3
c         write(6,*)' ix4,iy4,iz4=', ix4,iy4,iz4
c     endif
c
              do j=1,jlen
                jxf=mul2(1,j)+lb(1)
                jyf=mul2(2,j)+lb(2)
                jzf=mul2(3,j)+lb(3)
               if (jxf.ge.0.and.jyf.ge.0.and.jzf.ge.0) then
                xint(j,i,k12)=
     1          a2**2*xfc(jxf,ix1)*yfc(jyf,iy1)*zfc(jzf,iz1)
c
                if (ix2.ge.0.and.iy2.ge.0.and.iz2.ge.0) then
                  xint(j,i,k12)=xint(j,i,k12) -
     2         a2*xmult1*xfc(jxf,ix2)*yfc(jyf,iy2)*zfc(jzf,iz2)
                end if
c
                if (ix3.ge.0.and.iy3.ge.0.and.iz3.ge.0) then
                  xint(j,i,k12)=xint(j,i,k12) -
     2         a2*xmult2*xfc(jxf,ix3)*yfc(jyf,iy3)*zfc(jzf,iz3)
                end if
c
                if (ix4.ge.0.and.iy4.ge.0.and.iz4.ge.0) then
                  xint(j,i,k12)=xint(j,i,k12) + xmult1*xmult2*
     2            xfc(jxf,ix4)*yfc(jyf,iy4)*zfc(jzf,iz4)
                end if
c
                xint(j,i,k12)=xint(j,i,k12)*s0
c
               endif       ! (jxf.ge.0.and.jyf.ge.0.and.jzf.ge.0) then
              end do                            ! j
            end do                              ! i
          end if
   60   continue
   70 continue
      end
c======================================================================
c this is only for test of S2 matix
c
      subroutine intonST2all(natoms,  inx,  ifield,  ifgrad,  xfld,
     2                   datbas,  datnuc,ncs,     ncf,     ntri,
     3                   dn,      dw,    over2, key )
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c loop over all atoms including IAT NE JAT
c---------------------------------------------------------------------
c This routine calculates 2nd derivatives of the  overlap, kinetic,
c integrals and forms coresponding matrices
c derivative integrals are calculated using onehess called twice
c once with (ICS,JCS...) and then with (JCS,ICS,...)
c
c This routine is used only for testing purposes
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c                   (1 if present, 0 if not)
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix
c dw(ntri)        - weighted density matrix
c INTENT(OUT)
c athess(3,natoms,3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension dn(ntri),dw(ntri)   ! density & weight density
      dimension over2(ntri,3,natoms,3,natoms)
c
      parameter (mxsh=28,mxsh2=mxsh**2,mxsh6=mxsh2*6)
c
c  28 is the maximum number of angular momentum components in a shell
c  for the Cartesian i-type shell (i28). A given 1-electron
c  integral has only 6 independent second derivatives, say the upper
c  half of the second derivative matrix with respect to the coordinates
c  of the first center, d2/dAx2, d2/dAxdAy, d2/dAxdAz, d2/dAy2,
c  d2/dAudAz, d2/dAz2. The other second derivatives follow from the
c  symmetry of the derivation with respect to  interchange of the
c  coordinates: d2/dAydAx=d2/dAxdAy. The derivatives with respect to the
c  coordinates of the second center follow from the translational
c  invariance condition: d2/dAidBj=-d2/dAidAj, and d2/dBidBj=d2/dAidAj
clocal :
      dimension s3(mxsh6),t3(mxsh6),f3(mxsh6)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len6=len*6
               do 45 jgc=0,ngcjx
c222222222       if(iat.ne.jat) then
                    call zeroit(s3,len6)
                    call zeroit(t3,len6)
c
c................... overlap integral derivatives
c
c input parameter  key=1
                   call onehess(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
                   call onehess(j,i,jgc,igc,datbas,t3,inx,key,kk,k2)
c
c................... kinetic energy integral derivatives
c
c input parameter  key=2
c                  call onehess(i,j,igc,jgc,datbas,t3,inx,key,kk,k2)
c
c................... an external electric field & field gradient.....
c
c                  key=1
c                  if(ifield.eq.1) then
c ......................an external electric field
c                       (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)
c
c                     do k1=1,3
c                        field=xfld(k1)
c                        if(field.ne.zero) then
c                          call zeroit(f3,len6)
c
c                           call onehess(i, j,  igc, jgc, datbas,
c    *                                   f3,inx,key, k1,  k2)
c
c   replace with DAXPY
c                           do l=1,len6
c                              t3(l)=t3(l)+f3(l)*field
c                           enddo
c                        endif    ! field.ne.zero
c                     enddo       ! k1=1,3
c                  endif          ! if(ifield.eq.1) then
c
c                  if(ifgrad.eq.1) then
c ......................an external electric field gradient
c                       (i10 /xy/ j00) + (i00 /xy/ j10)
c                       where xy=xx,yx,yy,zx,zy,zz
c
c                     do k1=4,9
c                        field=xfld(k1)
c                        if(field.ne.zero) then
c                           call zeroit(f3,len6)
c
c                           call onehess(i, j,  igc, jgc, datbas,
c    *                                   f3,inx,key, k1,  k2)
c
c   replace with DAXPY
c                           do l=1,len6
c                              t3(l)=t3(l)+f3(l)*field
c                           enddo
c                        endif    ! field.ne.zero
c                     enddo       ! k1=4,9
c                  endif          ! if(ifgrad.eq.1) then
c222222222       endif            ! iat.ne.jat
c
c..................................................................
c reorder t3 like s3 is and put it into f3 :

           ijcr=0
           do icr=1,3
              do jcr=icr,3
                 ijcr=ijcr+1
                 ijcrl=(ijcr-1)*len
                 do j3=1,len2
                    do i3=1,len1
                       ji=(j3-1)*len1+i3
                       ij=(i3-1)*len2+j3
                       f3(ijcrl+ij)=t3(ijcrl+ji)
                    enddo
                 enddo
              enddo
           enddo
c..................................................................
c check out if S3 & F3 are the same :
           isame=1
           do iii=1,len6
              differ=s3(iii)-f3(iii)
              if(abs(differ).gt.1.d-6) then
                 isame=-1
                 exit
              endif
           enddo
c..................................................................
           if(isame.eq.-1) then
              write(6,*) ' S3 for & f3 iat=',iat,' jat=',jat,
     *                   ' shells=',i,j,' leni,j=',len1,len2
              write(6,99) (s3(iii),iii=1,len6)
              write(6,*)'  '
              write(6,99) (f3(iii),iii=1,len6)
 99           format(6(f10.6,1x))
           endif
c..................................................................
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
                     dnij=dn(ij)
                     dwij=dw(ij)
                     if(iff.ne.jff) then
                        dnij=2.d0*dnij
                        dwij=2.d0*dwij
                     endif
c
         if(isame.eq.-1) then
            write(6,*)' iij=',iij,' icf=',iff,' jcf=',jff,' ij=',ij
         endif
c
c222222222IF(IAT.NE.JAT) THEN
             icrij=iij-len
             do icr1=1,3
                do icr2=icr1,3
                   icrij=icrij+len
c.................................................................
c The contribution to the Hessian is -Tr[S"W], where W=2CeC+,
c and S" is the second derivative of the overlap matrix
c s3 is d2S/dAidAj  where i,j=1,2 or 3 (x,y,z), and Ai, Aj are the
c coordinates of the first atom. The derivatives for the second atom
c are obtained from the translational invariance conditions:
c d2S/dAidBj=-d2S/dAidBi; d2S/dBidBj=d2S/dAidAj
c...................................................................
          saaij=s3(icrij)
          sbbij=f3(icrij)
          sabij=-saaij
          sbaij=-sbbij
c
           if(isame.eq.-1) then
      write(6,*)'d2S/dAdA=',saaij,' at=',iat,' coor=',icr1,icr2
      write(6,*)'d2S/dBdB=',sbbij,' at=',jat,' coor=',icr1,icr2
           endif
c
      over2(ij,icr1,iat,icr2,iat)=over2(ij,icr1,iat,icr2,iat)+saaij
      over2(ij,icr1,jat,icr2,jat)=over2(ij,icr1,jat,icr2,jat)+sbbij
      over2(ij,icr1,iat,icr2,jat)=over2(ij,icr1,iat,icr2,jat)+sabij
      over2(ij,icr1,jat,icr2,iat)=over2(ij,icr1,jat,icr2,iat)+sbaij
c
      if(icr1.ne.icr2) then
         over2(ij,icr2,iat,icr1,iat)=over2(ij,icr2,iat,icr1,iat)+saaij
         over2(ij,icr2,jat,icr1,jat)=over2(ij,icr2,jat,icr1,jat)+sbbij
         over2(ij,icr2,jat,icr1,iat)=over2(ij,icr2,jat,icr1,iat)+sabij
         over2(ij,icr2,iat,icr1,jat)=over2(ij,icr2,iat,icr1,jat)+sbaij
      endif

                enddo
             enddo
c222222222ENDIF         !    (IAT.NE.JAT) THEN
c
   40             continue
                  jfu=jfu+len2
   45          continue          ! end of loop over gen.con.belonging to jcs
   50       continue             ! end of loop over jcs
            ifu=ifu+len1
   55    continue                ! end of loop over gen.con.belonging to ics
   60 continue                   ! end of loop over ics
c
c---------------------------------------------------------------------
c check if any of one center elements od S2 or T2 are non-zero :
c
c
         do ics=1,ncs
            iat=inx(2,ics)
            do jcs=1,ics
               jat=inx(2,jcs)
               if(iat.eq.jat) then
                  do icr=1,3
                  do jcr=icr,3
                     do icf=inx(11,ics)+1,inx(10,ics)
                     do jcf=inx(11,jcs)+1,inx(10,jcs)
                        if(icf.ge.jcf) then
                           ijcf=icf*(icf-1)/2 +jcf
                        else
                           ijcf=jcf*(jcf-1)/2 +icf
                        endif
                        s2ij=over2(ijcf,icr,iat,jcr,jat)
                        if(abs(s2ij).gt. 1.d-9) then
      write(6,*)
     * 'one center element of the 2nd derivative matrix is NOT zero'
                        write(6,*)'icf,jcf=',icf,jcf,' elem=',s2ij
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
               endif
            enddo
         enddo
c
c---------------------------------------------------------------------
c
      end
c======================================================================
c this is only for test of S2 matix
c
      subroutine intonST2neq(natoms,  inx,  ifield,  ifgrad,  xfld,
     2                   datbas,  datnuc,ncs,     ncf,     ntri,
     3                   dn,      dw,    over2, key )
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c loop over IAT NE JAT only
c---------------------------------------------------------------------
c This routine calculates 2nd derivatives of the  overlap, kinetic,
c integrals and forms coresponding matrices
c derivative integrals are calculated using onehess called twice
c once with (ICS,JCS...) and then with (JCS,ICS,...)
c
c This routine is used only for testing purposes
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c                   (1 if present, 0 if not)
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix
c dw(ntri)        - weighted density matrix
c INTENT(OUT)
c athess(3,natoms,3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension dn(ntri),dw(ntri)   ! density & weight density
      dimension over2(ntri,3,natoms,3,natoms)
c
      parameter (mxsh=28,mxsh2=mxsh**2,mxsh6=mxsh2*6)
c
c  28 is the maximum number of angular momentum components in a shell
c  for the Cartesian i-type shell (i28). A given 1-electron
c  integral has only 6 independent second derivatives, say the upper
c  half of the second derivative matrix with respect to the coordinates
c  of the first center, d2/dAx2, d2/dAxdAy, d2/dAxdAz, d2/dAy2,
c  d2/dAudAz, d2/dAz2. The other second derivatives follow from the
c  symmetry of the derivation with respect to  interchange of the
c  coordinates: d2/dAydAx=d2/dAxdAy. The derivatives with respect to the
c  coordinates of the second center follow from the translational
c  invariance condition: d2/dAidBj=-d2/dAidAj, and d2/dBidBj=d2/dAidAj
clocal :
      dimension s3(mxsh6),t3(mxsh6),f3(mxsh6)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len6=len*6
               do 45 jgc=0,ngcjx
      IF(IAT.NE.JAT) THEN
                    call zeroit(s3,len6)
                    call zeroit(t3,len6)
c
                    call onehess(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
                    call onehess(j,i,jgc,igc,datbas,t3,inx,key,kk,k2)
c
c..................................................................
c reorder t3 like s3 is and put it into f3 :

           ijcr=0
           do icr=1,3
              do jcr=icr,3
                 ijcr=ijcr+1
                 ijcrl=(ijcr-1)*len
                 do j3=1,len2
                    do i3=1,len1
                       ji=(j3-1)*len1+i3
                       ij=(i3-1)*len2+j3
                       f3(ijcrl+ij)=t3(ijcrl+ji)
                    enddo
                 enddo
              enddo
           enddo
c..................................................................
c check out if S3 & F3 are the same :
           isame=1
           do iii=1,len6
              differ=s3(iii)-f3(iii)
              if(abs(differ).gt.1.d-6) then
                 isame=-1
                 exit
              endif
           enddo
c..................................................................
           if(isame.eq.-1) then
              write(6,*) ' S3 for & f3 iat=',iat,' jat=',jat,
     *                   ' shells=',i,j,' leni,j=',len1,len2
              write(6,99) (s3(iii),iii=1,len6)
              write(6,*)'  '
              write(6,99) (f3(iii),iii=1,len6)
 99           format(6(f10.6,1x))
           endif
      ENDIF         !        (IAT.NE.JAT) THEN
c..................................................................
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
                     dnij=dn(ij)
                     dwij=dw(ij)
                     if(iff.ne.jff) then
                        dnij=2.d0*dnij
                        dwij=2.d0*dwij
                     endif
c
      IF(IAT.NE.JAT) THEN
         if(isame.eq.-1) then
            write(6,*)' iij=',iij,' icf=',iff,' jcf=',jff,' ij=',ij
         endif
             icrij=iij-len
             do icr1=1,3
                do icr2=icr1,3
                   icrij=icrij+len
c.................................................................
c The contribution to the Hessian is -Tr[S"W], where W=2CeC+,
c and S" is the second derivative of the overlap matrix
c s3 is d2S/dAidAj  where i,j=1,2 or 3 (x,y,z), and Ai, Aj are the
c coordinates of the first atom. The derivatives for the second atom
c are obtained from the translational invariance conditions:
c d2S/dAidBj=-d2S/dAidBi; d2S/dBidBj=d2S/dAidAj
c...................................................................
          saaij=s3(icrij)
          sbbij=f3(icrij)
          sabij=-saaij
          sbaij=-sbbij
c
           if(isame.eq.-1) then
      write(6,*)'d2S/dAdA=',saaij,' at=',iat,' coor=',icr1,icr2
      write(6,*)'d2S/dBdB=',sbbij,' at=',jat,' coor=',icr1,icr2
           endif
c
      over2(ij,icr1,iat,icr2,iat)=over2(ij,icr1,iat,icr2,iat)+saaij
      over2(ij,icr1,jat,icr2,jat)=over2(ij,icr1,jat,icr2,jat)+sbbij
      over2(ij,icr1,iat,icr2,jat)=over2(ij,icr1,iat,icr2,jat)+sabij
      over2(ij,icr1,jat,icr2,iat)=over2(ij,icr1,jat,icr2,iat)+sbaij
c
      if(icr1.ne.icr2) then
         over2(ij,icr2,iat,icr1,iat)=over2(ij,icr2,iat,icr1,iat)+saaij
         over2(ij,icr2,jat,icr1,jat)=over2(ij,icr2,jat,icr1,jat)+sbbij
         over2(ij,icr2,jat,icr1,iat)=over2(ij,icr2,jat,icr1,iat)+sabij
         over2(ij,icr2,iat,icr1,jat)=over2(ij,icr2,iat,icr1,jat)+sbaij
      endif

                enddo
             enddo
      ENDIF         !        (IAT.NE.JAT) THEN
   40             continue
                  jfu=jfu+len2
   45          continue          ! end of loop over gen.con.belonging to jcs
   50       continue             ! end of loop over jcs
            ifu=ifu+len1
   55    continue                ! end of loop over gen.con.belonging to ics
   60 continue                   ! end of loop over ics
c
c---------------------------------------------------------------------
c check if any of one center elements od S2 or T2 are non-zero :
c
c
         do ics=1,ncs
            iat=inx(2,ics)
            do jcs=1,ics
               jat=inx(2,jcs)
               if(iat.eq.jat) then
                  do icr=1,3
                  do jcr=icr,3
                     do icf=inx(11,ics)+1,inx(10,ics)
                     do jcf=inx(11,jcs)+1,inx(10,jcs)
                        if(icf.ge.jcf) then
                           ijcf=icf*(icf-1)/2 +jcf
                        else
                           ijcf=jcf*(jcf-1)/2 +icf
                        endif
                        s2ij=over2(ijcf,icr,iat,jcr,jat)
                        if(abs(s2ij).gt. 1.d-9) then
      write(6,*)
     * 'one center element of the 2nd derivative matrix is NOT zero'
                        write(6,*)'icf,jcf=',icf,jcf,' elem=',s2ij
                        endif
                     enddo
                     enddo
                  enddo
                  enddo
               endif
            enddo
         enddo
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine intonHall(natoms,  inx,  ifield,  ifgrad,  xfld,
     2                   datbas,  datnuc,ncs,     ncf,     ntri,
     3                   dn,      dw,    athess)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c This routine calculates 2nd derivatives of the  overlap, kinetic,
c dipole & quadrupole integrals
c
c last two for an external electric field & field gradient
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c                   (1 if present, 0 if not)
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix
c dw(ntri)        - weighted density matrix
c INTENT(OUT)
c athess(3,natoms,3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension dn(ntri),dw(ntri)   ! density & weight density
      dimension athess(3,natoms,3,natoms)
      parameter (mxsh=28,mxsh2=mxsh**2,mxsh6=mxsh2*6)
c  28 is the maximum number of angular momentum components in a shell
c  for the Cartesian i-type shell (i28). A given 1-electron
c  integral has only 6 independent second derivatives, say the upper
c  half of the second derivative matrix with respect to the coordinates
c  of the first center, d2/dAx2, d2/dAxdAy, d2/dAxdAz, d2/dAy2,
c  d2/dAudAz, d2/dAz2. The other second derivatives follow from the
c  symmetry of the derivation with respect to  interchange of the
c  coordinates: d2/dAydAx=d2/dAxdAy. The derivatives with respect to the
c  coordinates of the second center follow from the translational
c  invariance condition: d2/dAidBj=-d2/dAidAj, and d2/dBidBj=d2/dAidAj
clocal :
      dimension s3(mxsh6),t3(mxsh6),f3(mxsh6)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len6=len*6
               do 45 jgc=0,ngcjx
cccccccccccccc   if(iat.ne.jat) then
                    call zeroit(s3,len6)
                    call zeroit(t3,len6)
c
c................... overlap integral derivatives
c
                   key=1
                   call onehess(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
c
c................... kinetic energy integral derivatives
c
                   key=2
                   call onehess(i, j,  igc, jgc, datbas,
     1                          t3,inx,key, kk,  k2)
c
c................... an external electric field & field gradient.....
c
                   key=1
                   if(ifield.eq.1) then
c ......................an external electric field
c                       (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)
c
                      do k1=1,3
                         field=xfld(k1)
                         if(field.ne.zero) then
                           call zeroit(f3,len6)
c
                            call onehess(i, j,  igc, jgc, datbas,
     *                                   f3,inx,key, k1,  k2)
c
c   replace with DAXPY
                            do l=1,len6
                               t3(l)=t3(l)+f3(l)*field
                            enddo
                         endif    ! field.ne.zero
                      enddo       ! k1=1,3
                   endif          ! if(ifield.eq.1) then
c
                   if(ifgrad.eq.1) then
c ......................an external electric field gradient
c                       (i10 /xy/ j00) + (i00 /xy/ j10)
c                       where xy=xx,yx,yy,zx,zy,zz
c
                      do k1=4,9
                         field=xfld(k1)
                         if(field.ne.zero) then
                            call zeroit(f3,len6)
c
                            call onehess(i, j,  igc, jgc, datbas,
     *                                   f3,inx,key, k1,  k2)
c
c   replace with DAXPY
                            do l=1,len6
                               t3(l)=t3(l)+f3(l)*field
                            enddo
                         endif    ! field.ne.zero
                      enddo       ! k1=4,9
                   endif          ! if(ifgrad.eq.1) then
cccccccccccccc   endif            ! iat.ne.jat
c
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
c                    iijx=iij
c                    iijy=iijx+len
c                    iijz=iijx+len*2
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
c
                     dnij=dn(ij)
                     dwij=dw(ij)
                     if(iff.ne.jff) then
                        dnij=2.d0*dnij
                        dwij=2.d0*dwij
                     endif
c
cccccc    IF(IAT.NE.JAT) THEN
             icrij=iij-len
             do icr1=1,3
                do icr2=icr1,3
                   icrij=icrij+len
c.................................................................
c The contribution to the Hessian is -Tr[S"W], where W=2CeC+,
c and S" is the second derivative of the overlap matrix
c s3 is d2S/dAidAj  where i,j=1,2 or 3 (x,y,z), and Ai, Aj are the
c coordinates of the first atom. The derivatives for the second atom
c are obtained from the translational invariance conditions:
c d2S/dAidBj=-d2S/dAidBi; d2S/dBidBj=d2S/dAidAj
c.................................................................
c
          s3dw=t3(icrij)*dnij-s3(icrij)*dwij
c
c...................................................................
c      if(iat.gt.jat) then
c           athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+s3dw
c           athess(icr1,jat,icr2,iat)=athess(icr1,jat,icr2,iat)-s3dw
c           athess(icr2,jat,icr1,iat)=athess(icr2,jat,icr1,iat)-s3dw
c           athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+s3dw
c      endif
c
c      if(iat.lt.jat) then
c           athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+s3dw
c           athess(icr1,iat,icr2,jat)=athess(icr1,iat,icr2,jat)-s3dw
c           athess(icr2,iat,icr1,jat)=athess(icr2,iat,icr1,jat)-s3dw
c           athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+s3dw
c      end if
c...................................................................
          saaij=s3dw
          sbbij=s3dw
          sabij=-saaij
          sbaij=-sbbij
c
      athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+saaij
      athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+sbbij
      athess(icr1,iat,icr2,jat)=athess(icr1,iat,icr2,jat)+sabij
      athess(icr1,jat,icr2,iat)=athess(icr1,jat,icr2,iat)+sbaij
c
      if(icr1.ne.icr2) then
         athess(icr2,iat,icr1,iat)=athess(icr2,iat,icr1,iat)+saaij
         athess(icr2,jat,icr1,jat)=athess(icr2,jat,icr1,jat)+sbbij
         athess(icr2,jat,icr1,iat)=athess(icr2,jat,icr1,iat)+sabij
         athess(icr2,iat,icr1,jat)=athess(icr2,iat,icr1,jat)+sbaij
      endif
c...................................................................
                enddo
             enddo
cccccc    ENDIF         !    (IAT.NE.JAT) THEN
c
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine elenuc_hess(natoms,inx,datbas,datnuc,ncs,ncf,ntri,
     *                       dn,athess)
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical donra
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension dn(ntri)
      dimension athess(3,natoms,3,natoms)
c---------------------------------------------------------------------
c include charged dummies :
c
      call getival('ndum1',ndum1)
c---------------------------------------------------------------------
ccc   do nra=1,natoms
      do nra=1,natoms+ndum1
         donra=.true.
         if(nra.gt.natoms) donra=.false.
         call intoh7(donra, 
     *               nra,   inx,  datbas,  datnuc,  natoms,
     1               ncs,   ncf,  ntri,    dn,      athess)
      enddo
c
      end
c======================================================================
      subroutine intoh7(donra,
     *                  nra,   inx,  datbas,  datnuc,  natoms,
     1                  ncs,   ncf,  ntri,    dn,      athess)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c
c This routine calculates 2nd derivatives of the following integrals :
c
c             nuclear atraction & Hellmann-Feynman forces
c
c---------------------------------------------------------------------
c ARGUMENTS:
c nra             = reference atom
c inx(12,ncs)     = basis set data
c ifield & ifgrad = indicate presence of an ext.electric field & f.grad.
c datbas(13,nsh)  = basis set data
c datnuc(5,natoms)= nuclear data
c natoms          = the number of atoms
c ncs,ncf,ntri    = contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        = density matrix
c athess(3,natoms,3,natoms)= nuclear Hessian - INPUT/OUTPUT
c---------------------------------------------------------------------
      logical donra
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension dn(ntri)                         ! density
      dimension athess(3,natoms,3,natoms)        ! Hessian contribution
clocal :
      dimension sa(6*784),sb(6*784),sc(9*784) ! 6*(28x28) (i|i)
c
c  sa, sb and sab are implicitly dimensioned (unfortunately) as, say,
c  sa(j,i,6), sb(j,i,6) and sc(j,i,3,3)
c  Here i and j go over the shell components (e.g. s,x,y,z for an L
c  shell, and 6 is the number of perturbations wrt the coordinates :
c  xx,xy,xz,yy,yz,zz
c  sa gives the derivatives wrt the first atom, sb wrt the second,
c  and sc gives the mixed derivatives. In matrix form, the mixed
c  derivative block is defined as (Bj,Ai), i.e. the order is
c  d2/dBxdAx, dBydAx, dBzdAx, dBxdAy, dBydAy, dBzdAy, dBxdAz etc.
c  The derivatives with respect to the coordinates of C (The Hellmann-
c  Feynman Hessian, i.e. the derivatives of the Hellmann-Feynman forces)
c  are obtained from the translational invariance
c
      dimension xra(3)
      parameter (zero=0.0d0)
c---------------------------------------------------------------------
c     write(6,*)' density from elenuc_hess'
c     call drumh(dn,ncf, 6  ,'DENSITY ')
c---------------------------------------------------------------------
c note : for this contribution to the hessian the atomic charge must
c        be taken with a negative sign
c
      zza=-datnuc(1,nra)
c---------------------------------------------------------------------
c
      xra(1)=datnuc(2,nra)
      xra(2)=datnuc(3,nra)
      xra(3)=datnuc(4,nra)
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               do 45 jgc=0,ngcjx
c
c................... nuclear attraction derivatives
c
c
                     call oneh7(i,   j,   igc,   jgc,  datbas,
     1                          inx, xra, zza,    sa,  sb, sc)
c
c
c        write(6,*)
c    *' shells=',i,j,' (A|C|B)=',iat,nra,jat,' ij=',ij,' dij=',dij
c          write(6,*)' (Da|c|b) '
c          write(6,66) (sa(iii),iii=1,6)
c          write(6,*)' (a|c|Db) '
c          write(6,66) (sb(iii),iii=1,6)
c          write(6,*)' (Da|c|Db) '
c          write(6,99) (sc(iii),iii=1,9)
c 66  format(6(f12.6,2x))
c 99  format(3(f12.6,2x))
c.................................................................
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
                     dij=dn(ij)
                     if(iff.ne.jff) dij=dij+dij
c
c  iat,iat  and jat,jat terms
c----------------------------------------------------------
c
         dijzza=dij*zza
c
c----------------------------------------------------------
         icrij=iij-len
         do icr1=1,3
            do icr2=icr1,3
               icrij=icrij+len
c
               t1=sa(icrij)*dijzza
               t2=sb(icrij)*dijzza
c
        if(icr1.EQ.icr2) then
          athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+t1
          athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+t2
c
          if( donra ) then
c.......     contribution to mixed Hellmann-Feynman terms iat,nra
             athess(icr1,iat,icr2,nra)=athess(icr1,iat,icr2,nra)-t1
             athess(icr2,nra,icr1,iat)=athess(icr2,nra,icr1,iat)-t1
c
c.......     contribution to mixed Hellmann-Feynman terms jat,nra
             athess(icr1,jat,icr2,nra)=athess(icr1,jat,icr2,nra)-t2
             athess(icr2,nra,icr1,jat)=athess(icr2,nra,icr1,jat)-t2
c
c.......     contribution to the pure Hellmann-Feynman term
             athess(icr1,nra,icr2,nra)=athess(icr1,nra,icr2,nra)+t1+t2
          endif
c
        endif
c
        if(icr1.ne.icr2) then
          athess(icr1,iat,icr2,iat)=athess(icr1,iat,icr2,iat)+t1
          athess(icr2,iat,icr1,iat)=athess(icr2,iat,icr1,iat)+t1
c
          if( donra ) then
             athess(icr1,iat,icr2,nra)=athess(icr1,iat,icr2,nra)-t1
             athess(icr2,nra,icr1,iat)=athess(icr2,nra,icr1,iat)-t1
c
             athess(icr2,iat,icr1,nra)=athess(icr2,iat,icr1,nra)-t1
             athess(icr1,nra,icr2,iat)=athess(icr1,nra,icr2,iat)-t1
          endif
c
          athess(icr2,jat,icr1,jat)=athess(icr2,jat,icr1,jat)+t2
          athess(icr1,jat,icr2,jat)=athess(icr1,jat,icr2,jat)+t2
c
          if( donra ) then
             athess(icr1,jat,icr2,nra)=athess(icr1,jat,icr2,nra)-t2
             athess(icr2,nra,icr1,jat)=athess(icr2,nra,icr1,jat)-t2
c
             athess(icr2,jat,icr1,nra)=athess(icr2,jat,icr1,nra)-t2
             athess(icr1,nra,icr2,jat)=athess(icr1,nra,icr2,jat)-t2
c
             athess(icr1,nra,icr2,nra)=athess(icr1,nra,icr2,nra)+t1+t2
             athess(icr2,nra,icr1,nra)=athess(icr2,nra,icr1,nra)+t1+t2
          endif
c
        endif
c
            end do
         end do
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c mixed iat,jat terms
c
         icrij=iij-len
         do icr1=1,3
            do icr2=1,3
               icrij=icrij+len
c
          t1=sc(icrij)*dijzza
c
c
          athess(icr1,iat,icr2,jat)=athess(icr1,iat,icr2,jat)+t1
          athess(icr2,jat,icr1,iat)=athess(icr2,jat,icr1,iat)+t1
c
          if( donra ) then
c.......       contribution to mixed Hellmann-Feynman terms
             athess(icr1,iat,icr2,nra)=athess(icr1,iat,icr2,nra)-t1
             athess(icr2,nra,icr1,iat)=athess(icr2,nra,icr1,iat)-t1
             athess(icr2,jat,icr1,nra)=athess(icr2,jat,icr1,nra)-t1
             athess(icr1,nra,icr2,jat)=athess(icr1,nra,icr2,jat)-t1
c.......       contribution to the pure Hellmann-Feynman force
             athess(icr1,nra,icr2,nra)=athess(icr1,nra,icr2,nra)+t1
             athess(icr2,nra,icr1,nra)=athess(icr2,nra,icr1,nra)+t1
          endif
c
            enddo
         enddo
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
      end
c=======================================================================
      subroutine oneh7(ics,jcs,igc,jgc,datbas,inx,xra,zza,sa,sb,sc)
c  This routine calculates the second derivative of the electron-nucleus
c  attraction integral for a given contracted shell pair. The derivative
c  is with respect of the coordinates of the positions of the basis
c  functions A and B, corresponding to ics and jcs, respectively.
c  The result is put in the arrays sa, sb and sc. There is, of course,
c  another contribution, the second derivative of the Hellmann-Feynman
c  force, i.e. the second  derivative with respect to the coordinates
c  of the nucleus located by xra, with charge zra.
c  This is calculated from the the translational invariance criterion.
C  Arguments:
C  INPUT:
c  ics,jcs: contracted shell indices
c  igc,jgc: depth of general contaction (0 is unsegmented)
c  datbas: basis set data
c  inx: contraction data
c  xra(3): coordinates of the nucleus
c  zza: charge of the nucleus
C OUTPUT:
c  sa(6*maxshell) = 2nd derivatives of the 1-el. integrals with respect
c                   to the coordinates of the first center (ics)
c  sb(6*maxshell) = 2nd derivatives of the 1-el. integrals with respect
c                   to the coordinates of the second center (jcs)
c  only the derivatives xx,xy,xz,yy,yz,zz are given for the above two
c  sc(9*maxshell) = mixed 2nd derivatives of the 1-el. integral wrt the
c                   coordinates of the first nd second center
c  sa and sb are dimensioned implicitly as sa(jlen,ilen,6) and sb(jlen,ilen,6)
c  Here jlen and ilen are the shell sizes for the second and the first
c  function, respectively
c  sc is dimensioned implicitly as sc(j,i,l,k)  where i=1,ilen, j=1,jlen,
c  k,l=1,2,3   Using the notation I=<i|1/rC|j>, and the coordinates of
c  the center of |i>  Ak , coordinates of |j> Bl,
c  this is the derivative d2 I/(dAk dBl)
      implicit real*8 (a-h,o-z)
      common /forcdbl/ thre1,thre2,tchf
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension inx(12,*)
      dimension datbas(13,*)
      dimension xra(3),sa(*),sb(*),sc(*)
clocal:
      dimension ra(6*784),rb(6*784),rc(9*784)
      dimension xa(3),xb(3)
      data twopi/0.6366197723675d0/  !  this is 1/sqrt(12)
c----------------------------------------------------------------------
      ldim=6
      ldim1=9
c----------------------------------------------------------------------
c
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
c
      ilen=inx(3,ics)
      jlen=inx(3,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
      if (ityp.eq.6) ilen1=10
      if (jtyp.eq.6) jlen1=10
c
      if(ityp.eq.11) ilen1=15
      if(ityp.eq.12) ilen1=21
      if(ityp.eq.13) ilen1=28
c
      if(jtyp.eq.11) jlen1=15
      if(jtyp.eq.12) jlen1=21
      if(jtyp.eq.13) jlen1=28
c
      len1=ilen1*jlen1
c
      ityp1=ityp
      jtyp1=jtyp
      if(ityp.ge.5) ityp1=ityp-1
      if(jtyp.ge.5) jtyp1=jtyp-1
      if(ityp.ge.7) ityp1=ityp-2
      if(jtyp.ge.7) jtyp1=jtyp-2
      if(ityp.gt.10) ityp1=ityp-5
      if(jtyp.gt.10) jtyp1=jtyp-5
c
c----------------------------------------------------------------------
c zero out sa,sb,sc :
c
      do i=1,ldim*len1
         sa(i)=zero
         sb(i)=zero
      enddo
      do i=1,ldim1*len1
         sc(i)=zero
      enddo
c----------------------------------------------------------------------
c beg. and end of contraction :
c
      ia=inx(1,ics)+1
      ja=inx(1,jcs)+1
      ie=inx(5,ics)
      je=inx(5,jcs)
c----------------------------------------------------------------------
c
      do 600 i=ia,ie
         a=datbas(1,i)
         sqa=sqrt(a)
         csa=datbas(igc+2,i)
c    cpa is needed only for the l type
         cpa=datbas(3,i)
         do 200 l=1,3
         xa(l)=datbas(l+10,i)
  200    continue
c
         do 500 j=ja,je
            b=datbas(1,j)
            sqb=sqrt(b)
            csb=datbas(jgc+2,j)
            cpb=datbas(3,j)
c
            apb=a+b
            r=zero
            e=a*b/apb
            do 300 l=1,3
               xb(l)=datbas(l+10,j)
               r=r+(xa(l)-xb(l))**2
  300       continue
c
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
c
            estim=s0*abs(zza)*max(sqa,sqb)
            if(estim.lt.thre1 ) go to 500
c
             call nuchess(ityp1, jtyp1, ilen1, jlen1, a,
     1                    b,     s0,    xa,    xb,    xra,
     2                    ra,    rb,    rc)
c
c
           ij=0
           do ld=1,ldim            ! ldim=6
              do i1=1,ilen1
                 coefi=csa
                 if (ityp.eq.3.and.i1.gt.1) coefi=cpa
                 do j1=1,jlen1
                    coefj=csb
                    if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
                    ij=ij+1
                    coefij=coefi*coefj
                    sa(ij)=sa(ij)+ra(ij)*coefij
                    sb(ij)=sb(ij)+rb(ij)*coefij
                 enddo
              enddo
           enddo
c
           ij=0
           do ld=1,ldim1           ! ldim1=9
              do i1=1,ilen1
                 coefi=csa
                 if (ityp.eq.3.and.i1.gt.1) coefi=cpa
                 do j1=1,jlen1
                    coefj=csb
                    if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
                    ij=ij+1
                    coefij=coefi*coefj
                    sc(ij)=sc(ij)+rc(ij)*coefij
                 enddo
              enddo
           enddo
c
  500    continue
  600 continue
c----------------------------------------------------------------------
c transformation d6->d5, f10->f7
c
      incre=len1
      iadd=1-incre
      do  ld=1,ldim
      iadd=iadd+incre
c  transform ra and rb
        if(ityp.eq.4) then
           call dtran1a(ra,sa(iadd),jlen1)
           call dtran1a(rb,sb(iadd),jlen1)
        end if
        if(ityp.eq.6) then
           call ftran1a(ra,sa(iadd),jlen1)
           call ftran1a(rb,sb(iadd),jlen1)
        endif
      if(ityp.eq.11) then
         call gtran1a(ra,sa(iadd),jlen1)
         call gtran1a(rb,sb(iadd),jlen1)
      end if
      if(ityp.eq.12) then
           call htran1a(ra,sa(iadd),jlen1)
         call htran1a(rb,sb(iadd),jlen1)
      end if
      if(ityp.eq.13) then
c          call itran1a(ra,sa(iadd),jlen1)
c         call itran1a(rb,sb(iadd),jlen1)
      end if
c
c at this point the first function is transformed so the length is ilen
c
        if(jtyp.eq.4) then
           call dtran2a(ra,sa(iadd),ilen)
           call dtran2a(rb,sb(iadd),ilen)
        end if
        if(jtyp.eq.6) then
           call ftran2a(ra,sa(iadd),ilen)
           call ftran2a(rb,sb(iadd),ilen)
        end if
      if(jtyp.eq.11) then
         call gtran2a(ra,sa(iadd),ilen)
         call gtran2a(rb,sb(iadd),ilen)
      end if
      if(jtyp.eq.12) then
         call htran2a(ra,sa(iadd),ilen)
         call htran2a(rb,sb(iadd),ilen)
      end if
      if(jtyp.eq.13) then
c          call itran2a(ra,sa(iadd),ilen)
c          call itran2a(rb,sb(iadd),ilen)
      end if
c
      end do
c
c transform rc
c
      iadd=1-incre
      do  ld=1,ldim1       ! ldim1=9
      iadd=iadd+incre
c  transform rc
        if(ityp.eq.4) then
           call dtran1a(rc,sc(iadd),jlen1)
        end if
        if(ityp.eq.6) then
           call ftran1a(rc,sc(iadd),jlen1)
        endif
      if(ityp.eq.11) then
         call gtran1a(rc,sc(iadd),jlen1)
      end if
      if(ityp.eq.12) then
         call htran1a(rc,sc(iadd),jlen1)
      end if
      if(ityp.eq.13) then
c          call itran1a(rc,sc(iadd),jlen1)
      end if
c
        if(jtyp.eq.4) then
           call dtran2a(rc,sc(iadd),ilen)
        end if
        if(jtyp.eq.6) then
           call ftran2a(rc,sc(iadd),ilen)
        end if
      if(jtyp.eq.11) then
         call gtran2a(rc,sc(iadd),ilen)
      end if
      if(jtyp.eq.12) then
         call htran2a(rc,sc(iadd),ilen)
      end if
      if(jtyp.eq.13) then
c          call itran2a(rc,sc(iadd),ilen)
      end if
c
      end do

c
c end of transformation
c----------------------------------------------------------------------
c move to appropriate location
c
  730 continue
c
      incre=len1
      ij=0
      do ld=1,ldim
        iadd=(ld-1)*incre
        do ij1=1,len
          ij=ij+1
          sa(ij)=sa(iadd+ij1)
          sb(ij)=sb(iadd+ij1)
        enddo
      enddo
c
      incre=len1
      ij=0
      do ld=1,ldim1
        iadd=(ld-1)*incre
        do ij1=1,len
          ij=ij+1
          sc(ij)=sc(iadd+ij1)
        enddo
      enddo
c
c----------------------------------------------------------------------
      end
c======================================================================
      subroutine nuchess(n1,  n2,  ndim1,  ndim2,  a,
     1                   b,   sab, xa,     xb,     xc,
     2                   hessaa, hessbb, hessab)
      implicit real*8 (a-h,o-z)
c
C Arguments:
C INPUT:
c  n1      = type of the first function (1=S, 2=P,3=L,..)
c  n2      = type of the second function
c  a       = exponent of the 1st Gaussian
c  b       = expon. of the 2nd Gaussian
c  xa(3)   = center of the 1st Gaussian
c  xb(3)   = ditto for the 3nd Gaussian
c  xc(3)   = center coordinates of the nucleus
C OUTPUT:
c  hessaa  = 6 second detivatives wrt the coordiates of atom A
c  order: AxAx, AxAy,AxAz,AyAy,AyAz,AzAz
c  hessbb  = 6 second detivatives wrt the coordiates of atom B
c  hessab  = 9 mixed second detivatives wrt the coordiates of A and B
c  the latter are defined as AxBx, AxBy, AxBz, AyBx, AyBy, AyBz,
c   AzBx, AzBy, AzBz

      dimension xa(3),xb(3),xc(3),xp(3)
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
      dimension ndeg(8)
      dimension hessaa(ndim2,ndim1,6),
     *          hessbb(ndim2,ndim1,6),
     *          hessab(ndim2,ndim1,9)
c
c  xfch, yfch, zfch are the derivative "2D" (in this case 1D) integrals
c  The first is the undifferentiated value, the second is the derivative
c  wrt A, the third wrt B, the fourth is the second derivative wrt A.
c  the fifth is the mixed second deivative, the sixth is the second
c  derivative wrt B

      dimension xfch(7,7,7,6),yfch(7,7,7,6),zfch(7,7,7,6)
      data ndeg/0,1,1,2,3,4,5,6/
c
CPP
c  temparary - zero out xfch etc.
c

      ax=xa(1)
      ay=xa(2)
      az=xa(3)
      bx=xb(1)
      by=xb(2)
      bz=xb(3)
      cx=xc(1)
      cy=xc(2)
      cz=xc(3)
c
c     calculate number of roots and xrys for rys
c
      nroots=(ndeg(n1)+ndeg(n2)+4)/2
      p=a+b
      px=(a*ax+b*bx)/p
      py=(a*ay+b*by)/p
      pz=(a*az+b*bz)/p
      xp(1)=px
      xp(2)=py
      xp(3)=pz
      rpc2=(px-cx)**2+(py-cy)**2+(pz-cz)**2
      xrys=p*rpc2
c
c     calculate roots
c
      if (nroots.le.3) call rt123
      if (nroots.eq.4) call root4
      if (nroots.eq.5) call root5
      if (nroots.gt.5) call root6
c
c   calculate derivative "2D" integrals
c
      call fcaxyzh(n1,   n2,   sab,   p,   xp,
     1             xa,   xb,   xc,    a,   b,
     2             xfch, yfch, zfch)
c  Assemble the undifferentiated integrals
c      call frmcan(n1,  n2, xfch,  yfch,  zfch,
c     1            xint)
c  Assemble the first derivative integrals wrt A and B
c  Ax
c      call frmcan(n1,  n2, xfch(1,1,1,2),yfch,zfch,
c     1            xder(1,1,1))
c  Ay
c      call frmcan(n1,  n2, xfch,yfch(1,1,1,2),zfch,
c     1            xder(1,1,2))
c  Az
c      call frmcan(n1,  n2, xfch,yfch,zfch(1,1,1,2),
c     1            xder(1,1,3))
c  Bx
c      call frmcan(n1,  n2, xfch(1,1,1,3),yfch,zfch,
c     1            xder(1,1,4))
c  By
c      call frmcan(n1,  n2, xfch,yfch(1,1,1,3),zfch,
c     1            xder(1,1,5))
c  Bz
c      call frmcan(n1,  n2, xfch,yfch,zfch(1,1,1,3),
c     1            xder(1,1,6))
c  Assemble the second derivative integrals
c  AxAx
      call frmcan(n1,  n2, xfch(1,1,1,4),yfch,zfch,
     1            hessaa(1,1,1))
c  AxAy
      call frmcan(n1,  n2, xfch(1,1,1,2),yfch(1,1,1,2),zfch,
     1            hessaa(1,1,2))
c  AxAz
      call frmcan(n1,  n2, xfch(1,1,1,2),yfch,zfch(1,1,1,2),
     1            hessaa(1,1,3))
c  AyAy
      call frmcan(n1,  n2, xfch,yfch(1,1,1,4),zfch,
     1            hessaa(1,1,4))
c  AyAz
      call frmcan(n1,  n2, xfch,yfch(1,1,1,2),zfch(1,1,1,2),
     1            hessaa(1,1,5))
c  AzAz
      call frmcan(n1,  n2, xfch,yfch,zfch(1,1,1,4),
     1            hessaa(1,1,6))
c  AxBx
      call frmcan(n1,  n2, xfch(1,1,1,5),yfch,zfch,
     1            hessab(1,1,1))
c  AxBy
      call frmcan(n1,  n2, xfch(1,1,1,2),yfch(1,1,1,3),zfch,
     1            hessab(1,1,2))
c  AxBz
      call frmcan(n1,  n2, xfch(1,1,1,2),yfch,zfch(1,1,1,3),
     1            hessab(1,1,3))
c  AyBx
      call frmcan(n1,  n2, xfch(1,1,1,3),yfch(1,1,1,2),zfch,
     1            hessab(1,1,4))
c  AyBy
      call frmcan(n1,  n2, xfch,yfch(1,1,1,5),zfch,
     1            hessab(1,1,5))
c  AyBz
      call frmcan(n1,  n2, xfch,yfch(1,1,1,2),zfch(1,1,1,3),
     1            hessab(1,1,6))
c  AzBx
      call frmcan(n1,  n2, xfch(1,1,1,3),yfch,zfch(1,1,1,2),
     1            hessab(1,1,7))
c  AzBy
      call frmcan(n1,  n2, xfch,yfch(1,1,1,3),zfch(1,1,1,2),
     1            hessab(1,1,8))
c  AzBz
      call frmcan(n1,  n2, xfch,yfch,zfch(1,1,1,5),
     1            hessab(1,1,9))
c  BxBx
      call frmcan(n1,  n2, xfch(1,1,1,6),yfch,zfch,
     1            hessbb(1,1,1))
c  BxBy
      call frmcan(n1,  n2, xfch(1,1,1,3),yfch(1,1,1,3),zfch,
     1            hessbb(1,1,2))
c  BxBz
      call frmcan(n1,  n2, xfch(1,1,1,3),yfch,zfch(1,1,1,3),
     1            hessbb(1,1,3))
c  ByBy
      call frmcan(n1,  n2, xfch,yfch(1,1,1,6),zfch,
     1            hessbb(1,1,4))
c  ByBz
      call frmcan(n1,  n2, xfch,yfch(1,1,1,3),zfch(1,1,1,3),
     1            hessbb(1,1,5))
c  BzBz
      call frmcan(n1,  n2, xfch,yfch,zfch(1,1,1,6),
     1            hessbb(1,1,6))
c
      end
c
c=======================================================================
c
      subroutine fcaxyzh(n1,   n2,   sab,   p,   xp,
     1                   xa,   xb,   xc,    a,   b,
     2                   xfc,  yfc,  zfc)
      implicit real*8 (a-h,o-z)
c
c     calculation of the two-dimensional integrals, ix, iy, and iz
c     for the evaluation of the core attraction integral Hessian
c  Arguments:
c  INTENT(IN)
c  n1,n2       = the types of the two functions (S=1,P=2,L=3,D=4,..)
c  sab         = the overlap of the two s functions
c  p           = a+b where a and b are the orbital exponents
c  xp(3)       = coordinates of the of the overlap charge density
c                 P = (a*A + b*B)/(a+b)
c  A=xa(3)     = coordinates of the first center
c  B=xb(3)     = coordinates of the second center
c  C=xc(3)     = coordinates of the attractive center
c  a           = exponent of the first function
c  b           = exponent of the second function
c
c  INTENT(OUT)
c  Derivative "2D" (in reality 1D) integrals xfc, yfc, zfc
c  They are indexed as e.g. xfc(7,7,7,6) where the first index
c  is the degree of the first function, the second is the degree
c  of the second function, the third is the Rys root, and the last
c  one is the derivative (1=the undifferentiated factor,
c  second is d/dAx, the third is d/dBx, the fourth is d2/dAx**2,
c  the ffth is  d2/dAxdBx, the sixth is d2/dBx**2 )
c  Apparently, only a small part of the 3-dimensional array is used
c  xfc(1..jdeg+1,1..ideg+1,1..nroots). The maximum is 7,7,7

      PARAMETER (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0,
     1           five=5.0d0,six=6.0d0,pi=3.14159265358979323844d0)
      dimension xp(3),xa(3),xb(3),xc(3)
      dimension xfc(7,7,7,6),yfc(7,7,7,6),zfc(7,7,7,6)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      common /hermit/ h(55)
      common /wermit/ w(55)
      dimension ndege(8)
c      dimension itm4(7)
      dimension minh(10), maxh(10)
c  x1a1(i),...x1a6(i) are [(x-Ax)**1, (x-Ax)**2,... (x-Ax)**6]*w(i)
c  It looks like only part of these quantities are used (i=minh to maxh.
c  They are calculated very awkwardly
      dimension x1a1(36),y1a1(36),z1a1(36),x1b1(36),y1b1(36),z1b1(36)
      dimension x1a2(36),y1a2(36),z1a2(36),x1b2(36),y1b2(36),z1b2(36)
      dimension x1a3(36),y1a3(36),z1a3(36),x1b3(36),y1b3(36),z1b3(36)
      dimension x1a4(36),y1a4(36),z1a4(36),x1b4(36),y1b4(36),z1b4(36)
      dimension x1a5(36),y1a5(36),z1a5(36),x1b5(36),y1b5(36),z1b5(36)
      dimension x1a6(36),y1a6(36),z1a6(36),x1b6(36),y1b6(36),z1b6(36)
      dimension x1a7(36),y1a7(36),z1a7(36),x1b7(36),y1b7(36),z1b7(36)
      dimension x1a8(36),y1a8(36),z1a8(36),x1b8(36),y1b8(36),z1b8(36)
c  ndege is the maximum degree of the function+1
      data ndege/1,2,2,3,4,5,6,7/
c  minh and maxh are used in the Hermite integration over the x,y,z factors
c  For instance, if A and B are both d functions, ndege=3 and 3,
c  and the highest degree of x is 4 if there is no derivative, 5
c   with the derivative. For 5, we need at least a 3-point integration
c  We have mxh=(6+1)/2=3, mxh1=6, and the integration points will be
c  from 4 to 6
c  However, in this case other components of the d functions may have
c  lower degrees: 0 or 1, nd the program calculates all powers for these
c  lower order integration formulas
      data minh/1,2,4, 7,11,16,22,29,37,46/,
     1     maxh/1,3,6,10,15,21,28,36,45,55/
c
      px=xp(1)
      py=xp(2)
      pz=xp(3)
c
      ax=xa(1)
      ay=xa(2)
      az=xa(3)
c
      bx=xb(1)
      by=xb(2)
      bz=xb(3)
c
      cx=xc(1)
      cy=xc(2)
      cz=xc(3)
c
      ppx=p*px
      ppy=p*py
      ppz=p*pz
c
      a2=two*a
      b2=two*b
      const=two*sqrt(p/pi)*sab
c  ie1 and ie2 are the maximum degrees of the polynomial parts in the
c  two functions
      ie1=ndege(n1)
      ie2=ndege(n2)
c  quit if either ne1 or ie2 is too high. The limit is currently 5
c  (quartic, i.e. G functions) but could be extended to i functions
      if(ie1.gt.7.or.ie2.gt.7) then
        call nerror(1,'fcaxyzh',
     2  'second derivatives for j,k,. functions are not implemented',
     3   ie1,ie2)
      end if
c
c  we need to increase mxh1 by 2 because of the differentiation
      mxh=(ie1+ie2+2)/2
      mxh1=maxh(mxh)
      do 120 iroot=1,nroots
         u2=p*rysr(iroot)
         rw=rysw(iroot)*const
         ppu2=p+u2
         sqrt1=sqrt(ppu2)
         ttx=(ppx+u2*cx)/ppu2
         tty=(ppy+u2*cy)/ppu2
         ttz=(ppz+u2*cz)/ppu2
c  loop until the maximum power needed
c  the maximum power is 2 bigger than the power of the function
c  because of the differentiation
         do 20 i=1,mxh1
            hh1=h(i)/sqrt1
            x1=hh1+ttx
            y1=hh1+tty
            z1=hh1+ttz
c  the quantities x1ax,... x6ax are the powers of (x-Ax) from 1 to 6
c  times the Hermite weight w(i)
c  the quantities x1bx,... x6bx are the powers of (x-Bx) from 1 to 6
c  for 2nd derivatives, terms through quadratic are always needed
c  S functions
            x1ax=x1-ax
            y1ay=y1-ay
            z1az=z1-az
            x1a1(i)=x1ax*w(i)
            y1a1(i)=y1ay*w(i)
            z1a1(i)=z1az*w(i)
            x1a2(i)=x1a1(i)*x1ax
            y1a2(i)=y1a1(i)*y1ay
            z1a2(i)=z1a1(i)*z1az
            if (ie1.eq.1) go to 10
c  P and L
            x1a3(i)=x1a2(i)*x1ax
            y1a3(i)=y1a2(i)*y1ay
            z1a3(i)=z1a2(i)*z1az
            if (ie1.eq.2) go to 10
c  D  - needs quartic terms
            x1a4(i)=x1a3(i)*x1ax
            y1a4(i)=y1a3(i)*y1ay
            z1a4(i)=z1a3(i)*z1az
            if (ie1.eq.3) go to 10
c  F
            x1a5(i)=x1a4(i)*x1ax
            y1a5(i)=y1a4(i)*y1ay
            z1a5(i)=z1a4(i)*z1az
            if (ie1.eq.4) go to 10
c  G
            x1a6(i)=x1a5(i)*x1ax
            y1a6(i)=y1a5(i)*y1ay
            z1a6(i)=z1a5(i)*z1az
            if (ie1.eq.5) go to 10
c  H
            x1a7(i)=x1a6(i)*x1ax
            y1a7(i)=y1a6(i)*y1ay
            z1a7(i)=z1a6(i)*z1az
            if (ie1.eq.6) go to 10
c  I
            x1a8(i)=x1a7(i)*x1ax
            y1a8(i)=y1a7(i)*y1ay
            z1a8(i)=z1a7(i)*z1az
   10       continue
c  Calculate the powers of (x-Bx)
c  S
            x1bx=x1-bx
            y1by=y1-by
            z1bz=z1-bz
            x1b1(i)=x1bx
            y1b1(i)=y1by
            z1b1(i)=z1bz
            x1b2(i)=x1b1(i)*x1bx
            y1b2(i)=y1b1(i)*y1by
            z1b2(i)=z1b1(i)*z1bz
            if (ie2.eq.1) go to 20
c  P and L
            x1b3(i)=x1b2(i)*x1bx
            y1b3(i)=y1b2(i)*y1by
            z1b3(i)=z1b2(i)*z1bz
            if (ie2.eq.2) go to 20
c  D
            x1b4(i)=x1b3(i)*x1bx
            y1b4(i)=y1b3(i)*y1by
            z1b4(i)=z1b3(i)*z1bz
            if (ie2.eq.3) go to 20
c  F
            x1b5(i)=x1b4(i)*x1bx
            y1b5(i)=y1b4(i)*y1by
            z1b5(i)=z1b4(i)*z1bz
            if (ie2.eq.4) go to 20
c  G
            x1b6(i)=x1b5(i)*x1bx
            y1b6(i)=y1b5(i)*y1by
            z1b6(i)=z1b5(i)*z1bz
            if (ie2.eq.5) go to 20
c  H
            x1b7(i)=x1b6(i)*x1bx
            y1b7(i)=y1b6(i)*y1by
            z1b7(i)=z1b6(i)*z1bz
            if (ie2.eq.6) go to 20
c  G
            x1b8(i)=x1b7(i)*x1bx
            y1b8(i)=y1b7(i)*y1by
            z1b8(i)=z1b7(i)*z1bz
   20    continue
c  The above are the powers in the 2 functions: (x-Ax)**i .. (z-Bz)**j
c  j=1...8 max.
         do 110 i1=1,ie1
           do 110 i2=1,ie2
c  add 2 for 2 differentiations
            mpts=(i1+i2+2)/2
            minh1=minh(mpts)
            maxh1=maxh(mpts)
            fcx=zero
            fcy=zero
            fcz=zero
            fcxa=zero
            fcya=zero
            fcza=zero
            fcxb=zero
            fcyb=zero
            fczb=zero
            fcxaa=zero
            fcyaa=zero
            fczaa=zero
            fcxab=zero
            fcyab=zero
            fczab=zero
            fcxbb=zero
            fcybb=zero
            fczbb=zero
c  Quadrature using the Hermite factors
c  These are the actual integration points, not just the maximum
c  Add 2 because of the second derivative
            do 100 i=minh1,maxh1
               go to (50,40,30,25,23,22,21,19,18), i1+2
c  i1+2 because we have 2 higher power
c    m=6  (i) (i1=7)
   18          f1x=x1a6(i)
               f1y=y1a6(i)
               f1z=z1a6(i)
               f1xd=x1a7(i)*a2-x1a5(i)*6.0d0
               f1yd=y1a7(i)*a2-y1a5(i)*6.0d0
               f1zd=z1a7(i)*a2-z1a5(i)*6.0d0
               f1xh=(x1a8(i)*a2-13.0d0*f1x)*a2+30.0d0*x1a4(i)
               f1yh=(y1a8(i)*a2-13.0d0*f1y)*a2+30.0d0*y1a4(i)
               f1zh=(z1a8(i)*a2-13.0d0*f1z)*a2+30.0d0*z1a4(i)
               go to 60
c    m=5  (h) (i1=6)
   19          f1x=x1a5(i)
               f1y=y1a5(i)
               f1z=z1a5(i)
               f1xd=x1a6(i)*a2-x1a4(i)*five
               f1yd=y1a6(i)*a2-y1a4(i)*five
               f1zd=z1a6(i)*a2-z1a4(i)*five
               f1xh=(x1a7(i)*a2-11.0d0*f1x)*a2+20.0d0*x1a3(i)
               f1yh=(y1a7(i)*a2-11.0d0*f1y)*a2+20.0d0*y1a3(i)
               f1zh=(z1a7(i)*a2-11.0d0*f1z)*a2+20.0d0*z1a3(i)
               go to 60
c    m=4  (g) (i1=5)
   21          f1x=x1a4(i)
               f1y=y1a4(i)
               f1z=z1a4(i)
               f1xd=x1a5(i)*a2-x1a3(i)*four
               f1yd=y1a5(i)*a2-y1a3(i)*four
               f1zd=z1a5(i)*a2-z1a3(i)*four
               f1xh=(x1a6(i)*a2-9.0d0*f1x)*a2+12.0d0*x1a2(i)
               f1yh=(y1a6(i)*a2-9.0d0*f1y)*a2+12.0d0*y1a2(i)
               f1zh=(z1a6(i)*a2-9.0d0*f1z)*a2+12.0d0*z1a2(i)
               go to 60
c    m=3 (f)
   22          f1x=x1a3(i)
               f1y=y1a3(i)
               f1z=z1a3(i)
               f1xd=x1a4(i)*a2-x1a2(i)*three
               f1yd=y1a4(i)*a2-y1a2(i)*three
               f1zd=z1a4(i)*a2-z1a2(i)*three
               f1xh=(x1a5(i)*a2-7.0d0*f1x)*a2+6.0d0*x1a1(i)
               f1yh=(y1a5(i)*a2-7.0d0*f1y)*a2+6.0d0*y1a1(i)
               f1zh=(z1a5(i)*a2-7.0d0*f1z)*a2+6.0d0*z1a1(i)
               go to 60
c  m=2 (d)
   23          f1x=x1a2(i)
               f1y=y1a2(i)
               f1z=z1a2(i)
               f1xd=x1a3(i)*a2-x1a1(i)*two
               f1yd=y1a3(i)*a2-y1a1(i)*two
               f1zd=z1a3(i)*a2-z1a1(i)*two
               f1xh=(x1a4(i)*a2-five*f1x)*a2+two*w(i)
               f1yh=(y1a4(i)*a2-five*f1y)*a2+two*w(i)
               f1zh=(z1a4(i)*a2-five*f1z)*a2+two*w(i)
               go to 60
c  m=1  (p)
   25          f1x=x1a1(i)
               f1y=y1a1(i)
               f1z=z1a1(i)
               f1xd=x1a2(i)*a2-w(i)
               f1yd=y1a2(i)*a2-w(i)
               f1zd=z1a2(i)*a2-w(i)
               f1xh=(x1a3(i)*a2-three*f1x)*a2
               f1yh=(y1a3(i)*a2-three*f1y)*a2
               f1zh=(z1a3(i)*a2-three*f1z)*a2
               go to 60
c  m=0  s (i1=1)
   30          f1x=w(i)
               f1y=w(i)
               f1z=w(i)
               f1xd=x1a1(i)*a2
               f1yd=y1a1(i)*a2
               f1zd=z1a1(i)*a2
               f1xh=(x1a2(i)*a2-f1x)*a2
               f1yh=(y1a2(i)*a2-f1y)*a2
               f1zh=(z1a2(i)*a2-f1z)*a2
               go to 60
c  the next 2 labels are impossible
   40          continue
   50          continue
c
   60          go to (90,80,70,65,63,62,61,59,58), i2+2
c  n=6 i (i2=7)
   58          f1xa=f1xd*x1b6(i)
               f1ya=f1yd*y1b6(i)
               f1za=f1zd*z1b6(i)
               f1xb=(x1b7(i)*b2-x1b5(i)*six)
               f1yb=(y1b7(i)*b2-y1b5(i)*six)
               f1zb=(z1b7(i)*b2-z1b5(i)*six)
               f1xaa=f1xh*x1b6(i)
               f1yaa=f1yh*y1b6(i)
               f1zaa=f1zh*z1b6(i)
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*((x1b8(i)*b2-13.0d0*x1b6(i))*b2+30.0d0*x1b3(i))
               f1ybb=f1y*((y1b8(i)*b2-13.0d0*y1b6(i))*b2+30.0d0*y1b3(i))
               f1zbb=f1z*((z1b8(i)*b2-13.0d0*z1b6(i))*b2+30.0d0*z1b3(i))
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               f1x=f1x*x1b6(i)
               f1y=f1y*y1b6(i)
               f1z=f1z*z1b6(i)
c...
               go to 90
c  n=5 h (i2=6)
   59          f1xa=f1xd*x1b5(i)
               f1ya=f1yd*y1b5(i)
               f1za=f1zd*z1b5(i)
               f1xb=(x1b6(i)*b2-x1b4(i)*five)
               f1yb=(y1b6(i)*b2-y1b4(i)*five)
               f1zb=(z1b6(i)*b2-z1b4(i)*five)
               f1xaa=f1xh*x1b5(i)
               f1yaa=f1yh*y1b5(i)
               f1zaa=f1zh*z1b5(i)
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*((x1b7(i)*b2-11.0d0*x1b5(i))*b2+20.0d0*x1b3(i))
               f1ybb=f1y*((y1b7(i)*b2-11.0d0*y1b5(i))*b2+20.0d0*y1b3(i))
               f1zbb=f1z*((z1b7(i)*b2-11.0d0*z1b5(i))*b2+20.0d0*z1b3(i))
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               f1x=f1x*x1b5(i)
               f1y=f1y*y1b5(i)
               f1z=f1z*z1b5(i)
c...
               go to 90
c  n=4 g (i2=5)
   61          f1xa=f1xd*x1b4(i)
               f1ya=f1yd*y1b4(i)
               f1za=f1zd*z1b4(i)
               f1xb=(x1b5(i)*b2-x1b3(i)*four)
               f1yb=(y1b5(i)*b2-y1b3(i)*four)
               f1zb=(z1b5(i)*b2-z1b3(i)*four)
               f1xaa=f1xh*x1b4(i)
               f1yaa=f1yh*y1b4(i)
               f1zaa=f1zh*z1b4(i)
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*((x1b6(i)*b2-9.0d0*x1b4(i))*b2+12.0d0*x1b2(i))
               f1ybb=f1y*((y1b6(i)*b2-9.0d0*y1b4(i))*b2+12.0d0*y1b2(i))
               f1zbb=f1z*((z1b6(i)*b2-9.0d0*z1b4(i))*b2+12.0d0*z1b2(i))
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               f1x=f1x*x1b4(i)
               f1y=f1y*y1b4(i)
               f1z=f1z*z1b4(i)
c...
               go to 90
c  n=3  f
   62          f1xa=f1xd*x1b3(i)
               f1ya=f1yd*y1b3(i)
               f1za=f1zd*z1b3(i)
               f1xb=(x1b4(i)*b2-x1b2(i)*three)
               f1yb=(y1b4(i)*b2-y1b2(i)*three)
               f1zb=(z1b4(i)*b2-z1b2(i)*three)
               f1xaa=f1xh*x1b3(i)
               f1yaa=f1yh*y1b3(i)
               f1zaa=f1zh*z1b3(i)
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*((x1b5(i)*b2-7.0d0*x1b3(i))*b2+6.0d0*x1b1(i))
               f1ybb=f1y*((y1b5(i)*b2-7.0d0*y1b3(i))*b2+6.0d0*y1b1(i))
               f1zbb=f1z*((z1b5(i)*b2-7.0d0*z1b3(i))*b2+6.0d0*z1b1(i))
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               f1x=f1x*x1b3(i)
               f1y=f1y*y1b3(i)
               f1z=f1z*z1b3(i)
               go to 90
c  n=2  d
   63          f1xa=f1xd*x1b2(i)
               f1ya=f1yd*y1b2(i)
               f1za=f1zd*z1b2(i)
               f1xb=(x1b3(i)*b2-x1b1(i)*two)
               f1yb=(y1b3(i)*b2-y1b1(i)*two)
               f1zb=(z1b3(i)*b2-z1b1(i)*two)
               f1xaa=f1xh*x1b2(i)
               f1yaa=f1yh*y1b2(i)
               f1zaa=f1zh*z1b2(i)
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*((x1b4(i)*b2-five*x1b2(i))*b2+two)
               f1ybb=f1y*((y1b4(i)*b2-five*y1b2(i))*b2+two)
               f1zbb=f1z*((z1b4(i)*b2-five*z1b2(i))*b2+two)
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               f1x=f1x*x1b2(i)
               f1y=f1y*y1b2(i)
               f1z=f1z*z1b2(i)
               go to 90
c  n=1 p
   65          f1xa=f1xd*x1b1(i)
               f1ya=f1yd*y1b1(i)
               f1za=f1zd*z1b1(i)
               f1xb=x1b2(i)*b2-one
               f1yb=y1b2(i)*b2-one
               f1zb=z1b2(i)*b2-one
               f1xaa=f1xh*x1b1(i)
               f1yaa=f1yh*y1b1(i)
               f1zaa=f1zh*z1b1(i)
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*((x1b3(i)*b2-three*x1b1(i))*b2)
               f1ybb=f1y*((y1b3(i)*b2-three*y1b1(i))*b2)
               f1zbb=f1z*((z1b3(i)*b2-three*z1b1(i))*b2)
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               f1x=f1x*x1b1(i)
               f1y=f1y*y1b1(i)
               f1z=f1z*z1b1(i)
               go to 90
c  n=0  s  (i2=1)
   70          f1xa=f1xd
               f1ya=f1yd
               f1za=f1zd
               f1xb=x1b1(i)*b2
               f1yb=y1b1(i)*b2
               f1zb=z1b1(i)*b2
               f1xaa=f1xh
               f1yaa=f1yh
               f1zaa=f1zh
               f1xab=f1xd*f1xb
               f1yab=f1yd*f1yb
               f1zab=f1zd*f1zb
               f1xbb=f1x*(x1b2(i)*b2-one)*b2
               f1ybb=f1y*(y1b2(i)*b2-one)*b2
               f1zbb=f1z*(z1b2(i)*b2-one)*b2
c...
               f1xb=f1x*f1xb
               f1yb=f1y*f1yb
               f1zb=f1z*f1zb
               go to 90
c  The next 2 cases are impossible with i2>=1
   80          continue
   90          fcx=fcx+f1x
               fcy=fcy+f1y
               fcz=fcz+f1z
               fcxa=fcxa+f1xa
               fcya=fcya+f1ya
               fcza=fcza+f1za
               fcxb=fcxb+f1xb
               fcyb=fcyb+f1yb
               fczb=fczb+f1zb
               fcxaa=fcxaa+f1xaa
               fcyaa=fcyaa+f1yaa
               fczaa=fczaa+f1zaa
               fcxab=fcxab+f1xab
               fcyab=fcyab+f1yab
               fczab=fczab+f1zab
               fcxbb=fcxbb+f1xbb
               fcybb=fcybb+f1ybb
               fczbb=fczbb+f1zbb
  100       continue
             xfc(i2,i1,iroot,1)=fcx
             yfc(i2,i1,iroot,1)=fcy
             zfc(i2,i1,iroot,1)=fcz*rw
             xfc(i2,i1,iroot,2)=fcxa
             yfc(i2,i1,iroot,2)=fcya
             zfc(i2,i1,iroot,2)=fcza*rw
             xfc(i2,i1,iroot,3)=fcxb
             yfc(i2,i1,iroot,3)=fcyb
             zfc(i2,i1,iroot,3)=fczb*rw
             xfc(i2,i1,iroot,4)=fcxaa
             yfc(i2,i1,iroot,4)=fcyaa
             zfc(i2,i1,iroot,4)=fczaa*rw
             xfc(i2,i1,iroot,5)=fcxab
             yfc(i2,i1,iroot,5)=fcyab
             zfc(i2,i1,iroot,5)=fczab*rw
             xfc(i2,i1,iroot,6)=fcxbb
             yfc(i2,i1,iroot,6)=fcybb
             zfc(i2,i1,iroot,6)=fczbb*rw
  110    continue
  120 continue
      return
c
      end
c======================================================================
      subroutine frmcan(n1,  n2,  xfc,  yfc, zfc,  xint)
      implicit real*8 (a-h,o-z)
      PARAMETER (zero=0.0d0)
c     the array xfc,yfc,zfc should be (7,7,7), the indices
c     being xfc(ityp1,ityp2,iroot).
c      common /xyzfc/ xfc(343),yfc(343),zfc(343)
      dimension xfc(7,7,7),yfc(7,7,7),zfc(7,7,7)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      dimension nbeg(8), nend(8), nxtyp(84), nytyp(84), nztyp(84)
      dimension xint(784)
c      dimension itm4(7), xint(784)
      data nbeg/1,2,1,5,11,21,36,57/,nend/1,4,4,10,20,35,56,84/
      data nxtyp/1, 2,1,1, 3,1,1,2,2,1, 4,3,3,2,2,2,1,1,1,1,
     1           5,4,4,3,3,3,2,2,2,2,1,1,1,1,1,
     2           6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,
     3           7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,2,2,2,2,2,2,
     4           1,1,1,1,1,1,1/
      data nytyp/1, 1,2,1, 1,3,1,2,1,2, 1,2,1,3,2,1,4,3,2,1,
     1           1,2,1,3,2,1,4,3,2,1,5,4,3,2,1,
     2           1,2,1,3,2,1,4,3,2,1,5,4,3,2,1,6,5,4,3,2,1,
     3           1,2,1,3,2,1,4,3,2,1,5,4,3,2,1,6,5,4,3,2,1,
     4           7,6,5,4,3,2,1/
      data nztyp/1,1,1,2, 1,1,3,1,2,2, 1,1,2,1,2,3,1,2,3,4,
     1           1,1,2,1,2,3,1,2,3,4,1,2,3,4,5,
     2           1,1,2,1,2,3,1,2,3,4,1,2,3,4,5,1,2,3,4,5,6,
     3           1,1,2,1,2,3,1,2,3,4,1,2,3,4,5,1,2,3,4,5,6,
     4           1,2,3,4,5,6,7/
c     data itm4/0,4,8,12,16/
c      data itm4/0,7,14,21,28,35,42/
c
      ib1=nbeg(n1)
      ib2=nbeg(n2)
      ie1=nend(n1)
      ie2=nend(n2)
      ncount=0
      do 40 ityp1=ib1,ie1
         lax=nxtyp(ityp1)
         lay=nytyp(ityp1)
         laz=nztyp(ityp1)
c          indx1=itm4(lax)
c          indy1=itm4(lay)
c          indz1=itm4(laz)
      do 40 ityp2=ib2,ie2
c          indx=indx1+nxtyp(ityp2)
c          indy=indy1+nytyp(ityp2)
c          indz=indz1+nztyp(ityp2)
         lbx=nxtyp(ityp2)
         lby=nytyp(ityp2)
         lbz=nztyp(ityp2)
         sum=zero
c         go to (30,20,10,5,4,3,2), nroots
c    2    sum=sum+xfc(indx+294)*yfc(indy+294)*zfc(indz+294)
c    3    sum=sum+xfc(indx+245)*yfc(indy+245)*zfc(indz+245)
c    4    sum=sum+xfc(indx+196)*yfc(indy+196)*zfc(indz+196)
c    5    sum=sum+xfc(indx+147)*yfc(indy+147)*zfc(indz+147)
c   10    sum=sum+xfc(indx+98)*yfc(indy+98)*zfc(indz+98)
c   20    sum=sum+xfc(indx+49)*yfc(indy+49)*zfc(indz+49)
c   30    sum=sum+xfc(indx)*yfc(indy)*zfc(indz)
       do iroot=1,nroots
        sum=sum+xfc(lbx,lax,iroot)*yfc(lby,lay,iroot)*zfc(lbz,laz,iroot)
       end do
         ncount=ncount+1
         xint(ncount)=sum
   40 continue
      return
c
      end
c
c=======================================================================
      SubRoutine IntDDer(nAtoms,inx,datbas,datnuc,nCS,nCF,nTri,
     *                   Dens,DMDer)
c
      Implicit Real*8 (a-h,o-z)
      Common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c-----------------------------------------------------------------------
c
c This routine calculates 1st derivatives of the dipole integrals and
c with 0-th order density matrix forms the DipolMomentDerivative matrix
c needed for IR intensities
c
c This is a modification of the intof subroutine
c
c-----------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c Dens(ntri)      - zeroth order density matrix
c INTENT(OUT)
c DMDer(3,natoms,3) - Dipole moment derivatives
c-----------------------------------------------------------------------
      Dimension datbas(13,*),datnuc(5,*)
      Dimension inx(12,*)
      Dimension Dens(nTri)
      dimension DMDer(3,nAtoms,3)
c---- Local variables
      Dimension PLDer(1)                   !   only for WrDeriv call
      Dimension fff(3*784,3)
      Dimension f3a(3*784,3),f3b(3*784,3)
*
      Character*256 JobName
      Common /job/JobName,lenJ
c-----------------------------------------------------------------------
c
c---- Set the dipole moment integrals derivatives matrix to zero
      Call DCopy(3*3*nAtoms,0.d0,0,DMDer,1)
c
      k2=0
      key=1
      ifu=0
c
      Do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         Do 55 igc=0,ngci
            jfu=0
            Do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               If (j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               Do 45 jgc=0,ngcjx
c
                  Do k1=1,3
                     Call DCopy(len3,zero,0,f3a(1,k1),1)
                     Call DCopy(len3,zero,0,f3b(1,k1),1)
                     Call DCopy(len3,zero,0,fff(1,k1),1)
c
                     Call onef(i,j,igc,jgc,datbas,f3a(1,k1),
     &                         inx,key,k1,k2)
                     Call onef(j,i,jgc,igc,datbas,fff(1,k1),
     &                         inx,key,k1,k2)
                  End Do
c
                  Do k1=1,3
                     Do icr=1,3
                        icrl=(icr-1)*len
                        Do j3=1,len2
                           Do i3=1,len1
                              ji=(j3-1)*len1+i3
                              ij=(i3-1)*len2+j3
                              f3b(icrl+ij,k1)=fff(icrl+ji,k1)
                           End Do
                        End Do
                     End Do
                  End Do
c
                  iij=0
                  iff=ifu
                  Do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                     Do 40 j1=1,len2
                        jff=jff+1
                        jj=jff*(jff-1)/2
                        iij=iij+1
                        if (jff.gt.iff) go to 40
                        ij=ii+jff
c
                        dij=dens(ij)
                        if(iff.ne.jff) dij=2.d0*dij
c
                        Do icr=1,3
                           iij_icr=iij+(icr-1)*len
                           Do k1=1,3
c
                              dmder(icr,iat,k1)=dmder(icr,iat,k1)
     *                                         -f3a(iij_icr,k1)*dij
                              dmder(icr,jat,k1)=dmder(icr,jat,k1)
     *                                         -f3b(iij_icr,k1)*dij
c
                           End Do
                        End Do
   40             Continue
                  jfu=jfu+len2
   45          Continue
   50       Continue
            ifu=ifu+len1
   55    Continue
   60 Continue
c
c---- Save this part of the dipole moment derivatives on disc
      Open(Unit=40,File=jobname(1:lenJ)//'.deriv_temp',
     $     Form='formatted',Status='unknown')
      Do iAt=1,nAtoms
         Do i=1,3
            Write(40,910) DMDer(i,iAt,1),DMDer(i,iAt,2),DMDer(i,iAt,3)
         End Do
      End Do
      Close(Unit=40,Status='keep')
  910 Format(10X,3F20.14)
*     Call WrDeriv(3*nAtoms,DMDer,PLDer,.true.,.false.)
c
      End
c=======================================================================
      subroutine get_fromd(iunit,ntri3,natb,nate,f1)
      implicit real*8 (a-h,o-z)
      dimension f1(ntri3,natb:nate)
c
         do iat=natb,nate
            call read1mat(iunit,iat,ntri3,f1(1,iat))
         enddo
c
      end
c======================================================================
      subroutine trsp_fromd(iunit,ntri,natb,nate,f1,bl,listunq)
      implicit real*8 (a-h,o-z)
      dimension listunq(*)
      dimension bl(ntri,3)
      dimension f1(3,natb:nate,ntri)
c
      do iat=natb,nate
         iatr=listunq(iat)   !   real atom
cccccc   call read1mat(iunit,iatr,ntri*3,bl)
         read(unit=iunit,rec=iatr) bl
         do ij=1,ntri
            f1(1,iat,ij)=bl(ij,1)
            f1(2,iat,ij)=bl(ij,2)
            f1(3,iat,ij)=bl(ij,3)
         enddo
      enddo
c
      end
c======================================================================
      subroutine put_ondsk(iunit,irec,ntri3,iatr,f1)
      implicit real*8 (a-h,o-z)
      dimension f1(ntri3)
c
      call save1mat(iunit,irec+iatr,ntri3,f1)
c
      end
c======================================================================
      subroutine printest1(text,r1,ntri,iatu,iatr)
      implicit real*8 (a-h,o-z)
      character text*4
      dimension r1(ntri,3)
c
         sux=0.d0
         suy=0.d0
         suz=0.d0
         do ij=1,ntri
            sux=sux+r1(ij,1)
            suy=suy+r1(ij,2)
            suz=suz+r1(ij,3)
         enddo
         sum=sux+suy+suz
         write(6,66) text,iatu, sum, iatr
  66  format(a4,'unq at=',i3,' sx+sy+sz=',1(f12.6,1x),' real at=',i3)
      end
c======================================================================
      subroutine printest(text,iunit,r1,natonce,ntri)
      implicit real*8 (a-h,o-z)
      character text*4
      dimension r1(ntri,3)
c
      do iat=1,natonce
         call read1mat(iunit,iat,ntri*3,r1(1,1))
         sux=0.d0
         suy=0.d0
         suz=0.d0
         do ij=1,ntri
            sux=sux+r1(ij,1)
            suy=suy+r1(ij,2)
            suz=suz+r1(ij,3)
         enddo
         sum=sux+suy+suz
         write(6,66) text,iat, sum
      enddo
  66  format(a4,'  atom=',i3,' sx+sy+sz=',1(f12.6,1x))
      end
c======================================================================
      subroutine printest2(text,r1,ncs)
      implicit real*8 (a-h,o-z)
      character text*4
      dimension r1(ncs,ncs)
         sum=0.d0
         do i=1,ncs
         do j=1,ncs
            sum=sum+r1(i,j)
         enddo
         enddo
         write(6,66) text, sum
  66  format(a4,' sum =',f12.6,1x)
      end
c======================================================================
