C APR 1 98 PP  I have removed the thresholds on the overlap etc.
c  integrals. They are insignificant compared to the nuclear
c  attraction force. I ngelect the nucl. attraction force by
c  a different algorithm: S0*abs(zzz)*max(sqa,sqb), where
c  zza is the nucl. charge and sqa,sqb are square roots of the
c exponents. Note that one should use a nad b but this is
c consistent with the energy expression.
c===================================================================
c             One-electron part of the gradient
c
c                    July 97, K.Wolinski
c
c===================================================================
      subroutine force1(rhf,bl,inx,iprint)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      Logical rhf
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c external field info :
      common /fieldx/ xfield,xfgrad,elfiel(9)
c Gradient options :
c
      common /forcdbl/ thre1,thre2,tchf
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------------------------
c
c COSMO stuff
c
      character*100 cerm

C For CIM calculation. -NZG
      logical cim
      character*4 chopv
      character filefockder*256
      character scrfile*256
      real*8,allocatable::temp(:,:)
c
      call secund(tgrad1)
      call elapsec(elagrad1)
c------------------------------------------------------------------
      ifield=xfield
      ifgrad=xfgrad
c------------------------------------------------------------------
c Retrive memory addresses for alpha density (lden), beta density
c if any (mden), 1-electron forces (lforce1) and 2-el. forces (lforc2)
c
      call getival('ldensi',lden  )
      call getival('lforc1',lforc1)
      call getival('lforc3',lforc3)
      If(.not.rhf) call getival('ldensB',mden)
      call getival('nslv',nslv)
c------------------------------------------------------------------
c Allocate memory & calculate weighted density matrix :
c
      call mem_wdens(rhf,bl,ldew)
c
c------------------------------------------------------------------
c check-run only or lack of scf convergence :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c
c -- add alpha & beta densities if UHF
      If(.not.rhf) Call AddVEC(ntri,bl(lden),bl(mden),bl(lden))
c------------------------------------------------------------------
c TEST
c     call dipolex(bl,inx,ncf,ncs,na)
c------------------------------------------------------------------
c Calculate nuclear repulsion part of the gradient
c
      call nucrep_grad(na,bl(inuc),bl(lforc3),ifield,ifgrad,elfiel)
      if(iprint.gt.1) then
        call forcepr(na,bl(lforc3),'nucrep')
      end if
c------------------------------------------------------------------

c Check if CIM is the next step. If CIM force is to be calculated,
c the kinetic and electron-nuclear attraction contributions will be
c added up and stored.
c NZG_10/11/2017 @UARK
      cim=.false.
      read(inp,'(A4)',end=95) chopv
      call lowercas(chopv,4)
      if (chopv=='cim ') cim=.true.
      backspace inp
95    continue

      if (cim) then
         call getmem(na*3*ntri,lfock1)
         call setival('lfock1',lfock1)
         call zeroit(bl(lfock1),na*3*ntri)

         call getint(na,listreal)
         call setup_listreal(na,bl(listreal))

C For simple programming, I directly use the way of calculating the
C integrals in MP2 force. In this case, they will be calculated twice.
         call intonFder(1,na,inx,ifield,ifgrad,elfiel,bl(ibas),bl(inuc),
     &                  ncs,ncf,ntri,bl(listreal),bl(lfock1))
         call elenuc_Fder(na,1,na,inx,bl(ibas),bl(inuc),ncs,ncf,ntri,
     &                    bl(listreal),bl(lfock1))
C         call fder_print(bl(lfock1),na,ntri,1.d0)
         call getchval('scrf',scrfile)
         call rmblan(scrfile,256,len)
         filefockder=scrfile(1:len)//'.fock1'
         lrec=3*ntri*8
         nfockfile=723
         open(unit=nfockfile,file=filefockder,form='unformatted',
     &        access='direct',recl=lrec)
         allocate(temp(ntri,3))
         call trsp_disk(nfockfile,na,ntri,temp,bl(lfock1))
         close(unit=nfockfile,status='keep')
         deallocate(temp)
      endif

c Calculate S1 & T1 (overlap & kinetic) derivative integrals
c
      call intof(na,inx,ifield,ifgrad,elfiel,
     *           bl(ibas),bl(inuc),ncs,ncf,ntri,
     *           bl(lden),bl(ldew),bl(lforc1) )
c
      if(iprint.gt.1) then
        call forcepr(na,bl(lforc1),'s1+t1')
      end if
c------------------------------------------------------------------
c Calculate electron-nuclear attraction contributions to the forces
c
c ** This can now be done in parallel **
c
      If(nslv.eq.0) Then
        call elapsec(t1)
        call elenuc_grad(na,inx,bl(ibas),bl(inuc),ncs,ncf,ntri,
     *                   bl(lden),bl(lforc1) )
        call elapsec(t2)
cc        write(6,*) ' Time for 1-e e-nuclear forces is:',t2-t1,' sec'
c
        if(iprint.gt.1) then
          call forcepr(na,bl(lforc1),'ele-nuc')
        end if
      EndIf
c
c  pseudopotentials contribution to forces
c
      call getival('npsp',npsp)
      if(npsp.gt.0)then
cc        call secund(pspt0)
cc        call elapsec(pspet0)
        psptol=0.1d0*thre1
        call pspgrad(na,npsp,ncf,ncs,inx,bl(inuc),bl(ibas),
     $               bl(lden),psptol,bl(lforc1))
        if(iprint.gt.1) then
          call forcepr(na,bl(lforc1),'Pseudopotential contribution')
        end if
cc        call secund(pspt1)
cc        call elapsec(pspet1)
cc        psptime=pspt1-pspt0
cc        pspelap=pspet1-pspet0
      endif
c
c    check if the COSMO flag has been defined, otherwise
c    define it now.
c
      call tstival('cosmo',icosmo)
      If(icosmo.EQ.0) Then
        call setival('cosmo',0)
      Else
        call getival('cosmo',icosmo)
      EndIf
c
      tcosmo=0.0d0
      ecosmo=0.0d0
      if(icosmo.ne.0)then
c
c   COSMO contribution to forces:
c
        call secund(ct0)
        call elapsec(celt0)
c
c   obtain COSMO parameters from depository
c
        call getival('ndum',ndum)       ! no. of dummy atoms
        natom = na-ndum                 ! no. of real atoms
        call getival('c_nps',nps)       ! no. of surface segments
        call getival('c_nip',nip)       ! no. of intersection pairs
        call getival('c_lcavit',lcavity)! type of cavity
        call getival('nslv',nslv)
c
c   allocate memory for COSMO forces
c
        call mmark
c       call immark
        call getmem(3*na,icfor)
        call getmem(3*natom,icoxyz)
        call getmem(3*nps,icosurf)
        call getmem(nps,iar)
        call getmem(nps,iqcos)
        call getmem(4*nip,isude)
        call getint(nps,iiatsp)
        call getint(4*nip,iisude)
        call getint(natom,icharge)
c
c  read arrays from COSMO data file
c
        call cosmo_rdata(bl(icoxyz),bl(icosurf),bl(iar),bl(iqcos),
     $                   bl(iiatsp),bl(icharge),natom,nps)
c
c  read surface derivatives information
c
        ierr=0
        call reasude(nip,bl(isude),bl(iisude),ierr,cerm)
cc
        If(ierr.NE.0) Then
          write(iout,'(a)') cerm
          call nerror(7,'FORCE Module',
     $                   'Error in COSMO Forces',0,0)
        EndIf
c
c  contributions to forces from surface derivatives and
c  nuclear-surface potential.
c  cavnucgrad actually computes gradient contributions,
c  thus scale by -1.0 to obtain forces.
c
        call zeroit(bl(icfor),3*na)
        call cavnucgrd(natom,bl(icoxyz),nps,bl(icosurf),bl(iiatsp),
     $                 bl(iqcos),bl(iar),nip,bl(isude),bl(iisude),
     $                 lcavity,bl(icfor),bl(icharge))
        call VScal(3*na,-1.0d0,bl(icfor))
        if(iprint.gt.1) then
          call forcepr(na,bl(icfor),
     $                'COSMO surface + nuclear (unscaled)')
        end if
c
c  contribution to forces from electron-surface potential
c
        if(nslv.eq.0)then
        call cosmo_pot_forc(bl(lden),bl(icfor),bl(iqcos),bl(icosurf),
     $                      inx,bl(ibas),bl(iiatsp),nps,ncs,ncf,ntri)
        else
        call para_cosmo_forc(bl(lden),bl(icfor),bl(iqcos),bl(icosurf),
     $                       bl(iiatsp),na,nps,ntri)
        endif
        if(iprint.gt.1) then
          call forcepr(na,bl(icfor),'COSMO  + potential (unscaled)')
        end if
c
c  scale COSMO force contribution by fepsi
c
        call getrval('c_fepsi',fepsi)
        call VScal(3*na,fepsi,bl(icfor))
        if(iprint.gt.1) then
          call forcepr(na,bl(icfor),'COSMO total (scaled)')
        end if
c
c  add COSMO contribution to one electron forces
c
        call addvec(3*na,bl(lforc1),bl(icfor),bl(lforc1))
c
c  end of COSMO forces
c
        call secund(ct1)
        call elapsec(celt1)
        tcosmo=tcosmo+ct1-ct0
        ecosmo=ecosmo+celt1-celt0
      write(iout,
     $'(''Master CPU time for COSMO      gradient =''
     $,f8.2,'' Elapsed = '',f8.2,'' min'')') tcosmo/60.0d0,
     $ ecosmo/60.0d0
c       call retimark
        call retmark
      endif
c
c -- restore separate densities if UHF
      If(.not.rhf) Call MinusVEC(ntri,bl(lden),bl(mden),bl(lden))
c------------------------------------------------------------------
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) (total-tcosmo)/60.0d0,(elaps-ecosmo)/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for 1e part of gradient ='
     *,f8.2,' Elapsed = ',f8.2,' min')
c------------------------------------------------------------------
      end
c===================================================================
      subroutine forcepr(na,force,label)
      implicit real*8 (a-h,o-z)
      dimension force(3,na)
      character*(*) label
c  prints diagnostic forces
      iout=igetival('iout')
      write(iout,*) label
      do i=1,na
        write(iout,'(i4,2x,3f13.7)') i,(force(k,i),k=1,3)
      end do
      end
c===================================================================
      subroutine mem_wdens(rhf,bl,ldew)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c The purpose of this routine is to allocate memory for, and to construct
c the weighted density matrix (or matrices for UHF).
c  ARGUMENTS:
c  rhf     = true if restricted RHF, false for UHF
c  bl      = common area to allocate memory
c  ldew    = the address of the weighted density matrix
      Logical rhf
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension bl(*)
c------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c------------------------------------------------------------------
c check-run only or lack of scf convergance :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c
c------------------------------------------------------------------
c allocate memory for the weighted density
c
      call getmem(ntri,ldew)
c
c allocate temporary memory for eig_vec & eig_val
c
      call getmem(ncf*ncf,lvec)
      call getmem(ncf    ,lval)
      call getmem(2      ,lnoc)
c------------------------------------------------------------------
c read-in the eigenvectors & eigenvalues :
c
      np4=4
c
      If(rhf) Then
c -- read in closed-shell MOs and orbital energies
       call rea (bl(lvec),ncf*ncf,np4,'evec_rhf')
       call rea (bl(lval),ncf    ,np4,'eval_rhf')
       call rea (bl(lnoc),1      ,np4,'nocc_rhf')
       nocc=bl(lnoc)
c -- construct the weighted density matrix :
       call make_wdens(ncf,nocc,bl(lvec),bl(lval),bl(ldew),rhf)
      Else
c -- read in alpha MOs and orbital energies
       call rea (bl(lvec),ncf*ncf,np4,'evea_uhf')
       call rea (bl(lval),ncf    ,np4,'evaa_uhf')
       call rea (bl(lnoc),2      ,np4,'nocc_uhf')
       NAlpha = bl(lnoc)
       NBeta  = bl(lnoc+1)
c -- construct the alpha energy weighted density matrix
       call make_wdens(ncf,NAlpha,bl(lvec),bl(lval),bl(ldew),rhf)
c -- now read in beta MOs and orbital energies
       call rea (bl(lvec),ncf*ncf,np4,'eveb_uhf')
       call rea (bl(lval),ncf    ,np4,'evab_uhf')
c -- construct the total weighted density matrix
       call make_wdens(ncf,-NBeta,bl(lvec),bl(lval),bl(ldew),rhf)
      EndIf

c------------------------------------------------------------------
c release memory
c
      call retmem(3)  ! lvec, lval & lnoc
c------------------------------------------------------------------
      end
c======================================================================
      subroutine make_wdens(ncf,nocc,vec,val,dew,rhf)
      implicit real*8 (a-h,o-z)
c
c  construct energy-weighted density matrix
c  modified for UHF  April 98   (JB)
c
c  ARGUMENTS
c
c  ncf     -  number of contracted basis functions
c  nocc    -  number of occupied MOs
c             ** WARNING **  may be NEGATIVE for UHF indicating
c             that matrix is to be ADDED to existing (alpha) matrix
c  vec     -  MOs (closed-shell or alpha/beta)
c  val     -  corresponding orbital energies
c  dew     -  on exit contains energy-weighted density matrix
c  rhf     -  logical flag for closed-shell
c
      dimension vec(ncf,ncf),val(ncf)
      dimension dew(*)
      Logical rhf
c
      If(nocc.GT.0) Then
c -- alpha or closed-shell
        ij=0
        do 60 i=1,ncf
        do 50 j=1,i
        s=0.0d0
        ij=ij+1
         do 40 k=1,nocc
         s=s+vec(i,k)*vec(j,k)*val(k)
   40    continue
        dew(ij)=s
   50   continue
   60   continue
        If(rhf) call VScal((ncf*(ncf+1))/2,2.0d0,dew)
cc
      Else
c -- beta part (needs to be ADDED to existing matrix)
        ij=0
        do 61 i=1,ncf
        do 51 j=1,i
        s=0.0d0
        ij=ij+1
         do 41 k=1,-nocc
         s=s+vec(i,k)*vec(j,k)*val(k)
   41    continue
        dew(ij)=dew(ij)+s
   51   continue
   61   continue
      EndIf
c
      end
c======================================================================
      subroutine nucrep_grad(natoms,datnuc,atforces,ifield,ifgrad,xfld)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c This routine calculates nuclear repulsion terms of the gradient
c---------------------------------------------------------------------
c INPUT :
c natoms           - number of atoms
c datnuc(5,natoms) - nuclear date
c xfld(9)          - components of an external field (3)&field grad.(6)
c atforces(3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datnuc(5,*)
      dimension atforces(3,natoms)
      dimension xfld(9)
c---------------------------------------------------------------------
      do nra=1,natoms
         zza=datnuc(1,nra)
         xra=datnuc(2,nra)
         yra=datnuc(3,nra)
         zra=datnuc(4,nra)
c derivative of nuclear dipole interacting with an external electric field
         if(ifield.eq.1) then
            sx=-xfld(1)
            sy=-xfld(2)
            sz=-xfld(3)
         else
            sx=zero
            sy=zero
            sz=zero
         endif
         do nat=1,natoms
c ......... nuclear repulsion term
c
            if(nat.ne.nra) then
               zzb=datnuc(1,nat)
               xrb=datnuc(2,nat)
               yrb=datnuc(3,nat)
               zrb=datnuc(4,nat)
c
               xba=xrb-xra
               yba=yrb-yra
               zba=zrb-zra
c
               r=xba*xba+yba*yba+zba*zba
               r=sqrt(r)
               r=r*r*r
               sx=sx + zzb*xba/r
               sy=sy + zzb*yba/r
               sz=sz + zzb*zba/r
            endif
         enddo      ! nat=1,natoms
c
c....... accumulate nuclear repulsion forces
c
         atforces(1,nra)=atforces(1,nra)-zza*sx
         atforces(2,nra)=atforces(2,nra)-zza*sy
         atforces(3,nra)=atforces(3,nra)-zza*sz
c
      enddo         ! nra=1,natoms
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine elenuc_grad(natoms,inx,datbas,datnuc,ncs,ncf,ntri,
     *                       dn,atfor)
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension dn(ntri)
      dimension atfor(3,natoms)
c
      do nra=1,natoms
         call intof7(nra,inx,datbas,datnuc,ncs,ncf,ntri,dn,atfor)
      enddo
c
      end
c======================================================================
      subroutine intof7(nra,inx,datbas,datnuc,ncs,ncf,ntri,dn,atfor)
      implicit real*8 (a-h,o-z)
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
c dn(ntri)        - density matrix
c atfor(3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension dn(ntri)            ! density
      dimension atfor(3,*)          ! atfor(3,natoms)
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
                     dij=dn(ij)
                     if(iff.ne.jff) dij=2.d0*dij
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
                  atfor(icr,nra)=atfor(icr,nra)
     *                           -(sa(iij_icr)+sb(iij_icr))*dij*zza
                  atfor(icr,iat)=atfor(icr,iat)+sa(iij_icr)*dij*zza
                  atfor(icr,jat)=atfor(icr,jat)+sb(iij_icr)*dij*zza
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
      subroutine onef(ics,  jcs,  igc,  jgc,  datbas,
     1                sb,   inx,  key,  k1,   k2)
      implicit real*8 (a-h,o-z)
c this subroutine calculates the one-electron integral derivatives for a
c contracted pair of shells, and one particular general contraction
c (cf. onel which includes the general contraction loop)
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
c are in a separate routine intof7; 9=mixed GIAO and nuclear derivatives
c 10-ONE-SIDED nuclear derivative overlaps, 11=ONE-SIDED magnetic deriv.
c INTENT(OUT)
c sb holds the integral derivatives
c format is sb(jlen1,ilen1,3)
c  where ilen1 is the size of contracted shell ics (3 for P, 5 for D)
c  jlen1 is the size of contracted shell jcs,
c  and the last subscript denotes the derivatives (x,y,z)
c  for mixed GIAO-nuclear derivatives, the format is sb(jlen,ilen,3.3)
c  with the last index being the magnetic direction and the 3rd index
c  the nuclear direction OF THE FIRST CONTRACTION ONLY
c
      common /forcdbl/ thre1,thre2,tchf
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension inx(12,*)
      dimension datbas(13,*)
      dimension sb(*)
clocal:
      dimension rb(9*784)
      dimension xa(3),xb(3)
      data twopi/0.6366197723675d0/  !  this is 1/sqrt(12)
c----------------------------------------------------------------------
      ldim=3
      if(key.eq.9) ldim=9
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
      if(ityp.eq.11) ityp1=6
      if(ityp.eq.12) ityp1=7
      if(ityp.eq.13) ityp1=8
      if(jtyp.eq.11) jtyp1=6
      if(jtyp.eq.12) jtyp1=7
      if(jtyp.eq.13) jtyp1=8

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
c key.eq.1.or.key.eq.2 or key=9 or key=10 or 11 only
c
            call inconf(ityp1,jtyp1,key,k1,k2,a,b,s0,xa,xb,rb)
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
c transformation d6->d5, f10->f7   g15->g9
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
        end if
        if(ityp.eq.12) then
           call htran1a(rb,sb(iadd),jlen1)
        end if
        if(ityp.eq.13) then
c           call itran1a(rb,sb(iadd),jlen1)
        end if
c
c at this point i-functions are transformed so the length is ilen
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
      subroutine inconf(ityp,jtyp,key,k1,k2,a,b,s0,xa,xb,xint)
      implicit real*8 (a-h,o-z)
c    this routine constructs a set of one-electron integral derivatives
c    over two primitive shells
c   ARGUMENTS:
c
c  INTENT(IN):
c  ityp,jtyp: the types of the shells, s=1,p=2,l=3,d=4,f=5,g=6,h=7,i=8
c    note that only 6-component d and 10-component f are used here
c    transformation to d5 and f7 is later
c  key is the type of the integral: 1=overlap, 2=kinetic, 3=a dipole
c   component, 4=r**2,, 5=a second moment component
c  Currently only 1,2 and 9 are permitted
c  1= overlap, 2=kinetic energy, 9=mixed GIAO and nuclear derivative
c  overlap, 10=one-sided nuclear derivative overlap, 11=one-sided
c  GIAO magnetic derivative overlap
c  k1 and k2 are Cartesian components. They can range from 1=x to 3=z
c   they will be used for the dipole and 2nd moment integrals
c  a and b are the Gaussian exponents, sqa and sqb their square roots
c  s0 is the overlap integral of the spherical parts of the Gaussians
c  xa(3) and xb(3) are the center coordinates of the Gaussians
c  INTENT(OUT):
c  xint(1:lenj,1:leni,3) gives the integral derivatives.
c   leni and lenj are the sizes (lengths) of the shells, e.g.
c   3 for P, 4 for L, 6 for D6
c
c  For the mixed magnetic and nuclear derivatives xint is defined
c  (implicitly) as xint(1:lenj,1:leni,1:3,1:3); xint(j,i,l,k)
c  corresponds to element (i,j) differentiated by the k component of
c  the magnetic field and the l component of the nuclear coordinate of 
c  the atom corresponding to i, i.e. xa(l,i) if xa is (3,natoms)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      dimension xint(*)
      dimension xa(3),xb(3),len(8),la(3),lb(3),ll(3)
      dimension mulx(88),muly(88),mulz(88),nfu(9)
      dimension yint(3*784),zint2(784*9)
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
      do i=1,ln*3
        xint(i)=zero
        yint(i)=zero
      end do
      do 20 i=1,3
         la(i)=0
         lb(i)=0
         ll(i)=0
   20 continue
   25 go to(30,40,1000,1000,1000,1000,1000,1000,300,400,500) key
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
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      go to 1000
c
   40 call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 50 i=ia,ie
      do 50 j=ja,je
         ij=ij+1
         ft=dble(2*(mulx(j)+muly(j)+mulz(j))+3)*b
         xint(ij)=yint(ij)*ft
         xint(ij+ln)=yint(ij+ln)*ft
         xint(ij+2*ln)=yint(ij+2*ln)*ft
   50 continue
      lb(1)=2
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ft=-two*b**2
      do 60 i=1,ln*3
   60 xint(i)=xint(i)+yint(i)*ft
      lb(1)=0
      lb(2)=2
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      do 70 i=1,ln*3
   70 xint(i)=xint(i)+yint(i)*ft
      lb(2)=0
      lb(3)=2
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      lb(3)=0
      do 80 i=1,ln*3
   80 xint(i)=xint(i)+yint(i)*ft
      if (jtyp.lt.4) go to 1000
      lb(3)=0
      lb(1)=-2
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 90 i=ia,ie
        do 90 j=ja,je
         ij=ij+1
         ft=dble((mulx(j)*(mulx(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
         xint(ij+ln)=xint(ij+ln)-yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
   90 continue
      lb(1)=0
      lb(2)=-2
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 100 i=ia,ie
      do 100 j=ja,je
         ij=ij+1
         ft=dble((muly(j)*(muly(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
         xint(ij+ln)=xint(ij+ln)-yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
  100 continue
      lb(2)=0
      lb(3)=-2
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      lb(3)=0
      ij=0
      do 110 i=ia,ie
      do 110 j=ja,je
         ij=ij+1
         ft=dble((mulz(j)*(mulz(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
         xint(ij+ln)=xint(ij+ln)-yint(ij+ln)*ft
         xint(ij+2*ln)=xint(ij+2*ln)-yint(ij+2*ln)*ft
  110 continue
      go to 1000
  300 continue
c  Calculation of the mixed GIAO and nuclear derivatives
      do k=1,3
        ll(k)=1
        istart2=(k-1)*ln*3+1
c  integrals <d(mu)/dR|(x,y or z)(nu)>
        call intarf(ityp,  jtyp, a, b, s0,
     1              xa,    xb,   la,lb,ll,
     2              zint2(istart2))
c
        ll(k)=0
      end do
      do k=1,3          !magnetic field direction
        kk=ln*3*(k-1)
        kk1=k+1
        if(kk1.gt.3) kk1=kk1-3
        kk2=k+2
        if(kk2.gt.3) kk2=kk2-3
        da1=xb(kk1)
        da2=xb(kk2)
        do l=1,3        ! force direction
          kl=kk+ln*(l-1)
          iz21=((kk2-1)*3+l-1)*ln
          iz22=((kk1-1)*3+l-1)*ln
          do i=1,ln
c -- commented out due to apparent compiler bug   ! MM  April 2008
cc            call z2elem(ln,zint2,i,l,kk1,kk2,z21,z22)
            z21=zint2(iz21+i) !  z21=zint2(i,l,kk2)
            z22=zint2(iz22+i) !  z22=zint2(i,l,kk1)
            xint(i+kl)=(da1*z21-da2*z22)*half
          end do
        end do
      end do
      go to 1000
c
c  Left-sided nuclear derivatives
  400 continue
      call intarf(ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
CPP
      go to 1000
c  Right-side magnetic derivatives
  500 continue
      do k=1,3
        ll(k)=1
        istart1=(k-1)*ln+1
c  integrals <(mu)|(x,y or z)(nu)>, needed for the magnetic derivatives
        call intar(ityp,  jtyp, a, b, s0,
     1              xa,    xb,   la,lb,ll,
     2              zint2(istart1))
c  Construction of the GIAO derivatives from the dipole integrals
        ll(k)=0
      end do
      do k=1,3
        istart1=(k-1)*ln
        k1=k+1
        if(k1.gt.3) k1=k1-3
        k1start=(k1-1)*ln
        k2=k1+1
        if(k2.gt.3) k2=k2-3
        k2start=(k2-1)*ln
        do i=1,ln
          xint(istart1+i)=(xb(k1)*zint2(k2start+i)
     1                   -xb(k2)*zint2(k1start+i))*half 
        end do
      end do
      go to 1000
c
 1000 continue
c
      end
c========================================================================
      subroutine z2elem(ln,zint2,i,l,k1,k2,z21,z22)
      implicit real*8 (a-h,o-z)
c  This routine retrives elements z21=zint2(i,l,k2) and
c  zint22=zint2(i,l,k1)
      dimension zint2(ln,3,3)
      z21=zint2(i,l,k2)
      z22=zint2(i,l,k1)
c
      end
c=======================================================================
      subroutine intarf(ityp,  jtyp,  a,  b,  s0,
     1                  xa,    xb,    la, lb, l3,
     2                  xint)
      implicit real*8 (a-h,o-z)
c
c this routine calculates one-electron integral derivatives over primitive
c gaussian shells. It is assumed that the 1s (spherical) part of the gaussian
c is normalized (this is denoted by 'gaussian').
c the result is
c d/dxa(k){Int(x-xa(1))**la(1)*...(z-xb(3))**lb(3)*x**l3(1)..
c  * z**l3(3)*(x-xa(1))**mulx1*...(z-zb(3))**mulz2*gaussian1*gaussian2}
c  i.e. it can include, besides the x,y,z factors associated with the
c  angular momentum functions, extra (x-xa),...(z-zb) factors, as well
c  as extra x,y,z factors
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
c the results is built in (((xint(j,i,k),j=1,jsh),i=1,ish),k=1,3)
c where ish,jsh are the shell sizes: 1,3,4,6,10,15,21,28 for
c s,p,l,d,f,g,h,i
c it should be ok up to and including i functions, although it would
c be a good idea to test it for high ang. mom. functions befor using it
c it uses hermite-gaussian quadrature. note that this is not very
c economical but transparent.
c
c  Argument list:
c input:
c ityp,jtyp: function types, 1,2,3..8 for s,p,l,d,f,g,h,i
c a,b: Gaussian exponents
c s0: overlap between normalized 1s gaussians
c xa(3),xb(3): orbital centers
c la(3),lb(3),l3(3) : cartesian exponents, see above.
c output:
c  xint(jsh,ish,3) or xint(jsh,ish,3,3) where jsh=1,jlen
c  arguments; the last one for the mixed nuclear and magnetic
c  derivative overlap integrals
      dimension xa(3),xb(3),la(3),lb(3),l3(3),xint(*)
c xint should be dimensioned to 3*784 for i functions and 9*784
c  if mixed derivatives are calculated
c  data
      dimension idg(8),nfu(9),mul(3,88)
c  local variables
      dimension xfc(81),yfc(81),zfc(81)
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
      l1x=iadg+la(1)+1   ! the highest possible power of (x-Ax)
      l2x=ibdg+lb(1)     ! the highest possible power of (x-Bx)
      l1y=iadg+la(2)+1   ! the highest possible power of (y-Ay)
      l2y=ibdg+lb(2)     ! the highest possible power of (y-By)
      l1z=iadg+la(3)+1   ! the highest possible power of (z-Az)
      l2z=ibdg+lb(3)     ! the highest possible power of (z-Bz)
      istart=nfu(ityp)+1
      ilen=nfu(ityp+1)-istart+1
      jstart=nfu(jtyp)+1
      jlen=nfu(jtyp+1)-jstart+1
c
c     if(ityp.ge.7 .and. jtyp.ge.7) then
c     write(6,*)' ityp,jtyp=',ityp,jtyp,' l1x,l2x=',l1x,l2x
c     call f_lush(6)
c     endif
c
      call inta2f(ilen,  jlen,   mul(1,istart), mul(1,jstart), la,
     1            lb,    l3,     a,             b,             xa,
     3            xb,    l1x+1,  l1y+1,         l1z+1,         l2x,
     4            l2y,   l2z,    s0,            xfc,           yfc,
     5            zfc,   xint)
      end
c========================================================================
      subroutine inta2f(ilen,  jlen,   mul1,  mul2,   la,
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
      dimension xa(3),xb(3),xint(jlen,ilen,3),la(3),lb(3),l3(3)
      dimension mul1(3,ilen),mul2(3,jlen)
c  local variables
      dimension xfc(0:l2x,0:l1x),yfc(0:l2y,0:l1y),zfc(0:l2z,0:l1z)
      dimension ixf(3),iixf(3)
c
      PARAMETER (Zero=0.0d0,one=1.0d0,two=2.0d0)
c  call the routine xtor
c   in order to be correct for i functions and a fourth-degree
c   operator (say, hexadecapole), it needs the hermitian roots
c   and weights up to 9th degree; we have it up to the 10th,
c    so it could be expanded to j functions...
      a2=two*a
      call xtor (l1x,l2x,l3(1),a,b,xa(1),xb(1),xfc)
      call xtor (l1y,l2y,l3(2),a,b,xa(2),xb(2),yfc)
      call xtor (l1z,l2z,l3(3),a,b,xa(3),xb(3),zfc)
      do 60 k=1,3                   ! differentiation wrt ax(k)
        do 40 i=1,ilen
          ixf(1)=mul1(1,i)+la(1)
          ixf(2)=mul1(2,i)+la(2)
          ixf(3)=mul1(3,i)+la(3)
          iixf(1)=ixf(1)
          iixf(2)=ixf(2)
          iixf(3)=ixf(3)
          iixf(k)=ixf(k)+1
          xmult=ixf(k)
          ixf(k)=ixf(k)-1
          do 20 j=1,jlen
            jxf=mul2(1,j)+lb(1)
            jyf=mul2(2,j)+lb(2)
            jzf=mul2(3,j)+lb(3)
c  guard against an exponent being negative
            if (jxf.ge.0.and.jyf.ge.0.and.jzf.ge.0) then
            xint(j,i,k)=
     1      a2*xfc(jxf,iixf(1))*yfc(jyf,iixf(2))*zfc(jzf,iixf(3))
              if(ixf(1).ge.0.and.ixf(2).ge.0.and.ixf(3).ge.0) then
                xint(j,i,k)=xint(j,i,k) -
     2           xmult*xfc(jxf,ixf(1))*yfc(jyf,ixf(2))*zfc(jzf,ixf(3))
              end if
            end if
            xint(j,i,k)=xint(j,i,k)*s0
   20     continue
c  restore ixf(k) and iixf(k)
          ixf(k)=ixf(k)+1
          iixf(k)=iixf(k)-1
   40   continue
   60 continue
      end
c========================================================================
      subroutine onef7(ics,jcs,igc,jgc,datbas,sa,sb,inx,xra,zza)
C  Arguments:
C  INPUT:
c  ics,jcs: contracted shell indices
c  igc,jgc: depth of general contaction (0 is unsegmented)
c  datbas: basis set data
c  inx: contraction data
c  xra(3): coordinates of the nucleus
c  zza: charge of the nucleus
C OUTPUT:
c  sa(3*maxshell), sb(3*maxshell): derivatives of the 1-el. integrals
      implicit real*8 (a-h,o-z)
      common /forcdbl/ thre1,thre2,tchf
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension inx(12,*)
      dimension datbas(13,*)
      dimension xra(3),sa(*),sb(*)
clocal:
      dimension ra(3*784),rb(3*784)
      dimension xa(3),xb(3)
      data twopi/0.6366197723675d0/  !  this is 1/sqrt(12)
c----------------------------------------------------------------------
      ldim=3
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
c
      if(ityp.eq.11) ityp1=6
      if(ityp.eq.12) ityp1=7
      if(ityp.eq.13) ityp1=8
      if(jtyp.eq.11) jtyp1=6
      if(jtyp.eq.12) jtyp1=7
      if(jtyp.eq.13) jtyp1=8
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
      sa(i)=zero
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
            estim=s0*abs(zza)*max(sqa,sqb)
            if(estim.lt.thre1 ) go to 500
c
            call inconf7(ityp1,jtyp1,a,b,s0,xa,xb,ra,rb,xra)
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
                 coefij=coefi*coefj
                 sa(ij)=sa(ij)+ra(ij)*coefij
                 sb(ij)=sb(ij)+rb(ij)*coefij
  400         continue
  410       continue
c
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
c           call itran1a(ra,sa(iadd),jlen1)
c           call itran1a(rb,sb(iadd),jlen1)
        end if
c
c at this point i-functions are transformed so the length is ilen
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
      sa(ij)=sa(iadd+ij1)
      sb(ij)=sb(iadd+ij1)
  730 continue
  720 continue
c----------------------------------------------------------------------
      end
c====================================================================
      subroutine inconf7(ityp,jtyp,a,b,s0,xa,xb,yint,wint,xra)
      implicit real*8 (a-h,o-z)
c
      dimension yint(*),wint(*)
      dimension xa(3),xb(3),len(8)
      dimension mulx(88),muly(88),mulz(88),nfu(9)
clocal
      dimension xint(3*784)
      dimension xra(3)
c
      data zero,two /0.d0,2.d0/
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
      twoa=two*a
      twob=two*b
c
      ln=len(ityp)*len(jtyp)
c
      ia=nfu(ityp)+1
      ie=nfu(ityp+1)
      ja=nfu(jtyp)+1
      je=nfu(jtyp+1)
c
      do 10 i=1,ln*3
         xint(i)=zero
         yint(i)=zero
         wint(i)=zero
   10 continue
c
      kij=0
      call nucpof(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),xb(3),
clalb*            xra(1),xra(2),xra(3),la,lb,xint,wint,kij,ln)
     *            xra(1),xra(2),xra(3),      xint,wint,kij,ln)
c
      ij=0
      do 450 i=ia,ie
         ftx=-dble(mulx(i))
         fty=-dble(muly(i))
         ftz=-dble(mulz(i))
      do 450 j=ja,je
         ij=ij+1
         yint(ij+2*ln)=xint(ij+2*ln)*ftz + wint(ij+2*ln)*twoa
         yint(ij+ln)=xint(ij+ln)*fty +wint(ij+ln)*twoa
         yint(ij)=xint(ij)*ftx+wint(ij)*twoa
  450 continue
c
      kij=1
      call nucpof(ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),xb(3),
clalb*            xra(1),xra(2),xra(3),la,lb,xint,wint,kij,ln)
     *            xra(1),xra(2),xra(3),      xint,wint,kij,ln)
c
      ij=0
      do 451 i=ia,ie
      do 451 j=ja,je
         ftx=-dble(mulx(j))
         fty=-dble(muly(j))
         ftz=-dble(mulz(j))
         ij=ij+1
         wint(ij+2*ln)=xint(ij+2*ln)*ftz + wint(ij+2*ln)*twob
         wint(ij+ln)=xint(ij+ln)*fty + wint(ij+ln)*twob
         wint(ij)=xint(ij)*ftx + wint(ij)*twob
  451 continue
c
      end
c====================================================================
      subroutine nucpof(n1,n2,a,b,sab,ax,ay,az,bx,by,bz,cx,cy,cz,
     *                  xm1,xp1,kij,ln)
      implicit real*8 (a-h,o-z)
c this routine apparently calculates the components for the derivatives
c  of the nuclear integral with respect to the coordinates of the first atom
c
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
      dimension ndeg(8)
      dimension xm1(*),xp1(*)     ! dimension 3*ln each
c
      data ndeg/0,1,1,2,3,4,5,6/
c
c     calculate number of roots and xrys for rys
c
      nroots=(ndeg(n1)+ndeg(n2)+3)/2
c
      p=a+b
      px=(a*ax+b*bx)/p
      py=(a*ay+b*by)/p
      pz=(a*az+b*bz)/p
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
c calculate two dimensional derivative integrals (fcfxyz)
c and assemble core attraction integral derivatives (frmcf)
c
      do icr=1,3
         iadd=(icr-1)*ln
         call fcfxyz(n1,n2,sab,p,px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz,
     *               icr,kij)
         call frmcf(N1,N2,xm1(iadd+1),xp1(iadd+1))
      enddo
c
      end
c====================================================================
      subroutine fcfxyz(n1,n2,sab,p,px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz,
     *                  kk,kij)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c calculation of the two-dimensional integrals, ix, iy, and iz
c for the evaluation of the core attraction derivative integrals
c---------------------------------------------------------------------
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /xyzfc/ xfcm(343),yfcm(343),zfcm(343),
     *               xfcp(343),yfcp(343),zfcp(343)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      common /hermit/ h(55)
      common /wermit/ w(55)
      dimension itm4(7), ndege(8)
      dimension minh(8), maxh(8)
c
      dimension xyza(3,36),xyzb(3,36)
c
      dimension xa1m(3,36),xa1p(3,36),xb1(3,36)
      dimension xa2m(3,36),xa2p(3,36),xb2(3,36)
      dimension xa3m(3,36),xa3p(3,36),xb3(3,36)
      dimension xa4m(3,36),xa4p(3,36),xb4(3,36)
      dimension xa5m(3,36),xa5p(3,36),xb5(3,36)
      dimension xa6m(3,36),xa6p(3,36),xb6(3,36)
      dimension xa7m(3,36),xa7p(3,36),xb7(3,36)
c
      data itm4/0,7,14,21,28,35,42/
c     this is a relic and simulates a 3-index array xfc(i2,i1,iroot)
      data ndege/1,2,2,3,4,5,6,7/
      data minh/1,2,4,7,11,16,22,29/,maxh/1,3,6,10,15,21,28,36/
c
      ppx=p*px
      ppy=p*py
      ppz=p*pz
      const=two*sqrt(p/pi)*sab
      ie1=ndege(n1)
      ie2=ndege(n2)
      mxh1=(ie1+ie2+1)/2
c     if(mxh1.eq.0) mxh1=1
      mxh1=maxh(mxh1)
c
c     if(n1.ge.7.and.n2.ge.7) then
c     write(6,*)'n1,n2=',n1,n2,' ie1,ie2=',ie1,ie2,' Kij=',kij
c     endif
c
      IF(KIJ.eq.1) go to 4444
c
c
      indz=0
      do 120 iroot=1,nroots
         u2=p*rysr(iroot)
         rw=rysw(iroot)*const
         ppu2=p+u2
         sqrt1=sqrt(ppu2)
         ttx=(ppx+u2*cx)/ppu2
         tty=(ppy+u2*cy)/ppu2
         ttz=(ppz+u2*cz)/ppu2
         do 20 i=1,mxh1
            hh1=h(i)/sqrt1
            x1=hh1+ttx
            y1=hh1+tty
            z1=hh1+ttz
c
            xyza(1,i)=x1-ax
            xyza(2,i)=y1-ay
            xyza(3,i)=z1-az
c
            do 200 ki=1,3
            xa1m(ki,i)=one
            xa1p(ki,i)=one
  200       continue
            xa1m(kk,i)=zero
            xa1p(kk,i)=xyza(kk,i)
            if(ie1.eq.1) go to 10
c
            do 201 ki=1,3
            xa2m(ki,i)=xyza(ki,i)
            xa2p(ki,i)=xyza(ki,i)
  201       continue
            xa2m(kk,i)=one
            xa2p(kk,i)=xyza(kk,i)*xa2p(kk,i)
            if(ie1.eq.2) go to 10
c
            do 202 ki=1,3
            xa3m(ki,i)=xyza(ki,i)*xa2m(ki,i)
            xa3p(ki,i)=xyza(ki,i)*xa2p(ki,i)
  202       continue
            if(ie1.eq.3) go to 10
c
            do 203 ki=1,3
            xa4m(ki,i)=xyza(ki,i)*xa3m(ki,i)
            xa4p(ki,i)=xyza(ki,i)*xa3p(ki,i)
  203       continue
            if(ie1.eq.4) go to 10
c
            do 204 ki=1,3
            xa5m(ki,i)=xyza(ki,i)*xa4m(ki,i)
            xa5p(ki,i)=xyza(ki,i)*xa4p(ki,i)
  204       continue
            if(ie1.eq.5) go to 10
c
            do 205 ki=1,3
            xa6m(ki,i)=xyza(ki,i)*xa5m(ki,i)
            xa6p(ki,i)=xyza(ki,i)*xa5p(ki,i)
  205       continue
            if(ie1.eq.6) go to 10
c
            do 206 ki=1,3
            xa7m(ki,i)=xyza(ki,i)*xa6m(ki,i)
            xa7p(ki,i)=xyza(ki,i)*xa6p(ki,i)
  206       continue
c---
   10       continue
c---
c
            xyzb(1,i)=x1-bx
            xyzb(2,i)=y1-by
            xyzb(3,i)=z1-bz
c
            do 306 ki=1,3
            xb1(ki,i)=one
  306       continue
            if (ie2.eq.1) go to 20
c
            do 307 ki=1,3
            xb2(ki,i)=xyzb(ki,i)
  307       continue
            if (ie2.eq.2) go to 20
c
            do 308 ki=1,3
            xb3(ki,i)=xyzb(ki,i)*xb2(ki,i)
  308       continue
            if (ie2.eq.3) go to 20
c
            do 309 ki=1,3
            xb4(ki,i)=xyzb(ki,i)*xb3(ki,i)
  309       continue
            if (ie2.eq.4) go to 20
c
            do 310 ki=1,3
            xb5(ki,i)=xyzb(ki,i)*xb4(ki,i)
  310       continue
            if (ie2.eq.5) go to 20
c
            do 311 ki=1,3
            xb6(ki,i)=xyzb(ki,i)*xb5(ki,i)
  311       continue
            if (ie2.eq.6) go to 20
c
            do 312 ki=1,3
            xb7(ki,i)=xyzb(ki,i)*xb6(ki,i)
  312       continue
c----
   20    continue
c----
c
         do 110 i1=1,ie1
            ind1=itm4(i1)
         do 110 i2=1,ie2
            ind=ind1+i2
            mpts=(i1+i2+1)/2
            minh1=minh(mpts)
            maxh1=maxh(mpts)
            fcxm=zero
            fcym=zero
            fczm=zero
            fcxp=zero
            fcyp=zero
            fczp=zero
            do 100 i=minh1,maxh1
c
cccccc         go to (26,25,24,23,22,21), i1
               go to (27,26,25,24,23,22,21), i1
c
   21          f1xm=xa7m(1,i)
               f1ym=xa7m(2,i)
               f1zm=xa7m(3,i)
               f1xp=xa7p(1,i)
               f1yp=xa7p(2,i)
               f1zp=xa7p(3,i)
               go to 60
c
   22          f1xm=xa6m(1,i)
               f1ym=xa6m(2,i)
               f1zm=xa6m(3,i)
               f1xp=xa6p(1,i)
               f1yp=xa6p(2,i)
               f1zp=xa6p(3,i)
               go to 60
c
   23          f1xm=xa5m(1,i)
               f1ym=xa5m(2,i)
               f1zm=xa5m(3,i)
               f1xp=xa5p(1,i)
               f1yp=xa5p(2,i)
               f1zp=xa5p(3,i)
               go to 60
c
   24          f1xm=xa4m(1,i)
               f1ym=xa4m(2,i)
               f1zm=xa4m(3,i)
               f1xp=xa4p(1,i)
               f1yp=xa4p(2,i)
               f1zp=xa4p(3,i)
               go to 60

   25          f1xm=xa3m(1,i)
               f1ym=xa3m(2,i)
               f1zm=xa3m(3,i)
               f1xp=xa3p(1,i)
               f1yp=xa3p(2,i)
               f1zp=xa3p(3,i)
               go to 60
c
   26          f1xm=xa2m(1,i)
               f1ym=xa2m(2,i)
               f1zm=xa2m(3,i)
               f1xp=xa2p(1,i)
               f1yp=xa2p(2,i)
               f1zp=xa2p(3,i)
               go to 60
c
   27          f1xm=xa1m(1,i)
               f1ym=xa1m(2,i)
               f1zm=xa1m(3,i)
               f1xp=xa1p(1,i)
               f1yp=xa1p(2,i)
               f1zp=xa1p(3,i)
               go to 60

c----
   60          continue
c----
ccc            go to (80,70,65,63,62,61), i2
               go to (67,66,65,64,63,62,61), i2
c
   61          f1xm=f1xm*xb7(1,i)
               f1ym=f1ym*xb7(2,i)
               f1zm=f1zm*xb7(3,i)
               f1xp=f1xp*xb7(1,i)
               f1yp=f1yp*xb7(2,i)
               f1zp=f1zp*xb7(3,i)
               go to 90
c
   62          f1xm=f1xm*xb6(1,i)
               f1ym=f1ym*xb6(2,i)
               f1zm=f1zm*xb6(3,i)
               f1xp=f1xp*xb6(1,i)
               f1yp=f1yp*xb6(2,i)
               f1zp=f1zp*xb6(3,i)
               go to 90
c
   63          f1xm=f1xm*xb5(1,i)
               f1ym=f1ym*xb5(2,i)
               f1zm=f1zm*xb5(3,i)
               f1xp=f1xp*xb5(1,i)
               f1yp=f1yp*xb5(2,i)
               f1zp=f1zp*xb5(3,i)
               go to 90
c
   64          f1xm=f1xm*xb4(1,i)
               f1ym=f1ym*xb4(2,i)
               f1zm=f1zm*xb4(3,i)
               f1xp=f1xp*xb4(1,i)
               f1yp=f1yp*xb4(2,i)
               f1zp=f1zp*xb4(3,i)
               go to 90
c
   65          f1xm=f1xm*xb3(1,i)
               f1ym=f1ym*xb3(2,i)
               f1zm=f1zm*xb3(3,i)
               f1xp=f1xp*xb3(1,i)
               f1yp=f1yp*xb3(2,i)
               f1zp=f1zp*xb3(3,i)
               go to 90
c
   66          f1xm=f1xm*xb2(1,i)
               f1ym=f1ym*xb2(2,i)
               f1zm=f1zm*xb2(3,i)
               f1xp=f1xp*xb2(1,i)
               f1yp=f1yp*xb2(2,i)
               f1zp=f1zp*xb2(3,i)
               go to 90
c
   67          f1xm=f1xm*xb1(1,i)
               f1ym=f1ym*xb1(2,i)
               f1zm=f1zm*xb1(3,i)
               f1xp=f1xp*xb1(1,i)
               f1yp=f1yp*xb1(2,i)
               f1zp=f1zp*xb1(3,i)
c
   90          continue
c
               fcxm=fcxm+f1xm*w(i)
               fcym=fcym+f1ym*w(i)
               fczm=fczm+f1zm*w(i)
               fcxp=fcxp+f1xp*w(i)
               fcyp=fcyp+f1yp*w(i)
               fczp=fczp+f1zp*w(i)
  100       continue
            xfcm(ind+indz)=fcxm
            yfcm(ind+indz)=fcym
            zfcm(ind+indz)=fczm*rw
            xfcp(ind+indz)=fcxp
            yfcp(ind+indz)=fcyp
            zfcp(ind+indz)=fczp*rw
  110    continue
  120 indz=indz+49
c
      return
c
 4444 CONTINUE
c
      indz=0
      do 5120 iroot=1,nroots
         u2=p*rysr(iroot)
         rw=rysw(iroot)*const
         ppu2=p+u2
         sqrt1=sqrt(ppu2)
         ttx=(ppx+u2*cx)/ppu2
         tty=(ppy+u2*cy)/ppu2
         ttz=(ppz+u2*cz)/ppu2
         do 520 i=1,mxh1
            hh1=h(i)/sqrt1
            x1=hh1+ttx
            y1=hh1+tty
            z1=hh1+ttz
c
            xyza(1,i)=x1-ax
            xyza(2,i)=y1-ay
            xyza(3,i)=z1-az
c
            do 5200 ki=1,3
            xb1(ki,i)=one
 5200       continue
            if(ie1.eq.1) go to 510
c
            do 5201 ki=1,3
            xb2(ki,i)=xyza(ki,i)
 5201       continue
            if(ie1.eq.2) go to 510
c
            do 5202 ki=1,3
            xb3(ki,i)=xyza(ki,i)*xb2(ki,i)
 5202       continue
            if(ie1.eq.3) go to 510
c
            do 5203 ki=1,3
            xb4(ki,i)=xyza(ki,i)*xb3(ki,i)
 5203       continue
            if(ie1.eq.4) go to 510
c
            do 5204 ki=1,3
            xb5(ki,i)=xyza(ki,i)*xb4(ki,i)
 5204       continue
            if(ie1.eq.5) go to 510
c
            do 5205 ki=1,3
            xb6(ki,i)=xyza(ki,i)*xb5(ki,i)
 5205       continue
            if(ie1.eq.6) go to 510
c
            do 5206 ki=1,3
            xb7(ki,i)=xyza(ki,i)*xb6(ki,i)
 5206       continue
c---
  510       continue
c---
c
            xyzb(1,i)=x1-bx
            xyzb(2,i)=y1-by
            xyzb(3,i)=z1-bz
c
            do 6206 ki=1,3
            xa1m(ki,i)=one
            xa1p(ki,i)=one
 6206       continue
            xa1m(kk,i)=zero
            xa1p(kk,i)=xyzb(kk,i)
            if(ie2.eq.1) go to 520
c
            do 6207 ki=1,3
            xa2m(ki,i)=xyzb(ki,i)
            xa2p(ki,i)=xyzb(ki,i)
 6207       continue
            xa2m(kk,i)=one
            xa2p(kk,i)=xyzb(kk,i)*xyzb(kk,i)
            if(ie2.eq.2) go to 520
c
            do 6208 ki=1,3
            xa3m(ki,i)=xyzb(ki,i)*xa2m(ki,i)
            xa3p(ki,i)=xyzb(ki,i)*xa2p(ki,i)
 6208       continue
            if(ie2.eq.3) go to 520
c
            do 6209 ki=1,3
            xa4m(ki,i)=xyzb(ki,i)*xa3m(ki,i)
            xa4p(ki,i)=xyzb(ki,i)*xa3p(ki,i)
 6209       continue
            if(ie2.eq.4) go to 520
c
            do 6210 ki=1,3
            xa5m(ki,i)=xyzb(ki,i)*xa4m(ki,i)
            xa5p(ki,i)=xyzb(ki,i)*xa4p(ki,i)
 6210       continue
            if(ie2.eq.5) go to 520
c
            do 6211 ki=1,3
            xa6m(ki,i)=xyzb(ki,i)*xa5m(ki,i)
            xa6p(ki,i)=xyzb(ki,i)*xa5p(ki,i)
 6211       continue
            if(ie2.eq.6) go to 520
c
            do 6212 ki=1,3
            xa7m(ki,i)=xyzb(ki,i)*xa6m(ki,i)
            xa7p(ki,i)=xyzb(ki,i)*xa6p(ki,i)
 6212       continue
c----
  520    continue
c----
c
         do 5110 i1=1,ie1
            ind1=itm4(i1)
         do 5110 i2=1,ie2
            ind=ind1+i2
            mpts=(i1+i2+1)/2
            minh1=minh(mpts)
            maxh1=maxh(mpts)
            fcxm=zero
            fcym=zero
            fczm=zero
            fcxp=zero
            fcyp=zero
            fczp=zero
            do 5100 i=minh1,maxh1
c
ccccc          go to (540,530,525,523,522,521), i1
               go to (527,526,525,524,523,522,521), i1
c
  521          f1xm=xb7(1,i)
               f1ym=xb7(2,i)
               f1zm=xb7(3,i)
               f1xp=xb7(1,i)
               f1yp=xb7(2,i)
               f1zp=xb7(3,i)
               go to 560
c
  522          f1xm=xb6(1,i)
               f1ym=xb6(2,i)
               f1zm=xb6(3,i)
               f1xp=xb6(1,i)
               f1yp=xb6(2,i)
               f1zp=xb6(3,i)
               go to 560
c
  523          f1xm=xb5(1,i)
               f1ym=xb5(2,i)
               f1zm=xb5(3,i)
               f1xp=xb5(1,i)
               f1yp=xb5(2,i)
               f1zp=xb5(3,i)
               go to 560
c
  524          f1xm=xb4(1,i)
               f1ym=xb4(2,i)
               f1zm=xb4(3,i)
               f1xp=xb4(1,i)
               f1yp=xb4(2,i)
               f1zp=xb4(3,i)
               go to 560

  525          f1xm=xb3(1,i)
               f1ym=xb3(2,i)
               f1zm=xb3(3,i)
               f1xp=xb3(1,i)
               f1yp=xb3(2,i)
               f1zp=xb3(3,i)
               go to 560
c
  526          f1xm=xb2(1,i)
               f1ym=xb2(2,i)
               f1zm=xb2(3,i)
               f1xp=xb2(1,i)
               f1yp=xb2(2,i)
               f1zp=xb2(3,i)
               go to 560
c
  527          f1xm=xb1(1,i)
               f1ym=xb1(2,i)
               f1zm=xb1(3,i)
               f1xp=xb1(1,i)
               f1yp=xb1(2,i)
               f1zp=xb1(3,i)
               go to 560

c----
  560          continue
c----
ccccccc        go to (580,570,565,563,562,561), i2
               go to (567,566,565,564,563,562,561), i2
c
  561          f1xm=f1xm*xa7m(1,i)
               f1ym=f1ym*xa7m(2,i)
               f1zm=f1zm*xa7m(3,i)
               f1xp=f1xp*xa7p(1,i)
               f1yp=f1yp*xa7p(2,i)
               f1zp=f1zp*xa7p(3,i)
               go to 590
c
  562          f1xm=f1xm*xa6m(1,i)
               f1ym=f1ym*xa6m(2,i)
               f1zm=f1zm*xa6m(3,i)
               f1xp=f1xp*xa6p(1,i)
               f1yp=f1yp*xa6p(2,i)
               f1zp=f1zp*xa6p(3,i)
               go to 590
c
  563          f1xm=f1xm*xa5m(1,i)
               f1ym=f1ym*xa5m(2,i)
               f1zm=f1zm*xa5m(3,i)
               f1xp=f1xp*xa5p(1,i)
               f1yp=f1yp*xa5p(2,i)
               f1zp=f1zp*xa5p(3,i)
               go to 590
c
  564          f1xm=f1xm*xa4m(1,i)
               f1ym=f1ym*xa4m(2,i)
               f1zm=f1zm*xa4m(3,i)
               f1xp=f1xp*xa4p(1,i)
               f1yp=f1yp*xa4p(2,i)
               f1zp=f1zp*xa4p(3,i)
               go to 590
c
  565          f1xm=f1xm*xa3m(1,i)
               f1ym=f1ym*xa3m(2,i)
               f1zm=f1zm*xa3m(3,i)
               f1xp=f1xp*xa3p(1,i)
               f1yp=f1yp*xa3p(2,i)
               f1zp=f1zp*xa3p(3,i)
               go to 590
c
  566          f1xm=f1xm*xa2m(1,i)
               f1ym=f1ym*xa2m(2,i)
               f1zm=f1zm*xa2m(3,i)
               f1xp=f1xp*xa2p(1,i)
               f1yp=f1yp*xa2p(2,i)
               f1zp=f1zp*xa2p(3,i)
               go to 590
c
  567          f1xm=f1xm*xa1m(1,i)
               f1ym=f1ym*xa1m(2,i)
               f1zm=f1zm*xa1m(3,i)
               f1xp=f1xp*xa1p(1,i)
               f1yp=f1yp*xa1p(2,i)
               f1zp=f1zp*xa1p(3,i)
c
  590          continue
c
               fcxm=fcxm+f1xm*w(i)
               fcym=fcym+f1ym*w(i)
               fczm=fczm+f1zm*w(i)
               fcxp=fcxp+f1xp*w(i)
               fcyp=fcyp+f1yp*w(i)
               fczp=fczp+f1zp*w(i)
 5100       continue
            xfcm(ind+indz)=fcxm
            yfcm(ind+indz)=fcym
            zfcm(ind+indz)=fczm*rw
            xfcp(ind+indz)=fcxp
            yfcp(ind+indz)=fcyp
            zfcp(ind+indz)=fczp*rw
 5110    continue
 5120 indz=indz+49
c
      end
c====================================================================
      subroutine frmcf(n1,n2,xm1,xp1)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common /xyzfc/ xfcm(343),yfcm(343),zfcm(343),
     *               xfcp(343),yfcp(343),zfcp(343)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      dimension nbeg(8), nend(8), nxtyp(84), nytyp(84), nztyp(84)
      dimension itm4(7)
      dimension xm1(*),xp1(*)
c
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
c
      data itm4/0,7,14,21,28,35,42/
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
         indx1=itm4(lax)
         indy1=itm4(lay)
         indz1=itm4(laz)
      do 40 ityp2=ib2,ie2
         indx=indx1+nxtyp(ityp2)
         indy=indy1+nytyp(ityp2)
         indz=indz1+nztyp(ityp2)
         summ=zero
         sump=zero
         go to (30,20,10,5,4,3,2), nroots
    2    summ=summ+xfcm(indx+294)*yfcm(indy+294)*zfcm(indz+294)
         sump=sump+xfcp(indx+294)*yfcp(indy+294)*zfcp(indz+294)
    3    summ=summ+xfcm(indx+245)*yfcm(indy+245)*zfcm(indz+245)
         sump=sump+xfcp(indx+245)*yfcp(indy+245)*zfcp(indz+245)
    4    summ=summ+xfcm(indx+196)*yfcm(indy+196)*zfcm(indz+196)
         sump=sump+xfcp(indx+196)*yfcp(indy+196)*zfcp(indz+196)
    5    summ=summ+xfcm(indx+147)*yfcm(indy+147)*zfcm(indz+147)
         sump=sump+xfcp(indx+147)*yfcp(indy+147)*zfcp(indz+147)
   10    summ=summ+xfcm(indx+98)*yfcm(indy+98)*zfcm(indz+98)
         sump=sump+xfcp(indx+98)*yfcp(indy+98)*zfcp(indz+98)
   20    summ=summ+xfcm(indx+49)*yfcm(indy+49)*zfcm(indz+49)
         sump=sump+xfcp(indx+49)*yfcp(indy+49)*zfcp(indz+49)
   30    summ=summ+xfcm(indx)*yfcm(indy)*zfcm(indz)
         sump=sump+xfcp(indx)*yfcp(indy)*zfcp(indz)
         ncount=ncount+1
         xm1(ncount)=summ
         xp1(ncount)=sump
   40 continue
c
      end
c======================================================================
      subroutine dipolex(bl,inx,ncf,ncs,na)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension inx(12,*)
      dimension dip(3)
c
       ntri=ncf*(ncf+1)/2
c
      call getival('inuc',inuc)
      call getival('ibas',ibas)
      call getival('ldensi',lden)
c
c  define (temporary) dipole matrix
      call matdef('dipp','s',ncf,ncf)
      call matdef('dens','s',ncf,ncf)
c
      call dcopy(ntri,bl(lden),1,bl(mataddr('dens')),1)
c
      call mmark
c     call immark
c
c  3=dipole
      call zeroit(bl(mataddr('dipp')),ncf*(ncf+1)/2)
      call inton(3,na,bl(mataddr('dipp')),inx,1,0,bl(ibas),
     1           bl(inuc),ncs)
cxxx
      write(6,*) ' the dip X (i|X|j) matrix'
      call matprint('dipp', 6  )
cxxx
      call matprodtr('dens','dipp',dip(1))

      call zeroit(bl(mataddr('dipp')),ncf*(ncf+1)/2)
      call inton(3,na,bl(mataddr('dipp')),inx,2,0,bl(ibas),
     1           bl(inuc),ncs)
cxxx
      write(6,*) ' the dip Y (i|Y|j) matrix'
      call matprint('dipp', 6  )
cxxx
      call matprodtr('dens','dipp',dip(2))
c
      call zeroit(bl(mataddr('dipp')),ncf*(ncf+1)/2)
      call inton(3,na,bl(mataddr('dipp')),inx,3,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(3))
c
c  print (i|z|j) matrix
c
      write(6,*) ' the dip Z (i|Z|j) matrix'
      call matprint('dipp', 6  )
c
      call retmark
c     call retimark
c
      call matrem('dens')
      call matrem('dipp')
c
c  the electronic part is negative because of the negative charge
      dip(1) = -dip(1)
      dip(2) = -dip(2)
      dip(3) = -dip(3)
c
      write(6,66)' electronic dipole moment =',dip(1),dip(2),dip(3)
c  calculate nuclear dipole and add to dip
      call nucdipole(na,dip,bl(inuc))
c
      write(6,66)' total      dipole moment =',dip(1),dip(2),dip(3)
  66  format(a16,2x,3(f15.9,1x))
      call f_lush(6)
c
c
      end
c======================================================================
      subroutine intof(natoms,inx,ifield,ifgrad,xfld,
     *                 datbas,datnuc,ncs,ncf,ntri,dn,dw,atfor)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c This routine calculates 1st derivatives of the following integrals :
c
c             overlap , kinetic , dipole & quadrupole
c
c last two for an external electric field & field gradient
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c ifield & ifgrad - indicate presence of an ext.electric field & f.grad.
c xfld(9)         - components of an external field (3) &field grad.(6)
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix
c dw(ntri)        - weighted density matrix
c INTENT(OUT)
c atfor(3,natoms) - forces on atoms - INPUT/OUTPUT
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension xfld(9)
      dimension dn(ntri),dw(ntri)   ! density & weight density
      dimension atfor(3,natoms)
clocal :
      dimension s3(3*784),t3(3*784)
      dimension f3a(3*784,3),f3b(3*784,3),fff(3*784,3) ! for an ext.el.field
      dimension f31(3*784),f32(3*784)                  ! for an ext.el.field
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
            do 50 j=1,ncs
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
c over all        if(iat.ne.jat) then
                     do l=1,len3
                        s3(l)=zero
                        t3(l)=zero
                     enddo
               
c
c................... overlap integral derivatives
c
                     key=1
                     call onef(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
c
c................... kinetic energy integral derivatives
c
                     key=2
                     call onef(i,j,igc,jgc,datbas,t3,inx,key,kk,k2)
c
c................... an external electric field & field gradient.....
                     key=1
                     if(ifield.eq.1) then
c ......................an external electric field
c                       (i10 /x,y or z/ j00) + (i00 /x,y or z/ j10)
c
                        do k1=1,3
                           do l=1,len3
                              f3a(l,k1)=zero
                              f3b(l,k1)=zero
                              fff(l,k1)=zero
                           enddo
                        enddo
c
                        do k1=1,3
                           field=xfld(k1)
                           if(field.ne.zero) then
                              do l=1,len3
                                 f3a(l,k1)=zero
                                 f3b(l,k1)=zero
                                 fff(l,k1)=zero
                                 f31(l)   =zero
                                 f32(l)   =zero
                              enddo
c
                              call onef(i,j,igc,jgc,datbas,f31,
     *                                  inx,key,k1,k2)
                              call onef(j,i,jgc,igc,datbas,f32,
     *                                  inx,key,k1,k2)
c
                              do l=1,len3
                                 f3a(l,k1)=f31(l)*field
                                 fff(l,k1)=f32(l)*field
                              enddo
                           endif    ! field.ne.zero
                        enddo       ! k1=1,3
                     endif          ! if(ifield.eq.1) then
c
c..................................................................
c                    if(ifgrad.eq.1) then
c                       an external electric field gradient
c                       (i10 /xy/ j00) + (i00 /xy/ j10)
c                       where xy=xx,yx,yy,zx,zy,zz
c
c                       do k1=4,9
c                          field=xfld(k1)
c                          if(field.ne.zero) then
c                             do l=1,len3
c                                f3(l,k1)=zero
c                             enddo
c
c                             call onef(i,j,igc,jgc,datbas,f3,
c    *                                  inx,key,k1,k2)
c
c                             do l=1,len3
c wrong                          t3(l)=t3(l)+f3(l,k1)*field
c                             enddo
c                          endif    ! field.ne.zero
c                       enddo       ! k1=4,9
c                    endif          ! if(igrad.eq.1) then
c..................................................................
c over all        endif             ! if(iat.ne.jat) then
c..................................................................
c reorder f3 like f3a is and put it into f3b :

                 if(ifield.eq.1) then
                    do k1=1,3
                    do icr=1,3
                       icrl=(icr-1)*len
                       do j3=1,len2
                          do i3=1,len1
                             ji=(j3-1)*len1+i3
                             ij=(i3-1)*len2+j3
                             f3b(icrl+ij,k1)=fff(icrl+ji,k1)
                          enddo
                       enddo
                    enddo
                    enddo
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
                     jj=jff*(jff-1)/2
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
c over all     IF(IAT.NE.JAT) THEN
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
c.................................................................
c the correct formula for this term is : gradient contribution is
c -Tr[S'W] where W=2CeC+, and S' is the overlap mtx. derivative
c s3 is d/dXi(S), so the contribution to the FORCES is  positive
c for the first center. The contribution of the second center is
c obtained from translational invariance
c.................................................................
                  atfor(icr,iat)=atfor(icr,iat)+s3(iij_icr)*dwij
                  atfor(icr,jat)=atfor(icr,jat)-s3(iij_icr)*dwij
c
                  atfor(icr,iat)=atfor(icr,iat)-t3(iij_icr)*dnij
                  atfor(icr,jat)=atfor(icr,jat)+t3(iij_icr)*dnij
c
                if(ifield.eq.1) then
                  do k1=1,3
                     atfor(icr,iat)=atfor(icr,iat)-f3a(iij_icr,k1)*dnij
                     atfor(icr,jat)=atfor(icr,jat)-f3b(iij_icr,k1)*dnij
                  enddo
                endif
c
               enddo
c over all     ENDIF
c
   40             continue
                  jfu=jfu+len2
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine vcdint(natoms, inx,  datbas, datnuc,  ncs,
     1                  ncf,    ntri, dn,     iatom, atrot)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c  This routine calculates overlap integrals over derivative basis 
c  functions of the type
c      <d(mu)/dBa||d(nu)/dXib>
c  and contracts them with the density to generate the atrot array:
c  Here (mu),(nu) are atomic basis functions, Ba is the a-th compononent
c  of the magnetic induction, and Xib is the b-th Cartesian direction of
c  the i-th nucleus. NOTE THAT THIS IS NOT A SECOND DERIVATIVE. The 
c  basis fuctions are assumed to be GIAO (London) atomic orbitals and 
c  thus depend on B directly.
c  The output array is dimensioned as atrot(3,3). The first index
c  is the magnetic direction, the second is the nuclear coordinate.
c
c ARGUMENTS
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix  (D(mu,nu)
c iatom           - atom for nucl. pert.      
c INTENT(OUT)
c atrot(3,3) - mixed overlap derivatives for Vibrational 
c                     Circular Dichroism (VCD)
c       atrot(mag,nuc)=Sum(mu,nu) D(mu,nu)*<d mu/dX(nuc)|d nu/d B(mag)>
c  where mu,nu=basis functions, X(nuc) is the nuc component of the 
c  nuclear coordinate associated with the  basis function mu, and B is 
c  the magnetic induction. Note that this quantity is not invariant
c  translationally
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension dn(ntri)
      dimension atrot(3,3)
clocal :
      dimension s3(784*9)
c--------------------------------
c
      kk=0
      k2=0
c
      ifu=0
      call zeroit(atrot,9)
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,ncs
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               do 45 jgc=0,ngcjx
                 ij=0
                 do mag=1,3
                   do nuc=1,3
                     do l=1,len
                        ij=ij+1
                        s3(ij)=zero
                     enddo
                   end do
                end do
c
c......... mixed magnetic and coordinate overlap integral derivatives
c
                     key=9
                     if(iat.eq.iatom) then
                     call onef(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
                     endif
c
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
                     jj=jff*(jff-1)/2
                     iij=(i1-1)*len2+j1
                     ij=ii+jff
                     if (jff.gt.iff) then
                       ij=jj+iff
                     end if  
                     dnij=dn(ij)
c
               do mag=1,3
                 do nuc=1,3
c.................................................................
          call vcdsum(mag,nuc,len,iij,s3,dnij,atrot)
                enddo
              end do
c
   40     continue
          jfu=jfu+len2
   45   continue
   50 continue
      ifu=ifu+len1
   55 continue
   60 continue
c---------------------------------------------------------------------
c
      end
c======================================================================
      subroutine vcdsum(mag,nuc,len,iij,s3,den,atrot)
      implicit real*8 (a-h,o-z)
c  This routine is used only to avoid the cumbersome index manipulations 
c  (addressing a 3-index quantity by linear indices)
c  It sums up the rotational intensities for a mixed 2nd derivative integral
c  and the density
c  ARGUMENTS
c  INTENT(IN)
c  mag=magnetic direction (1,2 or 3)
c  nuc=nuclear direction (1,2 or 3)
c  len=the number of integrals for a shell pair, the product of the
c     lengths of the shells
c  s3 = derivative integral 
c  den=density matrix element
c   sign=+/- 1
c  atrot(3,3)=rotational matrix
c
      dimension s3(len,3,3),atrot(3,3)
c
      atrot(mag,nuc)=atrot(mag,nuc)-s3(iij,nuc,mag)*den
      end
c======================================================================
      subroutine vcdint12(natoms, inx,  datbas, datnuc,  ncs,
     1                    ncf,    inttp,  iatom,   sder)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c---------------------------------------------------------------------
c  This routine calculates overlap integrals over derivative basis 
c  functions of the type
c
c  <d(mu)/dXib|(nu)> or <(mu)|d(nu)/dBa> for itype=1 and itype=2, resp.
c
c  Here (mu),(nu) are atomic basis functions, Ba is the a-th compononent
c  of the magnetic induction, and Xib is the b-th Cartesian direction of
c  the i-th nucleus. NOTE THAT THIS IS NOT A DERIVATIVE OF THE INTEGRAL.
c  The basis fuctions are assumed to be GIAO (London) atomic orbitals 
c  and thus depend on B directly. Note that only the integrals are
c  calculated; contraction with density-like matrices is done in the
c  calling program. Note further that the first function is 
c  differentiated for the nuclei and the second for the magnetic field.
c  The output array is dimensioned as sder(ncf,ncf,3), a full matrix.
c
c ARGUMENTS
c---------------------------------------------------------------------
c INTENT(IN):
c natoms          - number of atoms
c inx(12,ncs)     - basis set data
c datbas(13,nsh)  - basis set data
c datnuc(5,natoms)- nuclear data
c ncs,ncf         - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c inttp           - 1 for nuclear perturbation, 2 for magnetic deriv.
c iatom           - the atom in question if itype=2, otherwise dummy
c INTENT(OUT)
c sder(ncf,ncf,3)   - mixed overlap derivatives for Vibrational 
c                     Circular Dichroism (VCD)
c       atrot(mag,nuc)=Sum(mu,nu) D(mu,nu)*<d mu/dX(nuc)|d nu/d B(mag)>
c  where mu,nu=basis functions, X(nuc) is the nuc component of the 
c  nuclear coordinate associated with the  basis function mu, and B is 
c  the magnetic induction. Note that this quantity is not invariant
c  translationally
c---------------------------------------------------------------------
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension sder(ncf,ncf,3)   
clocal :
      dimension s3(784*9)
c--------------------------------
c     iprint=igetival('nprint')
c
      call zeroit(sder,3*ncf**2)
      ifu=0
      do 60 i=1,ncs
        len1=inx(3,i)
        iat=inx(2,i)
        ngci=inx(4,i)
        do 55 igc=0,ngci
          jfu=0
          do 50 j=1,ncs
            len2=inx(3,j)
            jat=inx(2,j)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
            len=len1*len2
            do 45 jgc=0,ngcjx
              ij=0
              do mag=1,3 !direction of the nuclear coordinate
c    (or the magnetic field)
                do l=1,len
                  ij=ij+1
                  s3(ij)=zero
                enddo
              end do
c
              if(inttp.eq.1) then
                key=10    ! nuclear perturbation of (mu)
              else if(inttp.eq.2) then
                key=11    ! magnetic perturbation of (nu)
              else
                call nerror(1,'vcdint12','integral type must be 1 or 2',
     1               inttp,inttp)
              end if
CPP
c      print *,'Calling shells', i,j
c  in the nuclear perturbation case, derivatives are calculated for 
c  1 atom only
              if(inttp.eq.1.and.iat.eq.iatom.or.inttp.eq.2) then
                call onef(i,j,igc,jgc,datbas,s3,inx,key,kk,k2)
              end if
c
c..................................................................
c
              iij=0
              iff=ifu
              do 40 i1=1,len1
                iff=iff+1
                jff=jfu
                do 40 j1=1,len2
                  jff=jff+1
                  iij=iij+1
                  mag1=0
                  do mag=1,3  ! direction
                    sder(iff,jff,mag)=sder(iff,jff,mag)+s3(iij+mag1)
                    mag1=mag1+len
                  end do
c.................................................................
c
   40         continue
              jfu=jfu+len2
   45       continue
   50     continue
          ifu=ifu+len1
   55   continue
   60 continue
c---------------------------------------------------------------------
      
      end
c======================================================================
