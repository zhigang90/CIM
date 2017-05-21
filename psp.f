c...
c...  MM March-December 2004
c...
c...  thys file contains the driver routines for the calculation
c...  of integrals over ab initio pseudopotentials.
c...  The method of calculation is based on
c...
c...  L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44, 289 (1981)
c...
c...  Altough all the code has been originally written,
c..   this work has been grately facilitated by the the availability
c...  of Fortran routines from the ARGOS code that have been used as
c...  model and for testing.
c...
c...
c...       About the speed of the integral calculation
c...
c...   currently 2 types of threshold testing are implemented in order
c...   to speed up the calculation:
c...
c...   1) Test over the argument of the exponential prefactor
c...      (exparg) that is done in routines psp_type12xxx.
c...      The threshold for this term is currently set to 100.
c...      A lower threshold will furter speed up the calculation,
c...      but might start to affect the accuracy of the integrals
c...
c...   2) Test over the radial integrals. For this test, the one
c...      electron integral threshold is used (in many cases, this might
c...      have been further strengthened by a factor 0.01). This test
c...      is implemented in global form in routine psp_ar12, and in
c...      local form in the routines that compute the individual
c...      radial integrals. Note that the threshold used in the test
c...      is the logarithm of the actual integral threshold.
c...
c...   One possibility to further speed up the calculation might be
c...   to introduce a further test over the value of the radial
c...   integrals in subroutine sumar2. Currently, a test is done that
c...   will only discard integrals exactly equal to zero, and at the
c...   current stage this is mostly redundant because of the global
c...   test of point 2 above. Figuring out a threshold for this test
c...   might do some good.
c...
c======================================================================
      subroutine psph0(na,npsp,ncf,ncs,inx,xnuc,basdat,tol,hpsp)

      use memory

      implicit real*8 (a-h,o-z)
c...
c...  driver routine for the psp contribution to the one electron
c...  hamiltonian
c...
c...  na      number of atoms
c...  npsp    number of pseudopotentials
c...  ncf     number of contracted functions
c...  ncs     number of contracted shells
c...  inx     contraction information
c...  xnuc    nuclear data
c...  basdat  basis set data
c...  tol     tolerance threshold
c...  hpsp    in output will contain psp one electron contribution
c...
      dimension inx(12,*),basdat(13,*),xnuc(5,*),hpsp(*)
c     common /big/bl(1)
      dimension time(20)
c...
c     call secund(t0)
c     call zeroit(time,20)
c...
c...  get the psp parameters and pointers
c...
      call getival('nradp',nradp)       ! number of radial terms
      call getival('maxlpsp',maxlpsp)   ! maximum l value we can handle
      call getival('iiecpdat',iiecpdat) ! ptr. to psp info
      call getival('inrrad',inrrad)     ! ptr. to r exponents
      call getival('iarad',iarad)       ! ptr. to radial exponents
      call getival('icrad',icrad)       ! ptr. to radial coefficients
c...
c...  parameters for memory setup:
c...
c...  maxlp        maximum angular momentum of psp
c...  maxl         maximum angular momentum of Gaussians
c...  maxd2        maximum number of individual components for
c...               of a contracted shell (number of gaussians times
c...               contaction length or no. of general contractions)
c...
      call pspgetpar(npsp,maxlpsp,ncs,inx,bl(iiecpdat),maxlp,maxl,
     $               maxd2)
c...
c...  memory setup
c...
      mlp=max(maxl,maxlp)
      maxpw=max(maxl,2)+1
      call mmark
      call getmem(maxpw,icax)             ! powers of ca
      call getmem(maxpw,icay)             ! powers of ca
      call getmem(maxpw,icaz)             ! powers of ca
      call getmem(maxpw,icbx)             ! powers of cb
      call getmem(maxpw,icby)             ! powers of cb
      call getmem(maxpw,icbz)             ! powers of cb
      call getmem((maxl+mlp+1)**2,iyra)   ! real sperical harmonics
      call getmem((maxl+maxlp+1)**2,iyrb) ! real sperical harmonics
      call getmem((2*maxl+1)*(maxl+maxlp+1)**2,iq) ! radial integrals
      call getmem((4*mlp+1)**3,ixyzi)     ! precomputed angular terms
      call getmem((maxl+1)**6,iar) ! sums of radial and angular int.
      call getmem(maxlp+1,ilen)           ! array for psp contraction
      call getmem(maxd2,icna)             ! normalization constants
      call getmem(maxd2,icnb)             ! normalization constants
c...
      call getmem((maxd2)**2,ipspi)       ! integrals over primitives
      call getmem((maxd2)**2,ipspic)      ! integrals over contracted
c...
c...  zero out one electron matrix
c...
      ntri=ncf*(ncf+1)/2
      call zeroit(hpsp,ntri)
c...
c...  generate a table of double factorial products
c...  for angular integrals
c...
      mijk=4*mlp
      call preomega(mijk,bl(ixyzi))
c...
c...  loop over contracted shells
c...
      do ics=1,ncs
        do jcs=1,ics
c...
c...  compute integrals
c...
          call psp_int0(ics,       jcs,     npsp,    nradp,  maxlpsp,
     $                  ncs,       ncf,    maxlp,     maxl,     mijk,
     $         bl(iiecpdat),bl(inrrad),bl(iarad),bl(icrad), bl(ilen),
     $                 xnuc,       inx,   basdat, bl(icax), bl(icay),
     $             bl(icaz),  bl(icbx), bl(icby), bl(icbz), bl(iyra),
     $             bl(iyrb),    bl(iq),bl(ixyzi),  bl(iar),bl(ipspi),
     $           bl(ipspic),  bl(icna), bl(icnb),      tol,     hpsp,
     $                 time)
        enddo
      enddo
c...
c...  done. Release memory
c...
      call retmark
c     call secund(t1)
c...
c     call getival('iout',iout)
c     write(iout,*)
c     write(iout,*)'Total time for psp contribution       ',t1-t0
c     write(iout,*)
c     write(iout,*)'Overall integrals over primitives     ',time(8)
c     write(iout,*)'Zeroing integrals arrays              ',time(9)
c     write(iout,*)'Time for radial and angular integrals ',time(1)
c     write(iout,*)'Time for integrals over primitives    ',time(7)
c     write(iout,*)'       Time for summation loop        ',time(2)
c     write(iout,*)'       Time for assembling integrals  ',time(6)
c     write(iout,*)'Time for normalization/contraction    ',time(3)
c     write(iout,*)'Time for cartesian to spherical       ',time(4)
c     write(iout,*)'Time for fillup                       ',time(5)
c     write(iout,*)
c     write(iout,*)'Time for type 1 radial integrals      ',time(11)
c     write(iout,*)'Time for spherical harmonics          ',time(12)
c     write(iout,*)'Time for sumar1                       ',time(13)
c     write(iout,*)'Time for sumar0                       ',time(14)
c     write(iout,*)'Time for type 2 radial integrals      ',time(15)
c     write(iout,*)'Time for sumar2                       ',time(16)
c...
      end
c======================================================================
      subroutine pspgetpar(npsp,maxlpsp,ncs,inx,iecpdat,maxlp,maxl,
     $                     maxd2)
      implicit real*8 (a-h,o-z)
c...
c...  extract dimensional parameters for psp integral calculation
c...
c...  npsp    number of pseudopotentials
c...  maxlpsp maximum l value of psp we can handle
c...  ncs     number of contracted shells
c...  inx     contraction info
c...  iecpdat psp info
c...  maxlp   (out) maximum l value of psp
c...  maxl    (out) maximum l value of Gaussians
c...  maxd2   (out) maximum number of individual shell components
c...
      dimension inx(12,ncs),iecpdat(maxlpsp+6,npsp)
      dimension ltyp(13),ng(13)
      data ltyp/0,1,1,2,2,3,3,4,5,6,4,5,6/
      data ng/1,3,4,6,6,10,10,15,21,28,15,21,28/
c...
      maxlp=0
      do i=1,npsp
        maxlp=max(maxlp,iecpdat(1,i))
      enddo
      maxl=0
      maxd2=0
      do i=1,ncs
        lenc=inx(5,i)-inx(1,i)
        nng=ng(inx(12,i))
        ngc=inx(4,i)
        nd2=nng*max(lenc,ngc)
        maxd2=max(maxd2,nd2)
        maxl=max(maxl,ltyp(inx(12,i)))
      enddo
      end
c======================================================================
      subroutine psp_int0(ics,  jcs,   npsp,  nradp,maxlpsp,
     $                    ncs,  ncf,  maxlp,   maxl,   mijk,
     $                iecpdat,nrrad,  aradp,  cradp,    len,
     $                   xnuc,  inx, basdat,    cax,    cay,
     $                    caz,  cbx,    cby,    cbz,    yra,
     $                    yrb,    q,   xyzi,     ar,   pspi,
     $                  pspic,  cna,    cnb,    tol,   hpsp,
     $                   time)
      implicit real*8 (a-h,o-z)
c...
c...  compute psp contribution to one electron Hamiltonian for
c...  the contracted shell pair ics,jcs
c...
c...  ics        contracted shell i
c...  jcs        contracted shell j
c...  npsp       number of pseudopotentials
c...  nradp      number of radial terms
c...  maxlpsp    maximum angular momentum of psp we can handle
c...  ncs        number of contracted shells
c...  ncf        number of contracted basis functrions
c...  maxlp      maximum angular momentum of psp
c...  maxl       maximum angular momentum of Gaussians
c...  mijk       dimension for xyzi array
c...  iecpdat    array of psp info
c...  nrrad      array of r exponent of radial terms
c...  aradp      array of radial term exponents
c...  cradp      array of radial term coefficients
c...  len        scratch area for psp contraction length
c...  xnuc       array of nuclear info
c...  inx        array of contraction info
c...  basdat     array of  basis set data
c...  cax        scratch area for powers of cax
c...  cay        scratch area for powers of cay
c...  caz        scratch area for powers of caz
c...  cbx        scratch area for powers of cbx
c...  cby        scratch area for powers of cby
c...  cbz        scratch area for powers of cbz
c...  yra        scratch area for real spherical harmonics
c...  yrb        scratch area for real spherical harmonics
c...  q          scratch area for radial integrals
c...  xyzi       area for precomputed angular terms
c...  ar         scratch area for sums of radial and angular ints.
c...  pspi       scratch area for integrals over Gaussians primitives
c...  pspic      scratch area for integrals over Gaussians contracted
c...  cna        scratch are for normalization coefficients
c...  cnb        scratch are for normalization coefficients
c...  tol        tolerance threshold
c...  hpsp       psp contribution to one electron Hamiltonian
c...
      dimension iecpdat(maxlpsp+6,npsp),nrrad(nradp),aradp(nradp),
     $          cradp(nradp),len(0:maxlp)
      dimension xnuc(5,*),inx(12,*),basdat(13,*)
      dimension cax(0:maxl),cay(0:maxl),caz(0:maxl)
      dimension cbx(0:maxl),cby(0:maxl),cbz(0:maxl)
      dimension pspi(*),pspic(*),cna(*),cnb(*)
      dimension hpsp(*)
      dimension time(*)
c...
      dimension A(3),B(3),C(3),CA(3),CB(3),lima(3),limb(3)
c...
c...  the type of the Gaussian shell is given by ityp:
c...
c...  ityp  1  2  3   4   5   6    7    8   9    10   11  12   13
c...        s  p  l   d5  d6  f7  f10  g15  h21  i28  g9  h11  i13
c...
c...   the integrals are computed over the cartesian components
c...   (d6, f15, g15, h21, i28). The spherical harmonic components,
c...   if needed, will be obtained by linear combination.
c...
c...   data for number of components and angular momenta
c...
      dimension ityp1(13),ng(8),ns(8),lmax(8),mulx(88),muly(88),
     $          mulz(88)
      data ityp1/1,2,3,4,4,5,5,6,7,8,6,7,8/
      data ng/1,3,4,6,10,15,21,28/
      data ns/1,2,5,9,15,25,40,61/
      data lmax/0,1,1,2,3,4,5,6/
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
c...
      logical transa,transb,contra,contrb
      parameter (zero=0.0d0)
      parameter (argmax=100.0d0)
      tolp=tol
c...
c...  shell i:
c...
      ngci=inx(4,ics)+1          ! no. of general contractions
      ifi=inx(1,ics)             ! first primitive
      ifl=inx(5,ics)             ! last primitive
      ilen=ifl-ifi               ! contraction length
      iatom=inx(2,ics)           ! atomic center
      a(1)=xnuc(2,iatom)
      a(2)=xnuc(3,iatom)
      a(3)=xnuc(4,iatom)
      itypa=inx(12,ics)          ! type of shell i
      transa=itypa.eq.4.or.itypa.eq.6.or.itypa.ge.11
      ityp1a=ityp1(itypa)
      lmaxa=lmax(ityp1a)         ! angular momentum
      nga=ng(ityp1a)             ! no. of cartesian components
      nsa=ns(ityp1a)             ! ptr. for angular momenta vector
      nsha=inx(3,ics)            ! shell size
      ifci=inx(11,ics)           ! first contracted (minus one)
c...
c...  check if we can contract over shell i at the
c...  radial and angular integral level, and
c...  precompute normalization constants.
c...
      contra=itypa.ne.3.and.ngci.eq.1.and.ilen.gt.1
      if(contra)then
      do ip=1,ilen
         ipi=ifi+ip
         cna(ip)=basdat(2,ipi)*(2.0d0*basdat(1,ipi))**0.75d0
      enddo
      endif
c...
c...  shell j:
c...
      ngcj=inx(4,jcs)+1          ! no. of general contractions
      jfi=inx(1,jcs)             ! first primitive
      jfl=inx(5,jcs)             ! last primitive
      jlen=jfl-jfi               ! contraction length
      jatom=inx(2,jcs)           ! atomic center
      b(1)=xnuc(2,jatom)
      b(2)=xnuc(3,jatom)
      b(3)=xnuc(4,jatom)
      itypb=inx(12,jcs)          ! type of shell j
      transb=itypb.eq.4.or.itypb.eq.6.or.itypb.ge.11
      ityp1b=ityp1(itypb)
      lmaxb=lmax(ityp1b)         ! angular momentum
      ngb=ng(ityp1b)             ! no. of cartesian components
      nsb=ns(ityp1b)             ! ptr. for angular momenta vector
      nshb=inx(3,jcs)            ! shell size
      ifcj=inx(11,jcs)           ! first contracted (minus one)
      mxpw=max(2,lmaxa,lmaxb)
c...
c...  check if we can contract over shell j at the
c...  radial and angular integral level, and
c...  precompute normalization constants.
c...
      contrb=itypb.ne.3.and.ngcj.eq.1.and.jlen.gt.1
      if(contrb)then
      do jp=1,jlen
         jpj=jfi+jp
         cnb(jp)=basdat(2,jpj)*(2.0d0*basdat(1,jpj))**0.75d0
      enddo
      endif
c     if(debug)then
c         call getival('iout',iout)
c         write(iout,*)'iatom,itypa,lmaxa',iatom,itypa,lmaxa
c         write(iout,*)'jatom,itypb,lmaxb',jatom,itypb,lmaxb
c     endif
c...
c...  zero out integrals
c...
c     call secund(t0)
      call zeroit(pspic,nga*ngb*ngci*ngcj)
      call zeroit(pspi,nga*ngb*ilen*jlen)
c     call secund(tz)
c     time(9)=time(9)+tz-t0
c...
c...  loop over pseudopotentials
c...
      do ipsp=1,npsp
        icatom=iecpdat(3,ipsp)   ! atomic center of psp
        lpsp=iecpdat(1,ipsp)     ! maximum l value of psp
        irl=iecpdat(4,ipsp)
        irs=iecpdat(5,ipsp)
c...
c...  compute number of terms of local and non local parts
c...
        len0=irs-irl             ! number of terms in local part
        do l=0,lpsp-1,1
          len(l)=iecpdat(6+l,ipsp)-iecpdat(5+l,ipsp)
        enddo
        len(lpsp)=iecpdat(12,ipsp)-iecpdat(5+lpsp,ipsp)+1
c...
c...   distances CA and CB
c...
        c(1)=xnuc(2,icatom)
        c(2)=xnuc(3,icatom)
        c(3)=xnuc(4,icatom)
        CA(1)=C(1)-A(1)
        CA(2)=C(2)-A(2)
        CA(3)=C(3)-A(3)
        CB(1)=C(1)-B(1)
        CB(2)=C(2)-B(2)
        CB(3)=C(3)-B(3)
c       if(debug)then
c         call getival('iout',iout)
c         write(iout,*)'lpsp,len0,len',lpsp,len0,(len(ijz),ijz=0,lpsp)
c         write(iout,*)'nr',(nrrad(irl+ijz),ijz=0,len0-1)
c         write(iout,*)'coeff',(cradp(irl+ijz),ijz=0,len0-1)
c         write(iout,*)'exp',(aradp(irl+ijz),ijz=0,len0-1)
c         isk=0
c         do l=0,lpsp,1
c           write(iout,*)
c           write(iout,*)'l,len',l,len(l)
c           write(iout,*)'nr',(nrrad(irs+isk+ijz),ijz=0,len(l)-1)
c           write(iout,*)'coeff',(cradp(irs+isk+ijz),ijz=0,len(l)-1)
c           write(iout,*)'exp',(aradp(irs+isk+ijz),ijz=0,len(l)-1)
c           isk=isk+len(l)
c         enddo
c         write(iout,*)
c         write(iout,*)'a',a
c         write(iout,*)'b',b
c         write(iout,*)'c',c
c         write(iout,*)'ca',ca
c         write(iout,*)'cb',cb
c         write(iout,*)'mxpw',mxpw
c       endif
c...
c...  compute powers of distances ca and cb and do loops limits
c...
        call zeroit(cax,mxpw+1)
        call zeroit(cay,mxpw+1)
        call zeroit(caz,mxpw+1)
        call zeroit(cbx,mxpw+1)
        call zeroit(cby,mxpw+1)
        call zeroit(cbz,mxpw+1)
        call powcacb(CA(1),CA(2),CA(3),CB(1),CB(2),CB(3),
     $               mxpw,lima,limb,cax,cay,caz,cbx,cby,cbz)
        ca2=cax(2)+cay(2)+caz(2)
        cb2=cbx(2)+cby(2)+cbz(2)
c...
c...  compute the integrals using the appropriate routine
c...  according to the level of contraction that can
c...  be performed at the radial and angular integral level
c...
        if(contra.and.contrb)then ! contraction over shells a and b
        call psp_type12cab(
     $              basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $           muly(nsb), mulz(nsb),       nga,       ngb,      lpsp,
     $          nrrad(irl),aradp(irl),cradp(irl),      len0,       len,
     $                 cax,       cay,       caz,       cbx,       cby,
     $                 cbz,      lima,      limb,      mxpw,     lmaxa,
     $               lmaxb,      mijk,       yra,       yrb,      xyzi,
     $                   q,        ar,      tolp,    argmax,      ilen,
     $                 ifi,       ca2,      jlen,       jfi,       cb2,
     $                 cna,       cnb,     pspic,      time)
        else if(contrb) then   ! contraction over shell b
        call psp_type12cb(
     $              basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $           muly(nsb), mulz(nsb),       nga,       ngb,      lpsp,
     $          nrrad(irl),aradp(irl),cradp(irl),      len0,       len,
     $                 cax,       cay,       caz,       cbx,       cby,
     $                 cbz,      lima,      limb,      mxpw,     lmaxa,
     $               lmaxb,      mijk,       yra,       yrb,      xyzi,
     $                   q,        ar,      tolp,    argmax,      ilen,
     $                 ifi,       ca2,      jlen,       jfi,       cb2,
     $                 cna,       cnb,      pspi,      time)
        else if(contra) then  ! contraction over shell a
        call psp_type12ca(
     $              basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $           muly(nsb), mulz(nsb),       nga,       ngb,      lpsp,
     $          nrrad(irl),aradp(irl),cradp(irl),      len0,       len,
     $                 cax,       cay,       caz,       cbx,       cby,
     $                 cbz,      lima,      limb,      mxpw,     lmaxa,
     $               lmaxb,      mijk,       yra,       yrb,      xyzi,
     $                   q,        ar,      tolp,    argmax,      ilen,
     $                 ifi,       ca2,      jlen,       jfi,       cb2,
     $                 cna,       cnb,      pspi,      time)
        else   ! no contraction at this stage
        call psp_type12i(
     $              basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $           muly(nsb), mulz(nsb),       nga,       ngb,      lpsp,
     $          nrrad(irl),aradp(irl),cradp(irl),      len0,       len,
     $                 cax,       cay,       caz,       cbx,       cby,
     $                 cbz,      lima,      limb,      mxpw,     lmaxa,
     $               lmaxb,      mijk,       yra,       yrb,      xyzi,
     $                   q,        ar,      tolp,    argmax,      ilen,
     $                 ifi,       ca2,      jlen,       jfi,       cb2,
     $                pspi,      time)
        endif
      enddo   ! end loop over pseudopotentials
c     call secund(t1)
c     time(8)=time(8)+t1-t0
c...
c...  if we did not contract at the integral level, contract now
c...
c     call secund(t0)
      if(contra.and.contrb)then ! all is already contracted
      else if(contrb)then ! contract over a only
      call psp_contrna(ilen,ifi,ngci,nga,itypa,
     $                jlen,jfi,ngcj,ngb,itypb,
     $                basdat,1,pspi,pspic)
      else if(contra)then ! contract over b only
      call psp_contrnb(ilen,ifi,ngci,nga,itypa,
     $                jlen,jfi,ngcj,ngb,itypb,
     $                basdat,1,pspi,pspic)
      else ! contract over both a and b
      call psp_contrn(ilen,ifi,ngci,nga,itypa,
     $                jlen,jfi,ngcj,ngb,itypb,
     $                basdat,1,pspi,pspic)
      endif
c     call secund(t1)
c     time(3)=time(3)+t1-t0
c...
c...  eventually transform to spherical components.
c...  pspi is used as scratch space
c...
c     call secund(t0)
      if(transa.or.transb)then
        call cart2sphera(pspic,ngci*ngcj,itypa,transa,nga,nsha,
     $                           itypb,transb,ngb,nshb,pspi)
      endif
c     call secund(t1)
c     time(4)=time(4)+t1-t0
c...
c...  fill one electron matrix
c...
c     call secund(t0)
      ici=ifci
      nd3=ngci*ngcj*nga
      nd2=ngci*ngcj
      do igc=1,ngci
        do ica=1,nsha
          ici=ici+1
          ii=ici*(ici-1)/2
          icj=ifcj
          do jgc=1,ngcj
            do icb=1,nshb
              icj=icj+1
              if(ici.ge.icj)then
                ind=(icb-1)*nd3+(ica-1)*nd2+(jgc-1)*ngci+igc
                hpsp(ii+icj)=pspic(ind)
              endif
            enddo
          enddo
        enddo
      enddo
c     call secund(t1)
c     time(5)=time(5)+t1-t0
c...
      end
c=======================================================================
      subroutine psp_contrn(ilen,ifi,ngci,nga,itypa,
     $                      jlen,jfi,ngcj,ngb,itypb,
     $                       basdat,ndim,pspi,pspic)
      implicit real*8 (a-h,o-z)
c...
c...  normalizes and contracts a batch of psp integrals over
c...  primitives for both sells I and j
c...
c...  ilen        contraction length of shell i
c...  ifi         initial pointer to primitives of shell i
c...  ngci        number of general contractions on shell i
c...  nga         number of Gaussians of shell i (cartesian)
c...  itypa       type of shell i
c...  jlen        contraction length of shell j
c...  jfi         initial pointer to primitives of shell j
c...  ngcj        number of general contractions on shell j
c...  ngb         number of Gaussians of shell i (cartesian)
c...  itypb       type of shell j
c...  basdat      basis set info
c...  ndim        mumber of components of integrals
c...  pspi        integrals over primitives
c...  pspic       in exit will contain the integrals over contracted
c...
c...  a few words on normalization constants:
c...
c...  each Gaussian need to be multiplied by a normalization factor
c...
c...         N1=(2*a/%pi)**(3/4)*(4*a)**(l/2)
c...
c...   where a is the exponent and l is the angular momentum.
c...   (this assuming the integrals are normalized to
c...     (2*n_x-1)!!*(2*n_y-1)!!*(2*n_z-1)!!)
c...
c...  the factor (1/%pi)**(3/4) is already taken into account
c...  during the calculation of the radial integrals, while
c...  the factor (4*a)**(l/2) has already been included
c...  in the normalization of the primitives.
c...  Thus here the integrals need to be multipled only by the factor
c...
c...               NF=(2*a)**(3/4)
c...
c...        note that there are two such factors, one for each
c...        of the two Gaussians involved in the integral.
c...
      dimension basdat(13,*)
      dimension pspi(ndim,ilen,jlen,nga,ngb)
      dimension pspic(ndim,ngci,ngcj,nga,ngb)
c...
      do igc=1,ngci
        do jgc=1,ngcj
          do ip=1,ilen
            ipi=ifi+ip
            alpa=basdat(1,ipi)
            if(itypa.eq.3)then
              csa=basdat(igc+1,ipi)*(2.0d0*alpa)**0.75d0
              cpa=basdat(3,ipi)*(2.0d0*alpa)**0.75d0
            else
              csa=basdat(igc+1,ipi)*(2.0d0*alpa)**0.75d0
            endif
            do jp=1,jlen
              jpj=jfi+jp
              alpb=basdat(1,jpj)
              if(itypb.eq.3)then
                csb=basdat(jgc+1,jpj)*(2.0d0*alpb)**0.75d0
                cpb=basdat(3,jpj)*(2.0d0*alpb)**0.75d0
              else
                csb=basdat(jgc+1,jpj)*(2.0d0*alpb)**0.75d0
              endif
              do i=1,nga
                if(itypa.eq.3.and.i.gt.1)then
                  coefi=cpa
                else
                  coefi=csa
                endif
                do j=1,ngb
                  if(itypb.eq.3.and.j.gt.1)then
                    coefj=cpb
                  else
                    coefj=csb
                  endif
                  coefij=coefi*coefj
                  do ic=1,ndim
                    pspic(ic,igc,jgc,i,j)=pspic(ic,igc,jgc,i,j)+
     $                                 pspi(ic,ip,jp,i,j)*coefij
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      end
c=======================================================================
      subroutine psp_contrna(ilen,ifi,ngci,nga,itypa,
     $                      jlen,jfi,ngcj,ngb,itypb,
     $                       basdat,ndim,pspi,pspic)
      implicit real*8 (a-h,o-z)
c...
c...  normalizes and contracts a batch of psp integrals over
c...  primitives of shell a
c...
c...  ilen        contraction length of shell i
c...  ifi         initial pointer to primitives of shell i
c...  ngci        number of general contractions on shell i
c...  nga         number of Gaussians of shell i (cartesian)
c...  itypa       type of shell i
c...  jlen        contraction length of shell j
c...  jfi         initial pointer to primitives of shell j
c...  ngcj        number of general contractions on shell j
c...  ngb         number of Gaussians of shell i (cartesian)
c...  itypb       type of shell j
c...  basdat      basis set info
c...  ndim        mumber of components of integrals
c...  pspi        integrals over primitives
c...  pspic       in exit will contain the integrals over contracted
c...
c...  a few words on normalization constants:
c...
c...  each Gaussian need to be multiplied by a normalization factor
c...
c...         N1=(2*a/%pi)**(3/4)*(4*a)**(l/2)
c...
c...   where a is the exponent and l is the angular momentum.
c...   (this assuming the integrals are normalized to
c...     (2*n_x-1)!!*(2*n_y-1)!!*(2*n_z-1)!!)
c...
c...  the factor (1/%pi)**(3/4) is already taken into account
c...  during the calculation of the radial integrals, while
c...  the factor (4*a)**(l/2) has already been included
c...  in the normalization of the primitives.
c...  Thus here the integrals need to be multipled only by the factor
c...
c...               NF=(2*a)**(3/4)
c...
      dimension basdat(13,*)
      dimension pspi(ndim,ilen,1,nga,ngb)
      dimension pspic(ndim,ngci,1,nga,ngb)
c...
      do ip=1,ilen
        ipi=ifi+ip
          do igc=1,ngci
            alpa=basdat(1,ipi)
            if(itypa.eq.3)then
              csa=basdat(igc+1,ipi)*(2.0d0*alpa)**0.75d0
              cpa=basdat(3,ipi)*(2.0d0*alpa)**0.75d0
            else
              csa=basdat(igc+1,ipi)*(2.0d0*alpa)**0.75d0
            endif
              do i=1,nga
                if(itypa.eq.3.and.i.gt.1)then
                  coefi=cpa
                else
                  coefi=csa
                endif
                do j=1,ngb
                  do ic=1,ndim
                    pspic(ic,igc,1,i,j)=pspic(ic,igc,1,i,j)+
     $                                 pspi(ic,ip,1,i,j)*coefi
                  enddo
                enddo
              enddo
          enddo
      enddo
      end
c=======================================================================
      subroutine psp_contrnb(ilen,ifi,ngci,nga,itypa,
     $                      jlen,jfi,ngcj,ngb,itypb,
     $                       basdat,ndim,pspi,pspic)
      implicit real*8 (a-h,o-z)
c...
c...  normalizes and contracts a batch of psp integrals over
c...  primitives of shell j
c...
c...  ilen        contraction length of shell i
c...  ifi         initial pointer to primitives of shell i
c...  ngci        number of general contractions on shell i
c...  nga         number of Gaussians of shell i (cartesian)
c...  itypa       type of shell i
c...  jlen        contraction length of shell j
c...  jfi         initial pointer to primitives of shell j
c...  ngcj        number of general contractions on shell j
c...  ngb         number of Gaussians of shell i (cartesian)
c...  itypb       type of shell j
c...  basdat      basis set info
c...  ndim        mumber of components of integrals
c...  pspi        integrals over primitives
c...  pspic       in exit will contain the integrals over contracted
c...
c...  a few words on normalization constants:
c...
c...  each Gaussian need to be multiplied by a normalization factor
c...
c...         N1=(2*a/%pi)**(3/4)*(4*a)**(l/2)
c...
c...   where a is the exponent and l is the angular momentum.
c...   (this assuming the integrals are normalized to
c...     (2*n_x-1)!!*(2*n_y-1)!!*(2*n_z-1)!!)
c...
c...  the factor (1/%pi)**(3/4) is already taken into account
c...  during the calculation of the radial integrals, while
c...  the factor (4*a)**(l/2) has already been included
c...  in the normalization of the primitives.
c...  Thus here the integrals need to be multipled only by the factor
c...
c...               NF=(2*a)**(3/4)
c...
      dimension basdat(13,*)
      dimension pspi(ndim,1,jlen,nga,ngb)
      dimension pspic(ndim,1,ngcj,nga,ngb)
c...
        do jp=1,jlen
          jpj=jfi+jp
            do jgc=1,ngcj
              alpb=basdat(1,jpj)
              if(itypb.eq.3)then
                csb=basdat(jgc+1,jpj)*(2.0d0*alpb)**0.75d0
                cpb=basdat(3,jpj)*(2.0d0*alpb)**0.75d0
              else
                csb=basdat(jgc+1,jpj)*(2.0d0*alpb)**0.75d0
              endif
              do j=1,ngb
                if(itypb.eq.3.and.j.gt.1)then
                  coefj=cpb
                else
                  coefj=csb
                endif
                do i=1,nga
                  do ic=1,ndim
                    pspic(ic,1,jgc,i,j)=pspic(ic,1,jgc,i,j)+
     $                                 pspi(ic,1,jp,i,j)*coefj
                  enddo
                enddo
              enddo
            enddo
        enddo
      end
c=======================================================================
      subroutine cart2sphera(xc,ndim,itypa,transa,nga,nsha,
     $                         itypb,transb,ngb,nshb,xs)
      implicit real*8 (a-h,o-z)
c...
c...  transforms a batch of integrals from cartesian components
c...  to spherical components:
c...
c...      Cartesian                   Spherical
c...
c...               d  functions
c...
c...   1     xx                1/sqrt(12)*(2*zz-xx-yy)
c...   2     yy                  1/2*(xx-yy)
c...   3     zz                        xy
c...   4     xy                        xz
c...   5     xz                        yz
c...   6     yz
c...
c...               f  functions
c...
c...   1    xxx                   4*xxy-yyy-yyz
c...   2    xxy                   4*xxz-xxx-xzz
c...   3    xxz                   4*yyx-xxx-xzz
c...   4    xyy                   4*yyz-xxz-zzz
c...   5    xyz                   4*zzx-xxx-yyx
c...   6    xzz                   4*zzy-xxy-yyy
c...   7    yyy                    sqrt(40)*xyz
c...   8    yyz
c...   9    yzz
c...  10    zzz
c...
c...               g  functions
c...
c...   1   xxxx                sqrt(32/3)*(xxxx+yyyy-6*xxyy)
c...   2   xxxx                sqrt(32/3)*(xxxx+zzzz-6*xxzz)
c...   3   xxxx                sqrt(32/3)*(yyyy+zzzz-6*yyzz)
c...   4   xxxx                   sqrt(4/3)*(xxxy-3*xyzz)
c...   5   xxxx                   sqrt(4/3)*(xyyy-3*xyzz)
c...   6   xxxx                   sqrt(4/3)*(xxxz-3*xyyz)
c...   7   xxxx                   sqrt(4/3)*(xzzz-3*xyyz)
c...   8   xxxx                   sqrt(4/3)*(yyyz-3*xxyz)
c...   9   xxxx                   sqrt(4/3)*(yzzz-3*xxyz)
c...  10   xxxx
c...  11   xxxx
c...  12   xxxx
c...  13   xxxx
c...  14   xxxx
c...  15   xxxx
c...
c...               h  functions
c...
c...   1  xxxxx           4/sqrt(30)*(xxxxx-10*xxxyy+5*xyyyy)
c...   2  xxxxy           4/sqrt(30)*(xxxxx-10*xxxzz+5*xzzzz)
c...   3  xxxxz           4/sqrt(30)*(yyyyy-10*xxyyy+5*xxxxy)
c...   4  xxxyy           4/sqrt(30)*(yyyyy-10*yyyzz+5*yzzzz)
c...   5  xxxyz           4/sqrt(30)*(zzzzz-10*xxzzz+5*xxxxz)
c...   6  xxxzz           4/sqrt(30)*(zzzzz-10*yyzzz+5*yyyyz)
c...   7  xxyyy             4/sqrt(3)*(xxxxz+yyyyz-6*xxyyz)
c...   8  xxyyz             4/sqrt(3)*(xxxxy+yzzzz-6*xxyzz)
c...   9  xxyzz             4/sqrt(3)*(xxxxy+xzzzz-6*xyyzz)
c...  10  xxzzz                 16/sqrt(3)*(xxxyz-xyyyz)
c...  11  xyyyy                 16/3*(xxxyz-2*xyzzz+xyyyz)
c...  12  xyyyz
c...  13  xyyzz
c...  14  xyzzz
c...  15  xzzzz
c...  16  yyyyy
c...  17  yyyyz
c...  18  yyyzz
c...  19  yyzzz
c...  20  yzzzz
c...  21  zzzzz
c...
      dimension xc(ndim,nga,ngb)
      dimension xs(ndim,nga,ngb)
      logical transa,transb
      parameter (two=2.0d0,half=0.5d0,sqtw=2.886751345948129d-1)
      parameter (four=4.0d0,sqr40=6.324555320336759d0)
      parameter (three=3.0d0, six=6.0d0)
      parameter (yn=3.2659863237109d0,xn=1.15470053837925d0)
      parameter (five=5.0d0,ten=10.0d0)
      parameter (xnh=0.730296743340221d0,ynh=2.3094010767585d0,
     1           znh=9.23760430703401d0,snh=16.0d0/3.0d0)
c...
c...  first index
c...
      if(transa)then
        if(itypa.eq.4)then
          do j=1,ngb
          do ic=1,ndim
            xs(ic,1,j)=sqtw*(2*xc(ic,3,j)-xc(ic,1,j)-xc(ic,2,j))
            xs(ic,2,j)=half*(xc(ic,1,j)-xc(ic,2,j))
            xs(ic,3,j)=xc(ic,4,j)
            xs(ic,4,j)=xc(ic,5,j)
            xs(ic,5,j)=xc(ic,6,j)
          enddo
          enddo
        else if(itypa.eq.6)then
          do j=1,ngb
          do ic=1,ndim
            xs(ic,1,j)=four*xc(ic,2,j)-xc(ic,7,j)-xc(ic,9,j)
            xs(ic,2,j)=four*xc(ic,3,j)-xc(ic,8,j)-xc(ic,10,j)
            xs(ic,3,j)=four*xc(ic,4,j)-xc(ic,1,j)-xc(ic,6,j)
            xs(ic,4,j)=four*xc(ic,8,j)-xc(ic,3,j)-xc(ic,10,j)
            xs(ic,5,j)=four*xc(ic,6,j)-xc(ic,1,j)-xc(ic,4,j)
            xs(ic,6,j)=four*xc(ic,9,j)-xc(ic,2,j)-xc(ic,7,j)
            xs(ic,7,j)=sqr40*xc(ic,5,j)
          enddo
          enddo
        else if(itypa.eq.11)then
          do j=1,ngb
          do ic=1,ndim
            xs(ic,1,j)=xn*(xc(ic,1,j)+xc(ic,11,j)-six*xc(ic,4,j))
            xs(ic,2,j)=xn*(xc(ic,1,j)+xc(ic,15,j)-six*xc(ic,6,j))
            xs(ic,3,j)=xn*(xc(ic,11,j)+xc(ic,15,j)-six*xc(ic,13,j))
            xs(ic,4,j)=yn*(xc(ic,2,j)-three*xc(ic,9,j))
            xs(ic,5,j)=yn*(xc(ic,7,j)-three*xc(ic,9,j))
            xs(ic,6,j)=yn*(xc(ic,3,j)-three*xc(ic,8,j))
            xs(ic,7,j)=yn*(xc(ic,10,j)-three*xc(ic,8,j))
            xs(ic,8,j)=yn*(xc(ic,12,j)-three*xc(ic,5,j))
            xs(ic,9,j)=yn*(xc(ic,14,j)-three*xc(ic,5,j))
          enddo
          enddo
        else if(itypa.eq.12)then
          do j=1,ngb
          do ic=1,ndim
           xs(ic,1,j)=xnh*(xc(ic,1,j)-ten*xc(ic,4,j)+five*xc(ic,11,j))
           xs(ic,2,j)=xnh*(xc(ic,1,j)-ten*xc(ic,6,j)+five*xc(ic,15,j))
           xs(ic,3,j)=xnh*(xc(ic,16,j)-ten*xc(ic,7,j)+five*xc(ic,2,j))
           xs(ic,4,j)=xnh*(xc(ic,16,j)-ten*xc(ic,18,j)+five*xc(ic,20,j))
           xs(ic,5,j)=xnh*(xc(ic,21,j)-ten*xc(ic,10,j)+five*xc(ic,3,j))
           xs(ic,6,j)=xnh*(xc(ic,21,j)-ten*xc(ic,19,j)+five*xc(ic,17,j))
           xs(ic,7,j)=ynh*(xc(ic,3,j)+xc(ic,17,j)-six*xc(ic,8,j))
           xs(ic,8,j)=ynh*(xc(ic,2,j)+xc(ic,20,j)-six*xc(ic,9,j))
           xs(ic,9,j)=ynh*(xc(ic,11,j)+xc(ic,15,j)-six*xc(ic,13,j))
           xs(ic,10,j)=znh*(xc(ic,5,j)-xc(ic,12,j))
           xs(ic,11,j)=snh*(xc(ic,5,j)-two*xc(ic,14,j)+xc(ic,12,j))
          enddo
          enddo
        else if(itypa.eq.13)then
        endif
        if(itypa.ne.13)then
          do i=1,nsha
            do j=1,ngb
            do ic=1,ndim
              xc(ic,i,j)=xs(ic,i,j)
            enddo
            enddo
          enddo
        endif
      endif
c...
c...  second index
c...
      if(transb)then
        if(itypb.eq.4)then
          do i=1,nsha
          do ic=1,ndim
            xs(ic,i,1)=sqtw*(2*xc(ic,i,3)-xc(ic,i,1)-xc(ic,i,2))
            xs(ic,i,2)=half*(xc(ic,i,1)-xc(ic,i,2))
            xs(ic,i,3)=xc(ic,i,4)
            xs(ic,i,4)=xc(ic,i,5)
            xs(ic,i,5)=xc(ic,i,6)
          enddo
          enddo
        else if(itypb.eq.6)then
          do i=1,nsha
          do ic=1,ndim
            xs(ic,i,1)=four*xc(ic,i,2)-xc(ic,i,7)-xc(ic,i,9)
            xs(ic,i,2)=four*xc(ic,i,3)-xc(ic,i,8)-xc(ic,i,10)
            xs(ic,i,3)=four*xc(ic,i,4)-xc(ic,i,1)-xc(ic,i,6)
            xs(ic,i,4)=four*xc(ic,i,8)-xc(ic,i,3)-xc(ic,i,10)
            xs(ic,i,5)=four*xc(ic,i,6)-xc(ic,i,1)-xc(ic,i,4)
            xs(ic,i,6)=four*xc(ic,i,9)-xc(ic,i,2)-xc(ic,i,7)
            xs(ic,i,7)=sqr40*xc(ic,i,5)
          enddo
          enddo
        else if(itypb.eq.11)then
          do i=1,nsha
          do ic=1,ndim
            xs(ic,i,1)=xn*(xc(ic,i,1)+xc(ic,i,11)-six*xc(ic,i,4))
            xs(ic,i,2)=xn*(xc(ic,i,1)+xc(ic,i,15)-six*xc(ic,i,6))
            xs(ic,i,3)=xn*(xc(ic,i,11)+xc(ic,i,15)-six*xc(ic,i,13))
            xs(ic,i,4)=yn*(xc(ic,i,2)-three*xc(ic,i,9))
            xs(ic,i,5)=yn*(xc(ic,i,7)-three*xc(ic,i,9))
            xs(ic,i,6)=yn*(xc(ic,i,3)-three*xc(ic,i,8))
            xs(ic,i,7)=yn*(xc(ic,i,10)-three*xc(ic,i,8))
            xs(ic,i,8)=yn*(xc(ic,i,12)-three*xc(ic,i,5))
            xs(ic,i,9)=yn*(xc(ic,i,14)-three*xc(ic,i,5))
          enddo
          enddo
        else if(itypb.eq.12)then
          do i=1,nsha
          do ic=1,ndim
           xs(ic,i,1)=xnh*(xc(ic,i,1)-ten*xc(ic,i,4)+five*xc(ic,i,11))
           xs(ic,i,2)=xnh*(xc(ic,i,1)-ten*xc(ic,i,6)+five*xc(ic,i,15))
           xs(ic,i,3)=xnh*(xc(ic,i,16)-ten*xc(ic,i,7)+five*xc(ic,i,2))
           xs(ic,i,4)=xnh*(xc(ic,i,16)-ten*xc(ic,i,18)+five*xc(ic,i,20))
           xs(ic,i,5)=xnh*(xc(ic,i,21)-ten*xc(ic,i,10)+five*xc(ic,i,3))
           xs(ic,i,6)=xnh*(xc(ic,i,21)-ten*xc(ic,i,19)+five*xc(ic,i,17))
           xs(ic,i,7)=ynh*(xc(ic,i,3)+xc(ic,i,17)-six*xc(ic,i,8))
           xs(ic,i,8)=ynh*(xc(ic,i,2)+xc(ic,i,20)-six*xc(ic,i,9))
           xs(ic,i,9)=ynh*(xc(ic,i,11)+xc(ic,i,15)-six*xc(ic,i,13))
           xs(ic,i,10)=znh*(xc(ic,i,5)-xc(ic,i,12))
           xs(ic,i,11)=snh*(xc(ic,i,5)-two*xc(ic,i,14)+xc(ic,i,12))
          enddo
          enddo
        else if(itypb.eq.13)then
        endif
        if(itypb.ne.13)then
          do i=1,nsha
            do j=1,nshb
            do ic=1,ndim
              xc(ic,i,j)=xs(ic,i,j)
            enddo
            enddo
          enddo
        endif
      endif
      end
c=======================================================================
      subroutine psp_type12i(
     $                    basdat,   lla,    mma,   nna,   llb,
     $                       mmb,   nnb,    nga,   ngb,  lpsp,
     $                     nrpsp, alpha,     cc,  len0,   len,
     $                       cax,   cay,    caz,   cbx,   cby,
     $                       cbz,  lima,   limb,  mxpw, nmaxa,
     $                     nmaxb,  mijk,    yra,   yrb,  xyzi,
     $                         q,    ar,    tol,argmax,  ilen,
     $                       ifi,   ca2,   jlen,   jfi,   cb2,
     $                       pspi, time)
      implicit real*8 (a-h,o-z)
c...
c...   MM June-September 2004
c...
c...   pseudopotential integrals over primitives
c...
c...   Input parameters
c...
c...   basdat           basis set data
c...   lla,mma,nna,     arrays of angular momenta for A
c...   llb,mmb,nnb      arrays of angular momenta for B
c...   nga,ngb          number of Gaussian components for A and B
c...                    1=S, 3=P, 4=L, 6=D, 10=F, ...
c...   lpsp             maximum angular momentum of the pseudopotential
c...   nrpsp            array of powers of r of the psudopotential
c...   alpha            array of exponents of the pseudopotential
c...   cc               array of contraction coefficients of psp
c...   len0             contraction length of local part of psp
c...   len              array of contraction lengths for psp components
c...   cax,cay,caz      powers of CA
c...   cbx,cby,cbz      powers of CB
c...   lima,limb        do loop limits
c...   mxpw             maximum power of CA or CB
c...   nmaxa,nmaxb      maximum values of na and nb
c...   mijk             maximum exponent in angular integrals
c...   yra,yrb          scratch areas for real spherical harmonics
c...   xyzi             table of precomputed values for angular
c...                    integrals
c...   q                scratch area for radial integrals
c...   ar               scratch area for sums of angular and
c...                    radial integrals
c...   tol              tolerance threshold
c...   argmax           maximum argument of the exponential term
c...   ilen             contraction length of shell i
c...   ifi              pointer to first primitive of shell i
c...   ca2              distance CA squared
c...   jlen             contraction length of shell j
c...   jfi              pointer to first primitive of shell j
c...   cb2              distance CB squared
c...
c...   Output prameter
c...
c...   pspi             array of pseudopotential integrals
c...
      dimension basdat(13,*)
      dimension time(*)
c...
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      logical compute
c...
c...  loop over primitives
c...
      maxn=nmaxa+nmaxb
      do ip=1,ilen
        ipi=ifi+ip
        alpa=basdat(1,ipi)     ! exponent of shell a
        do jp=1,jlen
          jpj=jfi+jp
          alpb=basdat(1,jpj)     ! exponent of shell b
          exparg=alpa*ca2+alpb*cb2
          if(exparg.le.argmax)then ! test over the exponential
c...
c...  zero out ar array
c...
c           call secund(t0)
            call zeroit(ar,(nmaxa+1)**3*(nmaxb+1)**3)
            compute=.false.
c...
c...  compute radial and angular integrals
c...
            call psp_ar12(
     $                alpa,    one,  alpb,     one,  lpsp,
     $               nrpsp,  alpha,    cc,    len0,   len,
     $                 cax,    cay,    caz,    cbx,   cby,
     $                 cbz,   mxpw,  nmaxa,  nmaxb,  maxn,
     $                mijk,    yra,    yrb,   xyzi,     q,
     $              exparg,   tol,     ar,   compute,time)
c           call secund(t1)
c           time(1)=time(1)+t1-t0
c...
c...  assemble integrals over primitives
c...
            if(compute)then
c             call secund(tini)
              call psp_lg0(1,nga,1,ngb,
     $                     lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                     cax,cay,caz,cbx,cby,cbz,mxpw,
     $                     nmaxa,nmaxb,ar,ip,jp,ilen,jlen,nga,ngb,pspi)
c             call secund(tifi)
c             time(7)=time(7)+tifi-tini
            endif
          endif
        enddo
      enddo
      end
c=======================================================================
      subroutine psp_type12cab(
     $                    basdat,   lla,    mma,   nna,   llb,
     $                       mmb,   nnb,    nga,   ngb,  lpsp,
     $                     nrpsp, alpha,     cc,  len0,   len,
     $                       cax,   cay,    caz,   cbx,   cby,
     $                       cbz,  lima,   limb,  mxpw, nmaxa,
     $                     nmaxb,  mijk,    yra,   yrb,  xyzi,
     $                         q,    ar,    tol,argmax,  ilen,
     $                       ifi,   ca2,   jlen,   jfi,   cb2,
     $                       cna,   cnb,   pspc,  time)
      implicit real*8 (a-h,o-z)
c...
c...   MM June-September 2004
c...
c...   pseudopotential integrals over contracted
c...
c...   Input parameters
c...
c...   basdat           basis set data
c...   lla,mma,nna,     arrays of angular momenta for A
c...   llb,mmb,nnb      arrays of angular momenta for B
c...   nga,ngb          number of Gaussian components for A and B
c...                    1=S, 3=P, 4=L, 6=D, 10=F, ...
c...   lpsp             maximum angular momentum of the pseudopotential
c...   nrpsp            array of powers of r of the psudopotential
c...   alpha            array of exponents of the pseudopotential
c...   cc               array of contraction coefficients of psp
c...   len0             contraction length of local part of psp
c...   len              array of contraction lengths for psp components
c...   cax,cay,caz      powers of CA
c...   cbx,cby,cbz      powers of CB
c...   lima,limb        do loop limits
c...   mxpw             maximum power of CA or CB
c...   nmaxa,nmaxb      maximum values of na and nb
c...   maxn             maximum value of n
c...   mijk             maximum exponent in angular integrals
c...   yra,yrb          scratch areas for real spherical harmonics
c...   xyzi             table of precomputed values for angular
c...                    integrals
c...   q                scratch area for radial integrals
c...   ar               scratch area for sums of angular and
c...                    radial integrals
c...   tol              tolerance threshold
c...   argmax           maximum argument of the exponential term
c...   ilen             contraction length of shell i
c...   ifi              pointer to first primitive of shell i
c...   ca2              distance CA squared
c...   jlen             contraction length of shell j
c...   jfi              pointer to first primitive of shell j
c...   cb2              distance CB squared
c...   cna              array for normalization/contraction of shell a
c...   cnb              array for normalization/contraction of shell b
c...
c...   Output prameter
c...
c...   pspc              array of pseudopotential integrals
c...
      dimension basdat(13,*),cna(*),cnb(*)
      dimension time(*)
c...
      parameter (zero=0.0d0,two=2.0d0)
      logical compute
c...
c...  zero out ar array
c...
c...
c         call secund(t0)
          call zeroit(ar,(nmaxa+1)**3*(nmaxb+1)**3)
          compute=.false.
c...
c...  loop over primitives
c...
          maxn=nmaxa+nmaxb
          do ip=1,ilen
            ipi=ifi+ip
            alpa=basdat(1,ipi)          ! exponent of shell a
            ccna=cna(ip)                ! coefficient for shell a
            do jp=1,jlen
              jpj=jfi+jp
              alpb=basdat(1,jpj)        ! exponent of shell b
              ccnb=cnb(jp)              ! coefficient for shell b
              exparg=alpa*ca2+alpb*cb2
              if(exparg.le.argmax)then ! test over the exponential
              call psp_ar12(
     $                alpa,   ccna,  alpb,    ccnb,  lpsp,
     $               nrpsp,  alpha,    cc,    len0,   len,
     $                 cax,    cay,    caz,    cbx,   cby,
     $                 cbz,   mxpw,  nmaxa,  nmaxb,  maxn,
     $                mijk,    yra,    yrb,   xyzi,     q,
     $              exparg,   tol,     ar, compute,  time)
              endif
            enddo
          enddo
c         call secund(t1)
c         time(1)=time(1)+t1-t0
c...
c...  assemble integrals over contracted
c...
          if(compute)then
c         call secund(tini)
          call psp_lg0(1,nga,1,ngb,
     $                 lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                 cax,cay,caz,cbx,cby,cbz,mxpw,
     $                 nmaxa,nmaxb,ar,1,1,1,1,nga,ngb,pspc)
c         call secund(tifi)
c         time(7)=time(7)+tifi-tini
          endif
      end
c=======================================================================
      subroutine psp_type12cb(
     $                    basdat,   lla,    mma,   nna,    llb,
     $                       mmb,   nnb,    nga,   ngb,   lpsp,
     $                     nrpsp, alpha,     cc,  len0,    len,
     $                       cax,   cay,    caz,   cbx,    cby,
     $                       cbz,  lima,   limb,  mxpw,  nmaxa,
     $                     nmaxb,  mijk,    yra,   yrb,   xyzi,
     $                         q,    ar,    tol,argmax,   ilen,
     $                       ifi,   ca2,   jlen,   jfi,    cb2,
     $                       cna,   cnb,   pspi,  time)
      implicit real*8 (a-h,o-z)
c...
c...   MM June-September 2004
c...
c...   pseudopotential integrals over primitives of shell a and
c...   contracted of shell b
c...
c...   Input parameters
c...
c...   alpa,alphb       exponents of Gaussians A and B
c...   lla,mma,nna,     arrays of angular momenta for A
c...   llb,mmb,nnb      arrays of angular momenta for B
c...   nga,ngb          number of Gaussian components for A and B
c...                    1=S, 3=P, 4=L, 6=D, 10=F, ...
c...   lpsp             maximum angular momentum of the pseudopotential
c...   nrpsp            array of powers of r of the psudopotential
c...   alpha            array of exponents of the pseudopotential
c...   cc               array of contraction coefficients of psp
c...   len0             contraction length of local part of psp
c...   len              array of contraction lengths for psp components
c...   cax,cay,caz      powers of CA
c...   cbx,cby,cbz      powers of CB
c...   lima,limb        do loop limits
c...   mxpw             maximum power of CA or CB
c...   nmaxa,nmaxb      maximum values of na and nb
c...   maxn             maximum value of n
c...   mijk             maximum exponent in angular integrals
c...   yra,yrb          scratch areas for real spherical harmonics
c...   xyzi             table of precomputed values for angular
c...                    integrals
c...   q                scratch area for radial integrals
c...   ar               scratch area for sums of angular and
c...                    radial integrals
c...   tol              tolerance threshold
c...   argmax           maximum argument of the exponential term
c...   ilen             contraction length of shell i
c...   ifi              pointer to first primitive of shell i
c...   ca2              distance CA squared
c...   jlen             contraction length of shell j
c...   jfi              pointer to first primitive of shell j
c...   cb2              distance CB squared
c...   cna              array for normalization/contraction of shell a
c...   cnb              array for normalization/contraction of shell b
c...
c...   Output prameter
c...
c...   pspi             array of pseudopotential integrals
c...
      dimension basdat(13,*),cna(*),cnb(*)
      dimension time(*)
c...
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      logical compute
c...
c...  loop over primitives
c...
          maxn=nmaxa+nmaxb
          do ip=1,ilen
            ipi=ifi+ip
            alpa=basdat(1,ipi)     ! exponent of shell a
c...
c...  zero out ar array
c...
c...
c         call secund(t0)
            call zeroit(ar,(nmaxa+1)**3*(nmaxb+1)**3)
            compute=.false.
            do jp=1,jlen
              jpj=jfi+jp
              alpb=basdat(1,jpj)     ! exponent of shell b
              ccnb=cnb(jp)           ! coefficient for shell b
              exparg=alpa*ca2+alpb*cb2
              if(exparg.le.argmax)then ! test over the exponential
              loops=loops+1
              call psp_ar12(
     $                alpa,    one,  alpb,    ccnb,  lpsp,
     $               nrpsp,  alpha,    cc,    len0,   len,
     $                 cax,    cay,    caz,    cbx,   cby,
     $                 cbz,   mxpw,  nmaxa,  nmaxb,  maxn,
     $                mijk,    yra,    yrb,   xyzi,     q,
     $              exparg,   tol,     ar, compute,  time)
              endif
            enddo
c         call secund(t1)
c         time(1)=time(1)+t1-t0
c...
c...  assemble integrals over contracted
c...
          if(compute)then
c         call secund(tini)
          call psp_lg0(1,nga,1,ngb,
     $                 lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                 cax,cay,caz,cbx,cby,cbz,mxpw,
     $                 nmaxa,nmaxb,ar,ip,1,ilen,1,nga,ngb,pspi)
c         call secund(tifi)
c         time(7)=time(7)+tifi-tini
          endif
          enddo
      end
c=======================================================================
      subroutine psp_type12ca(
     $                    basdat,   lla,    mma,   nna,    llb,
     $                       mmb,   nnb,    nga,   ngb,   lpsp,
     $                     nrpsp, alpha,     cc,  len0,    len,
     $                       cax,   cay,    caz,   cbx,    cby,
     $                       cbz,  lima,   limb,  mxpw,  nmaxa,
     $                     nmaxb,  mijk,    yra,   yrb,   xyzi,
     $                         q,    ar,    tol,argmax,   ilen,
     $                       ifi,   ca2,   jlen,   jfi,    cb2,
     $                       cna,   cnb,   pspi,  time)
      implicit real*8 (a-h,o-z)
c...
c...   MM June-September 2004
c...
c...   pseudopotential integrals
c...
c...   Input parameters
c...
c...   alpa,alphb       exponents of Gaussians A and B
c...   lla,mma,nna,     arrays of angular momenta for A
c...   llb,mmb,nnb      arrays of angular momenta for B
c...   nga,ngb          number of Gaussian components for A and B
c...                    1=S, 3=P, 4=L, 6=D, 10=F, ...
c...   lpsp             maximum angular momentum of the pseudopotential
c...   nrpsp            array of powers of r of the psudopotential
c...   alpha            array of exponents of the pseudopotential
c...   cc               array of contraction coefficients of psp
c...   len0             contraction length of local part of psp
c...   len              array of contraction lengths for psp components
c...   cax,cay,caz      powers of CA
c...   cbx,cby,cbz      powers of CB
c...   lima,limb        do loop limits
c...   mxpw             maximum power of CA or CB
c...   nmaxa,nmaxb      maximum values of na and nb
c...   maxn             maximum value of n
c...   mijk             maximum exponent in angular integrals
c...   yra,yrb          scratch areas for real spherical harmonics
c...   xyzi             table of precomputed values for angular
c...                    integrals
c...   q                scratch area for radial integrals
c...   ar               scratch area for sums of angular and
c...                    radial integrals
c...   tol              tolerance threshold
c...   argmax           maximum argument of the exponential term
c...   ilen             contraction length of shell i
c...   ifi              pointer to first primitive of shell i
c...   ca2              distance CA squared
c...   jlen             contraction length of shell j
c...   jfi              pointer to first primitive of shell j
c...   cb2              distance CB squared
c...   cna              array for normalization/contraction of shell a
c...   cnb              array for normalization/contraction of shell b
c...
c...   Output prameter
c...
c...   pspi             array of pseudopotential integrals
c...
      dimension basdat(13,*),cna(*),cnb(*)
      dimension time(*)
c...
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      logical compute
c...
c...  loop over primitives
c...
        maxn=nmaxa+nmaxb
        do jp=1,jlen
          jpj=jfi+jp
          alpb=basdat(1,jpj)     ! exponent of shell b
c...
c...  zero out ar array
c...
c...
c         call secund(t0)
          call zeroit(ar,(nmaxa+1)**3*(nmaxb+1)**3)
          compute=.false.
          do ip=1,ilen
            ipi=ifi+ip
            alpa=basdat(1,ipi)     ! exponent of shell a
            ccna=cna(ip)           ! coefficient for shell a
            exparg=alpa*ca2+alpb*cb2
            if(exparg.le.argmax)then ! test over the exponential
              loops=loops+1
              call psp_ar12(
     $                alpa,   ccna,  alpb,     one,  lpsp,
     $               nrpsp,  alpha,    cc,    len0,   len,
     $                 cax,    cay,    caz,    cbx,   cby,
     $                 cbz,   mxpw,  nmaxa,  nmaxb,  maxn,
     $                mijk,    yra,    yrb,   xyzi,     q,
     $              exparg,   tol,     ar, compute,  time)
            endif
          enddo
c         call secund(t1)
c         time(1)=time(1)+t1-t0
c...
c...  assemble integrals over contracted
c...
          if(compute)then
c         call secund(tini)
          call psp_lg0(1,nga,1,ngb,
     $                 lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                 cax,cay,caz,cbx,cby,cbz,mxpw,
     $                 nmaxa,nmaxb,ar,1,jp,1,jlen,nga,ngb,pspi)
c         call secund(tifi)
c         time(7)=time(7)+tifi-tini
          endif
        enddo
      end
c=======================================================================
      subroutine psp_lg0(igas,igaf,igbs,igbf,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,i1,i2,nd1,nd2,nga,ngb,psp)
      implicit real*8 (a-h,o-z)
c...
c...  loop over the Gaussian component and assembling
c...  of the zero-order pseudopotential integrals
c...
      dimension lla(nga),mma(nga),nna(nga),llb(ngb),mmb(ngb),nnb(ngb)
      dimension lima(3),limb(3)
      dimension psp(nd1,nd2,nga,ngb)
c...
c...  loop over the Gaussian components
c...
      do iga=igas,igaf
        la=lla(iga)
        laf=min(la,lima(1))
        ma=mma(iga)
        maf=min(ma,lima(2))
        na=nna(iga)
        naf=min(na,lima(3))
        do igb=igbs,igbf
          lb=llb(igb)
          lbf=min(lb,limb(1))
          mb=mmb(igb)
          mbf=min(mb,limb(2))
          nb=nnb(igb)
          nbf=min(nb,limb(3))
c...
c...  summation
c...
          call psp_suml2(la,laf,ma,maf,na,naf,lb,lbf,mb,mbf,nb,nbf,
     $                  cax,cay,caz,cbx,cby,cbz,mxpw,
     $                nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppint)
c...
          psp(i1,i2,iga,igb)=psp(i1,i2,iga,igb)+ppint
        enddo
      enddo
      end
c=======================================================================
      subroutine psp_ar12(
     $                      alpa,  ccna,  alpb,  ccnb,  lpsp,
     $                     nrpsp, alpha,    cc,  len0,   len,
     $                       cax,   cay,   caz,   cbx,   cby,
     $                       cbz,  mxpw, nmaxa, nmaxb,  maxn,
     $                      mijk,   yra,   yrb,  xyzi,     q,
     $                    exparg,   tol,    ar,  compute,time)
      implicit real*8 (a-h,o-z)
c...
c...   MM June-September 2004
c...
c...   This subroutine computes a table of sums of angular and
c...   radial integral, accumulating type one and type 2 integrals
c...
c...   Input parameters
c...
c...   alpa             exponent of Gaussians A
c...   ccna             coefficient for Gaussian A
c...   alpa,alphb       exponents of Gaussian B
c...   ccnb             coefficient for Gaussian B
c...   lpsp             maximum angular momentum of the pseudopotential
c...   nrpsp            array of powers of r of the psudopotential
c...   alpha            array of exponents of the pseudopotential
c...   cc               array of contraction coefficients of psp
c...   len0             contraction length of local part of psp
c...   len              array of contraction lengths for psp components
c...   cax,cay,caz      powers of CA
c...   cbx,cby,cbz      powers of CB
c...   mxpw             maximum power of CA or CB
c...   nmaxa,nmaxb      maximum values of na and nb
c...   maxn             maximum value of n
c...   mijk             maximum exponent in angular integrals
c...   yra,yrb          scratch areas for real spherical harmonics
c...   xyzi             table of precomputed values for angular
c...                    integrals
c...   q                scratch area for radial integrals
c...   exparg           argument of the exponential term
c...   tol              tolerance threshold
c...
c...   Output parameter
c...
c...   compute          flag set to true if integrals have been computed
c...   ar               table of sums of angular and radial integrals
c...
      dimension alpha(*),cc(*),nrpsp(*),len(0:lpsp)
      dimension cax(0:mxpw),cay(0:mxpw),caz(0:mxpw)
      dimension cbx(0:mxpw),cby(0:mxpw),cbz(0:mxpw)
      dimension q(0:maxn,0:maxn)
      dimension time(20)
      real*8 k,kx,ky,kz
      real*8 ka,kax,kay,kaz,kb,kbx,kby,kbz
c...
c...   threshold under which we consider k to be zero
c...
      parameter (epsik=1.0d-12)
c...
      parameter (zero=0.0d0,two=2.0d0)
      parameter (pi=3.141592653589793d0,pim=1.570796326794897d0,
     $           tmpi=4.71238898038469d0)
      dimension ar(0:nmaxa,0:nmaxa,0:nmaxa,0:nmaxb,0:nmaxb,0:nmaxb)
      logical doint,compute
c...
c...  logarithm of the threshold
c...
      tolg=LOG(tol)
      aab=alpa+alpb
c...
c...  multiplicative factor for type 1 integrals.
c...
c...   f1=4/sqrt(%PI)*exp(-alpa*CA^2-alpb*CB^2)*ccna*ccnb
c...
c       call getival('iout',iout)
c       write(iout,*)'tol',tol
        f1=2.256758334191025D0*exp(-exparg)*ccna*ccnb
c       if(debug)write(iout,*)'exparg,f1',exparg,f1
c...
c...  type one integrals
c...
      if(.not.(len0.eq.1.and.cc(1).eq.zero))then
c...
c...  compute k
c...
        kx=-two*(alpa*cax(1)+alpb*cbx(1))
        ky=-two*(alpa*cay(1)+alpb*cby(1))
        kz=-two*(alpa*caz(1)+alpb*cbz(1))
        k=sqrt(kx**2+ky**2+kz**2)
c       if(debug)write(iout,*)'kx,ky,kz,k',kx,ky,kz,k
        if(k.lt.epsik)k=zero
c...
c...  global test for type 1 integrals
c...
        k2=k*k
        doint=.false.
        do nc=1,len0
          if(0.25d0*k2/(aab+alpha(nc))-exparg.gt.tolg)doint=.true.
        enddo
        if(doint)then
          compute=.true.
c...
c...  compute type 1 radial integrals
c...
c       call secund(t0)
        call domrad1(maxn,nrpsp,k,alpa+alpb,alpha,cc,len0,f1,exparg,
     $               tolg,maxn,q)
c       call secund(t1)
c       time(11)=time(11)+t1-t0
c...
c...  now compute sums of type 1 angular and radial integrals
c...
          if(k.ne.zero)then
c...
c...  spherical polar angle theta
c...
            theta=acos(kz/k)
            if(sqrt(kx**2+ky**2).ge.epsik)then
c...
c...  spherical polar angle phi can be defined
c...
c            call secund(t0)
               if(kx**2.lt.epsik)then
                 if(ky.gt.zero)then
                   phi=pim
                 else
                   phi=tmpi
                 endif
               else if(ky**2.lt.epsik)then
                 if(kx.gt.zero)then
                   phi=zero
                 else
                   phi=pi
                 endif
               else if(kx.gt.zero.and.ky.gt.zero)then
                 phi=atan(ky/kx)
               else if(kx.lt.zero.and.ky.gt.zero)then
                 phi=pi+atan(ky/kx)
               else if(kx.lt.zero.and.ky.lt.zero)then
                 phi=pi+atan(ky/kx)
               else if(kx.gt.zero.and.ky.lt.zero)then
                 phi=pi+pi+atan(ky/kx)
               endif
               call yrleg(maxn,theta,phi,yra)
c            call secund(t1)
               call sumar1(q,yra,xyzi,nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,
     $                     nmaxb,maxn,maxn,maxn,1,maxn,mijk,ar)
c            call secund(t2)
c            time(12)=time(12)+t1-t0
c            time(13)=time(13)+t2-t1
            else
c...
c...  spherical polar angle phi cannot be defined
c...
c           call secund(t0)
              call zeroit(yra,(maxn+1)**2)
            call yrleg0(maxn,theta,yra)
c           call secund(t1)
              call sumar1(q,yra,xyzi,nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,
     $                    nmaxb,maxn,maxn,maxn,0,maxn,mijk,ar)
c            call secund(t2)
c            time(12)=time(12)+t1-t0
c            time(13)=time(13)+t2-t1
            endif
          else
c           call secund(t0)
            call sumar0(q,xyzi,nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,
     $                  maxn,maxn,maxn,maxn,mijk,ar)
c           call secund(t1)
c            time(14)=time(14)+t1-t0
          endif
        endif
      endif
c...
c...  type 2 integrals
c...
c...
c...  compute ka and kb and the real spherical harmonics
c...
c           call secund(t0)
      kax=-two*alpa*cax(1)
      kay=-two*alpa*cay(1)
      kaz=-two*alpa*caz(1)
      ka=sqrt(kax**2+kay**2+kaz**2)
c       if(debug)write(iout,*)'kax,kay,kaz,ka',kax,kay,kaz,ka
      if(ka.lt.epsik) then
         ka=zero
         ika=0
      else
        thetaa=acos(kaz/ka)
        if(sqrt(kax**2+kay**2).ge.epsik)then
           if(kax**2.lt.epsik)then
             if(kay.gt.zero)then
               phia=pim
             else
               phia=tmpi
             endif
           else if(kay**2.lt.epsik)then
             if(kax.gt.zero)then
               phia=zero
             else
               phia=pi
             endif
           else if(kax.gt.zero.and.kay.gt.zero)then
             phia=atan(kay/kax)
           else if(kax.lt.zero.and.kay.gt.zero)then
             phia=pi+atan(kay/kax)
           else if(kax.lt.zero.and.kay.lt.zero)then
             phia=pi+atan(kay/kax)
           else if(kax.gt.zero.and.kay.lt.zero)then
             phia=pi+pi+atan(kay/kax)
           endif
           call yrleg(nmaxa+lpsp,thetaa,phia,yra)
           ika=2
        else
          call zeroit(yra,(nmaxa+lpsp+1)**2)
          call yrleg0(nmaxa+lpsp,thetaa,yra)
           ika=1
        endif
      endif
      kbx=-two*alpb*cbx(1)
      kby=-two*alpb*cby(1)
      kbz=-two*alpb*cbz(1)
      kb=sqrt(kbx**2+kby**2+kbz**2)
c       if(debug)write(iout,*)'kbx,kby,kbz,kb',kbx,kby,kbz,kb
      if(kb.lt.epsik) then
         kb=zero
         ikb=0
      else
        thetab=acos(kbz/kb)
        if(sqrt(kbx**2+kby**2).ge.epsik)then
           if(kbx**2.lt.epsik)then
             if(kby.gt.zero)then
               phib=pim
             else
               phib=tmpi
             endif
           else if(kby**2.lt.epsik)then
             if(kbx.gt.zero)then
               phib=zero
             else
               phib=pi
             endif
           else if(kbx.gt.zero.and.kby.gt.zero)then
             phib=atan(kby/kbx)
           else if(kbx.lt.zero.and.kby.gt.zero)then
             phib=pi+atan(kby/kbx)
           else if(kbx.lt.zero.and.kby.lt.zero)then
             phib=pi+atan(kby/kbx)
           else if(kbx.gt.zero.and.kby.lt.zero)then
             phib=pi+pi+atan(kby/kbx)
           endif
           call yrleg(nmaxb+lpsp,thetab,phib,yrb)
           ikb=2
        else
          call zeroit(yrb,(nmaxb+lpsp+1)**2)
          call yrleg0(nmaxb+lpsp,thetab,yrb)
          ikb=1
        endif
      endif
c           call secund(t1)
c            time(14)=time(14)+t1-t0
c...
c...  multiplicative factor for type 2 integrals.
c...
c...   4/sqrt(%PI)*4*%PI*ccna*ccnb
c...
      f2=28.35926161448825D0*ccna*ccnb
c...
c...  loop over the pseudopotential components
c...
      kab2=(ka+kb)*(ka+kb)
      argd=(ka*ka/(4.0d0*alpa)+kb*kb/(4.0d0*alpb))
      is=len0+1
      do lp=0,lpsp
        leni=len(lp)
        maxl=max(nmaxa,nmaxb)+lp
c...
c...  global test for type 2 integrals
c...
        doint=.false.
        do nc=1,leni
         if(0.25d0*kab2/(aab+alpha(is+nc-1))-argd.gt.tolg)doint=.true.
        enddo
        if(doint)then
          compute=.true.
c...
c...  compute type 2 radial integrals
c...
c           call secund(t0)
          call domrad2n(nmaxa,nmaxb,lp,nrpsp(is),ka,kb,alpa,alpb,
     $                 alpha(is),cc(is),leni,maxn,maxl,f2,tolg,q)
c           call secund(t1)
c...
c...  now accumulate the sums of type 2 angular and radial integrals
c...
          call sumar2n(q,yra,ika,yrb,ikb,xyzi,lp,nmaxa,nmaxa,nmaxa,
     $                nmaxb,nmaxb,nmaxb,nmaxa,nmaxb,maxn,maxl,mijk,ar)
c           call secund(t2)
c            time(15)=time(15)+t1-t0
c            time(16)=time(16)+t2-t1
        endif
        is=is+leni
      enddo
c...
      end
c=======================================================================
      subroutine psp_suml2(al,af,bl,bf,cl,cf,dl,df,el,ef,fl,ff,
     $                     cax,cay,caz,cbx,cby,cbz,mxpw,
     $                     am,bm,cm,dm,em,fm,ar2,ppint)
      implicit real*8 (a-h,o-z)
c...
c...  summation loop for pseudopotential integrals
c...
c...  al,af              limits for a
c...  bl,bf              limits for b
c...  cl,cf              limits for c
c...  dl,df              limits for d
c...  el,ef              limits for e
c...  fl,ff              limits for f
c...  cax,cay,caz        powers of distance CA
c...  cbx,cby,cbz        powers of distance CB
c...  mxpw               maximum power of distances
c...  am,bm,cm,dm,em,fm  dimensions of ar2 array
c...  ar2                table of sums of radial and angular integrals
c...  ppint              summation result
c...
      integer al,af,am,bl,bf,bm,cl,cf,cm,dl,df,dm,el,ef,em,fl,ff,fm
      dimension ar2(0:am,0:bm,0:cm,0:dm,0:em,0:fm)
      dimension cax(0:mxpw),cay(0:mxpw),caz(0:mxpw)
      dimension cbx(0:mxpw),cby(0:mxpw),cbz(0:mxpw)
      integer a,b,c,d,e,f
      parameter (zero=0.0d0)
c...
c...   table of binomial coefficients: BC(m,n)=binomial(m,n)
c...
      real*8 BC(0:16,0:16)
      save BC
      data BC/1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0
     $  D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,0.0D0,1.0D0,2.0D0,3
     $  .0D0,4.0D0,5.0D0,6.0D0,7.0D0,8.0D0,9.0D0,1.0D1,1.1D1,1.2D1,1.3D1
     $  ,1.4D1,1.5D1,1.6D1,0.0D0,0.0D0,1.0D0,3.0D0,6.0D0,1.0D1,1.5D1,2.1
     $  D1,2.8D1,3.6D1,4.5D1,5.5D1,6.6D1,7.8D1,9.1D1,1.05D2,1.2D2,0.0D0,
     $  0.0D0,0.0D0,1.0D0,4.0D0,1.0D1,2.0D1,3.5D1,5.6D1,8.4D1,1.2D2,1.65
     $  D2,2.2D2,2.86D2,3.64D2,4.55D2,5.6D2,0.0D0,0.0D0,0.0D0,0.0D0,1.0D
     $  0,5.0D0,1.5D1,3.5D1,7.0D1,1.26D2,2.1D2,3.3D2,4.95D2,7.15D2,1.001
     $  D3,1.365D3,1.82D3,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,6.0D0,2.1D
     $  1,5.6D1,1.26D2,2.52D2,4.62D2,7.92D2,1.287D3,2.002D3,3.003D3,4.36
     $  8D3,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,7.0D0,2.8D1,8.4D1,
     $  2.1D2,4.62D2,9.24D2,1.716D3,3.003D3,5.005D3,8.008D3,0.0D0,0.0D0,
     $  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,8.0D0,3.6D1,1.2D2,3.3D2,7.92
     $  D2,1.716D3,3.432D3,6.435D3,1.144D4,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
     $  ,0.0D0,0.0D0,0.0D0,1.0D0,9.0D0,4.5D1,1.65D2,4.95D2,1.287D3,3.003
     $  D3,6.435D3,1.287D4,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0
     $  D0,0.0D0,1.0D0,1.0D1,5.5D1,2.2D2,7.15D2,2.002D3,5.005D3,1.144D4,
     $  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D
     $  0,1.1D1,6.6D1,2.86D2,1.001D3,3.003D3,8.008D3,0.0D0,0.0D0,0.0D0,0
     $  .0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,1.2D1,7.8D1
     $  ,3.64D2,1.365D3,4.368D3,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D
     $  0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,1.3D1,9.1D1,4.55D2,1.82D3,
     $  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D
     $  0,0.0D0,0.0D0,1.0D0,1.4D1,1.05D2,5.6D2,0.0D0,0.0D0,0.0D0,0.0D0,0
     $  .0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0
     $  ,1.5D1,1.2D2,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0
     $  D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,1.6D1,0.0D0,0.0D0,0
     $  .0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
     $  ,0.0D0,0.0D0,0.0D0,1.0D0/
c...
      ppint=zero
c...
      do a=al,af,-1
        ax=BC(al,a)*cax(al-a)
        do b=bl,bf,-1
          abx=ax*BC(bl,b)*cay(bl-b)
          do c=cl,cf,-1
            abcx=abx*BC(cl,c)*caz(cl-c)
            do d=dl,df,-1
              abcdx=abcx*BC(dl,d)*cbx(dl-d)
              do e=el,ef,-1
                abcdex=abcdx*BC(el,e)*cby(el-e)
                do f=fl,ff,-1
                  ppint=ppint+abcdex*BC(fl,f)*cbz(fl-f)*ar2(a,b,c,d,e,f)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c...
      end
c======================================================================
      subroutine sumar1(q1,yr,xyzi,am,bm,cm,dm,em,fm,
     $                  im,jm,km,iphik,maxn,mijk,ar1)
      implicit real*8 (a-h,o-z)
c...
c...  MM june 2004
c...
c...  This subroutine computes sums of type one angular and radial
c...  integrals and returns them in ar1:
c...
c...                k + j + i
c...                ====
c...                \
c... ar1(i, j, k) =  >        OMEGA1(i, j, k, l) q1(k + j + i, l)
c...                /
c...                ====
c...                l = 0
c...
c...  where omega1 is a type one angular integral, and q1 is a type
c...  one radial integral. The summation has the additional
c...  constraint that (i+j+k)-l must be even.
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (16-19)
c...
c...  input parameters
c...
c...  q1            array of type 1 radial integrals
c...  yr            array of real spherical harmonics
c...  xyzi          table of double factorial products for angular
c...                integrals
c...  am,bm,cm,     maximum values of a,b and c
c...  dm,em,fm,     maximum values of d,e and f
c...  im,jm,km,     maximum values of i,j and k
c...  iphik         switch to tell which subroutine to use to compute
c...                the angular integrals omega1:
c...                       0      use omega0  (angle phik is undefined)
c...                       1      use omega1
c...  maxn          maximum allowed n value, for dimensioning
c...  mijk          maximum exponent for angular integrals,
c...                for dimensioning
c...
c...  output parameter
c...
c...  ar1    array of sums of angular and radial integrals
c...
      integer am,bm,cm,dm,em,fm
      dimension q1(0:maxn,0:maxn),ar1(0:am,0:bm,0:cm,0:dm,0:em,0:fm)
      integer a,b,c,d,e,f
      parameter (lmax=16)
      dimension omega(0:lmax)
      parameter (zero=0.0d0)
c...
      if(maxn.gt.lmax)
     $  call nerror(1,'Sumar1','maximum l value exceeded',maxn,lmax)
c...
      do i=0,im
        do j=0,jm
          do k=0,km
            ijk=i+j+k
            if(ijk.le.maxn)then
c...
c...  compute the batch of type 1 angular integrals
c...
              if(iphik.eq.1)then
                call omega1(i,j,k,yr,mijk,xyzi,omega)
              else
                call omega0(i,j,k,yr,mijk,xyzi,omega)
              endif
c...
c...  sum angular and radial integrals
c...
              sum=zero
              do l=mod(ijk,2),ijk,2
                sum=sum+omega(l)*q1(ijk,l)
              enddo
              if(sum.ne.zero)then
                do a=0,min(am,i)
                  d=i-a
                  if(d.le.dm)then
                  do b=0,min(bm,j)
                    e=j-b
                    if(e.le.em)then
                    do c=0,min(cm,k)
                      f=k-c
                      if(f.le.fm)then
                        ar1(a,b,c,d,e,f)=ar1(a,b,c,d,e,f)+sum
                      endif
                    enddo
                    endif
                  enddo
                  endif
                enddo
              endif
            endif
          enddo
        enddo
      enddo
      end
c======================================================================
      subroutine sumar0(q1,xyzi,am,bm,cm,dm,em,fm,
     $                  im,jm,km,maxn,mijk,ar1)
      implicit real*8 (a-h,o-z)
c...
c...  MM june 2004
c...
c...  This subroutine computes sums of type one angular and radial
c...  integrals and returns them in ar1:
c...
c... ar1(i, j, k) =   OMEGA1(i, j, k, 0) q1(k + j + i, 0)
c...
c...  where omega1 is a type one angular integral, and q1 is a type
c...  one radial integral.
c...
c...  this version is used when the vector k is zero, and thus only
c...  the l=0 terms of omega1 and q1 survive
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (16-19)
c...
c...  input parameters
c...
c...  q1            array of type 1 radial integrals
c...  xyzi          table of double factorial products for angular
c...                integrals
c...  am,bm,cm,     maximum values of a,b and c
c...  dm,em,fm,     maximum values of d,e and f
c...  im,jm,km,     maximum values of i,j and k
c...  maxn          maximum allowed n value, for dimensioning
c...  mijk          maximum exponent for angular integrals,
c...                for dimensioning
c...
c...  output parameter
c...
c...  ar1    array of sums of angular and radial integrals
c...
      integer am,bm,cm,dm,em,fm
      dimension q1(0:maxn,0:maxn),ar1(0:am,0:bm,0:cm,0:dm,0:em,0:fm)
      dimension xyzi(0:mijk,0:mijk,0:mijk)
      integer a,b,c,d,e,f
      parameter (zero=0.0d0)
c...
      do i=0,im
        do j=0,jm
          do k=0,km
            ijk=i+j+k
            if(ijk.le.maxn)then
c...
c...  compute the type 1 angular integral and sum with the appropriate
c...  radial integral
c...
              if(mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0)then
                sum=q1(ijk,0)*xyzi(i,j,k)
                if(sum.ne.zero)then
                do a=0,min(am,i)
                  d=i-a
                  if(d.le.dm)then
                  do b=0,min(bm,j)
                    e=j-b
                    if(e.le.em)then
                    do c=0,min(cm,k)
                      f=k-c
                      if(f.le.fm)then
                        ar1(a,b,c,d,e,f)=ar1(a,b,c,d,e,f)+sum
                      endif
                    enddo
                    endif
                  enddo
                  endif
                enddo
                endif
              endif
            endif
          enddo
        enddo
      enddo
      end
c======================================================================
      subroutine sumar2n(q2,yra,ika,yrb,ikb,xyzi,lp,am,bm,cm,dm,em,fm,
     $                  nmaxa,nmaxb,maxn,maxl,mijk,ar2)
      implicit none
c...
c...  MM august 2004
c...
c...  This subroutine computes sums of type 2 angular and radial
c...  integrals and accumulates them in ar2:
c...
c...
c...                      lp+c+b+a lp+f+e+d
c...                      ====     ====
c...                      \        \
c... ar2(lp,a,b,c,d,e,f) = >        >  Q2(f+e+d+c+b+a,la,lb) *
c...                      /        /
c...                      ====     ====
c...                      la=0     lb=0
c...
c...                 lp
c...                 ====
c...                 \
c...                ( > omega2(a,b,c,la,lp,m) omega2(d,e,f,lb,lp,m))
c...                 /
c...                 ====
c...                 m=-lp
c...
c...
c...  where omega2 are type 2 angular integrals, and q2 is a type
c...  2 radial integral. The summations have the additional
c...  constraints that (a+b+c+lp)-la and (d+e+f+lp)-lb must be even.
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (23-26)
c...
c...  this subroutine computes all ar2 sums for a given lp value
c...
c...  input parameters
c...
c...  q2            array of type 2 radial integrals
c...  ika           flag for polar angles of vector ka:
c...                    0    vector ka is zero (angles sre undefined)
c...                    1    angle phika is undefined
c...                   >1    general case (thetaka and phika defined)
c...  yra           array of real spherical harmonics computed for ka
c...  ikb           flag for polar angles of vector kb (same as ika)
c...  yrb           array of real spherical harmonics computed for kb
c...  xyzi          table of double factorial products for angular
c...                integrals
c...  lp            angular momentum of the projector
c...  am,bm,cm,     maximum values of a,b and c
c...  dm,em,fm,     maximum values of d,e and f
c...  nmaxa,nmaxb   maximum allowed values of na and nb
c...  maxn          maximum allowed n value, for dimensioning
c...  maxl          maximum allowed l value, for dimensioning
c...  mijk          maximum exponent for angular integrals,
c...                for dimensioning
c...
c...  output parameter
c...
c...  ar2    array of sums of type 2angular and radial integrals
c...
      integer lp,am,bm,cm,dm,em,fm,maxn,maxl,mijk,ika,ikb,nmaxa,nmaxb
      real*8 q2(0:maxn,0:maxl,0:maxl),ar2(0:am,0:bm,0:cm,0:dm,0:em,0:fm)
      real*8 yra,yrb,xyzi
c...
      integer a,b,c,d,e,f,abc,def,abcdef,la,lai,laf,lam,lb,lbi,lbf,lbm,m
      integer lan,lbn
      integer bc,ef,iabc,idef
      real*8  s,s1,qq2
c...
c...  scratch memory for angular integrals
c...
      integer lmax,lpmax
      parameter (lmax=16,lpmax=6)
      real*8 o1,o2a(-lpmax:lpmax),
     $                  o2b(-lpmax:lpmax)
c...
c...  data for do loops
c...
      integer ilo(0:10),ihi(0:10)
      save ilo,ihi
      data ilo/1,  2,  5, 11, 21, 36, 57, 85,121,166,221/
      data ihi/1,  4, 10, 20, 35, 56, 84,120,165,220,286/
      integer i(286),j(286),k(286)
      save i,j,k
      data i/0,0,0,1,0,0,0,1,1,2,0,0,0,0,1,1,1,2,2,3,0,0,0,0,0,1,1,1,1,
     $  2,2,2,3,3,4,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5,0,0,0,0,0,
     $  0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5,6,0,0,0,0,0,0,0,0,1,
     $  1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6,7,0,0,0,0,0,
     $  0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,5,5,
     $  5,5,6,6,6,7,7,8,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
     $  2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,8,8,9,0,
     $  0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,
     $  3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,7,7,7,7,8,8,8,9,9,
     $  10/
      data j/0,0,1,0,0,1,2,0,1,0,0,1,2,3,0,1,2,0,1,0,0,1,2,3,4,0,1,2,3,
     $  0,1,2,0,1,0,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,0,1,2,3,4,
     $  5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,0,1,2,3,4,5,6,7,0,
     $  1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,0,1,2,3,4,
     $  5,6,7,8,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,
     $  2,3,0,1,2,0,1,0,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,0,1,2,3,4,
     $  5,6,7,0,1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,0,
     $  1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,0,1,2
     $  ,3,4,5,6,7,0,1,2,3,4,5,6,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1
     $  ,0/
      data k/0,1,0,0,2,1,0,1,0,0,3,2,1,0,2,1,0,1,0,0,4,3,2,1,0,3,2,1,0,
     $  2,1,0,1,0,0,5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,6,5,4,3,2,
     $  1,0,5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,7,6,5,4,3,2,1,0,6,
     $  5,4,3,2,1,0,5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,8,7,6,5,4,
     $  3,2,1,0,7,6,5,4,3,2,1,0,6,5,4,3,2,1,0,5,4,3,2,1,0,4,3,2,1,0,3,2,
     $  1,0,2,1,0,1,0,0,9,8,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,7,6,5,4,3,
     $  2,1,0,6,5,4,3,2,1,0,5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,10
     $  ,9,8,7,6,5,4,3,2,1,0,9,8,7,6,5,4,3,2,1,0,8,7,6,5,4,3,2,1,0,7,6,5
     $  ,4,3,2,1,0,6,5,4,3,2,1,0,5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0
     $  ,0/
      real*8 zero,epsi
      parameter (zero=0.0d0,epsi=1.0d-17)
      real*8 st
c...
      if(nmaxa+nmaxb.gt.maxn)
     $  call nerror(1,'Sumar2','maximum n value exceeded',
     $              nmaxa+nmaxb,maxn)
      if(nmaxa.gt.10.or.nmaxb.gt.10)
     $  call nerror(2,'Sumar2','maximum n value exceeded',
     $              max(nmaxa,nmaxb),10)
      if(max(nmaxa,nmaxb)+lp.gt.lmax.or.max(nmaxa,nmaxb)+lp.gt.maxl)
     $  call nerror(3,'Sumar2','maximum l value exceeded',
     $              max(nmaxa,nmaxb)+lp,lmax)
      if(lp.gt.lpmax)
     $  call nerror(4,'Sumar2','maximum lp value exceeded',lp,lpmax)
c...
      if(ika.eq.0)then
        lam=0
      else
        lam=nmaxa+lp
      endif
      if(ikb.eq.0)then
        lbm=0
      else
        lbm=nmaxb+lp
      endif
      do abcdef=0,nmaxa+nmaxb
        do abc=0,min(abcdef,nmaxa)
           def=abcdef-abc
           if(def.le.nmaxb)then
             lai=max(lp-abc,0)
             lai=lai+mod(lp+abc+lai,2)
             laf=min(abc+lp,lam)
             lan=laf-lai+1
             lbi=max(lp-def,0)
             lbi=lbi+mod(lp+def+lbi,2)
             lbf=min(def+lp,lbm)
             lbn=lbf-lbi+1
             if(lai.le.laf.and.lbi.le.lbf)then
               do la=lai,laf,2
                 do lb=lbi,lbf,2
                   qq2=q2(abcdef,la,lb)
                   if(abs(qq2).gt.epsi)then
                     if(lan.lt.lbn)then
                       do iabc=ilo(abc),ihi(abc)
                         a=i(iabc)
                         b=j(iabc)
                         c=k(iabc)
                         call omega2one(a,b,c,lp,la,lpmax,ika,
     $                               yra,mijk,xyzi,o1,o2a)
                         st=zero
                         do m=-lp,lp
                           st=st+o2a(m)
                         enddo
                         if(st.ne.zero)then
                           do idef=ilo(def),ihi(def)
                             d=i(idef)
                             e=j(idef)
                             f=k(idef)
                             call omega2one(d,e,f,lp,lb,lpmax,ikb,yrb,
     $                                  mijk,xyzi,o1,o2b)
                             s1=zero
                             do m=-lp,lp
                               s1=s1+o2a(m)*o2b(m)
                             enddo
                             ar2(a,b,c,d,e,f)=ar2(a,b,c,d,e,f)+qq2*s1
                           enddo
                         endif
                       enddo
                     else
                       do idef=ilo(def),ihi(def)
                         d=i(idef)
                         e=j(idef)
                         f=k(idef)
                         call omega2one(d,e,f,lp,lb,lpmax,ikb,yrb,
     $                                mijk,xyzi,o1,o2b)
                         st=zero
                         do m=-lp,lp
                           st=st+o2b(m)
                         enddo
                         if(st.ne.zero)then
                           do iabc=ilo(abc),ihi(abc)
                             a=i(iabc)
                             b=j(iabc)
                             c=k(iabc)
                             call omega2one(a,b,c,lp,la,lpmax,ika,
     $                                 yra,mijk,xyzi,o1,o2a)
                             s1=zero
                             do m=-lp,lp
                               s1=s1+o2a(m)*o2b(m)
                             enddo
                             ar2(a,b,c,d,e,f)=ar2(a,b,c,d,e,f)+qq2*s1
                           enddo
                         endif
                       enddo
                     endif
                   endif
                 enddo
               enddo
             endif
           endif
         enddo
      enddo
      end
c=======================================================================
      subroutine powcacb(ax,ay,az,bx,by,bz,maxn,al,bl,cax,cay,caz,cbx,cb
     $  y,cbz)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 6/25/2004 19:20:58
c...
c...  Powers of ca and cb and do loop limits
c...
      real*8 ZERO,one
      parameter (ZERO=0.0D0,one=1.0D0)
c...
      integer maxn,al(3),bl(3),n
      real*8 ax,ay,az,bx,by,bz,cax(0:maxn),cbx(0:maxn),cay(0:maxn),cby(0
     $  :maxn),caz(0:maxn),cbz(0:maxn)
c...
      cax(0)=one
      cay(0)=one
      caz(0)=one
      cbx(0)=one
      cby(0)=one
      cbz(0)=one
      if(ax.ne.ZERO)then
        cax(1)=ax
        cax(2)=ax*ax
        do n=3,maxn,1
          cax(n)=ax*cax(n-1)
        enddo
        al(1)=0
      else
        cax(1)=zero
        cax(2)=zero
        al(1)=maxn
      endif
      if(ay.ne.ZERO)then
        cay(1)=ay
        cay(2)=ay*ay
        do n=3,maxn,1
          cay(n)=ay*cay(n-1)
        enddo
        al(2)=0
      else
        cay(1)=zero
        cay(2)=zero
        al(2)=maxn
      endif
      if(az.ne.ZERO)then
        caz(1)=az
        caz(2)=az*az
        do n=3,maxn,1
          caz(n)=az*caz(n-1)
        enddo
        al(3)=0
      else
        caz(1)=zero
        caz(2)=zero
        al(3)=maxn
      endif
      if(bx.ne.ZERO)then
        cbx(1)=bx
        cbx(2)=bx*bx
        do n=3,maxn,1
          cbx(n)=bx*cbx(n-1)
        enddo
        bl(1)=0
      else
        cbx(1)=zero
        cbx(2)=zero
        bl(1)=maxn
      endif
      if(by.ne.ZERO)then
        cby(1)=by
        cby(2)=by*by
        do n=3,maxn,1
          cby(n)=by*cby(n-1)
        enddo
        bl(2)=0
      else
        cby(1)=zero
        cby(2)=zero
        bl(2)=maxn
      endif
      if(bz.ne.ZERO)then
        cbz(1)=bz
        cbz(2)=bz*bz
        do n=1,maxn,1
          cbz(n)=bz*cbz(n-1)
        enddo
        bl(3)=0
      else
        cbz(1)=zero
        cbz(2)=zero
        bl(3)=maxn
      endif
c...
      end
