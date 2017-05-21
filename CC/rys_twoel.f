c *** data on the basis functions are stored in the array basdat
c *** beginning with the address ibas. 13 consecutive numbers are
c *** stored for each (uncontracted) shell.  these are the numbers..
c *** (1) orbital exponent
c *** (2) contraction coefficient (for the s function in an l shell)
c *** (3) contraction coefficient for p functions in an l shell or
c ***    contraction coeff. for the next general contraction
c *** (4)-(10): further contr. coeff. for general contractions
c *** (11)-(13): x, y and z coordinates of the center (not needed
c *** absolutely but convenient
c
c *** in the array inx 12 integers are stored for each contracted
c *** shell. the contraction goes from inx(1,k)+1 to inx(5,k)
c *** if inx is dimensioned as inx(12,ncs). inx(2,k) contains the
c *** atom corresponding to contraction k. inx(3,k) contains the
c *** shell size.. 1 for s,3 for p, 4 for l, 5 for d, 6 for d6,
c *** 7 for f, 10 for f10, 15 for g15. inx(4) contains the number of
c *** general contraction: note that in a segmented scheme it is 0
c *** inx(5) is the last primitive of the (general) contraction
c *** inx(6) to inx (nsym+6-1) (i.e. inx(6) or inx(6),inx(7) and (8)
c *** contain the symmetry pair of the contraction, i.e. another
c *** (general) contraction. this is not very useful, as it is the
c *** symmetry pairs of the contracted functions which is  needed.
c *** inx(9) is unused
c *** inx(11,k)+1 to inx(10,k) gives the indices of the
c *** contracted functions making up the contraction. note that
c *** the contracted functions for this gen. contraction go
c *** from inx(11,k)+1 to inx(10,k).
c *** inx(12,k) gives the function type for contracted shell k.
c
c
c ***
c      integrals are encoded the following way:
c       a 4-byte integer gives nns,npat,k and l.
c       nns is associated with symmetry. npat is a pattern
c       word. its 4 least significant bits (3 to 0) give
c       the following information: bit 3 is set if there
c       is a change in the values of the indices i or j. in
c       this case, these indices are given in a separate
c       4-byte integer  in the next word, see below. later.
c       bit 2 is set in npat if integral (ij kl) is present 
c       (zero is not stored). bit 1 is set if integral (ik jl)
c       is present. bit 0 is set if integral (il jk) is stored.
c       nns,npat,k and l are packed so that nns is the most
c       significant byte. note that with the present packing
c       scheme, the maximum index is 255. the coding can be easily
c       changed to allow longer indices, as nns and npat require
c       only 7 bits altogether. one could reserve 10 or 11 bits
c       for them.
c        if npat>7 (its 3rd bit is set) there follows another
c       4-byte word. it contains two dummy indices, followed by
c       i and j. the reason for this is the i and j change relatively
c       infrequently.
c                 ***    ***     ***     ***
c       next come the integrals in integer form. the unit for this
c       is stored permanently in par(26). there are either 3 integrals
c        (for npat=7), or 2 (for npat=3,5 or 6), or 1 (npat=1 or 2 or 4).
c                 ***    ***     ***     ***
c       if an integral is too big, it is encoded in multiple copies.
c       the limit is set by about 350 000 000, as the integrals are
c       used in a combination like 4(ij kl) -(ik jl) - (il jk) in the
c       closed-shell routine, and integer multiplication is used there.
c
      subroutine twoelinit
c
c  this routine initializes the Rys polynomial integral routines
c
      use memory
      implicit real*8 (a-h,o-z)
      parameter(itrian=10000)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,ictr
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c      common /index/ maxsh,ifp,inx(100)
      common /spherical/iexist,ispher(2000)
      save /ganz/,/number/,/spherical/
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      five=5.0d0
      ten=10.0d0
      ten6=1.0d6
      tenm9=1.0d-8
      pi=rgetrval('pi')
      acc=1.0d-14
      inuc=igetival('inuc')
      ibas=igetival('ibas')
      na=igetival('na')
      nbf=igetival('nbf')
      nsh=igetival('nsh')
      ncf=igetival('ncf')
      ncs=igetival('ncs')
      nsym=igetival('nsym')
c
      ictr=igetival('ictr')
      maxsh=1000   ! temporary
      ibas=igetival('ibas')
c
c   Precompute an integer array for spherical harmonics. see TRANS2SPHER
      call LocateSpher(ncs,bl(ictr),iexist,ispher)
c   Prime the incomplete gamma function engine
      call fprep
c
      end
c============================================================================
      subroutine LocateSpher(ncs,inx,iexist,ispher)
c  This subroutine determines whether there are any spherical harmonics basis
c  functions. If yes, it puts the angular momentum of the shell into the
c  corresponding element of ispher
c  Arguments
c  INTENT(IN)
c  ncs     = number of contracted shells
c  inx     = a (12,*) array, basis set contraction info
c  INTENT(OUT)
c  iexist  = 0 if there are no spherical harmonic basis functions
c  ispher  = ispher(ics) is 2,3,4,5 or 6 for D, F. G, H, I functions
c
      dimension inx(12,*),ispher(*)
      iexist=0
      ihigh=0   ! ihigh=2 for D, 3 for F, 4 for G, 5 for H
      do ics=1,ncs
        it=inx(12,ics)
        ispher(ics)=0
        if(it.eq.4.or.it.eq.6.or.it.eq.11.or.it.eq.12.or.it.eq.13) then
          iexist=1
          if(it.eq.4) ispher(ics)=2     ! D (5 comp.)
          if(it.eq.6) ispher(ics)=3     ! F (7 comp.)
          if(it.eq.11) ispher(ics)=4    ! G (9 comp.)
          if(it.eq.12) ispher(ics)=5    ! H(11 comp.Note: twoel is limited to G)
          if(ihigh.eq.13) ispher(ics)=6 ! I(13 comp. See note above)
        end if
      end do
      end
c==============================================================================
      subroutine twoelmain(ics,jcs,kcs,lcs,eps,
     1                     nonz,xint,buf2)
c
c  this is the calling routine of the rys polynomial integral program
c  it is just a wrapper, establishing connection with the basis set data
c  commons BASDAT and INX
c
c  arguments
c  INTENT(IN)
c  ics,jcs,kcs,lcs= four contracted shell indices
c  eps            = a threshold for an estimate on the uncontracted integrals
c  INTENT(OUT)
c  nonz=.true. if there are non-zero integrals calculated
c  xint= the integrals. They are ordered as in a Fortran array xint(lss,kss,jss,iss)
c    where lss=the shell size of the contracted shelllcs (3 for P, 6 for D6),
c    kss=shell size of kcs etc.
c  buf2= temporary location, should be at least as big as the largest Cartesian
c        shell block, i.e 50625 for (GG|GG) 
      use memory
      implicit real*8 (a-h,o-z)
        logical nonz
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,ictr
      common /spherical/iexist,ispher(2000)
      call twoel(ics,jcs,kcs,lcs,eps,bl(ibas),bl(ictr),nonz,xint)
c     
      if(iexist.ne.0) call trans2spher(ics,jcs,kcs,lcs,bl(ictr),
     1                                 ispher,xint,buf2)
      end
c==============================================================================
      subroutine twoel(ics,jcs,kcs,lcs,epsin,basdat,inx,nonz,buf)
      implicit real*8 (a-h,o-z)
c  This routine calculates the two-electron integrals over 4 contracted Cartesian
c    shells. The transformation to spherical harmonics is not yet included. For
c    SCF, DFT and perhaps even MP2, where each integral is used only once, 
c    it is probably a better idea to transfor the density to Cartesian basis,
c    build the Fock matrix in Cartesians, and the transform the Fock matrix to
c    spherical harmonics. For MP2, one should transform the orbitals, although
c    it becomes questionable whether this saves time or loses, as the dimensions
c    increase
c  Arguments
c  INTENT(IN)
c  ic,jc,kc,lc= contracted shell indices
c  epsin      = a threshold on the estimate of the uncontracted integrals
c  basdat(13,*)=primitive shell real data (exponent, location etc. see BASIN
c  inx(12,*)= contraction info (integers) See BASIN
c  INTENT(OUT)
c  nonz=.true. if there are non-zero integrals calculated
c  buf(*): integral buffer. Must be at least 50625 for (GG|GG).
c  The integrals are returned unscaled, as real numbers, BUT ARE MULTIPLIED BY A
c  PERMUTATION FACTOR (1-0.5*delta[i,j])*(1-0.5*delta[k,l])*(1-0.5*delta[ij,kl])
      dimension basdat(13,*)
c  the dimensioning of indxx,indyy,indzz goes up to g functions with 15 cartesian
c   components. the original code dimensioned them only 1296 (d functions)
c   lenn likewise goes up to g 
      logical nonz
      logical exsym,first
      dimension indxx(50625),indyy(50625),indzz(50625)
      save indxx,indyy,indzz
      dimension lenn(8), label(3)
c  maximum contraction length=20 with the dimensions below
      dimension sij(400),pij(1200)
      dimension inx(12,*),buf(*)
      dimension itypes(14),itypold(4)
      data itypold/0,0,0,0/
c  The COMMON /two1/ transfers data between the components of the Rys 2-electron
c  program (for supposedly greater efficiency)
c  ityp,jtyp,ktyp and ltyp are the Cartesian types, i.e they are the same for D 
c  and D6, F and F10, G and G15.
c  n1,n2,n3,n4 are important only if integrals are calculated in triplets
c     and are a permutation of ityp,jtyp,ktyp,ltyp
c  nall=  ?
c  nctn= it is usually an integral counter
c  nord= ordering parameter?, probably important only for triplets
c  maxl= is 0 for (SS|SS), otherwise counts the number of P or L type functions
c     in the integral. It is bigger than 4 if there  are any D,F or G functions
c  iwy=?
c  ilen1 etc. is the number of functions in the Cartesian shell (P=3, L=4, D=6)
c  ia= strating primitive SHELL of the first contracted shell
c  ie= lasyt primitive SHELL of the firts contracted shell
c  ja,je; ka,ke;la,le= same for the second, third, fourth contracted shells
c  lx=
c  eps,eps1= thresholds for integral neglect/approximation. The cut-off is 
c     smooth, by a cubic spline (Meyer's idea?). 
c  ca3,ca2,ca1,ca0 are the coefficients of the cubic spline 
      common /two1/ ityp,jtyp,ktyp,ltyp,
     1 n1,n2,n3,n4,nall,ncnt,nord,maxl,iwy,ilen1,
     1 jlen1,klen1,len1,ia,ie,ja,je,ka,ke,la,le,lx,eps,eps1,ca3,
     2 ca2,ca1,ca0
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(16)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      save /two1/
c      common /index/ maxsh,ifp,inx(12,1)
c      common bl(1),buf(1296),sij(100),pij(300),sik(100),pik(300),skj(100
c     1),pkj(300),indxxx(3888),buf1(1296)
c      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
c     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
c     20)
      data lenn/1,3,4,6,10,15,21,28/
      data label/4,2,1/
      data ca30,ca20,ca10,ca00/-.0150891632d0,0.304526749d0,-.5637860084
     1 d0, .2743484226d0/
c Conversion of true types (s,p,l,d,d6,f,f10,g15,h21,i28,g,h,i,j) to Cartesion type
c  (S,P,L,D,D,F,F,G,H,I,G,H,I,J). Note that the order is not logical because spherical
c   G,H and I (??) functions were added later. J functions are not yet implemented
       data itypes /1,2,3,4,4,5,5,6,7,8,6,7,8,9/
cpp      
c  check out basdat
c      write(6,*) 'basdat'
c      do i=1,nsh
c        write(6,*) (basdat(k,i),k=1,13)
c      end do
c      ca1=ca10
      eps=epsin
      exsym=nganz(1).eq.1
      pi23=(two/pi)**3
      two3=two/three
      three4=three/four
      rnsym=one/(one+float(nsym))
c      lprish=lopt(14)
c      if (lprish.eq.1) read (inp,10) iy,jy,ky,ly
c   10 format (4i2)
c      eps=par(1)
c      eps1=eps*ten
c      eps2=par(2)
c
cpp  print the inx array
c      write(6,*) 'inx'
c      do i=1,ncs
c        write(6,*) (inx(k,i),k=1,12)
c      end do
c      scfac=one/eps
      first=.true.
c      eps3=eps/ten
c
c this is the limit for the neglect of an integral
c
      iold=0
      jold=0
c **   these are used in deciding whether the indices i and j
c **   have to be stored
c
      iwx=0
c
c this is used for (ss,sp), (ss,ps), (sp,ss), (ps,ss) type integrals
c it indicates where the p function is situated
c this is zero if the block is (ss,ss), and 1 if there is only 1 p fct
c memory locations last+1to last+1296, last+1297 to last+2592,
c last+2593,last+3888 serve as temporary storages for thre
c integral blocks
c
   30 nci=0
      ifu=inx(11,ics)
c
c block count and blocksize
c
c      ca3=ca30/eps**2
c      ca2=ca20/eps
c      ca0=ca00*eps
c
         ia=inx(1,ics)+1
         ie=inx(5,ics)
         ityp1=inx(12,ics)
         ityp=itypes(ityp1)
c      print *, 'ics,ia,ie,type=',ics,ia,ie,ityp1
c   ityp1 is the actual type of the function: S,P,L,D,D6,F,F10,G,G15,(H,H21 are
c   not yet implemented). ityp is the Cartesian type. It is the same for D and 
c   D6, F and F10, G and G15, H and H21. 
         ilen1=lenn(ityp)
         llen=inx(3,ics)
         maxi=0
c  What is maxi?
         if (ityp.gt.1) maxi=1
         if (ityp.ge.4) maxi=5
         if (ityp.gt.1) iwx=1
         n1=ityp
         jfu=inx(11,jcs)
         je=inx(5,jcs)
         ja=inx(1,jcs)+1
         jlen=inx(3,jcs)
         jtyp1=inx(12,jcs)
         jtyp=itypes(jtyp1)
            jlen1=lenn(jtyp)
            ijlen=ilen1*jlen1
            maxj=maxi
            if (jtyp.gt.1) maxj=maxj+1
            if (jtyp.ge.4) maxj=maxj+4
            if (jtyp.gt.1.and.maxj.eq.1) iwx=2
            n2=jtyp
            kfu=inx(11,kcs)
c
               ka=inx(1,kcs)+1
               ke=inx(5,kcs)
               klen=inx(3,kcs)
               ktyp1=inx(12,kcs)
               ktyp=itypes(ktyp1)
               klen1=lenn(ktyp)
               ijklen=ijlen*klen1
               maxk=maxj
               if (ktyp.gt.1) maxk=maxk+1
c  Why is maxk increased by 4?? 
c  And why only k? and is this 4 sufficient for f and g functions as well?
c  Or should it be increased? Did the original code work for F functions?
               if (ktyp.ge.4) maxk=maxk+4
               if (ktyp.gt.1.and.maxk.eq.1) iwx=3
c
c calculate overlaps and orbital centers for the primitive shells
c ij,ik and kj
c  Changd to ij only - 2005
c
c
               call precal (ia,ie,ja,je,basdat,pi23,nij,sij,pij)
c permutation factors
      if(ics.eq.jcs) then
        do ij=1,nij
          sij(ij)=sij(ij)*half
        end do
      end if
      fact=one
      if(kcs.eq.lcs) fact=fact*half
      if((ics.eq.kcs.and.jcs.eq.lcs).or.(ics.eq.lcs.and.jcs.eq.kcs))
     1     fact=fact*half          
c
c due to the size of the arrays sij and pij, no more than 20 primitives
c can be defined within a contraction
c
               lfu=inx(11,lcs)
                  la=inx(1,lcs)+1
                  le=inx(5,lcs)
                  llen=inx(3,lcs)
                  ltyp1=inx(12,lcs)
                  ltyp=itypes(ltyp1)
                  len1=lenn(ltyp)
                  maxl=maxk
                  if (ltyp.gt.1) maxl=maxl+1
                  if (ltyp.ge.4) maxl=maxl+4
                  if (ltyp.gt.1.and.maxl.eq.1) iwx=4
                  n4=ltyp
c    \\Do not remove these! 
      n1=ityp
      n3=ktyp
                  if((ityp.ne.itypold(1)).or.(jtyp.ne.itypold(2)).or.
     1             (ktyp.ne.itypold(3)).or.(ltyp.ne.itypold(4))) then
                    call indxyz(indxx,indyy,indzz)
                    itypold(1)=ityp
                    itypold(2)=jtyp
                    itypold(3)=ktyp
                    itypold(4)=ltyp
                  end if

c
c      print *, 'before calling twoe'
      call twoe(sij,pij,basdat,indxx,indyy,indzz,fact,nonz,buf)
c
      end
c=======================================================================
      subroutine intprint(ics,jcs,kcs,lcs,buf)
      use memory
      implicit real*8 (a-h,o-z)
      dimension buf(*)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,ictr
      call intxprint(ics,jcs,kcs,lcs,buf,bl(ictr))
c
      end
c=======================================================================
      subroutine intxprint(ics,jcs,kcs,lcs,buf,inx)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*),buf(*)
      ilen=inx(3,ics)
      jlen=inx(3,jcs)
      klen=inx(3,kcs)
      llen=inx(3,lcs)
      intc=0
      write(6,*) 'K from ',inx(11,KCS)+1,' to ',inx(10,kcs)
      write(6,*) 'L from ',inx(11,LCS)+1,' to ',inx(10,lcs)
      do i=1,ilen
        do j=1,jlen
          write(6,*) 'I,J=',inx(11,ics)+i,inx(11,jcs)+j
          write(6,'(5e15.7)') (buf(intc+kl),kl=1,klen*llen)
          intc=intc+klen*llen
        end do
      end do
c
      end
c=======================================================================
