C APR 1 98 PP ADDED NEGLECT FOR THE NUCLEAR ATTRACTION INTEGRALS
c  the criterion is S0*abs(qch)*max(sqa,sqb)
c ....................................................................
c  routines <xtor> onwards modified slightly for potential
c  stand-alone use in SCF GUESS module.
c  routines <inton2> and <onel2>, used in GUESS module, REMOVED
c  and placed in new file - <inton2.f>       ! JB  July 25 1997
c ....................................................................
C PP 7/13/97 I have changed copydibl. It was not working correctly
c   for general contractions, and the symmetry was not working for them
C PP 6/24/97 changed the call of onel: basdat2 is eliminated
C   added inton2 and onel2 for matrix elements between different
C   basis sets
C   Renamed the file intcal2
C PP 5/19/97 set the value of the second threshold only if was really read in
C   lines changed or added: 20-25
c
c============================================================================
c
      subroutine integop(inp)
c  this routine reads the integral options - very few at present
      implicit real*8 (a-h,o-z)
      parameter (nopt=8)
      character*4 options(nopt)
      character*256 chopval(nopt)
      dimension ioptyp(nopt),iopval(3,nopt),ropval(3,nopt),ifound(nopt)
      common /intlim/ limxmem,limblks,limpair
      common /route/ iroute
      common /onethr/onethre
      data options /'thre','onel','ecp','prin','limi','rout','stab',
     *              'ncac'/
      data ioptyp /12,0,0,1,13,1,1,1/
c
c     the default integral threshold is set in basin
c     the value for instability is set in basin
c
      call readopt(inp,nopt,options,ioptyp,iopval,ropval,chopval,ifound)
      if (ifound(1).eq.1) then
c   the integral threshold is in ph notation
        thre=10.0d0**(-ropval(1,1))
c set the second threshold only if it was read in and is reasonable
        call setrval('ithr',thre)
        if(ropval(2,1).gt.5.0d0) then
          thre1=10.0d0**(-ropval(2,1))
          call setrval('ith1',thre1)
c if the second threshold is larger (sharper) than the first,
c  interchange them
          if(ropval(2,1).gt.ropval(1,1)) then
            call setrval('ithr',thre1)
            call setrval('ith1',thre)
          end if
        end if
      end if
      if (ifound(2).eq.1) then
        call setival('onel',1)
      end if
      if (ifound(3).eq.1) then
c    what do we do now?
      end if
      if(ifound(4).eq.1) then
        iprint=iopval(1,4)
      else
        iprint=0
      end if
      call setival('iprn',iprint)
      if (ifound(5).eq.1) then
        limxmem=ropval(1,5)
      limblks=ropval(2,5)
      limpair=ropval(3,5)
      call getival('iout',iout)
      call getival('icond',icond)
      write(iout,500) limxmem,limblks,limpair
cc      write(icond,500) limxmem,limblks,limpair
 500    format('Integral limits redefined, mem, blks, pairs=',i8,2i5)
      end if
      if(ifound(6).eq.1) then
        iroute=iopval(1,6)
        call setival('iroute',iroute)
      end if
      if (ifound(7).eq.1) then
c        numerical instability in the integral program
         istab=iopval(1,7)
          call setival('stab',istab)  ! meaning : +1 - stable, -1 unstable
          call getival('iout',iout)
          if(istab.lt.0) then
             write(iout,*)
     *       'instability overwritten : consider as unstable'
          else
             write(iout,*)
     *       'instability overwritten : consider as stable'
          endif
      end if
      if (ifound(8).eq.1) then
         ncache=iopval(1,8)
         call setival('ncache',ncache)
      end if
c-----------------------------------------------------------------------
C  set the 1-el. integral threshold
C  !!!this routine is not called in most jobs, setup is needed elsewhere
c  get the main integral threshold
      thre=rgetrval('ithr')
c  the nuclear potential threshold is 100 times less
      onethre=0.01d0*thre
c-----------------------------------------------------------------------
      end
c
c============================================================================
c
      subroutine sympar (ncs,nsh,ncf,nsym,nsy,datnuc,basdat,basdat2,
     1                   inx,ifup,ishp,eps,iout)
c     this subroutine determines the symmetry equivalents of the
c     contracted functions. It first makes a symmetry image of the contracted
c     basis functions, and calculates the overlap between this and the
c     original matrix. An overlap which is equal to the norm of the function
c     in abs. value is interpreted as a symmetry pair
c    PARAMETERS:
c  INPUT:
c      ncs: number of contracted shells
c      nsh: number of uncontracted shells
c      ncf: number of contracted basis functions
c      nsym: number of symmetry operations, only Abelian ones are permitted,
c        nsym<8
c      nsy is a 7 element array holding the symmetry operation types
c        1=X=mirror in the yz plane;
c        2=Y=mirror in the xz plane;
c        3=Z=mirror in the xy plane;
c        4=XY=C2 axis around z;
c        5=XZ=C2 axis around y;
c        6=YZ=C2 axis around x;
c        7=XYZ=inversion center;
c        The logic is that the coordinate which changes sign is named
c      datnuc: nuclear info array, see description of nuclear info
c      basdat: basis set info (see there)
c      eps: tolerance for numerical accuracy, usually 1.0E# times the
c           machine accuracy
c      iout: long output file
c  STORAGE:
c      basdat2 holds temporarily the symmetry image of the basis
c  OUTPUT
c      ifup(7,ncf) holds the symmetry images of the contracted functions
c      ishp(7,ncs) holds the symmetry images of the contracted shels
c

      use memory

      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*),basdat(13,*),basdat2(13,nsh)
      dimension inx(12,*), nsy(7), ifup(7,ncf),ishp(7,ncs),mirror(3,7),
     1s(784), ss(784)
c  784 is 28**2, the largest shell size (Cartesian i functions, i28)
c     common /big/bl(1)
      PARAMETER (zero=0.0d0,four=4.0d0)
      data mirror(1,1),mirror(2,1),mirror(3,1),mirror(1,2),mirror(2,2),m
     1irror(3,2),mirror(1,3),mirror(2,3),mirror(3,3),mirror(1,4),mirror(
     22,4),mirror(3,4),mirror(1,5),mirror(2,5),mirror(3,5),mirror(1,6),m
     3irror(2,6),mirror(3,6),mirror(1,7),mirror(2,7),mirror(3,7)/-1,1,1,
     41,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c   what is the significance of 63504? It is the number of gen. contr.
c (at most 9) squared times 28**2 This is quite big because we may have
c  up to i functions (28 components), i.e. 28**2=784, and a maximum of 9
c  general contractions for each function, i.e.  81*784=63504
c
      If(nsym.eq.0) return
      call getmem(63504,iss)
      call getmem(63504,is1)
      eps1=eps*four
      do 160 ns=1,nsym
         nsm=nsy(ns)
c  make a reflected copy of the basis set info in basdat2
c  basdat(11:13) contains the x,y,z coordinates of the basis function center
         do 30 i=1,nsh
            do 20 l=1,13
               basdat2(l,i)=basdat(l,i)
               if (l.gt.10.and.mirror(l-10,nsm).lt.0) then
                 basdat2(l,i)=-basdat(l,i)
               end if
   20       continue
   30    continue
         do 40 i=1,ncf
c        function pair
            ifup(ns,i)=0
   40    continue
c        inx(11,i)+1 is the first contracted function for this c. shell
         do 140 i=1,ncs
            i2=inx(11,i)
c           number of general contractions(0 if segmented)
            ilen=inx(3,i)
            ngc=inx(4,i)
            call onel(i,i,ilen,ilen,ngc,ngc,.true.,basdat,datnuc,
     1      bl(is1),inx,1,0,0,0)
          do 135 igc=0,ngc
c  due to the storage of the generally contracted one-electron integrals
c  as xx(jlen,ilen,0:ngcj,0:ngci), it is best to extract the diagonal
c  block (or the block we need) from the array bl(iss) in a separate
c  routine
            call copydibl(bl(is1),ss,ilen,ngc,igc,igc)
c  lstep=number of functions in a shell (3 for P, 5 for D,...)
            lstep=inx(3,i)
            lstep2=lstep**2
            lstep1=lstep+1
            ityp=inx(12,i)
            do 130 j=1,ncs
               j2=inx(11,j)
               if (inx(12,j).ne.ityp) go to 120
               jlen=inx(3,j)
               ngc2=inx(4,j)
               if(jlen.ne.ilen.or.ngc2.ne.ngc) go to 120
c  calculate the overlap matrix of the basis functiomns with themselves (ONEL, above)
c  and the mixed overlap between the original basis set and its symmetry-transformed
c  image in basdat2 (ONEL2, below). Do this only if the two functions have the same
c  type, contraction length and general contraction length
               call onel2(j,i,jlen,ilen,ngc2,ngc,.false.,basdat2,basdat,
     1          datnuc,bl(iss),inx,inx,1,0,0,0)
               do 125 jgc=0,ngc2
            call copydibl(bl(iss),s,jlen,ngc2,igc,jgc)
               sdif=zero
c  compare absolute values of the diagonal elements of the mixed overlap matrix with
c  the original values
               do 50 ii=1,lstep2,lstep1
                  sdif=sdif+abs(abs(s(ii))-abs(ss(ii)))
   50          continue
c  if the values agree, the two contracted shells form a symmetry pair
               if (sdif.lt.eps1) ishp(ns,i)=j
               i1=i2
               j1=j2
c  lstep is the shell size
               do 110 ii=1,lstep
                go to (60,70,80,90,92,94,96,98,102,104,106,108,109),ityp
c                      s, p, l, d, d6,f,f10,g15,h21,i28,g , h, i types
   60             imult=1
                  go to 100
   70             imult=mirror(ii,nsm)
                  go to 100
   80             if (ii.eq.1) imult=1
                  if (ii.gt.1) imult=mirror(ii-1,nsm)
                  go to 100
   90             if (ii.le.2) imult=1
c
c  z**2 and (x**2-y**2)
c
                  if (ii.eq.3) imult=mirror(1,nsm)*mirror(2,nsm)
c  xy
                  if (ii.eq.4) imult=mirror(1,nsm)*mirror(3,nsm)
c xz
                  if (ii.eq.5) imult=mirror(2,nsm)*mirror(3,nsm)
       go to 100
c      *** d6 type
92     if(ii.le.3) imult=1
c      *** xx,yy or zz
       if(ii.eq.4) imult=mirror(1,nsm)*mirror(2,nsm)
c      *** xy
       if (ii.eq.5) imult=mirror(1,nsm)*mirror(3,nsm)
c      *** xz
       if (ii.eq.6) imult=mirror(2,nsm)*mirror(3,nsm)
c      *** yz
       go to 100
c      *** f type
c      *** order of f functions xxy, xxz, yyx, yyz, zzx, zzy, xyz
 94    if (ii.eq.1. or.ii.eq.6) imult=mirror(2,nsm)
       if (ii.eq.2.or.ii. eq.4) imult=mirror(3,nsm)
       if (ii.eq.3.or.ii.eq.5)  imult=mirror(1,nsm)
       if (ii.eq.7) imult=mirror(1,nsm)*mirror(2,nsm)*mirror(3,nsm)
       go to 100
96     continue
c      *** f10 type: xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
       if(ii.eq.1.or.ii.eq.4.or.ii.eq.6) imult=mirror(1,nsm)
       if(ii.eq.2.or.ii.eq.7.or.ii.eq.9) imult=mirror(2,nsm)
       if(ii.eq.3.or.ii.eq.8.or.ii.eq.10) imult=mirror(3,nsm)
       if(ii.eq.5) imult=mirror(1,nsm)*mirror(2,nsm)*mirror(3,nsm)
       go to 100
 98    continue
c      *** g15 functions
c      *** order : xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz
c                   1     2    3    4    5    6    7    8   9    10
c
c                  yyyy,yyyz,yyzz,yzzz,zzzz
c                   11   12   13   14   15
       imult=1
       if(ii.eq.2.or.ii.eq.7.or.ii.eq.9) imult=mirror(1,nsm)*
     1   mirror(2,nsm)
       if(ii.eq.3.or.ii.eq.8.or.ii.eq.10) imult=mirror(1,nsm)*
     1   mirror(3,nsm)
       if(ii.eq.5.or.ii.eq.12.or.ii.eq.14) imult=mirror(2,nsm)*
     1   mirror(3,nsm)
       go to 100
 102   continue
c      *** h21 function
c      *** order x5,x4y,x4z,x3y2,x3yz,x3z2,x2y3,x2y2z,x2yz2,x2z3
c                1   2   3    4   5     6    7    8     9   10
c
c                xy4,xy3z,xy2z2,xyz3,xz4, y5,y4z,y3z2,y2z3,yz4
c                 11  12   13    14  15   16  17  18   19  20
c
c                z5
c                21
        if(ii.eq.1.or.ii.eq.4.or.ii.eq.6.or.ii.eq.11.or.ii.eq.13.or.
     1    ii.eq.15) imult=mirror(1,nsm)
        if(ii.eq.2.or.ii.eq.7.or.ii.eq.9.or.ii.eq.16.or.ii.eq.20.or.
     1    ii.eq.18) imult=mirror(2,nsm)
        if(ii.eq.3.or.ii.eq.8.or.ii.eq.10.or.ii.eq.17.or.ii.eq.19.or.
     1    ii.eq.21) imult=mirror(3,nsm)
        if(ii.eq.5.or.ii.eq.12.or.ii.eq.14) imult=mirror(1,nsm)*
     1    mirror(2,nsm)*mirror(3,nsm)
        go to 100
 104    continue
c       *** i28 functions
c       order  x6,x5y,x5z,x4y2,x4yz,x4z2,x3y3,x3y2z,x3yz2,x3z3
c               1  2   3   4     5    6    7     8    9    10
c
c              x2y4,x2y3z,x2y2z2,x2yz3,x2z4,xy5,xy4z,xy3z2,xy2z3,xyz4
c               11   12    13    14     15   16   17   18   19    20
c
c              xz5,y6,y5z,y4z2,y3z3,y2z4,yz5,z6
c              21  22  23  24  25   26   27  28
c
c          if(ii.eq.1.or.ii.eq.4.or.ii.eq.6.or.ii.eq.11.or.ii.eq.13
c    1    .or.ii.eq.15.or.ii.eq.22.or.ii.eq.24.or.ii.eq.26.or.ii.eq.28)
c    2     imult=1
        imult=1
        if(ii.eq.2.or.ii.eq.7.or.ii.eq.9.or.ii.eq.16.or.ii.eq.18
     1   .or.ii.eq.20) imult=mirror(1,nsm)*mirror(2,nsm)
        if(ii.eq.3.or.ii.eq.8.or.ii.eq.10.or.ii.eq.17.or.ii.eq.19
     1   .or.ii.eq.21) imult=mirror(1,nsm)*mirror(3,nsm)
        if(ii.eq.5.or.ii.eq.12.or.ii.eq.14.or.ii.eq.23.or.ii.eq.25
     1   .or.ii.eq.27) imult=mirror(2,nsm)*mirror(3,nsm)
        go to 100
  106 continue
c     *** spherical harmonics g functions (9 components)
c     order x4+y4-6x2y2, x4+z4-6z2z2, y4+z4-6y2z2,
c                1            2            3
c     x3y-3xyz2, xy3-3xyz2, x3z-3xy2z, xz3-3xy2z, y3z-3x2yz, yz3-3x2yz
c         4          5          6          7          8          9
         imult=1
         if(ii.eq.4.or.ii.eq.5) imult=mirror(1,nsm)*mirror(2,nsm)
         if(ii.eq.6.or.ii.eq.7) imult=mirror(1,nsm)*mirror(3,nsm)
         if(ii.eq.8.or.ii.eq.9) imult=mirror(2,nsm)*mirror(3,nsm)
         go to 100
  108    continue
c     *** spherical harmonics h functions (11 components)
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
         imult=1
           if(ii.eq.1.or.ii.eq.2.or.ii.eq.9) imult=mirror(1,nsm)
           if(ii.eq.3.or.ii.eq.4.or.ii.eq.8) imult=mirror(2,nsm)
           if(ii.eq.5.or.ii.eq.6.or.ii.eq.7) imult=mirror(3,nsm)
           if(ii.eq.10.or.ii.eq.11)
     1          imult=mirror(1,nsm)*mirror(2,nsm)*mirror(3,nsm)
           go to 100
  109    continue
c     *** spherical harmonics i functions - not yet done
         imult=1
         go to 100
  100             idi=ii*lstep1-lstep
                  s11=s(idi)*dble(imult)
                  j1=j1+1
                  i1=i1+1
cc                  if (abs(s11-ss(idi)).lt.eps) ifup(ns,i1)=j1
cc                  if (abs(s11+ss(idi)).lt.eps) ifup(ns,i1)=-j1
c -- replaced by looser threshold due to symmetry problems   ! JB Nov. 2005
                  if (abs(s11-ss(idi)).lt.eps1) ifup(ns,i1)=j1
                  if (abs(s11+ss(idi)).lt.eps1) ifup(ns,i1)=-j1
                  if (ishp(ns,i).eq.0) ifup(ns,i1)=0
  110          continue
               j2=j1
  125          continue
  120          continue
  130       continue
            i2=i1
  135       continue
cc            if (ishp(ns,i).eq.0) perfec=.false.
  140    continue
         do 150 i=1,ncs
            if (ishp(ns,i).lt.0) ishp(ns,i)=i
  150    continue
cc         write (iout,180) ns,(i,ishp(ns,i),i=1,ncs)
  180 format (/1x,40hpairs of shells for symmetry operation  ,i6,/,(10(2
     1i5,2x)))
cc         write (iout,190) ns,(i,ifup(ns,i),i=1,ncf)
  190 format (/1x,40hpairs of contractions for symm. operatn ,i6,/,(10(2
     1i4,2x)))
  160 continue
      call retmem(2)
      end
c
c============================================================================
c
      subroutine copydibl(x,s,ilen,ngc,igc,jgc)
      implicit real*8 (a-h,o-z)
      dimension x(ilen,ilen,0:ngc,0:ngc),s(*)
      ij=0
      do j=1,ilen
        do i=1,ilen
          ij=ij+1
          s(ij)=x(i,j,igc,jgc)
        end do
      end do
c      call tfer(x(1,1,igc,igc),s,ilen**2)
      end
c
c============================================================================
c
      subroutine inton (itype,nna,oneint,inx,kk1,kk2,
     1                  basdat,datnuc,ncs)

      use memory

      implicit real*8 (a-h,o-z)
c
c     This routine calculates a whole one-electron matrix.
c     It calls onel to calculate integrals over a contracted shell pair.
c     ONEL call incons for primitive integrals
c
c  ARGUMENTS
c     itype (in) is 0 for overlap, 1 for h0, 2,for kinetic energy, 3 for
c     dipole, 4 for r square, 5 for second moments, 6 for nuclear
c     potential (at nna), 7 for field at nna, 8 for field gradient nna
c     9 for the modified H0 matrix (for start), 10 for 3rd and 4th moments
c     nna (in) is the total number of nuclei if itype=1 (H0 matrix),
c     else a specified nucleus
c     oneint (out) is the one-electron integral matrix, stored as the
c     upper triangular in Fortran, i.e. ((h(i,j),i=1,j),j=1,n)
c     inx (in) holds the contraction info
c     kk1,kk2 give the Cartesian  components x,y,z; x',y',z'
c     (for 3rd and 4th moments, see the comments in incons)
c    name, iprint and iout removed
c     name is the code name of the integral. iprint=1 if integrals are
c     to be printed
c     iout is the output file (for printing)
c     basdat and datnuc (in) hold the basis set and nuclear data
c     ncs is the total number of contracted shells
c
ckwol dimension inx(12,*), tp(9), dir(4)
      dimension inx(12,*), tp(10), dir(4)
      dimension oneint(*),basdat(13,*),datnuc(5,*)
c     common /big/bl(1)
      common /onethr/onethre
      character*8 tp
      character*1 dir
      data tp/'ovrlap','h0-mtx','kinetic','dipole','rsquare','secmms',
     1 'nucpot','efield','fdgrad','h0-str'/
      data dir/ ' ','x','y','z'/
c
c     reserve memory for the array used to store the general contr.
c     integrals. This is quite big because we may have up to i functions
c     (28 components), i.e. 28**2=784, and a maximum of 9 general
c     contractions for each function, i.e.  81*784=63504
      call getmem(63504,is)
      call getmem(63504,iss)
      iprint=igetival('iprn')
      iout=igetival('iout')
      inuc=igetival('inuc')
      k1=kk1
      k2=kk2
      if (itype.le.2.or.itype.eq.4.or.itype.eq.6) k1=0
c     overlap (0),h0(1), r**2(4) or Vnuc(6) do not have Cartesian info
      if (iprint.eq.1) write (iout,80) tp(itype+1),nna,
     1    dir(k1+1),dir(k2+1)
   80 format (//,1x,a6,i6,3x,2a1/)
      key=itype
      if (itype.le.1) key=itype+1
ckwol
      if (itype.eq.9) key=2
ckwol
c     total number of 1-el. integrals
      iza=0
c     ifu=0
c  setup threshold for nuclear potential
c  it is 1/100th of the main int. threshold
c     if(itype.eq.6) onethre=0.01d0*rgetrval('ithr')
      onethre=0.01d0*rgetrval('ithr')
c
c     cycle through the contracted shells and calculate the integrals
      do 60 ics=1,ncs
c       number of gen. contractions
        ngci=inx(4,ics)
ckwol atoms=centers of contracted shells
        iatom=inx(2,ics)
c
         jfu=0
         len1=inx(3,ics)
         do 50 jcs=1,ics
           ngcj=inx(4,jcs)
           jatom=inx(2,jcs)
            len2=inx(3,jcs)
            len=len1*len2*(ngci+1)*(ngcj+1)
            call onel (ics,jcs,len2,len1,ngcj,ngci,.true.,basdat,
     1       datnuc,bl(is),inx,key,nna,k1,k2)
c           H0 matrix and H0-start matrix
ckwol       if (itype.eq.1) then
            if (itype.eq.1 .or. itype.eq.9) then
              do 20 nn=1,nna
                if(itype.eq.9) then
                  if(nn.ne.iatom .and. nn.ne.jatom) go to 20
                endif
                  call onel (ics,jcs,len2,len1,ngcj,ngci,.true.,basdat,
     1            datnuc,bl(iss),inx,6,nn,k1,k2)
                  ii=inuc+5*nn-4
                  za=-datnuc(1,nn)
                  call add1(bl(iss),za,bl(is),len)
   20         continue
            end if
c          end H0
c.....................
          iij=-1
          ifu=inx(11,ics)-len1
          do 45 igc=0,ngci
            ifu=ifu+len1
            jfu=inx(11,jcs)-len2
            do 45 jgc=0,ngcj
              jfu=jfu+len2
              iff=ifu
              do 40 i1=1,len1
                iff=iff+1
                jff=jfu
           ii=iff*(iff-1)/2
           do 40 j1=1,len2
             jff=jff+1
                  iij=iij+1
                  if (jff.le.iff) then
                    ij=ii+jff
                    oneint(ij)=bl(is+iij)
                    iza=iza+1
                  end if
 40           continue
 45       continue
c.....................
 50   continue
 60   continue
      n1=1
      if(iza.gt.0) then
c        call wri (oneint,iza,n1,0,name)
        if (iprint.eq.1) then
          write(iout,*)' Matrix  from inton'
          ncf=igetival('ncf')
          ii=0
          do i=1,ncf
            write(iout,90) (oneint(k),k=ii+1,ii+i)
            ii=ii+i
   90       format (2x,5f14.9)
          end do
        end if
      end if
c     return memory
      call retmem(2)
c
c
      end
c
c============================================================================
c
      subroutine onel (ics,jcs,jlen,ilen,ngcj,ngci,same,basdat,
     1  datnuc,ss,inx,key,nn,k1,k2)
      implicit real*8 (a-h,o-z)
      logical same
      dimension basdat(13,*),datnuc(5,*)
      dimension xa(3), xb(3), inx(12,*), s(784),xint(784), p(3)
      dimension ss(jlen,ilen,0:ngcj,0:ngci)
      common /onethr/onethre
      data sqtw/0.2886 7513 4594 8128 8d0/     ! 1/sqrt(12)
      PARAMETER (zero=0.0d0,two=2.0d0,pi=3.14159265358979323844d0)
c
c this subroutine calculates the one-electron integrals for a
c contracted pair of shells,including all general contractions
c ics,jcs are the indices of contracted shells
c jlen and ilen are the shell sizes (3 for p...) for the second
c   and first contractions - they are deliberately interchanged
c ngcj and ngci are the general contraction lengths, starting at 0
c if "same" is true then all data are taken from bbasdat and basdat
c is not used at all. This is important to avoid transmitting the
c  same argument twice which is not legal
c basdat and inx contain the shell info, datnuc the nucl. info
c ss holds the integrals, for all general contractions over the shells
c format is ss(jlen,ilen,0:ngci,0:ngcj)
c ics and jcs, in the order s(j,i,ngcj,ngci)
c inx is the contraction info
c key gives the kind of matrix element wanted (together with k1 and k2
c the latter can take the values 1,2 or 3 =x,y and z, or 11 for xx,...
c 33 for zz (only for the 3rd and 4th moments)
c key=1 for overlap, 2 for kinetic, 3 for dipole, 4 for r**2, 5 for
c second moments, 6 for potential at nucleus nn, 7 for electric
c field at nn, 8 for field gradient at nn, 10 for 3rd and 4th moments
c Note that 7 and 8 are not supported in incons, and 9 is left out
c (it is the special H0 start in inton)
c
      twopi=two/pi
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
c     ilen=inx(3,ics)
c     jlen=inx(3,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
c  f functions - 10 long before transformation
      if(ityp.eq.6) ilen1=10
      if(jtyp.eq.6) jlen1=10
c  g functions - 15 long before transformation
      if(ityp.eq.11) ilen1=15
      if(jtyp.eq.11) jlen1=15
c  h functions - 21 long before transformation
      if(ityp.eq.12) ilen1=21
      if(jtyp.eq.12) jlen1=21
c  i functions - 15 long before transformation
      if(ityp.eq.13) ilen1=28
      if(jtyp.eq.13) jlen1=28
c
c      len1=ilen1*jlen1
c
c for d functions, 6 spaces are needed prior to the transformation
c
c     zero out the buffer
      call zeroit(ss,len*(ngci+1)*(ngcj+1))
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      ja=inx(1,jcs)+1
      je=inx(5,jcs)
      do 60 i=ia,ie
         a=basdat(1,i)
         sqa=sqrt(a)
         do 20 l=1,3
           xa(l)=basdat(l+10,i)
   20    continue
         do 50 j=ja,je
            if (same) then
              b=basdat(1,j)
              do 24 l=1,3
               xb(l)=basdat(l+10,j)
  24          continue
            else
              b=basdat(1,j)
              do 26 l=1,3
               xb(l)=basdat(l+10,j)
  26          continue
            end if
            sqb=sqrt(b)
            apb=a+b
            r=zero
            e=a*b/apb
            do 30 l=1,3
c              xb(l)=bl(jj+l+4)
               p(l)=(a*xa(l)+b*xb(l))/apb
               r=r+(xa(l)-xb(l))**2
   30       continue
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
c
c     write(6,1000) ics,jcs,ityp,jtyp,i,j,ilen1,jlen1,csa,cpa,csb,cpb
c
      ityp1=ityp
      jtyp1=jtyp
c correspondance between ityp and ityp1:
c ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28) 11(g) 12(h) 13(i)
c ityp1 1    2    3    4    4     5    5      6      7       8        6     7    8
      if(ityp.ge.5) ityp1=ityp-1
      if(jtyp.ge.5) jtyp1=jtyp-1
      if(ityp.ge.7) ityp1=ityp-2
      if(jtyp.ge.7) jtyp1=jtyp-2
      if(ityp.eq.11) ityp1=6
      if(jtyp.eq.11) jtyp1=6
      if(ityp.eq.12) ityp1=7
      if(jtyp.eq.12) jtyp1=7
      if(ityp.eq.13) ityp1=8
      if(jtyp.eq.13) jtyp1=8
c
c     *** no special type is needed for d6 at this point
c
c Provide an estimate for the integrals
        if(key.eq.6) then
          qch=datnuc(1,nn)
          estim=s0*abs(qch)*max(sqa,sqb)
c         estim is the inverse square root of the larger exponent
c         times the nuclear charge times the primitive overlap s0
          if(estim.lt.onethre) go to 50
c  jump out of the loop (CYCLE)
        end if
            call incons(ityp1,jtyp1,key,nn,k1,k2,a,b,s0,xa,xb,
     1                  s,datnuc)
c
c     transform the first subscript if needed
c   d
        if (ityp.eq.4) then
           call dtran1a(xint,s,jlen1)
        end if
c  f
        if (ityp.eq.6) then
          call ftran1a(xint,s,jlen1)
        end if
c transform the 15 Cartesian g functions to 9 spherical harmonics components
        if(ityp.eq.11) then
          call gtran1a(xint,s,jlen1)
        end if
c transform the 21 Cartesian h functions to 11 spherical harmonics components
        if(ityp.eq.12) then
          call htran1a(xint,s,jlen1)
        end if
c transform the 28 Cartesian i functions to 13 spherical harmonics components
        if(ityp.eq.13) then
c          call itran1a(xint,s,jlen1)
        end if
c
c  transform the second subscript
c  d
        if(jtyp.eq.4) then
           call dtran2a(xint,s,ilen)
        end if
c  f
        if (jtyp.eq.6) then
          call ftran2a (xint,s,ilen)
        end if
c g
        if(jtyp.eq.11) then
          call gtran2a(xint,s,ilen)
        end if
c  h
        if(jtyp.eq.12) then
          call htran2a(xint,s,ilen)
        end if
c
        if(jtyp.eq.13) then
c          call itran2a(xint,s,ilen)
        end if
c
c  *** add the primitive integrals to the contracted ones
        do 45 igc=0,ngci
          csa=basdat(igc+2,i)
c         cpa is needed only for the l type
          cpa=basdat(3,i)
          do 45 jgc=0,ngcj
            if (same) then
                csb=basdat(jgc+2,j)
                cpb=basdat(3,j)
              else
                csb=basdat(jgc+2,j)
                cpb=basdat(3,j)
              end if
              ij=0
              do 40 i1=1,ilen
                coefi=csa
                if (ityp.eq.3.and.i1.gt.1) coefi=cpa
                do 40 j1=1,jlen
                  ij=ij+1
                  coefj=csb
                  if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
                  ss(j1,i1,jgc,igc)=ss(j1,i1,jgc,igc)+s(ij)*coefi*coefj
   40         continue
   45     continue
   50   continue
   60 continue
c
      end
c
c============================================================================
c
      subroutine incons (ityp,jtyp,key,nna,k1,k2,a,b,s0,xa,xb,
     1                   xint,datnuc)
c
c    this subroutine constructs a set of one-electron integrals  over
c    two primitive shells
c   ARGUMENTS:
c
c  INTENT(IN):
c  ityp,jtyp: the types of the shells, s=1,p=2,l=3,d=4,f=5,g=6,h=7,i=8
c    note that only 6-component d and 10-component f are used here
c    transformation to d5 and f7 is later
c  key is the type of the integral: 1=overlap, 2=kinetic, 3=a dipole
c   component, 4=r**2,, 5=a second moment component, 6=nuclear potential
c   7=a third or fourth moment component
c  nna is the nucleus for key=6
c  k1 and k2 are Cartesian components. They can range from 1=x to 3=z
c   they are used for the dipole (k1) and 2nd moment integrals (k1 and k2)
c  For 3rd and 4th moments, k1 and k2 can be bigger than 10, say  11 is x**2
c    For instance, the xyz component of the 3rd moment will have k1=12 and k2=3,
c    or k1=1,k2=23
c  a and b are the Gaussian exponents
c  s0 is the overlap integral of the spherical parts of the Gaussians
c  xa(3) and xb(3) are the center coordinates of the Gaussians
c  datnuc(5,*) stores the nuclear info in the usual way
c  INTENT(OUT):
c  xint(1:lenj,1:leni) gives the integrals; leni and lenj are the
c   sizes (lengths) of the shells, e.g. 3 for P, 4 for L, 6 for D6
c

      use memory

      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
      dimension xint(784),xa(3), xb(3), yint(784), len(8), la(3), lb(3),
     1ll(3)
      dimension mulx(88), muly(88), mulz(88), nfu(9), xc(3)
c     common /big/bl(1)
      PARAMETER (zero=0.0d0,two=2.0d0)
      common /tape/ inp,inp2,iout
      data len/1,3,4,6,10,15,21,28/
      data nfu/0,1,4,8,14,24,39,60,88/
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
c this threshold was not used at all fetching it took a lot of time
c  get the main integral threshold
c     thre=rgetrval('ithr')
c  the nuclear potential threshold is 10 times less
c     thre=0.1d0*thre
      ln=len(ityp)*len(jtyp)
      do 10 i=1,ln
         yint(i)=zero
   10 continue
      do 20 i=1,3
         la(i)=0
         lb(i)=0
         ll(i)=0
   20 continue
c            1  2   3   4   5   6   7   8   9   10
      go to (30,40,120,170,200,210,220,230,240,250), key
C  30=overlap 40=kinetic, 120=dipole,170=r**2, 200= 2nd moment,
C  210=nuclear pot
   30 call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      go to 240
   40 call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ia=nfu(ityp)+1
      ie=nfu(ityp+1)
      ja=nfu(jtyp)+1
      je=nfu(jtyp+1)
      ij=0
      do 50 i=ia,ie
      do 50 j=ja,je
         ij=ij+1
         ft=dble(2*(mulx(j)+muly(j)+mulz(j))+3)*b
         xint(ij)=yint(ij)*ft
   50 continue
      lb(1)=2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ft=-two*b**2
      do 60 i=1,ln
   60 xint(i)=xint(i)+yint(i)*ft
      lb(1)=0
      lb(2)=2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      do 70 i=1,ln
   70 xint(i)=xint(i)+yint(i)*ft
      lb(2)=0
      lb(3)=2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      do 80 i=1,ln
   80 xint(i)=xint(i)+yint(i)*ft
      if (jtyp.lt.4) go to 240
      lb(3)=0
      lb(1)=-2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 90 i=ia,ie
      do 90 j=ja,je
         ij=ij+1
         ft=dble((mulx(j)*(mulx(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
   90 continue
      lb(1)=0
      lb(2)=-2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 100 i=ia,ie
      do 100 j=ja,je
         ij=ij+1
         ft=dble((muly(j)*(muly(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
  100 continue
      lb(2)=0
      lb(3)=-2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 110 i=ia,ie
      do 110 j=ja,je
         ij=ij+1
         ft=dble((mulz(j)*(mulz(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
  110 continue
      go to 240
  120 go to (130,140,150), k1
  130 ll(1)=1
      go to 160
  140 ll(2)=1
      go to 160
  150 ll(3)=1
  160 call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      go to 240
  170 ll(1)=2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      ll(1)=0
      ll(2)=2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      do 180 i=1,ln
  180 xint(i)=xint(i)+yint(i)
      ll(2)=0
      ll(3)=2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      do 190 i=1,ln
  190 xint(i)=xint(i)+yint(i)
      go to 240
  200 ll(k1)=ll(k1)+1
      ll(k2)=ll(k2)+1
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      go to 240
  210 continue
c     ind=inuc+5*nna-4
c     qch=bl(ind)
c     xc(1)=bl(ind+1)
c     xc(2)=bl(ind+2)
c     xc(3)=bl(ind+3)
      xc(1)=datnuc(2,nna)
      xc(2)=datnuc(3,nna)
      xc(3)=datnuc(4,nna)
        call nucpot (ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),
     1  xb(1),xb(2),xb(3),xc(1),xc(2),xc(3),la,lb,0,xint)
      go to 240
  220 continue
      go to 240
  230 continue
      go to 240
  250 continue
c     print *,'key,k1,k2',key,k1,k2
      k11=0
      k22=0
      if(k1.gt.10) then
        k10=k1/10
        k11=k1-10*k10
      else
        k10=k1
      end if
      if(k2.gt.10) then
        k20=k2/10
        k22=k2-10*k20
      else
        k20=k2
      end if
CPP
c      print *, 'k10,k11,k20,k22',k10,k11,k20,k22
      if(k10.ge.1.and.k10.le.3) then
        ll(k10)=ll(k10)+1
      else
        call nerror(1,'incons',
     1  'Cartesian directions are wrong, <1 or >3',k10,0)
      end if
      if(k11.ge.1.and.k11.le.3) then
        ll(k11)=ll(k11)+1
      else
        call nerror(1,'incons',
     1  'Cartesian directions are wrong, <1 or >3',k11,0)
      end if
      if(k20.ge.1.and.k20.le.3) then
        ll(k20)=ll(k20)+1
      else
        call nerror(1,'incons',
     1  'Cartesian directions are wrong, <1 or >3',k20,0)
      end if
      if(k22.ge.0.and.k22.le.3) then
        ll(k22)=ll(k22)+1
      else
        call nerror(1,'incons',
     1  'Cartesian directions are wrong, <1 or >3',k22,0)
      end if
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,xint)
      go to 240
  240 continue
c
c     write(iout,1100) key,k1,k2,a,b
c
c
c     write(iout,1200) (xint(i),i=1,ln)
c
      return
c
      end
c
c============================================================================
c
      subroutine xtor (l1,l2,l3,a,b,ax,bx,fctor)
      implicit real*8 (a-h,o-z)
c
c generates the polynomial part of integral (x-ax)**k1*(x-bx)**k2*
c x**l3*exp(-a*(x-ax)**2-b*(x-bx)**2), ((fctor(k1,k2),k1=0,l1), k2=0,l2)
c
c results are produced in fctor. Although not dimensioned this way,
c fctor is filled as fctor(0:l2,0:l1)
c
c The maximum value of l1 and l2 is currently 8. This allows up to
c i type functions, taking into account that the kinetic energy term
c needs higher angular momentum.
c
c this is the hottest part if the 1-el code - to speed up:
c handcoded the exponentiation for small integers (library calls
c are required for variable exponents and are way too expensive)
c exchanged division by the same value for inverse multiplication
c
      dimension root(55), weig(55), fctor(49), ipoin(10), xx(8)
      common /hermit/ root
      common /wermit/ weig
      PARAMETER (zero=0.0d0,one=1.0d0)
      data ipoin/0,1,3,6,10,15,21,28,36,45/
      ideg=(l1+l2+l3)/2+1
      apb=a+b
      sqab=sqrt(apb)
      sqabi=1.d0/sqab
      p=(a*ax+b*bx)/apb
      ia=ipoin(ideg)+1
      ie=ipoin(ideg+1)
      i=0
      do 10 i1=ia,ie
         i=i+1
c         xtmp=root(i1)*sqabi
c         remain=(root(i1)-xtmp*sqab)*sqabi
c         xx(i)=xtmp+remain+p
cc these 3 lines are equiv. to xx(i)=root(i1)/sqrt(apb)+p, but faster???
          xx(i)=root(i1)*sqabi+p
   10 continue
      lab=0
      l11=l1+1
      l22=l2+1
      if (l11.lt.1) l11=1
      if (l22.lt.1) l22=1
      do 30 la=1,l11
      do 30 lb=1,l22
         lab=lab+1
         sum=zero
         i=0
         do 20 i1=ia,ie
            i=i+1
            if (la.eq.1) then
               x1=one
            elseif (la.eq.2) then
               x1=xx(i)-ax
            elseif (la.eq.3) then
               x1=(xx(i)-ax)**2
            else
               x1=(xx(i)-ax)**(la-1)
            endif
            if (lb.eq.1) then
               x2=one
            elseif (lb.eq.2) then
               x2=xx(i)-bx
            elseif (lb.eq.3) then
               x2=(xx(i)-bx)**2
            else
               x2=(xx(i)-bx)**(lb-1)
            endif
            if (l3.eq.0) then
               x3=one
            elseif (l3.eq.1) then
               x3=xx(i)
            elseif (l3.eq.2) then
               x3=xx(i)**2
            else
               x3=xx(i)**l3
            endif
            sum=sum+weig(i1)*x1*x2*x3
c           sum=sum+weig(i1)*x1**(la-1)*x2**(lb-1)*x3**l3
   20    continue
         fctor(lab)=sum
   30 continue
      return
c
      end
c============================================================================
      subroutine intar(ityp,  jtyp,  a,  b,  s0,
     1                  xa,    xb,    la, lb, l3,
     2                  xint)
      implicit real*8 (a-h,o-z)
c
c this routine calculates one-electron integrals over gaussian shells
c it is assumed that the 1s (spherical) part of the gaussian is
c normalized (this is denoted by 'gaussian').
c
c the integrand is (x-xa(1))**la(1)*...(z-xb(3))**lb(3)*x**l3(1)..
c  * z**l3(3)*(x-xa(1))**mulx1*...(z-zb(3))**mulz2*gaussian1*gaussian2
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
c gaussian exponents.
c s0 is the overlap integral between normalized
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
c  Argument list:
c input:
c ityp,jtyp: function types, 1,2,3..8 for s,p,l,d,f,g,h,i
c a,b: Gaussian exponents
c s0: overlap between normalized 1s gaussians
c xa(3),xb(3): orbital centers
c la(3),lb(3),l3(3) : cartesian exponents, see above.
c output:
c  xint(jsh,ish)
c  arguments
      dimension xa(3),xb(3),xint(784),la(3),lb(3),l3(3)
c  data
c      dimension idg(8),nfu(9),mulx(88),muly(88),mulz(88)
      dimension idg(8),nfu(9),mul(3,88)
c  local variables
      dimension xfc(63),yfc(63),zfc(63)
c
      data idg/0,1,1,2,3,4,5,6/
c  nfu stores the beginnings and endings of shells in the arrays
c  Shell type ityp goes from nfu(ityp)+1 to nfu(ityp)
      data nfu/0,1,4,8,14,24,39,60,88/
c  mulx(k,ifu) gives the power of the x,y,z (k=1..3) factors in the shells
      data mul/0,0,0,                                               ! s
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
      l1x=iadg+la(1)     ! the highest possible power of (x-Ax)
      l2x=ibdg+lb(1)     ! the highest possible power of (x-Bx)
      l1y=iadg+la(2)     ! the highest possible power of (y-Ay)
      l2y=ibdg+lb(2)     ! the highest possible power of (y-By)
      l1z=iadg+la(3)     ! the highest possible power of (z-Az)
      l2z=ibdg+lb(3)     ! the highest possible power of (z-Bz)
c
      istart=nfu(ityp)+1
      ilen=nfu(ityp+1)-istart+1
      jstart=nfu(jtyp)+1
      jlen=nfu(jtyp+1)-jstart+1
c
      call intar2(ilen, jlen, mul(1,istart), mul(1,jstart), la,
     2            lb,   l3,   a,             b,             xa,
     3            xb,   l1x,  l1y,           l1z,           l2x,
     4            l2y,  l2z,  s0,            xfc,           yfc,
     5            zfc,  xint)
      end
c========================================================================
      subroutine intar2(ilen,  jlen,   mul1,  mul2,   la,
     2                  lb,    l3,     a,     b,      xa,
     3                  xb,    l1x,    l1y,   l1z,    l2x,
     4                  l2y,   l2z,    s0,    xfc,   yfc,
     5                  zfc,   xint)
      implicit real*8 (a-h,o-z)
c This is the working routine to calculate 1-el. overlap, kinetic etc.
c integrals.
c x(i,j)= integral (x-xa(1))**[la(1)+mulx1(i)]*(x-xa(1))**[la(2)+muly1(i)]
c ...(z-xb(3))**[lb(3)+mulz2(j)] *x**l3(1)*y**l3(2)*z**l3(3) *
c  exp[-a(R-XA)**2]*[exp[-b(R-XB)**2] where R=(x,y,z),
c  XA=(xa(1),xa(2),xa(3);
c  i.e. it can include, besides the x,y,z factors associated with the
c  angular momentum functions, extra (x-xa),...(z-zb) factors, as well
c  as extra x,y,z factors
c  the loop over i,j is (xint(j,i),j=1,jlen),i=1,ilen)
c Arguments
c INTENT(IN)
c ilen,jlen         = the number of components in the shell, e.g. 3 for P
c mul1(3,ilen)      = array containing the power of (x-Ax),(y-Ay),(z-Az)
c  in the function. Here Ax,Ay,Az are the coordinates of the first atom,
c  see xa,ya,za as arguments
c mul2(3,jlen)      = the same for (x-Bx),(y-By),(z-Bz)
c la(3)             = extra powers of (x-ax), (y-ay), (z-az)
c lb(3)             = extra powers of (x-bx), (y-by), (z-bz)
c l3(3)             = extra powers of x,y,z
c a,b               = exponents of the Gaussians
c xa(3),xb(3)       = Cartesian coordinates of the centers, Ax..Bz
c l1x,l1y,l1z       = maximum possible power of (x-ax), (y-ay), (z-az)
c l2x,l2y,l2z       = same for (x-bx), (y-by), (z-bz)
c s0                = overlap of the spherical part of the 2 Gaussians
c INTENT - STORAGE
c xfc,yfc,zfc
c INTENT(OUT)
c xint(jlen,ilen)   = one-electron integrals
      dimension xa(3),xb(3),xint(jlen,ilen),la(3),lb(3),l3(3)
      dimension mul1(3,ilen),mul2(3,jlen)
      dimension xfc(0:l2x,0:l1x),yfc(0:l2y,0:l1y),zfc(0:l2z,0:l1z)
      PARAMETER (Zero=0.0d0)
c  call the routine xtor
c   in order to be correct for i functions and a fourth-degree
c   operator (say, hexadecapole), it needs the hermitian roots
c   and weights up to 9th degree; we have it up to the 10th,
c    so it could be expanded to j functions...
      call xtor (l1x,l2x,l3(1),a,b,xa(1),xb(1),xfc)
      call xtor (l1y,l2y,l3(2),a,b,xa(2),xb(2),yfc)
      call xtor (l1z,l2z,l3(3),a,b,xa(3),xb(3),zfc)
      do 40 i=1,ilen
        ixf=mul1(1,i)+la(1)
        iyf=mul1(2,i)+la(2)
        izf=mul1(3,i)+la(3)
        if(ixf.lt.0.or.iyf.lt.0.or.izf.lt.0) go to 30
        do 20 j=1,jlen
          jxf=mul2(1,j)+lb(1)
          jyf=mul2(2,j)+lb(2)
          jzf=mul2(3,j)+lb(3)
          if (jxf.lt.0.or.jyf.lt.0.or.jzf.lt.0) go to 10
          xint(j,i)=xfc(jxf,ixf)*yfc(jyf,iyf)*zfc(jzf,izf)*s0
          go to 20
   10     xint(j,i)=zero
   20   continue
        go to 40
   30   xint(j,i)=zero
   40 continue
      end
c========================================================================
c======================================================================
      subroutine nucpot (n1,n2,a,b,sab,ax,ay,az,bx,by,bz,cx,cy,cz,
     1 la,lb,kk,xint)
C Arguments:
C INPUT:
c  n1=type of the first function (1=S, 2=P,3=L,..)
c  n2=type of the second function
c  a=exponent of the 1st Gaussian
c  b=expon. of the 2nd Gaussian
c  ax,ay,az: center of the 1st Gaussian
c  bx,by,bz: ditto for the 3nd Gaussian
c  cx,cy,cz:  center coordinates of the nucleus
c  la,lb: ang. momenta of the functions
c  kk=
C OUTPUT:
c  xint= integrals, an array
      implicit real*8 (a-h,o-z)
      common /rys/ xrys,rysr(10),rysw(10),ftt(7),nroots
      dimension ndeg(8), la(3), lb(3), xint(784)
      data ndeg/0,1,1,2,3,4,5,6/
c
c     calculate number of roots and xrys for rys
c
      nroots=(ndeg(n1)+ndeg(n2)+2)/2
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
c      call rt123
c
c     calculate two dimensional integrals
c
      call fcaxyz (n1,n2,sab,p,px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz)
c
c     assemble core attraction integrals
c
      call frmca (n1,n2,xint)
      return
c
      end
c
c============================================================================
c
      subroutine frmca (n1,n2,xint)
      implicit real*8 (a-h,o-z)
      PARAMETER (zero=0.0d0)
c     the array xfc,yfc,zfc should be (7,7,7), the indices
c     being xfc(ityp1,ityp2,iroot). This is simulated here
c     by index inrements - a relic from old ages
      common /xyzfc/ xfc(343),yfc(343),zfc(343)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      dimension nbeg(8), nend(8), nxtyp(84), nytyp(84), nztyp(84)
      dimension itm4(7), xint(784)
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
         sum=zero
         go to (30,20,10,5,4,3,2), nroots
    2    sum=sum+xfc(indx+294)*yfc(indy+294)*zfc(indz+294)
    3    sum=sum+xfc(indx+245)*yfc(indy+245)*zfc(indz+245)
    4    sum=sum+xfc(indx+196)*yfc(indy+196)*zfc(indz+196)
    5    sum=sum+xfc(indx+147)*yfc(indy+147)*zfc(indz+147)
   10    sum=sum+xfc(indx+98)*yfc(indy+98)*zfc(indz+98)
   20    sum=sum+xfc(indx+49)*yfc(indy+49)*zfc(indz+49)
   30    sum=sum+xfc(indx)*yfc(indy)*zfc(indz)
         ncount=ncount+1
         xint(ncount)=sum
   40 continue
      return
c
      end
c
c============================================================================
c
      subroutine fcaxyz (n1,n2,sab,p,px,py,pz,ax,ay,az,bx,by,bz,
     $                   cx,cy,cz)
      implicit real*8 (a-h,o-z)
c
c     calculation of the two-dimensional integrals, ix, iy, and iz
c     for the evaluation of the core attraction integrals
c
      PARAMETER (zero=0.0d0,two=2.0d0,pi=3.14159265358979323844d0)
      common /xyzfc/ xfc(343),yfc(343),zfc(343)
      common /rys/ xrys,rysr(10),rysw(10),tff(7),nroots
      common /hermit/ h(55)
      common /wermit/ w(55)
      dimension itm4(7), ndege(8)
      dimension minh(8), maxh(8)
      dimension x1a1(36),y1a1(36),z1a1(36),x1b1(36),y1b1(36),z1b1(36)
      dimension x1a2(36),y1a2(36),z1a2(36),x1b2(36),y1b2(36),z1b2(36)
      dimension x1a3(36),y1a3(36),z1a3(36),x1b3(36),y1b3(36),z1b3(36)
      dimension x1a4(36),y1a4(36),z1a4(36),x1b4(36),y1b4(36),z1b4(36)
      dimension x1a5(36),y1a5(36),z1a5(36),x1b5(36),y1b5(36),z1b5(36)
      dimension x1a6(36),y1a6(36),z1a6(36),x1b6(36),y1b6(36),z1b6(36)
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
      mxh1=(ie1+ie2)/2
      mxh1=maxh(mxh1)
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
            if (ie1.eq.1) go to 10
            x1ax=x1-ax
            y1ay=y1-ay
            z1az=z1-az
            x1a1(i)=x1ax*w(i)
            y1a1(i)=y1ay*w(i)
            z1a1(i)=z1az*w(i)
            if (ie1.eq.2) go to 10
            x1a2(i)=x1a1(i)*x1ax
            y1a2(i)=y1a1(i)*y1ay
            z1a2(i)=z1a1(i)*z1az
            if (ie1.eq.3) go to 10
            x1a3(i)=x1a2(i)*x1ax
            y1a3(i)=y1a2(i)*y1ay
            z1a3(i)=z1a2(i)*z1az
            if (ie1.eq.4) go to 10
            x1a4(i)=x1a3(i)*x1ax
            y1a4(i)=y1a3(i)*y1ay
            z1a4(i)=z1a3(i)*z1az
            if (ie1.eq.5) go to 10
            x1a5(i)=x1a4(i)*x1ax
            y1a5(i)=y1a4(i)*y1ay
            z1a5(i)=z1a4(i)*z1az
            if (ie1.eq.6) go to 10
            x1a6(i)=x1a5(i)*x1ax
            y1a6(i)=y1a5(i)*y1ay
            z1a6(i)=z1a5(i)*z1az
   10       if (ie2.eq.1) go to 20
            x1bx=x1-bx
            y1by=y1-by
            z1bz=z1-bz
            x1b1(i)=x1bx
            y1b1(i)=y1by
            z1b1(i)=z1bz
            if (ie2.eq.2) go to 20
            x1b2(i)=x1b1(i)*x1bx
            y1b2(i)=y1b1(i)*y1by
            z1b2(i)=z1b1(i)*z1bz
            if (ie2.eq.3) go to 20
            x1b3(i)=x1b2(i)*x1bx
            y1b3(i)=y1b2(i)*y1by
            z1b3(i)=z1b2(i)*z1bz
            if (ie2.eq.4) go to 20
            x1b4(i)=x1b3(i)*x1bx
            y1b4(i)=y1b3(i)*y1by
            z1b4(i)=z1b3(i)*z1bz
            if (ie2.eq.5) go to 20
            x1b5(i)=x1b4(i)*x1bx
            y1b5(i)=y1b4(i)*y1by
            z1b5(i)=z1b4(i)*z1bz
            if (ie2.eq.6) go to 20
            x1b6(i)=x1b5(i)*x1bx
            y1b6(i)=y1b5(i)*y1by
            z1b6(i)=z1b5(i)*z1bz
   20    continue
c
         do 110 i1=1,ie1
            ind1=itm4(i1)
         do 110 i2=1,ie2
            ind=ind1+i2
            mpts=(i1+i2)/2
            minh1=minh(mpts)
            maxh1=maxh(mpts)
            fcx=zero
            fcy=zero
            fcz=zero
            do 100 i=minh1,maxh1
               go to (50,40,30,25,23,22,21), i1
   21          f1x=x1a6(i)
               f1y=y1a6(i)
               f1z=z1a6(i)
               go to 60
   22          f1x=x1a5(i)
               f1y=y1a5(i)
               f1z=z1a5(i)
               go to 60
   23          f1x=x1a4(i)
               f1y=y1a4(i)
               f1z=z1a4(i)
               go to 60
   25          f1x=x1a3(i)
               f1y=y1a3(i)
               f1z=z1a3(i)
               go to 60
   30          f1x=x1a2(i)
               f1y=y1a2(i)
               f1z=z1a2(i)
               go to 60
   40          f1x=x1a1(i)
               f1y=y1a1(i)
               f1z=z1a1(i)
               go to 60
   50          f1x=w(i)
               f1y=w(i)
               f1z=w(i)
   60          go to (90,80,70,65,63,62,61), i2
   61          f1x=f1x*x1b6(i)
               f1y=f1y*y1b6(i)
               f1z=f1z*z1b6(i)
               go to 90
   62          f1x=f1x*x1b5(i)
               f1y=f1y*y1b5(i)
               f1z=f1z*z1b5(i)
               go to 90
   63          f1x=f1x*x1b4(i)
               f1y=f1y*y1b4(i)
               f1z=f1z*z1b4(i)
               go to 90
   65          f1x=f1x*x1b3(i)
               f1y=f1y*y1b3(i)
               f1z=f1z*z1b3(i)
               go to 90
   70          f1x=f1x*x1b2(i)
               f1y=f1y*y1b2(i)
               f1z=f1z*z1b2(i)
               go to 90
   80          f1x=f1x*x1b1(i)
               f1y=f1y*y1b1(i)
               f1z=f1z*z1b1(i)
   90          fcx=fcx+f1x
               fcy=fcy+f1y
               fcz=fcz+f1z
  100       continue
            xfc(ind+indz)=fcx
            yfc(ind+indz)=fcy
            zfc(ind+indz)=fcz*rw
  110    continue
  120 indz=indz+49
      return
c
      end
c
c============================================================================
c
      subroutine dtran1a(xint,s,jlen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the 6*jlen1 raw integrals and the transformed ones
c     *** xint is only a temporary storage, it also contains the results
c     *** the order of the raw d functions xx,yy,zz,xy,xz,yz
c     *** the transformed functions are sqrt(1/12)(2zz-xx-yy),
c     *** (1/2)(xx-yy),xy,xz,yz
      dimension xint(*),s(*)
      PARAMETER (two=2.0d0,half=0.5d0)
      data sqtw/0.2886 7513 4594 8128 8d0/
      ij=0
      ij1=0
      do 90 i=1,6
        do 80 j=1,jlen1
         ij=ij+1
         if (i.eq.1) go to 80
         ij1=ij1+1
         if (i.gt.3) go to 70
         if (i.eq.2) xint(ij1)=(two*s(ij+jlen1)-s(ij)-s(ij-jlen1))*sqtw
         if (i.eq.3) xint(ij1)=(s(ij-jlen1-jlen1)-s(ij-jlen1))*half
         go to 80
   70    xint(ij1)=s(ij)
   80    continue
   90 continue
      call tfer(xint,s,ij1)
      end
c
c============================================================================
c
      subroutine dtran2a(xint,s,ilenx)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilenx*6 raw integrals and the transformed ones
c     *** xint holds temporarily the transformed integrals
c     *** the order of the raw d functions xx,yy,zz,xy,xz,yz
c     *** the transformed functions are (2zz-xx-yy),(xx-yy),xy,xz,yz
      dimension xint(*),s(*)
      PARAMETER (two=2.0d0,half=0.5d0)
      data sqtw/0.2886 7513 4594 8128 8d0/
      ij=0
      ij1=0
      do 140 i=1,ilenx
         do 130 j=1,6
            ij=ij+1
            if (j.eq.1) go to 130
            ij1=ij1+1
            if (j.gt.3) go to 120
            if (j.eq.2) xint(ij1)=(two*s(ij+1)-s(ij)-s(ij-1))*sqtw
            if (j.eq.3) xint(ij1)=(s(ij-2)-s(ij-1))*half
            go to 130
  120       xint(ij1)=s(ij)
  130    continue
  140 continue
      call tfer(xint,s,ij1)
      end
c
c============================================================================
c
      subroutine ftran1a(xint,s,jlen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilen1*jlen1 raw integrals
c     *** the results are in xint
c     *** the order of the raw f functions xxx,xxy,xxz,xyy,xyz,
c     *** xzz,yyy,yyz,yzz,zzz
c     *** the order of the transformed integrals (5xxy-rry), (5xxz-rrz),
c     *** (5yyx-rrx), (5yyz-rrz), (5zzx-rrx), (5zzy-rry), xyz
c     *** these are not orthogonal and are not normalized now
c
      dimension xint(jlen1,7),s(jlen1,10)
      PARAMETER (four=4.0d0)
      data sqr40/ 6.3245 5532 0336 7587d0 /
      do 100 j=1,jlen1
        xint(j,1)=four*s(j,2)-s(j,7)-s(j,9)
        xint(j,2)=four*s(j,3)-s(j,8)-s(j,10)
        xint(j,3)=four*s(j,4)-s(j,1)-s(j,6)
        xint(j,4)=four*s(j,8)-s(j,3)-s(j,10)
        xint(j,5)=four*s(j,6)-s(j,1)-s(j,4)
        xint(j,6)=four*s(j,9)-s(j,2)-s(j,7)
        xint(j,7)=sqr40*s(j,5)
 100  continue
      call tfer(xint,s,7*jlen1)
      end
c
c============================================================================
c
      subroutine ftran2a(xint,s,ilen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilen1*jlen1 raw integrals
c     *** the results are in xint
c     *** the order of the raw f functions xxx,xxy,xxz,xyy,xyz,
c     *** xzz,yyy,yyz,yzz,zzz
c     *** the order of the transformed integrals (5xxy-rry), (5xxz-rrz),
c     *** (5yyx-rrx), (5yyz-rrz), (5zzx-rrx), (5zzy-rry), xyz
c     *** these are not orthogonal and are not normalized now
c
c      implicit real*8 (a-h,o-z)
      dimension xint(7,ilen1),s(10,ilen1)
      PARAMETER (four=4.0d0)
      data sqr40/ 6.3245 5532 0336 7587d0 /
      do 100 i=1,ilen1
        xint(1,i)=four*s(2,i)-s(7,i)-s(9,i)
        xint(2,i)=four*s(3,i)-s(8,i)-s(10,i)
        xint(3,i)=four*s(4,i)-s(1,i)-s(6,i)
        xint(4,i)=four*s(8,i)-s(3,i)-s(10,i)
        xint(5,i)=four*s(6,i)-s(1,i)-s(4,i)
        xint(6,i)=four*s(9,i)-s(2,i)-s(7,i)
        xint(7,i)=sqr40*s(5,i)
 100  continue
      call tfer(xint,s,7*ilen1)
      end
c================================================================================
c
      subroutine gtran1a(xint,s,jlen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilen1*jlen1 raw integrals
c     *** the results are in xint
c     *** the order of the raw g functions is:
c     xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
c       1     2     3     4     5    6      7    8      9    10
c     yyyy,  yyyz, yyzz, yzzz,  zzzz
c      11     12    13    14     15
c     *** the order of the transformed integrals is
c     xxxx+yyyy-6xxyy, xxxx+zzzz-6xxzz, yyyy+zzzz-6yyzz
c            1                        2                        3
c     (xxxy-3xyzz), (xyyy-3xyzz), (xxxz-3xyyz), (xzzz-3xyyz).
c          4             5             6             7
c     (yyyz-3xxyz), (yzzz-3xxyz)
c          8             9
c     *** these are not orthogonal and are not normalized now
c
      dimension xint(jlen1,9),s(jlen1,15)
      PARAMETER (three=3.0d0, six=6.0d0)
      parameter (yn=3.2659863237109d0,xn=1.15470053837925d0)
c  3.265986 is sqrt(32/3);  1.1547005 is sqrt(4/3)
      do 100 j=1,jlen1
        xint(j,1)=xn*(s(j,1)+s(j,11)-six*s(j,4))
        xint(j,2)=xn*(s(j,1)+s(j,15)-six*s(j,6))
        xint(j,3)=xn*(s(j,11)+s(j,15)-six*s(j,13))
        xint(j,4)=yn*(s(j,2)-three*s(j,9))
        xint(j,5)=yn*(s(j,7)-three*s(j,9))
        xint(j,6)=yn*(s(j,3)-three*s(j,8))
        xint(j,7)=yn*(s(j,10)-three*s(j,8))
        xint(j,8)=yn*(s(j,12)-three*s(j,5))
        xint(j,9)=yn*(s(j,14)-three*s(j,5))
 100  continue
      call tfer(xint,s,9*jlen1)
c
      end
c================================================================================
c
      subroutine gtran2a(xint,s,ilen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilen1*jlen1 raw integrals
c     *** the results are in xint
c     *** the order of the raw g functions is:
c     xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
c       1     2     3     4     5    6      7    8      9    10
c     yyyy,  yyyz, yyzz, yzzz,  zzzz
c      11     12    13    14     15
c     *** the order of the transformed integrals is
c     *** the order of the transformed integrals is
c     xxxx+yyyy-6xxyy, xxxx+zzzz-6xxzz, yyyy+zzzz-6yyzz
c            1                        2                        3
c     (xxxy-3xyzz), (xyyy-3xyzz), (xxxz-3xyyz), (xzzz-3xyyz).
c          4             5             6             7
c     (yyyz-3xxyz), (yzzz-3xxyz)
c          8             9
c     *** these are not orthogonal and are not normalized now
c
c      implicit real*8 (a-h,o-z)
      dimension xint(9,ilen1),s(15,ilen1)
      PARAMETER (three=3.0d0, six=6.0d0)
      parameter (yn=3.2659863237109d0,xn=1.15470053837925d0)
c  3.265986 is sqrt(32/3);  1.1547005 is sqrt(4/3)
      do 100 i=1,ilen1
        xint(1,i)=xn*(s(1,i)+s(11,i)-six*s(4,i))
        xint(2,i)=xn*(s(1,i)+s(15,i)-six*s(6,i))
        xint(3,i)=xn*(s(11,i)+s(15,i)-six*s(13,i))
        xint(4,i)=yn*(s(2,i)-three*s(9,i))
        xint(5,i)=yn*(s(7,i)-three*s(9,i))
        xint(6,i)=yn*(s(3,i)-three*s(8,i))
        xint(7,i)=yn*(s(10,i)-three*s(8,i))
        xint(8,i)=yn*(s(12,i)-three*s(5,i))
        xint(9,i)=yn*(s(14,i)-three*s(5,i))
 100  continue
      call tfer(xint,s,9*ilen1)
      end
c================================================================================
c
      subroutine htran1a(xint,s,jlen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilen1*jlen1 raw integrals
c     *** the results are in xint
c     *** the order of the raw h functions is:
c     x5,    x4y,   x4z,   x3y2,  x3yz,  x3z2,  x2y3,  x2y2z, x2yz2, x2z3,
c     1       2      3      4      5      6      7       8      9     10
c     xy4,   xy3z,  xy3z2, xyz3,  xz4.   y5,    y4z,   y3z2.  y2z3.  yz4,    z5
c     11      12     13     14     15    16     17      18     19    20      21
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
c     *** these are not orthogonal but are normalized as computed below
c
      dimension xint(jlen1,11),s(jlen1,21)
      PARAMETER (two=2.0d0,three=3.0d0, five=5.0d0,six=6.0d0,ten=10.0d0)
      parameter (xn=0.730296743340221d0,yn=2.3094010767585d0,
     1           zn=9.23760430703401d0,sn=16.0d0/3.0d0)
c  xn=4/sqrt(30)  yn=4/sqrt(3)  zn=16/sqrt(3)   sn=16/3
      do 100 j=1,jlen1
        xint(j,1)=xn*(s(j,1)-ten*s(j,4)+five*s(j,11))
        xint(j,2)=xn*(s(j,1)-ten*s(j,6)+five*s(j,15))
        xint(j,3)=xn*(s(j,16)-ten*s(j,7)+five*s(j,2))
        xint(j,4)=xn*(s(j,16)-ten*s(j,18)+five*s(j,20))
        xint(j,5)=xn*(s(j,21)-ten*s(j,10)+five*s(j,3))
        xint(j,6)=xn*(s(j,21)-ten*s(j,19)+five*s(j,17))
        xint(j,7)=yn*(s(j,3)+s(j,17)-six*s(j,8))
        xint(j,8)=yn*(s(j,2)+s(j,20)-six*s(j,9))
        xint(j,9)=yn*(s(j,11)+s(j,15)-six*s(j,13))
        xint(j,10)=zn*(s(j,5)-s(j,12))
        xint(j,11)=sn*(s(j,5)-two*s(j,14)+s(j,12))
 100  continue
      call tfer(xint,s,11*jlen1)
c
      end
c================================================================================
c
      subroutine htran2a(xint,s,ilen1)
      implicit real*8 (a-h,o-z)
c     *** s contains the ilen1*jlen1 raw integrals
c     *** the results are in xint
c     *** the order of the raw h functions is:
c     x5,    x4y,   x4z,   x3y2,  x3yz,  x3z2,  x2y3,  x2y2z, x2yz2, x2z3,
c     1       2      3      4      5      6      7       8      9     10
c     xy4,   xy3z,  xy3z2, xyz3,  xz4.   y5,    y4z,   y3z2.  y2z3.  yz4,    z5
c     11      12     13     14     15    16     17      18     19    20      21
c     *** the order of the transformed integrals is
c     x5-10x3y2+5xy4.   x5-10x3z2+5xz4,    y5-10x2y3+5x4y,   y5-10y3z2+5yz4,
c            1                 2                 3                  4
c     z5-10x2z3+5x4z,   z5-10y2z3+5y4z,    x4z+y4z-6x2y2z,   x4y+yz4-6x2yz2
c          5                  6                  7                  8
c     xy4+xz4-6xy2z2,   x3yz-xy3z,        x3yz-2xyz3+xy3z
c           9              10                   11
c     *** these are not orthogonal but are normalized as computed below
c
c      implicit real*8 (a-h,o-z)
      dimension xint(11,ilen1),s(21,ilen1)
      PARAMETER (two=2.0d0,three=3.0d0, five=5.0d0,six=6.0d0,ten=10.0d0)
      parameter (xn=0.730296743340221d0,yn=2.3094010767585d0,
     1           zn=9.23760430703401d0,sn=16.0d0/3.0d0)
c  xn=4/sqrt(30)  yn=4/sqrt(3)  zn=16/sqrt(3)   sn=16/3
      do 100 i=1,ilen1
        xint(1,i)=xn*(s(1,i)-ten*s(4,i)+five*s(11,i))
        xint(2,i)=xn*(s(1,i)-ten*s(6,i)+five*s(15,i))
        xint(3,i)=xn*(s(16,i)-ten*s(7,i)+five*s(2,i))
        xint(4,i)=xn*(s(16,i)-ten*s(18,i)+five*s(20,i))
        xint(5,i)=xn*(s(21,i)-ten*s(10,i)+five*s(3,i))
        xint(6,i)=xn*(s(21,i)-ten*s(19,i)+five*s(17,i))
        xint(7,i)=yn*(s(3,i)+s(17,i)-six*s(8,i))
        xint(8,i)=yn*(s(2,i)+s(20,i)-six*s(9,i))
        xint(9,i)=yn*(s(11,i)+s(15,i)-six*s(13,i))
        xint(10,i)=zn*(s(5,i)-s(12,i))
        xint(11,i)=sn*(s(5,i)-two*s(14,i)+s(12,i))
 100  continue
      call tfer(xint,s,11*ilen1)
      end
c================================================================================
      block data herm
      implicit real*8 (a-h,o-z)
      common /hermit/r1,r2(2),r3(3),r4(4),r5(5),r6(6),r7(7),
     1  r8(8),r9(9),r10(10)
      common /wermit/w1,w2(2),w3(3),w4(4),w5(5),w6(6),w7(7),
     1  w8(8),w9(9),w10(10)
      data r1,r2,r3,r4,r5,r6,r7,r8,r9,r10
     1 /0.0d0,
     2 -0.7071067811865475d0 ,0.7071067811865475d0,
     3 -1.224744871391589d0 ,0.0d0 ,1.224744871391589d0,
     4 -1.6506801238857846d0 ,-0.5246476232752903d0 ,
     *  0.5246476232752903d0 , 1.6506801238857846d0,
     5 -2.020182870456086d0, -0.958572464613817d0 ,
     *  0.0d0, 0.958572464613817d0 , 2.020182870456086d0,
     6 -2.350604973674492d0 , -1.335849074013697d0 ,
     * -0.4360774119276165d0 , 0.4360774119276165d0,
     *  1.335849074013697d0 , 2.350604973674492d0,
     7 -2.651961356835233d0 , -1.673551628767471d0 ,
     * -0.816287882858965d0,  0.0d0,  0.816287882858965d0 ,
     *  1.673551628767471d0,  2.651961356835233d0,
     8 -2.930637420257244d0,   -1.981656756695843d0 ,
     * -1.15719371244678d0 , -0.3811869902073222d0 ,
     *  0.3811869902073222d0 , 1.15719371244678d0 ,
     *  1.981656756695843d0 , 2.930637420257244d0,
     9 -3.190993201781528d0 , -2.266580584531843d0,
     * -1.468553289216668d0 , -0.7235510187528376d0 ,
     *  0.0d0 , 0.7235510187528376d0 , 1.468553289216668d0,
     *  2.266580584531843d0 , 3.190993201781528d0,
     X -3.436159118837736d0 , -2.532731674232792d0,
     * -1.756683649299882d0 , -1.036610829789513d0 ,
     * -0.3429013272237051d0,  0.3429013272237051d0 ,
     *  1.036610829789513d0,  1.756683649299882d0 ,
     *  2.532731674232792d0,  3.436159118837736d0/
c
      data w1,w2,w3,w4,w5,w6,w7,w8,w9,w10
     1 /1.77245385090552d0,
     2 0.886226925452758d0, 0.886226925452758d0,
     3 0.295408975150919d0, 1.18163590060368d0, 0.295408975150919d0,
     4 0.0813128354472456d0, 0.804914090005512d0, 0.804914090005512d0,
     * 0.0813128354472456d0,
     5 0.0199532420590459d0, 0.393619323152242d0, 0.945308720482942d0,
     * 0.393619323152242d0, 0.0199532420590459d0,
     6 0.00453000990550885d0, 0.157067320322857d0, 0.724629595224393d0,
     * 0.724629595224393d0, 0.157067320322857d0, 0.00453000990550885d0,
     7 0.000971781245099521d0, 0.0545155828191271d0,0.425607252610128d0,
     * 0.810264617556807d0, 0.425607252610128d0, 0.0545155828191271d0,
     * 0.000971781245099521d0,
     8 0.000199604072211367d0, 0.0170779830074135d0,0.207802325814892d0,
     * 0.661147012558241d0, 0.661147012558241d0, 0.207802325814892d0,
     * 0.0170779830074135d0, 0.000199604072211367d0,
     9 0.0000396069772632635d0, 0.00494362427553720d0,
     * 0.088474527394376d0,  0.432651559002556d0, 0.720235215606051d0,
     * 0.432651559002556d0, 0.088474527394376d0, 0.00494362427553720d0,
     * 0.0000396069772632635d0,
     X 7.64043285523269d-6 ,0.00134364574678122d0, 0.0338743944554813d0,
     * 0.240138611082315d0, 0.610862633735326d0, 0.610862633735326d0,
     * 0.240138611082315d0, 0.0338743944554813d0, 0.00134364574678122d0,
     * 7.64043285523269d-6 /
c
      end
