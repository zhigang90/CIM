c =====================================================================
c  routines for determining overlap between two different basis sets.
c  modified from existing overlap routines by PP          June 1997
c  further modified for potential stand-alone use by JB   July 1997
c  used in SCF GUESS module
c =====================================================================
      subroutine inton2(itype,  nna,    oneint, inx,    inx2,
     $                  kk1,    kk2,    basdat, basdat2,XC,
     $                  IAN,    ncs,    ncs2,   ncf,    ncf2,
     $                  bl)
      implicit real*8 (a-h,o-z)
c
c This routine calculates one-elctron matrix elements between different
c basis sets. The nuclear info is considered common, though
c
c PARAMETERS:
c
c     itype (in) is 0 for overlap, 1 for h0, 2,for kinetic energy, 3 for
c     dipole, 4 for r square, 5 for second moments, 6 for nuclear
c     potential (at nna), 7 for field at nna, 8 for field gradient nna
c     nna (in) is the total number of nuclei if itype=1 (H0 matrix),
c     else a specified nucleus
c     oneint (out) is the one-electron integral matrix, stored as a
c     Fortran matrix, i.e. ((h(i,j),i=1,ncf),j=1,ncf2)
c     (ncf and ncf2 are the numbers of contr. basis functions in the
c     two sets, respectively).
c     inx (in) holds the contraction info for the first basis; inx2
c     does the same thing for the second basis set
c     kk1,kk2 give the Cartesian  components x,y,z; x',y',z'
c     basdat and basdat2 (in) hold the basis set
c     datnuc (in): nuclear data
c     ncs and ncs2 are the total number of contracted shells
c     ncf and ncf2 are the number of contracted functions in basis 1
c     and basis 2, resp.
c     bl is an integral work array
c  ------------------------------------------------------------------
      dimension inx(12,*), inx2(12,*)
      dimension oneint(ncf,ncf2),basdat(13,*),XC(3,*),IAN(*)
      dimension bl(127008)
cc      character*8 tp(10)
cc      character*1 dir(4)
cc      data tp/'ovrlap','h0-mtx','kinetic','dipole','rsquare','secmms',
cc     1 'nucpot','efield','fdgrad','h0-str'/
cc      data dir/ ' ','x','y','z'/
c     reserve memory for the array used to store the general contr.
c     integrals. This is quite big because we may have up to i functions
c     (28 components), i.e. 28**2=784, and a maximum of 9 general
c     contractions for each function, i.e.  81*784=63504
      is = 1
      iss = is + 63504
      k1=kk1
      k2=kk2
      if (itype.le.2.or.itype.eq.4.or.itype.eq.6) k1=0
c     overlap (0),h0(1), r**2(4) or Vnuc(6) do not have Cartesian info
cc      write (6,80) tp(itype+1),nna,dir(k1+1),dir(k2+1)
cc   80 format (//,1x,a6,i6,3x,2a1,/)
      key=itype
      if (itype.le.1) key=itype+1
ckwol
      if (itype.eq.9) key=2
ckwol
c     total number of 1-el. integrals
      iza=0
c     cycle through the contracted shells and calculate the integrals
      do 60 ics=1,ncs
c       number of gen. contractions
        ngci=inx(4,ics)
ckwol atoms=centers of contracted shells
        iatom=inx(2,ics)
c
         len1=inx(3,ics)
         jfu=0
         do 50 jcs=1,ncs2
           ngcj=inx2(4,jcs)
           jatom=inx2(2,jcs)
            len2=inx2(3,jcs)
            len=len1*len2*(ngci+1)*(ngcj+1)
            call onel2(ics,jcs,len2,len1,ngcj,ngci,.false.,basdat,
     1                 basdat2,XC,bl(is),inx,inx2,key,nna,k1,k2)
cc            write(6,*) ' back from <one12> bl array is:'
cc            write(6,*) (bl(lll),lll=1,len)
c           H0 matrix and H0-start matrix
            if (itype.eq.1 .or. itype.eq.9) then
              do 20 nn=1,nna
                if(itype.eq.9) then
                  if(nn.ne.iatom .and. nn.ne.jatom) go to 20
                endif
                  call onel2(ics,jcs,len2,len1,ngcj,ngci,.false.,basdat,
     1                       basdat2,XC,bl(iss),inx,inx2,6,nn,k1,k2)
                  za=-float(IAN(nn))
                  call add1(bl(iss),za,bl(is),len)
   20         continue
            end if
c          end H0
c.....................
          iij=0
          ifu=inx(11,ics)-len1
          do 45 igc=0,ngci
            ifu=ifu+len1
            jfu=inx2(11,jcs)-len2
            do 45 jgc=0,ngcj
              jfu=jfu+len2
              iff=ifu
              do 40 i1=1,len1
                iff=iff+1
                jff=jfu
                do 40 j1=1,len2
                  jff=jff+1
                       iij=iij+1
                  oneint(iff,jff)=bl(iij)
cc       write(6,*) ' just set oneint  iff:',iff,' jff:',jff,' iij:',iij
cc       write(6,*) ' len1:',len1,'len2:',len2,' oneint:',oneint(iff,jff)
                    iza=iza+1
 40           continue
 45       continue
c.....................
 50   continue
 60   continue
c     return memory
cc      call retmem90(2)
c
      return
      end
c

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine onel2(ics,jcs,jlen,ilen,ngcj,ngci,same,basdat,basdat2,
     1                 XC,ss,inx,inx2,key,nn,k1,k2)
      implicit real*8 (a-h,o-z)
      logical same
      dimension basdat(13,*),basdat2(13,*),XC(3,*)
      dimension xa(3), xb(3), inx(12,*),inx2(12,*), s(784),xint(784),
     1   p(3)
      dimension ss(jlen,ilen,0:ngcj,0:ngci)
      PARAMETER (zero=0.0d0,two=2.0d0,pi=3.14159265358979323844d0)
c -- this is 1/sqrt(12)
      data sqtw/0.2886 7513 4594 8128 8d0/
c
c this subroutine calculates the one-electron integrals for a
c contracted shell and all general contractions belonging to it
c BETWEEN TWO DIFFERENT BASIS SETS
c ics,jcs are the indices of contracted shells
c jlen and ilen are the shell sizes (3 for p...) for the second
c   and first contractions - they are deliberately interchanged
c ngcj and ngci are the general contraction lengths, starting at 0
c if "same" is true then all data are taken from basdat and inx;
c basdat2 is not used at all in this case. This is important to avoid
c transmitting the same argument twice which is not legal
c basdat and basdat2 contain the shell info, datnuc the nucl. info
c ss holds the integrals, for all general contractions over the shells
c format is ss(jlen,ilen,0:ngci,0:ngcj)
c ics and jcs, in the order s(j,i,ngcj,ngci)
c inx is the contraction info
c key gives the kind of matrix element wanted (together with k1 and k2
c the latter can take the values 1,2 or 3 =x,y and z
c key=1 for overlap, 2 for kinetic, 3 for dipole, 4 for r**2, 5 for
c second moments, 6 for nuclear potential at nucleus nn, 7 for electric
c field at nn, 8 for field gradient at nn
c
c   better than using /common/number:
c     call getrval('one',one)
c     call getrval('thre',three)
c     call getrval('four',four)
c     call getrval('two',two)
c     call getrval('pi',pi)
cc      call getival('inuc',inuc)
      twopi=two/pi
      ityp=inx(12,ics)
      jtyp=inx2(12,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
c
      if(ityp.eq.6) ilen1=10
      if(jtyp.eq.6) jlen1=10
c
      if(ityp.eq.11) ilen1=15
      if(ityp.eq.12) ilen1=21
      if(ityp.eq.13) ilen1=28
      if(jtyp.eq.11) jlen1=15
      if(jtyp.eq.12) jlen1=21
      if(jtyp.eq.13) jlen1=28
c
      len1=ilen1*jlen1
c
c for d functions, 6 spaces are needed prior to the transformation
c
c     zero out the buffer
      call zeroit(ss,len*(ngci+1)*(ngcj+1))
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      ja=inx2(1,jcs)+1
      je=inx2(5,jcs)
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
              b=basdat2(1,j)
              do 26 l=1,3
               xb(l)=basdat2(l+10,j)
  26          continue
            end if
            sqb=sqrt(b)
            apb=a+b
            r=zero
            e=a*b/apb
            do 30 l=1,3
               p(l)=(a*xa(l)+b*xb(l))/apb
               r=r+(xa(l)-xb(l))**2
   30       continue
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
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
            call incons2(ityp1,jtyp1,key,nn,k1,k2,a,b,sqa,sqb,s0,
     1                   xa,xb,s,XC)
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
c         call itran1a(xint,s,jlen1)
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
c	    call htran2a(xint,s,ilen)
        end if
c
c  *** add the primitive integrals to the contracted ones
cc        write(6,*) ' About to start Loop 45   s array is'
cc        write(6,*) (s(ll),ll=1,len)
        do 45 igc=0,ngci
          csa=basdat(igc+2,i)
c         cpa is needed only for the l type
          cpa=basdat(3,i)
          do 45 jgc=0,ngcj
            if (same) then
                csb=basdat(jgc+2,j)
                cpb=basdat(3,j)
              else
                csb=basdat2(jgc+2,j)
                cpb=basdat2(3,j)
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
cc        write(6,*) ' just set ss element     ss = ',ss(j1,i1,jgc,igc)
cc        write(6,*) ' j1:',j1,' i1:',i1,' jgc:',jgc,' igc:',igc
   40         continue
   44         continue
   45     continue
   50   continue
   60 continue
c
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine incons2(ityp,jtyp,key,nna,k1,k2,a,b,sqa,sqb,s0,
     1                   xa,xb,xint,XC)
c    this subroutine constructs a set of one-electron integrals  over
c    two primitive shells
c    parameters
c
c  input:
c  ityp,jtyp: the types of the shells, s=1,p=2,l=3,d=4,f=5,g=6,h=7,i=8
c    note that only 6-component d and 10-component f are used here
c    transformation to d5 and f7 is later
c  key is the type of the integral
      implicit real*8 (a-h,o-z)
      dimension XC(3,*)
      dimension xint(784),xa(3), xb(3), yint(784), len(8), la(3), lb(3),
     1ll(3)
      dimension mulx(88), muly(88), mulz(88), nfu(9)
      PARAMETER (zero=0.0d0,two=2.0d0,pi=3.14159265358979323844d0)
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
c
      ln=len(ityp)*len(jtyp)
      do 10 i=1,ln
         yint(i)=zero
   10 continue
      do 20 i=1,3
         la(i)=0
         lb(i)=0
         ll(i)=0
   20 continue
      go to (30,40,120,170,200,210,220,230), key
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
         ft=float(2*(mulx(j)+muly(j)+mulz(j))+3)*b
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
         ft=float((mulx(j)*(mulx(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
   90 continue
      lb(1)=0
      lb(2)=-2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 100 i=ia,ie
      do 100 j=ja,je
         ij=ij+1
         ft=float((muly(j)*(muly(j)-1))/2)
         xint(ij)=xint(ij)-yint(ij)*ft
  100 continue
      lb(2)=0
      lb(3)=-2
      call intar (ityp,jtyp,a,b,s0,xa,xb,la,lb,ll,yint)
      ij=0
      do 110 i=ia,ie
      do 110 j=ja,je
         ij=ij+1
         ft=float((mulz(j)*(mulz(j)-1))/2)
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
  210 continue
      call nucpot (ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),xb(1),xb(2),xb(3),
     1             XC(1,nna),XC(2,nna),XC(3,nna),la,lb,0,xint)
  220 continue
  230 continue
  240 continue
c
      return
      end
