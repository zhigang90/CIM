c JAN 11, 1995 PP: The definition of INV6 is now 1000/r**6
c SEP 13 PP CHANGED THE LIMIT DETERMINANT IN GDIIS TO 1E-10;
c  PUT IN THE RESTRICTION ON THE DIIS COEFFICIENTS WHICH WERE COMM. OUT
c JUN 28 PP CHANGED IT SO THAT A CARD BEGINNING WITH A BLANK SPACE
c    IS NO MORE MISTAKEN FOR A CONTINUATION LINE IN MACHB
      subroutine bread(na,inp,iout,nq,nprim,ibcode,maxcard,
     >                 ibcontr,mq,wri)
      implicit real*8 (a-h,o-z)
c
c   this routine reads the b matrix information and converts it to the
c   integer arrays ibcode(6,nprim) and ibcontr(20,nq). nprim is
c   the number of primitive internal coordinates, ibcode(1,ip) is the
c   type: 1=stretching, 2=invr (1/r), 3=inv6 (1/r**6),4=bending,
c   5=out of plane, 6=torsion, 7=linear1, 8=linear2
c   ibcode(2,ip) is the coefficient of the coordinate in the composite
c   internal coordinate. It is the product of the overall scale factor
c   for the coordinate and the (normalized coefficiient. It
c   is expressed as an integer, in units of 10**(-8) (to avoid
c   defining a structure containing both real and integer elements)
c   ibcode(3..6,ip) are the atoms participating in the coordinate
c   for stretchings only ibcode(3,ip) and ibcode(4,ip) are used etc.
c   ibcontr(nprim+1) is an array describing the contraction:
c   in goes from ibcontr(1,k-1)+1 to ibcontr(1,k), except when k=1.
c   ibcontr(2,k) holds the number of different atoms appearing in
c   the contracted internal coordinate (max. 18), and ibcontr(3..20,k)
c   gives these atoms.

c                                   INPUT DATA
c        EACH ELEMENTARY VALENCE COORDINATE ON A SEPARATE CARD
c        COL. 1 'K' OR' ' (BLANK). IF 'K' A NEW COORDINATE BEGINS, IF
c        BLANK THEN THE COMPOSITE INTERNAL COORDIMNATE BEGUN EARLIER IS
c         CONTINUED. ANY OTHER CHARACTER TERMINATES THE INPUT
c        COLS. 2-9 SCALE FACTOR FOR THE TOTAL COORDINATE (ONLY IF THERE
c        IS K IN COLUMN 1. BLANK OR ZERO IS INTERPRETED AS 1.0
c     0  COLS. 21-24 COORDINATE TYPE STRE,INVR,BEND,OUT ,TORS,LIN1,LIN2
c        COLS. 31-40,41-50,51-60,61-70 PARTICIPATING ATOMS A,B,C,D
c         FORMAT 4F10.X
c         A AND B ARE GIVEN FOR 'STRE' ORDER ARBITRARY
c         A AND B ARE GIVEN FOR -INVR- ORDER ARBITRARY
c        A,B,C FOR BEND - A AND B ARE END ATOMS, C IS THE APEX ATOM A-C-
c        ATOM A OUT OF THE BCD PLANE - C IS THE CENTRAL ATOM - COORDINA
c        TE POSITIVE IF A IS DISPLACED TOWARD THE VECTOR PRODUCT DB*DC
c        TORSION A-B-C-D, POSITIVE AS IN THE WILSON-DECIUS-CROSS BOOK
c        NOTE THAT THE VALUE OF THE COORDINATE VARIES BETWEEN -PI/2 TO
c        3PI/2   N O T  BETWEEN -PI/2 TO +PI/2
c        LIN1 L  COLLINEAR BENDING A-B-C DISTORTED IN THE PLANE OF ABD
c        POSITIVE IF A MOVES TOWARD D
c         LIN2 LINEAR BENDING A-C-B DISTORTED PERPENDICULAR TO THE PLANE
c        ABD - POSITIVE IF A MOVES TOWARD THE VECTOR CROSS PRODUCT CD*CA
c        THE LINEAR BENDINGS ARE A-C-B, I. E. THE APEX ATOM IS THIRD
c
      parameter(ntyp=8)
      logical wri,already
      character*1 wew(80),wort(3),we
      character*4 tipus(ntyp),blank,typ,tlast
      character*10 tencha,lblank,line*128
      dimension ibcode(6,maxcard),ibcontr(20,mq),a(4),ia(4)
      include "common/cmpp"
      equivalence (ka,ia(1)), (kb,ia(2)), (kc,ia(3)), (kd,ia(4))
      data tipus/'STRE','INVR','INV6','BEND','OUT ','TORS',
     1  'LIN1','LIN2'/
      data wort/'K',' ','C'/
      data anull/1.0d0/
      data blank/'    '/
      data lblank/'          '/
      icond=0
      ncard=0
      o=1.0d0
c  this is the (contracted) coordinate counter
      i=0
      c1=0.0d0
      ccx=0.0d0
      tlast=blank
      if (wri) write (iout,10)
   10 format (/,1x,34hdefinition of internal coordinates,/)
   20 if(iprocs.eq.0) then
        read(inp,'(a)',end=21) line
        goto 22
   21   line='---'
      end if
   22 continue
      if(nprocs.gt.1) call global_brdcst_str(line)
      if(line(1:3).eq.'---') goto 45
      read (line,30) wew
      we=wew(1)
   30 format (80a1)
      if(iprocs.eq.0) backspace inp
      do 40 k=1,3
         if (we.eq.wort(k)) go to 60
   40 continue
   45 continue
      if (i.eq.0) then
        write(iout,*) ' NO INTERNAL COORDINATES WERE FOUND ON',inp
        if(icond.ne.0) write(icond,*)
     1    ' NO INTERNAL COORDINATES WERE FOUND ON',inp
c        CALL ERROR
         call fehler
      end if
      c1=sqrt(1.0d0/c1)*ccx
c  put in the last value for ibcontr
      nq=i
      if(nq.gt.mq) call error('mq too small','bread')
      ibcontr(1,nq)=ncard
      if(i.eq.1) then
        ilow=1
      else
        ilow=ibcontr(1,i-1)+1
      end if
      ihigh=ibcontr(1,i)
      do 50 k=ilow,ihigh
         ibcode(2,k)=int(dble(ibcode(2,k))*c1)
   50 continue
      go to 370
   60 continue
c     CHECK IF THE BLANK IS REALLY BLANK
      if (we.eq.wort(2)) then
        if(iprocs.eq.0) read(inp,'(a)') line
        if(nprocs.gt.1) call global_brdcst_str(line)
         read (line,65) tencha
         if(iprocs.eq.0) backspace inp
 65      format(a10)
         if (tencha.ne.lblank) go to 45
      end if
      ncard=ncard+1
      if(iprocs.eq.0) read(inp,'(a)') line
      if(nprocs.gt.1) call global_brdcst_str(line)
      read (line,70) we,cc,c,typ,a
   70 format (a1,f9.5,f10.4,a4,6x,4f10.4)
      if(typ.eq.blank.and.ncard.eq.1) then
        write (iout,*) 'bread: error found'
        write (iout,*) 'first coordinate must be defined'
c        call error
      end if
      if (typ.eq.blank) typ=tlast
      tlast=typ
      if (cc.eq.0.0d0) cc=1.0d0
      if (we.eq.wort(1))then
          ccc=ccx
          ccx=cc
          if(i.gt.mq) call error('mq too small','bread')
          if(i.gt.0) ibcontr(1,i)=ncard-1
      end if
      if (c.eq.0.0d0) c=1.0d0
      if (we.eq.wort(2).or.we.eq.wort(3)) c1=c1+c**2
      if (we.eq.wort(1)) then
      if (i.ne.0) then
        if (wri) write (iout,80)
   80   format (1x)
        c1=sqrt(1.0d0/c1)*ccc
      if(i.eq.1) then
        ilow=1
      else
        ilow=ibcontr(1,i-1)+1
      end if
      ihigh=ibcontr(1,i)
      do 90 k=ilow,ihigh
         ibcode(2,k)=int(dble(ibcode(2,k))*c1)
   90 continue
      end if
      i=i+1
      c1=c**2
      end if
  130 do 140 k=1,4
         ia(k)=int(a(k)+0.01d0)
         ibcode(k+2,ncard)=ia(k)
  140 continue
      do 150 k=1,ntyp
         if (typ.eq.tipus(k)) go to 170
  150 continue
      ierror=7
      write (iout,160) i
      if (icond.ne.0) write (icond,160) i
  160 format (/,1x,38hundefined int.coordinate type at no.  ,i3,/,
     1 10(4h****))
      go to 380
  170 if (wri) write (iout,180) i,typ,ia,c,ccx
  180 format (1x,i3,1h.,a8,4i10,f12.5,f12.6)
      if (ka.lt.1.or.ka.gt.na.or.kb.lt.1.or.kb.gt.na) go to 350
      if (k.gt.3.and.(kc.lt.1.or.kc.gt.na)) go to 350
      if (k.gt.4.and.(kd.lt.1.or.kd.gt.na)) go to 350
      if(ncard.gt.maxcard) call error('maxcard too small','bread')
      ibcode(1,ncard)=k
      ibcode(2,ncard)=int(c*1.0d8)
      go to 20
  350 ierror=6
      write (iout,360) i
      write (iout,*) ' BMAT cannot deal with linear atoms that way.'
      write (iout,*) ' You have to include a dummy atom ( outside the',
     >' linear axis ) into your geometry input!'
      if (icond.ne.0) write (icond,360) i
      call error('Atoms erronously defined','bread')
  360 format (/,1x,41hAtoms erronously defined,coordinate no.  ,i3,/,1x,
     110(4h****))
  370 continue
c     fill in incontr
      do 500 i=1,nq
        if(i.eq.1) then
          ilow=1
        else
          ilow=ibcontr(1,i-1)+1
        end if
        ihigh=ibcontr(1,i)
c zero out ibcontr
        do 460 kk=3,12
          ibcontr(kk,i)=0
 460    continue
c  number of different primitives
        mprim=0
        do 490 k=ilow,ihigh
          do 480 l=3,6
            iatom=ibcode(l,k)
            if(iatom.eq.0) go to 480
            already=.false.
            do 470 kk=3,20
              if(iatom.eq.ibcontr(kk,i)) already=.true.
              if(ibcontr(kk+1,i).eq.0) go to 475
 470        continue
 475        if(.not.already) then
              mprim=mprim+1
c              if(mprim.gt.18) call error
              if(mprim.gt.18) then
                write (iout,*) 'bread: error found'
              end if
              ibcontr(mprim+2,i)=iatom
              ibcontr(2,i)=mprim
            end if
 480      continue
 490    continue
 500  continue
      return
c  380 CALL ERROR
 380  call fehler
      end
c
      subroutine machbnew(na,xa,nq,qonly,shftpt,ibcode,ibcontr,
     1   bmat,qq)
c     parameters: input
c     nek=3*na, na=number of nuclei
c     xa(3,*): nuclear coordinates (in Angstrom)
c     nq=number of internal coordinates
c     qonly: if .true., calculate only the coordinate values
c     shftpt: a constant which determines at which point does
c       a torsional coordinate change by 2*pi
c     ibcode(6,nprim): the encoding of primitive internal coordinates
c     ibcontr(20,nq): contains the contraction pattern, the
c     total number of atoms and the atoms participating in the
c     coordinate (max. 18)
c                 output
c     bmat(54,nq): contains the B matrix elements (for at most 18 atoms=
c         54 Cartesians). Note that this is in a way the transpose of
c         B since the SECOND subscript is the internal coordinate
c     qq(nq): the values of the internal coordinates
c
      implicit real*8 (a-h,o-z)
      logical qonly
c
      dimension xa(3,*),qq(nq),bmat(54,nq),ibcode(6,*),ibcontr(20,nq)
      dimension ia(4), u(3), v(3), w(3), z(3), x(3), uu(3), vv(3),
     1 ww(3), zz(3), uv(12)
      equivalence (ka,ia(1)), (kb,ia(2)), (kc,ia(3)), (kd,ia(4))
      equivalence (uv(1),uu(1)), (uv(4),vv(1)), (uv(7),ww(1)), (uv(10),z
     1z(1))
      parameter(one=1.0d0,pi=3.14159 26535 89793 24d0)
      nab=na
c      pi=dacos(-1.d0)
      call zeroit(qq,nq)
      do 1000 i=1,nq
        if(i.eq.1) then
          iprim1=1
        else
          iprim1=ibcontr(1,i-1)+1
        end if
        iprim2=ibcontr(1,i)
        iatoms=ibcontr(2,i)
        if(.not.qonly) then
          call zeroit(bmat(1,i),3*iatoms)
        end if
        do 900  ipr=iprim1,iprim2
          itype=ibcode(1,ipr)
          c=dble(ibcode(2,ipr))*1.0d-8
          ka=ibcode(3,ipr)
          kb=ibcode(4,ipr)
          kc=ibcode(5,ipr)
          kd=ibcode(6,ipr)
      go to (190,200,205,210,230,260,280,300), itype
c
c..... STRETCH
c
  190 call vektor1(uu,r1,ka,kb,xa)
      vv(1)=-uu(1)
      vv(2)=-uu(2)
      vv(3)=-uu(3)
      qq(i)=qq(i)+r1*c
      go to 320
c
c.....INVERSE
c
  200 call vektor1(uu,r1,ka,kb,xa)
      rm1=one/r1
      rm2=rm1**2
      uu(1)=-rm2*uu(1)
      uu(2)=-rm2*uu(2)
      uu(3)=-rm2*uu(3)
      vv(1)=-uu(1)
      vv(2)=-uu(2)
      vv(3)=-uu(3)
      ia(3)=0
      ia(4)=0
      qq(i)=qq(i)+rm1*c
      go to 320
c  ... inverse sixth power multiplied by 100
  205 continue
      uu(1)=xa(1,ka)-xa(1,kb)
      uu(2)=xa(2,ka)-xa(2,kb)
      uu(3)=xa(3,ka)-xa(3,kb)
      rm2=one/(uu(1)**2+uu(2)**2+uu(3)**2)
      rm6=rm2**3*1000.0d0
      rm8=-rm6*rm2*6.0d0
      uu(1)=uu(1)*rm8
      uu(2)=uu(2)*rm8
      uu(3)=uu(3)*rm8
      vv(1)=-uu(1)
      vv(2)=-uu(2)
      vv(3)=-uu(3)
      ia(3)=0
      ia(4)=0
      qq(i)=qq(i)+rm6*c
      go to 320
c.....BENDING
c
  210 call vektor1(u,r1,ka,kc,xa)
      call vektor1(v,r2,kb,kc,xa)
      co=scalar(u,v)
      si=s2(co)
      do 220 l=1,3
         uu(l)=(co*u(l)-v(l))/(si*r1)
         vv(l)=(co*v(l)-u(l))/(si*r2)
         ww(l)=-uu(l)-vv(l)
  220 continue
      ia(4)=0
      qq(i)=qq(i)+c*darcos(co)
      go to 320
c
c.....OUT OF PLANE
c
  230 call vektor1(u,r1,ka,kd,xa)
      call vektor1(v,r2,kb,kd,xa)
      call vektor1(w,r3,kc,kd,xa)
CPP
c     write(*,*) 'ka,kb,kc,kd',ka,kb,kc,kd
CPP
      call normalv (v,w,z)
      steta=scalar(u,z)
      cteta=s2(steta)
      cfi1=scalar(v,w)
      sfi1=s2(cfi1)
      cfi2=scalar(w,u)
      cfi3=scalar(v,u)
      den=cteta*sfi1**2
      st2=(cfi1*cfi2-cfi3)/(r2*den)
      st3=(cfi1*cfi3-cfi2)/(r3*den)
      do 240 l=1,3
         vv(l)=z(l)*st2
         ww(l)=z(l)*st3
  240 continue
      call normalv (z,u,x)
      call normalv (u,x,z)
      do 250 l=1,3
         uu(l)=z(l)/r1
         zz(l)=-uu(l)-vv(l)-ww(l)
  250 continue
CPP
c     write(*,'(4(3F9.5,2x))')(uu(k),k=1,3),(vv(k),k=1,3),(ww(k),k=1,3),
c    &     (zz(k),k=1,3)
cPP
      cx=-c
      if (steta.lt.0.0d0) cx=c
      qq(i)=qq(i)-cx*darcos(cteta)
      go to 320
c
c..... TORSION
c
  260 call vektor1(u,r1,ka,kb,xa)
      call vektor1(v,r2,kc,kb,xa)
      call vektor1(w,r3,kc,kd,xa)
      call normalv (u,v,z)
      call normalv (w,v,x)
      co=scalar(u,v)
      co2=scalar(v,w)
      si=s2(co)
      si2=s2(co2)
      do 270 l=1,3
         uu(l)=z(l)/(r1*si)
         zz(l)=x(l)/(r3*si2)
         vv(l)=(r1*co/r2-1.000000d0)*uu(l)-r3*co2/r2*zz(l)
         ww(l)=-uu(l)-vv(l)-zz(l)
  270 continue
      co=scalar(z,x)
      u(1)=z(2)*x(3)-z(3)*x(2)
      u(2)=z(3)*x(1)-z(1)*x(3)
      u(3)=z(1)*x(2)-z(2)*x(1)
      si3=sqrt(u(1)**2+u(2)**2+u(3)**2)
      co2=scalar(u,v)
      s=arc1(-co,si3)
      if (co2.lt.0.0d0) s=-s
      s1=s
      if(s.gt.shftpt) s=s-2.0d0*pi
      tau=-s
      qq(i)=qq(i)+c*tau
cpp
c      write(*,'(1x,"machb,atoms,co,s3,tau",4i5,6f12.6)') 
c     &   ka,kb,kc,kd,co,si3,tau,s1,s,shftpt
c
c.... REMEMBER THAT THE RANGE OF THIS COORDINATE IS -PI/2 TO 3*PI/2
c.... IN ORDER TO SHIFT THE DISCONTINUITY OFF THE PLANAR POSITION
c
      go to 320
c
c.....LINEAR COPLANAR BENDING
c
  280 call vektor1(u,r1,ka,kc,xa)
      call vektor1(v,r2,kd,kc,xa)
      call vektor1(x,r2,kb,kc,xa)
      co=scalar(v,u)
      co2=scalar(x,v)
      qq(i)=qq(i)+c*(pi-darcos(co)-darcos(co2))
      call normalv (v,u,w)
      call normalv (u,w,z)
      call normalv (x,v,w)
      call normalv (w,x,u)
c
c..... COORDINATE POSITIV IF ATOM A MOVES TOWARDS ATOM D
c
      do 290 l=1,3
         uu(l)=z(l)/r1
         vv(l)=u(l)/r2
         ww(l)=-uu(l)-vv(l)
  290 continue
      ia(4)=0
      go to 320
c
c.....LINEAR PERPENDICULAR BENDING
c
  300 call vektor1(u,r1,ka,kc,xa)
      call vektor1(v,r2,kd,kc,xa)
      call vektor1(z,r2,kb,kc,xa)
      call normalv (v,u,w)
      call normalv (z,v,x)
      do 310 l=1,3
         uu(l)=w(l)/r1
         vv(l)=x(l)/r2
         ww(l)=-uu(l)-vv(l)
  310 continue
      ia(4)=0
      co=scalar(u,w)
      co2=scalar(z,w)
      qq(i)=qq(i)+c*(pi-darcos(co)-darcos(co2))
  320 if (qonly) go to 900
      do 340 j=1,4
        m=ia(j)
        if (m.le.0) go to 340
        iatoms=ibcontr(2,i)
        do 330 kk=1,iatoms
c   search the atoms participating in the coordinate. Is one of them
c   (ikk) the same as the atom
          ikk=ibcontr(kk+2,i)
          if(ikk.eq.m) then
            imm=kk
            go to 335
          end if
 330    continue
c        print *, 'this should not happen, i,ibcontr,ia(1..4),j',
c     1 i,(ibcontr(kk,i),kk=1,iatoms+2),ia,j
        write(6,*) 'machbnew:'
        write(6,*) 'this should not happen: i,ibcontr,ia(1..4),j'
        write(6,*) i,(ibcontr(kk,i),kk=1,iatoms+2),ia,j
 335    continue
        imm3=imm*3
        j3=j*3
        bmat(imm*3-2,i)=bmat(imm*3-2,i)+uv(j3-2)*c
        bmat(imm*3-1,i)=bmat(imm*3-1,i)+uv(j3-1)*c
        bmat(imm*3,i)=bmat(imm*3,i)+uv(j3)*c
  340 continue
 900  continue
 1000 continue
c
      end
c=======================================================================
      double precision function arc1 (x,y)
      implicit real*8 (a-h,o-z)
      parameter(one=1.0d0,pi=3.14159 26535 89793 24d0)
c      pi=dacos(-1.d0)
      if (abs(x).lt.1.0d-9) go to 10
      s=atan(y/x)
      if (x.lt.0.0) s=s+pi
      arc1=s
      return
   10 cot=x/y
      arc1=0.5d0*pi-cot+cot**3/3.0d0
      return
c
      end
c=======================================================================
      double precision function darcos (x)
      implicit real*8 (a-h,o-z)
      pi=dacos(-1.d0)
      if (x.ge.1.0d0) go to 10
      if (x.le.-1.0d0) go to 20
      x1=sqrt(1.0d0-x**2)
      if (abs(x).lt.1.0d-11) go to 30
      s=atan(x1/x)
      if (x.lt.0.0d0) s=s+pi
      darcos=s
      return
   10 darcos=0.0d0
      return
   20 darcos=pi
      return
   30 darcos=0.5d0*pi
      return
c
      end
c=======================================================================
      subroutine vektor1(u,r,i,j,xa)
      implicit real*8 (a-h,o-z)
      dimension u(3), xa(3,*)
c
c       BILDET DEN NORMIERTEN ENTFERNUNGSVEKTOR VOM KERN J NACH KERN I
c        UND DIE ENTFERNUNG R
c
      u(1)=xa(1,i)-xa(1,j)
      u(2)=xa(2,i)-xa(2,j)
      u(3)=xa(3,i)-xa(3,j)
      r=sqrt(u(1)**2+u(2)**2+u(3)**2)
      rr=1.0d0/r
      u(1)=u(1)*rr
      u(2)=u(2)*rr
      u(3)=u(3)*rr
c
      end
c=======================================================================
      subroutine nom (u)
      implicit real*8 (a-h,o-z)
      dimension u(3)
      x=1.d0/sqrt(scalar(u,u))
      do 10 i=1,3
         u(i)=u(i)*x
   10 continue
      return
c
      end
c=======================================================================
      double precision function s2 (x)
      implicit real*8 (a-h,o-z)
      s2=sqrt(1.0d0-x**2)
      return
c
      end
c=======================================================================
      double precision function scalar (u,v)
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3)
      scalar=0.0d0
      do 10 i=1,3
         scalar=scalar+u(i)*v(i)
   10 continue
      return
c
      end
c=======================================================================
      subroutine normalv (u,v,w)
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3), w(3)
c
c     99999...  W WIRD EIN SENKRECHT AUF DIE EBENE(U,V) STEHENDER EINHEI
c      TOR
c
      w(1)=u(2)*v(3)-u(3)*v(2)
      w(2)=u(3)*v(1)-u(1)*v(3)
      w(3)=u(1)*v(2)-u(2)*v(1)
      call nom (w)
      return
c
      end
c
c=======================================================================
      subroutine btxv(nq,nek,b,ibcontr,v,w)
c
c  this routine multiplies a vector v by B transpose and puts
c  the result in w
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              v(nq): input vector
c              OUTPUT
c              w=B(transpose)*v
      implicit real*8 (a-h,o-z)
      dimension b(54,nq),v(nq),w(nek),ibcontr(20,nq)
      call zeroit(w,nek)
      do i=1,nq
        natom=ibcontr(2,i)
        k3=0
        do k=1,natom
          k3=k3+3
          iatom=ibcontr(k+2,i)
          iat3=iatom*3
          w(iat3-2)=w(iat3-2)+b(k3-2,i)*v(i)
          w(iat3-1)=w(iat3-1)+b(k3-1,i)*v(i)
          w(iat3)=w(iat3)+b(k3,i)*v(i)
        end do
      end do
      end
c
      subroutine bxv(nq,nek,b,ibcontr,w,v)
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              w(nek): input vector
c              OUTPUT
c              v=B*w
      implicit real*8 (a-h,o-z)
      parameter(zero=0.0d0)
      dimension b(54,nq),v(nq),w(nek),ibcontr(20,nq)
      do i=1,nq
        s=zero
        natom=ibcontr(2,i)
        k3=0
        do k=1,natom
          k3=k3+3
          iatom=ibcontr(k+2,i)
          iat3=iatom*3
          s=s+b(k3-2,i)*w(iat3-2)
          s=s+b(k3-1,i)*w(iat3-1)
          s=s+b(k3,i)*w(iat3)
        end do
        v(i)=s
      end do
      end
c
      subroutine bdiag(nq,nek,b,ibcontr,dm1)
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              OUTPUT
c              dm1=diag(BBtranspose)**-1 (inverse)
      implicit real*8 (a-h,o-z)
      parameter(zero=0.0d0,one=1.0d0)
      dimension b(54,nq),ibcontr(20,nq),dm1(nq)
      do i=1,nq
        s=zero
        natom=ibcontr(2,i)
        k3=0
        do k=1,natom
          k3=k3+3
          s=s+b(k3-2,i)**2
          s=s+b(k3-1,i)**2
          s=s+b(k3,i)**2
        end do
        dm1(i)=one/s
      end do
      end
c
      subroutine bdiag1(nq,nek,b,ibcontr,dm1)
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              OUTPUT
c              dm1=diag(BtransposeB)**-1 (inverse)
      implicit real*8 (a-h,o-z)
      parameter(one=1.0d0,eps=1.0d-9)
      dimension b(54,nq),ibcontr(20,nq),dm1(nek)
      call zeroit(dm1,nek)
      do i=1,nq
        natom=ibcontr(2,i)
        k3=0
        do k=1,natom
          k3=k3+3
          iatom=ibcontr(k+2,i)
          iat3=iatom*3
          dm1(iat3-2)=dm1(iat3-2)+b(k3-2,i)**2
          dm1(iat3-1)=dm1(iat3-1)+b(k3-1,i)**2
          dm1(iat3)=dm1(iat3)+b(k3,i)**2
        end do
      end do
      do k=1,nek
        if(abs(dm1(k)).gt.eps) then
          dm1(k)=one/dm1(k)
        else
          dm1(k)=eps
        end if
      end do
      end
c
      subroutine bbt(nq,nek,b,ibcontr,x,y,t)
c
c  This subroutine computes y= B(Bt)x. t is just storage
c  parameters: INPUT
c              nq=number of internal coordinates
c              nek=3*na, number of Cartesians
c              b(54,nq): contains the non-zero elements of B
c              ibcontr(20,nq): coding info for B
c              x(nq): input vector
c              OUTPUT
c              y: B(Bt)x
      implicit real*8 (a-h,o-z)
      dimension b(54,nq),ibcontr(20,nq),x(nq),y(nq),t(nek)
      call btxv(nq,nek,b,ibcontr,x,t)
      call bxv(nq,nek,b,ibcontr,t,y)
      end
