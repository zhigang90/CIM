      subroutine newintc(nat,inpf3,ian,xc,bndorder, nq,nprim)
      implicit real*8 (a-h,o-z)
c  This routine calls the internal coordinate generators
c  and also generates approximate diagonal force constants
c  Both are written to the file inpf3
c
c  integer, intent(in) ::  nat                    ! number of atomic centers
c  integer, intent(in) ::  inpf3                  ! file number 
c  integer, intent(in) ::  ian(nat)               ! atomic numbers
c  real*8,  intent(in) ::  xc(3,nat)              ! atomic coordinates
c  real*8,  intent(in) ::  bndorder(nat,nat)      ! bond order  matrix
c  integer, intent(out)::  nq                     ! number of composite internal coordinates
c  integer, intent(out)::  nprim                  ! number of primitive internal coordinates
      dimension ian(nat),xc(3,nat),bndorder(nat,nat)
c  calls the coordinate generators
c
c  Number of internal coordinates (composite)
      nq=0
      nprim=0
      call stretches(nat,inpf3,ian,xc,bndorder,nq,nprim) 
cpp
c      do i=1,nat
c        write(*,998) (bndorder(i,j),j=1,nat)
c      end do
c  998 format (10f8.3)
      call bendings(nat,inpf3,ian,xc,bndorder,nq,nprim) 
c  There is a "torsions" in INTC
      call torsion(nat,inpf3,ian,xc,bndorder,nq,nprim) 
      call outofplane(nat,inpf3,ian,xc,bndorder,nq,nprim) 
cpp
c      write(*,*) 'Number of internal coordinates=',nq 
c      write(*,*) 'Number of primitive coordinates=',nprim
c
cpp
      end
c=======================================================================
      subroutine bondlength(i,j,xc,r,rr)
      implicit real*8 (a-h,o-z)
c  Calculates the distance and the squared distance of two centers
c  INPUT: i,j= the two centers, xc(3,*) the coordinates in bohr
c  OUTPUT:  r=distance, rr=r**2
      dimension xc(3,*)
      rr=(xc(1,i)-xc(1,j))**2+(xc(2,i)-xc(2,j))**2+(xc(3,i)-xc(3,j))**2
      r=sqrt(rr)
c
      end
c=======================================================================
      subroutine stretches(nat,inpf3,ian,xc,bndorder,nq,nprim) 
      implicit real*8 (a-h,o-z)
c  Arguments
c  INTENT(IN)
c  nat               = number of centers
c  inpf3             = output input file number (formatted)
c  ian(nat)          = array of nuclear charges
c  xc(3,nat)         = Cartesian coordinates (bohr)
c  bndorder(nat,nat) = bond order matrix
c  INTENT(OUT)
c  nq                = current number of internal coordinates
c  nprim             = current number of primitive coordinates
c  In addition, the routine writes the definition of the internal
c    coordinates and approximate diagonal force constants on unit inpf3
c
      dimension ian(nat),xc(3,nat),bndorder(nat,nat)
      parameter (aa=4.0d0,bhh=-0.129d0,boh=0.186d0,boo=0.574d0,
     &   bsh=0.349d0,bso=0.804d0,bss=1.094d0,zero=0.0d0,bohr=0.529177d0)
c  This routine assigns  the stretching and inverse distance coordinates
c  and calculates approximate force constants
      do i=2,nat
        do j=1,i
          if(bndorder(i,j).eq.zero) cycle
          in=ian(i)
          jn=ian(j)
          if(in.lt.jn) then
            in=ian(j)
            jn=ian(i)
          end if
c
          if(in.le.2) then
            iper=1
          else if(in.le.10) then
            iper=2
          else
            iper=3
          end if
          if(jn.le.2) then
            jper=1
          else if(jn.le.10) then
            jper=2
          else
            jper=3
          end if
c
          call bondlength(i,j,xc,r,r2)
          r=r*bohr
          r2=r**2
c
          if(bndorder(i,j).gt.0.7d0) then
c  Standard force constants: H-first row
            if(iper.eq.1.and.jper.eq.1) then
              bb=bhh
            else if(iper.eq.2.and.jper.eq.1) then
              bb=boh
            else if(iper.eq.2.and.jper.eq.2) then
              bb=boo
            else if(iper.eq.3.and.jper.eq.1) then
              bb=bsh
            else if(iper.eq.3.and.jper.eq.2) then
              bb=bso
            else if(iper.eq.3.and.jper.eq.3) then
              bb=bss
            end if
c
            fc=aa/(r-bb)**3
            write(inpf3,100) i,j,fc,r
  100 format('K',19X,'STRE',6X,2(i7,'.',2x),20X,2F15.5)         
            nq=nq+1
            nprim=nprim+1
          else
            fc=bndorder(i,j)
c  Scale the inverse distance coordinate by the squared distance
            write(inpf3,200) r2,i,j,fc,r
            nq=nq+1
            nprim=nprim+1
  200 format('K',F9.2,10X,'INVR',6X,2(i7,'.',2x),20X,2F15.5)
          end if
        end do
      end do
c
      end
c=======================================================================
      subroutine bendings(nat,inpf3,ian,xc,bndorder,nq,nprim) 
      implicit real*8 (a-h,o-z)
c  Arguments
c  INTENT(IN)
c  nat               = number of centers
c  inpf3             = output input file number (formatted)
c  ian(nat)          = array of nuclear charges
c  xc(3,nat)         = Cartesian coordinates (bohr)
c  bndorder(nat,nat) = bond order matrix
c  INTENT(OUT)
c  nq                = current number of internal coordinates
c  nprim             = current number of primitive coordinates
c  In addition, the routine writes the definition of the internal
c    coordinates and approximate diagonal force constants on unit inpf3
c
       
c  This routine assigns  the bending coordinates, calculates
c  approximate force constants and writes them to inpf3
c
      dimension ian(nat),xc(3,nat),bndorder(nat,nat)
      dimension u(3),v(3),w(3),z(3)
      parameter(zero=0.0d0,half=0.5d0,bohr=0.529177d0,
     &          pi=3.14159 26535 89793 24d0)
cpp
c      write(*,*) 'In bendings, nat=',nat 
c      do i=1,nat
c        write(*,998) (bndorder(i,j),j=1,nat)
c      end do
c  998 format (10f8.3)
      do i=1,nat
        in=ian(i)
        do j=1,nat
          if(i.eq.j) cycle
          if(bndorder(i,j).eq.zero) cycle
          do k=1,nat
            if(bndorder(j,k).eq.zero) cycle
            if(k.le.i.or.k.eq.j) cycle
c  Eliminate 3-rings
            if(bndorder(i,k).gt.half) cycle
c Calculate the angle. generate linear bends if too close to 180 deg
            kn=ian(k)
c            write(*,*) 'i,j,k,in,kn',i,j,k,in,kn
c            call flush6()
            call anggrad(nat,i,j,k,xc,th,.false.,xjunk) 
            if(th.lt.zero) th=-th
c  Exclude acute angles
            if(th.lt.0.2d0) cycle
            if(pi-th.gt.0.1d0) then
c  Estimate force constant
              xmult=bndorder(i,j)*bndorder(j,k)
              if(xmult.lt.0.1d0) xmult=0.1d0
              if(in.eq.1.or.kn.eq.1) then
                 fc=0.7d0*xmult
              else
                fc=1.0d0*xmult
              end if
c  Note that the central atom is given LAST
              write(inpf3,100) i,k,j,fc,th
cpp
c      write(*,'("bend",3i4,2f7.3)') i,k,j,bndorder(i,j),bndorder(j,k)
              nq=nq+1
              nprim=nprim+1
  100 format('K',19X,'BEND',6X,3(i7'.',2x),10X,2F15.5)         
            else
c  Generate two linear bends. Find an atom not identical with i,j,k, and
c    not collinear or approximately collinear with the i-j or j-k bonds
              distance=1.0d6
              mm=0
              do m=1,nat
                if(m.eq.i.or.m.eq.j.or.m.eq.k) cycle
                call anggrad(nat,i,j,m,xc,th,.false.,xjunk)
                if(th.gt.0.5d0.and.th.lt.2.6425926d0) then
                  call bondlength(j,m,xc,r,r2)
                  if(r.lt.distance) then
                    distance=r
                    mm=m
                  end if
                end if
              end do
              if(mm.gt.0) then
                m=mm
                fc=0.8d0*xmult 
c  This theta value was totally wrong
c  Determine the LIN1 and LIN2 values
                call anggrad(nat,i,j,m,xc,th1,.false.,xjunk)
                call anggrad(nat,k,j,m,xc,th2,.false.,xjunk)
                th=th1+th2-pi
                write(inpf3,200) i,k,j,mm,fc,th
                call vektor1(u,r,i,j,xa)
                call vektor1(v,r,m,j,xa)
                call normalv(u,v,w)
                call vektor1(u,r,k,j,xa)
                call normalv(u,v,z)
                th=(w(2)*z(3)-w(3)*z(2))*v(1)+(w(3)*z(1)-w(1)*z(3))*v(2)
     &              +(w(1)*z(2)-w(2)*z(1))*v(3)
                write(inpf3,300) i,k,j,mm,fc,th
  200 format('K',19X,'LIN1',6X,4(i7,'.',2x),2F15.5)         
  300 format('K',19X,'LIN2',6X,4(i7,'.',2x),2F15.5)         
                nq=nq+2
                nprim=nprim+2
              end if
            end if
          end do
        end do
      end do
c
      end 
c=======================================================================
      subroutine torsion(nat,inpf3,ian,xc,bndorder,nq,nprim) 
      implicit real*8 (a-h,o-z)
c  Arguments
c  INTENT(IN)
c  nat               = number of centers
c  inpf3             = output input file number (formatted)
c  ian(nat)          = array of nuclear charges
c  xc(3,nat)         = Cartesian coordinates (bohr)
c  bndorder(nat,nat) = bond order matrix
c  INTENT(OUT)
c  nq                = current number of internal coordinates
c  nprim             = current number of primitive coordinates
c  In addition, the routine writes the definition of the internal
c    coordinates and approximate diagonal force constants on unit inpf3
c
       
c  This routine assigns  the torsional coordinates, calculates
c  approximate force constants and writes them to inpf3
c  The intended torsional coordinate is a single coordinate for a 
c  rotation around a bond, an for n prinitive torsions, it is 
c  1/n(Sum(torsions). As the program normlizes the internal coordinates,
c  this is achived iof the overall scale factor is 1/sqrt(n)
c
      dimension ian(nat),xc(3,nat),bndorder(nat,nat)
      dimension iother(10),jother(10)
      parameter(zero=0.0d0,one=1.0d0,bohr=0.529177d0,
     &                    pi=3.14159 26535 89793 24d0)
      do i=2,nat
        do j=1,i
          if(bndorder(i,j).eq.0.0d0) cycle
c  i-j is a bond. Find out whether i and j have other bonds.
          ibond=0
          do k=1,nat
            if(k.eq.i.or.k.eq.j) cycle
            if(bndorder(i,k).eq.0.0d0) cycle
c  Check whether the k-i-j angle is nearly linear (<8.594 or >171.4 deg.)
            call anggrad(nat,k,i,j,xc,theta,.false.,xjunk)
            if(theta.lt.0.15d0.or.theta.gt.2.9916d0) cycle
            ibond=ibond+1
            iother(ibond)=k
          end do
          jbond=0
          do k=1,nat
            if(k.eq.i.or.k.eq.j) cycle
            if(bndorder(j,k).eq.0.0d0) cycle
c  Check whether i-j-k is nearly linear
            call anggrad(nat,i,j,k,xc,theta,.false.,xjunk)
            if(theta.lt.0.15d0.or.theta.gt.2.9916d0) cycle
            jbond=jbond+1
            jother(jbond)=k
          end do
c  Stop if i or j have no other connections (are terminal atoms)
          if(ibond.eq.0.or.jbond.eq.0) cycle
cpp
c          write(*,*) 'ibond,jbond=',ibond,jbond
          coeff=one/(dble(ibond*jbond))
          call torsangle(iother(1),i,j,jother(1),xc, tau)
cpp
c      write(*,*) 'T',iother(1),i,j,jother(1)
          fc=0.1d0*bndorder(i,j)**4
c          write(inpf3,100) coeff,one, iother(1),i,j,jother(1),fc,tau
c          nq=nq+1
c          nprim=nprim+1
c 100 format('K',F9.6,f10.5,'TORS',6X,4(i7,'.',2x),2F15.5)         
          istart=0
          do i1=1,ibond
            ia=iother(i1)
            do j1=1,jbond
              ja=jother(j1)
              if(i1.eq.j1) cycle
              istart=istart+1
              call torsangle(ia,i,j,ja,xc,tau)
cpp
c      write(*,*) 'TT',i1,i,j,j1
               if(istart.eq.1) then
                 write(inpf3,100) coeff,one, ia,i,j,ja,fc,tau
                 nq=nq+1
               else
                 write(inpf3,200) coeff,one, ia,i,j,ja,fc,tau
               end if
              nprim=nprim+1
  100 format('K',F9.6,f10.5,'TORS',6X,4(i7,'.',2x),2F15.5)         
  200 format('C',F9.6,f10.5,'TORS',6X,4(i7,'.',2x),2F15.5)         
            end do
          end do
        end do
      end do
c
      end
c=======================================================================
      subroutine torsangle(i,j,k,m,xc, angle)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0,pi=3.14159 26535 89793 24d0,small=1.0d-11)
      dimension xc(3,*),u(3),v(3),w(3),z(3),x(3)
      shftpt=0.5d0*pi
      call vektor1(u,r1,i,j,xc)
c  make this counterintuitive way to match the old BMAT routine
      call vektor1(v,r2,k,j,xc)
      call vektor1(w,r3,k,m,xc)
      call normalv (u,v,z)
      call normalv (w,v,x)
      co=scalar(z,x)
      u(1)=z(2)*x(3)-z(3)*x(2)
      u(2)=z(3)*x(1)-z(1)*x(3)
      u(3)=z(1)*x(2)-z(2)*x(1)
      si3=sqrt(u(1)**2+u(2)**2+u(3)**2)
      co2=scalar(u,v)
      s=arc1(-co,si3)
      if (co2.lt.zero) s=-s
      if (s.gt.shftpt) s=s-2.0d0*pi
      angle=-s
c
cpp
c      write(*,'(4i4,"co,si3,angle",3f12.6)') i,j,k,m,co,si3,angle
      end
c=======================================================================
      subroutine outofplane(nat,inpf3,ian,xc,bndorder,nq,nprim) 
      implicit real*8 (a-h,o-z)
c  Arguments
c  INTENT(IN)
c  nat               = number of centers
c  inpf3             = output input file number (formatted)
c  ian(nat)          = array of nuclear charges
c  xc(3,nat)         = Cartesian coordinates (bohr)
c  bndorder(nat,nat) = bond order matrix
c  INTENT(OUT)
c  nq                = current number of internal coordinates
c  nprim             = current number of primitive coordinates
c  In addition, the routine writes the definition of the internal
c    coordinates and approximate diagonal force constants on unit inpf3
c
       
c  This routine assigns  the bending coordinates, calculates
c  approximate force constants and writes them to inpf3
c
      dimension ian(nat),xc(3,nat),bndorder(nat,nat)
      dimension jother(10),u(3),v(3),w(3),z(3)
      parameter(zero=0.0d0,one=1.0d0,bohr=0.529177d0,
     &      small=1.0d-11,pi=3.14159 26535 89793 24d0)
      do i=1,nat
        jj=0
        jjold=0
        do j=1,nat
          if(j.eq.i) cycle
          if(bndorder(i,j).gt.0.0d0) then
            jj=jj+1
            jother(jj)=j
          end if
        end do
c  jj=the number of bonded neighbors of i, iother(1..jj): these atoms
cpp
c        write(*,*) 'iat, bonded to',i,(jother(j),j=1,jj)
        call flush6()
        if(jj.eq.3) then
c  look at the local geometry
          call vektor1(u,r,jother(2),i,xc)
          call vektor1(v,r,jother(3),i,xc)
          call vektor1(z,r,jother(1),i,xc)
          w(1)=u(2)*v(3)-u(3)*v(2)
          w(2)=u(3)*v(1)-u(1)*v(3)
          w(3)=u(1)*v(2)-u(2)*v(1)
          d=scalar(w,z)
          d1=d/sqrt(scalar(w,w))
          theta=asin(d1)
cpp
c          write(*,*) 'Out of plane, d=',d
          if(abs(d).gt.0.9d0) cycle
c  divide by 3 - 3 similar values
          coeff=((one-abs(d))**3-0.001d0)/3.0d0
          fc=0.3d0
          write(inpf3,100) coeff,one, jother(2),jother(3),jother(1),i,
     &                    fc,theta
          nq=nq+1
          nprim=nprim+3
  100 format('K',F9.6,f10.5,'OUT ',6X,4(i7,'.',2x),2F15.5)         
c  Form symmetrical contribution
          do k=1,3
            tmp=u(k)
            u(k)=v(k)
            v(k)=z(k)
            z(k)=tmp
          end do
          w(1)=u(2)*v(3)-u(3)*v(2)
          w(2)=u(3)*v(1)-u(1)*v(3)
          w(3)=u(1)*v(2)-u(2)*v(1)
          d=scalar(w,z)
          d1=d/sqrt(scalar(w,w))
          theta=asin(d1)
          write(inpf3,200) coeff,one, jother(3),jother(1),jother(2),i,
     &                    fc,theta
  200 format('C',F9.6,f10.5,'OUT ',6X,4(i7,'.',2x),2F15.5)         
          do k=1,3
            tmp=u(k)
            u(k)=v(k)
            v(k)=z(k)
            z(k)=tmp
          end do
          w(1)=u(2)*v(3)-u(3)*v(2)
          w(2)=u(3)*v(1)-u(1)*v(3)
          w(3)=u(1)*v(2)-u(2)*v(1)
          d=scalar(w,z)
          d1=d/sqrt(scalar(w,w))
          theta=asin(d1)
          write(inpf3,200) coeff,one, jother(1),jother(2),jother(3),i,
     &                    fc,theta
        end if
c  Possibly add out-of-plane for more than 3 neighbors
      end do
c
      end
c=======================================================================
      subroutine intcvalues(nat,xc,shftpt,nq,qq)
      implicit real*8 (a-h,o-z)
c  This routine generates the internal coordinate values
c  BMAT calculates these, too. However, this routine also reads the previous
c  values of the PRIMITIVE internal coordinates, and corrects the torsions 
c  by multiples of 2*pi
      common /tapeneu/ inp,inp2,iout,ipun,ix,icond,inpf3,inpf4
      character*110 line
      character*1 start
      character*4 type
      dimension xc(3,nat),qq(*)
      dimension u(3),v(3),w(3),z(3)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,
     &           pi=3.14159 26535 89793 24d0)
      rewind inpf3
      rewind inpf4
      nq=0
cpp
c      write(*,*) 'In intcvalues'
c  This shoul not be neceessary
       call zeroit(qq,nq)
 100  continue
        read(inpf3,'(a110)',end=1000) line
        backspace inpf3
        read(inpf3,200,end=1000) start,scale,coeff,type,xia,xja,xka,xla,
     &                           fc,value
 200    format(a1,f9.5,f10.6,a4,6x,4f10.4,2f15.8)
        if(scale.eq.zero) scale=one
        if(coeff.eq.zero) coeff=one
        ia=nint(xia)
        ja=nint(xja)
        ka=nint(xka)
        la=nint(xla)
        if(type.eq.'STRE')  then
          call bondlength(ia,ja,xc,r,r2)
          primq=r
        else if(type.eq.'BEND') then
          call anggrad(nat,ia,ka,ja,xc,theta,.false.,xjunk)
          primq=theta
        else if (type.eq.'INVR') then
          call bondlength(ia,ja,xc,r,r2)
          primq=one/r
        else if(type.eq.'LIN1') then
c  Determine the LIN1 values
          i=ia
          k=ja
          j=ka
          m=la
          call anggrad(nat,i,j,m,xc,th1,.false.,xjunk)
          call anggrad(nat,k,j,m,xc,th2,.false.,xjunk)
          th=th1+th2-pi
          primq=th
        else if(type.eq.'LIN2') then
          i=ia
          k=ja
          j=ka
          m=la
          call vektor1(u,r,i,j,xa)
          call vektor1(v,r,m,j,xa)
          call normalv(u,v,w)
          call vektor1(u,r,k,j,xa)
          call normalv(u,v,z)
          th=(w(2)*z(3)-w(3)*z(2))*v(1)+(w(3)*z(1)-w(1)*z(3))*v(2)
     &        +(w(1)*z(2)-w(2)*z(1))*v(3)
          primq=th
        else if(type.eq.'TORS') then
          call torsangle(ia,ja,ka,la,xc,tau)
          primq=tau
        end if 
c  If the coordinate value changes by 2*Pi, correct it
        if(abs(primq-value).gt.abs(primq+two*pi-value))
     &                primq=primq+two*pi
        if(abs(primq-value).gt.abs(primq-two*pi-value))
     &                primq=primq-two*pi
        write(inpf4,200) start,scale,coeff,type,xia,xja,xka,xla,
     &                           fc,primq
c        write(*,200) start,scale,coeff,type,xia,xja,xka,xla,
c     &                           fc,primq
        if(start.eq.'K'.or.start.eq.'k') then
          if(nq.ge.1) then
c  Multiply the previous coordinate by the normalization factor
            snorm=one/sqrt(snorm)
            qq(nq)=qq(nq)*snorm
          end if
          nq=nq+1
          qq(nq)=scale*coeff*primq
          snorm=coeff**2
        else
          snorm=snorm+coeff**2
          qq(nq)=qq(nq)+scale*coeff*primq
        end if  
      go to 100
 1000 continue
c  Multiply the last coordinate by the normalization factor
c  This is needed only to keep compatibility with the old machb
      snorm=one/sqrt(snorm)
      qq(nq)=qq(nq)*snorm
c  Swap files (b/c it is not possible to write in the middle of a sequential file
      itmp=inpf4
      inpf4=inpf3
      inpf3=itmp
      rewind inpf3
      rewind inpf4
c
cpp
c      write(*,*) 'calculated values of the composite internal ',
c     &           'coordinates in intcvalues'
c      do k=1,nq
c        write(*,'(i5,f12.6)') k,qq(k)
c      end do
      end
c=======================================================================
