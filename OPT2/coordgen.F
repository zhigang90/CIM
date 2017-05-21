      subroutine coordgen(nat,iout,x,iuniq,atname,  xc,ian,nq,nprim)
      implicit real*8 (a-h,o-z)
      include "common/big"
      include "common/cpi"
      common /tapeneu/ inp,inp2,ioux,ipun,ix,icond,inpf3,inpf4
      dimension x(3,nat),iuniq(nat),xc(3,nat),ian(nat)
      character(*) atname(*)
      character*110 line
      logical writ
c
      istart=icorr(0)
!   MOLPRO nuclear info was converted to PQS (still in Bohr units!) 
!  in the calling program
c      call ConvertCart(nat,iout,x,iuniq,atname,  xc,ian)
      membo=icorr(nat**2)
      call bondorder(nat,ian,xc,q(membo),ierr)
      inpf3=33
      inpf4=34
c  The following routine opens an ASCII (units 1-9 or 22 and higher)
c  or binary file. If the file name=blank then it is a scratch file
c  It returns a file size ??
      call assgn(inpf3,'intcoors',isize,0)
      call assgn(inpf4,'intcoors1',isize,0)
c
      call newintc(nat,inpf3,ian,xc,q(membo),nq,nprim)
cpp
      rewind inpf3
 90   continue
       read(inpf3,'(a110)',end=99) line
       write(6,'(a110)') line 
       go to 90
 99   continue
       call corlsr(istart)
c
      end
c=======================================================================
      subroutine intforc1(nat,xc,grad,nq,nprim,
     &                    ibcode,ibcontr,qq,bmat,gmat,proj,phi)
c  Arguments:
c  INTENT(IN)
c  nat        = number of atoms (centers)
c  xc(3,nat)  = Cartesians (will be scaled!)
c  grad(3,nat)= analytic gradient (will be sign-reversed and scaled!)
c  nq         = number of (coimposite) internal coordinates
c  nprim      = number of primitive internal coordinates
c
c  INTENT(OUT)
c  ibcode     = encoding of the B matrix
c  ibcontr(*,i) gives the number of atoms and the atoms in the ith row
c               of the B matrix (see BREAD)
c  qq         = values of the composite internal coordinates
c  bmat       = the non-zero elements of the B matrix
c  gmat       = holds the G matrix and its inverse
c  proj       = projector to the non-null space of G
c  phi        = internal forces in aJ/A (stretches) or aJ/rad
      implicit real*8 (a-h,o-z)
      logical writ
      dimension xc(3,nat),grad(3,nat),ibcode(6,*),ibcontr(20,*),
     &          qq(nq),bmat(54,nq),gmat(nq,nq),proj(nq,nq),phi(nq)
      include "common/big"
      common /tapeneu/ inp,inp2,ioux,ipun,ix,icond,inpf3,inpf4
      parameter (zero=0.0d0,one=1.0d0,pi=3.14159 26535 89793 24d0)
c
      ibase=icorr(0)
      maxcard=nprim
      mq=nq
      writ=.false.
cpp
      iout= 6
      rewind inpf3
c  Convert the current coordinates to angstrom units for BMAT
      toang=get_scalar('TOANG')
      ajoule=get_scalar('JOULE')*1.0d18 ! contrary to its name, it is 1 Eh in aJ
      call dscal_x(3*nat,toang,xc,1) 
      call dscal_x(3*nat,-ajoule/toang,grad,1)
      call bread(nat,inpf3,6,nq,nprim,  ibcode,maxcard,ibcontr,mq,writ)
c  Reserve memory for the B matrix and for the internal coord. values
c      membmat=icorr(54*nq)
      shftpt=0.5d0*pi
      call machbnew(nat,xc,nq,.false.,shftpt,  
     &              ibcode,ibcontr,bmat,qq)
cpp
c       write(*,*) 'Composite internal coordinates in intforc1'
c       write(*,'(5f12.6)') (qq(k),k=1,nq)
c      call printbmat(nq,ibcontr,qq,bmat)
cpp
       call formgmat(nq,3*nat,bmat,ibcontr,gmat)
c       call prntmatn(6,nq,nq,nq,gmat,'G matrix')
       memw=icorr(nq)
       call bxv(nq,3*nat,bmat,ibcontr,grad,q(memw))
       tol=1.0d-6
       memu=icorr(nq**2)
       memx=icorr(nq)
       memy=icorr(nq)
       call geninv(nq,gmat,det,tol,q(memu),  q(memx),q(memy),ninv,proj) 
cpp
       write(6,*) 'Number of internal and independent coordinates',
     &             nq,ninv
c       call prntmatn(6,nq,nq,nq,gmat,'G inverse')
c       call prntmatn(6,nq,nq,nq,proj,'Projector')
c
c  Check if the coordinate system is complete
c  ninv is the number of non-zero eigenvalues of G
       if(ninv.lt.(3*nat-6).and.nat.gt.2) then
        write(*,*) 'Too few independent coordinates were found',
     & ' geomety optimization in internal coordinates may fail', ninv
c        call fehler
       end if
c       if(ninv.lt.(3*nat-6)) then
c         write(*,*) 'Warning: fewer than 3N-6 independent internals',
c     &  ' Geometry optimization may fail unless the  molecule is linear'
c       end if 
       call dgemv_x('N',nq,nq,one,gmat,nq,q(memw),1,zero,phi,1)
c       call prntmatn(6,1,nq,1,phi,'Internal forces')
       call corlsr(ibase)
cpp
c      stop 'intforc1'
      end
c=======================================================================
      subroutine printbmat(nq,ibcontr,qq,bmat)
      implicit real*8 (a-h,o-z)
      dimension ibcontr(20,nq),qq(nq),bmat(54,nq)
      write(*,*) 
     &    'Internal coordinates and non-zero elements of the B matrix'
      do i=1,nq
        na=ibcontr(2,i)
        write(*,130) i,qq(i),
     &  (ibcontr(k+2,i),bmat(3*k-2,i),bmat(3*k-1,i),bmat(3*k,i),k=1,na)
  130 format(i4,f10.5,i4,3f10.5,2(i4,3f10.5),(/,14x,3(i4,3f10.5)))
      end do
      end
c=======================================================================

      subroutine ConvertCart(nat,iout,x,iuniq,atname,  xc,ian)
      implicit real*8 (a-h,o-z)
c  aRGUMENTS
c integer,     intent(in) ::  nat        ! number of nuclei (centers)  
c integer,     intent(in) ::  iout       ! output file name
c real*8,      intent(in) ::  x(3,nat)   ! nuclear coordinates in bohr in the first step
c integer,     intent(in) ::  iuniq(nat) ! pointer to the names of atoms
c character*(*), intent(in) ::  atname(*) ! the names of unique atoms
c
c real*8,      intent(out)::  xc(3,nat)  ! Cartesians in atomic units
c integer,     intent(out)::  ian(nat)   ! atomic numbers
C  converts MOLPRO coordinates and atom names to PQS-like quantities
      dimension x(3,nat),xc(3,nat),iuniq(nat),ian(nat)
      character(*) atname(*)
      character*2 atn,elmt(0:99)
      DATA ELMT /'X ','H ','HE','LI','BE','B ','C ','N ',
     1           'O ','F ','NE','NA','MG','AL','SI',
     2           'P ','S ','CL','AR','K ',
     1           'CA','SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     2           'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR','NB',
     3           'MO','TC','RU','RH','PD','AG','CD','IN','SN','SB',
     4           'TE','I ','XE','CS','BA','LA','CE','PR','ND','PM','SM',
     5           'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF',
     6           'TA','W ','RE','OS','IR','PT','AU','HG','TL','PB','BI',
     7           'PO','AT','RN','FR','RA','AC','TH','PA','U ','NP','PU',
     8           'AM','CM','BK','CF','ES'/
c      write(*,*) 'Unique atom names=',atname(1),atname(2)
c      write(iout,*) 'Current atom names and coordinates'
      do iat=1,nat
        atn=atname(iuniq(iat))(1:2)
c        write(*,*) 'Atname=',atn
        do k=1,3
          xc(k,iat)=x(k,iat)
        end do
        natn=0
        do iel=0,99
          if(atn.eq.elmt(iel)) then
            natn=iel
            exit
          end if
        end do
        ian(iat)=natn
c        write(iout,100) atn,ian(iat),(xc(k,iat),k=1,3)
  100 format(a2,i3,2x,3f12.6)
       end do
!
       end
c======================================================================
c
      subroutine bondorder(natoms,ian,xc,bndorder,ierr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms a list of bond orders from the interatomic distances
C  and a list of "standard" single bond distances as contained
C  from the array BndLen.
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  cartesian coordinates
C             considered to be bonding
C  bndorder-  on exit contains the bond orders
C  IErr    -  error flag   0 - success
C                         -1 - no default bond length available
C
C
      REAL*8 XC(3,NAtoms),bndorder(NAtoms,NAtoms)
      INTEGER IAN(NAtoms)
C
      REAL*8 BndLen(92)      ! homonuclear bond distances in au
C
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0)
c                    H-H  He-He Li-Li Be-Be  B-B   C-C   N-N   O-O
      DATA BndLen / 1.42, 4.00, 5.10, 4.16, 3.21, 2.93, 2.66, 2.50,
c
c                    F-F  Ne-Ne Na-Na Mg-Mg Al-Al Si-Si  P-P   S-S
     $              2.12, 6.00, 5.86, 5.67, 5.01, 4.35, 4.16, 3.87,
c
c                   Cl-Cl Ar-Ar  K-K  Ca-Ca Sc-Sc Ti-Ti  V-V  Cr-Cr
     $              3.78, 8.00, 7.37, 6.58, 5.44, 5.00, 4.62, 4.42,
c
c                   Mn-Mn Fe-Fe Co-Co Ni-Ni Cu-Cu Zn-Zn Ga-Ga Ge-Ge
     $              4.40, 4.38, 4.38, 4.36, 4.35, 4.72, 4.72, 4.62,
c
c                   As-As Se-Se Br-Br Kr-Kr Rb-Rb Sr-Sr  Y-Y  Zr-Zr
     $              4.58, 4.44, 4.35, 9.00, 8.16, 7.22, 6.12, 5.50,
c
c                   Nb-Nb Mo-Mo Tc-Tc Ru-Ru Rh-Rh Pd-Pd Ag-Ag Cd-Cd
     $              5.06, 4.88, 4.80, 4.70, 4.72, 4.84, 5.06, 5.33,
c
c                   In-In Sn-Sn Sb-Sb Te-Te  I-I  Xe-Xe Cs-Cs Ba-Ba
     $              5.67, 5.30, 5.33, 5.18, 5.03, 9.00, 8.88, 7.48,
c
c                   La-La Ce-Ce Pr-Pr Nd-Nd Pm-Pm Sm-Sm Eu-Eu Gd-Gd
     $              6.40, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00,
c
c                   Tb-Tb Dy-Dy Ho-Ho Er-Er Tm-Tm Yb-Yb Lu-Lu Hf-Hf
     $              6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 5.44,
c
c                   Ta-Ta  W-W  Re-Re Os-Os Ir-Ir Pt-Pt Au-Au Hg-Hg
     $              5.06, 4.92, 4.84, 4.76, 4.76, 4.88, 5.06, 5.44,
c
c                   Tl-Tl Pb-Pb Bi-Bi Po-Po At-At Rn-Rn Fr-Fr Ra-Ra
     $              5.86, 5.82, 5.75, 6.00, 6.00, 6.00, 6.00, 6.00,
c
c                   Ac-Ac Th-Th Pa-Pa  U-U
     $              6.00, 6.00, 6.00, 6.00 /
C
C
C  set error flag
C
      IErr = -1
c
c      write(*,*) 'Cartesians'
c      do i=1,natoms
c        write(*,'(i3,3f12.6)') ian(i),(xc(k,i),k=1,3)
c      end do
c
      bndorder(1,1)=zero
      DO I=2,NAtoms
        bndorder(i,i)=zero
        DO J=1,I-1
          ii=i
          jj=j
c  order the atoms so that jj has the smaller atomic number
          if(ian(i).lt.ian(j)) then
            ii=j
            jj=i
          end if
c  Do not assign bond order to dummy atoms
          if(ian(ii).eq.0.or.ian(jj).eq.0) cycle
C  get interatomic distance I-J
          dist=sqrt((xc(1,i)-xc(1,j))**2+(xc(2,i)-xc(2,j))**2+
     &              (xc(3,i)-xc(3,j))**2)
C
C  get "standard" I-J bond length from BndLen table as geometric mean
C
          bondL = SQRT(BndLen(IAN(I))*BndLen(IAN(J)))
C
c  Estimate the bond order
cpp
c      write(*,999) i,j,ian(i),ian(j),bondL,dist
c  999 format('Atoms=',2i5,' Atomic numbers=',2i5,'bondL,dist',2f10.4)
cpp
          bo=(bondl/dist)**3
          bo1=bo
          if(bo.lt.0.4d0) bo1=zero
          if(ian(jj).eq.1.and.ian(ii).ge.7.and.ian(ii).le.9) then
            if(bo.gt.0.07d0) bo1=bo
          end if
          bndorder(i,j)=bo1
          bndorder(j,i)=bo1
c
        end do
      end do
C
C  Check for and remove "excess connectivity"
C
      CALL ChkXSC(NAtoms,XC,bndorder)
C
C  Now check that the connectivity matrix forms a "closed loop"
C  i.e., that all atoms are connected to the whole
C
      CALL ChkCMPLT(NAtoms,XC,bndorder,IErr)
C
cpp
c      do i=1,natoms
c        write(*,200) (bndorder(i,j),j=1,natoms)
c      end do
 200  format(5f8.3,2x,5f8.3)
c      stop 2013
      end
c =================================================================
c
      SUBROUTINE CheckLinear(NAtoms,XC,bndorder)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for and removes "excess connectivity"
C  Principally to prevent in the near-linear arrangement A--B--C
C  that if A is connected to B and B is connected to C then
C  A is not also connected to C
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  bndorder-  bond order matrix (may be modified on exit)
C
C
      DIMENSION XC(3,NAtoms),bndorder(NAtoms,NAtoms)
c
      PARAMETER (zero=0.0d0,TollZero=0.002d0)
C
C
C  First find a connectivity
C
      DO 20 I=2,NAtoms
      DO 20 J=1,I
      IF(bndorder(I,J).gt.zero) THEN
C
C  I and J are connected
C  What else is J connected to?
C
        DO 10 K=J+1,I
        IF(bndorder(K,J).gt.zero) THEN
C
C  J and K are connected
C  Is K also connected to I?
C
          If(bndorder(I,K).gt.zero) Then
C
C  check the angle I-J-K
C  if near-linear, remove longest bond
C
            CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,xjnk)
            If(Th.LT.TollZero) Then
              DIJ = SQRT( (XC(1,I) - XC(1,J))**2 +
     $                    (XC(2,I) - XC(2,J))**2 +
     $                    (XC(3,I) - XC(3,J))**2 )
              DIK = SQRT( (XC(1,I) - XC(1,K))**2 +
     $                    (XC(2,I) - XC(2,K))**2 +
     $                    (XC(3,I) - XC(3,K))**2 )
              If(DIJ.GT.DIK) Then
                bndorder(I,J) = zero
                bndorder(J,I) = zero
              Else
                bndorder(I,K) = zero
                bndorder(K,I) = zero
              EndIf
            EndIf
          EndIf
        ENDIF
 10     CONTINUE
cc
      ENDIF
 20   CONTINUE
C
      END
c =====================================================================
c
      SUBROUTINE ChkCMPLT(NAtoms,XC,bndorder,IErr)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks that the connectivity matrix forms a "closed loop"
C  with all atoms connected
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
c  bndorder- Bond order matrix
C  IErr    -  error flag on exit
C              0 - OK; -1 - problem
C
C
      DIMENSION XC(3,NAtoms),bndorder(NAtoms,NAtoms)
c .........................................................
c -- F90 dynamically allocated memory
      INTEGER I1(NAtoms),I2(NAtoms),I3(NAtoms)
c .........................................................
C
      Parameter (zero=0.0d0,Big=1000.0d0)
C
C
      IErr = 0      ! should not be any problems
c
      num = 0
      CALL IZeroIT(I1,NAtoms)
C
C  First find a connectivity
C
      DO 10 I=2,NAtoms
      DO 10 J=1,I
      If(bndorder(I,J).gt.zero) GO TO 11
 10   CONTINUE
C
C  should not get here unless NOTHING is bonded
C
      If(NAtoms.EQ.2) Then
        bndorder(1,2) = 0.1d0
        bndorder(2,1) = 0.1d0
      Else
        IErr = -1
      EndIf
      RETURN
C
 11   CONTINUE
      num = num+1
      I1(I) = 1
      I2(num) = I
      num = num+1
      I1(J) = 1
      I2(num) = J
C
C  now search among known connected atoms for further atoms
C  connected to any of the existing atoms
C
      num0 = 0
 20   CONTINUE
      numS = num        ! store number of atoms in group
      DO 30 inum=num0+1,numS
      I = I2(inum)
      DO 25 J=1,NAtoms
      If(bndorder(I,J).gt.zero.AND.I1(J).EQ.0) Then
c
c -- found a new atom not previously connected to existing group
        num = num+1
        I1(J) = 1
        I2(num) = J
      EndIf
 25   CONTINUE
 30   CONTINUE
C
C  check if num has changed
C
      num0 = numS
      If(num.GT.numS.AND.num.LT.NAtoms) GO TO 20
C
C  At this point we are finished
C  Either everything is connected, in which case we are done OR
C  not all atoms are connected
C
      If(num.EQ.NAtoms) RETURN
C
C  find which atoms are not connected and store them in I3
C
      num3 = 0
      DO 40 I=1,NAtoms
      If(I1(I).EQ.0) Then
        num3 = num3+1
        I3(num3) = I
      EndIf
 40   CONTINUE
C
C  now find the shortest distance between the atoms in the
C  two groups and connect them
C
      DistM = Big
      DO 50 inum=1,num
      I = I2(inum)
      DO 45 jnum=1,num3
      J = I3(jnum)
C
C  get interatomic distance I-J
C
      Dist = SQRT( (XC(1,I) - XC(1,J))**2 +
     $             (XC(2,I) - XC(2,J))**2 +
     $             (XC(3,I) - XC(3,J))**2 )
c
      If(Dist.LT.DistM) Then
        DistM = Dist
        IMin = I
        JMin = j
      EndIf
c
 45   CONTINUE
 50   CONTINUE
C
C  connect the new pair
      bndorder(IMin,JMin) = 0.05d0
      bndorder(JMin,IMin) = 0.051d0
C
C  check for other atoms (most likely symmetry-related) that
C  are the same distance as (IMin,JMin) pair and connect them
C
      DO 60 inum=1,num
      I = I2(inum)
      DO 55 jnum=1,num3
      J = I3(jnum)
C
C  get interatomic distance I-J
C
      Dist = SQRT( (XC(1,I) - XC(1,J))**2 +
     $             (XC(2,I) - XC(2,J))**2 +
     $             (XC(3,I) - XC(3,J))**2 )
c
      If(Dist.EQ.DistM) Then
        bndorder(I,J) = 0.05d0
        bndorder(J,I) = 0.05d0
      EndIf
c
 55   CONTINUE
 60   CONTINUE
C
C  now go back and continue the search
C
      num = num+1
      I2(num) = JMin
      I1(JMin) = 1
      GO TO 20
C
      END
c =====================================================================
c
      SUBROUTINE AngGRAD(NAtoms,I,J,K,XC,Th,grd,G)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculate the value and constraint normal for
C  the bond angle I-J-K
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  I       -  first atom in bond angle
C  J       -  central atom in bond angle
C  K       -  third atom in bond angle
C  XC      -  Cartesian coordinates
C  Th      -  on exit contains bond angle
C  grd     -  logical flag for calculating constraint normal
C              .true.  -  calculate angle gradient
C              .false. -  skip gradient calculation
C  G       -  on exit contains angle gradient (if calculated)
C
C
      REAL*8 XC(3,NAtoms),G(3,NAtoms)
      LOGICAL grd
C
      PARAMETER (One=1.0d0,Two=2.0d0,small=1.0d-6)
C
C
      XIJ = XC(1,I) - XC(1,J)
      XIK = XC(1,I) - XC(1,K)
      XJK = XC(1,J) - XC(1,K)
      YIJ = XC(2,I) - XC(2,J)
      YIK = XC(2,I) - XC(2,K)
      YJK = XC(2,J) - XC(2,K)
      ZIJ = XC(3,I) - XC(3,J)
      ZIK = XC(3,I) - XC(3,K)
      ZJK = XC(3,J) - XC(3,K)
      RIJ = SQRT(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ)
      RIK = SQRT(XIK*XIK + YIK*YIK + ZIK*ZIK)
      RJK = SQRT(XJK*XJK + YJK*YJK + ZJK*ZJK)
c
      D1 = Two*RIJ*RJK
      D2 = (RIJ*RIJ + RJK*RJK - RIK*RIK)
c
      CosTh = D2/D1
      If(Abs(CosTh).GT.One) CosTh = SIGN(One,CosTh)
      SinTh = SQRT(One - CosTh*CosTh)
      Th = ACOS(CosTh)
c
      If(.NOT.grd) RETURN
C
C  .............................................................
C    ** WARNING  **
C  If the three atoms are linear there are problems with the
C  angle bend derivative.  Skip derivative evaluation
C
      If(SinTh.LT.small) RETURN
C  .............................................................
C
      CALL ZeroIT(G,3*NAtoms)
c
      DCB = (Two/SinTh)/(D1*D1)
      RJKIJ = RJK/RIJ
      RIJJK = RIJ/RJK
c
      G(1,I) = -DCB*( D1*XJK + D2*XIJ*RJKIJ )
      G(1,J) =  DCB*( D1*(XJK-XIJ) - D2*(XJK*RIJJK - XIJ*RJKIJ) )
      G(1,K) =  DCB*( D1*XIJ + D2*XJK*RIJJK )
      G(2,I) = -DCB*( D1*YJK + D2*YIJ*RJKIJ )
      G(2,J) =  DCB*( D1*(YJK-YIJ) - D2*(YJK*RIJJK - YIJ*RJKIJ) )
      G(2,K) =  DCB*( D1*YIJ + D2*YJK*RIJJK )
      G(3,I) = -DCB*( D1*ZJK + D2*ZIJ*RJKIJ )
      G(3,J) =  DCB*( D1*(ZJK-ZIJ) - D2*(ZJK*RIJJK - ZIJ*RJKIJ) )
      G(3,K) =  DCB*( D1*ZIJ + D2*ZJK*RIJJK )
c
      RETURN
      END
c=============================================================
      subroutine izeroit (a,n)
      implicit integer(a-z)
      dimension a(n)
      do 100 i=1,n
       a(i) = 0
 100  continue
      return
      end
c==============================================================
c
      SUBROUTINE ChkXSC(NAtoms,XC,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Checks for and removes "excess connectivity"
C  Principally to prevent in the near-linear arrangement A--B--C
C  that if A is connected to B and B is connected to C then
C  A is not also connected to C
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XC      -  Cartesian coordinates
C  IC      -  connectivity matrix (may be modified on exit)
C
C
      DIMENSION XC(3,NAtoms),IC(NAtoms,NAtoms)
c
      PARAMETER (TollZero=0.002d0)
C
C
C  First find a connectivity
C
      DO 20 I=2,NAtoms
      DO 20 J=1,I
      IF(IC(I,J).EQ.1) THEN
C
C  I and J are connected
C  What else is J connected to?
C
        DO 10 K=J+1,I
        IF(IC(K,J).EQ.1) THEN
C
C  J and K are connected
C  Is K also connected to I?
C
          If(IC(I,K).EQ.1) Then
C
C  check the angle I-J-K
C  if near-linear, remove longest bond
C
            CALL AngGRAD(NAtoms,I,J,K,XC,Th,.false.,jnk)
            If(Th.LT.TollZero) Then
              DIJ = SQRT( (XC(1,I) - XC(1,J))**2 +
     $                    (XC(2,I) - XC(2,J))**2 +
     $                    (XC(3,I) - XC(3,J))**2 )
              DIK = SQRT( (XC(1,I) - XC(1,K))**2 +
     $                    (XC(2,I) - XC(2,K))**2 +
     $                    (XC(3,I) - XC(3,K))**2 )
              If(DIJ.GT.DIK) Then
                IC(I,J) = 0
                IC(J,I) = 0
              Else
                IC(I,K) = 0
                IC(K,I) = 0
              EndIf
            EndIf
          EndIf
        ENDIF
 10     CONTINUE
cc
      ENDIF
 20   CONTINUE
C
      RETURN
      END
c=============================================================
      subroutine diagfcinv(nq,hess,proj)
      implicit real*8 (a-h,o-z)
c  Reads the internal coordinate definition file and gets the initrial diagonal
c   force constants and puts them in a matrix H. Finally forms P(PHP)-P
      character*1 startc
      dimension hess(nq,nq),proj(nq,nq)
      include "common/big"
      common /tapeneu/ inp,inp2,ioux,ipun,ix,icond,inpf3,inpf4
      parameter(zero=0.0d0,one=1.0d0)
c  This routine reads the initial diagonal force constants from ifile
c  inverts them and transforms them to redundant internal coordinates
      ifile=inpf3
      rewind ifile
      ibase=icorr(0)
      memsq1=icorr(nq**2)
      memsq2=icorr(nq**2)
      call zeroit(q(memsq1),nq**2)
      k0=0
      kd=0
 100  continue
      read(ifile,200,end=500) startc,fc 
 200  format(a1,69x,f15.5)
      if(startc.eq.'K'.or.startc.eq.'k') then
        k0=k0+1
c  Put the force constants in the diagonal
        q(memsq1+kd)=fc
cpp
c        write(*,*) kd+1,fc
        kd=kd+nq+1
      end if
      go to 100
cpp
 500  continue
c      write(*,*) 'Number of f.c., nq=',k0,nq
      call dgemm('N','N',nq,nq,nq,   one,proj,nq,q(memsq1),nq,
     &            zero,q(memsq2),nq)
c      call prntmatn(6,nq,nq,nq,q(memsq2),'PD')
c      call prntmatn(6,nq,nq,nq,proj,'Projector')
      call dgemm('N','N',nq,nq,nq,   one,q(memsq2),nq,proj,nq,
     &           zero,q(memsq1),nq)
      memu=icorr(nq**2)
      memx=icorr(nq)
      memy=icorr(nq)
      tol=0.02d0
      call geninv(nq,q(memsq1),det,tol,q(memu),
     &            q(memx),q(memy),ninv,q(memsq2))
      call dgemm('N','N',nq,nq,nq,   one,proj,nq,q(memsq1),nq,
     &            zero,q(memsq2),nq)
      call dgemm('N','N',nq,nq,nq,   one,q(memsq2),nq,proj,nq,
     &            zero,hess,nq)
c  hess should contain the correctly projected Hessian, built from the
c  initial diagonal approximation
cpp
c      call prntmatn(6,nq,nq,nq,hess,'Projected diagonal Hessian')
      end
c=============================================================
      subroutine  relaxint(nq,phi,hinv,iout,iprint,
     &                     coordmax,dq,edecr,dismax,icmx)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0)
c  Applies H**(-1) to the internal forces
c  Arguments
c  INTENT(IN)
c  nq          = number (redundant) internal coordinates
c  phi(nq)     = internal forces
c  hinv(nq,nq) = inverse Hessian
c  iout        = output unit number
c  iprint      = print level. If 0, no printing, if 1, little
c  coordmax    = maximum displacement in internal coordinates, usually 0.1
c
c  INTENT(OUT)
c  dq(nq)      = change in the internal coordinates
      dimension hinv(nq,nq),intforc(nq),dq(nq)
      call dgemv_x('N',nq,nq,one,hinv,  nq,phi,1,zero,dq,  1)    
      edecraJ=-half*ddot(nq,phi,1,dq,1)
      attojoule=get_scalar('JOULE')*1.0d18  ! this is 1 Eh in aJ
      edecr=edecraJ/attojoule
      call absmax(nq,dq,icmx,dismax)
      if(iprint.gt.0) write(iout,100) dismax,icmx,edecraJ,edecr 
  100 format('Maximum internal displacement=',f10.6,' along',
     &       i5,'Energy lowering=',f10.6,' aJ',f11.7,' Eh')
      if(iprint.ge.2) then
      call prntmatn(iout,1,nq,1,dq,'Internal displacements')
c      write(iout,*) 'Internal displacements'
c      write(iout,100) (dq(k),k=1,nq)
c 100  format(5f10.6)
c
c  Reverse the step direction if it violates the descent condition
      scfac=one
      write(6,*) 'Relaxint: Estimated energy decrease in au=',edecr
      if(dismax.gt.coordmax) scfac=coordmax/dismax
      if(edecr.gt.zero) scfac=-one
      if(scfac.ne.one) call mult(dq,scfac,nq)
c      
      end if
c
      end
c=============================================================
      subroutine TransformHess(nq,nek,h,ibcontr,b,proj,hint)
      implicit real*8 (a-h,o-z)
      include "common/big"
      parameter(zero=0.0d0,one=1.0d0)
c  This routine takes the Cartesian Hessian matrix h, transforms it
c  to internal coordinates, and takes the projected generalized inverse
c  of  it, which is returned in hint. This is used if there is a Hessian
c  available
c  Arguments
c  INTENT(IN)
c  nq           = number of internal coordinates
c  nek          = 3* the number of atoms
c  h(nek,nek)   = cartesian Hessian
c  ibcontr(20,nq)=the encoding of the B matrix
c  b(54,nq)     = contains the non-zero elements of the B matrix
c  proj         = projector to the non-null subspace of G (unit matrix
c                 if there are no redundancies)
c  INTENT(OUT)
c  hint(nq,nq)  = the INVERSE of the internal Hessian, projected
      dimension h(nek,nek),ibcontr(20,nq),b(54,nq),ginv(nq,nq),
     &          proj(nq,nq),hint(nq,nq)
      ibase=icorr(0)
c  Reserve memory for the full (sparse) B matrix (a little wasteful)
      membb=icorr(nq*nek) 
c  Hes_bmat calculates its own generalized inverse of G - a little waste
      call hes_bmat(h,hint,q(membb),b,ibcontr,nq,nek)
c hes_bmat returns the Hessian transformed to internal coordinates without the
c gradient contributions
c Now take the generalized inverse of hint and project it
      memsq1=icorr(nq**2)
      memsq2=icorr(nq**2)
      memu=icorr(nq**2)
      memx=icorr(nq)
      memy=icorr(nq)
c  Eliminate very small eigenvalues of the Hessian (?)
      tol=1.0d-3
      call geninv(nq,hint,det,tol,q(memu), 
     & q(memx),q(memy),ninv,q(memsq2))
c  q(memsq2) contains a projector but it is not important. It should
c  be the unit matrix
c
      call dgemm('N','N',nq,nq,nq,   one,proj,nq,hint,nq,
     &            zero,q(memsq2),nq)
      call dgemm('N','N',nq,nq,nq,   one,q(memsq2),nq,proj,nq,
     &            zero,hess,nq)
c  hess should contain the correctly projected Hessian, built from the
c  initial Cartesian Hessian 
       call corlsr(ibase)
c
      end
c=============================================================
      subroutine bfgs(nq,qnew,qold,fnew,fold,  hint)
      implicit real*8 (a-h,o-z)
      include "common/big"   
      dimension qnew(nq),qold(nq),fnew(nq),fold(nq),hint(nq,nq)
      parameter(zero=0.0d0,one=1.0d0,small=0.02d0)
c  This routine makes a BFGS update step
c  Arguments
c  INTENT(IN)
c  nq        = nuber of internal coordinates
c  qnew      = current (composite) internal coordinates
c  qold      = old (composite) internal coordinates
c  fnew      = current internal forces
c  fold      = old internal forces
c
c  INTENT(INOUT)
c  hint      = current INVERSE Hessian
c
c  Formula: H(k)=(I-dq*dfT/dfT*dq)H(k-1)(I-df*dqT/dfT*dq) - dq*dqT/dfT*dq
c  Note that dfT*dq is a scalar while the dyadics dq*dfT and df*dqT are
c  nq x nq  matrices
cpp
c      call prntmatn(6,nq,nq,nq,hint,'Hessian in BFGS')
c      write(6,*) 'fold in BFGS'
c      write(6,'(5f12.6)') (fold(k),k=1,nq)
c      write(6,*) 'fnew in BFGS'
c      write(6,'(5f12.6)') (fnew(k),k=1,nq)
c      write(6,*) 'qold in BFGS'
c      write(6,'(5f12.6)') (qold(k),k=1,nq)
c      write(6,*) 'qnew in BFGS'
c      write(6,'(5f12.6)') (qnew(k),k=1,nq)
      ibase=icorr(0)
c  Transfer to qold and fold the differences
      call dscal_x(nq,-one,qold,1)
      call daxpy_x(nq,one,qnew,1,qold,1)
      call dscal_x(nq,-one,fold,1)
      call daxpy_x(nq,one,fnew,1,fold,1)
c      write(6,*) 'qdif in BFGS'
c      write(6,'(5f12.6)') (qold(k),k=1,nq)
c      write(6,*) 'fdif in BFGS'
c      write(6,'(5f12.6)') (fold(k),k=1,nq)
c  qold and fold now contain the changes in the coordinates and in the forces
      denom=ddot_x(nq,qold,1,fold,1)
c  This should not be very small and should be positive for optimization
      qnorm=ddot_x(nq,qold,1,qold,1)
      fnorm=ddot_x(nq,fold,1,fold,1)
      cangle=denom/sqrt(qnorm*fnorm)
      denom=one/denom
cpp
c      write(6,*) 'cangle,qnorm,fnorm',cangle,qnormn,fnorm
      if(cangle.gt.zero) then
        write(6,*) 'Warning in BFGS: cangle should be negative',cangle
      end if
      if(abs(cangle).lt.small) return
c  Calculate h*df 
      memx=icorr(nq)
      call dgemv_x('N',nq,nq,one,hint,nq,fold,1,zero,q(memx),1)
      call dscal_x(nq,denom,q(memx),1)
      do i=1,nq
        do j=1,nq
          hint(i,j)=hint(i,j)-q(memx+i-1)*qold(j)-qold(i)*q(memx+j-1)
        end do
      end do
      pr=(ddot_x(nq,q(memx),1,fold,1)-one)*denom
      do i=1,nq
        do j=1,nq
          hint(i,j)=hint(i,j)+pr*qold(i)*qold(j)
        end do
      end do
c  Test:  If everything is OK, the new Hint*delta f = -delta q
c      write(6,*) 'fold in BFGS before dgemv'
c      write(6,'(5f12.6)') (fold(k),k=1,nq)
c      call prntmatn(6,nq,nq,nq,hint,'Hint before dgemv')
c      call dgemv_x('N',nq,nq,one,hint,nq,fold,1,zero,q(memx),1)
c      write(6,*) 'Compare the two vectors'
c      write(6,'(5f12.6)') (q(memx+k-1),k=1,nq)
c      write(6,'(5f12.6)') (qold(k),k=1,nq)

      call corlsr(ibase)      
c
      end
c=============================================================
