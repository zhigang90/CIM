      SUBROUTINE SCFMAIN(inp)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  ............................................................
C    Wrapper for SCF module
C    Depending on number of alpha/beta electrons and
C    multiplicity calls Closed or Open-shell SCF
C  ............................................................
C
      Character*20 GUES,cdum
      Character*256 jobname
      PARAMETER (IUnit=1)          ! unit for checkpoint read
c
      Common /job/jobname,lenJ
C
C
C  read from <control> file
C  multiplicity
C  number of alpha electrons
C  number of beta electrons
C  SCF guess
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      call rdcntrl(IUnit,6,'$guess',3,idum,dum,GUES)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c -- set <s2> to zero in depository
      call setrval('<s2>',0.0d0)
c
c -- (re)initialize nnmr (# calls to NMR) for new wavefunction
      call setival('nnmr',0)
c
c -- check there are electrons!
      If(NAlpha.EQ.0) Call nerror(1,'SCF Module',
     $                     'There Are No Electrons!',0,0)
c
      If(NBeta.EQ.0.AND.IMult.EQ.1) Then
        Call RHFMAIN(inp,GUES)
      Else
        Call UHFMAIN(inp,GUES)
      EndIF
c
      RETURN
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fixdiag(diag,nmo,xlvsh)
      implicit real*8 (a-h,o-z)
      character*(*) diag
      do i=1,nmo
        call matelem(diag,i,i,xx)
        call mateset(diag,i,i,xx+xlvsh)
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ReadMOS(ncf,    cmo,    eorb,   rdeorb, len,
     $                   file,   itype,  IErr)
      implicit real*8(A-H,O-Z)
C
C  Attempts to read MOs from binary MO file
C  Also - optionally - reads orbital energies (Or their equivalent)
C  return IErr = 0 if successful
C
C  ARGUMENTS
C
C  ncf     -  number of basis functions
C             (actual number of MOs available read from file)
C  cmo     -  storage for MO coefficients
C  eorb    -  orbital energies (or fock matrix expectation values)
C  rdeorb  -  logical flag for reading of orbital energies
C  len     -  number of characters in <file>
C  file    -  filename
C  itype   -  wavefunction type found on MOs file
C              0 - do not check wavefunction type
C              1 - ab initio
C              2 - semiempirical
C              3 - density matrix (full)
C             if input value is greater than 0, then wavefunction type
C             will be checked and MOs will NOT be read if type is wrong
C  IErr    -  error flag on exit
C              0 - success
C             -1 - no MOs read
C
      Dimension cmo(ncf,*),eorb(*)
      Character*(*) file
      Character*35 Char
      Logical rdeorb
C
      IErr = -1
C
C  attempt to open MOS file
C
      OPEN (UNIT=40,FILE=file(1:len),FORM='UNFORMATTED',STATUS='OLD',
     $      ERR=95)
c
      READ(40,ERR=96) nmo,ichk
      write(6,*) "nmo read from mos file",nmo
      If(itype.gt.0.and.ichk.ne.itype) go to 94
      If(nmo.gt.ncf) Call nerror(13,'File IO routine <ReadMOS>',
     $     'Too many MOs on file for current system',0,0)
      CALL ReadBinary(40,ncf,nmo,cmo)
      If(rdeorb) CALL ReadBinary(40,nmo,1,eorb)
c
      itype = ichk
      IErr = 0
 94   CONTINUE
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
 95   CONTINUE
      RETURN
c
 96   CONTINUE
      Char = 'Error Reading Binary MOs file <'//file(1:len)//'>'
      Call nerror(13,'File IO routine <ReadMOS>',Char,0,0)
C
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ReadBinary(IUnit,n1,n2,A)
      real*8 A(n1,n2)
c
      READ(IUnit) A
c
      RETURN
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine WriteMOS(ncf,    nmo,    cmo,    eorb,   wreorb,
     $                    len,    file,   itype)
      implicit real*8(A-H,O-Z)
C
C  Writes orbital coefficients to binary <MOS> file
C  Also - optional - write orbital energies (or their equivalent)
C
C  ARGUMENTS
C
C  ncf     -  number of basis functions
C  nmo     -  number of orbitals to write
C  cmo     -  storage for MO coefficients
C  eorb    -  orbital energies (or fock matrix expectation values)
C  wreorb  -  logical flag for writing of orbital energies
C  len     -  number of characters in <file>
C  file    -  filename
C  itype   -  wavefunction type
C               1 - ab initio
C               2 - semiempirical
C
      Dimension cmo(ncf,nmo),eorb(nmo)
      Character*(*) file
      Logical wreorb
C
C  attempt to open MOS file
C
      OPEN (UNIT=40,FILE=file(1:len),FORM='UNFORMATTED',
     $      STATUS='UNKNOWN',ERR=95)
c
      write(6,*) 'nmo write to mos file',nmo
      WRITE(40,ERR=96) nmo,itype
      WRITE(40,ERR=96) cmo
      If(wreorb) WRITE(40,ERR=97) eorb
c
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
 95   CONTINUE
      Call message('**WARNING from routine <WriteMOS>',
     $ '  Unable to open MO coefficients file, unit =',40,0)
      RETURN
c
 96   Call message('**WARNING from routine <WriteMOS>',
     $ '  Unable to write to MOS file, unit =',40,0)
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
 97   Call message('**WARNING from routine <WriteMOS>',
     $ '  Unable to write orbital energies to MOS file, unit =',40,0)
      RETURN
C
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ChkThrsh(xlow,thint,thresh,iout)
      IMPLICIT REAL*8(A-H,O-Z)
c
c  Checks integral thresholds are compatible with lowest eigenvalue
c  of overlap matrix. Requested thresholds may be reduced.
c
c  ARGUMENTS
c
c  xlow    -  lowest eigienvalue of overlap matrix
c  thint   -  on exit contains rough integral threshold
c  thresh  -  on exit contains final integral threshold
c  iout    -  unit for printout
c
c -- get thresholds from depository
      call getrval('ithr',thresh)
      call getrval('ith1',firsthresh)  ! rough threshold
c
      if(thresh.lt.firsthresh) then
        thint=firsthresh
      else
        thint=thresh
      end if
c
c -- ensure integral thresholds compatible with lowest eigenvalue
c -- of overlap matrix (returned by <sminhalf> as xlow)
c -- if xlow is less than 1.0d-8, quit
cc      If(xlow.lt.1.0d-8) Call nerror(10,'SCF module',
cc     $      'Basis is Severely Linearly-Dependent - Change Basis',
cc     $       0,0)
c
      tratio = thint/xlow
      If(tratio.GT.1.0d-3) Then
       thint = xlow*0.001d0
      Endif
c
      ifact=0
 10   tratio = thresh/xlow
      If(tratio.GT.1.0d-4) Then
       thresh= thresh*0.1d0
       ifact=ifact+1
       go to 10
      EndIf
c
      if(ifact.gt.0) Write(iout,1000) ifact
c
      write(iout,1100) thresh,thint
c
      RETURN
c
 1000 Format('**WARNING** Near Linear Dependency in Basis',/,
     $       '  Sharpening final integral threshold by ',I1,
     $       ' orders of magnitude')
 1100 Format(' Integral thresholds are: Final=',d12.4,' Initial=',d12.4)
c
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine scfresu(NAlpha, NBeta,  ncf,    ncs,    na,
     $                   etot,   e1,     e2,     enuc,   iout,
     $                   ivirt,  inx,    ismpr,  ibas,   iprscf,
     $                   dip)

      use memory

      implicit real*8 (a-h,o-z)
c
c  Calculate and print final SCF results
c  **WARNING** This routine destroys Alpha density matrix if UHF
c
c  ARGUMENTS
c
c  NAlpha  -  number of occupied alpha/closed-shell MOs
c  NBeta   -  number of occupied beta/closed-shell MOs
c  ncf     -  number of contracted basis functions
c  ncs     -  number of contracted shells
c  na      -  number of atoms
c  etot    -  total SCF energy
c  e1      -  one-electron energy
c  e2      -  two-electron energy
c  enuc    -  nuclear repulsion energy
c  iout    -  unit for print out
c  ivirt   -  maximum number of virtual MOS to print
c  inx     -  basis set contraction data
c             (note NOT the same as in common /intbl/)
c  ismpr   -  symmetry equivalent basis function data
c  ibas    -  address of basis set definition block
c  iprscf  -  print flag
c  dip     -  on exit contains dipole moment
c
c  common bl used via <inton> to store 1-electron matrix elements
c
      dimension inx(*),ismpr(7,ncf)
      dimension dip(3),quad(6),octu(10),dipc(3),quadc(6),octuc(10)
      data np1/1/
c
c -- add alpha & beta density matrices
      If(NBeta.GT.0) call matadd('densB','dens')
c
c  define (temporary) dipole matrix
      call matdef('dipp','s',ncf,ncf)
      call matsub('tmatrix','dipp',1,ncf)
c  use dipp temporarily for the kinetic matrix
      call getival('inuc',inuc)
      call inton(2,na,bl(mataddr('tmatrix')),inx,0,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','tmatrix',ekin)
      call matrem('tmatrix')
c  3=dipole
      call inton(3,na,bl(mataddr('dipp')),inx,1,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(1))
      call matwrite('dipp',np1,0,'aoX     ')
c
      call inton(3,na,bl(mataddr('dipp')),inx,2,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(2))
      call matwrite('dipp',np1,0,'aoY     ')
c
      call inton(3,na,bl(mataddr('dipp')),inx,3,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(3))
      call matwrite('dipp',np1,0,'aoZ     ')
c
c  Calculate the components of the electronic quadrupole moment
c  X**2, Y**2, Z**2, XY, XZ, YZ
      call inton(5, na,       bl(mataddr('dipp')), inx, 1,
     1           1, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp',quad(1))
      call inton(5, na,       bl(mataddr('dipp')), inx, 2,
     1           2, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp',quad(2))
      call inton(5, na,       bl(mataddr('dipp')), inx, 3,
     1           3, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp',quad(3))
      call inton(5, na,       bl(mataddr('dipp')), inx, 1,
     1           2, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp',quad(4))
      call inton(5, na,       bl(mataddr('dipp')), inx, 1,
     1           3, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp',quad(5))
      call inton(5, na,       bl(mataddr('dipp')), inx, 2,
     1           3, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp',quad(6))
c
c  Calculate the electronic octupole moments XXX,XXY,XXZ,XYY,XYZ,
c  XZZ,YYY,YYZ,YZZ,ZZZ
      ioctu=0
      do k=1,3
        k0=10*k
        do l=k,3
          k1=k0+l
          do k2=l,3
            ioctu=ioctu+1
      call inton(10, na,       bl(mataddr('dipp')), inx, k1,
     1           k2, bl(ibas), bl(inuc),            ncs)
      call matprodtr('dens','dipp', octu(ioctu))
CPP
c      print *,'Electronic octupole in SCF',ioctu,octu(ioctu)
          end do
        end do
      end do
      call matrem('dipp')
c
c  the electronic part is negative because of the negative charge
      dip(1) = -dip(1)
      dip(2) = -dip(2)
      dip(3) = -dip(3)
c
      do i=1,6
        quad(i)=-quad(i)
      end do
c
      do i=1,10
        octu(i)=-octu(i)
      end do
CPP
c     print *,'Electronic octupole in SCF',octu
c  calculate nuclear dipole and add to dip; similarly the quadrupole
      call nucdipole(na,dip,bl(inuc))
      call nucquadru(na,bl(inuc),quad)
      call nucoctu(na,bl(inuc),octu)
c  Transform to an origin defined by the center of the nuclear charge
      call CenterOfCharge(na,bl(inuc),dip,quad,octu,
     1                    dipc,quadc,octuc)
c     print *,'Total octupole in SCF',octu
c
      call printscf(etot, e1,  e2,  ekin, enuc,
     1              dip,  quad,octu,dipc, quadc,
     2              octuc,iout,iout)
c
      mdiag=mataddr('diag')
      If(NBeta.gt.0) mdiab=mataddr('diagB')
c
      If(iprscf.ge.3) Then
        call druwf(bl(mdiag),'coef',etot,0,ivirt,NAlpha,iout)
        if(NBeta.gt.0) then
          call druwf(bl(mdiab),'coefB',etot,0,ivirt,NBeta,iout)
        endif
      EndIf
c
      call matdef('smat','s',ncf,ncf)
      call matread('smat',np1,'s matrix')
      ism=mataddr('smat')
      icof=mataddr('coef')
      nsym=igetival('nsym')
      Write(iout,1000)
 1000 FORMAT(/,'  Alpha spin')
      call PrintOrbSym(nsym,ncf,NAlpha,ivirt,bl(icof),bl(ism),ismpr,
     1                .true.,bl(mdiag))
      call absmax(NAlpha*ncf,bl(icof),ii,cfmax)
      if(cfmax.gt.5.0d0) then
        write(iout,1050) cfmax
 1050  format('Warning: largest alpha SCF coefficient exceeds 5.0',
     1        f12.5)
      end if
      If(NBeta.gt.0) Then
        icof=mataddr('coefB')
        Write(iout,1100)
 1100   FORMAT(/,'  Beta spin')
        call PrintOrbSym(nsym,ncf,NBeta,ivirt,bl(icof),bl(ism),ismpr,
     1                .true.,bl(mdiab))
      call absmax(NBeta*ncf,bl(icof),ii,cfmax)
      if(cfmax.gt.5.0d0) then
        write(iout,1150) cfmax
 1150  format('Warning: largest beta SCF coefficient exceeds 5.0',
     1   f12.5)
      end if
      EndIf
      call matrem('smat')
c
      end
c=======================================================================
      subroutine nucdipole(na,dip,xnuc)
      implicit real*8 (a-h,o-z)
      dimension dip(3),xnuc(5,*)
c  Transforms the dipole moment to the centroid of the nuclear charge
c  Changes only if the molecule has a net charge
c  ARGUMENTS
c  INTENT(IN)
c  na=number of nuclei
c  xnuc=nuclear data
c  INTENT(INOUT)
c  dip=on input: electronic dipole,
c     on output: dipole vector with nuclear contributions
      do i=1,na
c  nuclear charge
        z=xnuc(1,i)
        do k=1,3
          dip(k)=dip(k)+z*xnuc(k+1,i)
        end do
      end do
      end
c=======================================================================
      subroutine nucquadru(na,xnuc,quad)
      implicit real*8 (a-h,o-z)
      dimension quad(6),xnuc(5,*)
c  Transforms the quadrupole moment to the centroid of the nuclear charge
c  Changes only if the molecule has a net charge
c  ARGUMENTS
c  INTENT(IN)
c  na=number of nuclei
c  xnuc=nuclear data
c  INTENT(INOUT)
c  dip=on input: electronic quadrupole
c     on output: quadrupole with nuclear contributions
      do i=1,na
c  nuclear charge
        z=xnuc(1,i)
        do k=1,3
          quad(k)=quad(k)+z*xnuc(k+1,i)**2
        end do
        quad(4)=quad(4)+z*xnuc(2,i)*xnuc(3,i)   ! XY
        quad(5)=quad(5)+z*xnuc(2,i)*xnuc(4,i)   ! XZ
        quad(6)=quad(6)+z*xnuc(3,i)*xnuc(4,i)   ! YZ
      end do
CPP
c     print *, 'Total quad',quad
      end
c=======================================================================
      subroutine nucoctu(na,xnuc,octu)
      implicit real*8 (a-h,o-z)
c see comments in nucquadru
      dimension octu(10),xnuc(5,*)
      do i=1,na
c  nuclear charge
        z=xnuc(1,i)
        klm=0
        do k=1,3
          do l=k,3
            do m=l,3
              klm=klm+1
              octu(klm)=octu(klm)+z*xnuc(k+1,i)*xnuc(l+1,i)*xnuc(m+1,i)
CPP
c       print *, 'Nucoctu',i,klm,k,l,m,z,xnuc(k+1,i),xnuc(l+1,i),
c     1         xnuc(m+1,i),octu(klm)
C
            end do
          end do
        end do
      end do
      end
c=======================================================================
      subroutine CenterOfCharge(na,   xnuc, dip, quad, octu,
     1                          dic, quac,octc)
      implicit real*8 (a-h,o-z)
      dimension dip(3),quad(6),octu(10),xnuc(5,*),xn(3)
      dimension dic(3),quac(6),octc(10)
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0)
      charge=rgetrval('charge')
c  Calculate the centroid of nuclear charge
      chnuc=zero
      call zeroit(xn,3)
      do i=1,na
         q=xnuc(1,i)
        xx=xnuc(2,i)
        yy=xnuc(3,i)
        zz=xnuc(4,i)
        chnuc=chnuc+q
        xn(1)=xn(1)+q*xx
        xn(2)=xn(2)+q*yy
        xn(3)=xn(3)+q*zz
      end do
      do k=1,3
        xn(k)=xn(k)/chnuc
        dic(k)=dip(k)-xn(k)*charge
      end do
      quac(1)=quad(1)-two*xn(1)*dip(1)+xn(1)**2*charge                 ! xx
      quac(2)=quad(2)-two*xn(2)*dip(2)+xn(2)**2*charge                 ! yy
      quac(3)=quad(3)-two*xn(3)*dip(3)+xn(3)**2*charge                 ! zz
      quac(4)=quad(4)-xn(1)*dip(2)-xn(2)*dip(1)+xn(1)*xn(2)*charge     ! xy
      quac(5)=quad(5)-xn(1)*dip(3)-xn(3)*dip(1)+xn(1)*xn(3)*charge     ! xz
      quac(6)=quad(6)-xn(2)*dip(3)-xn(3)*dip(2)+xn(2)*xn(3)*charge     ! yz
      octc(1)=octu(1)-three*xn(1)*quac(1)-three*xn(1)**2*dic(1)
     1        -xn(1)**3*charge                                         ! xxx
      octc(2)=octu(2)-two*xn(1)*quac(4)-xn(2)*quac(1)-xn(1)**2*dic(2)
     1        -two*xn(1)*xn(2)*dic(1)-xn(1)**2*xn(2)*charge            ! xxy
      octc(3)=octu(3)-two*xn(1)*quac(5)-xn(3)*quac(1)-xn(1)**2*dic(3)
     1        -two*xn(1)*xn(3)*dic(1)-xn(1)**2*xn(3)*charge            ! xxz
      octc(4)=octu(4)-two*xn(2)*quac(4)-xn(1)*quac(2)-xn(2)**2*dic(1)
     1        -two*xn(2)*xn(1)*dic(2)-xn(1)*xn(2)**2*charge            ! xyy
      octc(5)=octu(5)-xn(1)*quac(6)-xn(2)*quac(5)-xn(3)*quac(4)
     1       -xn(1)*xn(2)*dic(3)-xn(1)*xn(3)*dic(2)-xn(2)*xn(3)*dic(1)
     2       -xn(1)*xn(2)*xn(3)*charge                                 ! xyz
      octc(6)=octu(6)-two*xn(3)*quac(5)-xn(1)*quac(3)-xn(3)**2*dic(1)
     1        -two*xn(3)*xn(1)*dic(3)-xn(1)*xn(3)**2*charge            ! xzz
      octc(7)=octu(7)-three*xn(2)*quac(2)-three*xn(2)**2*dic(2)
     1        -xn(2)**3*charge                                         ! yyy
      octc(8)=octu(8)-two*xn(2)*quac(6)-xn(3)*quac(2)-xn(2)**2*dic(3)
     1        -two*xn(2)*xn(3)*dic(2)-xn(2)**2*xn(3)*charge            ! yyz
      octc(9)=octu(9)-two*xn(3)*quac(6)-xn(2)*quac(3)-xn(3)**2*dic(2)
     1        -two*xn(3)*xn(2)*dic(2)-xn(2)*xn(3)**2*charge            ! yzz
      octc(10)=octu(10)-three*xn(3)*quac(3)-three*xn(3)**2*dic(3)
     1          -xn(3)**3*charge                                       ! zzz
c
      end
c=======================================================================
      subroutine printscf(etot,  e1,  e2,  ekin,  enuc,
     1                    dip,  quad, octu,dipc, quadc,
     2                    octuc,iout,icond)
      implicit real*8 (a-h,o-z)
      real*8 kcal,kelvin,kJ
c  prints final SCF results:
c  INPUT:
c etot,e1,e2,ekin,enuc=total,one-el.,two-el.,kinetic and nuclear energy
c  dip(3)=3 components of the dipole moment exp. value
c  quad(6)=6 components of the quadrupole exp. value (xx,yy,zz,xy,xz,yz)
c  iout,icond=printing files
c
      parameter (three=3.0d0, fifteen=15.0d0,xnine=9.0d0)
c
      dimension dip(3),quad(6),octu(10),deby(3),quadau(6),octuau(10)
      dimension dipc(3),quadc(6),octuc(10),octux(10)
C
      virial=-etot/ekin
      call getrval('one',one)
      call getrval('ajou',ajoule)
      call getrval('evol',evolt)
      call getrval('kcal/mol',kcal)
      call getrval('kJ/mol',kJ)
      call getrval('dkel',kelvin)
      call getrval('cmm1',cmm1)
      call getrval('Hertz',Hertz)
      call getrval('deby',debye)
      call getrval('angs',angs)
      debyeangs=debye*angs
      debyeangs2=debye*angs**2
      Eajou=etot/ajoule
      Eevolt=etot/evolt
      Ekj=etot/kJ
      Ekcal=etot/kcal
      Ekelvin=etot/kelvin
      Ecmm1=etot/cmm1
      Ehertz=etot/hertz
      call tfer(dip,deby,3)
      call dscal(3,one/debye,deby,1)
c  transfer the original (atomic unit) values to quadau
      call tfer(quad,quadau,6)
      call dscal(6,one/debyeangs,quad,1)
c  transfer the original (atomic unit) values to octuau
      call tfer(octu,octuau,6)
      call dscal(10,one/(debyeangs*angs),octu,1)
      diptot=sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
      dipdeb=diptot/debye
c  Build multipole moments: 3x**2-r**2,.., 3xy,... for quadrupole
c    15*x**3-9*r**2*x,..., 15x**2*y-3*r**2*y,...  15*x*y*z for octupole
c  Original order of octupoles: xxx, xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
      
      write(iout,100 )etot,Enuc,E1,E2,Ekin,virial,dip,diptot,deby,
     *               dipdeb,quad,octu,eajou,EkJ,Ekcal,Eevolt,Ekelvin,
     *                  Ehertz
      write(iout,*) 
     1  'Traceless electronic quadrupole tensor in au, 3xx-rr,3yy-rr,
     23zz-rr,  3xy,3xz,3yz'
      r2=quadau(1)+quadau(2)+quadau(3)
      quadau(1)=three*quadau(1)-r2 
      quadau(2)=three*quadau(2)-r2 
      quadau(3)=three*quadau(3)-r2 
      quadau(4)=three*quadau(4)
      quadau(5)=three*quadau(5)
      quadau(6)=three*quadau(6)
      write(iout,105) quadau
  105 format((6f12.6,2x))
      write(iout,*) 'Traceless electronic octupole tensor in au'
      write(iout,*)
     1' 15xxx-9rrx,  15yyy-9rry,  15zzz-rrz,    15xxy-3rry,  15xyy-3rrx'
      write(iout,*)
     1  ' 15xxz-3rrz,  15xzz-3rrx,  15yyz-3rrz,  15yzz-3rry,   15xyz'
      r2x=octuau(1)+octuau(4)+octuau(6)
      r2y=octuau(2)+octuau(7)+octuau(9)
      r2z=octuau(3)+octuau(8)+octuau(10)
      octux(1)=fifteen*octuau(1)-xnine*r2x 
      octux(7)=fifteen*octuau(7)-xnine*r2x 
      octux(10)=fifteen*octuau(10)-xnine*r2y 
      octux(2)=fifteen*octuau(2)-three*r2y
      octux(3)=fifteen*octuau(3)-three*r2z
      octux(4)=fifteen*octuau(4)-three*r2x
      octux(6)=fifteen*octuau(6)-three*r2x
      octux(8)=fifteen*octuau(8)-three*r2z
      octux(9)=fifteen*octuau(9)-three*r2y
      octux(5)=fifteen*octuau(5)
      octuau(1)=octux(1)   ! 15xxx-9rrx
      octuau(2)=octux(7)   ! 15yyy-9rry
      octuau(3)=octux(10)  ! 15zzz-9rrz
      octuau(4)=octux(2)   ! 15xxy-3rry
      octuau(5)=octux(4)   ! 15xyy-3rrx
      octuau(6)=octux(3)   ! 15xxz-3rrz
      octuau(7)=octux(6)   ! 15xzz-3rrx
      octuau(8)=octux(8)   ! 15yyz-3rrz
      octuau(9)=octux(9)   ! 15yzz-3rry
      octuau(10)=octux(5)  ! 15xyz

      write(iout,108) octuau
 108  format(3f14.6,2x,2f14.6/,4f14.6,2x,f14.6) 
      write(iout,*) ' '
c
cc      write(icond,100)etot,Enuc,E1,E2,Ekin,virial,dip,deby,quad,
cc     *               eajou,EkJ,Ekcal,Eevolt,Ekelvin,Ehertz
 100  format(/
     * '======================= SCF RESULTS ========================='/
     */' Total Energy =',F25.9,' Eh'/
     * ' nuclear ener.=',f25.9/
     * ' one-el.ener. =',f25.9/
     * ' two-el.ener. =',f25.9/
     * ' kinetic ener.=',f25.9/
     * ' virial coeff.=',f25.6/
     */' dipole/au=',2x,3f12.6, '  total=',f12.6,' au'/
     * ' dipole/D =',2x,3f12.6, '  total=',f12.6,' D'/
     * ' quadrupole/',3x,'  XX=',f13.4,'   YY=',f13.4,'   ZZ=',f13.4,/,
     * ' (Debye*Ang)',3x,'  XY=',f13.4,'   XZ=',f13.4,'   YZ=',f13.4,//,
     * ' octupole/  ',3x,' XXX=',f13.4,'  XXY=',f13.4,'  XXZ=',f13.4,/,
     * ' (Debye*A**2)',2x,' XYY=',f13.4,'  XYZ=',f13.4,'  XZZ=',f13.4,/,
     * '             ',2x,' YYY=',f13.4,'  YYZ=',f13.4,'  YZZ=',f13.4,/,
     * '             ',2x,' ZZZ=',f13.4,/,
     * '-------------------------------------------------------------'/
     * '    Total Energy in different units :'/
     *    /f16.6,'aJ',f16.3,'kJ/mol',f16.3,'kcal/mol',/
     *     f16.5, 'eV', f16.2, 'K ',e16.9,'Hz ',//
     * '======================= SCF RESULTS =========================')
c
      call tfer(dipc,deby,3)
      call dscal(3,one/debye,deby,1)
      call dscal(6,one/debyeangs,quadc,1)
      call dscal(10,one/(debyeangs*angs),octuc,1)
      write(iout,*)
     1 'Electric moments relative to the centroid of nuclear charge ',
     2 '(if different)'
      call absmaxdif(3,dip,dipc,i,dmax)
      if(dmax.gt.1.27d-7) then
        write(iout,110) deby
      end if
      call absmaxdif(6,quad,quadc,i,dmax)
      if(dmax.gt.0.5d-4) then
        write(iout,120) quadc
      end if
      call absmaxdif(10,octu,octuc,i,dmax)
      if(dmax.gt.0.5d-3) then
        write(iout,130) octuc
      end if
  110 format(
     * ' dipole/D =',2x,3f12.6)
  120 format(
     * ' quadrupole/',3x,'  XX=',f13.4,'   YY=',f13.4,'   ZZ=',f13.4,/,
     * ' (Debye*Ang)',3x,'  XY=',f13.4,'   XZ=',f13.4,'   YZ=',f13.4,/)
  130 format(
     * ' octupole/  ',3x,' XXX=',f13.4,'  XXY=',f13.4,'  XXZ=',f13.4,/,
     * ' (Debye*A**2)',2x,' XYY=',f13.4,'  XYZ=',f13.4,'  XZZ=',f13.4,/,
     * '             ',2x,' YYY=',f13.4,'  YYZ=',f13.4,'  YZZ=',f13.4,/,
     * '             ',2x,' ZZZ=',f13.4)
c
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine geneig(fock,u,uinv,diag,coef,dens,xlvsh,nroots,upper)

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) fock,u,uinv,diag,coef,dens,upper
      character*4 temp
      integer*4 info
c  this routine solves the generalized eigenvalue
c  equation FC = SCe   The matrices U and Uinv should satisfy
c  Univ(t)SUinv = I(identity) or S=U(t)U where (t) is transpose;
c  Uinv is the inverse of U.
c  It solves the equation (Uinv(t)FUinv(UC)=(UC)e
c  F =fock, u=U, uinv=Uinv, coef=C, diag=e
c  input data: fock,u,uinv,dens,xlvsh dens is needed for the level shift
c               nroots: number of roots desired. If nroots=0, all
c               roots are calculated, otherwise the lowest nroots
c              'upper' is a character variable, either 'U' or 'L'
c              it is 'u' if the matrix 'u' is an upper triangle,
c              and it is 'L' if it is a lower triangle. this has been
c              added so that it can handle the reverse case which
c              is needed somewhere in the NMR program
c  output: coef,diag
c  note that xlvsh must be half of the intended level shift
c  because it is multiplied by the density at double occupancy
c
c introduced triangular matrix multiplies for speed
c using lapack generalized eigensolver routines would be even faster
c just the levelshift should be applied in AO basis
c
c----------test------------
c      write(6,*)' from geneig :'
c      call matprint(fock, 6  )
c      call matprint(u   , 6  )
c      call matprint(uinv, 6  )
c----------test------------
      call getival('ncf',ncf)
      zero=rgetrval('zero')
      one=rgetrval('one')
      temp='temp'
      call matdef(temp,'s',ncf,ncf)
      iadra=mataddr(temp)
      iadrw=mataddr(diag)
      iadrc=mataddr(coef)
      iadru=mataddr(u)
      iadrui=mataddr(uinv)
c calculate Uinv(t)FUinv taking advantage of Uinv being triangular
c  use 'coef' as a temporary matrix
      call matcopy(fock,coef)
      call dtrmm('L',upper,'T','N',ncf,ncf,one,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
      call dtrmm('R',upper,'N','N',ncf,ncf,one,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
      call matcopy(coef,temp)
c     call matsimtr(fock,uinv,temp)
c
      if(xlvsh.ne.zero) then
c calculate U*D*U(t)  taking advantage of U being triangular
c  use 'coef' as a temporary matrix
        call matcopy(dens,coef)
        call dtrmm('L',upper,'N','N',ncf,ncf,one,bl(iadru),ncf,
     *             bl(iadrc),ncf)
        call matscal(coef,-xlvsh)
        call dtrmm('R',upper,'T','N',ncf,ncf,one,bl(iadru),ncf,
     *             bl(iadrc),ncf)
        call matadd(coef,temp)
c       call matmmult(u,dens,coef)
c       call matscal(coef,-xlvsh)
c       call matmmul2(coef,u,temp,'n','t','a')
c  the above operation subtracts U*D*U(t) from Uinv(t)*F*Uinv
      end if
c
c sdiag2 is a lot slower - use always LAPACK
c     if(nroots.eq.0) then
c       call matdiag(temp,diag,coef)
c     else
      if(nroots.eq.0) then        ! get all eigenvalues
         nroo=ncf
      else                        ! get first nroots eigenvalues
         nroo=nroots
      end if
      call getmem(8*ncf,iadrtmp)
      call getmem(5*ncf,iwork)
      call getmem(ncf,iadrfail)
      call dspevx('V','I','U',ncf,bl(iadra),zero,zero,1,nroo,
     1            1.0d-10,m,bl(iadrw),bl(iadrc),ncf,bl(iadrtmp),
     2            bl(iwork),bl(iadrfail),info)
      call retmem(3)
c calculate Uinv*Coeff taking advantage of Uinv being triangular
      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
c     call matmmult(uinv,coef,coef)
      call matrem(temp)
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine geneig2(dens,  u,     uinv,  diag,  coef,
     1                   nroots)

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) u,uinv,diag,coef,dens
      character*4 temp
      integer*4 info
c     common /big/bl(1000)
c  this routine solves the generalized eigenvalue
c  equation SDSC = SCe   The matrices U and Uinv should satisfy
c  Univ(t)SUinv = I(identity) or S=U(t)U where (t) is transpose;
c  Uinv is the inverse of U.
c  It solves the equation UDU(t)(UC)=(UC)e
c  D =density, u=U, uinv=Uinv, coef=C, diag=e
c  The sign of the density is reversed so that the eigenvectors come in 
c  descending order
c ARGUMENTS
c INTENT(IN)
c   dens = name of the charge density matrix
c   u    = name of the U matrix, usually triangular (will not work with general matrices)
c   uinv = name of the U inverse, both upper triangular
c   nroots= usually 0, meaning all eigenvalues are calculated. Can be set to
c           nocc+ a larger number (but less than ncf)
c   INTENT(OUT)
c   diag = occupation numbers of the NOs (name of the diagonal matrix)
c   coef = orbital coefficeints of the NOs (name of the matrix)
c
c  because it is multiplied by the density at double occupancy
c
c introduced triangular matrix multiplies for speed
c using lapack generalized eigensolver routines would be even faster
c
c----------test------------
c       write(6,*)' from geneig :'
c      call matprint(fock, 6  )
c      call matprint(u   , 6  )
c      call matprint(uinv, 6  )
c----------test------------
      call getival('ncf',ncf)
      zero=rgetrval('zero')
      one=rgetrval('one')
      xlvsh=one
      temp='temp'
      call matdef(temp,'s',ncf,ncf)
      iadra=mataddr(temp)
      iadrw=mataddr(diag)
      iadrc=mataddr(coef)
      iadru=mataddr(u)
      iadrui=mataddr(uinv)
c
c calculate U*D*U(t)  taking advantage of U being triangular
c  use 'coef' as a temporary matrix
        call matcopy(dens,coef)
        call dtrmm('L','U','N','N',ncf,ncf,one,bl(iadru),ncf,
     *             bl(iadrc),ncf)
        call matscal(coef,-xlvsh)
        call dtrmm('R','U','T','N',ncf,ncf,one,bl(iadru),ncf,
     *             bl(iadrc),ncf)
        call matcopy(coef,temp)
c
c sdiag2 is a lot slower - use always LAPACK
c     if(nroots.eq.0) then
c       call matdiag(temp,diag,coef)
c     else
      if(nroots.eq.0) then        ! get all eigenvalues
         nroo=ncf
      else                        ! get first nroots eigenvalues
         nroo=nroots
      end if
      call getmem(8*ncf,iadrtmp)
      call getmem(5*ncf,iwork)
      call getmem(ncf,iadrfail)
      call dspevx('V','I','U',ncf,bl(iadra),zero,zero,1,nroo,
     1            1.0d-10,m,bl(iadrw),bl(iadrc),ncf,bl(iadrtmp),
     2            bl(iwork),bl(iadrfail),info)
        call retmem(3)
c calculate Uinv*Coeff taking advantage of Uinv being triangular
      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
c     call matmmult(uinv,coef,coef)
      call matrem(temp)
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine densma(coef,dens,nmo,rhf)
      implicit real*8 (a-h,o-z)
      character*(*) coef,dens
      logical rhf
      parameter(two=2.0d0)
c
c -- if no occupied orbitals (for one-electron atom)
c -- simply zero density matrix
      If(nmo.eq.0) Then
       call matzero(dens)
       return
      EndIf
c
      call matsear('occu',iexist)
      if(iexist.eq.0) then
        call matsub('occu',coef,1,nmo)
      end if
      call matmmul2('occu','occu',dens,'n','t','n')
      if(rhf) call matscal(dens,two)
      call matrem('occu')
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine sminhalf(s,u,uinv,xlow)
c  this routine calculates the square root and inverse square root of
c  the overlap matrix s
c  INPUT:
c  s: the name of the (symmetric) matrix
c  the name of u=s**(1/2)
c  uinv the name of s**(-1/2)
c
c  OUTPUT:
c  upon output, u and uinv contain S**(1/2) and S**(-1/2)
c  xlow is the lowest eigenvalue of S

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) s,u,uinv
      character*8 eval,evec
c     common /big/bl(1000)
      parameter (one=1.0d0)
      call matinfo(s,mshape,idim1,idim2,iaddr,ilen)
      if(mshape.ne.3) then
        call nerror(1,'sminhalf','Overlap matrix is not symmetric',
     1   mshape,0)
      end if
      eval='XYZ87344'
      evec='YZX87344'
      call matdef(eval,'d',idim1,idim1)
      call matdef(evec,'q',idim1,idim1)
      call matdiag(s,eval,evec)
      call matpose(evec)
      iaddr=mataddr(eval)
      xlow=bl(iaddr)
      iaddr=iaddr-1
      do i=1,idim1
        bl(iaddr+i)=sqrt(bl(iaddr+i))
      end do
      call matsimtr(eval,evec,u)
      do i=1,idim1
        bl(iaddr+i)=one/bl(iaddr+i)
      end do
      call matsimtr(eval,evec,uinv)
      call matrem(evec)
      call matrem(eval)
      call getival('icond',icond)
      call getival('iout',iout)
      write(iout,100) xlow
cc      write(icond,100) xlow
 100  format(' Lowest eigenvalue of the overlap matrix ',e14.6)
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine cholesky(matrix,upper,n,xlow)

      use memory

      implicit real*8(a-h,o-z)
c This subroutine performs Cholesky factorization on matrix and returns
c     the upper triangular factor upper.  Matrix must be symmetric and
c     positive definite.
C  Arguments:
c     matrix: name of the  S matrix (character var.)
c     upper: name of the upper Cholesky factor (character var.)
c     matrix must be packed upper triangle styored by columns
c     order (11) (12) (22) (13) (23) (33) (14)...  or, what is the same,
c           (11) (21) (22) (31) (32) (33) ....
c     upper is a quadratic matrix currently because sminhalf yields a
c      quadratic one
c     n is the dimension
c     xlow is the smalles diagonal element of the factor
      character*(*) matrix, upper
c     common /big/bl(1000)
      parameter (zero=0.0d0)
      integer*4 i4nfo
c First create upper if necessary and get the addresses of the matrices.
      call matsear(upper,jcur)
      if(jcur.eq.0) then
         call matdef(upper,'q',n,n)
         iupperaddr=mataddr(upper)
      else
         call matinfo(upper,mshape,idim1,idim2,iupperaddr,ilen)
         if(mshape.ne.2) then
            call nerror(1,'Cholesky','Target matrix is not square
     1',mshape,0)
         end if
         if(ilen.ne.n*n) then
            call nerror(2,'Cholesky','Target matrix is not the right
     1      size',ilen, n*n)
         end if
      end if
      imatrixaddr=mataddr(matrix)
c     create a temporary matrix to hold intermediate results
      call matdef('temp','v',n*(n+1)/2,1)
      itempaddr=mataddr('temp')
c
c   calculate the lowest eigenvalue of the matrix by diagonalization
c    This is not bad, as the vectors are not needed
      call getmem(8*n,itempx)
      call getmem(5*n,itempi)
      call getmem(n,ifail)
      call getmem(n,ieval)
      call tfer(bl(imatrixaddr),bl(itempaddr),n*(n+1)/2)
      call dspevx('N','I','U',n,bl(itempaddr),VLdum,VUdum,1,1,
     1 1.0d-10,mfound,bl(ieval),Zdummy,n,bl(itempx),bl(itempi),
     1 bl(ifail), i4nfo)
      info=int(i4nfo)
c   xlow is the lowest eigenvalue of the S matrix
      xlow=bl(ieval)
      call retmem(4)
c
c     Transfer matrix into temp (again because dspevx destroys it)
      call tfer(bl(imatrixaddr),bl(itempaddr),n*(n+1)/2)
c     do the Cholesky factorization
      call dpptrf('U',n,bl(itempaddr),i4nfo)
      info=int(i4nfo)
c
      call getival('iout',iout)
      write(iout,100) xlow
 100  format(' Lowest eigenvalue of the overlap matrix',e14.6)
c     error checking
      if(info.lt.0) then
         call nerror(3,'Cholesky','Invalid argument in
     1 lapack routine dpptrf : number ',info,0)
      end if
      if(info.gt.0) then
         call nerror(4,'Cholesky','Matrix is not
     1positive definite',info,0)
      end if
c     transfer the result to the matrix upper
      call packed_to_upper(bl(itempaddr),bl(iupperaddr),n)
      call matrem('temp')
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine u_inverse(upper,upperinv,n)
c     Calculate the inverse of the upper triangular matrix

      use memory

      character*(*) upper, upperinv
      integer*4 i4nfo
c     real*8 bl
c     common /big/bl(1000)

c Find out about memory for upperinv

      call matsear(upperinv,jcur)
      if(jcur.eq.0) then
         call matdef(upperinv,'q',n,n)
         iupperinv_addr=mataddr(upperinv)
      else
         call matinfo(upperinv,msh,idim1,idim2,iupperinv_addr,ilen)
         if(msh.ne.2) then
            call nerror(1,'U_inverse','Target matrix is not square
     1',msh,0)
         end if
         if(ilen.ne.n*n) then
            call nerror(2,'U_inverse','Target matrix is not the right
     1 size',ilen, n*n)
         end if
      end if
      iupper_addr=mataddr(upper)

c make temporary matrix

      call matdef('temp','v',n*(n+1)/2,1)
      itemp_addr=mataddr('temp')

c transfer upper to packed format

      call upper_to_packed(bl(iupper_addr),bl(itemp_addr),n)

c find the inverse

      call dtptri('U','N',n,bl(itemp_addr),i4nfo)
      info=int(i4nfo)
      if(info.ne.0) then
        call nerror(1, 'U_inverse','S matrix is singular',
     1   info,info)
      end if
c transfer the result to upperinv

      call packed_to_upper(bl(itemp_addr),bl(iupperinv_addr),n)
      call matrem('temp')

      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine packed_to_upper(source,target,n)
c     transfer the packed upper triangular matrix source to the target
c     matrix in unpacked form
      real*8 source(n*(n+1)/2),target(n,n)

c first set the lower triangle of target to zero

      do j=1,n
         do i=j,n
            target(i,j)=0.0d0
         end do
      end do

c transfer source into the upper triangular portion of target

      jc=1
      do j=1,n
         do i=1,j
            target(i,j)=source(jc+i-1)
         end do
         jc=jc+j
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine upper_to_packed(source,target,n)
c     copies the upper triangular source matrix into packed format in the
c     target linear array.
      real*8 source(n,n),target(n*(n+1)/2)

      jc=1
      do j=1,n
         do i=1,j
            target(jc+i-1)=source(i,j)
         end do
         jc=jc+j
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine converged(u,uinv,density,fock,commut,error,errsq)

      use memory

      implicit real*8 (a-h,o-z)
c  ARGUMENTS:
c  INTENT(IN)
c  U and Uinv are decompositions of the overlap matrix S such that
c     U(t)U = S Uinv is the inverse of U; U may be S**(1/2),
c     though here it is the Cholesky decomposition of S
c  density is the name of the density matrix D(max. 8 char.)
c  fock is the name of the fock matrix F
c  commut is the commutator; it MUST be defined as antisymmetric,
c  otherwise the results will be wrong
c  INTENT(OUT)
c  error is the norm of the commutator UDFUinv - (UDFUinv)transpose
c
      character*(*) u,uinv,density,fock,commut
      character*4 imed
c     common /big/bl(10000)
c  Calculate UDFUinv - its transpose
      call getival('ncf ',ncf)
      ncf2=ncf*(ncf+1)/2
c
      iu=mataddr(u)
      iuinv=mataddr(uinv)
      icommut=mataddr(commut)
      imed='imed'
      call matdef(imed,'q',ncf,ncf)
      iimed=mataddr(imed)
c
c this is only needed for the intel blas library
c that does not conform to the standard and does not zero out the result
c problems surfaced in the MPI code that used uninitialized memory
      call matzero(imed)
c
c use triangular matrix multiply with U and Uinv
      call matmmult(density,fock,imed)
      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iu),ncf,
     *           bl(iimed),ncf)
      call dtrmm('R','U','N','N',ncf,ncf,1.0d0,bl(iuinv),ncf,
     *           bl(iimed),ncf)
c
      call matcopy(imed,commut)
      call matrem(imed)
      call absmax(ncf2,bl(icommut),iii,error)
      errsq=ddot(ncf2,bl(icommut),1,bl(icommut),1)
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine druwf (ew,coef,etot,lv,ndru,nel,ifi)

      use memory

      implicit double precision (a-h,o-z)
c  this interface gets some basic data for the printig routine druwff
c  INPUT:
c  ew(*) = orbital energies
c  coef  = name of the matrix of the wave function coefficients
c  etot  = total energy
c  lv    = number of virtuals to punch
c  ndru  = number of virtuals to print
c  nel   = number of occupied orbitals
c  ifi   = file to print on
c
      dimension ew(*)
      character*(*) coef
c     common /big/bl(10000)
c     common /intbl/ maxsh,inx(100)
      call getival('ncf',ncf)
      call getival('na',na)
      call getival('ncs',ncs)
      call getival('inuc',inuc)
      call getival('ictr',ictr)
      call getrval('zero',zero)
      call druwff(ncf, ew, etot, bl(mataddr(coef)),ndru,
     1            nel, ifi,na,   ncs,  inuc, bl(ictr))
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine WriteOrbEn(label,nocc,ntot,diag)
c  This routine is needed in Hyperchem to print the orbital energies
c  Arguments:
c   INTENT(IN)
c   label = the type of the matrix to be printed, 4 characters count,
c           can be rhf , alph, beta, loca, natu
c   nocc  = number of occupied orbitals to be printed
c   ntot  = number of total orbitals to be printed
c   diag  = name of the matrix holding the orbital energies
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) label,diag
      character*256 jobname
      parameter (ntypes=5)
      character*4 types(ntypes),type
      character*8 diag1
c     common /big/bl(10000)
      data types/'dble','alph','beta','loca','natu'/
c  needed for Hyperchem
c      return
c
c  determine the type of file to be printed
      do i=1,ntypes
        if(label(1:4).eq.types(i)) then
          type=types(i)
        end if
      end do
c  get the address of the matrix holding the orbital energies
      diag1=diag
      iew=mataddr(diag1)-1
c
c  get the jobname and open the file
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
      open(1,file=jobname(1:lenJ)//'.'//type(1:3))
c  write the number of occupied & virtual orbitals
      write(1,'(2i10)') nocc,ntot
c  write orbital energies
      write(1,'((f14.6))') (bl(iew+i),i=1,ntot)
c  close file
      close(1)
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine druwff(ncf,  ew,  etot,  wf,  ndru,
     1                  nel,  ifi, na,    ncs, nadr,
     2                  inx)
c .........
c     ncf   = number of ao basis functions
c     ew    = eigenvalues or similar quantities to be printed
c        for each orbital
c     etot  = total energy
c     wf is the array containing the lcao coefficients in a matrix
c       the columns of which are the mo coefficients
c     ndru and nel: typically, nel is the number of occupied
c       orbitals to be printed, and ndru is the number of virtuals,
c       so the first nel+ndru columns of wf (and the corresponding
c       energies or occupation numbers) are printed
c     ifi   = is the file to print on
c     na    = number of atoms
c     ncs   = number of contracted shells
c     nadr  = start of the nuclear info in bl (for atom symbols)
c     inx   = basis set info array
c..........

      use memory

      implicit double precision (a-h,o-z)
      character*4 hold(6)
c .................................................
c -- dynamic allocation of character arrays
      character symb(ncf)*4,atnum(ncf)*3,ftype(ncf)*16
c .................................................
      dimension ew(*), wf(ncf,*),inx(12,*)
c     common /big/bl(1000)
      norb=nel+ndru
c--------------------------- punch wavefunction -----------------------c
      if (norb.gt.ncf) norb=ncf
c---------- prepare vertical array of atomic orbital labels -----------c
      do i=1,ncf
        atnum(i)='   '
        symb(i)='    '
      end do
      ndx=0
      do 400 i=1,ncs
         ngc=inx(4,i)
      do 400 igc=1,ngc+1
         nt=ndx+1
         n=inx(2,i)
         iftyp=inx(12,i)
         write(atnum(nt),'(i3)') n
         write(symb(nt),'(a4)') bl(nadr+4+(5*(n-1)))
         if(iftyp.eq.1) then
           ftype(ndx+1)='s'
           ndx= ndx+1
           goto 400
         endif
         if(iftyp.eq.2) then
           ftype(ndx+1)='x'
           ftype(ndx+2)='y'
           ftype(ndx+3)='z'
           ndx= ndx+3
           goto 400
         endif
         if(iftyp.eq.3) then
           ftype(ndx+1)='s'
           ftype(ndx+2)='x'
           ftype(ndx+3)='y'
           ftype(ndx+4)='z'
           ndx= ndx+4
           goto 400
         endif
         if(iftyp.eq.4) then
           ftype(ndx+1)='3zz-rr'
           ftype(ndx+2)='xx-yy'
           ftype(ndx+3)='xy'
           ftype(ndx+4)='xz'
           ftype(ndx+5)='yz'
           ndx= ndx+5
           goto 400
         endif
         if(iftyp.eq.5) then
           ftype(ndx+1)='xx'
           ftype(ndx+2)='yy'
           ftype(ndx+3)='zz'
           ftype(ndx+4)='xy'
           ftype(ndx+5)='xz'
           ftype(ndx+6)='yz'
           ndx= ndx+6
           goto 400
         endif
         if(iftyp.eq.6) then
           ftype(ndx+1)='5xxy-rry'
           ftype(ndx+2)='5xxz-rrz'
           ftype(ndx+3)='5yyx-rrx'
           ftype(ndx+4)='5yyz-rrz'
           ftype(ndx+5)='5zzx-rrx'
           ftype(ndx+6)='5zzy-rry'
           ftype(ndx+7)='xyz'
           ndx= ndx+7
           goto 400
         endif
         if(iftyp.eq.7) then
           ftype(ndx+1)='xxx'
           ftype(ndx+2)='xxy'
           ftype(ndx+3)='xxz'
           ftype(ndx+4)='xyy'
           ftype(ndx+5)='xyz'
           ftype(ndx+6)='xzz'
           ftype(ndx+7)='yyy'
           ftype(ndx+8)='yyz'
           ftype(ndx+9)='yzz'
           ftype(ndx+10)='zzz'
           ndx= ndx+10
           goto 400
         end if
         if(iftyp.eq.8) then
           ftype(ndx+1)='xxxx'
           ftype(ndx+2)='xxxy'
           ftype(ndx+3)='xxxz'
           ftype(ndx+4)='xxyy'
           ftype(ndx+5)='xxyz'
           ftype(ndx+6)='xxzz'
           ftype(ndx+7)='xyyy'
           ftype(ndx+8)='xyyz'
           ftype(ndx+9)='xyzz'
           ftype(ndx+10)='xzzz'
           ftype(ndx+11)='yyyy'
           ftype(ndx+12)='yyyz'
           ftype(ndx+13)='yyzz'
           ftype(ndx+14)='yzzz'
           ftype(ndx+15)='zzzz'
           ndx= ndx+15
           goto 400
         endif
         if(iftyp.eq.9) then
           do kk=1,21
             ftype(ndx+kk)='h21'
           end do
           ndx=ndx+21
           go to 400
         end if
         if(iftyp.eq.10) then
           do kk=1,28
             ftype(ndx+kk)='i28'
           end do
           ndx=ndx+28
           go to 400
         end if
         if(iftyp.eq.11) then
           ftype(ndx+1)='x4+y4-6x2y2'
           ftype(ndx+2)='x4+z4-6x2z2'
           ftype(ndx+3)='y4+z4-6y2z2'
           ftype(ndx+4)='x3y-3xyz2'
           ftype(ndx+5)='xy3-3xyz2'
           ftype(ndx+6)='x3z-3xy2z'
           ftype(ndx+7)='xz3-3xy2z'
           ftype(ndx+8)='y3z-3x2yz'
           ftype(ndx+9)='yz3-3x2yz'
           ndx=ndx+9
         end if
         if(iftyp.eq.12) then
           ftype(ndx+1)='x5-10x3y2+5xy4'
           ftype(ndx+2)='x5-10x3z2+5xz4'
           ftype(ndx+3)='y5-10x2y3+5x4y'
           ftype(ndx+4)='y5-10y3z2+5yz4'
           ftype(ndx+5)='z5-10x2z3+5x4z'
           ftype(ndx+6)='z5-10y2z3+5y4z'
           ftype(ndx+7)='x4z+y4z-6z2y2z'
           ftype(ndx+8)='x4y+yz4-6x2yz2'
           ftype(ndx+9)='xy4+yz4-6xy2z2'
           ftype(ndx+10)='x3yz-xy3z'
           ftype(ndx+11)='x3yz-2xyz3+xy3z'
           ndx=ndx+11
         end if
         if(iftyp.eq.13) then
           do kk=1,13
             ftype(ndx+kk)='i'
           end do
           ndx=ndx+13
         end if
  400 continue
c-------------------- initialize printing variables -------------------c
      n=1
      m=3
      icnt=0
      mm=norb-(norb/4)*4 - 1
      do 500 i=1,4
         hold(i)=' occ'
  500 continue
c---------------------------- printing loop ---------------------------c
  600 do 700 i=1,4
         icnt=icnt+1
         if(icnt.gt.nel) hold(i)='virt'
  700 continue
      if(icnt.gt.norb) m=mm
      write(ifi,'(1x)')
      write(ifi,705) (i,i=n,n+m)
  705 format(1x,25x,4(7x,i4))
      write(ifi,710) (ew(i),i=n,n+m)
  710 format(' atom   type             no.',4(f11.6))
      write(ifi,720) (hold(i),i=1,1+m)
  720 format(1x,28x,4(4x,a4,3x))
      do 800 i=1,ncf
         write(ifi,900) atnum(i),symb(i),ftype(i),i,(wf(i,j),j=n,n+m)
  800 continue
      write(ifi,'(1x)')
  900 format(a3,1x,a4,a16,i4,4(f11.6))
      n=icnt+1
      if(icnt.lt.norb) goto 600
c----------------------------------------------------------------------c
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pdiag(fock,coef,d,nmo,ncf,cutoff,xlvsh)
c  this routine caries out a pseudo-diagonalization

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) fock, coef, d
      integer nmo, ncf
      call matsub('occu',coef,1,nmo)
      ioccu=mataddr('occu')
      call matinfo(coef,mshape,id1,id2,maddr,length)
      call matsub('virt',coef,nmo+1,id2)
      ivirt=mataddr('virt')
c
c     calculate FOV=occu(t)*F*virt
c
      nvi=id2-nmo
      call matdef('FOV','r',nmo,nvi)
      call matdef('temp','r',nmo,ncf)
      call matmmul2('occu',fock,'temp','t','n','n')
      call matmmult('temp','virt','FOV')
      call absmax(nmo*nvi,bl(mataddr('FOV')),ii,fovmax)
c     call matprint('FOV',6)
c
c     Prepare to calculate the matrix of rotation angles.
c     Declare matrices and get addresses of needed matrices.
c
      ifm=mataddr('FOV')
      idm1=mataddr(d)
      idm2=idm1+nmo
c
c     Calculate the rotation angles and simultaneously determine the
c     maximum rotation angle.
c
      call getangles(bl(ifm),bl(idm1),bl(idm2),nmo,nvi,ncf,smax,
     1  xlvsh)
c
c     calculate the rotation threshold
c
      tresh=cutoff*smax
c
c     perform the rotations that are above the threshold
c
      call rotate(bl(ioccu),bl(ivirt),bl(ifm),nmo,nvi,ncf,tresh)
      call matmmul2('occu',fock,'temp','t','n','n')
      call matmmult('temp','virt','FOV')
c     print *, 'after rotation'
c     call matprint('FOV',6)
cc      call absmax(nmo*nvi,bl(mataddr('FOV')),ii,fovmax1)
cc      call getival('iout',iout)
cc      write(iout,100) fovmax,fovmax1
  100 format('max. Brillouin violating element before and after',
     1 2e10.4)
      if(abs(xlvsh).lt.1.0d-3) then
        call getangles(bl(ifm),bl(idm1),bl(idm2),nmo,nvi,ncf,smax,
     1    xlvsh)
        tresh=cutoff*smax
        call rotate(bl(ioccu),bl(ivirt),bl(ifm),nmo,nvi,ncf,tresh)
      end if
      call matrem('temp')
      call matrem('FOV')
      call matrem('virt')
      call matrem('occu')
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine getangles(fov,d1,d2,nmo,nvi,ncf,smax,xlvsh)
      implicit real*8 (a-h,o-z)
c    this routine generates Jacobi rotational angles
c    INPUT
c    fov=Brillouin block of the orbital Fock matrix (destroyed)
c    d1= first (occcupied) set of approximate orbital energies
c    d2=second (virtual) approximate orbital energies
c    nmo,nvi= number of the orbitals in the first and second set
c      (occupied and virtual)
c    ncf=rotational angles (more accurately, sinuses)
c     smax= maximum rot. angle
c    xlvsh= level shift
c
c    OUTPUT
c    fov= rotational angles
      real*8 fov(nmo,nvi),d1(nmo),d2(nvi)
      parameter(zero=0.0d0,half=0.5d0)
      smax=zero
      angmax=0.1d0
      do j=1,nvi
         do i=1,nmo
           denom=d2(j)-d1(i)+xlvsh
           if(abs(denom).gt.0.0001d0) then
             s=fov(i,j)/(d2(j)-d1(i)+xlvsh)
           else
             s=sign(fov(i,j),angmax)
           end if
           if(abs(s).gt.angmax) s=sign(s,angmax)
           fov(i,j)=s
           if(abs(s).gt.smax) smax=abs(s)
         end do
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rotate(cocc,cvirt,x,nmo,nvi,ncf,tresh)
      implicit real*8 (a-h,o-z)
c  INPUT
c  cocc(ncf,nmo): first (occupied) set of orbital coefficients
c  cvirt(ncf,nvi): second (virtual) set of SCF coefficients
c  x(nmo,nvi): rotational angles
c  ncf: number of basis functions
c  tresh (note spelling): minimum rot. angle
      real*8 cocc(ncf,nmo),cvirt(ncf,nvi),x(nmo,nvi)
      parameter(one=1.0d0)
      do ii=1,nmo
         do ia=1,nvi
            s=-x(ii,ia)
            c=sqrt(one-s**2)
            if(abs(x(ii,ia)).gt.tresh) then
              call drot(ncf,cocc(1,ii),1,cvirt(1,ia),1,c,s)
            end if
         end do
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine combinedf(ntri,nfock,ifock,dm,dens,df)
      implicit real*8 (a-h,o-z)
c  this routine puts the ifock-th density matrix in the
c  combined matrix df, and zeroes out the second half of df if
c  ifock=1
c  INPUT:
c  ntri=lenghth of the Fock matrix, ncf*(ncf+1)/2
c  nfock2=number of fock + density matrices
c  ifock: the number corresponding to dm
c  df: large area  containing both the density & the Fock
      dimension dm(ntri),dens(nfock,ntri),df(nfock,ntri)
      data zero/0.0d0/

      do i=1,ntri
        dens(ifock,i)=dm(i)
      end do
      if(ifock.eq.1) then
        do i=1,ntri
          do j=1,nfock
            df(j,i)=zero
          end do
        end do
      end if
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine extractdf(ntri,nfock,ifock,dm,fm,df,dn)
      implicit real*8 (a-h,o-z)
c  this routine extracts the ifock-th density and fock matrices in the
c  combined matrix df
c  INPUT:
c  ntri=lenghth of the Fock matrix, ncf*(ncf+1)/2
c  nfock=number of fock matrices
c  ifock: the number corresponding to dma and fm
c  fm= fock matrix
c  df: large area  containing the Fock matrices
c  dn: contains the density matrices
      dimension dm(ntri),fm(ntri),df(nfock,ntri),dn(nfock,ntri)

      do i=1,ntri
        dm(i)=dn(ifock,i)
        fm(i)=df(ifock,i)
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GetIntThrsh
      implicit real*8(a-h,o-z)
C
C  if job is an optimization, check rms gradient
C  if below gthrsh, set initial and final integral thresholds
C  the same
C
      character*256 jobname
      logical isthere
      parameter (gtol=5.0d-3)
C
C  is there an <opt> file?
C
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
      inquire(file=jobname(1:lenJ)//'.opt',exist=isthere)
c
      IF(isthere) THEN
C
C  check for existence of 'rmsg'
C
        call tstrval('rmsg',ifound)
        If(ifound.eq.1) Then
          call getrval('rmsg',rmsg)
          if(rmsg.lt.gtol) then
            call getrval('ithr',thresh)
            call setrval('ith1',thresh)
          endif
        EndIf
      ENDIF
c
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE uhfNO(ncf,    NAlpha, NBeta,  dens,   densB,
     $                 olddA,  u,      uinv,   diag,   coef,
     $                 lprint, sexpv,  xmulti)

      use memory

      implicit real*8(a-h,o-z)
c
c  calculate UHF natural orbitals
c
c  ARGUMENTS
c
c  ncf     -  number of basis functions
c  NAlpha  -  number of occupied alpha MOs
c  NBeta   -  number of occupied Beta MOs
c  dens    -  alpha density matrix
c  densB   -  beta density matrix
c  olddA   -  work storage (for density sum)
c  u       -  S**1/2
c  uinv    -  S**-1/2
c  diag    -  storage for eigenvalues
c  coef    -  storage for UHF natural orbitals
c  lprint  -  print level; if > 1 then the fractionally
c               occupied NOs are printed
c
c  on exit
c
c  sexpv   -  expectation value <S**2>
c  xmulti  -  calculated spin multiplicity
c
c ** WARNING - This subroutine uses TEXAS matrix system **
      character*(*) dens,densB,olddA,u,uinv,diag,coef
c     common /big/bl(30000)
      parameter (Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0)
      parameter (np4=4)
c
c  Use the negative of the density matrix so that the largest
c  eigenvalues come out first
      call matcopy(dens,olddA)
      call matadd(densB,olddA)
      call geneig2(olddA,u,uinv,diag,coef,0)
c  restore the sign of diag
      call matscal(diag,-One)
c  write out these quantities
      call matwrite(coef,np4,0,'nato_uhf')
      call matwrite(diag,np4,0,'naoc_uhf')
c    idiag is the starting address of diag
      idiag=mataddr(diag)-1
      eps=0.005d0
      eps1=Two-eps
      nlow=0
      nhi=0
c  expectation value of <S**2>=[(na-nb)/2]**2 + n - Sum[de(i)**2]/2
c  The double of this quantity is calculated here first
      sexpv = Half*(NAlpha-Nbeta)**2 + Two*(NAlpha+Nbeta)
      do i=1,ncf
        if(bl(idiag+i).lt.eps1.and.nlow.eq.0) nlow=i
        if(bl(idiag+i).lt.eps.and.nhi.eq.0)   nhi=i-1
        sexpv=sexpv-bl(idiag+i)**2
      end do
c  The (fractional) multiplicity is sqrt(4*<S**2>+1); note that
c   at this point sexpv is twice its final value
      xmulti=sqrt(Two*sexpv+One)
      sexpv=sexpv*Half
c
      call setrval('<s2>',sexpv)     ! store in depository
c
      iout = igetival('iout')
      write(iout,*)
      write(iout,650) sexpv,xmulti
      If(nlow.gt.0.and.nlow.le.nhi) Then
        write(iout,700) eps,eps1
        do i=nlow,nhi
        write(iout,750) i,bl(idiag+i)
        enddo
        if(lprint.gt.0) then
          call matsub('natorb',coef,nlow,nhi)
          call druwf(bl(idiag+nlow),'natorb',0.0d0,0,0,nhi-nlow+1,
     1     iout)
          call matrem('natorb')
        end if
      EndIf
c
      return
c
 650  format(' Expectation value of S**2 and multiplicity=',2f11.7)
 700  format(' natural occupation numbers between ',f5.3,' and ',f5.3)
 750  format(1x,i4,2x,f10.6)
c
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine DensAnneal(coef,dens,nmo,rhf,temp)
      implicit real*8 (a-h,o-z)
c This routine constructs a finite temperature density matrix, based on the
c molecular orbital energies and the Fermi-Dirac statistics
c The latter gives the occupancy as n(i)=
      character*(*) coef,dens
      logical rhf
      parameter(two=2.0d0)
c
c -- if no occupied orbitals (for one-electron atom)
c -- simply zero density matrix
      If(nmo.eq.0) Then
       call matzero(dens)
       return
      EndIf
c

      call matmmul2('occu','occu',dens,'n','t','n')
      if(rhf) call matscal(dens,two)
      call matrem('occu')
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine FractOcc(coef,diag,temp,nmo,rhf,dens,Fermi)
c  This routine calculates a density matrix based on a finite temperature
c  occupancy and the Dirac-Fermi statistics.
c  Arguments
c  INTENT(IN)
c  coef      = the orbital coefficient matrix
c  diag      = the orbital energies
c  temp      = kT, the average thermal energy, in Hartrees
c  nmo       = number of electrons (for a given spin in UHF) or electron pairs in RHF
c  rhf       = logical, true if douby occupied orbitals
c  intent(OUT)
c  dens      = the density matrix
c  Fermi     = the Fermi level

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) coef,dens,diag
      logical rhf
      parameter(two=2.0d0)
c     common /big/bl(10000)
c
c -- if no occupied orbitals (e.g. for a one-electron atom)
c -- simply zero density matrix
      If(nmo.eq.0) Then
       call matzero(dens)
       return
      EndIf
c
c  Recall the number of contracted orbitals
      ncf=igetival('ncf')
c  Define the occupation number matrix
      call matdef('occnum','d',ncf,ncf)
c  Determine the address of the orbital energy matrix and the occupation no. matrix
      idiag=mataddr(diag)
      ioccnum=mataddr('occnum')
c  Determine the Fermi level and the occupation numbers
      call FermiLev(ncf,nmo,bl(idiag),temp,bl(ioccnum),Fermi)
cPP
c      print *, 'Fermi',Fermi
      call matdef('Cnmatrix','q',ncf,ncf)
      call matmmult(coef,'occnum','Cnmatrix')
      call matmmul2('Cnmatrix',coef,dens,'n','t','n')
c  Dens contains CnC(t) where n= the diagonal matrix of occupation numbers
c  Remove temporary matrices
      call matrem('Cnmatrix')
cPP
c      call matprint(diag,6)
c         call matprint('occnum',6)
      call matrem('occnum')
c  Multiply the density by 2 if RHF
      if(rhf) call matscal(dens,two)
       end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine fixdiagF(diag,n,xlvsh,Fermi,temp)
c  This routiner restores the orbital energies to their proper values in
c  a level shift calculation with fractional occupancy.
c  The level shift is implemented by adding to the Fock matrix, transformed
c  to orthonormal basis by Uinv(t)*F*Uinv, the quantity -lvsh*U*D*U(t)
c  Here F is the Fock matrix, U is the upper triangle of the Cholesky decomposition
c  of the overlap matris S: S=U(t)*U, U is upper triangular, D is the density
c  matrix (with unit occupancies; for a closed shell density, use lvsh/2),
c  and lvsh is the level shift,
c  This manipulation lowers occupied energy levels by lvsh*n(i), where n(i) is the
c  occupation number. If no fractional occupancy is used then simply add lvsh
c  to all occupied eigenvalues, as in fixdiag
c  To restore the values in the fractional occupation case, we have to add
c  lvsh*n(i) where ni=1/(1+exp((ei-F)/T), and ni= occupation number, ei=orbital
c  energy, F=Fermi energy, T=temperature
c  arguments
c  INTENT(INOUT)
c  diag    = symbolic NAME of the eigenvalue matrix
c  INTENT(IN)
c  n       = number of orbitals requiring correction (in general, all)
c  xlvsh   = level shift
c  Fermi   = Fermi level
c  temp    = temperature
      implicit real*8 (a-h,o-z)
      character*(*) diag
      parameter(one=1.0d0)
cPP
c      print *, 'Level Shift=',xlvsh
      do i=1,n
        call matelem(diag,i,i,epsi)
        if(epsi.lt.Fermi) then
           occnu=one/(one+exp((epsi-Fermi)/temp))
        else
c  Use alternative expression to avoid overflow if e>Fermi
           expo=exp((Fermi-epsi)/temp)
           occnu=expo/(expo+one)
        end if
cPP
c        print *, 'i,epsi,occnu',i,epsi,occnu
        call mateset(diag,i,i,epsi+xlvsh*occnu)
      end do
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine FermiLev(ncf,nmo,epsi,temp,occnum,Fermi)
      implicit real*8 (a-h,o-z)
c  This routine forms a set of occupation numbers for a finite temperature
c  Arguments
c  INTENT(IN)
c  ncf    = number of contracted basis functions
c  nmo    = number of electrons (for UHF) or electron pairs (for RHF)
c  epsi   = orbital energies, assumed to form an acending sequence
c  temp   = finite temperature (in Hartrees)
c  IMTENT(OUT)
c  occnum = Fermi-Dirac occupation numbers, given as ni=1/{1+exp[(ei-z)/kT]}
c           where z is the Fermi level. The latter has to be determined so that the
c           number of electrons comes out exactly
c  Fermi  = the Fermi level
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,tol=1.0d-10)
      dimension epsi(ncf),occnum(ncf)
c  Get a first approximation to z as (HOMO+LUMO)/2
      z=(epsi(nmo)+epsi(nmo+1))*half
      ncyc=0
  100 continue
      ncyc=ncyc+1
      totel=zero
      deriv=zero
      do i=1,ncf
        ei=epsi(i)
        if(ei.lt.z) then
           expo=exp((ei-z)/temp)
           occi=one/(one+expo)
           deri=occi**2*expo
         else
c  Use an alternative formula, exp[(z-ei)/kT]/{exp[(z-ei)/kT+1}
           expo=exp((z-ei)/temp)
           occi=expo/(expo+one)
           deri=occi/(expo+one)
         end if
         totel=totel+occi
         deriv=deriv+deri
         occnum(i)=occi
      end do
c  deriv is the derivative of the sum of the occupation numbers w.r.t. the Fermi
c  energy z
      deriv=-deriv/temp
      devi=totel-dble(nmo)
c      print *,'ncyc,totel',ncyc,totel
      if(abs(devi).gt.tol.and.ncyc.lt.100) then
        dz=devi/deriv
        z=z+dz
        go to 100
      end if
      if(ncyc.ge.25) then
        call nerror(1,'FermiLev',
     1  'Determination of the Fermi level does not converge',nmo,ncyc)
      end if
      Fermi=z
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine savebas4mp2(bl,datbas,inx,ncs,ncf,nsh,nbf)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension datbas(13,nsh)
      dimension inx(12,ncs)

c small basis set data (datbas & inx) will be written on a disk

c copy basis set info into name_sm :

      call setival('ncs_sm',ncs)
      call setival('ncf_sm',ncf)
      call setival('nsh_sm',nsh)
      call setival('nbf_sm',nbf)

c allocate double precision matrix for inx .

      call matdef('inx_sm','r',12,ncs)

      inxsm=mataddr('inx_sm')

      call set_inx_sm(inx,ncs,bl(inxsm))

      np4=4
      call wri(datbas,13*nsh,np4,0,'smallbas')
      call matwrite('inx_sm',np4,0,'smallinx')

      call matrem('inx_sm')

      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine set_inx_sm(inx,ncs,sm_inx)
      implicit real*8 (a-h,o-z)
      dimension inx(12,ncs),sm_inx(12,ncs)

      do ics=1,ncs
         do i=1,12
            sm_inx(i,ics)=dble( inx(i,ics) )
         enddo
      enddo

      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FormMaxDen(NBas,DA,DM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum density matrix element per column.
C  Used for density threshold in DFT code
C  (NOTE: Used for density AND difference-density matrices)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  DA      -  density matrix (lower triangle)
C  DM      -  on exit contains maximum density element per column
C
C
      DIMENSION DA(NBas*(NBas+1)/2),DM(NBas)
C
      DO 20 I=1,NBas
      DMx = 0.0d0
      II = (I*(I-1))/2
      DO 10 J=1,I
      IJ = II+J
      If(Abs(DA(IJ)).GT.DMx) DMx = Abs(DA(IJ))
 10   CONTINUE
      DO 11 J=I+1,NBas
      IJ = (J*(J-1))/2 + I
      If(Abs(DA(IJ)).GT.DMx) DMx = Abs(DA(IJ))
 11   CONTINUE
      DM(I) = DMx
 20   CONTINUE
C
      RETURN
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FormMaxDenU(NBas,DA,DB,DM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Forms vector of maximum density matrix element per column.
C  Used for density threshold in DFT code
C  (NOTE: Used for density AND difference-density matrices)
C
C  ARGUMENTS
C
C  NBas    -  number of basis functions
C  DA      -  alpha density matrix (lower triangle)
C  DB      -  beta density matrix (lower triangle)
C  DM      -  on exit contains maximum density element per column
C
C
      DIMENSION DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NBas)
C
      DO 20 I=1,NBas
      DMx = 0.0d0
      II = (I*(I-1))/2
      DO 10 J=1,I
      IJ = II+J
      If(Abs(DA(IJ))+Abs(DB(IJ)).GT.DMx) DMx = Abs(DA(IJ))+Abs(DB(IJ))
 10   CONTINUE
      DO 11 J=I+1,NBas
      IJ = (J*(J-1))/2 + I
      If(Abs(DA(IJ))+Abs(DB(IJ)).GT.DMx) DMx = Abs(DA(IJ))+Abs(DB(IJ))
 11   CONTINUE
      DM(I) = DMx
 20   CONTINUE
C
      RETURN
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SumSCF(etot,iter)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Writes the SCF energy, number of SCF iterations and basis set
C  to the summary LOG file
C
C  ARGUMENTS
C
C  etot    -  total SCF energy
C  iter    -  number of SCF iterations
C
C
      Character jobname*256,cdum*20
      Parameter (IUnit=1)
      Common /job/jobname,lenJ
c
      call getival('icond',icon)
C
C  open the <control file>
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
c
c -- recover basis set name
      call rdcntrl(IUnit,6,'$basis',3,idum,rdum,cdum)
      call RmBlank(20,cdum,len)
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      WRITE(icon,1000) etot,iter,cdum(1:len)
C
      RETURN
c
 1000 Format(/,' SCF Energy: ',f16.9,'   iterations: ',i3,
     $         '   basis: ',a16)
c
      END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine delete_ghosts(natom,XNUC,latom,XC,IAN)
      implicit real*8(a-h,o-z)
C
C  removes ghost atoms from the set of Cartesian coordinates
C  (originally written for DFT dispersion)
C
C  ARGUMENTS
C
C  natom   -  number of real atoms
C  XNUC    -  Cartesian coordinates, TEXAS format
C
C  on exit
C
C  latom   -  number of real atoms with ghost atoms removed
C  XC      -  Cartesian coordinates of non-ghost atoms
C  IAN     -  atomic numbers of non-ghost atoms
C
      REAL*8 XNUC(5,natom),XC(3,natom)
      INTEGER IAN(natom)
C
      latom = 0
      DO IAtm=1,natom
      If(XNUC(1,IAtm).GT.0.5d0) Then
        latom = latom+1
        IAN(latom) = NINT(XNUC(1,IAtm))
        XC(1,latom) = XNUC(2,IAtm)
        XC(2,latom) = XNUC(3,IAtm)
        XC(3,latom) = XNUC(4,IAtm)
      EndIf
      EndDO
C
      RETURN
      END
