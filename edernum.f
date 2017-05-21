c========================================================================
c
      subroutine edernum(inp,done)
      use memory
      implicit real*8(a-h,o-z)
      character*256 chopval
c
c Calculates the first- and the second drivatives of the total energy 
c with respect to an external electric field
c
c
c It also calculates first & secend derivatives of the gradient
c if the <grad> file exists
c
c This means the dipole moment and dipole polarizability
c and polarizability derivatives
c
c  reads the NUMElectric line in the input file 
c  and writes options (if any) to the <control> file
c
      parameter (nopt=5)
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),ifound(nopt)
      character*4 options(nopt)
      character cdum*20
      character*256 jobname
      character*80 wvfnc
c 
      dimension fiel0(3)   ! EXISTING electric field
      dimension field(3)   ! external electric field
      dimension dipSCF(3)  ! dipole moment from SCF 
      dimension dipole(3)  ! dipole moment from  e.g.MP2 - calculated here
      dimension polari(6)  ! dipole moment from  e.g.MP2 - calculated here
c      xx xy xz yy yz zz
c
      dimension ep(3),em(3)! energy with + and - field
      dimension epp(3),emm(3)! energy with +- field
      dimension epm(3),emp(3)! energy with +- field
c
      logical fflag
      logical lflag
      logical done
      logical diponly
      logical lgrad
c
      parameter (IUnit=1)
      data zero /0.0d0/
      data options/'nume','fdst','file','prin','dipo'/
      data ioptyp/0,11,21,1,0/
      data ncall /0/
      data delta /0.0d0/
      save ncall
      save e0,escf,dipSCF
      save ep,em, epp,epm,emp,emm
      save delta
      save fiel0
      save ntrans
      save diponly
      save ig0,igp,igm, igpp,igmm,igpm,igmp ! gradient in field pointers
      save natom
      save lgrad 
c---------------------------------------------------------------------
c           nacll filed: x       y      z
c           ------------------------------
c             1          0       0      0   !  no electric field
c     ep(1-3)
c             2         +h       0      0
c             3          0      +h      0
c             4          0       0     +h
c     em(1-3)
c             5         -h       0      0
c             6          0      -h      0
c             7          0       0     -h
c     epp(1-3)
c             8         +h      +h      0
c             9         +h       0     +h
c            10          0      +h     +h
c
c     epm(1-3)
c
c            11         +h      -h      0
c            12         +h       0     -h
c            13          0      +h     -h
c
c     emp(1-3)
c
c            14         -h      +h      0
c            15         -h       0     +h
c            16          0      -h     +h
c
c     emm(1-3)
c
c            17         -h      -h      0
c            18         -h       0     -h
c            19          0      -h     -h
c---------------------------------------------------------------------
      call getival('iout',iout)
      call getival('icon',icon)
c---------------------------------------------------------------------
c initializing
c
c     deal with jobname header
c
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
      call izeroit(iopval,3*nopt)
      call  zeroit(ropval,3*nopt)
      call readop1(inp,    nopt,   options,ioptyp, iopval,
     $             ropval, chopval,ifound)
c
c ...........................................................
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
c 
      if(ncall.eq.0) then
c
         call mmark
c
         call fdcntrl(iunit,9,'$wavefunc',iend)
         if(iend.eq.0)call rdcntrl(IUnit,9,'$wavefunc',3,idum,dum,wvfnc)
c        ....................................................
         call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
c
         call WrPOL(IUnit,field,  2  )
c
c        dummy atoms
c
         ndum1=0     ! number of charged dummies
         call fdcntrl(IUnit,7,'$ndummy',idum)
         If(idum.EQ.0) Then
          backspace IUnit
          READ(IUnit,900) Ndum1,Ndum2
          NAtom = NAtoms-Ndum1-Ndum2
         Else
          NAtom = NAtoms
         EndIf
  900 Format(9X,I4,2X,I4)
c
c ...........................................................
c        allocate memory for gradient
c
         natom3=natom*3
         call getmem(  natom3,ig0 )
         call getmem(3*natom3,igp )
         call getmem(3*natom3,igm )
         call getmem(3*natom3,igpp)
         call getmem(3*natom3,igmm)
         call getmem(3*natom3,igpm)
         call getmem(3*natom3,igmp)
c
         call zeroit(bl(ig0 ),  natom3)
         call zeroit(bl(igp ),3*natom3)
         call zeroit(bl(igm ),3*natom3)
         call zeroit(bl(igpp),3*natom3)
         call zeroit(bl(igmm),3*natom3)
         call zeroit(bl(igpm),3*natom3)
         call zeroit(bl(igmp),3*natom3)
c ...........................................................
c
c        get existing external electric field
c
         call getfield0(fiel0)
         write(iout,99) fiel0
 99      format('existing field=',3(f9.4,1x))
c
         call zeroit(field,3)
ccccc    call setfield(field)
c ...........................................................
         call zeroit(ep,3)
         call zeroit(em,3)
c ...........................................................
c
         if(ifound(2).eq.1) then
           delta = ropval(1,2)
         else
           delta = 0.005d0     ! default finite difference step size
         endif
c
         if(ifound(4).eq.1) then
           IPRNT = iopval(1,4)
         else
           IPRNT = 2          ! default print flag
         endif
         write(iout,1000)
         write(iout,1100) delta
c
         diponly=.false.
         if(ifound(5).eq.1) diponly=.true.
c
      endif
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
 1000 FORMAT(
     */' CALCULATING POLARIZABILITY BY CENTRAL DIFFERENCE ON ENERGY')
 1100 FORMAT('          Finite Difference Step Size: ',F9.6,/)
c---------------------------------------------------------------------
      ncall=ncall+1
c---------------------------------------------------------------------
c check if the <grad> file exist :
c
      if(ncall.eq.2) then
         inquire(file=jobname(1:lenJ)//'.grad',exist=lgrad) 
cccc     write(6,*)' gardient file :',lgrad
      endif
c---------------------------------------------------------------------
c set up the field
c
      if(ncall.eq.1) then
         field(1)=zero
         field(2)=zero
         field(3)=zero
      endif
      if(ncall.eq.2) then
         field(1)=delta+delta
         field(2)=zero
         field(3)=zero
      endif
      if(ncall.eq.3) then
         field(1)=zero
         field(2)=delta+delta
         field(3)=zero
      endif
      if(ncall.eq.4) then
         field(1)=zero
         field(2)=zero
         field(3)=delta+delta
      endif
      if(ncall.eq.5) then
         field(1)=-delta-delta
         field(2)=zero
         field(3)=zero
      endif
      if(ncall.eq.6) then
         field(1)=zero
         field(2)=-delta-delta
         field(3)=zero
      endif
      if(ncall.eq.7) then
         field(1)=zero
         field(2)=zero
         field(3)=-delta-delta
      endif
c
c for second derivatives
c             8         +h      +h      0
c             9         +h       0     +h
c            10          0      +h     +h
c
c            11         +h      -h      0
c            12         +h       0     -h
c            13          0      +h     -h
c
c            14         -h      +h      0
c            15         -h       0     +h
c            16          0      -h     +h
c
c            17         -h      -h      0
c            18         -h       0     -h
c            19          0      -h     -h
      if(ncall.eq.8) then
         field(1)=delta
         field(2)=delta
         field(3)=zero
      endif
      if(ncall.eq.9) then
         field(1)=delta
         field(2)= zero 
         field(3)=delta
      endif
      if(ncall.eq.10) then
         field(1)= zero
         field(2)=delta 
         field(3)=delta
      endif
      if(ncall.eq.11) then
         field(1)= delta
         field(2)=-delta 
         field(3)= zero
      endif
      if(ncall.eq.12) then
         field(1)= delta
         field(2)= zero 
         field(3)=-delta
      endif
      if(ncall.eq.13) then
         field(1)= zero 
         field(2)= delta
         field(3)=-delta
      endif
      if(ncall.eq.14) then
         field(1)=-delta
         field(2)= delta
         field(3)= zero 
      endif
      if(ncall.eq.15) then
         field(1)=-delta
         field(2)= zero 
         field(3)= delta
      endif
      if(ncall.eq.16) then
         field(1)= zero
         field(2)=-delta
         field(3)= delta
      endif
      if(ncall.eq.17) then
         field(1)=-delta
         field(2)=-delta
         field(3)= zero
      endif
      if(ncall.eq.18) then
         field(1)=-delta
         field(2)= zero
         field(3)=-delta
      endif
      if(ncall.eq.19) then
         field(1)= zero
         field(2)=-delta
         field(3)=-delta
      endif
c end of the story
      if(ncall.eq.20) then
         field(1)= zero
         field(2)= zero
         field(3)= zero
      endif
c
      call add2vec(3,fiel0,field,field)
c
      call setfield(field)
c
c---------------------------------------------------------------------
      if(ncall.lt.20 .or.(ncall.lt.8.and.diponly)) then
      write(iout,*)'  '
      write(iout,1234)ncall, field
 1234 format(' Step =',i3,' :  SCF with the field=',3(f8.5,2x))
      call flush(iout)
      endif
c---------------------------------------------------------------------
C  Read from the <control> file energy 
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
ctest
c     call RdField(IUnit,Field,FFlag)
c         write(6,*)' field flag=',fflag,' ncall=',ncall
ctest
c.....................................................................
         if(ncall.EQ.2) then
            call rdcntrl(IUnit,7,'$energy',2,idum,E0,cdum)
            call rdcntrl(IUnit,5,'$escf',2,idum,ESCF,cdum)
            call  rddip(iunit,dipSCF,lflag)
            if(lgrad) call getgrad(natom,  1  ,bl(ig0)) ! zero-field gradient
         endif
c.....................................................................
         if(ncall.gt.2.and.ncall.le.5) then
            icall=ncall-2
            call rdcntrl(IUnit,7,'$energy',2,idum,EP(icall),cdum)
            if(lgrad) call getgrad(natom,icall,bl(igp))
         endif
         if(ncall.gt.5.and.ncall.le.8) then
            icall=ncall-5
            call rdcntrl(IUnit,7,'$energy',2,idum,EM(icall),cdum)
            if(lgrad) call getgrad(natom,icall,bl(igm))
         endif
         if(ncall.gt.8.and.ncall.le.11)then
            icall=ncall-8
            call rdcntrl(IUnit,7,'$energy',2,idum,EPP(icall),cdum)
            if(lgrad) call getgrad(natom,icall,bl(igpp))
         endif
         if(ncall.gt.11.and.ncall.le.14)then
            icall=ncall-11
            call rdcntrl(IUnit,7,'$energy',2,idum,EPM(icall),cdum)
            if(lgrad) call getgrad(natom,icall,bl(igpm))
         endif
         if(ncall.gt.14.and.ncall.le.17)then
            icall=ncall-14
            call rdcntrl(IUnit,7,'$energy',2,idum,EMP(icall),cdum)
            if(lgrad) call getgrad(natom,icall,bl(igmp))
         endif
         if(ncall.gt.17.and.ncall.le.20)then
            icall=ncall-17
            call rdcntrl(IUnit,7,'$energy',2,idum,EMM(icall),cdum)
            if(lgrad) call getgrad(natom,icall,bl(igmm))
         endif
c.....................................................................
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
      done=.false.
c
      if(diponly) then
        if(ncall.eq. 8) done=.true.
      else
        if(ncall.eq.20) done=.true.
      endif
c     write(iout,*)' done=',done,'  ncall=',ncall
c---------------------------------------------------------------------
c     call matremark
c     call retmark
c     call retimark
c---------------------------------------------------------------------
      if(done) then
c    
calculate dipolemoment : dipole=(E(+2h)-E(-2h))/(4h)
c
            denom=2.d0*(delta+delta)
            denom=1.d0/denom
            do icall=1,3
               dipole(icall)=-(ep(icall)-em(icall))*denom 
            enddo
c
         if(diponly) go to 4321
c
c calculate polarizability :
c             8         +h      +h      0
c  epp        9         +h       0     +h
c            10          0      +h     +h
c
c            11         +h      -h      0
c epm        12         +h       0     -h
c            13          0      +h     -h
c
c            14         -h      +h      0
c emp        15         -h       0     +h
c            16          0      -h     +h
c
c            17         -h      -h      0
c epp        18         -h       0     -h
c            19          0      -h     -h
c
            denom=(2.0d0*delta)**2
            denom=1.d0/denom
c
            polari(1)= ep(1) -2.d0*e0 + em(1)  !  xx
            polari(4)= ep(2) -2.d0*e0 + em(2)  !  yy
            polari(6)= ep(3) -2.d0*e0 + em(3)  !  zz
c 
            polari(2)=epp(1)-epm(1)-emp(1)+emm(1)  !  xy
            polari(3)=epp(2)-epm(2)-emp(2)+emm(2)  !  xz
            polari(5)=epp(3)-epm(3)-emp(3)+emm(3)  !  yz
c
            call dscal(6,-denom,polari,1)
c
 4321 continue
c
c calculate dipole moment and polarizability derivatives :
c
            if(lgrad) then
               call getmem(natom*3*6,ipold)
               call getmem(natom*3*3,idipd)
               call calcDPD(natom,bl(idipd),bl(ipold),bl(ig0),bl(igp),
     *                     bl(igm),bl(igpp),bl(igpm),bl(igmp),bl(igmm),
     *                     diponly)
               dipd=0.25d0/delta
               call dscal(natom*3*3, dipd ,bl(idipd),1)
               if(.not.diponly) call dscal(natom*3*6, denom,bl(ipold),1)
            endif
ctest
c           if(lgrad) then
c              call getmem(natom*3*6,ipol1)
c              call calcP1(natom,bl(ipol1),bl(ig0),bl(igp),bl(igm),
c    *                     bl(igpp),bl(igpm),bl(igmp),bl(igmm))
c              call dscal(natom*3*6, denom,bl(ipol1),1)
c              call print_pold(iout,NAtom,bl(ipol1))
c              call retmem(1) ! ipol1
c           endif
c           if(lgrad) then
c              call getmem(natom*3*6,ipol2)
c              call calcP2(natom,bl(ipol2),bl(ig0),bl(igp),bl(igm),
c    *                     bl(igpp),bl(igpm),bl(igmp),bl(igmm))
c              call dscal(natom*3*6, denom,bl(ipol2),1)
c              call print_pold(iout,NAtom,bl(ipol2))
c              call retmem(1) ! ipol2
c           endif
ctest
c
c---------------------------------------------------------------------
c write original energy E0, SCF energy and calculated dipole momment
c to the <control> file
c
         OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $         FORM='FORMATTED',STATUS='OLD')
            call wrcntrl(IUnit,7,'$energy',2,idum,E0,cdum)
            call wrcntrl(IUnit,5,'$escf',2,idum,ESCF,cdum)
            call  wrdip(iunit,dipole)
            call fdcntrl(iunit,9,'$wavefunc',iend)
            wvfnc='          '
            if(iend.eq.0) then
                call rdcntrl(IUnit,9,'$wavefunc',3,idum,dum,wvfnc)
            endif
         CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c---------------------------------------------------------------------
         call setfield(fiel0)
         ncall=0
         if(diponly) then
            write(iout,2001)
         else
            write(iout,2000)
         endif
c---------------------------------------------------------------------
c print out energy, dipole moment and poliarizability
c
         call print_dipolar(iout,escf,e0,wvfnc, dipscf, dipole,polari,
     *                      diponly)
         call print_dipolar(icon,escf,e0,wvfnc, dipscf, dipole,polari,
     *                      diponly)
c---------------------------------------------------------------------
c print out poliarizability derivatives
c
         if(lgrad) then
cccc        call print_pold(iout,NAtom,bl(ipold))
cccc        call print_pold(icon,NAtom,bl(ipold))
            call print_dipold(iout,NAtom,bl(idipd),bl(ipold),diponly)
            call print_dipold(icon,NAtom,bl(idipd),bl(ipold),diponly)
            call retmem(1) ! idipd
            call retmem(1) ! ipold
         endif
c---------------------------------------------------------------------
         call retmark
      endif   ! done
c---------------------------------------------------------------------
 2000 FORMAT(/,
     $' -------------------------------------------------------',/,
     $' ** APPARENTLY SUCCESSFUL POLARIZABILITY CALCULATION  **',/,
     $' -------------------------------------------------------',/)
c
 2001 FORMAT(/,
     $' -------------------------------------------------------',/,
     $' ** APPARENTLY SUCCESSFUL DIPOLE MOMENT  CALCULATION  **',/,
     $' -------------------------------------------------------',/)
c
c---------------------------------------------------------------------
C
      END
c========================================================================
      subroutine print_dipolar(iout,escf,e0,wvfnc, dipscf, dipole,
     *                         polari,diponly)
      implicit real*8(a-h,o-z)
      dimension dipscf(3),dipole(3),polari(6)
      character*80 wvfnc
      character*80 scfwv
      logical diponly
c---------------------------------------------------------------------
c print out header and energy :
c
      write(iout,1200)
      write(iout,1150)
 1150 format(
     * '--            DIPOLE MOMENT and POLARIZABILITY',
     * '                        --',/,
     * '--       calculated by finite-difference on energy',
     * '                    --')
      write(iout,1155) wvfnc(1:10)
 1155 format(
     * '--                 The ',a10   ,' method      ',
     * '                        --')
      write(iout,1200)
      write(iout,777)'SCF       ' , escf
      write(iout,777) wvfnc(1:10),e0
  777 format(a10,' Energy  =',f20.9)
c
c print out dipole moment :
c
c total dipole moment
c
      dipTscf=sqrt( dipscf(1)**2+dipscf(2)**2+dipscf(3)**2)
      dipTwvf=sqrt( dipole(1)**2+dipole(2)**2+dipole(3)**2)
c
      write(iout,1200)
      write(iout,1250)
      write(iout,1200)
c
      write(iout,1175)
 1175 format(28x,' X', 7x,' Y', 7x,' Z',10x,' Total ')
      write(iout,88) 'SCF       ',dipSCF,dipTscf
      write(iout,88) wvfnc,       dipole,dipTwvf
   88 format(/,a10,' dipole    =',2x,3f9.4,4x,f9.4,' [au]')
c
      write(iout,1200)
c
 1250 format(
     * '--                        DIPOLE MOMENT ',18x,'            --',/
     * '--                       (atomic',
     * ' units)                               --')
c---------------------------------------------------------------------
      if(diponly) return
c---------------------------------------------------------------------
c Experimentally determined polarizability invariants are
c (P. Calaminici, K.Jug, A.M.Koster, J.C.P., 109,(18),7756 (1998)
c
c
c (1) mean polarizability : alpha=1/3(Axx+Ayy+Azz)
c (2) polarizabilty anisotropy: 1/2[(Axx-Ayy)^2+(Axx-Azz)^2+(Ayy-Azz)^2]
c
         amean=0.3333d0*(polari(1)+polari(4)+polari(6))
         aaniz=      (polari(1)-polari(4))**2
         aaniz=aaniz+(polari(1)-polari(6))**2
         aaniz=aaniz+(polari(4)-polari(6))**2
         aaniz=0.5d0*aaniz
c---------------------------------------------------------------------
c -- print out polarizabilites :
c
      write(iout,1200)
      write(iout,1300)
      write(iout,1200)
      write(iout,1400) polari(1),polari(2),polari(4),
     $                 polari(3),polari(5),polari(6)
      write(iout,1200)
      write(iout,1500) Amean, Aaniz
      write(iout,1200)
c---------------------------------------------------------------------
 1200 format(72('-'))
 1300 format(
     * '--                      POLARIZABILITY',
     * ' TENSOR                         --',/
     * '--                          (atomic',
     * ' units)                            --')
 1400 FORMAT(8X,'XX',10X,'XY',10X,'YY',10X,'XZ',10X,'YZ',10X,'ZZ',/,
     $         6F12.4)
 1500 FORMAT(6X,'Mean polarizability =',f12.4,'   Anisotropy =',f12.4)
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine getgrad(natom,icall,grad)
      implicit real*8(a-h,o-z)
      dimension grad(3*natom*3)
c 
       iadd=1+(icall-1)*natom*3
       call rdgrad(natom,grad(iadd),'save')
c
       thrsh=1.0d-14
ccc    call zerobelowt(natom*3,thrsh,grad(iadd))
       call cutoffshit(natom*3,thrsh,grad(iadd))
c
      end
c======================================================================
      subroutine calcPD(natom,pold,g0,gp,gm,gpp,gpm,gmp,gmm)
      implicit real*8(a-h,o-z)
      dimension pold(3*natom,6)
      dimension g0(3*natom), gp(3*natom,3), gm(3*natom,3)
      dimension gpp(3*natom,3), gpm(3*natom,3)
      dimension gmp(3*natom,3), gmm(3*natom,3)
c
c
         do i=1,natom*3
c
            pold(i,1)= gp(i,1) -2.d0*g0(i) + gm(i,1)  !  xx
            pold(i,4)= gp(i,2) -2.d0*g0(i) + gm(i,2)  !  yy
            pold(i,6)= gp(i,3) -2.d0*g0(i) + gm(i,3)  !  zz
c
            pold(i,2)=gpp(i,1)-gpm(i,1)-gmp(i,1)+gmm(i,1)  !  xy
            pold(i,3)=gpp(i,2)-gpm(i,2)-gmp(i,2)+gmm(i,2)  !  xz
            pold(i,5)=gpp(i,3)-gpm(i,3)-gmp(i,3)+gmm(i,3)  !  yz
         enddo
c
      end
c======================================================================
      subroutine calcP1(natom,pold,g0,gp,gm,gpp,gpm,gmp,gmm)
      implicit real*8(a-h,o-z)
      dimension pold(3*natom,6)
      dimension g0(3*natom), gp(3*natom,3), gm(3*natom,3)
      dimension gpp(3*natom,3), gpm(3*natom,3)
      dimension gmp(3*natom,3), gmm(3*natom,3)
c
c
         do i=1,natom*3
c
            pold(i,1)= gp(i,1) -2.d0*g0(i) + gm(i,1)  !  xx
            pold(i,4)= gp(i,2) -2.d0*g0(i) + gm(i,2)  !  yy
            pold(i,6)= gp(i,3) -2.d0*g0(i) + gm(i,3)  !  zz
c
            pold(i,2)=gpp(i,1)- 2.d0*g0(i)      +gmm(i,1)  !  xy
            pold(i,3)=gpp(i,2)- 2.d0*g0(i)      +gmm(i,2)  !  xz
            pold(i,5)=gpp(i,3)- 2.d0*g0(i)      +gmm(i,3)  !  yz
         enddo
c
      end
c======================================================================
      subroutine calcP2(natom,pold,g0,gp,gm,gpp,gpm,gmp,gmm)
      implicit real*8(a-h,o-z)
      dimension pold(3*natom,6)
      dimension g0(3*natom), gp(3*natom,3), gm(3*natom,3)
      dimension gpp(3*natom,3), gpm(3*natom,3)
      dimension gmp(3*natom,3), gmm(3*natom,3)
c
c
         do i=1,natom*3
c
            pold(i,1)= gp(i,1) -2.d0*g0(i) + gm(i,1)  !  xx
            pold(i,4)= gp(i,2) -2.d0*g0(i) + gm(i,2)  !  yy
            pold(i,6)= gp(i,3) -2.d0*g0(i) + gm(i,3)  !  zz
c
            xy=gp(i,1)+gp(i,2)+gm(i,1)+gm(i,2)
            xz=gp(i,1)+gp(i,3)+gm(i,1)+gm(i,3)
            yz=gp(i,2)+gp(i,3)+gm(i,2)+gm(i,3)
            pold(i,2)=gpp(i,1)+ 2.d0*g0(i) -xy             !  xy
            pold(i,3)=gpp(i,2)+ 2.d0*g0(i) -xz             !  xz
            pold(i,5)=gpp(i,3)+ 2.d0*g0(i) -yz             !  yz
         enddo
c
      end
c======================================================================
      subroutine calcDPD(natom,dipd,pold,g0,gp,gm,gpp,gpm,gmp,gmm,
     *                     diponly)
      implicit real*8(a-h,o-z)
      logical diponly
      dimension dipd(3*natom,3)
      dimension pold(3*natom,6)
      dimension g0(3*natom), gp(3*natom,3), gm(3*natom,3)
      dimension gpp(3*natom,3), gpm(3*natom,3)
      dimension gmp(3*natom,3), gmm(3*natom,3)
c
         do i=1,natom*3
c
            dipd(i,1)= gp(i,1)     -         gm(i,1)  !  x 
            dipd(i,2)= gp(i,2)     -         gm(i,2)  !  y 
            dipd(i,3)= gp(i,3)     -         gm(i,3)  !  z 
c
          if(.not.diponly) then
            pold(i,1)= gp(i,1) -2.d0*g0(i) + gm(i,1)  !  xx
            pold(i,4)= gp(i,2) -2.d0*g0(i) + gm(i,2)  !  yy
            pold(i,6)= gp(i,3) -2.d0*g0(i) + gm(i,3)  !  zz
c
            pold(i,2)=gpp(i,1)-gpm(i,1)-gmp(i,1)+gmm(i,1)  !  xy
            pold(i,3)=gpp(i,2)-gpm(i,2)-gmp(i,2)+gmm(i,2)  !  xz
            pold(i,5)=gpp(i,3)-gpm(i,3)-gmp(i,3)+gmm(i,3)  !  yz
          endif
         enddo
c
      end
c======================================================================
c
      SUBROUTINE RdGRADX(NATOMS,GC,FStat)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads gradient from <grad> file
C
C  ARGUMENTS
C
C  NATOMS  -  number of atoms
C  GC      -  on exit contains current gradient
C  FStat   -  status of file after read (save/delete)
C
C
      REAL*8 GC(3,NATOMS)
      CHARACTER*80 CHAR
      Character*4 FStat
c
      character*256 jobname
      Common /job/jobname,lenJ
C
C
C  open GRAD file
C
      OPEN (UNIT=40,FILE=jobname(1:lenJ)//'.grad',
     $      FORM='FORMATTED',STATUS='OLD',ERR=95)
C
C  locate $gradient section
C
 10   CONTINUE
      READ(40,900,END=96) CHAR
      IF(CHAR(1:9).NE.'$gradient') GO TO 10
C
C  read gradient
C
      DO 20 IAtm=1,NATOMS
      READ(40,910) GC(1,IAtm),GC(2,IAtm),GC(3,IAtm)
 20   CONTINUE
c
C
C  check for end of data
C  (just look for "$")
C
      READ(40,900) CHAR
      If(CHAR(1:1).NE.'$') GO TO 97
c
      IF(FStat.EQ.'save'.OR.FStat.EQ.'keep') THEN
        CLOSE (UNIT=40,STATUS='KEEP')
      ELSE
        CLOSE (UNIT=40,STATUS='DELETE')
      ENDIF
      RETURN
C  ..............................................
C    ERROR SECTION
C
 95   CONTINUE
      Call nerror(6,'File IO routine <RdGRAD>',
     $     'Unable to find <grad> file!',0,0)
c
 96   CONTINUE
      Call nerror(7,'File IO routine <RdGRAD>',
     $      'No $gradient section found on <grad> file!',0,0)
c
 97   CONTINUE
      Call nerror(8,'File IO routine <RdGRAD>',
     $  'No End-Of-File marker ($end) found on <grad> file',0,0)
c
 900  Format(A80)
 910  Format(10X,3F20.14)
c
      END
c =====================================================================
      subroutine print_pold(iout,natom,pold)
      implicit real*8(a-h,o-z)
      dimension pold(3,natom,6)
      write(iout,1000) 
      write(iout,1100) 
      write(iout,1000) 
      write(iout,1200) 
      do iat=1,natom
         write(iout,1300) iat
         write(iout,1400)
     *      ' X',pold( 1 ,iat,1),pold( 1 ,iat,2),pold( 1 ,iat,4),
     $           pold( 1 ,iat,3),pold( 1 ,iat,5),pold( 1 ,iat,6)
         write(iout,1400) 
     *      ' Y',pold( 2 ,iat,1),pold( 2 ,iat,2),pold( 2 ,iat,4),
     $           pold( 2 ,iat,3),pold( 2 ,iat,5),pold( 2 ,iat,6)
         write(iout,1400) 
     *      ' Z',pold( 3 ,iat,1),pold( 3 ,iat,2),pold( 3 ,iat,4),
     $           pold( 3 ,iat,3),pold( 3 ,iat,5),pold( 3 ,iat,6)
      enddo
c
      write(iout,1000) 
c
 1000 format(72('-'))
 1100 format(
     * '--                      POLARIZABILITY',
     * ' DERIVATIVES                    --',/
     * '--                          (atomic',
     * ' units)                            --')
 1200 FORMAT(8X,'XX',10X,'XY',10X,'YY',10X,'XZ',10X,'YZ',10X,'ZZ')
 1300 FORMAT('atom=',i3)
 1400 FORMAT(a2,6(   F10.4,2x))
c
      end
c =====================================================================
      subroutine print_dipold(iout,natom,dipd,pold,diponly)
      implicit real*8(a-h,o-z)
      logical diponly
      dimension dipd(3,natom,3)
      dimension pold(3,natom,6)
c
  100 format(
     * '--                      DIPOLE MOMENT',
     * ' DERIVATIVES                     --',/
     * '--                          (atomic',
     * ' units)                            --')
  200 FORMAT(8X,' X',10X,' Y',10X,' Z' )
      write(iout,1000) 
      write(iout,100 ) 
      write(iout,1000) 
      write(iout,200 ) 
c
      do iat=1,natom
         write(iout,1300) iat
         write(iout,1400)
     *      ' X',dipd( 1 ,iat,1),dipd( 1 ,iat,2),dipd( 1 ,iat,3)
         write(iout,1400) 
     *      ' Y',dipd( 2 ,iat,1),dipd( 2 ,iat,2),dipd( 2 ,iat,3)
         write(iout,1400) 
     *      ' Z',dipd( 3 ,iat,1),dipd( 3 ,iat,2),dipd( 3 ,iat,3)
      enddo
c
      if(diponly) return
c
      write(iout,1000) 
      write(iout,1100) 
      write(iout,1000) 
      write(iout,1200) 
      do iat=1,natom
         write(iout,1300) iat
         write(iout,1400)
     *      ' X',pold( 1 ,iat,1),pold( 1 ,iat,2),pold( 1 ,iat,4),
     $           pold( 1 ,iat,3),pold( 1 ,iat,5),pold( 1 ,iat,6)
         write(iout,1400) 
     *      ' Y',pold( 2 ,iat,1),pold( 2 ,iat,2),pold( 2 ,iat,4),
     $           pold( 2 ,iat,3),pold( 2 ,iat,5),pold( 2 ,iat,6)
         write(iout,1400) 
     *      ' Z',pold( 3 ,iat,1),pold( 3 ,iat,2),pold( 3 ,iat,4),
     $           pold( 3 ,iat,3),pold( 3 ,iat,5),pold( 3 ,iat,6)
      enddo
c
      write(iout,1000) 
c
 1000 format(72('-'))
 1100 format(
     * '--                      POLARIZABILITY',
     * ' DERIVATIVES                    --',/
     * '--                          (atomic',
     * ' units)                            --')
 1200 FORMAT(8X,'XX',10X,'XY',10X,'YY',10X,'XZ',10X,'YZ',10X,'ZZ')
 1300 FORMAT('atom=',i3)
 1400 FORMAT(a2,6(   F10.4,2x))
c
      end
c =====================================================================
      subroutine cutoffshit(ndim,thrsh,array) 
      implicit real*8(a-h,o-z)
      dimension array(ndim)
c
      t1=1.0d0/thrsh
c
      do i=1,ndim
         a=array(i)
         a=a*t1
         ia=int(a)
         a=ia
         array(i)=a*thrsh
      enddo
c
      end
c =====================================================================
