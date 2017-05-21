c=====================================================================
c
      subroutine cphfnum(inp,done)
      use memory
      implicit real*8(a-h,o-z)
      character*256 chopval
c 
c This routine calculates the Ist-order density matrix by numerical
c differentation of the Oth-order densities in the presence of an 
c external electric filed. This can be called :
c the numerical CPHF solver 
c
c can be used for polarizability and hyperpolarizability
c
c  reads the NUMCphf line in the input file and writes options
c  (if any) to the <control> file
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
      dimension field(3)   ! additional aplied external electric field
c
      logical lflag
      logical done
      logical rhf 
c
      parameter (IUnit=1)
      data zero /0.0d0/
      data options/'numc','fdst','file','prin','xxxx'/
      data ioptyp/0,11,21,1,0/
      data ncall /0/
      data delta /0.0d0/
      save ncall
      save delta
      save fiel0
      save wvfnc
c---------------------------------------------------------------------
      save ldens0, lsmat0,ldmp20
      dimension dipole(3)
c---------------------------------------------------------------------
c           nacll filed: x       y      z
c           ------------------------------
c for DX
c             1         +h       0      0
c             2         -h       0      0
c for DY
c             3          0      +h      0
c             4          0      -h      0
c for DZ
c             5          0       0     +h
c             6          0       0     -h
c---------------------------------------------------------------------
      call getival('iout',iout)
      call getival('icon',icon)
      call getival('ncf',ncf)
      call getival('ictr',ictr)
c---------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
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
c---------------------------------------------------------------------
         call open4numc 
c---------------------------------------------------------------------
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
c 
      if(ncall.eq.0) then
c
         call fdcntrl(iunit,9,'$wavefunc',iend)
         if(iend.eq.0)call rdcntrl(IUnit,9,'$wavefunc',3,idum,dum,wvfnc)
c
c get existing external electric field
c
         call getfield0(fiel0)
         write(iout,99) fiel0
 99      format('existing field=',3(f9.4,1x))
c
         call WrPOL(IUnit,field,  2  )
c
c        call zeroit(field,3)
c        call setfield(field)
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
         call setival('print',iprnt)
c
         write(iout,1000)
         write(iout,1155) wvfnc(1:10)
         write(iout,1100) delta
c
c        -------------------------------------------------------------
c        read and save 0-rder Density and Overlap (SCF)
c        needed for hyperpolarizability
c
         call matdef('sma0','s',ncf,ncf)   ! zero-order overlap
         call matdef('den0','s',ncf,ncf)   ! zero-order density
         ldens0= mataddr('den0')
         lsmat0= mataddr('sma0')
c
         np1=1
         np4=4
c
         call matread('sma0',np1,'s matrix')
         call matread('den0',np4,'den0_rhf')
c
         if(wvfnc(1:4).eq.'RMP2') then
            call matdef('dmp20','s',ncf,ncf)   ! zero-order MP2density
            ldmp20= mataddr('dmp20')
            call matread('dmp20',np4,'dens_mp2')
         endif
c        -------------------------------------------------------------
c
      endif
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
 1000 FORMAT(
     */' NUMERICAL CPHF SOLVER FOR AN ELECTRIC FIELD PERTURBATION')
 1100 FORMAT('          Finite Difference Step Size: ',F9.6,/)
 1155 format('                               ',a10)
c---------------------------------------------------------------------
      ncall=ncall+1
c---------------------------------------------------------------------
c set up the field
c
      if(ncall.eq.1) then
         field(1)=delta
         field(2)=zero
         field(3)=zero
      endif
      if(ncall.eq.2) then
         field(1)=-delta
         field(2)=zero
         field(3)=zero
      endif
      if(ncall.eq.3) then
         field(1)=zero
         field(2)=delta 
         field(3)=zero
      endif
      if(ncall.eq.4) then
         field(1)=zero
         field(2)=-delta
         field(3)=zero
      endif
      if(ncall.eq.5) then
         field(1)=zero
         field(2)=zero
         field(3)=delta
      endif
      if(ncall.eq.6) then
         field(1)=zero
         field(2)=zero
         field(3)=-delta
      endif
c end of the story
      if(ncall.eq.7) then
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
c Symmetry ???
c---------------------------------------------------------------------
      if(ncall.lt.7) then
      write(iout,*)'  '
      write(iout,1234) ncall, field
 1234 format(' Step =',i3,' :  SCF with the field=',3(f6.3,2x))
      call f_lush(iout)
      endif
c---------------------------------------------------------------------
c calculate D1 and F1(D1) (or G(D1)) 
c
      if(wvfnc(1:4).eq.'RMP2') then
      call calcD1num(ncf,ncall,delta,'dens_mp2','gock_mp2')
      else
      call calcD1num(ncf,ncall,delta,'den0_rhf','fock_rhf')
      endif
c
c---------------------------------------------------------------------
c the final D1 matrices are on disk : file nfile2
c  1 record  -    Dx
c  2 record  -    Dy
c  3 record  -    Dz
c---------------------------------------------------------------------
      done=.false.
      if(ncall.eq. 7) done=.true.
c---------------------------------------------------------------------
      if(done) then
         nfile1= 81 !   density
         nfile2= 82
         nfile3= 83 !   fock
         nfile4= 84
         call close4numc(nfile1)
         call close4numc(nfile3)
         call setfield(fiel0)
         ncall=0
         write(iout,2000) wvfnc(1:10)
cccc     write(iout,1155) wvfnc(1:10)
      endif
c---------------------------------------------------------------------
 2000 FORMAT(/,
     $' -------------------------------------------------------',/,
     $' ** APPARENTLY SUCCESSFUL Ist-order DENS CALCULATION  **',/,
     $' ----------------       ',a10,' ---------------------',/,
     $' -------------------------------------------------------',/)
c---------------------------------------------------------------------
      IF(done) THEN
c---------------------------------------------------------------------
c   Calculate polarizability and hyperpolarizability
c---------------------------------------------------------------------
c        get back saved 0-rder Density and Overlap 
c                 and write them back 
c            needed for hyperpolarizability
c
         call matconn('smat','s',ncf,ncf,lsmat0)
         call matconn('dens','s',ncf,ncf,ldens0)
c
         np1=1
         np4=4
c
         call matwrite('smat',np1,1,'s matrix')
         call matwrite('dens',np4,0,'den0_rhf')
c
         if(wvfnc(1:4).eq.'RMP2') then
            call matconn('dmp2','s',ncf,ncf,ldmp20)   ! zero-order MP2density
            call matwrite('dmp2',np4,0,'dens_mp2')
            call matrem('dmp2')
         endif
c
         call matrem('dens')
         call matrem('smat')
c---------------------------------------------------------------------
c---------------------------------------------------------------
c put down a memory mark
         call mmark
c---------------------------------------------------------------
c allocate memory for the final polarizability tensor polar(3,3)
c and zero it out
c
         call getmem(9,lpolar)
         call zeroit(bl(lpolar),9)
c
c and for hyperpolarizability beta_xyz
c
         call getmem(10,lhypol) ! 10 independent components
         call zeroit(bl(lhypol),10)
c---------------------------------------------------------------
c THE ONE-ELECTRON PART OF THE POLARIZABILITY CALCULATIONS:
c
c allocate memory needed in 1-el. part of polariz
c for Hx,Hy and Hz matrices (ntri,3) (dipole integrals)
c These are needed for the final polar. Tr(Hx*Dy) and
c If field-dependent basis set is used then for CPHF
c we need also (ix|h0|j)+(i|h0|jx) and S1ij=(i1|j)+(i|j1)
c
c So allocate Fxyz for CPHF and Hxyz for the final polariz.tensor
c For field-independent basis set Fxyz=Hxyz + G(Dxyz,g0)
c
         call getmem(3*ntri,lhxyz)
         call getmem(3*ntri,lfxyz)
         call getmem(3*ntri,lsxyz)
c
         call zeroit(bl(lhxyz),3*ntri)
         call zeroit(bl(lfxyz),3*ntri)
         call zeroit(bl(lsxyz),3*ntri)
c---------------------------------------------------------------
         CALL POLAR1(rhf,bl,bl(ictr),ntri,bl(lsxyz),bl(lhxyz),bl(lfxyz))
c
c        call drumh (bl(lhxyz),ncf,6,'  HX    ')
c        call drumh (bl(lhxyz+ntri),ncf,6,'  Hy    ')
c        call drumh (bl(lhxyz+2*ntri),ncf,6,'  Hz    ')
c---------------------------------------------------------------
c  THE CPHF PART replaced by Numerically obtained D1
c--------------------------------------------------------
         call getmem(ntri*3,ldxyz)
c--------------------------------------------------------
c read from a disk : D1 and F1 matrices 
c
         nfile2=82
         irec=1
         call read1mat(nfile2,irec,ntri,bl(ldxyz))  ! read DX
         irec=2
         call read1mat(nfile2,irec,ntri,bl(ldxyz+ntri))  ! read DY
         irec=3
         call read1mat(nfile2,irec,ntri,bl(ldxyz+2*ntri))  ! read DZ
c
         call zerobelowt(ntri*3,thrsh,bl(ldxyz))
c
         nfile4=84
         irec=1
         call read1mat(nfile4,irec,ntri,bl(lfxyz))  ! read FX
         irec=2
         call read1mat(nfile4,irec,ntri,bl(lfxyz+ntri))  ! read FY
         irec=3
         call read1mat(nfile4,irec,ntri,bl(lfxyz+2*ntri))  ! read fZ
c--------------------------------------------------------
c for MP2 add Hxyz to G(1st-order Dmp2)
c
         if(wvfnc(1:4).eq.'RMP2') then
            call add2vec(3*ntri,bl(lhxyz),bl(lfxyz),bl(lfxyz))
         endif
c--------------------------------------------------------
c At this point we are ready to calculate :
c dipole moment, polarizability and hyperpolarizability
c 
         call matdef('densX','s',ncf,ncf)
         ldensX= mataddr('densX')
         if(wvfnc(1:4).eq.'RMP2') then
            call matread('densX',np4,'dens_mp2')
         else
            call matread('densX',np4,'den0_rhf')
         endif
c--------------------------------------------------------
c
         call dipo_moment(bl(lhxyz),bl(ldensX),ncf,ntri,dipole,dipol)
c
c--------------------------------------------------------
         call matrem('densX')
c--------------------------------------------------------
c
         call pola_tensor(bl(lhxyz),bl(ldxyz),ncf,ntri,bl(lpolar),Am,Aa)
c
c--------------------------------------------------------
c
         call hypol_tensor(bl(lfxyz),bl(ldxyz),ncf,ntri,bl(lhypol),HPm,
     *                     wvfnc)
c
c--------------------------------------------------------
c zero out tensor components which are below threshold
c
         call zerobelowt( 3,thrsh,dipole)
         call zerobelowt( 9,thrsh,bl(lpolar))
         call zerobelowt(10,thrsh,bl(lhypol))
c--------------------------------------------------------
c print out the dipole moment
c
         call printDip(iout,icon,dipole,dipol)
c--------------------------------------------------------
c print out the polarizability and hyperpolarizability 
c
         call printPol(iout,icon,bl(lpolar),am,aa)
c--------------------------------------------------------
c print out the hyperpolarizability tensor
c
         if(wvfnc(1:4).eq.'RMP2') then
            write(iout,911)
            write(icon,911)
ccc11       format(/,18x,'********** APPROXIMATION **********')
  911       format(/,17x,'** APPROXIMATED HYPERPOLARIZABILITY **')
         endif
         call printHypol(iout,icon,bl(lhypol),hpm)
c---------------------------------------------------------------
c release memory
         call retmark
c---------------------------------------------------------------------
         call close4numc(nfile2)
         call close4numc(nfile4)
c---------------------------------------------------------------------
      ENDIF    !     (done) THEN
C
      END
c=======================================================================
      subroutine open4numc 
      implicit real*8(a-h,o-z)
      character*256 jobname,scrf,filename
      common /job/jobname,lenJ
c----------------------------------------------------
c open 4 direct access files 
c----------------------------------------------------
      call getival('ncf ',ncf)
      ntri=ncf*(ncf+1)/2
c----------------------------------------------------
      lrec = ntri*8          ! record length in bytes
      nfile = 80
c
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,len1)
c
      len = len1 + 6
c
      filename = scrf(1:len1)//'.den0f'
      open (unit=nfile+1,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.den1f'
      open (unit=nfile+2,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.foc0f'
      open (unit=nfile+3,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.foc1f'
      open (unit=nfile+4,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      end
c======================================================================
      subroutine close4numc(nfile)
         close (nfile,status='delete')
      end
c======================================================================
c======================================================================
      subroutine calcD1num(ncf,ncall,delta,dens,fock)
      use memory
      implicit real*8(a-h,o-z)
      character*8 dens,fock
c---------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c---------------------------------------------------------------------
c calculate :
c
      denom=2.0d0*delta
      denom=1.d0/denom
      np4=4
c
      nfile1= 81 !   density
      nfile2= 82
      nfile3= 83 !   fock
      nfile4= 84
c
      if(ncall.gt.1 .and. ncall.le.7) then
         call matdef('densp','s',ncf,ncf)
         call matdef('densm','s',ncf,ncf)
         ldensp= mataddr('densp')
         ldensm= mataddr('densm')
         call matdef('fockp','s',ncf,ncf)
         call matdef('fockm','s',ncf,ncf)
         lfockp= mataddr('fockp')
         lfockm= mataddr('fockm')
c
         if(ncall.eq.2) then 
            call matread ('densp',np4,dens) ! read D0(+x field) from previous Step
            call matread ('fockp',np4,fock) ! read F0(+x field) from previous Step
            irec=1
            call save1mat(nfile1,irec,ntri,bl(ldensp))
            call save1mat(nfile3,irec,ntri,bl(lfockp))
         endif
         if(ncall.eq.3) then
            call matread ('densm',np4,dens) ! read D0(-x field) from previous Step
            call matread ('fockm',np4,fock) ! read F0(-x field) from previous Step
            irec=1
            call read1mat(nfile1,irec,ntri,bl(ldensp)) ! read back densP
            call matadd1('densm',-1.0d0,'densp')       ! result in densp
            call matscal('densp',denom)     ! DX matrix
c
            call read1mat(nfile3,irec,ntri,bl(lfockp)) ! read back FockP
            call matadd1('fockm',-1.0d0,'fockp')       ! result in Fockp
            call matscal('fockp',denom)     ! FX matrix
            irec=1
            call save1mat(nfile2,irec,ntri,bl(ldensp)) ! final Dx rec=1 on file2
            call save1mat(nfile4,irec,ntri,bl(lfockp)) ! final Fx rec=1 on file4
         endif
         if(ncall.eq.4) then 
            call matread ('densp',np4,dens) ! read D0(+y field) from previous Step
            call matread ('fockp',np4,fock) ! read F0(+y field) from previous Step
            irec=1
            call save1mat(nfile1,irec,ntri,bl(ldensp))
            call save1mat(nfile3,irec,ntri,bl(lfockp))
         endif
         if(ncall.eq.5) then 
            call matread ('densm',np4,dens) ! read D0(-y field) from previous Step
            call matread ('fockm',np4,fock) ! read F0(-y field) from previous Step
            irec=1
            call read1mat(nfile1,irec,ntri,bl(ldensp))  ! read back densP
            call matadd1('densm',-1.0d0,'densp')        ! result in densp
            call matscal('densp',denom)     ! DY matrix
c
            call read1mat(nfile3,irec,ntri,bl(lfockp))  ! read back fockP
            call matadd1('fockm',-1.0d0,'fockp')        ! result in fockp
            call matscal('fockp',denom)     ! FY matrix
            irec=2
            call save1mat(nfile2,irec,ntri,bl(ldensp))  ! final Dy rec=2 on file2
            call save1mat(nfile4,irec,ntri,bl(lfockp))  ! final Fy rec=2 on file4
         endif
         if(ncall.eq.6) then 
            call matread ('densp',np4,dens) ! read D0(+z field) from previous Step
            call matread ('fockp',np4,fock) ! read F0(+z field) from previous Step
            irec=1
            call save1mat(nfile1,irec,ntri,bl(ldensp))
            call save1mat(nfile3,irec,ntri,bl(lfockp))
         endif
         if(ncall.eq.7) then 
            call matread ('densm',np4,dens) ! read D0(-z field) from previous Step
            call matread ('fockm',np4,fock) ! read F0(-z field) from previous Step
            irec=1
            call read1mat(nfile1,irec,ntri,bl(ldensp))  ! read back densP
            call matadd1('densm',-1.0d0,'densp')        ! result in densp
            call matscal('densp',denom)     ! DY matrix
c
            call read1mat(nfile3,irec,ntri,bl(lfockp))  ! read back fockP
            call matadd1('fockm',-1.0d0,'fockp')        ! result in fockp
            call matscal('fockp',denom)     ! DY matrix
            irec=3
            call save1mat(nfile2,irec,ntri,bl(ldensp))  ! final Dz rec=3 on file2
            call save1mat(nfile4,irec,ntri,bl(lfockp))  ! final Fz rec=3 on file4
         endif
         call matrem('fockm')
         call matrem('fockp')
         call matrem('densm')
         call matrem('densp')
      endif     !    (ncall.gt.1 .and. ncall.le.7)
c---------------------------------------------------------------------
      end
c=======================================================================
