c======================================================================
c the NMR driver
c
      subroutine chshift
c
      use memory
c
      implicit real*8 (a-h,o-z)
      character*256 jobname
      character*20 cdum,wvfnc
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c     common /big/ bl(1)
c     common /intbl/ maxsh,ifp,inx(1)
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
c-------------------------------------------
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
c-------------------------------------------
c put down a memory mark
      call mmark
c
c -- get value of ncall from the depository
c    (has NMR been done before for current wavefunction?)
      call tstival('nnmr',nnmr)
      if(nnmr.eq.1) then
         call getival('nnmr',ncall)
      else
         ncall = 0
      endif
      ncall = ncall + 1
c
c----------------------------------------------------------
c get original values of symmetry pair pointers
c
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c-------------------------------------------
c Restore commons big, index, and ener
c
      call readbl(ncall,ictr)
c
c output : ictr = new address in inx (texas95)
c-----------------------------------------------------------------------
c 2006 : because of segmentation of l shells into s,p
c allocate memory and readin density, eigenvectors and eigenvalues
c
      call read_deneigen
c
c-----------------------------------------------------------------------
c 2006 : segment l shells into s,p
c
c
      call getival('ncs ',ncs)
      call check4lgc(ncs,bl(ictr),lshells,lgcshell,ltype_gc,ldeep_gc)
c
c if l-shells are present then make S,P partitioning
c
      if(lshells.gt.0) then
        call trans_l_2_sp(rhf,bl,bl(ictr),ictr,lshell,'nmrshift')
      endif
c
c if there is a GC basis set but only as deep as 2 and only
c s-functions are GC then we DO NOT profit from GC, especially
c in the gradient and probably hessian. Thus, segment it
c
      if(lgcshell.gt.0.and.ltype_gc.eq.1.and.ldeep_gc.eq.1) then
        call trans_gc_2_sg(rhf,bl,bl(ictr),ictr,nogcsi,'nmrshift')
        write(iout,*)
        write(iout,*)
     * '   General Contracted shells have been segmented for NMR '
        write(iout,*)
     * '   because there were only s-doubly gen.contracted orbitals'
        write(iout,*) ' '
      else if(lgcshell.gt.0.and.iftc1.ne.0) then
        call nerror(8,'FORCE module',
     $        'Basis MUST be segmented for FTC forces',0,0)
      endif
c
c ictr is changed on return if l-shells or gc-shells are present
c
c-----------------------------------------------------------------------
c -- read wavefunction, dft flag and xlow from <control> file
      call getchval('jobname',jobname)
      call rmblan2(jobname,80,lenJ)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,cdum)
      CALL RdDFT(IUnit,idft)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c-------------------------------------------------------------
c
c -- NMR ONLY for closed-shell HF/DFT wavefunctions
      If(wvfnc(1:1).EQ.'U'.OR.wvfnc(2:3).EQ.'MP') Then
        call nerror(10,'NMR module',
     $    'NMR not available for open-shell or MP2',0,0)
      EndIf
c ............................................................
c
c make sure  xlow is in depository
c
      call setrval('xlows',xlow)
c
c-------------------------------------------
c Read-in the NMR options
c
      call nmr_options(ncall,idft)
c-------------------------------------------
c memory checking
c
      call nmr_memcheck(idft,ax)
c-------------------------------------------
c open files for NMR :
c
      call open4nmr(idft)
c-------------------------------------------
c translation
c
      call translat(bl,bl(ictr),'forw')
c-------------------------------------------
c if it is not a first call then D1 is stored
c and the only thing to do is to calculate
c H11n and H01n integrals for requested nuclei
c
CKW   If(ncall.gt.1) go to 9876
c
c WHAT about NOGIAO option ?
c-------------------------------------------
c one-electron part of the NMR calculations:
c
      call shift1(bl,bl(ictr))
c-------------------------------------------
c two-electron part of the NMR calculations:
c
c First construct the second screening density to be
c used in G(D0,g1) and G(D1,g0) (CPHF) calculations
c
      call make_d2screen(bl,bl(ictr))
c
c-------------------------------------------
c Calculate G(D0,g1)
c
      call shift2(idft,ax,bl,bl(ictr))
c-------------------------------------------
c Coupled Perturbed Hartree-Fock
c
      call get_d1nmr(idft,ax,bl,bl(ictr))
c-------------------------------------------
 9876 continue
c-------------------------------------------
c translation back
c
      call translat(bl,bl(ictr),'back')
c-------------------------------------------
c At this point the CPHF eq. have been solved
c and one needs to calculate H11n and H01n 1-el
c integrals and contract them with D0 and D1
c matrices to get finial NMR shielding tensors:
c calculate final shielding tensors :
c
      call nmr_tensor(bl,bl(ictr))
c-------------------------------------------
c -- save current ncall to depository
      call setival('nnmr',ncall)
c-------------------------------------------
c transfer back original L-shell or GC-shell basis set info
c
      if(lshells.gt.0) then
        call trans_sp_2_l(bl,bl(ictr),ictr,lshell)
      endif
      if(lgcshell.gt.0.and.ltype_gc.eq.1.and.ldeep_gc.eq.1) then
        call trans_sg_2_gc(bl,bl(ictr),ictr,lgcshell)
      endif
c
c returns the original value of ictr
c-------------------------------------------
c DO NOT CALL RETALL - annihilates everything !!!!release all memory
c
      call retmark
c-------------------------------------------
      call memory_status('end of nmr shift')
c-------------------------------------------
      call clos4nmr()
c-------------------------------------------
c restore original values of symmetry pair and atom pointers
c
      call setival('SymFunPr',ifp)
      call setival('SymShPr',ifp1)
c-------------------------------------------
      end
c=============================================
      subroutine readbl(ncall,ictr)

      use memory

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c------------------------------------------------
      data last_last /0/
      save last_last
c------------------------------------------------
clublin
      call getival('lcore',lcore)
      call getival('nsym',nsym)
clublin
c      call getival('icon',icond)
      call getival('iout',iout)
      call getival('iarc',iarc)
c
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('ictr',ictr)
c
c------------------------------------------------
c test:
c     write(8,*)'Readbl :ncall=',ncall,' inuc,ibas,na,nbf,nsh,ncf,ncs '
c     write(8,888) inuc,ibas,na,nbf,nsh,ncf,ncs
c 888 format(7(i5,2x))
c------------------------------------------------
c check status of bl(*) :
c
      if(ncall.eq.1) then
         call getmem(0,lastx)
         call retmem(1)
ctest
cc        write(6,*)'readbl; lastx=',lastx,' last_last=',last_last
ctest
         last=lastx
         last_last=last
      else if(ncall.gt.1) then
         call getmem(0,lastx)
         call retmem(1)
ctest
cc        write(6,*)'readbl; lastx=',lastx,' last_last=',last_last
         if(lastx .lt. last_last) then
            call getmem(last_last,lastx)
            lastx=last_last
         endif
         last=lastx
      endif
c------------------------------------------------
c
cc      write(6,88) last,lastx
  88  format('From CHSHIFT : On ENTRY                       :'/
     *       '               value of LAST is           =',i10/
     *       '               first free address in BL   =',i10/)
c------------------------------------------------
c
      end
c=============================================
      subroutine nmr_options(ncall,idft)
c
c read in nmr options
c

      use memory

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /aux/ nreq,nucnr(2000)
c     common /big/ bl(1)
c---
c options for integral calculations (from readin):
      common /intgop/ ncachx,maxpricx,iiii(2)
c---
      common /forcdbl/ thre1,thre2,tchf
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /nmrpri/ maxprice
      common /intlim/ limxmem,limblks,limpair
c---------------------------------------------------------------------
      parameter (nopt=35)
      character*4 word(nopt)
      parameter(ndft=27)
      character*6 dftype,dftxc(ndft)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      character*256 chopv(nopt)
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c
      data ioptyp/1,     2,    11,    12,     1,    11,    1,
     *            1,     0,     1,     1,     1,     3,   11,
     *            0,     0,     0,     0,     0,     0,    0,
     *            0,     0,     0,     0,     0,     0,    0,
     *            0,    21,     0,     3,    11,     1,   11/
c
      data word /'prin','for ','thr1','thr2','iter','thre','nogi',
     *           'ncac','noac','gaug','trsl','rese','limi','lvsh',
     *           'hydr','boro','carb','nitr','oxyg','fluo','sodi',
     *           'magn','alum','sili','phos','sulf','chlo','dumm',
     *           'nocp','malk','vcd ','nics','disp','grid','step'/
      data dftxc/'hfs   ','svwn  ','svwn5 ','hfb   ','bvwn  ','bvwn5 ',
     1           'bp86  ','bpw91 ','blyp  ','bvp86 ','optx  ','ovwn  ',
     2           'ovwn5 ','op86  ','opw91 ','olyp  ','pw91  ','pbe   ',
     3           'o3lyp ','b3lyp ','b3pw91','wah   ','user  ','b97   ',
     4           'b97-1 ','b97-2 ','hcth  '/

c
c---------------------------------------------------------------------
c maxprice (ipay) CAN NOT be changed :
      maxprice=maxpricx
c
c---------------------------------------------------------------------
ckw   if(ncall.eq.1) call infoprt(iout,icond)
                     call infoprt(iout,icond)
c---------------------------------------------------------------------
c thre1=1.0d-10   ! one-el. giao integral threshold
c thre2=1.0d-10   ! two-el. giao integral threshold and
c                   final two_el. integral threshold in direct-CPHF
c thre3=1.0d-8      loose two-el. integral threshold in direct-CPHF
c
c tchf =1.0d-5     accuracy for cphf convergance
c
c idelt=0 or 1 ; do not use or use delta D10 in CPHF
c
c xlvsh=0.0d0  ! levelshift in CPHF, added to all energy differences
c---------------------------------------------------------------------
c default values :
c
      nprint=0
      nreq=0
      thre1=1.0d-10
      thre2=1.0d-10
      thre3=1.0d-8
      ichf=30
      tchf=1.0d-5
      nogiao=0
      ncache=ncachx
      noacce=0
      ngauge=0
      ntrans=0
c2002 idelt =1                 ! delta density ALWAYS used with new CPHF
      ireset=20                ! reset CPHF after iter=20
      limxmem=800 000
c2002 limblks=300
      limblks=0
      limpair=100
      xlvsh=0.0d0
      ncphf=1                  ! cphf is enabled with hybrid dft
      if(idft.eq.22) ncphf=0   ! disable cphf for wah functional
      malk=0                   ! malkin correction is not done
c
c2008
      nics=0
c
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
c
      if(ifound(1).gt.0) nprint=iopv(1,1)
      if(ifound(2).gt.0) then
         iatom1=iopv(1,2)
         iatom2=iopv(2,2)
         if(iatom1.gt.na.or.iatom2.gt.na) go to 35
         do ii=iatom1,iatom2
           nreq=nreq+1
           nucnr(nreq)=ii
         end do
  35     continue
      endif
      if(ifound(3).gt.0) thre1=10.d0**(-ropv(1,3))
      if(ifound(4).gt.0) thre2=10.d0**(-ropv(1,4))
      if(ifound(4).gt.0) thre3=10.d0**(-ropv(2,4))
      if(ifound(5).gt.0) ichf =iopv(1,5)
      if(ifound(6).gt.0) tchf =10.d0**(-ropv(1,6))
      if(ifound(7).gt.0) nogiao=iopv(1,7)
      call setival('nogiao',nogiao)
c -- The nogiao option was fixed up. Needs documenting.    ! GM  Sep 2006
c -- Other options and a lot of crud should also be weeded.
      if(ifound(8).gt.0) ncache=iopv(1,8)
      if(ifound(9).gt.0) noacce=1
      if(ifound(10).gt.0)ngauge=iopv(1,10)
      if(ifound(11).gt.0)ntrans=iopv(1,11)
      if(ifound(12).gt.0)ireset=iopv(1,12)
      if(ifound(13).gt.0) then
         limxmem=iopv(1,13)
         limblks=iopv(2,13)
         limpair=iopv(3,13)
      endif
      if(ifound(14).gt.0) xlvsh=ropv(1,14)
c---------------------------------------------------------------------
c nuclei requested by type :
c
      if(ifound(15).gt.0) then     !     hydrogens
         qatom=1.d0
         call find_nuclei(qatom ,bl(inuc),na,nreq,nucnr)
      endif
c
c nuclei from Boron to Fluorine (5-9)
c
      do ii=16,20
         qatom=dble(ii-11)   ! 5,6,7,8,9
         if(ifound(ii).gt.0) then
            call find_nuclei(qatom ,bl(inuc),na,nreq,nucnr)
         endif
      enddo
c
c nuclei from Sodium to Chlorine (11-17)
c
      do ii=21,27
         qatom=dble(ii-10)   ! 11,12,13,14,15,16,17
         if(ifound(ii).gt.0) then
            call find_nuclei(qatom ,bl(inuc),na,nreq,nucnr)
         endif
      enddo
c
c dummy nuclei
c
         if(ifound(28).gt.0) then
            qatom=0.d0
            call find_nuclei(qatom ,bl(inuc),na,nreq,nucnr)
         endif
      if(ifound(29).gt.0) ncphf=0
      if(ifound(30).gt.0) then
          dftype=chopv(30)(1:6)
          call lowerca2(dftype,6)
          do i=1,ndft
            if(dftype.eq.dftxc(i)) then
              malk=i
              exit
            end if
            malk=ndft+1
          end do
          if(dftype.eq.'      ') then
             call getival('dft',malk)
          endif
          if(malk.gt.ndft) then
          call nerror(9,'nmropt','undefined dft XC potential',malk,malk)
          end if
      endif
      if(ifound(31).gt.0) call setival('vcd',1)   ! switch on VCD
c---------------------------------------------------------------------
c NICS= Nuclear Independet Chemical Shift
c
c NICS option;  3 atoms specifying a plane, if 2 or 3 of them are the same
c               then the program finds the plane that contains most atoms 
c NICS suboptions : 
c
c disp=disp -  displace a given plane by disp ( 0.0 au by default)
c grid=ngrid-  determines the number of NICS points on a plane
c              ngrid*ngrid
c
c    step -  step size for making NICS grid default is 1.0 au
c
      if(ifound(32).gt.0) then
         iatom1=iopv(1,32)
         iatom2=iopv(2,32)
         iatom3=iopv(3,32)
         disp=0.0d0  ! displace P plane parallelly by disp (P')
         ngrid=1     ! grid points on P'; ngrid*ngrid*ngrid + 1
         step=1.0d0  ! step size for grid NICS points
         if(ifound(33).gt.0) disp=ropv(1,33)
         if(ifound(34).gt.0) ngrid=iopv(1,34)
         if(ifound(35).gt.0) step=ropv(1,35)
c
c limit on number of NICS point 10 000 + 1
c
         if(ngrid.gt.100) ngrid=100
c
         nics=ngrid*ngrid+1
c
         call setrval('disp',disp)
         call setival('ngrid',ngrid)
         call setrval('step',step)
         call setival('iatom1',iatom1)
         call setival('iatom2',iatom2)
         call setival('iatom3',iatom3)
      endif
c---------------------------------------------------------------------
c check for linear dependencies in the basis set :
c
      call check4ld(thre2,thre3,tchf)
c---------------------------------------------------------------------
      call setrval('thr1',thre1)
      call setrval('thr2',thre2)
      call setrval('thr3',thre3)
c2002 call setival('delt',idelt)
      call setrval('xlvsh',xlvsh)
      call setival('ncphf',ncphf)
      call setival('malk',malk)
      call setival('nics',nics)
c
      if(ireset.gt.30) ireset=30
c
      call setival('reset',ireset)
c---------------------------------------------------------------------
c print options :
c
c
      write(iout,*) ' '
      do 100 iop=1,nopt
       if(ifound(iop).gt.0) then
          if(iop.eq. 1) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq. 2) write(iout,210) word(iop),iopv(1,2),iopv(2,2)
          if(iop.eq. 3) write(iout,220) word(iop),thre1
          if(iop.eq. 4) write(iout,220) word(iop),thre2,thre3
          if(iop.eq. 5) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq. 6) write(iout,220) word(iop),ropv(1,iop)
          if(iop.eq. 7) write(iout,240) word(iop),' GIAO not used' 
          if(iop.eq.11) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq.12) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq.13) write(iout,210) word(iop),(iopv(ii,iop),ii=1,3)
          if(iop.eq.14) write(iout,220) word(iop),ropv(1,iop)
c
          if(iop.eq.15) write(iout,230) word(iop),'hydrogen'
          if(iop.eq.16) write(iout,230) word(iop),'boron'
          if(iop.eq.17) write(iout,230) word(iop),'carbon'
          if(iop.eq.18) write(iout,230) word(iop),'nitrogen'
          if(iop.eq.19) write(iout,230) word(iop),'oxygen'
          if(iop.eq.20) write(iout,230) word(iop),'fluorine'
          if(iop.eq.21) write(iout,230) word(iop),'sodium'
          if(iop.eq.22) write(iout,230) word(iop),'magnesium'
          if(iop.eq.23) write(iout,230) word(iop),'aluminium'
          if(iop.eq.24) write(iout,230) word(iop),'silicon'
          if(iop.eq.25) write(iout,230) word(iop),'phosphorus'
          if(iop.eq.26) write(iout,230) word(iop),'sulfur'
          if(iop.eq.27) write(iout,230) word(iop),'chlorine'
          if(iop.eq.28) write(iout,230) word(iop),'dummy'
          if(iop.eq.29) write(iout,240) word(iop),' cphf turn off'
          if(iop.eq.30) write(iout,250) word(iop),dftxc(malk)
          if(iop.eq.31) write(iout,240) word(iop),' vcd'
          if(iop.eq.32) write(iout,210) word(iop),(iopv(ii,iop),ii=1,3)
       endif
c
  200 format(/58('-'))
  210 format(' NMR option = ',a4,' is on - value = ',3(i10,2x))
  220 format(' NMR option = ',a4,' is on - value = ',3(e12.5,1x))
  230 format(' NMR option = ',a4,' is on - value = ',a10,
     *                                       ' nuclei of interest' )
  240 format(' NMR option = ',a4,' is on;',a15)
  250 format(' NMR option = ',a4,' is on - XC functional: ',a6)
c
  100 continue
c     write(iout,200)
c---------------------------------------------------------------------
      lopt(1)=nprint
      if(nogiao.eq.1) ntrans=0
      if(nreq.eq.0) then
         nreq=na
         do 85 i=1,nreq
  85     nucnr(i)=i
      else
  186    continue
         do 86 i=1,nreq
         iatom=nucnr(i)
         do 87 j=1,nreq
         if(j.gt.i) then
            jatom=nucnr(j)
            if(jatom.eq.iatom) go to 88
         endif
  87     continue
         go to 89
  88     do 188 ii=j,nreq-1
 188     nucnr(ii)=nucnr(ii+1)
         nreq=nreq-1
         go to 186
  89     continue
  86     continue
      endif
c2002 print *, nreq,' atoms requested ',(nucnr(ii),ii=1,nreq)
c---------------------------------------------------------------------
      write(iout,260) nreq
  260 format(' NMR shift requested for ',i3,' nuclei '        )
c---------------------------------------------------------------------
      end
c=============================================
      subroutine infoprt(iout,icond)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /ganz/ lcore,iov,last,lflag(4),lrest(77)
      common /runtype/ scftype,where
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /fieldx/ xfield,xfgrad,elfiel(9)
c---------------------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
      write(iout,*)
     *  '                The NMR Chemical Shift Module '
      write(iout,*)' '
c
      call f_lush(iout)
c-----------------------------------------------------------
c
c electric field and electric field gradient :
c
      ifield=xfield
      ifgrad=xfgrad
c
ccc   write (iout,*)
      write (iout,223)
c-----------------------------------------------------------------------
      write(iout,301)
  301 format(
     * '----     Calculations of the Nuclear',
     * ' Magnetic Shielding Tensor      ----')
      write(iout,302)
  302 format(
     * '----                        in ppm units',
     * '                            ----')
c-----------------------------------------------------------------------
c----------------------
c     write (iout,280)
c     write (iout,281)
c     write (iout,282)
c----------------------
      if( ifield .eq.1) write (iout,2820) elfiel(1),elfiel(2),elfiel(3)
      if( ifgrad .eq.1) write (iout,2821) elfiel(4),elfiel(5),elfiel(6),
     *                                    elfiel(7),elfiel(8),elfiel(9)
      if(nogiao.eq.0) then
        write (iout,283)
      else
        if(ngauge.eq.0) then
           write (iout,284)
        else
           write (iout,285) ngauge
        endif
      endif
c     write (iout,223)
c------------------------------------------------
      if(lflag(1).eq.0 ) then
           if(noacce.gt.0) then
             write(iout ,110)
           endif
        write(iout ,223)
      endif
c
  110 format('----',18x,'without any acceleration',22x,'----')
c---------------------------------------------------------------------
  223 format(72(1H-))
c---------------------------------------------------------------------
c 280 format
c    *('----',15x,'       Calculations of the        ',15x,'----')
c 281 format
c    *('----',15x,'Nuclear Magnetic Shielding Tensor ',15x,'----')
c 282 format
c    *('----',15x,'          in ppm units            ',15x,'----')
c---------------------------------------------------------------------
 2820 format
     *('----',15x,'   external electric field is     ',15x,'----'/
     *'----',15x,'fx=',f6.4,2x,' fy=',f6.4,2x,' fz=',f6.4,
     *15x,' ----')
 2821 format
     *('----',13x,'external electric field gradient is   ',13x,'----'/
     *'----',15x,'xx=',f6.4,2x,' yx=',f6.4,2x,' yy=',f6.4,
     *15x,' ----'/
     *'----',15x,'zx=',f6.4,2x,' zy=',f6.4,2x,' zz=',f6.4,
     *15x,' ----')
  283 format
     *('----',15x,'      by the GIAO method          ',15x,'----')
  284 format
     *('----',15x,'the CPHF method - gauge orig. (0.0)',14x,'----')
  285 format
     *('----',12x,'the CPHF method - gauge orig. on atom ',i2,12x,'---')
c-----------------------------------------------------------------------
      end
c=============================================
      subroutine nmr_memcheck(idft,ax)

      use memory

      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c  Total memory (in bl(i)) needed for shielding's run is:
c
c  ( ntri=ncf*(ncf+1)/2 , ntri3=3*ntri )
c
c  (1) for  1-el. part ( in shift1.f) :
c
c     call getmem(ntri3,lfoc )  Fock matrices : Fx,Fy,Fz
c     call getmem(ntri3,love )  overlap matrices : Sx,Sy,Sz
c
c Total memory in shift1.f = 6*ntri
c------------------------------------------------------------
c  (2) for  2-el. GIAO ( in shift2.f) :
c
c     call getmem(maxlabels,labels)
c     call getmem(ntri3,lfxyz)  matrices G(d0,giaoX),G(d0,giaoY),G(d0,giaoZ)
c     call getmem(ntri ,lden )  matrices
c     call getmem(ntri3,lfoc )  matrices
c
c     Plus memory needed for 2-el.module (?)
c
c Total memory in shift2.f = 7*ntri
c------------------------------------------------------------
c  (3) for  CPHF       ( in copehf.f) :
c
c     in the main copehf.f    :
c        -------------------------
c     call getmem(ncf*ncf,lvec)
c     call getmem(ncf    ,lval)
c     call getmem(ntri3  ,lhfc)
c     call getmem(ntri3  ,ldn1)
c     call getmem(ntri3  ,lden)
c     call getmem(ntri3  ,lfoc)
c     call getmem(ntri3  ,love)
c     call getmem(12     ,iocc)
c        -------------------------
c        mem_part1=ncf*ncf + ncf + 15*ntri +12
c
c     in chfsol :
c        do CPHF
c        memory for w(3,ncf,ncf) in d1const for ds1d
c        mem_part3=3*ncf*ncf
c
c Total memory in copehf.f = 4*ncf*ncf + ncf +12 + 15*ntri
c-----------------------------------------------------------------------
      common /aux/ nreq,nucnr(2000)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c
      parameter (Zero=0.0d0)
c-----------------------------------------------------------------------
      call getmem(0,lastx)
      call retmem(1)
c
      ntri=ncf*(ncf+1)/2
c
c Total memory in shift1.f = 6*ntri
c Total memory in shift2.f = 7*ntri
c Total memory in copehf.f = 4*ncf*ncf + ncf +12 + 15*ntri
c
      memory_1= 6*ntri                      ! for shift1
      memory_2= 7*ntri                      ! for shift2
      memory_3= 4*ncf*ncf+ncf+ 15*ntri + 12  ! copehf
c
      memory_nmr=max(memory_1,memory_2,memory_3)
c
c add some for integrals
c
      memory_nmr=memory_nmr+100000  ! ??????
c-----------------------------------------------------------------------
      lbldim= lastx + memory_nmr
c
      if(lbldim.ge.lcore) then
         call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
         write(iout,210) lbldim-ioffset,lcore-ioffset
         call nerror(1,'Chshift','Memory: needed and available ',
     *               lbldim-ioffset,lcore-ioffset)
      endif
c
      if(lflag(2).eq.1.or.lflag(3).eq.1) then
         write(iout,209) lbldim-ioffset,lcore-ioffset
         return
      endif
c
  209 format (/1x,' Memory in BL for NMR shift calculations '/
     *            '  needed ',i10,'  available ',i10/)
  210 format (/1x,'common bl too small for shift run, required =',i10,
     *3x,' available =',i10/)
c-----------------------------------------------------------------------
      end
c=============================================
      subroutine translat(bl,inx,text)
c*
      implicit real*8 (a-h,o-z)
      character*4 text
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /translx/ trans(3)
c
      dimension bl(*)
      dimension inx(12,*)
c===========================
c  Check if it is "back" translation
c
      if(text.eq.'back') then
        if(ntrans.eq.1 .or. ntrans.eq.2) then
           trans(1)=-trans(1)
           trans(2)=-trans(2)
           trans(3)=-trans(3)
           call maketrsl(bl(inuc),bl(ibas),inx,trans,na,ncs)
        endif
        return
      endif
c===========================
      iprint=nprint
c===========================
      if(ntrans.eq.0) return
      if(nogiao.ne.0) return
c===========================
c  Find a translation according to two different criterions :
c
c  Translation 1 :
c
c  Find the translation T which makes
c  SUM(ij){ nij*[ RixRj + (Ri-Rj)xT ] } = minimum
c  Ri, Rj are the centers of the basis functions.
c  nij number of matrix elements a given pair Ri,Rj contributes to.
c  The minimum conditions may be written as
c
c  SUM(ij) nij*{ (RixRj)x(Ri-Rj) + [(Ri-Rj)xT]x(Ri-Rj) } = 0
c
c  in a vector form. (x - vector product)
c  This leads to the linear system of equations :
c
c [ (RixRj)xRij ]x = | Yij**2+Zij**2   -Xij*Yij     -Xij*Zij    | X
c [ (RixRj)xRij ]y =-|    -Xij*Yij   Xij**2+Zij**2  -Yij*Zij    | Y
c [ (RixRj)xRij ]z = |    -Xij*Zij     -Yij*Zij   Xij**2+Yij**2 | Z
c
c where the SUM(ij)nij  was omited
c
c  or  Translation 2 :
c
c  Find the translation T which makes RixRj + (Ri-Rj)xT = 0
C  for each pair (Ri,Rj). Than take the weighted average of these
c  Tij translations.
c  Ri, Rj are the centers of the basis functions.
c
c
      if(ntrans.eq.1 .or. ntrans.eq.3) then
        lastna=last
        last=last+na*na
        call findtrs1(bl(lastna+1),na,bl(inuc),
     *                trans,inx,ncs,ncf,ntrans)
c
        last=lastna
      endif
c============================
      if(ntrans.eq.2) then
        lastna=last
        last=last+na*na
        call findtrs2(bl(lastna+1),na,bl(inuc),
     *                trans,inx,ncs,ncf,ntrans)
c
        last=lastna
      endif
c===========================
c  Perform the translation of the molecule if it has been found
c  and if there is no symmetry specifeid
c
      if(ntrans.eq.0) then
        write(iout,*)' Requested translation has not been found '
        write(iout,223)
        return
      endif
      if(nsym.eq.0 .and. ntrans.ne.3) then
        call maketrsl(bl(inuc),bl(ibas),inx,trans,na,ncs)
        write(iout,223)
        write(iout,222)
        write(iout,224) trans
        write(iout,223)
        return
      endif
c
      if(nsym.ne.0 .or. ntrans.eq.3) then
c       use non-zero start for CPHF : (Ri-Rj)xT * D0
c       it will be calculated in the CHFSOL sub.
c
        ntrans=3
        write(iout,223)
        write(iout,222)
        write(iout,224) trans
        write(iout,225)
        write(iout,223)
      endif
c
c===========================
  222 format
     *(8x,'-    The molecule has been translated by the vector',6x,'-')
  223 format(8x,58(1H-))
  224 format(8x,'-    T = (',3(f11.6,1x),')  in a.u. -')
  225 format
     *(8x,'- Above translation used to generate non-zero CPHF start -')
      end
c======================================================================
      subroutine findtrs2(noij,na,datnuc,trans,inx,ncs,ncf,ntrans)
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datnuc(5,*)
      dimension noij(na,na)
      dimension inx(12,*)
      dimension trans(3)
c===========================
c  Find the translation T which makes RixRj + (Ri-Rj)xT = 0
C  for each pair (Ri,Rj). Than take the weighted average of these
c  Tij translations.
c  Ri, Rj are the centers of the basis functions.
c--
      trans(1)=zero
      trans(2)=zero
      trans(3)=zero
c--
c
      do 30 iat=1,na
      do 30 jat=1,iat
      noij(iat,jat)=0
   30 continue
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
c
         jfu=0
         do 50 j=1,i
            len2=inx(3,j)
            jat=inx(2,j)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
c
             len=len1*len2
c
c----
        do 45 jgc=0,ngcjx
c----
c*******************************
            ijel=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
            do 40 j1=1,len2
               jff=jff+1
               if (jff.gt.iff) go to 40
               ijel=ijel+1
   40       continue
            if(iat.gt.jat) then
              noij(iat,jat)=noij(iat,jat)+ijel
            endif
            if(jat.gt.iat) then
              noij(jat,iat)=noij(jat,iat)+ijel
            endif
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c*
c-----------------------
ctest
      ntot=0
      do 1001 iat=1,na
      do 1001 jat=1,iat
      ntot=ntot+noij(iat,jat)
c     write(8,88) iat,jat, noij(iat,jat)
 1001 continue
   88 format('iat,jat=',2i3,' no of elem=',i9)
c     write(8,89) ntot
   89 format(' total number of elements is =',i15)
c-----------------------
      ntotx=0
      ntoty=0
      ntotz=0
      do 100 iat=1,na
      xi=datnuc(2,iat)
      yi=datnuc(3,iat)
      zi=datnuc(4,iat)
        do 100 jat=1,iat
        xj=datnuc(2,jat)
        yj=datnuc(3,jat)
        zj=datnuc(4,jat)
c
        rxrx= yi*zj-zi*yj
        rxry= zi*xj-xi*zj
        rxrz= xi*yj-yi*xj
c
        rmrx= xi-xj
        rmry= yi-yj
        rmrz= zi-zj
c
        itransl=0
c
c    find the translation T that makes Ri x Rj + (Ri-Rj)xT =0
c    for each pair of atoms (centers)
c
        tx=zero
        ty=zero
        tz=zero
c
        nijx=0
        nijy=0
        nijz=0
c
        if(iat.eq.jat) go to 110

c
        arxrx=abs(rxrx)
        arxry=abs(rxry)
        arxrz=abs(rxrz)
c
        armrx=abs(rmrx)
        armry=abs(rmry)
        armrz=abs(rmrz)
c
c case 1:
        if(armrx.GT.acc.and.armry.le.acc.and.armrz.le.acc) then
          div=one/rmrx
          ty=-rxrz*div
          tz=+rxry*div
          itransl=1
          nijy=noij(iat,jat)
          nijz=noij(iat,jat)
          go to 110
        endif
c case 2:
        if(armrx.le.acc.and.armry.GT.acc.and.armrz.le.acc) then
          div=one/rmry
          tx=+rxrz*div
          tz=-rxrx*div
          itransl=1
          nijx=noij(iat,jat)
          nijz=noij(iat,jat)
          go to 110
        endif
c case 3:
        if(armrx.le.acc.and.armry.le.acc.and.armrz.GT.acc) then
          div=one/rmrz
          ty=+rxrx*div
          tx=-rxry*div
          itransl=1
          nijx=noij(iat,jat)
          nijy=noij(iat,jat)
          go to 110
        endif
c case 4:
        if(armrx.GT.acc.and.armry.GT.acc.and.armrz.le.acc) then
          if(arxrx.gt.acc) then
            div=one/rmry
            tx=+rxrz*div   ! +div*(Ty*rmrx) for arb. Ty (=0)
            tz=-rxrx*div
            itransl=1
            nijx=noij(iat,jat)
            nijz=noij(iat,jat)
            go to 110
          endif
          if(arxry.gt.acc) then
            div=one/rmrx
            ty=-rxrz*div   ! +div*(Tx*rmry) for arb. Tx (=0)
            tz=+rxry*div
            itransl=1
            nijy=noij(iat,jat)
            nijz=noij(iat,jat)
            go to 110
          endif
        endif
c case 5:
        if(armrx.GT.acc.and.armry.le.acc.and.armrz.GT.acc) then
          if(arxrx.gt.acc) then
            div=one/rmrz
            tx=-rxry*div   ! +div*(Tz*rmrx) for arb. Tz (=0)
            ty=+rxrx*div
            itransl=1
            nijx=noij(iat,jat)
            nijy=noij(iat,jat)
            go to 110
          endif
          if(arxrz.gt.acc) then
            div=one/rmrx
            tz=+rxry*div   ! +div*(Tx*rmrz) for arb. Tx (=0)
            ty=-rxrz*div
            itransl=1
            nijy=noij(iat,jat)
            nijz=noij(iat,jat)
            go to 110
          endif
        endif
c case 6:
        if(armrx.le.acc.and.armry.GT.acc.and.armrz.GT.acc) then
          if(arxry.gt.acc) then
            div=one/rmrz
            ty=+rxrx*div   ! +div*(Tz*rmry) for arb. Tz (=0)
            tx=-rxry*div
            itransl=1
            nijx=noij(iat,jat)
            nijy=noij(iat,jat)
            go to 110
          endif
          if(arxrz.gt.acc) then
            div=one/rmry
            tz=-rxrx*div   ! +div*(Ty*rmrz) for arb. Ty (=0)
            tx=+rxrz*div
            itransl=1
            nijx=noij(iat,jat)
            nijz=noij(iat,jat)
            go to 110
          endif
        endif
c case 7:
        if(armrx.GT.acc.and.armry.GT.acc.and.armrz.GT.acc) then
          if(arxrz.le.acc) then
            div=one/rmrz
            ty=+rxrx*div   ! +div*(Tz*rmry) for arb. Tz (=0)
            tx=-rxry*div   ! +div*(Tz*rmrx) for arb. Tz (=0)
            itransl=1
            nijx=noij(iat,jat)
            nijy=noij(iat,jat)
            go to 110
          endif
          if(arxry.le.acc) then
            div=one/rmry
            tx=+rxrz*div   ! +div*(Ty*rmrx) for arb. Ty (=0)
            tz=-rxrx*div   ! +div*(Ty*rmrz) for arb. Ty (=0)
            itransl=1
            nijx=noij(iat,jat)
            nijz=noij(iat,jat)
            go to 110
          endif
          if(arxrx.le.acc) then
            div=one/rmrx
            tz=+rxry*div   ! +div*(Tx*rmrz) for arb. Tx (=0)
            ty=-rxrz*div   ! +div*(Tx*rmry) for arb. Tx (=0)
            itransl=1
            nijy=noij(iat,jat)
            nijz=noij(iat,jat)
            go to 110
          endif
        endif
c
  110 continue
c
c average translation weighted by the number of elements
c
        ntotx=ntotx+nijx
        ntoty=ntoty+nijy
        ntotz=ntotz+nijz
c
        ntotal=ntotx+ntoty+ntotz
c
c  elements
c
        eijx=dble(nijx)
        eijy=dble(nijy)
        eijz=dble(nijz)
c
        trans(1)=trans(1)+tx*eijx
        trans(2)=trans(2)+ty*eijy
        trans(3)=trans(3)+tz*eijz
c
  100 continue
c
      if(ntotal.gt.0 ) then
        totx=zero
        toty=zero
        totz=zero
        if(ntotx.gt.0) totx=one/dble(ntotx)
        if(ntoty.gt.0) toty=one/dble(ntoty)
        if(ntotz.gt.0) totz=one/dble(ntotz)
c
        trans(1)=trans(1)*totx
        trans(2)=trans(2)*toty
        trans(3)=trans(3)*totz
      else
        ntrans=0
      endif
c------------------------
      trx=abs(trans(1))
      try=abs(trans(2))
      trz=abs(trans(3))
c
      if(trx.lt.0.1d0 .and. try.lt.0.1d0 .and. trz.lt.0.1d0) ntrans=0
c------------------------
c
c     write(8,881) trans
c 881 format(' translation vector T= ',3(f10.5,1x))
c
c     read(5,555) iread
c 555 format(i1)
c     if(iread.eq.1) then
c        read(5,556) trans(1),trans(2),trans(3)
c        write(8,881) trans
c     endif
c 556 format(3f10.6)
c
      end
c======================================================================
      subroutine findtrs1(noij,na,datnuc,trans,inx,ncs,ncf,ntrans)
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datnuc(5,*)
      dimension noij(na,na)
      dimension inx(12,*)
      dimension trans(3)
      dimension ri(3),rj(3),rxr(3),rij(3),aij(3),aa(3),bb(3,3)
      dimension lvec(3),mvec(3)
c===========================
c  Find the translation T which makes
c  SUM(ij){ nij*[ RixRj + (Ri-Rj)xT ] } = minimum
c  Ri, Rj are the centers of the basis functions.
c  nij number of matrix elements a given pair Ri,Rj contributes to.
c  The minimum conditions may be written as
c
c  SUM(ij) nij*{ (RixRj)x(Ri-Rj) + [(Ri-Rj)xT]x(Ri-Rj) } = 0
c
c  in a vector form. (x - vector product)
c  This leads to the linear system of equations :
c
c [ (RixRj)xRij ]x = | Yij**2+Zij**2   -Xij*Yij     -Xij*Zij    | X
c [ (RixRj)xRij ]y =-|    -Xij*Yij   Xij**2+Zij**2  -Yij*Zij    | Y
c [ (RixRj)xRij ]z = |    -Xij*Zij     -Yij*Zij   Xij**2+Yij**2 | Z
c
c where the SUM(ij)nij  was omited
c===========================
c
      do 30 iat=1,na
      do 30 jat=1,iat
      noij(iat,jat)=0
   30 continue
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
c
         jfu=0
         do 50 j=1,i
            len2=inx(3,j)
            jat=inx(2,j)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
c
             len=len1*len2
c
c----
        do 45 jgc=0,ngcjx
c----
c*******************************
            ijel=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
            do 40 j1=1,len2
               jff=jff+1
               if (jff.gt.iff) go to 40
               ijel=ijel+1
   40       continue
            if(iat.gt.jat) then
              noij(iat,jat)=noij(iat,jat)+ijel
            endif
            if(jat.gt.iat) then
              noij(jat,iat)=noij(jat,iat)+ijel
            endif
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c*
c-----------------------
ctest
      ntot=0
      ntot2=0
      do 1001 iat=1,na
      do 1001 jat=1,iat
      ntot=ntot+noij(iat,jat)
      ntot2=ntot2+noij(iat,jat)*noij(iat,jat)
cccc  write(8,88) iat,jat, noij(iat,jat)
 1001 continue
   88 format('iat,jat=',2i3,' no of elem=',i9)
cccc  write(8,89) ntot
   89 format(' total number of elements is =',i15)
c-----------------------
      do 99 i=1,3
      aa(i)=zero
      trans(i)=zero
      do 99 j=1,3
      bb(i,j)=zero
   99 continue
c-------------
      do 100 iat=1,na
      ri(1)=datnuc(2,iat)
      ri(2)=datnuc(3,iat)
      ri(3)=datnuc(4,iat)
        do 100 jat=1,iat
c
        eij=dble( noij(iat,jat) )
        eij2=eij*eij
c
        rj(1)=datnuc(2,jat)
        rj(2)=datnuc(3,jat)
        rj(3)=datnuc(4,jat)
c
        call vecprod(ri,rj,rxr)
c
        rij(1)=ri(1)-rj(1)
        rij(2)=ri(2)-rj(2)
        rij(3)=ri(3)-rj(3)
c
c-----> call vecprod(rxr,rij,aij)
        call vecprod(rij,rxr,aij)
c
        aa(1)=aa(1)+eij2*aij(1)
        aa(2)=aa(2)+eij2*aij(2)
        aa(3)=aa(3)+eij2*aij(3)
c
c setup the bb matrix
c
        xx=rij(1)*rij(1)
        xy=rij(1)*rij(2)
        xz=rij(1)*rij(3)
        yy=rij(2)*rij(2)
        yz=rij(2)*rij(3)
        zz=rij(3)*rij(3)
c
        bb(1,1)=bb(1,1)+eij2*(yy+zz)
        bb(2,2)=bb(2,2)+eij2*(xx+zz)
        bb(3,3)=bb(3,3)+eij2*(xx+yy)
        bb(1,2)=bb(1,2)-eij2*xy
        bb(1,3)=bb(1,3)-eij2*xz
        bb(2,3)=bb(2,3)-eij2*yz
        bb(2,1)=bb(1,2)
        bb(3,1)=bb(1,3)
        bb(3,2)=bb(2,3)
c
  100 continue
c
c re-normalize aa and bb
c
      renor=one/dble(ntot2)
      do 200 i=1,3
      aa(i)=aa(i)*renor
      do 200 j=1,3
      bb(i,j)=bb(i,j)*renor
  200 continue
c
c inverse of the bb matrix
c
      ndi3=3
      call osinv (bb,ndi3,det,acc,lvec,mvec)
c
cc      if(det.le.acc) write(8,*) 'determinant=',det
c
      call matvec(ndi3,bb,aa,trans)
c
c------------------------
      trx=abs(trans(1))
      try=abs(trans(2))
      trz=abs(trans(3))
c
      if(trx.lt.0.1d0 .and. try.lt.0.1d0 .and. trz.lt.0.1d0) ntrans=0
c------------------------
c
c     write(8,881) trans
c 881 format(' translation vector T2= ',3(f10.5,1x))
c
c     read(5,555) iread
c 555 format(i1)
c     if(iread.eq.1) then
c        read(5,556) trans(1),trans(2),trans(3)
c        write(8,881) trans
c     endif
c 556 format(3f10.6)
c
      end
c======================================================================
      subroutine vecprod(v1,v2,v3)
      implicit real*8 (a-h,o-z)
      dimension v1(3),v2(3),v3(3)
c
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      end
c======================================================================
      subroutine maketrsl(datnuc,datbas,inx,trans,na,ncs)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
      dimension trans(3)
c=======================
      do 100 iat=1,na
      datnuc(2,iat)= datnuc(2,iat) + trans(1)
      datnuc(3,iat)= datnuc(3,iat) + trans(2)
      datnuc(4,iat)= datnuc(4,iat) + trans(3)
  100 continue
c=======================
      do 200 ics=1,ncs
cccc  iat=inx(2,ics)
      ia=inx(1,ics)+1
      ie=inx(5,ics)
        do 250 i=ia,ie
        datbas(11,i)=datbas(11,i)+trans(1)
        datbas(12,i)=datbas(12,i)+trans(2)
        datbas(13,i)=datbas(13,i)+trans(3)
  250   continue
  200 continue
c=======================
      end
c======================================================================
      subroutine print_perf(iout,tot_cpu_m,tot_ela_m)
      implicit real*8 (a-h,o-z)
c print parallel performance results
      call getival('nslv',nslv)
      if(nslv.eq.0) return
c did step do anything in parallel?
      call getrval('cpus',tot_cpu_s)
      if(tot_cpu_s.eq.0.0d0) return
c efficiency = all cputime / elapsed time
      effic=(tot_cpu_m+tot_cpu_s)/tot_ela_m
      write(iout ,300) tot_cpu_s+tot_cpu_m, effic, nslv
  300 format('Total CPU time = ',f9.2,
     *       '        Efficiency = ',f7.3,' on ',i3,' processors')
c timbal is the cumulative wait time of slaves, could have finished
c timbal/nslv earlier
c tred is time spent in the reduce operations
      call getrval('imbal',timbal)
      call getrval('reduce',tred)
      call getrval('bcast',tbcast)
c if time lost is 1min or 1% print
      if(tred+timbal+tbcast.ge.min(tot_cpu_s*0.6,1.0d0))then
         write(iout,310) timbal/(nslv*60.0d0),tbcast/60.00d0,
     *                   tred/60.0d0
      endif
  310 format('Time lost due to imbalance = ',f7.2,
     *   ' Broadcast = ',f7.2,' Reduce = ',f6.2)
      call setrval('cpus',0.0d0)
      call setrval('imbal',0.0d0)
      call setrval('bcast',0.0d0)
      call setrval('reduce',0.0d0)
      call f_lush(iout)
      end
c======================================================================
      subroutine find_nuclei(qatom,datnuc,natoms,nreq,nucnr)
      implicit real*8 (a-h,o-z)
      dimension nucnr(*)
      dimension datnuc(5,*)
c
      do iat=1,natoms
         zza=datnuc(1,iat)
         if(zza.eq.qatom) then
              nreq=nreq+1
              nucnr(nreq)=iat
         endif
      enddo
c
      end
c======================================================================
      subroutine clos4nmr()
      implicit real*8(a-h,o-z)
c
      close (69,status='delete')
c
      close (61,status='delete',err=100)
  100 continue
      close (62,status='delete',err=200)
  200 continue
      close (63,status='delete',err=300)
  300 continue
      close (64,status='delete',err=400)
  400 continue
      close (65,status='delete',err=500)
  500 continue
      close (66,status='delete',err=600)
  600 return
      end
c========================================================
      subroutine open4nmr(idft)
      implicit real*8(a-h,o-z)
      character*256 jobname,scrf,filename
      common /job/jobname,lenJ
c----------------------------------------------------
c open 5 direct access files for NMR calculations
c
c Input
c idft - DFT flag; if >0 also open sequential scratch file
c ndim - record length
c----------------------------------------------------
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      ntri=ncf*(ncf+1)/2
      ndim=ntri*3
c----------------------------------------------------
      lrec = ndim*8          ! record length in bytes
      nfile = 60
c
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,len1)
c
      len = len1 + 6
c
      filename = scrf(1:len1)//'.fock1'
      open (unit=nfile+1,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.over1'
      open (unit=nfile+2,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.dens1'
      open (unit=nfile+3,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.resi1'
      open (unit=nfile+4,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.dcon1'
      open (unit=nfile+5,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
c new file for screening density
c
      ncs8=ncs*ncs*8
      filename = scrf(1:len1)//'.screen'
      open (unit=nfile+9,file=filename(1:len+1),
     *      form='unformatted',access='direct',recl=ncs8)
c
c
      If(idft.gt.0) Then
        filename = scrf(1:len1)//'.temp1'
        open (unit=nfile+6,file=filename(1:len),
     *        form='unformatted',status='unknown')
      EndIf
c
      end
c======================================================================
      subroutine make_d2screen(bl,inx)

      use memory, block => bl

      implicit real*8(a-h,o-z)
      common /cux/ lden,lfoc,lforc,love,lh01,lhfc,lval,lvec,lrem
      dimension bl(*)
      dimension inx(12,*)
c
c The second screening density should be the maximum of
c     Ci*Ca(T)/ei-ea
c over all occupied (i) and virtual (a) orbitals.
c
c
c
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
c
      call getival('lvec',lvec)
      call getival('lval',lval)
c
      call getmem(2,iocc)
c
c----------------------------------------------------
c Read-in no of occ orbitals
c
      np4=4
      call rea (bl(iocc),1 ,np4,'nocc_rhf')
      nocc=bl(iocc)
c------------------------------------------------------------
c calculate and store second screening density :
c
c     call secund(ct2d1 )
c     call elapsec(et2d1)
c
      call getmem(ncs*ncs,lden2)
      call d2screen(inx,ncs,ncf,nocc,bl(lvec),bl(lval),
     *              bl(lden2),bl)
      call save1mat(69, 1 ,ncs*ncs ,bl(lden2))
      call retmem(1)
c------------------------------------------------------------
      call retmem(1)   ! iocc
c------------------------------------------------------------
c     call secund(ct2d2 )
c     call elapsec(et2d2)
c     ct2d=(ct2d2-ct2d1)/60.0d0
c     et2d=(et2d2-et2d1)/60.0d0
c     write(6,*) 'master time for 2dscreen=',ct2d,et2d
c------------------------------------------------------------
c
      end
c======================================================================
      subroutine d2screen(inx,ncs,ncf,nocc,vec,val,den2,bl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension bl(*)
      dimension vec(*),val(*)
      dimension den2(ncs,ncs)
c--------------------------------------------------------------
      call getmem(ncf*ncf,lcica)     ! Ci(T)*Ca
c
      call getmem(ncf,lvocc)
      call getmem(ncf,lvirt)
      call cicamat(ncf,nocc,vec,val,bl(lvocc),bl(lvirt),bl(lcica))
      call retmem(2)
c
      call getmem(ncf,map_fs)
      call setup_densp3(inx,ncf,ncs,bl(lcica),den2,bl(map_fs))
      call retmem(1)   ! map_fs
      call retmem(1)   ! lcica
c--------------------------------------------------------------
      end
c=================================================================
      subroutine cicamat(ncf,nocc,vec,val,vocc,virt,den2)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension vec(ncf,ncf),val(ncf),vocc(ncf),virt(ncf)
      dimension den2(ncf,ncf)
c-----------------------------------------------------------------------
c This routine calculates :  { Cj*Ck(T)/(ej-ek) } matrix
c-----------------------------------------------------------------------
c
c  vec - eigenvectors
c  val - eigenvalues
c
c---------------------------------------------------
      do icf=1,ncf
         vimax=0.d0
         do ior=1,nocc
            vi=abs(vec(icf,ior))
            vimax=max(vimax,vi)
         enddo
         vocc(icf)=vimax
      enddo
c
      do jcf=1,ncf
         vjmax=0.d0
         do jor=nocc+1,ncf
            vj=abs(vec(jcf,jor))
            vjmax=max(vjmax,vj)
         enddo
         virt(jcf)=vjmax
      enddo
c
c
      do icf=1,ncf
         do jcf=1,ncf
            den2(icf,jcf)=vocc(icf)*virt(jcf)
         enddo
      enddo
c
c Rescale by minimum of 1/(ei-ea)
c
      ener=val(nocc)-val(nocc+1)
      ener=abs(ener)
      ener=1.d0/ener
c
      call dscal(ncf*ncf,ener,den2,1)
c
c----------------------------
      end
c====================================================================
      subroutine read_deneigen
      use memory
      implicit real*8 (a-h,o-z)
      common /cux/ lden,lfoc,lforc,love,lh01,lhfc,lval,lvec,lrem
c
      call getival('ncf ',ncf)
c
c memory allocation  :
c
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
c
      call getmem(ncf*ncf,lvec)
      call getmem(ncf    ,lval)
      call getmem(ntri3  ,lden)
c
      call setival('lden',lden)
      call setival('lvec',lvec)
      call setival('lval',lval)
c----------------------------------------------------
c Read-in :
c            eigen vectors and eigenvalues :
c            zero-order density
c
      np4=4
c
      call rea (bl(lvec),ncf*ncf,np4,'evec_rhf')
      call rea (bl(lval),ncf,np4,'eval_rhf')
      call rea (bl(lden),ntri,np4,'den0_rhf')
c
c----------------------------------------------------
c
      end
c====================================================================
