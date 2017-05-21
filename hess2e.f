c===================================================================
c            +-------------------------------------+
c            +                                     +
c            +      two-electron integral          +
c            +        first  derivatives           +
c            +           for hessian               +
c            +                                     +
c            +         Krzysztof Wolinski          +
c            +             July , 2001             +
c            +-------------------------------------+
c=======================================================================
      subroutine hess2e1(idft,ax,rhf,bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      logical do_allat
      character*11 scftype
      character*4 where
c------------------------------------------------
c     common /intbl/ifpp,inxx(100)
c------------------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /neglect/ eps,eps1,epsr
      common /lindvec/ lind,idensp
c------------------------------------------------
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c------------------------------------------------
c Calculates the first derivative 2-el. integrals
c and form a part of the Ist-order Fock matrices G(D,g1)
c
c deriv. integ. will never be re-calculated !!!
c------------------------------------------------
c Hessian options :
c
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension xintxx(9)
c----------------------------------------------------------------------
      dimension natodo(2,1000)
c----------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c----------------------------------------------------------------------
c remember the common block /neglect/ eps,eps1,epsr which has been used
c in the INTEG for SCF. It is needed because the integral threshold may
c be different for zero-order and first-order integrals.
c
      reps=eps
      reps1=eps1
      repsr=epsr
      icacher=icache
c----------------------------------------------------------------------
c prepare the integral system for gradient integrals (iforwhat=4) :
c (like twoint95 for ordinary integrals)
c Full list of parameters for twoint95 (13) is reduced here to only
c FOUR . The remaning parameters are taken from commons set up for
c the SCF . Hopefully, it will be always the case for post-scf calc.
c I want to use this interface for all post-scf applications (nmr,grad
c and second derivatives). Four input parameters will be taken from
c appropriate commons set up when a given application is called. Here
c four input parameters taken from HESS option commons :
c----------------------------------------------------------------------
c get integral threshold :
c
      call getrval('thref1',thref1)    ! int.thres for G(D0,g1)
      iforwhat=3       !  hessian  first derivative integrals
c
      thresh=thref1
c
      call post_scf(bl,inx,ncache,nprint,thresh,
     *              iforwhat,0,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
      if(scftype.ne.'full-direct')then
        call getival('incore',incore)
        if(incore.gt.0)then
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with core storage'
        else
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with disk storage'
        endif
      endif
c
c----------------------------------------------------------------------
c allocate memory for label's arrays :
c
      call getmem(maxlabels,labels)
c----------------------------------------------------------------------
c get addresses for density  and hessian matrix
c (memory allocated in hessana.f)
c
      call getival('ldensi',lden  )
      If(.not.rhf) then
        call getival('ldensB',mden)
      endif
c----------------------------------------------------------------------
      where='forc'
c----------------------------------------------------------------------
c transform ordinary density dens(ij) into denspar(ics,jcs)
c
      call getmem(2*ncs*ncs,idensp)
      call getmem(ncf,map_fs)
c
      If(rhf) Then
        call setup_densp2(inx,ncf,ncs,bl(lden),bl(idensp),bl(map_fs))
c           screening Dens NOT squared
c           read-in second screening density
        call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
      Else
c -- for UHF, allocate extra array for alpha+beta densities
        call getmem(ntri,idAB)
        Call AddVEC(ntri,bl(lden),bl(mden),bl(idAB))
        call setup_densp2(inx,ncf,ncs,bl(idAB),bl(idensp),bl(map_fs))
        call retmem(1)
c           read-in second screening density
        call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
      EndIf
c
       call retmem(1)     ! map_fs
c----------------------------------------------------------------
c allocate memory for mapping array : ncenter(icf)->iatom
c
      call getint(ncf,ncenter)
c
c integer allocation in bl
c----------------------------------------------------------------
c get centers of basis functions :
c
      call get_center(ncs,ncf,inx,bl(ncenter))
c----------------------------------------------------------------------
c Find out for how many atoms Fock-derivative matrices
c can be calculated at once.
c In case of symmetry we may construct Fock-der. matricies
c for symmetry unique atoms only (Abelian). However,
c symmtrization of these F1 matricies takes longer than
c symmetrization of skeleton F1 obtained for for all atoms.
c Thus, one can profit from calculating F1 for symmetry
c unique atoms ONLY if these can be calculated in ONE
c pass. Thus, see first how many passes will be required to
c calculated F1 for all atoms (and symmetrize it from skeleton)
c and how many if F1 are constracted for symmetry uniqe ones
c
      call getival('natunq',natunq)      ! number of unique atoms
      call getival('listunq',listunq)
      call getival('listrea',listreal)
c
      call fder4nat(idft,rhf,na    ,ntri,ntimes_all,natonce_all)
      call fder4nat(idft,rhf,natunq,ntri,ntimes_sym,natonce_sym)
c
c     write(6,*)
c    *' atoms=',na,' unique=',natunq,' ntimes=',ntimes_all,ntimes_sym
c
      do_allat=.true.
      ntimes=ntimes_all
      natonce=natonce_all
      natomsx=na
c
c it works for HF RHF ONLY :
c
      if(rhf.and. idft.eq.0) then
      if(ntimes_sym.lt.ntimes_all) then
         do_allat=.false.
         ntimes=ntimes_sym
         natonce=natonce_sym
         natomsx=natunq
      endif
      endif
c---------------------------------------------------------------------
ctest
c        do_allat=.false.
c        ntimes=ntimes_sym
c        natonce=natonce_sym
c        natomsx=natunq
ctest
c        do_allat=.false.
c        ntimes=2
c        natonce=3
c---------------------------------------------------------------------
c     write(6,*)' do_allat=',do_allat
c---------------------------------------------------------------------
      call getival('iout',iout)
      if(idft.eq.0) then
         if(ntimes.eq.1) write(iout,100)
  100    format(/'Fder will be calculated in 1 pass')
         if(ntimes.gt.1) write(iout,101) ntimes
  101    format(/'Fder will be calculated in ',i2,' passes')
         write(iout,*) '  '
      endif
      if(idft.GT.0) then
         if(ntimes.eq.1) write(iout,200)
  200    format(/
     *   'HF  contributions to Fder will be calculated in 1 pass')
         if(ntimes.gt.1) write(iout,201) ntimes
  201    format(/
     *   'HF  contributions to Fder will be calculated in ',i2,
     *   ' passes')
         write(iout,*) '  '
      endif
      call f_lush(iout)
c---------------------------------------------------------------------
c
      call make_natodo(natomsx,natonce,ntimes,natodo)
c
c     natodo(1,itime)-> nat_begining
c     natodo(2,itime)-> nat_end
c----------------------------------------------------------------------
c dimension of derivative fock   is (3*natonce,ntri )
c----------------------------------------------------------------------
      call secund(tgrad1)
      call elapsec(elagrad1)
c----------------------------------------------------------------------
c calculate two-electron integral derivatives for hessian
c (just like other int_name calls)
c
      do itime=1,ntimes
         natb=natodo(1,itime)
         nate=natodo(2,itime)
         natonce=nate-natb+1
c
         write(iout,1234) itime,natonce
 1234    format('  Constructing G(D0,g1) matrices ',I3,' time for ',I3,
     $          ' atoms')
         write(iout,*) '  '
         call f_lush(iout)
c
         call getmem(3*natonce*ntri,lf1)
         call zeroit(bl(lf1),3*natonce*ntri)
         if(.not.rhf) then
            call getmem(3*natonce*ntri,lf1B)
            call zeroit(bl(lf1B),3*natonce*ntri)
         endif
c
         call para_d0g1(idft,ax,rhf,nblocks,bl,inx,ntri,thref1,
     *                  natb,nate,bl(listreal),do_allat,
     *               bl(idensp),bl(ncenter),
     *             bl(lden),bl(mden),bl(lf1),bl(lf1B),bl(labels))
c
c----------------------------------------------------------------------
c  now we need to re-scale Fock derivative matrices.
c  divide the non-diagonal elements of the Fock matrix by 2
c  and scale back to the normal values from the large ones
c  After that put Fder matrices on a disk in trasposed form (ntri,3,nat)
c
      If(rhf) then
         call resc_fder(ncf,natonce,bl(lf1),thref1)
c
         call getmem(ntri*3,itmp)
        if(do_allat) then
         call trsp_ondsk(61,ntri,natb,nate,bl(lf1),bl(itmp) )
        else
         call trsp_ondsx(61,ntri,natb,nate,bl(lf1),bl(itmp),bl(listunq))
        endif
         call retmem(1)
      else
         call resu_fder(ncf,natonce,bl(lf1),bl(lf1B),thref1)
c
         call getmem(ntri*3,itmp)
        if(do_allat) then
         call trsp_ondsk(61,ntri,natb,nate,bl(lf1),bl(itmp) )
         call trsp_ondsk(71,ntri,natb,nate,bl(lf1b),bl(itmp) )
        else
         call trsp_ondsx(61,ntri,natb,nate,bl(lf1),bl(itmp),bl(listunq))
        call trsp_ondsx(71,ntri,natb,nate,bl(lf1b),bl(itmp),bl(listunq))
        endif
         call retmem(1)
      endif
c----------------------------------------------------------------------
         if(.not.rhf) call retmem(1)
         call retmem(1)
      enddo    ! do itime=1,ntimes
c----------------------------------------------------------------------
c release memory allocated for :
c
      call retmem(2) ! idensp and ncenter
      call retmem(1) ! and for labels
c----------------------------------------------------------------------
c  Return the SCF values of eps,eps1,epsr
c
      eps=reps
      eps1=reps1
      epsr=repsr
      icache=icacher
c----------------------------------------------------------------------
c SYMMETRIZATION of Fder matrices :
c
c If there was symmetry and the Fder1 matrices were constracted for
c all atoms then we symetrize them for exect symmetry from skeletons
c
      call getival('nsym',nsym)
c
      If(idft.EQ.0 .and. nsym.ne.0 .and. do_allat) Then
        call getival('ngener',ngener)
        call getival('SymFunPr',ifp)
        call getival('SymNuPr1',nupair)
        call getival('nsyo',nsyo)
c
        call getmem(6*ntri,lfock1)   ! for f1(ntri,3,2) 2 atoms
        iunit=61
        call fdersymm_hes1(iunit,ngener,ncf,ntri,na,bl(lfock1),
     *                    bl(nsyo),bl(ifp),bl(nupair))
        if(.not.rhf) then
           iunit=71
           call fdersymm_hes1(iunit,ngener,ncf,ntri,na,bl(lfock1),
     *                       bl(nsyo),bl(ifp),bl(nupair))
        endif
        call retmem(1)
      EndIf
c----------------------------------------------------------------------
c If there was symmetry and the Fder1 matrices were constracted for
c symmetry unique atoms only then they are already symmetrized
c (in int_d0g1). Here, generate Fder1 matrices for remaining atoms
c using these for symmetry unique ones. This should NOT be needed
c at all, I do it here to stay consistent with the rest of the program
c for a time being.
c
      If(idft.EQ.0 .and. nsym.ne.0 .and. .not.do_allat) Then
        call getival('SymFunPr',ifp)
        call getival('SymNuPr1',nupair)
        call getival('nsyo',nsyo)
c
        call getmem(6*ntri,lfock1)   ! for f1(ntri,3,2) 2 atoms
        iunit=61
        call fdersymm_hes2(iunit,0, nsym  ,ncf,ntri,na,bl(lfock1),
     *                    bl(nsyo),bl(ifp),bl(nupair),
     *                    natunq,bl(listunq))
        if(.not.rhf) then
           iunit=71
           call fdersymm_hes2(iunit,0, nsym  ,ncf,ntri,na,bl(lfock1),
     *                       bl(nsyo),bl(ifp),bl(nupair),
     *                       natunq,bl(listunq))
        endif
        call retmem(1)
      EndIf
c----------------------------------------------------------------------
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) total/60.0d0,elaps/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for 2e-1st deriv.integ. ='
     *,f8.2,' Elapsed = ',f8.2,' min')
  500 format(58('-'))
c----------------------------------------------------------------------
      end
c=======================================================================
c            +-------------------------------------+
c            +                                     +
c            +      two-electron integral          +
c            +        second derivatives           +
c            +           for hessian               +
c            +                                     +
c            +         Krzysztof Wolinski          +
c            +             July , 2001             +
c            +-------------------------------------+
c=======================================================================
      subroutine hess2e2(idft,ax,rhf,bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      character*11 scftype
      character*4 where
c------------------------------------------------
c     common /intbl/ifpp,inxx(100)
c------------------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /neglect/ eps,eps1,epsr
c------------------------------------------------
      common /lindvec/ lind,idensp
c------------------------------------------------
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c------------------------------------------------
c Hessian (Gradient) options :
c
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension xintxx(9)
c------------------------------------------------
c remember the common block /neglect/ eps,eps1,epsr which has been used
c in the INTEG for SCF. It is needed because the integral threshold may
c be different for zero-order and first-order integrals.
c
      reps=eps
      reps1=eps1
      repsr=epsr
      icacher=icache
c----------------------------------------------------------------------
c prepare the integral system for gradient integrals (iforwhat=4) :
c (like twoint95 for ordinary integrals)
c Full list of parameters for twoint95 (13) is reduced here to only
c FOUR . The remaning parameters are taken from commons set up for
c the SCF . Hopefully, it will be always the case for post-scf calc.
c I want to use this interface for all post-scf applications (nmr,grad
c and second derivatives). Four input parameters will be taken from
c appropriate commons set up when a given application is called. Here
c four input parameters taken from HESS option commons :
c----------------------------------------------------------------------
c get two integral thresholds :
c
      call getrval('threg2',threg2)    ! int.thres for G(D0,g2)
      iforwhat=4                       !  hessian  derivative integrals
c
      thresh=threg2
c
      call post_scf(bl,inx,ncache,nprint,thresh,
     *              iforwhat,0,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
      if(scftype.ne.'full-direct')then
        call getival('incore',incore)
        if(incore.gt.0)then
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with core storage'
        else
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with disk storage'
        endif
      endif
c
c The previous call sets up the local integral system
c the following initializes the slaves if we are running
c in parallel. Many input parameters are taken from commons
c or not sent (see comments in ptwoint.f)
c
c ***************************************************************
c     WARNING!!
c  The total number of atomic centers (na) is set to the number
c  of REAL atoms during Hessian evaluation. (This is because na
c  is added to various common blocks etc... and may also be
c  recovered from the depository in lower-level routines - God
c  alone knows what is going on.) HOWEVER na must correspond to
c  the actual number of atomic centers when the slaves are
c  initialized, else the symmetry-equivalent atoms array is
c  incorrect.  So we temporarily reset na here
c  **************************************************************
c
c -- reset na to total number of centers
      call getival('ndum',ndum)
      na = na+ndum
      call setival('na',na)
c
      call Para_HessInit(ncache,iforwhat,nblocks)
c
c -- now set it back to number of real atoms
      na = na-ndum
      call setival('na',na)
c
c output (see scf) : last line parameters
c----------------------------------------------------------------------
c
c Calculate ONLY :
c
c (1) the second derivative 2-el. integrals c  and contract
c     them with two-paritcle density to form hessian matrix  .
c
c
c  deriv. integ. will never be re-calculated !!!
c----------------------------------------------------------------------
c allocate memory for label's arrays :
c
      call getmem(maxlabels,labels)
c----------------------------------------------------------------------
c get addresses for density  and hessian matrix
c (memory allocated in hessana.f)
c
      call getival('ldensi',lden  )
      call getival('lhess',lhess)
      If(.not.rhf) call getival('ldensB',mden)
c
c dimension of hessian array is (3*natoms,3*natoms)
c----------------------------------------------------------------------
      where='hess'       !   nedded for setup_densp
c----------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c----------------------------------------------------------------------
c transform ordinary density dens(ij) into denspar(ics,jcs)
c
      call getmem(ncs*ncs,idensp)
      call getmem(ncf,map_fs)
c
      If(rhf) Then
        call setup_densp2(inx,ncf,ncs,bl(lden),bl(idensp),bl(map_fs))
      Else
c -- for UHF, allocate extra array for alpha+beta densities
        call getmem(ntri,idAB)
        Call AddVEC(ntri,bl(lden),bl(mden),bl(idAB))
        call setup_densp2(inx,ncf,ncs,bl(idAB),bl(idensp),bl(map_fs))
        call retmem(1)
      EndIf
c
       call retmem(1)     ! map_fs
c----------------------------------------------------------------
c allocate memory for mapping array : ncenter(icf)->iatom
c
      call getint(ncf,ncenter)
c
c integer allocation in bl
c----------------------------------------------------------------
c get centers of basis functions :
c
      call get_center(ncs,ncf,inx,bl(ncenter))
c----------------------------------------------------------------------
      call secund(tgrad1)
      call elapsec(elagrad1)
c----------------------------------------------------------------------
c calculate two-electron integral derivatives for hessian
c (just like other int_name calls)
c
      nat3=na*3
c
      call para_d0g2(idft,ax,rhf,nblocks,bl,inx,ntri,threg2,
     *               bl(idensp),bl(ncenter),
     *               bl(lden),bl(mden),bl(lhess),bl(labels))
c----------------------------------------------------------------------
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) total/60.0d0,elaps/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for 2e-2nd deriv.integ. ='
     *,f8.2,' Elapsed = ',f8.2,' min')
  500 format(58('-'))
c----------------------------------------------------------------------
c release memory allocated idensp and ncenter
c
      call retmem(2)
c
c and for labels
c
      call retmem(1)
c------------------------------------------------
c  Return the SCF values of eps,eps1,epsr
c
      eps=reps
      eps1=reps1
      epsr=repsr
      icache=icacher
c------------------------------------------------
c re-scale this 2-el part of the hessian matrix :
c
      lhess1=lhess-1
      do ii=1,nat3**2
         bl(lhess1+ii)=bl(lhess1+ii)*threg2
      enddo
c------------------------------------------------
      end
c===================================================================
c===========HESSIAN  INTEGRAL DERIVATIVE ROUTINES==================
c
c     It is called from Calcint2 when WHERE='hess'
c
      subroutine hessian_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
c only for first & second derivatives (for use in amshift):
      common /memor4b/ider0,ider1,ider2
c
c dimensions for assembling :
      common /dimasse/ lqijr,lqklr,lqmxr,lij3,lkl3,l3l,lsss
c
      dimension bl(*)
c----------------------------------------------------------
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=lqij
        if(lqkl.gt.lqij) lqmx=lqkl
c----------------------------------------------------------
      ndim=10   ! dimension for buf2(ndim,*) used in first_der
c----------------------------------------------------------
      call getmem(3*nbls,ixabq)
      call conv24r(nbls,npij,bl(idx1),bl(ixab),bl(ixabq))
c----------------------------------------------------------
c find these functions which have non-zero power of x,y,z on A and C

      ijdim=lnij-nfu(nqij)
      kldim=lnkl-nfu(nqkl)
c
      call getint(ijdim,ijvecx)
      call getint(ijdim,ijvecy)
      call getint(ijdim,ijvecz)
c
      call getint(kldim,klvecx)
      call getint(kldim,klvecy)
      call getint(kldim,klvecz)
c
      call find_non0(nqij,lnij,nqkl,lnkl,
     *               nijx,nijy,nijz,nklx,nkly,nklz,
     *               bl(ijvecx),bl(ijvecy),bl(ijvecz),
     *               bl(klvecx),bl(klvecy),bl(klvecz) )
c----------------------------------------------------------
      ibeg =ibuf
      ider2=ibuf2
      incr45=45*ngcd*nbls*lnij*lnkl
      ider1=ider2+incr45
      incr9 = 9*ngcd*nbls*lnij*lnkl
      ider0=ider1+incr9
c
      call first_der(ngcd,nbls,bl(ibeg),ndim,
     *               lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *               bl(ider1),bl(ider0),bl(ixabq),
     *               nijx,nijy,nijz,nklx,nkly,nklz,
     *               bl(ijvecx),bl(ijvecy),bl(ijvecz),
     *               bl(klvecx),bl(klvecy),bl(klvecz) )
c
      call secnd_der(ngcd,nbls,bl(ibeg),lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *               bl(ider2),bl(ixabq))
c
c--------------------------------------------------------------
c   NOTHING IS READY FOR L-shells and never will be !!!
c------
c     if(lshellt.eq.0) go to 100
c--------------------------------------------------------------
c
c 100 continue
c--------------------------------------------------------------
c
      call retmem(6)
      call retmem(1)
c--------------------------------------------------------------
      end
c=================================================================
      subroutine secnd_der(ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                     nqij,nqkl,der2,xab)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
cccc  dimension der0(nbls,lnij,lnkl,ngcd)     these two are constracted
cccc  dimension der1(9,nbls,lnij,lnkl,ngcd)   when first_der is called
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c----------------------------------------------------------------------
c This routine calculates UNFINISHED hessian derivatives (second order)
c of the type
c----------------------------------------------------------------------
c d2/dAxdAy (a+b,s|c+d,s)=    (2a)*(2a)*(a+b+1x+1y,s|c+d,s)
c                          -n_ab_x*(2a)*(a+b-1x+1y,s|c+d,s)
c                          -n_ab_y*(2a)*(a+b+1x-1y,s|c+d,s)
c                        +n_ab_x*n_ab_y*(a+b-1x-1y,s|c+d,s)
c
c d2/dAxdBy (a+b,s|c+d,s)=
c          (2a)*(2b)*[ (a+b+1x+1y,s|c+d,s) + (Ay-By)*(a+b+1x,s|c+d,s) ]
c       -n_ab_x*(2b)*[ (a+b-1x+1y,s|c+d,s) + (Ay-By)*(a+b-1x,s|c+d,s) ]
c
c d2/dAxdCy (a+b,s|c+d,s)=    (2a)*(2c)*(a+b+1x,s|c+d+1y,s)
c                          -n_ab_x*(2c)*(a+b-1x,s|c+d+1y,s)
c                          -n_cd_y*(2a)*(a+b+1x,s|c+d-1y,s)
c                        +n_ab_x*n_cd_y*(a+b-1x,s|c+d-1y,s)
c
c d2/dAxdDy (a+b,s|c+d,s)= - [ d2/dAxdAy (a+b,s|c+d,s)
c                             +d2/dAxdBy (a+b,s|c+d,s)
c                             +d2/dAxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c d2/dBxdAy (a+b,s|c+d,s)= d2/dAydBx (a+b,s|c+d,s) see above
c
c d2/dBxdBy (a+b,s|c+d,s)=          (2b)*(2b)*(a+b+1x+1y,s|c+d,s)
c                          +(Ay-By)*(2b)*(2b)*(a+b+1x,s|c+d,s)
c                          +(Ax-Bx)*(2b)*(2b)*(a+b+1y,s|c+d,s)
c                  +(Ax-Bx)*(Ay-By)*(2b)*(2b)*(a+b,s|c+d,s)
c                          +d(Ax-Bx)/dBy*(2b)*(a+b,s|c+d,s)
c
c d2/dBxdCy (a+b,s|c+d,s)=          (2b)*(2c)*(a+b+1x,s|c+d+1y,s)
c                           -n_cd_y*(2b)*     (a+b+1x,s|c+d-1y,s)
c                          +(Ax-Bx)*(2b)*(2c)*(a+b,s|c+d+1y,s)
c                        -(Ax-Bx)*n_cd_y*(2b)*(a+b,s|c+d-1y,s)
c
c d2/dBxdDy (a+b,s|c+d,s)= - [ d2/dBxdAy (a+b,s|c+d,s)
c                             +d2/dBxdBy (a+b,s|c+d,s)
c                             +d2/dBxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c d2/dCxdAy (a+b,s|c+d,s)= d2/dAydCx (a+b,s|c+d,s)  see above
c d2/dCxdBy (a+b,s|c+d,s)= d2/dBydCx (a+b,s|c+d,s)  see above
c
c d2/dCxdCy (a+b,s|c+d,s)=    (2c)*(2c)*(a+b,s|c+d+1x+1y,s)
c                          -n_cd_x*(2c)*(a+b,s|c+d-1x+1y,s)
c                          -n_cd_y*(2c)*(a+b,s|c+d+1x-1y,s)
c                        +n_cd_x*n_cd_y*(a+b,s|c+d-1x-1y,s)
c
c d2/dCxdDy (a+b,s|c+d,s)= - [ d2/dCxdAy (a+b,s|c+d,s)
c                             +d2/dCxdBy (a+b,s|c+d,s)
c                             +d2/dCxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c d2/dDxdAy (a+b,s|c+d,s)= d2/dAydDx (a+b,s|c+d,s)  see above
c d2/dDxdBy (a+b,s|c+d,s)= d2/dBydDx (a+b,s|c+d,s)  see above
c d2/dDxdCy (a+b,s|c+d,s)= d2/dCydDx (a+b,s|c+d,s)  see above
c
c d2/dDxdDy (a+b,s|c+d,s)= - [ d2/dDxdAy (a+b,s|c+d,s)
c                             +d2/dDxdBy (a+b,s|c+d,s)
c                             +d2/dDxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c Order for derivatives:       AA AB AC AD
c                                 BB BC BD
c                                    CC CD
c                                       DD
c
c       Block AA :         Block AB :          Block AC :     Block AD :
c  1 -  d2/dAxdAx     7 -  d2/dAxdBx     16 -  d2/dAxdCx      trans.inv.
c  2 -  d2/dAxdAy     8 -  d2/dAxdBy     17 -  d2/dAxdCy
c  3 -  d2/dAxdAz     9 -  d2/dAxdBz     18 -  d2/dAxdCz
c  4 -  d2/dAydAy    10 -  d2/dAydBx     19 -  d2/dAydCx
c  5 -  d2/dAydAz    11 -  d2/dAydBy     20 -  d2/dAydCy
c  6 -  d2/dAzdAz    12 -  d2/dAydBz     21 -  d2/dAydCz
c                    13 -  d2/dAzdBx     22 -  d2/dAzdCx
c                    14 -  d2/dAzdBy     23 -  d2/dAzdCy
c                    15 -  d2/dAzdBz     24 -  d2/dAzdCz
c
c       Block BA :         Block BB :          Block BC :     Block BD :
c                    25 -  d2/dBxdBx     31 -  d2/dBxdCx      trans.inv.
c                    26 -  d2/dBxdBy     32 -  d2/dBxdCy
c                    27 -  d2/dBxdBz     33 -  d2/dBxdCz
c                    28 -  d2/dBydBy     34 -  d2/dBydCx
c                    29 -  d2/dBydBz     35 -  d2/dBydCy
c                    30 -  d2/dBzdBz     36 -  d2/dBydCz
c                                        37 -  d2/dBzdCx
c                                        38 -  d2/dBzdCy
c                                        39 -  d2/dBzdCz
c
c       Block CA :         Block CB :          Block CC :     Block CD :
c                                        40 -  d2/dCxdCx      trans.inv.
c                                        41 -  d2/dCxdCy
c                                        42 -  d2/dCxdCz
c                                        43 -  d2/dCydCy
c                                        44 -  d2/dCydCz
c                                        45 -  d2/dCzdCz
c----------------------------------------------------------------------
c INPUT buf2(1,nbls,lnijr,lnklr,ngcd) - ordinary 2-el.integ.
c
c  buf2(2,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_a
c  buf2(3,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_b
c  buf2(4,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_c
c
c  buf2(5,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_a*exp_b
c  buf2(6,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_a*exp_c
c  buf2(7,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_b*exp_c
c
c  buf2(8,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_a*exp_a
c  buf2(9,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_b*exp_b
c  buf2(10nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_c*exp_c
c
c OUTPUT :
c    der000=der0(nbls,lnij,lnkl,ngcd) - ordinary integr.
c    needed in shifting for d/dA,d/dB (1->2) and d/dC (3->4)
c
c    der1AX=der1(1,nbls,lnij,lnkl,ngcd),
c    der1BX=der1(2,nbls,lnij,lnkl,ngcd)
c    der1CX=der1(3,nbls,lnij,lnkl,ngcd),
c
c    der1AY=der1(4,nbls,lnij,lnkl,ngcd),
c    der1BY=der1(5,nbls,lnij,lnkl,ngcd),
c    der1CY=der1(6,nbls,lnij,lnkl,ngcd),
c
c    der1AZ=der1(7,nbls,lnij,lnkl,ngcd),
c    der1BZ=der1(8,nbls,lnij,lnkl,ngcd),
c    der1CZ=der1(9,nbls,lnij,lnkl,ngcd),
c
c these above are not calculated here,
c they where calculated in first_der
c.................................................
c
c    der2_AxAx=der2(1,nbls,lnij,lnkl,ngcd),
c    der2_AxAy=der2(2,nbls,lnij,lnkl,ngcd),
c    der2_AxAz=der2(3,nbls,lnij,lnkl,ngcd),
c    der2_AyAy=der2(4,nbls,lnij,lnkl,ngcd),
c    der2_AyAz=der2(5,nbls,lnij,lnkl,ngcd),
c    der2_AzAz=der2(6,nbls,lnij,lnkl,ngcd),
c
c    der2_AxBx= der2(7,nbls,lnij,lnkl,ngcd),
c    der2_AxBy= der2(8,nbls,lnij,lnkl,ngcd),
c    der2_AxBz= der2(9,nbls,lnij,lnkl,ngcd),
c    der2_AyBx=der2(10,nbls,lnij,lnkl,ngcd),
c    der2_AyBy=der2(11,nbls,lnij,lnkl,ngcd),
c    der2_AyBz=der2(12,nbls,lnij,lnkl,ngcd),
c    der2_AzBx=der2(13,nbls,lnij,lnkl,ngcd),
c    der2_AzBy=der2(14,nbls,lnij,lnkl,ngcd),
c    der2_AzBz=der2(15,nbls,lnij,lnkl,ngcd),
c
c    der2_AxCx=der2(16,nbls,lnij,lnkl,ngcd),
c    der2_AxCy=der2(17,nbls,lnij,lnkl,ngcd),
c    der2_AxCz=der2(18,nbls,lnij,lnkl,ngcd),
c    der2_AyCx=der2(19,nbls,lnij,lnkl,ngcd),
c    der2_AyCy=der2(20,nbls,lnij,lnkl,ngcd),
c    der2_AyCz=der2(21,nbls,lnij,lnkl,ngcd),
c    der2_AzCx=der2(22,nbls,lnij,lnkl,ngcd),
c    der2_AzCy=der2(23,nbls,lnij,lnkl,ngcd),
c    der2_AzCz=der2(24,nbls,lnij,lnkl,ngcd),
c
c    der2_BxBx=der2(25,nbls,lnij,lnkl,ngcd),
c    der2_BxBy=der2(26,nbls,lnij,lnkl,ngcd),
c    der2_BxBz=der2(27,nbls,lnij,lnkl,ngcd),
c    der2_ByBy=der2(28,nbls,lnij,lnkl,ngcd),
c    der2_ByBz=der2(29,nbls,lnij,lnkl,ngcd),
c    der2_BzBz=der2(30,nbls,lnij,lnkl,ngcd),
c
c    der2_BxCx=der2(31,nbls,lnij,lnkl,ngcd),
c    der2_BxCy=der2(32,nbls,lnij,lnkl,ngcd),
c    der2_BxCz=der2(33,nbls,lnij,lnkl,ngcd),
c    der2_ByCx=der2(34,nbls,lnij,lnkl,ngcd),
c    der2_ByCy=der2(35,nbls,lnij,lnkl,ngcd),
c    der2_ByCz=der2(36,nbls,lnij,lnkl,ngcd),
c    der2_BzCx=der2(37,nbls,lnij,lnkl,ngcd),
c    der2_BzCy=der2(38,nbls,lnij,lnkl,ngcd),
c    der2_BzCz=der2(39,nbls,lnij,lnkl,ngcd),
c
c    der2_CxCx=der2(40,nbls,lnij,lnkl,ngcd),
c    der2_CxCy=der2(41,nbls,lnij,lnkl,ngcd),
c    der2_CxCz=der2(42,nbls,lnij,lnkl,ngcd),
c    der2_CyCy=der2(43,nbls,lnij,lnkl,ngcd),
c    der2_CyCz=der2(44,nbls,lnij,lnkl,ngcd),
c    der2_CzCz=der2(45,nbls,lnij,lnkl,ngcd),
c-------------------------------------------------------------------
c Second derivatives:
c
c Block AA (1-6)
c Block BB (25-30)
c Block CC (40-45)
c
        nder_aa=1
        nder_bb=25
        nder_cc=40
        do icart=1,3
           do jcart=icart,3
              call block_aa(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_aa,der2)
              call block_bb(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_bb,der2,xab)
              call block_cc(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_cc,der2)
              nder_aa=nder_aa+1
              nder_bb=nder_bb+1
              nder_cc=nder_cc+1
           enddo
        enddo
c
c Block AB (7-15)
c Block AC (16-24)
c Block BC (31-39)
c
        nder_ab=7
        nder_ac=16
        nder_bc=31
        do icart=1,3
           do jcart=1,3
              call block_ab(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_ab,der2,xab)
              call block_ac(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_ac,der2)
              call block_bc(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_bc,der2,xab)
              nder_ab=nder_ab+1
              nder_ac=nder_ac+1
              nder_bc=nder_bc+1
           enddo
        enddo
c
      end
c=================================================================
      subroutine block_aa(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_aa,der2)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      do 200 ij=nfu(nqij)+1,lnij
c
      ij_pi=npxyz(icart,ij)         ! funct. ( i+j + 1icart )
      ij_pj=npxyz(jcart,ij)         ! funct. ( i+j + 1jcart )
      ij_pi_pj=npxyz(jcart,ij_pi)   ! funct. ( i+j + 1ic + 1-jc)
c----------------------------------
c Three "powers" are needed :
c     n_ij_00_i=nia(icart,ij)
c     n_ij_pi_j=nia(jcart,ij_pi)
c     n_ij_mi_j=nia(jcart,ij_mi)
c----------------------------------
c
      n_ij_pi_j=nia(jcart,ij_pi)    ! jcart-power of (i+j+1ic |
c
      n_ij_00_i=nia(icart,ij)       ! icart-power of (i+j+0    |
      if(n_ij_00_i.gt.0) then
         ij_mi=nmxyz(icart,ij)      ! funct. (i+j-1ic |
      else
         ij_mi=0
      endif
c
      if(ij_mi.gt.0) then
         n_ij_mi_j=nia(jcart,ij_mi) ! jcart-power of (i+j-1ic
      else
         n_ij_mi_j=0
      endif
c
c first we added now we substract:
c
      ij_mi_pj=nmxyz(icart,ij_pj)   ! funct. (i+j -1ic +1jc)
      ij_pi_mj=nmxyz(jcart,ij_pi)   ! funct. (i+j +1ic -1jc)
c
      ij_mi_mj=0
      if(ij_mi.gt.0) ij_mi_mj=nmxyz(jcart,ij_mi)  ! funct. (i+j-1i-1j|
c                 or ij_mi_mj=nmxyz(icart,ij_mj)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_ab_pipj=buf2(8,ijkl,ij_pi_pj, kl,iqu)
c       two_a_mi_pj =buf2(2,ijkl,ij_mi_pj, kl,iqu)
c       two_a_pi_mj =buf2(2,ijkl,ij_pi_mj, kl,iqu)
c       two_0_mi_mj =buf2(1,ijkl,ij_mi_mj, kl,iqu)
c
      der=buf2(ijkl,ij_pi_pj, kl,iqu,8)
c
      if(ij_mi_pj.gt.0) then
         der=der-n_ij_00_i*buf2(ijkl,ij_mi_pj, kl,iqu,2)
      endif
      if(ij_pi_mj.gt.0) then
         der=der-n_ij_pi_j*buf2(ijkl,ij_pi_mj, kl,iqu,2)
      endif
      if(ij_mi_mj.gt.0) then
         der=der+n_ij_00_i*n_ij_mi_j*buf2(ijkl,ij_mi_mj, kl,iqu,1)
      endif
c
      der2(nder_aa,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_cc(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_cc,der2)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      kl_pi=npxyz(icart,kl)
      kl_pj=npxyz(jcart,kl)
      kl_pi_pj=npxyz(jcart,kl_pi)
c
      n_kl_pi_j=nia(jcart,kl_pi)
c
      n_kl_00_i=nia(icart,kl)
      if(n_kl_00_i.gt.0) then
         kl_mi=nmxyz(icart,kl)
      else
         kl_mi=0
      endif
c
      if(kl_mi.gt.0) then
         n_kl_mi_j=nia(jcart,kl_mi)
      else
         n_kl_mi_j=0
      endif
c
c first add than substract:
      kl_mi_pj=nmxyz(icart,kl_pj)
      kl_pi_mj=nmxyz(jcart,kl_pi)
c
      kl_mi_mj=0
      if(kl_mi.gt.0) kl_mi_mj=nmxyz(jcart,kl_mi)
c                 or kl_mi_mj=nmxyz(icart,kl_mj)
c
      do 200 ij=nfu(nqij)+1,lnij
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c
c       four_c_pipj=buf2(10,ijkl,ij,kl_pi_pj,iqu)
c       two_c_mi_pj=buf2( 4,ijkl,ij,kl_mi_pj,iqu)
c       two_c_pi_mj=buf2( 4,ijkl,ij,kl_pi_mj,iqu)
c       two_0_mi_mj=buf2( 1,ijkl,ij,kl_mi_mj,iqu)
c
      der=buf2(ijkl,ij,kl_pi_pj,iqu,10)
      if(kl_mi_pj.gt.0) then
        der=der-n_kl_00_i*buf2(ijkl,ij,kl_mi_pj,iqu,4)
      endif
      if(kl_pi_mj.gt.0) then
        der=der-n_kl_pi_j*buf2(ijkl,ij,kl_pi_mj,iqu,4)
      endif
      if(kl_mi_mj.gt.0) then
        der=der+n_kl_00_i*n_kl_mi_j*buf2(ijkl,ij,kl_mi_mj,iqu,1)
      endif
c
      der2(nder_cc,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_bb(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_bb,der2,xab)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
ctest
c     write(6,60) icart,jcart,nder_bb
c 60  format('from BB: icart,jcart=',2i3,' no=',i3)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      do 200 ij=nfu(nqij)+1,lnij
c
      ij_pi=npxyz(icart,ij)
      ij_pj=npxyz(jcart,ij)
      ij_pi_pj=npxyz(jcart,ij_pi)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_b_pipj=buf2(9,ijkl,ij_pi_pj, kl,iqu)
c       four_b_pi  =buf2(9,ijkl,ij_pi   , kl,iqu)
c       four_b_pj  =buf2(9,ijkl,ij_pj   , kl,iqu)
c       four_b_0   =buf2(9,ijkl,ij      , kl,iqu)
c       two_b_0    =buf2(3,ijkl,ij      , kl,iqu)
c
      der=                  buf2(ijkl,ij_pi_pj,kl,iqu,9)
     *              +xab(ijkl,jcart)*buf2(ijkl,ij_pi,kl,iqu,9)
     *              +xab(ijkl,icart)*buf2(ijkl,ij_pj,kl,iqu,9)
     *   +xab(ijkl,icart)*xab(ijkl,jcart)*buf2(ijkl,ij,kl,iqu,9)

      if(jcart.eq.icart) der=der-buf2(ijkl,ij,kl,iqu,3)
c
      der2(nder_bb,ijkl,ij,kl,iqu)=der
c
c        write(6,66) nder_bb, ij,kl
c  66    format('no=',i3,'ij,kl=',2i4)
c        write(6,67) '(2b)2 int. ++ =',buf2(9,ijkl,ij_pi_pj,kl,iqu)
c        write(6,67) '(2b)2 int. +0 =',buf2(9,ijkl,ij_pi   ,kl,iqu)
c        write(6,67) '(2b)2 int. 0+ =',buf2(9,ijkl,ij_pj   ,kl,iqu)
c        write(6,67) '(2b)2 int. 00 =',buf2(9,ijkl,ij      ,kl,iqu)
c        write(6,67) '(2b)1 int. 00 =',buf2(3,ijkl,ij      ,kl,iqu)
c        write(6,67) '  final deriv =',der
c 67     format(a15,f12.7)
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_ab(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_ab,der2,xab)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
ctest
c     write(6,60) icart,jcart,nder_ab
c 60  format('from AB: icart,jcart=',2i3,' no=',i3)
c
c
      do 200 kl=nfu(nqkl)+1,lnkl
      do 200 ij=nfu(nqij)+1,lnij
c
      ij_pi=npxyz(icart,ij)
      ij_pj=npxyz(jcart,ij)
      n_ij_i=nia(icart,ij)
c
      ij_mi=0
      if(n_ij_i.gt.0) ij_mi=nmxyz(icart,ij)
c
      ij_pi_pj=npxyz(jcart,ij_pi)
c first add than substract:
      ij_mi_pj=nmxyz(icart,ij_pj)
c
c     ij_mi_mj=0
c     if(ij_mi.gt.0) ij_mi_mj=nmxyz(jcart,ij_mi)
c                 or ij_mi_mj=nmxyz(icart,ij_mj)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_ab_pipj=buf2(5,ijkl,ij_pi_pj, kl,iqu)
c       four_ab_pi  =buf2(5,ijkl,ij_pi   , kl,iqu)
c
c       two_b_mi_pj=buf2(3,ijkl,ij_mi_pj, kl,iqu)
c       two_b_mi   =buf2(3,ijkl,ij_mi   , kl,iqu)
c
      der=                buf2(ijkl,ij_pi_pj, kl,iqu,5)
     *         +xab(ijkl,jcart)*buf2(ijkl,ij_pi , kl,iqu,5)
c
      if(n_ij_i.gt.0) then
        if(ij_mi_pj.gt.0) then
          der=der-n_ij_i*buf2(ijkl,ij_mi_pj, kl,iqu,3)
        endif
        if(ij_mi   .gt.0) then
          der=der -n_ij_i*xab(ijkl,jcart)*buf2(ijkl,ij_mi , kl,iqu,3)
        endif
      endif
c
c        write(6,66) nder_ab, ij,kl
c  66    format('no=',i3,'ij,kl=',2i4)
c        write(6,67) '(4ab  int. ++ =',buf2(5,ijkl,ij_pi_pj,kl,iqu)
c        write(6,67) '(4ab  int. +  =',buf2(5,ijkl,ij_pi   ,kl,iqu)
c        write(6,67) '(2b   int. -+ =',buf2(3,ijkl,ij_mi_pj,kl,iqu)
c        write(6,67) '(2b   int. -  =',buf2(3,ijkl,ij_mi   ,kl,iqu)
c        write(6,67) '  final deriv =',der
c 67     format(a15,f12.7)
c
c
      der2(nder_ab,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_ac(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_ac,der2)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
ctest
c     write(6,60) icart,jcart,nder_ac
c 60  format('from AC: icart,jcart=',2i3,' no=',i3)
c
c
      do 200 kl=nfu(nqkl)+1,lnkl
      kl_pj=npxyz(jcart,kl)
      n_kl_j=nia(jcart,kl)
      kl_mj=0
      if(n_kl_j.gt.0) kl_mj=nmxyz(jcart,kl)
c
      do 200 ij=nfu(nqij)+1,lnij
      ij_pi=npxyz(icart,ij)
      n_ij_i=nia(icart,ij)
      ij_mi=0
      if(n_ij_i.gt.0) ij_mi=nmxyz(icart,ij)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_ac_pipj=buf2(6,ijkl,ij_pi,kl_pj,iqu)
c
c       two_c_mi_pj=buf2(4,ijkl,ij_mi,kl_pj,iqu)
c       two_a_pi_mj=buf2(2,ijkl,ij_pi,kl_mj,iqu)
c       two_0_mi_mj=buf2(1,ijkl,ij_mi,kl_mj,iqu)
c
      der=buf2(ijkl,ij_pi,kl_pj,iqu,6)
      if(n_ij_i.gt.0 .and. ij_mi.gt.0) then
        der=der-n_ij_i*buf2(ijkl,ij_mi,kl_pj,iqu,4)
      endif
      if(n_kl_j.gt.0 .and. kl_mj.gt.0) then
        der=der-n_kl_j*buf2(ijkl,ij_pi,kl_mj,iqu,2)
      endif
      if(n_ij_i.gt.0 .and. n_kl_j.gt.0) then
        der=der+n_ij_i*n_kl_j*buf2(ijkl,ij_mi,kl_mj,iqu,1)
      endif
c
      der2(nder_ac,ijkl,ij,kl,iqu)=der
c
c        write(6,66) nder_ac, ij,kl
c        write(6,*)'ij   ,kl   =',ij   ,kl
c        write(6,*)'ij_pi,kl_pj=',ij_pi,kl_pj
c        write(6,*)'ij_mi,kl_pj=',ij_mi,kl_pj
c        write(6,*)'ij_pi,kl_mj=',ij_pi,kl_mj
c        write(6,*)'ij_mi,kl_mj=',ij_mi,kl_mj
c  66    format('no=',i3,'ij,kl=',2i4)
c        write(6,67) '(4ac  int. ++ =',buf2(6,ijkl,ij_pi,kl_pj,iqu)
c        write(6,67) '(2c   int. -+ =',buf2(4,ijkl,ij_mi,kl_pj,iqu)
c        write(6,67) '(2a   int. +- =',buf2(2,ijkl,ij_pi,kl_mj,iqu)
c        write(6,67) '(     int. -- =',buf2(1,ijkl,ij_mi,kl_mj,iqu)
c        write(6,67) '  final deriv =',der
c 67     format(a15,f12.7)
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_bc(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_bc,der2,xab)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
c2001 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      kl_pj=npxyz(jcart,kl)
      n_kl_j=nia(jcart,kl)
      kl_mj=0
      if(n_kl_j.gt.0) kl_mj=nmxyz(jcart,kl)
c
      do 200 ij=nfu(nqij)+1,lnij
      ij_pi=npxyz(icart,ij)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_bc_pipj=buf2(7,ijkl,ij_pi,kl_pj,iqu)
c       four_bc_0ipj=buf2(7,ijkl,ij   ,kl_pj,iqu)
c
c       two_b_pi_mj=buf2(3,ijkl,ij_pi,kl_mj,iqu)
c       two_b_0i_mj=buf2(3,ijkl,ij   ,kl_mj,iqu)
c
c
      der=                buf2(ijkl,ij_pi,kl_pj,iqu,7)
     *  + xab(ijkl,icart)*buf2(ijkl,ij   ,kl_pj,iqu,7)
c
      if(n_kl_j.gt.0 .and. kl_mj.gt.0) then
          der=der
     *    - n_kl_j*(                  buf2(ijkl,ij_pi,kl_mj,iqu,3)
     *              + xab(ijkl,icart)*buf2(ijkl,ij   ,kl_mj,iqu,3) )
     *
      endif
c
      der2(nder_bc,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine resc_fder(ncf,na,fock,thres)
      implicit real*8 (a-h,o-z)
      dimension fock(3,na,*)       ! ntri
      data half /0.5d0/
c------------------------------------------------------
c divide the non-diagonal elements of the Fder matrix by 2
c and scale back to the normal values from the large ones
c-------------------------------------------------------
      thinv2=half*thres
      ij=0
      do icf=1,ncf
         do jcf=1,icf
            ij=ij+1
            do iat=1,na
               do icr=1,3
                 fock(icr,iat,ij)=fock(icr,iat,ij)*thinv2
      if(icf.eq.jcf) fock(icr,iat,ij)=fock(icr,iat,ij)+fock(icr,iat,ij)
               enddo
            enddo
         enddo
      enddo
c-------------------------------------------------------
      end
c=================================================================
      subroutine resu_fder(ncf,na,fockA,fockB,thres)
      implicit real*8 (a-h,o-z)
      dimension fockA(3,na,*)       ! ntri
      dimension fockB(3,na,*)       ! ntri
c------------------------------------------------------
c for open shell need to double all Fock matrix elements
c so, compared to closed-shell, need to double only the
c diagonal elements (see <resc_fder> above)
c-------------------------------------------------------
      ij=0
      do icf=1,ncf
         do jcf=1,icf
            ij=ij+1
            do iat=1,na
               fockA(1,iat,ij)=fockA(1,iat,ij)*thres
               fockA(2,iat,ij)=fockA(2,iat,ij)*thres
               fockA(3,iat,ij)=fockA(3,iat,ij)*thres
               fockB(1,iat,ij)=fockB(1,iat,ij)*thres
               fockB(2,iat,ij)=fockB(2,iat,ij)*thres
               fockB(3,iat,ij)=fockB(3,iat,ij)*thres
            enddo
         enddo
         do iat=1,na
            fockA(1,iat,ij)=fockA(1,iat,ij)+fockA(1,iat,ij)
            fockA(2,iat,ij)=fockA(2,iat,ij)+fockA(2,iat,ij)
            fockA(3,iat,ij)=fockA(3,iat,ij)+fockA(3,iat,ij)
            fockB(1,iat,ij)=fockB(1,iat,ij)+fockB(1,iat,ij)
            fockB(2,iat,ij)=fockB(2,iat,ij)+fockB(2,iat,ij)
            fockB(3,iat,ij)=fockB(3,iat,ij)+fockB(3,iat,ij)
         enddo
      enddo
c-------------------------------------------------------
      end
c=====================================================================
      subroutine fdersymm_hes(ngener,ncf,na,fock,nsyo,ifp,nupair)
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------
c It symmetrizes Fock derivative matrices (hessian) for all atoms.
c These matrices are dimensioned as f(3,natoms,ntri).
c These matrices are defined as G(D0,g1) for x,y,z of all atoms
c-----------------------------------------------------------------
c This routine builds the whole Fder matrix from the skeleton one
c PARAMETERS:
c INPUT:
c    ngener= number of group generators
c    ncf=number of contracted basis functions
c    nfock=number of Fock matrices
c    fock=Fock matrices
c    ifp=function pairs, ifp(isymop,k) gives the symmetry
c       parter of basis function k under symm. op. isymop;
c       if ifp(isymop,k) is negative then the function changes
c       sign upon reflection or rotation
c OUTPUT: the full Fock matrices in fock
c-----------------------------------------------------------------
      common /negxyz/ negx(7),negy(7),negz(7)
      dimension nsyo(7)
      dimension mirror(3,7)
      dimension fock(3,na,*)           ! ntri
      dimension nupair(na,*)           ! (natoms,nsym)
      dimension ifp(7,ncf)
      dimension ngxyz(3)
      data one,half /1.d0 , 0.5d0/
      data mirror/-1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c-----------------------------------------
      do 10 ns=1,ngener
         nop=nsyo(ns)
         negx(ns)=mirror(1,nop)
         negy(ns)=mirror(2,nop)
         negz(ns)=mirror(3,nop)
   10 continue
c-----------------------------------------
c
      do 350 ns=1,ngener
         ngxyz(1)=negx(ns)
         ngxyz(2)=negy(ns)
         ngxyz(3)=negz(ns)
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if(i1.lt.0) then
               i1=-i1
               fct=-fct
            endif
            j1=ifp(ns,j)
            if(j1.lt.0) then
               j1=-j1
               fct=-fct
            endif
            ij1=i1*(i1-1)/2+j1
            if(j1.gt.i1) ij1=j1*(j1-1)/2+i1
                do iat=1,na
                   iat1=nupair(iat,ns)
                if(iat.eq.iat1) then
                   do icr=1,3
                     sign=fct*ngxyz(icr)
                     ff=fock(icr,iat,ij)+sign*fock(icr,iat1,ij1)
                     if(ij.gt.ij1) ff=ff*half
                     fock(icr,iat,ij)=ff
                     fock(icr,iat1,ij1)=sign*ff
                   enddo
                endif
                if(iat.ne.iat1) then
                   do icr=1,3
                     sign=fct*ngxyz(icr)
                     ff=fock(icr,iat,ij)+sign*fock(icr,iat1,ij1)
                     if(iat.gt.iat1) ff=ff*half
                     fock(icr,iat,ij)=ff
                     fock(icr,iat1,ij1)=sign*ff
                   enddo
                endif
                enddo
  340    continue
  350 continue
c
      end
c=====================================================================
      subroutine fder4nat(idft,rhf,natoms,ntri,ntimes,natonce)

      use memory

      logical rhf
      data ncall /0/
      save ncall
c-------------------------------------------------------------
      ncall=ncall+1
c-------------------------------------------------------------
ctest
c     natonce=6
c     ntimes=2
c     RETURN
c-------------------------------------------------------------
      call getival('iout',iout)
      call getival('lcore',lcore)
      call getmem(0,last0)
      call retmem(1)
c-------------------------------------------------------------
c the following allocations will be needed :
c
c     nmemory=ntri*3*natonce    for F1=(3,natonce,ntri)=G(d0,g1)
c
c     or twice as much if UHF
c
c     and one "3*ntri" for trsp_inplace
c-------------------------------------------------------------
c To do F1 for NATONCE atoms at once we need :
c
c   needed=(natonce*ntri*3)+ntri*3=(natonce+1)*3*ntri
c       or if UHF
c   needed=(2*natonce*ntri*3)+ntri*3=(2*ntonce+1)*3*ntri
c
c and also 6*ntri in the case of symmetry
c-------------------------------------------------------------
c     if(rhf) then
c        needed=(natonce+1)*ntri*3
c     else
c        needed=(2*natonce+1)*ntri*3
c     endif
c
c     call getival('nsym',nsym)
c     if(nsym.gt.0 .and.idft.gt.0) needed=needed+6*ntri
c-------------------------------------------------------------
c  memory we have :
c
      nmemory=lcore-last0
c leave some memory for integrals :
      nmemory=nmemory - 1000000
c-------------------------------------------------------------
c needed for one atom :
c
      if(rhf) then
         memory1=6*ntri
      else
         memory1=9*ntri
      endif
c
      call getival('nsym',nsym)
c
      if(nsym.gt.0 .and. idft.EQ.0) memory1=memory1+6*ntri
c-------------------------------------------------------------
      if(memory1.gt.nmemory) then
         write(iout,*)' not enough memory to handle one atom in Fder'
         write(iout,*)' needed =',memory1,' available =',nmemory
         write(iout,*)' increase declared memory by  =',memory1-nmemory
         call para_stop
         stop ' stopped because of memory problem in hess2e1'
      endif
c-------------------------------------------------------------
c calculate Natonce :    memory .GE. needed
c
      natonce=nmemory/(3*ntri)-1
c
      if(nsym.gt.0 .and. idft.EQ.0) then
         natonce=(nmemory-6*ntri)/(3*ntri)-1
      endif
c
      if(.not.rhf) natonce=natonce/2
c-------------------------------------------------------------
ctest
c       write(6,*)' memory=',nmemory,' atoms at once=',natonce
ctest
c
      if(natonce.lt.0) natonce=0
      if(natonce.gt.natoms) natonce=natoms
c-------------------------------------------------------------
      if(natonce.eq.0) then
         write(iout,*)' not enough memory to handle one atom in Fder'
         write(iout,*)' needed =',memory1,' available =',nmemory
         write(iout,*)' increase declared memory by  =',memory1-nmemory
         call para_stop
         stop ' stopped because of memory problem in hess2e1'
      else
         ntimes=natoms/natonce
         moreat=natoms-ntimes*natonce
         if(moreat.ne.0) ntimes=ntimes+1
c2003 to make more equal pieces :
         natonce=natoms/ntimes
         if(natoms.gt.ntimes*natonce) natonce=natonce+1
c2003
c
      endif
c-------------------------------------------------------------
      end
c======================================================================
      subroutine make_natodo(natoms,natonce,ntimes,natodo)
      dimension natodo(2,ntimes)
c
      if(ntimes.eq.1) then
        natodo(1,1)=1
        natodo(2,1)=natoms
        return
      endif
c
      nbeg=0
      do itime=1,ntimes-1
        natodo(1,itime)=nbeg+1
        natodo(2,itime)=natonce+nbeg
        nbeg=nbeg+natonce
      enddo
c
      natodo(1,ntimes)=nbeg+1
      natodo(2,ntimes)=natoms
c
      end
c======================================================================
      subroutine trsp_ondsk(iunit,ntri,natb,nate,f1,bl)
      implicit real*8 (a-h,o-z)
      dimension bl(ntri,3)
      dimension f1(3,natb:nate,ntri)
c
      do iat=natb,nate
         do ij=1,ntri
            bl(ij,1)=f1(1,iat,ij)
            bl(ij,2)=f1(2,iat,ij)
            bl(ij,3)=f1(3,iat,ij)
         enddo
ccccccc  call save1mat(iunit,iat,ntri*3,bl)
         write(unit=iunit,rec=iat) bl
      enddo
c
      end
c======================================================================
      subroutine fdersymm_hes1(iunit,ngener,ncf,ntri,na,fock,
     *                         nsyo,ifp,nupair)
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------
c It symmetrizes Fock derivative matrices (hessian) for all atoms.
c These matrices are dimensioned as f(ntri,3,natoms).
c These matrices are defined as G(D0,g1) for x,y,z of all atoms
c-----------------------------------------------------------------
c This routine builds the whole Fder matrix from the skeleton one
c PARAMETERS:
c INPUT:
c    ngener= number of group generators
c    ncf=number of contracted basis functions
c    fock=Fock matrices
c    ifp=function pairs, ifp(isymop,k) gives the symmetry
c       parter of basis function k under symm. op. isymop;
c       if ifp(isymop,k) is negative then the function changes
c       sign upon reflection or rotation
c OUTPUT: the full Fock matrices in fock
c-----------------------------------------------------------------
      common /negxyz/ negx(7),negy(7),negz(7)
      dimension nsyo(7)
      dimension mirror(3,7)
      dimension fock(ntri,3,2)
      dimension nupair(na,*)           ! (natoms,nsym)
      dimension ifp(7,ncf)
      dimension ngxyz(3)
      data one,half /1.d0 , 0.5d0/
      data mirror/-1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c-----------------------------------------
      do 10 ns=1,ngener
         nop=nsyo(ns)
         negx(ns)=mirror(1,nop)
         negy(ns)=mirror(2,nop)
         negz(ns)=mirror(3,nop)
   10 continue
c-----------------------------------------
c
      do 350 ns=1,ngener
         ngxyz(1)=negx(ns)
         ngxyz(2)=negy(ns)
         ngxyz(3)=negz(ns)
c
      do 345 iat=1,na
         iat1=nupair(iat,ns)
         call read1mat(iunit,iat ,ntri*3,fock(1,1,1))
       if(iat1.ne.iat) then
         call read1mat(iunit,iat1,ntri*3,fock(1,1,2))
       endif
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if(i1.lt.0) then
               i1=-i1
               fct=-fct
            endif
            j1=ifp(ns,j)
            if(j1.lt.0) then
               j1=-j1
               fct=-fct
            endif
            ij1=i1*(i1-1)/2+j1
            if(j1.gt.i1) ij1=j1*(j1-1)/2+i1
                if(iat.eq.iat1) then
                   do icr=1,3
                     sign=fct*ngxyz(icr)
                     ff=fock(ij,icr,1)+sign*fock(ij1,icr,1)
                     if(ij.gt.ij1) ff=ff*half
                     fock(ij ,icr,1)=ff
                     fock(ij1,icr,1)=sign*ff
                   enddo
                endif
                if(iat.ne.iat1) then
                   do icr=1,3
                     sign=fct*ngxyz(icr)
                     ff=fock(ij,icr,1)+sign*fock(ij1,icr,2)
                     if(iat.gt.iat1) ff=ff*half
                     fock(ij ,icr,1)=ff
                     fock(ij1,icr,2)=sign*ff
                   enddo
                endif
  340    continue
         call save1mat(iunit,iat ,ntri*3,fock(1,1,1))
        if(iat1.ne.iat) then
         call save1mat(iunit,iat1,ntri*3,fock(1,1,2))
        endif
  345    continue
  350 continue
c
      end
c=====================================================================
      subroutine trsp_ondsx(iunit,ntri,natb,nate,f1,bl,listunq)
      implicit real*8 (a-h,o-z)
      dimension listunq(*)
      dimension bl(ntri,3)
      dimension f1(3,natb:nate,ntri)
c
      do iat=natb,nate
         iatr=listunq(iat)   !   real atom
         do ij=1,ntri
            bl(ij,1)=f1(1,iat,ij)
            bl(ij,2)=f1(2,iat,ij)
            bl(ij,3)=f1(3,iat,ij)
         enddo
ccccccc  call save1mat(iunit,iatr,ntri*3,bl)
         write(unit=iunit,rec=iatr) bl
      enddo
c
      end
c======================================================================
      subroutine fdersymm_hes2(iunit,lrec,nsym  ,ncf,ntri,na,fock,
     *                         nsyo,ifp,nupair,natunq,listunq)
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------
c It symmetrizes Fock derivative matrices (hessian) for all atoms.
c These matrices are dimensioned as f(ntri,3,natoms).
c These matrices are defined as G(D0,g1) for x,y,z of all atoms
c-----------------------------------------------------------------
c This routine builds Fder matrices for all atoms
c from these for symmetry unique ones :
c PARAMETERS:
c INPUT:
c    ngener= number of group generators
c    ncf=number of contracted basis functions
c    fock=Fock matrices
c    ifp=function pairs, ifp(isymop,k) gives the symmetry
c       parter of basis function k under symm. op. isymop;
c       if ifp(isymop,k) is negative then the function changes
c       sign upon reflection or rotation
c OUTPUT: the full Fock matrices in fock
c-----------------------------------------------------------------
      common /negxyz/ negx(7),negy(7),negz(7)
      dimension nsyo(7)
      dimension mirror(3,7)
      dimension fock(ntri,3,2)
      dimension nupair(na,*)           ! (natoms,nsym)
      dimension ifp(7,ncf)
      dimension listunq(*)
      dimension ngxyz(3)
      data one,half /1.d0 , 0.5d0/
      data mirror/-1,1,1,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c-----------------------------------------
      do 10 ns=1,nsym
         nop=nsyo(ns)
         negx(ns)=mirror(1,nop)
         negy(ns)=mirror(2,nop)
         negz(ns)=mirror(3,nop)
   10 continue
c-----------------------------------------
      do 350 ns=1,nsym
         ngxyz(1)=negx(ns)
         ngxyz(2)=negy(ns)
         ngxyz(3)=negz(ns)
c
CCCC  do 345 iat=1,na
      do 345 iatu=1,natunq
         iat=listunq(iatu)
         iat1=nupair(iat,ns)
         if(iat.EQ.iat1) go to 345
         call read1mat(iunit,lrec+iat ,ntri*3,fock(1,1,1))
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if(i1.lt.0) then
               i1=-i1
               fct=-fct
            endif
            j1=ifp(ns,j)
            if(j1.lt.0) then
               j1=-j1
               fct=-fct
            endif
            ij1=i1*(i1-1)/2+j1
            if(j1.gt.i1) ij1=j1*(j1-1)/2+i1
               do icr=1,3
                  fock(ij1,icr,2)=fock(ij,icr,1)*fct*ngxyz(icr)
               enddo
  340    continue
         call save1mat(iunit,lrec+iat1,ntri*3,fock(1,1,2))
  345 continue
  350 continue
c
      end
c======================================================================
