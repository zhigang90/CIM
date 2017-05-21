c=======================================================================
c            +-------------------------------------+
c            +                                     +
c            +      two-electron integral          +
c            +       for polarizability            +
c            +   with field-independent basis set  +
c            +                                     +
c            +         Krzysztof Wolinski          +
c            +             May  , 2001             +
c            +-------------------------------------+
c=======================================================================
      subroutine polar2(idft,ax,rhf,bl,inx)

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
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c------------------------------------------------
c Hessian Gradient Polarizability  options :
c
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension xintxx(9)
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
c prepare the integral system for          integrals (iforwhat=4) :
c (like twoint95 for ordinary integrals)
c Full list of parameters for twoint95 (13) is reduced here to only
c FOUR . The remaning parameters are taken from commons set up for
c the SCF . Hopefully, it will be always the case for post-scf calc.
c I want to use this interface for all post-scf applications (nmr,grad
c and second derivatives). Four input parameters will be taken from
c appropriate commons set up when a given application is called. Here
c four input parameters taken from HESS option commons :
c
c
      call getrval('thr2',thresh)
      iforwhat=1       !   polarizability integrals (& derivatives)
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
c WARNING:
c MM 08/18/06
c as I understand it, currently the polarizability calcultion
c is not parallel. All the calculation is done by the Master!!!
c So all the call to parallel routines have the only effect
c of screwing up whatever parallel step is executed after!!!
c I am commenting them out.
c
c     call para_post_scf(ncache,iforwhat,nblocks)
c
c output (see scf) : last line parameters
c----------------------------------------------------------------------
      call getival('nopola',nopola)
      if(nopola.ne.0) return
c----------------------------------------------------------------------
c
c if filed-dependent basis set is used the Calculate
c
c (1) the first derivative 2-el. integrals and form a part
c     of the Ist-order Fock matreices G(D,g1)
c
c     will never be re-calculated !!!
c----------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c----------------------------------------------------------------------
c allocate memory for label's arrays & for G(D0,g1)
c
      call getmem(maxlabels,labels)
      call getmem(ntri*3,lfd0g1)
c
      call zeroit(bl(lfd0g1),ntri*3)
c----------------------------------------------------------------------
c get addresses for density , fock-der matrices and hessian matrix
c (memory allocated in hessana.f)
c
      call getival('ldensi',lden  )
      call getival('lfxyz',lfxyz)
      If(.not.rhf) call getival('ldensB',mden)
c
c dimension of Fock          is (ntri,3)
c----------------------------------------------------------------------
      call secund(tgrad1)
      call elapsec(elagrad1)
c----------------------------------------------------------------------
c calculate two-electron integral derivatives for hessian
c (just like other int_name calls)
c
CN    call para_pola(idft,ax,rhf,nblocks,bl,inx,ntri,thresh,
CN   *               bl(lden),bl(mden),bl(lfd0g1),bl(labels))
c----------------------------------------------------------------------
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) total/60.0d0,elaps/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for POLA 2e integrals = '
     *,f8.2,' Elapsed = ',f8.2,' min')
  500 format(58('-'))
c----------------------------------------------------------------------
c re-scale these parts of fock matrices :
c
        call dscal(ntri*3,thresh,bl(lfd0g1),1)
c----------------------------------------------------------------------
c calculates the POLA exchange matrix derivatives
c
      if(idft.ne.0) then
         call secund(tgiao2)
         call elapsec(etgi2)
c
         call getival('ndum',ndum)
CN       call dft_pola(idft,lfxyz,na-ndum,ncf,ibas,iprint)
c----------------------------------------------------------------------
         call secund(tgiao3)
         call elapsec(etgi3)
         total=(tgiao3-tgiao2)/60.0d0
         elaps=(etgi3-etgi2)/60.0d0
c        write(iout,500)
         write(iout,600) total,elaps
         call f_lush(iout)
c        write(iout,500)
  600 format('Master CPU time for POLA XC integrals = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
      endif
c----------------------------------------------------------------------
c test
c     write(8,*)'From POLAR2 : f(d0,g1)  '
c     lfx=lfxyz
c     lfy=lfx+ntri
c     lfz=lfy+ntri
c     call drumx (bl(lfx),ncf,8   ,'fock   X' )
c     call drumx (bl(lfy),ncf,8   ,'fock   Y' )
c     call drumx (bl(lfz),ncf,8   ,'fock   Z' )
c----------------------------------------------------------------------
c Symetrize these parts of the Fock matrices for exac symmetry :
c
      call getival('nsym',nsym)
      if(nsym.gt.0) then
         call getival('ngener',ngener)
         call getival('nsyo',nsyo)
         call getival('SymFunPr',ifp)
         call focksymm_nmr(ngener,ncf,bl(lfd0g1),bl(ifp),bl(nsyo))
      endif
c----------------------------------------------------------------------
c and add them to the Fock (lfxyz)
c
      call add(bl(lfd0g1),bl(lfxyz),ntri*3)
c----------------------------------------------------------------------
c
ctest
c     write(8,*)'From Shift2 : Full focks'
c     call drumx (bl(lfxyz+1),ncf,6   ,'fock   X' )
c     call drumx (bl(lfxyz+ntri),ncf,6   ,'fock   Y' )
c     call drumx (bl(lfxyz+ntri*2),ncf,6   ,'fock   Z' )
c----------------------------------------------------------------------
c release memory reserved for labels and fd0gX,Y,Z
c
      call retmem(2)
c----------------------------------------------------------------------
c  Return the SCF values of eps,eps1,epsr
c
      eps=reps
      eps1=reps1
      epsr=repsr
      icache=icacher
c----------------------------------------------------------------------
c
      end
c===================================================================
