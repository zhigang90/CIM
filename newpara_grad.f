c=============================================================
c
      subroutine para_post_scf(ncachex,iforwhatx,nblocks)
c
c=============================================================

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      Logical rhf
c------------------------------------------------
c NMR options :
c
      common /nmrint/ ichf,ncache,iprint,nogiao,noacce,ngauge,ntrans
c------------------------------------------------
c
c -- get data from the depository
      call getival('NAlpha',NAlpha)
      call getival('NBeta',NBeta)
      call getival('Multip',IMult)
c
      rhf = NBeta.EQ.0.AND.IMult.EQ.1
c
c Obtain and pack slave initialization data
c
      call para_initsend
c
c parameters for 2-electron integral derivatives
c
      call para_pack_int(ncachex,1)
      call para_pack_int(iforwhatx,1)
      call para_pack_int(NAlpha,1)
      call para_pack_int(NBeta,1)
      call para_pack_int(rhf,1)
c
      call getrval('thr1',thre1)
      call getrval('thr2',thre2)
c
      call para_pack_real(thre1,1)
      call para_pack_real(thre2,1)
c2002
      call getival('iroute',iroute)
      call getival('stab',istab)
c
      call para_pack_int(iroute,1)
      call para_pack_int(istab,1)
c2002
c2006
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc1)
      if(iftc1.NE.0) then
         call getival('ccdd',iccdd)
         call para_pack_int(iccdd,1)
      endif
c2006
c
c send data needed for NMR
c
      if (iforwhatx.eq.2) then
         call getrval('thr3',thre3)
         call para_pack_real(thre3,1)
         call para_pack_int(ichf,1)
         call para_pack_int(nogiao,1)
         call para_pack_int(noacce,1)
         call para_pack_int(ngauge,1)
         call para_pack_int(ntrans,1)
      end if
c
      call para_bcast_pack(TxPostInit)
c
      end
c=============================================================
c
      subroutine para_grad(idft,ax,rhf,nblocks,bl,inx,ntri,
     *                     thres,dens,denB,atforces,labels)

      use newpara

      implicit real*8 (a-h,o-z)
c
c=============================================================
c This routine calculates derivative two-el.integrals (for gradient)
c block by block and then calculates "2-el." forces on atoms .
c Actually, this routine handles the PVM calls and distributes tasks
c among nodes. If there is only one cpu then int_grad() is called.
c---------------------------------------------------------------------
c integral derivatives are calculated in int_grad()
c---------------------------------------------------------------------
c INPUT :
c  idft              - dft flag:  0 - no dft
c  ax                - factor to multiply exchange contribution for ACM
c  nblocks           - number of blocks of contracted shell quartets
c  bl                - scratch for everything
c  inx(12,ncs)       - basis set & nuclear data
c  natoms3           - naatom*3 , size of atforces array
c  thres             - integral treshould
c  dens()            - ordinary density matrix
c  labels            - labels buffer
C OUTPUT:
c atforces(natoms,3) - forces on atoms (two-electron part)
c---------------------------------------------------------------------
      character*11 scftype
      character*4 where
      Logical rhf
      common /runtype/ scftype,where
      COMMON /DFTCoeff/ aXX(18) ! contains 18 density functional coefficients
      dimension inx(12,*)
      dimension bl(*),dens(*),denB(*),atforces(*)
      dimension labels(*)
c
      call getival('gran',igran)
      call getival('na  ',na)
      nat3=na*3
c
      if(scftype.ne.'full-direct') then
         call nerror(53,'para_grad',
     *                  ' can only run full-direct mode',0,0)
      endif
c
c -- first get density matrices and other data
      call para_initsend              
      call para_pack_real(dens,ntri)
      If(.NOT.rhf) call para_pack_real(denB,ntri)
      call para_pack_int(idft,1)
      call para_pack_real(aXX,18)
      call para_bcast_pack(TxPostDens)
c
c *************************************************
c  one-electron integrals
c *************************************************
c
c -- first we are going to do the electron-nuclear attraction
c    1-e integrals that were skipped during the serial evaluation
c    These will be parallelized over atoms
c
      call elapsec(t1)
c
      If(nslv.GT.na) Then
c -- we are going to temporarily reduce the number of active slaves
        call message('**WARNING** in Force Module',
     $   'Temporarily reducing number of active slaves to ',na,0)
      EndIf
c
      DO IAtm=1,na
      call para_recv(imsg,islave,TxDftReq)
      call para_send(IAtm,islave,TxDftAssign)
      EndDO
c
c -- send termination to each slave
      DO IAtm=1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_send(0,islave,TxDftAssign)
      EndDO
c
      call elapsec(t2)
      t1 = (t2-t1)/60.0d0
cc      write(6,*) ' Time for 1-e e-nuclear forces is:',t1,' min'
c
c *************************************************
c  two-electron integrals proper
c *************************************************
c
c
c Give block assignments to slaves until all work is done
c
      call para_distr_blocks(nblocks,igran)      
c
      call para_reduce(atforces,nat3,TxPostFock)
c
c tell the slave that we are finished (argument=0)
c
      if(idft.eq.0) call para_next(0)
      end
