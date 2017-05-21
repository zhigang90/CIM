c=============================================================
c
      subroutine para_giao(idft,ax,nblocks,bl,inx,ntri,thres,
     *                     dscreen,dens,fock,labels)
c
c=============================================================
c This routine calls two-el. int. block by block
c and constructs closed-shell Fock matrix .
c---------------------------------------------------------------------

      use newpara

      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      COMMON /DFTCoeff/ aXX(18) ! contains 15 density functional coefficients
      dimension inx(12,*)
      dimension bl(*),dscreen(*),dens(*),fock(*)
      dimension labels(*)
c
      call getival('gran',igran)
      call getival('ncs ',ncs)
c
c     if(scftype.ne.'full-direct') then
c        call nerror(53,'para_grad',
c    *                  ' can only run full-direct mode',0,0)
c     endif
c
      ntri3=3*ntri
      ncs2=2*ncs*ncs
      call para_initsend              
      call para_pack_real(dens,ntri)
      call para_pack_real(dscreen,ncs2)
      call para_pack_int(idft,1)
      call para_pack_real(aXX,18)
      call para_bcast_pack(TxPostDens)
c
c Give work assignments to slaves until all work is done
c      
      call para_distr_blocks(nblocks,igran)
c
      call para_reduce(fock,ntri3,TxPostFock)
      end
c=====================================================================
c
      subroutine para_cphf(nblocks,bl,inx,ntri,thres,dscreen,dens,
     *                     fock,labels,myjob,igran,mgo)
c
c=====================================================================

      use newpara

      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension bl(*),dscreen(*),dens(ntri,3),fock(ntri,3)
      dimension labels(*)
c
      call getival('ncs ',ncs)
      ntri3=3*ntri
      ncs2=2*ncs*ncs
      call para_initsend
      call para_pack_real(dens,ntri3)
      call para_pack_real(dscreen,ncs2)
      call para_pack_int(mgo,1)
      call para_bcast_pack(TxPostDens)
c
c Give block assignments to slaves until all work is done
c
      call para_distr_blocks(nblocks,igran)
c      
      call para_reduce(fock,ntri3,TxPostFock)
      end
c=======================================================================
      subroutine para_shift_start(bl,lden,ldn1,ntri)
c set up shielding tensor calculations
c send over density and disturbed density matrices

      use newpara

      implicit real*8 (a-h,o-z)
      dimension bl(*)
      mgo=0
      call para_bcast(mgo,TxNext) !mgo=0
      call para_bcast_real(bl(lden),ntri,TxNext) !D00
      call para_bcast_real(bl(ldn1),3*ntri,TxNext) !D10
      end
c=======================================================================
      subroutine para_shi(na,bl,inx,last,iprint,iout,datbas,datnuc,
     *                ncs,ncf,ntri,nreq,den,dn1,shie,hna,nra,nrat)
c drive calculations of shielding tensors
c get request for work (and result if any)
c and send job

      use newpara

      implicit real*8 (a-h,o-z)
      dimension bl(1)
      dimension shie(3,3,*)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,1)
c
      call para_recv_pack(islave,TxShReq)
      call para_unpack_int(isres,1)
      if(isres.eq.1) then
         call para_unpack_int(nraold,1)
         call para_unpack_real(shie(1,1,nraold),9) !dia
         call para_unpack_real(shie(1,1,nraold+nreq),9) !para
      endif
      call para_initsend
      call para_pack_int(nrat,1)
      call para_pack_int(nra,1)
      call para_send_pack(islave,TxShAssign)
      end
c=======================================================================
      subroutine para_shift_end(shie,nreq)
c go over slaves, get work requests (and results if any)
c tell them to finish up

      use newpara

      implicit real*8 (a-h,o-z)
      dimension bl(1)
      dimension shie(3,3,*)
      do i=1,nslv
         call para_recv_pack(islave,TxShReq)
         call para_unpack_int(isres,1)
         if(isres.eq.1) then
            call para_unpack_int(nraold,1)
            call para_unpack_real(shie(1,1,nraold),9) !dia
            call para_unpack_real(shie(1,1,nraold+nreq),9) !para
         endif
         call para_initsend
         call para_pack_int(0,1) !nrat=0, we are over
         call para_send_pack(islave,TxShAssign)
      end do
      end
