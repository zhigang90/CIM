c =====================================================================
c
      subroutine get_storage(iwher)
c
c =====================================================================
c get the number of slave ready to store integrals
c put the number of already stored integrals and quartets in common datacore

      use newpara

      implicit real*8 (a-h,o-z)
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c
c
      integer*8 xntegrals,integ62,integ72,integx2,integx4,integx8
      common /datadisk/xntegrals,integ62,integ72,integx2,integx4,integx8
c
      call getival('indisk',indisk)
c
c check if the slaves are OK
c
      call para_check
c get the next available slave
      call para_recv_pack(iwher,TxStorReq)
      call para_unpack_int(integrals,1)
      call para_unpack_int(nqstore,1)
      call para_unpack_int(nblstore,1)
c
      if(indisk.gt.0) then
         call para_unpack_real(xntegrals,1)
         call para_unpack_real(integ62,1)
         call para_unpack_real(integ72,1)
         call para_unpack_real(integx2,1)
         call para_unpack_real(integx4,1)
         call para_unpack_real(integx8,1)
c        xntegrals=rntegrals
c        integ62=rnteg62
c        integ72=rnteg72
c        integx2=rntegx2
c        integx4=rntegx4
c        integx8=rntegx8
      endif
c
c     write(6,*)'     '
c     write(6,*)'PVM_STOR: IWHER=',iwher,' nblstored=',nblstore
c     write(6,*)'      xntegrals=',xntegrals
c     write(6,*)'      xntegra72=',integ72
c     write(6,*)'      xntegra62=',integ62
c     write(6,*)'      xntegrag2=',integx2
c     write(6,*)'      xntegrag4=',integx4
c     write(6,*)'      xntegrag8=',integx8
c
      end
c =====================================================================
c
      subroutine para_store(ifrom,ito,bl,inx,dens,thres,labels,
     *                      iwher,isecond)
c
c =====================================================================
c send the integral batch to store (superblocks between ifrom and ito)
c to the next slave, switch for doing positively priced blocks also sent

      use newpara

      implicit real*8 (a-h,o-z)
c     common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
c    *                          isiz,jsiz,ksiz,lsiz, nqstore
      dimension inx(12,*)
      dimension bl(*),dens(*)
c
      call getival('intpack',ipack)
      call getival('qrtpack',iqpac)
      call getival('blkpack',ibpac)
      call getival('maxsint',maxsi)
c
c       write(*,*)'pvm_stor packs=',ipack,iqpac,ibpac,' iwher=',iwher
c
c     call getival('indisk',indisk)
c
c     write(6,*)' pvm_store indisk=',indisk
c
c check if the slaves are OK
c
      call para_check
c send the message
      call para_initsend
      call para_pack_int(ifrom,1)
      call para_pack_int(ito,1)
      call para_pack_int(isecond,1)
ckw
      call para_pack_int(ipack,1)
      call para_pack_int(iqpac,1)
      call para_pack_int(ibpac,1)
      call para_pack_int(maxsi,1)
ckw
      call para_send_pack(iwher,TxStorAssign)
c
c     write(*,*)'pvm_stor after sent packs'
c
      end
c =====================================================================
c
      subroutine send_limits(isecond,ito)
c
c =====================================================================
c send out two parameters so that all slaves have the same info
c isecond shows if we have stored posively priced superblocks
c ito is the number of the last superblock stored

      use newpara

      implicit real*8 (a-h,o-z)
c
c check if the slaves are OK
c
      call para_check
c send the info to every slave
      do iproc=1,nslv
         call para_recv(imsg,islave,TxStore)
         call para_initsend
         call para_pack_int(ito,1)
         call para_pack_int(isecond,1)
         call para_send_pack(islave,TxStore)
      end do
      end
