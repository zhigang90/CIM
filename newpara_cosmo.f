c======================================================================
      subroutine para_cosmo_data(cosurf,xyz,natom,nps,npspher)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this subroutine send the COSMO data to the slaves
c
c  this routine is executed by the master
c
c   cosurf        cartesian coordinates of surface segments
c   xyz           atomic coordinates (real atoms only)
c   natom         number of real atoms
c   nps           number of surface segments
c   npspher       number of outer surface segments
c
      dimension cosurf(*),xyz(3,natom)
      character*256 conel
c
c   gather some parameters that need to be sent
c
      call getchval('c_onel',conel)
      call getrval('c_fepsi',fepsi)
      call getrval('ithr',thr)
c
c  check slaves
c
      call para_check
c
c  pack and send data
c
      call para_initsend
cc      call para_pack_int(natom,1)
      call para_pack_int(nps,1)
      call para_pack_int(npspher,1)
      call para_pack_string(conel,256)
      call para_pack_real(fepsi,1)
      call para_pack_real(thr,1)
cc      call para_pack_real(xyz,3*natom)
      call para_bcast_pack(TxCosmoDat)
c need to be sent separately (memory management and MPI)      
      call para_bcast_real(cosurf,3*(nps+npspher),TxCosmoDat)
c
      return
c
      end
c======================================================================
      subroutine para_cosmo_surfrep(nps)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this subroutine drives the calculation of
c  the electron-surface repulsion matrix elements
c
c  this routine is executed by the master
c
c   nps           number of surface segments
c
      character*256 conel
c
      call getchval('c_onel',conel)
c
      if(conel.ne.'direct')then
c
c  get granularity
c
        call getival('c_gran',icgran)
c
        call para_check
c
c  assign work
c
        do iw=1,nps,icgran
          ie=min(iw+icgran-1,nps)
          call para_recv(imsg,islave,TxCosmoReq)
          call para_initsend
          call para_pack_int(iw,1)
          call para_pack_int(ie,1)
          call para_send_pack(islave,TxCosmoAss)
        enddo
c
c  all done. stop slaves
c
        do i=1,nslv
          call para_recv(imsg,islave,TxCosmoReq)
          call para_initsend
          call para_pack_int(0,1)
          call para_pack_int(0,1)
          call para_send_pack(islave,TxCosmoAss)
        enddo
      endif
c
      return
c
      end
c=======================================================================
      subroutine slave_cosmo_surfrep(vi,cosurf,inx,basdat,myseg,
     $                               nps,ncs,ncf,ntri)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this subroutine computes and stores the matrix elements
c  of electron-surface repulsion.
c
c  this is the slave code. The corresponding master code is
c  executed by para_cosmo_surfrep
c
c   vi            area for matrix elements for a surface segment(ntri)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   inx           contraction info
c   basdat        basis set data
c   myseg         array to keep track of which segment have been
c                 computed
c   nps           number of surface segments
c   ncs           number of contracted shells
c   ncf           number of basis function
c   ntri          ncf*(ncf+1)/2
c
      dimension vi(ntri),cosurf(3,nps)
      dimension inx(12,*),basdat(13,*),myseg(nps)
      character*256 conel
      logical busy
c
      call getchval('c_onel',conel)
c
c  if conel='direct' the integrals will be recomputed each time
c  they are needed, thus nothing has to be done here
c
      if(conel.eq.'direct')return
c
c
c   look for a free I/O channel
c
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(1,'slave_cosmo_surfrep',
     $              'cannot find a free channel',0,0)
 10   continue
c
c   open the direct access file
c
      open(unit=ich,file=conel(1:len_trim(conel)),status='unknown',
     $     form='unformatted',access='direct',recl=ntri*8)
c
c   ready to start working
c
      call izeroit(myseg,nps)
111   continue
c
c  request work from master
c
      call para_send(MY_GID,0,TxCosmoReq)
      call para_recv_pack(isource,TxCosmoAss)
      call para_unpack_int(iw,1)
      call para_unpack_int(ie,1)
      if(iw.gt.0)then
        do ip=iw,ie
c
c   compute matrix elements for surface segment ip
c
          call cosmo_inton(ip,vi,inx,basdat,cosurf,ncs)
c
c   save the current batch of matrix elements
c
          write(unit=ich,rec=ip)vi
          myseg(ip)=1
        enddo
c
c  next work
c
        goto 111
      endif
c
c  all done
c
      close(ich,status='keep')
      return
c
      end
c======================================================================
      subroutine para_cosmo_pot(den,phin,phi,nps,ntri)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this routine coordinates the work for the calculation of
c  the cosmo surface potential
c
c  den     density matrix
c  phi     nuclear contribution to surface potential (nps)
c  phi     surface potential (nps)
c  nps     number of surface segments
c  ntri    dimension of lower triangular matrix
c
      dimension den(ntri),phin(nps),phi(nps)
      character*256 conel
      logical conven
c
c  we need to send the density matrix to tha slaves,
c  as COSMO needs the full density, while slaves might
c  be using delta density
c
      call para_check
      call para_bcast_real(den,ntri,TxCosmoDen)
c
      call getchval('c_onel',conel)
      conven=conel.ne.'direct'
c
c  if the matrix elements are to be recomputed, we have to
c  coordinate the slaves
c
      if(.not.conven)then
c
c  get granularity
c
        call getival('c_gran',icgran)
c
        call para_check
c
c  assign work
c
        do iw=1,nps,icgran
          ie=min(iw+icgran-1,nps)
          call para_recv(imsg,islave,TxCosmoReq)
          call para_initsend
          call para_pack_int(iw,1)
          call para_pack_int(ie,1)
          call para_send_pack(islave,TxCosmoAss)
        enddo
c
c  all done. stop slaves
c
        do i=1,nslv
          call para_recv(imsg,islave,TxCosmoReq)
          call para_initsend
          call para_pack_int(0,1)
          call para_pack_int(0,1)
          call para_send_pack(islave,TxCosmoAss)
        enddo
      endif
c
c  now sum up the results with the nuclear contribution
c
      call zeroit(phi,nps)
      call para_reduce(phi,nps,TxCosmoPot)
      do ip=1,nps
        phi(ip)=phi(ip)+phin(ip)
      enddo
c
      return
c
      end
c======================================================================
      subroutine slave_cosmo_pot(den,inx,basdat,phi,cosurf,vi,myseg,
     $                           nps,ncf,ncs,ntri)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this subroutine computes the electronic part of the potential
c  on the cosmo surface
c
c  this is the slave code. The corresponding master code is
c  executed by para_cosmo_pot
c
c   den           area for density matrix (ntri)
c   inx           contraction info
c   basdat        basis set data
c   phi           on exit will contain potential (nps)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   vi            area for matrix elements for a surface segment(ntri)
c   myseg         array to keep track of which segment have been
c                 computed (nps)
c   nps           number of surface segments
c   ncf           number of basis function
c   ncs           number of contracted shells
c   ntri          ncf*(ncf+1)/2
c
      dimension den(ntri),phi(nps)
      dimension cosurf(3,nps),vi(ntri),myseg(nps)
      dimension inx(12,*),basdat(13,*)
      character*256 conel
      logical busy,conven
      parameter (zero=0.0d0)
c
c   receive density matrix
c
      call para_recv_bcastreal(den,ntri,TxCosmoDen)
c
      call getchval('c_onel',conel)
      conven=conel.ne.'direct'
      if(conven)then
c
c  coventional calculation. We process our own matrix elements file
c
c   look for a free I/O channel
c
        do i=20,99
          inquire(i,opened=busy)
          if(.not.busy)then
             ich=i
             goto 10
           endif
        enddo
        call nerror(1,'slave_cosmo_pot',
     $              'cannot find a free channel',0,0)
 10     continue
c
c   open the direct access file containing the surface matix elements
c
        open(unit=ich,file=conel(1:len_trim(conel)),status='old',
     $       form='unformatted',access='direct',recl=ntri*8)
c
c  loop over surface segments
c
        do ip=1,nps
          if(myseg(ip).ne.0)then
c
c   read matrix elements for surface segment ip
c
            read(unit=ich,rec=ip)vi
c
c   take trace with density
c
            call spur(den,vi,ncf,trace)
c
c   fill phi array
c
            phi(ip)=-trace
          else
            phi(ip)=zero
          endif
        enddo
c
c  close the cosmo matrix elements file
c
        close(ich,status='keep')
      else
c
c  direct calculation.
c
        call zeroit(phi,nps)
111     continue
c
c  request work from master
c
        call para_send(MY_GID,0,TxCosmoReq)
        call para_recv_pack(isource,TxCosmoAss)
        call para_unpack_int(iw,1)
        call para_unpack_int(ie,1)
        if(iw.gt.0)then
          do ip=iw,ie
c
c   compute matrix elements for surface segment ip
c
            call cosmo_inton(ip,vi,inx,basdat,cosurf,ncs)
c
c   take trace with density
c
            call spur(den,vi,ncf,trace)
c
c   fill phi array
c
            phi(ip)=-trace
          enddo
c
c  next work
c
          goto 111
        endif
      endif
c
c  all done. Send results to master
c
      call para_reduce(phi,nps,TxCosmoPot)
c
      return
      end
c======================================================================
      subroutine para_cosmo_h0(hmat,qcos,fepsi,nps,ntri)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this routine coordinates the work for the calculation of
c  the cosmo contribution to the one electron Hamiltonian
c
c   hmat            one-electron Hamiltonian (ntri)
c   qcos            unscaled screening charges (nps)
c   fepsi           scaling factor
c   nps             number of surface segments
c   ntri            dimension of triangular arrays
c
      dimension hmat(ntri),qcos(nps)
      character*256 conel
      logical conven
c
c  send surface charges to slaves
c
      call para_check
      call para_bcast_real(qcos,nps,TxCosmoCha)
c
      call getchval('c_onel',conel)
      conven=conel.ne.'direct'
c
c  if the matrix elements are to be recomputed, we have to
c  coordinate the slaves
c
      if(.not.conven)then
c
c  get granularity
c
        call getival('c_gran',icgran)
c
        call para_check
c
c  assign work
c
        do iw=1,nps,icgran
          ie=min(iw+icgran-1,nps)
          call para_recv(imsg,islave,TxCosmoReq)
          call para_initsend
          call para_pack_int(iw,1)
          call para_pack_int(ie,1)
          call para_send_pack(islave,TxCosmoAss)
        enddo
c
c  all done. stop slaves
c
        do i=1,nslv
          call para_recv(imsg,islave,TxCosmoReq)
          call para_initsend
          call para_pack_int(0,1)
          call para_pack_int(0,1)
          call para_send_pack(islave,TxCosmoAss)
        enddo
      endif
c
c  now sum up the results
c
      call zeroit(hmat,ntri)
      call para_reduce(hmat,ntri,TxCosmoH0)
c
      end
c======================================================================
      subroutine slave_cosmo_h0(hmat,vi,inx,basdat,cosurf,qcos,myseg,
     $                    fepsi,nps,ncf,ncs,ntri)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  this subroutine computes the cosmo contribution to the
c  one electron Hamiltonian
c
c  this is the slave code. The corresponding master code is
c  executed by para_cosmo_h0
c
c   hmat          one-electron Hamiltonian (ntri)
c   vi            area for matrix elements for a surface segment(ntri)
c   inx           contraction info
c   basdat        basis set data
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   qcos          unscaled screening charges (nps)
c   myseg         array to keep track of which segment have been
c                 computed (nps)
c   nps           number of surface segments
c   ncf           number of basis function
c   ncs           number of contracted shells
c   ntri          ncf*(ncf+1)/2
c
      dimension hmat(ntri),qcos(nps)
      dimension cosurf(3,nps),vi(ntri),myseg(nps)
      dimension inx(12,*),basdat(13,*)
      character*256 conel
      logical busy,conven
      parameter (zero=0.0d0)
c
c   receive surface charges
c
      call para_recv_bcastreal(qcos,nps,TxCosmoCha)
c
      call getchval('c_onel',conel)
      conven=conel.ne.'direct'
      if(conven)then
c
c  coventional calculation. We process our own matrix elements file
c
c   look for a free I/O channel
c
        do i=20,99
          inquire(i,opened=busy)
          if(.not.busy)then
             ich=i
             goto 10
           endif
        enddo
        call nerror(1,'slave_cosmo_h0',
     $              'cannot find a free channel',0,0)
 10     continue
c
c   open the direct access file containing the surface matix elements
c
        open(unit=ich,file=conel(1:len_trim(conel)),status='old',
     $       form='unformatted',access='direct',recl=ntri*8)
c
c  loop over surface segments
c
        call zeroit(hmat,ntri)
        do ip=1,nps
          if(myseg(ip).ne.0)then
c
c   read matrix elements for surface segment ip
c
            read(unit=ich,rec=ip)vi
c
c  add contribution of surface segment to one-electron Hamiltonian
c
            fq=-fepsi*qcos(ip)
            do i=1,ntri
              hmat(i)=hmat(i)+fq*vi(i)
            enddo
          endif
        enddo
c
c  close the cosmo matrix elements file
c
        close(ich,status='keep')
      else
c
c  direct calculation.
c
        call zeroit(hmat,ntri)
111     continue
c
c  request work from master
c
        call para_send(MY_GID,0,TxCosmoReq)
        call para_recv_pack(isource,TxCosmoAss)
        call para_unpack_int(iw,1)
        call para_unpack_int(ie,1)
        if(iw.gt.0)then
          do ip=iw,ie
c
c   compute matrix elements for surface segment ip
c
            call cosmo_inton(ip,vi,inx,basdat,cosurf,ncs)
c
c  add contribution of surface segment to one-electron Hamiltonian
c
            fq=-fepsi*qcos(ip)
            do i=1,ntri
              hmat(i)=hmat(i)+fq*vi(i)
            enddo
          enddo
c
c  next work
c
          goto 111
        endif
      endif
c
c  all done. Send results to master
c
      call para_reduce(hmat,ntri,TxCosmoH0)
c
      return
      end
c======================================================================
      subroutine para_cosmo_forc(den,cfor,qcos,cosurf,iatsp,
     $                           na,nps,ntri)

      use newpara

      implicit real*8 (a-h,o-z)
c
c   master routine for the calculation of the COSMO contribution to
c   the forces
c
c   den        density matrix (ntri)
c   cfor       cosmo forces (3,na)
c   qcos       surface charges (nps)
c   cosurf     surface coordinates (3,nps)
c   iatsp      mapping surface segment -> atoms (nps)
c               NOTE: COSMO does not use ghost atoms for surface
c                     construction, thus the mapping of iatsp
c                     works only if the centers of the calculation
c                     are ordered with the real atoms first,
c                     followed by the ghost atoms.
c   na         number of atoms
c   nps        number of surface segments
c   ntri       dimension of triangular array
c
      dimension den(ntri),cfor(3,na),qcos(nps),cosurf(3,nps)
      dimension iatsp(nps)
c
c  tell slaves to do cosmo forces
c
      call para_check
      call para_bcast(TxDoCosmoForc,TxJobType)
c
c  send  general data
c
      call para_send_info
c
c  pack  specific data for cosmo forces
c
      call para_initsend
      call para_pack_int(nps,1)
      call para_bcast_pack(TxCosmoDat)
c  need to send separately (b/ of mem. allocation in MPI)      
      call para_initsend
      call para_pack_real(den,ntri)
      call para_pack_real(qcos,nps)
      call para_pack_real(cosurf,3*nps)
      call para_pack_int(iatsp,nps)
c
c  send
c
      call para_bcast_pack(TxCosmoDat)
c
c  coordinate the slave work
c  get granularity
c
      call getival('c_gran',icgran)
c
      call para_check
c
c  assign work
c
      do iw=1,nps,icgran
        ie=min(iw+icgran-1,nps)
        call para_recv(imsg,islave,TxCosmoReq)
        call para_initsend
        call para_pack_int(iw,1)
        call para_pack_int(ie,1)
        call para_send_pack(islave,TxCosmoAss)
      enddo
c
c  all done. stop slaves
c
      do i=1,nslv
        call para_recv(imsg,islave,TxCosmoReq)
        call para_initsend
        call para_pack_int(0,1)
        call para_pack_int(0,1)
        call para_send_pack(islave,TxCosmoAss)
      enddo
c
c  now sum up the results
c
      call para_reduce(cfor,3*na,TxCosmoForc)
c
c  get cumulative cputime from slaves
c
      cosmotime=0.0d0
      call para_reduce(cosmotime,1,TxCosmoTime)
      call setrval('cpucosmo',cosmotime)
c
      return
      end
c======================================================================
      subroutine do_cosmoforc

      use newpara
      use memory

      implicit real*8 (a-h,o-z)
c
c   slave routine for the calculation of the COSMO contribution to
c   the forces
c
c get the corresponding host name
c
      cosmotime=0.0d0
      call secund(tstart)
c
c -- receive the geometry,cpu,symmetry data
c
      call para_get_info
c
c  get values from depository
c
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      ntri=ncf*(ncf+1)/2
c
c  get specific cosmo data
c
      call para_recv_bcastpack(TxCosmoDat)
      call para_unpack_int(nps,1)
c
      call getmem(ntri,iden)
      call getmem(nps,iqcos)
      call getmem(3*nps,icosurf)
      call getmem(3*na,icfor)
      call getint(nps,iiatsp)
c two separate packages (to make memory allocation safe)      
      call para_recv_bcastpack(TxCosmoDat)
      call para_unpack_real(bl(iden),ntri)
      call para_unpack_real(bl(iqcos),nps)
      call para_unpack_real(bl(icosurf),3*nps)
      call para_unpack_int(bl(iiatsp),nps)
c
c  ready to go
c
      call zeroit(bl(icfor),3*na)
111   continue
c
c  request work from master
c
      call para_send(MY_GID,0,TxCosmoReq)
      call para_recv_pack(isource,TxCosmoAss)
      call para_unpack_int(iw,1)
      call para_unpack_int(ie,1)
      if(iw.gt.0)then
        do ip=iw,ie
        itipm1=3*(ip-1)
        call cosmo_getatsp(bl(iiatsp),ip,nps,natsp)
c
c   compute matrix elements for surface segment ip
c
        call cosmo_intof7(ip,bl(ictr),bl(ibas),bl(icosurf+itipm1),
     $                    natsp,ncs,ncf,ntri,bl(iden),
     $                    bl(icfor),bl(iqcos+ip-1))
        enddo
c
c  next work
c
        goto 111
      endif
c
c  all done. Send results to master
c
      call para_reduce(bl(icfor),3*na,TxCosmoForc)
      call secund(tend)
      cosmotime=cosmotime+tend-tstart
      call para_reduce(cosmotime,1,TxCosmoTime)
c
c  clear memory
c
cc      call retmark
cc      call retimark
c
      return
      end
c======================================================================
      subroutine cosmo_getatsp(iatsp,ip,nps,natsp)
      implicit real*8 (a-h,o-z)
c
c   this subroutine returns in natsp the atom to which surface segment
c   ip belongs to. The mapping surface segment --> atoms is stored
c   in the array iatsp (nps).
c
c   I use a subroutine for doing this because the memory for the
c   array iatsp is allocated in the real memory segment, thus for
c   the calling routine it is not easy to acces all the elements
c   of the array
c
      dimension iatsp(nps)
c
      natsp=iatsp(ip)
      end
