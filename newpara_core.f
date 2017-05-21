c This file contains the common routines for the paralel
c
c para_next:      tell slaves when iterative process over and get timings
c para_send_info: pack in the basic data for the slaves
c para_done:      get the slaves move on in the Hessian CPHF      
c para_distr_blocks:   handle integral block requests/assignments on master    
c para_get_block:       receive integral blocks on the slaves      
c para_bcast_real:      wrapper for big data chunks
c para_reduce:          ditto for reduce
c sortat:               sort atoms before distribution
c para_distr_atoms:     handle atom assignments on master      
c para_get_atom:        receive atom assignement on slaves      
c para_startmsg         compose parallel startup message
c para_affinity         compute CPU affinity mask
c
c====================================================================
c
      subroutine para_next(mgo)
c
c====================================================================
      
      use newpara

      character*60 hostname
      real*8 ttm,totcpu,zero,tccm,tcosmo
      parameter(zero=0.0d0)
      integer my_mgo
c
c  send next step mgo=0 - finish
c  1,2,3 start over scf with different thresholds
c  -1 - finish, gather timing data, but do NOT broadcast <TxNext>
c
      my_mgo = mgo  ! copy the argument into local storage,
                    ! because it needs to be writable in MPI calls
      If(my_mgo.NE.-1) call para_bcast(my_mgo,TxNext)
c
c  get timing messages if finished
c
      If(my_mgo.LE.0) Then
        totcpu=zero
        write(91,100)
        do iproc=1,nslv
        call para_recv_pack(ifrom,TxTime)
          call para_unpack_real(ttm,1)
          call para_unpack_string(hostname,60)
          ttmm=ttm/60.0d0
          call rmblan(hostname,60,lhostn)
          write(91,110) hostname(1:lhostn),ttm,ttmm
          totcpu=totcpu+ttm
        end do
c
c  there might be a cosmo step whose cpu time has not yet
c  been accounted for
c
        call tstrval('cpucosmo',icosmo)
        if(icosmo.ne.0)then
          call getrval('cpucosmo',tcosmo)
          totcpu=totcpu+tcosmo
          call setrval('cpucosmo',0.0d0)
        endif
c
c  cpus holds the total cpu time used by all slaves in the last step
c  cpuscum holds the cumulative cpu time used by the slaves
c
        call setrval('cpus',totcpu/60.0d0)
        call getrval('cpuscum',tccm)
        call setrval('cpuscum',tccm+totcpu/60.0d0)
      EndIf
 100  format('CPU time for daughter processes')
 110  format('process on ',a20,f10.2,' s',' (',f9.3,' m)')
      end
c
c==========================================================================
c
      subroutine para_send_info
c
c==========================================================================
c
c  send initialization data to the slaves. Updated to separate
c  the sending of single variables and of arrays, because in the MPI
c  version it is better not to allocate memory while packing and
c  unpacking data, as the memory array is used also as buffer for
c  send/receive
c
      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      common /intlim/ limxmem,limblks,limpair
c
c  get the data from the depository and pack it
c
c
c  initialization variables
c
      call para_initsend
c
c  output channel
c
      call tstival('iout',iout)
      if(iout.eq.0) then
        iout=6
      else
        call getival('iout',iout)
      endif
      call para_pack_int(iout,1)
c
c  geometry and basis set info
c
      call getival('na  ',na)
      call getival('ndum',ndum)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('ictr',ictr)
      call tstival('nonredun',iii)
      if(iii.eq.0) then
        nonredun=ncf
        call setival('nonredun',ncf)
      else
        call getival('nonredun',nonredun)
      endif
c
      call para_pack_int(na,1)
      call para_pack_int(ndum,1)
      call para_pack_int(ncf,1)
      call para_pack_int(ncs,1)
      call para_pack_int(nbf,1)
      call para_pack_int(nsh,1)
      call para_pack_int(nonredun,1)
c
c  COSMO flag
c
      call tstival('cosmo',icosmo)
      if(icosmo.eq.0) then
        call setival('cosmo',0)
      else
        call getival('cosmo',icosmo)
      endif
      call para_pack_int(icosmo,1)
c
c  CPU related information
c
      call getival('ints',intsize)
      call getival('accu',iacc)
      call getival('cach',icache)
      call getival('memr',memreal)
c
      call para_pack_int(intsize,1)
      call para_pack_int(iacc,1)
      call para_pack_int(icache,1)
      call para_pack_int(memreal,1)
c
c  these are not taken from the depository, it's not updated
c
      call para_pack_int(limxmem,1)
      call para_pack_int(limblks,1)
      call para_pack_int(limpair,1)
c
c  symmetry flag
c
      call getival('nsym',nsym)
      call para_pack_int(nsym,1)
c
      call para_bcast_pack(TxInitVar)
c
c  initialization arrays
c
      call para_initsend
      call para_pack_real(bl(ibas),13*nsh)
      call para_pack_real(bl(inuc),5*na)
      call para_pack_int(bl(ictr),12*ncs)
c
      call para_bcast_pack(TxInitArr)
c
c  Symmetry data
c
      if (nsym.gt.0) then
         call getival('ngener',ngener)
         call getival('nsyo',iadr)
         call getival('SymNuPr',nupair)
         call getival('SymNuPr1',nupair1)
         call getival('SymFunPr',ifp)
         call getival('SymShPr',ifp1)
         call para_initsend
         call para_pack_int(ngener,1)
         call para_pack_int(bl(iadr),7)
         call para_pack_int(bl(nupair),na*nsym)
         call para_pack_int(bl(ifp),7*ncf)
         call para_pack_int(bl(ifp1),7*ncs)
         if (ndum.gt.0) then
           natom = na-ndum
           call para_pack_int(bl(nupair1),natom*nsym)
         endif
c
        call para_bcast_pack(TxSymData)
      endif
      end
c====================================================================
c
      subroutine para_done(iarg)
c
c====================================================================

      use newpara
      implicit none
      integer iarg,imsg,iproc,islave

c send cphf termination to slaves
c
      call para_bcast(iarg,TxHessDat)
c
      If(iarg.EQ.-1) return
c
      do iproc=1,nslv
      call para_recv(imsg,islave,TxContinue)
      enddo
c
      call para_next(0)
c
      end
c====================================================================
c
      subroutine para_distr_blocks_o(nblocks,igran)
c
c====================================================================
c Give assignments to slaves requesting work
c Start with size igran blocks, then switch to single blocks
c to level out load imbalance.
c more slaves and rougher grain cause earlier switching
c slaves do the first 3*igran*nslv blocks w/o asking
c
      use newpara

      implicit real*8 (a-h,o-z)
c
      n85percent=nblocks-2*igran*nslv
      if(n85percent.le.3*igran*nslv) n85percent=3*igran*nslv+1
c
c make it an even multiple of igran
c
      n85percent=(n85percent/igran)*igran
      n85percent=nblocks-n85percent
c
      call para_check
c
      do iwork=nblocks-3*igran*nslv-igran+1,n85percent+1,-igran
         call para_recv(imsg,islave,TxBlockReq)
         call para_send(iwork,islave,TxBlockAssign)
      end do
c
c switch to single blocks
c the negative value of ijob signals the switch to the slaves
c send 0 if ready
c first check if the slaves are OK
c
      call para_check
      call elapsec(timo)
      call getrval('imbal',tlost)
      do iwork=n85percent,1-nslv,-1
         call para_recv(imsg,islave,TxBlockReq)
         if(iwork.lt.1) then
            call elapsec(timb)
            tlost=tlost-iwork*(timb-timo)
            timo=timb
            ijob=0
         else
            ijob=-iwork
         end if
         call para_send(ijob,islave,TxBlockAssign)
      end do
      call setrval('imbal',tlost)
      write (91,*) 'block imbal', tlost/60.0d0
      end
c====================================================================
c
      subroutine para_distr_blocks(nblocks,igran)
c
c====================================================================
c Give assignments to slaves requesting work
c Use the guided self-scheduling algorithm
c to level out load imbalance.
c Start from the last block and always send an interval
c depending on the no of remaining blocks and active slaves
c interval=work left/(nslv*variance)
c variance is the estimate for the ratio between the individual times
c one block might take to finish on a machine
c with more slave processors it was quite useful
c probably most cosmo and FTC calls can be converted to use this
c
c  MM 02/21/08: modified the while loop to stop when there are
c               no more blocks. Before it was possible to send
c               to a slave iwork=0, which is quite bad as this
c               might trigger the slave into computing the whole
c               set of integrals. Also, the way the imbalance time
c               was computed did not make any sense to me. Now,
c               after the loop is completed, the master send the
c               termination message to each slave, and at that point
c               it computes the time lost due to imbalamce, as now it
c               has to wait for each slave to finish.
c
      use newpara

      implicit none
      integer nblocks,igran
      integer ilast,ivari,imsg,islave,ilen,iwork,ic
      real*8 timo,tlost,timb
c
      ilast=nblocks
      ivari=igran      
c the default igran value of 20 was quite ok
c and eliminated imbalance in the integral jobs
c however the meaning of igran has changed
c the bigger value decreases the size of jobs sent out
c the documentation has to be updated!!!!      
c
      do while (ilast.gt.0)
        call para_recv(imsg,islave,TxBlockReq)
        ilen=max(1,ilast/(NSLV*ivari))  ! chunk length for this job
        iwork=ilast-ilen+1              ! start for this job
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(ilen,1)
        call para_send_pack(islave,TxBlockAssign) !sent
        ilast=ilast-ilen                ! last piece not done
      end do
c
c  all done, stop the slaves. After a first slave is done,
c  any time spent waiting for the other slaves to finish
c  is added to the imbalance time
c
      call para_recv(imsg,islave,TxBlockReq)
      call para_initsend
      call para_pack_int(-1,1)
      call para_pack_int(-1,1)
      call para_send_pack(islave,TxBlockAssign) ! stop first slave
c
      call elapsec(timo)
      call getrval('imbal',tlost)
      do ic=1,NSLV-1
        call para_recv(imsg,islave,TxBlockReq)
        call para_initsend
        call para_pack_int(-1,1)
        call para_pack_int(-1,1)
        call para_send_pack(islave,TxBlockAssign)
      enddo
      call elapsec(timb)
      tlost=tlost+(timb-timo)
      call setrval('imbal',tlost)
      end
c====================================================================
c
      subroutine para_get_block(iwork,igran)
c
c====================================================================
c Get block assignment from master
c iwork is the starting point of the interval, igran is the length
c
      use newpara

      implicit none
      integer iwork,igran,ifrom
c
c request more work from the master
c
      call para_send(MY_GID,0,TxBlockReq)
      call para_recv_pack(ifrom,TxBlockAssign)
      call para_unpack_int(iwork,1)
      call para_unpack_int(igran,1)
      end
c======================================================================      
c
      subroutine sortat(iunq,isortunq,nq)
      use memory
      implicit real*8(a-h,o-z)
c
c  Sorts pointer array to atoms 
c  so that the atoms pointed to iunq(x) follow with decreasing atomic no
c  useful to decrease parallel imbalance for stuff parallel over atoms
c  hopefully dummy atoms will not break it       
c  uses the shell sort algorithm      
c     
c  iunq:        pointer to equivalent atoms array
c  isortunq:    new pointer sorted by atomic charge      
c  nq:          number of equivalent atoms
      dimension iunq(nq),isortunq(nq),chrg(nq)
c
      parameter (aln2i=1.0d0/0.69314718d0,small=1.0d-5)
c obtain atomic charges and fill new pointer array
      call getival('inuc',inuc)
      do 10 i=1,nq
         isortunq(i)=iunq(i)
c charge is the first value in each row of nucdat(5,iat)         
         chrg(i)=bl(inuc+5*(iunq(i)-1))
 10   enddo
c
c  shell sort
c
      lognb2 = int( log(float(nq))*aln2i + small )
      m = nq
c
      do 40 nn=1,lognb2
         m = m/2
         k = nq-m
         do 30 j=1,k
            i = j
 20         continue
            l = i+m
            if(chrg(l).gt.chrg(i)) then
               charge=chrg(i)
               ind = isortunq(i)
               chrg(i) = chrg(l)
               isortunq(i) = isortunq(l)
               chrg(l) = charge
               isortunq(l) = ind
               i = i-m
               if(i.ge.1) go to 20
            endif
 30      continue
 40   continue
cc
c
      end subroutine sortat
c      
c====================================================================
c
      subroutine para_distr_atoms(iunq,nq)
c
c====================================================================
c Distribute atoms between slaves in a round-robin manner
c Only unique atoms are sent
c The atom list is sorted with decreasing nuclear charge       
c to level out load imbalance.
c
c iunq: pointer to symmetry unique atoms
c nq:   number of symmetry unique atoms      
c      
      use newpara

      implicit real*8 (a-h,o-z)
      dimension iunq(nq),isortunq(nq)
C
      If(nslv.GT.nq) Then
c -- inefficiency warning
        call message('**WARNING** in DFT',
     $  'Not possible to use slaves efficiently: Nslaves > Natoms ',
     *  nslv,NQ)
      EndIf
C sort the atoms
      call sortat(iunq,isortunq,nq)
c check slaves      
      call para_check
c
      DO 200 IAtom=1,NQ
      ICntr = IsortUNQ(IAtom)
      call para_recv(imsg,islave,TxDftReq)
      call para_send(ICntr,islave,TxDftAssign)
 200  CONTINUE
c
c -- all done  stop all slaves
      call elapsec(timo)
      call getrval('imbal',tlost)
      Do IAtom=1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_send(0,islave,TxDftAssign)
      call elapsec(timb)
      tlost=tlost+(IAtom-1)*(timb-timo)
      timo=timb
      EndDo
      call setrval('imbal',tlost)
c
      end subroutine para_distr_atoms
c      
c====================================================================
c
      subroutine para_get_atom(icntr)
c
c====================================================================
c Get atom assignment from master
c really simple
c
      use newpara

      implicit real*8 (a-h,o-z)
c
c -- request atom we are currently doing from master
      call para_send(MY_GID,0,TxDftReq)
      call para_recv(ICntr,isource,TxDftAssign)
c      
      end
c      
c====================================================================
c
      subroutine para_startmsg
c
c====================================================================
c  compose the parallel startup message

      use newpara
      use messages

      implicit none

      integer i, ii, j, nslvi, ncpui, jstrt, jend, k, nloop
      character*256 wline

      call addbuf('')
      wline=''
      write(wline,'(a,i0)')' Master process on '//trim(MY_HOSTNM)
     $               //'. Process ID: ',MY_PID
      call addbuf(trim(wline))
      call addbuf('')
      wline=''
      if( NHOST .gt. 1 )then
        write (wline,'(1x,2(i0,a))') NSLV,' slaves distributed over ',
     $   NHOST,' hosts are working on your job.'
      else
        write (wline,'(1x,i0,a)') NSLV,' slaves grouped on '//
     $   '1 host are working on your job.'
      endif
      call addbuf(trim(wline))

c     slaves on master host

      nslvi = NPROCH(1)-1
      ncpui=NCPUH(1)
      if( nslvi .gt. 0 ) then
        call addbuf('')
        wline=''
        write(wline,'(a,i0,a)')
     $    ' Slaves running on '//trim(HOSTNAMES(1))//' (',
     $     ncpui,' CPU):'
        call addbuf(trim(wline))
c
        nloop = nslvi/10 + 1
        jstrt = 1
        DO k=1,nloop
        jend = MIN(jstrt+9,nslvi)
        wline=''
        write(wline,'(1x,a,10(i8,1x))') 'Group ID     : ',
     $                        (IDH(1,j+1),j=jstrt,jend)
        call addbuf(trim(wline))
        wline=''
        write(wline,'(1x,a,10(i8,1x))') 'Process ID   : ',
     $                          (PIDSL(IDH(1,j+1)),j=jstrt,jend)
        call addbuf(trim(wline))
        wline='' 
        write(wline,'(a)') ' '
        call addbuf(trim(wline))
        jstrt = jstrt+10
        EndDO
c
        if( nslvi .gt. ncpui ) then
          call addbuf(' *** Warning: host '//trim(HOSTNAMES(1))
     $                //' is oversubscribed')
        endif
      endif

c     slaves on remote hosts

      if( NHOST .gt. 1 )then
        do i=2,NHOST
          call addbuf('')
          nslvi = NPROCH(i)
          ncpui=NCPUH(i)
          if( nslvi .gt. 0 ) then
            wline=''
            write(wline,'(a,i0,a)')
     $        ' Slaves running on '//trim(HOSTNAMES(i))//' (',
     $         ncpui,' CPU):'
            call addbuf(trim(wline))
c
            nloop = nslvi/10 + 1
            jstrt = 1
            DO k=1,nloop
            jend = MIN(jstrt+9,nslvi)
            wline=''
            write(wline,'(1x,a,10(i8,1x))') 'Group ID     : ',
     $                           (IDH(i,j),j=jstrt,jend)
            call addbuf(trim(wline))
            wline=''
            write(wline,'(1x,a,10(i8,1x))') 'Process ID   : ',
     $                         (PIDSL(IDH(i,j)),j=jstrt,jend)
            call addbuf(trim(wline))
            jstrt = jstrt+10
            wline='' 
            write(wline,'(a)') ' '
            call addbuf(trim(wline))
            EndDO
c
            if( nslvi .gt. ncpui ) then
              call addbuf(' *** Warning: host '//trim(HOSTNAMES(i))
     $                    //' is oversubscribed')
            endif
          endif
        enddo
      endif
      call addbuf('')

      contains

      character*8 function mask8(j,ncpu,imask)

        implicit none
        integer j, ncpu
        integer*4 imask(ncpu)
        integer i, ii

        mask8=''
        ii=(j-1)*8
        do i=1,8
          ii=ii+1
          if(ii.gt.ncpu) exit
          write(mask8(8-i+1:8-i+1),'(i1)')imask(ii)
        enddo

      end function mask8

      end

c      
c====================================================================
c
      subroutine para_affinity(gid,maxp,ncpus,idhost,imask)
c
c====================================================================
c compute the cpu affinity mask
c
c  gid    (input)  group id of the process for which the cpu mask
c                  has to be computed
c  maxp   (input)  maximum number of processes running on a host
c  ncpus  (input)  number of CPUs of the host
c  idhost (input)  array containing the group ids of the processes
c                  running on the host
c
c  imask  (input)  initial value of the CPU affinity mask (should
c                  show which CPU are available to the processes
c                  running on the chosen host)
c  imask  (output) will contain the mask for the chosen process
c
      use newpara

      implicit none

      integer gid, maxp, ncpus
      integer*4 imask(ncpus)
      integer idhost(maxp)

      integer*4 iret, getaffinity, setaffinity
      integer  ncpua, nproc, ncpup, ncpur, i, is, ig, il, igrp

               ! number of processes to accomodate
               ! (without counting the master)

      nproc=0
      do i=1,maxp
        if( idhost(i) .gt. 0 ) nproc=nproc+1
      enddo
      if( nproc .eq. 0) return ! nothing to be done

               ! number of available cpus

      ncpua=0
      do i=1,ncpus
        if( imask(i) .ne. 0 ) ncpua=ncpua+1
      enddo

      if( nproc .gt. ncpua ) return ! the cpus are oversubscribed

              ! number of cpus per process

      ncpup = ncpua / nproc
      ncpur = mod( ncpua, nproc )

              ! identify which group of cpus to allocate

      igrp=0
      if( gid .ne. 0 ) then
        if( idhost(1) .eq. 0 ) then
          is = 1
        else
          is = 0
        endif
        do i=1,nproc
          if( idhost(is+i) .eq. gid ) then
            igrp = i
            exit
          endif 
        enddo
      else
        igrp = 1               !  master uses first group
      endif
      if( igrp .eq. 0 ) return ! cannot identify group

      i = 0    ! cpu counter
      ig = 1   ! group counter

      if( ncpur .gt. 0 ) then  ! group length
        il = ncpup + 1
      else
        il = ncpup 
      endif

      is = 0   ! slot counter
      do
        i = i + 1
        if ( i .gt. ncpus )exit
        if( imask(i) .eq. 1 ) then
          is = is + 1
          if( is .gt. il ) then
            ig = ig + 1
            if( ig .le. ncpur ) then
              il = ncpup + 1
            else
              il = ncpup 
            endif
            is = 1
          endif
          if( ig .ne. igrp ) imask(i) = 0
        endif
      enddo

      end
