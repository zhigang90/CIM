
subroutine para_driver

! This file contains the PVM specific paralel routines:
!
! para_driver:    startup and forking between master and slave
! start_slaves    (function obvious)
! para_stop:      cleanup slaves in case of nonregular stop
! para_check:     check the slaves - and stop if any of them died

! set_pvm_def_int set the default integer type for PVM

! para_initsend:  initialize buffers
! para_pack_....  pack data(int,real,string,byte)      
! para_unpack_... unpack data(int,real,string,byte)      

! para_send_....  send msg w/ tag to host specified
! para_recv_....  recv msg w/ tag from any host, 
!                          source host is returned
!!!! para_recv_real receives only from host specified 

! para_bcast_.... msg w/ tag from master to all slaves
! para_recv_bcast.... recv msg w/ tag from master

! para_....()     msg is a single integer
! para_...._pack  msg is packed
! para_...._real  msg is a real array

! para_reduce     sum up array on master
! para_probe      check if message has arrived
! para_mem, para_slave_mem         send/receive memory limits

! para_..._..._flat is a version of the routine not using a 
!                   binary tree, but the slower native pvm calls
! para_printhelp  prints out the usage message
!
! para_sendbin    send bin to another slave during MP2 binsort


! here follows the code for para_driver:

! initialization, forking between slave and master      

  use pvm_defs
  use newpara
  use messages
  implicit none

  integer*4 info,mytid,itid,mygid,bufid,pidslv
  integer*4 get_pid,getaffinity,setaffinity,ncpu4
  integer iarg,pvmd,ncpu,i,grpsize,iret
  integer nhst,idslv,ncpus,ncpuslv,ii,j,ifrom
  character*256 hostfile,file,argument,prodstr,slname
  character*256 progname,inpname,wline
  external master_die, slave_die

  call set_pvm_def_int   ! set pvm type of default integer
!
!  set product string for license checking (override prodstr with parastr)
!
  call getchval('parastr',prodstr)
  call setchval('prodstr',prodstr)

  call parsearguments(progname,inpname) ! validate the program arguments


  call setchval('msgpass','PVM') ! set the message passing type

  call blankit(hostfile,256)
  call blankit(file,256)

          ! Temporarily supress the PVM errors

  call pvmfsetopt(PvmAutoErr,0,info)

          ! start the pvm demon
          ! the hostfile is the first argument
          ! if not found use the standard $HOME/.txhosts
          ! if the demon is already running, use the present configuration

  file='/.txhosts'

  call get_environment_variable('HOME',hostfile)
  call concat(hostfile,file)

  wbuf=''          ! zero out message buffer

  iarg=command_argument_count()
  if (iarg.gt.1) then
     i=1
     do
      call blankit(argument,256)
        call get_command_argument(i,argument)
        
        if (argument.eq.'-np'.and.i.lt.iarg) then
           call get_command_argument(i+1,argument)
           read (argument,*) ncpu
           NSLV = ncpu - 1 ! NSLV is part of newpara module data
           call setival('nslv',NSLV)
           i=i+2

        else if (argument.eq.'-f'.and.i.lt.iarg) then

             !   -f hostfile (machines participating in the 
             !                virt. machine listed there)

           call blankit(hostfile,256)
           call getarg(i+1,hostfile)
           i=i+2
        else
           i=i+1
        end if
        if (i.ge.iarg) exit
     enddo
  end if

  call pvmfstartpvmd(hostfile,1,info)
  if (info.eq.0) then

             ! This will tell if we need to shutdown the demon

     call setival('pvmd',1)
     call addbuf(' The virtual machine setup is based on: '//trim(hostfile))
  else if (info.ne.PvmDupHost) then

             ! Halt if the virtual machine did not start up

     call pvmfperror('Virtual machine did not start',info)
     call pvmfhalt(info)
  else
      call addbuf(' The existing virtual machine is used.')
  end if

             ! PVM errors enabled

  call pvmfsetopt(PvmAutoErr,1,info)

             ! setup direct routing

  call para_directroute
! call pvmfsetopt(PvmRoute,PvmRouteDirect,info)

             ! setup fragment size for direct routing (max. poss.)

  call pvmfsetopt(PvmFragSize,1048576,info)

             ! what is my tid? 

  call pvmfmytid(mytid)      

             ! get tid of parent process

  call pvmfparent(itid)

             ! set the data of the pvm_def module data block

  MY_ID = mytid
  MASTER_ID = itid
  MY_PID = get_pid()
  call gethstnm   ! this call will set MY_HOSTNM and PVM_NH

  if (itid.eq.PvmNoParent) then

             ! I am the master
             !
             ! supress the 'BEGIN' and 'END' strings
             ! but catch the output of the slaves (that should be only errors)

     call pvmfsetopt(PvmShowTids,0,info)
     call pvmfcatchout(1,info)

     call blankit(GROUP_NAME,79)
     write(GROUP_NAME,*) mytid

             ! join our calculational group, get group id

     call pvmfjoingroup(GROUP_NAME,mygid)
     MY_GID = mygid

             ! zero out the counter for slave's cpu time

     call setrval('cpuscum',0.0d0)
     call setrval('cpus',0.0d0)
     call setrval('reduce',0.0d0)
     call setrval('bcast',0.0d0)
     call setrval('imbal',0.0d0)

             ! let's start the slaves

     call start_slaves
 
     call para_gettids ! make sure slave ids are consistent

             ! gather id, names of slaves and number of CPU per host

     allocate( HOSTNAMES(NSLV+1), PIDSL(NSLV) )
     allocate( NPROCH(NSLV+1), NCPUH(NSLV+1) )
     NPROCH=0
     nhst=1
     HOSTNAMES(nhst)=MY_HOSTNM    ! the first host is the master's
     call getival('ncpus',ncpus)
     ncpu4=ncpus
     NCPUH(nhst)=ncpus
     NPROCH(nhst)=NPROCH(nhst)+1
     do i=1,NSLV
        call para_recv_pack(ifrom,TxId)
        call para_unpack_int(idslv,1)
        call para_unpack_int(ncpuslv,1)
        call para_unpack_int4(pidslv,1)
        call para_unpack_string(slname,256)
        HOSTSL(idslv)=slname
        PIDSL(idslv)=pidslv
        do ii=1,nhst
           if (slname.eq.HOSTNAMES(ii))then
              NPROCH(ii)=NPROCH(ii)+1
              exit
           endif
        enddo
        if (ii.gt.nhst)then
           nhst=nhst+1
           HOSTNAMES(nhst)=slname
           NCPUH(nhst)=ncpuslv
           NPROCH(nhst)=NPROCH(nhst)+1
        endif
     enddo

           !  number of hosts, maximum number of processes and CPU per host

     NHOST=nhst
     MAXPROCH=NPROCH(1)
     MAXCPUH=NCPUH(1)
     do i=1,NHOST
        if(NPROCH(i).gt.MAXPROCH)MAXPROCH=NPROCH(i)
        if(NCPUH(i).gt.MAXCPUH)MAXCPUH=NCPUH(i)
     enddo

       !  fill in the database of process id per host

     allocate( IDH(NHOST,MAXPROCH) )
     IDH=-1
     MY_HOSTID=1
     IDH(1,1)=0   ! the master is the first process of the first host
     do i=1,NSLV
        slname=HOSTSL(i)
        do ii =1,NHOST
           if (slname.eq.HOSTNAMES(ii))then
              do j=1,MAXPROCH
                 if (IDH(ii,j).eq.-1)then
                    IDH(ii,j)=i
                    exit
                 endif
              enddo
              exit
           endif
        enddo
     enddo
             ! broadcast the host database

     call para_bcast(NHOST,TxId)
     call para_bcast(MAXPROCH,TxId)
     call para_bcast(MAXCPUH,TxId)

     call para_initsend
     do i=1,NHOST
        call para_pack_string(HOSTNAMES(i),256)
     enddo
     call para_pack_int(NPROCH,NHOST)
     call para_pack_int(NCPUH,NHOST)
     call para_pack_int(IDH,NHOST*MAXPROCH)
     call para_pack_int4(PIDSL,NSLV)
     call para_bcast_pack(TxId)

             !  initialize cpu affinity mask

     allocate( ICPUMASK(ncpus), ICPU(ncpus) )
     allocate( MASKSL(MAXCPUH,NSLV) )
     MASKSL=0
     ICPUMASK=0
     info = getaffinity( ncpu4, ICPUMASK )
     ICPU=ICPUMASK

             !  set cpu affinity

     if (info.eq.0)then
        call para_affinity(MY_GID,MAXPROCH,ncpus,IDH(MY_HOSTID,1),ICPU)
        info = setaffinity( ncpu4, ICPU )
        if (info.ne.0)then
           call addbuf('')
           call addbuf(' *** Warning: Master cannot set cpu affinity')
        endif
     else
        call addbuf('')
        call addbuf(' *** Warning: Master cannot set cpu affinity')
     endif

             ! gather affinity masks from slaves

     do i=1,NSLV
        call para_recv_pack(ifrom,TxMask)
        call para_unpack_int(idslv,1)
        call para_unpack_int(ncpuslv,1)
        call para_unpack_int4(MASKSL(1,idslv),ncpuslv)
     enddo
             ! broadcast the mask database

     call para_initsend
     call para_pack_int4(MASKSL,NSLV*MAXCPUH)
     call para_bcast_pack(TxMask)

             ! compose the startup message

     call para_startmsg

             ! install message handler for master

     call pvmfaddmhf( -1, TxKissofDeath, -1, bufid, master_die, info )

             ! now do the work

     call driver

     call para_bcast(TxFinish,TxJobType) !stop slaves
  else

            ! I am a slave

     call para_recv_bcast(nhst,TxNslave)
     if (nhst.ge.0) then !normal slaves
        NSLV=nhst
        call setival('nslv',NSLV)
 
        allocate( SLVID(NSLV+1), HOSTSL(NSLV) )

        call blankit(GROUP_NAME,79)
        write(GROUP_NAME,*) itid

               ! join our calculational group, get group id

        call pvmfjoingroup(GROUP_NAME,mygid)
        MY_GID = mygid

        grpsize=NSLV+1
        call pvmffreezegroup(GROUP_NAME,grpsize,info)
        call para_gettids ! set up the slave ids

              ! send slave id, number of cpu and hostname

        call para_initsend
        call para_pack_int(MY_GID,1)
        call getival('ncpus',ncpus)
        ncpu4=ncpus
        call para_pack_int(ncpus,1)
        call para_pack_int4(MY_PID,1)
        call para_pack_string(MY_HOSTNM,256)
        call para_send_pack(0,TxId)

              ! receive the host database

        call para_recv_bcast(NHOST,TxId)
        call para_recv_bcast(MAXPROCH,TxId)
        call para_recv_bcast(MAXCPUH,TxId)

        allocate( PIDSL(NSLV) )
        allocate( HOSTNAMES(NHOST) )
        allocate( NPROCH(NHOST), NCPUH(NHOST) )
        allocate( IDH(NHOST,MAXPROCH) )

        call para_recv_bcastpack(TxId)
        do i=1,NHOST
           call para_unpack_string(HOSTNAMES(i),256)
        enddo
        call para_unpack_int(NPROCH,NHOST)
        call para_unpack_int(NCPUH,NHOST)
        call para_unpack_int(IDH,NHOST*MAXPROCH)
        call para_unpack_int4(PIDSL,NSLV)

                ! compute MY_HOSTID

        do i=1,NHOST
           if (HOSTNAMES(i).eq.MY_HOSTNM)then
              MY_HOSTID=i
              exit
           endif
        enddo

  !  initialize cpu affinity mask

        allocate( ICPUMASK(ncpus), ICPU(ncpus) )
        allocate( MASKSL(MAXCPUH,NSLV) )
        MASKSL=0
        ICPUMASK=0
        info = getaffinity( ncpu4, ICPUMASK )
        ICPU=ICPUMASK

  !  set cpu affinity

        if (info.eq.0)then
           call para_affinity( MY_GID, MAXPROCH, ncpus, IDH(MY_HOSTID,1), ICPU)
           info = setaffinity( ncpu4, ICPU )
           if (info.ne.0)then
              write(*,'(a,i0,a)')'*** Warning: Slave ',MY_GID,' cannot set cpu affinity'
           endif
        else
           write(*,'(a,i0,a)')'*** Warning: Slave ',MY_GID,' cannot set cpu affinity'
        endif

              ! send affinity mask

        call para_initsend
        call para_pack_int(MY_GID,1)
        call para_pack_int(ncpus,1)
        call para_pack_int4(ICPU,ncpus)
        call para_send_pack(0,TxMask)

              ! receive back the mask database

        call para_recv_bcastpack(TxMask)
        call para_unpack_int4(MASKSL,NSLV*MAXCPUH)

              ! install message handler for slave

        call pvmfaddmhf( -1, TxKissofDeath, -1, bufid, slave_die, info )

               ! now do the work

        call slave
     else       ! we are an IO daemon
        call afdx(MY_GID,-1*nhst)
        call pvmfexit(info)
        return
     end if

  end if

  deallocate( SLVID, HOSTSL, PIDSL )
  deallocate( HOSTNAMES )
  deallocate( NPROCH, NCPUH )
  deallocate( IDH )
  deallocate( ICPUMASK, ICPU)
  deallocate( MASKSL )
  call pvmfexit(info)

             ! dissolve virtual machine if we started it

  call tstival('pvmd',pvmd)
  if (pvmd.eq.1) call pvmfhalt(info)

end subroutine para_driver

!=================================================

subroutine start_slaves

  use pvm_defs
  use newpara
  use messages
  
  implicit none
  
  character*256 txprog,arch,wline
  integer*4 ncpu,narch,dtid,speed,gsize,info
  integer img,i

               ! get the name of the program running

  call get_command_argument(0,txprog)
  
               ! use no. of slaves from the command line, if available,
               ! othervise set no. of slaves equal to the number of hosts
               ! in the virtual machine

  if(NSLV.eq.0) then 
    call addbuf(' One slave on each host in the virtual machine.')
    NSLV = PVM_NH
    call setival('nslv',NSLV)
  endif

  if(NSLV.le.0) then
    call nerror(1,'Parallel Startup',                             &
                'No of slaves is 0: please use the serial version', &
                NSLV,0)
  endif

               ! ...............................................
               !  get the hostname of the other machines in PVM
               !  start a slave on every machine
               !  start slave on master LAST
               !  if anything goes wrong: die!
               ! ...............................................

  allocate( SLVID(NSLV+1), HOSTSL(NSLV) )

  do i=1,NSLV
    call pvmfconfig(ncpu,narch,dtid,HOSTSL(i),arch,speed,info)
    call pvmfspawn(txprog,PvmHost,HOSTSL(i),1,SLVID(i+1),info)
    if (info.ne.1) then
      if (info.eq.0) then
        if (SLVID(i+1).eq.PvmNoFile) then
            call addbuf('Your executable is: '//trim(txprog))
            call addbuf('It was not found on '//trim(HOSTSL(i)))
            call addbuf('Make sure you start it with its full path')
            call addbuf('or it is in the PVM binary dir')
            call addbuf('or its path is in the PVM hostfile with ep=path')
        else
           call pvmfperror('Slave startup ',SLVID(i+1))
        end if
      else
         call pvmfperror(HOSTSL(i),info)
      end if
      call setival('nslv',i-1)
      call para_stop
      stop
    end if
  end do

              ! request a message if any of the slaves dies

  if (NSLV.gt.0) then
    call pvmfnotify(PvmTaskExit,TxFailure,NSLV,SLVID(2),info)
  end if
  call setival('nslv',NSLV)

              ! send out no. of slaves, needed for freeze         

  do i=1,NSLV
     call para_send(NSLV,i,TxNSlave)
  end do
  gsize=NSLV+1
  call pvmffreezegroup(GROUP_NAME,gsize,info)

end subroutine start_slaves

!======================================================================      

subroutine para_stop
!
! This is the routine that cleans up the slaves
! This should be called before stopping
!
  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info, ibuf, mid4, sigterm
  parameter (sigterm=9)
  integer i,pvmd,itest,errid
  character*256 errm1,errm2

  ! pack the termination message

  call para_initsend
  call para_pack_int(MY_GID,1)
  call tstchval('errm1',itest)
  if(itest.eq.1)then
    call getchval('errm1',errm1)
  else
    errm1=''
  endif
  call para_pack_string(errm1,256)
  call tstchval('errm2',itest)
  if(itest.eq.1)then
    call getchval('errm2',errm2)
  else
    errm2=''
  endif
  call para_pack_string(errm2,256)

  ! Broadcast the "kiss of death"

  call pvmfbcast( GROUP_NAME,TxKissofDeath, info )

  call pvmfexit(info)

  ! if we are the master, dissolve virtual machine if we started it

  if(MY_GID.eq.0) then
    call tstival('pvmd',pvmd)
    if (pvmd.eq.1) call pvmfhalt(info)
  endif

end subroutine para_stop

!======================================================================      

subroutine master_die

! this routine is called when the master receives
! the "kiss of death" message.
! It prints out the termination message then dies

  implicit none

  integer errid
  character*256 errm1,errm2

  call para_unpack_int( errid, 1 )
  call para_unpack_string( errm1, 256 )
  call para_unpack_string( errm2, 256 )
  write(6,100)errid,errm1(1:len_trim(errm1)),errm2(1:len_trim(errm2))
  stop '*** PQS Error Termination ***'

100 format(/'Received termination message from slave',i3':'/a/a/)

end subroutine master_die

!======================================================================      

subroutine slave_die

! this routine is called when a slave receives
! the "kiss of death" message.

  implicit none

  stop '*** PQS Error Termination ***'

end subroutine slave_die

!======================================================================      

subroutine para_check

! checks on the slaves (mostly useful while debugging)      

  use pvm_defs
  use newpara
  implicit none
  
  integer i
  integer*4 ibuf,tid,info
  
  call pvmfprobe(-1,TxFailure,ibuf)
  if(ibuf.ne.0) then
     call pvmfrecv(-1,TxFailure,ibuf)
     call pvmfunpack(INTEGER4,tid,1,1,info)
     do i=1,NSLV
        if (tid.eq.SLVID(i+1)) write(6,*) 'Slave died on ', &
                                   HOSTSL(i),'tid=',tid
     end do
     call para_stop
     stop 'Slave failure'
  end if

end subroutine para_check

!======================================================================      

subroutine set_pvm_def_int

   ! sort out the default integer type

  use kinds
  use pvm_defs
  use newpara
  implicit none

  select case ( kind(0) )

    case ( i_1 )

      INTEGER_DEF = BYTE1

    case ( i_2 )

      INTEGER_DEF = INTEGER2

    case ( i_4 )

      INTEGER_DEF = INTEGER4

    case ( i_8 )

      INTEGER_DEF = INTEGER8

    case default

      write(6,*)'set_pvm_def_int: unknown integer kind'
      call para_stop

  end select

  return

end subroutine set_pvm_def_int

!======================================================================      

subroutine para_initsend

! initialize buffers for packing
! do not allocate and use new memory while packing      
! only use while setting up (for compatibility with MPI
! send big chunks without packing      

  use pvm_defs
  implicit none
  
  integer*4 info
  
  call pvmfinitsend(PvmDataRaw,info)      

end subroutine para_initsend

!======================================================================      

subroutine para_pack_int(inte,inum)

! pack integer data (PVM)
! data in inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer inte,inum
  
  call pvmfpack(INTEGER_DEF,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_pack_int

!======================================================================      

subroutine para_unpack_int(inte,inum)

! unpack integer data (PVM)
! data to inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer inte,inum
  
  call pvmfunpack(INTEGER_DEF,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_unpack_int

!======================================================================      

subroutine para_pack_int4(inte,inum)

! pack integer data (PVM)
! data in inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 inte
  integer inum
  
  call pvmfpack(INTEGER4,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_pack_int4

!======================================================================      

subroutine para_unpack_int4(inte,inum)

! unpack integer data (PVM)
! data to inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 inte
  integer inum
  
  call pvmfunpack(INTEGER4,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_unpack_int4

!======================================================================      

subroutine para_pack_int2(inte,inum)

! pack integer data (PVM)
! data in inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*2 inte
  integer inum
  
  call pvmfpack(INTEGER2,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_pack_int2

!======================================================================      

subroutine para_unpack_int2(inte,inum)

! unpack integer data (PVM)
! data to inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*2 inte
  integer inum
  
  call pvmfunpack(INTEGER2,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_unpack_int2

!======================================================================      

subroutine para_pack_real(dat,inum)

! pack integer data (PVM)
! data in inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  real*8 dat
  integer inum
  
  call pvmfpack(REAL8,dat,inum,1,info)
!  if (info.lt.0) call sleep(600)

end subroutine para_pack_real

!======================================================================      

subroutine para_unpack_real(dat,inum)

! unpack integer data (PVM)
! data to inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  real*8 dat
  integer inum
  
  call pvmfunpack(REAL8,dat,inum,1,info)
!  if (info.lt.0) call sleep(600)
  
  end subroutine para_unpack_real

!======================================================================      

subroutine para_pack_string(stri,inum)

! pack string (PVM), max. length of string is 256
! data in string 
! length in inum

  use pvm_defs
  implicit none
  
  integer*4 info
  integer inum
  character*256 stri
  
  call pvmfpack(STRING,stri,inum,1,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_pack_string

!======================================================================      

subroutine para_unpack_string(stri,inum)

! unpack string (PVM), max. length of string is 256
! data to string 
! length in inum

  use pvm_defs
  implicit none
  
  integer*4 info
  integer inum
  character*256 stri
  
  call pvmfunpack(STRING,stri,inum,1,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_unpack_string
    
!======================================================================      

subroutine para_pack_byte(inte,inum)

! pack byte data (PVM)
! data in inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  integer inum
  integer*1 inte
  
  call pvmfpack(BYTE1,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_pack_byte
    
!======================================================================      

subroutine para_unpack_byte(inte,inum)

! unpack byte data (PVM)
! data to inte (single variable or first member of array)
! array length in inum (1 for single, no bound checking!)

  use pvm_defs
  implicit none
  
  integer*4 info
  integer inum
  integer*1 inte
  
  call pvmfunpack(BYTE1,inte,inum,1,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_unpack_byte
    
!======================================================================      

subroutine para_send(inte,host,tag)

! send integer to host with message tag

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag
  integer host,inte
  
  call pvmfpsend(SLVID(host+1),tag,inte,1,INTEGER_DEF,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_send

!======================================================================      

subroutine para_recv(inte,ifrom,tag)

! receive integer from any host with message tag
! source is returned in ifrom

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer*4 tag,atid,atag,anum,i4from
  integer inte,ifrom
  
  call pvmfprecv(-1,tag,inte,1,INTEGER_DEF,atid,atag,anum,info)
!  if (info.lt.0) call sleep(600)
  call pvmfgetinst(GROUP_NAME,atid,i4from)
  ifrom=i4from
    
end subroutine para_recv

!======================================================================      

subroutine para_send_real(dat,inum,host,tag)

! send real(array) to host with message tag

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag
  integer host,inum
  real*8 dat
  
  call pvmfpsend(SLVID(host+1),tag,dat,inum,REAL8,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_send_real
    
!======================================================================      

subroutine para_recv_real(dat,inum,host,tag)

! receive real(array) from  host with message tag
! !!not from any host like the other recv routines!!)

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag,atid,atag,anum
  integer inum,host
  real*8 dat
  
  call pvmfprecv(SLVID(host+1),tag,dat,inum,REAL8, &
                 atid,atag,anum,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_recv_real

!======================================================================      

subroutine para_send_pack(host,tag)

! send packed data to host with message tag

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer host
  integer*4 tag
  
  call pvmfsend(SLVID(host+1),tag,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_send_pack

!======================================================================      

subroutine para_recv_pack(ifrom,tag)

! receive packed data with message tag from _ANYBODY_
! source id is _returned_ in ifrom      

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag,ibuf,ilen,itag,isrc,i4from
  integer ifrom
  
  call pvmfrecv(-1,tag,ibuf)
!  if (ibuf.le.0) call sleep(600)
  call pvmfbufinfo(ibuf,ilen,itag,isrc,info)
  call pvmfgetinst(GROUP_NAME,isrc,i4from)
  ifrom=i4from
!  if (info.lt.0) call sleep(600)
    
end subroutine para_recv_pack

!======================================================================      

subroutine para_bcast(inte,tag)

! broadcast int to all slaves with message tag

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer*4 tag
  integer inte
  
  call pvmfinitsend(PvmDataInPlace,info)
  call pvmfpack(INTEGER_DEF,inte,1,1,info)
  call pvmfbcast(GROUP_NAME,tag,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_bcast

!======================================================================      

subroutine para_recv_bcast(inte,tag)

! receive broadcast integer with message tag from master host

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag,atag,atid,anum,masterid4
  integer inte
  
  masterid4=MASTER_ID
  call pvmfprecv(masterid4,tag,inte,1,INTEGER_DEF, &
                 atid,atag,anum,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_recv_bcast

!======================================================================      

subroutine para_bcast_real(dat,inum,tag)

! broadcast real (array) to all slaves with message tag
!      
! this is just a short wrapper to chop up really big data blocks      
! if the array is bigger than a certain size
! the limit is 4,000,000 words, presently most linux kernels use
! this as maximum shared memory (cat /proc/sys/kernel/shmmax)
! MPI fails if larger chunks are sent on dual proc. machines
! PVM would benefit as well from the smaller buffers      

  use newpara
  implicit none

  integer inum,i,nchunk,ilen
  integer*4 tag
  real*8 dat(inum)

  nchunk=4000000
  do i=1,inum,nchunk
    ilen=min(inum-i+1,nchunk)
    call para_bcast_float(dat(i),ilen,tag)
  enddo

end subroutine para_bcast_real
        
!======================================================================      

subroutine para_recv_bcastreal(dat,inum,tag)

! receive broadcast real (array) with message tag from master host
! this is a wrapper to chop up big data blocks see above      

  use newpara
  implicit none

  integer inum,i,nchunk,ilen
  integer*4 tag
  real*8 dat(inum)

  nchunk=4000000
  do i=1,inum,nchunk
    ilen=min(inum-i+1,nchunk)
    call para_recv_bcastfloat(dat(i),ilen,tag)
  enddo

end subroutine para_recv_bcastreal

!======================================================================      

subroutine para_bcast_float(dat,inum,tag)

! broadcast using a binary tree
! master sends to slaves 1,2,4,8,.... (in reverse order)
! slaves distribute the data further       
! this is wrapped by para_bcast_real to cope with big data blocks

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info,tag
  integer nfact
  integer inum
  real*8 dat,tred,tafter,tbefore
  
  call getrval('bcast',tred)
  call elapsec(tbefore)
  nfact=1
  do while (nfact.lt.NSLV)
      nfact=nfact*2
  end do
  do while (nfact.ge.1)
    if ((MY_GID+nfact).le.NSLV) then
      call pvmfpsend(SLVID(MY_GID+nfact+1),tag,dat, &
                     inum,REAL8,info)
    endif
    nfact=nfact/2
  enddo
  call elapsec(tafter)
  call setrval('bcast',tred+tafter-tbefore)
    
end subroutine para_bcast_float
    
!======================================================================      

subroutine para_recv_bcastfloat(dat,inum,tag)

! receive broadcast using a binary tree
! receive data (from master on another slave)
! send it to other slaves
! source and denstinations are decided according to gid (1,2,3,nslv)      
! this is wrapped by para_recv_bcastfloat to cope with big data blocks

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer*4 nfact,atid,tag,atag,anum,ito
  integer inum
  real*8 dat
  
  nfact=1
  do while ((MY_GID-nfact).gt.0)
    nfact=nfact*2
  end do
  do while (nfact.ge.1)
    ito=MY_GID+nfact
    if (mod(ito*1,nfact*2).eq.0) then
      call pvmfprecv(SLVID(MY_GID-nfact+1),tag,dat, &
                     inum,REAL8,atid,atag,anum,info)
    else if ((mod(ito,nfact).eq.0).and.(ito.le.NSLV)) then 
      call pvmfpsend(SLVID(ito+1),tag,dat,inum,REAL8,info)
    endif
    nfact=nfact/2
  enddo
    
end subroutine para_recv_bcastfloat

!======================================================================      

subroutine para_bcast_pack(tag)

! broadcast using a binary tree
! master sends to slaves 1,2,4,8,.... (in reverse order)
! slaves distribute the data further       

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer*4 nfact,tag
  real*8 tred,tafter,tbefore
  
  call getrval('bcast',tred)
  call elapsec(tbefore)
  nfact=1
  do while (nfact.lt.NSLV)
      nfact=nfact*2
  end do
  do while (nfact.ge.1)
      if ((MY_GID+nfact).le.NSLV) then
         call pvmfsend(SLVID(MY_GID+nfact+1),tag,info)
      endif
      nfact=nfact/2
  enddo
  call elapsec(tafter)
  call setrval('bcast',tred+tafter-tbefore)
    
end subroutine para_bcast_pack
    
!====================================================================== 

subroutine para_barrier

  use pvm_defs
  use newpara
  implicit none
  integer*4 ierr,igrpsz

  igrpsz=NSLV+1
  call pvmfbarrier(GROUP_NAME,igrpsz,ierr)

end subroutine para_barrier

!======================================================================      

subroutine para_recv_bcastpack(tag)

! receive broadcast using a binary tree
! receive data (from master on another slave)
! send it to other slaves
! source and denstinations are decided according to gid (1,2,3,nslv)      

  use pvm_defs
  use newpara
  implicit none
  
  integer*4 info
  integer*4 nfact,tag,ibuf,oldbuf,ito
  
  nfact=1
  do while ((MY_GID-nfact).gt.0)
      nfact=nfact*2
  end do
  do while (nfact.ge.1)
      ito=MY_GID+nfact
      if (mod(ito*1,nfact*2).eq.0) then
         call pvmfrecv(SLVID(MY_GID-nfact+1),tag,ibuf)
      else if ((mod(ito,nfact).eq.0).and.(ito.le.NSLV)) then 
         call pvmfsetsbuf(ibuf,oldbuf)
         call pvmffreebuf(oldbuf,info) ! need to free the old buffer
         call pvmfsend(SLVID(ito+1),tag,info)
         call pvmfsetrbuf(ibuf, oldbuf)
         call pvmffreebuf(oldbuf,info) ! need to free the old buffer
      endif
      nfact=nfact/2
  enddo
    
end subroutine para_recv_bcastpack

!======================================================================      

subroutine para_reduce(dat,inum,tag)

! reduce data block
! this is a wrapper to chop up big data blocks
! on PVM this seems to be useful, big chunks seem to slow down trfer

  use newpara
  implicit none

  integer inum,i,nchunk,ilen,idum,ifrom
  integer*4 tag
  real*8 dat(inum)

  nchunk=4000000
  do i=1,inum,nchunk
     ilen=min(inum-i+1,nchunk)
     call para_reduce_block(dat(i),ilen,tag)
  enddo
        
end subroutine para_reduce

!======================================================================      

subroutine para_reduce_block(dat,inum,tag)

! sum up real (array) dat from all slaves on master
! binary tree version
! wrapped by para_reduce for big data blocks

  use pvm_defs
  use newpara
  use memory
  implicit none
  
  integer*4 info
  integer*4 tag,atid,atag,anum,nfact
  integer inum,ibfp
  real*8 dat,tred,tafter,tbefore
  
  call getmem(inum,ibfp)
  if (MY_GID.eq.0) then ! time on master
     call para_check
     call getrval('reduce',tred)
     call elapsec(tbefore)
  endif
  nfact=1
  do while (nfact.le.NSLV)
    if (mod(MY_GID*1,nfact*2).eq.0) then
      if (MY_GID+nfact.le.NSLV) then
        call pvmfprecv(SLVID(MY_GID+nfact+1),tag,bl(ibfp), &
                       inum,REAL8,atid,atag,anum,info)
        call daxpy(inum,1.0d0,bl(ibfp),1,dat,1)
      end if
    else 
      call pvmfpsend(SLVID(MY_GID-nfact+1),tag,dat, &
                     inum,REAL8,info)
      nfact=NSLV+1
    endif
    nfact=nfact*2
  end do
  call retmem(1)
  if (MY_GID.eq.0) then
     call elapsec(tafter)
     call setrval('reduce',tred+tafter-tbefore)
  endif
    
end subroutine para_reduce_block

!======================================================================      

subroutine para_probe(tag,isthere)

! isthere returns true if a message correponding to tag has 
! been received

  use pvm_defs
  implicit none
  
  integer*4 tag,bufid
  logical isthere
  
  call pvmfprobe(-1,tag,bufid)
  isthere=(bufid.gt.0)
    
end subroutine para_probe

!======================================================================      

subroutine para_mem(lcore,incore,indisk)

!
! Send memory info to slaves 
! send also the location of PQS_ROOT and scratch file root
!
  use newpara
  implicit none
  
  integer lcore,incore,indisk
  character*256 pqs_root,scrfile,scrdir,jobname
  
  call para_check
  call getchval('pqs_root',pqs_root)
  call getchval('scrf',scrfile)
  call getchval('scrdir',scrdir)
  call getchval('jobname',jobname)
  call para_initsend
  call para_pack_int(lcore,1)
  call para_pack_int(incore,1)
  call para_pack_int(indisk,1)
  call para_pack_string(pqs_root,256)
  call para_pack_string(scrfile,256)
  call para_pack_string(scrdir,256)
  call para_pack_string(jobname,256)
  call para_bcast_pack(TxMemSize)

end subroutine para_mem

!======================================================================      

subroutine para_slave_mem(lcore,incore,indisk,pqs_root,scrfile,scrdir,jobname)

!
! Get memory info (bl is not set up yet, need buffer)
!
  use newpara
  implicit none
  
  integer lcore,incore,indisk
  character*256 pqs_root,scrfile,scrdir,jobname
  
  call para_recv_bcastpack(TxMemSize)
  call para_unpack_int(lcore,1)
  call para_unpack_int(incore,1)
  call para_unpack_int(indisk,1)
  call para_unpack_string(pqs_root,256)
  call para_unpack_string(scrfile,256)
  call para_unpack_string(scrdir,256)
  call para_unpack_string(jobname,256)

end subroutine para_slave_mem

! native pvm version of calls implemented with binary tree
! pvm sends messages one by one internally


!======================================================================      

subroutine para_bcast_float_flat(dat,inum,tag)

! broadcast real (array) to all slaves with message tag
! this is better than packing in and out for large chunks      
! probably using pvm_psend and binary tree would be even better

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag
  integer inum
  real*8 dat
  real*8 tred,tbefore,tafter
  
  call getrval('bcast',tred)
  call elapsec(tbefore)
  call pvmfinitsend(PvmDataInPlace,info)
  call pvmfpack(REAL8,dat,inum,1,info)
!  if (info.lt.0) call sleep(600)
  call pvmfbcast(GROUP_NAME,tag,info)
!  if (info.lt.0) call sleep(600)
  call elapsec(tafter)
  call setrval('bcast',tred+tafter-tbefore)
    
end subroutine para_bcast_float_flat

!======================================================================      

subroutine para_recv_bcastfloat_flat(dat,inum,tag)

! receive broadcast real (array) with message tag from master host

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag,atag,atid,anum,masterid4
  integer inum
  real*8 dat
  
  masterid4=MASTER_ID
  call pvmfprecv(masterid4,tag,dat,inum,REAL8,atid,atag,anum,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_recv_bcastfloat_flat

!======================================================================      

subroutine para_bcast_pack_flat(tag)

! broadcast packed data to all slaves with message tag

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag
  real*8 tred,tafter,tbefore
  
  call getrval('bcast',tred)
  call elapsec(tbefore)
  call pvmfbcast(GROUP_NAME,tag,info)
  call elapsec(tafter)
  call setrval('bcast',tred+tafter-tbefore)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_bcast_pack_flat

!======================================================================      

subroutine para_recv_bcastpack_flat(tag)

! receive broadcast packed data with message tag from master host

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag
  
  call pvmfrecv(MASTER_ID,tag,info)
!  if (info.lt.0) call sleep(600)
    
end subroutine para_recv_bcastpack_flat

!======================================================================      

subroutine para_reduce_block_flat(dat,inum,tag)

! sum up real (array) dat from all slaves on master
! binary tree with summing and buffering proved better on PVM

  use newpara
  use pvm_defs
  implicit none
  
  external PvmSum
  integer*4 info
  integer*4 tag
  integer inum
  real*8 dat,tred,tafter,tbefore
  
  if (MY_GID.eq.0) then ! time on master
     call para_check
     call getrval('reduce',tred)
     call elapsec(tbefore)
  endif
  call pvmfreduce(PvmSum,dat,inum,REAL8,tag,GROUP_NAME,0,info)
!  if (info.lt.0) call sleep(600)
  if (MY_GID.eq.0) then
     call elapsec(tafter)
     call setrval('reduce',tred+tafter-tbefore)
  endif
    
end subroutine para_reduce_block_flat

!======================================================================

subroutine para_printhelp(ichan)
!
!  prints an help message (parallel MPI version)
!
  implicit none
  integer ichan
  write(ichan,100)
 100  format(/' USAGE: pqs_pvm.x [-v,-h] [-i] molecule-name'/)
end subroutine para_printhelp

!======================================================================

subroutine para_sendbin( ibins,       ibin1,   lbin,  icount,   &
                         isize,       igran,   idest, sendb,    &
                         igranulesize,ngran,   imod,  ndisk2,   &
                         nbin,        nrec,    islv)

  use newpara
  implicit none
!
!  Send the bin to another slave, PVM version.
!
!  This PVM version has noting special, just standard calls to the
!  newpara routines, because the PVM send is non-blocking by default
!
!  ARGUMENTS
!
!  ibins   -  integer storage for the bins
!  ibin1   -  precision overflow integer
!  lbin    -  kength of one bin
!  icount  -  number of non-zero elements in IJ bin
!  isize   -  size of present granule (in "ij units")
!  igran   -  current granule number
!  idest   -  destination slave
!  sendb   -  buffer for non-blocking send (used only by MPI version)
!  igranulesize - parameter for call to para_recv_bin (used only by MPI version)
!  ngran   - ditto
!  imo     - ditto
!  ndisk2  - ditto
!  nbin    - ditto
!  nrec,   - ditto
!  islv    - ditto
!
  integer lbin,isize,igran,idest
  INTEGER*4 ibins(2*lbin)
  INTEGER*1 ibin1(lbin)
  integer*4 icount(isize)
  real*8 sendb
  integer igranulesize,ngran,imod,ndisk2,nbin,nrec,islv
!
! pack bin and send it
!
  call para_initsend
  call para_pack_int(igran,1)       ! granule number
  call para_pack_int(isize,1)       ! granule size
  if(isize.gt.0)then
    call para_pack_int4(icount,isize) ! number of elements in bin
    call para_pack_int4(ibins,2*lbin*isize)  ! bin itself
    call para_pack_byte(ibin1,lbin*isize)  ! precision overflow
  endif
  call para_send_pack(idest,TxBinS1)

  return
end subroutine para_sendbin

!======================================================================

subroutine para_directroute

! enable diretc routing between processes

  use pvm_defs
  use newpara
  implicit none
  integer*4 info

  call pvmfsetopt(PvmRoute,PvmRouteDirect,info)

end subroutine para_directroute

!======================================================================

subroutine para_nodirectroute

! disable direct routing between processes.
! some small MP2 jobs show random errors if direct routing
! is enabled during the parallel binsort

  use pvm_defs
  use newpara
  implicit none
  integer*4 info

  call pvmfsetopt(PvmRoute,PvmDontRoute,info)

end subroutine para_nodirectroute
!======================================================================      

subroutine para_reduce_max(dat,inum,tag)

! gather max values of real (array) dat from all slaves on master
! binary tree version

  use pvm_defs
  use newpara
  use memory
  implicit none
  
  integer*4 info
  integer*4 tag,atid,atag,anum,nfact
  integer inum,ibfp,idt
  real*8 dat(inum),tred,tafter,tbefore
  
  call getmem(inum,ibfp)
  if (MY_GID.eq.0) then ! time on master
     call para_check
     call getrval('reduce',tred)
     call elapsec(tbefore)
  endif
  nfact=1
  do while (nfact.le.NSLV)
    if (mod(MY_GID*1,nfact*2).eq.0) then
      if (MY_GID+nfact.le.NSLV) then
        call pvmfprecv(SLVID(MY_GID+nfact+1),tag,bl(ibfp), &
                       inum,REAL8,atid,atag,anum,info)
        do idt=1,inum
           dat(idt)=max(dat(idt),bl(ibfp+idt-1))
        enddo
      end if
    else 
      call pvmfpsend(SLVID(MY_GID-nfact+1),tag,dat, &
                     inum,REAL8,info)
      nfact=NSLV+1
    endif
    nfact=nfact*2
  end do
  call retmem(1)
  if (MY_GID.eq.0) then
     call elapsec(tafter)
     call setrval('reduce',tred+tafter-tbefore)
  endif
    
end subroutine para_reduce_max

!======================================================================      

subroutine para_send_int(dat,inum,host,tag)

! send int(array) to host with message tag

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag
  integer host,inum
  integer dat
  
  call pvmfpsend(SLVID(host+1),tag,dat,inum,INTEGER_DEF,info)
    
end subroutine para_send_int
    
!======================================================================      

subroutine para_recv_int(dat,inum,host,tag)

! receive int(array) from  host with message tag

  use newpara
  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag,atid,atag,anum
  integer inum,host
  integer dat
  
  call pvmfprecv(SLVID(host+1),tag,dat,inum,INTEGER_DEF, &
                 atid,atag,anum,info)
    
end subroutine para_recv_int

!======================================================================      

subroutine para_recv_anypack(ifrom,tag)

! receive packed data with any message tag from _ANYBODY_
! source id is _returned_ in ifrom      
! tag is returned in tag

  use pvm_defs
  implicit none
  
  integer*4 info
  integer*4 tag,ibuf,ilen,itag,isrc,i4from
  integer ifrom
  
  call pvmfrecv(-1,-1,ibuf)
  call pvmfbufinfo(ibuf,ilen,itag,isrc,info)
  call pvmfgetinst(GROUP_NAME,isrc,i4from)
  ifrom=i4from
  tag=itag
    
end subroutine para_recv_anypack

!======================================================================      

  subroutine fafinit(iret)  ! collective!
! starts I/O daemons on master and transfers info to slaves
  use newpara
  use pvm_defs
  implicit none
  integer, intent(out) :: iret !are we I/O?
  integer ii
  integer*4 info

  allocate( DMNID(PVM_NH) )
  
  if (MY_GID .eq. 0) then
    call start_daemons
             ! broadcast the daemon tid database
    call para_initsend
    call para_pack_int4(DMNID,PVM_NH)
    call para_bcast_pack(TxMask)
  !  do ii=1,PVM_NH  !send daemon tids to all the daemons
  !    call pvmfsend(DMNID(ii),12345,info)
  !  enddo
  else    !on slaves
              ! receive the daemon tid database
    call para_recv_bcastpack(TxMask)
    call para_unpack_int4(DMNID,PVM_NH)
  endif
  call afinit(MY_GID,PVM_NH,DMNID,MASTER_ID,NSLV,SLVID)
  iret=0
  end subroutine fafinit

!======================================================================      

  subroutine fafterminate(info)
  use newpara
  use pvm_defs
  implicit none
  integer info
  call afterminate(info)
  deallocate(DMNID)
  end

!======================================================================

subroutine start_daemons

  use pvm_defs
  use newpara
  use messages
  
  implicit none
  
  character*256 txprog,arch,hostnm
  integer*4 ncpu,narch,dtid,speed,gsize,info
  integer i

               ! get the name of the program running

  call get_command_argument(0,txprog)

               ! ...............................................
               !  get the hostname of the other machines in PVM
               !  start a slave on every machine
               !  start slave on master LAST
               !  if anything goes wrong: die!
               ! ...............................................

  do i=1,PVM_NH
    call pvmfconfig(ncpu,narch,dtid,hostnm,arch,speed,info)
    call pvmfspawn(txprog,PvmHost,hostnm,1,DMNID(i),info)
    if (info.ne.1) then
      if (info.eq.0) then
        if (DMNID(i).eq.PvmNoFile) then
            call addbuf('Your executable is: '//trim(txprog))
            call addbuf('It was not found on '//trim(hostnm))
            call addbuf('Make sure you start it with its full path')
            call addbuf('or it is in the PVM binary dir')
            call addbuf('or its path is in the PVM hostfile with ep=path')
        else
           call pvmfperror('Slave startup ',DMNID(i))
        end if
      else
         call pvmfperror(HOSTSL(i),info)
      end if
      call para_stop
      stop
    end if
  end do

              ! send out negative no. of slaves, as a signal of being IO daemons

  do i=1,PVM_NH
     call pvmfpsend(DMNID(i),TxNSlave,-1*PVM_NH,1,INTEGER_DEF,info)
  end do

end subroutine start_daemons

!======================================================================      
      subroutine fafreopen
      end
!======================================================================
      subroutine sync_barrier
      end
!======================================================================
      subroutine array_files
      end
!======================================================================
      subroutine fafnbread
      end
!======================================================================
      subroutine fafrequest
      end
!======================================================================
      subroutine block_alarm
      end
!======================================================================
      subroutine unblock_alarm
      end
!======================================================================
