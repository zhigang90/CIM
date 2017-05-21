      subroutine para_driver

c This file contains the MPI specific parallel routines:
c
c para_driver:    startup and forking between master and slave
c start_slaves    (function obvious)
c para_stop:      cleanup slaves in case of nonregular stop
c para_check:     check the slaves - and stop if any of them died

c set_mpi_def_int set the default integer type for MPI

c para_initsend:  initialize buffers
c para_pack_....  pack data(int,real,string,byte)      
c para_unpack_... unpack data(int,real,string,byte)      

c para_send_....  send msg w/ tag to host specified
c para_recv_....  recv msg w/ tag from any host, 
c                          source host is returned
c!!! para_recv_real receives only from host specified 

c para_bcast_.... msg w/ tag from master to all slaves
c para_recv_bcast.... recv msg w/ tag from master

c para_....()     msg is a single integer
c para_...._pack  msg is packed
c para_...._real  msg is a real array

c para_reduce     sum up array on master
c para_probe      check if message has arrived
c para_mem, para_slave_mem         send/receive memory limits

c para_..._..._flat is a version of the routine not using a 
c                   binary tree, but the slower native pvm calls
c
c para_sendbin    send bin to another slave during MP2 binsort


c here follows the code for para_driver:

c initialization, forking between slave and master      

      use newpara
      use messages
      implicit none
      include 'mpif.h'

      integer*4 info,mytid,itid,mygid,ncpu,reslen,ipos,ncpu4,pidslv
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 get_pid,getaffinity,setaffinity, setbatch, ierr_handler
      integer img,i,ii,j,nhst,idslv,ncpus,ncpuslv
      logical ishost
      character*256 slname
      character*256 hstnm
      character*(MPI_MAX_PROCESSOR_NAME) mpi_hstnm
      character*256 progname,inpname,prodstr,wline
      character buffer(2048)
      external mpi_errorh
 
c  set product string for license checking (override prodstr with parastr)
 
      call getchval('parastr',prodstr)
      call setchval('prodstr',prodstr)

c validate the program arguments

      call parsearguments(progname,inpname)
        
      call mpi_init(info)

      call set_mpi_def_int   ! set mpi type of default integer
      NEW_WORLD=MPI_COMM_WORLD

#ifdef DEBUG
      call mpi_errhandler_create( mpi_errorh, ierr_handler, info )
      call mpi_errhandler_set(NEW_WORLD,ierr_handler,info)
#endif

c set MY_GID and NSLV

      call mpi_comm_rank(NEW_WORLD,mygid,info)
      call mpi_comm_size(NEW_WORLD,ncpu,info)
      MY_GID=mygid
      NSLV=ncpu-1
      call setival('nslv',NSLV)
      MY_PID = get_pid()

      allocate( SLVID(NSLV+1), HOSTSL(NSLV) )

c set the message passing type

      call setchval('msgpass','MPI')

c get my hostname (put in module data block)

      call blankit(hstnm,256)
      call mpi_get_processor_name(mpi_hstnm,reslen,info)
      hstnm=mpi_hstnm(1:min(len_trim(mpi_hstnm),256))
      MY_HOSTNM=hstnm
c     write (6,*), MY_GID,hstnm
c     call f_lush(6)

      wbuf=''          ! zero out message buffer

      if (mygid.eq.0) then

c I am the master
c gather id, names of slaves and number of CPU per host

        allocate( HOSTNAMES(NSLV+1), PIDSL(NSLV) )
        allocate( NPROCH(NSLV+1), NCPUH(NSLV+1) )
        NPROCH=0
        nhst=1
        HOSTNAMES(nhst)=MY_HOSTNM    ! the first host is the master's
        call getival('ncpus',ncpus)
        ncpu4=ncpus
        NCPUH(nhst)=ncpus
        NPROCH(nhst)=NPROCH(nhst)+1
        call mpi_pack_size(2048,MPI_CHARACTER,NEW_WORLD,
     &                    MPI_IBUF,info)
        do i=1,NSLV
          slname=''
          call mpi_recv(buffer,MPI_IBUF,MPI_PACKED,MPI_ANY_SOURCE, 
     &                  TxId,NEW_WORLD,stat,info)
          ipos=0
          call mpi_unpack(buffer,MPI_IBUF,ipos,idslv,1, 
     &               INTEGER_DEF,NEW_WORLD,info)
          call mpi_unpack(buffer,MPI_IBUF,ipos,ncpuslv,1, 
     &               INTEGER_DEF,NEW_WORLD,info)
          call mpi_unpack(buffer,MPI_IBUF,ipos,pidslv,1, 
     &               MPI_INTEGER4,NEW_WORLD,info)
          call mpi_unpack(buffer,MPI_IBUF,ipos,slname,256, 
     &                MPI_CHARACTER,NEW_WORLD,info)
          HOSTSL(idslv)=slname
          PIDSL(idslv)=pidslv
          do ii=1,nhst
            if(slname.eq.HOSTNAMES(ii))then
              NPROCH(ii)=NPROCH(ii)+1
              exit
            endif
          enddo
          if(ii.gt.nhst)then
            nhst=nhst+1
            HOSTNAMES(nhst)=slname
            NCPUH(nhst)=ncpuslv
            NPROCH(nhst)=NPROCH(nhst)+1
          endif
        enddo

c  number of hosts, maximum number of processes and CPU per host

        NHOST=nhst
        MAXPROCH=NPROCH(1)
        MAXCPUH=NCPUH(1)
        do i=1,NHOST
          if(NPROCH(i).gt.MAXPROCH)MAXPROCH=NPROCH(i)
          if(NCPUH(i).gt.MAXCPUH)MAXCPUH=NCPUH(i)
        enddo

c  fill in the database of process id per host

        allocate( IDH(NHOST,MAXPROCH) )
        IDH=-1
        MY_HOSTID=1
        IDH(1,1)=0   ! the master is the first process of the first host
        do i=1,NSLV
           slname=HOSTSL(i)
           do ii =1,NHOST
             if(slname.eq.HOSTNAMES(ii))then
               do j=1,MAXPROCH
                 if(IDH(ii,j).eq.-1)then
                   IDH(ii,j)=i
                   exit
                 endif
               enddo
               exit
             endif
           enddo
        enddo

c  broadcast the host database

        call MPI_Bcast(NHOST,1,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(MAXPROCH,1,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(MAXCPUH,1,INTEGER_DEF,0,NEW_WORLD,info)

        do i=1,NHOST
          call MPI_Bcast(HOSTNAMES(i),256,MPI_CHARACTER,0,
     &                 NEW_WORLD,info)
        enddo
        call MPI_Bcast(NPROCH,NHOST,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(NCPUH,NHOST,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(IDH,NHOST*MAXPROCH,INTEGER_DEF,0,
     &                 NEW_WORLD,info)
        call MPI_Bcast(PIDSL,NSLV,MPI_INTEGER4,0,NEW_WORLD,info)

c  initialize cpu affinity mask

        allocate( ICPUMASK(ncpus), ICPU(ncpus) )
        allocate( MASKSL(MAXCPUH,NSLV) )
        MASKSL=0
        ICPUMASK=0
        info = getaffinity( ncpu4, ICPUMASK )
        ICPU=ICPUMASK

c  set cpu affinity

        if(info.eq.0)then
          call para_affinity(MY_GID,MAXPROCH,ncpus,
     $                       IDH(MY_HOSTID,1),ICPU)
          info = setaffinity( ncpu4, ICPU )
          if(info.ne.0)then
            call addbuf('')
            call addbuf(' *** Warning: Master cannot set cpu affinity')
          endif
        else
          call addbuf('')
          call addbuf(' *** Warning: Master cannot set cpu affinity')
        endif
#ifdef SPIN

c  set scheduling policy 

        info = setbatch()
        if(info.ne.0)then
          call addbuf('')
          call addbuf(' *** Warning: Master cannot set scheduler')
        endif
#endif

c    gather affinity masks from slaves

        call mpi_pack_size(2048,MPI_CHARACTER,NEW_WORLD,
     &                    MPI_IBUF,info)
        do i=1,NSLV
          call mpi_recv(buffer,MPI_IBUF,MPI_PACKED,MPI_ANY_SOURCE, 
     &                  TxMask,NEW_WORLD,stat,info)
          ipos=0
          call mpi_unpack(buffer,MPI_IBUF,ipos,idslv,1, 
     &               INTEGER_DEF,NEW_WORLD,info)
          call mpi_unpack(buffer,MPI_IBUF,ipos,ncpuslv,1, 
     &               INTEGER_DEF,NEW_WORLD,info)
          call mpi_unpack(buffer,MPI_IBUF,ipos,MASKSL(1,idslv),ncpuslv, 
     &               MPI_INTEGER4,NEW_WORLD,info)
        enddo

c   broadcast the mask database

        call MPI_Bcast(MASKSL,NSLV*MAXCPUH,MPI_INTEGER4,0,
     &                 NEW_WORLD,info)

c compose the startup message

        call para_startmsg

c zero out the counter for slave's cpu time

        call setrval('cpuscum',0.0d0)
        call setrval('cpus',0.0d0)
        call setrval('reduce',0.0d0)
        call setrval('bcast',0.0d0)
        call setrval('imbal',0.0d0)

c now read the input and do the work

        call driver

c stop slaves         

        call para_bcast(TxFinish,TxJobType)
      else

c I am a slave
c
c Send my slave id, number of cpus and hostname

        ipos=0
        call mpi_pack(MY_GID,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call getival('ncpus',ncpus)
        ncpu4=ncpus
        call mpi_pack(ncpus,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call mpi_pack(MY_PID,1,MPI_INTEGER4,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call mpi_pack(MY_HOSTNM,256,MPI_CHARACTER,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call mpi_send(buffer,ipos,MPI_PACKED,0,TxId,NEW_WORLD,info)

c receive the broadcast of the host database

        call MPI_Bcast(NHOST,1,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(MAXPROCH,1,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(MAXCPUH,1,INTEGER_DEF,0,NEW_WORLD,info)

        allocate( PIDSL(NSLV) )
        allocate( HOSTNAMES(NHOST) )
        allocate( NPROCH(NHOST), NCPUH(NHOST) )
        allocate( IDH(NHOST,MAXPROCH) )

        do i=1,NHOST
          call MPI_Bcast(HOSTNAMES(i),256,MPI_CHARACTER,0,
     &                 NEW_WORLD,info)
        enddo
        call MPI_Bcast(NPROCH,NHOST,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(NCPUH,NHOST,INTEGER_DEF,0,NEW_WORLD,info)
        call MPI_Bcast(IDH,NHOST*MAXPROCH,INTEGER_DEF,0,
     &                 NEW_WORLD,info)
        call MPI_Bcast(PIDSL,NSLV,MPI_INTEGER4,0,NEW_WORLD,info)

c compute MY_HOSTID

        do i=1,NHOST
          if(HOSTNAMES(i).eq.MY_HOSTNM)then
            MY_HOSTID=i
            exit
          endif
        enddo

c  initialize cpu affinity mask

        allocate( ICPUMASK(ncpus), ICPU(ncpus) )
        allocate( MASKSL(MAXCPUH,NSLV) )
        ICPUMASK=0
        MASKSL=0
        info = getaffinity( ncpu4, ICPUMASK )
        ICPU=ICPUMASK

c  set cpu affinity

        if(info.eq.0)then
          call para_affinity(MY_GID,MAXPROCH,ncpus,
     $                       IDH(MY_HOSTID,1),ICPU)
          info = setaffinity( ncpu4, ICPU )
          if(info.ne.0)then
            write(*,'(a,i0,a)')'*** Warning: Slave ',MY_GID,
     $              ' cannot set cpu affinity'
          endif
        else
          write(*,'(a,i0,a)')'*** Warning: Slave ',MY_GID,
     $            ' cannot set cpu affinity'
        endif
#ifdef SPIN

c  set scheduling policy 

        info = setbatch( )
        if(info.ne.0)then
          write(*,'(a,i0,a)')'*** Warning: Slave ',MY_GID,
     $            ' cannot set scheduler'
        endif
#endif

c  send affinity mask

        ipos=0
        call mpi_pack(MY_GID,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call mpi_pack(ncpus,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call mpi_pack(ICPU,ncpus,MPI_INTEGER4,buffer,2048,ipos, 
     &              NEW_WORLD,info)
        call mpi_send(buffer,ipos,MPI_PACKED,0,TxMask,
     &              NEW_WORLD,info)

c receive back the mask database

        call MPI_Bcast(MASKSL,NSLV*MAXCPUH,MPI_INTEGER4,0,
     &                 NEW_WORLD,info)

c  now do the work

        call slave

      end if

      deallocate( SLVID, HOSTSL, PIDSL )
      deallocate( HOSTNAMES )
      deallocate( NPROCH, NCPUH )
      deallocate( IDH )
      deallocate( ICPUMASK, ICPU )
      deallocate( MASKSL )
      call mpi_finalize(info)

      end subroutine para_driver
        
c======================================================================      

      subroutine para_stop

c This is the routine that cleans up the slaves
c This should be called before stopping

      use newpara
      implicit none
      include 'mpif.h'
      integer*4 info

c     if (MY_GID.eq.0) then
        call mpi_abort(NEW_WORLD,20,info)
c     else
c       call mpi_finalize(info)
c     end if
      end subroutine para_stop

c======================================================================

#ifdef DEBUG
      subroutine mpi_errorh( mpi_comm, iercode )

c Error handling for MPI calls. Useful for debugging

      use newpara
#ifdef __INTEL_COMPILER
      use ifcore
#endif

      integer*4 mpi_comm, iercode
      integer*4 ierr,iercl,istrl
      character*256 erstr

      if(MY_GID.eq.0)then
        write(*,*)'*** MPI error on master'
      else
        write(*,'(a,i0)')' *** MPI error on slave ',MY_GID
      endif
      call mpi_error_class(iercode,iercl,ierr)
      erstr='none'
      call mpi_error_string(iercode,erstr,istrl,ierr)
      write(*,*)'Error class: ',iercl
      write(*,*)'Error string: '//trim(erstr)
      write(*,*)'MPI_IBUF=',MPI_IBUF
      write(*,*)'MPI_IPOS=',MPI_IPOS
      write(*,*)'MPI_IBFP=',MPI_IBFP
#ifdef __INTEL_COMPILER
      call tracebackqq(STRING='Internal MPI error',USER_EXIT_CODE=-1)
#endif
      call mpi_abort(NEW_WORLD,20,ierr)

      end subroutine mpi_errorh
#endif
        
c======================================================================      

      subroutine para_check

c no error tolerance in MPI, empty stub      
      end subroutine para_check

c======================================================================      

      subroutine set_mpi_def_int

c sort out the default integer type

      use kinds
      use newpara
      implicit none
      include 'mpif.h'

      select case ( kind(0) )

        case ( i_1 )

          INTEGER_DEF = MPI_INTEGER1

        case ( i_2 )

          INTEGER_DEF = MPI_INTEGER2

        case ( i_4 )

          INTEGER_DEF = MPI_INTEGER4

        case ( i_8 )

          INTEGER_DEF = MPI_INTEGER8

        case default

          write(6,*)'set_mpi_def_int: unknown integer kind'
          call para_stop

      end select

      return

      end subroutine set_mpi_def_int
        
c======================================================================      

      subroutine para_initsend

c initialize buffers for packing
c do not allocate and use new memory while packing and unpacking
c allocate before or after unpacking (for compatibility with MPI)
c send big chunks of data without packing

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer lcore
      integer*4 ibuflen,info

      call getival('lcore',lcore)
      call getmem(0,MPI_IBFP)
      call retmem(1)
      ibuflen=min(lcore-MPI_IBFP,250000000) ! limit to 2GB
      call mpi_pack_size(ibuflen,MPI_REAL8,NEW_WORLD,MPI_IBUF,info)
      MPI_IPOS=0

      end subroutine para_initsend

c======================================================================      

      subroutine para_pack_int(inte,inum)

c pack integer data (MPI)
c data in inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inte,inum
      integer*4 info

      call mpi_pack(inte,inum,INTEGER_DEF,bl(MPI_IBFP),MPI_IBUF, 
     &              MPI_IPOS,NEW_WORLD,info)

      end subroutine para_pack_int
        
c======================================================================      

      subroutine para_unpack_int(inte,inum)

c unpack integer data (MPI)
c data to inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inte,inum
      integer*4 info

      call mpi_unpack(bl(MPI_IBFP),MPI_IBUF,MPI_IPOS,inte,inum, 
     &               INTEGER_DEF,NEW_WORLD,info)

      end subroutine para_unpack_int
        
c======================================================================      

      subroutine para_pack_int4(inte,inum)

c pack integer data (MPI)
c data in inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*4 inte,info

      call mpi_pack(inte,inum,MPI_INTEGER4,bl(MPI_IBFP),MPI_IBUF, 
     &              MPI_IPOS,NEW_WORLD,info)

      end subroutine para_pack_int4
        
c======================================================================      
      subroutine para_unpack_int4(inte,inum)

c unpack integer data (MPI)
c data to inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*4 inte,info

      call mpi_unpack(bl(MPI_IBFP),MPI_IBUF,MPI_IPOS,inte,inum, 
     &                MPI_INTEGER4,NEW_WORLD,info)

      end subroutine para_unpack_int4
        
c======================================================================      

      subroutine para_pack_int2(inte,inum)

c pack integer data (MPI)
c data in inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*2 inte
      integer*4 info

      call mpi_pack(inte,inum,MPI_INTEGER2,bl(MPI_IBFP),MPI_IBUF, 
     &              MPI_IPOS,NEW_WORLD,info)

      end subroutine para_pack_int2
        
c======================================================================      
      subroutine para_unpack_int2(inte,inum)

c unpack integer data (MPI)
c data to inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*2 inte
      integer*4 info

      call mpi_unpack(bl(MPI_IBFP),MPI_IBUF,MPI_IPOS,inte,inum, 
     &                MPI_INTEGER2,NEW_WORLD,info)

      end subroutine para_unpack_int2
        
c======================================================================      

      subroutine para_pack_real(dat,inum)

c pack real data (MPI)
c data in dat (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      real*8 dat
      integer*4 info

      call mpi_pack(dat,inum,MPI_REAL8,bl(MPI_IBFP),MPI_IBUF, 
     &              MPI_IPOS,NEW_WORLD,info)

      end subroutine para_pack_real
        
c======================================================================      

      subroutine para_unpack_real(dat,inum)

c unpack real data (MPI)
c data to dat (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      real*8 dat
      integer*4 info

      call mpi_unpack(bl(MPI_IBFP),MPI_IBUF,MPI_IPOS,dat,inum, 
     &                MPI_REAL8,NEW_WORLD,info)

      end subroutine para_unpack_real
        
c======================================================================      

      subroutine para_pack_string(stri,inum)

c pack string (MPI), max. length of string is 256
c data in string 
c length in inum

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      character*256 stri
      integer*4 info

      call mpi_pack(stri,inum,MPI_CHARACTER,bl(MPI_IBFP),MPI_IBUF, 
     &              MPI_IPOS,NEW_WORLD,info)

      end subroutine para_pack_string

c======================================================================      

      subroutine para_unpack_string(stri,inum)

c unpack string (MPI), max. length of string is 256
c data to string 
c length in inum

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      character*256 stri
      integer*4 info

      call mpi_unpack(bl(MPI_IBFP),MPi_IBUF,MPI_IPOS,stri,inum, 
     &                MPI_CHARACTER,NEW_WORLD,info)

      end subroutine para_unpack_string
        
c======================================================================      

      subroutine para_pack_byte(inte,inum)

c pack byte data (MPI)
c data in inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*1 inte
      integer*4 info

      call mpi_pack(inte,inum,MPI_BYTE,bl(MPI_IBFP),MPI_IBUF, 
     &              MPI_IPOS,NEW_WORLD,info)

      end subroutine para_pack_byte
        
c======================================================================      

      subroutine para_unpack_byte(inte,inum)

c unpack byte data (MPI)
c data to inte (single variable or first member of array)
c array length in inum (1 for single, no bound checking!)

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*1 inte
      integer*4 info

      call mpi_unpack(bl(MPI_IBFP),MPI_IBUF,MPI_IPOS,inte,inum, 
     &                MPI_BYTE,NEW_WORLD,info)

      end subroutine para_unpack_byte
        
c======================================================================      

      subroutine para_send(inte,host,tag)

c send integer to host with message tag

      use newpara
      implicit none
      include 'mpif.h'

      integer host,inte
      integer*4 tag
      integer*4 info

      call mpi_send(inte,1,INTEGER_DEF,host,tag,NEW_WORLD,info)

      end subroutine para_send

c======================================================================      

      subroutine para_recv(inte,ifrom,tag)

c receive integer from any host with message tag
c source is returned in ifrom

      use newpara
      implicit none
      include 'mpif.h'

      integer inte,ifrom
      integer*4 tag
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info

      call mpi_recv(inte,1,INTEGER_DEF,MPI_ANY_SOURCE,tag, 
     &              NEW_WORLD,stat,info)

      ifrom=stat(MPI_SOURCE)

      end subroutine para_recv
        
c======================================================================      
      subroutine para_wait(tag)

c wait for the message identified by tag

      use newpara
      implicit none
      include 'mpif.h'

      integer*4 tag, yield
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info,iflag

      do
        call mpi_iprobe(MPI_ANY_SOURCE,tag,NEW_WORLD,
     $                  iflag,stat,info)
        if(iflag.ne.0)exit

        info=yield()

      enddo

      end subroutine para_wait
c======================================================================      

      subroutine para_send_real(dat,inum,host,tag)

c send real(array) to host with message tag

      use newpara
      implicit none
      include 'mpif.h'

      integer host,inum
      integer*4 tag
      integer*4 info
      real*8 dat

      call mpi_send(dat,inum,MPI_REAL8,host,tag,NEW_WORLD,info)

      end subroutine para_send_real
        
c======================================================================      

      subroutine para_recv_real(dat,inum,host,tag)

c receive real(array) from  host with message tag
c !!not from any host like the other recv routines!!)

      use newpara
      implicit none
      include 'mpif.h'

      integer inum,host
      integer*4 tag
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info
      real*8 dat

      call mpi_recv(dat,inum,MPI_REAL8,host,tag,NEW_WORLD, 
     &              stat,info)

      end subroutine para_recv_real

c======================================================================      

      subroutine para_send_pack(host,tag)

c send packed data to host with message tag

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer host
      integer*4 tag
      integer*4 info

      call mpi_send(bl(MPI_IBFP),MPI_IPOS,MPI_PACKED,host,tag, 
     &              NEW_WORLD,info)

      end subroutine para_send_pack

c======================================================================      

      subroutine para_recv_pack(ifrom,tag)

c receive packed data with message tag from _ANYBODY_
c data is put to unallocated bl memory
c do not allocate and use new memory while packing /unpacking     
c source id returned in ifrom      

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer ifrom
      integer*4 tag
      integer lcore
      integer*4 ibuflen
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info

      call getmem(0,MPI_IBFP)
      call retmem(1)
      call getival('lcore',lcore)
      ibuflen=min(lcore-MPI_IBFP,250000000) ! limit to 2GB
      call mpi_pack_size(ibuflen,MPI_REAL8,NEW_WORLD,MPI_IBUF,info)
      call mpi_recv(bl(MPI_IBFP),MPI_IBUF,MPI_PACKED,MPI_ANY_SOURCE, 
     &              tag,NEW_WORLD,stat,info)
      ifrom=stat(MPI_SOURCE)
      MPI_IPOS=0

      end subroutine para_recv_pack

c======================================================================      

      subroutine para_bcast(inte,tag)

c broadcast int to all slaves with message tag
c tag not used in mpi_bcast      

      use newpara
      implicit none
      include 'mpif.h'

      integer inte
      integer*4 tag
      integer*4 info
      integer idum,ifrom,i

#ifdef SPIN
      call para_barrier
#endif

      call mpi_bcast(inte,1,INTEGER_DEF,0,NEW_WORLD,info)

      end subroutine para_bcast

c======================================================================      

      subroutine para_recv_bcast(inte,tag)

c receive broadcast integer with message tag from master host
c tag not used in mpi_bcast      

      use newpara
      implicit none
      include 'mpif.h'

      integer inte
      integer*4 tag
      integer*4 info

#ifdef SPIN
      call para_barrier
#endif

      call mpi_bcast(inte,1,INTEGER_DEF,0,NEW_WORLD,info)

      end subroutine para_recv_bcast
        
c======================================================================      

      subroutine para_bcast_real(dat,inum,tag)

c broadcast real (array) to all slaves with message tag
c      
c this is just a short wrapper to chop up really big data blocks      
c if the array is bigger than a certain size
c the limit is 4,000,000 words, presently most linux kernels use
c this as maximum shared memory (cat /proc/sys/kernel/shmmax)
c MPI fails if larger chunks are sent on dual proc. machines
c PVM would benefit as well from the smaller buffers      
c internally the PVM/MPI routines are used

      use newpara
      implicit none

      integer inum,i,nchunk,ilen,idum,ifrom
      integer*4 tag
      real*8 dat(inum)

#ifdef SPIN
      call para_barrier
#endif

      nchunk=4000000
      do i=1,inum,nchunk
         ilen=min(inum-i+1,nchunk)
         call para_bcast_float(dat(i),ilen,tag)
      enddo
        
      end subroutine para_bcast_real
        
c======================================================================      

      subroutine para_recv_bcastreal(dat,inum,tag)

c receive broadcast real (array) with message tag from master host
c this is a wrapper to chop up big data blocks see above      

      use newpara
      implicit none

      integer inum,i,nchunk,ilen
      integer*4 tag
      real*8 dat(inum)

#ifdef SPIN
      call para_barrier
#endif

      nchunk=4000000
      do i=1,inum,nchunk
         ilen=min(inum-i+1,nchunk)
         call para_recv_bcastfloat(dat(i),ilen,tag)
      enddo
        
      end subroutine para_recv_bcastreal

c======================================================================      

      subroutine para_bcast_float(dat,inum,tag)

c broadcast real (array) to all slaves with message tag
c this is better than packing in and out for large chunks      
c tag not used in mpi_bcast      
c wrapped by para_bcast_real for big data blocks      

      use newpara
      implicit none
      include 'mpif.h'

      integer inum
      integer*4 tag
      real*8 dat,tred,tbefore,tafter
      integer*4 info

      call getrval('bcast',tred)
      call elapsec(tbefore)
      call mpi_bcast(dat,inum,MPI_REAL8,0,NEW_WORLD,info)
      call elapsec(tafter)
      call setrval('bcast',tred+tafter-tbefore)
        
      end subroutine para_bcast_float
        
c======================================================================      

      subroutine para_recv_bcastfloat(dat,inum,tag)

c receive broadcast real (array) with message tag from master host
c wrapped by para_recv_bcastreal for big data blocks      

      use newpara
      implicit none
      include 'mpif.h'

      integer inum
      integer*4 tag
      real*8 dat
      integer*4 info

      call mpi_bcast(dat,inum,MPI_REAL8,0,NEW_WORLD,info)

      end subroutine para_recv_bcastfloat

c======================================================================      

      subroutine para_bcast_pack(tag)

c broadcast packed data to all slaves with 

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer*4 tag
      real*8 tred,tbefore,tafter
      integer*4 info
      integer idum,ifrom,i

#ifdef SPIN
      call para_barrier
#endif

      call getrval('bcast',tred)
      call elapsec(tbefore)
      call mpi_bcast(MPI_IPOS,1,MPI_INTEGER4,0,NEW_WORLD,info)
      call mpi_bcast(bl(MPI_IBFP),MPI_IPOS,MPI_PACKED,0, 
     &               NEW_WORLD,info)
      call elapsec(tafter)
      call setrval('bcast',tred+tafter-tbefore)
        
      end subroutine para_bcast_pack
        
c====================================================================== 

      subroutine para_barrier

      use newpara
      implicit none
      include 'mpif.h'
      integer*4 ierr

      call MPI_Barrier(NEW_WORLD,ierr)

      end subroutine para_barrier

c====================================================================== 

      subroutine para_recv_bcastpack(tag)

c receive broadcast packed data with message tag from master host
c tag is not needed for broadcast in MPI
c data is put to unallocated bl memory
c do not allocate and use new memory while packing /unpacking     

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer*4 tag
      integer*4 info
      integer nmem

#ifdef SPIN
      call para_barrier
#endif

      call mpi_bcast(MPI_IBUF,1,MPI_INTEGER4,0,NEW_WORLD,info)
      nmem=MPI_IBUF/8+10     ! supposing 8 byte REALs in buffer      
      call getmem(nmem,MPI_IBFP)
      call retmem(1)

      call mpi_bcast(bl(MPI_IBFP),MPI_IBUF,MPI_PACKED,0, 
     &               NEW_WORLD,info)
      MPI_IPOS=0

      end subroutine para_recv_bcastpack

c======================================================================      

      subroutine para_reduce(dat,inum,tag)

c reduce data block
c this is a wrapper to chop up big data blocks

      use newpara
      implicit none

      integer inum,i,nchunk,ilen,idum,ifrom
      integer*4 tag
      real*8 dat(inum)

#ifdef SPIN
      call para_barrier
#endif

c  now do the actual reduce

      nchunk=4000000
      do i=1,inum,nchunk
         ilen=min(inum-i+1,nchunk)
         call para_reduce_block(dat(i),ilen,tag)
      enddo
        
      end subroutine para_reduce

c======================================================================      

      subroutine para_reduce_block(dat,inum,tag)

c sum up real (array) dat from all slaves on master
c wrapped by para_reduce for big data blocks      

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*4 tag
      real*8 dat,tred,tbefore,tafter
      integer*4 info

      call getmem(inum,MPi_IBFP)

c time on master

      if (MY_GID.eq.0) then
         call getrval('reduce',tred)
         call elapsec(tbefore)
         call zeroit(bl(MPi_IBFP),inum)
         call mpi_reduce(dat,bl(MPI_IBFP),inum,MPI_REAL8,MPI_SUM,0, 
     &                   NEW_WORLD,info)
         call dcopy(inum,bl(MPi_IBFP),1,dat,1)
         call elapsec(tafter)
         call setrval('reduce',tred+tafter-tbefore)
      else
         call mpi_reduce(dat,bl(MPI_IBFP),inum,MPI_REAL8,MPI_SUM,0, 
     &                   NEW_WORLD,info)
      endif
      call retmem(1)
        
      end subroutine para_reduce_block
        
c======================================================================      

      subroutine para_probe(tag,isthere)

c isthere returns true if a message correponding to tag has 
c been received

      use newpara
      implicit none
      include 'mpif.h'

      integer*4 tag
      logical isthere
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info,iflag

      call mpi_iprobe(MPI_ANY_SOURCE,tag,NEW_WORLD,iflag,stat,info)
      isthere=iflag.ne.0

      end subroutine para_probe
        
c======================================================================      

      subroutine para_mem(lcore,incore,indisk)

c
c Send memory info to slaves (bl is not set up yet, need buffer)
c
      use newpara
      implicit none
      include 'mpif.h'

      integer lcore,incore,indisk
      character*256 pqs_root,scrfile,scrdir,jobname
      character buffer(2048)
      integer*4 ipos,info

      call getchval('pqs_root',pqs_root)
      call getchval('scrf',scrfile)
      call getchval('scrdir',scrdir)
      call getchval('jobname',jobname)
      ipos=0
      call mpi_pack(lcore,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_pack(incore,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_pack(indisk,1,INTEGER_DEF,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_pack(pqs_root,256,MPI_CHARACTER,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_pack(scrfile,256,MPI_CHARACTER,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_pack(scrdir,256,MPI_CHARACTER,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_pack(jobname,256,MPI_CHARACTER,buffer,2048,ipos, 
     &              NEW_WORLD,info)
      call mpi_bcast(ipos,1,MPI_INTEGER4,0,NEW_WORLD,info)
      call mpi_bcast(buffer,ipos,MPI_PACKED,0,NEW_WORLD,info)

      end subroutine para_mem
        
c======================================================================      

      subroutine para_slave_mem(lcore,incore,indisk,pqs_root,scrfile,
     $                          scrdir,jobname)

c
c Get memory info (bl is not set up yet, need buffer)
c
      use newpara
      implicit none
      include 'mpif.h'

      integer lcore,incore,indisk
      character*256 pqs_root,scrfile,scrdir,jobname
      character buffer(2048)
      integer*4 ipos,info

      call mpi_bcast(ipos,1,MPI_INTEGER4,0,NEW_WORLD,info)
      call mpi_bcast(buffer,ipos,MPI_PACKED,0,NEW_WORLD,info)
      ipos=0
      call mpi_unpack(buffer,2048,ipos,lcore,1,INTEGER_DEF, 
     &                NEW_WORLD,info)
      call mpi_unpack(buffer,2048,ipos,incore,1,INTEGER_DEF, 
     &                NEW_WORLD,info)
      call mpi_unpack(buffer,2048,ipos,indisk,1,INTEGER_DEF, 
     &                NEW_WORLD,info)
      call mpi_unpack(buffer,2048,ipos,pqs_root,256,MPI_CHARACTER, 
     &                NEW_WORLD,info)
      call mpi_unpack(buffer,2048,ipos,scrfile,256,MPI_CHARACTER, 
     &                NEW_WORLD,info)
      call mpi_unpack(buffer,2048,ipos,scrdir,256,MPI_CHARACTER, 
     &                NEW_WORLD,info)
      call mpi_unpack(buffer,2048,ipos,jobname,256,MPI_CHARACTER, 
     &                NEW_WORLD,info)

      end subroutine para_slave_mem


c These needed to link with already compiled MPICH-1
      
        
      integer function mpir_iargc()
      mpir_iargc = iargc()
      return
      end
      subroutine mpir_getarg( i, s )
      integer       i
      character*(*) s
      call getarg(i,s)
      return
      end      
c======================================================================

      subroutine para_printhelp(ichan)
c
c  prints an help message (parallel MPI version)
c
      implicit none
      integer ichan
      write(ichan,100)
 100  format(/' USAGE: pqs_mpi.x [-v,-h] [-i] molecule-name'/)
      end

c======================================================================

      SUBROUTINE para_sendbin(
     $              ibins,       ibin1,   lbin,  icount,
     $              isize,       igran,   idest, sendb,
     $              igranulesize,ngran,   imod,  ndisk2,
     $              nbin,        nrec,    islv)

      use newpara
      implicit none
      include 'mpif.h'
C
C  Send the bin to another slave, MPI version
c
c  A separate MPI version is needed because we need to use
C  a non-blocking send routine, otherwise the MPI version 
c  of the bin sort might freeze if two slaves are trying
c  to send bins to each other at the same time.
C
C  Since we want a non-blocking send, and since MPI does not
c  use internal buffers (any buffer space has to be provided
C  by the application making the calls), we need a buffer space
c  for the send and need to make sure that the contents of the buffer
c  are not modified until the send is complete
C
C  ARGUMENTS
C
C  ibins   -  integer storage for the bins
C  ibin1   -  precision overflow integer
c  lbin    -  kength of one bin
C  icount  -  number of non-zero elements in IJ bin
C  isize   -  size of present granule (in "ij units")
C  igran   -  current granule number
c  idest   -  destination slave
c  sendb   -  buffer space for non-blocking send
C
      integer lbin,isize,igran,idest
      INTEGER*4 ibins(2*lbin)
      INTEGER*1 ibin1(lbin)
      integer*4 icount(isize)
      real*8 sendb
      integer igranulesize,ngran,imod,ndisk2,nbin,nrec,islv
      integer*4 ireq, stat(MPI_STATUS_SIZE), info, isbuf, isbp
      save ireq
      data ireq /-1/
      integer*4 lenb
C      integer*4 flag
      logical*4 flag ! In Fortran, this is a logical parameter. -Zhigang
C It is weird that, if I use logical, it doesn't work. When compiling,
C we use -i8 flag, in this case, the logical parameter and .true. are
C both 8 bytes. But in the program, they are not equal. When I use
C logical*4, flag equals .true.. We should test it using other mpi
C implementation and other compiler. -Zhigang 3/8/2017
c
c  check if previous send is finished, and wait if it is not
c
      if (ireq.ne.-1) then
         do
c        call mpi_wait( ireq, stat, info )
            call mpi_test( ireq, flag, stat, info )
c            write(*,*) "size .true.",bit_size(flag),bit_size(.true.)
            if (flag) then  !(flag.eq.1)=>(flag.eq..true.) or 
               ireq = -1    !In Fortran, .true. and 1 are different
               exit         !Zhigang,3/7/2017 @UARK 
            endif
#ifdef QLOGIC_MPI
         ! no busy_wait needed for Qlogic MPI
#else
            call busy_wait
#endif
            call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                         ndisk2 ,nbin, nrec, islv, .false.)
         enddo
      endif
c
c now pack bin and send it
c
      if (isize.gt.0)then
         lenb=16+4*isize+4*2*lbin*isize+lbin*isize
      else
         lenb=16
      endif

      call mpi_pack_size(lenb,MPI_BYTE,NEW_WORLD,isbuf,info)
      isbp=0
      call mpi_pack(igran,1,INTEGER_DEF,sendb,isbuf,isbp,
     $                                  NEW_WORLD,info) ! granule number
      call mpi_pack(isize,1,INTEGER_DEF,sendb,isbuf,isbp,
     $                                  NEW_WORLD,info) ! granule size
      if (isize.gt.0)then
         call mpi_pack(icount,isize,MPI_INTEGER4,sendb,isbuf,isbp,
     $                                  NEW_WORLD,info) ! no. of elements in bin
         call mpi_pack(ibins,2*lbin*isize,MPI_INTEGER4,sendb,isbuf,isbp,
     $                                  NEW_WORLD,info) ! bin itself
         call mpi_pack(ibin1,lbin*isize,MPI_BYTE,sendb,isbuf,isbp,
     $                                  NEW_WORLD,info) ! precision overfow
      endif
      call mpi_isend(sendb,isbp,MPI_PACKED,idest,TxBinS1,
     $                                  NEW_WORLD,ireq,info) ! non-blocking send

      return
      END
!======================================================================      

      subroutine para_bcast_float_flat(dat,inum,tag)

! broadcast real (array) to all slaves with message tag
! this is better than packing in and out for large chunks      
! probably using pvm_psend and binary tree would be even better
      use newpara
      implicit none
      include 'mpif.h'
      
      integer*4 tag
      integer inum
      integer*4 mpn,ier
      real*8 dat
      real*8 tred,tbefore,tafter
      
      call getrval('bcast',tred)
      call elapsec(tbefore)
      call MPI_Bcast(dat,inum,MPI_REAL8,0,NEW_WORLD,ier)
      if (ier.ne.MPI_SUCCESS) then
        print *, 'Error in MPI_Bcast'
        call f_lush(6)
      endif
      call elapsec(tafter)
      call setrval('bcast',tred+tafter-tbefore)
    
      end

!======================================================================      

      subroutine para_recv_bcastfloat_flat(dat,inum,tag)

! receive broadcast real (array) with message tag from master host
      use newpara
      implicit none
      include 'mpif.h'
      
      integer*4 ier,mpn
      integer*4 tag
      integer inum
      real*8 dat
      mpn=inum
      call MPI_Bcast(dat,inum,MPI_REAL8,0,NEW_WORLD,ier)
      if (ier.ne.MPI_SUCCESS) then
        print *, 'Error in MPI_Bcast'
        call f_lush(6)
      endif
      end subroutine para_recv_bcastfloat_flat

!======================================================================      

      subroutine para_send_int(dat,inum,host,tag)

! send int(array) to host with message tag

      use newpara
      implicit none
      include 'mpif.h'
      
      integer*4 info
      integer*4 tag
      integer host,inum
      integer dat
      
      call MPI_Send(dat,inum,INTEGER_DEF,host,tag,NEW_WORLD,info)
    
      end subroutine para_send_int
    
!======================================================================      

      subroutine para_recv_int(dat,inum,host,tag)

! receive integer(array) from  host with message tag
! !!not from any host like the other recv routines!!)

      use newpara
      implicit none
      include 'mpif.h'
      
      integer*4 info
      integer*4 tag,atid,atag,anum,stat(MPI_STATUS_SIZE)
      integer inum,host
      integer dat
      
      call para_wait(tag)

      call MPI_Recv(dat,inum,INTEGER_DEF,host,tag,NEW_WORLD,stat,
     *              info)
    
      end subroutine para_recv_int

!======================================================================      
      subroutine busy_wait
      implicit none
      integer i,j,generator,ndiv,k,itest
      integer loop_len
      save loop_len
      data loop_len/0/
      real*8 t0,t1
      if (loop_len.eq.0) then
c       print *,'Callibrating busy wait subrotine, wait interval 0.02 s'
        itest=2097152
        do k=1,10
          itest=2*itest
          call secund(t0)
          do i=1,itest
            j=4*generator()
          enddo
          call secund(t1)
          ndiv=nint((t1-t0)/0.02d0)
          if (ndiv.le.3) then
            continue
          else
            loop_len=itest/ndiv
            exit
          endif
        enddo
        if (loop_len.eq.0) loop_len=2147483647
c       print *, 'Busy wait callibrated to ',loop_len,'.'
      endif
c generator must be defined in a separate file, otherwise compiler may optimize
c this loop so that even 2147483647 repetitions will take no time. Linker
c cannot optimize generator calls, but Intel compiler can do it if
c the generator is defined in the same file where it is used. It probably
c simply detects this loop is doing nothing and perfoming the last
c repetition for i==loop_len gives identical result.
      do i=1,loop_len
        j=4*generator()
      enddo
      end
!======================================================================      

      subroutine para_directroute
      ! empty stub in MPI version
        implicit none
      end subroutine para_directroute
      
!======================================================================
      
      subroutine para_nodirectroute
      ! empty stub in MPI version
        implicit none
      end subroutine para_nodirectroute
c=================================================================CC CC

      subroutine para_reduce_max(dat,inum,tag)

c gather max values of real (array) dat from all slaves on master

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer inum
      integer*4 tag
      real*8 dat,tred,tbefore,tafter
      integer*4 info

      call block_alarm
      call getmem(inum,MPi_IBFP)

c time on master

      if (MY_GID.eq.0) then
         call getrval('reduce',tred)
         call elapsec(tbefore)
         call zeroit(bl(MPi_IBFP),inum)
         call mpi_reduce(dat,bl(MPI_IBFP),inum,MPI_REAL8,MPI_MAX,0, 
     &                   NEW_WORLD,info)
         call dcopy(inum,bl(MPi_IBFP),1,dat,1)
         call elapsec(tafter)
         call setrval('reduce',tred+tafter-tbefore)
      else
         call mpi_reduce(dat,bl(MPI_IBFP),inum,MPI_REAL8,MPI_MAX,0, 
     &                   NEW_WORLD,info)
      endif
      call retmem(1)
      call unblock_alarm
        
      end subroutine para_reduce_max
        
c======================================================================      

      subroutine int_probe(tag,slave)

c isthere returns true if a message correponding to tag has 
c been received

      use newpara
      implicit none
      include 'mpif.h'

      integer*4 tag
      integer slave
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info,iflag
      call block_alarm

      call mpi_probe(MPI_ANY_SOURCE,MPI_ANY_TAG,NEW_WORLD,stat,info)
      tag=stat(MPI_TAG)
      slave=stat(MPI_SOURCE)
      call unblock_alarm

      end subroutine int_probe
        
c======================================================================      

      subroutine int_recv_pack(ifrom,tag)

c receive packed data with message tag from ifrom
c data is put to unallocated bl memory
c do not allocate and use new memory while packing /unpacking     
c source id returned in ifrom      

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer ifrom
      integer*4 tag
      integer lcore
      integer*4 ibuflen
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info,ifrom1

      call block_alarm
      call getmem(0,MPI_IBFP)
      call retmem(1)
      call getival('lcore',lcore)
      ibuflen=lcore-MPi_IBFP
      ifrom1=ifrom
      call mpi_pack_size(ibuflen,MPI_REAL8,NEW_WORLD,MPI_IBUF,info)
      call mpi_recv(bl(MPI_IBFP),MPI_IBUF,MPI_PACKED,ifrom1, 
     &              tag,NEW_WORLD,stat,info)
      MPI_IPOS=0
      call unblock_alarm

      end subroutine int_recv_pack

c======================================================================      

      subroutine int_recv(ifrom,tag)

c receive dummy message from host ifrom with the message tag

      use newpara
      implicit none
      include 'mpif.h'

      integer ifrom
      real*8 inte
      integer*4 tag
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info,ifrom1
      ifrom1=ifrom

      call block_alarm
      call mpi_recv(inte,0,MPI_DOUBLE_PRECISION,ifrom1,tag,NEW_WORLD,
     1              stat,info)
      call unblock_alarm

      end subroutine int_recv
        
c======================================================================      

      subroutine para_recv_anypack(ifrom,tag)

c receive packed data with any message tag from _ANYBODY_
c data is put to unallocated bl memory
c do not allocate and use new memory while packing /unpacking     
c source id returned in ifrom  
c tag returned in tag    

      use newpara
      use memory
      implicit none
      include 'mpif.h'

      integer ifrom
      integer*4 tag
      integer lcore
      integer*4 ibuflen
      integer*4 stat(MPI_STATUS_SIZE)
      integer*4 info

      call getmem(0,MPI_IBFP)
      call retmem(1)
      call getival('lcore',lcore)
      ibuflen=min(lcore-MPI_IBFP,250000000) ! limit to 2GB
      call mpi_pack_size(ibuflen,MPI_REAL8,NEW_WORLD,MPI_IBUF,info)
      call mpi_recv(bl(MPI_IBFP),MPI_IBUF,MPI_PACKED,MPI_ANY_SOURCE, 
     &              MPI_ANY_TAG,NEW_WORLD,stat,info)
      ifrom=stat(MPI_SOURCE)
      tag=stat(MPI_TAG)
      MPI_IPOS=0

      end subroutine para_recv_anypack

c======================================================================      
c I modified this subroutine according to the source code from PQSv40 
c Now MPI version Coupled-Cluster calculation is OK!
c Zhigang 3/16/2017 @UARK
c======================================================================
      subroutine fafinit(info)  ! collective!
! splits into I/O daemons and master/slave group
      use newpara
      implicit none
      include 'mpif.h'
      integer, intent(out) :: info !are we I/O?
      integer color, firstio
      integer*4 ii,my_rank,io_size,InterComm
!      call set_color(color)
      color=1    ! normal processes have color 1
! the second process on each host will be an I/O daemon, color 2
      if (MY_GID .eq. IDH(MY_HOSTID,2)) then
          color=2
      endif
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ii)
      call MPI_Comm_split(MPI_COMM_WORLD,color,my_rank,NEW_WORLD,ii)
      call MPI_Comm_rank(NEW_WORLD,MY_GID,ii)
      call MPI_Comm_size(NEW_WORLD,NSLV,ii)
      NSLV=NSLV-1
      MY_ID=MY_GID
      call create_c_communicators(color)
      if (color.eq.2) then  ! I/O daemon
C        call MPI_Intercomm_create(NEW_WORLD,0,MPI_COMM_WORLD,0,
C     &                            TxSpinDat,InterComm,ii)
C        call MPI_Comm_size(InterComm,io_size,ii)
C        call MPI_Comm_rank(InterComm,my_rank,ii)
C        call afdx(my_rank,NEW_WORLD,InterComm,io_size,ii)
        call afdx
C        call MPI_Comm_free(InterComm,ii)
        call MPI_Comm_free(NEW_WORLD,ii)
        call MPI_Comm_rank(MPI_COMM_WORLD,MY_GID,ii)
        call MPI_Comm_size(MPI_COMM_WORLD,NSLV,ii)
        NEW_WORLD=MPI_COMM_WORLD
        NSLV=NSLV-1
        MY_ID=MY_GID
        info=1
        call remove_c_communicators(color)
      else  ! on normal processes
C        firstio=IDH(1,2)
C        if (firstio.ne.1) then
C!  find the lowest ranked I/O daemon, usually rank 1
C            do ii=2,NHOST
C               if(IDH(ii,2).lt.firstio) firstio=IDH(ii,2)
C            enddo
C        endif
C        call MPI_Intercomm_create(NEW_WORLD,0,MPI_COMM_WORLD,firstio,
C     &                            TxSpinDat,InterComm,ii)
C        call MPI_Comm_remote_size(InterComm,io_size,ii)
C        call afinit(my_rank,NEW_WORLD,Intercomm,io_size)
        info=0
      endif
      end

c======================================================================      

      subroutine fafterminate(info)
      use newpara
      implicit none
      include 'mpif.h'
      integer info
      integer*4 ii
      call afterminate(info)
      call MPI_Comm_rank(MPI_COMM_WORLD,MY_GID,ii)
      call MPI_Comm_size(MPI_COMM_WORLD,NSLV,ii)
      call MPI_Comm_free(NEW_WORLD,ii)
      NEW_WORLD=MPI_COMM_WORLD
      NSLV=NSLV-1
      MY_ID=MY_GID
      end

c======================================================================
      subroutine fafreopen
      end
c======================================================================
      subroutine sync_barrier
      end
c======================================================================
      subroutine array_files
      end
c======================================================================
      subroutine fafnbread
      end
c======================================================================
      subroutine fafrequest
      end
c======================================================================
      subroutine block_alarm
      end
c======================================================================
      subroutine unblock_alarm
      end
c======================================================================
