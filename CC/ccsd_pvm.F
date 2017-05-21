#include "maxi.h"
c==============================================================================
      subroutine 
     *   CoefInit_master(npairs, ncf,gauss_seidel,nodisk,   ccsd_inp,
     *                   vorbpar,nmo,ivalpair,    ivirpair, ifp,
     *                   nsym,   nval,nvirt,      ij_unique,isympairs,
     *                   restart,ndisk_dump,      nbf)
c This subroutine performs the same tasks as single processor version
c CoefInit. The additional task is that it sends to slaves the files unit
c numbers used internally. Also it is not supposed to work in 'nodisk' mode.
      use memory
      use newpara
      use kinds
      use ccounters
      implicit real*8 (a-h,o-z)
      character*256 filenamea,filenameb,scrfile,iternumbera,iternumberb
      character*256 filenamer
      character*2  ch_mygid
      logical gauss_seidel,nodisk,nodiskk,ccsd,ccsd_inp,fileopen
      logical vorb,vorbpar,af,restart
      common /ampldisk/ndiska,ndiskb,ncfa,internal_addr,nodiskk,ireturn,
     *                 icallno_propagate,npairs_propagate,ccsd,ndiskr,
     *                 ivirt,vorb,af
      common /amplsymm/ ivalpair1,ivirpair1,ifp1,nsym1,nval1,nvirt1,
     *                  ij_unique1,isympairs1
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      save /ampldisk/
      af=.true.
      ivirt=nbf-nmo
      vorb=vorbpar
      ivalpair1 = ivalpair
      ivirpair1 = ivirpair
      ifp1      = ifp
      nsym1     = nsym
      nval1     = nval
      nvirt1    = nvirt
      ij_unique1= ij_unique
      isympairs1= isympairs
      if (vorb) then
        idimen=ivirt
      else
        idimen=ncf
      endif
      ic_coefmast=ic_coefmast+1
      icallno_propagate=ic_coefmast
      npairs_propagate=npairs
      ccsd=ccsd_inp
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      if (ic_coefmast.eq.1) call dynamic_getmem(idimen*idimen,ireturn)
      if (nodisk) then
        STOP 'No support for nodisk tmp, should be added?'
        if (ic_coefmast.eq.1) then 
          nodiskk=nodisk
          ncfa=ncf
          call dynamic_getmem(npairs*idimen*idimen,ii)
          internal_addr=ii
          do i=1,npairs*idimen*idimen
          bl(internal_addr-1+i)=0.0d0
          enddo
        endif
        return
      endif
      if (ic_coefmast.eq.1)call dynamic_getmem(idimen*idimen,izerocoeff)
c
      internal_addr=0
      if (gauss_seidel) then
      if (ic_coefmast.eq.1) then
      ncfa=ncf
      nodiskk=nodisk
      write(iternumbera,'(I3.3)') ic_coefmast
      filenamea=scrfile(1:len)//'.ampl.'//iternumbera(1:3)
      len1=len+6+3
      call fafCreatem(filenamea(1:len1)//char(0),ndiska,idimen*idimen*8)
      ndiskb=ndiska
      else
      call fafCreatem(filenamea(1:len1)//char(0),ndiska,-1) ! dummy
      endif
      else
c=================
      write(iternumbera,'(I3.3)') ic_coefmast
      write(iternumberb,'(I3.3)') ic_coefmast+1
      filenamea=scrfile(1:len)//'.ampl.'//iternumbera(1:3)
      filenameb=scrfile(1:len)//'.ampl.'//iternumberb(1:3)
      filenamer=scrfile(1:len)//'.differen'
      len1=len+6+3
      if (ic_coefmast.eq.1) then
        ncfa=ncf
        ndiskr=62
        nodiskk=nodisk
      call fafCreatem(filenamea(1:len1)//char(0),ndiska,idimen*idimen*8)
      else
       call fafCreatem(filenamea(1:len1)//char(0),ndiska,-1) ! dummy
      endif
      if (ic_coefmast.gt.1) then
c Temporary switched off because no differences are needed now.
        call fafClosem(ndiska,0,info)
        ndiska=ndiskb
      else
        call fafClosem(ndiska,-1,info) ! dummy
      endif
      call fafCreatem(filenameb(1:len1)//char(0),ndiskb,idimen*idimen*8)
c=================
      endif
      call sync_barrier
      if (ic_coefmast.eq.1) then
        do j=1,idimen
          do i=1,idimen
            bl(izerocoeff+(i-1)+(j-1)*idimen)=0.0d0
          enddo
        enddo
        if (restart) then
          do ij=1,npairs
          call general_read1(ndisk_dump,bl(izerocoeff),idimen*idimen*8)
          call fafWrite(ndiska,bl(izerocoeff),8,idimen*idimen,1,ij,info)
          enddo
        else
          do ij=1,npairs
          call fafWrite(ndiska,bl(izerocoeff),8,idimen*idimen,1,ij,info)
          enddo
        endif
      endif
      if (ic_coefmast.eq.1)  call dynamic_retmem(1)
      call sync_barrier
      call para_initsend
      call para_pack_int(ndiska,1)
      call para_pack_int(ndiskb,1)
      call para_bcast_pack(TxCoefInit)
      end
c==============================================================================
      subroutine CCalphamast(ncf,   nval,   irecadrx, npairs,    ndiskx,
     *                       lbinx, thresh, doit,     byt8,      ccsd,
     *                       nmo,   vorb,   elaps,    ndiskalpha,af,
     *                       nbf)
c This is parallelized version of CCalphaonce subroutine. The parallelization
c is over external loop, driven by master. Master simply assignes work to
c slaves by sending them indices from external loop. Slaves have their external
c loop built as infinite loop and they work as long as master sends them
c iwork variable = 1. This is standard parallelization used in most places.
c For description of parameters see CCalphaonce.
      use memory
      use newpara
      use kinds
      implicit none
      integer ncf,nval,irecadrx(*),npairs,ndiskx,lbinx,nmo,ndiskalpha
      REAL*8 thresh,alpha,elaps(*)
      integer iexchAO,icoeffAO,mataddr,i,j,k,l,iaddress,icounter,nbf
      integer inval,icurraddr,nval2,nval3,initial_chunk
      logical doit,idoit,byt8,ccsd,vorb,af
      parameter (initial_chunk=25)
      REAL*8 elem,elaps0,elaps1
      common /CCAlphaCom/iaddress,inval,idoit
c new
      integer memavail,lcore,iad,ncf2,number_matrices,ichunk_T,ichunk_K
      integer iKstore,iTstore,iresult,iKpass,iTpass,iT,iK,iTstart,iTstop
      integer iTsize,iKstart,iKstop,iKsize,ij,kl,number_matrices1,idimen
      integer idimen2,ialphabuffer,imax_per_slave,iwork,islvid,info
      integer islgid
      logical directly
c
      save /CCAlphaCom/
c
c   For nslv divide:
c     if (nslv.gt.1) then
c       nslvdiv=nslv-1
c     else
c       nslvdiv=1
c     endif
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      idimen2=idimen*idimen
c
      inval=nval
      idoit=doit
      ncf2=ncf*ncf
      nval2=nval*nval
      nval3=nval*nval*nval
      call mmark
      call dynamic_mmark
      call matmark
c
c  Calculate space available for K and T blocks.
c
      call dynamic_show_free(memavail)
      memavail=memavail-300000 ! due to possible sparsity in memory allocator
c
      number_matrices=memavail/(idimen2) ! +1 for result matrix
c If all T matrices and initial_chunk K matrices will fit the memory
      looper: do
c round to nearest greater whole number:
        imax_per_slave=AInt(dble(npairs)/dble(nslv)+0.999999999999999d0)
        if ((imax_per_slave+initial_chunk).le.number_matrices) then
         ichunk_T=imax_per_slave
         number_matrices1=number_matrices-ichunk_T
         if (number_matrices1.gt.npairs) then 
           ichunk_K=npairs
         else
           ichunk_K=number_matrices1
         endif
        else
          ichunk_T=number_matrices-initial_chunk
          ichunk_K=initial_chunk
        endif
        if (ichunk_T.lt.initial_chunk.or.ichunk_K.lt.initial_chunk) then
          if (ichunk_T.ne.npairs.and.ichunk_K.ne.npairs) then
            if (ichunk_T.lt.imax_per_slave) then
            print *, 'ichunk_T: ',ichunk_T
            print *, 'ichunk_K: ',ichunk_K
            print *, 'idimen: ',idimen
            print *, 'Too little memory in alpha builder'
            STOP 'Too little memory in alpha builder'
            endif
          endif
        endif
        if ((idimen2*(ichunk_T+ichunk_K)+ichunk_T*ichunk_K+
     *                             ichunk_T*nval*nval).le.memavail)
     *    then
          exit
        else
         number_matrices=number_matrices-1
        endif
      enddo looper
c
      call para_initsend
      call para_pack_int(ichunk_K,1)
      call para_pack_int(ichunk_T,1)
      call para_pack_int(ndiskalpha,1)
      call para_bcast_pack(TxCCSDInit)
c
      iTpass=npairs/ichunk_T
      if (mod(npairs,ichunk_T).ne.0) iTpass=iTpass+1
c Loops over passes:
      iwork=1
      do iT=1,iTpass
        iTstart=ichunk_T*(iT-1)+1
        iTstop =ichunk_T*iT
        if (iTstop.gt.npairs) iTstop=npairs
        iTsize=iTstop-iTstart+1
        call para_recv(islvid,islgid,TxCCSDAReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(iTstart,1)
        call para_pack_int(iTstop,1)
        call para_pack_int(iTsize,1)
        call para_send_pack(islgid,TxCCSDAJob)
      enddo
      iwork=0
      do i=1,nslv
        call para_recv(islvid,islgid,TxCCSDAReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,TxCCSDAJob)
      enddo
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine B41a_G_Master(ncf,nval,norecalc,cc,ccsd,nfirst,F,
     *                        nmo,vorb,ndisk,af,ndiskalpha,betam,nbf)
c This is parallelized version of B41a_G subroutine. The parallelization
c is over external loop, driven by master. Master simply assignes work to
c slaves by sending them indices from external loop. Slaves have their external
c loop built as infinite loop and they work as long as master sends them
c iwork variable = 1. This is standard parallelization used in most places.
c For description of parameters see B41a_G.
      use memory
      use newpara
      use kinds
      implicit none
      integer ncf,nval,nfirst,nmo,ndisk,ndiskalpha,nbf
      logical norecalc,cc,ccsd,vorb,af
C
      integer iad,lcore,memavail,initial_chunk,number_matrices,ncf2
      integer npairs,ichunk_T,number_matrices1,ichunk_R,iTstore,iRstore
      integer ialphas,iT,iTpass,iTstart,iTstop,iTsize,iR,iRpass
      integer iRstart,iRstop,iRsize,itwo,idimen,idimen2,ialphabase
      integer imax_per_slave,info,iwork,islvid,i,islgid
      logical transposed,reset
      real*8 F(ncf,ncf),betam(nval,nval)
      parameter (initial_chunk=25)
C
c     if (.not.norecalc) return
c
c   For nslv divide:
c     if (nslv.gt.1) then
c       nslvdiv=nslv-1
c     else
c       nslvdiv=1
c     endif
c
      reset=.true. ! for storeAlpha
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      ncf2=ncf*ncf
      idimen2=idimen*idimen
      npairs=nval*(nval+1)/2
      call dynamic_show_free(memavail)
      memavail=memavail-300000 ! due to possible sparsity in memory allocator
c We need memory for alpha matrix and C matrix and result matrix.
c How much we can get?
      number_matrices=memavail/(idimen2)
c If all T matrices and initial_chunk Res matrices will fit the memory
      looper: do
c round to nearest greater whole number:
        imax_per_slave=AInt(dble(npairs)/dble(nslv)+0.999999999999999d0)
        if ((imax_per_slave+initial_chunk).le.number_matrices) then
         ichunk_R=imax_per_slave
         number_matrices1=number_matrices-ichunk_R
         if (number_matrices1.gt.npairs) then 
           ichunk_T=npairs
         else
           ichunk_T=number_matrices1
         endif
        else
          ichunk_R=number_matrices-initial_chunk
          ichunk_T=initial_chunk
        endif
        if (ichunk_R.lt.1.or.ichunk_T.lt.initial_chunk) then
          if (ichunk_R.ne.npairs.and.ichunk_T.ne.npairs) then
            if (ichunk_T.lt.imax_per_slave) then
            write(6,*)  'ichunk_R: ',ichunk_R
            write(6,*)  'ichunk_T: ',ichunk_T
            write(6,*)  'idimen: ',idimen
            write(6,*)  'Too little memory in 41a builder'
            STOP 'Too little memory in 41a builder'
            endif
          endif
        endif
        if ((idimen2*(ichunk_R+ichunk_T)+ichunk_R*ichunk_T+
     *        ichunk_R*nval*nval).le.memavail)
     *    then
          exit
        else
         number_matrices=number_matrices-1
c        write(6,*)  'Number of the matrices reduced by 1, looping...'
        endif
      enddo looper
c
c
c  Initialize file for 41a matrices
      call Fileinit41a(idimen,ndisk,af)
c
      call sync_barrier
      call para_initsend
      call para_pack_int(ichunk_R,1)
      call para_pack_int(ichunk_T,1)
      call para_pack_int(ndisk,1)
      call para_bcast_pack(TxCCSDInit)
c calculate passes:
      iRpass=npairs/ichunk_R
      if (mod(npairs,ichunk_R).ne.0) iRpass=iRpass+1
c Loops over passes:
      iwork=1
      do iR=1,iRpass
        iRstart=ichunk_R*(iR-1)+1
        iRstop =ichunk_R*iR
        if (iRstop.gt.npairs) iRstop=npairs
        iRsize=iRstop-iRstart+1
        call para_recv(islvid,islgid,TxCCSDAReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(iRstart,1)
        call para_pack_int(iRstop,1)
        call para_pack_int(iRsize,1)
        call para_send_pack(islgid,TxCCSDAJob)
      enddo
      iwork=0
      do i=1,nslv
        call para_recv(islvid,islgid,TxCCSDAReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,TxCCSDAJob)
      enddo
c Synchronize:
c
      call matremark
      call retmark
      call dynamic_retmark
      call para_barrier
      end
c==============================================================================
      subroutine EEO_vector_master(nval,irecadre,npairs,ndiske,ncf,
     *                             lbine,thresh,nfirst,nlast,byt8,
     *                             nmo,vorb,nslv1,iresult,nbf)
c This is parallelized version of EEO_vector subroutine. The parallelization
c is over external loop, driven by master. Master simply assignes work to
c slaves by sending them indices from external loop. Slaves have their external
c loop built as infinite loop and they work as long as master sends them
c iwork variable = 1. This is standard parallelization used in most places.
c For description of parameters see EEO_vector.
      use memory
      use newpara
      use kinds
      implicit none
      integer nval,npairs,ndiske,nbf
      integer ncf,lbine,nfirst,nlast,iresult,nmo,nslv1
      integer irecadre(npairs)
      real*8 thresh
      logical byt8,vorb
c
      integer i,k,mataddr,idimen,intermid,iwork,info,ibufid
      integer islvtid,islgid
c
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('intermid','r',ncf,nval)
      call matzero('intermid')
      intermid=mataddr('intermid')
      call matconn('t_result','r',idimen,nval,iresult)
c
      iwork=1
      do k=1,nval
        do i=1,k
        call para_recv(islvtid,islgid,TxCCSDReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(k,1)
        call para_pack_int(i,1)
        call para_send_pack(islgid,TxCCSDJob)
        enddo
      enddo
      iwork=0
      do i=1,nslv
        call para_recv(islvtid,islgid,TxCCSDReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,TxCCSDJob)
      enddo
c
      call sync_barrier
      call para_barrier
c
      iwork=1
      do k=1,nval
        call para_recv(islvtid,islgid,TxCCSDReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(k,1)
        call para_send_pack(islgid,TxCCSDJob)
      enddo
      iwork=0
      do i=1,nslv
        call para_recv(islvtid,islgid,TxCCSDReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,TxCCSDJob)
      enddo
      call sync_barrier
      call para_reduce(bl(iresult),idimen*nval,TxCCSDRes)
c     if (vorb) then
c       call matmmul2('genvirt','intermid','t_result','t','n','n')
c     else
c       call matcopy('intermid','t_result')
c     endif
      call matremark
      call retmark
      call dynamic_retmark
      call para_barrier
      end
c==============================================================================
      subroutine Generate_MP2_G_MAS(nval,ncf,ndisk,nmo,nfirst,vorb,af,F,
     *                              nslv1,nbf)
c This is parallelized version of Generate_MP2_G subroutine. The parallelization
c is over external loop, driven by master. Master simply assignes work to
c slaves by sending them indices from external loop. Slaves have their external
c loop built as infinite loop and they work as long as master sends them
c iwork variable = 1. This is standard parallelization used in most places.
c For description of parameters see Generate_MP2_G.
      use memory
      use newpara
      use kinds
      implicit none
      integer nval,ncf,ndisk,nmo,nfirst,nslv1,nbf
      logical vorb,af
      real*8 F(ncf,ncf)
c local:
c
      integer i,j,k,iresult,iamplit,ka,ir,ij,idimen,ifp,jfp,itmp,ifock
      integer istatus,info,iwork,islvtid,ibufid,islgid
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      iwork=1
      call para_barrier
      do j=1,nval   ! result column
        call para_recv(islvtid,islgid,MP2_G_GEN_req)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(j,1)
        call para_send_pack(islgid,MP2_G_GEN_job)
      enddo
      iwork=0
      do j=1,nslv   ! result column
        call para_recv(islvtid,islgid,MP2_G_GEN_req)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,MP2_G_GEN_job)
      enddo
      end
c==============================================================================
      subroutine Kijk_Vec_Mast(irecadrx,npairs,ndiskx,lbinx,ncf,thresh,
     *                     ifirst,last,byt8,nmo,vorb,af,ndiska,nbf)
c This is parallelized version of Kijk_Vec subroutine. The parallelization
c is over external loop, driven by master. Master simply assignes work to
c slaves by sending them indices from external loop. Slaves have their external
c loop built as infinite loop and they work as long as master sends them
c iwork variable = 1. This is standard parallelization used in most places.
c For description of parameters see Kijk_Vec.
      use memory
      use newpara
      use kinds
      implicit none
      integer irecadrx,npairs,ndiskx,lbinx,ncf,ifirst,last,nmo,ndiska
      integer nbf
      real*8 thresh
      logical byt8,vorb,vorba,af,afa
      dimension irecadrx(npairs)
c
      integer ndisk,nval,nval2,irecord,len,len1,mataddr,imatrix,iresult
      integer nrec,i,j,k,ncfa,idimen,iresultMO,info,ndisk_cache,iaddress
      integer nbfa
      common /kijk/ ndisk,ncfa,nval,idimen,vorba,afa,ndisk_cache,nbfa
      logical fileopen
      character*256 scrfile,filename
      character*2 ch_mygid
      integer mygid,ibufid,iwork,isl_tid,itempk,nrec1,islgid
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      save /kijk/
c
      vorba=vorb
      afa=af
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ndisk_cache=76
      if (af) then
        mygid=MY_GID
        write(ch_mygid,'(I2.2)') mygid
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
        filename=scrfile(1:len)//'.kijk--'//ch_mygid
        len1=len+9
        irecord=idimen*8
        OPEN(UNIT=ndisk_cache,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      endif
c
      nval=last-ifirst+1
      nval2=nval*nval
      ncfa=ncf
      nbfa=nbf
      ndisk=76
      irecord=idimen*8
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filename=scrfile(1:len)//'.kijk'
      len1=len+5
c
      if (af) then
         call fafCreatem(filename(1:len1)//char(0),ndisk,nval*idimen*8)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 76 was used by somebody'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      endif
c
      ndiska=ndisk
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('matrix','q',ncf,ncf)
      imatrix=mataddr('matrix')
      call dynamic_matdef('result','r',ncf,nval)
      iresult=mataddr('result')
      call matsub('occu','cano',ifirst,last)
      if (vorb) then
        call dynamic_matdef('resultMO','r',idimen,nval)
        iresultMO=mataddr('resultMO')
      endif
c to slave:
      call para_bcast(ndisk,TxCCSDInit)
      iwork=1
      do i=1,nval
        do j=1,nval
          call para_recv(isl_tid,islgid,TxCCSDReq)
          call para_initsend
          call para_pack_int(iwork,1)
          call para_pack_int(i,1)
          call para_pack_int(j,1)
          call para_send_pack(islgid,TxCCSDJob)
        enddo
      enddo
      iwork=0
ccs      call sleep(20)
      do i=1,nslv
        call para_recv(isl_tid,islgid,TxCCSDReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,TxCCSDJob)
      enddo
c Local cache:
      call para_barrier
      call dynamic_matdef('tempk','r',idimen,nval)
      itempk=mataddr('tempk')
      nrec=0
      nrec1=0
      do i=1,nval
        do j=1,nval
          nrec=nrec+1
          call fafRead(ndisk,bl(itempk),8,idimen*nval,1,nrec,info)
          if (info.lt.8*idimen*nval)
     *    ERROR
          do k=1,nval
            nrec1=nrec1+1
            iaddress=itempk+(k-1)*idimen
            call vectorwrite(ndisk_cache,nrec1,bl(iaddress),idimen)
          enddo
        enddo
      enddo
      call dynamic_matrem('tempk')
      call matremark
      call retmark
      call dynamic_retmark
      call sync_barrier
      end
c==============================================================================
      subroutine KijklInit_Master(irecadrx,npairs,ndiskx,lbinx,ncf,
     *                            thresh,ifirst,last,F,byt8,
     *                            nmo,vorb,ndisk,af,nslv1)
c This is parallelized version of KijklInit subroutine. The parallelization
c is over external loop, driven by master. Master simply assignes work to
c slaves by sending them indices from external loop. Slaves have their external
c loop built as infinite loop and they work as long as master sends them
c iwork variable = 1. This is standard parallelization used in most places.
c For description of parameters see KijklInit.
      use memory
      use newpara
      use kinds
      implicit none
      integer npairs,ndiskx,lbinx,ncf,ifirst,last,nmo,ndisk,nslv1
      integer irecadrx(npairs)
      real*8 thresh,F(ncf,ncf)  ! F-ock in MO basis
      logical byt8,vorb,af
c
      integer nval,imatrix,mataddr,ikij,ij,i,j,istatus,info,ibufid
      integer iwork,isl_tid,islgid
c
      nval=last-ifirst+1
      call para_barrier
      call para_bcast(ndisk,TxCCSDInit)
      ij=0
      iwork=1
      do i=1,nval
        do j=1,i
          ij=ij+1
          call para_recv(isl_tid,islgid,TxCCSDReq)
          call para_initsend
          call para_pack_int(iwork,1)
          call para_pack_int(i,1)
          call para_pack_int(j,1)
          call para_pack_int(ij,1)
          call para_send_pack(islgid,TxCCSDJob)
        enddo
      enddo
      iwork=0
      do i=1,nslv
        call para_recv(isl_tid,islgid,TxCCSDReq)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_send_pack(islgid,TxCCSDJob)
      enddo
      end
c==============================================================================
      subroutine par_diis_matrix(singles,mp2_cc,npairs,
     *                       irecord_tab,max_rec,ndiskr,
     *                       idiis_last,icalculate,nbf)
      use memory
      use newpara
      implicit none
      logical singles
      character*(*) mp2_cc
      integer npairs,irecord_tab(*),max_rec,nbf
      integer ndiskr,idiis_last,icalculate
c
      integer info,isl,iextra,ipass_size,ipass_reminder
      integer istart,istop,mp2_cc_len,isl_tid,islgid
      real*8  diis_row_matrix(icalculate)
c
        call sync_barrier
c
        mp2_cc_len=len(mp2_cc)
        call para_initsend
        call para_pack_int(mp2_cc_len,1)
        call para_pack_string(mp2_cc,mp2_cc_len)
        call para_pack_int(max_rec,1)
        call para_pack_int(irecord_tab,max_rec)
        call para_pack_int(ndiskr,1)
        call para_pack_int(icalculate,1)
        call para_bcast_pack(TxCCSDInit)
c
        if (singles) then 
          iextra=1
        else
          iextra=0
        endif
        ipass_size=(npairs+iextra)/nslv
        ipass_reminder=mod(npairs+iextra,nslv)
        istart=1
        istop=0
        do isl=1,nslv
          istart=istop +1
          istop =istart+ipass_size-1
          if (ipass_reminder.gt.0) istop=istop+1
          ipass_reminder=ipass_reminder-1
          if (ipass_reminder.lt.0.and.ipass_size.eq.0) then
            istart=-1
            istop =-1
          endif
          if (istop.gt.npairs+iextra) 
     *                STOP 'Cannot be true, par_diis_matrix'
          call para_recv(isl_tid,islgid,TxDiisReq)
          call para_initsend
          call para_pack_int(istart,1)
          call para_pack_int(istop ,1)
          call para_send_pack(islgid,TxDiisJob)
        enddo
        istart=0
        istop =0
        do isl=1,nslv
          call para_recv(isl_tid,islgid,TxDiisReq)
          call para_initsend
          call para_pack_int(istart,1)
          call para_pack_int(istop ,1)
          call para_send_pack(islgid,TxDiisJob)
        call para_recv_real(diis_row_matrix,icalculate,islgid,TxDiisRes)
          call add(diis_row_matrix,bl(idiis_last+1),icalculate)
        enddo
      end
C===============================================================================
      subroutine sla_diis_matrix(af,singles,npairs,list,ipairimages,
     *                           ncf,nmo,residVOcurr,residVOprev,
     *                           iMasID,mygid,mytid,nbf)
      use memory
      use newpara
      implicit none
      logical af,singles
      integer npairs,list(*),ipairimages(*),nbf
      integer ncf,nmo,iMasID,mygid,mytid
      real*8 residVOcurr(nbf-nmo,nbf-nmo)
      real*8 residVOprev(nbf-nmo,nbf-nmo)
c
      character*(80) mp2_cc
      integer irecord_tab(DIIS_LENGTH),max_rec,ndiskr
      integer ivirt,ij,irecord,info,i,k,icalculate,mp2_cc_len
      integer istart,istop,imastgid
      real*8 trace,diis_last(DIIS_LENGTH)
c
      ivirt=nbf-nmo
      call sync_barrier
c
      do i=1,80
        mp2_cc(i:i)=' '
      enddo
      do i=1,DIIS_LENGTH
        diis_last(i)=0.0d0
      enddo
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(mp2_cc_len,1)
      call para_unpack_string(mp2_cc,mp2_cc_len)
      call para_unpack_int(max_rec,1)
      call para_unpack_int(irecord_tab,max_rec)
      call para_unpack_int(ndiskr,1)
      call para_unpack_int(icalculate,1)
c
      do
        call para_send(mytid,0,TxDiisReq)
        call para_recv_pack(imastgid,TxDiisJob)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop ,1)
        if (istart.eq.-1) cycle
        if (istart.eq. 0) then
          call para_send_real(diis_last,icalculate,0,TxDiisRes)
          exit
        endif
        if (singles) then
        do ij=istart,istop
          if (ij.eq.(npairs+1).and. mp2_cc.eq.'mp2') exit
          if (ij.le.npairs) then
            if (list(ij).ne.1) cycle
            if (ipairimages(ij).eq.0) cycle
          endif
          irecord=(irecord_tab(max_rec)-1)*(npairs+1)+ij
          if (af) then
            call fafRead(ndiskr,residVOcurr,8,ivirt*ivirt,
     *                    1,irecord,info)
            if (info.lt.8*ivirt*ivirt)
     *    ERROR
          else
            call reader(ndiskr,ivirt,irecord,residVOcurr)
          endif
          do i=1,icalculate ! OPTIMIZE! LAST RECORD = residVOcurr
            k=i+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*(npairs+1)+ij
            if (af) then
              call fafRead(ndiskr,residVOprev,8,ivirt*ivirt,
     *                      1,irecord,info)
            if (info.lt.8*ivirt*ivirt)
     *    ERROR
            else
              call reader(ndiskr,ivirt,irecord,residVOprev)
            endif
              call tjmatprodtr1(residVOcurr,residVOprev,trace,ivirt)
            diis_last(i)=diis_last(i)+trace*ipairimages(ij)
          enddo
        enddo
        else !  if (singles) then
        do ij=istart,istop
          irecord=(irecord_tab(max_rec)-1)*npairs+ij
          if (ij.le.npairs) then
            if (list(ij).ne.1) cycle
            if (ipairimages(ij).eq.0) cycle
          endif
          if (af) then
            call fafRead(ndiskr,residVOcurr,8,ivirt*ivirt,
     *                    1,irecord,info)
            if (info.lt.8*ivirt*ivirt)
     *    ERROR
          else
            call reader(ndiskr,ivirt,irecord,residVOcurr)
          endif
          do i=1,icalculate
            k=i+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*npairs+ij
            if (af) then
              call fafRead(ndiskr,residVOprev,8,ivirt*ivirt,
     *                      1,irecord,info)
            if (info.lt.8*ivirt*ivirt)
     *    ERROR
            else
              call reader(ndiskr,ivirt,irecord,residVOprev)
            endif
            call tjmatprodtr1(residVOcurr,residVOprev,trace,ivirt)
            diis_last(i)=diis_last(i)+trace*ipairimages(ij)
          enddo
        enddo
        endif
      enddo
      end
C===============================================================================
      subroutine par_new_coeff(singles,af,mp2_cc,idimen,nval,list,
     *                          ipairimages,icalculate,max_rec,
     *                          irecord_tab,npairs,ndiskc,diiscoeff)
      use memory
      use newpara
      implicit none
      logical singles,af
      character mp2_cc*(*)
      integer idimen,nval,list(*),ipairimages(*),icalculate,max_rec
      integer irecord_tab(*),npairs,ndiskc
      real*8 diiscoeff(*)
c
      integer icoeff,ifinalcoeff,mataddr,ij,i,j,ii,k,irecord,info
      integer mp2_cc_len,ipass_size,ipass_reminder,istart,istop,isl
      integer isl_tid,islgid
c
c
      call sync_barrier
c
      mp2_cc_len=len(mp2_cc)
      call para_initsend
      call para_pack_int(icalculate,1)
      if (icalculate.le.1) then
        call para_bcast_pack(TxCCSDInit)
        return
      else
        call para_pack_int(mp2_cc_len,1)
        call para_pack_string(mp2_cc,mp2_cc_len)
        call para_pack_int(max_rec,1)
        call para_pack_int(irecord_tab,max_rec)
        call para_pack_int(ndiskc,1)
        call para_pack_real(diiscoeff,icalculate)
        call para_bcast_pack(TxCCSDInit)
      endif
      ipass_size=(npairs)/nslv
      ipass_reminder=mod(npairs,nslv)
      istart=1
      istop=0
        do isl=1,nslv
          istart=istop +1
          istop =istart+ipass_size-1
          if (ipass_reminder.gt.0) istop=istop+1
          ipass_reminder=ipass_reminder-1
          if (ipass_reminder.lt.0.and.ipass_size.eq.0) then
            istart=-1 ! means 'do nothing'
            istop =-1
          endif
          if (istop.gt.npairs) 
     *                STOP 'Cannot be true, par_diis_matrix'
          call para_recv(isl_tid,islgid,TxDiisReq)
          call para_initsend
          call para_pack_int(istart,1)
          call para_pack_int(istop ,1)
          call para_send_pack(islgid,TxDiisJob)
        enddo
        istart=0 ! means 'stop your work' (results are written to global AF)
        istop =0
        do isl=1,nslv
          call para_recv(isl_tid,islgid,TxDiisReq)
          call para_initsend
          call para_pack_int(istart,1)
          call para_pack_int(istop ,1)
          call para_send_pack(islgid,TxDiisJob)
        enddo
c
      if (singles.and.mp2_cc.eq.'ccc') then
        call dynamic_matdef('finalco','q',idimen,idimen)
        call dynamic_matdef('diiscoef','q',idimen,idimen)
        icoeff=mataddr('diiscoef')
        ifinalcoeff=mataddr('finalco')
        ij=npairs+1
        call matzero('finalco')
        do ii=1,icalculate
          k=ii+max_rec-icalculate
          irecord=(irecord_tab(k)-1)*(npairs+1)+ij
          if (af) then
            call fafRead(ndiskc,bl(icoeff),8,idimen*idimen,
     *                      1,irecord,info)
        if (info.lt.8*idimen*idimen)
     *  ERROR
          else
            call reader(ndiskc,idimen,irecord,bl(icoeff))
          endif
          call matadd1('diiscoef',diiscoeff(ii),'finalco')
        enddo
        do i=1,nval
          call writesingles(i,bl(ifinalcoeff+(i-1)*idimen))
        enddo
        call dynamic_matrem('diiscoef')
        call dynamic_matrem('finalco')
      endif
c
      end
C===============================================================================
      subroutine sla_new_coeff(singles,af,idimen,nval,list,
     *                          ipairimages,npairs,iMasID,mygid,mytid)
      use memory
      use newpara
      implicit none
      logical singles,af
      integer idimen,nval,list(*),ipairimages(*)
      integer npairs,iMasID,mygid,mytid
c
      integer icoeff,ifinalcoeff,mataddr,ij,i,j,ii,k,irecord,info
      character*256 mp2_cc
      integer icalculate,max_rec
      integer irecord_tab(DIIS_LENGTH),ndiskc
      real*8 diiscoeff(DIIS_LENGTH)
      integer mp2_cc_len,iextra,istart,istop,imastgid
c
      integer ix,iii
      real*8 xmax
c
      call sync_barrier
      do i=1,80
        mp2_cc(i:i)=' '
      enddo
      do i=1,DIIS_LENGTH
        diiscoeff(i)=0.0d0
      enddo
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(icalculate,1)
      if (icalculate.le.1) return
      call para_unpack_int(mp2_cc_len,1)
      call para_unpack_string(mp2_cc,mp2_cc_len)
      call para_unpack_int(max_rec,1)
      call para_unpack_int(irecord_tab,max_rec)
      call para_unpack_int(ndiskc,1)
      call para_unpack_real(diiscoeff,icalculate)
c
      call dynamic_matdef('finalco','q',idimen,idimen)
      call dynamic_matdef('diiscoef','q',idimen,idimen)
      icoeff=mataddr('diiscoef')
      ifinalcoeff=mataddr('finalco')
      if (singles) then
        iextra=1
      else
        iextra=0
      endif
      do
        call para_send(mytid,0,TxDiisReq)
        call para_recv_pack(imastgid,TxDiisJob)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop ,1)
        if (istart.eq.-1) cycle
        if (istart.eq. 0) exit
        ij=0
        do i=1,nval
          do j=1,i
            ij=ij+1
            if (ij.lt.istart) cycle
            if (ij.gt.istop ) exit
            if (list(ij).ne.1) cycle
            if (ipairimages(ij).eq.0) cycle
            call matzero('finalco')
            do ii=1,icalculate
              k=ii+max_rec-icalculate
              irecord=(irecord_tab(k)-1)*(npairs+iextra)+ij
              if (af) then
                call fafRead(ndiskc,bl(icoeff),8,idimen*idimen,
     *                          1,irecord,info)
            if (info.lt.8*idimen*idimen)
     *      ERROR
              else
                call reader(ndiskc,idimen,irecord,bl(icoeff))
              endif
              call matadd1('diiscoef',diiscoeff(ii),'finalco')
            enddo
c           do iii=0,idimen*idimen-1
c             bl(ifinalcoeff+iii)=ij
c           enddo
            call CoefWrite(i,j,ifinalcoeff)
c           call CoefRead('te',i,j,ix)
c           xmax=0d0
c           do iii=0,idimen*idimen-1
c           if (dabs(bl(ifinalcoeff+iii)-bl(ix+iii)).gt.xmax)
c    *              xmax=dabs(bl(ifinalcoeff+iii)-bl(ix+iii))
c           enddo
c           write(30+mygid,*) i,j,xmax
c           call flush(30+mygid)
c             if (i.eq.6.and.j.eq.6) then
c                call general_write1(97,bl(ifinalcoeff),idimen*idimen*8)
c                call flush(97)
c                call MPI_Abort(MPI_COMM_WORLD,1,i)
c             endif
          enddo
        enddo
      enddo
      call dynamic_matrem('diiscoef')
      call dynamic_matrem('finalco')
c
      end
C===============================================================================
      subroutine diis_slave_next(icode)
c icode=0, go ahead
c icode=1, step back
c icode=2, return
      use memory
      use newpara
      implicit none
      integer icode
      integer icode_loc
c
      integer info
c
      icode_loc=icode
      call sync_barrier
      call para_bcast(icode_loc,TxDIISFlow)
      end
C===============================================================================
      subroutine mas_trip(nval,idimen,npairs,ncf,nmo,vorb,epsi,res,ccsd,
     *                   irecadrx,ndiskx,lbinx,thresh,byt8,qcisd,af,
     *                   nslv1)
      use memory
      use newpara
      implicit none
      integer nval,idimen,npairs,ncf,nmo,irecadrx(*),ndiskx,lbinx
      integer nslv1
      logical vorb,ccsd,byt8,qcisd,af
      real*8 epsi(*),res,thresh
c
      real*8 :: total_energy
      real*8  energy,builder,abc,ene,t0,t1
      integer i,j,k,nrec,iWabc,itmp,len,len1,ndisk
      logical fileopen
      character*256 scrfile,filename
      real*8 :: time_i,time_o,time_mult,
     *                       etime_i,etime_o,etime_mult
      common /triples_times/  time_i, time_o, time_mult,
     *                       etime_i,etime_o,etime_mult
      integer :: mem,needmem
      integer i3ext_int,iamplit,iW_scratch,ij,iadr,ia,ic,irecords,n_int
      integer iamplits,ijk,istart,istop,ipass,itot_ijk,npass,ipass_size
      logical :: core,reset1,reset2,semicore
      real*8 ebuilder,eabc,eene,e0,e1,tot_e1,tot_e0,tot_t1,tot_t0
      integer ijk_loc,ivir_start,ivir_stop,ivirpass,nvirpass
      integer ivirpass_size,iscr,nrec_start,idisk_bytes,ii
      integer info,islvid,ireminder,i_add,idisk_space,islgid
      real*8 part_energy
      real*8 xbuilder, xabc, xene, xebuilder, xeabc, xeene, xtime_i
      real*8 xtime_o, xtime_mult, xetime_i, xetime_o, xetime_mult
      real*8 xtot_cpu,tot_cpu
c
      call secund(tot_t0)
      call elapsec(tot_e0)
      total_energy=0d0
      time_i=0d0
      time_o=0d0
      time_mult=0d0
      etime_i=0d0
      etime_o=0d0
      etime_mult=0d0
      ebuilder=0d0
      eabc=0d0
      eene=0d0
      itot_ijk=nval*(nval+1)*(nval+2)/6 - nval ! exclude (i,i,i)
c
c   Check if in-core algorithm is possible
c I need memory for all J(v,o)(v,v) integrals, all T(o,o)(v,v) amplitudes
c and two W(ijk)(abc) quantities
      needmem= (idimen*nval)*(idimen*idimen)
     *        +(npairs)*(idimen*idimen)
     *        +idimen*idimen*idimen
     *        +idimen*idimen*nval ! for amplitudes in second term from T
     *        +idimen*nval ! for three-internal integrals in the same term
      call dynamic_show_free(mem)
      core=.false.
      semicore=.false.
      if (mem.gt.needmem) core=.true.
c TEST
      core=.false.
      if (.not.core) print *, 'Would need ',needmem,
     *                        ' words of memory for in-core algorithm'
      if (core) print *, 'Use ',needmem,
     *                        ' words of memory for in-core algorithm'
      if (.not.core) then ! check if semi-core machanism is possible
        needmem= 3*idimen*idimen*idimen ! min 3 bunches of J(c,i) integrals
     *          +idimen*idimen*idimen
     *          +idimen*idimen*nval ! for amplitudes in second term from T
     *          +idimen*nval ! for three-internal integrals in the same term
      if (mem.gt.needmem) semicore=.true.
c TEST!
      semicore=.false.
      endif
      call flush(6)
      call para_initsend
      call para_pack_int(core,1)
      call para_pack_int(semicore,1)
      call para_bcast_pack(TxTriplesInit)
      do i=1,nslv
        call para_recv(idisk_space,islgid,TxTriplesRes)
        if (i.eq.1) idisk_bytes=idisk_space
        if (idisk_space.lt.idisk_bytes) idisk_bytes=idisk_space
      enddo
c Calculate how do we have to divide the virtual space
      call dynamic_mmark()
      call matmark()
c
c Triples loop:
c
c Divide ijk into passes:
c The total amount of ijk is 
      if (.not.core.and..not.semicore) then
c TEST!
        print *, 'Lowest free scratch space is: ', idisk_bytes
        idisk_bytes=idisk_bytes-idimen*idimen*idimen ! for one W
        npass=48*idimen*idimen*idimen*itot_ijk/idisk_bytes
        if (mod(48*idimen*idimen*idimen*itot_ijk,idisk_bytes).ne.0) then
          npass=npass+1
        endif
        do ! Make sure that for pass size (ii) resulting from npass the file
c            size will be smaller than idisk_bytes
          print *,'Npass: ',npass
          ii=itot_ijk/npass
          if (ii.eq.0) STOP 'To little disk space'
          if (mod(itot_ijk,npass).ne.0) ii=ii+1
          if (48*idimen*idimen*idimen*ii.lt.idisk_bytes) exit
          npass=npass+1
        enddo
      else
        npass=1
      endif
      ii=npass/nslv
      if (mod(npass,nslv).ne.0) ii=ii+1
      npass=nslv*ii
      ipass_size=itot_ijk/npass
      ireminder=mod(itot_ijk,npass)
      print *,'total ijk, ijk pass size, ijk passes: ',
     *         itot_ijk,ipass_size,npass
      call flush(6)
c The record formula for disk:
c (n-1)*6*idimen+(perm-1)*idimen+ic
c n is number of ijk in the present ipass, (ijk starts from 0 in a given pass)
c perm is number of a permutation, ic is a number of W column
      istop=0
      do ipass=1,npass
        if (ireminder.gt.0) then ! distribute reminder
          i_add=0                  ! normally we need to subtract 1 in
c                                    the istop = istart+ipass_size formula, 
c                                    if we won't do it it is equivalent to 
c                                    addition of one
          ireminder=ireminder-1
        else
          i_add=-1
        endif
        istart=istop+1
        istop = istart+ipass_size+i_add
        if (istop.gt.itot_ijk) istop=itot_ijk
        if (istart.gt.itot_ijk) STOP 'Wrong istart in triples'
        call para_recv(islvid,islgid,TxTriplesReq)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxTriplesJob)
        print *,'* * * * * * * * * * * IJK pass no.  ',
     *          ipass,',  ijk start, ijk stop: ',istart,istop
        call flush(6)
      enddo ! over ijk passes
      istart=0
      istop=0
      do i=1,nslv
        call para_recv(islvid,islgid,TxTriplesReq)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxTriplesJob)
        call para_recv_real(part_energy,1,islgid,TxTriplesRes)
        total_energy=total_energy+part_energy
      enddo
      print *,total_energy
      call secund(tot_t1)
      call elapsec(tot_e1)
c
c
      do i=1,nslv
        call para_recv_pack(islgid,TxTriplesRes)
        call para_unpack_real(xbuilder,1)
        call para_unpack_real(xabc,1)
        call para_unpack_real(xene,1)
        call para_unpack_real(xebuilder,1)
        call para_unpack_real(xeabc,1)
        call para_unpack_real(xeene,1)
        call para_unpack_real(xtime_i,1)
        call para_unpack_real(xtime_o,1)
        call para_unpack_real(xtime_mult,1)
        call para_unpack_real(xetime_i,1)
        call para_unpack_real(xetime_o,1)
        call para_unpack_real(xetime_mult,1)
        call para_unpack_real(xtot_cpu,1)
        builder=builder+xbuilder
        abc=abc+xabc
        ene=ene+xene
        ebuilder=ebuilder+xebuilder
        eabc=eabc+xeabc
        eene=eene+xeene
        time_i=time_i+xtime_i
        time_o=time_o+xtime_o
        time_mult=time_mult+xtime_mult
        etime_i=etime_i+xetime_i
        etime_o=etime_o+xetime_o
        etime_mult=etime_mult+xetime_mult
        tot_cpu=tot_cpu+xtot_cpu
      enddo
c
c
      print '(A)', '             Build amplitudes,            '//
     *             'Combine,             Energy:'
      print '(A10,3F20.1)','Sum CPU:',builder,abc,ene
      print '(A10,3F20.1)','Sum ELAPS:',ebuilder,eabc,eene
      print '(A)', 'WBuilder      Input,   Output,     Mult: '
      print '(A10,3F10.1)','Sum CPU:',time_i,time_o,time_mult
      print '(A10,3F10.1)','Sum ELAPS:',etime_i,etime_o,etime_mult
      print '(A21,F10.1)','Total triples CPU:  ',tot_cpu
      print '(A21,F10.1)','Total triples ELAPS:',tot_e1-tot_e0
      res=total_energy
      call matremark()
      call dynamic_retmark()
      end
C===============================================================================
      subroutine sla_trip(nval,idimen,npairs,ncf,nmo,vorb,epsi,ccsd,
     *                   irecadrx,ndiskx,lbinx,thresh,byt8,qcisd,af,
     *                   nslv1,iMasID,mygid,mytid)
      use memory
      use newpara
      implicit none
      integer nval,idimen,npairs,ncf,nmo,irecadrx(*),ndiskx,lbinx
      integer nslv1,iMasID,mygid,mytid
      logical vorb,ccsd,byt8,qcisd,af
      real*8 epsi(*),thresh
c
      real*8 :: total_energy
      real*8  energy,builder,abc,ene,t0,t1
      integer i,j,k,nrec,iWabc,itmp,len,len1,ndisk
      logical fileopen
      character*256 scrfile,filename
      real*8 :: time_i,time_o,time_mult,
     *                       etime_i,etime_o,etime_mult
      common /triples_times/  time_i, time_o, time_mult,
     *                       etime_i,etime_o,etime_mult
      integer :: irecadrtt,ndisktt,lbintt,mem,needmem
      integer i3ext_int,iamplit,iW_scratch,ij,iadr,ia,ic,irecords,n_int
      integer iamplits,ijk,istart,istop,ipass,npass,ipass_size
      logical :: core,reset1,reset2,semicore
      real*8 ebuilder,eabc,eene,e0,e1,tot_e1,tot_e0,tot_t1,tot_t0
      integer ijk_loc,ivir_start,ivir_stop,ivirpass,nvirpass
      integer ivirpass_size,iscr,nrec_start,idisk_bytes,ii
      integer ilast,ilastpass,ilast_start,ilast_stop,ilastpass_size
      integer k_ij,k_ik,k_jk,ndisk1,iaddr,ising
      character*2 ch_mygid
      integer info,idisk_space,imastgid
c
      irecadrtt=0
      ndisktt=0
      lbintt=0
c
      call secund(tot_t0)
      call elapsec(tot_e0)
      total_energy=0d0
      time_i=0d0
      time_o=0d0
      time_mult=0d0
      etime_i=0d0
      etime_o=0d0
      etime_mult=0d0
      ebuilder=0d0
      eabc=0d0
      eene=0d0
c
      call para_recv_bcastpack(TxTriplesInit)
      call para_unpack_int(core,1)
      call para_unpack_int(semicore,1)
      call getchval('scrf',scrfile)
      call partition_size(scrfile,idisk_space)
      call para_send(idisk_space,0,TxTriplesRes)
      nvirpass=1
      if (.not.core.and..not.semicore) then
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len)
        write(ch_mygid,'(I2.2)') mygid
        filename=scrfile(1:len)//'.triples-'//ch_mygid
        len1=len+11
        call s_open(filename(1:len1),ndisk,idimen*idimen*8)
c
        call dynamic_show_free(mem)
        mem=mem-3*idimen*idimen*nval   ! for amplitudes
        mem=mem-2*idimen*nval            ! 3-int. integrals & copy
        n_int=3
        do
          print *,'Try if it will go with ',nvirpass,' passes.'
          if (4*idimen*idimen*idimen/nvirpass.lt.mem) then
            exit
          endif
          if (nvirpass.gt.idimen) STOP 'Memory exhausted'
          nvirpass=nvirpass+1
        enddo
c TEST!!!
c       nvirpass=4
        print *,'Disk algorithm with ',nvirpass,' passes over virt.'
     *     //   ' space.'
      endif !  if (.not.core.and..not.semicore)
      ivirpass_size=idimen/nvirpass
      if (mod(idimen,nvirpass).ne.0) ivirpass_size=ivirpass_size+1
c
      call dynamic_mmark()
      call matmark()
      if (core) then
        call dynamic_getmem(idimen*nval*idimen*idimen,i3ext_int)
        call dynamic_getmem(npairs*idimen*idimen,iamplit)
        call dynamic_getmem(idimen*idimen*idimen,iW_scratch)
c Store integrals:
        ij=-1
        do ic=1,idimen
          do ia=1,nval
            ij=ij+1 ! starts from 0 !
            iadr=i3ext_int+ij*idimen*idimen
            call ExtrOne(ic,    ia,     irecadrtt,  npairs,    ndisktt,
     *                   ncf,  lbintt,thresh,     byt8,      'tt',
     *                   'mo', nmo,   vorb,       bl(iadr))
          enddo
        enddo
        call dynamic_lock(bl(i3ext_int),ij)
c Store amplitudes:
        ij=-1
        do i=1,nval
          do j=1,i
            ij=ij+1 ! starts from 0 !
            iadr=iamplit+ij*idimen*idimen
            call CoefRead('tt',i,j,itmp)
            call tfer(bl(itmp),bl(iadr),idimen*idimen)
          enddo
        enddo
        call dynamic_lock(bl(iamplit),ij)
      else if (semicore) then
        call dynamic_getmem(idimen*idimen*idimen,iW_scratch)
        needmem=idimen*nval ! for three-internal integrals in the same term
        call dynamic_getmem(nval,irecords) ! info about integrals location
        call dynamic_getmem(nval,iamplits) ! info about amplitudes
        call dynamic_show_free(mem)
        n_int=(mem-needmem)/(idimen*idimen*idimen+nval*idimen*idimen)
        if (n_int.lt.3) STOP 'Expected n_int>=3 in triples'
        if (n_int.gt.nval) n_int=nval
        call dynamic_getmem(n_int*idimen*idimen*idimen,i3ext_int)
        call dynamic_getmem(n_int*nval*idimen*idimen,iamplit)
        print *,'Semicore algorithm with cahing of ',n_int,' integrals.'
        call flush(6)
      endif
c
c Triples loop:
c
      reset1=.true.
      reset2=.true.
c Divide ijk into passes:
c The total amount of ijk is 
c The record formula for disk:
c (n-1)*6*idimen+(perm-1)*idimen+ic
c n is number of ijk in the present ipass, (ijk starts from 0 in a given pass)
c perm is number of a permutation, ic is a number of W column
      do ! over ijk passes, infinite, driven by master
      call para_send(mytid,0,TxTriplesReq)
      call para_recv_pack(imastgid,TxTriplesJob)
      call para_unpack_int(istart,1)
      call para_unpack_int(istop ,1)
      if (istart.eq.0.and.istop.eq.0) then
        call para_send_real(total_energy,1,0,TxTriplesRes)
        exit
      endif
c
      if (.not.core.and..not.semicore) then
        call dynamic_mmark()
        call dynamic_getmem(idimen*idimen*ivirpass_size,iW_scratch)
        call dynamic_getmem(3*ivirpass_size*idimen*idimen,i3ext_int)
        call dynamic_getmem(3*nval*idimen*idimen,iamplit)
        call dynamic_getmem(nval,irecords) ! info about integrals location
        call dynamic_getmem(nval,iamplits) ! info about amplitudes
      endif
c
      do ivirpass=1,nvirpass
        ivir_start=(ivirpass-1)*ivirpass_size+1
        ivir_stop =(ivirpass  )*ivirpass_size
        if (ivir_stop .gt.idimen) ivir_stop=idimen
        if (ivir_start.gt.idimen) STOP 'Suprise in triples ;)'
        print *,'Virt. passes. start,stop: ',ivir_start,ivir_stop
        if (.not.core.and..not.semicore) then
          reset1=.true.
          reset2=.true.
        endif
c
      ijk_loc=0
      ijk=0
      do i=1,nval
        do j=1,i
          do k=1,j
            if (i.eq.j.and.j.eq.k) cycle
            ijk=ijk+1
            if (ijk.lt.istart) cycle
            if (ijk.gt.istop) exit
            ijk_loc=ijk_loc+1
            print *,'Triples progress... ',i,j,k
            call flush(6)
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) then
              call zeroit(bl(iW_scratch),idimen*idimen*idimen)
            endif
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c
c  Store as much J(c,i) integrals as possible, max is n_int. I need to
c put priority on the integrals from i and j indices (they HAVE TO be
c kept, because they are used repeatedly for all "k"), then I need to
c store as much of low indices as possible (most frequently read) and also
c I HAVE to gather integrals with k index. So, I have to keep i,j,k and then
c as much as possible of all the lowest indices from 1 to nval range. I have
c bl(irecords) table which resolves index to the first record number from the
c "c" group (all "c" are kept for semidirect).
c This algorithm works so well, that I will probably remove "core"
c option in the future - semicore with n_int==nval gives the same
c performance, and locate_... subroutine simply returns as they do not
c have any work to do but only initial load of all the integrals and
c amplitudes
c
            if (.not.core) then
            call secund(t0)
            call elapsec(e0)
            call locate_integrals(i, j, k, i3ext_int, bl(irecords),
     *                           reset1,n_int, idimen,irecadrtt,npairs,
     *                           ndisktt,ncf,lbintt,thresh,byt8,nmo,
     *                           vorb,nval,ivir_start,ivir_stop)
            call locate_amplits(i, j, k, iamplit, bl(iamplits),
     *                          reset2,n_int, idimen,irecadrtt,npairs,
     *                          ndisktt,ncf,lbintt,thresh,byt8,nmo,
     *                          vorb,nval)
            call secund(t1)
            call elapsec(e1)
            time_i=time_i+t1-t0
            etime_i=etime_i+e1-e0
            endif
c Build \tilda W for six possible permutations:
c ijk
            nrec=(ijk_loc-1)*6*idimen
            call secund(t0)
            call elapsec(e0)
            call WTildaBuilderWriter(i,j,k,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                            iW_scratch,core,semicore,bl(irecords),
     *                               bl(iamplits),ivir_start,ivir_stop)
            call secund(t1)
            call elapsec(e1)
            builder=builder+t1-t0
            ebuilder=ebuilder+e1-e0
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) 
     *          call sorter(bl(iW_scratch),idimen,2,3)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c ikj
            nrec=(ijk_loc-1)*6*idimen+1*idimen
            call secund(t0)
            call elapsec(e0)
            call WTildaBuilderWriter(i,k,j,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                            iW_scratch,core,semicore,bl(irecords),
     *                               bl(iamplits),ivir_start,ivir_stop)
            call secund(t1)
            call elapsec(e1)
            builder=builder+t1-t0
            ebuilder=ebuilder+e1-e0
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) 
     *          call sorter(bl(iW_scratch),idimen,1,2)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c kij
            nrec=(ijk_loc-1)*6*idimen+2*idimen
            call secund(t0)
            call elapsec(e0)
            call WTildaBuilderWriter(k,i,j,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                            iW_scratch,core,semicore,bl(irecords),
     *                               bl(iamplits),ivir_start,ivir_stop)
            call secund(t1)
            call elapsec(e1)
            builder=builder+t1-t0
            ebuilder=ebuilder+e1-e0
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) 
     *          call sorter(bl(iW_scratch),idimen,2,3)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c kji
            nrec=(ijk_loc-1)*6*idimen+3*idimen
            call secund(t0)
            call elapsec(e0)
            call WTildaBuilderWriter(k,j,i,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                            iW_scratch,core,semicore,bl(irecords),
     *                               bl(iamplits),ivir_start,ivir_stop)
            call secund(t1)
            call elapsec(e1)
            builder=builder+t1-t0
            ebuilder=ebuilder+e1-e0
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) 
     *          call sorter(bl(iW_scratch),idimen,1,2)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c jki
            nrec=(ijk_loc-1)*6*idimen+4*idimen
            call secund(t0)
            call elapsec(e0)
            call WTildaBuilderWriter(j,k,i,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                            iW_scratch,core,semicore,bl(irecords),
     *                               bl(iamplits),ivir_start,ivir_stop)
            call secund(t1)
            call elapsec(e1)
            builder=builder+t1-t0
            ebuilder=ebuilder+e1-e0
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) 
     *          call sorter(bl(iW_scratch),idimen,2,3)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c jik
            nrec=(ijk_loc-1)*6*idimen+5*idimen
            call secund(t0)
            call elapsec(e0)
            call WTildaBuilderWriter(j,i,k,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                            iW_scratch,core,semicore,bl(irecords),
     *                               bl(iamplits),ivir_start,ivir_stop)
            call secund(t1)
            call elapsec(e1)
            builder=builder+t1-t0
            ebuilder=ebuilder+e1-e0
            call secund(t0)
            call elapsec(e0)
            if (core.or.semicore) 
     *          call sorter(bl(iW_scratch),idimen,1,2)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
            if (core.or.semicore) then
              iWabc=iW_scratch
c Calculate energy:
            call secund(t0)
            call elapsec(e0)
            call EnergyContribution(i,j,k,idimen,nval,bl(iWabc),epsi,
     *                              energy,ccsd,irecadrx,npairs,ndiskx,
     *                              ncf,lbinx,thresh,byt8,nmo,vorb,
     *                              qcisd)
            call secund(t1)
            call elapsec(e1)
            ene=ene+t1-t0
            eene=eene+e1-e0
            total_energy=total_energy+energy
            endif
          enddo
        enddo
      enddo
      enddo ! over virt passes
c
      if (.not.core.and..not.semicore) then
        call dynamic_retmark()
        call dynamic_mmark()
        call dynamic_getmem(idimen*idimen,iscr)
c Check if the whole W can be kept in memory
        call dynamic_show_free(mem)
        mem=mem-3*idimen*idimen ! singles
        ilastpass=0
        do ii=1,10
          ilastpass=ilastpass+1
          ilastpass_size=idimen/ilastpass
          if (mod(idimen,ilastpass).ne.0)ilastpass_size=ilastpass_size+1
          if (idimen*idimen*ilastpass_size.le.mem) exit
        enddo
c TEST!!!
c       ilastpass=2
        ilastpass_size=idimen/ilastpass
        if (mod(idimen,ilastpass).ne.0) ilastpass_size=ilastpass_size+1
        if (idimen*idimen*ilastpass_size.gt.mem) then
          STOP 'Too little memory to complete the triples calculation'
        endif
        print *,'''W'' assembly will be done in ',ilastpass,' passes'
        call flush(6)
c
        call dynamic_getmem(idimen*idimen*ilastpass_size,iWabc)
c
        if (ilastpass.gt.1) then
          call dynamic_getmem(idimen*idimen,k_ij)
          call dynamic_getmem(idimen*idimen,k_ik)
          call dynamic_getmem(idimen*idimen,k_jk)
          filename=scrfile(1:len)//'.Ws-'//ch_mygid
          len1=len+6
          call s_open(filename(1:len1),ndisk1,idimen*idimen*8)
        endif
c
      ijk_loc=0
      ijk=0
      do i=1,nval
        do j=1,i
          do k=1,j
           if (i.eq.j.and.j.eq.k) cycle
           ijk=ijk+1
           if (ijk.lt.istart) cycle
           if (ijk.gt.istop) exit
           ijk_loc=ijk_loc+1
           nrec_start=(ijk_loc-1)*6*idimen
           do ilast=1,ilastpass
            ilast_start=(ilast-1)*ilastpass_size+1
            ilast_stop =(ilast  )*ilastpass_size
            if (ilast_stop.gt.idimen) ilast_stop=idimen
            if (ilast_start.gt.idimen) STOP 'ilast_start.gt.idimen'
            call secund(t0)
            call elapsec(e0)
            call Wabc(ndisk,idimen,bl(iscr),bl(iWabc),nrec_start,
     *                ilast_start,ilast_stop)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
c
           call secund(t0)
           call elapsec(e0)
            if (ilastpass.gt.1) then
            do ii=ilast_start,ilast_stop
              iaddr=iWabc+(ii-ilast_start)*idimen*idimen
              call s_write(ndisk1,ii,bl(iaddr))
            enddo
            endif
           enddo
           call secund(t1)
           call elapsec(e1)
           time_o=time_o+t1-t0
           etime_o=etime_o+e1-e0
c
           if (ilastpass.eq.1) then
            call secund(t0)
            call elapsec(e0)
            call EnergyContribution(i,j,k,idimen,nval,bl(iWabc),epsi,
     *                              energy,ccsd,irecadrx,npairs,ndiskx,
     *                              ncf,lbinx,thresh,byt8,nmo,vorb,
     *                              qcisd)
            call secund(t1)
            call elapsec(e1)
            ene=ene+t1-t0
            eene=eene+e1-e0
            total_energy=total_energy+energy
           else
           if (ccsd.or.qcisd) then
             call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *                    ncf,  lbinx, thresh,     byt8,      'x',
     *                    'mo', nmo,   vorb,       bl(k_ij))
             call ExtrOne(i,    k,     irecadrx,   npairs,    ndiskx,
     *                    ncf,  lbinx, thresh,     byt8,      'x',
     *                    'mo', nmo,   vorb,       bl(k_ik))
             call ExtrOne(j,    k,     irecadrx,   npairs,    ndiskx,
     *                    ncf,  lbinx, thresh,     byt8,      'x',
     *                    'mo', nmo,   vorb,       bl(k_jk))
           call pointersingles(ising)
           endif
           do ilast=1,ilastpass
            ilast_start=(ilast-1)*ilastpass_size+1
            ilast_stop =(ilast  )*ilastpass_size
            if (ilast_stop.gt.idimen) ilast_stop=idimen
            if (ilast_start.gt.idimen) STOP 'ilast_start.gt.idimen'
            call secund(t0)
            call elapsec(e0)
            call EnergyContrPartial(i,j,k,idimen,nval,bl(iWabc),epsi,
     *                              energy,ccsd,irecadrx,npairs,ndiskx,
     *                              ncf,lbinx,thresh,byt8,nmo,vorb,
     *                              qcisd,ilast_start,ilast_stop,ndisk1,
     *                              bl(iscr),bl(k_ij),bl(k_ik),bl(k_jk),
     *                              bl(ising))
            call secund(t1)
            call elapsec(e1)
            ene=ene+t1-t0
            eene=eene+e1-e0
            total_energy=total_energy+energy
           enddo
           endif
          enddo
        enddo
      enddo
        if (ilastpass.gt.1) then
          call s_close(ndisk1)
        endif
        call dynamic_retmark()
      endif
c
      enddo ! over ijk passes
      print *,total_energy
      call secund(tot_t1)
      call elapsec(tot_e1)
c
      call para_initsend
      call para_pack_real(builder,1)
      call para_pack_real(abc,1)
      call para_pack_real(ene,1)
      call para_pack_real(ebuilder,1)
      call para_pack_real(eabc,1)
      call para_pack_real(eene,1)
      call para_pack_real(time_i,1)
      call para_pack_real(time_o,1)
      call para_pack_real(time_mult,1)
      call para_pack_real(etime_i,1)
      call para_pack_real(etime_o,1)
      call para_pack_real(etime_mult,1)
      call para_pack_real(tot_t1-tot_t0,1)
      call para_send_pack(0,TxTriplesRes)
      call matremark()
      call dynamic_retmark()
      call s_close(ndisk)
      end
C===============================================================================
      subroutine CCAOncemast(irecadrx, npairs, nval, ndiskx, ncf,
     *                   lbinx,    thresh, doit, byt8,   ioverlap,
     *                   nmo,      vorb,nslv1,nbf)
c This subroutine calculates A matrix, as it is defined in 
c J. Chem. Phys. 81, 1901 (1984) - formula no. (38)
c Parameters:
c   irecadrx - the record info for file containing Kext bins. Note, that 
c              this array has different size and different information is
c              contained in it in single processor mode and in parallel mode.
c   npairs   - number of correlated pairs
c   nval     - number of correlated orbitals
c   ndiskx   - file unit number where Kext bins are kept
c   ncf      - number of contracted functions
c   lbinx    - Kext bin length
c   thresh   - integral threshold
c   doit     - it is passed to common /CCAOnceCom/, but not used anywhere.
c   byt8     - .true.: double precision bins are used. .false.: int4+int1
c              bins are used
c   ioverlap - the overlap matrix bl index
c   nmo      - number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      use newpara
      use ccounters
      implicit none
      integer irecadrx(*),npairs,nval,ndiskx,ncf,lbinx,iA,ioverlap,nmo
      integer nslv1,nbf
      real*8 thresh
      real*8 half
      integer iexchAO,icoeffAO,mataddr,idimen,iwork33
      integer info,i,islv_tid,istart,istop,ipasssize,ipass,ireminder
      integer npasses,islgid
      logical doit,idoit,byt8,vorb
      parameter (half=0.5d0)
      common /CCAOnceCom/ iA,idoit
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ic_CCAOmas=ic_CCAOmas+1
      idoit=doit
      if (ic_CCAOmas.eq.1) call dynamic_getmem(idimen*idimen,iA)
      call mmark
      call dynamic_mmark
      call matmark
      call matconn('A','q',idimen,idimen,iA)
      call matzero('A')
      call dynamic_matdef('work33','q',idimen,idimen)
      iwork33=mataddr('work33')
c Divide npairs into groups, preferably the same as number of slaves.
      npasses=nslv
c let one slave to work on matrix mult. no more than 20 seconds
      do i=1,100
        ipasssize=npairs/npasses
        if (dble(ipasssize)*dble(idimen)**3d0/(1d9).gt.2d1) then
          npasses=nslv*i
        else
          exit
        endif
      enddo
c A is calculated according to the recipe taken from Pualy article. It is at
c the end transposed and multiplied by overlap in order to be identical
c with A from Werner article.
      ireminder=mod(npairs,npasses)
      ipasssize=npairs/npasses
      istart=1
      istop=0
      do ipass=1,npasses
        istart=istop+1
        istop=istart+ipasssize-1
        if (ireminder.gt.0) then
          istop=istop+1
          ireminder=ireminder-1
        else if (ipasssize.eq.0) then
          exit
        endif
        if (istop.gt.npairs.or.istart.gt.npairs) STOP 'CCAOnce, error'
        call para_recv(islv_tid,islgid,TxCCSDReq1)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob1)
      enddo
      istart=0
      istop=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq1)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob1)
        call para_recv_real(bl(iwork33),idimen*idimen,islgid,TxCCSDRes1)
        call matadd('work33','A')
      enddo
! A is changed in order to be the same as in Werner article
      call matpose('A')
      if (vorb) then
        continue
c       call matcopy('work33','A')
      else
        call matcopy('A','work33')
        call matconn('overl','q',ncf,ncf,ioverlap)
        call matmmult('work33','overl','A')
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
C==============================================================================
      subroutine CCAOnceslav(irecadrx, npairs, nval, ndiskx, ncf,
     *                   lbinx,    thresh, doit, byt8,   ioverlap,
     *                   nmo,      vorb,nslv1,iMasID,mygid,mytid,
     *                   nbf)
c This subroutine calculates A matrix, as it is defined in 
c J. Chem. Phys. 81, 1901 (1984) - formula no. (38)
c Parameters:
c   irecadrx - the record info for file containing Kext bins. Note, that 
c              this array has different size and different information is
c              contained in it in single processor mode and in parallel mode.
c   npairs   - number of correlated pairs
c   nval     - number of correlated orbitals
c   ndiskx   - file unit number where Kext bins are kept
c   ncf      - number of contracted functions
c   lbinx    - Kext bin length
c   thresh   - integral threshold
c   doit     - it is passed to common /CCAOnceCom/, but not used anywhere.
c   byt8     - .true.: double precision bins are used. .false.: int4+int1
c              bins are used
c   ioverlap - the overlap matrix bl index
c   nmo      - number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      use newpara
      implicit none
      integer irecadrx(*),npairs,nval,ndiskx,ncf,lbinx,iA,ioverlap,nmo
      integer nslv1,iMasID,mygid,mytid,imastgid,nbf
      real*8 thresh
      real*8 half
      integer iexchAO,icoeffAO,k,l,mataddr,idimen,info,istart,istop,kl
      logical doit,idoit,byt8,vorb
      parameter (half=0.5d0)
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('A','q',idimen,idimen)
      call matzero('A')
      iA=mataddr('A')
      call dynamic_matdef('exchAO','q',idimen,idimen)
      call dynamic_matdef('coetilda','q',idimen,idimen)
      call dynamic_matdef('work33','q',idimen,idimen)
      iexchAO=mataddr('exchAO')
c A is calculated according to the recipe taken from Pualy article. It is at
c the end transposed and multiplied by overlap in order to be identical
c with A from Werner article.
      do
        call para_initsend
        call para_send(mytid,0,TxCCSDReq1)
        call para_recv_pack(imastgid,TxCCSDJob1)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop,1)
        if (istart.eq.0.and.istop.eq.0) then
          call para_send_real(bl(iA),idimen*idimen,0,TxCCSDRes1)
          exit
        endif
        kl=0
        do k=1,nval
          do l=1,k
            kl=kl+1
            if (kl.lt.istart) cycle
            if (kl.gt.istop)  exit
            call ExtrOne(k,    l,     irecadrx,   npairs,    ndiskx,
     *                   ncf,  lbinx, thresh,     byt8,      'x',
     *                   'mo', nmo,   vorb,       bl(iexchAO))
            call CoefRead('tt',k,l,icoeffAO)
            call matconn('coeAO','q',idimen,idimen,icoeffAO)
            call matzero('coetilda')
            call matadd1('coeAO',+4.0d0,'coetilda')
            call matpose('coeAO')
            call matadd1('coeAO',-2.0d0,'coetilda')
            call matpose('coeAO')
            call matdisc('coeAO')
            if (k.eq.l) call matscal('coetilda',half)
            call matmmul2('coetilda','exchAO','work33','t','n','n')
            call matadd1('work33',half,'A')
            call matmmul2('coetilda','exchAO','work33','n','t','n')
            call matadd1('work33',half,'A')
          enddo
        enddo
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine Tl_gen_mast(ncf,nval,nmo,vorb,ioverl,imatrix,nslv1,nbf)
c This suroutine calculates the T x l term, according to 
c J. Chem. Phys. 104, 6286, (1996) it is part of formula no 44. 
c The output is written using writeS subroutine to storage which stores
c later the whole (44) term.
c Parameters:
cI  ncf    - the number of contracted functions
cI  nval   - the number of correlated orbitals
cI  nmo    - the number of occupied molecular orbitals
cI  vorb   - .true.: the dimension of virtual space is equal to no of virt.
c            orbitals, and MO form of integral matrices is used, otherwise,
c            the virtual space has ncf dimension and all calculations are 
c            performed for matrices in AO form.
      use memory
      use newpara
      implicit none
      integer ncf,nval,nmo,imatrix,ioverl,nslv1,nbf
      logical vorb
      integer k,l,i,idoubles,iwork11,iwork22,mataddr,idimen
      integer iresult,iresult1,iwork,npairs,ireminder,ipass
      integer info,islv_tid,istart,istop,npasses,ipasssize,islgid
c
      call mmark
      call dynamic_mmark
      call matmark
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      npairs=(nval*(nval+1))/2
      call matconn('fin_res','r',idimen,nval,imatrix)
      call matconn('xoverlap','q',ncf,ncf,ioverl)
      call dynamic_matdef('result','r',idimen,nval)
      iresult=mataddr('result')
      call matzero('result')
      call dynamic_matdef('result1','r',idimen,nval)
      iresult1=mataddr('result1')
c Divide npairs into groups, preferably the same as number of slaves.
      npasses=nslv
c let one slave to work on matrix mult. no more than 20 seconds
      do i=1,100
        ipasssize=npairs/npasses
        if (dble(ipasssize)*dble(idimen)**3d0/(1d9).gt.2d1) then
          npasses=nslv*i
        else
          exit
        endif
      enddo
c A is calculated according to the recipe taken from Pualy article. It is at
c the end transposed and multiplied by overlap in order to be identical
c with A from Werner article.
      ireminder=mod(npairs,npasses)
      ipasssize=npairs/npasses
      istart=1
      istop=0
      do ipass=1,npasses
        istart=istop+1
        istop=istart+ipasssize-1
        if (ireminder.gt.0) then
          istop=istop+1
          ireminder=ireminder-1
        else if (ipasssize.eq.0) then
          exit
        endif
        if (istop.gt.npairs.or.istart.gt.npairs) STOP 'Tl_gen, error'
        call para_recv(islv_tid,islgid,TxCCSDReq2)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob2)
      enddo
      istart=0
      istop=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq2)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob2)
        call para_recv_real(bl(iresult1),idimen*nval,islgid,TxCCSDRes2)
        call matadd('result1','result')
      enddo
      if (vorb) then
        call matadd1('result',-1d0,'fin_res')
      else
        call matscal('result',-1.0d0)
        call matmmul2('xoverlap','result','fin_res','n','n','a')
      endif
      call matremark
      call retmark
      call dynamic_retmark
      call para_barrier
      end
c==============================================================================
      subroutine Tl_gen_slav(ncf,nval,nmo,vorb,nslv1,iMasID,mygid,mytid,
     *                       nbf)
c This suroutine calculates the T x l term, according to 
c J. Chem. Phys. 104, 6286, (1996) it is part of formula no 44. 
c The output is written using writeS subroutine to storage which stores
c later the whole (44) term.
c Parameters:
cI  ncf    - the number of contracted functions
cI  nval   - the number of correlated orbitals
cI  nmo    - the number of occupied molecular orbitals
cI  vorb   - .true.: the dimension of virtual space is equal to no of virt.
c            orbitals, and MO form of integral matrices is used, otherwise,
c            the virtual space has ncf dimension and all calculations are 
c            performed for matrices in AO form.
      use memory
      use newpara
      implicit none
      integer ncf,nval,nmo,nslv1,iMasID,mygid,mytid,nbf
      logical vorb
      integer k,l,i,idoubles,iwork11,iwork22,mataddr,idimen
      integer iresult,iwork,info,istart,istop,kl,imastgid
c
      call mmark
      call dynamic_mmark
      call matmark
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call dynamic_matdef('result','r',idimen,nval)
      iresult=mataddr('result')
      call zeroit(bl(iresult),idimen*nval)
c
      call dynamic_matdef('work11','r',idimen,nval)
      iwork11=mataddr('work11')
      do
        call para_send(mytid,0,TxCCSDReq2)
        call para_recv_pack(imastgid,TxCCSDJob2)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop,1)
        if (istart.eq.0.and.istop.eq.0) then
          call para_send_real(bl(iresult),idimen*nval,0,TxCCSDRes2)
          exit
        endif
        kl=0
        do k=1,nval
          do l=1,k
            kl=kl+1
            if (kl.lt.istart) cycle
            if (kl.gt.istop)  exit
            call CoefRead('tt',l,k,idoubles)
            call matconn('doubl_am','q',idimen,idimen,idoubles)
            do i=1,nval
              iwork=iwork11+(i-1)*idimen
              call xlijk(k,l,i,bl(iwork))
            enddo
            call matmmul2('doubl_am','work11','result','n','n','a')
            if (k.ne.l) then
              do i=1,nval
                iwork=iwork11+(i-1)*idimen
                call xlijk(l,k,i,bl(iwork))
              enddo
              call matmmul2('doubl_am','work11','result','t','n','a')
            endif
            call matdisc('doubl_am')
          enddo
        enddo
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      call para_barrier
      end
C==============================================================================
      subroutine Lt_gen_mast(idimen,  nval,  npairs,  ncf,  vorb,
     *                       iresult,nslv1)
      use memory
      use newpara
      implicit none
c
      integer idimen,nval,npairs,ncf,iresult,nslv1
      logical vorb
c
      integer mataddr,iloc_res,npasses,ipasssize,i,info
      integer ireminder,istart,istop,ipass,islv_tid,islgid
c
      call dynamic_mmark
      call matmark
c
      call dynamic_matdef('loc_res','r',idimen,nval)
      iloc_res=mataddr('loc_res')
      call matconn('ext_res','r',idimen,nval,iresult)
c
c Divide npairs into groups, preferably the same as number of slaves.
      npasses=nslv
c let one slave to work on matrix mult. no more than 20 seconds
      do i=1,100
        ipasssize=npairs/npasses
        if (dble(ipasssize)*dble(idimen)**3d0/(1d9).gt.2d1) then
          npasses=nslv*i
        else
          exit
        endif
      enddo
c A is calculated according to the recipe taken from Pualy article. It is at
c the end transposed and multiplied by overlap in order to be identical
c with A from Werner article.
      ireminder=mod(npairs,npasses)
      ipasssize=npairs/npasses
      istart=1
      istop=0
      do ipass=1,npasses
        istart=istop+1
        istop=istart+ipasssize-1
        if (ireminder.gt.0) then
          istop=istop+1
          ireminder=ireminder-1
        else if (ipasssize.eq.0) then
          exit
        endif
        if (istop.gt.npairs.or.istart.gt.npairs) STOP 'Lt_gen error'
        call para_recv(islv_tid,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob3)
      enddo
      istart=0
      istop=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob3)
        call para_recv_real(bl(iloc_res),idimen*nval,islgid,TxCCSDRes3)
        call matadd('loc_res','ext_res')
      enddo
      call matremark
      call dynamic_retmark
      end
C==============================================================================
      subroutine Lt_gen_slav(idimen,nval,   irecadrx,npairs,ndiskx,
     *                       ncf,   lbinx,  thresh,  byt8,  nmo,
     *                       vorb,  nslv1,  iMasID,  mygid, mytid)
      use memory
      use newpara
      implicit none
c
      integer idimen,nval,irecadrx(*),npairs,ndiskx,ncf,lbinx
      integer nslv1,iMasID,mygid,mytid
      real*8 thresh
      logical byt8,vorb
      integer nmo
c
      integer isingles,iL,k,l,iwork11,ising,ires,mataddr,iloc_res
      integer info,istart,istop,kl,imastgid
c
      call dynamic_mmark
      call matmark
c
      call dynamic_matdef('L','q',idimen,idimen)
      iL=mataddr('L')
      call dynamic_matdef('work11','q',idimen,idimen)
      iwork11=mataddr('work11')
      call dynamic_matdef('loc_res','r',idimen,nval)
      call matzero('loc_res')
      iloc_res=mataddr('loc_res')
      call pointersingles(isingles)
      do
        call para_send(mytid,0,TxCCSDReq3)
        call para_recv_pack(imastgid,TxCCSDJob3)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop,1)
        if (istart.eq.0.and.istop.eq.0) then
          call para_send_real(bl(iloc_res),idimen*nval,0,TxCCSDRes3)
          exit
        endif
      kl=0
      do k=1,nval
        do l=1,k
            kl=kl+1
            if (kl.lt.istart) cycle
            if (kl.gt.istop)  exit
            call ExtrOne(k,    l,     irecadrx,     npairs,     ndiskx,
     *                   ncf,  lbinx, thresh,       byt8,       'x',
     *                   'mo', nmo,   vorb,         bl(iwork11))
            call CCLoperator(bl(iwork11),idimen,bl(iL))
c
            ising=isingles+(l-1)*idimen
            ires =iloc_res+(k-1)*idimen
            call matconn('singl1','r',idimen,1,ising)
            call matconn('res1','r',idimen,1,ires)
            call matmmul2('L','singl1','res1','n','n','a')
            call matdisc('res1')
            call matdisc('singl1')
            if (k.ne.l) then
              ising=isingles+(k-1)*idimen
              ires =iloc_res+(l-1)*idimen
              call matconn('singl1','r',idimen,1,ising)
              call matconn('res1','r',idimen,1,ires)
              call matmmul2('L','singl1','res1','t','n','a')
              call matdisc('res1')
              call matdisc('singl1')
            endif
          enddo
        enddo
      enddo
      call matremark
      call dynamic_retmark
      end
C==============================================================================
      subroutine Tr_gen_mast(idimen,nval,ncf,npairs,vorb,
     *                       ioverlap,nslv1,iresult)
c This subroutine calculates the second term of (38) formula taken from
c J. Chem. Phys. 104, 6286, (1996)
c Parameters:
cI  i        - the requested vector number
cI  nval     - number of correlated orbitals
cI  ncf      - number of contracted basis functions
cI  irecadrx - the records info for Kext operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - number of correlated pairs
cI  ndiskx   - the file unit number where Kext bins are kept.
cI  lbinx    - the length of Kext bin
cI  thresh   - the integral threshold
cI  byt8     - .true.: double precision bins, otherwise int4+int1
cI  cc       - .true.: the CCD contributions are calculated, otherwise CID or 
c              CISD
cO  iresult  - the result is deposited here
cI  nmo      - the number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      use newpara
      implicit none
      integer idimen,nval,ncf,iresult,ioverlap,nslv1,npairs
      logical vorb
c
      integer ires_tmp,i,isla_gath,npasses,ipasssize,ireminder
      integer mataddr,istart,istop,ipass,info,islv_tid,islgid
c
      call mmark
      call dynamic_mmark
      call matmark
c
      call matconn('fin_res','r',idimen,nval,iresult)
      call matconn('overl','q',ncf,ncf,ioverlap)
      call dynamic_matdef('res_tmp','r',idimen,nval)
      call matzero('res_tmp')
      ires_tmp=mataddr('res_tmp')
      call dynamic_matdef('sla_gath','r',idimen,nval)
      isla_gath=mataddr('sla_gath')
c Divide npairs into groups, preferably the same as number of slaves.
      npasses=nslv
c let one slave to work on matrix mult. no more than 20 seconds
      do i=1,100
        ipasssize=npairs/npasses
        if (dble(ipasssize)*dble(idimen)**3d0/(1d9).gt.2d1) then
          npasses=nslv*i
        else 
          exit
        endif
      enddo
      ireminder=mod(npairs,npasses)
      ipasssize=npairs/npasses
      istart=1
      istop=0
      do ipass=1,npasses
        istart=istop+1
        istop=istart+ipasssize-1
        if (ireminder.gt.0) then
          istop=istop+1
          ireminder=ireminder-1
        else if (ipasssize.eq.0) then
          exit
        endif
        if (istop.gt.npairs.or.istart.gt.npairs) STOP 'Tr_gen, error'
        call para_recv(islv_tid,islgid,TxCCSDReq4)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob4)
      enddo
      istart=0
      istop=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq4)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob4)
        call para_recv_real(bl(isla_gath),idimen*nval,islgid,TxCCSDRes4)
        call matadd('sla_gath','res_tmp')
      enddo
      if (vorb) then
        call matadd('res_tmp','fin_res')
      else
        call matmmul2('overl','res_tmp','fin_res','n','n','a')
      endif
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine Tr_gen_slav(idimen, nval,  ncf,  vorb, nslv1,
     *                       iMasID, mygid, mytid)
c This subroutine calculates the second term of (38) formula taken from
c J. Chem. Phys. 104, 6286, (1996)
c Parameters:
cI  i        - the requested vector number
cI  nval     - number of correlated orbitals
cI  ncf      - number of contracted basis functions
cI  irecadrx - the records info for Kext operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - number of correlated pairs
cI  ndiskx   - the file unit number where Kext bins are kept.
cI  lbinx    - the length of Kext bin
cI  thresh   - the integral threshold
cI  byt8     - .true.: double precision bins, otherwise int4+int1
cI  cc       - .true.: the CCD contributions are calculated, otherwise CID or 
c              CISD
cO  iresult  - the result is deposited here
cI  nmo      - the number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      use newpara
      implicit none
      integer idimen,nval,ncf,nslv1,iMasID,mygid,mytid
      logical vorb
c
      integer iiR,icoef_dif,ires_tmp,idoubles,ir_vec,ires_vec
      integer mataddr,info,istart,istop,i,k,ik,imastgid
c
      call mmark
      call dynamic_mmark
      call matmark
c
      call dynamic_matdef('coef_dif','q',idimen,idimen)
      icoef_dif=mataddr('coef_dif')
      call pointRR(iiR)
      call dynamic_matdef('res_tmp','r',idimen,nval)
      call matzero('res_tmp')
      ires_tmp=mataddr('res_tmp')
#if 1
      do
        call para_send(mytid,0,TxCCSDReq4)
        call para_recv_pack(imastgid,TxCCSDJob4)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop,1)
        if (istart.eq.0.and.istop.eq.0) then
          call para_send_real(bl(ires_tmp),idimen*nval,0,TxCCSDRes4)
          exit
        endif
        ik=0
        do i=1,nval
          do k=1,i
            ik=ik+1
            if (ik.lt.istart) cycle
            if (ik.gt.istop)  exit
            call CoefRead('tt',i,k,idoubles)
            call CCLoperator(bl(idoubles),idimen,bl(icoef_dif))
            ir_vec=iiR+(k-1)*idimen
            ires_vec=ires_tmp+(i-1)*idimen
            call matconn('RR','r',idimen,1,ir_vec)
            call matconn('result','r',idimen,1,ires_vec)
            call matmmul2('coef_dif','RR','result','n','n','a')
            call matdisc('result')
            call matdisc('RR')
            if (i.ne.k) then
              ir_vec=iiR+(i-1)*idimen
              ires_vec=ires_tmp+(k-1)*idimen
              call matconn('RR','r',idimen,1,ir_vec)
              call matconn('result','r',idimen,1,ires_vec)
              call matmmul2('coef_dif','RR','result','t','n','a')
              call matdisc('result')
              call matdisc('RR')
            endif
          enddo
        enddo
      enddo
#endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
c it is called but never used
      subroutine GenAmplSymm_Mas(npairs, ncf,  nval, nodisk,    nfirst,
     *                           nmo,    vorb, af,   mastslave, ndisks,
     *                           ndiska, scr,  list,isympairs,ij_unique,
     *                           nslv1)
c This subroutine constructs amplitudes needed for EEO generation. It takes
c doubles amplitudes and combines them with singles amplitudes, then it builds
c AO amplitudes and resolves each matrix into symmetric and antisymmetric part.
c If called on slave it sets up common block only.
c
c Input parameters:
c npairs     - number of correlated pairs
c ncf        - number of contracted functions, passed to common block below
c nval       - number of correlated orbitals
c nodisk     - if nodisk=.true. transformed amplitudes are kept in memory,
c              otherwise on disk
c nfirst     - first correlated orbital (orbitals are ordered according to 
c              SCF energy)
c nmo        - number of occupied SCF orbitals
c vorb       - if vorb=.true. all matrices are in transformed form. This has 
c              the following implication for this procedure that amplitudes 
c              need to be transformed from MO to AO basis
c af         - logical variable indicating whether we run single job with loca
c              I/O (.false.) or parallel job using Array Files.
c mastslave  - this subroutine is called on master and on slave. It performs
c              different tasks on slave than on master. This character variable
c              passes information whether we are on slave or on master.
c ndisks     - the disk unit number where the symmetrized amplitudes are written
c ndiska     - the disk unit number where antisymm. amplitudes are written
c Output parameters:
c scr        - here the screening matrix for EEO integarls is contructed. scr
c              contains, for given mu,nu max of T(ij)(mu,nu) for all ij.
c              To be more specific, first ncf*(ncf+1)/2 elements contains
c              max values for symmetrized T and the rest of scr contains max
c              values of antisymmetrized T.
c Common block /SymmAnti/:
c Output parameters:
c istartsymm - starting memory location of all symmetrized amplitudes, 
c              assigned and filled here if nodisk=.true.
c istartanti - as above for antisymmetrized
c incf       - number of contracted functions
c indisks    - disk unit number where symm ampl. are stored
c indiska    - as above for anti
c aff        - logical variable indicating whether we run single job with local
c              I/O (.false.) or parallel job using Array Files.
c Common block /SymmAnti/:
c Input parameters:
c cache      - if .true. it instructs transformed matrices I/O subsystem to 
c              cache amplitudes taken from AF on local disk.
c integerT   - if .true. it instructs this subroutine to convert double 
c              precision amplitudes into integer numbers (for space saving)
c              It turned out not to work, especially for diffuse basis sets 
c              because of large amount of big amplitudes. There is still
c              possible to store few large amplitudes in real*8 form, but 
c              experience shows that for good convergence we need integrals 
c              accurate to 12 decimal digits, the same applies PROBABLY to
c              amplitudes which have only 9 digits of accuracy in integer form
c
      use memory
      use newpara
      implicit none
c
      integer npairs,ncf,nval,nfirst,nmo,ndisks,ndiska,list(*)
      integer isympairs(*),ij_unique(*),nslv1
      logical nodisk,vorb,af
      character*(*) mastslave
      real*8 scr(*)
c
      integer istartsymm,istartanti,incf,indisks,indiska,npairsa
      logical aff
      common/SymmAnti/istartsymm,istartanti,incf,indisks,indiska,aff,
     *                npairsa
c
      character*256 scrfile,filenamea,filenames
      integer :: icountere = 0,i,k,info,idimen,len1,len,npasses,iscr
      integer ipasssize,ireminder,istart,istop,ipass,islv_tid,islgid
      real*8 xmax,xmax_sla
      integer :: icallno = 0
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
c
      xmax=0.0d0
      aff=af
      icountere=icountere+1
      npairsa=npairs
      icallno=icallno+1
      incf=ncf
      STOP 'Switched off'
c
      if (vorb) then
        idimen=ncf-nmo
      else
        idimen=ncf
      endif
c
      if (nodisk) then
        STOP 'nodisk not supported in symmanti'
      else   !!! * * * * * * * * * DISK * * * * * * * * *
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                 ! 80 is starting lenght
        filenames=scrfile(1:len)//'.ampl.'//'symm.'
        filenamea=scrfile(1:len)//'.ampl.'//'anti.'
        len1=len+6+5
        STOP 'Switched off'
        call fafCreatem(filenames(1:len1)//char(0),ndisks)
        call fafCreatem(filenamea(1:len1)//char(0),ndiska)
c
        call para_initsend
        call para_pack_int(ndisks,1)
        call para_pack_int(ndiska,1)
        call para_bcast_pack(TxCCSDInit)
c Divide npairs into groups, preferably the same as number of slaves.
      npasses=nslv
c let one slave to work on matrix mult. no more than 20 seconds
      do i=1,100
        ipasssize=npairs/npasses
        if (dble(ipasssize)*dble(idimen)**3d0/(1d9).gt.2d1) then
          npasses=nslv*i
        else
          exit
        endif
      enddo
      ireminder=mod(npairs,npasses)
      ipasssize=npairs/npasses
      istart=1
      istop=0
      do ipass=1,npasses
        istart=istop+1
        istop=istart+ipasssize-1
        if (ireminder.gt.0) then
          istop=istop+1
          ireminder=ireminder-1
        endif
        if (istop.gt.npairs.or.istart.gt.npairs) STOP 'CCAOnce, error'
        call para_recv(islv_tid,islgid,TxCCSDReq4)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob4)
      enddo
      istart=0
      istop=0
      call dynamic_getmem(ncf*ncf,iscr)
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq4)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop,1)
        call para_send_pack(islgid,TxCCSDJob4)
        call para_recv_real(xmax_sla,1,islgid,TxCCSDRes4)
        call para_recv_real(bl(iscr),ncf*ncf,islgid,TxCCSDRes4)
        if (dabs(xmax_sla).gt.xmax) xmax=dabs(xmax_sla)
        do k=0,ncf*ncf-1
          if (bl(iscr+k) .gt. scr(k+1)) scr(k+1)=bl(iscr+k)
        enddo
      enddo
      call dynamic_retmem(1)
      endif
      indisks=ndisks
      indiska=ndiska
      end
c==============================================================================
c called never used
      subroutine GenAmplSymm_Sla(npairs, ncf,  nval, nodisk,    nfirst,
     *                           nmo,    vorb, af,   mastslave, ndisks,
     *                           ndiska, scr,  list,isympairs,ij_unique,
     *                           nslv1,iMasID, mygid, mytid)
c This subroutine constructs amplitudes needed for EEO generation. It takes
c doubles amplitudes and combines them with singles amplitudes, then it builds
c AO amplitudes and resolves each matrix into symmetric and antisymmetric part.
c If called on slave it sets up common block only.
c
c Input parameters:
c npairs     - number of correlated pairs
c ncf        - number of contracted functions, passed to common block below
c nval       - number of correlated orbitals
c nodisk     - if nodisk=.true. transformed amplitudes are kept in memory,
c              otherwise on disk
c nfirst     - first correlated orbital (orbitals are ordered according to 
c              SCF energy)
c nmo        - number of occupied SCF orbitals
c vorb       - if vorb=.true. all matrices are in transformed form. This has 
c              the following implication for this procedure that amplitudes 
c              need to be transformed from MO to AO basis
c af         - logical variable indicating whether we run single job with loca
c              I/O (.false.) or parallel job using Array Files.
c mastslave  - this subroutine is called on master and on slave. It performs
c              different tasks on slave than on master. This character variable
c              passes information whether we are on slave or on master.
c ndisks     - the disk unit number where the symmetrized amplitudes are written
c ndiska     - the disk unit number where antisymm. amplitudes are written
c Output parameters:
c scr        - here the screening matrix for EEO integarls is contructed. scr
c              contains, for given mu,nu max of T(ij)(mu,nu) for all ij.
c              To be more specific, first ncf*(ncf+1)/2 elements contains
c              max values for symmetrized T and the rest of scr contains max
c              values of antisymmetrized T.
c Common block /SymmAnti/:
c Output parameters:
c istartsymm - starting memory location of all symmetrized amplitudes, 
c              assigned and filled here if nodisk=.true.
c istartanti - as above for antisymmetrized
c incf       - number of contracted functions
c indisks    - disk unit number where symm ampl. are stored
c indiska    - as above for anti
c aff        - logical variable indicating whether we run single job with local
c              I/O (.false.) or parallel job using Array Files.
c Common block /SymmAnti/:
c Input parameters:
c cache      - if .true. it instructs transformed matrices I/O subsystem to 
c              cache amplitudes taken from AF on local disk.
c integerT   - if .true. it instructs this subroutine to convert double 
c              precision amplitudes into integer numbers (for space saving)
c              It turned out not to work, especially for diffuse basis sets 
c              because of large amount of big amplitudes. There is still
c              possible to store few large amplitudes in real*8 form, but 
c              experience shows that for good convergence we need integrals 
c              accurate to 12 decimal digits, the same applies PROBABLY to
c              amplitudes which have only 9 digits of accuracy in integer form
c
      use memory
      use newpara
      implicit none
c
      integer npairs,ncf,nval,nfirst,nmo,ndisks,ndiska,list(*)
      integer isympairs(*),ij_unique(*),nslv1,iMasID, mygid, mytid
      logical nodisk,vorb,af
      character*(*) mastslave
      real*8 scr(*)
c
      integer istartsymm,istartanti,incf,indisks,indiska,npairsa
      logical aff
      common/SymmAnti/istartsymm,istartanti,incf,indisks,indiska,aff,
     *                npairsa
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
c
      integer :: icallno = 0,isingleMO,idimen,isinglvec,iworkhere,info
      integer irefsymm,irefanti,istart,istop,ij,i,j,ij_rec,icoeff
      integer imastgid
      integer :: icountere = 0,ntrisymm,ntrianti,mataddr
      real*8 xmax
c
      aff=af
      xmax=0.0d0
      icountere=icountere+1
c
      npairsa=npairs
      STOP 'Switched off'
      if (vorb) then
        idimen=ncf-nmo
      else
        idimen=ncf
      endif
c
      icallno=icallno+1
      incf=ncf
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      if (nodisk) then
        STOP 'nodisk not supported in symmanti'
      else
        if (vorb) then 
          call dynamic_matdef('worktran','r',ncf,idimen)
          call dynamic_matdef('singleMO','r',idimen,1)
          isingleMO=mataddr('singleMO')
        endif
        call dynamic_matdef('singlvec','r',ncf,1)
        isinglvec=mataddr('singlvec')
        call dynamic_matdef('workhere','q',ncf,ncf)
        iworkhere=mataddr('workhere')
c
c RECEIVE ndisks and ndiska from master
c
        call para_recv_bcastpack(TxCCSDInit)
        call para_unpack_int(ndisks,1)
        call para_unpack_int(ndiska,1)
c
        call dynamic_getmem(ntrisymm,irefsymm)
        call dynamic_getmem(ntrianti,irefanti)
      do
        call para_send(mytid,0,TxCCSDReq4)
        call para_recv_pack(imastgid,TxCCSDJob4)
        call para_unpack_int(istart,1)
        call para_unpack_int(istop,1)
        if (istart.eq.0.and.istop.eq.0) then
          call para_send_real(xmax,1,0,TxCCSDRes4)
          call para_send_real(scr,ncf*ncf,0,TxCCSDRes4)
          exit
        endif
        ij=0
        ij_rec=0
        do i=1,nval
          do j=1,i
            ij=ij+1
            if ((list(ij).eq.1.or.icountere.eq.1)
     *                        .and.isympairs(ij).ne.0) then
              continue
            else
              cycle
            endif
            ij_rec=ij_rec+1
            if (ij.lt.istart) cycle
            if (ij.gt.istop)  exit
            call CoefRead('tc',i,j,icoeff)
            if (vorb) then
            call matconn('coeff','q',idimen,idimen,icoeff)
            call matmmult('genvirt','coeff','worktran')
            call matmmul2('worktran','genvirt','workhere','n','t','n')
            call matdisc('coeff')
            else
            call matconn('coeff','q',ncf,ncf,icoeff)
            call matcopy('coeff','workhere')
            call matdisc('coeff')
            endif
c Singles part of EEO:
            if (list(ij).eq.1) then
            if (vorb) then
            call readsingles(j,bl(isingleMO))
            call matmmult('genvirt','singleMO','singlvec')
            else
            call readsingles(j,bl(isinglvec))
            endif
            call matsub('singlocc','occu',i,i)
            call matmmul2('singlocc','singlvec','workhere','n','t','a')
            call matdisc('singlocc')
c
            if (vorb) then
            call readsingles(i,bl(isingleMO))
            call matmmult('genvirt','singleMO','singlvec')
            else
            call readsingles(i,bl(isinglvec))
            endif
            call matsub('singlocc','occu',j,j)
            call matmmul2('singlvec','singlocc','workhere','n','t','a')
            call matdisc('singlocc')
            endif
c Singles part of EEO end
            call SymAntiBuild(bl(iworkhere),bl(irefsymm),
     *                        bl(irefanti),ncf,scr)
            call WriterCoefLin(ndisks,ij_rec,bl(irefsymm),ntrisymm,af,
     *                         xmax)
            call WriterCoefLin(ndiska,ij_rec,bl(irefanti),ntrianti,af,
     *                         xmax)
          enddo
        enddo
      enddo
        call dynamic_retmem(2)
      endif
      call dynamic_matrem('workhere')
      call dynamic_matrem('singlvec')
      if (vorb) call dynamic_matrem('singleMO')
      if (vorb) call dynamic_matrem('worktran')
      indisks=ndisks
      indiska=ndiska
      end
c==============================================================================
      subroutine new_tri_mas(nval,idimen,ndisk_a,ndisk_ie,ndisk_ii,
     *                       ndisk_ix,energy,eps,etimes,n,ccsd,
     *                       qcisd,af,nslv1,ext3_ratio,amp_ratio,npass,
     *                       isize,iprnt,energys,energyd,iabcstart)
      use memory
      use newpara
      implicit none
c
      integer nval,idimen,ndisk_a,ndisk_ie,ndisk_ii,ndisk_ix,n,nslv1
      integer iabcstart
      real*8 energy,eps(*),etimes(n),ext3_ratio,amp_ratio,energys
      real*8 energyd
      logical ccsd,qcisd,af
      integer npass,isize,iprnt
c
      integer iWabc,ia,ib,ic,i,irec_amp,n_ampl,iamplit,ip,ip_dim
      integer irec_int,n_int,i3int_int,irec_ext,n_ext,i3ext_int,isingles
      integer iap,ibp,icp,ia1,ia2,ib1,ib2,ic1,ic2,hb,hc
      integer irec,info,iK,irec_K,n_K,nvirpairs
      real*8 en,t0,t1,x0,x1,etimes1(n),energy_sl,percent,energy_sl1
      integer, parameter :: MAX_RES=10000
      real*8 energy_sl2,results(2,MAX_RES),pres
      integer :: itop
      logical reset1,reset2,reset3,reset4
      integer ip_array,fragmentation,iampl_store,i3int_store,i3ext_store
      integer iK_store,ialloc_amount,mem,ii,islv_tid,idone,idiffe
      integer ngranules,itrysize,ipass,irem,tot_work,ijk_done,iet,iet0
      integer islgid,iapr,ibpr,icpr,iabcr,iabc
      integer, dimension(8) :: dt,dt0
c
      itop=1
c
      call matmark()
      call dynamic_mmark()
      do i=1,MAX_RES
        results(1,i)=0d0
        results(2,i)=0d0
      enddo
      results(1,1)=iabcstart-1
      if (results(1,1).lt.0d0) results(1,1)=0d0
      results(2,1)=0d0
c
c etimes: 1 the time spend in W zeroing and sorting
c         2 W building, total time
c         3 in W: multiplication - over virt
c         4 in W: multiplication - over occ
c         5 reading
c         6 reading locate amplitudes
c         7 reading locate 3int integrals
c         8 reading locate 3ext integrals
      energy=0d0
      energys=0d0
      energyd=0d0
      do i=1,n
        etimes(i)=0d0
      enddo
      nvirpairs=idimen*(idimen+1)/2
c
      call dynamic_getmem(nval*nval*nval,iWabc)
c MEMORY MANAGEMENT:
      if (iprnt.gt.2)
     *write(6,*) 'Triples: calculated chunk size = ',isize,',  tail: ',
     *         mod(idimen,isize)
      call flush(6)
      call para_initsend
      call para_pack_int(isize,1)
      call para_pack_int(ndisk_a,1)
      call para_pack_int(ndisk_ie,1)
      call para_pack_int(ndisk_ii,1)
      call para_pack_int(ndisk_ix,1)
      call para_pack_int(n,1)
      call para_pack_real(ext3_ratio,1)
      call para_pack_real(amp_ratio,1)
      call para_bcast_pack(TxCCSDInit)
c END OF MEMORY MANAGEMENT
      tot_work=idimen*(idimen+1)*(idimen+2)/6-idimen
      ijk_done=0
      call date_and_time(values = dt0)
      call array_to_scalar(dt0,iet0)
      iabc=0
      do iap=1,npass
        ia1=(iap-1)*isize+1
        ia2=(iap  )*isize
        if (ia1.gt.idimen) STOP 'Error 054h32jvhsjdw;a'
        if (ia2.gt.idimen) ia2=idimen
        do ibp=1,iap
          ib1=(ibp-1)*isize+1
          ib2=(ibp  )*isize
          if (ib1.gt.idimen) STOP 'Error 054h32jvhsjdw;a'
          if (ib2.gt.idimen) ib2=idimen
          do icp=1,ibp
            ic1=(icp-1)*isize+1
            ic2=(icp  )*isize
            if (ic1.gt.idimen) STOP 'Error 054h32jvhsjdw;a'
            if (ic2.gt.idimen) ic2=idimen
              iabc=iabc+1
              if (iabc<iabcstart) cycle
              call para_recv_pack(islgid,TxTriplesReq)
              call para_unpack_int(islv_tid,1)
              call para_unpack_int(idone,1)
              if (idone.gt.0) then
                call para_unpack_int(iapr,1)
                call para_unpack_int(ibpr,1)
                call para_unpack_int(icpr,1)
                call para_unpack_real(pres,1)
                iabcr=((iapr-1)*iapr*(iapr+1))/6+((ibpr-1)*ibpr)/2+icpr
                itop=itop+1
                if (itop.ge.MAX_RES) STOP 'Condition unexpected'
                results(1,itop)=iabcr
                results(2,itop)=pres
                if (mod(iabc,20).eq.0) then
                  call sort_and_compress(results,itop,iprnt)
                endif
              endif
              ijk_done=ijk_done+idone
              call date_and_time(values = dt)
              call array_to_scalar(dt,iet)
              idiffe=iet-iet0
              percent=dble(ijk_done)/tot_work*1d2
              if (ijk_done.gt.0) then
                idiffe=dble(idiffe)/percent*1d2
                iet=iet0+idiffe
              endif
              if (iprnt.gt.2) then
              write(6,
     *        '(A,F6.2,A1,5X,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)',
     *        ADVANCE='NO')
     *              'Progress: ',
     *              percent,'%',dt(1),'-',
     *              dt(2),'-',dt(3),', ',dt(5),':',dt(6),':',dt(7)
              if (ijk_done.gt.0) then
                call scalar_to_array(iet,dt)
                write(6,'(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')
     *                 ",  Est. finish ",dt(1),'-',dt(2),'-',dt(3),', ',
     *                                   dt(5),':',dt(6),':',dt(7)
              else
                write(6,*)
              endif
              call flush(6)
              endif
              call para_initsend
              call para_pack_int(ia1,1)
              call para_pack_int(ia2,1)
              call para_pack_int(ib1,1)
              call para_pack_int(ib2,1)
              call para_pack_int(ic1,1)
              call para_pack_int(ic2,1)
              call para_pack_int(iap,1)
              call para_pack_int(ibp,1)
              call para_pack_int(icp,1)
              call para_send_pack(islgid,TxTriplesJob)
              call flush(6)
          enddo
        enddo
      enddo
      ia1=0
      do i=1,nslv
        call para_recv_pack(islgid,TxTriplesReq)
        call para_unpack_int(islv_tid,1)
        call para_unpack_int(idone,1)
              if (idone.gt.0) then
                call para_unpack_int(iapr,1)
                call para_unpack_int(ibpr,1)
                call para_unpack_int(icpr,1)
                call para_unpack_real(pres,1)
                iabcr=((iapr-1)*iapr*(iapr+1))/6+((ibpr-1)*ibpr)/2+icpr
                itop=itop+1
                if (itop.ge.MAX_RES) STOP 'Condition unexpected'
                results(1,itop)=iabcr
                results(2,itop)=pres
                call sort_and_compress(results,itop,iprnt)
              endif
        ijk_done=ijk_done+idone
        call date_and_time(values = dt)
        if (iprnt.gt.2) then
        write(6,
     *  '(A,F6.2,A1,10X,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')
     *        'Progress: ',
     *        dble(ijk_done)/tot_work*1d2,'%',dt(1),'-',
     *        dt(2),'-',dt(3),', ',dt(5),':',dt(6),':',dt(7)
        call flush(6)
        endif
        call para_send(ia1,islgid,TxTriplesJob)
        call para_recv_real(energy_sl,1,islgid,TxTriplesRes)
        call para_recv_real(energy_sl1,1,islgid,TxTriplesRes)
        call para_recv_real(energy_sl2,1,islgid,TxTriplesRes)
        call para_recv_real(etimes1,n,islgid,TxTriplesRes)
        energy=energy+energy_sl
        energys=energys+energy_sl1
        energyd=energyd+energy_sl2
        do ii=1,n
          etimes(ii)=etimes(ii)+etimes1(ii)
        enddo
      enddo
c
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine new_tri_sla(nval,idimen,eps,ccsd,qcisd,af,iMasID,
     *                     mygid, mytid,otable,vtable,vrevtable,m_table,
     *                     ijtable,nrep)
      use memory
      use newpara
      implicit none
c
      integer nval,idimen,iMasID,mygid,mytid
      real*8 eps(*)
      logical ccsd,qcisd,af
      integer otable(*),vtable(*),vrevtable(*),m_table(*),ijtable(*)
      integer nrep
c
      integer :: ndisk_a,ndisk_ie,ndisk_ii,ndisk_ix,n
      real*8 etimes(100),energy,enes,ened,en1
      integer iWabc,ia,ib,ic,i,irec_amp,n_ampl,iamplit,ip,ip_dim
      integer irec_int,n_int,i3int_int,irec_ext,n_ext,i3ext_int,isingles
      integer iap,ibp,icp,ia1,ia2,ib1,ib2,ic1,ic2,isize,npass,hb,hc
      integer iKab,iKac,iKbc,irec,info,iK,irec_K,n_K,nvirpairs
      real*8 en,t0,t1,x0,x1,ext3_ratio,amp_ratio,energys,energyd
      logical reset1,reset2,reset3,reset4
      integer ip_array,fragmentation,iampl_store,i3int_store,i3ext_store
      integer iK_store,ialloc_amount,mem,idone,ireal_size,ir_amp_size
      integer imastgid
c
      call matmark()
      call dynamic_mmark()
c
c etimes: 1 the time spend in W zeroing and sorting
c         2 W building, total time
c         3 in W: multiplication - over virt
c         4 in W: multiplication - over occ
c         5 reading
c         6 reading locate amplitudes
c         7 reading locate 3int integrals
c         8 reading locate 3ext integrals
      energy=0d0
      energys=0d0
      energyd=0d0
      nvirpairs=idimen*(idimen+1)/2
c
      call dynamic_getmem(nval*nval*nval,iWabc)
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(isize,1)
      call para_unpack_int(ndisk_a,1)
      call para_unpack_int(ndisk_ie,1)
      call para_unpack_int(ndisk_ii,1)
      call para_unpack_int(ndisk_ix,1)
      call para_unpack_int(n,1)
      call para_unpack_real(ext3_ratio,1)
      call para_unpack_real(amp_ratio,1)
      do i=1,n
        etimes(i)=0d0
      enddo
c
      ip_dim=6*isize*isize  ! is it woth to reduce to integer*4?
      call dynamic_getmem(ip_dim,ip)
      npass=idimen/isize
      if (mod(idimen,isize).ne.0) npass=npass+1
c locate_ampl setup:
      call dynamic_getmem(idimen,irec_amp)
      reset1=.true.
      n_ampl=isize*3
      ir_amp_size=int(amp_ratio*nval*nval*idimen+1)
      if (amp_ratio.gt.0.999999999999d0) ir_amp_size=nval*nval*idimen
      call dynamic_getmem(ir_amp_size*n_ampl,iamplit)
c locate_3int setup:
      call dynamic_getmem(idimen,irec_int)
      reset2=.true.
      n_int=isize*3
      call dynamic_getmem(nval*nval*nval*n_int,i3int_int)
c locate_3ext setup:
      call dynamic_getmem(idimen*idimen,irec_ext)
      reset3=.true.
      n_ext=3*isize*isize
      ireal_size=int(ext3_ratio*nval*idimen+1)
      if (ext3_ratio.gt.0.999999999999d0) ireal_size=nval*idimen
      call dynamic_getmem(ireal_size*n_ext*2,i3ext_int)
c locate_K setup:
      call dynamic_getmem(nvirpairs,irec_K)
      reset4=.true.
      n_K=3*isize*isize
      call dynamic_getmem(nval*nval*n_K,iK)
c
      call pointersingles(isingles)
c
      idone=0
      do
c       reset3=.true.
        reset4=.true.
        call para_initsend
        call para_pack_int(mytid,1)
        call para_pack_int(idone,1)
        if (idone.gt.0) then
          call para_pack_int(iap,1)
          call para_pack_int(ibp,1)
          call para_pack_int(icp,1)
          call para_pack_real(en1,1)
        endif
        call para_send_pack(0,TxTriplesReq)
        call para_recv_pack(imastgid,TxTriplesJob)
        call para_unpack_int(ia1,1)
        if (ia1.eq.0) then
          call para_send_real(energy,1,0,TxTriplesRes)
          call para_send_real(energys,1,0,TxTriplesRes)
          call para_send_real(energyd,1,0,TxTriplesRes)
          call para_send_real(etimes,n,0,TxTriplesRes)
          exit
        endif
        call para_unpack_int(ia2,1)
        call para_unpack_int(ib1,1)
        call para_unpack_int(ib2,1)
        call para_unpack_int(ic1,1)
        call para_unpack_int(ic2,1)
        call para_unpack_int(iap,1)
        call para_unpack_int(ibp,1)
        call para_unpack_int(icp,1)
            call elapsec(x0)
            call elapsec(t0)
            call locate_ampl(ia1,   ia2,    ib1,         ib2,   ic1,
     *                       ic2,   iamplit,bl(irec_amp),reset1,n_ampl,
     *                       idimen,nval,   ndisk_a,     af,    m_table,
     *                       ijtable,vtable,vrevtable, nrep,ir_amp_size)
            call elapsec(x1)
            etimes(6)=etimes(6)+x1-x0
            call elapsec(x0)
            call locate_3int(ia1,   ia2,      ib1,         ib2,   ic1,
     *                       ic2,   i3int_int,bl(irec_int),reset2,n_int,
     *                       idimen,nval,     ndisk_ii,    af)
            call elapsec(x1)
            etimes(7)=etimes(7)+x1-x0
            call elapsec(x0)
            call locate_3ext(ia1,   ia2,      ib1,         ib2,   ic1,
     *                       ic2,   i3ext_int,bl(irec_ext),reset3,n_ext,
     *                       idimen,nval,     ndisk_ie,    bl(ip),af,
     *                       nrep,vtable,m_table,otable,vrevtable,
     *                       ireal_size)
            call elapsec(x1)
            etimes(8)=etimes(8)+x1-x0
            call elapsec(x0)
            if (qcisd.or.ccsd) then
            call locate_K   (ia1,   ia2,      ib1,         ib2,   ic1,
     *                       ic2,   iK,       bl(irec_K),  reset4,n_K,
     *                       idimen,nval,     ndisk_ix,    bl(ip),af,
     *                       isize, npass,    iap,         ibp,   icp)
            endif
            call elapsec(x1)
            etimes(10)=etimes(10)+x1-x0
            call elapsec(t1)
            etimes(5)=etimes(5)+t1-t0
c
c
      idone=0
      en1=0d0
      do ia=ia1,ia2
        hb=min(ia,ib2)
        do ib=ib1,hb
          hc=min(ib,ic2)
          do ic=ic1,hc
            if (ia.eq.ib.and.ib.eq.ic) cycle
            idone=idone+1
c           print *,ia,ib,ic
            if (qcisd.or.ccsd) then
c             irec=ia*(ia-1)/2+ib
c             if (af) then
c             call fafread(ndisk_ix,bl(iKab),8,nval*nval,1,irec,info)
c             else
c             call general_read(ndisk_ix,irec,bl(iKab),nval*nval*8)
c             endif
c             irec=ia*(ia-1)/2+ic
c             if (af) then
c             call fafread(ndisk_ix,bl(iKac),8,nval*nval,1,irec,info)
c             else
c             call general_read(ndisk_ix,irec,bl(iKac),nval*nval*8)
c             endif
c             irec=ib*(ib-1)/2+ic
c             if (af) then
c             call fafread(ndisk_ix,bl(iKbc),8,nval*nval,1,irec,info)
c             else
c             call general_read(ndisk_ix,irec,bl(iKbc),nval*nval*8)
c             endif
              call return_addresses(ia,ib,ic,bl(irec_K),iKab,iKac,iKbc)
            endif
            call elapsec(t0)
            call zeroit(bl(iWabc),nval*nval*nval)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ia,      ib,       ic,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,2,3)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ia,      ic,       ib,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,1,2)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ic,      ia,       ib,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,2,3)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ic,      ib,       ia,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,1,2)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ib,      ic,       ia,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,2,3)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ib,      ia,       ic,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,1,2)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Energy_fromW(ia,  ib, ic,          nval,     bl(iWabc),
     *                        eps, en, bl(isingles),ccsd,     qcisd,
     *                        bl(iKab),bl(iKac),bl(iKbc),idimen,enes,
     *                        ened)
            call elapsec(t1)
            etimes(9)=etimes(9)+t1-t0
            energy=energy+en
            energys=energys+enes
            energyd=energyd+ened
            en1=en1+en
          enddo
        enddo
      enddo
      enddo ! infinite loop
c
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
c Integrals - 3ext J(d,i)(a,b) type -> K(b,d)(a,i) type
c===============================================================================
      subroutine sort_3ext_mas(nval,idimen,ndisk,irecadrtt,npairs,
     *                     ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,
     *                     af,igid,nslv1,irecord)
      use memory
      use kinds
      use newpara
      implicit none
c
      integer nval,idimen,ndisk,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nmo,igid,irecord(*),nslv1
      logical byt8,vorb
      real*8 thresh
      logical af
c
c
      integer memory1,ifill_level,lbin,ibin,iindex,len,len1
      integer ndisk_sort,ixmat,no_bins,info,nrec,isize,j,ij
      integer ipasssize,npass,ipass,istart,istop,islv_tid,i,iaddr,islgid
      character*256 scrfile,filename
      logical fileopen
c
      call matmark()
      call dynamic_mmark()
c
      call dynamic_getmem(idimen*idimen,ifill_level)
      call izeroit(bl(ifill_level),idimen*idimen*intsize)
c
c I need memory for bins, and indexes (2 indexes)
c All the bins will occuppy the amount of memory:
c idimen*idimen*lbin*(8+2*2)
c
      call dynamic_show_free(memory1)
      memory1=memory1-300000-idimen*idimen ! possible fragmentation
      memory1=memory1-idimen*idimen ! one memory reservation in ext3_sort
      lbin=8*memory1/(idimen*idimen*12)
      if (lbin.gt.nval*idimen) lbin=nval*idimen
      no_bins=nval*idimen/lbin
      if (mod(nval*idimen,lbin).ne.0) no_bins=no_bins+1
c
      call dynamic_mmark()
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
      filename=scrfile(1:len)//'.sorttmp'
      len1=len+8
      ndisk_sort=50
      STOP 'Switched off'
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk_sort)
      else
        inquire(ndisk_sort,opened=fileopen)
        if (fileopen) STOP 'Unit 50 opened'
        OPEN(UNIT=ndisk_sort,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lbin*12+8/intsize)
      endif
      filename=scrfile(1:len)//'.Kbc'
      len1=len+4
      ndisk=52
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 52 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*idimen*8)
      endif
c
c send lbin,ndisk_sort
      call para_initsend
      call para_pack_int(lbin,1)
      call para_pack_int(ndisk_sort,1)
      call para_pack_int(ndisk,1)
      call para_bcast_pack(TxCCSDInit)
c
      do isize=idimen*nval,1,-1
        npass=idimen*nval/isize
        if (mod(idimen*nval,isize).ne.0) npass=npass+1
        if (npass.ge.nslv) then
          ipasssize=isize
          exit
        endif
      enddo
      call optimal_isize(ipasssize,idimen*nval)
      npass=idimen*nval/ipasssize
      if (mod(idimen*nval,ipasssize).ne.0) npass=npass+1
      do ipass=1,npass
        istart=(ipass-1)*ipasssize+1
        istop =(ipass  )*ipasssize
        if (istart.gt.idimen*nval) STOP 'Error 392jfallaer'
        if (istop .gt.idimen*nval) istop=idimen*nval
        call para_recv(islv_tid,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop ,1)
        call para_send_pack(islgid,TxCCSDJob3)
      enddo
      istart=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(istart,1)
        call para_send_pack(islgid,TxCCSDJob3)
        iaddr=(i-1)*idimen*idimen+1
      call para_recv_int(irecord(iaddr),idimen*idimen,islgid,TxCCSDRes3)
      enddo
ccs      print *,'Entering sleep state for bin flush'
ccs      call flush(6)
ccs      call sleep(900)
ccs      print *,'Exiting sleep state'
      call flush(6)
c
      call para_initsend
      call para_pack_int(irecord,idimen*idimen*nslv)
      call para_bcast_pack(TxCCSDInit)
      do isize=idimen*idimen,1,-1
        npass=idimen*idimen/isize
        if (mod(idimen*idimen,isize).ne.0) npass=npass+1
        if (npass.ge.nslv) then
          ipasssize=isize
          exit
        endif
      enddo
      call optimal_isize(ipasssize,idimen*idimen)
      npass=idimen*idimen/ipasssize
      if (mod(idimen*idimen,ipasssize).ne.0) npass=npass+1
      print *,'Optimal pas. size: ',ipasssize,', npass: ',npass,
     *         ', tail: ',mod(idimen*idimen,ipasssize)
      call flush(6)
      do ipass=1,npass
        istart=(ipass-1)*ipasssize+1
        istop= (ipass  )*ipasssize
        if (istart.gt.idimen*idimen) STOP 'Error 392jfallaer'
        if (istop .gt.idimen*idimen) istop=idimen*idimen
        print *,'istart,istop: ',istart,istop
        call flush(6)
        call para_recv(islv_tid,islgid,TxCCSDReq2)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop ,1)
        call para_send_pack(TxCCSDJob2)
      enddo
      istart=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq2)
        call para_initsend
        call para_pack_int(istart,1)
        call para_send_pack(TxCCSDJob2)
      enddo
      if (af) then
        call fafClosem(ndisk_sort,0,info)
      else
        close(ndisk_sort,status='delete')
      endif
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine sort_3ext_sla(nval,idimen,irecadrtt,npairs,
     *                     ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,
     *                     af,iMasID,mygid,mytid,nslv1)
      use memory
      use kinds
      use newpara
      implicit none
c
      integer nval,idimen,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nmo,iMasID,mygid,mytid,nslv1
      logical byt8,vorb
      real*8 thresh
      logical af
c
c
      integer memory1,irecord,ifill_level,lbin,ibin,iindex,len,len1
      integer ndisk_sort,ixmat,no_bins,info,nrec,ndisk
      integer ipasssize,npass,ipass,istart,istop,imastgid
      character*256 scrfile,filename
      logical fileopen
c
      call matmark()
      call dynamic_mmark()
c
      call dynamic_getmem(idimen*idimen*nslv,irecord)
      call izeroit(bl(irecord),idimen*idimen*nslv*intsize)
      call dynamic_getmem(idimen*idimen,ifill_level)
      call izeroit(bl(ifill_level),idimen*idimen*intsize)
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(lbin,1)
      call para_unpack_int(ndisk_sort,1)
      call para_unpack_int(ndisk,1)
c
      call dynamic_mmark()
      call dynamic_getmem(idimen*idimen*lbin,ibin)
      call dynamic_getmem(idimen*idimen*lbin/2+1,iindex)
      call izeroit(bl(iindex),(idimen*idimen*lbin/2+1)*intsize)
c
      do
        call para_send(mytid,0,TxCCSDReq3)
        call para_recv_pack(imastgid,TxCCSDJob3)
        call para_unpack_int(istart,1)
        if (istart.eq.0) then
      call leftovers(idimen,   bl(ifill_level),lbin,    bl(iindex),nrec,
     *             bl(irecord),ndisk_sort,     bl(ibin),af,       mygid,
     *             idimen*idimen,2)
          call para_send_int(bl(irecord),idimen*idimen,0,TxCCSDRes3)
          exit
        endif
        call para_unpack_int(istop ,1)
      call ext3_sort(bl(ibin),bl(iindex),nval,bl(ifill_level),idimen,
     *               bl(irecord),ndisk_sort,lbin,irecadrtt,npairs,
     *               ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,af,mygid,
     *               nrec,istart,istop)
      enddo
c
      call dynamic_retmark()
c
      call dynamic_mmark()
      call dynamic_getmem(lbin,ibin)
      call dynamic_getmem(lbin/2+1,iindex)
      call dynamic_getmem(nval*idimen,ixmat)
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(bl(irecord),idimen*idimen*nslv)
c
      do
        call para_send(mytid,0,TxCCSDReq2)
        call para_recv_pack(imastgid,TxCCSDJob2)
        call para_unpack_int(istart,1)
        if (istart.eq.0) exit
        call para_unpack_int(istop, 1)
        call coll_ext3_bins(lbin,bl(ibin),bl(iindex),ndisk_sort,
     *                      ndisk,bl(irecord),idimen,bl(ixmat),
     *                      nval,af,istart,istop,nslv)
      enddo
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
c Integrals - 3ext J(d,i)(a,b) type -> K(b,d)(a,i) type, direct sort
c===============================================================================
      subroutine master_3ext(nval,idimen,ndisk,irecadrtt,npairs,ndisktt,
     *                       ncf,lbintt,thresh,byt8,nmo,vorb,af,igid,
     *                       nslv1,c_ratio)
      use memory
      use kinds
      use newpara
      implicit none
      integer nval,idimen,ndisk,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nmo,igid,nslv1
      real*8 thresh,c_ratio
      logical byt8,vorb,af
c
      character*256 scrfile,filename
      integer len,len1,ibigmat,ixmat,istart,istop,ipasssize,npass
      integer isize,ia,info,islv_tid,i,islgid
      logical fileopen
      real*8 c_ratio1
c
c
      call matmark()
      call dynamic_mmark()
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len)
      filename=scrfile(1:len)//'.Kbc'
      len1=len+4
      ndisk=52
      c_ratio=0d0
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk,nval*idimen*8)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 52 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*idimen*8)
      endif
c
      do isize=idimen,1,-1
        npass=idimen/isize
        if (mod(idimen,isize).ne.0) npass=npass+1
        if (npass.ge.nslv) then
          ipasssize=isize
          exit
        endif
      enddo
      call optimal_isize(ipasssize,idimen)
      npass=idimen/ipasssize
      if (mod(idimen,ipasssize).ne.0) npass=npass+1
      call para_bcast(ndisk,TxCCSDInit)
c
      do ia=1,npass
        istart=(ia-1)*ipasssize+1
        istop= (ia  )*ipasssize
        if (istart.gt.idimen) STOP 'Error 392jfallaer'
        if (istop .gt.idimen) istop=idimen
        call para_recv(islv_tid,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop ,1)
        call para_send_pack(islgid,TxCCSDJob3)
      enddo
      istart=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(istart,1)
        call para_send_pack(islgid,TxCCSDJob3)
        call para_recv_real(c_ratio1,1,islgid,TxCCSDRes3)
        if (c_ratio1.gt.c_ratio) c_ratio=c_ratio1
      enddo
      call para_bcast_float(c_ratio,1,TxCoefInit)
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine slave_3ext(nval,idimen,irecadrtt,npairs,ndisktt,
     *                      ncf,lbintt,thresh,byt8,nmo,vorb,af,
     *                      iMasID,mygid,mytid,nrep,vtable,m_table,
     *                      vrevtable,otable,c_ratio)
      use memory
      use kinds
      use newpara
      implicit none
      integer nval,idimen,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nmo,iMasID,mygid,mytid
      integer nrep,vtable(*),m_table(*),vrevtable(*),otable(*)
      real*8 thresh,c_ratio
      logical byt8,vorb,af
c
      character*256 scrfile,filename
      integer len,len1,ibigmat,ixmat,istart,istop,ipasssize,npass
      integer isize,ia,ndisk,info,imastgid
      logical fileopen
c
c
      c_ratio=0d0
      call matmark()
      call dynamic_mmark()
      call fafCreates(ndisk)
c
      call para_recv_bcast(ndisk,TxCCSDInit)
c
      call dynamic_getmem(idimen*idimen*nval,ibigmat)
      call dynamic_getmem(idimen*idimen,ixmat)
      do
        call para_send(mytid,0,TxCCSDReq3)
        call para_recv_pack(imastgid,TxCCSDJob3)
        call para_unpack_int(istart,1)
        if (istart.eq.0) then 
          call para_send_real(c_ratio,1,0,TxCCSDRes3)
          exit
        endif
        call para_unpack_int(istop,1)
        call direct_ext3_sort(idimen,nval,bl(ixmat),ndisk,irecadrtt,
     *                        npairs,ndisktt,ncf,lbintt,thresh,
     *                        byt8,nmo,vorb,bl(ibigmat),istart,
     *                        istop,af,nrep,vtable,m_table,vrevtable,
     *                        otable,c_ratio)
      enddo
      call para_recv_bcastfloat(c_ratio,1,TxCoefInit)
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine make_3ext_pairs_mas(idimen,nval,af,ndisk_s,ndisk_t,
     *                               nslv1)
      use memory
      use newpara
      implicit none
      integer idimen,nval,ndisk_s,ndisk_t,nslv1
      logical af
c
      character*256 scrfile,filename
      integer len,len1,ibase,iab,ia,ib,irec,iaddr,info,nvirpairs
      logical fileopen
      integer isize,npass,islv_tid,istart,istop,i,ipasssize,islgid
c
      nvirpairs=idimen*(idimen+1)/2
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len)
      filename=scrfile(1:len)//'.Kbc.1.'
      len1=len+7
      ndisk_t=55
      if (af) then
       call fafCreatem(filename(1:len1)//char(0),ndisk_t,nval*idimen*16)
      else
        inquire(ndisk_t,opened=fileopen)
        if (fileopen) STOP 'Unit 55 opened'
        OPEN(UNIT=ndisk_t,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*idimen*16)
      endif
c
      do isize=nvirpairs,1,-1
        npass=nvirpairs/isize
        if (mod(nvirpairs,isize).ne.0) npass=npass+1
        if (npass.ge.nslv) then
          ipasssize=isize
          exit
        endif
      enddo
      call optimal_isize(ipasssize,nvirpairs)
      npass=nvirpairs/ipasssize
      if (mod(nvirpairs,ipasssize).ne.0) npass=npass+1
c
      call para_initsend
      call para_pack_int(ndisk_s,1)
      call para_pack_int(ndisk_t,1)
      call para_bcast_pack(TxCCSDInit)
c
      do ia=1,npass
        istart=(ia-1)*ipasssize+1
        istop= (ia  )*ipasssize
        if (istart.gt.nvirpairs) STOP 'Error 392jfallaer'
        if (istop .gt.nvirpairs) istop=nvirpairs
        call para_recv(islv_tid,islgid,TxCCSDReq1)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(istop ,1)
        call para_send_pack(islgid,TxCCSDJob1)
      enddo
      istart=0
      do i=1,nslv
        call para_recv(islv_tid,islgid,TxCCSDReq1)
        call para_initsend
        call para_pack_int(istart,1)
        call para_send_pack(islgid,TxCCSDJob1)
      enddo
      if (af) then
        call fafClosem(ndisk_s,0,info)
      else
        close(ndisk_s,status='delete')
      endif
      end
c===============================================================================
      subroutine make_3ext_pairs_sla(idimen,nval,af,nrep,vtable,
     *                               m_table,otable,vrevtable,iMasID,
     *                               mygid,mytid)
      use memory
      use newpara
      implicit none
      integer idimen,nval,iMasID,mygid,mytid,vtable(*),m_table(*)
      integer otable(*),vrevtable(*),nrep
      logical af
c
      integer ibase,iab,ia,ib,irec,iaddr,info,istart,istop
      integer ndisk_s,ndisk_t,isize,imastgid
c
      call dynamic_getmem(2*idimen*nval,ibase)
c
      call fafCreates(ndisk_t)
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ndisk_s,1)
      call para_unpack_int(ndisk_t,1)
c
      do
        call para_send(mytid,0,TxCCSDReq1)
        call para_recv_pack(imastgid,TxCCSDJob1)
        call para_unpack_int(istart,1)
        if (istart.eq.0) exit
        call para_unpack_int(istop,1)
        iab=0
      x:do ia=1,idimen
          do ib=1,ia
            iab=iab+1
            if (iab.lt.istart) cycle
            if (iab.gt.istop ) exit x
            irec=(ia-1)*idimen+ib
            iaddr=ibase
            if (af) then
              call ext3_sizer(nrep,idimen,nval,vtable,m_table,
     *                        otable,vrevtable,ia,ib,isize)
              call fafread(ndisk_s,bl(iaddr),8,isize,1,irec,info)
          if (info.lt.8*isize) then
              ERROR
              call flush(6)
          endif
              iaddr=ibase+isize
              irec=(ib-1)*idimen+ia
              call fafread(ndisk_s,bl(iaddr),8,isize,1,irec,info)
          if (info.lt.8*isize) then
              ERROR
              call flush(6)
          endif
             call fafwrite(ndisk_t,bl(ibase),8,isize*2,1,iab,info)
            else
              call general_read(ndisk_s,irec,bl(iaddr),idimen*nval*8)
              iaddr=ibase+idimen*nval
              irec=(ib-1)*idimen+ia
              call general_read(ndisk_s,irec,bl(iaddr),idimen*nval*8)
              call general_write(ndisk_t,iab,bl(ibase),idimen*nval*16)
            endif
          enddo
        enddo x
      enddo
      call dynamic_retmem(1)
      call fafcloses
      end
c===============================================================================
