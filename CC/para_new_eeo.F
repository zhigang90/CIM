#include "maxi.h"
C==============================================================================
      subroutine MAS_EEO_INT(ncs,    ncf,     inx,   nval,  nmo,
     *                       nfirst, nlast,   thresh,
     *                       vorb,   ndisktr, npairs, nslv1,iprnt,
     *                       nbf)
      use memory
      use kinds
      use newpara
      use ccounters
      implicit none
c
      integer ncs,ncf,nval,nmo,nfirst,nlast,iprnt,nbf
      integer ndisktr,npairs,inx(12,*),nslv1
      logical vorb
      real*8  thresh
c
      integer ntrincf,idimen,iout,ischwarz,ndisk1,ics,jcs,kcs,lcs
      integer igetival,ijchunk1,ndisk2,ikcs,ibigress,ibigresa,ishpairs
      integer ijpasses1,iparts,iparta,ibuf1,ibuf2,ibuf3,isym,iant
      integer iisizeres,isize,jsize,ksize,lsize,iksize,jlsize,ijpass
      integer ijstart,ijstop,iress,iresa,irecs,ireca,jlcs,iresult
      integer ijsize,istart,istop,kstart,kstop,ij,ionemat,iikk,ii,kk
      integer, parameter :: igsize=15
      integer jlchunkmax,itot_int
      integer :: ikchunkmax,mult_memory,ikchunkfinal
      integer ikchunkno1,itail,memavail,ikchunktry
      integer jlarray(4,ncs*ncs),jlchunkno
      integer ikarray(4,ncs*ncs),ikchunkno
      integer iisizeamp,jlchunk,jlstart,jlstop,ikchunk,ikstart,ikstop
      integer ikcc,iskip,ifirst,jfirst,jlchunkfinal,jlchunktry
      integer jlchunkno1,iscr_s,iscr_a,nskip,nskip1,irowskip,irowskip1
      logical af
c     optimization common
      integer i,ist,info,isl_id,ijchunk2,ijpasses2,itail1
      integer inxold(12,ncs),inxnew(12,ncs),ish_tab_ik,ish_tab_jl
      integer iatpair,na
      real*8 timings(20),t0,t1,tt0,tt1,el_ik,el_ki,xmax
      integer list1(3,10),iip,kkp,iii,iorb,jorb,itr,islgid,itr1
c
c
      iout=igetival('iout')
      if (iprnt.gt.3) then
       write(iout,*) '                            The EEOij Module  '
       write(iout,*) ' EEOij integral thresh    = ',thresh
       write(iout,*)
       call flush(iout)
      endif
      af=.true.
c 1=int,2=mult,3=IO
      call dynamic_mmark()
      call matmark()
c allocate memory for integer*2 arrays:
      call dynamic_getmem((ncs*(ncs+1)/2+1)/2    +1,ish_tab_ik)
      call dynamic_getmem((ncs*(ncs+1)/2+1)/2    +1,ish_tab_jl)
c allocate memory for integer*4 array:
      call dynamic_getmem((ncs*(ncs+1)/2+1)*3/2  +1,iatpair)
c
      call dynamic_getmem(ncf*(ncf+1)/2,iscr_s)
      call dynamic_getmem(ncf*(ncf+1)/2,iscr_a)
      call zeroit(bl(iscr_s),ncf*(ncf+1)/2)
      call zeroit(bl(iscr_a),ncf*(ncf+1)/2)
c
      do i=1,20
        timings(i)=0d0
      enddo
      call generate_ikpairs(bl(ish_tab_ik),ncs)
      call getival('na',na)
      call generate_jlpairs(bl(ish_tab_jl),ncs,inx,na,bl(iatpair),iprnt)
      ijchunk1=min(150,npairs)
      ijchunk2=50
      ikchunkmax=207
      jlchunkmax=500
c Determine the max size of ikchunkmax
c 
      call elapsec(tt0)
      call dynamic_show_free(memavail)
c
c     Mam size of ij bunch of amplitudes, as limited by the prepare_ampl subr.
      itr=(memavail-4*ncf*ncf)/(ncf*(ncf+1)+2*jlchunkmax)
      itr1=npairs/nslv
      if (mod(npairs,nslv).ne.0) itr1=itr1+1
      if (itr1.gt.itr) then
        ijchunk1=itr
      else if (itr1.gt.ijchunk1) then
        ijchunk1=itr1
      endif
      if (ijchunk1.gt.itr) ijchunk1=itr
c     Mem size of ij bunch of amplitudes, as limited by the results collection
      itr=(memavail)/(2*ncf*ncf)
      if (itr1.gt.itr) then
        ijchunk2=itr
      else if (itr1.gt.ijchunk2) then
        ijchunk2=itr1
      endif
c
      ijpasses1=npairs/ijchunk1
      if (mod(npairs,ijchunk1).ne.0) ijpasses1=ijpasses1+1
c
      ijpasses2=npairs/ijchunk2
      if (mod(npairs,ijchunk2).ne.0) ijpasses2=ijpasses2+1
      ikchunkmax=ncf*(ncf+1)/2
c Calculate max. possible ikchunkmax
      do ikchunktry=ikchunkmax,1,-1
        mult_memory=2*ijchunk1*jlchunkmax+         ! chunks of amplit.
     *              3*igsize**4+                   ! integral buffers
     *              2*jlchunkmax*ikchunktry+       ! integral chunk
     *              2*ijchunk1*ikchunktry+         ! partial res for screening
     *              2*max(ijpasses1*ijchunk1*ikchunktry,
     *                    ijpasses2*ijchunk2*ikchunktry)+ ! results chunk
     *              500000                         ! possible fragmentation
        if (mult_memory.lt.memavail) then
          ikchunkmax=ikchunktry
          exit
        endif
      enddo
      if (iprnt.gt.2) write(6,*) 'Memory ikchunkmax: ',ikchunkmax
c
c Calculate max. possible ikchunkmax (but less than previous one,
c because of memory limitations) such that it gives the number of chunks
c which is multiple of number of slaves
c
      do ikchunktry=ikchunkmax,1,-1
        call int_chunks_iktst(inx,ncs,ikchunktry,itail,ikchunkno,
     *                        bl(ish_tab_ik))
        if (mod(ikchunkno,nslv).eq.0.and.ikchunkno/nslv.ge.1) then
          ikchunkmax=ikchunktry
          exit
        endif
      enddo
      if (iprnt.gt.2)
     *  write(6,*) 'max ik chunk, no of chunks, tail size: ',
     *             ikchunkmax,ikchunkno,itail
c
c Calculate the ikchunkmax, for which the last chunk is the biggest,
c but the total number of chunks is the same, as calculated previously
c
      ikchunkfinal=ikchunkmax
      do ikchunktry=ikchunkmax,1,-1
        call int_chunks_iktst(inx,ncs,ikchunktry,itail1,ikchunkno1,
     *                        bl(ish_tab_ik))
        if (ikchunkno1.ne.ikchunkno) exit
        if (itail1.gt.itail) then
          ikchunkfinal=ikchunktry
        endif
      enddo
      ikchunkmax=ikchunkfinal
      call int_chunks_iktst(inx,ncs,ikchunkmax,itail,ikchunkno,
     *                      bl(ish_tab_ik))
      if (iprnt.gt.2) then
         write(6,*) 'max ik chunk, no of chunks, tail size: ',
     *               ikchunkmax,ikchunkno,itail
c
      write(6,*) 'The total info for memory and file management:'
      write(6,*) 'ijchunk of amplitudes size, record size: ',
     *         ijchunk1,ijchunk1*jlchunkmax*8,' bytes'
      write(6,*) 'ijchunk of results    size, record size: ',
     *         ijchunk2,ijchunk2*ikchunkmax*8,' bytes'
      write(6,*) 'jlchunkmax,ikchunkmax: ',jlchunkmax,ikchunkmax
      call flush(6)
c
      endif
      iskip=0
      ic_MAS_EEO_INT=ic_MAS_EEO_INT+1
      ntrincf=ncf*(ncf+1)/2
      ishpairs=ncs*(ncs+1)/2
      iisizeres=ijchunk2*ikchunkmax*8
      iisizeamp=ijchunk1*jlchunkmax*8
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      call  matsub('occu','cano',nfirst,nlast)
c
c
c Perform schwartz integrals calculation:
      call dynamic_getmem(ncs*ncs,ischwarz)
      call schwarz_calc(inx,ncs,bl(ischwarz))
c
c  Prepare array with jlchunk info
      call int_chunks_jl(inx,ncs,jlchunkmax,jlarray,jlchunkno,
     *                   bl(ish_tab_jl),bl(iatpair),iprnt)
      call int_chunks_ik(inx,ncs,ikchunkmax,ikarray,ikchunkno,
     *                   bl(ish_tab_ik))
c
      call para_initsend
      call para_pack_int(ijchunk1,1)
      call para_pack_int(ijchunk2,1)
      call para_pack_int(ikchunkmax,1)
      call para_pack_int(jlchunkmax,1)
      call para_pack_int(jlchunkno,1)
      call para_pack_int(ikchunkno,1)
      call para_pack_int(jlarray,4*jlchunkno)
      call para_pack_int(ikarray,4*ikchunkno)
      call para_pack_int4(bl(ish_tab_ik),ncs*(ncs+1)/2+2)
      call para_pack_int4(bl(ish_tab_jl),ncs*(ncs+1)/2+2)
      call para_bcast_pack(TxCCSDInit)
c
      ijpasses1=npairs/ijchunk1
      if (mod(npairs,ijchunk1).ne.0) ijpasses1=ijpasses1+1
      call elapsec(t0)
      call prepare_ampl_mas(inx,ijchunk1,ncf,nmo,ndisk1,af,vorb,
     *                  npairs,nval,ncs,jlchunkmax,jlarray,jlchunkno,
     *                  xmax,nslv,bl(iscr_s),bl(iscr_a),nbf)
      call elapsec(t1)
      timings(1)=timings(1)+t1-t0
      if (iprnt.gt.2) write(6,*) 'The largest amplitude: ',xmax
c
      call dynamic_mmark()
c
      ndisk2=43
      call fafreopen(ndisk1)
      call Open_Scratch_File(ndisk2,ikchunkmax,ijchunk2,af,'EEOb')
      call sync_barrier
c
      call para_initsend
      call para_pack_int(ndisk1,1)
      call para_pack_int(ndisk2,1)
      call para_bcast_pack(TxCCSDInit)
ccs      call sleep(5)
      call elapsec(tt1)
      timings(8)=timings(8)+tt1-tt0
c    HERE CALCULATE ics, kcs
c
      call elapsec(tt0)
      do ikchunk=1,ikchunkno
          ikstart=ikarray(1,ikchunk)
          ikstop =ikarray(2,ikchunk)
          iksize =ikarray(3,ikchunk)
          if (iprnt.gt.2) write(6,'(4I8)') ikchunk,ikstart,ikstop,iksize
          call flush(6)
          call para_recv(isl_id,islgid,TxCCSDReq)
          call para_initsend
          call para_pack_int(ikstart,1)
          call para_pack_int(ikstop ,1)
          call para_pack_int(iksize ,1)
          call para_pack_int(ikchunk ,1)
          call para_send_pack(islgid,TxCCSDJob)
      enddo
      ikstart=0
      ikstop=0
      iksize=0
      nskip=0
      irowskip=0
      do i=1,nslv
          call para_recv(isl_id,islgid,TxCCSDReq)
          call para_initsend
          call para_pack_int(ikstart,1)
          call para_send_pack(islgid,TxCCSDJob)
          call para_recv(nskip1,islgid,TxCCSDRes)
          call para_recv(irowskip1,islgid,TxCCSDRes)
          nskip=nskip+nskip1
          irowskip=irowskip+irowskip1
      enddo
      itot_int=2*(ncf*(ncf+1)/2)**2
      if (iprnt.gt.2) then
      write(6,'(A,I5)')    "Total amount skipped: ",nskip
      write(6,'(A,F5.1)') "Percent skipped:      ",
     *                    dble(nskip)/itot_int*1d2
      write(6,'(A,F5.1)') "Row percent skipped:      ",
     *                    dble(irowskip)/itot_int*1d2
      endif
c
ccs      call sleep(5)
      call elapsec(tt1)
      timings(9)=timings(9)+tt1-tt0
c
      call elapsec(tt0)
c
      call fafClosem(ndisk1,0,ist)
      call sync_barrier
      if (ic_MAS_EEO_INT.gt.1) then 
        call fafClosem(ndisktr,0,ist)
      else
        call fafClosem(ndisktr,-1,ist)
      endif
      ndisktr=64
      call sync_barrier
      call Open_Scratch_File(ndisktr,ncf,ncf,af,'EEOc')
      call sync_barrier
c
      call dynamic_retmark()
c
      call results_collector_mas(npairs,ijchunk2,ndisk2,ndisktr)
c
      call fafClosem(ndisk2,0,ist)
      call sync_barrier
      call para_barrier
      call para_reduce(timings,20,TxCCSDRedu2)
c
      if (iprnt.gt.2) write(6,*) 'Shell pairs skipped: ',iskip
      call dynamic_retmark()
      call matremark()
      call elapsec(tt1)
      timings(10)=timings(10)+tt1-tt0
 33   FORMAT (A,F10.3)
      if (iprnt.gt.1) then
      write(6,33) 'Amplitudes transformation: ', timings(1)
      write(6,33) 'Integral generation:       ', timings(2)
      write(6,33) 'Multiplication:            ', timings(3)
      write(6,33) 'Amplitudes reads:          ', timings(4)
      write(6,33) 'Results write:             ', timings(5)
      write(6,33) 'Results collection IO      ', timings(6)
      write(6,33) 'Screening:                 ', timings(7)
      write(6,33) 'Master elaps. part 1       ', timings(8)
      write(6,33) 'Master elaps. part 2       ', timings(9)
      write(6,33) 'Master elaps. part 3       ', timings(10)
      write(6,33) 'Twoel:                     ', timings(11)
      endif
      call para_initsend
      call para_pack_int(ndisktr,1)
      call para_bcast_pack(TxCCSDInit)
      end
C==============================================================================
      subroutine SLA_EEO_INT(ncs,    ncf,     inx,   nval,  nmo,
     *                       nfirst, nlast,   thresh,
     *                       vorb,   ndisktr, npairs,
     *                       iMasID, mygid,   mytid,  nslv1,nbf)
      use memory
      use kinds
      use newpara
      implicit none
      integer ncs,ncf,inx(12,*),nval,nmo,nfirst,nlast,ndisktr,npairs
      integer iMasID,mygid,mytid,nslv1,nbf
      real*8 thresh
      logical vorb
c
      logical af
      integer jlarray(4,ncs*ncs),ish_tab_jl,ish_tab_ik,iscr_s,iscr_a
      real*8 timings(20),xt0,xt1,t1,t0,xinttime,ikarray(4,ncs*ncs)
      integer info,ijchunk1,ijchunk2,ikchunkmax,jlchunkmax,jlchunkno
      integer ndisk1,ndisk2,ikchunkno,ikstart,ikstop,ikchunk,iksize
      integer ischwarz,i,iskip,imastgid,nskip,irowskip,ixi
c
      af=.true.
      iskip=0
      xinttime=0d0
      do i=1,20
        timings(i)=0d0
      enddo
c     call itfer(inx,inxold,12*ncs)
      call dynamic_mmark()
      call matmark()
c allocate memory for integer*2 arrays:
      call dynamic_getmem((ncs*(ncs+1)/2+1)/2    +1,ish_tab_ik)
      call dynamic_getmem((ncs*(ncs+1)/2+1)/2    +1,ish_tab_jl)
c
      call dynamic_getmem(ncf*(ncf+1)/2,iscr_s)
      call dynamic_getmem(ncf*(ncf+1)/2,iscr_a)
      call zeroit(bl(iscr_s),ncf*(ncf+1)/2)
      call zeroit(bl(iscr_a),ncf*(ncf+1)/2)
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ijchunk1,1)
      call para_unpack_int(ijchunk2,1)
      call para_unpack_int(ikchunkmax,1)
      call para_unpack_int(jlchunkmax,1)
      call para_unpack_int(jlchunkno,1)
      call para_unpack_int(ikchunkno,1)
      call para_unpack_int(jlarray,4*jlchunkno)
      call para_unpack_int(ikarray,4*ikchunkno)
      call para_unpack_int4(bl(ish_tab_ik),ncs*(ncs+1)/2+2)
      call para_unpack_int4(bl(ish_tab_jl),ncs*(ncs+1)/2+2)
c 
      call  matsub('occu','cano',nfirst,nlast)
c 
      call prepare_ampl_sla(inx,ijchunk1,ncf,nmo,af,vorb,
     *                  npairs,nval,ncs,jlchunkmax,jlarray,jlchunkno,
     *                  iMasID,mygid,mytid,nslv,bl(iscr_s),bl(iscr_a),
     *                  bl(ish_tab_jl),nbf)
c 
      call fafCreates(ndisk2)
      call sync_barrier
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ndisk1,1)
      call para_unpack_int(ndisk2,1)
      call fafreopen(ndisk1)
c
c Perform schwartz integrals calculation:
      call dynamic_getmem(ncs*ncs,ischwarz)
      call schwarz_calc(inx,ncs,bl(ischwarz))
c
      timings(1)=timings(1)+t1-t0
c
      call dynamic_mmark()
c
      call elapsec(xt0)
      nskip=0
      irowskip=0
      do
        call para_send(mytid,0,TxCCSDReq)
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(ikstart,1)
        if (ikstart.eq.0) then 
          call para_send(nskip,0,TxCCSDRes)
          call para_send(irowskip,0,TxCCSDRes)
          exit
        endif
        call para_unpack_int(ikstop ,1)
        call para_unpack_int(iksize ,1)
        call para_unpack_int(ikchunk ,1)
c
        call jl_loop(bl(ish_tab_ik),inx,  ncs,        thresh,  npairs,
     *               ischwarz, jlchunkno, jlchunkmax, jlarray, ijchunk1,
     *               ijchunk2, ikchunkno, ikchunkmax, ikstart, ikstop,
     *               iksize,   ikchunk,   ndisk1 ,    ndisk2 , af,
     *               iskip,    xinttime,  timings,bl(iscr_s),bl(iscr_a),
     *               nskip,    irowskip,  bl(ish_tab_jl))
      enddo
      call elapsec(xt1)
      timings(11)=xinttime
      call fafCloses
      call sync_barrier
      call fafCloses
      call sync_barrier
      call fafCreates(ixi)
      call sync_barrier
      call results_collector_sla(ncf,af,ikchunkno,ikchunkmax,
     *                           timings,inx,bl(ish_tab_ik),ikarray)
      call fafcloses
      call sync_barrier
      call para_barrier
      call para_reduce(timings,20,TxCCSDRedu2)
c
      call dynamic_retmark()
c
      call dynamic_retmark()
      call matremark()
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ndisktr,1)
c     call itfer(inxold,inx,12*ncs)
c     call twoelinit
      end
C==============================================================================
      subroutine prepare_ampl_mas(inx,ijchunkmax,ncf,nmo,ndisk,
     *                        af,vorb,npairs,nval,ncs,
     *                        jlchunkmax,jlarray,jlchunkno,xmax,
     *                        nslv1,scr_s,scr_a)
      use memory
      use kinds
      use newpara
      implicit none
      integer inx(12,*),ijchunkmax,ncf,nmo,ndisk,npairs,nval,ncs
      integer jlchunkmax,jlarray(4,*),jlchunkno,nslv1
      logical af,vorb
      real*8 xmax,scr_s(*),scr_a(*)
c
      integer idimen,ijpasses,ijpass,iistart,iistop,ijsize,i,info
      integer isltid,islgid
      real*8 xmax1
c
      ndisk=78
      xmax=0d0
      call Open_Scratch_File(ndisk,jlchunkmax,ijchunkmax,af,'EEOa')
      call sync_barrier
      ijpasses=npairs/ijchunkmax
      if (mod(npairs,ijchunkmax).ne.0) ijpasses=ijpasses+1
c
      call para_initsend
      call para_pack_int(ndisk,1)
      call para_bcast_pack(TxTriplesInit)
c
      do ijpass=1,ijpasses
        iistart= (ijpass-1)*ijchunkmax+1
        iistop = (ijpass  )*ijchunkmax
        if (iistart.gt.npairs) STOP '4rdfgtgbcs'
        if (iistop .gt.npairs) iistop=npairs
        ijsize=iistop-iistart+1
        call para_recv(isltid,islgid,TxCCSDReq2)
c
        call para_initsend
        call para_pack_int(iistart,1)
        call para_pack_int(iistop,1)
        call para_pack_int(ijpass,1)
        call para_pack_int(ijsize,1)
        call para_send_pack(islgid,TxCCSDJob2)
c
      enddo
      iistart=0
      do i=1,nslv
        call para_recv(isltid,islgid,TxCCSDReq2)
c
        call para_initsend
        call para_pack_int(iistart,1)
        call para_send_pack(islgid,TxCCSDJob2)
c
      enddo
      call sync_barrier
      call para_reduce_max(xmax,1,TxCCSDRes3)
      call para_reduce_max(scr_s,ncf*(ncf+1)/2,TxCCSDRes3)
      call para_reduce_max(scr_a,ncf*(ncf+1)/2,TxCCSDRes3)
      call para_bcast_real(scr_s,ncf*(ncf+1)/2,TxCCSDRes2)
      call para_bcast_real(scr_a,ncf*(ncf+1)/2,TxCCSDRes2)
      end
c===============================================================================
      subroutine prepare_ampl_sla(inx,ijchunkmax,ncf,nmo,
     *                        af,vorb,npairs,nval,ncs,
     *                        jlchunkmax,jlarray,jlchunkno,
     *                        iMasID, mygid,   mytid,  nslv1,
     *                        scr_s,scr_a,ish_tab,nbf)
      use memory
      use kinds
      use newpara
      implicit none
      integer inx(12,*),ijchunkmax,ncf,nmo,npairs,nval,ncs,nbf
      integer jlchunkmax,jlarray(4,*),jlchunkno
      integer iMasID,mygid,mytid,nslv1
      logical af,vorb
      real*8 scr_s(*),scr_a(*)
      integer*2 ish_tab(2,*)
c
      integer idimen,ijpasses,ijpass,iistart,iistop,ijsize,ndisk,info
      integer imastgid
      real*8 xmax
c
      call dynamic_mmark()
      call matmark()
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      xmax=0d0
c
      call fafcreates(ndisk)
      call sync_barrier
      call para_recv_bcastpack(TxTriplesInit)
      call para_unpack_int(ndisk,1)
c
      do
        call para_send(mytid,0,TxCCSDReq2)
c
        call para_recv_pack(imastgid,TxCCSDJob2)
        call para_unpack_int(iistart,1)
        if (iistart.eq.0) then 
          exit
        endif
        call para_unpack_int(iistop,1)
        call para_unpack_int(ijpass,1)
        call para_unpack_int(ijsize,1)

        call singl_ij(nval,  ncf,   ncs,       idimen,   inx,
     *                af,    vorb,  ndisk,     iistart,  iistop,
     *                ijpass,ijsize,ijchunkmax,jlchunkno,jlchunkmax,
     *                jlarray, xmax,scr_s,     scr_a,    ish_tab)
      enddo
      call dynamic_retmark()
      call matremark()
      call sync_barrier
      call para_reduce_max(xmax,1,TxCCSDRes3)
      call para_reduce_max(scr_s,ncf*(ncf+1)/2,TxCCSDRes3)
      call para_reduce_max(scr_a,ncf*(ncf+1)/2,TxCCSDRes3)
      call para_recv_bcastreal(scr_s,ncf*(ncf+1)/2,TxCCSDRes2)
      call para_recv_bcastreal(scr_a,ncf*(ncf+1)/2,TxCCSDRes2)
      end
c===============================================================================
      subroutine results_collector_mas(npairs,ijchunk,ndisks,ndiskd)
      use memory
      use newpara
      implicit none
      integer npairs,ijchunk,ndisks,ndiskd
c
      integer ijpasses,ijpass,ijstart,ijstop,ijsize,i,islgid
c
c
      ijpasses=npairs/ijchunk
      if (mod(npairs,ijchunk).ne.0) ijpasses=ijpasses+1
      call para_bcast(ndisks,MP2file)
      call para_bcast(ndiskd,MP2file)
      call para_bcast(ijchunk,MP2file)
c
      do ijpass=1,ijpasses
        ijstart= (ijpass-1)*ijchunk+1
        ijstop = (ijpass  )*ijchunk
        if (ijstart.gt.npairs) STOP '4rdfgtgbcs'
        if (ijstop .gt.npairs) ijstop=npairs
        ijsize=ijstop-ijstart+1
        call para_recv(i,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(ijstart,1)
        call para_pack_int(ijstop ,1)
        call para_pack_int(ijsize ,1)
        call para_pack_int(ijpass ,1)
        call para_send_pack(islgid,TxCCSDJob3)
      enddo
      ijstart=0
      do ijpass=1,nslv
        call para_recv(i,islgid,TxCCSDReq3)
        call para_initsend
        call para_pack_int(ijstart,1)
        call para_send_pack(islgid,TxCCSDJob3)
      enddo
      end
c===============================================================================
      subroutine results_collector_sla(ncf,af,ikchunkno,ikchunkmax,
     *                                 timings,inx,ish_tab_ik,ikarray)
      use memory
      use newpara
      implicit none
      integer ncf,ikchunkno,ikchunkmax
      integer inx(12,*),ish_tab_ik(2,*),ikarray(4,*)
      logical af
      real*8 timings(*)
c
      integer ijstart,ijstop,ijsize,ndisks,ndiskd,ijpass,i,ijchunk
c
      call para_recv_bcast(ndisks,MP2file)
      call para_recv_bcast(ndiskd,MP2file)
      call para_recv_bcast(ijchunk,MP2file)
c
      i=0
      do
        call para_send(i,0,TxCCSDReq3)
        call para_recv_pack(i,TxCCSDJob3)
        if (i.ne.0) STOP 'Error 94y3jhlshgghfgdkfghghgh55g5'
        call para_unpack_int(ijstart,1)
        if (ijstart.eq.0) exit
        call para_unpack_int(ijstop ,1)
        call para_unpack_int(ijsize ,1)
        call para_unpack_int(ijpass ,1)
        call ij_res(ncf, af, inx, ish_tab_ik, ikchunkno,
     *              ikchunkmax, ikarray, ndisks, ndiskd, ijsize,
     *              ijpass, ijstart, ijstop, ijchunk,  timings)
      enddo
      end
