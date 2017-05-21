#include "maxi.h"
c     subroutine do_CCSD(iMasID,mygid,mytid)
      subroutine do_CCSD()
      use memory
      use kinds
      use newpara
      implicit none
      integer iMasID,mygid,mytid
c     
      logical nodisk,gauss_seidel,loca,omit,mp2,mp3,cep0,cep2,smal
      logical small,cc,norecalc,byt8,log_diis,singles,ccsd,vorb,nofr
      integer nmo,nvirt,nfirst,nlast,nval,idimen,npairs,icano,iepsi
      integer maxsh,na,ncf,ncs,nbfx,nsh,ibas,inuc,ictr,nsym,ngener
      integer iadr,nupair,nupair1,ifp,ifp1,mataddr,i
      integer ipass,IsWork,lbinc, lrecc, ndiskc, irecadrc
      integer              lbinx, lrecx, ndiskx, irecadrx
      integer              lbine, lrece, ndiske, irecadre
      integer              lbintx,lrectx,ndisktx,irecadrtx
      integer              lbintc,lrectc,ndisktc,irecadrtc
      integer              lbintt,lrectt,ndisktt,irecadrtt
      integer iterate,ndiskYZ,ilist,ising_res,iX,ndiskdiisr
      integer ndiskdiisc,ndiskalpha,length_diis
      logical af,l_triples,qcisd,restart
      real*8  thresh,xmaxdisk,core,shift,shiftini,erhf,tot_xmax
      character*20 chtid
      character*3 chgid,chite
      character*256 hostname
      character*2 int_kind
      integer islaves,ndisktrc,ndisktrx,ndisktre,ndisktrtx,ndisktrtc
      integer ndisktrtt
      integer mygidcp
      common /extrone_tmp/ af,islaves,ndisktrc,ndisktrx,ndisktre,
     *                     ndisktrtx,ndisktrtc,ndisktrtt,mygidcp,nbf
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      logical equal,fileopen
      integer iresAO,iexchAO,iwork1
      integer iwork2,iwork3,iresult,iwork,j,ij,icoeAO,iresult_iter,ibeta
      integer ifockAO,ifockMO,ioverlap,icorefock,icoeffAO,kijklndisk
      integer ndiskg41a,inum,iii,istart,istop,jstart,jstop,ndiskQ
      integer max_dyn_mem,iteration,iad,lcore,memavail,int_array
      real*8 epair,epairr,sqij,xmaxele,xepair,xepairr,xnormij
      real*8 s_ext_matstop,e_ext_matstop,e_ext_mat_tot1,e_ext_matstart
      real*8 eqstart,eqstop,sqstart,sqstop,s_ext_mat_tot1,eq_tot1
      real*8 sq_tot1,egstart,sgstart,egstop,sgstop,eg_tot1,sg_tot1
      real*8 sco_tot1,eco_tot1,ecostop,scostop,scostart,ecostart
      real*8 sen_tot1,een_tot1,t0,t1,eenstop,senstop,CPU_doubles_singles
      real*8 eenstart,senstart,xiterenergy
      real*8 s_ext_matstart,elaps(8),slave_cpu,slave_cpu0,slave_cpu1
      real*8 slave_tot_cpu0,slave_tot_cpu1,slave_tot_cpu
      character*256 home_dir,scrfile
      logical mp4,do_mp4
      integer len,istatus,ndisk_mp3,ndiskmp2,iamount
      real*8 strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
c TEST
      real*8 zeroed,zeroed_max,ext3_ratio
      common /zero_in_EEO/ zeroed,zeroed_max
      integer :: info, ibuffer, nstrong
      integer :: isympairs,ipairimages,ivalpair,ivirpair,ij_unique
      integer iad2,iad3,iinvcano,diis,itmp_mat1,itmp_mat2,icode
      integer ncore,iiR,iiS,iotable,ivtable,ivrevtable,im_table,ijtable
      integer ichar_count,iorevtable,imastgid,itrala,nbf
      real*8 det,sepair,sepairr,sxepair,sxepairr
      logical docansym
c
      zeroed=0.0d0
      zeroed_max=0.0d0
      info = 0
      ibuffer = 0
      ichar_count = 1
c
      call fafinit(iresult)
      call reset_counters
      if (iresult.eq.1) return ! I have been an I/O daemon
      if (iresult.ne.0.and.iresult.ne.1)
     * STOP 'fafinit error'
      iMasID=MASTER_ID
      mygid=MY_GID
      mytid=MY_ID
      call secund(slave_tot_cpu0)
      call setival('ccprint',0)
      af=.true.
      do_mp4=.false.
      mygidcp=mygid
c     call getenv('HOME',home_dir)
c     call rmblan(home_dir,80,len) !returns len of string without spaces,
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
      write (chtid,'(I10.10)') mytid
      open(99,FILE=scrfile(1:len)//'.debug_slaves'//chtid)
      write(chgid,'(I3.3)') mygid
c
c Memory markers:
      call matreset
      call mmark
      call matmark
c
      call fafCreates(itrala) ! peer of masters open_mp3_resid
c     Get basic information about geometry etc.
c     call para_recv_bcastpack(TxScfInit)
      call para_get_info
c
c Extract data from para_get_info
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('nbf',nbfx)
      call getival('nsh',nsh)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('ictr',ictr)
      call getival('nsym',nsym)
      if (nsym.gt.0) then
        call getival('ngener',ngener)
        call getival('nsyo',iadr)
        call getival('SymNuPr',nupair)
        call getival('SymNuPr1',nupair1)
        call getival('SymFunPr',ifp)
        call getival('SymShPr',ifp1)
      endif
      call setival('iout',99)
      call setival('icond',99)
c
c Receive symmetry info:
      
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(nbf,1)
      call para_unpack_int(max_dyn_mem,1)
      call dynamic_init(max_dyn_mem,bl(1))
      call para_unpack_int(npairs,1)
      call dynamic_getmem(npairs,isympairs)
      call dynamic_getmem(npairs+1,ipairimages)
      call para_unpack_int(bl(isympairs),npairs)
      call para_unpack_int(bl(ipairimages),npairs+1)
      call para_unpack_int(nval,1)
      call para_unpack_int(nvirt,1)
      call dynamic_getmem(7*nval, ivalpair)
      call dynamic_getmem(7*nvirt,ivirpair)
      call para_unpack_int(bl(ivalpair),7*nval)
      call para_unpack_int(bl(ivirpair),7*nvirt)
      call para_unpack_int(ij_unique,1)
      call para_unpack_int(loca,1)
      call para_unpack_int(docansym,1)
      if (.not.loca) then
        call dynamic_getmem(8*8,im_table)
        call dynamic_getmem(8*(nval+1),iotable)
        call dynamic_getmem(8*(nvirt+1),ivtable)
        call dynamic_getmem(nval,iorevtable)
        call dynamic_getmem(nvirt,ivrevtable)
        call dynamic_getmem((nval*nval+1)*8,ijtable)
        call para_unpack_int(bl(iotable),8*(nval+1))
        call para_unpack_int(bl(ivtable),8*(nvirt+1))
        call para_unpack_int(bl(iorevtable),nval)
        call para_unpack_int(bl(ivrevtable),nvirt)
        call para_unpack_int(bl(im_table),8*8)
        call para_unpack_int(bl(ijtable),(nval*nval+1)*8)
        call para_unpack_int(ichar_count,1)
      endif
c CCSD data:
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(nodisk,1)
      call para_unpack_int(gauss_seidel,1)
      call para_unpack_int(omit,1)
      call para_unpack_int(mp2,1)
      call para_unpack_int(mp3,1)
      call para_unpack_int(mp4,1)
      call para_unpack_int(cep0,1)
      call para_unpack_int(cep2,1)
      call para_unpack_int(smal,1)
      call para_unpack_int(small,1)
      call para_unpack_int(cc,1)
      call para_unpack_int(norecalc,1)
      call para_unpack_int(byt8,1)
      call para_unpack_int(log_diis,1)
      call para_unpack_int(singles,1)
      call para_unpack_int(ccsd,1)
      call para_unpack_int(qcisd,1)
      call para_unpack_int(vorb,1)
      call para_unpack_int(nofr,1)
      call para_unpack_int(cache,1)
      call para_unpack_int(integerT,1)
      call para_unpack_int(l_triples,1)
      call para_unpack_int(restart,1)
c int:
      call para_unpack_int(nmo,1)
      call para_unpack_int(nvirt,1)
      call para_unpack_int(nfirst,1)
      call para_unpack_int(nlast,1)
      call para_unpack_int(nval,1)
      call para_unpack_int(idimen,1)
      call para_unpack_int(npairs,1)
      call para_unpack_int(length_diis,1)
      call para_unpack_int(ndisk_mp3,1)
      call para_unpack_int(ncore,1)
      call dynamic_mmark
c real
      call para_unpack_real(erhf,1)
      call para_unpack_real(thresh,1)
      call para_unpack_real(xmaxdisk,1)
      call para_unpack_real(core,1)
      call para_unpack_real(shift,1)
      call para_unpack_real(shiftini,1)
      call dynamic_matdef('cano','q',ncf,ncf)
      call dynamic_matdef('epsi','d',ncf,ncf)
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
      call para_unpack_real(bl(icano),ncf*ncf)
      call para_unpack_real(bl(iepsi),ncf)
      islaves=nslv
      call setival('nslv',nslv)
      call matsub('genvirt','cano',nmo+1,nbf)
      call dynamic_lock(bl(icano),i)
      call dynamic_lock(bl(iepsi),i)
c
c Initialize symmetry commons
      call symm_memory_allocator(ncf)
      call pair_s_initializer(vorb,     nval,     nsym, ifp, ifp1,
     *                        ivalpair, ivirpair, loca,iorevtable,
     *                        ivrevtable,im_table)
c
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadrc)
      call izeroit(bl(irecadrc),npairs*nslv)    !zero out the pair records
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadrx)
      call izeroit(bl(irecadrx),npairs*nslv)    !zero out the pair records
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadre)
      call izeroit(bl(irecadre),npairs*nslv)         !zero out the pair records
      call dynamic_getmem((nval*nval/intsize+1)*nslv,irecadrtx)
      call izeroit(bl(irecadrtx),nval*nval*nslv)     !zero out the pair records
      call dynamic_getmem((nval*nval/intsize+1)*nslv,irecadrtc)
      call izeroit(bl(irecadrtc),nval*nval*nslv)   !zero out the pair records
      call dynamic_getmem((nval*(nbf-nmo)/intsize+1)*nslv,irecadrtt)
      call izeroit(bl(irecadrtt),nval*(nbf-nmo)*nslv)!zero out the pair records
c
      call twoelinit ! initalize callable two-electron program
      call getint(npairs,ilist)
      call filllist(bl(ilist),npairs)
c
      int_kind='c'
      ipass=0
      if (.not.mp2.and..not.(restart.and.l_triples)) then
      call slave_Gen(ncs,     ncf,          ictr,    nval,   nmo,
     *               nfirst,  nlast,        thresh,  core,   xmaxdisk,
     *               ndiskc,  bl(irecadrc), lrecc,   lbinc,  int_kind,
     *               nodisk,  byt8,         small,   vorb,   iMasID,
     *               mygid,   mytid,        ndisktrc,tot_xmax,bl(ilist),
     *               nstrong,bl(isympairs),ij_unique,nbf)
      endif
      call izeroit(bl(irecadrx),npairs*nslv)
      int_kind='x'
      ipass=0
      call slave_Gen(ncs,     ncf,          ictr,    nval,   nmo,
     *               nfirst,  nlast,        thresh,  core,   xmaxdisk,
     *               ndiskx,  bl(irecadrx), lrecx,   lbinx,  int_kind,
     *               nodisk,  byt8,         small,   vorb,   iMasID,
     *               mygid,   mytid,        ndisktrx,tot_xmax,bl(ilist),
     *               nstrong,bl(isympairs),ij_unique,nbf)
      if (l_triples) then
      call izeroit(bl(irecadrtt),nval*(nbf-nmo)*nslv)
      int_kind='tt'
      ipass=0
      call slave_Gen(ncs,    ncf,          ictr,     nval,    nmo,
     *               nfirst, nlast,        thresh,   core,    xmaxdisk,
     *               ndisktt,bl(irecadrtt),lrectt,   lbintt,  int_kind,
     *               nodisk, byt8,         small,    vorb,    iMasID,
     *               mygid,  mytid,        ndisktrtt,tot_xmax,bl(ilist),
     *               nstrong,bl(isympairs),ij_unique,nbf)
      endif
c
c  Initialization:
c
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
      call dynamic_matdef('fockAO','q',idimen,idimen)
      ifockAO=mataddr('fockAO')
      call dynamic_getmem(nbf*nbf,ifockMO)
      call dynamic_matdef('overlap','q',ncf,ncf)
      ioverlap=mataddr('overlap')
      call dynamic_matdef('corefock','q',idimen,idimen)
      icorefock=mataddr('corefock')
c
      if (.not.mp2) then
      call Kijk_Vec_Slave(bl(irecadrx),npairs,ndiskx,lbinx,thresh,
     *                          byt8,kijklndisk,ncf,nmo,nval,
     *                          vorb,af,mygid,nfirst,nlast,
     *                          nslv,iMasID,mytid,nbf)
      endif
      call sync_barrier
c
      call para_recv_bcastreal(bl(ifockAO),idimen*idimen,
     *                               TxCCSDInit)
      call para_recv_bcastreal(bl(ifockMO),nbf*nbf,TxCCSDInit)
      call para_recv_bcastreal(bl(ioverlap),ncf*ncf,TxCCSDInit)
      call para_recv_bcastreal(bl(icorefock),idimen*idimen,
     *                               TxCCSDInit)
      call dynamic_lock(bl(icorefock),iamount)
      call dynamic_lock(bl(ifockAO),iamount)
      call dynamic_lock(bl(ifockMO),iamount)
c
      if (.not.mp2) then
      call fafCreates(itrala)
      call sync_barrier
      call KijklInit_Slave(bl(irecadrx),npairs,ndiskx,lbinx,ncf,
     *                     thresh,nfirst,nlast,bl(ifockMO),byt8,
     *                     nmo,vorb,af,nslv,iMasID,mytid)
      call sync_barrier
      endif
c
      call Fock_Vector_Init(ncf,nfirst,nlast,ifockAO,ifockMO,nmo,vorb,
     *                      nbf)
      call initsingles_slave(nval,ncf,singles.or.mp4,nmo,vorb,iMasID,
     *                       nbf)
      call initS(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call initR(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call mmark
      call dynamic_mmark
      call matmark
c
c Space for one K_ij matrix
      call dynamic_matdef('resao','q',idimen,idimen)
      iresAO=mataddr('resao')
      call dynamic_matdef('exchao','q',idimen,idimen)
      iexchAO=mataddr('exchao')
      call dynamic_matdef('invcano','q',ncf,ncf)
      iinvcano=mataddr('invcano')
      call matzero('invcano')
      call matcopy('cano','invcano')
      call matmmul2('cano','overlap','invcano','t','n','n')
c
ctmp  call dynamic_getmem(ncf,iad2)
ctmp  call dynamic_getmem(ncf,iad3)
ctmp  call osinv(bl(iinvcano),ncf,det,1.0d-18,bl(iad2),bl(iad3))
ctmp  call dynamic_retmem(2)
ctmp  if (dabs(det).lt.1d-200) then
ctmp    print *,'DET: ',det
ctmp    STOP 'Determinant too small'
ctmp  endif
c
      call dynamic_matdef('rinvcano','q',ncf,ncf)
      call matcopy('invcano','rinvcano') !store r(eal) inverted cano for singles
      call dynamic_matdef('canotran','q',ncf,ncf)
      call matcopy('cano','canotran')
      call matpose('canotran')
c  This is necessary because of transf: T_MO=U-1 T_AO U-1(T)
c I use later matsimtr, which makes: U(T)A U
      call matpose('invcano')
c locking:
      i=mataddr('invcano')
      call dynamic_lock(bl(i),j)
      i=mataddr('rinvcano')
      call dynamic_lock(bl(i),j)
      i=mataddr('canotran')
      call dynamic_lock(bl(i),j)
c
      call dynamic_matdef('work1','q',idimen,idimen)
      iwork1=mataddr('work1')
      call dynamic_matdef('work2','q',idimen,idimen)
      iwork2=mataddr('work2')
      call dynamic_matdef('work3','q',idimen,idimen)
      iwork3=mataddr('work3')
      call dynamic_matdef('beta','q',nval,nval)
      ibeta=mataddr('beta')
      if (log_diis) then
      call fafcreates(ndiskdiisr)
      call fafcreates(ndiskdiisc)
      endif
      call sync_barrier
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ndiskdiisr,1)
      call para_unpack_int(ndiskdiisc,1)
      if (.not.log_diis) then
      call CCDiis(i,       j,          'nodiis',     ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      else
      i=length_diis
      call CCDiis(i,       j,          'diis',       ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      endif
c
c  Calculation:
c
      if (restart) goto 126
      iresult=0
      do
        call para_initsend
        call para_pack_int(mytid,1)
        call para_pack_int(iresult,1)
        if (iresult.eq.1) then
          call para_pack_real(epair,1)
          call para_pack_real(sepair,1)
          call para_pack_real(sqij,1)
          call para_pack_real(xmaxele,1)
        endif
        call para_send_pack(0,TxCCSDMP2Req)
        xmaxele=0.0d0
c
        call para_recv_pack(imastgid,TxCCSDMP2Job)
        call para_unpack_int(iwork,1)
        if (iwork.eq.0) exit
        call para_unpack_int(i,1)
        call para_unpack_int(j,1)
        call para_unpack_int(ij,1)
c
c calculation, the kernel from master:
c
        equal=.false.
        if (i.eq.j) equal=.true.
        call ExtrOne(i,    j,     bl(irecadrx),  npairs,     ndiskx,
     *               ncf,  lbinx, thresh,        byt8,       'x',
     *               'mo', nmo,   vorb,          bl(iexchAO))
        call CoefRead('tt',i,j,icoeAO)
        call UpdateCoef(i,       j,     icoeAO,  iexchAO,    bl(iepsi),
     *                  ncf,     nmo,   nfirst,  npairs,     shiftini,
     *                  'nodiis',vorb,  singles, ndiskdiisr, ndiskdiisc,
     *                  sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
        call CoefWrite(i,j,icoeAO)
        call PairEnergyMP2(icoeAO,iexchAO,iresAO,epairr,epair,
     *                     ncf,   equal,  nmo,   vorb,  sepairr,
     *                     sepair,nbf)
c       print '(A15,I5,5X,F20.10)', "Pair energy: ",ij,epair
        iresult=1
c
      enddo
      if (mp2.and..not.loca) then
        goto 999
      endif
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
c Local MP2:
 222  continue
      do             ! iterative loop
      call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
        iresult=0
        call fafCreates(ndiskmp2)
        call sync_barrier
        call para_recv_bcast(ndiskmp2,MP2file)
        call Generate_MP2_G_SLA(nval,ncf,ndiskmp2,nmo,nfirst,vorb,
     *                         .true.,bl(ifockMO),nslv,iMasID,mytid,nbf)
        call sync_barrier
        do           ! ij loop
          call para_initsend
          call para_pack_int(mytid,1)
          call para_pack_int(iresult,1)
          if (iresult.eq.1) then
            call para_pack_int(i,1)
            call para_pack_int(j,1)
            call para_pack_int(ij,1)
            call para_pack_real(epair,1)
            call para_pack_real(epairr,1)
            call para_pack_real(xepair,1)
            call para_pack_real(sepair,1)
            call para_pack_real(sepairr,1)
            call para_pack_real(sxepair,1)
            call para_pack_real(sqij,1)
            call para_pack_real(xmaxele,1)
          endif
          call para_send_pack(0,TxCCSDMP2Req)
          xmaxele=0.0d0
c
          call para_recv_pack(imastgid,TxCCSDMP2Job)
          call para_unpack_int(iwork,1)
          if (iwork.eq.0) exit
          call para_unpack_int(i,1)
          call para_unpack_int(j,1)
          call para_unpack_int(ij,1)
c
c calculation, the kernel from master:
c
          equal=.false.
          if (i.eq.j) equal=.true.
          if (do_mp4) then
            call fafread(ndisk_mp3,bl(iexchAO),8,idimen*idimen,1,ij,
     *                   istatus)
          if (istatus.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
          else
          call ExtrOne(i,    j,     bl(irecadrx),  npairs,     ndiskx,
     *                 ncf,  lbinx, thresh,        byt8,       'x',
     *                 'mo', nmo,   vorb,          bl(iexchAO))
          endif
          call matcopy('exchao','resao')
          call CoefRead('tt',i,j,icoeffAO)
          call matconn('coefAO','q',idimen,idimen,icoeffAO)
          if (vorb) then
          call matmmul2('fockAO','coefAO','resao','n','n','a')
          call matmmul2('coefAO','fockAO','resao','n','n','a')
          else
          call matmmul2('fockAO','coefAO','work1','n','n','n')
          call matmmul2('work1','overlap','resao','n','n','a')
          call matmmul2('coefAO','fockAO','work1','n','n','n')
          call matmmul2('overlap','work1','resao','n','n','a')
          endif
          call matdisc('coefAO')
c         call Ggen_MP2(i,j,ncf,nval,bl(ifockMO),nfirst,iwork1,nmo,vorb)
          call Read_MP2(i,j,nmo,ncf,nval,ndiskmp2,vorb,.true.,iwork1,
     *                  nbf)
          if (vorb) then
          call matadd('work1','resao')
          else
          call matsimtr('work1','overlap','work2')
          call matadd('work2','resao')
          endif
          call CoefRead('tt',i,j,icoeffAO)
          call PairEnergyMP2(icoeffAO,iexchAO,iresAO,epairr,epair,
     *                       ncf,     equal,  nmo,   vorb,  sepairr,
     *                       sepair,nbf)
          call UpdateCoef(i,       j,   icoeffAO,iresAO,     bl(iepsi),
     *                    ncf,     nmo, nfirst,  npairs,     shift,
     *                    'diis',vorb,singles, ndiskdiisr, ndiskdiisc,
     *                    sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
          call CoefWrite(i,j,icoeffAO)
          call PairEnergyMP2(icoeffAO,iexchAO,iresAO,xepairr,xepair,
     *                       ncf,     equal,  nmo,   vorb,   sxepairr,
     *                       sxepair, nbf)
c
          iresult=1
        enddo        ! ij loop
        call fafCloses
        call sync_barrier
        call sync_barrier
        call para_recv_bcast(diis,TxCCSDInit)
        if (diis.eq.1.and.log_diis) then
        call dynamic_getmem(nvirt*nvirt,itmp_mat1)
        call dynamic_getmem(nvirt*nvirt,itmp_mat2)
       call sla_diis_matrix(af,singles,npairs,bl(ilist),bl(ipairimages),
     *                      ncf,nmo,bl(itmp_mat1),bl(itmp_mat2),
     *                      iMasID,mygid,mytid,nbf)
        call dynamic_retmem(2)
 345    continue
        call sync_barrier
        call para_recv_bcast(icode,TxDIISFlow)
        if (icode.eq.2) goto 346
        if (icode.eq.1) goto 345
        if (icode.ne.0) STOP 'Cannot happen, on slave DIIS'
        call sla_new_coeff(singles,af,idimen,nval,bl(ilist),
     *                     bl(ipairimages),npairs,iMasID,mygid,mytid)

 346    continue
        endif
        call sync_barrier
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
      call para_barrier
        call para_recv_bcast(iterate,TxCCSDIter)
      call para_barrier
        if (iterate.eq.0) exit
      enddo          ! iterative loop
      if (mp2.or.do_mp4) then
        goto 999
      endif
 126  continue
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(bl(ilist),npairs)
      call para_unpack_int(nstrong,1)
        if (info.ne.0) then
          ERROR
          call flush(6)
        endif
      call dynamic_matdef('WernerX','q',idimen,idimen)
      iX=mataddr('WernerX')
      call dynamic_matdef('sing_res','r',idimen,1)
      ising_res=mataddr('sing_res')
      if ((mp4.or.restart).and.l_triples) then
c     call sla_trip(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),
c    *          ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
c    *          singles.and.cc.and..not.ccsd,af,nslv,iMasID,mygid,mytid)
c     call sort_3ext_sla(nval,idimen,bl(irecadrtt),npairs,
c    *                   ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,
c    *                   af,iMasID,mygid,mytid,nslv)
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
c     call fafCloses
      call slave_3ext(nval,idimen,bl(irecadrtt),npairs,ndisktt,
     *               ncf,lbintt,thresh,byt8,nmo,vorb,af,
     *               iMasID,mygid,mytid,ichar_count,bl(ivtable),
     *               bl(im_table),bl(ivrevtable),bl(iotable),ext3_ratio)
      call fafCloses
      call make_3ext_pairs_sla(idimen,nval,af,ichar_count,bl(ivtable),
     *                         bl(im_table),bl(iotable),bl(ivrevtable),
     *                         iMasID,mygid,mytid)
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
      call new_tri_sla(nval,idimen,bl(iepsi+ncore),ccsd,qcisd,af,iMasID,
     *              mygid, mytid,bl(iotable),bl(ivtable),bl(ivrevtable),
     *              bl(im_table),bl(ijtable),ichar_count)
      if (restart) goto 999
      endif
c
c The main CCSD loop:
c
      iteration=1
 666  continue   ! for MP4 only
      do
      call secund(slave_cpu0)
      iteration=iteration+1
        call para_barrier
        call para_recv_bcastpack(TxCCSDIter1)
        call para_unpack_real(tot_xmax,1)
        call para_unpack_int(iterate,1)
        call para_barrier
        if (iterate.eq.0) exit
        if (iteration.eq.2) tot_xmax=0.1d0
c
      call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
c
c       call izeroit(bl(irecadre),npairs*nslv)
        int_kind='e'
        ipass=0
        if (.not.do_mp4) then
      call SLA_EEO_INT(ncs,    ncf,     bl(ictr),   nval,  nmo,
     *                 nfirst, nlast,   thresh,
     *                 vorb,   ndisktre, npairs,
     *                 iMasID, mygid,   mytid,  nslv, nbf)
        ndiske=-1
        lrece =-1
        lbine =-1
        endif
        if (ccsd) then
c       call izeroit(bl(irecadrtx),nval*nval*nslv)
        int_kind='tx'
        ipass=0
      call slave_Gen(ncs,     ncf,          ictr,    nval,   nmo,
     *               nfirst,  nlast,        thresh,  core,   xmaxdisk,
     *               ndisktx, bl(irecadrtx),lrectx,  lbintx, int_kind,
     *               nodisk,  byt8,         small,   vorb,   iMasID,
     *              mygid,   mytid,        ndisktrtx,tot_xmax,bl(ilist),
     *               nstrong,bl(isympairs),ij_unique,nbf)
c       call izeroit(bl(irecadrtc),nval*nval*nslv)
        int_kind='tc'
        ipass=0
      call slave_Gen(ncs,     ncf,          ictr,    nval,   nmo,
     *               nfirst,  nlast,        thresh,  core,   xmaxdisk,
     *               ndisktc, bl(irecadrtc),lrectc,  lbintc, int_kind,
     *               nodisk,  byt8,         small,   vorb,   iMasID,
     *               mygid,   mytid,       ndisktrtc,tot_xmax,bl(ilist),
     *               nstrong,bl(isympairs),ij_unique,nbf)
        endif
c
        do i=1,8
          elaps(i)=0.0d0
        enddo
        call fafCreates(ndiskalpha)
        call sync_barrier
        call para_recv_bcast(ndiskalpha,TxCCSDInit)
        if (cc.or.ccsd.or.do_mp4) then
        call CCalphaslave(ncf,nval,bl(irecadrx),npairs,ndiskx,lbinx,
     *                       thresh,norecalc,byt8,ccsd,nmo,vorb,elaps,
     *                       ndiskalpha,af,iMasID,mytid,nbf)
        call sync_barrier
        call para_reduce(elaps,5,TxCCSDRedu)
        endif
        call sync_barrier
        call zeroit(bl(ibeta),nval*nval)
        call para_recv_bcastreal(bl(ibeta),nval*nval,TxCCSDInit)
        call B41a_G_Slave(ncf,nval,norecalc,cc,ccsd,nfirst,bl(ifockMO),
     *                nmo,vorb,ndiskG41a,.true.,ndiskalpha,bl(ibeta),
     *                iMasID,mytid,nbf)
        call fafCloses
        call sync_barrier
c A matrix
      if (cc.or.ccsd.or.do_mp4) then
      call CCAOnceslav(bl(irecadrx),npairs,nval,    ndiskx,ncf,
     *                 lbinx,       thresh,norecalc,byt8,  ioverlap,
     *                 nmo,         vorb,  nslv,    iMasID,mygid,
     *                 mytid,       nbf)
        call sync_barrier
      endif
c
c
      if (singles.or.ccsd.or.do_mp4) then
      call Tl_gen_slav(ncf,nval,nmo,vorb,nslv,iMasID,mygid,mytid,nbf)
      call sync_barrier
      call EEO_vector_slave(nval,bl(irecadre),ndiske,lbine,bl(irecadrx),
     *                      ndiskx,lbinx,bl(irecadrc),ndiskc,lbinc,
     *                      ncf,npairs,thresh,nfirst,nlast,
     *                      byt8,nmo,vorb,nslv,iMasID,
     *                      mytid,bl(ilist),nbf)
        call sync_barrier
        call unlockS
        call unlockR
        call pointS(iiS)
        call para_recv_bcastreal(bl(iiS),idimen*nval,TxCCSDInit)
        if (cc.or.ccsd) then
          call Lt_gen_slav(idimen,nval,   bl(irecadrx),npairs,ndiskx,
     *                     ncf,   lbinx,  thresh,      byt8,  nmo,
     *                     vorb,  nslv,   iMasID,      mygid, mytid)
          call sync_barrier
        endif
        call pointRR(iiR)
        call para_recv_bcastreal(bl(iiR),idimen*nval,TxCCSDInit)
        call lockS
        call lockR
        endif
c  X matrix:
        call sync_barrier
        call para_recv_bcastpack(TxCCSDInit)
        call para_unpack_real(bl(iX),idimen*idimen)
        call para_unpack_int(ndiskg41a,1)
        call para_unpack_real(xiterenergy,1)
c
c  YZ matrices:
c
        do i=1,5
          elaps(i)=0.0d0
        enddo
        iresult=0
        call fafCloses()
        call fafCreates(ndiskYZ)
        call sync_barrier
        call para_recv_bcast(ndiskYZ,TxCCSDInit)
        do
          call para_initsend
          call para_pack_int(mytid,1)
          call para_pack_int(iresult,1)
          if (iresult.eq.1) then
            call para_pack_int(istart,1)
            call para_pack_int(istop,1)
            call para_pack_int(jstart,1)
            call para_pack_int(jstop,1)
          endif
          call para_send_pack(0,TxCCSDReq)
c
          call para_recv_pack(imastgid,TxCCSDJob)
          call para_unpack_int(iwork,1)
          if (iwork.eq.0) exit
          call para_unpack_int(istart,1)
          call para_unpack_int(istop,1)
          call para_unpack_int(jstart,1)
          call para_unpack_int(jstop,1)
c
c calculation, the kernel from master:
c
        call CCYZ(istart,istop,jstart,jstop,ncf,nval,bl(irecadrx),
     *        bl(irecadrc),npairs,ndiskx,
     *        ndiskc,lbinx,lbinc,bl(irecadrtx),ndisktx,lbintx,
     *        bl(irecadrtc),ndisktc,lbintc,thresh,
     *        byt8,ioverlap,cc,ccsd,bl(ilist),nmo,vorb,ndiskYZ,.true.,
     *        elaps,do_mp4,
     *              ichar_count,bl(iorevtable),bl(ivtable),bl(im_table),
     *        nbf)
c
          iresult=1
        enddo
c
c  YZ matrices STOP
        call sync_barrier
        call para_reduce(elaps,5,TxCCSDRedu)
        call para_barrier
c  Q parts:
c
        iresult=0
        call fafCloses()
        call fafCreates(ndiskQ)
        call sync_barrier
        call para_recv_bcast(ndiskQ,TxCCSDInit)
        do
          call para_initsend
          call para_pack_int(mytid,1)
          call para_pack_int(iresult,1)
          if (iresult.eq.1) then
            call para_pack_int(istart,1)
            call para_pack_int(istop,1)
            call para_pack_int(jstart,1)
            call para_pack_int(jstop,1)
          endif
          call para_send_pack(0,TxCCSDReqQ)
c
          call para_recv_pack(imastgid,TxCCSDJobQ)
          call para_unpack_int(iwork,1)
          if (iwork.eq.0) exit
          call para_unpack_int(istart,1)
          call para_unpack_int(istop,1)
          call para_unpack_int(jstart,1)
          call para_unpack_int(jstop,1)
c
c calculation, the kernel from master:
        call Qparts(istart,istop,jstart,jstop,ncf,nval,npairs,byt8,
     *              cc,ccsd,bl(ilist),nmo,vorb,ndiskYZ,af,ndiskQ,
     *              ichar_count,bl(iorevtable),bl(ivtable),bl(im_table),
     *              nbf)
c
c
          iresult=1
        enddo
        call sync_barrier
c
c  YZ matrices STOP
c
c   The main CCSD doubles loop:
c
        iresult_iter=0
        xmaxele=0.0d0
        do           ! ij loop
          call para_initsend
          call para_pack_int(mytid,1)
          call para_pack_int(iresult_iter,1)
          if (iresult_iter.eq.1) then
            call para_pack_int(i,1)
            call para_pack_int(j,1)
            call para_pack_int(ij,1)
            call para_pack_real(epair,1)
            call para_pack_real(epairr,1)
            call para_pack_real(xepair,1)
            call para_pack_real(sqij,1)
            call para_pack_real(xmaxele,1)
            call para_pack_real(xnormij,1)
          endif
          call para_send_pack(0,TxCCSDReq1)
          if (iresult_iter.eq.1) xmaxele=0.0d0
c
          call para_recv_pack(imastgid,TxCCSDJob1)
          call para_unpack_int(iwork,1)
          if (iwork.eq.0) exit
          call para_unpack_int(i,1)
          call para_unpack_int(j,1)
          call para_unpack_int(ij,1)
          if (cep2) then
            call para_unpack_real(xiterenergy,1)
          endif
c
c  Kernel from master:
          equal=.false.
          if (i.eq.j) equal=.true.
          if (int_array(bl(ilist),ij).eq.1) then
c
c  Read exchange operators
          if (.not.do_mp4) then
          call ExtrOne(i,    j,     bl(irecadrx),  npairs,     ndiskx,
     *                 ncf,  lbinx, thresh,        byt8,       'x',
     *                 'mo', nmo,   vorb,          bl(iexchAO))
c  Read external echange operators
          call ExtrOne(i,    j,     bl(irecadre),  npairs,     ndiske,
     *                 ncf,  lbine, thresh,        byt8,       'e',
     *                 'mo', nmo,   vorb,          bl(iwork1))
c  Add it do residuum matrix
          call matcopy('exchao','resao')
          call matadd('work1','resao')
          else
          call matzero('resao')
          endif
c
          call secund(s_ext_matstop)
          call elapsec(e_ext_matstop)
c
          e_ext_mat_tot1=e_ext_mat_tot1+e_ext_matstop-e_ext_matstart
          s_ext_mat_tot1=s_ext_mat_tot1+s_ext_matstop-s_ext_matstart
c  Qgen, iwork1=Qij, iwork2=Qji(T), probably OK, order of en. OK
          call elapsec(eqstart)
          call secund(sqstart)
c
          call Q_read_build(i,j,iX,ndiskQ,af,vorb,ncf,nmo,nval,
     *                      iwork1,iwork2,nbf)
c
          call elapsec(eqstop)
          call secund(sqstop)
c
          eq_tot1=eq_tot1+eqstop-eqstart
          sq_tot1=sq_tot1+sqstop-sqstart
c   S*Qij and add to residuum
          call elapsec(e_ext_matstart)
          call secund(s_ext_matstart)
c
          if (vorb) then
          call matadd('work1','resao')
          call matadd('work2','resao')
          else
          call matmmult('overlap','work1','work3')
          call matadd('work3','resao')
          call matmmult('work2','overlap','work3')
          call matadd('work3','resao')
          endif
c   Qji(T)*S and add to residuum
c
          call secund(s_ext_matstop)
          call elapsec(e_ext_matstop)
c
          e_ext_mat_tot1=e_ext_mat_tot1+e_ext_matstop-e_ext_matstart
          s_ext_mat_tot1=s_ext_mat_tot1+s_ext_matstop-s_ext_matstart
c   Read Gij and Gji(T) with part of CCSD
          call elapsec(egstart)
          call secund(sgstart)
c
          call ReadCalcG_41a(i,j,norecalc,ncf,nmo,vorb,ndiskG41a,
     *                       iwork1,.true.,nbf)
c
          call elapsec(egstop)
          call secund(sgstop)
c
          eg_tot1=eg_tot1+egstop-egstart
          sg_tot1=sg_tot1+sgstop-sgstart
c
          if (cc.or.ccsd.or.do_mp4) then
            continue
          else
            call elapsec(e_ext_matstart)
            call secund(s_ext_matstart)
c
            call CoefRead('tt',i,j,icoeffAO)
            call matconn('coeffAO','q',idimen,idimen,icoeffAO)
            call matadd1('coeffAO',-xiterenergy,'work1')
            call matdisc('coeffAO')
          endif
c       S*Sum*S
          if (vorb) then
          call matadd('work1','resao')
          else
          call matsimtr('work1','overlap','work2')
          call matadd('work2','resao')
          endif
c
          call secund(s_ext_matstop)
          call elapsec(e_ext_matstop)
c
          e_ext_mat_tot1=e_ext_mat_tot1+e_ext_matstop-e_ext_matstart
          s_ext_mat_tot1=s_ext_mat_tot1+s_ext_matstop-s_ext_matstart
c
          call elapsec(eenstart)
          call secund(senstart)
c
c Singles part of CISD CCSD start
          call secund(t0)
          if (singles.or.ccsd) then
          call dynamic_matdef('CISD','q',idimen,idimen)
          iresult=mataddr('CISD')
          call singles_CISD(i,j,ncf,nval,ioverlap,bl(irecadre),npairs,
     *               nfirst,nlast,ndiske,lbine,thresh,byt8,ccsd,iresult,
     *               nmo,vorb,nbf)
          call matadd('CISD','resao')
          call dynamic_matrem('CISD')
          endif
          call secund(t1)
          CPU_doubles_singles=CPU_doubles_singles+t1-t0
c Singles part of CISD CCSD stop
          else ! if (int_array(bl(ilist),ij).eq.1) then
            call matzero('resao')
          call ExtrOne(i,    j,     bl(irecadrx),  npairs,     ndiskx,
     *                 ncf,  lbinx, thresh,        byt8,       'x',
     *                 'mo', nmo,   vorb,          bl(iexchAO))
          endif
          call CoefRead('tc',i,j,icoeffAO)
          if (mp3.or.do_mp4) then
          call PairEnergy(icoeffAO,iresAO,iresAO,
     *             epairr,epair,ncf,equal,nmo,vorb,
     *             nbf)
          sqij=0.0d0
          if (mp3.and.mp4) call fafWrite(ndisk_mp3,bl(iresAO),8,
     *                                   idimen*idimen,1,ij,istatus)
          else
          call PairEnergy(icoeffAO,iexchAO,iresAO,
     *             epairr,epair,ncf,equal,nmo,vorb,
     *             nbf)
          endif
          if (omit.and.int_array(bl(ilist),ij).ne.1) epairr=epair
c
          call elapsec(eenstop)
          call secund(senstop)
c
          een_tot1=een_tot1+eenstop-eenstart
          sen_tot1=sen_tot1+senstop-senstart
c
          call elapsec(ecostart)
          call secund(scostart)
c
          call CoefRead('tt',i,j,icoeffAO)
          if (.not.(mp3.or.do_mp4)) then
          call CINorm(icoeffAO,xnormij,ncf,equal,nmo,vorb,nbf)
          if (.not.omit.or.int_array(bl(ilist),ij).eq.1) then
          call UpdateCoef(i,     j,     icoeffAO,iresAO,     bl(iepsi),
     *                    ncf,   nmo,   nfirst,  npairs,     shift,
     *                    'diis',vorb,  singles, ndiskdiisr, ndiskdiisc,
     *                    sqij,  xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
          else
            sqij=0.0d0
            xmaxele=0.0d0
          endif
          endif
c
          call elapsec(ecostop)
          call secund(scostop)
c
          eco_tot1=eco_tot1+ecostop-ecostart
          sco_tot1=sco_tot1+scostop-scostart
c
          call CoefWrite(i,j,icoeffAO) ! Still I write "updated"
c          amplitudes in CoefWrite with "omit" switched on, because 
c          of non Gauss-Seidel algorithm
c
          call elapsec(eenstart)
          call secund(senstart)
c
          call CoefRead('ee',i,j,icoeffAO)
          call PairEnergy(icoeffAO,iexchAO,iresAO,
     *            xepairr,xepair,ncf,equal,nmo,vorb,
     *             nbf)
          if (omit.and.int_array(bl(ilist),ij).ne.1) xepairr=xepair
c
          call elapsec(eenstop)
          call secund(senstop)
c
          een_tot1=een_tot1+eenstop-eenstart
          sen_tot1=sen_tot1+senstop-senstart
c
          iresult_iter=1
        enddo ! of the loop over i,j pairs in main CCSD loop
      call sync_barrier
      if (singles.or.ccsd.or.do_mp4) then
      do
      call sync_barrier
      call para_recv_bcast(iwork,TxCCSDInit)
      if (iwork.eq.0) exit
      call Tr_gen_slav(idimen, nval,  ncf,  vorb, nslv,
     *                       iMasID, mygid, mytid)
      enddo
      endif
      call sync_barrier
c
      if (.not.(mp3.or.do_mp4)) then
        call sync_barrier
        call para_recv_bcast(diis,TxCCSDInit)
        if (diis.eq.1.and.log_diis) then
        call dynamic_getmem(nvirt*nvirt,itmp_mat1)
        call dynamic_getmem(nvirt*nvirt,itmp_mat2)
       call sla_diis_matrix(af,singles,npairs,bl(ilist),bl(ipairimages),
     *                      ncf,nmo,bl(itmp_mat1),bl(itmp_mat2),
     *                      iMasID,mygid,mytid,nbf)
        call dynamic_retmem(2)
 745    continue
        call sync_barrier
        call para_recv_bcast(icode,TxDIISFlow)
        if (icode.eq.2) goto 746
        if (icode.eq.1) goto 745
        if (icode.ne.0) STOP 'Cannot happen, on slave DIIS'
        call sla_new_coeff(singles,af,idimen,nval,bl(ilist),
     *                     bl(ipairimages),npairs,iMasID,mygid,mytid)

 746  continue
        endif
        endif
c
      call sync_barrier
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
        call initsingles_slave(nval,ncf,singles.or.mp4,nmo,vorb,iMasID,
     *                       nbf)
      call secund(slave_cpu1)
      slave_cpu=slave_cpu1-slave_cpu0
      call sync_barrier
      call para_reduce(slave_cpu,1,TxCCSDRedu)
      call para_barrier
c
      call sync_barrier
      enddo ! of the main ccsd loop
      if (do_mp4) then 
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
      goto 222
      endif
      if (mp4.and..not.do_mp4) then
        mp3=.false.
        do_mp4=.true.
        goto 666
      endif
c
      call dynamic_retmark ! Mark after initR
      call matremark
      call retmark
c
      if ((cc.and.singles).and.l_triples) then
c     call sla_trip(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),
c    *          ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
c    *          singles.and.cc.and..not.ccsd,af,nslv,iMasID,mygid,mytid)
c     call sort_3ext_sla(nval,idimen,bl(irecadrtt),npairs,
c    *                   ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,
c    *                   af,iMasID,mygid,mytid,nslv)
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
c     call fafCloses
      call slave_3ext(nval,idimen,bl(irecadrtt),npairs,ndisktt,
     *               ncf,lbintt,thresh,byt8,nmo,vorb,af,
     *               iMasID,mygid,mytid,ichar_count,bl(ivtable),
     *               bl(im_table),bl(ivrevtable),bl(iotable),ext3_ratio)
      call fafCloses
      call make_3ext_pairs_sla(idimen,nval,af,ichar_count,bl(ivtable),
     *                         bl(im_table),bl(iotable),bl(ivrevtable),
     *                         iMasID,mygid,mytid)
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
      call new_tri_sla(nval,idimen,bl(iepsi+ncore),ccsd,qcisd,af,iMasID,
     *              mygid, mytid,bl(iotable),bl(ivtable),bl(ivrevtable),
     *              bl(im_table),bl(ijtable),ichar_count)
      endif
 999  continue
      call secund(slave_tot_cpu1)
c -- send timings back to master then move on
      slave_tot_cpu=slave_tot_cpu1-slave_tot_cpu0
      call para_initsend
      call para_pack_real(slave_tot_cpu,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
      call fafTerminate(info)
      call allclose_and_delete(.true.)
c Release memory:
      call matremark
      call retmark
      call dynamic_retmark
      end
!=======================================================================
      subroutine slave_Gen(ncs,     ncf,     ictr,    nval,  nmo,
     *                     nfirst,  nlast,   thresh,  core,  xmaxdisk,
     *                     ndisk,   irecadr, lrec,    lbin,  int_kind,
     *                     nodisk,  byt8,    small,   vorb,  iMasID,
     *                     mygid,   mytid,   ndisktr,tot_xmax,list,
     *                     nstrong,isympairs,ij_unique,nbf)
c This subroutine calculates all necessary integrals needed for CCSD and
c transforms them.  The order of explanation is different that order of 
c parameters
c INPUT:
c ncs      - number of contracted shells
c ncf      - number of contracted basis functions
c ictr     - pointer do contraction info etc.
c nval     - number of correlating occupied orbitals.
c nmo      - number of occupied orbitals.
c nfirst   - first correlating orbital.
c nlast    - last  correlating orbital.
c thresh   - integral threshold
c core     - this parameter is used for determination of number of
c            correlating orbitals. Orbitals with energy below core are frozen
c xmaxdisk - not used at present, but it is kept for future use
c int_kind - it is character string indicating the kind of calculated
c            integrals. May be: 'c' coulomb, 'x' exchange, 'e' EEO, 
c            'tc' - three external coulomb, 'tx' - three external exchange
c nodisk   - the parameter indicating if the calculation is performed
c            with all amplitudes kept in RAM memory. It means that
c            symmetrized and antisymmetrized amplitudes will be kept in RAM in
c            addition to usual amplitudes (EEO only)
c byt8     - logical parameter, if .true. all bins are writen to disk in 
c            real*8 precision, otherwise in 5 byte compact format
c small    - how to perform calculations, using symm anti-method or
c            using full matrices of K integrals and amplitudes? (EEO only)
c            usually=.true., false was used for experimantal purposes,
c            may be useful for very lage molecules with sparsity
c vorb     - the logical indicating if intergals will transformed
c            from AO to MO - virt. orbitals (second transformation)
c iMasID   - the master TID
c mygid    - this process slave group ID
c mytid    - this process slave TID
c OUTPUT:
c irecadr  - irecadr contains the number of records for the ij bin,
c            written on each slave
c            where ij is a pair number. This is in fact two-dimensional
c            table, irecadr(npairs,nslv) where nslv is number of slaves
c lrec     - the length of bin record
c lbin     - the length of bin
c ndisktr  - contains integrals extracted from bins, if vorb=.true. the
c            integrals are already transformed, otherwise there are
c            rewritten only  from integrals in the bins. But integrals
c            in nodisk still are needed if vorb=.true. because some
c            parts of program requires integrals in AO form.
c            This should be changed, if vorb=.false. the bins should be
c            deleted because are not used at all.
c ndisk    - number of disk unit where bins are stored. This in really
c            output value
c tot_xmax - maximum residuum element in previous iteration, for
c            screening
      use memory
      use newpara
      use kinds
      use ccounters
      implicit real*8 (a-h,o-z)
      integer*1 i01
      logical LastSymPair,emp2only,nodisk,keepnodisk,byt8,small,vorb
      logical memory_success,af,split
      dimension xintxx(9)
      common /job/jobname,lenJ
      parameter (max_merg=(MAX_MERGE))  ! fragile!!
      common /indexsaver/ik,ikindex(2,max_merg) ! for EEOGenerator
      logical reset_s_a,cache,oldcache,integerT
      common /SymmAntiCache/reset_s_a
      common /GlobalCCSD/ cache,integerT
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*3 ch3
      character scftype*11
      character*256 jobname,scrfile,filename
      logical nofr,exst,restrt,dualbasis,smal,reset
      character*(*) int_kind
      dimension xnmo(2)
ctj header end
      dimension irecadr(*)
      character*256 filname1
      dimension Ipass(2,28), Jpass(2,28) !28 is 28 comp.of cartisian I-function
      parameter (i_orb_size=28)
      integer igran2pair,ipair2gran
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
ctest
      common /zero_in_EEO/ zeroed,zeroed_max
      character*7 int_kinder
      common /switch_off/ int_kinder
      parameter (ixint_time=20)
      real*8 xint_time(ixint_time)
      equivalence (tott,xint_time(ixint_time))
      equivalence (tot, xint_time(ixint_time-1))
      equivalence (skip,xint_time(ixint_time-2))
      character(len=100) :: err_msg
      integer(kind=i_4) :: inf
      integer list(*),nstrong
      integer, dimension(nval*nval) :: list_copy
      integer :: isympairs(*),ij_unique
      integer irec_no,irec_last
      irec_no=-1
      irec_last=-1
      int_kinder='xxxxxxx'
      xmax_x_int=0.0d0
      zeroed=0.0d0
      zeroed_max=0.0d0
      npairs=(nval**2+nval)/2
      if (int_kind.eq.'e') then
        icountere=icountere+1
        if (icountere.eq.1) then ! build list_copy with all strong pairs
          do i=1,npairs
            list_copy(i)=1
          enddo
          nstrong_copy=npairs
        else
          do i=1,npairs
            list_copy(i)=list(i)
          enddo
          nstrong_copy=nstrong
        endif
      endif
      if (int_kind.eq.'c'.or.int_kind.eq.'tc') int_kinder='coulomb'
c
ctj Int_kind=c coulomb, Int_kind=x exchange, Int_kind=e EEO
ctj -- BEGIN get some common variables
c     call getival('nslv',nslv)
      call secund(total_cpu0)
      keepnodisk=nodisk
      call getival('iout',iout)
      iprnt=0
      smal=.TRUE.
      strace=0.0d0
      sconstr=0.0d0
      sort=0.0d0
      extr=0.0d0
      const=0.0d0
      c1=0.0d0
      c2=0.0d0
      etrace=0.0d0
      econstr=0.0d0
      eort=0.0d0
      eetr=0.0d0
      eonst=0.0d0
      e1=0.0d0
      e2=0.0d0
ctj -- END   get some common variables
  46  format(/72('='))
      if (Int_kind.ne.'c'.and.Int_kind.ne.'x'.and.Int_kind.ne.'e'.and.
     *  Int_kind.ne.'tx'.and.Int_kind.ne.'tc'.and.Int_kind.ne.'tt') then
        call nerror(1,'GenCoulExInt',
     1  'Int_kind does not indicate any kind of integrals',0,0)
      endif
      write(iout,*)' '
c-----------------------------------------------------------
      ncs_sm=ncs
      ncf_sm=ncf
c-----------------------------------------------------------
c This is for use in EEO Generator and for bins construction when we
c decide how much memory we need to reserve
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ncf*(ncf-1)/2
      nvirt=nbf-nmo
c  put down a memory marker
c
      call mmark
      call matmark
      call dynamic_mmark 
c-----------------------------------------------------------
c zero out irrelevant dft stuff
      nfock=1
      idft=0
      ax=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
c-------------------------------------------------------
c allocate memory for an array mapping contr.func. to contr.shells :
c SLAVE!
#if 1
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
c
c get memory for screening density
c
      np4=4
      call dynamic_matdef('dsmx','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      idics=mataddr('dsmx')
      call dynamic_getmem(ncs*ncs,ischwarz)
      call schwarz_calc(bl(ictr),ncs,bl(ischwarz))
c-------------------------------------------------------
c initialize the two-el. integral program
c
ctmp  thint=thresh
ctmp  iforwhat=5
ctmp  call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
ctmp *             ax,     nrad,   nang,   lrad,   lang,
ctmp *             Iradq,  NBatch, .true., nmo,    0,
ctmp *             scftype,xintxx, nblocks,maxbuffer,maxlabels)
c                     ^       ^       ^       ^          ^
c OUTPUT:          NOT USED PRICES NO. BL. INT. BUF.  LABELS
c                          (ARRAY)           SIZE      SIZE
c-------------------------------------------------------
c
#endif
      call mmark
      call dynamic_mmark
c---------------------------------------------------------------------
c WARNING! Screening is turned off at present!
c SLAVE!
#if 1
      call  fill1(ncs,bl(idics))
c
c---------------------------------------------------------------------
c
      icano=mataddr('cano')
c.................................................
c check if there will be a split in integral calculations
c
      call check_sizes(bl(ictr),ncs,ncf,nval)
c.................................................
      if(ncore.gt.0) call matsub('occa','cano',1,nmo)
      call  matsub('occu','cano',nfirst,nlast)
      call  matsub('virt','cano',nmo+1,nbf)
      ioccu=mataddr('occu')
c.................................................
c   reserve space for one AO integral matrix
      call dynamic_matdef('xmat','q',ncf,ncf)
      ixadr=mataddr('xmat')
c space for one half-transformed coulomb or exchange operator
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
        call dynamic_matdef('halftra','q',nval,nval)
        ihalftra=mataddr('halftra')
      else if (int_kind.eq.'tt') then
        call dynamic_matdef('halftra','r',nvirt,nval)
        ihalftra=mataddr('halftra')
      else
        call dynamic_matdef('halfmulam','q',nval,nval)
        ihalfmulam=mataddr('halfmulam')
        call dynamic_matdef('halflammu','q',nval,nval)
        ihalflammu=mataddr('halflammu')
      endif
c  put down a memory marker
#endif
      call mmark
      call dynamic_mmark
c
c SLAVE!
      thint=thresh
      istat=igetival('istat')
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
c  nrec is the total number of records written on all files
c  irec is the current counter in the file
      nrec=0
      irec=0
c---------------------------------------------------
c
      call secund(tt)
      call secund(tt3)
      call elapsec(telap3)
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c  ijtot is the number of ics,jcs pairs really calculated & transformed
c   if an (ics,jcs) pair is skipped because there are no integrals in it
c   then it is NOT incremented
      ijtot=0
c  Count the actual number of (mu, lam) pairs, both the diagonal
c   and the non-diagonal ones
      mulam=0
      mulamd=0
      ENonZero=zero
ckw..........
c SLAVE!
#if 1
      call secund(txxx1)
      call init_info('fock')
ckw..........
c
c turn off symmetry for mp2 integrals :
      call symmoff
c
c  bl(icol) and bl(jcol) store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
      call getint(ncf,irow)
      call getint(ncf,icol)
      call getint(ncf,irow1)
      call getint(ncf,icol1)
      call getint(ncf,lzero)
      call matrix_int_zeroit(bl(lzero),ncf)
C
CTJ
#endif
#if 1
c SLAVE!
#if 1
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      call dynamic_getmem(npairs/intsize+1,icountr) !memory for bin counters
      call izeroit(bl(icountr),npairs)         !zero out the bin counters
      ioccvirt=npairs ! number of bins
      else if (int_kind.eq.'tt') then
      ioccvirt=nval*nvirt ! number of bins
      call dynamic_getmem(nval*nvirt/intsize,icountr) !memory for bin counters
      call izeroit(bl(icountr),nval*nvirt)         !zero out the bin counters
      else
      call dynamic_getmem(nval*nval/intsize+1,icountr) !memory for bin counters
      call izeroit(bl(icountr),nval*nval)         !zero out the bin counters
      endif
#endif
c
c Tell me master please how long bin should I use?
c
      call para_recv_bcast(lbin,TxCCSDJob)
c
c Here are memory allcations - slave!
c
#if 1
      if (Int_kind.eq.'e'.or.Int_kind.eq.'tt') then
      call dynamic_getmem((ioccvirt/i4size+1)*2,igran2pair)
      call dynamic_getmem((ioccvirt/i4size+1)*1,ipair2gran)
      endif
      call dynamic_mmark
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      indexmem=npairs*lbin/2+1
      if (byt8) then
      call dynamic_getmem(npairs*lbin,ibin8) !memory for bins 
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices 
      call zeroit(bl(ibin8),npairs*lbin)
      call izeroit(bl(indexadr),intsize*indexmem)        
      ibin4adr=1
      ibin1adr=1
c
      else
      int4mem=npairs*lbin/i4size+1
      int1mem=npairs*lbin/i1size+1
      call dynamic_getmem(int4mem,ibin4adr)  !memory for bins 
      call dynamic_getmem(int1mem,ibin1adr)  !memory for overflow 
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices 
      call zeroit(bl(ibin4adr),int4mem)
      call zeroit(bl(ibin1adr),int1mem)
      call izeroit(bl(indexadr),intsize*indexmem)
      ibin8=1
      endif
      call flush(6)
      else if (int_kind.eq.'tt') then
      indexmem=nval*(nbf-nmo)*lbin/2+1
      if (byt8) then
      call dynamic_getmem(nval*(nbf-nmo)*lbin,ibin8) !memory for bins 
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices 
      call zeroit(bl(ibin8),nval*(nbf-nmo)*lbin)
      call izeroit(bl(indexadr),intsize*indexmem)        
      ibin4adr=1
      ibin1adr=1
c
      else
      int4mem=nval*(nbf-nmo)*lbin/i4size+1
      int1mem=nval*(nbf-nmo)*lbin/i1size+1
      call dynamic_getmem(int4mem,ibin4adr)  !memory for bins 
      call dynamic_getmem(int1mem,ibin1adr)  !memory for overflow 
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices 
      call zeroit(bl(ibin4adr),int4mem)
      call zeroit(bl(ibin1adr),int1mem)
      call izeroit(bl(indexadr),intsize*indexmem)
      ibin8=1
      endif
      call flush(6)
      else !if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      indexmem=nval*nval*lbin/2+1
      if (byt8) then
      call dynamic_getmem(nval*nval*lbin,ibin8) !memory for bins 
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices 
      call zeroit(bl(ibin8),nval*nval*lbin)
      call izeroit(bl(indexadr),intsize*indexmem)        
      ibin4adr=1
      ibin1adr=1
c
      else
      int4mem=nval*nval*lbin/i4size+1
      int1mem=nval*nval*lbin/i1size+1
      call dynamic_getmem(int4mem,ibin4adr)  !memory for bins 
      call dynamic_getmem(int1mem,ibin1adr)  !memory for overflow 
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices 
      call zeroit(bl(ibin4adr),int4mem)
      call zeroit(bl(ibin1adr),int1mem)
      call izeroit(bl(indexadr),intsize*indexmem)
      ibin8=1
      endif
      endif !if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
#endif

c .................................................
c Tell me master please the number of ndisk
      call fafCloses
      call fafCreates(ndisk)
      call sync_barrier
      call para_recv_bcastpack(TxCCSDJob)
      call para_unpack_int(ndisk,1)
      if (Int_kind.eq.'e'.or.Int_kind.eq.'tt') then
        call para_unpack_int(igranules,1)
        call para_unpack_int(igranulesize,1)
        call para_unpack_int(ioccvirt,1)
        call para_unpack_int4(bl(igran2pair),2*ioccvirt)
        call para_unpack_int4(bl(ipair2gran),1*ioccvirt)
        call dynamic_lock(bl(igran2pair),ii)
        call dynamic_lock(bl(ipair2gran),ii)
      endif
#endif
CTJ
c--------------------------------------------------------------------
      icont=0
c--------------------------------------------------------------------
c 2*merged_size*npairs = matrix with results in EEO
c (ipass_size=30)*ncf*ncf
c
      tot_latency1=0.0d0
      tot_latency2=0.0d0
      if (Int_kind.eq.'e') then
        max_merg_shell_size=keepmax
        if (small) then
        call dynamic_matdef('screen','q',ncf,ncf)
        iscr=mataddr('screen') ! maybe xmat should be used (savings)
        call matzero('screen')
c
        call GenAmplSymm_Sla(npairs, ncf,       nval,  nodisk,   nfirst,
     *                       nmo,    vorb,     .true., 'slave', ndisks,
     *                       ndiska, bl(iscr), list,isympairs,ij_unique,
     *                       nslv,iMasID,mygid,mytid)
c                            ! Warning list is dummy! really!???? TJ
c Master give me the data!
c
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(max_merg_shell_size,1)
        ik_capacity=max_merg_shell_size
c       call para_unpack_int(ndisks,1)
c       call para_unpack_int(ndiska,1)
        call para_unpack_int(split,1)
        call para_unpack_real(bl(iscr),ncf*ncf)
        call para_unpack_real(bl(idics),ncs*ncs)
c tmp, because for symm-anti scr idics is not correct
        call  fill1(ncs,bl(idics))
c tmp END
        call dynamic_lock(bl(iscr),info)
c
        reset_s_a=.true.
        oldcache=cache
        cache=.false.
        call dynamic_getmem(ntrisymm*max_merg_shell_size,ibasesymm) 
        call dynamic_getmem(ntrianti*max_merg_shell_size,ibaseanti) 
        write(iout,*) 'max_merg_shell_size=  ',max_merg_shell_size
        if (.not.nodisk)
     *write(iout,*)'amplit. will be read',ncf*ncf/2/max_merg_shell_size,
     *         '  times'
        call flush(6)
      else ! if (small) then - no shell merging
        continue
      endif ! if (small) then
      endif ! if (Int_kind.eq.'e') then
c
      call dynamic_unlock(irecadr,iii)
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
        do iii=1,npairs*nslv
        irecadr(iii)=0
        enddo
      else if (int_kind.eq.'tt') then
        do iii=1,nval*(nbf-nmo)*nslv
        irecadr(iii)=0
        enddo
      else
        do iii=1,nval*nval*nslv
        irecadr(iii)=0
        enddo
      endif
C
      ncf2=ncf*ncf
      nskipped=0
      iare_zeros=0
      if (Int_kind.eq.'tx'.or.Int_kind.eq.'tc') then
        irecidx=nval*nval
      else if (Int_kind.eq.'e'.or.Int_kind.eq.'tt') then
        irecidx=igranules
      else
        irecidx=npairs
      endif
      do i=1,ixint_time
        xint_time(i)=0.0d0
      enddo
c     tot=0.0d0
c     skip=0.0d0
c     tott=0.0d0
c
      do ! ICS, JCS LOOP =================================================
c
c Send Master ready signal:
      iready=1
      call para_initsend
      call para_pack_int(iready,1)
      call para_pack_int(mytid,1)
      call para_send_pack(0,TxCCSDReq)
c
c iMasID,mygid,mytid
c
      call elapsec(xlatency0)
      call para_recv_pack(imastgid,TxCCSDJob)
      call para_unpack_int(iwork,1)
      call elapsec(xlatency1)
      tot_latency1=tot_latency1+xlatency1-xlatency0
      if (iwork.eq.0) then
#if 1
      if (int_kind.eq.'c'.or.int_kind.eq.'x'.or.int_kind.eq.'e') then
        if (int_kind.eq.'e') then
        reset=.true.
        call secund(ctrstart)
        call elapsec(etrstart)
        call EEOGenerator (ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),lmp2int,bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), nrow,
     3                    ncol,
     4                             lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .true.,reset,mygid,         ibasesymm,ibaseanti,
     8               max_merg_shell_size,nodisk,byt8,small,ndiskresid,
     9          max_amplitudes_size,ishellsaddr,iamplitaddr,iresultaddr,
     *       igranules,igranulesize,bl(igran2pair),bl(ipair2gran),split,
     1          bl(iscr),bl(mapf2s),ik_capacity,tot_xmax,list_copy,
     2          nstrong_copy,isympairs,ij_unique)
        call secund(ctrstop)
        call elapsec(etrstop)
        ttrans=ttrans+ctrstop-ctrstart
        elaptrans=elaptrans+etrstop-etrstart
        call WriteLeftOversGr(ndisk,npairs,lbin,bl(icountr),irecadr,
     1                      bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                      bl(ibin8),byt8,.true.,mygid,
     3          igranules,igranulesize,bl(igran2pair),bl(ipair2gran),
     4          irec_no,irec_last)
        else
      call WriteLeftOvers(ndisk,npairs,lbin,bl(icountr),irecadr,
     1                    bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                    bl(ibin8),byt8,.true.,mygid,irec_no,irec_last)
        endif
      else if (int_kind.eq.'tt') then
      call ext3_WriteLeftOvers(ndisk,npairs,lbin,bl(icountr),irecadr,
     1                    bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                    bl(ibin8),byt8,.true.,mygid,nval,ncf,nmo,
     3                    igranules,igranulesize,bl(igran2pair),
     4                    bl(ipair2gran),irec_no,irec_last,nbf)
      else
      call WriteLeftOversExt(ndisk,npairs,nval,lbin,bl(icountr),irecadr,
     1                    bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                    bl(ibin8),byt8,.true.,mygid,irec_no,irec_last)
      endif
#endif
c exit: no send data here
        exit
      endif
      call para_unpack_int(ics,1)
      call para_unpack_int(jcs,1)
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf*ncf*ics_size
#if 1
c inside if - slave!
c
           call get_shell_size(bl(ictr),jcs,jcs_size)
c
c
           lmp2_size=lmp2_siz1*jcs_size
c
c check if a size of the lmp2 integral buffer is not too big
c if it is  then split over kcs ( and possibly ics)
c
          call check_size2(lmp2_size,ncf,nval,max_merg_shell_size,
     *                     npairs,int_kind,ntimes,nodisk,split)
          ipassess=ipassess+ntimes
c
          call set_passes(bl(ictr),ics,jcs,ics_size,jcs_size,
     *                      ntimes,Ipass,Jpass,Itimes,Jtimes)
c
          do itime=1,itimes
             icf1=ipass(1,itime)
             icf2=ipass(2,itime)
             iatonce=icf2-icf1+1
          do jtime=1,jtimes
             jcf1=Jpass(1,jtime)
             jcf2=Jpass(2,jtime)
             jatonce=jcf2-jcf1+1
c
             lmp2_size=iatonce*jatonce*ncf2
             if (.not.split.and.int_kind.eq.'e'.and.
     *                 (iatonce*jatonce).gt.ik_capacity) then
               if (ik.eq.0) STOP 'Too little mem in EEO'
               reset=.true.
             print *,'* * * * * iatonce*jatonce,ik,ik_capacity flushed',
     *               iatonce*jatonce,ik,ik_capacity
               call secund(ctrstart)
               call elapsec(etrstart)
        call EEOGenerator (ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),lmp2int,bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), nrow,
     3                    ncol,
     4                             lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .true.,reset,mygid,         ibasesymm,ibaseanti,
     8               max_merg_shell_size,nodisk,byt8,small,ndiskresid,
     9          max_amplitudes_size,ishellsaddr,iamplitaddr,iresultaddr,
     *       igranules,igranulesize,bl(igran2pair),bl(ipair2gran),split,
     1          bl(iscr),bl(mapf2s),ik_capacity,tot_xmax,list_copy,
     2          nstrong_copy,isympairs,ij_unique)
               call secund(ctrstop)
               call elapsec(etrstop)
               ttrans=ttrans+ctrstop-ctrstart
               elaptrans=elaptrans+etrstop-etrstart
               reset=.false.
             endif
             if (.not.split.and.int_kind.eq.'e'.and.
     *                 (iatonce*jatonce).gt.ik_capacity) then
             print *, 'Program error, wrong memory assesment on master'
             STOP 'Program error, wrong memory assesment on master'
             endif
c
            call secund(tt2)
            call elapsec(telap2)
            call dynamic_getmem(lmp2_size,lmp2int) 
            if (Int_kind.eq.'e'.or.Int_kind.eq.'x') then
              do itj=0,lmp2_size-1
                bl(lmp2int+itj)=0.0d0
              enddo
            endif
            call mmark
            call dynamic_mmark
            if (Int_kind.eq.'c'.or.Int_kind.eq.'tc'.or.Int_kind.eq.'tt')
     *      then
c           call int_coul(bl, bl(ictr), thresh, ics, icf1,
c    1                    icf2, jcs, jcf1, jcf2, bl(mapf2s),
c    2                    bl(idics),iprnt,bl(lmp2int),nintotal,nrow,
c    3                    ncol, bl(irow),bl(icol),bl(lzero))
            call new_coul_int(bl(ictr), thresh,  ics,     icf1,   icf2,
     *                       jcs,  jcf1,   jcf2,  bl(mapf2s), bl(idics),
     *                        iprnt,bl(lmp2int), nintotal, nrow,   ncol,
     *                        bl(irow),  bl(icol), bl(lzero), ncf, ncs,
     *                        bl(ischwarz),skip,tot,xint_time)
            else if (Int_kind.eq.'x' .or.  Int_kind.eq.'e'
     *                                      .or.Int_kind.eq.'tx') then
c           call int_lmp2(bl,bl(ictr),thresh,
c    1                    ics,icf1,icf2,jcs,jcf1,jcf2,
c    2                 bl(mapf2s),bl(idics),iprnt,bl(lmp2int),nintotal,
c    3                    nrow,   ncol,    bl(irow),bl(icol),bl(lzero))
            call new_exch_int(bl(ictr), thresh,  ics,     icf1,   icf2,
     *                       jcs,  jcf1,   jcf2,  bl(mapf2s), bl(idics),
     *                        iprnt,bl(lmp2int), nintotal, nrow,   ncol,
     *                        bl(irow),  bl(icol), bl(lzero), ncf, ncs,
     *                        bl(ischwarz),skip,tot,xint_time)
            if (Int_kind.eq.'x') then
              do itj=0,lmp2_size-1
                if (dabs(bl(lmp2int+itj)).gt.xmax_x_int)
     *              xmax_x_int=dabs(bl(lmp2int+itj))
              enddo
            endif
            endif
            call retmark
            call dynamic_retmark
          if(nintotal.eq.0) then ! ought to calculate splitted shell again?
            nskipped=nskipped+1
            iare_zeros=iare_zeros+lmp2_size
            call dynamic_retmem(1) 
            call secund(tt3)
            call elapsec(telap3)
            tinteg=tinteg+tt3-tt2
            elapint=elapint+telap3-telap2
            cycle
          else
            ijtot=ijtot+1
            tott=tott+nintotal
          endif
      call secund(tt3)
      call elapsec(telap3)
      tinteg=tinteg+tt3-tt2
      elapint=elapint+telap3-telap2
c......................................................................
      if (Int_kind.eq.'e') then
       reset=.false.
       call EEOGenerator (ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),lmp2int,bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), nrow,
     3                    ncol,
     4                             lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .true.,reset,mygid,         ibasesymm,ibaseanti,
     8               max_merg_shell_size,nodisk,byt8,small,ndiskresid,
     9          max_amplitudes_size,ishellsaddr,iamplitaddr,iresultaddr,
     *       igranules,igranulesize,bl(igran2pair),bl(ipair2gran),split,
     1          bl(iscr),bl(mapf2s),ik_capacity,tot_xmax,list_copy,
     2          nstrong_copy,isympairs,ij_unique)
      else if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'tt')
     *     then

       call TransOneShell_tj(ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),bl(lmp2int),bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), bl(1),mulam, mulamd, 
     *                    nrow,
     3                    ncol, bl(irow), bl(icol), bl(irow1),bl(icol1),
     4                    ENonZero,smal,lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .true.,' ',mygid,int_kind,byt8,nmo,vorb,igranules,
     8               igranulesize,bl(igran2pair),bl(ipair2gran),irec_no,
     9               irec_last,nbf)
       else
       call TransOneShell_tj(ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),bl(lmp2int),bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalfmulam), bl(ihalflammu),
     *                    mulam, mulamd,  nrow,
     3                    ncol, bl(irow), bl(icol), bl(irow1),bl(icol1),
     4                    ENonZero,smal,lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .true.,' ',mygid,int_kind,byt8,nmo,vorb,igranules,
     8               igranulesize,bl(igran2pair),bl(ipair2gran),irec_no,
     9               irec_last,nbf)
      endif
c
      call secund(tt4)
      call elapsec(telap4)
      ttrans=ttrans+tt4-tt3
      elaptrans=elaptrans+telap4-telap3
c......................................................................
      call dynamic_retmem(1)          ! lmp2int 
c......................................................................
       enddo    ! over jtime (selected jcf belonging to jcs shell )
       enddo    ! over itime (selected icf belonging to ics shell )
c......................................................................
#endif
      end do       !  over ics,jcs shell
      call dynamic_retmark 
      if (Int_kind.eq.'e') then
      cache=oldcache
      endif
      call sync_barrier
      call para_initsend
      call para_pack_int(mygid,1)
      call para_pack_int(irecadr,irecidx)
      call para_send_pack(0,TxCCSDRes)
c......................................................................
c timing here is irrelevant :
      call secund(txxx2)
      txxx0=txxx1
      call term_info(thresh,txxx2-txxx1,txxx1-txxx0,'fock')
c==========================
      irecskip=irecidx
c     if (int_kind.eq.'tx'.or.int_kind.eq.'tc') irecskip=nval*nval
      islaves=nslv
      af=.true.
      call fafreopen(ndisk)
      call fafCloses
      call fafCreates(ndisktr)
      call sync_barrier
c Initialize PVM:
      call elapsec(xlatency0)
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(irecadr,irecskip*islaves)
      call dynamic_lock(irecadr,iii)
      call para_unpack_int(ndisktr,1)
      call elapsec(xlatency1)
      tot_latency2=xlatency1-xlatency0
c Do job:
      if (int_kind.ne.'e'.and.int_kind.ne.'tt') igranulesize=1
      call dynamic_getmem(lbin*igranulesize,ibin8) 
      call dynamic_getmem(lbin*igranulesize,ibin4adr) 
      call dynamic_getmem(lbin*igranulesize,ibin1adr) 
      call dynamic_getmem(lbin*igranulesize,indexadr) 
      do
        iready=1
        call para_initsend
        call para_pack_int(iready,1)
        call para_pack_int(mytid,1)
        call para_send_pack(0,TxCCSDReq)
c
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(iwork,1)
        if (iwork.eq.1) then
        call para_unpack_int(ij,1)
        if (Int_kind.eq.'e'.or.Int_kind.eq.'tt') then
       call transformatgr(ij,islaves,irecadr,irecskip,lbin,bl(indexadr),
     *        bl(ibin1adr),bl(ibin4adr),bl(ibin8),ndisk,ndisktr,npairs,
     *        ncf,bl(ixadr),thresh,af,byt8,nmo,vorb,int_kind,
     *   igranules,igranulesize,bl(igran2pair),bl(ipair2gran),list_copy,
     *   nbf)
        else
       call transformator(ij,islaves,irecadr,irecskip,lbin,bl(indexadr),
     *        bl(ibin1adr),bl(ibin4adr),bl(ibin8),ndisk,ndisktr,npairs,
     *        ncf,bl(ixadr),thresh,af,byt8,nmo,vorb,int_kind,nbf)
        endif
        else
          exit
        endif
      enddo
      call sync_barrier
      call fafreopen(ndisktr)
      info=0
      inf=0
      call dynamic_retmem(4) 
      call secund(total_cpu1)
      total_cpu=total_cpu1-total_cpu0
      call para_reduce(total_cpu,1,TxCCSDRedu1)
      call para_barrier
      call para_reduce(tot_latency1,1,TxCCSDRedu2)
      call para_barrier
      call para_reduce(tot_latency2,1,TxCCSDRedu3)
      call para_barrier
      call para_reduce(zeroed,2,TxCCSDRedu1)
      call para_barrier
      xiare_zeros=iare_zeros
      call para_reduce(xiare_zeros,1,TxCCSDRedu2)
      call para_barrier
      call para_reduce(strace,14,TxCCSDRedu3)
      call para_barrier
      call para_reduce(xint_time,ixint_time,TxCCSDRedu1)
      call para_barrier
      call para_recv_bcast(iwork,TxCCSDQuit)
        if (info.lt.0) print *,'Error: ', __LINE__
c==========================
c......................................................................
      nrec=nrec+irec
      if(iprnt.gt.2) then
         write(iout,*) ' Total number of records written=',nrec
      endif
c
      if(iprnt.gt.2) then
c        Calculate statistics
         PerCent=100.0d0*ENonZero/(dble(mulam+mulamd)*ncf**2)
         write(iout,40) PerCent
      endif
c
  40  format(
     *'  Percentage of the AO matrix used in the transformation=',f7.3/)
c......................................................................
  50  continue
c  release memory
      call retmark
      call dynamic_retmark ! mmark just before int_coul
c Slave!
#if 1
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e'
     *                   .or.int_kind.eq.'tt') then
      call dynamic_matrem('halftra')
      else
      call dynamic_matrem('halflammu')
      call dynamic_matrem('halfmulam')
      endif
#endif
      call retmark
      call dynamic_retmark ! mmark after irec allocation
c......................................................................
      write(iout,61)
   61 format('  CPU & Elapsed timings in the coulomb integrals module ')
      write(iout,*) '  '
c
      write(iout,62) tinteg/sixty, elapint/sixty
   62 format('  Coulomb Integrals =',f8.2,' and ',f8.2,' minutes')
      write(iout,64) ttrans/sixty, elaptrans/sixty
   64 format('  First half transf. =',f8.2,' and ',f8.2,' minutes')
      call flush(iout)
ckw---------------------------------------------------------
      if(.not.restrt .and. iprnt.ge.3) then
         call getrval('timeblks',timeblks)
         write(91,123) timeblks/sixty
  123    format('CPU time for new blocking =',f8.2)
      endif
ckw---------------------------------------------------------
      call secund(tt1)
      call elapsec(elaps1)
      if(iprnt.gt.2) then
        write(iout,*) ' CPU and elapsed time for integrals+first half ',
     1  'transformation ',(tt1-tt0)/sixty,(elaps1-elaps0)/sixty
      end if
c..................................................
c---------------------------------------------------------------------
      int_kinder='xxxxxxx'
      call symmon
c---------------------------------------------------------------------
      call dynamic_retmark 
      call matremark
      call retmark
      if (int_kind.eq.'tt') call fafCloses
      end
!=======================================================================
      subroutine initsingles_master(nval,ncf,singles,nmo,vorba,
     *                              restart,ndisk_dump,nbf)
c When called first time his subroutine allocates memory for internal storage.
c It also distributes the data to slaves via broadcast. When called second,
c third time and so on, it also distributes its data to slaves and copied
c singles from "writable" storage into "readable" storage.
c Parameters:
c   nval    - number of correlated orbitals
c   ncf     - number of contracted basis functions
c   singles - .true. singles contributions are calculated. depending on the
c             values of other variables like cc and CCSD different set of 
c             singles contributions will be calculated.
c   nmo     - number of occupied HF orbitals
cI  vorba   - .true.: the dimension of virtual space is equal to no of virt.
c            orbitals, and MO form of integral matrices is used, otherwise,
c            the virtual space has ncf dimension and all calculations are 
c            performed for matrices in AO form.
c This subroutine shares common block with readsingles & writesingles.
c commmon block /singleskeep/
c   iaddress0 - the memory location, where amplitudes to be read are kept
c   ncfa      - = ncf
c   singlesa  - = singles
c   iaddress1 - memory location for writing singles amplitudes
c   idimen    - the dimension of integral and amplitudes matrices
c   vorb      - = vorba
      use memory
      use newpara
      use kinds
      use ccounters
      implicit none
      integer nval,ncf,nmo,nbf
      logical singles,vorba,restart
      integer ndisk_dump
c
      logical singlesa
      integer iaddress0,iaddress1,ncfa,i,idimen,info
      integer icheck0,icheck1
      logical vorb
      common /singleskeep/ iaddress0,ncfa,singlesa,iaddress1,idimen,
     *                     vorb
      save /singleskeep/
c
      ic_initsingles_master=ic_initsingles_master+1
c
      vorb=vorba
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ncfa=ncf
      singlesa=singles
      if (singles.and.ic_initsingles_master.eq.1) then
        call dynamic_getmem(nval*idimen,iaddress0)
        call dynamic_getmem(nval*idimen,iaddress1)
        do i=0,nval*idimen-1
        bl(iaddress0+i)=0.0d0
c       bl(iaddress0+i)=1.d-10*dble(irand())
        enddo
      endif
      if (singles.and.ic_initsingles_master.gt.1) then
        do i=0,nval*idimen-1
        bl(iaddress0+i)=bl(iaddress1+i)
        bl(iaddress1+i)=0.0d0
        enddo
      endif
      if (singles) call symm2_s_ampl(idimen, bl(iaddress0))
      if (singles.and.restart.and.ic_initsingles_master.eq.1) then
        call general_read1(ndisk_dump,bl(iaddress0),nval*idimen*8)
      endif
      if (singles) then
        call para_initsend
        call para_pack_real(bl(iaddress0),nval*idimen)
        call para_bcast_pack(TxCCSDSing)
      endif
      call flush(6)
      end
C===============================================================================
      subroutine initsingles_slave(nval,ncf,singles,nmo,vorba,iMasID,
     *                             nbf)
c This subroutine is almost identical to initsingles_master, but it receives
c the amplitudes
      use memory
      use newpara
      use kinds
      use ccounters
      implicit none
      integer nval,ncf,nmo,nbf
      logical singles,vorba
c
      logical singlesa
      integer iaddress0,iaddress1,ncfa,i,idimen,info,iMasID
      integer icheck0,icheck1
      logical vorb
      common /singleskeep/ iaddress0,ncfa,singlesa,iaddress1,idimen,
     *                     vorb
      save /singleskeep/
c
      ic_initsingles_slave=ic_initsingles_slave+1
c
      vorb=vorba
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ncfa=ncf
      singlesa=singles
      if (singles.and.ic_initsingles_slave.eq.1) then
        call dynamic_getmem(nval*idimen,iaddress0)
        call dynamic_getmem(nval*idimen,iaddress1)
      else if (singles.and.ic_initsingles_slave.gt.1) then
        call dynamic_show_alloc(bl(iaddress0),icheck0)
        if (icheck0.ne.nval*idimen) 
     *    call nerror(1,'initsingles','Memory allocation error',
     *                ic_initsingles_slave,0)
        call dynamic_show_alloc(bl(iaddress1),icheck1)
        if (icheck1.ne.nval*idimen) 
     *    call nerror(1,'initsingles','Memory allocation error',
     *                ic_initsingles_slave,0)
      endif
c     if (ic_initsingles_slave.gt.1) then
c       do i=0,nval*idimen-1
c       bl(iaddress0+i)=bl(iaddress1+i)
c       bl(iaddress1+i)=0.0d0
c       enddo
c     endif
      if (singles) then
        call para_recv_bcastpack(TxCCSDSing)
        call para_unpack_real(bl(iaddress0),nval*idimen)
      endif
      end
c==============================================================================
      subroutine CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd_inp,
     *                          vorbpar,nmo,iMasID,     mygid,ivalpair,
     *                          ivirpair,ifp,nsym,      nval,  nvirt,
     *                          ij_unique,isympairs,    nbf)
c This subroutine initializes file with coefficients (amplitudes) on the
c slave side. This initialization encompasses receiving the file units
c numbers and preparation of local cache files if "cache" variable (passed
c via common block) is true. It also sets up common block 
c /ampldisk/ which is shared with coefwrite and coefread
c See CoefInit for single mode for common and parameters description.
c This subroutine assds its own:
c Parameters:
c   iMasID - the master TID
c   mygid  -  my (slave) TID
      use memory
      use newpara
      use kinds
      use ccounters
      implicit real*8 (a-h,o-z)
      character*256 filenamea,filenameb,scrfile,iternumbera,iternumberb
      character*256 filenamer
      character*2 ch_mygid
      logical gauss_seidel,nodisk,nodiskk,ccsd,ccsd_inp,fileopen
      logical vorb,vorbpar,af
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
      ic_CoefInit_slave=ic_CoefInit_slave+1
      icallno_propagate=ic_CoefInit_slave
      npairs_propagate=npairs
      ccsd=ccsd_inp
c     print *,'Coupled: ',ccsd
      if (ic_CoefInit_slave.eq.1) call dynamic_getmem(idimen*idimen,
     *                                                ireturn)
      if (nodisk) then
        STOP 'No support for nodisk tmp, should be added?'
        if (ic_CoefInit_slave.eq.1) then 
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
c
      internal_addr=0
      if (gauss_seidel) then
      if (ic_CoefInit_slave.eq.1) then
      ncfa=ncf
      nodiskk=nodisk
      endif
      call fafCreates(itrala)
      else
c=================
      if (ic_CoefInit_slave.eq.1) then
        ncfa=ncf
        ndiskr=62
        nodiskk=nodisk
      endif
      if (ic_CoefInit_slave.gt.1) then
      endif
      call fafCreates(itrala)
      call fafCloses
      call fafCreates(itrala)
c=================
      endif
      call sync_barrier
      call sync_barrier
      call para_recv_bcastpack(TxCoefInit)
      call para_unpack_int(ndiska,1)
      call para_unpack_int(ndiskb,1)
      call fafreopen(ndiska)
      end
C===============================================================================
      subroutine CCalphaslave(ncf,   nval,  irecadrx,npairs,    ndiskx,
     *                        lbinx, thresh,doit,    byt8,      ccsd,
     *                        nmo,   vorb,  elaps,   ndiskalpha,af,
     *                        iMasID,mytid, nbf)
c This is slave's part of CCalphaonce (see description there). It is 
c essentially the internal loop over K and external loop over T's is steered
c by master which sends the starting and closing indices.
c Additional parameters:
c   iMasID - the master TID
c   mygid  -  my (slave) TID
      use memory
      use newpara
      use kinds
      implicit none
      integer ncf,nval,irecadrx(*),npairs,ndiskx,lbinx,nmo,ndiskalpha
      integer iMasID,mytid,nbf
      REAL*8 thresh,alpha,elaps(*)
      integer iexchAO,icoeffAO,mataddr,i,j,k,l,iaddress
      integer inval,icurraddr,nval2,nval3,initial_chunk
      logical doit,idoit,byt8,ccsd,vorb,af
      parameter (initial_chunk=25)
      REAL*8 elem,elaps0,elaps1
      common /CCAlphaCom/iaddress,inval,idoit
c new
      integer memavail,lcore,iad,ncf2,number_matrices,ichunk_T,ichunk_K
      integer iKstore,iTstore,iresult,iKpass,iTpass,iT,iK,iTstart,iTstop
      integer iTsize,iKstart,iKstop,iKsize,ij,kl,number_matrices1,idimen
      integer idimen2,ialphabuffer,nslvdiv,info,iwork
      integer itransposition,ifree,imastgid
c
      save /CCAlphaCom/
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
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ichunk_K,1)
      call para_unpack_int(ichunk_T,1)
      call para_unpack_int(ndiskalpha,1)
      call dynamic_show_free(ifree)
c
c     call dynamic_getmem(idimen2*ichunk_K,iKstore)
c     call dynamic_getmem(idimen2*ichunk_T,iTstore)
c     call dynamic_getmem(ichunk_T*ichunk_K,iresult)
c     call dynamic_getmem(ichunk_T*nval*nval,ialphabuffer)
      call dynamic_getmem(idimen2*ichunk_K,iKstore)
      call dynamic_getmem(idimen2*ichunk_T,iTstore)
      call dynamic_getmem(ichunk_T*ichunk_K,iresult)
      call dynamic_getmem(ichunk_T*nval*nval,ialphabuffer)
      iKpass=npairs/ichunk_K
      if (mod(npairs,ichunk_K).ne.0) iKpass=iKpass+1
c
c Loops over passes:
      do
        call para_send(mytid,0,TxCCSDAReq)
        call para_recv_pack(imastgid,TxCCSDAJob)
        call para_unpack_int(iwork,1)
        if (iwork.eq.0) exit
        call para_unpack_int(iTstart,1)
        call para_unpack_int(iTstop,1)
        call para_unpack_int(iTsize,1)
c
        call elapsec(elaps0)
        call store(iTstart,iTstop,idimen,nval,iTstore)
        call elapsec(elaps1)
        elaps(1)=elaps(1)+elaps1-elaps0 ! amplitudes storing
        call matconn('BigT','r',idimen2,iTsize,iTstore)
        do iK=1,iKpass
          iKstart=ichunk_K*(iK-1)+1
          iKstop =ichunk_K*iK
          if (iKstop.gt.npairs) iKstop=npairs
          iKsize=iKstop-iKstart+1
          call elapsec(elaps0)
          call storeK(irecadrx,ndiskx,lbinx,thresh,iKstart,iKstop,
     *                ncf,nval,iKstore,byt8,nmo,vorb,nbf)
          call elapsec(elaps1)
          elaps(2)=elaps(2)+elaps1-elaps0 ! Kij storing
          call matconn('BigK','r',idimen2,iKsize,iKstore)
          do itransposition=1,2
          call elapsec(elaps0)
          if (itransposition.eq.2) then
            call BigTransposition(iKstart,iKstop,idimen,nval,iKstore)
          endif
          call elapsec(elaps1)
          elaps(3)=elaps(3)+elaps1-elaps0 ! transposition
          call elapsec(elaps0)
          call matconn('result','r',iTsize,iKsize,iresult)
          call matmmul2('BigT','BigK','result','t','n','n')
          call elapsec(elaps1)
          elaps(4)=elaps(4)+elaps1-elaps0 ! multiplication
          call elapsec(elaps0)
          call putresult(idimen,nval,bl(iresult),iTstart,iTstop,iKstart,
     *                   iKstop,ccsd,ialphabuffer,ndiskalpha,af,
     *                   itransposition)
          call elapsec(elaps1)
          elaps(5)=elaps(5)+elaps1-elaps0 ! putresult
          call matdisc('result')
          enddo
          call matdisc('BigK')
        enddo
        call matdisc('BigT')
      enddo
c
      call dynamic_retmem(4)
c
      call matremark
      call retmark
      call dynamic_retmark
      end
C===============================================================================
      subroutine B41a_G_Slave(ncf,   nval,      norecalc,cc,    ccsd,
     *                        nfirst,F,         nmo,     vorb,  ndisk,
     *                        af,    ndiskalpha,betam,   iMasID,mytid,
     *                        nbf)
c This is slave's part of B41a_G_Slave (see description there). It is 
c essentially the internal loop over T and external loop over R's is steered
c by master which sends the starting and closing indices.
c Additional parameters:
c   iMasID - the master TID
c   mygid  -  my (slave) TID
      use memory
      use newpara
      use kinds
      implicit none
      integer ncf,nval,nfirst,nmo,ndisk,ndiskalpha,iMasID,mytid,nbf
      logical norecalc,cc,ccsd,vorb,af
C
      integer iad,lcore,memavail,initial_chunk,number_matrices,ncf2
      integer npairs,ichunk_T,number_matrices1,ichunk_R,iTstore,iRstore
      integer ialphas,iT,iTpass,iTstart,iTstop,iTsize,iR,iRpass
      integer iRstart,iRstop,iRsize,itwo,idimen,idimen2,ialphabase
      integer info,iwork,imastgid
      logical transposed,reset
      real*8 F(ncf,ncf),betam(nval,nval)
      parameter (initial_chunk=25)
C
      reset=.true. ! for storeAlpha
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
c     call dynamic_mmark
      call matmark
      ncf2=ncf*ncf
      idimen2=idimen*idimen
      npairs=nval*(nval+1)/2
c
      call fafCloses
      call fafCreates(ndisk)
      call sync_barrier
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ichunk_R,1)
      call para_unpack_int(ichunk_T,1)
      call para_unpack_int(ndisk,1,1,info) ! file s
c
c     call dynamic_getmem(idimen2*ichunk_T,iTstore)
c     call dynamic_getmem(idimen2*ichunk_R,iRstore)
c     call dynamic_getmem(ichunk_T*ichunk_R,ialphas)
c     call dynamic_getmem(nval*nval*ichunk_R,ialphabase)
      call dynamic_getmem(idimen2*ichunk_T,iTstore)
      call dynamic_getmem(idimen2*ichunk_R,iRstore)
      call dynamic_getmem(ichunk_T*ichunk_R,ialphas)
      call dynamic_getmem(nval*nval*ichunk_R,ialphabase)
c calculate passes:
      iTpass=npairs/ichunk_T
      if (mod(npairs,ichunk_T).ne.0) iTpass=iTpass+1
c Loops over passes:
      iwork=1
      do
        call para_send(mytid,0,TxCCSDAReq)
        call para_recv_pack(imastgid,TxCCSDAJob)
        call para_unpack_int(iwork,1)
        if (iwork.eq.0) exit
        call para_unpack_int(iRstart,1)
        call para_unpack_int(iRstop,1)
        call para_unpack_int(iRsize,1)
c
        call matconn('BigPrR','r',idimen2,iRsize,iRstore)
        call matzero('BigPrR')
c
        do iT=1,iTpass
          iTstart=ichunk_T*(iT-1)+1
          iTstop =ichunk_T*iT
          if (iTstop.gt.npairs) iTstop=npairs
          iTsize=iTstop-iTstart+1
          call store(iTstart,iTstop,idimen,nval,iTstore)
          call matconn('BigT','r',idimen2,iTsize,iTstore)
c
          transposition: do itwo=1,2
          if (itwo.eq.1) then 
c Without transposition:
            transposed=.false.
          else
c With transposition:
          transposed=.true.
          call BigTransposition(iTstart,iTstop,idimen,nval,iTstore)
          endif
c
          call storeAlpha(ncf,nval,iTstart,iTstop,iRstart,
     *          iRstop,bl(ialphas),transposed,cc,ccsd,nfirst,F,nmo,vorb,
     *          ialphabase,ndiskalpha,betam,reset,af,nbf)
          call matconn('alpha','r',iTsize,iRsize,ialphas)
          call matmmul2('BigT','alpha','BigPrR','n','n','A')
c===========================
          call matdisc('alpha')
        enddo transposition
        call matdisc('BigT')
        enddo
c
        call writeresult(ndisk,iRstart,iRstop,idimen,nval,iRstore,af)
        call matdisc('BigPrR')
c
      enddo
c Synchronize:
c
      call dynamic_retmem(4)
c
      call matremark
      call retmark
c     call dynamic_retmark
      call para_barrier
      end
c==============================================================================
      subroutine 
     *    EEO_vector_slave(nval,  irecadre,ndiske,  lbine, irecadrx,
     *                     ndiskx,lbinx,   irecadrc,ndiskc,lbinc,
     *                     ncf,   npairs,  thresh,  nfirst,nlast,
     *                     byt8,  nmo,     vorb,    nslv1,  iMasID,
     *                     mytid, list,    nbf)
c This is slave's part of CCalphaonce (see description there). It is 
c essentially the internal loop over 'i' and external loop over k's is steered
c by master which sends the indices.
c Additional parameters:
c   iMasID - the master TID
c   mygid  -  my (slave) TID
      use memory
      use newpara
      use kinds
      implicit none
      integer nval,npairs,nbf
      integer ncf,nfirst,nlast,nmo,nslv1,iMasID,mytid
      integer irecadre(npairs),ndiske,lbine
      integer irecadrx(npairs),ndiskx,lbinx
      integer irecadrc(npairs),ndiskc,lbinc
      real*8 thresh
      logical byt8,vorb
      integer list(*)
c
      integer i,k,iEEOread,iEEOconv,mataddr,idimen,intermid,ires
      integer iwork,info,ibufid,imastgid
      integer iworker1,iworker2,icorrect,ising_vec,ik
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
      call matsub('occu','cano',nfirst,nlast)
c
      call dynamic_matdef('EEOread','q',ncf,ncf)
      iEEOread=mataddr('EEOread')
c
      call dynamic_matdef('EEOconv','q',ncf,ncf)
      iEEOconv=mataddr('EEOconv')
c
      call dynamic_matdef('intermid','r',ncf,nval)
      call matzero('intermid')
      intermid=mataddr('intermid')
c
      call dynamic_matdef('worker1','q',idimen,idimen)
      iworker1=mataddr('worker1')
      call dynamic_matdef('worker2','q',idimen,idimen)
      iworker2=mataddr('worker2')
      call dynamic_matdef('correct','r',idimen,nval)
      call matzero('correct')
      icorrect=mataddr('correct')
      call dynamic_matdef('sing_vec','r',idimen,1)
      ising_vec=mataddr('sing_vec')
c
      iwork=1
      do
        call para_send(mytid,0,TxCCSDReq)
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(iwork,1)
        if (iwork.eq.0) exit
        call para_unpack_int(k,1)
        call para_unpack_int(i,1)
          call ExtrOne(i,    k,     irecadre,   npairs,    ndiske,
     *                 ncf,  lbine, thresh,     byt8,      'e',
     *                 'ao', nmo,   vorb,       bl(iEEOread))
        call ijtriangle(k,i,ik)
c *****************************
        if (list(ik).ne.1) then
          call ExtrOne(i,    k,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'mo', nmo,   vorb,       bl(iworker1))
          call ExtrOne(i,    k,     irecadrc,   npairs,    ndiskc,
     *                 ncf,  lbinc, thresh,     byt8,      'c',
     *                 'mo', nmo,   vorb,       bl(iworker2))
          call matscal('worker1',2.0d0)
          call matadd1('worker2',-1.0d0,'worker1')
        endif
c *****************************
          call CCLoperator(bl(iEEOread),ncf,bl(iEEOconv))
c
          call matsub('onecol','occu',k,k)
          ires=intermid+(i-1)*ncf
          call matconn('result','r',ncf,1,ires)
          call matmmul2('EEOconv','onecol','result','n','n','a')
          call matdisc('result')
          call matdisc('onecol')
c
c *****************************
        if (list(ik).ne.1) then
          ires=icorrect+(i-1)*idimen
          call matconn('result','r',idimen,1,ires)
          call readsingles(k,bl(ising_vec))
          call matmmul2('worker1','sing_vec','result','n','n','a')
          call matdisc('result')
        endif
c *****************************
c
        if (i.ne.k) then
          call matsub('onecol','occu',i,i)
          ires=intermid+(k-1)*ncf
          call matconn('result','r',ncf,1,ires)
          call matmmul2('EEOconv','onecol','result','t','n','a')
          call matdisc('result')
          call matdisc('onecol')
c
c *****************************
        if (list(ik).ne.1) then
          ires=icorrect+(k-1)*idimen
          call matconn('result','r',idimen,1,ires)
          call readsingles(i,bl(ising_vec))
          call matmmul2('worker1','sing_vec','result','t','n','a')
          call matdisc('result')
        endif
c *****************************
        endif
      enddo
      call sync_barrier
      call para_barrier
      do
        call para_send(mytid,0,TxCCSDReq)
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(iwork,1)
        if (iwork.eq.0) exit
        call para_unpack_int(k,1)
        if (iwork.eq.0) exit
          call ExtrOne(k,    k,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'mo', nmo,   vorb,       bl(iworker1))
          call ExtrOne(k,    k,     irecadrc,   npairs,    ndiskc,
     *                 ncf,  lbinc, thresh,     byt8,      'c',
     *                 'mo', nmo,   vorb,       bl(iworker2))
        call matscal('worker2',2.0d0)
        call matadd1('worker1',-1.0d0,'worker2')
        do i=1,nval
          call ijtriangle(k,i,ik)
c *****************************
          if (list(ik).ne.1) then
          ires=icorrect+(i-1)*idimen
          call matconn('result','r',idimen,1,ires)
          call readsingles(i,bl(ising_vec))
          call matmmul2('worker2','sing_vec','result','n','n','a')
          call matdisc('result')
          endif
c *****************************
        enddo
      enddo
      if (vorb) then
        call matmmul2('genvirt','intermid','correct','t','n','a')
      else
        call matadd('intermid','correct')
      endif
      call sync_barrier
      call para_reduce(bl(icorrect),idimen*nval,TxCCSDRes)
c
      call matremark
      call retmark
      call dynamic_retmark
      call para_barrier
      end
c==============================================================================
      subroutine Generate_MP2_G_SLA(nval, ncf, ndisk, nmo,  nfirst,
     *                              vorb, af,  F,     nslv1, iMasID,
     *                              mytid,nbf)
c Paralelized equivalent of single processor Generate_MP2_G. Simply 
c parallelized external loop
c Additional parameters:
c   nslv   - total number of slaves
c   iMasID - the master TID
c   mygid  - my (slave) TID
      use memory
      use newpara
      use kinds
      implicit none
      integer nval,ncf,ndisk,nmo,nfirst,nslv1,iMasID,mytid,nbf
      logical vorb,af
      real*8 F(ncf,ncf)
c local:
c
      integer i,j,k,iresult,iamplit,ka,ir,ij,idimen,ifp,jfp,itmp,ifock
      integer istatus,info,ibufid,iwork,imastgid,memo,itmpfock,isubad
      integer istart1,istop1,imaxsize1,isize1,npass1,ipass1
      integer istart2,istop2,imaxsize2,isize2,npass2,ipass2
c
      call dynamic_mmark()
      call matmark()
      call dynamic_getmem(nval*nval,ifock)
      call matconn('fock_occ','r',nval,nval,ifock)
      do j=1,nval
        jfp=j+nfirst-1
        do i=1,nval
          ifp=i+nfirst-1
          bl(ifock+(i-1)+(j-1)*nval)=F(ifp,jfp)
        enddo
      enddo
      call matscal('fock_occ',-1.0d0)
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      call dynamic_getmem(nval*nval,itmpfock)
c
c Check multipassing1:
      call dynamic_show_free(memo)
      imaxsize1=(memo-5*idimen*idimen)/(idimen*idimen)
      if (imaxsize1.lt.1) STOP 'Insufficient amount of memory'
      if (imaxsize1.le.3) 
     *        print *,'Warning, low memory may cause efficiency problem'
      if (imaxsize1.gt.nval) imaxsize1=nval
c     imaxsize1=7
      npass1=nval/imaxsize1
      if (mod(nval,imaxsize1).ne.0) npass1=npass1+1
      call optimal_isize(imaxsize1,nval)
      call dynamic_getmem(idimen*idimen*imaxsize1,iamplit)
c
c Check multipassing2:
      call dynamic_show_free(memo)
      imaxsize2=memo/(idimen*idimen)
      if (imaxsize2.gt.nval) imaxsize2=nval
c     imaxsize2=8
      if (imaxsize2.lt.1) STOP 'Insufficient amount of memory'
      if (imaxsize2.le.3) 
     *        print *,'Warning, low memory may cause efficiency problem'
      npass2=nval/imaxsize2
      if (mod(nval,imaxsize2).ne.0) npass2=npass2+1
      call optimal_isize(imaxsize2,nval)
c
      call dynamic_getmem(idimen*idimen*imaxsize2,iresult)
      call para_barrier
      do            ! result column
        call para_send(mytid,0,MP2_G_GEN_req)
        call para_recv_pack(imastgid,MP2_G_GEN_job)
        call para_unpack_int(iwork,1)
        if (iwork.eq.0) exit
        call para_unpack_int(j,1)
        

        do ipass1=1,npass1
          istart1=(ipass1-1)*imaxsize1+1
          istop1= (ipass1  )*imaxsize1
          if (istop1 .gt.nval) istop1=nval
          if (istart1.gt.nval) STOP 'Error sdsdfu853kgbzxftgy'
          isize1=istop1-istart1+1
          do k=istart1,istop1
            ka=iamplit+(k-istart1)*idimen*idimen
            call CoefRead('tt',k,j,itmp)
            call tfer(bl(itmp),bl(ka),idimen*idimen)
          enddo
          do ipass2=1,npass2
            istart2=(ipass2-1)*imaxsize2+1
            istop2= (ipass2  )*imaxsize2
            if (istop2 .gt.nval) istop2=nval
            if (istart2.gt.nval) STOP 'Error sdsdfu853kgbzxftgy'
            isize2=istop2-istart2+1
c
            if (ipass1.eq.1) then
              call zeroit(bl(iresult),idimen*idimen*imaxsize2)
            else
            do i=istart2,istop2 ! result read
              ir=iresult+(i-istart2)*idimen*idimen
              ij=(j-1)*nval+i
              if (af) then
                call fafread(ndisk,bl(ir),8,idimen*idimen,1,ij,istatus)
              else
                call reader(ndisk,idimen,ij,bl(ir))
              endif
            enddo
            endif
c
            call matconn('result','r',idimen*idimen,isize2,iresult)
            call matconn('amplit','r',idimen*idimen,isize1,iamplit)
            call take_strip(bl(ifock),bl(itmpfock),nval,  nval, istart1,
     *                      istop1,   istart2,     istop2)
            call matconn('subfock','r',isize1,isize2,itmpfock)
            call matmmul2('amplit','subfock','result','n','n','a')
            call matdisc('subfock')
            call matdisc('amplit')
            call matdisc('result')
c
            do i=istart2,istop2 ! result row
              ir=iresult+(i-istart2)*idimen*idimen
              ij=(j-1)*nval+i
              if (af) then
                call fafwrite(ndisk,bl(ir),8,idimen*idimen,1,ij,istatus)
              else
              call writer(ndisk,idimen,ij,bl(ir))
              endif
            enddo
c
          enddo ! pass2
        enddo ! pass1
      enddo
      call matremark()
      call dynamic_retmark()
      end
c==============================================================================
      subroutine Kijk_Vec_Slave(irecadrx,npairs,ndiskx,lbinx, thresh,
     *                          byt8,    ndiska,ncfa,  nmo,   nvala,
     *                          vorba,   afa,   mygid, ifirst,last,
     *                          nslv1,    iMasID,mytid,nbfa)
c Paralelized equivalent of single processor Kijk_Vec. Simply 
c parallelized external loop. Refer single processor subroutine
c for description.
c Additional parameters:
c   nslv   - total number of slaves
c   iMasID - the master TID
c   mygid  - my (slave) TID
      use memory
      use newpara
      use kinds
      implicit none
      integer irecadrx(*),npairs,ndiskx,lbinx,nbfa
      real*8 thresh
      integer ndiska,ncfa,nmo,nvala,mygid,nslv1,iMasID,mytid,ifirst,last
      logical byt8,vorba,afa
c
      integer ndisk,ncf,nval,idimen,ndisk_cache,itempk,len,len1,irecord
      integer nrec,i,j,k,info,mataddr,imatrix,iresult,iresultMO,ibufid
      integer iwork,nval2,iaddress,nrec1,imastgid,nbf
      character*256 scrfile,filename
      character*3 ch_mygid
      logical vorb,af
      common /kijk/ ndisk,ncf,nval,idimen,vorb,af,ndisk_cache,nbf
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      save /kijk/
c
      call mmark
      call dynamic_mmark
      call matmark
cA
      nval2=nvala*nvala
      if (vorba) then
        idimen=nbfa-nmo
      else
        idimen=ncfa
      endif
c FRAGILE!
      ndisk_cache=76
c
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                   ! 80 is starting lenght
        write(ch_mygid,'(I3.3)') mygid
        filename=scrfile(1:len)//'.kijk--'//ch_mygid
        len1=len+10
        irecord=idimen*8
        OPEN(UNIT=ndisk_cache,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
c
      ncf=ncfa
      nbf=nbfa
      nval=nvala
      vorb=vorba
      af=afa
c
      call dynamic_matdef('matrix','q',ncf,ncf)
      imatrix=mataddr('matrix')
      call dynamic_matdef('result','r',ncf,nval)
      iresult=mataddr('result')
      call matsub('occu','cano',ifirst,last)
      if (vorb) then
        call dynamic_matdef('resultMO','r',idimen,nval)
        iresultMO=mataddr('resultMO')
      endif
c
      call fafCreates(ndiska)
      call sync_barrier
      call para_recv_bcast(ndiska,TxCCSDInit)
      ndisk=ndiska
      do
c PVM header:
        call para_send(mytid,0,TxCCSDReq)
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(iwork,1)
        if (iwork.ne.1) exit
        call para_unpack_int(i,1)
        call para_unpack_int(j,1)
c The work:
          call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'ao', nmo,   vorb,       bl(imatrix))
          call matmmult('matrix','occu','result')
          if (vorb) then
            call matmmul2('genvirt','result','resultMO','t','n','n')
          endif
          nrec=(i-1)*nval+j
          if (vorb) then
            iaddress=iresultMO
          else
            iaddress=iresult
          endif
          call fafWrite(ndiska,bl(iaddress),8,idimen*nval,1,nrec,info)
      enddo
      call para_barrier
      call dynamic_matdef('tempk','r',idimen,nval)
      itempk=mataddr('tempk')
      nrec=0
      nrec1=0
      do i=1,nval
        do j=1,nval
          nrec=nrec+1
          call fafRead(ndiska,bl(itempk),8,idimen*nval,1,nrec,info)
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
      call retmark
      call dynamic_retmark
      call matremark
      call sync_barrier
      end
c==============================================================================
      subroutine KijklInit_Slave(irecadrx,npairs,ndiskx,lbinx,ncf,
     *                           thresh,  ifirst,last,  F,    byt8,
     *                           nmo,     vorb,  af,    nslv1, iMasID,
     *                           mytid)
c Paralelized equivalent of single processor KijklInit. Simply 
c parallelized external loop. Refer single processor subroutine
c for description.
c Additional parameters:
c   nslv   - total number of slaves
c   iMasID - the master TID
c   mygid  - my (slave) TID
      use memory
      use newpara
      use kinds
      implicit none
      integer npairs,ndiskx,lbinx,ncf,ifirst,last,nmo
      integer irecadrx(npairs),nslv1,iMasID,mytid
      real*8 thresh,F(ncf,ncf)  ! F-ock in MO basis
      logical byt8,vorb,af
c
      integer nval,imatrix,mataddr,ikij,ij,i,j,istatus,ndisk,info
      integer iwork,ibufid,imastgid
c
      nval=last-ifirst+1
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('matrix','q',ncf,ncf)
      imatrix=mataddr('matrix')
      call dynamic_matdef('kij','q',nval,nval)
      ikij=mataddr('kij')
      call matsub('occu','cano',ifirst,last)
      call para_barrier
      call para_recv_bcast(ndisk,TxCCSDInit)
      do
C PVM header
        call para_send(mytid,0,TxCCSDReq)
        call para_recv_pack(imastgid,TxCCSDJob)
        call para_unpack_int(iwork,1)
        if (iwork.ne.1) exit
        call para_unpack_int(i,1)
        call para_unpack_int(j,1)
        call para_unpack_int(ij,1)
C job:
        call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *               ncf,  lbinx, thresh,     byt8,      'x',
     *               'ao', nmo,   vorb,       bl(imatrix))
        call matsimtr('matrix','occu','kij')
        call matpose('kij')
        call fafwrite(ndisk,bl(ikij),8,nval*nval,1,ij,istatus)
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine iprint(imat,idimen)
      implicit none
      integer imat(*),idimen,i
      do i=1,idimen
      print *,imat(i)
      enddo
      call flush(6)
      end
c==============================================================================
      subroutine take_strip(a,b,iar,iac,iar1,iar2,iac1,iac2)
      implicit none
      integer iar,iac,iar1,iar2,iac1,iac2
      real*8 a(iar,iac),b(iar2-iar1+1,iac2-iac1+1)
c
      integer ibr,ibc,i,j
      ibr=iar2-iar1+1
      ibc=iac2-iac1+1
      do j=1,ibc
        do i=1,ibr
          b(i,j)=a(iar1+i-1,iac1+j-1)
        enddo
      enddo
      end
