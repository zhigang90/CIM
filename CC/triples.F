      subroutine triples(nval,idimen,npairs,ncf,nmo,vorb,epsi,res,ccsd,
     *                   irecadrx,ndiskx,lbinx,thresh,byt8,qcisd,af,
     *                   nslv)
      use memory
      implicit none
      integer nval,idimen,npairs,ncf,nmo,irecadrx(*),ndiskx,lbinx
      integer nslv
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
      integer :: irecadrtt=0,ndisktt=0,lbintt=0,mem,needmem
      integer i3ext_int,iamplit,iW_scratch,ij,iadr,ia,ic,irecords,n_int
      integer iamplits,ijk,istart,istop,ipass,itot_ijk,npass,ipass_size
      logical :: core,reset1,reset2,semicore
      real*8 ebuilder,eabc,eene,e0,e1,tot_e1,tot_e0,tot_t1,tot_t0
      integer ijk_loc,ivir_start,ivir_stop,ivirpass,nvirpass
      integer ivirpass_size,iscr,nrec_start,idisk_bytes,ii
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
c     semicore=.false.
      endif
      call flush(6)
c 
      if (.not.core.and..not.semicore) then
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
        filename=scrfile(1:len)//'.triples'
        len1=len+8
        ndisk=99
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 99 was used by somebody'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=idimen*idimen*8)
c
c     Calculate how do we have to divide the virtual space
        call dynamic_show_free(mem)
        mem=mem-3*idimen*idimen*nval   ! for amplitudes
        mem=mem-2*idimen*nval            ! 3-int. integrals & copy
        nvirpass=1
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
      itot_ijk=nval*(nval+1)*(nval+2)/6 - nval ! exclude (i,i,i)
      if (.not.core.and..not.semicore) then
c TEST!
        idisk_bytes=200000000000
        ipass_size=idisk_bytes/(48*idimen*idimen*idimen)
        npass=48*idimen*idimen*idimen*itot_ijk/idisk_bytes
        if (mod(48*idimen*idimen*idimen*itot_ijk,idisk_bytes).ne.0) then
          npass=npass+1
        endif
        ivirpass_size=idimen/nvirpass
        if (mod(idimen,nvirpass).ne.0) ivirpass_size=ivirpass_size+1
      else
        nvirpass=1
        npass=1
        ivirpass_size=idimen
      endif
      ii=npass/nslv
      if (mod(npass,nslv).ne.0) ii=ii+1
      npass=nslv*ii
      ipass_size=itot_ijk/npass
      if (mod(itot_ijk,npass).ne.0) ipass_size=ipass_size+1
      print *,'total ijk, ijk pass size, ijk passes: ',
     *         itot_ijk,ipass_size,npass
      call flush(6)
c The record formula for disk:
c (n-1)*6*idimen+(perm-1)*idimen+ic
c n is number of ijk in the present ipass, (ijk starts from 0 in a given pass)
c perm is number of a permutation, ic is a number of W column
      do ipass=1,npass
      istart=(ipass-1)*ipass_size+1
      istop = ipass   *ipass_size
      if (istop.gt.itot_ijk) istop=itot_ijk
      if (istart.gt.itot_ijk) STOP 'Wrong istart in triples'
      print *,'* * * * * * * * * * * IJK pass no.  ',
     *        ipass,',  ijk start, ijk stop: ',istart,istop
      call flush(6)
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
        call dynamic_getmem(idimen*idimen*idimen,iWabc)
        call dynamic_getmem(idimen*idimen,iscr)
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
            call secund(t0)
            call elapsec(e0)
            call Wabc(ndisk,idimen,bl(iscr),bl(iWabc),nrec_start)
            call secund(t1)
            call elapsec(e1)
            abc=abc+t1-t0
            eabc=eabc+e1-e0
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
          enddo
        enddo
      enddo
        call dynamic_retmark()
      endif
c
      enddo ! over ijk passes
      print *,total_energy
      call secund(tot_t1)
      call elapsec(tot_e1)
      print '(A)', '             Build amplitudes,            '//
     *             'Combine,             Energy:'
      print '(A10,3F20.1)','CPU:',builder,abc,ene
      print '(A10,3F20.1)','ELAPS:',ebuilder,eabc,eene
      print '(A)', 'WBuilder      Input,   Output,     Mult: '
      print '(A10,3F10.1)','CPU:',time_i,time_o,time_mult
      print '(A10,3F10.1)','ELAPS:',etime_i,etime_o,etime_mult
      print '(A10,F10.1)','Total triples CPU:  ',tot_t1-tot_t0
      print '(A10,F10.1)','Total triples ELAPS:',tot_e1-tot_e0
      res=total_energy
      call matremark()
      call dynamic_retmark()
      end
C===============================================================================
      subroutine WTildaBuilderWriter(i,j,k,nrec,idimen,npairs,ncf,nmo,
     *                               vorb,ndisk,nval,i3ext_int,iamplit,
     *                               iW_scratch,core,semicore,irecords,
     *                               iamplits,icb,ice)
      use memory
      implicit none
      integer i,j,k,nrec,idimen,npairs,ncf,nmo,ndisk,nval,icb,ice
      integer i3ext_int,iamplit,iW_scratch,irecords(*),iamplits(*)
      logical vorb,core,semicore
      real*8 :: time_i,time_o,time_mult,
     *                       etime_i,etime_o,etime_mult
      common /triples_times/  time_i, time_o, time_mult,
     *                       etime_i,etime_o,etime_mult
c
      integer icoeff,iinteg,iWabc,ic
      integer :: irecadrtt=0,ndisktt=0,lbintt=0
      integer iaddr,iamplit_loc,ijkl,l,irec,mataddr,itmp,ij,il
      real*8 :: thresh = 1d-15,t0,t1,e0,e1
      logical :: byt8 = .true.
      integer :: ii,jj,ijkl_sm,irows,icount=0
c
      icount=icount+1
      call dynamic_mmark()
      call matmark()
c
      call secund(t0)
      call elapsec(e0)
      if (core) then
        if (i.ge.j) then
          ij=i*(i-1)/2+j
        else
          ij=j*(j-1)/2+i
        endif
        iaddr=iamplit+(ij-1)*idimen*idimen
        call matconn('coeff1','q',idimen,idimen,iaddr)
        icoeff=iaddr
        iWabc=iW_scratch
        call matconn('Wabc','r',idimen*idimen,idimen,iWabc)
      else if (semicore) then
        icoeff=iamplits(i)+(j-1)*idimen*idimen
        call matconn('coeff1','q',idimen,idimen,icoeff)
        iWabc=iW_scratch
        call matconn('Wabc','r',idimen*idimen,idimen,iWabc)
      else
        icoeff=iamplits(i)+(j-1)*idimen*idimen
        if (iamplits(i).eq.1) STOP 'Amplitudes for given i not loaded'
        call matconn('coeff1','q',idimen,idimen,icoeff)
        iWabc=iW_scratch
        call matconn('Wabc','r',idimen*idimen,ice-icb+1,iWabc)
      endif
      call elapsec(e1)
      call secund(t1)
      time_i=time_i+t1-t0
      etime_i=etime_i+e1-e0
      do ic=1,idimen
        call secund(t0)
        call elapsec(e0)
        if (core) then
          iaddr=i3ext_int+((ic-1)*nval+k-1)*idimen*idimen
          call matconn('integ1','q',idimen,idimen,iaddr)
        else if (semicore) then
          iaddr=irecords(k)+(ic-1)*idimen*idimen
          call matconn('integ1','q',idimen,idimen,iaddr)
        else
          if (ic.lt.icb) cycle
          if (ic.gt.ice) exit
          if (irecords(k).eq.1) STOP 'Amplitudes for given i not loaded'
          iaddr=irecords(k)+(ic-icb)*idimen*idimen
          call matconn('integ1','q',idimen,idimen,iaddr)
        endif
        call secund(t1)
        call elapsec(e1)
        time_i=time_i+t1-t0
        etime_i=etime_i+e1-e0
        iaddr=iWabc+(ic-icb)*idimen*idimen
        call secund(t0)
        call elapsec(e0)
        call matconn('resul1','q',idimen,idimen,iaddr)
        if (core) then
          if (i.lt.j) then
            call matmmul2('coeff1','integ1','resul1','t','n','a') !!
          else
            call matmmul2('coeff1','integ1','resul1','n','n','a') !!
          endif
        else if (semicore) then
          call matmmul2('coeff1','integ1','resul1','n','n','a') !!!add to ext. W
        else
          call matmmul2('coeff1','integ1','resul1','n','n','n') ! no add here
        endif
        call matdisc('resul1')
        call secund(t1)
        call elapsec(e1)
        time_mult=time_mult+t1-t0
        etime_mult=etime_mult+e1-e0
        call matdisc('integ1')
      enddo
c
c DEBUG:
#if 0
          ic=15
          if (ic.ge.icb.and.ic.le.ice) then
          iaddr=iWabc+(ic-icb)*idimen*idimen
          write(97,*) i,j,k,ic
          do ii=0,idimen*idimen-1
            write(97,'(F20.10)') bl(iaddr+ii)
          enddo
          STOP 'tmp'
          endif
#endif
c
      if (core) then
        call dynamic_matdef('amplit','r',idimen*idimen,nval)
        iamplit_loc = mataddr('amplit')
      endif
      call dynamic_matdef('jkl',   'r',idimen,nval)
      ijkl    = mataddr('jkl')
c Read amplitudes:
      if (.not.core) then
        if (iamplits(i).eq.1) STOP 'Amplitudes for given i not loaded'
        iamplit_loc=iamplits(i)
        call matconn('amplit','r',idimen*idimen,nval,iamplit_loc)
      endif
      call secund(t0)
      call elapsec(e0)
      if (core) then
        do l=1,nval
          iaddr=iamplit_loc+(l-1)*idimen*idimen
          if (i.ge.l) then
            il=i*(i-1)/2+l
          else
            il=l*(l-1)/2+i
          endif
          itmp=iamplit+(il-1)*idimen*idimen
          if (i.ge.l) then
            call tfer(bl(itmp),bl(iaddr),idimen*idimen)
          else
            call copy_with_trans(bl(itmp),bl(iaddr),idimen)
          endif
        enddo
      endif
c
c Read integrals:
      do l=1,nval
        iaddr=ijkl+(l-1)*idimen
        call xkijk(k,l,j,bl(iaddr))
      enddo
      if (.not.core.and..not.semicore) then
        irows=ice-icb+1
        call dynamic_matdef('jkl_sm',   'r',irows,nval)
        ijkl_sm=mataddr('jkl_sm')
        do jj=1,nval
          do ii=1,irows
            bl(ijkl_sm+(ii-1)+(jj-1)*irows)=
     *                               -bl(ijkl+(icb+ii-2)+(jj-1)*idimen)
          enddo
        enddo
      endif
      call secund(t1)
      call elapsec(e1)
      time_i=time_i+t1-t0
      etime_i=etime_i+e1-e0
c Multiply:
      call secund(t0)
      call elapsec(e0)
      if (.not.core.and..not.semicore) then
c     call matscal('jkl_sm',-1d0)
      call matmmul2('amplit','jkl_sm','Wabc','n','t','a')
      else
      call matscal('jkl',-1d0)
      call matmmul2('amplit','jkl','Wabc','n','t','a')
      endif
      call secund(t1)
      call elapsec(e1)
      time_mult=time_mult+t1-t0
      etime_mult=etime_mult+e1-e0
c Write results:
      call secund(t0)
      call elapsec(e0)
      if (.not.(core.or.semicore)) then
        do ic=icb,ice
          irec=nrec+ic
          iaddr=iWabc+(ic-icb)*idimen*idimen
          call s_write(ndisk,irec,bl(iaddr))
        enddo
      endif
      call secund(t1)
      call elapsec(e1)
      time_o=time_o+t1-t0
      etime_o=etime_o+e1-e0
c
      call matremark()
      call dynamic_retmark()
      end
C===============================================================================
      subroutine Wabc(ndisk,idimen,xmat,W,nrec_start,istart,istop)
      use memory
      implicit none
      integer idimen,ndisk,nrec_start,istart,istop
      real*8 W(idimen,idimen,istop-istart+1),xmat(idimen,idimen)
c
      integer ia,ib,ic,irec,nrec
c (ijk)(abc):
      nrec=nrec_start+0
      do ic=istart,istop
        irec=nrec+ic
        call s_read(ndisk,irec,xmat)
        do ia=1,idimen
          do ib=1,idimen
          W(ia,ib,ic-istart+1)=xmat(ia,ib)
          enddo
        enddo
      enddo
c (ikj)(acb):
      nrec=nrec_start+1*idimen
      do ib=1,idimen
        irec=nrec+ib
        call s_read(ndisk,irec,xmat)
        do ia=1,idimen
          do ic=istart,istop
          W(ia,ib,ic-istart+1)=W(ia,ib,ic-istart+1)+xmat(ia,ic)
          enddo
        enddo
      enddo
c (jik)(bac):
      nrec=nrec_start+5*idimen
      do ic=istart,istop
        irec=nrec+ic
        call s_read(ndisk,irec,xmat)
        do ib=1,idimen
          do ia=1,idimen
          W(ia,ib,ic-istart+1)=W(ia,ib,ic-istart+1)+xmat(ib,ia)
          enddo
        enddo
      enddo
c (jki)(bca):
      nrec=nrec_start+4*idimen
      do ia=1,idimen
        irec=nrec+ia
        call s_read(ndisk,irec,xmat)
        do ib=1,idimen
          do ic=istart,istop
          W(ia,ib,ic-istart+1)=W(ia,ib,ic-istart+1)+xmat(ib,ic)
          enddo
        enddo
      enddo
c (kij)(cab):
      nrec=nrec_start+2*idimen
      do ib=1,idimen
        irec=nrec+ib
        call s_read(ndisk,irec,xmat)
        do ic=istart,istop
          do ia=1,idimen
          W(ia,ib,ic-istart+1)=W(ia,ib,ic-istart+1)+xmat(ic,ia)
          enddo
        enddo
      enddo
c (kji)(cba):
      nrec=nrec_start+3*idimen
      do ia=1,idimen
        irec=nrec+ia
        call s_read(ndisk,irec,xmat)
        do ic=istart,istop
          do ib=1,idimen
          W(ia,ib,ic-istart+1)=W(ia,ib,ic-istart+1)+xmat(ic,ib)
          enddo
        enddo
      enddo
      end
C===============================================================================
      subroutine EnergyContribution(i,j,k,idimen,nval,W,epsi,energy,
     *                              ccsd,irecadrx,npairs,ndiskx,ncf,
     *                              lbinx,thresh,byt8,nmo,vorb,
     *                              qcisd)
      use memory
      implicit none
      integer i,j,k,idimen,nval,irecadrx(*),npairs,ndiskx,ncf,lbinx,nmo
      logical ccsd,byt8,vorb,qcisd
      real*8 W(idimen,idimen,idimen),epsi(*) ! epsi is assumed as starting
c                                              just from above frozen orbitals
      real*8 thresh
      real*8 energy
c
      real*8 occupied,denom,delta,singles
      integer ir,is,it,k_ij,k_ik,k_jk,isingles
c
      call dynamic_mmark()
      if (ccsd.or.qcisd) then
        call dynamic_getmem(idimen*idimen,k_ij)
        call dynamic_getmem(idimen*idimen,k_ik)
        call dynamic_getmem(idimen*idimen,k_jk)
        call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *               ncf,  lbinx, thresh,     byt8,      'x',
     *               'mo', nmo,   vorb,       bl(k_ij))
        call ExtrOne(i,    k,     irecadrx,   npairs,    ndiskx,
     *               ncf,  lbinx, thresh,     byt8,      'x',
     *               'mo', nmo,   vorb,       bl(k_ik))
        call ExtrOne(j,    k,     irecadrx,   npairs,    ndiskx,
     *               ncf,  lbinx, thresh,     byt8,      'x',
     *               'mo', nmo,   vorb,       bl(k_jk))
      call pointersingles(isingles)
      endif
      occupied=epsi(i)+epsi(j)+epsi(k)
      energy=0d0
      do ir=1,idimen
        do is=1,idimen
          do it=1,idimen
            if (ccsd.or.qcisd) then
            singles=
     *    bl(k_ij+(ir-1)+(is-1)*idimen)*bl(isingles+(it-1)+(k-1)*idimen)
     *   +bl(k_ik+(ir-1)+(it-1)*idimen)*bl(isingles+(is-1)+(j-1)*idimen)
     *   +bl(k_jk+(is-1)+(it-1)*idimen)*bl(isingles+(ir-1)+(i-1)*idimen)
            if (qcisd) singles=singles*2d0
            else
              singles=0d0
            endif
            denom=occupied-epsi(nval+ir)-epsi(nval+is)-epsi(nval+it)
            energy=energy+(W(ir,is,it)+singles)*
     *        (
     *        +4d0*W(ir,is,it)
     *        -2d0*W(ir,it,is)
     *        -2d0*W(it,is,ir)
     *        -2d0*W(is,ir,it)
     *        +    W(it,ir,is)
     *        +    W(is,it,ir)        )/denom
          enddo
        enddo
      enddo
      energy=energy*(2d0-delta(i,j)-delta(j,k))
      print '(F20.15)',energy
      call dynamic_retmark()
      end
C===============================================================================
      subroutine copy_with_trans(xsource,xtarget,idimen)
      implicit none
      integer i,j,idimen
      real*8 xsource(idimen,idimen),xtarget(idimen,idimen)
      do i=1,idimen
        do j=1,idimen
          xtarget(i,j)=xsource(j,i)
        enddo
      enddo
      end
C===============================================================================
      subroutine EnergyContrPartial(i,j,k,idimen,nval,W,epsi,energy,
     *                              ccsd,irecadrx,npairs,ndiskx,ncf,
     *                              lbinx,thresh,byt8,nmo,vorb,
     *                              qcisd,istart,istop,ndisk,xmat,
     *                              xk_ij,xk_ik,xk_jk,xsingl)
      use memory
      implicit none
      integer i,j,k,idimen,nval,irecadrx(*),npairs,ndiskx,ncf,lbinx,nmo
      integer istart,istop,ndisk
      logical ccsd,byt8,vorb,qcisd
      real*8 W(idimen,idimen,idimen),epsi(*) ! epsi is assumed as starting
c                                              just from above frozen orbitals
      real*8 xmat(idimen,idimen),xk_ij(idimen,idimen)
      real*8 xk_ik(idimen,idimen),xk_jk(idimen,idimen)
      real*8 xsingl(idimen,nval)
c
      real*8 thresh
      real*8 energy
c
      real*8 occupied,denom,delta,singles
      integer ia,ib,ic
c
      occupied=epsi(i)+epsi(j)+epsi(k)
      energy=0d0
      do ic=istart,istop
        do ib=1,idimen
          do ia=1,idimen
            W(ia,ib,ic-istart+1)= 0d0
          enddo
        enddo
      enddo
      do ic=1,idimen
        call s_read(ndisk,ic,xmat)
        if (ic.ge.istart.and.ic.le.istop) then
          do ib=1,idimen
            do ia=1,idimen
              W(ia,ib,ic-istart+1)= W(ia,ib,ic-istart+1)+4d0*xmat(ia,ib)
            enddo
          enddo
        endif
        do ib=1,idimen
          if (ib.lt.istart) cycle
          if (ib.gt.istop)  exit
          do ia=1,idimen
            W(ia,ic,ib-istart+1)=W(ia,ic,ib-istart+1)-2d0*xmat(ia,ib)
          enddo
        enddo
        do ia=1,idimen
          if (ia.lt.istart) cycle
          if (ia.gt.istop)  exit
          do ib=1,idimen
            W(ic,ib,ia-istart+1)=W(ic,ib,ia-istart+1)-2d0*xmat(ia,ib)
          enddo
        enddo
        if (ic.ge.istart.and.ic.le.istop) then
          do ib=1,idimen
            do ia=1,idimen
              W(ib,ia,ic-istart+1)=W(ib,ia,ic-istart+1)-2d0*xmat(ia,ib)
            enddo
          enddo
        endif
        do ia=1,idimen
          if (ia.lt.istart) cycle
          if (ia.gt.istop)  exit
          do ib=1,idimen
            W(ib,ic,ia-istart+1)=W(ib,ic,ia-istart+1)+xmat(ia,ib)
          enddo
        enddo
        do ib=1,idimen
          if (ib.lt.istart) cycle
          if (ib.gt.istop)  exit
          do ia=1,idimen
            W(ic,ia,ib-istart+1)=W(ic,ia,ib-istart+1)+xmat(ia,ib)
          enddo
        enddo
      enddo
c
      do ic=istart,istop
        call s_read(ndisk,ic,xmat)
        do ia=1,idimen
          do ib=1,idimen
            if (ccsd.or.qcisd) then
            singles=     xk_ij(ia,ib)*xsingl(ic,k)
     *                  +xk_ik(ia,ic)*xsingl(ib,j)
     *                  +xk_jk(ib,ic)*xsingl(ia,i)
              if (qcisd) singles=singles*2d0
            else
              singles=0d0
            endif
            denom=occupied-epsi(nval+ia)-epsi(nval+ib)-epsi(nval+ic)
            energy=energy+(xmat(ia,ib)+singles)*W(ia,ib,ic-istart+1)
     *                     /denom
          enddo
        enddo
      enddo
      energy=energy*(2d0-delta(i,j)-delta(j,k))
      print '(F20.15)',energy
      end
C===============================================================================
      subroutine Wabc_core(idimen,Wsourc,Wdest,typ)
      use memory
      implicit none
      integer idimen
      real*8 Wdest(idimen,idimen,idimen)
      real*8 Wsourc(idimen,idimen,idimen)
      character*3 typ
c
      integer ia,ib,ic
c (ijk)(abc):
      if (typ.eq.'ijk') then ! Warning! special case, W is also initialized!
      do ia=1,idimen
        do ib=1,idimen
          do ic=1,idimen
          Wdest(ia,ib,ic)=Wsourc(ia,ib,ic)
          enddo
        enddo
      enddo
      return
      endif
c (ikj)(acb):
      if (typ.eq.'ikj') then
      do ia=1,idimen
        do ib=1,idimen
          do ic=1,idimen
          Wdest(ia,ib,ic)=Wdest(ia,ib,ic)+Wsourc(ia,ic,ib)
          enddo
        enddo
      enddo
      return
      endif
c (jik)(bac):
      if (typ.eq.'jik') then
      do ia=1,idimen
        do ib=1,idimen
          do ic=1,idimen
          Wdest(ia,ib,ic)=Wdest(ia,ib,ic)+Wsourc(ib,ia,ic)
          enddo
        enddo
      enddo
      return
      endif
c (jki)(bca):
      if (typ.eq.'jki') then
      do ia=1,idimen
        do ib=1,idimen
          do ic=1,idimen
          Wdest(ia,ib,ic)=Wdest(ia,ib,ic)+Wsourc(ib,ic,ia)
          enddo
        enddo
      enddo
      return
      endif
c (kij)(cab):
      if (typ.eq.'kij') then
      do ia=1,idimen
        do ib=1,idimen
          do ic=1,idimen
          Wdest(ia,ib,ic)=Wdest(ia,ib,ic)+Wsourc(ic,ia,ib)
          enddo
        enddo
      enddo
      return
      endif
c (kji)(cba):
      if (typ.eq.'kji') then
      do ia=1,idimen
        do ib=1,idimen
          do ic=1,idimen
          Wdest(ia,ib,ic)=Wdest(ia,ib,ic)+Wsourc(ic,ib,ia)
          enddo
        enddo
      enddo
      return
      endif
      STOP 'Error in Wabc_core()'
      end
C===============================================================================
      subroutine sorter(W,idimen,ione,itwo)
      implicit none
      integer idimen,ione,itwo
      real*8 W(idimen,idimen,idimen)
      real*8 xtmp
      integer ia,ib,ic
      if (ione.eq.1.and.itwo.eq.2) then
        do ic=1,idimen
          do ia=1,idimen
            do ib=1,ia-1
              xtmp=W(ia,ib,ic)
              W(ia,ib,ic)=W(ib,ia,ic)
              W(ib,ia,ic)=xtmp
            enddo
          enddo
        enddo
        return
      endif
      if (ione.eq.2.and.itwo.eq.3) then
        do ic=1,idimen
          do ib=1,ic-1
            do ia=1,idimen
              xtmp=W(ia,ib,ic)
              W(ia,ib,ic)=W(ia,ic,ib)
              W(ia,ic,ib)=xtmp
            enddo
          enddo
        enddo
        return
      endif
      if (ione.eq.1.and.itwo.eq.3) then
        do ib=1,idimen
          do ic=1,idimen
            do ia=1,ic-1
              xtmp=W(ia,ib,ic)
              W(ia,ib,ic)=W(ic,ib,ia)
              W(ic,ib,ia)=xtmp
            enddo
          enddo
        enddo
        return
      endif
      STOP 'Error in sorter'
      end
C===============================================================================
      subroutine locate_integrals(i,    j, k,     i3ext_int,  irecords,
     *                           reset,n, idimen,irecadrtt,npairs,
     *                           ndisktt,ncf,lbintt,thresh,byt8,nmo,
     *                           vorb,nval,icb,ice)
      use memory
      implicit none
      integer i,j,k,i3ext_int,irecords(*),n,idimen,irecadrtt,npairs
      integer ndisktt,ncf,lbintt,nmo,nval,icb,ice
      real*8 thresh
      logical reset,byt8,vorb
c
      integer ia,iaddr,iaddr1,ic,ind1,ind2,ind3,irep
      integer addresses(n),iad_max
c
      if (reset) then
        reset=.false.
        do ia=1,n
          iaddr        =i3ext_int+(ia-1)*idimen*idimen*(ice-icb+1)
          addresses(ia)=iaddr
          if (iaddr.eq.1) 
     *  STOP 'Cannot be! I use assumption that the address = 1 is null!'
        enddo
        iad_max=n
        do ia=1,nval
          irecords(ia)=1
        enddo
        if (i.gt.n.or.j.gt.n.or.k.gt.n) then ! load i,j,k first
          do irep=1,3
          if (irep.eq.1) ind1=i
          if (irep.eq.2) ind1=j
          if (irep.eq.3) ind1=k
          if (irecords(ind1).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ind1)=iaddr
          do ic=icb,ice
            iaddr1=iaddr+(ic-icb)*idimen*idimen
            call ExtrOne(ic,    ind1,    irecadrtt,  npairs,    ndisktt,
     *                   ncf,  lbintt,thresh,     byt8,      'tt',
     *                   'mo', nmo,   vorb,       bl(iaddr1))
          enddo
          enddo
        endif
        do ia=1,n
          if (iad_max.le.0) exit
          if (irecords(ia).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ia)=iaddr
          do ic=icb,ice
            iaddr1=iaddr+(ic-icb)*idimen*idimen
            call ExtrOne(ic,    ia,     irecadrtt,  npairs,    ndisktt,
     *                   ncf,  lbintt,thresh,     byt8,      'tt',
     *                   'mo', nmo,   vorb,       bl(iaddr1))
          enddo
        enddo
      endif
c
      do irep=1,3
      if (irep.eq.1) then 
       ind1=i
       ind2=j
       ind3=k
      endif
      if (irep.eq.2) then 
       ind1=j
       ind2=i
       ind3=k
      endif
      if (irep.eq.3) then 
       ind1=k
       ind2=i
       ind3=j
      endif
c
      if (irecords(i).ne.1.and.irecords(j).ne.1.and.irecords(k).ne.1)
     *                                                              then
c DEBUG:
          do ia=1,nval
            if (irecords(ia).ne.1) write(*,'(I3)',ADVANCE='NO') ia
          enddo
          write(*,*)
c END DEBUG
          return
      endif
c
      if (irecords(ind1).eq.1) then ! I have to think about it
        do ia=nval,1,-1 ! search for largest record which can be replaced
          if (irecords(ia).ne.1) then
            if (ia.eq.ind1) STOP 'Cannot be: ia.eq.i'
            if (ia.eq.ind2.or.ia.eq.ind3) cycle
            iaddr=irecords(ia)
            irecords(ia)=1
            if (irecords(ind1).ne.1) STOP 'irecords(i).ne.1'
            irecords(ind1)=iaddr
            do ic=icb,ice
              iaddr1=iaddr+(ic-icb)*idimen*idimen
              call ExtrOne(ic,   ind1,   irecadrtt,  npairs,    ndisktt,
     *                     ncf,  lbintt,thresh,     byt8,      'tt',
     *                     'mo', nmo,   vorb,       bl(iaddr1))
            enddo
            exit
          endif
        enddo
      endif
c
      enddo
c
      do ia=1,nval
        if (irecords(ia).ne.1) write(*,'(I3)',ADVANCE='NO') ia
      enddo
      write(*,*)
c
      if (irecords(i).eq.1.or.irecords(j).eq.1.or.irecords(k).eq.1) then
c       write(*,'(F10.x)')
        STOP 'Cannot leave this procedure (int) with such condition .tr'
      endif
      end
C===============================================================================
      subroutine locate_amplits(i,    j, k,     i3ext_int,  irecords,
     *                           reset,n, idimen,irecadrtt,npairs,
     *                           ndisktt,ncf,lbintt,thresh,byt8,nmo,
     *                           vorb,nval)
      use memory
      implicit none
      integer i,j,k,i3ext_int,irecords(*),n,idimen,irecadrtt,npairs
      integer ndisktt,ncf,lbintt,nmo,nval
      real*8 thresh
      logical reset,byt8,vorb
c
      integer ia,iaddr,iaddr1,ib,ind1,ind2,ind3,irep,itmp
      integer addresses(n),iad_max
c
      if (reset) then
        reset=.false.
        do ia=1,n
          iaddr        =i3ext_int+(ia-1)*nval*idimen*idimen
          addresses(ia)=iaddr
          if (iaddr.eq.1)
     *  STOP 'Cannot be! I use assumption that the address = 1 is null!'
        enddo
        iad_max=n
        do ia=1,nval
          irecords(ia)=1
        enddo
        if (i.gt.n.or.j.gt.n.or.k.gt.n) then
          do irep=1,3
          if (irep.eq.1) ind1=i
          if (irep.eq.2) ind1=j
          if (irep.eq.3) ind1=k
          if (irecords(ind1).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ind1)=iaddr
          do ib=1,nval
            iaddr1=iaddr+(ib-1)*idimen*idimen
            call CoefRead('tt',ind1,ib,itmp)
            call tfer(bl(itmp),bl(iaddr1),idimen*idimen)
          enddo
          enddo
        endif
        do ia=1,n
          if (iad_max.le.0) exit
          if (irecords(ia).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ia)=iaddr
          do ib=1,nval
            iaddr1=iaddr+(ib-1)*idimen*idimen
            call CoefRead('tt',ia,ib,itmp)
            call tfer(bl(itmp),bl(iaddr1),idimen*idimen)
          enddo
        enddo
      endif
c
      do irep=1,3
      if (irep.eq.1) then 
       ind1=i
       ind2=j
       ind3=k
      endif
      if (irep.eq.2) then 
       ind1=j
       ind2=i
       ind3=k
      endif
      if (irep.eq.3) then 
       ind1=k
       ind2=i
       ind3=j
      endif
c
      if (irecords(i).ne.1.and.irecords(j).ne.1.and.irecords(k).ne.1)
     *                                                              then
          return
      endif
c
      if (irecords(ind1).eq.1) then
        do ia=nval,1,-1 ! search for largest record which can be replaced
          if (irecords(ia).ne.1) then
            if (ia.eq.ind1) STOP 'Cannot be: ia.eq.i'
            if (ia.eq.ind2.or.ia.eq.ind3) cycle
            iaddr=irecords(ia)
            irecords(ia)=1
            if (irecords(ind1).ne.1) STOP 'irecords(i).ne.1'
            irecords(ind1)=iaddr
            do ib=1,nval
              iaddr1=iaddr+(ib-1)*idimen*idimen
              call CoefRead('tt',ind1,ib,itmp)
              call tfer(bl(itmp),bl(iaddr1),idimen*idimen)
            enddo
            exit
          endif
        enddo
      endif
c
      enddo
c
      if (irecords(i).eq.1.or.irecords(j).eq.1.or.irecords(k).eq.1) 
     *  STOP 'Cannot leave this procedure with such condition .true.'
      end
