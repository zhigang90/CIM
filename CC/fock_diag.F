      subroutine FockRead(ncf,ifock)
      use memory
      implicit none
      integer ncf,ifock
c
      integer ni,np4,ixx,mataddr
c
      call dynamic_mmark()
      call matmark()
      call dynamic_matdef('work','s',ncf,ncf)
      ixx=mataddr('work')
c     call show_pointer("Work pointer ",bl(ixx))
      call matconn('F','q',ncf,ncf,ifock)
      np4=4
      call sudat(np4,'fock_rhf',ni)
      if(ni.gt.0) then
        call matread('work',np4,'fock_rhf')
      else
        call restart(np4,0)
        call sudat(np4,'fock_rhf',ni)
        if (ni.gt.0) then
          call matread('work',np4,'fock_rhf')
        else
          call nerror(2,'FockBuilder',
     1   'Program cannot find the fock matrix on <jobname>.14',
     2    ni,ni)
        end if
      end if
      call matcopy('work','F')
      call dynamic_retmark()
      call matremark()
      end
C===============================================================================
      subroutine Fock_Diag(ncf,ifock,ioverlap,icano,iepsi,ort_trans,
     *                     ibassym,nmo)
c
c Fock diagonalizer. It does not perform actual Fock matrix
c diagonalization, which is performed by Sub_Fock_Diag called from here,
c but transforms Fock and overlap matrix to symmetry adapted basis, then
c it extracts submatrices (this is facilitated by ibassym array). Then
c Sub_Fock_Diag diagonalizes the submatrices. The small block of
c coefficients is then stored in larger array and back transformation
c is performed at the end. The Fock matrix is purified too, for CCSD.
c This is done by putting back subblocks of Fock matrix into original,
c zeroed Fock matrix. This way off-diagonal parts are explicitly set to
c zero. Both Fock and coefficients are transformed back to AO at the end.
c
c Parameters:
c      Input:
c ncf        - number of contracted basis functions
c ioverlap   - bl index to overlap matrix
c ort_trans  - AO -> Sym. Adap. AO transformation matrix (SAAO)
c ibassym    - an array indicating first and last SAAO for a given irrep.
c nmo        - number of molecular orbitals
c      Both input and output (i.e. quantities are updated):
c ifock      - bl index to Fock matrix
c icano      - bl index canonical orbitals matrix
c iepsi      - bl index to orbital energies
c
c To be precise, the last three variables are not updated, but the
c quantities they point to are updated.
c
      use memory
      implicit none
      integer ncf,ifock,ioverlap,icano,iepsi,ibassym(2,8),nmo
      real*8 ort_trans(ncf,ncf)
c
      integer iort_trans,ifock_o,ioverlap_o,i,icano_o,j,k,ip,jp
      integer istart,istop,isize,isubfock,isuboverl,isubepsi,isubcano
      integer nepsi(ncf),isource,idest,iii,mataddr,inv_ort
      real*8 ep1,ep2,xepsisum2,xepsisum1,xnorm(ncf),xmax,ele
c
      do i=1,ncf
        nepsi(i)=i
      enddo
      xepsisum1=0d0
      do i=0,nmo-1
        xepsisum1=xepsisum1+bl(iepsi+i)
      enddo
c
      call dynamic_mmark()
      call matmark()
c
      call pointer2bl(ort_trans,iort_trans)
      call dynamic_getmem(ncf*ncf,inv_ort)
      call matconn('ortog1','q',ncf,ncf,iort_trans)
      call matconn('inv_ort','q',ncf,ncf,inv_ort)
      do i=1,ncf
        xnorm(i)=0d0
        do j=1,ncf
          xnorm(i)=xnorm(i)+ort_trans(j,i)*ort_trans(j,i)
        enddo
      enddo
      do i=1,ncf
        xnorm(i)=1d0/xnorm(i)
      enddo
      call matcopy('ortog1','inv_ort')
      do i=1,ncf
        do j=1,ncf
          bl(inv_ort+(j-1)+(i-1)*ncf)=ort_trans(i,j)*xnorm(j)
        enddo
      enddo
      call dynamic_getmem(ncf*ncf,ifock_o)
      call matconn('1Fock_O','q',ncf,ncf,ifock_o)
      call dynamic_getmem(ncf*ncf,icano_o)
      call zeroit(bl(icano_o),ncf*ncf)
      call matconn('1FockAO','q',ncf,ncf,ifock)
      call matsimtr('1FockAO','ortog1','1Fock_O')
      call dynamic_getmem(ncf*ncf,ioverlap_o)
      call matconn('over','q',ncf,ncf,ioverlap)
      call matconn('over_o','q',ncf,ncf,ioverlap_o)
      call matsimtr('over','ortog1','over_o')
      call matdisc('over_o')
      call matdisc('over')
      call dynamic_getmem(ncf*ncf,isubfock)
      call dynamic_getmem(ncf*ncf,isuboverl)
      call dynamic_getmem(ncf*ncf,isubcano)
      call dynamic_getmem(ncf,isubepsi)
      do k=1,8
        if (ibassym(1,k).eq.0) cycle
c
        istart=ibassym(1,k)
        istop =ibassym(2,k)
        isize=istop-istart+1
        do i=istart,istop
          ip=i-istart
          do j=istart,istop
            jp=j-istart
            bl(isubfock+ip+jp*isize)=bl(ifock_o+(i-1)+(j-1)*ncf)
            bl(isuboverl+ip+jp*isize)=bl(ioverlap_o+(i-1)+(j-1)*ncf)
          enddo
        enddo
        do i=istart,istop
          call zeroit(bl(ifock_o+(i-1)*ncf),ncf)
        enddo
        call Sub_Fock_Diag(isize,isubfock,isuboverl,isubcano,isubepsi)
        do i=istart,istop
          ip=i-istart
          do j=istart,istop
            jp=j-istart
            bl(icano_o+(i-1)+(j-1)*ncf)=bl(isubcano+ip+jp*isize)
            bl(ifock_o+(i-1)+(j-1)*ncf)=bl(isubfock+ip+jp*isize)
          enddo
          bl(iepsi+i-1)=bl(isubepsi+ip)
        enddo
      enddo
      call matsimtr('1Fock_O','inv_ort','1FockAO')
c Sort Epsilons:
      do i=1,ncf
        do j=2,ncf+1-i
          ep1=bl(iepsi+j-2)
          ep2=bl(iepsi+j-1)
          if (ep1.gt.ep2) then
            bl(iepsi+j-2)=ep2
            bl(iepsi+j-1)=ep1
            k=nepsi(j)
            nepsi(j)=nepsi(j-1)
            nepsi(j-1)=k
          endif
        enddo
      enddo
c
c Reorder orbitals according to epsilons:
c
      do i=1,ncf
        isource=icano_o  +(nepsi(i)-1)*ncf
        idest  =isubcano +(i-1)*ncf
        call tfer(bl(isource),bl(idest),ncf)
      enddo
      xepsisum2=0d0
      do i=0,nmo-1
        xepsisum2=xepsisum2+bl(iepsi+i)
      enddo
      if (dabs((xepsisum2-xepsisum1)/(xepsisum1)).gt.1d-7) then
        write(6,*) dabs(xepsisum2-xepsisum1),xepsisum2,xepsisum1
        call flush(6)
      STOP '43893453587349dfgdf'
      endif
c
      call matconn('cano_O','q',ncf,ncf,isubcano)
      call matconn('cano1','q',ncf,ncf,icano)
      call matmmult('ortog1','cano_O','cano1')
      call dynamic_retmark()
      call matremark()
      end
C===============================================================================
      subroutine Sub_Fock_Diag(ncf,ifock,ioverlap,icano,iepsi)
c Fock diagonalizer. This is more or less a copy of similar algorithm
c from the SCF program.
#define CHOLESKY
      use memory
      implicit none
      integer ncf,ifock,ioverlap,icano,iepsi
      real*8 ort_trans(ncf,ncf)
c
      integer ncf1,ncf2,ioverlap1,ivalues,ifailed,iortog
      integer iworksize1,iworksize2,iworkspace1,iworkspace2
      integer inew_can,i,j,ifockort,ivectors,ne,info,ij,iSpacked
      real*8 xmin,dlamch
c
c
      call dynamic_mmark()
      call matmark()
      ncf1=ncf
      ncf2=ncf
      info=0
      call dynamic_getmem(ncf,ivalues)
      call dynamic_getmem(ncf,ifailed)
      call dynamic_getmem(ncf*ncf,iortog)
      iworksize1=1+6*ncf+2*ncf*ncf
      iworksize2=3+5*ncf
      call dynamic_getmem(iworksize1,iworkspace1)
      call dynamic_getmem(iworksize2,iworkspace2)
      call dynamic_getmem(ncf*ncf,iortog)
c
#ifdef CHOLESKY
c Try Cholesky here:
      call dynamic_getmem(ncf*ncf,iSpacked)
c S in 'U' form:
      ij=-1
      do i=0,ncf-1
        do j=0,i
          ij=ij+1
          bl(iSpacked+ij)=bl(ioverlap+i+j*ncf)
        enddo
      enddo
      info=0
      call dpptrf('U',ncf,bl(iSpacked),info)
      if (info.ne.0) then 
        print *,info
        STOP 'Error in Cholesky'
      endif
      call dtptri('U','N',ncf,bl(iSpacked),info)
      ij=-1
      do i=0,ncf-1
        do j=0,i
          ij=ij+1
          bl(iortog+i+j*ncf)=0d0
          bl(iortog+j+i*ncf)=bl(iSpacked+ij)
        enddo
      enddo
#endif
      xmin=dlamch('S')*2d0
#ifndef CHOLESKY
      call dynamic_getmem(ncf*ncf,ioverlap1)
      call tfer(bl(ioverlap),bl(ioverlap1),ncf*ncf)
      call DSYEVX('V','A','L',ncf1,bl(ioverlap1),ncf,0d0,0d0,0d0,0d0,
     *             xmin,ne,bl(ivalues),bl(iortog),ncf2,bl(iworkspace1),
     *             iworksize1,bl(iworkspace2),bl(ifailed),info)
c
      if (info.ne.0) then
        print *,info
        STOP 'Diagonalization failed'
      endif
      print *,'Ten smallest eigenvalues of S: '
      do i=0,min(10,ncf-1)
        print '(E20.10)',bl(ivalues+i)
      enddo
      call flush(6)
      do i=0,ncf-1
        bl(ivalues+i)=1d0/sqrt(bl(ivalues+i))
      enddo
      do i=0,ncf-1
        do j=0,ncf-1
          bl(iortog+i+j*ncf)=bl(iortog+i+j*ncf)*bl(ivalues+j)
        enddo
      enddo
#endif
      call matconn('ortog','q',ncf,ncf,iortog)
c Transform Fock to ortogonal basis:
      call dynamic_getmem(ncf*ncf,ifockort)
      call matconn('FockOrt','q',ncf,ncf,ifockort)
      call matconn('FockAO','q',ncf,ncf,ifock)
      call matsimtr('FockAO','ortog','FockOrt')
c Diagnonalize ortogonal Fock:
      call dynamic_getmem(ncf*ncf,ivectors)
c     call sdiag2(ncf,ncf1,bl(ifockort),bl(ivalues),bl(ivectors))
      call DSYEVX('V','A','L',ncf1,bl(ifockort),ncf,0d0,0d0,0d0,0d0,
     *            xmin,ne,bl(ivalues),bl(ivectors),ncf2,bl(iworkspace1),
     *            iworksize1,bl(iworkspace2),bl(ifailed),info)
c
      if (info.ne.0) then
        print *,info
        STOP 'Fock diagonalization failed'
      endif
      call matconn('vec','q',ncf,ncf,ivectors)
      call dynamic_getmem(ncf*ncf,inew_can)
      call matconn('new_can','q',ncf,ncf,inew_can)
      call matmmult('ortog','vec','new_can')
      call tfer(bl(inew_can),bl(icano),ncf*ncf)
      call tfer(bl(ivalues),bl(iepsi),ncf)
      call dynamic_retmark()
      call matremark()
      end
C===============================================================================
