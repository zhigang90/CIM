C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C THE REAL PARA_CPHF ROUTINE WAS MOVED!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C July 07,97 KW: symmetrization of the F(D1,g0) matrices has been added:
c
c        call focksymm_nmr(ngener,ncf,fock,inxx(ifp),inxx(nsyo))
c
c        end of this routine ; line 100
c=====================================================================
c
      subroutine pcphf(idft,ax,nblocks,bl,inx,ntri,thres,mgo,
     *                 dens,fock,labels)
c
c=====================================================================
c This routine calls two-el. int. block by block
c and constructes closed-shell Fock matrix .
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      common /lindvec/ lind,idensp
      common /datstore/ thresh_stored,isto
      dimension inx(12,*)
      dimension bl(*),dens(*),fock(*)
      dimension labels(*)
      parameter (Zero=0.0d0)
c----------------------------------------------------------------
c make and get screening densities :
c
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getmem(2*ncs*ncs,idscree)
c
      call setup_screen_dens(inx,ncf,ntri,ncs,dens,bl(idscree),bl)
c----------------------------------------------------------------
      call getival('gran',igran)
      call zeroit(fock,ntri*3)
c setup flag for retrieval of stored integrals
      isto=1
      myjob=0
      call para_cphf(nblocks,bl,inx,ntri,thres,
     *               bl(idscree), dens,fock,
     *               labels,myjob,igran,mgo)
c
c----------------------------------------------------------------
      call retmem(1) ! idscree
c----------------------------------------------------------------
c
c introduce scaling factor in case of dft
c
      thresax=thres
      if(idft.gt.0.and.ax.NE.Zero) thresax=thresax*ax
c
c re-scale fock matrices :
c
      call rescale_fock(fock,ntri,bl(lind),thresax)
c
c symmetrize fock matrices for exact symmetry :
c
      call getival('nsym',nsym)
      if(nsym.gt.0) then
         call getival('ngener',ngener)
         call getival('nsyo',nsyo)
         call getival('SymFunPr',ifp)
         call getival('ncf',ncf)
         call focksymm_nmr(ngener,ncf,fock,bl(ifp),bl(nsyo))
      endif
      end
c=====================================================================
      subroutine setup_screen_dens(inx,ncf,ntri,ncs,dens,dscreen,bl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension bl(*),dens(*)
      dimension dscreen(ncs*ncs,2)
c----------------------------------------------------------------
c screening density (ites)
c----------------------------------------------------------------
c there are 3 density matrices; select one (max) for calcint2 call
c
      call getmem(ntri,lselect)
c
      call select_dens(dens,ntri, bl(lselect))
c
c transform selected density bl(lselect) dens(ij) into dscreen(ics,jcs)
c
      call getmem(ncf,map_fs)
c
c
      call setup_densp2(inx,ncf,ncs,bl(lselect),dscreen(1,1),bl(map_fs))
c
c      densp is NOT SQUARED
c
      call retmem(1)    ! map_fs
      call retmem(1)    ! lselect
c----------------------------------------------------------------
c read-in second screening density from a disk
c
         call read1mat(69, 1 ,ncs*ncs,dscreen(1,2))
c----------------------------------------------------------------
      end
c=====================================================================
