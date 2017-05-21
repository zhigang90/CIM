c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C The real PARA_FOCK routine was MOVED!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
C July 4,97 KW: number of generators for a symmetry group is
C               taken from options (getival('ngener',ngener)
C               Symmetrization of the Fock matrix is made bY :
C               call FockSymm_SCF(ngener,ncf,nfock,fock,inxx(ifp))
C    end of this routine
c---------------------------------------------------------------------
      subroutine pfock(idft,   ax,     nblocks,nfock,  rhf,
     $                 ncf,    bl,     inx,    thres,  mgo,
     $                 dens,   fock,   fockB,  bdn,    DenB,
     $                 labels, iforwhat)
c---------------------------------------------------------------------
c This routine is just a wrapper around para_fock
c calling the parallel integral calculation or int_fock
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /datstore/ thresh_stored,isto
c     common /intbl/ifpp,inxx(100)
      dimension inx(12,*)
      dimension bl(*),dens(*),fock(nfock,*),bdn(nfock,*)
      dimension fockB(*),DenB(*)
      dimension labels(*)
      integer igran
      logical rhf
c
      call getival('gran',igran)
c
c     if(scftype.ne.'full-direct') then
c        write(6,*)' cannot run job as ',scftype
c        stop 'cannot do it'
c     endif
c
c     if there is only one cpu then don't bother with PVM
c
c
c set flag for retrieval of stored integrals
c
      isto=1
      mywork=0
      call para_fock(idft,ax,nblocks,nfock,rhf,ncf,bl,inx,
     *               thres,dens,fock,fockB,bdn,DenB,
     *               labels,mywork,igran,mgo,iforwhat)
c
      If(iforwhat.EQ.11) RETURN     ! FTC code
c
c  now divide the non-diagonal elements of the Fock matrix by 2
c  and scale back to the normal values from the large ones
c
      If(rhf) call resc_nfock(ncf,nfock,fock,thres)
      If(.not.rhf) call resu_nfock(ncf,fock,fockB,thres)
c
c  symmetrize Fock matrix for exact symmetry
c
      call getival('nsym',nsym)
c
      if(nsym.eq.0) return
c
      call getival('ngener',ngener)
      call getival('SymFunPr',ifp)
c
      call FockSymm_SCF(ngener,ncf,nfock,fock,bl(ifp))
      If(.not.rhf) call FockSymm_SCF(ngener,ncf,nfock,fockB,bl(ifp))
c
      end
