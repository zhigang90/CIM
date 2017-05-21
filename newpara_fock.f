      subroutine para_JobInit(iforwhat)

      use newpara

      implicit real*8(a-h,o-z)
C
C  This little routine initializes parallelism on the slaves
C  Formerly this was done when the 2-electron integrals were
C  initialized, but this proved too restrictive as parallelism
C  was extended
C
c
c check if the slaves are OK
c
      call para_check
c
c -- tell slaves what calculations are next
c -- this is SCF (default), MP2 (if iforwhat=5)
c -- CIM input generation (if iforwhat=6) NZG_6/28/2016

      If(iforwhat.EQ.1.OR.iforwhat.EQ.11) Then
         call para_bcast(TxDoScf,TxJobType)
      Else If(iforwhat.EQ.5) Then
         call para_bcast(TxDoMP2,TxJobType)
      Else If(iforwhat.EQ.8) Then
         call para_bcast(TxCIMSubMP2,TxJobType)
      Else
         call para_bcast(TxDoPost,TxJobType)
      EndIf
c
c send the basic data necessary (geometry,cpu,symmetry)
c
      call para_send_info
c
c  send specific data
c
      call para_initsend
      call para_pack_int(iforwhat,1)
c
      call para_bcast_pack(TxJobInit)
C
      RETURN
      END
c=================================================================
c
      subroutine para_oneint(itype,  natom,  oneint, inx,    kk1,
     $                       kk2,    basdat, datnuc, ncs)

      use newpara

      implicit real*8 (a-h,o-z)
c
c  driver for parallel one-electron integral calculation
c
c  The only arguments to this routine that count are
c   oneint -  contains the final one-electron integrals
c             stored 1-dimensionally as an upper triangle
c   ncs    -  number of contracted shells
c
c
      dimension oneint(*),inx(12,ncs),basdat(13,*),datnuc(5,natom)
c
c
c -- zero out the oneint array
      call getival('ncf',ncf)
      ntri = (ncf*(ncf+1))/2
      call zeroit(oneint,ntri)
c
c -- we are going to send the shell indices over to the slaves,
c -- one at a time starting with the highest
c -- send indices on a round-robin basis
c
      do index=ncs,1,-1
      call para_recv(imsg,islv,TxDftReq)
      call para_send(index,islv,TxFTCAssign)
      enddo
c
c -- send termination flag
      do i=1,nslv
      call para_recv(imsg,islv,TxDftReq)
      call para_send(-1,islv,TxFTCAssign)
      enddo
c
c -- gather one-electron matrix from slaves
      call para_reduce(oneint,ntri,TxFTCBTM)
c
c  we are done!
c
      return
      end 
c=================================================================
c
      subroutine para_twoint(iforwhat,ncachex,threshx,threshy,nfock,
     *                       lsemi,   nblocks,idft,   ax,     nrad,
     *                       nang,    lrad,   lang,   IradQ,  NBatch,
     *                       rhf,     NAlpha, NBeta)
c
c=================================================================
c This is the routine that initializes the slaves for SCF
c
      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      Logical rhf
c
      COMMON /DFTCoeff/ aXX(18) ! contains 18 density functional coefficients
c
c2002....................................................
      call getival('iroute',iroute)
      call getival('stab',istab)
c2002....................................................
c
c Obtain and pack slave initialization data
c
      call para_initsend
c
c parameters for twoint and scf
c
      call para_pack_int(nfock,1)
      call para_pack_int(lsemi,1)
      call para_pack_int(iforwhat,1)
      call para_pack_int(ncachex,1)
      call para_pack_int(NAlpha,1)
      call para_pack_int(NBeta,1)
      call para_pack_int(rhf,1)
      call para_pack_real(threshx,1)
      call para_pack_real(threshy,1)
c
c     this is needed by blockint1
c     
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc1)
      if(iftc1.NE.0) then
        call getival('nbl4ftc',nbl4ftc)  ! number of FTC blocks
        nblocks=nbl4ftc
        call getival('ccdd',iccdd)
        call para_pack_int(iccdd,1)
      endif
      call para_pack_int(iroute,1)
      call para_pack_int(istab,1)
c .....................................................
      call para_pack_int(idft,1)
      If(idft.gt.0) Then
        call para_pack_real(aXX,18)
        call para_pack_int(nrad,1)
        call para_pack_int(nang,1)
        call para_pack_int(lrad,1)
        call para_pack_int(lang,1)
        call para_pack_int(IradQ,1)
        call para_pack_int(NBatch,1)
      EndIf
c ....................................................
c Finally, let's send the data
c
      call para_bcast_pack(TxScfInit)
      end
c=================================================================
c
      subroutine para_fock(idft,ax,nblocks,nfock,rhf,ncf,bl,inx,
     $                     thres,dens,fock,fockB,bdn,DenB,
     $                     labels,mywork,igran,mgo,iforwhat)
c
c=================================================================
c
c This routine drives the integral slaves
c constructing the closed-shell Fock matrix .
c
      use newpara

      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension bl(*),dens(*),fock(*),bdn(nfock,*)
      Dimension fockB(*),DenB(*)
      dimension labels(*)
      integer my_mgo
      Logical rhf
c
c check if the slaves are OK
c
      call para_check
      isize=ncf*(ncf+1)/2
      ifsize=nfock*isize
c
c send densities
c
c do not pack as ifsize can be large
c
      call para_bcast_real(dens,isize,TxScfDens)         
      call para_bcast_real(bdn,ifsize,TxScfDens)         
      If(.NOT.rhf) call para_bcast_real(DenB,isize,TxScfDens)  
      my_mgo=mgo
      call para_bcast(my_mgo,TxScfDens)
c
c set the fock matrices to zero
c
      call zeroit(fock,ifsize)
      If(.NOT.rhf) call zeroit(fockB,isize)
c
c Give block assignments to slaves requesting work
c
      call para_distr_blocks(nblocks,igran)      
c
      if(iforwhat.EQ.11) RETURN     ! FTC code
c      
      call para_reduce(fock,ifsize,TxScfFA)
c
      If(.NOT.rhf) call para_reduce(fockB,isize,TxScfFB)
      end
c=================================================================
c
      subroutine slave_oneint(itype,  natom,  oneint,
     $                        inx,    kk1,    kk2,    basdat, datnuc,
     $                        ncs)

      use newpara
      use memory

      implicit real*8(a-h,o-z)
c
c  This routine computes the one-electron integral matrix (H0)
c  on the slaves. It is modelled on subroutine <inton> and has
c  identical arguments, although several of these are redundant
c
c  ARGUMENTS
c
c  itype   -  which type of one-electron integral to calculate
c             there are many options (see subroutine <inton>) but
c             the only one that really counts is H0, for which
c             itype on entry should be 1
c  natom   -  number of real atoms (dummy atoms excluded)
c  oneint  -  linear array for one-electron integrals
c  inx     -  integer basis set information
c  kk1     -  don't really know what these are
c  kk2        usually have value 0, and should be this on entry
c  basdat  -  real basis set information
c  datnuc  -  nuclear data
c  ncs     -  number of contracted shells
c
c
      Dimension oneint(*),basdat(13,*),datnuc(5,natom)
      Integer inx(12,ncs)
c
c     common /big/bl(300)
      common /onethr/onethre
c
c
c -- zero out the oneint array
      call getival('ncf',ncf)      ! should be in argument list
      ntri = (ncf*(ncf+1))/2
      call zeroit(oneint,ntri)
c
c -- set up integral threshold for nuclear potential
c -- it is 1/100th of the main integral threshold (default 10**-10)
      onethre = 1.0d-12
c
c -- reserve memory for the array used to store the general contracted
c -- integrals. This is quite big because we have up to i-functions
c -- (28 components), i.e., 28**2=784, and a maximum of 9 possible
c -- general contractions for each function, i.e., 9*9*784=63504
      call getmem(63504,is)
      call getmem(63504,iss)
c
      k1 = kk1
      k2 = kk2
      if (itype.le.2.or.itype.eq.4.or.itype.eq.6) k1=0
      key = itype
      if (itype.le.1) key=itype+1
      if (itype.eq.9) key=2
c
c -- we are ready to go!
c -- get the shell index from the master
      do 
cc      write(6,*) ' SLAVE:',mytid,' Getting data from Master'
cc      call f_lush(6)
      call para_send(MY_GID,0,TxDftReq)
      call para_recv(ics,ifrom,TxFTCAssign)
cc      write(6,*) ' SLAVE:',mytid,' shell number:',ics
cc      call f_lush(6)
c
c -- check for termination
      If(ics.EQ.-1) exit     ! we are done!
c
c -- calculate the integrals!
      ngci=inx(4,ics)
      iatom=inx(2,ics)
      jfu=0
      len1=inx(3,ics)
      do 50 jcs=1,ics
      ngcj=inx(4,jcs)
      jatom=inx(2,jcs)
      len2=inx(3,jcs)
      len=len1*len2*(ngci+1)*(ngcj+1)
      call onel(ics,jcs,len2,len1,ngcj,ngci,.true.,basdat,
     1          datnuc,bl(is),inx,key,natom,k1,k2)
      if(itype.eq.1 .or. itype.eq.9) then
        do 20 nn=1,natom
        if(itype.eq.9) then
          if(nn.ne.iatom .and. nn.ne.jatom) go to 20
        endif
        call onel(ics,jcs,len2,len1,ngcj,ngci,.true.,basdat,
     1            datnuc,bl(iss),inx,6,nn,k1,k2)
        za=-datnuc(1,nn)
        call add1(bl(iss),za,bl(is),len)
   20   continue
      end if
c .............................
      iij=-1
      ifu=inx(11,ics)-len1
      do 45 igc=0,ngci
      ifu=ifu+len1
      jfu=inx(11,jcs)-len2
      do 45 jgc=0,ngcj
      jfu=jfu+len2
      iff=ifu
      do 40 i1=1,len1
      iff=iff+1
      jff=jfu
      ii=iff*(iff-1)/2
      do 40 j1=1,len2
      jff=jff+1
      iij=iij+1
      if(jff.le.iff) then
        ij=ii+jff
        oneint(ij)=bl(is+iij)
      end if
 40   continue
 45   continue
c .......................
 50   continue
cc      write(6,*) ' SLAVE:',mytid,' Getting another index'
cc      write(6,*) ' Oneint array is:'
cc      ii=0
cc      do i=1,ncf
cc      write(6,90) (oneint(k),k=ii+1,ii+i)
cc      ii=ii+i
cc 90   format (2x,5f14.9)
cc      end do
cc      call f_lush(6)
c
c -- get another shell index
      enddo
c
      call retmem(2)
c
c -- now accumulate final one-electron matrix on master
cc      write(6,*) ' SLAVE:',mytid,' About to accumulate data'
cc      call f_lush(6)
      call para_reduce(oneint,ntri,TxFTCBTM)
cc      write(6,*) ' SLAVE:',mytid,' Just sent data'
cc      call f_lush(6)
c
c  we are done!
c
      return
      end 
