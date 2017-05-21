      subroutine get_aat(natoms,nat_b,nat_e,npos,ncf,nocc,ntri,
     *                   listunq,lww,aat)
c-----------------------------------------------------------------------
c This routine calculates the GIAO atomic axial tensors
c for the atoms in this CPHF cycle      
c
c D1 = SUM(jocc,avir){Cj(T)*F(D1,g0)*Ca/(ej-ea-xlv)*[Cj*Ca(T)+Ca*Cj(T)]}
c
c  (T) = transpose
c-----------------------------------------------------------------------
c  INTENT(IN)
c  nat_b, nat_e  = first and last atom
c  npos          = offset in the atom list      
c  ncf           = number of contracted basis functions
c  nocc          = number of occupied orbitals
c  listunq       = pointer to atom position in full list      
c  STORAGE
c  lww 3*ncf*ncf
c  INTENT(OUT)
c  aat           = atomic axial tensors (for the symm. unique atoms)
c---------------------------------------------------
      use memory
      implicit real*8 (a-h,o-z)
      dimension w0(ncf*nocc,3)
      dimension aat(3,3,natoms),aaxt(3,3),listunq(*)
c
      ncoefile=79
      ncfmx=ncf*nocc
      call zeroit(aat,9*natoms)
      call getival('printh',nprint)    !print level
      call getival('nogiao',nogiao)    !giao basis
      call getival('iout',iout)
      call getival('na',natoms)
      call getival('ncs ',ncs)
      call getival('ictr',ictr)
      call getival('inuc',inuc)
      call getival('ibas',ibas)
      call getival('lsmat',lsmat)   ! overlap mx S0
      call getival('ldensi',ldens)   ! dens. mx
      call getival('lvec',lvec)   ! MO coefficients
      call matmark
      call mmark
      call matdef('cc','r',ncf,nocc)
      call matconn('s','s',ncf,ncf,lsmat)
      call matconn('c0','r',ncf,nocc,lvec)
      do iat=nat_b,nat_e
         iatr=listunq(iat)
c First contribution:
c Tr(Cgeom(T)*S0*Cmagn)
c     
         call getmem(ncfmx*3,lcm)  !bl(lcm)=dC/dX
         call read1mat(ncoefile,iat-npos,ncfmx*3,bl(lcm))
         call rea(bl(lww),ncfmx*3,2,'c1_magn')  !w0= dC/dB
         do i=1,3
            call matconn('c1_coo','r',ncf,nocc,lcm+(i-1)*ncfmx)
            do j=1,3
              call matconn('c1_magn','r',ncf,nocc,lww+(j-1)*ncfmx)
              call matmmul2('s','c1_magn','cc','n','n','n')
              call matprodtr('c1_coo','cc',aat(j,i,iatr))
              call matdisc('c1_magn')
            enddo
            call matdisc('c1_coo')
         enddo
         if(nprint.ge.2) write(iout,10) 'CmCx',iatr,
     &                         ((aat(j,i,iatr),j=1,3),i=1,3)
c second contribution
c Tr(Cgeom(T)*Smagn*C0)
         if(nogiao.eq.0) then
         call vcdint12(natoms,bl(ictr),bl(ibas),bl(inuc),ncs,ncf,2,0,
     &                 bl(lww))
         call dscal(ncfmx*3,-2.0d0,bl(lcm),1)
         do i=1,3
            call matconn('c1_coo','r',ncf,nocc,lcm+(i-1)*ncfmx)
            do j=1,3
               call matconn('sder','q',ncf,ncf,lww+(j-1)*ncf*ncf)
               call matmmul2('sder','c0','cc','n','n','n')
               call matprodtr('c1_coo','cc',aaxt(j,i))
               aat(j,i,iatr)=aat(j,i,iatr)+aaxt(j,i)
               call matdisc('sder')
            enddo
            call matdisc('c1_coo')
         enddo
         if(nprint.ge.2) write(iout,10) 'CCx',iatr,aaxt
         endif
c third contribution
c Tr(Cmagn*C0(T)*Sgeom)         
        call vcdint12(natoms,bl(ictr),bl(ibas),bl(inuc),ncs,ncf,1,iatr,
     &                bl(lww))
         call rea(bl(lcm),ncfmx*3,2,'c1_magn')  
         do i=1,3
            call matconn('sder','q',ncf,ncf,lww+(i-1)*ncf*ncf)
            do j=1,3
               call matconn('c1_magn','r',ncf,nocc,lcm+(j-1)*ncfmx)
               call matmmul2('sder','c1_magn','cc','n','n','n')
               call matprodtr('c0','cc',aaxt(j,i))
               call matdisc('c1_magn')
               aat(j,i,iatr)=aat(j,i,iatr)+aaxt(j,i)
            enddo
            call matdisc('sder')
         enddo
         if(nprint.ge.2) write(iout,10) 'CmC',iatr,aaxt
         call retmem(1)
c fourth contribution         
         if(nogiao.eq.0) then
         call vcdint(natoms,bl(ictr),bl(ibas),bl(inuc),ncs,ncf,ntri,
     &                bl(ldens),iatr,aaxt)
         do i=1,3
         do j=1,3
               aat(j,i,iatr)=aat(j,i,iatr)+aaxt(j,i)
         enddo
         enddo
         if(nprint.ge.2) write(iout,10) 'CC',iatr,aaxt
         endif
c nuclear aat         
         call aaxnuc(bl(inuc),natoms,iatr,aaxt)
         do i=1,3
         do j=1,3
               aat(j,i,iatr)=aat(j,i,iatr)+aaxt(j,i)
         enddo
         enddo
         if(nprint.ge.2) write(iout,10) 'nucl',iatr,aaxt
         if(nprint.ge.1) write(iout,10) 'Atomic Axial Tensor for atom',
     &                               iatr,((aat(j,i,iatr),j=1,3),i=1,3)
         
      enddo
      call retmark
      call matremark
 10   format(A,I3,/,(3(1F8.4)))
      end
c======================================================================
      subroutine aaxnuc(xnuc,natoms,iat,aaxt)
      implicit real*8 (a-h,o-z)
      parameter (fourth=0.25d0)
      dimension xnuc(5,natoms),aaxt(3,3)
c
c  calculates nuclear part of atomic axial tensor for atom iat
c  Arguments
c  INTENT(IN)
c  xnuc     = nuclear info
c  natoms   = no of atoms
c  iat      = no of atom used
c  INTENT(OUT)
c  aaxt     = nucl. a. a. tensor
      aaxt(1,1)=0.0d0
      aaxt(2,1)=fourth*xnuc(1,iat)*xnuc(4,iat)
      aaxt(3,1)=-fourth*xnuc(1,iat)*xnuc(3,iat)
      aaxt(1,2)=-fourth*xnuc(1,iat)*xnuc(4,iat)
      aaxt(2,2)=0.0d0
      aaxt(3,2)=fourth*xnuc(1,iat)*xnuc(2,iat)
      aaxt(1,3)=fourth*xnuc(1,iat)*xnuc(3,iat)
      aaxt(2,3)=-fourth*xnuc(1,iat)*xnuc(2,iat)
      aaxt(3,3)=0.0d0
      end
c ========================================================================
c
      SUBROUTINE SymAAT(NAtoms, NTrans, ISYM,   NEqATM, TRANS,
     $                  TRANSA, AAT,    V,      Idon)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Generates atomic axial tensors for all atoms
C  using the results for symmetry (full point group) unique atoms
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  NTrans  -  number of symmetry operations
C  ISYM    -  indicates which atoms are unique 
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  TRANS   -  transformation matrices for normal symmetry operations
C  TRANSA  -    ditto for axial pseudovectors
C  AAT     -  on input contains the unsymmetrized AATs
C  V       -  on exit contains symmetrized AATs
C  Idon    -  scratch (counting atoms already done)      
C
C
      DIMENSION NEqATM(NAtoms,NTrans),TRANS(3,3,NTrans),ISYM(Natoms),
     $          P(3,3),AAT(3,NAtoms,3),V(3,3,NAtoms),Idon(Natoms),
     $          TRANSA(3,3,NTrans)
C
      call zeroit(Idon,Natoms)
      DO 10 IAtm=1,NAtoms  ! over unique atoms
      If(ISYM(IAtm).EQ.0) GO TO 10
c
      DO 20 IOP=1,NTrans
      KAtm = NEqATM(IAtm,IOP)
      If(Idon(KAtm).EQ.1) GO TO 20
      Idon(KAtm)=1     
C
C  Form TRANSA * DipD * TRANS(t)
C
      CALL ZeroIT(P,9)
      DO 50 I=1,3
      DO 49 J=1,3
      DO 48 K=1,3
      P(I,J) = P(I,J) + TRANSA(I,K,IOP)*AAT(K,IAtm,J)
 48   CONTINUE
 49   CONTINUE
 50   CONTINUE
c
      DO 60 I=1,3
      DO 59 J=1,3
      DO 58 K=1,3
      V(J,I,KAtm) = V(J,I,KAtm) + P(I,K)*TRANS(J,K,IOP)
 58   CONTINUE
 59   CONTINUE
 60   CONTINUE
c
 20   CONTINUE
 10   CONTINUE
C
      RETURN
      END
c ========================================================================
c
      SUBROUTINE AssignSYM(NTrans, TRANS, TRANSA)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ARGUMENTS
C
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 matrices
C  TRANSA  -  on exit symmetry operations for axial pseudovectors
C
C
      DIMENSION TRANS(3,3,NTrans),TRANSA(3,3,NTrans)
c
      Parameter (one=1.0d0)
C
cccccc
cc      DO 10 I=1,NTrans
cc        WRITE(6,1100) I
cc 1100   FORMAT(/,' Symmetry operation: ',I3)
cc        CALL PrntMat(3,3,3,TRANS(1,1,I))
cc 10   CONTINUE
cccccc
c
      CALL ZeroIT(TRANSA,9*NTrans)
C
C  Loop over symmetry operations
C
      DO I=1,NTrans
c -- sort out the abelian operations
      If(trans(1,1,I).eq.one.and.trans(2,2,I).eq.one.and.
     $        trans(3,3,I).eq.one) Then
c -- identity
        transA(1,1,I) = one
        transA(2,2,I) = one
        transA(3,3,I) = one
      Else If(trans(1,1,I).eq.-one.and.trans(2,2,I).eq.one.and.
     $        trans(3,3,I).eq.one) Then
c -- reflection in YZ plane
        transA(1,1,I) = one
        transA(2,2,I) = -one
        transA(3,3,I) = -one
      Else If(trans(1,1,I).eq.one.and.trans(2,2,I).eq.-one.and.
     $        trans(3,3,I).eq.one) Then
c -- reflection in XZ plane
        transA(1,1,I) = -one
        transA(2,2,I) = one
        transA(3,3,I) = -one
      Else If(trans(1,1,I).eq.one.and.trans(2,2,I).eq.one.and.
     $        trans(3,3,I).eq.-one) Then
c -- reflection in XY plane
        transA(1,1,I) = -one
        transA(2,2,I) = -one
        transA(3,3,I) = one
      Else If(trans(1,1,I).eq.-one.and.trans(2,2,I).eq.-one.and.
     $        trans(3,3,I).eq.one) Then
c -- rotation about Z axis
        transA(1,1,I) = -one
        transA(2,2,I) = -one
        transA(3,3,I) = one
      Else If(trans(1,1,I).eq.-one.and.trans(2,2,I).eq.one.and.
     $        trans(3,3,I).eq.-one) Then
c -- rotation about Y axis
        transA(1,1,I) = -one
        transA(2,2,I) = one
        transA(3,3,I) = -one
      Else If(trans(1,1,I).eq.one.and.trans(2,2,I).eq.-one.and.
     $        trans(3,3,I).eq.-one) Then
c -- rotation about X axis
        transA(1,1,I) = one
        transA(2,2,I) = -one
        transA(3,3,I) = -one
      Else If(trans(1,1,I).eq.-one.and.trans(2,2,I).eq.-one.and.
     $        trans(3,3,I).eq.-one) Then
c -- inversion centre
        transA(1,1,I) = one
        transA(2,2,I) = one
        transA(3,3,I) = one
      EndIf
c
      EndDo
cccccc
cc      write(6,*) ' In AssignSYM - Pseudovector Symmetry Operations are:'
cc      DO 20 I=1,NTrans
cc        WRITE(6,1100) I
cc        CALL PrntMat(3,3,3,TRANSA(1,1,I))
cc 20   CONTINUE
cccccc
C
      RETURN
      END
