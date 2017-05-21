      Subroutine multipoles_main(
     &    ncs,           inx,           basdat,        ncspl,
     &    cssharps,      na,            ncf,           dens,
     &    sharpovs,      maxovs,        nsharpovs,     dist2mp,
     &    sharpness,     xnuc,          multipolecut,  scal,
     &    Fockmx)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*),xnuc(5,na)
      integer cssharps(ncspl)
      integer sharpness(ncs)
      real*8 dens(ncf,ncf)
      real*8 Fockmx(ncf*(ncf+1)/2)
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
      integer*1 multipolecut(na,na)
c     common /big/bl(300)
      data dzero/0.0d0/
c
c
cc      call secund(t1)
      scalfock=4.0d0/scal
c
c     allocate integer array ishellinfo(2,ncs)
c
      call mmark
      call getint(2*ncs,iishellinfo)
c
      call find_ncfp_descartes(ncs,inx,ncfp,bl(iishellinfo))
cc      call secund(t2)
cc      trest = t2-t1
      if(ncfp .gt. ncf) then
c
c      allocate array densp(ncfp,ncfp)
c
        call getmem(ncfp*ncfp,idensp)
c
c  allocate some work arrays
c
        call mmark
        call getmem(ncfp,iw1)
        call getmem(ncfp,iw2)
        call getmem(ncfp,iw3)
        call getmem(ncfp,iw4)
        call getmem(ncfp,iw5)
        call getmem(ncfp,iw6)
        call getmem(ncfp,iw7)
c
        call denstrafo_to_descartes(
     &    ncs,        inx,        ncf,        ncfp,       dens,
     &    bl(idensp), bl(iw1),    bl(iw2),    bl(iw3),    bl(iw4),
     &    bl(iw5),    bl(iw6),    bl(iw7))
        idescartestrafo=1
        call retmark ! deallocate work arrays
c
      else
        idescartestrafo=0
      end if
cc      call secund(t3)
cc      tdes = t3-t2
c
c    allocate array atomicpoles(35,na)
c
      call getmem(35*na,iatomicpoles)
      call zeroit(bl(iatomicpoles),35*na)
c
      if(idescartestrafo .eq. 0) then
        call collect_multipoles(
     &    ncs,           inx,           basdat,        ncspl,
     &    cssharps,      na,            ncf,           dens,
     &    sharpovs,      maxovs,        nsharpovs,     bl(iatomicpoles),
     &    sharpness,     bl(iishellinfo))
      else if(idescartestrafo .eq. 1) then
        call collect_multipoles(
     &    ncs,           inx,           basdat,        ncspl,
     &    cssharps,      na,            ncfp,          bl(idensp),
     &    sharpovs,      maxovs,        nsharpovs,     bl(iatomicpoles),
     &    sharpness,     bl(iishellinfo))
      end if
cc      call secund(t4)
cc      tmult = t4-t3
c
c     allocate array v(8,na)
c
      call getmem(8*na,iv)
      call zeroit(bl(iv),8*na)
c
c
cc      call secund(t1)
      DO ish=1,ncspl
        ics=cssharps(ish)
        jiatom=inx(2,ics)
        idim=inx(3,ics)
        if(idim .eq. 5) then
          idimp=6
        else if(idim .eq. 7) then
          idimp=10
        else
          idimp=idim
        end if
        ilastf=inx(10,ics)
        ifirstf=inx(11,ics)+1
        itype=inx(12,ics)
        Rxj=xnuc(2,jiatom)
        Ryj=xnuc(3,jiatom)
        Rzj=xnuc(4,jiatom)
        do iatom=1,na
          if(multipolecut(iatom,jiatom) .eq. 0) then
            Rxi=xnuc(2,iatom)              !they are far enough
            Ryi=xnuc(3,iatom)
            Rzi=xnuc(4,iatom)
            Rx=Rxj-Rxi
            Ry=Ryj-Ryi
            Rz=Rzj-Rzi
            call multipole_intermediates(bl(iatomicpoles+35*(iatom-1)),
     &                                     Rx,Ry,Rz,bl(iv+8*(iatom-1)))
          end if
        end do
        do jsh=1,nsharpovs(ish)
          jcs=sharpovs(jsh,ish)
          if(sharpness(jcs) .eq. 3) goto 150 !the function is sharp
          jdim=inx(3,jcs)
          if(jdim .eq. 5) then
            jdimp=6
          else if(jdim .eq. 7) then
            jdimp=10
          else
            jdimp=jdim
          end if
          jlastf=inx(10,jcs)
          jfirstf=inx(11,jcs)+1
          jtype=inx(12,jcs)
c
c    allocate arrays polesi(jdimp,idimp,35),polesit(35,jdimp,idimp)
c                    fcollect(jdimp,idimp),fcollect2(jdim,idim),
c                    poleshelp(jdimp,idimp,35)
c
          lenmp=idimp*jdimp
          lenmp35=35*lenmp
          call mmark
          call getmem(lenmp35,ipolesi)
          call zeroit(bl(ipolesi),lenmp35)
          call getmem(lenmp35,ipolesit)
          call getmem(lenmp,ifcollect)
          call zeroit(bl(ifcollect),lenmp)
          call getmem(jdim*idim,ifcollect2)
          call getmem(lenmp35,ipoleshelp)
c
          call mutipoles_jcs_ics(
     &    jcs,           ics,           jdimp,         idimp,
     &    jtype,         itype,         ncs,           inx,
     &    basdat,        bl(ipoleshelp),bl(ipolesi))
c
c   now transpose the array polesi. Do not ask me why
c
          call transpole(idimp,jdimp,bl(ipolesi),bl(ipolesit))
cc          do icomp=1,idimp
cc            i2=jdimp*(icomp-1)
cc            i3t=35*i2
cc            do ipole=1,35
cc              i3=lenmp*(ipole-1)
cc              i1t=ipole-1
cc              do jcomp=1,jdimp
cc                i1=jcomp-1
cc                i2t=35*(jcomp-1)
cc                bl(ipolesit+i1t+i2t+i3t)=bl(ipolesi+i1+i2+i3)
cc              end do
cc            end do
cc          end do
c
c
          do jatom=1,na
            if(multipolecut(jatom,jiatom) .eq. 0) then
              Rxi=xnuc(2,jatom)               !they are far enough
              Ryi=xnuc(3,jatom)
              Rzi=xnuc(4,jatom)
              Rx=Rxj-Rxi
              Ry=Ryj-Ryi
              Rz=Rzj-Rzi
              R=dsqrt(Rx**2+Ry**2+Rz**2)
c
c Now the main things are coming
c
c
              do icomp=1,idimp
                ifcoll=jdimp*(icomp-1)
                i3=35*ifcoll
                do jcomp=1,jdimp
                  i2=35*(jcomp-1)
                  call multipoles_interaction_R5(
     &        bl(ipolesit+i2+i3),bl(iatomicpoles+35*(jatom-1)),
     &        bl(iv+8*(jatom-1)),Rx,
     &        Ry,Rz,R,final)
c
                  bl(ifcollect+ifcoll+jcomp-1)=
     &               bl(ifcollect+ifcoll+jcomp-1)+final
                end do
              end do
c
            end if
          end do !for jatom
c
c Now transform it back if necessary
c
          if((idim .ne. idimp) .or. (jdim .ne. jdimp)) then
c
c  allocate array work(jdim,idimp)
c
            call getmem(jdim*idimp,iwork)
c
            call shellpair_trafo_from_descartes(idim,jdim,idimp,
     &           jdimp,bl(ifcollect2),bl(ifcollect),bl(iwork))
            call retmem(1) ! deallocate work
c
c
            if(ics .gt. jcs)then
              icomp=0
              do ifunc=ifirstf,ilastf
                indxi=ifunc*(ifunc-1)/2
                icomp=icomp+1
                jcomp=0
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect2+jdim*(icomp-1)+jcomp-1)
                end do
              end do
            else
              icomp=0
              do ifunc=ifirstf,ilastf
                icomp=icomp+1
                jcomp=0
                indxi=ifunc
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc*(jfunc-1)/2
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect2+jdim*(icomp-1)+jcomp-1)
                end do
              end do
            end if
c
          else
c
            if(ics .gt. jcs)then
              icomp=0
              do ifunc=ifirstf,ilastf
                indxi=ifunc*(ifunc-1)/2
                icomp=icomp+1
                jcomp=0
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect+jdimp*(icomp-1)+jcomp-1)
                end do
              end do
            else
              icomp=0
              do ifunc=ifirstf,ilastf
                indxi=ifunc
                icomp=icomp+1
                jcomp=0
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc*(jfunc-1)/2
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect+jdimp*(icomp-1)+jcomp-1)
                end do
              end do
            end if
c
          end if
c
          call retmark  ! deallocate polesi,polesit,fcollect,fcollect2,
                        ! poleshelp
 150      continue
        end do
      END DO
cc      call secund(t2)
cc      tloop = t2-t1
c
      call retmark  ! deallocate ishellinfo,densp,atomicpoles,v
c
cc      write(6,*) ' ** MULTIPOLE TIMINGS **'
cc      write(6,*) ' trest:',trest
cc      write(6,*) ' tdes:',tdes,' ncf:',ncf,' ncfp:',ncfp
cc      write(6,*) ' tmult:',tmult
cc      write(6,*) ' tloop:',tloop
      return
      end
c
c***********************************************************************
c
      Subroutine get_multipolecutmx(na,dist2mp,xnuc,multipolecut)
      IMPLICIT REAL*8(A-H,O-Z)
      integer*1 multipolecut(na,na)
      real*8 xnuc(5,na)
c
      do iatom=1,na
        pxi=xnuc(2,iatom)
        pyi=xnuc(3,iatom)
        pzi=xnuc(4,iatom)
        do jatom=1,na
          pxj=xnuc(2,jatom)
          pyj=xnuc(3,jatom)
          pzj=xnuc(4,jatom)
          dist=dsqrt((pxj-pxi)**2+(pyj-pyi)**2+(pzj-pzi)**2)
          if(dist .gt. dist2mp) then
            multipolecut(jatom,iatom)=0
          else
            multipolecut(jatom,iatom)=1
          end if
        end do
      end do

c
      return
      end
c
c***********************************************************************
c
      Subroutine shellpair_trafo_from_descartes(idim,jdim,idimp,
     &           jdimp,dens,densp,dens1)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 dens(jdim,idim)
      real*8 densp(jdimp,idimp)
      real*8 dens1(jdim,idimp)
c
c
      dsqrt3=dsqrt(3.0d0)
      c1=dsqrt3/6.0d0
      c2=-c1
      c3=2.0d0*c1
c
      cc1=-1.0d0
      cc2=4.0d0
      cc3=dsqrt(40.0d0)
c
      if(idim .eq. idimp) then
c
        if(jdimp .eq. 6) then
          do icount=1,idimp
            d1=densp(1,icount)
            d2=densp(2,icount)
            d3=densp(3,icount)
            d4=densp(4,icount)
            d5=densp(5,icount)
            d6=densp(6,icount)
c
            dens(1,icount)=c2*(d1+d2)+c3*d3
            dens(2,icount)=0.5d0*(d1-d2)
            dens(3,icount)=d4
            dens(4,icount)=d5
            dens(5,icount)=d6
          end do
c
        else if(jdimp .eq. 10) then
          do icount=1,idimp
            d1=densp(1,icount)
            d2=densp(2,icount)
            d3=densp(3,icount)
            d4=densp(4,icount)
            d5=densp(5,icount)
            d6=densp(6,icount)
            d7=densp(7,icount)
            d8=densp(8,icount)
            d9=densp(9,icount)
            d10=densp(10,icount)
c
            dens(1,icount)=cc2*d2+cc1*(d7+d9)
            dens(2,icount)=cc2*d3+cc1*(d8+d10)
            dens(3,icount)=cc2*d4+cc1*(d1+d6)
            dens(4,icount)=cc2*d8+cc1*(d3+d10)
            dens(5,icount)=cc2*d6+cc1*(d1+d4)
            dens(6,icount)=cc2*d9+cc1*(d2+d7)
            dens(7,icount)=cc3*d5
          end do
        end if
c
c
      else if(jdim .eq. jdimp) then
c
        if(idim .eq. 5) then
c
          do jcount=1,jdim
            dens(jcount,1)=c2*(densp(jcount,1)+densp(jcount,2))+
     &                 c3*densp(jcount,3)
            dens(jcount,2)=0.5d0*(densp(jcount,1)-densp(jcount,2))
            dens(jcount,3)=densp(jcount,4)
            dens(jcount,4)=densp(jcount,5)
            dens(jcount,5)=densp(jcount,6)
          end do
c
        else if(idim .eq. 7) then
          do jcount=1,jdim
c
            dens(jcount,1)=cc2*densp(jcount,2)+
     &      cc1*(densp(jcount,7)+densp(jcount,9))
            dens(jcount,2)=cc2*densp(jcount,3)+
     &      cc1*(densp(jcount,8)+densp(jcount,10))
            dens(jcount,3)=cc2*densp(jcount,4)+
     &      cc1*(densp(jcount,1)+densp(jcount,6))
            dens(jcount,4)=cc2*densp(jcount,8)+
     &      cc1*(densp(jcount,3)+densp(jcount,10))
            dens(jcount,5)=cc2*densp(jcount,6)+
     &      cc1*(densp(jcount,1)+densp(jcount,4))
            dens(jcount,6)=cc2*densp(jcount,9)+
     &      cc1*(densp(jcount,2)+densp(jcount,7))
            dens(jcount,7)=cc3*densp(jcount,5)
c
          end do
        end if
c
      else
        if(jdimp .eq. 6) then
          do icount=1,idimp
            d1=densp(1,icount)
            d2=densp(2,icount)
            d3=densp(3,icount)
            d4=densp(4,icount)
            d5=densp(5,icount)
            d6=densp(6,icount)
c
            dens1(1,icount)=c2*(d1+d2)+c3*d3
            dens1(2,icount)=0.5d0*(d1-d2)
            dens1(3,icount)=d4
            dens1(4,icount)=d5
            dens1(5,icount)=d6
          end do
c
        else if(jdimp .eq. 7) then
          do icount=1,idimp
            d1=densp(1,icount)
            d2=densp(2,icount)
            d3=densp(3,icount)
            d4=densp(4,icount)
            d5=densp(5,icount)
            d6=densp(6,icount)
            d7=densp(7,icount)
            d8=densp(8,icount)
            d9=densp(9,icount)
            d10=densp(10,icount)
c
            dens1(1,icount)=cc2*d2+cc1*(d7+d9)
            dens1(2,icount)=cc2*d3+cc1*(d8+d10)
            dens1(3,icount)=cc2*d4+cc1*(d1+d6)
            dens1(4,icount)=cc2*d8+cc1*(d3+d10)
            dens1(5,icount)=cc2*d6+cc1*(d1+d4)
            dens1(6,icount)=cc2*d9+cc1*(d2+d7)
            dens1(7,icount)=cc3*d5
          end do
        end if
c
c
c Now the j is already transformed. Trafo for i is coming.
c
        if(idim .eq. 5) then
c
          do jcount=1,jdim
            dens(jcount,1)=c2*(dens1(jcount,1)+dens1(jcount,2))+
     &                 c3*dens1(jcount,3)
            dens(jcount,2)=0.5d0*(dens1(jcount,1)-dens1(jcount,2))
            dens(jcount,3)=dens1(jcount,4)
            dens(jcount,4)=dens1(jcount,5)
            dens(jcount,5)=dens1(jcount,6)
          end do
c
        else if(idim .eq. 7) then
          do jcount=1,jdim
c
            dens(jcount,1)=cc2*dens1(jcount,2)+
     &      cc1*(dens1(jcount,7)+dens1(jcount,9))
            dens(jcount,2)=cc2*dens1(jcount,3)+
     &      cc1*(dens1(jcount,8)+dens1(jcount,10))
            dens(jcount,3)=cc2*dens1(jcount,4)+
     &      cc1*(dens1(jcount,1)+dens1(jcount,6))
            dens(jcount,4)=cc2*dens1(jcount,8)+
     &      cc1*(dens1(jcount,3)+dens1(jcount,10))
            dens(jcount,5)=cc2*dens1(jcount,6)+
     &      cc1*(dens1(jcount,1)+dens1(jcount,4))
            dens(jcount,6)=cc2*dens1(jcount,9)+
     &      cc1*(dens1(jcount,2)+dens1(jcount,7))
            dens(jcount,7)=cc3*dens1(jcount,5)
c
          end do
        end if
c
      end if
c
      return
      end
c
c
c************************************************************************
c
      Subroutine denstrafo_to_descartes(
     &    ncs,        inx,        ncf,        ncfp,      dens,
     &    densp,      v1,         v2,         v3,        v4,
     &    v5,         v6,         v7)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      real*8 dens(ncf,ncf)
      real*8 densp(ncfp,ncfp)
      real*8 v1(*)
      real*8 v2(*)
      real*8 v3(*)
      real*8 v4(*)
      real*8 v5(*)
      real*8 v6(*)
      real*8 v7(*)
c
      dsqrt3=dsqrt(3.0d0)
      c1=dsqrt3/6.0d0
      c2=-c1
      c3=2.0d0*c1
c
      cc1=-1.0d0
      cc2=4.0d0
      cc3=dsqrt(40.0d0)
c
      indexi=0
      indexip=0
      do ics=1,ncs
        idim=inx(3,ics)
        if(idim .eq. 5) then
          ijump=1
        else if(idim .eq. 7) then
          ijump=0
        else
          ijump=0
        end if
        indexip=indexip+ijump !starting from 2 for d and from 1 for f
        do icont=1,idim
          indexi=indexi+1
          indexip=indexip+1
          indexj=0
          indexjp=0
          do jcs=1,ncs
c
            jdim=inx(3,jcs)
            if(jdim .eq. 5) then
              indexj=indexj+1
              d1=dens(indexj,indexi)
              indexj=indexj+1
              d2=dens(indexj,indexi)
              indexj=indexj+1
              d3=dens(indexj,indexi)
              indexj=indexj+1
              d4=dens(indexj,indexi)
              indexj=indexj+1
              d5=dens(indexj,indexi)
c
              indexjp=indexjp+1
              densp(indexjp,indexip)=c2*d1+0.5d0*d2
              indexjp=indexjp+1
              densp(indexjp,indexip)=c2*d1-0.5d0*d2
              indexjp=indexjp+1
              densp(indexjp,indexip)=c3*d1
              indexjp=indexjp+1
              densp(indexjp,indexip)=d3
              indexjp=indexjp+1
              densp(indexjp,indexip)=d4
              indexjp=indexjp+1
              densp(indexjp,indexip)=d5
c
            else if(jdim .eq. 7) then
              indexj=indexj+1
              d1=dens(indexj,indexi)
              indexj=indexj+1
              d2=dens(indexj,indexi)
              indexj=indexj+1
              d3=dens(indexj,indexi)
              indexj=indexj+1
              d4=dens(indexj,indexi)
              indexj=indexj+1
              d5=dens(indexj,indexi)
              indexj=indexj+1
              d6=dens(indexj,indexi)
              indexj=indexj+1
              d7=dens(indexj,indexi)
c
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc1*(d3+d5)
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc2*d1+cc1*d6
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc2*d2+cc1*d4
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc2*d3+cc1*d5
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc3*d7
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc1*d3+cc2*d5
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc1*(d1+d6)
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc1*d2+cc2*d4
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc1*d1+cc2*d6
              indexjp=indexjp+1
              densp(indexjp,indexip)=cc1*(d2+d4)
c
            else
              do jcont=1,jdim
                indexjp=indexjp+1
                indexj=indexj+1
                densp(indexjp,indexip)=dens(indexj,indexi)
              end do
            end if
c
          end do
        end do
      end do
c
c Now the j is already transformed. Trafo for i is coming.
c
      indexip=0
      do ics=1,ncs
        idim=inx(3,ics)
        if(idim .eq. 5) then
c
          indexip=indexip+2
          do jcf=1,ncfp
            v1(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v2(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip-3
c
          indexip=indexip+1 !first
          do jcf=1,ncfp
            densp(jcf,indexip)=c2*v1(jcf)+0.5d0*v2(jcf)
          end do
          indexip=indexip+1 !second
          do jcf=1,ncfp
            densp(jcf,indexip)=c2*v1(jcf)-0.5d0*v2(jcf)
          end do
          indexip=indexip+1 !third
          do jcf=1,ncfp
            densp(jcf,indexip)=c3*v1(jcf)
          end do
          indexip=indexip+3 !All 6 is done
c
        else if(idim .eq. 7) then
c
          indexip=indexip+1
          do jcf=1,ncfp
            v1(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v2(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v3(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v4(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v5(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v6(jcf)=densp(jcf,indexip)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            v7(jcf)=densp(jcf,indexip)
          end do
c
          indexip=indexip-6
c
          do jcf=1,ncfp
            densp(jcf,indexip)=cc1*(v3(jcf)+v5(jcf))
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc2*v1(jcf)+cc1*v6(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc2*v2(jcf)+cc1*v4(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc2*v3(jcf)+cc1*v5(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc3*v7(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc1*v3(jcf)+cc2*v5(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc1*(v1(jcf)+v6(jcf))
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc1*v2(jcf)+cc2*v4(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc1*v1(jcf)+cc2*v6(jcf)
          end do
          indexip=indexip+1
          do jcf=1,ncfp
            densp(jcf,indexip)=cc1*(v2(jcf)+v4(jcf))
          end do
c
        else
          indexip=indexip+idim
        end if
      end do
c
      return
      end
c
c***********************************************************************
c
      Subroutine find_ncfp_descartes(ncs,inx,ncfp,ishellinfo)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      integer ishellinfo(2,ncs)
c
      ncfp=0
      icsprev=1
      do ics=1,ncs
        ishellinfo(1,ics)=icsprev
        idim=inx(3,ics)
        if(idim .eq. 5) then
          idimp=6
        else if(idim .eq. 7) then
          idimp=10
        else
          idimp=idim
        end if
        ncfp=ncfp+idimp
        icsprev=icsprev+idimp
        ishellinfo(2,ics)=icsprev-1
      end do
c
      return
      end
c
c*************************************************************************
c
      Subroutine multipole_intermediates(poles,Rx,Ry,Rz,v)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 poles(35) !precalculated multi-poles
      real*8 v(8) !output
c
      t5 = Rx*Rz
      t8 = Ry*Rz
      t11 = Rx*Ry
      t14 = Rx**2
      t16 = Rz**2
      t18 = Ry**2
      t37 = t14*Rx
      t39 = t14*Ry
      t45 = Rx*t18
      t48 = t18*Ry
      t53 = t16*Rz
      t86 = t18*poles(31)+2.D0*t8*poles(34)+t14*poles(24)+2.D0*t5*poles(
     #23)+t16*poles(33)+2.D0*t8*poles(25)+t14*poles(26)+t18*poles(33)+2.
     #D0*t11*poles(29)+t16*poles(26)+2.D0*t11*poles(27)+2.D0*t8*poles(32
     #)+2.D0*t5*poles(28)+2.D0*t5*poles(30)+t18*poles(24)+2.D0*t11*poles
     #(22)+t16*poles(35)+t14*poles(21)
      t99 = t16**2
      t119 = t14**2
      t121 = t18**2
      t129 = 4.D0*Rx*t53*poles(30)+4.D0*Ry*t53*poles(34)+4.D0*Rx*t48*pol
     #es(27)+12.D0*t11*t16*poles(29)+t99*poles(35)+4.D0*t37*Rz*poles(23)
     #+6.D0*t18*t16*poles(33)+12.D0*t39*Rz*poles(25)+12.D0*t45*Rz*poles(
     #28)+4.D0*t37*Ry*poles(22)+4.D0*t48*Rz*poles(32)+t119*poles(21)+t12
     #1*poles(31)+6.D0*t14*t18*poles(24)+6.D0*t14*t16*poles(26)
      v(1) = Rx*poles(2)+Ry*poles(3)+Rz*poles(4)
      v(2) = 2.D0*t5*poles(7)+2.D0*t8*poles(9)+2.D0*t11*poles(6)+t14*pol
     #es(5)+t16*poles(10)+t18*poles(8)
      v(3) = Rz*poles(18)+Ry*poles(19)+Rx*poles(14)+Rz*poles(13)+Rx*pole
     #s(16)+Ry*poles(12)+Rx*poles(11)+Rz*poles(20)+Ry*poles(17)
      v(4) = 3.D0*t18*Rz*poles(18)+3.D0*Ry*t16*poles(19)+t37*poles(11)+3
     #.D0*t39*poles(12)+3.D0*Rx*t16*poles(16)+3.D0*t45*poles(14)+t48*pol
     #es(17)+6.D0*t11*Rz*poles(15)+t53*poles(20)+3.D0*t14*Rz*poles(13)
      v(5) = t86
      v(6) = t129
      v(7) = poles(5)+poles(8)+poles(10)
      v(8) = poles(21)+2.D0*poles(24)+2.D0*poles(26)+poles(31)+2.D0*pole
     #s(33)+poles(35)

c
      return
      end
c
c*************************************************************************
c
      Subroutine collect_multipoles(
     &    ncs,           inx,           basdat,        ncspl,
     &    cssharps,      na,            ncf,           dens,
     &    sharpovs,      maxovs,        nsharpovs,     atomicpoles,
     &    sharpness,     ishellinfo)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer cssharps(ncspl)
      integer sharpness(ncs)
      real*8 dens(ncf,ncf)
      real*8 atomicpoles(35,na)!atomic poles times density mx
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
      integer ishellinfo(2,ncs)
c     common /big/bl(5000000)
      data dzero/0.0d0/
c
      do ish=1,ncspl
      ics=cssharps(ish)
      iatom=inx(2,ics)
      idim=inx(3,ics)
      itype=inx(12,ics)
      if(itype .eq. 4) then !d->d6
      idim=6
      else if(itype .eq. 6) then !f->f10
      idim=10
      end if
c      ilastf=inx(10,ics)
c      ifirstf=inx(11,ics)+1
      ifirstf=ishellinfo(1,ics)
      ilastf=ishellinfo(2,ics)
        do jsh=1,nsharpovs(ish)
        jcs=sharpovs(jsh,ish)
        if(sharpness(jcs) .eq. 3) goto 150 !the function is sharp
        jdim=inx(3,jcs)
        jtype=inx(12,jcs)
        if(jtype .eq. 4) then !d->d6
        jdim=6
        else if(jtype .eq. 6) then !f->f10
        jdim=10
        end if
c        jlastf=inx(10,jcs)
c        jfirstf=inx(11,jcs)+1
        jfirstf=ishellinfo(1,jcs)
        jlastf=ishellinfo(2,jcs)
c
c      allocate arrays poles(jdim,idim,35),poleshelp(jdim,idim,35)
c
        call mmark
        call getmem(35*idim*jdim,ipoles)
        call zeroit(bl(ipoles),35*idim*jdim)
        call getmem(35*idim*jdim,ipoleshelp)
c
        call mutipoles_jcs_ics(jcs,ics,jdim,idim,jtype,
     &                         itype,ncs,inx,basdat,
     &                         bl(ipoleshelp),bl(ipoles))
c
          do ipole=1,35
            i3=idim*jdim*(ipole-1)
            icomp=0
            do ifunc=ifirstf,ilastf
              icomp=icomp+1
              i2=jdim*(icomp-1)
              jcomp=0
              do jfunc=jfirstf,jlastf
              jcomp=jcomp+1
              atomicpoles(ipole,iatom)=atomicpoles(ipole,iatom)+
     &              dens(jfunc,ifunc)*bl(ipoles+jcomp-1+i2+i3)
              end do
            end do
          end do
c
        call retmark !  deallocate poles,poleshelp
 150    continue
        end do
      end do
c
      return
      end
c
c***********************************************************************
c
      Subroutine mutipoles_jcs_ics(jcs,ics,jdim,idim,jtype,
     &                             itype,ncs,inx,basdat,
     &                             poleshelp,poles)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 poles(jdim,idim,35) ! up to lmax=4
      real*8 poleshelp(jdim,idim,35) ! up to lmax=4
      data pi/3.1415926535897932384626433d0/
      data izero/0/
      dimension xa(3),xb(3)
c     common /big/bl(5000000)
c
c    allocate array xint(jdim,idim)
c
      call getmem(jdim*idim,ixint)
c
      const=2.0d0*dsqrt(2.0d0)*dsqrt(pi)**(-3)
ct      const=1.0d0
c
c
      ibeg=inx(1,ics)+1
      xa(1)=basdat(11,ibeg)
      xa(2)=basdat(12,ibeg)
      xa(3)=basdat(13,ibeg)
      iend=inx(5,ics)
c
      jbeg=inx(1,jcs)+1
      xb(1)=basdat(11,jbeg)
      xb(2)=basdat(12,jbeg)
      xb(3)=basdat(13,jbeg)
      jend=inx(5,jcs)
c
      dist2=(xa(1)-xb(1))**2+(xa(2)-xb(2))**2+(xa(3)-xb(3))**2
c
c
      do icont=ibeg,iend  ! loop over the contractions
        etaa=basdat(1,icont)
        cca1=basdat(2,icont)
        cca2=basdat(3,icont) !for itype=3 and icomp=1
c
        do jcont=jbeg,jend  ! loop over the contractions
          etab=basdat(1,jcont)
          ccb1=basdat(2,jcont)
          ccb2=basdat(3,jcont) !for jtype=3 and jcomp=1
c
          etasumm=etaa+etab
          etaprod=etaa*etab
          argexp=-etaprod*dist2/etasumm
          s0=const*dexp(argexp)*
     &    dsqrt(dsqrt(etaprod))**3/dsqrt(etasumm)**3
c
          call mutipoles_jsh_ish(jdim,idim,jtype,
     &                             itype,etaa,etab,s0,xa,
     &                             xb,bl(ixint),poleshelp)

c
c
c THIS PART GIVES YOU THE FINAL POLES USING CONTRACTION COEFS AND THE NORMALIZATION
C FACTOR.
c
          do ipole=1,35
            do icomp=1,idim
              if((itype .eq. 3) .and. (icomp .gt. 1)) then
                concoefa=cca2
              else
                concoefa=cca1
              end if
              do jcomp=1,jdim
                if((jtype .eq. 3) .and. (jcomp .gt. 1)) then
                  concoefb=ccb2
                else
                  concoefb=ccb1
                end if
c
ct            poles(jcomp,icomp,ipole)=poles(jcomp,icomp,ipole)+
ct     &        poleshelp(jcomp,icomp,ipole)*const*concoefa*concoefb
c
                poles(jcomp,icomp,ipole)=poles(jcomp,icomp,ipole)+
     &          poleshelp(jcomp,icomp,ipole)*concoefa*concoefb
c
              end do !jcomp
            end do !icomp
          end do !ipole
c
        end do !jcont
      end do !icont
c
      call retmem(1)  !  deallocate xint
c
      return
      end
c
c***********************************************************************
c
      Subroutine mutipoles_jsh_ish(jdim,idim,jtypeo,
     &                             itypeo,etaa,etab,s0,xa,
     &                             xb,xint,poleshelp)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 poleshelp(jdim,idim,35) ! up to lmax=4
      dimension xa(3),xb(3),xint(jdim*idim),la(3),lb(3),l3(3)
      data pi/3.1415926535897932384626433d0/
      data izero/0/
c
      if(itypeo .eq. 5) then
        itype=4
      else if((itypeo .eq. 6) .or. (itypeo .eq. 7)) then
        itype=5
      else
        itype=itypeo
      end if
      if(jtypeo .eq. 5) then
        jtype=4
      else if((jtypeo .eq. 6) .or. (jtypeo .eq. 7)) then
        jtype=5
      else
        jtype=jtypeo
      end if
c
      lb(1)=izero
      lb(2)=izero
      lb(3)=izero
      l3(1)=izero
      l3(2)=izero
      l3(3)=izero
c
c monopole
      la(1)=0
      la(2)=0
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,1)=xint(index)
        end do
      end do
c
c x
      la(1)=1
      la(2)=0
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,2)=xint(index)
        end do
      end do
c
c y
      la(1)=0
      la(2)=1
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,3)=xint(index)
        end do
      end do
c
c z
      la(1)=0
      la(2)=0
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,4)=xint(index)
        end do
      end do
c
c xx
      la(1)=2
      la(2)=0
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,5)=xint(index)
        end do
      end do
c
c xy
      la(1)=1
      la(2)=1
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,6)=xint(index)
        end do
      end do
c
c xz
      la(1)=1
      la(2)=0
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,7)=xint(index)
        end do
      end do
c
c y*y
      la(1)=0
      la(2)=2
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,8)=xint(index)
        end do
      end do
c
c y*z
      la(1)=0
      la(2)=1
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,9)=xint(index)
        end do
      end do
c
c z*z
      la(1)=0
      la(2)=0
      la(3)=2
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,10)=xint(index)
        end do
      end do
c
c x*x*x
      la(1)=3
      la(2)=0
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,11)=xint(index)
        end do
      end do
c
c x*x*y
      la(1)=2
      la(2)=1
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,12)=xint(index)
        end do
      end do
c
c x*x*z
      la(1)=2
      la(2)=0
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,13)=xint(index)
        end do
      end do
c
c x*y*y
      la(1)=1
      la(2)=2
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,14)=xint(index)
        end do
      end do
c
c x*y*z
      la(1)=1
      la(2)=1
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,15)=xint(index)
        end do
      end do
c
c x*z*z
      la(1)=1
      la(2)=0
      la(3)=2
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,16)=xint(index)
        end do
      end do
c
c y*y*y
      la(1)=0
      la(2)=3
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,17)=xint(index)
        end do
      end do
c
c y*y*z
      la(1)=0
      la(2)=2
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,18)=xint(index)
        end do
      end do
c
c y*z*z
      la(1)=0
      la(2)=1
      la(3)=2
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,19)=xint(index)
        end do
      end do
c
c z*z*z
      la(1)=0
      la(2)=0
      la(3)=3
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,20)=xint(index)
        end do
      end do
c
c x^4
      la(1)=4
      la(2)=0
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,21)=xint(index)
        end do
      end do
c
c x^3*y
      la(1)=3
      la(2)=1
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,22)=xint(index)
        end do
      end do
c
c x^3*z
      la(1)=3
      la(2)=0
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,23)=xint(index)
        end do
      end do
c
c x^2*y^2
      la(1)=2
      la(2)=2
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,24)=xint(index)
        end do
      end do
c
c x^2*y*z
      la(1)=2
      la(2)=1
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,25)=xint(index)
        end do
      end do
c
c x^2*z^2
      la(1)=2
      la(2)=0
      la(3)=2
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,26)=xint(index)
        end do
      end do
c
c x*y^3
      la(1)=1
      la(2)=3
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,27)=xint(index)
        end do
      end do
c
c x*y^2*z
      la(1)=1
      la(2)=2
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,28)=xint(index)
        end do
      end do
c
c x*y*z^2
      la(1)=1
      la(2)=1
      la(3)=2
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,29)=xint(index)
        end do
      end do
c
c x*z^3
      la(1)=1
      la(2)=0
      la(3)=3
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,30)=xint(index)
        end do
      end do
c
c y^4
      la(1)=0
      la(2)=4
      la(3)=0
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,31)=xint(index)
        end do
      end do
c
c y^3*z
      la(1)=0
      la(2)=3
      la(3)=1
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,32)=xint(index)
        end do
      end do
c
c y^2*z^2
      la(1)=0
      la(2)=2
      la(3)=2
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,33)=xint(index)
        end do
      end do
c
c y*z^3
      la(1)=0
      la(2)=1
      la(3)=3
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,34)=xint(index)
        end do
      end do
c
c z^4
      la(1)=0
      la(2)=0
      la(3)=4
      call intar_fml(itype,  jtype,  etaa,  etab,  s0,
     1           xa,    xb,    la, lb, l3,
     2           jdim,idim,xint)
      index=0
      do ish=1,idim
        do jsh=1,jdim
          index=index+1
          poleshelp(jsh,ish,35)=xint(index)
        end do
      end do
c
      return
      end
c
c***********************************************************************
c
      subroutine intar_fml(ityp,  jtyp,  a,  b,  s0,
     1                  xa,    xb,    la, lb, l3,
     2                  jdim,idim,xint)
      implicit real*8 (a-h,o-z)
c
c this routine calculates one-electron integrals over gaussian shells
c it is assumed that the 1s (spherical) part of the gaussian is
c normalized (this is denoted by 'gaussian').
C
C FML: YES IT IS NORMALIZED TO pi^(3/2) BUT NOT TO ONE !!! This subroutine
c      is modified by me using dynamical memory allocation for xint instead
c      of using always the biggest possible  xint(784) array.
c
c the integrand is (x-xa(1))**la(1)*...(z-xb(3))**lb(3)*x**l3(1)..
c  * z**l3(3)*(x-xa(1))**mulx1*...(z-zb(3))**mulz2*gaussian1*gaussian2
c  i.e. it can include, besides the x,y,z factors associated with the
c  angular momentum functions, extra (x-xa),...(z-zb) factors, as well
c  as extra x,y,z factors
c
c here xa(1), xa(2), xa(3) are the coordinates of center a
c xb(1),(2),(3) is the same for center b. a and b are the
c exponents of the two gaussians
c the arrays mulx, muly and mulz contain the x,y, z factors for
c whole shells of functions. the shell types available are (at this
c level) s, p, l (this is s and p with common exponent), d, f and g.
c these are cartesian functions and the transformation to spherical
c harmonics takes place later. for instance, the array elements
c nfu(ityp)+1 to nfu(ityp+1) refer to a shell of type itype.
c the types are 1-s,2-p,3-l,4-d,5-d,6-g. e.g. elements 9 to 14 of
c mulx,muly and mulz refer to d functions. mulx(9) is 2, showing that
c the first function in the d shell is x**2. muly and mulz(9) are 0.
c mulx(12) is 1, muly(12) is 1, mulz(12) is 0, i.e. the 4th function
c of a d shell (remember, we begin with 9) is a d(xy).
c gaussian exponents.
c s0 is the overlap integral between normalized
c 1s gaussians and it is explicitly evaluated in subroutine onel
c the vector la contains the powers of (x-xa), (y-ya),(z-za).
c lb contains (x-xb) etc. l3 contains the powers of x,y,z.
c the results is built in ((xint(j,i),j=1,jsh),i=1,ish)
c where ish,jsh are the shell sizes: 1,3,4,6,10,15,21,28 for
c s,p,l,d,f,g,h,i
c it should be ok up to and including i functions, although it would
c be a good idea to test it for high ang. mom. functions befor using it
c it uses hermite-gaussian quadrature. note that this is not very
c economical but transparent.
c
c  Argument list:
c input:
c ityp,jtyp: function types, 1,2,3..8 for s,p,l,d,f,g,h,i
c a,b: Gaussian exponents
c s0: overlap between normalized 1s gaussians
c xa(3),xb(3): orbital centers
c la(3),lb(3),l3(3) : cartesian exponents, see above.
c output:
c  xint(jsh,ish)
C
C
c  arguments
      dimension xa(3),xb(3),xint(jdim,idim),la(3),lb(3),l3(3)
c  data
c      dimension idg(8),nfu(9),mulx(88),muly(88),mulz(88)
      dimension idg(8),nfu(9),mul(3,88)
c  local variables
      dimension xfc(63),yfc(63),zfc(63)
c
      data idg/0,1,1,2,3,4,5,6/
c  nfu stores the beginnings and endings of shells in the arrays
c  Shell type ityp goes from nfu(ityp)+1 to nfu(ityp)
      data nfu/0,1,4,8,14,24,39,60,88/
c  mulx(k,ifu) gives the power of the x,y,z (k=1..3) factors in the shells
      data mul/0,0,0,                                               ! s
     1  1,0,0, 0,1,0, 0,0,1,                                        ! p
     2  0,0,0, 1,0,0, 0,1,0, 0,0,1,                                 ! l
     3  2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1,                   ! d6
     4  3,0,0, 2,1,0, 2,0,1, 1,2,0, 1,1,1, 1,0,2, 0,3,0, 0,2,1,
     5  0,1,2, 0,0,3,                                               ! f10
     6  4,0,0, 3,1,0, 3,0,1, 2,2,0, 2,1,1, 2,0,2, 1,3,0, 1,2,1,
     7  1,1,2, 1,0,3, 0,4,0, 0,3,1, 0,2,2, 0,1,3, 0,0,4,            ! g15
     8  5,0,0, 4,1,0, 4,0,1, 3,2,0, 3,1,1, 3,0,2, 2,3,0, 2,2,1,
     9  2,1,2, 2,0,3, 1,4,0, 1,3,1, 1,2,2, 1,1,3, 1,0,4, 0,5,0,
     &  0,4,1, 0,3,2, 0,2,3, 0,1,4, 0,0,5,                          ! h21
     1  6,0,0, 5,1,0, 5,0,1, 4,2,0, 4,1,1, 4,0,2, 3,3,0, 3,2,1,
     2  3,1,2, 3,0,3, 2,4,0, 2,3,1, 2,2,2, 2,1,3, 2,0,4, 1,5,0,
     3  1,4,1, 1,3,2, 1,2,3, 1,1,4, 1,0,5, 0,6,0, 0,5,1, 0,4,2,
     4  0,3,3, 0,2,4, 0,1,5, 0,0,6/                                 ! i28
      iadg=idg(ityp)
      ibdg=idg(jtyp)
      l1x=iadg+la(1)     ! the highest possible power of (x-Ax)
      l2x=ibdg+lb(1)     ! the highest possible power of (x-Bx)
      l1y=iadg+la(2)     ! the highest possible power of (y-Ay)
      l2y=ibdg+lb(2)     ! the highest possible power of (y-By)
      l1z=iadg+la(3)     ! the highest possible power of (z-Az)
      l2z=ibdg+lb(3)     ! the highest possible power of (z-Bz)
c
      istart=nfu(ityp)+1
      ilen=nfu(ityp+1)-istart+1
      jstart=nfu(jtyp)+1
      jlen=nfu(jtyp+1)-jstart+1
c
      call intar2(ilen, jlen, mul(1,istart), mul(1,jstart), la,
     2            lb,   l3,   a,             b,             xa,
     3            xb,   l1x,  l1y,           l1z,           l2x,
     4            l2y,  l2z,  s0,            xfc,           yfc,
     5            zfc,  xint)
      end
c
c========================================================================
c
c
      Subroutine multipoles_interaction_R5(polesi,polesj,v,Rx,Ry,
     &                                     Rz,R,final)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 polesi(35)
      real*8 polesj(35)
      real*8 v(8)
c
      t7 = Ry*polesi(3)+Rz*polesi(4)+Rx*polesi(2)
      t11 = R**2
      t12 = t11*R
      t13 = 1.D0/t12
      t15 = polesi(5)+polesi(8)+polesi(10)
      t25 = Ry*Rz
      t28 = Rz**2
      t29 = t28*polesi(10)
      t30 = Rx*Ry
      t33 = Ry**2
      t34 = t33*polesi(8)
      t35 = Rx*Rz
      t38 = Rx**2
      t39 = t38*polesi(5)
      t40 = 2.D0*t25*polesi(9)+t29+2.D0*t30*polesi(6)+t34+2.D0*t35*poles
     #i(7)+t39
      t48 = t11**2
      t50 = 1.D0/t48/R
      t61 = Rz*polesi(18)+Ry*polesi(19)+Rx*polesi(14)+Rz*polesi(13)+Rx*p
     #olesi(16)+Ry*polesi(12)+Rx*polesi(11)+Rz*polesi(20)+Ry*polesi(17)
      t114 = t38*Rx
      t116 = t38*Ry
      t122 = Rx*t33
      t125 = t33*Ry
      t130 = t28*Rz
      t135 = 3.D0*t33*Rz*polesi(18)+3.D0*Ry*t28*polesi(19)+t114*polesi(1
     #1)+3.D0*t116*polesi(12)+3.D0*Rx*t28*polesi(16)+3.D0*t122*polesi(14
     #)+t125*polesi(17)+6.D0*t30*Rz*polesi(15)+t130*polesi(20)+3.D0*t38*
     #Rz*polesi(13)
      t146 = 1.D0/t48/t12
      t178 = -0.15D1*polesi(14)*polesj(2)-0.15D1*polesi(12)*polesj(3)+0.
     #375D0*v(8)*polesi(1)+0.3D1*polesi(9)*polesj(9)-0.15D1*polesi(11)*p
     #olesj(2)-0.15D1*polesj(18)*polesi(4)-0.15D1*polesj(19)*polesi(3)-0
     #.15D1*polesj(12)*polesi(3)+0.375D0*polesj(1)*(polesi(21)+2.D0*pole
     #si(24)+2.D0*polesi(26)+polesi(31)+2.D0*polesi(33)+polesi(35))-0.15
     #D1*polesj(13)*polesi(4)-0.15D1*polesj(11)*polesi(2)-0.15D1*polesj(
     #17)*polesi(3)-0.15D1*polesj(14)*polesi(2)
      t207 = 0.15D1*polesi(10)*polesj(10)-0.15D1*polesj(16)*polesi(2)+0.
     #15D1*polesi(5)*polesj(5)-0.15D1*polesi(13)*polesj(4)-0.15D1*polesi
     #(19)*polesj(3)-0.15D1*polesi(20)*polesj(4)-0.15D1*polesi(17)*poles
     #j(3)-0.15D1*polesi(16)*polesj(2)-0.15D1*polesi(18)*polesj(4)+0.75D
     #0*t15*v(7)+0.3D1*polesi(6)*polesj(6)+0.15D1*polesi(8)*polesj(8)+0.
     #3D1*polesi(7)*polesj(7)-0.15D1*polesj(20)*polesi(4)
      t259 = -15.D0*t39*polesj(5)-15.D0*t35*polesi(7)*polesj(5)+0.75D1*v
     #(3)*t7+0.75D1*t33*polesj(17)*polesi(3)+0.75D1*t28*polesj(20)*poles
     #i(4)+0.75D1*t28*polesj(19)*polesi(3)+0.75D1*t38*polesi(12)*polesj(
     #3)+0.15D2*t30*polesj(14)*polesi(3)+0.15D2*t25*polesj(18)*polesi(3)
     #-15.D0*t25*polesi(9)*polesj(8)-15.D0*t30*polesi(6)*polesj(8)+0.15D
     #2*t35*polesj(15)*polesi(3)+0.15D2*t35*polesj(16)*polesi(4)+0.15D2*
     #t25*polesj(19)*polesi(4)-15.D0*t35*polesi(6)*polesj(9)-15.D0*t30*p
     #olesi(7)*polesj(9)+0.15D2*t30*polesj(15)*polesi(4)
      t328 = t38*polesi(21)+2.D0*t25*polesi(34)+2.D0*t35*polesi(28)+2.D0
     #*t25*polesi(25)+t28*polesi(35)+2.D0*t35*polesi(23)+t28*polesi(33)+
     #2.D0*t35*polesi(30)+2.D0*t30*polesi(29)+t33*polesi(33)+t28*polesi(
     #26)+t33*polesi(31)+2.D0*t30*polesi(22)+2.D0*t25*polesi(32)+t38*pol
     #esi(24)+t33*polesi(24)+t38*polesi(26)+2.D0*t30*polesi(27)
      t334 = -15.D0*t25*polesi(10)*polesj(9)+0.75D1*v(1)*t61-15.D0/4.D0*
     #polesi(1)*v(5)-15.D0*t30*polesi(8)*polesj(6)-0.375D1*v(7)*t40-15.D
     #0*t25*polesi(6)*polesj(7)-15.D0*t35*polesi(5)*polesj(7)-15.D0*t35*
     #polesi(9)*polesj(6)-15.D0*t25*polesi(8)*polesj(9)-15.D0*t30*polesi
     #(5)*polesj(6)-15.D0*t25*polesi(7)*polesj(6)-15.D0*t30*polesi(6)*po
     #lesj(5)-15.D0*t35*polesi(10)*polesj(7)-15.D0*t30*polesi(9)*polesj(
     #7)-0.375D1*v(2)*t15-0.375D1*polesj(1)*t328+0.75D1*t28*polesi(20)*p
     #olesj(4)
      t385 = 0.75D1*t33*polesi(17)*polesj(3)+0.75D1*t33*polesi(18)*poles
     #j(4)-15.D0*t35*polesi(7)*polesj(10)-15.D0*t38*polesi(7)*polesj(7)-
     #15.D0*t29*polesj(10)+0.75D1*t28*polesj(16)*polesi(2)+0.75D1*t38*po
     #lesj(11)*polesi(2)-15.D0*t28*polesi(7)*polesj(7)-15.D0*t34*polesj(
     #8)-15.D0*t38*polesi(6)*polesj(6)+0.75D1*t38*polesi(13)*polesj(4)+0
     #.75D1*t38*polesj(13)*polesi(4)-15.D0*t28*polesi(9)*polesj(9)+0.75D
     #1*t33*polesj(14)*polesi(2)+0.75D1*t33*polesi(14)*polesj(2)-15.D0*t
     #33*polesi(9)*polesj(9)-15.D0*t33*polesi(6)*polesj(6)
      t440 = 0.75D1*t38*polesi(11)*polesj(2)+0.15D2*t35*polesj(13)*poles
     #i(2)+0.15D2*t30*polesi(14)*polesj(3)+0.15D2*t25*polesi(19)*polesj(
     #4)+0.15D2*t35*polesi(16)*polesj(4)+0.15D2*t30*polesi(15)*polesj(4)
     #+0.15D2*t30*polesj(12)*polesi(2)+0.15D2*t25*polesj(15)*polesi(2)+0
     #.15D2*t35*polesi(13)*polesj(2)+0.15D2*t30*polesi(12)*polesj(2)-15.
     #D0*t25*polesi(9)*polesj(10)+0.15D2*t35*polesi(15)*polesj(3)+0.15D2
     #*t25*polesi(18)*polesj(3)+0.15D2*t25*polesi(15)*polesj(2)+0.75D1*t
     #28*polesi(16)*polesj(2)+0.75D1*t38*polesj(12)*polesi(3)+0.75D1*t33
     #*polesj(18)*polesi(4)+0.75D1*t28*polesi(19)*polesj(3)
      t447 = t28**2
      t464 = t38**2
      t481 = t33**2
      t486 = 4.D0*t125*Rz*polesi(32)+t447*polesi(35)+4.D0*Rx*t125*polesi
     #(27)+12.D0*t116*Rz*polesi(25)+12.D0*t122*Rz*polesi(28)+4.D0*Ry*t13
     #0*polesi(34)+4.D0*t114*Ry*polesi(22)+t464*polesi(21)+6.D0*t38*t33*
     #polesi(24)+4.D0*Rx*t130*polesi(30)+6.D0*t33*t28*polesi(33)+4.D0*t1
     #14*Rz*polesi(23)+6.D0*t38*t28*polesi(26)+t481*polesi(31)+12.D0*t30
     #*t28*polesi(29)
      t498 = t48**2
      s1 = polesi(1)*polesj(1)/R+(-polesj(1)*t7+polesi(1)*v(1))*t13+(-0.
     #5D0*polesj(1)*t15-0.5D0*polesi(1)*v(7)+polesi(2)*polesj(2)+polesi(
     #3)*polesj(3)+polesi(4)*polesj(4))*t13+(0.15D1*polesj(1)*t40-3.D0*v
     #(1)*t7+0.15D1*v(2)*polesi(1))*t50
      s2 = s1+(0.15D1*polesj(1)*t61-0.15D1*t15*v(1)+0.15D1*v(7)*t7-0.15D
     #1*polesi(1)*v(3)-3.D0*(polesi(9)*polesj(3)+polesi(10)*polesj(4)+po
     #lesi(7)*polesj(2))*Rz-3.D0*(polesi(8)*polesj(3)+polesi(6)*polesj(2
     #)+polesi(9)*polesj(4))*Ry-3.D0*(polesi(7)*polesj(4)+polesi(5)*pole
     #sj(2)+polesi(6)*polesj(3))*Rx+3.D0*(polesj(9)*polesi(3)+polesj(10)
     #*polesi(4)+polesj(7)*polesi(2))*Rz+3.D0*(polesj(8)*polesi(3)+poles
     #j(6)*polesi(2)+polesj(9)*polesi(4))*Ry+3.D0*(polesj(7)*polesi(4)+p
     #olesj(5)*polesi(2)+polesj(6)*polesi(3))*Rx)*t50+(-0.25D1*polesj(1)
     #*t135+0.75D1*t40*v(1)-0.75D1*v(2)*t7+0.25D1*v(4)*polesi(1))*t146
c
      final=s2+(t178+t207)*t50+(t259+t334+t385+t440)*t146+(0.4375D1*pol
     #esj(1)*t486-0.175D2*v(1)*t135+0.2625D2*v(2)*t40-0.175D2*v(4)*t7+0.
     #4375D1*polesi(1)*v(6))/t498/R
c
      return
      end
c
c***********************************************************************
c
      Subroutine get_i_atom(ncs,inx,ics,iatom)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      iatom=inx(2,ics)
      return
      end
c
c***********************************************************************
c
      subroutine transpole(idimp,jdimp,polesi,polesit)
      implicit real*8(a-h,o-z)
      dimension polesi(jdimp,idimp,35),polesit(35,jdimp,idimp)
c
      do icomp=1,idimp
      do ipole=1,35
      do jcomp=1,jdimp
      polesit(ipole,jcomp,icomp) = polesi(jcomp,icomp,ipole)
      enddo
      enddo
      enddo
c
      return
      end

