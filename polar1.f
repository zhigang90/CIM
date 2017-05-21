c===================================================================
c             One-electron part of the dipole polarizability
c
c                    May 2001 , K.Wolinski
c
c===================================================================
      subroutine polar1(rhf,bl,inx,ntri,sxyz,hxyz,fxyz)
      implicit real*8 (a-h,o-z)
      Logical rhf
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c external field info :
      common /fieldx/ xfield,xfgrad,elfiel(9)
c Hessian  options :
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension sxyz(ntri,3),hxyz(ntri,3),fxyz(ntri,3)  ! output
c------------------------------------------------------------------
      call secund(tgrad1)
      call elapsec(elagrad1)
c------------------------------------------------------------------
      ifield=xfield
      ifgrad=xfgrad
c------------------------------------------------------------------
c check-run only or lack of scf convergance :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c------------------------------------------------------------------
c                                  x
c calculate dipole integrals i.e. Hij=(i|x|j) for x=x,y,z
c
      call inton(3,na,hxyz(1,1),inx,1,0,bl(ibas),bl(inuc),ncs)
      call inton(3,na,hxyz(1,2),inx,2,0,bl(ibas),bl(inuc),ncs)
      call inton(3,na,hxyz(1,3),inx,3,0,bl(ibas),bl(inuc),ncs)
c
c     call druma(hxyz(1,1),ncf,iout,'Hxmatrix')
c     call druma(hxyz(1,2),ncf,iout,'Hymatrix')
c     call druma(hxyz(1,3),ncf,iout,'Hzmatrix')
c
c put h1 into f1 for cphf
c
      call dcopy(ntri*3,hxyz(1,1),1,fxyz(1,1),1)
c
ctest only.........................................................
c
c     call xyzdipole(bl,inx,ncf,ncs,na)
c
c------------------------------------------------------------------
c Calculate S1 & T1 (overlap & kinetic) integrals
c S1ij=(i1|j) + (i|j1) and T1ij=(i1|T|j)+(i|T|j1)
c
c------------------------------------------------------------------
c Calculate electron-nuclear attraction contributions
c to the derivative Fock
c V1ij=(i1|V|j)+(i|V|j1)
c
c------------------------------------------------------------------
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) total/60.0d0,elaps/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for 1e part of polariz  ='
     *,f8.2,' Elapsed = ',f8.2,' min')
c------------------------------------------------------------------
      end
c===================================================================
      subroutine xyzdipole(bl,inx,ncf,ncs,na)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension inx(12,*)
      dimension dip(3)
c
       ntri=ncf*(ncf+1)/2
c
      call getival('inuc',inuc)
      call getival('ibas',ibas)
      call getival('ldensi',lden)
c
c  define (temporary) dipole matrix
      call matdef('dipp','s',ncf,ncf)
      call matdef('dens','s',ncf,ncf)
c
      call dcopy(ntri,bl(lden),1,bl(mataddr('dens')),1)
c
      call mmark
c     call immark
c
c  3=dipole
      call inton(3,na,bl(mataddr('dipp')),inx,1,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(1))

      call inton(3,na,bl(mataddr('dipp')),inx,2,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(2))
c
      call inton(3,na,bl(mataddr('dipp')),inx,3,0,bl(ibas),
     1           bl(inuc),ncs)
      call matprodtr('dens','dipp',dip(3))
c
      call retmark
c     call retimark
c
      call matrem('dens')
      call matrem('dipp')
c
c  the electronic part is negative because of the negative charge
      dip(1) = -dip(1)
      dip(2) = -dip(2)
      dip(3) = -dip(3)
c
      write(6,66)' electronic dipole moment =',dip(1),dip(2),dip(3)
c  calculate nuclear dipole and add to dip
      call nucdipole(na,dip,bl(inuc))
c
      write(6,66)' total      dipole moment =',dip(1),dip(2),dip(3)
  66  format(a16,2x,3(f15.9,1x))
      call f_lush(6)
c
c
      end
