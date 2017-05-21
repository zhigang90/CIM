c...
c...  MM March-September 2004
c...
c...  thys file contains routines specific to the calculation
c...  of the radial part of integrals over ab initio
c...  pseudopotentials.
c...  The method of calculation is based on
c...
c...  L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44, 289 (1981)
c...
c...  Altough all the code has been originally written,
c..   this work has been grately facilitated by the the availability
c...  of Fortran routines from the ARGOS code that have been used as
c...  model and for testing.
c...
c=======================================================================
      subroutine DomRad1(nn,npseudo,k,aab,ALPHA,C,len,f1,expa,tolg,maxn,
     $  QQ)
      implicit none
c...
c...  Max2fort (MM) version 0.1 (September 2004)
c...  translated on 9/21/2004 18:12:23
c...
c...  This should be the final strategy for computing type 1 radial integrals.
c...
c...  This function computes all the radial integrals up to a chosen value of n.
c...  The integrals are ultimately generated using the recurrence relation based
c...  on the modified Spherical Bessell function of the first kind:
c...
c...  Q(n,l)=Q(n,l+2)+(2*l+3)/k*Q(n-1,l+1)
c...
c...  This recurrence is quite stable for all values of k^2/(4*alpha)
c...
c...  Using the above recurrence, and nothing than only integrals with n+l even
c...  are needed (as the angular integrals for n+l odd are zero), all the integr
c...  als
c...  can be generated from Q(i,i), i=0,n.
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  See also the comment lines of subroutine DoQ12f.
c...
c...  The starting values for the recurrence are computed as follow:
c...
c...  Case 1):  k^2/(4*alpha) <= 5.0
c...     1.a)   compute Q(n,n) and Q(n-1,n-1).
c...     1.b)   Apply recurrence F of McMurchie and Davidson paper to obtain the
c...   remaining
c...            integrals of type Q(i,i).
c...
c...  Case 2):  k^2/(4*alpha) > 5.0
c...     2.a)   compute Q(0,0) and Q(1,1).
c...     2.b)   Apply recurrence E of McMurchie and Davidson paper to obtain the
c...   remaining
c...            integrals of type Q(i,i).
c...
c...  Recurrence E is not stable for small values of K^2/(4*alpha) or large valu
c...  es of n,
c...  while recurrence F is not stable for large k^2/(4*alpha).
c...
c...  The Algorithm outlined above is used to generate the starting integrals fo
c...  r n up to and including 14.
c...  If n > 14 (which should never occur in routine calculations), the starting
c...   integrals
c...  Q(j,j) for 14 < j <= n are computed one by one.
c...
c...  The thresholds have been chosen a bit more conservatively than in McMurchi
c...  e and
c...  Davidson, and should ensure a global relative accuracy of 10^(-14)-10^(-15
c...  ).
c...  Input parameters:
c...
c...  nn          maximum value of the power of r coming from the cartesian term
c...  s
c...  npseudo     array of contributions to the power or r from the pseudopotent
c...  ial.
c...  k           positional parameter (depends only on the Gaussian part of the
c...   integral
c...              and on the center of the pseudopotential)
c...  aab         part of the exponential coefficient coming from
c...              the Gaussians
c...  alpha       array of coefficients of the exponential pseudopotential term
c...  c           array of contraction coefficients of the pseudopotential term
c...  len         length of the contraction of the pseudopotential term
c...  f1          multiplicative factor
c...  maxn        maximum allowed value of n (for dimensioning)
c...
c...  Output parameter:
c...
c...  qq          array of values of type one radial integrals
c...
      real*8 ZERO,one,two,three,four
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,four=4.0D0)
c...
      integer nn,len,maxn,npseudo(len),n,np,l,nr,nc
      real*8 k,k2,aab,ALPHA(len),C(len),QQ(0:maxn,0:maxn),f1,expa,tolg,z
     $  ,kk,kkm1,a,A2,a2k,a2m1,pl,pp,qq1
c...
      if(len.le.0)
     $  call nerror(1,'DomRad1','Improper value of parameter len',len,0)
c
      do np=0,nn,1
        do l=0,np,1
          QQ(np,l)=ZERO
        enddo
      enddo
c...
c...  if k=0, only the qq(np,0),np=0,nn,2 integrals survive
c...
      k2=k**2
      if(k.eq.ZERO)then
        if(len.eq.1)then
          z=2.5D-1*k2/(aab+ALPHA(1))
          if(z-expa.gt.tolg)then
            do np=0,nn,2
              call DoQ12f(np+npseudo(1),0,k,aab+alpha(1),qq1)
              qq(np,0)=qq1*c(1)*f1
            enddo
          endif
        else
          do nc=1,len,1
            z=2.5D-1*k2/(ALPHA(nc)+aab)
            if(z-expa.gt.tolg)then
              do np=0,nn,2
                call DoQ12f(np+npseudo(nc),0,k,aab+alpha(nc),qq1)
                qq(np,0)=qq(np,0)+qq1*c(nc)*f1
              enddo
            endif
          enddo
        endif
      else
        n=MIN(nn,14)
c...
c...  if nn < 2 we do not need any recurrence
c...
        if(n.lt.2)then
          if(len.eq.1)then
            z=2.5D-1*k2/(aab+ALPHA(1))
            if(z-expa.le.tolg)then
              goto 1
            endif
            call DoQ12f(npseudo(1),0,k,aab+alpha(1),qq1)
            qq(0,0)=qq1*c(1)*f1
            if(n.eq.1)then
              call DoQ12f(npseudo(1)+1,1,k,aab+alpha(1),qq1)
              qq(1,1)=qq1*c(1)*f1
            endif
 1          continue
          else
            do nc=1,len,1
              z=2.5D-1*k2/(ALPHA(nc)+aab)
              if(z-expa.le.tolg)then
                goto 2
              endif
              call DoQ12f(npseudo(nc),0,k,aab+alpha(nc),qq1)
              qq(0,0)=qq(0,0)+qq1*c(nc)*f1
              if(n.eq.1)then
                call DoQ12f(npseudo(nc)+1,1,k,aab+alpha(nc),qq1)
                qq(1,1)=qq(1,1)+qq1*c(nc)*f1
              endif
 2            continue
            enddo
          endif
        else
c...
c...  compute starting integrals for n up to 14
c...
          kk=k
          kkm1=one/kk
          if(len.eq.1)then
            nr=n+npseudo(1)
            a=aab+ALPHA(1)
            z=2.5D-1*k2/a
            if(z-expa.le.tolg)then
              goto 3
            endif
            A2=a*two
            a2k=A2/kk
            a2m1=one/A2
c...
c...  Case 1: power series plus recurrence F
c...
            if(z.le.5.0D0)then
c...
c...  Compute qq[n,n], and qq[n-1,n-1] by power series
c...
              call DoQ12f(nr,n,k,a,qq1)
              qq(n,n)=qq1*c(1)*f1
              call DoQ12f(nr-1,n-1,k,a,qq1)
              qq(n-1,n-1)=qq1*c(1)*f1
c...
c...  Apply recurrence F
c...
              pp=n*two-1.0D0*two
              do np=n-2,1,-1
                pp=pp-1.0D0*two
                QQ(np,np)=one*(A2*QQ(np+2,np+2)-1.0D0*QQ(np+1,np+1)*(kk-
     $            1.0D0*a2k*(three+pp)))/(pp+one+npseudo(1))
              enddo
              QQ(0,0)=one*(QQ(2,2)*A2-1.0D0*QQ(1,1)*(kk-1.0D0*a2k*three)
     $          )/(npseudo(1)+1.0D0)
c...
c...  Case 2: power series (or asymptotic series) plus recurrence E
c...
            else
c...
c...  Compute qq[0,0] and qq[1,1]
c...
              call DoQ12f(npseudo(1),0,k,a,qq1)
              qq(0,0)=qq1*c(1)*f1
              call DoQ12f(npseudo(1)+1,1,k,a,qq1)
              qq(1,1)=qq1*c(1)*f1
c...
c...  Apply recurrence E
c...
              QQ(2,2)=a2m1*(QQ(1,1)*(kk-1.0D0*a2k*three)+QQ(0,0)*(npseud
     $          o(1)+1.0D0))
              pp=four
              do np=3,n,1
                pp=two+pp
                QQ(np,np)=a2m1*(QQ(np-2,np-2)*(-1.0D0*three+pp+npseudo(1
     $            ))+QQ(np-1,np-1)*(kk-1.0D0*a2k*(pp-1.0D0*one)))
              enddo
            endif
 3          continue
          else
            do nc=1,len,1
              nr=npseudo(nc)+n
              a=ALPHA(nc)+aab
              z=2.5D-1*k2/a
              if(z-expa.le.tolg)then
                goto 4
              endif
              A2=a*two
              a2k=A2/kk
              a2m1=one/A2
c...
c...  Case 1: power series plus recurrence F
c...
              if(z.le.5.0D0)then
c...
c...  Compute qq[n,n], and qq[n-1,n-1] by power series
c...
                call DoQ12f(nr,n,k,a,qq1)
                qq(n,0)=qq1*c(nc)*f1
                call DoQ12f(nr-1,n-1,k,a,qq1)
                qq(n-1,0)=qq1*c(nc)*f1
c...
c...  Apply recurrence F
c...
                pp=n*two-1.0D0*two
                do np=n-2,1,-1
                  pp=pp-1.0D0*two
                  QQ(np,0)=one*(A2*QQ(np+2,0)-1.0D0*QQ(np+1,0)*(kk-1.0D0
     $              *a2k*(three+pp)))/(pp+one+npseudo(nc))
                enddo
                QQ(0,0)=one*(QQ(2,0)*A2-1.0D0*QQ(1,0)*(kk-1.0D0*a2k*thre
     $            e))/(npseudo(nc)+1.0D0)+QQ(0,0)
                do np=1,n,1
                  QQ(np,np)=QQ(np,np)+QQ(np,0)
                enddo
c...
c...  Case 2: power series (or asymptotic series) plus recurrence E
c...
              else
c...
c...  Compute qq[0,0] and qq[1,1]
c...
                call DoQ12f(npseudo(nc),0,k,a,qq1)
                qq(n,0)=qq1*c(nc)*f1
                call DoQ12f(npseudo(nc)+1,1,k,a,qq1)
                qq(n,1)=qq1*c(nc)*f1
c...
c...  Apply recurrence E
c...
                QQ(n,2)=a2m1*(QQ(n,1)*(kk-1.0D0*a2k*three)+QQ(n,0)*(npse
     $            udo(nc)+1.0D0))
                pp=four
                do np=3,n-1,1
                  pp=two+pp
                  QQ(n,np)=a2m1*(QQ(n,np-2)*(-1.0D0*three+pp+npseudo(nc)
     $              )+QQ(n,np-1)*(kk-1.0D0*a2k*(pp-1.0D0*one)))
                enddo
                pp=two+pp
                QQ(n,n)=a2m1*(QQ(n,n-2)*(-1.0D0*three+pp+npseudo(nc))+QQ
     $            (n,n-1)*(kk-1.0D0*a2k*(pp-1.0D0*one)))+QQ(n,n)
                do np=0,n-1,1
                  QQ(np,np)=QQ(np,np)+QQ(n,np)
                enddo
              endif
 4            continue
            enddo
          endif
c...
c...  if nn > 14, compute the remaining starting integrals one by one
c...
          if(nn.gt.14)then
            if(len.eq.1)then
              z=2.5D-1*k2/(aab+ALPHA(1))
              if(z-expa.gt.tolg)then
                do np=15,nn,1
                  call DoQ12f(npseudo(1)+np,np,k,aab+alpha(1),qq1)
                  qq(np,np)=qq1*c(1)*f1
                enddo
              endif
            else
              do nc=1,len,1
                z=2.5D-1*k2/(ALPHA(nc)+aab)
                if(z-expa.le.tolg)then
                  goto 5
                endif
                do np=15,nn,1
                  call DoQ12f(npseudo(nc)+np,np,k,aab+alpha(nc),qq1)
                  qq(np,np)=qq(np,np)+qq1*c(nc)*f1
                enddo
 5              continue
              enddo
            endif
          endif
c...
c...  All the starting integrals have been generated, now apply the final recurr
c...  ence
c...
          do np=2,nn,2
            pl=one
            do l=0,nn-np,1
              pl=two+pl
              QQ(np+l,l)=kkm1*QQ(np+l-1,l+1)*pl+QQ(np+l,l+2)
            enddo
          enddo
        endif
      endif
c...
      end
c=======================================================================
      subroutine DoQ12f(n,l,k,ALPHA,q1)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 8/5/2004 19:29:52
c...
c...  This subroutine computes a type 1 Radial integral.
c...  The target relative accuracy is 10^(-14)--10^(-15).
c...
c...                        INF
c...                       /                            2
c...                       [               n   - ALPHA r
c...  Q1(n, l, k, ALPHA) = I    M(l, k r) r  %E           dr
c...                       ]
c...                       /
c...                        0
c...
c...  Where M is a modified spherical Bessel function of the first kind
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (31-36)
c...
c...  **** This version works for n up to 30 ****
c...
c...  The method of computation depends on the input parameters.
c...
c...  zlist contains the minimum values of z=k^2/(4*alpha) after which,
c...  depending on the value of n, the asymptotic series provides the
c...  desired accuracy. This list is more conservative than in McMurchie
c...  and Davidson paper, where the target accuracy is 10^(-12).
c...
c...  if z >= zlist(min(n,15)), the integral is computed with the
c...  asymptotic series, eq (34)
c...
c...  if z < zlist(min(n+1,15)) and n+l even with n >= l+2 the alternating
c...  series is used, eq (36), as it truncates.
c...
c...  otherwise the power series is used, eq (33)
c...
c...  Input parameters:
c...
c...  n           exponent of r
c...  l           degree of the mofified spherical Bessel function of the
c...              first kind
c...  k           positional parameter
c...  alpha       coefficient of the exponential term
c...
c...  Output parameter:
c...
c...  q1          value of the type one radial integral
c...
      integer i2
      parameter (i2=2)
      real*8 ZERO,one,two,half
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,half=5.0D-1)
c...
      integer n,l
      real*8 k,ALPHA,q1,z,pre,pre1,SUM,sum1,term,termp,zppf,conv,a,b,p,d
     $  f1,df2,ddf1,ddf2,kk,NUM,den,ndf,ddf
c...
c...
      real*8 zlist(0:15)
      save zlist
      data zlist/3.6D1,3.6D1,3.30625D1,3.30625D1,3.025D1,2.75625D1,2.756
     $  25D1,2.5D1,2.5D1,2.25625D1,2.25625D1,2.025D1,2.025D1,1.80625D1,1
     $  .80625D1,1.6D1/
c...
      real*8 dbf(-1:61)
      save dbf
      data dbf/1.0D0,1.0D0,1.0D0,2.0D0,3.0D0,8.0D0,1.5D1,4.8D1,1.05D2,3.
     $  84D2,9.45D2,3.84D3,1.0395D4,4.608D4,1.35135D5,6.4512D5,2.027025D
     $  6,1.032192D7,3.4459425D7,1.8579456D8,6.54729075D8,3.7158912D9,1.
     $  3749310575D10,8.17496064D10,3.16234143225D11,1.9619905536D12,7.9
     $  05853580625D12,5.10117543936D13,2.13458046676875D14,1.4283291230
     $  208D15,6.190283353629375D15,4.2849873690624D16,1.918987839625106
     $  D17,1.371195958099968D18,6.332659870762851D18,4.662066257539891D
     $  19,2.216430954766998D20,1.678343852714361D21,8.200794532637891D2
     $  1,6.377706640314571D22,3.198309867728778D23,2.551082656125828D24
     $  ,1.311307045768799D25,1.071454715572848D26,5.638620296805835D26,
     $  4.714400748520531D27,2.537379133562626D28,2.168624344319444D29,1
     $  .192568192774434D30,1.040939685273333D31,5.843584144594727D31,5.
     $  204698426366666D32,2.980227913743311D33,2.706443181710666D34,1.5
     $  79520794283955D35,1.46147931812376D36,8.687364368561751D36,8.184
     $  284181493055D37,4.951797690080198D38,4.746884825265972D39,2.9215
     $  60637147317D40,2.848130895159583D41,1.782151988659863D42/
c...
      real*8 ff(-1:61)
      save ff
      data ff/1.0D0,1.0D0,1.0D0,2.0D0,6.0D0,2.4D1,1.2D2,7.2D2,5.04D3,4.0
     $  32D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,8.717829
     $  12D10,1.307674368D12,2.0922789888D13,3.55687428096D14,6.40237370
     $  5728D15,1.21645100408832D17,2.43290200817664D18,5.10909421717094
     $  4D19,1.124000727777608D21,2.585201673888498D22,6.204484017332394
     $  D23,1.551121004333099D25,4.032914611266056D26,1.088886945041835D
     $  28,3.048883446117139D29,8.841761993739702D30,2.652528598121911D3
     $  2,8.222838654177923D33,2.631308369336935D35,8.683317618811887D36
     $  ,2.952327990396041D38,1.033314796638614D40,3.719933267899012D41,
     $  1.376375309122634D43,5.230226174666011D44,2.039788208119744D46,8
     $  .159152832478977D47,3.345252661316381D49,1.40500611775288D51,6.0
     $  41526306337384D52,2.658271574788449D54,1.196222208654802D56,5.50
     $  2622159812089D57,2.586232415111682D59,1.241391559253607D61,6.082
     $  818640342676D62,3.041409320171338D64,1.551118753287382D66,8.0658
     $  17517094388D67,4.274883284060026D69,2.308436973392414D71,1.26964
     $  0335365828D73,7.109985878048635D74,4.052691950487722D76,2.350561
     $  331282879D78,1.386831185456898D80,8.32098711274139D81,5.07580213
     $  8772248D83/
c...
      if(n.gt.30.or.l.gt.30)
     $  call nerror(1,'DoQ12f','Maximum value of n or l exceeded',
     $              MAX(n,l),30)
c...
      z=2.5D-1*k**2/ALPHA
c...
c...  if k = 0 and l # 0, the integral is zero
c...
      if((k.eq.ZERO.and.l.ne.0))then
        q1=ZERO
      else
        if(z.ge.zlist(MIN(n,15)))then
c...
c...  Asymptotic series
c...
          zppf=one
          a=half*(n+l+1.0D0)
          b=half*(2.0D0*l+3.0D0)
          ddf1=one
          ddf2=one
          sum1=one
          if(mod(n+l,2).ne.0)then
            pre=one
            pre1=1.772453850905516D0*dbf(2*l+1)*two**(-l-1)/ff((n+l-1)/i
     $        2)
          else
            pre=1.2533141373155D0
            pre1=dbf(2*l+1)*two**(-one-half*(l-n))/dbf(n+l-1)
          endif
          pre=k**l*dbf(n+l-1)*pre/(dbf(2*l+1)*(ALPHA*two)**(half*(n+l+1)
     $      ))
          pre1=pre1*z**(a-b)*EXP(z)
          df1=-a
          df2=-one+b-a
          p=ZERO
          termp=one
 1        continue
            p=p+one
            zppf=zppf/(p*z)
            df1=one+df1
            df2=one+df2
            if((df1.eq.ZERO.or.df2.eq.ZERO))then
              goto 2
            endif
            ddf1=ddf1*df1
            ddf2=ddf2*df2
            term=ddf1*ddf2*zppf
            if(l+1.ge.n)then
              if(ABS(term).lt.1.0d-14.or.ABS(term).gt.ABS(termp))then
                goto 2
              else
                if(p.gt.500.0) call nerror(2,'DoQ12f',
     $                  'Asymptotic series Not converged',0,0)
              endif
            endif
            sum1=term+sum1
            termp=term
            goto 1
 2        continue
          q1=pre*pre1*sum1
        else
          if((mod(n+l,2).eq.0.and.n.gt.l+1))then
c...
c...  Alternating series (used only if it terminates)
c...
            zppf=one
            df1=half*(two-n+l)-one
            df2=half*(l*two+3.0D0)-one
            ddf1=one
            ddf2=one
            SUM=one
            if(k.ne.ZERO)then
              pre=1.2533141373155D0*k**l*dbf(n+l-1)*EXP(z)/(dbf(2*l+1)*(
     $          ALPHA*two)**(half*(n+l+1)))
            else
              pre=1.2533141373155D0*dbf(n+l-1)*EXP(z)/(dbf(2*l+1)*(ALPHA
     $          *two)**(half*(n+l+1)))
            endif
            p=ZERO
            if(z.ne.ZERO)then
 3            continue
                p=p+one
                zppf=-z*zppf/p
                df1=one+df1
                df2=one+df2
                if(df1.eq.ZERO)then
                  goto 4
                endif
                ddf1=ddf1*df1
                ddf2=ddf2*df2
                SUM=ddf1*zppf/ddf2+SUM
                goto 3
 4            continue
            endif
            q1=pre*SUM
          else
c...
c...  Power series
c...
            kk=k
            a=ALPHA
            zppf=one
            ndf=-one+n+l
            ddf=l*two+one
            NUM=one
            den=one
            SUM=one
            if(mod(n+l,2).ne.0)then
              pre=one
            else
              pre=1.2533141373155D0
            endif
            if(kk.ne.ZERO)then
              pre=kk**l*dbf(n+l-1)*pre/(dbf(2*l+1)*(a*two)**(half*(n+l+1
     $          )))
            else
              pre=dbf(n+l-1)*pre/(dbf(2*l+1)*(a*two)**(half*(n+l+1)))
            endif
            p=ZERO
            if(z.ne.ZERO)then
 5            continue
                p=p+one
                zppf=z*zppf/p
                ddf=two+ddf
                den=ddf*den
                ndf=two+ndf
                NUM=ndf*NUM
                term=NUM*zppf/den
                SUM=term+SUM
                conv=term/SUM
                if(ABS(conv).lt.1.0d-14)then
                  goto 6
                else
                  if(p.gt.500.0) call nerror(3,'DoQ12f',
     $                  'Power series Not converged',0,0)
                endif
                goto 5
 6            continue
            endif
            q1=pre*SUM
          endif
        endif
      endif
c...
      end
c=======================================================================
      subroutine recq1(ll,nps,nini,nfin,z,maxn,q1)
      implicit none
c...
c...  MM july 2004
c...
c...  (partially) fills a table of type one radial integrals
c...  using the upward recurrence for n
c...
c...  McMurchie and Davidson, J. Comp. Phys, 44, 289 (1981), eq. (52)
c...  (note that here alpha=1)
c...
c...  This is used to generate Q1 integrals for the single power
c...  series calculation of type 2 radial integrals
c...
c...  ll           maximum l value
c...  nps          contribution to exponent from pseudopotential
c...  nini         inital n for recurrence
c...  nfin         final n for recurrence
c...  z            positional parameter (k in the paper)
c...  maxn         maximum n value, for dimensioning
c...  q1           table of q1 values
c...
      integer ll,nini,nfin,maxn,nps
      real*8 z,q1(0:maxn,0:ll)
      integer n,l,nn
      real*8 z2fam1,nm4,nm3
c...
      real*8 ZERO,one,two,three,four,fourth
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,four=4.0D0)
      parameter (fourth=0.25d0)
c...
      if(nfin.gt.maxn)
     $  call nerror(1,'recq1','Maximum n value exceeded:',nfin,maxn)
c...
      if(nini.le.3)then
        do n=nini,3
          do l=0,ll
            if(mod(n+l,2).eq.0) call Doq12f(n+nps,l,z,one,q1(n,l))
          enddo
        enddo
        nn=4
      else
        nn=nini
      endif
      z2fam1=z*z*fourth+dfloat(nn+nps-1)-2.5d0
      do n=nn,nfin
        z2fam1=z2fam1+one
        nm4=dfloat(-(n+nps)+3)
        nm3=dfloat(n+nps-4)
        do l=0,ll
          nm4=nm4+one
          nm3=nm3+one
          if(mod(n+l,2).eq.0)then
            q1(n,l)=z2fam1*q1(n-2,l)+fourth*nm4*nm3*q1(n-4,l)
          endif
        enddo
      enddo
      end
c=======================================================================
      subroutine csumterm(ini,ifi,z,cst,maxst)
      implicit none
c...
c...  MM july 2004
c...
c...  precalculation of terms for the single power
c...  series calculation of type 2 radial integrals
c...
c...             2i
c...            z
c...  cst(i)= -------
c...             i
c...            2  i!
c...
      integer ini,ifi,maxst,i,ii
      real*8 z,cst(0:maxst),z2,p,st
c...
      if(ifi.gt.maxst)
     $  call nerror(1,'Csumterm','Maximum number of terms exceeded',
     $                ifi,maxst)
c
      if(ini.lt.0)
     $  call nerror(2,'Csumterm','Improper call',ini,0)
c
      z2=0.5d0*z*z
      if(ini.eq.0)then
        cst(0)=1.0d0
        ii=1
      else
        ii=ini
      endif
      st=cst(ii-1)
      p=dfloat(ii-1)
      do i=ii,ifi,1
        p=p+1.0d0
        st=st*z2/p
        cst(i)=st
      enddo
c...
      end
c=======================================================================
      subroutine SumSpn(n,nps,l1,l2,z1,z2,nadd,nst,maxst,cst,nfm,
     $                  lq1,maxnq1,q1,SUM)
      implicit none
c...
c...  MM July 2004
c...
c...  Single powers series summation for type 2 radial integrals.
c...
c...                          INF
c...                          ====       2p
c...                          \        z1   Q1(n+2p+l1,l2,z2,1)
c...  SumSp(n,l1,l2,z1,z2)=    >      --------------------------
c...                          /         p
c...                          ====     2  p! (2 p + 2 l1 + 1)!!
c...                          p = 0
c...
c...  Where Q1 is a type 1 radial integral that can be computed by recurrence:
c...
c...                                2
c...                     2 n - 5   k
c...  Q1(n, l, k, a) = ((------- + ---) Q1(n - 2, l, k, a)
c...                        2      4 a
c...
c...                    (- n + l + 4) (n + l - 3) Q1(n - 4, l, k, a)
c...                  + --------------------------------------------)/a
c...                                        4 a
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (51-52)
c...
c...  this version uses tables of precomputed summation terms and
c...  tupe one radial integrals
c...
c...  **** This version works for l1 up to 15 ****
c...
c...
      real*8 ZERO,one,two,three,four,half,fh,epsi
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,four=4.0D0,h
     $  alf=5.0D-1,fh=2.5D0,epsi=1.0D-14)
c...
      integer n,nps,l1,l2,nst,maxst,nfm,lq1,maxnq1,nadd
      real*8 z1,z2,SUM,cst(0:maxst),q1(0:maxnq1,0:lq1)
      integer nn,i,il,nmore,nterm
      real*8 p,z,ddf,qn,qnm2,qnm4,lmnp4,lpnm3,tnm52,NUM,den,term,conv
c...
      real*8 dbf(-1:41)
      save dbf
      data dbf/1.0D0,1.0D0,1.0D0,2.0D0,3.0D0,8.0D0,1.5D1,4.8D1,1.05D2,3.
     $  84D2,9.45D2,3.84D3,1.0395D4,4.608D4,1.35135D5,6.4512D5,2.027025D
     $  6,1.032192D7,3.4459425D7,1.8579456D8,6.54729075D8,3.7158912D9,1.
     $  3749310575D10,8.17496064D10,3.16234143225D11,1.9619905536D12,7.9
     $  05853580625D12,5.10117543936D13,2.13458046676875D14,1.4283291230
     $  208D15,6.190283353629375D15,4.2849873690624D16,1.918987839625106
     $  D17,1.371195958099968D18,6.332659870762851D18,4.662066257539891D
     $  19,2.216430954766998D20,1.678343852714361D21,8.200794532637891D2
     $  1,6.377706640314571D22,3.198309867728778D23,2.551082656125828D24
     $  ,1.311307045768799D25/
c...
      if(l1.gt.20)
     $  call nerror(1,'Sumsp','Maximum l1 value exceeded ',l1,20)
c
      nn=n+l1
c...
c...  the first two terms of the summation are computed outside the loop
c...
      den=dbf(2*l1+1)
      SUM=q1(nn,l2)/den
      nn=nn+2
      ddf=dfloat(l1)*two+three
      den=ddf*den
      SUM=SUM+cst(1)*q1(nn,l2)/den
      nterm=min(nst,int((nfm-nn)/2))
c...
c...  loop over the precomputed terms
c...
      do i=2,nterm
        nn=nn+2
        ddf=ddf+two
        den=den*ddf
        term=cst(i)*q1(nn,l2)/den
        SUM=term+SUM
        conv=term/SUM
        if(ABS(conv).lt.epsi) return
      enddo
c...
c...  not converged. See if we can precompute some more
c...
1     il=nterm
      nmore=min(nadd,maxst-il,int((maxnq1-nn)/2))
      if(nmore.gt.0)then
        if(il+nmore.gt.nst)then
           call csumterm(nst+1,il+nmore,z1,cst,maxst)
           nst=il+nmore
        endif
        if(nn+2*nmore.gt.nfm)then
           call recq1(lq1,nps,nfm+1,nn+2*nmore,z2,maxnq1,q1)
           nfm=nn+2*nmore
        endif
        nterm=il+nmore
        do i=il+1,nterm
          nn=nn+2
          ddf=ddf+two
          den=den*ddf
          term=cst(i)*q1(nn,l2)/den
          SUM=term+SUM
          conv=term/SUM
          if(ABS(conv).lt.epsi) return
        enddo
        goto 1
      endif
c...
c...  We cannot precompute anymore. Let's continue the hard way
c...
      z=half*z1**2
      qnm2=q1(nn,l2)
      qnm4=q1(nn-2,l2)
      lmnp4=dfloat(l2-(nn+nps)+4)*0.25d0
      lpnm3=dfloat(l2+nn+nps-3)
      tnm52=z2**2*0.25d0+dfloat(nn+nps)-fh
      p=dfloat(il)
      num=cst(il)
      do i=il+1,500
        p=p+one
        lmnp4=lmnp4-half
        lpnm3=two+lpnm3
        tnm52=two+tnm52
        ddf=two+ddf
        den=ddf*den
        num=num*z/p
c...
        qn=qnm2*tnm52+lmnp4*lpnm3*qnm4
        qnm4=qnm2
        qnm2=qn
c...
        term=NUM*qn/den
        SUM=term+SUM
        conv=term/SUM
        if(ABS(conv).lt.epsi) return
      enddo
c...
c...  if we made it to here, it is hopeless
c...
      call nerror(2,'Sumspn','Sum not converged after 500 terms',0,0)
c...
      end
c=======================================================================
      subroutine sbess1pn(l,x,n,sbess1p)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 9/15/2004 15:1:16
c...
c...  Modified spherical bessel function of the first kind
c...  divided by the exponential of the argument.
c...  This will be used for computing the radial integrals
c...  by Gauss-Hermite quadrature.
c...
c...  McMurchie and Davidson, eqs. (39,53-54) or
c...  Abramovitz and Stegun 10.2.5 and 10.2.9,10.2.11, pag 443
c...
c...  for x < 11.0 the ascending series is used (39 or 10.2.5),
c...  otherwise the sums 53-54 or 10.2.9,10.2.11 are used.
c...  this threshold is more conservative than in McMurchie and
c...  Davidson paper, and should ensure a relative accuracy
c...  greater than 1E-14
c...
      integer maxpt
      parameter (maxpt=20)
      real*8 ZERO,one,two,half,epsi
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,half=5.0D-1,epsi=1.0D-14
     $  )
c...
      integer l,k,n,i,nast,conv
      real*8 x(n),sbess1p(n),x2(maxpt),xp(maxpt),xm(maxpt),p,pf,pfm1,tp1
     $  ,sum1(maxpt),sum2(maxpt)
c...
c...
c...  dbfm1(i) = 1/i!!
c...
c...
      real*8 dbfm1(-1:41)
      save dbfm1
      data dbfm1/1.0D0,1.0D0,1.0D0,5.0D-1,3.333333333333333D-1,1.25D-1,6
     $  .666666666666667D-2,2.083333333333333D-2,9.523809523809524D-3,2.
     $  604166666666667D-3,1.058201058201058D-3,2.604166666666667D-4,9.6
     $  2000962000962D-5,2.170138888888889D-5,7.4000074000074D-6,1.55009
     $  9206349206D-6,4.9333382666716D-7,9.68812003968254D-8,2.901963686
     $  277412D-8,5.382288910934744D-9,1.527349308567059D-9,2.6911444554
     $  67372D-10,7.273091945557423D-11,1.223247479757896D-11,3.16221388
     $  9372793D-12,5.096864498991235D-13,1.264885555749117D-13,1.960332
     $  499612014D-14,4.684761317589322D-15,7.001187498614334D-16,1.6154
     $  34937099766D-16,2.333729166204778D-17,5.211080442257311D-18,7.29
     $  2903644389931D-19,1.579115285532518D-19,2.144971660114686D-20,4.
     $  511757958664338D-21,5.958254611429682D-22,1.219394042882254D-22,
     $  1.567961739849916D-23,3.126651392005778D-24,3.919904349624791D-2
     $  5,7.625979004892143D-26/
c...
c...  f(i) = i!
c...
c...
      real*8 f(-1:41)
      save f
      data f/1.0D0,1.0D0,1.0D0,2.0D0,6.0D0,2.4D1,1.2D2,7.2D2,5.04D3,4.03
     $  2D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,8.7178291
     $  2D10,1.307674368D12,2.0922789888D13,3.55687428096D14,6.402373705
     $  728D15,1.21645100408832D17,2.43290200817664D18,5.109094217170944
     $  D19,1.124000727777608D21,2.585201673888498D22,6.204484017332394D
     $  23,1.551121004333099D25,4.032914611266056D26,1.088886945041835D2
     $  8,3.048883446117139D29,8.841761993739702D30,2.652528598121911D32
     $  ,8.222838654177923D33,2.631308369336935D35,8.683317618811887D36,
     $  2.952327990396041D38,1.033314796638614D40,3.719933267899012D41,1
     $  .376375309122634D43,5.230226174666011D44,2.039788208119744D46,8.
     $  159152832478977D47,3.345252661316381D49/
c...
      if(l.gt.20)
     $  call nerror(1,'Sbess1pn','Maximum l value exceeded',l,20)
c
      if(n.gt.maxpt)
     $  call nerror(2,'Sbess1pn','Cannot compute more points than',
     $              maxpt,0)
c
      nast=1
      do i=1,n,1
        if(x(i).le.11.0)then
          nast=0
        endif
      enddo
      if(nast.eq.0)then
c...
c...  ascending series
c...
        do i=1,n,1
          x2(i)=half*x(i)**2
          xp(i)=one
          sum1(i)=one
        enddo
        pf=one
        tp1=l*two+one
        p=ZERO
 1      continue
          p=p+one
          tp1=two+tp1
          pf=p*pf*tp1
          pfm1=one/pf
          conv=1
          do i=1,n,1
            xp(i)=x2(i)*xp(i)
            sum2(i)=xp(i)*pfm1
            sum1(i)=sum2(i)+sum1(i)
            if(ABS(sum2(i))/ABS(sum1(i)).gt.epsi)then
              conv=0
            endif
          enddo
          if(conv.eq.1)then
            goto 2
          else
            if(p.gt.500.0) call nerror(1,'Sbess1pn',
     $              'Sum Not converged after 500 terms',0,0)
          endif
          goto 1
 2      continue
        do i=1,n,1
          sbess1p(i)=sum1(i)*x(i)**l*EXP(-1.0D0*x(i))*dbfm1(2*l+1)
        enddo
      else
c...
c...  truncated sum
c...
        do i=1,n,1
          x2(i)=half/x(i)
          xm(i)=one
          xp(i)=one
          sum1(i)=one
          sum2(i)=one
        enddo
        do k=1,l,1
          pf=f(l+k)/(f(k)*f(l-k))
          do i=1,n,1
            xm(i)=-1.0D0*x2(i)*xm(i)
            xp(i)=x2(i)*xp(i)
            sum1(i)=xm(i)*pf+sum1(i)
            sum2(i)=xp(i)*pf+sum2(i)
          enddo
        enddo
        if(mod(l,2).eq.0)then
          do i=1,n,1
            sbess1p(i)=x2(i)*(sum1(i)-1.0D0*sum2(i)*EXP(-1.0D0*x(i)*two)
     $        )
          enddo
        else
          do i=1,n,1
            sbess1p(i)=x2(i)*(sum2(i)*EXP(-1.0D0*x(i)*two)+sum1(i))
          enddo
        endif
      endif
c...
      end
c=======================================================================
      subroutine DomRad2n(nmaxa,nmaxb,lp,npseudo,ka,kb,aa,ab,ALPHA,C,
     $                    len,maxn,maxl,fspm1,tol,QQ)
      implicit none
c...
c...  MM july 2004
c...
c...  This subroutine computes a table of type 2 radial integrals.
c...
c...                          INF
c...                         /                                      2
c...                         [                         n   - alpha r
c...Q2(n,la,lb,ka,kb,alpha)= I M(la,ka r) M(lb,kb r)  r  %E          dr
c...                         ]
c...                         /
c...                          0
c...
c...  Where M are modified spherical Bessel functions of the first kind,
c...
c...  The table is filled using the recurrence relation for the modified
c...  spherical Bessel function of the first kind which, applied
c...  to type 2 radial integrals, become:
c...
c...  Q2(n,la,lb,ka,kb,a)=Q2(n,la+2,lb,ka,kb,a)+
c...                              (2*la+3)/ka*Q2(n-1,la+1,lb,ka,kb,a)
c...
c...  (the same applies to lb and kb)
c...
c...  the output values are scaled by
c...
c...  exp(-ka^2/(4*aa)-kb^2/(4*ab))=exp(-aa*CA^2-ab*CB^2)
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (25,39-60)
c...
c...  there can be the following cases:
c...
c...  Case 1: ka = 0.0 and kb = 0.0
c...
c...     only the integrals Q2(n,0,0)=Q1(n,0) n=0,nmaxa+nmaxb,
c...     n+lp=even survive. The integrals are computed as type 1.
c...     No recurrence has to be applied.
c...
c...  Case 2: ka # 0.0 and kb = 0.0
c...
c...     only the integrals Q2(na+nb,la,0)=Q1(na+nb,la),for na=0,nmaxa,
c...     la=0,na+lp, nb=0,nmaxb,2 and na+lp+la=even survive.
c...     the integrals are computed as type one. The table is
c...     filled by recurrence on a
c...
c...  Case 3: ka = 0.0 and kb # 0.0
c...
c...     same as case 2, but with a replaced by b.
c...
c...  Case 4: ka # 0.0 and kb # 0.0
c...                                                ka + kb
c...     the method depends on the value of rc = -------------
c...                                             2 sqrt(alpha)
c...              2
c...            rc  < 50     single power series (eq. 51)
c...
c...              2
c...      50 <  rc  < 500    20 points Gauss-Hermite quadrature (53-60)
c...
c...              2
c...      500 <  rc < 50000  10 points Gauss-Hermite quadrature
c...
c...                 2
c...      50000 <  rc         5 points Gauss-Hermite quadrature
c...
c...     after the starting integrals have been computed,
c...     the recurrence is applied. First on a alone, then on b alone,
c...     and finally on both a and b
c...
c...
c...  Input parameters:
c...
c...  nmaxa            maximum n value for Gaussian a
c...  nmaxb            maximum n value for Gaussian b
c...  lp               angular momentum of the pseudopontential (PSP)
c...                   projector
c...  npseudo          array of powers of r of PSP
c...  ka               positional parameter
c...  kb               positional parameter
c...  aa               exponent of Gaussian a
c...  ab               exponent of Gaussian b
c...  alpha            array of exponents of PSP
c...  C                array of contraction coefficients of PSP
c...  len              length of the PSP expansion
c...  maxn             maximum n value, for dimensioning
c...  maxl             maximum l value, for dimensioning
c...  fspm1            multiplicative factor for type 2 integrals
c...  tol              tolerance threshold
c...
c...  Output parameter:
c...
c...  qq               table of scaled type 2 radial integrals
c...
      integer nmaxa,nmaxb,lp,len,maxn,maxl,npseudo(len)
      real*8 ka,kb,aa,ab,alpha(len),c(len),qq(0:maxn,0:maxl,0:maxl)
      real*8 tol
c...
c...  scratch areas for type 1 radial integrals and power series terms
c...
      integer maxnq1,maxlq1,maxst,nadd
      parameter (maxnq1=110,maxlq1=15,maxst=50,nadd=10)
      real*8 q1(0:maxnq1,0:maxlq1),cst(0:maxst)
c...
c...  local variables
c...
      integer nn,n,np,nps,na,nap,nb,nbp,nc
      integer la,lam,lap,lap2,lap1,lal,lb,lbm,lbp,lbp2,lbp1,lbl,l
      integer mlp2,i
      integer lq1,l12,nst,nim,nfm
      real*8 a,cc,pla,plb,kam1,kbm1,plakam1,plbkbm1,kab2,argd,rc2,sf
      real*8 za,zb,z1,z2,pref,sqram1,sum,rc,q2,qq1
c...
      real*8 xa(20),xb(20),bpa(20),bpb(20),pre(20)
c...
c...  roots and weights for Gauss-Hermite quadrature,
c...  computed according to "Numerical Recipes in Fortran, Second
c...  Edition", page 147-148, relative precision EPS=1.0e-14
c...
      real*8 x5(1:5)
      save x5
      data x5/2.020182870456086D0,9.585724646138185D-1,0.0D0,-9.58572464
     $  6138185D-1,-2.020182870456086D0/
      real*8 w5(1:5)
      save w5
      data w5/1.995324205904591D-2,3.936193231522412D-1,9.45308720482941
     $  9D-1,3.936193231522412D-1,1.995324205904591D-2/
c...
      real*8 x10(1:10)
      save x10
      data x10/3.436159118837738D0,2.53273167423279D0,1.756683649299882D
     $  0,1.036610829789514D0,3.429013272237046D-1,-3.429013272237046D-1
     $  ,-1.036610829789514D0,-1.756683649299882D0,-2.53273167423279D0,-
     $  3.436159118837738D0/
      real*8 w10(1:10)
      save w10
      data w10/7.640432855232615D-6,1.343645746781233D-3,3.3874394455481
     $  06D-2,2.401386110823135D-1,6.108626337353253D-1,6.10862633735325
     $  3D-1,2.401386110823135D-1,3.387439445548106D-2,1.343645746781233
     $  D-3,7.640432855232615D-6/
c...
      real*8 x20(1:20)
      save x20
      data x20/5.387480890011233D0,4.603682449550744D0,3.944764040115625
     $  D0,3.347854567383216D0,2.788806058428131D0,2.254974002089276D0,1
     $  .738537712116586D0,1.234076215395323D0,7.374737285453944D-1,2.45
     $  3407083009013D-1,-2.453407083009013D-1,-7.374737285453944D-1,-1.
     $  234076215395323D0,-1.738537712116586D0,-2.254974002089276D0,-2.7
     $  88806058428131D0,-3.347854567383216D0,-3.944764040115625D0,-4.60
     $  3682449550744D0,-5.387480890011233D0/
      real*8 w20(1:20)
      save w20
      data w20/2.229393645534152D-13,4.399340992273174D-10,1.08606937076
     $  9282D-7,7.802556478532063D-6,2.283386360163525D-4,3.243773342237
     $  863D-3,2.481052088746361D-2,1.090172060200217D-1,2.8667550536283
     $  41D-1,4.622436696006101D-1,4.622436696006101D-1,2.86675505362834
     $  1D-1,1.090172060200217D-1,2.481052088746361D-2,3.243773342237863
     $  D-3,2.283386360163525D-4,7.802556478532063D-6,1.086069370769282D
     $  -7,4.399340992273174D-10,2.229393645534152D-13/
c...
      real*8 ZERO,one,two,three,four
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,four=4.0D0)
c...
c...  4/sqrt(%PI)*4*%PI
c...
c...  the factor 4/sqrt(%PI) is for comparison with Argo.
c...
      real*8 fspm1
c     parameter (fspm1=28.35926161448825D0)
c...
      nn=nmaxa+nmaxb
      if(nn.gt.maxn)
     $  call nerror(1,'DomRad2','Maximum n value exceeded:',nn,maxn)
c
      if(max(nmaxa,nmaxb)+lp.gt.maxl)
     $  call nerror(2,'DomRad2','Maximum l value exceeded:',
     $             max(nmaxa,nmaxb)+lp,maxl)
c
      if(len.le.0)
     $  call nerror(3,'DomRad2','Improper value of parameter len',len,0)
c...
      mlp2=mod(lp,2)
      kab2=(ka+kb)*(ka+kb)
      argd=(ka*ka/(four*aa)+kb*kb/(four*ab))
c...
      do np=0,nn,1
        do la=0,min(np+lp,nmaxa+lp),1
          do lb=0,min(np+lp,nmaxb+lp),1
            QQ(np,la,lb)=ZERO
          enddo
        enddo
      enddo
c...
c...  CASE 1
c...
      if(ka.eq.zero.and.kb.eq.zero)then
        do nc=1,len
          rc2=kab2/(4.0d0*(aa+ab+alpha(nc)))
          cc=c(nc)*fspm1
          if(rc2-argd.gt.tol)then
            do n=2*mlp2,nn,2
              call DoQ12f(n+npseudo(nc),0,zero,aa+ab+alpha(nc),qq1)
              qq(n,0,0)=qq(n,0,0)+qq1*cc
            enddo
          endif
        enddo
c...
c...  CASE 2
c...
      else if(ka.ne.zero.and.kb.eq.zero)then
        kam1=one/ka
        sf=exp(-argd)*fspm1
        do na=0,nmaxa
          if(na.eq.0)then
            lal=mlp2
          else
            lal=na+lp
          endif
          do nb=mlp2,nmaxb,2
            n=na+nb
            do la=lal,na+lp,2
              do nc=1,len
                rc2=kab2/(4.0d0*(aa+ab+alpha(nc)))
                if(rc2-argd.gt.tol)then
                  call DoQ12f(n+npseudo(nc),la,ka,aa+ab+alpha(nc),qq1)
                  qq(n,la,0)=qq(n,la,0)+qq1*c(nc)*sf
                endif
              enddo
            enddo
          enddo
        enddo
        lam=0
        do nap=-(lp-2),nmaxa,2
          lam=lam+2
          lal=max(1-nap,0)
          pla=dfloat(2*lal+1)
          do lap=lal,nmaxa+lp-lam,1
            lap1=lap+1
            lap2=lap+2
            pla=pla+two
            plakam1=pla*kam1
            na=nap+lap
            do nb=mlp2,nmaxb,2
              n=na+nb
              qq(n,lap,0)=qq(n,lap2,0)+plakam1*qq(n-1,lap1,0)
            enddo
          enddo
        enddo
c...
c...  CASE 3
c...
      else if(ka.eq.zero.and.kb.ne.zero)then
        kbm1=one/kb
        sf=exp(-argd)*fspm1
        do nb=0,nmaxb
          if(nb.eq.0)then
            lbl=mlp2
          else
            lbl=nb+lp
          endif
          do na=mlp2,nmaxa,2
            n=na+nb
            do lb=lbl,nb+lp,2
              do nc=1,len
                rc2=kab2/(4.0d0*(aa+ab+alpha(nc)))
                if(rc2-argd.gt.tol)then
                  call DoQ12f(n+npseudo(nc),lb,kb,aa+ab+alpha(nc),qq1)
                  qq(n,0,lb)=qq(n,0,lb)+qq1*c(nc)*sf
                endif
              enddo
            enddo
          enddo
        enddo
        lbm=0
        do nbp=-(lp-2),nmaxb,2
          lbm=lbm+2
          lbl=max(1-nbp,0)
          plb=dfloat(2*lbl+1)
          do lbp=lbl,nmaxb+lp-lbm,1
            lbp1=lbp+1
            lbp2=lbp+2
            plb=plb+two
            plbkbm1=plb*kbm1
            nb=nbp+lbp
            do na=mlp2,nmaxa,2
              n=na+nb
              qq(n,0,lbp)=qq(n,0,lbp2)+plbkbm1*qq(n-1,0,lbp1)
            enddo
          enddo
        enddo
c...
c...  CASE 4
c...
c...  compute the starting integrals
c...
      else
        kam1=one/ka
        kbm1=one/kb
        do nc=1,len
          a=aa+ab+alpha(nc)
          sqram1=one/sqrt(a)
          rc=0.5d0*(ka+kb)*sqram1
          rc2=rc*rc
          if(rc2-argd.le.tol)goto 1
          nps=npseudo(nc)
          za=ka*sqram1
          zb=kb*sqram1
          if(max(nmaxa,nmaxb)+lp.gt.maxlq1)
     $      call nerror(4,'DomRad2','Maximum lq1 value exceeded',
     $                  max(nmaxa,nmaxb)+lp,maxlq1)
          if (rc2.lt.50.0d0)then
c...
c...  single power series
c...
            sf=exp(-argd)
            cc=c(nc)*sf*fspm1
            if(za.le.zb)then
              l12=0
              z1=za
              z2=zb
              lq1=nmaxb+lp
            else
              l12=1
              z1=zb
              z2=za
              lq1=nmaxa+lp
            endif
            nim=nn+2*lp
            nfm=min(nim+nadd*2,maxnq1)
            call recq1(lq1,nps,0,nfm,z2,maxnq1,q1)
            nst=nadd
            call csumterm(0,nst,z1,cst,maxst)
            do na=0,nmaxa
              if(na.eq.0)then
                lal=mlp2
              else
                lal=na+lp
              endif
              do nb=0,nmaxb
                if(nb.eq.0)then
                  lbl=mlp2
                else
                  lbl=nb+lp
                endif
                n=na+nb
                if(l12.eq.0)then
                  do la=lal,na+lp,2
                    pref=cc*ka**la*sqram1**(la+n+nps+1)
                    do lb=lbl,nb+lp,2
                      call sumspn(n,nps,la,lb,za,zb,nadd,nst,
     +                            maxst,cst,nfm,lq1,maxnq1,q1,sum)
                      qq(n,la,lb)=qq(n,la,lb)+pref*sum
                    enddo
                  enddo
                else
                  do lb=lbl,nb+lp,2
                    pref=cc*kb**lb*sqram1**(lb+n+nps+1)
                    do la=lal,na+lp,2
                      call sumspn(n,nps,lb,la,zb,za,nadd,nst,
     +                            maxst,cst,nfm,lq1,maxnq1,q1,sum)
                      qq(n,la,lb)=qq(n,la,lb)+pref*sum
                    enddo
                  enddo
                endif
              enddo
            enddo
          else if(rc2.le.500.0d0)then
c...
c...  twenty points Gauss-Hermite quadrature
c...
            cc=c(nc)*exp(rc2-argd)*sqram1*fspm1
            do i=1,20
              xa(i)=za*(x20(i)+rc)
              xb(i)=zb*(x20(i)+rc)
            enddo
            do na=0,nmaxa
              if(na.eq.0)then
                lal=mlp2
              else
                lal=na+lp
              endif
              do nb=0,nmaxb
                if(nb.eq.0)then
                  lbl=mlp2
                else
                  lbl=nb+lp
                endif
                n=na+nb
                do i=1,20
                  pre(i)=((x20(i)+rc)*sqram1)**(n+nps)*w20(i)
                enddo
                do la=lal,na+lp,2
                  call sbess1pn(la,xa,20,bpa)
                  do lb=lbl,nb+lp,2
                    call sbess1pn(lb,xb,20,bpb)
                    q2=zero
                    do i=1,20
                      q2=q2+pre(i)*bpa(i)*bpb(i)
                    enddo
                    qq(n,la,lb)=qq(n,la,lb)+q2*cc
                  enddo
                enddo
              enddo
            enddo
          else if(rc2.le.50000.0d0)then
c...
c...  ten points Gauss-Hermite quadrature
c...
            cc=c(nc)*exp(rc2-argd)*sqram1*fspm1
            do i=1,10
              xa(i)=za*(x10(i)+rc)
              xb(i)=zb*(x10(i)+rc)
            enddo
            do na=0,nmaxa
              if(na.eq.0)then
                lal=mlp2
              else
                lal=na+lp
              endif
              do nb=0,nmaxb
                if(nb.eq.0)then
                  lbl=mlp2
                else
                  lbl=nb+lp
                endif
                n=na+nb
                do i=1,10
                  pre(i)=((x10(i)+rc)*sqram1)**(n+nps)*w10(i)
                enddo
                do la=lal,na+lp,2
                  call sbess1pn(la,xa,10,bpa)
                  do lb=lbl,nb+lp,2
                    call sbess1pn(lb,xb,10,bpb)
                    q2=zero
                    do i=1,10
                      q2=q2+pre(i)*bpa(i)*bpb(i)
                    enddo
                    qq(n,la,lb)=qq(n,la,lb)+q2*cc
                  enddo
                enddo
              enddo
            enddo
          else
c...
c...  five points Gauss-Hermite quadrature
c...
            cc=c(nc)*exp(rc2-argd)*sqram1*fspm1
            do i=1,5
              xa(i)=za*(x5(i)+rc)
              xb(i)=zb*(x5(i)+rc)
            enddo
            do na=0,nmaxa
              if(na.eq.0)then
                lal=mlp2
              else
                lal=na+lp
              endif
              do nb=0,nmaxb
                if(nb.eq.0)then
                  lbl=mlp2
                else
                  lbl=nb+lp
                endif
                n=na+nb
                do i=1,5
                  pre(i)=((x5(i)+rc)*sqram1)**(n+nps)*w5(i)
                enddo
                do la=lal,na+lp,2
                  call sbess1pn(la,xa,5,bpa)
                  do lb=lbl,nb+lp,2
                    call sbess1pn(lb,xb,5,bpb)
                    q2=zero
                    do i=1,5
                      q2=q2+pre(i)*bpa(i)*bpb(i)
                    enddo
                    qq(n,la,lb)=qq(n,la,lb)+q2*cc
                  enddo
                enddo
              enddo
            enddo
          endif
 1        continue
        enddo
c...
c...  the starting integrals have been computed.
c...  Now fill the table by recurrence
c...
c...  recurrence on a
c...
        lam=0
        do nap=-(lp-2),nmaxa,2
          lam=lam+2
          lal=max(1-nap,0)
          pla=dfloat(2*lal+1)
          do lap=lal,nmaxa+lp-lam,1
            lap1=lap+1
            lap2=lap+2
            pla=pla+two
            plakam1=pla*kam1
            na=nap+lap
            do nb=0,nmaxb,1
              n=na+nb
              if(nb.eq.0)then
                lbl=mlp2
              else
                lbl=nb+lp
              endif
              do lb=lbl,nb+lp,2
                if(qq(n,lap,lb).eq.zero)then    ! to avoid repetitions
                  qq(n,lap,lb)=qq(n,lap2,lb)+plakam1*qq(n-1,lap1,lb)
                endif
              enddo
            enddo
          enddo
        enddo
c...
c...  recurrence on b
c...
        lbm=0
        do nbp=-(lp-2),nmaxb,2
          lbm=lbm+2
          lbl=max(1-nbp,0)
          plb=dfloat(2*lbl+1)
          do lbp=lbl,nmaxb+lp-lbm,1
            lbp1=lbp+1
            lbp2=lbp+2
            plb=plb+two
            plbkbm1=plb*kbm1
            nb=nbp+lbp
            do na=0,nmaxa,1
              n=na+nb
              if(na.eq.0)then
                lal=mlp2
              else
                lal=na+lp
              endif
              do la=lal,na+lp,2
                if(qq(n,la,lbp).eq.zero)then    ! to avoid repetitions
                  qq(n,la,lbp)=qq(n,la,lbp2)+plbkbm1*qq(n-1,la,lbp1)
                endif
              enddo
            enddo
          enddo
        enddo
c...
c...  recurrence on a and b
c...
        lam=0
        do nap=-(lp-2),nmaxa,2
          lam=lam+2
          lal=max(1-nap,0)
          pla=dfloat(2*lal+1)
          do lap=lal,nmaxa+lp-lam,1
            lap1=lap+1
            lap2=lap+2
            pla=pla+two
            plakam1=pla*kam1
            na=nap+lap
            lbm=0
            do nbp=-(lp-2),nmaxb,2
              lbm=lbm+2
              lbl=max(1-nbp,0)
              plb=dfloat(2*lbl+1)
              do lbp=lbl,nmaxb+lp-lbm,1
                lbp1=lbp+1
                lbp2=lbp+2
                plb=plb+two
                plbkbm1=plb*kbm1
                nb=nbp+lbp
                n=na+nb
                if(qq(n,lap,lbp).eq.zero)then    ! to avoid repetitions
                  qq(n,lap,lbp)=
     $                       qq(n,lap2,lbp2)+plbkbm1*qq(n-1,lap2,lbp1)+
     $            plakam1*(qq(n-1,lap1,lbp2)+plbkbm1*qq(n-2,lap1,lbp1))
                endif
              enddo
            enddo
          enddo
        enddo
      endif
c...
      end
