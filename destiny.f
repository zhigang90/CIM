c===============================================================
c Dec. 97, KW :
c
c The general contraction loop goes ALWAYS over NGCD=ngci*ngcj*ngck*ngcl
c (the same for each quartet: it used to be  ngcq=lgenct(ijklp)  )
c.The diagonal cases are not distinguished anymore.
c Thus, the lgenct(nbls) vector does not have to be used anymore.
c It is still construced here because it is probably used in the DFT
c part of the program.
c
c Changes are made in the indexg procedure here.
c===============================================================
C March 26,96 KW : sub wriut has been commented out .
c==============================================================
      subroutine destret(ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *                   ipres,
     *                   ikbl,nbls,iis,jjs,ncs,inx,
     *                   icfg,jcfg,kcfg,lcfg,ngcd,indxr,
     *                   labels,length,lgenct )
c
c----------------------------------------------------------------
c This routine is called from calcint2 for each small blocks
c (if it has any intregrals)
c----------------------------------------------------------------
c This is for returning itegrals to the calling place .
c This routine constructs only label's arrays
c----------------------------------------------------------------
c
C???? implicit real*8 (a-h,o-z)
c
C???? common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
c
      dimension ijbl(nbl2,*)
      dimension ipres(*)
      dimension iis(*),jjs(*)
      dimension inx(12,*)
c
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension indxr(*)
c------------------------
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      length(1)=ilen
      length(2)=jlen
      length(3)=klen
      length(4)=llen
c----------------------------------------------------------------
ctest
c     write(8,*)' FROM DESTRET block no=',ikbl
c----------------------------------------------------------------
c
c  loop over quartets belonging to the block IKBL :
c
      ijkl =0
      do 100 ijp=nijbeg,nijend
         ijcs=ijbl(ibl,ijp)
         nklendx=nklend
         if(nklend.eq.0) nklendx=ijp
         do 100 klp=nklbeg,nklendx
            klcs=ijbl(kbl,klp)
            ijkl=ijkl+1
            if(ipres(ijkl).eq.0) go to 100
c
            ijklp=indxr(ijkl)
c                                   In the quartet IJKL :
c                                   --------------------
        ics=iis(ijcs)               ! contracted shell ICS
        jcs=jjs(ijcs)               !     -"-          JCS
        kcs=iis(klcs)               !     -"-          KCS
        lcs=jjs(klcs)               !     -"-          LCS
c
        if(ngcd.eq.1) then
           ngcq=1
           icfg(1)=inx(11,ics)
           jcfg(1)=inx(11,jcs)
           kcfg(1)=inx(11,kcs)
           lcfg(1)=inx(11,lcs)
        else
           call indexg(inx,ics,jcs,kcs,lcs,ijcs,klcs,
     *                 ilen,jlen,klen,llen, icfg,jcfg,kcfg,lcfg,ngcq)
        endif
c
          lgenct(ijklp)=ngcq
c
          do 150 iqu=1,ngcq
          icff=icfg(iqu)
          jcff=jcfg(iqu)
          kcff=kcfg(iqu)
          lcff=lcfg(iqu)
c
          labels(1,iqu,ijklp)=icff
          labels(2,iqu,ijklp)=jcff
          labels(3,iqu,ijklp)=kcff
          labels(4,iqu,ijklp)=lcff
c
  150   continue
  100 continue
c
      end
c==============================================================
      subroutine indexg(inx,ics,jcs,kcs,lcs,ijcs,klcs,
     *                  ilen,jlen,klen,llen,
     *                  icfg,jcfg,kcfg,lcfg,ngcq)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      dimension inx(12,*)
      dimension icfg(*),jcfg(*),kcfg(*),lcfg(*)
      dimension iix(100),jjx(100),kkx(100),llx(100)
c------------------------------------------------------------
c dim. 100 should be enough since the max. ge.con is 9, so 81
c is actually the max. for iix,jjx,kkx,llx
c------------------------------------------------------------
c
         icff=inx(11,ics)
         jcff=inx(11,jcs)
         kcff=inx(11,kcs)
         lcff=inx(11,lcs)
c
             ijpg=0
             icf=icff
             do 2041 igc=0,ngci
               jcf=jcff
               do 2042 jgc=0,ngcj
                  ijpg=ijpg+1
                  iix(ijpg)=icf
                  jjx(ijpg)=jcf
                  jcf=jcf+jlen
 2042          continue
               icf=icf+ilen
 2041        continue
c
             klpg=0
             kcf=kcff
             do 2043 kgc=0,ngck
               lcf=lcff
               do 2044 lgc=0,ngcl
                  klpg=klpg+1
                  kkx(klpg)=kcf
                  llx(klpg)=lcf
                  lcf=lcf+llen
 2044          continue
                kcf=kcf+klen
 2043        continue
c
             ijklg=0
             do 2045 ijp1=1,ijpg
             do 2045 klp1=1,klpg
                ijklg=ijklg+1
                icfg(ijklg)=iix(ijp1)
                jcfg(ijklg)=jjx(ijp1)
                kcfg(ijklg)=kkx(klp1)
                lcfg(ijklg)=llx(klp1)
 2045        continue
c
      ngcq=ijklg
c
      end
c==============================================================
ckw   subroutine wriut(isize,integ,last)
c-------------------------------------------------------------
c     this routine writes out a fully filled-up integral block
c     isize (inout): number of integrals
c     integ(about 3070 long) integrals and indices in nteger form
c     last=1 for the last block
c-------------------------------------------------------------
ckw   dimension integ(3070)
ckw   common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
ckw  1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
ckw  20),inpf,ioutf
ckw   np=npl(2)
ckw   integ(isize+1)=-1
ckw   if(last.eq.1) integ(isize+1)=-1 000 000 001
ckw   write(np) integ
ckw   isize=0
ckw   end
c==============================================================
      subroutine descend(i,j,k,l,iix)
      dimension iix(4)
c     this routine orders the 4 indices in canonical order
c    i>=j; k>=l; ij>=kl
      ij1=max0(i,j)
      ij0=min0(i,j)
      kl1=max0(k,l)
      kl0=min0(k,l)
      if(ij1.gt.kl1.or.(ij1.eq.kl1.and.ij0.ge.kl0)) then
        iix(1)=ij1
        iix(2)=ij0
        iix(3)=kl1
        iix(4)=kl0
      else
        iix(1)=kl1
        iix(2)=kl0
        iix(3)=ij1
        iix(4)=ij0
      end if
      end
c==============================================================
