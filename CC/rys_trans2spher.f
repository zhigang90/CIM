c===================================================================
      subroutine trans2spher(ics,jcs,kcs,lcs,inx,
     1                       ispher,buf,buf2)
c  Transforms a shell integral to spherical harmonics
c
c  Arguments
c  INTENT(IN)
c  ics,jcs,kcs,lcs  = 4 contracted shell indices
c  inx              = basis set contraction info
c  ispher           = array showing the presence of spherical harmonics functions
c                     This should be precomputed for the contracted shells
c                     The value is 0 if the function is Cartesian, 2 for D (D5),
c                     3 for F (F7), 4 for G (G9), 5 for H (H11)
c  INTENT(INOUT)
c  buf              = Cartesian integrals, replaced by spherical harmonics
c  buf2             = temporary storage for the transformed integrals
      implicit real*8 (a-h,o-z)
      dimension inx(12,*),ispher(*),buf(*),buf2(*),lcart(13)
      data lcart/1,3,4,6,6,10,10,15,21,28,15,21,28/
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
      ktyp=inx(12,kcs)
      ltyp=inx(12,lcs)
      ilen=lcart(ityp)
      jlen=lcart(jtyp)
      klen=lcart(ktyp)
      llen=lcart(ltyp)
      if(ispher(ics).gt.0) then
        jkln=jlen*klen*llen
        if(ispher(ics).eq.2) then
          call dtrans1(1,buf,buf2,jkln)
          ilen=5
        else if(ispher(ics).eq.3) then
c   Note that the call of the ftrans routines is different from the rest
          call ftrans1(1,buf,buf2,jlen,klen,llen)
          ilen=7
        else if(ispher(ics).eq.4) then
          call gtrans1(1,buf,buf2,jkln)
          ilen=9
        else if(ispher(ics).eq.5) then
          call htrans1(1,buf,buf2,jkln)
          ilen=11
        end if
      end if                  
c
      if(ispher(jcs).gt.0) then
        kln=klen*llen
        if(ispher(jcs).eq.2) then
          call dtrans2(1,buf,buf2,ilen,kln)
          jlen=5
        else if(ispher(jcs).eq.3) then
c   Note that the call of the ftrans routines is different from the rest
          call ftrans2(1,buf,buf2,ilen,klen,llen)
          jlen=7
        else if(ispher(jcs).eq.4) then
          call gtrans2(1,buf,buf2,ilen,kln)
          jlen=9
        else if(ispher(jcs).eq.5) then
          call htrans2(1,buf,buf2,ilen,kln)
          jlen=11
        end if
      end if                  
c
      if(ispher(kcs).gt.0) then
        ijlen=ilen*jlen
        if(ispher(kcs).eq.2) then
          call dtrans3(1,buf,buf2,ijlen,llen)
          klen=5
        else if(ispher(kcs).eq.3) then
c   Note that the call of the ftrans routines is different from the rest
          call ftrans3(1,buf,buf2,ilen,jlen,llen)
          klen=7
        else if(ispher(kcs).eq.4) then
          call gtrans3(1,buf,buf2,ijlen,llen)
          klen=9
        else if(ispher(kcs).eq.5) then
          call htrans3(1,buf,buf2,ijlen,llen)
          klen=11
        end if
      end if                  
c
      if(ispher(lcs).gt.0) then
        ijklen=ilen*jlen*klen
        if(ispher(lcs).eq.2) then
          call dtrans4(1,buf,buf2,ijklen)
        else if(ispher(lcs).eq.3) then
c   Note that the call of the ftrans routines is different from the rest
          call ftrans4(1,buf,buf2,ilen,jlen,klen)
        else if(ispher(lcs).eq.4) then
          call gtrans4(1,buf,buf2,ijklen)
        else if(ispher(lcs).eq.5) then
          call htrans4(1,buf,buf2,ijklen)
          klen=11
        end if
      end if                  
c
      end
c===================================================================
