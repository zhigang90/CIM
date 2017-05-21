c -- This collection of cipher routines serves a similar purpose to
c -- <para_single.f>. The post-SCF energy code has so far been implemented
c -- using MPI only. These are dummy routines for the PVM version
c
c======================================================================
      subroutine fafOpenm(filename,ifileID)
      character*(*) filename
      end
c======================================================================
      subroutine fafCreatem(filename,ifileID)
      character*(*) filename
      end
c======================================================================
      subroutine fafCreates(filename,ifileID)
      character*(*) filename
      end
c======================================================================
      subroutine fafwrite(array,length,num,filename,n,ir)
      end
c======================================================================
      subroutine fafread(array,length,num,filename,n,ir)
      end
c======================================================================
      subroutine fafClosem(ifileID,keep,info)
      end
c======================================================================
      subroutine fafCloses(ifileID,keep,info)
      end
c======================================================================
      subroutine fafrequest
      end
c======================================================================
      subroutine fafnbread
      end
c======================================================================
      subroutine afwritebin(ndisk,lbin,ibin4,ibin1,indxbin,irecord,ij,
     1                      npairs,afname,islvid)
      end
c======================================================================
      subroutine afwritebingr(ndisk,lbin,ibin4,ibin1,indxbin,irecord,ij,
     1                      npairs,afname,islvid,isize)
      end
c======================================================================
      subroutine do_ccsd
      end
c======================================================================
      subroutine para_coulomb
      call coulomb
      end
c======================================================================
      subroutine par_diis_matrix()
      end
c======================================================================
      subroutine par_new_coeff()
      end
c======================================================================
      subroutine diis_slave_next()
      end
c======================================================================
      subroutine coupled_cluster()
      call coulomb
      end
c======================================================================
      subroutine array_files()
      end
