#include "maxi.h"
module ccounters
 integer :: ic_CCAOnce=0,ic_Fileinit41a=0,ic_CCYZ_file=0,iRstartold=0
 integer :: ic_CCYZ_read=0,ic_CCAOmas=0,xxicountere=0,ic_coefmast=0
 integer :: ic_initsingles=0,ic_initsingles_master=0,ic_initsingles_slave=0
 integer :: ic_CoefInit_slave=0,ic_CoefInit=0,ic_ReaderCoefLin=0
 integer :: ic_prepare_Qparts_file=0,ic_NEW_EEO_INT=0,ic_dump_data=0
 integer :: ic_MAS_EEO_INT=0
! DIIS
 integer :: istoreno=0,icalculate=0,ijold=0
 integer :: ishift_record=0,indexer=0,idiis=0,idiis_last=0
 integer :: ifirst_rec=0,isiglresstore=0,isiglcoestore=0,max_rec=0
 integer, parameter :: abs_max_rec=(DIIS_LENGTH)
 integer :: irecord_tab(abs_max_rec)
 logical :: do_diis=.true.,diis_restart=.false.
!END DIIS
!GenCoulExInt and master_Gen
 integer :: ieeocounter=0,txcounter=0,tccounter=0,keepmax=30
 integer :: ttcounter=0,icounterx=0,icounterc=0,icountere=0,icountertx=0
 integer :: icountertc=0, icountertt=0
!END GenCoulExInt
!     acc_and_write_K
 integer :: curr_rec = 0, curr_pos = 0, max_rec_K = 0
! END acc_and_write_K
end module ccounters
subroutine reset_counters
use ccounters
implicit none
integer i
 ic_CCAOnce=0
 ic_Fileinit41a=0
 ic_CCYZ_file=0
 iRstartold=0
 ic_CCYZ_read=0
 ic_CCAOmas=0
 xxicountere=0
 ic_coefmast=0
 ic_initsingles=0
 ic_initsingles_master=0
 ic_initsingles_slave=0
 ic_CoefInit_slave=0
 ic_CoefInit=0
 ic_ReaderCoefLin=0
 ic_prepare_Qparts_file=0
 ic_NEW_EEO_INT=0
 ic_dump_data=0
 ic_MAS_EEO_INT=0
!
 istoreno=0
 icalculate=0
 ijold=0
 ishift_record=0
 indexer=0
 idiis=0
 idiis_last=0
 ifirst_rec=0
 isiglresstore=0
 isiglcoestore=0
 max_rec=0
 do i=1,abs_max_rec
   irecord_tab(i)=0
 enddo
 do_diis=.true.
 diis_restart=.false.
!END DIIS
!GenCoulExInt and master_Gen
 ieeocounter=0
 txcounter=0
 tccounter=0
 keepmax=30
 ttcounter=0
 icounterx=0
 icounterc=0
 icountere=0
 icountertx=0
 icountertc=0
 icountertt=0
!END GenCoulExInt
!     acc_and_write_K
 curr_rec = 0
 curr_pos = 0
 max_rec_K = 0
end subroutine reset_counters
