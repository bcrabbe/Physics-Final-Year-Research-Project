/********************************************************************/
/* COmPACT F77 user routine: USER_SUPERCMPEVENT(IB,QB,IE,QE)        */
/*                                                                  */
/* User F77 routine called everytime a SuperCOmPACT event is loaded.*/
/* A return value of greater than zero denotes that an error        */
/* condition has occured.                                           */
/*                                      BH 13/2/98      RWM 11/7/97 */
/********************************************************************/

#include "compactF77.h"

      FUNCTION USER_SUPERCMPEVENT(ISB,QSB,ISE,QSE)
      IMPLICIT NONE
      INTEGER USER_SUPERCMPEVENT
      INTEGER ISB(2),ISE(2)    /* fake array sizes since C structs */
      REAL*4  QSB(2),QSE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer nevt_read
      data nevt_read/0/

c      IF(nevt_read.LE.20) CALL PRINTSUPERCMPEVENT(SCMPEVT_PTR)
      nevt_read = nevt_read + 1
      USER_SUPERCMPEVENT=0



C ----------- End of user F77 code -----------
      RETURN
      END
