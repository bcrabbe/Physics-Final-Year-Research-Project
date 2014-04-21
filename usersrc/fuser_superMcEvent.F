/********************************************************************/
/* COmPACT F77 user routine: USER_SUPERMCEVENT(IB,QB,IME,QME)       */
/*                                                                  */
/* User F77 routine called everytime a SuperMC event is loaded.     */
/* A return value of greater than zero denotes that an error        */
/* condition has occured.                                           */ 
/* Modifications: 09-01-98 IWS correct IE, IQ for new MC            */
/*                13-02-98 BH  modified for use with SuperCOmPACT   */
/*                                                 RWM 11/7/97      */
/********************************************************************/

#include "compactF77.h"
/* These two defines alter the superCmpEvent macros so they can access the */
/* superCmpEvent struct inside the superMcEvent struct.                    */

#define ISE(x) IMSE(SUPERMCEVENT_SUPERCMPEVENT+(x)-1)
#define QSE(x) QMSE(SUPERMCEVENT_SUPERCMPEVENT+(x)-1)

      FUNCTION USER_SUPERMCEVENT(ISB,QSB,IMSE,QMSE)
      IMPLICIT NONE
      INTEGER USER_SUPERMCEVENT
      INTEGER ISB(2),IMSE(2)    /* fake array sizes since C structs */
      REAL*4  QSB(2),QMSE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer nevt_read
      data nevt_read/0/
      IF(nevt_read.LE.20) CALL PRINTSUPERMCEVENT(SMCEVT_PTR) 
      nevt_read = nevt_read + 1
      USER_SUPERMCEVENT=0

      print *, ' KTYPE, DECAY', SMCEVENT_DECAY_KTYPE, SMCEVENT_DECAY_DTYPE
      print *, ' Evt LKRenergy', SEVT_NEVT, SEVT_LKRENERGY


C ----------- End of user F77 code -----------
      RETURN
      END
