/***************************************************************/
/* COmPACT F77 user routine: USER_MCEVENT(IB,QB,IME,QME)       */
/*                                                             */
/* User F77 routine called everytime a COmPACT event is loaded.*/
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/* Modifications: 09-01-98 IWS correct IE, IQ for new MC       */
/*                                                 RWM 11/7/97 */
/***************************************************************/

#include "compactF77.h"
/* These two defines alter the cmpEvent macros so they can access the */
/* cmpEvent struct inside the mcEvent struct.                        */

#define IE(x) IME(MCEVENT_CMPEVT+(x)-1)
#define QE(x) QME(MCEVENT_CMPEVT+(x)-1)

      FUNCTION USER_MCEVENT(IB,QB,IME,QME)
      IMPLICIT NONE
      INTEGER USER_MCEVENT
      INTEGER IB(2),IME(2)    /* fake array sizes since C structs */
      REAL*4  QB(2),QME(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer nevt_read
      data nevt_read/0/
      IF(nevt_read.LE.20) CALL PRINTMCEVENT(MCEVT_PTR)
      nevt_read = nevt_read + 1
      USER_MCEVENT=0


C ----------- End of user F77 code -----------
      RETURN
      END
