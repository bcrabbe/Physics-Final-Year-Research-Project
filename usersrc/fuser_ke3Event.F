/***************************************************************/
/* COmPACT F77 user routine: USER_KE3EVENT(IB,QB,IE,QE)        */
/*                                                             */
/* User F77 routine called everytime a Ke3 event is loaded.    */
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 RWM  9/7/97 */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_KE3EVENT(IB,QB,IKE,QKE)
      IMPLICIT NONE
      INTEGER USER_KE3EVENT
      INTEGER IB(2),IKE(2)    /* fake array sizes since C structs */
      REAL*4  QB(2),QKE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      IF(KNEVENT.LE.20) CALL PRINTKE3EVENT(KE3EVT_PTR)
C      CALL PRINTKE3EVENT(KE3EVT_PTR)

      USER_KE3EVENT=0
C ----------- End of user F77 code -----------
      RETURN
      END
