/***************************************************************/
/* COmPACT F77 user routine: USER_KMU3FILTER(IB,QB,IKE,QKE)    */
/*                                                             */
/* User F77 routine called everytime a Kmu3 event is loaded    */
/* and there is an output stream opened as well. A return      */
/* of less than zero denotes that the event should be written  */
/* to the output stream. A return of greater than zero denotes */
/* an error has occured.                                       */
/*                                                 RWM  9/7/97 */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_KMU3FILTER(IB,QB,IKE,QKE)
      IMPLICIT NONE
      INTEGER USER_KMU3FILTER
      INTEGER IB(2),IKE(2)    /* fake array sizes since C structs */
      REAL*4  QB(2),QKE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      USER_KMU3FILTER=0
C ----------- End of user F77 code -----------
      RETURN
      END
