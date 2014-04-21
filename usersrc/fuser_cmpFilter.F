/***************************************************************/
/* COmPACT F77 user routine: USER_CMPFILTER(IB,QB,IE,QE)       */
/*                                                             */
/* User F77 routine called everytime a COmPACT event is loaded */
/* and there is an output stream opened as well. A return      */
/* of less than zero denotes that the event should be written  */
/* to the output stream. A return of greater than zero denotes */
/* an error has occured.                                       */
/*                                                 RWM  9/7/97 */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_CMPFILTER(IB,QB,IE,QE)
      IMPLICIT NONE
      INTEGER USER_CMPFILTER
      INTEGER IB(2),IE(2)     /* fake array sizes since C structs */
      REAL*4  QB(2),QE(2)     /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      USER_CMPFILTER=0
C ----------- End of user F77 code -----------
      RETURN
      END
