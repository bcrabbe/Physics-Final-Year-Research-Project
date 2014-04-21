/***************************************************************/
/* COmPACT F77 user routine: USER_SUPERCMPFILTER(IB,QB,IE,QE)  */
/*                                                             */
/* User F77 routine called everytime a COmPACT event is loaded */
/* and there is an output stream opened as well. A return      */
/* of less than zero denotes that the event should be written  */
/* to the output stream. A return of greater than zero denotes */
/* an error has occured.                                       */
/*                                  BH 10/2/98     RWM  9/7/97 */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_HYPERCMPFILTER(IHB,QHB,IHE,QHE)
      IMPLICIT NONE
      INTEGER USER_HYPERCMPFILTER
      INTEGER IHB(2),IHE(2)     /* fake array sizes since C structs */
      REAL*4  QHB(2),QHE(2)     /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      USER_HYPERCMPFILTER=-1
C ----------- End of user F77 code -----------
      RETURN
      END
