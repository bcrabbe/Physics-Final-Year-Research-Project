/***************************************************************/
/* COmPACT F77 user routine: USER_SUPERMCFILTER(IB,QB,IME,QME) */
/*                                                             */
/* User F77 routine called everytime a MC event is loaded      */
/* and there is an output stream opened as well. A return      */
/* of less than zero denotes that the event should be written  */
/* to the output stream. A return of greater than zero denotes */
/* an error has occured.                                       */
/*                                  BH 10/2/98    IWS 24/12/97 */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_SUPERMCFILTER(IB,QB,IME,QME,IHE,QHE)
      IMPLICIT NONE
      INTEGER USER_SUPERMCFILTER
      INTEGER IB(2),IME(2),IHE(2)    /* fake array sizes since C structs */
      REAL*4  QB(2),QME(2),QHE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      USER_SUPERMCFILTER=0
C ----------- End of user F77 code -----------
      RETURN
      END
