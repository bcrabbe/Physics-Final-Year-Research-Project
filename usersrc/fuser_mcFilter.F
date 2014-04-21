/***************************************************************/
/* COmPACT F77 user routine: USER_MCFILTER(IB,QB,IME,QME)      */
/*                                                             */
/* User F77 routine called everytime a MC event is loaded      */
/* and there is an output stream opened as well. A return      */
/* of less than zero denotes that the event should be written  */
/* to the output stream. A return of greater than zero denotes */
/* an error has occured.                                       */
/*                                                IWS 24/12/97 */
/* 10-02-2000 IWS add missing define (from Jean Duclos)        */
/***************************************************************/

#include "compactF77.h"
/* These two defines alter the cmpEvent macros so they can access the */
/* cmpEvent struct inside the mcEvent struct.                        */

#define IE(x) IME(MCEVENT_CMPEVT+(x)-1)
#define QE(x) QME(MCEVENT_CMPEVT+(x)-1)

      FUNCTION USER_MCFILTER(IB,QB,IME,QME)
      IMPLICIT NONE
      INTEGER USER_MCFILTER
      INTEGER IB(2),IME(2)    /* fake array sizes since C structs */
      REAL*4  QB(2),QME(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      USER_MCFILTER=0
C ----------- End of user F77 code -----------
      RETURN
      END
