/***************************************************************/
/* COmPACT F77 user routine: USER_KMU3EVENT(IB,QB,IE,QE)       */
/*                                                             */
/* User F77 routine called everytime a Kmu3 event is loaded.   */
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 RWM  9/7/97 */
/* Modified: 19-04-98 IWS Correct defines (IUE instead of IKE) */
/* ========                                                    */
/***************************************************************/

#include "compactF77.h"
/* These two defines alter the cmpEvent macros so they can access the */
/* cmpEvent struct inside the kmu3Event struct.                        */

#define IE(x) IUE(KMU3EVENT_CMPEVT+(x)-1)
#define QE(x) QUE(KMU3EVENT_CMPEVT+(x)-1)

      FUNCTION USER_KMU3EVENT(IB,QB,IUE,QUE)
      IMPLICIT NONE
      INTEGER USER_KMU3EVENT
      INTEGER IB(2),IUE(2)    ! fake array sizes since C structs 
      REAL*4  QB(2),QUE(2)    ! fake array sizes since C structs
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      INTEGER NEVENTREAD
      DATA NEVENTREAD/0/
      IF(NEVENTREAD.LE.20) CALL PRINTKMU3EVENT(IE(1),QE(1))
      NEVENTREAD = NEVENTREAD + 1
      USER_KMU3EVENT=0
C ----------- End of user F77 code -----------
      RETURN
      END
