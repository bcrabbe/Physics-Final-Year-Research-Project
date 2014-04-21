/*****************************************************************/
/* COmPACT F77 user routine: USER_HYPERCMPEVENT(IB,QB,IE,QE)     */
/*                                                               */
/* User F77 routine called everytime a hyperCmp event is loaded. */
/* A return value of greater than zero denotes that an error     */
/* condition has occured.                                        */
/*                         Marco.Clemencic@cern.ch  (11/7/2004)  */
/*****************************************************************/

#include "compactF77.h"
#include "CompactGeomF77.h"
#include "constants.h"

      FUNCTION USER_HYPERCMPEVENT(IHB,QHB,IHE,QHE)
      IMPLICIT NONE
      INTEGER USER_HYPERCMPEVENT
      INTEGER IHB(2),IHE(2)
      REAL*4  QHB(2),QHE(2)
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      
      USER_HYPERCMPEVENT=0
C ----------- End of user F77 code -----------
      RETURN
      END
