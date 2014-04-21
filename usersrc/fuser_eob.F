/***************************************************************/
/* COmPACT F77 user routine: USER_EOB(IEB,QEB)                 */
/*                                                             */
/* User F77 routine called at the end of every burst. It is    */
/* passed the end-of-burst structure for processing.           */
/*                                                 RWM  9/7/97 */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_EOB(IB,QB,IEB,QEB)
      IMPLICIT NONE
      INTEGER USER_EOB
      INTEGER IB(2),IEB(2)    /* fake array sizes since C struct */
      REAL*4  QB(2),QEB(2)    /* fake array sizes since C struct */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer NEOB
      data neob/0/
      USER_EOB=0
      if(neob .le. 5) call PRINTENDOFBURST(EOB_PTR)
      neob = neob + 1
C ----------- End of user F77 code -----------
      RETURN
      END

