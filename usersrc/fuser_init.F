/*******************************************************/
/* COmPACT F77 user routine: USER_INIT                 */
/*                                                     */
/* User F77 routine called upon program startup to     */
/* allow initialization of user files, variables etc.  */
/*                                         RWM  9/7/97 */
/*******************************************************/

#include "compactF77.h"
      FUNCTION USER_INIT()
      IMPLICIT NONE
      INTEGER USER_INIT
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      USER_INIT=0
C ----------- End of user F77 code -----------
      RETURN
      END
