/***************************************************************/
/* COmPACT F77 user routine: USER_SUPERBURST(IB,QB)            */
/*                                                             */
/* User F77 routine called everytime a SuperCOmPACT burst      */
/* '*sbur' is loaded.                                          */
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 BH 1/3/98   */
/***************************************************************/

#include "compactF77.h"
#include "cdbmapF77.h"

      FUNCTION USER_SUPERBURST(ISB,QSB)
      IMPLICIT NONE
      INTEGER USER_SUPERBURST
      INTEGER ISB(2)           /* fake array sizes since C struct */
      REAL*4  QSB(2)           /* fake array sizes since C struct */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      SBUR_BADB_SKIP=0 ! see fuser_superBurst.example.F to learn to use it

      USER_SUPERBURST=0
C ----------- End of user F77 code -----------
      RETURN
      END


