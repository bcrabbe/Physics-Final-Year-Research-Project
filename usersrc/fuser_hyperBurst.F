/*****************************************************************/
/* COmPACT F77 user routine: USER_HYPERBURST(IB,QB,IE,QE)        */
/*                                                               */
/* User F77 routine called everytime a hyperBurst is loaded.     */
/* A return value of greater than zero denotes that an error     */
/* condition has occured.                                        */
/*                         Marco.Clemencic@cern.ch  (11/7/2004)  */
/*****************************************************************/

#include "compactF77.h"
#include "CompactGeomF77.h"
#include "constants.h"

      FUNCTION USER_HYPERBURST(IHB,QHB)
      IMPLICIT NONE
      INTEGER USER_HYPERBURST
      INTEGER IHB(2)
      REAL*4  QHB(2)
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
C     FORTRAN interface not complete!!!!
CCC      HBUR_BADB_SKIP=0 ! see fuser_superBurst.example.F to learn to use it
      IHB(HYPERBURST_BADB + BADBURST_SKIP) = 0
      
      USER_HYPERBURST=0
C ----------- End of user F77 code -----------
      RETURN
      END
