/***************************************************************/
/* COmPACT F77 user routine: USER_SUPEREOB(IEB,QEB)            */
/*                                                             */
/* User F77 routine called at the end of every SuperCOmPACT    */
/* burst. It is passed the Super end-of-burst structure for    */
/* processing.                                                 */
/*                                                 BH 1/3/98   */
/***************************************************************/

#include "compactF77.h"

      FUNCTION USER_SUPEREOB(ISB,QSB,ISEB,QSEB)
      IMPLICIT NONE
      INTEGER USER_SUPEREOB
      INTEGER ISB(2),ISEB(2)    /* fake array sizes since C struct */
      REAL*4  QSB(2),QSEB(2)    /* fake array sizes since C struct */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer nsupereob
      data nsupereob/0/
      USER_SUPEREOB=0
C      if(nsupereob .le. 5) call PRINTSUPERENDOFBURST(SEOB_PTR)
C      nsupereob = nsupereob + 1
C ----------- End of user F77 code -----------
      RETURN
      END

