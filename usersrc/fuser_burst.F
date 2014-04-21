/***************************************************************/
/* COmPACT F77 user routine: USER_BURST(IB,QB)                 */
/*                                                             */
/* User F77 routine called everytime a burst '*bur' is loaded. */
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 RWM  9/7/97 */
/***************************************************************/

#include "compactF77.h"
#include "CompactGeomF77.h"
#include "cdbmapF77.h"
      FUNCTION USER_BURST(IB,QB)
      IMPLICIT NONE
      INTEGER USER_BURST
      INTEGER IB(2)           /* fake array sizes since C struct */
      REAL*4  QB(2)           /* fake array sizes since C struct */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer NBURST
      data NBURST/0/
      integer RUNCURRENT
      data RUNCURRENT/0/
      USER_BURST=0
      if(NBURST .eq. 0) call printcmpgeom(qgeom)
      if(NBURST .le. 5) call PRINTBURST(BURST_PTR)
      NBURST = NBURST + 1
C******************************************************************/
C* The following lines show how to call the printing routine      */
C* for the list of LKR dead cells, for the current run            */
C******************************************************************/
      IF (BUR_NRUN .NE. RUNCURRENT) THEN
         RUNCURRENT = BUR_NRUN
         CALL CMPLKRDEADPRINT()
         CALL CMPTIMEOFFSETPRINT(IB)
      ENDIF
C To print out some database data

c      PRINT *, ' =====>', RDB_NBUR, RDB_NBURBAD
c      PRINT *, ' =====>', BDB_NFILTCH, BDB_NGOODCH, BDB_NFILTNE, BDB_NGOODNE
C ----------- End of user F77 code -----------
      RETURN
      END


