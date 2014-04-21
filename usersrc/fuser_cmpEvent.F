/***************************************************************/
/* COmPACT F77 user routine: USER_CMPEVENT(IB,QB,IE,QE)        */
/*                                                             */
/* User F77 routine called everytime a COmPACT event is loaded.*/
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 RWM 11/7/97 */
/***************************************************************/

#include "compactF77.h"
#include "CompactGeomF77.h"
#include "constants.h"

      FUNCTION USER_CMPEVENT(IB,QB,IE,QE)
      IMPLICIT NONE
      INTEGER USER_CMPEVENT
      INTEGER IB(2),IE(2)    /* fake array sizes since C structs */
      REAL*4  QB(2),QE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      integer nevt_read
      integer icorr
      data nevt_read/0/
      real cogE, pt2E, phiE, EangleE, mlambdaE(2), DtargetE
      integer ktype, iv
      integer i
      real x,y,z

c      IF(nevt_read.LE.20) then
c         print *,'-----------------------------------------------------|'
c         print *,'| Selection flags: Charged/Neutral',ACHARG_IFLAG, ANEUT_IFLAG,
c     +            '               |'
c        call PRINTFLAGCORR(EVT_FLAGCORR)
c        CALL PRINTCMPEVENT(CMPEVT_PTR)
c      endif
      
      nevt_read = nevt_read + 1
      USER_CMPEVENT=0

C ----------- End of user F77 code -----------
      RETURN
      END

