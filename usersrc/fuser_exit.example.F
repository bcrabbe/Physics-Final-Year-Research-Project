/*******************************************************/
/* COmPACT F77 user routine: USER_EXIT                 */
/*                                                     */
/* User F77 routine called upon program finishing to   */
/* allow closing of user files etc.                    */
/*                                         RWM  9/7/97 */
/*******************************************************/

#include "compactF77.h"

      FUNCTION USER_EXIT()
      IMPLICIT NONE
      INTEGER USER_EXIT
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      INTEGER LREC,ISTAT,ICYCLE
      LREC=1024                 ! record length
      CALL HROPEN(1,'OUTPUT','compact.hbook',
     +     'N',LREC,ISTAT)      !open HBOOK file for output
      CALL HROUT(0,ICYCLE,' ')  ! write out all histograms
      CALL HREND('OUTPUT')      ! close HBOOK file
      USER_EXIT=0
C ----------- End of user F77 code -----------
      RETURN
      END
