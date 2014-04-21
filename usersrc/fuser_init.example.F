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
      INTEGER NWPAWC,I
      PARAMETER (NWPAWC=2000000)
      COMMON /PAWC/ I(NWPAWC)
      CALL HLIMIT(NWPAWC)
      CALL HBOOK1(100,'Particle Cuts',21,-0.5,20.5,0.0)
      CALL HBOOK1(101,'HAC Energy',50,0.0,150.0,0.0)
      CALL HBOOK1(102,'LKR Cluster Energies',50,0.0,150.0,0.0)
      CALL HBOOK1(103,'E/P for LKR',50,0.0,1.0,0.0)
      CALL HBOOK1(104,'E/P for 1 LKR',50,0.0,1.0,0.0)
      CALL HBOOK1(110,'Closest LKR Cluster',50,0.0,50.0,0.0)
      CALL HBOOK1(900,'No LKR clusters',6,-0.5,5.5,0.0)
      USER_INIT=0
C ----------- End of user F77 code -----------
      RETURN
      END
