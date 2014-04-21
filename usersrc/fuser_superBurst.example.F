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
      integer NSUPERBURST
      data NSUPERBURST/0/
      USER_SUPERBURST=0
C      if(NSUPERBURST .le. 5) call PRINTSUPERBURST(SBURST_PTR)
      NSUPERBURST = NSUPERBURST + 1
C*
C* Reject events for which relevant detectors are NOT OK
C* For epsilonPrime/epsilon analysis 
C*
      if( (SBUR_BADB_TAG.ne.0) .OR. 
     +(SBUR_BADB_AKS.ne.0) .OR.
     +(SBUR_BADB_DCH.ne.0) .OR.
     +(SBUR_BADB_MBX.ne.0) .OR.
     +(SBUR_BADB_HODC.ne.0) .OR.
     +(SBUR_BADB_LKR.ne.0) .OR.
     +(SBUR_BADB_NUT.ne.0) .OR.
     +(SBUR_BADB_MUV.ne.0) .OR.
     +(SBUR_BADB_PMB.ne.0) .OR.
     +(SBUR_BADB_CLK.ne.0) .OR.
     +(SBUR_BADB_SPARE2.ne.0) .OR.
     +(SBUR_BADB_NOEPS.ne.0) ) SBUR_BADB_SKIP=1

      if(SBUR_BADB_SKIP.gt.0) RETURN
C To print out some database data

c      PRINT *, ' =====>', RDB_NBUR, RDB_NBURBAD
c      PRINT *, ' =====>', BDB_NFILTCH, BDB_NGOODCH, BDB_NFILTNE, BDB_NGOODNE

C ----------- End of user F77 code -----------
      RETURN
      END


