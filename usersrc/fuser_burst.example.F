/***************************************************************/
/* COmPACT F77 user routine: USER_BURST(IB,QB)                 */
/*                                                             */
/* User F77 routine called everytime a burst '*bur' is loaded. */
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 RWM  9/7/97 */
/* Modified: 15-02-99 IWS change default in example            */
/* ========                                                    */
/***************************************************************/

#include "compactF77.h"
#include "CompactGeomF77.h"

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
      PRINT *,' XOFF DCH',
     + BUR_XOFFDCH_NTRANSITION, BUR_XOFFDCH_FRACTIMEON, BUR_XOFFDCH_FIRSTTIMEON
C*************************************************************/ 
C* The following lines show how to set up the selection      */
C* of analysis routines.                                     */
C* The default is set in compact to call SOME analysis       */
C* routines; if the user wants to change the default, here is*/
C* the way to go.                                            */
C* The code which calls the analysis routines is in          */
C* /afs/cern.ch/na48/offline2/compact/compact-4.3/compact/   */
C*       /rlib/src/cmpAnalysis.c                             */
C* The routine USER_TAGTIME is called by default             */
C* only for good neutral or charged events                   */
C* The following example is equivalent to setting            */
C* BUR_CALLANAROUTINE=1 i.e. letting compact decide which    */
C* routines to call.                                         */
C* If the user sets BUR_CALLANAROUTINE=0 NO routine          */
C* will be called from compact.                              */
C* Some routines are called ANYWAY:                          */
C* USER_BADBURST                                             */
C* USER_ESPY, USER_STD_EP,USER_STD_RCOG,USER_STD_PTSQ        */
C*************************************************************/
C*the user wants to choose the routines
      BUR_CALLANAROUTINE=-1    
C*to call selcharged
      BUR_TOCALL_SELCHARGED=0  
C*to call sel2pi0  
      BUR_TOCALL_SEL2PI0=1 
C*to call sel3pi0 - only if no 2pi0 cand.     
      BUR_TOCALL_SEL3PI0=1   
C*to call sel2gam
      BUR_TOCALL_SEL2GAM=1 
C*NOT to call correction for Bfield in blue tube (1997)     
      BUR_TOCALL_BLUEFIELD=0   
C*NOT to call lkr pedestal correction (1997)   
      BUR_TOCALL_LKRPEDCOR=0 
C*NOT to call lkr osition correction    
      BUR_TOCALL_LKRPOSCOR=0 
C*NOT to call lkr energy sharing correction (1997)   
      BUR_TOCALL_LKRSHARING=0 
C*to call hodotime    
      BUR_TOCALL_HODOTIME=1 
C*to call nhodtime       
      BUR_TOCALL_NHODTIME=1 
C*to call lkrtime only if there is a 2pi0/3pi0 cand.      
      BUR_TOCALL_LKRTIME=1      
C*NOT to call tagtime  
      BUR_TOCALL_TAGTIME=0  
C*to call flagging of good AKS events     
      BUR_TOCALL_AKSFLAG=1      
C*"new" muon reconstruction (1997) - NOT IN COMPACT ANYMORE (16-03-99)
      BUR_TOCALL_MUON_REC=0    
C*NOT to call muon rejection (1997) - NOT USEFULL ANYMORE     
      BUR_TOCALL_MUON_REJECT=0      
C*to call geometry correction routine (1997 and 1998 data)
      BUR_TOCALL_GEOMCOR=1
C*to call magnet correction routine (adjust the Kaon mass to be flat with time)
      BUR_TOCALL_MAGNETCOR=1
C*to call NEW charged event selection routine
      BUR_TOCALL_NEWCHARGED=1
C*to call track muon reconstruction - THE RECOMMANDED ONE
      BUR_TOCALL_MUON_TRACKREC=1
C*NOT to call vertex muon reconstruction (1997)     
      BUR_TOCALL_VTX_TRACKREC=0     
C*NOT to call lkr calibration correction 
      BUR_TOCALL_LKRCALCOR=0 
C*NOT to call lkr corrections (change of temperature) for KS HI 2000 data
      BUR_TOCALL_LKRCALHI2K=0 

C******************************************************************/
C* The following lines show how to call the routine user_badburst */
C* which fills the BUR_BADB_xxx variables to be used in the       */
C* user routines                                                  */
C******************************************************************/
      call printbadburst(IB)
C*************************************************************/ 
C* The following lines show how to reject events according   */
C* to the "bad burts" flags - The routine USER_BADBURST      */
C* SHOULD have been called before.                           */
C* From version 4.0 of compact (march 1998) the user_badburst*/
C* routine is called by compact.                             */
C* From version 4.1 of compact (may 1998) if the flag        */
C* BUR_BADB_SKIP > 0: the entire burst is skipped by compact */
C* This should NOT be used for filtered files where there are*/
C* more than ONE burst.                                      */
C* It is up to the user to decide on which criteria to reject*/
C* a burst; the example given here uses all the available    */
C* flags to reject a burst and also require the burst        */
C* to have epsPrime/eps data.                                */
C*************************************************************/ 
C*
C* Reject events for which relevant detectors are NOT OK
C* For epsilonPrime/epsilon analysis 
C*
      if( (BUR_BADB_TAG.lt.0) .OR. 
     +(BUR_BADB_AKS.lt.0) .OR.
     +(BUR_BADB_DCH.lt.0) .OR.
     +(BUR_BADB_MBX.lt.0) .OR.
     +(BUR_BADB_HODC.lt.0) .OR.
     +(BUR_BADB_LKR.lt.0) .OR.
     +(BUR_BADB_NUT.lt.0) .OR.
     +(BUR_BADB_MUV.lt.0) .OR.
     +(BUR_BADB_PMB.lt.0) .OR.
     +(BUR_BADB_CLK.lt.0) .OR.
     +(BUR_BADB_SPARE2.lt.0) .OR.
     +(BUR_BADB_NOEPS.lt.0) ) BUR_BADB_SKIP=1

      if(BUR_BADB_SKIP.gt.0) RETURN


C******************************************************************/
C* The following lines show how to call the printing routine      */
C* for the list of LKR dead cells, for the current run            */
C******************************************************************/
      IF (BUR_NRUN .NE. RUNCURRENT) THEN
         RUNCURRENT = BUR_NRUN
         CALL CMPLKRDEADPRINT()
         CALL CMPTIMEOFFSETPRINT(IB)
      ENDIF
C******************************************************************/
C The following lines show how to call the routine to empty part  */
C of the burst structure in order to save space when filtering    */
C events                                                          */
C******************************************************************/
      call EMPTY_BURST(IB,100)
C ----------- End of user F77 code -----------
      RETURN
      END


