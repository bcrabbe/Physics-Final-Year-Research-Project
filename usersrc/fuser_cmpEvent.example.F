/***************************************************************/
/* COmPACT F77 user routine: USER_CMPEVENT(IB,QB,IE,QE)        */
/*                                                             */
/* User F77 routine called everytime a COmPACT event is loaded.*/
/* A return value of greater than zero denotes that an error   */
/* condition has occured.                                      */
/*                                                 RWM  9/7/97 */
/* Modifications: 16-02-99 IWS correct ACHOD_TRACE(1),(2)      */
/* =============                                               */
/***************************************************************/

#include "compactF77.h"
#include "CompactGeomF77.h"  /* To get detector geometry */
/****************************************/
/* Definitions for cuts used in ntpread */
/****************************************/

#define PMIN       5.0      /* minimum track momentum */
#define CDA_CUT    1.4      /* CDA maximum in cm */
#define ZMIN       0.0      /* minimum vertex position (KStarget = 0cm) */
#define ZMAX    9000.0      /* maximum vertex position */
#define RMIN_LKR  17.0      /* minimum radius of LKR */
#define RMAX_LKR 110.0      /* maximum radius of LKR */
#define P_LAMBDA   3.0      /* Track momenta ratio cut to remove lambdas */
#define MAX_TRACKS 5        /* maximum number of tracks */
#define ASSOCLKR  10.0      /* max distance for LKR cluster association */
#define ASSOCHAC  15.0      /* max distance for HAC cluster association */


      FUNCTION USER_CMPEVENT(IB,QB,IE,QE)
      IMPLICIT NONE
      INTEGER USER_CMPEVENT
      INTEGER IB(2),IE(2)     /* fake array sizes since C structs */
      REAL*4  QB(2),QE(2)     /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C WARNING: do not alter things before this line
C ---------- Add user F77 code here ----------
      REAL RCLUS,XPRED,YPRED,DR,DRMIN
      INTEGER I,J,ICLUS
      INTEGER ILKR(MAX_TRACKS)       ! array of associated LKR clusters
      REAL EOP                      ! temporary E/p variable
      REAL NOCLUS
      integer iok, ICORR
      INTEGER ITP, ITN
      INTEGER IV
      INTEGER USER_HODOTIME, USER_TAGTIME
      EXTERNAL USER_HODOTIME, USER_TAGTIME
C************************************************************/ 
C* Please read the comments concerning the analysis routines*/
C* To call these routines, please add the following line    */
C* after the variable definitions                           */
C* #include "CommonGeomF77.h" /* Get the detector geomtry   */
C************************************************************/ 
#include "CommonGeomF77.h"    /* Detector Geometry */

      USER_CMPEVENT=0
C*************************************************************/ 
C* The following lines show how to reject events according   */
C* to the "bad burts" flags - The routine USER_BADBURST      */
C* SHOULD have been called from USER_BURST as shown in       */
C* fuser_burst.example.F file                                */
C* All the routines are not yet effective                    */
C* Please note that from compact-4.1 there is the possibility*/
C* to skip the ENTIRE burst (see in user_burst.example       */
C* ONLY valid when the input file contains only ONE burst    */
C*************************************************************/ 
C
C Reject events for which NUT or MBX or TAG or LKR or MUV not OK
C

      if(  BUR_BADB_NUT .lt. 0 .OR.
     +     BUR_BADB_MBX .lt. 0 .OR.
     +     BUR_BADB_LKR .lt. 0 .OR.
     +     BUR_BADB_MUV .lt. 0 .OR.
     +     BUR_BADB_TAG .lt. 0) RETURN
C 
C     For epsilonP/epsilon analysis 
C     if(  BUR_BADB_NOEPS .lt. 0) RETURN
C
      print *,'-----------------------------------------------------|'
      print *,'| Selection flags: Charged/Neutral',ACHARG_IFLAG, ANEUT_IFLAG,
     +     '               |'
      call PRINTFLAGCORR(EVT_FLAGCORR)
C************************************************************/ 
C* The following lines call the NEW MUV reconstruction      */
C* routine; if it is called, the content of the             */
C* MUV structure will be overwritten by the result of this  */
C* routine.                                                 */
C* Therefore SELCHARGED has to be called again              */
C************************************************************/ 

      if(BUR_NRUN.lt.6000) then
         call USER_NEW_MUON(IB, IE)
      endif
C 
C  Compare energies before and after LKR ped cor
C
      do i=1, EVT_NLKR
         print *, ' ================> ecorrke3(before) ',LKR_ECORRKE3(i),
     +     '--> (after)',LKR_ANAVAR(i,1)

      enddo
C 
C   !!!!!    !!!!!     !!!!!     !!!!!     !!!!!     !!!!! 
C   The following lines just show how to call the analysis 
C   routines; they dont need to be called (as they are already
C   called before the user routine) but users may want to
C   change parameters and re-run them
C   !!!!!    !!!!!     !!!!!     !!!!!     !!!!!     !!!!! 
C
C************************************************************/ 
C* The following lines call the standard routines           */ 
C* which compute E/p, the COG and the pt**2                 */
C* for charged tracks.                                      */ 
C* These routines fill data in the compact structure        */
C************************************************************/ 
      call USER_STD_EP(IB, IE, 0)
      call USER_STD_RCOG(IB, IE, 0)
      call USER_STD_PTSQ(IB, IE, 0)

C************************************************************/ 
C* The following line to call the routine                   */
C* which flags GOOD vertices                                */
C* One data needs to be filled before calling the routine   */
C* ACHARG_IVERTEX with the vertex index                     */
C* Two data are filled by the routine                       */
C* ACHARG_KTYPE with 1 for Ks and 2 for Kl                  */
C* ACHARG_IFLAG with 1 if OK, 0 not                         */
C************************************************************/

      do IV=1,EVT_NVERTEX
         ACHARG_IVERTEX=iv
         call USER_SELCHARGED(IB,QB,IE,QE)
C************************************************************/ 
C* The following lines to call the routine                  */
C* which computes time given by Ch. hodoscopes              */
C* Two data need to be filled before calling the routine    */
C* ACHOD_ITRACE(2) with the index of the first track        */ 
C* ACHOD_ITRACE(2) with the index of the second track       */ 
C* evt->achod.hodotime is compute by the routine            */
C************************************************************/

         ACHOD_ITRACE(1)=VERTEX_IPTRK(iv)
         ACHOD_ITRACE(2)=VERTEX_INTRK(iv)
         ICORR=0
         iok = USER_HODOTIME(IB,QB,IE,QE,ICORR)

      enddo


C************************************************************/ 
C* The following line to call the routine                   */
C* which flags GOOD 2pi0 events                             */
C* One data needs to be filled before calling the routine   */
C* ANEUT_ESC to set the energy scale                        */ 
C* ANEUT_IFLAG is filled by the routine: 1 OK, 0 not        */
C************************************************************/

       ANEUT_ESC=1.
       call USER_SEL2PI0(IB,QB,IE,QE)
       print *,' >>>>>>> iflag NEUTRAL - CHARGED',ANEUT_IFLAG, ACHARG_IFLAG

C************************************************************/ 
C* The following line to call the routine                   */
C* which computes the event time using the GOOD 2pi0 events */
C* selected by the SEL2PI0 routine                          */
C* In the routine, it is assumed that the 4 gamma clusters  */
C* are defined in ANEUT_INC(4).                             */
C************************************************************/

       if(ANEUT_IFLAG .gt. 0) THEN
          call USER_LKRTIME(IB,QB,IE,QE)
          print *,' >>>>>> LKRtime=',ANEUT_LKRTIME,' Ncells=',ANEUT_NTUSED
       endif

C************************************************************/ 
C* The following lines to call the routine                  */
C* which computes proton time given by the tagger.          */
C* The proton time is returned into ATAG_APROT_PROTIME(i)   */
C************************************************************/

      iok = USER_TAGTIME(IB,QB,IE,QE)
      print *,' >>>>> Nprotons=',ATAG_NAPROT

      if((ATAG_NAPROT).gt.0) then
         print *,' >>>>> 1st Proton Time / 1st Proton Nhits / hits',
     +        ATAG_APROT_PROTTIME(1), ATAG_APROT_NPROTHIT(1),
     +        ATAG_APROT_PROTHIT_INDEXHIT(1,1,1), ATAG_APROT_PROTHIT_INDEXHIT(1,1,2)
      endif




C******************************
C* Cuts applied to the events *
C******************************
      CALL HFILL(100,0.,0.,1.)      ! initialize cut histogram
      IF(EVT_NTRACK.NE.2) RETURN        ! check for 2 tracks
      CALL HFILL(100,1.,0.,1.)      ! store events with 2 tracks

C Check that there is one positive and one negative track

      IF((TRACK_Q(1)*TRACK_Q(2)).GT.0.) RETURN ! check for + and - tracks
      CALL HFILL(100,2.,0.,1.)      ! store events with + and - tracks

C Cut out lamdas

      IF(((TRACK_P(2)/TRACK_P(1)).LT.P_LAMBDA).AND.
     +   ((TRACK_P(1)/TRACK_P(2)).LT.P_LAMBDA)) RETURN
      CALL HFILL(100,3.,0.,1.)      ! non-lambda events

C Check that each track has momentum greater than pMin GeV

      IF((TRACK_P(1).LT.PMIN).OR.(TRACK_P(2).LT.PMIN)) RETURN
      CALL HFILL(100,4.,0.,1.)      ! track momentum > PMIN

C Check closest distance of approach of the tracks

      IF(VERTEX_CDA(1).GT.CDA_CUT) RETURN
      CALL HFILL(100,5.,0.,1.)      ! CDA < CDA_CUT

C Check Z position of vertex

      IF((VERTEX_Z(1).LT.ZMIN).OR.(VERTEX_Z(1).GT.ZMAX)) RETURN
      CALL HFILL(100,6.,0.,1.)      ! Z vertex in range
 
C Check for 1-2 clusters in LKr

      IF((EVT_NLKR.GT.2).OR.(EVT_NLKR.LT.1)) RETURN
      CALL HFILL(100,7.,0.,1.)      ! liquid krypton cluster check

C ----------- End of user F77 code -----------
      RETURN
      END
