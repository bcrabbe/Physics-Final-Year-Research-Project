/********************************************************/
/* COmPACT user routine: user_burst(Burst *bur)         */
/*                                                      */
/* User routine called everytime a burst `*bur' is      */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 20/6/97 */
/* Modified: 15-02-99 IWS change default in example     */
/* ========                                             */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"

int user_burst(Burst *bur) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  static int Nburst=0;
  static int RunCurrent=0;
  if(!Nburst) PrintCmpGeom(Geom, fprt);
  if(Nburst++ < 5) printBurst(bur,fprt);


/*************************************************************/ 
/* The following lines show how to set up the selection      */
/* of analysis routines.                                     */
/* The default is set in compact to call SOME analysis       */
/* routines; if the user wants to change the default, here is*/
/* the way to go.                                            */
/* The code which calls the analysis routines is in          */
/* /afs/cern.ch/na48/offline2/compact/compact-4.3/compact/   */
/*       /rlib/src/cmpAnalysis.c                             */
/* The routine USER_TAGTIME is called by default             */
/* only for good neutral or charged events                   */
/* The following example is equivalent to setting            */
/* bur->CallAnaRoutine=1 i.e. letting compact decide which   */
/* routines to call.                                         */
/* If the user sets bur->CallAnaRoutine=0 NO routine         */
/* will be called from compact.                              */
/* Some routines are called ANYWAY:                          */
/* USER_BADBURST                                             */
/* USER_ESPY, USER_STD_EP,USER_STD_RCOG,USER_STD_PTSQ        */
/*************************************************************/
  bur->CallAnaRoutine=-1;   /* the user wants to choose the routines */
  bur->tocall.selcharged=0; /* to call selcharged */
  bur->tocall.sel2pi0=1;    /* to call sel2pi0 */
  bur->tocall.sel3pi0=1;    /* to call sel3pi0 - only if no 2pi0 cand.*/
  bur->tocall.bluefield=0;  /* NOT to call correction for Bfield in blue tube (1997) */
  bur->tocall.lkrpedcor=0;  /* NOT to call lkr pedestal correction (1997) */
  bur->tocall.lkrposcor=1;  /* to call lkr position correction  */
  bur->tocall.lkrsharing=0; /* NOT to call lkr energy sharing correction */
  bur->tocall.hodotime=1;   /* to call hodotime */
  bur->tocall.nhodtime=1;   /* to call nhodtime */
  bur->tocall.lkrtime=1;    /* to call lkrtime only if there is a 2pi0 cand. */
  bur->tocall.tagtime=1;    /* to call tagtime */
  bur->tocall.aksflag=1;    /* to call flagging of good AKS events */
  bur->tocall.muon_rec=0;   /* NOT IN COMPACT anymore */
  bur->tocall.muon_reject=0;/* NOT to call muon Rejection routines (1997)-NOT USEFULL */
  bur->tocall.geomcor=1;    /* to call geometry correction routine (1997+1998) */
  bur->tocall.magnetcor=1;  /* to call magnet correction routine  (adjust the Kaon mass to be flat with time) */
  bur->tocall.newcharged=1; /* to call the new charged analysis routine */
  bur->tocall.muon_trackrec=1;   /* to call "track" muon reconstruction  - THE RECOMMANDED ONE */
  bur->tocall.muon_vtxrec=0;     /* NOT to call "vertex" muon reconstruction (1997) */
  bur->tocall.lkrcalcor=1;  /* to call lkr calibration correction */
  bur->tocall.sel2gam=1;    /* to call sel2gam - only if no 2 or 3 pi0 cand.*/
  bur->tocall.lkrcalhi2k=0; /*NOT to call lkr corrections (change of temperature) for KS HI 2000 data */

/*************************************************************/ 
/* The following lines show how to reject events according   */
/* to the "bad burts" flags - The routine USER_BADBURST      */
/* SHOULD have been called before.                           */
/* From version 4.0 of compact (march 1998) the user_badburst*/
/* routine is called by compact.                             */
/* From version 4.1 of compact (may 1998) if the flag        */
/* bur->BadB.Skip > 0: the entire burst is skipped by compact*/
/* This should NOT be used for filtered files where there are*/
/* more than ONE burst.                                      */
/* It is up to the user to decide on which criteria to reject*/
/* a burst; the example given here uses all the available    */
/* flags to reject a burst and also require the burst        */
/* to have eps'/eps data.                                    */
/*************************************************************/ 
/*
  Reject events for which relevant detector are not OK
  For epsilon'/epsilon analysis 
*/

  if( (bur->BadB.Tag<0) || 
      (bur->BadB.Aks<0) || 
      (bur->BadB.Dch<0) || 
      (bur->BadB.Mbx<0) || 
      (bur->BadB.HodC<0) || 
      (bur->BadB.Lkr<0) ||
      (bur->BadB.Nut<0) ||
      (bur->BadB.Muv<0) ||
      (bur->BadB.PMB<0) ||
      (bur->BadB.Clk<0) ||
      (bur->BadB.NoEps < 0) ||
      (bur->BadB.spare2<0) ) bur->BadB.Skip=1; /* spare2 makes equal eps' sample Lyon & CERN */

  if(bur->BadB.Skip>0) return 0;

/******************************************************************/
/* The following lines show how to call the routine               */
/* to print bad burst summary                                     */
/******************************************************************/
  printBadBurst(bur,fprt);
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of LKR dead cells, for the current run            */
/******************************************************************/
  if(bur->nrun != RunCurrent)
    {
      RunCurrent=bur->nrun;
      CmpLkrDeadPrint(fprt);
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of timing offsets, for the current run            */
/******************************************************************/
      CmpTimeOffsetPrint(bur,fprt);
    }
/*******************************************************************/
/* The following lines show how to call the routine to empty part  */
/* of the burst structure in order to save space when filtering    */
/* events                                                          */
/*******************************************************************/
      EMPTY_BURST(bur,100);
/*----------- End of user C code -----------*/
  return 0;
}
