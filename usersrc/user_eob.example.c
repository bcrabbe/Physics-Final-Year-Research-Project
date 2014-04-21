/***************************************************************/
/* COmPACT user routine: user_eob(Burst *bur, EndofBurst *eob) */
/*                                                             */
/* User routine called at the end of every burst and is passed */
/* the burst header and burst trailer.                         */
/*                                                             */
/*                                                 RWM 20/6/97 */
/***************************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"

int user_eob(Burst *bur,EndofBurst *eob) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  printEndofBurst(eob,fprt);
/* 
   Example on how to call the EOB error decoding routine
*/
  USER_EOBDEC(eob);
/*******************************************************************/
/* The following lines show how to call the routine to empty part  */
/* of the End of burst structure in order to save space when       */
/* filtering events                                                */
/*******************************************************************/
/*
    EMPTY_EOB(eob,0);
*/
/************************************************************/ 
/* The following line shows how to call the routine to move */
/* the TS+TW lists from ENDofBURST to BURST structure       */
/************************************************************/ 
      getEoBList(bur,eob);
/*----------- End of user C code -----------*/
  return 0;
}
