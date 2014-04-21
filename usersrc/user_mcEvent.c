/********************************************************/
/* COmPACT user routine: user_mcEvent(mcEvent *evt)     */
/*                                                      */
/* User routine called everytime an event `*evt' is     */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 11/7/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"  /* to get procedure calls for fortran routines */
#include "F77_ana.h" /* mapping of fortran common blocks for analysis routines */

int user_mcEvent(Burst *bur,mcEvent *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  static int nuserevt=0;
  int icorr, i;

  if(nuserevt<20)
    {
      printMcEvent(evt,fprt);
      nuserevt++;
    }

/* 
   Uncomment following lines to call the hodotime and tagtime functions 
*/
/*
  icorr=0; Dont use corrected variables 
  USER_HODOTIME(bur,evt->cmpevt,&icorr);
  USER_TAGTIME(bur,evt->cmpevt);
*/
/*----------- End of user C code -----------*/
  return 0;
}
