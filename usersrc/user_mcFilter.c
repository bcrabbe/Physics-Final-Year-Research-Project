/********************************************************/
/* COmPACT user routine: user_mcFilter(mcEvent *evt)    */
/*                                                      */
/* This routine acts as a filter for all mcEvent data.  */
/* If it returns a value <0 the event will be sent to   */
/* an output file if one has been opened.               */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"

int user_mcFilter(Burst *bur,mcEvent *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
/*----------- End of user C code -----------*/
  return -1;
}
