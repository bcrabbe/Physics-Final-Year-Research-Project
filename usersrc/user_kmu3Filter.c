/********************************************************/
/* COmPACT user routine: user_kmu3Filter(kmu3Event *evt)*/
/*                                                      */
/* This routine acts as a filter for all kmu3Event data.*/
/* If it returns a value <0 the event will be sent to   */
/* an output file if one has been opened.               */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "reader.h"

int user_kmu3Filter(Burst *bur,kmu3Event *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
/*----------- End of user C code -----------*/
  return 0;
}
