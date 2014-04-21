/********************************************************/
/* COmPACT user routine: user_cmpFilter(cmpEvent *evt)  */
/*                                                      */
/* This routine acts as a filter for all cmpEvent data. */
/* If it returns a value <0 the event will be sent to   */
/* an output file if one has been opened.               */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "reader.h"

#define OK  0
#define NOK 1

int user_cmpFilter(Burst *bur,cmpEvent *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/

  return -1;
}
