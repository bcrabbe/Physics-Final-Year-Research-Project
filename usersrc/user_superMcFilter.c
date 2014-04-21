/****************************************************************/
/* COmPACT user routine: user_superMcFilter(superMcEvent *evt)  */
/*                                                              */
/* This routine acts as a filter for all superMcEvent data.     */
/* If it returns a value <0 the event will be sent to           */
/* an output file if one has been opened.                       */
/*                                   BH 10/2/98    RWM 20/6/97  */
/****************************************************************/

#include "cmpio.h"
#include "reader.h"

int user_superMcFilter(superBurst *sbur,superMcEvent *sevt,hyperCmpEvent *hevt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
/* Remember to fill hyperCOmPACT structure if you want hcmp output */

/*----------- End of user C code -----------*/
  return -1;
}
