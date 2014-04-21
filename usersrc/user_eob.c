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
/*
  fprintf(fprt," ===> eob \n");
  printEndofBurst(eob,fprt);
  printTSTWList(bur,fprt);
*/

/*----------- End of user C code -----------*/
  return 0;
}
