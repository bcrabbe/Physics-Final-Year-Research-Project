/********************************************************/
/* COmPACT user routine: user_exit()                    */
/*                                                      */
/* User routine called once all data has been read. This*/
/* allows output files to be closed and any other user  */
/* resources to be tidied up.                           */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "cfortran.h"
#include "hbook.h"

int user_exit() {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  int istat,recl=4096;
  fprintf(fprt," -------> entering user_exit\n");
  HROPEN(1,"hbex","compact.hbook","N",recl,istat);
  HROUT(0,istat," ");
  HREND("hbex");
  fclose(fprt);
/*----------- End of user C code -----------*/
  return 0;
}
