/********************************************************/
/* COmPACT user routine: user_ke3Event(ke3Event *evt)   */
/*                                                      */
/* User routine called everytime an event `*evt' is     */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"

int user_ke3Event(Burst *bur,ke3Event *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
/*----------- End of user C code -----------*/
  static int nuserevt=0;
  
  if(nuserevt<20) {
    printKe3Event(evt,fprt);
    nuserevt++;
  }
  return 0;
}
