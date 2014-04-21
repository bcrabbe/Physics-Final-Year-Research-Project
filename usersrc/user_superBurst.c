/***********************************************************/
/* COmPACT user routine: user_superBurst(superBurst *sbur) */
/*                                                         */
/* User routine called everytime a SuperCOmPACT burst      */
/* `*bur' is loaded.                                       */
/* A return value of greater than zero denotes an error    */
/* condition has occured.                                  */
/*                                             BH 2/3/98   */
/***********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"

int user_superBurst(superBurst *sbur) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  sbur->BadB.Skip = 0; /* see user_superBurst.example.c to learn to use it */

/*----------- End of user C code -----------*/
  return 0;
}
