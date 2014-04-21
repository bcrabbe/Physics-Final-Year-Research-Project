#include "cmpio.h"
#include "reader.h"
#include "user.h"

int user_hyperBurst(hyperBurst *hbur){
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  hbur->BadB.Skip = 0; /* see user_superBurst.example.c to learn to use it */

/*----------- End of user C code -----------*/
  return 0;
}
