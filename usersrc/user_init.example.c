/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "cfortran.h"
#include "hbook.h"
#include "user.h"

#define NPAW 2000000

float pawc_[NPAW];

extern void hlimit_(int *);

int user_init() {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  int istat,recl=4096;
  int tempval=NPAW;           /* size of PAWC common block */
  hlimit_(&tempval);
  fprt=fopen("compact.txt","w");
  fprintf(fprt," -------> entering user_init\n");
  HBOOK1(100,"Particle Cuts",21,-0.5,20.5,0.0);
/*----------- End of user C code -----------*/
  return 0;
}
