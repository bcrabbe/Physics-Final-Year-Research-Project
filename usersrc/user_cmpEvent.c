/********************************************************/
/* COmPACT user routine: user_cmpEvent(cmpEvent *evt)   */
/*                                                      */
/* User routine called everytime an event `*evt' is     */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 20/6/97 */
/********************************************************/
#include "cmpio.h"
#include "reader.h"
#include "user.h"

int user_cmpEvent(Burst *bur,cmpEvent *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  int i, itrk,icorr;
  static int nuserevt=0;
  if( (nuserevt<20) )
    {
      fprintf(fprt,"-------------------------------------------------------\n");
      fprintf(fprt," Selection flags : Charged %1d| Neutral %1d               |\n",
              evt->acharg.iflag, evt->aneut.iflag);
      fprintf(fprt,"-------------------------------------------------------\n");
      PrintFlagCorr(evt->FlagCorr,fprt);
      printCmpEvent(evt,fprt);
    }
  nuserevt++;
    
/* Un-Comment following line to reject overlayed events */
/*
   if(evt->rndm.type > 0) return 0; 
*/

/*----------- End of user C code -----------*/
  return 0;
}
/*
  printf(" ===> ts: %9d/%9d/%9d tw=%8x/%8x/%8x Sim=%f %f \n",
         evt->spareInt[1], evt->timestamp, evt->spareInt[3],
         evt->spareInt[2], evt->trigword, evt->spareInt[4],
         evt->ovrflwSimBef,evt->ovrflwSimAft); 
*/
