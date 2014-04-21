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
  static int RunCurrent=0;
  static int Nburst=0;

  if(sbur->nrun != RunCurrent){
      RunCurrent=sbur->nrun;
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of timing offsets, for the current run            */
/******************************************************************/
      CmpTimeOffsetPrintSC(sbur,fprt); 
    }
/*
  Reject events for which relevant detector are not OK
  For epsilon'/epsilon analysis 
*/
  if (
      (sbur->BadB.Lkr!=0)
      || (sbur->BadB.Dch!=0)
      || (sbur->BadB.Nut!=0)
      || (sbur->BadB.Mbx!=0)
      || (sbur->BadB.Hac!=0)
      || (sbur->BadB.Muv!=0)
      || (sbur->BadB.HodC!=0)
      || (sbur->BadB.HodN!=0)
      || (sbur->BadB.PMB!=0)
      || (sbur->BadB.Akl!=0)
      || (sbur->BadB.Clk!=0)
      /*  || (sbur->BadB.Ksm!=0) */
      || (sbur->BadB.Kab!=0)
      || (sbur->BadB.Phys!=0)
      ) sbur->BadB.Skip = 1;

  if(sbur->BadB.Skip>0) return 0;

/* To print out some of the database data */
/*
  fprintf(fprt," **** SuperBurst: Run/time/intensity=%5d/%9d/%6d **** \n",
          sbur->nrun, sbur->time, sbur->intensity);
  if(rdb != NULL)
    {
      printf(" From superBurst nBur=%d nBurBad=%d \n",rdb->nBur, rdb->nBurBad);
    }
  if(bdb != NULL)
    {
      printf(" From superBurst BURST=%4d %4d %4d %4d \n",
             bdb->nFiltCh, bdb->nGoodCh, bdb->nFiltNe,  bdb->nGoodNe);   
    }
*/
/*----------- End of user C code -----------*/
  return 0;
}
