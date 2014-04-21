/********************************************************/
/* COmPACT user routine: user_burst(Burst *bur)         */
/*                                                      */
/* User routine called everytime a burst `*bur' is      */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"
int user_burst(Burst *bur) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  static int Nburst=0;
  static int RunCurrent=0;
  if(!Nburst) PrintCmpGeom(Geom, fprt);
  Nburst++;
/*  if(Nburst++ < 5) printBurst(bur,fprt);*/
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of LKR dead cells, for the current run            */
/******************************************************************/
  if(bur->nrun != RunCurrent)
    {
      RunCurrent=bur->nrun;
      CmpLkrDeadPrint(fprt);
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of timing offsets, for the current run            */
/******************************************************************/
      CmpTimeOffsetPrint(bur,fprt);
    }
/* To print out some of the database data */
/*
  printf(" From Burst nBur=%d nBurBad=%d \n",rdb->nBur, rdb->nBurBad);
  printf(" From Burst BURST=%4d %4d %4d %4d \n",
         bdb->nFiltCh, bdb->nGoodCh, bdb->nFiltNe,  bdb->nGoodNe);   
*/
/*----------- End of user C code -----------*/
  return 0;
}
