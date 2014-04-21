/********************************************************/
/* COmPACT user routine: user_cmpEvent(cmpEvent *evt)   */
/*                                                      */
/* User routine called everytime an event `*evt' is     */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 20/6/97 */
/********************************************************/
#include <math.h>
#include "cmpio.h"
#include "cfortran.h"
#include "hbook.h"
#include "user.h"
#include "reader.h"   /* function definitions */
/****************************************/
/* Definitions for cuts used in Example */
/****************************************/
#define PMIN       5.0      /* minimum track momentum */
#define CDA_CUT    1.4      /* CDA maximum in cm */
#define ZMIN       0.0      /* minimum vertex position (KStarget = 0cm) */
#define ZMAX    9000.0      /* maximum vertex position */
#define RMIN_LKR  17.0      /* minimum radius of LKR */
#define RMAX_LKR 110.0      /* maximum radius of LKR */
#define P_LAMBDA   3.0      /* Track momenta ratio cut to remove lambdas */
#define MAX_TRACKS 5        /* maximum number of tracks */
#define ASSOCLKR  10.0      /* max distance for LKR cluster association */
#define ASSOCHAC  15.0      /* max distance for HAC cluster association */

#define POW2(x) ((x)*(x))

int user_cmpEvent(Burst *bur,cmpEvent *evt) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  int i,j,iclus=0;
  static int nuserevt=0;
  int iv;                             /* index for vertex loop */
  int icorr;                          /* to select original or corrected structures */
/******************************/
/* Cuts applied to the events */
/******************************/

  printf(" entering user_cmpEvent evt->Nlkr=%3d\n",evt->Nlkr);
  fprintf(fprt," entering user_cmpEvent evt->Nlkr=%3d\n",evt->Nlkr);
  if(nuserevt<20) {
    
    fprintf(fprt,"-------------------------------------------------------\n");
    fprintf(fprt," Selection flags : Charged %1d| Neutral %1d               |\n",
            evt->acharg.iflag, evt->aneut.iflag);
    fprintf(fprt,"-------------------------------------------------------\n");
    PrintFlagCorr(evt->FlagCorr,fprt);
    printCmpEvent(evt,fprt);
  }
  nuserevt++;

/************************************************************/ 
/* Please read the comments concerning the analysis routines*/
/* To call these routines, please add the following line    */
/* at the begining of the file                              */
/* #include "reader.h"   function definitions               */
/************************************************************/ 

/************************************************************/ 
/* If the input file contains MORE than ONE burst, then the */
/* selction mecanism for "good" burst cannot be called from */
/* the burst routine but has to be called from the event    */
/* as shown below                                           */
/************************************************************/ 
/*
  Reject events for which NUT and MBX and TAG and LKR and MUV not OK
*/
  if( (bur->BadB.Nut<0) || 
      (bur->BadB.Mbx<0) || 
      (bur->BadB.Tag<0) ||
      (bur->BadB.Muv<0) ||
      (bur->BadB.Lkr<0) ) return;
/* 
For epsilon'/epsilon analysis 
*/
  if(bur->BadB.NoEps<0) return;

/************************************************************/ 
/* The following lines call the NEW MUV reconstruction      */
/* routine; if it is called, the content of the             */
/* MUV structure will be overwritten by the result of this  */
/* routine.                                                 */
/* Therefore SELCHARGED has to be called again              */
/************************************************************/ 

  if(bur->nrun < 6000)
    {
      USER_NEW_MUON(bur, evt);
    }

/* 
   Compare energies before and after LKR ped cor
*/
  for(i=0; i<evt->Nlkr; i++)
    {
      printf(" ================> ecorrke3=%5.2f (before) --> %5.2f (after)\n",
             evt->lkr[i].ecorrke3, evt->lkr[i].anavar[0]);
    }


/* 
   !!!!!    !!!!!     !!!!!     !!!!!     !!!!!     !!!!! 
   The following lines just show how to call the analysis 
   routines; they dont need to be called (as they are already
   called before the user routine) but users may want to
   change parameters and re-run them
   !!!!!    !!!!!     !!!!!     !!!!!     !!!!!     !!!!! 
*/

/************************************************************/ 
/* The following lines call the standard routines           */ 
/* which compute E/p, the COG and the pt**2                 */
/* for charged tracks.                                      */ 
/* These routines fill data in the compact structure        */
/************************************************************/ 

  /* 
     The third argument (0) is to use original tracks - 
                        (1) for corrected 
  */
  USER_STD_EP(bur, evt, 0);  /* Fills evt->track[i].EovP, LKRclu */
  USER_STD_RCOG(bur, evt, 0);/* Fills evt->vertex[i].rcog */
  USER_STD_PTSQ(bur, evt, 0);/* Fills evt->vertex[i].ptsqks, ptsqkl */


/************************************************************/ 
/* The following line to call the routine                   */
/* which flags GOOD vertices                                */
/* One data needs to be filled before calling the routine   */
/* evt->acharg.ivertex with the vertex index                */
/* Two data are filled by the routine                       */
/* evt->acharg.ktype with 1 for Ks and 2 for Kl             */
/* evt->acharg.iflag with 1 if OK, 0 not                    */
/************************************************************/

  for(iv=0; iv<evt->Nvertex; iv++)
    {
      evt->acharg.ivertex=iv+1; /* The vertex index is used in a fortran routine */
      USER_SELCHARGED(bur,evt); 
      printf(" vertex %2d is a %1d type vertex rcog=%f\n",
             evt->acharg.ivertex, evt->acharg.ktype, evt->vertex[iv].rcog);
/************************************************************/ 
/* The following lines to call the routine                  */
/* which computes time given by Ch. hodoscopes              */
/* Two data need to be filled before calling the routine    */
/* evt->achod.itrace[0] with the index of the first track   */ 
/* evt->achod.itrace[1] with the index of the second track  */ 
/* evt->achod.hodotime is compute by the routine            */
/************************************************************/
      evt->achod.itrace[0]=evt->vertex[iv].iptrk+1;
      evt->achod.itrace[1]=evt->vertex[iv].intrk+1;
      icorr=0;
      USER_HODOTIME(bur,evt,&icorr); /* Call HODOtime routine for original tracks */
      printf(" ---> hodotime=%f\n", evt->achod.hodotime);
    }


/************************************************************/ 
/* The following line to call the routine                   */
/* which flags GOOD 2pi0 events                             */
/* One data needs to be filled before calling the routine   */
/* evt->aneut.esc to set the energy scale                   */ 
/* evt->aneut.iflag is filled by the routine: 1 OK, 0 not   */
/************************************************************/

  evt->aneut.esc=1.;
  USER_SEL2PI0(bur, evt);

  printf(" ---> Selection flags neutral=%2d charged=%2d\n",
         evt->aneut.iflag, evt->acharg.iflag);

/************************************************************/ 
/* The following line to call the routine                   */
/* which computes the event time using the GOOD 2pi0 events */
/* selected by the SEL2PI0 routine                          */
/* In the routine, it is assumed that the 4 gamma clusters  */
/* are defined in ANEUT_INC(4).                             */
/************************************************************/

    if(evt->aneut.iflag > 0)
      {
        USER_LKRTIME(bur,evt);
        printf(" >>>> LKR time=%f Ncells=%d\n",
               evt->aneut.LKRtime,evt->aneut.ntUsed);
      }

/************************************************************/ 
/* The following lines to call the routine                  */
/* which computes proton time given by the tagger           */
/* This routine is not yet properly into compact and returns*/
/* data inside a common block which is accessible by the    */
/* C struture timtag_ ; file F77_ana.h needs to be included */
/************************************************************/
  USER_TAGTIME(bur,evt);
  printf(" ---> Nproton =%d\n", evt->atag.Naprot);
  if(evt->atag.Naprot > 0) 
    printf(" ---> time=%f hits=%d %d \n",
           evt->atag.aprot[0].protTime, 
           evt->atag.aprot[0].prothit[0].indexhit[0],
           evt->atag.aprot[0].prothit[0].indexhit[1]);



/*************************************************************/


  HFILL(100,0.,0.,1.);                /* initialize cut histogram */
  if(evt->Ntrack != 2)               /* check for 2 tracks */
    return 0;
  HFILL(100,1.,0.,1.);                /* store events with 2 tracks */

/* Check that there is one positive and one negative track */
  if((evt->track[0].p*evt->track[1].p)>0.) /* check for + and - tracks */
    return 0;
  HFILL(100,2.,0.,1.);                /* store events with + and - tracks */

/* Cut out lamdas */
  if(fabs(evt->track[1].p/evt->track[0].p) > P_LAMBDA ||
     fabs(evt->track[0].p/evt->track[1].p) > P_LAMBDA)
    return 0;
  HFILL(100,3.,0.,1.);                /* non-lambda events */

/* Check that each track has momentum greater than pMin GeV */
  if((evt->track[0].p < PMIN) || (evt->track[1].p < PMIN))
    return 0;
  HFILL(100,4.,0.,1.);                /* track momentum > PMIN */

/* Check closest distance of approach of the tracks */
  if(evt->vertex[0].cda > CDA_CUT)
    return 0;
  HFILL(100,5.,0.,1.);                /* CDA < CDA_CUT */

/* Check Z position of vertex */
  if((evt->vertex[0].z < ZMIN) || (evt->vertex[0].z > ZMAX)) 
    return 0;
  HFILL(100,6.,0.,1.);                /* Z vertex in range */
 
/* Check for 1-2 clusters in LKr */
  if((evt->Nlkr > 2) || (evt->Nlkr < 1))
    return 0;
  HFILL(100,7.,0.,1.);                /* liquid krypton cluster check */

/*----------- End of user C code -----------*/
  return 0;
}
