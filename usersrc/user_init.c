#define ENABLE_OUTPUT 1

/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"

#include "reader.h"

#include "user.h"
//char * dir = "/afs/cern.ch/user/b/bcrabbe/kl3-analysis/1/outputs/1/";
#if ENABLE_OUTPUT

//char * realFileName = "real-data.otxt";
//char * realFullPath = NULL;
FILE * realFP;

//char * ke3FileName = "ke3-data.otxt";
//char * ke3FullPath = NULL;
FILE * ke3FP;

//char * km3FileName = "km3-data.otxt";
//char * km3FullPath = NULL;
FILE * km3FP;

//char * k2piFileName = "k2pi-data.otxt";
//char * k2piFullPath = NULL;
FILE * k2piFP;

FILE *k3piFP;
FILE *km2FP;
FILE *k3pi0FP;
FILE *FP1;
FILE *FP2;
FILE *cutWhichKilledEventFP;
#endif

    int ke3Count =0;
    int km3Count = 0;
    int k2piCount = 0;
    int numUnidentified =0;
    int  nEvents =0;

    int MCke3Count =0;
    int MCkm3Count = 0;
    int MCk2piCount = 0;
    int MCk3piCount =0;
    int MCnumUnidentified =0;
    int  MCnEvents =0;

float CPDpos_leftDownCorner[256][2];
float CELLpos_leftDownCorner[256][64][2];
  /*) int ke3Count=0;
    int km3Count=0;
   int k2piCount=0;
    int nEvents=0;*/
//int   CPDindex, CELLindex;

//float CELLlength = 1.975;
//float CPDlength = 8 * 1.975;

// E/p corrections for each cell
FILE *EopCorrfile;
float EopCorr[256][64];

int user_init() {

/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
	fprt=fopen("compact.txt","w");
#if ENABLE_OUTPUT
	/*NOTE: asprintf isn't portable (if we stay with gcc all is good)*/
	/*
	asprintf(&realFullPath,"%s%s",dir,realFileName);
	realFP = fopen(realFullPath,"w");

	asprintf(&ke3FullPath,"%s%s",dir,ke3FileName);
	ke3FP = fopen(ke3FullPath,"w");

	asprintf(&km3FullPath,"%s%s",dir,km3FileName);
	km3FP = fopen(km3FullPath,"w");

	asprintf(&k2piFullPath,"%s%s",dir,k2piFileName);
	k2piFP = fopen(k2piFullPath,"w");
        */
	FP1 = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/test1.txt", "w");
	FP2 = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/test2.txt", "w");
	realFP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/real-data.otxt", "w");
	ke3FP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/ke3-data.otxt", "w");
	km3FP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/km3-data.otxt", "w");
	k2piFP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/k2pi-data.otxt","w");
	k3piFP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/k3pi-data.otxt","w");
	k3pi0FP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/k3pi0-data.otxt","w");
	km2FP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/km2-data.otxt","w");
    cutWhichKilledEventFP = fopen("/afs/cern.ch/user/b/bcrabbe/private/kl3-analysis/1/outputs/1/cutWhichKilledEvent.txt","w");

    /*int ke3Count =0;
    int km3Count = 0;
    int k2piCount = 0;
    int numUnidentified =0;
    int  nEvents =0;*/


	// CPD stuff:

    int i, j, k;
    int l, m, n;
    int cpd, cell;
    float CELLlength = 1.975;
    float CPDlength = 8 * 1.975;

    // Define the positions for CPDs and Cells and store them
    for (i=0; i<16; i++)
         for (j=0; j<16; j++)
    {
      k = i*16 + j;
      CPDpos_leftDownCorner[k][0] = (-1)*(7-i)*CPDlength;  // LKr RF is left-handed  --> the x sign has to be changed !!!
      CPDpos_leftDownCorner[k][1] = (7-j)*CPDlength;
      //printf ("CPD %d: position left down corner = %.2f, \t%.2f\n", k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1]);

      for (m=0; m<8; m++)
        for (n=0; n<8; n++)
          {
        l = m*8 + n;
        CELLpos_leftDownCorner[k][l][0] = CPDpos_leftDownCorner[k][0] - (7-m)*CELLlength;
        CELLpos_leftDownCorner[k][l][1] = CPDpos_leftDownCorner[k][1] + (7-n)*CELLlength;
        //printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n",
        //      l, k, CELLpos_leftDownCorner[k][l][0], CELLpos_leftDownCorner[k][l][1]);
          }
    }

    // Read E/p corrections for each cell from file
   // int    cpd, cell;

    //for first part of p5:
    EopCorrfile = fopen ("eopCorrfile_p5_20410-20415_v63.dat", "r");

    //2nd part:
    //EopCorrfile = fopen ("eopCorrfile_p5_20416-20485_v63.dat", "r");


    for (i=0; i<256; i++)
         for (j=0; j<64; j++)
    {
      fscanf (EopCorrfile, "%i %i %f\n", &cpd, &cell, &EopCorr[i][j]);
      //printf ("cpd%i\tcell%i\t%f\n", cpd, cell, EopCorr[i][j]);
    }

    fclose(EopCorrfile);
#endif
/*----------- End of user C code -----------*/
  return 0;
}



