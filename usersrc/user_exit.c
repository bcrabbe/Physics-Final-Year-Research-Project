#define ENABLE_OUTPUT 1
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
#include "reader.h"
/*
char * realDataK2piFileName = "realData-k2pi.otxt";
char * realDataK2piFullPath = NULL;
FILE * realDataK2piFP;

char * ke3DataK2piFileName = "realData-k2pi.otxt";
char * ke3DataK2piFullPath = NULL;
FILE * ke3DataK2piFP;

char * km3DataK2piFileName = "realData-k2pi.otxt";
char * km3DataK2piFullPath = NULL;
FILE * km3DataK2piFP;

char * k2piDataK2piFileName = "realData-k2pi.otxt";
char * k2piDataK2piFullPath = NULL;
FILE * k2piDataK2piFP;
*/

int user_exit() {
  /* WARNING: do not alter things before this line */
  /*---------- Add user C code here ----------*/
#if ENABLE_OUTPUT
    fclose(realFP);
    fclose(ke3FP);
    fclose(km3FP);
    fclose(k2piFP);

    fclose(km2FP);
    fclose(k3pi0FP);
    fclose(k3piFP);

    fclose(FP1);
    fclose(FP2);
    fclose(cutWhichKilledEventFP);

    float fracKe3 = (float)ke3Count/(float)nEvents;
    float fracKm3 = (float)km3Count/(float)nEvents;
    float fracK2pi = (float)k2piCount/(float)nEvents;
    float fracUnidentified = (float)numUnidentified/(float)nEvents;

    if(MCnEvents==0)
    {
        printf("\n\n%d signal Events \n",nEvents);
        printf("Ke3: %d Events, %f percent\n",ke3Count, fracKe3*100);
        printf("Km3: %d Events, %f percent\n",km3Count, fracKm3*100);
        printf("K2Pi: %d Events, %f percent\n",k2piCount,fracK2pi*100);
        printf("%d were unidentified, %f percent\n\n",numUnidentified,fracUnidentified);
    }
    if(MCnEvents>0)
    {
        printf("\n\n%d MC events\n",MCnEvents);
        printf("%d ke3, %d km3, %d k2pi and %d k3pi\n",MCke3Count,MCkm3Count,MCk2piCount,MCk3piCount);
        printf("Ke3: %d Events, %f percent \n",ke3Count, fracKe3*100);
        printf("Km3: %d Events, %f percent\n",km3Count, fracKm3*100);
        printf("K2Pi: %d Events, %f percent\n",k2piCount,fracK2pi*100);
        printf("%d were unidentified, %f percent\n\n",numUnidentified,fracUnidentified*100);
    }




  /*
  fclose(realMmFP);
  fclose(ke3MmFP);
  fclose(km3MmFP);
  fclose(k2piMmFP);
  */
  /*
  fclose(realEPRatioFP);
  fclose(ke3EPRatioFP);
  fclose(km3EPRatioFP);
  fclose(k2piEPRatioFP);
  */
/*
  fclose(mPionMkaonDeltaPFP);
  fclose(ke3mPionMkaonDeltaPFP);
  fclose(km3mPionMkaonDeltaPFP);
  fclose(k2pimPionMkaonDeltaPFP);

  asprintf(&realDataK2piFullPath,"%s%s",dir,realDataK2piFileName);
  printf("    >%s",realDataK2piFullPath);
  realDataK2piFP = fopen(realDataK2piFullPath,"w");

  asprintf(&ke3DataK2piFullPath,"%s%s",dir,ke3DataK2piFileName);
  printf("    >%s",ke3DataK2piFullPath);
  ke3DataK2piFP = fopen(ke3DataK2piFullPath,"w");

  asprintf(&km3DataK2piFullPath,"%s%s",dir,km3DataK2piFileName);
  printf("    >%s",km3DataK2piFullPath);
  km3DataK2piFP = fopen(km3DataK2piFullPath,"w");

  asprintf(&k2piDataK2piFullPath,"%s%s",dir,k2piDataK2piFileName);
  printf("    >%s",k2piDataK2piFullPath);
  k2piDataK2piFP = fopen(k2piDataK2piFullPath,"w");

  fprintf(realDataK2piFP,"%lf %lf\n",realDataNEvent,realDataNFailed);
  fprintf(ke3DataK2piFP,"%lf %lf\n", ke3DataNEvent, ke3DataNFailed);
  fprintf(km3DataK2piFP,"%lf %lf\n", km3DataNEvent, km3DataNFailed);
  fprintf(k2piDataK2piFP,"%lf %lf\n",k2piDataNEvent,k2piDataNFailed);

  fclose(realDataK2piFP);
  fclose(ke3DataK2piFP);
  fclose(km3DataK2piFP);
  fclose(k2piDataK2piFP);
  */
#endif
  /*----------- End of user C code -----------*/
  return 0;
}

