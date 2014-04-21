#ifndef __userh__
#define __userh__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmpio.h>
#include <reader.h>
#include <float.h>
#include "compactF77.h"

#define _POWER2(x) ((x)*(x))

#define ENABLE_OUTPUT 1
#define MC_ENABLE_OUTPUT 1

extern FILE *fprt;
extern char gString[50];
extern char * dir;

extern FILE * FP1;
extern FILE * FP2;
extern FILE * cutWhichKilledEventFP;

#if ENABLE_OUTPUT
extern FILE * realFP;
#endif

#if MC_ENABLE_OUTPUT
extern FILE * ke3FP;
extern FILE * km3FP;
extern FILE * k2piFP;
extern FILE * km2FP;
extern FILE * k3pi0FP;
extern FILE * k3piFP;
#endif

extern int ke3Count, km3Count, k2piCount, numUnidentified, nEvents;
extern int MCke3Count, MCkm3Count, MCk2piCount, MCk3piCount, MCnumUnidentified, MCnEvents;



int muvAccept(float x, float y);

float f3vmag2(float *v);
float f3vmag(float *v);
float f3vdot(float * v1, float * v2);

float f4vdot(float * v1, float * f2);

int crossProd(float * v1, float * v2, float * vOut);
int pointOfClosestApproach(float * point1, float * point2, float * v1, float * v2, float * distance, float * vertex);
void GetCpdCellIndex(double pos_x, double pos_y, int *cpd_index, int *cell_index);
//Get Cpd declarations:
extern float CPDpos_leftDownCorner[256][2];       // CPDpos_leftDownCorner[256CPDs][x,y]
extern float CELLpos_leftDownCorner[256][64][2];  // CELLpos_leftDownCorner[256CPDs][64Cells][x,y]
//extern int   CPDindex, CELLindex;
//extern float CPDlength, CELLlength;
// E/p corrections for each cell
extern FILE *EopCorrfile;                // correction file to be read
extern float EopCorr[256][64];           // correction for each cell

//void blue_1trk_(float *vx, int *nn,float *blue_hitbdxdz,float *blue_hitbdydz,void *ise,void *qse);
//void blue_1trk(float vx[3], int iTrack, float *blue_hitbdxdz, float*blue_hitbdydz, superCmpEvent *sevt);
void closap_(float * P1,float * P2,float * V1,float * V2,float* DMIN,float* VERTEX);
void blue_tack_(int*nchar,float *tmom,float*Vxyz,float*vpnt,float*vdir);

#define BREAK_ON_FAILED_CUT 1

#define ENABLE_BASIC_QUALITY_CUTS 1
#define ENABLE_TIMING_CUTS 1
#define ENABLE_CRAZY_LKR_ACC_CUT 1
#define ENABLE_MIN_BIAS_CUT 1
#define ENABLE_Z_COORD_CUT 1
#define ENABLE_ENERGY_SHARING_CUTS 1
#define ENABLE_DCH_LKR_VETO 1
#define ENABLE_EXCLUDE_BAD_BURSTS 1


#define MC_ENABLE_BASIC_QUALITY_CUTS 1
#define MC_ENABLE_TIMING_CUTS 0 //the timings in the MC are fucked
#define MC_ENABLE_CRAZY_LKR_ACC_CUT 1
#define MC_ENABLE_MIN_BIAS_CUT 0 //apparently all of the MC data fails this
#define MC_ENABLE_Z_COORD_CUT 1
#define MC_ENABLE_ENERGY_SHARING_CUTS 1
#define MC_ENABLE_DCH_LKR_VETO 1
#define MC_ENABLE_EXCLUDE_BAD_BURSTS 1


#define SURVIVED 0
#define BASIC_QUALITY_CUTS 1
#define TIMING_CUTS 2
#define CRAZY_LKR_ACC_CUT 3
#define MIN_BIAS_CUT 4
#define Z_COORD_CUT 5

#define ELECTRON_MASS 0.0005109
#define MUON_MASS 0.1056583715
#define PI0_MASS 0.1349766
#define PI1_MASS 0.13957018 //this is the pi+ mass but you can't use operators in c names
                            //it has a charge of 1 so the name seems fairly sensible
#define UNKNOWN_EVENT -1
#define KE3_EVENT 0
#define KM3_EVENT 1
#define K2P_EVENT 2
#define KM2_EVENT 3
#define K3P0_EVENT 4
#define K3P_EVENT 5
//for particle ID:
#define ELECTRON 1
#define MUON 2
#define PIPLUS 3
#define UNIDENTIFIED 4

int mcEventType(superMcEvent * evt);
/*
 *
 *
 * : user.h,v $
 * Revision 1.2  2003/10/31 12:32:30  andrew
 * Added the -string option to main
 *
 * An arbitrary string (length STRING_MAX) can be passed to compact which is
 * saved in a global variable gString (C), COMMON/GSTRING/GSTRING (FORTRAN)
 *
 *
 * made -ndb the default. For people needin the compact database the -db option
 * was created
 *
 *
 *
 */

#endif //user.h
