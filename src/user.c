/* File to define variables to be used in usersrc routines */
#include "user.h"
#include "reader.h"

FILE * fprt;
int muvAccept(float x, float y)
{
  if (fabs(x)<13 && fabs(y)<13) return 0;
  if (fabs(x)>130 || fabs(y)>130) return 0;
  return 1;
}
/* returns the dot product of two float 3-vectors */
float f3vdot(float *x,float *y)
{
  return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}

/* returns the square of a float 3-vector magnitude */
float f3vmag2(float *v)
{
  return ( _POWER2(v[0]) + _POWER2(v[1]) + _POWER2(v[2]) );
}

/* Returns the magnitude of a float 3-vector */
float f3vmag(float *v)
{
  return (float)sqrt( _POWER2(v[0]) + _POWER2(v[1]) + _POWER2(v[2]) );
}

/*Returns the (+---) dot product of two four vectors*/
float f4vdot(float * v1, float * v2)
{
  return +v1[0]*v2[0] -( v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3] );
}

int pointOfClosestApproach (float * point1, float * point2, float * v1, float * v2, float * dmin, float * vertex)
{

  //translate from track params to vector notation

  // printf("34 point1[0]: %f\n",point1[0]);
  int i;
  float r12[3];
  float det,q1,q2,t1,t2;
  float x;
  for(i=0;i<3;i++)
  {
    r12[i]=point1[i]-point2[i];
  }
  //printf("43 point1[0]: %f\n",point1[0]);

 /*
  12-08-98 IWS change pow to a simple product to avoid
               segmentation fault...
  det=pow(f3vdot(v1,v2),2)-f3vmag2(v1)*f3vmag2(v2);
  */
  x=f3vdot(v1,v2);
  det = _POWER2(x) - f3vmag2(v1)*f3vmag2(v2);
//  printf("Det: %f\n",det);
 // float  epsilon = 1e-9;
   if(det !=0.)
   {
      t1=(f3vmag2(v2)*f3vdot(r12,v1) - f3vdot(v1,v2) * f3vdot(r12,v2))/det;
      t2=(f3vdot(v1,v2)*f3vdot(r12,v1) - f3vmag2(v1) * f3vdot(r12,v2))/det;
      //printf("56 point1[0]: %f\n",point1[0]);
     for(i=0; i<3; i++)
     {
       q1 = point1[i] + (t1*v1[i]);
       q2 = point2[i] + (t2*v2[i]);
       vertex[i]=(q1+q2)/2.0;
       r12[i]=q1-q2;
     }
    *dmin = f3vmag(r12);

    return 0; //good
   }
   else
   {
    printf("ERROR: pointClosestApproach det too small\n");
	 return 1; //indicating parallel lines
   }
}


int mcEventType(superMcEvent * evt)
{
  /*
  if( (evt->Npart < 5) || (evt->part[0] != KAON) )
    {
      return UNKNOWN_EVENT;
    }
  */

  int pimiCount=0, pi0Count=0, piPlusCount=0, eCount=0, muCount=0, mumiCount=0, gammaCount=0,unknownCount=0,kCount=0;
  for(int i = 0;i<evt->Npart;++i)
    {
      switch(evt->part[i].type)
	{
	case 4:
	  ++pi0Count;
	  break;
	case 8:
	  ++piPlusCount;
	  break;
	case 16:
	  ++gammaCount;
	  break;
	case 32:
	  ++muCount;
	  break;
	case 64:
	  ++eCount;
	  break;
	case 512:
	  ++kCount;
	  break;
	case (-8):
	  ++pimiCount;
	  break;

	default:
	  ++unknownCount;
	  break;
	}
    }


  if( (pi0Count==1) && (eCount==1) && (gammaCount >=2)  &&
      (unknownCount==0) && (piPlusCount==0) && (muCount==0) )
    {
      return KE3_EVENT;
    }

  /****************************************/
  /*apparently the things I've decided are*/
  /*probably electrons randomly turn up in*/
  /*events in the file labeled kµ3        */
  /*                                      */
  /****************************************/
  if( (gammaCount>=2) && (pi0Count==1) && (muCount==1)
      /*&&(unknownCount==0) */&& (piPlusCount==0) /* && (eCount==0)*/ )
    {
      return KM3_EVENT;
    }
  if( (gammaCount>=2) && (pi0Count==1) && (piPlusCount==1)/* &&
      (unknownCount==0) && (eCount==0) && (muCount==0)*/ )
    {
      //printf("k2pi\n");
      return K2P_EVENT;
    }
if( (gammaCount == 1) && (muCount >=1) && (pi0Count==0) && (piPlusCount==0) )
    {
	return KM2_EVENT;
    }
if(  (gammaCount >=4) && (pi0Count==2) && (piPlusCount==1)   )
    {
	return K3P0_EVENT;
    }
if(  (piPlusCount == 2)  && ( pimiCount == 1) )
    {
	return K3P_EVENT;
    }

  printf("P0: %d\n",pi0Count);
  printf("P+: %d\n",piPlusCount);
  printf("P-: %d\n",pimiCount);
  printf("e+: %d\n",eCount);
  printf("m+: %d\n",muCount);
  printf("ga: %d\n",gammaCount);
  printf("k+: %d\n",kCount);
  printf("un: %d\n",unknownCount);
  printf("\n");

  for(int i = 0;i<evt->Npart;++i)
    {
      printf("Type: %d\n",evt->part[i].type);
    }
  printf("\n\n");
  //if we get here we are a bit clueless!

  return UNKNOWN_EVENT;

}

void GetCpdCellIndex(double pos_x, double pos_y, int *cpd_index, int *cell_index)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                              //
    // For a given x,y coordinate at the LKr front face, this routine returns the index of the CPD  //
    // and Cell at that position (official numbering). The information is needed e.g. for the       //
    // appropriate energy recalibration per cell and the electron ID correction for Ke2 MC which    //
    // can be different for single bad CPDs.                                                        //
    //                                                                                              //
    // To reduce computing time, the definition of the corner positions of CPDs and Cells has to be //
    // done before at the beginning of the job (e.g. in user_init.c) and then stored in memory.     //
    // An example how to do it can be found on the Ke2 analysis page, section Cell-by-cell E/p      //
    // recalibration.                                                                               //
    //                                                                                              //
    // A. Winhart  -  10.11.2009                                                                    //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////


    int   i, j, k;
    int   l, m, n;
    double CELLlength = 1.975;             // Cell length in cm
    double CPDlength = 8 * CELLlength;     // One CPD consists of 8x8 cells

    *cpd_index  = -1;
    *cell_index = -1;

    // Define CPD index
    for (i=0; i<16; i++)
    for (j=0; j<16; j++)
    {
        k = i*16 + j;
        if (pos_x <= CPDpos_leftDownCorner[k][0] && pos_x > (CPDpos_leftDownCorner[k][0]-CPDlength) &&
                pos_y >= CPDpos_leftDownCorner[k][1] && pos_y < (CPDpos_leftDownCorner[k][1]+CPDlength) )
        {
            //printf ("CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
            //  k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1], pos_x, pos_y);
            *cpd_index = k;
            break;
        }
    }

    // Define Cell index if CPD has been foundf
    if (*cpd_index >= 0)
    for (m=0; m<8; m++)
    for (n=0; n<8; n++)
    {
        l = m*8 + n;
        //printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n",
        //      l, *cpd_index, CELLpos_leftDownCorner[*cpd_index][l][0], CELLpos_leftDownCorner[*cpd_index][l][1]);

        if (pos_x <= CELLpos_leftDownCorner[*cpd_index][l][0] && pos_x > (CELLpos_leftDownCorner[*cpd_index][l][0]-CELLlength) &&
                pos_y >= CELLpos_leftDownCorner[*cpd_index][l][1] && pos_y < (CELLpos_leftDownCorner[*cpd_index][l][1]+CELLlength) )
        {
            //printf ("Cell %d in CPD %d: position left down corner = %.2f, \t%.2f \t track position = %.2f, \t %.2f\n",
            //  l, *cpd_index, CELLpos_leftDownCorner[*cpd_index][l][0], CELLpos_leftDownCorner[*cpd_index][l][1], pos_x, pos_y);
            *cell_index = l;
        }
    }
    //printf ("cpd_index = %3d \t cell_index = %3d\n", *cpd_index, *cell_index);

}


/*void blue_1trk(float vx[3], int iTrack, float *blue_hitbdxdz, float*blue_hitbdydz, superCmpEvent *sevt)
{
  /*
   * Routine to correct bluefield for single track
   * by G.L. (based on JBC Routine)
   * input: vx(3) vertex coordinate (0 in KS target)
   *        nn fortran inde of track
   * output: blue_hitbdxdz correct dx/dz before magnet
   *         blue_hitbdydz correct dy/dz before magnet
  */

 /* blue_1trk_(vx,&iTrack,blue_hitbdxdz,blue_hitbdydz,sevt,sevt);

}*/

//void closeap_(float * P1,float * P2,float * V1,float * V2,float* DMIN,float* VERTEX);
//void closeap_double_(double * P1,double * P2,double * V1,double * V2,double * DMIN,double* VERTEX);
