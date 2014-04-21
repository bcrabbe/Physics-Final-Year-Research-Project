/**********************************************************************/
/* C structure mapping the hodotime common block                      */
/**********************************************************************/
/*
#define NMAX_COUNTER 20
extern struct hodostruct{
  float hodotime;             Charged hod. event time 
  int itrace[2];              tracks index to be used by hodotime 
  float Pla[NMAX_COUNTER];    Planes indices used 
  float Cou[NMAX_COUNTER];    Counters indices used 
  float Tim[NMAX_COUNTER];    Time for counters used
  float sig[NMAX_COUNTER];   
  float hsigg;
  float vsigg;
}c_hodotime_;
*/

/**********************************************************************/
/* C structure mapping the tagtime common block                       */
/**********************************************************************/
#define MAXTAG 100
extern struct tagstruct{
  float prot_time[MAXTAG];  /* proton times */
  int prot_nhit[MAXTAG];    /* N. of tagger hits used to build the proton */
  int nprot;                /* Total number of protons */
  int prot_ihit[MAXTAG][MAXTAG]; /* First index of hit m of proton n [m][n] */
  int prot_jhit[MAXTAG][MAXTAG]; /* Snd index of hit m of proton n [m][n] */
  int prot_stat[MAXTAG];         /* proton code */
  int nprot_ladder;
  int nprot_monitor;
}timtag_;

/**********************************************************************/
/* C structure to hold the neutral analysis variables                 */
/**********************************************************************/
/*
struct ananeut{
  float esc;        Energy scale factor
  int iflag;        selection flag from routine user_sel2pi0 
  float cem[2];     2 pi0 masses 
  int inc[4];       clusters ids 
  float relli;      R ellipse 
  float ctau;       c*tau (lifetimes) 
  float ekaon;      Kaon energy 
}ananeut;
*/

