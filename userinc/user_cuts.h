#ifndef _usercuts_h
#define _usercuts_h

/* here all the cut values are defined */
/* set cut value to -999 if no cut     */
/* order of variables is order of cuts */


/* cuts for neutrals first .....*/ 

#define nngama       1.   /* 1 = exactly 4 gammas, 2 = 4 or 5 gammas */
#define nntrac       0.   /* maximum number of tracks allowed        */
#define nhacen       3.   /* maximum hac energy                      */
#define negesl      70.   /* minimum total energy in lkr             */
#define negesu     170.   /* maximum total energy in lkr             */
#define nbeamr      17.   /* minimum radius in lkr                   */
#define noktag     110.	  /* octagon cut for lkr                     */
#define negmin       3.   /* minimum photon energy                   */  
#define negmax     100.   /* maximum photon energy                   */  
#define ntimdi       5.   /* maximum time distance of showers        */     
#define ntimma    -999.   /* maximum time in lkr                     */
#define ndistm      20.   /* minimum shower distance                 */
#define nshowi     1.0    /* maximum shower width                    */
#define ncogra       8.   /* maximum center of gravity               */ 
#define nvertl       0.   /* neutral vertex lower bound              */ 
#define nvertu    5000.   /* neutral vertex upper bound              */ 
#define nellip       1.   /* ellipse number maximum                  */
#define ntimdt       2.   /* tagging window +-2ns                    */ 

/* now cuts for charged events ..... */
 
#define cntrac       2.   /* number of charged tracks                */
#define cpgesl      70.   /* minimum total momentum in dch           */
#define cpgesu     170.   /* maximum total momentum in dch           */
#define cptmin       3.   /* minimum track momentum                  */  
#define cptrat       3.   /* p ratio cut (against lambdas)           */ 
#define cptmax     100.   /* maximum track momentum                  */  
#define cepcut     0.85   /* E/p cut                                 */  
#define ccdacu       2.   /* cda cut                                 */
#define ccogks    -999.   /* maximum cog fuer ks                     */ 
#define ccogkl    -999.   /* maximum cog fuer kl                     */ 
#define ctimdt       2.   /* tagging window +-2ns                    */ 
#define cptks      0.0005 /* maximum pt**2 f"ur ks                   */ 
#define cptkl      0.0001 /* maximum pt**2 f"ur kl                   */ 


#endif


