/**************************************************************************/
/* Compact, OptiMized Program for Accessing Countless Terabytes           */
/*                                                                        */
/* COmPACT v4.3 - Main Program                                            */
/*                                                                        */
/* Modifications: 03-12-98 IWS try to read EoB when 2                     */
/*                         consecutives evts ahd same TS                  */
/*                         (pb found by G.U. + L.F.)                      */
/*                16-12-98 IWS uncomment mc2smc                           */ 
/*                05-01-99 IWS Dont call DB routine for MC events         */
/*                02-02-99 IWS from Guillaume: correct pb when            */
/*                             selecting SC from SC ({} missing)          */
/*                15-02-99 IWS trap EoB after double evts to              */
/*                             USER_EOB                                   */
/*                24-02-99 IWS option not to use CDB (-ndb)               */
/*                15-03-99 IWS call BuildEvtList when filtering           */
/*                         IWS Move test on BadB.Skip evt/evt             */
/*                             Do the same for super-compact              */
/*                24-03-99 IWS Call bad burst routine for SC              */
/*                10-05-99 IWS introduce cheat option                     */
/*                31-05-99 IWS more words into .compact-status file       */
/*                26-07-99 IWS 3 new arguments + call Ana routines for SC */
/*                05-08-99 IWS add sbur as argument to smc2float          */
/*                09-09-99 IWS call analysis routines for cheat option    */
/*                22-09-99 IWS test output file size (-b) for rfio        */
/*                29-09-99 IWS cmp2scmpBurst for option -SMO (M.Veltri)   */
/*                01-10-99 IWS option empty to empty bursts, EoB, lists   */
/*                             when filtering events                      */
/*                07-10-99 IWS update bursts counter when using -cheat    */
/*                21-10-99 IWS buffer for input events                    */
/**************************************************************************/

#include <string.h>
#include "cmpio.h"
#include "reader.h"
#include "cdbmap.h"
#define E_BUR 0x1
#define E_LIST 0x2
#define E_EOB 0x4
#define BUFFER_DEPTH 4
/* Define geometry structure + a pointer with a better name */
struct GeomCompact geom_;
struct GeomCompact * Geom = &geom_;
int main(int argc,char **argv) {
  int i,j;
  int userStatus;                           /* user routine return value */
  int maxBursts=0,readBursts=0;             /* burst limit and no. read */
  int readEvents=0;
  int ncmp=0,nke3=0,nkmu3=0,nmc=0,nscmp=0;  /* number of event type read */
  int nsmc=0;

  BurstID sBurst,eBurst,cBurst;              /* start, end and current */

  char SIFflag=FALSE,COflag=FALSE;           /* option flags */
  char KEOflag=FALSE,KMUOflag=FALSE;         /* option flags */
  char MOflag=FALSE,SOflag=FALSE;            /* SIF:SIngle File, CO:COmpact */
  char SMOflag=FALSE, ETAOflag=FALSE;        /* KEO:Ke3, KMUO:Kmu3, MO:MC */
  /* SO:SuperCOmPACT, SMO:SuperMC ETAO: Eta*/
  char SRWflag=FALSE;                        /* rewrite SuperCOmPACT INT event */
  char Bflag=FALSE;                          /* produce fixed size o/p files */
  char NDBflag=FALSE;                        /* Not to use CDB */
  char CHEATflag=FALSE;                      /* to metamorphose a compact event into SC */
  char EMPTYflag=FALSE;                      /* to empty bursts,... when filtering */
  long int Empty;                            /* Code to empty burst */
  char NEVTflag=FALSE;                       /* read NevtToRead */
  int  NevtToRead=0;                         /* Number of events to read in */  
  char RMINflag=FALSE;                       /* Only look at burst with Run>=RunMin */
  char RMAXflag=FALSE;                       /* Only look at burst with Run<=RunMax */
  int  RunMin, RunMax;                     
  long int Bsize;                            /* o/p file size (bytes) */
  int opfnum=1;                              /* number of o/p file */
  char fext[10];                             /* filename extension */

  char Xflag=FALSE;                           /* exec compact.pl flag */
  char xargstr[256]="runcompact ";            /* exec string for compact.pl */

  char COfname[81],KEOfname[81],MOfname[81];   /* names for I/O files */
  char KMUOfname[81],SOfname[81],SMOfname[81];
  char ETAOfname[81];

  char rundir[161],SIfname[256];               /* run directory and filename*/

  char buf[256];
  int LastFile=FALSE;                          /* last input file flag */
  int FLflag=FALSE;                            /* file list flag */
  Flist first;                                 /* first file in list */
  Flist *current=&first;                       /* current file in list */

  CMPstream *cmpin,*cmpout,*ke3out,*kmu3out;   /* COmPACT streams */
  CMPstream *mcout,*scmpout,*smcout;
  static Burst      bur;
  static superBurst sbur;
/* To buffer input events */
  static cmpEvent * pcmpevt[BUFFER_DEPTH];          /* pointers to buffered event structures */
  static cmpEvent cmpevt[BUFFER_DEPTH];             /* event structures for buffers */
  int ibuffer, ibuf, ibufIndex, ibufC;                     /* running variables */
  int bufIndex[BUFFER_DEPTH];                       /* position in buffer */
  static int NextStruct=0;                          /* To store NextStruct */
  static int EoBReached=0;                          /* EndofBurs reached */
  static int NextEvent=0;                           /* NUmber of events after the event */
  static ke3Event ke3evt;
  static kmu3Event kmu3evt;
  static mcEvent mcevt;

  static INTsuperCmpEvent supercmpevti;              /* integer versions of SuperCOmPACT */
  static INTsuperMcEvent supermcevti;                /* event structures - for i/o only  */

  static superCmpEvent supercmpevt;                  /* float versions of SuperCOmPACT */
  static superMcEvent supermcevt;                    /* event structures - for users   */

  static EndofBurst      eob;
  static superEndofBurst seob;

  int level=3;                               /* SuperCOmPACT level: 1/2/3 */

  int JumpToNextFile=0;                      /* Flag set when reading problem happens */
  int NotOpen=1;                             /* Problem when opening the file */
  bool_t zread;                              /* To check on rwRoutines */
  int *lpoint;                               /* pointer to zero Procerror + trig&tstamp list */
  gz_stream *s;                              /* pointer to o/p gz_stream */

  int numFiles=0;                            /* number of files in list */
  int iFile=0;                               /* index of current file */
  Flist *tmp=&first;
  FILE *fp;                                  /* pointer to .compact-status file */
  int statusUpd;
  int nByteOut;                              /* How many bytes output so far into CMPstream */
  int ComprFactor=4.0;                       /* compression factor */

  static int timestampPrev=0;                /* Keep timestamp of previous event */
  int LastEventOut=0;                        /* timestamp of last event written out */


#ifdef COMPACT_F77
#ifdef ALPHA
/* This call is required by the F77 runtime library on Alphas */
  for_rtl_init_(&argc,argv);         /* init F77 IO for alpha */
#endif   /* leave Alpha specific F77 */
#endif   /* leave F77 specific code */

  strcpy(rundir,".");
   

  for(i=0;i<argc;i++) {                /* loop over command line arguments */
    if(strcmp(argv[i],"-s") == 0)      /* search for starting burst */     
      sBurst=parseBurstID(argv[++i]);
    if(strcmp(argv[i],"-e") == 0)      /* search for end burst */
      eBurst=parseBurstID(argv[++i]);
    if(strcmp(argv[i],"-n") == 0)      /* search for max bursts */
      {  
        maxBursts=atoi(argv[++i]);
        eBurst.nrun=999999;            /* make sure eBurst doesn't limit */
        eBurst.nbur=999999;            /* the size of the run */
      }
    if(strcmp(argv[i],"-i") == 0)      /* search for input filename */
      {  
        strcpy(SIfname,argv[++i]);
        SIFflag=TRUE;                  /* set single input file flag */
      }
    if(strcmp(argv[i],"-d") == 0)      /* search for directory name */
      strcpy(rundir,argv[++i]);
    if(strcmp(argv[i],"-l") == 0)      /* search for file name list */
      {
        FLflag=readFlist(argv[++i],&first); /* read file list */
        /* count files */
        while (tmp != NULL) {
          numFiles++;
          tmp=tmp->next;
        }
        statusUpd=numFiles/10+1;
      }

    if(strcmp(argv[i],"-x") == 0)      /* execute compact.pl from compact */
      {
        strcat(xargstr,argv[++i]);
        Xflag=TRUE;
      }
     
    if(strcmp(argv[i],"-co") == 0)     /* search for output COmPACT filename */
      { 
        strcpy(COfname,argv[++i]);
        COflag=TRUE;                   /* set COmPACT output file flag */
      }
    if(strcmp(argv[i],"-ko") == 0)     /* search for output Ke3 filename */
      { 
        strcpy(KEOfname,argv[++i]);
        KEOflag=TRUE;                  /* set COmPACT output file flag */
      }
    if(strcmp(argv[i],"-kmo") == 0)    /* search for output Ke3 filename */
      { 
        strcpy(KMUOfname,argv[++i]);
        KMUOflag=TRUE;                 /* set COmPACT output file flag */
      }
    if(strcmp(argv[i],"-mo") == 0)     /* search for output mc filename */
      { 
        strcpy(MOfname,argv[++i]);
        MOflag=TRUE;                   /* set COmPACT output file flag */
      }
    if(strcmp(argv[i],"-so") == 0)     /* search for output SuperCOmPACT fname */
      { 
        strcpy(SOfname,argv[++i]);
        SOflag=TRUE;                   /* set COmPACT output file flag */
      }
    if(strcmp(argv[i],"-smo") == 0)    /* search for output SuperMC file name */
      { 
        strcpy(SMOfname,argv[++i]);
        SMOflag=TRUE;                  /* set COmPACT output file flag */
      }
    if(strcmp(argv[i],"-sl") == 0)     /* search for output SuperCOmPACT level */
      { 
        level=atoi(argv[++i]);
      }
    if(strcmp(argv[i],"-srw") ==0)     /* rewrite SuperCOmPACT variables */
      {
        SRWflag=TRUE;
      }
    if(strcmp(argv[i],"-ndb") ==0)     /* CDB not wanted */
      {
        NDBflag=TRUE;
      }
    if(strcmp(argv[i],"-cheat") ==0)   /* To metamorphose a compact evt in SC */
      {
        CHEATflag=TRUE;
      }
    if(strcmp(argv[i],"-empty") ==0)   /* To empty bursts,... when filtering (compact) */
      {
        EMPTYflag=TRUE;
        Empty=atoi(argv[++i]);
        if(Empty==0) Empty |= (E_BUR|E_LIST|E_EOB);
        sprintf(buf," =======> OK! Code to save disk space is %8x",Empty);
        cmpmsg(buf);
      }
    if(strcmp(argv[i],"-nevt") ==0)   /* Number of events to read in */
      {
        NEVTflag=TRUE;
        NevtToRead=atoi(argv[++i]);
        sprintf(buf," =======> OK! compact will stop after %5d events",NevtToRead);
        cmpmsg(buf);
      }
    if(strcmp(argv[i],"-rmin") ==0)   /* RunMin to look at */
      {
        RMINflag=TRUE;
        RunMin=atoi(argv[++i]);
        sprintf(buf," =======> OK! compact will look at Run>=%5d",RunMin);
        cmpmsg(buf);
      }
    if(strcmp(argv[i],"-rmax") ==0)   /* RunMin to look at */
      {
        RMAXflag=TRUE;
        RunMax=atoi(argv[++i]);
        sprintf(buf," =======> OK! compact will look at Run<=%5d",RunMax);
        cmpmsg(buf);
      }
    if((RMINflag) && (RMAXflag) && (RunMax<RunMin)) 
      {
        sprintf(buf," ======> Incompatible limits for run numbers - stop");
        cmperror(buf);
      }
    if (strcmp(argv[i],"-b") == 0)     /* generate big output files */
      {
        i++;
        if (strstr(argv[i],"KB"))
          {
            argv[i][(int)strlen(argv[i])-2]='\0';
            Bsize=1000*atoi(argv[i]);
          }
        else if (strstr(argv[i],"MB"))
          {
            argv[i][(int)strlen(argv[i])-2]='\0';
            Bsize=1000000*atoi(argv[i]);
          }
        else if (strstr(argv[i],"GB"))
          {
            argv[i][(int)strlen(argv[i])-2]='\0';
            Bsize=1000000000*atoi(argv[i]);
          }
        else
          {	     
            Bsize=atoi(argv[i]);
          }
        Bflag=TRUE;                    /* set output file size */
      }
    if(strcmp(argv[i],"-h")==0 || argc==1) 
      {
        printHelp();
        exit(0);
      }
  }
   
  /* consistency check for -x option */
  if (Xflag) {
    if (FLflag||SIFflag) cmperror("Please don't use -l or -i with -x option");
    if(system(xargstr)) cmperror("Failed to exec runcompact");
    FLflag=readFlist("compact.list",&first); /* read file list */
    /* count files */
    while (tmp != NULL) {
      numFiles++;
      tmp=tmp->next;
    }
    statusUpd=numFiles/10+1;
  }

  /* consistency check for -b option */
  if (Bflag)
    {
      if (!(COflag||KEOflag||KMUOflag||MOflag||SOflag||SMOflag))
        {
          cmpmsg("Use -co/-ko/-kmo/-mo/-so/-smo with -b option");
          printHelp();
          exit(0);
        }

    }
  /* Database initialisation routine */
  if(!NDBflag) getRecordCmpInit();
   
  NotOpen=1;                         /* set flag before opening file */
  if(FLflag)                         /* reading from input file list */
    {
      cmpin=openCompact(current->fname,"rb"); /* open current file */
      iFile++;
    }
  else if(SIFflag)                   /* reading from a single file */
    {                     
      cmpin=openCompact(SIfname,"rb"); /* open COmPACT file */
    } 
  /*
    23-06-98 IWS : dont understand this 
    else 
    {
    cmpmsg("Calling openNextBurst...");
    cBurst=sBurst;                                 set current burst to start burst 
    cmpin=openNextBurst(rundir,&cBurst,eBurst);    open next burst 
    }
  */
   
   
  if(cmpin==(CMPstream *)0)
    cmperror("Input file not found");
   
   
  if(COflag)                              /* writing an output file */
    { 
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(COfname,fext);
        }
      cmpout=openCompact(COfname,"wb");   /* open COmPACT file */
      if(cmpout == (CMPstream *)0)
        cmperror("compact: Failed to open COmPACT output stream");
    }
   
  if(KEOflag)                             /* writing an output file */
    {
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(KEOfname,fext);
        }
      ke3out=openCompact(KEOfname,"kwb"); /* open COmPACT Ke3 file */
      if(ke3out == (CMPstream *)0)
        cmperror("compact: Failed to open Ke3 COmPACT output stream");
    }
   
  if(KMUOflag)                            /* writing an output file */
    {
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(KMUOfname,fext);
        }
      kmu3out=openCompact(KMUOfname,"lwb");/* open COmPACT Kmu3 file */
      if(kmu3out == (CMPstream *)0)
        cmperror("compact: Failed to open Kmu3 COmPACT output stream");
    }
   
  if(MOflag)                              /* writing an output file */
    {                   
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(MOfname,fext);
        }
      mcout=openCompact(MOfname,"mwb");   /* open COmPACT MC file */
      if(mcout == (CMPstream *)0)
        cmperror("compact: Failed to open MC COmPACT output stream");
    }
   
  if(SOflag)                             /* writing an output file */
    {             
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(SOfname,fext);
        }
      scmpout=openCompact(SOfname,"swb"); /* open SuperCOmPACT file */
      if(scmpout == (CMPstream *)0)
        cmperror("compact: Failed to open SuperCompact output stream");
    }
   
  if(SMOflag)                             /* writing an output file */
    {             
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(SMOfname,fext);
        }
      smcout=openCompact(SMOfname,"nwb"); /* open SuperCOmPACT file */
      if(smcout == (CMPstream *)0)
        cmperror("compact: Failed to open SuperMC output stream");
    }
  
  cmpAnaInit();                          /* initialise analysis routines */
  superCmpAnaInit();                     /* initialise analysis routines */
  userStatus=USER_INIT;                  /* call user initialization */
  fflush(NULL);
   
  if(userStatus>0) 
    {
      sprintf(buf,"Error code %d returned from user_init",userStatus);
      cmpwarn(buf);
    }
   
   
/********************* Enter file loop **********************/
  while(!LastFile) 
    {
/*+++++++++++++++++++ Enter struct loop ++++++++++++++++++++*/
      while(((cmpin->nextStruct != FTLR_STRUCT) && (!JumpToNextFile)) || 
            ((cmpin->nextStruct == FTLR_STRUCT) && (NextEvent>0)))  
        {
          if((NEVTflag) && (readEvents>=NevtToRead)) 
            {
              LastFile=TRUE;
              break; /* Break Loop */
            }
          /* 
             Until EoB is reached, use current nextStruct
          */
          if(!EoBReached) NextStruct=cmpin->nextStruct; /* Until EoB read next event */

          /* 
             When meet EoB, read it, and then go on with the buffered events.
          */
          if((NextStruct == EOB_STRUCT) && (NextEvent==BUFFER_DEPTH))
            {
              zread=rwEndofBurst(cmpin,&eob);     /* read in end of burst to get lists */
              cmpAnaEoB(&bur,&eob);               /* Build lists */
              EoBReached=1;
              NextStruct=CMPEVENT_STRUCT;
            }

          /* 
             If EoB reached and last evt anal. => ana. EoB 
             This includes the buffers with depth==1;
          */
          if((EoBReached && (NextEvent==1))) 
            {
              NextStruct=EOB_STRUCT;
              NextEvent=0;
            }
          
          switch(NextStruct) 
            {
         
	
/*----------> Burst structure <----------*/

              case BURST_STRUCT:            /* burst structure */
                zread=rwBurst(cmpin,&bur);  /* read in burst header */
                if(!zread)                  /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                /* 21-10-99 - setup the buffering */
                NextEvent=0;
                EoBReached=0;
                ibuffer=0;
                for(ibuf=0; ibuf<BUFFER_DEPTH; ibuf++)
                  {
                    ibufIndex=(ibuffer+ibuf)%BUFFER_DEPTH; /* position in buffer */
                    bufIndex[ibuf]=ibufIndex;              /* store position in buffer */
                    pcmpevt[ibuf]=&cmpevt[ibufIndex];      /* store adress in buffer */
                  }
                if( (RMINflag) && (bur.nrun < RunMin) ) break;
                if( (RMAXflag) && (bur.nrun > RunMax) ) break;
/* 05-01-99 IWS Dont call DB routine for MC events */
                if ((bur.brtype == DATA_BURST) && (!NDBflag) )
                  {
                    bur.dberr=getRecordCmp(bur.nrun, bur.time);
                    /* 
                       25-01-99 IWS set the return code from cdb to a positive value 
                       so that it is not confused with fuser_database_decide return code
                    */
                    if(bur.dberr<0) bur.dberr = -1*bur.dberr;
                  }
                else
                  {
                    /* 10-05-99 IWS correct floowing line: it was sbur.dberr */
                    bur.dberr=1; /* set flag as if database had not data, when DB not wanted */
                  }

                cmpAnaPreBurst(&bur);                      /* Compact analysis for burst before user routine */
                readBursts++;                              /* increase number of bursts read */
#ifdef COMPACT_F77
                if(readBursts==1)                /* run check on first event read */
                  {         
                    cmpmsg("main: First Burst found - running C/F77 checks");
                    chkF77Burst(&bur);
                  }
#endif
                if (SOflag || SMOflag) cmp2scmpInitBurst(&bur);
                userStatus=USER_BURST(&bur); /* call user burst routine */
                cmpAnaPostBurst(&bur);       /* Compact analysis for burst after user routine */
                fflush(NULL);
             
/* 10-05-99 IWS introduce cheat option */
                if(CHEATflag)
                  {
                    cmp2scmpBurst(&bur,&sbur);  /* copy superBurst variables */
                    if ((sbur.brtype == DATA_BURST) && (!NDBflag) )
                      {
                        sbur.dberr=getRecordCmp(sbur.nrun, sbur.time);
                        if(sbur.dberr<0) sbur.dberr = -1*sbur.dberr;
                      }
                    else
                      {
                        sbur.dberr=1; /* set flag as if database had not data, when DB not wanted */
                      }
                    superCmpAnaPreBurst(&sbur);
                    userStatus=USER_SUPERBURST(&sbur);  /* call user burst routine */
                    superCmpAnaPostBurst(&sbur);
                    break;
                  } /* End test on CHEAT flag */

                if(userStatus>0) 
                  {
                    sprintf(buf,"Error code %d returned from user_burst",userStatus);
                    cmpwarn(buf);
                  }
                if((EMPTYflag) && (Empty&E_BUR))
                  {
                    EmptyBurst(&bur,100);
                  }
                if(COflag)
                  rwBurst(cmpout,&bur);
             
                if(KEOflag)
                  rwBurst(ke3out,&bur);
             
                if(KMUOflag)
                  rwBurst(kmu3out,&bur);
             
                if(MOflag)
                  rwBurst(mcout,&bur);
             
                if (SOflag) {
                  cmp2scmpBurst(&bur,&sbur);  /* copy superBurst variables */
                  rwSuperBurst(scmpout,&sbur);
                }
             
                if (SMOflag) {
                  cmp2scmpBurst(&bur,&sbur);  /* copy superBurst variables */
                  rwSuperBurst(smcout,&sbur);
                }
             
                break;
           

/*----------> SuperCOmPACT burst structure <----------*/

              case SBURST_STRUCT:                   /* SuperCOmPACT burst structure */
                zread=rwSuperBurst(cmpin,&sbur);    /* read in burst header */
                if(!zread)                          /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if( (RMINflag) && (sbur.nrun < RunMin) ) break;
                if( (RMAXflag) && (sbur.nrun > RunMax) ) break;

/* 05-01-99 IWS Dont call DB routine for MC events */
/* 05-02-99 IWS dont put in this protection yet
/*
  if ( (sbur.brtype == DATA_BURST) &&&& (!NDBflag)  )
  {
  }
*/
                if(!NDBflag) 
                  {
                    sbur.dberr=getRecordCmp(sbur.nrun, sbur.time);
                    /* 
                       25-01-99 IWS set the return code from cdb to a positive value 
                       so that it is not confused with fuser_database_decide return code
                    */
                    if(sbur.dberr<0) sbur.dberr = -1*sbur.dberr;
                  }
                else
                  {
                    sbur.dberr=1; /* set flag as if database had not data, when DB not wanted */
                  }
                superCmpAnaPreBurst(&sbur);             /* Compact analysis for burst before user routine */
                userStatus=USER_SUPERBURST(&sbur);  /* call user burst routine */
                superCmpAnaPostBurst(&sbur);        /* Compact analysis for burst before user routine */
                fflush(NULL);
                readBursts++;                       /* increase number of bursts read */
                 
                if(userStatus>0) 
                  {
                    sprintf(buf,"Error code %d returned from user_sburst",userStatus);
                    cmpwarn(buf);
                  }

                if (SOflag) {
                  rwSuperBurst(scmpout,&sbur);      /* write out superBurst */
                }
             
                if (SMOflag) {
                  rwSuperBurst(smcout,&sbur);
                }
           
                if (Bflag)
                  {
                    if (SOflag)
                      {
                        nByteOut=(int) ( (float) scmpout->bc/ComprFactor);
                        if (nByteOut>Bsize)   /* check COmPACT output file size */
                          {
                            closeCompact(scmpout);     /* close current file */
                            opfnum++;                 /* increment output file number */		    
                            i=(int)strlen(SOfname)-(int)strlen(fext);
                            SOfname[i]='\0';          /* strip off old file extension */
                            sprintf(fext,".%d",opfnum);
                            strcat(SOfname,fext);
                            scmpout=openCompact(SOfname,"swb"); /* open new file */
                          }
                      }
                    if (SMOflag)
                      {
                        nByteOut=(int) ( (float) smcout->bc/ComprFactor);
                        if (nByteOut>Bsize)   /* check COmPACT output file size */
                          {
                            closeCompact(smcout);  /* close current file */
                            opfnum++;              /* increment output file number */		    
                            i=(int)strlen(SMOfname)-(int)strlen(fext);
                            SMOfname[i]='\0';       /* strip off old file extension */
                            sprintf(fext,".%d",opfnum);
                            strcat(SMOfname,fext);
                            smcout=openCompact(SMOfname,"nwb"); /* open new file */
                          }
                      }
                  }
                break;
           
           
/*----------> Compact event structure <----------*/
           
              case CMPEVENT_STRUCT:                   /* normal COmPACT event */
                /* Read-in BUFFER_DEPTH events at begining of burst */
                ibufC=bufIndex[0];
                if(!NextEvent)
                  {
                    for(ibuf=0; ibuf<BUFFER_DEPTH; ibuf++)
                      {
                        if(cmpin->nextStruct == CMPEVENT_STRUCT)
                          {
                            zread=rwCmpEvent(cmpin,pcmpevt[ibuf],&bur); /* read event */
                            cmpAnaEvent(&bur,pcmpevt[ibuf]);            /* Compact analysis / event */
                          }
                        if( cmpin->nextStruct == EOB_STRUCT)                             
                          {
                            zread=rwEndofBurst(cmpin,&eob);     /* try to read in end of burst */
                            if(zread)                           /* End of Burst reached */
                              {
                                EoBReached=1;
                                printf(" breaking after only %1d\n", NextEvent);
                                break;
                              }
                            else
                              {
                                JumpToNextFile=1;              /* Problem with input file */
                                break;
                              }
                          }
                        else
                          {
                            NextEvent++;
                          }
                      }
                  }
                else if(!EoBReached)
                  {
                    ibuf=BUFFER_DEPTH-1;              /* Find out position in buffer for next evt */
                    zread=rwCmpEvent(cmpin,pcmpevt[ibuf],&bur); /* read event */
                    cmpAnaEvent(&bur,pcmpevt[ibuf]);            /* Compact analysis / event */
                  }
                else
                  {
                    NextEvent--;
                  }
                if(bur.BadB.Skip) break; /* Dont analyse events from bad bursts */
                if( (RMINflag) && (bur.nrun < RunMin) ) break;
                if( (RMAXflag) && (bur.nrun > RunMax) ) break;
                 
                getPrevNext(&bur,pcmpevt,bufIndex,NextEvent);
#ifdef COMPACT_F77
                if(ncmp==0)                           /* run check on first event read */
                  {         
                    cmpmsg("main: First cmpEvent found - running C/F77 checks");
                    chkF77CmpEvent(pcmpevt[0]);
                  }
#endif
/* 10-05-99 IWS introduce cheat option */
                if(CHEATflag)
                  {
                    cmp2scmp(&bur,pcmpevt[0],&supercmpevti,level); /* convert event */
                    scmp2float(&supercmpevti,&supercmpevt,&sbur);      /* convert user variables to floats */
                    calcScmpEvent(&supercmpevt);
                    superCmpAnaEvent(&sbur,&supercmpevt); /* Analysis routines for super-compact */
                    userStatus=USER_SUPERCMPEVENT(&sbur,&supercmpevt); /* <-- float version */
                    ncmp++;                     /* increment number of COmPACT events */
                    readEvents++;
                    /* 21-10-99 - setup the buffering */
                    ibuffer++;
                    for(ibuf=0; ibuf<BUFFER_DEPTH; ibuf++)
                      {
                        ibufIndex=(ibuffer+ibuf)%BUFFER_DEPTH; /* position in buffer */
                        bufIndex[ibuf]=ibufIndex;              /* store position in buffer */
                        pcmpevt[ibuf]=&cmpevt[ibufIndex];      /* store adress in buffer */
                      }
                    break;
                  } /* End test on CHEAT flag */
                userStatus=USER_CMPEVENT(&bur,pcmpevt[0]);
                fflush(NULL);
           
                if(userStatus>0)              /* check user return for error */
                  {      
                    sprintf(buf,"Error code %d returned from user_cmpEvent",
                            userStatus);
                    cmpwarn(buf);
                  }
             
                if(COflag)                    /* is the COmPACT output flag set? */
                  {               
                    userStatus=USER_CMPFILTER(&bur,pcmpevt[0]);
                    fflush(NULL);
                 
                    if(userStatus>0)         /* examine filter return for error */
                      {
                        sprintf(buf,"Error code %d returned from user_cmpFilter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)         /* check if event has passed filter */
                      {
                        if((EMPTYflag)&&(Empty&E_LIST))        /* Dont write out event lists */
                          {
                            pcmpevt[0]->NTrigWordL       = 0;
                            pcmpevt[0]->NTimeStampL      = 0;
                            pcmpevt[0]->NDCHDecErrorL    = 0;
                            pcmpevt[0]->NLKRHACDecErrorL = 0;
                            pcmpevt[0]->NChargedInfoL    = 0;
                            pcmpevt[0]->NNeutralInfoL    = 0;
                          }
                        else
                          {
                            BuildEvtList(&bur, pcmpevt[0], LastEventOut);
                            LastEventOut=pcmpevt[0]->timestamp;
                          }
                        rwCmpEvent(cmpout,pcmpevt[0],&bur); /* write out good event */
                      }
                  }
             
                if(SOflag)                  /* is the SuperCOmPACT output flag set? */
                  {               
                    userStatus=USER_CMPFILTER(&bur,pcmpevt[0]);
                    fflush(NULL);
                 
                    if(userStatus>0)        /* examine filter return for error */
                      {
                        sprintf(buf,"Error code %d returned from user_cmpFilter",
                                userStatus);
                        cmpwarn(buf);
                      }
                 
                    if(userStatus<0) {      /* check if event has passed filter */
                      cmp2scmp(&bur,pcmpevt[0],&supercmpevti,level); /* convert event */
                      rwSuperCmpEvent(scmpout,&supercmpevti);     /* write out event */
                    }
                  }
             
                ncmp++;                     /* increment number of COmPACT events */
                readEvents++;
                /* 21-10-99 - setup the buffering */
                ibuffer++;
                for(ibuf=0; ibuf<BUFFER_DEPTH; ibuf++)
                  {
                    ibufIndex=(ibuffer+ibuf)%BUFFER_DEPTH; /* position in buffer */
                    bufIndex[ibuf]=ibufIndex;              /* store position in buffer */
                    pcmpevt[ibuf]=&cmpevt[ibufIndex];      /* store adress in buffer */
                  }
                break;
           
           
           
/*----------> Monte Carlo event structure <----------*/
           
              case MCEVENT_STRUCT:                    /* MC event */
                zread=rwMcEvent(cmpin,&mcevt,&bur);   /* read event */
                if(!zread)                            /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if( (RMINflag) && (bur.nrun < RunMin) ) break;
                if( (RMAXflag) && (bur.nrun > RunMax) ) break;
                cmpAnaEvent(&bur,&mcevt.cmpevt);
#ifdef COMPACT_F77
                if(nmc==0) {                           /* run check on first event read */
                  cmpmsg("main: First mcEvent found - running C/F77 checks");
                  chkF77McEvent(&cmpevt);
                }
#endif
                userStatus=USER_MCEVENT(&bur,&mcevt);
                fflush(NULL);
             
                if(userStatus>0)                      /* check user return for error */
                  {
                    sprintf(buf,"Error code %d returned from user_McEvent",
                            userStatus);
                    cmpwarn(buf);
                  }
             
                if(MOflag)                            /* is the COmPACT output flag set? */
                  {         
                    userStatus=USER_MCFILTER(&bur,&mcevt);
                    fflush(NULL);
                 
                    if(userStatus>0)                  /* examine filter return for error */
                      { 
                        sprintf(buf,"Error code %d returned from user_mcFilter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)                 /* check if event has passed filter */
                      rwMcEvent(mcout,&mcevt,&bur);  /* write out good event */
                  }
             
                if(SMOflag)                          /* is the SuperMC output flag set? */
                  {               
                    userStatus=USER_MCFILTER(&bur,&mcevt);
                    fflush(NULL);
                 
                    if(userStatus>0)                /* examine filter return for error */
                      {
                        sprintf(buf,"Error code %d returned from user_mcFilter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)               /* check if event has passed filter */
                      {
                        mc2smc(&bur,&mcevt,&supermcevti,level);       /* convert event */
                        rwSuperMcEvent(smcout,&supermcevti);/* write out good event */
                      }
                  }
             
                nmc++;                             /* increment num of MC events read */
                readEvents++;
                break;
             
           
           
/*----------> Ke3 Compact event structure <----------*/
           
              case KE3EVENT_STRUCT:               /* Ke3 event */
                zread=rwKe3Event(cmpin,&ke3evt);  /* read event */
                if(!zread)                        /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if( (RMINflag) && (bur.nrun < RunMin) ) break;
                if( (RMAXflag) && (bur.nrun > RunMax) ) break;
#ifdef COMPACT_F77
                if(nke3==0) {                     /* run check on first event read */
                  cmpmsg("main: First ke3Event found - running C/F77 checks");
                  chkF77Ke3Event(&ke3evt);
                }
#endif
                userStatus=USER_KE3EVENT(&bur,&ke3evt);
                fflush(NULL);
             
                if(userStatus>0)                  /* check user return for error */
                  {
                    sprintf(buf,"Error code %d returned from user_ke3Event",
                            userStatus);
                    cmpwarn(buf);
                  }
             
                if(KEOflag)                      /* is the COmPACT output flag set? */
                  {         
                    userStatus=USER_KE3FILTER(&bur,&ke3evt);
                    fflush(NULL);
                 
                    if(userStatus>0)            /* examine filter return for error */
                      { 
                        sprintf(buf,"Error code %d returned from user_ke3Filter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)           /* check if event has passed filter */
                      rwKe3Event(ke3out,&ke3evt);/* write out good event */
                  }
                nke3++;                        /* increment num of Ke3 events read */
                readEvents++;
                break;
             
             
/*----------> Kmu3 Compact event structure <----------*/
           
              case KMU3EVENT_STRUCT:                    /* Kmu3 event */
                zread=rwKmu3Event(cmpin,&kmu3evt,&bur); /* read event */
                if(!zread)                              /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if( (RMINflag) && (bur.nrun < RunMin) ) break;
                if( (RMAXflag) && (bur.nrun > RunMax) ) break;
                cmpAnaEvent(&bur,&mcevt.cmpevt);
#ifdef COMPACT_F77
                if(nkmu3==0) {                            /* run check on first event read */
                  cmpmsg("main: First kmu3Event found - running C/F77 checks");
                  /* chkF77Kmu3Event(&cmpevt); */
                }
#endif
                userStatus=USER_KMU3EVENT(&bur,&kmu3evt);
                fflush(NULL);
             
                if(userStatus>0)                         /* check user return for error */
                  {
                    sprintf(buf,"Error code %d returned from user_kmu3Event",
                            userStatus);
                    cmpwarn(buf);
                  }
             
                if(KMUOflag)                            /* is the COmPACT output flag set? */
                  {         
                    userStatus=USER_KMU3FILTER(&bur,&kmu3evt);
                    fflush(NULL);
                 
                    if(userStatus>0)                    /* examine filter return for error */
                      { 
                        sprintf(buf,"Error code %d returned from user_kmu3Filter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)                    /* check if event has passed filter */
                      rwKmu3Event(kmu3out,&kmu3evt,&bur);/* write out good event */
                  }
                nkmu3++;                                 /* increment num of Ke3 events read */
                readEvents++;
                break;
           
           

/*----------> SuperCompact event structure <----------*/

              case SCMPEVENT_STRUCT:                        /* SuperCOmPACT event */
                zread=rwSuperCmpEvent(cmpin,&supercmpevti);     /* read event */
                if(!zread)                                  /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if(sbur.BadB.Skip) break; /* Dont analyse events from bad bursts 16-03-99 */
                if( (RMINflag) && (sbur.nrun < RunMin) ) break;
                if( (RMAXflag) && (sbur.nrun > RunMax) ) break;
                scmp2float(&supercmpevti,&supercmpevt,&sbur);      /* convert user variables to floats */
                calcScmpEvent(&supercmpevt);
#ifdef COMPACT_F77
                if(nscmp==0)                                /* run check on first event read */
                  {         
                    cmpmsg("main: First superCmpEvent found - running C/F77 checks");
                    chkF77SuperCmpEvent(&supercmpevt);
                  }
#endif
                superCmpAnaEvent(&sbur,&supercmpevt); /* Analysis routines for super-compact */
                userStatus=USER_SUPERCMPEVENT(&sbur,&supercmpevt); /* <-- float version */
                fflush(NULL);
           
                if(userStatus>0)                            /* check user return for error */
                  {      
                    sprintf(buf,"Error code %d returned from user_superCmpEvent",
                            userStatus);
                    cmpwarn(buf);
                  }
             
                if(COflag)                                 /* is the COmPACT output flag set? */
                  {
                    sprintf(buf,"Cannot write COmPACT o/p from SuperCOmPACT i/p file");
                    cmpwarn(buf);
                  }
             
                if(SOflag)                                 /* is the SuperCOmPACT output flag set? */
                  {               
                    userStatus=USER_SUPERCMPFILTER(&sbur,&supercmpevt);  /* <-- float version */
                    fflush(NULL);
                 
                    if(userStatus>0)                       /* examine filter return for error */
                      {
                        sprintf(buf,"Error code %d returned from user_superCmpFilter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)                       /* check if event has passed filter */
                      {
                        if (SRWflag) scmp2scmp(&supercmpevt,&supercmpevti);
                        reduceScmpLevel(level,&supercmpevti);
                        rwSuperCmpEvent(scmpout,&supercmpevti);/* write out good event */
                      }
                  }
             
                nscmp++;                                   /* increment number of SuperCOmPACT events */
                readEvents++;
                break;
             
           
           
/*----------> SuperMC event structure <----------*/

              case SMCEVENT_STRUCT:                        /* SuperMC event */
                zread=rwSuperMcEvent(cmpin,&supermcevti);  /* read event */
                if(!zread)                                 /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if( (RMINflag) && (sbur.nrun < RunMin) ) break;
                if( (RMAXflag) && (sbur.nrun > RunMax) ) break;
                smc2float(&supermcevti,&supermcevt,&sbur);   /* convert user variables to floats */
#ifdef COMPACT_F77
                if(nsmc==0)                               /* run check on first event read */
                  {         
                    cmpmsg("main: First superMcEvent found - running C/F77 checks");
                    chkF77SuperMcEvent(&supermcevt);
                  }
#endif
                userStatus=USER_SUPERMCEVENT(&sbur,&supermcevt); /* <-- float version */
                fflush(NULL);
             
                if(userStatus>0)                           /* check user return for error */
                  {      
                    sprintf(buf,"Error code %d returned from user_superMcEvent",
                            userStatus);
                    cmpwarn(buf);
                  }
             
                if(MOflag)                                /* is the MC output flag set? */
                  {
                    sprintf(buf,"Cannot write MC o/p from SuperMC i/p file");
                    cmpwarn(buf);
                  }
             
                if(SMOflag)                               /* is the SuperMC output flag set? */
                  {               
                    userStatus=USER_SUPERMCFILTER(&sbur,&supermcevt); /* <-- float version */
                    fflush(NULL);
                 
                    if(userStatus>0)                       /* examine filter return for error */
                      {
                        sprintf(buf,"Error code %d returned from user_superMcFilter",
                                userStatus);
                        cmpwarn(buf);
                      }
                    if(userStatus<0)                      /* check if event has passed filter */
                      rwSuperMcEvent(smcout,&supermcevti);/* write out good event */
                  }
             
                nsmc++;                                   /* increment number of SuperMC events */
                readEvents++;
                break;
             
             
             
/*----------> End-of-burst structure <----------*/

              case EOB_STRUCT:                      /* End-of-burst structure */
                if(!EoBReached)
                  {
                    zread=rwEndofBurst(cmpin,&eob);     /* read in end of burst */
                    if(!zread)                          /* skip burst if problem reading */
                      {
                        JumpToNextFile=1; 
                        break;
                      }
                    cmpAnaEoB(&bur,&eob);
                  }
                EoBReached=0;
                if( (RMINflag) && (bur.nrun < RunMin) ) break;
                if( (RMAXflag) && (bur.nrun > RunMax) ) break;


                userStatus=USER_EOB(&bur,&eob);
                fflush(NULL);
             
                if(userStatus>0)              /* examine filter return for error */
                  {
                    sprintf(buf,"Error code %d returned from user_eob",userStatus);
                    cmpwarn(buf);
                  }
                if(COflag)                    /* check for COmPACT output file */
                  {
                    if((EMPTYflag))
                      {
                        if(Empty&E_LIST)
                          {
                            bur.NTrigWord      = eob.NTrigWord      = 0;
                            bur.NTimeStamp     = eob.NTimeStamp     = 0;
                            bur.NDCHDecError   = eob.NDCHDecError   = 0;
                            bur.NLKRHACDecError= eob.NLKRHACDecError= 0;
                            bur.NChargedInfo   = eob.NChargedInfo   = 0;
                            bur.NNeutralInfo   = eob.NNeutralInfo   = 0;
                          }
                        if(Empty&E_EOB)
                          {
                            /* Empty EndofBurst structure */
                            memset(&eob,0,sizeof(EndofBurst));
                          }
                      }
                        /* Update event list for EoB structure */
                    BuildEoBList(&bur, &eob, LastEventOut);
                    rwEndofBurst(cmpout,&eob);  /* and write out end of burst struct */
                  }
                if(KEOflag)                   /* check for Ke3 output file */
                  rwEndofBurst(ke3out,&eob);  /* and write out end of burst struct */
                if(KMUOflag)                  /* check for Kmu3 output file */
                  rwEndofBurst(kmu3out,&eob); /* and write out end of burst struct */
                if(MOflag)                    /* check for Mc output file */
                  rwEndofBurst(mcout,&eob);   /* and write out end of burst struct */
             

                /* if fixed size output files are required check size now */
                if (Bflag)
                  {
                    if (COflag)
                      {
                        nByteOut=(int) ( (float) cmpout->bc/ComprFactor);
                        if (nByteOut>Bsize)   /* check COmPACT output file size */
                          {
                            closeCompact(cmpout);  /* close current file */
                            opfnum++;              /* increment output file number */		    
                            i=(int)strlen(COfname)-(int)strlen(fext);
                            COfname[i]='\0';       /* strip off old file extension */
                            sprintf(fext,".%d",opfnum);
                            strcat(COfname,fext);
                            cmpout=openCompact(COfname,"wb"); /* open new file */
                          }
                      }
                    if (KEOflag)
                      {
                        s=(gz_stream*)(ke3out->file);
                        if (ftell((FILE*)(s->file))>Bsize)   /* check Ke3 output file size */
                          {
                            closeCompact(ke3out);  /* close current file */
                            opfnum++;              /* increment output file number */		    
                            i=(int)strlen(KEOfname)-(int)strlen(fext);
                            KEOfname[i]='\0';      /* strip off old file extension */
                            sprintf(fext,".%d",opfnum);
                            strcat(KEOfname,fext);
                            ke3out=openCompact(KEOfname,"kwb"); /* open new file */
                          }
                      }
                    if (KMUOflag)
                      {
                        s=(gz_stream*)(kmu3out->file);
                        if (ftell((FILE*)(s->file))>Bsize)   /* check Kmu3 output file size */
                          {
                            closeCompact(kmu3out); /* close current file */
                            opfnum++;              /* increment output file number */		    
                            i=(int)strlen(KMUOfname)-(int)strlen(fext);
                            KMUOfname[i]='\0';     /* strip off old file extension */
                            sprintf(fext,".%d",opfnum);
                            strcat(KMUOfname,fext);
                            kmu3out=openCompact(KMUOfname,"lwb"); /* open new file */
                          }
                      }
                    if (MOflag)
                      {
                        s=(gz_stream*)(mcout->file);
                        if (ftell((FILE*)(s->file))>Bsize)       /* check MC output file size */
                          {
                            closeCompact(mcout);   /* close current file */
                            opfnum++;              /* increment output file number */		    
                            i=(int)strlen(MOfname)-(int)strlen(fext);
                            MOfname[i]='\0';       /* strip off old file extension */
                            sprintf(fext,".%d",opfnum);
                            strcat(MOfname,fext);
                            mcout=openCompact(MOfname,"mwb"); /* open new file */
                          }
                      }
                 
                  }
             
                break;


/*----------> SuperCOmPACT End-of-burst structure <----------*/

              case SEOB_STRUCT:                         /* End-of-burst structure */
                zread=rwSuperEndofBurst(cmpin,&seob);   /* read in end of burst */
                if(!zread)                              /* skip burst if problem reading */
                  {
                    JumpToNextFile=1; 
                    break;
                  }
                if( (RMINflag) && (sbur.nrun < RunMin) ) break;
                if( (RMAXflag) && (sbur.nrun > RunMax) ) break;
                userStatus=USER_SUPEREOB(&sbur,&seob);
                fflush(NULL);
             
                if(userStatus>0)                        /* examine filter return for error */
                  {
                    sprintf(buf,"Error code %d returned from user_seob",userStatus);
                    cmpwarn(buf);
                  }
                if(SOflag)                             /* check for SuperCOmPACT output file */
                  rwSuperEndofBurst(cmpout,&seob);     /* and write out end of burst struct */
                if (SMOflag)                           /* check for SuperMC output file */
                  rwSuperEndofBurst(cmpout,&seob);     /* and write out end of burst struct */

                /* if fixed size output files are required check size now */
                break;
             


/*----------> Unknown structure <----------*/

              default:
                if(readBursts)
                  {
                    JumpToNextFile=1;
                    printf(" 2: cmpin->nextStruct=%2d\n",cmpin->nextStruct);

                    cmpmsg("Unknown structure found in input file - try to jump to next file");
                 
                  }
                else
                  {
                    cmperror("Unknown structure found in input file");
                  }
            } /* End swith Structrure type */
        } /* End loop on structures */
/*+++++++++++++++++++ Leave struct loop ++++++++++++++++++++*/
/* 
   Call EndofBurst routine for files which 
   could not be read to the end.
   23-10-98 IWS (from Marco's suggestion)
*/
      if(JumpToNextFile) 
        {
          userStatus=USER_EOB(&bur,&eob);
          fflush(NULL);
           
          if(userStatus>0)              /* examine filter return for error */
            {
              sprintf(buf,"Error code %d returned from user_eob",userStatus);
              cmpwarn(buf);
            }
        }
      JumpToNextFile=0;
      NotOpen=1;                      /* No file open */
      closeCompact(cmpin);            /* close input file */
      while( NotOpen)
        {
          if(FLflag)                      /* reading from input file list */
            { 
              if(current->next)           /* is there another file in the list */
                {  
                  current=current->next;      /* update current pointer to next file */
                  cmpin=openCompact(current->fname,"rb"); /* open new file */
                  iFile++;
                  if (iFile%statusUpd == 0) 
                    {
                      if((fp=fopen(".compact-status","w")) == 0) 
                        {
                          cmpwarn("Failed to write .compact-status file");
                        } 
                      else 
                        {
                          fprintf(fp,"%d %s %d %d",iFile/statusUpd*10,current->fname,readBursts,readEvents);
                        }
                      fclose(fp);
                    } /* End test on iFile%statusUpd */
                }
              else
                {
                  LastFile=TRUE;
                }
            } /* End test on FLflag */
          else if(SIFflag)              /* single input file mode */
            {
              LastFile=TRUE;
            }
          else 
            {
              if(maxBursts==0 || readBursts<maxBursts)   /* check not over maxburst */
                {
                  cmpin=openNextBurst(rundir,&cBurst,eBurst); /* open next file */
                  if(!cmpin)                  /* check another file was found */
                    LastFile=TRUE;
                }
            } /* End test on type on input file */
          if(cmpin==(CMPstream *)0)
            {
              if(LastFile) break;
              cmpmsg(" Could not open file - try to skip to next one ");
            }
          else
            {
              NotOpen=0;
            }
        } /* End of test on not open */
    }/* End loop on files */


/********************* Leave file loop **********************/
  if((fp=fopen(".compact-status","w")) == 0) {
    cmpwarn("Failed to write .compact-status file");
  } else {
    fprintf(fp,"100");
  }
  fclose(fp);
  userStatus=USER_EXIT;              /* call user termination  */
  fflush(NULL);
   
  if(userStatus>0)                   /* examine filter return for error */
    { 
      sprintf(buf,"Error code %d returned from user_exit",userStatus);
      cmpwarn(buf);
    }
  if(COflag)                         /* check for COmPACT output file */
    closeCompact(cmpout);            /* and then close it */
  if(KEOflag)                        /* check for Ke3 output file */
    closeCompact(ke3out);            /* and then close it */
  if(KMUOflag)                       /* check for Ke3 output file */
    closeCompact(kmu3out);           /* and then close it */
  if(MOflag)                         /* check for Mc output file */
    closeCompact(mcout);             /* and then close it */
  if(SOflag)                         /* check for SuperCOmPACT output file */
    closeCompact(scmpout);           /* and then close it */
  if(SMOflag)                        /* check for SuperMC output file */
    closeCompact(smcout);            /* and then close it */
   
  cmpmsg("************ COmPACT finished ************");
  sprintf(buf,"Total number of  bursts read :%9d",readBursts);
  cmpmsg(buf);
  sprintf(buf,"Total number of  events read :%9d",readEvents);
  cmpmsg(buf);
  sprintf(buf,"  Number of CMP  events read :%9d",ncmp);
  cmpmsg(buf);
  sprintf(buf,"  Number of KE3  events read :%9d",nke3);
  cmpmsg(buf);
  sprintf(buf,"  Number of KMU3 events read :%9d",nkmu3);
  cmpmsg(buf);
  sprintf(buf,"  Number of MC   events read :%9d",nmc);
  cmpmsg(buf);
  sprintf(buf,"  Number of SCMP events read :%9d",nscmp);
  cmpmsg(buf);
  sprintf(buf,"  Number of SMC  events read :%9d",nsmc);
  cmpmsg(buf);
  cmpmsg("******************************************");
#if defined(COMPACT_F77) && defined(ALPHA)
/* This call is required by the F77 runtime library on Alphas */
  for_rtl_finish_();   /* shutdown F77 IO on Alpha */
#endif
}
