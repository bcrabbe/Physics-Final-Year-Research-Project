/**************************************************************************/
/* Compact, OptiMized Program for Accessing Countless Terabytes           */
/*                                                                        */
/* COmPACT v7.1 - Main Program                                            */
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
/*                27-04-2000 IWS use -sl to select type of summary        */
/*                31-05-2000 IWS cheat option for MC                      */
/*                22-11-2000 IWS introduce - combined option              */
/*                01-12-2000 IWS skip MC events (compact & SC) with Skip  */
/*                13-12-2000 IWS skip events from end of the bursts for   */
/*                               1998 data                                */
/*                05-03-2001 IWS add check for SMO, SO with opt. -b       */
/**************************************************************************/

#include <sys/types.h>
#include <libgen.h>

#include <string.h>
#include <errno.h>
#include "cmpio.h"
#include "reader.h"
#include "reader_options.h"
#include "database.h"
#include "constants.h"
#include "CompactDataBase.h"

#include "split_selection.h"
#include "scmp2hcmp.h"

#ifndef SEL3PIN
#ifdef __cplusplus 
extern  "C"
#endif
int hypSelNeutral (superBurst *sbur, superCmpEvent *sevt);
#endif

#ifdef USE_ROOT
#include "TROOT.h"
#endif
#define E_BUR 0x1
#define E_LIST 0x2
#define E_EOB 0x4
#define STRINGMAX 50
/* Define geometry structure + a pointer with a better name */
struct GeomCompact geom_;
struct GeomCompact * Geom = &geom_;
char gString[50];
#ifdef COMPACT_F77
struct {
  char gstring[STRINGMAX];
} gstring_;
#endif
#ifdef SCPROD
#  ifdef __cplusplus
extern "C" {
#  endif

  int charged_filter(superBurst *sbur, superCmpEvent *sevt);
  int filter_1tr_(superBurst *,superBurst *,superCmpEvent *,superCmpEvent *);
  int hypSelCharged(superBurst *sbur, superCmpEvent *sevt);
  int hypSelNeutral(superBurst *sbur, superCmpEvent *sevt);

#  ifdef __cplusplus
}
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void printHelp(void);
#ifdef __cplusplus
}
#endif



cmpReader_options_t cmpReader_options;

void cmpReader_options_clear(cmpReader_options_t *opts){
  /* super compact */
  /*
  opts->want_scmp = 0;
  opts->want_scmp_charged = 0;
  opts->want_scmp_neutral = 0;
  */
  opts->scmp_out_name = NULL;
  opts->scmp_out_name_charged = NULL;
  opts->scmp_out_name_neutral = NULL;
  opts->scmp_out_dir = 0;
  opts->scmp_out = FALSE;

  /* hyper compact */
  /*
  opts->want_hcmp = 0;
  opts->want_hcmp_charged = 0;
  opts->want_hcmp_neutral = 0;
  */
  opts->hcmp_out_name = NULL;
  opts->hcmp_out_name_charged = NULL;
  opts->hcmp_out_name_neutral = NULL;
  opts->hcmp_out_dir = 0;
  opts->hcmp_out = FALSE;

  /* MC */
  opts->mc_out_name = NULL;
  opts->mc_out_dir = 0;
  opts->mc_out = FALSE;

  /* super MC */
  opts->smc_out_name = NULL;
  opts->smc_out_dir = 0;
  opts->smc_out = FALSE;

  opts->filter = 0;
  opts->filter_kl2=TRUE;
  opts->abcog=FALSE;
  opts->calcbb=FALSE;
  opts->rewrite_timestamp=FALSE;

  //Behaviour options
  opts->use_mc_version=0;

#ifdef FIX_LKR_TIME
  opts->fix_lkr_toff = 0; /* no correction by default */
#endif
}

void log_skipped_file(char *file_skipped){
  char *filename_skipped = "compact.skipped";
  FILE *fskip;                                 /* Pointer to the list of skipped files. */
  char buf[256];
  static int first_time = 1;
  cmpmsg(" Could not open file - try to skip to next one ");
  //  cmpwarn("skipping file");
  errno = 0;
  if (first_time){ /* truncate the file */
    fskip = fopen(filename_skipped,"w");
    first_time = 0;
  } else { /* append to it */
    fskip = fopen(filename_skipped,"a");
  }
  if (!fskip) {
    sprintf(buf,"log_skipped_file, opening %s",filename_skipped);
    perror(buf);
  } else {
    fprintf(fskip,"%s\n",file_skipped);
    fclose(fskip);
  }
}

void writeExtendedHeader(CMPstream *cmps){
  if (cmps){
    cmps->exthdr.sql_db_ver = sql_db_version.version;
    cmps->exthdr.sql_db_upd = sql_db_version.update;

    cmps->exthdr.spare[0] = 
      cmps->exthdr.spare[1] = 
      cmps->exthdr.spare[2] = 
      cmps->exthdr.spare[3] = 0;

    xdrwarn(zxdr_FExtHeader(cmps),"writeExtendedHeader: Write of file extended header failed");
  }
}

int main(int argc,char **argv) {
  int i,j;
  int userStatus;                           /* user routine return value */
/* Marco.Clemencic@cern.ch === not used (-s -e -n -d options) ===  */
/*   int maxBursts=0;             /\* burst limit *\/ */
  int readBursts=0;                         /* number of bursts read */
  int readEoB=0;                         /* number of bursts read */
  int readEvents=0;
  int ncmp=0,nke3=0,nkmu3=0,nmc=0,nscmp=0;  /* number of event type read */
  int nsmc=0,nhcmp=0;

/* Marco.Clemencic@cern.ch === not used (-s -e -n -d options) ===  */
/*   BurstID sBurst,eBurst,cBurst;              /\* start, end and current *\/ */

//  char COflag=FALSE;           /* option flags */
  char KEOflag=FALSE,KMUOflag=FALSE;         /* option flags */
  //  char MOflag=FALSE;
//  char SOflag=FALSE;            /* CO:COmpact */
//  char SMOflag=FALSE;
  char ETAOflag=FALSE;        /* KEO:Ke3, KMUO:Kmu3, MO:MC */
//  char HOflag=FALSE;
  /* SO:SuperCOmPACT, SMO:SuperMC ETAO: Eta*/
  char SRWflag=FALSE;                        /* rewrite SuperCOmPACT INT event */
#if BFLAG
  char Bflag=FALSE;                          /* produce fixed size o/p files */
#endif
  char NDBflag=TRUE;                        /* Not to use CDB */
  char CHEATflag=FALSE;                      /* to metamorphose a compact event into SC */
  char EMPTYflag=FALSE;                      /* to empty bursts,... when filtering */
/* Marco Clemencic === Compact 7.x cannot read 1998 and 1999 data! */
/*   char COMBflag=FALSE;                       /\* to combined 1998 and 1999 data *\/ */
  long int Empty;                            /* Code to empty burst */
  char NEVTflag=FALSE;                       /* read NevtToRead */
  int  NevtToRead=0;                         /* Number of events to read in */  
  char RMINflag=FALSE;                       /* Only look at burst with Run>=RunMin */
  char RMAXflag=FALSE;                       /* Only look at burst with Run<=RunMax */
  int  RunMin, RunMax;                     
  long int Bsize;                            /* o/p file size (bytes) */
  int opfnum=1;                              /* number of o/p file */
  char fext[10];                             /* filename extension */

/* Marco.Clemencic@cern.ch === This part is not used since ages!!! (-x option) === */
/*   char Xflag=FALSE;                           /\* exec compact.pl flag *\/ */
/*   char xargstr[256]="runcompact ";            /\* exec string for compact.pl *\/ */

#ifndef SCPROD
  char KEOfname[81];   /* names for I/O files */
  char KMUOfname[81];
  char ETAOfname[81];
#else
  char COfname[200],KEOfname[200],MOfname[200];   /* names for I/O files */
  char KMUOfname[200],SOfname[200],SMOfname[200];
  char SOCfname[200],SONfname[200];
  char ETAOfname[200];
  char HOfname[200];
  char HOCfname[200],HONfname[200];
#endif

/* Marco.Clemencic@cern.ch === not used (-s -e -n -d options) ===  */
/*   char rundir[161];               /\* run directory*\/ */

  //#ifdef SCPROD
  char SIfsbase[256],SIfhbase[256];
  //#endif

  char buf[256];
  char *bufp;

  int LastFile=FALSE;                          /* last input file flag */
  int FLflag=FALSE;                            /* file list flag */
  Flist first;                                 /* first file in list */
  Flist *current=&first;                       /* current file in list */

  CMPstream *cmpin,*cmpout,*ke3out,*kmu3out;   /* COmPACT streams */
  CMPstream *mcout,*smcout;
  CMPstream *scmpout,*scmpoutn,*scmpoutc;
  CMPstream *hcmpout,*hcmpoutn,*hcmpoutc;

  /* Burst structures */
  static Burst      bur;
  static superBurst sbur;
  static hyperBurst hbur;

  /* event structures */
  static cmpEvent cmpevt;
  static ke3Event ke3evt;
  static kmu3Event kmu3evt;
  static mcEvent mcevt;

  /* SuperCompact structures */
  static INTsuperCmpEvent supercmpevti;              /* integer versions of SuperCOmPACT */
  static INTsuperMcEvent supermcevti;                /* event structures - for i/o only  */

  static superCmpEvent supercmpevt;                  /* float versions of SuperCOmPACT */
  static superMcEvent supermcevt;                    /* event structures - for users   */

  /* HyperCompact structure */
  static hyperCmpEvent hevt;
  
  /* End Of Burst structures */
  static EndofBurst      eob;
  static superEndofBurst seob;
  

  int level=EPSSUM3;                               /* SuperCOmPACT level: 1/2/3 */

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
  int ComprFactor=4;                       /* compression factor */

  static int timestampPrev=0;                /* Keep timestamp of previous event */
  int LastEventOut=0;                        /* timestamp of last event written out */

  struct stat stat_info;

  bool_t split_flag;
  bool_t split_flag_kl2;

#ifdef COMPACT_F77
#ifdef ALPHA
  /* This call is required by the F77 runtime library on Alphas */
  for_rtl_init_(&argc,argv);         /* init F77 IO for alpha */
#endif   /* leave Alpha specific F77 */
#endif   /* leave F77 specific code */

  /* initialize gString */
  memset(gString,0,50*sizeof(char));
  
/* Marco.Clemencic@cern.ch === not used (-s -e -n -d options) ===  */
/*  strcpy(rundir,".");  */
  if (getenv("COMPACT_SKIP_KABES_RECONSTRUCTION") != 0) {
    cmpwarn("env variable COMPACT_SKIP_KABES_RECONSTRUCTION is set: disabling Kabes reconstruction!");
    COMPACT_SKIP_KABES_RECONSTRUCTION = TRUE;
  }

  if(argc==1) { // no reason to enter the loop if no option is given
    printHelp();
    exit(0);
  }

  /* clear the cmdline options structure */
  cmpReader_options_clear(&cmpReader_options);
  /* initialize the first pointer of the files list */
  first.fname[0] = 0;
  first.next = NULL;

  for(i=1;i<argc;i++) {                /* loop over command line arguments */
/* Marco.Clemencic@cern.ch === they are not in the help file ===  */
/*     if(strcmp(argv[i],"-s") == 0){      /\* search for starting burst *\/      */
/*       sBurst=parseBurstID(argv[++i]); */
/*       continue; */
/*     } */
/*     if(strcmp(argv[i],"-e") == 0){      /\* search for end burst *\/ */
/*       eBurst=parseBurstID(argv[++i]); */
/*       continue; */
/*     } */
/*     if(strcmp(argv[i],"-n") == 0)      /\* search for max bursts *\/ */
/*       {   */
/*         maxBursts=atoi(argv[++i]); */
/*         eBurst.nrun=999999;            /\* make sure eBurst doesn't limit *\/ */
/*         eBurst.nbur=999999;            /\* the size of the run *\/ */
/*         continue; */
/*       } */
/*     if(strcmp(argv[i],"-d") == 0){      /\* search for directory name *\/ */
/*       strcpy(rundir,argv[++i]); */
/*       continue; */
/*     } */
    if(strcmp(argv[i],"-i") == 0){      /* search for input filename */
      tmp = &first; /* init tmp */
      if (first.fname[0] != 0){ /* we have files in the list */
        /* find the last entry */
        while (tmp->next != NULL) tmp=tmp->next;
        /* now tmp is the pointer to the last entry in the list */
        tmp->next = (Flist *)malloc(sizeof(Flist)); /* allocate next flist struct */
        if(!tmp->next) cmperror("Could not allocate memory for input file name!!!");
        tmp = tmp->next;
      }
      /* now tmp is the pointer to an empty entry in the list */
      strcpy(tmp->fname,argv[++i]);
      tmp->next = NULL;
      continue;
    }
    if(strcmp(argv[i],"-l") == 0)      /* search for file name list */
      {
        FLflag=readFlist(argv[++i],&first); /* read file list */
        /* count files */
        while (tmp != NULL) {
          numFiles++;
          tmp=tmp->next;
        }
        statusUpd=numFiles/10+1;
        continue;
      }

/* Marco.Clemencic@cern.ch === This part is not used since ages!!! (-x options) === */
/*     if(strcmp(argv[i],"-x") == 0)      /\* execute compact.pl from compact *\/ */
/*       { */
/*         strcat(xargstr,argv[++i]); */
/*         Xflag=TRUE; */
/*         continue; */
/*       } */
     
    if(strcmp(argv[i],"-ko") == 0)     /* search for output Ke3 filename */
      { 
        strcpy(KEOfname,argv[++i]);
        KEOflag=TRUE;                  /* set COmPACT output file flag */
        continue;
      }
    if(strcmp(argv[i],"-kmo") == 0)    /* search for output Ke3 filename */
      { 
        strcpy(KMUOfname,argv[++i]);
        KMUOflag=TRUE;                 /* set COmPACT output file flag */
        continue;
      }
    /************ OUTPUT OPTIONS **************/
    if(strcmp(argv[i],"-co") == 0){     /* search for output SuperCOmPACT fname */
      cmpReader_options.cmp_out_name = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-so") == 0){     /* search for output SuperCOmPACT fname */
      cmpReader_options.scmp_out_name = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-soc") == 0){     /* search for output SuperCOmPACT fname */
      cmpReader_options.scmp_out_name_charged = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-son") == 0){     /* search for output SuperCOmPACT fname */
      cmpReader_options.scmp_out_name_neutral = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-ho") == 0){     /* search for output hyperCOmPACT fname (charged)*/
      cmpReader_options.hcmp_out_name = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-hoc") == 0){     /* search for output hyperCOmPACT fname (charged)*/
      cmpReader_options.hcmp_out_name_charged = argv[++i];
      //      cmpReader_options.want_hcmp_charged=TRUE;
      continue;
    }
    if(strcmp(argv[i],"-hon") == 0){     /* search for output hyperCOmPACT fname (neutral)*/
      cmpReader_options.hcmp_out_name_neutral = argv[++i];
      //      cmpReader_options.want_hcmp_neutral=TRUE;
      continue;
    }
    if(strcmp(argv[i],"-ho-dir") == 0){     /* output file names for hyperCOmPACT are directories */
      cmpReader_options.hcmp_out_dir=TRUE;
      continue;
    }

    /* Filtering options */
    if(strcmp(argv[i],"-scprod-cuts") ==0){     /* using SCprod default cuts for sel3p[cn] */
      cmpReader_options.filter = 1;
      continue;
    }
    if(strcmp(argv[i],"-scprod-mode") ==0){     /* using SCprod default cuts for sel3p[cn] & split filter for scmp */
      cmpReader_options.filter = 2;
      continue;
    }
    if(strcmp(argv[i],"-nokl2")==0) /*Not to write kl2 events in the neutral split*/
      {
	cmpReader_options.filter_kl2 = FALSE;
	continue;
      }
    if(strcmp(argv[i],"-abcog")==0) { /* Produce kaon mass parameters and COGs */
      cmpReader_options.abcog=TRUE;
      continue;
    }
    if(strcmp(argv[i],"-calcbb")==0) { /* Perform data-quality checks */
      cmpReader_options.calcbb=TRUE;
      continue;
    }
    if(strcmp(argv[i],"-rwts")==0) { /* Rewrite the variable time-stamp of the previous event*/
      cmpReader_options.rewrite_timestamp=TRUE;
      continue;
    }
    if(strcmp(argv[i],"-mo") == 0){    /* search for output MC file name */
      cmpReader_options.mc_out_name = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-smo") == 0){    /* search for output SuperMC file name */
      cmpReader_options.smc_out_name = argv[++i];
      continue;
    }
    if(strcmp(argv[i],"-srw") ==0)     /* rewrite SuperCOmPACT variables */
      {
        SRWflag=TRUE;
        continue;
      }
    if(strcmp(argv[i],"-db") ==0)     /* CDB wanted */
      {
        NDBflag=FALSE;
        continue;
      }
    if(strcmp(argv[i],"-ndb") ==0)     /* CDB not wanted */
      {
        NDBflag=TRUE;
        continue;
      }
    if(strcmp(argv[i],"-cheat") ==0)   /* To metamorphose a compact evt in SC */
      {
        CHEATflag=TRUE;
        continue;
      }
    /********* BEHAVIOUR OPTIONS **********/
    if(strcmp(argv[i],"-skip-kab-rec") ==0)   /* Flag to skip kabes reconstruction in cmp */
      {
        COMPACT_SKIP_KABES_RECONSTRUCTION = TRUE;
        continue;
      }
#ifdef FIX_LKR_TIME
    if(strcmp(argv[i],"-fix-lkr-toff") ==0)   /* Flag to fix LKR cls time offset */
      {
	if ((i+1)<argc) {
	  cmpReader_options.fix_lkr_toff = atoi(argv[++i]);
	  switch ( cmpReader_options.fix_lkr_toff ){
	  case 1:
	    cmpmsg("Fix LKR Time Offset (run 12163 version)");
	    break;
	  case 2:
	    cmpmsg("Fix LKR Time Offset (right run version)");
	    break;
	  case 3:
	    cmpmsg("Fix LKR Time Offset (missing column version)");
	    break;
	  default:
	    cmperror("Option '-fix-lkr-toff' expects a parameter which must be 1, 2 or 3. Default is to not correct.");
	  }
	} else {
	  cmperror("You didn't give the parameter to the option '-fix-lkr-toff', use -h");
	}
        continue;
      }
#endif
    if(strcmp(argv[i],"-empty") ==0)   /* To empty bursts,... when filtering (compact) */
      {
        EMPTYflag=TRUE;
        Empty=atoi(argv[++i]);
        if(Empty==0) Empty |= (E_BUR|E_LIST|E_EOB);
        sprintf(buf," =======> OK! Code to save disk space is %8x",Empty);
        cmpmsg(buf);
        continue;
      }
/* Marco Clemencic === Compact 7.x cannot read 1998 and 1999 data! */
/*     if(strcmp(argv[i],"-combined") ==0)   /\* To combine 1998 and 1999 data *\/ */
/*       { */
/*         COMBflag=TRUE; */
/*         sprintf(buf," =======> You are analysing 1998 and 1999 data together"); */
/*         cmpmsg(buf); */
/*         continue; */
/*       } */

    if(strcmp(argv[i],"-mcver") ==0)   /* Number of events to read in */
      {
	cmpReader_options.use_mc_version=atoi(argv[++i]);
	sprintf(buf,"MC analysis level set to %d ",cmpReader_options.use_mc_version);
	cmpmsg(buf);
	continue;
      }

    if(strcmp(argv[i],"-nevt") ==0)   /* Number of events to read in */
      {
        NEVTflag=TRUE;
        NevtToRead=atoi(argv[++i]);
        sprintf(buf," =======> OK! compact will stop after %5d events",NevtToRead);
        cmpmsg(buf);
        continue;
      }
    if(strcmp(argv[i],"-rmin") ==0)   /* RunMin to look at */
      {
        RMINflag=TRUE;
        RunMin=atoi(argv[++i]);
        sprintf(buf," =======> OK! compact will look at Run>=%5d",RunMin);
        cmpmsg(buf);
        continue;
      }
    if(strcmp(argv[i],"-rmax") ==0)   /* RunMin to look at */
      {
        RMAXflag=TRUE;
        RunMax=atoi(argv[++i]);
        sprintf(buf," =======> OK! compact will look at Run<=%5d",RunMax);
        cmpmsg(buf);
        continue;
      }
    if(strcmp(argv[i],"-string")==0){
      strncpy(gString,argv[++i],49*sizeof(char));
#ifdef COMPACT_F77
      memcpy(gstring_.gstring,gString,50*sizeof(char));
#endif		
      continue;
    }
#ifdef BFLAG
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
        continue;
      }
#endif
    if(strcmp(argv[i],"-h")==0)
      {
        printHelp();
        exit(0);
      }
    /* If we arrived here, the option is not know */
    sprintf(buf," Option '%s' unknown!!! Try with -h",argv[i]);
    cmperror(buf);
  }

  /****************End of the loop on command line options****************/

  /* print the list fo files for debugging */
  /*
  if (first.fname[0] == 0) cmperror("no list of file given");
  for(tmp = &first; tmp != NULL; tmp = tmp->next){
    printf("%s\n",tmp->fname);
  }
  */

  /* consistency check for -rmin and -rmax options */
  if((RMINflag) && (RMAXflag) && (RunMax<RunMin)) 
    {
      sprintf(buf," ======> Incompatible limits for run numbers - stop");
      cmperror(buf);
    }
   
/* Marco.Clemencic@cern.ch === This part is not used since ages!!! (-x option) === */
/*   /\* consistency check for -x option *\/ */
/*   if (Xflag) { */
/*     if (FLflag||SIFflag) cmperror("Please don't use -l or -i with -x option"); */
/*     if(system(xargstr)) cmperror("Failed to exec runcompact"); */
/*     FLflag=readFlist("compact.list",&first); /\* read file list *\/ */
/*     /\* count files *\/ */
/*     while (tmp != NULL) { */
/*       numFiles++; */
/*       tmp=tmp->next; */
/*     } */
/*     statusUpd=numFiles/10+1; */
/*   } */

#ifdef BFLAG
  /* consistency check for -b option */
  if (Bflag)
    {
      if (!(KEOflag||KMUOflag
            ||cmpReader_options.cmp_out
            ||cmpReader_options.scmp_out
            ||cmpReader_options.hcmp_out
	    ||cmpReader_options.mc_out
	    ||cmpReader_options.smc_out
            ))
        {
          cmpmsg("Use -co/-ko/-kmo/-mo/-so/-smo with -b option");
          /* printHelp(); */
          exit(0);
        }

    }
#endif

  /* consistency check for -cheat */
  if (CHEATflag && (KEOflag||KMUOflag
                    ||cmpReader_options.cmp_out
                    ||cmpReader_options.scmp_out
                    ||cmpReader_options.hcmp_out
		    ||cmpReader_options.mc_out
		    ||cmpReader_options.smc_out)) {
    cmpmsg("-cheat is incompatible with any output flags");
    /* printHelp(); */
    exit(0);    
  }

  /* Database initialisation routine */
  if(!NDBflag) getRecordCmpInit();

#ifdef USE_ROOT
  cmpmsg("****************************************");
  cmpmsg("*                                      *");
  cmpmsg("* Setting up compact for use with ROOT *");
  cmpmsg("*                                      *");
  cmpmsg("****************************************");
  TROOT compact("compact","Compact Root Environment COmPACT++");
#endif
   
  NotOpen=1;                         /* set flag before opening file */
  cmpin=openCompact(current->fname,"rb"); /* open current file */
  iFile++;   
   
  if(cmpin==(CMPstream *)0) {
    log_skipped_file(current->fname);
    JumpToNextFile = TRUE;
  }

  if (cmpin) {
    printf("File written with SQL db version %d (upd. %d)\n",cmpin->exthdr.sql_db_ver,cmpin->exthdr.sql_db_upd);
  }
  
  sqlite_init(); /* open the database before the output files in order to set properly the version of the db in the header */

  if(KEOflag)                             /* writing an output file */
    {
#ifdef BFLAG
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(KEOfname,fext);
        }
#endif
      ke3out=openCompact(KEOfname,"kwb"); /* open COmPACT Ke3 file */
      if(ke3out == (CMPstream *)0)
        cmperror("compact: Failed to open Ke3 COmPACT output stream");
      writeExtendedHeader(ke3out);
    }
   
  if(KMUOflag)                            /* writing an output file */
    {
#ifdef BFLAG
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(KMUOfname,fext);
        }
#endif
      kmu3out=openCompact(KMUOfname,"lwb");/* open COmPACT Kmu3 file */
      if(kmu3out == (CMPstream *)0)
        cmperror("compact: Failed to open Kmu3 COmPACT output stream");
      writeExtendedHeader(kmu3out);
    }
   
#ifdef BFLAG
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
#endif
  
#ifdef BFLAG
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
#endif

  /***** COmPACT output files *****/
  cmpout = NULL;
  if (cmpReader_options.cmp_out_name) {
    cmpout=openCompact(cmpReader_options.cmp_out_name,"wb"); /* open COmPACT file (user filter) */
    if(cmpout == NULL)
      cmperror("compact: Failed to open Compact output stream");
    writeExtendedHeader(cmpout);
    cmpReader_options.cmp_out = TRUE;
  }
  /***** superCOmPACT output files *****/
  scmpout = scmpoutc = scmpoutn = NULL;
  if (cmpReader_options.scmp_out_name) {
    scmpout=openCompact(cmpReader_options.scmp_out_name,"swb"); /* open superCOmPACT file (user filter) */
    if(scmpout == NULL)
      cmperror("compact: Failed to open SuperCompact output stream");
    writeExtendedHeader(scmpout);
    cmpReader_options.scmp_out = TRUE;
  }
  if (cmpReader_options.scmp_out_name_charged) {
    scmpoutc=openCompact(cmpReader_options.scmp_out_name_charged,"swb"); /* open superCOmPACT file (charged filter) */
    if(scmpoutc == NULL)
      cmperror("compact: Failed to open SuperCompact output stream");
    writeExtendedHeader(scmpoutc);
    cmpReader_options.scmp_out = TRUE;
  }
  if (cmpReader_options.scmp_out_name_neutral) {
    scmpoutn=openCompact(cmpReader_options.scmp_out_name_neutral,"swb"); /* open superCOmPACT file (neutral filter) */
    if(scmpoutn == NULL)
      cmperror("compact: Failed to open SuperCompact output stream");
    writeExtendedHeader(scmpoutn);
    cmpReader_options.scmp_out = TRUE;
  }
  /***** COmPACT MC output files *****/
  mcout = NULL;
  if (cmpReader_options.mc_out_name) {
    mcout=openCompact(cmpReader_options.mc_out_name,"mwb"); /* open COmPACT MC file (conversion) */
    if(mcout == NULL)
      cmperror("compact: Failed to open Compact MC output stream");
    writeExtendedHeader(mcout);
    cmpReader_options.mc_out = TRUE;
  }
  /***** superCOmPACT MC output files *****/
  smcout = NULL;
  if (cmpReader_options.smc_out_name) {
    smcout=openCompact(cmpReader_options.smc_out_name,"nwb"); /* open superCOmPACT MC file (conversion) */
    if(smcout == NULL)
      cmperror("compact: Failed to open SuperCompact MC output stream");
    writeExtendedHeader(smcout);
    cmpReader_options.smc_out = TRUE;
  }
  /***** hyperCOmPACT output files *****/
  hcmpout = hcmpoutc = hcmpoutn = NULL;
  if (cmpReader_options.hcmp_out_name) {
    hcmpout=openCompact(cmpReader_options.hcmp_out_name,"hwb"); /* open superCOmPACT file (user filter) */
    if(hcmpout == NULL)
      cmperror("compact: Failed to open SuperCompact output stream");
    writeExtendedHeader(hcmpout);
    cmpReader_options.hcmp_out = TRUE;
  }
  if (cmpReader_options.hcmp_out_name_charged) {
    hcmpoutc=openCompact(cmpReader_options.hcmp_out_name_charged,"hwb"); /* open superCOmPACT file (charged filter) */
    if(hcmpoutc == NULL)
      cmperror("compact: Failed to open SuperCompact output stream");
    writeExtendedHeader(hcmpoutc);
    cmpReader_options.hcmp_out = TRUE;
  }
  if (cmpReader_options.hcmp_out_name_neutral) {
    hcmpoutn=openCompact(cmpReader_options.hcmp_out_name_neutral,"hwb"); /* open superCOmPACT file (neutral filter) */
    if(hcmpoutn == NULL)
      cmperror("compact: Failed to open SuperCompact output stream");
    writeExtendedHeader(hcmpoutn);
    cmpReader_options.hcmp_out = TRUE;
  }

#if 0
  if(SOflag)                             /* writing an output file */
    {
      if (Bflag) 
        {
          sprintf(fext,".%d",opfnum);
          strcat(SOfname,fext);
        }
#ifndef SCPROD
      scmpout=openCompact(SOfname,"swb"); /* open SuperCOmPACT file */
      if(scmpout == (CMPstream *)0)
        cmperror("compact: Failed to open SuperCompact output stream");
#else
      // check the presence of output directories
      sprintf(SOCfname,"%s/charged",SOfname);
      if (rfio_stat(SOCfname,&stat_info)) {
        bufp = basename(argv[0]);
        sprintf(buf,"%s: stat(\"%s\")",bufp,SOCfname);
        rfio_perror(buf);
        free(bufp);
        exit(1);
      } else {
        if (!S_ISDIR(stat_info.st_mode)){
          sprintf(buf,"I expected a directory called %s",SOCfname);
          cmperror(buf);
        }
      }
      // check the presence of output directories
      sprintf(SONfname,"%s/neutral",SOfname);
      if (rfio_stat(SONfname,&stat_info)) {
        bufp = basename(argv[0]);
        sprintf(buf,"%s: stat(\"%s\")",bufp,SONfname);
        rfio_perror(buf);
        free(bufp);
        exit(1);
      } else {
        if (!S_ISDIR(stat_info.st_mode)){
          sprintf(buf,"I expected a directory called %s",SONfname);
          cmperror(buf);
        }
      }
      /* build the outputname from the SOfname varibale */
      sprintf(SOCfname,"%s/charged/%s",SOfname,SIfsbase);
      sprintf(SONfname,"%s/neutral/%s",SOfname,SIfsbase);
      scmpoutc=openCompact(SOCfname,"swb"); /* open SuperCOmPACT file */
      if(scmpoutc == (CMPstream *)0){
        cmperror("compact: Failed to open charged filter output stream");
      }
      scmpoutn=openCompact(SONfname,"swb"); /* open SuperCOmPACT file */
      if(scmpoutn == (CMPstream *)0){
        cmperror("compact: Failed to open neutral filter output stream");
      }
#endif
    }
#endif

#ifdef BFLAG
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
#endif
 
#if 0
  if(cmpReader_options.want_hcmp_charged ||
     cmpReader_options.want_hcmp_neutral) {                             /* writing an output file */
    /*
   if (Bflag) 
      {
        sprintf(fext,".%d",opfnum);
        strcat(HOfname,fext);
      }
    */
  }
#endif
  
  cmpAnaInit();                          /* initialise analysis routines */
  superCmpAnaInit();                     /* initialise analysis routines */

  if ((cmpReader_options.filter >= 1) && (!cmpReader_options.abcog)) {
    cmpmsg("Initilizing ana3pic_cut and ana3pin_cut for SCProd filtering.");
    init_sel3pic_cuts_scprod();            /* re-initialize cuts for hypercmp production */
    init_sel3pin_cuts_scprod();            /* re-initialize cuts for hypercmp production */
  }
  if(cmpReader_options.abcog) {
    cmpmsg("Initilizing Alpha, Beta and COG calculation routine\n");
    cmpwarn("ABCOG cuts override SCPROD cuts for ana3pic !!! \n");
    abc_INIT;
  }

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
      while( (!JumpToNextFile) && cmpin && (cmpin->nextStruct != FTLR_STRUCT) ) 
        {
          if((NEVTflag) && (readEvents>=NevtToRead)) 
            {
              LastFile=TRUE;
              break; /* Break Loop */
            }
          switch(cmpin->nextStruct) 
            {
         
	
              /*----------> Burst structure <----------*/

            case BURST_STRUCT:            /* burst structure */
              zread=rwBurst(cmpin,&bur);  /* read in burst header */
              if(!zread)                  /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
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
              if (cmpReader_options.scmp_out
		  || cmpReader_options.smc_out
		  || CHEATflag) cmp2scmpInitBurst(&bur);
	      if (cmpReader_options.calcbb) { /*Bad bursts calculations necessary - call ana routines*/
		MUV_BURST(&bur); /* MUV bad burst */
	      }
              userStatus=USER_BURST(&bur); /* call user burst routine */
              cmpAnaPostBurst(&bur);       /* Compact analysis for burst after user routine */
              if(userStatus>0) 
                {
                  sprintf(buf,"Error code %d returned from user_burst",userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);
	      /*Let's have some information if the event is first in a burst*/
	      /* Here, to be compatible also with -cheat option */
	      if(cmpReader_options.rewrite_timestamp) 
		{
		  timestampPrev=-1;
		}
              
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
		  if(cmpReader_options.abcog) { /*need alpha, beta and COG constants*/
		    abc_SUPERBURST(&sbur);      
		  }
                  userStatus=USER_SUPERBURST(&sbur);  /* call user burst routine */
/* Marco Clemencic === Compact 7.x cannot read 1998 and 1999 data! */
/*                   if(COMBflag) sbur.dataset=1;       /\* Save info that data are combined *\/ */
                  superCmpAnaPostBurst(&sbur);
                  if(userStatus>0) 
                    {
                      sprintf(buf,"Error code %d returned from user_superBurst",userStatus);
                      cmpwarn(buf);
                    }
#ifdef HCHEAT
                  scmp2hcmpBurst(&sbur,&hbur);  /* copy superBurst variables */
                  userStatus=USER_HYPERBURST(&hbur);  /* call user burst routine */
                  if(userStatus>0) 
                    {
                      sprintf(buf,"Error code %d returned from user_hyperBurst",userStatus);
                      cmpwarn(buf);
                    }
#endif
                  break; /* CHEAT prevents output */
                } /* End test on CHEAT flag */

              if((EMPTYflag) && (Empty&E_BUR))
                {
                  EmptyBurst(&bur,100);
                }
              if (cmpReader_options.cmp_out) {
                if (cmpout) rwBurst(cmpout,&bur);      /* write out Burst */
              }
	      if (cmpReader_options.mc_out) {
                if (mcout) rwBurst(mcout,&bur);      /* write out Burst */
              }
	      if (cmpReader_options.scmp_out
		  || cmpReader_options.smc_out ) {
                cmp2scmpBurst(&bur,&sbur);  /* copy superBurst variables */
                if (scmpout) rwSuperBurst(scmpout,&sbur);      /* write out superBurst */
		if (smcout)   rwSuperBurst(smcout,&sbur);
		if(cmpReader_options.abcog) { /*need alpha, beta and COG constants*/
		  superCmpAnaPreBurst(&sbur);
		  abc_SUPERBURST(&sbur);      
		}

              }

              if(KEOflag)
                rwBurst(ke3out,&bur);
             
              if(KMUOflag)
                rwBurst(kmu3out,&bur);
             
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
              /* 05-02-99 IWS dont put in this protection yet */
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
	      if(cmpReader_options.abcog) { /*need alpha, beta and COG constants*/
		abc_SUPERBURST(&sbur);      
	      }
              readBursts++;                       /* increase number of bursts read */
              userStatus=USER_SUPERBURST(&sbur);  /* call user burst routine */
              superCmpAnaPostBurst(&sbur);        /* Compact analysis for burst before user routine */
              if(userStatus>0)
                {
                  sprintf(buf,"Error code %d returned from user_superBurst",userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);
#ifdef HCHEAT
              if(CHEATflag) {             
                scmp2hcmpBurst(&sbur,&hbur);  /* copy superBurst variables */
                userStatus=USER_HYPERBURST(&hbur);  /* call user burst routine */
                if(userStatus>0) 
                  {
                    sprintf(buf,"Error code %d returned from user_hyperBurst",userStatus);
                    cmpwarn(buf);
                  }
                break;
              }
#endif
	      /*Let's have some information if the event is first in a burst*/
	      if(cmpReader_options.rewrite_timestamp) 
		{
		  timestampPrev=-1;
		}
	                    
              if (cmpReader_options.scmp_out) {
                if (scmpout) rwSuperBurst(scmpout,&sbur);      /* write out superBurst */
                if (scmpoutc) rwSuperBurst(scmpoutc,&sbur);    /* write out superBurst (charged stream) */
                if (scmpoutn) rwSuperBurst(scmpoutn,&sbur);    /* write out superBurst (neutral stream) */
              }
	      if (cmpReader_options.smc_out) {
		if (smcout)   rwSuperBurst(smcout,&sbur);      /* write out superBurst (MC) */
              }
              if (cmpReader_options.hcmp_out) {
                scmp2hcmpBurst(&sbur,&hbur);
                if (hcmpout) rwHyperBurst(hcmpout,&hbur);      /* write out hyperBurst */
                if (hcmpoutc) rwHyperBurst(hcmpoutc,&hbur);    /* write out hyperBurst (charged stream) */
                if (hcmpoutn) rwHyperBurst(hcmpoutn,&hbur);    /* write out hyperBurst (neutral stream) */
              }

#ifdef BFLAG
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
#endif
              break;
           
            case HBURST_STRUCT:
              zread=rwHyperBurst(cmpin,&hbur);
              if(!zread){
                JumpToNextFile=1;
                break;
              }
              if( (RMINflag) && (hbur.nrun < RunMin) ) break;
              if( (RMAXflag) && (hbur.nrun > RunMax) ) break;
              
              hyperCmpAnaPreBurst(&hbur); /* Compact analysis for burst before user routine */
              readBursts++;                       /* increase number of bursts read */
              userStatus=USER_HYPERBURST(&hbur);  /* call user burst routine */
              //hyperCmpAnaPostBurst(&hbur);        /* Compact analysis for burst before user routine */
              if(userStatus>0) {
                sprintf(buf,"Error code %d returned from user_hyperBurst",userStatus);
                cmpwarn(buf);
              }
              fflush(NULL);

              if (cmpReader_options.hcmp_out) {
                if (hcmpout) rwHyperBurst(hcmpout,&hbur);      /* write out hyperBurst */
                /* no std filtering allowed */
              }
              break;
           
              /*----------> Compact event structure <----------*/
           
            case CMPEVENT_STRUCT:                   /* normal COmPACT event */
              zread=rwCmpEvent(cmpin,&cmpevt,&bur); /* read event */
              if(!zread)                            /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
                }
              if(bur.BadB.Skip) break; /* Dont analyse events from bad bursts */
              /* 13-12-2000 Remove events at end of bursts for 1998 data - 3pi0 */ 
	      /* 23/7/2004 Marco Clemencic: compact 7.x is not supposed to read them!! */
              /* if( (bur.nrun>RUN_MIN_1998) && (bur.nrun<RUN_MAX_1998) && (!NDBflag) && (cmpevt.n > bdb->nEvt)) break; */
              if( (RMINflag) && (bur.nrun < RunMin) ) break;
              if( (RMAXflag) && (bur.nrun > RunMax) ) break;
              if( (timestampPrev == cmpevt.timestamp) && (cmpevt.n != 0))
                {
                  zread=rwCmpEvent(cmpin,&cmpevt,&bur); /* read event */
                  if(!zread)                            /* skip burst if problem reading */
                    {
                      /* 03-12-98 May be EoB reached */
                      zread=rwEndofBurst(cmpin,&eob);     /* read in end of burst */
                      /* 15-02-99 to finish with this burst */
                      JumpToNextFile=1;
                      break;
                    }
                }
	      /*Hack to write the time-stamp of the previous event */
	      if(cmpReader_options.rewrite_timestamp) 
		{
		  /*Should be used only for main production of SC (scprod-mode)*/
		  cmpevt.spareInt[1]=timestampPrev;
		}
              timestampPrev = cmpevt.timestamp;
                 
              cmpAnaEvent(&bur,&cmpevt);            /* Compact analysis / event */
#ifdef COMPACT_F77
              if(ncmp==0)                           /* run check on first event read */
                {         
                  cmpmsg("main: First cmpEvent found - running C/F77 checks");
                  chkF77CmpEvent(&cmpevt);
                }
#endif
              ncmp++;                     /* increment number of COmPACT events */
              readEvents++;
	      if (cmpReader_options.calcbb) { /*Bad bursts calculations necessary - call ana routines*/
		MUV_EVENT(&bur,&cmpevt); /* MUV bad burst */
	      }
              userStatus=USER_CMPEVENT(&bur,&cmpevt);
              if(userStatus>0)              /* check user return for error */
                {      
                  sprintf(buf,"Error code %d returned from user_cmpEvent",
                          userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);
           
              if(CHEATflag)
                {
#ifdef COMPACT7
                  cmp2scmp(&bur,&cmpevt,&supercmpevti); /* convert event */
#else
                  cmp2scmp(&bur,&cmpevt,&supercmpevti,level); /* convert event */
#endif
                  scmp2float(&supercmpevti,&supercmpevt,&sbur);      /* convert user variables to floats */
                  calcScmpEvent(&sbur,&supercmpevt);
                  superCmpAnaEvent(&sbur,&supercmpevt); /* Analysis routines for super-compact */
		  if(cmpReader_options.abcog) { /*Need alpha, beta and COG constants*/
		    /*Before User Routine*/
		    abc_SUPERCMPEVENT(&sbur,&supercmpevt);
		  }		    
                  userStatus=USER_SUPERCMPEVENT(&sbur,&supercmpevt); /* <-- float version */

#ifdef HCHEAT
                  
#endif
                  break;
                } /* End test on CHEAT flag */
              
              
              /*************** user filtering required */
              if (cmpReader_options.cmp_out_name || cmpReader_options.scmp_out_name) {
                userStatus=USER_CMPFILTER(&bur,&cmpevt);
                fflush(NULL);
                if(userStatus>0){         /* examine filter return for error */
                  sprintf(buf,"Error code %d returned from user_cmpFilter",userStatus);
                  cmpwarn(buf);
                } else if(userStatus<0) {         /* check if event has passed filter */
                  if (cmpout){
                    if((EMPTYflag)&&(Empty&E_LIST)) {       /* Dont write out event lists */
                      cmpevt.NTrigWordL       = 0;
                      cmpevt.NTimeStampL      = 0;
                      cmpevt.NDCHDecErrorL    = 0;
                      cmpevt.NLKRHACDecErrorL = 0;
                      cmpevt.NChargedInfoL    = 0;
                      cmpevt.NNeutralInfoL    = 0;
                    } else {
                      BuildEvtList(&bur, &cmpevt, LastEventOut);
                      LastEventOut=cmpevt.timestamp;
                    }
                    rwCmpEvent(cmpout,&cmpevt,&bur); /* write out good event */
                  }
                  if (scmpout){
                    cmp2scmp(&bur,&cmpevt,&supercmpevti); /* convert event */
                    rwSuperCmpEvent(scmpout,&supercmpevti,&sbur);     /* write out event */
		    if(cmpReader_options.abcog) { /*Need alpha, beta and COG constants*/
		      scmp2float(&supercmpevti,&supercmpevt,&sbur);      /* convert user variables to floats */
		      calcScmpEvent(&sbur,&supercmpevt);
		      superCmpAnaEvent(&sbur,&supercmpevt); /* Analysis routines for super-compact */
		      abc_SUPERCMPEVENT(&sbur,&supercmpevt);
		    }		    
                  }
                }
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
              if(bur.BadB.Skip) break; /* Dont analyse events from bad bursts */
              if( (RMINflag) && (bur.nrun < RunMin) ) break;
              if( (RMAXflag) && (bur.nrun > RunMax) ) break;
              cmpAnaEvent(&bur,&mcevt.cmpevt);
#ifdef COMPACT_F77
              if(nmc==0) {                           /* run check on first event read */
                cmpmsg("main: First mcEvent found - running C/F77 checks");
                chkF77McEvent(&mcevt);
              }
#endif
              nmc++;                             /* increment num of MC events read */
              readEvents++;
              userStatus=USER_MCEVENT(&bur,&mcevt);
              if(userStatus>0)                      /* check user return for error */
                {
                  sprintf(buf,"Error code %d returned from user_McEvent",
                          userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);

              /* In case of CHEAT */
              if(CHEATflag)
                {
                  mc2smc(&bur,&mcevt,&supermcevti); /* convert event */
                  smc2float(&supermcevti,&supermcevt,&sbur);      /* convert user variables to floats */
                  calcScmpEvent(&sbur,&supermcevt.scmpevt);
                  superCmpAnaEvent(&sbur,&(supermcevt.scmpevt)); /* Analysis routines for super-compact */
		  if(cmpReader_options.abcog) { /*Need alpha, beta and COG constants*/
		    /*Before User Routine*/
		    abc_SUPERCMPEVENT(&sbur,&supermcevt.scmpevt);
		  }		    
                  userStatus=USER_SUPERMCEVENT(&sbur,&supermcevt); /* <-- float version */
                  break;
                } 
              /* End test on CHEAT flag */
             
              /*************** user filtering required */
              if (cmpReader_options.mc_out_name || cmpReader_options.smc_out_name) {
		userStatus=USER_MCFILTER(&bur,&mcevt);
		fflush(NULL);                 
		if(userStatus>0){                  /* examine filter return for error */
		  sprintf(buf,"Error code %d returned from user_mcFilter",
			  userStatus);
		  cmpwarn(buf);
		} else if(userStatus<0) {              /* check if event has passed filter */
		  if (mcout) rwMcEvent(mcout,&mcevt,&bur);  /* write out good event */
		  if (smcout) {
		    mc2smc(&bur,&mcevt,&supermcevti);       /* convert event */
		    rwSuperMcEvent(smcout,&supermcevti,&sbur);
		    if(cmpReader_options.abcog) { /*Need alpha, beta and COG constants*/
		      /*Before User Routine*/
		      smc2float(&supermcevti,&supermcevt,&sbur);      /* convert user variables to floats */
		      calcScmpEvent(&sbur,&supermcevt.scmpevt);
		      superCmpAnaEvent(&sbur,&(supermcevt.scmpevt)); /* Analysis routines for super-compact */
		      abc_SUPERCMPEVENT(&sbur,&supermcevt.scmpevt);
		    }		    
		  }
                }
              }
              break;
             
           
           
              /*----------> Ke3 Compact event structure <----------*/
           
            case KE3EVENT_STRUCT:               /* Ke3 event */
              zread=rwKe3Event(cmpin,&ke3evt);  /* read event */
              if(!zread)                        /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
                }
              if(bur.BadB.Skip) break; /* Dont analyse events from bad bursts */
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
              if(bur.BadB.Skip) break; /* Dont analyse events from bad bursts */
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
              zread=rwSuperCmpEvent(cmpin,&supercmpevti,&sbur);     /* read event */
              if(!zread)                                  /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
                }
              if(sbur.BadB.Skip) break; /* Dont analyse events from bad bursts 16-03-99 */
              if( (RMINflag) && (sbur.nrun < RunMin) ) break;
              if( (RMAXflag) && (sbur.nrun > RunMax) ) break;

#ifdef FIX_LKR_TIME
	      /* MarcoC: This fix can be written better... but I do not have time*/
	      switch ( cmpReader_options.fix_lkr_toff ){ /* FIX LKR time offset */
	      case 1: /* run 12163 */
		{
		  unsigned int i;
		  int corr = (int)((-(15.8359-40.5992) + sbur.tOffst.Lkr)*SCF_LKRTIME);
		  for (i=0; i<supercmpevti.Ncluster; i++){
		    supercmpevti.cluster[i].time += corr;
		  }
		  /*
		  register float t;
		  for (i=0; i<supercmpevti.Ncluster; i++){
		    t = ((float)supercmpevti->cluster[i].time)/SCF_LKRTIME;
		    t = t - (15.8359-40.5992) + sbur.tOffst.Lkr;
		    supercmpevti.cluster[i].time = (int) roundf(t*SCF_LKRTIME);
		  }
		  */
		}
		break;
	      case 2: /* right run */
		{
		  unsigned int i;
		  int corr = (int)((-(sbur.tOffst.Tag-sbur.tOffst.LkrTag) + sbur.tOffst.Lkr)*SCF_LKRTIME);
		  for (i=0; i<supercmpevti.Ncluster; i++){
		    supercmpevti.cluster[i].time += corr;
		  }
		}
		break;
	      case 3: /* missing column */
		{
		  unsigned int i;
		  int corr = (int)((-(sbur.tOffst.Tag-sbur.tOffst.Muv) + sbur.tOffst.Lkr)*SCF_LKRTIME);
		  for (i=0; i<supercmpevti.Ncluster; i++){
		    supercmpevti.cluster[i].time += corr;
		  }
		}
		break;
	      }
#endif
              scmp2float(&supercmpevti,&supercmpevt,&sbur);      /* convert user variables to floats */
              calcScmpEvent(&sbur,&supercmpevt);
#ifdef COMPACT_F77
              if(nscmp==0)                                /* run check on first event read */
                {         
                  cmpmsg("main: First superCmpEvent found - running C/F77 checks");
                  chkF77SuperCmpEvent(&supercmpevt);
                }
#endif
              superCmpAnaEvent(&sbur,&supercmpevt); /* Analysis routines for super-compact */
	      if(cmpReader_options.rewrite_timestamp) 
		{
		  /*Should be used only for main production of SC (scprod-mode)*/
		  supercmpevt.tsPrev=timestampPrev;
		}
              timestampPrev = supercmpevt.timeStamp;

	      if(cmpReader_options.abcog) { /*Need alpha, beta and COG constants*/
		abc_SUPERCMPEVENT(&sbur,&supercmpevt);
	      }
              nscmp++;                                   /* increment number of SuperCOmPACT events */
              readEvents++;
              userStatus=USER_SUPERCMPEVENT(&sbur,&supercmpevt); /* <-- float version */
              if(userStatus>0)                            /* check user return for error */
                {      
                  sprintf(buf,"Error code %d returned from user_superCmpEvent",
                          userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);
             
#ifdef HCHEAT
              if(CHEATflag){
                /* call selection routines for charged and neutral seperately */
                hevt.flag=-1; 
                if(hypSelCharged(&sbur,&supercmpevt)){
                  hevt.flag=0; 
                  scmp2hcmp(&sbur,&supercmpevt,&hevt);
                  userStatus = USER_HYPERCMPEVENT(&hbur,&hevt);
                  if(userStatus>0)                            /* check user return for error */
                    {      
                      sprintf(buf,"Error code %d returned from user_hyperCmpEvent (charged)",
                              userStatus);
                      cmpwarn(buf);
                    }
                  fflush(NULL);
                }

                hevt.flag=-1; 
                if(hypSelNeutral(&sbur,&supercmpevt)){
                  hevt.flag=1; /* make sure we call the neutral filling routine */
                  scmp2hcmp(&sbur,&supercmpevt,&hevt);
                  userStatus = USER_HYPERCMPEVENT(&hbur,&hevt);
                  if(userStatus>0)                            /* check user return for error */
                    {      
                      sprintf(buf,"Error code %d returned from user_hyperCmpEvent (neutral)",
                              userStatus);
                      cmpwarn(buf);
                    }
                  fflush(NULL);
                }
                break;
              }
#endif
              /*************** user filtering required */
              if (cmpReader_options.scmp_out_name || cmpReader_options.hcmp_out_name) {
                memset(&hevt,0,sizeof(hyperCmpEvent)); /* clean it -> the user must fill it */
                userStatus=USER_SUPERCMPFILTER(&sbur,&supercmpevt,&hevt);  /* <-- float version */
                fflush(NULL);
                if(userStatus>0){                       /* examine filter return for error */
                  sprintf(buf,"Error code %d returned from user_superCmpFilter", userStatus);
                  cmpwarn(buf);
                } else if(userStatus<0) { /* the user want this event on disk */
                  if (SRWflag) scmp2scmp(&supercmpevt,&supercmpevti); /* change integer structure */
                  if (scmpout) rwSuperCmpEvent(scmpout,&supercmpevti,&sbur);/* write out good event */
                  if (hcmpout) rwHyperCmpEvent(hcmpout,&hevt,&hbur); /* no conversion needed before writing */
                                                               /*     (the user should fill hevt)     */
                }
              }
              /*************** charged filtering required */
              if (cmpReader_options.scmp_out_name_charged || cmpReader_options.hcmp_out_name_charged) {
		split_flag = charged_split_selection(&sbur,&supercmpevt); /* this event should go to the split (C) */
		if (split_flag) {
		  /* run charged selection */
		  ana3pic.flag = 0;
		  user_superSel3pic(&sbur,&supercmpevt);

		  if ((scmpoutc) /* scmp out requested */
		      && ((cmpReader_options.filter == 2) /* want scmp split */
			  || (ana3pic.flag))              /* or is a good event */
		      )  rwSuperCmpEvent(scmpoutc,&supercmpevti,&sbur); /* write out scmp event */
		  if ((hcmpoutc) /* hcmp requested */
		      && (ana3pic.flag)){ /* the hcmp must be a good one (no split) */
		    FillHyperCharged(&sbur,&supercmpevt,&hevt);
		    rwHyperCmpEvent(hcmpoutc,&hevt,&hbur);
		  }
		}
              }
              /*************** neutral filtering required */
              if (cmpReader_options.scmp_out_name_neutral || cmpReader_options.hcmp_out_name_neutral) {
                /* check if this event should go to the split (N) */
		split_flag = neutral_split_selection(&sbur,&supercmpevt);  /* neutral 3pi */
                if ((!split_flag) && (cmpReader_options.filter == 2) && (cmpReader_options.filter_kl2)) {
                  split_flag_kl2 = ( kl2filter_(&sbur,&sbur,&supercmpevt,&supercmpevt) != 0 );/* filter by Luca Fiorini */
                } else {
                  split_flag_kl2 = 0;
                }

		if (split_flag||split_flag_kl2) {
		  /* run neutral selection */
		  ana3pin.flag = 0;
                  if (split_flag) { /* check neutral filter only for std split (not kl2)*/
#ifdef SEL3PIN
                    user_superSel3pin(&sbur,&supercmpevt);
#else
                    ana3pin.flag = hypSelNeutral(&sbur,&supercmpevt);
#endif
                  }
		  if ((scmpoutn) /* scmp out requested */
		      && ((cmpReader_options.filter == 2) /* want scmp split */
			  || (ana3pin.flag))              /* or is a good event */
		      )  rwSuperCmpEvent(scmpoutn,&supercmpevti,&sbur); /* write out scmp event */
		  if ((hcmpoutn) /* hcmp requested */
		      && (ana3pin.flag)){ /* the hcmp must be a good one (no split) */
		    FillHyperNeutral(&sbur,&supercmpevt,&hevt);
		    rwHyperCmpEvent(hcmpoutn,&hevt,&hbur);
		  }
		}
              }


              if(cmpReader_options.cmp_out)                                 /* is the COmPACT output flag set? */
                {
                  sprintf(buf,"Cannot write COmPACT o/p from SuperCOmPACT i/p file");
                  cmpwarn(buf);
                }
              break;
             
           
           
              /*----------> SuperMC event structure <----------*/

            case SMCEVENT_STRUCT:                        /* SuperMC event */
              zread=rwSuperMcEvent(cmpin,&supermcevti,&sbur);  /* read event */
              if(!zread)                                 /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
                }
              if(sbur.BadB.Skip) break; /* Dont analyse events from bad bursts */
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
              calcScmpEvent(&sbur,&(supermcevt.scmpevt));
              superCmpAnaEvent(&sbur,&(supermcevt.scmpevt)); /* Analysis routines for super-compact */
	      if(cmpReader_options.abcog) { /*Need alpha, beta and COG constants*/
		abc_SUPERCMPEVENT(&sbur,&(supermcevt.scmpevt));
	      }
              userStatus=USER_SUPERMCEVENT(&sbur,&supermcevt); /* <-- float version */
              nsmc++;                                   /* increment number of SuperMC events */
              readEvents++;
              fflush(NULL);
             
              if(userStatus>0)                           /* check user return for error */
                {      
                  sprintf(buf,"Error code %d returned from user_superMcEvent",
                          userStatus);
                  cmpwarn(buf);
                }
	      /*************** user filtering required */
              if(cmpReader_options.smc_out_name) {
		userStatus=USER_SUPERMCFILTER(&sbur,&supermcevt,&hevt);
		fflush(NULL);                 
		if(userStatus>0){                  /* examine filter return for error */
		  sprintf(buf,"Error code %d returned from user_superMcFilter",
			  userStatus);
		  cmpwarn(buf);
		} else if(userStatus<0) {              /* check if event has passed filter */
		  if (smcout) rwSuperMcEvent(smcout,&supermcevti,&sbur);
                }
              }
              /*************** charged filtering required */
              if (cmpReader_options.hcmp_out_name_charged) {
		split_flag = charged_split_selection(&sbur,&(supermcevt.scmpevt)); /* this event should go to the split (C) */
                ana3pic.flag = 0;
		if (split_flag) {
		  /* run charged selection */
		  user_superSel3pic(&sbur,&(supermcevt.scmpevt));
		}
                if (hcmpoutc) { /* hcmp requested */
                  if (ana3pic.flag){ /* the hcmp must be a good one (no split) */
                    FillHyperCharged(&sbur,&(supermcevt.scmpevt),&hevt);
                  } else { /* be sure that the structure is empty */
                    memset(&hevt,0,sizeof(hyperCmpEvent));
                  }
                  /* for MC events we need an output even if it is not passing the filter */
                  FillHyperMCParticles(&sbur,&supermcevt,&hevt);
                  rwHyperCmpEvent(hcmpoutc,&hevt,&hbur);
                }
              }
              /*************** neutral filtering required */
              if (cmpReader_options.hcmp_out_name_neutral) {
		split_flag = neutral_split_selection(&sbur,&(supermcevt.scmpevt)); /* this event should go to the split (N) */
                ana3pin.flag = 0;
		if (split_flag) {
		  /* run neutral selection */
#ifdef SEL3PIN
		  user_superSel3pin(&sbur,&(supermcevt.scmpevt));
#else
		  ana3pin.flag = hypSelNeutral(&sbur,&(supermcevt.scmpevt));
#endif
		}
                if (hcmpoutn) { /* hcmp requested */
                  if (ana3pin.flag){ /* the hcmp must be a good one (no split) */
                    FillHyperNeutral(&sbur,&(supermcevt.scmpevt),&hevt);
                  } else { /* be sure that the structure is empty */
                    memset(&hevt,0,sizeof(hyperCmpEvent));
                  }
                  /* for MC events we need an output even if it is not passing the filter */
                  FillHyperMCParticles(&sbur,&supermcevt,&hevt);
                  rwHyperCmpEvent(hcmpoutn,&hevt,&hbur);
                }
              }

	      /*User filtering required */
	      if (cmpReader_options.hcmp_out_name) {
		userStatus=USER_SUPERMCFILTER(&sbur,&supermcevt,&hevt);
		if(userStatus>0) {
		  sprintf(buf,"Error code %d returned from user_superMcFilter",
			  userStatus);
		  cmpwarn(buf);
		} else if (hcmpout) {
		  if(userStatus==0) { /*Be sure the event is empty */
		    memset(&hevt,0,sizeof(hyperCmpEvent));
		  }
		  FillHyperMCParticles(&sbur,&supermcevt,&hevt);
                  rwHyperCmpEvent(hcmpout,&hevt,&hbur);
		}
	      }
	      
              break;
             
              /*----------> Hyper event structure <----------*/
            case HCMPEVENT_STRUCT:
              zread=rwHyperCmpEvent(cmpin,&hevt,&hbur);
              if(!zread)                            /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
                }
              if(hbur.BadB.Skip) break; /* Dont analyse events from bad bursts */
              if( (RMINflag) && (hbur.nrun < RunMin) ) break;
              if( (RMAXflag) && (hbur.nrun > RunMax) ) break;

              //hyperCmpAnaEvent(&hbur,&hevt);            /* Compact analysis / event */
              nhcmp++;
              readEvents++;
              userStatus = USER_HYPERCMPEVENT(&hbur,&hevt);
              if(userStatus>0)                            /* check user return for error */
                {      
                  sprintf(buf,"Error code %d returned from user_hyperCmpEvent",
                          userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);

              /*************** user filtering required */
              if (cmpReader_options.hcmp_out_name) {
                userStatus=USER_HYPERCMPFILTER(&hbur,&hevt);
                fflush(NULL);
                if(userStatus>0){                       /* examine filter return for error */
                  sprintf(buf,"Error code %d returned from user_hyperCmpFilter", userStatus);
                  cmpwarn(buf);
                } else if(userStatus<0) { /* the user want this event on disk */
                  if (hcmpout) rwHyperCmpEvent(hcmpout,&hevt,&hbur);
                }
              }
              break;             
             
              /*----------> End-of-burst structure <----------*/

            case EOB_STRUCT:                      /* End-of-burst structure */
              zread=rwEndofBurst(cmpin,&eob);     /* read in end of burst */
              if(!zread)                          /* skip burst if problem reading */
                {
                  JumpToNextFile=1; 
                  break;
                }
              if( (RMINflag) && (bur.nrun < RunMin) ) break;
              if( (RMAXflag) && (bur.nrun > RunMax) ) break;
              /* 
                 ~15-04-98 BH+IWS: reset to zero unused words (for compression)
              */
              /* 
                 18-05-98 Remove from 4.1 
              */
              /*
                lpoint=&(eob.ProcError[0])+eob.NProcError;
                memset((void *)lpoint,0,sizeof(int)*(MAX_PROCERROR-eob.NProcError));

                lpoint=&(eob.TrigWord[0])+eob.NEventList;
                memset((void *)lpoint,0,sizeof(int)*(MAX_TRIGHIST-eob.NEventList));

                lpoint=&(eob.TimeStamp[0])+eob.NEventList;
                memset((void *)lpoint,0,sizeof(int)*(MAX_TRIGHIST-eob.NEventList));
              */

              readEoB++;                       /* increase number of "end of burst" */
#ifdef COMPACT_F77
              if(readEoB==1)                /* run check on first event read */
                {         
                  cmpmsg("main: First EoB found - running C/F77 checks");
                  chkF77EoB(&eob);
                }
#endif

              cmpAnaEoB(&bur,&eob);
	      if (cmpReader_options.calcbb) { /*Bad bursts calculations necessary - call ana routines*/
		MUV_EOB(&bur,&eob); /* MUV bad burst */
	      }
              userStatus=USER_EOB(&bur,&eob);
              if(userStatus>0)              /* examine filter return for error */
                {
                  sprintf(buf,"Error code %d returned from user_eob",userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);
              if(CHEATflag){
                cmp2scmpEndofBurst(&bur,&eob,&seob);
		if(cmpReader_options.abcog) { /* Need alpha, beta, COG production*/
		  abc_SUPEREOB(&sbur,&seob);
		}
                userStatus=USER_SUPEREOB(&sbur,&seob);
                fflush(NULL);
                if(userStatus>0){
                  sprintf(buf,"Error code %d returned from user_seob",userStatus);
                  cmpwarn(buf);
                }
                break;
              }

              if (cmpReader_options.cmp_out) {
                if (cmpout){
                  if((EMPTYflag)) {
                    if(Empty&E_LIST) {
                      bur.NTrigWord      = eob.NTrigWord      = 0;
                      bur.NTimeStamp     = eob.NTimeStamp     = 0;
                      bur.NDCHDecError   = eob.NDCHDecError   = 0;
                      bur.NLKRHACDecError= eob.NLKRHACDecError= 0;
                      bur.NChargedInfo   = eob.NChargedInfo   = 0;
                      bur.NNeutralInfo   = eob.NNeutralInfo   = 0;
                    }
                    if(Empty&E_EOB) {
                      /* Empty EndofBurst structure */
                      memset(&eob,0,sizeof(EndofBurst));
                    }
                  }
                  /* Update event list for EoB structure */
                  BuildEoBList(&bur, &eob, LastEventOut);
                  rwEndofBurst(cmpout,&eob);  /* and write out end of burst struct */
                }
              }
              if (cmpReader_options.mc_out_name) {
		if (mcout) rwEndofBurst(mcout,&eob);
	      }
	      if (cmpReader_options.scmp_out || cmpReader_options.smc_out_name) {
                cmp2scmpEndofBurst(&bur, &eob, &seob);  /* copy superEoB variables */
                if (scmpout) rwSuperEndofBurst(scmpout,&seob);      /* write out superBurst */
		if (smcout)  rwSuperEndofBurst(smcout,&seob);
		if(cmpReader_options.abcog) { /* Need alpha, beta, COG production*/
		  abc_SUPEREOB(&sbur,&seob);
		}
              }
              

              if(KEOflag)                   /* check for Ke3 output file */
                rwEndofBurst(ke3out,&eob);  /* and write out end of burst struct */
              if(KMUOflag)                  /* check for Kmu3 output file */
                rwEndofBurst(kmu3out,&eob); /* and write out end of burst struct */

#ifdef BFLAG
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
#endif
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
	      if(cmpReader_options.abcog) { /* Need alpha, beta, COG production*/
		abc_SUPEREOB(&sbur,&seob);
	      }
              userStatus=USER_SUPEREOB(&sbur,&seob);
              if(userStatus>0)                        /* examine filter return for error */
                {
                  sprintf(buf,"Error code %d returned from user_seob",userStatus);
                  cmpwarn(buf);
                }
              fflush(NULL);

              if (cmpReader_options.scmp_out) {
                if (scmpout) rwSuperEndofBurst(scmpout,&seob);      /* write out superBurst */
                if (scmpoutc) rwSuperEndofBurst(scmpoutc,&seob);    /* write out superBurst (charged stream) */
                if (scmpoutn) rwSuperEndofBurst(scmpoutn,&seob);    /* write out superBurst (neutral stream) */
              }
              if (cmpReader_options.smc_out) {
		if (smcout)   rwSuperEndofBurst(smcout,&seob);      /* write out superBurst (MC) */
              }
	      /* no hyperCOmPACT write because hyperCOmPACT has not an EOB*/

#ifdef BFLAG
              /* if fixed size output files are required check size now */
              if(Bflag)
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
                }
#endif
              break;
             


              /*----------> Unknown structure <----------*/

            default:
              if(readBursts)
                {
                  JumpToNextFile=1;
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
      if (cmpin) closeCompact(cmpin);            /* close input file */
      while( NotOpen) {
        if(current->next) {           /* is there another file in the list */
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
                  fclose(fp);
                }
            } /* End test on iFile%statusUpd */
        }
        else
          {
            LastFile=TRUE;
          }
/* Marco.Clemencic@cern.ch === This part is not used since ages!!! (-s -e -n -d options) === */
/*           else  */
/*             { */
/*               if(maxBursts==0 || readBursts<maxBursts)   /\* check not over maxburst *\/ */
/*                 { */
/*                   cmpin=openNextBurst(rundir,&cBurst,eBurst); /\* open next file *\/ */
/*                   if(!cmpin)                  /\* check another file was found *\/ */
/*                     LastFile=TRUE; */
/*                 } */
/*             } */
          /* End test on type on input file */
          if(cmpin==(CMPstream *)0)
            {
              if(LastFile) break;
              log_skipped_file(current->fname);
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
    fclose(fp);
  }
  if(cmpReader_options.abcog) {/* Alpha, beta and COG production - EXIT */
    abc_EXIT;
  }
  userStatus=USER_EXIT;              /* call user termination  */
  fflush(NULL);
  sql_close();
   
  if(userStatus>0)                   /* examine filter return for error */
    { 
      sprintf(buf,"Error code %d returned from user_exit",userStatus);
      cmpwarn(buf);
    }
  if(KEOflag)                        /* check for Ke3 output file */
    closeCompact(ke3out);            /* and then close it */
  if(KMUOflag)                       /* check for Ke3 output file */
    closeCompact(kmu3out);           /* and then close it */

  /* check for COmPACT output files */
  if (cmpout) closeCompact(cmpout);
  /* check for superCOmPACT output files */
  if (scmpout) closeCompact(scmpout);
  if (scmpoutc) closeCompact(scmpoutc);
  if (scmpoutn) closeCompact(scmpoutn);
  /* check for hyperCOmPACT output files */
  if (hcmpout) closeCompact(hcmpout);
  if (hcmpoutc) closeCompact(hcmpoutc);
  if (hcmpoutn) closeCompact(hcmpoutn);
  /* check for COmPACT MC output files */
  if (mcout) closeCompact(mcout);
  /* check for superCOmPACT MC output files */
  if (smcout) closeCompact(smcout);
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
  sprintf(buf,"  Number of HCMP events read :%9d",nhcmp);
  cmpmsg(buf);
  cmpmsg("******************************************");
#if defined(COMPACT_F77) && defined(ALPHA)
  /* This call is required by the F77 runtime library on Alphas */
  for_rtl_finish_();   /* shutdown F77 IO on Alpha */
#endif
}
