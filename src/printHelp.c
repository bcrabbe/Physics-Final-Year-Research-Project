/********************************************************/
/* COmPACT routine: printHelp(void)                     */
/*                                                      */
/* This routine prints the command li e options for the */
/* COmPACT reader program. It can be invoked by either  */
/* the command line option `-h' or by incompatible      */
/* command line options being chosen.                   */
/*                                                      */
/*                                         RWM 19/6/97  */
/********************************************************/

#include "cmpio.h"

#ifdef __cplusplus
extern "C" 
#endif
void printHelp(void) {
  cmpmsg("Compact command arguments:");
/* 
   cmpmsg("   -s   A:B      start at run A, burst B");
   cmpmsg("   -e   C:D      end at run A, burst B");
   cmpmsg("   -n   #        run over a maximum of # bursts");
   cmpmsg("   -d   <dir>    set data directory");
*/
  cmpmsg("   -i   <file>      declare input file name (replaces s/e/n options)");
  cmpmsg("   -l   <file>      give input file list (replaces s/e/n/i options)");
/*
  cmpmsg("   -x   \"args\"   run \"runcompact\" from within compact");
  cmpmsg("                 use \"runcompact -d\" and don't use i/l options");
*/
  cmpmsg("   -co  <file>      declare COmPACT output file name");
  cmpmsg("   -ko  <file>      declare Ke3 output file name");
  cmpmsg("   -kmo <file>      declare Kmu3 output file name");
  cmpmsg("   -mo  <file>      declare MC output file name");
  cmpmsg("   -so  <file>      declare SuperCOmPACT output file name (user filter)");
  cmpmsg("   -soc <file>      declare SuperCOmPACT output file name (charged filter)");
  cmpmsg("   -son <file>      declare SuperCOmPACT output file name (neutral filter)");
  cmpmsg("   -ho  <file>      declare HyperCOmPACT output file name (user filter)");
  cmpmsg("   -hoc <file>      declare HyperCOmPACT output file name (charged filter from scmp)");
  cmpmsg("   -hon <file>      declare HyperCOmPACT output file name (neutral filter from scmp)");
  cmpmsg("   -scprod-cuts     set default for looser filter cuts (SCProd)");
  cmpmsg("   -nokl2           disable kl2 filter in the neutral split (enabled by default)");
  cmpmsg("   -mcver           Set MC version to ...");
  cmpmsg("                    1: - use mk=0.9*mk+0.0493377 for a,b corrections");
  cmpmsg("                    2: - use mk=0.9*mk+0.0493377-0.000035 for a,b corrections");
  cmpmsg("                    not set: - use mk=0.9*mk+0.0493377-0.000035 for a,b corrections");
  cmpmsg("   -smo <file>      declare MC output file name");
  cmpmsg("   -sl  #           set SuperCOmPACT output level");
  cmpmsg("                    default (no option): epsilon standard super-compact");
  cmpmsg("                    #=16+iflag (iflag=1,2,4) for 2pi0, 3pi0, 2 gamma respectively");
  cmpmsg("                    #=32 for Dalitz summary");
  cmpmsg("                    #=64 for rare decay summary");
  cmpmsg("                    #=128 for four-trackssummary");
  cmpmsg("                    #=a combination of 16+iflag, 32, 64, 128 is also accepted");
  cmpmsg("   -srw             rewrite SuperCOmPACT integer event structure");
  cmpmsg("                    from float version");
#ifdef BFLAG
  cmpmsg("   -b   #           produce o/p files of fixed size");
  cmpmsg("                    # can be nnn for bytes, nnnMB or nnnGB");
#endif
  cmpmsg("   -db              access cdb (compact database)"); 
  cmpmsg("   -ndb             do not access cdb (compact database) This is the default"); 
  cmpmsg("   -cheat           transform compact events into super-compact");
  cmpmsg("   -skip-kab-rec    turn off Kabes reconstruction (when reading cmp events)");
#ifdef FIX_LKR_TIME
  cmpmsg("   -fix-lkr-toff #  fix LKR Clusters time offset: 1->use run 12163, 2->use right run");
  cmpmsg("                                                  3->use muv instead of lkr-tag");
#endif
  cmpmsg("   -rmin <RMIN>     Only analysis bursts for which run>=RMIN");
  cmpmsg("   -rmax <RMAX>     Only analysis bursts for which run<=RMAX");
  cmpmsg("   -nevt <NEVT>     Only analysis NEVT events and then stops");
  cmpmsg("   -empty #         Empty Burst, EoB, evt lists structures:");
  cmpmsg("                    to save disk space when filtering events");
  cmpmsg("                    #=0: empty burst, evt lists, EoB");
  cmpmsg("                    #=1: empty burst");
  cmpmsg("                    #=2: empty evt lists");
  cmpmsg("                    #=4: empty EoB");
  cmpmsg("                    #=6,3,.. a combination of above is also accepted");
  cmpmsg("  -string <string>  Pass an arbitrary string. This string is available as global");
  cmpmsg("  -h                print this help message");
}


/*
 *
 *
 * $Log: printHelp.c,v $
 * Revision 1.9  2006/01/17 13:25:14  venelin
 * Added integer option "-mcver"
 * 	1: - use mk=0.9*mk+0.0493377 for a,b corrections
 * 	2: - use mk=0.9*mk+0.0493377-0.000035 for a,b corrections
 * 	not set: - use mk=0.9*mk+0.0493377-0.000035 for a,b corrections
 *
 * Revision 1.8  2004/10/12 12:08:20  clemenci
 * Added the fix for a 3rd type of lkr time offset bug
 *
 * Revision 1.7  2004/10/07 08:28:53  venelin
 * Added switch for kl2 filter in the neutral split
 *
 * Revision 1.6  2004/09/21 13:35:19  clemenci
 * Fixed the bug on LKR time offset in SCompact.
 * The code to correct the wrong value is in the reader and has to be activated in the Makefile.
 *
 * Revision 1.5  2004/07/21 17:54:17  clemenci
 * removed (commented) the -b option and implemented the writing of EoB structures, which was not done yesterday
 *
 * Revision 1.4  2004/07/21 01:02:57  clemenci
 * redesigned the code to write cmp, scmp and hcmp (filtering and conversions)
 *
 * Revision 1.3  2004/07/19 17:08:00  clemenci
 * Kabes reconstruction can be turned off by a command line option or an env
 * variable
 *
 * Revision 1.2  2004/07/10 22:43:25  clemenci
 * fixed c++ compilation error
 *
 * Revision 1.1  2004/06/23 15:31:19  clemenci
 * moved printHelp.c from libreader to the reader main tree (why was it there?)
 *
 * Revision 1.3  2003/10/31 12:32:30  andrew
 * Added the -string option to main
 *
 * An arbitrary string can be passed to compact which is
 * saved in a global variable gString (C), COMMON/GSTRING/GSTRING (FORTRAN)
 *
 *
 * made -ndb the default. For people needin the compact database the -db option
 * was created
 *
 *
 *
 */
