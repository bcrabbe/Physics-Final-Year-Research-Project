#
# Makefile for generating the COmPACT reader - requires GNUmake (gmake)
#
# Get OS type from uname shell command
SHELL   = /bin/sh
OSTYPE = $(shell uname)
HOSTNAME = $(shell hostname)
#-----------------------User Selections-----------------------------
# Uncomment the following line for F77 mode compilation
#F77DEF  = -DCOMPACT_F77
# set the following line to DEBUG=yes for debug compilation of reader
DEBUG = yes
# Is zlib installed on your system?
ZLIBSYS = no
# To use rfio, in particular on the PC farm, set USE_RFIO=yes
USE_RFIO=yes
#Do you want to use the libraries stored on public afs ?
# If not, (PUBLIC_AFS = no) libraries in subdirectories from compact 
# is assumed	
# Change following line to use "local" libraries
PUBLIC_AFS = yes

# Uncomment the following line if you want to enable the run-time fix
# for LKR time offsets
# FIX_LKR_TIME = yes

#
# If you want to use Root instead of Paw set USE_ROOT to yes
#
#
USE_ROOT=no

.SUFFIXES =
.SUFFIXES = .o .c .cc .C .h .f .F .a .d .Fd

# CFortran include directory
CFORINC  = -I/cern/pro/include/cfortran

# Setup ZLIB variables to point to correct location

AR       = ar
RM       = rm -f
MV       = mv
TOUCH    = touch
RANLIB   = ranlib
ifeq ($(OSTYPE),Linux)
MKDIR    = mkdir -p -v
else
MKDIR    = mkdir -p
endif

VERSION  = 7.3
SUBVERSION = .1

ifeq ($(USE_RFIO),yes)
ZLIBTYPE=-rfio
else
ZLIBTYPE=
endif

ifeq ($(OSTYPE),Linux)
#Checking which version of Linux we have
ISSLC3 = $(shell cat /etc/redhat-release |grep -ic "Scientific Linux CERN release 3")
ISSLC4 = $(shell cat /etc/redhat-release |grep -ic "Scientific Linux CERN SLC release 4")
ISSLC5 = $(shell cat /etc/redhat-release |grep -ic "Scientific Linux CERN SLC release 5")
ISSLC6 = $(shell cat /etc/redhat-release |grep -ic "Scientific Linux CERN SLC release 6")

ifeq ($(ISSLC6),1)
LINVER = -slc6
#GFORTRAN needed on lxplus/lxbatch since CERNLIB is compiled with gfortran
GFORLIBS = -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4/32 -lgfortran
#32 bit libshift library not linked to libshift.so, use the binary
LIBSHIFT = -L. $(shell ls /usr/lib/libshift.* | head -1)
else
ifeq ($(ISSLC5),1)
LINVER = -slc5
#GFORTRAN needed on lxplus/lxbatch since CERNLIB is compiled with gfortran
GFORLIBS = -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/32 -lgfortran
#32 bit libshift library not linked to libshift.so, use the binary
LIBSHIFT = $(shell ls /usr/lib/libshift.* | head -1)
else
ifeq ($(ISSLC3),1)
LINVER = -slc3
else 
ifeq ($(ISSLC4),1)
LINVER = -slc4
else
LINVER = -rh7.3
endif #($(ISSLC4),1)
endif #ifeq ($(ISSLC3),1)
LIBSHIFT = -lshift
GFORLIBS = 
endif #ifeq ($(ISSLC5),1)
endif #ifeq ($(ISSLC6),1)
SQLITETYPE =
else
SQLITETYPE = -$(OSTYPE)
endif

TOPDIR  = .
INCDIR  = $(TOPDIR)/inc
SRCDIR  = $(TOPDIR)/src
UINCDIR = $(TOPDIR)/userinc
USRCDIR = $(TOPDIR)/usersrc
USRADIR = $(TOPDIR)/userana
OBJDIR  = $(TOPDIR)/obj
DEPDIR  = $(TOPDIR)/depends
DOCDIR  = $(TOPDIR)/doc
ifeq ($(USE_ROOT),yes)
ROOTCFLAGS= $(shell $(ROOTSYS)/bin/root-config --cflags) -DUSE_ROOT
ROOTLIBS=$(shell $(ROOTSYS)/bin/root-config --libs)
else
ROOTCFLAGS=
ROOTLIBS=
endif

4SQRT =
# to run axna48mon, uncomment the 3 follwing lines
#ifeq ($(OSTYPE),OSF1)
#    4SQRT = -lm_4sqrt
#endif
ifeq ($(PUBLIC_AFS),yes)
   CINCDIR = /afs/cern.ch/na48/offline2/compact/compact-$(VERSION)/compact/lib/inc
   ZINCDIR = /afs/cern.ch/na48/offline2/compact/compact-$(VERSION)/compact/zlib-1.1.4
   RINCDIR = /afs/cern.ch/na48/offline2/compact/compact-$(VERSION)/compact/rlib/inc
   RAINCDIR = /afs/cern.ch/na48/offline2/compact/compact-$(VERSION)/compact/rlib/anainc
   ifeq ($(OSTYPE),OSF1)
      LIBDIR = /afs/cern.ch/na48/offline2/compact/.alpha_osf20/lib
   endif
   ifeq ($(OSTYPE),HP-UX)
      LIBDIR = /afs/cern.ch/na48/offline2/compact/.hp700_ux90/lib
   endif
   ifeq ($(OSTYPE),SunOS)
      LIBDIR = /afs/cern.ch/na48/offline2/compact/.sun4m_53/lib
   endif
   ifeq ($(OSTYPE),Linux)
      LIBDIR = /afs/cern.ch/na48/offline2/compact/.linux$(LINVER)/lib
   endif
# for the database
   DBINCDIR = /afs/cern.ch/user/c/cdbase1/public/cdbase/code/inc
ifeq ($(findstring cerncs2,$(HOSTNAME)),cerncs2)
   DBLIBDIR = /afs/cern.ch/user/c/cdbase1/public/cdbase/code/.sun4x_55/lib   
else
   DBLIBDIR = /afs/cern.ch/user/c/cdbase1/public/cdbase/code/lib
endif
   DBLIB = -lcdb

   ZLIB = -lz$(ZLIBTYPE)
   SQLITEINCDIR = /afs/cern.ch/na48/offline2/compact/compact-$(VERSION)/compact/sqlite/include
   SQLITELIB = -lsqlite$(SQLITETYPE)
ifeq ($(DEBUG),yes)
   RLIB = -lreader-gdb.$(VERSION)$(SUBVERSION)
   CLIB =  $(LIBSHIFT) -lcompact-gdb.$(VERSION)$(SUBVERSION)
else   
   RLIB = -lreader.$(VERSION)$(SUBVERSION)
   CLIB = $(LIBSHIFT) -lcompact.$(VERSION)$(SUBVERSION)
endif   
   LDIRS    = -L$(LIBDIR) -L$(DBLIBDIR)

# If libraries not taken from afs
else
   WORKDIR = $(shell pwd)
   LIBDIR  = $(WORKDIR)/../lib
   CINCDIR = $(LIBDIR)/inc
   ZLIBDIR = $(WORKDIR)/../zlib-1.1.4
   SQLITEDIR = $(WORKDIR)/../sqlite
   SQLITEINCDIR = $(WORKDIR)/../sqlite/include
   ZINCDIR = $(ZLIBDIR)
   RLIBDIR = $(WORKDIR)/../rlib
   RINCDIR = $(RLIBDIR)/inc
   RAINCDIR = $(RLIBDIR)/anainc
   ZLIB = -lz$(ZLIBTYPE)-$(OSTYPE)
   RLIB = -lreader.$(VERSION)
   CLIB = $(LIBSHIFT) -lcompact.$(VERSION)
   SQLITELIB = -lsqlite$(SQLITETYPE)

# for the database
   DBINCDIR = $(WORKDIR)/../cdbase/code/inc
   DBLIBDIR = $(WORKDIR)/../cdbase/code/lib
   DBLIB = -lcdb

   LDIRS = -L$(LIBDIR) -L$(RLIBDIR) -L$(ZLIBDIR) -L$(DBLIBDIR) -L$(SQLITEDIR)/lib
endif

LIBS     =   $(CLIB) $(SYSLIBS) $(RLIB) $(ZLIB) $(4SQRT) $(DBLIB) -lm $(SQLITELIB)
DEFINES  = 
CPPFLAGS = -I. -I$(UINCDIR) -I$(CINCDIR) -I$(RINCDIR) -I$(RAINCDIR)  -I$(ZINCDIR) \
           -I$(DBINCDIR) -I$(SQLITEINCDIR) $(ROOTCFLAGS) $(F77DEF) $(SYSDEFS) $(DEFINES)
CPPFLAGS += -I/usr/local/include/shift -I/usr/include/shift $(CFORINC)
#CERNLIBS = `/cern/pro/bin/cernlib mathlib packlib`
CERNLIBS = -L/cern/pro/lib -lmathlib -lpacklib -lnsl -lcrypt -ldl

ifeq ($(FIX_LKR_TIME),yes)
DEFINES += -DFIX_LKR_TIME
endif

# Set this to be the GNU C compiler
ifeq ($(USE_ROOT),yes)
CC     = g++  
else
CC     = gcc 
endif
# Set any GCC compiler flags here.
ifeq ($(DEBUG),yes)
CFLAGS = -m32 -DCOMPACT7 -std=c99 $(ROOTCFLAGS)
else
CFLAGS = -m32 -O -DCOMPACT7 -std=c99 $(ROOTCFLAGS)
endif
# Set this to be your F77 compiler
FC     = g77
#FFLAGS = -m32 -w -ffixed-line-length-none -fno-second-underscore -g -DCOMPACT7 -lg2c -lm
#----------------------------------------------------------------------
#---------------------System Specific configurations-------------------
#----------------------------------------------------------------------
# The option needed to enable extended source mode (>72 characters/line)
# must be included in FFLAGS
# CPP must be set to the location of the GNU C pre-processor. The
# current locations work for 'standard' CERN machines.
#----------------------------Linux-------------------------------------
ifeq ($(OSTYPE),Linux)
# Following 3 lines are compiler options when using Portland Group 
# fortran compiler: pgf77
# Can be used on pcna48df3 for instance
# Uncomment these 3 lines and comment the other 3
#F77LIBS = -L /usr/local/pgi/linux86/lib -lpgftnrtl -lpgc -lf2c -lm -lc
#FC    = pgf77 
#FFLAGS = -O -Mextend -f77libs 

# Following 3 lines are compiler options when using g77 
# fortran compiler: f77 is in fact calling g77.
# Can be used on pcna48df1 for instance
# Uncomment these 3 lines and comment the other 3 (above)
# 
#F77LIBS = 
#FC    = f77 
# for offline PC farm, use g77
#FC = g77
ifeq ($(DEBUG),yes)
FFLAGS =  -m32 -w -ffixed-line-length-none -fno-second-underscore -g -DCOMPACT7
else
FFLAGS = -m32 -O -w -ffixed-line-length-none -fno-second-underscore -DCOMPACT7
endif
# End of compiler options for Linux
CFORINC := -Df2cFortran $(CFORINC)
#CLDFLAGS  = -Wl,-t
CLDFLAGS  = -Wl,-V
# IWS (following Steffen's advice) 07-01-98 CPP = /lib/cpp
CPP = gcc -E
# redhat 5.1 
#SYSLIBS = -L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.90.27
# redhat 6.1 
SYSLIBS = -L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66 $(GFORLIBS)
endif
#----------------------------DEC Alpha---------------------------------
# For Alpha require "-DALPHA" and "-DLONG64", for other 64 bit systems
# with 64bit long integers "-DLONG64" is required due to a bug in XDR
ifeq ($(OSTYPE),OSF1)
F77LIBS = -lfor -lots
ifeq ($(DEBUG),yes)
FFLAGS = -extend_source -fpe2 -g -DCOMPACT7
else
FFLAGS = -extend_source -fpe2  -DCOMPACT7
endif
SYSDEFS = -DALPHA -DLONG64
#CPP=/usr/local/lib/gcc-lib/alpha-dec-osf4.0b/egcs-2.91.60/cpp
#CPP     = /usr/local/lib/gcc-lib/alpha-dec-osf3.2/2.7.2.1/cpp
CPP      =  cpp
endif
#------------------------------Sun-------------------------------------
# Suns need '-lnsl' library to access XDR functions and socket to use
# CERNlibs
ifeq ($(OSTYPE),SunOS)
F77LIBS   = -L/opt/SUNWspro/SC4.0/lib -lF77 -lM77 -lsunmath
# Turn off stupid linker warning using -Wl,-t
CLDFLAGS  = -Wl,-t
# 25-04-2001 switch to g77
ifeq ($(DEBUG),yes)
FFLAGS = -g -w -ffixed-line-length-none -fno-second-underscore -DCOMPACT7
else
FFLAGS = -O -w -ffixed-line-length-none -fno-second-underscore -DCOMPACT7
endif
#ifeq ($(DEBUG),yes)
#FFLAGS    = -e -g
#else
#FFLAGS    = -e  
#endif
SYSLIBS   = -lnsl
CERNLIBS += -lsocket -lposix4
# cpp on cerncs2
ifeq ($(findstring cerncs2,$(HOSTNAME)),cerncs2)
CPP       = /usr/local/lib/gcc-lib/sparc-sun-solaris2.3/2.7.2/cpp
else
# cpp on na48 work group server
CPP       = /usr/local/lib/gcc-lib/sparc-sun-solaris2.6/egcs-2.91.66/cpp
endif
endif
#-------------------------------HP-------------------------------------
ifeq ($(OSTYPE),HP-UX)
F77LIBS = -lf
# Have to use f2c for CFortran because of "+ppe" option used for CERNLIB
CFORINC := -Df2cFortran $(CFORINC)
# 25-04-2001 switch to g77
ifeq ($(DEBUG),yes)
FFLAGS = -g -w -ffixed-line-length-none -fno-second-underscore -DCOMPACT7
else
FFLAGS = -O -w -ffixed-line-length-none -fno-second-underscore -DCOMPACT7
endif
#ifeq ($(DEBUG),yes)
#FFLAGS  = +es +ppu -g
#else
#FFLAGS  = +es +ppu
#endif
SYSLIBS   = -L/lib/pa1.1
#25-04-2001 switch to cpp
#CPP   = /usr/local/lib/gcc-lib/hppa1.1-hp-hpux10.20/2.7.2.2/cpp
CPP = gcc -E
endif
#-----------------------------AIX--------------------------------------
ifeq ($(OSTYPE),AIX)
F77LIBS = 
FC = g77
ifeq ($(DEBUG),yes)
#FFLAGS = -k -P -g -DCOMPACT7
FFLAGS = -g -w -ffixed-line-length-none -fno-second-underscore  -DCOMPACT7
else
#FFLAGS = -k -P -DCOMPACT7
FFLAGS = -w -ffixed-line-length-none -fno-second-underscore  -DCOMPACT7
endif
#CPP = /usr/local/lib/gcc-lib/rs6000-ibm-aix4.1.4.0/2.7.2.2/cpp
CPP = gcc -E
endif
#-----------------------------SGI-----------------------------------------
ifeq ($(OSTYPE),IRIX64)
CC     = cc
CFLAGS += -32 -common -Olimit 1070 -Umips -DCOMPACT7
#FC     = f77
F77LIBS = -lF77 -lI77 -lisam
ifeq ($(DEBUG),yes)
FFLAGS = -O -static -32 -extend_source -g -DCOMPACT7
else
FFLAGS = -O -static -32 -extend_source -DCOMPACT7
endif

#FFLAGS = -O -static -n32 -extend_source

endif
#----------------------------------------------------------------------

# Define search paths for the various files
vpath %.c $(SRCDIR):$(USRCDIR):$(USRADIR)
vpath %.F $(USRCDIR):$(USRADIR)
vpath %.h $(INCDIR):$(UINCDIR):$(DBINCDIR)
vpath %.o $(OBJDIR)
vpath %.d $(DEPDIR)

# Define main COmPACT sources and objects
SRCS = compact_main.c \
       printHelp.c \
       user.c \

OBJS = $(SRCS:.c=.o)

# Define user sources and objects
UFSRCS = fuser_init.F \
         fuser_burst.F \
         fuser_superBurst.F \
         fuser_ke3Event.F \
         fuser_kmu3Event.F \
         fuser_mcEvent.F \
         fuser_superCmpEvent.F \
         fuser_superMcEvent.F \
         fuser_cmpFilter.F \
         fuser_ke3Filter.F \
         fuser_kmu3Filter.F \
         fuser_mcFilter.F \
         fuser_superCmpFilter.F \
         fuser_hyperCmpFilter.F \
         fuser_superMcFilter.F \
         fuser_eob.F \
         fuser_superEob.F \
      	 fuser_cmpEvent.F \
         fuser_exit.F \
		 fuser_hyperBurst.F\
		 fuser_hyperCmpEvent.F

UFOBJS = $(UFSRCS:.F=.o)

UCSRCS = user_init.c \
         user_burst.c \
         user_superBurst.c \
         user_cmpEvent.c \
       	user_ke3Event.c \
         user_kmu3Event.c \
         user_mcEvent.c \
         user_superCmpEvent.c \
         user_superMcEvent.c \
         user_cmpFilter.c \
         user_ke3Filter.c \
         user_kmu3Filter.c \
         user_mcFilter.c \
         user_superCmpFilter.c \
         user_hyperCmpFilter.c \
         user_superMcFilter.c \
         user_eob.c \
         user_superEob.c \
         user_exit.c \
		 user_hyperBurst.c\
		 user_hyperCmpEvent.c

UCOBJS = $(UCSRCS:.c=.o)

UCASRCS = \

UCAOBJS = $(UCASRCS:.c=.o)

UFASRCS =     	Blue_Tack.F\

UFAOBJS = $(UFASRCS:.F=.o)

UAOBJS = $(UFAOBJS) $(UCAOBJS) 

DEPENDS   = $(SRCS:%.c=$(DEPDIR)/%.d)
UFDEPENDS = $(UFSRCS:%.F=$(DEPDIR)/%.d)
UCDEPENDS = $(UCSRCS:%.c=$(DEPDIR)/%.d)
UFADEPENDS = $(UFASRCS:%.F=$(DEPDIR)/%.d)
UCADEPENDS = $(UCASRCS:%.c=$(DEPDIR)/%.d) 

ifeq ($(F77DEF),-DCOMPACT_F77)
 UOBJS=$(UFOBJS)
 LASTCMP=.last.f77.compile.h
 OTHERCMP=.last.c.compile.h
else
 UOBJS=$(UCOBJS)
 LASTCMP=.last.c.compile.h
 OTHERCMP=.last.f77.compile.h
endif

default:
	gmake compact

compact: $(OBJS) $(UOBJS) $(UAOBJS)
	$(FC) $(CFLAGS) $(CLDFLAGS) -o compact \
	$(OBJS:%=$(OBJDIR)/%) \
	$(UOBJS:%=$(OBJDIR)/%) \
	$(UAOBJS:%=$(OBJDIR)/%) \
	$(LDIRS) $(LIBS) $(CERNLIBS) $(F77LIBS) $(LIBS) $(ROOTLIBS)

clean:
	$(RM) -r $(OBJDIR) $(DEPDIR)
	$(RM) compact
	$(RM) .last.c.compile.h
	$(RM) .last.f77.compile.h

purge:
	$(RM) $(SRCDIR)/*~
	$(RM) $(USRCDIR)/*~
	$(RM) $(USRADIR)/*~
	$(RM) $(UINCDIR)/*~

ifneq ("$(findstring clean,$(MAKECMDGOALS))","clean")
ifneq ("$(findstring purge,$(MAKECMDGOALS))","purge")
include $(DEPENDS) $(UCDEPENDS) $(UFDEPENDS) $(UCADEPENDS) $(UFADEPENDS)
endif
endif

$(DEPDIR)/%.d: %.c
	@ [ ! -d $(DEPDIR) ] && $(MKDIR) $(DEPDIR) || true
	$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $(CFORINC) $< \
                      | sed s/$*.o/$*.o\ $*.d/ > $@'

$(DEPDIR)/%.d: %.F
	@ [ ! -d $(DEPDIR) ] && $(MKDIR) $(DEPDIR) || true
	$(SHELL) -ec '$(FC) -MM $(CPPFLAGS) $< \
                      | sed s/$(*F).F.o/$*.o\ $*.d/ > $@'

$(SRCS) $(UCSRCS) $(UFSRCS) $(UCASRCS) $(UFASRCS):

$(OBJS):
	@ [ ! -d $(OBJDIR) ] && $(MKDIR) $(OBJDIR) || true
	$(CC) $(CFLAGS) $(CFORINC) $(CPPFLAGS) -c $< -o $(OBJDIR)/$(*F).o

$(UCOBJS):
	@ [ ! -d $(OBJDIR) ] && $(MKDIR) $(OBJDIR) || true
	$(CC) $(CFLAGS) $(CFORINC) $(CPPFLAGS) -c $< -o $(OBJDIR)/$(*F).o

$(UFOBJS):
	@ [ ! -d $(OBJDIR) ] && $(MKDIR) $(OBJDIR) || true
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)/$(*F).o

$(UCAOBJS):
	@ [ ! -d $(OBJDIR) ] && $(MKDIR) $(OBJDIR) || true
	$(CC) $(CFLAGS) $(CFORINC) $(CPPFLAGS) -c $< -o $(OBJDIR)/$(*F).o 

$(UFAOBJS):
	@ [ ! -d $(OBJDIR) ] && $(MKDIR) $(OBJDIR) || true
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)/$(*F).o
# For Linux	mv *.o $(OBJDIR)/.

### This targets are needed to correctly recompile when changeing from C to F and back
#$(OBJS): $(LASTCMP)
$(LASTCMP):
	$(RM) $(OTHERCMP)
	$(TOUCH) $(LASTCMP)

# This dummy dependency forces others to be recompiled everytime. This
# is used to ensure the COmPACT library is always uptodate.
FORCE:

# $Log: Makefile,v $
# Revision 1.35  2010/07/16 12:25:21  venelin
# Disable generation of functions time profile - safe of 1MB of extra disk-space per compact process
# 	(forgotten from previous debugging purposes)
#
# Revision 1.34  2009/11/29 17:51:36  venelin
# Make compact compile on SLC5
# Changes affect the usage of:
# 	libncurses  -> needed by sqlite
# 	libshift    -> compact
# 	libgfortran -> CERNLIB on SLC5 compiled with gfortran instead of g77
#
# Revision 1.33  2007/12/17 17:05:31  venelin
# Changing the version of compact to 7.3.0
#
# Revision 1.32  2007/01/18 15:15:05  venelin
# Temporary hack to make the CERN reader working on a x86_64 bit
# architecture (using -m32 option for the compiler)
#
# Revision 1.31  2006/03/31 16:58:20  venelin
# Adding a check for SLC4 operating system and using the proper libs for it
#
# Revision 1.30  2006/02/01 17:21:53  venelin
# Adding a switch also to the subversion of compact
#
# Revision 1.29  2005/03/18 12:06:17  venelin
# Added "grep -i" to ignore case sensitivity
# reason - different config files on lxplus machines
#
# Revision 1.28  2005/03/01 13:36:14  venelin
# Changes in the paths when compiling in stand-alone mode
# ( PUBLIC_AFS=no )
#
# Revision 1.27  2004/12/06 02:33:48  venelin
# Modified in order to decide which libraries to take according to the platform we run on.
# Makes reader -  independent of the Linux version
#
# Revision 1.26  2004/10/27 10:49:29  clemenci
# fixed dependencies to make it possible to compile using -j option of gmake
#
# Revision 1.25  2004/10/18 18:56:21  clemenci
# Changed COmPACT version: 7.1->7.2
#
# Revision 1.24  2004/10/04 13:55:58  clemenci
# other fix for portability
#
# Revision 1.23  2004/10/04 12:27:04  clemenci
# same as before... missed one line
#
# Revision 1.22  2004/10/04 12:19:16  clemenci
# fixed non-Linux compilation issues (AIX)
#
# Revision 1.21  2004/09/21 13:35:18  clemenci
# Fixed the bug on LKR time offset in SCompact.
# The code to correct the wrong value is in the reader and has to be activated in the Makefile.
#
# Revision 1.20  2004/07/21 00:59:12  clemenci
# added hyperCOmPACT filtering routines
#
# Revision 1.19  2004/07/15 22:33:08  clemenci
# the path to cfortran.h was missing
#
# Revision 1.18  2004/07/11 03:12:53  clemenci
# cosmetics, improvements and bug fixes
#
# Revision 1.17  2004/06/23 15:41:48  clemenci
# added pritHelp.c, and a trick to avoid the generation of the dependencied when cleaning
#
# Revision 1.16  2004/05/28 12:37:32  maiera
# Updated the Makefiles in reader and rlib for compact-7.1
#
# Revision 1.15  2004/04/28 17:50:58  andrew
# Added calls to open and close the database
# added the sqlite lib to the link step
#
# Revision 1.14  2004/04/02 13:29:32  andrew
# Added the root include dir to depends
#
# Revision 1.13  2004/02/17 16:12:26  andrew
# Added hypercompact support
#
# Revision 1.12  2004/01/29 17:02:55  andrew
# Added -lshift to the debug version
#
# changed the search order of includes
#
# Revision 1.11  2003/11/25 13:13:51  andrew
# Fix Makefile to USE_RFIO=yes (as before)
#
# Revision 1.10  2003/10/24 11:09:26  andrew
# set up compact for use with ROOT. Set USE_ROOT to yes in the Makefile to use
# ROOT
#
# Revision 1.9  2003/09/02 08:00:11  andrew
# fixed a bug which prevented saving the burst scaler
#
# fixed the reader makefile to point to zlib-1.1.4
#
# Revision 1.8  2003/04/16 13:26:36  andrew
# Added for new cmp2scmp logic
#
# Revision 1.7  2003/04/16 13:19:02  andrew
# Changed to accomodate for compact 7
#
# Revision 1.6  2002/10/22 14:52:55  mmartini
# *** empty log message ***
#
# Revision 1.5  2002/07/17 00:19:50  andrew
# Removed a comment which prevented linkiing
#







