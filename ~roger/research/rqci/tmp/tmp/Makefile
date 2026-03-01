# Makefile created by S-PLUS utility CHAPTER, version 3.24
# S-PLUS Version 3.2 Release 1 for Sun SPARC, SunOS 4.x : 1993 

SHOME=/home/splus3.2

#WHICH_LOAD=static.load
WHICH_LOAD=dyn.load

EXTRA_OBJ_FILES=
EXTRA_C_NAMES=
EXTRA_F_NAMES=

# This is a Makefile produced by the S-PLUS utility CHAPTER. It guides
# compilation of C, Ratfor and Fortran source code, loading of the resulting
# object code, and installation of functions and helpfiles.
# 
# Overview: 
# ========
# 
# You will need to attend to a few 'make' macro settings at the top of this
# file first.  See the Detailed Instructions below for these. Then, if you
# want to compile and load the default set of C and Fortran routines, type
# (1) 'make install' to install functions on .Data and thereby enable
# building of the default list of compiled routine names, then (2) 'make
# load' to compile the object code, build the set of names to load, and load
# this object code.  Other, secondary, 'make' targets are 'install.funs' and
# 'install.help' ('install' makes both of those), 'clean' (remove only
# intermediate object code) and 'virgin' (remove object code, the local
# standalone S-PLUS executable if you built it, and the directory
# .Data).
# 
# Detailed Instructions:
# =====================
# 
# 1. On the line near the top that starts with "SHOME=", change the
# directory pathname on the right of the equals sign to the pathname of the
# top S-PLUS directory on your machine; use Splus SHOME to get this
# pathname.
# 
# 2. Near the top, the macro WHICH_LOAD specifies whether 'make' targets
# 'load' and 'all' should make a standalone S-PLUS binary executable named
# "local.Sqpe" which includes the local C and Fortran code, or a
# dyn.loadable file "rq_l.o" containing only this local code.  To
# make a complete standalone executable, WHICH_LOAD should be set to
# "static.load"; the only other possible value is "dyn.load".  Make sure a
# comment symbol (#) appears in front of the setting you DON'T want.
# The default has been set to "static.load".
# 
# Incidentally, you can always make targets 'static.load' or 'dyn.load'
# regardless of the value of WHICH_LOAD. The macro just specifies what the
# generic targets 'load' and 'all' should make.
# 
# 3. Near the top, the macro EXTRA_OBJ_FILES names object files (ending in
# ".o") which should be loaded in addition to those which are part of this
# chapter.  Add any additional object file names you want to load, separated
# by spaces, on the right hand side of the equals sign.  There MUST be
# corresponding C (".c"), Ratfor (".r") or Fortran (".f") files present in
# this directory.
# 
# If you add file names to EXTRA_OBJ_FILES, you must also add the names of
# the corresponding C functions or Fortran subroutines to the macros
# EXTRA_C_NAMES and EXTRA_F_NAMES respectively, underneath EXTRA_OBJ_FILES.
# The names added should be as they appear at the source-code level; do not
# add leading or trailing underscores unless these are present in the
# source-level name.
# 
# You won't generally need to change any more 'make' macros, or indeed
# anything else in this Makefile.
# 
# 4. The file rq_i.c is a C language source file which will ensure
# that the make target 'load' will load C and Fortran routines referenced
# only through .C or .Fortran calls in S functions. If rq_i.c does
# not exist, it will be made during 'make load' based on the objects in
# .Data. If rq_i.c exists already and you have since added new C
# or Fortran routine names to EXTRA_C_NAMES and EXTRA_F_NAMES, you must
# remove rq_i.c so that it will get remade using this new
# information. Do this now.
# 
# 5. Type 'make install.funs' to create the S functions on .Data.  This
# must precede 'make load' in Step 6. (You can install the helpfiles at
# this point as well by typing making 'install' instead of 'install.funs'.)
# 
# 6. Type 'make load'. This will build the list of routines to load, then
# compile the C, Ratfor and Fortran source files specified by this Makefile,
# and lastly load the routines named on the list into either a standalone
# S-PLUS executable or a dyn.loadable file.
# 
# 7. To remove unnecessary object files, type 'make clean'. This will not
# remove the dyn.loadable file rq_l.o or the standalone local.Sqpe.
# To clean out everything including .Data and start over, use 'make
# virgin'.
# 
# 8. See the helpfile for the "library" function in S-PLUS for hints on
# how to install the new code as a library section.


# ========================= End instructions. =============================

default : all
include $(SHOME)/newfun/lib/S_makefile
chapters=../rq
SRC= rq.r
OBJ= rq.o $(EXTRA_OBJ_FILES)
CFLAGS=-O2
FFLAGS=-O2
FUNS=force_ld.q rq.s
HELPS= ranks.d rrs.test.d rq.d trq.d qrq.d trq.print.d
FLAGS=
RM=-rm

all load: $(WHICH_LOAD)

static.load: rq.a
	Splus LOAD $(FLAGS) CHAPTERS='"$(chapters)"'

dyn.load rq_l.o: rq.a
	ld -r -o rq_l.o rq_i.o rq.a $(FLAGS)
	@echo dynamically loadable file in rq_l.o

rq.a: rq_i.o $(OBJ) 
	Splus LIBRARY rq.a rq_i.o $(OBJ)

rq_i.c:
	-mkdir .Data
	Splus make.init rq .Data

install : install.funs install.help

funs install.funs : $(FUNS) .Data
	Splus QINSTALL .Data $(FUNS)

force_ld.q :
	( echo "force.loading <- function(){" ;\
	set $(EXTRA_C_NAMES) terminator ;\
	while test $$1 != terminator ;\
	do \
		echo $$1 | sed 's:.*:.C("&"):' ;\
		shift ;\
	done ;\
	set $(EXTRA_F_NAMES) terminator ;\
	while test $$1 != terminator ;\
	do \
		echo $$1 | sed 's:.*:.Fortran("&"):' ;\
		shift ;\
	done ;\
	echo } \
	) > $@


# force_ld.q :
# 	echo "force.loading <- function(){ X\
# 		stop(\"should not be executed\") X\
# 		# Add .C and .Fortran calls here to force routines to  X\
# 		# be loaded. For example, .C(\"fun1\"), .Fortran(\"fun2\"). X\
# 	}" | tr X \\012  > $@

help install.help : $(HELPS) .Data/.Help
	Splus HINSTALL .Data/.Help $(HELPS)
	Splus help.findsum .Data

.Data :
	-mkdir $@
.Data/.Help : .Data
	-mkdir $@

virgin : clean virgin.std
# add additional cleanup rules/targets above here for target 'virgin'
clean :
	$(RM) -f $(OBJ) rq_i.o S_load_time.[oc] core
virgin.std :
	$(RM) -rf .Data
	$(RM) -f rq.a local.Sqpe rq_l.o rq_i.c force_ld.q
