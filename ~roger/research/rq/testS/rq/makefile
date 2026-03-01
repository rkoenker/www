# makefile for local CHAPTER
SHELL=/bin/sh

SRC= rqbr.r rqfn.r xlapack.f
OBJ= rqbr.o rqfn.o xlapack.o

# Use LOCAL_CFLAGS to add arguments for the C compiler
LOCAL_CFLAGS=
# Use LOCAL_CXXFLAGS to add arguments for the C++ compiler
LOCAL_CXXFLAGS=
# Use LOCAL_FFLAGS to add arguments for the FORTRAN compiler
LOCAL_FFLAGS=

# Use LOCAL_LIBS to add arguments or additional libraries to the linker
# This next line is needed to avoid an unresolved __d_sign for the rq package
LOCAL_LIBS="-lF77"
#LOCAL_LIBS=


include $(SHOME)/library/S_FLAGS

S.so:  $(OBJ) 
	$(SHOME)/cmd/Splus5 LIBRARY S.so  $(OBJ) $(LOCAL_LIBS)

dump:
	@if test -d ./.Data; then S dumpChapter $(SRC); fi

boot:
	@if test -s all.Sdata; \
	  then (BOOTING_S="TRUE" export BOOTING_S;  echo "terminate(should have been booting S)"| $(SHOME)/cmd/S); \
	fi


clean:
	-rm $(OBJ)

