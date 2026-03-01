S=Splus

CFLAGS=
FFLAGS= -sun4 -c
FORTRAN= f77

OBJS = drqssbr.o drqssbc.o dblas1.o
FUNCS = qss qss.fit
DATAS = data.s

FUNC = $(FUNCS:%=.Data/%)

.SUFFIXES:
.SUFFIXES: .o .c .c~ .r .r~ .f .f~ 

.f.o:	
	$(FORTRAN) $(FFLAGS) $*.f 

.r.o:
	$(FORTRAN) $(FFLAGS) $*.r 

DYNAMIC: local.o $(FUNC) DATA
	mv local.o /tmp/local.o
	make -f Makefile clean.dynamic
	mv /tmp/local.o .
	touch DYNAMIC

STATIC: local.o $(FUNC) DATA
	$S LOAD local.o
	make -f Makefile clean.static
	touch STATIC

local.o: $(OBJS)
	ld -r -d -o local.o $(OBJS) 

.Data/%: %.s
	$S <  $<

DATA:  data.s
	$S < data.s
	touch DATA

clean.static:
	rm -rf *.o DATA DYNAMIC

clean.dynamic:
	rm -rf *.o DATA STATIC
