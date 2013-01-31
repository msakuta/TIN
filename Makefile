# Time-stamp: </wrf/c/tin/Dist/Makefile, Thu,  1 Nov 2001, 17:46:42 EST, wrf@benvolio.ecse.rpi.edu>

CCC=g++
CPPFLAGS=-O3   -W -felide-constructors -fstrength-reduce -finline-functions 
LDFLAGS=-lm

DISTFILES = \
	README tin.cc tin.help doadir doadirz Makefile adirondacks.elevs.gz

all: tin


clean:
	rm centers.all* edges.all* triangles.all* stats*


# Make a distribution package
# TODO:  add sample files, README, fns

dist:  $(DISTFILES) 
	tar zcf tin.tgz $(DISTFILES)
