# makefile of the MFE and Annealed MAP simulations studies (using GIMP and Linux binaries for libDAI)
# run with make (or make rebuild)
# ELF 64 bit by default (run with "make MACHFLAG=-m32 DAILIB=-ldai32" for ELF 32 bit to downgrade)
# see LICENSE for license information

CC = g++
MACHFLAG ?= -m64
DAILIB ?= -ldai
CFLAGS = -O0 -DNDEBUG $(MACHFLAG) -ffast-math -Wall -g -fPIC -std=c++11 -I./include
AR = ar
ARFLAGS = -rv

OBJECT = ./object
LIBDIR = ./lib
SOURCE = ./source
RELEASE = ./release

REDIRC = >errors_gcc.txt 2>&1
REDIRL = >errors_ld.txt 2>&1

# look in source directory for source files
VPATH = $(OBJECT);$(SOURCE);$(INCLUDE)

.DEFAULT_GOAL := simulate

objs = mfesim_main.o mfe.o ann.o rel.o util.o map_indep.o

# make rebuild cleans and rebuilds all targets
rebuild: clean simulate bif2fg

.PHONY: clean
clean:
	rm -f $(OBJECT)/*.o

# builds a linux executable for the simulation software
simulate: $(addprefix $(OBJECT)/,$(objs))
	$(CC) -static $(CFLAGS) -o $(RELEASE)/mfesim $(addprefix $(OBJECT)/,$(objs)) -L $(LIBDIR) $(DAILIB) -lgmpxx -lgmp $(REDIRL)

# builds a helper executable for transferring .bif network descriptions into libDAI .fg factor graphs
bif2fg: $(OBJECT)/bif2fg.o
	$(CC) -static $(CFLAGS) -o $(RELEASE)/bif2fg $(OBJECT)/bif2fg.o -L $(LIBDIR) $(DAILIB) -lgmpxx -lgmp $(REDIRL)

# rules for individual objects
$(OBJECT)/mfesim_main.o : $(SOURCE)/mfesim_main.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/mfesim_main.cpp -o $(OBJECT)/mfesim_main.o $(REDIRC)

$(OBJECT)/mfe.o : $(SOURCE)/mfe.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/mfe.cpp -o $(OBJECT)/mfe.o $(REDIRC)

$(OBJECT)/ann.o : $(SOURCE)/ann.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/ann.cpp -o $(OBJECT)/ann.o $(REDIRC)

$(OBJECT)/rel.o : $(SOURCE)/rel.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/rel.cpp -o $(OBJECT)/rel.o $(REDIRC)

$(OBJECT)/map_indep.o : $(SOURCE)/map_indep.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/map_indep.cpp -o $(OBJECT)/map_indep.o $(REDIRC)

$(OBJECT)/util.o : $(SOURCE)/util.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/util.cpp -o $(OBJECT)/util.o $(REDIRC)

$(OBJECT)/bif2fg.o : $(SOURCE)/bif2fg.cpp
	$(CC) $(CFLAGS) -c $(SOURCE)/bif2fg.cpp -o $(OBJECT)/bif2fg.o $(REDIRC)

