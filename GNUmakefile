# $Id: GNUmakefile,v 1.3 2014/10/13 18:43:32 swesterd Exp $
# --------------------------------------------------------------
# If you want to use the graphics (visualisation) in the programm
# just comment the line : export G4VIS_NONE=1
# and recompile the program by 2 commands:
# make & make clean
# ---------------------------------------------------------------

G4BIN := ./

name := g4ds

G4TARGET := $(name)
G4EXLIB := true

G4VIS_USE := 0
BX_USE_OPERA := 0

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

HOST = $(shell hostname)
ARGUS = argus.Princeton.EDU
DEATHSTAR = deathstar
YODA = phy-yoda.Princeton.EDU
PALPATINE = palpatine.Princeton.EDU
CPPFLAGS += -I$(ROOTSYS)/include -DBX_USE_OPERA 

ifeq ($(HOST),$(ARGUS))
EXTRALIBS = $(shell root-config --glibs) /sw/lib/libxerces-c.27.dylib
else
ifeq ($(HOST),$(DEATHSTAR))
EXTRALIBS = -lxerces-c -L/usr/local/bin/root5.34/lib -lMathCore -lCore -lCint -lHist -lNet -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lThread -lGui -pthread -rdynamic
CPPFLAGS += -I/usr/local/root_v5.32.04/include/root/
else
ifeq ($(HOST),$(PALPATINE))
EXTRALIBS = $(shell root-config --glibs --cflags)
CPPFLAGS += -I$(ROOTSYS)/include/root
else
ifeq ($(HOST),$(YODA))
EXTRALIBS = $(shell root-config --glibs --cflags) -L/usr/local/root_v5.32.04/lib/root
CPPFLAGS += -I/usr/local/root_v5.32.04/include/root
else
EXTRALIBS = $(shell root-config --glibs) -lxerces-c
endif
endif
endif
endif

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps
	rm -f ./Darwin-g++/*.wrl
	rm -f .DAWN_*
	rm -f ./Linux-g++/*.wrl
logclean:
	rm -f ./Linux-g++/*.fil
	rm -f ./Linux-g++/*.log
	rm -f ./Linux-g++/*.out
#	rm -f ./Linux-g++/*.root
	rm -f ./Darwin-g++/*.fil
	rm -f ./Darwin-g++/*.log
	rm -f ./Darwin-g++/*.out
#	rm -f ./Darwin-g++/*.root
