# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linuxx8664gcc

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 
DictSuf       = Dict.o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


ifeq ($(ARCH),linuxx8664gcc)
# Linux 
CXX           = g++ 
CXXFLAGS      = -pipe -m64 -Wall -W -fPIC
LD            = g++ 
LDFLAGS       = -m64 -O -Wl,--no-as-needed
SOFLAGS       = -shared -Wl,-soname,
endif


ifeq ($(ARCH),linux)
# Linux with gcc 2.7.2.x
CXX           = g++ 
CXXFLAGS      = -O -Wall -fPIC
LD            = g++ 
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxrh42)
# Linux with gcc 2.7.2.x (RedHat 4.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxdeb)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxdeb2)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxsuse6)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxegcs)
# Linux with egcs (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxkcc)
# Linux with the KAI compiler
CXX           = KCC
CXXFLAGS      = -fPIC +K0
LD            = KCC
LDFLAGS       = -O
SOFLAGS       =
endif

ifeq ($(ARCH),linuxppcegcs)
# MkLinux with egcs/glibc
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared -Wl,-soname,
endif

ifeq ($(ARCH),linuxia64gcc)
# Itanium Linux with gcc 2.9x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxia64sgi)
# Itanium Linux with sgiCC
CXX           = sgiCC
CXXFLAGS      = -O -Wall -fPIC
LD            = gsgiCC
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxalphaegcs)
# Alpha Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxarm)
# ARM Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

CXXFLAGS     += $(ROOTCFLAGS) 
LIBS          = $(ROOTLIBS) $(SYSLIBS) -lTreePlayer
GLIBS         = $(ROOTGLIBS) 
#GLIBS         = $(ROOTGLIBS) -lrfstream -L/usr/local/lib/ -lshift $(SYSLIBS) 
#GLIBS         = $(ROOTGLIBS) -lRFIO -L/usr/local/lib/ -lshift $(SYSLIBS) 

#------------------------------------------------------------------------------

ana_clas12O          = ana_clas12.$(ObjSuf) \
                    ana_clas12Dict.$(ObjSuf)

ana_clas12S          = ana_clas12.$(SrcSuf) \
                      ana_clas12Dict.$(SrcSuf)

ana_clas12SO         = libana_clas12.$(DllSuf)

OBJS          = $(ana_clas12O)

SLIB          = $(ana_clas12SO) \
              $(ROOTSYS)/lib/libPhysics.so $(ROOTSYS)/lib/libTreePlayer.so

MAINO		= ana.$(ObjSuf)
MAINS		= ana.$(SrcSuf)
MAIN		= ana$(ExeSuf)
PROGRAMS	= $(MAIN)
#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(SLIB) $(PROGRAMS)

lib:            $(SLIB) 

doc:
	root.exe -b -q make_html.C

%.$(ObjSuf):%.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<

%Dict.$(SrcSuf):%.h
	rootcint -f $@ -c $<

lib%.$(DllSuf): %Dict.$(ObjSuf) %.$(ObjSuf)
	$(LD) $(SOFLAGS) $(CXXFLAGS) $^ -o $@

$(MAIN): $(MAINO) $(SLIB)
	$(LD) $(LDFLAGS) $(GLIBS) $^ -o $@

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(SLIB) *Dict.* *.def
