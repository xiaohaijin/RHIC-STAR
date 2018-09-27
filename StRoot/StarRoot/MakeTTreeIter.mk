# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linuxegcs

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 


EVENTLIB      = $(EVENTSO)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


ifeq ($(ARCH),hpux)
# HP-UX with CC
CXX           = CC
CXXFLAGS      = -O +Z
LD            = CC
LDFLAGS       = -O +a1 -z
SOFLAGS       = -b
DllSuf        = sl
endif

ifeq ($(ARCH),hpuxacc)
# HP-UX 10.x with aCC
CXX           = aCC
CXXFLAGS      = -O +Z
LD            = aCC
LDFLAGS       = -O -z
SOFLAGS       = -b
endif

ifeq ($(ARCH),hpuxegcs)
# HP-UX 10.x with g++
CXXFLAGS      = -O -fPIC
CXX           = g++
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -fPIC -shared
endif

ifeq ($(ARCH),aix)
# IBM AIX
CXX           = xlC
CXXFLAGS      = -O
LD            = xlC
LDFLAGS       = -O
SOFLAGS       =
ROOTLIBS     := $(shell root-config --nonew --libs)
ROOTGLIBS    := $(shell root-config --nonew --glibs)
endif

ifeq ($(ARCH),aixegcs)
# IBM AIX with GCC
CXX           = g++
CXXFLAGS      = -O
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),solaris)
# Solaris CC
CXX           = /opt/SUNWspro/bin/CC
CXXFLAGS      = -O -KPIC
LD            = /opt/SUNWspro/bin/CC
LDFLAGS       = -O
SOFLAGS       = -G
endif

ifeq ($(ARCH),solarisCC5)
# Solaris CC 5.0
CXX           = CC
CXXFLAGS      = -O -KPIC
LD            = CC
LDFLAGS       = -O
SOFLAGS       = -G
endif

ifeq ($(ARCH),solarisegcs)
# Solaris egcs
CXX           = g++
CXXFLAGS      = -O -fPIC
LD            = CC
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),solarisgcc)
# Solaris gcc
CXX           = g++
CXXFLAGS      = -O -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),solariskcc)
# Solaris kcc
CXX           = KCC
CXXFLAGS      = -O4 -KPIC
LD            = KCC
LDFLAGS       = -O4
SOFLAGS       =
endif

ifeq ($(ARCH),solarisx86)
# Solaris CC on Intel
CXX           = CC
CXXFLAGS      = -O -KPIC
LD            = CC
LDFLAGS       = -O
SOFLAGS       = -G
endif

ifeq ($(ARCH),sgicc)
# SGI
CXX           = CC -n32  -I/usr/include/CC.sgi
CXXFLAGS      = -O
LD            = CC -n32  -I/usr/include/CC.sgi
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),sgiegcs)
# SGI 6.x with EGCS
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O -Wl,-u,__builtin_new -Wl,-u,__builtin_delete -Wl,-u,__nw__FUiPv
SOFLAGS       = -shared
endif

ifeq ($(ARCH),sgin32egcs)
# SGI 6.x with EGCS for n32 ABI
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O -L/usr/lib32 -Wl,-woff,134
SOFLAGS       = -shared
endif

ifeq ($(ARCH),sgigcc)
# SGI with GCC
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O -Wl,-u,__builtin_new -Wl,-u,__builtin_delete -Wl,-u,__nw__FUiPv
SOFLAGS       = -shared
endif

ifeq ($(ARCH),sgikcc)
# SGI with KCC
CXX           = KCC -n32 --no_exceptions
CXXFLAGS      = -O
LD            = KCC -n32 --no_exceptions
LDFLAGS       = -O
SOFLAGS       =
endif

ifeq ($(ARCH),alphagcc)
# Alpha/OSF with g++
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -Wl,-expect_unresolved,* -shared
endif

ifeq ($(ARCH),alphaegcs)
# Alpha/OSF with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -Wl,-expect_unresolved,* -shared
endif

ifeq ($(ARCH),alphakcc)
# Alpha/OSF with kai compiler (not yet valid)
CXX           = g++
CXXFLAGS      = -O -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -Wl,-expect_unresolved,* -shared
endif

ifeq ($(ARCH),alphacxx6)
# Alpha/OSF with cxx6
CXX           = cxx
CXXFLAGS      = -O0
LD            = cxx
LDFLAGS       = -O
SOFLAGS       = -Wl,-expect_unresolved,* -shared
endif

ifeq ($(ARCH),alphacxx)
# Alpha/OSF with cxx5
CXX           = cxx
CXXFLAGS      = -O
LD            = cxx
LDFLAGS       = -O
SOFLAGS       = -Wl,-expect_unresolved,* -call_shared
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

ifeq ($(ARCH),mklinux)
# MkLinux with libc5
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared -Wl,-soname,
endif

ifeq ($(ARCH),freebsd)
# FreeBSD with libc5
CXX           = g++
CXXFLAGS      = -O -pipe -W -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared -Wl,-x
endif

ifeq ($(ARCH),freebsd4)
# FreeBSD with glibc
CXX           = g++
CXXFLAGS      = -O -pipe -W -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared -Wl,-x
endif

ifeq ($(ARCH),hiux)
# Hitachi HIUX
CXX           = g++
CXXFLAGS      = -O2 -fPIC
LD            = g++
LDFLAGS       = -Wl,+s
SOFLAGS       = -Wl,-b,-E -nostdlib -nostartfiles
DllSuf        = sl
endif

ifeq ($(ARCH),win32)
# Windows with the VC++ compiler
ObjSuf        = obj
SrcSuf        = cxx
ExeSuf        = .exe
DllSuf        = dll
OutPutOpt     = -out:
CXX           = cl
CXXOPT        = -O2
#CXXOPT        = -Z7
CXXFLAGS      = $(CXXOPT) -G5 -GR -MD -DWIN32 -D_WINDOWS -nologo \
                -DVISUAL_CPLUSPLUS -D_X86_=1 -D_DLL
LD            = link
LDOPT         = -opt:ref
#LDOPT         = -debug
LDFLAGS       = $(LDOPT) -pdb:none -nologo -nodefaultlib -incremental:no
SOFLAGS       = -DLL
SYSLIBS       = msvcrt.lib oldnames.lib kernel32.lib  ws2_32.lib mswsock.lib \
                advapi32.lib  user32.lib gdi32.lib comdlg32.lib winspool.lib \
                msvcirt.lib
EVENTLIB      = libTTreeIter.lib

ROOTLIBS     := $(shell root-config --nonew --libs)
ROOTGLIBS    := $(shell root-config --nonew --glibs)
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------

EVENTO        = TTreeIter.$(ObjSuf) \
                TTreeIterDict.$(ObjSuf)

EVENTS        = TTreeIter.$(SrcSuf) \
                TTreeIterDict.$(SrcSuf)


EVENT         = TTreeIter$(ExeSuf)
EVENTSO       = libTTreeIter.$(DllSuf)





OBJS          = $(EVENTO)


#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)


$(EVENTSO):     $(EVENTO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $(EVENTSO) \
		   $(LIBS) -p 0 $(EVENTO)
else
ifeq ($(ARCH),alphacxx)
# due to a bug in cxx/ld under osf3.xx, one cannot use cxx to generate
# a shared library. One must use ld instead.
		ld -L/usr/lib/cmplrs/cxx -rpath /usr/lib/cmplrs/cxx \
		   -expect_unresolved "*" -g0 -O1 -shared \
		   /usr/lib/cmplrs/cc/crt0.o /usr/lib/cmplrs/cxx/_main.o \
		   -o libTTreeIter.so TTreeIter.o TTreeIterDict.o -lcxxstd -lcxx \
		   -lexc -lots -lc
else
ifeq ($(ARCH),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $(EVENTO) $(OutPutOpt) $(EVENTSO)
endif
endif
endif
		@echo "$@ done"

$(EVENT):       $(EVENTSO) 
		$(LD) $(LDFLAGS)  $(EVENTLIB) $(LIBS) \
		   $(OutPutOpt)$(EVENT)
		@echo "$@ done"





clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations

.SUFFIXES: .$(SrcSuf)

###

TTreeIter.$(ObjSuf): TTreeIter.h

TTreeIterDict.$(SrcSuf): TTreeIter.h TTreeIterLinkDef.h
	@echo "Generating dictionary TTreeIterDict..."
	@rootcint -f TTreeIterDict.$(SrcSuf) -c TTreeIter.h TTreeIterLinkDef.h


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
