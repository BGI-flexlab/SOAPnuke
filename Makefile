cc := g++
src := ./src
obj := ./obj

MIN_GCC_VERSION = 4.7
GCC_VERSION := $(shell gcc -dumpversion)
IS_GCC_ABOVE_MIN_VERSION := $(shell echo `gcc -dumpversion | cut -f1-2 -d.` \>= $(MIN_GCC_VERSION) | bc )
ifeq "$(IS_GCC_ABOVE_MIN_VERSION)" "1"
    GCC_VERSION_STRING := "GCC version Passes, $(GCC_VERSION) >= $(MIN_GCC_VERSION)"
else
    GCC_VERSION_STRING := "Warning: GCC version $(GCC_VERSION) is lower than $(MIN_GCC_VERSION)."
endif
$(info $(GCC_VERSION_STRING))



all=SOAPnuke
exe=SOAPnuke
MIN_ZLIB_VERSION = 1.2.3.5
#ZLIB_VERSION := $(shell readlink $(whereis libz) | sed 's/libz.so.\(.*\)$/\1/g')
ZLIBLOCATION:=$(shell whereis libz)
ZLIB_VERSION=$(shell readlink $(ZLIBLOCATION) | awk -F 'so.' '{print $$NF}')
#$(info $(ZLIB_VERSION_STRING))
IS_ZLIB_ABOVE_MIN_VERSION := $(shell expr "$(ZLIB_VERSION)" ">=" "$(MIN_ZLIB_VERSION)")
ifeq "$(IS_ZLIB_ABOVE_MIN_VERSION)" "1"
	ZLIB_VERSION_STRING := "ZLIB version Passes, $(ZLIB_VERSION) >= $(MIN_ZLIB_VERSION)"
else
	ZLIB_VERSION_STRING := "Warning: ZLIB version $(ZLIB_VERSION) is lower than $(MIN_ZLIB_VERSION)."
endif
$(info $(ZLIB_VERSION_STRING))

#set USEHTS true if you want use filterHts module
USEHTS=false

source := $(wildcard ${src}/*.cpp)
source := $(filter-out ./src/mGzip.cpp, $(source))
LIBS := -lz -lpthread -labc
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)
$(warning $(LD_FLAGS))
DFLAG=-lz -lpthread
CXXFLAGS=-std=c++11 -g -O3
ifeq "$(USEHTS)" "true"
	CXXFLAGS+=-D _PROCESSHTS
	DFLAG+= -lhts
else
	source := $(filter-out ./src/processHts.cpp, $(source))
endif
#$(warning $(source))
objfile := $(patsubst %.cpp,${obj}/%.o,$(notdir ${source}))
objfile:=$(filter-out processHts.o,$(objfile))
$(exe):${objfile}
	$(cc) $(objfile) -o $@ $(DFLAG)
${obj}/%.o:${src}/%.cpp mk_dir
	$(cc) $(CXXFLAGS) -c $< -o $@
mk_dir:
	@if test ! -d $(obj);\
	then\
		mkdir $(obj);\
	fi
.PHONY:clean
clean:
	rm -f obj/*.o
	rm -f SOAPnuke
