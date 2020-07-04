cc := g++
src := ./src
obj := ./obj

INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

source := $(wildcard ${src}/*.cpp)
objfile := $(patsubst %.cpp,${obj}/%.o,$(notdir ${source}))
all=SOAPnuke
exe=SOAPnuke

CXXFLAGS := -std=c++11 -g -O3 $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) $(CXXFLAGS)
LIBS := -lz -lpthread -lhts
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)

MIN_GCC_VERSION = "4.7"
GCC_VERSION := $(shell gcc -dumpversion)
IS_GCC_ABOVE_MIN_VERSION := $(shell expr "$(GCC_VERSION)" ">=" "$(MIN_GCC_VERSION)")
ifeq "$(IS_GCC_ABOVE_MIN_VERSION)" "1"
    GCC_VERSION_STRING := "GCC version Passes, $(GCC_VERSION) >= $(MIN_GCC_VERSION)"
else
    GCC_VERSION_STRING := "Warning: GCC version $(GCC_VERSION) is lower than $(MIN_GCC_VERSION)."
endif
$(info $(GCC_VERSION_STRING))

MIN_ZLIB_VERSION = "1.2.3.5"
ZLIB_LOCATION := $(shell ldconfig -p | grep libz | head -n 1 | cut -d ">" -f 2) 
ZLIB_VERSION := $(shell readlink $(ZLIB_LOCATION) | sed 's/\(.*\)libz\.so\.\(.*\)$/\2/g/'| sed 's/.//')
IS_ZLIB_ABOVE_MIN_VERSION := $(shell expr "$(ZLIB_VERSION)" ">=" "$(MIN_ZLIB_VERSION)")
ifeq "$(IS_ZLIB_ABOVE_MIN_VERSION)" "1"
    ZLIB_VERSION_STRING := "ZLIB version Passes, $(ZLIB_VERSION) >= $(MIN_ZLIB_VERSION)"
else
    ZLIB_VERSION_STRING := "Warning: ZLIB version $(ZLIB_VERSION) is lower than $(MIN_ZLIB_VERSION)."
endif
$(info $(ZLIB_VERSION_STRING))


$(exe):${objfile}
	$(cc) $(objfile) -o $@ $(LD_FLAGS)

${obj}/%.o:${src}/%.cpp mk_dir
	$(cc) -c $< -o $@ $(CXXFLAGS)

mk_dir:
	@if test ! -d $(obj) ; \
	then \
		mkdir $(obj) ; \
	fi
