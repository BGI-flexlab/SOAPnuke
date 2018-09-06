cc := g++
src := ./src
obj := ./obj

source := $(wildcard ${src}/*.cpp)
objfile := $(patsubst %.cpp,${obj}/%.o,$(notdir ${source}))
all=SOAPnuke
exe=SOAPnuke


$(exe):${objfile}
	$(cc) $(objfile) -o $@ -L/hwfssz1/BIGDATA_COMPUTING/gongchun/zlib/lib -lz -lpthread

${obj}/%.o:${src}/%.cpp mk_dir
	$(cc) -std=c++11 -g -O3 -c $< -o $@ -I/hwfssz1/BIGDATA_COMPUTING/gongchun/zlib/include 
mk_dir:
	@if test ! -d $(obj) ; \
	then \
		mkdir $(obj) ; \
	fi
