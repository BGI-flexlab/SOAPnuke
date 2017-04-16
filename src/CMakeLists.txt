set(CMAKE_C_FLAGS " -O2 ")
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,8388608")
set(CMAKE_CXX_STANDARD 98)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(SOAPnuke)

ADD_DEFINITIONS(-Wall -O2)#附加编译选项

include_directories(/usr/local/include)

add_library(system STATIC IMPORTED)
    set_property(TARGET system PROPERTY IMPORTED_LOCATION
		/usr/local/lib/libboost_system.a)

add_library(thread STATIC IMPORTED)
    set_property(TARGET thread PROPERTY IMPORTED_LOCATION
		/usr/local/lib/libboost_thread.a)

add_library(zlib STATIC IMPORTED)
    set_property(TARGET zlib PROPERTY IMPORTED_LOCATION
		/usr/local/lib/libz.a)

#log4cplus version 1.0.3
#wget https://sourceforge.net/projects/log4cplus/files/log4cplus-stable/1.0.3/log4cplus-1.0.3.tar.bz2
add_library(log4cplus STATIC IMPORTED)
    set_property(TARGET log4cplus PROPERTY IMPORTED_LOCATION
		/usr/local/lib/liblog4cplus.a)

add_library(openssl STATIC IMPORTED)
    set_property(TARGET openssl PROPERTY IMPORTED_LOCATION
		/usr/local/lib/libcrypto.a)


ADD_EXECUTABLE(SOAPnuke Main.cpp Common.cpp FilterProcessor.cpp FqBuffer.cpp
		PeBuffer.cpp Logger.cpp SRNACommon.cpp SRNAProcessor.cpp
		SRNACleanFA.cpp
		SRNACleanFQ.cpp DGECommon.cpp DGEProcessor.cpp FqFile.cpp
		MetaCommon.cpp MetaProcessor.cpp)

TARGET_LINK_LIBRARIES(SOAPnuke system openssl thread log4cplus zlib pthread)
