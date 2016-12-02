set(CMAKE_C_FLAGS " -O2 ") #-mno-cygwin
set(CMAKE_CXX_FLAGS " -O2 ")
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,8388608")

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(SOAPnuke)

ADD_DEFINITIONS(-Wall -O2)#附加编译选项 

		 
add_library(system STATIC IMPORTED)
    set_property(TARGET system PROPERTY IMPORTED_LOCATION
         /usr/lib/libboost_system.dll.a)

add_library(thread STATIC IMPORTED)
    set_property(TARGET thread PROPERTY IMPORTED_LOCATION
         /usr/lib/libboost_thread.dll.a)

add_library(zlib STATIC IMPORTED)
    set_property(TARGET zlib PROPERTY IMPORTED_LOCATION
		/usr/lib/libz.a)

add_library(log4cplus STATIC IMPORTED)
    set_property(TARGET log4cplus PROPERTY IMPORTED_LOCATION
		/usr/local/lib/liblog4cplus.dll.a)

add_library(openssl STATIC IMPORTED)
    set_property(TARGET openssl PROPERTY IMPORTED_LOCATION
         /usr/lib/libcrypto.a) 

#add_library(cygwin STATIC IMPORTED)
#    set_property(TARGET cygwin PROPERTY IMPORTED_LOCATION
#        /usr/lib/libcygwin.a)	
		
#add_library(libc STATIC IMPORTED)
#    set_property(TARGET libc PROPERTY IMPORTED_LOCATION
#        /usr/lib/libc.a)

ADD_EXECUTABLE(SOAPnuke Main.cpp Common.cpp FilterProcessor.cpp FqBuffer.cpp
        	PeBuffer.cpp Logger.cpp SRNACommon.cpp SRNAProcessor.cpp
		SRNACleanFA.cpp
		SRNACleanFQ.cpp DGECommon.cpp DGEProcessor.cpp FqFile.cpp
		MetaCommon.cpp MetaProcessor.cpp)

TARGET_LINK_LIBRARIES(SOAPnuke system openssl thread log4cplus zlib pthread)# cygwin libc) 


