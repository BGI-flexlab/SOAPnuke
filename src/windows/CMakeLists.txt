CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(SOAPnuke)

INCLUDE_DIRECTORIES(C:/Users/guxiang/Desktop/Windows_lib/zlib128-dll/include
    C:/Users/guxiang/Desktop/Windows_lib/log4cplus-1.2.0-rc3/include
	F:/OpenSSL-Win32/include
	C:/Users/guxiang/Desktop/Windows_lib/getopt9/include
	#F:/cygwin64/usr/include
	C:/Users/guxiang/Desktop/Windows_lib/boost_1_59_0)

    
# add_library(thread STATIC IMPORTED)
    # set_property(TARGET thread PROPERTY IMPORTED_LOCATION
        # C:/Users/guxiang/Desktop/SoapNuke-Windows/boost_1_59_0/stage/lib/libboost_thread-vc120-mt-gd-1_59.lib)

# add_library(zlib STATIC IMPORTED)
    # set_property(TARGET zlib PROPERTY IMPORTED_LOCATION
		# C:/Users/guxiang/Desktop/SoapNuke-Windows/zlib128-dll/lib/zdll.lib)

# add_library(log4cplus STATIC IMPORTED)
    # set_property(TARGET log4cplus PROPERTY IMPORTED_LOCATION
		# C:/Users/guxiang/Desktop/SoapNuke-Windows/log4cplus-1.2.0-rc3/msvc10/Win32/bin.Debug/log4cplusD.lib)

# add_library(openssl STATIC IMPORTED)
    # set_property(TARGET openssl PROPERTY IMPORTED_LOCATION
        # F:/OpenSSL/lib/VC/static/libeay32MTd.lib) 
		
# add_library(pthread STATIC IMPORTED)
    # set_property(TARGET openssl PROPERTY IMPORTED_LOCATION
        # C:/Users/guxiang/Desktop/SoapNuke-Windows/pthread/Pre-built.1/lib/pthreadVC1.lib)

# add_library(getopt STATIC IMPORTED)
    # set_property(TARGET getopt PROPERTY IMPORTED_LOCATION
        # C:/Users/guxiang/Desktop/SoapNuke-Windows/getopt9/lib/getopt.lib)

ADD_EXECUTABLE(SOAPnuke Main.cpp Common.cpp FilterProcessor.cpp FqBuffer.cpp
        PeBuffer.cpp Logger.cpp SRNACommon.cpp SRNAProcessor.cpp SRNACleanFA.cpp
		SRNACleanFQ.cpp DGECommon.cpp DGEProcessor.cpp FqFile.cpp
		MetaCommon.cpp MetaProcessor.cpp)

#TARGET_LINK_LIBRARIES(SOAPnuke getopt openssl thread log4cplus zlib pthread)


