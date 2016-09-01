CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(SOAPnuke)

ADD_DEFINITIONS(-Wall -O2)#附加编译选项 

    INCLUDE_DIRECTORIES(/ifs4/PROTON_TEST/software/zlib/include/
            /ifs4/PROTON_TEST/software/log4cplus/include/
			/ifs4/PROTON_TEST/software/openssl/include/openssl
			/ifs4/ISDC_BD/zhangyong2/lib/boost_1_49/include/ /usr/include)

#LINK_DIRECTORIES(/usr/local/lib/)
    
add_library(thread STATIC IMPORTED)
    set_property(TARGET thread PROPERTY IMPORTED_LOCATION
            /ifs4/ISDC_BD/zhangyong2/lib/boost_1_49/lib/libboost_thread.a)

add_library(zlib STATIC IMPORTED)
    set_property(TARGET zlib PROPERTY IMPORTED_LOCATION
		/ifs4/PROTON_TEST/software/zlib/lib/libz.a)

add_library(log4cplus STATIC IMPORTED)
    set_property(TARGET log4cplus PROPERTY IMPORTED_LOCATION
            /ifs4/PROTON_TEST/software/log4cplus/lib/liblog4cplus.a)

add_library(openssl STATIC IMPORTED)
    set_property(TARGET openssl PROPERTY IMPORTED_LOCATION
           /ifs4/PROTON_TEST/software/openssl/lib/libcrypto.a) 


ADD_EXECUTABLE(SOAPnuke Main.cpp Common.cpp FilterProcessor.cpp FqBuffer.cpp
        	PeBuffer.cpp Logger.cpp SRNACommon.cpp SRNAProcessor.cpp
		SRNACleanFA.cpp
		SRNACleanFQ.cpp DGECommon.cpp DGEProcessor.cpp FqFile.cpp
		MetaCommon.cpp MetaProcessor.cpp)

TARGET_LINK_LIBRARIES(SOAPnuke openssl thread log4cplus zlib pthread) 


