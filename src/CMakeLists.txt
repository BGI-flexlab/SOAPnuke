#     SOAPnuke - A MapReduce Acceleration supported Software for integrated Quality Control and Preprocessing of High-Throughput Sequencing Data
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(SOAPnuke)

set(CMAKE_C_FLAGS " -O2 ")
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
# set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,8388608")
set(CMAKE_CXX_STANDARD 11)

ADD_DEFINITIONS(-Wall -O2)  #附加编译选项

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(RT REQUIRED)
include_directories(${RT_INCLUDE_DIR})
message(STATUS "${RT_INCLUDE_DIR}")
    
find_library(temp1 libboost_system.a PATHS ENV LD_LIBRARY_PATH)
add_library(system STATIC IMPORTED)
        set_property(TARGET system PROPERTY IMPORTED_LOCATION
                ${temp1})
find_library(temp2 libboost_thread.a PATHS ENV LD_LIBRARY_PATH)
add_library(thread STATIC IMPORTED)
        set_property(TARGET thread PROPERTY IMPORTED_LOCATION
                ${temp2})
find_library(temp3 libcrypto.a PATHS ENV LD_LIBRARY_PATH)
add_library(openssl STATIC IMPORTED)
        set_property(TARGET openssl PROPERTY IMPORTED_LOCATION
                ${temp3})
find_library(temp4 libz.a PATHS ENV LD_LIBRARY_PATH)
add_library(zlib STATIC IMPORTED)
        set_property(TARGET zlib PROPERTY IMPORTED_LOCATION
                ${temp4})
find_library(temp5 liblog4cplus.a PATHS ENV LD_LIBRARY_PATH)
add_library(log4cplus STATIC IMPORTED)
        set_property(TARGET log4cplus PROPERTY IMPORTED_LOCATION
                ${temp5})

ADD_EXECUTABLE(SOAPnuke Main.cpp Common.cpp FilterProcessor.cpp FqBuffer.cpp
        PeBuffer.cpp Logger.cpp SRNACommon.cpp SRNAProcessor.cpp
        SRNACleanFA.cpp
        SRNACleanFQ.cpp DGECommon.cpp DGEProcessor.cpp FqFile.cpp
        MetaCommon.cpp MetaProcessor.cpp)

TARGET_LINK_LIBRARIES(SOAPnuke openssl thread system log4cplus zlib pthread "${RT_LIBRARY}")
