cmake_minimum_required(VERSION 3.11)

project(imcqd)

#set(CMAKE_C_COMPILER "gcc-11")
#set(CMAKE_CXX_COMPILER "g++-11")
#set(CMAKE_MACOSX_RPATH OFF)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED true)
#set(CMAKE_CXX_FLAGS "-Wall -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
#set(CMAKE_EXE_LINKER_FLAGS "-static")
set(CMAKE_CXX_FLAGS_DEBUG "-DNDEBUG")

set(Boost_USE_STATIC_LIBS ON)
set(BOOST_ROOT_DIR /usr/local/opt/boost)
set(OPENSSL_ROOT_DIR /usr/local/opt/openssl)
set(ZLIB_ROOT_DIR /usr/local/opt/zlib)

find_package(Boost REQUIRED COMPONENTS date_time system regex iostreams)
find_package(GSL REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(ZLIB)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${OPENSSL_INCLUDE_DIR})

include_directories(".")
include_directories(helper)
include_directories(graph)

add_subdirectory(helper)
add_subdirectory(graph)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
    helper
    graph
    ${Boost_DATE_TIME_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_IOSTREAMS_LIBRARY}
    ${ZLIB_LIBRARIES}
)


