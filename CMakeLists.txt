cmake_minimum_required(VERSION 3.7)
project(RaptRanker CXX)

set(CMAKE_CXX_STANDARD 14)

IF(NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

#set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 ")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
