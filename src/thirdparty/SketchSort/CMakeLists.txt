cmake_minimum_required(VERSION 3.7)

find_package(Boost)
IF (Boost_FOUND)
    include_directories(SketchSort ${Boost_INCLUDE_DIR})
endif()

add_library(SketchSort STATIC
        SketchSort.cpp
        )



