cmake_minimum_required(VERSION 3.7)


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/CapR)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/SketchSort)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty/json11)


add_executable(RaptRanker
        main.cpp
        RaptRanker.cpp
        )


target_link_libraries(RaptRanker
        CapR
        SketchSort
        SQLiteCpp
        sqlite3
        pthread
        dl
        json11
        )

target_include_directories(RaptRanker PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/CapR
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/SketchSort
        ${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp/include/SQLiteCpp
        ${CMAKE_CURRENT_LIST_DIR}/thirdparty/SQLiteCpp/include/json11
        )

find_package(Boost)
IF (Boost_FOUND)
   target_include_directories(RaptRanker PUBLIC ${Boost_INCLUDE_DIR})
endif()