cmake_minimum_required( VERSION 3.2 FATAL_ERROR )
set(CMAKE_CXX_STANDARD 17)
project(Significance)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package( ROOT COMPONENTS MathCore Core RooStats )
message(STATUS "Using ROOT From: ${ROOT_INCLUDE_DIRS}")
include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__Significance Significance.h LINKDEF LinkDef.h)

add_library(Significance SHARED Significance.cpp G__Significance.cxx)
target_link_libraries(Significance ${ROOT_LIBRARIES})



install(FILES Significance.h DESTINATION include)

install( FILES ${CMAKE_BINARY_DIR}/libSignificance.rootmap
        ${CMAKE_BINARY_DIR}/libSignificance_rdict.pcm
        DESTINATION lib
        COMPONENT libraries)

install(TARGETS Significance
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        COMPONENT libraries)

