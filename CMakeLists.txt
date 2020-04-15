cmake_minimum_required (VERSION 2.8.11)

macro(use_cxx11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(use_cxx11)

project(dvm)
use_cxx11()

include_directories(include)
include_directories(logging)

set(SOURCES src/Cell.cpp 
            src/CollisionPartnersGenerator.cpp
            src/distribution.cpp 
            src/Grid.cpp 
            src/interpolation.cpp 
            src/PostCollisionVelocitiesGenerator.cpp
            src/collisions.cpp
            src/simulate.cpp
            logging/easylogging++.cc)

add_library(dvmLibrary ${SOURCES})
add_executable(dvm src/main.cpp)
target_link_libraries(dvm dvmLibrary)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_executable(dvmTest test/collisionPartnersGenerationTest.cpp
                        test/interpolationTest.cpp)
target_link_libraries(dvmTest dvmLibrary)
add_definitions(-DBOOST_TEST_DYN_LINK) 
target_link_libraries(dvmTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

