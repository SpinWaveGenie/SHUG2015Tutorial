cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#add CFLAGS if not already present.
if(CMAKE_COMPILER_IS_GNUCXX AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+0x")
    message(STATUS "compiler detected as: GNU compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse4")
elseif(APPLE AND NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+11")
    message(STATUS "Setting up Xcode for C++11 with libc++.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -msse4")
endif()

if(NOT TARGET SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()

add_executable(VillainSC2Dispersion SC2Dispersion.cpp)
target_link_libraries( VillainSC2Dispersion SpinWaveGenie)
