cmake_minimum_required(VERSION 3.0)

#add CFLAGS if not already present.

if(${CMAKE_VERSION} VERSION_GREATER 3.1.0 OR ${CMAKE_VERSION} VERSION_EQUAL 3.1.0)
  set(CMAKE_CXX_STANDARD 14)
  set(CXX_STANDARD_REQUIRED 11)
else()
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+[01][xyz147]")
    message(STATUS "compiler detected as: GNU compiler")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    if(COMPILER_SUPPORTS_CXX14)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    elseif(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
  endif()
endif()

if(NOT TARGET SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()

add_executable(FMDispersion Dispersion.cpp)
target_link_libraries( FMDispersion SpinWaveGenie)
add_executable(FMcut TwoDimensionalCut.cpp)
target_link_libraries( FMcut SpinWaveGenie)
add_executable(FMPowderAverage PowderAverage.cpp)
target_link_libraries( FMPowderAverage SpinWaveGenie)

file( COPY "FM.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file( COPY "FMavg.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
