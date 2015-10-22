cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#add CFLAGS if not already present.
if(${CMAKE_VERSION} VERSION_GREATER 3.1.0 OR ${CMAKE_VERSION} VERSION_EQUAL 3.1.0)
  set(CMAKE_CXX_STANDARD 14)
  set(CXX_STANDARD_REQUIRED 11)
else()
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "c\\+\\+[01][xyz147]")
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
# XCode isn't picking up the standard set above.
if(CMAKE_GENERATOR STREQUAL Xcode)
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
endif()

find_package(NLOPT)
if(NOT TARGET SpinWaveGenie)
  find_package(SpinWaveGenie)
endif()

add_executable(VillainSC1Dispersion SC1Dispersion.cpp)
target_link_libraries( VillainSC1Dispersion SpinWaveGenie ${NLOPT_LIBRARIES})
target_include_directories(VillainSC1Dispersion SYSTEM PRIVATE ${NLOPT_INCLUDE_DIRECTORY})

add_executable(VillainFMDispersion FMDispersion.cpp)
target_link_libraries( VillainFMDispersion SpinWaveGenie)

file( COPY "SC1Dispersion.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

