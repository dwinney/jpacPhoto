cmake_minimum_required (VERSION 3.0)

project (jpacPhoto)

set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# Set output lib and bin locations
set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_BINDIR})

# Make sure executables are always built with the install RPATH 
# else exes get linked to the build library and not the installed one
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_CURRENT_SOURCE_DIR}/lib")
endif("${isSystemDir}" STREQUAL "-1")

##-----------------------------------------------------------------------
## Finding and connecting ROOT

# Make sure we can find ROOT 
execute_process(COMMAND root-config --prefix 
                COMMAND tr -d \\n
                RESULT_VARIABLE RC_RESULT 
                OUTPUT_VARIABLE ROOT_PREFIX)
list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS MathMore)
if (ROOT_FOUND)
    message("-- ROOT found! (${ROOT_USE_FILE})")
    #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
    include(${ROOT_USE_FILE})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
else()
    message(SEND_ERROR "-- Root not found!")
endif()

##-----------------------------------------------------------------------
## Look for BOOST needed for jpacBox

set(BOOST_MIN_VERSION 1.68.0) # Need complex Gauss-Kronrod integration
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
if (Boost_FOUND)
    message("-- Building jpacBox")
else()
    message("-- BOOST not found! jpacBox will not be available.")
endif()

# Build library 
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(src/core)