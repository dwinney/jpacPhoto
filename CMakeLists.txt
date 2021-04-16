cmake_minimum_required (VERSION 2.9)

project (jpacPhoto)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# Make sure gcc version is atleast 5!
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    message(FATAL_ERROR "GCC version must be atleast 5.0!")
  endif()
endif()

# INSTALLATION SETTINGS
set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

# Make sure we can find ROOT
execute_process(COMMAND root-config --prefix 
                COMMAND tr -d \\n
                RESULT_VARIABLE RC_RESULT 
                OUTPUT_VARIABLE ROOT_PREFIX)

list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")

##-----------------------------------------------------------------------
## Build the base library jpacPhoto

## Connect ROOT
find_package(ROOT REQUIRED MathMore)
if (ROOT_FOUND)
    message("-- ROOT found!")
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "-- ROOT not found!")
endif()

# BUILD LIBRARY FROM LOCAL FiLES
include_directories("include")
include_directories("src")

file(GLOB INC "include/*.hpp" "include/amplitudes/*.hpp")
file(GLOB SRC "src/*.cpp"     "src/amplitudes/*.cpp")

add_library( jpacPhoto SHARED ${INC} ${SRC} )
target_link_libraries( jpacPhoto ${ROOT_LIBRARIES})

##-----------------------------------------------------------------------
## Look for BOOSt and if found, build the auxiliary library jpacBox

## Connect BOOST
find_package(Boost)
if (Boost_FOUND)
    message("-- Building jpacBox")
    include_directories(${BOOST_INCLUDE_DIRS})
    link_directories(${BOOST_LIBRARY_DIRS})
    file(GLOB BOX_INC "include/box/*.hpp")
    file(GLOB BOX_SRC "src/box/*.cpp")

    add_library( jpacBox SHARED ${BOX_INC} ${BOX_SRC} )
    target_link_libraries( jpacBox ${ROOT_LIBRARIES})
    target_link_libraries( jpacBox jpacPhoto)
    target_link_libraries( jpacBox ${Boost_LIBRARIES})
else()
    message("-- BOOST not found! jpacBox will not be available.")
endif()

##-----------------------------------------------------------------------
## Installation

if (Boost_FOUND)
    install(TARGETS jpacBox jpacPhoto
            LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}" )
else()
    install(TARGETS jpacPhoto
            LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}" )
endif()

##-----------------------------------------------------------------------
## Look for jpacStyle if found, build all the executables

# Find the jpacStyle library
find_library(JSTYLELIB NAMES jpacStyle libjpacStyle 
                       HINTS "$ENV{JPACSTYLE}/lib")
if (JSTYLELIB)
    message("-- jpacStyle found! (${JSTYLELIB})")
    message("-- Building executables")
    include_directories("$ENV{JPACSTYLE}/include")
    link_libraries(${JSTYLELIB})

    include_directories("executables")
    if (Boost_FOUND)
        file(GLOB EXE_FILES "executables/pentaquark/*.cpp" "executables/XYZ_plots/*.cpp" "executables/open_charm/*.cpp")
    else()
        file(GLOB EXE_FILES "executables/pentaquark/*.cpp" "executables/XYZ_plots/*.cpp")
    endif()

    foreach( exefile ${EXE_FILES} )
        get_filename_component( exename ${exefile} NAME_WE)
        add_executable( ${exename} ${exefile} )
        target_link_libraries( ${exename} jpacPhoto)
        if (Boost_FOUND)
            target_link_libraries( ${exename} jpacBox)
        endif()
        target_link_libraries( ${exename} ${ROOT_LIBRARIES})
        target_link_libraries( ${exename} ${JSTYLELIB})
    endforeach( exefile ${EXE_FILES} )
else()
    message("-- jpacStyle not found! Executables will not be installed.")
endif()

