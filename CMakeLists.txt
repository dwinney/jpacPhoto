cmake_minimum_required (VERSION 2.9)

project (jpacPhoto)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# Make sure gcc version is atleast 5!
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
    message(FATAL_ERROR "GCC version must be atleast 7.0!")
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
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
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
## build the auxiliary library jpacInclusive

file(GLOB INC_SI  "include/inclusive/*.hpp")
file(GLOB SRC_SI  "src/inclusive/*.cpp")

add_library( jpacInclusive SHARED ${INC_SI} ${SRC_SI} )
target_link_libraries( jpacInclusive ${ROOT_LIBRARIES})
target_link_libraries( jpacInclusive jpacPhoto)

##-----------------------------------------------------------------------
## Installation

if (Boost_FOUND)
    install(TARGETS jpacBox jpacPhoto jpacInclusive
            LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}" )
else()
    install(TARGETS jpacPhoto jpacInclusive
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

    ## Sort different executables into three categories:
    ## Exclusive XYZP plots
    ## Inclusive XYZP plots
    ## Box amplitude plots

    include_directories("executables")
    file(GLOB EXC_EXE_FILES "executables/pentaquark/*.cpp" "executables/exclusive_XYZ/*.cpp")
    file(GLOB INC_EXE_FILES "executables/inclusive_XYZ/*.cpp")
    file(GLOB BOX_EXE_FILES "executables/open_charm/*.cpp")

    foreach( exefile ${EXC_EXE_FILES} )
        get_filename_component( exename ${exefile} NAME_WE)
        add_executable( ${exename} ${exefile} )
        target_link_libraries( ${exename} ${ROOT_LIBRARIES})
        target_link_libraries( ${exename} jpacPhoto)
        target_link_libraries( ${exename} ${JSTYLELIB})
    endforeach( exefile ${EXE_FILES} )

    foreach( exefile ${INC_EXE_FILES} )
        get_filename_component( exename ${exefile} NAME_WE)
        add_executable( ${exename} ${exefile} )
        target_link_libraries( ${exename} ${ROOT_LIBRARIES})
        target_link_libraries( ${exename} jpacInclusive)
        target_link_libraries( ${exename} ${JSTYLELIB})
    endforeach( exefile ${EXE_FILES} )

    if (Boost_FOUND)
        foreach( exefile ${BOX_EXE_FILES} )
            get_filename_component( exename ${exefile} NAME_WE)
            add_executable( ${exename} ${exefile} )
            target_link_libraries( ${exename} ${ROOT_LIBRARIES})
            target_link_libraries( ${exename} jpacPhoto)
            target_link_libraries( ${exename} jpacBox)
            target_link_libraries( ${exename} ${JSTYLELIB})
        endforeach( exefile ${EXE_FILES} )
    endif()

else()
    message("-- jpacStyle not found! Executables will not be installed.")
endif()

