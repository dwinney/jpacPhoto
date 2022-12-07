# BUILD LIBRARY FROM LOCAL FiLES FOR THE MAIN LIBRARY
file(GLOB SRC "core/*.cpp" #[[ "amplitudes/*.cpp" ]]) 
file(GLOB INC "core/*.hpp" #[[ "amplitudes/*.hpp" ]]) 

include_directories(core)
# include_directories(amplitudes)
# include_directories(covariants)

add_library(JPACPHOTO SHARED ${SRC})

target_link_libraries( JPACPHOTO ROOT::Core ROOT::MathMore )

install(TARGETS JPACPHOTO
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )


# Build the executable to integrate with root cling
add_executable( jpacPhoto "cling/jpacPhoto.cpp" )
target_link_libraries( jpacPhoto ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( jpacPhoto JPACPHOTO )
target_link_libraries( jpacPhoto ${JSTYLELIB})

# if (BUILD_INCLUSIVE)
#     target_link_libraries( ${exename} jpacInclusive)
# endif()

# if (BUILD_BOX)
#     target_link_libraries( ${exename} jpacBox)
# endif()