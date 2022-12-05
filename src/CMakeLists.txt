# BUILD LIBRARY FROM LOCAL FiLES FOR THE MAIN LIBRARY
file(GLOB SRC "core/*.cpp" "amplitudes/*.cpp" "covariants/*.cpp") 
file(GLOB INC "core/*.hpp" "amplitudes/*.hpp" "covariants/*.hpp") 

include_directories(core)
include_directories(amplitudes)
include_directories(covariants)

add_library(JPACPHOTO SHARED ${SRC})

target_link_libraries( JPACPHOTO ROOT::Core ROOT::MathMore )

install(TARGETS JPACPHOTO
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )