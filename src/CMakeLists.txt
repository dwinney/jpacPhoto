# BUILD LIBRARY FROM LOCAL FiLES FOR THE MAIN LIBRARY
file(GLOB SRC "core/*.cpp" "covariants/*.cpp" "plotting/*.cpp") 
file(GLOB JPACTOOLS_FILES "AmpTools/*.cpp")

include_directories(core)
include_directories(covariants)
include_directories(plotting)

if (JPACTOOLS)
    include_directories(AmpTools)
    include_directories(${AMPTOOLS_SRC})
    # AmpTools by default is built as a static library
    # This will build a local copy in a dynamic library
    add_library(AMPTOOLS SHARED ${STATIC_AMPTOOLS})
    set_target_properties(AMPTOOLS PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(AMPTOOLS PUBLIC
      -Wl,--whole-archive
      ${STATIC_AMPTOOLS}
      -Wl,--no-whole-archive
    )
    install(TARGETS AMPTOOLS
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
endif()

ROOT_GENERATE_DICTIONARY(G__JPACPHOTO
    amplitude.hpp
    partial_wave.hpp
    kinematics.hpp
    dirac_spinor.hpp
    dirac_matrix.hpp
    lorentz_tensor.hpp
    fitter.hpp
    data_set.hpp
    inclusive_process.hpp
    LINKDEF JPACPHOTOLinkDef.h
)
add_library(JPACPHOTO SHARED ${SRC} G__JPACPHOTO.cxx )

target_link_libraries( JPACPHOTO ROOT::Core ROOT::MathMore )
target_link_libraries( JPACPHOTO Boost::boost)

install(TARGETS JPACPHOTO
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libJPACPHOTO_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Build the executable to integrate with root cling
add_executable( jpacPhoto "cling/jpacPhoto.cpp" )
target_link_libraries( jpacPhoto ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( jpacPhoto JPACPHOTO )
if (JPACTOOLS)
    target_link_libraries( jpacPhoto ${AMPTOOLS} )
endif()