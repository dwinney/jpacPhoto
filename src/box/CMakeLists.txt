# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

ROOT_GENERATE_DICTIONARY(G__jpacBox
        box_amplitude.hpp
        box_discontinuity.hpp
        projected_amplitude.hpp
        LINKDEF jpacBoxLinkDef.h
        )

add_library(jpacBox SHARED ${SRC} G__jpacBox.cxx)
target_link_libraries( jpacBox ROOT::Core ROOT::MathMore )
target_link_libraries( jpacBox jpacPhoto )
target_link_libraries( jpacBox Boost::boost)

install(TARGETS jpacBox
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
        
install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libjpacBox_rdict.pcm    
DESTINATION ${CMAKE_INSTALL_LIBDIR})