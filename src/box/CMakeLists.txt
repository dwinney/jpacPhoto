# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

if (BUILD_DICTIONARY)
    ROOT_GENERATE_DICTIONARY(G__jpacBox
            box_amplitude.hpp
            box_discontinuity.hpp
            projected_amplitude.hpp
            helicity_PWA.hpp
            interpolation_2D.hpp
            interpolated_discontinuity.hpp
            LINKDEF jpacBoxLinkDef.h
            )
    add_library(jpacBox SHARED ${SRC} G__jpacBox.cxx)
else()
    add_library(jpacBox SHARED ${SRC})
endif()

target_link_libraries( jpacBox ROOT::Core ROOT::MathMore )
target_link_libraries( jpacBox jpacPhoto )
target_link_libraries( jpacBox Boost::boost)

install(TARGETS jpacBox
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
if (BUILD_DICTIONARY)
install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libjpacBox_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()