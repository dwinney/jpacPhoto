# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

ROOT_GENERATE_DICTIONARY(G__jpacPhoto
        amplitude.hpp
        amplitude_sum.hpp
        baryon_resonance.hpp
        dirac_exchange.hpp
        pomeron_exchange.hpp
        primakoff_effect.hpp
        pseudoscalar_exchange.hpp
        vector_exchange.hpp
        LINKDEF jpacPhotoLinkDef.h
        )

add_library(jpacPhoto SHARED ${SRC} G__jpacPhoto.cxx )
target_link_libraries( jpacPhoto ROOT::Core ROOT::MathMore )

install(TARGETS jpacPhoto
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libjpacPhoto_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})