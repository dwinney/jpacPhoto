# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

ROOT_GENERATE_DICTIONARY(G__jpacInclusive
        inclusive_kinematics.hpp
        inclusive_production.hpp
        sigma_JPAC.hpp
        sigma_PDG.hpp
        triple_regge.hpp
        LINKDEF jpacInclusiveLinkDef.h
        )

add_library(jpacInclusive SHARED ${SRC} G__jpacInclusive.cxx )
target_link_libraries( jpacInclusive ROOT::Core ROOT::MathMore )
target_link_libraries( jpacInclusive jpacPhoto )

install(TARGETS jpacInclusive
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libjpacInclusive_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})