# Build the executable to integrate with root cling
file(GLOB EXE_FILES "*.cpp")
foreach( exefile ${EXE_FILES} )
    get_filename_component( exename ${exefile} NAME_WE)
    add_executable( ${exename} ${exefile} )
    target_link_libraries( ${exename} ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
    target_link_libraries( ${exename} jpacPhoto)
    target_link_libraries( ${exename} ${JSTYLELIB})

    if (BUILD_INCLUSIVE)
        target_link_libraries( ${exename} jpacInclusive)
    endif()

    if (BUILD_BOX)
        target_link_libraries( ${exename} jpacBox)
    endif()

endforeach( exefile ${EXE_FILES} )