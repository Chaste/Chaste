project(cell_based)
#The compiler flag /bigobj deals with the "fatal error C1128": 'number of sections exceeded object file format limit'
if (MSVC)
    add_definitions(-bigobj)
endif()

chaste_do_component(cell_based)
