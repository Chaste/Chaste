find_package(Chaste COMPONENTS ${Chaste_DEPENDS_mesh})
#Define TETLIBRARY: prevents Tetgen from including main() method
add_definitions(-DTETLIBRARY)
add_definitions(-DANSI_DECLARATORS -DTRILIBRARY)

chaste_do_component(mesh)
set_source_files_properties(src/3rdparty/tetgen1.4.2/predicates.cpp 
    PROPERTIES COMPILE_FLAGS "-Wno-error -O0")
set_source_files_properties(src/3rdparty/tetgen1.4.2/tetgen.cpp 
    PROPERTIES COMPILE_FLAGS "-Wno-error -O0")
set_source_files_properties(src/3rdparty/triangle/triangle.cpp
        PROPERTIES COMPILE_FLAGS "-Wno-error -O0")
