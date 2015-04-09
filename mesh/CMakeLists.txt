find_package(Chaste COMPONENTS ${CHASTE_DEPENDS_mesh})
#Define TETLIBRARY: prevents Tetgen from including main() method
add_definitions(-DTETLIBRARY)
add_definitions(-DANSI_DECLARATORS -DTRILIBRARY)

chaste_do_component(mesh)
