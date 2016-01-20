list(APPEND Chaste_LINK_LIBRARIES "fftw3")
list(APPEND Chaste_LINK_LIBRARIES "fftw3_threads")
list(APPEND Chaste_LINK_LIBRARIES "boost_program_options")

find_package(Chaste COMPONENTS cell_based)
chaste_do_project(ImmersedBoundary)