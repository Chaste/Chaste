list(APPEND Chaste_LINK_LIBRARIES "fftw3")
list(APPEND Chaste_LINK_LIBRARIES "fftw3_threads")

find_package(Chaste COMPONENTS cell_based)
chaste_do_project(ImmersedBoundary)