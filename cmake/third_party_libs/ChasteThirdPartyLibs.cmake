#URLs to Third party libraries needed by Chaste

# Specify the urls of the libraries you want to build separated by spaces and/or newlines, or as separate strings.
# Note that the URLS of PARMETIS, METIS, F2CBLAS, F2CLAPACK are all automatically obtained from the
# PETSc distribution once it has been downloaded and unzipped. So, there is no need to manually
# specify the URLs for these libraries 
set(PETSC_URLS "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p6.tar.gz")
set(HDF5_URLS "http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.10-patch1/src/hdf5-1.8.10-patch1.tar.gz")
set(SUNDIALS_URLS "https://computation.llnl.gov/casc/sundials/download/code/sundials-2.5.0.tar.gz")
set(BOOST_URLS "http://kent.dl.sourceforge.net/project/boost/boost/1.53.0/boost_1_53_0.tar.gz" 
               #"http://kent.dl.sourceforge.net/project/boost/boost/1.52.0/boost_1_52_0.tar.gz" 
            )
set(VTK_URLS "http://www.vtk.org/files/release/5.8/vtk-5.8.0.tar.gz")
