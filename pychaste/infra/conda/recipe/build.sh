#!/bin/bash -ex

mkdir -p ${PREFIX}/build
cd ${PREFIX}/build || exit

# CMakeLists.txt creates a venv and pip-installs chaste_codegen and cppwg from PyPI.
# Unset conda-build's pip isolation flags so those installs can reach the internet.
unset PIP_NO_INDEX
unset PIP_NO_DEPENDENCIES

# Configure
cmake ${CMAKE_ARGS} \
  -DChaste_ENABLE_PYCHASTE=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_C_COMPILER_EXECUTABLE="${PREFIX}/bin/h5pcc" \
  -DPETSC_DIR="${PREFIX}" \
  -DPython3_EXECUTABLE="${PYTHON}" \
  -DXERCESC_INCLUDE="${PREFIX}/include" \
  -DXERCESC_LIBRARY="${PREFIX}/lib/libxerces-c.so" \
  -DXSD_EXECUTABLE="${PREFIX}/bin/xsd" \
  $SRC_DIR

# Build
make -j ${CPU_COUNT} pychaste

# Install
${PYTHON} -m pip install -v pychaste/package --prefix="${PREFIX}"

# Remove build artifacts so conda-build's rpath scanner doesn't trip over
# object files, cmake probe binaries, and test executables.
find . -type d -name "CMakeFiles" -prune -exec rm -rf {} +
find . -type d -name "test" -prune -exec rm -rf {} +
find . -type f -name "*.o" -delete
rm -rf \
  chaste_python3_venv \
  python \
  pychaste/package \
  pychaste/wrappers
