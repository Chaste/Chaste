#!/bin/bash -ex

mkdir -p ${PREFIX}/build
cd ${PREFIX}/build || exit

# Modify pip settings for internal Chaste Python env
export PIP_NO_DEPENDENCIES="False"
export PIP_NO_INDEX="False"

# Configure
cmake \
  -DChaste_ENABLE_PYCHASTE=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_LIBRARY_PATH="${PREFIX}/lib" \
  -DBUILD_SHARED_LIBS=ON \
  -DBOOST_ROOT="${PREFIX}" \
  -DHDF5_C_COMPILER_EXECUTABLE="${PREFIX}/bin/h5pcc" \
  -DPETSC_DIR="${PREFIX}" \
  -DPYTHON_EXECUTABLE=${PYTHON} \
  -DVTK_DIR=${PREFIX} \
  -DXERCESC_INCLUDE="${PREFIX}/include" \
  -DXERCESC_LIBRARY="${PREFIX}/lib/libxerces-c.so" \
  -DXSD_EXECUTABLE="${PREFIX}/bin/xsd" \
  $SRC_DIR

# Revert pip settings
export PIP_NO_DEPENDENCIES="True"
export PIP_NO_INDEX="True"

# Build
make -j ${CPU_COUNT} pychaste

# Install
${PYTHON} -m pip install -v pychaste/package --prefix="${PREFIX}"

# Cleanup
rm -rf \
  cell_based/CMakeFiles \
  global/CMakeFiles \
  io/CMakeFiles \
  linalg/CMakeFiles \
  mesh/CMakeFiles \
  ode/CMakeFiles \
  pde/CMakeFiles \
  python \
  pychaste/CMakeFiles \
  pychaste/package
