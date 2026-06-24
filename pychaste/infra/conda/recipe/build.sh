#!/bin/bash -ex

mkdir -p ${PREFIX}/build
cd ${PREFIX}/build || exit

# CMakeLists.txt creates a venv and pip-installs chaste_codegen and cppwg from PyPI.
# Unset conda-build's pip isolation flags so those installs can reach the internet.
unset PIP_NO_INDEX
unset PIP_NO_DEPENDENCIES

# Configure
cmake ${CMAKE_ARGS} \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -w" \
  -DX11_X11_INCLUDE_PATH="${PREFIX}/include" \
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

# Build. The pychaste wrapper translation units and their LTO link are very
# memory-heavy, so CPU_COUNT is kept low on memory-constrained machines (set by
# build-package.sh and forwarded through recipe.yaml) to avoid exhausting RAM.
make -j "${CPU_COUNT}" pychaste

# Install. --no-deps: the runtime dependencies (matplotlib, numpy, xvfbwrapper)
# are provided as conda run requirements, so pip must not pull them from PyPI
# (where e.g. building matplotlib from an sdist needs meson-python).
${PYTHON} -m pip install -v --no-deps pychaste/package --prefix="${PREFIX}"

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
