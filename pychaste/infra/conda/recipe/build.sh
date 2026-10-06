#!/bin/bash -ex

mkdir -p ${PREFIX}/build
cd ${PREFIX}/build || exit

# CMakeLists.txt creates a venv and pip-installs chaste_codegen and cppwg from PyPI.
# Unset conda-build's pip isolation flags so those installs can reach the internet.
unset PIP_NO_INDEX
unset PIP_NO_DEPENDENCIES

# Platform-specific configuration. On Linux the shared library suffix is .so and
# VTK's find_package(X11) needs the X11 headers; on macOS the suffix is .dylib
# and VTK renders through Cocoa, so no X11 include path is required.
if [[ "$(uname)" == "Darwin" ]]; then
  shlib_ext="dylib"
  platform_args=()
  # The Fortran compiler (compiler('fortran')) pulls in conda's gcc, whose
  # activation sets CC to gcc while the C++ compiler is clang. Force the matching
  # clang C wrapper so the whole configure uses a consistent clang toolchain.
  platform_args+=(-DCMAKE_C_COMPILER="${HOST}-clang")
  # Point cmake straight at conda's VTK config so find_package locates it
  # irrespective of macOS find-root behaviour.
  for _vtk_dir in "${PREFIX}"/lib/cmake/vtk-*; do
    [ -d "${_vtk_dir}" ] && platform_args+=(-DVTK_DIR="${_vtk_dir}")
  done
else
  shlib_ext="so"
  platform_args=(-DX11_X11_INCLUDE_PATH="${PREFIX}/include")
fi

# Configure.
# Disable error-on-warning for the package build (both linux-64 and osx-64).
# conda's CFLAGS carry gcc-only options (e.g. -fno-merge-constants) that clang
# only warns about; combined with Chaste's default -Werror, such a warning
# becomes an error and breaks configure-time C probes - notably find_package
# (Threads), which VTK depends on, yielding a spurious "Did not find VTK" on
# macOS. CXX warnings are silenced with -w; turning off -Werror also covers the
# C probes.
cmake ${CMAKE_ARGS} \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -w" \
  -DChaste_ERROR_ON_WARNING=OFF \
  "${platform_args[@]}" \
  -DChaste_ENABLE_PYCHASTE=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_C_COMPILER_EXECUTABLE="${PREFIX}/bin/h5pcc" \
  -DPETSC_DIR="${PREFIX}" \
  -DPython3_EXECUTABLE="${PYTHON}" \
  -DXERCESC_INCLUDE="${PREFIX}/include" \
  -DXERCESC_LIBRARY="${PREFIX}/lib/libxerces-c.${shlib_ext}" \
  -DXSD_EXECUTABLE="${PREFIX}/bin/xsd" \
  $SRC_DIR

# Build. The pychaste wrapper translation units and their LTO link are very
# memory-heavy, so CPU_COUNT is kept low on memory-constrained machines (set by
# the build-package-* scripts and forwarded through recipe.yaml) to avoid exhausting RAM.
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
