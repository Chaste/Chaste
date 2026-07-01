#!/bin/bash -e

# Build a PyChaste conda package for macOS (osx-64) natively with rattler-build.
#
# Unlike build-package-linux.sh, which builds linux-64 inside a conda-forge Docker
# container, this script runs directly on a macOS host: cross-building macOS
# conda packages requires the macOS SDK and the host toolchain, so there is no
# container to run it in.
#
# Example usage (local source):
#   ./build-package-osx.sh --variant=osx_64_python3.12_cpython \
#       --source-dir=/path/to/Chaste --cpu-count=4
#
# Example usage (clone from GitHub):
#   ./build-package-osx.sh --variant=osx_64_python3.12_cpython --branch=2026.1

# Parse args
variant=
branch=develop
cpu_count=
source_dir=

for option; do
  case $option in
  --variant=*)
    variant=$(expr "x$option" : "x--variant=\(.*\)")
    ;;
  --source-dir=*)
    source_dir=$(expr "x$option" : "x--source-dir=\(.*\)")
    ;;
  --branch=*)
    branch=$(expr "x$option" : "x--branch=\(.*\)")
    ;;
  --cpu-count=*)
    cpu_count=$(expr "x$option" : "x--cpu-count=\(.*\)")
    ;;
  *)
    echo "Unknown option: $option" 1>&2
    exit 1
    ;;
  esac
done

if [ -z "${variant}" ]; then
  echo "Error: --variant is required" 1>&2
  echo "Usage: $(basename "$0") --variant=<name> [--branch=<branch>] [--cpu-count=<n>] [--source-dir=<path>]" 1>&2
  exit 1
fi

if [ "$(uname)" != "Darwin" ]; then
  echo "Error: build-package-osx.sh must be run on macOS; use build-package-linux.sh for linux-64." 1>&2
  exit 1
fi

set -x

# Configure environment
export FEEDSTOCK_ROOT="$(cd "$(dirname "$0")" && pwd)"
export RECIPE_ROOT="${FEEDSTOCK_ROOT}/recipe"
export CONFIG_FILE="${FEEDSTOCK_ROOT}/variants/${variant}.yaml"
# CPU_COUNT bounds parallelism for both rattler-build and the C++ compile
# (make -j). Keep it low on memory-constrained machines: the pychaste wrapper
# translation units and their LTO link consume roughly 3 GB each. It is
# forwarded into the sanitized build script environment via recipe.yaml.
export CPU_COUNT="${cpu_count:-$(sysctl -n hw.ncpu)}"
export PYTHONUNBUFFERED=1

OUTPUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/rattler_output.XXXXXX")"

# Install rattler-build as a self-contained static binary (no conda solve).
RATTLER_BUILD=/tmp/rattler-build-osx
curl -fsSL \
  https://github.com/prefix-dev/rattler-build/releases/latest/download/rattler-build-x86_64-apple-darwin \
  -o "${RATTLER_BUILD}"
chmod +x "${RATTLER_BUILD}"

# macOS SDK. conda-forge's clang activation expects CONDA_BUILD_SYSROOT to point
# at the SDK matching MACOSX_SDK_VERSION (10.13, pinned in the variant). The host
# only ships current SDKs, so download and cache the 10.13 SDK if it is absent.
# Both variables are exported here (for rattler-build's environment activation)
# and forwarded into the build script via recipe.yaml.
SDK_CACHE="${FEEDSTOCK_ROOT}/.osx-sdk"
export CONDA_BUILD_SYSROOT="${SDK_CACHE}/MacOSX10.15.sdk"
export MACOSX_DEPLOYMENT_TARGET=10.15
if [ ! -d "${CONDA_BUILD_SYSROOT}" ]; then
  mkdir -p "${SDK_CACHE}"
  curl -fsSL \
    https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX10.15.sdk.tar.xz \
    | tar -xJ -C "${SDK_CACHE}"
fi

# Show configuration
cat "${CONFIG_FILE}"
"${RATTLER_BUILD}" --version

# Get source code as a shallow (depth 1) clone so the full git history is not
# dragged into the build.
if [ -n "${source_dir}" ]; then
  # The local source tree may be owned differently; mark it safe so git works.
  git config --global --add safe.directory "${source_dir}"
  # Determine the version from the full local repo (it has the tags and history
  # the shallow clone lacks) before making the shallow clone. file:// is required
  # for git to honour --depth when cloning a local path.
  export CHASTE_VERSION=$(bash "${FEEDSTOCK_ROOT}/package-version.sh" "${source_dir}")
  rm -rf /tmp/Chaste
  git clone --recurse-submodules --depth 1 "file://${source_dir}" /tmp/Chaste
else
  rm -rf /tmp/Chaste
  git clone --recursive --branch "${branch}" --depth 1 --tags https://github.com/Chaste/Chaste.git /tmp/Chaste
  export CHASTE_VERSION=$(bash "${FEEDSTOCK_ROOT}/package-version.sh" /tmp/Chaste)
fi
export CHASTE_SOURCE_DIR=/tmp/Chaste

"${RATTLER_BUILD}" build \
  --recipe "${RECIPE_ROOT}/recipe.yaml" \
  --variant-config "${CONFIG_FILE}" \
  --target-platform osx-64 \
  --output-dir "${OUTPUT_DIR}" \
  --channel conda-forge \
  --channel pychaste

# Copy the built package to the host build_artifacts directory
mkdir -p "${FEEDSTOCK_ROOT}/build_artifacts/osx-64"
cp "${OUTPUT_DIR}"/osx-64/chaste-*.conda "${FEEDSTOCK_ROOT}/build_artifacts/osx-64/"
