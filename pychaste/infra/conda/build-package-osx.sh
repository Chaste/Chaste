#!/bin/bash -e

# Build a PyChaste conda package for macOS natively with rattler-build.
#
# Unlike build-package-linux.sh, which builds linux-64 inside a conda-forge Docker
# container, this script runs directly on a macOS host: building macOS conda
# packages requires the macOS SDK and the host toolchain, so there is no
# container to run it in.
#
# --target selects the macOS architecture and must match the host: use
# --target=osx-64 on an Intel Mac and --target=osx-arm64 on an Apple Silicon Mac.
#
# Example usage (local source, osx-64 on an Intel Mac):
#   ./build-package-osx.sh --variant=osx_64_python3.12_cpython \
#       --target=osx-64 --source-dir=/path/to/Chaste --cpu-count=4
#
# Example usage (local source, osx-arm64 on an Apple Silicon Mac):
#   ./build-package-osx.sh --variant=osx_arm64_python3.12_cpython \
#       --target=osx-arm64 --source-dir=/path/to/Chaste --cpu-count=4
#
# Example usage (clone from GitHub):
#   ./build-package-osx.sh --variant=osx_64_python3.12_cpython --branch=2026.1

# Parse args
variant=
branch=develop
cpu_count=
source_dir=
target=osx-64

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
  --target=*)
    target=$(expr "x$option" : "x--target=\(.*\)")
    ;;
  *)
    echo "Unknown option: $option" 1>&2
    exit 1
    ;;
  esac
done

if [ -z "${variant}" ]; then
  echo "Error: --variant is required" 1>&2
  echo "Usage: $(basename "$0") --variant=<name> [--target=osx-64|osx-arm64] [--branch=<branch>] [--cpu-count=<n>] [--source-dir=<path>]" 1>&2
  exit 1
fi

case "${target}" in
  osx-64 | osx-arm64) ;;
  *)
    echo "Error: --target must be osx-64 or osx-arm64 (got '${target}')" 1>&2
    exit 1
    ;;
esac

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

# Install rattler-build as a self-contained static binary (no conda solve),
# matching the host architecture (arm64 on Apple Silicon, x86_64 on Intel).
case "$(uname -m)" in
  arm64 | aarch64) rb_arch=aarch64 ;;
  *)               rb_arch=x86_64  ;;
esac
RATTLER_BUILD=/tmp/rattler-build-osx
curl -fsSL \
  "https://github.com/prefix-dev/rattler-build/releases/latest/download/rattler-build-${rb_arch}-apple-darwin" \
  -o "${RATTLER_BUILD}"
chmod +x "${RATTLER_BUILD}"

# macOS SDK. conda-forge's clang activation expects CONDA_BUILD_SYSROOT to point
# at an SDK matching the deployment target pinned in the variant. The host only
# ships current SDKs, so download and cache the required SDK if it is absent.
# osx-64 targets 10.15 (the minimum providing aligned_alloc, used by pocketfft);
# osx-arm64 targets 11.0, the floor for Apple Silicon. Both variables are
# exported here (for rattler-build's environment activation) and forwarded into
# the build script via recipe.yaml.
case "${target}" in
  osx-64)    sdk_version=10.15; deployment=10.15 ;;
  osx-arm64) sdk_version=11.3;  deployment=11.0  ;;
esac
SDK_CACHE="${FEEDSTOCK_ROOT}/.osx-sdk"
export CONDA_BUILD_SYSROOT="${SDK_CACHE}/MacOSX${sdk_version}.sdk"
export MACOSX_DEPLOYMENT_TARGET="${deployment}"
if [ ! -d "${CONDA_BUILD_SYSROOT}" ]; then
  mkdir -p "${SDK_CACHE}"
  curl -fsSL \
    "https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX${sdk_version}.sdk.tar.xz" \
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
  --target-platform "${target}" \
  --output-dir "${OUTPUT_DIR}" \
  --channel conda-forge \
  --channel pychaste

# Copy the built package to the host build_artifacts directory
mkdir -p "${FEEDSTOCK_ROOT}/build_artifacts/${target}"
cp "${OUTPUT_DIR}/${target}"/chaste-*.conda "${FEEDSTOCK_ROOT}/build_artifacts/${target}/"
