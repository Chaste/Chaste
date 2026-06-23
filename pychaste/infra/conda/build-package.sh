#!/bin/bash -e

# Build a PyChaste conda package inside a conda-forge build container.
#
# Example usage (clone from GitHub):
#   ./build-package.sh --variant=linux_64_python3.10_cpython --branch=2026.1 --parallel=4
#
# Example usage (local source):
#   ./build-package.sh --variant=linux_64_python3.10_cpython --source-dir=/chaste
#
# Intended for use in a build container e.g.:
#   docker run --rm -it \
#     -v $(pwd):/home/conda \
#     -e HOST_USER_ID="$(id -u)" \
#     quay.io/condaforge/linux-anvil-cos7-x86_64 \
#     ./build-package.sh --variant=linux_64_python3.10_cpython
#
# To use a local source tree, also mount it and pass --source-dir:
#   docker run --rm -it \
#     -v $(pwd):/home/conda \
#     -v /path/to/Chaste:/chaste:ro \
#     -e HOST_USER_ID="$(id -u)" \
#     quay.io/condaforge/linux-anvil-cos7-x86_64 \
#     ./build-package.sh --variant=linux_64_python3.10_cpython --source-dir=/chaste

# Parse args
variant=
branch=develop
parallel=
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
  --parallel=*)
    parallel=$(expr "x$option" : "x--parallel=\(.*\)")
    ;;
  *)
    echo "Unknown option: $option" 1>&2
    exit 1
    ;;
  esac
done

if [ -z "${variant}" ]; then
  echo "Error: --variant is required" 1>&2
  echo "Usage: $(basename "$0") --variant=<name> [--branch=<branch>] [--parallel=<n>] [--source-dir=<path>]" 1>&2
  exit 1
fi

set -x

# Configure environment
export FEEDSTOCK_ROOT="$(pwd)"
export RECIPE_ROOT="${FEEDSTOCK_ROOT}/recipe"
export CONFIG_FILE="${FEEDSTOCK_ROOT}/variants/${variant}.yaml"
export CONDA_BLD_PATH="${FEEDSTOCK_ROOT}/build_artifacts"
export CPU_COUNT="${parallel:-$(nproc)}"
export PYTHONUNBUFFERED=1

mkdir -p "${CONDA_BLD_PATH}"

# Configure conda
cat >~/.condarc <<CONDARC
conda-build:
  root-dir: ${CONDA_BLD_PATH}
pkgs_dirs:
  - ${CONDA_BLD_PATH}/pkg_cache
  - /opt/conda/pkgs
channels:
  - conda-forge
  - pychaste
channel_priority: strict
solver: libmamba
CONDARC

# Install conda-build
mamba install --update-specs --yes --quiet conda-build

# Show configuration
cat "${CONFIG_FILE}"
conda info
conda list --show-channel-urls

# Install system dependencies
/usr/bin/sudo -n yum install -y libXt-devel mesa-libGLU-devel patch

# Get source code
if [ -n "${source_dir}" ]; then
  export CHASTE_SOURCE_DIR="${source_dir}"
else
  git clone --recursive --branch "${branch}" --depth 1 --tags https://github.com/Chaste/Chaste.git /tmp/Chaste
  export CHASTE_SOURCE_DIR=/tmp/Chaste
fi

# Determine package version from git tag or commit hash
export CHASTE_VERSION=$(bash "${FEEDSTOCK_ROOT}/package-version.sh" "${CHASTE_SOURCE_DIR}")

conda build "${RECIPE_ROOT}" --variant-config-files "${CONFIG_FILE}"
