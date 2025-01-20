#!/bin/bash -e

# Example usage:
#   ./build-package.sh --variant=linux_64_python3.8_cpython --branch=2024.1 --parallel=4
#
# Intended for use in a build container e.g.:
#   docker run --rm -it quay.io/condaforge/linux-anvil-cos7-x86_64 /bin/bash

# Parse args
variant=
branch=
parallel=

for option; do
  case $option in
  --variant=*)
    variant=$(expr "x$option" : "x--variant=\(.*\)")
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

if [ -z "${variant}" ]; then usage; fi
if [ -z "${branch}" ]; then branch=develop; fi

set -x

# Configure environment
export FEEDSTOCK_ROOT="$(pwd)"
export RECIPE_ROOT="${FEEDSTOCK_ROOT}/recipe"
export CONFIG_FILE="${FEEDSTOCK_ROOT}/variants/${variant}.yaml"
export CONDA_BLD_PATH="${FEEDSTOCK_ROOT}/build_artifacts"

export CPU_COUNT="${parallel:-$(nproc)}"

export PYTHONUNBUFFERED=1

# Configure conda build path
mkdir -p "${CONDA_BLD_PATH}"

cat >~/.condarc <<CONDARC
conda-build:
  root-dir: ${CONDA_BLD_PATH}
pkgs_dirs:
  - ${CONDA_BLD_PATH}/pkg_cache
  - /opt/conda/pkgs
CONDARC

# Install conda build tools
mamba install --update-specs --yes --quiet --channel conda-forge \
  conda-build pip boa liblief=0.11.5 conda-forge-ci-setup=3
mamba update --update-specs --yes --quiet --channel conda-forge \
  conda-build pip boa liblief=0.11.5 conda-forge-ci-setup=3

# Configure conda channels
conda config --add channels conda-forge
conda config --add channels pychaste
conda config --env --set show_channel_urls true
conda config --env --set auto_update_conda false
conda config --env --set add_pip_as_python_dependency false
conda config --env --append aggressive_update_packages ca-certificates
conda config --env --remove-key aggressive_update_packages
conda config --env --append aggressive_update_packages ca-certificates
conda config --env --append aggressive_update_packages certifi
conda config --env --set channel_priority strict

# Configure conda activation
mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d"
cat >"${CONDA_PREFIX}"/etc/conda/activate.d/conda-forge-ci-setup-activate.sh <<CONDAACTIVATE
export CONDA_BLD_PATH='${CONDA_BLD_PATH}'
export CPU_COUNT='${CPU_COUNT:-}'
export PYTHONUNBUFFERED='1'
CONDAACTIVATE

# Show conda build configuration
cat "${CONFIG_FILE}"

conda info
conda config --env --show-sources
conda list --show-channel-urls

# Install system dependencies
/usr/bin/sudo -n yum install -y libXt-devel mesa-libGLU-devel patch

# Get source code
git clone --recursive --branch "${branch}" --depth 1 https://github.com/Chaste/Chaste.git /tmp/Chaste

# Build conda package
conda mambabuild "${RECIPE_ROOT}" --variant-config-files "${CONFIG_FILE}"
