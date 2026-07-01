#!/bin/bash -e

# Build a PyChaste conda package for linux-64 with rattler-build inside a
# conda-forge build container. For macOS (osx-64), use build-package-osx.sh.
#
# Example usage (clone from GitHub):
#   ./build-package-linux.sh --variant=linux_64_python3.10_cpython --branch=2026.1 --cpu-count=4
#
# Example usage (local source):
#   ./build-package-linux.sh --variant=linux_64_python3.10_cpython --source-dir=/chaste
#
# Intended for use in a build container e.g.:
#   docker run --rm -it \
#     -v $(pwd):/home/conda \
#     -e HOST_USER_ID="$(id -u)" \
#     quay.io/condaforge/linux-anvil-cos7-x86_64 \
#     ./build-package-linux.sh --variant=linux_64_python3.10_cpython
#
# To use a local source tree, also mount it and pass --source-dir:
#   docker run --rm -it \
#     -v $(pwd):/home/conda \
#     -v /path/to/Chaste:/chaste:ro \
#     -e HOST_USER_ID="$(id -u)" \
#     quay.io/condaforge/linux-anvil-cos7-x86_64 \
#     ./build-package-linux.sh --variant=linux_64_python3.10_cpython --source-dir=/chaste

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

set -x

# The container entrypoint runs useradd with the mounted source tree as the home
# directory, which drops bash profile files into it. Move them and redirect HOME
# to a temp path so subsequent writes also stay off the source tree.
OLD_HOME="${HOME}"
export HOME=/tmp/conda_home
mkdir -p "${HOME}"
mv -f "${OLD_HOME}/.bash_logout" "${OLD_HOME}/.bash_profile" "${OLD_HOME}/.bashrc" \
      "${OLD_HOME}/.profile" "${OLD_HOME}/.condarc" "${HOME}/" 2>/dev/null || true

# Configure environment
export FEEDSTOCK_ROOT="$(pwd)"
export RECIPE_ROOT="${FEEDSTOCK_ROOT}/recipe"
export CONFIG_FILE="${FEEDSTOCK_ROOT}/variants/${variant}.yaml"
# CPU_COUNT bounds parallelism for both rattler-build and the C++ compile
# (make -j). Keep it low on memory-constrained machines: the pychaste wrapper
# translation units and their LTO link consume roughly 3 GB each, so building
# with all cores can exhaust RAM. It is forwarded into the sanitized build
# script environment via recipe.yaml.
export CPU_COUNT="${cpu_count:-$(nproc)}"
export PYTHONUNBUFFERED=1

# Keep all build artefacts on the container's own filesystem to avoid
# case-sensitivity issues on macOS hosts (e.g. ncurses terminfo entries that
# differ only by case are collapsed on APFS, breaking package verification).
export OUTPUT_DIR=/tmp/rattler_output
mkdir -p "${OUTPUT_DIR}"

# Install rattler-build as a self-contained static binary (no conda solve).
curl -fsSL \
  https://github.com/prefix-dev/rattler-build/releases/latest/download/rattler-build-x86_64-unknown-linux-musl \
  -o /tmp/rattler-build
chmod +x /tmp/rattler-build
export PATH="/tmp:${PATH}"

# Show configuration
cat "${CONFIG_FILE}"
rattler-build --version

# Install system dependencies
/usr/bin/sudo -n yum install -y libXt-devel mesa-libGLU-devel patch

# Get source code as a shallow (depth 1) clone so the full git history is not
# dragged into the build.
if [ -n "${source_dir}" ]; then
  # The mounted source tree is owned by the host user; mark it safe so git works.
  git config --global --add safe.directory "${source_dir}"
  # Determine the version from the full local repo (it has the tags and history
  # the shallow clone lacks) before making the shallow clone. file:// is required
  # for git to honour --depth when cloning a local path.
  export CHASTE_VERSION=$(bash "${FEEDSTOCK_ROOT}/package-version.sh" "${source_dir}")
  git clone --recurse-submodules --depth 1 "file://${source_dir}" /tmp/Chaste
else
  git clone --recursive --branch "${branch}" --depth 1 --tags https://github.com/Chaste/Chaste.git /tmp/Chaste
  export CHASTE_VERSION=$(bash "${FEEDSTOCK_ROOT}/package-version.sh" /tmp/Chaste)
fi
export CHASTE_SOURCE_DIR=/tmp/Chaste

rattler-build build \
  --recipe "${RECIPE_ROOT}/recipe.yaml" \
  --variant-config "${CONFIG_FILE}" \
  --target-platform linux-64 \
  --output-dir "${OUTPUT_DIR}" \
  --channel conda-forge \
  --channel pychaste

# Copy the built package to the mounted host directory
mkdir -p "${FEEDSTOCK_ROOT}/build_artifacts/linux-64"
cp "${OUTPUT_DIR}"/linux-64/chaste-*.conda "${FEEDSTOCK_ROOT}/build_artifacts/linux-64/"
