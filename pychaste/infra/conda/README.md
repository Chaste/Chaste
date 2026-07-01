# PyChaste conda packages

The conda [pychaste](https://anaconda.org/pychaste) channel hosts the `chaste` and `xsd` conda packages. All other dependencies are installed from the [conda-forge](https://anaconda.org/conda-forge) channel.

## Building packages via GitHub Actions (recommended)

The easiest way to build new conda packages is to run the [PyChaste build conda package](../../../.github/workflows/pychaste-build-conda-package.yml) workflow manually (`workflow_dispatch`) from the GitHub Actions tab, with the following inputs:

- **python-version**: The Python version to build for (3.10–3.14).
- **chaste-ref**: The Chaste branch, tag, or commit ref to build (default: `develop`).
- **upload**: Whether to upload the built packages to the `pychaste` Anaconda channel (default: `false`).

The workflow builds a matrix of platforms, each natively on a matching-architecture runner — `linux-64` (`ubuntu-24.04`), `linux-aarch64` (`ubuntu-24.04-arm`), and `osx-arm64` (`macos-15`) — running a smoke test on each and uploading each package as a GitHub Actions artifact. If `upload` is checked, it also pushes them to the Anaconda channel. Intel macOS (`osx-64`) is not built in CI (it needs a paid larger runner); build it locally instead (see below).

## Building packages locally with Docker

This directory also contains scripts for building the `chaste` conda package locally inside a Docker container.

The following instructions assume that [mamba](https://mamba.readthedocs.io) and [docker](https://docs.docker.com/get-docker/) have already been installed.

Launch a docker container to build the package:

```bash
cd /path/to/Chaste/pychaste/infra/conda

docker run -it --rm \
  -v $(pwd):/home/conda \
  -e HOST_USER_ID="$(id -u)" \
  quay.io/condaforge/linux-anvil-cos7-x86_64 \
  ./build-package-linux.sh --variant=<package-variant> --branch=<chaste-branch> --cpu-count=<num-cpus>
```

The `--variant` argument specifies a build variant name. A list of variant files can be found in the `variants` subdirectory. The variant name is the file name without the `.yaml` extension.

The optional `--branch` argument specifies the Chaste branch or tag to build; the default is `develop`.

The optional `--cpu-count` argument specifies the maximum number of parallel build processes.

To build from a local Chaste checkout instead of cloning from GitHub, mount it into the container and pass `--source-dir` (see the comments at the top of `build-package-linux.sh`).

To build for `linux-aarch64`, run the aarch64 container image (`quay.io/condaforge/linux-anvil-cos7-aarch64`, on an arm64 host) and pass `--target=linux-aarch64` with the matching `linux_aarch64_*` variant.

After the build completes, verify that the package was created:

```bash
ls ./build_artifacts/linux-64
```

There should be a `chaste-<version>-<hash>.conda` file in the directory.

To upload the package to the [pychaste](https://anaconda.org/pychaste) Anaconda channel:

```bash
conda install anaconda-client
anaconda login
anaconda upload -u pychaste ./build_artifacts/linux-64/chaste-<version>-<hash>.conda
```

## Building macOS packages locally

macOS packages cannot be built in the Linux container; they must be built
natively on a macOS host with `build-package-osx.sh`, which runs `rattler-build`
directly (no Docker). It downloads the macOS `rattler-build` binary and, because
conda-forge pins an older macOS SDK than recent Xcode ships, downloads and caches
the matching SDK (into `.osx-sdk/`) to use as `CONDA_BUILD_SYSROOT`.

`--target` must match the host architecture: `osx-64` on an Intel Mac (SDK 10.15,
the minimum providing `aligned_alloc`) or `osx-arm64` on an Apple Silicon Mac
(SDK 11.0, the floor for Apple Silicon).

```bash
cd /path/to/Chaste/pychaste/infra/conda

# On an Intel Mac:
./build-package-osx.sh \
  --variant=osx_64_python3.12_cpython \
  --target=osx-64 \
  --source-dir=/path/to/Chaste \
  --cpu-count=<num-cpus>

# On an Apple Silicon Mac:
./build-package-osx.sh \
  --variant=osx_arm64_python3.12_cpython \
  --target=osx-arm64 \
  --source-dir=/path/to/Chaste \
  --cpu-count=<num-cpus>
```

The other arguments match `build-package-linux.sh` (`--variant`, `--branch`,
`--source-dir`, `--cpu-count`). The built package is placed in
`./build_artifacts/<target>`.

## Directory structure

```
├── build-package-linux.sh  Docker-based local build script (linux-64/linux-aarch64)
├── build-package-osx.sh    Native local build script (osx-64/osx-arm64, run on macOS)
├── package-version.sh      Derives the package version from the source tree's git tags/history
├── envs/                   Conda environment specs for each Python version (used by CI tests)
│   ├── env_python3.10.yaml
│   ├── env_python3.11.yaml
│   ├── env_python3.12.yaml
│   ├── env_python3.13.yaml
│   └── env_python3.14.yaml
├── recipe/
│   ├── build.sh            Build script executed by rattler-build (configure + make + pip install)
│   └── recipe.yaml         Package metadata (dependencies, version, source path)
└── variants/               Per-platform/per-Python build configuration for rattler-build
    ├── linux_64_python3.10_cpython.yaml
    ├── linux_64_python3.11_cpython.yaml
    ├── linux_64_python3.12_cpython.yaml
    ├── linux_64_python3.13_cpython.yaml
    ├── linux_64_python3.14_cpython.yaml
    ├── linux_aarch64_python3.10_cpython.yaml
    ├── linux_aarch64_python3.11_cpython.yaml
    ├── linux_aarch64_python3.12_cpython.yaml
    ├── linux_aarch64_python3.13_cpython.yaml
    ├── linux_aarch64_python3.14_cpython.yaml
    ├── osx_64_python3.10_cpython.yaml
    ├── osx_64_python3.11_cpython.yaml
    ├── osx_64_python3.12_cpython.yaml
    ├── osx_64_python3.13_cpython.yaml
    ├── osx_64_python3.14_cpython.yaml
    ├── osx_arm64_python3.10_cpython.yaml
    ├── osx_arm64_python3.11_cpython.yaml
    ├── osx_arm64_python3.12_cpython.yaml
    ├── osx_arm64_python3.13_cpython.yaml
    └── osx_arm64_python3.14_cpython.yaml
```

- `build-package-linux.sh`: Sets up the build environment and runs `rattler-build` to create a Linux package (`linux-64` or `linux-aarch64`, selected with `--target`) inside a Docker container.
- `build-package-osx.sh`: Builds a macOS package (`osx-64` or `osx-arm64`, selected with `--target`) natively on a macOS host (no Docker); fetches `rattler-build` and the matching macOS SDK.
- `package-version.sh`: Prints the package version derived from the source tree's git tags and commit history (used by the `build-package-*` scripts to set `CHASTE_VERSION`).
- `recipe/build.sh`: Used by `rattler-build` to compile the source. See the [rattler-build script documentation](https://rattler.build/latest/build_script/).
- `recipe/recipe.yaml`: Package metadata used by `rattler-build`. The source path is read from the `CHASTE_SOURCE_DIR` environment variable (defaults to `/tmp/Chaste` for Docker builds). See the [rattler-build recipe documentation](https://rattler.build/latest/reference/recipe_file/).
- `variants/`: Per-platform, per-Python dependency pinning added on top of `recipe.yaml`. See the [rattler-build variant documentation](https://rattler.build/latest/variants/).
- `envs/`: Conda environment files used by the [PyChaste conda tests](../../../.github/workflows/pychaste-conda-tests.yml) CI workflow.
