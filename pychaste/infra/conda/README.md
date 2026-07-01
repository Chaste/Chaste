# PyChaste conda packages

The conda [pychaste](https://anaconda.org/pychaste) channel hosts the `chaste` and `xsd` conda packages. All other dependencies are installed from the [conda-forge](https://anaconda.org/conda-forge) channel.

## Building packages via GitHub Actions (recommended)

The easiest way to build a new conda package is to use the [PyChaste build conda package](../../../.github/workflows/pychaste-build-conda-package.yml) workflow. This can be triggered from the GitHub Actions tab with the following inputs:

- **python-version**: The Python version to build for (3.10–3.14).
- **chaste-branch**: The Chaste branch or tag to build (default: `develop`).
- **upload**: Whether to upload the built package to the `pychaste` Anaconda channel (default: `false`).

The workflow builds the package, runs a smoke test to verify that imports work, and uploads the built package as a GitHub Actions artifact. If `upload` is checked, it also pushes the package to the Anaconda channel.

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

## Building macOS (osx-64) packages locally

macOS packages cannot be built in the Linux container; they must be built
natively on a macOS host with `build-package-osx.sh`, which runs `rattler-build`
directly (no Docker). It downloads the macOS `rattler-build` binary and, because
conda-forge pins an older macOS SDK than recent Xcode ships, downloads and caches
the matching 10.13 SDK (into `.osx-sdk/`) to use as `CONDA_BUILD_SYSROOT`.

```bash
cd /path/to/Chaste/pychaste/infra/conda

./build-package-osx.sh \
  --variant=osx_64_python3.12_cpython \
  --source-dir=/path/to/Chaste \
  --cpu-count=<num-cpus>
```

The arguments match `build-package-linux.sh` (`--variant`, `--branch`, `--source-dir`,
`--cpu-count`). The built package is placed in `./build_artifacts/osx-64`.

## Directory structure

```
├── build-package-linux.sh  Docker-based local build script (linux-64)
├── build-package-osx.sh    Native local build script (osx-64, run on macOS)
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
    ├── osx_64_python3.10_cpython.yaml
    ├── osx_64_python3.11_cpython.yaml
    ├── osx_64_python3.12_cpython.yaml
    ├── osx_64_python3.13_cpython.yaml
    └── osx_64_python3.14_cpython.yaml
```

- `build-package-linux.sh`: Sets up the build environment and runs `rattler-build` to create the linux-64 package inside a Docker container.
- `build-package-osx.sh`: Builds the osx-64 package natively on a macOS host (no Docker); fetches `rattler-build` and the macOS 10.13 SDK.
- `package-version.sh`: Prints the package version derived from the source tree's git tags and commit history (used by the `build-package-*` scripts to set `CHASTE_VERSION`).
- `recipe/build.sh`: Used by `rattler-build` to compile the source. See the [rattler-build script documentation](https://rattler.build/latest/build_script/).
- `recipe/recipe.yaml`: Package metadata used by `rattler-build`. The source path is read from the `CHASTE_SOURCE_DIR` environment variable (defaults to `/tmp/Chaste` for Docker builds). See the [rattler-build recipe documentation](https://rattler.build/latest/reference/recipe_file/).
- `variants/`: Per-Python-version dependency pinning added on top of `recipe.yaml`. See the [rattler-build variant documentation](https://rattler.build/latest/variants/).
- `envs/`: Conda environment files used by the [PyChaste conda tests](../../../.github/workflows/pychaste-conda-tests.yml) CI workflow.
