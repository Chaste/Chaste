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
  ./build-package.sh --variant=<package-variant> --branch=<chaste-branch> --parallel=<num-cpus>
```

The `--variant` argument specifies a build variant name. A list of variant files can be found in the `variants` subdirectory. The variant name is the file name without the `.yaml` extension.

The optional `--branch` argument specifies the Chaste branch or tag to build; the default is `develop`.

The optional `--parallel` argument specifies the maximum number of parallel build processes.

After the build completes, verify that the package was created:

```bash
ls ./build_artifacts/linux-64
```

There should be a `chaste-<version>-<hash>.tar.bz2` file in the directory.

To upload the package to the [pychaste](https://anaconda.org/pychaste) Anaconda channel:

```bash
conda install anaconda-client
anaconda login
anaconda upload -u pychaste ./build_artifacts/linux-64/chaste-<version>-<hash>.tar.bz2
```

## Directory structure

```
├── build-package.sh        Docker-based local build script
├── envs/                   Conda environment specs for each Python version (used by CI tests)
│   ├── env_python3.10.yaml
│   ├── env_python3.11.yaml
│   ├── env_python3.12.yaml
│   ├── env_python3.13.yaml
│   └── env_python3.14.yaml
├── recipe/
│   ├── build.sh            Build script executed by conda-build (configure + make + pip install)
│   └── meta.yaml           Package metadata (dependencies, version, source path)
└── variants/               Per-Python-version build configuration for conda-build
    ├── linux_64_python3.10_cpython.yaml
    ├── linux_64_python3.11_cpython.yaml
    ├── linux_64_python3.12_cpython.yaml
    ├── linux_64_python3.13_cpython.yaml
    └── linux_64_python3.14_cpython.yaml
```

- `build-package.sh`: Sets up the build environment and runs `conda build` to create the package.
- `recipe/build.sh`: Used by `conda build` to compile the source. See the [conda-build script documentation](https://docs.conda.io/projects/conda-build/en/latest/resources/build-scripts.html).
- `recipe/meta.yaml`: Package metadata used by `conda build`. The source path is read from the `CHASTE_SOURCE_DIR` environment variable (defaults to `/tmp/Chaste` for Docker builds). See the [conda-build metadata documentation](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html).
- `variants/`: Per-Python-version dependency pinning added on top of `meta.yaml`. See the [conda-build variant documentation](https://docs.conda.io/projects/conda-build/en/stable/resources/variants.html).
- `envs/`: Conda environment files used by the [PyChaste conda tests](../../../.github/workflows/pychaste-conda-tests.yml) CI workflow.
