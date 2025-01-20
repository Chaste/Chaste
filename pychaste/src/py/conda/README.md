# PyChaste conda packages

The conda [pychaste](https://anaconda.org/pychaste) channel hosts the `chaste` and `xsd` conda packages. All other dependencies are installed from the [conda-forge](https://anaconda.org/conda-forge) channel.

## How to build the conda packages

This directory contains scripts for building the `chaste` conda package and uploading it to the conda [pychaste](https://anaconda.org/pychaste) channel.

The following instructions assume that [mamba](https://mamba.readthedocs.io) and [docker](https://docs.docker.com/get-docker/) have already been installed.

Launch a docker container to build the package:

```bash
cd /path/to/Chaste/pychaste/src/py/conda

docker run -it --rm \
  -v $(pwd):/home/conda \
  -e HOST_USER_ID="$(id -u)" \
  quay.io/condaforge/linux-anvil-cos7-x86_64 \
  ./build-package.sh --variant=<package-variant> --branch=<chaste-branch> --parallel=<num-cpus>
```

The `--variant` argument specifies a build variant name. A list of variant files can be found in the `variants` subdirectory. The variant name is the name of the file without the `.yaml` extenstion.

The optional `--build` argument specifies the Chaste branch/tag to build from; the default is `develop`.

The optional `--parallel` argument specifies the maximum number of parallel build processes.

After the build is complete, verify that the package has been created:

```bash
ls ./build_artifacts/linux-64
```

There should be a `<chaste>-<version>-<hash>.tar.bz2` file in the directory.

Test the newly built package by installing it in a conda environment:

```bash
conda create -n py38-env python=3.8
conda activate py38-env
conda install -c ./build_artifacts <chaste>-<version>-<hash>.tar.bz2
```

Run tests on the package:

```bash
python -m unittest discover /path/to/Chaste/pychaste/test
```

To upload the package to the conda [pychaste](https://anaconda.org/pychaste) channel, first install `anaconda-client`:

```bash
conda install anaconda-client
```

Login and upload the package:

```bash
anaconda login
anaconda upload -u 'pychaste' ./build_artifacts/<chaste>-<version>-<hash>.tar.bz2
```

## Directory structure

The conda scripts directory has this structure:

```
├── build-package.sh
├── recipe
│   ├── build.sh
│   └── meta.yaml
└── variants
    ├── <variant-01>.yaml
    ├── <variant-02>.yaml
    └── <variant-03>.yaml
```

- `build-package.sh`: This sets up the environment for the build and runs `conda mambabuild` to create the package.
- `recipe/build.sh`: This script is used by `conda mambabuild` to build the source code. Typically, this performs the "configure" and "make" steps. See the [conda-build script documentation](https://docs.conda.io/projects/conda-build/en/latest/resources/build-scripts.html) for more information.
- `recipe/meta.yaml`: This contains metadata used by `conda mambabuild` such as build dependencies. For more information, see the [conda-build metadata documentation](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html).
- `variants`: This contains additional metadata for different variants of the build, which are added on top of the metadata provided in `meta.yaml`. Typically, a variant file adds specific versions of dependencies required for that variant. For further information, see the [conda-build variant documentation](https://docs.conda.io/projects/conda-build/en/stable/resources/variants.html).
