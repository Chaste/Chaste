PyChaste Notebook Dev Container
===============================

A lightweight environment for **running PyChaste notebooks** — not for building
Chaste from source. It uses the prebuilt [`chaste/pychaste`](https://hub.docker.com/r/chaste/pychaste)
image (~672 MB: a `micromamba` base with the `chaste` conda package and
JupyterLab installed), so it is far smaller than `chaste/develop` (~24 GB).

Usage
-----

**GitHub Codespaces:** create a Codespace on the repository and select the
**Chaste/PyChaste** configuration when prompted (or set it as the default dev
container). It boots ready to run notebooks in a couple of minutes.

**Local VS Code:** run *Dev Containers: Reopen in Container* and choose
**Chaste/PyChaste** (requires Docker + the Dev Containers extension; see the
[main devcontainer README](../README.md)).

Then open a notebook (for example under `pychaste/src/py/doc/tutorial/`), select
the `/opt/conda/bin/python` kernel, and run the cells.

What the configuration sets up
------------------------------

- **Image:** `chaste/pychaste`, running as the image's `mambauser`.
- **Workspace:** the repository is mounted at `/home/mambauser/src`.
- **Display:** `DISPLAY=:99` with `Xvfb` started on attach, so VTK
  visualisation renders without a physical display.
- **Extensions:** the VS Code Python and Jupyter extensions, with the
  interpreter pinned to `/opt/conda/bin/python`.

Notes
-----

- The image ships a `pychaste` conda build from the `pychaste` channel.
- This image does **not** include a C++ toolchain or the Chaste sources built;
  use the default `.devcontainer` (or `chaste/develop`) if you need C++.
