VS Code Development Container
=============================

Using this development environment requires:
* [Docker](https://www.docker.com/get-started) we recommend installing Docker Desktop (unless you want to use a remote server)
* [VS Code](https://code.visualstudio.com/)
* [Remote development extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack)

Cloning the Chaste repository and opening it in VS Code should result in a prompt asking if you want to reopen the project in a container. Clicking `Yes` will start the build process. This will take a few minutes for the first time but will be faster on subsequent rebuilds.

You now have an isolated development environment (which will not conflict with packages installed elsewhere on your system) with all the dependencies needed for Chaste already installed, laidout as described in the main repository [README](https://github.com/Chaste/chaste-docker/blob/master/README.md).

The container environment can be customised in many ways, such as with [dotfiles](https://code.visualstudio.com/docs/remote/containers#_personalizing-with-dotfile-repositories) if hosted in a public repository. Further documentation for development in containers can be found here: https://code.visualstudio.com/docs/devcontainers/containers.
