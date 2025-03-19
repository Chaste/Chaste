[![DOI](https://joss.theoj.org/papers/10.21105/joss.01848/status.svg)](https://doi.org/10.21105/joss.01848)[![Build chaste/develop](https://github.com/Chaste/Chaste/actions/workflows/docker-develop-image.yml/badge.svg)](https://github.com/Chaste/Chaste/actions/workflows/docker-develop-image.yml)

# Welcome to Chaste

If you are new to Chaste please see [our Getting Started wiki page](https://chaste.github.io/docs/).

The files you have downloaded contain the source code for all Chaste functionality. 
Chaste makes use of a variety of external libraries and packages that need to be installed on your machine. 
The [Install Guide webpage](https://chaste.github.io/docs/installguides/) 
provides a comprehensive guide on how to install these external tools.

Chaste is distributed as open source software under the [3-clause BSD licence](https://opensource.org/licenses/BSD-3-Clause). 
For full details see the file [Copying.pdf](docs/licencing/Copying.pdf).
Chaste uses various third party libraries which have their own licences. 
For details of these licences and the impact they may have on your use of 
Chaste please see [Licences.html](docs/licencing/Licences.html).

Chaste includes a complete test suite covering all the source code. 
The easiest way to use existing source codes is to create a test file 
which can call upon any of the source files.  
The Chaste build system can build this file for you and handle 
all of the dependencies and library calls.

We suggest you use the projects directory in this manner to store your own 
source and test files if you do not wish to modify the Chaste source code. 
For more information, see the [User Projects webpage](https://chaste.github.io/docs/user-guides/user-projects/).

For more information please refer to the Chaste website at: https://chaste.github.io

Information on changes in this release can be found on the [release notes webpage](https://chaste.github.io/docs/release-notes/release-notes/).

Tutorial examples for this release are available at:
https://chaste.github.io/releases/2024.2/user-tutorials/

API documentation generated from the code by Doxygen is available at:
https://chaste.github.io/doxygen-releases/release_2024.2

Chaste welcomes contributions from the community.
For information on how to contribute to Chaste, and for support and bug reports, please see the file [CONTRIBUTING.md](docs/CONTRIBUTING.md).

A number of external libraries have been created that build on the Chaste trunk code. These include the following:
 * Microvessel Chaste (https://jmsgrogan.github.io/MicrovesselChaste/)
 * ChemChaste (https://github.com/OSS-Lab/ChemChaste)

Note that, while Chaste developers may have been contributed to the development of these external libraries, we are unable to offer any support in their maintenance, testing or usage. If you have any questions about one of these external libraries, please contact that library's lead developer directly.
