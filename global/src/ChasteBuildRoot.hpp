/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CHASTEBUILDROOT_HPP_
#define CHASTEBUILDROOT_HPP_

/**
 * @file
 * A collection of functions providing information about the
 * filesystem layout of this Chaste source tree.
 */

#include <string>

/**
 * @return the path to the root directory of the Chaste build tree.
 * Will always give you the absolute path with a trailing slash.
 */
const char* ChasteBuildRootDir();

/**
 * @return the path to the root directory of the Chaste source tree.
 * Will always give you the absolute path with a trailing slash.
 */
const char* ChasteSourceRootDir();

/**
 * @return the folder in which compiled files are placed for the given
 * Chaste component.
 * Will always give you the absolute path with a trailing slash.
 *
 * @param rComponent  e.g. global, heart, pde, ...
 */
std::string ChasteComponentBuildDir(const std::string& rComponent);

/**
 * @return the name of the folder within the 'build' dir of a component
 * that contains the compiled files.
 */
std::string ChasteBuildDirName();

/**
 * @return the build type string used in building Chaste.
 */
std::string ChasteBuildType();

#endif /*CHASTEBUILDROOT_HPP_*/
