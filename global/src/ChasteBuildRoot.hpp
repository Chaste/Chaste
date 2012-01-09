/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
 * Get the path to the root directory of the Chaste source tree.
 * Will always give you the absolute path with a trailing slash.
 */
const char* ChasteBuildRootDir();

/**
 * Get the folder in which compiled files are placed for the given
 * Chaste component.
 * Will always give you the absolute path with a trailing slash.
 *
 * @param rComponent  e.g. global, heart, pde, ...
 */
std::string ChasteComponentBuildDir(const std::string& rComponent);

/**
 * Get the name of the folder within the 'build' dir of a component
 * that contains the compiled files.
 */
std::string ChasteBuildDirName();

/**
 * Get the build type string used in building Chaste.
 */
std::string ChasteBuildType();

#endif /*CHASTEBUILDROOT_HPP_*/
