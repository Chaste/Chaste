/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef VERSION_HPP_
#define VERSION_HPP_

#include <string>

/**
 * A class with static methods providing various information about this build of Chaste.
 */
class ChasteBuildInfo
{
public:

    /**
     * @return The path to the root directory of the Chaste source tree.
     */
    static const char* GetRootDir();

    /**
     * Get a string representation of the current Chaste version. This combines the
     * information from GetMajorReleaseNumber, GetMinorReleaseNumber, and GetRevisionNumber.
     */
    static std::string GetVersionString();

    /**
     * Get the major number of the "current" Chaste release.
     * If this is a development build, this will be the number of the last release.
     *
     * @note This must be set manually by modifying Version.cpp.in.
     */
    static unsigned GetMajorReleaseNumber();

    /**
     * Get the minor number of the "current" Chaste release.
     * If this is a development build, this will be the number of the last release.
     *
     * @note This must be set manually by modifying Version.cpp.in.
     */
    static unsigned GetMinorReleaseNumber();

    /**
     * Get the subversion revision number of the Chaste source tree.
     *
     * If the file ReleaseVersion.txt exists in the directory given by GetRootDir, then
     * we assume this is not a working copy, and read the version information from there.
     *
     * Otherwise, we assume this is a checked-out tree, and call svnversion
     * during the build.  If it returns a range of versions, the upper end of this range
     * is used.  Whether the working copy is modified is ignored by this method; use
     * IsWorkingCopyModified to test that.
     */
    static unsigned GetRevisionNumber();

    /**
     * If this Chaste was built from a subversion working copy, then return whether there
     * were local modifications.  If it's not a working copy, return false.
     */
    static bool IsWorkingCopyModified();

    /**
     * @return The date and time at which Chaste was built.
     */
    static const char* GetBuildTime();

    /**
     * Get the current date and time, in the same format as GetBuildTime.
     * The returned 'string' is statically allocated, so you don't need to free the memory.
     * However, if you call this method again, the contents will be overwritten.
     */
    static const char* GetCurrentTime();

    /**
     * @return The output of "uname -a" on the machine that built Chaste.
     */
    static const char* GetBuilderUnameInfo();

    /**
     * Get information about this build of Chaste: the build type used, whether libraries
     * were used, and if so what kind.
     */
    static const char* GetBuildInformation();

    /**
     * Get the compiler type used to build (must be either 'intel' or 'gcc').
     */
     static const char* GetCompilerType();

    /**
     * Get the compiler version number.
     */
     static const char* GetCompilerVersion();

    /**
     * Get the compiler flags.
     */
     static const char* GetCompilerFlags();

     /**
      * Get the XSD binary version number.
      */
     static const char* GetXsdVersion();

    /**
     * Get a single-line string representation of the provenance information to be attached
     * to any files we generate.  This includes the version of the Chaste code used, how and
     * when it was built, and the current date and time.
     */
    static std::string GetProvenanceString();
};

#endif // VERSION_HPP_
