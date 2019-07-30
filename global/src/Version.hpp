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

#ifndef VERSION_HPP_
#define VERSION_HPP_

#include <string>
#include <map>

/**
 * A class with static methods providing various information about this build of Chaste.
 */
class ChasteBuildInfo
{
public:
    /**
     * @return The licence notice for Chaste.
     */
    static std::string GetLicenceText();

    /**
     * @return The path to the root directory of the Chaste source tree.
     */
    static const char* GetRootDir();

    /**
     * @return A string representation of the current Chaste version. This combines the
     * information from GetMajorReleaseNumber, GetMinorReleaseNumber, and GetRevisionNumber.
     */
    static std::string GetVersionString();

    /**
     * @return The major number of the "current" Chaste release.
     * If this is a development build, this will be the number of the last release.
     *
     * @note This must be set manually by modifying Version.cpp.in.
     */
    static unsigned GetMajorReleaseNumber();

    /**
     * @return The minor number of the "current" Chaste release.
     * If this is a development build, this will be the number of the last release.
     *
     * @note This must be set manually by modifying Version.cpp.in.
     */
    static unsigned GetMinorReleaseNumber();

    /**
     * @return  Get the subversion revision number of the Chaste source tree.
     *
     * If the file ReleaseVersion.txt exists in the directory given by GetRootDir, then
     * we assume this is not a working copy, and read the version information from there.
     *
     * Otherwise, we assume this is a checked-out tree, and call svnversion
     * during the build.  If it returns a range of versions, the upper end of this range
     * is used.  Whether the working copy is modified is ignored by this method; use
     * IsWorkingCopyModified to test that.
     */
    static unsigned long long GetRevisionNumber();

    /**
     * @return  If this Chaste was built from a subversion working copy, then return whether there
     * were local modifications.  If it's not a working copy, return false.
     */
    static bool IsWorkingCopyModified();

    /**
     * @return The date and time at which Chaste was built.
     */
    static const char* GetBuildTime();

    /**
     * @return Get the current date and time, in the same format as GetBuildTime.
     * The returned 'string' is statically allocated, so you don't need to free the memory.
     * However, if you call this method again, the contents will be overwritten.
     */
    static const char* GetCurrentTime();

    /**
     * @return The output of "uname -a" on the machine that built Chaste.
     */
    static const char* GetBuilderUnameInfo();

    /**
     * @return  Information about this build of Chaste: the build type used, whether libraries
     * were used, and if so what kind.
     */
    static const char* GetBuildInformation();

    /**
     * @return The compiler type used to build (must be either 'intel' or 'gcc').
     */
    static const char* GetCompilerType();

    /**
     * @return  The compiler version number.
     */
    static const char* GetCompilerVersion();

    /**
     * @return  The compiler flags.
     */
    static const char* GetCompilerFlags();

    /**
     * @return  The XSD binary version number.
     */
    static const char* GetXsdVersion();

    /**
     * @return  The version numbers (i.e. revisions) of any checked-out projects.
     */
    static const std::map<std::string, std::string>& rGetProjectVersions();

    /**
     * @return  Whether any checked-out projects have uncommitted revisions.
     */
    static const std::map<std::string, std::string>& rGetIfProjectsModified();

    /**
     * @return  A single-line string representation of the provenance information to be attached
     * to any files we generate.  This includes the version of the Chaste code used, how and
     * when it was built, and the current date and time.
     */
    static std::string GetProvenanceString();
};

#endif // VERSION_HPP_
