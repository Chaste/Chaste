/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTCHASTEBUILDINFO_HPP_
#define TESTCHASTEBUILDINFO_HPP_

#include <cxxtest/TestSuite.h>

#include "ExecutableSupport.hpp"
#include "Version.hpp"
#include <sstream>
#include "FakePetscSetup.hpp"

/**
 * Most of what ChasteBuildInfo reports depends on the build, so we mainly print out all the
 * information available, and check the parts that should hold for any build.
 */
class TestChasteBuildInfo : public CxxTest::TestSuite
{
private:
    /**
     * @return whether rTime ends in a UTC offset of the form +hhmm or -hhmm
     * @param rTime  a time string
     */
    bool EndsWithUtcOffset(const std::string& rTime)
    {
        if (rTime.size() < 5)
        {
            return false;
        }
        const std::string offset = rTime.substr(rTime.size() - 5);
        return (offset[0] == '+' || offset[0] == '-') && offset.find_first_not_of("0123456789", 1) == std::string::npos;
    }

    /**
     * @return whether rRevision is either unknown, or a full commit hash (SHA-1 or SHA-256) whose
     * short form is its first SourceRevision::SHORT_COMMIT_LENGTH characters
     * @param rRevision  a source revision
     */
    bool IsWellFormed(const SourceRevision& rRevision)
    {
        if (!rRevision.IsKnown())
        {
            return rRevision.GetShortCommit() == "unknown";
        }
        const std::string& r_commit = rRevision.rGetCommit();
        return (r_commit.size() == 40u || r_commit.size() == 64u)
               && r_commit.find_first_not_of("0123456789abcdef") == std::string::npos
               && rRevision.GetShortCommit() == r_commit.substr(0, SourceRevision::SHORT_COMMIT_LENGTH);
    }

public:
    void TestShowInfo()
    {
        std::string info;
        ExecutableSupport::GetBuildInfo(info);
        std::cout << info << std::flush;
    }

    /**
     * Checks that hold for any build. Kept to a single test, since this suite is mostly run just
     * to print the information above.
     */
    void TestChecks()
    {
        // ChasteBuildInfo reads its build-time provenance from plain files under provenance/ (see
        // global/src/Version.cpp.in), falling back on fixed defaults if they cannot be found, so
        // check a broken lookup isn't passing unnoticed behind those defaults. Both times end in
        // their offset from UTC, e.g. "+0100".
        const std::string build_time = ChasteBuildInfo::GetBuildTime();
        TS_ASSERT_DIFFERS(build_time, "unknown");
        TS_ASSERT(EndsWithUtcOffset(build_time));
        TS_ASSERT(EndsWithUtcOffset(ChasteBuildInfo::GetCurrentTime()));

        // The banner keeps the licence terms, without the LICENSE file's heading, the source
        // header line, or a trailing newline
        const std::string licence = ChasteBuildInfo::GetLicenceText();
        TS_ASSERT(licence.find("University of Oxford") != std::string::npos);
        TS_ASSERT(licence.find("BSD 3-Clause License.") == std::string::npos);
        TS_ASSERT(licence.find("This file is part of Chaste.") == std::string::npos);
        TS_ASSERT(!licence.empty() && licence.back() != '\n');

        // Revisions are a full commit hash, shown shortened to a fixed length, or unknown; built
        // from a git working copy, Chaste's own commit must be known
        const SourceRevision chaste_revision = ChasteBuildInfo::GetChasteRevision();
        TS_ASSERT(IsWellFormed(chaste_revision));
        if (FileFinder(".git", RelativeTo::ChasteSourceRoot).Exists())
        {
            TS_ASSERT(chaste_revision.IsKnown());
        }
        for (const auto& r_project : ChasteBuildInfo::GetProjectRevisions())
        {
            TS_ASSERT(IsWellFormed(r_project.second));
        }

        // The version string carries the short commit as build metadata, when it is known
        std::stringstream expected_version;
        expected_version << ChasteBuildInfo::GetMajorReleaseNumber() << "." << ChasteBuildInfo::GetMinorReleaseNumber();
        if (chaste_revision.IsKnown())
        {
            expected_version << "+" << chaste_revision.GetShortCommit();
        }
        TS_ASSERT_EQUALS(ChasteBuildInfo::GetVersionString(), expected_version.str());

        // Nothing inside Chaste tests these macros any more (VTK, CVODE and Xerces are required
        // dependencies), but downstream code may still guard on them, so they must remain defined.
        // See cmake/Modules/ChasteLegacyDefinitions.cmake.
#ifndef CHASTE_VTK
        TS_FAIL("CHASTE_VTK must remain defined, see cmake/Modules/ChasteLegacyDefinitions.cmake");
#endif
#ifndef CHASTE_CVODE
        TS_FAIL("CHASTE_CVODE must remain defined, see cmake/Modules/ChasteLegacyDefinitions.cmake");
#endif
#ifndef CHASTE_XERCES
        TS_FAIL("CHASTE_XERCES must remain defined, see cmake/Modules/ChasteLegacyDefinitions.cmake");
#endif
    }
};

#endif /* TESTCHASTEBUILDINFO_HPP_ */
