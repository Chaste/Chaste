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

#ifndef TESTSOURCEREVISION_HPP_
#define TESTSOURCEREVISION_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include "SourceRevision.hpp"
#include "FakePetscSetup.hpp"

class TestSourceRevision : public CxxTest::TestSuite
{
public:

    void TestUnknownRevision()
    {
        // The default is what ChasteBuildInfo reports when a commit could not be determined
        const SourceRevision revision;
        TS_ASSERT_EQUALS(revision.rGetCommit(), "unknown");
        TS_ASSERT_EQUALS(revision.GetShortCommit(), "unknown");
        TS_ASSERT(!revision.IsKnown());
        TS_ASSERT(!revision.IsModified());
    }

    void TestKnownRevision()
    {
        const std::string commit = "0123456789abcdef0123456789abcdef01234567";
        const SourceRevision revision(commit, true);
        TS_ASSERT_EQUALS(revision.rGetCommit(), commit);
        TS_ASSERT(revision.IsKnown());
        TS_ASSERT(revision.IsModified());

        // Shortened to a fixed length for display
        TS_ASSERT_EQUALS(SourceRevision::SHORT_COMMIT_LENGTH, 12u);
        TS_ASSERT_EQUALS(revision.GetShortCommit(), "0123456789ab");

        // A SHA-256 repository's longer hashes shorten to the same length
        const SourceRevision sha256(commit + "89abcdef0123456789abcdef", false);
        TS_ASSERT_EQUALS(sha256.GetShortCommit(), "0123456789ab");
        TS_ASSERT(!sha256.IsModified());

        // A commit already shorter than that is shown in full
        TS_ASSERT_EQUALS(SourceRevision("abc123").GetShortCommit(), "abc123");
    }
};

#endif /* TESTSOURCEREVISION_HPP_ */
