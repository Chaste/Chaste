/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef TESTCITATIONS_HPP_
#define TESTCITATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "Citations.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "FileComparison.hpp"
#include "OutputFileHandler.hpp"
#include "PetscAndChasteCitations.hpp"
#include "PetscException.hpp"
#include "PetscSetupUtils.hpp"

class TestCitations : public CxxTest::TestSuite
{
public:
    void TestChasteCitation() throw(Exception)
    {
        /*
         * Get location of output file to pass through as -citation argument.
         * At this point we haven't called PetscSetup, so the number of threads
         * running at this point is implementation-dependent, and it's not really
         * safe to do anything collective! Fortunately FileFinder isn't.
         */
        FileFinder output_file("TestCitations/citations.txt", RelativeTo::ChasteTestOutput);
        // Turn on citations with argument
        CommandLineArgumentsMocker mocker("-citations " + output_file.GetAbsolutePath());

        // First test Chaste implementation on its own, as PETSc isn't yet initialised
        {
            Citations::Register(PetscCitation1, &PetscCite1);
            Citations::Register(PetscCitation2, &PetscCite2);
            Citations::Register(ChasteCitation, &ChasteCite);

            Citations::Print();

            FileFinder reference_citations("global/test/data/citations.txt", RelativeTo::ChasteSourceRoot);
            FileComparison check_files(output_file, reference_citations, false); // Not collective
            check_files.CompareFiles();

            TS_ASSERT_EQUALS(Citations::mUseChasteImplementation, true);
        }

        // Reset to default
        Citations::mUseChasteImplementation = false;

        // Now test as part of PETSc, PETSc implementation will be used if PETSc is new enough...
        // if not this is something of a duplicate of the above test!
        {

            PetscSetupUtils::CommonSetup(); // This automatically includes some citations
            /*
			 * Make empty directory now that Petsc is set up, because this must be done
			 * collectively.
			 */
            OutputFileHandler handler("TestCitations");
            PetscSetupUtils::CommonFinalize(); // This prints the citations to disk

// Check PETSc version - this is just because they reformatted their BibTex between versions, no change to function!
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7) //PETSc 3.7
            FileFinder reference_citations("global/test/data/citations-2016.txt", RelativeTo::ChasteSourceRoot);
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 6) //PETSc 3.6
            FileFinder reference_citations("global/test/data/citations-2015.txt", RelativeTo::ChasteSourceRoot);
#else
            // Use PETSc 3.4 older version (or mocked up version that matches 3.4)
            FileFinder reference_citations("global/test/data/citations.txt", RelativeTo::ChasteSourceRoot);
#endif
            FileComparison check_files(output_file, reference_citations, false); // Not collective
            check_files.CompareFiles();
        }
    }
};

#endif // TESTCITATIONS_HPP_
