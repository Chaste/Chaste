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

#ifndef TESTCITATIONS_HPP_
#define TESTCITATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "Citations.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "FileComparison.hpp"
#include "OutputFileHandler.hpp"
#if ((PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) || PETSC_VERSION_MAJOR > 3)
#include "PetscAndChasteCitations.hpp"
#endif
#include "PetscException.hpp"
#include "PetscSetupUtils.hpp"

// Important that this test DOESN'T include a PetscSetupAndFinalise.hpp !

class TestCitations : public CxxTest::TestSuite
{
public:
    void TestChasteCitation()
    {
        /*
         * Get location of output file to pass through as -citation argument.
         * At this point we haven't called PetscSetup, so the number of threads
         * running at this point is implementation-dependent, and it's not really
         * safe to do anything collective! Fortunately FileFinder isn't.
         */
        FileFinder output_petsc_file("TestCitations/petsc_citations.txt", RelativeTo::ChasteTestOutput);

        // First part of PETSc citation testing.
        {
            // Turn on citations with argument - N.B. we direct PETSc to a different file to Chaste implementation.
            // N.B. This is needed twice - PETSc implementation reads the first one,
            // Chaste implementation masquerading as PETSc reads the second!
            CommandLineArgumentsMocker mocker("-citations " + output_petsc_file.GetAbsolutePath());
            // Setup PETSc (with command line for citations pointing to PETSc output file)
            PetscSetupUtils::CommonSetup();
        }

        /*
         * Make empty directory, has to be done after PETSc is set up to be safe in parallel.
         */
        OutputFileHandler handler("TestCitations");

// If we need to check the Chaste version manually (PETSc is new)
// Otherwise on old PETSc the test below covers what is in this guard.
#if ((PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5) || PETSC_VERSION_MAJOR > 3)
        FileFinder output_chaste_file("TestCitations/chaste_citations.txt", RelativeTo::ChasteTestOutput);

        // Test Chaste implementation on its own (for PETSc <= 3.4, or no PETSc setup)
        Citations::mUseChasteImplementation = true;
        {
            std::cout << "Testing Chaste implementation" << std::endl;
            // Turn on citations with argument
            CommandLineArgumentsMocker mocker("-citations " + output_chaste_file.GetAbsolutePath());

            Citations::Register(PetscCitation1, &PetscCite1);
            Citations::Register(PetscCitation2, &PetscCite2);
            Citations::Register(ChasteCitation, &ChasteCite);

            Citations::Print(); // Writes the citations file.

            PetscTools::Barrier("Make sure the master process has finished writing the citations file.");

            FileFinder reference_citations("global/test/data/citations.txt", RelativeTo::ChasteSourceRoot);
            FileComparison check_files(output_chaste_file, reference_citations); // Collective call (default behaviour)
            check_files.CompareFiles();

            TS_ASSERT_EQUALS(Citations::mUseChasteImplementation, true);
        }

        // Reset to default
        Citations::mUseChasteImplementation = false;
        Citations::mCitations.clear();
#endif

        // Now test as part of PETSc,
        // PETSc implementation will be used if PETSc is >= 3.5
        // if not, this is something of a duplicate of the above test!
        {
            std::cout << "Testing PETSc implementation" << std::endl;

            // Turn on citations with argument - N.B. we direct PETSc to a different file to Chaste implementation.
            // N.B. This is needed twice - PETSc implementation reads the first one,
            // Chaste implementation masquerading as PETSc reads the second!
            CommandLineArgumentsMocker mocker("-citations " + output_petsc_file.GetAbsolutePath());

            PetscSetupUtils::CommonFinalize(); // This prints the citations to disk

// Check PETSc version - this is just because they reformatted their BibTex between versions, no change to function!
#if (PETSC_VERSION_MAJOR == 3 && (PETSC_VERSION_MINOR == 9 || PETSC_VERSION_MINOR == 10 || PETSC_VERSION_MINOR == 11)) //PETSc 3.9, 3.10 and 3.11
            FileFinder reference_citations("global/test/data/citations-2018.txt", RelativeTo::ChasteSourceRoot);
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 8) //PETSc 3.8
            FileFinder reference_citations("global/test/data/citations-2017.txt", RelativeTo::ChasteSourceRoot);
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7) //PETSc 3.7
            FileFinder reference_citations("global/test/data/citations-2016.txt", RelativeTo::ChasteSourceRoot);
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 6) //PETSc 3.6
            FileFinder reference_citations("global/test/data/citations-2015.txt", RelativeTo::ChasteSourceRoot);
#else
            // Use PETSc 3.4 older version (or mocked up version that matches 3.4)
            FileFinder reference_citations("global/test/data/citations.txt", RelativeTo::ChasteSourceRoot);
#endif
            FileComparison check_files(output_petsc_file, reference_citations, false); // false = not collective (this is after Finalize)
            check_files.CompareFiles();
        }
    }
};

#endif // TESTCITATIONS_HPP_
