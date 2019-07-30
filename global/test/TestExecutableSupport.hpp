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

#ifndef TESTEXECUTABLESUPPORT_HPP_
#define TESTEXECUTABLESUPPORT_HPP_

#include <cxxtest/TestSuite.h>
#include "ExecutableSupport.hpp"
#include "CommandLineArguments.hpp"

// NB: PETSc is set up in the test itself, so don't include a setup header!

/**
 * This doesn't really do any testing - it's just here for coverage.
 */
class TestExecutableSupport : public CxxTest::TestSuite
{
public:
    void TestStaticMethods()
    {
        const std::string output_dir("TestExecutableSupport");
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        ExecutableSupport::StandardStartup(p_args->p_argc, p_args->p_argv);
        ExecutableSupport::StartupWithoutShowingCopyright(p_args->p_argc, p_args->p_argv);
        ExecutableSupport::SetOutputDirectory(output_dir);
        std::string msg("This is not an error, it's just for coverage.");
        ExecutableSupport::PrintError(msg, true);
        ExecutableSupport::PrintError(msg);
        ExecutableSupport::Print(msg);
        ExecutableSupport::WriteProvenanceInfoFile();
        ExecutableSupport::WriteMachineInfoFile("write_test");

        OutputFileHandler handler(output_dir, false);
        ExecutableSupport::SetOutputDirectory(handler.GetOutputDirectoryFullPath());

        ExecutableSupport::FinalizePetsc();

        TS_ASSERT_EQUALS(ExecutableSupport::EXIT_OK, 0);
        TS_ASSERT_EQUALS(ExecutableSupport::EXIT_ERROR, 1);
        TS_ASSERT_EQUALS(ExecutableSupport::EXIT_BAD_ARGUMENTS, 2);
    }
};

#endif /* TESTEXECUTABLESUPPORT_HPP_ */
