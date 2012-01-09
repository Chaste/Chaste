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

#ifndef TESTEXECUTABLESUPPORT_HPP_
#define TESTEXECUTABLESUPPORT_HPP_

#include <cxxtest/TestSuite.h>
#include "ExecutableSupport.hpp"
#include "CommandLineArguments.hpp"

/**
 * This doesn't really do any testing - it's just here for coverage.
 */
class TestExecutableSupport : public CxxTest::TestSuite
{
public:
    void TestStaticMethods() throw(Exception)
    {
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        ExecutableSupport::StandardStartup(p_args->p_argc, p_args->p_argv);
        ExecutableSupport::SetOutputDirectory("TestExecutableSupport");
        std::string msg("This is not an error, it's just for coverage.");
        ExecutableSupport::PrintError(msg, true);
        ExecutableSupport::PrintError(msg);
        ExecutableSupport::Print(msg);
        ExecutableSupport::WriteProvenanceInfoFile();
        ExecutableSupport::WriteMachineInfoFile("write_test");
        ExecutableSupport::FinalizePetsc();

        TS_ASSERT_EQUALS(ExecutableSupport::EXIT_OK, 0);
        TS_ASSERT_EQUALS(ExecutableSupport::EXIT_ERROR, 1);
        TS_ASSERT_EQUALS(ExecutableSupport::EXIT_BAD_ARGUMENTS, 2);
    }
};

#endif /* TESTEXECUTABLESUPPORT_HPP_ */
