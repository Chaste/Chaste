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

#ifndef _PETSCSETUPANDFINALIZE_HPP_
#define _PETSCSETUPANDFINALIZE_HPP_

/**
 * This file is designed to be included by any test suites that use PETSc.
 * It controls the PETSc initialisation and finalisation.
 */

#ifdef TEST_FOR_FPE
#include <fenv.h>
#include <signal.h>
#endif

#include <cxxtest/GlobalFixture.h>
#include <petsc.h>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <unistd.h>
#include <iostream>

#include "Exception.hpp"
#include "PetscException.hpp"
#include "CommandLineArguments.hpp"
#include "ChasteBuildRoot.hpp"
#include "GetCurrentWorkingDirectory.hpp"

#ifdef TEST_FOR_FPE
void FpeSignalToAbort(int sig_num, siginfo_t* info, void* context )
{
       if ( info->si_code == FPE_FLTDIV)
       {
           std::cerr<<"SIGFPE: floating point exception was divide by zero.\n";
       }
       else if ( info->si_code == FPE_FLTINV)
       {
           std::cerr<<"SIGFPE: floating point exception was an invalid operation (like 0.0/0.0).\n";
       }
       else
       {
           std::cerr<<"SIGFPE: unexpected error code.\n";
       }
}
#endif

class PetscSetup : public CxxTest::GlobalFixture
{
public:

    /** Standard setup method for PETSc. */
    bool setUpWorld()
    {
        /**
         * The cxxtest_argc and cxxtest_argv variables are global, and filled in
         * from the arguments passed to the cxxtest test suite runner.
         */
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        PETSCEXCEPT(PetscInitialize(p_args->p_argc, p_args->p_argv,
                                    PETSC_NULL, PETSC_NULL) );

        // Check that the working directory is correct, or many tests will fail
        std::string cwd = GetCurrentWorkingDirectory() + "/";
        if (strcmp(cwd.c_str(), ChasteBuildRootDir()) != 0)
        {
#define COVERAGE_IGNORE
            // Change directory
            std::cout << std::endl << "Changing directory from '" << cwd << "' to '" << ChasteBuildRootDir() << "'." << std::endl;
            EXPECT0(chdir, ChasteBuildRootDir());
            std::cout << "CWD now: " << GetCurrentWorkingDirectory() << std::endl;
#undef COVERAGE_IGNORE
        }

#ifdef TEST_FOR_FPE
        // Give all PETSc enabled tests the ability to trap for divide-by-zero
        feenableexcept(FE_DIVBYZERO | FE_INVALID );
        // Catch all SIGFPE signals and convert them to exceptions (before PETSc gets to them)
        struct sigaction sa;
        sa.sa_sigaction = FpeSignalToAbort;
        sa.sa_flags = SA_RESETHAND|SA_SIGINFO;
        sa.sa_restorer = 0;
        sigaction(SIGFPE, &sa, NULL);
#endif
        return true;
    }
    /** Clean up PETSc after running all tests. */
    bool tearDownWorld()
    {
        PETSCEXCEPT(PetscFinalize());
        return true;
    }
};

static PetscSetup thisSetup;

#endif //_PETSCSETUPANDFINALIZE_HPP_
