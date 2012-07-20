/*

Copyright (c) 2005-2012, University of Oxford.
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
        // The CommandLineArguments instance is filled in by the cxxtest test suite runner.
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
        ///Causes memory failure (and seg fault) in PETSc 3.2 with MPICH-1
        PETSCEXCEPT(PetscFinalize());
        return true;
    }
};

static PetscSetup thisSetup;

#endif //_PETSCSETUPANDFINALIZE_HPP_
