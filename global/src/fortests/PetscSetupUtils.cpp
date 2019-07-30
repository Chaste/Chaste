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

#include "PetscSetupUtils.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <petsc.h>

#include "ChasteBuildRoot.hpp"
#include "ChasteSyscalls.hpp"
#include "Citations.hpp"
#include "CommandLineArguments.hpp"
#include "Exception.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "PetscException.hpp"
#include "PetscTools.hpp"

#include "PetscAndChasteCitations.hpp"

#ifdef TEST_FOR_FPE
#include <fenv.h>
#include <signal.h>

void FpeSignalToAbort(int sig_num, siginfo_t* info, void* context)
{
    if (info->si_code == FPE_FLTDIV)
    {
        std::cerr << "SIGFPE: floating point exception was divide by zero.\n";
    }
    else if (info->si_code == FPE_FLTINV)
    {
        std::cerr << "SIGFPE: floating point exception was an invalid operation (like 0.0/0.0).\n";
    }
    else
    {
        std::cerr << "SIGFPE: unexpected error code.\n";
    }
}
#endif

void PetscSetupUtils::InitialisePetsc()
{
    // The CommandLineArguments instance is filled in by the cxxtest test suite runner.
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    PETSCEXCEPT(PetscInitialize(p_args->p_argc, p_args->p_argv, PETSC_NULL, PETSC_NULL));
    // Work around what seems to be an Intel compiler bug/quirk that makes the cache stale,
    // by using an explicit reset to ensure all code is aware we're running in parallel.
    PetscTools::ResetCache();
}

void PetscSetupUtils::CommonSetup()
{
    InitialisePetsc();

    if ((PETSC_VERSION_MAJOR < 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5)) // PETSc 3.4 or earlier
    {
        // Add some PETSc citations if we're not using the built-in citation mechanisms
        Citations::Register(PetscCitation1, &PetscCite1);
        Citations::Register(PetscCitation2, &PetscCite2);
    }
    // Add main Chaste citation
    Citations::Register(ChasteCitation, &ChasteCite);

    // Check that the working directory is correct, or many tests will fail
    std::string cwd = GetCurrentWorkingDirectory() + "/";
    if (strcmp(cwd.c_str(), ChasteSourceRootDir()) != 0)
    {
        // Change directory
        std::cout << std::endl
                  << "Changing directory from '" << cwd << "' to '" << ChasteSourceRootDir() << "'." << std::endl;
        EXPECT0(chdir, ChasteSourceRootDir());
        std::cout << "CWD now: " << GetCurrentWorkingDirectory() << std::endl;
    }

#ifdef TEST_FOR_FPE
    // Give all PETSc enabled tests the ability to trap for divide-by-zero
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
    // Catch all SIGFPE signals and convert them to exceptions (before PETSc gets to them)
    struct sigaction sa;
    sa.sa_sigaction = FpeSignalToAbort;
    sa.sa_flags = SA_RESETHAND | SA_SIGINFO;
    sa.sa_restorer = 0;
    sigaction(SIGFPE, &sa, NULL);
#endif
}

void PetscSetupUtils::CommonFinalize()
{
    // This does nothing if we are on a new PETSc, just allows Chaste to print citations instead in this case
    Citations::Print();

    PETSCEXCEPT(PetscFinalize());
}

void PetscSetupUtils::ResetStatusCache()
{
    PetscTools::ResetCache();
}
