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

//Note that this name deliberately clashes.  We want to avoid a test being both sequential and parallel.
#ifndef _PETSCSETUPANDFINALIZE_HPP_
#define _PETSCSETUPANDFINALIZE_HPP_

/**
 * This file is designed to be included by any test suites that DO NOT use PETSc.
 * It does the PETSc initialisation and finalisation at the top of a test suite, in order to make
 * sure that only one process is left to run the actual tests.
 */

#include <cxxtest/GlobalFixture.h>

#include <cstdlib>
#include <petsc.h>

#include "PetscSetupUtils.hpp"

#include "PetscException.hpp"

class PetscSetup : public CxxTest::GlobalFixture
{
public:

    /**
     * Run the standard setup method for PETSc, but then fake being a single process.
     * Only the process with PETSc rank zero continues beyond setUpWorld(); the others
     * exit gracefully.
     *
     * Note that we need to keep MPI initialized so we can use MPI_Wtime etc.
     *
     * @return true (by CxxTest convention)
     */
    bool setUpWorld()
    {
        PetscSetupUtils::CommonSetup();

        // Get rank
        PetscInt my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        // Make sure that only one process proceeds into the test itself
        if (my_rank != 0)
        {
            PETSCEXCEPT(PetscFinalize());
            exit(0);
        }

        // Fool PETSc into thinking it was only run on one process.
        // This ensures any barriers etc. don't cause deadlock.
        PETSC_COMM_WORLD = MPI_COMM_SELF;
        PetscSetupUtils::ResetStatusCache();

        return true;
    }

    /**
     * Clean up PETSc on the master process after running tests.
     * @return true (by CxxTest convention)
     */
    bool tearDownWorld()
    {
        // Remind PETSc there were originally more processes so MPI is finalized properly.
        PETSC_COMM_WORLD = MPI_COMM_WORLD;
        PetscSetupUtils::CommonFinalize();
        return true;
    }
};

static PetscSetup thisSetup;

#endif //_PETSCSETUPANDFINALIZE_HPP_
