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

#ifndef TESTPETSCEVENTS_HPP_
#define TESTPETSCEVENTS_HPP_

#include "PetscSetupAndFinalize.hpp"
#include <petsc.h>

class TestPetscEvents : public CxxTest::TestSuite
{
public:
    void TestEvents()
    {
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        PetscLogEvent my_event;
        PetscLogEventRegister("My first event", 0, &my_event);
#else
        PetscEvent my_event;
        PetscLogEventRegister(&my_event, "My first event", 0);

#endif
        (void)PetscLogEventBegin(my_event,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(my_event,0,0,0,0);

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        PetscLogEvent my_event2;
        PetscLogEventRegister("My second event", 0, &my_event2);
#else
        PetscEvent my_event2;
        PetscLogEventRegister(&my_event2, "My second event", 0);
#endif

        (void)PetscLogEventBegin(my_event2,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(my_event2,0,0,0,0);

        (void)PetscLogEventBegin(24,0,0,0,0);
        for (unsigned i=0; i<1000000; i++);
        (void)PetscLogEventEnd(24,0,0,0,0);

        //PetscLogPrintDetailed(MPI_COMM_WORLD, filename);
    }
    // this test should be run on the command line with -log_summary
    // to check that a summary is given
};

#endif /*TESTPETSCEVENTS_HPP_*/
