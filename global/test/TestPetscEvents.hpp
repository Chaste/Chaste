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
