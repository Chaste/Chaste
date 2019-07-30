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

#ifndef TESTTIMEADAPTIVITYCONTROLLER_HPP_
#define TESTTIMEADAPTIVITYCONTROLLER_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractTimeAdaptivityController.hpp"

#include "PetscSetupAndFinalize.hpp"

// Very simple toy time-adaptivity controller
class ToyController : public AbstractTimeAdaptivityController
{
    double ComputeTimeStep(double currentTime, Vec currentSolution)
    {
        if (currentTime < 1)
        {
            return 0.1;
        }
        else if (currentTime < 2)
        {
            return 0.5;
        }
        else
        {
            return 1.5;
        }
    }

public:
    ToyController(double minDt, double maxDt)
        : AbstractTimeAdaptivityController(minDt, maxDt)
    {
    }
};

class TestTimeAdaptivityController : public CxxTest::TestSuite
{
public:
    void TestToyController()
    {
        ToyController controller(0.2,1.0);
        TS_ASSERT_EQUALS(controller.GetNextTimeStep(0.5,NULL), 0.2);
        TS_ASSERT_EQUALS(controller.GetNextTimeStep(1.5,NULL), 0.5);
        TS_ASSERT_EQUALS(controller.GetNextTimeStep(10 ,NULL), 1.0);
    }
};

#endif /*TESTTIMEADAPTIVITYCONTROLLER_HPP_*/
