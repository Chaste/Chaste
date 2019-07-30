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

#ifndef TESTPDESIMULATIONTIME_HPP_
#define TESTPDESIMULATIONTIME_HPP_

#include <cxxtest/TestSuite.h>
#include "PdeSimulationTime.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPdeSimulationTime : public CxxTest::TestSuite
{
public:
    void TestTime()
    {
        PdeSimulationTime::SetTime(1.0453346);
        TS_ASSERT_EQUALS(PdeSimulationTime::GetTime(), 1.0453346);

        PdeSimulationTime::SetPdeTimeStepAndNextTime(0.025, 1.0453346+0.025);
        TS_ASSERT_EQUALS(PdeSimulationTime::GetPdeTimeStep(), 0.025);
        TS_ASSERT_EQUALS(PdeSimulationTime::GetPdeTimeStepInverse(), 1.0/0.025);

        TS_ASSERT_EQUALS(PdeSimulationTime::GetNextTime(), 1.0453346+0.025);
    }
    void TestAssertionWithNonRationalTimeStep()
    {
        PdeSimulationTime::SetTime(5.0);
        for (unsigned i=1; i<100; i++)
        {
            double end_time = 5.0+0.1*i;
            /*
             *  With non-rational floating point time step we expect a small "error" in
             *  time-step (or end-points).  We simulate the worse possible situation (or something
             *  like it but slightly worse since neither 0.1 or 0.1+0.1*epsilon are exactly representable).
             *
             */

            PdeSimulationTime::SetPdeTimeStepAndNextTime((1.0+DBL_EPSILON)*0.1, end_time);
            PdeSimulationTime::SetTime(end_time);
        }
    }
};

#endif /*TESTPDESIMULATIONTIME_HPP_*/
