/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef TESTSSACHSE_HPP_
#define TESTSSACHSE_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>
#include "RegularStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

#include "CellProperties.hpp"

#include "sachse_moreno_abildskov_2008_b.hpp"
#include "sachse_moreno_abildskov_2008_bBackwardEuler.hpp"

class TestSsachse: public CxxTest::TestSuite
{
public:
    void TestCompareNormalBackwarwdEuler(){
        double magnitude=-40.0;
        double duration = 2.0; // ms
        double when = 50.0; // ms
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        
	Cellsachse_moreno_abildskov_2008_bFromCellML normal(p_solver, p_stimulus);
	Cellsachse_moreno_abildskov_2008_bFromCellMLBackwardEuler backwardEuler(p_solver, p_stimulus);

        double start_time = 0.0;
        double end_time = 1000.0; //One second in milliseconds
	double sampling_interval = 1.0;

	OdeSolution normalSolution = normal.Compute(start_time, end_time, sampling_interval);
	OdeSolution backwardEulerSolution = backwardEuler.Compute(start_time, end_time, sampling_interval);

	std::vector<double> normalVoltages = normalSolution.GetVariableAtIndex(normal.GetVoltageIndex());
	std::vector<double> backwardEulerVoltages = backwardEulerSolution.GetVariableAtIndex(backwardEuler.GetVoltageIndex());

        TS_ASSERT_EQUALS(normalVoltages.size(), backwardEulerVoltages.size());
        for (unsigned i = 0; i < normalVoltages.size(); i++)
        {
            TS_ASSERT_DELTA(normalVoltages[i], backwardEulerVoltages[i], 0.1);
        }

    }
};

#endif // TESTSSACHSE_HPP_