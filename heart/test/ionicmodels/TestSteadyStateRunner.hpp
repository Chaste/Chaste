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

#ifndef _TESTSTEADYSTATERUNNER_HPP_
#define _TESTSTEADYSTATERUNNER_HPP_

#include <cxxtest/TestSuite.h>

#include "CellProperties.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "VectorHelperFunctions.hpp" // For CopyToStdVector and DeleteVector

#include "Shannon2004Cvode.hpp"
#include "SteadyStateRunner.hpp"
#include "ZeroStimulus.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSteadyStateRunner : public CxxTest::TestSuite
{

public:
    void TestSteadyStateRunnerConverges(void)
    {
#ifdef CHASTE_CVODE
        //////////// DEFINE PARAMETERS ///////////////
        // Get the frequency
        double hertz = 1.0;

        ///////// END DEFINE PARAMETERS ////////////////////////

        // Setup a CVODE model that has empty solver and stimulus
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellShannon2004FromCellMLCvode(p_solver, p_stimulus));

        // Get it to use the default stimulus from CellML
        boost::shared_ptr<RegularStimulus> p_reg_stim = p_model->UseCellMLDefaultStimulus();

        { // Test that the steady state analysis throws a nice error if we try to run it with a non-RegularStimulus

            boost::shared_ptr<ZeroStimulus> p_stim(new ZeroStimulus());
            p_model->SetStimulusFunction(p_stim);

            SteadyStateRunner bad_steady_runner(p_model);
            TS_ASSERT_THROWS_THIS(bad_steady_runner.RunToSteadyState(),
                                  "Steady State approximations only work for models with RegularStimulus objects.");

            // Reset to a sensible stimulus function.
            p_model->SetStimulusFunction(p_reg_stim);
        }

        p_reg_stim->SetPeriod(1000.0 / hertz);

        // Note that increasing the strictness of the tolerances here (default is 1e-5, 1e-7)
        // actually leads to a faster convergence to steady state, as a model solved in
        // a sloppy way might never actually get close to its true limit cycle!
        p_model->SetTolerances(1e-6, 1e-8);

        /**
         * STEADY STATE PACING EXPERIMENT
         */
        /*
         * HOW_TO_TAG Cardiac/Cell Models
         * Get a cardiac cell model to (roughly) a steady state, given a regular stimulus, using the `SteadyStateRunner` class.
         */
        SteadyStateRunner steady_runner(p_model);

        bool result;

        // Here we don't reach steady state by max num paces
        steady_runner.SetMaxNumPaces(1u);
        result = steady_runner.RunToSteadyState();

        TS_ASSERT_EQUALS(result, false);

        // Here we do reach the steady state OK.
        steady_runner.SetMaxNumPaces(10000u);
        result = steady_runner.RunToSteadyState();

        TS_ASSERT_EQUALS(result, true);
        // Your mileage may vary. 32-bit machine, default build gives 520 evaluations, 484 on recent 64-bit CVODE.
        TS_ASSERT_LESS_THAN(steady_runner.GetNumEvaluations(), 550u);

        // For coverage
        TS_ASSERT_THROWS_THIS(steady_runner.SetMaxNumPaces(0u),
                              "Please set a maximum number of paces that is positive");

        // For coverage
        steady_runner.SuppressOutput();
#else
        std::cout << "CVODE must be enabled for the steady state runner to work." << std::endl;
#endif //_CHASTE_CVODE
    }
};

#endif //_TESTSTEADYSTATERUNNER_HPP_
