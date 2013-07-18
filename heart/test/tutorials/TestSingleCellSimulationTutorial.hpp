/*

Copyright (c) 2005-2013, University of Oxford.
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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 */
#ifndef TESTSINGLECELLSIMULATIONTUTORIAL_HPP_
#define TESTSINGLECELLSIMULATIONTUTORIAL_HPP_
/*
 * = An example showing how to run a single cell simulation =
 *
 * == Introduction ==
 *
 * In this tutorial we run a single cell simulation,
 * showing:
 *  * how to load a cardiac cell (using CVODE - best solver to use for single cell simulations).
 *  * how to define the stimulus (using the default from CellML).
 *  * run the model to steady state using the {{{SteadyStateRunner}}}.
 *  * solve the model for the number of paces of interest.
 *  * Write to voltage data to a file.
 *  * Use {{{CellProperties}}} to output information such as APD and upstroke velocity.
 *
 *
 *
 * The first thing to do is to include the headers.
 */
#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
/* This test is always run sequentially (never in parallel)*/
#include "FakePetscSetup.hpp"


/* Now we define the test class, which must inherit from {{{CxxTest::TestSuite}}}
 * as usual, and the (public) test method
 */
class TestSingleCellTutorial : public CxxTest::TestSuite
{
public:
    void TestLuoRudySimulation() throw(Exception)
    {
/* Make sure this code is only run if CVODE is installed/enabled on the computer */
#ifdef CHASTE_CVODE
        /*
         * == Defining a CVODE model ==
         *
         * Setup a CVODE model that has empty solver and stimulus
         * This is necessary to initialise the cell model.
         *
         * If you want to define your own stimulus without using the default one,
         * you can define it here instead of giving it an empty stimulus:
         * {{{boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-25.5,2.0,50.0,500);}}}
         * the parameters are magnitude, duration, start time and period of stimulus
         */
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellShannon2004FromCellMLCvode(p_solver, p_stimulus));

        /* Once the model is set up we can get the the default stimulus from CellML */
        p_model->UseCellMLDefaultStimulus();

        /* And retrieve a pointer to this. We need to cast it to a {{{RegularStimulus}}} */
        boost::shared_ptr<RegularStimulus> p_reg_stim =
                             boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction());

        /* Now you can modify certain parameters of the stimulus function, such as the period
         * {{{p_reg_stim->SetPeriod(500.0);}}}
         */

        /*
         * You can change the absolute and relative tolerances of the solver, the default being (1e-5,1e-7)
         * {{{p_model->SetTolerances(1e-6,1e-8);}}}
         */


        /*
         * == Running model to steady state ==
         *
         * Now we run the model to steady state.
         * You can detect for steady state alternans by giving it true as a second parameter
         * {{{SteadyStateRunner steady_runner(p_model, true);}}}
         * You may change the number of maximum paces the runner takes. The default is 1e5.
         */
        SteadyStateRunner steady_runner(p_model);
        steady_runner.SetMaxNumPaces(100u);
        bool result;
        result = steady_runner.RunToSteadyState();

        /*
         * Check that the model has NOT reached steady state.
         * The model needs more than a 100 paces to reach steady state.
         */
        TS_ASSERT_EQUALS(result,false);

        /*
         * == Solving model for paces of interest ==
         *
         * Now we solve for the number of paces we are interested in
         * max_timestep and sampling time step are the same for CVODE
         * The start time and end time are only relevant for the stimulus
         *
         */
        double max_timestep = 1.0;
        double sampling_time = max_timestep;
        double start_time = 0.0;
        double end_time = 1000.0;
        OdeSolution solution = p_model->Solve(start_time, end_time, max_timestep, sampling_time);

        /*
         * Write the data out to a file.
         */
        solution.WriteToFile("TestCvodeCells","Shannon2004Cvode","ms");

        /*
         * == Calculating APD and Upstroke Velocity ==
         *
         * Calculate APD and upstroke velocity using {{{CellProperties}}}
         */
        unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
        CellProperties cell_props(voltages, solution.rGetTimes());

        double apd = cell_props.GetLastActionPotentialDuration(90);
        double upstroke_velocity = cell_props.GetLastMaxUpstrokeVelocity();
        /*
         * Here we just check that the values are equal to the ones we expect with 1e-2 precision.
         */
        TS_ASSERT_DELTA(apd,212.42,1e-2);
        TS_ASSERT_DELTA(upstroke_velocity,94.85,1e-2);

#else
        /* CVODE is not enable or installed*/
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#endif /*TESTANOTHERBIDOMAINSIMULATIONTUTORIAL_HPP_*/
