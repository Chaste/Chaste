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
 * [[PageOutline]]
 *
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
#include "AbstractCvodeCell.hpp"
#include "CellProperties.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"
#include "Shannon2004Cvode.hpp"
#include "SteadyStateRunner.hpp"
/* This test is always run sequentially (never in parallel)*/
#include "FakePetscSetup.hpp"

/* Now we define the test class, which must inherit from {{{CxxTest::TestSuite}}}
 * as usual, and the (public) test method
 */
class TestSingleCellSimulationTutorial : public CxxTest::TestSuite
{
public:
    void TestShannonSimulation()
    {
/* CVODE is still an optional Chaste dependency, but it is highly recommended for
 * working with single cell simulations. This tutorial code will only run if CVODE is installed and enabled
 * (see InstallCvode and ChasteGuides/CmakeBuildGuide). */
#ifdef CHASTE_CVODE
        /*
         * == Defining a CVODE model ==
         *
         * Setup a CVODE model that has empty solver and stimulus
         * This is necessary to initialise the cell model.
         *
         * If you want to define your own stimulus without using the default one,
         * you can define it here instead of giving it an empty stimulus:
         *
         * {{{boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-25.5,2.0,50.0,500));}}}
         *
         * the parameters are magnitude, duration, period, and start time of stimulus.
         */
        boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new CellShannon2004FromCellMLCvode(p_solver, p_stimulus));

        /*
         * Once the model is set up we can tell it to use the the default stimulus from CellML,
         * (if one has been labelled, you get an exception if not), and return it.
         *
         * NB. You could automatically check whether one is available with:
         *
         * {{{p_model->HasCellMLDefaultStimulus()}}}
         *
         */
        boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();

        /*
         * Now you can modify certain parameters of the stimulus function, such as the period
         */
        p_regular_stim->SetPeriod(1000.0);

        /*
         * == Numerical Considerations ==
         *
         * Cardiac cell models can be pretty tricky to deal with, as they are very stiff and sometimes full
         * of singularities.
         *
         * One of the first things you want to ensure is that the maximum timestep CVODE can take is less
         * than or equal to the duration of the stimulus. Otherwise CVODE could evaluate the right-hand side
         * before and after the stimulus, and never see it (giving you a cell that never does anything).
         * This can be done using something like:
         *
         * {{{double max_timestep = p_regular_stim->GetDuration();}}}
         *
         * instead of the declaration of `max_timestep` below. In this tutorial we want an answer that is
         * refined in time to give an accurate upstroke velocity. So we make the maximum timestep even
         * smaller, to match the printing timestep. A rough rule of thumb would be to use the smaller of
         * stimulus duration and printing time step as your CVODE maximum timestep. But note CVODE will still
         * give sensible answers if printing timestep is less than maximum timestep, it might just decide
         * to interpolate the output instead of evaluate it directly, if it thinks it will be
         * accurate enough to meet your tolerances.
         *
         *
         * A common error from CVODE is '''TOO_MUCH_WORK''', this means CVODE tried to exceed the maximum number of
         * internal time steps it is allowed to do. You can try using the method `SetMaxSteps` to change
         * the default (500) to a larger value, with a command like:
         *
         * {{{p_model->SetMaxSteps(1e5);}}}
         *
         * We have found that 1e5 should be enough for a single pace of all models we've tried so far,
         * but if you were running for a long time (e.g. 1000 paces in one Solve call) you would need to increase this.
         *
         * Another common error from CVODE is:
         * '''the error test failed repeatedly or with |h| = hmin.'''
         *
         * Since we don't change hmin (and it defaults to a very small value), this generally means the
         * ODE system has got to a situation where refining the timestep is not helping the convergence.
         *
         * This generally indicates that you are hitting some sort of singularity, or divide by zero, in
         * the model. Unfortunately cardiac models are full of these, they can sometimes be manually edited out
         * by changing the cellML file, for instance using [http://en.wikipedia.org/wiki/L%27H%C3%B4pital%27s_rule L'Hopital's rule].
         *
         * In this case, one other thing you can try is to change the absolute and relative
         * tolerances of the CVODE solver, the default being (1e-5,1e-7), although it isn't clear whether
         * refining sometimes makes things worse for models with singularities,
         * as CVODE goes to look for trouble in areas with steep gradients.
         *
         * For this particular test, we are going to specify quite strict tolerances, so that the test gets the same results
         * on different versions of CVODE and different compilers.
         */
        p_model->SetTolerances(1e-8, 1e-8);

        /*
         * By default we use an analytic Jacobian for CVODE cells (where available - see [wiki:ChasteGuides/CodeGenerationFromCellML]
         * for instructions on how to provide one using Maple). In some cases (the Hund-Rudy model particularly being one) the
         * analytic Jacobian contains effectively divide-by-zero entries, even at resting potential. If you observe
         * CVODE errors when trying to run simulations, it can be worth switching off the analytic Jacobian and resorting
         * to a numerical approximation (as happens by default if no analytic Jacobian is available). This can be done with the
         * following command:
         *
         * {{{p_model->ForceUseOfNumericalJacobian();}}}
         *
         */

        /*
         * == Changing Parameters in the Cell Model ==
         *
         * You can also change any parameters that are labelled in the cell model.
         *
         * Instructions for annotating parameters can be found at [wiki:ChasteGuides/CodeGenerationFromCellML]
         *
         * Here we show how to change the parameter dictating the maximal conductance of the IKs current.
         * Note this call actually leaves it unchanged from the default,
         * you can experiment with changing it and examine the impact on APD.
         */
        p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.07);

        /*
         * == Running model to steady state ==
         *
         * Now we run the model to steady state.
         * You can detect for steady state alternans by giving it true as a second parameter
         * {{{SteadyStateRunner steady_runner(p_model, true);}}}
         *
         * You may change the number of maximum paces the runner takes. The default is 1e5.
         */
        SteadyStateRunner steady_runner(p_model);
        steady_runner.SetMaxNumPaces(100u);
        bool result;
        result = steady_runner.RunToSteadyState();

        /*
         * Check that the model has NOT reached steady state
         * (this model needs more than 100 paces to reach steady state).
         *
         */
        TS_ASSERT_EQUALS(result, false);

        /*
         * == Getting detail for paces of interest ==
         *
         * Now we solve for the number of paces we are interested in.
         *
         * The absolute values of start time and end time are typically only relevant for the stimulus, in general
         * nothing else on the right-hand side of the equations uses time directly.
         *
         * i.e. if you have a `RegularStimulus` of period 1000ms then you would get exactly the same results
         * calling Solve(0,1000,...) twice, as you would calling Solve(0,1000,...) and Solve(1000,2000,...).
         *
         * Single cell results can be very sensitive to the sampling time step, because of the steepness of the upstroke.
         *
         * For example, try changing the line below to 1 ms. The upstroke velocity that is detected will change
         * from 339 mV/ms to around 95 mV/ms. APD calculations will only ever be accurate to sampling timestep
         * for the same reason.
         */
        double max_timestep = 0.1;
        p_model->SetMaxTimestep(max_timestep);

        double sampling_timestep = max_timestep;
        double start_time = 0.0;
        double end_time = 1000.0;
        OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

        /*
         * `p_model` retains the state variables at the end of `Solve`, if you call `Solve` again the state
         * variables will evolve from their new state, not the original initial conditions.
         *
         * Write the data out to a file.
         */
        solution.WriteToFile("TestCvodeCells", "Shannon2004Cvode", "ms");

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

        std::cout << "APD = " << apd << "ms" << std::endl;
        std::cout << "Upstroke velocity = " << upstroke_velocity << "mV/ms" << std::endl;

        /*
         * Here we just check that the values are equal to the ones we expect,
         * with appropriate precision to pass on different versions of CVODE.
         *
         * (These reference values were generated with tolerances of Abs=1e-12, Rel=1e-12.)
         */
        TS_ASSERT_DELTA(apd, 212.411, 1e-2);
        TS_ASSERT_DELTA(upstroke_velocity, 338.704, 1.25);

        /* CVODE is still an optional dependency for Chaste, but is required for this tutorial.
         * If CVODE is not installed this tutorial will
         * not do anything, but we can at least alert the user to this.*/
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#endif /*TESTSINGLECELLSIMULATIONTUTORIAL_HPP_*/
