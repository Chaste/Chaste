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

#ifndef TESTMONODOMAINWITHTIMEADAPTIVITY_HPP_
#define TESTMONODOMAINWITHTIMEADAPTIVITY_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "ReplicatableVector.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "Warnings.hpp"


/* HOW_TO_TAG Cardiac/Solver
 * Run using (simple, user-defined) time-adaptivity
 */

// Toy controller which just goes alters the timestep from 0.01ms to 1ms after a given
// threshold time.
class FixedTimeAdaptivityController : public AbstractTimeAdaptivityController
{
private:
    double mThresholdTime;

    double ComputeTimeStep(double currentTime, Vec currentSolution)
    {
        if (currentTime < mThresholdTime)
        {
            return 0.01; // ms
        }
        else
        {
            return 1;
        }
    }

public:
    FixedTimeAdaptivityController(double thresholdTime)
      : AbstractTimeAdaptivityController(0.01, 1.0),
        mThresholdTime(thresholdTime)
    {
        assert(thresholdTime > 0);
    }
};


class TestMonodomainWithTimeAdaptivity : public CxxTest::TestSuite
{
public:
    void TestWithCube()
    {
        HeartConfig::Instance()->SetPrintingTimeStep(1.0);
        HeartConfig::Instance()->SetSimulationDuration(3); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-600.0*1000);

        //////////////////////////////////////////////////////////////////////////
        // run original simulation - no adaptivity, dt=0.01 all the way through
        //////////////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("MonoWithTimeAdaptivity/OrigNoAdapt");
        MonodomainProblem<3> problem(&cell_factory);
        problem.Initialise();
        problem.Solve();

        //HeartEventHandler::Headings();
        //HeartEventHandler::Report();

        Vec solution = problem.GetSolution();
        int index; //dummy
        double min_non_adaptive;
        double max_non_adaptive;
        VecMin(solution, &index, &min_non_adaptive);
        VecMax(solution, &index, &max_non_adaptive);
        //std::cout << "Non adaptive: range at final time: " << min_non_adaptive << "mV to " << max_non_adaptive << "mV\n";

        //////////////////////////////////////////////////////////////////////////
        // run adaptive simulation - dt=0.01 for first 2ms, then dt=1
        //////////////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("MonoWithTimeAdaptivity/SimpleAdapt");
        MonodomainProblem<3> adaptive_problem(&cell_factory);
        FixedTimeAdaptivityController controller(2.0);
        adaptive_problem.SetUseTimeAdaptivityController(true, &controller);
        adaptive_problem.Initialise();
        adaptive_problem.Solve();
        adaptive_problem.SetUseTimeAdaptivityController(false);

        //HeartEventHandler::Headings();
        //HeartEventHandler::Report();   // note: adaptive is slower in this short sim due to second matrix assemble

        Vec adaptive_solution = adaptive_problem.GetSolution();
        double min_adaptive;
        double max_adaptive;
        VecMin(adaptive_solution, &index, &min_adaptive);
        VecMax(adaptive_solution, &index, &max_adaptive);
        ///std::cout << "Adaptive:     range at final time: " << min_adaptive << "mV to " << max_adaptive << "mV\n";


        // compare

        TS_ASSERT_DELTA(min_non_adaptive, 21.9062, 1e-3);
        TS_ASSERT_DELTA(max_non_adaptive, 28.8345, 1e-3);
        TS_ASSERT_DELTA(min_adaptive, 19.9010, 1e-3);
        TS_ASSERT_DELTA(max_adaptive, 25.6083, 1e-3);
    }

    void TestWithChebyshevAndFixedIterations()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetPrintingTimeStep(1.0);
        HeartConfig::Instance()->SetSimulationDuration(3); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");

        HeartConfig::Instance()->SetKSPSolver("chebychev");
        HeartConfig::Instance()->SetUseFixedNumberIterationsLinearSolver(true, 30);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-600.0*1000);

        double min_non_adaptive;
        double max_non_adaptive;
        double min_adaptive;
        double max_adaptive;

        try
        {
            //////////////////////////////////////////////////////////////////////////
            // run original simulation - no adaptivity, dt=0.01 all the way through
            //////////////////////////////////////////////////////////////////////////
            {
                HeartConfig::Instance()->SetOutputDirectory("MonoWithTimeAdaptivityOrigNoAdapt");
                MonodomainProblem<3> problem(&cell_factory);
                problem.Initialise();
                problem.Solve();

                HeartEventHandler::Headings();
                HeartEventHandler::Report();

                Vec solution = problem.GetSolution();
                int index; //dummy
                VecMin(solution, &index, &min_non_adaptive);
                VecMax(solution, &index, &max_non_adaptive);
                //std::cout << "Non adaptive: range at final time: " << min_non_adaptive << "mV to " << max_non_adaptive << "mV\n";

                //This test assumes no exception.  Requires PETSc > 2.3.2
                TS_ASSERT_DELTA(min_non_adaptive, 21.9062, 1e-3);
                TS_ASSERT_DELTA(max_non_adaptive, 28.8346, 1e-3);
            }


            //////////////////////////////////////////////////////////////////////////
            // run adaptive simulation - dt=0.01 for first 2ms, then dt=1
            //////////////////////////////////////////////////////////////////////////
            {
                HeartConfig::Instance()->SetOutputDirectory("MonoWithTimeAdaptivitySimpleAdapt");
                MonodomainProblem<3> adaptive_problem(&cell_factory);
                FixedTimeAdaptivityController controller(2.0);
                adaptive_problem.SetUseTimeAdaptivityController(true, &controller);
                adaptive_problem.Initialise();
                adaptive_problem.Solve();

                HeartEventHandler::Headings();
                HeartEventHandler::Report();   // note: adaptive is slower in this short sim due to second matrix assemble

                Vec adaptive_solution = adaptive_problem.GetSolution();
                int index; //dummy
                VecMin(adaptive_solution, &index, &min_adaptive);
                VecMax(adaptive_solution, &index, &max_adaptive);
                ///std::cout << "Adaptive:     range at final time: " << min_adaptive << "mV to " << max_adaptive << "mV\n";

                // compare
                //This test assumes no exception.  Requires PETSc > 2.3.2
                TS_ASSERT_DELTA(min_adaptive, 19.9009, 1e-3);
                TS_ASSERT_DELTA(max_adaptive, 25.6082, 1.25e-3);
                // Note that the threshold was increased to 1.25e-3 due to a minor change in PETSc. See #2997.
            }
        }
        catch (Exception& e)
        {
            if (e.GetShortMessage() == "Chebyshev with fixed number of iterations is known to be broken in PETSc <= 2.3.2")
            {
                WARNING(e.GetShortMessage());
            }
            else
            {
                TS_FAIL(e.GetShortMessage());
            }
        }
    }
};

#endif /*TESTMONODOMAINWITHTIMEADAPTIVITY_HPP_*/
