/*

Copyright (C) University of Oxford, 2005-2011

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


// Toy controller which just goes alters the timestep from 0.01ms to 1ms after a given
// threshold time.
class FixedTimeAdaptivityController : public AbstractTimeAdaptivityController
{
private:
    double mThresholdTime;

    double ComputeTimeStep(double currentTime, Vec currentSolution)
    {
        if(currentTime < mThresholdTime)
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
    void TestWithCube() throw(Exception)
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

        TS_ASSERT_DELTA(min_non_adaptive, 22.0383, 1e-3);
        TS_ASSERT_DELTA(max_non_adaptive, 29.0697, 1e-3);
        TS_ASSERT_DELTA(min_adaptive, 19.8749, 1e-3);
        TS_ASSERT_DELTA(max_adaptive, 25.0398, 1e-3);
    }

    void TestWithChebyshevAndFixedIterations() throw(Exception)
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

                TS_ASSERT_DELTA(min_non_adaptive, 22.0383, 1e-3);
                TS_ASSERT_DELTA(max_non_adaptive, 29.0697, 1e-3);
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
                TS_ASSERT_DELTA(min_adaptive, 19.8749, 1e-3);
                TS_ASSERT_DELTA(max_adaptive, 25.0398, 1e-3);
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
