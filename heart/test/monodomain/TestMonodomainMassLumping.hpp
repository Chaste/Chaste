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


#ifndef TESTMONODOMAINMASSLUMPING_HPP_
#define TESTMONODOMAINMASSLUMPING_HPP_

#include <cxxtest/TestSuite.h>
#include "LuoRudy1991BackwardEuler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestMonodomainMassLumping : public CxxTest::TestSuite
{

public:

    void TestCompareCubePlaneStimulus() throw(Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(20); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01,0.1,0.1);
        double spatial_step = 0.05;
        HeartConfig::Instance()->SetSlabDimensions(0.3, 0.3, 0.3, spatial_step); // 3mm edge cube meshed at 500um

        /*
         *  Standard solve
         */
        HeartConfig::Instance()->SetOutputDirectory("CompareCubeStandard");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CompareCubeStandard");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellMLBackwardEuler,3> cell_factory(-3e5, 1.0);

        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        DistributedVector standard_solution = monodomain_problem.GetSolutionDistributedVector();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();


        /*
         *  Mass lumping solve
         */
        HeartEventHandler::Reset();
        HeartConfig::Instance()->SetOutputDirectory("CompareCubeMassLumping");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CompareCubeMassLumping");
        HeartConfig::Instance()->SetUseMassLumping();

        MonodomainProblem<3> monodomain_problem_ml( &cell_factory );

        monodomain_problem_ml.Initialise();
        monodomain_problem_ml.Solve();

        DistributedVector mass_lumping_solution = monodomain_problem_ml.GetSolutionDistributedVector();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        // The idea is to check that the error stays O(h)
        double tolerance = 30*spatial_step;
        for (DistributedVector::Iterator index = standard_solution.Begin();
             index != standard_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(standard_solution[index], mass_lumping_solution[index], tolerance);
        }

    }

    void TestCompareCubePlaneStimulusOnlyPrecondLumping() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(20); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01,0.1,0.1);
        double spatial_step = 0.05;
        HeartConfig::Instance()->SetSlabDimensions(0.3, 0.3, 0.3, spatial_step); // 3mm edge cube meshed at 500um

        /*
         *  Standard solve
         */
        HeartConfig::Instance()->SetOutputDirectory("CompareCubeStandard");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CompareCubeStandard");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellMLBackwardEuler,3> cell_factory(-3e5, 1.0);

        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        DistributedVector standard_solution = monodomain_problem.GetSolutionDistributedVector();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();


        /*
         *  Preconditioning mass lumping solve
         */
        HeartEventHandler::Reset();
        HeartConfig::Instance()->SetOutputDirectory("CompareCubeMassLumping");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CompareCubeMassLumping");
        HeartConfig::Instance()->SetUseMassLumpingForPrecond();

        MonodomainProblem<3> monodomain_problem_ml( &cell_factory );

        monodomain_problem_ml.Initialise();
        monodomain_problem_ml.Solve();

        DistributedVector mass_lumping_solution = monodomain_problem_ml.GetSolutionDistributedVector();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        // Solutions should agree closely
        double tolerance = 1e-6;
        for (DistributedVector::Iterator index = standard_solution.Begin();
             index != standard_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(standard_solution[index], mass_lumping_solution[index], tolerance);
        }

    }


};

#endif /* TESTMONODOMAINMASSLUMPING_HPP_ */
