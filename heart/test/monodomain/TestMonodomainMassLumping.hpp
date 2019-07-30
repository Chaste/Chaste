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

    void TestCompareCubePlaneStimulus()
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

    void TestCompareCubePlaneStimulusOnlyPrecondLumping()
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
