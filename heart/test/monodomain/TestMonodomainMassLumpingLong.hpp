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

#ifndef TESTMONODOMAINMASSLUMPINGLONG_HPP_
#define TESTMONODOMAINMASSLUMPINGLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "LuoRudy1991BackwardEuler.hpp"
#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"


template<class CELL>
class SmallBenchmarkStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{

private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    SmallBenchmarkStimulusHeartCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1000.0*500, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        // Stimulate the apex
        if (GetMesh()->GetNode(node)->rGetLocation()[0] > 0.94)
        {
            return new CELL(mpSolver, mpStimulus);
        }
        else
        {
            return new CELL(mpSolver, mpZeroStimulus);
        }
    }
};


class TestMonodomainMassLumping : public CxxTest::TestSuite
{

public:

    void TestCompareRealisticGeometry() throw(Exception)
    {
        HeartConfig::Reset();
        HeartConfig::Instance()->SetSimulationDuration(50); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01,0.1,0.1);

        double spatial_step = 0.05;
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/heart");

        /*
         *  Standard solve
         */
        HeartConfig::Instance()->SetOutputDirectory("CompareRealisticGeometryStandard");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CompareRealisticGeometryStandard");

        SmallBenchmarkStimulusHeartCellFactory<CellLuoRudy1991FromCellMLBackwardEuler> cell_factory;

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
        HeartConfig::Instance()->SetOutputDirectory("CompareRealisticGeometryMassLumping");
        HeartConfig::Instance()->SetOutputFilenamePrefix("CompareRealisticGeometryMassLumping");
        HeartConfig::Instance()->SetUseMassLumping();

        MonodomainProblem<3> monodomain_problem_ml( &cell_factory );

        monodomain_problem_ml.Initialise();
        monodomain_problem_ml.Solve();

        DistributedVector mass_lumping_solution = monodomain_problem_ml.GetSolutionDistributedVector();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        // The idea is to check that the error stays O(h)
        double tolerance = 100*spatial_step;
        for (DistributedVector::Iterator index = standard_solution.Begin();
             index != standard_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(standard_solution[index], mass_lumping_solution[index], tolerance);
        }

    }

};


#endif /* TESTMONODOMAINMASSLUMPINGLONG_HPP_ */
