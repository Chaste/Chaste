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

#ifndef TESTMONODOMAINMASSLUMPINGLONG_HPP_
#define TESTMONODOMAINMASSLUMPINGLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "LuoRudy1991BackwardEuler.hpp"
#include "MonodomainProblem.hpp"
#include "SimpleStimulus.hpp"
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        // Stimulate the apex
        if (pNode->rGetLocation()[0] > 0.94)
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

    void TestCompareRealisticGeometry()
    {
        HeartConfig::Reset();
        HeartConfig::Instance()->SetSimulationDuration(50); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01,0.1,0.1);

        double spatial_step = 0.05;
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/UCSD_heart");

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
