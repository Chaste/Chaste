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

#ifndef TESTCHECKPOINTING_HPP_
#define TESTCHECKPOINTING_HPP_

#include <cxxtest/TestSuite.h>

#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "BidomainProblem.hpp"
#include "CardiacSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "OutputFileHandler.hpp"

class TestCheckpointing : public CxxTest::TestSuite
{
public:

    /*
     *  This test makes sure that a simulation of x seconds can be run in multiple steps (SetSimulationTime()
     * plus Solve()) and return the same results that a single call to Solve()
     */
    void TestMultipleCallsToProblemSolve()
    {
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain3d");

        ///////////////////////////////////////////////////////////////////
        // Multiples calls to Solve()
        ///////////////////////////////////////////////////////////////////
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-600.0*1000);
        BidomainProblem<3> bidomain_problem_multiple( &cell_factory );

        bidomain_problem_multiple.Initialise();
        HeartConfig::Instance()->SetSimulationDuration(0.1);
        bidomain_problem_multiple.Solve();
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        bidomain_problem_multiple.Solve();

        ///////////////////////////////////////////////////////////////////
        // Single call to Solve()
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("Bidomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain3d");

        BidomainProblem<3> bidomain_problem_single( &cell_factory );
        bidomain_problem_single.Initialise();
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        bidomain_problem_single.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector single_solution = bidomain_problem_single.GetSolutionDistributedVector();
        DistributedVector multiple_solution = bidomain_problem_multiple.GetSolutionDistributedVector();
        DistributedVector::Stripe single_vm(single_solution,0);
        DistributedVector::Stripe single_phie(single_solution,1);
        DistributedVector::Stripe multiple_vm(multiple_solution,0);
        DistributedVector::Stripe multiple_phie(multiple_solution,1);
        for (DistributedVector::Iterator index = single_solution.Begin();
             index != single_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(single_vm[index], multiple_vm[index], 1e-8);
            TS_ASSERT_DELTA(single_phie[index], multiple_phie[index], 1e-8);
        }

    }

    void TestCheckpointingGeneratesMultipleDirectories() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_checkpoint.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 0.2);
        TS_ASSERT(HeartConfig::Instance()->GetCheckpointSimulation());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCheckpointTimestep(), 0.1);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "SaveBi2DCheckPoint");

        // Test that two directories have been created
        OutputFileHandler handler("");
        EXPECT0(system, "test -d " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.1ms/SaveBi2DCheckPoint");
        EXPECT0(system, "test -d " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.1ms/SaveBi2DCheckPoint_0.1ms");
        EXPECT0(system, "test -d " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint");
        EXPECT0(system, "test -d " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms");

        // Test the content of one of them
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/AbstractCardiacProblem_mSolution.h5");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/ChasteParameters.xml");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/mesh.ele");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/mesh.edge");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/mesh.node");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "/SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/archive.arch");
    }

    void TestCheckpointingGeneratesSameResults()
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_checkpoint.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "SaveBi2DCheckPoint");

        CardiacSimulation simulation_compare("heart/test/data/xml/bidomain2d_compare_with_checkpoint.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "SaveBi2DCheckPointCompare");

        TS_ASSERT( CompareFilesViaHdf5DataReader("SaveBi2DCheckPoint", "bidomain2d", true,
                                                 "SaveBi2DCheckPointCompare", "bidomain2d", true));
    }
};

#endif /*TESTCHECKPOINTING_HPP_*/
