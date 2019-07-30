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
#include "FileFinder.hpp"

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

    void TestCheckpointingGeneratesMultipleDirectories()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_checkpoint.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 0.2);
        TS_ASSERT(HeartConfig::Instance()->GetCheckpointSimulation());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetCheckpointTimestep(), 0.1);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "SaveBi2DCheckPoint");

        // Test that two directories have been created
        OutputFileHandler handler("");
        FileFinder directory_finder("SaveBi2DCheckPoint_checkpoints/0.1ms/SaveBi2DCheckPoint",RelativeTo::ChasteTestOutput);
        TS_ASSERT(directory_finder.IsDir());
        directory_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.1ms/SaveBi2DCheckPoint_0.1ms",RelativeTo::ChasteTestOutput);
        TS_ASSERT(directory_finder.IsDir());
        directory_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint",RelativeTo::ChasteTestOutput);
        TS_ASSERT(directory_finder.IsDir());
        directory_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms",RelativeTo::ChasteTestOutput);
        TS_ASSERT(directory_finder.IsDir());

        FileFinder file_finder("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/AbstractCardiacProblem_mSolution.h5",RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.IsFile());
        file_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/ChasteParameters.xml",RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.IsFile());
        file_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/mesh.ele",RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.IsFile());
        file_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/mesh.edge",RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.IsFile());
        file_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/mesh.node",RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.IsFile());
        file_finder.SetPath("SaveBi2DCheckPoint_checkpoints/0.2ms/SaveBi2DCheckPoint_0.2ms/archive.arch",RelativeTo::ChasteTestOutput);
        TS_ASSERT(file_finder.IsFile());
    }

    void TestCheckpointingGeneratesSameResults()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_checkpoint.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "SaveBi2DCheckPoint");

        CardiacSimulation simulation_compare("heart/test/data/xml/bidomain2d_compare_with_checkpoint.xml");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputDirectory(), "SaveBi2DCheckPointCompare");

        TS_ASSERT( CompareFilesViaHdf5DataReader("SaveBi2DCheckPoint", "bidomain2d", true,
                                                 "SaveBi2DCheckPointCompare", "bidomain2d", true));
    }
};

#endif /*TESTCHECKPOINTING_HPP_*/
