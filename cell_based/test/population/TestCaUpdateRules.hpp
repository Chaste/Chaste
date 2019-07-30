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
#ifndef TESTCAUPDATERULES_HPP_
#define TESTCAUPDATERULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCaUpdateRule.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "AbstractCaSwitchingUpdateRule.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "PottsMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCaUpdateRules : public AbstractCellBasedTestSuite
{
public:

    void TestDiffusionCaUpdateRuleIn2d()
    {
        // Set the timestep and size of domain to let us calculate the probabilities of movement
        double delta_t = 1;
        double delta_x = 1;
        double diffusion_parameter = 0.1;

        // Create an update law system
        DiffusionCaUpdateRule<2> diffusion_update_rule;

        // Test get/set methods
        TS_ASSERT_DELTA(diffusion_update_rule.GetDiffusionParameter(), 0.5, 1e-12);

        diffusion_update_rule.SetDiffusionParameter(1.0);

        TS_ASSERT_DELTA(diffusion_update_rule.GetDiffusionParameter(), 1.0, 1e-12);

        diffusion_update_rule.SetDiffusionParameter(diffusion_parameter);

        // Test EvaluateProbability()

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 1u, p_diff_type);

        // Specify where cells lie here we have one cell on the bottom left site
        std::vector<unsigned> location_indices;
        location_indices.push_back(0u);

        // Create cell population
        CaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices);

        // Note we just pass a pointer to the only cell as DiffusionUpdateRule is independent of cell type.
        CellPtr p_cell = cell_population.rGetCells().front();
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(0,1,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(0,6,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(0,5,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);

        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,1,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,6,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,10,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);

        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(6,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(7,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(10,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(12,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(15,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(16,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(17,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);

        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(24,19,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(24,18,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(24,23,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
    }

    void TestDiffusionCaUpdateRuleIn2dWithMultipleCells()
    {
        // Set the timestep and size of domain to let us calculate the probabilities of movement
        double delta_t = 1;
        double delta_x = 1;
        double diffusion_parameter = 0.1;

        // Create an update law system
        DiffusionCaUpdateRule<2> diffusion_update_rule;

        // Test get/set methods
        TS_ASSERT_DELTA(diffusion_update_rule.GetDiffusionParameter(), 0.5, 1e-12);

        diffusion_update_rule.SetDiffusionParameter(1.0);

        TS_ASSERT_DELTA(diffusion_update_rule.GetDiffusionParameter(), 1.0, 1e-12);

        diffusion_update_rule.SetDiffusionParameter(diffusion_parameter);

        // Test EvaluateProbability()

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 10u, p_diff_type);

        // Specify where cells lie here we have one cell on the bottom left site
        std::vector<unsigned> location_indices;

        // Two cells will be added in locations 0 and 1
        location_indices.push_back(0u);
        location_indices.push_back(0u);
        location_indices.push_back(1u);
        location_indices.push_back(1u);

        // Then other cells are added to the lattice to check if the probabilities are still the same
        for (unsigned i=4; i<10; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        CaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices, 2);

        // Note we just pass a pointer to the first cell as DiffusionUpdateRule is independent of cell type.
        CellPtr p_cell = cell_population.rGetCells().front();
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(0,1,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(0,6,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(0,5,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);

        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,1,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,6,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,10,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(5,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);

        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(6,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(7,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(10,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(12,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(15,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(16,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(17,11,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);

        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(24,19,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(24,18,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/4.0,1e-6);
        TS_ASSERT_DELTA(diffusion_update_rule.EvaluateProbability(24,23,cell_population, delta_t, delta_x, p_cell),diffusion_parameter*delta_t/delta_x/delta_x/2.0,1e-6);
    }

    void TestArchiveDiffusionCaUpdateRule()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DiffusionCaUpdateRule.arch";

        {
            DiffusionCaUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            update_rule.SetDiffusionParameter(1.0);

            // Serialize via pointer to most abstract class possible
            AbstractCaUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractCaUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            TS_ASSERT_DELTA((static_cast<DiffusionCaUpdateRule<2>*>(p_update_rule))->GetDiffusionParameter(), 1.0, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestUpdateRuleOutputUpdateRuleInfo()
    {
        std::string output_directory = "TestCaUpdateRulesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with VolumeConstraintPottsUpdateRule
        DiffusionCaUpdateRule<2> diffusion_update_rule;
        diffusion_update_rule.SetDiffusionParameter(1.0);

        TS_ASSERT_EQUALS(diffusion_update_rule.GetIdentifier(), "DiffusionCaUpdateRule-2");

        out_stream diffusion_update_rule_parameter_file = output_file_handler.OpenOutputFile("diffusion_update_rule_results.parameters");
        diffusion_update_rule.OutputUpdateRuleInfo(diffusion_update_rule_parameter_file);
        diffusion_update_rule_parameter_file->close();

        // Compare the generated file in test output with a reference copy in the source code.
        FileFinder generated = output_file_handler.FindFile("diffusion_update_rule_results.parameters");
        FileFinder reference("cell_based/test/data/TestCaUpdateRules/diffusion_update_rule_results.parameters",
                RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }

    /*
     * Now test the switching rules.
     */
    void TestRandomCaSwitchingUpdateRuleIn2d()
    {
        // Set the timestep and size of domain to let us calculate the probabilities of movement
        double delta_t = 0.1;
        double delta_x = 1;
        double switching_parameter = 0.1;

        // Create an update law system
        RandomCaSwitchingUpdateRule<2> random_switching_update_rule;

        // Test get/set methods
        TS_ASSERT_DELTA(random_switching_update_rule.GetSwitchingParameter(), 0.5, 1e-12);

        random_switching_update_rule.SetSwitchingParameter(1.0);

        TS_ASSERT_DELTA(random_switching_update_rule.GetSwitchingParameter(), 1.0, 1e-12);

        random_switching_update_rule.SetSwitchingParameter(switching_parameter);

        // Test EvaluateProbability()

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 6u, p_diff_type);

        // Specify where cells lie here we have cells on the bottom two rows
        std::vector<unsigned> location_indices;

        for (unsigned i=0; i<6; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        CaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices);

        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(0,1,cell_population, delta_t, delta_x),switching_parameter*delta_t,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(0,6,cell_population, delta_t, delta_x),switching_parameter*delta_t,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(0,5,cell_population, delta_t, delta_x),switching_parameter*delta_t,1e-6);
        // Note this is independent of node index and population so even returns for nodes not in the mesh
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(UNSIGNED_UNSET,UNSIGNED_UNSET,cell_population, delta_t, delta_x),switching_parameter*delta_t,1e-6);
    }

    void TestArchiveRandomCaSwitchingUpdateRule()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "RandomCaSwitchingUpdateRule.arch";

        {
            RandomCaSwitchingUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            update_rule.SetSwitchingParameter(1.0);

            // Serialize via pointer to most abstract class possible
            AbstractCaSwitchingUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractCaSwitchingUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            TS_ASSERT_DELTA((static_cast<RandomCaSwitchingUpdateRule<2>*>(p_update_rule))->GetSwitchingParameter(), 1.0, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestSwitchingUpdateRuleOutputUpdateRuleInfo()
    {
        std::string output_directory = "TestCaSwitchingUpdateRulesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with RandomCaSwitchingUpdateRule
        RandomCaSwitchingUpdateRule<2> random_switching_update_rule;
        random_switching_update_rule.SetSwitchingParameter(1.0);

        TS_ASSERT_EQUALS(random_switching_update_rule.GetIdentifier(), "RandomCaSwitchingUpdateRule-2");

        out_stream random_switching_update_rule_parameter_file = output_file_handler.OpenOutputFile("random_switching_update_rule_results.parameters");
        random_switching_update_rule.OutputUpdateRuleInfo(random_switching_update_rule_parameter_file);
        random_switching_update_rule_parameter_file->close();

        // Compare the generated file in test output with a reference copy in the source code.
        FileFinder generated = output_file_handler.FindFile("random_switching_update_rule_results.parameters");
        FileFinder reference("cell_based/test/data/TestCaUpdateRules/random_switching_update_rule_results.parameters",
                RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }
};

#endif /*TESTCAUPDATERULES_HPP_*/
