/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTLATTICEBASEDUPDATERULES_HPP_
#define TESTLATTICEBASEDUPDATERULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "AdvectionCaUpdateRule.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestCaBasedUpdateRules : public AbstractCellBasedTestSuite
{
public:

    void TestDiffusionCaUpdateRuleConstructor()
    {
        // Create update rule using default input argument
        DiffusionCaUpdateRule<2> update_rule1;

        // Check the diffusion constant has been set correctly
        TS_ASSERT_DELTA(update_rule1.GetDiffusionConstant(), 1.0000, 1e-6);

        // Create update rule using different input argument
        DiffusionCaUpdateRule<3> update_rule2(3.1412);

        // Check the diffusion constant has been set correctly
        TS_ASSERT_DELTA(update_rule2.GetDiffusionConstant(), 3.1412, 1e-6);
    }

    void TestDiffusionCaUpdateRuleGetNewLocationOfCell() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create a cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_cell(new Cell(p_state, p_model));

        // Create lattice-based cell population
        std::vector<CellPtr> cells;
        cells.push_back(p_cell);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Create update rule
        DiffusionCaUpdateRule<2> update_rule(2.0);

        // Set the time step high enough that movement is guaranteed
        double dt = 1.0;

        // Test an exception is thrown when using a location index that does not correspond to a cell
        TS_ASSERT_THROWS_THIS(update_rule.GetNewLocationOfCell(0, cell_population, dt),
                              "There is no cell at the current location.");

        // Check the cell moves to the correct new location (random, but reproducible)
        unsigned new_location_index = update_rule.GetNewLocationOfCell(4, cell_population, dt);

        TS_ASSERT_EQUALS(new_location_index, 7u);
    }

    void TestArchiveDiffusionCaUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DiffusionCaUpdateRule.arch";

        {
            DiffusionCaUpdateRule<2> update_rule(1.25);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

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
            TS_ASSERT_DELTA((static_cast<DiffusionCaUpdateRule<2>*>(p_update_rule))->GetDiffusionConstant(), 1.25, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestAdvectionCaUpdateRuleGetNewLocationOfCell() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(5, 5, true); // 6*6 nodes

        // Create a line of cells along the bottom of the mesh

        // Create location indices
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<6; i++)
        {
            real_node_indices.push_back(i);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Create advection update rule: impose a flow 'north' with (mean) speed 2
        unsigned flow_direction = 0; // north
        AdvectionCaUpdateRule<2> update_rule(flow_direction, 2.0);

        // Take the time step large enough for flow-induced movement to be guaranteed
        double dt = 1.0;

        // Test an exception is thrown when using a location index that does not correspond to a cell
        TS_ASSERT_THROWS_THIS(update_rule.GetNewLocationOfCell(10, cell_population, dt),
                              "There is no cell at the current location.");

        // Call GetNewLocationOfCell() on each cell
        for (unsigned i=0; i<6; i++)
        {
            unsigned new_location_index = update_rule.GetNewLocationOfCell(i, cell_population, dt);
            TS_ASSERT_EQUALS(new_location_index, i+6);
        }
    }

    void TestArchiveAdvectionCaUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "AdvectionCaUpdateRule.arch";

        {
            AdvectionCaUpdateRule<2> update_rule(3, 1.25);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

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
            TS_ASSERT_DELTA((static_cast<AdvectionCaUpdateRule<2>*>(p_update_rule))->GetAdvectionSpeed(), 1.25, 1e-6);
            TS_ASSERT_EQUALS((static_cast<AdvectionCaUpdateRule<2>*>(p_update_rule))->GetAdvectionDirection(), 3u);

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestUpdateRuleOutputUpdateRuleInfo()
    {
        std::string output_directory = "TestCaUpdateRulesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with AdvectionCaUpdateRule
        AdvectionCaUpdateRule<2> advection_update_rule(3, 1.6);
        TS_ASSERT_EQUALS(advection_update_rule.GetIdentifier(), "AdvectionCaUpdateRule-2");

        out_stream advection_update_rule_parameter_file = output_file_handler.OpenOutputFile("advection_update_rule_results.parameters");
        advection_update_rule.OutputUpdateRuleInfo(advection_update_rule_parameter_file);
        advection_update_rule_parameter_file->close();

        std::string advection_update_rule_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + advection_update_rule_results_dir + "advection_update_rule_results.parameters cell_based/test/data/TestCaUpdateRules/advection_update_rule_results.parameters").c_str()), 0);

        // Test with DiffusionCaUpdateRule
        DiffusionCaUpdateRule<2> diffusion_update_rule(1.4);
        TS_ASSERT_EQUALS(diffusion_update_rule.GetIdentifier(), "DiffusionCaUpdateRule-2");

        out_stream diffusion_update_rule_parameter_file = output_file_handler.OpenOutputFile("diffusion_update_rule_results.parameters");
        diffusion_update_rule.OutputUpdateRuleInfo(diffusion_update_rule_parameter_file);
        diffusion_update_rule_parameter_file->close();

        std::string diffusion_update_rule_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + diffusion_update_rule_results_dir + "diffusion_update_rule_results.parameters cell_based/test/data/TestCaUpdateRules/diffusion_update_rule_results.parameters").c_str()), 0);
    }
};

#endif /*TESTLATTICEBASEDUPDATERULES_HPP_*/
